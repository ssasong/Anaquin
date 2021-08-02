#include <assert.h>
#include <unistd.h>
#include <iostream>
#include <sys/types.h>
#include <htslib/bgzf.h>
#include "KmerIndex.h"
#include "Kallisto.hpp"
#include "ProcessReads.h"
#include "tools/tools.hpp"
#include "KmerIterator.hpp"
#include "tools/errors.hpp"
#include <boost/format.hpp>
#include "tools/samtools.hpp"
#include "sequins/rna/rna.hpp"
#include "sequins/meta/meta.hpp"
#include "writers/bam_writer.hpp"
#include "sequins/genomics/genomics.hpp"

using namespace Anaquin;

// Cache for combined BAM
static std::map<Thread, std::unique_ptr<char>> s1_, s2_, q1_, q2_;

// Extremly unlikely but it's only 1MB
static const int MAX_READ_LENGTH = 1000000;

class Partition;

struct FAIndex
{
    // Partition for this index
    std::map<Thread, std::shared_ptr<Partition>> part;
    
    // Running statistics by threads
    std::map<Thread, KStats> stats;
    
    // Kallisto index
    std::shared_ptr<KmerIndex> index;
};

struct KData
{
    Library lib;
    
    FAIndex __R__;
    
    // Temp directory for the current exectution
    Path __tmp__;
    
    KOptions o;
};

class Partition
{
    public:
        Partition(const Path &, Thread, KData *);
        ~Partition() { close(); }

        void close()
        {
            for (auto &i : o1)   { i.second->close();    }
            for (auto &i : o2)   { i.second->close();    }
            for (auto &i : bamW) { bgzf_close(i.second); }
            bamW.clear();
        }
    
        void writeB(Bin k, const bam1_t *x)
        {
            if (bam_write1(bamW.at(k), x) == -1)
            {
                throw std::runtime_error("bam_write1 failed");
            }
        }
    
        void write(KData *,
                   Bin, const char *, const char *, const char *, void *,
                        const char *, const char *, const char *, void *, bool);
    
        /*
         * Partitioned reads to FASTQ
         */
    
        std::map<Bin, FileName> c1, c2;
        std::map<Bin, std::shared_ptr<std::ofstream>> o1, o2;

        /*
         * Partitioned reads to BAM
         */
    
        std::map<Bin, BGZF *>   bamW;
        std::map<Bin, FileName> bamF;
};

typedef std::map<unsigned, std::shared_ptr<Partition>> PartitionReads;

void Partition::write(KData *d,
                      Bin k, const char *r1, const char *s1, const char *q1, void *b1,
                             const char *r2, const char *s2, const char *q2, void *b2, bool rev)
{
    auto writeF = [&](const char *r, const char *s, const char *q, std::shared_ptr<std::ofstream> o, bool rev)
    {
        assert(o);
        std::string ss = s;
        if (rev) { complement(ss); }
        
        std::stringstream str;
        str << "@" + std::string(r) + "\n";
        str << ss + "\n";
        str << "+\n";
        str << std::string(q) + "\n";
        
        const auto x = compressGZ(str.str());
        o->write(x.data(), x.size());
    };

    if (d->o.writeBAM())
    {
        if (!d->o.onlySeqLad || k == GR || k == LD)
        {
            writeB(k, (bam1_t *) b1);
            writeB(k, (bam1_t *) b2);
        }
    }
    else
    {
        writeF(r1, s1, q1, o1[k], rev);
        writeF(r2, s2, q2, o2[k], rev);
    }
}

Partition::Partition(const Path &p, Thread tID, KData *d)
{
    forBin([&](Bin c)
    {
        bool x;
        
        switch (d->o.prod)
        {
            case Product::RNA:      { x = RValid(c); break; }
            case Product::Meta:     { x = MValid(c); break; }
            case Product::Genomics: { x = GValid(c); break; }
        }
        
        const auto f1 = p + "/partition_B" + std::to_string(c) + "_" + S0(tID) + "_1.fq.gz.tmp";
        const auto f2 = p + "/partition_B" + std::to_string(c) + "_" + S0(tID) + "_2.fq.gz.tmp";
        o1[c] = std::shared_ptr<std::ofstream>(new std::ofstream(c1[c] = f1, std::ios::binary | std::ios::out));
        o2[c] = std::shared_ptr<std::ofstream>(new std::ofstream(c2[c] = f2, std::ios::binary | std::ios::out));
    });
}

static std::shared_ptr<KmerIndex> KIndex(const FileName &file, KmerLen k)
{
    ProgramOptions o;
    o.k = k;
    o.index = file;
    
    auto x = std::shared_ptr<KmerIndex>(new KmerIndex(o));
    x->load(o);
    return x;
}

// Build Kallisto index from a FASTA file
static FileName KBuildIndex(const FileName &file, unsigned k)
{
    ProgramOptions opt;
    
    opt.k = ::Kmer::k = k;
    opt.index = tmpFile();
    opt.transfasta.push_back(file);
    
    KmerIndex index(opt);
    index.BuildTranscripts(opt);
    index.write(opt.index);
    
    if (index.dbGraph.contigs.empty() || !index.kmap.size())
    {
        throw std::runtime_error("Failed to build index for " + file);
    }
    
    return opt.index;
}

// Initialise Kallisto. Required before running anything else
static KData KInit(const KOptions &o)
{
    KData x;
    
    // Always initialize in case it has been used (like ladder recalibration)
    x.__R__ = FAIndex();

    x.o = o;
    assert(x.o.prod == Product::RNA || x.o.prod == Product::Meta || x.o.prod == Product::Genomics);

    x.__tmp__ = tmpPath();
    removeD(x.__tmp__);
    createD(x.__tmp__);

    x.__R__.index = KIndex(KBuildIndex(o.index, o.k), o.k);

    if (x.__R__.index->target_names_.empty())
    {
        throw std::runtime_error("No sequence found in reference index: " + o.index);
    }

    for (auto i = 0u; i < x.__R__.index->target_names_.size(); i++)
    {
        const auto &seq = x.__R__.index->target_names_[i];
        
        StandardID std;
        
        switch (x.o.prod)
        {
            case Product::RNA:      { std = RSeq2Std(seq); break; }
            case Product::Meta:     { std = MSeq2Std(seq); break; }
            case Product::Genomics: { std = GSeq2Std(seq); break; }
        }

        assert(!std.empty());

        for (auto j = 0; j < o.thr; j++)
        {
            x.__R__.stats[j].stds.insert(std);
            x.__R__.stats[j].seqs.insert(seq);
            
            if (x.o.prod == Product::Genomics)
            {
                auto isRSeq = [&](const SequinID &x)
                {
                    return isEnd(x, "_R");
                };
                
                if (isRSeq(seq))
                {
                    x.__R__.stats[j].rSeqs[std] = seq;
                }
                else
                {
                    x.__R__.stats[j].aSeqs[std].insert(seq);
                }
            }
        }
    }
    
    return x;
}

struct EdgeMatch
{
    // n-th k-mer
    Base i;
    
    // Position matching in the reference
    Base p;
};

static bool matchIndex(Count nMatch, Count nNMatch, const KData *d)
{
    const auto mp = static_cast<float>(nMatch) / (nMatch + nNMatch);
    return nMatch > nNMatch || mp >= d->o.rule;
}

// Match a read by indexing
static Bin KMatch(Bin (*f)(const std::string &),
                  Thread tID,
                  const char *r,
                  const char *s,
                  void *b,
                  std::map<SequinID, Count> &ms,
                  bool &isForw,
                  Base &dist,
                  KData *d)
{
    Base start, end;
    
    // Our only index
    auto &x = d->__R__;

    if (!x.part.count(tID))
    {
        d->__R__.part.insert(std::pair<unsigned, std::shared_ptr<Partition>>(
                          tID, std::shared_ptr<Partition>(
                          new Partition(d->__tmp__, tID, d))));
    }
    
    std::string tmp;
    
    if (d->o.flipBefore)
    {
        tmp = s;
        complement(tmp);
        s = tmp.data();
    }

    // Number of k-mers in the read
    const auto n = strlen(s) - d->o.k + 1;

    auto __match__ = [&](std::shared_ptr<KmerIndex> index)
    {
        ms.clear();
        KmerIterator ki(s), ke;
        
        Count nMatch  = 0;
        Count nNMatch = 0;
        
        SequinKM shared, uniqs;
        
        auto j = 0;
        bool isLastMatched = false;
        
        // Edge matching for first and last
        std::shared_ptr<EdgeMatch> fm, lm;

        // Matched sequins?
        std::set<SequinID> seqs;

        for (int i = 0; ki != ke; ++i, ++ki)
        {
            // Never skip the first k-mer
            if (i && j != d->o.skipKM)
            {
                j++;
                continue;
            }
            
            j = 0;
            auto search = index->kmap.find(ki->first.rep());

            // Where the k-mer starts
            const auto l = i;
            
            // Where the k-mer finishes
            //const auto u = l + d->o.k - 1;

            // Can we find this k-mer?
            if (search != index->kmap.end())
            {
                isLastMatched = true;
                
                nMatch++;
                const auto km = search->first.toString();

                for (const auto &tran : index->dbGraph.contigs[search->second.contig].transcripts)
                {
                    seqs.insert(index->target_names_[tran.trid]);
                    
                    auto findPos = [&]()
                    {
                        auto x = index->findPosition(tran.trid, search->first, search->second);
                        return x.second ? x.first : x.first - d->o.k;
                    };
                    
                    if (!fm)
                    {
                        fm = std::shared_ptr<EdgeMatch>(new EdgeMatch());
                        fm->i = i; fm->p = findPos();
                    }
                    else
                    {
                        if (!lm) { lm = std::shared_ptr<EdgeMatch>(new EdgeMatch()); }
                        lm->i = i; lm->p = findPos();
                    }
                }

                // Unique matching?
                const auto isUniq = seqs.size() == 1;
                
                for (const auto &seq : seqs)
                {
                    ms[seq]++;

                    if (!isUniq) { shared[seq][km]++; }
                    else         { uniqs [seq][km]++;  }
                }
            }
            else
            {
                isLastMatched = false;
                nNMatch++;
            }
        }

        auto copy = [&](const SequinKM &x1, SequinKM &x2)
        {
            for (const auto &i : x1)
            {
                for (const auto &j : i.second)
                {
                    x2[i.first][j.first] += j.second;
                }
            }
        };
        
        // Only matched a single k-mer?
        if (fm && !lm) { lm = fm; }

        // Match this index?
        if (matchIndex(nMatch, nNMatch, d))
        {
            copy(uniqs,  x.stats[tID].uniqs);
            copy(shared, x.stats[tID].shared);
            
            x.stats[tID].nMatch++;
            
            if (!fm || !lm)
            {
                // Assume the whole sequence
                start = 0;
                
                // Assume the whole sequence
                end = strlen(s);
            }
            else
            {
                lm = !lm ? fm : lm;
                
                // Forward strand?
                isForw = fm->p < lm->p;
                
                // Estimated beginning of the read relative to sequin (like on IGV)
                start = isForw ? (fm->p - fm->i) : lm->p - (n - lm->i - 1);
                
                // Estimated ending of the read relative to sequin (like on IGV)
                end = start + strlen(s);
            }
            
            extern std::map<Anaquin::SequinID, long> __KMSize__;

            /*
             * d2 can be negative. For example, the code would count the 8 skips in 128M8S.
             */
            
            const auto d1 = start;
            const auto d2 = __KMSize__[*(seqs.begin())] - end;
            
            // Closest distance to sequin edges
            dist = std::max((Base) 0, std::min(d1, d2));
            
            return true;
        }
        
        // Doesn't match this index
        else
        {
            ms.clear();
            x.stats[tID].nNMatch++;
            return false;
        }
    };
    
    return (__match__(x.index)) ? (ms.empty() ? GR: f(max(ms))) : ES;
}

void KPartition(void *p,
                Thread tID,
                const char *r1,
                const char *s1,
                const char *q1,
                void       *b1,
                bool       rc1,
                const char *r2,
                const char *s2,
                const char *q2,
                void       *b2,
                bool       rc2)
{
    auto d = (KData *) p;
    
    if ((tID == 0) && !d->lib.heads())
    {
        d->lib.addInfo(r1, s1);
    }
    
    std::string tmp1, tmp2, tmp3, tmp4;
    
    if (false)
    {
        assert(!s1 && !q1 && !s2 && !q2);
        assert(s1_.count(tID) && q1_.count(tID));
        assert(s2_.count(tID) && q2_.count(tID));

        bam2seq (s1_.at(tID).get(), static_cast<const bam1_t *>(b1));
        bam2qual(q1_.at(tID).get(), static_cast<const bam1_t *>(b1));
        bam2seq (s2_.at(tID).get(), static_cast<const bam1_t *>(b2));
        bam2qual(q2_.at(tID).get(), static_cast<const bam1_t *>(b2));
        
        s1 = s1_[tID].get(); q1 = q1_[tID].get();
        s2 = s2_[tID].get(); q2 = q2_[tID].get();

        if (rc1)
        {
            tmp1 = std::string(q1);
            std::reverse(tmp1.begin(), tmp1.end());
            q1 = tmp1.c_str();
            tmp2 = s1;
            tmp2 = revcomp(tmp2);
            s1 = tmp2.c_str();
        }
        
        if (rc2)
        {
            tmp3 = std::string(q2);
            std::reverse(tmp3.begin(), tmp3.end());
            q2 = tmp3.c_str();
            tmp4 = s2;
            tmp4 = revcomp(tmp4);
            s2 = tmp4.c_str();
        }

        assert(strlen(s1) && strlen(s2));
        assert(strlen(s1) == strlen(q1));
        assert(strlen(s2) == strlen(q2));
    }

    // Counting for all matches (all k-mers)
    std::map<SequinID, Count> R1, R2, F1, F2;

    // Forward orientation?
    bool isF1, isF2;
    
    // Position for trimming
    Base d1, d2;
    
    auto match = [&](FAIndex &i, std::map<SequinID, Count> &m1, std::map<SequinID, Count> &m2)
    {
        std::pair<Bin, Bin> p;
        
        switch (d->o.prod)
        {
            case Product::Genomics:
            {
                p = std::pair<Bin, Bin>(KMatch(GBin, tID, r1, s1, b1, m1, isF1, d1, d),
                                        KMatch(GBin, tID, r2, s2, b2, m2, isF2, d2, d));
                break;
            }
                
            case Product::RNA:
            {
                p = std::pair<Bin, Bin>(KMatch(RBin, tID, r1, s1, b1, m1, isF1, d1, d),
                                        KMatch(RBin, tID, r2, s2, b2, m2, isF2, d2, d));
                break;
            }
                
            case Product::Meta:
            {
                p = std::pair<Bin, Bin>(KMatch(MBin, tID, r1, s1, b1, m1, isF1, d1, d),
                                        KMatch(MBin, tID, r2, s2, b2, m2, isF2, d2, d));
                break;
            }
        }
        
        auto shouldTrim = [](Bin x)
        {
            switch (x)
            {
                case Bin::MT:
                case Bin::MS:
                case Bin::HL:
                case Bin::HP:
                case Bin::SV:
                case Bin::IM:
                case Bin::GR:
                case Bin::SO:
                case Bin::MI: { return true;  }
                default:      { return false; }
            }
        };
        
        const auto shouldTrim1 = shouldTrim(p.first);
        const auto shouldTrim2 = shouldTrim(p.second);

        if (shouldTrim1 && shouldTrim2 && m1.size() == 1 && m2.size() == 1 && (d1 <= 15 || d2 <= 15))
        {
            return std::pair<Bin, Bin>(VC, VC);
        }
        
        return p;
    };

    auto fill = [&](Bin k, FAIndex &i, const std::map<SequinID, Count> &m1, const std::map<SequinID, Count> &m2)
    {
        const auto x1 = !m1.empty() ? max(m1) : "-";
        const auto x2 = !m2.empty() ? max(m2) : "-";
        
        auto addRead = [&]()
        {
            auto r1_ = first(r1, " ");
                 r1_ = isSubstr(r1_, "/") ? noLast(r1_, "/") : r1_;
            auto r2_ = first(r2, " ");
                 r2_ = isSubstr(r2_, "/") ? noLast(r2_, "/") : r2_;

            assert(!r1_.empty() && !r2_.empty());
            
            i.stats[tID].r1[k][r1_] = x1;
            i.stats[tID].r2[k][r2_] = x2;
        };
        
        // Always required for ladder calibration
        //if (d->o.writeReads || k == Bin::LD)
        if (k == Bin::LD)
        {
            addRead();
        }

        i.stats[tID].sqc[x1]++;
        i.stats[tID].sqc[x2]++;
        i.stats[tID].c1[k]++;
        i.stats[tID].c2[k]++;

        auto flip = !(k == ES || k == VC || k == IF || k == LD);
       
        // Flipping is only on genomics
        if (d->o.prod == Product::RNA || d->o.prod == Product::Meta)
        {
           flip = false;
        }
 
        // Force flipping off?
        if (!d->o.flip)
        {
            flip = false;
        }

        if (i.part.count(tID))
        {
            i.part[tID]->write(d, k, r1, s1, q1, b1, r2, s2, q2, b2, flip);
        }
    };

    auto write = [&](Bin x1, Bin x2)
    {
        if ((x1 != Bin::SV && x2 == Bin::SV) || (x1 == Bin::SV && x2 != Bin::SV))
        {
            fill(ES, d->__R__, R1, R2);
        }
        else if (x1 == Bin::VC || x2 == Bin::VC)
        {
            fill(VC, d->__R__, R1, R2);
        }
        else
        {
            if (x1 != x2)
            {
                // Write ambigious reads to the bin for sample reads
                fill(ES, d->__R__, R1, R2);
            }
            else
            {
                fill(x1, d->__R__, R1, R2);
            }
        }
    };
    
    const auto R = match(d->__R__, R1, R2);
    write(R.first, R.second);
}

// Merge the internal Kallisto's results to more accessible KStats
static KStats KMerge(KData &d)
{
    auto &x = d.__R__;
    auto r = x.stats[0];
    
    // Close all connections
    for (const auto &i : d.__R__.part) { i.second->close(); }

    // For each thread, merge with all previous threads...
    for (auto i = 1; i < x.stats.size(); i++)
    {
        auto addX1X2 = [](KStats &x1, KStats &x2)
        {
            assert(x1.aSeqs.size() == x2.aSeqs.size());
            
            // We'll just need to add x2
            KStats x = x1;
            
            auto f = [&](KStats &x1, KStats &x2)
            {
                x1.nMatch  += x2.nMatch;
                x1.nNMatch += x2.nNMatch;
                
                x1.uniqs  = add(x1.uniqs, x2.uniqs);
                x1.shared = add(x1.shared, x2.shared);
                
                x1.c1  = add(x1.c1,  x2.c1);
                x1.c2  = add(x1.c2,  x2.c2);
                x1.sqc = add(x1.sqc, x2.sqc);
                
                for (auto i = 0; i <= (int)VC; i++)
                {
                    x1.r1[(Bin)i] = add(x1.r1[(Bin)i], x2.r1[(Bin)i]);
                    x1.r2[(Bin)i] = add(x1.r2[(Bin)i], x2.r2[(Bin)i]);
                }
            };
            
            f(x, x2);
            return x;
        };
        
        r = addX1X2(r, x.stats[i]);
    }
    
    for (const auto &i : x.part)
    {
        forBin([&](Bin c)
        {
            if (d.o.writeBAM())
            {
                if (c != ES)
                {
                    if (!d.o.onlySeqLad || c == GR || c == LD)
                    {
                        r.f1[c].push_back(i.second->bamF.at(c));
                    }
                }
            }
            else
            {
                if (r.c1.count(c) && r.c2.count(c))
                {
                    r.f1[c].push_back(i.second->c1.at(c));
                    r.f2[c].push_back(i.second->c2.at(c));
                }
            }
        });
    }
    
    s1_.clear(); s2_.clear(); q1_.clear(); q2_.clear();
    r.work = d.__tmp__; d.__tmp__.clear();
    r.lib  = d.lib;

    if (!d.o.tsv.empty())
    {
        std::ofstream w;
        w.open(d.o.tsv);
        w << "NAME\tCOUNT" << std::endl;
        
        for (const auto &i : r.seqs)
        {
            w << i << "\t" << (r.sqc.count(i) ? r.sqc.at(i) : 0) << std::endl;
        }
        
        w.close();
    }
    
    return r;
};

KStats Kallisto(const FileName &f1, const FileName &f2, const KOptions &o)
{
    ProgramOptions opt;
    auto x = KInit(o);
    
    opt.d = (void *) (&x);
    opt.k = ::Kmer::k = o.k;
    opt.index = o.index;
    opt.threads = o.thr;
    
    opt.files.push_back(f1);
    opt.files.push_back(f2);

    // Required for feteching read names
    opt.fusion = true;
    
    KmerIndex index(opt);
    MinCollector collection(index, opt);
    ProcessReads(index, opt, collection);
    
    return KMerge(x);
}
