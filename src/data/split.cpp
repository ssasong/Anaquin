#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "data/split.hpp"
#include "tools/tools.hpp"
#include "data/reader.hpp"
#include "data/resources.hpp"
#include "tools/samtools.hpp"
#include "tools/calibrator.hpp"
#include "writers/fq_writer.hpp"
#include "writers/bam_writer.hpp"
#include "writers/sam_writer.hpp"
#include "parsers/parser_bam.hpp"

extern std::map<std::string, std::map<int, Anaquin::KMInfo>> __KMInfo__;

using namespace Anaquin;

// Defined in Kallisto
extern KStats Kallisto(const FileName &, const FileName &, const KOptions &);

template <typename T1, typename T2> void SKFreq(const T1 &x, T2 &m)
{
    for (const auto &i : x) // For all sequins
    {
        for (const auto &j : i.second) // For all k-mers in the sequin
        {
            assert(j.second);
            m[i.first].push_back(j.second);
        }
    }
};

std::string SAlignStats::SInsert(Bin b, Base min, Base max) const
{
    assert(ins.count(b));
    std::vector<double> x;
    
    for (auto i : ins.at(b))
    {
        if (i.first >= min && i.first <= max)
        {
            for (auto j = 0; j < ins.at(b).at(i.first); j++)
            {
                x.push_back(i.first);
            }
        }
    }
    
    return S2(SS::mean(x)) + " +- " + S2(SS::SD(x));
}

template <typename T> void SQuantKM(T &x)
{
    // Descriptive statistics for unique k-mers
    for (auto &i : x.s2u)
    {
        std::sort(i.second.begin(), i.second.end());
        x.d2u.mus[i.first]  = SS::mean(i.second);
        x.d2u.q25[i.first]  = quant(i.second, 0.25);
        x.d2u.q75[i.first]  = quant(i.second, 0.75);
        x.d2u.sds[i.first]  = SS::SD(i.second);
        x.d2u.mins[i.first] = *(i.second.begin());
        x.d2u.meds[i.first] = quant(i.second, 0.50);
        x.d2u.maxs[i.first] = *(i.second.rbegin());
    }
    
    // Descriptive statistics for shared k-mers
    for (auto &i : x.s2s)
    {
        std::sort(i.second.begin(), i.second.end());
        x.d2u.mus[i.first]  = SS::mean(i.second);
        x.d2s.q25[i.first]  = quant(i.second, 0.25);
        x.d2s.q75[i.first]  = quant(i.second, 0.75);
        x.d2s.sds[i.first]  = SS::SD(i.second);
        x.d2s.mins[i.first] = *(i.second.begin());
        x.d2s.meds[i.first] = quant(i.second, 0.50);
        x.d2s.maxs[i.first] = *(i.second.rbegin());
    }
}

void Anaquin::writeSTable_1(const FileName &src, const FileName &dst, const SOptions &o, Count nExp, Count nObs, Count nCV, const Label &expL, const Label &obsL, float d)
{
    const auto tmp = tmpFile();
    RMeanCV(src, tmp, expL, obsL, nExp, nObs, nCV, d);
    auto t = readFile(tmp);
         t = replace(t, "NAME",  "EXPECTED_ABUNDANCE");
         t = replace(t, "COUNT", "SEQUIN_DETECTED");
         t = replace(t, "MEAN",  "ABUNDANCE_(MEAN_READ_COUNT)");
    o.generate(dst);
    o.writer->open(dst);
    o.writer->write(t);
    o.writer->close();
}

void Anaquin::SWriteReads(Product mode, const FileName &file, const FQStats &stats, const SOptions &o)
{
    if (!o.debug)
    {
        return;
    }
    
    const auto format = "%1%\t%2%\t%3%\t%4%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "NAME" % "READ1" % "READ2" % "LABEL").str());
    
    auto write = [&](const std::map<ReadName, SequinID> &x, const std::map<ReadName, SequinID> &y, const Label &l)
    {
        for (const auto &i : x)
        {
            o.writer->write((boost::format(format) % i.first % i.second % y.at(i.first) % l).str());
        }
    };
    
    forBin([&](Bin c)
    {
        std::string s;
        
        switch (mode)
        {
            case Product::Genomics:
            {
                switch (c)
                {
                    case HP: { s = "HP";         break; }
                    case IF: { s = "Info";       break; }
                    case MT: { s = "Mito";       break; }
                    case HL: { s = "HLA";        break; }
                    case LD: { s = "Ladder";     break; }
                    case IM: { s = "Immune";     break; }
                    case ES: { s = "Sample";     break; }
                    case VC: { s = "Vector";     break; }
                    case SV: { s = "Structural"; break; }
                    case GR: { s = "Germline";   break; }
                    case SO: { s = "Somatic";    break; }
                    case MS: { s = "Micro";      break; }
                    case MI: { s = "MSI";        break; }
                    default: { break; }
                }
                
                break;
            }

            case Product::Meta:
            {
                switch (c)
                {
                    case IF: { s = "Info";    break; }
                    case LD: { s = "Ladder";  break; }
                    case ES: { s = "Sample";  break; }
                    case VC: { s = "Vector";  break; }
                    case GR: { s = "Sequins"; break; }
                    default: { break; }
                }

                break;
            }

            case Product::RNA:
            {
                switch (c)
                {
                    case IF: { s = "Info";    break; }
                    case ES: { s = "Sample";  break; }
                    case GR: { s = "Sequins"; break; }
                    default: { break; }
                }
                
                break;
            }

            default: { break; }
        }

        if (stats.K.r1.count(c) || stats.K.r2.count(c))
        {
            write(stats.K.r1.at(c), stats.K.r2.at(c), s);
        }
    });

    o.writer->close();
}

/*
 * Post-calculation from completed Kallisto statistics. This function should be called after Kallisto. It
 * is responsible for constructing various ladders.
 */

static void SPostCal(FQStats &stats, const SOptions &o)
{
    SKFreq(stats.K.uniqs, stats.R.s2u);
    SKFreq(stats.K.shared, stats.R.s2s);
    SQuantKM(stats.R);
     
    if (stats.R.d2s.mins.empty() && stats.R.d2u.mins.empty())
    {
        o.warn("Sequin not found. Please check your input files.");
    }
    
    std::set<Bin> only;
    
    switch (o.prod)
    {
        case Product::RNA:  { only = std::set<Bin> { ES, GR };         break; }
        case Product::Meta: { only = std::set<Bin> { ES, GR, LD, VC }; break; }
        case Product::Genomics:
        {
            only = std::set<Bin> { MT, MS, HL, HP, LD, SV, IM, ES, GR, SO, MI, VC };
            break;
        }
    }

    // Merge and generate FASTQ for each bin
    SMerge(stats, o, only);

    assert(stats.dil() >= 0 && stats.dil() <= 1.0);
}

struct AbstractMerge
{
    virtual void init(const SOptions &o, Bin m, const Label &x) = 0;
    virtual void close() = 0;
    virtual void merge(FQStats &, Bin, Thread, const SOptions &) = 0;
};

struct FQMerge : public AbstractMerge
{
    void init(const SOptions &o, Bin m, const Label &x) override
    {
        const auto f1_ = o.work + "/" + o.name + "_" + x + "_1.fq.gz";
        const auto f2_ = o.work + "/" + o.name + "_" + x + "_2.fq.gz";

        auto i1 = std::find_if(f1.begin(), f1.end(), [&](const std::pair<Bin, FileName> &i)
        {
            return i.second == f1_;
        });

        auto i2 = std::find_if(f2.begin(), f2.end(), [&](const std::pair<Bin, FileName> &i)
        {
            return i.second == f2_;
        });

        f1[m] = f1_;
        f2[m] = f2_;

        if (i1 != f1.end())
        {
            assert(i2 != f2.end());
            o1[m] = o1[i1->first];
            o2[m] = o2[i2->first];
            prim[m] = false;
        }
        else
        {
            prim[m] = true;
            o1[m] = std::shared_ptr<std::ofstream>(new std::ofstream(f1[m], std::ios::binary | std::ios::out));
            o2[m] = std::shared_ptr<std::ofstream>(new std::ofstream(f2[m], std::ios::binary | std::ios::out));
        }
    }
    
    void close() override
    {
        for (auto &i : o1) { if (prim[i.first]) { i.second->close(); } }
        for (auto &i : o2) { if (prim[i.first]) { i.second->close(); } }
    }

    void merge(FQStats &stats, Bin b, Thread i, const SOptions &o) override
    {
        assert(f1.count(b) && f2.count(b));
        assert(stats.K.f1[b][i].size() == stats.K.f2[b][i].size());
        
        auto myCopyGZ = [&](const FileName &src, std::shared_ptr<std::ofstream> w)
        {
            switch (copyGZ(src, w))
            {
                case CopyGZStatus::Failed:
                {
                    o.logWarn(src + " malformed");
                    return false;
                }
                    
                case CopyGZStatus::Corrected:
                {
                    o.logWarn(src + " malformed but corrected");
                    return true;
                }

                case CopyGZStatus::Success: { return true; }
            }
        };

        if (!myCopyGZ(stats.K.f1[b][i], o1[b])) { return; }
        if (!myCopyGZ(stats.K.f2[b][i], o2[b])) { return; }

        rm(stats.K.f1[b][i]); // Remove the old file
        rm(stats.K.f2[b][i]); // Remove the old file
    }
    
    std::map<Bin, bool> prim;
    std::map<Bin, FileName> f1, f2;
    std::map<Bin, std::shared_ptr<std::ofstream>> o1, o2;
};

void Anaquin::SMerge(FQStats &stats, const SOptions &o, const std::set<Bin> &only)
{
    assert(!o.bam);
    auto xs = std::vector<std::shared_ptr<AbstractMerge>>();
    xs.push_back(std::shared_ptr<AbstractMerge>(new FQMerge()));
    
    auto bin2Str = [&](Bin x)
    {
        switch (x)
        {
            case SV: { return "sv";     }
            case IF: { return "info";   }
            case HL: { return "hla";    }
            case MT: { return "mito";   }
            case LD: { return "ladder"; }
            case IM: { return "immune"; }
            case ES: { return "sample"; }
            case VC: { return "vector"; }
            case HP:
            case MS:
            case GR:
            case MI:
            case SO: { return "sequin"; }
        }
    };
    
    forBin([&](Bin x)
    {
        if (!only.count(x))
        {
            return;
        }
        else if (o.onlySeqLad && x != GR && x != LD)
        {
            return;
        }
        
        for (auto &i : xs)
        {
            i->init(o, x, bin2Str(x));
        }
    });
    
    // For each category...
    for (auto c : stats.K.f1)
    {
        auto bin = c.first;

        if (!only.count(bin))
        {
            continue;
        }
        else if (o.onlySeqLad && bin != GR && bin != LD)
        {
//            continue;
        }

        auto m = std::map<Bin, std::string>
        {
            { ES, "ES" },
            { IF, "IF" },
            { MT, "MT" },
            { MS, "MS" },
            { HL, "HL" },
            { HP, "HP" },
            { LD, "LD" },
            { SV, "SV" },
            { IM, "IM" },
            { ES, "ES" },
            { GR, "GR" },
            { SO, "SO" },
            { VC, "VC" },
        };
        
        // For each thread...
        for (auto i = 0u; i < stats.K.f1[bin].size(); i++)
        {
            o.info("Merging thread " + S0(i));

            for (auto &j : xs)
            {
                j->merge(stats, bin, i, o);
            }
        }
    }
    
    for (auto &i : xs) { i->close(); }
    
    assert(!stats.K.work.empty());
    
    // Remove the directory keeping Kallisto files
    removeD(stats.K.work);
}

SCStats Anaquin::SCalibrateP(const std::vector<Bin> &b, Proportion p, const FQStats &stats, const SOptions &o, const SCalibratePFiles &f)
{
    assert(!o.bam);
    SCStats x;
    
    // Number of sample paired reads (equal before and after calibration)
    x.bSam = x.aSam = stats.K.binN(ES);

    // Number of paired reads before calibration
    x.bSeq = 0; for (auto &i : b) { x.bSeq += stats.K.binN(i); }
    
    if (x.bSam == 0) { o.warn("Number of sample reads is zero. Scaling factor set to 1."); }
    if (x.bSeq == 0) { o.warn("Number of sequin reads is zero. Scaling factor set to 0."); }

    if (p == NO_CALIBRATION) // No calibration?
    {
        o.info("Calibration is NA. Skipped.");
        x.p  = NO_CALIBRATION;
        x.o1 = f.i1(o); // No calibration
        x.o2 = f.i2(o); // No calibration
        x.aSeq = x.bSeq;
        return x;
    }
    else if (p > 1.0) // Calibrate by reads? (p is the number of target paired reads)
    {
        // Convert to target paired reads (very important)
        //p = int(2 * p);

        // Number of target sequin reads after calibration
        x.tar = x.bSeq < p ? x.bSeq : p;

        // Scaling factor for calibration (0.5 because of paired-ends)
        x.p = x.bSeq ? 0.5 * x.tar / x.bSeq : 0.0;

        o.info("Calibration by reads: " + S0(p));
        o.info("Sequin: " + std::to_string(x.bSeq));
        o.info("Target: " + std::to_string(x.tar));
        o.info("Scaling: " + std::to_string(x.p));
    }
    else // Calibrate by percentage?
    {
        // Number of target sequin reads after calibration
        x.tar = (p / (1.0 - p)) * x.bSam;
        
        // Make sure the target doesn't goto zero if sample reads is non-zero
        if (x.bSam && !x.tar)
        {
            x.tar = x.bSam;
        }
        
        // Scaling factor for calibration
        x.p = (x.bSam == 0) ? 1.0 : (x.bSeq == 0) ? 0.0 : (x.tar >= x.bSeq ? 1.0 : ((float) x.tar) / x.bSeq);

        o.info("Calibration by percentage: " + std::to_string(p));
        o.info("Sample: " + std::to_string(x.bSam));
        o.info("Sequin: " + std::to_string(x.bSeq));
        o.info("Target: " + std::to_string(x.tar));
        o.info("Scaling: " + std::to_string(x.p));
    }
    
    assert(x.p >= 0.0 && x.p <= 1.0);
    auto r = SelectionCalibrator::createFQ(f.i1(o), f.i2(o), f.o1(o), f.o2(o))->calibrate(x.p, o);
    
    x.o1   = r.o1;
    x.o2   = r.o2;
    x.aSeq = r.nSel;
    assert(x.bSeq >= x.aSeq);

    x.bSam *= 2; // Convert from number of pair-ended
    x.aSam *= 2; // Convert from number of pair-ended
    x.bSeq *= 2; // Convert from number of pair-ended
    x.aSeq *= 2; // Convert from number of pair-ended

    return x;
}

void Anaquin::SWriteLadder(std::shared_ptr<Ladder> l3, const FileName &file, const FQStats &stats, const SOptions &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "NAME" % "STOICHOMETRY" % "COPY" % "Q50" % "READ").str());
    
    for (const auto &seq : stats.K.seqs)
    {
        if (isSubstr(seq, "LD_"))
        {
            auto write = [&](const SequinID &seq)
            {
                const auto u = stats.R.d2u.meds.count(seq) ? stats.R.d2u.meds.at(seq) : NAN;
                const auto m = u;
                const auto x = stats.K.sqc.count(seq) ? stats.K.sqc.at(seq) : 0;
                
                if (l3->contains(seq))
                {
                    o.writer->write((boost::format(format) % seq
                                                           % 1
                                                           % l3->input(seq)
                                                           % (std::isnan(m) ? MISSING : S2(m))
                                                           % x).str());
                }
                else
                {
                    o.writer->write((boost::format(format) % seq
                                                           % MISSING
                                                           % MISSING
                                                           % (std::isnan(m) ? MISSING : S2(m))
                                                           % x).str());
                }
            };
            
            write(seq);
        }
    }
    
    o.writer->close();
}

void Anaquin::SKallisto(FQStats &stats, const FileName &f1, const FileName &f2, const SOptions &o)
{
    auto o_ = o;
    o_.bam = false;
    o_.writeReads = o.debug;
    o_.tsv = stats.tsv = tmpFile(); // Always request for TSV
    
    // Running Kallisto
    stats.K = Kallisto(f1, f2, o_);
    
    // Post-calculation
    SPostCal(stats, o_);
}

std::string Anaquin::SVersion(const Reference &r, const KStats &x)
{
    std::map<SequinID, Count> c;
    
    for (const auto &i : x.sqc)
    {
        if (r.t1()->count(i.first))
        {
            c[i.first] = i.second;
        }
    }
    
    return !c.empty() ? r.t1()->at(max(c)) : MISSING;
}

void Anaquin::writeKmers(const FileName &file, const FQStats &stats, const SOptions &o)
{
    // Consturct a convenient structure for mapping
    std::map<std::string, std::map<Kmer, KMInfo *>> m;

    for (auto &i : __KMInfo__)
    {
        for (auto &j : i.second)
        {
            m[i.first][j.second.kmer] = &j.second;
        }
    }

    auto add = [&](const SequinKM &x)
    {
        for (const auto &i : x)
        {
            for (const auto &j : i.second)
            {
                const auto k = m[i.first].count(j.first) ? j.first : revcomp(j.first);
                assert(m[i.first][k]);
                m[i.first][k]->abund += j.second;
            }
        }
    };
    
    add(stats.K.uniqs);
    add(stats.K.shared);
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "NAME"
                                           % "TYPE"
                                           % "SEQUENCE"
                                           % "POSITION"
                                           % "COUNT").str());

    for (const auto &i : __KMInfo__)
    {
        const auto &seq = i.first;
        
        for (const auto &j : i.second)
        {
            if (isSubstr(seq, "LD_"))
            {
                const auto isUniq = stats.K.uniqs.count(seq) && stats.K.uniqs.at(seq).count(j.second.kmer);
                
                o.writer->write((boost::format(format) % seq
                                                       % (isUniq ? "Unique" : "Shared")
                                                       % j.second.kmer
                                                       % j.first
                                                       % j.second.abund).str());
            }
        }
    }

    o.writer->close();
}
