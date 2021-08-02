#include "sequins/sequins.hpp"
#include "tools/calibrator.hpp"
#include "sequins/meta/meta.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_csv.hpp"

using namespace Anaquin;

typedef LadderCalibrator::Method Method;
typedef LadderCalibrator::Options::PoolData PoolData;

static void writeL(const FileName &file, const DecoyAnalyzer::Results &r, const MDecoyOptions &o)
{
    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    const auto tmp = tmpFile();
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % "NAME"
                                      % "SEQUIN"
                                      % "STANDARD"
                                      % "STOICHOMETRY"
                                      % "COPY"
                                      % "READ"
                                      % "MIN"
                                      % "MEAN"
                                      % "MAX").str());
    
    auto r1 = o.r1.inters();
    
    for (auto &i : r1)
    {
        for (auto &j : i.second.data())
        {
            const auto &name = j.second.id();
            
            if (!isSubstr(name, "_LD"))
            {
                continue;
            }
            
            // Synthetic doesn't need any calibration. No need for r.after.
            const auto stats = r.decoy.r1.find(name)->stats();
            
            o.writer->write((boost::format(f) % name
                                              % noLast(name, "_")
                                              % MSeq2Std(name)
                                              % 1
                                              % S0(o.l3->input(noLast(name, "_")))
                                              % stats.n
                                              % stats.min
                                              % stats.mean
                                              % stats.max).str());
        }
    }
    
    o.writer->close();

    // Always generate the merged version
    mergeL(o.work + "/" + file, o.work + "/" + file);
}

void MDecoyResults::writeR(const MDecoyResults &r, const MDecoyOptions &o)
{
    if (o.tsvR.empty())
    {
        return;
    }
    
    const auto f = boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%");
    
    o.generate(o.tsvR);
    o.writer->open(o.tsvR);
    o.writer->write((boost::format(f) % "NAME"
                                      % "MIX"
                                      % "CHROM"
                                      % "START"
                                      % "END"
                                      % "EDGE"
                                      % "MEAN_READ_LENGTH"
                                      % "READ"
                                      % "COVERAGE").str());

    auto r1 = o.r1.inters();

    for (auto &i : r1)
    {
        for (auto &j : i.second.data())
        {
            const auto &name = j.second.id();
            const auto stats = r.decoy.r1.find(name)->stats();
            
            const auto l1 = o.l1->findBySub(noPID(name));
            o.writer->write((boost::format(f) % name
                                              % (std::isnan(l1) ? MISSING : S6(l1))
                                              % i.first
                                              % j.second.l().start
                                              % j.second.l().end
                                              % o.edge
                                              % r.lib.meanRL()
                                              % stats.n
                                              % stats.mean).str());
        }
    }

    o.writer->close();
}

// Used in reportE()
struct SearchResults
{
    SearchResults & operator+=(const SearchResults &x)
    {
        this->ts += x.ts;
        this->tv += x.tv;
        this->in += x.in;
        this->dl += x.dl;
        return *this;
    }
    
    inline Count total() const { return ts + tv + in + dl; }
    
    Count ts = 0; // Transition
    Count tv = 0; // Transversion
    Count in = 0; // Insertion
    Count dl = 0; // Deletion
};

MDecoyResults::ErrorReport MDecoyResults::reportE(const FileName &tsvE, const FileName &tsvR)
{
    MDecoyResults::ErrorReport r;

    auto tmp1 = tmpFile();
    auto tmp2 = tmpFile();
    auto tmp3 = tmpFile();

    // Keep only mismatches
    RGrep(tsvE, tmp1, "LABEL", "Match", false);
    
    // Only calibrated (that's what sequins can inform about the sample)
    RGrep(tmp1, tmp2, "CALIBRATED", "true");

    // Only chrQB is being used for error calculation
    RGrep(tsvR, tmp3, "CHROM", MDecoyChrQB);

    struct Data
    {
        SequinID name;
        Label label;
        std::string data1, data2, data3;
    };
    
    std::map<Base, Data> m;
    
    ParserCSV::parse(Reader(tmp2), [&](const ParserCSV::Data &x, Progress i)
    {
        if (!i) { return; }
        
        assert(x[0] == "true" || x[0] == "false");
        assert(x[1] == MDecoyChrQB);
        assert(x[4] == "All");

        Data d;
        d.name  = x[2];
        d.label = x[3];
        d.data1 = x[5];
        d.data2 = x[6];
        d.data3 = x[7];

        m[stod(x[5])] = d;
    }, "\t");

    // Mean read length
    const auto rl = RMean(tmp3, "MEAN_READ_LENGTH");

    // Total number of reads in all chrQB regions
    const auto allReads = RSum(tmp3, "READ");

    // Number of errors, number of alignments and length for a stratification
    auto f1 = [&](Count x, Count rl, Base reads, Base size)
    {
        // (Number of errors / number of reads) / (mean read length)
        const auto base = S6(((Proportion) x / reads) / rl);

        // COUNT ; PER_BASE ; REGION_SIZE
        return S0(x) + " ; " + base + " ; " + S4(size);
    };
    
    // Count SNP, insertions and deletions from tsvE
    auto f2 = [&](Base start, Base end)
    {
        SearchResults r;
     
        for (auto i = start; i < end; i++)
        {
            if (m.count(i))
            {
                if (m[i].label == "SNP")
                {
                    const auto tv = std::set<std::string> { "AG", "GA", "CT", "TC" };

                    // Transversion?
                    if (tv.count(m[i].data2))
                    {
                        r.tv += stod(m[i].data3);
                    }
                    else
                    {
                        r.ts += stod(m[i].data3);
                    }
                }
                else if (m[i].label == "Insertion")
                {
                    r.in += stod(m[i].data3);
                }
                else
                {
                    r.dl += stod(m[i].data3);
                }
            }
        }
        
        return r;
    };

    // Total size in chrQB
    auto size = (Base) sum(RSubtract(tmp3, "START", "END"));
    
    // Edge for each region in chrQB
    const auto edge = numeric<Base>(RReadTSV(tmp3, "EDGE"));

    // Assume edges are all equal (they must)
    size -= edge.size() * (2 * edge.front());

    // Total number of sequencing errors on calibrated sample
    r.text[0] = f1(RCount(tmp2, "CHROM", MDecoyChrQB), (Count) rl, allReads, size);
    
    // Total number of mismatches on calibrated sample
    r.text[1] = f1(RCount(tmp2, "LABEL", "SNP"), (Count) rl, allReads, size);

    // Total number of insertions on calibrated sample
    r.text[3] = f1(RCount(tmp2, "LABEL", "Insertion"), (Count) rl, allReads, size);
    
    // Total number of deletions on calibrated sample
    r.text[4] = f1(RCount(tmp2, "LABEL", "Deletion"), (Count) rl, allReads, size);

    if (m.empty())
    {
        r.text[2] = MISSING;
    }
    else
    {
        // Break down SNP into transversions and transitions
        const auto x = f2(m.begin()->first, m.rbegin()->first);
        
        // Transition to Transversion ratio
        r.text[2] = S4((Proportion) x.ts / x.tv);
    }

    return r;
}

MDecoyResults Anaquin::MDecoyAnalysis(const FileName &src, const MDecoyOptions &o_)
{
    MDecoyResults r;
    
    auto o = o_;
    const auto writeD = !o.writeD.empty() ? o.writeD : tmpFile(); // Always required
    const auto writeT = !o.writeT.empty() ? o.writeT : tmpFile(); // Always required
    
    ParserFA::parse(Reader(o.index), [&](const ParserFA::Data &x)
    {
        o.seqs[x.id] = x.seq;
    });
    
    logRun([&]()
    {
        // Never ask DecoyAnalyzer to calibrate because it'd mess up with ladder reads
        auto tmp = o; tmp.writeD = writeD; tmp.writeT = writeT;
        
        auto sq = std::shared_ptr<BAMWriter>(new BAMWriter()); sq->open(o.originalW + "/meta_sequin.bam");
        auto ld = std::shared_ptr<BAMWriter>(new BAMWriter()); ld->open(o.originalW + "/meta_ladder.bam");
        auto vc = std::shared_ptr<BAMWriter>(new BAMWriter()); vc->open(o.originalW + "/meta_vector.bam");

        auto pc_ = ParitionCount();
        
        r.B1 = DecoyAnalyzer::analyze(src, "", tmp, [&](ParserBAM::Data &x, const DInter *r1, const DInter *r2, bool trimmed)
        {
            // Skip the read if trimmed (Tim's request to make trimmed going to VC)
            if (trimmed) { vc->write(x); pc_[Bin::VC]++; return; }
            
            if (!isDecoy(x.cID))      { pc_[Bin::ES]++; return; }
            if (x.cID == MDecoyChrQL) { ld->write(x); pc_[Bin::LD]++; return; }
            
            if (r2)
            {
                Bin b;
                
                switch (b = MBin(r2->name()))
                {
                    case IF:
                    case VC: { vc->write(x); break; }
                    default: { sq->write(x); break; }
                }
                
                pc_[b]++;
            }
            else
            {
                sq->write(x);
                pc_[Bin::GR]++;
            }
        }, [&]() { sq->close(); ld->close(); vc->close(); }); r.B1._pc = pc_;
        
        // Check calibrated BAM?
        if (o.seqC != NO_CALIBRATION)
        {
            auto _o = cloneO(o);
            
            _o.ladC   = NO_CALIBRATION;
            _o.seqC   = NO_CALIBRATION;
            _o.writeC = ""; // Only statistics
            _o.writeS = ""; // Only statistics
            _o.writeD = ""; // Only statistics
            _o.writeT = ""; // Only statistics
            
            r.B2 = DecoyAnalyzer::analyze(o_.originalW + "/meta_calibrated.bam", "", _o, [&](ParserBAM::Data &x, const DInter *r1, const DInter *r2, bool) {});
        }
    }, o, "Started initial analysis", "Completed initial analysis");

    r.lib   = r.B1.bamLib;
    r.samp  = r.B1.samp;
    r.decoy = r.B1.decoy;
    
    // Ladder before ladder calibration
    writeL(o.tsvL1, r.B1, o);
    
    /*
     * Ladder calibration. Sampling coverage unavailable in metagenomics, thus Method::Pooling is the
     * only option.
     */
    
    FileName tmp;

    // Metagenomics doesn't have mirror regions or sampling coverage
    if (o.ladC != NO_CALIBRATION)
    {
        tmp = o.writeL; // Calibrated ladder alignment
        assert(!tmp.empty());
        
        if (o.ladC <= 1.0)
        {
            LadderCalibrator::Options o_;
            
            o_.d2 = std::shared_ptr<PoolData>(new PoolData());
            o_.d2->p = o.ladC; // How much to calibrate
            o_.d2->tsvL = o.work + "/" + o.tsvL1;
            
            o_.r1 = o.r1;
            o_.meth = Method::Pooling;
            o_.logger = o.logger;
            o_.output = o.output;
            o_.writer = o.writer;
            o_.merged = true;

            logRun([&]()
            {
                LadderCalibrator::createBAM(writeT, tmp)->calibrate(o_);
            }, o, "Started synthetic ladder calibration", "Completed synthetic ladder calibration");
        }
        else
        {
            const auto bSeq = r.B1._pc[Bin::LD];
            
            // Number of target sequin reads after calibration
            const auto tar = bSeq < o.ladC ? bSeq : o.ladC;

            // Scaling factor for calibration (0.5 because of paired-ends)
            const auto p = bSeq ? (Proportion) tar / bSeq : 0.0;

            o.logInfo("Ladder reads before calibration: " + std::to_string(bSeq));
            o.logInfo("Target: " + std::to_string(tar));
            o.logInfo("Scaling factor: " + std::to_string(p));

            // Calibrate and return the number of reads after calibration
            SelectionCalibrator::createBAM(writeT, tmp)->calibrate(p, o);
        }
        
        auto _o = cloneO(o);
        
        _o.writeC = ""; // Only statistics
        _o.writeS = ""; // Only statistics
        _o.writeD = ""; // Only statistics
        _o.writeT = ""; // Only statistics

        // Analysis on post ladder calibration
        r.B3 = DecoyAnalyzer::analyze(tmp, "", _o, [&](ParserBAM::Data &x, const DInter *r1, const DInter *r2, bool) {});
        
        writeL(o.tsvL2, r.B3, o); // After synthetic normalisation
    }
    else
    {
        // No ladder calibration, the output file is also the input file
        tmp = writeT;
    }
    
    //GDecoyResults::writeF(r, o);
    
    DecoyAnalyzer::ErrorOptions eo(o.tsvE, o.r1, o.r2, MDecoyChrQB, o.v1, true, o);
    eo.data[&r.decoy] = false; eo.data[&r.decoy] = true;
    
    // Generating error profiles
    DecoyAnalyzer::writeE(eo);
    
    MDecoyResults::writeR(r, o);

    // Merging possibly calibrated sequin reads?
    if (r.B1.wM)
    {
        logRun([&]()
        {
            assert(r.B1.wM);
            assert(exists(writeT));
            
            ParserBAM::parse(writeT, [&](ParserBAM::Data &x, const ParserBAM::Info &)
            {
                r.B1.wM->write(x);
            });

            r.B1.wM->close();
        }, o, "Started merging sample and calibrated sequin reads",
              "Completed merging sample and calibrated sequin reads");
    }
    
    return r;
}
