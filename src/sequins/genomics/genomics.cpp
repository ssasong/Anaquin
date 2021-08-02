#include "tools/tools.hpp"
#include "tools/random.hpp"
#include "data/standard.hpp"
#include "data/resources.hpp"
#include "sequins/sequins.hpp"
#include "tools/calibrator.hpp"
#include "parsers/parser_fa.hpp"
#include "writers/vcf_writer.hpp"
#include "writers/bam_writer.hpp"
#include "writers/file_writer.hpp"
#include "parsers/parser_bambed.hpp"
#include "sequins/genomics/genomics.hpp"

using namespace Anaquin;
using namespace std::chrono;

typedef DecoyAnalyzer::Results::MData MData;
typedef DecoyAnalyzer::IDData IDData;
typedef DecoyAnalyzer::Results DResults;
typedef ParserBAMBED::Response Response;
typedef LadderCalibrator::Method Method;
typedef LadderCalibrator::Options::PoolData PoolData;

/*
 * Capture regions may break a sequin region into non-overlapping smaller regions. This function merges them and return
 * a data structure that resembles as if a single region.
 */

static std::shared_ptr<DInter::Stats> mergeStats(const Chr2DInters &r, const std::string &x, const WriterOptions &o)
{
    extern bool __HACK_IS_CAPTURE__;
    auto stats = std::shared_ptr<DInter::Stats>(new DInter::Stats());
    
    // Only if capture and original sequin name failed to match
    if (__HACK_IS_CAPTURE__ && !r.find(x))
    {
        std::vector<Coverage> p25, p50, p75, mean, n;

        for (auto i = 1;; i++) // Start from "1"
        {
            const auto name = x + "_" + std::to_string(i); // Capture
            const DInter *m = nullptr;
            
            if ((m = r.find(name)))
            {
                const auto tmp = m->stats();
                p25.push_back(tmp.p25);
                p50.push_back(tmp.p50);
                p75.push_back(tmp.p75);
                n.push_back(tmp.n);
                mean.push_back(tmp.mean);
                continue;
            }
            
            break;
        }
        
        if (p25.empty())
        {
            return nullptr;
        }
        
        stats->n    = SS::sum(n);
        stats->p25  = SS::mean(p25);
        stats->p50  = SS::mean(p50);
        stats->p75  = SS::mean(p75);
        stats->mean = SS::mean(mean);
        
        // Only a few selected fields should be sufficient for coverage calculation
        return stats;
    }
    else
    {
        if (r.find(x))
        {
            *stats = r.find(x)->stats();
            return stats;
        }
        
        return nullptr;
    }
}

void GDecoyResults::writeR(const GDecoyResults &r, const GDecoyOptions &o)
{
    if (o.tsvR.empty())
    {
        return;
    }
    
    const auto f = boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10$.2f\t%11$.2f\t%12$.2f\t%13$.2f");
    
    o.generate(o.tsvR);
    o.writer->open(o.tsvR);
    o.writer->write((boost::format(f) % "NAME"
                                      % "CHROM"
                                      % "START"
                                      % "END"
                                      % "EDGE"
                                      % "MEAN_READ_LENGTH"
                                      % "SAMPLE_READ"
                                      % "PRE_READ"
                                      % "POST_READ"
                                      % "SAMPLE_COVERAGE"
                                      % "PRE_COVERAGE"
                                      % "POST_COVERAGE"
                                      % "SCALE").str());

    o.writer->write((boost::format(f) % "All"
                                      % MISSING
                                      % MISSING
                                      % MISSING
                                      % MISSING
                                      % r.lib.meanRL()
                                      % r.samp.total()
                                      % r.before.total()
                                      % r.after.total()
                                      % MISSING
                                      % MISSING
                                      % MISSING
                                      % MISSING).str());

    auto r1 = o.r1.inters();

    for (auto &i : r1)
    {
        for (auto &j : i.second.data())
        {
            const auto &name = j.second.id();
            const auto norm = o.seqC == NO_CALIBRATION ? getNAFromD(r.c1.norms, name, 4) : MISSING;
            
            const auto afterN = r.after.r1.find(name) ? S0(r.after.r1.find(name)->stats().n) : MISSING;
            
            // Sample coverage before and after calibration (r2 to make it consistent)
            const auto sampS = mergeStats(r.samp.r2, name, o);
            
            // Sequin coverage before calibration (merging because it could be intersecting capture regions)
            const auto beforeS = mergeStats(r.before.r2, name, o);

            // Sequin coverage after calibration (merging because it could be intersecting capture regions)
            const auto afterS = mergeStats(r.after.r2, name, o);
            
            o.writer->write((boost::format(f) % name
                                              % i.first
                                              % j.second.l().start
                                              % j.second.l().end
                                              % o.edge
                                              % r.lib.meanRL()
                                              % (sampS ? S0(sampS->n) : MISSING)
                                              % r.before.r1.find(name)->stats().n
                                              % afterN
                                              % (sampS ? S4(sampS->mean) : MISSING)
                                              % (beforeS ? S4(beforeS->mean) : MISSING)
                                              % (afterS  ? S4(afterS->mean)  : MISSING)
                                              % norm).str());
        }
    }

    o.writer->close();
}

void GDecoyResults::writeF(const GDecoyResults &r, const GDecoyOptions &o)
{
    if (o.tsvF.empty())
    {
        return;
    }
    
    const auto hasA = !o.writeMC.empty(); // Calibration?
    const auto f = boost::format("%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7$.2f\t%8$.2f");
    
    o.generate(o.tsvF);
    o.writer->open(o.tsvF);
    o.writer->write((boost::format(f) % "NAME"
                                      % "CHROM"
                                      % "START"
                                      % "END"
                                      % "PRE_READ"
                                      % "POST_READ"
                                      % "PRE_COVERAGE"
                                      % "POST_COVERAGE").str());
    
    auto r5 = o.attr.inters();
    
    for (auto &i : r5)
    {
        for (auto &j : i.second.data())
        {
            const auto &name = j.second.id();
            const auto  samp = r.samp.r1.find(name);
            
            o.writer->write((boost::format(f) % name
                                              % i.first
                                              % j.second.l().start
                                              % j.second.l().end
                                              % r.before.r5.find(name)->stats().n
                                              % (hasA ? S0(r.after.r5.find(name)->stats().n) : MISSING)
                                              % S4(r.before.r5.find(name)->stats().mean)
                                              % (hasA ? S4(r.after.r5.find(name)->stats().mean) : MISSING)).str());
        }
    }
    
    o.writer->close();
}

void GDecoyResults::writeL(const FileName &tsvL, const DecoyAnalyzer::Results &r, const DecoyAnalyzer::Options &o)
{
    if (tsvL.empty())
    {
        return;
    }

    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%";
    
    o.generate(tsvL);
    o.writer->open(tsvL);
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
                                              % GSeq2Std(noLast(name, "_"))
                                              % 1
                                              % S0(o.l3->input(noLast(name, "_")))
                                              % stats.n
                                              % stats.min
                                              % stats.mean
                                              % stats.max).str());
        }
    }

    o.writer->close();
}

void GDecoyResults::writeV(const GDecoyResults &r, const GDecoyOptions &o)
{
    if (o.tsvA.empty())
    {
        return;
    }

    o.generate(o.tsvA);
    o.writer->open(o.tsvA);

    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";
    o.writer->write("NAME\tLABEL\tCHROM\tPOSITION\tEXP_FREQ\tOBS_FREQ_CALIB\tREF_COUNT_CALIB\tVAR_COUNT_CALIB\tOBS_FREQ\tREF_COUNT\tVAR_COUNT");

    for (auto &i : o.v1->data.vars())
    {
        if (i.gt == Genotype::MSI)
        {
            continue;
        }
        else if (r.bAF.find(noPID(i.name)))
        {
            const auto X1 = r.bAF.find(noPID(i.name));
            const auto X2 = r.aAF.find(noPID(i.name));

            const auto &R1 = X1 ? X1->R : 0; // Reference counts before calibration
            const auto &V1 = X1 ? X1->V : 0; // Reference counts before calibration
            const auto &R2 = X2 ? X2->R : 0; // Variant counts after calibration
            const auto &V2 = X2 ? X2->V : 0; // Variant counts after calibration

            assert(R1 >= 0);
            assert(R2 >= 0);
            assert(V1 >= 0);
            assert(V2 >= 0);

            o.writer->write((boost::format(f) % i.name
                                              % bin2Label(GBin(i.name))
                                              % i.cID
                                              % i.l.start
                                              % S6(o.l1->input(i.name))
                                              % replaceNA((Proportion) V2 / (R2+V2), 6)
                                              % R2
                                              % V2
                                              % replaceNA((Proportion) V1 / (R1+V1), 6)
                                              % R1
                                              % V1).str());
        }
    }
    
    o.writer->close();
}

static void calibrateBAM(const FileName &src,
                         const FileName &dst,
                         GDecoyResults::SingleCalibrate &c,
                         const GDecoyOptions &o)
{
    assert(c.p >= 0.0 && c.p <= 1.0);
    c.after = 0;
    auto r = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - c.p));
    
    BAMWriter w;
    w.open(dst);
    
    ParserBAM::parse(src, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
    {
        if (i.p && !(i.p % 1000000)) { o.logInfo(S0(i.p)); }
        if (i.p == 0) { w.writeH(x); }
        x.lName();

        if (r->select(x.name))
        {
            c.after++;
            w.write(x);
            return Response::OK;
        }
        
        return Response::SKIP_EVERYTHING;
    });
    
    w.close();
}

static void calibrateBAM(const FileName &src,
                         const FileName &dst,
                         const Norms &norms,
                         const GDecoyOptions &o)
{
    typedef std::map<SequinID, std::shared_ptr<RandomSelection>> Selection;
    
    Selection select;
    
    auto r1 = o.r1.inters();
    auto r2 = o.r2.inters();
    
    for (const auto &i : norms)
    {
        // Guard against no sequin coverage
        const auto p = std::isnan(i.second) ? 1.0 : i.second;
        
        assert(p >= 0 && p <= 1.0 && !std::isnan(p));
        
        // Create independent random generator for each calibration region
        select[i.first] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - p));
    }
    
    assert(select.size() == norms.size());
    
    BAMWriter w;
    w.open(dst);

    ParserBAMBED::parse(src, r1, [&](ParserBAM::Data &x, const ParserBAM::Info &i, const DInter *)
    {
        if (i.p && !(i.p % 1000000)) { o.logInfo(S0(i.p)); }
        if (i.p == 0) { w.writeH(x); }
        
        x.lName();
        
        auto mustKept = !x.mapped || x.isDuplicate;
        
        if (!mustKept)
        {
            DInter *inter = nullptr;
            
            /*
             * Should that be contains or overlap? We prefer overlaps because any read that is overlapped
             * into the regions still give valuable information and sequencing depth.
             */
            
            if (r1.count(x.cID) && (inter = r1.at(x.cID).overlap(x.l)))
            {
                // Calibrate both paired-ends because they have the same name
                if (!select.count(inter->name()) || select.at(inter->name())->select(x.name))
                {
                    mustKept = true;
                }
            }
            else
            {
                // Never throw away reads outside the regions
                mustKept = true;
            }
        }
        
        if (mustKept)
        {
            w.write(x);
            
            if (r2.count(x.cID) && r2.at(x.cID).overlap(x.l) && !x.isDuplicate)
            {
                r2.at(x.cID).overlap(x.l)->map(x.l);
            }
            
            return Response::OK;
        }
        
        return Response::SKIP_EVERYTHING;
    });
    
    w.close();
}

std::map<SequinID, GDecoyResults::VariantData> GDecoyResults::buildAF(const DecoyAnalyzer::Results::Metrics &x, bool isChrQ)
{
    const auto &r = Standard::instance().gen;
    std::map<SequinID, GDecoyResults::VariantData> mm;
    
    auto r2 = r.r2()->inters();
    
    for (auto &v : (isChrQ ? r.v4() : r.v1())->data.vars())
    {
        const auto std = GSeq2Std(v.name);
        
        // Relative to chrQ or hg38
        auto start = v.l.start - 1; // 0-based
        
        // Convert from hg38 to sequin relative
        if (!isChrQ)
        {
            const auto tmp = r2.find(noLast(noPID(v.name), "_"), false);
            
            // Unable to construct AF, for example, custom BED file
            if (!tmp)
            {
                continue;
            }
            
            start -= (tmp->l().start - 1);
        }

        switch (v.type())
        {
            case Variation::SNP:
            {
                auto f1 = [&](const MData &x)
                {
                    if (!x.count(std) || !x.at(std).count(start))
                    {
                        return (Count) 0;
                    }
                    
                    return x.at(std).at(start);
                };
                
                auto f2 = [&](const DecoyAnalyzer::SData &x)
                {
                    if (!x.count(std) || !x.at(std).count(start) || !x.at(std).at(start).count(v.snpType()))
                    {
                        return (Count) 0;
                    }
                    
                    return x.at(std).at(start).at(v.snpType());
                };
                
                if (x.sd.count("All"))
                {
                    mm[v.name].R = f1(x.sd.at("All").match);
                    mm[v.name].V = f2(x.sd.at("All").snps);
                    assert(mm[v.name].R >= 0);
                    assert(mm[v.name].V >= 0);
                }

                break;
            }
                
            case Variation::Deletion:
            case Variation::Insertion:
            {
                auto f1 = [&](const MData &x, Base base)
                {
                    if (!x.count(std) || !x.at(std).count(base))
                    {
                        return (Count) 0;
                    }
                    
                    return x.at(std).at(base);
                };
                
                auto f2 = [&](const IDData &x, Base base)
                {
                    if (!x.count(std) || !x.at(std).count(base))
                    {
                        return (Count) 0;
                    }
                    
                    return sum(x.at(std).at(base));
                };
                
                // AF for insertion
                auto I = [&](const MData &d1, const IDData &d2)
                {
                    // Number of insertions is simply the number reported in CIGAR (IGV has the information)
                    const auto V = f2(d2, v.l.start);

                    /*
                     * Inserted sequences are not in the reference, but we can infer the
                     * number of references by subtracing all reads by inserted reads. Note
                     * each inserted read is also a contribution to the total reads.
                     */
                    
                    const auto R1 = f1(d1, v.l.start);
                    const auto R2 = R1 - V;
                    
                    assert(R2 >= 0);
                    return std::pair<Coverage, Coverage>(f1(d1, v.l.start) - V, V);
                };

                // AF for Deletion
                auto D = [&](const MData &d1, const IDData &d2)
                {
                    // Number of deletions is simply the number reported in CIGAR
                    const auto V = f2(d2, v.l.start);

                    /*
                       * When we say deletion, we mean something like "TGTC" to "T".
                       * The actual sequences being deleted is "G" not "T".
                       */
                    
                    const auto R = f1(d1, v.l.start + 1);
                    
                    return std::pair<Coverage, Coverage>(R, V);
                };

                if (x.sd.count("All"))
                {
                    const auto b = v.type() == Variation::Deletion ? D(x.sd.at("All").match,
                                                                       x.sd.at("All").dls) :
                                                                     I(x.sd.at("All").match,
                                                                       x.sd.at("All").ins);
                    
                    mm[v.name].R = b.first;
                    mm[v.name].V = b.second;
                    
                    assert(mm[v.name].R >= 0);
                    assert(mm[v.name].V >= 0);
                }

                break;
            }
                
            default:
            {
                throw std::runtime_error("Only short mutations are supported");
            }
        }
    }
    
    return mm;
}

template <typename T> Coverage getLocalCoverage(const T &x, CalibrateMethod m)
{
    switch (m)
    {
        case CalibrateMethod::Mean:
        case CalibrateMethod::Percent: { return x.mean; }
        case CalibrateMethod::Median:  { return x.p50; }
        default: { throw std::runtime_error("Only local coverage methods supported"); }
    }
}

Coverage Anaquin::meanSamp(const FileName &file, const Label &key)
{
    const auto tmp = tmpFile();
    RGrep(file, tmp, "CHROM", key);
    return RMean(tmp, "SAMPLE_COVERAGE");
}

GDecoyOptions GDecoyOptions::create(const FileName &index, const AnalyzerOptions &o)
{
    GDecoyOptions tmp;

    if (o.debug) { tmp.work = o.work;    } // Always write to working directory if debug
    else         { tmp.work = tmpPath(); }

    tmp.edge = 550;
    tmp.showGen = false;
    tmp.logger  = o.logger;
    tmp.output  = o.output;
    tmp.writer = std::shared_ptr<Writer<>>(new FileWriter(tmp.work));

    const auto &r = Standard::instance().gen;
    tmp.h1 = *(r.r1()); tmp.h2 = *(r.r3());
    tmp.r1 = *(r.r2()); tmp.r2 = *(r.r4());
    tmp.attr = *(r.r5());

    tmp.l1 = r.l1();
    tmp.l3 = r.l3();
    tmp.v1 = r.v4();
    tmp.a3 = r.a3();

    tmp.index = index;
    tmp.errors.insert(GDecoyChrQS);

    return tmp;
}

GDecoyOptions GDecoyOptions::create(const FileName &writeS,
                                    const FileName &writeD,
                                    const FileName &writeM,
                                    const FileName &writeT,
                                    const FileName &writeL1,
                                    const FileName &writeL2,
                                    const FileName &calibF,
                                    const FileName &tsvE,
                                    const FileName &tsvR,
                                    const FileName &tsvF,
                                    const FileName &tsvA,
                                    const FileName &tsvL1,
                                    const FileName &tsvL2,
                                    const FileName &index,
                                    Proportion seqC,
                                    Proportion ladC,
                                    const AnalyzerOptions &o)
{
    GDecoyOptions o2;

    if (o.debug) { o2.work = o.work;    } // Always write to working directory if debug
    else         { o2.work = tmpPath(); }

    o2.seqC    = seqC;
    o2.ladC    = ladC;
    o2.writer  = o.writer;
    o2.logger  = o.logger;
    o2.output  = o.output;
    o2.showGen = false;

    o2.writeS  = writeS;
    o2.writeD  = writeD;
    o2.writeM  = writeM;
    o2.writeT  = writeT;
    o2.writeMC = calibF;
    o2.writeL1 = writeL1;
    o2.writeL2 = writeL2;

    o2.tsvE  = tsvE;
    o2.tsvR  = tsvR;
    o2.tsvF  = tsvF;
    o2.tsvA  = tsvA;
    o2.tsvL1 = tsvL1;
    o2.tsvL2 = tsvL2;

    const auto &r = Standard::instance().gen;

    o2.l1   = r.l1();
    o2.l3   = r.l3();
    o2.v1   = r.v4();
    o2.a3   = r.a3();
    o2.edge = 550;

    o2.h1 = *(r.r1()); o2.h2 = *(r.r3());
    o2.r1 = *(r.r2()); o2.r2 = *(r.r4());
    o2.attr = *(r.r5());

    o2.index = index;
    o2.writer = std::shared_ptr<Writer<>>(new FileWriter(o2.work));
    o2.errors.insert(GDecoyChrQS);

    return o2;
}

GDecoyResults Anaquin::GDecoyAnalysis(const FileName &f1, const FileName &f2, const GDecoyOptions &_o_, G1 g1, G2 g2)
{
    assert(!_o_.r1.inters().empty());
    assert(!_o_.r2.inters().empty());
    
    assert(!_o_.writeL1.empty() && !_o_.inputL2.empty());
    
    // Either decoy calibration or sample/sequin calibration or nothing
    assert(!(_o_.seqC != NO_CALIBRATION && !_o_.writeMC.empty()));

    GDecoyResults r;
    
    const auto twoBAM = !f2.empty();
    
    auto o = _o_; o.showGen = _o_.showGen;
    const auto isChrQ = f2.empty();
    const auto writeD = !o.writeD.empty() ? o.writeD : tmpFile(); // Always required (all decoy reads)
    const auto writeT = !o.writeT.empty() ? o.writeT : tmpFile(); // Always required (only trimmed decoy reads)
    
    if (o.seqs.empty())
    {
        ParserFA::parse(Reader(o.index), [&](const ParserFA::Data &x)
        {
            auto seq = x.seq;
            
            if (!isChrQ)
            {
                switch (GBin(x.id))
                {
                    case Bin::GR:
                    case Bin::SO:
                    {
                        seq = reverse(seq); // Sequins on hg19/hg38 are reversed
                        break;
                    }
                        
                    default: { break; }
                }
            }
            
            o.seqs[x.id] = seq;
        });
    }
    
    Count nLad = 0;
    auto ld = std::shared_ptr<BAMWriter>(new BAMWriter()); ld->open(_o_.writeL1);

    assert(!o.seqs.empty());
    logRun([&]()
    {
        DecoyAnalyzer::Options tmp = o; tmp.writeD = writeD; tmp.writeT = writeT;
        assert(tmp.seqC == o.seqC);
        
        r.B1 = DecoyAnalyzer::analyze(f1, f2, tmp, [&](const ParserBAM::Data &x,
                                                       const DInter *r1,
                                                       const DInter *r2,
                                                       bool trimmed)
        {
            g1(GDecoyStatus::Before, x, r1, r2, trimmed);
            if (x.cID == GDecoyChrQL && !trimmed) { nLad++; ld->write(x); }
        }); ld->close(); g2(GDecoyStatus::Before);
    }, o, "Started initial analysis", "Completed initial analysis");

    r.lib    = r.B1.bamLib; // Library diagnostic
    r.samp   = r.B1.samp;   // Sampling diagnostic
    r.before = r.B1.decoy;  // Decoy diagnostic

    // Ladder before calibration
    GDecoyResults::writeL(o.tsvL1, r.B1, o);
    
    o.tsvL1NotMerged = "notMerged.tsv";
    mv(o.work + "/" + o.tsvL1, o.work + "/" + o.tsvL1NotMerged);
    
    // TODO: ???? o.tsvL2 = "notMerged.tsv";
    // TODO: ???? mv(o.work + "/" + o.tsvL1, o.work + "/" + o.tsvL2);
    
    // Note writeL() writes out to like "S0015_LD_016_D_1", "S0015_LD_016_D_2". Let's combine.
    mergeL(o.work + "/" + o.tsvL1NotMerged, o.work + "/" + o.tsvL1);
    
    // We'll do again, but this is needed for calculating sample coverage
    GDecoyResults::writeR(r, o);
    
    /*
     * Ladder calibration
     */
    
    typedef LadderCalibrator::Options::SampleCNV2Data SampleCNV2Data;
    LadderCalibrator::Options o_;
    
    o_.d1 = std::shared_ptr<SampleCNV2Data>(new SampleCNV2Data());
    o_.r1 = o.r1;
    o_.d1->tsvL = o.work + "/" + o.tsvL1NotMerged; // Unmerged for calibration
    
    // Making it "chrQ" will also read in "chrQL". Let's assume "chrQR" is not causing problems here...
    o_.d1->sampC = meanSamp(o.work + "/" + o.tsvR, f2.empty() ? "chrQS" : "chr");

    o_.logger = o.logger;
    o_.output = o.output;
    o_.writer = o.writer;
    o.info("Sampling coverage for ladder calibration: " + std::to_string(o_.d1->sampC));
    
    const auto writeL2 = _o_.writeL2.empty() ? tmpFile() : _o_.writeL2;

    if (o.ladC == NO_CALIBRATION)
    {
        logRun([&]()
        {
            LadderCalibrator::createBAM(_o_.inputL2, writeL2)->calibrate(o_);
        }, o, "Started ladder calibration", "Completed ladder calibration");
    }
    else if (o.ladC <= 1.0)
    {
        LadderCalibrator::Options o_;
        
        o_.d2 = std::shared_ptr<PoolData>(new PoolData());
        o_.d2->p = o.ladC; // How much to calibrate
        o_.d2->tsvL = o.work + "/" + o.tsvL1NotMerged;
        o_.r1 = o.r1;
        o_.meth = Method::Pooling;
        o_.logger = o.logger;
        o_.output = o.output;
        o_.writer = o.writer;

        logRun([&]()
        {
            assert(!o_.d2->tsvL.empty());
            LadderCalibrator::createBAM(_o_.inputL2, writeL2)->calibrate(o_);
        }, o, "Started synthetic ladder calibration", "Completed synthetic ladder calibration");
    }
    else
    {
        const auto bSeq = nLad;
        
        // Number of target sequin reads after calibration
        const auto tar = bSeq < o.ladC ? bSeq : o.ladC;

        // Scaling factor for calibration
        const auto p = bSeq ? (Proportion) tar / bSeq : 0.0;

        o.logInfo("Ladder reads before calibration: " + std::to_string(bSeq));
        o.logInfo("Target: " + std::to_string(tar));
        o.logInfo("Scaling factor: " + std::to_string(p));

        // Calibrate and return the number of reads after calibration
        SelectionCalibrator::createBAM(_o_.inputL2, writeL2)->calibrate(p, o);
    }

    // Before calibration, but after synthetic ladder normalization
    GDecoyOptions ol = o; ol.showGen = o.showGen;
    
    if (!writeL2.empty())
    {
        logRun([&]()
        {
            // Don't override what we've already done
            ol.writeS = ol.writeD = ol.writeT = ol.writeM = ""; ol.seqC = NO_CALIBRATION;

            // Analyze for ladder calibration
            r.B2 = DecoyAnalyzer::analyze(writeL2, "", ol, [&](const ParserBAM::Data &x, const DInter *r1, const DInter *r2, bool trimmed){
                g1(GDecoyStatus::AfterLC, x, r1, r2, trimmed);
            }); g2(GDecoyStatus::AfterLC);
        }, o, "Started analysis on post-ladder normalization",
              "Completed analysis on post-ladder normalization");
    }

    /*
     * Sequin calibration
     */

    auto getNormFactors = [&]()
    {
        o.info("Calculating normalisation factors");
        const auto inters = o.r1.inters();

        if (o.seqC == NO_CALIBRATION)
        {
            /*
             * Calibration by region to region
             */
            
            // Global statistics for the entire sample
            const auto samp = r.samp.r2.stats();
            
            if (o.debug)
            {
                o.logInfo("Inside getNormFactors(). sample.mean = " + std::to_string(samp.mean));
                o.logInfo("Inside getNormFactors(). sample.median = " + std::to_string(samp.p50));
            }
                        
            for (auto &i : inters)
            {
                if (isChrQ && i.first != GDecoyChrQS && i.first != GDecoyChrQR)
                {
                    continue;
                }
                
                for (auto &j : i.second.data())
                {
                    const auto &name = j.second.id();
                    const auto before = r.before.r2.find(name);
                    
                    const auto m1 = o.meth != CalibrateMethod::Custom ? o.meth : CalibrateMethod::Mean;
                    const auto m2 = o.meth != CalibrateMethod::Custom ? o.meth : CalibrateMethod::Mean;

                    // Sequin coverage before calibration (merging because it could be intersecting capture regions)
                    const auto beforeS = mergeStats(r.before.r2, name, o);

                    // Sequin coverage
                    const auto seqC = beforeS ? getLocalCoverage(*beforeS, m1) : NAN;

                    // Sequin median (only used in custom)
                    const auto samM = samp.p50;
                    
                    // Sample coverage (merging because it could be intersecting capture regions)
                    const auto sampS = mergeStats(r.samp.r2, name, o);
                                        
                    // Provisional local sample coverage
                    const auto samC1 = sampS ? getLocalCoverage(*sampS, m2) : NAN;

                    // Final sample coverage
                    auto samC2 = samC1;
                    
                    /*
                     * Custom calibration. Possible scenarios:
                     *
                     *   1. Sequin coverage > sample coverage
                     *
                     *      Anaquin was designed to calibrate region by region. However, the sequin reads would be
                     *      all lost if sample coverage is too low. What's too low? That is set by
                     *      customSequinThreshold. If it's too low, calibrate to the sample median.
                     *
                     *   2. Sample coverage > sequin coveage
                     *
                     *      Anaquin can't calibrate locally. Calibrate to the sample median.
                     */
                    
                    if (o.debug) { o.logInfo("Calibrating: " + name); }

                    if (o.meth == CalibrateMethod::Custom)
                    {
                        // Scenario 1
                        if (seqC > samC1)
                        {
                            // Sample coverage is high enough?
                            if (samC1 > o.customSequinThreshold)
                            {
                                samC2 = samC1; // Calibrate locally
                                if (o.debug) { o.logInfo("Custom 1.1"); }
                            }
                            else
                            {
                                samC2 = samM; // Don't want to lose too much sequin reads. Calibrate to sample median
                                if (o.debug) { o.logInfo("Custom 1.2"); }
                            }
                        }
                        
                        // Scenario 2
                        else if (samC1 > seqC)
                        {
                            samC2 = seqC; // No calibration
                            if (o.debug) { o.logInfo("Custom 2.1"); }
                        }
                    }
                    
                    if (o.debug)
                    {
                        o.logInfo("Sample: " + std::to_string(samC2) + " for " + name);
                        o.logInfo("Sequin: " + std::to_string(seqC)  + " for " + name);
                    }

                    // Calibrate by sample coverage
                    r.c1.norms[name] = seqC ? std::min(samC2 / seqC, 1.0) : NAN;
                }
            }
            
            assert(!r.c1.norms.empty());
        }
        else if (o.seqC > 1.0)
        {
            // What's the number of sequin reads to chrQS before calibration?
            auto sum = 0;
            
            for (auto &i : inters.at(GDecoyChrQS).data())
            {
                sum += r.before.r2.find(i.second.id())->stats().n;
            }

            const auto tar = o.seqC > sum ? sum : o.seqC;
            r.c2.p = sum ? tar / sum : 0.0;
            
            o.info("Target: "  + S0(tar));
            o.info("Scaling: " + S4(r.c2.p));
        }
        else
        {
            auto sam = 0;
            auto seq = 0;
            
            for (auto &i : inters.at(GDecoyChrQS).data())
            {
                sam += r.samp.r2.find(i.second.id())->stats().n;
                seq += r.before.r2.find(i.second.id())->stats().n;
            }

            // Number of target sequin reads after calibration
            auto tar = (o.seqC / (1.0 - o.seqC)) * sam;
            
            // Make sure our target doesn't goto zero if the sample is non-zero
            if (seq && !tar)
            {
                tar = sam;
            }
            
            // Scaling factor for calibration
            r.c2.p = (sam == 0) ? 1.0 : (seq == 0) ? 0.0 : (tar >= seq ? 1.0 : ((float) tar) / seq);

            o.info("Target: "  + S0(tar));
            o.info("Scaling: " + S4(r.c2.p));
        }
    };
    
    r.bAF = GDecoyResults::buildAF(r.before, isChrQ); // Allele frequency before calibration

    // Only perform mirror calibration if no percentage/absolue calibration
    if (o.seqC == NO_CALIBRATION && !o.inputMC.empty() && !o.writeMC.empty())
    {
        // Build normalization for each sequin
        getNormFactors();
        
        logRun([&]()
        {
            // Calibrate BAM (only on the decoy BAM)
            if (!o.inputMC.empty() && !o.writeMC.empty())
            {
                o.info("Calibrating to " + o.writeMC);
                
                if (o.seqC == NO_CALIBRATION)
                {
                    // Calibrate by sampling coverage
                    calibrateBAM(o.inputMC, o.writeMC, r.c1.norms, o);
                }
                else
                {
                    calibrateBAM(o.inputMC, o.writeMC, r.c2, o);
                }
                
                auto o_ = o;
                o_.tsvA = o_.tsvE = o_.tsvR = "";
                o_.writeM = o_.writeS = o_.writeD = o_.writeT = "";
                
                const auto f1 = twoBAM ? "" : o.writeMC;
                const auto f2 = twoBAM ? o.writeMC : "";
                
                // Analyze the calibrated BAM
                r.B3 = DecoyAnalyzer::analyze(f1, f2, o_, [&](const ParserBAM::Data &x, const DInter *r1, const DInter *r2, bool trimmed) {
                    g1(GDecoyStatus::AfterSC, x, r1, r2, trimmed);
                }); g2(GDecoyStatus::AfterSC); r.after = r.B3.decoy;
            }
        }, o, "Started mirror sequin calibration", "Completed mirror sequin calibration");

        r.aAF = GDecoyResults::buildAF(r.after, isChrQ);  // Allele frequency after calibration
    }
    else if (o.seqC != NO_CALIBRATION)
    {
        // Make sure this is the calibrated BAM
         assert(!_o_.writeC.empty());

        auto tmp = cloneO(o);
        
        tmp.ladC   = NO_CALIBRATION;
        tmp.seqC   = NO_CALIBRATION;
        tmp.writeC = ""; // Only statistics
        tmp.writeS = ""; // Only statistics
        tmp.writeD = ""; // Only statistics
        tmp.writeT = ""; // Only statistics

        r.B3 = DecoyAnalyzer::analyze(_o_.writeC, "", tmp, [&](const ParserBAM::Data &, const DInter *, const DInter *, bool) {});
    }
    
    GDecoyResults::writeF(r, o);
    GDecoyResults::writeR(r, o);
    GDecoyResults::writeV(r, o);

    DecoyAnalyzer::ErrorOptions eo(o.tsvE, o.r1, o.r2, GDecoyChrQS, o.v1, isChrQ, o);
    eo.data[&r.before] = false; eo.data[&r.after] = true;
    //assert(eo.isDecoy == isChrQ);
    
    // Generating error profiles
    DecoyAnalyzer::writeE(eo);

    // After ladder normalisation
    GDecoyResults::writeL(ol.tsvL2, r.B2, ol);

    if (r.B1.wM)
    {
        logRun([&]()
        {
            // Calibrated sequins or uncalibrated (but trimmed)
            const auto src = !o.writeMC.empty() ? o.writeMC : writeT;
            
            // Write out calibrated reads to merged
            ParserBAM::parse(src, [&](ParserBAM::Data &x, const ParserBAM::Info &)
            {
                r.B1.wM->write(x);
            });

            r.B1.wM->close();
        }, o, "Started merging sample and calibrated sequin reads",
              "Completed merging sample and calibrated sequin reads");
    }
    
    return r;
}

GResource::GResource(const Path &path, const FileName &file, const FileName &ext, Build x)
{
    const auto s1 = x == Build::chrQ ? "chrQ" : toString(x);
    const auto s2 = toString(x);

    Resource::path = (x == Build::None) ? Bundle::latest(path + "/" + file + "_", ext) :
                                                  Bundle::latest(path + "/" + s1 + "/" + file + "_" + s2, ext);
}

LResource::LResource(const Path &path, const FileName &file, const FileName &ext, Build x)
{
    Resource::path = Bundle::latest(path + "/" + toString(x) + "/" + file + "_" + toString(x), ext);
}

Bin Anaquin::GBin(const SequinID &x, std::shared_ptr<VCFLadder> v1)
{
    v1 = v1 ? v1 : Standard::instance().gen.v1();
    
    auto f = [&](const std::string &s, Bin &r)
    {
        /*
         * Sequin names can be classified directly without VCF
         */
        
             if (s == "IF") { r = Bin::IF; return true; }
        else if (s == "HP") { r = Bin::HP; return true; }
        else if (s == "MS") { r = Bin::MS; return true; }
        else if (s == "MT") { r = Bin::MT; return true; }
        else if (s == "HL") { r = Bin::HL; return true; }
        else if (s == "LD") { r = Bin::LD; return true; }
        else if (s == "SV") { r = Bin::SV; return true; }
        else if (s == "PV") { r = Bin::SV; return true; }
        else if (s == "TL") { r = Bin::SV; return true; }
        else if (s == "IM") { r = Bin::IM; return true; }
        else if (s == "VC") { r = Bin::VC; return true; }
        
        const auto v = v1 ? v1->data.find(x) : nullptr;
        
        if (v)
        {
            switch (v->gt)
            {
                case Genotype::MSI:         { r = Bin::MI; return true; }
                case Genotype::Somatic:     { r = Bin::SO; return true; }
                case Genotype::Homozygous:
                case Genotype::Heterzygous: { r = Bin::GR; return true; }
            }
        }
        
             if (s == "CM") { r = Bin::SO; return true; }
        else if (s == "CV") { r = Bin::SO; return true; }
        else if (s == "CL") { r = Bin::GR; return true; }
        else if (s == "DV") { r = Bin::GR; return true; }
        else if (s == "GV") { r = Bin::GR; return true; }
        
        return false;
    };
    
    Bin r;
    
    if (f(first(noPID(x), "_"), r))
    {
        return r;
    }

    throw std::runtime_error("Unknown binning for: " + x);
}

bool Anaquin::isHP    (const Bin x) { return x == Bin::HP; }
bool Anaquin::isMS    (const Bin x) { return x == Bin::MS; }
bool Anaquin::isSoma  (const Bin x) { return x == Bin::SO; }
bool Anaquin::isVector(const Bin x) { return x == Bin::VC; }

bool Anaquin::isHP    (const SequinID &x) { return isHP(GBin(x));     }
bool Anaquin::isMS    (const SequinID &x) { return isMS(GBin(x));     }
bool Anaquin::isSoma  (const SequinID &x) { return isSoma(GBin(x));   }
bool Anaquin::isVector(const SequinID &x) { return isVector(GBin(x)); }
