#include "data/data.hpp"
#include "stats/stats.hpp"
#include "data/resources.hpp"
#include "sequins/sequins.hpp"
#include "sequins/meta/meta.hpp"
#include "parsers/parser_csv.hpp"
#include "writers/file_writer.hpp"
#include "sequins/meta/m_split.hpp"
#include "sequins/genomics/genomics.hpp"
#include "sequins/genomics/g_broad_bam.hpp"

using namespace Anaquin;

typedef MSplit::Stats Stats;
typedef MSplit::Options Options;

static MDecoyOptions create(const FileName &writeS,
                            const FileName &writeD,
                            const FileName &writeM,
                            const FileName &writeT,
                            const FileName &tsvE,
                            const FileName &tsvR,
                            const FileName &tsvF,
                            const FileName &tsvA,
                            const FileName &tsvL1,
                            const FileName &tsvL2,
                            const FileName &index,
                            Proportion seqC,
                            Proportion ladC,
                            const Options &o)
{
    MDecoyOptions o2;
    
    if (o.debug) { o2.work = o.work;   } // Always write to working directory if debug
    else         { o2.work = tmpPath(); }

    o2.seqC    = seqC;
    o2.ladC    = ladC;
    o2.writer  = o.writer;
    o2.logger  = o.logger;
    o2.output  = o.output;
    o2.showGen = false;
    
    o2.writeS = writeS;
    o2.writeD = writeD;
    o2.writeM = writeM;
    o2.writeT = writeT;
    o2.inputC = o.seqC != NO_CALIBRATION ? o.work + "/meta_sequin.bam" : "";
    o2.writeC = o.seqC != NO_CALIBRATION ? o.work + "/meta_calibrated.bam" : "";
    o2.writeL = o.work + "/meta_ladder_calibrated.bam";

    o2.tsvE = tsvE;
    o2.tsvR = tsvR;
    o2.tsvF = tsvF;

    o2.tsvL1 = tsvL1;
    o2.tsvL2 = tsvL2;
    
    const auto &r = Standard::instance().meta;
    
    o2.l1   = r.l1();
    o2.l3   = r.l3();
    o2.edge = o.edge;
    o2.r1   = *(r.r1()); o2.r2 = *(r.r2());
    
    o2.index  = index;
    o2.writer = std::shared_ptr<Writer<>>(new FileWriter(o2.work));
    o2.errors.insert(MDecoyChrQB);
    
    o2.originalW = o.work;
    return o2;
}

static std::string errorSummary(const Stats &stats, const Options &o)
{
    const auto f = "\n\n7. SEQUENCING ERRORS            COUNT; PER_BASE; REGION_SIZE\n"
                   "Total:                          %1%\n"
                   "Mismatch:                       %2%\n"
                   "Transitions/Transversions:      %3%\n"
                   "Insertion:                      %4%\n"
                   "Deletion:                       %5%";

    if (!o.bam)
    {
        return (boost::format(f) % MISSING
                                 % MISSING
                                 % MISSING
                                 % MISSING
                                 % MISSING).str();
    }

    const auto E = MDecoyResults::reportE(stats.tsvE, stats.tsvR);
    return (boost::format(f) % E.text.at(0)
                             % E.text.at(1)
                             % E.text.at(2)
                             % E.text.at(3)
                             % E.text.at(4)).str();
}

Stats MSplit::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    o.info("Index: "   + o.index);
    o.info("Threads: " + S0(o.thr));
    o.info("Mixture: " + mixToStr(o.mix));

    Stats stats;

    if (o.bam)
    {
        const auto o_ = create(o.work + "/meta_sample.bam",
                               "",
                              !o.skipMerge ? o.work + "/meta_merged.bam" : "",
                               "",
                               "meta_errors.tsv",
                               "meta_sequin.tsv",
                               "meta_features.tsv",
                               "meta_attritures.tsv",
                               "meta_ladder.tsv",
                               "meta_ladder_calibrated.tsv",
                               o.index,
                               o.seqC,
                               o.ladC,
                               o);

        const auto r = MDecoyAnalysis(f1, o_);
        
        stats.B1 = r.B1;
        stats.B2 = r.B2;
        stats.B3 = r.B3;
        stats.tsvE = o_.work + "/" + o_.tsvE;
        stats.tsvF = o_.work + "/" + o_.tsvF;
        stats.tsvR = o_.work + "/" + o_.tsvR;

        stats.tsvL1 = o.work + "/meta_ladder.tsv";
        stats.tsvL2 = o.work + "/meta_ladder_calibrated.tsv";

        copy(o_.work + "/" + o_.tsvR,  stats.tsvR);
        copy(o_.work + "/" + o_.tsvL1, stats.tsvL1);
        copy(o_.work + "/" + o_.tsvL2, stats.tsvL2);
    }
    else
    {
        // Kallisto before calibration
        SKallisto(stats.S1, f1, f2, o);

        const auto o_ = SCalibrateDefault("sequin");

        // Calibration on sequins
        stats.S1.C = SCalibrateP(GR, o.seqC, stats.S1, o, o_);

        // Generating "_ladder.tsv"
        SWriteLadderPostCalib(Standard::instance().meta.l3(), stats.S1, o, false);

        // Generating meta_reads.tsv
        SWriteReads(Product::Meta, "meta_reads.tsv", stats.S1, o);

        if (o.seqC != NO_CALIBRATION)
        {
            // Make sure the working directory is unaffected
            auto tmp = cloneO(o); tmp.seqC = NO_CALIBRATION; tmp.ladC = NO_CALIBRATION;
            
            // Kallisto after sequin calibration
            SKallisto(stats.S2, o_.o1(o), o_.o2(o), tmp); removeD(tmp.work);
        }
        
        if (o.ladC != NO_CALIBRATION)
        {
            const auto f1 = o.work + "/meta_ladder_1.fq.gz";
            const auto f2 = o.work + "/meta_ladder_2.fq.gz";
            const auto o1 = o.work + "/meta_ladder_calibrated_1.fq.gz";
            const auto o2 = o.work + "/meta_ladder_calibrated_2.fq.gz";

            if (o.ladC <= 1.0)
            {
                LadderCalibrator::Options o_;
                o_.meth = LadderCalibrator::Method::Pooling;
                o_.d2 = std::shared_ptr<LadderCalibrator::Options::PoolData>(new LadderCalibrator::Options::PoolData());
                o_.d2->tsvL = stats.S1.tsv;
                o_.d2->p = o.ladC;

                if (stats.S1.K.r1.count(Bin::LD))
                {
                    o_.d2->reads = stats.S1.K.r1.at(Bin::LD);
                }

                // Calibration by pooling (because sample coverage not available)
                LadderCalibrator::createFQ(f1, f2, o1, o2)->calibrate(o_);
            }
            else
            {
                const auto bSeq = stats.S1.binN(Bin::LD);
                
                // Number of target sequin reads after calibration
                const auto tar = bSeq < o.ladC ? bSeq : o.ladC;

                // Scaling factor for calibration
                const auto p = bSeq ? (Proportion) tar / bSeq : 0.0;

                o.logInfo("Sequin reads: " + std::to_string(bSeq));
                o.logInfo("Target: " + std::to_string(tar));
                o.logInfo("Scaling factor: " + std::to_string(p));

                // Calibrate and return the number of reads after calibration
                SelectionCalibrator::createFQ(f1, f2, o1, o2)->calibrate(p, o);
            }
            
            // Make sure the working directory is unaffected
            auto tmp = cloneO(o); tmp.seqC = NO_CALIBRATION; tmp.ladC = NO_CALIBRATION;

            // Kallisto after ladder calibration
            SKallisto(stats.S3, o1, o2, tmp); removeD(tmp.work);
            
            // Generating "_ladder_calibrated.tsv"
            SWriteLadderPostCalib(Standard::instance().meta.l3(), stats.S3, o, true);
        }
        
        if (!o.skipMerge)
        {
            o.info("Merging sample, sequin and ladder reads");
            
            const auto seq_1 = (!o.isSCalib()) ? o.work + "/meta_sequin_1.fq.gz" : o.work + "/meta_sequin_calibrated_1.fq.gz";
            const auto seq_2 = (!o.isSCalib()) ? o.work + "/meta_sample_2.fq.gz" : o.work + "/meta_sequin_calibrated_2.fq.gz";
            const auto lad_1 = (!o.isLCalib()) ? o.work + "/meta_ladder_1.fq.gz" : o.work + "/meta_ladder_calibrated_1.fq.gz";
            const auto lad_2 = (!o.isLCalib()) ? o.work + "/meta_ladder_2.fq.gz" : o.work + "/meta_ladder_calibrated_2.fq.gz";
            mergeFQ(std::vector<FileName> { o.work + "/meta_sample_1.fq.gz", seq_1, lad_1 },
                    std::vector<FileName> { o.work + "/meta_sample_2.fq.gz", seq_2, lad_2 },
                    o.work + "/meta_merged_1.fq.gz", o.work + "/meta_merged_2.fq.gz");
        }
    }

    return stats;
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const FileName &tsv1, const FileName &tsv2, const Stats &stats, const Options &o)
{
    o.generate(file);
    o.writer->open(file);
    
    const auto f = "SEQUIN META REPORT\n"
                   "1. ANALYSIS\n"
                   "Date:                           %1%\n"
                   "Anaquin version:                %2%\n"
                   "Command:                        %3%\n"
                   "Sequin mixture:                 %4%\n"
                   "Mixture type:                   %5%\n"
                   "Reference index:                %6%\n"
                   "Reference regions:              %7%\n"
                   "Sample reads path:              %8%\n"
                   "Sequin reads path:              %9%\n"
                   "Calibrated reads path:          %10%\n"
                   "Anaquin k-mer length:           %11%\n"
                   "Threshold:                      %12%\n\n"
                   "2. LIBRARY DETAILS\n"
                   "Instrument ID:                  %13%\n"
                   "Run number:                     %14%\n"
                   "Flowcell ID:                    %15%\n"
                   "Lane:                           %16%\n\n"
                   "3. LIBRARY FRACTIONS\n"
                   "Sample reads; fraction:         %17%\n"
                   "Sequin reads; fraction:         %18%\n"
                   "Ladder reads; fraction:         %19%\n"
                   "Vector reads; fraction:         %20%\n"
                   "Total reads:                    %21%\n"
                   "Sequin dilution:                %22%\n\n"
                   "4. CALIBRATION SUMMARY\n"
                   "Sequin calibration:             %23%\n"
                   "Sequin calibration factor:      %24%\n"
                   "Sequin reads after calibration: %25%\n"
                   "Ladder calibration:             %26%\n"
                   "Ladder reads after calibration: %27%\n"
                   "Total reads after calibration:  %28%\n\n"
                   "5. LADDER - LIBRARY QUALITY\n"
                   "Slope:                          %29%\n"
                   "R2:                             %30%\n"
                   "Mean Ratio:                     %31%\n"
                   "Ladder table:                   %32%\n\n"
                   "6. SEQUIN QUANTIFICATION\n"
                   "Slope:                          %33%\n"
                   "R2:                             %34%\n"
                   "Sequin table:                   %35%";

    const auto tmp = tmpFile();
    
    // Only ladder sequins
    RGrep(tsv1, tmp, "NAME", "LD_"); const auto l1 = RLinear(tmp, "NAME", "COPY", "READ").linear();
    
    // Only metagenomics sequins
    RGrep(tsv2, tmp, "NAME", "MG_"); const auto l2 = RLinear(tmp, "NAME", "MIX", "TPM").linear(true, true);

    const auto p1 = o.bam ? (CommonResults *) &stats.B1 : (CommonResults *) &stats.S1;
    const auto p2 = o.bam ? (CommonResults *) &stats.B2 : (CommonResults *) &stats.S2;
    const auto p3 = o.bam ? (CommonResults *) &stats.B3 : (CommonResults *) &stats.S3;

    // Parition before calibration
    const auto pc = p1->pc();
    
    // Number of calibrated sequin reads
    const auto cs = o.isSCalib() ? p2->total() : 0;
    
    // Number of calibrated ladder reads
    const auto ls = o.isLCalib() ? p3->total() : 0;
    
    o.writer->write((boost::format(f) % date()
                                      % o.version
                                      % o.cmd
                                      % "META v3 Mix"
                                      % mixToStr(o.mix)
                                      % o.index
                                      % Standard::instance().meta.r1()->src
                                      % (o.work + "/meta_sample*")         // 8
                                      % (o.work + "/meta_sequin*")         // 9
                                      % (o.seqC == NO_CALIBRATION ? MISSING : o.work + "/meta_sequin_calibrated*")      // 10
                                      % (o.bam ? MISSING : S0(o.k))        // 11
                                      % (o.bam ? MISSING : S2(o.rule))     // 12
                                      % p1->lib().inst(p1->lib().format()) // 13
                                      % p1->lib().run(p1->lib().format())  // 14
                                      % p1->lib().flow(p1->lib().format()) // 15
                                      % p1->lib().lane(p1->lib().format()) // 16
                                      % (S0(pc.binN(ES)) + " ; " + S2(pc.binP(ES)))
                                      % (S0(pc.binN(GR)) + " ; " + S2(pc.binP(GR)))
                                      % (S0(pc.binN(LD)) + " ; " + S2(pc.binP(LD)))
                                      % (S0(pc.binN(VC)) + " ; " + S2(pc.binP(VC)))
                                      % p1->total()                         // 21
                                      % S2(p1->dil())                       // 22
                                      % calib2str(o.seqC)                   // 23
                                      % (!o.isSCalib() ? MISSING : S2(p1->calibF(Calibration::Sequin)))
                                      % (!o.isSCalib() ? MISSING : S2(cs))  // 25
                                      % calib2str(o.ladC)                   // 26
                                      % (o.isLCalib() ? S0(ls) : MISSING)   // 27
                                      % ((!o.isSCalib() && !o.isLCalib()) ? MISSING : S0(pc.binN(ES) + cs + ls))
                                      % replaceNA(l1.m)                     // 29
                                      % replaceNA(l1.R2)                    // 30
                                      % replaceNA(RLadTable(tsv1, tmp, "NAME", o.bam ? "MEAN" : "Q50")) // 31
                                      % (o.work + "/meta_ladder_table.tsv") // 32
                                      % replaceNA(l2.m)                     // 33
                                      % replaceNA(l2.R2)                    // 34
                                      % (o.work + "/meta_sequin_table.tsv") // 35
                     ).str() + errorSummary(stats, o));
    o.writer->close();
}

static void writeQuin(const FileName &file, const CustomMap<SequinID, Read> &r1, const Options &o)
{
    const auto f = "%1%\t%2%\t%3%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % "NAME"
                                      % "MIX"
                                      % "READ").str());
    
    auto l1 = Standard::instance().meta.l1();
    
    for (const auto &seq : l1->seqs)
    {
        const auto mix = (MBin(seq) == LD || MBin(seq) == IF || MBin(seq) == VC) ? MISSING :
                        (l1->contains(seq, o.mix) ? S6(l1->input(seq, o.mix)) : MISSING);
        
        o.writer->write((boost::format(f) % seq
                                          % mix
                                          % (r1.count(seq) ? r1.at(seq) : 0)).str());
    }
    
    o.writer->close();
}

static void addTPM(const FileName &file, const Options &o)
{
    const auto scale = RSum(o.work + "/" + file, "READ") / 1000000;
    auto j = -1;

    const auto tmp = file + ".tmp";
    o.writer->open(tmp);
    
    // Assume only the READ column
    ParserCSV::parse(Reader(o.work + "/" + file), [&](const ParserCSV::Data &x, Progress i)
    {
        if (!i)
        {
            auto p = std::find(x.begin(), x.end(), "READ");
            assert(p != x.end());
            j = p - x.begin();
            
            for (auto k = 0; k < x.size(); k++)
            {
                o.writer->write(x[k] + "\t", false);
            }

            o.writer->write("TPM", false);
        }
        else
        {
            assert(j >= 0);
            o.writer->write("\n", false);
            
            for (auto k = 0; k < x.size(); k++)
            {
                o.writer->write(x[k] + "\t", false);

                if (k == x.size() - 1)
                {
                    o.writer->write(S4(stod(x[j]) / scale), false);
                }
            }
        }
    }, "\t");

    o.writer->write("");
    o.writer->close();
    mv(o.work + "/" + tmp, o.work + "/" + file);
}

void MSplit::report(const FileName &f1, const FileName &f2, const Options &o)
{
    const auto stats = analyze(f1, f2, o);
    
    // Statistics before calibration
    const auto p1 = o.bam ? (CommonResults *) &stats.B1 : (CommonResults *) &stats.S1;

    writeQuin("meta_sequin.tsv", p1->rn(), o); addTPM("meta_sequin.tsv", o);

    // Generating meta_sequin_calibrated.tsv
    if (o.isSCalib())
    {
        // Statistics after calibration
        const auto p2 = o.bam ? (CommonResults *) &stats.B2 : (CommonResults *) &stats.S2;
        
        writeQuin("meta_sequin_calibrated.tsv", p2->rn(), o); addTPM("meta_sequin_calibrated.tsv", o);
    }

    const auto tsv1 = o.isLCalib() ? o.work + "/meta_ladder_calibrated.tsv" : o.work + "/meta_ladder.tsv";
    const auto tsv2 = o.isSCalib() ? o.work + "/meta_sequin_calibrated.tsv" : o.work + "/meta_sequin.tsv";
    
    // Generating meta_ladder.tsv
    writeLTable(o.work + "/meta_ladder.tsv", "meta_ladder_table.tsv", o);

    // Generating sequin abundance table
    writeSTable_1(tsv2, "meta_sequin_table.tsv", o, 6, 6, 6, "MIX", "READ");

    // Generating meta_report.txt
    writeSummary("meta_report.txt", f1, f2, tsv1, tsv2, stats, o);
}
