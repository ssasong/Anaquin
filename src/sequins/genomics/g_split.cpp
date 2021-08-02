#include "sequins/sequins.hpp"
#include "tools/calibrator.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_bam.hpp"
#include "writers/file_writer.hpp"
#include "sequins/genomics/g_split.hpp"
#include "sequins/genomics/genomics.hpp"
#include "sequins/genomics/g_broad_bam.hpp"

using namespace Anaquin;

typedef GSplit::Stats Stats;
typedef GSplit::Options Options;

Stats GSplit::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    o.info("Index: "   + o.index);
    o.info("Threads: " + S0(o.thr));

    Stats stats;
    
    if (o.bam)
    {
        auto tmp = GDecoyOptions::create(o.index, o);

        tmp.seqC  = o.seqC;
        tmp.ladC  = o.ladC;

        tmp.tsvE  = "split_errors.tsv";
        tmp.tsvR  = "split_region.tsv";
        tmp.tsvF  = "split_features.tsv";
        tmp.tsvA  = "split_variant.tsv";
        tmp.tsvL1 = "split_ladder.tsv";
        tmp.tsvL2 = "split_ladder_calibrated.tsv";

        tmp.inputC  = o.seqC != NO_CALIBRATION ? o.work + "/split_sequin.bam" : "";
        tmp.writeC  = o.work + "/split_calibrated.bam";
        tmp.writeS  = o.work + "/split_sample.bam";
        tmp.writeD  = "";
        tmp.writeT  = "";
        tmp.writeM  = !o.skipMerge ? o.work + "/split_merged.bam" : "";
        tmp.writeMC = "";
        tmp.inputL2 = o.work + "/split_ladder.bam";
        tmp.writeL1 = o.work + "/split_ladder.bam";
        tmp.writeL2 = o.work + "/split_ladder_calibrated.bam"; //o.ladC != NO_CALIBRATION ? o.work + "/split_ladder_calibrated.bam" : "";

        // Never allow mirror calibation
        assert(tmp.writeMC.empty());
        
        std::shared_ptr<BAMWriter> sqW, svW, imW, miW, mtW, hlW, ifW, vcW;
        
        sqW = std::shared_ptr<BAMWriter>(new BAMWriter()); sqW->open(o.work + "/split_sequin.bam");
        svW = std::shared_ptr<BAMWriter>(new BAMWriter()); svW->open(o.work + "/split_structural.bam");
        imW = std::shared_ptr<BAMWriter>(new BAMWriter()); imW->open(o.work + "/split_immune.bam");
        hlW = std::shared_ptr<BAMWriter>(new BAMWriter()); hlW->open(o.work + "/split_hla.bam");
        ifW = std::shared_ptr<BAMWriter>(new BAMWriter()); ifW->open(o.work + "/split_info.bam");
        mtW = std::shared_ptr<BAMWriter>(new BAMWriter()); mtW->open(o.work + "/split_mito.bam");
        vcW = std::shared_ptr<BAMWriter>(new BAMWriter()); vcW->open(o.work + "/split_vector.bam");

        auto pc_ = ParitionCount();
        auto r = GDecoyAnalysis(f1, f2, tmp, [&](GDecoyStatus s,
                                                 const ParserBAM::Data &x,
                                                 const DInter *r1,
                                                 const DInter *r2,
                                                 bool trimmed)
        {
            if (s != GDecoyStatus::Before)
            {
                return;
            }
           
            // Skip the read if trimmed (Tim's request to make trimmed reads going to VC)
            if (trimmed) { vcW->write(x); pc_[Bin::VC]++; return; }

            if (!isDecoy(x.cID))      { pc_[Bin::ES]++; return; }
            if (x.cID == GDecoyChrQL) { pc_[Bin::LD]++; return; }

            if (r2)
            {
                Bin b;
                
                switch (b = GBin(r2->name()))
                {
                    case HL: { hlW->write(x); break; }
                    case LD: { break; }
                    case SV: { svW->write(x); break; }
                    case IM: { imW->write(x); break; }
                    case MT: { mtW->write(x); break; }
                    case IF: { ifW->write(x); break; }
                    case VC: { vcW->write(x); break; }
                    default: { sqW->write(x); break; }
                }
                
                pc_[b]++;
            }
            else
            {
                sqW->write(x);
                pc_[Bin::GR]++;
            }            
        }, [&](GDecoyStatus s) {
            if (s == GDecoyStatus::Before)
            {
                sqW->close(); svW->close(); imW->close(); hlW->close(); ifW->close(); mtW->close(); vcW->close();
                sqW = svW = imW = hlW = ifW = mtW = vcW = nullptr;
            }
        }); r.B1._pc = pc_;

        stats.B1 = r.B1;
        stats.B2 = r.B3;
        stats.B3 = r.B2;
        
        stats.tsvV = o.work + "/split_variant.tsv";
        stats.tsvE = tmp.work + "/" + tmp.tsvE;
        stats.tsvR = tmp.work + "/" + tmp.tsvR;
        stats.tsvF = tmp.work + "/" + tmp.tsvF;
        
        copy(tmp.work + "/" + tmp.tsvA, o.work + "/split_variant.tsv");
        copy(tmp.work + "/split_ladder.tsv", o.work + "/split_ladder.tsv");
        
        // TODO: This should have been done in the workflow
        if (o.ladC == NO_CALIBRATION) { rm(tmp.writeL2); }
    }
    else
    {
        // Kallisto before calibration
        SKallisto(stats.S1, f1, f2, o);
        
        const auto o_ = SCalibrateDefault("sequin");

        // Calibration on sequins formed by several categories
        stats.S1.C = SCalibrateP(std::vector<Bin> { GR, SO, MI, MS, HP }, o.seqC, stats.S1, o, o_);

        // Generating "_ladder.tsv"
        SWriteLadderPostCalib(Standard::instance().gen.l3(), stats.S1, o, false);

        // Generating split_reads.tsv
        SWriteReads(Product::Genomics, "split_reads.tsv", stats.S1, o);

        if (o.seqC != NO_CALIBRATION)
        {
            // Make sure the working directory is unaffected
            auto tmp = cloneO(o); tmp.seqC = NO_CALIBRATION; tmp.ladC = NO_CALIBRATION;
            
            // Kallisto after sequin calibration
            SKallisto(stats.S2, o_.o1(o), o_.o2(o), tmp); removeD(tmp.work);
        }

        if (o.ladC != NO_CALIBRATION)
        {
            const auto f1 = o.work + "/split_ladder_1.fq.gz";
            const auto f2 = o.work + "/split_ladder_2.fq.gz";
            const auto o1 = o.work + "/split_ladder_calibrated_1.fq.gz";
            const auto o2 = o.work + "/split_ladder_calibrated_2.fq.gz";

            if (o.ladC <= 1.0)
            {
                LadderCalibrator::Options o_;
                o_.meth = LadderCalibrator::Method::Pooling;
                o_.d2 = std::shared_ptr<LadderCalibrator::Options::PoolData>(new LadderCalibrator::Options::PoolData());
                o_.d2->tsvL = stats.S1.tsv;
                o_.d2->p = o.ladC;
                o_.d2->reads = stats.S1.K.r1.at(Bin::LD);

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
            SWriteLadderPostCalib(Standard::instance().gen.l3(), stats.S3, o, true);
        }
        
        if (!o.skipMerge)
        {
            o.info("Merging sample and sequin");
            
            const auto seq_1 = !o.isSCalib() ? o.work + "/split_sequin_1.fq.gz" : o.work + "/split_sequin_calibrated_1.fq.gz";
            const auto seq_2 = !o.isSCalib() ? o.work + "/split_sequin_2.fq.gz" : o.work + "/split_sequin_calibrated_2.fq.gz";
            
            mergeFQ(std::vector<FileName> { o.work + "/split_sample_1.fq.gz", seq_1  },
                    std::vector<FileName> { o.work + "/split_sample_2.fq.gz", seq_2  },
                    o.work + "/split_merged_1.fq.gz", o.work + "/split_merged_2.fq.gz");
        }
    }
    
    return stats;
}

static std::string errorSummary(const Stats &stats, const Options &o)
{
    if (!o.bam)
    {
        return "";
    }
    
    const auto E = GBroadBam::reportE(stats.tsvE, stats.tsvR, stats.tsvF, false);

    const auto f =  "\n\nSEQUENCING ERRORS               COUNT; PER_BASE; REGION_SIZE\n"
                    "Total:                          %1%\n"
                    "Mismatch:                       %2%\n"
                    "Transitions/Transversions:      %3%\n"
                    "Insertion:                      %4%\n"
                    "Deletion:                       %5%\n"
                    "Coding:                         %6%\n"
                    "Noncoding:                      %7%";

    return (boost::format(f) % E.text.at(0)
                             % E.text.at(1)
                             % E.text.at(2)
                             % E.text.at(3)
                               % E.text.at(4)
                             % E.text.at(5)
                             % E.text.at(6)).str();
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const Stats &stats, const Options &o)
{
    o.generate(file);
    o.writer->open(file);

    const auto f = "SEQUIN SPLIT REPORT\n\n"
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
                   "Ladder reads path:              %10%\n"
                   "Structural variant reads path:  %11%\n"
                   "Immune reads path:              %12%\n"
                   "Mitochondria reads path:        %13%\n"
                   "HLA reads path:                 %14%\n"
                   "Information reads path:         %15%\n"
                   "Vector reads path:              %16%\n"
                   "Calibrated reads path:          %17%\n"
                   "Anaquin k-mer length:           %18%\n"
                   "Threshold:                      %19%\n\n"
                   "2. LIBRARY DETAILS\n"
                   "Instrument ID:                  %20%\n"
                   "Run number:                     %21%\n"
                   "Flowcell ID:                    %22%\n"
                   "Lane:                           %23%\n\n"
                   "3. LIBRARY FRACTIONS\n"
                   "Sample reads; fraction:         %24%\n"
                   "Sequin reads; fraction:         %25%\n"
                   "Ladder reads; fraction:         %26%\n"
                   "Structural reads; fraction:     %27%\n"
                   "Immune reads; fraction:         %28%\n"
                   "Mitochondria; fraction:         %29%\n"
                   "HLA reads; fraction:            %30%\n"
                   "Information reads; fraction:    %31%\n"
                   "Vector reads; fraction:         %32%\n"
                   "Total reads:                    %33%\n"
                   "Sequin dilution:                %34%\n\n"
                   "4. CALIBRATION SUMMARY\n"
                   "Sequin calibration:             %35%\n"
                   "Sequin calibration factor:      %36%\n"
                   "Sequin reads after calibration: %37%\n"
                   "Ladder calibration:             %38%\n"
                   "Ladder reads after calibration: %39%\n"
                   "Total reads after calibration:  %40%\n\n"
                   "5. LADDER - LIBRARY QUALITY\n"
                   "Slope:                          %41%\n"
                   "R2:                             %42%\n"
                   "Mean Ratio:                     %43%\n"
                   "Ladder table:                   %44%\n\n"
                   "6. SEQUIN QUANTIFICATION\n"
                   "Slope:                          %45%\n"
                   "R2:                             %46%\n"
                   "Sequin table:                   %47%";

    const auto tmp = tmpFile();

    RGrep(o.work + "/split_ladder.tsv", tmp, "NAME", "LD_"); const auto l1 = RLinear(tmp, "NAME", "COPY", "READ").linear(true, true);
    RGrep(o.work + "/split_variant.tsv", tmp, "LABEL", "Somatic"); const auto l2 = RLinear(tmp, "NAME", "EXP_FREQ", "OBS_FREQ").linear(true, true);

    const auto &r = Standard::instance().gen;
    
    const auto p1 = o.bam ? (CommonResults *) &stats.B1 : (CommonResults *) &stats.S1;
    const auto p2 = o.bam ? (CommonResults *) &stats.B2 : (CommonResults *) &stats.S2;
    const auto p3 = o.bam ? (CommonResults *) &stats.B3 : (CommonResults *) &stats.S3;

    // Parition before calibration
    const auto pc = p1->pc();
    
    // Number of calibrated sequin reads
    const auto cs = o.isSCalib() ? p2->total() : 0;
    
    // Number of calibrated ladder reads
    const auto ls = o.isLCalib() ? p3->total() : 0;

    const auto gn = pc.binN(GR) + pc.binN(SO) + pc.binN(MI) + pc.binN(MS) + pc.binN(HP);
    const auto gp = pc.binP(GR) + pc.binP(SO) + pc.binP(MI) + pc.binP(MS) + pc.binP(HP);
    assert(gp >= 0 && gp <= 100);
    
    o.writer->write((boost::format(f) % date()
                                      % o.version
                                      % o.cmd
                                      % "Genome Mix"
                                      % "A"
                                      % o.index
                                      % Standard::instance().gen.r1()->src
                                      % (o.work + "/split_sample*") // 8
                                      % (o.work + "/split_sequin*") // 9
                                      % (o.work + "/split_ladder*") // 10
                                      % (o.work + "/split_sv*")     // 11
                                      % (o.work + "/split_immune*") // 12
                                      % (o.work + "/split_mito*")   // 13
                                      % (o.work + "/split_hla*")    // 14
                                      % (o.work + "/split_info*")   // 15
                                      % (o.work + "/split_vector*") // 16
                                      % (o.seqC == NO_CALIBRATION ? MISSING : o.work + "/split_sequin_calibrated*")
                                      % (o.bam ? MISSING : S0(o.k))        // 18
                                      % (o.bam ? MISSING : S2(o.rule))     // 19
                                      % p1->lib().inst(p1->lib().format()) // 20
                                      % p1->lib().run(p1->lib().format())  // 21
                                      % p1->lib().flow(p1->lib().format()) // 22
                                      % p1->lib().lane(p1->lib().format()) // 23
                                      % (S0(pc.binN(ES)) + " ; " + S2(pc.binP(ES)))
                                      % (S0(gn) + " ; " + S2(gp)) // 25
                                      % (S0(pc.binN(LD)) + " ; " + S2(pc.binP(LD)))
                                      % (S0(pc.binN(SV)) + " ; " + S2(pc.binP(SV)))
                                      % (S0(pc.binN(IM)) + " ; " + S2(pc.binP(IM)))
                                      % (S0(pc.binN(MT)) + " ; " + S2(pc.binP(MT)))
                                      % (S0(pc.binN(HL)) + " ; " + S2(pc.binP(HL)))
                                      % (S0(pc.binN(IF)) + " ; " + S2(pc.binP(IF)))
                                      % (S0(pc.binN(VC)) + " ; " + S2(pc.binP(VC)))
                                      % p1->total()                         // 33
                                      % S2(p1->dil())                       // 34
                                      % calib2str(o.seqC)                   // 35
                                      % (!o.isSCalib() ? MISSING : S2(p1->calibF(Calibration::Sequin)))
                                      % (!o.isSCalib() ? MISSING : S2(cs))  // 37
                                      % calib2str(o.ladC)                   // 38
                                      % (o.isLCalib() ? S0(ls) : MISSING)   // 39
                                      % ((!o.isSCalib() && !o.isLCalib()) ? MISSING : S0(pc.binN(ES) + cs + ls))
                                      % replaceNA(l1.m)                      // 41
                                      % replaceNA(l1.R2)                     // 42
                                      % S2(RLadTable(o.work + "/split_ladder.tsv", tmp, "NAME", o.bam ? "MEAN" : "Q50"))
                                      % (o.work + "/split_ladder_table.tsv") // 44
                                      % replaceNA(l2.m)                      // 45
                                      % replaceNA(l2.R2)                     // 46
                                      % (o.work + "/split_sequin_table.tsv") // 47
            ).str() + errorSummary(stats, o));
    o.writer->close();
}

GVariantResults GSplit::Stats::V1(const AnalyzerOptions &o) const
{
    const auto &r = Standard::instance().gen;
    GVariantResults v;
    
    if (o.bam)
    {
        // Allele frequency from alignments
        const auto m = GDecoyResults::buildAF(this->B1.decoy, false);

        for (auto &i : r.v4()->data.vars())
        {
            v.R[i.name] = m.count(i.name) ? m.at(i.name).R : 0;
            v.V[i.name] = m.count(i.name) ? m.at(i.name).V : 0;
        }
    }
    else
    {
        for (const auto &i : S1.K.rSeqs)
        {
            const auto std = i.first;
            const auto R_  = S1.K.rSeqs.at(std);
            
            if (!S1.K.aSeqs.count(std) || (GBin(i.first) != Bin::GR && GBin(i.first) != Bin::SO))
            {
                continue;
            }
            
            for (const auto &j : S1.K.aSeqs.at(std))
            {
                const auto V_  = j;
                const auto exp = r.l1()->findBySub(noPID(V_));
                
                if (!exp)
                {
                    continue;
                }
                
                assert(isSubstr(V_, "_A"));
                
                const auto hasR_ = !R_.empty() ? findBySub(this->S1.K.uniqs, noPID(R_)) : nullptr;
                const auto hasV_ = findBySub(this->S1.K.uniqs, noPID(V_));
                const auto hasR  = hasR_ && hasR_->size() >= 3;
                const auto hasV  = hasV_ && hasV_->size() >= 3;
                
                const auto R = v.R[R_] = hasR ? med(toVector(*hasR_)) : 0;
                const auto V = v.V[V_] = hasV ? med(toVector(*hasV_)) : 0;
            }
        }
    }
    
    return v;
}

void GSplit::writeV(const FileName &file, const Stats &stats, const SOptions &o)
{
    if (o.bam)
    {
        return; // Copied from DecoyAnalyzer in analyze()
    }
    
    const auto &r = Standard::instance().gen;
    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%";

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % "NAME"
                                      % "LABEL"
                                      % "TYPE"
                                      % "GENOTYPE"
                                      % "CHROM"
                                      % "POSITION"
                                      % "EXP_FREQ"
                                      % "OBS_FREQ"
                                      % "REF_COUNT"
                                      % "VAR_COUNT").str());
    
    const auto V1 = stats.V1(o);
    assert(V1.R.size() == V1.V.size());

    // hg19/hg38? or chrQ?
    const auto vr = o.bam ? r.v4() : r.v1();

    for (const auto &i : V1.V)
    {
        const auto &x = o.bam ? i.first : toStand(i.first);
        
        const auto R = V1.R.find(noPID(x)); // Reference counts
        const auto V = V1.V.find(noPID(x)); // Variant counts
        assert(R && V);
        
        const auto exp = r.l1()->findBySub(noPID(x));
        const auto obs = *R + *V ? *V / (*R + *V) : 0;

        const auto v = vr->data.find(x, true);
        
        if (v)
        {
            o.writer->write((boost::format(f) % i.first
                                            % bin2Label(GBin(x))
                                            % var2str(v->type())
                                            % gt2str(v->gt)
                                            % v->cID
                                            % v->l.start
                                            % S6(exp)
                                            % S6(obs)
                                            % (R ? S0(*R) : MISSING)
                                            % (V ? S0(*V) : MISSING)).str());
        }
    }
    
    o.writer->close();
}

static void writeSTable(const FileName &src, const FileName &dst, const SOptions &o, Count nExp, Count nObs, Count nCV)
{
    const auto tmp = tmpFile();
    RMeanCV(src, tmp, "EXP_FREQ", "OBS_FREQ", nExp, nObs, nCV);
    auto t = readFile(tmp);
         t = replace(t, "NAME",  "EXPECTED_FREQ");
         t = replace(t, "COUNT", "SEQUIN_DETECTED");
         t = replace(t, "MEAN",  "MEAN_OBS_FREQ");
    o.generate(dst);
    o.writer->open(dst);
    o.writer->write(t);
    o.writer->close();
}

void GSplit::report(const FileName &f1, const FileName &f2, const Options &o)
{
    const auto stats = analyze(f1, f2, o);
    
    // Generating split_variant.tsv
    GSplit::writeV("split_variant.tsv", stats, o);
    
    // Generating split_ladder.tsv
    writeLTable(o.work + "/split_ladder.tsv", "split_ladder_table.tsv", o);
    
    const auto tmp = tmpFile();
    
    RGrep(o.work + "/split_variant.tsv", tmp, "LABEL", "Somatic");
    
    // Generating split_variant_table.tsv
    writeSTable(tmp, "split_variant_table.tsv", o, 6, 6, 6);

    // Generating split_report.txt
    writeSummary("split_report.txt", f1, f2, stats, o);
}
