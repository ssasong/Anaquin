#include "data/resources.hpp"
#include "sequins/rna/rna.hpp"
#include "parsers/parser_fa.hpp"
#include "sequins/rna/r_split.hpp"
#include "writers/file_writer.hpp"

using namespace Anaquin;

typedef RSplit::Stats Stats;
typedef RSplit::Options Options;

Stats RSplit::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    o.info("Index: "   + o.index);
    o.info("Threads: " + S0(o.thr));
    o.info("Mixture: " + mixToStr(o.mix));
    
    Stats stats;

    if (o.bam)
    {
        DecoyAnalyzer::Options o2;
        
        o2.edge    = 550;
        o2.showGen = false;
        o2.logger  = o.logger;
        o2.writer  = o.writer;
        o2.output  = o.output;
        o2.shouldTrim = false; // Never trim in RNA

        o2.writeS = o.work + "/rna_sample.bam";
        o2.writeD = o.work + "/rna_sequin.bam";
        o2.writeM = o.skipMerge ? o.work + "/rna_merged.bam" : "";
        o2.inputC = o.work + "/rna_sequin.bam"; // Calibrate o2.writeD
        o2.writeC = o.work + "/rna_calibrated.bam";
        
        const auto &r = Standard::instance().rna;

        o2.r1 = *(r.r3()); // No trimming at the gene level
        o2.r2 = *(r.r4()); // With trimming at the gene level
        
        o2.seqC   = o.seqC;
        o2.index  = o.index;
        o2.writer = std::shared_ptr<Writer<>>(new FileWriter(o.work));
        o2.errors.insert(GDecoyChrIS);
        o2.shouldError = false; // Skip all error calculations
        
        ParserFA::parse(Reader(o.index), [&](const ParserFA::Data &x)
        {
            o2.seqs[x.id] = x.seq;
        });
        
        assert(!o2.seqs.empty());
        assert(o2.writeT.empty());      // RNA not supporting trimming
        assert(o2.h1.inters().empty()); // Entire human chromosomes
        assert(o2.h2.inters().empty()); // Entire human chromosomes

        // Statistics before calibration
        stats.B1 = DecoyAnalyzer::analyze(f1, "", o2);

        if (stats.B1.wM)
        {
            assert(exists(o2.writeD));
            
            ParserBAM::parse(o2.writeD, [&](ParserBAM::Data &x, const ParserBAM::Info &)
            {
                stats.B1.wM->write(x);
            });

            stats.B1.wM->close();
        }

        if (o2.seqC != NO_CALIBRATION)
        {
            o2.seqC   = NO_CALIBRATION; // No need to do it again
            o2.writeS = ""; // Only statistics
            o2.writeD = ""; // Only statistics
            o2.writeM = ""; // Only statistics

            stats.B2 = DecoyAnalyzer::analyze(o2.writeC, "", o2);
        }
    }
    else
    {
        // Kallisto before calibration
        SKallisto(stats.S1, f1, f2, o);
        
        // Calibration on sequins
        stats.S1.C = SCalibrateP(GR, o.seqC, stats.S1, o, SCalibrateDefault("sequin"));
        
        // Kallisto after calibration
        auto o_ = cloneO(o); SKallisto(stats.S2, stats.S1.C.o1, stats.S1.C.o2, o_); removeD(o_.work);
        
        if (!o.skipMerge)
        {
            o.info("Merging sample and sequin");
            
            const auto seq_1 = !o.isSCalib() ? o.work + "/rna_sequin_1.fq.gz" : o.work + "/rna_sequin_calibrated_1.fq.gz";
            const auto seq_2 = !o.isSCalib() ? o.work + "/rna_sequin_2.fq.gz" : o.work + "/rna_sequin_calibrated_2.fq.gz";
         
            mergeFQ(std::vector<FileName> { o.work + "/rna_sample_1.fq.gz", seq_1  },
                    std::vector<FileName> { o.work + "/rna_sample_2.fq.gz", seq_2  },
                    o.work + "/rna_merged_1.fq.gz", o.work + "/rna_merged_2.fq.gz");
        }
    }
    
    return stats;
}

static void createTables(const FileName &tsv, const Options &o, LinearModel &l1, LinearModel &l2)
{
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    const auto tmp3 = tmpFile();

    // Sequins at the isoform level
    RGrep(tsv, tmp1, "NAME", "R"); l1 = RLinear(tmp1, "NAME", "MIX", "TPM").linear();
    
    if (!o.bam)
    {
        // Generating sequin abundance table
        writeSTable_1(tmp1, "rna_sequin_isoform_table.tsv", o, 6, 6, 6, "MIX", "READ", 0);
    }
    
    // Gene expression threshold
    const auto p = 0.05;
    
    // Sequins at the gene level by TPM
    RFilterC(tmp1, tmp2, std::set<Label> { "GENE", "MIX", "TPM" }, true); RAggregateSum(tmp2, tmp3, "GENE", p, Imputation::ToZero);

    // Regression at the gene level
    l2 = RLinear(tmp3, "GENE", "MIX", "TPM").linear();

    // Sequins at the gene level by READ
    RFilterC(tmp1, tmp2, std::set<Label> { "GENE", "MIX", "READ" }, true); RAggregateSum(tmp2, tmp3, "GENE", p, Imputation::ToZero);

    // Generating sequin abundance table
    writeSTable_1(tmp3, "rna_sequin_gene_table.tsv", o, 6, 6, 6, "MIX", "READ", p);
}

static void writeSummary(const FileName &file, const FileName &f1, const FileName &f2, const Stats &stats, const FileName &tsv, const Options &o)
{
    LinearModel l1, l2;
    
    // Generating "rna_sequin_isoform_table.tsv" and "rna_sequin_gene_table.tsv"
    createTables(tsv, o, l1, l2);

    const auto f = "SEQUIN RNA REPORT\n"
                   "1. ANALYSIS\n"
                   "Date:                           %1%\n"
                   "Anaquin version:                %2%\n"
                   "Command:                        %3%\n"
                   "Sequin mixture:                 RNA v2 Mix\n"
                   "Mixture type:                   %4%\n"
                   "Reference index:                %5%\n"
                   "Reference regions:              %6%\n"
                   "Sample reads path:              %7%\n"
                   "Sequin reads path:              %8%\n"
                   "Calibrated reads path:          %9%\n"
                   "Anaquin k-mer length:           %10%\n"
                   "Threshold:                      %11%\n\n"
                   "2. LIBRARY DETAILS\n"
                   "Instrument ID:                  %12%\n"
                   "Run number:                     %13%\n"
                   "Flowcell ID:                    %14%\n"
                   "Lane:                           %15%\n\n"
                   "3. LIBRARY FRACTIONS\n"
                   "Sample reads; fraction:         %16%\n"
                   "Sequin reads; fraction:         %17%\n"
                   "Vector reads; fraction:         %18%\n"
                   "Total reads:                    %19%\n"
                   "Sequin dilution:                %20%\n\n"
                   "4. CALIBRATION SUMMARY\n"
                   "Sequin calibration:             %21%\n"
                   "Sequin calibration factor:      %22%\n"
                   "Sequin reads after calibration: %23%\n"
                   "Total reads after calibration:  %24%\n\n"
                   "5. ISOFORM QUANTIFICATION\n"
                   "Slope:                          %25%\n"
                   "R2:                             %26%\n"
                   "Isoform expression table:       %27%\n\n"
                   "6. GENE QUANTIFICATION\n"
                   "Slope:                          %28%\n"
                   "R2:                             %29%\n"
                   "Gene expression table:          %30%";
    
    const auto p1 = o.bam ? (CommonResults *) &stats.B1 : (CommonResults *) &stats.S1;
    const auto p2 = o.bam ? (CommonResults *) &stats.B2 : (CommonResults *) &stats.S2;

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % date()                             // 1
                                      % o.version                          // 2
                                      % o.cmd                              // 3
                                      % mixToStr(o.mix)                    // 4
                                      % o.index                            // 5
                                      % Standard::instance().rna.r1()->src // 6
                                      % (o.work + "/rna_sample*")          // 7
                                      % (o.work + "/rna_sequin*")          // 8
                                      % (o.seqC == NO_CALIBRATION ? MISSING : o.work + "/rna_sequin_calibrated*")       // 9
                                      % (o.bam ? MISSING : S0(o.k))        // 10
                                      % (o.bam ? MISSING : S2(o.rule))     // 11
                                      % p1->lib().inst(p1->lib().format()) // 12
                                      % p1->lib().run(p1->lib().format())  // 13
                                      % p1->lib().flow(p1->lib().format()) // 14
                                      % p1->lib().lane(p1->lib().format()) // 15
                                      % (S0(p1->binN(ES)) + " ; " + S2(p1->binP(ES)))
                                      % (S0(p1->total() - p1->binN(ES)) + " ; " + S2(1.0 - p1->binP(ES)))
                                      % ("0 ; 0.00") // RNA has no vector currently...
                                      % p1->total()                        // 19
                                      % S2(p1->dil())                      // 20
                                      % calib2str(o.seqC)                  // 21
                                      % (o.seqC == NO_CALIBRATION ? MISSING : S2(p1->calibF(Calibration::Sequin)))
                                      % (o.seqC == NO_CALIBRATION ? MISSING : S0(p2->total()))
                                      % (o.seqC == NO_CALIBRATION ? MISSING : S0(p1->binN(ES) + p2->total()))
                                      % (o.bam ? MISSING : S2(l1.m))            // 25
                                      % (o.bam ? MISSING : S2(l1.R2))           // 26
                                      % (o.bam ? MISSING : o.work + "/rna_sequin_isoform_table.tsv") // 27
                                      % S2(l2.m)                                // 28
                                      % S2(l2.R2)                               // 29
                                      % (o.work + "/rna_sequin_gene_table.tsv") // 30
                     ).str());
    o.writer->close();
}

static void writeQuin(const FileName &file, CommonResults *stats, const Options &o)
{
    const auto l1 = Standard::instance().rna.l1(); assert(l1); // Sequin mixture
    const auto l2 = Standard::instance().rna.l2(); assert(l2); // Gene mixture
    const auto l3 = Standard::instance().rna.l3(); assert(l3); // Sequin length

    /*
     * While it's possible deriving abundance on each individual sequin transcript, alignments don't
     * give us the information. Reads from "R1_22_1" and "R1_22_2" will both aligned to the same region.
     */

    // Sequins (R1_22_1) for FASTQ or standard (R1_22) otherwise
    std::vector<Label> names, genes;
    
    // Lengths
    std::vector<Base> lens;
    
    // Mixture
    std::vector<double> mixs;
    
    for (const auto &i : l1->m1)
    {
        if (!o.bam)
        {
            names.push_back(i.first);
            genes.push_back(RSeq2Std(i.first));
            lens.push_back(l3->input(i.first));
            mixs.push_back(l1->input(i.first, o.mix));
        }
        else
        {
            const auto std = RSeq2Std(i.first);
            
            if (std::find(names.begin(), names.end(), std) == names.end())
            {
                names.push_back(std);
                genes.push_back(std);
                lens.push_back(l3->input(i.first));
                mixs.push_back(l2->input(std, o.mix));
            }
        }
    }

    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % "NAME"
                                      % "GENE"
                                      % "LENGTH"
                                      % "MIX"
                                      % "READ"
                                      % "TPM").str());
    
    // "per million" scaling factor
    auto scale = 0.0;

    /*
     * Sum up all the reads to derive scaling factor
     */
    
    for (auto i = 0; i < names.size(); i++)
    {
        // Raw read count
        const auto read = stats->count(names[i]);
        
        // Reads per kilobase (RPK)
        const auto rpk = (Proportion) read / lens[i];

        scale += rpk;
    }
    
    // Divide by a million
    scale /= 1000000;
    
    for (auto i = 0; i < names.size(); i++)
    {
        const auto n = stats->count(names[i]);
        
        // Reads per kilobase (RPK)
        const auto rpk = (double) n / lens[i];
        
        // Transcripts per million
        const auto tpm = scale ? S2(rpk / scale) : MISSING;

        o.writer->write((boost::format(f) % names[i]
                                          % genes[i]
                                          % lens[i]
                                          % mixs[i]
                                          % n
                                          % tpm).str());
    }
    
    o.writer->close();
}

void RSplit::report(const FileName &f1, const FileName &f2, const Options &o)
{
    assert(Standard::instance().rna.l1());
    assert(Standard::instance().rna.l2());
    
    auto stats = analyze(f1, f2, o);
    
    // Generating rna_sequin.tsv
    writeQuin("rna_sequin.tsv", o.bam ? (CommonResults *) &stats.B1 : (CommonResults *) &stats.S1, o);
    
    if (o.seqC != NO_CALIBRATION)
    {
        // Generating rna_sequin_calibrated.tsv
        writeQuin("rna_sequin_calibrated.tsv", o.bam ? (CommonResults *) &stats.B2 : (CommonResults *) &stats.S2, o);
    }

    // Generating rna_reads.tsv
    if (!o.bam) { SWriteReads(Product::RNA, "rna_reads.tsv", stats.S1, o); }

    const auto tsv = o.isSCalib() ? o.work + "/rna_sequin_calibrated.tsv" : o.work + "/rna_sequin.tsv";
    
    // Generating rna_report.txt
    writeSummary("rna_report.txt", f1, f2, stats, tsv, o);
}
