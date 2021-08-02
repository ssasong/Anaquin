#include "tools/tools.hpp"
#include "parsers/parser_fa.hpp"
#include "parsers/parser_csv.hpp"
#include "writers/file_writer.hpp"
#include "sequins/genomics/g_broad_bam.hpp"

using namespace Anaquin;

GBroadBam::SomaticReport GBroadBam::reportS(const FileName &file)
{
    GBroadBam::SomaticReport r;
    
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    const auto tmp3 = tmpFile();
    
    // Construct log-linear model
    RGrep(file, tmp1, "LABEL", "Somatic"); r.lm = RLinear(tmp1, "NAME", "EXP_FREQ", "OBS_FREQ_CALIB").linear(false);
    
    const auto exps = RReadTSV(tmp1, "EXP_FREQ", "EXP_FREQ");
    
    // TODO: hg38 calibration not supported for now
    if (exps.empty())
    {
        return r;
    }
    
    assert(!exps.empty());
    
    struct Data
    {
        // Alignment Count for reference and alternative Count
        std::vector<Count> ref, var;
        
        // Observed allele frequency
        std::vector<Proportion> obs;
    };
    
    std::map<double, Data> m;
    
    // For each allele frequency group...
    for (const auto &exp : exps)
    {
        RGrep(tmp1, tmp2, "EXP_FREQ", exp.first);
        m[stod(exp.first)].ref = numeric<Count>(RList(tmp2, "REF_COUNT_CALIB"));
        m[stod(exp.first)].var = numeric<Count>(RList(tmp2, "VAR_COUNT_CALIB"));
        m[stod(exp.first)].obs = numeric<Proportion>(RList(tmp2, "OBS_FREQ_CALIB"));
    }
    
    std::stringstream ss;
    
    for (auto i = m.rbegin(); i != m.rend(); i++)
    {
        ss << (i != m.rbegin() ? "\n" : "")
           << std::to_string(i->first) << "\t"
           << (!i->second.obs.empty() ? S4(SS::mean(i->second.obs)) : "NA") << " ; "
           << (!i->second.ref.empty() ? S1(SS::mean(i->second.ref)) : "NA") << " ; "
           << (!i->second.var.empty() ? S1(SS::mean(i->second.var)) : "NA");
    }
    
    r.table = ss.str();
    return r;
}

GBroadBam::MixtureReport GBroadBam::reportM(const FileName &file, std::shared_ptr<Translation> t)
{
    GBroadBam::MixtureReport r;
    std::map<SequinID, Coverage> m;

    const auto x = RReadTSV(file, "NAME", "PRE_COVERAGE");
    
    for (auto &i : *(t))
    {
        if (x.count(i.first) && x.at(i.first) != MISSING)
        {
            m[i.first] = stod(x.at(i.first));
        }
    }
    
    typedef std::pair<SequinID, Coverage> Pair;
    
    const auto max = std::max_element(m.begin(), m.end(), [](const Pair &p1, const Pair &p2)
    {
        return p2.second > p1.second;
    });

    r.text = m.empty() ? MISSING : (*t)[max->first];
    return r;
}

GBroadBam::SyntheticReport GBroadBam::reportL(const FileName &file)
{
    GBroadBam::SyntheticReport r;

    const auto tmpA = tmpFile();
    const auto tmpB = tmpFile();
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    const auto tmp4 = tmpFile();
    const auto tmp8 = tmpFile();

    RFilterC(file, tmpA, std::set<Label> { "NAME", "STANDARD" }, false);
    RAggregateSum(tmpA, tmpB, "SEQUIN");

    r.lm = RLinear(tmpB, "SEQUIN", "COPY", "MEAN").linear();

    r.text[0] = S4(r.lm.m);
    r.text[1] = S4(r.lm.R2);
    
    RGrep(tmpB, tmpA, "COPY", "1"); RGrep(tmpA, tmp1, "COPY", "16", false);
    RGrep(tmpB, tmpA, "COPY", "4"); RGrep(tmpA, tmp2, "COPY", "64", false);
    RGrep(tmpB, tmp4, "COPY", "16"); // CP = 4
    RGrep(tmpB, tmp8, "COPY", "64"); // CP = 8

    const auto x1 = numeric<Coverage>(RReadTSV(tmp1, "MEAN"));
    const auto x2 = numeric<Coverage>(RReadTSV(tmp2, "MEAN"));
    const auto x4 = numeric<Coverage>(RReadTSV(tmp4, "MEAN"));
    const auto x8 = numeric<Coverage>(RReadTSV(tmp8, "MEAN"));

    r.text[2] = S4(SS::mean(x1)) + " ; " + S4(SS::CV(x1));
    r.text[3] = S4(SS::mean(x2)) + " ; " + S4(SS::CV(x2));
    r.text[4] = S4(SS::mean(x4)) + " ; " + S4(SS::CV(x4));
    r.text[5] = S4(SS::mean(x8)) + " ; " + S4(SS::CV(x8));
    
    const auto mean12 = (!x1.empty() && !x2.empty()) ? SS::mean(x2) / SS::mean(x1) : NAN;
    const auto mean24 = (!x2.empty() && !x4.empty()) ? SS::mean(x4) / SS::mean(x2) : NAN;
    const auto mean48 = (!x4.empty() && !x8.empty()) ? SS::mean(x8) / SS::mean(x4) : NAN;

    const auto meanR = removeNA(std::vector<double> { mean12, mean24, mean48 });
    r.text[6] = !meanR.empty() ? S4(SS::mean(meanR)) : MISSING;

    return r;
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

GBroadBam::ErrorReport GBroadBam::reportE(const FileName &tsvE, const FileName &tsvR, const FileName &tsvF, bool onlyC, bool isChrQ)
{
    GBroadBam::ErrorReport r;

    auto tmp1 = tmpFile();
    auto tmp2 = tmpFile();
    auto tmp3 = tmpFile();

    // Keep only mismatches
    RGrep(tsvE, tmp1, "LABEL", "Match", false);
    
    // Only calibrated (that's what sequins can inform about the sample)
    if (onlyC) { RGrep(tmp1, tmp2, "CALIBRATED", "true"); } else { tmp2 = tmp1; }

    // Only chrQS is being used for error calculation
    if (isChrQ) { RGrep(tsvR, tmp3, "CHROM", "chrQS"); } else { tmp3 = tsvR; }

    struct Data
    {
        SequinID name;
        Label label;
        std::string data1, data2, data3;
    };
    
    // Data structure for all mutations
    std::map<ChrID, std::map<Base, Data>> m;
    
    ParserCSV::parse(Reader(tmp2), [&](const ParserCSV::Data &x, Progress i)
    {
        if (!i) { return; }
        
        assert(x[0] == "true" || x[0] == "false");
        assert(x[4] == "All");

        Data d;
        d.name  = x[2];
        d.label = x[3];
        d.data1 = x[5];
        d.data2 = x[6];
        d.data3 = x[7];

        m[x[1]][stod(x[5])] = d;
    }, "\t");

    // Mean read length
    const auto rl = RMean(tmp3, "MEAN_READ_LENGTH");

    // Total number of reads in all chrQS regions
    const auto allReads = RSum(tmp3, onlyC ? "POST_READ" : "PRE_READ");

    // Number of errors, number of alignments and length for a stratification
    auto f1 = [&](Count x, Count rl, Base reads, Base size)
    {
        // (Number of errors / number of reads) / (mean read length)
        const auto base = S6(((Proportion) x / reads) / rl);

        // COUNT ; PER_BASE ; REGION_SIZE
        return S0(x) + " ; " + base + " ; " + S4(size);
    };
    
    // Count SNP, insertions and deletions from tsvE
    auto f2 = [&](const ChrID &cID, Base start, Base end)
    {
        SearchResults r;
        
        // Do we have any mutations in the chromosome?
        if (!m.count(cID))
        {
            return r;
        }

        for (auto i = start; i < end; i++)
        {
            if (m[cID].count(i))
            {
                if (m[cID][i].label == "SNP")
                {
                    const auto tv = std::set<std::string> { "AG", "GA", "CT", "TC" };

                    // Transversion?
                    if (tv.count(m[cID][i].data2))
                    {
                        r.tv += stod(m[cID][i].data3);
                    }
                    else
                    {
                        r.ts += stod(m[cID][i].data3);
                    }
                }
                else if (m[cID][i].label == "Insertion")
                {
                    r.in += stod(m[cID][i].data3);
                }
                else
                {
                    r.dl += stod(m[cID][i].data3);
                }
            }
        }
        
        return r;
    };

    auto f3 = [&](const Label &name, bool exact = true)
    {
        const auto tmp = tmpFile();

        // Grep to only the stratificaion
        RGrep(tsvF, tmp, "NAME", name, true, false, exact);
     
        // Total size
        const auto size = (Base) sum(RSubtract(tmp, "START", "END"));
        
        auto chr = RReadTSV(tmp, "CHROM");
        auto str = numeric<Base>(RReadTSV(tmp, "START"));
        auto end = numeric<Base>(RReadTSV(tmp, "END"));
        assert(chr.size() == str.size());
        assert(str.size() == end.size());

        SearchResults rr;
        
        for (auto i = 0; i < str.size(); i++)
        {
            rr += f2(chr[i], str[i], end[i]);
        }
        
        return f1(rr.total(), rl, RSum(tmp, onlyC ? "POST_READ" : "PRE_READ"), size);
    };
    
    // Total size in chrQS
    auto size = (Base) sum(RSubtract(tmp3, "START", "END"));
    
    // Edge for each region in chrQS
    const auto edge = numeric<Base>(RReadTSV(tmp3, "EDGE"));

    // Assume edges are all equal (they must)
    size -= !edge.empty() ? edge.size() * (2 * edge.front()) : 0;
    
    // Total number of sequencing errors on calibrated sample
    r.text[0] = f1(RCount(tmp2, "CHROM", isChrQ ? "chrQS" : "chr", false), (Count) rl, allReads, size);
    
    // Total number of mismatches on calibrated sample
    r.text[1] = f1(RCount(tmp2, "LABEL", "SNP"), (Count) rl, allReads, size);

    // Total number of insertions on calibrated sample
    r.text[3] = f1(RCount(tmp2, "LABEL", "Insertion"), (Count) rl, allReads, size);
    
    // Total number of deletions on calibrated sample
    r.text[4] = f1(RCount(tmp2, "LABEL", "Deletion"), (Count) rl, allReads, size);

    r.text[5] = f3("GeneContext_CodingRegion", false);      // All coding including ACGM and PGx
    r.text[6] = f3("GeneContext_NoncodingRegion", false);   // Noncoding
    r.text[7] = f3("GeneContext_CodingRegion_ACMG", false); // ACMG
    r.text[8] = f3("GeneContext_CodingRegion_PGx", false);  // Pharmacogenes
    
    if (m.empty())
    {
        r.text[2] = MISSING;
    }
    else
    {
        SearchResults x;
        
        // Add up all regions in all chromosomes
        for (auto &i : m)
        {
            x += f2(i.first, i.second.begin()->first, i.second.rbegin()->first);
        }
        
        // Transition to Transversion ratio
        r.text[2] = S4((Proportion) x.ts / x.tv);
    }

    return r;
}

GBroadBam::CoverageReport GBroadBam::reportC(const FileName &file, bool isChrQ)
{
    GBroadBam::CoverageReport r;
    
    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    auto tmp3 = tmpFile();

    RGrep(file, tmp1, "NAME", "All");
    RGrep(file, tmp2, "NAME", "All", false);
    
    // Making it "chrQ" will also read in "chrQL". Let's assume "chrQR" is not causing problems here...
    if (isChrQ) { RGrep(tmp2, tmp3, "CHROM", "chrQS"); } else { tmp3 = tmp2; }
    
    const auto x4 = numeric<Coverage>(RReadTSV(tmp3, "SAMPLE_COVERAGE"));
    const auto x5 = numeric<Coverage>(RReadTSV(tmp3, "PRE_COVERAGE"));
    const auto x6 = numeric<Coverage>(RReadTSV(tmp3, "POST_COVERAGE"));
    const auto x7 = numeric<Coverage>(RReadTSV(tmp1, "SAMPLE_READ")); // All sample reads
    const auto x8 = numeric<Coverage>(RReadTSV(tmp1, "PRE_READ"));    // All sequin reads before
    const auto x9 = numeric<Coverage>(RReadTSV(tmp1, "PRE_READ"));    // All sequin reads after

    const auto s1  = x7[0];   // All sample reads
    const auto s2  = x8[0];   // All sequin reads before calibration
    const auto s3  = x9[0];   // All sequin reads after calibration
    const auto s12 = s1 + s2; // Total before calibration
    const auto s13 = s1 + s3; // Total after calibration
    const auto dil = ((Proportion) s2 / s12); // Before calibration

    // After calibration
    RGrep(file, tmp1, "CHROM", "chrQL"); const auto QL = numeric<Coverage>(RReadTSV(tmp1, "POST_READ"));
    
    // After calibration
    RGrep(file, tmp2, "CHROM", "chrQV"); const auto QV = numeric<Coverage>(RReadTSV(tmp2, "POST_READ"));

    r.text[0] = S0(s12);
    r.text[1] = S0(s1) + " ; " + S4((Proportion) s1 / s12);
    r.text[2] = S4(SS::mean(x4)) + " ; " + S4(SS::CV(x4));
    r.text[3] = S0(s2) + " ; " + S4((Proportion) s2 / s12);
    r.text[4] = S4(SS::mean(x5)) + " ; " + S4(SS::CV(x5));
    r.text[5] = S4(SS::mean(x6)) + " ; " + S4(SS::CV(x6));
    
    // Synthetic ladder alignments after calibraation
    r.text[6] = S0(SS::sum(QL)) + " ; " + S4((Proportion) SS::sum(QL) / s13);
    
    // Vector alignments after calibration
    r.text[7] = S0(SS::sum(QV)) + " ; " + S4((Proportion) SS::sum(QV) / s13);

    // Keep it as fraction
    r.text[8] = S4(dil);

    return r;
}

static void writeSummary(const FileName &file,
                         const FileName &f1,
                         const FileName &f2,
                         const GDecoyResults &rr,
                         const GBroadBam::Options &o1,
                         const GDecoyOptions &o2,
                         const FileName &origW)
{
    const auto tmp = tmpFile();
    const auto &r = Standard::instance().gen;
    const auto isChrQ = f2.empty();
    
    const auto f = "SEQUIN CALIBRATE REPORT:\n\n"
                   "1. ANALYSIS\n"
                   "Date:                                        %1%\n"
                   "Anaquin version:                             %2%\n"
                   "Command:                                     %3%\n"
                   "Sequin mixture:                              %4%\n"
                   "Reference index:                             %5%\n"
                   "Reference regions:                           %6%\n"
                   "Calibration regions:                         %7%\n"
                   "Input alignment file:                        %8%\n"
                   "Sequin reads path:                           %9%\n"
                   "Calibrated reads path:                       %10%\n"
                   "Sample & sequin reads path:                  %11%\n\n"
                   "2. LIBRARY DETAILS\n"
                   "Instrument ID:                               %12%\n"
                   "Run number:                                  %13%\n"
                   "Flowcell ID:                                 %14%\n"
                   "Lane:                                        %15%\n\n"
                   "3. SEQUENCING ERRORS                         COUNT; PER_BASE; REGION_SIZE\n"
                   "Total:                                       %16%\n"
                   "Mismatch:                                    %17%\n"
                   "Transitions/Transversions:                   %18%\n"
                   "Insertion:                                   %19%\n"
                   "Deletion:                                    %20%\n"
                   "Coding:                                      %21%\n"
                   "Noncoding:                                   %22%\n"
                   "ACMG genes:                                  %23%\n"
                   "Pharmacogenes                                %24%\n\n"
                   "4. SEQUIN CALIBRATION\n"
                   "Total alignments:                            %25%\n"
                   "Sample alignments; fraction:                 %26%\n"
                   "Sample mean coverage; CV:                    %27%\n"
                   "Sequin alignments; fraction:                 %28%\n"
                   "Sequin mean coverage before calibration; CV: %29%\n"
                   "Sequin mean coverage after calibration; CV:  %30%\n"
                   "Synthetic ladder alignments; fraction:       %31%\n"
                   "Vector alignments; fraction:                 %32%\n"
                   "Sequin dilution:                             %33%\n\n"
                   "5. COPY-NUMBER VARIATION\n"
                   "Slope:                                       %34%\n"
                   "R2:                                          %35%\n"
                   "1n (deletion) coverage; CV:                  %36%\n"
                   "2n (diploid) coverage; CV:                   %37%\n"
                   "4n (amplification) coverage; CV:             %38%\n"
                   "8n (amplification) coverage; CV:             %39%\n"
                   "Mean ratio between 1,2,4,8n:                 %40%\n\n"
                   "6. SOMATIC VARIANTS\n"
                   "Quantitative accuracy (slope):               %41%\n"
                   "Quantitative accuracy (correlation):         %42%\n"
                   "Mutation allele frequency (AF; REF; VAR)\n"
                   "%43%";
    
    const auto S = GBroadBam::reportS(o2.work + "/" + o2.tsvA);
    const auto C = GBroadBam::reportC(o2.work + "/" + o2.tsvR, isChrQ);
    const auto L = GBroadBam::reportL(o2.work + "/" + o2.tsvL2);
    const auto M = GBroadBam::reportM(o2.work + "/" + o2.tsvR, r.t1());
    
    // Only after calibration. Neat mixture won't work.
    const auto E = GBroadBam::reportE(o2.work + "/" + o2.tsvE,
                                      o2.work + "/" + o2.tsvR,
                                      o2.work + "/" + o2.tsvF,
                                      true,
                                      isChrQ);

    extern bool __HACK_IS_CANCER__;
    
    o1.generate(file);
    o1.writer->open(file);
    o1.writer->write((boost::format(f) % date()                                 // 1
                                       % o1.version                             // 2
                                       % o1.cmd                                 // 3
                                       % (__HACK_IS_CANCER__ ? "v3" : M.text)   // 4
                                       % o1.index                               // 5
                                       % r.r2()->src                            // 6
                                       % r.r1()->src                            // 7
                                       % (f2.empty() ? f1 : f1 + " and " + f2)  // 8
                                       % (origW + "/sequin.bam")                // 9
                                       % (origW + "/sequin_calibrated.bam")     // 10
                                       % (origW + "/merged.bam")                // 11
                                       % rr.lib.inst(rr.lib.format())           // 12
                                       % rr.lib.run(rr.lib.format())            // 13
                                       % rr.lib.flow(rr.lib.format())           // 14
                                       % rr.lib.lane(rr.lib.format())           // 15
                                       % E.text.at(0)      // 16
                                       % E.text.at(1)      // 17
                                       % E.text.at(2)      // 18
                                       % E.text.at(3)      // 19
                                       % E.text.at(4)      // 20
                                       % E.text.at(5)      // 21
                                       % E.text.at(6)      // 22
                                       % E.text.at(7)      // 23
                                       % E.text.at(8)      // 24
                                       % C.text.at(0)      // 25
                                       % C.text.at(1)      // 26
                                       % C.text.at(2)      // 27
                                       % C.text.at(3)      // 28
                                       % C.text.at(4)      // 29
                                       % C.text.at(5)      // 30
                                       % C.text.at(6)      // 31
                                       % C.text.at(7)      // 32
                                       % C.text.at(8)      // 33
                                       % L.text.at(0)      // 34
                                       % L.text.at(1)      // 35
                                       % L.text.at(2)      // 36
                                       % L.text.at(3)      // 37
                                       % L.text.at(4)      // 38
                                       % L.text.at(5)      // 39
                                       % L.text.at(6)      // 40
                                       % replaceNA(S.lm.m) // 41
                                       % replaceNA(S.lm.r) // 42
                                       % S.table           // 43
                     ).str());
    o1.writer->close();
}

void GBroadBam::report(const FileName &f1, const FileName &f2, const Options &o1)
{    
    GDecoyOptions o2;
    
    if (o1.debug) { o2.work = o1.work;   } // Always write to working directory if debug
    else          { o2.work = tmpPath(); }
    
    o2.meth    = o1.meth;
    o2.seqC    = o1.seqC; // (BroadBAM doesn't use seqL because it uses sampling coverage)
    o2.writer  = o1.writer;
    o2.logger  = o1.logger;
    o2.output  = o1.output;
    o2.showGen = false;
    o2.debug   = o1.debug;
    o2.meth    = o1.meth;
    o2.customSequinThreshold = o1.customSequinThreshold;

    o2.writeS  = o1.work + "/sample.bam";
    o2.writeD  = o1.work + "/sequin.bam";
    o2.writeM  = o1.work + "/merged.bam";
    o2.writeT  = o2.work + "/trimmed.bam";    
    o2.writeL1 = tmpFile();
    o2.inputL2 = o2.work + "/trimmed.bam";
    o2.writeL2 = tmpFile();
    o2.inputMC = o2.writeL2;
    o2.writeMC = o1.work + "/sequin_calibrated.bam";

    o2.tsvE  = "broadBAM_errors.tsv";
    o2.tsvR  = "broadBAM_regions.tsv";
    o2.tsvF  = "broadBAM_features.tsv";
    o2.tsvA  = "broadBAM_variants.tsv";
    o2.tsvL1 = "broadBAM_synthetic_L1.tsv";
    o2.tsvL2 = "broadBAM_synthetic.tsv";

    const auto &r = Standard::instance().gen;
    const auto isChrQ = f2.empty();
    
    o2.l1   = r.l1();
    o2.l3   = r.l3();
    o2.v1   = isChrQ ? r.v4() : r.v1();
    o2.a3   = r.a3();
    o2.edge = 550;

    o2.h1 = *(r.r1()); o2.h2 = *(r.r3());
    o2.r1 = *(r.r2()); o2.r2 = *(r.r4());
    o2.attr = *(r.r5());

    o2.index = o1.index;
    o2.writer = std::shared_ptr<Writer<>>(new FileWriter(o2.work));

    assert(!o2.r1.empty());
    assert(!o2.r2.empty());

    if (f2.empty())
    {
        o2.errors.insert(GDecoyChrQS);
    }
    else
    {
        for (auto &i : o2.r1)
        {
            o2.errors.insert(i.first);
        }
    }
    
    // Generating broad_bam.txt
    writeSummary("broad_bam.txt", f1, f2, GDecoyAnalysis(f1, f2, o2,
        [&](GDecoyStatus,
            const ParserBAM::Data &,
            const DInter *,
            const DInter *,
            bool)
        {}, [&](GDecoyStatus) {}), o1, o2, o1.origW);
}
