#include <cmath>
#include "writers/file_writer.hpp"
#include "sequins/genomics/g_somatic.hpp"

using namespace Anaquin;

typedef GSomatic::Stats Stats;
typedef GVariant::Options Options;

template <typename T, typename O> void writeQualFilter(const T &, const O &o)
{
    o.generate("somatic_qualFilter.R");
    o.writer->open("report_files/somatic_qualFilter.R");
    o.writer->write((boost::format(PlotQualFilter()) % date()
                                                     % o.cmd
                                                     % o.work
                                                     % "somatic_sequin.tsv").str());
    o.writer->close();
}

template <typename O> void writeAllele(const O &o)
{
    o.generate("somatic_ladder.R");
    o.writer->open("report_files/somatic_ladder.R");
    o.writer->write((boost::format(PlotAllele()) % date()
                                                 % o.cmd
                                                 % o.work
                                                 % "somatic_sequin.tsv").str());
    o.writer->close();
}

static void writeTable(const FileName &src, const FileName &file, const GSomatic::Stats &stats, const GVariant::Options &o)
{
    const auto tmp1  = tmpFile();
    const auto tmp2  = tmpFile();
    const auto src1  = tmpFile(); // Filtered source including all detected sequins
    const auto src1_ = tmpFile(); // src1 but only the filtered columns
    const auto tmpM  = tmpFile(); // Temporary for means
    
    RGrep(src,  tmp1, "LABEL", "SV", false); // Remove sample variants not giving allele frequency
    RGrep(tmp1, tmp2, "LABEL", "FP", false); // Remove false positives not giving allele frequency
    RGrep(tmp2, src1, "LABEL", "FN", false); // Remove false negatives
 
    RFilterC(src1, src1_, std::set<Label> { "EXP_FREQ", "REF_DEPTH_NORMAL", "VAR_DEPTH_NORMAL", "OBS_FREQ_NORMAL", "REF_DEPTH_TUMOR", "VAR_DEPTH_TUMOR", "OBS_FREQ_TUMOR", "QUAL" }, true);
 
    // Aggregated mean from all detected sequins
    RAggregateMean(src1_, tmpM, "EXP_FREQ", -1, Imputation::Remove);
    
    auto read = [&](const Label &l, std::map<std::string, std::string> &m, std::map<std::string, std::string> &s)
    {
        // Average for each allele frequency group (all detected sequins)
        RFilterC(tmpM, tmp1, std::set<Label> { "EXP_FREQ", l}, true); m = RReadTSV(tmp1, "EXP_FREQ", l);

        // Standard deviation for each allele frequency group
        RFilterC(tmpM, tmp1, std::set<Label> { "EXP_FREQ", l}, true); s = RReadTSV(tmp1, "EXP_FREQ", l);
    };
    
    std::map<std::string, std::string> m1, m3, m4, m5, m6;
    std::map<std::string, std::string> s1, s3, s4, s5, s6;

    read("REF_DEPTH_NORMAL", m1, s1);
    read("REF_DEPTH_TUMOR",  m3, s3);
    read("VAR_DEPTH_TUMOR",  m4, s4);
    read("OBS_FREQ_TUMOR",   m5, s5);
    read("QUAL",             m6, s6);

    o.generate("somatic_variant_table.tsv");
    o.writer->open("somatic_variant_table.tsv");
    o.writer->write("EXP_FREQ\tOBS_FREQ\tTP\tFN\tREF_DEPTH_TUMOR\tVAR_DEPTH_TUMOR\tDEPTH_NORMAL\tQUAL_SCORE");

    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";
    
    // For each allele frequency group ...
    for (auto &i : m1)
    {
        const auto key = i.first;
        
        auto ms = [&](const std::map<std::string, std::string> &m, const std::map<std::string, std::string> &s)
        {
            if (!m.count(key))
            {
                return std::string("NA");
            }
            else
            {
                return m.at(key) + " +- " + (s.count(i.first) ? s.at(i.first) : "NA");
            }
        };
        
        auto count = [&](const Label &x)
        {
            const auto tmp1 = tmpFile();
            const auto tmp2 = tmpFile();
            const auto tmp3 = tmpFile();

            RGrep(src,  tmp1, "LABEL", "SV", false); // Remove sample variants not giving allele frequency
            RGrep(tmp1, tmp2, "LABEL", "FP", false); // Remove false positives not giving allele frequency
            RGrep(tmp2, tmp3, "EXP_FREQ", key, true, true);
            
            return RCount(tmp3, "LABEL", x);
        };
        
        o.writer->write((boost::format(f) % key
                                          % ms(m5, s5)
                                          % count("TP")
                                          % count("FN")
                                          % ms(m3, s3)
                                          % ms(m4, s4)
                                          % ms(m1, s1)
                                          % ms(m6, s6)).str());
    }
    
    o.writer->close();
}

template <typename O> void writeROC(const FileName &file, const O &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotStrelkaROC()) % date()
                                                     % o.cmd
                                                     % o.work
                                                     % "somatic_sequin.tsv").str());
    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &f1,
                         const FileName &f2,
                         const Stats &stats,
                         const GVariant::Options &o)
{
    const auto &r = Standard::instance().gen;
    const auto l1 = stats.oa.linear();
    
    // Total sequin size for this analysis
    const auto sl = GVariant::countSizeForG(o);
    
    const auto f = "1. ANALYSIS\n"
                   "Date:                            %1%\n"
                   "Anaquin version:                 %2%\n"
                   "Command:                         %3%\n"
                   "Reference variants:              %4%\n"
                   "Reference regions:               %5%\n"
                   "Reference stratification:        %6%\n"
                   "Input variant file:              %7%\n"
                   "TP sequin variants:              %8%\n"
                   "FP sequin variants:              %9%\n"
                   "FN sequin variants:              %10%\n\n"
                   "2. SOMATIC VARIANTS\n"
                   "Detected sequin variants (TP):   %11%\n"
                   "Undetected sequin variants (FN): %12%\n"
                   "Erroneous sequin variants (FP):  %13%\n"
                   "Sensitivity (SN):                %14%\n"
                   "Precision (PC):                  %15%\n"
                   "Analysed region size:            %16%\n"
                   "FP per KB:                       %17%\n"
                   "Median FP observed frequency:    %18%\n"
                   "Median FP depth (AP):            %19%\n"
                   "Median FP quality (QUAL):        %20%\n\n"
                   "3. VARIANT ALLELE FREQUENCY\n"
                   "Slope:                           %21%\n"
                   "R2:                              %22%\n"
                   "Sequin table:                    %23%";

    auto S = [&](const GVariant::TableRow &x)
    {
        return !x.valid ? MISSING : x.broad();
    };

    const auto fpq  = GVariant::getFP(o).qs;
    const auto fpd  = GVariant::getFP(o).ds;
    const auto fpa  = GVariant::getFP(o).af;
    const auto size = GVariant::countSizeForG(o);
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % date()
                                      % o.version
                                      % o.cmd
                                      % r.v1()->src
                                      % r.r3()->src
                                      % r.a1()->src
                                      % (f1 + (f1.empty() ? "" : " and " + f2))
                                      % (o.work + "/sequins_TP.vcf")
                                      % (o.work + "/sequins_FP.vcf")
                                      % (o.work + "/sequins_FN.vcf")
                                      % stats.ds.oc.tp()     // 11
                                      % stats.ds.oc.fn()     // 12
                                      % stats.ds.oc.fp()     // 13
                                      % S3(stats.ds.oc.sn()) // 14
                                      % S3(stats.ds.oc.pc()) // 15
                                      % size                 // 16
                                      % S2(1000.0 * ((double) stats.ds.oc.fp() / size))
                                      % (fpa.empty() ? MISSING : S2(med(fpa))) // 18
                                      % (fpd.empty() ? MISSING : S2(med(fpd))) // 19
                                      % (fpq.empty() ? MISSING : S2(med(fpq))) // 20
                                      % l1.m
                                      % l1.R2
                                      % (o.work + "/somatic_variant_table.tsv")).str());
    o.writer->close();
}

void GSomatic::report(const FileName &f1, const FileName &f2, const Options &o)
{
    const auto &r = Standard::instance().gen;

    auto o1 = cloneO(o);
    o1.es = tmpFile();
    o1.tp = o.work + "/sequins_TP.vcf";
    o1.fp = o.work + "/sequins_FP.vcf";
    o1.fn = o.work + "/sequins_FN.vcf";
    o1.tsvS = tmpFile();

    Stats stats;
    GVariant::analyze(stats, GSomatic(), f1, f2, o1, Standard::instance().gen.v3());

    for (const auto &m : stats.ds.tps)
    {
        auto oAF = [&]()
        {
            /*
             * Since this tool works with both germline and somatic callers.
             * We need to guess what the observed allele frequency.
             */
            
            // Tumor allele fequency?
            const auto tAF = tumorAF(m.qry);
            
            return std::isnan(tAF) ? m.qry.allF : tAF;
        };

        const auto exp = r.af(m.var->name);
        const auto obs = oAF();
        
        assert(!std::isnan(exp));
        
        // Eg: 2821292107
        const auto x = S2(m.var->key());
        
        // Add for all variants
        stats.oa.add(x, exp, obs);
        
        // Add for mutation type
        stats.m2a[m.qry.type()].add(x, exp, obs);
    }
    
    // Determine quantification limits
    stats.oa.limit = stats.oa.limitQuant();
    
    for (auto &i : stats.m2a)
    {
        i.second.limit = i.second.limitQuant();
    }
    
    /*
     * Performance by allele frequency
     */
    
    for (const auto &i : r.v1()->data.vars())
    {
        if (isSoma(i.name))
        {
            const auto af = r.af(i.name);
            stats.af[af].nr()++;
            
            if (stats.ds.findTP(i.name))
            {
                stats.af[af].tp()++;
            }
            else
            {
                stats.af[af].fn()++;
            }
        }
    }

    auto o2 = o.debug ? o : cloneO(o);
    o2.tsvS = o2.work + "/somatic_sequin.tsv";
    o2.writer = std::shared_ptr<Writer<>>(new FileWriter(o2.work));

    // Generating somatic_sequin.tsv
    writeSomatic("somatic_sequin.tsv", stats, o2);

    // Generating somatic_variant_table.tsv
    writeTable(o2.work + "/somatic_sequin.tsv", "somatic_variant_table.tsv", stats, o);

    auto o3 = o;
    o3.tsvS = o2.tsvS;

    // Generating somatic_report.txt
    writeSummary("somatic_report.txt", f1, f2, stats, o3);

    if (o.debug)
    {
        // Generating somatic_ROC.R
        writeROC("report_files/somatic_ROC.R", o);
        
        // Generating somatic_qualFilter.R
        writeQualFilter(stats, o);
        
        // Generating somatic_ladder.R
        writeAllele(o);
    }
}
