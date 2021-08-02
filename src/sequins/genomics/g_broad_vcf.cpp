#include "tools/tools.hpp"
#include "writers/file_writer.hpp"
#include "sequins/genomics/g_broad_vcf.hpp"

using namespace Anaquin;

typedef GBroadVCF::Options Options;

static void writeSummary(const FileName &file,
                         const FileName &src,
                         const GGerm::Stats &stats,
                         const GVariant::Options &o,
                         const GVariant::Options &o_)
{
    const auto &r = Standard::instance().gen;

    // Total sequin size for this analysis
    const auto sl = GVariant::countSizeForG(o);
    
    o.info("Total sequin size: " + std::to_string(sl));
    
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
                   "2. GERMLINE VARIANTS\n"
                   "Detected sequin variants (TP):   %11%\n"
                   "Undetected sequin variants (FN): %12%\n"
                   "Erroneous sequin variants (FP):  %13%\n"
                   "Sensitivity (SN):                %14%\n"
                   "Precision (PC):                  %15%\n"
                   "TP median depth:                 %16%\n"
                   "TP median QUAL score (TPMQ):     %17%\n"
                   "FP median depth:                 %18%\n"
                   "FP median QUAL score (FPMQ):     %19%\n"
                   "Analysed region size:            %20%\n"
                   "False positives per KB:          %21%\n\n"
                   "3. SNVs                          TP; FN; FP; SN; PC; TPMQ; FPMQ\n"
                   "Total:                           %22%\n"
                   "Homozygous:                      %23%\n"
                   "Heterozygous:                    %24%\n"
                   "High GC (>70%%):                  %25%\n"
                   "Moderate GC (30-70%%):            %26%\n"
                   "Low GC (<30%%):                   %27%\n"
                   "High confidence:                 %28%\n"
                   "Low confidence:                  %29%\n"
                   "Coding:                          %30%\n"
                   "Noncoding:                       %31%\n"
                   "ACMG genes:                      %32%\n"
                   "Pharmacogenes:                   %33%\n"
                   "Simple repeats:                  %34%\n"
                   "5-12nt homopolymer:              %35%\n"
                   ">12nt homopolymer:               %36%\n\n"
                   "4. INDELS                        TP; FN; FP; SN; PC; TPMQ; FPMQ\n"
                   "Total:                           %37%\n"
                   "Homozygous:                      %38%\n"
                   "Heterozygous:                    %39%\n"
                   "High GC (>70%%):                  %40%\n"
                   "Moderate GC (30-70%%):            %41%\n"
                   "Low GC (<30%%):                   %42%\n"
                   "High confidence:                 %43%\n"
                   "Low confidence:                  %44%\n"
                   "Coding:                          %45%\n"
                   "Noncoding:                       %46%\n"
                   "ACMG genes:                      %47%\n"
                   "Pharmacogenes:                   %48%\n"
                   "Simple repeats:                  %49%\n"
                   "5-12nt homopolymer:              %50%\n"
                   ">12nt homopolymer:               %51%";
    
    auto S = [&](const GVariant::TableRow &x)
    {
        return !x.valid ? MISSING : x.broad();
    };

    const auto tpq = GVariant::getTP(o_).qs;
    const auto tpd = GVariant::getTP(o_).ds;
    const auto fpq = GVariant::getFP(o_).qs;
    const auto fpd = GVariant::getFP(o_).ds;

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % date()
                                      % o.version
                                      % o.cmd
                                      % r.v1()->src
                                      % r.r3()->src
                                      % r.a1()->src
                                      % src
                                      % "output_files/broad_vcf_TP.vcf"
                                      % "output_files/broad_vcf_FP.vcf"
                                      % "output_files/broad_vcf_FN.vcf"
                                      % stats.ds.oc.tp()     // 11
                                      % stats.ds.oc.fn()     // 12
                                      % stats.ds.oc.fp()     // 13
                                      % S3(stats.ds.oc.sn()) // 14
                                      % S3(stats.ds.oc.pc()) // 15
                                      % (tpd.empty() ? MISSING : S2(med(tpd))) // 16
                                      % (tpq.empty() ? MISSING : S2(med(tpq))) // 17
                                      % (fpd.empty() ? MISSING : S2(med(fpd))) // 18
                                      % (fpq.empty() ? MISSING : S2(med(fpq))) // 19
                                      % sl                                     // 20
                                      % S4(1000.0 * ((Proportion) stats.ds.oc.fp() / sl)) // 21
                                      % S(GVariant::getSNV(o_))                  // 22
                                      % S(GVariant::getHom(o_, "SNP"))           // 23
                                      % S(GVariant::getHet(o_, "SNP"))           // 24
                                      % S(GVariant::getHighGC(o_, "SNP"))        // 25
                                      % S(GVariant::getModGC(o_, "SNP"))         // 26
                                      % S(GVariant::getLowGC(o_, "SNP"))         // 27
                                      % S(GVariant::getHighConf(o_, "SNP"))      // 28
                                      % S(GVariant::getLowConf(o_, "SNP"))       // 29
                                      % S(GVariant::getCode(o_, "SNP"))          // 30
                                      % S(GVariant::getNCode(o_, "SNP"))         // 31
                                      % S(GVariant::getACMG(o_, "SNP"))          // 32
                                      % S(GVariant::getPGX(o_, "SNP"))           // 33
                                      % S(GVariant::getRepeat(o_, "SNP"))        // 34
                                      % S(GVariant::getShortRepeat(o_, "SNP"))   // 35
                                      % S(GVariant::getLongRepeat(o_, "SNP"))    // 36
                                      % S(GVariant::getIndel(o_))                // 37
                                      % S(GVariant::getHom(o_, "Indel"))         // 38
                                      % S(GVariant::getHet(o_, "Indel"))         // 39
                                      % S(GVariant::getHighGC(o_, "Indel"))      // 40
                                      % S(GVariant::getModGC(o_, "Indel"))       // 41
                                      % S(GVariant::getLowGC(o_, "Indel"))       // 42
                                      % S(GVariant::getHighConf(o_, "Indel"))    // 43
                                      % S(GVariant::getLowConf(o_, "Indel"))     // 44
                                      % S(GVariant::getCode(o_, "Indel"))        // 45
                                      % S(GVariant::getNCode(o_, "Indel"))       // 46
                                      % S(GVariant::getACMG(o_, "Indel"))        // 47
                                      % S(GVariant::getPGX(o_, "Indel"))         // 48
                                      % S(GVariant::getRepeat(o_, "Indel"))      // 49
                                      % S(GVariant::getShortRepeat(o_, "Indel")) // 50
                                      % S(GVariant::getLongRepeat(o_, "Indel"))  // 51
                     ).str());
    o.writer->close();
}

void GBroadVCF::report(const FileName &file, const Options &o)
{
    auto o_  = o.debug ? o : cloneO(o);
    o_.base  = "germline";
    o_.name  = "germline";
    o_.writer = std::shared_ptr<FileWriter>(new FileWriter(o_.work));
    
    const auto stats = GGerm::analyze("", file, o_);
    
    createD(o.work + "/output_files");
    copy(o_.work + "/germline_TP.vcf", o.work + "/output_files/broad_vcf_TP.vcf");
    copy(o_.work + "/germline_FP.vcf", o.work + "/output_files/broad_vcf_FP.vcf");
    copy(o_.work + "/germline_FN.vcf", o.work + "/output_files/broad_vcf_FN.vcf");

    // Generating summary
    writeSummary("report_vcf.txt", file, stats, o, o_);
}
