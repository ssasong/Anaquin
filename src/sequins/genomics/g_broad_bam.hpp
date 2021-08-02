#ifndef G_BROAD_BAM_HPP
#define G_BROAD_BAM_HPP

#include "writers/file_writer.hpp"
#include "sequins/genomics/genomics.hpp"

namespace Anaquin
{
    struct GBroadBam
    {
        struct Options : public AnalyzerOptions
        {
            Options() : seqC(NO_CALIBRATION) {}
            Options(const AnalyzerOptions &o) : AnalyzerOptions(o), seqC(NO_CALIBRATION)
            {
                writer = std::shared_ptr<Writer<>>(new FileWriter(o.work));
            }
            
            FileName origW;
            FileName index;
            Proportion seqC;
            CalibrateMethod meth = CalibrateMethod::Mean;
            
            // Only if meth == CalibrationMethod::Custom
            double customSequinThreshold = 0;
        };
        
        struct SomaticReport
        {
            LinearModel lm;
            std::string table;
        };
        
        struct ErrorReport
        {
            std::map<unsigned, std::string> text;
        };

        struct CoverageReport
        {
            std::map<unsigned, std::string> text;
        };
        
        struct SyntheticReport
        {
            LinearModel lm;
            std::map<unsigned, std::string> text;
        };
        
        struct MixtureReport
        {
            std::string text;
        };
        
        static MixtureReport reportM(const FileName &, std::shared_ptr<Translation>);
        static SyntheticReport reportL(const FileName &);
        static ErrorReport reportE(const FileName &,
                                   const FileName &,
                                   const FileName &,
                                   bool onlyC = true,
                                   bool isChrQ = true);
        static CoverageReport reportC(const FileName &, bool isChrQ = true);
        static SomaticReport reportS(const FileName &);
        
        static void report(const FileName &, const FileName &, const Options &o = Options());
        static void report(const FileName &f1, const Options &o = Options())
        {
            GBroadBam::report(f1, "", o);
        }
    };
}

#endif
