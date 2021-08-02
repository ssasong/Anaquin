#ifndef G_CALIBRATE_HPP
#define G_CALIBRATE_HPP

#include "sequins/genomics/g_broad_bam.hpp"

namespace Anaquin
{
    struct GCalibrate
    {
        typedef CalibrateMethod Method;
        
        struct Stats
        {
            // Empty implementation
        };
        
        enum class CalibrateMode
        {
            TwoBAM,
            Combined
        };

        struct Options : SOptions
        {
            Options() {}

            // How to calculate coverage? Defined only if seqC == NO_CALIBRATION
            Method meth = Method::Mean;
            
            // Only for Method::Custom
            double customSequinThreshold = 0;

            Base edge;

            // Calibrate to sampling if not defined
            Proportion seqC = NO_CALIBRATION;
            
            // Defined only if meth==Reads
            Count reads = 0;
            
            // The specialized cancer tool?
            bool isCancer = false;
            
            bool writeS = false;
            bool writeD = false;
            bool writeC = false;
        };

        static Stats analyze(const FileName &, const FileName &, const Options &);

        static void report(const FileName &, const Options &o = Options());
        static void report(const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
