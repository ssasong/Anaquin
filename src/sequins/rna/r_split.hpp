#ifndef R_SPLIT_HPP
#define R_SPLIT_HPP

#include "data/split.hpp"
#include "sequins/decoy_analyzer.hpp"

namespace Anaquin
{
    struct RSplit
    {
        struct Options : public SOptions
        {
            Options()
            {
                prod = Product::RNA;
            }
            
            Mixture mix;
        };
        
        struct Stats
        {
            // (FASTQ) Before calibration
            FQStats S1;
            
            // (FASTQ) After calibration
            FQStats S2;
            
            // (BAM) Alignment before calibration
            DecoyAnalyzer::Results B1;
            
            // (BAM) Alignment after calibration
            DecoyAnalyzer::Results B2;

            // Isoform and gene-level expression
            SequinStats l1, l2;
        };

        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
