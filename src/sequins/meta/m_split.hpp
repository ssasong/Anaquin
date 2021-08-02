#ifndef M_SPLIT_HPP
#define M_SPLIT_HPP

#include "data/split.hpp"
#include "sequins/decoy_analyzer.hpp"

namespace Anaquin
{
    struct MSplit
    {
        struct Stats
        {
            // FASTQ results (before calibration)
            FQStats S1;
            
            // FASTQ results (after sequin calibration)
            FQStats S2;
            
            // FASTQ results (after ladder calibration)
            FQStats S3;
            
            // BAM results (before calibration)
            DecoyAnalyzer::Results B1;
            
            // BAM results (after sequin calibration)
            DecoyAnalyzer::Results B2;            

            // BAM results (after ladder calibration)
            DecoyAnalyzer::Results B3;

            FileName tsvE;  // Only for BAM
            FileName tsvR;  // Only for BAM
            FileName tsvF;  // Only for BAM
            FileName tsvL1; // Only for BAM (before ladder calibration)
            FileName tsvL2; // Only for BAM (after ladder calibration)
        };

        struct Options : public SOptions
        {
            Options()
            {
                prod = Product::Meta;
            }
            
            Base edge;
            
            Mixture mix;
        };
        
        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
