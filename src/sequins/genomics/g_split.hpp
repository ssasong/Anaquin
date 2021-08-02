#ifndef G_SPLIT_HPP
#define G_SPLIT_HPP

#include "sequins/decoy_analyzer.hpp"
#include "sequins/genomics/genomics.hpp"

namespace Anaquin
{
    struct GSplit
    {
        struct Options : public SOptions
        {
            Options() : SOptions() { prod = Product::Genomics; }
        };

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

            FileName tsvE; // Only for BAM
            FileName tsvR; // Only for BAM
            FileName tsvF; // Only for BAM
            FileName tsvV; // Only for BAM
            
            // Variant results from FASTQ or BAM (before calibration)
            GVariantResults V1(const AnalyzerOptions &) const;
        };
        
        static void writeV(const FileName &, const Stats &, const SOptions &);
        
        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void  report (const FileName &, const FileName &, const Options &o = Options());
    };
}

#endif
