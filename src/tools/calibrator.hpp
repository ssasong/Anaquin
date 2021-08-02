#ifndef CALIBRATOR_HPP
#define CALIBRATOR_HPP

#include "tools/random.hpp"
#include "data/analyzer.hpp"

namespace Anaquin
{
    struct SelectionCalibrator
    {
        struct Result
        {
            // Number of reads selected
            Count n = 0;
            
            // Number of paired reads selected
            Count nSel = 0;
            
            // Number of paired reads skipped
            Count nSkip = 0;
            
            FileName o1, o2;
        };
        
        /*
         * Calibrate by specifying what the reads belong to. If a read is not in the map,
         * it will be written out.
         */
        
        virtual Result calibrate(const std::map<ReadName, SequinID> &,
                                 const Selection &,
                                 const WriterOptions &) = 0;

        /*
         * Calibrate by a constant probability of selection
         */
        
        virtual Result calibrate(Probability, const WriterOptions &) = 0;
        
        // Calibrator for BAM
        static std::shared_ptr<SelectionCalibrator> createBAM(const FileName &, const FileName &);
        
        // Calibrator for FASTQ
        static std::shared_ptr<SelectionCalibrator> createFQ(const FileName &, const FileName &,
                                                    const FileName &, const FileName &);
    };
    
    struct LadderCalibrator
    {
        struct Results
        {
            // Calibrated results
            SelectionCalibrator::Result cr;
            
            // Scaling factors
            std::vector<Proportion> scales;
            
            // Probability of selection
            Selection rnd;
        };
        
        enum class Method
        {
            SampleCNV2,
            Pooling
        };
        
        struct Options : public AnalyzerOptions
        {
            Method meth = Method::SampleCNV2;
            
            struct SampleCNV2Data
            {
                // Sample coverage which the CNV=2 calibrated to
                Coverage sampC;

                FileName tsvL;
            };
            
            struct PoolData
            {
                Proportion p = NO_CALIBRATION;
                
                // Unmerged ladder TSV
                FileName tsvL;
                
                std::map<ReadName, SequinID> reads;
            };
            
            bool merged = false;
            
            // Sequin regions, can be unmapped
            BedData r1;
            
            // Data for Method::SampleCNV2
            std::shared_ptr<SampleCNV2Data> d1;

            // Data for Method::PoolData
            std::shared_ptr<PoolData> d2;
        };
        
        static void getNorms(const Options &, Results &);
        
        virtual Results calibrate(const Options &) = 0;
        
        static std::shared_ptr<LadderCalibrator> createFQ(const FileName &, const FileName &,
                                                          const FileName &, const FileName &);
        static std::shared_ptr<LadderCalibrator> createBAM(const FileName &, const FileName &);
    };}

#endif
