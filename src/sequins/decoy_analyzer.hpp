#ifndef DECOY_ANALYZER_HPP
#define DECOY_ANALYZER_HPP

#include "data/bData.hpp"
#include "data/library.hpp"
#include "data/analyzer.hpp"
#include "data/reference.hpp"
#include "writers/bam_writer.hpp"

namespace Anaquin
{
    struct DecoyAnalyzer
    {
        struct Options : public WriterOptions
        {
            // Sequin calibration relative to sample?
            Proportion seqC = NO_CALIBRATION;
            
            BedData h1, h2;
            BedData r1, r2, attr;
            
            std::shared_ptr<Ladder> l1, l3;
            std::shared_ptr<VCFLadder> v1;
            std::shared_ptr<AttributeBed> a3;

            Base edge;
            FileName index;
            
            // File for merged alignments
            FileName writeM;
            
            // File for sample alignments
            FileName writeS;
            
            // File for non-trimmed decoy alignments
            FileName writeD;
            
            // File for trimmed decoy alignments
            FileName writeT;
            
            // Input alignment file writeC (usually writeT or writeD)
            FileName inputC;

            // Output file for calibrated alignments
            FileName writeC;
                        
            // Chromosomes will be reported for error profiles
            std::set<ChrID> errors;
            
            Base trim = 1;
            
            bool shouldTrim  = true;
            bool shouldError = true;
            
            // Sequin sequences
            CustomMap<ChrID, Sequence> seqs;
        };
        
        // 0-based (required conversion to 1-based for reporting)
        typedef std::map<SequinID, std::map<Base, std::map<SNPType, Count>>> SData;
        
        // 0-based (required conversion to 1-based for reporting)
        typedef std::map<SequinID, std::map<Base, std::map<Base, Count>>> IDData;
        
        struct Results : public CommonResults
        {
            // 0-based (required conversion to 1-based for reporting)
            typedef std::map<SequinID, std::map<Base, Count>> MData;
            
            struct StratData
            {
                // Data for matching
                MData match;

                // Data for insertion
                SData snps;

                // Data for insertion and deletion
                IDData ins, dls;
                
                inline Count totalS(SNPType *type = nullptr) const
                {
                    Count n = 0;
                    
                    for (auto &i : snps)
                    {
                        for (auto &j : i.second)
                        {
                            for (auto &k : j.second)
                            {
                                if (!type || *type == k.first)
                                {
                                    n += k.second;
                                }
                            }
                        }
                    }
                    
                    return n;
                }
                
                inline Count totalIDData(const IDData &x) const
                {
                    Count n = 0;
                    
                    for (auto &i : x)
                    {
                        for (auto &j : i.second)
                        {
                            for (auto &k : j.second)
                            {
                                n += k.second;
                            }
                        }
                    }
                    
                    return n;
                }
                
                inline Count totalD() const { return totalIDData(dls); }
                inline Count totalI() const { return totalIDData(ins); }
                
                inline Count total() const { return totalS() + totalI() + totalD(); }
            };
            
            struct Metrics
            {
                // Number of reads not overlapping trimmed regions
                Count out = 0;

                // Number of reads overlapping trimmed regions
                Count in = 0;
                
                std::map<AttKey, StratData> sd;
                
                // Regional performance with and without trimming
                Chr2DInters r1, r2, r5;
                
                // Total number of reads
                inline Count total() const { return in + out; }
            };
            
            // Library details from alignments
            Library bamLib;
            
            // Metrics for samples, before and after calibration
            Metrics samp, decoy;
            
            // BAM writer for merging
            std::shared_ptr<BAMWriter> wM;
            
            // This needs to be populated ...
            ParitionCount _pc;
            
            CustomMap<SequinID, Read> rn() const override;
            
            ParitionCount pc() const override { return _pc; }
            
            const Library &lib() const override { return bamLib; }
            
            Count binN(Bin x) const override;

            Count total() const override { return samp.total() + decoy.total(); }
            
            // Counts for a given sequin
            Count count(const SequinID &x) const override;
            
            // Calibration factor
            Proportion calibF(Calibration) const override { return _p; };
                   
            // Number of reads after calibration
            Count calibA(Calibration) const override { return _aSeq; }
                   
            // Dilution after calibration
            Proportion calibD(Calibration) const override { return (Proportion) _aSeq / (_aSam + _aSeq); }
            
            /*
             * Metrics only defined if calibration
             */
            
            Proportion _p;
            Count _tar;
            Count _bSam; // Sample reads before calibration
            Count _aSam; // Sample reads after calibration
            Count _bSeq; // Sequin reads before calibration
            Count _aSeq; // Sequin reads after calibration
        };
        
        typedef std::function<void (ParserBAM::Data &x, const DInter *, const DInter *, bool trimmed)> D1;
        typedef std::function<void ()> D2;

        static Results analyze(const FileName &, const FileName &, const Options &, D1, D2);

        static Results analyze(const FileName &f1, const FileName &f2, const Options &o, D1 d1)
        {
            return DecoyAnalyzer::analyze(f1, f2, o, d1, []() {});
        }
        
        static Results analyze(const FileName &f1, const FileName &f2, const Options &o)
        {
            return DecoyAnalyzer::analyze(f1, f2, o, [](ParserBAM::Data &, const DInter *, const DInter *, bool) {}, []() {});
        }
        
        /*
         * Options to define writing error outputs from alignment results. Single decoy chromosome assumed.
         */
        
        struct ErrorOptions : WriterOptions
        {
            ErrorOptions(const FileName &tsvE,
                         const BedData &r1,
                         const BedData &r2,
                         const ChrID &chr,
                         std::shared_ptr<VCFLadder> v1,
                         bool isDecoy,
                         WriterOptions o) : WriterOptions(o), isDecoy(isDecoy), debug(false), r1(r1), r2(r2), chr(chr), tsvE(tsvE), v1(v1) {}
            
            typedef bool Calibrated;
            
            bool debug;
            
            // Metrics data to write
            std::map<DecoyAnalyzer::Results::Metrics *, Calibrated> data;
            
            ChrID chr;
            BedData r1, r2;
            bool isDecoy;
            FileName tsvE;
            std::shared_ptr<VCFLadder> v1;
        };
        
        static void writeE(const ErrorOptions &);
    };
}

#endif
