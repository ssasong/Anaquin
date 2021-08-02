#ifndef SPLIT_HPP
#define SPLIT_HPP

#include <fstream>
#include "Kallisto.hpp"
#include "stats/linear.hpp"
#include "data/library.hpp"
#include "data/analyzer.hpp"
#include "stats/ss/stats.hpp"
#include "data/resources.hpp"
#include "sequins/sequins.hpp"
#include "tools/calibrator.hpp"
#include "parsers/parser_bam.hpp"

namespace Anaquin
{
    struct SOptions : public KOptions, public AnalyzerOptions
    {
        SOptions() : skipMerge(true), seqC(NO_CALIBRATION), ladC(NO_CALIBRATION) {}
        
        bool skipMerge;
        
        // How much to calibrate for sequins and ladder in stage one?
        Proportion seqC;

        // Ladder calibration
        Proportion ladC;

        // Are we doing sequin calibration?
        inline bool isSCalib() const
        {
            return seqC != NO_CALIBRATION;
        }

        // Are we doing ladder calibration?
        inline bool isLCalib() const
        {
            return ladC != NO_CALIBRATION;
        }
    };
    
    // Statistics for calibration
    struct SCStats
    {
        Proportion p = NAN;
       
        // Target number of reads
        Count tar;

        Count bSam = 0; // Sample reads before calibration
        Count aSam = 0; // Sample reads after calibration
        Count bSeq = 0; // Sequin reads before calibration
        Count aSeq = 0; // Sequin reads after calibration
       
        // Output file names from calibration
        FileName o1, o2;
    };

    struct FQStats : public CommonResults
    {
        // Calibration statistics for sequins
        SCStats C;
        
        // Calibration statistics for ladders
        SCStats L;
        
        FileName tsv;
        
        // Kallisto statistics
        KStats K;

        typedef std::map<SequinID, std::vector<double>> KMCount;
        
        struct Abundance
        {
            // Shared k-mer Count
            KMCount s2s;
            
            // Unique k-mer Count
            KMCount s2u;
            
            struct Descriptive
            {
                // Descriptive statistics (unique k-mers)
                std::map<SequinID, Count> mins, maxs;
                
                // Descriptive statistics (unique k-mers)
                std::map<SequinID, Coverage> mus, q25, q75, meds, sds;
            };
            
            // Descriptive statistics for shared and unique k-mers
            Descriptive d2s, d2u;
        };
        
        Abundance R;
        
        /*
         * Implementation for CommonResults
         */

        ParitionCount pc() const override
        {
            ParitionCount x;
            
            forBin([&](Bin b)
            {
                x[b] = this->binN(b);
            });

            return x;
        }

        CustomMap<SequinID, Read> rn() const override
        {
            CustomMap<SequinID, Read> m;
            
            for (const auto &x : K.seqs)
            {
                auto p = K.sqc.find(noPID(x));
                m[x] = p ? *p : 0;
            }
            
            return m;
        }
        
        const Library &lib() const override { return K.lib; }

        // Multiply by two to reach number of reads (not paired)
        Count binN(Bin x) const override { return 2 * K.binN(x); }

        Proportion binP(Bin x) const override { return K.binP(x); }
        
        Count count(const SequinID &x) const override
        {
            return K.sqc.count(x) ? K.sqc.at(x) : 0;
        }
        
        // Scaling factor
        Proportion calibF(Calibration) const override { return C.p; };
                 
        // Number of reads after calibration
        Count calibA(Calibration) const override { return C.aSeq; }
                 
        // Dilution after calibration
        Proportion calibD(Calibration) const override { return (Proportion) C.aSeq / (C.aSam + C.aSeq); }
    };

    // Version number by translation
    std::string SVersion(const Reference &, const KStats &);
    
    struct SAlignStats
    {
        SAlignStats()
        {
            ins[Bin::ES]; // Make sure we always have an entry
            ins[Bin::GR]; // Make sure we always have an entry
        }

        // Insertion size in mean and SD
        std::string SInsert(Bin b = Bin::GR, Base min = 10, Base max = 1000) const;
        
        // Frequency table of insertion size for each bin
        std::map<Bin, std::map<Base, Count>> ins;
    };
    
    std::string SLibraryRun();
    std::string SLibraryInst();
    std::string SLibraryFlow();
    std::string SLibraryLane();

    void writeSTable_1(const FileName &, const FileName &, const SOptions &, Count nExp, Count nObs, Count nCV, const Label &, const Label &, float d = 0);
    
    void writeKmers(const FileName &, const FQStats &, const SOptions &);

    void SMerge(FQStats &stats, const SOptions &o, const std::set<Bin> &only);

    void SWriteLadder(std::shared_ptr<Ladder>, const FileName &, const FQStats &, const SOptions &);

    template <typename Stats> void SWriteLadderPostCalib(std::shared_ptr<Ladder> l3,
                                                         const Stats &stats,
                                                         const SOptions &o,
                                                         bool isCalib)
    {
        const auto l2 = (!isCalib ? "_ladder.tsv" : "_ladder_calibrated.tsv");
        
        // Generating "_ladder.tsv" or "_ladder_calibrated.tsv"
        SWriteLadder(l3, o.name + l2, stats, o);
    }

    // Writing reads to TSV file
    void SWriteReads(Product, const FileName &, const FQStats &, const SOptions &);

    void SKallisto(FQStats &, const FileName &, const FileName &, const SOptions &);

    /*
     * -------------------- Calibration by percentage --------------------
     */

    struct SCalibratePFiles
    {
        virtual FileName i1(const SOptions &) const = 0;
        virtual FileName i2(const SOptions &) const = 0;
        virtual FileName o1(const SOptions &) const = 0;
        virtual FileName o2(const SOptions &) const = 0;
    };
    
    /*
     * Default implementation for SCalibratePFiles
     */
    
    struct SCalibrateDefault : public SCalibratePFiles
    {
        SCalibrateDefault(const Label &x, const Path &f1 = "", const Path &f2 = "") : x(x), f1(f1), f2(f2) {}
        
        FileName i1(const SOptions &o) const override
        {
            const auto f1_ = f1.empty() ? o.work : f1;
            return o.writeBAM() ? f1_ + "/" + o.name + "_" + x + ".bam" :
                                  f1_ + "/" + o.name + "_" + x + "_1.fq.gz";
        }
        
        FileName i2(const SOptions &o) const override
        {
            const auto f1_ = f1.empty() ? o.work : f1;
            return o.writeBAM() ? "" : f1_ + "/" + o.name + "_" + x + "_2.fq.gz";
        }
        
        FileName o1(const SOptions &o) const override
        {
            const auto f2_ = f2.empty() ? o.work : f2;
            return o.writeBAM() ? f2_ + "/" + o.name + "_" + x + "_calibrated.bam" :
                                  f2_ + "/" + o.name + "_" + x + "_calibrated_1.fq.gz";
        }
        
        FileName o2(const SOptions &o) const override
        {
            const auto f2_ = f2.empty() ? o.work : f2;
            return o.writeBAM() ? "" : f2_ + "/" + o.name + "_" + x + "_calibrated_2.fq.gz";
        }
        
        const Label x;
        const FileName f1, f2;
    };

    SCStats SCalibrateP(const std::vector<Bin> &, Proportion, const FQStats &, const SOptions &, const SCalibratePFiles &);

    inline SCStats SCalibrateP(Bin b, Proportion p, const FQStats &stats, const SOptions &o, const SCalibratePFiles &f)
    {
        return SCalibrateP(std::vector<Bin> { b }, p, stats, o, f);
    }
}

#endif
