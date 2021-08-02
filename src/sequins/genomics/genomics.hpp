#ifndef GENOMICS_HPP
#define GENOMICS_HPP

#include "data/data.hpp"
#include "data/split.hpp"
#include "data/bundle.hpp"
#include "tools/tools.hpp"
#include "sequins/decoy_analyzer.hpp"

namespace Anaquin
{
    const auto GDecoyChrQS = "chrQS";
    const auto GDecoyChrQR = "chrQR";
    const auto GDecoyChrQL = "chrQL";
    
    struct GResource : public Resource
    {
        GResource(const Path &, const FileName &, const FileName &, Build);
    };

    struct LResource : public Resource
    {
        LResource(const Path &, const FileName &, const FileName &, Build);
    };

    inline Resource GInfoCode(const Path &p)
    {
        return LResource(p + "/synthetic", "sequin_barcodes", ".tsv", Build::None);
    }

    inline Resource GSeqDecoy(const Path &p)
    {
        return GResource(p + "/genome/chrQ", "genome_chrQ_decoys", ".fa", Build::None);
    }

    inline Resource CSeqDecoy(const Path &p)
    {
        return GResource(p + "/cancer/chrQ", "genome_chrQ_decoys", ".fa", Build::None);
    }

    inline Resource CSeqFA(const Path &p)
    {
        return GResource(p + "/cancer", "sequin_sequences", ".fa", Build::None);
    }

    inline Resource GSeqFA(const Path &p)
    {
        return GResource(p + "/genome", "sequin_sequences", ".fa", Build::None);
    }

    inline Resource CRegionBED(const Path &p, Build x)
    {
        return GResource(p + "/cancer", "sequin_regions", ".bed", x);
    }

    inline Resource GRegionBED(const Path &p, Build x)
    {
        return GResource(p + "/genome", "sequin_regions", ".bed", x);
    }

    inline Resource CVarVCF(const Path &p, Build x)
    {
        return GResource(p + "/cancer", "sequin_smallvariants", ".vcf", x);
    }

    inline Resource GVarVCF(const Path &p, Build x)
    {
        return GResource(p + "/genome", "sequin_smallvariants", ".vcf", x);
    }
    
    inline Resource GFeatBED(const Path &p, Build x)
    {
        return GResource(p + "/genome", "sequin_features", ".bed", x);
    }

    inline Resource GAttrBED(const Path &p)
    {
        return GResource(p + "/genome", "sequin_attributes", ".tsv", Build::None);
    }

    inline Resource GSynTSV(const Path &p)
    {
        return GResource(p + "/synthetic", "synthetic_ladder", ".tsv", Build::None);
    }

    inline std::string bin2Label(Bin x)
    {
        switch (x)
        {
            case MI: { return "MSI";        }
            case HP: { return "HP";         }
            case IF: { return "Info";       }
            case MT: { return "Mito";       }
            case HL: { return "HLA";        }
            case LD: { return "Ladder";     }
            case IM: { return "Immune";     }
            case ES: { return "Sample";     }
            case VC: { return "Vector";     }
            case SV: { return "Structural"; }
            case GR: { return "Germline";   }
            case SO: { return "Somatic";    }
            case MS: { return "Micro";      }
        }
    }

    Bin GBin(const SequinID &, std::shared_ptr<VCFLadder> v1);
    inline Bin GBin(const SequinID &x) { return GBin(x, nullptr); }
    
    struct GVariantResults
    {
        // Measured coverage for reference and variants
        CustomMap<SequinID, Measured> R, V;
    };

    typedef std::map<SequinID, Proportion> Norms;

    // Mean sampling coverage from tsvR
    Coverage meanSamp(const FileName &, const Label &);

    inline std::string gt2str(Genotype x)
    {
        switch (x)
        {
            case Genotype::MSI:         { return "MSI";          }
            case Genotype::Somatic:     { return "Somatic";      }
            case Genotype::Homozygous:  { return "Homozygous";   }
            case Genotype::Heterzygous: { return "Heterozygous"; }
        }
    }

    inline std::string var2str(Variation x)
    {
        switch (x)
        {
            case Variation::SNP:         { return "SNP";             }
            case Variation::Deletion:    { return "Indel_Deletion";  }
            case Variation::Insertion:   { return "Indel_Insertion"; }
            case Variation::Duplication: { return "Duplication";     }
            case Variation::Inversion:   { return "Inversion";       }
        }
    }

    // Convert genomics sequin name to it's standard
    inline StandardID GSeq2Std(const SequinID &x)
    {
        return noLast(noFirst(x, "_"), "_");
    }
    
    inline bool GValid(Bin) { return true; }

    bool isHP(const Bin x);
    bool isHP(const SequinID &x);
    
    bool isMS(const Bin x);
    bool isMS(const SequinID &x);

    bool isSoma(const Bin x);
    bool isSoma(const SequinID &x);

    bool isVector(const Bin x);
    bool isVector(const SequinID &x);

    // Valid for germline analysis?
    inline bool isGermA(Bin x)
    {
        return x != Bin::SO && x != Bin::SV;
    }

    // Valid for germline analysis?
    inline bool isGermA(const SequinID &x)
    {
        return isGermA(GBin(x));
    }
    
    enum class GDecoyStatus
    {
        Before,
        AfterLC,
        AfterSC
    };
    
    struct GDecoyOptions : public DecoyAnalyzer::Options
    {
        Proportion ladC = NO_CALIBRATION;
        
        /*
         * 1. Decoy analysis
         *
         *     - Perform standard sequin and ladder percentage/absolute calibration
         *     - writeL1 is the output ladder alignment file
         *
         * 2. Ladder calibration
         *
         *     - writeL1 calibrated to writeL2
         *
         * 3. Sequin mirror calibration
         *
         *     - writeMC is the output calibration file
         */
        
        bool debug = false;
        
        // Inout alignment file for writeL2
        FileName inputL2;
        
        // Output alignment file for ladder alignments before and after calibration
        FileName writeL1, writeL2;
        
        FileName inputMC;
        
        // Output alignment file for sequin mirror calibration
        FileName writeMC;
        
        // TSV output files
        FileName tsvR, tsvA, tsvE, tsvF;
        
        // Ladder outputs before and after ladder calibration
        FileName tsvL1, tsvL2;
        
        // tsvL but not merged
        FileName tsvL1NotMerged;
        
        // How sequin mirror calibration is performed
        CalibrateMethod meth = CalibrateMethod::Mean;
        
        // Only for Method::Custom
        double customSequinThreshold = 0;

        static GDecoyOptions create(const FileName &, const AnalyzerOptions &);
        static GDecoyOptions create(const FileName &writeS,
                                    const FileName &writeD,
                                    const FileName &writeM,
                                    const FileName &writeT,
                                    const FileName &writeL1,
                                    const FileName &writeL2,
                                    const FileName &calibF,
                                    const FileName &tsvE,
                                    const FileName &tsvR,
                                    const FileName &tsvF,
                                    const FileName &tsvA,
                                    const FileName &tsvL1,
                                    const FileName &tsvL2,
                                    const FileName &index,
                                    Proportion seqC,
                                    Proportion ladC,
                                    const AnalyzerOptions &o);
    };

    struct GDecoyResults
    {
        // Directly from DecoyAnalyzer (before calibration)
        DecoyAnalyzer::Results B1;
        
        // Directly from DecoyAnalyzer (after ladder calibration)
        DecoyAnalyzer::Results B2;

        // Directly from DecoyAnalyzer (after sequin calibration)
        DecoyAnalyzer::Results B3;
        
        Library lib;
        DecoyAnalyzer::Results::Metrics samp, before, after;
        
        struct VariantData
        {
            Count R, V;
        };
        
        // Allele frequency for each sequin (only sequins with a variant)
        CustomMap<SequinID, VariantData> bAF, aAF;
        
        struct MultiCalibrate
        {
            // Calibration factors (only valid if calib == -1)
            Norms norms;
        };
        
        struct SingleCalibrate
        {
            Proportion p;
            
            // Number of reads after calibration
            Count after;
        };
        
        VariantData v;
        MultiCalibrate  c1;
        SingleCalibrate c2;
        
        // Construct allele frequeny ladder
        static std::map<SequinID, VariantData> buildAF(const DecoyAnalyzer::Results::Metrics &, bool);
        
        // Write TSV for variants
        static void writeV(const GDecoyResults &, const GDecoyOptions &);
        
        // Write TSV for synthetic ladders
        static void writeL(const FileName &, const DecoyAnalyzer::Results &, const DecoyAnalyzer::Options &);
        
        // Write TSV for each region (for calibration)
        static void writeR(const GDecoyResults &, const GDecoyOptions &);
        
        // Write TSV for each attribute (or feature) region
        static void writeF(const GDecoyResults &, const GDecoyOptions &);
    };

    typedef std::function<void (GDecoyStatus, const ParserBAM::Data &x, const DInter *, const DInter *, bool)> G1;
    typedef std::function<void (GDecoyStatus)> G2; // Invoked when completed

    GDecoyResults GDecoyAnalysis(const FileName &, const FileName &, const GDecoyOptions &, G1, G2);
}

#endif
