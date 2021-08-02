#ifndef META_HPP
#define META_HPP

#include "data/bundle.hpp"
#include "sequins/decoy_analyzer.hpp"

namespace Anaquin
{
    const auto MDecoyChrQB = "chrQB";
    const auto MDecoyChrQL = "chrQL";

    struct MResource : public Resource
    {
        MResource(const Path &path, const FileName &file, const FileName &ext)
        {
            Resource::path = Bundle::latest(path + "/" + file, ext);
        }
    };

    inline Resource MetaDecoy(const Path &p)
    {
        return MResource(p + "/metagenome/chrQ", "metagenome_chrQ_decoys_", ".fa");
    }

    inline Resource MetaFA(const Path &p)
    {
        return MResource(p + "/metagenome", "metasequin_sequences_", ".fa");
    }

    inline Resource MetaMix(const Path &p)
    {
        return MResource(p + "/metagenome", "metasequin_abundance_", ".tsv");
    }

    inline Resource MetaBED(const Path &p)
    {
        return MResource(p + "/metagenome", "metasequin_regions_", ".bed");
    }

    inline Resource MetaDBED(const Path &p)
    {
        return MResource(p + "/metagenome/chrQ", "metagenome_regions_chrQ_", ".bed");
    }

    inline Bin MBin(const SequinID &x)
    {
        if (second(x, "_") == "VC")
        {
            return VC;
        }
        else if (first(x, "_") == "LD" || second(x, "_") == "LD")
        {
            return LD;
        }
        else if ((first(x, "_") == "IF" || second(x, "_") == "IF") ||
                 (first(x, "_") == "UQ" || second(x, "_") == "UQ"))
        {
            return IF;
        }
        
        return GR;
    }
        
    inline bool MValid(Bin x) { return x == GR || x == ES || x == LD || x == IF || x == VC; }

    inline StandardID MSeq2Std(const SequinID &x) { return isSubstr(x, "LD_") ? noLast(noFirst(x, "_"), "_") : x; }

    struct MDecoyOptions : public DecoyAnalyzer::Options
    {
        Proportion ladC = NO_CALIBRATION;
        
        // Output files
        FileName tsvR, tsvE, tsvF;
        
        // Ladder outputs before and after ladder calibration
        FileName tsvL1, tsvL2;

        FileName writeL;
        
        Path originalW; // TOOD...
    };

    struct MDecoyResults
    {
        struct ErrorReport
        {
            std::map<unsigned, std::string> text;
        };

        // Directly from DecoyAnalyzer (before calibration)
        DecoyAnalyzer::Results B1;
        
        // Directly from DecoyAnalyzer (after sequin calibration)
        DecoyAnalyzer::Results B2;
        
        // Directly from DecoyAnalyzer (after ladder calibration)
        DecoyAnalyzer::Results B3;
        
        Library lib;
        DecoyAnalyzer::Results::Metrics samp, decoy;
        
        // Write TSV for synthetic ladders
        static void writeL(const DecoyAnalyzer::Results &, const MDecoyOptions &);
        
        // Write TSV for each region (for calibration)
        static void writeR(const MDecoyResults &, const MDecoyOptions &);
        
        // Write TSV for each attribute (or feature) region
        static void writeF(const MDecoyResults &, const MDecoyOptions &);

        static ErrorReport reportE(const FileName &, const FileName &);
    };
    
    MDecoyResults MDecoyAnalysis(const FileName &, const MDecoyOptions &);
}

#endif
