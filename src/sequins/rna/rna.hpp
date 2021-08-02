#ifndef RNA_HPP
#define RNA_HPP

#include <map>
#include <fstream>
#include <iostream>
#include "data/bundle.hpp"
#include "parsers/parser_csv.hpp"

namespace Anaquin
{
    const auto GDecoyChrIS = "chrIS";
    
    struct RResource : public Resource
    {
        RResource(const Path &path) { Resource::path = path; }
        RResource(const Path &path, const FileName &file, const FileName &ext)
        {
            Resource::path = Bundle::latest(path + "/" + file, ext);
        }
    };

    inline StandardID RSeq2Std(const SequinID &x) { return noLast(x, "_"); }

    inline Resource RNAGBed(const Path &p)
    {
        const auto x = RResource(p + "/transcriptome", "rnasequin_regions", ".bed");
        std::map<StandardID, Base> m1, m2;
        
        const auto tmp = tmpFile();
        
        ParserCSV::parse(Reader(x.path), [&](const ParserCSV::Data &x, Progress i)
        {
            const auto std = RSeq2Std(x[3]);
            const auto x1  = std::stoi(x[1]);
            const auto x2  = std::stoi(x[2]);
            
            if (!m1.count(std) || x1 < m1[std]) { m1[std] = x1; }
            if (!m2.count(std) || x2 > m2[std]) { m2[std] = x2; }
        }, "\t");
        
        std::ofstream w;
        w.open(tmp);
        
        for (const auto &i : m1)
        {
            w << GDecoyChrIS << "\t" << m1[i.first] << "\t" << m2[i.first] << "\t" << i.first << std::endl;
        }
        
        w.close();
        return RResource(tmp);
    }

    inline Resource RNAMix(const Path &p)
    {
        return RResource(p + "/transcriptome", "rnasequin_isoforms_", ".tsv");
    }
    
    inline Resource RNATBed(const Path &p)
    {
        return RResource(p + "/transcriptome", "rnasequin_regions", ".bed");
    }
    
    inline Resource RNAFA(const Path &p)
    {
        return RResource(p + "/transcriptome", "rnasequin_sequences_", ".fa");
    }

    inline Resource RNADecoy(const Path &p)
    {
        return RResource(p + "/transcriptome", "rnasequin_decoychr_", ".fa");
    }

    inline Bin RBin(const SequinID &) { return GR; }
    
    inline bool RValid(Bin x) { return x == GR || x == ES; }
}

#endif
