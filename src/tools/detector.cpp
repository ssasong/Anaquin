#include <fstream>
#include "tools/tools.hpp"
#include "tools/detector.hpp"
#include "parsers/parser_bam.hpp"
#include "parsers/parser_vcf.hpp"
#include "writers/file_writer.hpp"
#include "writers/terminal_writer.hpp"

using namespace Anaquin;

static bool fromVCFRef(const FileName &file, Build &g)
{
    std::fstream f;
    f.open(file, std::fstream::in);
    
    std::string line;
    while (std::getline(f, line))
    {
        if (!line.empty() && line[0] == '#')
        {
            // Eg: "reference=file:///share/NatureP/Homo_sapiens_assembly38.fasta"
            if (isSubstr(line, "reference"))
            {
                if (isSubstr(line, "Homo_sapiens_assembly38.fasta"))
                {
                    g = Build::hg38;
                    return true;
                }
            }
            
            continue;
        }
        
        break;
    }
    
    f.close();
    return false;
}

Build Detector::fromBAM(const FileName &)
{
    return Build::hg38;
}

Build Detector::fromVCF(const FileName &file)
{
    Build hg;
    
    if (fromVCFRef(file, hg))
    {
        return hg;
    }
    else
    {
        return Build::hg38;
    }
}
