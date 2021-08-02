#ifndef SEQUINS_HPP
#define SEQUINS_HPP

#include "stats/stats.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    inline void mergeL(const FileName &src, const FileName &dst)
    {
        const auto tmp1 = tmpFile();
        const auto tmp2 = tmpFile();

        RFilterC(src, tmp1, std::set<Label> { "NAME", "STANDARD" }, false);
        RAggregateSum(tmp1, tmp2, "SEQUIN");
        RFilterC(tmp2, tmp1, std::set<Label> { "SEQUIN", "STOICHOMETRY", "COPY", "MEAN", "READ" }, true);
        RRenameC(tmp1, tmp2, std::vector<Column> { "NAME", "STOICHOMETRY", "COPY", "MEAN", "READ" });
    
        RApply(tmp2, tmp1, "STOICHOMETRY", [&](const std::string &x)
        {
            return "1";
        });

        RApply(tmp1, dst, "COPY", [&](const std::string &x)
        {
            if (x == "4")  { return std::string("2"); }
            if (x == "16") { return std::string("4"); }
            if (x == "64") { return std::string("8"); }
            return x;
        });
    }

    template <typename O> void writeLTable(const FileName &src, const FileName &dst, const O &o)
    {
        const auto tmp = tmpFile();
        RLadTable(src, tmp, "NAME", o.bam ? "MEAN" : "Q50");
        o.generate(dst);
        o.writer->open(dst);
        o.writer->write(readFile(tmp));
        o.writer->close();
    }

    inline SequinID toSeqID(const Tok &x)
    {
        auto toks = std::vector<std::string>();
        split(x, "_", toks);
        assert(toks.size() >= 2 && toks.size() <= 4);
        
        // Eg: S0552_CM_001_A
        if (toks.size() == 4)
        {
            return noFirst(x, "_"); // Eg: CM_001_A
        }
        
        return x;
    }
    
    inline StandardID toStand(const Tok &x)
    {
        const auto seq = toSeqID(x);
        auto toks = std::vector<std::string>();
        split(seq, "_", toks);
        return toks.size() == 3 ? noLast(seq, "_") : seq;
    }
    
    inline Label noPID(const Label &x)
    {
        return isBegin(x, "S") && countT(x, "_") > 2 ? noFirst(x, "_") : x;
    }
}

#endif
