#ifndef G_GERM_HPP
#define G_GERM_HPP

#include "sequins/genomics/g_variant.hpp"

namespace Anaquin
{
    struct GGerm : public GVariant
    {
        bool isValid(const SequinID &x) const override { return isGermA(x); }

        typedef GVariant::Stats Stats;
        
        static Stats analyze(const FileName &, const FileName &, const Options &);
        static void report (const FileName &, const FileName &, const Options &);
    };
}

#endif
