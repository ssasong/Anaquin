#ifndef G_BROAD_VCF_HPP
#define G_BROAD_VCF_HPP

#include "sequins/genomics/g_germ.hpp"

namespace Anaquin
{
    struct GBroadVCF : public GGerm
    {
        typedef GGerm::Options Options;
        
        bool isValid(const SequinID &x) const override { return true; }

        static void report(const FileName &, const Options &);
    };
}

#endif
