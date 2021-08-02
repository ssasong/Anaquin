#ifndef G_TMM_HPP
#define G_TMM_HPP

#include "sequins/genomics/g_split.hpp"

namespace Anaquin
{
    struct GTMM
    {
        struct Stats
        {
            std::vector<Sequence> seqs;            
            SS::Matrix dt, dt_n;
            
            // Scaling factor for the two samples
            Proportion p1, p2;
            
            inline bool isFirstDeeper() const
            {
                return p1 < p2;
            }
            
            // Scalling probability for the deeper sample
            inline double scale() const
            {
                return isFirstDeeper() ? (p1 / p2) : (p2 / p1);
            }
        };

        typedef GSplit::Options Options;
        
//        static void writeNorm(const FileName &, const FileName &, GNorm::Stats &, const GNorm::Options &);
        
        static Stats analyze(const FileName &, const Options &);
        static void  report (const FileName &, const Options &o = Options());
    };
}

#endif
