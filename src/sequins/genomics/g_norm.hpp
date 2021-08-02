#ifndef G_NORM_HPP
#define G_NORM_HPP

#include "data/data.hpp"
#include "sequins/genomics/g_split.hpp"

namespace Anaquin
{
    struct GNorm
    {
        struct Stats
        {
            // Split statistics for all samples
            std::vector<GSplit::Stats> S;
            
            // Scaling factor for all samples
            std::vector<Proportion> P;
            
            // Linear model before normalization
            std::vector<SequinStats> l1;
            
            // Linear model after normalization
            std::vector<SequinStats> l2;
            
            inline Index shallowest() const
            {
                return std::distance(P.begin(), std::max_element(P.begin(), P.end()));
            }
            
            inline double scale(Index i) const
            {
                return P[i] / P[shallowest()];
            }
        };
        
        enum class Method
        {
            Simple
        };
        
        struct Options : public GSplit::Options
        {
            Options() : meth(Method::Simple) {}
            
            Method meth;
        };
        
        /*
         * Similarity to Andre's recommendation, apply normalisation on ladder sequins and write
         * new Count dst.
         */
        
        static void writeNorm(const FileName &, const FileName &, GNorm::Stats &, const GNorm::Options &);
        
        static void writeQuin(const FileName &, const Stats &, const SOptions &);
        
        static Stats analyze(const std::vector<FileName> &, const std::vector<FileName> &, const Options  &);
        static void  report (const std::vector<FileName> &, const std::vector<FileName> &, const Options  &o = Options());
    };
}

#endif
