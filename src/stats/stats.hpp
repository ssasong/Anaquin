#ifndef STATS_HPP
#define STATS_HPP

#include <set>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "data/data.hpp"
#include "tools/tools.hpp"
#include "stats/linear.hpp"
#include "stats/ss/matrix.hpp"

namespace Anaquin
{
    /*
     * Implementation of edgeR's TMM method. Adopted from Andre's method.
     */

    SS::Matrix TMM(const SS::Matrix &x);
    
    template <typename T1, typename T2> typename T1::value_type quant(const T1 &x, T2 q)
    {
        assert(q >= 0.0 && q <= 1.0);
        
        auto y = std::vector<typename T1::value_type>(x.begin(), x.end());
        auto const i = q * x.size();
        
        std::nth_element(y.begin(), y.begin() + i, y.end());
        return (y.size() % 2 || y.size() == 1) ? y.at(i) : 0.5 * (y.at(i-1) + y.at(i));
    }

    template <typename T> typename T::value_type med(const T &x)
    {
        return quant(x, 0.5);
    }

    template <typename T> typename T::value_type min(const T &x)
    {
        return *(std::min_element(x.begin(), x.end()));
    }
    
    template <typename T> typename T::value_type max(const T &x)
    {
        return *(std::max_element(x.begin(), x.end()));
    }
    
    Count RFilterR(const FileName &, const FileName &, const std::set<Column> &);
    
    void RReplaceC(const FileName &, const FileName &, const Column &, const std::map<std::string, std::string> &);
    void RRenameC(const FileName &, const FileName &, const std::vector<Column> &);

    void RFilterC(const FileName &, const FileName &, const std::set<Column> &, bool keep = false);
    
    void RFilterC(const FileName &, const FileName &, const std::set<std::size_t> &, bool keep = false);
    
    typedef std::function<double (const std::vector<double> &)> Apply;

    std::map<double, double> RBinaryTSV(const FileName &, const Column &, const Column &);

    // Read TSV and return it as a dictionary
    std::map<std::string, std::string> RReadTSV(const FileName &, const Column &, const Column &);

    // Read TSV and return it as a vector
    std::vector<std::string> RReadTSV(const FileName &, const Column &);
    
    enum class Imputation
    {
        None,
        ToZero,
        Remove
    };
    
    void RAggregate(const FileName &, const FileName &, const Column &, Apply, float d = 0, Imputation impute = Imputation::None);
    void RAggregateSD(const FileName &, const FileName &, const Column &, float d = 0, Imputation impute = Imputation::None);
    void RAggregateCount(const FileName &, const FileName &, const Column &, float d = 0);
    void RAggregateMean(const FileName &, const FileName &, const Column &, float d = 0, Imputation impute = Imputation::None);
    void RAggregateSum(const FileName &, const FileName &, const Column &, float d = 0, Imputation impute = Imputation::None);
    void RAggregateMax(const FileName &, const FileName &, const Column &, float  = 0);
    void RAggregateMin(const FileName &, const FileName &, const Column &, float d = 0);
    void RAggregateQ0(const FileName &, const FileName &, const Column &, float n = 0);
    void RAggregateQ25(const FileName &, const FileName &, const Column &, float n = 0);
    void RAggregateQ50(const FileName &, const FileName &, const Column &, float n = 0);
    void RAggregateQ75(const FileName &, const FileName &, const Column &, float n = 0);
    void RAggregateQ100(const FileName &, const FileName &, const Column &, float n = 0);

    typedef std::function<std::string (const std::string &)> RApplyF;
    void RApply(const FileName &, const FileName &, const Column &, RApplyF);

    enum class Arithmetic
    {
        Add,
        Subtract,
        Multipy,
        Divide,
    };
    
    // Apply arithmetic operations to two selected columns
    std::vector<double> RArith(const FileName &, const Column &, const Column &, Arithmetic, bool skipNA = true);
    
    inline std::vector<double> RSubtract(const FileName &src, const Column &x, const Column &y)
    {
        return RArith(src, x, y, Arithmetic::Subtract);
    }

    inline std::vector<double> RAdd(const FileName &src, const Column &x, const Column &y)
    {
        return RArith(src, x, y, Arithmetic::Add);
    }

    void RGrep(const FileName &, const FileName &, const Column &, const std::set<std::string> &, bool keep = true, bool isNum = false, bool exact = false);

    void RNoLast(const FileName &, const FileName &, const Column &, const Tok &);
    
    inline void RGrep(const FileName &src,
                      const FileName &dst,
                      const Column   &lab,
                      const std::string &val,
                      bool keep = true,
                      bool isNum = false,
                      bool exact = false)
    {
        RGrep(src, dst, lab, std::set<std::string> { val }, keep, isNum, exact);
    }

    double RSum(const FileName &, const Column &);
    double RMean(const FileName &, const Column &);
    Count  RCount(const FileName &, const Column &, const std::string &, bool exact = true);

    inline Count RCount(const FileName &file, const Column &c)
    {
        return RCount(file, c, "", false);
    }
    
    std::vector<std::string> RList(const FileName &, const Column &);
    
    // Return number of columns
    Count RCountCols(const FileName &);
    
    bool RHead(const FileName &, const Column &);
    
    void RMeanCV(const FileName &, const FileName &, const Column &, const Column &, Count nExp = 2, Count Obs = 2, Count nCV = 2, float d = 0);

    SequinStats RLinear(const FileName &src, const Column &, const Column &, const Column &);

    // Returns mean ratios for a ladder table
    double RLadTable(const FileName &, const FileName &, const Column &, const Column &);

    struct LadderResult
    {
        double meanR;
        
        struct CopyNumber
        {
            Coverage mean, sd;
            Coverage q0, q25, q50, q75, q100;
            
            // Coefficent of variation
            double cv;
            
            // Mean ratios
            double ratio;
            
            // 1n (deletion) coverage; CV
            inline std::string broad() const
            {
                return S2(mean) + " ; " + S2(cv);
            }
        };
        
        std::map<unsigned, CopyNumber> cn;
    };
    
    /*
     * The source file must have a column for the median ("Q50"), and also "NAME".
     */
    
    LadderResult RLadder(const FileName &);
}

#endif
