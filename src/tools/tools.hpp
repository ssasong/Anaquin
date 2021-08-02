#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <set>
#include <chrono>
#include <cstdio>
#include <vector>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "data/data.hpp"
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace Anaquin
{
    typedef std::string Tok;
    typedef std::vector<Tok> Toks;

    void createD(const Path &);
    void removeD(const Path &);
    
    // Empty file?
    bool isEmpty(const FileName &);
    
    void runCmd(const std::string &);
    
    void clearAllTmp();
    Path tmpPath();

    void copy(const FileName &, const FileName &);
    
    template <typename O> O cloneO(const O &o)
    {
        auto o_ = o;
        o_.work = tmpPath();
        return o_;
    }
    
    class Compressor
    {
        public:
        
            ~Compressor() { flush(); }
             Compressor(std::shared_ptr<std::ofstream> w, Count n = 10000) : _n(n), _i(0), _w(w) {}
        
            inline void flush() { write(); }

            void write(const std::string &);

        private:

            void write();
        
            // Number of reads to buffer
            const Count _n;
        
            Count _i;
            std::stringstream _buf;
            std::shared_ptr<std::ofstream> _w;
    };
    
    // Compress a string to GZ format
    std::string compressGZ(const std::string &);
    
    // Can we read the head of a FQ GZ file?
    bool headFQGZ(const FileName &);
    
    enum class CopyGZStatus
    {
        Failed,
        Success,
        Corrected,
    };
    
    CopyGZStatus copyGZ(const FileName &, std::shared_ptr<std::ofstream>);
    
    template <typename F, typename O> void logRun(F f,
                                                  const O &o,
                                                  const std::string &s1,
                                                  const std::string &s2)
    {
        o.info(s1);
        using namespace std::chrono;
        auto t1 = high_resolution_clock::now();
        f();
        auto t2 = high_resolution_clock::now();
        o.info(s2);
        o.info((boost::format("%1% seconds") % duration_cast<seconds>(t2 - t1).count()).str());
    }

    inline std::string trim(const std::string &str)
    {
        std::stringstream ss;
        
        auto seenLastClose = false;
        
        for (auto iter = str.rbegin(); iter != str.rend(); iter++)
        {
            if (*iter == ')')
            {
                seenLastClose = true;
            }
            
            if (!seenLastClose)
            {
                continue;
            }
            
            ss << *iter;
        }
        
        auto tmp = ss.str();
        std::reverse(tmp.begin(), tmp.end());
        return tmp;
    }

    inline void reverse(char *s)
    {
        int l = strlen(s);
        int c, i, j;
        
        for (i = 0, j = l - 1; i < j; i++, j--)
        {
            c = s[i];
            s[i] = s[j];
            s[j] = c;
        }
    }

    inline std::string reverse(const std::string &x)
    {
        auto tmp = x;
        std::reverse(tmp.begin(), tmp.end());
        return tmp;
    }
    
    bool exists(const FileName &);

    void mv(const FileName &, const FileName &);
    void rm(const FileName &);
    void removeD(const Path &);

    std::string date();
    
    inline std::string random(std::size_t l, const std::string &x)
    {        
        auto randchar = [&]() -> char
        {
            const auto max_index = (sizeof(x) - 1);
            return x.at(rand() % max_index);
        };

        std::string str(l,0);
        std::generate_n(str.begin(), l, randchar);

        return str;
    }

    template <typename T> static void split(const Tok &x, const Tok &d, T &r)
    {
        r.clear();
        boost::split(r, x, boost::is_any_of(d));
    }
    
    inline std::string extend(const std::string &x, unsigned size)
    {
        auto tmp = x;
        while (tmp.size() < size) { tmp += " "; }    
        return tmp;
    }
    
    inline Tok join(const Toks &x, const std::string &d)
    {
        return boost::algorithm::join(x, d);
    }
    
    inline bool isEnd(const std::string &x, const std::string &y)
    {
        return boost::algorithm::ends_with(x, y);
    }

    inline bool isBegin(const std::string &x, const std::string &y)
    {
        return boost::algorithm::starts_with(x, y);
    }
    
    // Decoy chromosome?
    inline bool isDecoy(const ChrID &x) { return isBegin(x, "chrQ") || isBegin(x, "chrIS"); }

    FileName hackUZero(const FileName &);
    
    template <typename T1, typename T2> T2 getBestOverlap(const T1 &l, std::vector<T2> &v)
    {
        return *std::max_element(v.begin(), v.end(), [&](const T2 x, const T2 y)
        {
            return x->l().overlap(l) < y->l().overlap(l);
        });
    }

    template <typename T> void complement(T &x)
    {
        std::map<char, char> m = { { 'A', 'T' },
                                   { 'T', 'A' },
                                   { 'G', 'C' },
                                   { 'C', 'G' },
                                   { 'N', 'N' } };
        
        std::transform(x.begin(), x.end(), x.begin(), [&](char c)
        {
            return m.at(c);
        });
    }
    
    template <typename T> T revcomp(const T &x)
    {
        T r(x);
        std::transform(x.rbegin(), x.rend(), r.begin(), [](char c)
        {
            switch(c)
            {
                case 'A': return 'T';
                case 'C': return 'G';
                case 'G': return 'C';
                case 'T': return 'A';
                default: return 'N';
            }
            return 'N';
        });

        return r;
    }
    
    template <typename T> T CHROM(const T &x)
    {
        return !x.empty() && (isdigit(x.front()) || x == "x" || x == "X" || x == "y" || x == "Y" || x == "m" || x == "M") ? "chr" + x : x;
    }

    inline std::string toString(Build x)
    {
        switch (x)
        {
            case Build::None: { return "";     }
            case Build::chrQ: { return "chrQ"; }
            case Build::hseq: { return "hseq"; }
            case Build::gr37: { return "gr37"; }
            case Build::gr38: { return "gr38"; }
            case Build::hg19: { return "hg19"; }
            case Build::hg38: { return "hg38"; }
        }
    }
    
    FileName tmpFile();
    FileName script2File(const Scripts &);

    std::vector<FileName> listFiles(const Path &, const std::string &pat = "", const std::string &ext = "");
    
    // Run an external program and return the output
    std::string run(const Command &);

    void runScript(const Scripts &, const std::string &);

    inline Tok tokN(const Label &x, const Tok &d, unsigned i)
    {
        Toks toks;
        toks.clear();
        split(x, d, toks);
        return toks[i];
    }

    inline Tok first(const Label &x, const Tok &d)
    {
        Toks toks;
        toks.clear();
        split(x, d, toks);
        return toks.front();
    }

    inline Tok second(const Label &x, const Tok &d)
    {
        Toks toks;
        toks.clear();
        split(x, d, toks);
        return toks[1];
    }

    inline Count countT(const Label &x, const Tok &d)
    {
        Toks toks;
        toks.clear();
        split(x, d, toks);
        return toks.size();
    }
    
    inline Tok last(const Label &x, const Tok &d)
    {
        Toks toks;
        toks.clear();
        split(x, d, toks);
        return toks.back();
    }
    
    inline Tok noLast(const Label &x, const Tok &d)
    {
        Toks toks;
        split(x, d, toks);
        toks.pop_back();
        return join(toks, d);
    }

    inline Tok noFirst(const Label &x, const Tok &d)
    {
        Toks toks;
        split(x, d, toks);
        toks.erase(toks.begin());
        return join(toks, d);
    }

    inline Tok replace(const Label &x, const Tok &s1, const Tok &s2)
    {
        auto t = x;
        boost::replace_all(t, s1, s2);
        return t;
    }
    
    inline std::string remove(const Label &s1, const Label &s2)
    {
        return replace(s1, s2, "");
    }

    int parseChrID(const ChrID &);

    template <typename T1, typename T2, typename T3> std::map<T1, std::map<T2, T3>>
            add(const std::map<T1, std::map<T2, T3>> &x1,
                const std::map<T1, std::map<T2, T3>> &x2)
    {
        std::map<T1, std::map<T2, T3>> x = x1;
        
        for (const auto &i : x2)
        {
            if (!x.count(i.first))
            {
                x[i.first] = i.second;
            }
            else
            {
                for (const auto &j : i.second)
                {
                    x[i.first][j.first] += j.second;
                }
            }
        }
        
        return x;
    }
    
    template <typename T> T add(const T &x1, const T &x2)
    {
        T x = x1;
        
        for (const auto &i : x2)
        {
            x[i.first] += i.second;
        }
        
        return x;
    }

    template <typename T1, typename T2> const T1 &max(const std::map<T1, T2> &x)
    {
        return std::max_element(std::begin(x), std::end(x),
                                [&] (const typename std::map<T1, T2>::value_type &t1,
                                     const typename std::map<T1, T2>::value_type &t2)
        {
            return t1.second < t2.second;
        })->first;
    }
    
    inline bool isNumber(const std::string &s)
    {
        return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) { return !std::isdigit(c) && c != '.'; }) == s.end();
    };

    inline bool isSubstr(const std::string &x, const std::string &y)
    {
        return x.find(y) != std::string::npos;
    }

    template <typename K, typename V> const V* findBySub(const std::map<K, V> &x, const K &key)
    {
        for (auto &i : x)
        {
            if (isSubstr(i.first, key))
            {
                return &(i.second);
            }
        }

        return nullptr;
    }
    
    template <typename K, typename V> struct CustomMap : public std::map<K, V>
    {
        CustomMap<K, V> &operator=(const std::map<K, V> &x)
        {
            this->insert(x.begin(), x.end());
            return *this;
        }

        inline V* find(const K &key)
        {
            return findBySub(*this, key);
        }

        inline const V* find(const K &key) const
        {
            return findBySub(*this, key);
        }
    };

    std::string readFile(const FileName &);

    inline ReadName trimRName(const ReadName &x)
    {
        auto x_ = isSubstr(x, "/") ? noLast(x, "/") : x;
        return isSubstr(x_, " ") ? first(x_, " ") : x_;
    }

    template <typename T> std::vector<T> numeric(const std::vector<std::string> &x)
    {
        std::vector<T> r;
        
        for (const auto &i : x)
        {
            if (i != MISSING)
            {
                r.push_back((T) stod(i));
            }
        }
        
        return r;
    }

    template <typename T> std::vector<T> logV(const std::vector<T> &x)
    {
        auto tmp = x;
        for (auto i = 0; i < tmp.size(); i++)
        {
            tmp[i] = log(tmp[i]);
        }
        
        return tmp;
    }
    
    std::vector<double> rank(const std::vector<double> &);
    
    bool isFloat(const std::string &);

    template <typename T> std::vector<T> removeNA(const std::vector<T> &x)
    {
        std::vector<T> tmp;
        
        for (auto &i : x)
        {
            if (!std::isnan(i) && std::isfinite(i))
            {
                tmp.push_back(i);
            }
        }
        
        return tmp;
    }
    
    inline std::string mixToStr(Mixture m)
    {
        switch (m)
        {
            case Mix_1: { return "A"; }
            case Mix_2: { return "B"; }
            case Mix_3: { return "C"; }
        }
    }
    
    template <typename T> std::string toString(const T &x, unsigned n = 2, bool naTrailZero = false)
    {
        if (std::isnan(x) || !std::isfinite(x))
        {
            return MISSING;
        }
        
        std::ostringstream out;
        out << std::fixed << std::setprecision(n) << x;

        auto str = out.str();
        
        if (isEnd(str, "0") && isSubstr(str, ".") && naTrailZero)
        {
            int offset{1};
            if (str.find_last_not_of('0') == str.find('.')) { offset = 0; }
            str.erase(str.find_last_not_of('0') + offset, std::string::npos);
        }

        return str;
    }

    template <typename T1, typename T2> std::string getNAFromD(const T1 &x, const T2 &key, unsigned n = 2)
    {
        return (!x.count(key) || std::isnan(x.count(key))) ? MISSING : toString(x.at(key), n);
    }

    #define S0(x) toString(x,0)
    #define S1(x) toString(x,1)
    #define S2(x) toString(x,2)
    #define S3(x) toString(x,3)
    #define S4(x) toString(x,4)
    #define S5(x) toString(x,5)
    #define S6(x) toString(x,6)

    template <typename T> std::string replaceNA(T x, unsigned n = 2)
    {
        return (std::isnan(x) ? "NA" : toString(x, n));
    }

    template <typename T1, typename T2> std::vector<T2> toVector(const std::map<T1, T2> &m)
    {
        std::vector<T2> x;
        
        for (const auto &i : m)
        {
            x.push_back(i.second);
        }
        
        return x;
    }

    template <typename T> unsigned count(const std::map<T, unsigned> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), 0, [](unsigned c, const std::pair<T, unsigned>& p)
        {
            return c + (p.second ? 1 : 0);
        });
    }
    
    template <typename X, typename F> unsigned countMap(const X &x, F f)
    {
        unsigned n = 0;
        
        for (const auto &i : x)
        {
            n += f(i.first, i.second);
        }
        
        return n;
    }
    
    template <typename T> T sum(const std::vector<T> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), T(), [](T c, const T &p)
        {
            return c + p;
        });
    }

    inline Label calib2str(Proportion p)
    {
        return p == NO_CALIBRATION ? MISSING : p > 1 ? S0(p) : toString(p, 10, true);
    }

    inline FileName path2file(const Path &path)
    {
        auto tmp = path;
        const auto last = path.find_last_of("\\/");
        
        if (std::string::npos != last)
        {
            tmp.erase(0, last + 1);
        }
        
        return tmp;
    }
    
    Path execPath();    
    
    void mergeFQ (const std::vector<FileName> &,
                  const std::vector<FileName> &,
                  const FileName &,
                  const FileName &);

    inline std::vector<FileName> path2file(const std::vector<FileName> &files)
    {
        auto tmp = files;
        
        for (auto i = 0u; i < tmp.size(); i++)
        {
            tmp[i] = path2file(tmp[i]);
        }
        
        return tmp;
    }

    template <typename T1, typename T2> T2 sum(const std::map<T1, T2> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), T2(), [](T2 c, const std::pair<T1, T2>& p)
        {
            return c + p.second;
        });
    }

    template <typename Key, typename Value> std::set<Key> keys(const std::map<Key, Value> &x)
    {
        std::set<Key> keys;
        
        for (const auto i: x)
        {
            keys.insert(i.first);
        }
        
        return keys;
    }
    
    struct ParitionCount : public std::map<Bin, Count>
    {
        inline Count total() const { return sum(*this); }
        
        // Count for each partition
        inline Count binN(Bin x) const { return this->count(x) ? this->at(x) : 0; }
        
        // Proportion for each partition
        inline Proportion binP(Bin x) const { return this->count(x) ? ((Proportion) this->at(x) / total()) : 0.0; }
        
        inline Proportion dil() const { return 1.0 - binP(Bin::ES); }
    };
}

#endif
