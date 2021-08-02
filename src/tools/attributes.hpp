#ifndef ATTRIBUTES_HPP
#define ATTRIBUTES_HPP

#include <set>
#include <map>
#include <sstream>
#include "data/bData.hpp"
#include "tools/tools.hpp"

namespace Anaquin
{
    typedef std::string AttKey;
    typedef std::string AttGrp;
    typedef std::string AttStr;

    inline AttKey attrKey(const AttStr &x) { return first(x, "_"); }
    inline AttGrp attrGrp(const AttStr &x) { return last(x, "_");  }
    
    class AttributeBed : public std::map<AttKey, std::shared_ptr<BedData>>
    {
        public:

            inline const std::set<AttKey> &keys() const
            {
                if (_keys.empty())
                {
                    for (const auto &i : (*this))
                    {
                        _keys.insert(i.first);
                    }
                }
            
                return _keys;
            }

            inline std::string strForKeys() const
            {
                std::stringstream ss;
                for (const auto &i : (*this)) { ss << "\t" << i.first; }            
                return ss.str();
            }

            inline std::string strForNulls() const
            {
                if (_nulls.empty())
                {
                    std::stringstream ss;
                    for (auto i = 0u; i < keys().size(); i++) { ss << "\t"; ss << MISSING; }
                    _nulls = ss.str();
                }
                
                return _nulls;
            }

            static std::string strForVals(const std::map<AttKey, AttGrp> &am)
            {
                std::stringstream ss;
                for (const auto &i : am) { ss << "\t"; ss << i.second; }
                return ss.str();
            }
        
            inline std::string strForLocus(const ChrID &cID, const Locus &l)
            {
                std::stringstream ss;
                
                for (const auto &i : keys())
                {
                    ss << "\t"; ss << strForLocus(i, cID, l);
                }
                
                return ss.str();
            }
        
            inline bool overlap(const AttKey &key, const ChrID &cID, const Locus &l)
            {
                initCache(key);
                return _cache[key].count(cID) && _cache[key].at(cID).overlap(l);
            }
        
            inline std::string strForLocus(const AttKey &key, const ChrID &cID, const Locus &l)
            {
                initCache(key);
                
                const auto &m = _cache[key];
                std::vector<GInterval *> o;
            
                if (m.count(cID) && m.at(cID).overlap(l, &o))
                {
                    std::vector<std::string> os;
                
                    std::transform(o.begin(), o.end(), std::back_inserter(os), [&](GInterval *i)
                    {
                        return noFirst(i->name(), "_");
                    });

                    return os.back();
                }
            
                return MISSING;
            }
        
            FileName src;
            std::map<AttKey, std::set<AttGrp>> vals;

        private:
        
            void initCache(const AttKey &key)
            {
                if (!_cache.count(key))
                {
                    _cache[key] = (*this)[key]->ginters();
                }
            }
        
            mutable std::string _nulls;
            mutable std::set<AttKey> _keys;
            mutable std::map<AttKey, std::map<ChrID, GIntervals<>>> _cache;
    };
    
    AttributeBed readAttrib(const FileName &);
}

#endif