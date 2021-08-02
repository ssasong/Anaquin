#ifndef B_DATA_HPP
#define B_DATA_HPP

#include "tools/tools.hpp"
#include "data/ginters.hpp"
#include "data/dinters.hpp"
#include "parsers/parser_bed.hpp"

namespace Anaquin
{
    typedef ParserBed::Data Data;

    struct BedChrData
    {
        std::map<Locus, Data> l2d;
    };
    
    struct BedData : public std::map<ChrID, BedChrData>
    {
        inline const ParserBed::Data * find(const Name &x, bool exact = true) const
        {
            for (const auto &i : *this)
            {
                for (const auto &j : i.second.l2d)
                {
                    if ((exact && j.second.name == x) || (!exact && isSubstr(j.second.name, x)))
                    {
                        return &j.second;
                    }
                }
            }
            
            return nullptr;
        }

        inline DIntervals<> inters(const ChrID &x) const
        {
            DIntervals<> r;
            
            for (const auto &i : at(x).l2d)
            {
                r.add(DInter(i.second.name, i.second.l));
            }
            
            r.build();
            return r;
        }

        inline Chr2DInters inters() const
        {
            Chr2DInters r;
            
            for (const auto &i : *this)
            {
                r[i.first] = inters(i.first);
            }
            
            return r;
        }
        
        inline const ParserBed::Data * overlap(const ChrID &x, const Locus &l)
        {
            for (const auto &i : *this)
            {
                if (i.first == x)
                {
                    for (const auto &j : i.second.l2d)
                    {
                        if (j.first.overlap(l))
                        {
                            return &j.second;
                        }
                    }
                }
            }
            
            return nullptr;
        }

        inline GIntervals<> ginters(const ChrID &x) const
        {
            GIntervals<> r;
            
            for (const auto &i : at(x).l2d)
            {
                r.add(GInterval(x, i.second.name, i.second.l));
            }
            
            r.build();
            return r;
        }
        
        inline std::map<ChrID, GIntervals<>> ginters() const
        {
            std::map<ChrID, GIntervals<>> r;
            
            for (const auto &i : *this)
            {
                r[i.first] = ginters(i.first);
            }
            
            return r;
        }
        
        // Source for regions (if available)
        FileName src;
    };
    
    struct RegionOptions
    {
        RegionOptions(Base edge = 0) : edge(edge) {}
        
        Base edge = 0;
        
        // Only those chromosomes?
        std::set<ChrID> onlyC;
        
        // Only those sequins?
        std::set<SequinID> only;
        
        // Any sequin containing the substrings?
        std::set<std::string> conts;
    };
    
    template <typename F> BedData readRegions(const Reader &r, F f, RegionOptions o = RegionOptions())
    {
        BedData c2d;
        
        ParserBed::parse(r, [&](ParserBed::Data &x, Progress i)
        {
            auto edge = o.edge;
            
            if (x.l.length() < 2.0 * o.edge)
            {
                edge = 0.0; // Ignore edge ...
            }
            
            /*
             * Filter out unwanted regions if we have the information (fourth column in BED).
             */
            
            if (!o.only.empty() && !o.only.count(x.name) && !std::any_of(o.conts.begin(), o.conts.end(), [&](const SequinID &i)
            {
                return isSubstr(x.name, i);
            }))
            {
                return;
            }            
            else if (!o.onlyC.empty() && !o.onlyC.count(x.cID))
            {
                return;
            }
            
            x.l.end   -= edge;
            x.l.start += edge;
            
            if (x.l.end < x.l.start)
            {
                throw std::runtime_error("Region: " + x.name + " has an error. The end position is equal or before the start position. The start position is: " + std::to_string(x.l.start) + " and the end position is: " + std::to_string(x.l.end));
            }
            
            c2d[x.cID].l2d[x.l] = x;
            f(x, i);
        });

        c2d.src = r.src();
        assert(!c2d.src.empty());
        
        return c2d;
    }
    
    inline BedData readRegions(const Reader &r, RegionOptions o = RegionOptions())
    {
        return readRegions(r, [](const ParserBed::Data &, Progress) {}, o);
    }
}

#endif
