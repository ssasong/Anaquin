#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <map>
#include <cmath>
#include "data/locus.hpp"
#include <boost/format.hpp>

namespace Anaquin
{
    inline long var2hash(const SequinID &id, Variation type, const Locus &l)
    {
        const auto str = (boost::format("%1%_%2%_%3%_%4%") % id
                                                           % type
                                                           % l.start
                                                           % l.end).str();
        return std::hash<std::string>{}(str);
    }

    struct Variant
    {
        inline bool operator<(const Variant &x) const
        {
            if (parseChrID(cID) < parseChrID(x.cID))
            {
                return true;
            }
            else if (parseChrID(cID) > parseChrID(x.cID))
            {
                return false;
            }
            
            return l < x.l;
        }
        
        inline bool isSV() const
        {
            return alt == "<DUP>" || alt == "<INV>" || alt == "<INS>" || alt == "<DEL>";
        }
        
        inline Variation type() const
        {
            if (alt == "<DUP>" || alt[0] == '=')
            {
                return Variation::Duplication;
            }
            else if (alt == "<INV>" || alt[0] == '^')
            {
                return Variation::Inversion;
            }
            else if (alt == "<INS>" || alt[0] == '+')
            {
                return Variation::Insertion;
            }
            else if (alt == "<DEL>" || alt[0] == '-')
            {
                return Variation::Deletion;
            }
            else if (ref.size() == alt.size())
            {
                return Variation::SNP;
            }
            else if (ref.size() > alt.size())
            {
                return Variation::Deletion;
            }
            else
            {
                return Variation::Insertion;
            }
        }

        inline long key() const
        {
            return var2hash(name, type(), l);
        }
        
        /*
         * Convenient functions
         */
        
        // Depth for a specified allele
        inline std::string dp(unsigned i) const
        {
            return AD[i] >= 0 ? S2(AD[i]) : MISSING;
        }
        
        // Observed allele frequency
        inline std::string obsAF() const
        {
            if (AD[0] == -1 && AD[1] >= 0)
            {
                return S2(1.00);
            }
            else if (AD[0] >= 0 && AD[1] == -1)
            {
                return S2(0.00);
            }
            
            return S2((Proportion) AD[1] / (AD[0] + AD[1]));
        }
        
        // SNP type. Undefined if not a SNP.
        inline SNPType snpType() const
        {
            assert(type() == Variation::SNP);
            
            if (ref == "A" && alt == "C") { return SNPType::AC; }
            if (ref == "A" && alt == "T") { return SNPType::AT; }
            if (ref == "A" && alt == "G") { return SNPType::AG; }

            if (ref == "C" && alt == "A") { return SNPType::CA; }
            if (ref == "C" && alt == "T") { return SNPType::CT; }
            if (ref == "C" && alt == "G") { return SNPType::CG; }

            if (ref == "T" && alt == "A") { return SNPType::TA; }
            if (ref == "T" && alt == "G") { return SNPType::TG; }
            if (ref == "T" && alt == "C") { return SNPType::TC; }

            if (ref == "G" && alt == "A") { return SNPType::GA; }
            if (ref == "G" && alt == "C") { return SNPType::GC; }
            if (ref == "G" && alt == "T") { return SNPType::GT; }

            if (ref == alt) { return SNPType::RF; }

            throw std::runtime_error("Unknown SNP");
        }

        ChrID cID;

        // Eg: GI_005
        SequinID name;

        // The reference position, with the 1st base having position 1
        Locus l;
        
        Genotype gt;
        
        // Reference and alternative allele
        Sequence ref, alt;

        // Allelle frequency
        Proportion allF = NAN;
        
        // Quality score
        float qual[2] = { NAN, NAN };
        
        // Depth for reference and allele (e.g. normal and tumor)
        int AD[4] = {-1,-1,-1,-1};

        std::map<std::string, int> ifi;
        std::map<std::string, float> iff;
        std::map<std::string, std::string> ifs;
        
        std::map<std::string, int> fi;
        std::map<std::string, float> ff;
        std::map<std::string, std::string> fs;

        void *hdr, *line;
    };
}

#endif
