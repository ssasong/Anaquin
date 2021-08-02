#ifndef G_SOMATIC_HPP
#define G_SOMATIC_HPP

#include "stats/linear.hpp"
#include "data/standard.hpp"
#include "data/analyzer.hpp"
#include "data/resources.hpp"
#include "sequins/genomics/g_variant.hpp"

namespace Anaquin
{
    template <typename T> std::string extra(const std::vector<std::string> &hs, const T &x, long key)
    {
        std::stringstream ss;
        
        for (auto &h : hs)
        {
            if (x.si.count(h))
            {
                if (x.si.at(h).count(key)) { ss << ("\t" + S2(x.si.at(h).at(key))); }
                else                       { ss << "\t-"; }
            }
            else if (x.sf.count(h))
            {
                if (x.sf.at(h).count(key)) { ss << ("\t" + S2(x.sf.at(h).at(key))); }
                else                       { ss << "\t-"; }
            }
            else
            {
                ss << "\t-";
            }
        }
        
        return ss.str();
    }

    template <typename S, typename T, typename O> void writeFP(const S &,
                                                               const std::vector<T> &x,
                                                               const std::string &label,
                                                               const std::string &format,
                                                               const O &o)
    {
        const auto &r = Standard::instance().gen;

        for (const auto &i : x)
        {
            if (i.rID.empty())
            {
                continue;
            }
                
            auto sID = (i.var && i.alt && i.ref ? i.var->name : MISSING);
            
            #define _FI_(x) (i.qry.fi.count(x) ? S2(i.qry.fi.at(x)) : MISSING)
            #define _FF_(x) (i.qry.ff.count(x) ? S2(i.qry.ff.at(x)) : MISSING)
            
            o.writer->write((boost::format(format) % i.rID
                                                   % i.qry.cID
                                                   % i.qry.l.start
                                                   % label
                                                   % i.qry.alt
                                                   % var2str(i.qry.type())
                                                   % (sID != MISSING ? std::to_string(r.af(sID)) : MISSING)
                                                   % normalDPR(i.qry)
                                                   % normalDPV(i.qry)
                                                   % (std::isnan(normalAF(i.qry)) ? MISSING : S4(normalAF(i.qry)))
                                                   % tumorDPR(i.qry)
                                                   % tumorDPV(i.qry)
                                                   % (std::isnan(tumorAF(i.qry)) ? MISSING : S4(tumorAF(i.qry)))
                                                   % MISSING
                                                   % S2(i.qry.qual[0])
                                                   % r.a1()->strForLocus(i.qry.cID, i.qry.l)).str());
        }
    }
    
    template <typename T, typename O> void writeSD(const std::string &format, const T &x, const O &o)
    {
        writeFP(x.ds, x.ds.fps, "FP", format, o);
    }

    template <typename T, typename O> void writeSS(const std::string &format, const T &stats, const O &o)
    {
        const auto &r = Standard::instance().gen;

        for (const auto &i : stats.hs.vs)
        {
            o.writer->write((boost::format(format) % i.name
                                                   % i.cID
                                                   % i.l.start
                                                   % "SV"
                                                   % i.alt
                                                   % var2str(i.type())
                                                   % MISSING
                                                   % normalDPR(i)
                                                   % normalDPV(i)
                                                   % (std::isnan(normalAF(i)) ? MISSING : S4(normalAF(i)))
                                                   % tumorDPR(i)
                                                   % tumorDPV(i)
                                                   % (std::isnan(tumorAF(i))  ? MISSING : S4(tumorAF(i)))
                                                   % MISSING
                                                   % S2(i.qual[0])
                                                   % r.a1()->strForLocus(i.cID, i.l)).str());
        }
        
        o.writer->close();
    }

    template <typename T, typename O> void writeSQ(const std::string &format, const T &stats, const O &o)
    {
        const auto &r = Standard::instance().gen;
        const auto r4 = r.r4()->ginters();
        const auto r5 = r.r5()->ginters();

        for (const auto &i : r.v1()->data.vars())
        {
            if (isSoma(i.name))
            {
                // Can we find this sequin?
                const auto isTP = stats.ds.findTP(i.name);
    
                // Locus for the sequin
                const auto l5 = r5.at(i.cID).overlap(i.l)->l();

                #define FORMAT_I(x) (c.fi.count(x) ? S2(c.fi.at(x)) : MISSING)
                #define FORMAT_F(x) (c.ff.count(x) ? S2(c.ff.at(x)) : MISSING)
    
                if (isTP)
                {
                    // Called variant (if found)
                    const auto &c = isTP->qry;
                    
                    o.writer->write((boost::format(format) % i.name
                                                           % i.cID
                                                           % i.l.start
                                                           % "TP"
                                                           % i.alt
                                                           % var2str(i.type())
                                                           % r.af(i.name)
                                                           % normalDPR(c)
                                                           % normalDPV(c)
                                                           % (std::isnan(normalAF(c)) ? MISSING : S4(normalAF(c)))
                                                           % tumorDPR(c)
                                                           % tumorDPV(c)
                                                           % (std::isnan(tumorAF(c))  ? MISSING : S4(tumorAF(c)))
                                                           % r4.at(i.cID).length(l5)
                                                           % S2(isTP->qry.qual[0])
                                                           % r.a1()->strForLocus(i.cID, i.l)).str());
                }
                else
                {
                    o.writer->write((boost::format(format) % i.name
                                                           % i.cID
                                                           % i.l.start
                                                           % "FN"
                                                           % MISSING
                                                           % var2str(i.type())
                                                           % r.af(i.name)
                                                           % MISSING
                                                           % MISSING
                                                           % MISSING
                                                           % MISSING
                                                           % MISSING
                                                           % MISSING
                                                           % r4.at(i.cID).length(l5)
                                                           % MISSING
                                                           % r.a1()->strForLocus(i.cID, i.l)).str());
                }
            }
        }
    }

    template <typename T, typename O> void writeSomatic(const FileName &file, const T &stats, const O &o)
    {
        const auto &r = Standard::instance().gen;
        const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9$.5f\t%10%\t%11%\t%12%\t%13%\t%14%\t%15%%16%";
        
        o.writer->open(file);
        o.writer->write((boost::format(format) % "NAME"
                                               % "CHROM"
                                               % "POSITION"
                                               % "LABEL"
                                               % "ALLELE"
                                               % "TYPE"
                                               % "EXP_FREQ"
                                               % "REF_DEPTH_NORMAL"
                                               % "VAR_DEPTH_NORMAL"
                                               % "OBS_FREQ_NORMAL"
                                               % "REF_DEPTH_TUMOR"
                                               % "VAR_DEPTH_TUMOR"
                                               % "OBS_FREQ_TUMOR"
                                               % "SIZE"
                                               % "QUAL"
                                               % r.a1()->strForKeys()).str());
        
        writeSQ(format, stats, o);
        writeSD(format, stats, o);
        writeSS(format, stats, o);
        
        o.writer->close();
    }

    /*
     * Strelka implementation
     */
    
    inline bool isStrelka(const Variant &x)
    {
        auto isSNP = [&]()
        {
            return x.fi.count("AU_2_1") &&
                   x.fi.count("CU_2_1") &&
                   x.fi.count("GU_2_1") &&
                   x.fi.count("TU_2_1") &&
                   x.fi.count("AU_2_2") &&
                   x.fi.count("CU_2_2") &&
                   x.fi.count("GU_2_2") &&
                   x.fi.count("TU_2_2");
        };
        
        auto isInd = [&]()
        {
            return x.fi.count("TAR_1_1") &&
                   x.fi.count("TAR_1_2") &&
                   x.fi.count("TAR_2_1") &&
                   x.fi.count("TAR_2_2") &&
                   x.fi.count("TIR_1_1") &&
                   x.fi.count("TIR_1_2") &&
                   x.fi.count("TIR_2_1") &&
                   x.fi.count("TIR_2_2");
        };
        
        return isSNP() || isInd();
    }
    
    inline Count strelkaNormalT(const Variant &x, const Sequence &a)
    {
        switch (x.type())
        {
            case SNP:
            {
                if (a == "A") { return x.fi.at("AU_1_1"); }
                if (a == "C") { return x.fi.at("CU_1_1"); }
                if (a == "G") { return x.fi.at("GU_1_1"); }
                else          { return x.fi.at("TU_1_1"); }
                break;
            }
                
            case Deletion:
            case Insertion:
            {
                if (a == x.ref)      { return x.fi.at("TAR_1_1"); }
                else if (a == x.alt) { return x.fi.at("TIR_1_1"); }
                else { throw std::runtime_error("Unknown: " + a); }
            }

            default: { throw std::runtime_error("Not supported in strelkaNormalT()"); }
        }
    }
    
    inline Count strelkaTumorT(const Variant &x, const Sequence &a)
    {
        switch (x.type())
        {
            case SNP:
            {
                if (a == "A") { return x.fi.at("AU_2_1"); }
                if (a == "C") { return x.fi.at("CU_2_1"); }
                if (a == "G") { return x.fi.at("GU_2_1"); }
                else          { return x.fi.at("TU_2_1"); }
                break;
            }
                
            case Deletion:
            case Insertion:
            {
                if (a == x.ref)      { return x.fi.at("TAR_2_1"); }
                else if (a == x.alt) { return x.fi.at("TIR_2_1"); }
                else { throw std::runtime_error("Unknown: " + a); }
            }

            default: { throw std::runtime_error("Not supported in strelkaTumorT()"); }
        }
    }
    
    // Depth for reference allele in normal sample
    inline Count strelkaNormalDPR(const Variant &x)
    {
        return strelkaNormalT(x, x.ref);
    }
    
    // Depth for variant allele in normal sample
    inline Count strelkaNormalDPV(const Variant &x)
    {
        return strelkaNormalT(x, x.alt);
    }
    
    // Depth for reference allele in tumor sample
    inline Count strelkaTumorDPR(const Variant &x)
    {
        return strelkaTumorT(x, x.ref);
    }
    
    // Depth for variant allele in tumor sample
    inline Count strelkaTumorDPV(const Variant &x)
    {
        return strelkaTumorT(x, x.alt);
    }
    
    inline Proportion strelkaNormalAF(const Variant &x)
    {
        const auto r = strelkaNormalDPR(x);
        const auto v = strelkaNormalDPV(x);
        return ((Proportion) v) / (r + v);
    }
    
    inline Proportion strelkaTumorAF(const Variant &x)
    {
        const auto r = strelkaTumorDPR(x);
        const auto v = strelkaTumorDPV(x);
        return ((Proportion) v) / (r + v);
    }
    
    // Depth for reference allele in normal sample
    inline Proportion normalDPR(const Variant &x)
    {
        if (isStrelka(x)) { return strelkaNormalDPR(x); }
        else              { return x.AD[0] == -1 ? -1 : x.AD[0]; }
    }
    
    // Depth for variant allele in normal sample
    inline Proportion normalDPV(const Variant &x)
    {
        if (isStrelka(x)) { return strelkaNormalDPV(x); }
        else              { return x.AD[1] == -1 ? -1 : x.AD[1]; }
    }

    // Depth for reference allele in tumor sample
    inline Proportion tumorDPR(const Variant &x)
    {
        if (isStrelka(x)) { return strelkaTumorDPR(x); }
        else              { return x.AD[2] == -1 ? -1 : x.AD[2]; }
    }
    
    // Depth for variant allele in tumor sample
    inline Proportion tumorDPV(const Variant &x)
    {
        if (isStrelka(x)) { return strelkaTumorDPV(x); }
        else              { return x.AD[3] == -1 ? -1 : x.AD[3]; }
    }
    
    // Measured allele frequency for normal
    inline Proportion normalAF(const Variant &x)
    {
        if (x.ff.count("AF_1"))    { return x.ff.at("AF_1");    }
        else if (isStrelka(x))     { return strelkaNormalAF(x); }
        else
        {
            const auto r = normalDPR(x);
            const auto v = normalDPV(x);
            return ((Proportion) v) / (r + v);
        }
    }
    
    // Measured allele frequency for tumor
    inline Proportion tumorAF(const Variant &x)
    {
        if (x.ff.count("AF_2"))    { return x.ff.at("AF_2");   }
        else if (x.ff.count("AF")) { return x.ff.at("AF");     }
        else if (isStrelka(x))     { return strelkaTumorAF(x); }
        else
        {
            const auto r = tumorDPR(x);
            const auto v = tumorDPV(x);
            return ((Proportion) v) / (r + v);
        }
    }

    struct GSomatic : public GVariant
    {
        struct Stats : public GVariant::Stats
        {
            /*
             * Statistics for allele frequency
             */
            
            struct AlleleStats : public SequinStats, public LimitStats {};
            
            // Performance for each variation
            std::map<Variation, AlleleStats> m2a;
            
            // Overall performance
            AlleleStats oa;

            /*
             * Performance by allele frequency
             */
            
            std::map<float, Confusion> af;
        };
        
        bool isValid(const SequinID &x) const override { return isSoma(x); }

        static void report(const FileName &, const FileName &, const GVariant::Options &);
    };
}

#endif
