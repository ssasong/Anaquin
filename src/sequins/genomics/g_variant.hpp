#ifndef G_VARIANT_HPP
#define G_VARIANT_HPP

#include "data/analyzer.hpp"
#include "writers/vcf_writer.hpp"
#include "sequins/genomics/genomics.hpp"

namespace Anaquin
{
    struct VCFMatch
    {
        // The called variant
        Variant qry;
        
        // Sequin matched by position?
        const Variant *var = nullptr;
        
        // Matched by variant allele? Only if position is matched.
        bool alt;
        
        // Matched by reference allele? Only if position is matched.
        bool ref;
        
        // Does the variant fall into one of the reference regions?
        SequinID rID;
    };

    template <typename T, typename O> void writeFN(const FileName &file, const T &g, const O &o)
    {
        o.generate(file);
        
        auto wFN = VCFWriter();
        wFN.open(file);
        
        std::set<long> keys;
        for (const auto &i : g) { keys.insert(i.var->key()); }
        
        ParserVCF::parse(Standard::instance().gen.v1()->src, [&](const Variant &x)
        {
            if (keys.count(x.key()))
            {
                wFN.write(x.hdr, x.line);
            }
        });
    }

    struct GVariant
    {
        struct HStats
        {
            std::set<Variant> vs;
            
            /*
             * Caller specific fields
             */
            
            std::map<std::string, std::map<long, int>>   si;
            std::map<std::string, std::map<long, float>> sf;
        };
        
        struct DStats
        {
            std::vector<VCFMatch> tps, fns, fps;
            
            // Overall performance
            Confusion oc;
            
            inline const VCFMatch * findTP(const SequinID &id) const
            {
                for (auto &i : tps)
                {
                    if (i.var->name == id)
                    {
                        return &i;
                    }
                }
                
                return nullptr;
            }
            
            /*
             * Caller specific fields
             */
            
            std::map<std::string, std::map<long, int>>   si;
            std::map<std::string, std::map<long, float>> sf;
        };
        
        struct Stats
        {
            HStats hs;
            DStats ds;
        };
        
        struct Options : public AnalyzerOptions
        {
            Base edge;
            
            FileName uBED;
          
            bool decoy;
            
            // Generating HTML report?
            bool report;
            
            // Required for common operations
            Label base;
            
            FileName tsvS;
            
            // Files written for sample, TP, FP etc
            FileName es, tp, fp, fn;
        };
        
        virtual bool isValid(const SequinID &) const = 0;
        
        static Stats analyze(Stats &, const GVariant &, const FileName &, const FileName &, const Options &, std::shared_ptr<VCFLadder>);
        
        struct TableRow
        {
            TableRow & operator+(const TableRow &x)
            {
                assert(valid && x.valid);
                this->nr += x.nr;
                this->tp += x.tp;
                this->fp += x.fp;
                this->fn += x.fn;
                this->depth += x.depth;
                this->sample += x.sample;
                
                for (auto &i : x.ds) { this->ds.push_back(i); }
                for (auto &i : x.qs) { this->qs.push_back(i); }

                return *this;
            }
            
            bool valid = false;
            
            Count nr, tp, fp, fn;
            long size;
            
            inline Proportion sn() const { return (Proportion) tp / nr; }
            inline Proportion pc() const { return (Proportion) tp / (tp + fp); }

            inline double fpKB() const
            {
                return size ? (double) fp / size : NAN;
            }
            
            Count depth;  // Combined depth
            Count sample; // Number of samples found
            
            // Quality scores
            std::vector<double> tpq, fpq, qs;
            
            // Allele frequency
            std::vector<double> af;
            
            // Depth coverage
            std::vector<Coverage> ds;
            
            // Generate a report in Broad's style
            std::string broad() const;
        };
        
        // Sequin region size for germline analysis. Excluding edges.
        static Base countSizeForG(const Options &);
        
        static TableRow getTRow(const std::vector<Label> &,
                                const std::vector<Label> &,
                                const GVariant::Options &, bool keep = true);

        static TableRow getTRow(const Label &c1,
                                const Label &x1,
                                const Label &c2,
                                const Label &x2,
                                const GVariant::Options &o, bool keep = true)
        {
            return GVariant::getTRow(std::vector<Label> { c1, c2 }, std::vector<Label> { x1, x2 }, o, keep);
        }

        static TableRow getTRow(const Label &c,
                                const Label &x,
                                const GVariant::Options &o,
                                bool keep = true)
        {
            return GVariant::getTRow(std::vector<Label> { c }, std::vector<Label> { x }, o, keep);
        }

        static TableRow getTotal(const Options &o) { return getTRow("NAME", "", o); }
        static TableRow getTP(const Options &o)    { return getTRow("LABEL", "TP", o); }
        static TableRow getFP(const Options &o)    { return getTRow("LABEL", "FP", o); }

        static TableRow getSNV(const Options &o)   { return getTRow("TYPE", "SNP", o); }
        static TableRow getIndel(const Options &o) { return getTRow("TYPE", "SNP", o, false); }

        static TableRow getHom(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GENOTYPE", "Homozygous", o);
        }
        
        static TableRow getHet(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GENOTYPE", "Heterozygous", o);
        }
        
        static TableRow getCode(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GeneContext", "CodingRegion", o);
        }
        
        static TableRow getNCode(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GeneContext", "NoncodingRegion", o);
        }
        
        static TableRow getACMG(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GeneContext", "CodingRegion_ACMG", o);
        }
                
        static TableRow getPGX(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GeneContext", "CodingRegion_PGx", o);
        }
       
        static TableRow getHighGC(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GCcontent", "GCrich", o);
        }
                
        static TableRow getModGC(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GCcontent", "GCmoderate", o);
        }
                
        static TableRow getLowGC(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "GCcontent", "ATrich", o);
        }
                
        static TableRow getLowConf(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "NistRegions", "LowConf", o);
        }
                
        static TableRow getHighConf(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "NistRegions", "HighConf", o);
        }

        static TableRow getRepeat(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "SimpleRepeat", "NA", o, false);
        }

        static TableRow getShortRepeat(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "SimpleRepeat", "Short", o);
        }

        static TableRow getLongRepeat(const Options &o, const std::string &x = "")
        {
            return getTRow("TYPE", x, "SimpleRepeat", "Long", o);
        }

        static TableRow getSine(const Options &o)   { return getTRow("MobileElement", "SINE", o); }
        static TableRow getDNA(const Options &o)    { return getTRow("MobileElement", "DNA", o);  }
        static TableRow getLine(const Options &o)   { return getTRow("MobileElement", "LINE", o); }
        static TableRow getLTR(const Options &o)    { return getTRow("MobileElement", "LTR", o);  }

        static TableRow getDi(const Options &o)     { return getTRow("SimpleRepeat", "Di", o);   }
        static TableRow getTri(const Options &o)    { return getTRow("SimpleRepeat", "Tri", o);  }
        static TableRow getMono(const Options &o)   { return getTRow("SimpleRepeat", "Mono", o); }
        static TableRow getQuad(const Options &o)   { return getTRow("SimpleRepeat", "Quad", o); }
    };
}

#endif
