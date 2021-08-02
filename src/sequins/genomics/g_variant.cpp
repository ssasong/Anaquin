#include "writers/vcf_writer.hpp"
#include "sequins/genomics/g_variant.hpp"

using namespace Anaquin;

enum Mode
{
    HumanVCF = 1,
    DecoyVCF = 2,
    HumanAndDecoy = 3
};

const auto muts = std::set<Variation>
{
    Variation::SNP,
    Variation::Deletion,
    Variation::Insertion,
};

inline bool isShort(const SequinID &x)
{
    return isSubstr(x, "GV_") ||
           isSubstr(x, "DV_") ||
           isSubstr(x, "HP_") ||
           isSubstr(x, "CV_") ||
           isSubstr(x, "CM_") ||
           isSubstr(x, "CL_") ||
           isSubstr(x, "MS_");
}

static void analyzeMode(Mode mode, const GVariant &v, GVariant::Stats &stats, const FileName &file, const GVariant::Options &o, std::shared_ptr<VCFLadder> vl)
{
    auto w1 = (HumanVCF & mode) ? std::shared_ptr<VCFWriter>(new VCFWriter()) : nullptr;
    auto wT = (DecoyVCF & mode) ? std::shared_ptr<VCFWriter>(new VCFWriter()) : nullptr;
    auto wF = (DecoyVCF & mode) ? std::shared_ptr<VCFWriter>(new VCFWriter()) : nullptr;

    if (w1) { w1->open(o.es); }
    if (wT) { wT->open(o.tp); }
    if (wF) { wF->open(o.fp); }

    const auto &r = Standard::instance().gen;

    const auto r2 = r.r2()->inters(); // Human regions including edges
    const auto r3 = r.r3()->inters(); // Decoy regions without edges
    const auto r4 = r.r4()->inters(); // Decoy regions including edges

    // We'll need it for false negatives
    std::set<long> keys;

    o.analyze(file);
    
    ParserVCF::parse(file, [&](const Variant &x)
    {
        DInter *d = nullptr;

        if (mode & HumanVCF)
        {
            assert(w1);
            
            if ((d = contains(r2, x.cID, x.l)) && isShort(d->name()))
            {
                w1->write(x.hdr, x.line);
                auto t = x;
                t.name = d->name();
                stats.hs.vs.insert(t);
            }
        }
        
        if (!d && (mode & DecoyVCF) && (d = contains(r4, x.cID, x.l)))
        {
            auto findMatch = [&](const Variant &q)
            {
                VCFMatch m;
                
                m.qry = q;
                m.var = nullptr;
                
                // Can we match by position?
                if ((m.var = r.v1()->data.findVar(q.cID, q.l)))
                {
                    // Match by reference allele?
                    m.ref = m.var->ref == q.ref;
                    
                    // Match by alternative allele?
                    m.alt = m.var->alt == q.alt;
                }
                
                if (m.var)
                {
                    m.rID = m.var->name;
                    assert(!m.rID.empty());
                }
                else
                {
                    GIntervals<> inters;
                    
                    try
                    {
                        // Search where the FPs are (finding shoudl ignore edges)
                        inters = r.r3()->ginters(x.cID);
                        
                        const auto m2 = inters.contains(x.l);
                        
                        // Can we find the corresponding region for the FP?
                        if (m2)
                        {
                            m.rID = m2->name();
                            assert(!m.rID.empty());
                        }
                    }
                    catch (...) {}
                }
                
                return m;
            };
            
            auto m = findMatch(x);
            
            // No sequin variant found? It must be a false positive.
            if (!m.var)
            {
                wF->write(x.hdr, x.line);
                stats.ds.fps.push_back(m);
                return;
            }
            
            assert(m.var);
            
            if (!m.rID.empty() && !v.isValid(m.rID))
            {
                return;
            }
            
            // Matched if both the position and alleles agree
            const auto matched = m.var && m.ref && m.alt;
            
            if (matched)
            {
                keys.insert(x.key());
                wT->write(x.hdr, x.line);
                stats.ds.tps.push_back(m);
                assert(!std::isnan(r.af(m.var->name)));
            }
            else
            {
                wF->write(x.hdr, x.line);
                stats.ds.fps.push_back(m);
            }
        }
    });
    
    if (DecoyVCF & mode)
    {
        auto forTP = [&]()
        {
            for (auto i = 0u; i < stats.ds.tps.size(); i++)
            {
                // Overall performance
                stats.ds.oc.tp()++;
            }
        };
        
        auto forFP = [&]()
        {
            for (auto i = 0u; i < stats.ds.fps.size(); i++)
            {
                // Overall performance
                stats.ds.oc.fp()++;
            }
        };
        
        forTP();
        forFP();
        
        for (auto &mut : muts)
        {
            stats.ds.oc.nr() += r.nType(vl, mut);
        }
        
        stats.ds.oc.fn() = stats.ds.oc.nr() - stats.ds.oc.tp();
        assert(stats.ds.oc.nr() >= stats.ds.oc.fn());
        
        /*
         * Work out false negatives
         */
        
        for (const auto &i : r.v1()->data.vars())
        {
            if (v.isValid(i.name) && !stats.ds.findTP(i.name))
            {
                VCFMatch m;
                
                m.var = r.v1()->data.findVar(i.cID, i.l);
                m.rID = i.name;
                assert(m.var);
                
                if (!keys.count(m.var->key()))
                {
                    stats.ds.fns.push_back(m);
                }
            }
        }
    }
    
    writeFN(o.fn, stats.ds.fns, o);
}

std::string GVariant::TableRow::broad() const
{
    std::stringstream ss;
    
    // TP; FN; FP; SN; PC; TPMQ; FPMQ
    ss << tp << "; "
       << fn << "; "
       << fp << "; "
       << S3(sn()) << "; "
       << S3(pc()) << "; "
       << (tpq.empty() ? MISSING : S0(med(tpq))) << "; "
       << (fpq.empty() ? MISSING : S0(med(fpq)));

    return ss.str();
}

GVariant::TableRow GVariant::getTRow(const std::vector<Label> &c,
                                     const std::vector<Label> &x,
                                     const GVariant::Options &o,
                                     bool keep)
{
    assert(c.size() == x.size());
    
    const auto src = o.tsvS;
    assert(!src.empty());    
    auto dst = src;
    
    TableRow t;

    for (auto i = 0; i < c.size(); i++)
    {
        if (!RHead(src, c[i]))
        {
            t.valid = false;
            return t;
        }

        const auto tmp = tmpFile();
        RGrep(dst, tmp, c[i], x[i], keep);
        
        // Prepare for next iteration
        dst = tmp;
    }

    t.tp = RCount(dst, "LABEL", "TP");
    t.fp = RCount(dst, "LABEL", "FP");
    t.fn = RCount(dst, "LABEL", "FN");
    t.nr = t.tp + t.fn;

    if (t.nr + t.fp == 0)
    {
        t.valid = false;
        return t;
    }
    
    t.valid = true;
    
    t.size = RSum(dst, "SIZE");
    t.sample = RCount(dst, "LABEL", "SV");
    t.depth  = (RHead(dst, "REF_DEPTH")       ? RSum(dst, "REF_DEPTH")       : 0) +
               (RHead(dst, "VAR_DEPTH")       ? RSum(dst, "VAR_DEPTH")       : 0) +
               (RHead(dst, "REF_DEPTH_TUMOR") ? RSum(dst, "REF_DEPTH_TUMOR") : 0) +
               (RHead(dst, "VAR_DEPTH_TUMOR") ? RSum(dst, "VAR_DEPTH_TUMOR") : 0);
    
    t.qs = numeric<double>(RReadTSV(dst, "QUAL"));
    
    if (RHead(dst, "REF_DEPTH_TUMOR"))
    {
        t.ds = RAdd(dst, "REF_DEPTH_TUMOR", "VAR_DEPTH_TUMOR");
    }
    else
    {
        t.ds = RAdd(dst, "REF_DEPTH", "VAR_DEPTH");
    }

    if (RHead(dst, "OBS_FREQ_TUMOR"))
    {
        t.af = numeric<double>(RReadTSV(dst, "OBS_FREQ_TUMOR"));
    }
    else if (RHead(dst, "OBS_FREQ"))
    {
        t.af = numeric<double>(RReadTSV(dst, "OBS_FREQ"));
    }

    const auto tmp1 = tmpFile();
    const auto tmp2 = tmpFile();
    
    RGrep(dst, tmp1, "LABEL", "TP"); // TP
    RGrep(dst, tmp2, "LABEL", "FP"); // FP

    t.tpq = numeric<double>(RReadTSV(tmp1, "QUAL"));
    t.fpq = numeric<double>(RReadTSV(tmp2, "QUAL"));
    
    return t;
}

Base GVariant::countSizeForG(const Options &o)
{
    Base x = 0;
    
    // Trimmed (for decoy and non-decoy)
    auto r = o.decoy ? *(Standard::instance().gen.r4()) : *(Standard::instance().gen.r2());
    
    for (auto &i : r)
    {
        for (auto &j : i.second.l2d)
        {
            x += j.second.l.length(); // Everything...
        }
    }
    
    assert(x);
    return x;
}

GVariant::Stats GVariant::analyze(Stats &stats, const GVariant &x, const FileName &f1, const FileName &f2, const Options &o, std::shared_ptr<VCFLadder> vl)
{
    assert(!o.es.empty());
    assert(!o.fp.empty());
    assert(!o.tp.empty());
    assert(!o.fn.empty());

    o.info("Edge: " + S0(o.edge));
    o.info("Variant: " + Standard::instance().gen.v1()->src);

    if (o.decoy)
    {
        assert(f2.empty());
        analyzeMode(HumanAndDecoy, x, stats, f1, o, vl);
    }
    else
    {
        if (!f1.empty()) { analyzeMode(HumanVCF, x, stats, f1, o, vl); }
        analyzeMode(DecoyVCF, x, stats, f2, o, vl);
    }

    o.info("TP: " + S0(stats.ds.oc.tp()));
    o.info("FP: " + S0(stats.ds.oc.fp()));
    o.info("FN: " + S0(stats.ds.fns.size()));

    return stats;
}
