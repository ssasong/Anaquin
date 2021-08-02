#include <htslib/bgzf.h>
#include "tools/calibrator.hpp"
#include "writers/fq_writer.hpp"
#include "parsers/parser_fq.hpp"
#include "writers/bam_writer.hpp"
#include "parsers/parser_bam.hpp"
#include "sequins/genomics/genomics.hpp"

using namespace Anaquin;

typedef LadderCalibrator::Method  Method;
typedef LadderCalibrator::Options Options;
typedef LadderCalibrator::Results Results;

void LadderCalibrator::getNorms(const Options &o, Results &r)
{
    switch (o.meth)
    {
        case Method::Pooling:
        {
            assert(o.d2 && !o.d2->tsvL.empty() && o.d2->p != NO_CALIBRATION);

            const auto tmp1 = tmpFile();
            const auto tmp2 = tmpFile();

            RGrep(o.d2->tsvL, tmp1, "NAME", "LD_");
            RFilterC(tmp1, tmp2, std::set<Column> { "NAME", "READ" }, true);
            
            std::map<StandardID, Count> m;
            const auto x = RReadTSV(tmp2, "NAME", "READ");
            
            for (const auto &i : x)
            {
                m[i.first] += stoi(i.second);
            }
            
            // Pool size for all ladder reads
            const auto pool = sum(m);
            
            // Proportion of the pool
            const auto tar = o.d2->p * pool;
            
            o.logInfo("Pool: " + std::to_string(pool));
            o.logInfo("Target: " + std::to_string(tar));
            
            for (const auto &i : m)
            {
                // Discard if less than the target
                const auto scale = i.second < tar ? 1.0 : 1.0 - (tar / i.second);
                
                r.scales.push_back(scale);
                
                // Scale each individual ladder to the target
                r.rnd[i.first] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - scale));
            }
            
            assert(!r.rnd.empty());
            break;
        }
            
        case Method::SampleCNV2:
        {
            o.info("Sample coverage: " + std::to_string(o.d1->sampC));
            assert(!o.d1->tsvL.empty());
            const auto tmp1 = tmpFile();
            const auto tmp2 = tmpFile();

            RFilterC(o.d1->tsvL, tmp1, std::set<Column> { "NAME", "MEAN" }, true);

            for (const auto &name : RReadTSV(tmp1, "NAME"))
            {
                const auto tmp2 = tmpFile();
                RGrep(tmp1, tmp2, "NAME", name);
                const auto seqC = RMean(tmp2, "MEAN");

                // 0.5 because of paired-end reads
                auto p = seqC ? std::min((0.5 * o.d1->sampC) / seqC, 1.0) : 1.0;

                if (std::isnan(p) || std::isinf(p))
                {
                    p = 1; // Keep all reads
                }

                r.scales.push_back(p);
                
                // Scale each individual standard to the target (same factor)
                r.rnd[name] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - p));
            }

            break;
        }
    }
}

struct FQLadderCalibrator : public LadderCalibrator
{
    FQLadderCalibrator(const FileName &f1, const FileName &f2,
                       const FileName &o1, const FileName &o2) : f1(f1), f2(f2), o1(o1), o2(o2) {}
    
    Results calibrate(const Options &o) override
    {
        Results r;        
        assert(o.d2 && !o.d2->tsvL.empty() && !o.d2->reads.empty());
        
        if (o.meth == Method::SampleCNV2)
        {
            throw std::runtime_error("FASTQ ladder calibration doesn't support SampleCNV2");
        }
        else if (o.d2->p == NO_CALIBRATION)            
        {
            return r;
        }
        
        const auto tmp1 = tmpFile();
        const auto tmp2 = tmpFile();

        RGrep(o.d2->tsvL, tmp1, "NAME", "LD_");
        RFilterC(tmp1, tmp2, std::set<Column> { "NAME", "COUNT" }, true);

        std::map<StandardID, Count> m;
        const auto x = RReadTSV(tmp2, "NAME", "COUNT");
        
        for (const auto &i : x)
        {
            m[i.first] += stoi(i.second);
        }
        
        // Pool size for all ladder reads
        const auto pool = sum(m);
        
        // Proportion of the pool
        const auto targetC = o.d2->p * pool;
        
        o.logInfo("Ladder calibration (pool size): " + std::to_string(pool));
        o.logInfo("Ladder calibration (target): " + std::to_string(targetC));

        for (const auto &i : m)
        {
            // Discard if less than the target
            const auto scale = i.second < targetC ? 1.0 : 1.0 - (targetC / i.second);
            
            r.scales.push_back(scale);
            
            // Scale each individual ladder to the target
            r.rnd[i.first] = std::shared_ptr<RandomSelection>(new RandomSelection(1.0 - scale));
            
            o.logInfo("Ladder calibration (standard): " + i.first);
            o.logInfo("Ladder calibration (scale): " + std::to_string(scale));
        }
        
        // Calibrate ladder reads. Anything else will be written out.
        SelectionCalibrator::createFQ(f1, f2, o1, o2)->calibrate(o.d2->reads, r.rnd, o);
        
        return r;
    }
    
    FileName f1, f2, o1, o2;
};

struct BAMLadderCalibrator : public LadderCalibrator
{
    BAMLadderCalibrator(const FileName &src, const FileName &dst) : src(src), dst(dst) {}
    
    Results calibrate(const Options &o) override
    {
        Results r;

        // Calculate normalisation for each standard
        LadderCalibrator::getNorms(o, r);
        
        // Standard for each synthetic read (required by Calibrator)
        std::map<ReadName, StandardID> m;
        
        auto tmp = o.r1.inters();

        ParserBAM::parse(src, [&](ParserBAM::Data &x, const ParserBAM::Info &)
        {
            if (x.cID == GDecoyChrQL)
            {
                std::vector<DInter *> v;
                               
                if (tmp.at(x.cID).overlap(x.l, &v))
                {
                    auto name = getBestOverlap(x.l, v)->name();
                    
                    if (o.merged)
                    {
                        x.lName(); m[x.name] = noLast(name, "_");
                    }
                    else
                    {
                        x.lName(); m[x.name] = name;
                    }
                }
            }
        });

        // Calibrate ladder reads. Anything else will be written out.
        r.cr = SelectionCalibrator::createBAM(src, dst)->calibrate(m, r.rnd, o);

        return r;
    }
    
    FileName src, dst;
};

std::shared_ptr<LadderCalibrator> LadderCalibrator::createBAM(const FileName &src, const FileName &dst)
{
    return std::shared_ptr<LadderCalibrator>(new BAMLadderCalibrator(src, dst));
}

std::shared_ptr<LadderCalibrator> LadderCalibrator::createFQ(const FileName &f1, const FileName &f2,
                                                             const FileName &o1, const FileName &o2)
{
    return std::shared_ptr<LadderCalibrator>(new FQLadderCalibrator(f1, f2, o1, o2));
}

struct BAMSelectionCalibrator : public SelectionCalibrator
{
    virtual ~BAMSelectionCalibrator() {}
    
    Result calibrate(const std::map<ReadName, SequinID> &m, const Selection &r, const WriterOptions &wo) override
    {
        Result rr; rr.o1 = o;
        BAMWriter w; w.open(o);
     
        ParserBAM::parse(f, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
        {
            if (i.p && !(i.p % 10000)) { wo.wait(S0(i.p)); }
            x.lName();
            
            if (m.count(x.name)) { assert(r.count(m.at(x.name))); }
            if (!m.count(x.name) || r.at(m.at(x.name))->select(x.name))
            {
                rr.nSel++;
                w.write(x);
            }
            else
            {
                rr.nSkip++;
            }
        });
        
        w.close();
        return rr;
    }
    
    Result calibrate(Probability p, const WriterOptions &wo) override
    {
        RandomSelection r(1.0 - p);

        BAMWriter w; w.open(o);
        Result rr; rr.o1 = o;
        
        std::set<ReadName> ns;
        
        ParserBAM::parse(f, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
        {
            if (i.p && !(i.p % 10000)) { wo.wait(S0(i.p)); }
            x.lName();
            
            // Trimmed name
            const auto rn = trimRName(x.name);
            
            // Calibrate both pair-ended as they use the same trimmed read name
            if (r.select(rn))
            {
                rr.n++;
                
                // Only increment for a paired (not each read)
                if (!ns.count(rn))
                {
                    rr.nSel++;
                }
                else
                {
                    rr.nSkip++;
                }
                
                ns.insert(rn);
                w.write(x);
            }
        });

        w.close();
        return rr;
    }

    FileName f, o;
};

struct FQSelectionCalibrator : public SelectionCalibrator
{
    virtual ~FQSelectionCalibrator() {}

    Result calibrate(Probability p, const WriterOptions &wo) override
    {
        RandomSelection r(1.0 - p);
        
        Result rr;
        rr.o1 = o1; rr.o2 = o2;

        auto w1 = FQWriter(o1);
        auto w2 = FQWriter(o2);
        auto i = 0;
        
        ParserFQ::parse(Reader(f1), Reader(f2), [&](const ParserFQ::Data &x)
        {
            if (i && !(i % 10000)) { wo.wait(S0(i)); } i++;

            auto n1 = trimRName(x.name1);
            auto n2 = trimRName(x.name2);
            assert(n1 == n2);
            
            if (r.select(n1))
            {
                rr.nSel++;
                w1.write(n1, x.seq1, x.qual1);
                w2.write(n2, x.seq2, x.qual2);
            }
        });
        
        w1.close(); w2.close();
        return rr;
    }
    
    Result calibrate(const std::map<ReadName, SequinID> &m, const Selection &r, const WriterOptions &wo) override
    {
        FQWriter w1(o1);
        FQWriter w2(o2);

        ParserFQ::parse(Reader(f1), Reader(f2), [&](const ParserFQ::Data &x)
        {
            const auto n1 = first(replace(x.name1, "/1", ""), " ");
            const auto n2 = first(replace(x.name2, "/2", ""), " ");
            assert(n1 == n2);
            
            if (!m.count(n1) || r.at(m.at(n1))->select(n1))
            {
                w1.write(x.name1, x.seq1, x.qual1);
                w2.write(x.name2, x.seq2, x.qual2);
            }
        });

        w1.close();
        w2.close();

        Result rr; rr.o1 = o1; rr.o2 = o2;
        return rr;
    }

    FileName f1, f2, o1, o2;
};

std::shared_ptr<SelectionCalibrator> SelectionCalibrator::createBAM(const FileName &f, const FileName &o)
{
    auto x = std::shared_ptr<BAMSelectionCalibrator>(new BAMSelectionCalibrator());
    x->f = f; x->o = o;
    return x;
}

std::shared_ptr<SelectionCalibrator> SelectionCalibrator::createFQ(const FileName &f1, const FileName &f2,
                                                 const FileName &o1, const FileName &o2)
{
    auto x = std::shared_ptr<FQSelectionCalibrator>(new FQSelectionCalibrator());
    x->f1 = f1; x->f2 = f2; x->o1 = o1; x->o2 = o2;
    return x;
}
