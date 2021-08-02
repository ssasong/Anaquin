#include "sequins/sequins.hpp"
#include "tools/calibrator.hpp"
#include "writers/bam_writer.hpp"
#include "parsers/parser_bam.hpp"
#include "sequins/decoy_analyzer.hpp"

using namespace Anaquin;

typedef Anaquin::DecoyAnalyzer::D1 D1;
typedef Anaquin::DecoyAnalyzer::D2 D2;
typedef DecoyAnalyzer::Options Options;
typedef DecoyAnalyzer::Results Results;
typedef DecoyAnalyzer::Results::Metrics Metrics;

struct AnalyzeBAMResults
{
    std::shared_ptr<BAMWriter> wM;
    
    /*
     * Metrics only defined if calibration
     */

    Proportion p; // Scaling factor
    
    Count tar;
    Count bSam; // Sample reads before calibration
    Count aSam; // Sample reads after calibration
    Count bSeq; // Sequin reads before calibration
    Count aSeq; // Sequin reads after calibration
};

Count DecoyAnalyzer::Results::binN(Bin x) const
{
    switch (x)
    {
        case Bin::ES: { return samp.total();  }
        default:      { return decoy.total(); }
    }
}

CustomMap<SequinID, Read> DecoyAnalyzer::Results::rn() const
{
    CustomMap<SequinID, Read> x;
    
    for (const auto &i : decoy.r2)
    {
        for (const auto &j : i.second.data())
        {
            x[j.first] = j.second.stats().n;
        }
    }
    
    return x;
}

Count DecoyAnalyzer::Results::count(const SequinID &x) const
{
    assert(decoy.r2.find(x));
    return decoy.r2.find(x)->stats().n;
}

static AnalyzeBAMResults analyzeBAM(const FileName &f1,
                                    const FileName &f2,
                                    Metrics &samp,
                                    Metrics &decoy,
                                    Library &lib,
                                    const CustomMap<ChrID, Sequence> &chrs,
                                    const Options &o,
                                    D1 d1,
                                    D2 d2)
{
    assert(!chrs.empty());
    const auto isChrQ = f2.empty();

    assert(!decoy.r2.empty());
    
    o.analyze(f1);
    if (!f2.empty()) { o.analyze(f2); }
    std::shared_ptr<BAMWriter> wS, wD, wT, wM;
    
    if (!o.writeS.empty()) { wS = std::shared_ptr<BAMWriter>(new BAMWriter()); wS->open(o.writeS); }
    if (!o.writeD.empty()) { wD = std::shared_ptr<BAMWriter>(new BAMWriter()); wD->open(o.writeD); }
    if (!o.writeT.empty()) { wT = std::shared_ptr<BAMWriter>(new BAMWriter()); wT->open(o.writeT); }
    if (!o.writeM.empty()) { wM = std::shared_ptr<BAMWriter>(new BAMWriter()); wM->open(o.writeM); }

    // Trimmed read names
    std::set<ReadName> trimmed;
    
    AnalyzeBAMResults r;

    auto analyze = [&](const FileName &file, bool isSam, bool isSeq)
    {
        ParserBAM::parse(file, [&](ParserBAM::Data &x, const ParserBAM::Info &i)
        {
            if (i.p && !(i.p % 10000)) { o.wait(S2(i.p)); }
            
            if (!i.p)
            {
                x.lName(); x.lSeq(); lib.addInfo(x.name, x.seq);
            }
            
            // Only primary alignments are considered
            if (!x.isPrimary) { return; }
            
            std::vector<DInter *> v1, v2, v3;

            // Overlapped with human regions (mirrored with sequin regions)?
            DInter *h1 = isSam && samp.r1.count(x.cID) ? samp.r1.at(x.cID).overlap(x.l) : nullptr;

            // Overlapped with trimmed human regions?
            DInter *h2 = isSam && samp.r2.count(x.cID) ? samp.r2.at(x.cID).overlap(x.l) : nullptr;

            // It's possible a read overlaps two regions
            v1.clear(); v2.clear();
            
            // Overlapped with non-trimmed decoy regions?
            DInter *r1 = isSeq && !h1 && decoy.r1.count(x.cID) ? decoy.r1.at(x.cID).overlap(x.l, &v1) : nullptr;
            
            // Make sure we get the best overlap
            r1 = r1 ? getBestOverlap(x.l, v1) : nullptr;
            
            // Overlapped with trimmed decoy regions?
            DInter *r2 = isSeq && !h1 && decoy.r2.count(x.cID) ? decoy.r2.at(x.cID).overlap(x.l, &v2) : nullptr;
            
            // Make sure we get the best overlap
            r2 = r2 ? getBestOverlap(x.l, v2) : nullptr;
            
            auto f = [&](Results::Metrics &rm, const CustomMap<ChrID, Sequence> &ind, bool shouldA)
            {
                const auto isR1Empty = rm.r1.empty();
                
                if (!isR1Empty && (!rm.r1.count(x.cID) || !rm.r1.at(x.cID).overlap(x.l)))
                {
                    rm.out++;
                }
                else
                {
                    rm.in++;
                    
                    if (isR1Empty)
                    {
                        return;
                    }

                    const auto m1 = rm.r1.at(x.cID).overlap(x.l);
                    const auto m2 = rm.r2.count(x.cID) ? rm.r2.at(x.cID).overlap(x.l) : nullptr; // Trimmed
                    assert(m1);
                    
                    auto name = m1->name();
                    
                    if (countT(name, "_") != 2)
                    {
                        name = noFirst(name, "_");
                    }
                    
                    // How to search the index?
                    const auto key = isChrQ ? x.cID : noPID(m1->name() + "_R");
                    
                    if (!ind.find(key) || name.empty())
                    {
                        return;
                    }
                    
                    // Don't run the time-consuming cigar unless it's required
                    if (!shouldA)
                    {
                        return;
                    }
                    
                    // Only the requested chromosomes will be reported
                    else if (!o.errors.count(x.cID))
                    {
                        return;
                    }
                    
                    x.lCigar(); x.lSeq();

                    const auto seq = *(ind.find(key));

                    auto s1 = (int) x.l.start; // Temporary starting
                    auto s2 = (int) x.l.end;   // Temporary ending
                    
                    /*
                     * If hg19/hg38, convert to sequin relative just like chrQ
                     */
                    
                    if (!isChrQ)
                    {
                        s1 -= (m1->l().start - 1);
                        s2 -= (m1->l().start - 1);
                    }
                    
                    // Read extend outside sequin region? Don't use it.
                    if (s1 <= 0 || s2 >= seq.size())
                    {
                        return;
                    }
                    
                    auto tmp = x;
                    tmp.l.start = s1; tmp.l.end = s2;
                    assert(tmp.l.start); // 1-based, so zero is impossible
                    
                    forCigar(tmp, seq, [&](Base i, Base j, char r, char q) // Matching
                    {
                        rm.sd["All"].match[name][j]++;
                    }, [&](Base i, Base j, char r, char q) // SNP
                    {
                        auto cb = SNPType::AC;
                        
                        if (r == 'A' && q == 'C') { cb = SNPType::AC; }
                        if (r == 'A' && q == 'T') { cb = SNPType::AT; }
                        if (r == 'A' && q == 'G') { cb = SNPType::AG; }
                        if (r == 'C' && q == 'A') { cb = SNPType::CA; }
                        if (r == 'C' && q == 'T') { cb = SNPType::CT; }
                        if (r == 'C' && q == 'G') { cb = SNPType::CG; }
                        if (r == 'T' && q == 'A') { cb = SNPType::TA; }
                        if (r == 'T' && q == 'G') { cb = SNPType::TG; }
                        if (r == 'T' && q == 'C') { cb = SNPType::TC; }
                        if (r == 'G' && q == 'A') { cb = SNPType::GA; }
                        if (r == 'G' && q == 'C') { cb = SNPType::GC; }
                        if (r == 'G' && q == 'T') { cb = SNPType::GT; }

                        rm.sd["All"].snps[name][j][cb]++;
                    }, [&](Base i, Base j, Base l) // Insertion
                    {
                        rm.sd["All"].ins[name][j][l]++;
                    }, [&](Base i, Base j, Base l) // Deletion
                    {
                        rm.sd["All"].dls[name][j][l]++;
                    }, [&](Base, Base, Base) // Skip and clips
                    {
                        // Ignored
                    });
                }
            };
            
            const auto mustSam = isSam && !isSeq;
            const auto mustSeq = isSeq && !isSam;
            const auto both    = isSam && isSam;
            
            // Decoy chromosome?
            const auto isDecoyC = decoy.r2.count(x.cID);
            
            // Only trimmed if it's decoy and satisfy conditions
            auto isTrimmed = false;
            
            // Sample read?
            if (mustSam || (both && !isDecoyC))
            {
                if (wS) { wS->write(x); }
                if (wM) { wM->write(x); }

                f(samp, chrs, false); // Update general statistics

                if (h1) { h1->map(x.l); }
                if (h2) { h2->map(x.l); } // Ignore trimming, just set h1 == h2
            }
            else
            {
                if (wD) { wD->write(x); } // Writing for decoy
                
                f(decoy, chrs, o.shouldError); // Update general statistics

                if (r1)
                {
                    if (o.shouldTrim && shouldTrim(x.l, r1->l().start, r1->l().end, o.trim))
                    {
                        x.lName();
                        assert(!x.name.empty());
                        trimmed.insert(x.name);
                        isTrimmed = true;
                    }
                    else if (r2)
                    {
                        // Update coverage on trimmed regions
                        r2->map(x.l);
                    }
                    
                    r1->map(x.l);
                }
                
                // Overlapped with feature regions?
                if (decoy.r5.count(x.cID) && decoy.r5.at(x.cID).overlap(x.l, &v3))
                {
                    for (auto &i : v3) { i->map(x.l); }
                }
            }
            
            d1(x, r1, r2, isTrimmed);
        });
    };
    
    if (f2.empty())
    {
        analyze(f1, true, true);
    }
    else
    {
        if (!f1.empty()) { analyze(f1, true, false); }
        if (!f2.empty()) { analyze(f2, false, true); }
    }

    d2(); // Completed
    if (wS) { wS->close(); }
    if (wD) { wD->close(); }

    if (wD && wT) // wD is not longer valid but it's still not NULL
    {
        ParserBAM::parse(o.writeD, [&](ParserBAM::Data &x, const ParserBAM::Info &)
        {
            x.lName();
            
            if (!trimmed.count(x.name))
            {
                wT->write(x);
            }
        });
    }
    
    o.info(f1 + " completed");
    if (!f2.empty()) { o.info(f2 + " completed"); }
    
    if (wT) { wT->close(); }
    
    // Return the writer without closing it, so calibrated reads can be added
    r.wM = wM;
    
    // Calibrate sequin reads?
    if (o.seqC != NO_CALIBRATION)
    {
        assert(!o.inputC.empty());
        assert(!o.writeC.empty());
        
        r.bSam = r.aSam = samp.total();
        r.bSeq = 0; ParserBAM::parse(o.inputC, [&](ParserBAM::Data &x, const ParserBAM::Info &) {
            r.bSeq++;
        }); // r.bSeq could be just decoy.total() but doesn't have to be

        auto tmp = o.seqC;
        
        // Convert to target paired reads (very important)
        //tmp = 2.0 * tmp;

        if (tmp > 1.0) // Absolute?
        {
            o.logInfo("Absolute calibration");
            tmp = int(tmp);
            
            // Number of target sequin reads after calibration
            r.tar = r.bSeq < tmp ? r.bSeq : tmp;

            // Scaling factor for calibration (0.5 because of paired-ends)
            r.p = r.bSeq ? (Proportion) r.tar / r.bSeq : 0.0;
        }
        else // Percentage?
        {
            o.logInfo("Percentage calibration");

            // Number of target sequin reads after calibration
            r.tar = (tmp / (1.0 - tmp)) * r.bSam;
            
            // Make sure the target doesn't goto zero if sample reads is non-zero
            if (r.bSam && !r.tar)
            {
                o.logInfo("Set target to sample");
                r.tar = r.bSam;
            }
            
            // Scaling factor for calibration
            r.p = (r.bSam == 0) ? 1.0 : (r.bSeq == 0) ? 0.0 : (r.tar >= r.bSeq ? 1.0 : ((float) r.tar) / r.bSeq);
        }

        assert(!std::isnan(r.p));
        assert(!std::isnan(r.tar));
        
        o.logInfo("Input: "   + o.inputC);
        o.logInfo("Output: "  + o.writeC);
        o.logInfo("Sample: "  + S0(r.bSam));
        o.logInfo("Sequin: "  + S0(r.bSeq));
        o.logInfo("Target: "  + S0(r.tar));
        o.logInfo("Scaling: " + std::to_string(r.p));
        
        // Calibrate and return the number of reads after calibration
        r.aSeq = SelectionCalibrator::createBAM(o.inputC, o.writeC)->calibrate(r.p, o).nSel;
    }
    
    return r;
}

Results DecoyAnalyzer::analyze(const FileName &f1, const FileName &f2, const Options &o, D1 d1, D2 d2)
{
    assert(exists(o.index));
    assert(!o.errors.empty());
    assert(!o.seqs.empty());
    
    assert(!o.r1.inters().empty());
    assert(!o.r2.inters().empty());

    Results r;

    r.samp.r1  = o.h1.inters();   // No trimming
    r.samp.r2  = o.h2.inters();   // With trimming
    r.decoy.r1 = o.r1.inters();   // No trimming
    r.decoy.r2 = o.r2.inters();   // With trimming
    r.decoy.r5 = o.attr.inters(); // Attributes

    // Regions without trimming must always be there
    assert(!r.decoy.r1.empty());
    
    const auto tmp = analyzeBAM(f1, f2, r.samp, r.decoy, r.bamLib, o.seqs, o, d1, d2);

    r._p    = tmp.p;
    r.wM    = tmp.wM;
    r._tar  = tmp.tar;
    r._bSam = tmp.bSam;
    r._aSam = tmp.aSam;
    r._bSeq = tmp.bSeq;
    r._aSeq = tmp.aSeq;

    return r;
}

void DecoyAnalyzer::writeE(const ErrorOptions &o)
{
    if (o.tsvE.empty())
    {
        return;
    }
    
    auto m2 = std::map<SNPType, Label>
    {
        {SNPType::AC, "AC"}, {SNPType::AT, "AT"}, {SNPType::AG, "AG"},
        {SNPType::CA, "CA"}, {SNPType::CT, "CT"}, {SNPType::CG, "CG"},
        {SNPType::TA, "TA"}, {SNPType::TG, "TG"}, {SNPType::TC, "TC"},
        {SNPType::GA, "GA"}, {SNPType::GC, "GC"}, {SNPType::GT, "GT"},
        {SNPType::RF, "Ref"}
    };
    
    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%";
    
    o.generate(o.tsvE);
    o.writer->open(o.tsvE);
    o.writer->write("CALIBRATED\tCHROM\tNAME\tLABEL\tSTRAT\tDATA_1\tDATA_2\tDATA_3");
    
    auto r2 = o.r2.ginters();
    
    auto shouldReport = [r2, o](Base x)
    {
        return !o.isDecoy || r2.at(o.chr).overlap(Locus(x, x));
    };
    
    auto run = [&](const Results::Metrics &m, const std::string &calib)
    {
        for (auto &i : m.sd)
        {
            for (auto &j : i.second.match)
            {
                for (auto &k : j.second)
                {
                    o.writer->write((boost::format(f) % calib
                                                      % o.chr
                                                      % j.first
                                                      % "Match"
                                                      % i.first       // Stratification
                                                      % (k.first + 1) // 1-based Base
                                                      % k.second      // Count
                                                      % MISSING).str());
                }
            }
            
            for (auto &j : i.second.snps)
            {
                for (auto &k : j.second)
                {
                    if (!shouldReport(k.first + 1))
                    {
                        continue;
                    }
                    
                    // For each SNP type ...
                    for (auto &l : k.second)
                    {
                        auto v = o.v1 ? o.v1->data.findVar(o.chr, Locus(k.first+1, k.first+1)) : nullptr;
                        
                        // Skip sequin SNP variant?
                        if (v && v->type() == Variation::SNP && v->snpType() == l.first)
                        {
                            continue;
                        }
                        
                        const auto &name = j.first;
                        const auto r = !o.isDecoy ? o.r1.find(noPID(name), false) : nullptr;
                        assert(o.isDecoy || r);
                        
                        // 1-based base
                        const auto ll = o.isDecoy ? (k.first + 1) : (r->l.start + k.first);
                        
                        o.writer->write((boost::format(f) % calib
                                                          % (o.isDecoy ? o.chr : r->cID)
                                                          % name
                                                          % "SNP"
                                                          % i.first        // Stratification
                                                          % ll
                                                          % m2.at(l.first) // Type
                                                          % l.second       // How many?
                                         ).str());
                    }
                }
            }
            
            for (auto &j : i.second.ins)
            {
                for (auto &k : j.second)
                {
                    if (!shouldReport(k.first + 1))
                    {
                        continue;
                    }
                    
                    for (auto &l : k.second)
                    {
                        if (!o.v1 || !o.v1->data.findVar(o.chr, Locus(k.first+1, k.first+1)))
                        {
                            const auto &name = j.first;
                            const auto r = !o.isDecoy ? o.r1.find(noPID(name), false) : nullptr;
                            assert(o.isDecoy || r);

                            // 1-based base
                            const auto ll = o.isDecoy ? k.first : (r->l.start + k.first);

                            o.writer->write((boost::format(f) % calib
                                                              % (o.isDecoy ? o.chr : r->cID)
                                                              % name
                                                              % "Insertion"
                                                              % i.first // Stratification
                                                              % ll
                                                              % l.first
                                                              % l.second).str());
                            
                        }
                    }
                }
            }
            
            for (auto &j : i.second.dls)
            {
                for (auto &k : j.second)
                {
                    if (!shouldReport(k.first + 1))
                    {
                        continue;
                    }
                    
                    for (auto &l : k.second)
                    {
                        if (!o.v1 || !o.v1->data.findVar(o.chr, Locus(k.first+1, k.first+1)))
                        {
                            const auto &name = j.first;
                            const auto r = !o.isDecoy ? o.r1.find(noPID(name), false) : nullptr;
                            assert(o.isDecoy || r);

                            // 1-based base
                            const auto ll = o.isDecoy ? k.first : (r->l.start + k.first);

                            o.writer->write((boost::format(f) % calib
                                                              % (o.isDecoy ? o.chr : r->cID)
                                                              % name
                                                              % "Deletion"
                                                              % i.first // Stratification
                                                              % ll
                                                              % l.first
                                                              % l.second).str());
                        }
                    }
                }
            }
        }
    };
    
    // Combing error profiles
    for (auto &i : o.data) { run(*i.first, i.second ? "true" : "false"); }
    
    o.writer->close();
    
    if (o.debug)
    {
        o.writer->open(replace(o.tsvE, ".tsv", ".bed"));
        
        auto run2 = [&](const Results::Metrics &m, const std::string &calib)
        {
            for (auto &i : m.sd)
            {
                if (i.first != "All")
                {
                    continue;
                }
                
                for (auto &j : i.second.snps) // For each SNP
                {
                    for (auto &k : j.second) // For each reported position
                    {
                        if (!shouldReport(k.first + 1))
                        {
                            continue;
                        }
                        
                        for (auto &l : k.second) // For each SNP type
                        {
                            auto v = o.v1 ? o.v1->data.findVar(o.chr, Locus(k.first+1, k.first+1)) : nullptr;
                            
                            if (v && v->type() == Variation::SNP && v->snpType() == l.first)
                            {
                                continue;
                            }
                            
                            const auto m = k.first + 1; // 1-based Base
                            const auto s = S0(m) + "_" + m2.at(l.first) + "_" + S0(l.second);
                            
                            const auto f = "%1%\t%2%\t%3%\t%4%";
                            o.writer->write((boost::format(f) % o.chr
                                                              % m
                                                              % m
                                                              % s).str());
                        }
                    }
                }
            }
        };
        
        // BED for SNP before calibration
        for (auto &i : o.data) { run2(*i.first, i.second ? "true" : "false"); }
        
        o.writer->close();
    }
}
