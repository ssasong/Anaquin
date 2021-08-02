#include <future>
#include <algorithm>
#include "stats/stats.hpp"
#include "tools/calibrator.hpp"
#include "parsers/parser_csv.hpp"
#include "sequins/genomics/g_norm.hpp"
#include "stats/ss/regression/linear.hpp"

using namespace Anaquin;

typedef GNorm::Stats Stats;
typedef GNorm::Options Options;

extern std::map<std::string, std::map<int, Anaquin::KMInfo>> __KMInfo__;

static void split(const FileName &f1)
{
    GSplit::analyze(f1, "", GSplit::Options());
}

void GNorm::writeNorm(const FileName &src, const FileName &dst, GNorm::Stats &stats, const Options &o)
{
    std::set<Label> labs = { "SAMP_A", "SAMP_B", "S1", "S2" };
    if (stats.S.size() > 2) { labs.insert("SAMP_C"); }
    if (stats.S.size() > 3) { labs.insert("SAMP_D"); }
    if (stats.S.size() > 4) { labs.insert("SAMP_E"); }
    
    const auto tmp = tmpFile();
    RFilterC(src, tmp, labs, true);

    // Number of samples
    const auto n = stats.S.size();
    
    const auto dt = SS::MatrixTools::readMatrix(tmp, RCount(tmp, "", ""), n, true);
    const auto dt_norm = TMM(dt);

    for (auto i = 0; i < n; i++)
    {
        const auto S  = SS::MatrixTools::toSTDVect(dt.col(i));
        const auto SN = SS::MatrixTools::toSTDVect(dt_norm.col(i));
        stats.P.push_back(SS::linearModel(false, SN, S).coeffs[0].est);
        o.info("Scaling Factor: " + std::to_string(stats.P[i]));
    }
    
    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%";

    std::stringstream ss;
    ss << (boost::format(f) % "NAME"
                            % "SEQUENCE"
                            % "SAMP_A"
                            % "SAMP_B"
                            % "SAMP_C"
                            % "SAMP_D"
                            % "SAMP_E"
                            % "NORM_A"
                            % "NORM_B"
                            % "NORM_C"
                            % "NORM_D"
                            % "NORM_E").str() << std::endl;

    ParserCSV::parse(Reader(src), [&](const ParserCSV::Data &x, Progress i)
    {
        if (!i || x.size() == 1) { return; }

        const auto SA = stod(x[2]);
        const auto SB = stod(x[3]);
        const auto SC = stats.S.size() > 2 ? stod(x[4]) : NAN;
        const auto SD = stats.S.size() > 3 ? stod(x[5]) : NAN;
        const auto SE = stats.S.size() > 4 ? stod(x[6]) : NAN;

        const auto NA = (stats.P.size() > 0 ? stats.P[0] * SA : NAN);
        const auto NB = (stats.P.size() > 1 ? stats.P[1] * SA : NAN);
        const auto NC = (stats.P.size() > 2 ? stats.P[2] * SA : NAN);
        const auto ND = (stats.P.size() > 3 ? stats.P[3] * SA : NAN);
        const auto NE = (stats.P.size() > 4 ? stats.P[4] * SA : NAN);

        ss << (boost::format(f) % x[0]
                                % x[1]
                                % SA
                                % SB
                                % SC
                                % SD
                                % SE
                                % NA
                                % NB
                                % NC
                                % ND
                                % NE) << std::endl;
    }, "\t");
    
    o.generate(dst);
    o.writer->open(dst);
    o.writer->write(ss.str());
    o.writer->close();
}

static void writeKM(const FileName &file, const std::vector<const KStats *> &x, const Options &o)
{
    std::map<Kmer, std::map<std::size_t, Count>> m;
    
    for (auto &i : __KMInfo__)
    {
        if (isSubstr(i.first, "LD_"))
        {
            for (auto &j : i.second)
            {
                for (auto k = 0; k < x.size(); k++)
                {
                    m[j.second.kmer][k] = 0;
                }
            }
        }
    }

    auto add = [&](size_t id, const SequinKM &x)
    {
        for (const auto &i : x)
        {
            if (isSubstr(i.first, "LD_"))
            {
                for (const auto &j : i.second)
                {
                    m[j.first][id] += j.second;
                }
            }
        }
    };
    
    for (auto i = 0; i < x.size(); i++)
    {
        add(i, x.at(i)->uniqs);
        add(i, x.at(i)->shared);
    }
    
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "NAME"
                                           % "SEQUENCE"
                                           % "SAMP_A"
                                           % "SAMP_B"
                                           % "SAMP_C"
                                           % "SAMP_D"
                                           % "SAMP_E").str());
    
    for (const auto &i : __KMInfo__)
    {
        if (isSubstr(i.first, "LD_"))
        {
            for (const auto &j : i.second)
            {
                assert(m.count(j.second.kmer));
                o.writer->write((boost::format(format) % i.first
                                                       % j.second.kmer
                                                       % m[j.second.kmer][0]
                                                       % m[j.second.kmer][1]
                                                       % (m.size() > 2 ? S0(m[j.second.kmer][2]) : MISSING)
                                                       % (m.size() > 3 ? S0(m[j.second.kmer][3]) : MISSING)
                                                       % (m.size() > 4 ? S0(m[j.second.kmer][4]) : MISSING)).str());
            }
        }
    }
    
    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const std::vector<FileName> &f1,
                         const std::vector<FileName> &f2,
                         const Stats &stats,
                         const Options &o)
 {
     o.generate(file);
     o.writer->open(file);
 
      const auto f = "SEQUIN REPORT:                 %1%\n\n"
                     "REFERENCE FILES\n"
                     "Reference index:               %2%\n\n"
                     "LIBRARY INFORMATION\n"
                     "Version:                       %3%\n"
                     "Instrument ID:                 %4%\n"
                     "Run number:                    %5%\n"
                     "Flowcell ID:                   %6%\n"
                     "Lane:                          %7%\n\n"
                     "USER-SUPPLIED FILES\n"
                     "Input files (first):           %8%\n"
                     "Input files (second):          %9%\n"
                     "Input files (third):           %10%\n"
                     "Input files (fourth):          %11%\n"
                     "Input files (fifth):           %12%\n\n"
                     "ANAQUIN PARAMETERS\n"
                     "K-mer length:                  %13%\n"
                     "Threshold:                     %14%\n\n"
                     "NORMALIZATION SUMMARY\n"
                     "Scaling factor (sample A):     %15%\n"
                     "Scaling factor (sample B):     %16%\n"
                     "Scaling factor (sample C):     %17%\n"
                     "Scaling factor (sample D):     %18%\n"
                     "Scaling factor (sample E):     %19%\n"
                     "K-mers table:                  %20%\n";

     // Statistics before calibration
     const auto p1 = (CommonResults *) &stats.S;
     
     // Libray format
     const auto fo = p1->lib().format();

     //const auto l1 = stats.l1.linear(false);
     //const auto l2 = stats.l2.linear(false);

     o.writer->write((boost::format(f) % date()
                                       % o.index
                                       % SVersion(Standard::instance().gen, stats.S[0].S1.K)
                                       % p1->lib().inst(fo) // 4
                                       % p1->lib().run(fo)  // 5
                                       % p1->lib().flow(fo) // 6
                                       % p1->lib().lane(fo) // 7
                                       % (f1.size() > 0 ? f1[0] + " and " + f1[0] : "") // 8
                                       % (f1.size() > 1 ? f1[1] + " and " + f1[1] : "") // 9
                                       % (f1.size() > 2 ? f1[2] + " and " + f1[2] : "") // 10
                                       % (f1.size() > 3 ? f1[3] + " and " + f1[3] : "") // 11
                                       % (f1.size() > 4 ? f1[4] + " and " + f1[4] : "") // 12
                                       % o.k                  // 13
                                       % o.rule               // 14
                                       % (stats.S.size() > 0 ? S2(stats.P[0]) : MISSING) // 15
                                       % (stats.S.size() > 1 ? S2(stats.P[1]) : MISSING) // 16
                                       % (stats.S.size() > 2 ? S2(stats.P[2]) : MISSING) // 17
                                       % (stats.S.size() > 3 ? S2(stats.P[3]) : MISSING) // 18
                                       % (stats.S.size() > 4 ? S2(stats.P[4]) : MISSING) // 19
                                       % (o.work + "/norm_kmers.tsv")                    // 20
                      ).str());

     o.writer->close();
}

Stats GNorm::analyze(const std::vector<FileName> &f1, const std::vector<FileName> &f2, const Options &o)
{
    const auto o2 = cloneO(o);
    
    std::vector<std::future<GSplit::Stats>> p;

    for (auto i = 0; i < f1.size(); i++)
    {
        p.push_back(std::async(std::launch::deferred, [i, f1, f2, o2]() {
            return GSplit::analyze(f1[i], f2[i], o2);
        }));
    }
    
    Stats stats;
    
    for (auto i = 0; i < p.size(); i++)
    {
        stats.S.push_back(p[i].get());
    }
                    
    return stats;
}

static void writeRNorm(const FileName &file,
                       const FileName &src,
                       const Label &x,
                       const Label &y,
                       const Label &title,
                       const Label &xl,
                       const Label &yl,
                       const SOptions &o)
{
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(PlotNorm()) % date()
                                               % o.cmd
                                               % o.work
                                               % src
                                               % x
                                               % y
                                               % title
                                               % xl
                                               % yl).str());
    o.writer->close();
}

void GNorm::report(const std::vector<FileName> &f1, const std::vector<FileName> &f2, const Options &o)
{
    auto stats = analyze(f1, f2, o);
    
    std::vector<const KStats *> ks;
    for (auto &stat : stats.S) { ks.push_back(&stat.S1.K); }
    
    // Generating norm_kmers.tsv
    writeKM("old_norm_kmers.tsv", ks, o);
    
    // Determining scaling factors (from already written "norm_kmers.tsv"). Also determining "P"
    GNorm::writeNorm(o.work + "/old_norm_kmers.tsv", "norm_kmers.tsv", stats, o);
    
    // Calibrator for each sample
    std::vector<std::shared_ptr<SelectionCalibrator>> calib;
    
    for (auto i = 0; i < stats.S.size(); i++)
    {
        const auto o1 = o.work + "/" + "norm_library" + std::to_string(i+1) + "_1.fq.gz";
        const auto o2 = o.work + "/" + "norm_library" + std::to_string(i+1) + "_2.fq.gz";

        const auto scale = stats.scale(i);
        
        // Don't calibrate if there's no ladder reads or TMM fail...
        const auto shouldCalib = !std::isnan(scale) && !std::isinf(scale) && scale >= 0.0 && scale <= 1.0;

        SelectionCalibrator::createFQ(f1[i], f2[i], o1, o2)->calibrate(shouldCalib ? scale : 1.0, o);
    }
    
    // Generating norm_summary.txt
    writeSummary("norm_summary.txt", f1, f2, stats, o);

    if (o.debug)
    {
        // Generating norm_before.R
        writeRNorm("report_files/norm_before.R", "norm_kmers.tsv", "data$S1", "data$S2", "Ladder Sequins (before normalisation)", "Sample 1", "Sample 2", o);

        // Generating norm_after.R
        writeRNorm("report_files/norm_after.R", "norm_kmers.tsv", "data$NORM_S1", "data$NORM_S2", "Ladder Sequins (after normalisation)", "Sample 1 (Normalized)", "Sample 2 (Normalized)", o);
    }
}
