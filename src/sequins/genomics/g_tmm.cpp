#include "stats/stats.hpp"
#include "sequins/genomics/g_tmm.hpp"
#include "parsers/parser_csv.hpp"

using namespace Anaquin;

typedef GTMM::Stats Stats;
typedef GTMM::Options Options;

static void writeKM(const FileName &file, const Stats &stats, const Options &o)
{
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(format) % "NAME"
                                           % "S1"
                                           % "S2"
                                           % "S3"
                                           % "S4"
                                           % "S5"
                                           % "NORM_S1"
                                           % "NORM_S2"
                                           % "NORM_S3"
                                           % "NORM_S4"
                                           % "NORM_S5").str());
    
    for (auto i = 0; i < stats.seqs.size(); i++)
    {
        o.writer->write((boost::format(format) % stats.seqs[i]
                                               % stats.dt(i,0)
                                               % stats.dt(i,1)
                                               % (stats.dt.cols() > 2 ? S0(stats.dt(i,2)) : MISSING)
                                               % (stats.dt.cols() > 3 ? S0(stats.dt(i,3)) : MISSING)
                                               % (stats.dt.cols() > 4 ? S0(stats.dt(i,4)) : MISSING)
                                               % stats.dt_n(i,0)
                                               % stats.dt_n(i,1)
                                               % (stats.dt.cols() > 2 ? S0(stats.dt_n(i,2)) : MISSING)
                                               % (stats.dt.cols() > 3 ? S0(stats.dt_n(i,3)) : MISSING)
                                               % (stats.dt.cols() > 4 ? S0(stats.dt_n(i,4)) : MISSING)).str());
    }

    o.writer->close();
}

Stats GTMM::analyze(const FileName &file, const Options &o)
{
    const auto tmp = tmpFile();
    RFilterC(file, tmp, std::set<std::size_t> { 1, 2, 3, 4, 5, 6, 7, 8, 9 }, true);
    
    Stats stats;
    
    stats.dt   = SS::MatrixTools::readMatrix(tmp, RCount(tmp, "", ""), RCountCols(tmp), true);
    stats.dt_n = TMM(stats.dt);
    
    RFilterC(file, tmp, std::set<std::size_t> { 0 }, true);

    ParserCSV::parse(tmp, [&](const ParserCSV::Data &x, Progress i)
    {
        if (i)
        {
            assert(x.size() == 1);
            stats.seqs.push_back(x[0]);
        }
    });
    
    assert(stats.seqs.size() == stats.dt.rows());
    assert(stats.seqs.size() == stats.dt_n.rows());
    
    return stats;
}

void GTMM::report(const FileName &file, const Options &o)
{
    writeKM("tmm_kmers.tsv", analyze(file, o), o);
}
