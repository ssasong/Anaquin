#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include "data/data.hpp"
#include "data/library.hpp"
#include "data/standard.hpp"
#include "writers/writer.hpp"

namespace Anaquin
{
    struct LimitStats
    {
        // Absolute detection limit
        Limit limit;
    };

    /*
     * Common results between FASTQ and BAM implementation
     */

    struct CommonResults
    {
        virtual const Library &lib() const = 0;
        
        virtual Count binN(Bin) const = 0;

        virtual Proportion binP(Bin x) const
        {
            return (Proportion) binN(x) / total();
        }
        
        // Partition for each bin
        virtual ParitionCount pc() const = 0;

        // Read for each sequin
        virtual CustomMap<SequinID, Read> rn() const = 0;

        // Counts for a given sequin
        virtual Count count(const SequinID &) const = 0;

        virtual Count total() const
        {
            Count c = 0;
            forBin([&](Bin x) { c += binN(x); });
            return c;
        }

        virtual Proportion dil() const { return 1.0 - binP(ES); }
        
        /*
         * Only defined if calibrating
         */
        
        // Scaling factor
        virtual Proportion calibF(Calibration) const = 0;
        
        // Number of reads after calibration
        virtual Count calibA(Calibration) const = 0;
        
        // Dilution after calibration
        virtual Proportion calibD(Calibration) const = 0;
    };

    class WriterOptions
    {
        public:
            WriterOptions() : showGen(true),
                              showWarn(true),
                              writer(std::make_shared<MockWriter>()),
                              logger(std::make_shared<MockWriter>()),
                              output(std::make_shared<MockWriter>()) {}
        
            WriterOptions(const WriterOptions &x)
            {
                this->work = x.work;
                this->writer = x.writer;
                this->logger = x.logger;
                this->output = x.output;
            }
        
            // Working directory
            Path work;

            std::shared_ptr<Writer<>> writer, logger, output;

            inline void wait(const std::string &s) const
            {
                log("[WAIT]: " + s);
                out("[WAIT]: " + s);
            }

            inline void logWait(const std::string &s) const
            {
                log("[WAIT]: " + s);
            }
              
            inline void warn(const std::string &s) const
            {
                if (showWarn)
                {
                    log("[WARN]: " + s);
                    out("[WARN]: " + s);
                }
            }

            inline void logWarn(const std::string &s) const
            {
                if (showWarn)
                {
                    log("[WARN]: " + s);
                }
            }

            inline void analyze(const std::string &s) const
            {
                info("Analyzing: " + s);
            }
        
            inline void generate(const FileName &f) const
            {
                if (showGen) { info("Generating " + f); }
            }
        
            inline void info(const std::string &s) const
            {
                log("[INFO]: " + s);
                out("[INFO]: " + s);
            }
        
            inline void logInfo(const std::string &s) const
            {
                log("[INFO]: " + s);
            }

            inline void error(const std::string &s) const
            {
                info("[ERROR]: " + s);
            }
        
            bool showWarn, showGen;
        
        private:
        
            inline void out(const std::string &s) const { if (output) { output->write(s); } }
            inline void log(const std::string &s) const { if (logger) { logger->write(s); } }
    };

    struct AnalyzerOptions : public WriterOptions
    {
        AnalyzerOptions() {}
        AnalyzerOptions(const AnalyzerOptions &x) : WriterOptions(x), version(x.version), name(x.name), cmd(x.cmd), debug(false) {}
        
        bool debug;
        
        // Eg: "split"_summary.txt
        std::string name;
        
        // Full comamnd
        std::string cmd;
        
        std::string version;
        
        // BAM mode?
        bool bam;
    };
    
    inline bool shouldTrim(const Locus &l, Base start, Base end, Base trim)
    {
        const auto lTrim = std::abs(l.start - start) <= trim;
        const auto rTrim = std::abs(l.end - end) <= trim;
        return lTrim || rTrim;
    }
}

#endif
