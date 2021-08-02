#include <ctime>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <strings.h>
#include <execinfo.h>
#include <sys/stat.h>

#include "sequins/rna/rna.hpp"
#include "sequins/meta/meta.hpp"
#include "parsers/parser_csv.hpp"
#include "sequins/rna/r_split.hpp"
#include "sequins/meta/m_split.hpp"
#include "sequins/genomics/g_tmm.hpp"
#include "sequins/genomics/g_germ.hpp"
#include "sequins/genomics/g_norm.hpp"
#include "sequins/genomics/g_split.hpp"
#include "sequins/genomics/g_somatic.hpp"
#include "sequins/genomics/g_calibrate.hpp"
#include "sequins/genomics/g_broad_bam.hpp"
#include "sequins/genomics/g_broad_vcf.hpp"

#include "parsers/parser_csv.hpp"
#include "parsers/parser_vcf.hpp"

#include "tools/tools.hpp"
#include "data/resources.hpp"
#include "tools/detector.hpp"
#include "tools/bedtools.hpp"
#include "tools/checkPackage.h"
#include "tools/attributes.hpp"
#include "writers/file_writer.hpp"
#include "writers/terminal_writer.hpp"

#ifdef UNIT_TEST
#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#endif

#define DEFAULT_EDGE 550
#define DEFAULT_LADDER_POOL_CALIBRATION 0.1

using namespace Anaquin;

typedef int Option;

typedef std::string Value;

static Build __ha__;

static std::string version() { return "3.22.0"; }

/*
 * Options specified in the command line
 */

#define OPT_TEST    320
#define OPT_TOOL    321
#define OPT_PATH    325

#define OPT_R_BED    801
#define OPT_METHOD   802
#define OPT_R_VCF    803
#define OPT_COMBINE  804
#define OPT_MIXTURE  808
#define OPT_R_HUMAN  809
#define OPT_R_DECOY  810
#define OPT_R_REGS   811
#define OPT_RESOURCE 812
#define OPT_BUILD    813
#define OPT_SYNC     814
#define OPT_ATTTSV   815
#define OPT_FASTA    816
#define OPT_KMER     817
#define OPT_RULE     818
#define OPT_THREAD   819
#define OPT_EDGE     820
#define OPT_SKIP     821
#define OPT_MERGE    822
#define OPT_WINDOW   823
#define OPT_1        824
#define OPT_2        825
#define OPT_3        826
#define OPT_4        827
#define OPT_SAMPLE   829
#define OPT_SEQUIN   830
#define OPT_LCALIB   831
#define OPT_FASTQ    832
#define OPT_DEBUG    833
#define OPT_VCF      834
#define OPT_A1       835
#define OPT_A2       836
#define OPT_B1       837
#define OPT_B2       838
#define OPT_C1       839
#define OPT_C2       840
#define OPT_D1       841
#define OPT_D2       842
#define OPT_E1       843
#define OPT_E2       844
#define OPT_ONLY_S   845
#define OPT_ONLY_D   846
#define OPT_ONLY_C   847
#define OPT_CALIB    848
#define OPT_FIXED_L1 849
#define OPT_CUSTOM_SEQUIN_THRESHOLD 850
#define OPT_MMIX     851

// Shared with other modules
std::string __full_command__;

static Path __working__;

// Shared with other modules
Path __output__;

std::shared_ptr<FileWriter> __loggerW__;
std::shared_ptr<TerminalWriter> __outputW__;

static std::map<Value, Tool> _tools =
{
    { "rna",       Tool::RNA       },
    { "meta",      Tool::Meta      },
    { "norm",      Tool::Norm      },
    { "tmm",       Tool::TMM       },
    { "split",     Tool::Split     },
    { "germline",  Tool::Germline  },
    { "somatic",   Tool::Somatic   },
    { "cancer",    Tool::Cancer    },
    { "calibrate", Tool::Calibrate },
    { "broad_bam", Tool::BroadBAM  },
    { "broad_vcf", Tool::BroadVCF  }
};

struct Parsing
{
    Path path = "output";

    // Specific options
    std::map<Option, std::string> opts;
    
    std::map<Option, double> od;
    
    // How Anaquin is invoked
    Command  cmd;

    // Mixture A(1) or mixture B(2)?
    Mixture mix = Mix_1;
    
    Tool tool;
};

// Wrap the variables so that it'll be easier to reset them
static Parsing _p;

static Scripts fixManual(const Scripts &str)
{
    auto x = str;
    
    boost::replace_all(x, "<b>", "\e[1m");
    boost::replace_all(x, "</b>", "\e[0m");
    boost::replace_all(x, "<i>", "\e[3m");
    boost::replace_all(x, "</i>", "\e[0m");
    
    return x;
}

struct InvalidUsageException : public std::exception {};

struct InvalidDependency : public std::runtime_error
{
    InvalidDependency(const std::string &x) : std::runtime_error(x) {}
};

struct InvalidOptionException : public std::exception
{
    InvalidOptionException(const std::string &o) : o(o) {}

    const std::string o;
};

struct InvalidValueException : public std::exception
{
    InvalidValueException(const std::string &o, const std::string &v) : o(o), v(v) {}

    const std::string o, v;
};

struct InvalidToolError : public InvalidValueException
{
    InvalidToolError(const std::string &v) : InvalidValueException("-t", v) {}
};

struct MissingOptionError : public std::exception
{
    MissingOptionError(const std::string &o) : o(o) {}

    // Option that is missing
    const std::string o;
    
    // Possible values for the missing option
    const std::string r;
};

struct MissingBundleError : public std::runtime_error
{
    MissingBundleError(const Path &path) : std::runtime_error(path), path(path) {}
    const Path path;
};

struct InvalidBundleError : public std::runtime_error
{
    InvalidBundleError(const FileName &file) : std::runtime_error(file), file(file) {}
    const FileName file;
};

static const char *short_opts = ":";

static const struct option long_opts[] =
{
    { "1", required_argument, 0, OPT_1 },
    { "2", required_argument, 0, OPT_2 },
    { "3", required_argument, 0, OPT_3 },
    { "4", required_argument, 0, OPT_4 },

    { "A1", required_argument, 0, OPT_A1 },
    { "A2", required_argument, 0, OPT_A2 },
    { "B1", required_argument, 0, OPT_B1 },
    { "B2", required_argument, 0, OPT_B2 },
    { "C1", required_argument, 0, OPT_C1 },
    { "C2", required_argument, 0, OPT_C2 },
    { "D1", required_argument, 0, OPT_D1 },
    { "D2", required_argument, 0, OPT_D2 },
    { "E1", required_argument, 0, OPT_E1 },
    { "E2", required_argument, 0, OPT_E2 },
    
    { "v",   required_argument, 0, OPT_VCF },
    { "vcf", required_argument, 0, OPT_VCF },

    { "t",       required_argument, 0, OPT_THREAD },
    { "threads", required_argument, 0, OPT_THREAD },

    { "sample", required_argument, 0, OPT_SAMPLE },
    { "sequin", required_argument, 0, OPT_SEQUIN },
    { "sample_alignment", required_argument, 0, OPT_SAMPLE },
    { "sequin_alignment", required_argument, 0, OPT_SEQUIN },

    { "b",        required_argument, 0, OPT_COMBINE },
    { "combined", required_argument, 0, OPT_COMBINE },

    { "fq",      no_argument, 0, OPT_FASTQ },
    { "debug",   no_argument, 0, OPT_DEBUG },
    { "merge",   no_argument, 0, OPT_MERGE },

    { "sample_bam",            no_argument, 0, OPT_ONLY_S },
    { "sequin_bam",            no_argument, 0, OPT_ONLY_D },
    { "calibrated_sequin_bam", no_argument, 0, OPT_ONLY_C },

    { "human_regions",    required_argument, 0, OPT_R_HUMAN },
    { "decoy_regions",    required_argument, 0, OPT_R_DECOY },
    { "restrict_regions", required_argument, 0, OPT_R_REGS  },

    { "manual_mix",      required_argument, 0, OPT_MMIX  },
    { "manual_bed",      required_argument, 0, OPT_R_BED },
    { "manual_fasta",    required_argument, 0, OPT_FASTA },

    { "sample_variants", required_argument, 0, OPT_SAMPLE },
    { "vcf_hg38", required_argument, 0, OPT_SEQUIN },

    { "kmer",      required_argument, 0, OPT_KMER   },
    { "threshold", required_argument, 0, OPT_RULE   },
    { "skip",      required_argument, 0, OPT_SKIP   },
    
    { "r",            required_argument, 0, OPT_RESOURCE },
    { "resource_dir", required_argument, 0, OPT_RESOURCE },

    { "build", required_argument, 0, OPT_BUILD }, // Alternative genome assembly

    { "mix",    required_argument, 0, OPT_MIXTURE },
    { "method", required_argument, 0, OPT_METHOD  },
    { "calibration_method", required_argument, 0, OPT_METHOD },
    { "custom_sequin_threshold",  required_argument, 0, OPT_CUSTOM_SEQUIN_THRESHOLD },

    { "calibrate",        required_argument, 0, OPT_CALIB  },
    { "ladder_calibrate", required_argument, 0, OPT_LCALIB },
    
    { "edge",   required_argument, 0, OPT_EDGE   },
    { "window", required_argument, 0, OPT_WINDOW },

    { "o",      required_argument, 0, OPT_PATH },
    { "output", required_argument, 0, OPT_PATH },

    {0, 0, 0, 0 }
};

static Build parseBuild(const std::string &x)
{
    if      (x == "hg38") { return Build::hg38; }
    else if (x == "hg19") { return Build::hg19; }
    else if (x == "gr37") { return Build::gr37; }
    else if (x == "gr38") { return Build::gr38; }
    else
    {
        throw InvalidValueException("--build", "Invalid assembly build. Must be \"hg38\", \"hg19\", \"gr37\" or \"gr38\"");
    }
}

static std::string optToStr(int opt)
{
    for (const auto &o : long_opts)
    {
        if (o.val == opt)
        {
            return o.name;
        }
    }
    
    throw std::runtime_error("Invalid option: " + std::to_string(opt));
}

static void printUsage()
{
    extern Scripts Manual();
    std::cout << std::endl << fixManual(replace(Manual(), "%VERSION%", version())) << std::endl << std::endl;
}

static Scripts manual(Tool tool)
{
    switch (tool)
    {
        case Tool::RNA:       { return rna();       }
        case Tool::Meta:      { return meta();      }
        case Tool::Norm:      { return norm();      }
        case Tool::Split:     { return split();     }
        case Tool::Cancer:    { return cancer();    }
        case Tool::Somatic:   { return somatic();   }
        case Tool::Germline:  { return germline();  }
        case Tool::Calibrate: { return calibrate(); }
        case Tool::BroadBAM:  { return broadBAM();  }
        case Tool::BroadVCF:  { return broadVCF();  }
        default:              { return cancer();    }
    }
}

std::string __HACK_PRINT_CON__ = "";
bool __HACK_IS_CANCER__ = false;
bool __HACK_IS_CAPTURE__ = false;

inline std::string option(const Option &key, const Scripts &x = "")
{
    return _p.opts.count(key) ? _p.opts[key] : x;
}

template <typename F> std::shared_ptr<Ladder> readL(F f, Option key, UserReference &, const Scripts &x = "")
{
    return _p.opts.count(key) ? std::shared_ptr<Ladder>(new Ladder(f(Reader(_p.opts[key])))) :
                                std::shared_ptr<Ladder>(new Ladder(f(Reader(x))));
}

static std::shared_ptr<VCFLadder> readV(Option opt, UserReference &, std::shared_ptr<BedData> rb, const FileName &def)
{
    auto rr = (_p.opts.count(opt)) ? Reader(_p.opts[opt]) : Reader(def);
    return std::shared_ptr<VCFLadder>(new VCFLadder(Standard::addVCF(rr, rb)));
}

static std::shared_ptr<VCFLadder> readGV(Option opt, UserReference &, std::shared_ptr<BedData> rb, const FileName &def)
{
    auto rr = (_p.opts.count(opt)) ? Reader(_p.opts[opt]) : Reader(def);
    return std::shared_ptr<VCFLadder>(new VCFLadder(Standard::addGVCF(rr, rb)));
}

static std::shared_ptr<VCFLadder> readSV(Option opt, UserReference &, std::shared_ptr<BedData> rb, const FileName &def)
{
    auto rr = (_p.opts.count(opt)) ? Reader(_p.opts[opt]) : Reader(def);
    return std::shared_ptr<VCFLadder>(new VCFLadder(Standard::addSVCF(rr, rb)));
}

static std::shared_ptr<BedData> readR(const FileName &file, const RegionOptions &o = RegionOptions())
{
    return std::shared_ptr<BedData>(new BedData(Standard::readBED(Reader(file), o)));
}

static WriterOptions *__o__;

void info(const std::string &s) { if (__o__) { __o__->info(s); } }
void wait(const std::string &s) { if (__o__) { __o__->wait(s); } }

template <typename Analyzer, typename O, typename F> void start(const std::string &name, F f, O &o)
{
    o.info("-----------------------------------------");
    o.info("------------- Sequin Analysis -----------");
    o.info("-----------------------------------------");

    const auto path = _p.path;

    // Required for logging
    __o__ = &o;

    o.output = __outputW__;
    o.logger = __loggerW__;
    o.writer = std::shared_ptr<FileWriter>(new FileWriter(path));
    
    if (system(("mkdir -p " + path).c_str()))
    {
        throw std::runtime_error("Failed to create output directory");
    }
    
    o.name = name;
    o.debug = _p.opts.count(OPT_DEBUG);
    o.cmd = __full_command__ = _p.cmd;
    o.work = path;
    o.version = version();
    
    assert(!o.cmd.empty());
    assert(!o.name.empty());
    assert(!o.work.empty());

    o.info("Version: " + version());
    o.info(_p.cmd);
    o.info(date());
    o.info("Path: " + path);
    o.info("Resources: " + _p.opts[OPT_RESOURCE]);

    if (!__HACK_PRINT_CON__.empty())
    {
        o.info("Conversion from hg38 to chrQ completed for " + __HACK_PRINT_CON__);
    }
    
    using namespace std::chrono;
    
    auto begin = high_resolution_clock::now();

    f(o);
    
    auto end = high_resolution_clock::now();

    // Remove all working directories
    clearAllTmp();
    
    const auto elapsed = (boost::format("Completed. %1% seconds.") % duration_cast<seconds>(end - begin).count()).str();
    o.info(elapsed);

    o.logger->close();
}

template <typename Analyzer> void analyze_1(const std::string &p, Option x, typename Analyzer::Options &o = typename Analyzer::Options())
{
    return start<Analyzer>(p, [&](const typename Analyzer::Options &o)
    {
        Analyzer::report(_p.opts.at(x), o);
    }, o);
}

template <typename Analyzer> void analyze_2(const std::string &p, Option x1, Option x2, typename Analyzer::Options &o = typename Analyzer::Options())
{
    return start<Analyzer>(p, [&](const typename Analyzer::Options &o)
    {
        Analyzer::report(_p.opts.count(x1) ? _p.opts[x1] : "", _p.opts.count(x2) ? _p.opts[x2] : "", o);
    }, o);
}

template <typename T> std::shared_ptr<Ladder> readTSV(const Reader &r, T &, int con = 1)
{
    Ladder l = Ladder();
    
    ParserCSV::parse(r, [&](const ParserCSV::Data &x, Progress)
    {
        if (x.size() < 2)
        {
            throw std::runtime_error("Invalid format. Two or more columns expected");
        }
        else if (x[con] == MISSING || !isNumber(x[con]))
        {
            return;
        }
        
        l.add(x[0], Mix_1, stof(x[con]));
    }, "\t");

    return std::shared_ptr<Ladder>(new Ladder(l));
}

void parse(int argc, char ** argv)
{
    auto tmp = new char*[argc+1];
    
    for (auto i = 0; i < argc; i++)
    {
        tmp[i] = (char *) malloc((strlen(argv[i]) + 1) * sizeof(char));
        strcpy(tmp[i], argv[i]);
    }
    
    _p = Parsing();

    if ((argc <= 1) || (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")))
    {
        printUsage();
        return;
    }

    // Required for unit-testing
    optind = 1;

    /*
     * Reconstruct the overall command
     */
    
    for (auto i = 0; i < argc; i++)
    {
        _p.cmd += std::string(argv[i]) + " ";
    }

    assert(!_p.cmd.empty());

    int next, index;

    auto checkPath = [&](const Path &path)
    {
        if (path[0] == '/')
        {
            return path;
        }
        else
        {
            return __working__ + "/" + path;
        }
    };
    
    auto checkBundle = [&](const Path &path)
    {
        if (!exists(path))
        {
            throw MissingBundleError(path);
        }
        
        const auto files = listFiles(path);
        
        auto check = [&](const FileName &x)
        {
            if (std::find_if(files.begin(), files.end(), [&](const FileName &i) {
                return isBegin(i, x);
            }) == files.end()) { throw InvalidBundleError(x); }
        };
        
        check("genome");
        check("synthetic");
        check("metagenome");
        check("transcriptome");

        return path;
    };

    auto checkFile = [&](const FileName &file)
    {
        if (!std::ifstream(file).good())
        {
            throw InvalidFileError(file);
        }
        
        return file;
    };

    /*
     * Pre-process arguments. This way, we can examine the options in whatever order we'd like to
     */

    std::vector<Value>  vals;
    std::vector<Option> opts;

    if (strcmp(argv[1], "-v") == 0)
    {
        std::cout << version() << std::endl;
        return;
    }
    else if (strcmp(argv[1], "-t") == 0)
    {
#ifdef UNIT_TEST
        Catch::Session().run(1, argv);
#else
        A_THROW("UNIT_TEST is undefined");
#endif
        return;
    }
    else if (!_tools.count(argv[1]))
    {
        throw InvalidToolError(argv[1]);
    }

    _p.tool = _tools[argv[1]];
    const auto isHelp = argc >= 3 && (!strcmp(argv[2], "-h") || !strcmp(argv[2], "--help"));

    if (isHelp)
    {
        if (argc != 3)
        {
            throw std::runtime_error("Too many arguments for help usage. Usage: anaquin <tool> -h or anaquin <tool> --help");
        }

        std::cout << std::endl << fixManual(manual(_p.tool)) << std::endl << std::endl;
        return;
    }
    
    unsigned n = 2;

    while ((next = getopt_long_only(argc, argv, short_opts, long_opts, &index)) != -1)
    {
        if (next < OPT_TOOL)
        {
            throw InvalidOptionException(argv[n]);
        }
        
        opts.push_back(next);
        
        // Whether this option has an value
        const auto hasValue = optarg;
        
        n += hasValue ? 2 : 1;
        
        vals.push_back(hasValue ? std::string(optarg) : "");
    }

    for (auto i = 0u; i < opts.size(); i++)
    {
        auto opt = opts[i];
        auto val = vals[i];

        switch (opt)
        {
            case OPT_RULE:
            case OPT_CUSTOM_SEQUIN_THRESHOLD:
            {
                try
                {
                    stof(val);
                    _p.opts[opt] = val;
                }
                catch (...)
                {
                    throw std::runtime_error(val + " is not an. Please check and try again.");
                }
                
                break;
            }

            case OPT_CALIB:
            case OPT_LCALIB:
            {
                try
                {
                    _p.od[opt] = stod(val);
                }
                catch (...)
                {
                    throw std::runtime_error(val + " is not a floating number");
                }

                break;
            }
            
            case OPT_EDGE:
            case OPT_KMER:
            case OPT_SKIP:
            case OPT_WINDOW:
            {
                try
                {
                    stoi(val);
                    _p.opts[opt] = val;
                }
                catch (...)
                {
                    throw std::runtime_error(val + " is not an integer. Please check and try again.");
                }

                break;
            }

            case OPT_METHOD:
            {
                switch (_p.tool)
                {
                    case Tool::Cancer:
                    case Tool::Somatic:
                    case Tool::Germline:
                    case Tool::Calibrate: { _p.opts[opt] = val; break; }
                    default : { throw InvalidOptionException("Invalid usage for --method"); }
                }
                
                break;
            }

            case OPT_FASTQ:
            case OPT_FASTA:
            case OPT_BUILD:
            case OPT_DEBUG:
            case OPT_MERGE:
            case OPT_THREAD:
            case OPT_ONLY_S:
            case OPT_ONLY_D:
            case OPT_ONLY_C:
            {
                _p.opts[opt] = val;
                break;
            }

            case OPT_MIXTURE:
            {
                if      (val == "A") { _p.mix = Mixture::Mix_1; }
                else if (val == "B") { _p.mix = Mixture::Mix_2; }
                else if (val == "C") { _p.mix = Mixture::Mix_3; }
                else                 { throw InvalidValueException("-mix", val); }
                break;
            }

            case OPT_1:
            case OPT_2:
            case OPT_3:
            case OPT_4:
            case OPT_A1:
            case OPT_A2:
            case OPT_B1:
            case OPT_B2:
            case OPT_C1:
            case OPT_C2:
            case OPT_D1:
            case OPT_D2:
            case OPT_E1:
            case OPT_E2:
            case OPT_VCF:
            case OPT_SAMPLE:
            case OPT_SEQUIN:
            case OPT_COMBINE: { checkFile(_p.opts[opt] = val); break; }

            case OPT_RESOURCE: { checkBundle(_p.opts[opt] = val); break; }
                
            case OPT_MMIX:
            case OPT_R_VCF:
            case OPT_R_BED:
            case OPT_R_REGS:
            case OPT_R_HUMAN:
            case OPT_R_DECOY:
            {
                checkFile(_p.opts[opt] = val);
                break;
            }

            case OPT_PATH: { _p.path = val; break; }

            default: { throw InvalidUsageException(); }
        }
    }
    
    if (!_p.opts.count(OPT_RESOURCE))
    {
        _p.opts[OPT_RESOURCE] = checkBundle(execPath() + "/resources");
    }

    __output__ = _p.path = checkPath(_p.path);
    
    if (opts.empty())
    {
        std::cout << std::endl << fixManual(manual(_p.tool)) << std::endl << std::endl;
        return;
    }
    
    UserReference r;
    
    auto &s = Standard::instance();
    
    __loggerW__ = std::shared_ptr<FileWriter>(new FileWriter(_p.path));
    __outputW__ = std::shared_ptr<TerminalWriter>(new TerminalWriter());
    __loggerW__->open("anaquin.log");
    
    #define CSFA()          (CSeqFA(_p.opts.at(OPT_RESOURCE)).path)
    #define GSFA()          (GSeqFA(_p.opts.at(OPT_RESOURCE)).path)
    #define GDFA()          (GSeqDecoy(_p.opts.at(OPT_RESOURCE)).path)
    #define CDFA()          (CSeqDecoy(_p.opts.at(OPT_RESOURCE)).path)
    #define GFeatBED_(x)    (GFeatBED(_p.opts.at(OPT_RESOURCE), x).path)
    #define GLTSV(x)        (GSynTSV(_p.opts.at(OPT_RESOURCE)).path)
    #define GAttrBED_(x)    (GAttrBED(_p.opts.at(OPT_RESOURCE)).path)
    #define CRegionBED_(x)  (CRegionBED(_p.opts.at(OPT_RESOURCE), x).path)
    #define GRegionBED_(x)  (GRegionBED(_p.opts.at(OPT_RESOURCE), x).path)
    #define CVCF(x)         (CVarVCF(_p.opts.at(OPT_RESOURCE), x).path)
    #define GVCF(x)         (GVarVCF(_p.opts.at(OPT_RESOURCE), x).path)
    #define RSFA()          (RNAFA(_p.opts.at(OPT_RESOURCE)).path)
    #define RDFA()          (RNADecoy(_p.opts.at(OPT_RESOURCE)).path)

    auto initAR = [&](UserReference &r)
    {
        const auto a1 = readAttrib(GFeatBED_(Build::hseq));
        const auto a2 = readAttrib(GFeatBED_(Build::hg38));
        const auto a3 = readAttrib(GFeatBED_(Build::chrQ));

        r.a1 = std::shared_ptr<AttributeBed>(new AttributeBed(a1));
        r.a2 = std::shared_ptr<AttributeBed>(new AttributeBed(a2));
        r.a3 = std::shared_ptr<AttributeBed>(new AttributeBed(a3));
    };

    const auto rpath = _p.opts.at(OPT_RESOURCE);
    
    switch (_p.tool)
    {
        case Tool::RNA:
        case Tool::TMM:
        case Tool::Meta:
        case Tool::Norm:
        case Tool::Split:
        case Tool::Cancer:
        case Tool::Somatic:
        case Tool::Germline:
        case Tool::BroadBAM:
        case Tool::BroadVCF:
        case Tool::Calibrate:
        {
            switch (_p.tool)
            {
                case Tool::Split:
                case Tool::Cancer:
                case Tool::BroadBAM:
                case Tool::Calibrate:
                {
                    const auto edge  = _p.opts.count(OPT_EDGE)  ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    const auto build = _p.opts.count(OPT_BUILD) ? parseBuild(_p.opts.at(OPT_BUILD)) :  (_p.opts.count(OPT_COMBINE) ? Detector::fromBAM(_p.opts.at(OPT_COMBINE)) : Build::hg38);
                    
                    const auto bothHR = _p.opts.count(OPT_SAMPLE) && _p.opts.count(OPT_SEQUIN);
                    const bool isCancer = __HACK_IS_CANCER__ = _p.tool == Tool::Cancer;

                    const auto gb1 = isCancer ? CRegionBED_(build) : GRegionBED_(build);
                    const auto gb2 = isCancer ? CRegionBED_(Build::chrQ) : GRegionBED_(Build::chrQ);
                    
                    const auto hr = _p.opts.count(OPT_R_BED)   ? _p.opts[OPT_R_BED] :
                                    _p.opts.count(OPT_R_HUMAN) ? _p.opts[OPT_R_HUMAN] : gb1;
                    const auto dr = _p.opts.count(OPT_R_BED)   ? _p.opts[OPT_R_BED] :
                                    _p.opts.count(OPT_R_DECOY) ? _p.opts[OPT_R_DECOY] : !bothHR ? gb2 : gb1;

                    // Restricted regions
                    auto rr = _p.opts.count(OPT_R_REGS) ? _p.opts[OPT_R_REGS] : dr;

                    if (dr == rr)
                    {
                        RegionOptions o;
                        r.r1 = readR(hr, o); // No trimming
                        r.r2 = readR(dr, o); // No trimming
                        
                        o.edge = edge;
                        r.r3 = readR(hr, o); // Trimmed
                        r.r4 = readR(dr, o); // Trimmed
                    }
                    else // Restricted regions provided
                    {
                        auto oldRR  = rr;
                        auto isHG38 = false;
                        
                        // Decoy mode? If so, maybe we should try to convert to chrQ?
                        if (_p.opts.count(OPT_COMBINE))
                        {
                            const auto tmp = tmpFile();
                            runScript(hg382chrQ(), hr + " " + dr + " " + rr + " > " + tmp);
                            
                            auto nChrQ = 0;
                            auto nHG38 = 0;
                            
                            ParserCSV::parse(Reader(rr), [&](const ParserCSV::Data &x, Progress i)
                            {
                                if (isSubstr(x[0], "chrQ"))
                                {
                                    nChrQ++;
                                }
                                else
                                {
                                    nHG38++;
                                }
                            }, "\t");
                            
                            // Is this really a hg38 lift-over file?
                            if (nHG38 > nChrQ)
                            {
                                isHG38 = true;
                            }
                            
                            // Use the chrQ version instead of hg38 (for lift-over)
                            __HACK_PRINT_CON__ = rr = tmp;
                            
                            __HACK_IS_CAPTURE__ = true;
                            __HACK_PRINT_CON__ += (" to " + rr);
                        }
                        
                        RegionOptions o;
                        r.r1 = readR(hr, o); // No trimming
                        r.r2 = readR(dr, o); // No trimming
                        
                        o.edge = edge;
                        
                        if (isHG38)
                        {
                            r.r3 = readR(BedTools::intersect(hr, oldRR, edge), o); // Trimmed
                        }
                        else
                        {
                            r.r3 = readR(hr, o); // Trimmed
                        }
                        
                        r.r4 = readR(BedTools::intersect(dr, rr, edge), o); // Trimmed
                    }
                    
                    assert(!r.r1->inters().empty());
                    assert(!r.r2->inters().empty());
                    assert(!r.r3->inters().empty());
                    assert(!r.r4->inters().empty());

                    initAR(r);
                    r.t1 = readTrans(Reader(GInfoCode(_p.opts.at(OPT_RESOURCE)).path));
                    
                    // Attributes on chrQS (hg19/hg38 not supported)
                    r.r5 = readR(GFeatBED_(!bothHR ? chrQ : build), RegionOptions()); // Never trimming
                    
                    const auto hg38V = isCancer ? CVCF(Build::hg38) : GVCF(Build::hg38);
                    const auto chrQV = isCancer ? CVCF(Build::chrQ) : GVCF(Build::chrQ);
                    
                    /*
                     * Allele frequency ladder, it's assumed that chrQ and hg38 encode the same AF so only the hg38
                     * file is used. Locations will not be read.
                     */
                    
                    r.l1 = readL(std::bind(&Standard::addAF, &s, std::placeholders::_1), OPT_R_VCF, r, hg38V);
                    r.l2 = readTSV(Reader(option(OPT_ATTTSV, GAttrBED_())), r); // Attributes
                    r.l3 = readTSV(Reader(option(OPT_SYNC, GLTSV())), r, 2); // Synthetuc ladder
                    
                    r.v1 = readV(OPT_R_VCF,  r, nullptr, hg38V); // All variants
                    r.v2 = readGV(OPT_R_VCF, r, nullptr, hg38V); // Germline variants
                    r.v3 = readSV(OPT_R_VCF, r, nullptr, hg38V); // Somatic variants
                    r.v4 = readV(OPT_R_VCF,  r, nullptr, chrQV); // All decoy variants

                    break;
                }

                case Tool::Somatic:
                case Tool::Germline:
                case Tool::BroadVCF:
                {
                    const auto useDecoy = _p.opts.count(OPT_VCF);

                    // Key for checking automatic detection
                    const auto key = useDecoy ? OPT_VCF : OPT_SEQUIN;
                    
                    // File for automatic detection
                    const auto file = _p.opts.at(key);

                    // Attempt to detect build from sequin VCF
                    __ha__ = _p.opts.count(OPT_BUILD) ? parseBuild(_p.opts.at(OPT_BUILD)) : Detector::fromVCF(file);

                    // Human BED regions (not used if decoy)
                    const auto hr = _p.opts.count(OPT_R_HUMAN) ? _p.opts[OPT_R_HUMAN] : GRegionBED_(__ha__);

                    // Decoy BED regions
                    const auto dr = !useDecoy ? hr : (_p.opts.count(OPT_R_DECOY) ? _p.opts[OPT_R_DECOY] : GRegionBED_(Build::chrQ));

                    const auto rr   = _p.opts.count(OPT_R_REGS) ? _p.opts[OPT_R_REGS] : dr;
                    const auto edge = _p.opts.count(OPT_EDGE)   ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;

                    /*
                     * For combined,     a1 == chrQS and a2 == hg38
                     * For non-combined, a1 == hg38  and a2 == null
                     */
                    
                    try
                    {
                        r.a1 = std::shared_ptr<AttributeBed>(new AttributeBed(readAttrib(
                                        BedTools::intersect2(BedTools::intersect(dr, rr, 0),
                                              useDecoy ? GFeatBED_(Build::chrQ) : GFeatBED_(__ha__)))));
                    }
                    catch (const ZeroIntersectionError &)
                    {
                        r.a1 = std::shared_ptr<AttributeBed>(new AttributeBed());
                    }
                    
                    // TODO: Why are we intersecting the attributes?
                    r.a1->src = useDecoy ? GFeatBED_(Build::chrQ) : GFeatBED_(__ha__);
                    
                    if (useDecoy)
                    {
                        try
                        {
                            r.a2 = std::shared_ptr<AttributeBed>(new AttributeBed(
                                        readAttrib(BedTools::intersect2(hr, GFeatBED_(__ha__)))));
                        }
                        catch (const ZeroIntersectionError &)
                        {
                            r.a2 = std::shared_ptr<AttributeBed>(new AttributeBed());
                        }
                    }

                    /*
                     * r1 = human regions without edge
                     * r2 = human regions with edge
                     * r3 = decoy regions without edge
                     * r4 = decoy regions with edge
                     */
                    
                    RegionOptions o1, o2, o3, o4;
                    o2.edge = edge;
                    o4.edge = edge;

                    if (useDecoy)
                    {
                        o3.onlyC.insert(GDecoyChrQS); // No region other than chrQS should be included
                        o4.onlyC.insert(GDecoyChrQS); // No region other than chrQS should be included
                    }

                    if (dr == rr)
                    {
                        r.r1 = readR(hr, o1); // Human regions without edge
                        r.r2 = readR(hr, o2); // Human regions with edge
                        r.r3 = readR(dr, o3); // Decoy regions without edge
                        r.r4 = readR(dr, o4); // Decoy regions with edge
                    }
                    else
                    {
                        const auto f1 = BedTools::intersect(dr, rr, 0);
                        const auto f2 = BedTools::intersect(dr, rr, edge);
                        
                        if (useDecoy)
                        {
                            r.r1 = readR(hr, o1); // Human regions without edge
                            r.r2 = readR(hr, o2); // Human regions with edge
                        }
                        else
                        {
                            r.r1 = readR(f1, o1); // Human regions without edge
                            r.r2 = readR(f2, o1); // Human regions with edge
                        }
                        
                        r.r3 = readR(f1, o3); // Decoy regions without edge
                        r.r4 = readR(f2, o4); // Decoy regions with edge
                    }
                    
                    r.r5 = readR(dr, o4); // Decoy regions without intersection with edge

                    const auto vcf = useDecoy ? GVCF(Build::chrQ) : GVCF(__ha__);
                    
                    r.v1 = readV (OPT_R_VCF, r, r.r4, vcf); // All variants
                    r.v2 = readGV(OPT_R_VCF, r, r.r4, vcf); // Germline variants
                    r.v3 = readSV(OPT_R_VCF, r, r.r4, vcf); // Somatic variants

                    break;
                }

                case Tool::Meta:
                {
                    #define MMix()  (MetaMix(_p.opts.at(OPT_RESOURCE)).path)
                    #define MQBed() (MetaDBED(_p.opts.at(OPT_RESOURCE)).path)
                    
                    r.l1 = readL(std::bind(&Standard::readMMix, &s, std::placeholders::_1), OPT_MMIX, r, MMix());
                    r.l3 = readTSV(Reader(option(OPT_SYNC, GLTSV())), r, 2); // Unit
                    r.t1 = readTrans(Reader(GInfoCode(_p.opts.at(OPT_RESOURCE)).path));

                    const auto edge = _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    RegionOptions o2;
                    o2.edge = edge;

                    r.r1 = readR(MQBed());
                    r.r2 = readR(MQBed(), o2);

                    break;
                }

                case Tool::RNA:
                {
                      #define RMix()  (RNAMix(_p.opts.at(OPT_RESOURCE)).path)
                      #define RTBed() (RNATBed(_p.opts.at(OPT_RESOURCE)).path)
                      #define RGBed() (RNAGBed(_p.opts.at(OPT_RESOURCE)).path)

                    r.l1 = readL(std::bind(&Standard::readRMix,  &s, std::placeholders::_1), OPT_MMIX, r, RMix());
                    r.l2 = readL(std::bind(&Standard::readRGMix, &s, std::placeholders::_1), OPT_MMIX, r, RMix());
                    r.l3 = readL(std::bind(&Standard::readRLen,  &s, std::placeholders::_1), OPT_MMIX, r, RMix());

                    const auto edge  = _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    RegionOptions o_;
                    o_.edge = edge;

                    r.r1 = readR(RTBed());     // Non-trimming transcripts
                    r.r2 = readR(RTBed(), o_); // Trimmed transcripts
                    r.r3 = readR(RGBed());     // Non-trimming genes
                    r.r4 = readR(RGBed(), o_); // Trimmed genes

                    break;
                }

                case Tool::TMM:
                case Tool::Norm:
                {
                    RegionOptions o2;
                    o2.edge = DEFAULT_EDGE;
                    
                    const auto build = _p.opts.count(OPT_BUILD) ? parseBuild(_p.opts.at(OPT_BUILD)) : Build::hg38;
                    
                    // Attribute regions
                    initAR(r);
                    
                    r.l1 = readL(std::bind(&Standard::addAF, &s, std::placeholders::_1), OPT_R_VCF, r, GVCF(build));
                    r.l2 = readTSV(Reader(option(OPT_ATTTSV, GAttrBED_())), r);
                    r.l3 = readTSV(Reader(option(OPT_SYNC, GLTSV())), r, 2); // Unit
                    r.t1 = readTrans(Reader(GInfoCode(_p.opts.at(OPT_RESOURCE)).path));

                    r.v1 = readV (OPT_R_VCF, r, r.r2, GVCF(build)); // All variants
                    r.v2 = readGV(OPT_R_VCF, r, r.r2, GVCF(build)); // Germline variants
                    r.v3 = readSV(OPT_R_VCF, r, r.r2, GVCF(build)); // Somatic variants

                    /*
                     * l1 is always the base quantitative ladder. We should check if our mixture
                     * is compatible or not.
                     */
                    
                    if (_p.mix == Mix_2 && !r.l1->hasMix2())
                    {
                        throw InvalidOptionException("Mixture B specified, but only a single mixture found in reference file: " + r.l1->src);
                    }
                    else if (_p.mix == Mix_3 && !r.l1->hasMix3())
                    {
                        throw InvalidOptionException("Mixture C specified, but only a single mixture found in reference file: " + r.l1->src);
                    }

                    break;
                }

                default: { break; }
            }

            Standard::instance().gen.finalize(r);
            Standard::instance().rna.finalize(r);
            Standard::instance().meta.finalize(r);

            auto initSOptions = [&](SOptions &o, const FileName &i1, const FileName &i2)
            {
                assert(o.flip && o.skipMerge);

                o.bam = _p.opts.count(OPT_COMBINE);

                // How many k-mers to skip?
                o.skipKM = _p.opts.count(OPT_SKIP) ? stoi(_p.opts[OPT_SKIP]) : 5;
                
                // Number of threads
                o.thr = _p.opts.count(OPT_THREAD) ? stoi(_p.opts[OPT_THREAD]) : 1;
                
                // K-mer length
                o.k = _p.opts.count(OPT_KMER) ? stoi(_p.opts[OPT_KMER]) : K_DEFAULT_L;
                
                // For debugging?
                o.debug = _p.opts.count(OPT_DEBUG);
                
                // Classification rule
                o.rule = _p.opts.count(K_DEFAULT_R) ? stoi(_p.opts[K_DEFAULT_R]) : K_DEFAULT_R;
                
                o.index = _p.opts.count(OPT_FASTA) ? _p.opts[OPT_FASTA] : (!o.bam ? i1 : i2);
            };
            
            auto initCalib = [&](int key, double &x)
            {
                if (!_p.od.count(key))
                {
                    return;
                }
                else if (_p.od.at(key) <= 0)
                {
                    throw std::runtime_error("Calibration value must be greater than zero");
                }
                
                x = _p.od.at(key);
            };
            
            switch (_p.tool)
            {
                case Tool::Meta:
                {
                    MSplit::Options o;
                    initSOptions(o, MetaFA(_p.opts.at(OPT_RESOURCE)).path, MetaFA(_p.opts.at(OPT_RESOURCE)).path);
                    if (o.bam) { o.index = MetaDecoy(_p.opts.at(OPT_RESOURCE)).path; }

                    o.mix = _p.mix;
                    assert(o.seqC == NO_CALIBRATION && o.ladC == NO_CALIBRATION);

                    double tmp = -1;
                    
                    initCalib(OPT_CALIB,  o.seqC);
                    initCalib(OPT_LCALIB, o.ladC);
                    
                    o.edge = _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    
                    start<MSplit>("meta", [&](const MSplit::Options &)
                    {
                        if (o.bam) { MSplit::report(_p.opts[OPT_COMBINE], "", o); }
                        else       { MSplit::report(_p.opts[OPT_1], _p.opts.count(OPT_2) ? _p.opts[OPT_2] : _p.opts[OPT_1], o); }
                    }, o);

                    break;
                }
                    
                case Tool::RNA:
                {
                    RSplit::Options o;
                    initSOptions(o, RSFA(), RDFA());
                    
                    o.mix = _p.mix;
                    assert(o.seqC == NO_CALIBRATION && o.ladC == NO_CALIBRATION);
                    
                    initCalib(OPT_CALIB, o.seqC);
                    
                    start<RSplit>("rna", [&](const RSplit::Options &)
                    {
                        if (o.bam) { RSplit::report(_p.opts[OPT_COMBINE], "", o); }
                        else       { RSplit::report(_p.opts[OPT_1], _p.opts.count(OPT_2) ? _p.opts[OPT_2] : _p.opts[OPT_1], o); }
                    }, o);

                    break;
                }
                    
                case Tool::Split:
                {
                    GSplit::Options o;

                    initSOptions(o, GSFA(), GDFA());
                    initCalib(OPT_CALIB,  o.seqC);
                    initCalib(OPT_LCALIB, o.ladC);

                    start<GSplit>("split", [&](const GSplit::Options &)
                    {
                        if (o.bam) { GSplit::report(_p.opts[OPT_COMBINE], "", o); }
                        else       { GSplit::report(_p.opts[OPT_1], _p.opts.count(OPT_2) ? _p.opts[OPT_2] : _p.opts[OPT_1], o); }
                    }, o);

                    break;
                }

                case Tool::TMM:
                {
                    GTMM::Options o;
                    
                    start<GTMM>("tmm", [&](const GSplit::Options &)
                    {
                        GTMM::report(_p.opts[OPT_1], o);
                    }, o);

                    break;
                }

                case Tool::Norm:
                {
                    GNorm::Options o;
                    initSOptions(o, GSFA(), GSFA());
                    assert(!o.bam);
                    
                    std::vector<FileName> f1, f2;
                    
                    if (_p.opts.count(OPT_A1)) { f1.push_back(_p.opts.at(OPT_A1)); }
                    if (_p.opts.count(OPT_B1)) { f1.push_back(_p.opts.at(OPT_B1)); }
                    if (_p.opts.count(OPT_C1)) { f1.push_back(_p.opts.at(OPT_C1)); }
                    if (_p.opts.count(OPT_D1)) { f1.push_back(_p.opts.at(OPT_D1)); }
                    if (_p.opts.count(OPT_E1)) { f1.push_back(_p.opts.at(OPT_E1)); }
                    if (_p.opts.count(OPT_A2)) { f2.push_back(_p.opts.at(OPT_A2)); }
                    if (_p.opts.count(OPT_B2)) { f2.push_back(_p.opts.at(OPT_B2)); }
                    if (_p.opts.count(OPT_C2)) { f2.push_back(_p.opts.at(OPT_C2)); }
                    if (_p.opts.count(OPT_D2)) { f2.push_back(_p.opts.at(OPT_D2)); }
                    if (_p.opts.count(OPT_E2)) { f2.push_back(_p.opts.at(OPT_E2)); }
                    
                    start<GNorm>("norm", [&](const GSplit::Options &)
                    {
                        GNorm::report(f1, f2, o);
                    }, o);

                    break;
                }

                case Tool::BroadVCF:
                {
                    GVariant::Options o;
                    
                    o.edge  = _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    o.uBED  = _p.opts.count(OPT_R_REGS) ? _p.opts[OPT_R_REGS] : MISSING;
                    o.decoy = true;
                    
                    start<GBroadVCF>("broadVCF", [&](const GVariant::Options &)
                    {
                        GBroadVCF::report(_p.opts.at(OPT_VCF), o);
                    }, o);

                    break;
                }
                    
                case Tool::Somatic:
                case Tool::Germline:
                {
                    GVariant::Options o;
                    
                    o.decoy = _p.opts.count(OPT_VCF);
                    o.edge  = _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    o.uBED  = _p.opts.count(OPT_R_REGS) ? _p.opts[OPT_R_REGS] : MISSING;
                    
                    const auto f1 = o.decoy ? _p.opts.at(OPT_VCF) : _p.opts.count(OPT_SAMPLE) ? _p.opts.at(OPT_SAMPLE) : "";
                    const auto f2 = o.decoy ? "" : _p.opts.at(OPT_SEQUIN);
                    
                    if (_p.tool == Tool::Germline)
                    {
                        start<GGerm>(o.base = "germline", [&](const GVariant::Options &)
                        {
                            GGerm::report(f1, f2, o);
                        }, o);
                    }
                    else
                    {
                        start<GSomatic>(o.base = "somatic", [&](const GVariant::Options &)
                        {
                            GSomatic::report(f1, f2, o);
                        }, o);
                    }

                    break;
                }

                case Tool::BroadBAM:
                {
                    GBroadBam::Options o;
                    
                    o.seqC  = _p.od.count(OPT_CALIB) ? _p.od.at(OPT_CALIB) : -1;
                    o.index = GSeqDecoy(_p.opts.at(OPT_RESOURCE)).path;
                    analyze_1<GBroadBam>("broadBAM", OPT_COMBINE, o);

                    break;
                }

                case Tool::Cancer:
                case Tool::Calibrate:
                {
                    typedef GCalibrate::Method Method;
                    const auto isChrQ = _p.opts.count(OPT_COMBINE);
                    const auto isCancer = _p.tool == Tool::Cancer;
                    
                    const auto f1 = isCancer ? CDFA() : GDFA();
                    const auto f2 = isCancer ? CSFA() : GSFA();
                    
                    GCalibrate::Options o;
                    initSOptions(o, isChrQ ? f1 : f2, isChrQ ? f1 : f2);
                    o.edge = _p.opts.count(OPT_EDGE) ? stoi(_p.opts[OPT_EDGE]) : DEFAULT_EDGE;
                    
                    o.writeS = _p.opts.count(OPT_ONLY_S);
                    o.writeD = _p.opts.count(OPT_ONLY_D);
                    o.writeC = _p.opts.count(OPT_ONLY_C);
                    o.isCancer = isCancer;

                    // How to calibrate?
                    const auto meth = _p.opts.count(OPT_METHOD) ? _p.opts[OPT_METHOD] : "mean";
                    
                    if      (meth == "mean")   { o.meth = Method::Mean;   }
                    else if (meth == "median") { o.meth = Method::Median; }
                    else if (meth == "custom")
                    {
                        if (_p.opts.count(OPT_CUSTOM_SEQUIN_THRESHOLD) == 0)
                        {
                            throw MissingOptionError("--custom_sequin_threshold");
                        }
                        
                        o.meth = Method::Custom;
                        o.customSequinThreshold = stof(_p.opts.at(OPT_CUSTOM_SEQUIN_THRESHOLD));
                    }
                    else if (isFloat(meth)) { o.meth = Method::Percent; o.seqC = stof(meth); }
                    else
                    {
                        throw std::runtime_error("Unknown method: " + meth);
                    }
                    
                    if (o.bam)
                    {
                        analyze_1<GCalibrate>("calibrate", OPT_COMBINE, o);
                    }
                    else
                    {
                        if (o.meth != Method::Percent && !_p.opts.count(OPT_SAMPLE))
                        {
                            throw MissingOptionError("--sample");
                        }
                        
                        analyze_2<GCalibrate>("calibrate", o.meth != Method::Percent ? OPT_SAMPLE : OPT_SEQUIN, OPT_SEQUIN, o);
                    }
                    
                    break;
                }

                default : { break; }
            }

            break;
        }

        default : { break; }
    }
}

extern int parse_options(int argc, char ** argv)
{
    char cwd[1024];
    
    auto printError = [&](const std::string &x)
    {
        std::cerr << "***********************" << std::endl;
        std::cerr << "[ERRO]: " << x << std::endl;
        std::cerr << "***********************" << std::endl << std::endl;
    };
    
    if (getcwd(cwd, sizeof(cwd)))
    {
        __working__ = cwd;
    }
    
    try
    {
        parse(argc, argv);
        return 0;
    }
    catch (const FailedCommandException &ex)
    {
        printError(std::string(ex.what()));
    }
    catch (const InvalidFormatException &ex)
    {
        printError("Invalid file format: " + std::string(ex.what()));
    }
    catch (const InvalidUsageException &)
    {
        printError("Invalid usage. Please check and try again.");
    }
    catch (const InvalidToolError &ex)
    {
        printError("Invalid command. Unknown tool: " + ex.v + ". Please check your usage and try again.");
    }
    catch (const InvalidOptionException &ex)
    {
        printError((boost::format("Invalid usage. Unknown option: %1%") % ex.o).str());
    }
    catch (const InvalidValueException &ex)
    {
        printError((boost::format("Invalid command. %1% not expected for %2%.") % ex.v % ex.o).str());
    }
    catch (const MissingOptionError &ex)
    {
        const auto format = "Invalid command. Mandatory option is missing. Please specify %1%.";
        printError((boost::format(format) % ex.o).str());
    }
    catch (const MissingBundleError &ex)
    {
        printError((boost::format("%1%%2%") % "Failed to find bundle at: " % ex.path).str());
    }
    catch (const InvalidBundleError &ex)
    {
        printError((boost::format("%1%%2%") % "Invalid bundle. Missing bundle folder: " % ex.file).str());
    }
    catch (const InvalidFileError &ex)
    {
        printError((boost::format("%1%%2%") % "Invalid command. File is invalid: " % ex.file).str());
    }
    catch (const InvalidDependency &ex)
    {
        printError(ex.what());
    }
    catch (const std::runtime_error &ex)
    {
        printError(ex.what());
    }

    return 1;
}

int main(int argc, char ** argv)
{
    srand(time(0));
    return parse_options(argc, argv);
}
