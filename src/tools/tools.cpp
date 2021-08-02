#include <mutex>
#include <fstream>
#include <iostream>
#include <dirent.h>
#include <libgen.h>
#include <sys/stat.h>
#include <gzip/utils.hpp>
#include "data/reader.hpp"
#include "tools/tools.hpp"
#include "tools/whereami.h"
#include "tools/errors.hpp"
#include <gzip/compress.hpp>
#include <gzip/decompress.hpp>
#include "parsers/parser_fq.hpp"
#include "writers/fq_writer.hpp"
#include "writers/bam_writer.hpp"
#include <boost/iostreams/copy.hpp>
#include <boost/algorithm/string.hpp>

using namespace Anaquin;

bool Anaquin::isEmpty(const FileName& file)
{
    std::ifstream ss(file);
    return ss.peek() == std::ifstream::traits_type::eof();
}

void Anaquin::createD(const Path &x)
{
    mkdir(x.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

void Anaquin::removeD(const Path &x)
{
    const auto cmd = "rm -rf " + x;
    
    if (system(cmd.c_str()))
    {
        throw FailedCommandException("Failed: " + cmd);
    }
}

void Anaquin::Compressor::write()
{
    if (_i) { *_w << compressGZ(_buf.str()); _i = 0; _buf.clear(); }
}

void Anaquin::Compressor::write(const std::string &s)
{
    _buf << s;
    if (++_i >= _n) { write(); assert(_i == 0); }
}

static std::vector<FileName> __TMP_PATHS__;

void Anaquin::clearAllTmp()
{
    for (const auto &i : __TMP_PATHS__) { removeD(i); }
    __TMP_PATHS__.clear();
}

Path Anaquin::tmpPath()
{
    static int i = 0;
    const auto path = "tmp_anaquin_" + S0(getpid()) + std::to_string(i++);
    createD(path);
    __TMP_PATHS__.push_back(path);
    return path;
}

void Anaquin::mergeFQ(const std::vector<FileName> &s1,
                      const std::vector<FileName> &s2,
                      const FileName &o1,
                      const FileName &o2)
{
    assert(s1.size() == s2.size());
    
    FQWriter w1(o1);
    FQWriter w2(o2);
    
    for (auto i = 0; i < s1.size(); i++)
    {
        if (Reader::valid(s1[i]) && Reader::valid(s2[i]))
        {
            ParserFQ::parse(Reader(s1[i]), Reader(s2[i]), [&](const ParserFQ::Data &x)
            {
                w1.write(x.name1, x.seq1, x.qual1);
                w2.write(x.name2, x.seq2, x.qual2);
            });
        }
    }
    
    w1.close();
    w2.close();
}

bool Anaquin::headFQGZ(const FileName &file)
{
    Reader r(file, true);
    
    for (auto i = 0; i < 16; i++)
    {
        auto onlyASCII = [&](const std::string &x)
        {
            for (auto c: x)
            {
                if (static_cast<char>(c) < 0 || static_cast<char>(c) > 127)
                {
                    return false;
                }
            }
            return true;
        };
        
        std::string line;
        
        if (r.nextLine(line))
        {
            if ((i % 4) == 0 && (!onlyASCII(line) || line[0] != '@'))
            {
                return false;
            }
        }
    }
    
    return true;
}

void Anaquin::mv(const FileName &src, const FileName &dst)
{
    const auto cmd = "mv " + src + " " + dst;
    
    if (system(cmd.c_str()))
    {
        throw FailedCommandException("Failed: " + cmd);
    }
}

CopyGZStatus Anaquin::copyGZ(const FileName &src_, std::shared_ptr<std::ofstream> dst)
{
    assert(dst);
    auto src = src_;
    
    if (src.empty() || !Reader::valid(src)) { return CopyGZStatus::Success; }
    
    // Need to autocorrect?
    auto hack = false;
    
    if (!headFQGZ(src))
    {
        const auto tmp = hackUZero(src);
        
        if (headFQGZ(tmp))
        {
            hack = true;
            
            // This is the new source...
            src = hack;
        }
        else
        {
            // Even correction failed...
            return CopyGZStatus::Failed;
        }
    }
    
    std::ifstream r(src, std::ios::in | std::ios::binary);
    
    const static int BUF_SIZE = 4096;
    char buf[BUF_SIZE];
    
    do {
        r.read(&buf[0], BUF_SIZE);
        dst->write(&buf[0], r.gcount());
    } while (r.gcount() > 0);
    
    r.close();
    
    return hack ? CopyGZStatus::Corrected : CopyGZStatus::Success;
};

std::string Anaquin::date()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    
    strftime(buffer, 80, "%d-%m-%Y %H:%M:%S", timeinfo);
    std::string str(buffer);
    
    return str;
}

std::string Anaquin::compressGZ(const std::string &x)
{
    return gzip::compress(x.data(), x.size());
}

void Anaquin::copy(const FileName &x, const FileName &y)
{
    if (x != y)
    {
        std::ifstream src(x);
        std::ofstream dst(y);
        dst << src.rdbuf();
    }
}

void Anaquin::runCmd(const Scripts &cmd)
{
    if (system(cmd.c_str()))
    {
        throw FailedCommandException("Failed: " + cmd);
    }
}

FileName Anaquin::hackUZero(const FileName &file)
{
    const auto file2 = file + "_";
    
    std::ofstream w;
    w.open(file2, std::ios::binary | std::ios::out);
    
    std::ifstream r(file, std::ios::in | std::ios::binary);
    bool hack = true;
    
    char c;
    while(r.read(&c, sizeof(char)))
    {
        if (hack && c == 0) { continue; }
        hack = false;
        w.write(&c, sizeof(c));
    }
    
    w.close();
    return file2;
}

int Anaquin::parseChrID(const ChrID &x)
{
    const auto t = remove(x, "chr");
    
    // Eg: chr1, chr2 ...
    if (isNumber(t))
    {
        return stoi(t);
    }
    else
    {
        std::hash<std::string> hasher;
        auto hashed = hasher(t);
        return 22 + abs((int) hashed);
    }
}

std::vector<FileName> Anaquin::listFiles(const Path &path, const std::string &r, const std::string &ext)
{
    std::vector<FileName> files;
    struct dirent *dp;
    auto dirp = opendir(path.c_str());
    while ((dp = readdir(dirp))) { files.push_back(dp->d_name); }
    
    files.erase(std::remove_if(files.begin(), files.end(),
                               [&](const FileName& x) { return !isSubstr(x, r) || !isSubstr(x, ext); }), files.end());
    
    return files;
}

void Anaquin::runScript(const Scripts &script, const std::string &args)
{
    const auto f = [&](const std::string &cmd)
    {
        if (system(cmd.c_str()))
        {
            throw FailedCommandException("Failed: " + cmd);
        }
    };
    
    // Create a copy of the script
    const auto tmp = tmpFile();
    
    std::ofstream out(tmp);
    out << script;
    out.close();
    
    // Run the script with given arguments
    f("python " + tmp + " " + args);
}

FileName Anaquin::script2File(const Scripts &x)
{
    const auto tmp = tmpFile();
    std::ofstream out(tmp);
    out << x;
    out.close();
    return tmp;
}

FileName Anaquin::tmpFile()
{
    static std::mutex mtx;
    mtx.lock();
    
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    auto f = std::string(tmpnam(NULL));
#pragma clang diagnostic pop
    
    auto random = [&](std::size_t len)
    {
        auto randchar = []() -> char
        {
            const char charset[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
            const size_t max_index = (sizeof(charset) - 1);
            return charset[rand() % max_index];
        };
        
        std::string str(len, 0);
        std::generate_n(str.begin(), len, randchar);
        
        return str;
    };
    
    // Avoid filename crashing...
    f = f + random(3) + random(3);
    
    mtx.unlock();
    
    return f;
}

Path Anaquin::execPath()
{
    int length = wai_getExecutablePath(NULL, 0, NULL);
    char *path = (char*)malloc(length + 1);
    int dirname_length;
    wai_getExecutablePath(path, length, &dirname_length);
    path[length] = '\0';
    
    return dirname(path);
}

bool Anaquin::exists(const FileName &file)
{
    struct stat buffer;
    return (stat(file.c_str(), &buffer) == 0);
}

std::string Anaquin::readFile(const FileName &file)
{
    std::ifstream x(file);
    return std::string((std::istreambuf_iterator<char>(x)), std::istreambuf_iterator<char>());
}

std::vector<double> Anaquin::rank(const std::vector<double> &x)
{
    auto tmp1 = x;
    std::sort(tmp1.begin(), tmp1.end());
    auto tmp2 = std::vector<double>();
    
    for (auto i = 0; i < tmp1.size(); i++)
    {
        auto j = i;
        
        // What's the last element equal?
        for (; (j+1) < tmp1.size() && (tmp1[j] == tmp1[j+1]); j++) {}
        
        // Delta
        const auto n = 1+(j-i);
        
        // Sum of the ranks (startin from 1)
        auto sum = 0.0;
        
        for (auto k = 0; k < n; k++) { sum += (i+1+k); }
        
        // Average the ranks
        const auto m = sum / n;
        
        for (auto k = 0; k < n; k++)
        {
            tmp2.push_back(m);
        }
        
        i += (j-i);
    }
    
    std::map<double, double> tmp3;
    for (auto i = 0; i < tmp1.size(); i++)
    {
        tmp3[tmp1[i]] = tmp2[i];
    }
    
    std::vector<double> r;
    for (auto i = 0; i < x.size(); i++) {
        r.push_back(tmp3.at(x[i]));
    }
    
    return r;
}

bool Anaquin::isFloat(const std::string &x)
{
    std::istringstream iss(x);
    float f;
    iss >> std::noskipws >> f;
    return iss.eof() && !iss.fail();
}

void Anaquin::rm(const FileName &file)
{
    if (exists(file))
    {
        std::remove(file.c_str());
    }
}

std::string Anaquin::run(const Command &x)
{
    std::string data;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    
    auto cmd = x;
    cmd.append(" 2>&1");
    
    stream = popen(cmd.c_str(), "r");
    if (stream)
    {
        while (!feof(stream)) { if (fgets(buffer, max_buffer, stream)) data.append(buffer); }
        pclose(stream);
    }
    
    return data;
}
