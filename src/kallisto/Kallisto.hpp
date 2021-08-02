#ifndef KALLISTO_HPP
#define KALLISTO_HPP

#include <map>
#include <memory>
#include "tools/tools.hpp"
#include "data/library.hpp"

namespace Anaquin
{
    typedef std::map<SequinID, std::map<Kmer, Count>> SequinKM;

    // Default k-mer length
    const KmerLen K_DEFAULT_L = 23;
    
    // Default classification rule
    const Proportion K_DEFAULT_R  = 0.20;

    /*
     * Statistics from Kallisto
     */
    
    struct KStats
    {
        Path work;
        
        // Count for the categories
        std::map<Bin, Count> c1, c2;
        
        // Output file names (f2 might not be defined)
        std::map<Bin, std::vector<FileName>> f1, f2;
 
        // Total number of reads
        inline Count total() const { return sum(c1); }

        // Number of reads for a bin
        inline Count binN(Bin x) const { return c1.count(x) ? c1.at(x) : 0; }

        // Percentage of reads for a bin
        inline Proportion binP(Bin x) const { return c1.count(x) ? ((Proportion) c1.at(x) / total()) : 0.0; }

        Library lib;
        std::set<StandardID> stds, seqs;
        std::map<StandardID, SequinID> rSeqs;
        std::map<StandardID, std::set<SequinID>> aSeqs;
        
        // Measured Count for unique k-mers
        SequinKM uniqs;
        
        // Measured Count for shared k-mers
        SequinKM shared;
        
        // Number of reads matching this index (e.g. sequins)
        Count nMatch = 0;
        
        // Number of reads not matching this index (e.g. genome)
        Count nNMatch = 0;

        // Number of reads for each sequin (2x for paired-end)
        CustomMap<SequinID, Count> sqc;
        
        // Information about reads
        std::map<Bin, std::map<ReadName, SequinID>> r1, r2;
    };
    
    struct KOptions
    {
        KOptions() : flipBefore(false), flip(true), onlySeqLad(false) {}
        
        inline bool writeBAM() const
        {
            return false;
        }
        
        bool writeReads = false;
        
        // Writing detailed statistics?
        FileName tsv;
        
        FileName index;
        
        // Flip before matching?
        bool flipBefore;
        
        // Flip sequin reads?
        bool flip;
        
        // K-mer length
        KmerLen k;
        
        // How to run Kallisto?
        Product prod;
        
        // How many k-mers to skip?
        Count skipKM;
        
        // Classification rule (minimum percentage matching)
        Proportion rule;
        
        // Number of threads
        Count thr;
        
        // Only "sequin" and "ladder" bins?
        bool onlySeqLad;
    };
}

#endif
