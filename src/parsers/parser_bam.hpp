#ifndef PARSER_BAM_HPP
#define PARSER_BAM_HPP

#include <memory>
#include "data/alignment.hpp"
#include "data/analyzer.hpp"
#include "parsers/parser_bam.hpp"

namespace Anaquin
{
    struct ParserBAM
    {
        struct Info
        {
            Progress p = 0;

            // Whether this is a multi-alignment
            bool multi;
            
            // Whether there is insertion
            bool ins;
            
            // Whether there is deletion
            bool del;
            
            // Whether there is skipped region
            bool skip;
            
            // Whether there is clipping
            bool clip;
            
            // Size of the chromosome of the alignment
            Base length;
            
            void *b;
            void *h;
        };

        class Data : public Alignment
        {
            friend struct ParserBAM;
            
            public:
            
                /*
                 * Optional fields
                 */

                void lSeq();     // Lazy loadding of sequence
                void lQual();    // Lazy loading of quality
                void lName();    // Lazy loading of reads
                void lMateID();  // Lazy loading of mate ID
                void lMatePos(); // Lazy loading of mate position
                void lCigar();

                inline void *b() const { return _b; } // bam1_t
                inline void *h() const { return _h; } // bam_hdr_t

                void *copyH() const;
                void *copyB() const;
            
            private:
            
                Base getEnd();

                void *_b;
                void *_h;
        };
        
        typedef std::function<void (Data &, const Info &)> Functor;
        
        static std::map<ChrID, Base> header(const FileName &);
        
        /*
         * Parse alignment file and call a function for each read. The locus of the read is determined
         * by all cigars. Use forCigar() if cigar-level resolution is necessary.
         */

        static void parse(const FileName &, Functor);
        
        struct ParseResult
        {
            // Sequin read?
            bool isSeq = false;
            
            /*
             * Should we reverse complement this read sequence and
             * also reverse the quality scores?
             *
             *     std::reverse(x.qual.begin(), x.qual.end());
             *     x.seq = revcomp(x.seq);
             */
            
            bool rc = false;
        };
        
        // Functor for parsing alignments selectively
        typedef std::function<ParseResult (Data &, const Info &)> F2;

        static void parse_(const FileName &, AlignedReads &r, F2);
    };
    
    /*
     * Define a template for looping alignment cigars. The bases are all in 0-based, like all
     * other programming. The caller should convert them to 1-based by an increment.
     */
    
    typedef std::function<void (Base i, Base j, char, char)> F1;
    typedef std::function<void (Base i, Base j, Base l)> F2;

    void forCigar(const Alignment &, const Sequence &, F1, F1, F2, F2, F2);
}

#endif
