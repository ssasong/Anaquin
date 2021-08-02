#include <htslib/sam.h>
#include "tools/samtools.hpp"
#include "parsers/parser_bam.hpp"
#include <boost/algorithm/string/predicate.hpp>

using namespace Anaquin;

void Anaquin::forCigar(const Alignment &x, const Sequence &seq,
                       F1 M, F1 S, F2 I, F2 D, F2 SC)
{
    // Relative to the query sequence
    Base i = 0;
    
    // Relative to the decoy chromosome (x.l.start is 1-based)
    Base j = x.l.start-1;
    
    for (auto ci: x.cigars)
    {
        const auto val = ci.second;
     
        switch (ci.first)
        {
            case BAM_CMATCH:
            {
                // For each position...
                for (auto k = 0; k < val; k++)
                {
                    assert(j < seq.size());
                    
                    const char r = toupper(seq.at(j));
                    const char q = toupper(x.seq[i]);

                    if (seq.at(j) == x.seq[i])
                    {
                        M(i, j, r, q); // j is 0-based (it needs to be converted to 1-based)
                    }
                    else
                    {
                        assert(i < x.seq.size());
                        
                        if (r != 'N' && q != 'N')
                        {
                            assert(r == 'A' || r == 'C' || r == 'T' || r == 'G');
                            assert(q == 'A' || q == 'C' || q == 'T' || q == 'G');
                            S(i, j, r, q); // j is 0-based (it needs to be converted to 1-based)
                        }
                    }
                    
                    i++;
                    j++;
                }
                
                break;
            }
                
            case BAM_CINS: { I(i, j, val); i += val; break; }
            case BAM_CDEL: { D(i, j, val); j += val; break; }
                
            case BAM_CREF_SKIP:  { SC(i, j, val); j += val; break; }
            case BAM_CSOFT_CLIP: { SC(i, j, val); i += val; break; }
            case BAM_CHARD_CLIP: { SC(i, j, val); break; }
                
            case BAM_CPAD:
            case BAM_CDIFF:
            case BAM_CBACK:
            case BAM_CEQUAL: { assert(false); break; }
                
            default:
            {
                throw std::runtime_error("Unknown " + std::to_string(ci.first));
            }
        }
    }
};

void ParserBAM::Data::lMateID()
{
    mID = std::string(static_cast<bam_hdr_t *>(_h)->target_name[static_cast<bam1_t *>(_b)->core.mtid]);
}

void ParserBAM::Data::lMatePos()
{
    mPos = static_cast<bam1_t *>(_b)->core.mpos;
}

void ParserBAM::Data::lName()
{
    name = bam_get_qname(static_cast<bam1_t *>(_b));
}

void ParserBAM::Data::lSeq()
{
    seq = bam2seq(static_cast<bam1_t *>(_b));
}

void ParserBAM::Data::lQual()
{
    qual = bam2qual(static_cast<bam1_t *>(_b));
}

void * ParserBAM::Data::copyH() const
{
    assert(_h);
    return bam_hdr_dup(static_cast<bam_hdr_t *>(_h));
}

void * ParserBAM::Data::copyB() const
{
    assert(_b);
    return bam_dup1(static_cast<bam1_t *>(_b));
}

void ParserBAM::Data::lCigar()
{
    const auto x = static_cast<bam1_t *>(_b);
    const auto cig = bam_get_cigar(x);
    cigars.clear();
    
    for (auto i = 0u; i < x->core.n_cigar; i++)
    {
        cigars.push_back(std::pair<Cigar, Base>(bam_cigar_op(cig[i]), bam_cigar_oplen(cig[i])));
    }
}

Base ParserBAM::Data::getEnd()
{
    assert(_h && _b);
    
    auto t = static_cast<bam1_t *>(_b);
    const auto cig = bam_get_cigar(t);
    auto p = t->core.pos;
    
    for (auto i = 0; i < t->core.n_cigar; i++)
    {
        const auto ol = bam_cigar_oplen(cig[i]);
        
        switch (bam_cigar_op(cig[i]))
        {
            case BAM_CDEL:   { p += ol; break; }
            case BAM_CMATCH: { p += ol; break; }
            case BAM_CINS:
            case BAM_CPAD:
            case BAM_CBACK:
            case BAM_CDIFF:
            case BAM_CEQUAL:
            case BAM_CREF_SKIP:
            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP: { break; }
        }
    }

    return p;
}

std::map<ChrID, Base> ParserBAM::header(const FileName &file)
{
    auto f = sam_open(file.c_str(), "r");
    
    if (!f)
    {
        throw std::runtime_error("Failed to open: " + file);
    }
    
    auto h = sam_hdr_read(f);

    std::map<ChrID, Base> c2b;
    
    for (auto i = 0; i < h->n_targets; i++)
    {
        c2b[std::string(h->target_name[i])] = h->target_len[i];
    }
    
    return c2b;
}

void ParserBAM::parse_(const FileName &file, AlignedReads &r, F2 f)
{
    r.clear();
    
    ParserBAM::parse(file, [&](Data &x, const Info &i)
    {
        const auto rr = f(x, i);
        
        if (rr.isSeq)
        {
            /*
             * Assume f() lazy initalize name, qual and seq
             */
            
            if (x.isFirstPair)
            {
                r[x.name].first = std::shared_ptr<AlignedRead>(new AlignedRead());
                r[x.name].first->b = x.copyB();
                r[x.name].first->rc = rr.rc;
            }
            else
            {
                r[x.name].second = std::shared_ptr<AlignedRead>(new AlignedRead());
                r[x.name].second->b = x.copyB();
                r[x.name].second->rc = rr.rc;
            }
        }
    });
}

void ParserBAM::parse(const FileName &file, Functor x)
{
    auto f = sam_open(file.c_str(), "r");
    
    if (!f)
    {
        throw std::runtime_error("Failed to open: " + file);
    }

    auto t = bam_init1();
    auto h = sam_hdr_read(f);

    Info info;
    Data align;

    while (sam_read1(f, h, t) >= 0)
    {
        info.length = h->target_len[t->core.tid];
        align.mapped = false;
        
        info.b = t;
        info.h = h;

        align._b = t;
        align._h = h;

        align.mapq = t->core.qual;
        align.flag = t->core.flag;
        
        const auto hasCID = t->core.tid >= 0;
        
        #define isPairedEnd(b)    (((b)->core.flag&0x1)   != 0)
        #define isAllAligned(b)   (((b)->core.flag&0x2)   != 0)
        #define isUnmapped(b)     (((b)->core.flag&0x4)   != 0)
        #define isMateUnmapped(b) (((b)->core.flag&0x8)   != 0)
        #define isReversed(b)     (((b)->core.flag&0x10)  != 0)
        #define isMateReversed(b) (((b)->core.flag&0x20)  != 0)
        #define isFirstPair(b)    (((b)->core.flag&0x40)  != 0)
        #define isLastPair(b)     (((b)->core.flag&0x80)  != 0)
        #define isSecondary(b)    (((b)->core.flag&0x100) != 0)
        #define isFailed(b)       (((b)->core.flag&0x200) != 0)
        #define isDuplicate(b)    (((b)->core.flag&0x400) != 0)
        #define isSupplement(b)   (((b)->core.flag&0x800) != 0)
        #define isPrimary(b)      (((b)->core.flag&0x900) == 0)

        align.isPaired      = isPairedEnd(t);
        align.isAllAligned  = isAllAligned(t);
        align.isAligned     = !isUnmapped(t);
        align.isMateAligned = !isMateUnmapped(t);
        align.isReverseC    = isReversed(t);
        align.isMateReverse = isMateReversed(t);
        align.isFirstPair   = isFirstPair(t);
        align.isSecondPair  = isLastPair(t);
        align.isPassed      = !isFailed(t);
        align.isDuplicate   = isDuplicate(t);
        align.isSupplement  = isSupplement(t);
        align.isPrimary     = isPrimary(t);
        align.isSecondary   = isSecondary(t);

        if (hasCID)
        {
            align.cID = std::string(h->target_name[t->core.tid]);
        }
        else
        {
            align.cID = "*";
            align.l.start = 0;
            align.l.end = 0;
        }

        align.tlen = hasCID ? t->core.isize : 0;
        align.mapped = hasCID && !(t->core.flag & BAM_FUNMAP);

        if (align.mapped)
        {
            const auto cigar = bam_get_cigar(t);

            // Is this a multi alignment?
            info.multi = t->core.n_cigar > 1;

            info.ins  = false;
            info.del  = false;
            info.clip = false;
            info.skip = false;
            
            for (auto i = 0u; i < t->core.n_cigar; i++)
            {
                switch (bam_cigar_op(cigar[i]))
                {
                    case BAM_CINS:       { info.ins  = true; break; }
                    case BAM_CDEL:       { info.del  = true; break; }
                    case BAM_CREF_SKIP:  { info.skip = true; break; }
                    case BAM_CSOFT_CLIP:
                    case BAM_CHARD_CLIP:
                    case BAM_CPAD:       { info.clip  = true; break; }
                    default: { break; }
                }
            }
            
            align.l.start = t->core.pos + 1; // 1-based position
            align.l.end   = align.getEnd();
        }
        
        x(align, info);
        info.p++;
    }

    bam_destroy1(t);
    sam_close(f);
}
