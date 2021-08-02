#ifndef LIBRARY_HPP
#define LIBRARY_HPP

#include <string>
#include <memory>
#include "data/data.hpp"

namespace Anaquin
{
    struct Library
    {
        enum Format
        {
            None,
            Illumina_V2
        };
        
        struct Impl;
        Library();
        
        // Add a FQ/BAM header info
        void addInfo(const ReadName &, const Sequence &);
        
        // Number of headers already added
        unsigned heads() const;
        
        // Mean read length
        Base meanRL() const;

        Format format() const;
        
        std::string run(Format)  const;
        std::string inst(Format) const;
        std::string flow(Format) const;
        std::string lane(Format) const;

        std::shared_ptr<Impl> _impl;
    };
}

#endif
