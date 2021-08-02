#ifndef DETECTOR_HPP
#define DETECTOR_HPP

#include "data/data.hpp"

namespace Anaquin
{
    struct Detector
    {
        static Build fromBAM(const FileName &);
        static Build fromVCF(const FileName &);
    };
}

#endif
