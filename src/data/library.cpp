#include <vector>
#include "stats/stats.hpp"
#include "tools/tools.hpp"
#include "data/library.hpp"
#include "stats/ss/stats.hpp"

using namespace Anaquin;

struct Library::Impl
{
    std::vector<std::string> hs;
    std::vector<Base> ls; // Read lengths
};

Library::Library()
{
    _impl = std::shared_ptr<Library::Impl>(new Library::Impl());
}

unsigned Library::heads() const { return _impl->hs.size(); }

Base Library::meanRL() const
{
    return SS::mean(_impl->ls);
}

std::string Library::inst(Format f) const
{
    switch (f)
    {
        case Illumina_V2: { return tokN(_impl->hs.front(), ":", 0); }
        default : { return MISSING; }
    }
}

std::string Library::run(Format f) const
{
    switch (f)
    {
        case Illumina_V2: { return tokN(_impl->hs.front(), ":", 1); }
        default : { return MISSING; }
    }
}

std::string Library::flow(Format f) const
{
    switch (f)
    {
        case Illumina_V2: { return tokN(_impl->hs.front(), ":", 2); }
        default : { return MISSING; }
    }
}

std::string Library::lane(Format f) const
{
    switch (f)
    {
        case Illumina_V2: { return tokN(_impl->hs.front(), ":", 3); }
        default : { return MISSING; }
    }
}

Library::Format Library::format() const
{
    if (_impl->hs.empty())
    {
        return Format::None;
    }

    auto tmp1 = std::vector<std::string>();
    split(_impl->hs.front(), " ", tmp1);

    auto tmp2 = std::vector<std::string>();
    split(tmp1.front(), ":", tmp2);
    
    if (tmp2.size() == 7)
    {
        return Format::Illumina_V2;
    }

    return Format::None;
}

void Library::addInfo(const std::string &head, const Sequence &seq)
{
    _impl->hs.push_back(head);
    _impl->ls.push_back(seq.size());
}
