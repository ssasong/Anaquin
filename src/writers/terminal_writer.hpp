#ifndef TERMINAL_WRITER_HPP
#define TERMINAL_WRITER_HPP

#include <iostream>
#include "writers/writer.hpp"

namespace Anaquin
{
    class TerminalWriter : public Writer<>
    {
        public:

            inline void close() override {}

            inline void open(const FileName &) override {}

            inline void write(const std::string &str, bool newLine = true) override
            {
                std::cout << str << std::endl;
            }
    };
}

#endif
