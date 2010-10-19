#include "Util.hh"
#include <iostream>

namespace BASim {

void _error(const char* file, const char* function,
            int line, const char* message)
{
    std::cerr << file << ':' << function << ':' << line << ": " << message << std::endl;
    exit(1);
}

} // namespace BASim
