#pragma once

#include <array>
#include <string>
#include <iosfwd>       // forward-declares std::ostream without pulling in <ostream>

struct UniformityResult
{
   std::array<int, 10> bins{};
   double chi_square = 0.0;
   double uniformity_p_value = 0.0;

   std::string toString () const;
};

std::ostream &operator<< (std::ostream &os, const UniformityResult &r);

