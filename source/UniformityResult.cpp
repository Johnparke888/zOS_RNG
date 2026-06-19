#include "UniformityResult.h"

#include <sstream>
#include <iomanip>
#include <ostream>

std::string UniformityResult::toString () const
{
   std::ostringstream os;
   os << "bins=[";
   for (std::size_t i = 0; i < bins.size (); ++i)
   {
      if (i != 0)
      {
         os << ", ";
      }
      os << bins[i];
   }
   os << "], chi_square=" << std::fixed << std::setprecision (4) << chi_square << ", pvalue=" << uniformity_p_value << "";
   return os.str ();
}

std::ostream &operator<< (std::ostream &os, const UniformityResult &r)
{
   return os << r.toString ();
}
