#include <iomanip>  // necessary for setprecision
#include <iostream> // necessary for 'cout'

#include "simpsons.hpp"

int main() {
  long double a = 0, b = 1;
  // Actual value = 2 for integral of sin(pi * x)/pi over [a, b]
  long double result_ldbl = simpsons(a, b);
  double result_dbl = simpsons<double>(a, b);
  float result_flt = simpsons<float>(a, b);

  std::cout << std::setprecision(15)
            << "Simpsons with long double: " << result_ldbl
            << " error: " << std::abs(2.0 - result_ldbl) << "\n"
            << "Simspons with double: " << result_dbl
            << " error: " << std::abs(2.0 - result_dbl) << "\n"
            << "Simspons with float: " << result_flt
            << " error: " << std::abs(2.0 - result_flt);
}
