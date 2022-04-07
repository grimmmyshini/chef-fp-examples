// Clad header
#include "clad/Differentiator/Differentiator.h"

#include <iostream> // necessary for 'cout'
#include <cmath> // necessary for 'sin' and other math functions
#include <iomanip> // necessary for setprecision

// defines sin(pi * x)
template <typename prec> 
prec f(prec x) {
 prec pi = M_PI;
 prec tmp = x * pi;
 prec tmp2 = sin(tmp);
 return tmp2;
}

// defines the simpsons rule for integral estimation.
template <typename prec>
prec simpsons(prec a, prec b) {
 int n = 1000000;
 prec h = (b - a) / (2.0 * n);
 prec x = a;
 prec tmp;
 prec fa = f(a), fb = f(b);
 prec s1 = fa + fb;

 for(int l = 0; l < n; l++) {
   x = x + h;
   s1 = s1 + 4.0 * f(x);
   x = x + h;
   s1 = s1 + 2.0 * f(x);
 }

 s1 = s1 - fb;
 tmp = h * M_PI / 3.0;
 s1 = s1 * tmp;
 return s1;
}

int main() {
  long double a = 0, b = 1;
  // Actual value = 2 for integral of sin(pi * x)/pi over [a, b]
  long double result_ldbl = simpsons(a, b);
  double result_dbl = simpsons<double>(a, b);
  float result_flt = simpsons<float>(a, b);

  std::cout << std::setprecision(15)
            << "Simpsons with long double: " << result_ldbl << " error: " << std::abs(2.0 - result_ldbl) << "\n"
            << "Simspons with double: " << result_dbl << " error: " << std::abs(2.0 - result_dbl) << "\n"
            << "Simspons with float: " << result_flt << " error: " << std::abs(2.0 - result_flt);
}
