#include <iomanip>  // necessary for setprecision
#include <iostream> // necessary for 'cout'

#include "simpsons.hpp"

// Clad header
#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#define ITERATIONS 1000

#include "Derivative.hpp"

int main() {
  // clad::estimate_error(simpsons<double>);
  double a = 0, b = 1, d_a = 0, d_b = 0, final_error = 0;
  
  double result_dbl = simpsons<double>(a, b, ITERATIONS);

  clad::simpsons_grad(a, b, &d_a, &d_b, final_error);

  clad::printErrorReport();
}
