#ifndef SIMPSONS_HPP
#define SIMPSONS_HPP

#include <cmath>

// defines the simpsons rule for integral estimation of sin( x * pi ) in [a,b]
template <typename prec>
prec simpsons(prec a, prec b) {
  int n = ITERATIONS;
  prec pi = M_PI;
  prec h = (b - a) / (2.0 * n);
  prec x = a;
  prec tmp, x_pi, sin_x_pi;
  x_pi = a * pi;
  sin_x_pi = sin(x_pi);
  prec fa = sin_x_pi;
  
  x_pi = b * pi;
  sin_x_pi = sin(x_pi);
  prec fb = sin_x_pi;
  prec s1 = fa + fb;

  for (int l = 0; l < n; l++) {
    x = x + h;
    x_pi = x * pi;
    sin_x_pi = sin(x_pi);
    s1 = s1 + 4.0 * sin_x_pi;
    x = x + h;
    x_pi = x * pi;
    sin_x_pi = sin(x_pi);
    s1 = s1 + 2.0 * sin_x_pi;
  }

  s1 = s1 - fb;
  tmp = h * M_PI / 3.0;
  s1 = s1 * tmp;
  return s1;
}

#endif