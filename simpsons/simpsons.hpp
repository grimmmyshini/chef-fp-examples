#ifndef SIMPSONS_HPP
#define SIMPSONS_HPP

#include <cmath>

// defines the simpsons rule for integral estimation of sin( x * pi ) in [a,b]
template <typename highprec, typename lowPrec = highprec>
highprec simpsons(lowPrec a, lowPrec b, int n) {
  highprec pi = M_PI;
  lowPrec h = (b - a) / (2.0 * n);
  highprec x = a;
  lowPrec tmp, x_pi, sin_x_pi;
  x_pi = a * pi;
  sin_x_pi = pi * sin(x_pi);
  lowPrec fa = sin_x_pi;
  
  x_pi = b * pi;
  sin_x_pi = pi * sin(x_pi);
  lowPrec fb = sin_x_pi;
  highprec s1 = fa + fb;

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
