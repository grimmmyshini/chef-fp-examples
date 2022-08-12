#ifndef SIMPSONS_ADAPT_HPP
#define SIMPSONS_ADAPT_HPP

#include <cmath>

#include "adapt.h"
#include "adapt-impl.cpp"

namespace adapt
{
  // defines the simpsons rule for integral estimation.
  AD_real simpsons(AD_real a, AD_real b)
  {
    int n = ITERATIONS;

    AD_real pi = M_PI;

    AD_real h = (b - a) / (2.0 * n);
    AD_INTERMEDIATE(h, "h");
    AD_real x = a;
    AD_INTERMEDIATE(x, "x");

    AD_real tmp, x_pi, sin_x_pi;
    x_pi = a * pi;
    AD_INTERMEDIATE(x_pi, "x_pi");
    sin_x_pi = sin(x_pi);
    AD_INTERMEDIATE(sin_x_pi, "sin_x_pi");
    AD_real fa = sin_x_pi;

    x_pi = b * pi;
    AD_INTERMEDIATE(x_pi, "x_pi");
    sin_x_pi = sin(x_pi);
    AD_INTERMEDIATE(sin_x_pi, "sin_x_pi");
    AD_real fb = sin_x_pi;

    AD_real s1 = fa + fb;
    AD_INTERMEDIATE(s1, "s1");

    for (int l = 0; l < n; l++)
    {
      x = x + h;
      AD_INTERMEDIATE(x, "x");
      x_pi = x * pi;
      AD_INTERMEDIATE(x_pi, "x_pi");
      sin_x_pi = sin(x_pi);
      AD_INTERMEDIATE(sin_x_pi, "sin_x_pi");
      s1 = s1 + 4.0 * sin_x_pi;
      AD_INTERMEDIATE(s1, "s1");
      x = x + h;
      AD_INTERMEDIATE(x, "x");
      x_pi = x * pi;
      AD_INTERMEDIATE(x_pi, "x_pi");
      sin_x_pi = sin(x_pi);
      AD_INTERMEDIATE(sin_x_pi, "sin_x_pi");
      s1 = s1 + 2.0 * sin_x_pi;
      AD_INTERMEDIATE(s1, "s1");
    }

    s1 = s1 - fb;
    AD_INTERMEDIATE(s1, "s1");
    tmp = h * M_PI / 3.0;
    AD_INTERMEDIATE(tmp, "tmp");
    s1 = s1 * tmp;
    AD_INTERMEDIATE(s1, "s1");
    return s1;
  }
}

#endif