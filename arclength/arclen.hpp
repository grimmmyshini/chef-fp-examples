#ifndef ARCLEN_HPP
#define ARCLEN_HPP

#include <cmath>
#include <iostream>
#define ABS(x) (((x) < 0.0) ? (-(x)) : (x))

#define PI 3.1415926535897932L
#define ANS 5.795776322412856L
#define EPS 1e-10

template <typename T, typename L = T> T fun(T x) {
  L d1 = 1.0;
  L d2 = d1; // also d1 in original
  T t3 = x;  // also t1 in original
  int k;
  for (k = 1; k <= 5; k += 1) {
    d2 = 2.0 * d2;
    t3 = t3 + sin(d2 * x) / d2;
  }
  return t3;
}

template <typename T, typename L = T> T do_fun(int iterations) {
  T t2;
  T h = PI/iterations;
  T t1 = 0;
  T s1 = 0;
  int i;
  for (i = 1; i <= iterations; i += 1) {
    t2 = fun(i * h);
    s1 = s1 + sqrt(h * h + (t2 - t1) * (t2 - t1));
    t1 = +t2;
  }
  return s1;
}

#endif
