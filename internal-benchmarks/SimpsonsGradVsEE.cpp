#include "benchmark/benchmark.h"

// FIXME: If we move this before benchmark.h we have tons of errors due to a bug
#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include <cmath>

// Simpsons function, approximates the integral of sin(x * pi) over the interval
// [a,b] for n iterations.
long double simpsons(long double a, long double b, int n) {
  long double h = (b - a) / (2.0 * n);
  long double x = a;
  long double tmp;
  long double pi = M_PI;
  long double fa = pi * std::sin(pi * a), fb = pi * std::sin(pi * b);
  long double s1 = fa;

  for (int l = 0; l < n; l++) {
    x = x + h;
    s1 = s1 + 4.0 * pi * std::sin(pi * x);
    x = x + h;
    s1 = s1 + 2.0 * pi * std::sin(pi * x);
  }

  s1 = s1 - fb;
  tmp = h / 3.0;
  s1 = s1 * tmp;
  return s1;
}

// Error estimation
static void ErrorEstimateCalc(benchmark::State& state) {
  auto df = clad::estimate_error(simpsons);
  int n = state.range(0);
  for (auto _ : state) {
    long double diffa = 0, diffb = 0, diffn = 0;
    double finError = 0;
    df.execute(0, 1, n, &diffa, &diffb, &diffn, finError);
  }
}
BENCHMARK(ErrorEstimateCalc)->RangeMultiplier(10)->Range(1E2, 1E6);

// Gradient generation
static void GradientCalc(benchmark::State& state) {
  auto df = clad::gradient(simpsons);
  int n = state.range(0);
  for (auto _ : state) {
    long double diffa = 0, diffb = 0, diffn = 0;
    df.execute(0, 1, n, &diffa, &diffb, &diffn);
  }
}
BENCHMARK(GradientCalc)->RangeMultiplier(10)->Range(1E2, 1E6);

// Define our main.
BENCHMARK_MAIN();
