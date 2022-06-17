#include "benchmark/benchmark.h"

#include "clad/Differentiator/Differentiator.h"

#include <cmath>
#include <fstream> // Necessary for saving error estimates.
#include <iostream>
#include <math.h>
#include <stdio.h>

/*
 * Example adapted from ADAPT-FP: https://github.com/LLNL/adapt-fp
 * Precimonious output labels located in comments beside variable declarations.
 *
 * Typical invocation of the benchmark looks like this:
 * clang -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang /path/to/clad.so -Ipath/to/clad/include -x c++ -lstdc++ 4
 * -lm -isystem path/to/Google/benchmark/include  -Lpath/to/benchmark/build/src -O2 ./arclenBenchmark.cpp -lbenchmark -lpthread
 * 
 */

#define ABS(x) (((x) < 0.0) ? (-(x)) : (x))

#define PI 3.1415926535897932L
#define ANS 5.795776322412856L
#define EPS 1e-10

#define N 1000000

double fun(double x) {
  double d2 = 1.0; // also d1 in original
  double t3 = x;   // also t1 in original

  int k;
  for (k = 1; k <= 5; k += 1) {
    d2 = 2.0 * d2;
    t3 = t3 + sin(d2 * x) / d2;
  }
  return t3;
}

double do_fun(double s1, double t1) {
  double t2, h = PI / (double)N;

  for (int i = 1; i <= N; i += 1) {
    t2 = fun(i * h);
    s1 = s1 + sqrt(h * h + (t2 - t1) * (t2 - t1));
    t1 = t2;
  }

  return s1;
}

static void ErrorEstimateArcLen(benchmark::State& state) {
  double s1 = do_fun(0.0, 0.0);
  auto df = clad::estimate_error<true>(do_fun);
  // Enable before using std::cout for measurements.
  // std::ios_base::sync_with_stdio(false); 
  std::ofstream fil("err_file.txt");
  for (auto _ : state) {
    double ds1 = 0, dt1 = 0;
    double fin_error = 0;
    df.execute(0, 0, &ds1, &dt1, fin_error, fil);
  }
}

BENCHMARK(ErrorEstimateArcLen)->Unit(benchmark::kSecond);

// Define our main.
BENCHMARK_MAIN();
