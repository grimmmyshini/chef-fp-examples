/* 
 * This benchmark outputs all errors calculated by clad.
 * Both clad and adapt output to a stringstream buffer.
 * 
 * The outputs are redirected to a stringstream buffer so as to separate
 * them from the benchmark results.
 * 
 * Example adapted from ADAPT-FP: https://github.com/LLNL/adapt-fp
 *
 * Typical invocation of the benchmark looks like this:
 * clang -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang /path/to/clad.so 
 * -Ipath/to/clad/include -Ipath/to/adapt-fp -x c++ -lstdc++ -O2
 * -lm -isystem path/to/Google/benchmark/include  -Lpath/to/benchmark/build/src
 * ./arclen_full_buf.cpp -lbenchmark -lpthread -DCODI_ZeroAdjointReverse=0
 * 
 */

#include "benchmark/benchmark.h"

#include <cmath>
#include <fstream> // Necessary for saving error estimates.
#include <iostream>
#include <stdio.h>
#include <math.h>

#include "adapt.h"
#include "adapt-impl.cpp"

#include "clad/Differentiator/Differentiator.h"

/*
 * Precimonious output labels located in comments beside variable declarations
 */

#define ABS(x) ( ((x) < 0.0) ? (-(x)) : (x) )

#define PI  3.1415926535897932L
#define ANS 5.795776322412856L
#define EPS 1e-10

#define N 1000000

AD_real h  = PI / (double)N;     // double

AD_real t1 = 0.0;
AD_real t2;
AD_real t3;                      // double

AD_real s1 = 0.0;                // double

AD_real d1 = 1.0;
AD_real d2;                      // float

AD_real fun (AD_real x)
{
    d2 = d1;    // also d1 in original
    AD_INTERMEDIATE(d2, "d2");
    t3 = x;     // also t1 in original
    AD_INTERMEDIATE(t3, "t3");

    int k;
    for (k = 1; k <= 5; k+=1)
    {
        d2 = 2.0 * d2;
        AD_INTERMEDIATE(d2, "d2");
        t3 = t3 + sin (d2 * x) / d2;
        AD_INTERMEDIATE(t3, "t3");
    }
    return t3;
}

void do_fun ()
{
    int i;
    for (i = 1; i <= N; i+=1)
    {
        t2 = fun (i * h);
        AD_INTERMEDIATE(t2, "t2");
        s1 = s1 + sqrt (h * h + (t2 - t1) * (t2 - t1));
        AD_INTERMEDIATE(s1, "s1");
        t1 = t2;
        AD_INTERMEDIATE(t1, "t1");
    }
}

struct cout_redirect {
    cout_redirect( std::streambuf * new_buffer ) 
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

private:
    std::streambuf * old;
};

static void ErrorEstimateArcLenAdapt(benchmark::State& state) {
    std::stringstream buffer;
    cout_redirect output(buffer.rdbuf());
    for (auto _ : state) {
        AD_begin();
        AD_INTERMEDIATE(h, "h");
        AD_INTERMEDIATE(t1, "t1");
        AD_INTERMEDIATE(s1, "s1");
        AD_INTERMEDIATE(d1, "d1");

        do_fun();

        AD_DEPENDENT(s1, "s1", EPS);
        AD_report();
    }
}

double fun_clad(double x) {
  double d2 = 1.0; // also d1 in original
  double t3 = x;   // also t1 in original

  int k;
  for (k = 1; k <= 5; k += 1) {
    d2 = 2.0 * d2;
    t3 = t3 + sin(d2 * x) / d2;
  }
  return t3;
}

double do_fun_clad(double s1, double t1) {
  double t2, h = PI / (double)N;

  for (int i = 1; i <= N; i += 1) {
    t2 = fun_clad(i * h);
    s1 = s1 + sqrt(h * h + (t2 - t1) * (t2 - t1));
    t1 = t2;
  }

  return s1;
}

void do_fun_clad_grad(double, double, clad::array_ref<double>, clad::array_ref<double>, double &, std::ostream &);

static void ErrorEstimateArcLenClad(benchmark::State& state) {
  // double s1 = do_fun_clad(0.0, 0.0);
  auto df = clad::estimate_error<true>(do_fun_clad);
  // df.dump();
  // Enable before using std::cout for measurements.
  // std::ios_base::sync_with_stdio(false);
  std::stringstream buffer;
  cout_redirect output(buffer.rdbuf());
  for (auto _ : state) {
    double ds1 = 0, dt1 = 0;
    double fin_error = 0;
    do_fun_clad_grad(0, 0, &ds1, &dt1, fin_error, std::cout);
  }
}

BENCHMARK(ErrorEstimateArcLenAdapt)->Unit(benchmark::kSecond)->Iterations(10);

BENCHMARK(ErrorEstimateArcLenClad)->Unit(benchmark::kSecond)->Iterations(10);

// Define our main
BENCHMARK_MAIN();
