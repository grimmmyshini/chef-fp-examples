#include "benchmark/benchmark.h"

#include "arclen.hpp"

namespace lower_prec
{
  double do_fun(double s1, double t1)
  {
    int i, k;
    double t2, h = PI / ITERATIONS, x, t3;
    float d2;

        for (i = 1; i <= ITERATIONS; i += 1)
    {
      // t2 = fun_clad(i * h);
      x = i * h;
      d2 = 1.0; // also d1 in original
      t3 = x;   // also t1 in original

      for (k = 1; k <= 5; k += 1)
      {
        d2 = 2.0 * d2;
        t3 = t3 + sin(d2 * x) / d2;
      }

      t2 = t3;

      s1 = s1 + sqrt(h * h + (t2 - t1) * (t2 - t1));
      t1 = t2;
    }

    return s1;
  }
}

static void ArcLenLowerPrec(benchmark::State &state)
{
  double result;

  for (auto _ : state)
  {
    result = lower_prec::do_fun(0, 0);

    benchmark::DoNotOptimize(result);
  }
}

static void ArcLenHighPrec(benchmark::State &state)
{
  double result;

  for (auto _ : state)
  {
    result = do_fun(0, 0);

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(ArcLenLowerPrec)->Unit(benchmark::kSecond);
BENCHMARK(ArcLenHighPrec)->Unit(benchmark::kSecond);

// Define our main
BENCHMARK_MAIN();