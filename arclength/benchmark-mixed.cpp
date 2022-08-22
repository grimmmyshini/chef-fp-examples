#include "benchmark/benchmark.h"

int ITERATIONS = 50000000;

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

  double do_fun2(double s1, float t1)
  {
    int i, k;
    double t2, h = PI / ITERATIONS, x, t3;
    double d2;

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
static void ArcLenLowerPrec2(benchmark::State &state)
{
  double result;
  ITERATIONS = state.range(0);

  for (auto _ : state)
  {
    result = lower_prec::do_fun2(0, 0);

    benchmark::DoNotOptimize(result);
  }
}

static void ArcLenLowerPrec(benchmark::State &state)
{
  double result;
  ITERATIONS = state.range(0);

  for (auto _ : state)
  {
    result = lower_prec::do_fun(0, 0);

    benchmark::DoNotOptimize(result);
  }
}

static void ArcLenHighPrec(benchmark::State &state)
{
  double result;
  ITERATIONS = state.range(0);

  for (auto _ : state)
  {
    result = do_fun(0, 0);

    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(ArcLenLowerPrec2)->RangeMultiplier(10)->Range(10000000, 100000000);
BENCHMARK(ArcLenLowerPrec)->RangeMultiplier(10)->Range(10000000, 100000000);
BENCHMARK(ArcLenHighPrec)->RangeMultiplier(10)->Range(10000000, 100000000);

// Define our main
// BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    ::benchmark::RegisterMemoryManager(mm.get());
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}