#include "benchmark/benchmark.h"

#include <iostream>

#include "arclen.hpp"

static void ArcLength_MixedPrecision(benchmark::State &state)
{
  double result;
  int iters = state.range(0);
  for (auto _ : state) {
    result = do_fun<double, float>(iters);
    benchmark::DoNotOptimize(result);
  }
}

static void ArcLength_HighPrecision(benchmark::State &state)
{
  double result;
  int iters = state.range(0);

  for (auto _ : state)
  {
    result = do_fun<long double>(iters);
    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(ArcLength_HighPrecision)->RangeMultiplier(10)->Range(10000, 100000000);
BENCHMARK(ArcLength_HighPrecision)->RangeMultiplier(10)->Range(10000, 100000000);

// Define our main
BENCHMARK_MAIN();
