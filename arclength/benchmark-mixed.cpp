#include "benchmark/benchmark.h"

#include <iostream>

#include "arclen.hpp"

static void ArcLenLowerPrec(benchmark::State &state)
{
  double result;
  int iters = state.range(0);
  for (auto _ : state) {
    result = do_fun<double, float>(iters);
    benchmark::DoNotOptimize(result);
  }
}

static void ArcLenHighPrec(benchmark::State &state)
{
  double result;
  int iters = state.range(0);

  for (auto _ : state)
  {
    result = do_fun<double>(iters);
    benchmark::DoNotOptimize(result);
  }
}

static void ArcLenLongPrec(benchmark::State &state)
{
  double result;
  int iters = state.range(0);

  for (auto _ : state)
  {
    result = do_fun<long double>(iters);
    benchmark::DoNotOptimize(result);
  }
}

BENCHMARK(ArcLenLowerPrec)->Unit(benchmark::kMillisecond)->RangeMultiplier(10)->Range(10000, 100000000);
BENCHMARK(ArcLenHighPrec)->Unit(benchmark::kMillisecond)->RangeMultiplier(10)->Range(10000, 100000000);
BENCHMARK(ArcLenLongPrec)->Unit(benchmark::kMillisecond)->RangeMultiplier(10)->Range(10000, 100000000);

// Define our main
BENCHMARK_MAIN();
