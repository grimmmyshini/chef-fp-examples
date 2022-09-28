

#include "benchmark/benchmark.h"
#include <math.h>
#ifdef VERBOSE 
#include <iostream>
#endif

#include "simpsons.hpp"

static void SimpLowerPrec(benchmark::State &state)
{
    int n = state.range(0);
    float a = 0, b = 1;
    #ifdef VERBOSE 
    std::cout << "Difference in lower vs higher precision: "
              << std::fabs(simpsons<double, float>(a, b, n) -
                           simpsons<double>(a, b, n))
              << "\n"
              << "Difference in actual vs lower: "
              << std::fabs(2 - simpsons<double, float>(a, b, n)) << "\n"
              << "Difference in actual vs higher: "
              << std::fabs(2 - simpsons<double>(a, b, n)) << "\n";
    #endif
    for (auto _ : state) {
      double result_dbl = simpsons<double, float>(a, b, n);

      benchmark::DoNotOptimize(result_dbl);
    }
}

static void SimpHighPrec(benchmark::State &state)
{
    int n = state.range(0);
    double a = 0, b = 1;

    for (auto _ : state)
    {

        double result_dbl = simpsons<double>(a, b, n);

        benchmark::DoNotOptimize(result_dbl);
    }
}

BENCHMARK(SimpLowerPrec)->Unit(benchmark::kSecond)->RangeMultiplier(10)->Range(10000, 100000000);
BENCHMARK(SimpHighPrec)->Unit(benchmark::kSecond)->RangeMultiplier(10)->Range(10000, 100000000);

BENCHMARK_MAIN();
