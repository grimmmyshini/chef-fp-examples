#include "benchmark/benchmark.h"
#include "../benchmark-utils/memory-manager.hpp"
#include <math.h>
#ifdef VERBOSE 
#include <iostream>
#endif

#include "simpsons.hpp"

static void Simpsons_MixedPrecision(benchmark::State &state)
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

static void Simpsons_HighPrecision(benchmark::State &state)
{
    int n = state.range(0);
    double a = 0, b = 1;

    for (auto _ : state)
    {

        double result_dbl = simpsons<double>(a, b, n);

        benchmark::DoNotOptimize(result_dbl);
    }
}

BENCHMARK(Simpsons_MixedPrecision)->RangeMultiplier(10)->Range(10000, 100000000);
BENCHMARK(Simpsons_HighPrecision)->RangeMultiplier(10)->Range(10000, 100000000);

// Define our main
// BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    ::benchmark::RegisterMemoryManager(mm.get());
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}
