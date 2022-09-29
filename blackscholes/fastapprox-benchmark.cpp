#include "benchmark/benchmark.h"

#include <random>
#include <cmath>
#include "fastapprox/fastapprox/src/fastonebigheader.h"

void exp_bench(benchmark::State &state) {
    double lower_bound = 0;
    double upper_bound = 10000;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    int iterations = state.range(0);

    for (auto _ : state)
    {
        for (int i = 0; i < iterations; i++) {
            double x = unif(re);
            double e = exp(x);
            benchmark::DoNotOptimize(e);
        }
    }
}

void fastexp_bench(benchmark::State &state) {
    double lower_bound = 0;
    double upper_bound = 10000;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    int iterations = state.range(0);

    for (auto _ : state)
    {
        for (int i = 0; i < iterations; i++) {
            double x = unif(re);
            double e = fastexp(x);
            benchmark::DoNotOptimize(e);
        }
    }
}

BENCHMARK(exp_bench)->RangeMultiplier(10)->Range(1, 100000000);
BENCHMARK(fastexp_bench)->RangeMultiplier(10)->Range(1, 100000000);

BENCHMARK_MAIN();