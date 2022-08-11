#include <iomanip>  // necessary for setprecision
#include <iostream> // necessary for 'cout'

#include "benchmark/benchmark.h"

#define ITERATIONS 100000000

#include "simpsons.hpp"

namespace lower_prec
{
    double simpsons(float a, float b)
    {
        int n = ITERATIONS;
        double pi = M_PI;
        float h = (b - a) / (2.0 * n);
        double x = a;
        double tmp;
        double x_pi, sin_x_pi;
        x_pi = a * pi;
        sin_x_pi = sin(x_pi);
        double fa = sin_x_pi;

        x_pi = b * pi;
        sin_x_pi = sin(x_pi);
        double fb = sin_x_pi;
        double s1 = fa + fb;

        for (int l = 0; l < n; l++)
        {
            x = x + h;
            x_pi = x * pi;
            sin_x_pi = sin(x_pi);
            s1 = s1 + 4.0 * sin_x_pi;
            x = x + h;
            x_pi = x * pi;
            sin_x_pi = sin(x_pi);
            s1 = s1 + 2.0 * sin_x_pi;
        }

        s1 = s1 - fb;
        tmp = h * M_PI / 3.0;
        s1 = s1 * tmp;
        return s1;
    }
}

static void SimpLowerPrec(benchmark::State &state)
{
    float a = 0, b = 1;

    for (auto _ : state)
    {
        double result_dbl = lower_prec::simpsons(a, b);

        benchmark::DoNotOptimize(result_dbl);
    }
}

static void SimpHighPrec(benchmark::State &state)
{
    double a = 0, b = 1;

    for (auto _ : state)
    {

        double result_dbl = simpsons(a, b);

        benchmark::DoNotOptimize(result_dbl);
    }
}

BENCHMARK(SimpLowerPrec)->Unit(benchmark::kSecond);
BENCHMARK(SimpHighPrec)->Unit(benchmark::kSecond);

BENCHMARK_MAIN();
