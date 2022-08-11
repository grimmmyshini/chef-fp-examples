#include <sstream>

#include "benchmark/benchmark.h"

#include "simpsons.hpp"

#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include "Derivative.hpp"

struct cout_suppressor
{
  cout_suppressor()
      : buffer(), old(std::cout.rdbuf(buffer.rdbuf()))
  {
  }

  ~cout_suppressor()
  {
    std::cout.rdbuf(old);
  }

private:
  std::stringstream buffer;
  std::streambuf *old;
};

static void Simp(benchmark::State &state)
{
    for (auto _ : state)
    {
        double a = 0, b = 1, d_a = 0, d_b = 0, final_error = 0;

        double result_dbl = simpsons(a, b);
        benchmark::DoNotOptimize(result_dbl);
    }
}

static void SimpWithClad(benchmark::State &state)
{
    cout_suppressor suppress;
    
    for (auto _ : state)
    {
        double a = 0, b = 1, d_a = 0, d_b = 0, final_error = 0;

        double result_dbl = simpsons(a, b);
        benchmark::DoNotOptimize(result_dbl);

        clad::resetErrors();
        clad::simpsons_grad(a, b, &d_a, &d_b, final_error);
        clad::printErrorReport();
    }
}

BENCHMARK(Simp)->Unit(benchmark::kSecond);
BENCHMARK(SimpWithClad)->Unit(benchmark::kSecond);

// Define our main
BENCHMARK_MAIN();