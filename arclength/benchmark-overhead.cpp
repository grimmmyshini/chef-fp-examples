#include <sstream>

#include "benchmark/benchmark.h"

#include "arclen.hpp"

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

static void ArcLen(benchmark::State &state)
{
    for (auto _ : state)
    {
        double result = do_fun(0, 0);

        benchmark::DoNotOptimize(result);
    }
}

static void ArcLenWithClad(benchmark::State &state)
{
    cout_suppressor suppressor;

    for (auto _ : state)
    {
        double ds1 = 0, dt1 = 0;
        double fin_error = 0;

        double result = do_fun(0, 0);

        clad::resetErrors();
        clad::do_fun_grad(0, 0, &ds1, &dt1, fin_error);
        clad::printErrorReport();

        benchmark::DoNotOptimize(result);
    }
}

BENCHMARK(ArcLen)->Unit(benchmark::kSecond);
BENCHMARK(ArcLenWithClad)->Unit(benchmark::kSecond);

// Define our main
BENCHMARK_MAIN();