#include "benchmark/benchmark.h"
#include "../benchmark-utils/memory-manager.hpp"

#include "arclen.hpp"
#include "arclen-adapt.hpp"

#include "adapt.h"

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

static void ErrorEstimateArcLen(benchmark::State &state)
{
    int iters = state.range(0);
    for (auto _ : state)
    {
        double result = do_fun<double>(iters);

        benchmark::DoNotOptimize(result);
    }
}

static void ErrorEstimateArcLenAdapt(benchmark::State& state) {
    cout_suppressor suppressor;

    using namespace adapt;
    ITERATIONS = state.range(0);

    for (auto _ : state) {
        AD_begin();
        AD_INTERMEDIATE(h, "h");
        AD_INTERMEDIATE(t1, "t1");
        AD_INTERMEDIATE(s1, "s1");
        AD_INTERMEDIATE(d1, "d1");

        do_fun();

        AD_DEPENDENT(s1, "s1", EPS);
        AD_report();
        t1 = 0.0;
        s1 = 0.0;
        d1 = 1.0;
    }
}

static void ErrorEstimateArcLenClad(benchmark::State& state) {
//   auto df = clad::estimate_error(clad::do_fun);
  
  cout_suppressor suppressor;
  double result;
  double fin_error;
  int iters = state.range(0), diters  = 0;

  for (auto _ : state) {
    ds1 = 0, dt1 = 0;
    fin_error = 0;

    clad::resetErrors();

    // result = do_fun(0, 0);
    clad::do_fun_grad(iters, &diters, fin_error);

    benchmark::DoNotOptimize(fin_error);
    clad::printErrorReport();
  }
}

BENCHMARK(ErrorEstimateArcLen)->Unit(benchmark::kMillisecond)->RangeMultiplier(10)->Range(10000, 100000000);
BENCHMARK(ErrorEstimateArcLenClad)->Unit(benchmark::kMillisecond)->RangeMultiplier(10)->Range(10000, 100000000);
BENCHMARK(ErrorEstimateArcLenAdapt)->Unit(benchmark::kMillisecond)->RangeMultiplier(10)->Range(10000, 10000000);

// Define our main
// BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    ::benchmark::RegisterMemoryManager(mm.get());
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}
