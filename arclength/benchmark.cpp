#include "benchmark/benchmark.h"

#include "arclen.hpp"


#include "adapt.h"
#include "adapt-impl.cpp"

#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include "Derivative.hpp"

struct cout_redirect {
    cout_redirect( std::streambuf * new_buffer ) 
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

private:
    std::streambuf * old;
};

static void ErrorEstimateArcLenAdapt(benchmark::State& state) {
    std::stringstream buffer;
    cout_redirect output(buffer.rdbuf());
    using namespace adapt;

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
  
  std::stringstream buffer;
  cout_redirect output(buffer.rdbuf());

  for (auto _ : state) {
    double ds1 = 0, dt1 = 0;
    double fin_error = 0;

    clad::resetErrors();

    clad::do_fun_grad(0, 0, &ds1, &dt1, fin_error);

    benchmark::DoNotOptimize(fin_error);
    clad::printErrorReport();
  }
}

BENCHMARK(ErrorEstimateArcLenClad)->Unit(benchmark::kSecond)->Iterations(1);
BENCHMARK(ErrorEstimateArcLenAdapt)->Unit(benchmark::kSecond)->Iterations(1);

// Define our main
BENCHMARK_MAIN();