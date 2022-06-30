/* 
 * This benchmark outputs only the final error contribution calculated by clad.
 * (See lines 201, 211, 219, 223, 298, 314, 324, 327 for ommitted outputs) 
 * Clad outputs to a file while adapt outputs to a stringstream buffer.
 * 
 * The outputs are redirected to a stringstream buffer so as to separate
 * them from the benchmark results.
 * 
 * Example adapted from ADAPT-FP: https://github.com/LLNL/adapt-fp
 *
 * Typical invocation of the benchmark looks like this:
 * clang -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang /path/to/clad.so 
 * -Ipath/to/clad/include -Ipath/to/adapt-fp -x c++ -lstdc++ -O2
 * -lm -isystem path/to/Google/benchmark/include  -Lpath/to/benchmark/build/src
 * ./arclen_full_buf.cpp -lbenchmark -lpthread -DCODI_ZeroAdjointReverse=0
 * 
 */

#include "benchmark/benchmark.h"

#include <cmath>
#include <fstream> // Necessary for saving error estimates.
#include <iostream>
#include <stdio.h>
#include <math.h>

#include "adapt.h"
#include "adapt-impl.cpp"

#include "clad/Differentiator/Differentiator.h"
/*
 * Precimonious output labels located in comments beside variable declarations
 */

#define ABS(x) ( ((x) < 0.0) ? (-(x)) : (x) )

#define PI  3.1415926535897932L
#define ANS 5.795776322412856L
#define EPS 1e-10

#define N 1000000

AD_real h  = PI / (double)N;     // double

AD_real t1 = 0.0;
AD_real t2;
AD_real t3;                      // double

AD_real s1 = 0.0;                // double

AD_real d1 = 1.0;
AD_real d2;                      // float

AD_real fun (AD_real x)
{
    d2 = d1;    // also d1 in original
    AD_INTERMEDIATE(d2, "d2");
    t3 = x;     // also t1 in original
    AD_INTERMEDIATE(t3, "t3");

    int k;
    for (k = 1; k <= 5; k+=1)
    {
        d2 = 2.0 * d2;
        AD_INTERMEDIATE(d2, "d2");
        t3 = t3 + sin (d2 * x) / d2;
        AD_INTERMEDIATE(t3, "t3");
    }
    return t3;
}

void do_fun ()
{
    int i;
    for (i = 1; i <= N; i+=1)
    {
        t2 = fun (i * h);
        AD_INTERMEDIATE(t2, "t2");
        s1 = s1 + sqrt (h * h + (t2 - t1) * (t2 - t1));
        AD_INTERMEDIATE(s1, "s1");
        t1 = t2;
        AD_INTERMEDIATE(t1, "t1");
    }
}

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
    for (auto _ : state) {
        AD_begin();
        AD_INTERMEDIATE(h, "h");
        AD_INTERMEDIATE(t1, "t1");
        AD_INTERMEDIATE(s1, "s1");
        AD_INTERMEDIATE(d1, "d1");

        do_fun();

        AD_DEPENDENT(s1, "s1", EPS);
        AD_report();
    }
}

BENCHMARK(ErrorEstimateArcLenAdapt)->Unit(benchmark::kSecond)->Iterations(10);

double fun_clad(double x) {
  double d2 = 1.0; // also d1 in original
  double t3 = x;   // also t1 in original

  int k;
  for (k = 1; k <= 5; k += 1) {
    d2 = 2.0 * d2;
    t3 = t3 + sin(d2 * x) / d2;
  }
  return t3;
}

double do_fun_clad(double s1, double t1) {
  double t2, h = PI / (double)N;

  for (int i = 1; i <= N; i += 1) {
    t2 = fun_clad(i * h);
    s1 = s1 + sqrt(h * h + (t2 - t1) * (t2 - t1));
    t1 = t2;
  }

  return s1;
}

// void do_fun_clad_grad(double, double, clad::array_ref<double>, clad::array_ref<double>, double &, std::ostream &);

void fun_clad_pullback(double x, double _d_y, clad::array_ref<double> _d_x, double &_final_error, std::ostream &_error_stream) {
    double _d_d2 = 0;
    double _delta_d20 = 0;
    double _EERepl_d200;
    double _d_t3 = 0;
    double _delta_t30 = 0;
    double _EERepl_t300;
    int _d_k = 0;
    unsigned long _t0;
    clad::tape<double> _t1 = {};
    clad::tape<double> _EERepl_d201 = {};
    clad::tape<double> _t2 = {};
    clad::tape<double> _t3 = {};
    clad::tape<double> _t4 = {};
    clad::tape<double> _t5 = {};
    clad::tape<double> _t6 = {};
    clad::tape<double> _EERepl_t301 = {};
    double d20 = 1.;
    _EERepl_d200 = d20;
    double t30 = x;
    _EERepl_t300 = t30;
    int k;
    int _EERepl_k0 = 1;
    int _EERepl_k1 = 1;
    _t0 = 0;
    for (k = 1; k <= 5; k += 1) {
        _t0++;
        d20 = 2. * clad::push(_t1, d20);
        clad::push(_EERepl_d201, d20);
        t30 = t30 + clad::push(_t6, sin(clad::push(_t5, clad::push(_t4, d20) * clad::push(_t3, x)))) / clad::push(_t2, d20);
        clad::push(_EERepl_t301, t30);
    }
    double fun_clad_return = t30;
    goto _label0;
  _label0:
    _d_t3 += _d_y;
    for (; _t0; _t0--) {
        {
            int _r_d1 = _d_k;
            _d_k += _r_d1;
            _d_k -= _r_d1;
        }
        {
            {
                double _r_d3 = _d_t3;
                _d_t3 += _r_d3;
                double _r3 = clad::pop(_t2);
                double _r4 = _r_d3 / _r3;
                double _r5 = _r4 * clad::custom_derivatives::sin_pushforward(clad::pop(_t5), 1.).pushforward;
                double _r6 = _r5 * clad::pop(_t3);
                _d_d2 += _r6;
                double _r7 = clad::pop(_t4) * _r5;
                * _d_x += _r7;
                double _r8 = _r_d3 * -clad::pop(_t6) / (_r3 * _r3);
                _d_d2 += _r8;
                double _r9 = clad::pop(_EERepl_t301);
                _delta_t30 += std::abs(_r_d3 * _r9 * 1.1920928955078125E-7);
                // _error_stream << "t30" << " : " << std::abs(_r_d3 * _r9 * 1.1920928955078125E-7) << "\n";
                _d_t3 -= _r_d3;
            }
            {
                double _r_d2 = _d_d2;
                double _r0 = _r_d2 * clad::pop(_t1);
                double _r1 = 2. * _r_d2;
                _d_d2 += _r1;
                double _r2 = clad::pop(_EERepl_d201);
                _delta_d20 += std::abs(_r_d2 * _r2 * 1.1920928955078125E-7);
                // _error_stream << "d20" << " : " << std::abs(_r_d2 * _r2 * 1.1920928955078125E-7) << "\n";
                _d_d2 -= _r_d2;
            }
        }
    }
    * _d_x += _d_t3;
    {
        _delta_d20 += std::abs(_d_d2 * _EERepl_d200 * 1.1920928955078125E-7);
        // _error_stream << "d20" << " : " << std::abs(_d_d2 * _EERepl_d200 * 1.1920928955078125E-7) << "\n";
    }
    double _delta_x = 0;
    _delta_x += std::abs(* _d_x * x * 1.1920928955078125E-7);
    // _error_stream << "x" << " : " << std::abs(* _d_x * x * 1.1920928955078125E-7) << "\n";
    _final_error += _delta_x + _delta_t30 + _delta_d20;
    _error_stream << "\nFinal error contribution by x = " << _delta_x << "\n";
    _error_stream << "\nFinal error contribution by t30 = " << _delta_t30 << "\n";
    _error_stream << "\nFinal error contribution by d20 = " << _delta_d20 << "\n";
}

void do_fun_clad_grad(double s1, double t1, clad::array_ref<double> _d_s1, clad::array_ref<double> _d_t1, double &_final_error, std::ostream &_error_stream) {
    long double _t0;
    double _d_t2 = 0, _d_h = 0;
    double _delta_t20 = 0;
    double _EERepl_t200;
    double _delta_t200 = 0;
    double _EERepl_t201;
    unsigned long _t1;
    int _d_i = 0;
    int _EERepl_i0 = 1;
    clad::tape<double> _t2 = {};
    clad::tape<int> _t3 = {};
    clad::tape<double> _t4 = {};
    clad::tape<double> _EERepl_t202 = {};
    double _delta_s1 = 0;
    double _EERepl_s10 = s1;
    clad::tape<double> _t6 = {};
    clad::tape<double> _t7 = {};
    clad::tape<double> _t8 = {};
    clad::tape<double> _t9 = {};
    clad::tape<double> _t10 = {};
    clad::tape<double> _EERepl_s11 = {};
    _t0 = (double)1000000;
    double t20, h0 = 3.14159265358979319992L / _t0;
    _EERepl_t201 = t20;
    _EERepl_t200 = t20;
    _t1 = 0;
    for (int i = 1; i <= 1000000; i += 1) {
        _t1++;
        t20 = fun_clad(clad::push(_t4, clad::push(_t3, i) * clad::push(_t2, h0)));
        clad::push(_EERepl_t202, t20);
        s1 = s1 + sqrt(clad::push(_t10, clad::push(_t7, h0) * clad::push(_t6, h0) + clad::push(_t9, (t20 - t1)) * clad::push(_t8, (t20 - t1))));
        clad::push(_EERepl_s11, s1);
        t1 = t20;
    }
    double do_fun_clad_return = s1;
    goto _label0;
  _label0:
    * _d_s1 += 1;
    for (; _t1; _t1--) {
        {
            int _r_d0 = _d_i;
            _d_i += _r_d0;
            _d_i -= _r_d0;
        }
        {
            {
                double _r_d3 = * _d_t1;
                _d_t2 += _r_d3;
                * _d_t1 -= _r_d3;
                * _d_t1;
            }
            {
                double _r_d2 = * _d_s1;
                * _d_s1 += _r_d2;
                double _r6 = _r_d2 * clad::custom_derivatives::sqrt_pushforward(clad::pop(_t10), 1.).pushforward;
                double _r7 = _r6 * clad::pop(_t6);
                _d_h += _r7;
                double _r8 = clad::pop(_t7) * _r6;
                _d_h += _r8;
                double _r9 = _r6 * clad::pop(_t8);
                _d_t2 += _r9;
                * _d_t1 += -_r9;
                double _r10 = clad::pop(_t9) * _r6;
                _d_t2 += _r10;
                * _d_t1 += -_r10;
                double _r11 = clad::pop(_EERepl_s11);
                _delta_s1 += std::abs(_r_d2 * _r11 * 1.1920928955078125E-7);
                // _error_stream << "s1" << " : " << std::abs(_r_d2 * _r11 * 1.1920928955078125E-7) << "\n";
                * _d_s1 -= _r_d2;
                * _d_s1;
            }
            {
                double _r_d1 = _d_t2;
                double _grad0 = 0.;
                double _t5 = 0;
                fun_clad_pullback(clad::pop(_t4), _r_d1, &_grad0, _t5, _error_stream);
                double _r2 = _grad0;
                double _r3 = _r2 * clad::pop(_t2);
                _d_i += _r3;
                double _r4 = clad::pop(_t3) * _r2;
                _d_h += _r4;
                double _r5 = clad::pop(_EERepl_t202);
                _delta_t20 += _t5;
                // _error_stream << "t20" << " : " << _t5 << "\n";
                _d_t2 -= _r_d1;
            }
        }
    }
    {
        long double _r0 = _d_h / _t0;
        long double _r1 = _d_h * -3.14159265358979319992L / (_t0 * _t0);
    }
    _delta_s1 += std::abs(* _d_s1 * _EERepl_s10 * 1.1920928955078125E-7);
    // _error_stream << "s1" << " : " << std::abs(* _d_s1 * _EERepl_s10 * 1.1920928955078125E-7) << "\n";
    double _delta_t1 = 0;
    _delta_t1 += std::abs(* _d_t1 * t1 * 1.1920928955078125E-7);
    // _error_stream << "t1" << " : " << std::abs(* _d_t1 * t1 * 1.1920928955078125E-7) << "\n";
    _final_error += _delta_t1 + _delta_s1 + _delta_t20;
    _error_stream << "\nFinal error contribution by t1 = " << _delta_t1 << "\n";
    _error_stream << "\nFinal error contribution by s1 = " << _delta_s1 << "\n";
    _error_stream << "\nFinal error contribution by t20 = " << _delta_t20 << "\n";
}

static void ErrorEstimateArcLenClad(benchmark::State& state) {
  // double s1 = do_fun_clad(0.0, 0.0);
  // auto df = clad::estimate_error<true>(do_fun_clad);
  // df.dump();
  // Enable before using std::cout for measurements.
  // std::ios_base::sync_with_stdio(false);
  std::ofstream fil("part_err_file.txt");
  // std::stringstream buffer;
  // cout_redirect output(buffer.rdbuf());
  for (auto _ : state) {
    double ds1 = 0, dt1 = 0;
    double fin_error = 0;
    do_fun_clad_grad(0, 0, &ds1, &dt1, fin_error, fil);
  }
}

BENCHMARK(ErrorEstimateArcLenClad)->Unit(benchmark::kSecond)->Iterations(10);

// Define our main
BENCHMARK_MAIN();
