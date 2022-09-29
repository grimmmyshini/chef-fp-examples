#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include <string>
#include <cmath>
#include <iomanip>

#include "benchmark/benchmark.h"
#include "../benchmark-utils/memory-manager.hpp"

#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include "generate_matrix.hpp"
#include "read_HPC_row.hpp"
#include "HPC_sparsemv.hpp"
#include "compute_residual.hpp"
#include "HPCCG-adapt.hpp"
#include "HPCCG-clad.hpp"
#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"
#include "ddot.hpp"

#undef DEBUG

#include "adapt.h"


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

static void HPCCG(benchmark::State &state)
{
  double *x, *b, *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;
  double times[7];
  double t6 = 0.0;

  int nx = state.range(0);
  int ny = state.range(1);
  int nz = state.range(2);
  
  generate_matrix(nx, ny, nz, clad::A, &x, &b, &xexact);

  bool dump_matrix = false;
  if (dump_matrix)
    dump_matlab_matrix(clad::A);

  int nrow = clad::A.local_nrow;
  int ncol = clad::A.local_ncol;

  double *r = new double[nrow]();
  double *p = new double[ncol]();
  double *Ap = new double[nrow]();
  int max_iter = 100, niters = 0;
  double tolerance = 0.0, normr = 0.0;

  double residual;

  for (auto _ : state)
  {
    residual = clad::HPCCG(b, x, max_iter, tolerance, niters, normr, r, p, Ap, xexact);
    
    for (int i = 0; i < nrow; i++) {
      x[i] = 0;
      r[i] = 0;
      Ap[i] = 0;
    }
    for (int i = 0; i < ncol; i++) {
      p[i] = 0;
    }
    niters = 0;
    normr = 0;
  }

  delete[] p;
  delete[] Ap;
  delete[] r;

  delete[] x;
  delete[] xexact;
  delete[] b;
}

static void HPCCG_Adapt(benchmark::State &state)
{
  HPC_Sparse_Matrix *A;
  AD_real *x, *b;
  double *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;
  double t6 = 0.0;

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  int nx = state.range(0);
  int ny = state.range(1);
  int nz = state.range(2);

  adapt::generate_matrix(nx, ny, nz, &A, &x, &b, &xexact);

  bool dump_matrix = false;
  if (dump_matrix && size <= 4)
    adapt::dump_matlab_matrix(A, rank);

  int niters = 0;
  AD_real normr = 0.0;
  int max_iter = 100;
  double tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations

  int nrow = A->local_nrow, ncol = A->local_ncol;

  AD_real *r = new AD_real[nrow];
  AD_real *p = new AD_real[ncol]; // In parallel case, A is rectangular
  AD_real *Ap = new AD_real[nrow];

  cout_suppressor suppressor;

  for (auto _ : state)
  {
    AD_begin();
    //  AD_enable_absolute_value_error();
    AD_enable_source_aggregation();
    for (int i = 0; i < A->total_nrow; i++)
    {
      AD_INDEPENDENT(x[i], "x");
      AD_INDEPENDENT(b[i], "b");
    }

    adapt::HPCCG(A, b, x, max_iter, tolerance, niters, normr, r, p, Ap);

    // Compute difference between known exact solution and computed solution
    // All processors are needed here.

    AD_real residual = 0;
    adapt::compute_residual(A->local_nrow, x, xexact, &residual);

    if (rank == 0)
      cout << std::setprecision(5) << "Difference between computed and exact (residual)  = "
           << AD_value(residual) << ".\n"
           << endl;

    AD_DEPENDENT(residual, "residual", 0.0);
    AD_report();
    AD_end();

    for (int i = 0; i < nrow; i++) {
      x[i] = 0;
      r[i] = 0;
      Ap[i] = 0;
    }
    for (int i = 0; i < ncol; i++) {
      p[i] = 0;
    }
    niters = 0;
    normr = 0;
  }

  delete[] p;
  delete[] Ap;
  delete[] r;

  delete[] x;
  delete[] xexact;
  delete[] b;
  delete A;
}

static void HPCCG_Clad(benchmark::State &state)
{
  double *x, *b, *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;
  double t6 = 0.0;

  int nx = state.range(0);
  int ny = state.range(1);
  int nz = state.range(2);

  generate_matrix(nx, ny, nz, clad::A, &x, &b, &xexact);

  bool dump_matrix = false;
  if (dump_matrix)
    dump_matlab_matrix(clad::A);

  int nrow = clad::A.local_nrow;
  int ncol = clad::A.local_ncol;

  // cout << "Float Result: " << executefunction<float>(nrow, ncol, x, b, xexact) << std::endl;
  // cout << "Double Result: " << executefunction<double>(nrow, ncol, x, b, xexact);

  // executeGradient(nrow, ncol, x, b, xexact);

  double *x_diff = new double[nrow]();
  clad::array_ref<double> d_x(x_diff, nrow);

  double *b_diff = new double[nrow]();
  clad::array_ref<double> d_b(b_diff, nrow);

  double *xexact_diff = new double[nrow]();
  clad::array_ref<double> d_xexact(xexact_diff, nrow);

  double *r = new double[nrow]();
  double *r_diff = new double[nrow]();
  clad::array_ref<double> d_r(r_diff, nrow);

  double *p = new double[ncol]();
  double *p_diff = new double[ncol]();
  clad::array_ref<double> d_p(p_diff, ncol);

  double *Ap = new double[nrow]();
  double *Ap_diff = new double[nrow]();
  clad::array_ref<double> d_Ap(Ap_diff, nrow);

  double _final_error = 0, residual;
  int max_iter = 100, d_max_iter, niters = 0, d_niters;
  double tolerance =  0.0, d_tolerance, normr = 0, d_normr;

  // cout << "b: ";
  // printVals(b, nrow);
  // cout << "x: ";
  // printVals(x, nrow);
  // cout << "x exact: ";
  // printVals(xexact, nrow);

  cout_suppressor suppressor;

  for (auto _ : state)
  {
    clad::resetErrors();

    clad::HPCCG_grad(b, x, max_iter, tolerance, niters, normr, r, p, Ap, xexact, d_b, d_x, &d_max_iter, &d_tolerance, &d_niters, &d_normr, d_r, d_p, d_Ap, d_xexact, _final_error);

    cout << "\nFinal error in HPCCG =" << _final_error << endl;

    clad::printErrorReport();

    for (int i = 0; i < nrow; i++) {
      x[i] = 0;
      r[i] = 0;
      Ap[i] = 0;
    }
    for (int i = 0; i < ncol; i++) {
      p[i] = 0;
    }
    niters = 0;
    normr = 0;
  }

  // cout << "Gradients are: " << endl;
  // cout << "b: ";
  // printVals(d_b.ptr(), nrow);
  // cout << "x: ";
  // printVals(d_x.ptr(), nrow);

  delete[] b_diff;
  delete[] x_diff;
  delete[] xexact_diff;
  delete[] p_diff;
  delete[] Ap_diff;
  delete[] r_diff;
  delete[] p;
  delete[] Ap;
  delete[] r;

  delete[] x;
  delete[] xexact;
  delete[] b;
}

BENCHMARK(HPCCG)->Args({20, 30, 10})->Args({20, 30, 20})->Args({20, 30, 40})->Args({20, 30, 80})->Args({20, 30, 160})->Args({20, 30, 320});
BENCHMARK(HPCCG_Clad)->Args({20, 30, 10})->Args({20, 30, 20})->Args({20, 30, 40})->Args({20, 30, 80})->Args({20, 30, 160})->Args({20, 30, 320});
BENCHMARK(HPCCG_Adapt)->Args({20, 30, 10})->Args({20, 30, 20})->Args({20, 30, 40})->Args({20, 30, 80})->Args({20, 30, 160});


// BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    ::benchmark::RegisterMemoryManager(mm.get());
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}