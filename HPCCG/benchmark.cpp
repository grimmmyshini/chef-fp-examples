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

#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include "generate_matrix.hpp"
#include "read_HPC_row.hpp"
#include "HPC_sparsemv.hpp"
#include "compute_residual.hpp"
#include "HPCCG.hpp"
#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"
#include "ddot.hpp"

#undef DEBUG

#include "adapt.h"

#define nx 20
#define ny 30
#define nz 160

struct cout_redirect
{
  cout_redirect(std::streambuf *new_buffer)
      : old(std::cout.rdbuf(new_buffer))
  {
  }

  ~cout_redirect()
  {
    std::cout.rdbuf(old);
  }

private:
  std::streambuf *old;
};

static void ErrorEstimateHPCCGAdapt(benchmark::State &state)
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

  std::stringstream buffer;
  cout_redirect output(buffer.rdbuf());

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
  }

  delete[] p;
  delete[] Ap;
  delete[] r;
}

static void ErrorEstimateHPCCGClad(benchmark::State &state)
{
  double *x, *b, *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;
  double t6 = 0.0;

  clad::generate_matrix(nx, ny, nz, clad::A, &x, &b, &xexact);

  bool dump_matrix = false;
  if (dump_matrix)
    clad::dump_matlab_matrix(clad::A);

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

  double _final_error = 0;

  // cout << "b: ";
  // printVals(b, nrow);
  // cout << "x: ";
  // printVals(x, nrow);
  // cout << "x exact: ";
  // printVals(xexact, nrow);

  std::stringstream buffer;
  cout_redirect output(buffer.rdbuf());

  for (auto _ : state)
  {
    HPCCG_residual_grad(b, x, xexact, r, p, Ap, d_b, d_x, d_xexact, d_r, d_p, d_Ap, _final_error);

    cout << "\nFinal error in HPCCG =" << _final_error << endl;

    clad::printErrorReport();
  }

  // cout << "Gradients are: " << endl;
  // cout << "b: ";
  // printVals(d_b.ptr(), nrow);
  // cout << "x: ";
  // printVals(d_x.ptr(), nrow);

  delete[] b_diff;
  delete[] x_diff;
  delete[] xexact_diff;
  delete[] p;
  delete[] p_diff;
  delete[] Ap;
  delete[] Ap_diff;
  delete[] r;
  delete[] r_diff;

  delete[] x;
  delete[] xexact;
  delete[] b;
}

BENCHMARK(ErrorEstimateHPCCGAdapt)->Unit(benchmark::kSecond)->Iterations(1);

BENCHMARK(ErrorEstimateHPCCGClad)->Unit(benchmark::kSecond)->Iterations(1);

BENCHMARK_MAIN();