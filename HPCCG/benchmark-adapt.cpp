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

#include "generate_matrix.hpp"
#include "read_HPC_row.hpp"
#include "HPC_sparsemv.hpp"
#include "compute_residual.hpp"
#include "HPCCG-adapt.hpp"

#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"
#include "ddot.hpp"

#undef DEBUG

#include "adapt.h"

#define nx 20
#define ny 30
#define nz 160

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


BENCHMARK(ErrorEstimateHPCCGAdapt)->Unit(benchmark::kSecond);


BENCHMARK_MAIN();
