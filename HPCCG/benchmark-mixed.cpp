#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>

#include "benchmark/benchmark.h"
#include "../benchmark-utils/memory-manager.hpp"

#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"
#include "generate_matrix.hpp"

#undef DEBUG

struct cout_suppressor {
  cout_suppressor() : buffer(), old(std::cout.rdbuf(buffer.rdbuf())) {}

  ~cout_suppressor() { std::cout.rdbuf(old); }

private:
  std::stringstream buffer;
  std::streambuf* old;
};

HPC_Sparse_Matrix A;

namespace lower_loop_prec {
double HPCCG(double const* b, double* x, const int max_iter,
             const double tolerance, int& niters, double& normr, double* r,
             double* p, double* Ap, double const* xexact) {
  // ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr);
  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  normr = 0.0;
  double rtrans = 0.0;
  double oldrtrans = 0.0;
  double alp, bta;

  int rank = 0; // Serial case (not using MPI)

  // // p is of length ncols, copy x to p for sparse MV operation
  // waxpby(nrow, 1.0, x, 0.0, x, p);
  alp = 1.0;
  bta = 0.0;
  if (alp == 1.0) {
    for (int i = 0; i < nrow; i++)
      p[i] = x[i] + bta * x[i];
  } else if (bta == 1.0) {
    for (int i = 0; i < nrow; i++)
      p[i] = alp * x[i] + x[i];
  } else {
    for (int i = 0; i < nrow; i++)
      p[i] = alp * x[i] + bta * x[i];
  }
  // waxpy end

  // HPC_sparsemv(A, p, Ap);
  for (int i = 0; i < nrow; i++) {
    Ap[i] = 0;
    for (int j = 0; j < A.nnz_in_row[i]; j++)
      Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
  }
  // HPC_sparsemv end

  // waxpby(nrow, 1.0, b, -1.0, Ap, r);
  alp = 1.0;
  bta = -1.0;
  if (alp == 1.0) {
    for (int i = 0; i < nrow; i++)
      r[i] = b[i] + bta * Ap[i];
  } else if (bta == 1.0) {
    for (int i = 0; i < nrow; i++)
      r[i] = alp * b[i] + Ap[i];
  } else {
    for (int i = 0; i < nrow; i++)
      r[i] = alp * b[i] + bta * Ap[i];
  }
  // waxpby end

  // ddot(nrow, r, r, &rtrans);
  rtrans = 0.0;
  for (int i = 0; i < nrow; i++)
    rtrans += r[i] * r[i];
  // ddot end

  normr = sqrt(rtrans);

  int k = 1;
  int high_prec_iters = 60;
  for (; k < high_prec_iters; k++) {
    if (k == 1) {
      // waxpby(nrow, 1.0, r, 0.0, r, p);
      alp = 1.0;
      bta = 0.0;
      if (alp == 1.0) {
        for (int i = 0; i < nrow; i++)
          p[i] = r[i] + bta * r[i];
      } else if (bta == 1.0) {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + r[i];
      } else {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + bta * r[i];
      }
      // waxpby end
    } else {
      oldrtrans = rtrans;

      // ddot (nrow, r, r, &rtrans); // 2*nrow ops
      rtrans = 0.0;
      for (int i = 0; i < nrow; i++)
        rtrans += r[i] * r[i];
      // ddot end

      double beta = rtrans / oldrtrans;

      // waxpby(nrow, 1.0, r, beta, p, p); // 2*nrow ops
      alp = 1.0;
      bta = beta;
      if (alp == 1.0) {
        for (int i = 0; i < nrow; i++)
          p[i] = r[i] + bta * p[i];
      } else if (bta == 1.0) {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + p[i];
      } else {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + bta * p[i];
      }
      // waxpby end
    }

    normr = sqrt(rtrans);

    // HPC_sparsemv(A, p, Ap); // 2*nnz ops
    for (int i = 0; i < nrow; i++) {
      Ap[i] = 0;
      for (int j = 0; j < A.nnz_in_row[i]; j++)
        Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
    }
    // HPC_sparsemv end

    double alpha = 0.0;

    // ddot(nrow, p, Ap, &alpha); // 2*nrow ops
    alpha = 0.0;
    if (Ap == p)
      for (int i = 0; i < nrow; i++)
        alpha += p[i] * p[i];
    else
      for (int i = 0; i < nrow; i++)
        alpha += p[i] * Ap[i];
    // ddot end

    alpha = rtrans / alpha;

    // waxpby(nrow, 1.0, x, alpha, p, x);   // 2*nrow ops
    alp = 1.0;
    bta = alpha;
    if (alp == 1.0) {
      for (int i = 0; i < nrow; i++)
        x[i] = x[i] + bta * p[i];
    } else if (bta == 1.0) {
      for (int i = 0; i < nrow; i++)
        x[i] = alp * x[i] + p[i];
    } else {
      for (int i = 0; i < nrow; i++)
        x[i] = alp * x[i] + bta * p[i];
    }
    // waxpby end

    // waxpby(nrow, 1.0, r, -alpha, Ap, r); // 2*nrow ops
    alp = 1.0;
    bta = -alpha;
    if (alp == 1.0) {
      for (int i = 0; i < nrow; i++)
        r[i] = r[i] + bta * Ap[i];
    } else if (bta == 1.0) {
      for (int i = 0; i < nrow; i++)
        r[i] = alp * r[i] + Ap[i];
    } else {
      for (int i = 0; i < nrow; i++)
        r[i] = alp * r[i] + bta * Ap[i];
    }
    // waxpby end

    niters = k;
  }

  float* rf = new float[nrow];
  float* pf = new float[ncol]; // In parallel case, A is rectangular
  float* Apf = new float[nrow];
  float* xf = new float[nrow];
  float oldrtransf = oldrtrans;
  float rtransf = rtrans;
  float alpf = alp, btaf = bta;

  for (int i = 0; i < ncol; i++) {
    pf[i] = p[i];
  }
  for (int i = 0; i < nrow; i++) {
    rf[i] = r[i];
    Apf[i] = Ap[i];
    xf[i] = x[i];
  }

  for (; k < max_iter; k++) {
    oldrtransf = rtransf;

    // ddot (nrow, r, r, &rtrans); // 2*nrow ops
    rtransf = 0.0;
    for (int i = 0; i < nrow; i++)
      rtransf += rf[i] * rf[i];
    // ddot end

    float beta = rtransf / oldrtransf;

    // waxpby(nrow, 1.0, r, beta, p, p); // 2*nrow ops
    alpf = 1.0;
    btaf = beta;
    if (alpf == 1.0) {
      for (int i = 0; i < nrow; i++)
        pf[i] = rf[i] + btaf * pf[i];
    } else if (btaf == 1.0) {
      for (int i = 0; i < nrow; i++)
        pf[i] = alpf * rf[i] + pf[i];
    } else {
      for (int i = 0; i < nrow; i++)
        pf[i] = alpf * rf[i] + btaf * pf[i];
    }
    // waxpby end

    normr = sqrt(rtransf);

    // HPC_sparsemv(A, p, Ap); // 2*nnz ops
    for (int i = 0; i < nrow; i++) {
      Apf[i] = 0;
      for (int j = 0; j < A.nnz_in_row[i]; j++)
        Apf[i] += A.ptr_to_vals_in_row_f[i][j] * pf[A.ptr_to_inds_in_row[i][j]];
    }
    // HPC_sparsemv end

    float alpha = 0.0;

    // ddot(nrow, p, Ap, &alpha); // 2*nrow ops
    alpha = 0.0;
    if (Apf == pf)
      for (int i = 0; i < nrow; i++)
        alpha += pf[i] * pf[i];
    else
      for (int i = 0; i < nrow; i++)
        alpha += pf[i] * Apf[i];
    // ddot end

    alpha = rtransf / alpha;

    // waxpby(nrow, 1.0, x, alpha, p, x);   // 2*nrow ops
    alpf = 1.0;
    btaf = alpha;
    if (alpf == 1.0) {
      for (int i = 0; i < nrow; i++)
        xf[i] = xf[i] + btaf * pf[i];
    } else if (btaf == 1.0) {
      for (int i = 0; i < nrow; i++)
        xf[i] = alpf * xf[i] + pf[i];
    } else {
      for (int i = 0; i < nrow; i++)
        xf[i] = alpf * xf[i] + btaf * pf[i];
    }
    // waxpby end

    // waxpby(nrow, 1.0, r, -alpha, Ap, r); // 2*nrow ops
    alpf = 1.0;
    btaf = -alpha;
    if (alpf == 1.0) {
      for (int i = 0; i < nrow; i++)
        rf[i] = rf[i] + btaf * Apf[i];
    } else if (btaf == 1.0) {
      for (int i = 0; i < nrow; i++)
        rf[i] = alpf * rf[i] + Apf[i];
    } else {
      for (int i = 0; i < nrow; i++)
        rf[i] = alpf * rf[i] + btaf * Apf[i];
    }
    // waxpby end

    niters = k;
  }

  // Compute difference between known exact solution and computed solution
  // All processors are needed here.

  double residual = 0;

  oldrtrans = oldrtransf;
  rtrans = rtransf;

  for (int i = 0; i < ncol; i++) {
    p[i] = pf[i];
  }
  for (int i = 0; i < nrow; i++) {
    r[i] = rf[i];
    Ap[i] = Apf[i];
    x[i] = xf[i];
  }

  delete[] rf;
  delete[] pf;
  delete[] Apf;
  delete[] xf;

  // compute residual(A.local_nrow, x, xexact, &residual);
  for (int i = 0; i < nrow; i++) {
    double diff = fabs(x[i] - xexact[i]);
    if (diff > residual)
      residual = diff;
  }

  return residual;
}
} // namespace lower_loop_prec

namespace high_prec {
double HPCCG(double const* b, double* x, const int max_iter,
             const double tolerance, int& niters, double& normr, double* r,
             double* p, double* Ap, double const* xexact) {
  // ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr);
  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  normr = 0.0;
  double rtrans = 0.0;
  double oldrtrans = 0.0;
  double alp, bta;

  int rank = 0; // Serial case (not using MPI)

  // // p is of length ncols, copy x to p for sparse MV operation
  // waxpby(nrow, 1.0, x, 0.0, x, p);
  alp = 1.0;
  bta = 0.0;
  if (alp == 1.0) {
    for (int i = 0; i < nrow; i++)
      p[i] = x[i] + bta * x[i];
  } else if (bta == 1.0) {
    for (int i = 0; i < nrow; i++)
      p[i] = alp * x[i] + x[i];
  } else {
    for (int i = 0; i < nrow; i++)
      p[i] = alp * x[i] + bta * x[i];
  }
  // waxpy end

  // HPC_sparsemv(A, p, Ap);
  for (int i = 0; i < nrow; i++) {
    Ap[i] = 0;
    for (int j = 0; j < A.nnz_in_row[i]; j++)
      Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
  }
  // HPC_sparsemv end

  // waxpby(nrow, 1.0, b, -1.0, Ap, r);
  alp = 1.0;
  bta = -1.0;
  if (alp == 1.0) {
    for (int i = 0; i < nrow; i++)
      r[i] = b[i] + bta * Ap[i];
  } else if (bta == 1.0) {
    for (int i = 0; i < nrow; i++)
      r[i] = alp * b[i] + Ap[i];
  } else {
    for (int i = 0; i < nrow; i++)
      r[i] = alp * b[i] + bta * Ap[i];
  }
  // waxpby end

  // ddot(nrow, r, r, &rtrans);
  rtrans = 0.0;
  for (int i = 0; i < nrow; i++)
    rtrans += r[i] * r[i];
  // ddot end

  normr = sqrt(rtrans);

  for (int k = 1; k < max_iter; k++) {
    if (k == 1) {
      // waxpby(nrow, 1.0, r, 0.0, r, p);
      alp = 1.0;
      bta = 0.0;
      if (alp == 1.0) {
        for (int i = 0; i < nrow; i++)
          p[i] = r[i] + bta * r[i];
      } else if (bta == 1.0) {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + r[i];
      } else {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + bta * r[i];
      }
      // waxpby end
    } else {
      oldrtrans = rtrans;

      // ddot (nrow, r, r, &rtrans); // 2*nrow ops
      rtrans = 0.0;
      for (int i = 0; i < nrow; i++)
        rtrans += r[i] * r[i];
      // ddot end

      double beta = rtrans / oldrtrans;

      // waxpby(nrow, 1.0, r, beta, p, p); // 2*nrow ops
      alp = 1.0;
      bta = beta;
      if (alp == 1.0) {
        for (int i = 0; i < nrow; i++)
          p[i] = r[i] + bta * p[i];
      } else if (bta == 1.0) {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + p[i];
      } else {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + bta * p[i];
      }
      // waxpby end
    }

    normr = sqrt(rtrans);

    // HPC_sparsemv(A, p, Ap); // 2*nnz ops
    for (int i = 0; i < nrow; i++) {
      Ap[i] = 0;
      for (int j = 0; j < A.nnz_in_row[i]; j++)
        Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
    }
    // HPC_sparsemv end

    double alpha = 0.0;

    // ddot(nrow, p, Ap, &alpha); // 2*nrow ops
    alpha = 0.0;
    if (Ap == p)
      for (int i = 0; i < nrow; i++)
        alpha += p[i] * p[i];
    else
      for (int i = 0; i < nrow; i++)
        alpha += p[i] * Ap[i];
    // ddot end

    alpha = rtrans / alpha;

    // waxpby(nrow, 1.0, x, alpha, p, x);   // 2*nrow ops
    alp = 1.0;
    bta = alpha;
    if (alp == 1.0) {
      for (int i = 0; i < nrow; i++)
        x[i] = x[i] + bta * p[i];
    } else if (bta == 1.0) {
      for (int i = 0; i < nrow; i++)
        x[i] = alp * x[i] + p[i];
    } else {
      for (int i = 0; i < nrow; i++)
        x[i] = alp * x[i] + bta * p[i];
    }
    // waxpby end

    // waxpby(nrow, 1.0, r, -alpha, Ap, r); // 2*nrow ops
    alp = 1.0;
    bta = -alpha;
    if (alp == 1.0) {
      for (int i = 0; i < nrow; i++)
        r[i] = r[i] + bta * Ap[i];
    } else if (bta == 1.0) {
      for (int i = 0; i < nrow; i++)
        r[i] = alp * r[i] + Ap[i];
    } else {
      for (int i = 0; i < nrow; i++)
        r[i] = alp * r[i] + bta * Ap[i];
    }
    // waxpby end

    niters = k;
  }

  // Compute difference between known exact solution and computed solution
  // All processors are needed here.

  double residual = 0;

  // compute residual(A.local_nrow, x, xexact, &residual);
  for (int i = 0; i < nrow; i++) {
    double diff = fabs(x[i] - xexact[i]);
    if (diff > residual)
      residual = diff;
  }

  return residual;
}

} // namespace high_prec

void HPCCG_MixedPrecision(benchmark::State& state) {
  double *x, *b, *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;

  int nx = state.range(0);
  int ny = state.range(1);
  int nz = state.range(2);

  int size = 1; // Serial case (not using MPI)

  generate_matrix(nx, ny, nz, A, &x, &b, &xexact);

  bool dump_matrix = false;
  if (dump_matrix && size <= 4)
    dump_matlab_matrix(A);

  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  double* r = new double[nrow];
  double* p = new double[ncol]; // In parallel case, A is rectangular
  double* Ap = new double[nrow];
  double residual;

  int niters = 0;
  double normr = 0.0;
  int max_iter = 100;
  double tolerance =
      0.0; // Set tolerance to zero to make all runs do max_iter iterations

  cout_suppressor suppressor;
  for (auto _ : state) {
    residual = lower_loop_prec::HPCCG(b, x, max_iter, tolerance, niters,
                                      normr, r, p, Ap, xexact);

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
  delete[] xexact;
  delete[] b;
  delete[] x;
  delete[] r;

  cout << std::setprecision(5)
       << "Difference between computed and exact (residual)  = " << residual
       << ".\n"
       << endl;
}

static void HPCCG_HighPrecision(benchmark::State& state) {
  double *x, *b, *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;

  int nx = state.range(0);
  int ny = state.range(1);
  int nz = state.range(2);

  int size = 1; // Serial case (not using MPI)

  generate_matrix(nx, ny, nz, A, &x, &b, &xexact);

  bool dump_matrix = false;
  if (dump_matrix && size <= 4)
    dump_matlab_matrix(A);

  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  double* r = new double[nrow];
  double* p = new double[ncol]; // In parallel case, A is rectangular
  double* Ap = new double[nrow];
  double residual;

  int niters = 0;
  double normr = 0.0;
  int max_iter = 100;
  double tolerance =
      0.0; // Set tolerance to zero to make all runs do max_iter iterations

  cout_suppressor suppressor;
  for (auto _ : state) {
    residual = high_prec::HPCCG(b, x, max_iter, tolerance, niters, normr, r,
                                p, Ap, xexact);

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
  delete[] xexact;
  delete[] b;
  delete[] x;
  delete[] r;

  cout << std::setprecision(5)
       << "Difference between computed and exact (residual)  = " << residual
       << ".\n"
       << endl;
}

BENCHMARK(HPCCG_MixedPrecision)->Args({20, 30, 10})->Args({20, 30, 20})->Args({20, 30, 40})->Args({20, 30, 80})->Args({20, 30, 160})->Args({20, 30, 320});
BENCHMARK(HPCCG_HighPrecision)->Args({20, 30, 10})->Args({20, 30, 20})->Args({20, 30, 40})->Args({20, 30, 80})->Args({20, 30, 160})->Args({20, 30, 320});

// BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    ::benchmark::RegisterMemoryManager(mm.get());
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}