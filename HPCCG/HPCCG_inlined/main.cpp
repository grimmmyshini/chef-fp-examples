
//@HEADER
// ************************************************************************
//
//               HPCCG: Simple Conjugate Gradient Benchmark Code
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

// Changelog
//
// Version 0.3
// - Added timing of setup time for sparse MV
// - Corrected percentages reported for sparse MV with overhead
//
/////////////////////////////////////////////////////////////////////////

// Main routine of a program that reads a sparse matrix, right side
// vector, solution vector and initial guess from a file  in HPC
// format.  This program then calls the HPCCG conjugate gradient
// solver to solve the problem, and then prints results.

// Calling sequence:

// test_HPCCG linear_system_file

// Routines called:

// read_HPC_row - Reads in linear system

// mytimer - Timing routine (compile with -DWALL to get wall clock
//           times

// HPCCG - CG Solver

// compute_residual - Compares HPCCG solution to known solution.

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

#include "generate_matrix.hpp"
#include "read_HPC_row.hpp"
#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"

#undef DEBUG

HPC_Sparse_Matrix A;

double HPCCG_residual(double *b, double *x, int max_iter, double tolerance,
                      int &niters, double &normr, double *r, double *p,
                      double *Ap, double *xexact)
{
  // ierr = HPCCG(A, b, x, max_iter, tolerance, niters, normr);
  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  normr = 0.0;
  double rtrans = 0.0;
  double oldrtrans = 0.0;
  double alp, bta;

  int rank = 0; // Serial case (not using MPI)

  int print_freq = max_iter / 10;
  if (print_freq > 50)
    print_freq = 50;
  if (print_freq < 1)
    print_freq = 1;

  // // p is of length ncols, copy x to p for sparse MV operation
  // waxpby(nrow, 1.0, x, 0.0, x, p);
  alp = 1.0;
  bta = 0.0;
  if (alp == 1.0)
  {
    for (int i = 0; i < nrow; i++)
      p[i] = x[i] + bta * x[i];
  }
  else if (bta == 1.0)
  {
    for (int i = 0; i < nrow; i++)
      p[i] = alp * x[i] + x[i];
  }
  else
  {
    for (int i = 0; i < nrow; i++)
      p[i] = alp * x[i] + bta * x[i];
  }
  // waxpy end

  // HPC_sparsemv(A, p, Ap);
  for (int i = 0; i < nrow; i++)
  {
    double sum = 0.0;

    for (int j = 0; j < A.nnz_in_row[i]; j++)
      sum += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
    Ap[i] = sum;
  }
  // HPC_sparsemv end

  // waxpby(nrow, 1.0, b, -1.0, Ap, r);
  alp = 1.0;
  bta = -1.0;
  if (alp == 1.0)
  {
    for (int i = 0; i < nrow; i++)
      r[i] = b[i] + bta * Ap[i];
  }
  else if (bta == 1.0)
  {
    for (int i = 0; i < nrow; i++)
      r[i] = alp * b[i] + Ap[i];
  }
  else
  {
    for (int i = 0; i < nrow; i++)
      r[i] = alp * b[i] + bta * Ap[i];
  }
  // waxpby end

  // ddot(nrow, r, r, &rtrans);
  rtrans = 0.0;
  if (r == r)
    for (int i = 0; i < nrow; i++)
      rtrans += r[i] * r[i];
  else
    for (int i = 0; i < nrow; i++)
      rtrans += r[i] * r[i];
  // ddot end

  normr = sqrt(rtrans);

  if (rank == 0)
    cout << "Initial Residual = " << normr << endl;

  for (int k = 1; k < max_iter && normr > tolerance; k++)
  {
    if (k == 1)
    {
      // waxpby(nrow, 1.0, r, 0.0, r, p);
      alp = 1.0;
      bta = 0.0;
      if (alp == 1.0)
      {
        for (int i = 0; i < nrow; i++)
          p[i] = r[i] + bta * r[i];
      }
      else if (bta == 1.0)
      {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + r[i];
      }
      else
      {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + bta * r[i];
      }
      // waxpby end
    }
    else
    {
      oldrtrans = rtrans;

      // ddot (nrow, r, r, &rtrans); // 2*nrow ops
      rtrans = 0.0;
      if (r == r)
        for (int i = 0; i < nrow; i++)
          rtrans += r[i] * r[i];
      else
        for (int i = 0; i < nrow; i++)
          rtrans += r[i] * r[i];
      // ddot end

      double beta = rtrans / oldrtrans;

      // waxpby(nrow, 1.0, r, beta, p, p); // 2*nrow ops
      alp = 1.0;
      bta = beta;
      if (alp == 1.0)
      {
        for (int i = 0; i < nrow; i++)
          p[i] = r[i] + bta * p[i];
      }
      else if (bta == 1.0)
      {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + p[i];
      }
      else
      {
        for (int i = 0; i < nrow; i++)
          p[i] = alp * r[i] + bta * p[i];
      }
      // waxpby end
    }

    normr = sqrt(rtrans);
    if (rank == 0 && (k % print_freq == 0 || k + 1 == max_iter))
      cout << "Iteration = " << k << "   Residual = " << normr << endl;

    // HPC_sparsemv(A, p, Ap); // 2*nnz ops
    for (int i = 0; i < nrow; i++)
    {
      double sum = 0.0;

      for (int j = 0; j < A.nnz_in_row[i]; j++)
        sum += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
      Ap[i] = sum;
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
    if (alp == 1.0)
    {
      for (int i = 0; i < nrow; i++)
        x[i] = x[i] + bta * p[i];
    }
    else if (bta == 1.0)
    {
      for (int i = 0; i < nrow; i++)
        x[i] = alp * x[i] + p[i];
    }
    else
    {
      for (int i = 0; i < nrow; i++)
        x[i] = alp * x[i] + bta * p[i];
    }
    // waxpby end

    // waxpby(nrow, 1.0, r, -alpha, Ap, r); // 2*nrow ops
    alp = 1.0;
    bta = -alpha;
    if (alp == 1.0)
    {
      for (int i = 0; i < nrow; i++)
        r[i] = r[i] + bta * Ap[i];
    }
    else if (bta == 1.0)
    {
      for (int i = 0; i < nrow; i++)
        r[i] = alp * r[i] + Ap[i];
    }
    else
    {
      for (int i = 0; i < nrow; i++)
        r[i] = alp * r[i] + bta * Ap[i];
    }
    // waxpby end

    niters = k;
  }

  // Compute difference between known exact solution and computed solution
  // All processors are needed here.

  double residual = 0;

  // compute_residual(A.local_nrow, x, xexact, &residual);
  for (int i = 0; i < nrow; i++)
  {
    double diff = fabs(x[i] - xexact[i]);
    if (diff > residual)
      residual = diff;
  }

  return residual;
}

int main(int argc, char *argv[])
{
  double *x, *b, *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;

  int nx, ny, nz;

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  if (argc != 2 && argc != 4)
  {
    if (rank == 0)
      cerr << "Usage:" << endl
           << "Mode 1: " << argv[0] << " nx ny nz" << endl
           << "     where nx, ny and nz are the local sub-block dimensions, or" << endl
           << "Mode 2: " << argv[0] << " HPC_data_file " << endl
           << "     where HPC_data_file is a globally accessible file containing matrix data." << endl;
    exit(1);
  }

  if (argc == 4)
  {
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nz = atoi(argv[3]);
    generate_matrix(nx, ny, nz, A, &x, &b, &xexact);
  }
  else
  {
    read_HPC_row(argv[1], A, &x, &b, &xexact);
  }

  bool dump_matrix = false;
  if (dump_matrix && size <= 4)
    dump_matlab_matrix(A, rank);

  int niters = 0;
  double normr = 0.0;
  int max_iter = 100;
  double tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations

  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  double *r = new double[nrow];
  double *p = new double[ncol]; // In parallel case, A is rectangular
  double *Ap = new double[nrow];

  double residual =
      HPCCG_residual(b, x, max_iter, tolerance, niters, normr, r, p, Ap, xexact);

  delete[] p;
  delete[] Ap;
  delete[] r;

  if (rank == 0)
    cout << std::setprecision(5) << "Difference between computed and exact (residual)  = "
         << residual << ".\n"
         << endl;

  return 0;
}
