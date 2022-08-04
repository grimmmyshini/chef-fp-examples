
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
#include "HPC_sparsemv.hpp"
#include "compute_residual.hpp"
#include "HPCCG.hpp"
#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"

#undef DEBUG

#include "adapt.h"

int main(int argc, char *argv[])
{
  HPC_Sparse_Matrix *A;
  AD_real *x, *b;
  double *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;
  double t6 = 0.0;
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
    adapt::generate_matrix(nx, ny, nz, &A, &x, &b, &xexact);
  }
  else
  {
    adapt::read_HPC_row(argv[1], &A, &x, &b, &xexact);
  }

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

  delete[] p;
  delete[] Ap;
  delete[] r;
  return 0;
}
