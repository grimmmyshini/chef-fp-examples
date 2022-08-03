
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

#include "clad/Differentiator/Differentiator.h"
#include "../PrintModel/ErrorFunc.h"

#include "generate_matrix.hpp"
#include "read_HPC_row.hpp"
#include "HPC_Sparse_Matrix.hpp"
#include "dump_matlab_matrix.hpp"
#include "ddot.hpp"

#undef DEBUG

HPC_Sparse_Matrix A;

double HPCCG_residual(double *b, double *x,
                      double *xexact, double *r, double *p,
                      double *Ap);

template <typename precision>
double HPCCG_residual(double *b, double *x,
                      double *xexact, precision *r, precision *p,
                      precision *Ap)
{
  int niters = 0;
  precision normr = 0.0;
  int max_iter = 100;
  precision tolerance = 0.0; // Set tolerance to zero to make all runs do max_iter iterations
  int cur_nnz;

  ////// HPCCG
  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  normr = 0.0;
  precision rtrans = 0.0;
  precision oldrtrans = 0.0;

  precision beta = 0.0;

  for (int i = 0; i < nrow; i++)
    p[i] = x[i] + beta * x[i];

  // HPC_sparsemv(A, p, Ap);
  for (int i = 0; i < nrow; i++)
  {
    cur_nnz = A.nnz_in_row[i];
    Ap[i] = 0;
    for (int j = 0; j < cur_nnz; j++)
      Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
  }

  beta = -1.0;
  for (int i = 0; i < nrow; i++)
    r[i] = b[i] + beta * Ap[i];

  rtrans = ddot(nrow, r, r);
  normr = sqrt(rtrans);

  for (int k = 1; k < max_iter && normr > tolerance; k++)
  {
    if (k == 1)
    {
      //// waxpby(nrow, 1.0, r, 0.0, r, p);
      beta = 0.0;
      for (int i = 0; i < nrow; i++)
        p[i] = r[i] + beta * r[i];
    }
    else
    {
      oldrtrans = rtrans;
      rtrans = ddot(nrow, r, r);
      beta = rtrans / oldrtrans;

      //// waxpby(nrow, 1.0, r, beta, p, p);
      for (int i = 0; i < nrow; i++)
        p[i] = r[i] + beta * p[i];
    }
    normr = sqrt(rtrans);

    //// HPC_sparsemv(A, p, Ap);
    for (int i = 0; i < nrow; i++)
    {
      cur_nnz = A.nnz_in_row[i];
      Ap[i] = 0;
      for (int j = 0; j < cur_nnz; j++)
        Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
      
    }

    // waxpby(nrow, 1.0, x, alpha, p, x);
    // beta = 0.0;
    beta = ddot(nrow, p, Ap);
    beta = rtrans / beta;

    for (int i = 0; i < nrow; i++)
      x[i] = x[i] + beta * p[i];

    // waxpby(nrow, 1.0, r, -alpha, Ap, r);
    beta = -beta;
    for (int i = 0; i < nrow; i++)
      r[i] = r[i] + beta * Ap[i];

    niters = k;
  }

  // Compute difference between known exact solution and computed solution
  // All processors are needed here.

  precision residual = 0.0;
  precision diff = 0;

  for (int i = 0; i < nrow; i++)
  {
    diff = fabs(x[i] - xexact[i]);
    if(diff > residual) residual = diff;
  }

  return residual;
}

template <typename T>
void printVals(T* arr, int n) {
  for(int i = 0; i < n; i++ )
    cout << arr[i] << " ";
  cout << endl;
}

#include "Derivative.h"

void executeGradient(int nrow, int ncol, double* x, double* b, double* xexact) {
  // auto df = clad::gradient(HPCCG_residual<double>);
  auto df = clad::estimate_error(HPCCG_residual<double>);


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

  HPCCG_residual_grad(b, x, xexact, r, p, Ap, d_b, d_x, d_xexact, d_r, d_p, d_Ap, _final_error);

  cout << "\nFinal error in HPCCG =" << _final_error << endl;

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
}

template <typename precision = double>
precision executefunction(int nrow, int ncol, double* x, double* b, double* xexact) {
  precision *r = new precision[nrow];
  precision *p = new precision[ncol]; // In parallel case, A is rectangular
  precision *Ap = new precision[nrow];

  precision residual = HPCCG_residual<precision>(b, x, xexact, r, p, Ap);

  delete[] p;
  delete[] Ap;
  delete[] r;

  return residual;
} 

int main(int argc, char *argv[])
{

  double *x, *b, *xexact;
  double *x_diff, *b_diff, *xexact_diff;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;
  double times[7];
  double t6 = 0.0;
  int nx, ny, nz;

  if (argc != 2 && argc != 4)
  {
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
  if (dump_matrix)
    dump_matlab_matrix(A);

  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  
  cout << "Actual Error: " << executefunction<double>(nrow, ncol, x, b, xexact);

  // executeGradient(nrow, ncol, x, b, xexact);

  delete[] x;
  delete[] xexact;
  delete[] b;

  return 0;
}
