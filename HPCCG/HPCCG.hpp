
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

#ifndef HPCCG_H
#define HPCCG_H
#include "HPC_sparsemv.hpp"
#include "ddot.hpp"
#include "waxpby.hpp"
#include "HPC_Sparse_Matrix.hpp"
#include "HPCCG.hpp"

#include "Derivative.hpp"

#include "adapt.h"
#include "adapt-impl.cpp"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <cmath>

// this function will compute the Conjugate Gradient...
// A <=> Matrix
// b <=> constant
// xnot <=> initial guess
// max_iter <=> how many times we iterate
// tolerance <=> specifies how "good"of a value we would like
// x <=> used for return value

// A is known
// x is unknown vector
// b is known vector
// xnot = 0
// niters is the number of iterations

namespace adapt {

int HPCCG(HPC_Sparse_Matrix *A,
          const AD_real *const b, AD_real *const x,
          const int max_iter, const double tolerance, int &niters,
          AD_real &normr, AD_real *r, AD_real *p, AD_real *Ap)

{
  int nrow = A->local_nrow;
  int ncol = A->local_ncol;

  normr = 0.0;
  AD_real rtrans = 0.0;
  AD_real oldrtrans = 0.0;

  int rank = 0; // Serial case (not using MPI)

  int print_freq = max_iter / 10;
  if (print_freq > 50)
    print_freq = 50;
  if (print_freq < 1)
    print_freq = 1;

  // p is of length ncols, copy x to p for sparse MV operation
  waxpby(nrow, 1.0, x, 0.0, x, p);

  HPC_sparsemv(A, p, Ap);
  
  waxpby(nrow, 1.0, b, -1.0, Ap, r);
  
  ddot(nrow, r, r, &rtrans);
  
  normr = sqrt(AD_value(rtrans));
  AD_INTERMEDIATE(normr, "normr");

  if (rank == 0)
    cout << "Initial Residual = " << AD_value(normr) << endl;

  for (int k = 1; k < max_iter && normr > tolerance; k++)
  {
    if (k == 1)
    {
      waxpby(nrow, 1.0, r, 0.0, r, p);
    }
    else
    {
      oldrtrans = rtrans;
	  
      ddot(nrow, r, r, &rtrans);
	  
      AD_real beta = rtrans / oldrtrans;
      AD_INTERMEDIATE(beta, "beta");
	  
      waxpby(nrow, 1.0, r, beta, p, p);
	  
      // for(int ii=0; ii < nrow; ii++) {AD_INTERMEDIATE(p[ii], "p"+std::to_string(k));};
    }
    normr = sqrt(rtrans);
    AD_INTERMEDIATE(normr, "normr");
    if (rank == 0 && (k % print_freq == 0 || k + 1 == max_iter))
      cout << "Iteration = " << k << "   Residual = " << AD_value(normr) << endl;

    HPC_sparsemv(A, p, Ap);

    for (int ii = 0; ii < nrow; ii++)
    {
      AD_INTERMEDIATE(Ap[ii], "Ap");
    };
    AD_real alpha = 0.0;
	
    ddot(nrow, p, Ap, &alpha);
	
    AD_INTERMEDIATE(alpha, "alpha");
    alpha = rtrans / alpha;
    AD_INTERMEDIATE(alpha, "alpha");
	
    waxpby(nrow, 1.0, x, alpha, p, x); // 2*nrow ops
    for (int ii = 0; ii < nrow; ii++)
    {
      AD_INTERMEDIATE(x[ii], "x");
    };
    waxpby(nrow, 1.0, r, -alpha, Ap, r);
	
    for (int ii = 0; ii < nrow; ii++)
    {
      AD_INTERMEDIATE(r[ii], "r");
    };
    niters = k;
  }

  return (0);
}

} // namespace adapt


namespace clad {

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
  precision sum, temp;

  ////// HPCCG
  int nrow = A.local_nrow;
  int ncol = A.local_ncol;

  normr = 0.0;
  precision rtrans = 0.0;
  precision oldrtrans = 0.0;

  precision beta = 0.0;

  for (int i = 0; i < nrow; i++) {
    temp = x[i];
    p[i] = temp + beta * temp;
  }

  // HPC_sparsemv(A, p, Ap);
  for (int i = 0; i < nrow; i++)
  {
    cur_nnz = A.nnz_in_row[i];
    sum = 0;
    for (int j = 0; j < cur_nnz; j++)
      sum += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
    Ap[i] = sum;
  }

  beta = -1.0;
  for (int i = 0; i < nrow; i++)
    r[i] = b[i] + beta * Ap[i];

  // ddot(nrow, r, r, &rtrans);
  rtrans = 0.0;
  for (int i = 0; i < nrow; i++) {
    temp = r[i];
    rtrans += temp * temp;
  }
  // ddot end

  normr = sqrt(rtrans);

  for (int k = 1; k < max_iter && normr > tolerance; k++)
  {
    if (k == 1)
    {
      //// waxpby(nrow, 1.0, r, 0.0, r, p);
      beta = 0.0;
      for (int i = 0; i < nrow; i++) {
        temp = r[i];
        p[i] = temp + beta * temp;
      }
    }
    else
    {
      oldrtrans = rtrans;
      // ddot (nrow, r, r, &rtrans); // 2*nrow ops
      rtrans = 0.0;
      for (int i = 0; i < nrow; i++) {
        temp = r[i];
        rtrans += temp * temp;
      }
      // ddot end

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
      sum = 0;
      for (int j = 0; j < cur_nnz; j++)
        sum += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
      Ap[i] = sum;
      
    }

    // waxpby(nrow, 1.0, x, alpha, p, x);
    // beta = 0.0;
    // beta = ddot(nrow, p, Ap);
    // ddot(nrow, p, Ap, &alpha); // 2*nrow ops
    beta = 0.0;
    for (int i = 0; i < nrow; i++)
      beta += p[i] * Ap[i];
    // ddot end
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

void executeGradient(int nrow, int ncol, double* x, double* b, double* xexact) {
  // auto df = clad::gradient(HPCCG_residual<double>);
  // auto df = clad::estimate_error(HPCCG_residual<double>);


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

} // clad

#endif
