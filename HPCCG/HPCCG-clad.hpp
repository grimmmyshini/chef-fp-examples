
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

#ifndef HPCCG_CLAD_H
#define HPCCG_CLAD_H
#include "HPC_sparsemv.hpp"
#include "ddot.hpp"
#include "waxpby.hpp"
#include "HPC_Sparse_Matrix.hpp"

#include "Derivative.hpp"

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


namespace clad {

double HPCCG(double *b, double *x, int max_iter, double tolerance,
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
    for (int j = 0; j < A.nnz_in_row[i]; j++)
      Ap[i] += A.ptr_to_vals_in_row[i][j] * p[A.ptr_to_inds_in_row[i][j]];
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
    for (int i = 0; i < nrow; i++)
      rtrans += r[i] * r[i];
  // ddot end

  normr = sqrt(rtrans);

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

    // HPC_sparsemv(A, p, Ap); // 2*nnz ops
    for (int i = 0; i < nrow; i++)
    {
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


} // clad

#endif
