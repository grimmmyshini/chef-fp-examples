
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

#ifndef HPCCG_ADAPT_H
#define HPCCG_ADAPT_H
#include "HPC_sparsemv.hpp"
#include "ddot.hpp"
#include "waxpby.hpp"
#include "HPC_Sparse_Matrix.hpp"

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

#endif
