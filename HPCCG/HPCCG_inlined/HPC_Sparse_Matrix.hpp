
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

#ifndef HPC_SPARSE_MATRIX_H
#define HPC_SPARSE_MATRIX_H

// These constants are upper bounds that might need to be changes for 
// pathological matrices, e.g., those with nearly dense rows/columns.

const int max_external = 100000;
const int max_num_messages = 500;
const int max_num_neighbors = max_num_messages;


struct HPC_Sparse_Matrix_STRUCT {
  char   *title;
  int start_row;
  int stop_row;
  int total_nrow;
  long long total_nnz;
  int local_nrow;
  int local_ncol;  // Must be defined in make_local_matrix
  int local_nnz;
  int  * nnz_in_row;
  double ** ptr_to_vals_in_row;
  int ** ptr_to_inds_in_row;
  double ** ptr_to_diags;

  double *list_of_vals;   //needed for cleaning up memory
  int *list_of_inds;      //needed for cleaning up memory

};
typedef struct HPC_Sparse_Matrix_STRUCT HPC_Sparse_Matrix;


void destroyMatrix(HPC_Sparse_Matrix * &A)
{
  if(A->title)
  {
    delete [] A->title;
  }
  if(A->nnz_in_row)
  {
    delete [] A->nnz_in_row;
  }
  if(A->list_of_vals)
  {
    delete [] A->list_of_vals;
  }
  if(A->ptr_to_vals_in_row !=0)
  {
    delete [] A->ptr_to_vals_in_row;
  }
  if(A->list_of_inds)
  {
    delete [] A->list_of_inds;
  }
  if(A->ptr_to_inds_in_row !=0)
  {
    delete [] A->ptr_to_inds_in_row;
  }
  if(A->ptr_to_diags)
  {
    delete [] A->ptr_to_diags;
  }

  delete A;
  A = 0;
}

#endif

