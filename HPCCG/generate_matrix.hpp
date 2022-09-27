
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

#ifndef GENERATE_MATRIX_H
#define GENERATE_MATRIX_H

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cstdlib>
#include <cstdio>
#include <cassert>

#include "HPC_Sparse_Matrix.hpp"
#include "adapt.h"

namespace adapt {

void generate_matrix(int nx, int ny, int nz, HPC_Sparse_Matrix **A, double **x, double **b, double **xexact)

{
  int debug = 0;
  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it
  (*A)->title = 0;


  // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
  bool use_7pt_stencil = false;

  int local_nrow = nx*ny*nz; // This is the size of our subblock
  assert(local_nrow>0); // Must have something to work with
  int local_nnz = 27*local_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)

  int total_nrow = local_nrow*size; // Total number of grid points in mesh
  long long total_nnz = 27* (long long) total_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)

  int start_row = local_nrow*rank; // Each processor gets a section of a chimney stack domain
  int stop_row = start_row+local_nrow-1;
  

  // Allocate arrays that are of length local_nrow
  (*A)->nnz_in_row = new int[local_nrow];
  (*A)->ptr_to_vals_in_row = new double*[local_nrow];
  (*A)->ptr_to_vals_in_row_f = new float*[local_nrow];
  (*A)->ptr_to_inds_in_row = new int   *[local_nrow];
  (*A)->ptr_to_diags       = new double*[local_nrow];
  (*A)->ptr_to_diags_f       = new float*[local_nrow];

  *x = new double[local_nrow];
  *b = new double[local_nrow];
  *xexact = new double[local_nrow];


  // Allocate arrays that are of length local_nnz
  (*A)->list_of_vals = new double[local_nnz];
  (*A)->list_of_inds = new int   [local_nnz];
  (*A)->list_of_vals_f = new float[local_nnz];

  double * curvalptr = (*A)->list_of_vals;
  int * curindptr = (*A)->list_of_inds;
  float * curvalptr_f = (*A)->list_of_vals_f;

  long long nnzglobal = 0;
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int curlocalrow = iz*nx*ny+iy*nx+ix;
	int currow = start_row+iz*nx*ny+iy*nx+ix;
	int nnzrow = 0;
  (*A)->ptr_to_vals_in_row_f[curlocalrow] = curvalptr_f;
	(*A)->ptr_to_vals_in_row[curlocalrow] = curvalptr;
	(*A)->ptr_to_inds_in_row[curlocalrow] = curindptr;
	for (int sz=-1; sz<=1; sz++) {
	  for (int sy=-1; sy<=1; sy++) {
	    for (int sx=-1; sx<=1; sx++) {
	      int curcol = currow+sz*nx*ny+sy*nx+sx;
//            Since we have a stack of nx by ny by nz domains , stacking in the z direction, we check to see
//            if sx and sy are reaching outside of the domain, while the check for the curcol being valid
//            is sufficient to check the z values
              if ((ix+sx>=0) && (ix+sx<nx) && (iy+sy>=0) && (iy+sy<ny) && (curcol>=0 && curcol<total_nrow)) {
                if (!use_7pt_stencil || (sz*sz+sy*sy+sx*sx<=1)) { // This logic will skip over point that are not part of a 7-pt stencil
                  if (curcol==currow) {
		    (*A)->ptr_to_diags[curlocalrow] = curvalptr;
		    (*A)->ptr_to_diags_f[curlocalrow] = curvalptr_f;
		    *curvalptr++ = 27.0;
		    *curvalptr_f++ = 27.0;
		  }
		  else {
		    *curvalptr++ = -1.0;
        *curvalptr_f++ = -1.0;
                  }
		  *curindptr++ = curcol;
		  nnzrow++;
	        } 
              }
	    } // end sx loop
          } // end sy loop
        } // end sz loop
	(*A)->nnz_in_row[curlocalrow] = nnzrow;
	nnzglobal += nnzrow;
	(*x)[curlocalrow] = 0.0;
	(*b)[curlocalrow] = 27.0 - ((double) (nnzrow-1));
	(*xexact)[curlocalrow] = 1.0;
      } // end ix loop
     } // end iy loop
  } // end iz loop  
  if (debug) cout << "Process "<<rank<<" of "<<size<<" has "<<local_nrow;
  
  if (debug) cout << " rows. Global rows "<< start_row
		  <<" through "<< stop_row <<endl;
  
  if (debug) cout << "Process "<<rank<<" of "<<size
		  <<" has "<<local_nnz<<" nonzeros."<<endl;

  (*A)->start_row = start_row ; 
  (*A)->stop_row = stop_row;
  (*A)->total_nrow = total_nrow;
  (*A)->total_nnz = total_nnz;
  (*A)->local_nrow = local_nrow;
  (*A)->local_ncol = local_nrow;
  (*A)->local_nnz = local_nnz;

  return;
}

void generate_matrix(int nx, int ny, int nz, HPC_Sparse_Matrix **A, AD_real **x, AD_real **b, double **xexact)

{
  int debug = 0;
  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it
  (*A)->title = 0;


  // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
  bool use_7pt_stencil = false;

  int local_nrow = nx*ny*nz; // This is the size of our subblock
  assert(local_nrow>0); // Must have something to work with
  int local_nnz = 27*local_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)

  int total_nrow = local_nrow*size; // Total number of grid points in mesh
  long long total_nnz = 27* (long long) total_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)

  int start_row = local_nrow*rank; // Each processor gets a section of a chimney stack domain
  int stop_row = start_row+local_nrow-1;
  

  // Allocate arrays that are of length local_nrow
  (*A)->nnz_in_row = new int[local_nrow];
  (*A)->ptr_to_vals_in_row = new double*[local_nrow];
  (*A)->ptr_to_inds_in_row = new int   *[local_nrow];
  (*A)->ptr_to_diags       = new double*[local_nrow];

  *x = new AD_real[local_nrow];
  *b = new AD_real[local_nrow];
  *xexact = new double[local_nrow];


  // Allocate arrays that are of length local_nnz
  (*A)->list_of_vals = new double[local_nnz];
  (*A)->list_of_inds = new int   [local_nnz];

  double * curvalptr = (*A)->list_of_vals;
  int * curindptr = (*A)->list_of_inds;

  long long nnzglobal = 0;
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int curlocalrow = iz*nx*ny+iy*nx+ix;
	int currow = start_row+iz*nx*ny+iy*nx+ix;
	int nnzrow = 0;
	(*A)->ptr_to_vals_in_row[curlocalrow] = curvalptr;
	(*A)->ptr_to_inds_in_row[curlocalrow] = curindptr;
	for (int sz=-1; sz<=1; sz++) {
	  for (int sy=-1; sy<=1; sy++) {
	    for (int sx=-1; sx<=1; sx++) {
	      int curcol = currow+sz*nx*ny+sy*nx+sx;
//            Since we have a stack of nx by ny by nz domains , stacking in the z direction, we check to see
//            if sx and sy are reaching outside of the domain, while the check for the curcol being valid
//            is sufficient to check the z values
              if ((ix+sx>=0) && (ix+sx<nx) && (iy+sy>=0) && (iy+sy<ny) && (curcol>=0 && curcol<total_nrow)) {
                if (!use_7pt_stencil || (sz*sz+sy*sy+sx*sx<=1)) { // This logic will skip over point that are not part of a 7-pt stencil
                  if (curcol==currow) {
		    (*A)->ptr_to_diags[curlocalrow] = curvalptr;
		    *curvalptr++ = 27.0;
		  }
		  else {
		    *curvalptr++ = -1.0;
                  }
		  *curindptr++ = curcol;
		  nnzrow++;
	        } 
              }
	    } // end sx loop
          } // end sy loop
        } // end sz loop
	(*A)->nnz_in_row[curlocalrow] = nnzrow;
	nnzglobal += nnzrow;
	(*x)[curlocalrow] = 0.0;
	(*b)[curlocalrow] = 27.0 - ((double) (nnzrow-1));
	(*xexact)[curlocalrow] = 1.0;
      } // end ix loop
     } // end iy loop
  } // end iz loop  
  if (debug) cout << "Process "<<rank<<" of "<<size<<" has "<<local_nrow;
  
  if (debug) cout << " rows. Global rows "<< start_row
		  <<" through "<< stop_row <<endl;
  
  if (debug) cout << "Process "<<rank<<" of "<<size
		  <<" has "<<local_nnz<<" nonzeros."<<endl;

  (*A)->start_row = start_row ; 
  (*A)->stop_row = stop_row;
  (*A)->total_nrow = total_nrow;
  (*A)->total_nnz = total_nnz;
  (*A)->local_nrow = local_nrow;
  (*A)->local_ncol = local_nrow;
  (*A)->local_nnz = local_nnz;

  return;
}

} // adapt

void generate_matrix(int nx, int ny, int nz, HPC_Sparse_Matrix& A, double** x,
                     double** b, double** xexact)

{
  int debug = 0;

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  // *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it
  A.title = 0;

  // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
  bool use_7pt_stencil = false;

  int local_nrow = nx * ny * nz;   // This is the size of our subblock
  assert(local_nrow > 0);          // Must have something to work with
  int local_nnz = 27 * local_nrow; // Approximately 27 nonzeros per row (except
                                   // for boundary nodes)

  int total_nrow = local_nrow * size; // Total number of grid points in mesh
  long long total_nnz =
      27 * (long long)total_nrow; // Approximately 27 nonzeros per row (except
                                  // for boundary nodes)

  int start_row =
      local_nrow *
      rank; // Each processor gets a section of a chimney stack domain
  int stop_row = start_row + local_nrow - 1;

  // Allocate arrays that are of length local_nrow
  A.nnz_in_row = new int[local_nrow];
  A.ptr_to_vals_in_row = new double*[local_nrow];
  A.ptr_to_vals_in_row_f = new float*[local_nrow];
  A.ptr_to_inds_in_row = new int*[local_nrow];
  A.ptr_to_diags = new double*[local_nrow];
  A.ptr_to_diags_f = new float*[local_nrow];

  *x = new double[local_nrow];
  *b = new double[local_nrow];
  *xexact = new double[local_nrow];

  // Allocate arrays that are of length local_nnz
  A.list_of_vals = new double[local_nnz];
  A.list_of_vals_f = new float[local_nnz];
  A.list_of_inds = new int[local_nnz];

  double* curvalptr = A.list_of_vals;
  float* curvalptr_f = A.list_of_vals_f;
  int* curindptr = A.list_of_inds;

  long long nnzglobal = 0;
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int curlocalrow = iz * nx * ny + iy * nx + ix;
        int currow = start_row + iz * nx * ny + iy * nx + ix;
        int nnzrow = 0;
        A.ptr_to_vals_in_row[curlocalrow] = curvalptr;
        A.ptr_to_vals_in_row_f[curlocalrow] = curvalptr_f;
        A.ptr_to_inds_in_row[curlocalrow] = curindptr;
        for (int sz = -1; sz <= 1; sz++) {
          for (int sy = -1; sy <= 1; sy++) {
            for (int sx = -1; sx <= 1; sx++) {
              int curcol = currow + sz * nx * ny + sy * nx + sx;
              //            Since we have a stack of nx by ny by nz domains ,
              //            stacking in the z direction, we check to see if sx
              //            and sy are reaching outside of the domain, while the
              //            check for the curcol being valid is sufficient to
              //            check the z values
              if ((ix + sx >= 0) && (ix + sx < nx) && (iy + sy >= 0) &&
                  (iy + sy < ny) && (curcol >= 0 && curcol < total_nrow)) {
                if (!use_7pt_stencil ||
                    (sz * sz + sy * sy + sx * sx <=
                     1)) { // This logic will skip over point that are not part
                           // of a 7-pt stencil
                  if (curcol == currow) {
                    A.ptr_to_diags[curlocalrow] = curvalptr;
                    A.ptr_to_diags_f[curlocalrow] = curvalptr_f;
                    *curvalptr++ = 27.0;
                    *curvalptr_f++ = 27.0;
                  } else {
                    *curvalptr++ = -1.0;
                    *curvalptr_f++ = -1.0;
                  }
                  *curindptr++ = curcol;
                  nnzrow++;
                }
              }
            } // end sx loop
          }   // end sy loop
        }     // end sz loop
        A.nnz_in_row[curlocalrow] = nnzrow;
        nnzglobal += nnzrow;
        (*x)[curlocalrow] = 0.0;
        (*b)[curlocalrow] = 27.0 - ((double)(nnzrow - 1));
        (*xexact)[curlocalrow] = 1.0;
      } // end ix loop
    }   // end iy loop
  }     // end iz loop
  if (debug)
    cout << "Process " << rank << " of " << size << " has " << local_nrow;

  if (debug)
    cout << " rows. Global rows " << start_row << " through " << stop_row
         << endl;

  if (debug)
    cout << "Process " << rank << " of " << size << " has " << local_nnz
         << " nonzeros." << endl;

  A.start_row = start_row;
  A.stop_row = stop_row;
  A.total_nrow = total_nrow;
  A.total_nnz = total_nnz;
  A.local_nrow = local_nrow;
  A.local_ncol = local_nrow;
  A.local_nnz = local_nnz;

  return;
}

#endif
