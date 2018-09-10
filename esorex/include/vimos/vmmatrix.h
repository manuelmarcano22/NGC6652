/* $Id: vmmatrix.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2004 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-03-25 11:43:04 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VM_MATRIX_H
#define VM_MATRIX_H

#include <pilmacros.h>

#include <vmtypes.h>


PIL_BEGIN_DECLS

#define NullMatrix NULL
#define _(b,i,j) (*((b)->data+(i)*(b)->nc+(j))) /* b(i,j)*/


typedef struct _VIMOS_FLOAT_ARRAY_ {
  float *data;
  int len;
} VimosFloatArray;


typedef struct _VIMOS_MATRIX_ {
  double *data;
  int nr;
  int nc;
} VimosMatrix;


VimosFloatArray *newFloatArray(int nels);
void deleteFloatArray(VimosFloatArray *array);


typedef struct _VIMOS_INT_ARRAY_ {
  int *data;
  int len;
} VimosIntArray;


VimosIntArray *newIntArray(int nels);
void deleteIntArray(VimosIntArray *array);

/*---------------------------------------------------------------------------
 * Function	:	create_mx()
 * In 		:	number of rows, number of columns
 * Out 		:	Matrix
 * Job		:	allocates memory for a matrix
 * Notice	:	returns Nulmat if cannot allocate
 *--------------------------------------------------------------------------*/
VimosMatrix *newMatrix(int nr, int nc) ;


/*---------------------------------------------------------------------------
 * Function	:	copy_mx()
 * In 		:	Matrix
 * Out 		:	Matrix
 * Job		:	copy a matrix
 * Notice	:	returns Nulmat if cannot allocate
 *--------------------------------------------------------------------------*/
VimosMatrix *copyMatrix(VimosMatrix *m) ;


/*---------------------------------------------------------------------------
 * Function	:	close_mx()
 * In 		:	Matrix
 * Out 		:	void
 * Job		:	free memory associated to a matrix
 * Notice	:
 *--------------------------------------------------------------------------*/
void deleteMatrix(VimosMatrix *m) ;


/*---------------------------------------------------------------------------
 * Function	:	mul_mx()
 * In 		:	2 Matrix
 * Out 		:	Matrix
 * Job		:	Multiplies 2 Matrixes
 * Notice	:
 *--------------------------------------------------------------------------*/
VimosMatrix *mulMatrix(VimosMatrix *a, VimosMatrix *b) ;


/*---------------------------------------------------------------------------
 * Function	:	invert_mx()
 * In 		:	(square) Matrix
 * Out 		:	Matrix
 * Job		:	invert a matrix
 * Notice	:	Hardcoded for 1x1, 2x2, 3x3 matrices
 *--------------------------------------------------------------------------*/
VimosMatrix *invertMatrix(VimosMatrix *aa) ;


/*---------------------------------------------------------------------------
 * Function	:	transp_mx()
 * In 		:	Matrix
 * Out 		:	Matrix
 * Job		:	transposition of a Matrix
 * Notice	:
 *--------------------------------------------------------------------------*/

VimosMatrix *transpMatrix(VimosMatrix *a) ;


/*---------------------------------------------------------------------------
 * Function	:	gauss_pivot()
 * In 		:	2 matrix lines, number of rows in line
 * Out 		:	error code: 1 if Ok, 0 else
 * Job		:	line simplification with Gauss method
 * Notice	:	should be private to invert_mx()
 *--------------------------------------------------------------------------*/
int gaussPivot(double *ptra, double *ptrc, int n) ;


/*---------------------------------------------------------------------------
   Function :   least_sq_mx()
   In       :   2 matrices A, B
   Out      :   1 matrix
   Job      :   from the set XA = B, compute the matrix P defined by
                P = B.tA.inv(A.tA)
                P solves X by a least-squares criterion
   Notice   :
 ---------------------------------------------------------------------------*/

VimosMatrix *lsqMatrix(VimosMatrix  *A, VimosMatrix *B) ;



/*
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Description:
  Some mathematical functions for "compatibility" with NR-like routines.
  These can handle different ranges for array indexes (not necessarilily
  starting from 0) 

  Updates:
  15 Jun 00: Created (AZ)

------------------------------------------------------------------------------
*/


/* 
  allocate and free a vector of integers with subscript range v[nlow..nhigh]
*/
int *intVector(int nl, int nh);
void freeIntVector(int *v, int nl, int nh);


/* 
  allocate and free a vector of floats with subscript range v[nlow..nhigh]
*/
float *floatVector(int nl, int nh);
void freeFloatVector(float *v, int nl, int nh);

/* 
  allocate and free a vector of doubles with subscript range v[nlow..nhigh]
*/
double *doubleVector(int nl, int nh);
void freeDoubleVector(double *v, int nl, int nh);

/*
  allocate and free a floating point  matrix with subscript range 
  m[nrl..nrh][ncl..nch]
*/
float **Matrix(int nrl, int nrh, int ncl, int nch);
void freeMatrix(float **m, int nrl, int nrh, int ncl, int nch);

/*
  allocate and free a double   matrix with subscript range 
  m[nrl..nrh][ncl..nch]
*/
double **doubleMatrix(long nrl, long nrh, long ncl, long nch);
void freeDoubleMatrix(double **m, long nrl, long nrh, long ncl, long nch);

/*
  allocate a floating point matrix m[nrl..nrh][ncl..nch] that points to the 
  matrix declared in the standard C manner as a[nrow][ncol], where 
  nrow=nrh-nrl+1 and ncol=nch-ncl+1. The routine should be called with the 
  address &a[0][0] as the first argument.
*/
float **convertMatrix(float *a, int nrl, int nrh, int ncl, int nch);
void freeConvertedMatrix(float **b, int nrl, int nrh, int ncl, int nch);




/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 These functions define the VimosFloat2DArray type

   THESE HAVE BEEN ADDED TO VIMOSMATRIX 

------------------------------------------------------------------------------
*/
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  structure VimosFloat2DArray

  Description:
  Structure similar to VimosFloatArray, but 2-dimanesional. It is used
  to compute histograms, for instance.

  Layout:
  float x;
  float y;
  int len;

  Updates:
  01 Aug 00: Created (AZ)

------------------------------------------------------------------------------
*/

typedef struct _VIMOS_FLOAT_2D_ARRAY_ {
  float *x;
  float *y;
  int len;
} VimosFloat2DArray;


VimosFloat2DArray *newFloat2DArray(int len);

void deleteFloat2DArray(VimosFloat2DArray *array);

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 End of functions to define the VimosFloat2DArray type

------------------------------------------------------------------------------
*/

PIL_END_DECLS

#endif /* VM_MATRIX_H */
