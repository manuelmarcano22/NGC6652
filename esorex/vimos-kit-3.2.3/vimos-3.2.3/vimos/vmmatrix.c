/* $Id: vmmatrix.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmmatrix.h"
#include "cpl.h"


#define dtiny(a) ((a)<0?(a)> -1.e-30:(a)<1.e-30)


VimosFloatArray *newFloatArray(int len) 
{
  VimosFloatArray *tArray;
  
  tArray = (VimosFloatArray *) cpl_malloc(sizeof(VimosFloatArray));

  /* check if space was allocated */
  if (tArray == NULL) {
    cpl_msg_error("newFloatArray","Allocation Error");
    return(NULL);
  }
  
  tArray->data = (float *) cpl_calloc(len, sizeof(float));
  if (tArray->data == NULL) {
    /* cleanup */
    cpl_free(tArray);
    cpl_msg_error("newFloatArray","Allocation Error");
    return(NULL);
  }
  
  tArray->len = len;
  
  return(tArray);
}


void deleteFloatArray(VimosFloatArray *array) 
{
  if (array == NULL) {
    return;
  }
  cpl_free(array->data);
  cpl_free(array);
}


VimosIntArray *newIntArray(int len) 
{
  VimosIntArray *tArray;
  
  tArray = (VimosIntArray *) cpl_malloc(sizeof(VimosIntArray));

  /* check if space was allocated */
  if (tArray == NULL) {
    cpl_msg_error("newIntArray","Allocation Error");
    return(NULL);
  }
  
  tArray->data = (int *) cpl_calloc(len, sizeof(int));
  if (tArray->data == NULL) {
    /* cleanup */
    cpl_free(tArray);
    cpl_msg_error("newIntArray","Allocation Error");
    return(NULL);
  }
  
  tArray->len = len;
  
  return(tArray);
}


void deleteIntArray(VimosIntArray *array) 
{
  if (array == NULL) {
    return;
  }
  cpl_free(array->data);
  cpl_free(array);
}


/*---------------------------------------------------------------------------
 * Function	:	create_mx()
 * In 		:	number of rows, number of columns
 * Out 		:	Matrix
 * Job		:	allocates memory for a matrix
 * Notice	:	returns Nulmat if cannot allocate
 *--------------------------------------------------------------------------*/
VimosMatrix *newMatrix(int nr, int nc)
{
  VimosMatrix *b;
  b = (VimosMatrix *) cpl_calloc(1,sizeof(VimosMatrix));

  /* check if space was allocated */
  if (b == NULL) {
    cpl_msg_error("newMatrix","Allocation Error");
    return(NULL);
  }

  b->data = (double *) cpl_calloc(nr*nc, sizeof(double));
  if (b->data == NULL) {
    /* cleanup */
    cpl_free(b);
    cpl_msg_error("newMatrix","Allocation Error");
    return(NULL);
  }

  b->nr = nr;
  b->nc = nc;
  return(b);
}

/*---------------------------------------------------------------------------
 * Function	:	copy_mx()
 * In 		:	Matrix
 * Out 		:	Matrix
 * Job		:	copy a matrix
 * Notice	:	returns Nulmat if cannot allocate
 *--------------------------------------------------------------------------*/
VimosMatrix *copyMatrix(VimosMatrix *a)
{
  VimosMatrix *b;
  b = newMatrix(a->nr,a->nc);
  if (b == NULL) {
    cpl_msg_error("copyMatrix","The function newMatrix has returned NULL");
    return(NULL);
  }
  {
    register long s = a->nr*a->nc;
    register double *mm = b->data+s;
    register double *am = a->data+s;
    while (s--) {
      *--mm = *--am;
    }
    
  }
  return(b);
}


/*---------------------------------------------------------------------------
 * Function	:	close_mx()
 * In 		:	Matrix
 * Out 		:	void
 * Job		:	free memory associated to a matrix
 * Notice	:
 *--------------------------------------------------------------------------*/
void deleteMatrix(VimosMatrix *a)
{
  if (a == NULL) { 
    return;
  }
  
  if (a->data != NULL) {
    cpl_free(a->data);
  }
  
  cpl_free(a);
  return;
}


/*---------------------------------------------------------------------------
 * Function	:	mul_mx()
 * In 		:	2 Matrix
 * Out 		:	Matrix
 * Job		:	Multiplies 2 Matrixes
 * Notice	:
 *--------------------------------------------------------------------------*/
VimosMatrix *mulMatrix(VimosMatrix *a, VimosMatrix *b)
{
  VimosMatrix *c, *d;
  int n1 = a->nr;
  int n2 = a->nc;
  int n3 = b->nc;

  register double *a0;
  register double *c0;
  register double *d0;
  register int i,j,k;
  
  if (n2 != b->nr) {
    cpl_msg_error("mulMatrix","Number of row has to be equal to number of column");
    return(NULL);
  }
  c = newMatrix(n1, n3);
  if (c == NULL) {
    cpl_msg_error("mulMatrix","The function newMatrix has returned NULL");
    return(NULL);
  }

  d = transpMatrix(b);
  if (d == NULL) {
    cpl_msg_error("mulMatrix","The function transpMatrix has returned NULL");
    return(NULL);
  }

  for (i = 0, c0 = c->data; i < n1; i++)
    for (j = 0, d0 = d->data; j<n3; j++, c0++)
      for (k = 0, *c0 = 0, a0 = a->data+i*n2; k<n2; k++)
        *c0 += *a0++ * *d0++;

  deleteMatrix(d);
  return(c);
}


/*---------------------------------------------------------------------------
 * Function	:	invert_mx()
 * In 		:	(square) Matrix
 * Out 		:	Matrix
 * Job		:	invert a matrix
 * Notice	:	Hardcoded for 1x1, 2x2, 3x3 matrices
 *--------------------------------------------------------------------------*/
VimosMatrix *invertMatrix(VimosMatrix *aa)
{
  VimosMatrix *bb;
  int test = 1;

  if (aa->nr != aa->nc) {
    cpl_msg_error("invertMatrix","The matrix has to be a square matrix");
    return(NULL);
  }
  
  bb = newMatrix(aa->nr,aa->nc);

  if (bb == NULL) {
    cpl_msg_error("invertMatrix","The function newMatrix has returned NULL");
    return(NULL);
  }
  
  if (aa->nr == 1) {
    double det;
    register double ted;
    det= *(aa->data);
    if (dtiny(det)) {
      test = 0;
    }
    
    ted=1./det;
    *(bb->data)=ted;
  } 
  else {
    if (aa->nr == 2) {
      double det;
      register double ted;
      register double *mm = aa->data;
      double a = *(mm++),b = *(mm++);
      double c = *(mm++),d = *(mm);
      det = a*d-b*c;
      if (dtiny(det)) {
        test = 0;
      }
      
      ted = 1./det;
      mm = bb->data;
      *(mm++) = d*ted,*(mm++) = -b*ted;
      *(mm++) = -c*ted,*(mm) = a*ted; 
    }
    else {
      if (aa->nr == 3) {
        double det;
        register double ted;
        register double *mm = aa->data;
        double a = *(mm++), b = *(mm++), c = *(mm++);
        double d = *(mm++), e = *(mm++), f = *(mm++);
        double g = *(mm++), h = *(mm++), i = *(mm);
        det = a*e*i-a*h*f-b*d*i+b*g*f+c*d*h-c*g*e;
        if (dtiny(det)) {
          test = 0;
        }
        
        ted = 1./det;
        mm = bb->data;
        *(mm++) = (e*i-f*h)*ted, *(mm++) =(c*h-b*i)*ted, *(mm++) = (b*f-e*c)*ted;
        *(mm++) = (f*g-d*i)*ted, *(mm++) =(a*i-g*c)*ted, *(mm++) = (d*c-a*f)*ted;
        *(mm++) = (d*h-g*e)*ted, *(mm++) =(g*b-a*h)*ted, *(mm) = (a*e-d*b)*ted;
      } 
      else {
        VimosMatrix *temp = copyMatrix(aa);
	if (temp == NULL) {
	  cpl_msg_error("invertMatrix","The function copyMatrix has returned NULL");
	  return(NULL);
	}

        if (gaussPivot(temp->data, bb->data, aa->nr) == 0) {
          test = 0;
        }
        deleteMatrix(temp);
      }
    }
  }
  
  if (test == 0) {
    cpl_msg_error("invertMatrix","matrix::invert: not invertible, aborting inversion");
    return(NULL);
  }
  return(bb);
}


/*---------------------------------------------------------------------------
 * Function	:	transp_mx()
 * In 		:	Matrix
 * Out 		:	Matrix
 * Job		:	transposition of a Matrix
 * Notice	:
 *--------------------------------------------------------------------------*/

VimosMatrix *transpMatrix(VimosMatrix *a)
{
  register int nc = a->nc, nr = a->nr;
  register double *a0;
  register double *b0;
  register int i,j;
  VimosMatrix *b = newMatrix(nc,nr);
  
  if (b == NULL) {
    cpl_msg_error("transpMatrix","The function newMatrix has returned NULL");
    return(NULL);
  }
  
  for (i = 0, b0 = b->data; i < nc; i++)
    for (j = 0, a0 = a->data+i; j < nr; j++, a0 += nc, b0++)
      *b0 = *a0;
  return(b);
}


/*---------------------------------------------------------------------------
 * Function	:	gauss_pivot()
 * In 		:	2 matrix lines, number of rows in line
 * Out 		:	error code: 1 if Ok, 0 else
 * Job		:	line simplification with Gauss method
 * Notice	:	should be private to invert_mx()
 *--------------------------------------------------------------------------*/
int gaussPivot(double *ptra, double *ptrc, int n)
     /* c(n,n) = a(n,n)^-1 */
{
  
  register int i,j,k,l;
  int maj;
  double max,r,t;
  double *ptrb;
  
  ptrb = (double *)cpl_calloc(n*n,sizeof(double));

  /* check if space was allocated */
  if (ptrb == NULL) {
    cpl_msg_error("gaussPivot","Allocation Error");
    return(0);
  }

  for (i = 0; i < n; i++) {
    ptrb[i*n+i]= 1.0;
  }
  
  
  for (i = 1; i <= n; i++) {
    /* Search max in current column  */
    max = ABS(*(ptra + n*i-n));
    maj = i;
    for (j = i;j <= n;j++) {
      if (ABS(*(ptra+n*j+i-n-1)) > max) {
        maj = j;
        max = ABS(*(ptra+n*j+i-n-1));
      }
    }
    
    /* swap lines i and maj */
    if (maj != i) {
      for (j = i;j <= n; j++) {
        r = *(ptra+n*maj+j-n-1);
        *(ptra+n*maj+j-n-1) = *(ptra+n*i+j-n-1);
        *(ptra+n*i+j-n-1) = r;
      }
      for(l = 0; l < n; l++) {
        r = *(ptrb+l*n+maj-1);
        *(ptrb+l*n+maj-1) = *(ptrb+l*n+i-1);
        *(ptrb+l*n+i-1) = r;
      }
    }
      
    /* Subtract line by line */
    for (j = i + 1;j <= n;j++) {
      t = (*(ptra+(n+1)*i-n-1));
      if (dtiny(t)) {
        return(0);
      }
      
      r = (*(ptra+n*j+i-n-1)) / t;
      for (l = 0; l < n;l++) {
        *(ptrb+l*n+j-1) -= r * (*(ptrb+l*n+i-1));
      }
      for (k = i; k <= n; k++) {
        *(ptra+n*j+k-n-1) -= r * (*(ptra+n*i+k-n-1));
      }
    }
  }
  
  /* Triangular system resolution	*/
  for (l = 0; l < n; l++) {
    for (i = n; i >= 1; i--) {
      t = (*(ptra+(n+1)*i-n-1));
      if (dtiny(t)) {
        return(0);
      }
      
      *(ptrc+l+(i-1)*n) = (*(ptrb+l*n+i-1)) / t;
      if (i > 1) {
        for (j = i - 1;j > 0;j--) {
          *(ptrb+l*n+j-1) -= (*(ptra+n*j+i-n-1)) * (*(ptrc+l+(i-1)*n));
        }
      }
      
    }
  }
  
  cpl_free(ptrb);
  return(1);

#undef ABS
}



/*---------------------------------------------------------------------------
   Function	:	least_sq_mx()
   In 		:	2 matrices A, B
   Out 		:	1 matrix
   Job		:	from the set XA = B, compute the matrix P defined by
   				P = B.tA.inv(A.tA) 
				P solves X by a least-squares criterion
   Notice	:
 ---------------------------------------------------------------------------*/

VimosMatrix *lsqMatrix(VimosMatrix *A, VimosMatrix *B)
{
  VimosMatrix	*m1, *m2, *m3, *m4, *m5 ;
  
  m1 = transpMatrix(A);
  if (m1 == NULL) {
    cpl_msg_error("lsqMatrix","The function transpMatrix has returned NULL");
    return(NULL);
  }

  m2 = mulMatrix(A, m1);
  if (m2 == NULL) {
    cpl_msg_error("lsqMatrix","The function mulMatrix has returned NULL");
    return(NULL);
  }

  m3 = invertMatrix(m2);
  if (m3 == NULL) {
    cpl_msg_error("lsqMatrix","The function invertMatrix has returned NULL");
    return(NULL);
  }

  m4 = mulMatrix(B, m1);
  if (m4 == NULL) {
    cpl_msg_error("lsqMatrix","The function mulMatrix has returned NULL");
    return(NULL);
  }

  m5 = mulMatrix(m4, m3);
  if (m5 == NULL) {
    cpl_msg_error("lsqMatrix","The function mulMatrix has returned NULL");
    return(NULL);
  }
  
  deleteMatrix(m1);
  deleteMatrix(m2);
  deleteMatrix(m3);
  deleteMatrix(m4);
  
  return(m5);
}

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 These functions are for "compatibility with NR-like routines

------------------------------------------------------------------------------
*/


#define RANGE_END 1
#define FREE_ARG char*

/* 
  allocate a vector of integers with subscript range v[nlow..nhigh]
*/

int *intVector(int nlow, int nhigh)
{
  int *v;
 
  v = (int *)cpl_malloc((size_t) ((nhigh-nlow+1+RANGE_END)*sizeof(int)));
  if (!v)
    {
     abort();
    }

  return v-nlow+RANGE_END;
}

/* 
  free a vector of integers
*/

void freeIntVector(int *v, int nlow, int nhigh)
{
  nhigh = nhigh ;
  cpl_free((FREE_ARG) (v+nlow-RANGE_END));
}

/* 
  allocate a vector of floating points with subscript range v[nlow..nhigh]
*/

float *floatVector(int nlow, int nhigh)
{
  float *v;
 
  v = (float *)cpl_malloc((size_t) ((nhigh-nlow+1+RANGE_END)*sizeof(float)));
  if (!v)
    {
     abort();
    }
  return v-nlow+RANGE_END;
}


/* 
  allocate a vector of doubles with subscript range v[nlow..nhigh]
*/

double *doubleVector(int nlow, int nhigh)
{
  double *v;
 
  v = (double *)cpl_malloc((size_t) ((nhigh-nlow+1+RANGE_END)*sizeof(double)));
  if (!v)
    {
     abort();
    }
  return v-nlow+RANGE_END;
}

/*
  free a vector of floating points
*/

void freeFloatVector(float *v, int nlow, int nhigh)
{
  nhigh = nhigh ;
  cpl_free((FREE_ARG) (v+nlow-RANGE_END));
}

/*
  free a vector of double 
*/

void freeDoubleVector(double *v, int nlow, int nhigh)
{
  nhigh = nhigh ;
  cpl_free((FREE_ARG) (v+nlow-RANGE_END));
}

/*
  allocate a floating point  matrix with subscript range m[nrl..nrh][ncl..nch]
*/

float **Matrix(int nrl, int nrh, int ncl, int nch)
{
  int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  float **m;
 
  /* allocate pointers to rows */
  m = (float **) cpl_malloc((size_t)((nrow+RANGE_END)*sizeof(float*)));
        
  if (!m)
    {
     abort();
    }
  m += RANGE_END;
  m -= nrl;
  /* allocate rows and set pointers to them */
  m[nrl] = (float *) cpl_malloc((size_t)((nrow*ncol+RANGE_END)*sizeof(float)));

  if (!m[nrl])
    {
     abort();
    }

  m[nrl] += RANGE_END;
  m[nrl] -= ncl;

  for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;
 
  /* return pointer to array of pointers to rows */
  return m;
}
 
/*
  free a floating point matrix
*/

void freeMatrix(float **m, int nrl, int nrh, int ncl, int nch)
{
  nrh = nrh;
  nch = nch;
  cpl_free((FREE_ARG) (m[nrl]+ncl-RANGE_END));
  cpl_free((FREE_ARG) (m+nrl-RANGE_END));
}

double **doubleMatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) cpl_malloc((size_t)((nrow+RANGE_END)*sizeof(double*)));
	if (!m) {
	  exit(-1) ;
	}
	m += RANGE_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) cpl_malloc((size_t)((nrow*ncol+RANGE_END)*sizeof(double)));
	if (!m[nrl]) exit(-1);
	m[nrl] += RANGE_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void freeDoubleMatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by doubleMatrix() */
{
  nrh = nrh;
  nch = nch;
  cpl_free((FREE_ARG) (m[nrl]+ncl-RANGE_END));
  cpl_free((FREE_ARG) (m+nrl-RANGE_END));
}


/*
  allocate a floating point matrix m[nrl..nrh][ncl..nch] that points to the 
  matrix declared in the standard C manner as a[nrow][ncol], where 
  nrow=nrh-nrl+1 and ncol=nch-ncl+1. The routine should be called with the 
  address &a[0][0] as the first argument.
*/

float **convertMatrix(float *a, int nrl, int nrh, int ncl, int nch)
{
  int i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m = (float **) cpl_malloc((size_t) ((nrow+RANGE_END)*sizeof(float*)));

  if (!m)
    {
     abort();
    }

  m += RANGE_END;
  m -= nrl;
  /* set pointers to rows */
  m[nrl] = a - ncl;
  for(i=1,j=nrl+1; i<nrow; i++,j++) m[j] = m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

/*
  free a matrix allocated by convertMatrix
*/

void freeConvertedMatrix(float **b, int nrl, int nrh, int ncl, int nch)
{
  nrh = nrh ;
  ncl = ncl ;
  nch = nch ;
  cpl_free((FREE_ARG) (b+nrl-RANGE_END));
}

/* commented out to avoid warning message sin compilation
#undef RANGE_END 1
#undef FREE_ARG char*
*/

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 End of functions for "compatibility with NR-like routines

------------------------------------------------------------------------------
*/


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 These functions define the VimosFloat2DArray type

------------------------------------------------------------------------------
*/

/* 
 Create a VimosFloat2DArray
*/

VimosFloat2DArray *newFloat2DArray(int len) 
{
  VimosFloat2DArray *tArray;
  
  tArray = (VimosFloat2DArray *) cpl_malloc(sizeof(VimosFloat2DArray));
  if (tArray == NULL) {
    exit(-2);
  }
  
  tArray->x = (float *) cpl_calloc(len, sizeof(float));
  tArray->y = (float *) cpl_calloc(len, sizeof(float));
  if ( (tArray->x == NULL) || (tArray->y == NULL) ) {
    cpl_free(tArray);
    exit(-1);
  }
  
  tArray->len = len;
  
  return(tArray);
}


/* 
 Delete a VimosFloat2DArray
*/

void deleteFloat2DArray(VimosFloat2DArray *array) 
{
  if (array == NULL) {
    return;
  }
  cpl_free(array->x);
  cpl_free(array->y);
  cpl_free(array);
}

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 End of functions to define the VimosFloat2DArray type

------------------------------------------------------------------------------
*/




/*

 $Id: vmmatrix.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
 $Author: cgarcia $
 $Date: 2013-03-25 11:43:04 $
 $Revision: 1.2 $

 */
