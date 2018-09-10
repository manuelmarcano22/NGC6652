/* $Id: vmfit.c,v 1.3 2013-08-07 18:30:28 cgarcia Exp $
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
 * $Date: 2013-08-07 18:30:28 $
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilstrutils.h>
#include <pilerrno.h>

#include "vmtypes.h"
#include "vmmatrix.h"
#include "vmmath.h"
#include "vmfit.h"
#include "cpl.h"

/**
 * @name vimosFit
 *
 * @doc
 *   The module vimosFit collects functions
 *   for fitting in the DRS.
 */

/*@{*/

static VimosBool robustLinearFit(double*,double*,int,double*,double*,double*) ;

/*
---------------------------------------------------------------------------
   Function	:  fit_1d_poly()
   In 		:  requested polynomial degree
                   a list of pixel positions + number of pixels in the list
                   (out) mean squared error, set to NULL if you do not want
                   to compute it.
   Out 		:  newly allocated array containing fit coefficients
   Job		:  fit a polynomial to a list of pixel positions
   Notice	:  The fitted polynomial is such that:
                   y = c[0] + c[1].x + c[2].x^2 + ... + c[n].x^n
                   So requesting a polynomial of degree n will return n+1
                   coefficients. Beware that with such polynomials, two
                   input points shall never be on the same vertical!
---------------------------------------------------------------------------
*/

double *fit1DPoly(int polyDeg, VimosDpoint *list, int np, 
                  double *meanSquaredError)
{
  const char    modName[] = "fit1DPoly";
  int           i, k;
  VimosMatrix   *mA, *mB, *mX;
  double        *c;
  double        err;
  double        xp, y;

  if (np < polyDeg+1) {
    cpl_msg_debug(modName, 
      "The number of pixel in the list is less then polynomial degree");
    return NULL;
  }
  
  mA = newMatrix(polyDeg+1, np) ;
  if (mA == NULL) {
    cpl_msg_debug(modName, "The function newMatrix has returned NULL");
    return(NULL);
  }

  mB = newMatrix(1, np) ;
  
  if (mB == NULL) {
    cpl_msg_debug(modName, "The function newMatrix has returned NULL");
    return(NULL);
  }

/*  c = cpl_malloc((polyDeg+1)*sizeof(double)) ; */
/* */
/* check if space was allocated */
/*  if (c == NULL) { */
/*    cpl_msg_debug(modName, "Allocation Error"); */
/*    return(NULL); */
/*  } */
/* */
/*  for (i = 0; i < (polyDeg+1); i++) { */
/*    c[i] = 0.0; */
/*  } */

  for (i = 0; i < np; i++) {
    mA->data[i] = 1.0;
    for (k = 1; k <= polyDeg; k++) {
      mA->data[i+k*np] = ipow(list[i].x, k);
    }
    mB->data[i] = list[i].y;
  }
  
  /*
   * Solve XA=B by a least-square solution (aka pseudo-inverse).
   */
  mX = lsqMatrix(mA, mB);
  /*
   * Delete input matrices
   */
  deleteMatrix(mA);
  deleteMatrix(mB);
  /*
   * Examine result
   */
  if (mX == NULL) {
    cpl_msg_debug(modName, "The function lsqMatrix has returned NULL");
    return(NULL);
  }
  
  c = cpl_malloc((polyDeg+1)*sizeof(double)) ;

  /* check if space was allocated */
  if (c == NULL) {
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }

  for (i = 0; i < (polyDeg+1); i++) {
    c[i] = mX->data[i];
  }
  deleteMatrix(mX);
  
  /*
   * If requested, compute mean squared error
   */
  if (meanSquaredError != NULL) {
    err = 0.00 ;
    for (i = 0 ; i < np ; i++) {
      y = c[0] ;
      /*
       * Compute the value obtained through the fit
       */
      for (k = 1; k <= polyDeg; k++) {
        xp = ipow(list[i].x, k) ;
        y += c[k] * xp ; 
      }
      /*
       * Subtract from the true value, square, accumulate
       */
      xp   = ipow(list[i].y - y, 2) ;
      err += xp ; 
    }
    /* Average the error term */
    err /= (double)np ;
    *meanSquaredError = err ;
  }
  return(c);
}



/*---------------------------------------------------------------------------
   Function	:   fit_surface_polynomial()
   In 		:   list of pixels, # of pixels in the list.
                    character string indicating which coefficients should be
                    taken into account, maximum polynomial degree.
   Out 		:   double * (table of fitted coefficients)
                    number of coefficients returned
                    mean squared error for the fit.
   Job		:   fit a 2d surface with a polynomial in (x,y).
   Notice	:   To define which coefficients should be computed, either
		    provide NULL for the control string and the maximal
                    polynomial degree, or fill up the control string as
                    follows:

                    The control string contains (int,int) couples. The first
                    integer specifies the degree for X, the second one the
                    degree for Y. Couples are given in parentheses, integers
                    separated by a comma, with no blanks within the
                    parentheses. Couples are separated from other couples by
                    one blank character. Example: to compute the fit for an
                    equation of type:
                    
                    P(x,y) = c[0] + c[1].x + c[2].x^2 + c[3].x.y
                    
                    You would provide the following control string:
                    
                    "(0,0) (1,0) (2,0) (1,1)"
                    (0,0) is degx=0 and degy=0 -> constant term	c[0]
                    (1,0) is degx=1 and degy=0 -> term in x 	c[1]
                    (2,0) is degx=2 and degy=0 -> term in x^2 	c[2]
                    (1,1) is degx=1 and degy=1 -> term in x.y 	c[3]
                    
                    The maximal polynomial degree indicates the highest sum
                    for X and Y degrees. Example: for poly_deg=3, only the
                    following terms can be computed:
                    
                    1		x		x^2		x^3
                    y		x.y		x^2.y
                    y^2		x.y^2
                    y^3
                    
                    If you do not provide any control string, use NULL as
                    argument and set the polynomial degree to what you wish,
                    all terms satisfying the condition (degx+degy<=polydeg)
                    will be taken into account for the fit.
				
 ---------------------------------------------------------------------------*/

double *fitSurfacePolynomial(VimosPixel *surface, int np, char *controlString,
	int polyDeg, int *ncoeffs, double *meanSquaredError)
{
  const char  modName[] = "fitSurfacePolynomial";
  int         i, j;
  int         degx, degy;
  VimosMatrix *mA, *mB, *mX;
  double      x, y, z;
  int         nc;
  double      *c;
  double      err;
  int         *degxTab;
  int         *degyTab;
  
 /*
  * Fill up look-up table for coefficients to compute
  */
  
  nc = (1+polyDeg)*(2+polyDeg) / 2 ;
  degxTab = cpl_malloc(nc * sizeof(int)) ;

 /* check if space was allocated */
  if (degxTab == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  degyTab = cpl_malloc(nc * sizeof(int)) ;
  
  /* check if space was allocated */
  if (degyTab == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  if (controlString == NULL) {
    i = 0 ;
    for (degy = 0; degy <= polyDeg; degy++) {
      for (degx = 0; degx <= polyDeg; degx++) {
        if (degx+degy <= polyDeg) {
          degxTab[i] = degx ;
          degyTab[i] = degy ;
          i++ ;
        }
      }
    }
  } 
  else {
    nc = buildupPolytabFromString(controlString, polyDeg, degxTab, degyTab) ;
    if (nc == -1) {
      cpl_msg_error(modName, "function buildupPolytabFromString returned error");
      return(NULL);
    }  
  }
  
 /*
  * Initialize matrices
  * mA contains the polynomial terms in the order described
  * above in each column, for each input point.
  * mB contains the intensity (z-axis) values in a single line
  */
  
  mA = newMatrix(nc, np) ;
  if (mA == NULL) {
    cpl_msg_error(modName, "The function newMatrix has returned NULL");
    return(NULL);
  }  

  mB = newMatrix(1, np) ;
  if (mB == NULL) {
    cpl_msg_error(modName, "The function newMatrix has returned NULL");
    return(NULL);
  }  

 /*
  * Fill up matrices
  */
  for (i = 0; i < np; i++) { /* loop over surface points */
   /* Get x and y value for current point */
    x    = surface[i].x ;
    y    = surface[i].y ;
    
    for (j = 0; j < nc; j++) { /* loop over coefficients */
     /*
      * compute jth coefficient
      */
      mA->data[i+j*np] = ipow(x, degxTab[j]) * ipow(y, degyTab[j]) ; 
    }
    /* mB contains surface values (z-axis) */
    mB->data[i] = surface[i].i ;
  }
  
 /*
  * Solve XA=B by a least-square solution (aka pseudo-inverse).
  */
  mX = lsqMatrix(mA,mB) ;

 /*
  * Delete input matrices
  */
  deleteMatrix(mA) ;
  deleteMatrix(mB) ;

 /*
  * Examine result
  */
  if (mX == NULL) {
    cpl_msg_error(modName, "The function lsqMatrix has returned NULL"); 
    return(NULL);
  }
 /*
  * Store coefficients for output in a single array
  */
  c = (double*) cpl_malloc(nc * sizeof(double)) ;
  /* check if space was allocated */
  if (c == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  for (i = 0; i < nc; i++) {
    c[i] = mX->data[i] ;
  }
  deleteMatrix(mX) ;
  *ncoeffs = nc ;
  
 /*
  * If requested, compute mean squared error
  */
  if (meanSquaredError != NULL) {
    err = 0.00 ;
    for (i = 0; i < np; i++) {
      z = 0.00 ; 
      for (j = 0; j < nc; j++) {
        z += c[j] *
          ipow(surface[i].x, degxTab[j]) *
          ipow(surface[i].y, degyTab[j]);
      }
      /*
        Subtract from the true value, square, accumulate
      */
      err += ipow(surface[i].i-z, 2) ;
    }
    /* Average the error term */
    err /= (double)np ;
    *meanSquaredError = err ;
  }
  cpl_free(degxTab) ;
  cpl_free(degyTab) ;
  return c ;
}


/*------------------------------------------------------------------------
   Function :   fitSurPolErrors
   In       :   list of pixels, # of pixels in the list.
   Out      :   fitted coefficients
                error on each coefficient

   Notice   :   This routine uses fitSurfacePolynomial to fine coefficients
                and then also computes the error on each coefficient
		using a jacknife algorithm
                Added by BG and MS, based on an MS program
---------------------------------------------------------------------------*/
void fitSurPolErrors(VimosPixel *poly,int ndat,double *a_all,double *b_all,double *c_all,double *sa,double *sb,double *sc)
{
   double *ajpar,*bjpar,*cjpar;
   double *ap,*bp,*cp,aasum2,bbsum2,ccsum2;
   double a1,b1,c1,asum,bsum,csum,aasum1,bbsum1,ccsum1;
   double *coeff;
   int i,j,ncoef;
   VimosPixel  *polywork        = NULL;

   ajpar = doubleVector(0,ndat);
   bjpar = doubleVector(0,ndat);
   cjpar = doubleVector(0,ndat);
   ap = doubleVector(0,ndat);
   bp = doubleVector(0,ndat);
   cp = doubleVector(0,ndat);

   polywork      = newPixel(ndat);

   coeff = fitSurfacePolynomial(poly, ndat, NULL, 1, &ncoef, NULL);
   *a_all =  coeff[0];
   *b_all =  coeff[1];
   *c_all =  coeff[2];

	
   for(i = 0; i < ndat; i++)
   {
	for(j = 0; j < i; j++)
        {
		polywork[j].x =  poly[j].x;
        	polywork[j].y =  poly[j].y;
        	polywork[j].i =  poly[j].i;
	}
	for(j = i+1; j < ndat; j++)
        {
		polywork[j-1].x =  poly[j].x;
        	polywork[j-1].y =  poly[j].y;
        	polywork[j-1].i =  poly[j].i;
	}
	coeff = fitSurfacePolynomial(polywork, ndat-1, NULL, 1, &ncoef, NULL);
	ajpar[i] =  coeff[0];
	bjpar[i] =  coeff[1];
	cjpar[i] =  coeff[2];

        ap[i] = ndat * (*a_all) - (ndat-1) * ajpar[i];
        bp[i] = ndat * (*b_all) - (ndat-1) * bjpar[i];
        cp[i] = ndat * (*c_all) - (ndat-1) * cjpar[i];
   }
   asum = 0.;
   bsum = 0.;
   csum = 0.;
   aasum1 = 0.;
   bbsum1 = 0.;
   ccsum1 = 0.;
   for(i = 0; i < ndat; i++)
   {
	asum += ap[i];
	bsum += bp[i];
	csum += cp[i];
	aasum1 += ap[i]*ap[i];
	bbsum1 += bp[i]*bp[i];
	ccsum1 += cp[i]*cp[i];
   }
   aasum2 = (asum * asum) / ndat;
   bbsum2 = (bsum * bsum) / ndat;
   ccsum2 = (csum * csum) / ndat;
   a1 = ((aasum1-aasum2) >= 0.) ? (aasum1-aasum2) : -(aasum1-aasum2);
   b1 = ((bbsum1-bbsum2) >= 0.) ? (bbsum1-bbsum2) : -(bbsum1-bbsum2);
   c1 = ((ccsum1-ccsum2) >= 0.) ? (ccsum1-ccsum2) : -(ccsum1-ccsum2);
/*   a = asum / ndat;
   b = bsum / ndat;
   c = csum / ndat; */
   *sa = sqrt( a1 / (ndat*(ndat-1)) );
   *sb = sqrt( b1 / (ndat*(ndat-1)) );
   *sc = sqrt( c1 / (ndat*(ndat-1)) );

   freeDoubleVector(ajpar,0,ndat);
   freeDoubleVector(bjpar,0,ndat);
   freeDoubleVector(cjpar,0,ndat);
   freeDoubleVector(ap,0,ndat);
   freeDoubleVector(bp,0,ndat);
   freeDoubleVector(cp,0,ndat);
}


/*------------------------------------------------------------------------
   Function :   fit_surface_matrix()
   In       :   list of pixels, # of pixels in the list.
                character string indicating which coefficients should be
                taken into account, maximum polynomial degree.
   Out      :   double * (table of fitted coefficients)
                number of coefficients returned
                mean squared error for the fit.

   Notice   :   This is very similar to the previous function, but with an
                important difference: the definition of "degree of the fit"
		is made according to the size of the matrix used in the fit,
		and not to a polynomial degree. Thus 2nd degree fit means
		that all terms in the series "x^2, x^2 y, X^2 y^2" will
		be used.
---------------------------------------------------------------------------*/

double *fitSurfaceMatrix(VimosPixel *surface, int np, char *controlString,
	int polyDeg, int *ncoeffs, double *meanSquaredError)
{
  int         i, j;
  int         degx, degy;
  VimosMatrix *mA, *mB, *mX;
  double      x, y, z;
  int         nc;
  double      *c;
  double      err;
  int         *degxTab;
  int         *degyTab;
  
  
  /*
   * Fill up look-up table for coefficients to compute
   */
  
  nc = (1+polyDeg)*(1+polyDeg);
  degxTab = cpl_malloc(nc * sizeof(int)) ;
  degyTab = cpl_malloc(nc * sizeof(int)) ;
  
  if (controlString == NULL) {
    i = 0 ;
    for (degy = 0; degy <= polyDeg; degy++) {
      for (degx = 0; degx <= polyDeg; degx++) {
	degxTab[i] = degx ;
	degyTab[i] = degy ;
	i++ ;
      }
    }
  } 
  else {
    nc = buildupPolytabFromString(controlString, polyDeg, degxTab, degyTab) ;
  }
  
  /*
   * Initialize matrices
   * mA contains the polynomial terms in the order described
   * above in each column, for each input point.
   * mB contains the intensity (z-axis) values in a single line
   */
  
  mA = newMatrix(nc, np) ;
  mB = newMatrix(1, np) ;
  
  /*
    Fill up matrices
  */
  for (i = 0; i < np; i++) { /* loop over surface points */
    /* Get x and y value for current point */
    x    = surface[i].x ;
    y    = surface[i].y ;
    
    for (j = 0; j < nc; j++) { /* loop over coefficients */
      /*
        compute jth coefficient
      */
      mA->data[i+j*np] = ipow(x, degxTab[j]) * ipow(y, degyTab[j]) ; 
    }
    /* mB contains surface values (z-axis) */
    mB->data[i] = surface[i].i ;
  }
  
  /*
    Solve XA=B by a least-square solution (aka pseudo-inverse).
  */
  mX = lsqMatrix(mA,mB) ;
  /*
    Delete input matrices
  */
  deleteMatrix(mA) ;
  deleteMatrix(mB) ;
  /*
    Examine result
  */
  if (mX == NULL) {
    return NULL ;
  }
  /*
    Store coefficients for output in a single array
  */
  c = (double*) cpl_malloc(nc * sizeof(double)) ;
  for (i = 0; i < nc; i++) {
    c[i] = mX->data[i] ;
  }
  deleteMatrix(mX) ;
  *ncoeffs = nc ;
  
  /*
    If requested, compute mean squared error
  */
  if (meanSquaredError != NULL) {
    err = 0.00 ;
    for (i = 0; i < np; i++) {
      z = 0.00 ; 
      for (j = 0; j < nc; j++) {
        z += c[j] *
          ipow(surface[i].x, degxTab[j]) *
          ipow(surface[i].y, degyTab[j]);
      }
      /*
        Subtract from the true value, square, accumulate
      */
      err += ipow(surface[i].i-z, 2) ;
    }
    /* Average the error term */
    err /= (double)np ;
    *meanSquaredError = err ;
  }


  cpl_free(degxTab) ;
  cpl_free(degyTab) ;
  return c ;
}


/*---------------------------------------------------------------------------
   Function : buildup_polytab_from_string()
   In       : control string,
              polynomial degree,
              allocated table to fill in for x degrees,
              allocated table to fill in for y degrees,
              -> degx_tab and degy_tab must have allocated at least 
              (1+poly_deg)*(2+poly_deg)/2 integers
   Out      : number of coefficients found
   Job      : translates a control string into a list of polynomial
              degrees for x and y.
   Notice   : returns -1 in case of error.
              A control string is given as:

              "(int,int) (int,int) [...] (int,int)"
              
              each couple (int,int) represents the degree in x and y
              to be computed for the fit. Couples are given in
              parentheses and separated by commas, without any space
              between the parentheses.
              
              Couples are separated from each other by any number of
              blank characters (at least one is required).
              
              The following is a valid control string:
              "(0,0) (1,2) (2,1) (1,1)"
              
              The following are invalid control strings:
              "(0, 0)"		blanks in parentheses
              "( 0 , 0 )"		blanks in parentheses
              "(0,0)(1,2)"	no blank between couples

 ---------------------------------------------------------------------------
*/

int buildupPolytabFromString(char *s, int polyDeg, int * degxTab, int *degyTab)
{
  const char modName[] = "buildupPolytabFromString";
  char *s2;
  char *token;
  int  degx, degy;
  int  nc;
  int  i, j;
  int  ret;

  /*
    Count number of coefficients provided
  */

  if (s == NULL) {
    cpl_msg_error(modName, "Invalid input string");
    pilErrno = 1;
    return(-1);
  }
  
  if (polyDeg < 0) {
    cpl_msg_error(modName, "Invalid input polynomial degree");
    pilErrno = 1;
    return(-1);
  }
  
  if ( (degxTab == NULL) || 
       (degyTab == NULL) ) {
    cpl_msg_error(modName, "Invalid input");
    pilErrno = 1;
    return(-1);
  }
  
  nc = 0 ;
  for (i = 0; i < (int)strlen(s); i++) {
    if (s[i] == ',') {
      nc++ ;
    }
  }

  /*
    Cut the string into tokens, get degrees for x and y
  */
  s2 = cpl_strdup(s) ;
  token = strtok(s2," ") ;
  if (token == NULL) {
    cpl_free(s2) ;
    cpl_msg_error(modName, "No tokens have been found");
    pilErrno = 1;
    return(-1);
  }
  ret = sscanf(token, "(%d,%d)", &degx, &degy) ;
  if (ret != 2) {
    cpl_free(s2) ;
    cpl_msg_error(modName, "Not enough tokens have been found");
    pilErrno = 1;
    return(-1);
  }
  degxTab[0] = degx ;
  degyTab[0] = degy ;

  for (i = 1; i < nc; i++) {
    token = strtok(NULL, " ") ;
    sscanf(token, "(%d,%d)", &degx, &degy) ;
    if (ret != 2) {
      cpl_free(s2) ;
      cpl_msg_error(modName, "Not enough tokens have been found");
      pilErrno = 1;
      return(-1);
    }
    if (degx+degy > polyDeg) {
      cpl_free(s2) ;
      cpl_msg_error(modName, 
        "The sum of degrees of x and y is greater then polynomial degree"); 
      pilErrno = 1;
      return(-1);
    }
    /*
      check for duplicates
    */
    for (j = 0; j < i; j++) {
      if ((degxTab[j] == degx) && 
          (degyTab[j] == degy) ) {
        cpl_free(s2) ;
	cpl_msg_error(modName, "Duplicates have been found");
	pilErrno = 1;
        return(-1);
      }
    }
    degxTab[i] = degx ; 
    degyTab[i] = degy ;
  }
  cpl_free(s2) ;
  return(nc);
}

/*---------------------------------------------------------------------------
   Function	:	fit_slope_robust()
   In 		:	list of dpoints, # of points in the list
   Out 		:	pointer to (newly allocated) 3 doubles
   				y = c[0] + c[1] * x
				c[2] is the median squared error 
   Job		:	fit a slope to a list of points
   Notice	:	very robust - up to 50% outliers in input
 ---------------------------------------------------------------------------*/

double *fitSlopeRobust(VimosDpoint *list, int np)
{
  const char modName[] = "fitSlopeRobust";
  double *c ;
  double *x ;
  double *y ;
  int    i ;
  
  x = cpl_malloc(np * sizeof(double)) ;

  /* check if space was allocated */
  if (x == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  /* check if space was allocated */
  y = cpl_malloc(np * sizeof(double)) ;
  if (y == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  
  for (i = 0; i < np; i++) {
    x[i] = list[i].x ;
    y[i] = list[i].y ;
  }
  
  c = cpl_malloc(3 * sizeof(double)) ;

  /* check if space was allocated */
  if (c == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  if (!robustLinearFit(x, y, np, c, c+1, c+2)) {
    cpl_msg_error(modName, "The function robustLinearFit has returned an error");
    return(NULL);
  }
 
  cpl_free(x) ;
  cpl_free(y) ;
  return(c);
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX_ITERATE             30

static VimosBool robustLinearFit(double *x, double *y, int np, 
                                 double *a, double *b, double *abdev)
{
	int 	i;
	double 	aa, bb, bcomp, b1, b2, del, abdevt,
			f, f1, f2,
			sigb, temp, d, sum ;
	double	sx, sy,
			sxy, sxx,
			chisq ;
	double* arr ;
	double	aa_ls, bb_ls ;
	int		iter ;

	sx = sy = sxx = sxy = 0.00 ;
	for (i=0 ; i<np ; i++) {
		sx  += x[i];
		sy  += y[i];
		sxy += x[i] * y[i];
		sxx += x[i] * x[i];
	}

	del = np * sxx - sx * sx;
	aa_ls = aa  = (sxx * sy - sx * sxy) / del;
	bb_ls = bb  = (np * sxy - sx * sy) / del;

	chisq = 0.00 ;
	for (i=0;i<np;i++) {
		temp = y[i] - (aa+bb*x[i]) ;
		temp *= temp ;
		chisq += temp ;
	}

	arr = cpl_malloc(np * sizeof(double)) ;
	sigb = sqrt(chisq/del);
	b1   = bb ;

	bcomp = b1 ;
	sum = 0.00 ;
	for (i=0 ; i<np ; i++) {
		arr[i] = y[i] - bcomp * x[i];
	}
	aa = medianDouble(arr, np) ;
	abdevt = 0.0;
	for (i=0 ; i<np ; i++) {
		d = y[i] - (bcomp * x[i] + aa);
		abdevt += fabs(d);
		if (y[i] != 0.0) d /= fabs(y[i]);
		if (fabs(d) > 1e-7) sum += (d >= 0.0 ? x[i] : -x[i]);
	}
	f1 = sum ;

	b2   = bb + SIGN(3.0 * sigb, f1);

	bcomp = b2 ;
	sum = 0.00 ;
	for (i=0 ; i<np ; i++) {
		arr[i] = y[i] - bcomp * x[i];
	}
	aa = medianDouble(arr, np) ;
	abdevt = 0.0;
	for (i=0 ; i<np ; i++) {
		d = y[i] - (bcomp * x[i] + aa);
		abdevt += fabs(d);
		if (y[i] != 0.0) d /= fabs(y[i]);
		if (fabs(d) > 1e-7) sum += (d >= 0.0 ? x[i] : -x[i]);
	}
	f2 = sum ;

	if (fabs(b2-b1)<1e-7) {
		*a = aa ;
		*b = bb ;
		*abdev = abdevt / (double)np;
		cpl_free(arr);
		return(VM_TRUE);
	}

	iter = 0 ;
	while (f1*f2 > 0.0) {
		bb = 2.0*b2-b1;
		b1 = b2;
		f1 = f2;
		b2 = bb;

		bcomp = b2 ;
		sum = 0.00 ;
		for (i=0 ; i<np ; i++) {
			arr[i] = y[i] - bcomp * x[i];
		}
		aa = medianDouble(arr, np) ;
		abdevt = 0.0;
		for (i=0 ; i<np ; i++) {
			d = y[i] - (bcomp * x[i] + aa);
			abdevt += fabs(d);
			if (y[i] != 0.0) d /= fabs(y[i]);
			if (fabs(d) > 1e-7) sum += (d >= 0.0 ? x[i] : -x[i]);
		}
		f2 = sum ;
		iter++;
		if (iter>=MAX_ITERATE) break ;
	}
	if (iter>=MAX_ITERATE) {
		*a = aa_ls ;
		*b = bb_ls ;
		*abdev = -1.0 ;
		cpl_free(arr);
		return(VM_FALSE);
	}

	sigb = 0.01 * sigb;
	while (fabs(b2-b1) > sigb) {
		bb = 0.5 * (b1 + b2) ;
		if ((fabs(bb-b1)<1e-7) || (fabs(bb-b2)<1e-7)) break;

		bcomp = bb ;
		sum = 0.00 ;
		for (i=0 ; i<np ; i++) {
			arr[i] = y[i] - bcomp * x[i];
		}
		aa = medianDouble(arr, np) ;
		abdevt = 0.0;
		for (i=0 ; i<np ; i++) {
			d = y[i] - (bcomp * x[i] + aa);
			abdevt += fabs(d);
			if (y[i] != 0.0) d /= fabs(y[i]);
			if (fabs(d) > 1e-7) sum += (d >= 0.0 ? x[i] : -x[i]);
		}
		f = sum ;

		if (f*f1 >= 0.0) {
			f1=f;
			b1=bb;
		} else {
			f2=f;
			b2=bb;
		}
	}
	cpl_free(arr) ;
	*a=aa;
	*b=bb;
	*abdev=abdevt/np;
        return(VM_TRUE);
}
#undef MAX_ITERATE
#undef SIGN


VimosBool stupidLinearFit(double *x, double *y, int np, 
                                 double *a, double *b, double *adev, double *bdev)
{
  int 	 i;
  double aa, bb, b1,  del,   sigb, siga, temp ;
  double sx, sy, sxy, sxx, chisq ;

  sx = sy = sxx = sxy = 0.00 ;
  for (i = 0; i < np; i++) {
    sx  += x[i];
    sy  += y[i];
    sxy += x[i] * y[i];
    sxx += x[i] * x[i];
  }

  del = np * sxx - sx * sx;
  aa  = (sxx *  sy - sx * sxy) / del;
  bb  = (np * sxy - sx * sy) / del;
  
  chisq = 0.00 ;
  for (i = 0; i < np;i++) {
    temp = y[i] - (aa+bb*x[i]) ;
    temp *= temp ;
    chisq += temp ;
  }

  sigb = sqrt((chisq/del)*(np/(np-2)));
  siga = sqrt ((chisq/del)*sxx/(np-2));
  b1   = bb ;
  
  *a = aa;
  *b = bb;
  *bdev = sigb;
  *adev = siga;
  return VM_TRUE;
}
/**
 * @memo
 *   Generating a control string for the function fitSurfacePolynomial,
 *   listing all the coefficients of an incomplete polynomial here named
 *   conventionally "VIMOS Model". The VIMOS Model is defined by giving
 *   the max allowed power for, respectively, the x and the y variables 
 *   in the polynomial. Then, all terms containing a power of x and y 
 *   less than or equal to their max allowed power, are inserted in 
 *   the model. In this way the polynomial has (xord+1)*(yord+1) free 
 *   coefficients, and has a degree xord+yord (differently from the 
 *   definition of a complete polinomial of order n, having (n+1)*(n+2)/2 
 *   coefficents).
 *
 * @return Pointer to the generated string.
 *
 * @param xord    max degree of variable x
 * @param yord    max degree of variable y
 *
 * @doc
 *   Enough memory is allocated for the result control string,
 *   then the control string is generated. Memory should be 
 *   deallocated afterwards.
 *
 * @author C. Izzo
 */

char *createVimosCtrlStr(int xord, int yord) 
{
  char      *controlString;
  char      *pos;
  int        nCoeff = (xord+1)*(yord+1);
  int        powx = 0;
  int        powy = 0;
  int        valx = 1;
  int        valy = 1;
  int        spacex, spacey;
  int        i, j;

  if (xord < 0 || yord < 0) return NULL;

 /*
  *  What follows is just to determine how much memory to allocate
  *  for the resulting control string
  */

 /* 
  *  Highest power of 10 less then xord and yord: exponent, and
  *  figure itself.
  */
  if (xord > 0) powx = log10((double) xord);
  if (yord > 0) powy = log10((double) yord);

  for (i=0; i<powx; i++) valx *= 10;
  for (i=0; i<powy; i++) valy *= 10;

 /*
  *  Compute number of digits necessary for the x and the y indexes
  */
  spacex = (xord+1)*(powx+1);
  while (powx) {
    powx--;
    spacex -= valx;
    valx /=10;
  }
  spacex *= yord+1;

  spacey = (yord+1)*(powy+1);
  while (powy) {
    powy--;
    spacey -= valy;
    valy /=10;
  }
  spacey *= xord+1;

 /*
  *  Generate control string
  */
  pos = controlString = (char *) cpl_malloc((4*nCoeff+spacex+spacey)*
                                            sizeof(char));

  for (i=0; i<=xord; i++) {
    for (j=0; j<=yord; j++) {
      if (i || j) sprintf(pos, " (%d,%d)", i, j);
      else        sprintf(pos, "(%d,%d)", i, j);
      pos += strlen(pos);
    }
  }

  return(controlString);
}



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  These are functions for fitting & c. (NR- or IDL- like)

------------------------------------------------------------------------------
*/

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void gaussFunc(float anX, float coeffs[], float *yFit, float derY[], 
               int nTerms)

  Description:
  Computes the sum of a gaussian + variable continuum 
  as described by the function:

     F(x)= A0 * EXP(-z^2 / 2) + A3 + A4*x + A5*x^2
     with:
     z = (x - A1) / A2

  where:
  A0 = height of exponential
  A1 = center of exonential
  A2 = sigma (the width)
  A3 = constant term
  A4 = linear term
  A5 = quadratic term


  Input:
  float anX
  x value

  float coeffs[]
  Array of A* coefficients

  float *yFit
  Computed f(x)

  float derY[]
  Partial derivatives with respect to the A*

  int nTerms
  Number of A coefficients one wants to fit. Default to 6.

  Updates:
  15 Jun 00: Created (AZ)

------------------------------------------------------------------------------
*/

void gaussFunc(float anX, float coeffs[], float *yFit, float derY[], 
               int nTerms)
{
  float z, ez;

  if (coeffs[3] != 0.0)
    {
     /* get z */
     z = (anX - coeffs[2]) / coeffs[3];
     /* gaussian part */
     ez = exp( -(pow(z,2.)) / 2.);
    }
  else
    {
     z = 100.;
     ez = 0.0;
    }
  switch (nTerms)
    {
     case 3:
      {
       *yFit = coeffs[1] * ez;
       break;
      }
     case 4:
      {
       *yFit = coeffs[1]*ez + coeffs[4];
       break;
      }
     case 5:
      {
       *yFit = coeffs[1]*ez + coeffs[4] + coeffs[5]*anX;
       break;
      }
     case 6:
      {
       *yFit = coeffs[1]*ez + coeffs[4] + coeffs[5]*anX + 
              coeffs[6]*(pow(anX,2.));
       break;
      }
     default:
      {
      }
    }
 
  derY[1] = ez;

  if (coeffs[3] != 0.)
    {
     derY[2] = coeffs[1] * ez * (z/coeffs[3]);
    }
  else
    {
      derY[2] = 0.;
    }
  derY[3] = derY[2] * z;
  if (nTerms > 3) derY[4] =  1.;
  if (nTerms > 4) derY[5] = anX;
  if (nTerms > 5) derY[6] = pow(anX, 2.);
}



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void gaussJordan(float **a, int n, float **b, int m)

  Description:
  Gauss-Jordan elimination (gaussj from NR). Used by LevenMar
 
  Input:
  float **a
  Input matrix

  int n
  Dimension of the input matrix [1...n][1...n]

  float **b
  Input m right-hand sided vectors

  int m
  Dimension of the m right-hand sided vectors [1...n][1...m]

  Updates:
  16 Jun 00: Created (AZ)

------------------------------------------------------------------------------
*/

#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussJordan(float **a, int n, float **b, int m)
{
  int *indexc, *indexr, *iPivot;
  int i, j, k, l, ll;
  int iCol = 0;
  int iRow = 0;
  float big, dum, pivotInv, temp;

  indexc = intVector(1, n);
  indexr = intVector(1, n);
  iPivot = intVector(1, n);
  for (j=1; j<=n; j++) iPivot[j] = 0;
  for (i=1; i<=n; i++)
     {
      big = 0.0;
      for (j=1; j<=n; j++)
         if (iPivot[j] != 1)
           for (k=1; k<=n; k++)
              {
	       if (iPivot[k] == 0)
                 {
	          if (fabs(a[j][k]) >= big)
                    {
                     big = fabs(a[j][k]);
        	     iRow = j;
		     iCol = k;
		    }
	         } 
               else if (iPivot[k] > 1)
                 {
                  puts("gaussJordan: Singular Matrix-1");
                  abort();
                 }
	      }
      ++(iPivot[iCol]);
      if (iRow != iCol)
        {
	 for (l=1; l<=n; l++) SWAP(a[iRow][l],a[iCol][l])
         for (l=1; l<=m; l++) SWAP(b[iRow][l],b[iCol][l])
        }
      indexr[i] = iRow;
      indexc[i] = iCol;
      if (a[iCol][iCol] == 0.0)
        {
         puts("gaussJordan: Singular Matrix-2");
         abort();
        }
      pivotInv = 1.0 / a[iCol][iCol];
      a[iCol][iCol] = 1.0;
      for (l=1; l<=n; l++) a[iCol][l] *= pivotInv;
      for (l=1; l<=m; l++) b[iCol][l] *= pivotInv;
      for (ll=1; ll<=n; ll++)
         if (ll != iCol)
           {
            dum = a[ll][iCol];
            a[ll][iCol] = 0.0;
            for (l=1; l<=n; l++) a[ll][l] -= a[iCol][l]*dum;
            for (l=1; l<=m; l++) b[ll][l] -= b[iCol][l]*dum;
           }
     }
  for (l=n; l>=1; l--)
     {
      if (indexr[l] != indexc[l])
	for (k=1; k<=n; k++) SWAP(a[k][indexr[l]],a[k][indexc[l]]);
     }
  freeIntVector(iPivot, 1, n);
  freeIntVector(indexr, 1, n);
  freeIntVector(indexc, 1, n);
}
#undef SWAP


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void minimizeChisq(float x[], float y[], float sig[], int n, float coeffs[],
                   int nToFit[], int nTerms, float **alpha, float beta[], 
                   float *chisq, 
                   void (*gaussFunc)(float, float [], float *, float [], int))

  Description:
  Evaluates the linearized fitting matrix and vector beta and calculates Chi 
  square (mrqcof from NR). Used by levenMar
 
  Input:
  float x[]
  Input x data points

  float y[]
  Input y data points

  float sig[]
  Individual standard deviations on data points

  int n
  Dimension of x and y arrays (number of data points)

  float coeffs[]
  Coefficients of the non-linear function, to be best-fitted

  int nToFit[]
  Array to select which coefficients are to be fitted (nToFit[i]=1 fit, 
  nToFit[i]=0 don't fit)

  int nTerms
  Number of coefficients of the non-linear function

  float **alpha
  Fitting matrix

  float beta[]
  Beta vector

  float *chisq
  Chi square

  void (*gaussFunc)(float, float [], float *, float [], int)

  Updates:
  16 Jun 00: Created (AZ)

------------------------------------------------------------------------------
*/

void minimizeChisq(float x[], float y[], float sig[], int n, float coeffs[], 
            int nToFit[], int nTerms, float **alpha, float beta[], 
            float *chisq, 
            void (*gaussFunc)(float, float [], float *, float [], int))
{
  int i, j, k, l, m, nFit=0;
  float ymod, wt, sig2i, dy, *dyda;

  dyda = floatVector(1, nTerms);
  for (j=1; j<=nTerms; j++) if (nToFit[j]) nFit++;

  for (j=1; j<=nFit; j++)
     {
      for (k=1; k<=j; k++) alpha[j][k] = 0.0;
      beta[j] = 0.0;
     }
  *chisq = 0.0;

  for (i=1; i<=n; i++)
     {
      (*gaussFunc)(x[i], coeffs, &ymod, dyda, nTerms);
      sig2i = 1.0 / (sig[i]*sig[i]);
      dy = y[i] - ymod;
      for (j=0, l=1; l<=nTerms; l++)
         {
	  if (nToFit[l])
            {
	     wt = dyda[l] * sig2i;
	     for (j++, k=0, m=1; m<=l; m++)
		if (nToFit[m]) alpha[j][++k] += wt*dyda[m];
	     beta[j] += dy*wt;
	    }
	 }
      *chisq += dy*dy*sig2i;
     }
  for (j=2; j<=nFit; j++) for (k=1; k<j; k++) alpha[k][j] = alpha[j][k];

  freeFloatVector(dyda, 1, nTerms);
}


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void expandCovar(float **covar, int nTerms, int nToFit[], int nFit)

  Description:
  Expand the covariance matrix (covsrt from NR). Used by levenMar

  Input:
  float **covar
  Covariance matrix

  int nTerms
  Number of coefficients of the function to be fitted

  int nToFit[]
  Array to select which coefficients are to be fitted (nToFit[i]=1 fit, 
  nToFit[i]=0 don't fit)

  int nFit
  Number of nToFit terms different from zero (i.e. how many to be fitted)

  Updates:
  16 Jun 00: Created (AZ)

------------------------------------------------------------------------------
*/

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void expandCovar(float **covar, int nTerms, int nToFit[], int nFit)
{
  int i, j, k;
  float swap;

  for (i=nFit+1; i<=nTerms; i++)
      for (j=1; j<=i; j++) covar[i][j] = covar[j][i] = 0.0;
  k = nFit;
  for (j=nTerms; j>=1; j--)
     {
      if (nToFit[j])
        {
	 for (i=1; i<=nTerms; i++) SWAP(covar[i][k], covar[i][j])
	 for (i=1; i<=nTerms; i++) SWAP(covar[k][i], covar[j][i])
	 k--;
	}
     }
}
#undef SWAP



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void levenMar(float x[], float y[], float sig[], int n, float coeffs[], 
              int nToFit[], int nTerms, float **covar, float **alpha,
              float *chisq, float *aLambda)

  Description:
  Levenberg-Marquardt method. (modified mrqmin from NR)
  It attempts to reduce the chi-square of a fit between a set of data points 
  (x[i],y[i]) and a non-linear function that depends on NTERMS coefficients. 
  The non-linear function is provided by the function gaussFunc. 
  Recursively called until convergence is achieved.
 
  Input:
  float x[]
  Input x data points

  float y[]
  Input y data points

  float sig[]
  Individual standard deviations on data points

  int n
  Dimension of x and y arrays (number of data points)

  float coeffs[]
  Coefficients of the non-linear function, to be best-fitted

  int nToFit[]
  Array to select which coefficients are to be fitted (nToFit[i]=1 fit, 
  nToFit[i]=0 don't fit)

  int nTerms
  Number of coefficients of the non-linear function

  float **covar
  Covariance matrix

  float **alpha
  Fitting matrix

  float *chisq
  Chi square

  float *aLambda
  Convergence parameter

  Updates:
  16 Jun 00: Created (AZ)

------------------------------------------------------------------------------
*/

void levenMar(float x[], float y[], float sig[], int n, float coeffs[], 
            int nToFit[], int nTerms, float **covar, float **alpha, 
            float *chisq, float *aLambda)
{
  int j,k,l,m;
  static int nFit;
  static float outChisq,*tryCoeffs,*beta,*da,**oneda;

  if (*aLambda < 0.0)
    {
     tryCoeffs = floatVector(1,nTerms);
     beta = floatVector(1,nTerms);
     da = floatVector(1,nTerms);
     for (nFit=0, j=1; j<=nTerms; j++) if (nToFit[j]) nFit++;
     oneda = Matrix(1,nFit,1,1);
     *aLambda = 0.001;
     minimizeChisq(x, y, sig, n, coeffs, nToFit, nTerms, alpha, beta, chisq, 
            gaussFunc);
     outChisq = (*chisq);
     for (j=1; j<=nTerms; j++) tryCoeffs[j]=coeffs[j];
    }
  for (j=0,l=1;l<=nTerms;l++)
     {
      if (nToFit[l])
        {
	 for (j++, k=0, m=1; m<=nTerms; m++)
            {
	     if (nToFit[m])
               {
		k++;
		covar[j][k] = alpha[j][k];
	       }
	    }
	 covar[j][j] = alpha[j][j] * (1.0 + (*aLambda));
	 oneda[j][1] = beta[j];
	}
     }
  gaussJordan(covar, nFit, oneda, 1);
  for (j=1; j<=nFit; j++) da[j]=oneda[j][1];
  if (*aLambda == 0.0)
    {
     expandCovar(covar, nTerms, nToFit, nFit);
     freeMatrix(oneda, 1, nFit, 1, 1);
     freeFloatVector(da, 1, nTerms);
     freeFloatVector(beta, 1, nTerms);
     freeFloatVector(tryCoeffs, 1, nTerms);
     return;
    }
  for (j=0, l=1; l<=nTerms; l++) 
     if (nToFit[l]) tryCoeffs[l] = coeffs[l]+da[++j];
  minimizeChisq(x, y, sig, n, tryCoeffs, nToFit, nTerms, covar, da, chisq, 
              gaussFunc);
  if (*chisq < outChisq)
    {
     *aLambda *= 0.1;
     outChisq = (*chisq);
     for (j=0, l=1; l<=nTerms; l++)
        {
	 if (nToFit[l])
           {
	    for (j++, k=0, m=1; m<=nTerms; m++)
               {
		if (nToFit[m])
                  {
		   k++;
		   alpha[j][k] = covar[j][k];
		  }
	       }
	    beta[j] = da[j];
	    coeffs[l] = tryCoeffs[l];
           }
	}
    }
  else
    {

      /*ALEX:*/
      /*      printf("CHI: %12.9f DIFF: %12.9f\n",(*chisq), 
	      (outChisq -(*chisq))/outChisq);
      */

     *aLambda *= 10.0;
     *chisq = outChisq;

    }
}
#undef NRANSI


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void fit1DGauss(VimosFloatArray *xArray, VimosFloatArray *yArray, 
                float coeffArray[], int nTerms)

  Description:
  Non-linear least squares fit with a function made of gaussian+variable
  continuum. Uses Levenberg-Marquard minimization method for chi square.
  See the function gaussFunc for details on the fitted function.

  Input:
  VimosFloatArray *xArray
  Array of independent variables

  VimosFloatArray *yArray
  Array containing the dependent variable

  float coeffArray[]
  Array of NTERMS coefficients A, with initial estimates

  int nTerms
  number of coefficients A to be fitted

  Updates:
  15 Jun 00: Created (AZ)
  23 Jun 00: Working (AZ)
  28 Feb 01: Ongoing tests on fit (AZ)

------------------------------------------------------------------------------
*/
/***************************************************************/
/* Fit the equation y=f(x) where:                              */
/*                                                             */
/*     F(x)= A0 * EXP(-z^2 / 2) + A3 + A4*x + A5*x^2           */
/*     and:                                                    */
/*     z = (x - A1) / A2                                       */
/*                                                             */
/*  A0 = height of exponential                                 */
/*  A1 = center of exonential                                  */
/*  A2 = sigma (the width)                                     */
/*  A3 = constant term                                         */
/*  A4 = linear term                                           */
/*  A5 = quadratic term                                        */
/*                                                             */
/*  See gaussFunc function                                     */
/***************************************************************/

#define ITER_MAX 60
#define TOLERANCE 0.001

void fit1DGauss(VimosFloatArray *xArray, VimosFloatArray *yArray, 
                float coeffArray[], int nTerms)
{
  int i, l, iMin, iMax, i0, nPoints;
  int *iCoeffs;
  float Chisq, aLambda;
  float *sig;
  float chiDiff, prevChisq, prevALambda;
  float **alpha;
  float **covar;
  float *x, *y;
  double dy, factor;
  double yMin, yMax;
  double xMin, xMax;
  double *coeffs1DPoly;

  VimosFloatArray *yFit;
  VimosFloatArray *yDiff;
  VimosDpoint *toPolyFit;

  /* get number of points */
  nPoints = xArray->len;

  yFit = newFloatArray(nPoints);
  yDiff = newFloatArray(nPoints);
  sig = floatVector(1,nPoints);
  iCoeffs = intVector(1,nTerms);

  /*  puts("CONTROLLA IL CASTING DA FLOAT A DOUBLE");
  puts("                                       ");
  puts("WEIGHTS! SIGMAS NOW SET TO POISSONIAN NOISE");
  */

  if ((nTerms < 3) || (nTerms > 6))
    {
     puts(" fit1DGauss: nTerms at least 3, at most 6");
     abort();
    }

  /* to work with curvefit it is better to have arrays. I transform to    */
  /* VimosDPoint so to use TAO's fit1DPoly function (already in vimos lib */

  toPolyFit = newDpoint(nPoints);

  /* fill toPolyFit after casting */
  for (l=0; l<nPoints; l++)
    {
      toPolyFit[l].x = (double)xArray->data[l];
      toPolyFit[l].y = (double)yArray->data[l];
    }

  /* fit a straight line to get A3, A4 */
  coeffs1DPoly = fit1DPoly(1, toPolyFit, nPoints, NULL);

  /*  printf("COEFF_POLYFIT(0): %10.7f COEFF_POLYFIT(1): %10.7f\n" , 
      coeffs1DPoly[0],coeffs1DPoly[1]);*/

  /* evaluate y from fit and difference with data points */
  for (i=0; i<nPoints; i++)
    {
     yFit->data[i] = (float)coeffs1DPoly[0] +
                      (float)coeffs1DPoly[1]*(xArray->data[i]);
     yDiff->data[i] = yArray->data[i] - yFit->data[i];
    }

  /*  for (i=0; i<nPoints; i++) printf("YFIT: %10.7f\n" , yFit->data[i]);*/

  /*  Get minimum and maximum of difference between linear fit and data */
  /* and save the correspondent index */
  /*  We do not treat the edges, so we start at pixel 1 and end at n-1 */
  /* CHECK THIS! */
  yMin = yMax = (double)yDiff->data[0];
  xMin = xMax = (double)xArray->data[0];
  iMin = iMax = 0;

  for (i = 1; i < (nPoints-1); i++)
     {
      yMin = MIN(yMin, (double)yDiff->data[i]);
      if (yMin == (double)yDiff->data[i]) iMin = i;

      yMax = MAX(yMax, (double)yDiff->data[i]);
      if (yMax == (double)yDiff->data[i]) iMax = i;
     }
  xMax = (double)xArray->data[iMin];
  xMin = (double)xArray->data[iMax];

  /* this is to check if it is an absorption or emission line */
  if (ABS(yMax) > ABS(yMin))
    {
     i0 = iMax;
    }
  else
    {
     i0=iMin;
    }
 
  /* never take edges: */
  if (i0 < 1) i0 = 1;
  if (i0 > (nPoints-2)) i0 = (nPoints-2);

  /* this is the difference between extreme and mean */
  dy = (double)yDiff->data[i0];

  /* now look for the coefficient A0=width. The width is found by searching */
  /* out from the extrema until a point is found less than the 1/e value */
  /* save also the index at which this happens */
  factor = dy / exp(1.);
  i = 0;

  while( ((i0+i+1) < nPoints) && ((i0-i) > 0) && 
         (ABS((double)yDiff->data[i0+i]) > ABS(factor)) &&
         (ABS((double)yDiff->data[i0-i]) > ABS(factor)) ) i++;

  /* now set the estimates for the coefficient array, to be passed */
  /* to levenMar */
  coeffArray[1] = yDiff->data[i0];
  coeffArray[2] = xArray->data[i0];
  coeffArray[3] = ABS(xArray->data[i0] - xArray->data[i0+i]);
  if (nTerms >= 4) coeffArray[4] = (float)coeffs1DPoly[0];
  if (nTerms >= 5) coeffArray[5] = (float)coeffs1DPoly[1];
  if (nTerms == 6) coeffArray[6] = 0.;

  /* always fit all the parameters */

  for (l=1; l<=nTerms; l++) iCoeffs[l]=1;

  /* create matrix alpha and delta for levenMar */
  alpha = Matrix(1,nTerms,1,nTerms);
  covar = Matrix(1,nTerms,1,nTerms); 

  /* create and fill floatVectors for passing to levenMar */
  x = floatVector(1,nPoints);
  y = floatVector(1,nPoints);
  for (i=0; i<nPoints; i++)
    {
     x[i+1]=xArray->data[i];
     y[i+1]=yArray->data[i];
    }


  /* initialize sigmas: we do not know the individual standard deviations */

  /*ALEX!*/ /* SET THEM TO 1 SINCE THEY ARE USED TO COMPUTE WEIGHTS */
  /*ALEX!  for (l=1; l<=nPoints; l++) sig[l] = 1.0; */
  /*ALEX! SET THEM TO POISSONIAN NOISE */
  for (l=1; l<=nPoints; l++) sig[l] = pow(y[l],0.5);
  /*for (l=1; l<=nPoints; l++)printf("Y: %15.7f SIG: %15.7f\n",y[l],sig[l]);*/
  /*  printf("NPOINTS: %4d\n",nPoints);*/


  /*  for (i=1;i<=nTerms;i++) printf("A BEFORE:�%10.7f\n",coeffArray[i]);*/

  /* initialize aLambda */
  aLambda = -1.0;

  /* first call to levenMar with aLambda=-1.0, to initialize parameters */
  levenMar(x, y, sig, nPoints, coeffArray, iCoeffs, nTerms, covar, 
         alpha, &Chisq, &aLambda);

  for (l=0; l<ITER_MAX; l++)
    {
     prevChisq = Chisq;
     prevALambda = aLambda;

     levenMar(x, y, sig, nPoints, coeffArray, iCoeffs, nTerms, covar, 
              alpha, &Chisq, &aLambda);

     /* chiDiff is the fractional amount of variation in ChuSquare value */
     /* between two subsequent iterations. It is checked against TOLERANCE */

     chiDiff = (prevChisq - Chisq) / prevChisq;

     /*printf("IT: %4d C:�%10.7f PC:�%10.7f D: %12.9f L:�%12.9f PL: %12.9f \n",
       l, Chisq,prevChisq,chiDiff,aLambda,prevALambda); */

     /* if Chisq is decreasing and changes very little from previous iter */
     if ((aLambda < prevALambda) && (chiDiff <= TOLERANCE)) 
       {
	 /*for (li=1;li<=nTerms; li++) 
	   printf("A END: %16.12f\n",coeffArray[li]);*/

        /* last call to levenMar with aLambda = 0 to solve alpha and covar */
        aLambda = 0.0;
        levenMar(x, y, sig, nPoints, coeffArray, iCoeffs, nTerms, covar, 
                 alpha, &Chisq, &aLambda);
         
	/*        printf("Convergence achieved at iteration: %4d\n", l+1);*/
        return;
       }

    }


}


#undef MIN
#undef MAX
#undef ABS
#undef ITER_MAX
#undef TOLERANCE



/*
  Given the coefficients A from fit1DGauss, given a value for theX,
  it computes f(x) by using the same formula as in gaussFunc (i.e. 
  gaussian + variable continuum).
*/

float evalYFit(float theX, float coeff[])
{
  float z, ez;
  float Yfit;
  if (coeff[3] != 0.0)
    {
     /* get z */
     z = (theX - coeff[2]) / coeff[3];
     /* gaussian part */
     ez = exp( -(pow(z,2.)) / 2.);
    }
  else
    {
     z = 100.;
     ez = 0.0;
    }
  Yfit = coeff[1]*ez + coeff[4] + 
         coeff[5]*theX + coeff[6]*(pow(theX,2.));
          
  return(Yfit);

}    



/*
  Computes line flux. First it makes a fit with fit1DGauss (gaussian
  + variable continuum). Then it integrates the fitted function twice, one
  for continuum+gaussian (using all coefficients from fit1DGauss) and one for
  continuum alone (use only some of the coefficients). the difference between
  integrated_continuum+gaussian and integrated_continuum gives the flux in 
  that line.
*/

float evalLineFlux(VimosFloatArray *anX, VimosFloatArray *anY, 
                   float coeffs[], int nT)
{
  int i, lun;

  float integratedTotal, integratedContinuum;
  float start, end, lineF;
  float *selectedCoeffs;

  selectedCoeffs = floatVector(1, nT);

  for (i=1; i<=nT; i++) coeffs[i]=0.;

  /* perform fit gaussian+"variable continuum" */
  fit1DGauss(anX, anY, coeffs, nT);

  /*  for (l=1; l<=nT; l++) printf(" %20.7f",coeffs[l]);*/

  /* start to evaluate integrals for the gaussian+"variable continuum" */
  /* and "variable continuum" alone */
  integratedTotal = 0.0;
  integratedContinuum = 0.0;

  lun = anX->len;
  start=anX->data[0];
  end = anX->data[lun-1];

  /* define first which coeffs are to be considered: we want */
  /* gaussian+"variable continuum" alone, so selectedCoeffs[*] = coeffs[*] */

  for (i=1; i<=nT; i++) selectedCoeffs[i]=coeffs[i];

  /* integrate the whole profile */
  integratedTotal = rombergInt(evalYFit, selectedCoeffs, start, end);

  /*  printf("TOTAL: %20.10f",integratedTotal);*/

  /* define first which coeffs are to be considered: we want */
  /* "variable continuum" alone, so selectedCoeffs[1,2,3] = 0. */
  /* and selectedCoeffs[4,5,6] = coeffs[4,5,6] */

  for (i=1; i<=3; i++) selectedCoeffs[i]=0.;
  for (i=4; i<=nT; i++) selectedCoeffs[i]=coeffs[i];

  /* integrate the continuum part alone */
  integratedContinuum = rombergInt(evalYFit, selectedCoeffs, start, end);
 
  /*  printf("CONTINUUM: %20.10f",integratedContinuum);*/

  /* line flux */
  lineF = integratedTotal - integratedContinuum;

  return (lineF);

}

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  End of functions for fitting & c. (NR- or IDL- like)

------------------------------------------------------------------------------
*/
