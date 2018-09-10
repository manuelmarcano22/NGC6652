/* $Id: vmfit.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_FIT_H
#define VM_FIT_H

#include <pilmacros.h>

#include <vmtypes.h>
#include <vmmatrix.h>


PIL_BEGIN_DECLS

/*---------------------------------------------------------------------------
   Function :   fit_1d_poly()
   In       :   requested polynomial degree
                a list of pixel positions + number of pixels in the list
                (out) mean squared error, set to NULL if you do not want
                to compute it.
   Out      :   newly allocated array containing fit coefficients
   Job      :   fit a polynomial to a list of pixel positions
   Notice   :
                The fitted polynomial is such that:
                y = c[0] + c[1].x + c[2].x^2 + ... + c[n].x^n
                So requesting a polynomial of degree n will return n+1
                coefficients. Beware that with such polynomials, two
                input points shall never be on the same vertical!
 ---------------------------------------------------------------------------*/

double *fit1DPoly(int polyDeg, VimosDpoint *list, int np, 
                  double *meanSquaredError);



/*---------------------------------------------------------------------------
   Function :   fit_surface_polynomial()
   In       :   list of pixels, # of pixels in the list.
                character string indicating which coefficients should be
                taken into account, maximum polynomial degree.
   Out      :   double * (table of fitted coefficients)
                number of coefficients returned
                mean squared error for the fit.
   Job      :   fit a 2d surface with a polynomial in (x,y).
   Notice   :   To define which coefficients should be computed, either
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
                (0,0) is degx=0 and degy=0 -> constant term c[0]
                (1,0) is degx=1 and degy=0 -> term in x     c[1]
                (2,0) is degx=2 and degy=0 -> term in x^2   c[2]
                (1,1) is degx=1 and degy=1 -> term in x.y   c[3]

                The maximal polynomial degree indicates the highest sum
                for X and Y degrees. Example: for poly_deg=3, only the
                following terms can be computed:

                1       x       x^2     x^3
                y       x.y     x^2.y
                y^2     x.y^2
                y^3

                If you do not provide any control string, use NULL as
                argument and set the polynomial degree to what you wish,
                all terms satisfying the condition (degx+degy<=polydeg)
                will be taken into account for the fit.
                
 ---------------------------------------------------------------------------*/

double *fitSurfacePolynomial(VimosPixel *surface, int np,
                             char *controlString, int polyDeg, int *ncoeffs,
                             double *meanSquaredError);


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
void fitSurPolErrors(VimosPixel *poly,int ndat,double *a_all,double *b_all,double *c_all,double *sa,double *sb,double *sc);


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

double *fitSurfaceMatrix(VimosPixel *surface, int np,
                             char *controlString, int polyDeg, int *ncoeffs,
                             double *meanSquaredError);

/*---------------------------------------------------------------------------
   Function :   buildup_polytab_from_string()
   In       :   control string,
                polynomial degree,
                allocated table to fill in for x degrees,
                allocated table to fill in for y degrees,
                -> degx_tab and degy_tab must have allocated at least 
                   (1+poly_deg)*(2+poly_deg)/2 integers
   Out      :   number of coefficients found
   Job      :   translates a control string into a list of polynomial
                degrees for x and y.
   Notice   :   returns -1 in case of error.
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
                "(0, 0)"        blanks in parentheses
                "( 0 , 0 )"     blanks in parentheses
                "(0,0)(1,2)"    no blank between couples

 ---------------------------------------------------------------------------*/

int buildupPolytabFromString(char *s, int polyDeg, int *degxTab, int *degyTab);



/*---------------------------------------------------------------------------
   Function :   fit_slope_robust()
   In       :   list of dpoints, # of points in the list
   Out      :   pointer to (newly allocated) 3 doubles
                y = c[0] + c[1] * x
                c[2] is the median squared error 
   Job      :   fit a slope to a list of points
   Notice   :   very robust - up to 50% outliers in input
 ---------------------------------------------------------------------------*/



double *fitSlopeRobust(VimosDpoint *list, int np);

char *createVimosCtrlStr(int xord, int yord);

VimosBool stupidLinearFit(double *x, double *y, int np, 
                                 double *a, double *b, double *adev, double *bdev);

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
               int nTerms);



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

void gaussJordan(float **a, int n, float **b, int m);


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
                  void (*gaussFunc)(float, float [], float *, float [], int));

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

void expandCovar(float **covar, int nTerms, int nToFit[], int nFit);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void levenMar(float x[], float y[], float sig[], int n, float coeffs[], 
              int nToFit[], int nTerms, float **covar, float **alpha,
              float *chisq, float *aLambda);

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
              float *chisq, float *aLambda);




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

------------------------------------------------------------------------------
*/

void fit1DGauss(VimosFloatArray *xArray, VimosFloatArray *yArray, 
                float coeffArray[], int nTerms);




/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  float evalYFit(float theX, float coeff[])

  Description:
  Given the coefficients A[] from fit1DGauss, given a value for theX,
  it computes f(x) by using the same formula as in gaussFunc (i.e. 
  gaussian + variable continuum).
  Needed by Romberg integration (rombergInt function).

  Input:
  float theX
  X value where to compute the y according to the same formula as in gaussFunc

  float coeff[]
  Coefficients defining the function to be evaluated

  Updates:
  03 Jul 00: Created (AZ)

------------------------------------------------------------------------------
*/

float evalYFit(float theX, float coeff[]);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  float evalLineFlux(VimosFloatArray *anX, VimosFloatArray *anY, 
                     float coeffs[], int nT)

  Description:
  Computes line flux. First it makes a fit with fit1DGauss (gaussian
  + variable continuum). Then it integrates the fitted function twice, one
  for continuum+gaussian (using all coefficients from fit1DGauss) and one for
  continuum alone (use only some of the coefficients). the difference between
  integrated_continuum+gaussian and integrated_continuum gives the flux in 
  that line.
  (Used to compute sky line intensity in VmIfuCalRel)
 
  Input:
  VimosFloatArray *anX
  Array containing the independent variable X

  VimosFloatArray *anY
  Array containing the dependent variable Y

  float coeffs[]
  Array of coefficients from the fit in fit1DGauss

  int nT
  Size  of coeffs[] = number of terms to be fitted by fit1DGauss

  Updates:
  03 Jul 00: Created (AZ)
  14 Jul 00: Working (AZ)

------------------------------------------------------------------------------
*/

float evalLineFlux(VimosFloatArray *anX, VimosFloatArray *anY, 
                   float coeffs[], int nT);

PIL_END_DECLS

#endif /* VM_FIT_H */
