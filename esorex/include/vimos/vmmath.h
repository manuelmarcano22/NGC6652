/* $Id: vmmath.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_MATH_H
#define VM_MATH_H

#include <limits.h>

#include <pilmacros.h>

#include <vmmatrix.h>


PIL_BEGIN_DECLS

/* Cell size for surface fit */
#define BOXSIZE  100

/* defines for resampling code */
/* Number of tabulations in kernel  */
#define TABSPERPIX      (1000)
#define KERNEL_WIDTH    (2.0)
#define KERNEL_SAMPLES  (1+(int)(TABSPERPIX * KERNEL_WIDTH))

#define TANH_STEEPNESS	(5.0)
#define MIN_DIVISOR     ((double)1e-10)

#define MEANDEV_TO_SIGMA (1.25)

#ifndef PI_NUMB
#define PI_NUMB     (3.1415926535897932384626433832795)
#endif

/* Arithmetic operations */

typedef enum _VIMOS_OPER_TYPE_
{
  VM_OPER_ADD,
  VM_OPER_SUB,
  VM_OPER_MUL,
  VM_OPER_DIV
} VimosOperator;

typedef enum _NORM_METHOD
{
  MEAN,
  MEDIAN,
  MODE
} Method;

/*---------------------------------------------------------------------------
   Function	:	ipow()
   In 		:	double, int
   Out 		:	double
   Job		:	same as pow(x,y) but for integer values of y
   Notice	:	much faster than the math function due to the integer.
   				On some compilers, this optimization is already done in
				the pow() function, we do not rely on that fact.
 ---------------------------------------------------------------------------*/


double ipow(double x, int i);



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function :   kth_smallest()
In       :   array of elements, # of elements in the array, rank k
Out      :   one element
Job      :   find the kth smallest element in the array
Notice   :   use the medianWirth() macro to get the median. 

Reference:

Author: Wirth, Niklaus 
Title: Algorithms + data structures = programs 
Publisher: Englewood Cliffs: Prentice-Hall, 1976 
Physical description: 366 p. 
Series: Prentice-Hall Series in Automatic Computation 
   
--------------------------------------------------------------------------------
*/


float kthSmallest(float a[], int n, int k);
#define medianWirth(a,n) kthSmallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

double kthSmallestDouble(double a[], int n, int k) ;
#define medianDouble(a,n) \
kthSmallestDouble(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))



VimosBool setupInterpolation(double **kernel, VimosLong32 **leaps, 
                             int imageXlen);


/*----------------------------------------------------------------------------
   Function :   opt_med3()
   In       :   pointer to array of 3 pixel values
   Out      :   a pixelvalue
   Job      :   optimized search of the median of 3 pixel values
   Notice   :   found on sci.image.processing
                cannot go faster unless assumptions are made
                on the nature of the input signal.
				MODIFIES THE INPUT ARRAY
 ---------------------------------------------------------------------------*/

float opt_med3(float  *p) ;
float opt_med5(float  *p) ;
float opt_med7(float  *p) ;
float opt_med9(float  *p) ;

/*---------------------------------------------------------------------------
   Function	:	medianPixelvalue() (extracted from ECLIPSE)
   In 		:	allocated array of pixelvalues, # of pixels in the array
   Out 		:	1 pixel value
   Job		:	compute the median pixel value out of an array
   Notice	:	calls the fastest method depending on the number of
			elements in input.
				MODIFIES THE INPUT ARRAY
 ---------------------------------------------------------------------------*/



float medianPixelvalue(float * a, int n);

double medianPixelvalueDouble(double * a, int n);

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void sort(int arrayLen, float *array)

  Description:
  Sorts an array of floats.
  I took this code from Numerical Recipes. NR has 1-based arrays and for
  Virmos we use 0-based. I modified the code for this brute force (i.e. not
  elegantly...). I tested the routine, seems to work...

  Input:
  int arrayLen
  Length of the input array
  
  Input/Output:
  float *array
  Array of floats to sort. Sorting is inplace.

  Return Value:
  void

  Updates:
  18 Apr 99: Created (TAO)

--------------------------------------------------------------------------------
*/
void sort(int n, float *ra);

/*
 * Calculate average of n values (avoiding overflows)
 */

double computeAverageInt(int a[], int n);
double computeAverageFloat(float a[], int n);
double computeAverageDouble(double a[], int n);

/*
 * Calculate variance avoiding trends in image
 */

double computeVarianceFloat2D(float a[], int nx, int ny);

/*
 * Calculate variance avoiding trends in image
 */

double computeVarianceDouble2D(double a[], int nx, int ny);

/*
 * This is for extracting a submatrix from a matrix of sizex X sizey
 * elements, starting from element x,y till element x+nx-1,y+ny-1
 */

float *extractFloatImage
         (float a[], int sizex, int sizey, int x, int y, int nx, int ny);

/*
 * This is for inserting a submatrix into a matrix of sizex
 * times sizey elements, starting from element x,y till element
 * x+nx-1,y+ny-1
 */

VimosBool insertFloatImage
(float a[], int sizex, int sizey, int x, int y, int nx, int ny, float region[]);

/*
 * Compute histogram of a double array: input the array, the array dimension, 
 * and the number of histogram bins
 */

VimosDpoint *darrayHistogram(double *darray, unsigned int arDim, 
			     unsigned int nbins);
/*
 * Compute histogram of double array giving the value of the starting and
 * of the ending bins, and the bin size
 */

VimosDpoint *darrayHistoStartEnd(double *darray, int arrDim, double start,
				 double end, double bin_size);

 /* find the position of the peak and the FWHM of an histogram. 
  * The peak position
  * depends on the  bin size (i.e on the number of bins "npoints")
  * of the histogram. 
  */

double histogramPeak(VimosDpoint *histogram,double *fwhm, unsigned int nbins);
/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  int  waterShed(float *profile, int numPoints, int numLevels, 
                 float fluxLimit, int width, int *mask);

  Description:
  Performs watershed detection on *profile. Idea taken from SExtractor. 

  Input:
  float *profile
  Pointer to float array to process
  
  int numPoints
  length of pointer[] and of mask[]

  int numLevels
  Number of levels to use in watershed thresholding

  float fluxLimit
  Flux limit to use in detection. Only objects with flux larger than
  fluxFraction * total flux will be considered as objects.
  
  int width
  Width of boxcar smoothing applied to profile[]

  Input/Output:
  int *mask
  Pointer to array that will contain the mask. The 

  Return Value (succes):
  number of objects detected

  Return Value (error):
  -1
  
  Updates:
  15 Jun 00: Return -1 on error (Maura)
  21 Apr 99: Created (TAO)

--------------------------------------------------------------------------------
*/
int  waterShed(float *profile, int numPoints, int numLevels, 
               float fluxLimit, int width, int *mask);

VimosBool findPeak1D(float *, int size, float *position, int minPoints);
VimosBool findPeak2D(float *, int sizex, int sizey, float *posx, float *posy, 
                     int minPoints);
VimosBool findDip1D(float *, int size, float *position, int minPoints);
VimosBool findUpJump(float *, int size, float *position, int minPoints);
VimosBool findJump(float *, int size, float *position, int minPoints);
VimosBool findDownJump(float *, int size, float *position, int minPoints);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  These are functions taken from the Rostat program written by T.Beers
  to compute robust statistics on a data set. See Beers, Flynn & Gebhardt
  1990, AJ 100, 32
-------------------------------------------------------------------------*/
float median(float x[],int n);

void xbiwt(float xdata[],int n,float *xlbiwt,float *xsbiwt,float *xlbiwt1,
	   float *xsbiwt1);

float xmad(float xdata[],int n,float xmed);

float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);




/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  These are functions for integrating & c. (NR- like)

------------------------------------------------------------------------------
*/

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  float trapezeInt(float (*theFunc)(float, float[]), float coeff[], 
                   float a, float b, int n)

  Description:
  Trapeze integration (modified trapzd from NR, addded coeff[]). 
  When called with n = 1 the function returns the crudest estimate of the 
  integral in[a,b] of f(x)dx. Subsequent calls with n=2,3,... (in that 
  sequential order) will improve the accuracy by adding  2^(n-2)  additional 
  interior points.

  Input:
  float (*theFunc)(float, float[])] 
  Function that gives, for each x, an y value computed by means of some 
  coefficients (see e.g. evalYFit)to be integrated

  float coeff[]
  Coefficients to be used to define the function to be integrated (e.g. linear 
  alone, or linear + gaussian, see  evalLineFlux)

  float a
  Lower limit of the integration

  float b
  Upper limit of the integration

  int n
  Stage of refinement of an extended trapezoidal rule.

  Updates:
  03 Jul 00: Created (AZ)

------------------------------------------------------------------------------
*/

float trapezeInt(float (*theFunc)(float, float[]), float coeff[], 
                 float a, float b, int n);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void polint(float xa[], float ya[], int n, float x, float *y, float *dy)

  Description:
  Polinomial interpolation (polint from NR). Used by rombergInt

  Input:
  float xa[]
  Input x array

  float ya[]
  Input y array

  int n
  Dimension of arrays x and y

  float x
  The x at which one wants to interpolate

  float *y
  the returned interpolated y corresponding to x

  float *dy
  Error estimate

  Updates:
  03 Jul 00: Created (AZ)

------------------------------------------------------------------------------
*/

void polint(float xa[], float ya[], int n, float x, float *y, float *dy);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  float rombergInt(float (*theFunc)(float, float[]), float coeff[], 
                   float a, float b)

  Description:
  Romberg integration method (modified qromb from NR, added coeffs[]).

  Input:
  float (*theFunc)(float, float[])
  Function that gives, for each x, an y value computed by means of some 
  coefficients (see e.g. evalYFit)

  float coeff[]
  Coefficients to be used to define the function to be integrated (e.g. linear 
  alone, or linear + gaussian, see evalLineFlux)

  float a
  Lower limit for integration

  float b
  Upper limit for integration

  Updates:
  03 Jul 00: Created (AZ)

------------------------------------------------------------------------------
*/

float rombergInt(float (*theFunc)(float, float[]), float coeff[], 
                 float a, float b);



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  float integrateSpectrum(VimosFloatArray *tmpSpec, float wlenStep)

  Description:
  This function performs a VERY simple integration, by summing all data
  using the formula:
  integrFlux += (tmpSpec->data[i] * wlenStep)
  To be improved!

  Input:
  VimosFloatArray *tmpSpec
  Array of input data to be integrated

  float wlenStep
  Integration step. Usually: wavelength increment in spectra taken from 
  descriptors

  Updates:
  25 Jul 00: Created (AZ)

------------------------------------------------------------------------------
*/

float integrateSpectrum(VimosFloatArray *tmpSpec, float wlenStep);



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  End of functions for integrating & c. (NR- like)

------------------------------------------------------------------------------
*/



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void Indexx(int n, float *arrin, int *indx)

   Description:
   Construct an index table to sort an array. Taken from Num. Rec. , changing 
   arrays from 1-based to 0-based.
   Needed for IFU: when sorting in fiber flux, to avoid losing memory of 
   which  fiber corresponds to which flux. Gives an array of subscripts to 
   sort the fiber numbers too.

   Input: 
   int n
   Length of input array

   float arrin
   array for which you want the index table
   
   Output:
   int indx
   Array containing the subscripts to sort arrin in ascendin order

   Return Value:
   void
   
   Updates:
   09 Feb 00: Created (AZ)
------------------------------------------------------------------------------
*/

void Indexx(int n, float *arrin, int *indx);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  float computeRMS(float *theData, int n)

  Description:
  Computes root mean square of theData. NB: variance is computed as the
  corrected, not-distorted one, i.e. divide by (n-1)

  Input:
  float *theData
  Array of data for which to compute RMS

  int n
  Number of points in theData

  Updates:
  25 Jul 00: Created (AZ)

------------------------------------------------------------------------------
*/

float computeRMS(float *theData, int n);


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  int  computeHistogram(VimosFloatArray *theData, int theXaxLen,
                        VimosFloatArray *theHisto, float mini, 
			float maxi, float step)

  Description:
  Computes the histogram of values in a given input array, between a minimum 
  and maximum and with a given bin size.

  Input:
  VimosFloatArray *theData
  Array of input values.

  int theXaxLen
  Length of Array of x axis (computed between mini and maxi, with step) .

  VimosFloatArray *theHisto
  Array of histogram values (defined outside, filled here)

  float mini
  Minimum value for computing the histogram

  float maxi
  Maximum value for computing the histogram

  float step
  Histogram bin size

  Updates:
  25 Jul 00: Created (AZ)

------------------------------------------------------------------------------
*/


int computeHistogram(VimosFloatArray *theData, int theXaxLen,
		     VimosFloatArray *theHisto, float mini, 
		     float maxi, float step);

PIL_END_DECLS

#endif /* VM_MATH_H */
