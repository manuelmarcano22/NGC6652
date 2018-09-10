/* $Id: vmmath.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmmath.h"
#include "cpl.h"


/* private functions */
static double *generateInterpolationKernel(char *kernelType) ;
static double *generateTanhKernel(double steep) ;
static void reverseTanhKernel(double *data, int nn) ;
static double sinc(double x) ;


/*---------------------------------------------------------------------------
   Function	:	ipow()
   In 		:	double, int
   Out 		:	double
   Job		:	same as pow(x,y) but for integer values of y
   Notice	:	much faster than the math function due to the integer.
   				On some compilers, this optimization is already done in
				the pow() function, we do not rely on that fact.
 ---------------------------------------------------------------------------*/


double ipow(double x, int i)
{
  double ret ;
  
  if (i == 0) {
    ret = 1.00;
  } 
  else {
    if (i == 1) {
      ret = x ;
    } 
    else {
      ret = x ;
      while (--i) {
        ret *= x;
      }
    }
  }
  
  return(ret);
}


/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :  	use the median_WIRTH() macro to get the median. 
   				MODIFIES THE INPUT ARRAY

   				Reference:

		          Author: Wirth, Niklaus 
				   Title: Algorithms + data structures = programs 
			   Publisher: Englewood Cliffs: Prentice-Hall, 1976 
	Physical description: 366 p. 
				  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/

#define PIX_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

float kthSmallest(float a[], int n, int k)
{
  register int i,j,l,m ;
  register float x ;
  
  l = 0; 
  m = n-1;
  while (l < m) {
    x = a[k];
    i = l;
    j = m;
    do {
      while (a[i] < x) {
        i++;
      }
      while (x < a[j]) {
        j--;
      }
      if (i <= j) {
        PIX_SWAP(a[i],a[j]);
        i++; 
        j--;
      }
    } while (i <= j);

    if (j < k) {
      l = i;
    }
    if (k < i) {
      m = j;
    }
    
  }
  return(a[k]);
}

#undef PIX_SWAP



/*---------------------------------------------------------------------------
   Function :   kth_smallest_double()
   In       :   array of doubles, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :  	use the median_double() macro to get the median. 
 ---------------------------------------------------------------------------*/


#define DBL_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double kthSmallestDouble(double a[], int n, int k)
{
    register int i,j,l,m;
    register double x;
    
    l = 0; 
    m = n-1;
    while (l < m) {
      x = a[k];
      i = l;
      j =m;
      do {
        while (a[i] < x) {
          i++ ;
        }
        while (x < a[j]) {
          j--;
        }
        if (i <= j) {
          DBL_SWAP(a[i],a[j]);
          i++; 
          j--;
        }
      } while (i <= j);
      if (j<k) {
        l = i;
      }
      if (k < i) {
        m = j;
      }
    }
    return(a[k]);
}

#undef DBL_SWAP



static double *generateInterpolationKernel(char *kernel_type)
{
  double  *tab;
  int     i;
  double  x;
  double  alpha;
  double  inv_norm;
  int     samples = KERNEL_SAMPLES ;

  if (!strcmp(kernel_type, "default")) {
    tab = generateInterpolationKernel("tanh") ;
    if (tab == NULL) {
      cpl_msg_error("generateInterpolationKernel","The function generateInterpolationKernel has returned NULL");
      return(NULL);
    }
    return(tab) ;
  } 
  
  if (!strcmp(kernel_type, "tanh")) {
    tab = generateTanhKernel(TANH_STEEPNESS) ;
    if (tab == NULL) {
      cpl_msg_error("generateInterpolationKernel","The function generateTanhKernel has returned NULL");
      return(NULL); 
    }
    return(tab) ;
  } 

  if (!strcmp(kernel_type, "sinc")) {
    tab = (double*) cpl_malloc(samples * sizeof(double)) ;

  /* check if space was allocated */
    if (tab == NULL) {
      cpl_msg_error("generateInterpolationKernel","Allocation Error");
      return(NULL);
    }

    tab[0] = 1.0 ;
    tab[samples-1] = 0.0 ;
    for (i = 1 ; i < samples ; i++) {
      x = (double)KERNEL_WIDTH * (double)i/(double)(samples-1) ;
      tab[i] = sinc(x) ;
    }
    return(tab) ;
  }  
  
  if (!strcmp(kernel_type, "sinc2")) {
    tab = (double*) cpl_malloc(samples * sizeof(double)) ;

  /* check if space was allocated */
    if (tab == NULL) {
      cpl_msg_error("generateInterpolationKernel","Allocation Error");
      return(NULL);
    }

    tab[0] = 1.0 ;
    tab[samples-1] = 0.0 ;
    for (i = 1 ; i < samples ; i++) {
      x = 2.0 * (double)i/(double)(samples-1) ;
      tab[i] = sinc(x) ;
      tab[i] *= tab[i] ;
    }
    return(tab) ;
  }  
  
  if (!strcmp(kernel_type, "lanczos")) {
    tab = (double*) cpl_malloc(samples * sizeof(double)) ;

  /* check if space was allocated */
    if (tab == NULL) {
      cpl_msg_error("generateInterpolationKernel","Allocation Error");
      return(NULL);
    }

    for (i = 0 ; i < samples ; i++) {
      x = (double)KERNEL_WIDTH * (double)i/(double)(samples-1) ;
      if (fabs(x)<2) {
        tab[i] = sinc(x) * sinc(x/2) ;
      } else {
        tab[i] = 0.00 ;
      }
    }
    return(tab) ;
  } 

  if (!strcmp(kernel_type, "hamming")) {
    tab = (double*) cpl_malloc(samples * sizeof(double)) ;

  /* check if space was allocated */
    if (tab == NULL) {
      cpl_msg_error("generateInterpolationKernel","Allocation Error");
      return(NULL);
    }

    alpha = 0.54 ;
    inv_norm  = 1.00 / (double)(samples - 1) ;
    for (i = 0 ; i < samples ; i++) {
      x = (double)i ;
      if (i < (samples-1)/2) {
        tab[i] = alpha + (1-alpha) * cos(2.0*PI_NUMB*x*inv_norm) ;
      } else {
        tab[i] = 0.0 ;
      }
    }
    return(tab) ;
  } 

  if (!strcmp(kernel_type, "hann")) {
    tab = (double*) cpl_malloc(samples * sizeof(double)) ;

  /* check if space was allocated */
    if (tab == NULL) {
      cpl_msg_error("generateInterpolationKernel","Allocation Error");
      return(NULL);
    }

    alpha = 0.50 ;
    inv_norm  = 1.00 / (double)(samples - 1) ;
    for (i = 0 ; i < samples ; i++) {
      x = (double)i ;
      if (i < (samples-1)/2) {
        tab[i] = alpha + (1-alpha) * cos(2.0*PI_NUMB*x*inv_norm) ;
      } else {
        tab[i] = 0.0 ;
      }
    }
    return(tab) ;
  } 



  cpl_msg_error("generateInterpolationKernel","Unrecognized kernel type [%s]: aborting generation", kernel_type) ;
  return(NULL);
}
  



/*---------------------------------------------------------------------------
 * Function :   sinc()
 * In       :   double
 * Out      :   double
 * Job      :   cardinal sine
 * Notice   :   rescaled by PI
 *--------------------------------------------------------------------------*/

static double sinc(double x)
{
if (fabs(x) < 1e-4)
     return (double)1.00 ;
     else
     return ((sin(x * (double)PI_NUMB)) / (x * (double)PI_NUMB)) ;
}



/*
 * The following function is a good approximation of a box filter,
 * built up from a product of hyperbolic tangents. It has the following
 * properties:
 * 1. It converges very quickly towards +/- 1
 * 2. The converging transition is sharp
 * 3. It is infinitely differentiable everywhere (smooth)
 * 4. The transition sharpness is scalable
 *
 * It is defined as a macro here (yes, it's dirty) to optimize the code
 */

#define hkGen(x,s) (((tanh(s*(x+0.5))+1)/2)*((tanh(s*(-x+0.5))+1)/2))


/*---------------------------------------------------------------------------
   Function	:	generate_tanh_kernel()
   In 		:	double: steepness of the hyperbolic tangent
   Out 		:	pointer to (samples) doubles
   Job		:	generate an hyperbolic tangent kernel
   Notice	:	don't touch
 ---------------------------------------------------------------------------*/


static double *generateTanhKernel(double steep)
{
  double  *   kernel ;
  double  *   x ;
  double      width ;
  double      inv_np ;
  double      ind ;
  int         i ;
  int         np ;
  int         samples ;
  
  width   = (double)TABSPERPIX / 2.0 ; 
  samples = KERNEL_SAMPLES ;
  np      = 32768 ; /* Hardcoded: should never be changed */
  inv_np  = 1.00 / (double)np ;
  
  /*
   * Generate the kernel expression in Fourier space
   * with a correct frequency ordering to allow standard FT
   */
  x = (double*) cpl_malloc((2*np+1)*sizeof(double)) ;

  /* check if space was allocated */
    if (x == NULL) {
      cpl_msg_error("generateTanhKernel","Allocation Error");
      return(NULL);
    }

  for (i = 0 ; i < np/2 ; i++) {
    ind      = (double)i * 2.0 * width * inv_np ;
    x[2*i]   = hkGen(ind, steep) ;
    x[2*i+1] = 0.00 ;
  }
  for (i = np/2 ; i < np ; i++) {
    ind      = (double)(i-np) * 2.0 * width * inv_np ;
    x[2*i]   = hkGen(ind, steep) ;
    x[2*i+1] = 0.00 ;
  }
  
  /* 
   * Reverse Fourier to come back to image space
   */
  reverseTanhKernel(x, np) ;
  
  /*
   * Allocate and fill in returned array
   */
  kernel = (double*) cpl_malloc(samples * sizeof(double)) ;

  /* check if space was allocated */
    if (kernel == NULL) {
      cpl_msg_error("generateTanhKernel","Allocation Error");
      return(NULL);
    }

  for (i = 0 ; i < samples ; i++) {
    kernel[i] = 2.0 * width * x[2*i] * inv_np ;
  }
  cpl_free(x) ;
  return kernel ;
}


/*---------------------------------------------------------------------------
   Function	:	reverse_tanh_kernel()
   In 		:	a tanh generated kernel in Fourier space
   Out 		:	a reversed kernel in image space
   Job		:	transforms from Fourier to image space a tanh kernel
   Notice	:	optimized, only useful for generate_tanh_kernel()
 ---------------------------------------------------------------------------*/


#define KERNEL_SW(a,b) tempr=(a);(a)=(b);(b)=tempr
static void reverseTanhKernel(double *data, int nn)
{
  VimosUlong32   n,
    mmax,
    m,
    i, j,
    istep ;
  double  wtemp,
    wr,
    wpr,
    wpi,
    wi,
    theta;
  double  tempr,
    tempi;
  
  n = (VimosUlong32)nn << 1;
  j = 1;
  for (i=1 ; i<n ; i+=2) {
    if (j > i) {
      KERNEL_SW(data[j-1],data[i-1]);
      KERNEL_SW(data[j],data[i]);
    }
    m = n >> 1;
    while (m>=2 && j>m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = 2 * PI_NUMB / mmax;
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr  = 1.0;
    wi  = 0.0;
    for (m=1 ; m<mmax ; m+=2) {
      for (i=m ; i<=n ; i+=istep) {
        j = i + mmax;
        tempr = wr * data[j-1] - wi * data[j];
        tempi = wr * data[j]   + wi * data[j-1];
        data[j-1] = data[i-1] - tempr;
        data[j]   = data[i]   - tempi;
        data[i-1] += tempr;
        data[i]   += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}
#undef KERNEL_SW


VimosBool setupInterpolation(double **kernel, VimosLong32 **leaps, 
                             int imageXlen)
{
  VimosLong32 *tmpLong;
  
  *kernel = generateInterpolationKernel("default");
  if (*kernel == NULL) {
    cpl_msg_error("setupInterpolation","The function generateInterpolationKernel has returned NULL");
    return(VM_FALSE);
  }
  
  tmpLong = (VimosLong32 *) cpl_malloc(16*sizeof(VimosLong32));

  /* check if space was allocated */
    if (tmpLong == NULL) {
      cpl_msg_error("setupInterpolation","Allocation Error");
      return(VM_FALSE);
    }
  
  tmpLong[0] = -1 - imageXlen ;
  tmpLong[1] =    - imageXlen ;
  tmpLong[2] =  1 - imageXlen ;
  tmpLong[3] =  2 - imageXlen ;
  
  tmpLong[4] = -1 ;
  tmpLong[5] =  0 ;
  tmpLong[6] =  1 ;
  tmpLong[7] =  2 ;
  
  tmpLong[8] = -1 + imageXlen ;
  tmpLong[9] =      imageXlen ;
  tmpLong[10]=  1 + imageXlen ;
  tmpLong[11]=  2 + imageXlen ;
  
  tmpLong[12]= -1 + 2*imageXlen ;
  tmpLong[13]=      2*imageXlen ;
  tmpLong[14]=  1 + 2*imageXlen ;
  tmpLong[15]=  2 + 2*imageXlen ;

  *leaps = tmpLong;
  
  return(VM_TRUE);
}


void sort(int n, float *ra)
{
  int l,j,ir,i;
  float rra;
  
  if (n == 1)
      return;
    
  l=(n >> 1)+1;
  ir=n;
  for (;;) {  
    if (l > 1) {
      --l;
      rra=ra[l-1];
    } else {   
      rra=ra[ir-1];
      ra[ir-1]=ra[0];
      if (--ir == 1) {
        ra[0]=rra;
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir) {
      if (j < ir && ra[j-1] < ra[j])
        ++j;  
      if (rra < ra[j-1]) {
        ra[i-1]=ra[j-1];
        j += (i=j);
      } else 
        j=ir+1;
    }
    ra[i-1]=rra;
  }

}

double computeAverageInt(int a[], int n)
{
  int    i;
  double ave = 0.;
  char   modName[] = "computeAverageInt";

  if (n > 0) {
    for (i=0; i<n; i++) {
     /*
      * This is a bit slow but safer against overflows
      */
      ave = (((double) i)/(i+1))*ave + ((double) a[i])/(i+1);
    }
  }
  else {
    cpl_msg_error(modName, "Array size must be positive");
  }
  return(ave);
}

double computeAverageFloat(float a[], int n)
{
  int    i;
  double ave = 0.;
  char   modName[] = "computeAverageFloat";

  if (n > 0) {
    for (i=0; i<n; i++) {
     /*
      * This is a bit slow but safer against overflows
      */
      ave = (((double) i)/(i+1))*ave + ((double) a[i])/(i+1);
    }
  }
  else {
    cpl_msg_error(modName, "Array size must be positive");
  }
  return(ave);
}

double computeAverageDouble(double a[], int n)
{
  int    i;
  double ave = 0.;
  char   modName[] = "computeAverageDouble";

  if (n > 0) {
    for (i=0; i<n; i++) {
     /*
      * This is a bit slow but safer against overflows
      */
      ave = (((double) i)/(i+1))*ave + a[i]/(i+1);
    }
  }
  else {
    cpl_msg_error(modName, "Array size must be positive");
  }
  return(ave);
}

double computeVarianceFloat2D(float a[], int nx, int ny)
{

  int    x, y;
  int    n = 0;
  int    dx = 1;
  int    dy = 1;
  double delta;
  double variance = 0.0;

  if (nx < 4 || ny < 4)
    return 0;

  nx--;
  ny--;

  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      delta = a[x + y * nx] - a[x + dx + (y + dy) * nx];
      variance = (((double)n)/(n+1))*variance + ((double)delta*delta)/(n+1);
      n++;
    }
  }

  return variance / 2.0;
  
}

double computeVarianceDouble2D(double a[], int nx, int ny)
{

  int    x, y;
  int    n = 0;
  int    dx = 1;
  int    dy = 1;
  double delta;
  double variance = 0.0;

  if (nx < 4 || ny < 4)
    return 0;

  nx--;
  ny--;

  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      delta = a[x + y * nx] - a[x + dx + (y + dy) * nx];
      variance = (((double)n)/(n+1))*variance + ((double)delta*delta)/(n+1);
      n++;
    }
  }

  return variance / 2.0;
  
}


/*************************************************************************/

float *extractFloatImage
  (float a[], int sizex, int sizey, int x, int y, int nx, int ny)
{
  float *extractedImage;
  int    i,j;
  char   modName[] = "extractFloatImage";

  if (x >= 0 && y >= 0 && sizex >= x + nx && sizey >= y + ny) {
    extractedImage = (float *) cpl_malloc(nx*ny*sizeof(float));
    for (i=0; i<ny; i++) {
      for (j=0; j<nx; j++) {
        extractedImage[nx*i+j] = a[sizex*(y+i)+x+j];
      }
    }
    return(extractedImage);
  }
  else {
    cpl_msg_error(modName, "Extracted image is not contained in source image");
    return NULL;
  }
}

VimosBool insertFloatImage(float a[], int sizex, int sizey,
int x, int y, int nx, int ny, float region[])
{
  int    i,j;
  char   modName[] = "insertFloatImage";

  if (x >= 0 && y >= 0 && sizex >= x + nx && sizey >= y + ny)
  {
    for (i=0; i<ny; i++)
    {
      for (j=0; j<nx; j++)
      {
        a[sizex*(y+i)+x+j] = region[nx*i+j];
      }
    }
  }
  else
  {
    cpl_msg_error(modName, "Extracted image is not contained in source image");
    return(VM_FALSE);
  }
  return(VM_TRUE);
}

VimosDpoint *darrayHistogram(double *darray, unsigned int arDim,
			     unsigned int nbins)

{
  float           max,min;
  unsigned int    i ;
  unsigned int   *h ;      /* frequency of each bin_val */
  int             bin_val ;
  double          bin_size ;
  VimosDpoint    *histogram ;

  max = darray[0];
  min = darray[0];
  for(i=0; i<arDim; i++) {
    if(darray[i] < min)  min = darray[i];
    if(darray[i] > max)  max = darray[i];
  }

  bin_size = (double)(max - min)/(double)nbins;
  h = (unsigned int*)cpl_calloc(nbins,sizeof(unsigned int)) ;
 
  for (i=0 ; i<arDim ; i++) { 
    if (darray[i] >= max) {
      bin_val = nbins - 1;
    }
    else {
     bin_val = (int)((darray[i] - min) /bin_size) ;
    }
    h[bin_val]++ ;
  }
  histogram = newDpoint(nbins) ;
  for (i=0 ; i<nbins ; i++) {
     histogram[i].x = min + (double)i * bin_size ;
     histogram[i].y = (double)h[i] ;
  }
  cpl_free(h) ;
  return histogram ;
}

VimosDpoint *darrayHistoStartEnd(double *darray, int arrDim, double start,
				double end, double bin_size)
{
  int    i ;
  unsigned int   *h ;      /* frequency of each bin_val */
  int             nbins, bin_val=0;
  VimosDpoint    *histogram ;
  char            modName[] = "farrayHistoStartEnd";

  if(start > end) {
    cpl_msg_error(modName, "start point must be lower than end point");
    return NULL;
  }

  nbins = (floor)((start-end)/bin_size);
  h = (unsigned int*)cpl_calloc(nbins,sizeof(unsigned int)) ;
 
  for (i=0 ; i<arrDim ; i++) { 
    if(((int)(darray[i] - start) > 0) && ((int)(end - darray[i]) > 0))
      bin_val = (int)((darray[i] - start) /bin_size) ;
    h[bin_val]++ ;
  }
  histogram = newDpoint(nbins) ;
  for (i=0 ; i<nbins ; i++) {
    histogram[i].x = start + (double)i * bin_size ;
    histogram[i].y = (double)h[i] ;
  }
  cpl_free(h) ;
  return histogram ;
}

double histogramPeak(VimosDpoint *histogram, double *fwhm, unsigned int nbins)
{
  unsigned int          i, left, right, peakind;
  double       peakval,peakpos,peak,hpeak,left_x,right_x;
  char         modName[] = "histogramPeak";

  peakind=0;
  peak   = histogram[0].y;
  for (i=1; i<nbins; i++) {
    if (histogram[i].y > peak) {
      peak    = histogram[i].y;
      peakind = i;
    }
  }
  peakpos= histogram[peakind].x;
  peakval= histogram[peakind].y;

  hpeak = peakval * 0.5;
  left  = peakind -1 ;
  right = peakind +1 ;
  while (histogram[left].y > hpeak) left-- ;
  while ((histogram[right].y > hpeak) && (right < (nbins-1))) right++ ;
  if ((left == 0) || 
      (right == nbins-1)) {
    cpl_msg_error(modName, "Cannot compute histogram FWHM") ;
    return -1. ;
  }
  left_x = histogram[left].x + (hpeak - histogram[left].y) *
    (histogram[left+1].x - histogram[left].x) /
    (histogram[left+1].y - histogram[left].y) ;

  right_x = histogram[right-1].x +
    (hpeak - histogram[right-1].y) *
    (histogram[right].x - histogram[right-1].x) /
    (histogram[right].y - histogram[right-1].y) ;

  *fwhm = fabs(right_x - left_x) ;
  return peakpos;
}

int  waterShed(float *profile, int numPoints, int numLevels, 
               float fluxLimit, int width, int *mask)
{

  float *tmpProf;
  int level;
  int i, j;
  int numCand;
  int ind1, ind2;
  int candStart, candEnd;
  int numObjects;
  
  float totalFlux;
  float objectFlux;
  
  float minimum, maximum;
  float sum, num;

  /* copy input profile and box smooth */
  tmpProf = (float *) cpl_malloc(numPoints*sizeof(float));

  /* check if space was allocated */
    if (tmpProf == NULL) {
      cpl_msg_error("waterShed","Allocation Error");
      return(-1);
    }

  for (i = 0; i < numPoints; i++) {
    sum = 0.0;
    num = 0.0;
    
    for (j = MAX(i-width, 0); j < MIN(i+width, numPoints-1); j++) {
      sum += profile[j];
      num = num+1.0;
    }
    tmpProf[i] = sum/num;
  }

  /* 
     Get minimum and maximum.
     We do not treat the edges, so we start at pixel 1
  */
  minimum = maximum = tmpProf[1];
  for (i = 1; i < numPoints-1; i++) {
    minimum = MIN(minimum, tmpProf[i]);
    maximum = MAX(maximum, tmpProf[i]);
  }

  /* in the unlikely case that the input is constant, the whole profile is
     object */
  if (minimum == maximum) {
    for (i = 0; i < numPoints; i++) {
      mask[i] = 1;
      return(1);
    }
  }

  /* scale profile into levels range
     get total flux to use in detection criterion */
  totalFlux = 0.0;
  for (i = 1; i < numPoints-1; i++) {
    totalFlux += tmpProf[i]/(maximum-minimum)*numLevels;
    tmpProf[i] = (tmpProf[i]-minimum)/(maximum-minimum) * numLevels;
  }
  /* initialise output mask */
  for (i = 0; i < numPoints; i++) {
    mask[i] = 0;
  }
  /* initialise counter of number of objects */
  numObjects = 0;
  
  /* loop through levels, top->down */
  for (level = numLevels-1; level > 0; level--) {
    /* initialise number of candidate new objects */
    numCand = -1;
    /* loop through profile */
    for (i = 1; i < numPoints-1; i++) {
      /* if not yet object.. */
      if (mask[i] == 0) {
        /* check if value is above current level */
        if (tmpProf[i] > level) {
          /* if it is, check if this pixel is contiguous with previous pixel */
          if (mask[i-1] != numCand) {
            /* if not, this is a new candidate, so increment counter */
            numCand--;
          }
          /* 
             set mask for this pixel to candidate number.
             numCand is negative to flag that his is only a candidate, not yet
             an accepted object 
          */
          mask[i] = numCand;
        }
      }
    }

    /* loop over new candidates */
    for (j = 2; j <= -numCand; j++) {
      /* set index to first pixel in profile that we consider */
      i = 1;
      /* find start of candidate j (actually candidate j-1...) */
      while (mask[i] != -j) {
        i++;
      }
      candStart = i;


      /* initialise flux of candidate */
      objectFlux = 0.0;      
      /* find end of range of candidate j, and integrate */
      while ( (i < numPoints-1) && (mask[i] == -j) ) {
        objectFlux += (tmpProf[i]-level)/totalFlux;
        i++;
      }
      candEnd = i-1;
      /* check if flux of this candidate is above threshold
       AND not touching any other real object. That is what I think a
       watershed algorithm requires and should work if the number of levels is
       high enough. The critetion of objectFlux > fluxLimit fails for low S/N
       profiles. Obviously the criterion should be S/N dependent. 
       Note that a candidate can only touch real objects, not candidates 
      */
      if ( (objectFlux > fluxLimit) &&  
           (mask[candStart-1] == 0) &&
           (mask[candEnd + 1] == 0) ) {
        /* seems we have a new object */
        numObjects++;
        /* flag it as such */
        for (i = candStart; i <= candEnd; i++) {
          mask[i] = numObjects;
        }
      } 
      /* not strong enough, should we merge? */
      else {
        /* real object  to the left? */
        if (mask[candStart-1] > 0) {
          /* only object to the left? */
          if (mask[candEnd+1] == 0) {
            /* merge with object to the left */
            for (i = candStart; i <= candEnd; i++) {
              mask[i] = mask[candStart-1];
            }
          } 
          else {
            /* in between two real objects.
               object to the right has to be real, since two candidates cannot
               touch */
            ind1 = candStart;
            ind2 = candEnd;
            /* loop from left and right */
            while(ind1 < ind2) {
              /* merge left part with left object, right part with right object */
              mask[ind1] = mask[candStart-1];
              mask[ind2] = mask[candEnd+1];
              /* increment left index, decrement right index */
              ind1++;
              ind2--;
            }
            /* odd pixel stuff. Assign central pixel to "heaviest" */
            if (ind1 == ind2) {
              if (tmpProf[ind1-1] > tmpProf[ind1+1]) {
                mask[ind1] = mask[candStart-1];
              }
              else {
                mask[ind2] = mask[candEnd+1];
              }
            }
          }
        }
        else {
          /* only real object to the right ? */
          if (mask[candEnd+1] > 0) {
            /* merge with object to the right */
            for (i = candStart; i <= candEnd; i++) {
              mask[i] = mask[candEnd+1];
            }
          }
          else {
            for (i = candStart; i <= candEnd; i++) {
              mask[i] = 0;
            }
          }

        }
      }
        
    }
      
  }
  
  return(numObjects);
  
  
  cpl_free(tmpProf);
  
}




#define PIX_SORT(a,b) { if ((a)>(b)) PIX_SWAP((a),(b)); }
#define PIX_SWAP(a,b) { float temp=(a);(a)=(b);(b)=temp; }

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

float opt_med3(float  *p)
{
    PIX_SORT(p[0],p[1]) ; PIX_SORT(p[1],p[2]) ; PIX_SORT(p[0],p[1]) ;
    return(p[1]) ;
}

/*----------------------------------------------------------------------------
   Function :   opt_med5()
   In       :   pointer to array of 5 pixel values
   Out      :   a pixelvalue
   Job      :   optimized search of the median of 5 pixel values
   Notice   :   found on sci.image.processing
                cannot go faster unless assumptions are made
                on the nature of the input signal.
				MODIFIES THE INPUT ARRAY
 ---------------------------------------------------------------------------*/

float opt_med5(float  *p)
{
    PIX_SORT(p[0],p[1]) ; PIX_SORT(p[3],p[4]) ; PIX_SORT(p[0],p[3]) ;
    PIX_SORT(p[1],p[4]) ; PIX_SORT(p[1],p[2]) ; PIX_SORT(p[2],p[3]) ;
    PIX_SORT(p[1],p[2]) ; return(p[2]) ;
}

/*----------------------------------------------------------------------------
   Function :   opt_med7()
   In       :   pointer to array of 7 pixel values
   Out      :   a pixelvalue
   Job      :   optimized search of the median of 7 pixel values
   Notice   :   found on sci.image.processing
                cannot go faster unless assumptions are made
                on the nature of the input signal.
				MODIFIES THE INPUT ARRAY
 ---------------------------------------------------------------------------*/

float opt_med7(float  *p)
{
    PIX_SORT(p[0], p[5]) ; PIX_SORT(p[0], p[3]) ; PIX_SORT(p[1], p[6]) ;
    PIX_SORT(p[2], p[4]) ; PIX_SORT(p[0], p[1]) ; PIX_SORT(p[3], p[5]) ;
    PIX_SORT(p[2], p[6]) ; PIX_SORT(p[2], p[3]) ; PIX_SORT(p[3], p[6]) ;
    PIX_SORT(p[4], p[5]) ; PIX_SORT(p[1], p[4]) ; PIX_SORT(p[1], p[3]) ;
    PIX_SORT(p[3], p[4]) ; return (p[3]) ;
}

/*----------------------------------------------------------------------------
   Function :   opt_med9()
   In       :   pointer to an array of 9 pixelvalues
   Out      :   a pixelvalue
   Job      :   optimized search of the median of 9 pixelvalues
   Notice   :   in theory, cannot go faster without assumptions on the
                signal.
                Formula from:
                XILINX XCELL magazine, vol. 23 by John L. Smith
  
                The input array is modified in the process
                The result array is guaranteed to contain the median
                value
                in middle position, but other elements are NOT sorted.

				MODIFIES THE INPUT ARRAY
 ---------------------------------------------------------------------------*/

float opt_med9(float  *p)
{
    PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ;
    PIX_SORT(p[0], p[1]) ; PIX_SORT(p[3], p[4]) ; PIX_SORT(p[6], p[7]) ;
    PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ;
    PIX_SORT(p[0], p[3]) ; PIX_SORT(p[5], p[8]) ; PIX_SORT(p[4], p[7]) ;
    PIX_SORT(p[3], p[6]) ; PIX_SORT(p[1], p[4]) ; PIX_SORT(p[2], p[5]) ;
    PIX_SORT(p[4], p[7]) ; PIX_SORT(p[4], p[2]) ; PIX_SORT(p[6], p[4]) ;
    PIX_SORT(p[4], p[2]) ; return(p[4]) ;
}


#undef PIX_SORT
#undef PIX_SWAP


/*---------------------------------------------------------------------------
   Function :   medianPixelvalue()
   In       :   allocated array of pixelvalues, # of pixels in the array
   Out      :   1 pixel value
   Job      :   compute the median pixel value out of an array
   Notice   :   calls the fastest method depending on the number of
                elements in input.
                 
 ---------------------------------------------------------------------------*/



float medianPixelvalue(float * a, int n)
{
  float  medianVal ;
  float *a_copy;
  int    i;

  if(n==1) return a[0];

  a_copy = (float*) cpl_malloc (n*sizeof(float));
  for (i=0; i<n; i++) a_copy[i] = a[i];

    switch(n) {
    case 3:
      medianVal = opt_med3(a_copy);
      break ;

    case 5:
      medianVal = opt_med5(a_copy);
      break ;

    case 7:
      medianVal = opt_med7(a_copy);
      break ;

    case 9:
      medianVal = opt_med9(a_copy);
      break ;

    default:
      medianVal = ((n>1000) ? medianWirth(a_copy,n): median(a_copy,n));
      break ;
    }
    cpl_free(a_copy);
    return medianVal;
}

/*---------------------------------------------------------------------------
   Function :   medianPixelvalue()
   In       :   allocated array of pixelvalues, # of pixels in the array
   Out      :   1 pixel value
   Job      :   compute the median pixel value out of an array
   Notice   :   calls the fastest method depending on the number of
                elements in input.
                 
 ---------------------------------------------------------------------------*/



double medianPixelvalueDouble(double * a, int n)
{
  float  medianVal ;
  float *a_copy;
  int    i;

  if(n==1) return a[0];

  a_copy = (float*) cpl_malloc (n*sizeof(float));
  for (i=0; i<n; i++) a_copy[i] = a[i];

    switch(n) {
    case 3:
      medianVal = opt_med3(a_copy);
      break ;

    case 5:
      medianVal = opt_med5(a_copy);
      break ;

    case 7:
      medianVal = opt_med7(a_copy);
      break ;

    case 9:
      medianVal = opt_med9(a_copy);
      break ;

    default:
      medianVal = ((n>1000) ? medianWirth(a_copy,n): median(a_copy,n));
      break ;
    }
    cpl_free(a_copy);
    return medianVal;
}

VimosBool findPeak1D(float *data, int size, float *position, int minPoints) 
{
  int    i;
  int    count = 0;
  float *copy;
  float  max, median, level, pos, variance, uniformVariance;
  double sum, weights;

  if (data == NULL) return VM_FALSE;
  if (size < 5) return VM_FALSE;         /* Hardcoded, I know... */
 /*
  *  Find median
  */
  copy = (float *) cpl_malloc(size*sizeof(float));
  for (i=0; i<size; i++) copy[i] = data[i];
  median = medianWirth(copy, size);
  cpl_free(copy);

 /*
  *  Find max
  */
  max = data[0];
  for (i=1; i<size; i++) if (data[i] > max) max = data[i];

 /*
  *  If the max equals the median we have a flat input, therefore
  *  no peak is found.
  */

  if (max-median < MIN_DIVISOR) return VM_FALSE;

 /*
  *  Discrimination level: only pixels with values above this
  *  level are considered in baricenter calculation.
  */
  level = (max+median)/2;

 /*
  *  Of the values above this level compute the baricenter and
  *  then the variance of the positions used. Note that the weights
  *  are taken as the difference between the pixels values and 
  *  the median level (supposedly the background).
  */
  count = 0;
  for (i=0, sum=0., weights=0.; i<size; i++) {
    if (data[i] > level) {
      count++;
      weights += (data[i] - median);
      sum     += i*(data[i] - median);
    }
  }
 /*
  *  If too few values are above threshold, refuse the position
  *  as insignificant
  */
  if (count < minPoints) return VM_FALSE;

  pos = sum/weights;
  for (i=0, sum=0., weights=0.; i<size; i++) {
    if (data[i] > level) {
      weights++;
      sum += (i-pos)*(i-pos);
    }
  }
  variance = sqrt(sum/weights);

 /*
  *  The "uniform variance" is the variance that should be obtained
  *  in the case of uniform distribution of the points positions in
  *  the selected interval. If the real variance is comparable with
  *  this value, the peak is considered not found.
  */

  uniformVariance = sqrt(size*size/3 - pos*size + pos*pos);

  if (variance > 0.8 * uniformVariance) return VM_FALSE;
  *position = pos;
  return VM_TRUE;
}

VimosBool findPeak2D(float *data, int sizex, int sizey, 
                     float *posx, float *posy, int minPoints)
{
  int    i, j;
  int    count = 0;
  int    size = sizex*sizey;
  float *copy;
  float  max, median, level, xpos, ypos;
  float  levelSigma, levelFraction;
  float  variancex, variancey, uniformVariancex, uniformVariancey;
  float  value;
  double sumx, sumy, sum, weights;

  if (data == NULL) return VM_FALSE;
  if (sizex < 5 || sizey < 5) return VM_FALSE; /* Hardcoded, I know... */

 /*
  *  Find median
  */
  copy = (float *) cpl_malloc(size*sizeof(float));
  for (i=0; i<size; i++) copy[i] = data[i];
  median = medianWirth(copy, size);
  cpl_free(copy);

 /*
  *  Find max
  */
  max = data[0];
  for (i=1; i<size; i++) if (data[i] > max) max = data[i];

 /*
  *  If the max equals the median we have a flat input, therefore
  *  no peak is found.
  */

  if (max-median < MIN_DIVISOR) {
    return VM_FALSE;
  }

 /*
  *  Discrimination level: only pixels with values above this
  *  level are considered in baricenter calculation.
  */
  levelFraction = (max+3*median)/4; 

  sum = count = 0;
  for (j=0; j<sizey; j++) {
    for (i=0; i<sizex; i++) {
      value = data[i+j*sizex];
      if ((median - value) > 0) {
        count++;
        sum += (median-value)*(median-value);
      }
    }
  }
  levelSigma = median + 3*sqrt(sum/count);
  level = MAX(levelFraction, levelSigma);
 /*
  *  Of the values above this level compute the baricenter and
  *  then the variance of the positions used. Note that the weights
  *  are taken as the difference between the pixels values and
  *  the median level (supposedly the background).
  */
  count = 0;
  sumx = sumy = weights = 0.;
  for (j=0; j<sizey; j++) {
    for (i=0; i<sizex; i++) {
      value = data[i+j*sizex];
      if (value > level) {
        count++;
        weights += (value - median);
        sumx    += i*(value - median);
        sumy    += j*(value - median);
      }
    }
  }
 /*
  *  If too few values are above threshold, refuse the position
  *  as insignificant
  */
  if (count < minPoints) {
    return VM_FALSE;
  }

  xpos = sumx/weights;
  ypos = sumy/weights;

  sumx = sumy = weights = 0.;
  for (j=0; j<sizey; j++) {
    for (i=0; i<sizex; i++) {
      if (data[i+j*sizex] > level) {
        weights++;
        sumx += (i-xpos)*(i-xpos);
        sumy += (j-ypos)*(j-ypos);
      }
    }
  }
  variancex = sqrt(sumx/weights);
  variancey = sqrt(sumy/weights);
  

 /*
  *  The "uniform variance" is the variance that should be obtained
  *  in the case of uniform distribution of the points positions in
  *  the selected interval. If the real variance is comparable with
  *  this value, the peak is considered not found.
  */

  uniformVariancex = sqrt(sizex*sizex/3 - xpos*sizex + xpos*xpos);
  uniformVariancey = sqrt(sizey*sizey/3 - ypos*sizey + ypos*ypos);

  if (variancex > uniformVariancex/2.) {
    return VM_FALSE;
  }
  if (variancey > uniformVariancey/2.) {
    return VM_FALSE;
  }
  *posx = xpos;
  *posy = ypos;

  return VM_TRUE;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  These are functions taken from the Rostat program written by T.Beers
  to compute robust statistics on a data set. See Beers, Flynn & Gebhardt
  1990, AJ 100, 32
-------------------------------------------------------------------------*/

/****************************************************************/
/**  This is just to compute the median of a set of points     **/
/****************************************************************/
float median(float x[],int n)
{
   float *xwork,medVal;
   int i,n2;
  
   xwork = (float *) cpl_malloc(n*sizeof(float));
  
   for (i = 0; i < n; i++)
     xwork[i] = x[i];
  
   sort(n, xwork);
   n2 = (int) (n/2);
   if (2*n2 == n)
     medVal = (0.5*(xwork[n2-1] + xwork[n2]));
   else
     medVal = xwork[n2];

   cpl_free(xwork);
   return medVal;
}


/*************************************************************************/
/**  The subroutine XBIWT provides an estimator of the location and     **/
/**  scale of the data set XDATA.  The scale uses the Biweight function **/
/**  in the general formula of "A-estimators." This formula is given    **/
/**  on page of 416 in UREDA (formula 4). The BIWEIGHT scale estimate   **/
/**  is returned as the value XSBIWT. The BIWEIGHT function is given    **/
/**  by:                                                                **/
/**                                                                     **/
/**                                u((1-u*u)**2)     abs(u) <= 1        **/
/**                       f(u) =                                        **/
/**                                0                 abs(u) >  1        **/
/**                                                                     **/
/**  where u is defined by                                              **/
/**                                                                     **/
/**                       u = (XDATA(I) - M) / c*MAD  .                 **/
/**                                                                     **/
/**  M, MAD, and c are the median, the median absolute deviation from   **/
/**  the median, and the tuning constant respectively. The tuning       **/
/**  constant is a parameter which is chosen depending on the sample    **/
/**  size and the specific function being used for the scale estimate.  **/
/**  (See page 417 in UREDA).  Here we take c = 9.0.                    **/
/**                                                                     **/
/**  The biweght location is found using the formula:                   **/
/**                                                                     **/
/**                       T = M + (sums)                                **/
/**                                                                     **/
/**                       where M is the sample median and sums are     **/
/**                       as given on page 421 in UREDA                 **/
/**                                                                     **/
/**               the tuning constant c is set to 6.0 for calculation   **/
/**               of the location as reccommended by Tukey ()           **/
/**                                                                     **/
/**  NOTE that the biweight is meant to be an iterated estimator, but   **/
/**  one commonly only takes the first step of the iteration.  Here we  **/
/**  report both the one-step estimators (XLBIWT1, XSBIWT1) and the     **/
/**  preferred fully iterated versions (XLBIWT, XSBIWT).                **/
/**                                                                     **/
/**  This function was taken from the library of ROSTAT functions       **/
/**  provided by T. Beers for his ROSTAT code.                          **/
/**  Ported to C by Marco Scodeggio, sometimes in 1998.                 **/ 
/*************************************************************************/
#define D1 1.0
#define D5 5.0
#define D6 6.0
#define D9 9.0

void xbiwt(float xdata[],int n,float *xlbiwt,float *xsbiwt,float *xlbiwt1,
       float *xsbiwt1)
{
   float *u1,*u2,*xlb,*xsb,xm,xmadm,c1,c2,xmm,s1,s2,s3,s4,*xwork;
   int i,j;

   u1 = (float *) cpl_malloc(n*sizeof(float));
   u2 = (float *) cpl_malloc(n*sizeof(float));
   xlb = (float *) cpl_malloc(11*sizeof(float));
   xsb = (float *) cpl_malloc(11*sizeof(float));
   xwork = (float *) cpl_malloc(n*sizeof(float));
  
   for (i = 0; i < n; i++)
     xwork[i] = xdata[i];

/**  find the median **/

   xm = median(xwork,n);

/**     call xmad to find the median absolute deviation **/

   xmadm = xmad(xwork,n,xm);

/**     must choose value of the tuning constant "c"   **/
/**     here c = 6.0 for the location estimator and    **/
/**     9.0 for the scale estimator                    **/

   c1 = D6;
   c2 = D9;

   if (xmadm <= .0001)
   {
         *xlbiwt=xm;
         *xlbiwt1=xm;
         *xsbiwt=xmadm;
         *xsbiwt1=xmadm;
         return;
   }

   for(i = 0; i < n; i++)
   {
         u1[i] = (xwork[i] - xm)/(c1*xmadm);
         u2[i] = (xwork[i] - xm)/(c2*xmadm);
   }

   s1 = 0.;
   s2 = 0.;
   s3 = 0.;
   s4 = 0.;

   for(i = 0; i < n; i++)
   {
        if (fabs(u2[i]) < D1)
        {
               s1 += (pow((xwork[i]-xm),2.)*pow((D1-(u2[i]*u2[i])),4.));
               s2 += ((D1-u2[i]*u2[i])*(D1-(D5*u2[i]*u2[i])));
        }
        if (fabs(u1[i]) < D1)
        {
               s3 += (xwork[i]-xm)*pow((D1-(u1[i]*u1[i])),2.);
               s4 += pow((D1-(u1[i]*u1[i])),2.);
        }
   }

/** here are the one-step estimators **/

   *xlbiwt1 = xm+s3/s4;
   *xsbiwt1 = ((float) n /sqrt((float) (n-1))) * sqrt(s1) / fabs(s2);

/** now obtain the fully-iterated versions **/

/** solve for new estimates of u1 and u2 **/

   xlb[0] = *xlbiwt1;
   xsb[0] = *xsbiwt1;

   for(i = 1; i <= 10; i++)
   {
         xmm = xlb[i-1];
         for(j = 0; j < n; j++)
         {
                 u1[j] = (xwork[j] - xmm)/(c1*xmadm);
                 u2[j] = (xwork[j] - xmm)/(c2*xmadm);
         }

         s1 = 0.;
         s2 = 0.;
         s3 = 0.;
         s4 = 0.;

         for(j = 0; j < n; j++)
         {
                if (fabs(u2[j]) < D1)
                {
                    s1 += (pow((xwork[j]-xmm),2.)*pow((D1-(u2[j]*u2[j])),4.));
                    s2 += ((D1-u2[j]*u2[j])*(D1-(D5*u2[j]*u2[j])));
                }
                if (fabs(u1[j]) < D1)
                {
                    s3 += (xwork[j]-xmm)*pow((D1-(u1[j]*u1[j])),2.);
                    s4 += pow((D1-(u1[j]*u1[j])),2.);
                }
         }

         xlb[i] = xlb[i-1]+s3/s4;
         xsb[i] = ((float) n /sqrt((float) (n-1))) * sqrt(s1) / fabs(s2);
   }

   *xlbiwt = xlb[10];
   *xsbiwt = xsb[10];
   cpl_free(xwork);
   cpl_free(u1);
   cpl_free(u2);
   cpl_free(xlb);
   cpl_free(xsb);
}


/**************************************************************************/
/**  The XMAD subroutine calculates the Median Absolute Deviation from   **/
/**  the sample median. The median, M , is subtracted from each          **/
/**  ORDERED statistic and then the absolute value is taken. This new    **/
/**  set of of statistics is then resorted so that they are ORDERED      **/
/**  statistics. The MAD is then defined to be the median of this        **/
/**  new set of statistics and is returned as XMADM. The MAD can         **/
/**  be defined:                                                         **/
/**                                                                      **/
/**                 XMADM = median{ abs(x(i) - M) }                      **/
/**                                                                      **/
/**  where the x(i) are the values passed in the array XDATA, and        **/
/**  the median, M, is passed in the array XLETTER. The set of stats     **/
/**  in the brackets is assumed to be resorted. For more information     **/
/**  see page 408 in UREDA.                                              **/
/**                                                                      **/
/**  This function was taken from the library of ROSTAT functions        **/
/**  provided by T. Beers for his ROSTAT code.                           **/
/**  Adapted to C by Marco Scodeggio, sometimes in 1998.                  **/ 
/**************************************************************************/

float xmad(float xdata[],int n,float xmed)
{
   float *xdata2,madVal;
   int i;

   xdata2 = (float *) cpl_malloc(n*sizeof(float));

   for(i = 0; i < n; i++)
        xdata2[i] = fabs(xdata[i] - xmed);
        
   madVal = median(xdata2,n); 
   cpl_free(xdata2);
   return madVal;
}

/*
 * These two are taken from Numerical Recipes, as far as I can tell *
 */

#define NR_END 1
#define FREE_ARG char*

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;
 
        v=(float *)cpl_malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  nh = nh ;
        cpl_free((FREE_ARG) (v+nl-NR_END));
}


/* end of section taken from the Rostat code */

VimosBool findDip1D(float *data, int size, float *position, int minPoints) 
{
  int       i;
  float    *neg;
  VimosBool status;

  neg = (float *) cpl_malloc(size*sizeof(float));
  for (i=0; i<size; i++) neg[i] = -data[i];
  status = findPeak1D(neg, size, position, minPoints);
  cpl_free(neg);
  return status;
}

VimosBool findUpJump(float *data, int size, float *position, int minPoints) 
{
  int       i;
  int       derSize;
  float    *der;
  VimosBool status;

  derSize = size - 1;
  der = (float *) cpl_malloc(derSize*sizeof(float));
  for (i=0; i<derSize; i++) der[i] = MAX(data[i+1]-data[i],0.);
  status = findPeak1D(der, derSize, position, minPoints);
  cpl_free(der);
  if (status == VM_TRUE) *position += 0.5;
  return status;
}

VimosBool findJump(float *data, int size, float *position, int minPoints)
{
  int       i;
  int       derSize;
  float    *der;
  VimosBool status;

  derSize = size - 1;
  der = (float *) cpl_malloc(derSize*sizeof(float));
  for (i=0; i<derSize; i++) der[i] = fabs(data[i+1]-data[i]);
  status = findPeak1D(der, derSize, position, minPoints);
  cpl_free(der);
  if (status == VM_TRUE) *position += 0.5;
  return status;
}

VimosBool findDownJump(float *data, int size, float *position, int minPoints) 
{
  int       i;
  float    *neg;
  VimosBool status;

  neg = (float *) cpl_malloc(size*sizeof(float));
  for (i=0; i<size; i++) neg[i] = -data[i];
  status = findUpJump(neg, size, position, minPoints);
  cpl_free(neg);
  return status;
}




/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  These are functions for integrating & c. (NR- like)

------------------------------------------------------------------------------
*/

/*
  Trapeze integration
*/

#define FUNC(x) ((*theFunc)(x,coeff))

float trapezeInt(float (*theFunc)(float, float[]), float coeff[], 
                 float a, float b, int n)
{
  float x, tnm, sum, del;
  static float s;
  int it, j;

  if (n == 1)
    {
     return (s = 0.5 * (b-a) * (FUNC(a) + FUNC(b)));
    }
  else
    {
     for (it=1,j=1; j<n-1; j++) it <<= 1;
     tnm = it;
     del = (b - a) / tnm;
     x = a + 0.5 * del;
     for (sum=0.0,j=1; j<=it; j++,x+=del) sum += FUNC(x);
     s = 0.5 * (s + (b - a) * sum / tnm);
     return s;
    }
}
#undef FUNC

/*
  Polinomial interpolation
*/

#define NRANSI

void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
  int i, m, ns=1;
  float den, dif, dift, ho, hp, w;
  float *c, *d;

  dif = fabs(x - xa[1]);

  c = floatVector(1,n);
  d = floatVector(1,n);

  for (i=1; i<=n; i++)
     {
      if ((dift = fabs(x - xa[i])) < dif)
        {
         ns = i;
         dif = dift;
        }
      c[i] = ya[i];
      d[i] = ya[i];
     }
  *y = ya[ns--];

  for (m=1; m<n; m++)
     {
      for (i=1; i<=n-m; i++)
         {
          ho = xa[i] - x;
          hp = xa[i+m] - x;
          w = c[i+1] - d[i];
          if ((den = ho-hp) == 0.0) puts("Error in routine polint");
          den = w / den;
          d[i] = hp * den;
          c[i] = ho * den;
         }
      *y += (*dy = (2*ns < (n-m) ? c[ns+1] : d[ns--]));
     }

  freeFloatVector(d,1,n);
  freeFloatVector(c,1,n);

}

#undef NRANSI

/*
  Romberg integration method
*/

#define EPS 1.0e-6
#define JMAX 30
#define JMAXP (JMAX+1)
#define K 5

float rombergInt(float (*theFunc)(float, float[]), float coeff[], 
                 float a, float b)
{
  float ss, dss;
  float s[JMAXP+1], h[JMAXP+1];
  int j;

  h[1] = 1.0;

  for (j=1; j<=JMAX; j++)
     {
      printf("J: %2d",j);
      /* ALEX:     s[j] = trapezeInt(func, a, b, j); */
      s[j] = trapezeInt(theFunc, coeff, a, b, j);
      if (j >= K)
        {
         polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
         if (fabs(dss) < EPS*fabs(ss)) return ss;
        }
      s[j+1] = s[j];
      h[j+1] = 0.25 * h[j];
     }
  puts("Too many steps in routine rombergInt");
  return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K



/*
 Integrate a spectrum - to be improved!
*/

float integrateSpectrum(VimosFloatArray *tmpSpec, float wlenStep)
{
  int i, nstep;
  float integrFlux;

  nstep = tmpSpec->len;
  integrFlux = 0.0;

  /* simple integration by summing all data */
  for (i=0; i<nstep; i++) integrFlux += (tmpSpec->data[i] * wlenStep);

  return(integrFlux);
}

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  End of functions for integrating & c. (NR- like)

------------------------------------------------------------------------------
*/



/*
 Construct an index table to sort an array.
*/

void Indexx(int n, float *arrin, int *indx)
 {
  int l,j,ir,indxt,i;
  float q;

  for (j=0;j<n;j++) indx[j]=j;

  /* this is like n / (2^1) + 1 i.e. (l = half n +1) */
  l=(n >> 1) + 1;

  ir=n;
  for (;;)
     {
      if (l > 1)
        {/* ADDED */
	 /* ALEX: */
	 --l; /* ADDED */
	 q=arrin[(indxt=indx[l-1])];
	}/* ADDED */

      /*ORIGINAL:  q=arrin[(indxt=indx[--l])];*/

      else
        {
	 /* ALEX: */
         q=arrin[(indxt=indx[ir-1])];
         indx[ir-1]=indx[0];

	 /*ORIGINAL:  q=arrin[(indxt=indx[ir])];*/
	 /*ORIGINAL:  indx[ir]=indx[1];*/

         if (--ir == 1)
           {
            indx[0]=indxt;

	    /*ORIGINAL:  indx[1]=indxt;*/

            return;
	   }
	}
      i=l;
      j=l << 1;
      while (j <= ir)
        {
         if (j < ir && arrin[indx[j-1]] < arrin[indx[j]]) j++;

	 /*ORIGINAL:   if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;*/

	 /*ORIGINAL:   if (q < arrin[indx[j]]) {*/
	 /*ORIGINAL:   indx[i]=indx[j];*/

         if (q < arrin[indx[j-1]])
           {
	    indx[i-1]=indx[j-1];
	    j += (i=j);
	   }
	 else j=ir+1;
	}

      /*ORIGINAL: indx[i]=indxt;*/

      indx[i-1]=indxt;
     }
 }


/* 
 Compute RMS of a data set (corrected, not-distorted estimate) 
*/

float computeRMS(float *theData, int n)
{
  int i;
  float med = 0.0;
  float sig = 0.0;
  float variance = 0.0;
  float dif, quaddif;

  /* compute mean */
  for (i=0; i<n; i++)
    {
     med += theData[i];
    }

  med /= (float)n;

  /* compute variance */
  for (i=0; i<n; i++)
     {
      dif=(theData[i] - med);
      quaddif = pow(dif,2);
      variance += quaddif;
     }

  /* corrected, not-distorted estimate: divide by (n-1) instead of by n */
  variance /= (float)(n-1);

  /* compute r.m.s. */
  sig = sqrt(variance);

  return((float)sig);

}




/*
 Compute a 2D array containing the histogram of a data set
*/

int computeHistogram(VimosFloatArray *theData, 
		     int theXaxLen, 
		     VimosFloatArray *histo, float xmin,
		     float xmax, float xstep)
{
  int maxbin, bin;
  int i, k;

  /* largest integer value */
  maxbin = (int)((xmax - xmin) / xstep + 1.);
  if(maxbin != theXaxLen) puts("ERROR!!!");

  for (i=0; i<maxbin; i++)
    histo->data[i]=0.;

  for (k=0; k<theData->len; k++)
    {
      bin=(int)(((theData->data[k] - xmin) / xstep));
      if(bin < 0) bin = 0;
      if(bin >= maxbin) bin = maxbin-1;

      histo->data[bin] += 1.;
    }

  return VM_TRUE;

}
