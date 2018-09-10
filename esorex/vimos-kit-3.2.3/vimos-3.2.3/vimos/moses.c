/* $Id: moses.c,v 1.48 2013-08-22 16:57:33 cgarcia Exp $
 *
 * This file is part of the MOSES library
 * Copyright (C) 2002-2010 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-08-22 16:57:33 $
 * $Revision: 1.48 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include <fors_tools.h>
#include <moses.h>

/* Prototypes */
static cpl_polynomial *read_global_distortion(cpl_table *global, cpl_size row);


/* Cheating: moving here the cpl_tools_get_median_float() prototype,
 * even if cpl_tool.h is not public. It should be removed as soon as 
 * an image median filtering with generic kernel will be implemented
 * in the CPL, or as soon as this module will be moved into the CPL. */

float cpl_tools_get_median_float(float *, cpl_size);

#define MAX_COLNAME      (80)
#define STRETCH_FACTOR   (1.20)

// Related to mos_identify_peaks(), used in multiplex mode

static int mos_multiplex   = -1;
static int mos_region_size = 800;

static double default_lines_hi[] = {   /* Default sky line catalog */
                    5577.338,          /* for high res data        */
                    5889.953,
                    5895.923,
                    5915.301,
                    5932.862,
                    5953.420,
                    6257.961,
                    6287.434,
                    6300.304,
                    6306.869,
                    6363.780,
                    6498.729,
                    6533.044,
                    6553.617,
                    6841.945,
                    6863.955,
                    6870.994,
                    6889.288,
                    6900.833,
                    6912.623,
                    6923.220,
                    6939.521,
                    6969.930,
                    7003.858,
                    7244.907,
                    7276.405,
                    7284.439,
                    7316.282,
                    7329.148,
                    7340.885,
                    7358.659,
                    7571.746,
                    7750.640,
                    7759.996,
                    7794.112,
                    7808.467,
                    7821.503,
                    7841.266,
                    7913.708,
                    7949.204,
                    7964.650,
                    7993.332,
                    8014.059,
                    8310.719,
                    8344.602,
                    8382.392,
                    8399.170,
                    8415.231,
                    8430.174,
                    8452.250,
                    8493.389,
                    8791.186,
                    8827.096,
                    8885.850,
                    8903.114,
                    8943.395,
                    8988.366
                    };

static double default_lines_lo[] = {   /* Default sky line catalog */
                    5577.338,          /* for low res data         */
                    6300.304,
                    6863.955,
                    7571.746,
                    7964.650,
                    7993.332
                    };


/**
 * @defgroup moses MOS data reduction library
 *
 *   The module moses collects low/medium level functions related to 
 *   MOS data reduction.
 */

/**@{*/

/*
 * The following macros and function for finding the k-th smallest
 * value on a float array will be accessible from cpl_tools once
 * this module will be moved into the CPL.
 */

/****

#define PIX_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

static float kthSmallest(float a[], int n, int k)
{
  register int i,j,l,m;
  register float x;

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

#define medianWirth(a, n) kthSmallest(a, n, (((n)&1) ? ((n)/2) : (((n)/2)-1)))

****/

/* 
 * Return random number with gaussian distribution (mean = 0, variance = 1)
 * (Box-Mueller method). The mos_randg() argument is either true or false, 
 * indicating whether to "seed" or not the sequence of generated random 
 * numbers. The "seeding" is performed just at the first mos_randg(1) call, 
 * and at further calls the input argument is ignored. This function
 * generates two random numbers at each call, returning the first one
 * at odd calls, and the second one at even calls.
 */

static void mos_seed(void)
{
    srand((unsigned int)time((time_t *)0));
}

static double mos_randg(int seme)
{
    static int doit = 1;
    static int gotit = 1;
    double x1, x2, w, y1;
    static double y2;

    if (gotit && seme) {
        mos_seed();
        gotit = 0;
    }

    if (doit) {
        doit = 0;
        do {
            x1 = 2.0 * (double)rand() / RAND_MAX - 1.0;
            x2 = 2.0 * (double)rand() / RAND_MAX - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0 || w == 0.0);
    
        w = sqrt( (-2.0 * log(w)) / w);
    
        y1 = x1 * w;
        y2 = x2 * w;
        return y1;
    }

    doit = 1;
    return y2;
}

/* 
 * This function contained a dependency on the VIMOS library
 * (function medianPixelvalue()): it should be removed as soon as an 
 * image median filtering with generic kernel will be implemented
 * in the CPL. Currently it has been solved by a direct call to
 * a cpl_tool function.
 */

static cpl_image *mos_image_vertical_median_filter(cpl_image *ima_in,
                                            int filtsizey, int refrow,
                                            int above, int below, int step)
{

  const char *func = "mos_image_general_median_filter";

  cpl_image  *filt_img = NULL;
  int         col, row;
  float      *buf = NULL;
  float      *data;
  float      *fdata;
  int         upright_y, loleft_y;
  int         j;
  int         yIsEven = !(filtsizey - (filtsizey/2)*2);
  int         f2y;
  int         nx = cpl_image_get_size_x(ima_in);
  int         ny = cpl_image_get_size_y(ima_in);
  int         firstRow;


  if (yIsEven) filtsizey++;

  if (ny <= filtsizey) {
    cpl_msg_error(func, 
                  "Median filter size: %d, image size: %d", filtsizey, ny);
    return NULL;
  }

  f2y = filtsizey / 2;

  filt_img = cpl_image_duplicate(ima_in);
  buf = cpl_malloc(filtsizey * sizeof(float));
  data = cpl_image_get_data(ima_in);
  fdata = cpl_image_get_data(filt_img);

  firstRow = refrow - step * (below / step);
  if (firstRow < f2y)
    firstRow += step;

  for (col = 0; col < nx; col++) {
    for (row = firstRow; row < refrow + above; row += step) {
      if (row >= ny - f2y)
        break;
      loleft_y = row - f2y;
      upright_y = row + f2y + 1;
      for (j = loleft_y; j < upright_y; j++)
        buf[j - loleft_y] = data[col + j * nx];

      fdata[col + row * nx] = cpl_tools_get_median_float(buf, filtsizey);
    }
  }

  cpl_free(buf);

  return filt_img;

}


/*
 * The following static function is used to find an accurate position
 * of a peak within a short interval (however at least 5 pixels long). 
 * The basic idea is to find the baricenter of all the pixel values 
 * that pass a threshold level between the median value and the maximum
 * value within the examined interval (in case such levels are equal,
 * the input is considered flat and no position is returned). At least
 * minPoints must pass this threshold, or no position is computed. To
 * evaluate the significance of the computed baricenter, the variance 
 * of the contributing positions (relative to the found baricenter) is 
 * also evaluated, and compared with the expected variance for a uniform 
 * distribution of positions. If the observed variance is greater than 
 * 80% of the variance of the uniform distribution, the found position 
 * is rejected.
 */

static int peakPosition(const float *data, int size, float *position,
                        int minPoints)
{
  int    i;
  int    count = 0;
  float *copy;
  float  max, median, level, pos, variance, uniformVariance;
  double sum, weights;


  if (data == NULL)
      return 1;

  if (size < 5)         /* Hardcoded, I know... */
      return 1;


  /*
   *  Find median level
   */

  copy = (float *) cpl_malloc(size*sizeof(float));
  for (i = 0; i < size; i++)
      copy[i] = data[i];
  median = cpl_tools_get_median_float(copy, size);
  cpl_free(copy);


  /*
   *  Find max
   */

  max = data[0];
  for (i = 1; i < size; i++)
      if (data[i] > max)
          max = data[i];


  /*
   *  If the max equals the median we have a flat input, therefore
   *  no peak is found.
   */

  if (max-median < 0.00001)
      return 1;


  /*
   *  Discrimination level: only pixels with values above this
   *  level are considered in baricenter calculation.
   */

  level = (max + median) / 2;


  /*
   *  Of the values above this level compute the baricenter and
   *  then the variance of the positions used. Note that the weights
   *  are taken as the difference between the pixels values and
   *  the median level (supposedly the background).
   */

  count = 0;
  for (i = 0, sum = 0., weights = 0.; i < size; i++) {
      if (data[i] > level) {
          count++;
          weights += (data[i] - median);
          sum     += i * (data[i] - median);
      }
  }


  /*
   *  If too few values are above threshold, refuse the position
   *  as insignificant
   */

  if (count < minPoints)
      return 1;

  pos = sum / weights;
  for (i = 0, sum = 0., weights = 0.; i < size; i++) {
      if (data[i] > level) {
          weights++;
          sum += (i - pos) * (i - pos);
      }
  }
  variance = sqrt(sum / weights);


 /*
  *  The "uniform variance" is the variance that should be obtained
  *  in the case of uniform distribution of the points positions in
  *  the selected interval. If the real variance is comparable with
  *  this value, the peak is considered not found.
  */

  uniformVariance = sqrt(size*size/3 - pos*size + pos*pos);

  if (variance > 0.8 * uniformVariance)
      return 1;

  *position = pos + 0.5;

  return 0;
}


/*
 *  The following static function determines the quantity dx to be
 *  added to the position of the highest pixel of a fiber profile,
 *  to get the true position of the profile maximum. All is needed
 *  is the maximum observed value v2 in the profile, and the observed
 *  values v1 and v3 of the previous and the next pixels in the profile.
 *  
 *  The following ratio is defined:
 *  
 *      R = 0.5 (v3 - v1) / (2*v2 - v3 - v1)
 *      
 *  This is a conventional ratio that wouldn't diverge for any set of
 *  pixel values, and that would not depend on the presence of background
 *  (with the assumption that the background level is the same for the 
 *  three pixels). R has also been chosen in such a way that its value
 *  is already quite close to the real dx. It should be noted that the
 *  following condition should be fulfilled:
 *
 *           v1  <= v2   and   v3  <  v2
 *  or
 *           v1  <  v2   and   v3  <=  v2
 *
 *  This implies that dx varies between -0.5 and 0.5 pixels. In such
 *  boundary cases, one has:
 *
 *           v2 = v1   and   R = dx = -0.5
 *           v2 = v3   and   R = dx =  0.5
 *
 *  Another special case is when the observed pixel values are perfectly
 *  symmetrical:
 *
 *           v1 = v3   and   R = dx =  0.0
 *
 *  In all the intermediate cases the relation between R and dx depends
 *  on the shape of the fiber profile, that has been determined elsewhere.
 *  Using the accurate reconstruction of the fiber profile obtained by 
 *  the *  functions ifuProfile() and rebinProfile(), it can be shown 
 *  that R differs from dx always less than 0.01 pixels. If the condition
 *
 *           v1  <= v2   and   v3  <  v2
 *  or
 *           v1  <  v2   and   v3  <=  v2
 *
 *  is not fulfilled, then this function returns the value 2.0.
 */

static double values_to_dx(double v1, double v2, double v3)
{

  static double epsilon = 0.00000001;
  double        r       = 2.0;


  if (v1 > v2 || v3 > v2)
    return r;

  if (2 * v2 - v1 - v3 < epsilon)
    return r;

  r = 0.5 * (v3 - v1) / (2 * v2 - v3 - v1);

  return r;

}


/*
 * The following static function passes a min filter of given box
 * size on the data buffer. The box size must be a positive odd integer.
 */

static float *min_filter(float *buffer, int length, int size)
{
    float *minf  = cpl_calloc(length, sizeof(float));
    float  min;
    int    start = size / 2;
    int    end   = length - size / 2;
    int    i, j;


    for (i = start; i < end; i++) {
        min = buffer[i-start];
        for (j = i - start + 1; j <= i + start; j++)
            if (min > buffer[j])
                min = buffer[j];
        minf[i] = min;
    }

    for (i = 0; i < start; i++)
        minf[i] = minf[start];

    for (i = end; i < length; i++)
        minf[i] = minf[end-1];

    return minf;
}


/*
 * The following static function passes a max filter of given box
 * size on the data buffer. The box size must be a positive odd integer.
 */
 
static float *max_filter(float *buffer, int length, int size)
{
    float *maxf  = cpl_calloc(length, sizeof(float));
    float  max;
    int    start = size / 2;
    int    end   = length - size / 2;
    int    i, j;


    for (i = start; i < end; i++) {
        max = buffer[i-start];
        for (j = i - start + 1; j <= i + start; j++)
            if (max < buffer[j])
                max = buffer[j];
        maxf[i] = max;
    }

    for (i = 0; i < start; i++)
        maxf[i] = maxf[start];

    for (i = end; i < length; i++)
        maxf[i] = maxf[end-1];

    return maxf;
}


/*
 * The following static function passes a running average of given box
 * size on the data buffer. The box size must be a positive odd integer.
 */
 
static float *smo_filter(float *buffer, int length, int size)
{
    float *smof  = cpl_calloc(length, sizeof(float));
    double sum;
    int    start = size / 2;
    int    end   = length - size / 2;
    int    i, j;


    for (i = start; i < end; i++) {
        sum = 0.0;
        for (j = i - start; j <= i + start; j++)
            sum += buffer[j];
        smof[i] = sum / size;
    }

    for (i = 0; i < start; i++)
        smof[i] = smof[start];

    for (i = end; i < length; i++)
        smof[i] = smof[end-1];

    return smof;
}

/*
 * The following two static functions are used to read and write from the 
 * global distortion table the different model components. Conventionally
 * the table consists of 6 columns and 10 rows. Each row is just ordered 
 * storage for model coefficients, and these functions guarantee that the
 * coefficients are read in and written out correctly, independent on their
 * physical meaning. The first 6 table rows are a description of the IDS
 * coefficients, followed by a row containing only the used reference 
 * wavelength. The remaining 3 are a description of the spectral curvature.
 * The first row is a description of coefficient c0, the second of coefficient
 * c1, etc., of the IDS. The 8th row is a description of coefficient c0,
 * the 9th of coefficient c1, etc., of the spectral curvature. All are
 * bivariate polynomialx on x,y mask coordinates. If the input table
 * to the write routine is NULL, it is allocated and initialised. Also
 * the input polynomial could be NULL, and nothing would be written to 
 * the table. If both pointers are NULL the function is basically a
 * constructor of the global distortion table.
 */

static cpl_polynomial *read_global_distortion(cpl_table *global, cpl_size row)
{
    cpl_polynomial *poly = NULL;
    cpl_size        p[2];
    cpl_size        degree = 2;
    int             null;
    double          coeff;

    char   name[MAX_COLNAME];


    for (p[0] = 0; p[0] <= degree; p[0]++) {
        for (p[1] = 0; p[1] <= degree - p[0]; p[1]++) {
            snprintf(name, MAX_COLNAME, "a%"CPL_SIZE_FORMAT"%"CPL_SIZE_FORMAT"", p[0], p[1]);
            coeff = cpl_table_get_double(global, name, row, &null);
            if (null)
                continue;
            if (poly == NULL)
                poly = cpl_polynomial_new(2);
            cpl_polynomial_set_coeff(poly, p, coeff);
        }
    }

    return poly;
}

static cpl_table *write_global_distortion(cpl_table *global, int row, 
                                          cpl_polynomial *poly)
{
    cpl_table *table;
    cpl_size   p[2];
    cpl_size   degree = 2;
    int        nrow = 10;

    char       name[MAX_COLNAME];


    if (global) {
        table = global;
    }
    else {
        table = cpl_table_new(nrow);
        for (p[0] = 0; p[0] <= degree; p[0]++) {
            for (p[1] = 0; p[1] <= degree - p[0]; p[1]++) {
                snprintf(name, MAX_COLNAME, "a%"CPL_SIZE_FORMAT"%"CPL_SIZE_FORMAT"", p[0], p[1]);
                cpl_table_new_column(table, name, CPL_TYPE_DOUBLE);
            }
        }
    }

    if (poly) {
        for (p[0] = 0; p[0] <= degree; p[0]++) {
            for (p[1] = 0; p[1] <= degree - p[0]; p[1]++) {
                snprintf(name, MAX_COLNAME, "a%"CPL_SIZE_FORMAT"%"CPL_SIZE_FORMAT"", p[0], p[1]);
                cpl_table_set_double(table, name, row, 
                                     cpl_polynomial_get_coeff(poly, p));
            }
        }
    }

    return table;
}


/*
 * The following static function is performing a robust linear fit
 * (drawn from the VIMOS library, and originally from ESO-Eclipse).
 *
 *  ----> y = a + b * x
 *
 * This function return 0 on success.
 */

#define SEGNO(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static int robustLinearFit(cpl_bivector *list, double *a, double *b, 
                           double *abdev)
{
    cpl_vector *vx;
    cpl_vector *vy;
    cpl_vector *va;

    double  aa, bb, bcomp, b1, b2, del, abdevt, f, f1, f2, sigb, temp, d, sum;
    double  sx, sy, sxy, sxx, chisq;
    double *arr;
    double  aa_ls, bb_ls;
    double *x;
    double *y;
    int     np;
    int     iter;
    int     max_iterate = 30;
    int     i;


    np = cpl_bivector_get_size(list);
    vx = cpl_bivector_get_x(list);
    vy = cpl_bivector_get_y(list);
    x = cpl_vector_get_data(vx);
    y = cpl_vector_get_data(vy);

    sx = sy = sxx = sxy = 0.00;
    for (i = 0; i < np; i++) {
        sx  += x[i];
        sy  += y[i];
        sxy += x[i] * y[i];
        sxx += x[i] * x[i];
    }

    del = np * sxx - sx * sx;
    aa_ls = aa = (sxx * sy - sx * sxy) / del;
    bb_ls = bb = (np * sxy - sx * sy) / del;

    chisq = 0.00;
    for (i = 0; i < np; i++) {
        temp = y[i] - (aa+bb*x[i]);
        temp *= temp;
        chisq += temp;
    }

    va = cpl_vector_new(np);
    arr = cpl_vector_get_data(va);
    sigb = sqrt(chisq/del);
    b1 = bb;

    bcomp = b1;
    sum = 0.00;
    for (i = 0; i < np; i++) {
        arr[i] = y[i] - bcomp * x[i];
    }
    aa = cpl_vector_get_median_const(va);
    abdevt = 0.0;
    for (i = 0; i < np; i++) {
        d = y[i] - (bcomp * x[i] + aa);
        abdevt += fabs(d);
        if (y[i] != 0.0) 
            d /= fabs(y[i]);
        if (fabs(d) > 1e-7) 
            sum += (d >= 0.0 ? x[i] : -x[i]);
    }
    f1 = sum;

    b2 = bb + SEGNO(3.0 * sigb, f1);

    bcomp = b2;
    sum = 0.00;
    for (i = 0; i < np; i++) {
        arr[i] = y[i] - bcomp * x[i];
    }
    aa = cpl_vector_get_median_const(va);
    abdevt = 0.0;
    for (i = 0; i < np; i++) {
        d = y[i] - (bcomp * x[i] + aa);
        abdevt += fabs(d);
        if (y[i] != 0.0) 
            d /= fabs(y[i]);
        if (fabs(d) > 1e-7) 
            sum += (d >= 0.0 ? x[i] : -x[i]);
    }
    f2 = sum;

    if (fabs(b2-b1)<1e-7) {
        *a = aa;
        *b = bb;
        *abdev = abdevt / (double)np;
        cpl_vector_delete(va);
        return 0;
    }

    iter = 0;
    while (f1*f2 > 0.0) {
        bb = 2.0*b2-b1;
        b1 = b2;
        f1 = f2;
        b2 = bb;

        bcomp = b2;
        sum = 0.00;
        for (i = 0; i < np; i++) {
            arr[i] = y[i] - bcomp * x[i];
        }
        aa = cpl_vector_get_median_const(va);
        abdevt = 0.0;
        for (i = 0; i < np; i++) {
            d = y[i] - (bcomp * x[i] + aa);
            abdevt += fabs(d);
            if (y[i] != 0.0) 
                d /= fabs(y[i]);
            if (fabs(d) > 1e-7) 
                sum += (d >= 0.0 ? x[i] : -x[i]);
        }
        f2 = sum;
        iter++;
        if (iter >= max_iterate) 
            break;
    }
    if (iter >= max_iterate) {
        *a = aa_ls;
        *b = bb_ls;
        *abdev = -1.0;
        cpl_vector_delete(va);
        return 1;
    }

    sigb = 0.01 * sigb;
    while (fabs(b2-b1) > sigb) {
        bb = 0.5 * (b1 + b2);
        if ((fabs(bb-b1) < 1e-7) || (fabs(bb-b2) < 1e-7)) 
            break;
        bcomp = bb;
        sum = 0.0;
        for (i = 0; i < np; i++) {
            arr[i] = y[i] - bcomp * x[i];
        }
        aa = cpl_vector_get_median_const(va);
        abdevt = 0.0;
        for (i = 0; i < np; i++) {
            d = y[i] - (bcomp * x[i] + aa);
            abdevt += fabs(d);
            if (y[i] != 0.0) 
                d /= fabs(y[i]);
            if (fabs(d) > 1e-7) 
                sum += (d >= 0.0 ? x[i] : -x[i]);
        }
        f = sum;

        if (f*f1 >= 0.0) {
            f1=f;
            b1=bb;
        } 
        else {
            f2=f;
            b2=bb;
        }
    }
    cpl_vector_delete(va);
    *a = aa;
    *b = bb;
    *abdev = abdevt / np;
    return 0;
}
#undef SEGNO

/*      
 * The following static function applies the Hough transform from a table
 * of points to another table of points. Given the points p_i = (x_i,y_i)
 * and p_j = (x_j,y_j), the point (X,Y) with X = (y_i - y_j)/(x_i - x_j)
 * and Y = y_i - X*x_i is computed and added to the output table for each
 * p_i, p_j pair. This means that if the input table has N points, the
 * output table has N*(N-1)/2 points.
 */
    
/* static */
cpl_table *mos_hough_table(cpl_table *table, const char *x, const char *y)
{
    cpl_table *output;
    double    *xdata;
    double    *ydata;
    double    *xodata;
    double    *yodata;
    double     max;
    int        npoints;
    int        opoints;
    int        i, j, k;


    if (!cpl_table_has_valid(table, x))
        return NULL;

    npoints = cpl_table_get_nrow(table);
    opoints = npoints*(npoints-1)/2;

    output = cpl_table_new(opoints);
    cpl_table_new_column(output, "m", CPL_TYPE_DOUBLE);
    cpl_table_new_column(output, "q", CPL_TYPE_DOUBLE);

    xodata = cpl_table_get_data_double(output, "m");
    yodata = cpl_table_get_data_double(output, "q");
    
    cpl_table_cast_column(table, x, "x", CPL_TYPE_DOUBLE);
    cpl_table_cast_column(table, y, "y", CPL_TYPE_DOUBLE);

    max = cpl_table_get_column_max(table, "x");
    cpl_table_fill_invalid_double(table, "x", max + 1.0);
    max += 0.5;

    xdata = cpl_table_get_data_double(table, "x");
    ydata = cpl_table_get_data_double(table, "y");

    k = 0;
    for (i = 0; i < npoints; i++) {
        if (xdata[i] < max) {                   // Element i is valid
            for (j = i+1; j < npoints; j++) {
                if (xdata[j] < max) {           // Element j is valid
                    cpl_table_set_double(output, "m", k, 
                                (ydata[i]-ydata[j])/(xdata[i]-xdata[j]));
                    cpl_table_set_double(output, "q", k, 
                                ydata[i] - xodata[k] * xdata[i]);
                    k++;
                }
            }
        }
    }

    if (k != opoints)
        printf("Assert k = %d, expected %d\n", k, opoints);

    cpl_table_erase_column(table, "x");
    cpl_table_erase_column(table, "y");

    return output;
}


/*
 * The following static function is performing the spectral
 * extraction for the function mos_extract_objects()
 */

static void mos_extraction(cpl_image *sciwin, cpl_image *sci_var_win, 
                           cpl_image *skywin, 
                           cpl_image *extracted, cpl_image *sky, 
                           cpl_image *error, int nobjects, int extraction, 
                           double ron, double conad, int ncomb)
{

  cpl_vector *vprofile;
  cpl_image  *smowin;

  int i, j;
  int specLen;
  int numRows;
  int index;
  int iter;
  int maxIter   = 2;         /* Not less than 2 !!! */
  int smoothBox = 31;        /* Not less than 5 !!! */
  double nsigma = 5.0;

  double sumWeight, sum, sumSky, sumProf, sumVar, variance, weight;
  double *profile;
  double *buffer;
  float  *edata;
  float  *ekdata;
  float  *endata;
  float  *sdata;
  float  *kdata;
  float  *fdata;
  float  *vardata;

  double value;


  specLen = cpl_image_get_size_x(sciwin);
  numRows = cpl_image_get_size_y(sciwin);

  edata = cpl_image_get_data(extracted);
  edata += nobjects * specLen;

  ekdata = cpl_image_get_data(sky);
  ekdata += nobjects * specLen;

  endata = cpl_image_get_data(error);
  endata += nobjects * specLen;

  sdata   = cpl_image_get_data(sciwin);
  kdata   = cpl_image_get_data(skywin);
  if(sci_var_win != NULL)
      vardata = cpl_image_get_data(sci_var_win);

  /*
   * Initial spectrum estimate
      if (sdata[i + j * specLen] > 0.0)
   */

  if (extraction && numRows > 5) {
      smowin = mos_image_filter_median(sciwin, 3, 3);
      fdata = cpl_image_get_data(smowin);
      for (i = 0; i < specLen; i++)
        for (j = 0, edata[i] = 0.0; j < numRows; j++)
            edata[i] += fdata[i + j * specLen];
      cpl_image_delete(smowin);
  }
  else {
      for (i = 0; i < specLen; i++)
        for (j = 0, edata[i] = 0.0; j < numRows; j++)
            edata[i] += sdata[i + j * specLen];
  }

  if (extraction) {

    profile = cpl_calloc(specLen * numRows, sizeof(double));
    buffer  = cpl_calloc(specLen, sizeof(double));

    for (iter = 0; iter < maxIter; iter++) {

      /*
       * Normalised spatial profile
       */

      for (i = 0; i < specLen; i++) {
        for (j = 0; j < numRows; j++) {
          index = i + j * specLen;
/*          if (sdata[index] > 0.0 && edata[i] > 0.00001)     */
          if (fabs(edata[i]) > 0.00001)
            profile[index] = sdata[index] / edata[i];
          else
            profile[index] = 0.0;
        }
      }

      for (j = 0; j < numRows; j++) {

        /*
         * Smooth each row in the dispersion direction, and enforce positivity
         */

        for (i = 0; i < specLen - smoothBox; i++) {
          vprofile = cpl_vector_wrap(smoothBox, profile + i + j*specLen);
          value = cpl_vector_get_median_const(vprofile);
          cpl_vector_unwrap(vprofile);
          if (value < 0)
            value = 0.0;
          buffer[i + smoothBox / 2] = value;
        }

        /*
         * Replace the end portions (i.e., not median filtered) with a mean
         */

        vprofile = cpl_vector_wrap(smoothBox / 2, profile + j*specLen);
        value = cpl_vector_get_mean(vprofile);
        cpl_vector_unwrap(vprofile);

        if (value < 0)
            value = 0.0;

        for (i = 0; i < smoothBox / 2; i++)
          buffer[i] = value;

        vprofile = cpl_vector_wrap(smoothBox / 2, 
                                   profile + specLen - smoothBox/2 + j*specLen);
        value = cpl_vector_get_mean(vprofile);
        cpl_vector_unwrap(vprofile);

        if (value < 0)
            value = 0.0;

        for (i = 0; i < smoothBox / 2; i++)
          buffer[i + specLen - smoothBox / 2] = value;

        for (i = 0; i < specLen; i++)
          profile[i + j * specLen] = buffer[i];

      }

      /*
       * Enforce normalization of spatial profile after smoothing
       */

      for (i = 0; i < specLen; i++) {
        for (j = 0, value = 0.0; j < numRows; j++)
          value += profile[i + j * specLen];
        if (value > 0.00001)
          for (j = 0; j < numRows; j++)
            profile[i + j * specLen] /= value;
        else
          for (j = 0; j < numRows; j++)
            profile[i + j * specLen] = 0.0;
      }


      /*
       * Optimal extraction
       */

      for (i = 0; i < specLen; i++) {
        sum = 0.0;
        sumSky = 0.0;
        sumWeight = 0.0;
        sumProf = 0.0;
        sumVar = 0;
        for (j = 0; j < numRows; j++) {
            index = i + j * specLen;
            /*        
            if (sdata[index] > 0.0) {
             */
            //This is the theoretical estimated variance. In principle, since we
            //have the propagated variance, we could use that one, but I leave
            //this as this is the original algorithm (cgarcia)
            variance = ron*ron + fabs(edata[i] * profile[index] + kdata[index])
                                     / conad;
            variance /= ncomb;  /* If input dataset is sum of ncomb images */
            value = sdata[index] - edata[i] * profile[index];
            if (fabs(value) / sqrt(variance) < nsigma) {
                weight = 1000000 * profile[index] / variance;
                sum += weight * sdata[index];
                sumSky += weight * kdata[index];
                sumWeight += weight * profile[index];
                sumProf += profile[index];
                //This is how we propagated the variance. We assume that the
                //weigth has no error, although in has been computed from the
                //profile and the theoretical variance (which also includes the data)
                if(sci_var_win != NULL)
                    sumVar += weight * weight * vardata[index];
            }
        }

        if (sumWeight > 0.00001) {
          edata[i] = sum / sumWeight;
          ekdata[i] = sumSky / sumWeight;
          if(sci_var_win != NULL)
              endata[i] = sqrt(sumVar / sumWeight / sumWeight); //This is the error, not the variance.
          else              
              endata[i] = 1000 * sqrt(sumProf / sumWeight); //This was the old formula, which is not a real error propagation
        }
        else {
/*
          edata[i] = 0.0;
          ekdata[i] = 0.0;
          endata[i] = 0.0;
*/
          //endata[i] = sqrt(ron*ron + fabs(edata[i] + ekdata[i]) / conad);
        }
      }
    }
    cpl_free(profile);
    cpl_free(buffer);
  }
  else {

    /*
     * Add sky estimation for the simple aperture extraction.
        if (kdata[i + j * specLen] > 0.0)
     */

    for (i = 0; i < specLen; i++)
      for (j = 0, ekdata[i] = 0.0; j < numRows; j++)
          ekdata[i] += kdata[i + j * specLen];

    /*
     * Add error estimation for the simple aperture extraction.
     */
    for (i = 0; i < specLen; i++)
    {
        if(sci_var_win != NULL)
        {
            //We propagate the variance of a simple addition
            for (j = 0, endata[i] = 0.0; j < numRows; j++)
                endata[i] += vardata[i + j * specLen];
            endata[i] = sqrt(endata[i]); //We return the error, not the variance
        }
        else 
            endata[i] = sqrt(ron*ron + fabs(edata[i] + ekdata[i]) / conad);
    }
  }

}


/**
 * @brief
 *   Determine all global distortions models
 *
 * @param slits       Table with slits positions on CCD
 * @param maskslits   Table with slits positions on mask
 * @param ids         IDS coefficients table
 * @param crv         Spectral curvature coefficients table
 * @param reference   Reference wavelength
 *
 * @return The global distortions table
 *
 * The @em ids table refers to the spatially rectified spectra, and
 * therefore it has as many rows as the rectified image. The x coordinate
 * of this image corresponds to the x coordinate on the CCD. From the
 * "position" and "length" columns of the @em slits table it is possible
 * to select from the @em ids table the rows belonging to the same slit.
 * For each slit the median value of each coefficient is computed, and
 * associated to the CCD coordinates of the slit of the @em slits table
 * having the same "slit_id". Each coefficient can then be modeled as a
 * function of the CCD coordinates by fitting a second order bivariate 
 * polynomial. The only exception to this scheme is with the modeling
 * of the IDS constant term (the zeropoint of the wavelength calibration),
 * that is expressed as a function of the mask coordinates listed in the
 * @em maskslits table.
 *
 * The @em crv table lists the values of the curvature polynomial for
 * each end of the detected slits. Such values are associated to the 
 * mask coordinates listed in the @em maskslits table having the same
 * "slit_id". Each coefficient can then be modeled as a function of the 
 * mask coordinates by fitting a second order bivariate polynomial.
 *
 * All the coefficients of the obtained bivariate polynomial are written
 * to the newly created global distortions table. Conventionally this 
 * table consists of 6 columns and 10 rows. Each row corresponds to the
 * modeling of one coefficient of the original polynomial coefficients
 * belonging to the local solutions. The first 6 table rows are a 
 * description of the IDS coefficients, up to the fifth polynomial
 * degree; these rows are followed by a row where just the first
 * element is assigned the value of the reference wavelength for 
 * the given IDS model. The remaining 3 rows are a description of 
 * the spectral curvature, up to the second polynomial degree.
 *
 * At least 12 valid slits must be listed in the slits tables.
 */

cpl_table *mos_global_distortion(cpl_table *slits, cpl_table *maskslits,
                                 cpl_table *ids, cpl_table *crv, 
                                 double reference)
{
    const char *func = "mos_global_distortion";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};

    cpl_table      *global = NULL;
    cpl_table      *coeff;
    cpl_table      *dummy;
    cpl_vector     *ci;
    cpl_vector     *xmask;
    cpl_vector     *ymask;
    cpl_bivector   *mask;
    cpl_vector     *xccd;
    cpl_vector     *yccd;
    cpl_bivector   *ccd;
    cpl_polynomial *poly;
    double         *xtop;
    double         *ytop;
    double         *xbottom;
    double         *ybottom;
    double         *mxtop;
    double         *mytop;
    double         *mxbottom;
    double         *mybottom;
    int            *position;
    int            *length;
    int            *slit_id;
    int            *mslit_id;
    int             nslits, nmaskslits, npoints;
    int             order;
    int             i, j;
    int             minslit = 6;    // 12;


    if (slits == NULL || maskslits == NULL || ids == NULL || crv == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nslits = cpl_table_get_nrow(slits);

    if (nslits < minslit) {
        cpl_msg_warning(func, "Too few slits (%d < %d) for global "
                        "distortion model determination", nslits, minslit);
        return NULL;
    }

    nmaskslits = cpl_table_get_nrow(maskslits);

    length   = cpl_table_get_data_int(slits, "length");
    position = cpl_table_get_data_int(slits, "position");
    slit_id  = cpl_table_get_data_int(slits, "slit_id");
    mslit_id = cpl_table_get_data_int(maskslits, "slit_id");
    xtop     = cpl_table_get_data_double(slits, "xtop");
    ytop     = cpl_table_get_data_double(slits, "ytop");
    xbottom  = cpl_table_get_data_double(slits, "xbottom");
    ybottom  = cpl_table_get_data_double(slits, "ybottom");
    mxtop    = cpl_table_get_data_double(maskslits, "xtop");
    mytop    = cpl_table_get_data_double(maskslits, "ytop");
    mxbottom = cpl_table_get_data_double(maskslits, "xbottom");
    mybottom = cpl_table_get_data_double(maskslits, "ybottom");


    /*
     * Global IDS
     */

    coeff = cpl_table_new(nslits);
    cpl_table_copy_structure(coeff, ids);
    cpl_table_new_column(coeff, "xccd", CPL_TYPE_DOUBLE);
    cpl_table_new_column(coeff, "yccd", CPL_TYPE_DOUBLE);
    cpl_table_new_column(coeff, "xmask", CPL_TYPE_DOUBLE);
    cpl_table_new_column(coeff, "ymask", CPL_TYPE_DOUBLE);

    for (i = 0; i < nslits; i++) {
        for (j = 0; j < nmaskslits; j++) {
            if (slit_id[i] == mslit_id[j]) {
                cpl_table_set_double(coeff, "xmask", i,
                                     (mxtop[j] + mxbottom[j]) / 2);
                cpl_table_set_double(coeff, "ymask", i,
                                     (mytop[j] + mybottom[j]) / 2);
            }
        }
    }

    if (cpl_table_has_invalid(coeff, "xmask")) {
        cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
        cpl_table_delete(coeff);
        return NULL;
    }

    for (i = 0; i < nslits; i++) {
        cpl_table_set_double(coeff, "xccd", i, (xtop[i] + xbottom[i]) / 2);
        cpl_table_set_double(coeff, "yccd", i, (ytop[i] + ybottom[i]) / 2);
    }

    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        cpl_table_and_selected_window(ids, position[i], length[i]);
        dummy = cpl_table_extract_selected(ids);
        for (j = 0; j < 6; j++) {
            if (cpl_table_has_column(dummy, clab[j])) {
                if (length[i] - cpl_table_count_invalid(dummy, clab[j]) > 10) {
                    cpl_table_set_double(coeff, clab[j], i, 
                         cpl_table_get_column_median(dummy, clab[j]));
                }
            }
        }

        cpl_table_delete(dummy);
        cpl_table_select_all(ids);
            
    }

    for (j = 0; j < 6; j++) {
        if (cpl_table_has_column(coeff, clab[j])) {
            cpl_table_and_selected_invalid(coeff, clab[j]);

            if (cpl_table_not_selected(coeff))
                dummy = cpl_table_extract_selected(coeff);
            else
                break;

            npoints = cpl_table_get_nrow(dummy);

            if (npoints >= 6) {

                if (npoints >= 12)
                    order = 2;
                else
                    order = 1;
                   
                ci = cpl_vector_wrap(npoints,
                                     cpl_table_get_data_double(dummy, clab[j]));
                if (j) {
                    xccd = cpl_vector_wrap(npoints,
                                     cpl_table_get_data_double(dummy, "xccd"));
                    yccd = cpl_vector_wrap(npoints,
                                     cpl_table_get_data_double(dummy, "yccd"));
                    ccd = cpl_bivector_wrap_vectors(xccd, yccd);

/* %%% */
                    poly = cpl_polynomial_fit_2d_create(ccd, ci, order, NULL);

                    cpl_bivector_unwrap_vectors(ccd);
                    cpl_vector_unwrap(xccd);
                    cpl_vector_unwrap(yccd);
                    cpl_vector_unwrap(ci);
                }
                else {
                    xmask = cpl_vector_wrap(npoints,
                                     cpl_table_get_data_double(dummy, "xmask"));
                    ymask = cpl_vector_wrap(npoints,
                                     cpl_table_get_data_double(dummy, "ymask"));
                    mask = cpl_bivector_wrap_vectors(xmask, ymask);

/* %%% */
                    poly = cpl_polynomial_fit_2d_create(mask, ci, order, NULL);

                    cpl_bivector_unwrap_vectors(mask);
                    cpl_vector_unwrap(xmask);
                    cpl_vector_unwrap(ymask);
                    cpl_vector_unwrap(ci);
                }
            }
            else {
                cpl_size p[2] = {0, 0};
                poly = cpl_polynomial_new(2);
                cpl_polynomial_set_coeff(poly, p, 
                               cpl_table_get_column_median(dummy, clab[j]));
            }

            cpl_table_delete(dummy);

            global = write_global_distortion(global, j, poly);

            cpl_polynomial_delete(poly);

            cpl_table_select_all(coeff);
        }
    }

    cpl_table_delete(coeff);


    /*
     * Add model's reference wavelength
     */

    cpl_table_set_double(global, "a00", 6, reference);


    /*
     * Global curvature model
     */

    coeff = cpl_table_duplicate(crv);
    cpl_table_new_column(coeff, "xmask", CPL_TYPE_DOUBLE);
    cpl_table_new_column(coeff, "ymask", CPL_TYPE_DOUBLE);
    slit_id = cpl_table_get_data_int(coeff, "slit_id");
    npoints = cpl_table_get_nrow(coeff);

    for (i = 0; i < npoints; i++) {
        for (j = 0; j < nmaskslits; j++) {
            if (slit_id[i] == mslit_id[j]) {
                if (i%2) {
                    cpl_table_set_double(coeff, "xmask", i, mxbottom[j]);
                    cpl_table_set_double(coeff, "ymask", i, mybottom[j]);
                }
                else {
                    cpl_table_set_double(coeff, "xmask", i, mxtop[j]);
                    cpl_table_set_double(coeff, "ymask", i, mytop[j]);
                }
            }
        }
    }

    if (cpl_table_has_invalid(coeff, "xmask")) {
        cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
        cpl_table_delete(coeff);
        return NULL;
    }

    for (j = 0; j < 3; j++) {
        if (cpl_table_has_column(coeff, clab[j])) {
            cpl_table_and_selected_invalid(coeff, clab[j]);

            if (cpl_table_not_selected(coeff))
                dummy = cpl_table_extract_selected(coeff);
            else
                break;

            npoints = cpl_table_get_nrow(dummy);

            if (npoints >= 6) {

                if (npoints >= 12)
                    order = 2;
                else
                    order = 1;

                ci = cpl_vector_wrap(npoints,
                                     cpl_table_get_data_double(dummy, clab[j]));
                xmask = cpl_vector_wrap(npoints,
                                     cpl_table_get_data_double(dummy, "xmask"));
                ymask = cpl_vector_wrap(npoints,
                                     cpl_table_get_data_double(dummy, "ymask"));
                mask = cpl_bivector_wrap_vectors(xmask, ymask);

                poly = cpl_polynomial_fit_2d_create(mask, ci, order, NULL);

                cpl_bivector_unwrap_vectors(mask);
                cpl_vector_unwrap(ci);
                cpl_vector_unwrap(xmask);
                cpl_vector_unwrap(ymask);
            }
            else {
                cpl_size p[2] = {0, 0};
                poly = cpl_polynomial_new(2);
                cpl_polynomial_set_coeff(poly, p,
                               cpl_table_get_column_median(dummy, clab[j]));
            }

            cpl_table_delete(dummy);

            global = write_global_distortion(global, j + 7, poly);

            cpl_polynomial_delete(poly);
            cpl_table_select_all(coeff);
        }
    }

    cpl_table_delete(coeff);

    return global;

}

/**
 * @brief
 *   Average global distortion tables.
 *
 * @param global      Array of global distortion tables to average.
 * @param nglobal     Number of global distortion tables to average.
 * @param scale       Expected scale.
 * @param tolerance   Tolerance around expected scale.
 *
 * @return The averaged global distortions table
 *
 * The @em scale is the value of the first element of column a10 of
 * a global distortion table. If @em scale and @em tolerance are 
 * positive, only the tables with a scale within tolerance are
 * averaged.
 */

cpl_table *mos_average_global_distortion(cpl_table **global, int nglobal, 
                                         double scale, double tolerance)
{
    cpl_table *table  = NULL;
    cpl_array *column = NULL;
    int       *good   = NULL;
    int        i, j, ngood, ncolumn, first;
    double     value;


    if (nglobal <= 0) {
        return NULL;
    }

    good = cpl_calloc(nglobal, sizeof(int));
    ngood = 0;

    if (scale > 0.0 && tolerance > 0.0) {
        for (i = 0; i < nglobal; i++) {
            value = cpl_table_get_double(global[i], "a20", 0, NULL);
            if (value != 0.0) {
                value = cpl_table_get_double(global[i], "a10", 0, NULL);
                if (fabs(scale - value) < tolerance) {
                    good[i] = 1;
                    ngood++;
                }
            }
        }
    }
    else {
        for (i = 0; i < nglobal; i++) {
            value = cpl_table_get_double(global[i], "a20", 0, NULL);
            if (value != 0.0) {
                good[i] = 1;
            }
        }
    }

    if (ngood == 0)
        return NULL;

    first = 1;
    for (i = 0; i < nglobal; i++) {
        if (good[i]) {
            if (first) {
                first = 0;
                table = cpl_table_duplicate(global[i]);
                column = cpl_table_get_column_names(table);
                ncolumn = cpl_array_get_size(column);
            }
            else {
                const char *name;
                for (j = 0; j < ncolumn; j++) {
                    name = cpl_array_get_string(column, j);
                    cpl_table_duplicate_column(table, "tmp", global[i], name);
                    cpl_table_add_columns(table, name, "tmp");
                    cpl_table_erase_column(table, "tmp");
                }
            }
        }
    }

    cpl_free(good);

    if (ngood > 1) {
        for (j = 0; j < ncolumn; j++) {
            cpl_table_divide_scalar(table, 
                                    cpl_array_get_string(column, j), ngood);
        }
    }

    cpl_array_delete(column);

    return table;
}


/**
 * @brief
 *   Build the slit location table from a global distortions table
 *
 * @param global      Global distortions table
 * @param maskslits   Table with slits positions on mask
 * @param ysize       Y size of the CCD
 *
 * @return Slits location table
 *
 * The output slits location table has the same structure of the
 * output of the function @c mos_identify_slits(), i.e. the "position"
 * and "length" columns are still missing (such columns would be
 * added by the @c mos_spatial_calibration() function). The column
 * "slit_id" is obtained from the "slit_id" column of the input
 * @em maskslits table, while the "xtop", "ytop", "xbottom", "ybottom"
 * columns are obtained as
 * @code
 *               xtop    = poly0(mxtop, mytop)
 *               xbottom = poly0(mxbottom, mybottom)
 *               ytop    = poly7(mxtop, mytop)
 *                       + poly8(mxtop, mytop) * xtop;
 *                       + poly9(mxtop, mytop) * xtop^2
 *               ybottom = poly7(mxbottom, mybottom)
 *                       + poly8(mxbottom, mybottom) * xbottom;
 *                       + poly9(mxbottom, mybottom) * xbottom^2
 * @endcode
 * where polyX is the polynomial obtained from row X
 * of the input global distortions table, and mxtop, mytop, mxbottom,
 * mybottom are the coordinates of the slits ends listed in the input
 * @em maskslits table. The resulting slits location table is finally
 * sorted according to the descending value of "ytop", ordering the
 * slits from top to bottom of the image starting from the first table 
 * row. The slits that are completely outside the CCD are excluded
 * from the table. The argument @em ysize is used for this purpose.
 */

cpl_table *mos_build_slit_location(cpl_table *global, cpl_table *maskslits,
                                   int ysize)
{
    const char *func = "mos_build_slit_location";

    cpl_propertylist *sort_col;
    cpl_polynomial   *ids0;
    cpl_polynomial   *crv[3];
    cpl_polynomial   *loc_crv;
    cpl_vector       *point;
    cpl_table        *slits;
    cpl_size         nslits;
    int              *slit_id;
    double           *dpoint;
    double           *xtop;
    double           *ytop;
    double           *xbottom;
    double           *ybottom;
    double           *mxtop;
    double           *mytop;
    double           *mxbottom;
    double           *mybottom;
    cpl_size          i;
    cpl_size          j;


    if (global == NULL || maskslits == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nslits   = cpl_table_get_nrow(maskslits);
    slit_id  = cpl_table_get_data_int(maskslits, "slit_id");
    mxtop    = cpl_table_get_data_double(maskslits, "xtop");
    mytop    = cpl_table_get_data_double(maskslits, "ytop");
    mxbottom = cpl_table_get_data_double(maskslits, "xbottom");
    mybottom = cpl_table_get_data_double(maskslits, "ybottom");

    slits = cpl_table_duplicate(maskslits);

    xtop    = cpl_table_get_data_double(slits, "xtop");
    ytop    = cpl_table_get_data_double(slits, "ytop");
    xbottom = cpl_table_get_data_double(slits, "xbottom");
    ybottom = cpl_table_get_data_double(slits, "ybottom");

    ids0 = read_global_distortion(global, 0);
    crv[0] = read_global_distortion(global, 7);
    crv[1] = read_global_distortion(global, 8);
    crv[2] = read_global_distortion(global, 9);

    loc_crv = cpl_polynomial_new(1);

    point = cpl_vector_new(2);
    dpoint = cpl_vector_get_data(point);

    for (i = 0; i < nslits; i++) {
        dpoint[0] = mxtop[i];
        dpoint[1] = mytop[i];

        xtop[i] = cpl_polynomial_eval(ids0, point);

        for (j = 0; j < 3; j++)
            if (crv[j])
                cpl_polynomial_set_coeff(loc_crv, &j, 
                                         cpl_polynomial_eval(crv[j], point));

        ytop[i] = cpl_polynomial_eval_1d(loc_crv, xtop[i], NULL);

        dpoint[0] = mxbottom[i];
        dpoint[1] = mybottom[i];
        xbottom[i] = cpl_polynomial_eval(ids0, point);

        for (j = 0; j < 3; j++)
            if (crv[j])
                cpl_polynomial_set_coeff(loc_crv, &j,
                                         cpl_polynomial_eval(crv[j], point));

        ybottom[i] = cpl_polynomial_eval_1d(loc_crv, xbottom[i], NULL);
    }

    cpl_vector_delete(point);
    cpl_polynomial_delete(ids0);
    cpl_polynomial_delete(loc_crv);
    for (j = 0; j < 3; j++)
        cpl_polynomial_delete(crv[j]);

    sort_col = cpl_propertylist_new();
    cpl_propertylist_append_bool(sort_col, "ytop", 1);
    cpl_table_sort(slits, sort_col);
    cpl_table_sort(maskslits, sort_col);
    cpl_propertylist_delete(sort_col);

    /*
     * Eliminate slits which are _entirely_ outside the CCD
     */

    cpl_table_and_selected_double(slits, "ybottom", CPL_GREATER_THAN, ysize-1);
    cpl_table_or_selected_double(slits, "ytop", CPL_LESS_THAN, 0);
    cpl_table_erase_selected(slits);

    nslits = cpl_table_get_nrow(slits);

    if (nslits == 0) {
        cpl_msg_warning(func, "No slits found on the CCD");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        cpl_table_delete(slits);
        return NULL;
    }

    if (nslits > 1)
        cpl_msg_info(func, "Slit location: %"CPL_SIZE_FORMAT" slits are entirely or partially "
                     "contained in CCD", nslits);
    else
        cpl_msg_info(func, "Slit location: %"CPL_SIZE_FORMAT" slit is entirely or partially "
                     "contained in CCD", nslits);

    return slits;

}


/**
 * @brief
 *   Build the curvature coefficients table from a global distortions table
 *
 * @param global      Global distortions table
 * @param maskslits   Table with slits positions on mask
 * @param slits       Table with slits positions on CCD
 *
 * @return Curvature coefficients table
 *
 * The output curvature coefficients table has the same structure of the
 * output of the function @c mos_poly_trace(). The column "slit_id" is 
 * obtained from the "slit_id" column of the input @em maskslits table.
 * The coefficients columns are obtained as
 * @code
 *               c0 = poly7(mx, my)
 *               c1 = poly8(mx, my)
 *               c2 = poly9(mx, my)
 * @endcode
 * where polyX is the polynomial obtained from row X
 * of the input global distortions table, and (mx, my) are the 
 * coordinates of the slits ends listed in the input @em maskslits 
 * table. The slits that are completely outside the CCD are excluded
 * from the table.
 */

cpl_table *mos_build_curv_coeff(cpl_table *global, cpl_table *maskslits,
                                cpl_table *slits)
{
    const char *func = "mos_build_curv_coeff";

    const char     *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_polynomial *crv[3];
    cpl_vector     *point;
    cpl_table      *polytraces;
    double         *dpoint;
    double         *xtop;
    double         *ytop;
    double         *xbottom;
    double         *ybottom;
    int            *slit_id;
    int            *valid_id;
    int             nslits, nvalid;
    int             found;
    int             i, j, k;


    if (global == NULL || slits == NULL || maskslits == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nslits  = cpl_table_get_nrow(maskslits);
    slit_id = cpl_table_get_data_int(maskslits, "slit_id");
    xtop    = cpl_table_get_data_double(maskslits, "xtop");
    ytop    = cpl_table_get_data_double(maskslits, "ytop");
    xbottom = cpl_table_get_data_double(maskslits, "xbottom");
    ybottom = cpl_table_get_data_double(maskslits, "ybottom");

    polytraces = cpl_table_new(2*nslits);
    cpl_table_new_column(polytraces, "slit_id", CPL_TYPE_INT);
    for (i = 0; i < 3; i++)
        cpl_table_new_column(polytraces, clab[i], CPL_TYPE_DOUBLE);

    crv[0] = read_global_distortion(global, 7);
    crv[1] = read_global_distortion(global, 8);
    crv[2] = read_global_distortion(global, 9);

    point = cpl_vector_new(2);
    dpoint = cpl_vector_get_data(point);

    for (i = 0; i < nslits; i++) {
        for (j = 0; j < 2; j++) {  /* For top and bottom trace of each slit */

            cpl_table_set_int(polytraces, "slit_id", 2*i+j, slit_id[i]);

            if (j) {
                dpoint[0] = xbottom[i];
                dpoint[1] = ybottom[i];                
            }
            else {
                dpoint[0] = xtop[i];
                dpoint[1] = ytop[i];                
            }

            for (k = 0; k < 3; k++)
                if (crv[j])
                    cpl_table_set_double(polytraces, clab[k], 2*i+j,
                                         cpl_polynomial_eval(crv[k], point));
        }
    }

    cpl_vector_delete(point);
    for (j = 0; j < 3; j++)
        cpl_polynomial_delete(crv[j]);

    /*
     * Eliminate slits which are _entirely_ outside the CCD
     */
 
    nvalid  = cpl_table_get_nrow(slits);
    valid_id = cpl_table_get_data_int(slits, "slit_id");
    cpl_table_unselect_all(polytraces);
    for (i = 0; i < nslits; i++) {
        found = 0;
        for (j = 0; j < nvalid; j++) {
            if (slit_id[i] == valid_id[j]) {
                found = 1;
                break;
            }
        }
        if (!found) {
            cpl_table_select_row(polytraces, 2*i);
            cpl_table_select_row(polytraces, 2*i + 1);
        }
    }
    cpl_table_erase_selected(polytraces);
 
    nslits = cpl_table_get_nrow(polytraces);

    if (nslits == 0) {
        cpl_msg_warning(func, "No slits found on the CCD");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        cpl_table_delete(polytraces);
        return NULL;
    }

    if (nslits > 2) 
        cpl_msg_info(func, "Curvature model: %d slits are entirely or "
                     "partially contained in CCD", nslits / 2);
    else
        cpl_msg_info(func, "Curvature model: %d slit is entirely or "
                     "partially contained in CCD", nslits / 2);

    return polytraces;
}


/**
 * @brief
 *   Build the IDS coefficients table from a global distortions table
 *   
 * @param global      Global distortions table
 * @param slits       Table with slits positions on CCD
 * 
 * @return IDS coefficients table
 * 
 * The input @em slits table should be already processed by the function
 * @c mos_spatial_calibration(), i.e., it should already have the columns
 * "position" and "length". The output IDS coefficients table will have the 
 * same structure of the one created by @c mos_wavelength_calibration_final(), 
 * but without the "error" and "nlines" columns. This output table will have
 * as many rows as the sum of the "length" values in the input @em slits 
 * table, corresponding to the number of spatial pseudo-pixels of a 
 * corresponding rectified data image. If a given slit extends between CCD
 * coordinates ytop and ybottom (as listed in the input @em slits table), 
 * the corresponding number of spatial pseudo-pixels is conventionally set
 * to length = ceil(ytop-ybottom)+1 (corresponding to the content of the
 * "length" column, as computed by the @c mos_spatial_calibration() function). 
 * The spatial pseudo-pixels p are counted from top to bottom, starting from 0,
 * and their corresponding y coordinate on the CCD is therefore given by 
 * y = ytop - p*(ytop-ybottom)/length . The corresponding x coordinate 
 * is computed in the same way, as x = xtop - p*(xtop-xbottom)/length .
 * The coefficients columns are obtained as
 * @code
 *               c0 = poly0(xmask, ymask)
 *               c1 = poly1(x, y)
 *               c2 = poly2(x, y)
 *               c3 = poly3(x, y)
 *               c4 = poly4(x, y)
 *               c5 = poly5(x, y)
 * @endcode
 * where polyX is the polynomial obtained from row X
 * of the input global distortions table, and the (xmask, ymask)
 * are the slit coordinate on the mask. The value of poly0(xmask, ymask)
 * is already contained in the input @em slits table (columns "xtop" and
 * "xbottom").
 */

cpl_table *mos_build_disp_coeff(cpl_table *global, cpl_table *slits)
{
    const char *func = "mos_build_disp_coeff";

    const char     *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};

    cpl_polynomial *ids[6];
    cpl_vector     *point;
    cpl_table      *idscoeff;
    double         *dpoint;
    double         *xtop;
    double         *ytop;
    double         *xbottom;
    double         *ybottom;
    int            *position;
    int            *length;
    int             nslits;
    int             nrows;
    int             order;
    int             ylow, yhig;
    int             i, j, k;


    if (global == NULL || slits == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }
    
    nslits   = cpl_table_get_nrow(slits);
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");
    xtop     = cpl_table_get_data_double(slits, "xtop");
    ytop     = cpl_table_get_data_double(slits, "ytop");
    xbottom  = cpl_table_get_data_double(slits, "xbottom");
    ybottom  = cpl_table_get_data_double(slits, "ybottom");

    for (i = 0; i < 6; i++)
        ids[i] = read_global_distortion(global, i);

    for (i = 0; i < 6; i++)
        if (ids[i] == NULL)
            break;

    order = i - 1;

    if (order < 1) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    nrows = 0;
    for (i = 0; i < nslits; i++)
        nrows += length[i]; 

    idscoeff = cpl_table_new(nrows);

    for (j = 0; j <= order; j++)
        cpl_table_new_column(idscoeff, clab[j], CPL_TYPE_DOUBLE);

    cpl_table_new_column(idscoeff, "error", CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(idscoeff, "error", 0, nrows, 0.0);
    cpl_table_new_column(idscoeff, "nlines", CPL_TYPE_INT);
    cpl_table_fill_column_window_int(idscoeff, "nlines", 0, nrows, 0);

    point = cpl_vector_new(2);
    dpoint = cpl_vector_get_data(point);

    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        ylow = position[i];
        yhig = ylow + length[i];

        for (j = 0; j <= order; j++) {
            if (j) {
                for (k = 0; k < length[i]; k++) {
                    dpoint[0] = xbottom[i] + k*(xtop[i]-xbottom[i])/length[i];
                    dpoint[1] = ybottom[i] + k*(ytop[i]-ybottom[i])/length[i];
                    cpl_table_set_double(idscoeff, clab[j], ylow + k,
                                         cpl_polynomial_eval(ids[j], point));
                }
            }
            else {
                for (k = 0; k < length[i]; k++) {
                    cpl_table_set_double(idscoeff, clab[0], ylow + k,
                                xbottom[i] + k*(xtop[i]-xbottom[i])/length[i]);
                }
            }
        }
    }

    cpl_vector_delete(point);
    for (j = 0; j < 6; j++)
        cpl_polynomial_delete(ids[j]);

    return idscoeff;

}


/**
 * @brief 
 *   Subtract the sky from the scientific CCD exposure
 *  
 * @param science     Image containing the CCD scientific spectra
 * @param slits       Table with slits positions
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param reference   Reference wavelength
 * @param blue        Start lambda to process
 * @param red         End lambda to process
 * @param dispersion  Mean spectral dispersion
 *  
 * @return The sky subtracted scientific image
 *
 * The input @em science frame should be already bias subtracted, and should
 * be oriented so that the dispersion direction is horizontal with @em blue
 * on the left and @em red on the right. The sky level is determined by
 * robust linear fitting of flux along the spatial direction. The input 
 * @em science image is sky subtracted IN PLACE using the obtained sky 
 * map, that is also returned. 
 */
 
cpl_image *mos_subtract_sky(cpl_image *science, cpl_table *slits, 
                            cpl_table *polytraces, double reference, 
                            double blue, double red, double dispersion)
{
    const char     *func = "mos_subtract_sky";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_image      *sky;
    cpl_bivector   *list;
    cpl_vector     *listx;
    cpl_vector     *listy;
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;
    cpl_polynomial *trend;

    int            *slit_id;
    double         *dlistx;
    double         *dlisty;
    float          *sdata;
    float          *kdata;
    double          top, bot;
    int             itop, ibot;
    double          coeff;
    double          ytop, ybot;
    double          m, q, err;
    int             npix;

    int             pixel_above, pixel_below, refpixel, start_pixel, end_pixel;
    int             nx, ny;
    int             nslits;
    int            *length;
    int             missing_top, missing_bot;
    int             order;
    int             null;
    int             window = 50;  /* Longer slits have polynomial sky model */
    int             count;
    int             i, j;
    cpl_size        k;


    if (science == NULL || slits == NULL || polytraces == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }
 
    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (red - blue < dispersion) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    nx = cpl_image_get_size_x(science);
    ny = cpl_image_get_size_y(science);

    sky = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);

    sdata = cpl_image_get_data(science);
    kdata = cpl_image_get_data(sky);

    nslits   = cpl_table_get_nrow(slits);
    order    = cpl_table_get_ncol(polytraces) - 2;
    length   = cpl_table_get_data_int(slits, "length");
    slit_id  = cpl_table_get_data_int(slits, "slit_id");

    /*
     * The spatial resampling is performed for a certain number of
     * pixels above and below the position of the reference wavelength:
     */
    
    pixel_above = STRETCH_FACTOR * (red - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - blue) / dispersion;

    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        
        /*
         * Recover from the table of spectral curvature coefficients
         * the curvature polynomials.
         */

        refpixel = cpl_table_get_double(slits, "xtop", i, NULL);

        start_pixel = refpixel - pixel_below;
        if (start_pixel < 0)
            start_pixel = 0;

        end_pixel = refpixel + pixel_above;
        if (end_pixel > nx)
            end_pixel = nx;

        missing_top = 0;
        polytop = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i, &null);
            if (null) {
                cpl_polynomial_delete(polytop);
                missing_top = 1;
                break;
            }
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        missing_bot = 0;
        polybot = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i+1, &null);
            if (null) {
                cpl_polynomial_delete(polybot);
                missing_bot = 1;
                break;
            }
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        if (missing_top && missing_bot) {
            cpl_msg_debug(func, "Slit %d was not traced: no extraction!",
                          slit_id[i]);
            continue;
        }

        /*
         * In case just one of the two edges was not traced, the other
         * edge curvature model is duplicated and shifted to the other
         * end of the slit: better than nothing!
         */

        if (missing_top) {
            cpl_msg_debug(func, "Upper edge of slit %d was not traced: "
                          "the spectral curvature of the lower edge "
                          "is used instead.", slit_id[i]);
            polytop = cpl_polynomial_duplicate(polybot);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polybot, &k);
            coeff += ytop - ybot;
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        if (missing_bot) {
            cpl_msg_debug(func, "Lower edge of slit %d was not traced: "
                          "the spectral curvature of the upper edge "
                          "is used instead.", slit_id[i]);
            polybot = cpl_polynomial_duplicate(polytop);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polytop, &k);
            coeff -= ytop - ybot;
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }


        /*
         * Now read pixel values along spatial direction, and fit them.
         */

        for (j = start_pixel; j < end_pixel; j++) {
            top = cpl_polynomial_eval_1d(polytop, j, NULL);
            bot = cpl_polynomial_eval_1d(polybot, j, NULL);
            itop = floor(top + 0.5) + 1;
            ibot = floor(bot + 0.5);
            if (itop > ny)
                itop = ny;
            if (ibot < 0)
                ibot = 0;
            npix = itop - ibot;
            if (npix < 5)
                break;

            list = cpl_bivector_new(npix);
            listx = cpl_bivector_get_x(list);
            listy = cpl_bivector_get_y(list);
            dlistx = cpl_vector_get_data(listx);
            dlisty = cpl_vector_get_data(listy);

            for (k = 0; k < npix; k++) {
                dlistx[k] = k;
                dlisty[k] = sdata[j + (ibot + k)*nx];
            }

            if (robustLinearFit(list, &q, &m, &err)) {
                cpl_bivector_delete(list);
                continue;
            }

            cpl_bivector_delete(list);

            for (k = 0; k < npix; k++) {
                kdata[j + (ibot + k)*nx] = m*k + q;
            }

            if (npix > window) {

                /*
                 * Polynomial iteration
                 */

                err = 3*sqrt(err);

                count = 0;
                for (k = 0; k < npix; k++)
                    if (fabs(sdata[j + (ibot + k)*nx] - m*k - q) < err)
                        count++;

                if (count < 10)
                    continue;

                list = cpl_bivector_new(count);
                listx = cpl_bivector_get_x(list);
                listy = cpl_bivector_get_y(list);
                dlistx = cpl_vector_get_data(listx);
                dlisty = cpl_vector_get_data(listy);

                count = 0;
                for (k = 0; k < npix; k++) {
                    if (fabs(sdata[j + (ibot + k)*nx] - m*k - q) < err) {
                        dlistx[count] = k;
                        dlisty[count] = sdata[j + (ibot + k)*nx];
                        count++;
                    }
                }

                trend = cpl_polynomial_fit_1d_create(listx, listy, 2, &err);
 
                cpl_bivector_delete(list);

                err = 3*sqrt(err);

                count = 0;
                for (k = 0; k < npix; k++)
                    if (fabs(sdata[j + (ibot + k)*nx] 
                             - cpl_polynomial_eval_1d(trend, k, NULL)) < err)
                        count++;

                if (count < 10) {
                    cpl_polynomial_delete(trend);
                    continue;
                }

                list = cpl_bivector_new(count);
                listx = cpl_bivector_get_x(list);
                listy = cpl_bivector_get_y(list);
                dlistx = cpl_vector_get_data(listx);
                dlisty = cpl_vector_get_data(listy);

                count = 0;
                for (k = 0; k < npix; k++) {
                    if (fabs(sdata[j + (ibot + k)*nx] 
                             - cpl_polynomial_eval_1d(trend, k, NULL)) < err) {
                        dlistx[count] = k;
                        dlisty[count] = sdata[j + (ibot + k)*nx];
                        count++;
                    }
                }

                cpl_polynomial_delete(trend);

                trend = cpl_polynomial_fit_1d_create(listx, listy, 3, &err);

                cpl_bivector_delete(list);
 
                for (k = 0; k < npix; k++) {
                    kdata[j + (ibot + k)*nx] = cpl_polynomial_eval_1d(trend, 
                                               k, NULL);
                }

                cpl_polynomial_delete(trend);
            }
        }
        cpl_polynomial_delete(polytop);
        cpl_polynomial_delete(polybot);
    }

    cpl_image_subtract(science, sky);

    return sky;
}


/**
 * @brief 
 *   Select wavelengths from a single spectrum.
 *  
 * @param wavemap     Wavelength map
 * @param slit        Slit number
 * @param slits       Table with slits positions
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param reference   Reference wavelength
 * @param blue        Start lambda to process
 * @param red         End lambda to process
 * @param dispersion  Mean spectral dispersion
 *  
 * @return @c CPL_ERROR_NONE on success.
 *
 * The input wavelength map is filled with zeroes, but for the specified
 * slit.
 */
 
cpl_error_code mos_slit_wavemap(cpl_image *wavemap, int slit, cpl_table *slits, 
                                cpl_table *polytraces, double reference, 
                                double blue, double red, double dispersion)
{
    const char     *func = "mos_slit_wavemap";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_image      *mask;
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;

    int            *slit_id;
    float          *sdata;
    double          top, bot;
    int             itop, ibot;
    double          coeff;
    double          ytop, ybot;
    int             npix;

    int             pixel_above, pixel_below, refpixel, start_pixel, end_pixel;
    int             nx, ny;
    int             nslits;
    int            *length;
    int             missing_top, missing_bot;
    int             order;
    int             null;
    int             j;
    cpl_size        k;


    if (wavemap == NULL || slits == NULL || polytraces == NULL)
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);
 
    if (dispersion <= 0.0)
        return cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);

    if (red - blue < dispersion)
        return cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);

    nslits = cpl_table_get_nrow(slits);

    if (slit < 0 || slit >= nslits)
        return cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);

    nx = cpl_image_get_size_x(wavemap);
    ny = cpl_image_get_size_y(wavemap);

    mask = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
    sdata = cpl_image_get_data(mask);

    order    = cpl_table_get_ncol(polytraces) - 2;
    length   = cpl_table_get_data_int(slits, "length");
    slit_id  = cpl_table_get_data_int(slits, "slit_id");

    /*
     * The spatial reading is performed for a certain number of
     * pixels above and below the position of the reference wavelength:
     */
    
    pixel_above = STRETCH_FACTOR * (red - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - blue) / dispersion;

    if (length[slit] == 0)
        return cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        
        
    /*
     * Recover from the table of spectral curvature coefficients
     * the curvature polynomials.
     */

    refpixel = cpl_table_get_double(slits, "xtop", slit, NULL);

    start_pixel = refpixel - pixel_below;
    if (start_pixel < 0)
        start_pixel = 0;

    end_pixel = refpixel + pixel_above;
    if (end_pixel > nx)
        end_pixel = nx;

    missing_top = 0;
    polytop = cpl_polynomial_new(1);
    for (k = 0; k <= order; k++) {
        coeff = cpl_table_get_double(polytraces, clab[k], 2*slit, &null);
        if (null) {
            cpl_polynomial_delete(polytop);
            missing_top = 1;
            break;
        }
        cpl_polynomial_set_coeff(polytop, &k, coeff);
    }

    missing_bot = 0;
    polybot = cpl_polynomial_new(1);
    for (k = 0; k <= order; k++) {
        coeff = cpl_table_get_double(polytraces, clab[k], 2*slit+1, &null);
        if (null) {
            cpl_polynomial_delete(polybot);
            missing_bot = 1;
            break;
        }
        cpl_polynomial_set_coeff(polybot, &k, coeff);
    }

    if (missing_top && missing_bot) {
        cpl_msg_debug(func, "Slit %d was not traced: no extraction!",
                      slit_id[slit]);
        return cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
    }

    /*
     * In case just one of the two edges was not traced, the other
     * edge curvature model is duplicated and shifted to the other
     * end of the slit: better than nothing!
     */

    if (missing_top) {
        cpl_msg_debug(func, "Upper edge of slit %d was not traced: "
                      "the spectral curvature of the lower edge "
                      "is used instead.", slit_id[slit]);
        polytop = cpl_polynomial_duplicate(polybot);
        ytop = cpl_table_get_double(slits, "ytop", slit, NULL);
        ybot = cpl_table_get_double(slits, "ybottom", slit, NULL);
        k = 0;
        coeff = cpl_polynomial_get_coeff(polybot, &k);
        coeff += ytop - ybot;
        cpl_polynomial_set_coeff(polytop, &k, coeff);
    }

    if (missing_bot) {
        cpl_msg_debug(func, "Lower edge of slit %d was not traced: "
                      "the spectral curvature of the upper edge "
                      "is used instead.", slit_id[slit]);
        polybot = cpl_polynomial_duplicate(polytop);
        ytop = cpl_table_get_double(slits, "ytop", slit, NULL);
        ybot = cpl_table_get_double(slits, "ybottom", slit, NULL);
        k = 0;
        coeff = cpl_polynomial_get_coeff(polytop, &k);
        coeff -= ytop - ybot;
        cpl_polynomial_set_coeff(polybot, &k, coeff);
    }


    /*
     * Now read pixel values along spatial direction, and fit them.
     */

    for (j = start_pixel; j < end_pixel; j++) {
        top = cpl_polynomial_eval_1d(polytop, j, NULL);
        bot = cpl_polynomial_eval_1d(polybot, j, NULL);
        itop = floor(top + 0.5) + 1;
        ibot = floor(bot + 0.5);
        if (itop > ny)
            itop = ny;
        if (ibot < 0)
            ibot = 0;
        npix = itop - ibot;

        for (k = 0; k < npix; k++) {
            sdata[j + (ibot + k)*nx] = 1.0;
        }
    }

    cpl_polynomial_delete(polytop);
    cpl_polynomial_delete(polybot);

    cpl_image_multiply(wavemap, mask);

    cpl_image_delete(mask);

    return CPL_ERROR_NONE;
}


/**
 * @brief 
 *   Normalise a flat field exposure
 *  
 * @param flat        Image containing the original flat field spectra
 * @param spatial     Spatial calibration image
 * @param slits       Table with slits positions
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param reference   Reference wavelength
 * @param blue        Start lambda to process
 * @param red         End lambda to process
 * @param dispersion  Mean spectral dispersion
 * @param sradius     Radius of smoothing box along the dispersion direction
 * @param polyorder   Order of fitting polynomial along the dispersion
 *  
 * @return The smoothed flat field exposure used for normalisation
 *
 * The input @em flat frame should be already bias subtracted, and should
 * be oriented so that the dispersion direction is horizontal with @em blue
 * on the left and @em red on the right. The flat field spectra are spatially 
 * rectified, heavily smoothed, and then mapped back on the CCD. The original 
 * @em flat image is divided IN PLACE by its smoothed counterpart, which is 
 * also returned. If the polynomial @em polyorder is set to a negative number 
 * the smoothing consists of a linear fit along the spatial direction 
 * (excluding 3+3 pixels at the spectral edges), and by a median filtering 
 * along the dispersion direction using a window with the specified 
 * @em sradius; alternatively, if @em polyorder is not negative, the smoothing
 * will consist of a polynomial fitting of the illumination profile along
 * the dispersion direction, performed independently for each row of the
 * spatially remapped spectra.
 */
 
cpl_image *mos_normalise_flat(cpl_image *flat, cpl_image *spatial, 
                              cpl_table *slits, cpl_table *polytraces, 
                              double reference, double blue, double red, 
                              double dispersion, int sradius, int polyorder)
{
    const char     *func = "mos_normalise_flat";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_image      *rectified;
    cpl_image      *smo_flat;
    cpl_image      *exslit;
    cpl_vector     *positions;
    cpl_vector     *flux;
    cpl_vector     *smo_flux;
    cpl_polynomial *trend;
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;

    int            *slit_id;
    float          *p;
    float          *data;
    double         *fdata;
    double         *pdata;
    float          *sdata;
    float          *xdata;
    float          *wdata;
    double          vtop, vbot, value;
    double          top, bot;
    double          coeff;
    double          ytop, ybot;
    double          ypos;
    double          fvalue;
    int             ivalue;
    int             yint, yprev;
    int             npseudo;

    int             pixel_above, pixel_below, refpixel, start_pixel, end_pixel;
    int             nx, ny, nsubx, nsuby;
    int             xlow, ylow, xhig, yhig;
    int             nslits;
    int            *position;
    int            *length;
    int             missing_top, missing_bot;
    int             order;
    int             npoints;
    int             uradius;
    int             null;
    int             i, j;
    cpl_size        k;

/*    int             exclude = 5;     Number of excluded pixels at edges */

    /* For debug puposes only: cpl_image      *smo_rectified; */


    if (flat == NULL || slits == NULL || polytraces == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }
 
    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (red - blue < dispersion) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    rectified = mos_spatial_calibration(flat, slits, polytraces, reference,
                                        blue, red, dispersion, 0, NULL);

    nx = cpl_image_get_size_x(rectified);
    ny = cpl_image_get_size_y(rectified);

    smo_flat = cpl_image_new(cpl_image_get_size_x(spatial), 
                             cpl_image_get_size_y(spatial), CPL_TYPE_FLOAT);
    wdata = cpl_image_get_data(smo_flat);

    nslits   = cpl_table_get_nrow(slits);
    order    = cpl_table_get_ncol(polytraces) - 2;
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");
    slit_id  = cpl_table_get_data_int(slits, "slit_id");

    /*
     * The spatial resampling is performed for a certain number of
     * pixels above and below the position of the reference wavelength:
     */
    
    pixel_above = STRETCH_FACTOR * (red - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - blue) / dispersion;

    xlow = 1;
    xhig = nx;
    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        /*
         * We DON'T write:
         *
         * ylow = position[i];
         * yhig = ylow + length[i];
         *
         * because the cpl_image pixels are counted from 1, and because in 
         * cpl_image_extract() the coordinates of the last pixel are inclusive.
         */

        ylow = position[i] + 1;
        yhig = ylow + length[i] - 1;

        exslit = cpl_image_extract(rectified, xlow, ylow, xhig, yhig);

        if (polyorder < 0) {

            cpl_image_turn(exslit, -1);   /* For faster memory access */
    
            nsubx = cpl_image_get_size_x(exslit);
            nsuby = cpl_image_get_size_y(exslit);
            data = cpl_image_get_data(exslit);
            flux = cpl_vector_new(nsubx);

            uradius = nsubx / 2;
            if (uradius > sradius)
                uradius = sradius;

            for (j = 0; j < nsuby; j++) {
                fdata = cpl_vector_get_data(flux);
                p = data;
                for (k = 0; k < nsubx; k++)
                    *fdata++ = *p++;
                smo_flux = cpl_vector_filter_median_create(flux, uradius);
                fdata = cpl_vector_get_data(smo_flux);
                p = data;
                for (k = 0; k < nsubx; k++)
                    *p++ = *fdata++;
                cpl_vector_delete(smo_flux);
                data += nsubx;
            }

            cpl_vector_delete(flux);


            /*
             * First fit fluxes along the spatial direction with a low-degree
             * polynomial (excluding the first and the last pixels, close to
             * the edges)
             */
/*
            if (nsubx-2*exclude > 10) {
                flux = cpl_vector_new(nsubx-2*exclude);
                fdata = cpl_vector_get_data(flux);
                positions = cpl_vector_new(nsubx-2*exclude);
                for (j = 0; j < nsubx-2*exclude; j++)
                    cpl_vector_set(positions, j, j+exclude);
        
                for (k = 0; k < nsuby; k++) {
                    for (j = 0; j < nsubx-2*exclude; j++)
                        fdata[j] = data[j+exclude];
                    trend = cpl_polynomial_fit_1d_create(positions, flux, 
                                                         1, NULL);
                    for (j = 0; j < nsubx; j++)
                        data[j] = cpl_polynomial_eval_1d(trend, j, NULL);
                    cpl_polynomial_delete(trend);
                    data += nsubx;
                }

                cpl_vector_delete(flux);
                cpl_vector_delete(positions);
            }
*/

            /*
             * Now smooth along the dispersion direction 
             */

            cpl_image_turn(exslit, 1);   /* For faster memory access */
            nsubx = cpl_image_get_size_x(exslit);
            nsuby = cpl_image_get_size_y(exslit);
            data = cpl_image_get_data(exslit);

            for (j = 0; j < nsuby; j++) {
                flux = cpl_vector_new(nsubx);
                fdata = cpl_vector_get_data(flux);
                p = data;
                for (k = 0; k < nsubx; k++)
                    *fdata++ = *p++;
                smo_flux = cpl_vector_filter_median_create(flux, sradius);
                cpl_vector_delete(flux);
                fdata = cpl_vector_get_data(smo_flux);
                p = data;
                for (k = 0; k < nsubx; k++)
                    *p++ = *fdata++;
                cpl_vector_delete(smo_flux);
                data += nsubx;
            }
        }
        else {

            /*
             * Fit with a polynomial the flat field trend row by row.
             */

            nsubx = cpl_image_get_size_x(exslit);
            nsuby = cpl_image_get_size_y(exslit);
            data = cpl_image_get_data(exslit);

            for (j = 0; j < nsuby; j++) {

                /*
                 * First get the size of the vectors to allocate:
                 */

                npoints = 0;
                p = data + j*nsubx;
                for (k = 0; k < nsubx; k++)
                    if (p[k] > 1.0)
                        npoints++;

                if (npoints > polyorder + 1) {

                    /*
                     * Fill the vectors for the fitting
                     */

                    flux = cpl_vector_new(npoints);
                    fdata = cpl_vector_get_data(flux);
                    positions = cpl_vector_new(npoints);
                    pdata = cpl_vector_get_data(positions);

                    npoints = 0;
                    p = data + j*nsubx;
                    for (k = 0; k < nsubx; k++) {
                        if (p[k] > 1.0) {
                            fdata[npoints] = p[k];
                            pdata[npoints] = k;
                            npoints++;
                        }
                    }
    
                    trend = cpl_polynomial_fit_1d_create(positions, flux, 
                                                         polyorder, NULL);

                    cpl_vector_delete(flux);
                    cpl_vector_delete(positions);

                    if (trend) {
                        p = data + j*nsubx;
                        for (k = 0; k < nsubx; k++)
                            if (p[k] > 1.0)
                                p[k] = cpl_polynomial_eval_1d(trend, k, NULL);
                        cpl_polynomial_delete(trend);
                    }
                    else {
                        cpl_msg_warning(func, "Invalid flat field flux fit "
                                        "(ignored)");
                    }
                }
            }
        }

        
        /*
         * Recover from the table of spectral curvature coefficients
         * the curvature polynomials.
         */

        refpixel = cpl_table_get_double(slits, "xtop", i, NULL);

        start_pixel = refpixel - pixel_below;
        if (start_pixel < 0)
            start_pixel = 0;

        end_pixel = refpixel + pixel_above;
        if (end_pixel > nx)
            end_pixel = nx;

        missing_top = 0;
        polytop = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i, &null);
            if (null) {
                cpl_polynomial_delete(polytop);
                missing_top = 1;
                break;
            }
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        missing_bot = 0;
        polybot = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i+1, &null);
            if (null) {
                cpl_polynomial_delete(polybot);
                missing_bot = 1;
                break;
            }
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        if (missing_top && missing_bot) {
            cpl_msg_debug(func, "Slit %d was not traced: no extraction!",
                          slit_id[i]);
            continue;
        }

        /*
         * In case just one of the two edges was not traced, the other
         * edge curvature model is duplicated and shifted to the other
         * end of the slit: better than nothing!
         */

        if (missing_top) {
            cpl_msg_debug(func, "Upper edge of slit %d was not traced: "
                          "the spectral curvature of the lower edge "
                          "is used instead.", slit_id[i]);
            polytop = cpl_polynomial_duplicate(polybot);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polybot, &k);
            coeff += ytop - ybot;
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        if (missing_bot) {
            cpl_msg_debug(func, "Lower edge of slit %d was not traced: "
                          "the spectral curvature of the upper edge "
                          "is used instead.", slit_id[i]);
            polybot = cpl_polynomial_duplicate(polytop);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polytop, &k);
            coeff -= ytop - ybot;
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }


        /*
         * Now map smoothed image to CCD.
         * Note that the npseudo value related to this slit is equal
         * to the number of spatial pseudo-pixels decreased by 1
         * (compare with function mos_spatial_calibration()).
         */

        nx = cpl_image_get_size_x(flat);
        ny = cpl_image_get_size_y(flat);

        sdata = cpl_image_get_data(spatial);
        xdata = cpl_image_get_data(exslit);
        npseudo = cpl_image_get_size_y(exslit) - 1;

        /*
         * Write interpolated smoothed values to CCD image
         */

        for (j = start_pixel; j < end_pixel; j++) {
            top = cpl_polynomial_eval_1d(polytop, j, NULL);
            bot = cpl_polynomial_eval_1d(polybot, j, NULL);
            for (k = 0; k <= npseudo; k++) {
                ypos = top - k*(top-bot)/npseudo;
                yint = ypos;

                /*
                 * The line:
                 *     value = sdata[j + nx*yint];
                 * should be equivalent to:
                 *     value = npseudo*(top-yint)/(top-bot);
                 */

                if (yint < 0 || yint >= ny-1) {
                    yprev = yint;
                    continue;
                }

                value = sdata[j + nx*yint];   /* Spatial coordinate on CCD */
                ivalue = value;               /* Nearest spatial pixels:   */
                fvalue = value - ivalue;      /* ivalue and ivalue+1       */
                if (ivalue < npseudo && ivalue >= 0) {
                    vtop = xdata[j + nx*(npseudo-ivalue)];
                    vbot = xdata[j + nx*(npseudo-ivalue-1)];
                    wdata[j + nx*yint] = vtop*(1-fvalue) + vbot*fvalue;

                    if (k) {

                        /*
                         * This is added to recover lost pixels on
                         * the CCD image (pixels are lost because
                         * the CCD pixels are less than npseudo+1).
                         */

                        if (yprev - yint > 1) {
                            value = sdata[j + nx*(yint+1)];
                            ivalue = value;
                            fvalue = value - ivalue;
                            if (ivalue < npseudo && ivalue >= 0) {
                                vtop = xdata[j + nx*(npseudo-ivalue)];
                                vbot = xdata[j + nx*(npseudo-ivalue-1)];
                                wdata[j + nx*(yint+1)] = vtop*(1-fvalue) 
                                                       + vbot*fvalue;
                            }
                        }
                    }
                }
                yprev = yint;
            }
        }
        cpl_polynomial_delete(polytop);
        cpl_polynomial_delete(polybot);
        cpl_image_delete(exslit);
    }

    cpl_image_delete(rectified);

    cpl_image_divide(flat, smo_flat);

    return smo_flat;
}


/**
 * @brief 
 *   Normalise a long slit flat field exposure
 *  
 * @param flat        Image containing the original flat field spectra
 * @param sradius     Radius of smoothing box along the spatial direction
 * @param dradius     Radius of smoothing box along the dispersion direction
 * @param polyorder   Order of fitting polynomial along the spatial direction
 *  
 * @return The smoothed flat field exposure used for normalisation
 *
 * The input @em flat frame should be already bias subtracted, and should
 * be oriented so that the dispersion direction is horizontal with @em blue
 * on the left and @em red on the right. The original @em flat image is 
 * divided IN PLACE by its smoothed counterpart, which is also returned. 
 * If the polynomial @em polyorder is set to a negative number the smoothing 
 * consists of a median filtering box of specified sizes. Alternatively, 
 * if @em polyorder is not negative, the smoothing will consist of a 
 * polynomial fitting of the illumination profile along the spatial 
 * direction (and not along the dispersion direction, as for the case 
 * of shorter slits), performed independently for each column of the 
 * spectrum.
 */
 
cpl_image *mos_normalise_longflat(cpl_image *flat, int sradius, int dradius, 
                                  int polyorder)
{
    const char     *func = "mos_normalise_longflat";

    cpl_image      *smo_flat;
    cpl_image      *profile;
    cpl_vector     *flux;
    cpl_vector     *smo_flux;
    cpl_vector     *positions;
    cpl_polynomial *trend;

    float          *level;
    float          *p;
    float          *data;
    double         *fdata;
    double         *pdata;

    int             nx, ny;
    int             npoints;
    int             i, j;


    if (flat == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }
 
    if (sradius < 1 || dradius < 1) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    smo_flat = cpl_image_duplicate(flat);

    if (polyorder < 0) {

        /*
         * First smooth along the spatial direction
         */

        cpl_image_turn(smo_flat, -1);   /* For faster memory access */

        nx = cpl_image_get_size_x(smo_flat);
        ny = cpl_image_get_size_y(smo_flat);
        data = cpl_image_get_data(smo_flat);
    
        for (i = 0; i < ny; i++) {
            flux = cpl_vector_new(nx);
            fdata = cpl_vector_get_data(flux);
            p = data;
            for (j = 0; j < nx; j++)
                *fdata++ = *p++;
            smo_flux = cpl_vector_filter_median_create(flux, sradius);
            cpl_vector_delete(flux);
            fdata = cpl_vector_get_data(smo_flux);
            p = data;
            for (j = 0; j < nx; j++)
                *p++ = *fdata++;
            cpl_vector_delete(smo_flux);
            data += nx;
        }

        /*
         * Second smooth along the dispersion direction
         */

        cpl_image_turn(smo_flat, 1);   /* For faster memory access */

        nx = cpl_image_get_size_x(smo_flat);
        ny = cpl_image_get_size_y(smo_flat);
        data = cpl_image_get_data(smo_flat);

        for (i = 0; i < ny; i++) {
            flux = cpl_vector_new(nx);
            fdata = cpl_vector_get_data(flux);
            p = data;
            for (j = 0; j < nx; j++)
                *fdata++ = *p++;
            smo_flux = cpl_vector_filter_median_create(flux, sradius);
            cpl_vector_delete(flux);
            fdata = cpl_vector_get_data(smo_flux);
            p = data;
            for (j = 0; j < nx; j++)
                *p++ = *fdata++;
            cpl_vector_delete(smo_flux);
            data += nx;
        }
    }
    else {

        /*
         * Fit with a polynomial the flat field trend column by column.
         */

        cpl_image_turn(smo_flat, -1);   /* For faster memory access */

        nx = cpl_image_get_size_x(smo_flat);
        ny = cpl_image_get_size_y(smo_flat);
        data = cpl_image_get_data(smo_flat);

        profile = cpl_image_collapse_median_create(smo_flat, 1, 0, 0);
        level = cpl_image_get_data(profile);

        for (i = 0; i < ny; i++) {

            /*
             * First get the size of the vectors to allocate:
             * eliminate from fit any value more than 20% away
             * from median level in current column.
             */

            npoints = 0;
            p = data + i*nx;
            for (j = 0; j < nx; j++)
                if (fabs(p[j]/level[i] - 1) < 0.20)
                    npoints++;

            if (npoints > polyorder + 1) {

                /*
                 * Fill the vectors for the fitting
                 */

                flux = cpl_vector_new(npoints);
                fdata = cpl_vector_get_data(flux);
                positions = cpl_vector_new(npoints);
                pdata = cpl_vector_get_data(positions);

                npoints = 0;
                p = data + i*nx;
                for (j = 0; j < nx; j++) {
                    if (fabs(p[j]/level[i] - 1) < 0.20) {
                        fdata[npoints] = p[j];
                        pdata[npoints] = j;
                        npoints++;
                    }
                }
    
                trend = cpl_polynomial_fit_1d_create(positions, flux, 
                                                     polyorder, NULL);

                cpl_vector_delete(flux);
                cpl_vector_delete(positions);

                if (trend) {
                    p = data + i*nx;
                    for (j = 0; j < nx; j++)
                        p[j] = cpl_polynomial_eval_1d(trend, j, NULL);
                    cpl_polynomial_delete(trend);
                }
                else {
                    cpl_msg_warning(func, 
                                    "Invalid flat field flux fit (ignored)");
                }
            }
        }

        cpl_image_delete(profile);
        cpl_image_turn(smo_flat, 1);

    }

    cpl_image_divide(flat, smo_flat);

    return smo_flat;
}


/**
 * @brief 
 *   Interpolate MOS wavelength calibration
 *  
 * @param idscoeff    Table with IDS polynomials
 * @param slits       Table with slits positions
 * @param order       0 = median, > 0 order of fitting polynomial.
 *  
 * @return @c CPL_ERROR_NONE on success.
 *
 * This function is used on MOS data, to interpolate the wavelength
 * calibration obtained with the function @c mos_wavelength_calibration_final()
 * also in those image rows where the calibration failed. The input 
 * @em idscoeff table is reprocessed and modified in-place according 
 * to what indicated by the @em mode argument. The @em idscoeff table 
 * coefficients are modeled by low degree polynomials: if @em mode is 1, 
 * the model values are used to fill the gaps in the input solutions, 
 * otherwise all the inputs are replaced with the fitted model values. 
 * The corresponding wavelength map could be calculated using the 
 * function @c mos_map_idscoeff().
 */

cpl_error_code mos_interpolate_wavecalib_slit(cpl_table *idscoeff,
                                              cpl_table *slits, 
                                              int order, int global)
{
    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    int nrow = cpl_table_get_nrow(slits);
    int i, j;

    
    if (order < 0)
        return CPL_ERROR_NONE;

    cpl_table_new_column(idscoeff, "x", CPL_TYPE_DOUBLE);
    cpl_table_new_column(idscoeff, "y", CPL_TYPE_DOUBLE);

    for (i = 0; i < nrow; i++) {
        int        position = cpl_table_get_int   (slits, "position", i, NULL);
        int        length   = cpl_table_get_int   (slits, "length",   i, NULL);
        double     xtop     = cpl_table_get_double(slits, "xtop",     i, NULL);
        double     xbot     = cpl_table_get_double(slits, "xbottom",  i, NULL);
        double     ytop     = cpl_table_get_double(slits, "ytop",     i, NULL);
        double     ybot     = cpl_table_get_double(slits, "ybottom",  i, NULL);
        double     dx       = xtop - xbot;
        double     dy       = ytop - ybot;
        cpl_table *table    = cpl_table_extract(idscoeff, position, length);

        if (mos_interpolate_wavecalib(table, NULL, 2, order))
            continue;

        cpl_table_erase_window(idscoeff, position, length);
        cpl_table_insert(idscoeff, table, position);

        cpl_table_delete(table);

        for (j = 0; j < length; j++) {
            cpl_table_set_double(idscoeff, "x", j + position,
                                 xbot + j*(dx/length));
            cpl_table_set_double(idscoeff, "y", j + position,
                                 ybot + j*(dy/length));
        }
    }

    if (global) {

        /*
         * Now fit a global solution
         */

        nrow = cpl_table_get_nrow(idscoeff);

        for (i = 0; i < 6; i++) {
            cpl_table      *dummy;
            cpl_vector     *x;
            cpl_vector     *y;
            cpl_bivector   *z;
            cpl_vector     *c;
            cpl_polynomial *p;
            cpl_vector     *point;
            double         *dpoint;
            int             npoints;

            if (!cpl_table_has_column(idscoeff, clab[i]))
                break;

            npoints = nrow - cpl_table_count_invalid(idscoeff, clab[i]);
            if (npoints < 18)
                break;

            dummy = cpl_table_new(nrow);
            cpl_table_duplicate_column(dummy, "x", idscoeff, "x");
            cpl_table_duplicate_column(dummy, "y", idscoeff, "y");
            cpl_table_duplicate_column(dummy, clab[i], idscoeff, clab[i]);
            cpl_table_erase_invalid(dummy);

            x = cpl_vector_wrap(npoints, cpl_table_get_data_double(dummy, "x"));
            y = cpl_vector_wrap(npoints, cpl_table_get_data_double(dummy, "y"));
            z = cpl_bivector_wrap_vectors(x, y);
            c = cpl_vector_wrap(npoints, cpl_table_get_data_double(dummy, 
                                                                   clab[i]));
            p = cpl_polynomial_fit_2d_create(z, c, 2, NULL);
            cpl_bivector_unwrap_vectors(z);
            cpl_vector_unwrap(x);
            cpl_vector_unwrap(y);
            cpl_vector_unwrap(c);
            cpl_table_delete(dummy);

            point  = cpl_vector_new(2);
            dpoint = cpl_vector_get_data(point);
            for (j = 0; j < nrow; j++) {
                dpoint[0] = cpl_table_get_double(idscoeff, "x", j, NULL);
                dpoint[1] = cpl_table_get_double(idscoeff, "y", j, NULL);
                cpl_table_set_double(idscoeff, clab[i], j,
                                     cpl_polynomial_eval(p, point));
            }
            cpl_vector_delete(point);
            cpl_polynomial_delete(p);
        }
    }

    return CPL_ERROR_NONE;
}


/**
 * @brief 
 *   Interpolate LSS wavelength calibration
 *  
 * @param idscoeff    Table with IDS polynomials
 * @param wavemap     Wavelength calibration image
 * @param mode        0 = do nothing, 1 = fill gaps, 2 = global model
 * @param degree      0 = median, > 0 order of fitting polynomial.
 *  
 * @return @c CPL_ERROR_NONE on success.
 *
 * This function is used on LSS data, to interpolate the wavelength
 * calibration obtained with the function @c mos_wavelength_calibration_raw()
 * also in those image rows where the calibration failed. The input 
 * @em idscoeff table and @em wavemap image are reprocessed and modified
 * in-place according to what indicated by the @em mode argument. 
 * The @em idscoeff table coefficients and the wavelengths in @em wavemap
 * are modeled by low degree polynomials: if @em mode is 1, the model
 * values are used to fill the gaps in the input solutions, otherwise
 * all the inputs are replaced with the fitted model values. If the
 * @em wavemap is not given, just the @em idscoeff table will be
 * interpolated: the corresponding wavelength map could still be 
 * calculated using the function @c mos_map_idscoeff().
 */
 
cpl_error_code mos_interpolate_wavecalib(cpl_table *idscoeff, 
                                         cpl_image *wavemap, int mode,
                                         int degree)
{
    const char *func = "mos_interpolate_wavecalib";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_vector     *wave;
    cpl_vector     *positions;
    cpl_polynomial *trend;

    float          *p;
    float          *data;
    double         *wdata;
    double         *pdata;

    double          c;
    double          mse, ksigma;

    int             order;
    int             nrows, first_row, last_row;
    int             nx, ny;
    int             npoints, rpoints;
    int             null;
    int             i, j, k;

    int             polyorder = 4;  /* Candidate input argument */


    if (idscoeff == NULL)
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);

    if (mode < 0 || mode > 2)
        return cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);

    if (mode == 0 || degree < 0)
        return CPL_ERROR_NONE;

    if (wavemap) {

        /*
         * Fit with a polynomial the wavelength trend column by column.
         */

        cpl_image_turn(wavemap, -1);   /* For faster memory access */

        nx = cpl_image_get_size_x(wavemap);
        ny = cpl_image_get_size_y(wavemap);
        data = cpl_image_get_data(wavemap);

        for (i = 0; i < ny; i++) {

            /*
             * First get the size of the vectors to allocate:
             * eliminate from fit any value with "impossible" wavelength.
             */

            npoints = 0;
            p = data + i*nx;
            for (j = 0; j < nx; j++)
                if (p[j] > 1.0)
                    npoints++;

            if (npoints > polyorder + 1) {

                /*
                 * Fill the vectors for the fitting
                 */

                wave = cpl_vector_new(npoints);
                wdata = cpl_vector_get_data(wave);
                positions = cpl_vector_new(npoints);
                pdata = cpl_vector_get_data(positions);

                npoints = 0;
                p = data + i*nx;
                for (j = 0; j < nx; j++) {
                    if (p[j] > 1.0) {
                        wdata[npoints] = p[j];
                        pdata[npoints] = j;
                        npoints++;
                    }
                }
    
                trend = cpl_polynomial_fit_1d_create(positions, wave, 
                                                     polyorder, &mse);

                ksigma = 3*sqrt(mse);

                cpl_vector_delete(wave);
                cpl_vector_delete(positions);

                if (trend) {

                    /*
                     * Apply 3-sigma rejection
                     */

                    rpoints = 0;
                    p = data + i*nx;
                    for (j = 0; j < nx; j++)
                        if (p[j] > 1.0)
                            if (fabs(cpl_polynomial_eval_1d(trend, j, NULL) 
                                                    - p[j]) < ksigma)
                                rpoints++;

                    if (rpoints < npoints && rpoints > polyorder + 1) {

                        wave = cpl_vector_new(rpoints);
                        wdata = cpl_vector_get_data(wave);
                        positions = cpl_vector_new(rpoints);
                        pdata = cpl_vector_get_data(positions);

                        npoints = 0;
                        p = data + i*nx;
                        for (j = 0; j < nx; j++) {
                            if (p[j] > 1.0) {
                                if (fabs(cpl_polynomial_eval_1d(trend, 
                                                                j, NULL) - p[j])
                                                                < ksigma) {
                                    wdata[npoints] = p[j];
                                    pdata[npoints] = j;
                                    npoints++;
                                }
                            }
                        }
        
                        cpl_polynomial_delete(trend);
                        trend = cpl_polynomial_fit_1d_create(positions, wave,
                                                             polyorder, NULL);

                        cpl_vector_delete(wave);
                        cpl_vector_delete(positions);
                    }
                }

                if (trend) {
                    p = data + i*nx;
                    if (mode == 1) {
                        for (j = 0; j < nx; j++)
                            if (p[j] < 1.0)
                                p[j] = cpl_polynomial_eval_1d(trend, j, NULL);
                    }
                    else if (mode == 2) {
                        for (j = 0; j < nx; j++)
                            p[j] = cpl_polynomial_eval_1d(trend, j, NULL);
                    }
                    cpl_polynomial_delete(trend);
                }
                else {
                    cpl_msg_warning(func, 
                                    "Invalid wavelength field fit (ignored)");
                }
            }
    
        }

        cpl_image_turn(wavemap, 1);

    }


    /*
     * Interpolating the IDS coefficients
     */

    nrows = cpl_table_get_nrow(idscoeff);

    order = 0;
    while (order < 6 && cpl_table_has_column(idscoeff, clab[order]))
        ++order;
    --order;

    if (degree == 0) {
        for (k = 0; k <= order; k++) {
            double m;
            if (cpl_table_has_column(idscoeff, clab[k])) {
                m = cpl_table_get_column_median(idscoeff, clab[k]);
                cpl_table_fill_column_window_double(idscoeff, clab[k], 
                                                    0, nrows, m);
            }
        }

        return CPL_ERROR_NONE;
    }

    first_row = 0;
//    while (!cpl_table_is_valid(idscoeff, clab[0], first_row))
//        first_row++;

    last_row = nrows - 1;
//    while (!cpl_table_is_valid(idscoeff, clab[0], last_row))
//        last_row--;

// The lines above are commented to enable extrapolation of the solution

//{
//bla++;
//char *tablename = cpl_sprintf("before%d.fits", bla);
//cpl_table_save(idscoeff, NULL, NULL, tablename, CPL_IO_DEFAULT);
//cpl_free(tablename);
//}
//
//    for (k = 0; k <= order; k++) {
//
//        /*
//         * Eliminate obvious outliers from fit. 
//         */
//
//        mos_clean_outliers(idscoeff, clab[k]);
//    }
//
//{
//char *tablename = cpl_sprintf("after%d.fits", bla);
//cpl_table_save(idscoeff, NULL, NULL, tablename, CPL_IO_DEFAULT);
//cpl_free(tablename);
//}

    for (k = 0; k <= order; k++) {

        npoints = nrows - cpl_table_count_invalid(idscoeff, clab[k]);
        wave = cpl_vector_new(npoints);
        wdata = cpl_vector_get_data(wave);
        positions = cpl_vector_new(npoints);
        pdata = cpl_vector_get_data(positions);

        npoints = 0;
        for (i = first_row; i <= last_row; i++) {
            c = cpl_table_get_double(idscoeff, clab[k], i, &null);
            if (null == 0) {
                wdata[npoints] = c;
                pdata[npoints] = i;
                npoints++;
            }
        }

        if (degree == 1) {
            cpl_table *points = cpl_table_new(npoints);
            cpl_table *hough;
            double     q, m;

            cpl_table_wrap_double(points, pdata, "p");
            cpl_table_wrap_double(points, wdata, "w");
            hough = mos_hough_table(points, "p", "w");
            cpl_table_unwrap(points, "p");
            cpl_table_unwrap(points, "w");
            cpl_table_delete(points);

            if (hough == NULL)
                continue;

            m = cpl_table_get_column_median(hough, "m");
            q = cpl_table_get_column_median(hough, "q");

            cpl_table_delete(hough);
            
            /*
             * If degree = 1, dare an extrapolation to the whole 
             * slit. Before it was:
             *
             *  for (i = first_row; i <= last_row; i++)
             */

            for (i = 0; i < nrows; i++)
                 cpl_table_set_double(idscoeff, clab[k], i, q + m*i);

            continue;
        }

// This doesn't seem to provide good results, I have not understood why.
// Restore for robust linear fitting.
//
//        if (degree == 1) {
//            cpl_vector   *p;
//            cpl_vector   *w;
//            cpl_bivector *list;
//            double        q, m;
//
//            if (npoints > 4) {
//                p = cpl_vector_extract(positions, 2, npoints - 2, 1);
//                w = cpl_vector_extract(wave, 2, npoints - 2, 1);
//            }
//            else {
//                p = positions;
//                w = wave;
//            }
//
//            list = cpl_bivector_wrap_vectors(p, w);
//
//            robustLinearFit(list, &q, &m, &mse);
//            cpl_bivector_unwrap_vectors(list);
//            for (i = first_row; i <= last_row; i++)
//                 cpl_table_set_double(idscoeff, clab[k], i, q + m*i);
//
//            if (npoints > 4) {
//                cpl_vector_delete(p);
//                cpl_vector_delete(w);
//            }
//
//            continue;
//        }

        trend = cpl_polynomial_fit_1d_create(positions, wave, degree, &mse);

        ksigma = 3*sqrt(mse);

        cpl_vector_delete(wave);
        cpl_vector_delete(positions);

        /*
         * Iteration
         */

        if (trend) {
            rpoints = 0;
            for (i = first_row; i <= last_row; i++) {
                c = cpl_table_get_double(idscoeff, clab[k], i, &null);
                if (null == 0) {
                    if (fabs(cpl_polynomial_eval_1d(trend, i, NULL) - c) 
                                                                 < ksigma) {
                        rpoints++;
                    }
                }
            }

            if (rpoints > 0 && rpoints < npoints) {
                cpl_msg_debug(func, "%d points rejected from "
                              "wavelength calibration fit", 
                              npoints - rpoints);

                wave = cpl_vector_new(rpoints);
                wdata = cpl_vector_get_data(wave);
                positions = cpl_vector_new(rpoints);
                pdata = cpl_vector_get_data(positions);

                npoints = 0;
                for (i = first_row; i <= last_row; i++) {
                    c = cpl_table_get_double(idscoeff, clab[k], i, &null);
                    if (null == 0) {
                        if (fabs(cpl_polynomial_eval_1d(trend, i, NULL) - c)
                                                                 < ksigma) {
                            wdata[npoints] = c;
                            pdata[npoints] = i;
                            npoints++;
                        }
                    }
                }

                if (npoints) {
                    cpl_polynomial_delete(trend);
                    trend = cpl_polynomial_fit_1d_create(positions, 
                                                         wave, degree, NULL);
                }

                cpl_vector_delete(wave);
                cpl_vector_delete(positions);

            }
        }

        if (trend) {
            for (i = first_row; i <= last_row; i++) {
                if (mode == 1) {
                    if (!cpl_table_is_valid(idscoeff, clab[k], i)) {
                        cpl_table_set_double(idscoeff, clab[k], i, 
                                             cpl_polynomial_eval_1d(trend, i,
                                                                    NULL));
                    }
                }
                else if (mode == 2) {
                    cpl_table_set_double(idscoeff, clab[k], i, 
                                     cpl_polynomial_eval_1d(trend, i, NULL));
                }
            }
            cpl_polynomial_delete(trend);
        }
        else {
            cpl_msg_warning(func, "Invalid IDS coefficient fit (ignored)");
        }

    }

    return CPL_ERROR_NONE;
}

/**
 * @brief 
 *   Interpolate wavelength calibration for a single MOS slit
 *  
 * @param idscoeff    Table with IDS polynomials
 * @param wavemap     Wavelength calibration image
 * @param mode        0 = do nothing, 1 = fill gaps, 2 = global model
 * @param degree      0 = median, > 0 order of fitting polynomial. This applies
 *                    only to the constant coefficient
 *  
 * @return @c CPL_ERROR_NONE on success.
 *
 * This function is a copy of mos_interpolate_wavecalib which does the
 * fitting only for the first coefficient, while for the rest it just
 * computes the median
 */
 
cpl_error_code mos_interpolate_wavecalib_mos(cpl_table *idscoeff, 
                                             int mode,
                                             int degree)
{
    const char *func = "mos_interpolate_wavecalib";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_vector     *wave;
    cpl_vector     *positions;
    cpl_polynomial *trend;

    double         *wdata;
    double         *pdata;

    double          c;
    double          mse, ksigma;

    int             order;
    int             nrows, first_row, last_row;
    int             npoints, rpoints;
    int             null;
    int             i;


    if (idscoeff == NULL)
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);

    if (mode < 0 || mode > 2)
        return cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);

    if (mode == 0 || degree < 0)
        return CPL_ERROR_NONE;

    /*
     * Interpolating the IDS coefficients
     */

    nrows = cpl_table_get_nrow(idscoeff);

    order = 0;
    while (order < 6 && cpl_table_has_column(idscoeff, clab[order]))
        ++order;
    --order;

    //For coefficients > 0 always do a median 
    
    int kinit = 1;
    if (degree == 0) 
        kinit = 0;
    for (int k = kinit; k <= order; k++) {
        double m;
        if (cpl_table_has_column(idscoeff, clab[k])) {
            m = cpl_table_get_column_median(idscoeff, clab[k]);
            cpl_table_fill_column_window_double(idscoeff, clab[k], 
                                                0, nrows, m);
        }
    }

    if(degree > 0)
    {
        first_row = 0;
        while (!cpl_table_is_valid(idscoeff, clab[0], first_row))
            first_row++;

        last_row = nrows - 1;
        while (!cpl_table_is_valid(idscoeff, clab[0], last_row))
            last_row--;

        int korder = 0;

        npoints = nrows - cpl_table_count_invalid(idscoeff, clab[korder]);
        wave = cpl_vector_new(npoints);
        wdata = cpl_vector_get_data(wave);
        positions = cpl_vector_new(npoints);
        pdata = cpl_vector_get_data(positions);

        npoints = 0;
        for (i = first_row; i <= last_row; i++) {
            c = cpl_table_get_double(idscoeff, clab[korder], i, &null);
            if (null == 0) {
                wdata[npoints] = c;
                pdata[npoints] = i;
                npoints++;
            }
        }

        // This doesn't seem to provide good results, I have not understood why.
        // Restore for robust linear fitting.
        //
        if (degree == 1) {
            cpl_vector   *p;
            cpl_vector   *w;
            cpl_bivector *list;
            double        q, m;

            if (npoints > 4) {
                p = cpl_vector_extract(positions, 2, npoints - 2, 1);
                w = cpl_vector_extract(wave, 2, npoints - 2, 1);
            }
            else {
                p = positions;
                w = wave;
            }

            list = cpl_bivector_wrap_vectors(p, w);

            robustLinearFit(list, &q, &m, &mse);
            cpl_bivector_unwrap_vectors(list);
            for (i = first_row; i <= last_row; i++)
                cpl_table_set_double(idscoeff, clab[korder], i, q + m*i);

            if (npoints > 4) {
                cpl_vector_delete(p);
                cpl_vector_delete(w);
            }

            cpl_vector_delete(wave);
            cpl_vector_delete(positions);

        }
        else
        {

            // End robust linear fitting

            trend = cpl_polynomial_fit_1d_create(positions, wave, degree, &mse);

            ksigma = 3*sqrt(mse);

            cpl_vector_delete(wave);
            cpl_vector_delete(positions);

            /*
             * Iteration
             */

            if (trend) {
                rpoints = 0;
                for (i = first_row; i <= last_row; i++) {
                    c = cpl_table_get_double(idscoeff, clab[korder], i, &null);
                    if (null == 0) {
                        if (fabs(cpl_polynomial_eval_1d(trend, i, NULL) - c) 
                                < ksigma) {
                            rpoints++;
                        }
                    }
                }

                if (rpoints > 0 && rpoints < npoints) {
                    cpl_msg_debug(func, "%d points rejected from "
                            "wavelength calibration fit", 
                            npoints - rpoints);

                    wave = cpl_vector_new(rpoints);
                    wdata = cpl_vector_get_data(wave);
                    positions = cpl_vector_new(rpoints);
                    pdata = cpl_vector_get_data(positions);

                    npoints = 0;
                    for (i = first_row; i <= last_row; i++) {
                        c = cpl_table_get_double(idscoeff, clab[korder], i, &null);
                        if (null == 0) {
                            if (fabs(cpl_polynomial_eval_1d(trend, i, NULL) - c)
                                    < ksigma) {
                                wdata[npoints] = c;
                                pdata[npoints] = i;
                                npoints++;
                            }
                        }
                    }

                    if (npoints) {
                        cpl_polynomial_delete(trend);
                        trend = cpl_polynomial_fit_1d_create(positions, 
                                wave, degree, NULL);
                    }

                    cpl_vector_delete(wave);
                    cpl_vector_delete(positions);

                }
            }

            if (trend) {
                for (i = first_row; i <= last_row; i++) {
                    if (mode == 1) {
                        if (!cpl_table_is_valid(idscoeff, clab[korder], i)) {
                            cpl_table_set_double(idscoeff, clab[korder], i, 
                                    cpl_polynomial_eval_1d(trend, i,
                                            NULL));
                        }
                    }
                    else if (mode == 2) {
                        cpl_table_set_double(idscoeff, clab[korder], i, 
                                cpl_polynomial_eval_1d(trend, i, NULL));
                    }
                }
                cpl_polynomial_delete(trend);
            }
            else {
                cpl_msg_warning(func, "Invalid IDS coefficient fit (ignored)");
            }
        }
    }

    return CPL_ERROR_NONE;
}


/**
 * @brief
 *   Subtract the bias from a CCD exposure
 *
 * @param image       Image containing the data to correct
 * @param bias        Master bias
 * @param overscans   Table with the overscan information
 *
 * @return A newly allocated bias subtracted image
 *
 * If the master @em bias has the same sizes of the input @em image it is 
 * simply subtracted from the @em image. If it has different sizes, they 
 * should be compatible with the description in the @em overscans table 
 * produced by one of the functions @em mos_load_overscans_<instrument>(). 
 * In this case the mean level of the master bias is compared with the 
 * mean level of the overscan regions of the input image. The master 
 * @em bias and the difference of mean levels are then subtracted from 
 * the valid region of the input @em image. The output image will always 
 * have the overscan regions trimmed, having in this way the same sizes 
 * of the input master @em bias. If the @em bias frame is not specified,
 * only the mean level of the overscan regions is subtracted from the
 * input image. Only in case no overscan regions are present, nothing
 * is done.
 */

cpl_image *mos_remove_bias(cpl_image *image, cpl_image *bias, 
                           cpl_table *overscans)
{
    const char *func = "mos_remove_bias";

    cpl_image *unbiased;
    cpl_image *overscan;
    double     mean_bias_level;
    double     mean_overscans_level;
    int        count;
    int        nrows;
    int        xlow, ylow, xhig, yhig;
    int        i;


    if (image == NULL || overscans == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nrows = cpl_table_get_nrow(overscans);

    if (nrows == 0) {
        cpl_msg_error(func, "Empty overscan table");
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    if (bias) {
        if (nrows == 1) {
            unbiased = cpl_image_subtract_create(image, bias);
            if (unbiased == NULL) {
                cpl_msg_error(func, "Incompatible master bias");
                cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
            }
            return unbiased;
        }
        mean_bias_level = cpl_image_get_mean(bias);
    }
    else {
        if (nrows == 1) {
            cpl_msg_error(func, "No master bias in input, and no overscan "
                          "regions in input image: bias subtraction "
                          "cannot be performed!");
            cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
            return NULL;
        }
        mean_bias_level = 0.0;
    }

    mean_overscans_level = 0.0;
    count = 0;
    for (i = 0; i < nrows; i++) {
        xlow = cpl_table_get_int(overscans, "xlow", i, NULL);
        ylow = cpl_table_get_int(overscans, "ylow", i, NULL);
        xhig = cpl_table_get_int(overscans, "xhig", i, NULL);
        yhig = cpl_table_get_int(overscans, "yhig", i, NULL);

        if (i == 0) {
            unbiased = cpl_image_extract(image, xlow+1, ylow+1, xhig, yhig);
            if (unbiased == NULL) {
                cpl_msg_error(func, "Incompatible overscan table");
                cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
                return NULL;
            }
            if (bias) {
                if (cpl_image_subtract(unbiased, bias)) {
                    cpl_msg_error(func, "Incompatible master bias");
                    cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
                    cpl_image_delete(unbiased);
                    return NULL;
                }
            }
        }
        else {
            overscan = cpl_image_extract(image, xlow+1, ylow+1, xhig, yhig);
            if (overscan == NULL) {
                cpl_msg_error(func, "Incompatible overscan table");
                cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
                cpl_image_delete(unbiased);
                return NULL;
            }

            mean_overscans_level += cpl_image_get_median(overscan);
            count++;

/***
 * Here the mean level was used: not very robust...

            mean_overscans_level += cpl_image_get_flux(overscan);
            count += cpl_image_get_size_x(overscan)
                   * cpl_image_get_size_y(overscan);
***/
            cpl_image_delete(overscan);
        }
    }

    /*
     * Overscan correction
     */

    mean_overscans_level /= count;

    cpl_image_subtract_scalar(unbiased, mean_overscans_level - mean_bias_level);

    cpl_msg_info(cpl_func, 
                 "Difference between mean overscans level "
                 "and mean bias level: %.2f",
                 mean_overscans_level - mean_bias_level);

    return unbiased;

}


/**
 * @brief
 *   Background determination on 1D emission line spectrum (arc)
 *
 * @param spectrum A 1D emission line spectrum
 * @param back     A pre-allocated buffer where the background will be written
 * @param length   Length of spectrum
 * @param msize    Size of min-filter
 * @param fsize    Size of running-average-filter
 *
 * @return CPL_ERROR_NONE in case of success
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The input spectrum is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_ILLEGAL_INPUT</td>
 *       <td class="ecr">
 *         Either <i>msize</i> is less than 3, or <i>fsize</i> is less 
 *         than <i>msize</i>, or <i>fsize</i> is greater than <i>length/2</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * This function fills the array @em back with the estimated values
 * of the background along the input spectrum. The algorithm is
 * based on the assumption that there is at least one background 
 * value at any position of the min-filter box running along the 
 * spectrum. A min-filter is passed on the spectrum, and the result 
 * is smoothed by averaging on a running box of size @em fsize.
 * The min-filter is run between the positions @em msize / 2
 * and @em length - @em msize / 2, and the min value found at
 * such positions is then repeated up to the spectrum ends. Similarly, 
 * the running average is limited to the interval from @em fsize / 2 
 * and @em length - @em fsize / 2, leaving the most external values 
 * untouched. After this, a max filter and a smoothing using boxes 
 * with double the specified sizes are run, as a way to eliminate 
 * the contamination from occasional cold pixels. Finally, the 
 * min filter and the smoothing are applied again to obviate the
 * slight background over-estimation introduced by the max filter.
 *
 * It is required that the @em back array is at least long as the
 * array @em spectrum.  Moreover @em msize must be greater than 1, 
 * and @em fsize greater than, or equal to, @em msize. Likewise, 
 * @em length must be greater than twice @em fsize. If such conditions
 * are not met, or if the input arrays are @c NULL pointers, this 
 * function will set an error code, and leave the @em back array 
 * untouched. If either @em msize or @em fsize are even numbers, 
 * they are made odd by adding 1. Suggested values for @em msize 
 * and @em fsize are 15 pixels for typical arc lamp spectra.
 */

cpl_error_code mos_arc_background_1D(float *spectrum, float *back, 
                                     int length, int msize, int fsize) 
{
    const char *func = "mos_arc_background_1D";

    float  *minf;
    float  *maxf;
    float  *smof;
    int     i;


    if (spectrum == NULL || back == NULL)
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);

    if (msize % 2 == 0)
        msize++;

    if (fsize % 2 == 0)
        fsize++;

    if (msize < 3 || fsize < msize || length < 2*fsize)
        return cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);


    minf = min_filter(spectrum, length, msize);
    smof = smo_filter(minf, length, fsize);
    cpl_free(minf);
    maxf = max_filter(smof, length, 2*msize+1);
    cpl_free(smof);
    smof = smo_filter(maxf, length, 2*fsize+1);
    cpl_free(maxf);
    minf = min_filter(smof, length, 2*msize+1);
    cpl_free(smof);
    smof = smo_filter(minf, length, 2*fsize+1);
    cpl_free(minf);

    for (i = 0; i < length; i++)
        back[i] = smof[i];

    cpl_free(smof);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Background determination on emission line spectrum (arc)
 *
 * @param image    An emission line spectrum exposure
 * @param msize    Size of min-filter
 * @param fsize    Size of running-average-filter
 *
 * @return Background image.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The input image is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_ILLEGAL_INPUT</td>
 *       <td class="ecr">
 *         Either <i>msize</i> is less than 3, or <i>fsize</i> is less 
 *         than <i>msize</i>, or <i>fsize</i> is greater than half the 
 *         input image size in the X direction.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The input @em image is assumed to be of type @c float and to contain 
 * MOS arc lamp line spectra dispersed (roughly) along the X direction. 
 * The background is estimated independently for each image row. The 
 * algorithm is based on the assumption that there is at least one 
 * background value at any position of the 1D-min-filter box running 
 * along the image row. A min-filter is passed along the row, and the 
 * result is smoothed by averaging on a running box of size @em fsize. 
 * The min-filter is run between the positions @em msize / 2 and 
 * @em length - @em msize / 2, and the min value found at such positions 
 * is then repeated up to the spectrum ends. Similarly, the running average 
 * is limited to the interval from @em fsize / 2 and @em length - @em fsize / 2,
 * leaving the most external values untouched. After this, a max filter 
 * and a smoothing using boxes with double the specified sizes are run, 
 * as a way to eliminate the contamination from occasional cold pixels. 
 * Finally, the min filter and the smoothing are applied again to obviate 
 * the slight background over-estimation introduced by the max filter.
 *
 * It is required that @em msize is greater than 1, and @em fsize greater 
 * than, or equal to, @em msize. Likewise, the image size along the X
 * direction must be greater than @em fsize / 2. If such conditions are not 
 * met, or if the input image is a @c NULL pointer, this function will 
 * set an error code and return a @c NULL pointer. If either @em msize 
 * or @em fsize are even numbers, they are made odd by adding 1.
 * Suggested values for @em msize and @em fsize are 15 pixels for 
 * typical arc lamp spectra.
 */

cpl_image *mos_arc_background(cpl_image *image, int msize, int fsize) 
{
    const char *func = "mos_arc_background";

    cpl_image  *fimage;
    cpl_image  *bimage;
    float      *data;
    float      *bdata;
    float      *row;
    float      *brow;
    int         nx, ny;
    int         i;


    if (image == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (msize % 2 == 0)
        msize++;

    if (fsize % 2 == 0)
        fsize++;

    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);

    bimage = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);

    fimage = mos_image_filter_median(image, 3, 3);

    data = cpl_image_get_data_float(fimage);
    bdata = cpl_image_get_data_float(bimage);

    for (i = 0; i < ny; i++) {
        row = data + i * nx;
        brow = bdata + i * nx;
        if (mos_arc_background_1D(row, brow, nx, msize, fsize)) {
            cpl_error_set_where(func); 
            cpl_image_delete(fimage);
            cpl_image_delete(bimage);
            return NULL;
        }
    }

    cpl_image_delete(fimage);

    return bimage;
}


/**
 * @brief
 *   Estimate lines widths (in pixel) in arc lamp spectrum
 *
 * @param spectrum  A 1D emission line spectrum
 * @param length    Length of spectrum
 *
 * @return Mean lines width
 *
 * This function only works with emission lines spectra. The derivative
 * of the input spectrum is calculated. The result is processed once by
 * setting all its negative values to zero, and once by setting all its
 * positive values to zero and then making positive its negative values.
 * The two profiles are cross-correlated, finding in this way the mean
 * width of all the emission lines. Note that the returned width is the
 * width of the lines measured at the point of max steepness of the
 * signal (i.e., not necessarily the FWHM). In case of pure noise a
 * width 1 is returned. The maximum returned width is 20.
 */

int mos_lines_width(const float *spectrum, int length)
{

  const char *func = "mos_lines_width";

  double *profile1 = cpl_calloc(length - 1, sizeof(double));
  double *profile2 = cpl_calloc(length - 1, sizeof(double));

  double  norm, value, max;
  int     radius = 20;
  int     short_length = length - 2*radius - 1;
  int     width;
  int     i, j, k;


  /*
   * Derivative, and separation of positive and negative derivatives
   */

  for (j = 0, i = 1; i < length; j++, i++) {
      profile1[j] = profile2[j] = spectrum[i] - spectrum[j];
      if (profile1[j] < 0)
          profile1[j] = 0;
      if (profile2[j] > 0)
          profile2[j] = 0;
      else
          profile2[j] = -profile2[j];
  }


  /*
   * Profiles normalisation
   */

  length--;

  norm = 0;
  for (i = 0; i < length; i++)
      if (norm < profile1[i])
          norm = profile1[i];

  for (i = 0; i < length; i++) {
      profile1[i] /= norm;
      profile2[i] /= norm;
  }


  /*
   * Cross-correlation
   */

  max = -1;
  for (i = 0; i <= radius; i++) {
      value = 0;
      for (j = 0; j < short_length; j++) {
          k = radius+j;
          value += profile1[k] * profile2[k+i];
      }
      if (max < value) {
          max = value;
          width = i;
      }
  }

  cpl_free(profile1);
  cpl_free(profile2);

  if (max < 0.0) {
      cpl_msg_debug(func, "Cannot estimate line width");
      width = 1;
  }

  return width;

}


/**
 * @brief
 *   Find positions of peaks candidates.
 *
 * @param spectrum  A 1D emission line spectrum
 * @param length    Length of spectrum
 * @param level     Significance level
 * @param exp_width Expected lines FWHM (in pixels)
 *
 * @return List of peaks candidates positions
 *
 * A peak candidate corresponds to any pixel value above @em level
 * that is preceded and followed by a pixel with lower values. The
 * peak candidate position is determined by parabolic interpolation
 * of the three pixel values. A @c NULL pointer is returned in case
 * no peak candidate is found. No error code is set in this case.
 *
 * If very very broad and flat-topped peaks are expected, i.e.,
 * if the expected lines FWHM is more than 5 pixels, the @em spectrum
 * is slightly smoothed before peaks are searched. If the expected 
 * lines FWHM is more than 20 pixels, the spectrum is preliminary 
 * sampled at a step of half @em exp_width. These operations are
 * applied to ensure that the top of an emission line profile is 
 * never flat.
 */

cpl_vector *mos_peak_candidates(const float *spectrum, 
                                int length, float level, 
                                float exp_width)
{ 

  const char *func = "mos_peak_candidates";

  int     i, j;
  int     nint   = length - 1;
  int     n      = 0;
  int     width  = 2 * ceil(exp_width / 2) + 1;
  int     start  = width / 2;
  int     end    = length - width / 2;
  int     step;
  float  *smo;
  double *data   = cpl_calloc(length/2, sizeof(double));


  if (spectrum == NULL) {
      cpl_error_set(func, CPL_ERROR_NULL_INPUT);
      return NULL;
  }


  /*
   * If lines have a flat top (as in the case of broad slit), smooth
   * before determining the max.
   */

  if (width > 7) {
    smo = cpl_calloc(length, sizeof(float));
    start = width / 2;
    end = length - width / 2;
    for (i = 0; i < start; i++)
      smo[i] = spectrum[i];
    for (i = start; i < end; i++) {
      for (j = i - start; j <= i + start; j++)
        smo[i] += spectrum[j];
      smo[i] /= width;
    }
    for (i = end; i < length; i++)
      smo[i] = spectrum[i];
  }
  else {
      smo = (float *)spectrum;
  }

  /*
   * Collect all relative maxima along spectrum, that are higher than the
   * specified level.
   */

  if (width > 20)
    step = width / 2;
  else
    step = 1;

  for (i = step; i < nint - step + 1; i += step) {
    if (smo[i] > level) {
      if (smo[i] >= smo[i-step] && smo[i] > smo[i+step]) {
        if (smo[i-step] != 0.0 && smo[i+step] != 0.0) {
          data[n] = i + step * values_to_dx(smo[i-step], smo[i], smo[i+step]);
          ++n;
        }
      }
    }
  }

  if (width > 7) {
    cpl_free(smo);
  }

  if (n == 0) {
    cpl_free(data);
    return NULL;
  }

  return cpl_vector_wrap(n, data);

}


/**
 * @brief
 *   Improve (when possible) accuracy of peaks candidates positions.
 *
 * @param spectrum  A 1D emission line spectrum
 * @param length    Length of spectrum
 * @param peaks     List of peaks candidates
 * @param sradius   Search radius for expected peaks
 *  
 * @return Vector with refined peak positions
 *
 * The list of input peaks candidates (obtained with the function
 * @c mos_peak_candidates() ) is made more accurate, whenever possible, 
 * by applying a line baricenter determination method. In case the 
 * baricentric method fails, the corresponding peak position is not
 * rejected, but just kept unchanged. The final peaks candidates
 * list is finally cleaned from peak positions that are closer
 * than 0.5 pixel (to avoid duplications). The input vector is
 * destroyed, and a newly allocated vector is returned.
 */

cpl_vector *mos_refine_peaks(const float *spectrum, int length, 
                             cpl_vector *peaks, int sradius)
{

    const char *func = "mos_refine_peaks";

    double *data;
    float   pos;
    int     npeaks;
    int     startPos, endPos;
    int     window = 2*sradius+1;
    int     i, j;


    if (peaks == NULL || spectrum == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    npeaks = cpl_vector_get_size(peaks);
    data = cpl_vector_unwrap(peaks);

    for (i = 0; i < npeaks; i++) {
        startPos = data[i] - window/2;
        endPos   = startPos + window;
        if (startPos < 0 || endPos >= length)
            continue;

        if (0 == peakPosition(spectrum + startPos, window, &pos, 1)) {
            pos += startPos;
            data[i] = pos;
        }
    }

    for (i = 1; i < npeaks; i++)
        if (data[i] - data[i-1] < 0.5)
            data[i-1] = -1.0;

    for (i = 0, j = 0; i < npeaks; i++) {
        if (data[i] > 0.0) {
            if (i != j)
                data[j] = data[i];
            j++;
        }
    }

    return cpl_vector_wrap(j, data);

}


void mos_set_multiplex(int multiplex)
{
    mos_multiplex = multiplex;
}

/**
 * @brief
 *   Identify peak candidates
 *   
 * @param peaks     List of peaks candidates
 * @param lines     List of wavelengths
 * @param min_disp  Min expected spectral dispersion (Angstrom/pixel)
 * @param max_disp  Max expected spectral dispersion (Angstrom/pixel)
 * @param tolerance Tolerance for interval ratio comparison
 *  
 * @return List of pixel positions and wavelengths of all identified peaks
 * 
 * The list of input peaks candidates (obtained with the functions
 * @c mos_peak_candidates(), and possibly @c mos_refine_peaks() ) 
 * is compared with a list of expected emission lines wavelengths.
 * The algorithm is based on pattern recognition, where the pattern
 * is contained in the vector @em lines, and is searched in the list
 * @em peak. 
 *
 * In order to work, this method just requires a rough expectation
 * value of the spectral dispersion (in Angstrom/pixel), and a line
 * catalog. The line catalog @em lines should just include lines that 
 * are expected somewhere in the CCD exposure of the calibration lamp
 * (note, however, that a catalog including extra lines at its blue 
 * and/or red ends is still allowed).
 * 
 * Typically, the arc lamp lines candidates @em peak will include 
 * light contaminations, hot pixels, and other unwanted signal, 
 * but only in extreme cases this prevents the pattern-recognition 
 * algorithm from identifying all the spectral lines. The pattern 
 * is detected even in the case @em peak contained more arc lamp 
 * lines than actually listed in the input line catalog.
 * 
 * This method is based on the assumption that the relation between
 * wavelengths and CCD positions is with good approximation locally
 * linear (this is always true, for any existing spectrograph).
 * 
 * The ratio between consecutive intervals pairs in wavelength and in 
 * pixel is invariant to linear transformations, and therefore this 
 * quantity can be used in the recognition of local portions of the 
 * searched pattern. All the examined sub-patterns will overlap, leading 
 * to the final identification of the whole pattern, notwithstanding the 
 * overall non-linearity of the relation between pixels and wavelengths.
 *
 * Ambiguous cases, caused by exceptional regularities in the pattern,
 * or by a number of undetected (but expected) peaks that disrupt the
 * pattern on the data, are solved by linear interpolation and extrapolation
 * of the safe identifications. 
 *
 * More details about the applied algorithm can be found in the comments
 * to the function code.
 */

cpl_bivector *mos_identify_peaks(cpl_vector *peaks, cpl_vector *lines,
                                 double min_disp, double max_disp,
                                 double tolerance)
{

  int      i, j, k, l;
  int      nlint, npint;
  int      minpos;
  float    min;
  double   lratio, pratio;
  double   lo_start, lo_end, hi_start, hi_end, denom;
  double   disp, variation, prev_variation;
  int      max, maxpos, minl, mink;
  int      ambiguous;
  int      npeaks_lo, npeaks_hi;
  int     *peak_lo;
  int     *peak_hi;
  int    **ident;
  int     *nident;
  int     *lident;

  double  *peak;
  double  *line;
  int      npeaks, nlines;

  double  *xpos;
  double  *lambda;
  int     *ilambda;
  double  *tmp_xpos;
  double  *tmp_lambda;
  int     *tmp_ilambda;
  int     *flag;
  int      n = 0;
  int      nn;
  int      nseq = 0;
  int      gap;
  int     *seq_length;
  int      found;

  peak        = cpl_vector_get_data(peaks);
  npeaks      = cpl_vector_get_size(peaks);
  line        = cpl_vector_get_data(lines);
  nlines      = cpl_vector_get_size(lines);

  if (npeaks < 4)
      return NULL;

  peak_lo     = cpl_malloc(npeaks * sizeof(int));
  peak_hi     = cpl_malloc(npeaks * sizeof(int));
  nident      = cpl_calloc(npeaks, sizeof(int));
  lident      = cpl_calloc(nlines, sizeof(int));
  xpos        = cpl_calloc(npeaks, sizeof(double));
  lambda      = cpl_calloc(npeaks, sizeof(double));
  ilambda     = cpl_calloc(npeaks, sizeof(int));
  tmp_xpos    = cpl_calloc(npeaks, sizeof(double));
  tmp_lambda  = cpl_calloc(npeaks, sizeof(double));
  tmp_ilambda = cpl_calloc(npeaks, sizeof(int));
  flag        = cpl_calloc(npeaks, sizeof(int));
  seq_length  = cpl_calloc(npeaks, sizeof(int));
  ident       = cpl_malloc(npeaks * sizeof(int *));
  //In the very worst case scenario, a given peak is always assigned
  //a possible identification (within the tolerance). This is actually
  //very unlikely but not impossible.
  for (i = 0; i < npeaks; i++)
    ident[i]  = cpl_malloc(npeaks * npeaks * sizeof(int));

  /*
   * This is just the number of intervals - one less than the number
   * of points (catalog wavelengths, or detected peaks).
   */

  nlint = nlines - 1;
  npint = npeaks - 1;


  /*
   * Here the big loops on catalog lines begins.
   */

  for (i = 1; i < nlint; i++) {


    /*
     * For each catalog wavelength I take the previous and the next one, 
     * and compute the ratio of the corresponding wavelength intervals.
     * This ratio will be compared to all the ratios obtained doing the
     * same with all the detected peaks positions.
     */

    lratio = (line[i+1] - line[i]) / (line[i] - line[i-1]);


    /*
     * Here the loop on detected peaks positions begins.
     */

    for (j = 1; j < npint; j++) {

      /*
       * Not all peaks are used for computing ratios: just the ones
       * that are compatible with the expected spectral dispersion
       * are taken into consideration. Therefore, I define the pixel
       * intervals before and after any peak that are compatible with
       * the specified dispersion interval, and select just the peaks
       * within such intervals. If either of the two intervals doesn't
       * contain any peak, then I skip the current peak and continue
       * with the next.
       */

      lo_start = peak[j] - (line[i] - line[i-1]) / min_disp;
      lo_end   = peak[j] - (line[i] - line[i-1]) / max_disp;
      hi_start = peak[j] + (line[i+1] - line[i]) / max_disp;
      hi_end   = peak[j] + (line[i+1] - line[i]) / min_disp;

      for (npeaks_lo = 0, k = 0; k < npeaks; k++) {
        if (peak[k] > lo_end)
          break;
        if (peak[k] > lo_start) {
          peak_lo[npeaks_lo] = k;
          ++npeaks_lo;
        }
      }

      if (npeaks_lo == 0)
        continue;

      for (npeaks_hi = 0, k = 0; k < npeaks; k++) {
        if (peak[k] > hi_end)
          break;
        if (peak[k] > hi_start) {
          peak_hi[npeaks_hi] = k;
          ++npeaks_hi;
        }
      }

      if (npeaks_hi == 0)
        continue;


      /*
       * Now I have all peaks that may help for a local identification.
       * peak_lo[k] is the sequence number of the k-th peak of the lower
       * interval; peak_hi[l] is the sequence number of the l-th peak of
       * the higher interval. j is, of course, the sequence number of the
       * current peak (second big loop).
       */

      prev_variation = 1000.0;
      minl = mink = 0;

      for (k = 0; k < npeaks_lo; k++) {
        denom = peak[j] - peak[peak_lo[k]];
        for (l = 0; l < npeaks_hi; l++) {

          /*
           * For any pair of peaks - one from the lower and the other
           * from the higher interval - I compute the same ratio that
           * was computed with the current line catalog wavelength.
           */

          pratio = (peak[peak_hi[l]] - peak[j]) / denom;

          /*
           * If the two ratios are compatible within the specified
           * tolerance, we have a preliminary identification. This
           * will be marked in the matrix ident[][], where the first
           * index corresponds to a peak sequence number, and the second
           * index is the counter of the identifications made during
           * this whole process. The array of counters is nident[].
           * If more than one interval pair fulfills the specified
           * tolerance, the closest to the expected ratio is selected.
           */

          variation = fabs(lratio-pratio) / pratio;

          if (variation < tolerance) {
            if (variation < prev_variation) {
              prev_variation = variation;
              minl = l;
              mink = k;
            }
          }
        }
      }
      if (prev_variation < tolerance) {
        ident[j][nident[j]]                         = i;
        ident[peak_hi[minl]][nident[peak_hi[minl]]] = i + 1;
        ident[peak_lo[mink]][nident[peak_lo[mink]]] = i - 1;
        ++nident[j];
        ++nident[peak_hi[minl]];
        ++nident[peak_lo[mink]];
      }
    }   /* End loop on positions */
  }    /* End loop on lines     */


  /*
   * At this point I have filled the ident matrix with all my preliminary
   * identifications. Ambiguous identifications must be eliminated.
   */


  for (i = 0; i < npeaks; i++) {


    /*
     * I don't take into consideration peaks that were never identified.
     * They are likely contaminations, or emission lines that were not
     * listed in the input wavelength catalog.
     */

    if (nident[i] > 1) {


      /*
       * Initialise the histogram of wavelengths assigned to the i-th peak.
       */

      for (j = 0; j < nlines; j++)
        lident[j] = 0;


      /*
       * Count how many times each catalog wavelength was assigned
       * to the i-th peak.
       */

      for (j = 0; j < nident[i]; j++)
        ++lident[ident[i][j]];


      /*
       * What wavelength was most frequently assigned to the i-th peak?
       */

      max = 0;
      maxpos = 0;
      for (j = 0; j < nlines; j++) {
        if (max < lident[j]) {
          max = lident[j];
          maxpos = j;
        }
      }


      /*
       * Were there other wavelengths assigned with the same frequency?
       * This would be the case of an ambiguous identification. It is
       * safer to reject this peak...
       */

      ambiguous = 0;

      for (k = maxpos + 1; k < nlines; k++) {
        if (lident[k] == max) {
          ambiguous = 1;
          break;
        }
      }

      if (ambiguous)
        continue;


      /*
       * Otherwise, I assign to the i-th peak the wavelength that was
       * most often assigned to it.
       */

      tmp_xpos[n]   = peak[i];
      tmp_lambda[n] = line[maxpos];
      tmp_ilambda[n] = maxpos;

      ++n;

    }

  }


  /*
   * Check on identified peaks. Contaminations from other spectra might 
   * be present and should be excluded: this type of contamination 
   * consists of peaks that have been _correctly_ identified! The non-
   * spectral type of light contamination should have been almost all 
   * removed already in the previous steps, but it may still be present.
   * Here, the self-consistent sequences of identified peaks are
   * separated one from the other. At the moment, just the longest of
   * such sequences is selected (in other words, spectral multiplexing
   * is ignored).
   */

  if (n > 1) {
    nn = 0;                  /* Number of peaks in the list of sequences */
    nseq = 0;                /* Current sequence */
    for (k = 0; k < n; k++) {
      if (flag[k] == 0) {    /* Was peak k already assigned to a sequence? */
        flag[k] = 1;
        xpos[nn] = tmp_xpos[k];       /* Begin the nseq-th sequence */
        lambda[nn] = tmp_lambda[k];
        ilambda[nn] = tmp_ilambda[k];
        ++seq_length[nseq];
        ++nn;

        /*
         * Now look for all the following peaks that are compatible
         * with the expected spectral dispersion, and add them in 
         * sequence to xpos. Note that missing peaks are not a problem...
         */
         
        i = k;
        while (i < n - 1) {
          found = 0;
          for (j = i + 1; j < n; j++) {
            if (flag[j] == 0) {
              disp = (tmp_lambda[j] - tmp_lambda[i])
                   / (tmp_xpos[j] - tmp_xpos[i]);
              if (disp >= min_disp && disp <= max_disp) {
                flag[j] = 1;
                xpos[nn] = tmp_xpos[j];
                lambda[nn] = tmp_lambda[j];
                ilambda[nn] = tmp_ilambda[j];
                ++seq_length[nseq];
                ++nn;
                i = j;
                found = 1;
                break;
              }
            }
          }
          if (!found)
            break;
        }

        /*
         * Current sequence is completed: begin new sequence on the
         * excluded peaks, starting the loop on peaks again.
         */

        ++nseq;
        k = 0;
      }
    }


    /*
     * Find the longest sequence of self-consistent peaks.
     */

    maxpos = max = 0;

    if (mos_multiplex < 0) {
      for (i = 0; i < nseq; i++) {
        if (seq_length[i] > max) {
          max = seq_length[i];
          maxpos = i;
        }
      }
    }
    else {

      /*
       * Now consider the sequence which lays in the specified 
       * CCD region (indicated by mos_multiplex): that is, _most_ 
       * of its lines (more than half) must be in that region...
       */

      nn = 0;
      found = 0;

      for (i = 0; i < nseq; i++) {
        n = seq_length[i];
        if (n > 5) {
          cpl_array *regions = cpl_array_new(n, CPL_TYPE_INT);
          int        region;

          for (j = 0; j < n; j++)
            cpl_array_set_int(regions, j, 
                              ((int)floor(xpos[nn + j])) / mos_region_size);

          region = (int)cpl_array_get_median(regions);
          cpl_array_delete(regions);

          if (mos_multiplex == region) {
            if (found) {
              cpl_msg_debug(cpl_func, "More than one spectrum found in "
                            "region %d (only the first one is extracted)", 
                            mos_multiplex);
              break;
            }
            found = 1;
            max = seq_length[i];
            maxpos = i;
          }
        }
        nn += seq_length[i];
      }
    }

    /*
     * Find where this sequence starts in the whole peak position
     * storage.
     */

    nn = 0;
    for (i = 0; i < maxpos; i++)
      nn += seq_length[i];

    /*
     * Move the longest sequence at the beginning of the returned lists
     */

    n = max;
    for (i = 0; i < n; i++, nn++) {
      xpos[i] = xpos[nn];
      lambda[i] = lambda[nn];
      ilambda[i] = ilambda[nn];
    }


    /*
     * Are some wavelengths missing? Recover them.
     */

    for (i = 1; i < n; i++) {
      gap = ilambda[i] - ilambda[i-1];
      for (j = 1; j < gap; j++) {

        if (j == 1) {

          /*
           * Determine the local dispersion from the current pair of peaks
           */
  
          disp = (lambda[i] - lambda[i-1]) / (xpos[i] - xpos[i-1]);
        }

        /*
         * With this, find the expected position of the missing
         * peak by linear interpolation.
         */

        hi_start = xpos[i-1] + (line[ilambda[i-1] + j] - lambda[i-1]) / disp;

        /*
         * Is there a peak at that position? Here a peak from the
         * original list is searched, that is closer than 2 pixels
         * to the expected position. If it is found, insert it at
         * the current position on the list of identified peaks,
         * and leave immediately the loop (taking the new position
         * for the following linear interpolation, in case more
         * than one peak is missing in the current interval).
         * If it is not found, stay in the loop, looking for 
         * the following missing peaks in this interval.
         */

        found = 0;
        for (k = 0; k < npeaks; k++) {
          if (fabs(peak[k] - hi_start) < 2) {
            for (l = n; l > i; l--) {
              xpos[l] = xpos[l-1];
              lambda[l] = lambda[l-1];
              ilambda[l] = ilambda[l-1];
            }
            xpos[i] = peak[k];
            lambda[i] = line[ilambda[i-1] + j];
            ilambda[i] = ilambda[i-1] + j;
            ++n;
            found = 1;
            break;
          }
        }
        if (found)
          break;
      }
    }


    /*
     * Try to extrapolate forward
     */

    found = 1;

    if (n > 0) {

        while (ilambda[n-1] < nlines - 1 && found) {

          /*
           * Determine the local dispersion from the last pair of 
           * identified peaks
           */

          if (n > 1)
              disp = (lambda[n-1] - lambda[n-2]) / (xpos[n-1] - xpos[n-2]);
          else
              disp = 0.0;

          if (disp > max_disp || disp < min_disp)
            break;


          /*
           * With this, find the expected position of the missing
           * peak by linear interpolation.
           */

          hi_start = xpos[n-1] + (line[ilambda[n-1] + 1] - lambda[n-1]) / disp;

          /*
           * Is there a peak at that position? Here a peak from the
           * original list is searched, that is closer than 6 pixels
           * to the expected position. If it is found, insert it at
           * the end of the list of identified peaks. If it is not
           * found, leave the loop.
           */

          found = 0;
          min = fabs(peak[0] - hi_start);
          minpos = 0;
          for (k = 1; k < npeaks; k++) {
            if (min > fabs(peak[k] - hi_start)) {
                min = fabs(peak[k] - hi_start);
                minpos = k;
            }
          }
          if (min < 6 && fabs(peak[minpos] - xpos[n-1]) > 1.0) {
            xpos[n] = peak[minpos];
            lambda[n] = line[ilambda[n-1] + 1];
            ilambda[n] = ilambda[n-1] + 1;
            ++n;
            found = 1;
          }
        }
    }


    /*
     * Try to extrapolate backward
     */

    found = 1;
    while (ilambda[0] > 0 && found) {

      /*
       * Determine the local dispersion from the first pair of
       * identified peaks
       */

      disp = (lambda[1] - lambda[0]) / (xpos[1] - xpos[0]);

      if (disp > max_disp || disp < min_disp)
        break;


      /*
       * With this, find the expected position of the missing
       * peak by linear interpolation.
       */

      hi_start = xpos[0] - (lambda[0] - line[ilambda[0] - 1]) / disp;


      /*
       * Is there a peak at that position? Here a peak from the
       * original list is searched, that is closer than 6 pixels
       * to the expected position. If it is found, insert it at
       * the beginning of the list of identified peaks. If it is not
       * found, leave the loop.
       */

      found = 0;
      min = fabs(peak[0] - hi_start);
      minpos = 0;
      for (k = 1; k < npeaks; k++) {
        if (min > fabs(peak[k] - hi_start)) {
            min = fabs(peak[k] - hi_start);
            minpos = k;
        }
      }
      if (min < 6 && fabs(peak[minpos] - xpos[0]) > 1.0) {
        for (j = n; j > 0; j--) {
          xpos[j] = xpos[j-1];
          lambda[j] = lambda[j-1];
          ilambda[j] = ilambda[j-1];
        }
        xpos[0] = peak[minpos];
        lambda[0] = line[ilambda[0] - 1];
        ilambda[0] = ilambda[0] - 1;
        ++n;
        found = 1;
      }
    }
  }


  /*
   * At this point all peaks are processed. Free memory, and return
   * the result.
   */

/************************************************+
  for (i = 0; i < npeaks; i++) {
    printf("Peak %d:\n   ", i);
    for (j = 0; j < nident[i]; j++)
      printf("%.2f, ", line[ident[i][j]]);
    printf("\n");
  }

  printf("\n");

  for (i = 0; i < n; i++)
    printf("%.2f, %.2f\n", xpos[i], lambda[i]);
+************************************************/
  for (i = 0; i < npeaks; i++)
    cpl_free(ident[i]);
  cpl_free(ident);
  cpl_free(nident);
  cpl_free(lident);
  cpl_free(ilambda);
  cpl_free(tmp_xpos);
  cpl_free(tmp_lambda);
  cpl_free(tmp_ilambda);
  cpl_free(peak_lo);
  cpl_free(flag);
  cpl_free(seq_length);
  cpl_free(peak_hi);

  if (n == 0) {
    cpl_free(xpos);
    cpl_free(lambda);
    return NULL;
  }

  return cpl_bivector_wrap_vectors(cpl_vector_wrap(n, xpos), 
                                   cpl_vector_wrap(n, lambda));
}

/* 
 * This is an attempt to interface ppm to the new ppm function
 * of the CPL. Very slow, very inaccurate, probably the radius
 * should be decreased to max 4 pixel.
 */

cpl_bivector *mos_identify_peaks_new(cpl_vector *peaks, cpl_vector *lines,
                                     double min_disp, double max_disp,
                                     double tolerance)
{
    cpl_matrix   *data;
    cpl_matrix   *pattern;
    cpl_matrix   *mdata;
    cpl_matrix   *mpattern;
    cpl_matrix   *subdata;
    int           use_data    = cpl_vector_get_size(peaks);
    int           use_pattern = cpl_vector_get_size(lines);
    const double  err_data    = 0.3;
    const double  err_pattern = 0.01;
    const double  radius      = 20;

    cpl_vector   *x;
    cpl_vector   *l;
    cpl_bivector *xl;

    //Avoid compiler warning
    (void)max_disp;
    (void)min_disp;

    data    = cpl_matrix_new(2, use_data);
    subdata = cpl_matrix_wrap(1, use_data, cpl_vector_get_data(peaks));
    cpl_matrix_copy(data, subdata, 0, 0);
    cpl_matrix_unwrap(subdata);

    pattern = cpl_matrix_new(2, use_pattern);
    subdata = cpl_matrix_wrap(1, use_pattern, cpl_vector_get_data(lines));
    cpl_matrix_copy(pattern, subdata, 0, 0);
    cpl_matrix_unwrap(subdata);

printf("input data:\n");
cpl_matrix_dump(data, NULL);
printf("input pattern:\n");
cpl_matrix_dump(pattern, NULL);

    cpl_array_delete(
                     cpl_ppm_match_points(data, use_data, err_data, pattern, 
                                          use_pattern, err_pattern, tolerance, 
                                          radius, &mdata, &mpattern, 
                                          NULL, NULL)
                    );

    cpl_matrix_delete(data);
    cpl_matrix_delete(pattern);

    if (mdata == NULL)
        return NULL;

    cpl_matrix_sort_columns(mdata, 1);
    cpl_matrix_sort_columns(mpattern, 1);

printf("RISULTATO:\n");
printf("data:\n");
cpl_matrix_dump(mdata, NULL);
printf("pattern:\n");
cpl_matrix_dump(mpattern, NULL);

    use_data = cpl_matrix_get_ncol(mdata);

    x = cpl_vector_wrap(use_data, cpl_matrix_get_data(mdata));
    l = cpl_vector_wrap(use_data, cpl_matrix_get_data(mpattern));

    xl = cpl_bivector_wrap_vectors(cpl_vector_duplicate(x),
                                   cpl_vector_duplicate(l));

    cpl_vector_unwrap(x);
    cpl_vector_unwrap(l);
    cpl_matrix_delete(mdata);
    cpl_matrix_delete(mpattern);

    return xl;

}


/**
 * @brief
 *   Evaluate the wavelength of a pixel position
 *
 * @param ids      Inverse dispersion relation (from wave to pixel)
 * @param blue     Start wavelength of ids validity
 * @param red      End wavelength of ids validity
 * @param refwave  Reference wavelength
 * @param pixel    Pixel position
 *
 * @return Wavelength of pixel
 *
 * The accuracy of the returned wavelength is guaranteed to be better 
 * than 0.02 pixels (converted from wavelength units to pixel). If @em pixel
 * is outside the @em ids validity range, the wavelength 0.0 is returned.
 */

/*
double mos_eval_dds(cpl_polynomial *ids, double blue, double red, 
                    double refwave, double pixel)
{
    double yellow;
    double cpixel;
    double tolerance = 0.02;
    int    max_iter = 20;
    int    iter = 0;

    
    if (cpl_polynomial_eval_1d(ids, blue-refwave, NULL) > pixel)
        return 0.0;
    
    if (cpl_polynomial_eval_1d(ids, red-refwave, NULL) < pixel)
        return 0.0;

    yellow = (blue + red) / 2;

    cpixel = cpl_polynomial_eval_1d(ids, yellow-refwave, NULL);

    while (fabs(cpixel - pixel) > tolerance && iter < max_iter) {

        if (cpixel > pixel)
            red = yellow;
        else
            blue = yellow;

        yellow = (blue + red) / 2;
        cpixel = cpl_polynomial_eval_1d(ids, yellow-refwave, NULL);

        iter++;

    }

    return yellow;

}
*/

double mos_eval_dds(cpl_polynomial *ids, double blue, double red,
                    double refwave, double pixel)
{
    double     yellow;
    double     coeff;
    cpl_size   zero = 0;

    if (cpl_polynomial_eval_1d(ids, blue-refwave, NULL) > pixel)
        return 0.0;

    if (cpl_polynomial_eval_1d(ids, red-refwave, NULL) < pixel)
        return 0.0;

    yellow = (blue + red) / 2 - refwave;

    coeff = cpl_polynomial_get_coeff(ids, &zero);
    cpl_polynomial_set_coeff(ids, &zero, coeff - pixel);

    cpl_polynomial_solve_1d(ids, yellow, &yellow, 1);
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_error_reset();
        return 0.0;
    }

    cpl_polynomial_set_coeff(ids, &zero, coeff);

    return yellow + refwave;

}

/**
 * @brief
 *   Fit polynomial relation from wavelengths to pixels
 *
 * @param pixwav   List of pixel positions and associated wavelengths
 * @param order    Order of the fitting polynomial
 * @param reject   Max residual tolerated for line rejection (in pixel)
 * @param minlines Min number of lines to be used in the fit
 * @param nlines   Output number of lines actually used in the fit
 * @param err      Output RMS of the fit residuals
 *
 * @return Polynomial relation from wavelengths to pixels
 *
 * The list of identified peaks and their wavelengths (obtained either with 
 * the function @c mos_identify_peaks() or the function @c mos_find_peaks() ) 
 * must contain at least @em minlines entries. A fit is tried, and all the 
 * peaks farther than the @em reject threshold from the model are rejected. 
 * This process is iterated until either the number of surviving lines is 
 * less than @em minlines or there are no more outliers (whichever comes 
 * first). Only in the latter case a fit is returned. If the @em reject 
 * threshold is negative, no outliers rejection is applied, and the first 
 * fit is accepted. In case of error, a @c NULL pointer is returned, and 
 * @em nlines and @em err are set to zero.
 */  

cpl_polynomial *mos_poly_wav2pix(cpl_bivector *pixwav, int order, 
                                 double reject, int minlines, 
                                 int *nlines, double *err,
                                 cpl_bivector **pixwav_used)
{
    const char   *func = "mos_poly_wav2pix";

    cpl_bivector *pixwav2;
    cpl_vector   *wavel;
    cpl_vector   *pixel;
    double       *d_wavel;
    double       *d_pixel;
    double        pixpos;
    int           fitlines;
    int           rejection = 0;
    int           i, j;

    cpl_polynomial *ids;


    *nlines = 0;
    *err = 0;

    if (pixwav == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    fitlines = cpl_bivector_get_size(pixwav);

    if (fitlines < minlines) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }


    /*
     * If outliers rejection was requested, allocate a working
     * vector (that can be modified as soon as outliers are removed)
     */

    if (reject > 0.0)
        rejection = 1;

    if (rejection)
        pixwav2 = cpl_bivector_duplicate(pixwav);
    else
        pixwav2 = pixwav;


    /*
     * The single vectors are extracted just because the fitting routine
     * requires it
     */

    pixel = cpl_bivector_get_x(pixwav2);
    wavel = cpl_bivector_get_y(pixwav2);


    /*
     * Get rid of the wrapper, in case of duplication
     */

    if (rejection)
        cpl_bivector_unwrap_vectors(pixwav2);


    /*
     * Here begins the iterative fit of identified lines
     */

    while (fitlines >= minlines) {

        ids = cpl_polynomial_fit_1d_create(wavel, pixel, order, err);
        *err = sqrt(*err);
    
        if (ids == NULL) {
            cpl_msg_debug(cpl_error_get_where(), "%s", cpl_error_get_message());
            cpl_msg_debug(func, "Fitting IDS");
            cpl_error_set_where(func);
            if (rejection) {
                cpl_vector_delete(wavel);
                cpl_vector_delete(pixel);
            }
            return NULL;
        }

        if (rejection) {
            cpl_vector * wavel_used = cpl_vector_duplicate(wavel);
            cpl_vector * pixel_used = cpl_vector_duplicate(pixel);


            /*
             * Now work directly with the vector data buffers...
             */

            d_pixel = cpl_vector_unwrap(pixel);
            d_wavel = cpl_vector_unwrap(wavel);

            for (i = 0, j = 0; i < fitlines; i++) {
                pixpos = cpl_polynomial_eval_1d(ids, d_wavel[i], NULL);
                if (fabs(pixpos - d_pixel[i]) < reject) {
                    d_pixel[j] = d_pixel[i];
                    d_wavel[j] = d_wavel[i];
                    j++;
                }
            }
    
            if (j == fitlines) {       /* No rejection in last iteration */
                cpl_bivector * pixwav_used_temp = 
                        cpl_bivector_wrap_vectors(pixel_used, wavel_used);
                *pixwav_used = cpl_bivector_duplicate(pixwav_used_temp);
                cpl_bivector_unwrap_vectors(pixwav_used_temp);
                cpl_vector_delete(wavel_used);
                cpl_vector_delete(pixel_used);
                cpl_free(d_wavel);
                cpl_free(d_pixel);
                *nlines = fitlines;
                return ids;
            }
            else {                     /* Some lines were rejected       */
                fitlines = j;
                cpl_polynomial_delete(ids);
                if (fitlines >= minlines) {
                    pixel = cpl_vector_wrap(fitlines, d_pixel);
                    wavel = cpl_vector_wrap(fitlines, d_wavel);
                }
                else {                 /* Too few lines: failure         */
                    cpl_free(d_wavel);
                    cpl_free(d_pixel);
                    cpl_error_set(func, CPL_ERROR_CONTINUE);
                    return NULL;
                }
            }
            cpl_vector_delete(wavel_used);
            cpl_vector_delete(pixel_used);
        }
        else {
            *nlines = fitlines;
            *pixwav_used = cpl_bivector_duplicate(pixwav2);
            return ids;       /* Exit at first iteration if no rejection */
        }
    }
    
    return ids;               /* To avoid compiler warnings */
}


/**
 * @brief
 *   Fit polynomial relation from pixels to wavelengths
 *
 * @param pixwav   List of pixel positions and associated wavelengths
 * @param order    Order of the fitting polynomial
 * @param reject   Max residual tolerated for line rejection (in wave units)
 * @param minlines Min number of lines to be used in the fit
 * @param nlines   Output number of lines actually used in the fit
 * @param err      Output RMS of the fit residuals
 *
 * @return Polynomial relation from wavelengths to pixels
 *
 * The list of identified peaks and their wavelengths (obtained with the
 * function @c mos_identify_peaks() ) must contain at least @em minlines
 * entries. A fit is tried, and all the peaks farther than the @em reject
 * threshold from the model are rejected. This process is iterated until
 * either the number of surviving lines is less than @em minlines or there
 * are no more outliers (whichever comes first). Only in the latter case
 * a fit is returned. If the @em reject threshold is negative, no outliers
 * rejection is applied, and the first fit is accepted. In case of error,
 * a @c NULL pointer is returned, and @em nlines and @em err are set to zero.
 */

cpl_polynomial *mos_poly_pix2wav(cpl_bivector *pixwav, int order,
                                 double reject, int minlines, 
                                 int *nlines, double *err)
{

    cpl_bivector *wavpix;
    cpl_vector   *wavel;
    cpl_vector   *pixel;

    cpl_polynomial *dds;
    
    cpl_bivector *wavepix_used;


    /*
     * Swap vectors in bivector, in order to reuse mos_poly_wav2pix()
     */

    pixel = cpl_bivector_get_x(pixwav);
    wavel = cpl_bivector_get_y(pixwav);

    wavpix = cpl_bivector_wrap_vectors(wavel, pixel);

    dds = mos_poly_wav2pix(wavpix, order, reject, minlines, nlines, err,
                           &wavepix_used);

    cpl_bivector_unwrap_vectors(wavpix);
    
    cpl_bivector_delete(wavepix_used);

    return dds;

}


/**
 * @brief
 *   Find the reference lines peaks using a polynomial first-guess
 *
 * @param spectrum  A 1D emission line spectrum
 * @param length    Length of spectrum
 * @param lines     List of wavelengths
 * @param ids       Polynomial conversion from wavelengths to pixel
 * @param refwave   Zero wavelength used in ids determination
 * @param sradius   Search radius for expected peaks
 *
 * @return List of pixel positions and wavelengths of all identified peaks
 *
 * The input polynomial @em ids is applied to the input wavelengths
 * to find the expected position of the corresponding peak along the
 * input @em spectrum. The expected peak is searched within a window
 * of radius @em sradius. A list is returned, with the positions of
 * the detected peaks with their associated wavelengths. The @em sradius
 * must be at least 1 pixel, and the input spectrum must be at least 
 * twice + 1 @em sradius. In case of error, a @c NULL pointer is returned.
 */

cpl_bivector *mos_find_peaks(const float *spectrum, int length, 
                             cpl_vector *lines, cpl_polynomial *ids, 
                             double refwave, int sradius)
{
    const char   *func = "mos_find_peaks";

    double       *data;
    double       *d_pixel;
    double       *d_wavel;
    float         pos;
    int           nlines;
    int           pixel;
    int           i, j;


    if (spectrum == NULL || lines == NULL || ids == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nlines = cpl_vector_get_size(lines);

    if (sradius < 1 || length < 2*sradius+1 || nlines < 1) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    d_wavel = cpl_malloc(nlines * sizeof(double));
    d_pixel = cpl_malloc(nlines * sizeof(double));

    data = cpl_vector_get_data(lines);

    for (i = 0, j = 0; i < nlines; i++) {
        pixel = cpl_polynomial_eval_1d(ids, data[i]-refwave, NULL) + 0.5;
        if (pixel < 0 || pixel - sradius < 0 || pixel + sradius >= length)
            continue;
        if (peakPosition(spectrum+pixel-sradius, 2*sradius+1, &pos, 1) == 0) {
            pos += pixel - sradius;
            d_pixel[j] = pos;
            d_wavel[j] = data[i];
            j++;
        }
    }

    if (j > 0) {
        return cpl_bivector_wrap_vectors(cpl_vector_wrap(j, d_pixel),
                                         cpl_vector_wrap(j, d_wavel));
    }
    else {
        cpl_free(d_wavel);
        cpl_free(d_pixel);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_OUTPUT);
        return NULL;
    }
}


/**
 * @brief
 *   Derive wavelength calibration from a raw arc lamp or sky exposure
 *
 * @param image       An arc lamp or sky exposure
 * @param lines       List of reference lines wavelengths
 * @param dispersion  Expected value of the dispersion (wavelength units/pixel)
 * @param level       Threshold for peak detection
 * @param sradius     Search radius for expected peaks (pixels)
 * @param order       Degree of fitting polynomial for wavelength calibration
 * @param reject      Max residual tolerated for line rejection (pixels)
 * @param refwave     Reference wavelength
 * @param wavestart   I/O wavelength of first pixel of resampled image
 * @param waveend     I/O wavelength of last pixel of resampled image
 * @param nlines      Returned array of number of lines used for each fit
 * @param error       Returned array of mean accuracies obtained for each fit
 * @param idscoeff    Returned table with IDS polynomials
 * @param calibration Returned wavelength calibration image
 * @param residuals   Returned residuals image
 * @param restable    Returned residuals table
 * @param refmask     Returned monochromatic image at reference wavelength
 *
 * @return Input exposure resampled at constant wavelength step.
 *
 * The input @em image must contain spectra with the dispersion direction
 * approximately horizontal, with blue on the left and red on the right.
 * 
 * The rows of the input @em image are independently calibrated one by 
 * one. First, the function @c mos_arc_background() is used to subtract 
 * the continuum from the input spectrum. Second, the two functions 
 * @c mos_peak_candidates() and @c mos_refine_peaks() are used to create 
 * a list of positions of reference lines candidates along each @em image 
 * row. Only peaks above @em level are selected. Third, the function 
 * @c mos_identify_peaks() is applied to select from the found peaks 
 * the ones corresponding to reference lines, associating them to the 
 * appropriate wavelengths from the line catalog @em lines. The ensuing 
 * polynomial fit is adaptive: it is performed using the specified 
 * polynomial @em order in case the resulting list of positions in 
 * pixels and wavelengths contains at least twice the degrees of freedom 
 * of the polynomial (i.e., @em order + 1). If this is not the case, 
 * the order of the fitting polynomial is adapted to the number of 
 * actual entries. The number of points to fit, however, should never 
 * be less than 4. The fit is performed both from wavelength to pixel 
 * and from pixel to wavelength, using the functions @c mos_poly_wav2pix() 
 * and @c mos_poly_pix2wav(). If a @em reject threshold greater than 
 * zero is specified, outlaying entries are removed iteratively from 
 * the list. If the @em reject threshold is negative, no outliers 
 * rejection is applied, and the first fit is accepted.
 *
 * Optionally, in case @em sradius is positive, the obtained polynomial 
 * solution is passed to the function @c mos_find_peaks(), that searches 
 * again along each image row the reference lines candidates 
 * around their expected positions, within the specified search radius;
 * The polynomial fitting is then repeated with the new found positions. 
 * This option can be useful for recovering very faint (i.e., below 
 * @em level) reference lines, or reference lines that were lost by
 * @c mos_identify_peaks() because of a partially wrong input @em lines
 * list. 
 *
 * An array @em nlines, containing the number of lines used for each 
 * fit, and an array @em error, containing the mean error of the 
 * polynomial models (in pixels), are returned. A fit failure is 
 * indicated with the corresponding element of @em nlines set to 
 * zero. The mean error of the polynomial model is evaluated by 
 * dividing the RMS of the fit residuals by the square root of the 
 * number of fitted points divided the degrees of freedom of the model:
 *
 *               mean error = RMS / sqrt(N / (@em order + 1))
 * 
 * The arrays @em nlines and @em error must be pre-allocated, and 
 * should all have as many elements as the number of rows in the input 
 * @em image. If @c NULL pointers are passed, they are not computed.
 * In the same way a preallocated @em idscoeff table may or may not
 * be passed to this function, and it will be filled with the dispersion
 * relation coefficients only in the former case. This table should be
 * preallocated with the same number of rows as the input @em image. 
 * No columns should be defined in this table: they will be created 
 * automatically by this function, and will be labeled c0, c1, c2, ... 
 * up to the specified @em order of the fitting polynomial.
 *
 * As a by-product of the wavelength calibration, the input @em image 
 * is resampled at a constant wavelength step, @em dispersion, and is 
 * returned by this function. In case of error a @c NULL pointer is 
 * returned. If the input arguments @em wavestart and @em waveend
 * are greater than 1.0, they are taken as the spectral interval
 * where the spectra are resampled; alternatively, the wavelength
 * range covered by the resampled image is equal to the wavelength
 * range of the input reference @em lines catalog, extended by 10
 * percent on the blue and the red sides: such used spectral interval
 * is then returned via the same variables, @em wavestart and @em waveend.
 * Note that limiting the spectral range doesn't prevents this function
 * to perform a wavelength calibration based on all the reference 
 * wavelengths listed in the input line catalog, including those 
 * outside the specified range (if they are found somewhere on the
 * detector).
 *
 * Optionally, an image of the wavelength calibrated input exposure, 
 * @em calibration, an image of the fit residuals, @em residuals,
 * and a monochromatic mask image, @em refmask, that is obtained at a 
 * given wavelength @em refwave, can be returned. These images must be 
 * pre-allocated, and should all have the same size of the input 
 * @em image. If @c NULL pointers are passed, they are not computed. 
 * The @em calibration image consists of pixels having the value of 
 * their wavelength, or the value zero if this is not available. The 
 * @em residuals image has, at the positions of the reference lines,
 * the value of the corresponding distance from the polynomial model.
 * The @em refmask image is used to flag the pixels containing the 
 * specified reference wavelength. All the other pixels are set to 
 * zero. This mask is the monochromatic image on the CCD of the slits 
 * located on the telescope focal plane, and can be used in the
 * determination of an optical distortion model. In order to clean
 * this image from occasional bad fits contributions and fill possible 
 * small gaps in the wavelength calibration, a morphological closing
 * (dilation + erosion), followed by a morphological opening 
 * (erosion + dilation), is applied.
 *
 * If detected_lines table is an input, then a table with the position of
 * the detected lines, its peak flux, wavelength identified and final xpos
 * position in the wavelength calibrated image using the wavelength solution
 *  is returned. The table has to be allocated but empty.
 */

cpl_image *mos_wavelength_calibration_raw(const cpl_image *image,
                                          cpl_vector *lines,
                                          double dispersion, float level,
                                          int sradius, int order,
                                          double reject, double refwave, 
                                          double *wavestart, double *waveend,
                                          int *nlines, double *error, 
                                          cpl_table *idscoeff,
                                          cpl_image *calibration,
                                          cpl_image *residuals, 
                                          cpl_table *restable,
                                          cpl_mask *refmask,
                                          cpl_table *detected_lines, 
                                          double disp_tolerance,
                                          double ratio_tolerance)
{

    const char *func = "mos_wavelength_calibration_raw";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    double  tolerance = disp_tolerance;     
    double  r_toleran = ratio_tolerance;     
    int     step      = 10;       /* Compute restable every "step" rows */

    char            name[MAX_COLNAME];
    cpl_image      *resampled;
    cpl_bivector   *output;
    cpl_bivector   *new_output;
    cpl_vector     *peaks;
    cpl_vector     *wavel;
    cpl_polynomial *ids;
    cpl_polynomial *lin;
    cpl_matrix     *kernel;
    double          ids_err;
    double          max_disp, min_disp;
    double         *line;
    double          firstLambda, lastLambda, lambda;
    double          value, wave, pixe;
    cpl_binary     *mdata;
    const float    *sdata;
    float          *rdata;
    float          *idata;
    float          *ddata;
    float           v1, v2, vi;
    float           fpixel;
    int            *have_it;
    int             pixstart, pixend;
    int             extrapolation;
    int             nref;
    int             nl, nx, ny, pixel;
    int             countLines, usedLines;
    int             uorder;
    int             in, first, last;
    int             width, uradius;
    int             i, j;
    int             null;
    cpl_size        k;


    if (dispersion == 0.0) {
        cpl_msg_error(func, "The expected dispersion (A/pixel) must be given");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (dispersion < 0.0) {
        cpl_msg_error(func, "The expected dispersion must be positive");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    max_disp = dispersion + dispersion * tolerance;
    min_disp = dispersion - dispersion * tolerance;

    if (order < 1) {
        cpl_msg_error(func, "The order of the fitting polynomial "
                      "must be at least 1");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (image == NULL || lines == NULL) {
        cpl_msg_error(func, "Both spectral exposure and reference line "
                      "catalog are required in input");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);
    sdata = cpl_image_get_data_float_const(image);

    nref = cpl_vector_get_size(lines);
    line = cpl_vector_get_data(lines);

    if (*wavestart < 1.0 && *waveend < 1.0) {
        firstLambda = line[0];
        lastLambda = line[nref-1];
        extrapolation = (lastLambda - firstLambda) / 10;
        firstLambda -= extrapolation;
        lastLambda += extrapolation;
        *wavestart = firstLambda;
        *waveend = lastLambda;
    }
    else {
        firstLambda = *wavestart;
        lastLambda = *waveend;
    }

    nl = (lastLambda - firstLambda) / dispersion;
    resampled = cpl_image_new(nl, ny, CPL_TYPE_FLOAT);
    rdata = cpl_image_get_data_float(resampled);

    if (calibration)
        idata = cpl_image_get_data_float(calibration);

    if (residuals)
        ddata = cpl_image_get_data_float(residuals);

    if (idscoeff)
        for (j = 0; j <= order; j++)
            cpl_table_new_column(idscoeff, clab[j], CPL_TYPE_DOUBLE);

    if (restable) {
        cpl_table_set_size(restable, nref);
        cpl_table_new_column(restable, "wavelength", CPL_TYPE_DOUBLE);
        cpl_table_copy_data_double(restable, "wavelength", line);
        for (i = 0; i < ny; i += step) {
             snprintf(name, MAX_COLNAME, "r%d", i);
             cpl_table_new_column(restable, name, CPL_TYPE_DOUBLE);
             snprintf(name, MAX_COLNAME, "d%d", i);
             cpl_table_new_column(restable, name, CPL_TYPE_DOUBLE);
             snprintf(name, MAX_COLNAME, "p%d", i);
             cpl_table_new_column(restable, name, CPL_TYPE_DOUBLE);
        }
    }

    if (detected_lines) {
        cpl_table_set_size(detected_lines, 0);
        cpl_table_new_column(detected_lines, "xpos", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "ypos", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "xpos_iter", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "ypos_iter", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "peak_flux", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "wave_ident", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "wave_ident_iter", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "xpos_fit_rect_wavecal", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "res_xpos", CPL_TYPE_DOUBLE);
    }

    /*
     * Here is the real thing: detecting and identifying peaks,
     * and then fit the transformation from wavelength to pixel
     * and from pixel to wavelength.
     */

    for (i = 0; i < ny; i++) {
        width = mos_lines_width(sdata + i*nx, nx);
        if (sradius > 0) {
            if (width > sradius) {
                uradius = width;
            }
            else {
                uradius = sradius;
            }
        }
        if (width < 5)
            width = 5;
        peaks = mos_peak_candidates(sdata + i*nx, nx, level, width);
        if (peaks) {
            peaks = mos_refine_peaks(sdata + i*nx, nx, peaks, width);
        }
        if (peaks) {
            output = mos_identify_peaks(peaks, lines, 
                                        min_disp, max_disp, r_toleran);
            if (output) {
                cpl_bivector * peaks_ident_used_fit;
                countLines = cpl_bivector_get_size(output);
                if (countLines < 4) {
                    cpl_bivector_delete(output);
                    cpl_vector_delete(peaks);
                    if (nlines)
                        nlines[i] = 0;
                    if (error)
                        error[i] = 0.0;
                    continue;
                }

                /*
                 * Set reference wavelength as zero point
                 */

                wavel = cpl_bivector_get_y(output);
                cpl_vector_subtract_scalar(wavel, refwave);

                uorder = countLines / 2 - 1;
                if (uorder > order)
                    uorder = order;

/* This part is now commented out. In case the first-guess iteration
 * was requested, the first fit was made with a lower polynomial degree:
 * more robust, and still accurate enough to be used as a first-guess.

                if (sradius > 0 && uorder > 2)
                    --uorder;

 * End of commented part */

                ids = mos_poly_wav2pix(output, uorder, reject,
                                       2 * (uorder + 1), &usedLines,
                                       &ids_err, &peaks_ident_used_fit);

                if (ids == NULL) {
                    cpl_bivector_delete(output);
                    cpl_vector_delete(peaks);
                    if (nlines)
                        nlines[i] = 0;
                    if (error)
                        error[i] = 0.0;
                    cpl_error_reset();
                    continue;
                }

                if (idscoeff) {

                    /*
                     * Write it anyway, even in case a first-guess based
                     * solution will be searched afterwards: in case of
                     * failure, the "blind" solution is kept.
                     */

                    for (k = 0; k <= order; k++) {
                        if (k > uorder) {
                            cpl_table_set_double(idscoeff, clab[k], i, 0.0);
                        }
                        else {
                            cpl_table_set_double(idscoeff, clab[k], i,
                                      cpl_polynomial_get_coeff(ids, &k));
                        }
                    }
                }

                if(detected_lines)
                {
                    cpl_size newlines = cpl_vector_get_size(peaks); 
                    cpl_size oldsize = cpl_table_get_nrow(detected_lines); 
                    cpl_table_set_size(detected_lines, oldsize + newlines);
                    for(cpl_size iline = 0; iline < newlines; ++iline)
                    {
                        cpl_table_set_double(detected_lines, "xpos",
                             oldsize + iline, cpl_vector_get(peaks, iline) + 1);
                        cpl_table_set_double(detected_lines, "ypos",
                             oldsize + iline, (double)i + 1);
                        cpl_table_set_double(detected_lines, "peak_flux",
                             oldsize + iline, 
                             sdata[i*nx+(int)(cpl_vector_get(peaks, iline)+0.5)]);
                    }
                }

                //Fill the line identification information in 
                //the detected_lines table
                if(detected_lines)
                {
                    cpl_size nidentlines = cpl_bivector_get_size(output); 
                    cpl_size ndetectlines = cpl_vector_get_size(peaks); 
                    cpl_size totalsize = cpl_table_get_nrow(detected_lines);
                    for(cpl_size idline = 0; idline < nidentlines; ++idline)
                    {
                        for(cpl_size detline = 0; detline < ndetectlines; ++detline)
                        {
                            if(cpl_vector_get(peaks, detline) == 
                               cpl_bivector_get_x_data(output)[idline])
                            {
                                cpl_size table_pos = totalsize - ndetectlines + detline;
                                double wave_ident = cpl_bivector_get_y_data(output)[idline] + refwave;
                                double xpix_fit = cpl_polynomial_eval_1d(ids,
                                        wave_ident - refwave, NULL);
                                double xpos_det = cpl_table_get_double(detected_lines,
                                        "xpos",
                                        table_pos, &null);
                                cpl_table_set_double(detected_lines,
                                                     "wave_ident",
                                                     table_pos,
                                                     wave_ident);
                                cpl_table_set_double(detected_lines,
                                                     "xpos_fit_rect_wavecal",
                                                     table_pos,
                                                     xpix_fit + 1);
                                cpl_table_set_double(detected_lines,
                                                     "res_xpos",
                                                     table_pos,
                                                     xpos_det - xpix_fit - 1);
                               
                            }
                        }
                    }
                }

                if (sradius > 0) {
                    cpl_bivector * peaks_ident_used_fit;

                    /*
                     * Use ids as a first-guess
                     */

                    new_output = mos_find_peaks(sdata + i*nx, nx, lines, 
                                                ids, refwave, uradius);

                    if (new_output) {
                        cpl_bivector_delete(output);
                        output = new_output;
                    }
                    else
                        cpl_error_reset();


                    cpl_polynomial_delete(ids);

                    countLines = cpl_bivector_get_size(output);

                    if (countLines < 4) {
                        cpl_bivector_delete(output);
                        cpl_vector_delete(peaks);

                        /* 
                         * With the following code a decision is taken:
                         * if using the first-guess gives no results,
                         * then also the "blind" solution is rejected.
                         */

                        if (nlines)
                            nlines[i] = 0;
                        if (error)
                            error[i] = 0.0;
                        if (idscoeff)
                            for (k = 0; k <= order; k++)
                                cpl_table_set_invalid(idscoeff, clab[k], i);
                        continue;
                    }

                    wavel = cpl_bivector_get_y(output);
                    cpl_vector_subtract_scalar(wavel, refwave);

                    uorder = countLines / 2 - 1;
                    if (uorder > order)
                        uorder = order;

                    ids = mos_poly_wav2pix(output, uorder, reject,
                                           2 * (uorder + 1), &usedLines,
                                           &ids_err, &peaks_ident_used_fit);

                    if (ids == NULL) {
                        cpl_bivector_delete(output);
                        cpl_vector_delete(peaks);

                        /* 
                         * With the following code a decision is taken:
                         * if using the first-guess gives no results,
                         * then also the "blind" solution is rejected.
                         */

                        if (nlines)
                            nlines[i] = 0;
                        if (error)
                            error[i] = 0.0;
                        if (idscoeff)
                            for (k = 0; k <= order; k++)
                                cpl_table_set_invalid(idscoeff, clab[k], i);
                        cpl_error_reset();
                        continue;
                    }

                    if (idscoeff) {
                        for (k = 0; k <= order; k++) {
                            if (k > uorder) {
                                cpl_table_set_double(idscoeff, clab[k], i, 0.0);
                            }
                            else {
                                cpl_table_set_double(idscoeff, clab[k], i,
                                            cpl_polynomial_get_coeff(ids, &k));
                            }
                        }
                    }

                    
                    if(detected_lines)
                    {
                        cpl_size oldsize = cpl_table_get_nrow(detected_lines); 
                        cpl_size nidentlines = cpl_bivector_get_size(output); 
                        cpl_table_set_size(detected_lines, oldsize + nidentlines);
                        for(cpl_size idline = 0; idline < nidentlines ; ++idline)
                        {
                            double wave_ident = cpl_bivector_get_y_data(output)[idline] + refwave;
                            double xpix_fit = cpl_polynomial_eval_1d(ids,
                                    wave_ident - refwave, NULL);
                            cpl_table_set_double(detected_lines, "xpos_iter",
                                 oldsize + idline, cpl_bivector_get_x_data(output)[idline] + 1);
                            cpl_table_set_double(detected_lines, "ypos_iter",
                                 oldsize + idline, (double)i + 1);
                            cpl_table_set_double(detected_lines, "peak_flux",
                                 oldsize + idline, 
                                 sdata[i*nx+(int)(cpl_bivector_get_x_data(output)[idline]+0.5)]);
                            cpl_table_set_double(detected_lines, "wave_ident_iter",
                                 oldsize + idline, wave_ident);
                            cpl_table_set_double(detected_lines, "xpos_fit_rect_wavecal",
                                 oldsize + idline, xpix_fit + 1);
                        }
                    }

                } /* End of "use ids as a first-guess" */

                if (nlines)
                    nlines[i] = usedLines;
                if (error)
                    error[i] = ids_err / sqrt(usedLines/(uorder + 1));

                pixstart = cpl_polynomial_eval_1d(ids, 
                    cpl_bivector_get_y_data(output)[0], NULL);
                pixend = cpl_polynomial_eval_1d(ids,
                    cpl_bivector_get_y_data(output)[countLines-1], NULL);
                extrapolation = (pixend - pixstart) / 5;
                pixstart -= extrapolation;
                pixend += extrapolation;
                if (pixstart < 0)
                    pixstart = 0;
                if (pixend > nx)
                    pixend = nx;

                /*
                 * Wavelength calibrated image (if requested):
                 */

                if (calibration) {
                    for (j = pixstart; j < pixend; j++) {
                        (idata + i*nx)[j] = mos_eval_dds(ids, firstLambda, 
                                                         lastLambda, refwave, 
                                                         j);
                    }
                }

                /*
                 * Resampled image:
                 */

                for (j = 0; j < nl; j++) {
                    lambda = firstLambda + j * dispersion;
                    fpixel = cpl_polynomial_eval_1d(ids, lambda - refwave, 
                                                    NULL);
                    pixel = fpixel;
                    if (pixel >= 0 && pixel < nx-1) {
                        v1 = (sdata + i*nx)[pixel];
                        v2 = (sdata + i*nx)[pixel+1];
                        vi = v1 + (v2-v1)*(fpixel-pixel);
                        (rdata + i*nl)[j] = vi;
                    }
                }

                /*
                 * Residuals image
                 */

                if (residuals || (restable && !(i%step))) {
                    if (restable && !(i%step)) {
                        lin = cpl_polynomial_new(1);
                        for (k = 0; k < 2; k++)
                            cpl_polynomial_set_coeff(lin, &k, 
                                          cpl_polynomial_get_coeff(ids, &k));
                    }
                    for (j = 0; j < countLines; j++) {
                        pixe = cpl_bivector_get_x_data(output)[j];
                        wave = cpl_bivector_get_y_data(output)[j];
                        value = pixe - cpl_polynomial_eval_1d(ids, wave, NULL);
                        if (residuals) {
                            pixel = pixe + 0.5;
                            (ddata + i*nx)[pixel] = value;
                        }
                        if (restable && !(i%step)) {
                            for (k = 0; k < nref; k++) {
                                if (fabs(line[k] - refwave - wave) < 0.1) {
                                    snprintf(name, MAX_COLNAME, "r%d", i);
                                    cpl_table_set_double(restable, name, 
                                                         k, value);
                                    value = pixe
                                          - cpl_polynomial_eval_1d(lin, wave,
                                                                   NULL);
                                    snprintf(name, MAX_COLNAME, "d%d", i);
                                    cpl_table_set_double(restable, name, 
                                                         k, value);
                                    snprintf(name, MAX_COLNAME, "p%d", i);
                                    cpl_table_set_double(restable, name,
                                                         k, pixe);
                                    break;
                                }
                            }
                        }
                    }
                    if (restable && !(i%step)) {
                        cpl_polynomial_delete(lin);
                    }
                }

                /*
                 * Mask at reference wavelength
                 */

                if (refmask) {
                    mdata = cpl_mask_get_data(refmask);
                    pixel = cpl_polynomial_eval_1d(ids, 0.0, NULL) + 0.5;
                    if (pixel - 1 >= 0 && pixel + 1 < nx) {
                        mdata[pixel-1 + i*nx] = CPL_BINARY_1;
                        mdata[pixel + i*nx] = CPL_BINARY_1;
                        mdata[pixel+1 + i*nx] = CPL_BINARY_1;
                    }
                }

                cpl_polynomial_delete(ids);
                cpl_bivector_delete(output);
            }
            cpl_vector_delete(peaks);
        }
    }

    if (refmask) {
        kernel = cpl_matrix_new(3, 3);
        cpl_matrix_set(kernel, 0, 1, 1.0);
        cpl_matrix_set(kernel, 1, 1, 1.0);
        cpl_matrix_set(kernel, 2, 1, 1.0);

        cpl_mask_dilation(refmask, kernel);
        cpl_mask_erosion(refmask, kernel);
        cpl_mask_erosion(refmask, kernel);
        cpl_mask_dilation(refmask, kernel);

        cpl_matrix_delete(kernel);

        /*
         *  Fill possible gaps
         */

        mdata = cpl_mask_get_data(refmask);
        have_it = cpl_calloc(ny, sizeof(int));

        for (i = 0; i < ny; i++, mdata += nx) {
            for (j = 0; j < nx; j++) {
                if (mdata[j] == CPL_BINARY_1) {
                    have_it[i] = j;
                    break;
                }
            }
        }

        mdata = cpl_mask_get_data(refmask);
        in = 0;
        first = last = 0;

        for (i = 0; i < ny; i++) {
            if (have_it[i]) {
                if (!in) {
                    in = 1;
                    if (first) {
                        last = i;
                        if (abs(have_it[first] - have_it[last]) < 4) {
                            for (j = first; j < last; j++) {
                                mdata[have_it[first] + nx*j + 0] = CPL_BINARY_1;
                                mdata[have_it[first] + nx*j + 1] = CPL_BINARY_1;
                                mdata[have_it[first] + nx*j + 2] = CPL_BINARY_1;
                            }
                        }
                    }
                }
            }
            else {
                if (in) {
                    in = 0;
                    first = i - 1;
                }
            }
        }

        cpl_free(have_it);

    }

/*
    for (i = 0; i < ny; i++) {
        if (nlines[i] == 0) {
            for (k = 0; k <= order; k++) {
                cpl_table_set_invalid(idscoeff, clab[k], i);
            }
        }
    }
*/

    return resampled;
}


/**
 * @brief
 *   Find the location of detected spectra on the CCD
 *
 * @param mask   A reference mask at a given reference wavelength
 *
 * @return Table with characteristics of the detected spectra
 *
 * The input @em mask is the one obtained with the function
 * @c mos_wavelength_calibration_raw().
 * The output table contains the start and end image coordinates
 * of the slit on the input @em mask at reference wavelength.
 * The slits are ordered from top to bottom of the image starting
 * from the first table row.
 *
 * Note that possible gaps within the images of the slits will
 * result in splitting the same slit into two or more sub-slits.
 * This kind of problem is solved within the slit identification
 * task, performed by the function @c mos_identify_slits().
 */

/*
cpl_table *mos_locate_spectra_bis(cpl_mask *mask)
{
    const char *func = "mos_locate_spectra_bis";

    cpl_apertures    *slits;
    cpl_image        *labimage;
    cpl_image        *refimage;
    cpl_binary       *mdata;
    cpl_table        *slitpos;
    cpl_propertylist *sort_col;
    int               nslits;
    int              *have_it;
    int               in, first, last;
    int               i, j;


    if (mask == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nx = cpl_mask_get_size_x(mask);
    ny = cpl_mask_get_size_y(mask);

    mdata = cpl_mask_get_data(refmask);
    have_it = cpl_calloc(ny, sizeof(int));

    for (i = 0; i < ny; i++, mdata += nx) {
        for (j = 0; j < nx; j++) {
            if (mdata[j] == CPL_BINARY_1) {
                have_it[i] = j + 1;
                break;
            }
        }
    }

    mdata = cpl_mask_get_data(refmask);
    in = 0;
    first = last = 0;
    nslits = 0;

    for (i = 0; i < ny; i++) {
        if (have_it[i]) {
            if (in) {
                if (i) {
                    if (abs(have_it[i] - have_it[i-1]) > 3) {
                        nslits++;
                    }
                }
            }
            else {
                in = 1;
                nslits++;
            }
        }
        else {
            if (in) {
                in = 0;
            }
        }
    }
}
*/


/**
 * @brief
 *   Find the location of detected spectra on the CCD
 *
 * @param mask   A reference mask at a given reference wavelength
 *
 * @return Table with characteristics of the detected spectra
 *
 * The input @em mask is the one obtained with the function
 * @c mos_wavelength_calibration_raw().
 * The output table contains the start and end image coordinates 
 * of the slit on the input @em mask at reference wavelength. 
 * The slits are ordered from top to bottom of the image starting 
 * from the first table row.
 *
 * Note that possible gaps within the images of the slits will 
 * result in splitting the same slit into two or more sub-slits.
 * This kind of problem is solved within the slit identification
 * task, performed by the function @c mos_identify_slits().
 */

cpl_table *mos_locate_spectra(cpl_mask *mask)
{
    const char *func = "mos_locate_spectra";

    cpl_apertures    *slits;
    cpl_image        *labimage;
    cpl_image        *refimage;
    cpl_table        *slitpos;
    cpl_propertylist *sort_col;
    cpl_size          nslits;
    int               i;


    if (mask == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    labimage = cpl_image_labelise_mask_create(mask, &nslits);

    if (nslits < 1) {
        cpl_image_delete(labimage);
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    refimage = cpl_image_new_from_mask(mask);

    slits = cpl_apertures_new_from_image(refimage, labimage);

    cpl_image_delete(labimage);
    cpl_image_delete(refimage);

    nslits = cpl_apertures_get_size(slits);  /* Overwriting nslits - safer! */
    if (nslits < 1) {
        cpl_apertures_delete(slits);
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    slitpos = cpl_table_new(nslits);
    cpl_table_new_column(slitpos, "xtop", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slitpos, "ytop", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slitpos, "xbottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slitpos, "ybottom", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(slitpos, "xtop", "pixel");
    cpl_table_set_column_unit(slitpos, "ytop", "pixel");
    cpl_table_set_column_unit(slitpos, "xbottom", "pixel");
    cpl_table_set_column_unit(slitpos, "ybottom", "pixel");

    for (i = 0; i < nslits; i++) {
        cpl_table_set_double(slitpos, "xtop", i, 
                             cpl_apertures_get_top_x(slits, i+1) - 1);
        cpl_table_set_double(slitpos, "ytop", i, 
                             cpl_apertures_get_top(slits, i+1));
        cpl_table_set_double(slitpos, "xbottom", i, 
                             cpl_apertures_get_bottom_x(slits, i+1) - 1);
        cpl_table_set_double(slitpos, "ybottom", i, 
                             cpl_apertures_get_bottom(slits, i+1));
    }

    cpl_apertures_delete(slits);

    sort_col = cpl_propertylist_new();
    cpl_propertylist_append_bool(sort_col, "ytop", 1);
    cpl_table_sort(slitpos, sort_col);
    cpl_propertylist_delete(sort_col);

    return slitpos;

}


/**
 * @brief
 *   Check validity of a slit location table
 *
 * @param slits    Slit location table to validate
 *
 * @return CPL_ERROR_NONE in case of success
 *
 * The input @em slits table is the one obtained with the functions
 * @c mos_locate_spectra() and @c mos_rotate_slits. This table is 
 * expected to contain double precision columns labeled @em xtop, 
 * @em ytop, @em xbottom, and @em ybottom, and this is the check 
 * that is performed here.
 */

cpl_error_code mos_validate_slits(cpl_table *slits) 
{
    const char *func = "mos_validate_slits";


    if (slits == NULL)
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);

    if (1 != cpl_table_has_column(slits, "xtop"))
        return cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);

    if (1 != cpl_table_has_column(slits, "ytop"))
        return cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);

    if (1 != cpl_table_has_column(slits, "xbottom"))
        return cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);

    if (1 != cpl_table_has_column(slits, "ybottom"))
        return cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);

    if (CPL_TYPE_DOUBLE != cpl_table_get_column_type(slits, "xtop"))
        return cpl_error_set(func, CPL_ERROR_INVALID_TYPE);

    if (CPL_TYPE_DOUBLE != cpl_table_get_column_type(slits, "ytop"))
        return cpl_error_set(func, CPL_ERROR_INVALID_TYPE);

    if (CPL_TYPE_DOUBLE != cpl_table_get_column_type(slits, "xbottom"))
        return cpl_error_set(func, CPL_ERROR_INVALID_TYPE);

    if (CPL_TYPE_DOUBLE != cpl_table_get_column_type(slits, "ybottom"))
        return cpl_error_set(func, CPL_ERROR_INVALID_TYPE);

    return CPL_ERROR_NONE;
}


/**
 * @brief
 *   Rotate a slit location table
 *
 * @param slits    Slit location table to rotate
 * @param rotation Rotation angle in multiples of 90 degrees (counterclockwise)
 * @param nx       X size of reference image
 * @param ny       Y size of reference image
 *
 * @return CPL_ERROR_NONE in case of success
 *
 * This function is meant to align the input @em slits table to the mask 
 * coordinates convention used for a given instrument. The input @em slits 
 * table is the one obtained with the function @c mos_locate_spectra(), or 
 * by this self. This table is expected to contain double precision columns 
 * labeled @em xtop, @em ytop, @em xbottom, and @em ybottom, containing 
 * the start and end image coordinates of the slits on the mask obtained 
 * at reference wavelength. The transformation will affect these columns, 
 * other columns are ignored. The reference wavelength image was expected 
 * to be oriented with horizontal dispersion direction and red wavelengths 
 * on the right side, but the @em slits columns are not renamed after a 
 * rotation to reflect the change of convention. If the sizes of the 
 * reference mask, @em nx and @em ny, are greater than 0, the new coordinates
 * will be related to the image rotated in the same way, otherwise a simple 
 * geometrical rotation is applied. The sizes @em nx and @em ny must refer 
 * to the reference image @em after it is rotated in the same way.
 */

cpl_error_code mos_rotate_slits(cpl_table *slits, int rotation, int nx, int ny)
{
    const char *func = "mos_rotate_slits";

    cpl_error_code error;
    char aux_name[] = "_0";
    int i;


    rotation %= 4;
    if (rotation < 0)
        rotation += 4;

    if (rotation == 0)
        return CPL_ERROR_NONE;

    error = mos_validate_slits(slits);
    if (error)
        return cpl_error_set(func, error);

    if (rotation == 1 || rotation == 3) {

        /*
         * Swap x and y column names
         */

        for (i = 0; i < 77; i++)
            if (1 == cpl_table_has_column(slits, aux_name))
                aux_name[1]++;
        if (1 == cpl_table_has_column(slits, aux_name))
            return cpl_error_set(func, CPL_ERROR_CONTINUE);
        cpl_table_name_column(slits, "xtop", aux_name);
        cpl_table_name_column(slits, "ytop", "xtop");
        cpl_table_name_column(slits, aux_name, "ytop");
        cpl_table_name_column(slits, "xbottom", aux_name);
        cpl_table_name_column(slits, "ybottom", "xbottom");
        cpl_table_name_column(slits, aux_name, "ybottom");
    }

    if (rotation == 1 || rotation == 2) {
        cpl_table_multiply_scalar(slits, "xtop", -1.0);
        cpl_table_multiply_scalar(slits, "xbottom", -1.0);
        cpl_table_add_scalar(slits, "xtop", nx);
        cpl_table_add_scalar(slits, "xbottom", nx);
    }

    if (rotation == 3 || rotation == 2) {
        cpl_table_multiply_scalar(slits, "ytop", -1.0);
        cpl_table_multiply_scalar(slits, "ybottom", -1.0);
        cpl_table_add_scalar(slits, "ytop", ny);
        cpl_table_add_scalar(slits, "ybottom", ny);
    }

    return CPL_ERROR_NONE;
}


/**
 * @brief
 *   Identify slits listed in a slit location table
 *
 * @param slits     Slit location on the camera focal plane (CCD pixels).
 * @param maskslits Slit location on the telescope focal plane (mask).
 * @param global    Global distortion table.
 *
 * @return New slit location table 
 *
 * This function is meant to assign to the slits positions listed in 
 * the input @em slits table the slit identifiers contained in the 
 * input @em maskslits table. At least 3 slits should be listed in
 * both tables. The input @em slits table is the one obtained with 
 * the function @c mos_locate_spectra(), with no rotation applied. 
 * This table is expected to contain the double precision columns 
 * labeled @em xtop, @em ytop, @em xbottom, and @em ybottom, 
 * containing the start and end image coordinates of the slits on 
 * the CCD obtained at reference wavelength. The table @em maskslits 
 * is expected to contain the same columns, but with the start 
 * and end coordinates of the slits on the telescope focal plane.
 * The coordinate system should have approximately the same orientation
 * of the input @em slits table, i.e., with horizontal dispersion 
 * direction and red wavelengths dispersed toward the right side. 
 * In addition to the standard columns listed above, the input
 * @em maskslits table should also have a slit identifying integer 
 * column, labeled "slit_id", containing the unique slit identifiers
 * that will be assigned to the identified slits in the input @em slits 
 * table. The construction of this table is instrument dependent, and
 * should be provided at instrument recipe level. 
 * 
 * The output slit location table will contain the same columns as
 * the input tables, with the CCD positions of all the slits listed
 * in the @em maskslits table: such positions are not necessarily
 * all contained in the CCD. The new positions are obtained by mean
 * of a low degree bivariate polynomial model converting from mask
 * positions to CCD positions. This model is derived from a subset
 * of safely identified slits positions. The preliminary identification
 * is performed by matching similar triangles constructed both on the
 * mask and on the CCD of the slits taken three-by-three. Recomputing
 * all positions will remove false detections, join slits containing 
 * gaps, and separate slits that were accidentally joined together.
 *
 * The slit identification may fail in case of masks containing a 
 * regular spacing of slits: such masks would invariably lead to 
 * ambiguous pattern matching, that would not be processed. This
 * would not prevent that data reduction in itself: simply, the
 * reduced spectra would miss their identification. Note that this 
 * is not a real problem, since ambiguous masks are typically masks 
 * used for calibration, and not for scientific observations.
 *
 * In case a @em global distortion table is specified in input, 
 * the coefficients of the bivariate polynomials describing the 
 * relation between mask and CCD coordinates (at reference wavelength)
 * are written to rows 0 and 7.
 */

cpl_table *mos_identify_slits(cpl_table *slits, cpl_table *maskslits,
                              cpl_table *global)
{
    cpl_array        *top_ident = NULL;;
    cpl_array        *bot_ident = NULL;;
    cpl_matrix       *mdata;
    cpl_matrix       *mpattern;
    cpl_matrix       *top_data;
    cpl_matrix       *top_pattern;
    cpl_matrix       *top_mdata;
    cpl_matrix       *top_mpattern;
    cpl_matrix       *bot_data;
    cpl_matrix       *bot_pattern;
    cpl_matrix       *bot_mdata;
    cpl_matrix       *bot_mpattern;
    cpl_propertylist *sort_col;
    double           *xtop;
    double           *ytop;
    double           *xmtop;
    double           *ymtop;
    double           *xbot;
    double           *ybot;
    double           *xmbot;
    double           *ymbot;
    double            top_scale, bot_scale;
    double            angle, top_angle, bot_angle;
    double            xmse, ymse;
    double            xrms, top_xrms, bot_xrms;
    double            yrms, top_yrms, bot_yrms;
    int               nslits, use_data;
    int               nmaskslits, use_pattern;
    int               found_slits, found_slits_top, found_slits_bot;
    int               degree;
    int               i;
    cpl_table        *positions;
    cpl_error_code    error;

    cpl_vector       *point;
    double           *dpoint;
    cpl_vector       *xpos;
    cpl_vector       *ypos;
    cpl_vector       *xmpos;
    cpl_vector       *ympos;
    cpl_bivector     *mpos;
    cpl_polynomial   *xpoly     = NULL;
    cpl_polynomial   *ypoly     = NULL;
    cpl_polynomial   *top_xpoly = NULL;
    cpl_polynomial   *top_ypoly = NULL;
    cpl_polynomial   *bot_xpoly = NULL;
    cpl_polynomial   *bot_ypoly = NULL;

    char *msg_multiplex = " ";



    error = mos_validate_slits(slits);
    if (error) {
        cpl_msg_error(cpl_func, "CCD slits table validation: %s",
                      cpl_error_get_message());
        cpl_error_set(cpl_func, error);
        return NULL;
    }

    error = mos_validate_slits(maskslits);
    if (error) {
        cpl_msg_error(cpl_func, "Mask slits table validation: %s",
                      cpl_error_get_message());
        cpl_error_set(cpl_func, error);
        return NULL;
    }

    if (1 != cpl_table_has_column(maskslits, "slit_id")) {
        cpl_msg_error(cpl_func, "Missing slits identifiers");
        cpl_error_set(cpl_func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    if (CPL_TYPE_INT != cpl_table_get_column_type(maskslits, "slit_id")) {
        cpl_msg_error(cpl_func, "Wrong type used for slits identifiers");
        cpl_error_set(cpl_func, CPL_ERROR_INVALID_TYPE);
        return NULL;
    }

    nslits = cpl_table_get_nrow(slits);
    nmaskslits = cpl_table_get_nrow(maskslits);

    if (nslits == 0 || nmaskslits == 0) {
        cpl_msg_error(cpl_func, "Empty slits table");
        cpl_error_set(cpl_func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (nslits > 200 && mos_multiplex < 0) {
        cpl_msg_info(cpl_func, "Many slits: using 'fast' pattern matching...");
        positions = mos_identify_slits_fast(slits, maskslits, global);
        if (positions == NULL)
            cpl_error_set_where(cpl_func);
        return positions;
    }

    /*
     * Guarantee that both input tables are sorted in the same way
     */

    sort_col = cpl_propertylist_new();

    if (mos_multiplex < 0) {
        cpl_propertylist_append_bool(sort_col, "ytop", 1);
        cpl_table_sort(slits, sort_col);
        cpl_table_sort(maskslits, sort_col);
    }
    else if (nslits > 25) {
        double xmeans = cpl_table_get_column_mean(slits, "xtop");
        double ymeans = cpl_table_get_column_mean(slits, "ytop");
        double xmeanm = cpl_table_get_column_mean(maskslits, "xtop");
        double ymeanm = cpl_table_get_column_mean(maskslits, "ytop");

        cpl_table_duplicate_column(slits, "x", slits, "xtop");
        cpl_table_subtract_scalar(slits, "x", xmeans); // Relative to baricenter
        cpl_table_multiply_columns(slits, "x", "x");   // Squared
    
        cpl_table_duplicate_column(slits, "y", slits, "ytop");
        cpl_table_subtract_scalar(slits, "y", ymeans); // Relative to baricenter
        cpl_table_multiply_columns(slits, "y", "y");   // Squared
    
        cpl_table_add_columns(slits, "x", "y");        // Dist from baricenter

        cpl_table_duplicate_column(maskslits, "x", maskslits, "xtop");
        cpl_table_subtract_scalar(maskslits, "x", xmeanm);
        cpl_table_multiply_columns(maskslits, "x", "x");
    
        cpl_table_duplicate_column(maskslits, "y", maskslits, "ytop");
        cpl_table_subtract_scalar(maskslits, "y", ymeanm);
        cpl_table_multiply_columns(maskslits, "y", "y");
    
        cpl_table_add_columns(maskslits, "x", "y");

        cpl_propertylist_append_bool(sort_col, "x", 0);
        cpl_table_sort(slits, sort_col);
        cpl_table_sort(maskslits, sort_col);
    
        cpl_table_erase_column(slits, "x");
        cpl_table_erase_column(slits, "y");
        cpl_table_erase_column(maskslits, "x");
        cpl_table_erase_column(maskslits, "y");
    }

    cpl_propertylist_delete(sort_col);


    /*
     * First we handle all the special cases (too few slits...)
     */

    if (nslits < 3) {

        /*
         * If there are just 1 or 2 slits on the CCD, and more on the
         * mask, the ambiguity cannot be solved, and an error is returned.
         * This is a case that must be solved with a first-guess relation
         * between mask and CCD.
         */

        if (nslits > 1)
            cpl_msg_warning(cpl_func, "Cannot match the %d found CCD slits "
                            "with the %d mask slits: process will continue "
                            "using the detected CCD slits positions", nslits,
                            nmaskslits);
        else
            cpl_msg_warning(cpl_func, "Cannot match the found CCD slit with "
                            "the %d mask slits: process will continue using "
                            "the detected CCD slit position", nmaskslits);
        return NULL;
    }

    if (nmaskslits < 3) {

        /*
         * If there are less than 3 slits on the mask the ambiguity cannot
         * be solved, and an error is returned. This is a case that must
         * be solved with a first-guess relation between mask and CCD.
         */

        cpl_msg_warning(cpl_func, "Cannot match the %d found CCD slits with "
                        "the %d mask slits: process will continue using "
                        "the detected CCD slits positions", nslits,
                        nmaskslits);
        return NULL;
    }

    /*
     * Pattern matching related operations begin here. Two pattern
     * matching will be run, one based on the "top" and another one
     * based on the "bottom" slit coordinates. The one with the
     * smallest rms in the Y coordinate will be chosen.
     */

    xtop  = cpl_table_get_data_double(slits, "xtop");
    ytop  = cpl_table_get_data_double(slits, "ytop");
    xmtop = cpl_table_get_data_double(maskslits, "xtop");
    ymtop = cpl_table_get_data_double(maskslits, "ytop");

    xbot  = cpl_table_get_data_double(slits, "xbottom");
    ybot  = cpl_table_get_data_double(slits, "ybottom");
    xmbot = cpl_table_get_data_double(maskslits, "xbottom");
    ymbot = cpl_table_get_data_double(maskslits, "ybottom");

    top_data    = cpl_matrix_new(2, nslits);
    top_pattern = cpl_matrix_new(2, nmaskslits);
    bot_data    = cpl_matrix_new(2, nslits);
    bot_pattern = cpl_matrix_new(2, nmaskslits);

    for (i = 0; i < nslits; i++)
        cpl_matrix_set(top_data, 0, i, xtop[i]);

    for (i = 0; i < nslits; i++)
        cpl_matrix_set(top_data, 1, i, ytop[i]);

    for (i = 0; i < nmaskslits; i++)
        cpl_matrix_set(top_pattern, 0, i, xmtop[i]);

    for (i = 0; i < nmaskslits; i++)
        cpl_matrix_set(top_pattern, 1, i, ymtop[i]);

    for (i = 0; i < nslits; i++)
        cpl_matrix_set(bot_data, 0, i, xbot[i]);

    for (i = 0; i < nslits; i++)
        cpl_matrix_set(bot_data, 1, i, ybot[i]);

    for (i = 0; i < nmaskslits; i++)
        cpl_matrix_set(bot_pattern, 0, i, xmbot[i]);

    for (i = 0; i < nmaskslits; i++)
        cpl_matrix_set(bot_pattern, 1, i, ymbot[i]);

    if (mos_multiplex < 0) {
        if (nmaskslits > nslits) { // Safety no longer necessary with CPL 6.+
            use_pattern = nslits;
            use_data = nslits;
        }
        else {
            use_pattern = nmaskslits;
            use_data = nslits;
        }
    }
    else if (nslits > 25) {
        if (nmaskslits > 25) { // Safety no longer necessary with CPL 6.+
            use_pattern = 25;
            if (nslits > 50)
                use_data = 50;
            else
                use_data = nslits;
        }
        else {
            use_pattern = nmaskslits;
            if (nslits > 50)
                use_data = 50;
            else
                use_data = nslits;
        }
    }
    else {
        if (nmaskslits > nslits) { // Safety no longer necessary with CPL 6.+
            use_pattern = nslits;
            use_data = nslits;
        }
        else {
            use_pattern = nmaskslits;
            use_data = nslits;
        }
    }

// printf("use_data = %d/%d, use_pattern = %d/%d\n", use_data, nslits, 
//                                                   use_pattern, nmaskslits);

// cpl_table_save(slits, NULL, NULL, "slits.fits", CPL_IO_DEFAULT);
// cpl_table_save(maskslits, NULL, NULL, "maskslits.fits", CPL_IO_DEFAULT);
// cpl_table_dump(slits, 0, nslits, NULL);
// cpl_table_dump(maskslits, 0, use_pattern, NULL);

    if (mos_multiplex >= 0) {
        sort_col = cpl_propertylist_new();
        cpl_propertylist_append_bool(sort_col, "ytop", 1);
        cpl_table_sort(slits, sort_col);
        cpl_table_sort(maskslits, sort_col);
        cpl_propertylist_delete(sort_col);
    }

    top_ident = cpl_ppm_match_points(top_data, use_data, 3.0, top_pattern,
                                     use_pattern, 0.0, 0.1, 10, &top_mdata,
                                     &top_mpattern, &top_scale, &top_angle);

    bot_ident = cpl_ppm_match_points(bot_data, use_data, 3.0, bot_pattern,
                                     use_pattern, 0.0, 0.1, 10, &bot_mdata,
                                     &bot_mpattern, &bot_scale, &bot_angle);

    if (top_ident == NULL && bot_ident == NULL) {
        cpl_msg_warning(cpl_func, "Pattern matching failure: cannot match "
                        "the %d found CCD slits with the %d mask slits: "
                        "process will continue using the detected CCD "
                        "slits positions", nslits, nmaskslits);
        return NULL;
    }

    found_slits_top = 0;
    found_slits_bot = 0;
    if (top_ident && bot_ident) {
        cpl_msg_info(cpl_func, "Median platescale: %f +/- %f pixel/mm",
                     (top_scale + bot_scale) / 2, fabs(top_scale - bot_scale));
        cpl_msg_info(cpl_func, "Median rotation: %f +/- %f degrees",
                     (top_angle + bot_angle) / 2, fabs(top_angle - bot_angle));
        if (fabs(top_angle) < fabs(bot_angle))
            angle = fabs(top_angle);
        else
            angle = fabs(bot_angle);
        found_slits_top = cpl_matrix_get_ncol(top_mdata);
        found_slits_bot = cpl_matrix_get_ncol(bot_mdata);
    }
    else if (top_ident) {
        cpl_msg_info(cpl_func, "Median platescale: %f pixel/mm", top_scale);
        cpl_msg_info(cpl_func, "Median rotation: %f degrees", top_angle);
        angle = fabs(top_angle);
        found_slits_top = cpl_matrix_get_ncol(top_mdata);
    }
    else {
        cpl_msg_info(cpl_func, "Median platescale: %f pixel/mm", bot_scale);
        cpl_msg_info(cpl_func, "Median rotation: %f degrees", bot_angle);
        angle = fabs(bot_angle);
        found_slits_bot = cpl_matrix_get_ncol(bot_mdata);
    }

    cpl_array_delete(top_ident);
    cpl_array_delete(bot_ident);

    if (angle > 4.0) {
        cpl_msg_warning(cpl_func, "Uncertain pattern matching: the rotation "
                        "angle is expected to be around zero. This match is "
                        "rejected: the process will continue using the %d "
                        "detected CCD slits positions", nslits);
        return NULL;
    }

    found_slits = found_slits_top;
    if (found_slits < found_slits_bot)
        found_slits = found_slits_bot;     /* Max value */

    if (found_slits < 4) {
        cpl_msg_warning(cpl_func,
                        "Too few safely identified slits: %d out of %d "
                        "candidates (%d expected). Process will continue "
                        "using the detected CCD slits positions", found_slits,
                        nslits, nmaskslits);
        return NULL;
    }

    cpl_msg_info(cpl_func, "Preliminary identified slits: %d out of %d "
                 "candidates\n(%d expected)", found_slits, nslits,
                 nmaskslits);

    if (found_slits_top < 4)
        found_slits_top = 0;

    if (found_slits_bot < 4)
        found_slits_bot = 0;

    /*
     * Now for each set select the points of the identified slits, and fit
     * two bivariate polynomials to determine a first approximate relation
     * between positions on the mask and positions on the CCD.
     */

    for (i = 0; i < 2; i++) {
        if (i) {
            found_slits = found_slits_top;
            mdata = top_mdata;
            mpattern = top_mpattern;
        }
        else {
            found_slits = found_slits_bot;
            mdata = bot_mdata;
            mpattern = bot_mpattern;
        }

        if (found_slits == 0)
            continue;
        else if (found_slits < 10)
            degree = 1;
        else
            degree = 2;

        xpos  = cpl_vector_wrap(found_slits,
                                cpl_matrix_get_data(mdata)              );
        ypos  = cpl_vector_wrap(found_slits,
                                cpl_matrix_get_data(mdata) + found_slits);
        xmpos = cpl_vector_wrap(found_slits,
                                cpl_matrix_get_data(mpattern)              );
        ympos = cpl_vector_wrap(found_slits,
                                cpl_matrix_get_data(mpattern) + found_slits);
        mpos  = cpl_bivector_wrap_vectors(xmpos, ympos);
        xpoly = cpl_polynomial_fit_2d_create(mpos, xpos, degree, &xmse);
        ypoly = cpl_polynomial_fit_2d_create(mpos, ypos, degree, &ymse);

        cpl_bivector_unwrap_vectors(mpos);
        cpl_vector_unwrap(xpos);
        cpl_vector_unwrap(ypos);
        cpl_vector_unwrap(xmpos);
        cpl_vector_unwrap(ympos);
        cpl_matrix_delete(mdata);
        cpl_matrix_delete(mpattern);

        if (i) {
            top_xpoly = xpoly;
            top_ypoly = ypoly;
            top_xrms = sqrt(xmse*2*degree/(found_slits - 1));
            top_yrms = sqrt(ymse*2*degree/(found_slits - 1));
        }
        else {
            bot_xpoly = xpoly;
            bot_ypoly = ypoly;
            bot_xrms = sqrt(xmse*2*degree/(found_slits - 1));
            bot_yrms = sqrt(ymse*2*degree/(found_slits - 1));
        }
    }

    if (top_xpoly && bot_xpoly) {
        if (top_xrms < bot_xrms) {  /* top X solution wins... */
            xrms = top_xrms;
            xpoly = top_xpoly;
            cpl_polynomial_delete(bot_xpoly);
        }
        else {                      /* bottom X solution wins... */
            xrms = bot_xrms;
            xpoly = bot_xpoly;
            cpl_polynomial_delete(top_xpoly);
        }
    }
    else if (top_xpoly) {
        xrms = top_xrms;
        xpoly = top_xpoly;
    }
    else {
        xrms = bot_xrms;
        xpoly = bot_xpoly;
    }

    if (top_ypoly && bot_ypoly) {
        if (top_yrms < bot_yrms) {  /* top Y solution wins... */
            yrms = top_yrms;
            ypoly = top_ypoly;
            cpl_polynomial_delete(bot_ypoly);
        }
        else {                      /* bottom Y solution wins... */
            yrms = bot_yrms;
            ypoly = bot_ypoly;
            cpl_polynomial_delete(top_ypoly);
        }
    }
    else if (top_ypoly) {
        yrms = top_yrms;
        ypoly = top_ypoly;
    }
    else {
        yrms = bot_yrms;
        ypoly = bot_ypoly;
    }

    if (xpoly == NULL || ypoly == NULL) {
        cpl_msg_warning(cpl_func, "Fit failure: the accuracy of the "
                        "identified slits positions cannot be improved.");
        cpl_polynomial_delete(xpoly);
        cpl_polynomial_delete(ypoly);
        cpl_error_reset();
        return NULL;
    }

    cpl_msg_info(cpl_func,
                 "Fit successful: X rms = %.3g, Y rms = %.3g (pixel)",
                 xrms, yrms);

    if (global) {
        write_global_distortion(global, 0, xpoly);
        write_global_distortion(global, 7, ypoly);
    }

    /*
     * The fit was successful: use the polynomials to obtain a new
     * position table with the improved positions of the slits
     */

    positions = cpl_table_duplicate(maskslits);
    cpl_table_duplicate_column(positions, "xmtop", positions, "xtop");
    cpl_table_duplicate_column(positions, "ymtop", positions, "ytop");
    cpl_table_duplicate_column(positions, "xmbottom", positions, "xbottom");
    cpl_table_duplicate_column(positions, "ymbottom", positions, "ybottom");

    point = cpl_vector_new(2);
    dpoint = cpl_vector_get_data(point);

    for (i = 0; i < nmaskslits; i++) {
        double position;

        dpoint[0] = cpl_table_get_double(positions, "xmtop", i, NULL);
        dpoint[1] = cpl_table_get_double(positions, "ymtop", i, NULL);
        position  = cpl_polynomial_eval(xpoly, point);
//        if (mos_multiplex >= 0) {
//            if (mos_multiplex != ((int)floor(position)) / mos_region_size) {
//                cpl_table_unselect_row(positions, i);
//                continue;
//            }
//        }
        cpl_table_set_double(positions, "xtop", i, position);
        position  = cpl_polynomial_eval(ypoly, point);
        cpl_table_set_double(positions, "ytop", i, position);
        dpoint[0] = cpl_table_get_double(positions, "xmbottom", i, NULL);
        dpoint[1] = cpl_table_get_double(positions, "ymbottom", i, NULL);
        position  = cpl_polynomial_eval(xpoly, point);
        cpl_table_set_double(positions, "xbottom", i, position);
        position  = cpl_polynomial_eval(ypoly, point);
        cpl_table_set_double(positions, "ybottom", i, position);
    }

//    if (mos_multiplex >= 0) {
//        cpl_table_not_selected(positions);
//        cpl_table_erase_selected(positions);
//        nmaskslits = cpl_table_get_nrow(positions);
//    }

    cpl_vector_delete(point);
    cpl_polynomial_delete(xpoly);
    cpl_polynomial_delete(ypoly);

    cpl_table_erase_column(positions, "xmtop");
    cpl_table_erase_column(positions, "ymtop");
    cpl_table_erase_column(positions, "xmbottom");
    cpl_table_erase_column(positions, "ymbottom");

//    if (mos_multiplex >= 0) {
//        msg_multiplex = 
//        cpl_sprintf("in the CCD section between %d and %d pixel", 
//                    mos_multiplex * mos_region_size, 
//                    (mos_multiplex + 1) * mos_region_size);
//    }

    if (nmaskslits > nslits)
        cpl_msg_info(cpl_func,
                     "Finally identified slits: %d out of %d expected %s\n"
                     "(%d recovered)", nmaskslits, nmaskslits, msg_multiplex,
                     nmaskslits - nslits);
    else if (nmaskslits < nslits)
        cpl_msg_info(cpl_func,
                     "Finally identified slits: %d out of %d expected %s\n"
                     "(%d rejected)", nmaskslits, nmaskslits, msg_multiplex,
                     nslits - nmaskslits);
    else
        cpl_msg_info(cpl_func,
                     "Finally identified slits: %d out of %d expected %s",
                     nmaskslits, nmaskslits, msg_multiplex);

//    if (mos_multiplex >= 0) {
//        cpl_free(msg_multiplex);
//    }

    return positions;

}


cpl_table *mos_identify_slits_fast(cpl_table *slits, cpl_table *maskslits,
                                   cpl_table *global)
{
    const char *func = "mos_identify_slits_fast";

    cpl_propertylist *sort_col;
    cpl_table        *positions;
    cpl_vector       *scales;
    cpl_vector       *angles;
    cpl_vector       *point;
    cpl_vector       *xpos;
    cpl_vector       *ypos;
    cpl_vector       *xmpos;
    cpl_vector       *ympos;
    cpl_bivector     *mpos;
    cpl_polynomial   *xpoly = NULL;
    cpl_polynomial   *ypoly = NULL;
    cpl_error_code    error;
    int nslits;
    int nmaskslits;
    int found_slits;
    int i, j, k;

    double  dist1, dist2, dist3, dist, mindist;
    double  scale, minscale, maxscale;
    double  angle, minangle, maxangle;
    double *dscale;
    double *dangle;
    double *dpoint;
    double *xtop;
    double *ytop;
    double *xbottom;
    double *ybottom;
    double *xcenter;
    double *ycenter;
    double *xpseudo;
    double *ypseudo;
    int    *slit_id;
    double *xmtop;
    double *ymtop;
    double *xmbottom;
    double *ymbottom;
    double *xmcenter;
    double *ymcenter;
    double *xmpseudo;
    double *ympseudo;
    double  xmse, ymse;
    int    *mslit_id;
    int    *good;
    int     minpos;
    int     degree;

    double  sradius = 0.01;   /* Candidate input argument... */
    int     in_sradius;

    double pi = 3.14159265358979323846;


    error = mos_validate_slits(slits);
    if (error) {
        cpl_msg_error(func, "CCD slits table validation: %s", 
                      cpl_error_get_message());
        cpl_error_set(func, error);
        return NULL;
    }

    error = mos_validate_slits(maskslits);
    if (error) {
        cpl_msg_error(func, "Mask slits table validation: %s", 
                      cpl_error_get_message());
        cpl_error_set(func, error);
        return NULL;
    }

    if (1 != cpl_table_has_column(maskslits, "slit_id")) {
        cpl_msg_error(func, "Missing slits identifiers");
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    if (CPL_TYPE_INT != cpl_table_get_column_type(maskslits, "slit_id")) {
        cpl_msg_error(func, "Wrong type used for slits identifiers");
        cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
        return NULL;
    }

    nslits = cpl_table_get_nrow(slits);
    nmaskslits = cpl_table_get_nrow(maskslits);

    if (nslits == 0 || nmaskslits == 0) {
        cpl_msg_error(func, "Empty slits table");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }


    /*
     * Compute middle point coordinates for each slit listed in both
     * input tables.
     */

    if (cpl_table_has_column(slits, "xcenter"))
        cpl_table_erase_column(slits, "xcenter");

    if (cpl_table_has_column(slits, "ycenter"))
        cpl_table_erase_column(slits, "ycenter");

    if (cpl_table_has_column(maskslits, "xcenter"))
        cpl_table_erase_column(maskslits, "xcenter");

    if (cpl_table_has_column(maskslits, "ycenter"))
        cpl_table_erase_column(maskslits, "ycenter");

    cpl_table_duplicate_column(slits, "xcenter", slits, "xtop");
    cpl_table_add_columns(slits, "xcenter", "xbottom");
    cpl_table_divide_scalar(slits, "xcenter", 2.0);
    cpl_table_duplicate_column(slits, "ycenter", slits, "ytop");
    cpl_table_add_columns(slits, "ycenter", "ybottom");
    cpl_table_divide_scalar(slits, "ycenter", 2.0);

    cpl_table_duplicate_column(maskslits, "xcenter", maskslits, "xtop");
    cpl_table_add_columns(maskslits, "xcenter", "xbottom");
    cpl_table_divide_scalar(maskslits, "xcenter", 2.0);
    cpl_table_duplicate_column(maskslits, "ycenter", maskslits, "ytop");
    cpl_table_add_columns(maskslits, "ycenter", "ybottom");
    cpl_table_divide_scalar(maskslits, "ycenter", 2.0);


    /*
     * Guarantee that both input tables are sorted in the same way
     */

    sort_col = cpl_propertylist_new();
    cpl_propertylist_append_bool(sort_col, "ycenter", 1);
    cpl_table_sort(slits, sort_col);
    cpl_table_sort(maskslits, sort_col);
    cpl_propertylist_delete(sort_col);


    /*
     * First we handle all the special cases (too few slits...)
     */

    if (nslits < 3 && nmaskslits > nslits) {

        /*
         * If there are just 1 or 2 slits on the CCD, and more on the
         * mask, the ambiguity cannot be solved, and an error is returned.
         * This is a case that must be solved with a first-guess relation
         * between mask and CCD.
         */

        if (nslits > 1)
            cpl_msg_warning(func, "Cannot match the found CCD slit with the "
                            "%d mask slits: process will continue using the "
                            "detected CCD slit position", nmaskslits);
        else
            cpl_msg_warning(func, "Cannot match the %d found CCD slits with "
                            "the %d mask slits: process will continue using "
                            "the detected CCD slits positions", nslits, 
                            nmaskslits);
        return NULL;
    }

    if (nslits <= 3 && nslits == nmaskslits) {

        cpl_msg_warning(func, "Too few slits (%d) on mask and CCD", nslits);
        cpl_msg_warning(func, "Their detected positions are left unchanged");

        /*
         * If there are just up to 3 slits, both on the mask and on the CCD,
         * we can reasonably hope that those slits were found, and accept 
         * that their positions on the CCD cannot be improved. We prepare
         * therefore an output position table containing the slits with
         * their original positions. We can however give an estimate of
         * the platescale if there is more than one slit.
         */

        positions = cpl_table_duplicate(slits);
        cpl_table_erase_column(slits, "xcenter");
        cpl_table_erase_column(slits, "ycenter");
        cpl_table_duplicate_column(positions, "xmtop", maskslits, "xtop");
        cpl_table_duplicate_column(positions, "ymtop", maskslits, "ytop");
        cpl_table_duplicate_column(positions, "xmbottom", maskslits, "xbottom");
        cpl_table_duplicate_column(positions, "ymbottom", maskslits, "ybottom");
        cpl_table_duplicate_column(positions, "xmcenter", maskslits, "xcenter");
        cpl_table_duplicate_column(positions, "ymcenter", maskslits, "ycenter");
        cpl_table_duplicate_column(positions, "slit_id", maskslits, "slit_id");
        cpl_table_erase_column(maskslits, "xcenter");
        cpl_table_erase_column(maskslits, "ycenter");

        if (nslits > 1) {
            xcenter = cpl_table_get_data_double(positions, "xcenter");
            ycenter = cpl_table_get_data_double(positions, "ycenter");
            xmcenter = cpl_table_get_data_double(positions, "xmcenter");
            ymcenter = cpl_table_get_data_double(positions, "ymcenter");

            dist1 = (xcenter[0] - xcenter[1])*(xcenter[0] - xcenter[1])
                  + (ycenter[0] - ycenter[1])*(ycenter[0] - ycenter[1]);
            dist2 = (xmcenter[0] - xmcenter[1])*(xmcenter[0] - xmcenter[1])
                  + (ymcenter[0] - ymcenter[1])*(ymcenter[0] - ymcenter[1]);
            scale = sqrt(dist1/dist2);

            if (nslits == 3) {
                dist1 = (xcenter[1] - xcenter[2])*(xcenter[1] - xcenter[2])
                      + (ycenter[1] - ycenter[2])*(ycenter[1] - ycenter[2]);
                dist2 = (xmcenter[1] - xmcenter[2])*(xmcenter[1] - xmcenter[2])
                      + (ymcenter[1] - ymcenter[2])*(ymcenter[1] - ymcenter[2]);
                scale += sqrt(dist1/dist2);
                scale /= 2;
            }

            cpl_msg_info(func, "Platescale: %f pixel/mm", scale);
        }

        return positions;
    }

    if (nmaskslits < 3 && nslits > nmaskslits) {

        /*
         * If there are less than 3 slits on the mask the ambiguity cannot 
         * be solved, and an error is returned. This is a case that must 
         * be solved with a first-guess relation between mask and CCD.
         */

        cpl_msg_warning(func, "Cannot match the %d found CCD slits with "
                        "the %d mask slits: process will continue using "
                        "the detected CCD slits positions", nslits, 
                        nmaskslits);
        return NULL;
    }


    /*
     * At this point of the program all the region of the plane
     * (nslits, nmaskslits) where either or both mask and CCD display
     * less than 3 slits are handled in some way. It would be better
     * to add in this place a special handling for identifying slits
     * in case of a very reduced number of slits (say, below 6).
     * It is also clear that if there are many more slits on the
     * mask than on the CCD, or many more on the CCD than on the
     * mask, something went deeply wrong with the preliminary 
     * wavelength calibration. Such cases should be handled with
     * a _complete_ pattern-recognition algorithm based on the
     * construction of all possible triangles. For the moment, 
     * we go directly to the limited pattern-recognition applied
     * below, based on triangles build only for consecutive slits.
     * This is reasonably safe, since the preliminary wavelength
     * calibration performed by mos_identify_peaks() is generally
     * robust.
     */


    /*
     * Compute (X, Y) coordinates on pseudo-plane describing the
     * different position ratios of successive slits, in both
     * input tables.
     */

    if (cpl_table_has_column(slits, "xpseudo"))
        cpl_table_erase_column(slits, "xpseudo");

    if (cpl_table_has_column(slits, "ypseudo"))
        cpl_table_erase_column(slits, "ypseudo");

    if (cpl_table_has_column(maskslits, "xpseudo"))
        cpl_table_erase_column(maskslits, "xpseudo");

    if (cpl_table_has_column(maskslits, "ypseudo"))
        cpl_table_erase_column(maskslits, "ypseudo");

    cpl_table_duplicate_column(slits, "xpseudo", slits, "xcenter");
    cpl_table_duplicate_column(slits, "ypseudo", slits, "ycenter");

    xcenter = cpl_table_get_data_double(slits, "xcenter");
    ycenter = cpl_table_get_data_double(slits, "ycenter");
    xpseudo = cpl_table_get_data_double(slits, "xpseudo");
    ypseudo = cpl_table_get_data_double(slits, "ypseudo");

    for (i = 1; i < nslits - 1; i++) {
        dist1 = (xcenter[i-1] - xcenter[i]) * (xcenter[i-1] - xcenter[i])
              + (ycenter[i-1] - ycenter[i]) * (ycenter[i-1] - ycenter[i]);
        dist2 = (xcenter[i-1] - xcenter[i+1]) * (xcenter[i-1] - xcenter[i+1])
              + (ycenter[i-1] - ycenter[i+1]) * (ycenter[i-1] - ycenter[i+1]);
        dist3 = (xcenter[i] - xcenter[i+1]) * (xcenter[i] - xcenter[i+1])
              + (ycenter[i] - ycenter[i+1]) * (ycenter[i] - ycenter[i+1]);
        xpseudo[i] = sqrt(dist1/dist2);
        ypseudo[i] = sqrt(dist3/dist2);
    }

    cpl_table_set_invalid(slits, "xpseudo", 0);
    cpl_table_set_invalid(slits, "xpseudo", nslits-1);
    cpl_table_set_invalid(slits, "ypseudo", 0);
    cpl_table_set_invalid(slits, "ypseudo", nslits-1);

    cpl_table_duplicate_column(maskslits, "xpseudo", maskslits, "xcenter");
    cpl_table_duplicate_column(maskslits, "ypseudo", maskslits, "ycenter");

    xcenter = cpl_table_get_data_double(maskslits, "xcenter");
    ycenter = cpl_table_get_data_double(maskslits, "ycenter");
    xmpseudo = cpl_table_get_data_double(maskslits, "xpseudo");
    ympseudo = cpl_table_get_data_double(maskslits, "ypseudo");

    for (i = 1; i < nmaskslits - 1; i++) {
        dist1 = (xcenter[i-1] - xcenter[i])*(xcenter[i-1] - xcenter[i])
              + (ycenter[i-1] - ycenter[i])*(ycenter[i-1] - ycenter[i]);
        dist2 = (xcenter[i-1] - xcenter[i+1])*(xcenter[i-1] - xcenter[i+1])
              + (ycenter[i-1] - ycenter[i+1])*(ycenter[i-1] - ycenter[i+1]);
        dist3 = (xcenter[i] - xcenter[i+1])*(xcenter[i] - xcenter[i+1])
              + (ycenter[i] - ycenter[i+1])*(ycenter[i] - ycenter[i+1]);
        xmpseudo[i] = sqrt(dist1/dist2);
        ympseudo[i] = sqrt(dist3/dist2);
    }
    
    cpl_table_set_invalid(maskslits, "xpseudo", 0);
    cpl_table_set_invalid(maskslits, "xpseudo", nmaskslits-1);
    cpl_table_set_invalid(maskslits, "ypseudo", 0);
    cpl_table_set_invalid(maskslits, "ypseudo", nmaskslits-1);


    /*
     * For each (X, Y) on the pseudo-plane related to the mask positions,
     * find the closest (X, Y) on the pseudo-plane related to the CCD
     * positions. If the closest point is closer than a given search
     * radius, a triangle has been found, and 3 slits are identified.
     * However, if more than one point is found within the search
     * radius, we have an ambiguity and this is rejected.
     */

    if (cpl_table_has_column(slits, "slit_id"))
        cpl_table_erase_column(slits, "slit_id");
    cpl_table_new_column(slits, "slit_id", CPL_TYPE_INT);
    slit_id = cpl_table_get_data_int(maskslits, "slit_id");

    for (i = 1; i < nmaskslits - 1; i++) {
        in_sradius = 0;
        mindist = (xmpseudo[i] - xpseudo[1]) * (xmpseudo[i] - xpseudo[1])
                + (ympseudo[i] - ypseudo[1]) * (ympseudo[i] - ypseudo[1]);
        minpos = 1;
        if (mindist < sradius*sradius)
            in_sradius++;
        for (j = 2; j < nslits - 1; j++) {
            dist = (xmpseudo[i] - xpseudo[j]) * (xmpseudo[i] - xpseudo[j])
                 + (ympseudo[i] - ypseudo[j]) * (ympseudo[i] - ypseudo[j]);
            if (dist < sradius*sradius)
                in_sradius++;
            if (in_sradius > 1)    /* More than one triangle within radius */
                break;
            if (mindist > dist) {
                mindist = dist;
                minpos = j;
            }
        }

        mindist = sqrt(mindist);

        if (mindist < sradius && in_sradius == 1) {
            cpl_table_set_int(slits, "slit_id", minpos-1, slit_id[i-1]);
            cpl_table_set_int(slits, "slit_id", minpos, slit_id[i]);
            cpl_table_set_int(slits, "slit_id", minpos+1, slit_id[i+1]);
        }
    }


    /*
     * At this point, the slit_id column contains invalid elements 
     * corresponding to unidentified slits.
     */

    found_slits = nslits - cpl_table_count_invalid(slits, "slit_id");

    if (found_slits < 3) {
        cpl_msg_warning(func, "Too few preliminarily identified slits: "
                        "%d out of %d", found_slits, nslits);
        if (nslits == nmaskslits) {
            cpl_msg_warning(func, "(this is not an error, it could be caused "
                            "by a mask with regularly located slits)");
            cpl_msg_warning(func, "The detected slits positions are left "
                            "unchanged");

            /*
             * If there are less than 3 identified slits, this is probably 
             * a mask with regularly spaced slits (leading to an ambiguous
             * pattern). Only in the case all expected slits appear to have 
             * been found on the CCD we can proceed...
             */

            cpl_table_erase_column(slits, "slit_id");
            cpl_table_erase_column(slits, "xpseudo");
            cpl_table_erase_column(slits, "ypseudo");
            positions = cpl_table_duplicate(slits);
            cpl_table_erase_column(slits, "xcenter");
            cpl_table_erase_column(slits, "ycenter");

            cpl_table_erase_column(maskslits, "xpseudo");
            cpl_table_erase_column(maskslits, "ypseudo");
            cpl_table_duplicate_column(positions, "xmtop", 
                                       maskslits, "xtop");
            cpl_table_duplicate_column(positions, "ymtop", 
                                       maskslits, "ytop");
            cpl_table_duplicate_column(positions, "xmbottom", 
                                       maskslits, "xbottom");
            cpl_table_duplicate_column(positions, "ymbottom", 
                                       maskslits, "ybottom");
            cpl_table_duplicate_column(positions, "xmcenter", 
                                       maskslits, "xcenter");
            cpl_table_duplicate_column(positions, "ymcenter", 
                                       maskslits, "ycenter");
            cpl_table_duplicate_column(positions, "slit_id", 
                                       maskslits, "slit_id");
            cpl_table_erase_column(maskslits, "xcenter");
            cpl_table_erase_column(maskslits, "ycenter");
            return positions;
        }
        else {
            cpl_table_erase_column(slits, "slit_id");
            cpl_table_erase_column(slits, "xpseudo");
            cpl_table_erase_column(slits, "ypseudo");
            positions = cpl_table_duplicate(slits);
            cpl_table_erase_column(slits, "xcenter");
            cpl_table_erase_column(slits, "ycenter");
            cpl_msg_warning(func, "(the failure could be caused "
                            "by a mask with regularly located slits)");
            return NULL;
        }
    }
    else {
        cpl_msg_info(func, "Preliminarily identified slits: %d out of %d "
                     "candidates (%d expected)", found_slits, nslits, 
                     nmaskslits);
    }


    /*
     * Create a table with the coordinates of the preliminarily identified 
     * slits, both on CCD and mask. The original order of the slits positions 
     * is preserved.
     */

    positions = cpl_table_new(found_slits);
    cpl_table_new_column(positions, "slit_id", CPL_TYPE_INT);
    cpl_table_new_column(positions, "xtop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "ytop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "xbottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "ybottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "xcenter", CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "ycenter", CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "xmtop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "ymtop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "xmbottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "ymbottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "xmcenter", CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "ymcenter", CPL_TYPE_DOUBLE);
    cpl_table_new_column(positions, "good", CPL_TYPE_INT);
    cpl_table_fill_column_window_int(positions, "good", 0, found_slits, 0);

    slit_id = cpl_table_get_data_int   (slits, "slit_id");
    xtop    = cpl_table_get_data_double(slits, "xtop");
    ytop    = cpl_table_get_data_double(slits, "ytop");
    xbottom = cpl_table_get_data_double(slits, "xbottom");
    ybottom = cpl_table_get_data_double(slits, "ybottom");
    xcenter = cpl_table_get_data_double(slits, "xcenter");
    ycenter = cpl_table_get_data_double(slits, "ycenter");

    mslit_id = cpl_table_get_data_int   (maskslits, "slit_id");
    xmtop    = cpl_table_get_data_double(maskslits, "xtop");
    ymtop    = cpl_table_get_data_double(maskslits, "ytop");
    xmbottom = cpl_table_get_data_double(maskslits, "xbottom");
    ymbottom = cpl_table_get_data_double(maskslits, "ybottom");
    xmcenter = cpl_table_get_data_double(maskslits, "xcenter");
    ymcenter = cpl_table_get_data_double(maskslits, "ycenter");


    /*
     * Transferring the valid slits information to the new table.
     * Note that invalid elements are coded as 0 in the internal
     * buffer, and this is the way they are recognised and excluded.
     */

    k = 0;
    cpl_table_fill_invalid_int(slits, "slit_id", 0);
    for (i = 0; i < nmaskslits; i++) {
        for (j = 0; j < nslits; j++) {
            if (slit_id[j] == 0)
                continue; /* Skip invalid slit */
            if (mslit_id[i] == slit_id[j]) {
                cpl_table_set_int   (positions, "slit_id",  k, slit_id[j]);

                cpl_table_set_double(positions, "xtop",     k, xtop[j]);
                cpl_table_set_double(positions, "ytop",     k, ytop[j]);
                cpl_table_set_double(positions, "xbottom",  k, xbottom[j]);
                cpl_table_set_double(positions, "ybottom",  k, ybottom[j]);
                cpl_table_set_double(positions, "xcenter",  k, xcenter[j]);
                cpl_table_set_double(positions, "ycenter",  k, ycenter[j]);

                cpl_table_set_double(positions, "xmtop",    k, xmtop[i]);
                cpl_table_set_double(positions, "ymtop",    k, ymtop[i]);
                cpl_table_set_double(positions, "xmbottom", k, xmbottom[i]);
                cpl_table_set_double(positions, "ymbottom", k, ymbottom[i]);
                cpl_table_set_double(positions, "xmcenter", k, xmcenter[i]);
                cpl_table_set_double(positions, "ymcenter", k, ymcenter[i]);

                k++;

                break;
            }
        }
    }

    found_slits = k;

    cpl_table_erase_column(slits, "slit_id");
    cpl_table_erase_column(slits, "xpseudo");
    cpl_table_erase_column(slits, "ypseudo");
    cpl_table_erase_column(slits, "xcenter");
    cpl_table_erase_column(slits, "ycenter");
    cpl_table_erase_column(maskslits, "xpseudo");
    cpl_table_erase_column(maskslits, "ypseudo");
    cpl_table_erase_column(maskslits, "xcenter");
    cpl_table_erase_column(maskslits, "ycenter");


    /*
     * Find the median platescale and rotation angle from the identified 
     * slits, and then exclude slits outlaying more than 10% from the 
     * median platescale, and more than 2 degrees from the median
     * rotation angle.
     */

    ytop    = cpl_table_get_data_double(positions, "ytop");
    ybottom = cpl_table_get_data_double(positions, "ybottom");
    xcenter = cpl_table_get_data_double(positions, "xcenter");
    ycenter = cpl_table_get_data_double(positions, "ycenter");
    xmcenter = cpl_table_get_data_double(positions, "xmcenter");
    ymcenter = cpl_table_get_data_double(positions, "ymcenter");

    scales = cpl_vector_new(found_slits - 1);
    dscale = cpl_vector_get_data(scales);
    angles = cpl_vector_new(found_slits - 1);
    dangle = cpl_vector_get_data(angles);

    for (i = 1; i < found_slits; i++) {
        dist1 = (xcenter[i-1] - xcenter[i]) * (xcenter[i-1] - xcenter[i])
              + (ycenter[i-1] - ycenter[i]) * (ycenter[i-1] - ycenter[i]);
        dist2 = (xmcenter[i-1] - xmcenter[i]) * (xmcenter[i-1] - xmcenter[i])
              + (ymcenter[i-1] - ymcenter[i]) * (ymcenter[i-1] - ymcenter[i]);
        dscale[i-1] = sqrt(dist1/dist2);
        dangle[i-1] = atan2(ycenter[i-1] - ycenter[i], 
                            xcenter[i-1] - xcenter[i])
                    - atan2(ymcenter[i-1] - ymcenter[i], 
                            xmcenter[i-1] - xmcenter[i]);
        dangle[i-1] *= 180;
        dangle[i-1] /= pi;
    }

    minscale = cpl_vector_get_min(scales);
    scale = cpl_vector_get_median_const(scales);
    maxscale = cpl_vector_get_max(scales);

    minangle = cpl_vector_get_min(angles);
    angle = cpl_vector_get_median_const(angles);
    maxangle = cpl_vector_get_max(angles);

    cpl_msg_info(func, "Median platescale: %f pixel/mm", scale);
    cpl_msg_info(func, "Minmax platescale: %f, %f pixel/mm", 
                 minscale, maxscale);

    cpl_msg_info(func, "Median rotation: %f degrees", angle);
    cpl_msg_info(func, "Minmax rotation: %f, %f degrees", 
                 minangle, maxangle);

    good = cpl_table_get_data_int(positions, "good");

    good[0] = good[found_slits - 1] = 1;
    for (i = 1; i < found_slits; i++) {
        if (fabs((dscale[i-1] - scale)/scale) < 0.10
         && fabs(dangle[i-1] - angle) < 2) {
            good[i-1]++;
            good[i]++;
        }
    }

    for (i = 0; i < found_slits; i++) {
        if (good[i] < 2)
            good[i] = 0;
        else
            good[i] = 1;
    }

/*
    for (i = 1; i < found_slits; i++)
        if (fabs((dscale[i-1] - scale)/scale) < 0.10)
            good[i-1] = good[i] = 1;
*/

/* DEBUG ************+
    for (i = 0; i < found_slits; i++) {
        if (good[i]) {
            if (i == found_slits - 1)
                printf("include slit %d, prev = %f, %f\n", 
                       i, dscale[i-1], dangle[i-1]);
            else if (i == 0)
                printf("include slit %d, next %f, %f\n", 
                       i, dscale[i], dangle[i]);
            else
                printf("include slit %d, prev = %f, %f, next %f, %f\n", i, 
                       dscale[i-1], dangle[i-1], dscale[i], dangle[i]);
        }
        else {
            if (i == found_slits - 1)
                printf("EXclude slit %d, prev = %f, %f\n", 
                       i, dscale[i-1], dangle[i-1]);
            else if (i == 0)
                printf("EXclude slit %d, next %f, %f\n", 
                       i, dscale[i], dangle[i]);
            else
                printf("EXclude slit %d, prev = %f, %f, next %f, %f\n", i,    
                       dscale[i-1], dangle[i-1], dscale[i], dangle[i]);
        }
    }
+*********** DEBUG */

    cpl_vector_delete(scales);
    cpl_vector_delete(angles);

    cpl_table_and_selected_int(positions, "good", CPL_EQUAL_TO, 0);
    cpl_table_erase_selected(positions);
    cpl_table_erase_column(positions, "good");
    found_slits = cpl_table_get_nrow(positions);

    if (found_slits < 4) {

        /*
         * If the self-consistency check gives such a poor result,
         * something must have gone really wrong in the preliminary
         * wavelength calibration... Nothing can be done.
         */

        cpl_msg_warning(func, "Too few safely identified slits: %d out of %d "
                        "candidates (%d expected). Process will continue "
                        "using the detected CCD slits positions", found_slits, 
                        nslits, nmaskslits);
        cpl_table_delete(positions);
        return NULL;
    }
    else {
        cpl_msg_info(func, "Safely identified slits: %d out of %d "
                     "candidates\n(%d expected)", found_slits, nslits,
                     nmaskslits);
    }


    /*
     * Now select the central points of the identified slits, and
     * fit two bivariate polynomials to determine a first approximate
     * relation between positions on the mask and positions on the CCD.
     */

    xpos = cpl_vector_wrap(found_slits, 
                           cpl_table_get_data_double(positions, "xcenter"));
    ypos = cpl_vector_wrap(found_slits, 
                           cpl_table_get_data_double(positions, "ycenter"));
    xmpos = cpl_vector_wrap(found_slits, 
                            cpl_table_get_data_double(positions, "xmcenter"));
    ympos = cpl_vector_wrap(found_slits, 
                            cpl_table_get_data_double(positions, "ymcenter"));
    mpos = cpl_bivector_wrap_vectors(xmpos, ympos);

    if (found_slits < 10)
        degree = 1;
    else
        degree = 2;

    xpoly = cpl_polynomial_fit_2d_create(mpos, xpos, degree, &xmse);
    if (xpoly != NULL)
        ypoly = cpl_polynomial_fit_2d_create(mpos, ypos, degree, &ymse);
    cpl_bivector_unwrap_vectors(mpos);
    cpl_vector_unwrap(xpos);
    cpl_vector_unwrap(ypos);
    cpl_vector_unwrap(xmpos);
    cpl_vector_unwrap(ympos);
    if (ypoly == NULL) {
        if (found_slits == nmaskslits) {
            cpl_msg_warning(func, "Fit failure: the accuracy of the "
                            "identified slits positions is not improved.");

            /*
             * The determination of the transformation from mask to CCD
             * failed, but since all slits have been already identified
             * this is not a problem: data can be reduced also without
             * such transformation. Simply the accuracy of the slits 
             * positions cannot be improved.
             */ 

        } else {
            cpl_msg_info(func, "Fit failure: not all slits have been "
                         "identified. Process will continue using "
                         "the detected CCD slits positions");
        }

        cpl_polynomial_delete(xpoly);
        return positions;
    }

    cpl_msg_info(func, "Fit successful: X rms = %.3g, Y rms = %.3g (pixel)", 
                 sqrt(xmse), sqrt(ymse));

    if (global) {
        write_global_distortion(global, 0, xpoly);
        write_global_distortion(global, 7, ypoly);
    }

    /*
     * The fit was successful: use the polynomials to obtain a new 
     * position table with the improved positions of the slits
     */

    cpl_table_delete(positions);

    positions = cpl_table_duplicate(maskslits);
    cpl_table_duplicate_column(positions, "xmtop", positions, "xtop");
    cpl_table_duplicate_column(positions, "ymtop", positions, "ytop");
    cpl_table_duplicate_column(positions, "xmbottom", positions, "xbottom");
    cpl_table_duplicate_column(positions, "ymbottom", positions, "ybottom");

    point = cpl_vector_new(2);
    dpoint = cpl_vector_get_data(point);

    for (i = 0; i < nmaskslits; i++) {
        dpoint[0] = cpl_table_get_double(positions, "xmtop", i, NULL);
        dpoint[1] = cpl_table_get_double(positions, "ymtop", i, NULL);
        cpl_table_set_double(positions, "xtop", i, 
                             cpl_polynomial_eval(xpoly, point));
        cpl_table_set_double(positions, "ytop", i, 
                             cpl_polynomial_eval(ypoly, point));
        dpoint[0] = cpl_table_get_double(positions, "xmbottom", i, NULL);
        dpoint[1] = cpl_table_get_double(positions, "ymbottom", i, NULL);
        cpl_table_set_double(positions, "xbottom", i, 
                             cpl_polynomial_eval(xpoly, point));
        cpl_table_set_double(positions, "ybottom", i, 
                             cpl_polynomial_eval(ypoly, point));
    }

    cpl_vector_delete(point);
    cpl_polynomial_delete(xpoly);
    cpl_polynomial_delete(ypoly);

    cpl_table_erase_column(positions, "xmtop");
    cpl_table_erase_column(positions, "ymtop");
    cpl_table_erase_column(positions, "xmbottom");
    cpl_table_erase_column(positions, "ymbottom");

    if (nmaskslits > nslits)
        cpl_msg_info(func, "Finally identified slits: %d out of %d expected\n"
                 "(%d recovered)", nmaskslits, nmaskslits, nmaskslits - nslits);
    else if (nmaskslits < nslits)
        cpl_msg_info(func, "Finally identified slits: %d out of %d expected\n"
                 "(%d rejected)", nmaskslits, nmaskslits, nslits - nmaskslits);
    else
        cpl_msg_info(func, "Finally identified slits: %d out of %d expected",
                 nmaskslits, nmaskslits);

    return positions;
}

/**
 * @brief
 *   Identify slits listed in a slit location table
 *
 * @param slits     Slit location on the camera focal plane (CCD pixels).
 * @param maskslits Slit location on the telescope focal plane (mask).
 *
 * @return New slit location table 
 *
 * This function is meant to assign to the slits positions listed in 
 * the input @em slits table the slit identifiers contained in the 
 * input @em maskslits table. This is a special case of @c mos_ident_slits
 * for LSS-like masks, i.e., all the slits have the same Y position 
 * (X after rotation).
 * The input @em slits table is the one obtained with 
 * the function @c mos_locate_spectra(), with no rotation applied. 
 * This table is expected to contain the double precision columns 
 * labeled @em xtop, @em ytop, @em xbottom, and @em ybottom, 
 * containing the start and end image coordinates of the slits on 
 * the CCD obtained at reference wavelength. The table @em maskslits 
 * is expected to contain the same columns, but with the start 
 * and end coordinates of the slits on the telescope focal plane.
 * The coordinate system should have approximately the same orientation
 * of the input @em slits table, i.e., with horizontal dispersion 
 * direction and red wavelengths dispersed toward the right side. 
 * In addition to the standard columns listed above, the input
 * @em maskslits table should also have a slit identifying integer 
 * column, labeled "slit_id", containing the unique slit identifiers
 * that will be assigned to the identified slits in the input @em slits 
 * table. The construction of this table is instrument dependent, and
 * should be provided at instrument recipe level. 
 * 
 * The output slit location table will contain the same columns as
 * the input tables, with the CCD positions of all the slits listed
 * in the @em maskslits table: such positions are not necessarily
 * all contained in the CCD. The new positions are obtained by mean
 * of a low degree bivariate polynomial model converting from mask
 * positions to CCD positions. This model is derived from a subset
 * of safely identified slits positions. The preliminary identification
 * is performed by matching similar triangles constructed both on the
 * mask and on the CCD of the slits taken three-by-three. Recomputing
 * all positions will remove false detections, join slits containing 
 * gaps, and separate slits that were accidentally joined together.
 *
 */

cpl_table *mos_identify_slits_linear(cpl_table *slits, cpl_table *maskslits)
{
    cpl_propertylist *sort_col;
    int               nslits;
    int               nmaskslits;
    int               i;
    cpl_table        *positions;
    cpl_error_code    error;

    error = mos_validate_slits(slits);
    if (error) {
        cpl_msg_error(cpl_func, "CCD slits table validation: %s",
                      cpl_error_get_message());
        cpl_error_set(cpl_func, error);
        return NULL;
    }

    error = mos_validate_slits(maskslits);
    if (error) {
        cpl_msg_error(cpl_func, "Mask slits table validation: %s",
                      cpl_error_get_message());
        cpl_error_set(cpl_func, error);
        return NULL;
    }

    if (1 != cpl_table_has_column(maskslits, "slit_id")) {
        cpl_msg_error(cpl_func, "Missing slits identifiers");
        cpl_error_set(cpl_func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    if (CPL_TYPE_INT != cpl_table_get_column_type(maskslits, "slit_id")) {
        cpl_msg_error(cpl_func, "Wrong type used for slits identifiers");
        cpl_error_set(cpl_func, CPL_ERROR_INVALID_TYPE);
        return NULL;
    }

    nslits = cpl_table_get_nrow(slits);
    nmaskslits = cpl_table_get_nrow(maskslits);

    if (nslits == 0 || nmaskslits == 0) {
        cpl_msg_error(cpl_func, "Empty slits table");
        cpl_error_set(cpl_func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (nslits != nmaskslits) {
        cpl_msg_error(cpl_func, "Number of detected and nominal slits do not match. "
                "Cannot identify slits");
        cpl_error_set(cpl_func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    /*
     * Guarantee that both input tables are sorted in the same way
     */

    sort_col = cpl_propertylist_new();

    cpl_propertylist_append_bool(sort_col, "ytop", 1);
    cpl_table_sort(slits, sort_col);
    cpl_table_sort(maskslits, sort_col);
    cpl_propertylist_delete(sort_col);

    /*
     * We assume that we can match the slits based on the order of ytop.
     * TODO: Not that this might not be completely true if there are missing slts
     * detected.
     * TODO: Maybe fit a polynomial to get the refined slit positions, 
     * similar to the pattern matching method.
     */

    positions = cpl_table_duplicate(maskslits);
    cpl_table_duplicate_column(positions, "xmtop", positions, "xtop");
    cpl_table_duplicate_column(positions, "ymtop", positions, "ytop");
    cpl_table_duplicate_column(positions, "xmbottom", positions, "xbottom");
    cpl_table_duplicate_column(positions, "ymbottom", positions, "ybottom");

    for (i = 0; i < nmaskslits; i++) {
        double position;

        cpl_table_set_double(positions, "xtop", i, cpl_table_get_double(slits, "xtop", i, NULL));
        cpl_table_set_double(positions, "ytop", i, cpl_table_get_double(slits, "ytop", i, NULL));
        cpl_table_set_double(positions, "xbottom", i, cpl_table_get_double(slits, "xbottom", i, NULL));
        cpl_table_set_double(positions, "ybottom", i, cpl_table_get_double(slits, "ybottom", i, NULL));
    }

    cpl_table_erase_column(positions, "xmtop");
    cpl_table_erase_column(positions, "ymtop");
    cpl_table_erase_column(positions, "xmbottom");
    cpl_table_erase_column(positions, "ymbottom");

    return positions;

}

/**
 * @brief
 *   Trace flat field spectra
 *
 * @param flat       Flat field image
 * @param slits      Slits positions on the CCD
 * @param reference  Reference wavelength
 * @param blue       Start lambda for tracing
 * @param red        End lambda for tracing
 * @param dispersion Mean spectral dispersion
 *
 * @return A table with the tracings for each spectrum
 *
 * The input @em flat field image is expected to be oriented with 
 * horizontal dispersion direction and red wavelengths on the right 
 * side, and it should be of type @c CPL_TYPE_FLOAT. The @em slits
 * table should be the output of either the function @c mos_identify_slits()
 * (if available) or the function @c mos_locate_spectra().
 * 
 * The flat image is shifted one pixel down and is subtracted from
 * the original image. The result is a vertical gradient map. Next,
 * the negative values are forced positive, to obtain an absolute 
 * gradient map. The map is passed with a horizontal median filter, 
 * and after that the gradient peaks are traced starting from the 
 * slits positions listed in the input @em slits table. The number 
 * of pixels to the left and to the right of the reference pixel is 
 * trivially derived from the specified spectral range @em blue to 
 * @em red and @em dispersion.
 *
 * The output table contains the traced spectral edges positions 
 * in CCD (Y) coordinates for each spectrum. The columns are named 
 * after the "slit_id" listed in the input @em slits table, and are 
 * preceded by a "t" for upper edges, and by a "b" for bottom edges: 
 * for instance, the trace corresponding to the upper edge of a spectrum 
 * with id = 123 will be found in a column named "t123". If the "slit_id" 
 * column is missing in the input table, one will be created with 
 * conventional (unique) numbers. One more column will be created, 
 * named "x", listing the corresponding CCD (X) coordinates of the 
 * trace.
 */

cpl_table *mos_trace_flat(cpl_image *flat, cpl_table *slits, double reference,
                          double blue, double red, double dispersion)
{

    const char  *func = "mos_trace_flat";

    cpl_image   *gradient;
    cpl_image   *sgradient;
    float       *dgradient;
    float        level = 500;   /* It was 250... */
    cpl_vector  *row;
    cpl_vector  *srow;
    cpl_vector **peaks;
    double      *peak;
    int         *slit_id;
    float       *g;
    double      *r;
    double      *xtop;
    double      *ytop;
    double      *xbottom;
    double      *ybottom;
    double       min, dist;
    double       sradius;
    double       tolerance;
    double       start_y, prev_y;
    int          minpos;
    int          nslits;
    int          nrows;
    int          step = 10;
    int          filtbox = 15;
    int          nx, ny, npix;
    int          pos, ypos;
    int          npeaks;
    int          pixel_above, pixel_below;
    int          i, j, k, l;
    char         trace_id[MAX_COLNAME];

    cpl_table   *traces;


    if (flat == NULL || slits == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (red - blue < dispersion) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    /*
     * Create a dummy slit_id column if it is missing in the
     * input slits table
     */

    nslits  = cpl_table_get_nrow(slits);
    if (1 != cpl_table_has_column(slits, "slit_id")) {
        cpl_table_new_column(slits, "slit_id", CPL_TYPE_INT);
        for (i = 0; i < nslits; i++)
            cpl_table_set_int(slits, "slit_id", i, -(i+1));  /* it was (i+1) */
    }

    slit_id = cpl_table_get_data_int(slits, "slit_id");

    nx = cpl_image_get_size_x(flat);
    ny = cpl_image_get_size_y(flat);
    npix = nx * ny;

    gradient = cpl_image_duplicate(flat);
    dgradient = cpl_image_get_data_float(gradient);

    for (i = 0; i < ny - 1; i++) {
        k = i * nx;
        for (j = 0; j < nx; j++) {
            l = k + j;
            dgradient[l] = fabs(dgradient[l] - dgradient[l + nx]);
        }
    }

    npix--;
    for (j = 0; j < nx; j++)
        dgradient[npix - j] = 0.0;

    cpl_image_turn(gradient, -1);
    nx = cpl_image_get_size_x(gradient);
    ny = cpl_image_get_size_y(gradient);
    sgradient = mos_image_vertical_median_filter(gradient, 
                                                 filtbox, 0, ny, 0, step);
    cpl_image_delete(gradient);


    /*
     * Remove background from processed image rows
     */

    dgradient = cpl_image_get_data_float(sgradient);

    for (i = 1; i <= ny; i += step) {
        row = cpl_vector_new_from_image_row(sgradient, i);
        srow = cpl_vector_filter_median_create(row, filtbox);
        cpl_vector_subtract(row, srow);
        cpl_vector_delete(srow);
        g = dgradient + (i-1)*nx;
        r = cpl_vector_get_data(row);
        for (j = 0; j < nx; j++)
            g[j] = r[j];
        cpl_vector_delete(row);
    }


    /*
     * Rotate (temporarily) the input slits table, to get coordinates
     * compatible with the rotated gradient image.
     */

    mos_rotate_slits(slits, 1, nx, ny);
    xtop    = cpl_table_get_data_double(slits, "xtop");
    ytop    = cpl_table_get_data_double(slits, "ytop");
    xbottom = cpl_table_get_data_double(slits, "xbottom");
    ybottom = cpl_table_get_data_double(slits, "ybottom");


    /*
     * Get positions of peaks candidates for each processed gradient
     * image row
     */

    peaks = cpl_calloc(ny, sizeof(cpl_vector *));

    for (i = 0; i < ny; i += step) {
        g = dgradient + i*nx;
        peaks[i] = mos_peak_candidates(g, nx, level, 1.0);

        /* I thought this would be required, but apparently I was wrong.
         * Check twice... */
        if (peaks[i])
            cpl_vector_subtract_scalar(peaks[i], 0.5);
         /**/
    }

    cpl_image_delete(sgradient);


    /*
     * Tracing the flat field spectra edges, starting from the
     * slits positions obtained at reference wavelength. The 
     * gradient maximum closest to each slits ends is found and
     * accepted only within a given search radius:
     */

    sradius = 7.0;  /* Pixel  - it was 5.0 */

    /*
     * The tracing proceeds along the processed gradient image rows,
     * above and below the start position, finding the closest peak
     * to the previous obtained position, accepting the new position
     * only if it is closer than a given tolerance:
     */

    tolerance = 0.9;  /* Pixel */

    /*
     * The trace is attempted for a certain number of pixels above
     * and below the position of the reference wavelength:
     */

    pixel_above = STRETCH_FACTOR * (red - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - blue) / dispersion;


    /*
     * Prepare the structure of the output table:
     */

    nrows = (ny-1)/step + 1;
    traces = cpl_table_new(nrows);
    cpl_table_new_column(traces, "x", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(traces, "x", "pixel");
    for (i = 0, j = 0; i < ny; i += step, j++)
        cpl_table_set(traces, "x", j, i);

    for (i = 0; i < nslits; i++) {

        /*
         * Find the closest processed gradient image row
         */

        ypos = ytop[i];

        if (ypos < 0)
            ypos = 0;
        if (ypos >= ny)
            ypos = ny - 1;

        pos = ypos / step;
        pos *= step;

        /*
         * Find the peak in that row that is closest to xtop[i]
         */

        if (peaks[pos]) {
            peak = cpl_vector_get_data(peaks[pos]);
            npeaks = cpl_vector_get_size(peaks[pos]);

            min = fabs(peak[0] - xtop[i]);
            minpos = 0;
            for (j = 1; j < npeaks; j++) {
                dist = fabs(peak[j] - xtop[i]);
                if (min > dist) {
                    min = dist;
                    minpos = j;
                }
            }
        }
        else {
            npeaks = 0;
        }

        snprintf(trace_id, MAX_COLNAME, "t%d", slit_id[i]);
        if (cpl_table_has_column(traces, trace_id)) {
            cpl_table_delete(traces);
            cpl_msg_error(func, 
               "There are slits with the same ID in header!");
            cpl_msg_error(func, 
               "Something is wrong with the header information: either "
               "both slits with ID=%d exist physically, and in that case "
               "a new ID should be assigned to one of the two slits; or "
               "one of the two slits is actually invalid and should be "
               "removed (but I cannot know which one).\n", slit_id[i]);
            cpl_msg_error(func,
               "This is very likely a problem with VMMPS, which cannot "
               "be fixed here, and as such it should be reported "
               "either to Paranal or to usd-help@eso.org\n");
            return NULL;
        }
        cpl_table_new_column(traces, trace_id, CPL_TYPE_DOUBLE);

        if (min > sradius || npeaks == 0) {
            cpl_msg_warning(func, "Cannot find spectrum edge for "
                            "top (or left) end of slit %d", slit_id[i]);
        }
        else {

            /*
             * Add to output table the start y position. Note that
             * y positions are written in coordinates of the unrotated
             * image, i.e., with horizontal dispersion. Currently nx
             * is the x size of the temporarily rotated image (for
             * faster memory access), but it has the sense of a y size.
             */

            cpl_table_set(traces, trace_id, pos/step, nx - peak[minpos]);
            start_y = peak[minpos];

            /*
             * Perform the tracing of current edge. Above:
             */

            prev_y = start_y;

            for (j = pos + step; j < ny; j += step) {
                if (j - pos > pixel_above)
                    break;
                if (peaks[j]) {
                    peak = cpl_vector_get_data(peaks[j]);
                    npeaks = cpl_vector_get_size(peaks[j]);
                    min = fabs(peak[0] - prev_y);
                    minpos = 0;
                    for (k = 1; k < npeaks; k++) {
                        dist = fabs(peak[k] - prev_y);
                        if (min > dist) {
                            min = dist;
                            minpos = k;
                        }
                    }
                    if (min < tolerance) {
                        cpl_table_set(traces, trace_id, j/step, 
                                      nx - peak[minpos]);
                        prev_y = peak[minpos];
                    }
                }
            }

            /*
             * Perform the tracing of current edge. Below:
             */

            prev_y = start_y;

            for (j = pos - step; j >= 0; j -= step) {
                if (pos - j > pixel_below)
                    break;
                if (peaks[j]) {
                    peak = cpl_vector_get_data(peaks[j]);
                    npeaks = cpl_vector_get_size(peaks[j]);
                    min = fabs(peak[0] - prev_y);
                    minpos = 0;
                    for (k = 1; k < npeaks; k++) {
                        dist = fabs(peak[k] - prev_y);
                        if (min > dist) {
                            min = dist;
                            minpos = k;
                        }
                    }
                    if (min < tolerance) {
                        cpl_table_set(traces, trace_id, j/step, 
                                      nx - peak[minpos]);
                        prev_y = peak[minpos];
                    }
                }
            }
        }


        /*
         * Find the peak in that row that is closest to xbottom[i]
         */

        if (peaks[pos]) {
            peak = cpl_vector_get_data(peaks[pos]);
            npeaks = cpl_vector_get_size(peaks[pos]);
    
            min = fabs(peak[0] - xbottom[i]);
            minpos = 0;
            for (j = 1; j < npeaks; j++) {
                dist = fabs(peak[j] - xbottom[i]);
                if (min > dist) {
                    min = dist;
                    minpos = j;
                }
            }
        }
        else {
            npeaks = 0;
        }

        snprintf(trace_id, MAX_COLNAME, "b%d", slit_id[i]);
        cpl_table_new_column(traces, trace_id, CPL_TYPE_DOUBLE);

        if (min > sradius || npeaks == 0) {
            cpl_msg_warning(func, "Cannot find spectrum edge for "
                            "bottom (or right) end of slit %d", slit_id[i]);
        }
        else {

            cpl_table_set(traces, trace_id, pos/step, nx - peak[minpos]);
            start_y = peak[minpos]; 

            /*
             * Perform the tracing of current edge. Above:
             */

            prev_y = start_y;

            for (j = pos + step; j < ny; j += step) {
                if (j - pos > pixel_above)
                    break;
                if (peaks[j]) {
                    peak = cpl_vector_get_data(peaks[j]);
                    npeaks = cpl_vector_get_size(peaks[j]);
                    min = fabs(peak[0] - prev_y);
                    minpos = 0;
                    for (k = 1; k < npeaks; k++) {
                        dist = fabs(peak[k] - prev_y);
                        if (min > dist) {
                            min = dist;
                            minpos = k;
                        }
                    }
                    if (min < tolerance) {
                        cpl_table_set(traces, trace_id, j/step, 
                                      nx - peak[minpos]);
                        prev_y = peak[minpos];
                    }
                }
            }

            /*
             * Perform the tracing of current edge. Below:
             */

            prev_y = start_y;

            for (j = pos - step; j >= 0; j -= step) {
                if (pos - j > pixel_below)
                    break;
                if (peaks[j]) {
                    peak = cpl_vector_get_data(peaks[j]);
                    npeaks = cpl_vector_get_size(peaks[j]);
                    min = fabs(peak[0] - prev_y);
                    minpos = 0;
                    for (k = 1; k < npeaks; k++) {
                        dist = fabs(peak[k] - prev_y);
                        if (min > dist) {
                            min = dist;
                            minpos = k;
                        }
                    }
                    if (min < tolerance) {
                        cpl_table_set(traces, trace_id, j/step, 
                                      nx - peak[minpos]);
                        prev_y = peak[minpos];
                    }
                }
            }
        }

    }   /* End of loop on slits */

    for (i = 0; i < ny; i += step)
        cpl_vector_delete(peaks[i]);
    cpl_free(peaks);

    /*
     * Restore original orientation of slits positions table
     */

    mos_rotate_slits(slits, -1, ny, nx);

    return traces;

}


/**
 * @brief
 *   Fit spectral traces
 *
 * @param slits      Slits positions on the CCD
 * @param traces     Spectral traces
 * @param order      Order of fitting polynomial
 *
 * @return Table with tracing polynomials coefficients
 *
 * The @em traces table should be the product of the function
 * @c mos_trace_flat(), and the @em slits table should be the
 * same used and processed there. The @em order of the fitting polynomial
 * should not be greater than 5. If the fit of a trace is successful,
 * a column containing the residuals of the fit is added to the
 * @em traces table. If the fitted data are contained in a column
 * named "xyz" (in the convention used by the function @c mos_trace_flat()),
 * the column containing the residuals will be named "xyz_res".
 */

cpl_table *mos_poly_trace(cpl_table *slits, cpl_table *traces, int order)
{
    const char *func = "mos_poly_trace";

    cpl_table      *polytraces;
    cpl_table      *dummy;
    cpl_vector     *x;
    cpl_vector     *trace;
    cpl_polynomial *polytrace;
    char            trace_id[MAX_COLNAME];
    char            trace_res[MAX_COLNAME];
    char            trace_mod[MAX_COLNAME];
    const char     *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    double         *xdata;
    int            *slit_id;
    int             nslits;
    int             nrows;
    int             npoints;
    int             i, j;
    cpl_size        k;


    if (traces == NULL || slits == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (order > 5) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    nrows   = cpl_table_get_nrow(traces);
    xdata   = cpl_table_get_data_double(traces, "x");
    nslits  = cpl_table_get_nrow(slits);
    slit_id = cpl_table_get_data_int(slits, "slit_id");

    polytraces = cpl_table_new(2*nslits);
    cpl_table_new_column(polytraces, "slit_id", CPL_TYPE_INT);
    for (i = 0; i <= order; i++)
        cpl_table_new_column(polytraces, clab[i], CPL_TYPE_DOUBLE);

    for (i = 0; i < nslits; i++) {
        for (j = 0; j < 2; j++) {  /* For top and bottom trace of each slit */

            if (j) {
                snprintf(trace_id, MAX_COLNAME, "b%d", slit_id[i]);
                snprintf(trace_res, MAX_COLNAME, "b%d_res", slit_id[i]);
                snprintf(trace_mod, MAX_COLNAME, "b%d_mod", slit_id[i]);
            }
            else {
                snprintf(trace_id, MAX_COLNAME, "t%d", slit_id[i]);
                snprintf(trace_res, MAX_COLNAME, "t%d_res", slit_id[i]);
                snprintf(trace_mod, MAX_COLNAME, "t%d_mod", slit_id[i]);
            }

            cpl_table_set_int(polytraces, "slit_id", 2*i+j, slit_id[i]);

            /*
             * The "dummy" table is just a tool for eliminating invalid
             * points from the vectors to be fitted.
             */

            dummy = cpl_table_new(nrows);
            cpl_table_duplicate_column(dummy, "x", traces, "x");
            cpl_table_duplicate_column(dummy, trace_id, traces, trace_id);
            npoints = nrows - cpl_table_count_invalid(dummy, trace_id);
            if (npoints < 2 * order) {
                cpl_table_delete(dummy);
                continue;
            }
            cpl_table_erase_invalid(dummy);
            x     = cpl_vector_wrap(npoints, 
                                    cpl_table_get_data_double(dummy, "x"));
            trace = cpl_vector_wrap(npoints, 
                                    cpl_table_get_data_double(dummy, trace_id));
            polytrace = cpl_polynomial_fit_1d_create(x, trace, order, NULL);
            cpl_vector_unwrap(x);
            cpl_vector_unwrap(trace);
            cpl_table_delete(dummy);

            /*
             * Screen bad solutions. At the moment, a primitive screening
             * consists in excluding solutions displaying excessive
             * curvature (larger than 1E-5 / pixel)
             */

            k = 2;
            if (fabs(cpl_polynomial_get_coeff(polytrace, &k)) >  1.E-4) {
                cpl_polynomial_delete(polytrace);
                cpl_table_new_column(traces, trace_mod, CPL_TYPE_DOUBLE);
                cpl_table_duplicate_column(traces, trace_res, traces, 
                                           trace_mod);
                if (j) 
                    cpl_msg_warning(func, "Exclude bad curvature solution "
                           "for bottom (right) edge of slit %d", slit_id[i]);
                else
                    cpl_msg_warning(func, "Exclude bad curvature solution "
                                "for top (left) edge of slit %d", slit_id[i]);
                continue;
            }

            /*
             * Write polynomial coefficients to the output table,
             * tagged with the appropriate slit_id
             */

            for (k = 0; k <= order; k++)
                cpl_table_set_double(polytraces, clab[k], 2*i+j, 
                                     cpl_polynomial_get_coeff(polytrace, &k));

            /*
             * Add column of residuals to input traces table
             */

            cpl_table_new_column(traces, trace_mod, CPL_TYPE_DOUBLE);
            cpl_table_set_column_unit(traces, trace_mod, "pixel");

            for (k = 0; k < nrows; k++) {
                cpl_table_set_double(traces, trace_mod, k,
                          cpl_polynomial_eval_1d(polytrace, xdata[k], NULL));
            }

            cpl_polynomial_delete(polytrace);

            cpl_table_duplicate_column(traces, trace_res, traces, trace_mod);
            cpl_table_subtract_columns(traces, trace_res, trace_id);
            cpl_table_multiply_scalar(traces, trace_res, -1.0);

        }
    }

    return polytraces;

}


/**
 * @brief
 *   Recompute tracing coefficients globally
 *
 * @param slits       Slits positions on the CCD
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param mode        0 = do nothing, 1 = fill gaps, 2 = global model
 *
 * @return @c CPL_ERROR_NONE on success
 *
 * The @em polytraces table should be the product of the function
 * @c mos_trace_flat(), and the @em slits table should be the same
 * used and processed there. The trend of the tracing coefficients 
 * as a function of the first (offset) coefficient c0 is modelled
 * by a linear fit. If @em mode is 2, all the coefficients are 
 * recomputed according to this model. If @em mode is 1, just 
 * missing solutions are found by interpolation. If some tracings 
 * are missing from the @em polytraces table, the value of the 
 * coefficient c0 is drawn from the @em y coordinate of the 
 * corresponding slit edge in the @em slits table, and the rest 
 * of the coefficients are derived from it.
 */

cpl_error_code mos_global_trace(cpl_table *slits, cpl_table *polytraces,
                                int mode)
{
    const char *func = "mos_global_trace";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    cpl_table      *table;
    cpl_vector     *c0;
    cpl_vector     *cn;
    cpl_bivector   *list;
/* alternative (not robust)
    cpl_polynomial *poly;
*/

    double *offset;
    double  rms, q, m;

    int order, nrows, nslits;
    int i, j;


    if (polytraces == NULL) {
        cpl_msg_error(func, "Missing spectral curvature table");
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);
    }

    if (slits == NULL) {
        cpl_msg_error(func, "Missing slits positions table");
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);
    }

    nslits = cpl_table_get_nrow(slits);

    table = cpl_table_duplicate(polytraces);
    cpl_table_erase_invalid(table);

    nrows = cpl_table_get_nrow(table);

    if (nrows < 4) {
        cpl_msg_warning(func, "Too few successful spectral curvature tracings "
                      "(%d): the determination of a global curvature model "
                      "failed", nrows);
        return CPL_ERROR_NONE;
    }
    
    order = cpl_table_get_ncol(polytraces) - 2;

    for (i = 0; i <= order; i++) {
        if (!cpl_table_has_column(table, clab[i])) {
            cpl_msg_error(func, "Wrong spectral curvature table");
            return cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        }
    }


    /*
     * Fill in advance the missing offset terms
     */

    for (i = 0; i < nslits; i++) {
        if (!cpl_table_is_valid(polytraces, clab[0], 2*i)) {
            cpl_table_set_double(polytraces, clab[0], 2*i, 
                                cpl_table_get_double(slits, "ytop", i, NULL));
        }
        if (!cpl_table_is_valid(polytraces, clab[0], 2*i+1)) {
            cpl_table_set_double(polytraces, clab[0], 2*i+1,
                             cpl_table_get_double(slits, "ybottom", i, NULL));
        }
    }

    offset = cpl_table_get_data_double(polytraces, clab[0]);


    /*
     * Fit the global model and modify polytraces table accordingly
     */

    c0 = cpl_vector_wrap(nrows, cpl_table_get_data_double(table, clab[0]));

    for (i = 1; i <= order; i++) {
        cn = cpl_vector_wrap(nrows, cpl_table_get_data_double(table, clab[i]));
        list = cpl_bivector_wrap_vectors(c0, cn);
        robustLinearFit(list, &q, &m, &rms);
/* alternative (not robust)
        poly = cpl_polynomial_fit_1d_create(c0, cn, 1, NULL);
*/
        for (j = 0; j < 2*nslits; j++) {
            if (mode == 1)
                if (cpl_table_is_valid(polytraces, clab[i], j))
                    continue;
            cpl_table_set_double(polytraces, clab[i], j, offset[j]*m + q);
/* alternative (not robust)
            cpl_table_set_double(polytraces, clab[i], j, 
                                 cpl_polynomial_eval_1d(poly, offset[j], NULL));
*/
        }
        cpl_bivector_unwrap_vectors(list);
/* alternative (not robust)
        cpl_polynomial_delete(poly);
*/
        cpl_vector_unwrap(cn);
    }

    cpl_vector_unwrap(c0);
    cpl_table_delete(table);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Spatial remapping of CCD spectra eliminating the spectral curvature
 *
 * @param spectra     CCD image of spectra
 * @param slits       Slits positions on the CCD
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param reference   Reference wavelength
 * @param blue        Start lambda for spatial remapping
 * @param red         End lambda for spatial remapping
 * @param dispersion  Mean spectral dispersion
 * @param flux        If zero, flux conservation is not applied
 * @param calibration Returned spatial calibration image
 *
 * @return Image with the spatially resampled spectra
 * 
 * The input @em spectra image is expected to be oriented with
 * horizontal dispersion direction and red wavelengths on the right
 * side, and it should be of type @c CPL_TYPE_FLOAT. The @em slits
 * table should be the output of either the function @c mos_identify_slits()
 * (if available) or the function @c mos_locate_spectra(). The @em polytraces 
 * table is the output of the function @c mos_poly_trace().
 *
 * The spectra are spatially remapped (with a slight oversampling)
 * starting from the pixel containing the reference wavelength
 * (as reported for each entry of the @em slits table). The number 
 * of pixels to the left and to the right of the reference pixel 
 * is trivially derived from the specified spectral range @em blue 
 * to @em red and @em dispersion.
 *
 * At each @em x pixel position the interval between the top and the
 * bottom edges of each spectrum is linearly remapped into N = ceil(t-b)+1
 * spatial pseudo-pixels, where t and b are its edges positions at the CCD 
 * reference pixel (i.e., if the slit spectrum is widening or narrowing
 * along the wavelength range, it is mapped always to the same number 
 * of pseudo-pixels). The returned image will have the same x size of 
 * the input @em spectra image, and a y size equal to the sum of the
 * N spatial sizes of the resampled spectra.
 *
 * Optionally, an image of the spatially calibrated input exposure,
 * @em calibration, can be returned. This image must be pre-allocated, 
 * and should have the same size of the input @em spectra image. If a
 * @c NULL pointer is passed, it is not computed. The @em calibration 
 * image will consist of pixels having the value of their spatial
 * coordinate along the slit they belong to, or the value zero if 
 * the pixels do not belong to any spectrum.
 *
 * In case @em flux is set to a value different from zero, the remapping
 * is performed applying a correction for flux conservation.
 *
 * Here is how the spatial resampling is applied: the value of a spatial
 * pseudo-pixel p (counted from 0, corresponding to the top spectral
 * trace at a given @em x position on the CCD of a given spectrum),
 * is computed by linear interpolation between the two values of the
 * consecutive CCD pixels closest to the position y = t - p*(t-b)/N,
 * with the same meaning for N, b, and t, described above, and with 
 * 0 <= p <= N.
 *
 * On the @em calibration image the original y pixel positions (at
 * a given x) are assigned the value of their distance from the top 
 * spectral trace, computed as p = N*(t-y)/(t-b), where p is now a
 * float, and y an integer varying from ceil(t) to floor(b).
 *
 * The input @em slits table is added two new columns, labeled "position"
 * and "length", reporting the number of the first (bottom) row belonging 
 * to the corresponding slit spectrum, and its extension in pseudo-pixels
 * on the returned spatially rectified image.
 */

cpl_image *mos_spatial_calibration(cpl_image *spectra, cpl_table *slits,
                                   cpl_table *polytraces, double reference,
                                   double blue, double red, double dispersion,
                                   int flux, cpl_image *calibration)
{
    const char *func = "mos_spatial_calibration";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;
    cpl_image     **exslit;
    cpl_image      *resampled;
    float          *data;
    float          *sdata;
    float          *xdata;
    double          vtop, vbot, value;
    double          top, bot;
    double          coeff;
    double          ytop, ybot;
    double          ypos, yfra;
    double          factor;
    int             yint, ysize, yprev;
    int             nslits;
    int             npseudo;
    int            *slit_id;
    int            *length;
    int             nx, ny;
    int             pixel_above, pixel_below, refpixel, start_pixel, end_pixel;
    int             missing_top, missing_bot;
    int             null;
    int             order;
    int             i, j; 
    cpl_size        k;

    int             create_position = 1;


    if (spectra == NULL || slits == NULL || polytraces == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (red - blue < dispersion) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);
    sdata = cpl_image_get_data(spectra);
    if (calibration)
        data = cpl_image_get_data(calibration);

    if (cpl_table_has_column(slits, "position"))
        create_position = 0;

    if (create_position) {
        cpl_table_new_column(slits, "position", CPL_TYPE_INT);
        cpl_table_new_column(slits, "length", CPL_TYPE_INT);
        cpl_table_set_column_unit(slits, "position", "pixel");
        cpl_table_set_column_unit(slits, "length", "pixel");
    }
    else
        length = cpl_table_get_data_int(slits, "length");

    nslits   = cpl_table_get_nrow(slits);
    slit_id  = cpl_table_get_data_int(slits, "slit_id");
    order    = cpl_table_get_ncol(polytraces) - 2;

    /*
     * The spatial resampling is performed for a certain number of 
     * pixels above and below the position of the reference wavelength:
     */

    pixel_above = STRETCH_FACTOR * (red - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - blue) / dispersion;

    exslit = cpl_calloc(nslits, sizeof(cpl_image *));

    for (i = 0; i < nslits; i++) {
        
        if (create_position == 0)
            if (length[i] == 0)
                continue;

        /*
         * Note that the x coordinate of the reference pixels on the CCD
         * is taken arbitrarily at the top end of each slit. This wouldn't
         * be entirely correct in case of curved slits, or in presence of
         * heavy distortions: in such cases the spatial resampling is
         * really performed across a wide range of wavelengths. But
         * the lag between top and bottom spectral curvature models 
         * would introduce even in such cases negligible effects on
         * the spectral spatial resampling.
         */

        refpixel = cpl_table_get_double(slits, "xtop", i, NULL);

        start_pixel = refpixel - pixel_below;
        if (start_pixel < 0)
            start_pixel = 0;

        end_pixel = refpixel + pixel_above;
        if (end_pixel > nx)
            end_pixel = nx;

        /*
         * Recover from the table of spectral curvature coefficients
         * the curvature polynomials.
         */

        missing_top = 0;
        polytop = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i, &null);
            if (null) {
                cpl_polynomial_delete(polytop);
                missing_top = 1;
                break;
            }
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        missing_bot = 0;
        polybot = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i+1, &null);
            if (null) {
                cpl_polynomial_delete(polybot);
                missing_bot = 1;
                break;
            }
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        if (missing_top && missing_bot) {
            cpl_msg_warning(func, "Spatial calibration, slit %d was not "
                            "traced: no extraction!", 
                            slit_id[i]);
            continue;
        }

        /*
         * In case just one of the two edges was not traced, the other
         * edge curvature model is duplicated and shifted to the other
         * end of the slit: better than nothing!
         */

        if (missing_top) {
            cpl_msg_warning(func, "Upper edge of slit %d was not traced: "
                            "the spectral curvature of the lower edge "
                            "is used instead.", slit_id[i]);
            polytop = cpl_polynomial_duplicate(polybot);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polybot, &k);
            coeff += ytop - ybot;
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        if (missing_bot) {
            cpl_msg_warning(func, "Lower edge of slit %d was not traced: "
                            "the spectral curvature of the upper edge "
                            "is used instead.", slit_id[i]);
            polybot = cpl_polynomial_duplicate(polytop);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polytop, &k);
            coeff -= ytop - ybot;
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        /*
         * Allocate image for current extracted slit
         */

        top = cpl_polynomial_eval_1d(polytop, refpixel, NULL);
        bot = cpl_polynomial_eval_1d(polybot, refpixel, NULL);
        npseudo = ceil(top-bot) + 1;

        if (npseudo < 1) {
            cpl_polynomial_delete(polytop);
            cpl_polynomial_delete(polybot);
            cpl_msg_warning(func, "Slit %d was badly traced: no extraction!",
                            slit_id[i]);
            continue;
        }

        exslit[i] = cpl_image_new(nx, npseudo+1, CPL_TYPE_FLOAT);
        xdata = cpl_image_get_data(exslit[i]);

        /*
         * Write interpolated values to slit image.
         */

        for (j = start_pixel; j < end_pixel; j++) {
            top = cpl_polynomial_eval_1d(polytop, j, NULL);
            bot = cpl_polynomial_eval_1d(polybot, j, NULL);
            factor = (top-bot)/npseudo;
            for (k = 0; k <= npseudo; k++) {
                ypos = top - k*factor;
                yint = floor(ypos);
                yfra = ypos - yint;
                if (yint >= 0 && yint < ny-1) {
                    vtop = sdata[j + nx*yint];
                    vbot = sdata[j + nx*(yint+1)];

                    //This means that the top and bottom traces are crossing,
                    //which is physically impossible, so let's set it to 0. 
                    if(factor <= 0 )  
                        value = 0;
                    else if(vtop == FLT_MAX || vbot == FLT_MAX)
                        value = FLT_MAX;
                    else
                    {
                        value = vtop*(1-yfra) + vbot*yfra;
                        if (flux)
                            value *= factor;
                    }
                    xdata[j + nx*(npseudo-k)] = value;
                    if (calibration) {
                        data[j + nx*yint] = (top-yint)/factor;
                        if (k) {

                            /*
                             * This is added to recover lost pixels on
                             * the CCD image (pixels are lost because
                             * the CCD pixels are less than npseudo+1).
                             */

                            if (yprev - yint > 1) {
                                data[j + nx*(yint+1)] = (top-yint-1)/factor;
                            }
                        }
                    }
                }
                yprev = yint;
            }
        }
        cpl_polynomial_delete(polytop);
        cpl_polynomial_delete(polybot);
    }

    /*
     * Now all the slits images are copied to a single image
     */

    ysize = 0;
    for (i = 0; i < nslits; i++)
        if (exslit[i])
            ysize += cpl_image_get_size_y(exslit[i]);

    if(ysize == 0)
        return NULL;
    
    resampled = cpl_image_new(nx, ysize, CPL_TYPE_FLOAT);

    yint = -1;
    for (i = 0; i < nslits; i++) {
        if (exslit[i]) {
            yint += cpl_image_get_size_y(exslit[i]);
            cpl_image_copy(resampled, exslit[i], 1, ysize - yint);
            if (create_position) {
                cpl_table_set_int(slits, "position", i, ysize - yint - 1);
                cpl_table_set_int(slits, "length", i, 
                                  cpl_image_get_size_y(exslit[i]));
            }
            cpl_image_delete(exslit[i]);
        }
        else if (create_position) {
            cpl_table_set_int(slits, "position", i, -1);
            cpl_table_set_int(slits, "length", i, 0);
        }
    }


    /*
     * Elimination of non-traced slits from slit position table: we cannot do
     * it because we would lose sync with polytraces and other slit-oriented
     * tables. COMMENTED OUT.
     * 

    if (create_position) {

        if (cpl_table_and_selected_int(slits, "position", CPL_EQUAL_TO, -1))
            cpl_table_erase_selected(slits);

    }

    */

    cpl_free(exslit);

    return resampled;

}


/**
 * @brief
 *   Derive wavelength calibration from a rectified arc lamp or sky exposure
 *
 * @param image       A rectified arc lamp or sky exposure
 * @param slits       Slits positions table
 * @param lines       List of reference lines wavelengths
 * @param dispersion  Expected value of the dispersion (wavelength units/pixel)
 * @param level       Threshold for peak detection
 * @param sradius     Search radius for expected peaks (pixels)
 * @param order       Degree of fitting polynomial for wavelength calibration
 * @param reject      Max residual tolerated for line rejection (pixels)
 * @param refwave     Reference wavelength
 * @param wavestart   I/O wavelength of first pixel of resampled image
 * @param waveend     I/O wavelength of last pixel of resampled image
 * @param nlines      Returned array of number of lines used for each fit
 * @param error       Returned array of mean accuracies obtained for each fit
 * @param idscoeff    Returned table with IDS polynomials
 * @param calibration Returned wavelength calibration image
 * @param residuals   Returned residuals image
 * @param restable    Returned residuals table
 * @param detected_lines Returned info on lines detected and used for the fit
 * 
 * @return Input image resampled at constant wavelength step.
 *
 * This function applies the same algorithm for line identification
 * that is applied by the function @c mos_wavelength_calibration_raw().
 * The fundamental difference is that the algorithm is here applied 
 * to an image containing just spectra where the spectral curvature 
 * was eliminated, as produced by the function @c mos_spatial_calibration().
 * The input @em slits table should be the same used and processed by the 
 * function @c mos_spatial_calibration(), containing a column named
 * "position" listing the position of the individual spectra in the
 * rectified image.
 *
 * The rows of the input @em image are independently calibrated one by
 * one. The spectral continuum is assumed to have been already removed
 * (and it should be so, if the function @c mos_wavelength_calibration_raw()
 * was earlier called).
 *
 * Optionally, in case @em sradius is positive, the polynomial solutions 
 * obtained for all rows of each individual slit are averaged into a 
 * single solution, that is used as a first-guess that is passed to the 
 * function @c mos_find_peaks(). This function will search again along 
 * each row belonging to that slit the reference lines candidates around 
 * their expected positions, within the specified search radius; the 
 * polynomial fitting is then repeated using the new found positions.
 * This option can be useful for recovering very faint (i.e., below
 * @em level) reference lines, or reference lines that were lost by
 * @c mos_identify_peaks() because of a partially wrong input @em lines
 * list. The first-guess polynomial solution is derived from the initial 
 * polynomial solutions by determining their median coefficients.
 *
 * An array @em nlines, containing the number of lines used for each
 * fit, and an array @em error, containing the mean error of the
 * polynomial models (in pixels), are returned. A fit failure is
 * indicated with the corresponding element of @em nlines set to
 * zero. The mean error of the polynomial model is evaluated by
 * dividing the RMS of the fit residuals by the square root of the
 * number of fitted points divided the degrees of freedom of the model:
 *
 *               mean error = RMS / sqrt(N / (@em order + 1))
 *
 * The arrays @em nlines and @em error must be pre-allocated, and
 * should all have as many elements as the number of rows in the input
 * @em image. If @c NULL pointers are passed, they are not computed.
 *
 * In the table @em idscoeff will be written the polynomial fits coefficients 
 * obtained for each input @em image row: this table must therefore be 
 * preallocated, before calling this function, with the same number of 
 * rows as the input @em image. No columns should be defined in this
 * table: they will be created automatically by this function, and will
 * be labeled c0, c1, c2, ... up to the specified @em order of the
 * fitting polynomial. 
 *
 * As a by-product of the wavelength calibration, the input @em image
 * is resampled at a constant wavelength step, @em dispersion, and is
 * returned by this function. In case of error a @c NULL pointer is
 * returned. If the input arguments @em wavestart and @em waveend
 * are greater than 1.0, they are taken as the spectral interval 
 * where the spectra are resampled; alternatively, the wavelength 
 * range covered by the resampled image is equal to the wavelength 
 * range of the input reference @em lines catalog, extended by 10 
 * percent on the blue and the red sides: such used spectral interval 
 * is then returned via the same variables, @em wavestart and @em waveend.
 * Note that limiting the spectral range doesn't prevents this function
 * to perform a wavelength calibration based on all the reference 
 * wavelengths listed in the input line catalog, including those
 * outside the specified range (if they are found somewhere on the
 * detector).
 *
 * Optionally, an image of the wavelength calibrated input exposure,
 * @em calibration, might be filled with the value of the wavelength 
 * for each pixel. This image may be significantly more accurate than 
 * the @em calibration image obtained with the preliminary wavelength 
 * calibration, returned by the function @c mos_wavelength_calibration_raw(). 
 * For this reason it would also be appropriate to map it back to the
 * original CCD coordinate system, using the function @c mos_map_wavelengths().
 * The image @em calibration should be pre-allocated in order to have 
 * it computed, and it should have the same size of the input @em image 
 * that was used for the calibration. Furthermore, an image of the fit 
 * residuals, @em residuals, can be returned. This image must also be
 * be pre-allocated for having it computed, and should as well have the 
 * same size of the input rectified @em image. If @c NULL pointers are 
 * passed, the @em calibration and @em residuals images are not computed.
 *
 * If detected_lines table is an input, then a table with the position of
 * the detected lines, its peak flux, wavelength identified and final xpos
 * position in the wavelength calibrated image using the wavelength solution
 *  is returned. The table has to be allocated but empty.
 */

cpl_image *mos_wavelength_calibration_final(cpl_image *image, cpl_table *slits,
                                            cpl_vector *lines,
                                            double dispersion, float level,
                                            int sradius, int order,
                                            double reject, double refwave,
                                            double *wavestart, double *waveend,
                                            int *nlines, double *error, 
                                            cpl_table *idscoeff,
                                            cpl_image *calibration,
                                            cpl_image *residuals,
                                            cpl_table *restable,
                                            cpl_table *detected_lines, 
                                            double disp_tolerance,
                                            double ratio_tolerance)
{

    const char *func = "mos_wavelength_calibration_final";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    double  tolerance = disp_tolerance;    
    double  r_toleran = ratio_tolerance;    
    int     step = 10;           /* Compute restable every "step" rows */

    char            name[MAX_COLNAME];

    cpl_image      *resampled;
    cpl_bivector   *peaks_ident;
    cpl_vector     *wavel;
    cpl_vector     *peaks;
    cpl_polynomial *ids;
    cpl_polynomial *lin;
    cpl_polynomial *fguess;
    cpl_table      *coeff;
    double          ids_err;
    double          max_disp, min_disp;
    double         *line;
    double          firstLambda, lastLambda, lambda;
    double          wave, pixe, value;
    double          c;
    float          *sdata;
    float          *rdata;
    float          *idata;
    float          *ddata;
    float           v1, v2, vi;
    float           fpixel;
    int            *length;
    int             pixstart, pixend;
    int             row_top, row_bot;
    int             extrapolation;
    int             nref;
    int             nslits;
    int             nfits;
    int             nl, nx, ny, pixel;
    int             countLines, usedLines;
    int             uorder;
    int             missing;
    int             null;
    int             width, uradius;
    int             i, j, s;
    cpl_size        k;


    if (dispersion == 0.0) {
        cpl_msg_error(func, "The expected dispersion (A/pixel) must be given");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (dispersion < 0.0) {
        cpl_msg_error(func, "The expected dispersion must be positive");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (idscoeff == NULL) {
        cpl_msg_error(func, "A preallocated IDS coeff table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    max_disp = dispersion + dispersion * tolerance;
    min_disp = dispersion - dispersion * tolerance;

    if (order < 1) {
        cpl_msg_error(func, "The order of the fitting polynomial "
                      "must be at least 1");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (image == NULL || lines == NULL) {
        cpl_msg_error(func, "Both spectral exposure and reference line "
                      "catalog are required in input");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);
    sdata = cpl_image_get_data_float(image);

    nref = cpl_vector_get_size(lines);
    line = cpl_vector_get_data(lines);

    if (*wavestart < 1.0 && *waveend < 1.0) {
        firstLambda = line[0];
        lastLambda = line[nref-1];
        extrapolation = (lastLambda - firstLambda) / 10;
        firstLambda -= extrapolation;
        lastLambda += extrapolation;
        *wavestart = firstLambda;
        *waveend = lastLambda;
    }
    else {
        firstLambda = *wavestart;
        lastLambda = *waveend;
    }

    nl = (lastLambda - firstLambda) / dispersion;
    resampled = cpl_image_new(nl, ny, CPL_TYPE_FLOAT);
    rdata = cpl_image_get_data_float(resampled);

    /*
     * Allocate total output table of IDS coefficients
     */

    for (j = 0; j <= order; j++)
        cpl_table_new_column(idscoeff, clab[j], CPL_TYPE_DOUBLE);

    if (calibration)
        idata = cpl_image_get_data_float(calibration);

    if (residuals)
        ddata = cpl_image_get_data_float(residuals);

    if (restable) {
        cpl_table_set_size(restable, nref);
        cpl_table_new_column(restable, "wavelength", CPL_TYPE_DOUBLE);
        cpl_table_copy_data_double(restable, "wavelength", line);
        for (i = 0; i < ny; i += step) {
             snprintf(name, MAX_COLNAME, "r%d", i);
             cpl_table_new_column(restable, name, CPL_TYPE_DOUBLE);
             snprintf(name, MAX_COLNAME, "d%d", i);
             cpl_table_new_column(restable, name, CPL_TYPE_DOUBLE);
             snprintf(name, MAX_COLNAME, "p%d", i);
             cpl_table_new_column(restable, name, CPL_TYPE_DOUBLE);
        }
    }

    if (detected_lines) {
        cpl_table_set_size(detected_lines, 0);
        cpl_table_new_column(detected_lines, "slit_id", CPL_TYPE_INT);
        cpl_table_new_column(detected_lines, "xpos_rectified", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "ypos_rectified", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "xpos_rectified_iter", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "ypos_rectified_iter", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "peak_flux", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "wave_ident", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "wave_ident_iter", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "xpos_fit_rect_wavecal", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "res_xpos", CPL_TYPE_DOUBLE);
        cpl_table_new_column(detected_lines, "fit_used", CPL_TYPE_INT);
    }


    /*
     * Process all slits separately.
     */

    nslits   = cpl_table_get_nrow(slits);
    length   = cpl_table_get_data_int(slits, "length");

    row_top = ny;
    for (s = 0; s < nslits; s++) {

        int slit_id;
        slit_id = cpl_table_get_int(slits, "slit_id", s, NULL);

        if (length[s] == 0)
            continue;

        /*
         * row_top and row_bot define the boundaries of the current slit.
         * Here we begin (arbitrarily...) from the top slit.
         */

        row_bot = cpl_table_get_int(slits, "position", s, NULL);

        if (sradius > 0) {

            /*
             * If a search radius was defined, allocate the table of
             * the fitting polynomials coefficients. This table is
             * just used to generate the first-guess polynomial made
             * of the median coefficients of all polynomials found
             * for this slit.
             */

            coeff = cpl_table_new(row_top - row_bot);
            for (j = 0; j <= order; j++)
                cpl_table_new_column(coeff, clab[j], CPL_TYPE_DOUBLE);
        }

        /*
         * Here is the loop on all rows of the current slit. They are
         * wavelength calibrated one by one.
         */

        for (i = row_bot; i < row_top; i++) {
            width = mos_lines_width(sdata + i*nx, nx);
            if (width < 5)
                width = 5;
            peaks = mos_peak_candidates(sdata + i*nx, nx, level, width);
            if (peaks) {
                peaks = mos_refine_peaks(sdata + i*nx, nx, peaks, width);
            }
            if (peaks) {
                int keep_multiplex = mos_multiplex;
                mos_multiplex = -1;
                if(detected_lines)
                {
                    cpl_size newlines = cpl_vector_get_size(peaks); 
                    cpl_size oldsize = cpl_table_get_nrow(detected_lines); 
                    cpl_table_set_size(detected_lines, oldsize + newlines);
                    for(cpl_size iline = 0; iline < newlines; ++iline)
                    {
                        cpl_table_set_int(detected_lines, "slit_id",
                             oldsize + iline, slit_id);
                        cpl_table_set_double(detected_lines, "xpos_rectified",
                             oldsize + iline, cpl_vector_get(peaks, iline) + 1);
                        cpl_table_set_double(detected_lines, "ypos_rectified",
                             oldsize + iline, (double)i + 1);
                        cpl_table_set_double(detected_lines, "peak_flux",
                             oldsize + iline, 
                             sdata[i*nx+(int)(cpl_vector_get(peaks, iline)+0.5)]);
                        cpl_table_set_int(detected_lines,
                                          "fit_used",
                                          oldsize + iline, 0);
                    }
                }
                peaks_ident = mos_identify_peaks(peaks, lines, 
                                            min_disp, max_disp, r_toleran);
                mos_multiplex = keep_multiplex;
                if (peaks_ident) {
                    cpl_bivector * peaks_ident_used_fit;
                    countLines = cpl_bivector_get_size(peaks_ident);
                    if (countLines < 4) {
                        cpl_bivector_delete(peaks_ident);
                        cpl_vector_delete(peaks);
                        if (nlines)
                            nlines[i] = 0;
                        if (error)
                            error[i] = 0.0;
                        continue;
                    }

                    /*
                     * Set reference wavelength as zero point
                     */

                    wavel = cpl_bivector_get_y(peaks_ident);
                    cpl_vector_subtract_scalar(wavel, refwave);

                    uorder = countLines / 2 - 1;
                    if (uorder > order)
                        uorder = order;

                    ids = mos_poly_wav2pix(peaks_ident, uorder, reject,
                                           2 * (uorder + 1), &usedLines,
                                           &ids_err, &peaks_ident_used_fit);

                    if (ids == NULL) {
                        cpl_bivector_delete(peaks_ident);
                        cpl_vector_delete(peaks);
                        if (nlines)
                            nlines[i] = 0;
                        if (error)
                            error[i] = 0.0;
                        cpl_error_reset();
                        continue;
                    }

                    if (sradius > 0) {
                        for (k = 0; k <= order; k++) {
                            if (k > uorder) {
                                cpl_table_set_double(coeff, clab[k], 
                                i - row_bot, 0.0);
                            }
                            else {
                                cpl_table_set_double(coeff, clab[k], 
                                i - row_bot, cpl_polynomial_get_coeff(ids, &k));
                            }
                        }
                    }
               /*   else {   */
                        if (calibration) {
                            pixstart = cpl_polynomial_eval_1d(ids,
                              cpl_bivector_get_y_data(peaks_ident)[0], 
                              NULL);
                            pixend = cpl_polynomial_eval_1d(ids,
                              cpl_bivector_get_y_data(peaks_ident)[countLines-1],
                              NULL);
                            extrapolation = (pixend - pixstart) / 5;
                            pixstart -= extrapolation;
                            pixend += extrapolation;
                            if (pixstart < 0)
                                pixstart = 0;
                            if (pixend > nx)
                                pixend = nx;
   
                            for (j = pixstart; j < pixend; j++) {
                                (idata + i*nx)[j] = mos_eval_dds(ids, 
                                     firstLambda, lastLambda, refwave, j);
                            }
                        }

                        /*
                         * Residuals image
                         */
        
                        if (residuals || (restable && !(i%step))) {
                            if (restable && !(i%step)) {
                                lin = cpl_polynomial_new(1);
                                for (k = 0; k < 2; k++)
                                    cpl_polynomial_set_coeff(lin, &k,
                                          cpl_polynomial_get_coeff(ids, &k));
                            }
                            for (j = 0; j < countLines; j++) {
                                pixe = cpl_bivector_get_x_data(peaks_ident)[j];
                                wave = cpl_bivector_get_y_data(peaks_ident)[j];
                                value = pixe 
                                      - cpl_polynomial_eval_1d(ids, wave, NULL);
                                if (residuals) {
                                    pixel = pixe + 0.5;
                                    (ddata + i*nx)[pixel] = value;
                                }
                                if (restable && !(i%step)) {
                                    for (k = 0; k < nref; k++) {
                                        if (fabs(line[k]-refwave-wave) < 0.1) {
                                            snprintf(name, MAX_COLNAME, 
                                                     "r%d", i);
                                            cpl_table_set_double(restable, name,
                                                                 k, value);
                                            value = pixe
                                                  - cpl_polynomial_eval_1d(lin,
                                                              wave, NULL);
                                            snprintf(name, MAX_COLNAME, 
                                                     "d%d", i);
                                            cpl_table_set_double(restable, name,
                                                                 k, value);
                                            snprintf(name, MAX_COLNAME,
                                                     "p%d", i);
                                            cpl_table_set_double(restable, name,
                                                                 k, pixe);
                                            break;
                                        }
                                    }
                                }
                            }
                            if (restable && !(i%step)) {
                                cpl_polynomial_delete(lin);
                            }
/***
                            for (j = 0; j < countLines; j++) {
                                pixel = cpl_bivector_get_x_data(output)[j] 
                                      + 0.5;
                                (ddata + i*nx)[pixel] =
                                cpl_bivector_get_x_data(output)[j]
                              - cpl_polynomial_eval_1d(ids,
                                cpl_bivector_get_y_data(output)[j], 
                                NULL);
                            }
***/
                            //Fill the line identification information in 
                            //the detected_lines table
                            if(detected_lines)
                            {
                                cpl_size nidentlines = cpl_bivector_get_size(peaks_ident); 
                                cpl_size ndetectlines = cpl_vector_get_size(peaks); 
                                cpl_size totalsize = cpl_table_get_nrow(detected_lines);
                                for(cpl_size idline = 0; idline < nidentlines; ++idline)
                                {
                                    for(cpl_size detline = 0; detline < ndetectlines; ++detline)
                                    {
                                        if(cpl_vector_get(peaks, detline) == 
                                           cpl_bivector_get_x_data(peaks_ident)[idline])
                                        {
                                            cpl_size table_pos = totalsize - ndetectlines + detline;
                                            double wave_ident = cpl_bivector_get_y_data(peaks_ident)[idline] + refwave;
                                            double xpix_fit = cpl_polynomial_eval_1d(ids,
                                                    wave_ident - refwave, NULL);
                                            double xpos_det = cpl_table_get_double(detected_lines,
                                                    "xpos_rectified",
                                                    table_pos, &null);
                                            cpl_table_set_double(detected_lines,
                                                                 "wave_ident",
                                                                 table_pos,
                                                                 wave_ident);
                                            cpl_table_set_double(detected_lines,
                                                                 "xpos_fit_rect_wavecal",
                                                                 table_pos,
                                                                 xpix_fit + 1);
                                            cpl_table_set_double(detected_lines,
                                                                 "res_xpos",
                                                                 table_pos,
                                                                  xpos_det - xpix_fit - 1);
                                            cpl_table_set_int(detected_lines,
                                                              "fit_used",
                                                              table_pos, 0);
                                            for(cpl_size i_used = 0; i_used < cpl_bivector_get_size(peaks_ident_used_fit); ++i_used)
                                            {
                                                if(cpl_bivector_get_x_data(peaks_ident)[idline] == cpl_bivector_get_x_data(peaks_ident_used_fit)[i_used])
                                                    cpl_table_set_int(detected_lines,
                                                                     "fit_used",
                                                                     table_pos, 1);
                                            }
                                        }                                            
                                    }
                                }
                            }
                        }
                /*  }   */

                    /*
                     * Write it anyway, even in case a first-guess based
                     * solution will be searched afterwards: in case of
                     * failure, the "blind" solution is kept.
                     */

                    if (nlines)
                        nlines[i] = usedLines;
                    if (error)
                        error[i] = ids_err / sqrt(usedLines/(uorder + 1));

                    for (k = 0; k <= order; k++) {
                        if (k > uorder) {
                            cpl_table_set_double(idscoeff, clab[k], i, 0.0);
                        }
                        else {
                            cpl_table_set_double(idscoeff, clab[k], i,
                                      cpl_polynomial_get_coeff(ids, &k));
                        }
                    }

                    cpl_polynomial_delete(ids);
                    cpl_bivector_delete(peaks_ident);
                }
                cpl_vector_delete(peaks);
            }
        }       /* End of loop on current slit's rows */


        if (sradius > 0) {

            /*
             * See whether there are valid fits at all...
             */
    
            nfits = row_top - row_bot - cpl_table_count_invalid(coeff, clab[0]);
    
            if (nfits) {
                int slope = 0;

                fguess = cpl_polynomial_new(1);

                if (mos_interpolate_wavecalib_mos(coeff, 2, 1)) {

                    slope = 0;

                    /*
                     * Compute a median IDS polynomial for the current slit
                     */

                    for (k = 0; k <= order; k++) {
                        c = cpl_table_get_column_median(coeff, clab[k]);
                        cpl_polynomial_set_coeff(fguess, &k, c);
                    }
                }
                else {
                    slope = 1;
                }

                for (i = row_bot; i < row_top; i++) {
                    cpl_bivector * peaks_ident_used_fit;

                    /*
                     * Use first-guess to find the reference lines again
                     */

                    width = mos_lines_width(sdata + i*nx, nx);
                    if (width > sradius) {
                        uradius = width; 
                    }
                    else {
                        uradius = sradius;
                    }

                    if (slope) {
                        for (k = 0; k <= order; k++) {
                            c = cpl_table_get_double(coeff, clab[k], 
                                                     i - row_bot, NULL);
                            cpl_polynomial_set_coeff(fguess, &k, c);
                        }
                    }

                    peaks_ident = mos_find_peaks(sdata + i*nx, nx, lines, 
                                            fguess, refwave, uradius);

                    if (peaks_ident == NULL) {
                        cpl_error_reset();
                        continue;
                    }

                    countLines = cpl_bivector_get_size(peaks_ident);

                    if (countLines < 4) {
                        cpl_bivector_delete(peaks_ident);
                        continue;
                    }

                    /*
                     * Set reference wavelength as zero point
                     */

                    wavel = cpl_bivector_get_y(peaks_ident);
                    cpl_vector_subtract_scalar(wavel, refwave);

                    uorder = countLines / 2 - 1;
                    if (uorder > order)
                        uorder = order;

                    ids = mos_poly_wav2pix(peaks_ident, uorder, reject,
                                           2 * (uorder + 1), &usedLines,
                                           &ids_err, &peaks_ident_used_fit);

                    if (ids == NULL) {
                        cpl_error_reset();
                        cpl_bivector_delete(peaks_ident);
                        continue;
                    }

                    if (nlines)
                        nlines[i] = usedLines;
                    if (error)
                        error[i] = ids_err / sqrt(usedLines/(uorder + 1));

                    if (calibration) {
                        pixstart = cpl_polynomial_eval_1d(ids,
                           cpl_bivector_get_y_data(peaks_ident)[0], 
                           NULL);
                        pixend = cpl_polynomial_eval_1d(ids,
                           cpl_bivector_get_y_data(peaks_ident)[countLines-1], 
                           NULL);
                        extrapolation = (pixend - pixstart) / 5;
                        pixstart -= extrapolation;
                        pixend += extrapolation;
                        if (pixstart < 0)
                            pixstart = 0;
                        if (pixend > nx)
                            pixend = nx;

                        for (j = pixstart; j < pixend; j++) {
                            (idata + i*nx)[j] = mos_eval_dds(ids,
                                     firstLambda, lastLambda, refwave, j);
                        }
                    }

                    /*
                     * Residuals image
                     */

                    if (residuals || (restable && !(i%step))) {
                        if (restable && !(i%step)) {
                            lin = cpl_polynomial_new(1);
                            for (k = 0; k < 2; k++)
                                cpl_polynomial_set_coeff(lin, &k,
                                      cpl_polynomial_get_coeff(ids, &k));
                        }
                        for (j = 0; j < countLines; j++) {
                            pixe = cpl_bivector_get_x_data(peaks_ident)[j];
                            wave = cpl_bivector_get_y_data(peaks_ident)[j];
                            value = pixe
                                  - cpl_polynomial_eval_1d(ids, wave, NULL);
                            if (residuals) {
                                pixel = pixe + 0.5;
                                (ddata + i*nx)[pixel] = value;
                            }
                            if (restable && !(i%step)) {
                                for (k = 0; k < nref; k++) {
                                    if (fabs(line[k]-refwave-wave) < 0.1) {
                                        snprintf(name, MAX_COLNAME,
                                                 "r%d", i);
                                        cpl_table_set_double(restable, name,
                                                             k, value);
                                        value = pixe
                                              - cpl_polynomial_eval_1d(lin,
                                                          wave, NULL);
                                        snprintf(name, MAX_COLNAME,
                                                 "d%d", i);
                                        cpl_table_set_double(restable, name,
                                                             k, value);
                                        snprintf(name, MAX_COLNAME,
                                                 "p%d", i);
                                        cpl_table_set_double(restable, name,
                                                             k, pixe);
                                        break;
                                    }
                                }
                            }
                        }
                        if (restable && !(i%step)) {
                            cpl_polynomial_delete(lin);
                        }
/***
                        for (j = 0; j < countLines; j++) {
                            pixel = cpl_bivector_get_x_data(output)[j]
                                  + 0.5; 
                            (ddata + i*nx)[pixel] =
                            cpl_bivector_get_x_data(output)[j]
                          - cpl_polynomial_eval_1d(ids,
                            cpl_bivector_get_y_data(output)[j], 
                            NULL);
                        }
***/
                    }

                    if(detected_lines)
                    {
                        cpl_size oldsize = cpl_table_get_nrow(detected_lines); 
                        cpl_size nidentlines = cpl_bivector_get_size(peaks_ident); 
                        cpl_table_set_size(detected_lines, oldsize + nidentlines);
                        for(cpl_size idline = 0; idline < nidentlines ; ++idline)
                        {
                            double wave_ident = cpl_bivector_get_y_data(peaks_ident)[idline] + refwave;
                            double xpix_fit = cpl_polynomial_eval_1d(ids,
                                    wave_ident - refwave, NULL);
                            cpl_table_set_int(detected_lines, "slit_id",
                                 oldsize + idline, slit_id);
                            cpl_table_set_double(detected_lines, "xpos_rectified_iter",
                                 oldsize + idline, cpl_bivector_get_x_data(peaks_ident)[idline] + 1);
                            cpl_table_set_double(detected_lines, "ypos_rectified_iter",
                                 oldsize + idline, (double)i + 1);
                            cpl_table_set_double(detected_lines, "peak_flux",
                                 oldsize + idline, 
                                 sdata[i*nx+(int)(cpl_bivector_get_x_data(peaks_ident)[idline]+0.5)]);
                            cpl_table_set_double(detected_lines, "wave_ident_iter",
                                 oldsize + idline, wave_ident);
                            cpl_table_set_double(detected_lines, "xpos_fit_rect_wavecal",
                                 oldsize + idline, xpix_fit + 1);
                            cpl_table_set_int(detected_lines,
                                              "fit_used",
                                              oldsize + idline, 0);
                            for(cpl_size i_used = 0; i_used < cpl_bivector_get_size(peaks_ident_used_fit); ++i_used)
                            {
                                if(cpl_bivector_get_x_data(peaks_ident)[idline] == cpl_bivector_get_x_data(peaks_ident_used_fit)[i_used])
                                    cpl_table_set_int(detected_lines,
                                                     "fit_used",
                                                     oldsize + idline, 1);
                            }
                        }
                    }

                    
                    for (k = 0; k <= order; k++) {
                        if (k > uorder) {
                            cpl_table_set_double(idscoeff, clab[k], i, 0.0);
                        }
                        else {
                            cpl_table_set_double(idscoeff, clab[k], i,
                                        cpl_polynomial_get_coeff(ids, &k));
                        }
                    }

                    cpl_bivector_delete(peaks_ident);
                    cpl_polynomial_delete(ids);

                } /* End of loop "use ids as a first-guess" */

                cpl_polynomial_delete(fguess);
            }

            cpl_table_delete(coeff);

        }

        row_top = row_bot;

    } /* End of loop on slits */


    /*
     * At this point the idscoeff table has been filled with all the 
     * fits coefficients obtained for all the rows of the input image.
     * Now we apply these coefficients to resample the input image
     * at constant wavelength step.
     */

    for (i = 0; i < ny; i++) {

        missing = 0;
        ids = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            c = cpl_table_get_double(idscoeff, clab[k], i, &null);
            if (null) {
                cpl_polynomial_delete(ids);
                missing = 1;
                break;
            }
            cpl_polynomial_set_coeff(ids, &k, c);
        }
        if (missing)
            continue;

        pixstart = cpl_polynomial_eval_1d(ids, firstLambda - refwave, NULL);
        pixend = cpl_polynomial_eval_1d(ids, lastLambda - refwave, NULL);
        if (pixstart < 0)
            pixstart = 0;
        if (pixend > nx)
            pixend = nx;

        /*
         * Resampled image:
         */

        for (j = 0; j < nl; j++) {
            lambda = firstLambda + j * dispersion;
            fpixel = cpl_polynomial_eval_1d(ids, lambda - refwave, NULL);
            pixel = fpixel;
            if (pixel >= 0 && pixel < nx-1) {
                v1 = (sdata + i*nx)[pixel];
                v2 = (sdata + i*nx)[pixel+1];
                vi = v1 + (v2-v1)*(fpixel-pixel);
                (rdata + i*nl)[j] = vi;
            }
        }

        cpl_polynomial_delete(ids);
    }

    /* Set the invalid flag for the integer columns. Otherwise cpl_save_table()
       can crash (PIPE-6839) due to uninitialise values. For double columns this is 
       just fine. See documentation for cpl_table_save()  */
    if(detected_lines)
        cpl_table_fill_invalid_int(detected_lines, "fit_used", -1);

    return resampled;
}


/**
 * @brief
 *   Remap at constant wavelength step an image of rectified scientific spectra
 *
 * @param image       Image with rectified scientific spectra
 * @param refwave     Reference wavelength
 * @param firstLambda Wavelength of first pixel of resampled image
 * @param lastLambda  Wavelength of last pixel of resampled image
 * @param dispersion  Resampling step (wavelength units/pixel)
 * @param idscoeff    Table with IDS polynomials
 * @param flux        If zero, flux conservation correction is not applied.
 * 
 * @return Input image resampled at constant wavelength step.
 *
 * The input @em image is the image of the scientific spectra with
 * the spectral curvature removed, as produced by the function 
 * @c mos_spatial_calibration(). The table @em idscoeff contains 
 * the polynomial fits coefficients obtained for each input @em image 
 * row, as produced by the function @c mos_wavelength_calibration_final(),
 * and possibly modified by the function @c mos_wavelength_align(),
 * used to align the wavelength solution to available sky lines.
 * The input @em image is resampled at a constant wavelength step, 
 * @em dispersion, and is returned by this function. In case of error 
 * a @c NULL pointer is returned.
 */

cpl_image *mos_wavelength_calibration(cpl_image *image, double refwave, 
                                      double firstLambda, double lastLambda, 
                                      double dispersion, cpl_table *idscoeff, 
                                      int flux)
{

    const char *func = "mos_wavelength_calibration";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_image      *resampled;
    cpl_polynomial *ids;
    double          pixel_per_lambda;
    double          lambda;
    double          c;
    float          *sdata;
    float          *rdata;
    float           v0, v1, v2, v3, vi;
    float           fpixel;
    int             order;
    int             pixstart, pixend;
    int             nl, nx, ny, pixel;
    int             missing;
    int             null;
    int             i, j;
    cpl_size        k;


    if (dispersion <= 0.0) {
        cpl_msg_error(func, "The resampling step must be positive");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (lastLambda - firstLambda < dispersion) {
        cpl_msg_error(func, "Invalid spectral range: %.2f to %.2f", 
                      firstLambda, lastLambda);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (idscoeff == NULL) {
        cpl_msg_error(func, "An IDS coeff table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (image == NULL) {
        cpl_msg_error(func, "A scientific spectral image must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);
    sdata = cpl_image_get_data_float(image);

    nl = (lastLambda - firstLambda) / dispersion;
    resampled = cpl_image_new(nl, ny, CPL_TYPE_FLOAT);
    rdata = cpl_image_get_data_float(resampled);

    order = 0;
    while (order < 6 && cpl_table_has_column(idscoeff, clab[order]))
        ++order;
    --order;

    for (i = 0; i < ny; i++) {

        missing = 0;
        ids = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            c = cpl_table_get_double(idscoeff, clab[k], i, &null);
            if (null) {
                cpl_polynomial_delete(ids);
                missing = 1;
                break;
            }
            cpl_polynomial_set_coeff(ids, &k, c);
        }
        if (missing)
            continue;

        pixstart = cpl_polynomial_eval_1d(ids, firstLambda - refwave, NULL);
        pixend = cpl_polynomial_eval_1d(ids, lastLambda - refwave, NULL);
        if (pixstart < 0)
            pixstart = 0;
        if (pixend > nx)
            pixend = nx;

        /*
         * Resampled image:
         */

        for (j = 0; j < nl; j++) {
            lambda = firstLambda + j * dispersion;
            fpixel = cpl_polynomial_eval_1d(ids, lambda - refwave, 
                                            &pixel_per_lambda);

            /*
             * The local dispersion is 1 / pixel_per_lambda
             * and this factor is used for applying the flux
             * conservation correction (if requested).
             */

            pixel = fpixel;

         // if (dispersion * pixel_per_lambda < 2.0) {
            if (1) {  /* Old behaviour: this is safe. */

                /*
                 * In this case we just sample interpolating the
                 * signal of nearby pixels.
                 */

                //the wave calibration should be a monotonically increasing 
                //function. If the detivative is negative, we are outside of
                //the wavelength solution domain.
                if(pixel_per_lambda <= 0) 
                    vi = 0;
                else if (fpixel < 0)
                    vi = 0;
                else if (pixel >= 1 && pixel < nx-2) {
                   v0 = (sdata + i*nx)[pixel-1];
                   v1 = (sdata + i*nx)[pixel];
                   v2 = (sdata + i*nx)[pixel+1];
                   v3 = (sdata + i*nx)[pixel+2];
                   vi = (fpixel-pixel)*(fpixel-pixel)*(v3 - v2 - v1 + v0)
                      + (fpixel-pixel)*(3*v2 - v3 - v1 - v0)
                      + 2*v1;
                   vi /= 2;
                   if (v1 > v2) {
                       if (vi > v1) { 
                           vi = v1;
                       }
                       else if (vi < v2) {
                           vi = v2;
                       }
                   }
                   else {
                       if (vi > v2) { 
                           vi = v2;
                       }
                       else if (vi < v1) {
                           vi = v1;
                       }
                   }
                   if (flux)
                       vi *= dispersion * pixel_per_lambda;
               }
               else if (pixel >= 0 && pixel < nx-1) {
                   v1 = (sdata + i*nx)[pixel];
                   v2 = (sdata + i*nx)[pixel+1];
                   vi = v1 + (v2-v1)*(fpixel-pixel);
                   if (flux)
                       vi *= dispersion * pixel_per_lambda;
               }
               else
                   vi = 0;
               (rdata + i*nl)[j] = vi;
           }
           else {

               /*
                * Here instead we integrate the pixel values in
                * the interval centered at the interpolation point.
                * This interval is long dispersion * pixel_per_lambda
                * of the original pixels, and is centered at fpixel.
                * So it starts at fpixel - dispersion * pixel_per_lambda / 2,
                * and it ends at fpixel + dispersion * pixel_per_lambda / 2.
                */

               double spos = fpixel - dispersion * pixel_per_lambda / 2;
               double epos = fpixel + dispersion * pixel_per_lambda / 2;

               /*
                * Brutal sum over all involved pixels
                */

               int spix = spos;
               int epix = epos + 1;

               if (spix < 0)
                   spix = 0;

               if (epix > nx)
                   epix = nx;

               vi = 0.0;
               for (k = spix; k < epix; k++) {
                   if (pixel >= 0 && pixel < nx) {
                       vi += (sdata + i*nx)[k];
                   }
               }

               /*
                * Correct integrated flux by true length
                * of interval. This is clearly an approximation,
                * but it's good enough if the PSF is much larger
                * than the original pix.
                */

               vi *= dispersion * pixel_per_lambda / (epix - spix);

                /*
                 * Flux conservation is a geometric factor that is applied
                 * always in the same way...
                 */

               if (flux)
                   vi *= dispersion * pixel_per_lambda;

               (rdata + i*nl)[j] = vi;
            }
        }

        cpl_polynomial_delete(ids);
    }

    return resampled;
}


/**
 * @brief
 *   Modify the input wavelength solution to match reference sky lines.
 *   
 * @param image       Image with rectified scientific spectra
 * @param slits       Slits positions table
 * @param refwave     Reference wavelength 
 * @param firstLambda Wavelength range start
 * @param lastLambda  Wavelength range end
 * @param idscoeff    Table with IDS polynomials
 * @param skylines    Vector with wavelengths of reference skylines
 * @param highres     1 = high resolution data, 0 = low resolution data
 * @param order       Order of sky lines offsets fitting polynomial
 * @param calibration Rectified wavelength calibration image
 * @param sradius     Search radius for sky lines
 *
 * @return Table with measured sky lines offsets (in pixel)
 *
 * The input @em image is the image of the scientific spectra with
 * the spectral curvature removed, as produced by the function
 * @c mos_spatial_calibration(). The input @em slits table should be
 * the same used and processed by the same function, containing a
 * column named "position" listing the position of the individual
 * spectra in the rectified image. The table @em idscoeff contains
 * the polynomial fits coefficients obtained for each input @em image
 * row, as produced by the function @c mos_wavelength_calibration_final().
 * The vector @em skylines should list a number of reference sky-lines
 * wavelengths. The calibrating polynomials will be used to get the
 * expected positions of the reference skylines for each slit, for
 * measuring their offset from this position as a function of wavelength. 
 * If the @em skylines vector is not provided, an internal list of
 * sky lines is used instead. Only in this case the argument @em highres
 * is checked in order to select a sky lines list for high or low
 * resolution data (by high resolution data we mean here R > 700). 
 * In all cases, only the sky lines within
 * the specified range @em firstLambda @em lastLambda will be used.
 * If the rectified wavelengths @em calibration image produced by the
 * function @c mos_wavelength_calibration_final() is provided, the 
 * mapped wavelengths will be upgraded to reflect the alignment of 
 * the old solution to the position of the detected sky lines.
 * Subsequently, such map may be transferred to the CCD image
 * using the function @c mos_map_wavelengths().
 * 
 * More in detail: for each slit, a "sky" spectrum is derived by computing 
 * the median slit spectrum, reducing in this way possible contamination 
 * by bright scientific object(s) contained in the slit. All the IDS
 * polynomials modelling the slit are used to determine the sky lines 
 * offsets: with N slit rows, N independent offsets for each reference 
 * sky line are found; a new table is created, with columns containing
 * the median offsets obtained for each reference sky line and for each
 * slit. The columns are named after the "slit_id" listed in the input 
 * @em slits table, preceded by "offset". For instance, the median
 * offsets for all sky lines observed on a slit with id = 123 will be 
 * written to a column named "offset123". A column "wave" listing the
 * reference wavelengths used will be also included to the table.
 *
 * For each slit, the median offsets are fitted by a polynomial relation 
 * that is then summed to all the IDS polynomials pertaining to that slit,
 * modifying in this way the input @em idscoeff table. It is of course 
 * advisable not to go beyond a fitting polynomial of 0th or 1st order, 
 * unless many and well distributed reference sky lines are available 
 * within the specified range. Note that this way of proceeding is applied 
 * to avoid destroying the information about irregularities in the slit 
 * ("slit geometry").
 */

cpl_table *mos_wavelength_align(cpl_image *image, cpl_table *slits, 
                                double refwave, double firstLambda, 
                                double lastLambda, cpl_table *idscoeff,
                                cpl_vector *skylines, int highres, int order,
                                cpl_image *calibration, int sradius)
{
    const char *func = "mos_wavelength_align";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    double         *line;
    double         *data;
    double          expPos, offset;
    double          c;
    double          lambda1, lambda2;
    double          rms;
    float           pos;
    float          *sdata;
    float          *cdata;
    int            *idata;
    int             startPos, endPos;
    int             window = 2*sradius + 1;
    int             nlines;
    int             nslits;
    int             npoints;
    int             nrows;
    int             nx, ny;
    int             xlow, ylow, xhig, yhig;
    int             idsorder, uorder;
    int            *slit_id;
    int            *position;
    int            *length;
    int             missing;
    int             null;
    int             i;
    cpl_size        j, k;

    char            offname[MAX_COLNAME];
    char            name[MAX_COLNAME];

    cpl_polynomial *ids;
    cpl_polynomial *polycorr;
    cpl_image      *exslit;
    cpl_image      *sky;
    cpl_table      *offsets;
    cpl_table      *dummy;
    cpl_vector     *wave;
    cpl_vector     *offs;
    

    if (idscoeff == NULL) {
        cpl_msg_error(func, "An IDS coeff table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (image == NULL) {
        cpl_msg_error(func, "A scientific spectral image must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (slits == NULL) {
        cpl_msg_error(func, "A slit position table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (skylines) {
        line = cpl_vector_get_data(skylines);
        nlines = cpl_vector_get_size(skylines);
    }
    else {
        cpl_msg_warning(func, "A catalog of sky lines wavelengths was not "
                        "given: using internal list of reference sky lines");
        if (highres) {
           line = default_lines_hi;
           nlines = sizeof(default_lines_hi) / sizeof(double);
        }
        else {
           line = default_lines_lo;
           nlines = sizeof(default_lines_lo) / sizeof(double);
        }
    }

    if (calibration)
        cdata = cpl_image_get_data(calibration);

    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);

    nslits   = cpl_table_get_nrow(slits);
    slit_id  = cpl_table_get_data_int(slits, "slit_id");
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");


    /*
     * Define the output table of offsets
     */

    nrows = 0;
    for (i = 0; i < nlines; i++)
        if (line[i] > firstLambda && line[i] < lastLambda)
            nrows++;

    offsets = cpl_table_new(nrows);
    cpl_table_new_column(offsets, "wave", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(offsets, "wave", "Angstrom");

    nrows = 0;
    for (i = 0; i < nlines; i++) {
        if (line[i] > firstLambda && line[i] < lastLambda) {
            cpl_table_set_double(offsets, "wave", nrows, line[i]);
            nrows++;
        }
    }

    /*
     * Here "line" is made to point to the new list of selected wavelengths
     */

    line = cpl_table_get_data_double(offsets, "wave");
    nlines = nrows;

    idsorder = 0;
    while (idsorder < 6 && cpl_table_has_column(idscoeff, clab[idsorder]))
        ++idsorder;
    --idsorder;

    xlow = 1;
    xhig = nx;
    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        snprintf(offname, MAX_COLNAME, "offset%d", slit_id[i]);
        cpl_table_new_column(offsets, offname, CPL_TYPE_DOUBLE);

        /* 
         * Define the extraction boundaries. We DON'T write:
         *
         * ylow = position[i];
         * yhig = ylow + length[i];
         *
         * because the cpl_image pixels are counted from 1, and because in
         * cpl_image_extract() the coordinates of the last pixel are inclusive.
         */

        ylow = position[i] + 1;
        yhig = ylow + length[i] - 1;

        exslit = cpl_image_extract(image, xlow, ylow, xhig, yhig);
        sky    = cpl_image_collapse_median_create(exslit, 0, 0, 1);
        sdata  = cpl_image_get_data(sky);

        cpl_image_delete(exslit);

        /* 
         * Return here to a decent way of counting pixels (i.e., starting
         * from 0)
         */
         
        ylow--;

        /*
         * Allocate a dummy table for collecting all the offsets
         * for all the lines: this is only needed for the computation
         * of the median offset for each sky line
         */

        dummy = cpl_table_new(yhig - ylow);
        for (j = 0; j < nlines; j++) {
            snprintf(name, MAX_COLNAME, "%"CPL_SIZE_FORMAT, j);
            cpl_table_new_column(dummy, name, CPL_TYPE_DOUBLE);
        }

        for (j = ylow; j < yhig; j++) {

            /*
             * Get the IDS polynomial for the current slit row
             */

            missing = 0;
            ids = cpl_polynomial_new(1);
            for (k = 0; k <= idsorder; k++) {
                c = cpl_table_get_double(idscoeff, clab[k], j, &null);
                if (null) {
                    cpl_polynomial_delete(ids);
                    missing = 1;
                    break;
                }
                cpl_polynomial_set_coeff(ids, &k, c);
            }
            if (missing)
                continue;

            for (k = 0; k < nlines; k++) {
                expPos = cpl_polynomial_eval_1d(ids, line[k] - refwave, NULL);
                startPos = expPos - sradius;
                endPos   = startPos + window;
                if (startPos < 0 || endPos >= nx)
                    continue;
           
                if (0 == peakPosition(sdata + startPos, window, &pos, 1)) {
                    pos += startPos;
                    offset = pos - expPos;
                    snprintf(name, MAX_COLNAME, "%"CPL_SIZE_FORMAT, k);
                    cpl_table_set_double(dummy, name, j - ylow, offset);
                }
            }

            cpl_polynomial_delete(ids);
        }

        cpl_image_delete(sky);

        for (j = 0; j < nlines; j++) {
            snprintf(name, MAX_COLNAME, "%"CPL_SIZE_FORMAT, j);
            if (cpl_table_has_valid(dummy, name)) {
                offset = cpl_table_get_column_median(dummy, name);
                cpl_table_set_double(offsets, offname, j, offset);
            }
        }

        cpl_table_delete(dummy);

    }


    /*
     * In the following the input idscoeff table is modified by simply
     * adding the coefficients of the polynomial used to fit the sky
     * line residuals to the coefficients of the IDS polynomials.
     */

    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        snprintf(offname, MAX_COLNAME, "offset%d", slit_id[i]);

        /*
         * In the following, the "dummy" table is just a tool for
         * eliminating invalid points from the vectors to be fitted.
         */

        dummy = cpl_table_new(nlines);
        cpl_table_duplicate_column(dummy, "wave", offsets, "wave");
        cpl_table_duplicate_column(dummy, "offset", offsets, offname);

        npoints = nlines - cpl_table_count_invalid(dummy, "offset");
        if (npoints == 0) {
            cpl_msg_warning(func, "No sky lines alignment was possible "
                            "for slit ID=%d: no sky line found", slit_id[i]);
            cpl_table_delete(dummy);
            continue;
        }

        uorder = order;
        if (npoints <= uorder) {
            uorder = npoints - 1;
            if (uorder) {
                cpl_msg_warning(func, "Just %d sky lines detected for slit "
                                "ID=%d, while a polynomial order %d was "
                                "requested. Using polynomial order %d for "
                                "this slit!", npoints, slit_id[i], order, 
                                uorder);
            }
            else {
                cpl_msg_warning(func, "Just %d sky lines detected for slit "
                                "ID=%d, while a polynomial order %d was "
                                "requested. Computing a median offset for "
                                "this slit!", npoints, slit_id[i], order);
            }
        }

        cpl_table_erase_invalid(dummy);

        if (uorder > 1) {

            /*
             * Model offsets with polynomial fitting
             */

            wave = cpl_vector_wrap(npoints,
                                   cpl_table_get_data_double(dummy, "wave"));
            offs = cpl_vector_wrap(npoints,
                                   cpl_table_get_data_double(dummy, "offset"));

            /*
             * Set reference wavelength as zero point
             */

            cpl_vector_subtract_scalar(wave, refwave);

            polycorr = cpl_polynomial_fit_1d_create(wave, offs, uorder, &rms);

            rms = sqrt(rms * (uorder + 1) / npoints);

            cpl_vector_unwrap(wave);
            cpl_vector_unwrap(offs);
            cpl_table_delete(dummy);

            /*
             * Now correct the coefficients of the corresponding IDS
             * polynomials related to this slit:
             */

            ylow = position[i];
            yhig = ylow + length[i];

            for (j = 0; j <= uorder; j++) {
                data = cpl_table_get_data_double(idscoeff, clab[j]);
                c = cpl_polynomial_get_coeff(polycorr, &j);
                for (k = ylow; k < yhig; k++)
                    data[k] += c;
            }

            data = cpl_table_get_data_double(idscoeff, "error");
            for (k = ylow; k < yhig; k++)
                 data[k] = sqrt(data[k]*data[k] + rms*rms);

            idata = cpl_table_get_data_int(idscoeff, "nlines");
            for (k = ylow; k < yhig; k++)
                 idata[k] = npoints;

            /*
             * If a wavelengths map was provided, correct it to keep
             * into account the alignment to skylines:
             */

            if (calibration) {
                for (j = ylow; j < yhig; j++) {
                    for (k = 1; k < nx; k++) {
                        lambda1 = cdata[k - 1 + j*nx];
                        lambda2 = cdata[k + j*nx];
                        if (lambda1 < 1.0 || lambda2 < 1.0)
                            continue;
                        offset = cpl_polynomial_eval_1d(polycorr, 
                                                        lambda1-refwave, NULL);
                        cdata[k - 1 + j*nx] -= offset * (lambda2-lambda1);
                    }
                }
            }
    
            cpl_polynomial_delete(polycorr);
        }
        else if (uorder == 1) {

            /*
             * Model offsets with robust linear fitting
             */

            double        q, m;
            cpl_bivector *list;


            wave = cpl_vector_wrap(npoints,
                                   cpl_table_get_data_double(dummy, "wave"));
            offs = cpl_vector_wrap(npoints,
                                   cpl_table_get_data_double(dummy, "offset"));

            list = cpl_bivector_wrap_vectors(wave, offs);

            /*
             * Set reference wavelength as zero point
             */

            cpl_vector_subtract_scalar(wave, refwave);

            robustLinearFit(list, &q, &m, &rms);

            rms = sqrt(rms * (uorder + 1) / npoints);

            cpl_bivector_unwrap_vectors(list);
            cpl_vector_unwrap(wave);
            cpl_vector_unwrap(offs);
            cpl_table_delete(dummy);

            /*
             * Now correct the coefficients of the corresponding IDS
             * polynomials related to this slit:
             */

            ylow = position[i];
            yhig = ylow + length[i];

            for (j = 0; j <= uorder; j++) {
                data = cpl_table_get_data_double(idscoeff, clab[j]);
                if (j)
                    c = m;
                else
                    c = q;
                for (k = ylow; k < yhig; k++)
                    data[k] += c;
            }

            data = cpl_table_get_data_double(idscoeff, "error");
            for (k = ylow; k < yhig; k++)
                 data[k] = sqrt(data[k]*data[k] + rms*rms);

            idata = cpl_table_get_data_int(idscoeff, "nlines");
            for (k = ylow; k < yhig; k++)
                 idata[k] = npoints;

            /*
             * If a wavelengths map was provided, correct it to keep
             * into account the alignment to skylines:
             */

            if (calibration) {
                for (j = ylow; j < yhig; j++) {
                    for (k = 1; k < nx; k++) {
                        lambda1 = cdata[k - 1 + j*nx];
                        lambda2 = cdata[k + j*nx];
                        if (lambda1 < 1.0 || lambda2 < 1.0)
                            continue;
                        offset = q + m*(lambda1-refwave);
                        cdata[k - 1 + j*nx] -= offset * (lambda2-lambda1);
                    }
                }
            }
        }
        else {

            /*
             * Just compute median offset
             */

            offs = cpl_vector_wrap(npoints,
                                   cpl_table_get_data_double(dummy, "offset"));

            offset = cpl_vector_get_median_const(offs);

            if (npoints > 1)
                rms = cpl_table_get_column_stdev(dummy, "offset");
            else
                rms = 0.0;

            rms /= sqrt(npoints);

            cpl_vector_unwrap(offs);
            cpl_table_delete(dummy);

            /*
             * Now correct the constant term of the corresponding IDS
             * polynomials related to this slit:
             */

            ylow = position[i];
            yhig = ylow + length[i];

            data = cpl_table_get_data_double(idscoeff, clab[0]);
            for (k = ylow; k < yhig; k++)
                data[k] += offset;

            data = cpl_table_get_data_double(idscoeff, "error");
            for (k = ylow; k < yhig; k++)
                 data[k] = sqrt(data[k]*data[k] + rms*rms);

            idata = cpl_table_get_data_int(idscoeff, "nlines");
            for (k = ylow; k < yhig; k++)
                 idata[k] = npoints;

            /*
             * If a wavelengths map was provided, correct it to keep
             * into account the alignment to skylines. Note that 
             * the offset must be converted from pixels to wavelengths.
             */

            if (calibration) {
                for (j = ylow; j < yhig; j++) {
                    for (k = 1; k < nx; k++) {
                        lambda1 = cdata[k - 1 + j*nx];
                        lambda2 = cdata[k + j*nx];
                        if (lambda1 < 1.0 || lambda2 < 1.0)
                            continue; 
                        cdata[k - 1 + j*nx] -= offset * (lambda2-lambda1);
                    }
                }
            }
        }
    }

    return offsets;

}


/**
 * @brief
 *   Modify the input wavelength solution to match reference sky lines (LSS).
 *   
 * @param image       Image with rectified scientific spectra
 * @param refwave     Reference wavelength 
 * @param firstLambda Wavelength range start
 * @param lastLambda  Wavelength range end
 * @param idscoeff    Table with IDS polynomials
 * @param skylines    Vector with wavelengths of reference skylines
 * @param highres     1 = High resolution, 0 = Low resolution
 * @param order       Order of sky lines offsets fitting polynomial
 * @param calibration Rectified wavelength calibration image
 * @param sradius     Search radius for sky lines
 *
 * @return Table with measured sky lines offsets (in pixel)
 *
 * The input @em image is the image of the scientific spectra with
 * the spectral curvature removed, as produced by the function
 * @c mos_spatial_calibration(). The table @em idscoeff contains
 * the polynomial fits coefficients obtained for each input @em image
 * row, as produced by the function @c mos_wavelength_calibration_final().
 * The vector @em skylines should list a number of reference sky-lines
 * wavelengths. The calibrating polynomials will be used to get the
 * expected positions of the reference skylines for each slit, for
 * measuring their offset from this position as a function of wavelength. 
 * If the @em skylines vector is not provided, an internal list of
 * sky lines is used instead. Only in this case the argument @em highres
 * is checked in order to select a sky lines list for high or low
 * resolution data (by high resolution data we mean here R > 700).
 * In all cases, only the sky lines within
 * the specified range @em firstLambda @em lastLambda will be used.
 * If the rectified wavelengths @em calibration image produced by the
 * function @c mos_wavelength_calibration_final() is provided, the 
 * mapped wavelengths will be upgraded to reflect the alignment of 
 * the old solution to the position of the detected sky lines.
 * Subsequently, such map may be transferred to the CCD image
 * using the function @c mos_map_wavelengths().
 * 
 * More in detail: a "sky" spectrum is derived by computing the
 * median spectrum, reducing in this way possible contamination 
 * by bright scientific object(s) contained in the image. All the IDS
 * polynomials modelling the image are used to determine the sky lines 
 * offsets: with N image rows, N independent offsets for each reference 
 * sky line are found; a new table is created, with a column containing
 * the median offsets obtained for each reference sky line and for each
 * image row. This column is named "offset". A column "wave" listing the
 * reference wavelengths used will be also included to the table.
 *
 * For each image row, the measured offsets are fitted by a polynomial 
 * relation of given @em order. The coefficients of all the obtained 
 * polynomial relations are then modeled as a function of the image row,
 * and their modeled values are then summed to all the IDS polynomials 
 * modifying in this way the input @em idscoeff table. It is of course 
 * advisable not to go beyond a fitting polynomial of 0th or 1st order, 
 * unless many and well distributed reference sky lines are available 
 * within the specified range. Note that this way of proceeding is applied 
 * to avoid destroying the information about irregularities in the slit 
 * ("slit geometry").
 */

cpl_table *mos_wavelength_align_lss(cpl_image *image, double refwave, 
                                    double firstLambda, double lastLambda, 
                                    cpl_table *idscoeff, cpl_vector *skylines, 
                                    int highres, int order, 
                                    cpl_image *calibration, int sradius)
{
    const char *func = "mos_wavelength_align_lss";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    double         *line;
    double         *data;
    double         *wdata;
    double         *odata;
    double          expPos, offset;
    double          c;
    double          lambda1, lambda2;
    double          rms;
    float           pos;
    float          *sdata;
    float          *cdata;
    int            *idata;
    int             startPos, endPos;
    int             window = 2*sradius + 1;
    int             nlines;
    int             npoints;
    int             nrows;
    int             nx, ny;
    int             idsorder, uorder;
    int             missing;
    int             i;
    cpl_size        j, k;

    char            name[MAX_COLNAME];
    char            fname[MAX_COLNAME];

    cpl_polynomial *ids;
    cpl_polynomial *polycorr;
    cpl_table      *offsets;
    cpl_table      *fittable;
    cpl_table      *dummy;
    cpl_vector     *wave;
    cpl_vector     *offs;
    cpl_vector     *row;
    

    if (idscoeff == NULL) {
        cpl_msg_error(func, "An IDS coeff table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (image == NULL) {
        cpl_msg_error(func, "A scientific spectral image must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (skylines) {
        line = cpl_vector_get_data(skylines);
        nlines = cpl_vector_get_size(skylines);
    }
    else {
        cpl_msg_warning(func, "A catalog of sky lines wavelengths was not "
                        "given: using internal list of reference sky lines");
        if (highres) {
           line = default_lines_hi;
           nlines = sizeof(default_lines_hi) / sizeof(double);
        }
        else {
           line = default_lines_lo;
           nlines = sizeof(default_lines_lo) / sizeof(double);
        }
    }

    if (calibration)
        cdata = cpl_image_get_data(calibration);

    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);

    sdata = cpl_image_get_data(image);

    if (ny != cpl_table_get_nrow(idscoeff)) {
        cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
        return NULL;
    }
    

    /*FIXME: This is a remnant of the adaptation of the function
     * mos_wavelength_align(), where an offset table was created.
     * I leave it here because I am in a hurry, it is just used to
     * hold the list of selected sky lines.
     *
     * Define table of wavelengths
     */

    nrows = 0;
    for (i = 0; i < nlines; i++)
        if (line[i] > firstLambda && line[i] < lastLambda)
            nrows++;

    offsets = cpl_table_new(nrows);
    cpl_table_new_column(offsets, "wave", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(offsets, "wave", "Angstrom");

    nrows = 0;
    for (i = 0; i < nlines; i++) {
        if (line[i] > firstLambda && line[i] < lastLambda) {
            cpl_table_set_double(offsets, "wave", nrows, line[i]);
            nrows++;
        }
    }

    /*
     * Here "line" is made to point to the new list of selected wavelengths
     */

    line = cpl_table_get_data_double(offsets, "wave");
    nlines = nrows;

    idsorder = 0;
    while (idsorder < 6 && cpl_table_has_column(idscoeff, clab[idsorder]))
        ++idsorder;
    --idsorder;


    /*
     * Allocate a dummy table for collecting all the offsets
     * for all the lines
     */

    dummy = cpl_table_new(ny);
    for (j = 0; j < nlines; j++) {
        snprintf(name, MAX_COLNAME, "off_%d", (int)line[j]);
        snprintf(fname, MAX_COLNAME, "fit_%d", (int)line[j]);
        cpl_table_new_column(dummy, name, CPL_TYPE_DOUBLE);
        cpl_table_new_column(dummy, fname, CPL_TYPE_DOUBLE);
    }

    for (j = 0; j < ny; j++, sdata += nx) {

        /*
         * Get the IDS polynomial for the current slit row
         */

        missing = 0;
        ids = cpl_polynomial_new(1);
        for (k = 0; k <= idsorder; k++) {
            c = cpl_table_get_double(idscoeff, clab[k], j, &missing);
            if (missing) {
                cpl_polynomial_delete(ids);
                break;
            }
            cpl_polynomial_set_coeff(ids, &k, c);
        }
        if (missing)
            continue;

        for (k = 0; k < nlines; k++) {
            expPos = cpl_polynomial_eval_1d(ids, line[k] - refwave, NULL);
            startPos = expPos - sradius;
            endPos   = startPos + window;
            if (startPos < 0 || endPos >= nx)
                continue;
           
            if (0 == peakPosition(sdata + startPos, window, &pos, 1)) {
                pos += startPos;
                offset = pos - expPos;
                snprintf(name, MAX_COLNAME, "off_%d", (int)line[k]);
                cpl_table_set_double(dummy, name, j, offset);
            }
        }

        cpl_polynomial_delete(ids);
    }


    /*
     * At this point for each sky line we model its offset along
     * the image rows using a robust linear fitting
     */

    for (j = 0; j < nlines; j++) {
        snprintf(name, MAX_COLNAME, "off_%d", (int)line[j]);
        snprintf(fname, MAX_COLNAME, "fit_%d", (int)line[j]);
        if (cpl_table_has_valid(dummy, name)) {

            /*
             * In the following, the "fittable" is just a tool for
             * eliminating invalid points from the vectors to be fitted.
             */

            double        q, m;
            cpl_bivector *list;

            fittable = cpl_table_new(ny);
            cpl_table_new_column(fittable, "row", CPL_TYPE_DOUBLE);
            cpl_table_set_column_unit(fittable, "row", "pixel");
            for (k = 0; k < ny; k++)
                 cpl_table_set_double(fittable, "row", k, k);
            cpl_table_duplicate_column(fittable, "offset", dummy, name);
            npoints = ny - cpl_table_count_invalid(fittable, "offset");
            cpl_table_erase_invalid(fittable);
            row = cpl_vector_wrap(npoints,
                               cpl_table_get_data_double(fittable, "row"));
            offs = cpl_vector_wrap(npoints,
                               cpl_table_get_data_double(fittable, "offset"));
            list = cpl_bivector_wrap_vectors(row, offs);
            robustLinearFit(list, &q, &m, &rms);
            cpl_bivector_unwrap_vectors(list);
            cpl_vector_unwrap(row);
            cpl_vector_unwrap(offs);
            cpl_table_delete(fittable);
            for (k = 0; k < ny; k++)
                 cpl_table_set_double(dummy, fname, k, q + m*k);
        }
    }


    /*
     * Now each dummy table row consists of a sequence of offsets,
     * one for each wavelength. A table row corresponds to an image row.
     * We must fit a polynomial to each one of these rows, in order to
     * express the offsets as a function of wavelength. The obtained 
     * polynomial coefficients are used to correct the IDS coefficients.
     */

    for (i = 0; i < ny; i++) {

        if (!cpl_table_is_valid(idscoeff, clab[0], i))
            continue;

        npoints = 0;
        for (j = 0; j < nlines; j++) {
            snprintf(name, MAX_COLNAME, "fit_%d", (int)line[j]);
            if (cpl_table_is_valid(dummy, name, i))
                npoints++;
        }

        if (npoints == 0)
            continue;

        uorder = order;
        if (npoints <= uorder)
            uorder = npoints - 1;

        if (uorder > 1) {

            /*
             * Model offsets with polynomial fitting
             */

            wave = cpl_vector_new(npoints);
            wdata = cpl_vector_get_data(wave);
            offs = cpl_vector_new(npoints);
            odata = cpl_vector_get_data(offs);

            npoints = 0;
            for (j = 0; j < nlines; j++) {
                snprintf(name, MAX_COLNAME, "fit_%d", (int)line[j]);
                if (cpl_table_is_valid(dummy, name, i)) {
                    wdata[npoints] = line[j] - refwave;
                    odata[npoints] = cpl_table_get_double(dummy, name, i, NULL);
                    npoints++;
                }
            }

            polycorr = cpl_polynomial_fit_1d_create(wave, offs, uorder, &rms);

            rms = sqrt(rms * (uorder + 1) / npoints);

            cpl_vector_delete(wave);
            cpl_vector_delete(offs);

            /*
             * Now correct the coefficients of the corresponding IDS
             * polynomials related to this slit:
             */

            for (j = 0; j <= uorder; j++) {
                data = cpl_table_get_data_double(idscoeff, clab[j]);
                c = cpl_polynomial_get_coeff(polycorr, &j);
                data[i] += c;
            }

            data = cpl_table_get_data_double(idscoeff, "error");
            data[i] = sqrt(data[i]*data[i] + rms*rms);

            idata = cpl_table_get_data_int(idscoeff, "nlines");
            idata[i] = npoints;

            /*
             * If a wavelengths map was provided, correct it to keep
             * into account the alignment to skylines:
             */

            if (calibration) {
                for (k = 1; k < nx; k++) {
                    lambda1 = cdata[k - 1 + i*nx];
                    lambda2 = cdata[k + i*nx];
                    if (lambda1 < 1.0 || lambda2 < 1.0)
                        continue;
                    offset = cpl_polynomial_eval_1d(polycorr,
                                                    lambda1-refwave, NULL);
                    cdata[k - 1 + i*nx] -= offset * (lambda2-lambda1);
                }
            }

            cpl_polynomial_delete(polycorr);

        }
        else if (uorder == 1) {

            /*
             * Model offsets with robust linear fitting
             */

            cpl_bivector *list;
            double        q, m;

            wave = cpl_vector_new(npoints);
            wdata = cpl_vector_get_data(wave);
            offs = cpl_vector_new(npoints);
            odata = cpl_vector_get_data(offs);

            npoints = 0;
            for (j = 0; j < nlines; j++) {
                snprintf(name, MAX_COLNAME, "fit_%d", (int)line[j]);
                if (cpl_table_is_valid(dummy, name, i)) {
                    wdata[npoints] = line[j] - refwave;
                    odata[npoints] = cpl_table_get_double(dummy, name, i, NULL);
                    npoints++;
                }
            }

            list = cpl_bivector_wrap_vectors(wave, offs);
            robustLinearFit(list, &q, &m, &rms);

            rms = sqrt(rms * (uorder + 1) / npoints);

            cpl_bivector_unwrap_vectors(list);
            cpl_vector_delete(wave);
            cpl_vector_delete(offs);

            /*
             * Now correct the coefficients of the corresponding IDS
             * polynomials related to this row:
             */

            for (j = 0; j <= uorder; j++) {
                data = cpl_table_get_data_double(idscoeff, clab[j]);
                if (j)
                    c = m;
                else
                    c = q;
                data[i] += c;
            }

            data = cpl_table_get_data_double(idscoeff, "error");
            data[i] = sqrt(data[i]*data[i] + rms*rms);

            idata = cpl_table_get_data_int(idscoeff, "nlines");
            idata[i] = npoints;

            /*
             * If a wavelengths map was provided, correct it to keep
             * into account the alignment to skylines:
             */

            if (calibration) {
                for (k = 1; k < nx; k++) {
                    lambda1 = cdata[k - 1 + i*nx];
                    lambda2 = cdata[k + i*nx];
                    if (lambda1 < 1.0 || lambda2 < 1.0)
                        continue;
                    offset = q + m*(lambda1-refwave);
                    cdata[k - 1 + i*nx] -= offset * (lambda2-lambda1);
                }
            }
        }
        else {

            /*
             * Just compute median offset
             */

            offs = cpl_vector_new(npoints);
            odata = cpl_vector_get_data(offs);

            npoints = 0;
            for (j = 0; j < nlines; j++) {
                snprintf(name, MAX_COLNAME, "fit_%d", (int)line[j]);
                if (cpl_table_is_valid(dummy, name, i)) {
                    odata[npoints] = cpl_table_get_double(dummy, name, i, NULL);
                    npoints++;
                }
            }

            offset = cpl_vector_get_median_const(offs);

            if (npoints > 1) {
                rms = cpl_vector_get_stdev(offs);
            }
            else if (npoints == 1) {
                snprintf(name, MAX_COLNAME, "off_%d", (int)line[0]);
                if (cpl_table_has_valid(dummy, name)) {
                    rms = cpl_table_get_column_stdev(dummy, name);
                    rms /= sqrt(ny - cpl_table_count_invalid(dummy, name));
                }
                else {
                    rms = 0.0;
                }
            }
            else {
                rms = 0.0;
            }

            rms /= sqrt(npoints);

            cpl_vector_delete(offs);

            /*
             * Now correct the constant term of the corresponding IDS
             * polynomials related to this slit:
             */

            data = cpl_table_get_data_double(idscoeff, clab[0]);
            data[i] += offset;

            data = cpl_table_get_data_double(idscoeff, "error");
            data[i] = sqrt(data[i]*data[i] + rms*rms);

            idata = cpl_table_get_data_int(idscoeff, "nlines");
            idata[i] = npoints;

            /*
             * If a wavelengths map was provided, correct it to keep
             * into account the alignment to skylines. Note that
             * the offset must be converted from pixels to wavelengths.
             */

            if (calibration) {
                for (k = 1; k < nx; k++) {
                    lambda1 = cdata[k - 1 + i*nx];
                    lambda2 = cdata[k + i*nx];
                    if (lambda1 < 1.0 || lambda2 < 1.0)
                        continue;
                    cdata[k - 1 + i*nx] -= offset * (lambda2-lambda1);
                }
            }
        }
    }

    missing = 1;
    for (j = 0; j < nlines; j++) {
        snprintf(name, MAX_COLNAME, "off_%d", (int)line[j]);
        if (cpl_table_has_valid(dummy, name)) {
            missing = 0;
            offset = cpl_table_get_column_median(dummy, name);
            cpl_msg_info(func, "Median offset for %.3f: %.3f pixel",
                         line[j], offset);
        }
        else {
            cpl_msg_info(func, 
                         "Median offset for %.2f: not available", line[j]);
        }
    }

    cpl_table_delete(offsets);

    if (missing) {
        cpl_table_delete(dummy);
        dummy = NULL;
    }

    return dummy;

}


/**
 * @brief
 *   Estimate the spectral distortion modeling goodness.
 *
 * @param rectified  Calibrated calibration image
 * @param lines      Reference wavelengths (line catalog)
 * @param wavestart  Wavelength of bluest (left) pixel in input image
 * @param dispersion Angstrom per pixel of input image
 * @param radius     Search radius (in pixels of input image)
 * @param highres    1 = high resolution data, 0 = low resolution data
 *
 * @return Mean RMS of residuals in pixels
 *
 * The input @em rectified image is the product of either the functions 
 * @c mos_wavelength_calibration_raw() and @c mos_wavelength_calibration_final()
 * The RMS of the residuals between the expected positions of the arc 
 * lamp lines and the actual position of a detected peak is computed. 
 * The peak is searched within the search radius specified by the last 
 * argument, which typically would depend on the expected lines FWHM. 
 * This function prints also the computed RMS for each wavelength
 * listed in the line catalog to screen. If the @em lines vector is 
 * not provided, an internal list of sky lines wavelengths is used 
 * instead. Only in this case the argument @em highres is checked in 
 * order to select a sky lines list for high or low resolution data 
 * (by high resolution data we mean here R > 700).
 */

double mos_distortions_rms(cpl_image *rectified, cpl_vector *lines, 
                           double wavestart, double dispersion, int radius,
                           int highres)
{

    const char *func = "mos_distortions_rms";

    int xlen;
    int ylen;
    int numLines;
    int cpix, npix, nzero;
    int sp, ep;
    int i, j, k;
    int npeaks, allPeaks;

    float *profile;
    float  peak, expectPeak, offset;
    double lambda;

    double  average;
    double  rms, oneRms;

    float  *sdata;
    double *wdata;

  
    xlen = cpl_image_get_size_x(rectified);
    ylen = cpl_image_get_size_y(rectified);
    sdata = cpl_image_get_data(rectified);

    if (lines) {
        wdata = cpl_vector_get_data(lines);
        numLines = cpl_vector_get_size(lines);
    }
    else {
        cpl_msg_warning(func, "A catalog of sky lines wavelengths was not "
                        "given: using internal list of reference sky lines");
        if (highres) {
           wdata = default_lines_hi;
           numLines = sizeof(default_lines_hi) / sizeof(double);
        }
        else {
           wdata = default_lines_lo;
           numLines = sizeof(default_lines_lo) / sizeof(double);
        }
    }

    npix = 2 * radius + 1;
    profile = cpl_calloc(npix, sizeof(float));

    rms = 0.0;
    allPeaks = 0;

    for (i = 0; i < numLines; i++) {

        /*
         *  Expected peak and closest pixel to specified wavelength.
         */

        lambda = wdata[i];
        expectPeak = (lambda - wavestart) / dispersion;
        cpix = floor(expectPeak + 0.5);

        /*
         *  Search interval for peak. Abort if too close to image border.
         */

        sp = cpix - radius;
        ep = cpix + radius;

        if (sp < 0 || ep > xlen)
            continue;

        average = 0.0;
        npeaks = 0;
        oneRms = 0.0;

        for (j = 0; j < ylen; j++) {    /*  For each row of each slit  */
            nzero = 0;
            for (k = 0; k < npix; k++) {
                profile[k] = sdata[sp + k + j * xlen];
                if (fabs(profile[k]) < 0.0001)
                    nzero++; /* Count number of 0 pixels (spectrum truncated) */
            }
            if (nzero > 0)
                continue;

            if (peakPosition(profile, npix, &peak, 1) == 0) {
                offset = (sp + peak) - expectPeak;
                average += offset;
                rms += fabs(offset);
                oneRms += fabs(offset);
                npeaks++;
                allPeaks++;
            }
        }

        if (npeaks)
            cpl_msg_info(func, "RMS for %.2f: %.3f pixel (%d points)",
                         lambda, oneRms / npeaks * 1.25, npeaks);
        else
            cpl_msg_info(func, "RMS for %.2f: line not available", lambda);
    }

    cpl_free(profile);

    if (allPeaks < 10)
        return 0.0;

    rms /= allPeaks;
    rms *= 1.25;       /* Factor to convert average deviation to sigma */

    return rms;

}


/**
 * @brief
 *   Create a pixel map from an IDS coefficients table
 *
 * @param idscoeff    Table with IDS polynomials
 * @param reference   Reference wavelength
 * @param blue        Start lambda for spatial remapping
 * @param red         End lambda for spatial remapping
 * @param dispersion  Mean spectral dispersion
 * @param trend       Trend to remove from IDS polynomials
 *
 * @return Pixel map
 * 
 * The output pixel map will have the Y size equal to the
 * number of rows in the @em idscoeff table. The X size is
 * computed as (@em red - @em blue) / @em dispersion. With @em trend
 * the number of polynomial coefficients to ignore in the computation
 * can be specified.
 */

cpl_image *mos_map_pixel(cpl_table *idscoeff, double reference,
                         double blue, double red, double dispersion, int trend)
{
    const char *func = "mos_map_pixel";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_polynomial *ids;
    cpl_image      *map;
    float          *mdata;
    double          lambda;
    double          c;
    int             order;
    int             xsize, ysize;
    int             missing;
    int             i, j;
    cpl_size        k;


    if (idscoeff == NULL) {
        cpl_msg_error(func, "An IDS coeff table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    xsize = (red - blue) / dispersion;
    ysize = cpl_table_get_nrow(idscoeff);
    map = cpl_image_new(xsize, ysize, CPL_TYPE_FLOAT);
    mdata = cpl_image_get_data(map);

    order = 0;
    while (order < 6 && cpl_table_has_column(idscoeff, clab[order]))
        ++order;
    --order;

    for (i = 0; i < ysize; i++, mdata += xsize) {

        missing = 0;
        ids = cpl_polynomial_new(1);
        for (k = trend; k <= order; k++) {
            c = cpl_table_get_double(idscoeff, clab[k], i, &missing);
            if (missing) {
                cpl_polynomial_delete(ids);
                break;
            }
            cpl_polynomial_set_coeff(ids, &k, c);
        }
        if (missing)
            continue;

        for (j = 0; j < xsize; j++) {
            lambda = blue + j*dispersion;
            mdata[j] = cpl_polynomial_eval_1d(ids, lambda-reference, NULL);
        }

        cpl_polynomial_delete(ids);
    }

    return map;

}


/**
 * @brief
 *   Create a wavelengths map from an IDS coefficients table
 *
 * @param idscoeff    Table with IDS polynomials
 * @param xsize       X size of used CCD
 * @param reference   Reference wavelength
 * @param blue        Start lambda for spatial remapping
 * @param red         End lambda for spatial remapping
 *
 * @return Wavelengths map
 * 
 * The output wavelength map will have the Y size equal to the
 * number of rows in the @em idscoeff table. The input @em xsize
 * will typically be the X size of the used detector. In general
 * this should be identical to the X size of the arc lamp image
 * from where the @em idscoeff table was derived. The wavelengths
 * map is created by applying the inverses of the IDS polynomials
 * to all the pixel positions of each output image row.
 */

cpl_image *mos_map_idscoeff(cpl_table *idscoeff, int xsize, double reference,
                            double blue, double red)
{
    const char *func = "mos_map_idscoeff";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_polynomial *ids;
    cpl_image      *map;
    float          *mdata;
    double          lambda;
    double          c;
    int             order;
    int             ysize;
    int             missing;
    int             i, j;
    cpl_size        k;


    if (idscoeff == NULL) {
        cpl_msg_error(func, "An IDS coeff table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (xsize < 1) {
        cpl_msg_error(func, "Invalid image size");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (xsize < 20 || xsize > 5000) {
        cpl_msg_warning(func, "Do you really have a detector %d pixels long?",
                        xsize);
    }

    ysize = cpl_table_get_nrow(idscoeff);
    map = cpl_image_new(xsize, ysize, CPL_TYPE_FLOAT);
    mdata = cpl_image_get_data(map);

    order = 0;
    while (order < 6 && cpl_table_has_column(idscoeff, clab[order]))
        ++order;
    --order;

    for (i = 0; i < ysize; i++, mdata += xsize) {

        missing = 0;
        ids = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            c = cpl_table_get_double(idscoeff, clab[k], i, &missing);
            if (missing) {
                cpl_polynomial_delete(ids);
                break;
            }
            cpl_polynomial_set_coeff(ids, &k, c);
        }
        if (missing)
            continue;

        for (j = 0; j < xsize; j++) {
            lambda = mos_eval_dds(ids, blue, red, reference, j);

            if (lambda >= blue && lambda <= red) {
                mdata[j] = lambda;
            }
        }

        cpl_polynomial_delete(ids);
    }

    return map;

}


/**
 * @brief
 *   Remapping of spatially rectified wavelengths to original CCD pixels
 *
 * @param spatial     CCD image of spatial coordinates
 * @param calibration Spatially rectified image of wavelengths
 * @param slits       Slits positions on the CCD
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param reference   Reference wavelength
 * @param blue        Start lambda for spatial remapping
 * @param red         End lambda for spatial remapping
 * @param dispersion  Mean spectral dispersion
 *
 * @return Wavelengths mapping on CCD.
 *
 * The input @em spatial image is the one produced by the function
 * @c mos_spatial_calibration() (argument: @em calibration). It is 
 * expected to be oriented with horizontal dispersion direction and 
 * red wavelengths on the right side, and it should have the same x-length 
 * of the input @em calibration image. The @em calibration image is the
 * one produced by the function @c mos_wavelength_calibration_final() 
 * (argument: @em calibration). The @em slits table should be the 
 * same processed by the function @c mos_spatial_calibration(). The 
 * @em polytraces table is the output of the function @c mos_poly_trace().
 * The other arguments should be in principle the same specified for
 * the function @c mos_spatial_calibration().
 *
 * For each slit, each (x,y) coordinate on the @em spatial image is 
 * characterised by a spatial coordinate p. The wavelength at that
 * spatial coordinate is derived by linear interpolation of the two
 * wavelengths with the same coordinate x that are closer to p on
 * the @em calibration image.
 */

cpl_image *mos_map_wavelengths(cpl_image *spatial, cpl_image *calibration,
                               cpl_table *slits, cpl_table *polytraces, 
                               double reference, double blue, double red, 
                               double dispersion)
{
    const char *func = "mos_map_wavelengths";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;
    cpl_image      *remapped;
    float          *data;
    float          *wdata;
    float          *sdata;
    float          *xdata;
    double          vtop, vbot, value;
    double          top, bot;
    double          coeff;
    double          ytop, ybot;
    double          ypos;
    double          fvalue;
    int             ivalue;
    int             yint, ysize, yprev;
    int             nslits;
    int             npseudo;
    int            *slit_id;
    int            *position;
    int            *length;
    int             nx, ny;
    int             pixel_above, pixel_below, refpixel, start_pixel, end_pixel;
    int             missing_top, missing_bot;
    int             null;
    int             order;
    int             i, j;
    cpl_size        k;


    if (spatial == NULL || calibration == NULL || 
        slits == NULL || polytraces == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (red - blue < dispersion) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    nx = cpl_image_get_size_x(spatial);
    ny = cpl_image_get_size_y(spatial);
    ysize = cpl_image_get_size_y(calibration);
    remapped = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
    data  = cpl_image_get_data(remapped);
    sdata = cpl_image_get_data(spatial);
    wdata = cpl_image_get_data(calibration);

    nslits   = cpl_table_get_nrow(slits);
    slit_id  = cpl_table_get_data_int(slits, "slit_id");
    order    = cpl_table_get_ncol(polytraces) - 2;
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");

    /*
     * The spatial resampling is performed for a certain number of 
     * pixels above and below the position of the reference wavelength:
     */

    pixel_above = STRETCH_FACTOR * (red - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - blue) / dispersion;

    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        /*
         * Note that the x coordinate of the reference pixels on the CCD
         * is taken arbitrarily at the top end of each slit. This wouldn't
         * be entirely correct in case of curved slits, or in presence of
         * heavy distortions: in such cases the spatial resampling is
         * really performed across a wide range of wavelengths. But
         * the lag between top and bottom spectral curvature models 
         * would introduce even in such cases negligible effects on
         * the spectral spatial resampling.
         */

        refpixel = cpl_table_get_double(slits, "xtop", i, NULL);

        start_pixel = refpixel - pixel_below;
        if (start_pixel < 0)
            start_pixel = 0;

        end_pixel = refpixel + pixel_above;
        if (end_pixel > nx)
            end_pixel = nx;

        /*
         * Recover from the table of spectral curvature coefficients
         * the curvature polynomials.
         */

        missing_top = 0;
        polytop = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i, &null);
            if (null) {
                cpl_polynomial_delete(polytop);
                missing_top = 1;
                break;
            }
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        missing_bot = 0;
        polybot = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i+1, &null);
            if (null) {
                cpl_polynomial_delete(polybot);
                missing_bot = 1;
                break;
            }
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        if (missing_top && missing_bot) {
            cpl_msg_debug(func, "Slit %d was not traced: no extraction!", 
                          slit_id[i]);
            continue;
        }

        /*
         * In case just one of the two edges was not traced, the other
         * edge curvature model is duplicated and shifted to the other
         * end of the slit: better than nothing!
         */

        if (missing_top) {
            cpl_msg_debug(func, "Upper edge of slit %d was not traced: "
                          "the spectral curvature of the lower edge "
                          "is used instead.", slit_id[i]);
            polytop = cpl_polynomial_duplicate(polybot);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polybot, &k);
            coeff += ytop - ybot;
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        if (missing_bot) {
            cpl_msg_debug(func, "Lower edge of slit %d was not traced: "
                          "the spectral curvature of the upper edge "
                          "is used instead.", slit_id[i]);
            polybot = cpl_polynomial_duplicate(polytop);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polytop, &k);
            coeff -= ytop - ybot;
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        /*
         * Point to current slit on wavelength calibration image.
         * Note that the npseudo value related to this slit is equal 
         * to the number of spatial pseudo-pixels decreased by 1 
         * (compare with function mos_spatial_calibration()).
         */

        xdata = wdata + nx*position[i];
        npseudo = length[i] - 1;

        /*
         * Write interpolated wavelengths to CCD image
         */

        for (j = start_pixel; j < end_pixel; j++) {
            top = cpl_polynomial_eval_1d(polytop, j, NULL);
            bot = cpl_polynomial_eval_1d(polybot, j, NULL);
            for (k = 0; k <= npseudo; k++) {
                ypos = top - k*(top-bot)/npseudo;
                yint = ypos;

                /* 
                 * The line:
                 *     value = sdata[j + nx*yint];
                 * should be equivalent to:
                 *     value = npseudo*(top-yint)/(top-bot);
                 */

                if (yint < 0 || yint >= ny-1) {
                    yprev = yint;
                    continue;
                }

                value = sdata[j + nx*yint];   /* Spatial coordinate on CCD */
                ivalue = value;               /* Nearest spatial pixels:   */
                fvalue = value - ivalue;      /* ivalue and ivalue+1       */
                if (ivalue < npseudo && ivalue >= 0) {
                    vtop = xdata[j + nx*(npseudo-ivalue)];
                    vbot = xdata[j + nx*(npseudo-ivalue-1)];
                    if (vtop < 1.0) {  /* Impossible wavelength */
                        if (vbot < 1.0) {
                            value = 0.0;
                        }
                        else {
                            value = vbot;
                        }
                    }
                    else if (vbot < 1.0) {
                        if (k)
                            value = vtop;
                        else
                            value = 0.0;
                    }
                    else if (fabs(vbot-vtop) > 10*dispersion) {
                        value = 0.0;
                    }
                    else {
                        value = vtop*(1-fvalue) + vbot*fvalue;
                    }
                    data[j + nx*yint] = value;

                    if (k) {

                        /*
                         * This is added to recover lost pixels on
                         * the CCD image (pixels are lost because
                         * the CCD pixels are less than npseudo+1).
                         */

                        if (yprev - yint > 1) {
                            value = sdata[j + nx*(yint+1)];
                            ivalue = value;
                            fvalue = value - ivalue;
                            if (ivalue < npseudo && ivalue >= 0) {
                                vtop = xdata[j + nx*(npseudo-ivalue)];
                                vbot = xdata[j + nx*(npseudo-ivalue-1)];
                                if (vtop < 1.0) {
                                    if (vbot < 1.0) {
                                        value = data[j + nx*(yint+1)];
                                    }
                                    else {
                                        value = vbot;
                                    }
                                }
                                else if (vbot < 1.0) {
                                    value = vtop;
                                }
                                else if (fabs(vbot-vtop) > 2*dispersion) {
                                    value = vtop;
                                }
                                else {
                                    value = vtop*(1-fvalue) + vbot*fvalue;
                                }
                                data[j + nx*(yint+1)] = value;
                            }
                        }
                    }
                }
                yprev = yint;
            }
        }
        cpl_polynomial_delete(polytop);
        cpl_polynomial_delete(polybot);
    }

    return remapped;
}

/**
 * @brief
 *   Remapping of slit spectra into a grid of lambda-space coordinates
 *   
 * @param spectra     CCD image containing the observed slit spectra
 * @param wavecalib   CCD image of wavelengths
 * @param spatial     CCD image of spatial coordinates
 * @param slits       Slits positions on the CCD
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param reference   Reference wavelength
 * @param blue        Start lambda for remapping
 * @param red         End lambda for remapping
 * @param dispersion  Mean spectral dispersion
 * @param flux        flux = 0 means no flux conservation correction applied
 * 
 * @return Extracted slit spectra
 *
 * DO NOT USE THIS FUNCTION, IT GIVES BAD RESULTS (and it should be
 * eventually removed).
 * 
 * The input @em wavecalib image is the one returned by the function
 * @c mos_map_wavelengths() or, if this is not available, the approximate
 * one returned by the function @em mos_wavelength_calibration_raw().
 * The input @em spatial image is the one produced by the function
 * @c mos_spatial_calibration() (argument: @em calibration). Both images
 * are expected to be oriented with horizontal dispersion direction and
 * red wavelengths on the right side, and they should have the same sizes
 * of the input @em spectra image. The @em slits table should be the 
 * same processed by the function @c mos_spatial_calibration(). The
 * @em polytraces table is the output of the function @c mos_poly_trace().
 * If @em flux is different from 0, the factors to be applied for flux 
 * conservation are calculated and applied to the remapped image.
 * The other arguments should be in principle the same specified for
 * the function @c mos_spatial_calibration(), even if it would be
 * conceivable (and perhaps even reasonable) to specify shorter
 * spectral ranges, and/or a smaller value of the dispersion (A/pixel)
 * for supersampling the signal. 
 * 
 * The target grid in lambda and pseudo spatial coordinates is defined
 * on the basis of the indicated spectral range and dispersion, while
 * the number of spatial pixels assigned to each slits are compatible
 * with the positions listed in the @em slits table (column "position")
 * following the convention used in the rest of the spectral reduction,
 * in particular the functions @em mos_spatial_calibration() and 
 * @c mos_map_wavelengths().
 *
 * Indicating with (x,y) a pixel position on the CCD, and with (L,S) a
 * pixel position on one extracted slit, the mapping from (x,y) to (L,S)
 * is done in the following way (for each slit): for each (x,y) read 
 * the corresponding value of lambda and space from the input images 
 * @em wavecalib and @em spatial; find the pixel (L,S) on the output 
 * (remapped) slit that has the highest lambda and space that are less 
 * than the lambda and space obtained above; read the values of the 
 * following CCD pixels from the input spectrum: (x,y), (x-1,y), (x+N,y+1), 
 * (x+N-1,y+1), where N is an offset pointing to the pixel on the row 
 * y+1 having the closest wavelength to the wavelength of pixel (x,y). 
 * This is a way to avoid possible discontinuities in the wavelength 
 * calibration along the cross-dispersion direction (due to slit 
 * irregularities, also commonly indicated with the expression "slit 
 * geometry"). Typically, it will always be N = 0.
 * 
 * The pixel value to assign to the pixel (L,S) is computed by linear 
 * interpolation to its corresponding (x',y') position of the 4 pixel 
 * values (horizontal interpolation followed by vertical interpolation 
 * of the interpolated values). The corresponding flux-conservation 
 * factor, if requested, is computed as (dL/dx)*(dS/dy). The interpolated 
 * value is multiplied by this factor before being written to the pixel 
 * (L,S).
 *
 * DO NOT USE THIS FUNCTION, IT GIVES BAD RESULTS (and it should be
 * eventually removed).
 */

cpl_image *mos_map_spectrum(cpl_image *spectra, cpl_image *wavecalib, 
                            cpl_image *spatial, cpl_table *slits,
                            cpl_table *polytraces, double reference,
                            double blue, double red, double dispersion,
                            int flux)
{
    const char *func = "mos_map_spectrum";
    
    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;
    cpl_image      *remapped;
    cpl_image     **exslit;
    float          *data;
    float          *wdata;
    float          *sdata;
    float          *xdata;
    double          lambda00, lambda01, lambda10, lambda11, lambda;
    double          space00, space01, space10, space11, space;
    double          value00, value01, value10, value11, value0, value1, value;
    double          dL, dS;
    double          top, bot;
    double          coeff;
    double          ytop, ybot;
    double          xfrac, yfrac;
    int             yint, ysize;
    int             itop, ibot;
    int             shift;
    int             L, S;
    int             nslits;
    int             npseudo;
    int            *slit_id;
    int            *position;
    int            *length;
    int             nx, ny;
    int             x, y;
    int             nlambda;
    int             pixel_above, pixel_below, refpixel, start_pixel, end_pixel;
    int             missing_top, missing_bot; 
    int             null;
    int             order;
    int             i; 
    cpl_size        k;
    

    flux += flux;

    if (spectra == NULL || spatial == NULL || wavecalib == NULL ||
        slits == NULL || polytraces == NULL) { 
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (red - blue < dispersion) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }
    
    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);

    if (nx != cpl_image_get_size_x(spatial) ||
        ny != cpl_image_get_size_y(spatial) ||
        nx != cpl_image_get_size_x(wavecalib) ||
        ny != cpl_image_get_size_y(wavecalib)) {
        cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
        return NULL;
    }

    nlambda     = STRETCH_FACTOR * (red - blue) / dispersion;
    pixel_above = STRETCH_FACTOR * (red - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - blue) / dispersion;

    data  = cpl_image_get_data(spectra);
    sdata = cpl_image_get_data(spatial);
    wdata = cpl_image_get_data(wavecalib);
    
    nslits   = cpl_table_get_nrow(slits);
    slit_id  = cpl_table_get_data_int(slits, "slit_id");
    order    = cpl_table_get_ncol(polytraces) - 2;
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");
    
    exslit = cpl_calloc(nslits, sizeof(cpl_image *));

    for (i = 0; i < nslits; i++) {

         if (length == 0)
             continue;

        /*
         * Note that the x coordinate of the reference pixels on the CCD
         * is taken arbitrarily at the top end of each slit. This wouldn't
         * be entirely correct in case of curved slits, or in presence of
         * heavy distortions: in such cases the spatial resampling is
         * really performed across a wide range of wavelengths. But
         * the lag between top and bottom spectral curvature models
         * would introduce even in such cases negligible effects on
         * the spectral spatial resampling.
         */

        refpixel = cpl_table_get_double(slits, "xtop", i, NULL);

        start_pixel = refpixel - pixel_below;
        if (start_pixel < 1)
            start_pixel = 1;

        end_pixel = refpixel + pixel_above;
        if (end_pixel > nx)
            end_pixel = nx;

        /*
         * Recover from the table of spectral curvature coefficients
         * the curvature polynomials.
         */

        missing_top = 0;
        polytop = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i, &null);
            if (null) {
                cpl_polynomial_delete(polytop);
                missing_top = 1;
                break;
            }
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        missing_bot = 0;
        polybot = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i+1, &null);
            if (null) {
                cpl_polynomial_delete(polybot);
                missing_bot = 1;
                break;
            }
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        if (missing_top && missing_bot) {
            cpl_msg_debug(func, "Slit %d was not traced: no extraction!",
                          slit_id[i]);
            continue;
        }

        /*
         * In case just one of the two edges was not traced, the other
         * edge curvature model is duplicated and shifted to the other
         * end of the slit: better than nothing!
         */

        if (missing_top) {
            cpl_msg_debug(func, "Upper edge of slit %d was not traced: "
                          "the spectral curvature of the lower edge "
                          "is used instead.", slit_id[i]);
            polytop = cpl_polynomial_duplicate(polybot);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polybot, &k);
            coeff += ytop - ybot;
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        if (missing_bot) {
            cpl_msg_debug(func, "Lower edge of slit %d was not traced: "
                          "the spectral curvature of the upper edge "
                          "is used instead.", slit_id[i]);
            polybot = cpl_polynomial_duplicate(polytop);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polytop, &k);
            coeff -= ytop - ybot;
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        /*
         * Allocate image for current extracted slit
         */

        top = cpl_polynomial_eval_1d(polytop, refpixel, NULL);
        bot = cpl_polynomial_eval_1d(polybot, refpixel, NULL);
        npseudo = ceil(top-bot) + 1;

        if (npseudo < 1) {
            cpl_polynomial_delete(polytop);
            cpl_polynomial_delete(polybot);
            cpl_msg_debug(func, "Slit %d was badly traced: no extraction!",
                          slit_id[i]);
            continue;
        }

        exslit[i] = cpl_image_new(nlambda, npseudo+1, CPL_TYPE_FLOAT);
        xdata = cpl_image_get_data(exslit[i]);

        /*
         * Write interpolated spectral values to remapped slit spectrum.
         */

        for (x = start_pixel; x < end_pixel; x++) {
            top = cpl_polynomial_eval_1d(polytop, x, NULL);
            bot = cpl_polynomial_eval_1d(polybot, x, NULL);
            itop = top + 1;
            ibot = bot;
            if (itop < 0)
                itop = 0;
            if (itop > ny - 1)
                itop = ny - 1;
            if (ibot < 0)
                ibot = 0;
            if (ibot > ny - 1)
                ibot = ny - 1;
            for (y = ibot; y < itop; y++) {
                 lambda11 = wdata[x + y*nx];
                 if (lambda11 < 1.0)        /* Impossible wavelength */
                     continue;
                 space11 = sdata[x + y*nx];
                 if (space11 < 0.0)         /* Impossible spatial coordinate */
                     continue;
                 lambda01 = wdata[x - 1 + y*nx];
                 if (lambda01 < 1.0)        /* Impossible wavelength */
                     continue;
                 space01 = sdata[x - 1 + y*nx];
                 if (space01 < 0.0)         /* Impossible spatial coordinate */
                     continue;

                 shift = 0;

/****+
                 if (wdata[x + (y+1)*nx] > 1.0) {
                     if (wdata[x + (y+1)*nx] - lambda11 > 0) {
                         shift = -1;
                         while (wdata[x + shift + (y+1)*nx] - lambda11 > 0)
                             shift--;
                         if (lambda11 - wdata[x + shift + (y+1)*nx] > 
                             wdata[x + shift + 1 + (y+1)*nx] - lambda11) {
                             shift++;
                         }
                     }
                     else {
                         shift = 1;
                         while (wdata[x + shift + (y+1)*nx] - lambda11 < 0)
                             shift++;
                         if (wdata[x + shift + (y+1)*nx] - lambda11 >
                             lambda11 - wdata[x + shift + 1 + (y+1)*nx]) {
                             shift--;
                         }
                     }
                 }
****/

/****
printf("y = %d, shift = %d\n", y, shift);
****/

                 lambda10 = wdata[x + shift + (y+1)*nx];
                 if (lambda10 < 1.0)        /* Impossible wavelength */
                     continue;
                 space10 = sdata[x + shift + (y+1)*nx];
                 if (space10 < 0.0)         /* Impossible spatial coordinate */
                     continue;
                 lambda00 = wdata[x - 1 + shift + (y+1)*nx];
                 if (lambda00 < 1.0)        /* Impossible wavelength */
                     continue;
                 space00 = sdata[x - 1 + shift + (y+1)*nx];
                 if (space00 < 0.0)         /* Impossible spatial coordinate */
                     continue;
                 
                 /*
                  * Find the variation in lambda and space in this
                  * position for each CCD pixel (both quantities are 
                  * expected to be positive).
                  */

                 dL = lambda11 - lambda01;
                 dS = space11 - space10;

                 /*
                  * Find the position (L,S) of the output pixel 
                  * (by integer truncation).
                  */

                 L = (lambda11 - blue)/dispersion + 0.5;
                 S = space11 + 0.5;                   /* Counted from top! */

                 if (L < 0 || L >= nlambda)
                     continue;
                 if (S < 0 || S > npseudo)
                     continue;

                 /*
                  * Find the coordinate of pixel (L,S)
                  */

                 lambda = blue + L*dispersion;
                 space  = S;

                 /*
                  * Find the interpolation point on the CCD: it is
                  * defined as the (positive) distance from current
                  * CCD pixel (x,y) of the interpolation point (x',y'),
                  * measured in CCD pixels. The interpolation point
                  * is located between the four CCD pixels selected
                  * above.
                  */

                 xfrac = (lambda11-lambda)/dL;
                 yfrac = (space11-space)/dS;

/*
if (xfrac < 0.0 || xfrac > 1.0 || yfrac < 0.0 || yfrac > 1.0)
printf("xyfrac = %f, %f\n", xfrac, yfrac);
*/

                 /*
                  * Get the four values to interpolate
                  */

                 value11 = data[x + y*nx];
                 value01 = data[x - 1 + y*nx];
                 value10 = data[x + shift + (y+1)*nx];
                 value00 = data[x + shift - 1 + (y+1)*nx];

                 /*
                  * Interpolation
                  */

                 value1 = (1-xfrac)*value11 + xfrac*value01;
                 value0 = (1-xfrac)*value10 + xfrac*value00;
                 value  = (1-yfrac)*value1  + yfrac*value0;

                 /*
                  * Write this value to the appropriate (L,S) coordinate
                  * on output slit
                  */

                 xdata[L + nlambda*(npseudo-S)] = value;
                 
            }
        }
        cpl_polynomial_delete(polytop);
        cpl_polynomial_delete(polybot);
    }

    /*
     * Now all the slits images are copied to a single image
     */

    ysize = 0;
    for (i = 0; i < nslits; i++)
        if (exslit[i])
            ysize += cpl_image_get_size_y(exslit[i]);

    remapped = cpl_image_new(nlambda, ysize, CPL_TYPE_FLOAT);

    yint = -1;
    for (i = 0; i < nslits; i++) {
        if (exslit[i]) {
            yint += cpl_image_get_size_y(exslit[i]);
            cpl_image_copy(remapped, exslit[i], 1, ysize - yint);
            cpl_image_delete(exslit[i]);
            cpl_table_set_int(slits, "position", i, ysize - yint - 1);
        }
    }

    cpl_free(exslit);

    return remapped;

}


/** 
 * @brief
 *   Create a CCD median sky map
 *  
 * @param spectra     CCD image of spectra
 * @param wavemap     CCD image of wavelengths
 * @param dispersion  Mean spectral dispersion
 * @param factor      Supersampling factor
 * @param minpoints   Minimum points required per supersampled spectrum bin
 * @param skymap      Returned CCD image of sky
 *  
 * @return Table with median sky spectrum
 *  
 * The spectra contained in the input image, @em spectra, taken all
 * together are assumed to contain at least 50% of their pixels on 
 * the sky. Moreover, all slits are assumed to have the same width.
 * The CCD image of wavelengths, @em wavemap, is the one produced 
 * by the function @c mos_wavelength_calibration_final(), possibly 
 * modified by the function @c mos_wavelength_align() used to align 
 * the wavelength solution to available sky lines. The sky spectrum 
 * is computed as the median level of all the pixel values of all 
 * the CCD spectra in the wavelength interval corresponding to that 
 * bin, on bins containing at least @em minpoints; otherwise it is 
 * computed by linear interpolation from the nearest valid bins.
 * To each bin is assigned in the first case the median of the 
 * contributing wavelengths, and in the second case its central
 * wavelength. The @em skymap image should be preallocated: each 
 * one of its pixels, corresponding to a wavelength read from @em wavemap, 
 * will be assigned a value obtained by linear interpolation of the two 
 * pixels of the supersampled spectrum that are closest to its wavelength.
 */

cpl_table *mos_sky_map_super(cpl_image *spectra, cpl_image *wavemap,
                             double dispersion, double factor, int minpoints,
                             cpl_image *skymap)
{
    const char *func = "mos_sky_map_super";

    cpl_vector **vector;
    cpl_vector **wvector;
    double       firstLambda, lastLambda;
    double       lambda, lambda1, lambda2;
    double       value, value1, value2;
    double       frac;
    float        min, max;
    int         *count;
    int          nbin, bin;
    int          nx, ny, npix;
    int          first_valid, valid_bins;
    int          i, j;

    cpl_table   *sky;
    double      *sky_spectrum;
    double      *sky_wave;
    float       *data;
    float       *sdata;
    float       *kdata;


    if (spectra == NULL || wavemap == NULL || skymap == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }
    
    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        cpl_msg_error(func, "Negative dispersion: %s", cpl_error_get_message());
        return NULL;
    }

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);
    npix = nx * ny;

    if (nx != cpl_image_get_size_x(wavemap) ||
        ny != cpl_image_get_size_y(wavemap) ||
        nx != cpl_image_get_size_x(skymap) ||
        ny != cpl_image_get_size_y(skymap)) {
        cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
        cpl_msg_error(func, "Image sizes: %s", cpl_error_get_message());
        return NULL;
    }

    if (factor < 1.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        cpl_msg_error(func, "Undersampling (%f): %s", factor, 
                      cpl_error_get_message());
        return NULL;
    }

    if (minpoints < 0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        cpl_msg_error(func, "Negative threshold: %s", cpl_error_get_message());
        return NULL;
    }

    dispersion /= factor;


    /*
     * Find bluest and reddest wavelengths in the whole image
     */

    data = cpl_image_get_data(wavemap);

    j = -1;
    for (i = 0; i < npix; i++) {
        if (data[i] > 1.0) {
            min = max = data[i];
            j = i+1;
            break;
        }
    }
    
    if(j == -1) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        cpl_msg_warning(func, "Wavelength map has no valid values: %s", 
                      cpl_error_get_message());
        return NULL;
    }

    for (i = j; i < npix; i++) {
        if (data[i] < 1.0)      /* Impossible wavelength */
            continue;
        if (min > data[i])
            min = data[i];
        if (max < data[i])
            max = data[i];
    }

    firstLambda = min;
    lastLambda = max;


    /*
     * Determine length of median spectrum
     */

    nbin = (lastLambda - firstLambda) / dispersion;

    /*
     * Count how many values will be found for each spectral bin.
     * The ith bin has a wavelength range from firstLambda + i*dispersion
     * (inclusive) to firstLambda + (i+1)*dispersion (exclusive), and
     * it is assigned to its central wavelength.
     */

    count = cpl_calloc(nbin, sizeof(int));

    data = cpl_image_get_data(wavemap);

    for (i = 0; i < npix; i++) {
        if (data[i] < 1.0)
            continue;
        bin = (data[i] - firstLambda) / dispersion;
        if (bin < nbin)                               /* Safer */
            count[bin]++;
    }

    valid_bins = 0;
    for (i = 0; i < nbin; i++)
        if (count[i] >= minpoints)
            valid_bins++;

    if (valid_bins < nbin/3) {
        cpl_msg_warning(func, "Cannot determine a good global sky "
                        "spectrum from input data");
        return NULL;
    }


    /*
     * Allocate an array of vectors with the appropriate size, to
     * contain a list of all the spectral pixels values. At the same
     * time, reset the array of counters (because we will have to
     * count again...).
     */

    vector = cpl_calloc(nbin, sizeof(cpl_vector *));
    wvector = cpl_calloc(nbin, sizeof(cpl_vector *));
    for (i = 0; i < nbin; i++) {
        if (count[i] >= minpoints) {
            vector[i] = cpl_vector_new(count[i]);
            wvector[i] = cpl_vector_new(count[i]);
        }
        count[i] = 0;
    }


    /*
     * Read the wavemap and the spectral images, and add the data values
     * to the appropriate wavelength bins
     */

    data  = cpl_image_get_data(wavemap);
    sdata = cpl_image_get_data(spectra);

    for (i = 0; i < npix; i++) {
        if (data[i] < 1.0)
            continue;
        bin = (data[i] - firstLambda) / dispersion;
        if (bin < nbin) {                             /* Safer */
            if (vector[bin]) {
                cpl_vector_set(vector[bin], count[bin], sdata[i]);
                cpl_vector_set(wvector[bin], count[bin], data[i]);
            }
            count[bin]++;
        }
    }


    /*
     * Compute the median flux for each wavelength bin, and destroy
     * at the same time the used vectors
     */

    sky_spectrum = cpl_calloc(nbin, sizeof(double));
    sky_wave = cpl_calloc(nbin, sizeof(double));
    for (i = 0; i < nbin; i++) {
        if (vector[i]) {
            sky_spectrum[i] = cpl_vector_get_median_const(vector[i]);
            sky_wave[i] = cpl_vector_get_median_const(wvector[i]);
            cpl_vector_delete(vector[i]);
            cpl_vector_delete(wvector[i]);
        }
    }

    cpl_free(vector);
    cpl_free(wvector);


    /*
     * Here possible gaps in the final spectrum are filled by interpolation
     */

    for (i = 0; i < nbin; i++) {
        if (count[i] >= minpoints) {
            first_valid = i;
            break;
        }
    }
    
    for (i = first_valid; i < nbin; i++) {
        if (count[i] < minpoints) {
            sky_wave[i] = firstLambda + (i+0.5)*dispersion;
            for (j = i+1; j < nbin; j++) {
                if (count[j] >= minpoints) {
                    if (sky_wave[j] - sky_wave[i-1] < 0.1) {
                        sky_spectrum[i] = (sky_spectrum[j] + sky_spectrum[i-1])
                                        / 2;
                    }
                    else {
                        frac = (sky_wave[i] - sky_wave[i-1]) 
                             / (sky_wave[j] - sky_wave[i-1]);
                        sky_spectrum[i] = frac * sky_spectrum[j]
                                        + (1 - frac) * sky_spectrum[i-1];
                    }
                }
            }
        }
    }


    /*
     * Create the output table
     */

    sky = cpl_table_new(nbin);
    cpl_table_wrap_double(sky, sky_wave, "wavelength");
    cpl_table_wrap_double(sky, sky_spectrum, "sky");
    cpl_table_wrap_int(sky, count, "npoints");


    /*
     * Fill the sky map
     */

    data  = cpl_image_get_data(wavemap);
    sdata = cpl_image_get_data(spectra);
    kdata = cpl_image_get_data(skymap);

    for (i = 0; i < npix; i++) {

        /*
         * Currently based on linear interpolation
         */

        lambda = data[i];
        if (lambda < 1.0)
            continue;
        bin = (lambda - firstLambda) / dispersion;
        if (bin >= nbin)                               /* Safer */
            continue;
        lambda1 = sky_wave[bin];
        value1 = sky_spectrum[bin];
        if (lambda1 < lambda) {
            bin++;
            if (bin < nbin) {
                lambda2 = sky_wave[bin];
                value2  = sky_spectrum[bin];
                if (lambda2 - lambda1 < 0.1) {
                    value = (value1 + value2) / 2;
                }
                else {
                    frac = (lambda - lambda1) / (lambda2 - lambda1);
                    value = frac * value2 + (1 - frac) * value1;
                }
            }
            else {
                value = value1;
            }
        }
        else {
            if (bin > 0) {
                bin--;
                lambda2 = lambda1;
                value2  = value1;
                lambda1 = sky_wave[bin];
                value1  = sky_spectrum[bin];
                if (lambda2 - lambda1 < 0.1) {
                    value = (value1 + value2) / 2;
                }
                else {
                    frac = (lambda - lambda1) / (lambda2 - lambda1);
                    value = frac * value2 + (1 - frac) * value1;
                }
            }
            else {
                value = value1;
            }
        }
        kdata[i] = value;
    }

    if (first_valid)
        cpl_table_erase_window(sky, 0, first_valid);

    return sky;

}


/** 
 * @brief
 *   Create a CCD median sky map
 *  
 * @param spectra     CCD image of spectra
 * @param wavemap     CCD image of wavelengths
 * @param dispersion  Mean spectral dispersion
 * @param skymap      Returned CCD image of sky
 *  
 * @return Table with median sky spectrum
 *  
 * The spectra contained in the input image, @em spectra, taken all
 * together are assumed to contain at least 50% of their pixels on 
 * the sky. Moreover, all slits are assumed to have the same width.
 * The CCD image of wavelengths, @em wavemap, is the one produced 
 * by the function @c mos_wavelength_calibration_final(), possibly 
 * modified by the function @c mos_wavelength_align() used to align 
 * the wavelength solution to available sky lines. The computed sky 
 * spectrum will consist of pixels of size @em dispersion: typically
 * the passed dispersion should be slightly smaller than the dispersion 
 * of the original spectra. How much smaller will affect the level of
 * supersampling of the computed sky spectrum: a too small value will
 * increase the noise of the output sky spectrum, leading to no advantage
 * in the determination of the sky map to be subtracted from the scientific
 * data. Each pixel of the resulting sky spectrum is computed as the 
 * median of all the pixel values of all the CCD spectra in the wavelength 
 * interval corresponding to that pixel. The @em skymap image should be 
 * preallocated: each one of its pixels, corresponding to a wavelength 
 * read from @em wavemap, will be assigned a value obtained by linear 
 * interpolation of the two pixels of the supersampled spectrum that are 
 * closest to its wavelength.
 */

cpl_table *mos_sky_map(cpl_image *spectra, cpl_image *wavemap,
                       double dispersion, cpl_image *skymap)
{
    const char *func = "mos_sky_map";

    cpl_vector **vector;
    double       firstLambda, lastLambda;
    double       lambda, lambda1, lambda2;
    double       value, value1, value2;
    float        min, max;
    int         *count;
    int          nbin, bin;
    int          nx, ny, npix;
    int          i, j;

    cpl_table   *sky;
    double      *sky_spectrum;
    float       *data;
    float       *sdata;
    float       *kdata;
    double      *wdata;


    if (spectra == NULL || wavemap == NULL || skymap == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }
    
    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);
    npix = nx * ny;

    if (nx != cpl_image_get_size_x(wavemap) ||
        ny != cpl_image_get_size_y(wavemap) ||
        nx != cpl_image_get_size_x(skymap) ||
        ny != cpl_image_get_size_y(skymap)) {
        cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
        return NULL;
    }


    /*
     * Find bluest and reddest wavelengths in the whole image
     */

    data = cpl_image_get_data(wavemap);

    for (i = 0; i < npix; i++) {
        if (data[i] > 1.0) {
            min = max = data[i];
            j = i+1;
            break;
        }
    }

    for (i = j; i < npix; i++) {
        if (data[i] < 1.0)      /* Impossible wavelength */
            continue;
        if (min > data[i])
            min = data[i];
        if (max < data[i])
            max = data[i];
    }

    firstLambda = min;
    lastLambda = max;


    /*
     * Determine length of median spectrum
     */

    nbin = (lastLambda - firstLambda) / dispersion;

    /*
     * Count how many values will be found for each spectral bin.
     * The ith bin has a wavelength range from firstLambda + i*dispersion
     * (inclusive) to firstLambda + (i+1)*dispersion (exclusive), and
     * it is assigned to its central wavelength.
     */

    count = cpl_calloc(nbin, sizeof(int));

    data = cpl_image_get_data(wavemap);

    for (i = 0; i < npix; i++) {
        if (data[i] < 1.0)
            continue;
        bin = (data[i] - firstLambda) / dispersion;
        if (bin < nbin)                               /* Safer */
            count[bin]++;
    }


    /*
     * Allocate an array of vectors with the appropriate size, to
     * contain a list of all the spectral pixels values. At the same
     * time, reset the array of counters (because we will have to
     * count again...).
     */

    vector = cpl_calloc(nbin, sizeof(cpl_vector *));
    for (i = 0; i < nbin; i++) {
        if (count[i])
            vector[i] = cpl_vector_new(count[i]);
        else
            vector[i] = NULL;
        count[i] = 0;
    }


    /*
     * Read the wavemap and the spectral images, and add the data values
     * to the appropriate wavelength bins
     */

    data  = cpl_image_get_data(wavemap);
    sdata = cpl_image_get_data(spectra);

    for (i = 0; i < npix; i++) {
        if (data[i] < 1.0)
            continue;
        bin = (data[i] - firstLambda) / dispersion;
        if (bin < nbin) {                             /* Safer */
            cpl_vector_set(vector[bin], count[bin], sdata[i]);
            count[bin]++;
        }
    }


    /*
     * Compute the median flux for each wavelength bin, and destroy
     * at the same time the used vectors
     */

    sky_spectrum = cpl_calloc(nbin, sizeof(double));
    for (i = 0; i < nbin; i++) {
        if (vector[i]) {
            sky_spectrum[i] = cpl_vector_get_median_const(vector[i]);
            cpl_vector_delete(vector[i]);
        }
    }

    cpl_free(vector);


    /*
     * Here possible gaps in the final spectrum should be filled
     * by interpolation
     */

    /* ... */

    /*
     * Create the output table
     */

    sky = cpl_table_new(nbin);
    cpl_table_new_column(sky, "wavelength", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(sky, "wavelength", "pixel");
    cpl_table_wrap_double(sky, sky_spectrum, "sky");
    cpl_table_wrap_int(sky, count, "npoints");
    for (i = 0; i < nbin; i++)
        cpl_table_set_double(sky, "wavelength", i, 
                             firstLambda + (i+0.5)*dispersion);


    /*
     * Fill the sky map
     */

    data  = cpl_image_get_data(wavemap);
    sdata = cpl_image_get_data(spectra);
    kdata = cpl_image_get_data(skymap);
    wdata = cpl_table_get_data_double(sky, "wavelength");

    for (i = 0; i < npix; i++) {

        /*
         * Currently based on linear interpolation
         */

        lambda = data[i];
        if (lambda < 1.0)
            continue;
        bin = (lambda - firstLambda) / dispersion;
        lambda1 = wdata[bin];
        value1 = sky_spectrum[bin];
        if (lambda1 < lambda) {
            bin++;
            if (bin < nbin) {
                lambda2 = wdata[bin];
                value2  = sky_spectrum[bin];
                value   = ((lambda2 - lambda)*value1 
                        +  (lambda - lambda1)*value2) / dispersion;
            }
            else {
                value = value1;
            }
        }
        else {
            if (bin > 0) {
                bin--;
                lambda2 = lambda1;
                value2  = value1;
                lambda1 = wdata[bin];
                value1  = sky_spectrum[bin];
                value   = ((lambda2 - lambda)*value1 
                        +  (lambda - lambda1)*value2)/dispersion;
            }
            else {
                value = value1;
            }
        }
        kdata[i] = value;
    }

    return sky;

}


/**
 * @brief
 *   Local determination of sky
 *   
 * @param spectra     Rectified image of scientific spectra
 * @param slits       Slits positions table
 *  
 * @return Image of (rectified) sky
 *  
 * Each slit spectrum contained in the input image, @em spectra, 
 * is assumed to contain at least 50% of its pixels on the sky. 
 * The median level at each wavelength is determined, and the
 * sky map generated in this way is returned.
 */

cpl_image *mos_sky_local_old(cpl_image *spectra, cpl_table *slits)
{
    const char *func = "mos_sky_local_old";

    cpl_image *exslit;
    cpl_image *sky;
    cpl_image *skymap;
    float     *data;
    float     *sdata;
    int        nx, ny;
    int        xlow, ylow, xhig, yhig;
    int        nslits;
    int       *slit_id;
    int       *position;
    int       *length;
    int        i, j, k;


    if (spectra == NULL) {
        cpl_msg_error(func, 
                      "A scientific rectified spectral image must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (slits == NULL) {
        cpl_msg_error(func, "A slits position table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nslits   = cpl_table_get_nrow(slits);
    slit_id  = cpl_table_get_data_int(slits, "slit_id");
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);

    skymap = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);

    xlow = 1;
    xhig = nx;
    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        /*
         * Define the extraction boundaries. We DON'T write:
         *
         * ylow = position[i];
         * yhig = ylow + length[i];
         *
         * because the cpl_image pixels are counted from 1, and because in
         * cpl_image_extract() the coordinates of the last pixel are inclusive.
         */

        ylow = position[i] + 1;
        yhig = ylow + length[i] - 1;

        exslit = cpl_image_extract(spectra, xlow, ylow, xhig, yhig);
        sky    = cpl_image_collapse_median_create(exslit, 0, 0, 1);
        cpl_image_delete(exslit);

        data   = cpl_image_get_data(skymap);
        data  += nx * position[i];

        for (j = 0; j < length[i]; j++) {
            sdata  = cpl_image_get_data(sky);
            for (k = 0; k < nx; k++) {
                *data++ = *sdata++;
            }
        }

        cpl_image_delete(sky);
    }

    return skymap;

}


/**
 * @brief
 *   Local determination of sky
 *   
 * @param spectra     Rectified image of scientific spectra
 * @param slits       Slits positions table
 * @param order       Order of the polynomial fitting the sky.
 *  
 * @return Image of (rectified) sky
 *  
 * The median level at each wavelength is determined and subtracted
 * from the data. The position and extension of the objects is detected
 * and the rest of the pixels are flagged as sky, and used for the
 * final determination of the sky level at each wavelength. If the
 * @em order of the sky-fitting polynomial is zero the median sky
 * level is determined, otherwise a fit with outliers rejection is
 * applied. The sky map generated in this way is returned.
 */

cpl_image *mos_sky_local(cpl_image *spectra, cpl_table *slits, int order)
{
    const char *func = "mos_sky_local";

    char        name[MAX_COLNAME];

    cpl_polynomial *fit;
    cpl_vector     *points;
    cpl_vector     *values;
    cpl_vector     *keep_points;
    cpl_vector     *keep_values;
    cpl_image      *exslit;
    cpl_image      *sky;
    cpl_image      *subtracted;
    cpl_image      *profile;
    cpl_image      *skymap;
    cpl_table      *objects;
    float          *data;
    float          *sdata;
    float          *xdata;
    double         *vdata;
    double         *pdata;
    double          median;
    int             nx, ny;
    int             xlow, ylow, xhig, yhig;
    int             nslits;
    int            *slit_id;
    int            *position;
    int            *length;
    int            *is_sky;
    int             nsky, nbad;
    int             maxobjects;
    int             margin = 3;
    int             radius = 6;
    int             i, j, k;


    if (spectra == NULL) {
        cpl_msg_error(func, 
                      "A scientific rectified spectral image must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (slits == NULL) {
        cpl_msg_error(func, "A slits position table must be given");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (order < 0) {
        cpl_msg_error(func, "Invalid fit order");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nslits   = cpl_table_get_nrow(slits);
    slit_id  = cpl_table_get_data_int(slits, "slit_id");
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);

    skymap = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);

    xlow = 1;
    xhig = nx;
    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        /*
         * Define the extraction boundaries. We DON'T write:
         *
         * ylow = position[i];
         * yhig = ylow + length[i];
         *
         * because the cpl_image pixels are counted from 1, and because in
         * cpl_image_extract() the coordinates of the last pixel are inclusive.
         */

        ylow = position[i] + 1;
        yhig = ylow + length[i] - 1;

        exslit = cpl_image_extract(spectra, xlow, ylow, xhig, yhig);
        sky    = cpl_image_collapse_median_create(exslit, 0, 0, 1);
        cpl_image_delete(exslit);

        data   = cpl_image_get_data(skymap);
        data  += nx * position[i];

        for (j = 0; j < length[i]; j++) {
            sdata  = cpl_image_get_data(sky);
            for (k = 0; k < nx; k++) {
                *data++ = *sdata++;
            }
        }

        cpl_image_delete(sky);
    }


    /*
     * Preliminary sky-subtracted image
     */

    subtracted = cpl_image_duplicate(spectra);
    cpl_image_subtract(subtracted, skymap);
    cpl_image_delete(skymap);


    /*
     * Detect objects positions in all slits
     */

    objects = cpl_table_duplicate(slits);
    profile = mos_detect_objects(subtracted, objects, margin, radius, 0, 5);
    cpl_image_delete(profile);
    cpl_image_delete(subtracted);


    /*
     * Flag the sky pixels. Note that maxobjects is intentionally 
     * the max number of objects increased by one.
     */

    maxobjects = 1;
    snprintf(name, MAX_COLNAME, "object_%d", maxobjects);
    while (cpl_table_has_column(objects, name)) {
        maxobjects++;
        snprintf(name, MAX_COLNAME, "object_%d", maxobjects);
    }

    is_sky = cpl_calloc(ny, sizeof(int));

    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        ylow = position[i] + margin;
        yhig = position[i] + length[i] - margin;

        for (j = ylow; j < yhig; j++)
            is_sky[j] = 1;

        for (j = 1; j < maxobjects; j++) {
            snprintf(name, MAX_COLNAME, "object_%d", j);
            if (cpl_table_is_valid(objects, name, i)) {
                snprintf(name, MAX_COLNAME, "start_%d", j);
                ylow = cpl_table_get_int(objects, name, i, NULL);
                snprintf(name, MAX_COLNAME, "end_%d", j);
                yhig = cpl_table_get_int(objects, name, i, NULL);
                for (k = ylow; k <= yhig; k++)
                    is_sky[k] = 0;
            }
        }


        /*
         * Eliminate isolated sky points
         */

        ylow = position[i] + margin + 1;
        yhig = position[i] + length[i] - margin - 1;

        for (j = ylow; j < yhig; j++)
            if (is_sky[j])
                if (is_sky[j-1] == 0 && is_sky[j+1] == 0)
                    is_sky[j] = 0;

    }


    /*
     * Determination of the sky map
     */

    skymap = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);

    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        ylow = position[i];
        yhig = ylow + length[i];

        nsky = 0;
        for (j = ylow; j < yhig; j++)
            if (is_sky[j])
                nsky++;

        if (nsky > order + 1) {
            if (order) {
                points = cpl_vector_new(nsky);
                nsky = 0;
                for (j = ylow; j < yhig; j++) {
                    if (is_sky[j]) {
                        cpl_vector_set(points, nsky, j);
                        nsky++;
                    }
                }

                exslit = cpl_image_extract(spectra, 1, ylow+1, nx, yhig);
                xdata = cpl_image_get_data(exslit);
                values = cpl_vector_new(nsky);

                for (j = 0; j < nx; j++) {
                    nsky = 0;
                    for (k = ylow; k < yhig; k++) {
                        if (is_sky[k]) {
                            cpl_vector_set(values, nsky, xdata[j+(k-ylow)*nx]);
                            nsky++;
                        }
                    }

                    /*
                     * Eliminate obvious outliers
                     */

                    median = cpl_vector_get_median_const(values);
                    vdata = cpl_vector_get_data(values);
                    pdata = cpl_vector_get_data(points);
                    nbad = 0;
                    for (k = 0; k < nsky; k++) {
                        if (fabs(vdata[k] - median) < 100) {
                            if (nbad) {
                                vdata[k-nbad] = vdata[k];
                                pdata[k-nbad] = pdata[k];
                            }
                        }
                        else
                            nbad++;
                    }

                    if (nsky == nbad)
                        continue;

                    if (nbad && nsky - nbad > order + 1) {
                        keep_values = values;
                        keep_points = points;
                        values = cpl_vector_wrap(nsky-nbad, vdata);
                        points = cpl_vector_wrap(nsky-nbad, pdata);
                    }

                    if (nsky - nbad > order + 1) {

                        fit = cpl_polynomial_fit_1d_create(points, values, 
                                                           order, NULL);

                        if (fit) {
                            for (k = ylow; k < yhig; k++) {
                                xdata[j+(k-ylow)*nx] = 
                                         cpl_polynomial_eval_1d(fit, k, NULL);
                            }

                            cpl_polynomial_delete(fit);
                        }
                        else
                            cpl_error_reset();
                    }
                    else {
                        for (k = 0; k < nsky; k++) {
                            xdata[j+k*nx] = median;
                        }
                    }

                    if (nbad && nsky - nbad > order + 1) {
                        cpl_vector_unwrap(values);
                        cpl_vector_unwrap(points);
                        values = keep_values;
                        points = keep_points;
                    }

                    if (nbad) {
                        nsky = 0;
                        for (k = ylow; k < yhig; k++) {
                            if (is_sky[k]) {
                                cpl_vector_set(points, nsky, k);
                                nsky++;
                            }
                        }
                    }

                }

                cpl_vector_delete(values);
                cpl_vector_delete(points);

                cpl_image_copy(skymap, exslit, 1, ylow+1);
                cpl_image_delete(exslit);

            }
            else {
                exslit = cpl_image_extract(spectra, 1, ylow+1, nx, yhig);
                xdata = cpl_image_get_data(exslit);
                values = cpl_vector_new(nsky);

                for (j = 0; j < nx; j++) {
                    nsky = 0;
                    for (k = ylow; k < yhig; k++) {
                        if (is_sky[k]) {
                            cpl_vector_set(values, nsky, xdata[j+(k-ylow)*nx]);
                            nsky++;
                        }
                    }

                    median = cpl_vector_get_median_const(values);

                    for (k = ylow; k < yhig; k++)
                        xdata[j+(k-ylow)*nx] = median;

                }

                cpl_vector_delete(values);

                cpl_image_copy(skymap, exslit, 1, ylow+1);
                cpl_image_delete(exslit);
            }
        }
        else
            cpl_msg_warning(func, "Too few sky points in slit %d", i + 1);
    }

    cpl_free(is_sky);

    return skymap;

}


/**
 * @brief
 *   Remove cosmic rays from sky-subtracted CCD spectral exposure.
 *
 * @param image      Input image to be cleaned (in ADU)
 * @param gain       Inverse gain factor (e-/ADU)
 * @param threshold  Threshold for cosmics detection, given in noise sigmas
 * @param ratio      Ratio for discrimination between objects and cosmics
 *
 * @return @c CPL_ERROR_NONE in case of success
 *
 * If @em treshold is negative, it is assigned the value 4.0. 
 * If @em ratio is negative, it is assigned the value 2.0.
 * The algorithm used is the same of the MIDAS command FILTER/COSMIC.
 * This function should be used to process not-rectified spectral 
 * exposures from where the sky spectrum was already subtracted.
 * Cosmics are not cleaned if either x or y size of the image is
 * not greater than 3: in this case the function returns without
 * setting an error.
 */

cpl_error_code mos_clean_cosmics(cpl_image *image, float gain,
                                 float threshold, float ratio)
{
    const char *func = "mos_clean_cosmics";

    cpl_image  *smoothImage;
    cpl_table  *table;
    cpl_matrix *kernel;
    int        *xdata;
    int        *ydata;
    float      *idata;
    float      *sdata;
    float       sigma, sum, value, smoothValue;
    double      noise;
    int         count;
    float       fMax;
    int         iMin, iMax, jMin, jMax, iPosMax, jPosMax;
    int         xLen;
    int         yLen;
    int         nPix;
    int         first = 1;  /* position of first cosmic ray candidate
                               encountered while scanning the image */
    int         pos, i, j, k, l, ii, jj, iii = 0, jjj = 0;
    int         numCosmic = 0;
    int         found, foundContiguousCandidate;
    int        *cosmic;
  

    if (image == NULL)
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);


    /*
     *  "cosmic" is a flags holder (initialized to zero):
     *
     *           -1 = candidate for cosmic ray
     *            0 = not a cosmic
     *            1 = a cosmic ray
     *            2 = member of current group of contiguous candidates
     *            3 = examined member of current group
     */

    xLen = cpl_image_get_size_x(image);
    yLen = cpl_image_get_size_y(image);

    if (xLen < 4 || yLen < 4)
        return CPL_ERROR_NONE;

    nPix = xLen * yLen;

    /*
     * Noise estimation from negative offsets in image. Note that this
     * assumes that the background level (skyLevel) has already been 
     * subtracted from the data. In this way we estimate the noise due 
     * to detector readout and to the background signal level (before 
     * it were removed). Theoretically this is given by 
     *
     *        noise = sqrt(ron^2 + skyLevel/gain)
     *
     * where ron is the read-out-noise. To this we will sum the noise 
     * contribution due to any increase of the signal above the background
     * by an amount scienceLevel. Theoretically the total noise is given by
     *
     *        totalNoise = sqrt(ron^2 + (skyLevel+scienceLevel)/gain)
     *
     * that is
     *
     *        totalNoise = sqrt(noise^2 + scienceLevel/gain)
     *
     */

    idata = cpl_image_get_data(image);
    noise = 0.0;
    count = 0;

    for (i = 0; i < nPix; i++) {
        if (idata[i] < -0.00001) {
            noise -= idata[i];
            count++;
        }
    }

    noise /= count;
    noise *= 1.25;       /* Factor to convert average deviation to sigma */

    cosmic = cpl_calloc(nPix, sizeof(int));

    if (threshold < 0.)
        threshold = 4.0;
    if (ratio < 0.)
        ratio = 2.0;

    kernel = cpl_matrix_new(3, 3);
    cpl_matrix_fill(kernel, 1.0);
    cpl_matrix_set(kernel, 1, 1, 0.0);
    smoothImage = cpl_image_filter_median(image, kernel);
    cpl_matrix_delete(kernel);
    
    /*
     *  Loop on images pixels, searching for cosmic rays candidates.
     *  Border pixels are currently excluded (they cannot contain
     *  candidates), to avoid that the search for groups of contiguous
     *  pixels would ever go out of image boundaries. In future we may
     *  overcome this limit, adding an appropriate check when contiguous
     *  pixels are searched.
     */

    sdata = cpl_image_get_data(smoothImage);

    for (j = 1; j < yLen - 1; j++) {
        for (i = 1; i < xLen - 1; i++) {
            value = idata[i + j * xLen];
            smoothValue = sdata[i + j * xLen];
            if (smoothValue < 1.0)
                smoothValue = 1.0;
            sigma = sqrt(noise * noise + smoothValue / gain);
            if (value - smoothValue >= threshold * sigma) 
                cosmic[i + j * xLen] = -1;
        }
    }

    cpl_image_delete(smoothImage);


    /*
     *  Search for groups of contiguous cosmic rays candidates.
     */

    do {
        found = 0;
        for (pos = first; pos < nPix; pos++) {
            if (cosmic[pos] == -1) {
                cosmic[pos] = 2;         /*  Candidate found.  */
                i = pos % xLen;          /*  Its coordinates.  */
                j = pos / xLen;
                first = pos;
                first++;      /* ???  really necessary? */
                found = 1;
                break;
            }
        }

        if (found) {

            /*
             *  Determine new group of contiguous cosmic rays candidates.
             *  Initialize the working box boundaries, iMin, iMax, jMin, jMax, 
             *  and the value of the max pixel and its position, fMax, iPosMax,
             *  jPosMax.
             */

            iMin = iMax = iPosMax = i;
            jMin = jMax = jPosMax = j;
            fMax = idata[i + j * xLen];

            do {
                foundContiguousCandidate = 0;
                for (l = 0; l <= 1; l++) {
                    for (k = 0; k <= 1; k++) {

                        /*
                         *  Looping on 4 pixels to North, East, South and West
                         */

                        ii = i + k - l;
                        jj = j + k + l - 1;
                        if (cosmic[ii + jj * xLen] == -1) {
                            foundContiguousCandidate = 1;
                            cosmic[ii + jj * xLen] = 2;
                                        /* Candidate belongs to current group */
                            iii = ii;   /* Keep its position */
                            jjj = jj;

                            /*
                             * Upgrade search box
                             */

                            if (ii < iMin)
                                iMin = ii;
                            if (ii > iMax)
                                iMax = ii;
                            if (jj < jMin)
                                jMin = jj;
                            if (jj > jMax)
                                jMax = jj;

                            if (idata[ii + jj * xLen] > fMax) {
                                fMax = idata[ii + jj * xLen];
                                iPosMax = ii;
                                jPosMax = jj;
                            }
                        }
                    }
                }

                /*
                 *  We are done exploring the "cross". Now mark as "examined"
                 *  the current candidate (at the center of the cross):
                 */

                cosmic[i + j * xLen] = 3; /* It may probably be set to 1 now */

                if (foundContiguousCandidate) {

                    /*
                     * Pass (arbitrarily) the coordinates of the LAST found 
                     * candidate
                     */

                    i = iii;
                    j = jjj;

                    /* 
                     * Skip the rest, continue loop on new candidate 
                     */

                    continue; 
                }


                /*
                 *  Look for leftovers in the (growing!) search box
                 */

                for (l = jMin; l <= jMax; l++) {
                    for (k = iMin; k <= iMax; k++) {
                        if (cosmic[k + l * xLen] == 2) {
                            i = k;
                            j = l;
                            foundContiguousCandidate = 1;
                            break;
                        }
                    }
                    if (foundContiguousCandidate) 
                        break;
                }
            } while (foundContiguousCandidate);


            /*
             *  No more contiguous candidates are found. Decide now
             *  whether the current group is a cosmic ray or not.
             */

            sum = 0.;                /* Sum of 8 pixels around max position */
            for (l = -1; l <= 1; l++) {
                for (k = -1; k <= 1; k++) {
                    if (l != 0 || k != 0) {
                        sum += idata[iPosMax + k + (jPosMax + l) * xLen];
                    }
                }
            }

            sum /= 8.;
            if (fMax > ratio * sum) {
                for (l = jMin - 1; l <= jMax + 1; l++) {
                    for (k = iMin - 1; k <= iMax + 1; k++) {
                        if (cosmic[k + l * xLen] == 3) {
                            cosmic[k + l * xLen] = 1;
                            numCosmic++;
                        }
                    }
                }
            }
            else {
                for (l = jMin - 1; l <= jMax + 1; l++) {
                    for (k = iMin - 1; k <= iMax + 1; k++) {
                        if (cosmic[k + l * xLen] != -1) {
                            if (cosmic[k + l * xLen] == 1) 
                                numCosmic--;
                            cosmic[k + l * xLen] = 0;
                        }
                    }
                }
            }
        }
    } while (found);


    /*
     *  Prepare table containing cosmic rays coordinates. 
     */

    table = cpl_table_new(numCosmic);
    cpl_table_new_column(table, "x", CPL_TYPE_INT);
    cpl_table_new_column(table, "y", CPL_TYPE_INT);
    cpl_table_set_column_unit(table, "x", "pixel");
    cpl_table_set_column_unit(table, "y", "pixel");
    xdata = cpl_table_get_data_int(table, "x");
    ydata = cpl_table_get_data_int(table, "y");

    for (pos = 0, i = 0; pos < nPix; pos++) {
        if (cosmic[pos] == 1) {
            xdata[i] = (pos % xLen);
            ydata[i] = (pos / xLen);
            i++;
        }
    }

    mos_clean_bad_pixels(image, table, 1);

    cpl_free(cosmic);
    cpl_table_delete(table);

    return CPL_ERROR_NONE;

}


cpl_error_code mos_clean_bad_pixels(cpl_image *image, cpl_table *table,
                                    int spectral)
{
    const char *func = "mos_clean_cosmics";
 
    float       *idata;
    int         *isBadPix;
    int          i, j, k, d;
    int          xlen, ylen, totPix;
    int          nBadPixels = 0;
    int          sign, foundFirst;
    int         *xValue = NULL;
    int         *yValue = NULL;
    float        save = 0.;
    double       sumd;
    int          cx, cy;
    int          nPairs;
    float        estimate[4];
    int          sx[] = {0, 1, 1, 1};
    int          sy[] = {1,-1, 0, 1};
    int          searchHorizon = 100;
    int          percent = 15;


    if (image == NULL || table == NULL)
        return cpl_error_set(func, CPL_ERROR_NULL_INPUT);

    if (1 != cpl_table_has_column(table, "x"))
        return cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);

    if (1 != cpl_table_has_column(table, "y"))
        return cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);

    if (CPL_TYPE_INT != cpl_table_get_column_type(table, "x"))
        return cpl_error_set(func, CPL_ERROR_INVALID_TYPE);

    if (CPL_TYPE_INT != cpl_table_get_column_type(table, "y"))
        return cpl_error_set(func, CPL_ERROR_INVALID_TYPE);

    nBadPixels = cpl_table_get_nrow(table);

    if (nBadPixels) {
        xlen = cpl_image_get_size_x(image);
        ylen = cpl_image_get_size_y(image);
        idata = cpl_image_get_data(image);
        totPix = xlen * ylen;
        if (((float) nBadPixels) / ((float) totPix) < percent/100.) {
            isBadPix = cpl_calloc(totPix, sizeof(int));
        }
        else {
            cpl_msg_warning(func, "Too many bad pixels (> %d%%): "
                            "skip bad pixel correction", percent);
            return CPL_ERROR_ILLEGAL_INPUT;
        }
    }
    else {
        cpl_msg_debug(func, "No pixel values to interpolate");
        return CPL_ERROR_NONE;
    }

    xValue = cpl_table_get_data_int(table, "x");
    yValue = cpl_table_get_data_int(table, "y");

    for (i = 0; i < nBadPixels; i++)
        isBadPix[xValue[i] + yValue[i] * xlen] = 1;

    for (i = 0; i < nBadPixels; i++) {

        /*
         *  Search for the closest good pixel along the 4 fundamental 
         *  directions (in both senses):
         *                            \ | /
         *                             \|/
         *                           --- ---
         *                             /|\
         *                            / | \
         *
         *  Then collect pairs of values to interpolate linearly.
         */

        nPairs = 0;
        for (j = 0; j < 4; j++) {

            if (spectral) /* Just horizontal interpolation for spectral data */
                if (j != 2)
                    continue;

            estimate[nPairs] = 0.;  /* Pairs interpolations are stored here */
            sumd = 0.;
            foundFirst = 0;
            for (k = 0; k < 2; k++) {
                sign = 2 * k - 1;
                d = 0;
                cx = xValue[i];
                cy = yValue[i];
                do {
                    cx += sign * sx[j];
                    cy += sign * sy[j];
                    if (cx < 0 || cx >= xlen || cy < 0 || cy >= ylen) 
                        break;
                    d++;
                } while (isBadPix[cx + cy * xlen] && d < searchHorizon);

                if (cx >= 0 && cx < xlen && 
                    cy >= 0 && cy < ylen && d < searchHorizon) {

                    /*
                     *  In this block is cripted the linear interpolation...
                     */

                    save = idata[cx + cy * xlen];
                    estimate[nPairs] += save / d;
                    sumd += 1. / (double) d;
                    if (k) {
                        estimate[nPairs] /= sumd;
                        nPairs++;
                    }
                    else {
                        foundFirst = 1;
                    }
                }
                else {

                    /*
                     * Image borders were crossed, incomplete pair of values
                     */

                    if (k) {
                        if (foundFirst) {
                            estimate[nPairs] = save;
                            nPairs++;
                        }
                    }
                }
            }
        }

        /*
         * Replace pixel value of the input image, corresponding to
         * the current bad pixel, with the median of the estimates
         * resulted from the 4 linear interpolations.
         */

        if (nPairs > 2) {
            idata[xValue[i] + yValue[i] * xlen] = 
                               cpl_tools_get_median_float(estimate, nPairs);
        }
        else if (nPairs == 2) {
            idata[xValue[i] + yValue[i] * xlen] =
                                            (estimate[0] + estimate[1]) / 2.;
        }
        else if (nPairs == 1) {
            idata[xValue[i] + yValue[i] * xlen] = estimate[0];
        }
        else {
            cpl_msg_debug(func, "Cannot correct bad pixel %d,%d\n",
                          xValue[i], yValue[i]);
        }
    }

    cpl_free(isBadPix);

    return CPL_ERROR_NONE;
}


/**
 * @brief
 *   Create coordinate map from spectral curvature table
 *
 * @param spectra     CCD image of spectra
 * @param slits       Slits positions on the CCD
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param reference   Reference wavelength
 * @param blue        Start lambda for spatial remapping
 * @param red         End lambda for spatial remapping
 * @param dispersion  Mean spectral dispersion
 *
 * @return Coordinate map
 * 
 * The input @em spectra image is expected to be oriented with
 * horizontal dispersion direction and red wavelengths on the right
 * side, and it should be of type @c CPL_TYPE_FLOAT. The @em slits
 * table should be the output of either the function @c mos_identify_slits()
 * (if available) or the function @c mos_locate_spectra(). The @em polytraces 
 * table is the output of the function @c mos_poly_trace().
 *
 * The returned image image has the same size of the input @em spectra 
 * image, and it will consist of pixels having the value of their spatial
 * coordinate along the slit they belong to, or the value zero if the
 * pixels do not belong to any spectrum. The original y pixel positions 
 * (at a given x) are assigned the value of their distance from the top 
 * spectral trace measured in CCD pixels.
 */

cpl_image *mos_spatial_map(cpl_image *spectra, cpl_table *slits,
                           cpl_table *polytraces, double reference,
                           double blue, double red, double dispersion)
{
    const char *func = "mos_spatial_map";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;
    cpl_image      *calibration;
    float          *data;
    double          top, bot;
    double          coeff;
    double          ytop, ybot;
    double          ypos, yfra;
    double          factor;
    int             yint, yprev;
    int             nslits;
    int             npseudo;
    int            *slit_id;
    int            *length;
    int             nx, ny;
    int             pixel_above, pixel_below, refpixel, start_pixel, end_pixel;
    int             missing_top, missing_bot;
    int             null;
    int             order;
    int             i, j;
    cpl_size        k;


    if (spectra == NULL || slits == NULL || polytraces == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (red - blue < dispersion) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);

    calibration = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
    data = cpl_image_get_data(calibration);

    length  = cpl_table_get_data_int(slits, "length");
    nslits  = cpl_table_get_nrow(slits);
    slit_id = cpl_table_get_data_int(slits, "slit_id");
    order   = cpl_table_get_ncol(polytraces) - 2;

    /*
     * The spatial resampling is performed for a certain number of 
     * pixels above and below the position of the reference wavelength:
     */

    pixel_above = STRETCH_FACTOR * (red - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - blue) / dispersion;

    for (i = 0; i < nslits; i++) {
        
        if (length[i] == 0)
            continue;

        /*
         * Note that the x coordinate of the reference pixels on the CCD
         * is taken arbitrarily at the top end of each slit. This wouldn't
         * be entirely correct in case of curved slits, or in presence of
         * heavy distortions: in such cases the spatial resampling is
         * really performed across a wide range of wavelengths. But
         * the lag between top and bottom spectral curvature models 
         * would introduce even in such cases negligible effects on
         * the spectral spatial resampling.
         */

        refpixel = cpl_table_get_double(slits, "xtop", i, NULL);

        start_pixel = refpixel - pixel_below;
        if (start_pixel < 0)
            start_pixel = 0;

        end_pixel = refpixel + pixel_above;
        if (end_pixel > nx)
            end_pixel = nx;

        /*
         * Recover from the table of spectral curvature coefficients
         * the curvature polynomials.
         */

        missing_top = 0;
        polytop = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i, &null);
            if (null) {
                cpl_polynomial_delete(polytop);
                missing_top = 1;
                break;
            }
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        missing_bot = 0;
        polybot = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i+1, &null);
            if (null) {
                cpl_polynomial_delete(polybot);
                missing_bot = 1;
                break;
            }
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        if (missing_top && missing_bot) {
            cpl_msg_warning(func, "Spatial map, slit %d was not traced!", 
                            slit_id[i]);
            continue;
        }

        /*
         * In case just one of the two edges was not traced, the other
         * edge curvature model is duplicated and shifted to the other
         * end of the slit: better than nothing!
         */

        if (missing_top) {
            cpl_msg_warning(func, "Upper edge of slit %d was not traced: "
                            "the spectral curvature of the lower edge "
                            "is used instead.", slit_id[i]);
            polytop = cpl_polynomial_duplicate(polybot);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polybot, &k);
            coeff += ytop - ybot;
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        if (missing_bot) {
            cpl_msg_warning(func, "Lower edge of slit %d was not traced: "
                            "the spectral curvature of the upper edge "
                            "is used instead.", slit_id[i]);
            polybot = cpl_polynomial_duplicate(polytop);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polytop, &k);
            coeff -= ytop - ybot;
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        top = cpl_polynomial_eval_1d(polytop, refpixel, NULL);
        bot = cpl_polynomial_eval_1d(polybot, refpixel, NULL);
        npseudo = ceil(top-bot) + 1;

        if (npseudo < 1) {
            cpl_polynomial_delete(polytop);
            cpl_polynomial_delete(polybot);
            cpl_msg_warning(func, "Slit %d was badly traced: no extraction!",
                            slit_id[i]);
            continue;
        }

        for (j = start_pixel; j < end_pixel; j++) {
            top = cpl_polynomial_eval_1d(polytop, j, NULL);
            bot = cpl_polynomial_eval_1d(polybot, j, NULL);
            factor = (top-bot)/npseudo;
            for (k = 0; k <= npseudo; k++) {
                ypos = top - k*factor;
                yint = ypos;
                yfra = ypos - yint;
                if (yint >= 0 && yint < ny-1) {
                    data[j + nx*yint] = (top-yint)/factor;
                    if (k) {

                        /*
                         * This is added to recover lost pixels on
                         * the CCD image (pixels are lost because
                         * the CCD pixels are less than npseudo+1).
                         */

                        if (yprev - yint > 1) {
                            data[j + nx*(yint+1)] = (top-yint-1)/factor;
                        }
                    }
                }
                yprev = yint;
            }
        }
        cpl_polynomial_delete(polytop);
        cpl_polynomial_delete(polybot);
    }

    return calibration;
}


/**
 * @brief
 *   Detect objects in rectified scientific frame
 *
 * @param image      Rectified image of scientific spectra
 * @param slits      Table with slits positions
 * @param margin     Number of pixels to exclude at slits edges
 * @param maxradius  Maximum extraction radius
 * @param conradius  Contamination radius
 *
 * @return @c Slits mean flux spatial profile
 *
 * This function upgrades the input @em slits position table with the
 * positions of the objects detected within each slit, and their 
 * corresponding extraction intervals. The object positions are listed 
 * in columns labeled "object_1", "object_2", etc., while the start and
 * end positions of the extraction interval are marked with "start_1",
 * "start_2", ... and "end_1", "end_2", etc., where the pixel coordinate
 * is counted starting from the bottom position of the slit. Objects
 * are not searched in slits having less than 10 valid pixels.
 *
 * The algorithm applied is based on the relative peak intensity of
 * each detected object. First of all, a peak is identified by a
 * positive value that is preceded and followed by two lower positive 
 * values that decrease with distance. Also the first and the last
 * (valid) spatial pixel of a slit is considered a peak, if followed
 * or preceded by decreasing pixels, to keep it into account in the
 * computation of the contaminations, but it is never extracted in
 * the end. Each peak is compared with all the other peaks to determine 
 * if this peak is contaminated by any of the others. Indicating with 
 * L_o the peak value of the examined peak and with L the peak value 
 * of another peak, the quantity
 * @code
 *            S = C * (L / L_o)
 * @endcode
 * is computed, where C is the indicated contamination radius, 
 * @em conradius, that represents the minimum distance at which 
 * two point-like object of equal luminosity can stay without
 * contaminating each others (a typical value for this parameter
 * may be 16 pixels). If the distance between the two peaks is
 * less than S, the examined peak is flagged as contaminated,
 * and is excluded from the @em final list of detected objects.
 * This empirical formula has the effect of assigning a larger
 * contamination radius to relatively brighter objects with
 * respect to dimmer ones. With the final list of object positions
 * the extraction intervals are determined in the following way:
 * for each pair of consecutive peaks, an intermediate positions
 * is determined with the inverse baricenter formula (defining
 * the point of minimal reciprocal contamination):
 * @code
 *            B_i = (P_i * L_j + P_j * L_i) / (L_i + L_j)
 * @endcode
 * where P_i is the position of the i-th peak, L_i its peak value,
 * and j = i + 1. The position of the lower limit of the first object 
 * is set at the top border of the slit, excluding the number of
 * pixels indicated by the @em margin argument. Analogously, the 
 * position of the upper limit of the last object is set at the 
 * bottom border of the slit, excluding the same number of pixels.
 * Finally, the extraction borders exceeding @em maxradius are
 * corrected accordingly.
 */

cpl_image *mos_detect_objects(cpl_image *image, cpl_table *slits, int margin,
                              int maxradius, int conradius, double noise)
{
    const char *func = "mos_detect_objects";

    cpl_image  *profile;
    float      *pdata;
    float      *p;

    cpl_image  *valid;
    int        *ivalid;

    float      *data;

    char        name[MAX_COLNAME];

    int         nslits;
    int         npeaks;
    int         nobjects, objpos, totobj;
    int         maxobjects;
    int        *position;
    int        *length;
    int        *reject;
    double     *place;
    double     *bright;
    double      mindistance;
    int         pos, count;
    int         up;
    int         low, hig;
    int         row;
    int         nx, ny, npix;
    int         i, j, k;

    const int   min_pixels = 10;

    
    if (cpl_error_get_code() != CPL_ERROR_NONE)
        return NULL;
 
    if (image == NULL || slits == NULL) { 
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (margin < 0)
        margin = 0;

    if (maxradius < 0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (conradius < 0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    nslits   = cpl_table_get_nrow(slits);
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");
    nx       = cpl_image_get_size_x(image);
    ny       = cpl_image_get_size_y(image);
    npix     = nx * ny;

    profile = cpl_image_collapse_create(image, 1);

    valid  = cpl_image_new(1, ny, CPL_TYPE_INT);
    ivalid = cpl_image_get_data_int(valid);
    data   = cpl_image_get_data_float(image);

    for (i = 0; i < npix; i++) {
        if (data[i] != 0.0) {
            ivalid[i/nx]++;
        }
    }

    cpl_image_divide(profile, valid);

/*
    cpl_image_save(profile, "profile.fits", 
                   CPL_BPP_IEEE_FLOAT, NULL, CPL_IO_DEFAULT);
    cpl_image_save(valid, "valid.fits", 
                   CPL_BPP_IEEE_FLOAT, NULL, CPL_IO_DEFAULT);
*/

    cpl_image_delete(valid);

    pdata = cpl_image_get_data(profile);

    row = 1;
    maxobjects = 0;
    totobj = 0;
    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
            continue;

        pos = position[i] + margin;
        count = length[i] - 2*margin;

        if (count < min_pixels)
            continue;

        p = pdata + pos;


        /*
         * Count peaks candidates
         */

        npeaks = 0;
        if (p[0] > p[1] && p[1] > p[2] && p[2] > p[3] && p[3] > 0) {
            npeaks++;
        }

        up = 0;
        for (j = 0; j < count - 3; j++) {
            if (p[j] > 0) {
                if (p[j+1] > p[j]) {
                    up++;
                }
                else {
                    if (up > 2) {
                        if (p[j+1] > p[j+2] && p[j+2] > 0) {
                            if (p[j] > noise) {
                                npeaks++;
                            }
                        }
                    }
                    else if (up > 1) {
                        if (p[j+1] > p[j+2] && p[j+2] > p[j+3] && p[j+3] > 0) {
                            if (p[j] > noise) {
                                npeaks++;
                            }
                        }
                    }
                    up = 0;
                }
            }
            else {
                up = 0;
            }
        }

        if (p[count-1] > p[count-2] && p[count-2] > p[count-3] 
            && p[count-3] > p[count-4] && p[count-4] > 0) {
            npeaks++;
        }

        if (npeaks == 0)
            continue;


        /*
         * Get candidates parameters
         */

        reject = cpl_calloc(npeaks, sizeof(int));
        bright = cpl_calloc(npeaks, sizeof(double));
        place  = cpl_calloc(npeaks, sizeof(double));

        npeaks = 0;
        if (p[0] > p[1] && p[1] > p[2] && p[2] > p[3] && p[3] > 0) {
            bright[0] = p[0];
            place[0] = position[i] + margin;
            npeaks++;
        }

        up = 0;
        for (j = 0; j < count - 3; j++) {
            if (p[j] > 0) {
                if (p[j+1] > p[j]) {
                    up++;
                }
                else {
                    if (up > 2) {
                        if (p[j+1] > p[j+2] && p[j+2] > 0) {
                            if (p[j] > noise) {
                                bright[npeaks] = p[j];
                                place[npeaks] = position[i] + margin + j + 1
                                       + values_to_dx(p[j-1], p[j], p[j+1]);
                                npeaks++;
                            }
                        }
                    }
                    else if (up > 1) {
                        if (p[j+1] > p[j+2] && p[j+2] > p[j+3] && p[j+3] > 0) {
                            if (p[j] > noise) {
                                bright[npeaks] = p[j];
                                place[npeaks] = position[i] + margin + j + 1
                                       + values_to_dx(p[j-1], p[j], p[j+1]);
                                npeaks++;
                            }
                        }
                    }
                    up = 0;
                }
            }
            else {
                up = 0;
            }
        }

        if (p[count-1] > p[count-2] && p[count-2] > p[count-3]
            && p[count-3] > p[count-4] && p[count-4] > 0) {
            bright[npeaks] = p[count-1];
            place[npeaks] = position[i] + count;
            npeaks++;
        }


        /*
         * Now select the uncontaminated peaks
         */

        if (fabs(place[0] - pos) < 1.0)
            reject[0] = 1;
        if (fabs(place[npeaks-1] - pos - count) < 1.0)
            reject[npeaks-1] = 1;
        for (j = 0; j < npeaks; j++) {
            for (k = 0; k < npeaks; k++) {
                if (k == j)
                    continue;
                mindistance = conradius * bright[k] / bright[j] 
                                        * bright[k] / bright[j];
                if (fabs(place[j] - place[k]) < mindistance)
                    reject[j] = 1;
            }
        }

/* new part */
        for (j = 0; j < npeaks; j++) {
            if (reject[j])
                continue;
            if (j) {
                low = (place[j-1]*bright[j] + place[j]*bright[j-1])
                    / (bright[j-1] + bright[j]) + 1;
            }
            else {
                low = pos;
            }
            if (j < npeaks - 1) {
                hig = (place[j+1]*bright[j] + place[j]*bright[j+1])
                    / (bright[j+1] + bright[j]) + 1;
            }
            else {
                hig = pos + count;
            }

            if (low < pos)
                low = pos;
            if (hig > pos + count)
                hig = pos + count;
            if (place[j] - low > maxradius)
                low = place[j] - maxradius;
            if (hig - place[j] > maxradius)
                hig = place[j] + maxradius;
            if (hig == low)
                reject[j] = 1;
        }
/* end new part */

        nobjects = npeaks;
        for (j = 0; j < npeaks; j++)
            if (reject[j])
                nobjects--;

        for (j = 0; j < nobjects; j++) {
            snprintf(name, MAX_COLNAME, "object_%d", j+1);
            if (cpl_table_has_column(slits, name))
                continue;
            cpl_table_new_column(slits, name, CPL_TYPE_DOUBLE);
            snprintf(name, MAX_COLNAME, "start_%d", j+1);
            cpl_table_new_column(slits, name, CPL_TYPE_INT);
            cpl_table_set_column_unit(slits, name, "pixel");
            snprintf(name, MAX_COLNAME, "end_%d", j+1);
            cpl_table_new_column(slits, name, CPL_TYPE_INT);
            cpl_table_set_column_unit(slits, name, "pixel");
            snprintf(name, MAX_COLNAME, "row_%d", j+1);
            cpl_table_new_column(slits, name, CPL_TYPE_INT);
            cpl_table_set_column_unit(slits, name, "pixel");
        }

        objpos = nobjects;
        for (j = 0; j < npeaks; j++) {
            if (reject[j])
                continue;
            if (j) {
                low = (place[j-1]*bright[j] + place[j]*bright[j-1])
                    / (bright[j-1] + bright[j]) + 1;
            }
            else {
               low = pos;
            }
            if (j < npeaks - 1) {
                hig = (place[j+1]*bright[j] + place[j]*bright[j+1])
                    / (bright[j+1] + bright[j]) + 1;
            }
            else {
                hig = pos + count;
            }

            if (low < pos)
                low = pos;
            if (hig > pos + count)
                hig = pos + count;
            if (place[j] - low > maxradius)
                low = place[j] - maxradius;
            if (hig - place[j] > maxradius)
                hig = place[j] + maxradius;

            snprintf(name, MAX_COLNAME, "object_%d", objpos);
            cpl_table_set_double(slits, name, i, place[j]);
            snprintf(name, MAX_COLNAME, "start_%d", objpos);
            cpl_table_set_int(slits, name, i, low);
            snprintf(name, MAX_COLNAME, "end_%d", objpos);
            cpl_table_set_int(slits, name, i, hig);
            snprintf(name, MAX_COLNAME, "row_%d", objpos);
            cpl_table_set_int(slits, name, i, row + objpos - 1);
            totobj++;
            objpos--;
        }

        row += nobjects;

        if (maxobjects < nobjects)
            maxobjects = nobjects;

        cpl_free(reject);
        cpl_free(bright);
        cpl_free(place);

    }

/*    nobjects = row - nobjects;     A bug, I think... */
    row = cpl_table_get_nrow(slits);

    for (i = 0; i < row; i++) {
        for (j = 0; j < maxobjects; j++) {
            snprintf(name, MAX_COLNAME, "row_%d", j+1);
            if (cpl_table_is_valid(slits, name, i))
                cpl_table_set_int(slits, name, i, totobj -
                                  cpl_table_get_int(slits, name, i, NULL));
        }
    }

    for (i = 0; i < maxobjects; i++) {
        snprintf(name, MAX_COLNAME, "start_%d", i+1);
        cpl_table_fill_invalid_int(slits, name, -1);
        snprintf(name, MAX_COLNAME, "end_%d", i+1);
        cpl_table_fill_invalid_int(slits, name, -1);
        snprintf(name, MAX_COLNAME, "row_%d", i+1);
        cpl_table_fill_invalid_int(slits, name, -1);
    }

    return profile;
}


/**
 * @brief
 *   Extract detected objects from rectified scientific frame
 *   
 * @param science     Rectified and sky subtracted scientific spectra
 * @param sky         Rectified sky spectra 
 * @param objects     Spectra and objects position table
 * @param extraction  Extraction mode: 0 = aperture, 1 = optimal (Horne)
 * @param ron         Read-out-noise in ADU
 * @param gain        Conversion from ADU to electrons (e-/ADU)
 * @param ncombined   Number of combined scientific frames
 * 
 * @return @c Images with extracted science, sky, and error spectra
 * 
 * The objects spatial extraction intervals are those listed in the
 * input @em objects table produced by the function @c mos_detect_objects().
 * The arguments @em ron, @em gain, and @em ncombined are used only in
 * case the @em extraction mode is set to 1 (optimal extraction). The
 * optimal extraction is based on Horne, K., (1986), PASP, 98, 609.
 * If the @em science and @em science_sky frames are the result of
 * the combination of different frames, the value of @em ron will be
 * divided by the square root of @em ncombined. 
 */

cpl_image **mos_extract_objects(cpl_image *science, cpl_image *science_var, 
                                cpl_image *sky,
                                cpl_table *objects, int extraction, double ron,
                                double gain, int ncombined)
{
    const char *func = "mos_extract_objects";

    char        name[MAX_COLNAME];

    cpl_image **output;
    cpl_image  *extracted;
    cpl_image  *extr_sky;
    cpl_image  *error;
    cpl_image  *sciwin;
    cpl_image  *sci_var_win = NULL;
    cpl_image  *skywin;
    int         nslits;
    int         nobjects;
    int         maxobjects;
    int         nx;
    int         ylow, yhig;
    int         i, j;


    if (science == NULL || sky == NULL) {
        cpl_msg_error(func, "Both scientific exposures are required in input");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (objects == NULL) {
        cpl_msg_error(func, "An object table is required in input");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (extraction < 0 || extraction > 1) {
        cpl_msg_error(func, "Invalid extraction mode (%d): it should be "
                      "either 0 or 1", extraction); 
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (ron < 0.0) {
        cpl_msg_error(func, "Invalid read-out-noise (%f ADU)", ron);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (gain < 0.1) {
        cpl_msg_error(func, "Invalid gain factor (%f e-/ADU)", gain);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (ncombined < 1) {
        cpl_msg_error(func, "Invalid number of combined frames (%d): "
                      "it should be at least 1", ncombined);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }


    /*
     * Count the max number of objects per slit. Note that maxobjects 
     * is intentionally the max number of objects increased by one.
     */

    maxobjects = 1;
    snprintf(name, MAX_COLNAME, "object_%d", maxobjects);
    while (cpl_table_has_column(objects, name)) {
        maxobjects++;
        snprintf(name, MAX_COLNAME, "object_%d", maxobjects);
    }


    /*
     * Count objects to extract
     */

    nobjects = 0;
    nslits = cpl_table_get_nrow(objects);

    for (i = 0; i < nslits; i++) {
        for (j = 1; j < maxobjects; j++) {
            snprintf(name, MAX_COLNAME, "object_%d", j);
            if (cpl_table_is_valid(objects, name, i))
                nobjects++;
        }
    }

    if (nobjects == 0)
        return NULL;

    nx = cpl_image_get_size_x(science);

    output = cpl_calloc(3, sizeof(cpl_image *));
    extracted = output[0] = cpl_image_new(nx, nobjects, CPL_TYPE_FLOAT);
    extr_sky  = output[1] = cpl_image_new(nx, nobjects, CPL_TYPE_FLOAT);
    error     = output[2] = cpl_image_new(nx, nobjects, CPL_TYPE_FLOAT);


    /*
     * Extract objects
     */

    nobjects = 0;
    for (i = 0; i < nslits; i++) {
        for (j = 1; j < maxobjects; j++) {
            snprintf(name, MAX_COLNAME, "object_%d", j);
            if (cpl_table_is_valid(objects, name, i)) {
                snprintf(name, MAX_COLNAME, "start_%d", j);
                ylow = cpl_table_get_int(objects, name, i, NULL);
                snprintf(name, MAX_COLNAME, "end_%d", j);
                yhig = cpl_table_get_int(objects, name, i, NULL);
                snprintf(name, MAX_COLNAME, "row_%d", j);
                nobjects = cpl_table_get_int(objects, name, i, NULL);
                sciwin = cpl_image_extract(science, 1, ylow+1, nx, yhig);
                if(science_var != NULL)
                    sci_var_win = cpl_image_extract(science_var, 1, ylow+1, nx, yhig);
                skywin = cpl_image_extract(sky, 1, ylow+1, nx, yhig);
/*
 * Cleaning the cosmics locally was really NOT a good idea...
 * I leave it here, commented out, to never forget this mistake!

                if (extraction) {
                    mos_clean_cosmics(sciwin, gain, -1., -1.);
                }
 */
                mos_extraction(sciwin, sci_var_win, skywin, extracted, extr_sky, error, 
                               nobjects, extraction, ron, gain, ncombined);

                /*
                 * Hidden check whether the spectrum was saturated or not
                 */

                {
                    cpl_image *total = cpl_image_add_create(sciwin, skywin);
                    float     *data  = cpl_image_get_data_float(total);
                    int        size  = cpl_image_get_size_x(total)
                                     * cpl_image_get_size_y(total);
                    int        k;
                    char      *saturation_level = getenv("SATURATION_LEVEL");
                    float      saturation = 62000.0;
                    char      *max_saturated = getenv("MAX_SATURATED");
                    int        max_satur = 10;
                    int        saturated;

                    if (saturation_level)
                        saturation = atof(saturation_level);

                    if (max_saturated)
                        max_satur = atoi(max_saturated);

                    saturated = 0;
                    for (k = 0; k < size; k++) {
                        if (data[k] > saturation) {
                            saturated++;
                            if (saturated > max_satur) {
                                break;
                            }
                        }
                    }

                    if (saturated > max_satur)
                        saturated = 1;
                    else
                        saturated = 0;

                    data = cpl_image_get_data(extracted);
                    data[nobjects * nx] = saturated;
                }

                cpl_image_delete(sciwin);
                cpl_image_delete(skywin);
                nobjects++;
            }
        }
    }

    return output;

}


/**
 * @brief
 *   Compute mean spectral resolution at a given arc lamp line
 * 
 * @param image       Rectified and wavelength calibrated arc lamp image
 * @param lambda      Wavelength to examine
 * @param startwave   Shortest wavelength in image.
 * @param dispersion  Wavelength units per image pixel
 * @param saturation  Saturation value
 * @param mfwhm       Returned: median FWHM
 * @param rmsfwhm     Returned: RMS of median FWHM
 * @param resolution  Returned: spectral resolution
 * @param rmsres      Returned: RMS of spectral resolution
 * @param nlines      Returned: Number of examined line profiles
 * 
 * @return 1 in case of success, 0 in case of failure
 * 
 * Given a wavelength, determine spectral resolution from a given arc
 * lamp line. A high S/N is assumed (virtually no noise). Working
 * with 2D extracted images!
 */

int mos_spectral_resolution(cpl_image *image, double lambda, double startwave, 
                            double dispersion, int saturation, 
                            double *mfwhm, double *rmsfwhm,
                            double *resolution, double *rmsres, int *nlines)
{
    cpl_vector *vector;

    int     i, j, n, m;
    int     position, maxpos;
    int     xlen, ylen;
    int     sp, ep;
    int     radius;
    int     sradius = 40;
    int     threshold = 250;    /* Peak must be so many ADUs above min */

    int     ifwhm;
    double  fwhm;
    double *buffer;
    double  min, max, halfmax;
    double  cut = 1.5;         /* To cut outliers from FWHM values (pixel) */
    double  value, rms;

    float  *data;


    *resolution = 0.0;
    *rmsres = 0.0;
    *nlines = 0;

    xlen = cpl_image_get_size_x(image);
    ylen = cpl_image_get_size_y(image);
    data = cpl_image_get_data(image);

    buffer = cpl_malloc(ylen * sizeof(double));

    /*
     *  Closest pixel to specified wavelength.
     */

    position = floor((lambda - startwave) / dispersion + 0.5);

    sp = position - sradius;
    ep = position + sradius;

    if (sp < 0 || ep > xlen) {
        cpl_free(buffer);
        return 0;
    }

    for (i = 0, n = 0; i < ylen; i++) {    /*  For each row of each slit  */

        /*
         *  Search interval for peak. Abort if too close to image border.
         */

        radius = mos_lines_width(data + i*xlen + position - sradius, 
                                 2*sradius + 1);
        if (radius < 5)
            radius = 5;

        sp = position - radius;
        ep = position + radius;

        if (sp < 0 || ep > xlen) {
            cpl_free(buffer);
            return 0;
        }


        /*
         *  Determine min-max value and position.
         */

        maxpos = sp;
        min = max = data[sp + i * xlen];
        for (j = sp; j < ep; j++) {
            if (data[j + i * xlen] > max) {
                max = data[j + i * xlen];
                maxpos = j;
            }
            if (data[j + i * xlen] < min) {
                min = data[j + i * xlen];
            }
        }

        if (fabs(min) < 0.0000001)        /* Truncated spectrum */
            continue;

        if (max - min < threshold)        /* Low signal... */
            continue;

        if (max > saturation)             /* Saturation */
            continue;

        /*
         *  Determine FWHM counting pixels with value greater than
         *  half of the max value, to the right and to the left of
         *  the max. Linear interpolation between the pixels where
         *  the transition happens.
         */

        halfmax = (max + min)/ 2.0;

        fwhm = 0.0;
        ifwhm = 0;
        for (j = maxpos; j < maxpos + radius; j++) {
            if (j < xlen) {
                if (data[j + i * xlen] < halfmax) {
                    fwhm = ifwhm + (data[j - 1 + i * xlen] - halfmax)
                         / (data[j - 1 + i * xlen] - data[j + i * xlen]);
                    break;
                }
                ifwhm++;
            }
        }

        ifwhm = 0;
        for (j = maxpos; j > maxpos - radius; j--) {
            if (j >= 0) {
                if (data[j + i * xlen] < halfmax) {
                    fwhm += ifwhm + (data[j + 1 + i * xlen] - halfmax)
                          / (data[j + 1 + i * xlen] - data[j + i * xlen]);
                    break;
                }
                ifwhm++;
            }
        }

        if (fwhm > 3.0) {
            buffer[n] = fwhm - 2.0;
            n++;
        }

    }

    if (n == 0) {
        cpl_free(buffer);
        return 0;
    }

    vector = cpl_vector_wrap(n, buffer);
    value = cpl_vector_get_median_const(vector);
    cpl_vector_unwrap(vector);

    rms = 0.0;
    for (i = 0, m = 0; i < n; i++) {
        if (fabs(buffer[i] - value) < cut) {
            rms += fabs(buffer[i] - value);
            m++;
        }
    }

    cpl_free(buffer);

    if (m < 3)
        return 0;

    rms /= m;
    rms *= 1.25;       /* Factor to convert average deviation to sigma */

    value *= dispersion;
    rms *= dispersion;

    *mfwhm = value;
    *rmsfwhm = rms;

    *resolution = lambda / value;
    *rmsres = *resolution * rms / value;

    *nlines = m;

    return 1;
}


/** 
 * @brief
 *   Compute mean spectral resolution at a given arc lamp line
 *  
 * @param image       Rectified and wavelength calibrated arc lamp image
 * @param startwave   Shortest wavelength in image.
 * @param dispersion  Wavelength units per image pixel
 * @param saturation  Saturation value
 * @param lines       Line catalog
 *
 * @return Spectral resolution table
 *
 * This function loops on the function @c mos_spectral_resolution() for
 * each line of the input line catalog, and creates a spectral resolution
 * table. This table has a column "wavelength" listing the line catalog 
 * wavelengths, a column "resolution" containing the corresponding mean
 * spectral resolutions, a column "rms" with the population standard
 * deviation of the measured resolutions, and a column "nlines" reporting
 * the number of measurement for each wavelength.
 */

cpl_table *mos_resolution_table(cpl_image *image, double startwave, 
                                double dispersion, int saturation, 
                                cpl_vector *lines)
{

    cpl_table *table;
    double    *line;
    double     fwhm;
    double     rmsfwhm;
    double     resolution;
    double     rmsres;
    int        nref;
    int        nlines;
    int        i;


    nref = cpl_vector_get_size(lines);
    line = cpl_vector_get_data(lines);

    table = cpl_table_new(nref);
    cpl_table_new_column(table, "wavelength", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(table, "wavelength", "Angstrom");
    cpl_table_new_column(table, "fwhm", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(table, "fwhm", "Angstrom");
    cpl_table_new_column(table, "fwhm_rms", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(table, "fwhm_rms", "Angstrom");
    cpl_table_new_column(table, "resolution", CPL_TYPE_DOUBLE);
    cpl_table_new_column(table, "resolution_rms", CPL_TYPE_DOUBLE);
    cpl_table_new_column(table, "nlines", CPL_TYPE_INT);

    for (i = 0; i < nref; i++) {
        if (mos_spectral_resolution(image, line[i], startwave, dispersion, 
                                    saturation, &fwhm, &rmsfwhm, 
                                    &resolution, &rmsres, &nlines)) {
            cpl_table_set_double(table, "wavelength", i, line[i]);
            cpl_table_set_double(table, "fwhm", i, fwhm);
            cpl_table_set_double(table, "fwhm_rms", i, rmsfwhm);
            cpl_table_set_double(table, "resolution", i, resolution);
            cpl_table_set_double(table, "resolution_rms", i, rmsres);
            cpl_table_set_int(table, "nlines", i, nlines);
        }
        else {
            cpl_table_set_int(table, "nlines", i, 0);
            cpl_table_set_double(table, "wavelength", i, line[i]);
        }
    }

    if (cpl_table_has_valid(table, "wavelength"))
        return table;

    cpl_table_delete(table);

    return NULL;
    
}


/**
 * @brief
 *   Integrate signal from wavelength and spatial interval
 *
 * @param image    CCD exposure
 * @param wavemap  Wavelengths map of CCD exposure
 * @param ystart   Start Y-pixel coordinate on CCD
 * @param yend     End Y-pixel coordinate on CCD
 * @param wstart   Start wavelength
 * @param wend     End wavelength
 *
 * @return Integrated signal
 *
 * This function sum the signal in the specified interval. @em ystart
 * is inclusive and @em wend exclusive.
 */

double mos_integrate_signal(cpl_image *image, cpl_image *wavemap,
                            int ystart, int yend, double wstart, double wend)
{
    const char *func = "mos_integrate_signal";

    double sum;
    float *sdata;
    float *wdata;
    int    nx, ny;
    int    x, y;
    

    if (image == NULL || wavemap == NULL) { 
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0.0;
    }

    if (ystart > yend || wstart >= wend) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return 0.0;
    }

    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);

    if (!(nx == cpl_image_get_size_x(wavemap) 
        && ny == cpl_image_get_size_y(wavemap))) {
        cpl_error_set(func, CPL_ERROR_INCOMPATIBLE_INPUT);
        return 0.0;
    }

    if (ystart < 0 || yend > ny) {
        cpl_error_set(func, CPL_ERROR_ACCESS_OUT_OF_RANGE);
        return 0.0;
    }

    sdata = cpl_image_get_data(image);
    wdata = cpl_image_get_data(wavemap);

    sdata += ystart*nx;
    wdata += ystart*nx;

    sum = 0.0;
    for (y = ystart; y < yend; y++) {
        for (x = 0; x < nx; x++) {
            if (wdata[x] < wstart || wdata[x] > wend)
                continue;
            sum += sdata[x];
        }
        sdata += nx;
        wdata += nx;
    }

    return sum;

}

/****************************************************************************
 * From this point on, the instrument dependent functions are added:
 * they are functions that retrieve information that is stored in
 * the data headers in some instrument specific way, such as the
 * location of overscans, the slits positions on the telescope
 * focal plane, the gain factor, etc.
 */


/**
 * @brief
 *   Create slit location table from FITS header of FORS2-MXU data
 *
 * @param header    FITS header of FORS data containing slits information
 *
 * @return A slit position table
 *
 * This function is meant to convert the information contained in 
 * FORS2 MXU data FITS header into a slit position table. This table 
 * will contain the double precision columns labeled @em xtop, 
 * @em ytop, @em xbottom, and @em ybottom, containing the start 
 * and end coordinates of the slits on the telescope focal plane
 * (mask). The coordinates are expected to have a horizontal 
 * dispersion direction and red wavelengths on the right side.
 * A flip to the Y coordinate is applied to match the increasing
 * Y CCD pixel coordinate. A slit identifying integer column, 
 * labeled "slit_id", containing unique slit identifiers, will
 * also be added. In case of FORS2 data only the slits pertaining
 * to the used chip are loaded.
 */

cpl_table *mos_load_slits_fors_mxu(cpl_propertylist *header)
{
    const char *func = "mos_load_slits_fors_mxu";

    cpl_table  *slits;
    char        keyname[MAX_COLNAME];
    const char *instrume;
    const char *target_name;
    float       slit_x;
    float       slit_y;
    float       length;
/*    double      arc2mm = 0.53316;         */
    double      arc2mm = 0.528;
    int         nslits;
    int         slit_id;
    int         fors;
    int         chip;
    int         found;

    /*
     * The limits below are used to exclude from the loaded slit list
     * any slit that surely doesn't belong to the used chip. This is
     * a way to reduce the chance of ambiguous slit identification.
     */

    float      low_limit1 = 10.0;
    float      hig_limit2 = 30.0;


    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        return NULL;
    }

    if (header == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }


    /*
     * See if this is FORS1 or FORS2;
     */

    instrume = cpl_propertylist_get_string(header, "INSTRUME");

    fors = 0;
    if (instrume[4] == '1')
        fors = 1;
    if (instrume[4] == '2')
        fors = 2;

    if (fors != 2) {
        cpl_msg_error(func, "Wrong instrument: %s\n"
                      "FORS2 is expected for MXU data", instrume);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }


    /*
     * The master and slave chips can be identified by their positions
     * in the chip array in the case of FORS2 data (with fors1 the chip
     * is always 1). chip = 2 is the master, chip = 1 is the slave.
     */

    chip = cpl_propertylist_get_int(header, "ESO DET CHIP1 Y");

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(func, "Missing keyword ESO DET CHIP1 Y "
                      "in FITS header");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (chip != 1 && chip != 2) {
        cpl_msg_error(func, "Unexpected chip position in keyword "
                      "ESO DET CHIP1 Y: %d", chip);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }


    /*
     * Count slits in header (excluding reference slits, and the slits
     * that _surely_ belong to the other chip)
     */

    nslits = 0;
    slit_id = 0;
    found = 1;

    while (found) {
        slit_id++;
        snprintf(keyname, MAX_COLNAME, "ESO INS MOS%d YPOS", slit_id + 100);
        if (cpl_propertylist_has(header, keyname)) {
            slit_y = cpl_propertylist_get_double(header, keyname);

            if (chip == 1)
                if (slit_y < low_limit1)
                    continue;
            if (chip == 2)
                if (slit_y > hig_limit2)
                    continue;
                
            snprintf(keyname, MAX_COLNAME, "ESO INS TARG%d NAME", 
                     slit_id + 100);
            if (cpl_propertylist_has(header, keyname)) {
                target_name = cpl_propertylist_get_string(header, keyname);
                if (strncmp(target_name, "refslit", 7))
                    nslits++;
            }
            else
                nslits++;
        }
        else
            found = 0;
    }

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(func, "%s while loading slits coordinates from "
                      "FITS header", cpl_error_get_message());
        cpl_error_set_where(func);
        return NULL;
    }

    if (nslits == 0)  {
        cpl_msg_error(func, "No slits coordinates found in header");
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    slits = cpl_table_new(nslits);
    cpl_table_new_column(slits, "slit_id", CPL_TYPE_INT);
    cpl_table_new_column(slits, "xtop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ytop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "xbottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ybottom", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(slits, "xtop",    "pixel");
    cpl_table_set_column_unit(slits, "ytop",    "pixel");
    cpl_table_set_column_unit(slits, "xbottom", "pixel");
    cpl_table_set_column_unit(slits, "ybottom", "pixel");

    nslits = 0;
    slit_id = 0; 
    found = 1;
    while (found) {
        slit_id++;
        snprintf(keyname, MAX_COLNAME, "ESO INS MOS%d YPOS", slit_id + 100);
        if (cpl_propertylist_has(header, keyname)) {
            slit_y = cpl_propertylist_get_double(header, keyname);

            if (chip == 1) 
                if (slit_y < low_limit1)
                    continue;
            if (chip == 2)
                if (slit_y > hig_limit2)
                    continue;

            /*
             * Y-flip the slit position, to match CCD pixel coordinate
             * convention
             */

            slit_y = -slit_y;

            snprintf(keyname, MAX_COLNAME, "ESO INS MOS%d XPOS", slit_id + 100);
            slit_x = cpl_propertylist_get_double(header, keyname);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_table_delete(slits);
                cpl_msg_error(func, "Missing keyword %s in FITS header", 
                              keyname);
                cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
                return NULL;
            }

            snprintf(keyname, MAX_COLNAME, "ESO INS MOS%d LEN", slit_id + 100);
            length = cpl_propertylist_get_double(header, keyname);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_table_delete(slits);
                cpl_msg_error(func, "Missing keyword %s in FITS header", 
                              keyname);
                cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
                return NULL;
            }

            length *= arc2mm;

            snprintf(keyname, MAX_COLNAME, "ESO INS TARG%d NAME", 
                     slit_id + 100);
            if (cpl_propertylist_has(header, keyname)) {
                target_name = cpl_propertylist_get_string(header, keyname);
                if (strncmp(target_name, "refslit", 7)) {
                    cpl_table_set_int(slits, "slit_id", nslits, slit_id);
                    cpl_table_set(slits, "xtop", nslits, slit_x);
                    cpl_table_set(slits, "ytop", nslits, slit_y + length/2);
                    cpl_table_set(slits, "xbottom", nslits, slit_x);
                    cpl_table_set(slits, "ybottom", nslits, slit_y - length/2);
                    nslits++;
                }
            }
            else {
                cpl_table_set_int(slits, "slit_id", nslits, slit_id);
                cpl_table_set(slits, "xtop", nslits, slit_x);
                cpl_table_set(slits, "ytop", nslits, slit_y + length/2);
                cpl_table_set(slits, "xbottom", nslits, slit_x);
                cpl_table_set(slits, "ybottom", nslits, slit_y - length/2);
                nslits++;
            }
        }
        else
            found = 0;
    }

    return slits;
}


/**
 * @brief
 *   Create slit location table from FITS header of FORS1/2 MOS data
 *
 * @param header    FITS header of FORS data containing slits information
 *
 * @return A slit position table
 *
 * This function is meant to convert the information contained in 
 * FORS1 and FORS2 MOS data FITS header into a slit position table. 
 * This table will contain the double precision columns labeled 
 * @em xtop, @em ytop, @em xbottom, and @em ybottom, containing the 
 * start and end coordinates of the slits on the telescope focal 
 * plane (mask). The coordinates are expected to have a horizontal 
 * dispersion direction and red wavelengths on the right side.
 * A flip to the Y coordinate is applied to match the increasing
 * Y CCD pixel coordinate. A slit identifying integer column, 
 * labeled "slit_id", containing unique slit identifiers, will 
 * also be added. In case of FORS2 data only the slits pertaining
 * to the used chip are loaded.
 */

cpl_table *mos_load_slits_fors_mos(cpl_propertylist *header)
{
    const char *func = "mos_load_slits_fors_mos";

    cpl_table  *slits;
    char        keyname[MAX_COLNAME];
    const char *instrume;
    const char *chipname;
    float       slit_x;
    int         first_slit, last_slit;
    cpl_size    nslits;
    int         slit_id;
    int         fors;
    int         chip;
    int         fors_is_old;

    /*
     * The Y coordinates of the slits are fixed
     */

    float       ytop[19]    = { 113.9, 101.3,  89.9,  77.3,  65.9,  53.3, 
                                 41.9,  29.3,  17.9,   5.3,  -6.1, -18.7, 
                                -30.1, -42.7, -54.1, -66.7, -78.1, -90.7, 
                               -102.1 };
    float       ybottom[19] = { 102.1,  90.7,  78.1,  66.7,  54.1,  42.7,
                                 30.1,  18.7,   6.1,  -5.3, -17.9, -29.3,
                                -41.9, -53.3, -65.9, -77.3, -89.9, -101.3,
                               -113.9 };


    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        return NULL;
    }

    if (header == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }


    /*
     * See if this is FORS1 or FORS2;
     */

    instrume = cpl_propertylist_get_string(header, "INSTRUME");

    fors = 0;
    if (instrume[4] == '1')
        fors = 1;
    if (instrume[4] == '2')
        fors = 2;

    if (fors == 0) {
        cpl_msg_error(func, "Wrong instrument found in FITS header: %s", 
                      instrume);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    /* FIXME:
     * This is the way FORS1 data belong to the upgraded chips,
     * named "Marlene" and "Norma III". It's a quick solution,
     * there are hardcoded values here!!!
     */

    chipname = cpl_propertylist_get_string(header, "ESO DET CHIP1 ID");

    if (chipname[0] == 'M' || chipname[0] == 'N')
        fors_is_old = 0;
    else
        fors_is_old = 1;

    if (fors == 1 && fors_is_old) {
        first_slit = 1;
        last_slit = 19;
    }
    else {

        /*
         * The master and slave chips can be identified by their positions
         * in the chip array in the case of FORS2 data: chip = 2 is the 
         * master, chip = 1 is the slave.
         */

        chip = cpl_propertylist_get_int(header, "ESO DET CHIP1 Y");

        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_msg_error(func, "Missing keyword ESO DET CHIP1 Y "
                          "in FITS header");
            cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
            return NULL;
        }

        if (chip != 1 && chip != 2) {
            cpl_msg_error(func, "Unexpected chip position in keyword "
                          "ESO DET CHIP1 Y: %d", chip);
            cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
            return NULL;
        }

        if (chip == 1) {
            first_slit = 12;
            last_slit = 19;
        }
        else {
            first_slit = 1;
            last_slit = 11;
        }
    }


    /*
     * Count slits in header (excluding closed slits - i.e. those with
     * offsets greater than 115 mm - and the slits that do not belong 
     * to this chip)
     */

    nslits = 0;
    slit_id = 0;
    for (slit_id = first_slit; slit_id <= last_slit; slit_id++) {
        snprintf(keyname, MAX_COLNAME, "ESO INS MOS%d POS", slit_id);
        if (cpl_propertylist_has(header, keyname)) {
            slit_x = cpl_propertylist_get_double(header, keyname);
            if (fabs(slit_x) < 115.0)
                nslits++;
        }
        else {
            cpl_msg_error(func, "Missing keyword %s in FITS header", keyname);
            cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
            return NULL;
        }
    }

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(func, "%s while loading slits coordinates from "
                      "FITS header", cpl_error_get_message());
        cpl_error_set_where(func);
        return NULL;
    }

    if (nslits == 0)  {
        cpl_msg_error(func, "No slits coordinates found in header");
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    slits = cpl_table_new(nslits);
    cpl_table_new_column(slits, "slit_id", CPL_TYPE_INT);
    cpl_table_new_column(slits, "xtop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ytop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "xbottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ybottom", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(slits, "xtop",    "pixel");
    cpl_table_set_column_unit(slits, "ytop",    "pixel");
    cpl_table_set_column_unit(slits, "xbottom", "pixel");
    cpl_table_set_column_unit(slits, "ybottom", "pixel");

    nslits = 0;
    slit_id = 0;
    for (slit_id = first_slit; slit_id <= last_slit; slit_id++) {
        snprintf(keyname, MAX_COLNAME, "ESO INS MOS%d POS", slit_id);
        slit_x = cpl_propertylist_get_double(header, keyname);
        if (fabs(slit_x) < 115.0) {
            cpl_table_set_int(slits, "slit_id", nslits, slit_id);
            cpl_table_set(slits, "xtop", nslits, slit_x);
            cpl_table_set(slits, "ytop", nslits, ytop[slit_id-1]);
            cpl_table_set(slits, "xbottom", nslits, slit_x);
            cpl_table_set(slits, "ybottom", nslits, ybottom[slit_id-1]);
            nslits++;
        }
    }

    return slits;
}


/**
 * @brief
 *   Create slit location table from FITS header of FORS1/2 LSS data
 *
 * @param header    FITS header of FORS data containing slits information
 *
 * @return A slit position table
 *
 * This function is meant to convert the information contained in 
 * FORS1 and FORS2 LSS data FITS header into a slit position table. 
 * This table will contain the double precision columns labeled 
 * @em xtop, @em ytop, @em xbottom, and @em ybottom, containing the 
 * start and end coordinates of the slits on the telescope focal 
 * plane (mask). The coordinates are expected to have a horizontal 
 * dispersion direction and red wavelengths on the right side.
 * A flip to the Y coordinate is applied to match the increasing
 * Y CCD pixel coordinate. A slit identifying integer column, 
 * labeled "slit_id", containing unique slit identifiers, will 
 * also be added. This identifier will be set to 1 for lSlit0_3arcsec,
 * 2 for lSlit0_4arcsec, up to 9 for lSlit2_5arcsec. In case of FORS2 
 * data only the part of the slit pertaining to the used chip is loaded.
 */

cpl_table *mos_load_slits_fors_lss(cpl_propertylist *header)
{
    const char *func = "mos_load_slits_fors_lss";

    cpl_table  *slits;
    char       *slit_name;
    const char *instrume;
    int         fors;
    int         chip;
    float       ytop;
    float       ybottom;

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        return NULL;
    }

    if (header == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }


    /*
     * See if this is FORS1 or FORS2;
     */

    instrume = cpl_propertylist_get_string(header, "INSTRUME");

    fors = 0;
    if (instrume[4] == '1')
        fors = 1;
    if (instrume[4] == '2')
        fors = 2;

    if (fors == 0) {
        cpl_msg_error(func, "Wrong instrument found in FITS header: %s", 
                      instrume);
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (fors == 1) {
        ytop = 109.94;
        ybottom = -109.94;
    }
    else {

        /*
         * The master and slave chips can be identified by their positions
         * in the chip array in the case of FORS2 data: chip = 2 is the 
         * master, chip = 1 is the slave.
         */

        chip = cpl_propertylist_get_int(header, "ESO DET CHIP1 Y");

        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_msg_error(func, "Missing keyword ESO DET CHIP1 Y "
                          "in FITS header");
            cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
            return NULL;
        }

        if (chip != 1 && chip != 2) {
            cpl_msg_error(func, "Unexpected chip position in keyword "
                          "ESO DET CHIP1 Y: %d", chip);
            cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
            return NULL;
        }

        if (chip == 1) {
            ytop = 30.0;
            ybottom = -109.94;
        }
        else {
            ytop = 109.94;
            ybottom = -20.0;
        }
    }


    slits = cpl_table_new(1);
    cpl_table_new_column(slits, "slit_id", CPL_TYPE_INT);
    cpl_table_new_column(slits, "xtop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ytop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "xbottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ybottom", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(slits, "xtop",    "pixel");
    cpl_table_set_column_unit(slits, "ytop",    "pixel");
    cpl_table_set_column_unit(slits, "xbottom", "pixel");
    cpl_table_set_column_unit(slits, "ybottom", "pixel");

    slit_name = (char *)cpl_propertylist_get_string(header, 
                                                    "ESO INS SLIT NAME");

    cpl_table_set(slits, "ytop", 0, ytop);
    cpl_table_set(slits, "ybottom", 0, ybottom);

    if (!strncmp(slit_name, "lSlit0_3arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 1);
        cpl_table_set(slits, "xbottom", 0, -0.075);
        cpl_table_set(slits, "xtop", 0, 0.075);
    }
    else if (!strncmp(slit_name, "lSlit0_4arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 2);
        cpl_table_set(slits, "xbottom", 0, 5.895);
        cpl_table_set(slits, "xtop", 0, 6.105);
    }
    else if (!strncmp(slit_name, "lSlit0_5arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 3);
        cpl_table_set(slits, "xbottom", 0, -6.135);
        cpl_table_set(slits, "xtop", 0, -5.865);
    }
    else if (!strncmp(slit_name, "lSlit0_7arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 4);
        cpl_table_set(slits, "xbottom", 0, 11.815);
        cpl_table_set(slits, "xtop", 0, 12.185);
    }
    else if (!strncmp(slit_name, "lSlit1_0arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 5);
        cpl_table_set(slits, "xbottom", 0, -12.265);
        cpl_table_set(slits, "xtop", 0, -11.735);
    }
    else if (!strncmp(slit_name, "lSlit1_3arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 6);
        cpl_table_set(slits, "xbottom", 0, 17.655);
        cpl_table_set(slits, "xtop", 0, 18.345);
    }
    else if (!strncmp(slit_name, "lSlit1_6arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 7);
        cpl_table_set(slits, "xbottom", 0, -18.425);
        cpl_table_set(slits, "xtop", 0, -17.575);
    }
    else if (!strncmp(slit_name, "lSlit2_0arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 8);
        cpl_table_set(slits, "xbottom", 0, 23.475);
        cpl_table_set(slits, "xtop", 0, 24.525);
    }
    else if (!strncmp(slit_name, "lSlit2_5arcsec", 14)) {
        cpl_table_set_int(slits, "slit_id", 0, 9);
        cpl_table_set(slits, "xbottom", 0, -24.66);
        cpl_table_set(slits, "xtop", 0, -23.34);
    }
    else {
        cpl_msg_error(func, "Invalid slit %s in keyword ESO INS SLIT NAME",
                      slit_name);
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        cpl_table_delete(slits);
        return NULL;
    }

    return slits;
}


/**
 * @brief
 *   Return gain factor for a VIMOS exposure
 *   
 * @param header  FITS header of VIMOS data containing information
 * 
 * @return Gain factor (e-/ADU)
 * 
 * This function is meant to read the gain factor from VIMOS data 
 * FITS headers. In the specific case of VIMOS the keyword used is 
 * ESO DET OUT1 CONAD. If no keyword is found a negative gain factor
 * is returned.
 */

double mos_get_gain_vimos(cpl_propertylist *header)
{
    const char *func = "mos_get_gain_vimos";

    double gain = -1.0;


    if (cpl_error_get_code() != CPL_ERROR_NONE)
        return gain;

    if (header == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return gain;
    }

    gain = cpl_propertylist_get_double(header, "ESO DET OUT1 CONAD");
    if (cpl_error_get_code()) {
        cpl_error_set_where(func);
        gain = -1.0;
    }

    return gain;

}


/**
 * @brief
 *   Create slit location table from FITS header of VIMOS data
 *
 * @param header    FITS header of VIMOS data containing slits information
 * @param include_ref Include reference slits
 *
 * @return A slit position table
 *
 * This function is meant to convert the information contained in 
 * VIMOS data FITS header into a slit position table. This table 
 * will contain the double precision columns labeled @em xtop, 
 * @em ytop, @em xbottom, and @em ybottom, containing the start 
 * and end coordinates of the slits on the telescope focal plane
 * (mask). The coordinates will be rotated to have a horizontal 
 * dispersion direction and red wavelengths on the right side.
 * A slit identifying integer column, labeled "slit_id", containing 
 * unique slit identifiers, will also be added.
 */

cpl_table *mos_load_slits_vimos(cpl_propertylist *header, int include_ref)
{
    const char *func = "mos_load_slits_vimos";

    cpl_table *slits;
    char       keyname[MAX_COLNAME];
    float      tolerance = 1.0;
    float      slit_x;
    float      slit_y;
    float      slit_x_next;
    float      slit_y_next;
    float      dim_x;
    float      dim_y;
    int        nslits;
    int        nrefslits;
    int        duplicated;
    int        refduplicated;
    int        slit_id;
    int        slit_id_next;
    int        curved;
    int        i, j, k;


    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        return NULL;
    }

    if (header == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    nslits = cpl_propertylist_get_int(header, "ESO INS SLIT NO");

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_error_set_where(func);
        return NULL;
    }

    /*
     * Check there are no duplicated slits
     */

    duplicated = 0;
    for (i = 0; i < nslits; i++) {
        sprintf(keyname, "ESO INS SLIT%d X", i+1);
        slit_x = cpl_propertylist_get_double(header, keyname);
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_error_set_where(func);
            return NULL;
        }
        sprintf(keyname, "ESO INS SLIT%d Y", i+1);
        slit_y = cpl_propertylist_get_double(header, keyname);
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_error_set_where(func);
            return NULL;
        }
        for (j = i + 1; j < nslits; j++) {
            sprintf(keyname, "ESO INS SLIT%d X", j+1);
            slit_x_next = cpl_propertylist_get_double(header, keyname);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_error_set_where(func);
                return NULL;
            }
            sprintf(keyname, "ESO INS SLIT%d Y", j+1);
            slit_y_next = cpl_propertylist_get_double(header, keyname);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_error_set_where(func);
                return NULL;
            }

            if ((slit_x - slit_x_next) * (slit_x - slit_x_next) +
                (slit_y - slit_y_next) * (slit_y - slit_y_next) < tolerance) {
                sprintf(keyname, "ESO INS SLIT%d ID", i+1);
                slit_id = cpl_propertylist_get_int(header, keyname);
                if (cpl_error_get_code() != CPL_ERROR_NONE) {
                    cpl_error_set_where(func);
                    return NULL;
                }
                sprintf(keyname, "ESO INS SLIT%d ID", j+1);
                slit_id_next = cpl_propertylist_get_int(header, keyname);
                if (cpl_error_get_code() != CPL_ERROR_NONE) {
                    cpl_error_set_where(func);
                    return NULL;
                }
                cpl_msg_warning(func, "Slit %d (ID=%d) is a duplicate of "
                                "slit %d (ID=%d), and it is therefore ignored", 
                                j+1, slit_id_next, i+1, slit_id);
                ++duplicated;
                break; //We want only the first slit that is duplicated.
                       //If both slits b and c are duplications of a, 
                       //first b will be assigned as a duplication of a,
                       //this break will be executed and the outer loop
                       //will eventually reach b, and then c will be catched
                       //as a duplicate of c. In total: two duplications.
            }
        }
    }
    
    nrefslits = 0;
    refduplicated = 0;
    if(include_ref && cpl_propertylist_has(header, "ESO INS REF NO"))
    {
        nrefslits = cpl_propertylist_get_int(header, "ESO INS REF NO");
        for (i = 0; i < nrefslits; i++) {
            sprintf(keyname, "ESO INS REF%d X", i+1);
            slit_x = cpl_propertylist_get_double(header, keyname);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_error_set_where(func);
                return NULL;
            }
            sprintf(keyname, "ESO INS REF%d Y", i+1);
            slit_y = cpl_propertylist_get_double(header, keyname);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_error_set_where(func);
                return NULL;
            }
            for (j = i + 1; j < nrefslits; j++) {
                sprintf(keyname, "ESO INS REF%d X", j+1);
                slit_x_next = cpl_propertylist_get_double(header, keyname);
                if (cpl_error_get_code() != CPL_ERROR_NONE) {
                    cpl_error_set_where(func);
                    return NULL;
                }
                sprintf(keyname, "ESO INS REF%d Y", j+1);
                slit_y_next = cpl_propertylist_get_double(header, keyname);
                if (cpl_error_get_code() != CPL_ERROR_NONE) {
                    cpl_error_set_where(func);
                    return NULL;
                }

                if ((slit_x - slit_x_next) * (slit_x - slit_x_next) +
                    (slit_y - slit_y_next) * (slit_y - slit_y_next) < tolerance) {
                    sprintf(keyname, "ESO INS REF%d ID", i+1);
                    slit_id = cpl_propertylist_get_int(header, keyname);
                    if (cpl_error_get_code() != CPL_ERROR_NONE) {
                        cpl_error_set_where(func);
                        return NULL;
                    }
                    sprintf(keyname, "ESO INS REF%d ID", j+1);
                    slit_id_next = cpl_propertylist_get_int(header, keyname);
                    if (cpl_error_get_code() != CPL_ERROR_NONE) {
                        cpl_error_set_where(func);
                        return NULL;
                    }
                    cpl_msg_warning(func, "Reference slit %d (ID=%d) is a duplicate of "
                                    "reference slit %d (ID=%d), and it is therefore ignored", 
                                    j+1, slit_id_next, i+1, slit_id);
                    ++refduplicated;
                    break; //We want only the first slit that is duplicated.
                           //If both slits b and c are duplications of a, 
                           //first b will be assigned as a duplication of a,
                           //this break will be executed and the outer loop
                           //will eventually reach b, and then c will be catched
                           //as a duplicate of c. In total: two duplications.
                }
            }
        }
    }

    /*
     * Create table with valid slits
     */

    slits = cpl_table_new(nslits + nrefslits - duplicated - refduplicated);
    cpl_table_new_column(slits, "slit_id", CPL_TYPE_INT);
    cpl_table_new_column(slits, "xtop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ytop",    CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "xbottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ybottom", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "xwidth", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "ywidth", CPL_TYPE_DOUBLE);
    cpl_table_new_column(slits, "curved", CPL_TYPE_INT);
    cpl_table_set_column_unit(slits, "xtop",    "pixel");
    cpl_table_set_column_unit(slits, "ytop",    "pixel");
    cpl_table_set_column_unit(slits, "xbottom", "pixel");
    cpl_table_set_column_unit(slits, "ybottom", "pixel");
    cpl_table_set_column_unit(slits, "xwidth", "mm");
    cpl_table_set_column_unit(slits, "ywidth", "mm");

    for (k = 0, i = 0; i < nslits; i++) {
        sprintf(keyname, "ESO INS SLIT%d X", i+1);
        slit_x = cpl_propertylist_get_double(header, keyname);
        sprintf(keyname, "ESO INS SLIT%d Y", i+1);
        slit_y = cpl_propertylist_get_double(header, keyname);
        duplicated = 0;
        for (j = i + 1; j < nslits; j++) {
            sprintf(keyname, "ESO INS SLIT%d X", j+1);
            slit_x_next = cpl_propertylist_get_double(header, keyname);
            sprintf(keyname, "ESO INS SLIT%d Y", j+1);
            slit_y_next = cpl_propertylist_get_double(header, keyname);

            if ((slit_x - slit_x_next) * (slit_x - slit_x_next) +
                (slit_y - slit_y_next) * (slit_y - slit_y_next) < tolerance) {
                duplicated = 1; //If there are several duplicates is ok,
                                //with the first one we already have the
                                //condition.
            }
        }

        if (duplicated)
            continue;

        sprintf(keyname, "ESO INS SLIT%d ID", i+1);
        slit_id = cpl_propertylist_get_int(header, keyname);

        sprintf(keyname, "ESO INS SLIT%d DIMX", i+1);
        dim_x = cpl_propertylist_get_double(header, keyname);
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_error_set_where(func);
            return NULL;
        }

        sprintf(keyname, "ESO INS SLIT%d BEZIER DY", i+1);
        if (cpl_propertylist_has(header, keyname)) {
            curved = 1;
        }
        else {
            sprintf(keyname, "ESO INS SLIT%d DIMY", i+1);
            curved = 0;
        }
        dim_y = cpl_propertylist_get_double(header, keyname);
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_error_set_where(func);
            return NULL;
        }

        cpl_table_set_int(slits, "slit_id", k, slit_id);
        cpl_table_set(slits, "xtop", k, slit_x - dim_x/2);
        cpl_table_set(slits, "ytop", k, slit_y);
        cpl_table_set(slits, "xbottom", k, slit_x + dim_x/2);
        cpl_table_set(slits, "ybottom", k, slit_y);
        cpl_table_set(slits, "xwidth", k, dim_x);
        cpl_table_set(slits, "ywidth", k, dim_y);
        cpl_table_set_int(slits, "curved", k, curved);
        k++;
    }

    if(include_ref && cpl_propertylist_has(header, "ESO INS REF NO"))
    {
        for (i = 0; i < nrefslits; i++) {
            sprintf(keyname, "ESO INS REF%d X", i+1);
            slit_x = cpl_propertylist_get_double(header, keyname);
            sprintf(keyname, "ESO INS REF%d Y", i+1);
            slit_y = cpl_propertylist_get_double(header, keyname);
            refduplicated = 0;
            for (j = i + 1; j < nrefslits; j++) {
                sprintf(keyname, "ESO INS REF%d X", j+1);
                slit_x_next = cpl_propertylist_get_double(header, keyname);
                sprintf(keyname, "ESO INS REF%d Y", j+1);
                slit_y_next = cpl_propertylist_get_double(header, keyname);

                if ((slit_x - slit_x_next) * (slit_x - slit_x_next) +
                        (slit_y - slit_y_next) * (slit_y - slit_y_next) < tolerance) {
                    refduplicated = 1; //If there are several duplicates is ok,
                    //with the first one we already have the
                    //condition.
                }
            }

            if (refduplicated)
                continue;

            sprintf(keyname, "ESO INS REF%d ID", i+1);
            slit_id = cpl_propertylist_get_int(header, keyname);

            sprintf(keyname, "ESO INS REF%d DIMX", i+1);
            dim_x = cpl_propertylist_get_double(header, keyname);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_error_set_where(func);
                return NULL;
            }

            sprintf(keyname, "ESO INS REF%d DIMY", i+1);
            curved = 0; //There are no reference curved slits (Burkhard dixit)

            dim_y = cpl_propertylist_get_double(header, keyname);
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_error_set_where(func);
                return NULL;
            }

            cpl_table_set_int(slits, "slit_id", k, slit_id);
            cpl_table_set(slits, "xtop", k, slit_x - dim_x/2);
            cpl_table_set(slits, "ytop", k, slit_y);
            cpl_table_set(slits, "xbottom", k, slit_x + dim_x/2);
            cpl_table_set(slits, "ybottom", k, slit_y);
            cpl_table_set(slits, "xwidth", k, dim_x);
            cpl_table_set(slits, "ywidth", k, dim_y);
            cpl_table_set_int(slits, "curved", k, curved);
            k++;
        }
    }
    return slits;
}

/**
 * @brief
 *   Determining whether a VIMOS mask has spectral multplexing or not.
 *
 * @param slits   Slit positions produced by mos_load_slits_vimos()
 *
 * @return Max observed spectral multiplexing.
 */

int mos_check_multiplex(cpl_table *slits)
{
    cpl_propertylist *sort;
    int               nrow = cpl_table_get_nrow(slits);
    int               i, j, multiplex;
    double            xtop, xtopnext, xbot;
    double            tolerance = 1.0; // About spatially aligned slits (mm)


    /*
     * Sort according to increasing xtop positions, and when those are
     * equal according to the increasing ytop position.
     */

    sort = cpl_propertylist_new();
    cpl_propertylist_append_bool(sort, "xtop", 0);
    cpl_propertylist_append_bool(sort, "ytop", 0);
    cpl_table_sort(slits, sort);
    cpl_propertylist_delete(sort);

    /*
     * Now assign to each slit its multiplex order.
     */

    if (!cpl_table_has_column(slits, "multiplex")) {
        cpl_table_new_column(slits, "multiplex", CPL_TYPE_INT);
        cpl_table_fill_column_window_int(slits, "multiplex", 0, nrow, 0);
    }

    for (i = 0; i < nrow; i++) {
        multiplex = cpl_table_get_int(slits, "multiplex", i, NULL);
        xtop = cpl_table_get_double(slits, "xtop", i, NULL);
        xbot = cpl_table_get_double(slits, "xbottom", i, NULL);
        for (j = i + 1; j < nrow; j++) {
            xtopnext = cpl_table_get_double(slits, "xtop", j, NULL);
            if (xbot - xtopnext < tolerance)
                break;
            cpl_table_set_int(slits, "multiplex", j, multiplex + 1);
        }
    }

// debug
//    cpl_table_save(slits, NULL, NULL, "multiplex.fits", CPL_IO_DEFAULT);

    return 1 + cpl_table_get_column_max(slits, "multiplex");

}


/**
 * @brief
 *   Determining whether a VIMOS mask has spectral multplexing or not.
 *
 * @param slits   Slit positions produced by mos_load_slits_vimos()
 *
 * @return Max observed spectral multiplexing.
 */

int mos_check_multiplex_old(cpl_table *slits)
{
    cpl_propertylist *sort;
    int               nrow;
    int               i, multiplex, xprev, xcur;
    double            prev, cur;
    double            tolerance = 1.0; // About spatially aligned slits (mm)


    /*
     * Create an auxiliary column containing a sort of integer
     * x coordinate of the slit, to guarantee that slits at the
     * same spatial offset are recognised immediately as in spectral 
     * multiplexing.
     */

    sort = cpl_propertylist_new();
    cpl_propertylist_append_bool(sort, "xtop", 0);
    cpl_table_sort(slits, sort);
    cpl_propertylist_delete(sort);

    prev = cpl_table_get_double(slits, "xtop", 0, NULL);
    cpl_table_new_column(slits, "xind", CPL_TYPE_INT);
    cpl_table_set_int(slits, "xind", 0, prev);   // cast to int is intentional
    nrow = cpl_table_get_nrow(slits);
    for (i = 1; i < nrow; i++) {
        cur = cpl_table_get_double(slits, "xtop", i, NULL);
        if (fabs(prev - cur) > tolerance)
            prev = cur;
        cpl_table_set_int(slits, "xind", i, prev);
    }

    /*
     * Now sort according to increasing (integer) x positions, and when
     * those are equal (multiplexed) according to the increasing y position.
     */

    sort = cpl_propertylist_new();
    cpl_propertylist_append_bool(sort, "xind", 0);
    cpl_propertylist_append_bool(sort, "ytop", 0);
    cpl_table_sort(slits, sort);
    cpl_propertylist_delete(sort);

    /*
     * Now assign to each slit its multiplex order.
     */

    multiplex = 0;
    if (!cpl_table_has_column(slits, "multiplex"))
        cpl_table_new_column(slits, "multiplex", CPL_TYPE_INT);
    xprev = cpl_table_get_int(slits, "xind", 0, NULL);
    cpl_table_set_int(slits, "multiplex", 0, multiplex);
    nrow = cpl_table_get_nrow(slits);
    for (i = 1; i < nrow; i++) {
        xcur = cpl_table_get_int(slits, "xind", i, NULL);
        if (xcur == xprev) {
            multiplex++;
        }
        else {
            xprev = xcur;
            multiplex = 0;
        }
        cpl_table_set_int(slits, "multiplex", i, multiplex);
    }

// debug
//    cpl_table_save(slits, NULL, NULL, "multiplex.fits", CPL_IO_DEFAULT);

    cpl_table_erase_column(slits, "xind");

    return 1 + cpl_table_get_column_max(slits, "multiplex");

}


/**
 * @brief
 *   Assign to each slit a group index, where a group has no multiplexing
 *
 * @param slits   Slit positions produced by mos_load_slits_vimos()
 *
 * @return Number of groups
 */

int mos_assign_multiplex_group(cpl_table *slits)
{
    int    group = 0;
    int    nsel;
    double jump;
    double ymin, ymax, ykeep, y;
    double tolerance = 2.0;  // millimeters
    int    i, nrows;


    ymin = cpl_table_get_column_min(slits, "ytop") - 2*tolerance;
    ymax = cpl_table_get_column_max(slits, "ytop") + 2*tolerance;

    cpl_table_new_column(slits, "group", CPL_TYPE_INT);

    while (ymax - ymin > tolerance) {
        y     = ymax;
        ykeep = ymin;
        jump  = ymin - ymax;             // Initially negative on purpose
//printf("Start search between %.2f and %.2f\n", ymin, ymax);
        while (abs(jump) > tolerance) {
            jump /= 2;
            y += jump;
            if (jump > 0.0 && ymax - y <= tolerance) {
//printf("    It HAPPENED, y from %.2f to %.2f\n", y, ymax);
                y = ymax;
            }
//printf("    Search between %.2f and %.2f\n", ymin, y);
            cpl_table_select_all(slits);
            cpl_table_and_selected_double(slits, "ytop",
                                          CPL_NOT_LESS_THAN, ymin);
            cpl_table_and_selected_double(slits, "ytop", CPL_LESS_THAN, y);
            nsel = cpl_table_and_selected_invalid(slits, "group");
            if (nsel > 0) {
                cpl_table *subslits  = cpl_table_extract_selected(slits);
                int        multiplex;

                if (cpl_table_has_column(subslits, "multiplex"))
                    cpl_table_erase_column(subslits, "multiplex");

                multiplex = mos_check_multiplex(subslits);

                if (multiplex > 1) {
//printf("        MULTIPLEXED\n");
                    jump  = -fabs(jump);
                }
                else {
//printf("        Not multiplexed\n");
                    ykeep = y;
                    jump  = fabs(jump);
                }

                cpl_table_delete(subslits);
            }
            else {
//printf("        No slits found\n");
                jump  = fabs(jump);
            }
        }

        if (ykeep <= ymin) {
//printf("Leave loop, with group = %d\n", group);
//printf("ykeep <= ymin, %.2f <= %.2f\n", ykeep, ymin);
            break;
        }
//printf("-----------\n");

        cpl_table_select_all(slits);
        cpl_table_and_selected_double(slits, "ytop", CPL_NOT_LESS_THAN, ymin);
        cpl_table_and_selected_double(slits, "ytop", CPL_LESS_THAN, ykeep);
        nsel = cpl_table_and_selected_invalid(slits, "group");

        if (nsel <= 0) {
//printf("Leave loop, with no slits found\n");
            break;
        }

        nrows = cpl_table_get_nrow(slits);

        for (i = 0; i < nrows; i++) {
            if (cpl_table_is_selected(slits, i)) {
//printf("Write group %d to row %d\n", group, i);
                cpl_table_set_int(slits, "group", i, group);
            }
        }

        ymin = ykeep;

        ++group;
    }

    cpl_table_select_all(slits);

    if (cpl_table_has_invalid(slits, "group")) {
printf("SOMETHING'S WRONG\n");
cpl_table_dump_structure(slits, NULL);
cpl_table_dump(slits, 0, nrows, NULL);
        return 0;
    }

    return group;
}


/**
 * @brief
 *   Get the overscan positions from FITS header of VIMOS data
 *
 * @param header              FITS header containing overscan information
 * @param check_consistency   If true (non-zero), this function fails if
 *                            the sum of prescan, overscan and detector size
 *                            is different from the actual FITS image size.
 *                            (this does not hold true for old FORS data)
 *
 * @return A overscan position table
 *
 * This function is meant to convert the information contained in
 * VIMOS data FITS header into an overscan position table. This table
 * will contain the integer columns labeled @em xlow, @em ylow, @em xhig, 
 * and @em yhig, containing the pixel coordinates of opposite corners of 
 * the overscan regions, and at its first row the corners coordinates
 * of the valid region (i.e., the CCD proper) within the image. This is 
 * a standard table valid for all instruments that should be used by 
 * the function @c mos_remove_bias(). In the specific case of VIMOS the 
 * keywords used are NAXIS1, NAXIS2, ESO DET OUT1 PRSCX, ESO DET OUT1 OVSCX, 
 * ESO DET OUT1 PRSCY, and ESO DET OUT1 OVSCY. The input header should come 
 * from a raw, unprocessed image (typically the image from which a master 
 * bias was not yet subtracted). 
 */

cpl_table *mos_load_overscans_vimos(const cpl_propertylist *header, 
                                    int check_consistency)
{
    const char *func = "mos_load_overscans_vimos";

    int        nx = 0;
    int        ny = 0;
    int        px = 0;
    int        py = 0;
    int        ox = 0;
    int        oy = 0;
    int        vx = 0;
    int        vy = 0;
    int        nrows;
    cpl_table *overscans;


    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(func, "Reset your error: %s", cpl_error_get_message());
        return NULL;
    }

    if (header == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (cpl_propertylist_has(header, "NAXIS1"))
        nx = cpl_propertylist_get_int(header, "NAXIS1");
    if (cpl_propertylist_has(header, "NAXIS2"))
        ny = cpl_propertylist_get_int(header, "NAXIS2");
    if (cpl_propertylist_has(header, "ESO DET OUT1 PRSCX"))
        px = cpl_propertylist_get_int(header, "ESO DET OUT1 PRSCX");
    if (cpl_propertylist_has(header, "ESO DET OUT1 PRSCY"))
        py = cpl_propertylist_get_int(header, "ESO DET OUT1 PRSCY");
    if (cpl_propertylist_has(header, "ESO DET OUT1 OVSCX"))
        ox = cpl_propertylist_get_int(header, "ESO DET OUT1 OVSCX");
    if (cpl_propertylist_has(header, "ESO DET OUT1 OVSCY"))
        oy = cpl_propertylist_get_int(header, "ESO DET OUT1 OVSCY");
    if (cpl_propertylist_has(header, "ESO DET OUT1 NX"))
        vx = cpl_propertylist_get_int(header, "ESO DET OUT1 NX");
    if (cpl_propertylist_has(header, "ESO DET OUT1 NY"))
        vy = cpl_propertylist_get_int(header, "ESO DET OUT1 NY");

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(func, "Missing overscan keywords in header");
        cpl_error_set_where(func);
        return NULL;
    }

    if (px < 0 || py < 0 || ox < 0 || oy < 0) {
        cpl_msg_error(func, "Missing overscan keywords in header");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if ((px + vx + ox != nx) || (py + vy + oy != ny)) {
        if (check_consistency) {
            cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
            return NULL;
        }
        else {
            cpl_msg_debug(func, "Overscans description conflicts with "
                          "reported image sizes, "
                          "%d + %d + %d != %d or "
                          "%d + %d + %d != %d",
                          px, vx, ox, nx,
                          py, vy, oy, ny);
        }
    }

    nrows = 0;
    if (px > 0)
        nrows++;
    if (ox > 0)
        nrows++;
    if (py > 0)
        nrows++;
    if (oy > 0)
        nrows++;

    if (nrows > 2) {
        cpl_msg_error(func, "Unexpected overscan regions "
                      "(both in X and Y direction)");
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }


    /*
     * A row is added for the description of the valid region of the
     * exposure the input header belongs to.
     */

    nrows++;

    overscans = cpl_table_new(nrows);
    cpl_table_new_column(overscans, "xlow", CPL_TYPE_INT);
    cpl_table_new_column(overscans, "ylow", CPL_TYPE_INT);
    cpl_table_new_column(overscans, "xhig", CPL_TYPE_INT);
    cpl_table_new_column(overscans, "yhig", CPL_TYPE_INT);

    nrows = 0;

    cpl_table_set_int(overscans, "xlow", nrows, px);
    cpl_table_set_int(overscans, "ylow", nrows, py);
    cpl_table_set_int(overscans, "xhig", nrows, nx - ox);
    cpl_table_set_int(overscans, "yhig", nrows, ny - oy);
    nrows++;

    if (px > 0) {
        cpl_table_set_int(overscans, "xlow", nrows, 0);
        cpl_table_set_int(overscans, "ylow", nrows, 0);
        cpl_table_set_int(overscans, "xhig", nrows, px);
        cpl_table_set_int(overscans, "yhig", nrows, ny);
        nrows++;
    }

    if (ox > 0) {
        cpl_table_set_int(overscans, "xlow", nrows, nx - ox);
        cpl_table_set_int(overscans, "ylow", nrows, 0);
        cpl_table_set_int(overscans, "xhig", nrows, nx);
        cpl_table_set_int(overscans, "yhig", nrows, ny);
        nrows++;
    }

    if (py > 0) {
        cpl_table_set_int(overscans, "xlow", nrows, 0);
        cpl_table_set_int(overscans, "ylow", nrows, 0);
        cpl_table_set_int(overscans, "xhig", nrows, nx);
        cpl_table_set_int(overscans, "yhig", nrows, py);
        nrows++;
    }

    if (oy > 0) {
        cpl_table_set_int(overscans, "xlow", nrows, 0);
        cpl_table_set_int(overscans, "ylow", nrows, ny - oy);
        cpl_table_set_int(overscans, "xhig", nrows, nx);
        cpl_table_set_int(overscans, "yhig", nrows, ny);
        nrows++;
    }

    return overscans;

}


cpl_table *mos_load_overscans_fors(const cpl_propertylist *header)
{
    const char *func = "mos_load_overscans_fors";

    int        nports;
    int        nx = 0;
    int        ny = 0;
    int        px = 0;
    int        py = 0;
    int        ox = 0;
    int        oy = 0;
    int        rebin;
    int        nrows;
    cpl_table *overscans;


    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(func, "Reset your error: %s", cpl_error_get_message());
        return NULL;
    }

    if (header == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (cpl_propertylist_has(header, "ESO DET OUTPUTS"))
        nports = cpl_propertylist_get_int(header, "ESO DET OUTPUTS");

    if (nports == 4                                        && 
        cpl_propertylist_has(header, "ESO DET OUT1 PRSCX") &&
        cpl_propertylist_has(header, "ESO DET WIN1 BINX")) {

        rebin = cpl_propertylist_get_int(header, "ESO DET WIN1 BINX");

        overscans = cpl_table_new(3);
        cpl_table_new_column(overscans, "xlow", CPL_TYPE_INT);
        cpl_table_new_column(overscans, "ylow", CPL_TYPE_INT);
        cpl_table_new_column(overscans, "xhig", CPL_TYPE_INT);
        cpl_table_new_column(overscans, "yhig", CPL_TYPE_INT);

        px = 16 / rebin;
        ox = 16 / rebin;
        nx = 2080 / rebin;
        ny = 2048 / rebin;
        nrows = 0;

        cpl_table_set_int(overscans, "xlow", nrows, px);
        cpl_table_set_int(overscans, "ylow", nrows, py);
        cpl_table_set_int(overscans, "xhig", nrows, nx - ox);
        cpl_table_set_int(overscans, "yhig", nrows, ny - oy);
        nrows++;

        cpl_table_set_int(overscans, "xlow", nrows, 0);
        cpl_table_set_int(overscans, "ylow", nrows, 0);
        cpl_table_set_int(overscans, "xhig", nrows, px);
        cpl_table_set_int(overscans, "yhig", nrows, ny);
        nrows++;

        cpl_table_set_int(overscans, "xlow", nrows, nx - ox);
        cpl_table_set_int(overscans, "ylow", nrows, 0);
        cpl_table_set_int(overscans, "xhig", nrows, nx);
        cpl_table_set_int(overscans, "yhig", nrows, ny);
        nrows++;
    }
    else {
        overscans = mos_load_overscans_vimos(header, 0);
    }

    return overscans;

}

/**
 * @brief
 *   Montecarlo simulation to evaluate error on polynomial fit.
 *
 * @param points   Table with (x,y) coordinates and uncertainties on y.
 * @param evaluate Table with x coordinates where to evaluate model variance.
 * @param samples  Number of simulations.
 * @param order    Degree of the fitted polynomial.
 *
 * @return Polynomial fit
 *
 * This function evaluates the effects of variations of coordinates on
 * a fitted polynomial model. The @em points table must contain two 
 * columns labeled @em x and @em y, and optionally a third column
 * labeled @em y_err. @em x is assumed errorless, and @em y_err must 
 * be given at 1-sigma error. All columns must be of type CPL_TYPE_DOUBLE. 
 * A first (reference) polynomial fit p(x) is made. If the column @em y_err 
 * was not given, one is created and filled with the RMS of the fit residuals. 
 * The fit is then repeated @em samples times, varying at each time the set 
 * y coordinates according to a gaussian random distribution. The varied y 
 * coordinates to fit are taken as 
 *
 *    y = p(x) + y_err * mos_randg(1)
 *
 * where p(x) is the fit on the original points, and @em x are taken 
 * from the @em points table together with the corresponding @em y_err. 
 * The @em evaluate table should contain a column labeled @em x, and 
 * is returned with a new column labeled @em sigma, listing the sigma 
 * of the model variation at each of the indicated @em x coordinates.
 */

#define READY 1
#ifdef READY

cpl_polynomial *mos_montecarlo_polyfit(cpl_table *points, cpl_table *evaluate, 
                                       int samples, int order)
{

    const char *func = "mos_montecarlo_polyfit";

    cpl_polynomial *p;
    cpl_polynomial *q;
    cpl_vector     *listx;
    cpl_vector     *listy;
    double          err;
    double         *x;
    double         *px;
    double         *x_eval;
    double         *px_eval;
    double         *sigma;
    double         *vy;
    double         *dy;
    int             npoints, nevaluate;
    int             i, j;


    if (points == NULL || evaluate == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (!cpl_table_has_column(points, "x")) {
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    if (cpl_table_get_column_type(points, "x") != CPL_TYPE_DOUBLE) {
        cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
        return NULL;
    }

    if (cpl_table_has_invalid(points, "x")) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (!cpl_table_has_column(points, "y")) {
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    if (cpl_table_get_column_type(points, "y") != CPL_TYPE_DOUBLE) {
        cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
        return NULL;
    }

    if (cpl_table_has_invalid(points, "y")) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (cpl_table_has_column(points, "y_err")) {

        if (cpl_table_get_column_type(points, "y_err") != CPL_TYPE_DOUBLE) {
            cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
            return NULL;
        }
    
        if (cpl_table_has_invalid(points, "y_err")) {
            cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
            return NULL;
        }
    }

    if (!cpl_table_has_column(evaluate, "x")) {
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    if (cpl_table_get_column_type(evaluate, "x") != CPL_TYPE_DOUBLE) {
        cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
        return NULL;
    }

    if (cpl_table_has_invalid(evaluate, "x")) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (samples < 2 || order < 0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    npoints = cpl_table_get_nrow(points);
    listx = cpl_vector_wrap(npoints, cpl_table_get_data_double(points, "x"));
    listy = cpl_vector_wrap(npoints, cpl_table_get_data_double(points, "y"));

    p = cpl_polynomial_fit_1d_create(listx, listy, order, &err);

    if (!cpl_table_has_column(points, "y_err")) {
        err = sqrt(err);
        cpl_table_new_column(points, "y_err", CPL_TYPE_DOUBLE);
        cpl_table_fill_column_window_double(points, "y_err", 0, npoints, err);
        cpl_msg_info(func, "Error column not found - set to %f\n", err);
    }

    /*
     * Create columns containing modeled values at each x
     */

    if (cpl_table_has_column(points, "px"))
        cpl_table_erase_column(points, "px");
    cpl_table_new_column(points, "px", CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(points, "px", 0, npoints, 0);
    x = cpl_table_get_data_double(points, "x");
    px = cpl_table_get_data_double(points, "px");
    for (i = 0; i < npoints; i++)
        px[i] = cpl_polynomial_eval_1d(p, x[i], NULL);

    nevaluate = cpl_table_get_nrow(evaluate);

    if (cpl_table_has_column(evaluate, "px"))
        cpl_table_erase_column(evaluate, "px");
    cpl_table_new_column(evaluate, "px", CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(evaluate, "px", 0, nevaluate, 0);
    x_eval = cpl_table_get_data_double(evaluate, "x");
    px_eval = cpl_table_get_data_double(evaluate, "px");
    for (i = 0; i < nevaluate; i++)
        px_eval[i] = cpl_polynomial_eval_1d(p, x_eval[i], NULL);

    /*
     * Initialise column with sigma
     */

    if (cpl_table_has_column(evaluate, "sigma"))
        cpl_table_erase_column(evaluate, "sigma");
    cpl_table_new_column(evaluate, "sigma", CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(evaluate, "sigma", 0, nevaluate, 0);
    sigma = cpl_table_get_data_double(evaluate, "sigma");

    /*
     * Compute varied y cordinates to fit
     */

    if (cpl_table_has_column(points, "vy"))
        cpl_table_erase_column(points, "vy");
    cpl_table_new_column(points, "vy", CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(points, "vy", 0, npoints, 0);
    vy = cpl_table_get_data_double(points, "vy");
    dy = cpl_table_get_data_double(points, "y_err");
    cpl_vector_unwrap(listy);
    listy = cpl_vector_wrap(npoints, vy);

    for (i = 0; i < samples; i++) {
        for (j = 0; j < npoints; j++)
            vy[j] = px[j] + dy[j] * mos_randg(1);
        q = cpl_polynomial_fit_1d_create(listx, listy, order, NULL);
        for (j = 0; j < nevaluate; j++)
            sigma[j] += fabs(px_eval[j] 
                      - cpl_polynomial_eval_1d(q, x_eval[j], NULL));
        cpl_polynomial_delete(q);
    }

    /* 
     * Factor 1.25 to convert average deviation to sigma 
     */

    cpl_table_multiply_scalar(evaluate, "sigma", 1.25);
    cpl_table_divide_scalar(evaluate, "sigma", samples);

    cpl_vector_unwrap(listx);
    cpl_vector_unwrap(listy);

    return p;
}

#endif

/**
 * @brief Randomise image
 *
 * @param image        Image to randomise (in place)
 * @param ron          Read out noise (ADU)
 * @param gain         Gain (electrons/ADU)
 * @param bias         Bias level (ADU)
 *
 * @return CPL_ERROR_NONE or corresponding cpl_error_code on error.
 *
 * For each pixel the 1-sigma error is evaluated as the square root of
 * the variance:
 * @code
                V = ron*ron + (S - bias) / gain
 * @endcode
 * where S is the pixel value.
 * The pixel value is randomised as
 * @code
                S = S + sqrt(V) * mos_randg(1)
 * @endcode
 */

cpl_error_code mos_randomise_image(cpl_image *image, double ron, 
                                   double gain, double bias)
{
    float *data;
    int    npix, i;


    if (image == NULL)
        return cpl_error_set(cpl_func, CPL_ERROR_NULL_INPUT);

    if (ron < 0.0 || gain <= FLT_EPSILON)
        return cpl_error_set(cpl_func, CPL_ERROR_ILLEGAL_INPUT);

    data = cpl_image_get_data_float(image);
    npix = cpl_image_get_size_x(image) * cpl_image_get_size_y(image);
    ron *= ron;

    for (i = 0; i < npix; i++) {
        if (data[i] < bias) {
            data[i] += sqrt(ron) * mos_randg(1);
        }
        else {
            data[i] += sqrt(ron + (data[i] - bias) / gain) * mos_randg(1);
        }
    }

    return CPL_ERROR_NONE;
}


/**
 * @brief Reconstruct the gaps required for slit location
 *
 * @param refmask      Reference mask 
 * @param master_flat  Masterflat. The gaps are detected and inserted based 
 *                     upon the values found in this image.
 *
 * @return CPL_ERROR_NONE of corresponding cpl_error_code on error.
 *
 * Deviation larger than 1 sigma from median are killed from map.
 * If @em level is greater than zero, deviations less than @em level from
 * median are killed.
 */

cpl_error_code mos_refmask_find_gaps(cpl_mask  *refmask,
                                     cpl_image *master_flat,
                                     double     level)
{
    int          nx     = cpl_mask_get_size_x(refmask);
    int          ny     = cpl_mask_get_size_y(refmask);

    int        * xpos   = cpl_calloc(sizeof(int), ny);

    cpl_image  * filtered = cpl_image_duplicate(master_flat);
    cpl_mask   * kernel = cpl_mask_new(9, 9);
    cpl_vector * v      = cpl_vector_new(ny);
    cpl_vector * truev;
    int          nvalid = 0;
    double     * flats  = cpl_vector_get_data(v);

    double       median, stdev, delta;

    int          i, kill;


    cpl_mask_not(kernel);
    cpl_image_filter_mask(filtered, master_flat, kernel, 
                          CPL_FILTER_MEDIAN, CPL_BORDER_COPY);
    cpl_mask_delete(kernel);

    for (i = 1; i <= ny; i++) {
        int j = 0;

        do j++;
        while (!cpl_mask_get(refmask, j, i) && j < nx);

        if (j < nx) {
            int rejected;

            xpos[i - 1] = j;
            flats[nvalid] = cpl_image_get(filtered, j, i, &rejected);
            nvalid++;
        }
        else {
            xpos[i - 1] = -1;
        }
    }

    if (nvalid == 0)
        return cpl_error_set(cpl_func, CPL_ERROR_DATA_NOT_FOUND);

    truev = cpl_vector_wrap(nvalid, flats);

    median = cpl_vector_get_median(truev);

    if (level < 0.0)
       stdev = cpl_vector_get_stdev(truev);

    cpl_vector_unwrap(truev);
    cpl_vector_delete(v);

    for (i = 1; i <= ny; i++) {
	if (xpos[i - 1] > 0) {
	    int    rejected;
	    double kappa = 1.0;

            delta = cpl_image_get(filtered, xpos[i - 1], i, &rejected) - median;

            if (level < 0.0)
                kill = fabs(delta) > stdev * kappa;
            else
                kill = delta < level;

            if (kill) {
                int j = 0;
            
                while (cpl_mask_get(refmask, xpos[i - 1] + j, i)) {
                    cpl_mask_set(refmask, xpos[i - 1] + j, i, CPL_BINARY_0);
                    j++;
                }
            }
        }
    }

    cpl_image_delete(filtered);
    cpl_free(xpos);

    return cpl_error_get_code();
}

/**
 * @brief Process saturation
 *
 * @param image        Image to process saturation
 *
 * @return CPL_ERROR_NONE of corresponding cpl_error_code on error.
 */
cpl_error_code mos_saturation_process(cpl_image * image)
{
    int     nx    = cpl_image_get_size_x(image);
    int     ny    = cpl_image_get_size_y(image);
    int     npix  = nx * ny;
    float * sdata = cpl_image_get_data_float(image);

    int count, i, j, k;

    /*
     * This is used to avoid saturation level coded with pixel value zero
     * To make it more robust against random 0.0 values, check that also
     * next pixel along the spatial direction is 0.0.
     */

    //This could be applied only to raw images, but it was being applied
    //to already bias/overscan processed images, which doesn't make sense.
//    for (i = 0; i < npix - nx; i++)
//        if (sdata[i] == 0.0 && sdata[i + nx] == 0.0)
//            sdata[i] = 65535.0;


//    for (i = npix - nx; i < npix; i++)
//        if (sdata[i] == 0.0) 
//            sdata[i] = 65535.0;

    /*
     * This is a dirty trick to overcome saturations (making up a false
     * tip on their flat tops). This should be useless with a better
     * peak detection algorithm.
     */

    for (i = 0; i < npix; i++) {
        if (sdata[i] >= 65535.0) {
            count = 0;
            for (j = i; j < npix; j++) {
                if (sdata[j] < 65535.0) {
                    break;
                }
                else {
                    count++;
                }
            }
            if (count < 30 && count > 2) {
                for (j = i; j < i + count/2; j++)
                    sdata[j] = sdata[i] + 1000.0 * (j - i);
                if (count % 2 != 0) {
                    sdata[j] = sdata[j-1] + 1000.0;
                    j++;
                }
                for (k = j; k <= i + count; k++)
                    sdata[k] = sdata[i] - 1000.0 * (k - i - count);
                i = k;
            }
        }
    }

    return cpl_error_get_code();
}


/**
 * @brief Subtract the background
 *
 * @param image        Image to subtract the background
 *
 * @return CPL_ERROR_NONE of corresponding cpl_error_code on error.
 */

cpl_error_code mos_subtract_background(cpl_image * image)
{
    /*
     * Create and subtract background
     */

    cpl_image * bimage = mos_arc_background(image, 15, 15);
    cpl_image_subtract(image, bimage);
    cpl_image_delete(bimage);

    return cpl_error_get_code();
}


/**
 * @brief Intersect a number of slit tables.
 *
 * @param slitss         Pointer to the list of slit tables
 * @param origslits      Pointer to output intersected table
 * @param nscience       Number of tables in the list
 * @param tolerance      Tolerance in object position (pixel)
 *
 * @return CPL_ERROR_NONE of corresponding cpl_error_code on error.

 * The intersected table will only contain polarimetric objects which have 
 * been detected in all tables in both beams in the input table list 
 * (that is, are present at all angles so polarimetric computation 
 * is possible)
 */

cpl_error_code mos_object_intersect(cpl_table **slitss, cpl_table *origslits, 
                                    int nscience, float tolerance)
{
    int i, j;

    cpl_table *summary;
    int summary_nobjs = 0;
 
    int nobjs;

    int nmatches;
    int nslits = cpl_table_get_nrow(slitss[0]);

    int maxobjs;
    int k, m;
    int nstokes, sstokes;

    cpl_table **work;

    work = (cpl_table **)cpl_malloc(sizeof(cpl_table *) * nscience);


    /* 
     * First we build a table listing the offset of each detected
     * object at each angle and each beam, from the bottom of each 
     * slit spectrum, and the pair that slit spectrum belongs to.
     * This summary table will have as many rows as objects found 
     * in total at all angles.
     */

    for (j = 0; j < nscience; j++) {
        int c_nobjs = mos_get_nobjects(slitss[j]);
        if (!c_nobjs) 
            return cpl_error_set(cpl_func, CPL_ERROR_DATA_NOT_FOUND);
        summary_nobjs += c_nobjs;
    }

    summary = cpl_table_new(summary_nobjs);

    cpl_table_new_column(summary, "offset", CPL_TYPE_DOUBLE);
    cpl_table_new_column(summary, "pair",   CPL_TYPE_INT);
    cpl_table_new_column(summary, "absolute", CPL_TYPE_DOUBLE);
    cpl_table_new_column(summary, "pos", CPL_TYPE_DOUBLE);

    /*
     * Fill the summary table with data from all objects:
     */

    nobjs = 0;

    /* Loop on all object tables (one for each angle) */
    for (j = 0; j < nscience; j++) {
        int c_maxobjs = mos_get_maxobjs_per_slit(slitss[j]);

        /* Loop on all slits found on first - i.e., ALL - object table */
        for (k = 0; k < nslits; k++) {

            /* Loop on all objects found on each object table */
            for (m = 0; m < c_maxobjs; m++) {
                int null;
                char *name = cpl_sprintf("object_%d", m + 1);
                double obj = cpl_table_get_double(slitss[j], name, k, &null);
                int pos;
                int pair;

                cpl_free(name);

                if (null) 
                    break;  /* No object #m+1 in this slit - go to next slit */

                /*
                 * Copy necessary object data to summary table. Note 
                 * that the absolute object position (row) in the
                 * rectified image is made relative to the bottom
                 * position (row) of the current slit.
                 */ 
        
                pos  = cpl_table_get_int(slitss[j], "position", k, &null);
                pair = cpl_table_get_int(slitss[j], "pair_id", k, &null);
                cpl_table_set(summary, "absolute", nobjs, obj);
                cpl_table_set(summary, "pos", nobjs, pos);
                cpl_table_set(summary, "offset", nobjs, obj - pos);
                cpl_table_set(summary, "pair", nobjs, pair);

                nobjs++;
            }
        }
    }

//    cpl_table_save(summary, NULL, NULL, "susu.fits", CPL_IO_DEFAULT);

    /* 
     * Perform the intersection: what are the objects belonging
     * to the same slit (same pair ordinary + extraordinary) which 
     * are observed at the same offset at all angles? Those are
     * the polarimetric objects.
     */

    nmatches = 0;
    maxobjs = mos_get_maxobjs_per_slit(slitss[0]);

    /*
     * We loop on the objects of the first-angle object table as 
     * reference, and check whether those objects are present also
     * at *all* other angles. Note that the loop advances by pairs.
     * If the top (k = 0) slit spectrum is not an ordinary beam,
     * it is ignored. The loop advances by pairs, starting at the
     * first complete pair. It is implicitely assumed that the 
     * slit spectrum on top is always from the ordinary beam, and 
     * the spectrum below (k+1) its extraordinary match.
     */

    for (k = 0; k < nslits; k+=2) {
        int slitmatches = 0;

        if (k == 0) {
            if (cpl_table_get_int(slitss[0], "pair_id",  0, NULL) !=
                cpl_table_get_int(slitss[0], "pair_id",  1, NULL)) {

                /*
                 * This is not an ordinary beam - advance to next slit.
                 */

                k++;
                continue;
            }
        }

	for (m = 0; m < maxobjs; m++) {
	    int null;
	    char *name = cpl_sprintf("object_%d", m + 1);
	    double obj = cpl_table_get_double(slitss[0], name, k, &null);
	    double pos;
            int pair;

	    char *name_obj;
	    char *name_start;
	    char *name_end;
            char *name_row;
            char *name_row_s;

	    char *name_start_o;
	    char *name_end_o;
            char *name_row_o;

            int start, end;
            int length;
 
            int selected;


	    cpl_free(name);

	    if (null) 
                break;

            /*
             * Each object of the first object table belongs to a
             * slit spectrum (k). This slit spectrum has a position
             * in the rectified image, and it belongs to a given 
             * ordinary + extraordinary pair.
             */
     
            pos  = cpl_table_get_int(slitss[0], "position", k, &null);
            pair = cpl_table_get_int(slitss[0], "pair_id",  k, &null);

            /*
             * Now from the summary table we can select all objects
             * which have the same offset (obj - pos) within all slit
             * spectra belonging to the same ordinary + extraordinary 
             * pair (at all angles).
             */

	    cpl_table_select_all(summary);  /* Reset selection */

            cpl_table_and_selected_int(summary, "pair", CPL_EQUAL_TO, pair);
            cpl_table_and_selected_double(summary, "offset", CPL_LESS_THAN,
                                          obj - pos + tolerance);
	    selected = 
            cpl_table_and_selected_double(summary, "offset", CPL_GREATER_THAN,
                                          obj - pos - tolerance);


            /*
             * If this object were observed at all angles (nscience) and 
             * at all beams (2), we should have selected exactly 2*nscience
             * objects. If not, this is not a polarimetric object, and it
             * is discarded from the intersection.
             */
	    
	    if (selected != nscience * 2) 
                continue;

            /*
             * If we reach this point we have found one valid polarimetric
             * object, that must be inserted in the intersection object
             * table.
             */
 
            slitmatches++;

            /*
             * Names of the columns of the output table where the
             * object information needs to be copied. Note that a
             * new column is created, the "row_stokes_#", where the
             * row number of the extracted polarimetric signal is
             * also computed. For the moment this column will be 
             * left empty - it will be filled only when all matches 
             * are collected.
             */

	    name_obj   = cpl_sprintf("object_%d",     slitmatches);
	    name_start = cpl_sprintf("start_%d",      slitmatches);
	    name_end   = cpl_sprintf("end_%d",        slitmatches);
	    name_row   = cpl_sprintf("row_%d",        slitmatches);
	    name_row_s = cpl_sprintf("row_stokes_%d", slitmatches);

            /*
             * Names of the columns of the input table where the
             * object information is available.
             */

	    name_start_o = cpl_sprintf("start_%d",  m + 1);
	    name_end_o   = cpl_sprintf("end_%d",    m + 1);
	    name_row_o   = cpl_sprintf("row_%d",    m + 1);

            /*
             * If the output columns do not exist yet, create them.
             */
 
	    if (!cpl_table_has_column(origslits, name_obj)) {
		cpl_table_new_column(origslits, name_obj, CPL_TYPE_DOUBLE);
         	cpl_table_new_column(origslits, name_start, CPL_TYPE_INT);
		cpl_table_new_column(origslits, name_end,   CPL_TYPE_INT);
		cpl_table_new_column(origslits, name_row,   CPL_TYPE_INT);
		cpl_table_new_column(origslits, name_row_s, CPL_TYPE_INT);
            }

            /*
             * The current slit spectrum is k. The slit spectrum immediately
             * below (in the rectified image) is k+1. We need the length of
             * the spectrum below for computing the _absolute_ coordinates
             * of the objects in the rectified image in both beams.
             */
 
	    length = cpl_table_get_int(origslits, "length", k + 1, &null);

            /*
             * Read from the first input object table (first angle)
             * the spatial window enclosing the object.
             */

	    start = cpl_table_get_int(slitss[0], name_start_o, k, &null);
	    end   = cpl_table_get_int(slitss[0], name_end_o,   k, &null);

            /*
             * Write the object coordinates in the same slit, and in the
             * slit below. Note that here we assume that all slits were
             * traced perfectly, and we compute the theoretical coords
             * (obj - length) within the next slit spectrum (k + 1). In
             * principle we should read them as well from the input
             * table!
             */

            cpl_table_set_double(origslits, name_obj,   k,     obj);
            cpl_table_set_double(origslits, name_obj,   k + 1, obj - length);

            cpl_table_set_int(origslits,    name_start, k,     start);
            cpl_table_set_int(origslits,    name_start, k + 1, start - length);

            cpl_table_set_int(origslits,    name_end,   k,     end);
            cpl_table_set_int(origslits,    name_end,   k + 1, end - length);

            /*
             * "nmatches" is counting at what "reduced" image row the
             * extracted spectra are. Note that this is s preliminary
             * numbering - which is wrong: other objects may be found
             * in the same slit, and then the indeces would not be in
             * sequence. What is important is that at the end of this
             * loop "nmatches" would be the total number of matching 
             * objects. The two cpl_table_set_int() calls made here
             * cannot be removed - they "validate" those table elements
             * (see ahead). 
             */

            cpl_table_set_int(origslits,    name_row,   k,     nmatches);
	    nmatches++;
            cpl_table_set_int(origslits,    name_row,   k + 1, nmatches);
	    nmatches++;

	    cpl_free(name_obj);
	    cpl_free(name_start);
	    cpl_free(name_end);
            cpl_free(name_row);
            cpl_free(name_row_s);

	    cpl_free(name_start_o);
	    cpl_free(name_end_o);
            cpl_free(name_row_o);
        }
    }

    /*
     * The summary table has fulfilled its function. If no matching 
     * objects are found, the function returns with an error.
     */

    cpl_table_delete(summary);

    if (!nmatches)
        return cpl_error_set(cpl_func, CPL_ERROR_DATA_NOT_FOUND); 

    /*
     * Now we consider the resulting intersection object table,
     * listing all matches. As seen, the image row number reported
     * in the columns "row_#" was not really performed sequentially.
     * We need to renumber sequentially...
     * We need also to fill the "row_stokes_#" column the way the
     * extracted polarimetric signal will be stored in the 
     * reduced_pol_images...
     */
 
    maxobjs = mos_get_maxobjs_per_slit(origslits);
    nstokes = nmatches / 2;         /* nmatches is always an even number     */

    for (k = 0; k < nslits; k++) {
        if (k % 2) { /* Extraordinary beam */
            nstokes = sstokes;      /* Use same start value as for ordinary  */
        }
        else {       /* Ordinary beam      */
            sstokes = nstokes;      /* Memorise start value at ordinary beam */
        }

	for (m = 0; m < maxobjs; m++) {
	    char *name       = cpl_sprintf("row_%d",        m + 1);
            char *namestokes = cpl_sprintf("row_stokes_%d", m + 1);

	    if (!cpl_table_is_valid(origslits, name, k)) {
	        cpl_free(name);
                cpl_free(namestokes);
	        break;
	    }
            else { 
                nmatches--;
                nstokes--;
	        cpl_table_set_int(origslits, name, k, nmatches);
                cpl_table_set_int(origslits, namestokes, k, nstokes);
	    }

	    cpl_free(name);
            cpl_free(namestokes);
	}
    }


    /*
     * This is done to avoid the NULL value is zero (it would invalidate
     * also the row_# = 0 or start_# = 0 for an object), and to enable 
     * working directly with the column data buffers, when using this 
     * table afterwards.
     */

    for (j = 0; j < maxobjs; j++) {
	char *name = cpl_sprintf("object_%d", j + 1);
	cpl_table_fill_invalid_double(origslits, name, -1);
	cpl_free(name);

	name       = cpl_sprintf("start_%d", j + 1);
	cpl_table_fill_invalid_int(origslits, name, -1);
	cpl_free(name);

	name       = cpl_sprintf("end_%d", j + 1);
	cpl_table_fill_invalid_int(origslits, name, -1);
	cpl_free(name);

	name       = cpl_sprintf("row_%d", j + 1);
	cpl_table_fill_invalid_int(origslits, name, -1);
	cpl_free(name);

	name       = cpl_sprintf("row_stokes_%d", j + 1);
	cpl_table_fill_invalid_int(origslits, name, -1);
	cpl_free(name);
    }

    /*********************************************************************
     * This tail has been added to propagate the selection of valid
     * objects also to the input slitss[] tables. Just eliminate all
     * this final part to suppress this behaviour.
     */

    /*
     * First of all, make a working copy and remove all columns related 
     * to objects from the input object tables. 
     */

    for (i = 0; i < nscience; i++) {
        int c_maxobjs = mos_get_maxobjs_per_slit(slitss[i]);

        work[i] = cpl_table_duplicate(slitss[i]);

        for (m = 0; m < c_maxobjs; m++) {
            char *object_o = cpl_sprintf("object_%d", m + 1);
            char *start_o  = cpl_sprintf("start_%d",  m + 1);
            char *end_o    = cpl_sprintf("end_%d",    m + 1);
            char *row_o    = cpl_sprintf("row_%d",    m + 1);

            cpl_table_erase_column(slitss[i], object_o);
            cpl_table_erase_column(slitss[i], start_o);
            cpl_table_erase_column(slitss[i], end_o);
            cpl_table_erase_column(slitss[i], row_o);
        }
    }

    /* 
     * Now just consider all the objects in the intersection table.
     */

    for (k = 0; k < nslits; k++) {
        for (j = 0; j < maxobjs; j++) {
            double object_w, object_r;
            int    row_w;

	    char  *object_i = cpl_sprintf("object_%d", j + 1);
	    char  *start_i  = cpl_sprintf("start_%d",  j + 1);
	    char  *end_i    = cpl_sprintf("end_%d",    j + 1);
	    char  *row_i    = cpl_sprintf("row_%d",    j + 1);


            if (!cpl_table_is_valid(origslits, object_i, k))
                break;

            /* 
             * We have found a valid object (valid because it belongs
             * to the intersection). Now we look for this object in each
             * one of the original tables, we get its parameters, and
             * copy them at the right position (i.e., same position as
             * in intersection table). The object will be the one closest
             * to the object position (column object_i) in the intersection
             * table. Note that we examine the same row, k, in all tables.
             */

            object_w = cpl_table_get_double(origslits, object_i, k, NULL);
            row_w    = cpl_table_get_int   (origslits, row_i,    k, NULL);

            for (i = 0; i < nscience; i++) {
                int        c_maxobjs = mos_get_maxobjs_per_slit(work[i]);
                int        minpos;
                double     mindiff, diff;
	        char      *object_o;
	        char      *start_o;
	        char      *end_o;
	        char      *row_o;

                for (m = 0; m < c_maxobjs; m++) {
	            object_o = cpl_sprintf("object_%d", m + 1);
	            start_o  = cpl_sprintf("start_%d",  m + 1);
	            end_o    = cpl_sprintf("end_%d",    m + 1);
	            row_o    = cpl_sprintf("row_%d",    m + 1);

                    if (!cpl_table_is_valid(work[i], object_o, k))
                        break;

                    object_r = cpl_table_get_double(work[i], object_o, k, NULL);

                    diff = fabs(object_w - object_r);
                    if (m) {
                        if (mindiff > diff) {
                            mindiff = diff;
                            minpos = m;
                        }
                    }
                    else {
                        mindiff = diff;
                        minpos = 0;
                    }

	            cpl_free(object_o);
	            cpl_free(start_o);
	            cpl_free(end_o);
	            cpl_free(row_o);
                }

                object_o = cpl_sprintf("object_%d", minpos + 1);
                start_o  = cpl_sprintf("start_%d",  minpos + 1);
                end_o    = cpl_sprintf("end_%d",    minpos + 1);
                row_o    = cpl_sprintf("row_%d",    minpos + 1);

                if (!cpl_table_has_column(slitss[i], object_i)) {
                    cpl_table_new_column(slitss[i], object_i, CPL_TYPE_DOUBLE);
                    cpl_table_new_column(slitss[i], start_i,  CPL_TYPE_INT);
                    cpl_table_new_column(slitss[i], end_i,    CPL_TYPE_INT);
                    cpl_table_new_column(slitss[i], row_i,    CPL_TYPE_INT);
	            cpl_table_fill_invalid_double(slitss[i], object_i, -1);
	            cpl_table_fill_invalid_int   (slitss[i], start_i,  -1);
	            cpl_table_fill_invalid_int   (slitss[i], end_i,    -1);
	            cpl_table_fill_invalid_int   (slitss[i], row_i,    -1);
                }

                cpl_table_set_double(slitss[i], object_i, k,
                                     cpl_table_get_double(work[i], object_o, 
                                                          k, NULL));
                cpl_table_set_int(slitss[i], start_i , k,
                                  cpl_table_get_int(work[i], start_o, k, NULL));
                cpl_table_set_int(slitss[i], end_i , k,
                                  cpl_table_get_int(work[i], end_o, k, NULL));
                cpl_table_set_int(slitss[i], row_i , k, row_w);

	        cpl_free(object_o);
	        cpl_free(start_o);
	        cpl_free(end_o);
	        cpl_free(row_o);
            }

	    cpl_free(object_i);
	    cpl_free(start_i);
	    cpl_free(end_i);
	    cpl_free(row_i);
        }
    }

    for (i = 0; i < nscience; i++)
        cpl_table_delete(work[i]);

    cpl_free(work);


    return cpl_error_get_code();
}


/**
 * @brief Get the maximum possible number of objects in a slit
 *
 * @param slits   Slits table
 *
 * @return maximum possible number of objects in a slit
 */
int mos_get_maxobjs_per_slit(cpl_table * slits)
{
    int maxobjs = 1;

    char * colname = cpl_sprintf("object_%d", maxobjs);
    
    while (cpl_table_has_column(slits, colname)) {
        maxobjs++;
        cpl_free(colname);
        colname = cpl_sprintf("object_%d", maxobjs);
    }
    
    cpl_free(colname);

    maxobjs--;

    return maxobjs;
}

/**
 * @brief Get the total number of objects detected in a slits table
 *
 * @param slits   Slits table
 *
 * @return the number of objects in the table
 */
int mos_get_nobjects(cpl_table * slits)
{
    int nobjs = 0;

    int nslits  = cpl_table_get_nrow(slits);
    int maxobjs = mos_get_maxobjs_per_slit(slits);

    int k, m;

    for (k = 0; k < nslits; k++) {
        for (m = 0; m < maxobjs; m++) {
            char * name = cpl_sprintf("object_%d", m + 1);
            int    null = !cpl_table_is_valid(slits, name, k);

            cpl_free(name);

            if (null)  break;
            else nobjs++;
        }
    }

    return nobjs;
}

/**
 * @brief Check that all slit have been detected, insert them if not
 *
 * @param slits   Slits table
 *
 * @return 0 or -1 on error.
 */
int mos_check_slits(cpl_table *slits, float rescale)
{

    cpl_propertylist *sort;

    int nslits  = cpl_table_get_nrow(slits);

    int k, null;

    const float interval = 90.0 * rescale;
    const float offset   = (90.0 - 5) * rescale;


    for (k = 0; k < nslits; k++) {
        double ytop    = cpl_table_get_double(slits, "ytop",    k, &null);
        double ybottom = cpl_table_get_double(slits, "ybottom", k, &null);

        double xtop    = cpl_table_get_double(slits, "xtop",    k, &null);
        double xbottom = cpl_table_get_double(slits, "xbottom", k, &null);

        int nmiss = (int)((ytop - ybottom) / interval + 0.5);

        if (nmiss > 1) {
            cpl_msg_warning(cpl_func, 
                            "Some slits could not be properly detected. "
                            "There might be accountable inaccuracies.");
            while (nmiss > 1) {
                cpl_table_set_size(slits, nslits + 1);

                /* Fill in new slit 'cut' */

                /* x coordinates be the same (acceptable approximation) */
                cpl_table_set_double(slits, "xtop",    nslits, xtop);
                cpl_table_set_double(slits, "xbottom", nslits, xbottom);

                /* Cut */
                if (k == 0) {
                    cpl_table_set_double(slits, "ybottom", nslits, ybottom); 
                    cpl_table_set_double(slits, "ytop",    nslits, ybottom
                                                                   + offset);
                    ybottom += interval;
                    cpl_table_set_double(slits, "ybottom", k,      ybottom);
                } else {
                    cpl_table_set_double(slits, "ytop",    nslits, ytop);
                    cpl_table_set_double(slits, "ybottom", nslits, ytop 
                                                                   - offset);
                    ytop -= interval;
                    cpl_table_set_double(slits, "ytop",     k,     ytop);
                }

                nslits++; nmiss--;
            }
        }
    }

    sort = cpl_propertylist_new();
    cpl_propertylist_append_bool(sort, "ytop", 1);
    cpl_table_sort(slits, sort);
    cpl_propertylist_delete(sort);

    /*
     * Add here an ad hoc check on the last slit: is it too long 
     * (by more than 10%)? Then shorten it...
     */

    k = cpl_table_get_nrow(slits) - 1;

    {
        double ytop    = cpl_table_get_double(slits, "ytop",    k, &null);
        double ybottom = cpl_table_get_double(slits, "ybottom", k, &null);
        double length  = (ytop - ybottom) / interval;

        if (length > 1.1) {
            cpl_table_set_double(slits, "ybottom", k, ytop - offset);
        }
  
    }

    return 0;
}

/**
 * @brief
 *   Create PMOS slit location table from FITS header of FORS1/2 MOS data
 *
 * @param header    FITS header of FORS data containing slits information
 *
 * @return A slit position table
 *
 * This function is meant to convert the information contained in 
 * FORS1 and FORS2 MOS data FITS header into a slit position table. 
 * This table will contain the double precision columns labeled 
 * @em xtop, @em ytop, @em xbottom, and @em ybottom, containing the 
 * start and end coordinates of the slits on the telescope focal 
 * plane (mask). The coordinates are expected to have a horizontal 
 * dispersion direction and red wavelengths on the right side.
 * A flip to the Y coordinate is applied to match the increasing
 * Y CCD pixel coordinate. A slit identifying integer column, 
 * labeled "slit_id", containing unique slit identifiers, will 
 * also be added. In case of FORS2 data only the slits pertaining
 * to the used chip are loaded.
 */

cpl_table *mos_load_slits_fors_pmos(cpl_propertylist *header)
{
    int m, null;
    int halfsize;

    cpl_propertylist * sort;
    cpl_table        * slits; 

    slits    = mos_load_slits_fors_mos(header);
    halfsize = cpl_table_get_nrow(slits);

    cpl_table_set_size(slits, 2 * halfsize);

    for (m = 0; m < halfsize; m++) {

        double gap = 1.4;

        double length = 
            cpl_table_get(slits, "ytop",    m, &null) -
            cpl_table_get(slits, "ybottom", m, &null);

        if (m) {
            double interval = 
                cpl_table_get(slits, "ybottom", m - 1, &null) -
                cpl_table_get(slits, "ytop",    m,     &null);

            gap = (interval - length) / 2;
        }

        cpl_table_set(slits, "slit_id", m + halfsize,
                      cpl_table_get(slits, "slit_id", m, &null) - 1);

        cpl_table_set(slits, "xtop",    m + halfsize,
                      cpl_table_get(slits, "xtop",    m, &null));

        cpl_table_set(slits, "xbottom", m + halfsize,
                      cpl_table_get(slits, "xbottom", m, &null));

        cpl_table_set(slits, "ytop",    m + halfsize, 
                      cpl_table_get(slits, "ytop", m, &null) + gap + length);

        cpl_table_set(slits, "ybottom", m + halfsize,
                      cpl_table_get(slits, "ytop", m, &null) + gap);
    }

    for (m = 0; m < 2 * halfsize; m++) {
        cpl_table_set(slits, "ytop",    m, 
                      cpl_table_get(slits, "ytop",    m, &null) - 5.3);

        cpl_table_set(slits, "ybottom", m,
                      cpl_table_get(slits, "ybottom", m, &null) - 5.3);

    }

    sort = cpl_propertylist_new();
    cpl_propertylist_append_bool(sort, "ytop", 1);
    cpl_table_sort(slits, sort);

    cpl_propertylist_delete(sort);

    return slits;
}

int * fors_get_nobjs_perslit(cpl_table * slits)
{
    int nslits  = cpl_table_get_nrow(slits);
    int maxobjs = mos_get_maxobjs_per_slit(slits);

    int * nobjs_per_slit = cpl_malloc(sizeof(int) * nslits);

    int k, m;

    for (k = 0; k < nslits; k++) {
        int nobjs = 0;
        for (m = 0; m < maxobjs; m++) {
            char * name = cpl_sprintf("object_%d", m + 1);
            int    null = !cpl_table_is_valid(slits, name, k);

            cpl_free(name);

            if (null)  break;
            else nobjs++;
        }
        
        nobjs_per_slit[k] = nobjs;
    }

    return nobjs_per_slit;
}

double fors_get_object_position(cpl_table *slits, int slit, int object)
{
    char   *name = cpl_sprintf("object_%d", object);
    double  position;

    position = cpl_table_get_double(slits, name, slit, NULL)
             - cpl_table_get_int(slits, "position", slit, NULL);

    cpl_free(name);

    return position;
}

int mos_rebin_signal(cpl_image **image, int rebin)
{
    cpl_image *rebinned;


    if (*image == NULL)
        return 1;

    if (rebin == 1)
        return 0;

    rebinned = cpl_image_rebin(*image, 1, 1, rebin, 1);

    cpl_image_delete(*image);

    *image = rebinned;

    return 0;
}

int mos_rebin_error(cpl_image **image, int rebin)
{
    if (*image == NULL)
        return 1;

    if (rebin == 1)
        return 0;

    cpl_image_power(*image, 2);
    mos_rebin_signal(image, rebin);
    cpl_image_power(*image, 0.5);

    return 0;
}

/*
 * @brief
 *   Map table values into a 1D image
 *  
 * @param image       Target image
 * @param start       Coordinate of first pixel in image
 * @param step        Coordinate step for one pixel in image
 * @param table       Source table
 * @param xname       Name of coordinate column
 * @param yname       Name of values column
 *
 * @return 0 on success
 *
 * The values in @em yname are linearly interpolated at the @em image 
 * pixel coordinates. The @em image must have Nx1 size.
 */

int map_table(cpl_image *image, double start, double step,
              cpl_table *table, const char *xname, const char *yname)
{
    int      length = cpl_image_get_size_x(image);
    int      nrows  = cpl_table_get_nrow(table);
    float   *data   = cpl_image_get_data_float(image);
    float   *fdata  = NULL;
    double  *xdata  = NULL;
    double  *ydata  = NULL;
    cpl_type xtype  = cpl_table_get_column_type(table, xname);
    cpl_type ytype  = cpl_table_get_column_type(table, yname);
    double   xzero, pos;
    int      i, j, n;


    /*
     * Initialization of output image at 0.0 - this value is left 
     * on non-overlapping portions.
     */

    for (i = 0; i < length; i++)
        data[i] = 0.0;


    /*
     * Do everything in double precision
     */

    if (xtype == CPL_TYPE_FLOAT) {
        fdata = cpl_table_get_data_float(table, xname);
        xdata = cpl_malloc(nrows * sizeof(double));
        for (i = 0; i < nrows; i++) {
           xdata[i] = fdata[i];
        }
    }
    else {
        xdata = cpl_table_get_data_double(table, xname);
    }

    if (ytype == CPL_TYPE_FLOAT) {
        fdata = cpl_table_get_data_float(table, yname);
        ydata = cpl_malloc(nrows * sizeof(double));
        for (i = 0; i < nrows; i++) {
           ydata[i] = fdata[i];
        }
    }
    else {
        ydata = cpl_table_get_data_double(table, yname);
    }

    /*
     * Mapping
     */

    n = 0;
    xzero = xdata[n];

    for (i = 0; i < length; i++) {
        pos = start + step * i;
        if (pos < xzero)
            continue;
        for (j = n; j < nrows; j++) {
            if (xdata[j] > pos) {
                n = j;
                data[i] = ydata[j-1]
                        + (ydata[j] - ydata[j-1])
                        * (pos - xdata[j-1]) / (xdata[j] - xdata[j-1]);
                break;
            }
        }
    }

    if (xtype == CPL_TYPE_FLOAT)
        cpl_free(xdata);

    if (ytype == CPL_TYPE_FLOAT)
        cpl_free(ydata);

    return 0;
}


/*
 * @brief
 *   Fit overall trend of a Nx1 image
 *  
 * @param image       Values to smooth
 * @param order       Order of fitting polynomial
 * @param hw          Half width of smoothing window
 *
 * @return Smoothed image, or NULL on failure.
 *
 * Heavily smooth and fit data in the input Nx1 size @em image.
 */

static cpl_image *polysmooth(cpl_image *image, int order, int hw)
{
    int             npoints;
    cpl_vector     *x;
    cpl_vector     *y;
    double         *xdata;
    double         *ydata;
    cpl_polynomial *poly;
    cpl_vector     *ysmooth;
    cpl_image      *smoothed;
    float          *sdata;
    int             i;


    npoints = cpl_image_get_size_x(image);

    if (2 * hw + 1 > npoints)
        return NULL;

    x       = cpl_vector_new(npoints);
    y       = cpl_vector_new(npoints);
    xdata   = cpl_vector_get_data(x);
    ydata   = cpl_vector_get_data(y);

    smoothed = cpl_image_duplicate(image);
    sdata = cpl_image_get_data_float(smoothed);

    for (i = 0; i < npoints; i++) {
        xdata[i] = i;
        ydata[i] = sdata[i];
    }

    ysmooth = cpl_vector_filter_median_create(y, hw);
    cpl_vector_delete(y);

    poly = cpl_polynomial_fit_1d_create(x, ysmooth, order, NULL);
    cpl_vector_delete(x);
    cpl_vector_delete(ysmooth);

    if (poly) {
        for (i = 0; i < npoints; i++)
            sdata[i] = cpl_polynomial_eval_1d(poly, i, NULL);
    
        cpl_polynomial_delete(poly);
    }
    else {
        cpl_image_delete(smoothed);
        return NULL;
    }

    return smoothed;
}

#undef cleanup
#define cleanup                       \
do {                                  \
    cpl_image_delete(spectrum);       \
    cpl_image_delete(flux);           \
    cpl_image_delete(efficiency);     \
    cpl_image_delete(smo_efficiency); \
    cpl_image_delete(extinction);     \
    cpl_image_delete(response);       \
    cpl_image_delete(smo_response);   \
    cpl_image_delete(physical);       \
} while (0)

/**
 * @brief
 *   Produce instrument response curve, with some ancillary information
 *
 * @param spectra     Image with extracted spectra
 * @param startwave   Shortest wavelength in input image (Angstrom)
 * @param dispersion  Angstrom per pixel of input image
 * @param gain        Gain factor (e-/ADU)
 * @param exptime     Exposure time (seconds)
 * @param ext_table   Atmospheric extinction table
 * @param airmass     Airmass of observation
 * @param flux_table  Standard star catalog flux
 * @param order       Order of polynomial to model the instrument response
 *
 * @return A slit position table
 *
 * This function is meant to convert the information contained in 
 *
 * The brightest extracted spectrum in @em spectra is assumed to be
 * the standard star spectrum.
 * The @em order of the polynomial fitting the instrument must be at least 2.
 */

cpl_table *mos_photometric_calibration(cpl_image *spectra, double startwave, 
                                 double dispersion, double gain,
                                 double exptime, cpl_table *ext_table,
                                 double airmass, cpl_table *flux_table,
                                 int order)
{

    cpl_image *spectrum       = NULL; // Extracted standard star spectrum
    float     *data;
    cpl_image *extinction     = NULL; // Extinction binned as "spectrum"
    float     *ext_data;
    cpl_image *flux           = NULL; // Standard star flux binned as "spectrum"
    float     *flux_data;
    cpl_image *physical       = NULL; // Physical units of above
    float     *phys_data;
    cpl_image *efficiency     = NULL; // Raw efficiency curve
    float     *eff_data;
    cpl_image *smo_efficiency = NULL; // Smoothed efficiency curve
    float     *smo_eff_data;
    cpl_image *response       = NULL; // Raw response curve
    float     *res_data;
    cpl_image *smo_response   = NULL; // Smoothed response curve
    float     *smo_res_data;
    cpl_image *image;
    cpl_image *smo_image;
    cpl_table *table;
    float      lambda;
    int        nx, ny;
    int        ext_count, ext_pos;
    int        eff_count, eff_pos;
    int        flux_count, flux_pos;
    int        start, end;
    int        i;


    if (spectra == NULL || ext_table == NULL || flux_table == NULL) {
        cpl_error_set(cpl_func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (!cpl_table_has_column(ext_table, "WAVE")) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                              "Column WAVE in atmospheric extinction table");
        return NULL;
    }

    if (!cpl_table_has_column(ext_table, "EXTINCTION")) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                        "Column EXTINCTION in atmospheric extinction table");
        return NULL;
    }

    if (!cpl_table_has_column(flux_table, "WAVE")) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                              "Column WAVE in standard star flux table");
        return NULL;
    }

    if (!cpl_table_has_column(flux_table, "FLUX")) {
        cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                              "Column FLUX in standard star flux table");
        return NULL;
    }

    if (gain < 0.1) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Invalid gain factor (%.2f)", gain);
        return NULL;
    }

    if (exptime < 0.001) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Invalid exposure time (%.2f)", exptime);
        return NULL;
    }

    if (dispersion < 0.001) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Invalid dispersion (%.2f)", dispersion);
        return NULL;
    }

    if (order < 2) {
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Order of the polynomial fitting the "
                              "instrument response must be at least 2");
        return NULL;
    }

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);


    /*
     * Find brightest spectrum and duplicate it.
     */

    if (ny == 1) {
        spectrum = cpl_image_duplicate(spectra);
    }
    else {
        cpl_size        x, y;
        cpl_image *brights = cpl_image_collapse_create(spectra, 1);

        cpl_image_get_maxpos(brights, &x, &y);
        cpl_image_delete(brights);
        spectrum = cpl_image_extract(spectra, 1, y, nx, y);
    }


    /*
     * Convert standard star spectrum in electrons per second per Angstrom.
     */

    cpl_image_multiply_scalar(spectrum, gain / exptime / dispersion);


    /*
     * Map the atmospheric extinction factors to the same lambda sampling
     * of the extracted spectrum.
     */

    extinction = cpl_image_duplicate(spectrum);
    map_table(extinction, startwave + dispersion/2, dispersion, 
              ext_table, "WAVE", "EXTINCTION");


    /*
     * Convert from magnitudes to actual flux loss.
     */

    cpl_image_multiply_scalar(extinction, 0.4 * airmass);
    cpl_image_exponential(extinction, 10.);


    /*
     * Correct the scientific spectrum to airmass 0
     */

    cpl_image_multiply(spectrum, extinction);


    /*
     * Find in what pixel interval (start at "ext_pos", for "ext_count" 
     * pixels) the atmospheric extinction is available.
     */
    
    ext_data = cpl_image_get_data_float(extinction);

    ext_count = 0;
    ext_pos = 0;
    for (i = 0; i < nx; i++) {
        if (ext_data[i] > 0.0) {
            if (ext_count == 0) {
                ext_pos = i;
            }
            ext_count++;
        }
        else {
            if (ext_count) {
                break;
            }
        }
    }

    cpl_image_delete(extinction); extinction = NULL;


    /*
     * Map the standard star catalog flux to the same lambda sampling
     * of the extracted spectrum.
     */

    flux = cpl_image_duplicate(spectrum);
    map_table(flux, startwave + dispersion/2, dispersion, 
              flux_table, "WAVE", "FLUX");


    /*
     * Find in what pixel interval (start at "flux_pos", for "flux_count" 
     * pixels) the standard star flux is available.
     */
    
    flux_data = cpl_image_get_data_float(flux);

    flux_count = 0;
    flux_pos = 0;
    for (i = 0; i < nx; i++) {
        if (flux_data[i] > 0.0) {
            if (flux_count == 0) {
                flux_pos = i;
            }
            flux_count++;
        }
        else {
            if (flux_count) {
                break;
            }
        }
    }


    /*
     * Intersection with previous selection
     */

    start      = ext_pos > flux_pos ? ext_pos : flux_pos;
    end        = (ext_pos + ext_count) < (flux_pos + flux_count) ?
                 (ext_pos + ext_count) : (flux_pos + flux_count);
    flux_pos   = start;
    flux_count = end - start;


    /*
     * Convert the flux to photons (per second per Angstrom).
     * std_flux is in units of erg / cm^2 / s / Angstrom. This
     * must be multiplied by the efficient area of the telescope,
     * 5.18E+5 cm^2, and divided by hv (v = frequency). With 
     * hc = 1.98E-8 erg*Angstrom one obtains the following:
     */

    physical = cpl_image_duplicate(spectrum);
    phys_data = cpl_image_get_data_float(physical);

    for (i = 0; i < nx; i++) {
        lambda = startwave + dispersion * (i + 0.5);
        phys_data[i] = 0.0026 * lambda * flux_data[i];
    }

    efficiency = cpl_image_duplicate(spectrum);
    eff_data = cpl_image_get_data_float(efficiency);
    data = cpl_image_get_data_float(spectrum);

    for (i = 0; i < nx; i++) {
        if (phys_data[i] > 0.0)
            eff_data[i] = data[i] / phys_data[i];
        else
            eff_data[i] = 0.0;
    }

    cpl_image_delete(physical); physical = NULL;


    /*
     * Find interval (longer than 300 pixels) where efficiency is 
     * greater than 1%
     */

    eff_count = 0;
    eff_pos = 0;
    for (i = 0; i < nx; i++) {
        if (eff_data[i] > 0.01) {
            if (eff_count == 0) {
                eff_pos = i; 
            }
            eff_count++;
        }
        else {
            if (eff_count > 300) {
                break;
            }
        }
    }


    /*
     * Intersection with previous selection
     */

    start      = eff_pos > flux_pos ? eff_pos : flux_pos;
    end        = (eff_pos + eff_count) < (flux_pos + flux_count) ?
                 (eff_pos + eff_count) : (flux_pos + flux_count);
    eff_pos    = start;
    eff_count  = end - start;

    if (eff_count < 1) {
        cpl_error_set_message(cpl_func, CPL_ERROR_INCOMPATIBLE_INPUT,
                              "No overlap between catalog and spectrum");
        cleanup;
        return NULL;
    }


    /*
     * Extract only data to fit, i.e., where the efficiency is available.
     */

    image = cpl_image_extract(efficiency, eff_pos + 1, 1, 
                              eff_pos + eff_count, 1);

    smo_image = polysmooth(image, order, 50);
    cpl_image_delete(image);

    smo_efficiency = cpl_image_duplicate(efficiency);
    smo_eff_data = cpl_image_get_data_float(smo_efficiency);
    cpl_image_copy(smo_efficiency, smo_image, eff_pos + 1, 1);

    cpl_image_delete(smo_image);


    /*
     * Compute instrument response as the ratio between the catalog
     * spectrum and the observed spectrum (converted in physical units).
     * The polynomial smoothing, however, is performed on the inverse
     * of this ration, for obvious reasons (i.e., no divergence at zero
     * efficiency).
     */

    response = cpl_image_duplicate(spectrum);
    res_data = cpl_image_get_data_float(response);

    for (i = 0; i < nx; i++) {
        if (eff_data[i] > 0.01 && flux_data[i] > 0.0)
            res_data[i] = data[i] / flux_data[i];
        else
            res_data[i] = 0.0;
    }


    /*
     * Extract only data to fit, i.e., where the response is available.
     */

    image = cpl_image_extract(response, eff_pos + 1, 1, eff_pos + eff_count, 1);

    smo_image = polysmooth(image, order, 50);
    cpl_image_delete(image);

    smo_response = cpl_image_duplicate(response);
    smo_res_data = cpl_image_get_data_float(smo_response);
    cpl_image_copy(smo_response, smo_image, eff_pos + 1, 1);

    cpl_image_delete(smo_image);

    for (i = 0; i < nx; i++) {
        if (eff_data[i] > 0.01) {
            res_data[i] = 1 / res_data[i];
            smo_res_data[i] = 1 / smo_res_data[i];
        }
        else {
            res_data[i] = 0.0;
            smo_res_data[i] = 0.0;
        }
    }


    /*
     * Assemble the product spectrophotometric table.
     */

    table = cpl_table_new(nx);

    cpl_table_new_column(table, "WAVE", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(table, "WAVE", "Angstrom");

    for (i = 0; i < nx; i++)
        cpl_table_set_float(table, "WAVE", i, startwave + dispersion*(i+0.5));

    cpl_table_new_column(table, "STD_FLUX", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(table, "STD_FLUX", 
                              "10^(-16) erg/(cm^2 s Angstrom)");
    cpl_table_copy_data_float(table, "STD_FLUX", flux_data);
    cpl_image_delete(flux); flux = NULL;

    cpl_table_new_column(table, "OBS_FLUX", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(table, "OBS_FLUX", "electron/(s Angstrom)");
    cpl_table_copy_data_float(table, "OBS_FLUX", data);
    cpl_image_delete(spectrum); spectrum = NULL;

    cpl_table_new_column(table, "RAW_EFFICIENCY", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(table, "RAW_EFFICIENCY", "electron/photon");
    cpl_table_copy_data_float(table, "RAW_EFFICIENCY", eff_data);
    cpl_image_delete(efficiency); efficiency = NULL;

    cpl_table_new_column(table, "EFFICIENCY", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(table, "EFFICIENCY", "electron/photon");
    cpl_table_copy_data_float(table, "EFFICIENCY", smo_eff_data);
    cpl_image_delete(smo_efficiency); smo_efficiency = NULL;

    cpl_table_new_column(table, "RAW_RESPONSE", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(table, "RAW_RESPONSE", 
                              "10^(-16) erg/(cm^2 electron)");
    cpl_table_copy_data_float(table, "RAW_RESPONSE", res_data);
    cpl_image_delete(response); response = NULL;

    cpl_table_new_column(table, "RESPONSE", CPL_TYPE_FLOAT);
    cpl_table_set_column_unit(table, 
                              "RESPONSE", "10^(-16) erg/(cm^2 electron)");
    cpl_table_copy_data_float(table, "RESPONSE", smo_res_data);
    cpl_image_delete(smo_response); smo_response = NULL;

    cleanup;

    return table;
}

static double ksigma_vector(cpl_vector *values, 
                            double klow, double khigh, int kiter, int *good)
{
    cpl_vector *accepted;
    double  mean  = 0.0;
    double  sigma = 0.0;
    double *data  = cpl_vector_get_data(values);
    int     n     = cpl_vector_get_size(values);
    int     ngood = n;
    int     count = 0;
    int     i;


    /*
     * At first iteration the mean is taken as the median, and the
     * standard deviation relative to this value is computed.
     */

    mean = cpl_vector_get_median(values);

    for (i = 0; i < n; i++) 
        sigma += (mean - data[i]) * (mean - data[i]);

    sigma = sqrt(sigma / (n - 1));

    while (kiter) {
        count = 0;
        for (i = 0; i < ngood; i++) {
            if (data[i]-mean < khigh*sigma && mean-data[i] < klow*sigma) {
                data[count] = data[i];
                ++count;
            }
        }

        if (count == 0) // This cannot happen at first iteration.
            break;      // So we can break: we have already computed a mean.

        /*
         * The mean must be computed even if no element was rejected
         * (count == ngood), because at first iteration median instead 
         * of mean was computed.
         */

        accepted = cpl_vector_wrap(count, data);
        mean = cpl_vector_get_mean(accepted);
        if (count > 1)
            sigma = cpl_vector_get_stdev(accepted);
        cpl_vector_unwrap(accepted);

        if (count == ngood || count == 1)
            break;

        ngood = count;
        --kiter;
    }

    if (good)
        *good = ngood;

    return mean;
}


/**
 * @brief
 *   Stack images using k-sigma clipping
 *
 * @param imlist      List of images to stack
 * @param klow        Number of sigmas for rejection of lowest values
 * @param khigh       Number of sigmas for rejection of highest values
 * @param kiter       Max number of iterations
 *
 * @return Stacked image.
 *
 * At the first iteration the value of sigma is computed relatively to 
 * the median value of all pixels at a given image position. For the 
 * next iterations the sigma is computed in the standard way. If 
 * at some iteration all points would be rejected, the mean computed
 * at the previous iteration is returned.
 */

cpl_image *mos_ksigma_stack(cpl_imagelist *imlist, 
                            double klow, double khigh, int kiter,
                            cpl_image **good)
{
    int         ni, nx, ny, npix;
    cpl_image  *out_ima;
    float      *pout_ima;
    float      *good_ima;
    cpl_image  *image;
    float     **data;
    cpl_vector *time_line;
    double     *ptime_line;
    int         ngood;
    int         i, j;


    ni         = cpl_imagelist_get_size(imlist);

    image      = cpl_imagelist_get(imlist, 0);
    nx         = cpl_image_get_size_x(image);
    ny         = cpl_image_get_size_y(image);
    npix       = nx * ny;
    
    out_ima    = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
    pout_ima   = cpl_image_get_data_float(out_ima);

    if (good) {
        *good = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
        good_ima = cpl_image_get_data_float(*good);
    }

    time_line  = cpl_vector_new(ni);
    ptime_line = cpl_vector_get_data(time_line);

    data = cpl_calloc(sizeof(float *), ni);
    
    for (i = 0; i < ni; i++) {
        image = cpl_imagelist_get(imlist, i);
        data[i] = cpl_image_get_data_float(image);
    }

    for (i = 0; i < npix; i++) {
        for (j = 0; j < ni; j++) {
            ptime_line[j] = data[j][i];
        }
        pout_ima[i] = ksigma_vector(time_line, klow, khigh, kiter, &ngood);
        if (good) {
            good_ima[i] = ngood;
        }
    }

    cpl_free(data);
    cpl_vector_delete(time_line);

    return out_ima;

}


/**
 * @brief
 *   Apply response curve to extracted spectra.
 *
 * @param spectra     Image containing extracted spectra (slits or objects)
 * @param response    Table including the response curve
 * @param ext_table   Atmospheric extinction table
 * @param startwave   Start wavelength
 * @param dispersion  Angstrom per pixel of input image
 * @param gain        Gain factor (e-/ADU)
 * @param exptime     Exposure time (seconds)
 * @param airmass     Airmass of observation
 *
 * The response function is resampled, trough interpolation, to the target
 * spectrum resolution. For the wavelengths outside of the response range,
 * the target spectrum is set to -1. 
 *
 * @return Photometrically calibrated image at airmass 0.
 */

cpl_image *mos_apply_photometry(cpl_image *spectra, cpl_table *response,
                                cpl_table *ext_table, double startwave,
                                double dispersion, double gain,
                                double exptime, double airmass)
{
    cpl_image *extinction;
    cpl_image *outspectra;
    cpl_image *mapresponse;
    float     *res_data;
    float     *out_data;
    float     *ext_data;
    int        tlength, xlength, ylength;
    int        i, j, k;
    double     resp_startwave;
    double     resp_endwave;
    int        null;


    if (spectra == NULL || ext_table == NULL || response == NULL) {
        cpl_error_set(cpl_func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }
    
    /* Use the normal response if available. If not, use the flat corrected
     * response. Usually only one of the two is available in the interpolated
     * response.
     */
    if(cpl_table_has_column(response, "RESPONSE"))
        cpl_table_cast_column(response, "RESPONSE", "RESPONSE_F", CPL_TYPE_FLOAT);
    else if(cpl_table_has_column(response, "RESPONSE_FFSED"))
        cpl_table_cast_column(response, "RESPONSE_FFSED", "RESPONSE_F", CPL_TYPE_FLOAT);
    else
        return NULL;
    
    res_data = cpl_table_get_data_float(response, "RESPONSE_F");

    if (res_data == NULL) {
        cpl_error_set(cpl_func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    tlength = cpl_table_get_nrow(response);
    xlength = cpl_image_get_size_x(spectra);
    ylength = cpl_image_get_size_y(spectra);

    /* Map the response */
    mapresponse = cpl_image_new(xlength, 1, CPL_TYPE_FLOAT);
    map_table(mapresponse, startwave + dispersion/2, dispersion,
              response, "WAVE", "RESPONSE_F");
    res_data = cpl_image_get_data_float(mapresponse);

    /*
     * Map the atmospheric extinction factors to the same lambda sampling
     * of the extracted spectrum.
     */

    extinction = cpl_image_new(xlength, 1, CPL_TYPE_FLOAT);
    map_table(extinction, startwave + dispersion/2, dispersion,
              ext_table, "WAVE", "EXTINCTION");


    /*
     * Convert from magnitudes to actual flux loss.
     */

    cpl_image_multiply_scalar(extinction, 0.4 * airmass);
    cpl_image_exponential(extinction, 10.);

    outspectra = cpl_image_duplicate(spectra);

    ext_data = cpl_image_get_data_float(extinction);
    out_data = cpl_image_get_data_float(outspectra);

    for (k = 0, i = 0; i < ylength; i++) {
        for (j = 0; j < xlength; j++, k++)
            out_data[k] *= ext_data[j] * res_data[j];
    }

    cpl_image_delete(extinction);
    cpl_image_delete(mapresponse);

    cpl_image_multiply_scalar(outspectra, gain / exptime / dispersion);

    /* 
     * Set to -1 the extrapolated values
     */
    resp_startwave = cpl_table_get(response, "WAVE", 0, &null);
    resp_endwave = cpl_table_get(response, "WAVE", 
                                 cpl_table_get_nrow(response) -1, &null);
    for (j = 0; j < xlength; j++) {
        double this_wave = startwave + j * dispersion;
        if(this_wave < resp_startwave ||this_wave > resp_endwave)
        {
            for (i = 0; i < ylength; i++) 
                out_data[j + xlength * i] = -1;
        }
    }
    
    cpl_table_erase_column(response, "RESPONSE_F");

    return outspectra;
}


/**
 * @brief
 *   Propagate errors from response curve and extracted spectra
 *
 * @param errors      Image containing errors on extracted spectra (objects)
 * @param response    Table including the response curve errors
 * @param ext_table   Atmospheric extinction table
 * @param startwave   Start wavelength
 * @param dispersion  Angstrom per pixel of input image
 * @param gain        Gain factor (e-/ADU)
 * @param exptime     Exposure time (seconds)
 * @param airmass     Airmass of observation
 *
 * @return Error on photometrically calibrated spectrum at airmass 0.
 */

cpl_image *mos_propagate_photometry_error(cpl_image *spectra, 
                                          cpl_image *errors, 
                                          cpl_table *response,
                                          cpl_table *ext_table, 
                                          double startwave,
                                          double dispersion, double gain,
                                          double exptime, double airmass)
{
    cpl_image *extinction;
    cpl_image *outerrors;
    cpl_image *mapresponse;
    cpl_image *maperror;
    float     *err_data;
    float     *out_data;
    float     *ext_data;
    float     *res_data;
    float     *spe_data;
    int        tlength, xlength, ylength;
    int        i, j, k;


    if (errors == NULL || ext_table == NULL || response == NULL) {
        cpl_error_set(cpl_func, CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (!cpl_table_has_column(response, "ERROR")) {
        return mos_apply_photometry(errors, response, ext_table, startwave,
                                    dispersion, gain, exptime, airmass);
    }

    cpl_table_cast_column(response, "RESPONSE", "RESPONSE_F", CPL_TYPE_FLOAT);
    res_data = cpl_table_get_data_float(response, "RESPONSE_F");

    if (res_data == NULL) {
        cpl_error_set(cpl_func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    err_data = cpl_table_get_data_float(response, "ERROR");

    if (err_data == NULL) {
        cpl_error_set(cpl_func, CPL_ERROR_DATA_NOT_FOUND);
        return NULL;
    }

    tlength = cpl_table_get_nrow(response);
    xlength = cpl_image_get_size_x(errors);
    ylength = cpl_image_get_size_y(errors);

    if (xlength != tlength) {
        mapresponse = cpl_image_new(xlength, 1, CPL_TYPE_FLOAT);
        map_table(mapresponse, startwave + dispersion/2, dispersion,
                  response, "WAVE", "RESPONSE_F");
        res_data = cpl_image_get_data_float(mapresponse);

        maperror = cpl_image_new(xlength, 1, CPL_TYPE_FLOAT);
        map_table(maperror, startwave + dispersion/2, dispersion,
                  response, "WAVE", "ERROR");
        err_data = cpl_image_get_data_float(maperror);
    }

    /*
     * Map the atmospheric extinction factors to the same lambda sampling
     * of the extracted spectrum.
     */

    extinction = cpl_image_new(xlength, 1, CPL_TYPE_FLOAT);
    map_table(extinction, startwave + dispersion/2, dispersion,
              ext_table, "WAVE", "EXTINCTION");


    /*
     * Convert from magnitudes to actual flux loss.
     */

    cpl_image_multiply_scalar(extinction, 0.4 * airmass);
    cpl_image_exponential(extinction, 10.);

    outerrors = cpl_image_duplicate(errors);

    ext_data = cpl_image_get_data_float(extinction);
    out_data = cpl_image_get_data_float(outerrors);
    spe_data = cpl_image_get_data_float(spectra);

    for (k = 0, i = 0; i < ylength; i++) {
        for (j = 0; j < xlength; j++, k++) {
            out_data[k] = ext_data[j] * 
              sqrt(err_data[j] * err_data[j] * spe_data[k] * spe_data[k] +
                   res_data[j] * res_data[j] * out_data[k] * out_data[k]);
        }
    }

    cpl_image_delete(extinction);
    if (xlength != tlength) {
        cpl_image_delete(maperror);
    }

    cpl_image_multiply_scalar(outerrors, gain / exptime / dispersion);

    cpl_table_erase_column(response, "RESPONSE_F");
    return outerrors;
}


/**
 * @brief
 *   Estimate linear polarisation parameters on spectral interval
 *
 * @param q_image      Image whose rows correspond to extracted spectra 
 * @param q_error      Image with errors of @em q_image
 * @param u_image      Image whose rows correspond to extracted spectra 
 * @param u_error      Image with errors of @em u_image
 * @param startwave    Start wavelength
 * @param dispersion   Angstrom per pixel of input image
 * @param band         Width of band where the signal is averaged (Angstrom)
 * @param pol_sta      Table with polarimetric standard stars parameters
 * @param ra           Right Ascension of center of field of view (degrees)
 * @param dec          Declination of center of field of view (degrees)
 * @param filter       Returned filter (UBVRI) to which the measurement refers
 * @param polarisation Returned polarisation flag (0 = star is not polarised)
 * @param p_offset     Returned relative variation of P in @em filter
 * @param p_error      Returned error on @em p_offset
 * @param a_offset     Returned variation of angle in @em filter (degrees)
 * @param a_error      Returned error on @em a_offset (degrees)
 *
 * @return 0 on success.
 *
 * This function is used to estimate the deviation of the linear 
 * polarization of observed polarimetric standard stars from 
 * expected (catalog) values.
 *
 * The parameter @em p_offset is defined as observed minus expected
 * polarisation, divided by expected polarisation. Only in case the
 * expected polarisation is zero, @em p_offset is simply defined as
 * the observed polarisation.
 *
 * The parameter @em p_error is obtained by propagating the measurement 
 * and the catalog errors.
 *
 * The parameter @em a_offset is defined as observed minus expected
 * polarisation angle. In case the expected polarisation is zero, 
 * @em a_offset is simply set to zero.
 *
 * The parameter @em a_error is obtained by propagating the measurement 
 * and the catalog errors. In case the expected polarisation is zero, 
 * @em a_error is simply set to zero.
 *
 * The applied algorithm is the following:
 *
 * Each row of the input images corresponds to a different object.
 * The row corresponding to the standard star is selected as the one
 * where the S/N ratio is the greatest (this check is performed
 * conventionally on @em q_error, which is minimal for high S/N
 * objects). The selected row is extracted from all input images.
 *
 * From the table @em pol_sta the standard star closest to the specified
 * @em ra and @em dec coordinates is found, and if this is closer than
 * half a degree from the input @em ra and @em dec is selected. From the 
 * available expected values at the different spectral bands (U, B, V, 
 * etc.) the one which is closest to the center of the valid spectral 
 * interval for the standard star observation is selected (the valid 
 * spectral interval is defined here as the region where the error is 
 * greater than zero). The following wavelengths are associated to each 
 * band:
 *
 * U -> 3650 Angstrom
 * B -> 4450 Angstrom
 * V -> 5510 Angstrom
 * R -> 6580 Angstrom
 * I -> 8060 Angstrom
 *
 * If the interval @em band centered at the wavelength closest to center
 * of the valid spectral range is not entirely contained in the valid
 * spectral range, this function fails and 1 is returned. Otherwise,
 * within the defined band, the median values of Q and U, together 
 * with their median errors, are determined. From those the values 
 * for the output parameters are computed in the usual way.
 */

int mos_check_polarisation(cpl_image *q_image, cpl_image *q_error,
                           cpl_image *u_image, cpl_image *u_error,
                           double startwave, double dispersion,
                           double band, cpl_table *pol_sta,
                           double ra, double dec, char *filter, 
                           int *polarisation,
                           double *p_offset, double *p_error,
                           double *a_offset, double *a_error)
{
    cpl_table *standard;
    cpl_image *q_noise;
    cpl_image *q_signal;
    cpl_image *u_noise;
    cpl_image *u_signal;
    cpl_image *noise;
    double    *q_ndata;
    double     arctol = 0.5;  /* Arc tolerance in degrees */
    double     mindist;
    double     cwave;
    double     bwave[] = {3650., 4450., 5510., 6580., 8060};
    char      *bands = "UBVRI";
    char       p_label[] = {' ', 'p', '\0'};
    char       dp_label[] = {' ', 'd', 'p', '\0'};
    char       a_label[] = {' ', 'a', '\0'};
    char       da_label[] = {' ', 'd', 'a', '\0'};
    int        nbands = strlen(bands);
    int        selected;
    int        first, last, count, center;
    int        nx;
    cpl_size   col, row;
    int        i, found, closest;
    int        pband;
    int        polarised;
    double     q_obs;
    double     q_err;
    double     u_obs;
    double     u_err;
    double     p_obs;
    double     p_err;
    double     p_ref;
    double     dp_ref;
    double     a_obs;
    double     a_err;
    double     a_ref;
    double     da_ref;


    *filter       = '\0';
    *polarisation = 0;
    *p_offset     = 0.0;
    *p_error      = 0.0;
    *a_offset     = 0.0;
    *a_error      = 0.0;

    /*
     * Select reference standard star
     */

    cpl_table_select_all(pol_sta);
    cpl_table_and_selected_double(pol_sta, "RA",  CPL_GREATER_THAN, ra-arctol);
    cpl_table_and_selected_double(pol_sta, "RA",  CPL_LESS_THAN,    ra+arctol);
    cpl_table_and_selected_double(pol_sta, "DEC", CPL_GREATER_THAN, dec-arctol);
    selected =
    cpl_table_and_selected_double(pol_sta, "DEC", CPL_LESS_THAN,    dec+arctol);

    if (selected == 0) {
        cpl_msg_warning(cpl_func, "No standard star found in FOV");
        return 1;
    }

    if (selected > 1) {
        cpl_msg_warning(cpl_func, 
                        "Ambiguity: %d standard stars found in FOV", selected);
        return 1;
    }

    standard = cpl_table_extract_selected(pol_sta);

    cpl_msg_info(cpl_func, "Standard star: %s", 
                 cpl_table_get_string(standard, "name", 0));

    /*
     * Check whether the star is polarised or not
     */

    polarised = cpl_table_get_int(standard,  "polarised", 0, NULL);

    cpl_msg_info(cpl_func, "This star is%sexpected to be polarised",
                 polarised ? " " : " not ");


    /*
     * Determine the image row with the smallest median noise: this 
     * row is assumed to refer to the standard star.
     * (note: the higher the S/N ratio of the original spectra, the 
     * smaller the noise of the Stokes parameters Q and U).
     */

    nx = cpl_image_get_size_x(q_error);

    noise = cpl_image_collapse_median_create(q_error, 1, 0, 0);
    cpl_image_get_minpos(noise, &col, &row);

    cpl_image_delete(noise);

    if (col != 1) {
        cpl_table_delete(standard);
        cpl_msg_error(cpl_func, 
                      "Assertion failure!!! col = %"CPL_SIZE_FORMAT" (it should be 1)", col);
        return 1;
    }

    q_signal = cpl_image_extract(q_image, 1, row, nx, row);
    q_noise  = cpl_image_extract(q_error, 1, row, nx, row);
    u_signal = cpl_image_extract(u_image, 1, row, nx, row);
    u_noise  = cpl_image_extract(u_error, 1, row, nx, row);

    q_ndata = cpl_image_get_data_double(q_noise);


    /*
     * Determine valid interval in input images (where error is positive).
     */

    first = -1;
    last = nx = cpl_image_get_size_x(q_signal);
    for (i = 0; i < nx; i++) {
        if (first < 0) {
            if (q_ndata[i] > 0.0) {
                first = i;
            }
        }
        else {
            if (q_ndata[i] <= 0.0) {
                last = i - 1;
                break;
            }
        }
    }

    count = last - first + 1;

    if (first < 0 || count < band) {
        cpl_table_delete(standard);
        cpl_image_delete(q_signal);
        cpl_image_delete(q_noise);
        cpl_image_delete(u_signal);
        cpl_image_delete(u_noise);
        cpl_msg_warning(cpl_func, "Too short spectrum (%d pixels)", count);
        return 1;
    }

    center = (first + last) / 2;              // Center of valid spectrum
    cwave = startwave + dispersion * center;  // Corresponding wavelength


    /*
     * Find the band UBVRI closest to the central wavelength.
     */

    found = 0;
    for (i = 0; i < nbands; i++) {
        p_label[0] = bands[i];
        if (cpl_table_is_valid(standard, p_label, 0)) {
            if (found == 0) {
                found = 1;
                mindist = fabs(bwave[i] - cwave);
                closest = i;
            }
            else if (mindist > fabs(bwave[i] - cwave)) {
                mindist = fabs(bwave[i] - cwave);
                closest = i;
            }
        }
    }

    if (!found) {
        cpl_table_delete(standard);
        cpl_image_delete(q_signal);
        cpl_image_delete(q_noise);
        cpl_image_delete(u_signal);
        cpl_image_delete(u_noise);
        cpl_msg_warning(cpl_func, "No reference value available");
        return 1;
    }

    center = (bwave[closest] - startwave) / dispersion; // Center of band (pix)
    cwave  =  bwave[closest];                           // Wavelength of band


    /*
     * Check that the integration interval is entirely contained
     * in the valid interval, or give it up.
     */

    pband = floor(band / dispersion);  // Band width in pixels

    if (center - pband/2 < first || center + pband/2 > last) {
        cpl_table_delete(standard);
        cpl_image_delete(q_signal);
        cpl_image_delete(q_noise);
        cpl_image_delete(u_signal);
        cpl_image_delete(u_noise);
        cpl_msg_warning(cpl_func, "No reference value available");
        return 1;
    }

    first = center - pband/2;
    last  = center + pband/2;

    /*
     * Collect reference values. Note that if angle info is not available,
     * angle stuff is set automaticaly to zero.
     */

     p_label[0] = bands[closest];
    dp_label[0] = bands[closest];
     a_label[0] = bands[closest];
    da_label[0] = bands[closest];

     p_ref = cpl_table_get(standard,  p_label, 0, NULL);
    dp_ref = cpl_table_get(standard, dp_label, 0, NULL);
     a_ref = cpl_table_get(standard,  a_label, 0, NULL);
    da_ref = cpl_table_get(standard, da_label, 0, NULL);

    cpl_msg_info(cpl_func, 
                 "The expected polarisation is %.2f +- %.2f %%", 
                 p_ref, dp_ref);

    if (polarised) {
        cpl_msg_info(cpl_func, 
                     "The expected polarisation angle is %.2f +- %.2f degrees", 
                     a_ref, da_ref);
    }

    /*
     * Find median signal and median error.
     */

    q_obs = cpl_image_get_median_window(q_image, first, 1, last, 1);
    q_err = cpl_image_get_median_window(q_error, first, 1, last, 1);
    u_obs = cpl_image_get_median_window(u_image, first, 1, last, 1);
    u_err = cpl_image_get_median_window(u_error, first, 1, last, 1);

    /*
     * Measured linear polarisation and its error
     */

    p_obs = sqrt(q_obs * q_obs + u_obs * u_obs);
    p_err = CPL_MATH_SQRT1_2 * 0.5 * (q_err + u_err);

    /*
     * Measured polarisation angle
     */

    a_obs = 0.0;
    if (polarised) {
        if (fabs(q_obs) < 0.00001) {
            if (u_obs > 0.0) {
                a_obs = 45.0;
            }
            else {
                a_obs = 135.0;
            }
        }
        else {
            a_obs = 0.5 * atan(u_obs / q_obs) * 180 / CPL_MATH_PI;
            if (q_obs > 0.0) {
                if (u_obs < 0.0) {
                    a_obs += 180.;
                }
            }
            else {
                a_obs += 90.;
            }
        }
    }

    /*
     * Error on polarisation angle
     */

    a_err = 0.0;
    if (polarised) {
        a_err = sqrt(q_obs*q_obs*u_err*u_err + u_obs*u_obs*q_err*q_err)
              / (p_obs * p_obs) 
              * 90 / CPL_MATH_PI;
    }

    p_obs *= 100;
    p_err *= 100;
    cpl_msg_info(cpl_func, 
                 "The measured polarisation is %.2f +- %.2f %%", 
                 p_obs, p_err);

    if (polarised) {
        cpl_msg_info(cpl_func, 
                     "The measured polarisation angle is %.2f +- %.2f degrees", 
                     a_obs, a_err);
    }

    *filter       = bands[closest];
    *polarisation = polarised;

    if (polarised) {
        *p_offset = (p_obs - p_ref) / p_ref;
        *p_error  = sqrt(p_err * p_err + dp_ref * dp_ref) / p_ref;
    }
    else {
        *p_offset = p_obs - p_ref;
        *p_error  = sqrt(p_err * p_err + dp_ref * dp_ref);
    }

    *a_offset     = a_obs - a_ref;
    *a_error      = sqrt(a_err*a_err + da_ref*da_ref);

    return 0;

}


/**
 * @brief
 *   Estimate offset between two object tables.
 *
 * @param reference    Reference object table.
 * @param objects      Object table from offset frame.
 * @param offset       Returned offset in CCD pixels.
 *
 * @return CPL_ERROR_NONE in case of success.
 *
 * Given two object tables, derived from two scientific exposures obtained
 * with the same mask, this function determines the offset between the
 * two tables.
 *
 * The procedure is the following:
 * For each slit (corresponding to one row of the object tables), two 
 * integer arrays of length "length" (read from the object table) are 
 * allocated. The position of objects in that slit, as given in both 
 * the reference and the offset object tables, are flagged with 1, all 
 * the rest is left to zero. The two arrays are correlated, finding
 * a preliminary (integer) offset in pixel. This preliminary offset
 * is just used to match corresponding objects. At a second step, the
 * median offset between matching objects is computed. This offset is
 * converted to CCD pixels according to the
 * @code
 * CCD_offset = Map_offset * (t - b) / length
 * @endcode
 * in the same convention used in function @c mos_spatial_calibration().
 * The returned offset is the median offset obtained from all slits.
 */

int mos_compute_offset(cpl_table *reference, cpl_table *objects, double *offset)
{
    cpl_array *offsets;
    int        noffset;
    int        nslits = cpl_table_get_nrow(reference);
    int       *nref;
    int       *nobj;
    int        corr, maxcorr;
    int        best_shift;
    int        i, j, k;

    cpl_error_code status = CPL_ERROR_NONE;


    *offset = 0.0;

    if (objects == NULL)
        return CPL_ERROR_NULL_INPUT;

    if (nslits != cpl_table_get_nrow(objects))
        return CPL_ERROR_INCOMPATIBLE_INPUT;

    nref = fors_get_nobjs_perslit(reference);
    nobj = fors_get_nobjs_perslit(objects);

    noffset = 0;
    for (i = 0; i < nslits; i++)
        noffset += nobj[i];

    if (noffset == 0) {
        cpl_free(nref);
        cpl_free(nobj);
        return CPL_ERROR_DATA_NOT_FOUND;
    }

    noffset = 0;
    for (i = 0; i < nslits; i++)
        noffset += nref[i];

    if (noffset == 0) {
        cpl_free(nref);
        cpl_free(nobj);
        return CPL_ERROR_DATA_NOT_FOUND;
    }

    offsets = cpl_array_new(noffset, CPL_TYPE_DOUBLE);

    noffset = 0;    // The real number of offsets found will be counted.

    for (i = 0; i < nslits; i++) {
        if (nref[i] > 0 && nobj[i] > 0) {
            double shift;
            int    length  = cpl_table_get_int(objects, "length", i, NULL);
            double ytop    = cpl_table_get_double(objects, "xtop", i, NULL);
            double ybottom = cpl_table_get_double(objects, "xbottom", i, NULL);
            int   *aref    = cpl_calloc(length, sizeof(int));
            int   *aobj    = cpl_calloc(length, sizeof(int));
            float *pref    = cpl_calloc(nref[i], sizeof(float));
            float *pobj    = cpl_calloc(nobj[i], sizeof(float));
            
            for (j = 0; j < nref[i]; j++) {
                pref[j] = fors_get_object_position(reference, i, j + 1);
                aref[(int)pref[j]] = 1;
            }

            for (j = 0; j < nobj[i]; j++) {
                pobj[j] = fors_get_object_position(objects, i, j + 1);
                aobj[(int)pobj[j]] = 1;
            }

            /*
             * Do not consider objects at border
             */

            aref[0] = 0;
            aref[length - 1] = 0;
            aobj[0] = 0;
            aobj[length - 1] = 0;

//for (j = 0; j < nref[i]; j++)
//printf("references: %f, ", pref[j]);
//printf("\n");
//for (j = 0; j < nref[i]; j++)
//printf("objects   : %f, ", pobj[j]);
//printf("\n");
//for (j = 0; j < length; j++)
//printf("%d", aref[j]);
//printf("\n");
//for (j = 0; j < length; j++)
//printf("%d", aobj[j]);
//printf("\n");

            /*
             * Object matching by correlation
             */

            maxcorr = 0;
            best_shift = length;

            for (shift = length/2, j = 0; j <= length; shift--, j++) {
                int rstart, ostart, count;

                if (shift > 0) {
                   rstart = shift;
                   ostart = 0;
                   count  = length - shift;
                }
                else {
                   rstart = 0;
                   ostart = -shift;
                   count  = length + shift;
                }

                corr = 0;
                for (k = 0; k < count; k++) {
                    corr += aref[rstart + k] * aobj[ostart + k];
                }

                if (maxcorr < corr) {
                    maxcorr = corr;
                    best_shift = shift;
                }
            }

            if (best_shift == length) { // No shift found
//printf("%d: No shift found\n", i);
                cpl_free(aref);
                cpl_free(aobj);
                cpl_free(pref);
                cpl_free(pobj);
                continue;
            }
//printf("%d: Integer shift found = %d\n", i, best_shift);

            for (j = 0; j < nref[i]; j++) {
                for (k = 0; k < nobj[i]; k++) {
                    if (fabs(pref[j] - pobj[k] - best_shift) < 2) {
                       double ccd_offset = (pref[j] - pobj[k]) 
                                         * (ytop - ybottom)
                                         / length;

//printf("%d: Match found: %f\n", i, ccd_offset);
                       /* 
                        * The matching object is found, store the
                        * corresponding offset
                        */

                       cpl_array_set(offsets, noffset, ccd_offset);
                       noffset++;
                       break;
                    }
                }
            }

            cpl_free(aref);
            cpl_free(aobj);
            cpl_free(pref);
            cpl_free(pobj);
        }
//else
//printf("%d: No object found\n", i);
    }

    cpl_free(nref);
    cpl_free(nobj);

//printf("%d offsets found in total\n", noffset);
    if (noffset > 0) {
        if (noffset % 2) {
            *offset = cpl_array_get_median(offsets);
        }
        else {
            double *a = cpl_malloc(sizeof(double) * noffset);
            for (i = 0; i < noffset; i++) {
                a[i] = cpl_array_get_double(offsets, i, NULL);
            }
            *offset = (fors_tools_get_kth_double(a, noffset, (noffset-1)/2) +
                       fors_tools_get_kth_double(a, noffset, (noffset/2))) / 2.0;
            cpl_free(a);
        }
    }
    else
        status = CPL_ERROR_DATA_NOT_FOUND;
//printf("Median offset: %f\n", *offset);

    cpl_array_delete(offsets);

    return status;

}


/**
 * @brief
 *   Shift values in an image
 *
 * @param image  Input image
 * @param dx     Shift in x
 * @param dy     Shift in y
 *
 * @return @c CPL_ERROR_NONE on success.
 */

cpl_error_code mos_image_shift(cpl_image *image, double dx, double dy)
{
    cpl_image *source;
    int        nx = cpl_image_get_size_x(image);
    int        ny = cpl_image_get_size_y(image);
    float     *idata;
    float     *sdata;
    int        i, j, pos;
    double     xpos, ypos, xfrac, yfrac;
    int        xint, yint;


    if (fabs(dx) >= nx || fabs(dy) >= ny)
        return CPL_ERROR_ACCESS_OUT_OF_RANGE;

    source = cpl_image_duplicate(image);
    idata = cpl_image_get_data_float(image);
    sdata = cpl_image_get_data_float(source);

    /*
     * Shift in y
     */

    yfrac = - dy - floor(- dy);
    xfrac = - dx - floor(- dx);

    for (pos = 0, j = 0; j < ny; j++) {
        ypos = j - dy;
        yint = floor(ypos);
        for (i = 0; i < nx; i++) {
            xpos = i - dx;
            xint = floor(xpos);
            if (xint < 0 || yint < 0 || xint > nx - 2 || yint > ny - 2) {
                idata[pos] = 0.0;
            }
            else {
                idata[pos] = sdata[xint + nx*yint] * (1 - xfrac) * (1 - yfrac)
                           + sdata[xint + 1 + nx*yint] * xfrac * (1 - yfrac)
                           + sdata[xint + nx*(yint + 1)] * (1 - xfrac) * yfrac
                           + sdata[xint + 1 + nx*(yint + 1)] * xfrac * yfrac;
            }
            pos++;
        }
    }

    cpl_image_delete(source);

    return CPL_ERROR_NONE;
}

/**
 * @brief
 *   Return slit closest to CCD center.
 *
 * @param slits  Slit table
 * @param nx     X size of CCD
 * @param ny     Y size of CCD
 *
 * @return Table row with the slit closest to center.
 */

int mos_slit_closest_to_center(cpl_table *slits, int nx, int ny)
{
#ifdef CPL_SIZE_FORMAT
    cpl_size row;
#else
    int row;
#endif

    cpl_table_duplicate_column(slits, "x", slits, "xtop");
    cpl_table_add_columns(slits, "x", "xbottom");
    cpl_table_divide_scalar(slits, "x", 2);         // Mean x position
    cpl_table_subtract_scalar(slits, "x", nx/2);    // Relative to CCD center
    cpl_table_multiply_columns(slits, "x", "x");    // Squared

    cpl_table_duplicate_column(slits, "y", slits, "ytop");
    cpl_table_add_columns(slits, "y", "ybottom");
    cpl_table_divide_scalar(slits, "y", 2);         // Mean y position
    cpl_table_subtract_scalar(slits, "y", ny/2);    // Relative to CCD center
    cpl_table_multiply_columns(slits, "y", "y");    // Squared

    cpl_table_add_columns(slits, "x", "y");         // Distance from center
    cpl_table_get_column_minpos(slits, "x", &row);  // Min distance from center

    cpl_table_erase_column(slits, "x");
    cpl_table_erase_column(slits, "y");

    return (int)row;
}

/**
 * @brief
 *   Measure flux from spectral interval on CCD
 *
 * @param image  Image containing raw spectra.
 * @param slits  Table with slits properties.
 * @param dx     Pixels to extract along the dispersion direction (radius).
 * @param gain   In electrons/ADU, used for error computation.
 * @param o_flux Returned integrated flux, in ADU/mm^2
 * @param o_err  Returned error on integrated flux.
 *
 * @return Status
 *
 * This one should integrate counts on a rectangle around a given 
 * wavelength on a spectrum corresponding to a given slit. Then
 * the counts are normalized to the corresponding physical area on
 * the slit.
 */

cpl_error_code mos_extract_flux(cpl_image *image, cpl_table *slits, 
                     double xwidth, double ywidth,
                     int dx, double gain, double *o_flux, double *o_err)
{
    int    nx      = cpl_image_get_size_x(image);
    int    ny      = cpl_image_get_size_y(image);
    int    slit    = mos_slit_closest_to_center(slits, nx, ny);
    int    ytop    = (int)cpl_table_get(slits, "ytop", slit, NULL);
    int    ybottom = (int)cpl_table_get(slits, "ybottom", slit, NULL);
    int    dy      = ytop - ybottom;
    int    xcenter = (int)((cpl_table_get(slits, "xtop", slit, NULL) +
                            cpl_table_get(slits, "xbottom", slit, NULL)) / 2);
    int    xleft   = xcenter - dx;
    int    xright  = xcenter + dx + 1;
    double area    = xwidth * ywidth; // squared mm
    int    npix    = (2*dx + 1) * dy;
    int    count   = 0;
    float *data    = cpl_image_get_data_float(image);
    double flux    = 0.0;
    double error   = 0.0;
    int    satur   = 60000;
    int    x, y;


    if (cpl_table_has_column(slits, "ywidth")) {
        area    = cpl_table_get(slits, "xwidth", slit, NULL)
                * cpl_table_get(slits, "ywidth", slit, NULL);
    }

    *o_flux = 0.0;
    *o_err = 0.0;

    if (xleft < 0)
        xleft = 0;

    if (xleft > nx)
        xleft = nx;

    if (xright < 0)
        xright = 0;

    if (xright > nx)
        xright = nx;

    if (ytop < 0)
        ytop = 0;

    if (ytop > ny)
        ytop = ny;

    if (ybottom < 0)
        ybottom = 0;

    if (ybottom > ny)
        ybottom = ny;

    count = (xright - xleft) * (ytop - ybottom);

    if (count == 0)
        return CPL_ERROR_ACCESS_OUT_OF_RANGE;

    count = 0;

    for (y = ybottom; y < ytop; y++) {
        for (x = xleft; x < xright; x++) {
            double value = data[x + y * nx];
            if (value < satur) {
                flux += value;
                count++;
            }
        }
    }

    if (count == 0)
        return CPL_ERROR_DIVISION_BY_ZERO;

    error = sqrt(flux/gain);

    /*
     * Flux correction for lost pixels
     */

    flux *= (float)npix / count;
    error *= (float)npix / count;

    flux /= area;
    error /= area;

    *o_flux = flux;
    *o_err = error;

    return CPL_ERROR_NONE;
}


/**
 * @brief
 *   Measure flux from spectral interval on remapped frame
 *
 * @param image       Image containing remapped spectra.
 * @param slits       Table with slits properties.
 * @param lambda      Wavelength to examine
 * @param startwave   Shortest wavelength in image.
 * @param dispersion  Wavelength units per image pixel
 * @param dx          Pixels to extract along the dispersion direction (radius).
 * @param gain        In electrons/ADU, used for error computation.
 * @param o_flux      Returned integrated flux, in ADU/mm^2
 * @param o_err       Returned error on integrated flux.
 *
 * @return Status
 *
 * This one should integrate counts on a rectangle around a given 
 * wavelength on a spectrum corresponding to the slit closest to
 * the CCD center. Then the counts are normalized to the corresponding 
 * physical area on the slit.
 */

cpl_error_code mos_extract_flux_mapped(cpl_image *image, cpl_table *slits,
                                       double xwidth, double ywidth,
                                       double lambda, double startwave, 
                                       double dispersion, int dx, double gain, 
                                       double *o_flux, double *o_err)
{
    int    nx      = cpl_image_get_size_x(image);
    int    ny      = cpl_image_get_size_y(image);
    int    slit    = mos_slit_closest_to_center(slits, nx, ny);
    int    dy      = (int)cpl_table_get(slits, "length", slit, NULL);
    int    ybottom = (int)cpl_table_get(slits, "position", slit, NULL);
    int    ytop    = ybottom + dy;
    int    xcenter = (int)floor((lambda - startwave) / dispersion + 0.5);
    int    xleft   = xcenter - dx;
    int    xright  = xcenter + dx + 1;
    double area    = xwidth * ywidth;
    int    npix    = (2*dx + 1) * dy;
    int    count   = 0;
    float *data    = cpl_image_get_data_float(image);
    double flux    = 0.0;
    double error   = 0.0;
    int    satur   = 60000;
    int    x, y;


    if (cpl_table_has_column(slits, "ywidth")) {
        area    = cpl_table_get(slits, "xwidth", slit, NULL)
                * cpl_table_get(slits, "ywidth", slit, NULL);
    }

    *o_flux = 0.0;
    *o_err = 0.0;

    if (xleft < 0)
        xleft = 0;

    if (xleft > nx)
        xleft = nx;

    if (xright < 0)
        xright = 0;

    if (xright > nx)
        xright = nx;

    if (ytop < 0)
        ytop = 0;

    if (ytop > ny)
        ytop = ny;

    if (ybottom < 0)
        ybottom = 0;

    if (ybottom > ny)
        ybottom = ny;

    count = (xright - xleft) * (ytop - ybottom);

    if (count == 0)
        return CPL_ERROR_ACCESS_OUT_OF_RANGE;

    count = 0;

    for (y = ybottom; y < ytop; y++) {
        for (x = xleft; x < xright; x++) {
            double value = data[x + y * nx];
            if (value < satur) {
                flux += value;
                count++;
            }
        }
    }

    if (count == 0)
        return CPL_ERROR_DIVISION_BY_ZERO;

    if (flux >= 0.0)
        error = sqrt(flux/gain);
    else
        error = sqrt(1/gain);

    /*
     * Flux correction for lost pixels
     */

    flux *= (float)npix / count;
    error *= (float)npix / count;
    
    flux /= area;  
    error /= area; 
    
    *o_flux = flux;
    *o_err = error;

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Compute median from a table column section corresponding to a slit
 *
 * @param table  Table with as many rows as rectified images
 * @param slits  Table with slits properties.
 * @param slit   Row in slits corresponding to slit to examine
 * @param label  Name of column to examine
 * @param mvalue Returned median value
 *
 * @return 0 in case of success.
 */

int mos_median_in_slit(cpl_table *table, cpl_table *slits, int slit, 
                       char *label, double *mvalue)
{
    int        position   = cpl_table_get_int(slits, "position", slit, NULL);
    int        length     = cpl_table_get_int(slits, "length", slit, NULL);
    if(position>=0)
    {
        cpl_table *tmp        = cpl_table_extract(table, position, length);

        *mvalue = cpl_table_get_column_median(tmp, label);
        cpl_table_delete(tmp);
    }
    else
        return 1;

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        return 1;

    return 0;
}


/**
 * @brief
 *   Convenience function for standard median filtering
 *
 * @param image  Image to smooth
 * @param nx     Filter size in x.
 * @param ny     Filter size in y.
 *
 * @return Filtered image
 */

cpl_image *mos_image_filter_median(cpl_image *image, int nx, int ny)
{
      cpl_mask  *kernel   = cpl_mask_new(nx, ny);
      cpl_image *filtered = cpl_image_new(cpl_image_get_size_x(image),
                                          cpl_image_get_size_y(image),
                                          cpl_image_get_type(image));

      cpl_mask_not(kernel);
      cpl_image_filter_mask(filtered, image, kernel,
                            CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
      cpl_mask_delete(kernel);

      return filtered;
}

/**
 * @brief
 *   Eliminate obvious outliers from table column values
 *
 * @param table       Input table
 * @param name        Column name
 *
 * @return 0 on success.
 *
 * The method is base on finite differences. First, the median
 * finite difference M is found. Then the median of the absolute
 * second order differences S is also found, and taken as a typical
 * deviation of the first order differences. Then values having both 
 * forward and backward differences more than 2S are invalidated.
 */

int mos_clean_outliers(cpl_table *table, const char *name)
{
    int    nrows = cpl_table_get_nrow(table);
    int    first_valid, last_valid;
    int    i, j, k;
    double trend, sigma;


    /*
     * Do something only if there are at least 10 valid values
     */

    if (nrows - cpl_table_count_invalid(table, name) < 10)
        return 0;

    /*
     * Start to compute differences at first valid element
     */

    for (i = 0; i < nrows; i++) {
        if (cpl_table_is_valid(table, name, i)) {
            first_valid = i;
            break;
        }
    }

    /*
     * Stop to compute differences at last valid element
     */

    for (i = nrows - 1; i >= 0; i--) {
        if (cpl_table_is_valid(table, name, i)) {
            last_valid = i;
            break;
        }
    }

    /*
     * Compute forward differences
     */

    cpl_table_new_column(table, "_D1_", CPL_TYPE_DOUBLE);

    for (i = first_valid; i < last_valid; i++) {
        if (cpl_table_is_valid(table, name, i)) {
            for (j = i + 1; j <= last_valid; j++) {
                if (cpl_table_is_valid(table, name, j)) {
                    double diff = cpl_table_get(table, name, j, NULL)
                                - cpl_table_get(table, name, i, NULL);

                    diff /= j - i;

                    for (k = i; k < j; k++) {
                        cpl_table_set(table, "_D1_", k, diff);
                    }
                    break;
                }
            }
        }
    }

    /*
     * Compute absolute forward differences of _D1_
     */

    cpl_table_new_column(table, "_D2_", CPL_TYPE_DOUBLE);

    for (i = first_valid; i < last_valid; i++) {
        if (cpl_table_is_valid(table, "_D1_", i)) {
            for (j = i + 1; j <= last_valid; j++) {
                if (cpl_table_is_valid(table, "_D1_", j)) {
                    double diff = cpl_table_get(table, "_D1_", j, NULL)
                                - cpl_table_get(table, "_D1_", i, NULL);

                    diff /= j - i;

                    for (k = i; k < j; k++) {
                        cpl_table_set(table, "_D2_", k, fabs(diff));
                    }
                    break;
                }
            }
        }
    }

    /*
     * Invalidate column elements where both forward and backward
     * differences are either invalid or deviate from the median
     * of _D1_ more than three times the median of _D2_.
     */

    trend = cpl_table_get_column_median(table, "_D1_");
    sigma = cpl_table_get_column_median(table, "_D2_");

    if (sigma == 0.0) {
        sigma = cpl_table_get_column_mean(table, "_D2_");
    }

    if (sigma == 0.0) {
        return 0.0;          // No outliers, not even noise...
    }

    for (i = first_valid; i <= last_valid; i++) {
        int    null;
        double t = cpl_table_get(table, "_D1_", i, &null);

        if (fabs(trend - t) > 2*sigma || null) {
            if (i > 0) {
                t = cpl_table_get(table, "_D1_", i - 1, &null);
                if (fabs(trend - t) > 2*sigma || null) {
                    cpl_table_set_invalid(table, name, i);
                }
            }
            else {
                cpl_table_set_invalid(table, name, i);
            }
        }
    }

    cpl_table_erase_column(table, "_D1_");
    cpl_table_erase_column(table, "_D2_");

//{
//char *savename = cpl_sprintf("D1_%s", name);
//cpl_table_name_column(table, "_D1_", savename);
//cpl_free(savename);
//savename = cpl_sprintf("D2_%s", name);
//cpl_table_name_column(table, "_D2_", savename);
//cpl_free(savename);
//}

    return 0;
}
    
