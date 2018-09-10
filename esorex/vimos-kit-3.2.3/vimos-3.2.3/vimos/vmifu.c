/* $Id: vmifu.c,v 1.13 2013-08-23 10:13:51 cgarcia Exp $
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
 * $Date: 2013-08-23 10:13:51 $
 * $Revision: 1.13 $
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
#include <float.h> /* for DBL_EPSILON, added by Peter Weilbacher */

#include <vmmath.h>
#include <vmfit.h>
#include <vmifu.h>

/* #include <qfits.h> */
#include "vmextractiontable.h"
#include <vmmoswavecalib.h>
#include <cpl_type.h>
#include <cpl_propertylist.h>
#include <cpl_error.h>
#include <cpl_memory.h>
#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>


#define N_SLITS          (4)
#define N_BLOCKS         (5)
#define FIBERS_PER_BLOCK (80)
#define FIBERS_STEP      (5)

#define MAX_COLNAME      (15)

/**
 * @name vmifu
 *
 * @doc
 *   The module vmifu collects low/medium level functions related to 
 *   IFU data reduction.
 */

/*@{*/

static double modelValue1D(double *c, int order, double position)
{

  double value  = 0.0;
  double factor = 1.0;
  int    j;

  
  for (j = 0; j <= order; j++) {
    value += c[j] * factor;
    factor *= position;
  }

  return value;

}

/***  Commented to avoid warning from the compiler (because unused)

static double modelValue2D(double *c, int order, double xpos, double ypos)
{

  double value  = 0.0;
  int    i, j, k;


  for (k = 0, j = 0; j <= order; j++)
    for (i = 0; i <= order - j; i++, k++)
      value += c[k] * ipow(xpos, i) * ipow(ypos, j);

  return value;

}

***/


static void drawModel(cpl_table *model, 
                      const char *colName, double *c, int order) 
{

  int    i;
  int    range;
  int   *idata;
  float *fdata;

  fdata = cpl_table_get_data_float(model, colName);
  idata = cpl_table_get_data_int(model, "y");
  range = cpl_table_get_nrow(model);
  cpl_table_fill_column_window_float(model, colName, 0, range, 0.0);

  for (i = 0; i < range; i++)
    fdata[i] = modelValue1D(c, order, idata[i]);

}


static int countRejections(VimosDpoint *list, 
                           int npix, double *c, int order, float tolerance)
{

  int    i, j;
  int    rejected = 0;
  double value;


  j = 0;
  for (i = 0; i < npix; i++) {

    value = modelValue1D(c, order, list[i].x);

    if (fabs(list[i].y - value) > tolerance)
      rejected++;
    else {
      if (j < i) {
        list[j].x = list[i].x;
        list[j].y = list[i].y;
      }
      j++;
    }

  }

  return rejected;

}


/*** 
 *** This is just a version of VmFrMedFil() adapted to a cpl_image.
 *** It should be removed as soon as a general median filtering for
 *** images is implemented in CPL.
 ***/

cpl_image *cpl_image_general_median_filter(cpl_image *ima_in, 
                                           int filtsizex, 
                                           int filtsizey,
                                           int excludeCenter) 
{

  char          task[] = "cpl_image_general_median_filter";

  cpl_image    *filt_img = NULL;
  int           col, row;
  float        *buf = NULL;
  float        *data;
  float        *fdata;
  int           medsize, upright_x, loleft_x, upright_y, loleft_y;
  int           uprightuse_x, loleftuse_x;
  int           i, j;
  int           xIsEven = !(filtsizex - (filtsizex/2)*2);
  int           yIsEven = !(filtsizey - (filtsizey/2)*2);
  float        *inpt;
  float        *outpt;
  int           f2x, f2y;
  int           nx = cpl_image_get_size_x(ima_in);
  int           ny = cpl_image_get_size_y(ima_in);


  if (xIsEven) filtsizex++;
  if (yIsEven) filtsizey++;

  if (nx <= filtsizex || ny <= filtsizey) {
    cpl_msg_error(task, "Median filter size: %dx%d, image size: %d,%d",
                filtsizex, filtsizey, nx, ny);
    return NULL;
  }

  if (excludeCenter) excludeCenter = 1;

  f2x = filtsizex / 2;
  f2y = filtsizey / 2;

  filt_img = cpl_image_duplicate(ima_in);
  buf = cpl_malloc(filtsizex * filtsizey * sizeof(float));
  data = cpl_image_get_data(ima_in);
  fdata = cpl_image_get_data(filt_img);

  for (row = 0; row < ny; row++) {
    loleft_y = row - f2y;
    upright_y = row + f2y + 1;

    for (col = 0; col < nx; col++) {
      loleft_x = col - f2x;
      loleftuse_x = MAX(loleft_x, 0);    /* Lowest x-value on image */
      medsize = filtsizex * filtsizey - excludeCenter;
      upright_x = col + f2x + 1;
      uprightuse_x = MIN(upright_x, nx); /* Highest x-value on image */

      /* Optimized extraction */

      outpt = buf;
      if (excludeCenter) {
        for (j = loleft_y; j < upright_y; j++) {

          if (j < 0)
            inpt = data + loleftuse_x;
          else if (j > ny - 1)
            inpt = data + loleftuse_x + (ny - 1) * nx;
          else
            inpt = data + loleftuse_x + j * nx;

          for (i = loleft_x; i < loleftuse_x; i++)
            *outpt++ = *inpt;

          for (i = loleftuse_x; i < uprightuse_x; i++) {
            if (i == col && j == row)
              inpt++;                /*** Skip "central" pixel value ***/
            else
              *outpt++ = *inpt++;
          }

          for (i = uprightuse_x; i < upright_x; i++)
            *outpt++ = *inpt;
        }
      }
      else {
        for (j = loleft_y; j < upright_y; j++) {

          if (j < 0)
            inpt = data + loleftuse_x;
          else if (j > ny - 1)
            inpt = data + loleftuse_x + (ny - 1) * nx;
          else
            inpt = data + loleftuse_x + j * nx;

          for (i = loleft_x; i < loleftuse_x; i++)
            *outpt++ = *inpt;

          for (i = loleftuse_x; i < uprightuse_x; i++)
            *outpt++ = *inpt++;

          for (i = uprightuse_x; i < upright_x; i++)
            *outpt++ = *inpt;
        }
      }

      fdata[col + row * nx] = medianPixelvalue(buf, medsize);
    }
  }

  cpl_free(buf);

  return filt_img;

}


/*** 
 *** This is a variation of cpl_image_general_median_filter,
 *** applying a median filter only in the vertical direction.
 ***/

cpl_image *cpl_image_vertical_median_filter(cpl_image *ima_in, 
                                            int filtsizey, int refrow, 
                                            int above, int below, int step)
{

  char       task[] = "cpl_image_general_median_filter";

  cpl_image *filt_img = NULL;
  int        col, row;
  float     *buf = NULL;
  float     *data;
  float     *fdata;
  int        upright_y, loleft_y;
  int        j;
  int        yIsEven = !(filtsizey - (filtsizey/2)*2);
  int        f2y;
  int        nx = cpl_image_get_size_x(ima_in);
  int        ny = cpl_image_get_size_y(ima_in);
  int        firstRow;


  if (yIsEven) filtsizey++;

  if (ny <= filtsizey) {
    cpl_msg_error(task, "Median filter size: %d, image size: %d", filtsizey, ny);
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
      fdata[col + row * nx] = medianPixelvalue(buf, filtsizey);
    }
  }

  cpl_free(buf);

  return filt_img;

}


/*** 
 ***  This should go to CPL some day... Currently works just with
 ***  float columns, it should work also for other numerical column.
 ***  The radius of the smooth box is hw (half width). No NULL elements
 ***  are allowed! Also this should be fixed.
 ***/

static int cpl_table_median_filter_column(cpl_table *table, 
                                          char *colName, int hw)
{

  int    box  = 2 * hw + 1;
  int    nrow = cpl_table_get_nrow(table);
  float *data = cpl_table_get_data_float(table, colName);
  float *smoo = cpl_malloc(nrow * sizeof(float));
  float *row  = cpl_malloc(box * sizeof(float));

  int    i, j;


  /* 
   * Copy first hw and last hw items 
   */

  for (i = 0; i < hw; i++) 
    smoo[i] = data[i];

  for (i = nrow - hw; i < nrow; i++) 
    smoo[i] = data[i];

  /* 
   * Median filtering
   */

  for (i = hw; i < nrow - hw; i++) {
    for (j = -hw; j <= hw; j++) 
      row[j + hw] = data[i + j];
    smoo[i] = median(row, box);
  }

  cpl_free(row);

  cpl_table_copy_data_float(table, colName, smoo);

  cpl_free(smoo);

  return 0;

}


/***
 ***  Also this should go to CPL some day... Currently works just with 
 ***  float columns, it should work also for other numerical columns. 
 ***  The radius of the filter box is hw (half width). No NULL elements
 ***  are allowed! Also this should be fixed.
 ***/

static int cpl_table_min_filter_column(cpl_table *table, 
                                       char *in_col, char *out_col, int hw)
{

  int    nrow = cpl_table_get_nrow(table);
  float *data = cpl_table_get_data_float(table, in_col);
  float *smoo;

  float  min;
  int    i, j;


  cpl_table_duplicate_column(table, out_col, table, in_col);
  smoo = cpl_table_get_data_float(table, out_col);


  /*
   * Min filtering. First and last part of the buffer are left untouched.
   */

  for (i = hw; i < nrow - hw; i++) {
    min = data[i];
    for (j = -hw; j <= hw; j++) {
      if (min > data[i + j])
        min = data[i + j];
    }
    smoo[i] = min;
  }

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
 *  Using the accurate reconstruction of the fiber profile obtained by the 
 *  functions ifuProfile() and rebinProfile(), an empirical relation has 
 *  been obtained. This function computes the value of R from the three 
 *  input values, and derives the corresponding value of dx by linear 
 *  interpolation on the tabulated values. If the condition 
 *
 *           v1  <= v2   and   v3  <  v2
 *  or
 *           v1  <  v2   and   v3  <=  v2
 *
 *  is not fulfilled, then this function returns the value 2.0.
 */

/*** Commented out because the correction r to dx is so small to be
     negligible.

static double values_to_dx(double v1, double v2, double v3)
{

  static double rt[] = {
    -0.50000000,
    -0.47308246,
    -0.42050084,
    -0.36780536,
    -0.32124185,
    -0.26728730,
    -0.22053576,
    -0.17183302,
    -0.12233649,
    -0.07436058,
    -0.02070184,
     0.00000000,
     0.02070184,
     0.07436058,
     0.12233649,
     0.17183302,
     0.22053576,
     0.26728730,
     0.32124185,
     0.36780536,
     0.42050084,
     0.47308246,
     0.50000000
  };

  static double dxt[] = {
    -0.500,
    -0.475,
    -0.425,
    -0.375,
    -0.325,
    -0.275,
    -0.225,
    -0.175,
    -0.125,
    -0.075,
    -0.025,
     0.000,
     0.025,
     0.075,
     0.125,
     0.175,
     0.225,
     0.275,
     0.325,
     0.375,
     0.425,
     0.475,
     0.500,
  };

  static int    n = sizeof(rt) / sizeof(double);
  static double epsilon = 0.00000001;

  double dx = 2.0;
  double r;
  int    i;


  if (v1 > v2 || v3 > v2) 
    return dx;

  if (2 * v2 - v1 - v3 < epsilon)
    return dx;

  r = 0.5 * (v3 - v1) / (2 * v2 - v3 - v1);

  for (i = 0; i < n; i++) {
    if (r > rt[i]) 
      continue;
    dx = dxt[i - 1] 
       + (r - rt[i - 1]) / (rt[i] - rt[i - 1]) * (dxt[i] - dxt[i - 1]);
    break;
  }

  return dx;

}

***/

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
 *  For a pixel at a distance dx from the spectral profile center, this
 *  function returns its theoretical value based on the modelled fiber
 *  profile. This value is normalized to the profile maximum value, i.e.,
 *  for dx = 0 this function returns 1. There are no limits to the value
 *  given for dx, but it's clear that beyond a certain distance this 
 *  function will just return 0.0. The returned value is computed from 
 *  an empirical table, linearly interpolating between values.
 */

static double dx_to_value(double dx)
{

/*  static double norm   = 3.017532;   */
  static double offset = 0.025;
  static double step   = 0.05;
  static double p[]    = {
    1.00000,
    0.99996,
    0.99600,
    0.99363,
    0.99024,
    0.98371,
    0.97754,
    0.96946,
    0.95954,
    0.94931,
    0.93719,
    0.92335,
    0.90941,
    0.89560,
    0.87707,
    0.86130,
    0.84158,
    0.82178,
    0.80195,
    0.77927,
    0.76022,
    0.73276,
    0.70892,
    0.68236,
    0.65778,
    0.63209,
    0.60109,
    0.57527,
    0.54045,
    0.51107,
    0.48302,
    0.45345,
    0.42497,
    0.39587,
    0.36674,
    0.34252,
    0.31339,
    0.28750,
    0.26153,
    0.24057,
    0.21730,
    0.19604,
    0.17451,
    0.15626,
    0.13760,
    0.12041,
    0.10540,
    0.09260,
    0.07954,
    0.06847,
    0.05968,
    0.05272,
    0.04530,
    0.04062,
    0.03373,
    0.02871,
    0.02387,
    0.02167,
    0.01669,
    0.01503,
    0.01247,
    0.00993,
    0.00758,
    0.00669,
    0.00526,
    0.00301,
    0.00289,
    0.00123,
    0.00065
  };

  static int n = sizeof(p) / sizeof(double);

  int   i;


  dx = fabs(dx);

  if (dx < offset)
    return 1.0;

  i = (dx - offset) / step;

  if (i < 0)
    i = 0;                            /* Safety belt */

  if (i >= n)
    return 0.0;

  return (p[i] + (dx - offset - step * i) / step * (p[i + 1] - p[i])); 
                         /* / norm;  */

}


/*
 * Internal utility, to derive the correct normalization factor
 * of the fiber profile model used by dx_to_value(). The sum of
 * the values of this model taken from a regular grid of step 1
 * should be 1. The actual constant is returned.
 */

void flux_constant()
{

  int     i, j, count;
  int     nvalues = 10;
  double  dx;
  double  step = 1./nvalues;
  double  total, subtotal;


  count = 0;
  total = 0.0;
  for (i = 0; i < nvalues; i++) {
    subtotal = 0.0;
    for (j = -4; j < 5; j++) {
      dx = j + i * step;
      subtotal += dx_to_value(dx);
    }
    printf("Subtotal = %f\n", subtotal);
    total += subtotal;
    count++;
  }
  printf("Total = %f\n", total / 3.017532 / count);

}


/*
 * alternative findPeak() using Gaussian fitting
 * (added by Peter Weilbacher)
 *
 * buffer: data buffer
 * npix:   length of data buffer
 * level:  min peak value
 * pos:    return peak center
 * Returns 0 if a failure occurs, 1 for success
 */

static int findPeakGaussian(double *buffer, int npix, double level, 
                            double *pos)
{
  cpl_vector *v = cpl_vector_wrap(npix, buffer),
             *p = cpl_vector_new(npix);
  int i, rc = CPL_ERROR_NONE;
  double center, sigma, area, bglevel, mse;

#ifdef DEBUG_PEAK
  cpl_bivector* biv;
#endif

  /* If at least one pixel of the profile is equal to zero,
   * a failure is returned.
   */
  for (i = 0; i < npix; i++) {
    if (fabs(buffer[i]) < DBL_EPSILON) {
      return 0;
    }
  }

  /* If two contiguous pixels of the profile
   * are identical, a failure is returned (suspected saturated line).
   */
  for (i = 1; i < npix; i++) {
    if (fabs(buffer[i] - buffer[i - 1]) < DBL_EPSILON) {
      return 0;
    }
  }

  /* fill positions into vector */
  for (i = 0; i < npix; i++) {
     cpl_vector_set(p, i, (double)i);
  }
#ifdef DEBUG_PEAK
  biv = cpl_bivector_wrap_vectors(p, v);
  cpl_bivector_dump(biv, stdout);
  cpl_bivector_unwrap_vectors(biv);
#endif

  rc = cpl_vector_fit_gaussian(p, NULL, v, NULL,
                               CPL_FIT_ALL,
                               &center, &sigma, &area, &bglevel, &mse,
                               NULL, NULL);
  buffer = cpl_vector_unwrap(v);

  if (rc == CPL_ERROR_NONE || rc == CPL_ERROR_CONTINUE ||
      rc == CPL_ERROR_SINGULAR_MATRIX) {
    cpl_error_reset();
    if (buffer[(int)floor(center)] < level &&
        buffer[(int)ceil(center)] < level) {
      /* Insignificant peak */
      return 0;
    }
    *pos = center; /* save only the measured center */
    return 1; /* it worked */
  } else {
    cpl_error_reset();
    return 0;
  }

}


/*
 *  This function can be used just to find the peak of an arc lamp line.
 *  It works under the assumption that the S/N is typically high. The
 *  peak is defined as a pixel having a value that is higher then the
 *  previous and the next pixels. The idea is to look for the peak
 *  with the highest value in the interval, and determine its position
 *  using a parabolic interpolation. The same will be repeated with
 *  findPeakGaussian(), and if the two positions will differ more
 *  than one pixel the parabolic interpolation result will be preferred
 *  (suspecting a blend between lines). This function returns 0 in case
 *  of failure. If at least one pixel of the profile is equal to zero, 
 *  a failure is returned. If two contiguous pixels of the profile
 *  are identical, a failure is returned (suspected saturated line).
 */

static int findPeak(double *buffer, int npix, double level, double *pos)
{

  double diffprev, diffnext;
  double max = 0;
  int    maxpos = 0;
  double gpos;
  int    i;

  npix--;

  for (i = 0; i < npix; i++)
    if (fabs(buffer[i]) < 0.000001)
      return 0;

  for (i = 1; i < npix; i++)
    if (fabs(buffer[i] - buffer[i - 1]) < 0.000001)
      return 0;

  for (i = 1; i < npix; i++) {
#ifdef DEBUG_PEAK
    printf("i=%d, buffer=%f\n",i,buffer[i]);
#endif
    if (buffer[i] > buffer[i - 1]) {
      if (buffer[i] > buffer[i + 1]) {
        if (buffer[i] > max || maxpos == 0) {
          max = buffer[i];
          maxpos = i;
        }
      }
    }
  }

#ifdef DEBUG_PEAK
  fflush(stdout);
#endif

  if (maxpos == 0)            /* No peaks - monotonic sequence of values */
    return 0;

  if (max < level)            /* Unsignificant peak                      */
    return 0;

  diffnext = buffer[maxpos + 1] - buffer[maxpos];
  diffprev = buffer[maxpos - 1] - buffer[maxpos];
  *pos = maxpos + 0.5 * (diffprev - diffnext) / (diffnext + diffprev);

  /*
   * Try to improve accuracy
   */

/* */
  npix++;
  if (findPeakGaussian(buffer, npix, level, &gpos)) {
    if (fabs(*pos - gpos) < 1.0) {
       *pos = gpos;
    }
  }
/* */

  return 1;

}


/*
 * This function just returns the position of the maximum value along 
 * the buffer.
 */

static int whereMax(double *buffer, int npix)
{

  double max    = buffer[0];
  int    maxpos = 0;
  int    i;


  for (i = 1; i < npix; i++) {
    if (buffer[i] > max) {
      max = buffer[i];
      maxpos = i;
    }
  }

  return maxpos;
  
}


/**
 * @memo
 *   Get IFU spectra extraction parameters.
 *
 * @return First guess fiber positions.
 *
 * @param grism    Grism identifier.
 * @param quadrant Quadrant number [1-4].
 * @param slit     IFU slit number [0-3].
 * @param mode     CCD or wavelength oriented extraction parameters.
 * @param row      Returned reference row.
 * @param above    Returned pixels to extract above reference row.
 * @param below    Returned pixels to extract below reference row.
 * @param zero     Returned expected position of zero order contamination.
 *
 * @doc
 *   This function returns the system default spectral extraction parameters,
 *   to be used when no configuration file is found. The returned parameters
 *   corresponding to HR grisms can be either CCD oriented or wavelength 
 *   oriented. When @em mode is set to 0, the returned values are CCD 
 *   oriented: this means that the pixels above and below the reference 
 *   row are chosen so that any spectrum is extracted in its whole extension 
 *   on the CCD. Alternatively, when @em mode is set to 1, the returned 
 *   values are wavelength oriented: this means that the pixels above and 
 *   below the reference row are chosen so that roughly the same wavelength 
 *   interval is covered in the four quadrants for a given grism, i.e., 
 *   only the overlapping spectral ranges in the four quadrants would be 
 *   extracted. For MR and LR grisms the @em mode is ignored, and the 
 *   parameters are always wavelength oriented.
 * 
 *   The grism identifier is a number with the following meaning:
 *
 *     0 = LR_red
 *     1 = LR_blue
 *     2 = MR
 *     3 = HR_red
 *     4 = HR_orange
 *     5 = HR_blue
 *     6 = HR_red holographic
 *     7 = HR_blue holographic
 *
 *   The reference row is chosen roughly within a region that is free from
 *   sky emission lines.
 *
 * @author C. Izzo
 */

int ifuExtractionParameters(int grism, int quadrant, int slit, int mode, 
                            int *row, int *above, int *below, int *zero)
{

  char       task[] = "ifuExtractionParameters";

  int        rows[4];
  int        aboves[4];
  int        zeros[4];
  int        exception;
  int        min, max;
  int        i;


  switch(grism) {
  case 0:
    switch(slit) {
    case 0:
      rows[0] = 1030;           /* LR red, quadrant 1, slit 0 (lowest)  */
      rows[1] = 1019;           /* LR red, quadrant 2, slit 0 (lowest)  */
      rows[2] = 3097;           /* LR red, quadrant 3, slit 0 (highest) */
      rows[3] = 2947;           /* LR red, quadrant 4, slit 0 (highest) */
      aboves[0] = 304;          /* 312; */
      aboves[1] = 304;          /* 311; */
      aboves[2] = 304;          /* 307; */
      aboves[3] = 303;
      zeros[0] = 1210;
      zeros[1] = 1196;
      zeros[2] = 0;
      zeros[3] = 0;
      *above = aboves[quadrant - 1];
      *below = 500 - *above;
      *row = rows[quadrant - 1];
      *zero = zeros[quadrant - 1];
      break;

    case 1:
      rows[0] = 2157;           /* LR red, quadrant 1, slit 1 */
      rows[1] = 2144;           /* LR red, quadrant 2, slit 1 */
      rows[2] = 1974;           /* LR red, quadrant 3, slit 1 */
      rows[3] = 1819;           /* LR red, quadrant 4, slit 1 */
      aboves[0] = 308;
      aboves[1] = 306;
      aboves[2] = 308;
      aboves[3] = 309;
      zeros[0] = 2356;
      zeros[1] = 2340;
      zeros[2] = 2159;
      zeros[3] = 2007;
      *above = aboves[quadrant - 1];
      *below = 500 - *above;
      *row = rows[quadrant - 1];
      *zero = zeros[quadrant - 1];
      break;

    case 2:
      rows[0] = 2726;           /* LR red, quadrant 1, slit 2 */
      rows[1] = 2710;           /* LR red, quadrant 2, slit 2 */
      rows[2] = 1405;           /* LR red, quadrant 3, slit 2 */
      rows[3] = 1252;           /* LR red, quadrant 4, slit 2 */
      aboves[0] = 304;
      aboves[1] = 304;
      aboves[2] = 310;
      aboves[3] = 310;
      zeros[0] = 0;
      zeros[1] = 0;
      zeros[2] = 0;
      zeros[3] = 0;
      *above = aboves[quadrant - 1];
      *below = 500 - *above;
      *row = rows[quadrant - 1];
      *zero = zeros[quadrant - 1];
      break;

    case 3:
      rows[0] = 3292;           /* LR red, quadrant 1, slit 3 (highest)  */
      rows[1] = 3276;           /* LR red, quadrant 2, slit 3 (highest)  */
      rows[2] = 835;            /* LR red, quadrant 3, slit 3 (lowest) */
      rows[3] = 679;            /* LR red, quadrant 4, slit 3 (lowest) */
      aboves[0] = 303;
      aboves[1] = 304;
      aboves[2] = 305;     /* 315 */
      aboves[3] = 305;     /* 315 */
      zeros[0] = 0;
      zeros[1] = 0;
      zeros[2] = 1021;
      zeros[3] = 866;
      *above = aboves[quadrant - 1];
      *below = 500 - *above;
      *row = rows[quadrant - 1];
      *zero = zeros[quadrant - 1];
      break;

    default:
      cpl_msg_error(task, "Wrong slit number");
      return 1;
    }
    break;

  case 1:
    switch(slit) {
    case 0:
      rows[0] = 1191;           /* LR blue, quadrant 1, slit 0 (lowest)  */
      rows[1] = 1181;           /* LR blue, quadrant 2, slit 0 (lowest)  */
      rows[2] = 3250;           /* LR blue, quadrant 3, slit 0 (highest) */
      rows[3] = 3098;           /* LR blue, quadrant 4, slit 0 (highest) */
      zeros[0] = 1361;
      zeros[1] = 1356;
      zeros[2] = 0;
      zeros[3] = 0;
      *above = 269;
      *below = 269;
      *row = rows[quadrant - 1];
      *zero = zeros[quadrant - 1];
      break;

    case 1:
      rows[0] = 2314;           /* LR blue, quadrant 1, slit 1 */
      rows[1] = 2304;           /* LR blue, quadrant 2, slit 1 */
      rows[2] = 2130;           /* LR blue, quadrant 3, slit 1 */
      rows[3] = 1974;           /* LR blue, quadrant 4, slit 1 */
      zeros[0] = 2504;
      zeros[1] = 2497;
      zeros[2] = 2318;
      zeros[3] = 2162;
      *above = 269;
      *below = 269;
      *row = rows[quadrant - 1];
      *zero = zeros[quadrant - 1];
      break;

    case 2:
      rows[0] = 2877;           /* LR blue, quadrant 1, slit 2 */
      rows[1] = 2867;           /* LR blue, quadrant 2, slit 2 */
      rows[2] = 1565;           /* LR blue, quadrant 3, slit 2 */
      rows[3] = 1411;           /* LR blue, quadrant 4, slit 2 */
      zeros[0] = 0;
      zeros[1] = 0;
      zeros[2] = 0;
      zeros[3] = 0;
      *above = 269;
      *below = 269;
      *row = rows[quadrant - 1];
      *zero = zeros[quadrant - 1];
      break;

    case 3:
      rows[0] = 3440;           /* LR blue, quadrant 1, slit 3 (highest)  */
      rows[1] = 3430;           /* LR blue, quadrant 2, slit 3 (highest)  */
      rows[2] = 1001;           /* LR blue, quadrant 3, slit 3 (lowest) */
      rows[3] = 848;            /* LR blue, quadrant 4, slit 3 (lowest) */
      zeros[0] = 0;
      zeros[1] = 0;
      zeros[2] = 1182;
      zeros[3] = 1026;
      *above = 269;
      *below = 269;
      *row = rows[quadrant - 1];
      *zero = zeros[quadrant - 1];
      break;

    default:
      cpl_msg_error(task, "Wrong slit number");
      return 1;
    }
    break;

  case 2:
    rows[0] = 2244;           /* MR, quadrant 1, slit 1 */
    rows[1] = 2234;           /* MR, quadrant 2, slit 1 */
    rows[2] = 2058;           /* MR, quadrant 3, slit 1 */
    rows[3] = 1895;           /* MR, quadrant 4, slit 1 */
    zeros[0] = 0;
    zeros[1] = 0;
    zeros[2] = 0;
    zeros[3] = 0;
    *above = 1175;
    *below = 918;
    *row = rows[quadrant - 1];
    *zero = zeros[quadrant - 1];
    break;
    
  case 6:
  case 3:
    rows[0] = 1718;           /* HR red, quadrant 1, slit 1 */
    rows[1] = 1712;           /* HR red, quadrant 2, slit 1 */
    rows[2] = 1514;           /* HR red, quadrant 3, slit 1 */
    rows[3] = 1560;           /* HR !!!orange!!!, quadrant 4, slit 1 */
    zeros[0] = 0;
    zeros[1] = 0;
    zeros[2] = 0;
    zeros[3] = 0;
    *row = rows[quadrant - 1];
    *zero = zeros[quadrant - 1];
    break;

  case 4:
    rows[0] = 1900;           /* HR orange, quadrant 1, slit 1 */
    rows[1] = 1893;           /* HR orange, quadrant 2, slit 1 */
    rows[2] = 1691;           /* HR orange, quadrant 3, slit 1 */
    rows[3] = 1515;           /* HR orange, quadrant 4, slit 1 */
    zeros[0] = 0;
    zeros[1] = 0;
    zeros[2] = 0;
    zeros[3] = 0;
    *row = rows[quadrant - 1];
    *zero = zeros[quadrant - 1];
    break;

  case 5:
    rows[0] = 3398;           /* HR blue, quadrant 1, slit 1 */
    rows[1] = 3388;           /* HR blue, quadrant 2, slit 1 */
    rows[2] = 3236;           /* HR blue, quadrant 3, slit 1 */
    rows[3] = 3080;           /* HR blue, quadrant 4, slit 1 */
    zeros[0] = 0;
    zeros[1] = 0;
    zeros[2] = 0;
    zeros[3] = 0;
    *row = rows[quadrant - 1];
    *zero = zeros[quadrant - 1];
    break;
  case 7:
    rows[0] = 3398;           /* HR blue holog, quadrant 1, slit 1 */
    rows[1] = 3450;           /* HR blue holog, quadrant 2, slit 1 */
    rows[2] = 3228;           /* HR blue holog, quadrant 3, slit 1 */
    rows[3] = 3046;           /* HR blue holog, quadrant 4, slit 1 */
    zeros[0] = 0;
    zeros[1] = 0;
    zeros[2] = 0;
    zeros[3] = 0;
    *above = 640;
    *below = 1900;
    *row = rows[quadrant - 1];
    *zero = zeros[quadrant - 1];
    break;

  default:
    cpl_msg_error(task, "Wrong grism");
    return 1;
  }

  /* For grism = 7 the spectra do not cover the whole CCD */
  if (grism > 2 && grism < 7) {

    /*
     *  The exception consists in having an HR orange grism on quadrant 4
     *  when HR red grism is used in all other quadrants. In this case the 
     *  definition for a common wavelength interval has not much sense, 
     *  and the whole CCD range is always returned, even if mode == 1. 
     *  In addition to that, the computation of the common wavelength 
     *  interval is limited to the first three quadrants.
     */

    exception = (grism == 3 && quadrant == 4);

    if (mode == 0 || exception) {
      *above = 4096 - *row - 5;
      *below = *row - 5;
    }
    else {
      min = max = rows[0];
      for (i = 1; i < 4; i++) {
        if (rows[i] < min)
          min = rows[i];
        if (rows[i] > max)
          max = rows[i];
        if (i == 2 && exception) /* Leave at quadrant 3 in HR orange case */
          break;
      }
      *above = 4096 - max - 5;
      *below = min - 5;
    }
  }

  return 0;

}


/**
 * @memo
 *   Get IFU wavelength calibration first guess.
 *
 * @return Coefficients of IDS polynomial.
 *
 * @param grism    Grism identifier.
 * @param quadrant Quadrant number [1-4].
 * @param slit     IFU slit number [0-3].
 * @param order    Returned order of the IDS polynomial.
 * @param lambda   Returned reference wavelength.
 *
 * @doc
 *   This function returns the system default spectral IDS coefficients,
 *   to be used when no configuration file is found. The same set of
 *   coefficients is returned for a given grism, quadrant, and pseudo-slit.
 *
 *   The grism identifier is a number with the following meaning:
 *
 *     0 = LR_red
 *     1 = LR_blue
 *     2 = MR
 *     3 = HR_red
 *     4 = HR_orange
 *     5 = HR_blue
 *     6 = HR_red holographic
 *     7 = HR_blue holographic
 *
 *   The pseudo-slit 0 is the one that is most separated from the others.
 *   For MR anf HR grisms the slit number has no effect, as it is always
 *   assumed to be 1. The returned polynomial transforms the wavelength
 *   difference from the reference wavelength to a CCD pixel position
 *   along the dispersion direction.
 *
 * @author C. Izzo
 */

double *ifuFirstIds(int grism, int quadrant, int slit, 
                    int *order, double *lambda)
{

  double *c;


  if (grism > 1 && grism != 7)
    *order = 4;
  else
    *order = 3;

  c = cpl_malloc((*order + 1) * sizeof(double));

  switch (grism) {
  case 0:                                /* LR_red */
    *lambda = 7635.105;
    switch (quadrant) {
    case 1:
      switch (slit) {
      case 0:
        c[0] = 1106;
        c[1] = 1.42721E-01;
        c[2] = -1.39587E-06;
        c[3] = 2.05069E-10;
        break;
      case 1:
        c[0] = 2232;
        c[1] = 1.40079E-01;
        c[2] = -1.41532E-06;
        c[3] = 1.91742E-10;
        break;
      case 2:
        c[0] = 2802;
        c[1] = 1.39292E-01;
        c[2] = -1.36808E-06;
        c[3] = 1.82394E-10;
        break;
      case 3:
        c[0] = 3366;
        c[1] = 1.38695E-01;
        c[2] = -1.52763E-06;
        c[3] = 1.73369E-10;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
    case 2:
      switch (slit) {
      case 0:
        c[0] = 1092;
        c[1] = 1.42768E-01;
        c[2] = -1.50645E-06;
        c[3] = 2.16611E-10;
        break;
      case 1:
        c[0] = 2214;
        c[1] = 1.40184E-01;
        c[2] = -1.44645E-06;
        c[3] = 1.89765E-10;
        break;
      case 2:
        c[0] = 2780;
        c[1] = 1.39336E-01;
        c[2] = -1.32686E-06;
        c[3] = 2.39379E-10;
        break;
      case 3:
        c[0] = 3350;
        c[1] = 1.38752E-01;
        c[2] = -1.45132E-06;
        c[3] = 2.13249E-10;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
    case 3:
      switch (slit) {
      case 0:
        c[0] = 3168;
        c[1] = 1.39195E-01;
        c[2] = -1.42261E-06;
        c[3] = 1.69776E-10;
        break;
      case 1:
        c[0] = 2046;
        c[1] = 1.40507E-01;
        c[2] = -1.36907E-06;
        c[3] = 3.15178E-10;
        break;
      case 2:
        c[0] = 1480;
        c[1] = 1.42087E-01;
        c[2] = -1.50725E-06;
        c[3] = 1.90814E-10;
        break;
      case 3:
        c[0] = 912;
        c[1] = 1.43513E-01;
        c[2] = -1.44775E-06;
        c[3] = 1.90265E-10;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
    case 4:
      switch (slit) {
      case 0:
        c[0] = 3018;
        c[1] = 1.39187E-01;
        c[2] = -1.39253E-06;
        c[3] = 2.44681E-10;
        break;
      case 1:
        c[0] = 1894;
        c[1] = 1.40789E-01;
        c[2] = -1.43792E-06;
        c[3] = 2.71668E-10;
        break;
      case 2:
        c[0] = 1326;
        c[1] = 1.42259E-01;
        c[2] = -1.49931E-06;
        c[3] = 2.32267E-10;
        break;
      case 3:
        c[0] = 754;
        c[1] = 1.43676E-01;
        c[2] = -1.37310E-06;
        c[3] = 2.04999E-10;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
    default:
      cpl_free(c);
      return NULL;
    }
    break;
  case 1:                                /* LR_blue */
    *lambda = 5015.675;
    switch (quadrant) {
    case 1:
      switch (slit) {
      case 0:
        c[0] = 1128;
        c[1] = 1.92205E-01;
        c[2] = -3.71742E-06;
        c[3] = 6.87774E-10;
        break;
      case 1:
        c[0] = 2252;
        c[1] = 1.88831E-01;
        c[2] = -3.38905E-06;
        c[3] = 6.43705E-10;
        break;
      case 2:
        c[0] = 2820;
        c[1] = 1.87606E-01;
        c[2] = -2.69320E-06;
        c[3] = 3.23760E-10;
        break;
      case 3:
        c[0] = 3386;
        c[1] = 1.87251E-01;
        c[2] = -3.24139E-06;
        c[3] = 5.89868E-10;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
    case 2:
      switch (slit) {
      case 0:
        c[0] = 1120;
        c[1] = 1.92105E-01;
        c[2] = -3.74334E-06;
        c[3] = 7.80496E-10;
        break;
      case 1:
        c[0] = 2244;
        c[1] = 1.88911E-01;
        c[2] = -3.54657E-06;
        c[3] = 7.73775E-10;
        break;
      case 2:
        c[0] = 2810;
        c[1] = 1.88185E-01;
        c[2] = -3.83450E-06;
        c[3] = 7.98968E-10;
        break;
      case 3:
        c[0] = 3372;
        c[1] = 1.87408E-01;
        c[2] = -3.24542E-06;
        c[3] = 5.89129E-10;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
    case 3:
      switch (slit) {
      case 0:
        c[0] = 3196;
        c[1] = 1.87613E-01;
        c[2] = -3.64721E-06;
        c[3] = 6.93420E-10;
        break;
      case 1:
        c[0] = 2070;
        c[1] = 1.89155E-01;
        c[2] = -3.48750E-06;
        c[3] = 6.66328E-10;
        break;
      case 2:
        c[0] = 1500;
        c[1] = 1.90859E-01;
        c[2] = -4.25190E-06;
        c[3] = 9.84393E-10;
        break;
      case 3:
        c[0] = 936;
        c[1] = 1.92296E-01;
        c[2] = -3.22008E-06;
        c[3] = 5.27671E-10;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
    case 4:
      switch (slit) {
      case 0:
        c[0] = 3040;
        c[1] = 1.87776E-01;
        c[2] = -3.55411E-06;
        c[3] = 6.93461E-10;
        break;
      case 1:
        c[0] = 1916;
        c[1] = 1.89409E-01;
        c[2] = -3.64742E-06;
        c[3] = 7.44002E-10;
        break;
      case 2:
        c[0] = 1350;
        c[1] = 1.91058E-01;
        c[2] = -3.59109E-06;
        c[3] = 6.58365E-10;
        break;
      case 3:
        c[0] = 782;
        c[1] = 1.92951E-01;
        c[2] = -3.81789E-06;
        c[3] = 7.85007E-10;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
    default:
      cpl_free(c);
      return NULL;
    }
    break;
  case 2:                                       /* MR */
    *lambda = 7635.105;
    switch (quadrant) {
    case 1:
      c[0] = 2428;
      c[1] = 3.85953E-01;
      c[2] = -2.93889E-06;
      c[3] = 4.44543E-10;
      c[4] = -9.78740E-14;
      break;
    case 2:
      c[0] = 2420;
      c[1] = 3.87053E-01;
      c[2] = -2.99889E-06;
      c[3] = 4.76732E-10;
      c[4] = -8.91030E-14;
      break;
    case 3:
      c[0] = 2252;
      c[1] = 3.89297E-01;
      c[2] = -3.27377E-06;
      c[3] = 4.68983E-10;
      c[4] = -7.18225E-14;
      break;
    case 4:
      c[0] = 2094;
      c[1] = 3.89315E-01;
      c[2] = -3.15185E-06;
      c[3] = 5.11989E-10;
      c[4] = -8.92864E-14;
      break;
    default:
      cpl_free(c);
      return NULL;
    }
    break;
  case 3:                                /* HR_red */
    *lambda = 7245.167;
    switch (quadrant) {
    case 1:
      c[0] = 1944;
      c[1] = 1.69035;
      c[2] = -1.12423E-04;
      c[3] = 2.11259E-08;
      c[4] = -3.69459E-12;
      break;
    case 2:
      c[0] = 1942;
      c[1] = 1.68526;
      c[2] = -1.10435E-04;
      c[3] = 2.10878E-08;
      c[4] = -4.26251E-12;
      break;
    case 3:
      c[0] = 1762;
      c[1] = 1.71871;
      c[2] = -1.18593E-04;
      c[3] = 2.14324E-08;
      c[4] = -3.37848E-12;
      break;
    case 4:                     /* This is actually an HR_orange */
      *lambda = 6598.953;
      c[0] = 2434;
      c[1] = 1.60816;
      c[2] = -8.77165E-05;
      c[3] = 1.52295E-08;
      c[4] = -5.61383E-12;
      break;
    default:
      cpl_free(c);
      return NULL;
    }
    break;
  case 4:                              /* HR_orange */
    *lambda = 6598.953;
    switch (quadrant) {
    case 1:
      c[0] = 2744;
      c[1] = 1.60;                     /* 1.59241; */
      c[2] = -8.51278E-05;
      c[3] = 1.30080E-08;
      c[4] = -5.93688E-12;
      break;
    case 2:
      c[0] = 2748;
      c[1] = 1.58803;
      c[2] = -8.39732E-05;
      c[3] = 1.27792E-08;
      c[4] = -5.26719E-12;
      break;
    case 3:
      c[0] = 2574;
      c[1] = 1.60;                     /* 1.59546;  */
      c[2] = -8.65212E-05;
      c[3] = 1.41364E-08;
      c[4] = -4.32910E-12;
      break;
    case 4:
      c[0] = 2434;
      c[1] = 1.60816;
      c[2] = -8.77165E-05;
      c[3] = 1.52295E-08;
      c[4] = -5.61383E-12;
      break;
    default:
      cpl_free(c);
      return NULL;
    }
    break;
  case 5:                             /* HR_blue */
    *lambda = 5015.675;
    switch (quadrant) {
    case 1:
      c[0] = 2056;
      c[1] = 1.79965;
      c[2] = -1.14937E-04;
      c[3] = 2.95681E-08;
      c[4] = -8.40936E-12;
      break;
    case 2:
      c[0] = 2046;
      c[1] = 1.80073;
      c[2] = -1.14859E-04;
      c[3] = 2.72775E-08;
      c[4] = -6.89394E-12;
      break;
    case 3:
      c[0] = 1876;
      c[1] = 1.83359;
      c[2] = -1.23169E-04;
      c[3] = 2.94921E-08;
      c[4] = -6.77784E-12;
      break;
    case 4:
      c[0] = 1726;
      c[1] = 1.83450;
      c[2] = -1.24277E-04;
      c[3] = 3.10028E-08;
      c[4] = -6.91772E-12;
      break;
    default:
      cpl_free(c);
      return NULL;
    }
    break;
  case 6:                           /* HR_red holographic */
    *lambda = 7245.167;
    switch (quadrant) {
    case 1:
      c[0] = 1804;
      c[1] = 1.67;
      c[2] = -0.00004015;
      c[3] = 1.0E-08;
      c[4] = -1.75E-12;
      break;
    case 2:
      c[0] = 1777;
      c[1] = 1.672;
      c[2] = -0.00004011;
      c[3] = 1.0E-08;
      c[4] = -1.8E-12;
      break;
    case 3:
      c[0] = 1597;
      c[1] = 1.677;
      c[2] = -0.0000426;
      c[3] = 9.0E-09;
      c[4] = -1.0E-12;
      break;
    case 4:
      c[0] = 1447;
      c[1] = 1.6765;
      c[2] = -0.0000425;
      c[3] = 1.0E-08;
      c[4] = -1.5E-12;
      break;
    default:
      cpl_free(c);
      return NULL;
    }
    break;
  case 7:                             /* HR_blue holographic */
      *lambda = 5015.675;
      switch (quadrant) {
      case 1:
        c[0] = 3312;
        c[1] = 1.3561;
        c[2] = -4.4e-5;
        c[3] = 1.7e-8;
        break;
      case 2:
        c[0] = 3312;
        c[1] = 1.3561;
        c[2] = 0;
        c[3] = 0;
        break;
      case 3:
        c[0] = 3312;
        c[1] = 1.3561;
        c[2] = 0;
        c[3] = 0;
        break;
      case 4:
        c[0] = 3312;
        c[1] = 1.3561;
        c[2] = 0;
        c[3] = 0;
        break;
      default:
        cpl_free(c);
        return NULL;
      }
      break;
  }

  return c;

}


/**
 * @memo
 *   Get default IFU spectra resampling parameters.
 *
 * @return 0 on success.
 *
 * @param grism       Grism identifier.
 * @param startLambda Returned conventional start wavelength.
 * @param endLambda   Returned conventional end wavelength.
 * @param stepLambda  Returned constant sampling step.
 *
 * @doc
 *   This function returns the system default spectral resampling
 *   parameters, i.e., the wavelength range, and the constant
 *   resampling wavelength step for a given grism.
 *
 *   The grism identifier is a number with the following meaning:
 *
 *     0 = LR_red
 *     1 = LR_blue
 *     2 = MR
 *     3 = HR_red
 *     4 = HR_orange
 *     5 = HR_blue
 *     6 = HR_red holographic
 *     7 = HR_blue holographic
 *
 * @author C. Izzo
 */

int ifuRange(int grism, 
             double *startLambda, double *endLambda, double *stepLambda)
{

  switch (grism) {
  case 0:                          /* LR_red */
    *startLambda = 5500.;
    *endLambda = 10000.;
    *stepLambda = 7.0;
    break;
  case 1:                          /* LR_blue */
    *startLambda = 3500.;
    *endLambda = 7000.;
    *stepLambda = 5.2;
    break;
  case 2:                          /* MR */
    *startLambda = 4000.;
    *endLambda = 11000.;
    *stepLambda = 2.6;
    break;
  case 6:                          /* HR_red holographic */
  case 3:                          /* HR_red */
    *startLambda = 6100.;
    *endLambda = 8900.;
    *stepLambda = 0.58;
    break;
  case 4:                          /* HR_orange */
    *startLambda = 5100.;
    *endLambda = 7700.;
    *stepLambda = 0.62;
    break;
  case 5:                          /* HR_blue */
    *startLambda = 4000.;
    *endLambda = 6300.;
    *stepLambda = 0.54;
    break;
  case 7:                          /* HR_blue holographic */ 
    *startLambda = 3450.;
    *endLambda = 5350.;
    *stepLambda = 0.71;
    break;
  default:
    return 1;
  }
  return 0;

}


/**
 * @memo
 *   Get default IFU spectra range where transmission is determined.
 *
 * @return 0 on success.
 *
 * @param grism       Grism identifier.
 * @param startLambda Returned start wavelength.
 * @param endLambda   Returned end wavelength.
 *
 * @doc
 *   This function returns the system default spectral interval where
 *   the fiber-to-fiber transmission correction is determined on a
 *   flat field exposure for a given grism. Such wavelength intervals
 *   are meant to select the brightest part of the spectra, avoiding
 *   possible zero order contaminations.
 *
 *   The grism identifier is a number with the following meaning:
 *
 *     0 = LR_red
 *     1 = LR_blue
 *     2 = MR
 *     3 = HR_red
 *     4 = HR_orange
 *     5 = HR_blue
 *     6 = HR_red holographic 
 *     7 = HR_blue holographic 
 *
 * @author C. Izzo
 */

int ifuRangeTransmission(int grism, double *startLambda, double *endLambda)
{

  switch (grism) {
  case 0:                          /* LR_red */
    *startLambda = 6500.;
    *endLambda = 8000.;
    break;
  case 1:                          /* LR_blue */
    *startLambda = 4600.;
    *endLambda = 6100.;
    break;
  case 2:                          /* MR */
    *startLambda = 6500.;
    *endLambda = 8000.;
    break;
  case 6:                          /* HR_red holographic */
  case 3:                          /* HR_red */
    *startLambda = 7000.;
    *endLambda = 8000.;
    break;
  case 4:                          /* HR_orange */
    *startLambda = 6000.;
    *endLambda = 7000.;
    break;
  case 7:                          /* HR_blue holographic TODO: Review*/
    *startLambda = 4000.;
    *endLambda = 5000.;
    break;
  case 5:                          /* HR_blue */
    *startLambda = 4700.;
    *endLambda = 5700.;
    break;
  default:
    return 1;
  }
  return 0;

}


/**
 * @memo
 *   Bias subtraction.
 *
 * @return Pointer to the bias subtracted image.
 *
 * @param image   Image of an IFU exposure.
 * @param bias    Master bias.
 *
 * @doc
 *   This is just a quick replacement of the VIMOS DRS bias
 *   subtraction function, adapted to CPL images. This function
 *   assumes that the master bias has the overscan regions trimmed,
 *   so they are grown back before subtraction. After this, an
 *   overscan correction is applied to the result, and finally 
 *   the overscan regions are trimmed. The new image is returned.
 *
 * @author C. Izzo
 */

cpl_image *removeBias(cpl_image *image, cpl_image *mbias)
{

  int     nx    = cpl_image_get_size_x(image);
  int     ny    = cpl_image_get_size_y(image);
  float  *data  = cpl_image_get_data(image);
  int     bnx   = cpl_image_get_size_x(mbias);
  int     bny   = cpl_image_get_size_y(mbias);
  float  *bdata = cpl_image_get_data(mbias);

  cpl_image *grown;
  float     *gdata;
  float     *edata;

  int     dx;
  float   residual;


  if (mbias) {

    /*
     *  Overscans are assumed to lay on the right and left sides of image,
     *  and to have equal width. Then this is the width of one overscan 
     *  region:
     */

    dx = (nx - bnx) / 2;

    grown = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
    gdata = cpl_image_get_data(grown);

    insertFloatImage(gdata, nx, ny, dx, 0, bnx, bny, bdata);

    edata = extractFloatImage(bdata, bnx, bny, 0, 0, dx, bny);
    insertFloatImage(gdata, nx, ny, 0, 0, dx, ny, edata);
    free(edata);

    edata = extractFloatImage(bdata, bnx, bny, bnx - dx - 1, 0, dx, bny);
    insertFloatImage(gdata, nx, ny, nx - dx - 1, 0, dx, ny, edata);
    free(edata);

    cpl_image_subtract(image, grown);

    cpl_image_delete(grown);

  }
  else
    dx = 50;  /*  It should be read from the header, but what the hell...  */

  edata = extractFloatImage(data, nx, ny, 0, 0, dx, ny);
  residual = medianPixelvalue(edata, dx * ny);
  free(edata);

  edata = extractFloatImage(data, nx, ny, nx - dx - 1, 0, dx, ny);
  residual += medianPixelvalue(edata, dx * ny);
  free(edata);

  residual /= 2;

  cpl_image_subtract_scalar(image, residual);

  return cpl_image_extract(image, dx + 1, 1, nx - dx, ny);

}


/**
 * @memo
 *   Remove the bias constant from an image.
 * 
 * @return Pointer to the bias subtracted image.
 * 
 * @param image   Image of an IFU exposure.
 * 
 * @doc
 *   The mean bias level is evaluated from the overscan regions of
 *   the input frame. Currently the values are hardcoded - 50 pixels
 *   in the X direction for each overscan are taken. After subtraction
 *   the overscan regions are trimmed and the new image is returned.
 *   
 * @author C. Izzo
 */

cpl_image *removeBiasLevel(cpl_image *image)
{

  int     nx    = cpl_image_get_size_x(image);
  int     ny    = cpl_image_get_size_y(image);
  float  *data  = cpl_image_get_data(image);

  float     *edata;

  int     dx = 50; 
  float   residual;


  edata = extractFloatImage(data, nx, ny, 0, 0, dx, ny);
  residual = medianPixelvalue(edata, dx * ny);
  free(edata);

  edata = extractFloatImage(data, nx, ny, nx - dx - 1, 0, dx, ny);
  residual += medianPixelvalue(edata, dx * ny);
  free(edata);

  residual /= 2;

  cpl_image_subtract_scalar(image, residual);

  return cpl_image_extract(image, dx + 1, 1, nx - dx, ny);

}


/**
 * @memo
 *   Identify fibers. 
 *
 * @return 0 on success.
 *
 * @param image    Image of an IFU flat field exposure.
 * @param refrow   Image row where the fibers are to be identified.
 * @param refdata  Reference image row with identified fiber peaks.
 * @param ident    Identified fibers positions table.
 * @param radius   Correlation radius.
 * @param wradius  Correlation half-window.
 *
 * @doc
 *   Identifying the IFU fibers here means to give an approximate
 *   position for each IFU fiber along the specified image row.
 *   This means to give the expected positions for dead fibers too.
 *   @em refdata and @em ident are typically produced from a task
 *   external to the automatic pipeline, where the fibers are identified
 *   interactively, through a trial-and-error procedure.
 *
 *   The algorithm is based on the cross-correlation between the
 *   flat field image at the reference row, and the reference
 *   image consisting of a row of peaks where the corresponding
 *   fiber spectra were already identified. The flat field image 
 *   at the reference row is divided into 5 contiguous sections 
 *   of 409 pixels each. The first section begins at the second
 *   pixel on the left (pixel 1). The correlation window is twice
 *   the value of @em wradius plus 1. The correlation is made on
 *   an interval of twice the value of @em radius plus 1. For this
 *   reason, 2 * (@em radius + @em wradius ) + 1 must be less than
 *   409, or an error will be returned. The correlation is tried
 *   for each one of the sections. The median of the 5 offsets
 *   obtained from the correlations is applied to the positions
 *   listed in the @em ident table.
 *
 * @author C. Izzo
 */

int ifuIdentifyUpgrade(cpl_image *image, int refrow, float *refdata,
                       cpl_table *ident, int radius, int wradius)
{

  char    task[]   = "ifuIdentifyUpgrade";

  int     nx       = cpl_image_get_size_x(image);
  float  *data     = cpl_image_get_data(image);
  float  *line     = data + refrow * nx;
  double *normdata = NULL;
  double *normref  = NULL;
  double *cross    = NULL;

  int     countFibers = N_BLOCKS * FIBERS_PER_BLOCK;
  int     count = 0;
  int     firstpix = 1;
  int     secSize  = 680;
  int     secCount = 3;
  float   offset[5];                   /* Size equal to secCount */
  float   shift;
  float   fpos, max;
  float   dx;
  float  *posdata;
  int     maxpos;
  double  sum;
  int     length   = 2 * radius + 1;
  int     wlength  = 2 * wradius + 1;
  int     start;
  int     i, j, k;


  normdata = cpl_malloc(secSize * sizeof(double));
  normref = cpl_malloc(wlength * sizeof(double));
  cross = cpl_malloc(length * sizeof(double));

  for (i = 0; i < secCount; i++) {

    /*
     *  Load data from current section, and normalize them.
     */

    start = firstpix + i * secSize;
    max = line[start];
    for (j = 0; j < secSize; j++) {
      normdata[j] = line[start + j];
      if (normdata[j] > max)
        max = normdata[j];
    }

    if (fabs(max) < 0.000001) {
      cpl_free(normdata);
      cpl_free(normref);
      cpl_free(cross);
      return 1;
    }

    for (j = 0; j < secSize; j++)
      normdata[j] /= max;

    /*
     *  Load data from the corresponding reference image section,
     *  and normalize them.
     */

    start += secSize / 2 - wradius;
    max = refdata[start];
    for (j = 0; j < wlength; j++) {
      normref[j] = refdata[start + j];
      if (normref[j] > max)
        max = normref[j];
    }

    if (fabs(max) < 0.000001) {
      cpl_free(normdata);
      cpl_free(normref);
      cpl_free(cross);
      return 1;
    }

    for (j = 0; j < wlength; j++) 
      normref[j] /= max;

    /*
     *  Cross-correlation of normalized data.
     */

    start = secSize / 2 - wradius - radius;
    for (j = 0; j < length; j++) {
      sum = 0.0;
      for (k = 0; k < wlength; k++) 
        sum += normref[k] * normdata[start + j + k];
      cross[j] = sum;
    }

    max = cross[0];
    maxpos = 0;
    for (j = 1; j < length; j++) {
      if (cross[j] > max) {
        max = cross[j];
        maxpos = j;
      }
    }

    /*
     * The offset is set to a value greater than the correlation
     * radius in case of error.
     */

    offset[i] = radius + 1;
    if (maxpos != 0 && maxpos != length - 1) {
      dx = values_to_dx(cross[maxpos - 1], cross[maxpos], cross[maxpos + 1]);
      if (dx < 1.0)
        offset[i] = maxpos - radius + dx;
    }

  }


  /*
   *  Find median offset among good offsets
   */

  j = 0;
  for (i = 0; i < secCount; i++) {
    if (offset[i] < radius) {
      if (i > j)
        offset[j] = offset[i];
      ++j;
    }
  }

  if (j == 0) {
    cpl_free(normdata);
    cpl_free(normref);
    cpl_free(cross);
    return 1;
  }

  shift = median(offset, j);

  /*
   *  Ensure that missing data (marked by zeroes) are not brought in
   *  by a positive shift.
   */

  posdata = cpl_table_get_data_float(ident, "Position");
  for (i = 0; i < countFibers; i++) {
    if (posdata[i] < 0.0001) {
      posdata[i] -= length;
      count++;
    }
  }

  if (count == countFibers) {
    cpl_free(normdata);
    cpl_free(normref);
    cpl_free(cross);
    return 1;
  }

  cpl_msg_info(task, "Cross-correlation offset "
             "with reference identification: %f", shift);

  cpl_table_add_scalar(ident, "Position", shift);

  for (i = 0; i < countFibers; i++) {
    fpos = cpl_table_get_float(ident, "Position", i, NULL);
    if (!fiberPeak(image, refrow, &fpos, NULL))
      cpl_table_set_float(ident, "Position", i, fpos);
  }

  cpl_free(normdata);
  cpl_free(normref);
  cpl_free(cross);

  return 0;

}


/**
 * @memo
 *   Identify fibers.
 *
 * @return Pointer to the identification table.
 *
 * @param image   Image of an IFU flat field exposure.
 * @param row     Image row where the fibers are to be identified.
 *
 * @doc
 *   Identifying the IFU fibers here means to give an approximate
 *   position for each IFU fiber along the specified image row.
 *   This means to give the expected positions for dead fibers too.
 *   The algorithm is implemented in two basic steps: 1) determination
 *   of the position of the gaps between blocks, and 2) determination
 *   of the fibers positions within each block. The first step is 
 *   carried out with a folding analysis, and the second by correlating
 *   the fiber signal of each block with a "comb" of 80 equally spaced 
 *   teeth. Ambiguities may arise for the first and the last blocks, 
 *   in case they are vignetted: for a bullet-proof identification a 
 *   first-guess identification should be prepared, and passed to the 
 *   function ifuIdentifyUpgrade(). The advantage of ifuIdentify() lays 
 *   in its generality, that permits to reduce IFU data also when a 
 *   first-guess fiber identification is missing, turning out to be 
 *   especially useful in the preparation of such first guesses.
 *
 * @author C. Izzo
 */

cpl_table *ifuIdentify(cpl_image *image, int refrow)
{

  char    task[] = "ifuIdentify";

  int     nx   = cpl_image_get_size_x(image);
  int     ny   = cpl_image_get_size_y(image);
  float  *data = cpl_image_get_data(image);
  float  *line = data + refrow * nx;

  cpl_table *ident;
  cpl_table *folds;
  cpl_table *stats;
  float     *folded;
  int       *position;
  int       *count;
  int        p, i, j, k, n, pos;
#ifdef CPL_SIZE_FORMAT
  cpl_size   row;
#else
  int        row;
#endif
  int        startp     = 402;
  int        endp       = 412;
  int        interval;
  int        ntrials;
  int        ups, downs;
  float      period;
  float      gap[N_BLOCKS];
  float      min, frow;
  float      meanLevel1, meanLevel2, meanLevel3;
  double     level, plevel;
  char       colName[MAX_COLNAME];

  float      fpos;
  float      max;
  float      candidate[20];
  int        peak[20];


  cpl_msg_debug(task, "Identify fibers in image row %d", refrow);

  if (refrow < 0 || refrow >= ny) {
    cpl_msg_error(task, "Image row %d out of bounds.", refrow);
    return NULL;
  }

  /*
   *  A row of the input image is pointed by 'line'. The signal in
   *  this line is folded with different periods, starting from
   *  'startp' and ending with 'endp' (included). This is done to
   *  obtain the mean pattern for each IFU block, and determine
   *  the positions of the gaps between blocks. The pattern of a
   *  block is expected to be as long as the distance between 80 - 1 
   *  fibers plus the width of a gap. Being the fiber-to-fiber 
   *  distance about 5 pixels, and a gap width about 10 pixels,
   *  the pattern of a block should be around 405 pixels long. 
   *  The choice of 'startp' and 'endp' is constrained by this.
   *  A table containing all the obtained folded profiles is
   *  created. The folded profile containing the minimum value
   *  is found, and the position of this minimum is determined.
   *  This position roughly corresponds to the X coordinate of 
   *  the first gap in the image. 
   */

  folds = cpl_table_new(endp);
  cpl_table_new_column(folds, "count", CPL_TYPE_INT);
  count = cpl_table_get_data_int(folds, "count");

  /*
   *  The "count" column is used to normalize the folded profile.
   */

  for (i = 0, p = startp; p <= endp; i++, p++) {
    snprintf(colName, MAX_COLNAME, "p%d", p);
    cpl_table_fill_column_window_int(folds, "count", 0, p, 0);
    cpl_table_new_column(folds, colName, CPL_TYPE_FLOAT);
    cpl_table_fill_column_window_float(folds, colName, 0, p, 0.0);
    folded = cpl_table_get_data_float(folds, colName);
    for (j = 0; j < nx; j++) {
      folded[j % p] += line[j];
      count[j % p]++;
    }
    cpl_table_divide_columns(folds, colName, "count");
  }

  cpl_table_erase_column(folds, "count");

  /*
   *  For each folded profile i the position M_i of the minimum value
   *  is determined. Then the folded profile j containing the minimum 
   *  of all the minima is selected. The position of the first gap in
   *  the image line is given either by the position M_j of this lowest
   *  minimum value, or by the sum (M_j + P_j) (i.e. of this value with
   *  the period of the selected folded profile). To decide which one 
   *  is the value to prefer, it is checked whether the mean level of 
   *  a pixel interval before the first position is significantly above 
   *  the background level. Currently it is just checked whether the 
   *  signal is above a conventional level of 500 ADU.
   */

  stats = cpl_table_new(endp - startp + 1);
  cpl_table_new_column(stats, "Period", CPL_TYPE_INT);
  cpl_table_new_column(stats, "Minimum", CPL_TYPE_FLOAT);
  cpl_table_new_column(stats, "Position", CPL_TYPE_INT);

  for (i = 0, p = startp; p <= endp; i++, p++) {

    snprintf(colName, MAX_COLNAME, "p%d", p);
    cpl_table_get_column_minpos(folds, colName, &row);
    min = cpl_table_get_float(folds, colName, row, NULL);

    cpl_table_set_int(stats, "Period", i, p);
    cpl_table_set_int(stats, "Position", i, row);
    cpl_table_set_float(stats, "Minimum", i, min);

  }

  /*
   *  Fit the trend of the minima as a function of the folding period, 
   *  to find the best folding period. If the fit fails, take the
   *  folding period producing the profile with the deepest minimum.
   */

  data = cpl_table_get_data_float(stats, "Minimum");

  if (findDip1D(data, endp - startp + 1, &frow, 1) == VM_TRUE)
    row = frow;
  else {
    cpl_table_get_column_minpos(stats, "Minimum", &row);
    frow = row;
  }

  /*
   *  Best period, and beginning of first 80-fibers block.
   */

  period = cpl_table_get_int(stats, "Period", row, NULL);
  period += frow - row;

  cpl_msg_debug(task, "Best period: %f", period);

  gap[0] = cpl_table_get_int(stats, "Position", row, NULL);

  /*
   *  Here we try to eliminate the ambiguity inherent to the folding
   *  analysis, that gives just the phase of the gap positions along 
   *  the folding period. The position of the first gap along the 
   *  image line may coincide with the phase, or with the phase plus 
   *  or minus the folding period. To determine the right case, average
   *  the signal external to all fiber blocks (5 periods), considering
   *  the three different starting positions. The start position returning 
   *  the lowest value for the average is the right one.
   */

  n = 0;
  meanLevel1 = 0.0;
  interval = gap[0] - period;
  for (i = 0; i < interval; i++) {
    meanLevel1 += line[i];
    n++;
  }
  for (i = interval + N_BLOCKS * period; i < nx; i++) {
    meanLevel1 += line[i];
    n++;
  }

  meanLevel1 /= n;

  n = 0;
  meanLevel2 = 0.0;
  interval = gap[0];
  for (i = 0; i < interval; i++) {
    meanLevel2 += line[i];
    n++;
  }
  for (i = interval + N_BLOCKS * period; i < nx; i++) {
    meanLevel2 += line[i];
    n++;
  }

  meanLevel2 /= n;

  n = 0;
  meanLevel3 = 0.0;
  interval = gap[0] + period;
  for (i = 0; i < interval; i++) {
    meanLevel3 += line[i];
    n++;
  }
  for (i = interval + N_BLOCKS * period; i < nx; i++) {
    meanLevel3 += line[i];
    n++;
  }

  meanLevel3 /= n;

  if (meanLevel1 < meanLevel2 && meanLevel1 < meanLevel3)
    gap[0] -= period;
  else if (meanLevel2 > meanLevel3)
    gap[0] += period;

  gap[0] -= 1.0;  /* Start search position (fuzz parameter). */

  /*
   *  Now the 80 fibers of each block are "combed", i.e., they are
   *  correlated with a 80-teeth "comb" starting for the assigned
   *  position for the gaps. The number of correlation steps are
   *  defined as the folding period minus the distance between the
   *  first and the last peaks of a block. This is exactly the
   *  length where it makes sense to look for the max correlation
   *  factor without invading the correlation intervals of the other
   *  blocks.
   */

  ntrials = period - FIBERS_STEP * (FIBERS_PER_BLOCK - 1);

  for (i = 0; i < N_BLOCKS; i++) {

    /*
     *  After the first block position is determined, the start search
     *  of the next block is set at the last peak of the previous block
     *  (minus one pixel tolerance).
     */

    if (i)
      gap[i] = gap[i - 1] + FIBERS_STEP * FIBERS_PER_BLOCK - 1.0;

    /*
     *  The following code is nothing more than the search for the
     *  highest correlation factor. The ups and downs counters are
     *  just used for screening false maxima, coming from descending
     *  tails at the beginning of the interval (i.e., only maxima
     *  found after an ascending trend are considered).
     */

    n = 0;
    ups = 1;  /* Setting this to 1 (instead of 0) eliminates the screening */
    downs = 0;
    for (j = 0; j < ntrials; j++) {
      level = 0.0;
      for (k = 0; k < FIBERS_PER_BLOCK; k++) {
        pos = gap[i] + j + FIBERS_STEP * k;
        if (pos >= 0 && pos < nx)
          level += line[pos];
      }
      if (j) {
        if (level < plevel) {
          if (ups > 0) {
            max = plevel;
            downs++;
          }
          else {
            ups = 0;
          }
 
          if (downs > 0) {
            candidate[n] = max;
            peak[n] = j - 1;
            n++;
            ups = 0;
            downs = 0;
            continue;
          }
        }
        else {
          ups++;
        }
      }
      plevel = level;
    }

    row = peak[0];
    max = 0;
    for (j = 1; j < n; j++) {
      if (candidate[j] > max) {
        max = candidate[j];
        row = peak[j];
      }
    }

    /*
     *  To the start search position is added the offset where the
     *  max correlation was found.
     */

    gap[i] += row;

    cpl_msg_debug(task, "First peak of block %d: %f", i, gap[i]);

  }

  ident = cpl_table_new(N_BLOCKS * FIBERS_PER_BLOCK);
  cpl_table_new_column(ident, "Position", CPL_TYPE_INT);
  cpl_table_fill_column_window_int(ident, "Position", 0, 
                            N_BLOCKS * FIBERS_PER_BLOCK, 0);
  position = cpl_table_get_data_int(ident, "Position");

  for (k = 0, i = 0; i < N_BLOCKS; i++)
    for (j = 0; j < FIBERS_PER_BLOCK; j++, k++)
      position[k] = gap[i] + 5 * j;

  /* 
   *  Tables are written to disk just for debug purposes - keep the 
   *  destructors... 
   */

/*  cpl_table_save(folds, NULL, NULL, "folds.fits", CPL_IO_CREATE);  */
  cpl_table_delete(folds);

/*  cpl_table_save(stats, NULL, NULL, "stats.fits", CPL_IO_CREATE);  */
  cpl_table_delete(stats);

  /*
   *  Refine positions
   */

  cpl_table_name_column(ident, "Position", "FirstGuess");
  cpl_table_cast_column(ident, "FirstGuess", "Position", CPL_TYPE_FLOAT);

  for (i = 0; i < N_BLOCKS * FIBERS_PER_BLOCK; i++) {
    fpos = cpl_table_get_float(ident, "Position", i, NULL);
    if (!fiberPeak(image, refrow, &fpos, NULL))
      cpl_table_set_float(ident, "Position", i, fpos);
  }

  return ident;

}


/**
 * @memo
 *   Improve peak position.
 *
 * @return 0 on success.
 *
 * @param image    Image of an IFU flat field exposure.
 * @param row      Image row where the peak is.
 * @param position First-guess position of the peak (pixel).
 *
 * @doc
 *   In the assumption of high S/N, the peak is obtained by parabolic
 *   interpolation of the three highest pixels values within the peak.
 *
 * @author C. Izzo
 */

int fiberPeak(cpl_image *image, int row, float *position, float *max)
{

  char    task[] = "fiberPeak";

  int     nx   = cpl_image_get_size_x(image);
  int     ny   = cpl_image_get_size_y(image);
  float  *data = cpl_image_get_data(image);
  float  *line = data + row * nx;

  int     pos = *position + 0.5;    /* Nearest integer */
  int     ipos = pos;
  int     step;
  float   fpos;
  float   rpos;


  if (row < 0 || row >= ny) {
    cpl_msg_debug(task, "Image row %d out of bounds.", row);
    return 1;
  }
  
  if (!(pos > 0 && pos < nx - 1)) {
    cpl_msg_debug(task, "Peak position %f out of bounds.", *position);
    return 1;
  }

  /*
   *  Follow the gradient to find the highest peak.
   */

  if (line[pos] < line[pos - 1] && line[pos] > line[pos + 1])
    step = -1;
  else if (line[pos] > line[pos - 1] && line[pos] < line[pos + 1])
    step = 1;
  else if (line[pos] < line[pos - 1] && line[pos] < line[pos + 1])
    return 1;
  else
    step = 0;

  if (step) {
    while (line[pos] < line[pos + step]) {
      pos += step;
      if (!(pos > 0 && pos < nx - 1)) {
        cpl_msg_debug(task, "Peak position %f out of bounds.", *position);
        return 1;
      }
      if (abs(pos - ipos) > 2) {
        cpl_msg_debug(task, "Dead fiber at position %f.", *position);
        return 1;
      }
    }
  }

/*  if (line[pos] < 300) {    FIXME:  300 = too-low-signal threshold */
/*    cpl_msg_error(task, "Dead fiber at position %f.", *position);  */
/*    return 1;                                                      */
/*  }                                                                */

  /*
   * The peak position, and the value that a pixel would have at that
   * position, are determined on the basis of a mean fiber profile
   * model.
   */

  rpos = values_to_dx(line[pos - 1], line[pos], line[pos + 1]); 
  fpos = pos + rpos;

  if (fabs(*position - fpos) > 1.9)
    return 1;

  *position = fpos;

  if (max) {
    rpos = dx_to_value(rpos);
    *max = line[pos] / rpos;
  }

  return 0;
  
}


/**
 * @memo
 *   Trace all fibers.
 * 
 * @return Tables of fibers positions and fluxes along the dispersion direction.
 *
 * @param image    Image of an IFU flat field exposure.
 * @param row      Image reference row.
 * @param above    Pixels to trace above reference row.
 * @param below    Pixels to trace below reference row.
 * @param ident    Fibers identification table.
 *
 * @doc
 *   A table having a column for each fiber and as long as the spectral
 *   extraction interval (above + below + 1) is produced. The columns
 *   contain the fibers X CCD positions for each Y CCD pixel starting
 *   from (row - below). Another table containing the peak value at each
 *   computed fiber position on the same range is also produced. The fiber 
 *   identification table in input should be the product of either the 
 *   ifuIdentify() or the ifuIdentifyUpgrade() functions.
 *
 * @author C. Izzo
 */

cpl_table **ifuTrace(cpl_image *image, int row, int above, int below, 
                     int step, cpl_table *ident)
{

  char        task[] = "ifuTrace";

  cpl_table  *trace;
  cpl_table  *signal;
  cpl_table **tables;
  float      *pdata;
  float      *fdata;
  float      *mdata;
  int        *idata;
  int         ny     = cpl_image_get_size_y(image);

  int         i, j, k;
  int         range, offset;
  float       position, prePosition;
  float       max;
  char        colName[MAX_COLNAME];

  
  if (row + above > ny - 1 || row - below < 0) {
    cpl_msg_error(task, "Spectral extraction interval out of bounds.");
    return NULL;
  }

  range = above + below + 1;
  offset = row - below;

  trace = cpl_table_new(range);

  cpl_table_new_column(trace, "y", CPL_TYPE_INT);
  cpl_table_fill_column_window_int(trace, "y", 0, range, 1);
  idata = cpl_table_get_data_int(trace, "y");

  for (i = 0; i < range; i++)
    idata[i] = i;

  cpl_table_add_scalar(trace, "y", offset);

  signal = cpl_table_new(range);

  cpl_table_new_column(signal, "y", CPL_TYPE_INT);
  cpl_table_fill_column_window_int(signal, "y", 0, range, 1);
  idata = cpl_table_get_data_int(signal, "y");

  for (i = 0; i < range; i++)
    idata[i] = i;

  cpl_table_add_scalar(signal, "y", offset);

  pdata = cpl_table_get_data_float(ident, "Position");

  for (i = 0; i < N_BLOCKS * FIBERS_PER_BLOCK; i++) {
    snprintf(colName, MAX_COLNAME, "f%d", i + 1);
    cpl_table_new_column(trace, colName, CPL_TYPE_FLOAT);

    if (step > 1) {
      for (j = 0, k = row; j <= above; j += step, k += step)
        cpl_table_set_float(trace, colName, k - offset, 0.0); 
      for (j = step, k = row - step; j <= below; j += step, k -= step)
        cpl_table_set_float(trace, colName, k - offset, 0.0); 
    }
    else 
      cpl_table_fill_column_window_float(trace, colName, 0, range, 0.0);

    fdata = cpl_table_get_data_float(trace, colName);

    cpl_table_new_column(signal, colName, CPL_TYPE_FLOAT);

    if (step > 1) {
      for (j = 0, k = row; j <= above; j += step, k += step)
        cpl_table_set_float(signal, colName, k - offset, 0.0); 
      for (j = step, k = row - step; j <= below; j += step, k -= step)
        cpl_table_set_float(signal, colName, k - offset, 0.0); 
    }
    else 
      cpl_table_fill_column_window_float(signal, colName, 0, range, 0.0);

    mdata = cpl_table_get_data_float(signal, colName);
    position = pdata[i];
    for (j = 0, k = row; j <= above; j += step, k += step) {
      prePosition = position;
      if (fiberPeak(image, k, &position, &max)) {
        cpl_table_set_invalid(trace, colName, k - offset);
        cpl_table_set_invalid(signal, colName, k - offset);
      }
      else {
        if (fabs(prePosition - position) < 0.9) {    /*  WAS 0.2, poi 0.4 */
          fdata[k - offset] = position;
          mdata[k - offset] = max;
        }
        else {
          cpl_table_set_invalid(trace, colName, k - offset);
          cpl_table_set_invalid(signal, colName, k - offset);
          position = prePosition;
        }
      }
    }
    position = pdata[i];
    for (j = step, k = row - step; j <= below; j += step, k -= step) {
      prePosition = position;
      if (fiberPeak(image, k, &position, &max)) {
        cpl_table_set_invalid(trace, colName, k - offset);
        cpl_table_set_invalid(signal, colName, k - offset);
      }
      else {
        if (fabs(prePosition - position) < 0.9) {     /*  ERA 0.2, poi 0.4 */
          fdata[k - offset] = position;
          mdata[k - offset] = max;
        }
        else {
          cpl_table_set_invalid(trace, colName, k - offset);
          cpl_table_set_invalid(signal, colName, k - offset);
          position = prePosition;
        }
      }
    }
  }

/*  cpl_image_save(image, "fibers.fits", -32, NULL);  */

  tables = cpl_malloc(2 * sizeof(cpl_table *));

  tables[0] = trace;
  tables[1] = signal;

  return tables;

}


/**
 * @memo
 *   Smooth fibers intensities along the traces.
 *
 * @return 0 on success.
 *
 * @param signal      Table with all fibers signals along the traces.
 * @param order       Polynomial order.
 * @param maxNulls    Max tolerated number of NULL values in fiber signal.
 *
 * @doc
 *   A polynomial fit is made to each fiber signal contained in the
 *   table produced by the function ifuTrace(). The fit is performed
 *   only if the number of NULL points in the fiber signal is less 
 *   than @em maxNulls. The data values are replaced by the model 
 *   values, otherwise the whole column is replaced by NULLs.
 *
 * @author C. Izzo
 */

int ifuSignal(cpl_table *signal, int order, int maxNulls) 
{

  VimosDpoint *list;
  cpl_table   *oneSignal;
  float       *fdata;
  int         *idata = cpl_table_get_data_int(signal, "y");
  int          range = cpl_table_get_nrow(signal);
  int          npix;
  double      *c = NULL;

  int          rejected;
  int          i, j;
  char         colName[MAX_COLNAME];


  for (i = 0; i < N_BLOCKS * FIBERS_PER_BLOCK; i++) {

    snprintf(colName, MAX_COLNAME, "f%d", i + 1);
    rejected = cpl_table_count_invalid(signal, colName);
    if (rejected > maxNulls) {
      cpl_table_set_column_invalid(signal, colName, 0, range);
      continue;
    }

    if (rejected) {
      oneSignal = cpl_table_new(range);
      cpl_table_duplicate_column(oneSignal, "y", signal, "y");
      cpl_table_duplicate_column(oneSignal, colName, signal, colName);
      cpl_table_erase_invalid(oneSignal);
      fdata = cpl_table_get_data_float(oneSignal, colName);
      idata = cpl_table_get_data_int(oneSignal, "y");
      npix = cpl_table_get_nrow(oneSignal);
    }
    else {
      fdata = cpl_table_get_data_float(signal, colName);
      idata = cpl_table_get_data_int(signal, "y");
      npix = range;
    }

    list = newDpoint(npix);

    for (j = 0; j < npix; j++) {
      list[j].x = idata[j];
      list[j].y = fdata[j];
    }

    if (rejected)
      cpl_table_delete(oneSignal);

    c = fit1DPoly(order, list, npix, NULL);

    if (c) {
      drawModel(signal, colName, c, order);
      free(c);
      c = NULL;
    }

    deleteDpoint(list);

  }

  return 0;

}


/**
 * @memo
 *   Fit all traces.
 *
 * @return Table with fit coefficients for each fiber.
 *
 * @param trace       Table with all fiber traces.
 * @param order       Polynomial order.
 * @param tolerance   Max residual (in pixels) for point rejection.
 * @param maxReject   Number of point rejections to flag a dead fiber.
 *
 * @doc
 *   A polynomial fit is made to each fiber trace contained in the 
 *   table produced by the function ifuTrace(). The fit is performed
 *   only if the number of NULL points in the tracing is less than
 *   @em maxReject. After the fit, all points deviating from the
 *   model more than @em tolerance are also excluded. If the total
 *   number of NULLs and rejected points is less than @em maxReject
 *   the fiber is flagged as "dead", and the corresponding row of
 *   the table of coefficients is left empty (NULL).
 *
 * @author C. Izzo
 */

cpl_table **ifuFit(cpl_table *trace, int order, float tolerance, int maxReject)
{

  char task[] = "ifuFit";

  cpl_table   *coeff;
  cpl_table   *model;
  cpl_table  **tables;
  VimosDpoint *list;
  double      *c = NULL;
  double       rms;
  float       *fdata;
  int         *idata;
  int          range, npix;
  int          rejected, moreRejected;
  int          i, j, k;
  char         colName[MAX_COLNAME];


  range = cpl_table_get_nrow(trace);
  model = cpl_table_new(range);
  cpl_table_copy_structure(model, trace);
  idata = cpl_table_get_data_int(trace, "y");
  cpl_table_copy_data_int(model, "y", idata);

  coeff = cpl_table_new(N_BLOCKS * FIBERS_PER_BLOCK);
  for (i = 0; i <= order; i++) {
    snprintf(colName, MAX_COLNAME, "c%d", i);
    cpl_table_new_column(coeff, colName, CPL_TYPE_DOUBLE);
  }
  cpl_table_new_column(coeff, "rms", CPL_TYPE_DOUBLE);

  list = newDpoint(range);

  idata = cpl_table_get_data_int(trace, "y");

  for (i = 0; i < N_BLOCKS * FIBERS_PER_BLOCK; i++) {
    snprintf(colName, MAX_COLNAME, "f%d", i + 1);
    rejected = cpl_table_count_invalid(trace, colName);
    if (rejected > maxReject) {
      cpl_msg_debug(task, "Rejected fiber: %d (%d NULLs)", i + 1, rejected);
      continue;
    }

    fdata = cpl_table_get_data_float(trace, colName);

    if (rejected) {
      cpl_table_fill_invalid_float(trace, colName, -1);
      npix = 0;
      for (j = 0; j < range; j++) {
        if (fdata[j] < 0.)
          continue;
        list[npix].x = idata[j];
        list[npix].y = fdata[j];
        npix++;
      }
    }
    else {
      npix = range;
      for (j = 0; j < npix; j++) {
        list[j].x = idata[j];
        list[j].y = fdata[j];
      }
    }

    c = fit1DPoly(order, list, npix, &rms);
    if (c) {
      moreRejected = countRejections(list, npix, c, order, tolerance);
      if (rejected + moreRejected > maxReject) {
        cpl_msg_debug(task, "Rejected fiber: %d (%d bad values)", 
                   i + 1, rejected + moreRejected);
        free(c);
        c = NULL;
        continue;
      }
      else if (moreRejected) {    /*  Iteration  */
        free(c);
        c = NULL;
        c = fit1DPoly(order, list, npix - moreRejected, &rms);
      }
      if (c) {
        drawModel(model, colName, c, order);
        for (k = 0; k <= order; k++) {
          snprintf(colName, MAX_COLNAME, "c%d", k);
          cpl_table_set_double(coeff, colName, i, c[k]);
        }
        cpl_table_set_double(coeff, "rms", i, sqrt(rms));
        free(c);
        c = NULL;
      }
    }

  }

  deleteDpoint(list);

  tables = cpl_malloc(2 * sizeof(cpl_table *));

  tables[0] = coeff;
  tables[1] = model;

  return tables;

}


/*
 * @memo
 *   Fill missing traces.
 *
 * @return 0 on success.
 *
 * @param coeff       Table with coefficients of all fitted traces.
 *
 * @doc
 *   The ifuFit() function may reject a number of fibers, flagging
 *   them as "dead" as soon as the fit of the tracing fails. The
 *   reason of the failure may be a bad CCD column, or a cosmic
 *   ray with a trace about parallel to the fiber traces. In both
 *   such cases a fiber is not really "dead", and it may be used
 *   in the extraction of science spectra. For this reason it makes
 *   sense that in the table carrying the coefficients of the
 *   successful traces, all the missing solutions are restored by 
 *   interpolation. The interpolation is based on a polynomial
 *   fit of all the available coefficients, as a function of the
 *   refined start position of each fiber...  (HOLD ON...)
 *
 * @author C. Izzo
 */


/**
 * @memo
 *   Create a background image.
 *   
 * @return Background image.
 * 
 * @param image   Input image.
 * @param hr      0 = low resolution, else middle or high resolution exposure.
 * @param coeffs  Tables with coefficients of all fitted traces.
 * @param start   Start Y positions of validity ranges for traces.
 * @param end     End Y positions of validity ranges for traces.
 * @param order   Order of the bivariate polynomial for background modeling.
 * 
 * @doc
 *   An input image and all its available tracing solutions are
 *   specified. For HR and MR observations just one tracing table
 *   is available, while for LR observation four different tracing 
 *   tables must be specified, one for each IFU slit containing 
 *   400 fibers. It is unimportant the order in which such tables
 *   are given. For each table also the start and the end positions 
 *   of the tracings validity ranges must be specified.
 *   The algorithm goes this way: The tracings of the first and 
 *   the last fibers of each block are used to delimit the regions 
 *   where the background signal will be evaluated. If any of these 
 *   tracings are missing (because of lost or damaged fibers), they 
 *   are reconstructed by shifting the closest tracing to their 
 *   expected positions. The background values are taken from all 
 *   image pixels that are distant from the border tracings more 
 *   than 4 pixels. The pixel values are fitted by a low degree
 *   bivariate polynomial, that is then used to construct the
 *   background image.
 *   
 * @author C. Izzo
 */

cpl_image *ifuBack(cpl_image *image, int hr, cpl_table **coeffs, 
                   int *start, int *end, int order)
{

  float *auxildata;   
  float *data        = cpl_image_get_data(image);
  int    xlen        = cpl_image_get_size_x(image);
  int    ylen        = cpl_image_get_size_y(image);
  int    npix        = xlen * ylen;
  int    nc          = cpl_table_get_ncol(coeffs[0]);
  int    null;
  float  value;

  int    tableCount;
 
  int    limit[]     = {0, 79, 80, 159, 160, 239, 240, 319, 320, 399};
  int    limitCount  = 10;

  char   colName[MAX_COLNAME];

/*  VimosPixel  *pixel;  */
  VimosDpoint *list;
  double      *bra;       /* Tracing coefficients for left limit  */
  double      *ket;       /* Tracing coefficients for right limit */
  double      *c = NULL;  /* Generic tracing coefficients         */
  double      *b = NULL;  /* Background fit surface coefficients  */

  int dir;
  int radius = 5;         /* In number of fibers                  */
  int N = 4;              /* Distance of background from centroid */
  int x1, x2;
  int i, j, k, m, n;

  cpl_image *background;
  cpl_image *sbackground;
  cpl_image *auxil;   


  if (hr)
    tableCount = 1;
  else
    tableCount = 4;

  /*
   *  The coefficients of the last and the first fibers of one block and
   *  the next are bracketing a background region, and are therefore
   *  quantistically called "bra" and "ket". Only the tracings of the
   *  first fiber of the first block, and the last fiber of the last
   *  block, are a ket without a bra, and a bra without a ket.
   *  (This is clear, uh?).
   */

  bra = (double *)cpl_malloc(nc * sizeof(double));
  ket = (double *)cpl_malloc(nc * sizeof(double));

  n = 0;

  auxil = cpl_image_duplicate(image);   
  auxildata = cpl_image_get_data(auxil);   

  for (i = 0; i < npix; i++)   
    auxildata[i] = -1.0;   

  for (i = 0; i < tableCount; i++) {
    for (j = 0; j < limitCount; j++) {

      if (j % 2)                            /* First fiber of block (bra)   */
        c = bra;
      else                                  /* Last fiber of block (ket)    */
        c = ket;

      k = 0;

      c[0] = cpl_table_get_double(coeffs[i], "c0", limit[j], &null);

      if (null) {                           /* Lost fiber                   */

        if (c == bra)
          dir = -1;                         /* Backward alternatives        */
        else
          dir = 1;                          /* Forward alternatives         */

        /*
         *  Since a fiber trace is missing, we try to replace it with
         *  a shifted trace from a nearby fiber. To to this it is enough
         *  to shift the nearest tracing to the expected position of the
         *  lost fiber. Nearby tracings are searched within the value
         *  specified in 'radius' (in number of fibers).
         */

        for (k = 1; k < radius; k++) {
          c[0] = cpl_table_get_double(coeffs[i], "c0", 
                                      limit[j] + k * dir, &null);
          if (!null) {                      /* Alternative found            */
            c[0] -= FIBERS_STEP * k * dir;  /* The only changed coefficient */
            break;
          }
        }
      }

      if (null) {                           /* No alternative was found     */
        if (c == bra)                       /* Skip also the associated ket */
          j++;
        c = NULL;
        continue;
      }

      if (c[0] < 10 || c[0] > xlen - 10) {  /* 10 = safety marging          */
        if (c == bra)
          j++;
        c = NULL;
        continue;
      }

      for (m = 1; m < nc; m++)  {           /* Other coeffs are just copied */

        snprintf(colName, MAX_COLNAME, "c%d", m);
        c[m] = cpl_table_get_double(coeffs[i], colName, 
                                    limit[j] + k * dir, NULL);
      }

      if (c == ket) {                       /* We have closed a bracket     */
        for (k = start[i]; k < end[i]; k++) {
          if (j == 0)                       /* First fiber of first block   */
            x1 = 0;
          else
            x1 = modelValue1D(bra, nc - 1, k) + N;

          x2 = modelValue1D(ket, nc - 1, k) - N;

          if (x1 < 0)
            x1 = 0;

          if (x2 > xlen)
            x2 = xlen;

          for (m = x1; m < x2; m++, n++)
            auxildata[k * xlen + m] = data[k * xlen + m];

        }
      }
    }

    if (c == bra) {
      for (k = start[i]; k < end[i]; k++) {
        x1 = modelValue1D(bra, nc - 1, k) + N;
  
        for (m = x1; m < xlen; m++, n++)
          auxildata[k * xlen + m] = data[k * xlen + m];

      }
    }

  }

  /*
   * One fit for each image row. The alternative is to fit a bivariate 
   * polynomial (see commented part). If that is used, all code from
   * HERE to THERE should be commented out.
   */

  /*** HERE ***/

  background = cpl_image_duplicate(image);
  data = cpl_image_get_data(background);

  list = newDpoint(xlen);

  for (j = 0; j < ylen; j++) {
    for (i = 0, m = 0; i < xlen; i++) {
      value = auxildata[j * xlen + i];
      if (value > 0.0) {
        list[m].x = i;
        list[m].y = value;
        m++;
      }
    }

    if (m > order + 1) {
      b = fit1DPoly(order, list, m, NULL);
      if (b) {
        for (i = 0; i < xlen; i++)
          data[j * xlen + i] = modelValue1D(b, order, i);
        free(b);
        b = NULL;
      }
    }

  }

  deleteDpoint(list);

  sbackground = cpl_image_general_median_filter(background, 1, 15, 0);

  cpl_image_delete(background);

/*** THERE ***/

/*** This part (up to FINIS) does the bivariate fitting:

  pixel = newPixel(n);

  for (j = 0, m = 0; j < ylen; j++) {
    for (i = 0; i < xlen; i++) {
      value = auxildata[j * xlen + i];
      if (value > 0.0) {
        pixel[m].x = i;
        pixel[m].y = j;
        pixel[m].i = value;
        m++;
      }
    }
  }

  cpl_image_save(auxil, "debug.fits", -32, NULL);   
  cpl_image_delete(auxil);   

  b = fitSurfacePolynomial(pixel, n, NULL, order, &n, NULL);

  deletePixel(pixel);

  background = cpl_image_duplicate(image);
  data = cpl_image_get_data(background);

  for (j = 0; j < ylen; j++)
    for (i = 0; i < xlen; i++)
      data[j * xlen + i] = modelValue2D(b, order, i, j);

 *** FINIS ***/

  cpl_image_delete(auxil);   
  return sbackground;

}


/**
 * @memo
 *   Model gap background.
 *   
 * @return Background table.
 * 
 * @param image   Input image.
 * @param coeffs  Table with coefficients of all fitted traces.
 * @param start   Start Y position of validity range for traces.
 * @param end     End Y position of validity range for traces.
 * @param sbox    Smooth box length (in pixel).
 * 
 * @doc
 *   An input image and the available tracing solutions for a given
 *   range (corresponding to one IFU slit) are specified. The start 
 *   and the end positions of the tracings validity ranges must be 
 *   specified. The algorithm goes this way: The tracings of the first 
 *   and the last fibers of each block are used to delimit the regions 
 *   where the background signal will be evaluated. If any of these 
 *   tracings are missing (because of lost or damaged fibers), they 
 *   are reconstructed by shifting the closest tracing to their 
 *   expected positions. The mean value of the pixels between the
 *   the two tracings of a gap (excluding a margin of N pixels from
 *   the tracings) at a given Y position is stored in an output table 
 *   column that is then median filtered. The table columns are named 
 *   b0 to b5, where b1 to b4 refers to the gaps between fibers blocks, 
 *   b0 to the background before the first block, and b5 to the background 
 *   after the last block. Some of the columns may be missing, either 
 *   because the background layed outside the image boundaries, or 
 *   because the appropriate tracings could not be found.
 *   
 * @author C. Izzo
 */

cpl_table *ifuGap(cpl_image *image, cpl_table *coeffs, 
                  int start, int end, int sbox)
{

  cpl_table *table;
  float     *data        = cpl_image_get_data(image);
  int        xlen        = cpl_image_get_size_x(image);
  int        nc          = cpl_table_get_ncol(coeffs);
  int        null;
  double     mean;
 
  int        limit[]     = {0, 79, 80, 159, 160, 239, 240, 319, 320, 399};
  int        limitCount  = 10;

  char       colName[MAX_COLNAME];

  double     *bra;       /* Tracing coefficients for left limit  */
  double     *ket;       /* Tracing coefficients for right limit */
  double     *c;         /* Generic tracing coefficients         */

  int dir;
  int radius = 5;        /* Given in number of fibers            */
  int N = 4;             /* Distance of background from centroid */
  int x1, x2;
  int j, k, m, n;
  int nrow;


  /*
   *  The coefficients of the last and the first fibers of one block and
   *  the next are bracketing a background region, and are therefore
   *  quantistically called "bra" and "ket". Only the tracings of the
   *  first fiber of the first block, and the last fiber of the last
   *  block, are a ket without a bra, and a bra without a ket.
   *  (This is clear, uh?).
   */

  bra = (double *)cpl_malloc(nc * sizeof(double));
  ket = (double *)cpl_malloc(nc * sizeof(double));

  nrow = end - start;
  table = cpl_table_new(nrow);
  cpl_table_new_column(table, "y", CPL_TYPE_INT);

  for (j = start; j < end; j++)
    cpl_table_set_int(table, "y", j - start, j);

  for (j = 0; j < limitCount; j++) {

    if (j % 2)                            /* First fiber of block (bra)   */
      c = bra;
    else                                  /* Last fiber of block (ket)    */
      c = ket;

    k = 0;

    c[0] = cpl_table_get_double(coeffs, "c0", limit[j], &null);

    if (null) {                           /* Lost fiber                   */

      if (c == bra)
        dir = -1;                         /* Backward alternatives        */
      else
        dir = 1;                          /* Forward alternatives         */

      /*
       *  Since a fiber trace is missing, we try to replace it with
       *  a shifted trace from a nearby fiber. To to this it is enough
       *  to shift the nearest tracing to the expected position of the
       *  lost fiber. Nearby tracings are searched within the value
       *  specified in 'radius' (in number of fibers).
       */

      for (k = 1; k < radius; k++) {
        c[0] = cpl_table_get_double(coeffs, "c0", limit[j] + k * dir, &null);
        if (!null) {                      /* Alternative found            */
          c[0] -= FIBERS_STEP * k * dir;  /* The only changed coefficient */
          break;
        }
      }
    }

    if (null) {                           /* No alternative was found     */
      if (c == bra)                       /* Skip also the associated ket */
        j++;
      c = NULL;
      continue;
    }

    if (c[0] < 10 || c[0] > xlen - 10) {  /* 10 = safety marging          */
      if (c == bra)
        j++;
      c = NULL;
      continue;
    }

    for (m = 1; m < nc; m++)  {           /* Other coeffs are just copied */

      snprintf(colName, MAX_COLNAME, "c%d", m);
      c[m] = cpl_table_get_double(coeffs, colName, limit[j] + k * dir, NULL);
    }

    if (c == ket) {                       /* We have closed a bracket     */

      snprintf(colName, MAX_COLNAME, "b%d", j / 2);
      cpl_table_new_column(table, colName, CPL_TYPE_FLOAT);

      for (k = start; k < end; k++) {

        x2 = modelValue1D(ket, nc - 1, k) - N;

        if (j == 0)                       /* First fiber of first block   */
          x1 = x2 - 10;
        else
          x1 = modelValue1D(bra, nc - 1, k) + N;

        if (x1 < 0)
          x1 = 0;

        if (x2 > xlen)
          x2 = xlen;

        mean = 0.0;
        for (m = x1, n = 0; m < x2; m++, n++)
          mean += data[k * xlen + m];

/* To exclude max value from mean: eliminated, because leads to background
   underestimation.

        if (n > 1) {
          max = data[k * xlen + x1];
          for (m = x1; m < x2; m++)
            if (data[k * xlen + m] > max)
              max = data[k * xlen + m];
          mean -= max;
          --n;
        }
*/

        mean /= n;
  
        cpl_table_set_float(table, colName, k - start, mean);

      }

      cpl_table_median_filter_column(table, colName, sbox / 2);

    }
  }

  if (c == bra) {

    snprintf(colName, MAX_COLNAME, "b%d", 5);
    cpl_table_new_column(table, colName, CPL_TYPE_FLOAT);

    for (k = start; k < end; k++) {
      x1 = modelValue1D(bra, nc - 1, k) + N;
      x2 = x1 + 10;

      if (x2 > xlen)
        x2 = xlen;

      mean = 0.0;
      for (m = x1, n = 0; m < x2; m++, n++)
        mean += data[k * xlen + m];

/* To exclude max value from mean: eliminated, because leads to background
   underestimation.

      if (n > 1) {
        max = data[k * xlen + x1];
        for (m = x1; m < x2; m++)
          if (data[k * xlen + m] > max)
            max = data[k * xlen + m];
        mean -= max;
        --n;
      }
*/

      mean /= n;
  
      cpl_table_set_float(table, colName, k - start, mean);

    }

    cpl_table_median_filter_column(table, colName, sbox / 2);

  }

  return table;

}


/**
 * @memo
 *   Produce table of samplings of the fiber profile.
 *      
 * @return Table of samplings of the fiber profile.
 *
 * @param image   Input image.
 * @param model   Table with all fitted traces.
 * @param signal  Table with all fibers signal along the traces.
 * @param backs   Table with background values along each fiber block gap.
 * 
 * @doc
 *   The 10 half profiles of the 10 spectra siding the background
 *   regions are used to construct a table consisting of pixel
 *   values and their distances from the profile centroid (derived
 *   from the tracing solution), at each Y CCD coordinate within
 *   the specified range. The profile maximum contained in the
 *   signal table and the pixel values are corrected for the local 
 *   background contained in table backs, produced by the function
 *   ifuGap()). The pixel values are then normalized to the profile 
 *   maximum. 
 *
 * @author C. Izzo
 */

cpl_table *ifuProfile(cpl_image *image, cpl_table *model, cpl_table *signal,
                      cpl_table *backs)
{

  char task[] = "ifuProfile";

  cpl_table *table;

  float     *data   = cpl_image_get_data(image);
  int        xlen   = cpl_image_get_size_x(image);

  int       *idata  = cpl_table_get_data_int(model, "y");
  int       *iodata;                              /* Y coordinate on image  */
  float     *mdata;                               /* Centroid positions     */
  float     *sdata;                               /* Normalization factor   */
  float     *odata;                               /* Output pixel values    */
  float     *cdata;                               /* Distance from centroid */
  float     *bdata;                               /* Background data        */
  int        nrow   = cpl_table_get_nrow(model);

  int        i, j, k, m;
  int        x, y;
  float      fx;
  float      max;

  int        nsamples = 6;    /*  Points for each PSF profile  */

  int        limit[]     = {1, 80, 81, 160, 161, 240, 241, 320, 321, 400};
  int        limitCount  = 10;
  int        dir;

  double     (*myfloor)(double);

  char       colName[MAX_COLNAME];
  char       bakName[MAX_COLNAME];


  table = cpl_table_new(nrow * nsamples);
  cpl_table_new_column(table, "y", CPL_TYPE_INT);
  cpl_table_fill_column_window_int(table, "y", 0, nrow * nsamples, 0.0);
  iodata = cpl_table_get_data_int(table, "y");
  for (j = 0, m = 0; j < nrow; j++)
    for (k = 0; k < nsamples; k++, m++)
      iodata[m] = idata[j];

  for (i = 0; i < limitCount; i++) {

    if (i % 2) {
      dir = 1;
      myfloor = floor;
    }
    else {
      dir = -1;
      myfloor = ceil;
    }

    snprintf(colName, MAX_COLNAME, "f%d", limit[i]);

    if (cpl_table_has_invalid(model, colName)) {
      cpl_msg_debug(task, "Cannot build profile of fiber %d", limit[i]);
      continue;
    }

    /*
     *  The background values
     */

    snprintf(bakName, MAX_COLNAME, "b%d", (i + 1) / 2);
    bdata = cpl_table_get_data_float(backs, bakName);

    if (!bdata) {
      cpl_msg_debug(task, "Cannot build profile of fiber %d", limit[i]);
      continue;
    }

    /*
     *  The signal
     */

    cpl_table_fill_invalid_float(signal, colName, -1.0);
    sdata = cpl_table_get_data_float(signal, colName);

    /*
     *  The modeled traces
     */

    mdata = cpl_table_get_data_float(model, colName);

    /*
     *  The output data
     */

    cpl_table_new_column(table, colName, CPL_TYPE_INT);
    cpl_table_fill_column_window_float(table, colName, 0, nrow * nsamples, 0.0);
    odata = cpl_table_get_data_float(table, colName);

    /*
     *  The distances from the centroid
     */

    snprintf(colName, MAX_COLNAME, "d%d", limit[i]);
    cpl_table_new_column(table, colName, CPL_TYPE_FLOAT);
    cpl_table_fill_column_window_float(table, colName, 0, nrow * nsamples, 0.0);
    cdata = cpl_table_get_data_float(table, colName);

    /*
     *  Writing the result
     */

    for (j = 0, m = 0; j < nrow; j++) {
      y = idata[j];
      fx = mdata[j];
      x = myfloor(fx);
      max = sdata[j] - bdata[j];
      for (k = 0; k < nsamples; k++, x += dir, m++) {
        if (x > 0 && x < xlen && max > 0.0) {
          odata[m] = (data[y * xlen + x] - bdata[j]) / max;
          cdata[m] = fabs(x - fx);
        }
        else
          cpl_table_set_invalid(table, colName, m);
      }
    }
  }

  if (cpl_table_get_ncol(table) > 1)
    return table;

  cpl_msg_warning(task, "Table of fiber profiles not created!");
  cpl_table_delete(table);

  return NULL;

}


/**
 * @memo
 *   Build model profile table from table created by ifuProfile().
 *
 * @return Model profile table.
 *
 * @param profiles   Table created by ifuProfile().
 * @param start      Start Y position to select profile points (inclusive).
 * @param end        End Y position to select profile points (exlcusive).
 * @param length     Length of rebin interval (in pixel).
 * @param bin        Bin size (in pixel).
 *
 * @doc
 *   Within the specified range along the dispersion direction the
 *   pixel values are averaged. The obtained values are written to
 *   a table having the same column names of the input table (but 
 *   with column y removed). The "distance" column contains the
 *   distance of the center of each bin from the centroid of the 
 *   fiber profile. The number of rows in the output table is the
 *   integer part of length/bin.
 *
 * @author C. Izzo
 */

cpl_table *rebinProfile(cpl_table *profiles, int start, int end, 
                        double length, double bin)
{

  char       task[]  = "rebinProfile";
  int        nbin       = length / bin;
  cpl_table *table      = cpl_table_new(nbin);
  cpl_table *selected;

  double     xvalue, yvalue;
  double    *buffer;
  int       *count;
  int        null;

  int        limit[]    = {1, 80, 81, 160, 161, 240, 241, 320, 321, 400};
  int        limitCount = 10;

  char       distance[MAX_COLNAME];
  char       flux[MAX_COLNAME];

  int        nrow;
  int        i, j, pos;


  cpl_table_copy_structure(table, profiles);


  /*
   *  Extract just the indicated interval along the dispersion direction
   */

  cpl_table_and_selected_int(profiles, "y", CPL_NOT_LESS_THAN, start);
  nrow = cpl_table_and_selected_int(profiles, "y", CPL_LESS_THAN, end);

  selected = cpl_table_extract_selected(profiles);
  cpl_table_select_all(profiles);

  cpl_table_erase_column(table, "y");

  /*
   *  Initialize the "distance" column of the rebinned profiles.
   *  It contains the distance of the midpoint of each bin from the
   *  fiber profile centroid.
   */

  cpl_table_new_column(table, "distance", CPL_TYPE_FLOAT);
  for (j = 0; j < nbin; j++)
    cpl_table_set_float(table, "distance", j, bin * (j + 0.5));

  /*
   *  Auxiliary buffers
   */

  buffer = cpl_malloc(nbin * sizeof(double));
  count = cpl_malloc(nbin * sizeof(int));

  /*
   *  Loop on available profiles.
   */

  for (i = 0; i < limitCount; i++) {
  
    snprintf(distance, MAX_COLNAME, "d%d", limit[i]);
    snprintf(flux, MAX_COLNAME, "f%d", limit[i]);
  
    cpl_error_reset();
  
    if (!cpl_table_has_valid(selected, distance)) {
      if (CPL_ERROR_DATA_NOT_FOUND == cpl_error_get_code()) 
        cpl_msg_debug(task, "Missing profile for fiber %d", limit[i]);
      else
        cpl_msg_debug(task, "Cannot rebin profile of fiber %d", limit[i]);
      continue;
    }    

    /*
     *  Cleaning from output table useless columns as we go...
     */

    cpl_table_erase_column(table, distance);

    for (j = 0; j < nbin; j++) {
      buffer[j] = 0.0;
      count[j] = 0;
    }
  
    for (j = 0; j < nrow; j++) {
      xvalue = cpl_table_get_float(selected, distance, j, &null);
      yvalue = cpl_table_get_float(selected, flux, j, NULL);
      if (!null) {
        pos = floor(xvalue / bin);
        if (pos < nbin) {
          buffer[pos] += yvalue;
          count[pos]++;
        }
      }
    }

    for (j = 0; j < nbin; j++)
      if (count[j] > 0)
        cpl_table_set_float(table, flux, j, buffer[j] / count[j]);

    /*
     * Elements left at NULL are linearly interpolated:
     */

    /* %%%%%%%%%%%%%%%%%%%%%%%%%% DA FARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  
  }

  cpl_free(buffer);
  cpl_free(count);

  return table;

}


/**
 * @memo
 *   Fit a gaussian model to IFU fiber empirical profiles.
 *
 * @return Table with model coefficients (scale, mean, sigma).
 *
 * @param profiles   Table created by ifuProfile().
 * @param start      Start Y position to select profile points (inclusive).
 * @param end        End Y position to select profile points (excusive).
 *
 * @doc
 *   Within the specified range along the dispersion direction the
 *   pixel values are used for a gaussian fit. The model coefficients
 *   when available are stored at a row of the output table. In theory
 *   the scale should always be 1.0 and the mean 0.0, so any deviation
 *   from these values indicates a bad fit.
 *
 * @author C. Izzo
 */

cpl_table *ifuGauss(cpl_table *table, int start, int end)
{

  char task[] = "ifuGauss";

  cpl_table       *models;
  cpl_table       *selected;
  VimosFloatArray *x;
  VimosFloatArray *y;
  float            coeff[3];

  int    limit[]     = {1, 80, 81, 160, 161, 240, 241, 320, 321, 400};
  int    limitCount  = 10;

  int    minCount = 100;   /* Minimum number of points for fit */
  int    count;

  char   distance[MAX_COLNAME];
  char   flux[MAX_COLNAME];

  int    nrow;
  int    null;
  int    i, j;

  float  xvalue, yvalue;


  cpl_table_and_selected_int(table, "y", CPL_NOT_LESS_THAN, start);
  nrow = cpl_table_and_selected_int(table, "y", CPL_LESS_THAN, end);

  if (nrow < minCount)
    return NULL;

  selected = cpl_table_extract_selected(table);
  cpl_table_select_all(table);

  models = cpl_table_new(limitCount);
  cpl_table_new_column(models, "max", CPL_TYPE_FLOAT);
  cpl_table_new_column(models, "mean", CPL_TYPE_FLOAT);
  cpl_table_new_column(models, "sigma", CPL_TYPE_FLOAT);

  x = newFloatArray(nrow);
  y = newFloatArray(nrow);

  for (i = 0; i < limitCount; i++) {

    snprintf(distance, MAX_COLNAME, "d%d", limit[i]);
    snprintf(flux, MAX_COLNAME, "f%d", limit[i]);

    cpl_error_reset();

    if (!cpl_table_has_valid(selected, distance)) {
      cpl_msg_debug(task, "Cannot fit profile of fiber %d", limit[i]);
      continue;
    }

    if (CPL_ERROR_DATA_NOT_FOUND == cpl_error_get_code()) {
      cpl_msg_debug(task, "Missing profile for fiber %d", limit[i]);
      continue;
    }

    count = cpl_table_count_invalid(selected, distance);
    count = nrow - count;

    if (count < minCount)
      continue;

    x->len = y->len = count;  /* Use just a part of the allocated array */
    
    count = 0;
    for (j = 0; j < nrow; j++) {
      xvalue = cpl_table_get_float(selected, distance, j, &null);
      yvalue = cpl_table_get_float(selected, flux, j, NULL);
      if (!null) {
        x->data[count] = xvalue;
        y->data[count] = yvalue;
        count++;
      }
    }

    fit1DGauss(x, y, coeff, 3);

    cpl_msg_debug(task, "Profile %d: max = %f, mean = %f, sigma = %f",
               i, coeff[0], coeff[1], coeff[2]);

    cpl_table_set_float(models, "max", i, coeff[0]);
    cpl_table_set_float(models, "mean", i, coeff[1]);
    cpl_table_set_float(models, "sigma", i, coeff[2]);

  }

  deleteFloatArray(x);
  deleteFloatArray(y);

  cpl_table_delete(selected);

  return models;

}


/**
 * @memo
 *   Fit a gaussian model to IFU fiber empirical profiles.
 *
 * @return Table with model coefficients (scale, mean, sigma).
 *
 * @param profiles   Table created by ifuProfile().
 * @param start      Start Y position to select profile points (inclusive).
 * @param end        End Y position to select profile points (excusive).
 *
 * @doc
 *   Within the specified range along the dispersion direction the
 *   pixel values are used for a gaussian fit. The model coefficients
 *   when available are stored at a row of the output table. In theory
 *   the scale should always be 1.0 and the mean 0.0, so any deviation
 *   from these values indicates a bad fit.
 *
 * @author C. Izzo
 */

cpl_table *ifuGauss2(cpl_table *table, int start, int end)
{

  char task[] = "ifuGauss";

  cpl_table       *models;
  cpl_table       *selected;
  VimosFloatArray *x;
  VimosFloatArray *y;
  float            coeff[3];

  int    limit[]     = {1, 80, 81, 160, 161, 240, 241, 320, 321, 400};
  int    limitCount  = 10;

  int    minCount = 100;   /* Minimum number of points for fit */
  int    count;

  char   distance[MAX_COLNAME];
  char   flux[MAX_COLNAME];

  int    nrow;
  int    null;
  int    i, j;

  float  xvalue, yvalue;


  cpl_table_and_selected_int(table, "y", CPL_NOT_LESS_THAN, start);
  nrow = cpl_table_and_selected_int(table, "y", CPL_LESS_THAN, end);

  if (nrow < minCount)
    return NULL;

  selected = cpl_table_extract_selected(table);
  cpl_table_select_all(table);

  models = cpl_table_new(limitCount);
  cpl_table_new_column(models, "max", CPL_TYPE_FLOAT);
  cpl_table_new_column(models, "mean", CPL_TYPE_FLOAT);
  cpl_table_new_column(models, "sigma", CPL_TYPE_FLOAT);

  x = newFloatArray(2 * nrow);
  y = newFloatArray(2 * nrow);

  for (i = 0; i < limitCount; i++) {

    snprintf(distance, MAX_COLNAME, "d%d", limit[i]);
    snprintf(flux, MAX_COLNAME, "f%d", limit[i]);

    cpl_error_reset();

    if (!cpl_table_has_valid(selected, distance)) {
      cpl_msg_debug(task, "Cannot fit profile of fiber %d", limit[i]);
      continue;
    }

    if (CPL_ERROR_DATA_NOT_FOUND == cpl_error_get_code()) {
      cpl_msg_debug(task, "Missing profile for fiber %d", limit[i]);
      continue;
    }

    count = cpl_table_count_invalid(selected, distance);
    count = nrow - count;

    if (count < minCount)
      continue;

    x->len = y->len = 2 * count;  /* Use just a part of the allocated array */
    
    count = 0;
    for (j = 0; j < nrow; j++) {
      xvalue = cpl_table_get_float(selected, distance, j, &null);
      yvalue = cpl_table_get_float(selected, flux, j, NULL);
      if (!null) {
        x->data[count] = xvalue;
        y->data[count] = yvalue;
        count++;
        x->data[count] = -xvalue;
        y->data[count] = yvalue;
        count++;
      }
    }

    fit1DGauss(x, y, coeff, 3);

    cpl_msg_debug(task, "Profile %d: max = %f, mean = %f, sigma = %f",
               i, coeff[0], coeff[1], coeff[2]);

    cpl_table_set_float(models, "max", i, coeff[0]);
    cpl_table_set_float(models, "mean", i, coeff[1]);
    cpl_table_set_float(models, "sigma", i, coeff[2]);

  }

  deleteFloatArray(x);
  deleteFloatArray(y);

  cpl_table_delete(selected);

  return models;

}


/**
 * @memo 
 *   Detect traceable spectra on image.
 *
 * @return Table with cross-dispersion positions of traceable spectra.
 *
 * @param image     Bias subtracted image containing IFU spectra.
 * @param row       Image row where to run peak detection.
 * @param minSignal Min peak signal above mean minimum.
 *  
 * @doc
 *   Along the specified image row significant signal peaks are detected
 *   and analysed. A peak is considered significant when it rises above 
 *   the background level more than a given amount (currently, 70 ADU),
 *   and it is a sequence of 3 increasing values and 4 decreasing, or 
 *   a sequence of 4 increasing values and 3 decreasing. This is a secure
 *   protection against the effects of bad pixels, hot CCD columns and
 *   cosmic rays. The backgound level is taken as the mean of the output 
 *   of a min-filter with running box size of 7 pixels applied to the 
 *   image row data.
 *
 * @author C. Izzo
 */
  
cpl_table *ifuDetect(cpl_image *image, int row, float minSignal)
{

  int        nx     = cpl_image_get_size_x(image);
  float     *data   = cpl_image_get_data(image);
  float     *line   = data + row * nx;
  float     *fdata;
  int       *idata;
  int       *mask1;
  int       *mask2;

  cpl_table *imageRow;
  cpl_table *ident;

  char       colXpos[]   = "x";
  char       colSvalue[] = "svalue";
  char       colValue[]  = "value";
  float      mean, fpos;
  /* float      minSignal   = 70.0;  Removed by Peter Weilbacher */
  int        hw          = 3;
  int        count, countTraceable;
  int        i;

  cpl_propertylist *reflist;

  /*
   *  Create a table containing the selected image row data.
   *  A light median filter is passed to eliminate possible hot 
   *  columns. One column contains the x positions and another 
   *  the values of the corresponding pixels.
   */

  imageRow = cpl_table_new(nx);

  cpl_table_new_column(imageRow, colValue, CPL_TYPE_FLOAT);
  cpl_table_copy_data_float(imageRow, colValue, line);

  /*
   *  Originally a median filtering of the image row was added
   *  here, to eliminate contributes from hot CCD columns. This 
   *  was not a good idea, because it caused all peaks to have
   *  flat maxima, that introduced a bias on the peak positions.
   *  The call was:
   *
   *    cpl_table_median_filter_column(imageRow, colValue, 1);
   */

  cpl_table_new_column(imageRow, colXpos, CPL_TYPE_INT);
  idata = cpl_table_get_data_int(imageRow, colXpos);
  cpl_table_fill_column_window_int(imageRow, colXpos, 0, nx, 0);
  for (i = 0; i < nx; i++)
    idata[i] = i;

  /*
   *  Apply the min filter to estimate the mean backgound level.
   */

  cpl_table_min_filter_column(imageRow, colValue, colSvalue, hw);
  mean = cpl_table_get_column_mean(imageRow, colSvalue);
  cpl_table_subtract_scalar(imageRow, colValue, mean);
  cpl_table_erase_column(imageRow, colSvalue);

  /*
   *  Sort table according to pixel values, highest values on top.
   */

  reflist = cpl_propertylist_new();
  cpl_propertylist_append_bool(reflist, colValue, TRUE);
  cpl_table_sort(imageRow, reflist);
  cpl_propertylist_delete(reflist);

  idata = cpl_table_get_data_int(imageRow, colXpos);
  fdata = cpl_table_get_data_float(imageRow, colValue);

  /*
   *  Apply search method (first part, filling masks).
   */

  mask1 = cpl_calloc(nx, sizeof(int));
  mask2 = cpl_calloc(nx, sizeof(int));

  for (i = 0; i < nx; i++) {
    if (fdata[i] < minSignal)
      break;
    if (idata[i] > hw && idata[i] < nx - hw) {
      mask2[idata[i]] = 1;
      if (mask2[idata[i] - 1] == 0 && mask2[idata[i] + 1] == 0) 
        mask1[idata[i]] = 1;
    }
  }

  cpl_table_delete(imageRow);
  cpl_free(mask2);

  /*
   *  Search method (second part, selecting peaks having the right size).
   *  Among all the local maxima found, select just those that are the
   *  tip of a 5 pixels wide bell profile - i.e., if the peak is at
   *  position N and its value is S(N), then it must be:
   *
   *                S(N) > S(N+1) > S(N+2)
   *  and
   *                S(N) > S(N-1) > S(N-2)
   *
   *  Moreover, it must be either S(N+2) > S(N+3) or S(N-2) > S(N-3).
   *  It is not necessary to repeat in the code the comparison of S(N)
   *  with its immediate two neighbours, because this is ensured by the
   *  previous search.
   */

  countTraceable = 0;

  for (i = 0; i < nx; i++) {
    if (mask1[i]) {
      mask1[i] = 0;
      if (line[i + 1] > line[i + 2]) {
        if (line[i - 1] > line[i - 2]) {
          if (line[i + 2] > line[i + 3] || line[i - 2] > line[i - 3]) {
            countTraceable++;
            mask1[i] = 1;
          }
        }
      }
    }
  }

  if (countTraceable == 0) {
    cpl_free(mask1);
    return NULL;
  }

  ident = cpl_table_new(countTraceable);
  cpl_table_new_column(ident, "Position", CPL_TYPE_INT);
  for (count = 0, i = 0; i < nx; i++) {
    if (mask1[i] == 1) {
      cpl_table_set_int(ident, "Position", count, i);
      count++;
    }
  }

  cpl_free(mask1);

  /*
   *  Refine positions
   */

  cpl_table_name_column(ident, "Position", "FirstGuess");
  cpl_table_cast_column(ident, "FirstGuess", "Position", CPL_TYPE_FLOAT);

  for (i = 0; i < countTraceable; i++) {
    fpos = cpl_table_get_float(ident, "Position", i, NULL);
    if (!fiberPeak(image, row, &fpos, NULL))
      cpl_table_set_float(ident, "Position", i, fpos);
  }

  return ident;

}


/**
 * @memo
 *   Trace all detected fibers.
 * 
 * @return Table of fibers positions along the dispersion direction.
 *
 * @param image    Image of an IFU exposure (science).
 * @param row      Image reference row.
 * @param above    Pixels to trace above reference row.
 * @param below    Pixels to trace below reference row.
 * @param ident    Fibers detection table.
 *
 * @doc
 *   This function is pretty similar to the ifuTrace(), with the
 *   difference that it tries to trace just starting from positions
 *   of detected spectra. In practice, this function is used to
 *   trace spectra on science frames. A table having a column for 
 *   each fiber and as long as the spectral extraction interval 
 *   (above + below + 1) is produced. The columns contain the 
 *   fibers X CCD positions for each Y CCD pixel starting from 
 *   (row - below). The detected spectra table in input should 
 *   be the product of the ifuDetect() function.
 *
 * @author C. Izzo
 */

cpl_table *ifuTraceDetected(cpl_image *image, int row, int above, int below, 
                            int step, cpl_table *ident)
{

  char        task[] = "ifuTraceDetected";

  cpl_table  *trace;
  float      *pdata;
  float      *fdata;
  int        *idata;
  int         ny = cpl_image_get_size_y(image);
  int         countTraceable = cpl_table_get_nrow(ident);

  int         i, j, k;
  int         range, offset;
  float       position, prePosition;
  float       max;
  char        colName[MAX_COLNAME];

  
  if (row + above > ny - 1 || row - below < 0) {
    cpl_msg_error(task, "Spectral extraction interval out of bounds.");
    return NULL;
  }

  range = above + below + 1;
  offset = row - below;

  trace = cpl_table_new(range);

  cpl_table_new_column(trace, "y", CPL_TYPE_INT);
  cpl_table_fill_column_window_int(trace, "y", 0, range, 1);
  idata = cpl_table_get_data_int(trace, "y");

  for (i = 0; i < range; i++)
    idata[i] = i;

  cpl_table_add_scalar(trace, "y", offset);

  pdata = cpl_table_get_data_float(ident, "Position");

  for (i = 0; i < countTraceable; i++) {
    snprintf(colName, MAX_COLNAME, "t%d", i + 1);
    cpl_table_new_column(trace, colName, CPL_TYPE_FLOAT);

    if (step > 1) {
      for (j = 0, k = row; j <= above; j += step, k += step)
        cpl_table_set_float(trace, colName, k - offset, 0.0);
      for (j = step, k = row - step; j <= below; j += step, k -= step)
        cpl_table_set_float(trace, colName, k - offset, 0.0);
    }
    else
      cpl_table_fill_column_window_float(trace, colName, 0, range, 0.0);

    fdata = cpl_table_get_data_float(trace, colName);

    position = pdata[i];
    for (j = 0, k = row; j <= above; j += step, k += step) {
      prePosition = position;
      if (fiberPeak(image, k, &position, &max)) {
        cpl_table_set_invalid(trace, colName, k - offset);
      }
      else {
        if (fabs(prePosition - position) < 0.4) {    /*  ERA 0.2 */
          fdata[k - offset] = position;
        }
        else {
          cpl_table_set_invalid(trace, colName, k - offset);
          position = prePosition;
        }
      }
    }
    position = pdata[i];
    for (j = step, k = row - step; j <= below; j += step, k -= step) {
      prePosition = position;
      if (fiberPeak(image, k, &position, &max)) {
        cpl_table_set_invalid(trace, colName, k - offset);
      }
      else {
        if (fabs(prePosition - position) < 0.4) {     /*  ERA 0.2 */
          fdata[k - offset] = position;
        }
        else {
          cpl_table_set_invalid(trace, colName, k - offset);
          position = prePosition;
        }
      }
    }
  }

  return trace;

}


/**
 * @memo
 *   Fit all fiber tracings on a science exposure.
 *
 * @return Table with fit coefficients for each fiber.
 *
 * @param trace       Table with detected fiber traces ifuTraceDetected().
 * @param order       Polynomial order.
 * @param tolerance   Max residual (in pixels) for point rejection.
 * @param maxReject   Number of point rejections to flag a dead fiber.
 *
 * @doc
 *   This function is pretty similar to the ifuFit(), with the
 *   difference that it tries to fit the unidentified spectra
 *   traced on science frames. A polynomial fit is made to each 
 *   fiber trace contained in the table produced by the function 
 *   ifuTraceDetected(). 
 *
 * @author C. Izzo
 */

cpl_table *ifuFitDetected(cpl_table *trace, int order, 
                          float tolerance, int maxReject)
{

  char task[] = "ifuFitDetected";

  cpl_table   *oneTrace;
  cpl_table   *coeff;
  VimosDpoint *list;
  double      *c = NULL;
  float       *fdata;
  int         *idata;
  int          range, npix;
  int          rejected, moreRejected;
  int          countTraced = cpl_table_get_ncol(trace) - 1;
  int          i, j, k;
  char         colName[MAX_COLNAME];


  coeff = cpl_table_new(countTraced);
  for (i = 0; i <= order; i++) {
    snprintf(colName, MAX_COLNAME, "c%d", i);
    cpl_table_new_column(coeff, colName, CPL_TYPE_DOUBLE);
  }

  range = cpl_table_get_nrow(trace);
  list = newDpoint(range);

  for (i = 0; i < countTraced; i++) {
    snprintf(colName, MAX_COLNAME, "t%d", i + 1);
    rejected = cpl_table_count_invalid(trace, colName);
    if (rejected > maxReject) {
      cpl_msg_debug(task, "Rejected fiber: %d (%d NULLs)", i + 1, rejected);
      continue;
    }

    if (rejected) {
      oneTrace = cpl_table_new(range);
      cpl_table_duplicate_column(oneTrace, "y", trace, "y");
      cpl_table_duplicate_column(oneTrace, colName, trace, colName);
      cpl_table_erase_invalid(oneTrace);
      fdata = cpl_table_get_data_float(oneTrace, colName);
      idata = cpl_table_get_data_int(oneTrace, "y");
      npix = cpl_table_get_nrow(oneTrace);
    }
    else {
      fdata = cpl_table_get_data_float(trace, colName);
      idata = cpl_table_get_data_int(trace, "y");
      npix = range;
    }

    for (j = 0; j < npix; j++) {
      list[j].x = idata[j];
      list[j].y = fdata[j];
    }

    if (rejected)
      cpl_table_delete(oneTrace);

    c = fit1DPoly(order, list, npix, NULL);
    if (c) {
      moreRejected = countRejections(list, npix, c, order, tolerance);
      if (rejected + moreRejected > maxReject) {
        cpl_msg_debug(task, "Rejected fiber: %d (%d bad values)", 
                   i + 1, rejected + moreRejected);
        free(c);
        c = NULL;
        continue;
      }
      else if (moreRejected) {    /*  Iteration  */
        free(c);
        c = NULL;
        c = fit1DPoly(order, list, npix - moreRejected, NULL);
      }
      if (c) {
        for (k = 0; k <= order; k++) {
          snprintf(colName, MAX_COLNAME, "c%d", k);
          cpl_table_set_double(coeff, colName, i, c[k]);
        }
        free(c);
        c = NULL;
      }
    }

  }

  deleteDpoint(list);

  return coeff;

}


/**
 * @memo
 *   Match fibers traced on science with fibers traced on flat
 *
 * @return Table with matches
 *
 * @param short_fcoeff Table with fit to the short tracings on flats
 * @param short_coeff  Table with fit to the short tracings on science
 * @param row          Reference row
 * @param dc0          Output median offset between matching fibers
 * @param dc1          Output median slope difference between matching fibers
 *
 * @doc
 *   The X values at the reference row are computed for each set of
 *   coefficients, both for science and for flat tracings. The two
 *   sets of values are compared: for each science fiber X position,
 *   the closest flat fiber X position is searched. If this is at
 *   a distance greater than 2.5 pixels it is rejected. The set of
 *   accepted matches is further used, for computing the median
 *   difference between slopes. The returned table contains the
 *   tracing sequence number, the corresponding identified fiber 
 *   number, the position of the flat fibers at the reference row,
 *   and its offset position with the matching science tracing.
 *
 * @author C. Izzo
 */

cpl_table *ifuMatch(cpl_table *short_fcoeff, cpl_table *short_coeff, 
                    int row, double *dc0, double *dc1)
{
  
  cpl_table *matches;
  char       colName[MAX_COLNAME];
  double     fpos[N_BLOCKS * FIBERS_PER_BLOCK];
  double     fslope[N_BLOCKS * FIBERS_PER_BLOCK];
  double     position;
  double     slope;
  int        countTraced = cpl_table_get_nrow(short_coeff);
  int        countFibers = cpl_table_get_nrow(short_fcoeff);
  int        order       = cpl_table_get_ncol(short_fcoeff) - 2;
  int        fiber;
  int        null        = 0;
  int        i, j;
  double    *c;


  if (countFibers != N_BLOCKS * FIBERS_PER_BLOCK)
    return NULL;

  matches = cpl_table_new(countTraced);
  cpl_table_new_column(matches, "science", CPL_TYPE_INT);
  cpl_table_new_column(matches, "flat", CPL_TYPE_INT);
  cpl_table_new_column(matches, "position", CPL_TYPE_DOUBLE);
  cpl_table_new_column(matches, "offset", CPL_TYPE_DOUBLE);
  cpl_table_new_column(matches, "dslope", CPL_TYPE_DOUBLE);


  /*
   *  Compute all X positions of flat tracings at the reference row.
   */

  c = cpl_malloc((order + 1) * sizeof(double));

  for (i = 0; i < countFibers; i++) {
    for (j = 0; j <= order; j++) {
      snprintf(colName, MAX_COLNAME, "c%d", j);
      c[j] = cpl_table_get_double(short_fcoeff, colName, i, &null);
      if (null)
        break;
    }
    if (null) {
      null = 0;
      fpos[i] = -1.0;
      continue;
    }

    fpos[i] = modelValue1D(c, order, row);
    fslope[i] = c[1];
  }


  /*
   *  Loop on science tracings, finding their X position at the
   *  reference row. Then find the closest flat tracing, and accept
   *  it only if its offset is less than 2.5. Write matching fiber
   *  number to the output table, slope difference, and position
   *  of the science trace.
   */

  for (i = 0; i < countTraced; i++) {
    cpl_table_set_int(matches, "science", i, i + 1);
    for (j = 0; j <= order; j++) {
      snprintf(colName, MAX_COLNAME, "c%d", j);
      c[j] = cpl_table_get_double(short_coeff, colName, i, &null);
      if (null)
        break;
    }
    if (null) {
      null = 0;
      continue;
    }

    position = modelValue1D(c, order, row);
    slope = c[1];

    for (j = 0; j < countFibers; j++) {
      if (fpos[j] > 0.0) {
        if (fabs(fpos[j] - position) < 2.5) {
          cpl_table_set_int(matches, "flat", i, j + 1);
          cpl_table_set_double(matches, "position", i, position);
          cpl_table_set_double(matches, "dslope", i, slope - fslope[j]);
/*          cpl_table_set_double(matches, "offset", i, position);  */
          break;
        }
      }
    }
  }


  /*
   *  The median slope is simple.
   */

  *dc1 = cpl_table_get_column_median(matches, "dslope");

  /*
   *  The median X-shift between traces is a bit more complicated,
   *  because the rotation must be applied before to measure the offset.
   *  Recompute all X positions of flat tracings at the reference row.
   */

  for (i = 0; i < countFibers; i++) {
    for (j = 0; j <= order; j++) {
      snprintf(colName, MAX_COLNAME, "c%d", j);
      c[j] = cpl_table_get_double(short_fcoeff, colName, i, &null);
      if (null)
        break;
    }
    if (null) {
      null = 0;
      fpos[i] = -1.0;
      continue;
    }

    c[1] += *dc1;

    fpos[i] = modelValue1D(c, order, row);
  }

  cpl_free(c);


  /*
   *  Loop again on science tracings, finding their offset from the
   *  new (rotated) flat field tracings.
   */

  for (i = 0; i < countTraced; i++) {
    position = cpl_table_get_double(matches, "position", i, &null);
    if (null)
      continue;

    fiber = cpl_table_get_int(matches, "flat", i, &null);
    if (null)
      continue;

    position -= fpos[fiber - 1];
    cpl_table_set_double(matches, "offset", i, position);
  }

  /*
   *  Median offset.
   */

  *dc0 = cpl_table_get_column_median(matches, "offset");


  return matches;

}


/**
 * @memo
 *   Interpolate missing tracings.
 *
 * @return 0 on success.
 *
 * @param coeff Table with fits to the long tracings on flats
 *
 * @doc
 *   The input table is the main output of the function ifuFit(), the
 *   table containing the coefficients of all the fits to the flat field
 *   tracings. Some of these fits may be missing: their coefficients
 *   are added as an interpolation of the available nearby tracings.
 *   However, no interpolated solution is computed if the interpolation
 *   region includes one or more gaps between different fibers blocks.
 *
 * @author C. Izzo
 */

int ifuFillTracings(cpl_table *coeff, cpl_table *model)
{

  double *c1;
  double *c2;
  double *c3;
  int     countFibers = cpl_table_get_nrow(coeff);
  int     order       = cpl_table_get_ncol(coeff) - 2;
  int     null        = 0;
  int     skip        = 0;
  int     lastGood    = -1;
  int     firstGood   = 0;
  int     inNull      = 1;
  int     limit[]     = {0, 79, 80, 159, 160, 239, 240, 319, 320, 399};
  int     limitCount  = 10;
  int     i, j, k;
  char    colName[MAX_COLNAME];


  if (countFibers != N_BLOCKS * FIBERS_PER_BLOCK)
    return 1;

  c1 = cpl_malloc((order + 1) * sizeof(double));
  c2 = cpl_malloc((order + 1) * sizeof(double));
  c3 = cpl_malloc((order + 1) * sizeof(double));

  lastGood = -1;
  for (i = 0; i < countFibers; i++) {
    null = !cpl_table_is_valid(coeff, "c0", i);
    if (inNull) {
      if (!null) {
        inNull = 0;
        firstGood = i;

        /*
         *  Avoid interpolations across gaps
         */

        if (lastGood < 0)
          continue;

        for (j = 0, skip = 0; j < limitCount; j++) {
          if (limit[j] > lastGood && limit[j] < firstGood) {
            skip = 1;
            break;
          }
        }
        if (skip)
          continue;

        /*
         *  Interpolation of missing solutions
         */

        for (j = 0; j <= order; j++) {
          snprintf(colName, MAX_COLNAME, "c%d", j);
          c1[j] = cpl_table_get_double(coeff, colName, lastGood, NULL);
          c2[j] = cpl_table_get_double(coeff, colName, firstGood, NULL);
        }

        for (k = lastGood + 1; k < firstGood; k++) {
          for (j = 0; j <= order; j++) {
            snprintf(colName, MAX_COLNAME, "c%d", j);
            c3[j] = (k - lastGood) * c2[j] + (firstGood - k) * c1[j];
            c3[j] /= firstGood - lastGood;
            cpl_table_set_double(coeff, colName, k, c3[j]);
          }
          snprintf(colName, MAX_COLNAME, "f%d", k + 1);
          drawModel(model, colName, c3, order);
        }

      }
    }
    else {
      if (null) {
        inNull = 1;
        lastGood = i - 1;
      }
    }

  }

  cpl_free(c1);
  cpl_free(c2);
  cpl_free(c3);

  return 0;

}


/**
 * @memo
 *   Create table with tracing models.
 *
 * @return Table with tracing models.
 *
 * @param coeff Table with fit to the long tracings on flats
 *
 * @doc
 *   The input table is the first output of the function ifuFit().
 *   The coefficients are used to generate the new output model table
 *   with the computed tracings.
 *
 * @author C. Izzo
 */

cpl_table *ifuComputeTraces(cpl_table *coeff, int row, int above, int below)
{

  cpl_table *model;
  double    *c;
  int       *idata;
  int        countFibers = cpl_table_get_nrow(coeff);
  int        order       = cpl_table_get_ncol(coeff) - 2;
  int        null        = 0;
  int        range, offset;
  int        i, j;
  char       colName[MAX_COLNAME];


  if (countFibers != N_BLOCKS * FIBERS_PER_BLOCK)
    return NULL;

  range = above + below + 1;
  offset = row - below;

  model = cpl_table_new(range);

  cpl_table_new_column(model, "y", CPL_TYPE_INT);
  cpl_table_fill_column_window_int(model, "y", 0, range, 1);
  idata = cpl_table_get_data_int(model, "y");

  for (i = 0; i < range; i++)
    idata[i] = i;

  cpl_table_add_scalar(model, "y", offset);

  c = cpl_malloc((order + 1) * sizeof(double));

  for (i = 0; i < countFibers; i++) {
    for (j = 0; j <= order; j++) {
      snprintf(colName, MAX_COLNAME, "c%d", j);
      c[j] = cpl_table_get_double(coeff, colName, i, &null);
      if (null)
        break;
    }

    snprintf(colName, MAX_COLNAME, "f%d", i + 1);
    cpl_table_new_column(model, colName, CPL_TYPE_FLOAT);

    if (null) {
      null = 0;
      continue;
    }

    drawModel(model, colName, c, order);
  }

  cpl_free(c);

  return model;

}


/**
 * @memo
 *   Create table with shifted and rotated tracing models.
 *   
 * @return Table with shifted and rotated tracing models.
 * 
 * @param coeff Table with fit to the long tracings on flats
 * @param model Table with model of the long tracings on flats
 * @param dc0   Output median offset between matching fibers
 * @param dc1   Output median slope difference between matching fibers
 * 
 * @doc
 *   The two input tables are the outputs of the function ifuFit(). The 
 *   two input number are the offsets computed by the function ifuMatch(). 
 *   The @em model table is used just as a template for the output table, 
 *   and is left untouched. The alignment is accomplished by adding the 
 *   offset @em dc0 to the column "c0" of the @em coeff table, and the 
 *   median slope difference @em dc1 to the column "c1". After this, the 
 *   new tracings coefficients are used to generate the new output model 
 *   table with the flat tracings recomputed for the science exposure. 
 *
 * @author C. Izzo
 */

cpl_table *ifuAlign(cpl_table *coeff, cpl_table *model, double dc0, double dc1)
{

  cpl_table *amodel;
  double    *c;
  int        countFibers = cpl_table_get_nrow(coeff);
  int        order       = cpl_table_get_ncol(coeff) - 2;
  int        null        = 0;
  int        i, j;
  char       colName[MAX_COLNAME];


  if (countFibers != N_BLOCKS * FIBERS_PER_BLOCK)
    return NULL;

  amodel = cpl_table_duplicate(model);

  cpl_table_add_scalar(coeff, "c0", dc0);
  cpl_table_add_scalar(coeff, "c1", dc1);

  c = cpl_malloc((order + 1) * sizeof(double));

  for (i = 0; i < countFibers; i++) {
    for (j = 0; j <= order; j++) {
      snprintf(colName, MAX_COLNAME, "c%d", j);
      c[j] = cpl_table_get_double(coeff, colName, i, &null);
      if (null)
        break;
    }
    if (null) {
      null = 0;
      continue;
    }

    // c[0] += dc0;
    // c[1] += dc1;

    snprintf(colName, MAX_COLNAME, "f%d", i + 1);
    drawModel(amodel, colName, c, order);
  }

  cpl_free(c);

  return amodel;

}


/**
 * @memo
 *   Extraction of IFU spectra.
 *
 * @return Table with extracted spectra.
 *
 * @param image Bias subtracted IFU raw image.
 * @param model Table with model of tracings for all fibers.
 *
 * @doc
 *   For each fiber having a tracing model, a spectrum is extracted
 *   and written to the output table as a function of the CCD pixel
 *   along the dispersion direction (Y). For a given fiber, the 
 *   extracted value for each CCD pixel is obtained from the values 
 *   of the 3 pixels along the cross-dispersion direction (X) that 
 *   are closest to the position of the centroid obtained from the 
 *   tracing model. These values are divided by the normalized fiber
 *   profile values obtained with the function dx_to_value(). The model
 *   of the spatial profile is tuned (by horizontal stretching/shrinking)
 *   to an estimate of the true spatial profile of each fiber spectrum.
 *   The three values will require typically no cross-talk correction
 *   (with the only possible exception of extremely bright spectra
 *   close to extremely dim ones, that never happens), and would 
 *   provide 3 independent estimates of the total flux. In the 
 *   approximation of this function, a simple average of these three 
 *   values is taken as the extracted value.
 *
 * @note
 *   Typically, if the input image is an IFU science exposure, then
 *   the input table is the product of the function ifuAlign(); if
 *   the input image is an IFU flat field exposure, then the input
 *   table is the second product of the function ifuFit().
 *
 * @author C. Izzo
 */

cpl_table *ifuExtraction(cpl_image *image, cpl_table *model)
{

  char       task[]      = "ifuExtraction";

  int        nx          = cpl_image_get_size_x(image);
  float     *data        = cpl_image_get_data(image);
  float     *line;

  int        countFibers = cpl_table_get_ncol(model) - 1;
  int        range       = cpl_table_get_nrow(model);
  double    *spectrum;
  double    *fdata;
  float     *cdata;

  double  c;             /* Running position of the centroid                */
  int     x1, x2, x3;    /* Positions of the pixels closest to the centroid */
  double  d1, d2, d3;    /* Distances of the above pixels from the centroid */
  double  v1, v2, v3;    /* Values of the above pixels                      */
  double  c1, c2, c3;    /* Normalized values of the above pixels           */
  double  value;         /* Estimated spectral flux                         */

  double  norm = 3.017532; /* This is the normalization for transforming
                            * the normalized pixel value to an estimate of
                            * the total spectral flux, obtained from the
                            * flux_constant() function.
                            */

  int        i, j, k;
  int        offset      = cpl_table_get_int(model, "y", 0, NULL);
  double     variance, factor, step, f, min, sum_weights;
  char       colTrace[MAX_COLNAME];
  char       colSpectrum[MAX_COLNAME];


  if (countFibers != N_BLOCKS * FIBERS_PER_BLOCK)
    return NULL;

  cpl_table *spectra     = cpl_table_new(range);
  cpl_table_duplicate_column(spectra, "y", model, "y");

  cpl_table *factors = cpl_table_new(range);
  cpl_table_duplicate_column(factors, "y", model, "y");

  for (i = 0; i < countFibers; i++) {                   /* Loop on fibers */

    snprintf(colTrace, MAX_COLNAME, "f%d", i + 1);
    snprintf(colSpectrum, MAX_COLNAME, "s%d", i + 1);

    if (cpl_table_has_invalid(model, colTrace)) {
      cpl_msg_debug(task, "Trace not available for spectrum %d\n", i + 1);
      continue;
    }
    cdata = cpl_table_get_data_float(model, colTrace);

    cpl_table_new_column(factors, colSpectrum, CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(factors, colSpectrum, 0, range, 0.0);
    fdata = cpl_table_get_data_double(factors, colSpectrum);

    /*
     * First iteration, looking for the best "stretching" of the
     * standard spatial profile adapting to the current fibre.
     */

    for (j = 0; j < range; j++) {    /* Loop along the dispersion */

      line = data + (j + offset) * nx;
/* Too slow:
      c = cpl_table_get_float(model, colTrace, j, NULL);
*/
      c = cdata[j];                  /* Centroid                */
      x1 = c - 0.5;                  /* Intentional truncation  */
      x2 = x1 + 1;
      x3 = x2 + 1;
      if (x1 > 0 && x3 < nx) {
        d1 = c - x1;                 /* Distances from centroid */ 
        d2 = c - x2;
        d3 = c - x3;
        v1 = line[x1];               /* Pixel values            */
        v2 = line[x2];
        v3 = line[x3];
        
        step = 0.05;
        for (k = 0; k < 22; k++) {
          f = 0.5 + k*step;
          c1 = v1 / dx_to_value(f*d1);        /* Normalization  */
          c2 = v2 / dx_to_value(f*d2);
          c3 = v3 / dx_to_value(f*d3);
          value = (c1 + c2 + c3) / 3;         /* Mean value     */
          variance = (value - c1)*(value - c1)
                   + (value - c2)*(value - c2)
                   + (value - c3)*(value - c3);
          if (k) {
            if (min > variance) {
              min = variance;
              factor = f;
            }
          }
          else {
            min = variance;
            factor = f;
          }
        }
        if (factor < 0.55 || factor > 1.5)
            cpl_table_set_invalid(factors, colSpectrum, j);
        else
            fdata[j] = factor;
/* Too slow:
            cpl_table_set_double(factors, colSpectrum, j, factor);
*/
      }

    } /* End loop along the dispersion */

    cpl_table_new_column(spectra, colSpectrum, CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(spectra, colSpectrum, 0, range, 0.0);
    spectrum = cpl_table_get_data_double(spectra, colSpectrum);
    f = cpl_table_get_column_median(factors, colSpectrum);
/* *+
    printf("median of fiber %d = %f\n", i, f);
+* */

    for (j = 0; j < range; j++) {            /* Loop along the dispersion */

      line = data + (j + offset) * nx;
/* Too slow:
      c = cpl_table_get_float(model, colTrace, j, NULL);
*/
      c = cdata[j];                  /* Centroid                */
      x1 = c - 0.5;                  /* Intentional truncation */
      x2 = x1 + 1;
      x3 = x2 + 1;
      if (x1 > 0 && x3 < nx) {
        d1 = c - x1;                 /* Distances from centroid */
        d2 = c - x2;
        d3 = c - x3;
        v1 = line[x1];               /* Pixel values            */
        v2 = line[x2];
        v3 = line[x3];
/*
 * Old, non-optimal code:

        c1 = v1 / dx_to_value(f*d1);   \* Normalization           *\
        c2 = v2 / dx_to_value(f*d2);
        c3 = v3 / dx_to_value(f*d3);
        value = (c1 + c2 + c3) / 3;  \* Mean value              *\
 End of old, non-optimal code */

 /* New, optimal code: */

        sum_weights = dx_to_value(f*d1) * dx_to_value(f*d1)
                    + dx_to_value(f*d2) * dx_to_value(f*d2)
                    + dx_to_value(f*d3) * dx_to_value(f*d3);

        /*
         * Here come the optimal estimate (assuming that the variance is
         * background dominated, and therefore constant, that is true for
         * weak objects, while for bright objects S/N is good anyway).
         * Note that the weights are dx_to_value(f*di) squared, and
         * IT IS CORRECT NOT TO SQUARE the profile in the sum of the
         * measured profile values - do the algebra if you want to see 
         * why...
         */

        c1 = v1 * dx_to_value(f*d1);
        c2 = v2 * dx_to_value(f*d2);
        c3 = v3 * dx_to_value(f*d3);
        value = (c1 + c2 + c3) / sum_weights;  /* Flux estimator */

        value *= norm*f;               /* ADU for each Y CCD pixel */
        spectrum[j] = value;
/* Too slow:
        cpl_table_set_double(spectra, colSpectrum, j, value);
*/
      }

    } /* End loop along the dispersion */

  } /* End loop on fibers */

/*
  cpl_table_save(factors, NULL, NULL, "factors.fits", CPL_IO_CREATE);
*/
  cpl_table_delete(factors);

  return spectra;

}


/**
 * @memo
 *   Simple extraction of IFU spectra.
 *
 * @return Table with extracted spectra.
 *
 * @param image Bias subtracted IFU raw image.
 * @param model Table with model of tracings for all fibers.
 *
 * @doc
 *   For each fiber having a tracing model, a spectrum is extracted
 *   and written to the output table as a function of the CCD pixel
 *   along the dispersion direction (Y). For a given fiber, the 
 *   extracted value for each CCD pixel is obtained from the values 
 *   of the 3 pixels along the cross-dispersion direction (X) that 
 *   are closest to the position of the centroid obtained from the 
 *   tracing model. These values are divided by the normalized fiber
 *   profile values obtained with the function dx_to_value(). These
 *   three values will require typically no cross-talk correction
 *   (with the only possible exception of extremely bright spectra
 *   close to extremely dim ones), and would provide 3 independent
 *   estimates of the total flux. In the simple approximation of
 *   this function, a simple average of these three values is taken 
 *   as the extracted value.
 *
 * @note
 *   Typically, if the input image is an IFU science exposure, then
 *   the input table is the product of the function ifuAlign(); if
 *   the input image is an IFU flat field exposure, then the input
 *   table is the second product of the function ifuFit().
 *
 * @author C. Izzo
 */

cpl_table *ifuSimpleExtraction(cpl_image *image, cpl_table *model)
{

  char       task[]      = "ifuSimpleExtraction";

  int        nx          = cpl_image_get_size_x(image);
  float     *data        = cpl_image_get_data(image);
  float     *line;

  int        countFibers = cpl_table_get_ncol(model) - 1;
  int        range       = cpl_table_get_nrow(model);
  cpl_table *spectra     = cpl_table_new(range);

  double  c;             /* Running position of the centroid                */
  int     x1, x2, x3;    /* Positions of the pixels closest to the centroid */
  double  d1, d2, d3;    /* Distances of the above pixels from the centroid */
  double  v1, v2, v3;    /* Values of the above pixels                      */
  double  c1, c2, c3;    /* Normalized values of the above pixels           */
  double  value;         /* Estimated spectral flux                         */

  double  norm = 3.017532; /* This is the normalization for transforming
                            * the normalized pixel value to an estimate of
                            * the total spectral flux, obtained from the
                            * flux_constant() function.
                            */

  int        i, j;
  int        offset      = cpl_table_get_int(model, "y", 0, NULL);
  char       colTrace[MAX_COLNAME];
  char       colSpectrum[MAX_COLNAME];


  if (countFibers != N_BLOCKS * FIBERS_PER_BLOCK)
    return NULL;

  cpl_table_duplicate_column(spectra, "y", model, "y");

  for (i = 0; i < countFibers; i++) {                   /* Loop on fibers */

    snprintf(colTrace, MAX_COLNAME, "f%d", i + 1);
    snprintf(colSpectrum, MAX_COLNAME, "s%d", i + 1);
    if (cpl_table_has_invalid(model, colTrace)) {
      cpl_msg_debug(task, "Trace not available for spectrum %d\n", i + 1);
      continue;
    }
    cpl_table_new_column(spectra, colSpectrum, CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(spectra, colSpectrum, 0, range, 0.0);

    for (j = 0; j < range; j++) {            /* Loop along the dispersion */

      line = data + (j + offset) * nx;
      c = cpl_table_get_float(model, colTrace, j, NULL);
      x1 = c - 0.5;                  /* Intentional truncation */
      x2 = x1 + 1;
      x3 = x2 + 1;
      if (x1 > 0 && x3 < nx) {
        d1 = c - x1;                 /* Distances from centroid */
        d2 = c - x2;
        d3 = c - x3;
        v1 = line[x1];               /* Pixel values            */
        v2 = line[x2];
        v3 = line[x3];
        c1 = v1 / dx_to_value(d1);   /* Normalization           */
        c2 = v2 / dx_to_value(d2);
        c3 = v3 / dx_to_value(d3);
        value = (c1 + c2 + c3) / 3;  /* Mean value              */
        value *= norm;               /* ADU for each Y CCD pixel */
        cpl_table_set_double(spectra, colSpectrum, j, value);
      }

    } /* End loop along the dispersion */

  } /* End loop on fibers */

  return spectra;

}


/**
 * @memo
 *   Very simple extraction of IFU spectra.
 *
 * @return Table with extracted spectra.
 *
 * @param image Bias subtracted IFU raw image.
 * @param model Table with model of tracings for all fibers.
 *
 * @doc
 *   For each fiber having a tracing model, a spectrum is extracted
 *   and written to the output table as a function of the CCD pixel
 *   along the dispersion direction (Y). For a given fiber, the 
 *   extracted value for each CCD pixel is obtained from the value 
 *   of the pixel along the cross-dispersion direction (X) that 
 *   is closest to the position of the centroid obtained from the 
 *   tracing model. This value is divided by the normalized fiber
 *   profile values obtained with the function dx_to_value(). This
 *   value will require typically no cross-talk correction (with 
 *   the only possible exception of extremely bright spectra close
 *   to extremely dim ones).
 *
 * @note
 *   Typically, if the input image is an IFU science exposure, then
 *   the input table is the product of the function ifuAlign(); if
 *   the input image is an IFU flat field exposure, then the input
 *   table is the second product of the function ifuFit().
 *
 * @author C. Izzo
 */

cpl_table *ifuVerySimpleExtraction(cpl_image *image, cpl_table *model)
{

  char       task[]      = "ifuVerySimpleExtraction";

  int        nx          = cpl_image_get_size_x(image);
  float     *data        = cpl_image_get_data(image);
  float     *line;

  int        countFibers = cpl_table_get_ncol(model) - 1;
  int        range       = cpl_table_get_nrow(model);
  cpl_table *spectra     = cpl_table_new(range);
  double    *spectrum;

  double     c;          /* Running position of the centroid              */
  int        x;          /* Position of the pixel closest to the centroid */
  double     d;          /* Distance of the above pixel from the centroid */
  double     v;          /* Value of the above pixel                      */
  double     value;      /* Estimated spectral flux                       */

  int        i, j;
  int        offset      = cpl_table_get_int(model, "y", 0, NULL);
  char       colTrace[MAX_COLNAME];
  char       colSpectrum[MAX_COLNAME];


  if (countFibers != N_BLOCKS * FIBERS_PER_BLOCK)
    return NULL;

  cpl_table_duplicate_column(spectra, "y", model, "y");

  for (i = 0; i < countFibers; i++) {                   /* Loop on fibers */

    snprintf(colTrace, MAX_COLNAME, "f%d", i + 1);
    snprintf(colSpectrum, MAX_COLNAME, "s%d", i + 1);
    if (cpl_table_has_invalid(model, colTrace)) {
      cpl_msg_debug(task, "Trace not available for spectrum %d\n", i + 1);
      continue;
    }
    cpl_table_new_column(spectra, colSpectrum, CPL_TYPE_DOUBLE);
    cpl_table_fill_column_window_double(spectra, colSpectrum, 0, range, 0.0);
    spectrum = cpl_table_get_data_double(spectra, colSpectrum);

    for (j = 0; j < range; j++) {            /* Loop along the dispersion */

      line = data + (j + offset) * nx;
      c = cpl_table_get_float(model, colTrace, j, NULL);
      x = c + 0.5;                           /* Intentional truncation    */
      if (x > 0 && x < nx) {
        d = c - x;                           /* Distances from centroid   */
        v = line[x];                         /* Pixel values              */
        value = v / dx_to_value(d);          /* Normalization             */
        cpl_table_set_double(spectra, colSpectrum, j, value);
      }

    } /* End loop along the dispersion */

  } /* End loop on fibers */

  return spectra;

}


/**
 * @memo
 *   Generate transmission correction from flat field spectra.
 *
 * @return 0 on success.
 *
 * @param image     Image containing wavelength calibrated flat field spectra.
 * @param startPix  Start integration pixel.
 * @param endPix    End integration pixel.
 * @param norm      Returned normalization factor (median value).
 *
 * @doc
 *   The input image contains the extracted and wavelength calibrated
 *   flat field spectra for each fiber. If the input image contains
 *   1600 spectra (LR grisms case), the first (bottom) 400 spectra
 *   refer to pseudo-slit 0, the next 400 to pseudo-slit 1, ..., and
 *   the last (top) 400 to pseudo-slit 3. If the input image contains
 *   400 spectra (MR and HR grisms case), the spectra refer to pseudo-
 *   slit 1 (the only one used). In all cases, the spectra of a pseudo-
 *   slit starting from the left side of the chip (lowest X CCD ccordinates)
 *   are stored from bottom to top in the output image.
 *
 *   The output table contains the fluxes integrated from each spectrum
 *   on the specified interval, and normalized to their median value.
 *
 * @author C. Izzo
 */

cpl_table *ifuTransmission(cpl_image *image, 
                           int startPix, int endPix, double *norm, double *err)
{

/*  char task[]      = "ifuTransmission";   */

  int        nx          = cpl_image_get_size_x(image);
  int        ny          = cpl_image_get_size_y(image);
  float     *data        = cpl_image_get_data(image);
  float     *line        = data;
  cpl_table *table       = cpl_table_new(ny);
  int        i, j;
  double     sum;
  double     level;


  cpl_table_new_column(table, "trans", CPL_TYPE_DOUBLE);

  for (i = 0; i < ny; i++, line += nx) {
    sum = 0.0;
    for (j = startPix; j < endPix; j++)
      sum += line[j];
    if (sum > 0.00001)
      cpl_table_set_double(table, "trans", i, sum);
  }

  level = cpl_table_get_column_median(table, "trans");

  cpl_table_divide_scalar(table, "trans", level);

  *norm = level;

  /* FIXME:
   * The statistical error on the computed total flux is extremely
   * small, since the flux is integrated on a wide spectral range
   * from all the available spectra. In order to compute the error 
   * rigorously, an error image should be compute at extraction
   * time, and passed to this function. Temporarily, a worst
   * scenario is adopted for an estimation of the error, setting
   * a very high value for the gain (3 e-/ADU), and assuming that
   * the spectral extraction operation did not improve the S/N ratio.
   */

  *err = sqrt(3*level);

  return table;

}


/**
 * @memo
 *   Apply transmission correction to extracted and wave-calibrated spectra
 *
 * @return 0 on success.
 *
 * @param image     Image containing wavelength calibrated spectra.
 * @param table     Table containing transmission correction factors.
 *
 * @doc
 *   The input image contains extracted and wavelength calibrated
 *   spectra for each fiber. Each image row is multiplied by the
 *   corresponding transmission factor.
 *
 * @author C. Izzo
 */

int ifuApplyTransmission(cpl_image *image, cpl_table *table)
{

  int        nx          = cpl_image_get_size_x(image);
  int        ny          = cpl_image_get_size_y(image);
  float     *data        = cpl_image_get_data(image);
  float     *line        = data;
  double     factor;
  int        null;
  int        i, j;


  for (i = 0; i < ny; i++, line += nx) {
    factor = cpl_table_get_double(table, "trans", i, &null);
    if (null)
      continue;
    if (factor < 0.00001)
      continue;
    for (j = 0; j < nx; j++)
      line[j] /= factor;
  }

  return 0;

}


/**
 * @memo
 *   Sum spectral signal.
 *
 * @return Array of integrated signals from each fiber spectrum.
 *
 * @param spectra Table with extracted spectra.
 * @param zero    Expected position of contamination.
 * @param skip    Number of pixels to skip around contamination (radius).
 *
 * @doc
 *   The input table may be the product of any IFU spectral extraction 
 *   function, from ifuSimpleExtraction() to [ADD FUNCTION NAME HERE].
 *   The only requirement is that the table contains 400 columns named
 *   "s1", "s2", "s3", ... , "s400", and one column "y" for the coordinate
 *   along the dispersion direction. Other columns may exist, and would
 *   be ignored. The values of each column are averaged and written to output. 
 *   Any column containing invalid values would be assigned a zero sum. 
 *   If the expected position of the contamination is zero, or outside 
 *   the extraction range, all the spectral values are integrated.
 *
 * @author C. Izzo
 */

double *ifuIntegrateSpectra(cpl_table *spectra, int zero, int skip)
{

/*  char    task[]      = "ifuIntegrateSpectra";   */

  char    colSpectrum[MAX_COLNAME];

  double *buffer;
  double *spectrum;
  double  sum;
  int    *y           = cpl_table_get_data_int(spectra, "y");
  int     countFibers = N_BLOCKS * FIBERS_PER_BLOCK;
  int     range       = cpl_table_get_nrow(spectra);
  int     dist;
  int     count;
  int     i, j;


  buffer = cpl_malloc(countFibers * sizeof(double));

  for (i = 0; i < countFibers; i++) {
    snprintf(colSpectrum, MAX_COLNAME, "s%d", i + 1);
    if (cpl_table_has_column(spectra, colSpectrum)) {
      if (cpl_table_has_invalid(spectra, colSpectrum)) {
        buffer[i] = 0.0;
      }
      else {
        spectrum = cpl_table_get_data_double(spectra, colSpectrum);
        for (j = 0, sum = 0.0, count = 0; j < range; j++) {
          dist = abs(y[j] - zero);
          if (dist > skip) {
            sum += spectrum[j];
            count++;
          }
        }
        buffer[i] = sum / count;
      }
    }
    else
      buffer[i] = 0.0;
  }

  return buffer;
  
}


/**
 * @memo
 *   Fill IFU reconstructed field.
 *
 * @return 0 on success.
 *
 * @param image     IFU field.
 * @param integrals Array with integrated spectral signals.
 * @param quadrant  Quadrant number [1-4].
 * @param slit      IFU slit number [0-3].
 *
 * @doc
 *   The input buffer, @em integrals, may be the product of any IFU 
 *   spectral signal integration function, from ifuIntegrateSpectra() 
 *   to [ADD FUNCTION NAME HERE]. The buffer must include 400 values, 
 *   each one of them is written to the appropriate position of the
 *   allocated 80x80 @em image given in input, according to the
 *   specified quadrant and IFU slit number.
 *
 * @author C. Izzo
 */

int ifuImage(cpl_image *image, double *integrals, int quadrant, int slit)
{

  char   task[] = "ifuImage";

  float *data   = cpl_image_get_data_float(image);
  int    lpos, mpos;
  int    startl[N_SLITS], startm[N_SLITS], dm[N_SLITS], jump[N_SLITS];
  int    dl;
  int    i, j, k;


  /*
   * Coding of the IFU table:
   *
   * Each block consists of 80 fibers.
   * Each IFU slit consists of 5 blocks.
   * Each quadrant accomodates 4 IFU slits.
   * Here all is counted starting from zero.
   * IFU slit 1 is the central one (the only one used with HR grisms).
   *
   * startl[L] indicates the start l coordinate of block 0 in IFU slit L.
   * startm[L] indicates the start m coordinate of block 0 in IFU slit L.
   * dm[L]     indicates whether the m coordinate increases or decreases
   *           during wrapping of IFU slit L.
   * jump[L]   indicates the m gap between m start positions of blocks
   *           of IFU slit L.
   */

  switch (quadrant) {
  case 1:
    dl = 1;
    startl[0] = 79;
    startl[1] = 59;
    startl[2] = 59;
    startl[3] = 79;
    startm[0] = 60;
    startm[1] = 43;             /* 43; CONSORTIUM,  56; Isabelle */
    startm[2] = 63;
    startm[3] = 43;
    dm[0] =  1;
    dm[1] = -1;                 /* -1; CONSORTIUM,   1; Isabelle */
    dm[2] = -1;
    dm[3] = -1;
    jump[0] =  4;
    jump[1] =  4;               /*  4; CONSORTIUM,  -4; Isabelle */
    jump[2] =  4;
    jump[3] =  4;
    break;
  case 2:
/*  Original setting!
    dl =  1;
    startl[0] = 19;
    startl[1] = 39;
    startl[2] = 39;
    startl[3] = 19;
    startm[0] = 76;
    startm[1] = 59;
    startm[2] = 79;
    startm[3] = 56;
    dm[0] =  1;
    dm[1] = -1;
    dm[2] = -1;
    dm[3] =  1;
    jump[0] = -4;
    jump[1] = -4;
    jump[2] = -4;
    jump[3] = -4;
    break;
*/
    dl =  1;
    startl[0] = 19;
    startl[1] = 39;
    startl[2] = 39;
    startl[3] = 19;
    startm[0] = 76;
    startm[1] = 59;
    startm[2] = 79;
    startm[3] = 59;
    dm[0] =  1;
    dm[1] = -1;
    dm[2] = -1;
    dm[3] = -1;
    jump[0] = -4;
    jump[1] = -4;
    jump[2] = -4;
    jump[3] = -4;
    break;
  case 3:
    dl = 1;
    startl[0] = 19;
    startl[1] = 39;
    startl[2] = 39;
    startl[3] = 19;
    startm[0] = 3;
    startm[1] = 20;
    startm[2] = 0;
    startm[3] = 20;
    dm[0] = -1;
    dm[1] =  1;
    dm[2] =  1;
    dm[3] =  1;
    jump[0] =  4;
    jump[1] =  4;
    jump[2] =  4;
    jump[3] =  4;
    break;
  case 4:
    dl = 1;
    startl[0] = 79;
    startl[1] = 59;
    startl[2] = 59;
    startl[3] = 79;
    startm[0] = 19;
    startm[1] = 36;             /* 36; CONSORTIUM, 23; Isabelle */
    startm[2] = 16;
    startm[3] = 36;
    dm[0] = -1;
    dm[1] =  1;                 /* 1; CONSORTIUM, -1; Isabelle */
    dm[2] =  1;
    dm[3] =  1;
    jump[0] = -4;
    jump[1] = -4;                /* -4; CONSORTIUM, 4; Isabelle */
    jump[2] = -4;
    jump[3] = -4;
    break;
  default:
    cpl_msg_error(task, "Wrong quadrant number (you should never get here!)");
    return 1;
  }

  /*
   *  Write values to image
   */

  lpos = startl[slit];

  k = 0;
  for (i = 0; i < N_BLOCKS; i++) {
    if (quadrant == 2 && slit == 3) {
      if (i == 4) {
        startm[slit] = 47;
        jump[slit] = 0;
        dm[slit] = -1;
      }
      if (i == 3) {
        startm[slit] = 43;
        jump[slit] = 0;
        dm[slit] = -1;
      }
    }
    mpos = startm[slit] + i * jump[slit];
    for (j = 0; j < FIBERS_PER_BLOCK / 4; j++) {
      data[lpos + mpos * FIBERS_PER_BLOCK] = integrals[k];
      k++;
/* printf("%d %d %d\n", k + slit * 400, lpos + 1, mpos + 1); */
/* printf("%d %d %d\n", k, lpos + 1, mpos + 1); */
      lpos -= dl;
    }
    mpos += dm[slit];
    for (j = 0; j < FIBERS_PER_BLOCK / 4; j++) {
      lpos += dl;
      data[lpos + mpos * FIBERS_PER_BLOCK] = integrals[k];
      k++;
/* printf("%d %d %d\n", k + slit * 400, lpos + 1, mpos + 1); */
/* printf("%d %d %d\n", k, lpos + 1, mpos + 1); */
    }
    mpos += dm[slit];
    for (j = 0; j < FIBERS_PER_BLOCK / 4; j++) {
      data[lpos + mpos * FIBERS_PER_BLOCK] = integrals[k];
      k++;
/* printf("%d %d %d\n", k + slit * 400, lpos + 1, mpos + 1); */
/* printf("%d %d %d\n", k, lpos + 1, mpos + 1); */
      lpos -= dl;
    }
    mpos += dm[slit];
    for (j = 0; j < FIBERS_PER_BLOCK / 4; j++) {
      lpos += dl;
      data[lpos + mpos * FIBERS_PER_BLOCK] = integrals[k];
      k++;
/* printf("%d %d %d\n", k + slit * 400, lpos + 1, mpos + 1); */
/* printf("%d %d %d\n", k, lpos + 1, mpos + 1); */
    }
  }

  return 0;

}


/**
 * @memo
 *   Compute wavelength calibration for each extracted arc lamp spectrum.
 *
 * @return Table with wavelength calibrations.
 *
 * @param spectra  Extracted arc lamp spectra for each fiber.
 * @param linecat  Line catalog.
 * @param coeff    First guess wavelength calibration.
 * @param order    Order of the IDS polynomial.
 * @param lambda   Reference wavelength.
 * @param zero     Expected position of the zero order contamination.
 *
 * @doc
 *   This function returns the inverse dispersion solution for each
 *   extracted arc lamp spectrum. The first guess wavelength calibration
 *   may come from the appropriate calibration file, or from the function 
 *   ifuFirstIds().
 *
 * @author C. Izzo
 */

cpl_table *ifuComputeIds(cpl_table *spectra, cpl_table *linecat, 
                         double *coeff, int order, double lambda, 
                         int zero, double maxRms)
{

/*  char       task[] = "ifuComputeIds";  */

  char         colName[MAX_COLNAME];
  int          countFibers = N_BLOCKS * FIBERS_PER_BLOCK;
  int          nrow        = cpl_table_get_nrow(spectra);
  int          nlines      = cpl_table_get_nrow(linecat);
  int          offset      = cpl_table_get_int(spectra, "y", 0, NULL);
  cpl_table   *coeffTable  = cpl_table_new(countFibers);
  float       *wdata       = cpl_table_get_data_float(linecat, "WLEN");
  VimosDpoint *list        = newDpoint(nlines);
  double      *data;
  double      *c = NULL;
  double       rms;
  int          contamRadius = 15;
  int          bigSearchRadius = 30;
  int          startSearchRadius = 3;
  int          endSearchRadius = 15;
  int          searchRadius;
  int          countLines;

  double       level;             /* Current detection level               */
  double       aboveLevel = 120.; /* Start detection level                 */
  double       thresRms = 1.0;    /* This is used by the iteration process */
  double       pos;
  int          ipos;
  int          maxIter = 4;
  int          yStart, yEnd, length;
  int          yStartShort, yEndShort, lengthShort;
  int          found, foundFirst;
  int          skipFirst = 1;
  int          i, j, k;


  for (i = 0; i <= order; i++) {
    snprintf(colName, MAX_COLNAME, "c%d", i);
    cpl_table_new_column(coeffTable, colName, CPL_TYPE_DOUBLE);
  }
  cpl_table_new_column(coeffTable, "rms", CPL_TYPE_DOUBLE);
  cpl_table_new_column(coeffTable, "nlines", CPL_TYPE_INT);
  cpl_table_fill_column_window_int(coeffTable, "nlines", 0, countFibers, 0);

repeat:

  /*
   *  Determine interval where to look for the reference arc line
   */

  yStart = coeff[0] - offset - bigSearchRadius;
  yEnd = coeff[0] - offset + bigSearchRadius;
  if (yStart < 0)
    yStart = 0;
  if (yEnd > nrow)
    yEnd = nrow;
  if (yEnd < 0 || yStart >= nrow)
    return NULL;

  length = yEnd - yStart;


  /*
   *  Lines identification
   */

  for (i = 0; i < countFibers; i++) {
    cpl_msg_debug(cpl_func,"Computing wavelength solution for fiber %d", i+1);
    snprintf(colName, MAX_COLNAME, "s%d", i + 1);

    if (!cpl_table_has_column(spectra, colName))
      continue;

    if (cpl_table_has_invalid(spectra, colName))
      continue;

    if (skipFirst) {
        skipFirst = 0;
        continue;
    }

    /*
     * Determine significant level for line detection at 5 * sigma.
     * Determine standard deviation from the median level using only
     * the negative deviations.
     */

    data = cpl_table_get_data_double(spectra, colName);
    level = cpl_table_get_column_median(spectra, colName);
    cpl_msg_debug(cpl_func,"The median level of fiber %d is %f", i+1, level);

    ipos = whereMax(data + yStart, length);
    coeff[0] = yStart + offset + ipos;

    for (k = 0; k < maxIter; k++) {
      level += aboveLevel;
      countLines = 0;
      for (j = 0; j < nlines; j++) {
        ipos = modelValue1D(coeff, order, wdata[j] - lambda);

        /*
         * Skip line if too close to zero order contamination
         */

        if (abs(ipos - zero) < contamRadius + endSearchRadius)
          continue;

        found = 0;
        for (searchRadius = startSearchRadius; searchRadius <= endSearchRadius;
             searchRadius++) {

          /*
           *  Determine interval where to look for the current arc line
           */

          yStartShort = ipos - offset - searchRadius;
          yEndShort = ipos - offset + searchRadius;
          if (yStartShort < 0)
            break;                              /* yStartShort = 0; */
          if (yEndShort > nrow)
            break;                              /* yEndShort = nrow; */
          if (yEndShort < 0 || yStartShort >= nrow)
            break;
          lengthShort = yEndShort - yStartShort;
    
          if (findPeak(data + yStartShort, lengthShort, level, &pos)) {
            found = 1;
            break;
          }
        }

        if (found) {
          list[countLines].x = wdata[j] - lambda;
          list[countLines].y = yStartShort + offset + pos;
          countLines++;
        }

      }

      /*
       *  The number of identified lines should be at least twice
       *  the degrees of freedom.
       */

      if (countLines < 2 * (order + 1))
        continue;

      c = fit1DPoly(order, list, countLines, &rms);
      if (c) {
        if (rms < thresRms)
          break;
        if (k == maxIter - 1)
          break;         /* If last, don't free(c) (keep the bad solution) */
        free(c);
        c = NULL;
      }
    }

    if (c) {
      if (rms < maxRms) {
        for (j = 0; j <= order; j++) {
          snprintf(colName, MAX_COLNAME, "c%d", j);
          cpl_table_set_double(coeffTable, colName, i, c[j]);
        }
        cpl_table_set_double(coeffTable, "rms", i, sqrt(rms));
        cpl_table_set_int(coeffTable, "nlines", i, countLines);
        for (j = 0; j <= order; j++)
          coeff[j] = c[j];
        if (!foundFirst) {
          foundFirst = 1;
          free(c);
          c = NULL;
          goto repeat;
        }
      }
      free(c);
      c = NULL;
    }
    else
    {
      cpl_msg_debug(cpl_func,
              "Wavelenght calibration failed for fiber %d", i+1);
      continue;
    }

  }

  deleteDpoint(list);

  return coeffTable;

}


/**
 * @memo
 *   Compute wavelength calibration for each extracted arc lamp spectrum.
 *
 * @return Table with wavelength calibrations.
 *
 * @param spectra    Extracted arc lamp spectra for each fiber.
 * @param linecat    Line catalog.
 * @param lambda2pix Rough expected inverse-dispersion (pixel/A).
 * @param order      Order of the IDS polynomial.
 * @param lambda     Reference wavelength.
 * @param zero       Expected position of the zero order contamination.
 *
 * @doc
 *   This function returns the inverse dispersion solution for each
 *   extracted arc lamp spectrum. There is no need of first-guess
 *   model - just a rough inverse-dispersion should be provided.
 *
 * @author C. Izzo
 */

double *ifuComputeIdsBlind(cpl_table *spectra, cpl_table *linecat, 
                           double lambda2pix, int order, double lambda, 
                           double maxRms)
{

/*  char       task[] = "ifuComputeIdsBlind";  */

  char         colName[MAX_COLNAME];
  int          countFibers = N_BLOCKS * FIBERS_PER_BLOCK;
  int          nrow        = cpl_table_get_nrow(spectra);
  int          nlines      = cpl_table_get_nrow(linecat);
  int          offset      = cpl_table_get_int(spectra, "y", 0, NULL);
  cpl_table   *coeffTable  = cpl_table_new(countFibers);
  float       *wdata       = cpl_table_get_data_float(linecat, "WLEN");
  double      *lines;
  VimosDpoint *list        = newDpoint(nlines);
  double      *peaks;
  double     **output;
  double      *data;
  double      *c = NULL;
  double      *coeff = NULL;
  double       rms;
  int          countLines;
  int          npeaks;

  double       level;             /* Spectral background level             */
  double       aboveLevel = 120.; /* Above background level                */
  double       max_disp, min_disp;
  
  int          i, j;


  for (i = 0; i <= order; i++) {
    snprintf(colName, MAX_COLNAME, "c%d", i);
    cpl_table_new_column(coeffTable, colName, CPL_TYPE_DOUBLE);
  }
  cpl_table_new_column(coeffTable, "rms", CPL_TYPE_DOUBLE);
  cpl_table_new_column(coeffTable, "nlines", CPL_TYPE_INT);
  cpl_table_fill_column_window_int(coeffTable, "nlines", 0, countFibers, 0);

  lines = cpl_malloc(nlines * sizeof(double));

  for (i = 0; i < nlines; i++)
    lines[i] = wdata[i];

  max_disp = min_disp = 1 / lambda2pix;
  max_disp += max_disp / 5.5;
  min_disp -= min_disp / 5.5;

  /*
   *  Peak detection and lines identification
   */

  for (i = 0; i < countFibers; i++) {
    snprintf(colName, MAX_COLNAME, "s%d", i + 1);

    if (!cpl_table_has_column(spectra, colName))
      continue;

    if (cpl_table_has_invalid(spectra, colName))
      continue;

    data = cpl_table_get_data_double(spectra, colName);
    level = cpl_table_get_column_median(spectra, colName);
    level += aboveLevel;

    peaks = collectPeaks_double(data, nrow, level, 1.0, &npeaks);
    cpl_msg_debug(cpl_func,"Found %d peaks for fiber %d", 
                  npeaks, i+1);

    if (peaks) {
      output = identPeaks(peaks, npeaks, lines, nlines, min_disp, max_disp,
                          0.07, &countLines);
      if (output) {
        for (j = 0; j < countLines; j++) {
          list[j].x = output[1][j] - lambda;
          list[j].y = output[0][j] + offset;
        }
        cpl_free(output[0]);
        cpl_free(output[1]);
        cpl_free(output);
      }
      cpl_free(peaks);
    }
    else
      countLines = 0;

    /*
     *  The number of identified lines should be at least twice
     *  the degrees of freedom.
     */
    cpl_msg_debug(cpl_func,"Number of identified lines for fiber %d: %d order %d", 
                  i+1, countLines, order);

    if (countLines < 2 * (order + 1))
    {
      cpl_msg_debug(cpl_func,
              "Number of lines (%d) not enough for blind fitting (%d)",
              countLines, 2 * (order + 1));
      continue;
    }

    c = fit1DPoly(order, list, countLines, &rms);

    if (c) {
      if (rms < maxRms) {
        for (j = 0; j <= order; j++) {
          snprintf(colName, MAX_COLNAME, "c%d", j);
          cpl_table_set_double(coeffTable, colName, i, c[j]);
        }
        cpl_table_set_double(coeffTable, "rms", i, sqrt(rms));
        cpl_table_set_int(coeffTable, "nlines", i, countLines);
      }
      free(c);
      c = NULL;
    }
    else
    {
        cpl_msg_debug(cpl_func, 
                "Fitting of wavelength polynomial failed for fiber %d",i+1);
    }

  }

  cpl_free(lines);
  deleteDpoint(list);

/*
 *  return coeffTable;
 */

  coeff = cpl_malloc((order + 1) * sizeof(double));
  cpl_msg_debug(cpl_func,"Mean blind wavelength ids:");
  for (i = 0; i <= order; i++) {
    snprintf(colName, MAX_COLNAME, "c%d", i);
    coeff[i] = cpl_table_get_column_median(coeffTable, colName);
    cpl_msg_debug(cpl_func,"  c[%d]: %f",i, coeff[i]);
  }

  cpl_table_delete(coeffTable);

  return coeff;

}


/**
 * @memo
 *   Resample all extracted spectra at a constant wavelength step.
 *   
 * @return Image with resampled spectra.
 * 
 * @param image       Allocated image to be filled with the resampled spectra.
 * @param spectra     Extracted arc lamp spectra for each fiber.
 * @param ids         IDS table.
 * @param slit        IFU pseudo-slit (from 0 to 3, with 0 the isolated slit).
 * @param lambda      Reference wavelength.
 * @param startLambda First wavelength.
 * @param stepLambda  Step of wavelength grid.
 * 
 * @doc
 *   This function fills an image with all the extracted fiber spectra 
 *   resampled at a constant wavelength step. The resampling is made 
 *   conserving the flux locally. Signal excessive undersampling is not 
 *   allowed: this function returns a NULL if the wavelength step is
 *   more than twice the maximum value of the inverse of the first order 
 *   coefficient in the IDS for all fibers. No limits are posed to signal 
 *   oversampling: however, it is recommended to choose a wavelength step 
 *   close to, or slightly less than, the mean spectral dispersion (in 
 *   A/pixel). The resampled spectra from fiber 1 to 400 are stored in
 *   the output image from bottom to top, with the 400 spectra coming 
 *   from the pseudo-slit 0 at the bottom, and those coming from the
 *   pseudo-slit 3 at the top.
 *
 * @author C. Izzo
 */

int ifuResampleSpectra(cpl_image *image, cpl_table *spectra, cpl_table *ids, 
                       int slit, double lambda, double startLambda, 
                       double stepLambda)
{

/*  char       task[] = "ifuResampleSpectra";  */

  int        countFibers = N_BLOCKS * FIBERS_PER_BLOCK;
  int        nx          = cpl_image_get_size_x(image);
  int        ny          = cpl_image_get_size_y(image);
  float     *data        = cpl_image_get_data(image);
  float     *line        = data;
  int        order       = cpl_table_get_ncol(ids) - 3;
  int        offset      = cpl_table_get_int(spectra, "y", 0, NULL);
  int        npix        = cpl_table_get_nrow(spectra);
  int        null        = 0;
  int        intPixel;
  double    *spectrum;
  double     pixel;
  double    *p;
  double     value;
  double    *v;
  double    *c;
  cpl_table *table;
  char       colName[MAX_COLNAME];
  int        i, j, k;


  /*
   *  Allocate the work table were each resampled spectrum is constructed.
   */

  table = cpl_table_new(nx);
  /* Pixels at each lambda    */
  cpl_table_new_column(table, "pixel", CPL_TYPE_DOUBLE);  

  /* Derivative of the above  */
  cpl_table_new_column(table, "dpixel", CPL_TYPE_DOUBLE); 

  /* Interpolated spec values */
  cpl_table_new_column(table, "values", CPL_TYPE_DOUBLE); 

  p = cpl_table_get_data_double(table, "pixel");
  v = cpl_table_get_data_double(table, "values");
  c = cpl_malloc((order + 1) * sizeof(double));

  if (ny > countFibers)
    line += slit * countFibers * nx;

  for (i = 0; i < countFibers; i++, line += nx) {

    /*
     *  Get the IDS for the current fiber
     */

    for (j = 0; j <= order; j++) {
      snprintf(colName, MAX_COLNAME, "c%d", j);
      c[j] = cpl_table_get_double(ids, colName, i, &null);
      if (null)
        break;
    }
    if (null) {
      null = 0;
      continue;
    }

    /*
     *  Determine the pixel positions corresponding to each wavelength
     *  of the sampling, and linearly interpolate the signal at those 
     *  positions.
     */

    snprintf(colName, MAX_COLNAME, "s%d", i + 1);
    spectrum = cpl_table_get_data_double(spectra, colName);

    if (spectrum == NULL) {
      cpl_error_reset();
      continue;
    }

    for (k = 0; k < nx; k++) {
      pixel = modelValue1D(c, order, startLambda + k * stepLambda - lambda);
      pixel -= offset;
      p[k] = pixel;

      /*
       *  Linear interpolation of the signal
       */

      intPixel = pixel;
      if (intPixel <= 0 || intPixel > npix - 2) {
        v[k] = 0.0;
        continue;
      }

      v[k] = spectrum[intPixel] * (1 - pixel + intPixel)
           + spectrum[intPixel + 1] * (pixel - intPixel);
    }


    /*
     *  Apply flux conservation correction. Compute the pixel positions 
     *  derivative, and multiply the interpolated values by it.
     */

    cpl_table_copy_data_double(table, "dpixel", p);
    cpl_table_shift_column(table, "dpixel", -1);
    cpl_table_subtract_columns(table, "dpixel", "pixel");
    value = cpl_table_get_double(table, "dpixel", nx - 2, NULL);
    cpl_table_set_double(table, "dpixel", nx - 1, value);
    cpl_table_multiply_columns(table, "values", "dpixel");

    for (k = 0; k < nx; k++)
      line[k] = v[k];
    
  }

  cpl_table_delete(table);
  cpl_free(c);

  return 0;

}


/**
 * @memo
 *   Read all extracted spectra from table to image
 *
 * @return Image with spectra.
 *
 * @param image       Allocated image to be filled with the spectra.
 * @param spectra     Extracted arc lamp spectra for each fiber.
 * @param slit        IFU pseudo-slit (from 0 to 3, with 0 the isolated slit).
 *
 * @doc
 *   This function fills an image with all the extracted fiber spectra.
 *   The extracted spectra from fiber 1 to 400 are stored in the
 *   output image from bottom to top, with the 400 spectra coming
 *   from the pseudo-slit 0 at the bottom, and those coming from the
 *   pseudo-slit 3 at the top.
 *
 * @author C. Izzo
 */


int ifuReadSpectra(cpl_image *image, cpl_table *spectra, int slit)
{ 
    
/*  char       task[] = "ifuReadSpectra";  */
  
  int        countFibers = N_BLOCKS * FIBERS_PER_BLOCK;
  int        nx          = cpl_image_get_size_x(image);
  int        ny          = cpl_image_get_size_y(image);
  int        nrow        = cpl_table_get_nrow(spectra);
  float     *line        = cpl_image_get_data(image);
  double    *spectrum;
  char       colName[MAX_COLNAME];
  int        i, k;


  if (ny > countFibers)
    line += slit * countFibers * nx;

  for (i = 0; i < countFibers; i++, line += nx) {

    snprintf(colName, MAX_COLNAME, "s%d", i + 1);
    spectrum = cpl_table_get_data_double(spectra, colName);

    if (spectrum == NULL) {
      cpl_error_reset();
      continue;
    }

    for (k = 0; k < nrow; k++)
      line[k] = spectrum[k];

  }

  return 0;

}


/**
 * @memo
 *   Correct arc calibration with offset of sky lines.
 *
 * @return Applied shift.
 *
 * @param spectra     Extracted scientific spectra for each fiber.
 * @param ids         IDS table.
 * @param lambda      Reference wavelength.
 *
 * @doc
 *   This function determine on extracted (but still uncalibrated)
 *   spectra the distance of the skylines from their expected
 *   position according to the given IDS. The median offset is
 *   then added to the constant term of the IDS. Currently, only
 *   skylines 5577.338, 6300.304, 6363.780, and 8344.602 are used. 
 *
 * @author C. Izzo
 */
    
double ifuAlignSkylines(cpl_table *spectra, cpl_table *ids, double lambda,
                        int individual)
{

  int        countFibers  = N_BLOCKS * FIBERS_PER_BLOCK;
  int        order        = cpl_table_get_ncol(ids) - 3;
  int        offset       = cpl_table_get_int(spectra, "y", 0, NULL);
  int        npix         = cpl_table_get_nrow(spectra);
  int        null         = 0;
  int        searchRadius = 7;  /* It was 4, Peter Weilbacher */
  int        nLines       = 4;  /* Number of skylines */
  double     line[]       = {5577.338, 6300.304, 6363.780, 8344.602};
  double     expected;
  double     detected;
  double     shift;
  int        iexpected;
  double    *spectrum;
  cpl_table *table;
  double    *c;
  /* double     maxShift = 2000.0;   Removed by Peter Weilbacher */

  /* Added by Peter Weilbacher: */
  double     maxShift = 30.0,
             maxResidShift = 2.0; /* everything beyond a 2px individual shift *
                                   * means that the wavelength calibration    *
                                   * was totally wrong (cosmic ray?)          */
  char       colName[MAX_COLNAME];
  int        found;
  int        startBin, endBin, length;
  int        countNulls;
  int        i, j;

#if DEBUG_SHIFTS
  char       tablename[1024];
#endif


  /*
   *  Allocate the work table an offset is determined for each spectrum.
   */

  table = cpl_table_new(countFibers); 
  cpl_table_new_column(table, "shift", CPL_TYPE_DOUBLE);

  c = cpl_malloc((order + 1) * sizeof(double));


  /*
   *  Compute expected position for skylines, and compare it to their
   *  real position.
   */

  for (i = 0; i < countFibers; i++) {

    /*
     *  Get the IDS for the current fiber
     */
    
    for (j = 0; j <= order; j++) {
      snprintf(colName, MAX_COLNAME, "c%d", j);
      c[j] = cpl_table_get_double(ids, colName, i, &null);
      if (null)
        break;
    }
    if (null) {
      null = 0;
      continue;
    }
    
    snprintf(colName, MAX_COLNAME, "s%d", i + 1);
    spectrum = cpl_table_get_data_double(spectra, colName);

    if (spectrum == NULL) {
      cpl_error_reset();
      continue;
    }

    found = 0;
    shift = 0.0;
    for (j = 0; j < nLines; j++) {
      expected = modelValue1D(c, order, line[j] - lambda);
      iexpected = expected;


      /*
       *  Determine interval where to look for the current sky line
       */

      startBin = iexpected - offset - searchRadius;
      endBin = iexpected - offset + searchRadius;
      if (startBin < 0)
        continue;
      if (endBin > npix)
        continue;
      length = endBin - startBin;

      if (findPeak(spectrum + startBin, length, 0.0, &detected)) {
        shift += startBin + offset + detected - expected;
        found++;
      }
    }

    if (found) {
      shift /= found;
      if (shift < maxShift) {
        cpl_table_set_double(table, "shift", i, shift);
      }
    }
  }

  cpl_free((void*)c);

  countNulls = cpl_table_count_invalid(table, "shift");

  if (countNulls == countFibers)
    return 0.0;

/* cpl_table_save(table, NULL, NULL, "shifts.fits", CPL_IO_CREATE); */

#if DEBUG_SHIFTS
  /* save for debugging */
  sprintf(tablename, "align_shifts1_%s.fits", individual ? "indi" : "norm");
  cpl_table_save(table, NULL, NULL, tablename, CPL_IO_CREATE);
  sprintf(tablename, "align_ids1_%s.fits", individual ? "indi" : "norm");
  cpl_table_save(ids, NULL, NULL, tablename, CPL_IO_CREATE);
#endif

  /*  
   *  Median shift.
   *  Apply this in any case, so that even fibers with invalid values
   *  are shifted by the extra zeropoint.
   */

  shift = cpl_table_get_column_median(table, "shift");
  cpl_msg_info("ifuAlignSkylines", "Applying median shift of %f px", shift);
  cpl_table_add_scalar(ids, "c0", shift); /* apply the median shift */

  /* Added by Peter Weilbacher: */

  if (individual) {
    cpl_msg_info("ifuAlignSkylines", "Now applying individual shifts...");

    /* subtract it from shift column to make those shift into residual shifts */
    cpl_table_subtract_scalar(table, "shift", shift);
#if DEBUG_SHIFTS
    cpl_table_save(table, NULL, NULL, "align_shifts2_indi.fits", CPL_IO_CREATE);
#endif

    /* now apply individual shift if available, one by one */
    for (i = 0; i < countFibers; i++) {
      if (cpl_table_is_valid(table, "shift", i) == 1 &&
          cpl_table_is_valid(ids, "c0", i) == 1) {
        double residualshift, c0;
        int rscheck = 0, c0check = 0;

        residualshift = cpl_table_get_double(table, "shift", i, &rscheck);
#if DEBUG_SHIFTS
        cpl_msg_debug("ifuAlignSkylines()", "%d: %f", i+1, residualshift);
#endif
        if (fabs(residualshift) <= maxResidShift) {
          c0 = cpl_table_get_double(ids, "c0", i, &c0check);
          if (rscheck || c0check) {
            continue;
          }
          cpl_table_set_double(ids, "c0", i, c0 + residualshift);
        }
      }
    }
  } else {
    cpl_msg_info("ifuAlignSkylines", "NOT applying individual shifts");
#if DEBUG_SHIFTS
    cpl_table_save(ids, NULL, NULL, "align_ids2_norm.fits", CPL_IO_CREATE);
#endif
  }
  /* End of addition by Peter Weilbacher */

  cpl_table_delete(table);

  return shift;

}


/**
 * @memo
 *   Return sequence number of active fiber closest to CCD center.
 *   
 * @return Sequence number of active fiber closest to CCD center.
 * 
 * @param short_coeff Coefficients of linear tracing around reference row.
 * @param row         Reference row.
 * 
 * @doc
 *   This function finds the active fiber that is closest to the CCD 
 *   coordinate X = 1024, and returns its sequence number (counted 
 *   starting from 0). The input table contains the coefficients of 
 *   a linear fit to the spectral tracings around the reference row, 
 *   i.e., at the CCD coordinate Y = @em row. In case of error, a 
 *   negative number is returned.
 *
 * @author C. Izzo
 */

int findCentralFiber(cpl_table *short_coeff, int row)
{

  int     countFibers = cpl_table_get_nrow(short_coeff);
  int     null        = 0;
  int     fiber       = -1;
  int     i;
  double  c[2];
  double  pos, prepos;


  if (countFibers != N_BLOCKS * FIBERS_PER_BLOCK)
    return fiber;

  /*
   *  Compute all X positions of tracings at the reference row.
   */

  for (i = 0; i < countFibers; i++) {
    c[0] = cpl_table_get_double(short_coeff, "c0", i, &null);
    if (null) {
      null = 0;
      continue;
    }
    c[1] = cpl_table_get_double(short_coeff, "c1", i, NULL);
    pos = modelValue1D(c, 1, row);

    if (pos > 1024.) {
      if ((pos - 1024.) < (1024. - prepos))
        fiber = i;
      else
        fiber = i - 1;
      break;
    }
    else
      prepos = pos;

  }

  return fiber;

}


/**
 * @memo
 *   Determine sky spectrum and subtract it from the data.
 *
 * @return Sky spectrum
 *
 * @param extracted Image of extracted spectra.
 *
 * @doc
 *   This function finds the sky as the median value along the cross
 *   dispersion direction (y coordinate), and subtract it from the data.
 *   This method is appropriately applied on data where more than 50%
 *   of the extracted spectra are coming from the sky, typically on a
 *   standard star exposure.
 *
 * @author C. Izzo
 */

cpl_image *ifuSubtractSky(cpl_image *extracted)
{

  float     *data = cpl_image_get_data(extracted);
  int        nx   = cpl_image_get_size_x(extracted);
  int        ny   = cpl_image_get_size_y(extracted);
  cpl_image *sky;
  float     *skydata;
  float     *column;
  int        i, j;


  sky = cpl_image_new(nx, 1, CPL_TYPE_FLOAT);
  skydata = cpl_image_get_data(sky);

  column = cpl_malloc(ny * sizeof(float));

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++)
      column[j] = data[i + j * nx];
    skydata[i] = median(column, ny);
    for (j = 0; j < ny; j++)
      data[i + j * nx] -= skydata[i];
  }

  cpl_free(column);

  return sky;

}


/**
 * @memo
 *   Obtain total spectrum from IFU observation.
 *
 * @return Total spectrum
 *
 * @param extracted   Image of extracted spectra.
 *
 * @doc
 *   This function simply sum all the spectra and returns the sum
 *   spectrum.
 *
 * @author C. Izzo
 */

cpl_image *ifuSumSpectrum(cpl_image *extracted)
{

  float     *data = cpl_image_get_data(extracted);
  int        nx   = cpl_image_get_size_x(extracted);
  int        ny   = cpl_image_get_size_y(extracted);
  cpl_image *sum;
  float     *sumdata;
  double     value;
  int        i, j;


  sum = cpl_image_new(nx, 1, CPL_TYPE_FLOAT);
  sumdata = cpl_image_get_data(sum);

  for (i = 0; i < nx; i++) {
    value = 0.0;
    for (j = 0; j < ny; j++)
      value += data[i + j * nx];
    sumdata[i] = value;
  }

  return sum;

}


/**
 * @memo
 *   Obtain mean flux per extracted fiber per given wavelength.
 *   
 * @return 0 = success
 * 
 * @param extracted   Image of extracted spectra.
 * @param lambda      Wavelength of line to evaluate.
 * @param start       Wavelength of first pixel in image.
 * @param step        Wavelength interval per pixel.
 * @param flux        Mean flux per fiber.
 * @param flux_err    Error on mean flux per fiber.
 * 
 * @doc
 *   This function simply sums along the whole Y image dimension all ADUs 
 *   on a 11 pixel interval centered on the given wavelength along the 
 *   whole Y image dimension, and divides it by the number of image rows
 *   containing signal (value > 0.0).
 *   
 * @author C. Izzo
 */

int extractIfuFlux(cpl_image *extracted, double lambda, double start, 
                   double step, double *flux, double *flux_err)
{

  float      *data = cpl_image_get_data(extracted);
  int         nx   = cpl_image_get_size_x(extracted);
  int         ny   = cpl_image_get_size_y(extracted);
  int         pixcen, pixsta, pixend;
  int         pos;
  int         nfib;
  double     *buffer;
  double      value, rms;
  int         i, j;
  cpl_vector *vect;


  *flux = 0.0;
  *flux_err = 0.0;

  pixcen = (lambda - start) / step;

  pixsta = pixcen - 5;
  pixend = pixcen + 6;

  if (pixsta < 0 || pixend > nx)
    return 1;

  buffer = cpl_calloc(ny, sizeof(double));

  nfib = 0;
  for (i = 0; i < ny; i++) {
    value = 0.0;
    pos = i*nx + pixsta;
    for (j = pixsta; j < pixend; j++, pos++)
      value += data[pos];
    if (value > 0.0) {
      buffer[nfib] = value;
      nfib++;
    }
  }

  if (nfib < 3) {
    cpl_free(buffer);
    return 2;
  }

  vect = cpl_vector_wrap(nfib, buffer);
  value = cpl_vector_get_median_const(vect);
  *flux = cpl_vector_get_mean(vect);
  cpl_vector_unwrap(vect);

  rms = 0.0;
  for (i = 0; i < nfib; i++)
    rms += fabs(buffer[i] - value);

  cpl_free(buffer);

  rms /= nfib;
  rms *= MEANDEV_TO_SIGMA;

  *flux_err = rms;

  return 0;

}



/**
 * @memo
 *   Set to zero all background spectra.
 * 
 * @return 0 on success.
 * 
 * @param extracted   Image of extracted spectra.
 *
 * @doc
 *   This function simply set to zero all spectra that can be
 *   considered background spectra. A background spectrum is 
 *   detected under the assumption that the background level
 *   has already been subtracted: in such a case, a background
 *   spectrum is one that has enough negative values (more than
 *   a given percentage).
 *
 * @author C. Izzo
 */
 
int ifuSetZeroLevel(cpl_image *extracted)
{
 
  float     *data = cpl_image_get_data(extracted);
  int        nx   = cpl_image_get_size_x(extracted);
  int        ny   = cpl_image_get_size_y(extracted);
  int        count;
  int        i, j;
  double     ratio;


  for (j = 0; j < ny; j++) {
    count = 0;
    for (i = 0; i < nx; i++)
      if (data[i + j * nx] < 0.0)
        count++;
    ratio = ((double)count) / nx;
    if (ratio > 0.2)
    for (i = 0; i < nx; i++)
      data[i + j * nx] = 0.0;
  }

  return 0;

}

/*@}*/
