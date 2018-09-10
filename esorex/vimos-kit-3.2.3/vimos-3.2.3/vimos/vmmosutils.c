/* $Id: vmmosutils.c,v 1.4 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>
#include <math.h>

#include <pilmemory.h>
#include <pildate.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include <cpl_table.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmextractiontable.h"
#include "vmmoswavecalib.h"
#include "vmadf.h"
#include "vmfit.h"
#include "vmutils.h"
#include "vmmosutils.h"
#include "cpl.h"


/*
 * C.Izzo, ESO
 *
 * Local spectral distortion
 */

int
modelWavcal(VimosExtractionTable *extractionTable, int lorder)
{
  
  VimosExtractionSlit *slit;
  VimosDpoint         *row_coefs;     /* Storage of row,coef points */
  int                  order;         /* Order of local wavecal solution */
  double              *coefs = NULL;  /* Coefs of local fit on coefficients */
  int                  nrows;         /* Number of rows in slit */
  int                  valid;         /* Number of valid wavecal solutions */
  double               value;
  float               *buffer;
  int                  i, j, k;

  int                  exrows = 3;    /* Rows excluded at start and end */
  

  slit = extractionTable->slits;
  order = slit->invDis[0]->order;
  
  while (slit) {

    /* 
     * Get number of rows in current slit 
     */

    nrows = slit->numRows;
    row_coefs = newDpoint(nrows);

    if (lorder == 0) {
      buffer = cpl_malloc(nrows * sizeof(float));

      for (i = 0; i <= order; i++) {

        for (j = exrows, valid = 0; j < nrows - exrows; j++) {

          /*
           * Check that this row has a wavelength solution
           */

          if (slit->invDisQuality->data[j]) {
            buffer[valid] = slit->invDis[j]->coefs[i];
            valid++;
          }

        }

        if (!valid)
          break;

        value = medianPixelvalue(buffer, valid);

        for (j = 0; j < nrows; j++) {

          /*
           * If a row had not a wavelength solution, we assign the
           * model solution anyway.
           */

          if (i == order)
            slit->invDisQuality->data[j] = 1;
          slit->invDis[j]->coefs[i] = value;

        }
      }
      cpl_free(buffer);
    }
    else {

     /* 
      * Loop over all orders and all rows of current slit,
      * excluding first and last rows.
      */

      for (i = 0; i <= order; i++) {

        for (j = exrows, valid = 0; j < nrows - exrows; j++) {

         /* 
          * Check that this row has a wavelength solution 
          */

          if (slit->invDisQuality->data[j]) {
            row_coefs[valid].x = j;
            row_coefs[valid].y = slit->invDis[j]->coefs[i];
            valid++;
          }

        }

        if (!valid)
          break;

        coefs = fit1DPoly(lorder, row_coefs, valid, NULL);

        if (coefs) {
          for (j = 0; j < nrows; j++) {

           /*
            * If a row had not a wavelength solution, we assign the
            * model solution anyway.
            */

            if (i == order)
              slit->invDisQuality->data[j] = 1;

            for (k = 0, value = 0.0; k <= lorder; k++)
              value += coefs[k] * ipow(j, k);

            slit->invDis[j]->coefs[i] = value;

          }

          cpl_free(coefs);
          coefs = NULL;
        }

      }

      cpl_free(row_coefs);
    }
    slit = slit->next;
  }

  return EXIT_SUCCESS;
}

/*
 * This one should compute the mean level (in ADU) on a rectangle around 
 * a given wavelength on a spectrum corresponding to a given slit.
 * (C.Izzo, September 2005).
 *
 *   image   = image containing raw spectra.
 *   descs   = Header containing grism information
 *   slit    = slit of chosen spectrum.
 *   wave    = wavelength around which to extract.
 *   dy      = pixels to extract along the dispersion direction (radius).
 *             If dy = n, then 2n+1 pixels are extracted (n >= 0).
 *   o_flux  = Returned mean level in ADU.
 *
 * Returning EXIT_SUCCESS or EXIT_FAILURE.
 */

int
extractSpecLevel(VimosImage *image, VimosExtractionSlit *slit, 
                    double wave, int dy, double *o_flux)
{

  const char modName[] = "extractSpecLevel";

  int    minRows = 3;   /* Min number of rows required in slit */
  int    dx, sx, ex, sy, ey, sr, er;
  int    npix, xlen, ylen;
  int    r, x, y, yCcd;
  float  yCen, yOff;
  double flux;


  *o_flux = 0.0;

  if (!image)
    return EXIT_FAILURE;

  if (!slit)
    return EXIT_FAILURE;

  if (slit->numRows < minRows)
    return EXIT_FAILURE;

  if (dy < 0)
    return EXIT_FAILURE;

  dx = slit->numRows / 2;             /* Number of rows from middle slit */
  sr = dx - dx / 2;                   /* Start row                       */
  er = sr + dx;                       /* End row                         */
  sx = slit->ccdX->data[0] + sr;      /* Start row on CCD                */
  ex = sx + dx;                       /* End row on CCD                  */

  cpl_msg_debug(modName, "Extract %d rows from %d to %d", dx, sx, ex);

  xlen = image->xlen;
  ylen = image->ylen;

  flux = 0.0;
  npix = 0;


  for (x = sx, r = sr;  r < er; x++, r++) {
    if (x >= 0 && x < xlen) {
      if (slit->invDisQuality->data[r]) {
        yCen = slit->ccdY->data[r];
        yOff = computeDistModel1D(slit->invDis[r], wave);
        yCcd = floor(yCen + yOff + 0.5);
        sy = yCcd - dy;
        ey = yCcd + dy + 1;

        cpl_msg_debug(modName, "  %d: yCen = %.2f, yOff = %.2f, yCcd = %d\n", 
                    r, yCen, yOff, yCcd);

        for (y = sy; y < ey; y++) {
          if (y >= 0 && y < ylen) {
            flux += image->data[x + y * xlen];
            npix++;
          }
        }
      }
    }
  }

  if (!npix)
    return EXIT_FAILURE;

  flux /= npix;

  *o_flux = flux;

  return EXIT_SUCCESS;

}


/*
 * This one should integrate counts on a rectangle around a given 
 * wavelength on a spectrum corresponding to a given slit. Then
 * the counts are normalized to the corresponding physical area on
 * the slit (C.Izzo, March 2003).
 *
 *   image   = image containing raw spectra.
 *   descs   = Header containing grism information
 *   slit    = slit of chosen spectrum.
 *   wave    = wavelength around which to extract.
 *   dy      = pixels to extract along the dispersion direction (radius).
 *             If dy = n, then 2n+1 pixels are extracted (n >= 0).
 *   o_flux  = Returned integrated flux, in ADU/cm^2
 *   o_err   = Returned error on integrated flux.
 *
 * Returning EXIT_SUCCESS or EXIT_FAILURE.
 */

int
extractSpecFlux(VimosImage *image, VimosExtractionSlit *slit, 
                    double wave, int dy, double *o_flux, double *o_err)
{

  const char modName[] = "extractSpecFlux";

  int    minRows = 7;   /* Min number of rows required in slit */
  int    dx, sx, ex, sy, ey, sr, er;
  int    npix, xlen, ylen;
  int    r, x, y, yCcd;
  int    saturation = 60000; /* Used values must be less than so many ADUs */
  float  yCen, yOff;
  double flux, error, area;


  *o_flux = 0.0;
  *o_err = 0.0;

  if (!image)
    return EXIT_FAILURE;

  if (!slit)
    return EXIT_FAILURE;

  if (slit->numRows < minRows)
    return EXIT_FAILURE;

  if (dy < 0)
    return EXIT_FAILURE;

  dx = slit->numRows / 2;             /* Number of rows from middle slit */
  sr = dx - dx / 2;                   /* Start row                       */
  er = sr + dx;                       /* End row                         */
  sx = slit->ccdX->data[0] + sr;      /* Start row on CCD                */
  ex = sx + dx;                       /* End row on CCD                  */

  cpl_msg_debug(modName, "Extract %d rows from %d to %d", dx, sx, ex);

  xlen = image->xlen;
  ylen = image->ylen;

  flux = 0.0;
  npix = 0;


  for (x = sx, r = sr;  r < er; x++, r++) {
    if (x >= 0 && x < xlen) {
      if (slit->invDisQuality->data[r]) {
        yCen = slit->ccdY->data[r];
        yOff = computeDistModel1D(slit->invDis[r], wave);
        yCcd = floor(yCen + yOff + 0.5);
        sy = yCcd - dy;
        ey = yCcd + dy + 1;

        cpl_msg_debug(modName, "  %d: yCen = %.2f, yOff = %.2f, yCcd = %d\n", 
                    r, yCen, yOff, yCcd);

        for (y = sy; y < ey; y++) {
          if (y >= 0 && y < ylen) {
            if (image->data[x + y * xlen] < saturation) {
              flux += image->data[x + y * xlen];
              npix++;
            }
          }
        }
      }
    }
  }

  if (!npix)
    return EXIT_FAILURE;

  error = sqrt(flux);


  /*
   * Flux correction for lost pixels
   */

  flux *= dx * (2*dy + 1) / (float)npix;
  error *= dx * (2*dy + 1) / (float)npix;

  area = (slit->maskX->data[er] - slit->maskX->data[sr]) * slit->width;

  flux /= area;
  error /= area;

  *o_flux = flux;
  *o_err = error;

  return EXIT_SUCCESS;

}


/*
 * Given wavelength, determine spectral resolution from a given arc 
 * lamp line. A high S/N is assumed (virtually no noise). Working 
 * with 2D extracted images!
 */

int
spectralResolution(VimosImage *image, float lambda, double *meanValue,
                   double *rmsValue, int saturation)
{

  int     i, j, n, m;
  int     position, maxpos;
  int     xlen, ylen;
  int     sp, ep;
  int     radius = 5;         /* Search radius for max */
  int     threshold = 500;    /* Peak must be at least so many ADUs above min */

  int     ifwhm;
  float   fwhm;
  float  *buffer;
  float   min, max, halfmax;
  float   cut = 1.5;         /* To cut outliers from FWHM values (pixel) */

  double  start, step;
  double  value, rms;


  *meanValue = 0.0;
  *rmsValue = 0.0;
  
  xlen = image->xlen;
  ylen = image->ylen;

  buffer = cpl_malloc(ylen * sizeof(double));

  readDoubleDescriptor(image->descs, pilTrnGetKeyword("Crval", 1), &start, 0);
  readDoubleDescriptor(image->descs, pilTrnGetKeyword("Cdelt", 1), &step, 0);


  /*
   *  Closest pixel to specified wavelength.
   */

  position = floor((lambda - start) / step + 0.5);


  /*
   *  Search interval for peak. Abort if too close to image border.
   */

  sp = position - radius;
  ep = position + radius;

  if (sp < 0 || ep > xlen)
    return EXIT_FAILURE;

  for (i = 0, n = 0; i < ylen; i++) {    /*  For each row of each slit  */

    /*
     *  Determine min-max value and position.
     */

    maxpos = sp;
    min = max = image->data[sp + i * xlen];
    for (j = sp; j < ep; j++) {
      if (image->data[j + i * xlen] > max) {
        max = image->data[j + i * xlen];
        maxpos = j;
      }
      if (image->data[j + i * xlen] < min) {
        min = image->data[j + i * xlen];
      }
    }

    if (fabs(min) < 0.0000001)
      continue;                       /* Truncated spectrum */

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
        if (image->data[j + i * xlen] < halfmax) {
          fwhm = ifwhm + (image->data[j - 1 + i * xlen] - halfmax)
               / (image->data[j - 1 + i * xlen] - image->data[j + i * xlen]);
          break;
        }
        ifwhm++;
      }
    }

    ifwhm = 0;
    for (j = maxpos; j > maxpos - radius; j--) {
      if (j >= 0) {
        if (image->data[j + i * xlen] < halfmax) {
          fwhm += ifwhm + (image->data[j + 1 + i * xlen] - halfmax)
               / (image->data[j + 1 + i * xlen] - image->data[j + i * xlen]);
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
    return EXIT_FAILURE;
  }

  value = medianPixelvalue(buffer, n);

  rms = 0.0;
  for (i = 0, m = 0; i < n; i++) {
    if (fabs(buffer[i] - value) < cut) {
      rms += fabs(buffer[i] - value);
      m++;
    }
  }

  cpl_free(buffer);

  if (m < 3)
    return EXIT_FAILURE;

  rms /= m;
  rms *= MEANDEV_TO_SIGMA;

  value *= step;
  rms *= step;

  *meanValue = lambda / value;
  *rmsValue = *meanValue * rms / value;

  return EXIT_SUCCESS;

}


/*
 * Around each line expected position check the presence of saturated 
 * pixels on all spectra. Allow a small percentage of saturated pixels 
 * for cosmics, zero order contamination, hot pixels, and the kind.
 * If this threshold is NOT passed for any of the arc lamp lines, this
 * function returns zero.
 */

int testLineSaturation(VimosImage *image, VimosTable *lineCat)
{

  const char modName[] = "testLineSaturation";

  int xlen     = image->xlen;
  int ylen     = image->ylen;
  int numLines = lineCat->cols->len;
  int radius   = 3;
  int found    = 0;
  int saturated, totPix;
  int cpix, npix;
  int sp, ep;
  int i, j, k;

  float *profile;
  float  expectPeak;
  float  lambda;

  float  saturationLevel = 65000;
  double tolerance = 0.2;

  double start, step;

  VimosColumn *wLen = findColInTab(lineCat, "WLEN");


  readDoubleDescriptor(image->descs, pilTrnGetKeyword("Crval", 1), &start, 0);
  readDoubleDescriptor(image->descs, pilTrnGetKeyword("Cdelt", 1), &step, 0);

  npix = 2 * radius + 1;
  profile = cpl_calloc(npix, sizeof(float));

  for (i = 0; i < numLines; i++) {

    totPix = 0;
    saturated = 0;

    /*
     *  Expected peak and closest pixel to specified wavelength.
     */

    lambda = wLen->colValue->fArray[i];
    expectPeak = (lambda - start) / step;
    cpix = floor(expectPeak + 0.5);

    /*
     *  Search interval for peak. Abort if too close to image border.
     */

    sp = cpix - radius;
    ep = cpix + radius;

    if (sp < 0 || ep > xlen)
      continue;

    for (j = 0; j < ylen; j++) {    /*  For each row of each slit  */
      for (k = 0; k < npix; k++) {
        if (image->data[sp + k + j * xlen] > MIN_DIVISOR)
          totPix++;
        if (image->data[sp + k + j * xlen] > saturationLevel)
          saturated++;
      }
    }

    if (totPix) {
      if ((double)saturated / totPix > tolerance) {
        cpl_msg_info(modName, "Line %7.2f (X = %d): SATURATED", lambda, cpix);
        found = 1;
      }
      else
        cpl_msg_debug(modName, "Line %7.2f (X = %d): ok", lambda, cpix);
    }
    else {
      cpl_msg_debug(modName, "Line %7.2f (X = %d): ok (not in spectral range)", 
                   lambda, cpix);
    }
  }

  return found;

}


/*
 * Estimate the spectral distortion modelling goodness. This is done
 * computing the rms of the residuals between the expected positions
 * of the arc lamp lines (according to the WCS of the 2D extracted arc
 * lamp image), and the actual position of a detected peak. The peak
 * is searched within the range specified by the last argument (given
 * in Angstrom), which typically would depend on the expected lines
 * FWHM.
 */

double
distortionsRms(VimosImage *image, VimosTable *lineCat, double fwhm) 
{

  const char modName[] = "distortionsRms";

  int xlen     = image->xlen;
  int ylen     = image->ylen;
  int numLines = lineCat->cols->len;
  int radius;
  int cpix, npix, nzero;
  int sp, ep; 
  int i, j, k;
  int npeaks, allPeaks;

  float *profile;
  float  peak, expectPeak, offset;
  float  lambda;

  double  start, step;
  double  average;
  double  rms, oneRms;

  VimosColumn *wLen = findColInTab(lineCat, "WLEN");


  readDoubleDescriptor(image->descs, pilTrnGetKeyword("Crval", 1), &start, 0);
  readDoubleDescriptor(image->descs, pilTrnGetKeyword("Cdelt", 1), &step, 0);

  radius = ceil(fwhm / step);
  npix = 2 * radius + 1;
  profile = cpl_calloc(npix, sizeof(float));

  rms = 0.0;
  allPeaks = 0;
  
  for (i = 0; i < numLines; i++) {

    /*
     *  Expected peak and closest pixel to specified wavelength.
     */ 
   
    lambda = wLen->colValue->fArray[i];
    expectPeak = (lambda - start) / step;
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
        profile[k] = image->data[sp + k + j * xlen];
        if (fabs(profile[k]) < MIN_DIVISOR)
          nzero++;    /* Count number of zero pixels (spectrum truncated) */
      }
      if (nzero > 0)
        continue;
      if (findPeak1D(profile, npix, &peak, 2) == VM_TRUE) {
        offset = (sp + peak) - expectPeak - 0.5;
        average += offset;
        rms += fabs(offset);
        oneRms += fabs(offset);
        npeaks++;
        allPeaks++;
      }
    }

    if (npeaks)
      cpl_msg_info(modName, "RMS for %.2f: %.3f", 
                  lambda, oneRms / npeaks * MEANDEV_TO_SIGMA);
    else
      cpl_msg_info(modName, "RMS for %.2f: line not available", lambda);
/*
    printf("XXXX %.2f %.4f\n", lambda, average / npeaks);
*/
  }

  cpl_free(profile);

  if (allPeaks < 10)
    return 0.0;

  rms /= allPeaks;
  rms *= MEANDEV_TO_SIGMA;

  return rms;

}


/*
 * This is the same as distortionsRms(), but the input catalog is a
 * CPL table. This function is just a temporary solution in these
 * intermediate pipeline status between the VIMOS libraries and the
 * CPL libraries.
 *
 * Estimate the spectral distortion modelling goodness. This is done
 * computing the rms of the residuals between the expected positions
 * of the arc lamp lines (according to the WCS of the 2D extracted arc
 * lamp image), and the actual position of a detected peak. The peak
 * is searched within the range specified by the last argument (given
 * in Angstrom), which typically would depend on the expected lines
 * FWHM. If no peak is detected, the offset is set in that case to
 * a "punitive" value corresponding to twice the last argument value. 
 */

double
distortionsRms_CPL(VimosImage *image, cpl_table *lineCat, double fwhm) 
{

  const char modName[] = "distortionsRms";

  int xlen     = image->xlen;
  int ylen     = image->ylen;
  int numLines = cpl_table_get_nrow(lineCat);
  int radius;
  int cpix, npix, nzero;
  int sp, ep; 
  int i, j, k;
  int npeaks, allPeaks;

  float *profile;
  float  peak, expectPeak, offset;
  float  lambda;

  double  start, step;
  double  average;
  double  rms, oneRms;

  float  *wdata = cpl_table_get_data_float(lineCat, "WLEN");


  readDoubleDescriptor(image->descs, pilTrnGetKeyword("Crval", 1), &start, 0);
  readDoubleDescriptor(image->descs, pilTrnGetKeyword("Cdelt", 1), &step, 0);

  radius = ceil(fwhm / step);
  npix = 2 * radius + 1;
  profile = cpl_calloc(npix, sizeof(float));

  rms = 0.0;
  allPeaks = 0;

  for (i = 0; i < numLines; i++) {

    /*
     *  Expected peak and closest pixel to specified wavelength.
     */ 
   
    lambda = wdata[i];
    expectPeak = (lambda - start) / step;
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
        profile[k] = image->data[sp + k + j * xlen];
        if (fabs(profile[k]) < MIN_DIVISOR)
          nzero++;    /* Count number of zero pixels (spectrum truncated) */
      }
      if (nzero > 0)
        continue;
      if (findPeak1D(profile, npix, &peak, 2) == VM_TRUE) {
        offset = (sp + peak) - expectPeak;   /*    - 0.5;     */
        average += offset;
        rms += fabs(offset);
        oneRms += fabs(offset);
        npeaks++;
        allPeaks++;
      }
    }

    if (npeaks)
      cpl_msg_info(modName, "RMS for %.2f: %.3f", 
                  lambda, oneRms / npeaks * MEANDEV_TO_SIGMA);
    else
      cpl_msg_info(modName, "RMS for %.2f: line not available", lambda);
/* +/
    printf("XXXX %.2f %.4f\n", lambda, average / npeaks);
/+ */
  }

  cpl_free(profile);

  if (allPeaks < 10)
    return 0.0;

  rms /= allPeaks;
  rms *= MEANDEV_TO_SIGMA;

  return rms;

}


/*
 * Just identify the grism used, and return a code:
 *  0 = LR_red
 *  1 = LR_blue
 *  2 = MR
 *  3 = HR_red
 *  4 = HR_orange
 *  5 = HR_blue
 *  6 = HR_red holographic
 *  7 = HR_blue holographic
 * -1 = Unidentified grism
 */

int
getGrism(VimosImage *image) 
{
  
  int  quadrant;  
  char grismName[80];
  char grismId[80];

  
  readIntDescriptor(image->descs, pilTrnGetKeyword("Quadrant"),
                    &quadrant, NULL);
  
  readStringDescriptor(image->descs,
                       pilTrnGetKeyword("GrismName", quadrant),
                       grismName, NULL);

  readStringDescriptor(image->descs,
                       pilTrnGetKeyword("GrismId", quadrant),
                       grismId, NULL);

  
  if (grismName[0] == 'L') 
  {  
    if (grismName[3] == 'r')  
      return 0;                     /* LR_red    */
    if (grismName[3] == 'b')
      return 1;                     /* LR_blue   */
  }
  
  if (grismName[0] == 'M')
    return 2;                       /* MR        */
  
  if (grismName[0] == 'H') 
  {
    if (grismName[3] == 'r') 
    {
      if (grismId[8] == 'H')
        return 6;                   /* HR_red holographic   */
      else
        return 3;                   /* HR_red    */
    }
    if (grismName[3] == 'o')
      return 4;                     /* HR_orange */
    if (grismName[3] == 'b')
    {
      if (grismId[9] == 'H')
        return 7;                    /* HR_blue  holographic */
      else
        return 5;                    /* HR_blue   */
    }
  }

  return -1;

}


int
getGrismAgain(VimosTable *table)
{

  int  quadrant;
  char grismName[10];
  char grismId[80];


  readIntDescriptor(table->descs, pilTrnGetKeyword("Quadrant"),
                    &quadrant, NULL);

  readStringDescriptor(table->descs,
                       pilTrnGetKeyword("GrismName", quadrant),
                       grismName, NULL);

  readStringDescriptor(table->descs,
                       pilTrnGetKeyword("GrismId", quadrant),
                       grismId, NULL);

  if (grismName[0] == 'L') {
    if (grismName[3] == 'r')
      return 0;                     /* LR_red    */
    if (grismName[3] == 'b')
      return 1;                     /* LR_blue   */
  }

  if (grismName[0] == 'M')
    return 2;                       /* MR        */

  if (grismName[0] == 'H') {
    if (grismName[3] == 'r') {
      if (grismId[8] == 'H')
        return 6;                   /* HR_red holographic   */
      else
        return 3;                   /* HR_red    */
    }
    if (grismName[3] == 'o')
      return 4;                     /* HR_orange */
    if (grismName[3] == 'b')
      return 5;                     /* HR_blue   */
  }

  return -1;

}


/*
 * Derive the X and Y components of the Optical Distortion Model using
 * an arc lamp frame.
 */

int
findCentralPosition(VimosImage *image, VimosDescriptor *desc, 
                    double X, double Y, double slitLength, float width,
                    VimosTable *lineCat,
                    double *deltaX, double *deltaY)
{

  int      nlines;
  double  *lines;
  int      npeaks;
  double  *peaks;
  int      length = slitLength;
  float    a2pix;
  float    mm2pix;
  float    pixWidth;
  float   *region;
  double   min_disp, max_disp;
  double   max, value;
  double   xTol = 1000;
  double   yTol = 1000;
  double **out;
  float    refWave;
  int      maxpos;
  int      nident;
  int      found;
  int      i, j;

  int      pixelAbove, pixelBelow;

  int      xStart;
  int      nx;

  int      yStart;
  int      yOffset;
  int      ny;

  int      c_length;
  int      c_yStart;
  int      c_yOffset;
  int      c_ny;

  VimosColumn         *wLen;


  /*
   * The spatial extension of the slit must be entirely contained in the CCD.
   */

  xStart = X - length;
  nx = 3 * length;

  if (xStart < 0 || xStart + nx >= image->xlen || Y < 0 || Y >= image->ylen) {
    *deltaX = 0;
    *deltaY = 0;
    return 1;
  }


  /*
   * Define the collapse window, and collapse along dispersion to obtain
   * spatial profile.
   */

  c_length  = 400;
  c_yStart  = Y - c_length / 2;
  c_yOffset = 0;
  c_ny      = c_length;

  if (c_yStart < 0) {
    c_ny += c_yStart;
    c_yStart = 0;
  }

  if (c_yStart + c_ny >= image->ylen)
    c_ny = image->ylen - c_yStart;

  region = collapse2Dto1D(image, xStart, c_yStart, nx, c_ny, ROW);


  /*
   * Normalisation of the spatial profile.
   */

  max = region[0];
  for (i = 1; i < nx; i++)
    if (max < region[i])
      max = region[i];

  for (i = 0; i < nx; i++) {
    region[i] /= max;
  }


  /*
   * Find the rough position of the couple of edges.
   */

  max = 0;
  maxpos = 0;
  for (i = 1, j = length + 1; i < 2 * length; i++, j++) {
    value = (region[i] - region[i-1]) * (region[j-1] - region[j]);
    if (max < value) {
      max = value;
      maxpos = i - 1;
    }
  }

  cpl_free(region);


  /*
   * Here is the offset from the expected X:
   */

  if (fabs(maxpos - length) > xTol) {
    *deltaX = 0;
    *deltaY = 0;
    return 1;
  }

  *deltaX = maxpos - length;


  /*
   * This is the center of the newly positioned spatial profile:
   */

  xStart = X + maxpos - length / 2;


  /*
   * Now define the interval of extraction
   */

  readIntDescriptor(desc, pilTrnGetKeyword("NumPixBelow"), &pixelBelow, NULL);
  readIntDescriptor(desc, pilTrnGetKeyword("NumPixAbove"), &pixelAbove, NULL);

  yStart = Y - pixelBelow;
  yOffset = 0;
  ny = pixelAbove + pixelBelow + 1;

  if (yStart < 0) {
    yOffset = -yStart;
    ny += yStart;
    yStart = 0;
  }

  if (yStart + ny >= image->ylen)
    ny = image->ylen - yStart;

  region = extractFloatImage(image->data, image->xlen, image->ylen,
                             xStart, yStart, 1, ny);

  readFloatDescriptor(desc, pilTrnGetKeyword("WlenCen"), &refWave, NULL);
  readFloatDescriptor(desc, pilTrnGetKeyword("OptDistY", 0, 1), &mm2pix, NULL);
  readFloatDescriptor(desc, pilTrnGetKeyword("Dispersion", 1, 0, 0), &a2pix, 
                      NULL);

  max_disp = min_disp = 1 / a2pix;
  max_disp += max_disp / 6;
  min_disp -= min_disp / 6;

  pixWidth = mm2pix * width;

  peaks = collectPeaks(region, ny, 200, pixWidth, &npeaks);

  cpl_free(region);

  if (npeaks == 0) {
    *deltaX = 0;
    *deltaY = 0;
    return 1;
  }

  nlines = lineCat->cols->len;
  wLen = findColInTab(lineCat, "WLEN");
  lines = cpl_malloc(nlines * sizeof(double));

  for (j = 0; j < nlines; j++)
    lines[j] = wLen->colValue->fArray[j];

  out = identPeaks(peaks, npeaks, lines, nlines, 
                   min_disp, max_disp, .1, &nident);

  cpl_free(peaks);
  cpl_free(lines);

  if (out == NULL) {
    *deltaX = 0;
    *deltaY = 0;
    return 1;
  }

  found = 0;
  for (i = 0; i < nident; i++) {
    if (fabs(out[1][i] - refWave) < 1.) {
      *deltaY = yStart + out[0][i] - Y;
      if (fabs(*deltaY) > yTol) {
        *deltaX = 0;
        *deltaY = 0;
        return 1;
      }
      found = 1;
      break;
    }
  }

  cpl_free(out[0]);
  cpl_free(out[1]);
  cpl_free(out);

  if (found)
    return 0;

  *deltaX = 0;
  *deltaY = 0;

  return 1;

}

int
alignWavePattern(VimosImage *image, double X, double Y, double slitLength, 
                 double *deltaX, double *deltaY) 
{

  int     i, j, k, n, m;
  int     grism;
  int     xlen    = image->xlen;
  int     ylen    = image->ylen;
  int     xLength = slitLength + 1;
  double  dX      = 0.0;
  double  dY      = 0.0;

  double  d, disp;
  double  refWave;
  double *listWave;
  double *listCcd;
  int     nWave;

  double  x, y, xFrac, yFrac;
  int     xStart, xInt, yInt;

  int     nDisp   = 20;      /* 100 Number of different dispersions to try  */
  double  dRange  = 0.2;      /* 0.1 Range of dispersions to try             */
  double  dStep   = dRange / nDisp;

  int     nOff    = 50;       /* 50 Number of X, and then Y, offsets to try */
  double  oRange  = 20.0;     /* 20.0 Range of offsets                        */
  double  oStep   = oRange / nOff;

  double  value[7];
  double  index, oldIndex;


  grism = getGrism(image);

  switch (grism) {
  case 0:                                                   /* LR_red */
    nWave = 4;
    listCcd = (double *)malloc(nWave * sizeof(double));
    listWave = (double *)malloc(nWave * sizeof(double));
    listWave[0] = 7383.980;
    listWave[1] = 7507.000;
    listWave[2] = refWave = 7635.105;
    listWave[3] = 7723.800;
    disp = 0.141;  /* pix/Angstrom  */
    break;
  case 1:                                                   /* LR_blue */
    nWave = 3;
    listCcd = (double *)malloc(nWave * sizeof(double));
    listWave = (double *)malloc(nWave * sizeof(double));
    listWave[0] = 4713.143;
    listWave[1] = 4921.929;
    listWave[2] = refWave = 5015.675;
    disp = 0.190;
    break;
  case 2:                                                   /* MR */
    nWave = 3;
    listCcd = (double *)malloc(nWave * sizeof(double));
    listWave = (double *)malloc(nWave * sizeof(double));
    listWave[0] = 7383.980;
    listWave[1] = refWave = 7635.105;
    listWave[2] = 7723.800;
    disp = 0.390;
    break;
  case 3:                                                   /* HR_red */
    nWave = 3;
    listCcd = (double *)malloc(nWave * sizeof(double));
    listWave = (double *)malloc(nWave * sizeof(double));
    listWave[0] = refWave = 7948.175;
    listWave[1] = 8006.156;
    listWave[2] = 8014.786;
    disp = 1.576;
    break;
  case 4:                                                   /* HR_orange */
    nWave = 5;
    listCcd = (double *)malloc(nWave * sizeof(double));
    listWave = (double *)malloc(nWave * sizeof(double));
    listWave[0] = 6678.200;
    listWave[1] = 6717.043;
    listWave[2] = refWave = 6929.468;
    listWave[3] = 6965.430;
    listWave[4] = 7032.413;
    disp = 1.550;
    break;
  case 5:                                                   /* HR_blue */
    nWave = 3;
    listCcd = (double *)malloc(nWave * sizeof(double));
    listWave = (double *)malloc(nWave * sizeof(double));
    listWave[0] = 5852.488;
    listWave[1] = refWave = 5875.618;
    listWave[2] = 5944.834;
    disp = 1.670;
    break;
  case 6:                                        /* HR_red holographic */
    nWave = 3;
    listCcd = (double *)malloc(nWave * sizeof(double));
    listWave = (double *)malloc(nWave * sizeof(double));
    listWave[0] = refWave = 7948.175;
    listWave[1] = 8006.156;
    listWave[2] = 8014.786;
    disp = 1.627;
    break;
  default:
    return EXIT_FAILURE;
  }

  /*
   * We begin with the wavelengths converted to CCD positions.
   * The reference wavelength position corresponds to the input X, Y.
   * (X,Y) are the coordinates of the beginning of the slit. We try
   * a grid of (X + dX, Y + dY) for the slit start, at different
   * dispersions.
   */

  oldIndex = 0.0;
  for (i = 0; i < nDisp; i++) {
    d = disp - dRange / 2 + i * dStep;              /* Try this dispersion  */
    for (n = 0; n < nWave; n++)
      listCcd[n] = Y + d * (listWave[n] - refWave); /* From wave to Y pixel */
    for (j = 0; j < nOff; j++) {
      dX = j * oStep - oRange / 2.0;
      x = X + dX;
      xStart = floor(x);
      xFrac = x - xStart;
      for (k = 0; k < nOff; k++) {
        dY = k * oStep - oRange / 2.0;
        index = 0.0;
        for (n = 0; n < nWave; n++) {
          y = listCcd[n] + dY;
          yInt = floor(y);
          yFrac = y - yInt;
          if (yInt < 0 || (yInt + 1) >= ylen || 
            xStart < 0 || xStart + xLength >= xlen) {
            free(listWave);
            free(listCcd);
            *deltaX = 0.0;
            *deltaY = 0.0;
            return EXIT_FAILURE;
          }
          for (m = 0, xInt = xStart; m < xLength; m++, xInt++) {
            value[1] = image->data[xInt + yInt * xlen];
            value[2] = image->data[xInt + 1 + yInt * xlen];
            value[3] = image->data[xInt + (yInt + 1) * xlen];
            value[4] = image->data[xInt + 1 + (yInt + 1) * xlen];
            value[5] = (1 - xFrac) * value[1] + xFrac * value[2];
            value[6] = (1 - xFrac) * value[3] + xFrac * value[4];
            value[0] = (1 - yFrac) * value[5] + yFrac * value[6];
            index += value[0];
          }
        }
        if (index > oldIndex) {
          oldIndex = index;
          *deltaX = dX;
          *deltaY = dY;
        }
      }
    }
  }

  free(listWave);
  free(listCcd);

  return EXIT_SUCCESS;

}


void
findSpectrumBorders(VimosFloatArray *profile, double *upper, double *lower,
                    int fuzz)
{
  int i;
  float *window;
  int    windowSize;
  float  fupper, flower;

 /*
  *  A positive edge is expected at the beginning, and a negative
  *  edge at the end of the profile.
  */

  windowSize = 2*fuzz + 1;
  window = (float *) cpl_malloc(windowSize*sizeof(float));

 /*  Look for lower edge */

  for (i=0; i<windowSize; i++) 
    window[i] = profile->data[i];

  if (findUpJump(window, windowSize, &flower, 1) == VM_FALSE) {
/*  if (findJump(window, windowSize, &flower, 1) == VM_FALSE) { */
    *lower = -999.0;
/* printf("0. "); */
  }
  else {
    *lower = flower;
/* printf("%f ",flower - fuzz); */
  }

 /*  Look for upper edge */

  for (i=0; i<windowSize; i++) 
    window[i] = profile->data[profile->len - windowSize + i];

  if (findDownJump(window, windowSize, &fupper, 1) == VM_FALSE) {
/*  if (findJump(window, windowSize, &fupper, 1) == VM_FALSE) { */
    *upper = -999.0;
/* printf("0.\n"); */
  }
  else {
    *upper = profile->len - windowSize + fupper;
/* printf("%f\n",fupper - fuzz); */
  }
  cpl_free(window);
}


/**
 * @memo
 *   Sort a list of spectroscopic frames into groups, according to 
 *   shutter position.
 *
 * @return Number of images in each images group of the rearranged
 *   input list.
 *
 * @param image        Array of input images.
 * @param imageCount   Number of input images.
 * @param shutPosCount Number of different shutter positions found (returned).
 *
 * @doc 
 *   The pointers to image in the input list are rearranged within
 *   the input list. The input frames are assumed to be all from the 
 *   same quadrant.
 *
 * @author C.Izzo
 *
 */

int *
sortByShutterPosition(VimosImage **image, int imageCount, 
                      int *shutPosCount)
{
  const char modName[] = "sortByShutterPosition";

  int          i, j, pos, error;
  int          quadrant;
  int         *shutPositionGroup = NULL;

  float       *shutLow = NULL;
  float       *shutHigh = NULL;

  double       tolerance = 1.;   /* Tolerance whithin which two shutter *
                                  * positions are considered equal      */
  int          nGroups = 0;
  int         *imageCountInGroup = NULL;

  VimosImage **holder;
  VimosMshuMode mshuMode;

  char         mshuString[80];


 /*
  * Prepare containers for high and low shutter position, and for
  * the codes used to identify identical shutter positions.
  */

  shutLow  = (float *)cpl_malloc(imageCount * sizeof(float));
  shutHigh = (float *)cpl_malloc(imageCount * sizeof(float));
  shutPositionGroup = (int *)cpl_malloc(imageCount * sizeof(int));

  if (!(shutLow && shutHigh && shutPositionGroup)) {
    cpl_msg_debug(modName, "Not enough memory");
    cpl_free(shutLow);
    cpl_free(shutHigh);
    cpl_free(shutPositionGroup);
    return NULL;
  }

  for (i = 0; i < imageCount; i++) {
    shutPositionGroup[i] = -1;
  }

  if (readIntDescriptor(image[0]->descs, pilTrnGetKeyword("Quadrant"), 
                        &quadrant, NULL) == VM_FALSE) {
    cpl_msg_debug(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("Quadrant"));
    cpl_free(shutLow);
    cpl_free(shutHigh);
    cpl_free(shutPositionGroup);
    return NULL;
  }

  /* find out if mask shutters are being used or not */
  mshuMode = VM_MSHU_MODE_UNDEF;
  if (readStringDescriptor(image[0]->descs,
                           pilTrnGetKeyword("MshuMode", quadrant), 
                           mshuString, NULL) == VM_FALSE) {
    cpl_msg_debug(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("MshuMode", quadrant));
    cpl_free(shutLow);
    cpl_free(shutHigh);
    cpl_free(shutPositionGroup);
    return NULL;
  }
  if (!strncmp(mshuString, "ON", 2))
      mshuMode = VM_MSHU_MODE_ON;
  if (!strncmp(mshuString, "OFF", 3))
      mshuMode = VM_MSHU_MODE_OFF;

  if (mshuMode == VM_MSHU_MODE_ON)
  {
    /*
     * In case they are being used, 
     * Read in start and end shutter positions for each input image
     */

    for (i = 0; i < imageCount; i++) {
      error = 1;
      if (readFloatDescriptor(image[i]->descs, 
                              pilTrnGetKeyword("MshuPosL", quadrant), 
                              shutLow + i, NULL) == VM_TRUE) {

        if (readFloatDescriptor(image[i]->descs, 
                                pilTrnGetKeyword("MshuPosH", quadrant), 
                                shutHigh + i, NULL) == VM_TRUE) {
          error = 0;
        }
        else 
          cpl_msg_debug(modName, "Cannot read descriptor %s", 
                      pilTrnGetKeyword("MshuPosH", quadrant));
      }
      else
        cpl_msg_debug(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("MshuPosL", quadrant));

      if (error) {
        cpl_free(shutLow);
        cpl_free(shutHigh);
        cpl_free(shutPositionGroup);
        return NULL;
      }
    }
    

 /*
  * Determine number of different shutter positions, and flag each
  * one of them with an internal code (i.e., and integer number 
  * starting from 0).
  */

    for (i = 0, nGroups = 0; i < imageCount; i++) {

      if (shutPositionGroup[i] < 0) {
        shutPositionGroup[i] = nGroups;
        for (j = i + 1; j < imageCount; j++) {
          if (fabs(shutLow[i] - shutLow[j]) < tolerance) {
            if (fabs(shutHigh[i] - shutHigh[j]) < tolerance) {
              shutPositionGroup[j] = nGroups;
            }
          }
        }
        nGroups++;
      }

    }

 /*
  * Count how many images are belonging to each group
  */

    imageCountInGroup = (int *)cpl_calloc(nGroups, sizeof(int));

    if (!imageCountInGroup) {
      cpl_msg_debug(modName, "Not enough memory");
      cpl_free(shutLow);
      cpl_free(shutHigh);
      cpl_free(shutPositionGroup);
      return NULL;
    }

    for (i = 0; i < imageCount; i++) {
      for (j = 0; j < nGroups; j++) {
        if (shutPositionGroup[i] == j) {
          imageCountInGroup[j]++;
          break;
        }
      }
    }

    /*
     * Finally, resort array of input images according to group
     */

    holder = (VimosImage **)cpl_malloc(imageCount * sizeof(VimosImage *));

    if (!holder) {
      cpl_msg_debug(modName, "Not enough memory");
      cpl_free(shutLow);
      cpl_free(shutHigh);
      cpl_free(shutPositionGroup);
      cpl_free(imageCountInGroup);
      return NULL;
    }

    for (i = 0, pos = 0; i < nGroups; i++) {
      for (j = 0; j < imageCount; j++) {
        if (shutPositionGroup[j] == i) {
          holder[pos] = image[j];
          pos++;
        }
      }
    }

    for (i = 0; i < imageCount; i++)
      image[i] = holder[i];
    
    cpl_free(holder);
  }
  else {
    /*
     * This is the case when mshuMode == VM_MSHU_MODE_OFF, i.e.
     * no mask shutters
     */
    nGroups = 1;
    imageCountInGroup = (int *)cpl_calloc(nGroups, sizeof(int));
    imageCountInGroup[0] = imageCount;
  }

  /*
   * Cleanup
   */
  
  cpl_free(shutLow);
  cpl_free(shutHigh);
  cpl_free(shutPositionGroup);

 /*
  * Returned:
  */

  *shutPosCount = nGroups;

  return imageCountInGroup;
  
}


/**
 * @memo
 *   Create the PAF file containing the spectral curvature and the optical 
 *   distorsion models for a given quadrant.
 *
 * @return Name of created PAF file, or NULL in case of failure.
 *
 * @param descs      List of descriptors containing the spectral curvature
 *                   and the optical distorsion coefficients.
 * @param namePAF    PAF file root name.
 *
 * @doc
 *   The PAF file name depends on the quadrant.
 *
 * @author C. Izzo
 */

char *
createSpectralDistModelsPAF(VimosDescriptor *descs, char *namePAF)
{
  const char       modName[] = "createSpectralDistModelsPAF";

  char            *pafName;
  int              ord, xOrd, yOrd;
  int              i, j, k, len, quad;
  FILE            *fp;
  VimosDescriptor *desc;

  cpl_msg_debug(modName, "Write spectral distorsion models into PAF file");

  readIntDescriptor(descs, pilTrnGetKeyword("Quadrant"), &quad, NULL);

  len = sizeof(char) * (strlen(namePAF)+7);
  pafName = (char *)cpl_malloc(len);
  if (!pafName)
    return NULL;

  sprintf(pafName, "%s-%d.paf", namePAF, quad);
  fp = fopen(pafName, "w");
  if (!fp) {
    cpl_free(pafName);
    return NULL;
  }


 /*
  * Write PAF header
  */

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafHeaderStart"),NULL);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafType"),"Configuration");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafId"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafName"), pafName);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafDesc"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafCrteName"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafCrteDaytim"),
                      pilDateGetISO8601());
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafLchgName"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafLchgDaytim"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChckName"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChckDaytim"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChecksum"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafHeaderEnd"),NULL);


 /*
  * Date
  */

  desc = findDescriptor(descs, pilTrnGetKeyword("DateObs"));
  if (!desc) {
    cpl_free(pafName);
    return NULL;
  }

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvOptDate"), 
                      desc->descValue->s);
 /*
  * temp
  */

  desc = findDescriptor(descs, pilTrnGetKeyword("BeamTemperature", quad));
  if (desc == NULL) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvOptTemp"), 
                      desc->descValue->d);


 /*
  *  Write curvature model coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("CurvatureOrd"), 
                        &ord, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvModOrd"), ord);

  if (readIntDescriptor(descs, pilTrnGetKeyword("CurvatureOrdX"), 
                        &xOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvModXord"), xOrd);

  if (readIntDescriptor(descs, pilTrnGetKeyword("CurvatureOrdY"), 
                        &yOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvModYord"), yOrd);

  for (i = 0; i <= ord; i++) {
    for (j = 0; j <= xOrd; j++) {
      for (k = 0; k <= yOrd; k++) {

        if ((desc = findDescriptor(descs, 
            pilTrnGetKeyword("Curvature", i, j, k))) == NULL) {
          cpl_msg_error(modName, "Cannot read descriptor %s", 
                      pilTrnGetKeyword("Curvature", i, j, k));
          cpl_free(pafName);
          return NULL;
        }

        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvMod", i, j, k), 
                            desc->descValue->s);

      }
    }
  }


 /*
  *  Write optical distorsion model coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("OptDistOrdX"),
                        &xOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisXord"), xOrd);

  for (i = 0; i <= xOrd; i++) {
    for (j = 0; j <= xOrd; j++) {

      if ((desc = findDescriptor(descs,
          pilTrnGetKeyword("OptDistX", i, j))) == NULL) {
        cpl_msg_error(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("OptDistX", i, j));
        cpl_free(pafName);
        return NULL;
      }

      writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisX", i, j),
                          desc->descValue->s);

    }
  }

  if (readIntDescriptor(descs, pilTrnGetKeyword("OptDistOrdY"),
                        &yOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisYord"), yOrd);

  for (i = 0; i <= yOrd; i++) {
    for (j = 0; j <= yOrd; j++) {

      if ((desc = findDescriptor(descs,
          pilTrnGetKeyword("OptDistY", i, j))) == NULL) {
        cpl_msg_error(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("OptDistY", i, j));
        cpl_free(pafName);
        return NULL;
      }

      writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisY", i, j),
                          desc->descValue->s);

    }
  }

  fclose(fp);

  return pafName;

}


/**
 * @memo
 *   Create the PAF file containing the IDS coefficients for a given quadrant.
 *
 * @return Name of created PAF file, or NULL in case of failure.
 *
 * @param descs      List of descriptors containing the IDS.
 * @param namePAF    PAF file root name.
 *
 * @doc
 *   The PAF file name depends on the quadrant.
 *
 * @author C. Izzo
 */

char *
createIdsPAF(VimosDescriptor *descs, char *namePAF)
{
  const char       modName[] = "createIdsPAF";
  char            *pafName;
  int              ord, xOrd, yOrd;
  int              i, j, k, len, quad;
  FILE            *fp;
  VimosDescriptor *desc;

  cpl_msg_debug(modName, "Write IDS into PAF file");

  readIntDescriptor(descs, pilTrnGetKeyword("Quadrant"), &quad, NULL);

  len = sizeof(char) * (strlen(namePAF)+7);
  pafName = (char *)cpl_malloc(len);
  if (!pafName)
    return NULL;

  sprintf(pafName, "%s-%d.paf", namePAF, quad);
  fp = fopen(pafName, "w");
  if (!fp) {
    cpl_free(pafName);
    return NULL;
  }


 /*
  * Write PAF header
  */

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafHeaderStart"),NULL);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafType"),"Configuration");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafId"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafName"), pafName);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafDesc"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafCrteName"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafCrteDaytim"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafLchgName"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafLchgDaytim"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChckName"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChckDaytim"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChecksum"),"");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafHeaderEnd"),NULL);


 /*
  * Date
  */

  if ((desc = findDescriptor(descs, "DATE-OBS")) == NULL) {
    cpl_free(pafName);
    return NULL;
  }

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsDate"), 
                      desc->descValue->s);
 /*
  * temp
  */

  desc = findDescriptor(descs, pilTrnGetKeyword("BeamTemperature", quad));
  if (desc == NULL) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsTemp"), 
                      desc->descValue->d);


 /*
  *  Write IDS coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("DispersionOrd"), 
                        &ord, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsModOrd"), ord);

  if (readIntDescriptor(descs, pilTrnGetKeyword("DispersionOrdX"), 
                        &xOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsModXord"), xOrd);

  if (readIntDescriptor(descs, pilTrnGetKeyword("DispersionOrdY"), 
                        &yOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsModYord"), yOrd);

  for (i = 0; i <= ord; i++) {
    for (j = 0; j <= xOrd; j++) {
      for (k = 0; k <= yOrd; k++) {

        if ((desc = findDescriptor(descs, 
            pilTrnGetKeyword("Dispersion", i, j, k))) == NULL) {
          cpl_msg_error(modName, "Cannot read descriptor %s", 
                      pilTrnGetKeyword("Dispersion", i, j, k));
          cpl_free(pafName);
          return NULL;
        }

        writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsMod", i, j, k),
                            desc->descValue->d);

      }
    }
  }

  fclose(fp);

  return pafName;

}


/**
 * @memo
 *   Create the PAF file containing the coefficients for all spectral 
 *   distortions for a given quadrant.
 *
 * @return Name of created PAF file, or NULL in case of failure.
 *
 * @param descs      List of descriptors containing the distortion models.
 * @param namePAF    PAF file root name.
 *
 * @doc
 *   The PAF file name depends on the quadrant and the grism used.
 *
 * @author C. Izzo
 */

char *
createSpectralDistPAF(VimosDescriptor *descs, char *namePAF)
{
  const char       modName[] = "createIdsPAF";
  char            *pafName_noext;
  char             grismName[80];
  char             filename[PATHNAME_MAX + 1];
  char            *pafName;
  int              ord, xOrd, yOrd;
  int              i, j, k, len, quad;
  double           rms, pWidth;
  FILE            *fp;
  VimosDescriptor *desc;

  cpl_msg_debug(modName, "Write spectral distortion models into PAF file");

  readIntDescriptor(descs, pilTrnGetKeyword("Quadrant"), &quad, NULL);

  readStringDescriptor(descs, pilTrnGetKeyword("GrismName", quad),
                       grismName, NULL);
  sprintf(filename, "%s_%s_%d.cmf", namePAF, grismName, quad);

  len = strlen(filename) + 1;
  pafName = cpl_malloc(len * sizeof(char));
  strcpy(pafName, filename);    /* Need this, because pafName is returned */

  fp = fopen(pafName, "w");
  if (!fp) {
    cpl_free(pafName);
    return NULL;
  }


 /*
  * Write PAF header
  */

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafHeaderStart"), NULL);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafType"),
                      "Configuration");
  pafName_noext = (char *)cpl_malloc(len - 4);
  sprintf(pafName_noext, "%s_%s_%d", namePAF, grismName, quad);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafId"), pafName_noext);
  cpl_free(pafName_noext);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafName"), pafName);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafDesc"), "");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafCrteName"), "");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafCrteDaytim"), 
                      pilDateGetISO8601());
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafLchgName"), "");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafLchgDaytim"), "");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChckName"), "");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChckDaytim"), "");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafChecksum"), "");
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafHeaderEnd"), NULL);


 /*
  * IDS date
  */

  desc = findDescriptor(descs, pilTrnGetKeyword("DateObs"));
  if (!desc) {
    cpl_free(pafName);
    return NULL;
  }

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsDate"), 
                      desc->descValue->s);

 /*
  * IDS temperature
  */

  desc = findDescriptor(descs, pilTrnGetKeyword("BeamTemperature", quad));
  if (!desc) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsTemp"), 
                      desc->descValue->d);


 /*
  *  IDS coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("DispersionOrd"), 
                        &ord, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsModOrd"), ord);

  if (readIntDescriptor(descs, pilTrnGetKeyword("DispersionOrdX"), 
                        &xOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsModXord"), xOrd);

  if (readIntDescriptor(descs, pilTrnGetKeyword("DispersionOrdY"), 
                        &yOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsModYord"), yOrd);

  for (i = 0; i <= ord; i++) {
    for (j = 0; j <= xOrd; j++) {
      for (k = 0; k <= yOrd; k++) {

        if ((desc = findDescriptor(descs, 
            pilTrnGetKeyword("Dispersion", i, j, k))) == NULL) {
          cpl_msg_error(modName, "Cannot read descriptor %s", 
                      pilTrnGetKeyword("Dispersion", i, j, k));
          cpl_free(pafName);
          return NULL;
        }

        writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsMod", i, j, k),
                            desc->descValue->d);

      }
    }
  }

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("IdsXrms"),
                           &rms, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsXrms"), rms);

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("IdsYrms"),
                           &rms, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFIdsYrms"), rms);


 /*
  * Date of optical and curvature models
  */

  desc = findDescriptor(descs, pilTrnGetKeyword("DateObs"));
  if (!desc) {
    cpl_free(pafName);
    return NULL;
  }

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDate"), 
                      desc->descValue->s);


 /*
  * Optical distortion temperature (valid also for curvature?)
  */

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("OptTemp", quad))) == 0) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptTemp"), 
                      desc->descValue->d);


 /*
  *  Write optical distorsion model coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("OptDistOrdX"),
                        &xOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisXord"), xOrd);

  for (i = 0; i <= xOrd; i++) {
    for (j = 0; j <= xOrd; j++) {

      if ((desc = findDescriptor(descs,
          pilTrnGetKeyword("OptDistX", i, j))) == NULL) {
        cpl_msg_error(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("OptDistX", i, j));
        cpl_free(pafName);
        return NULL;
      }

      writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisX", i, j),
                          desc->descValue->d);

    }
  }

  if (readIntDescriptor(descs, pilTrnGetKeyword("OptDistOrdY"),
                        &yOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisYord"), yOrd);

  for (i = 0; i <= yOrd; i++) {
    for (j = 0; j <= yOrd; j++) {

      if ((desc = findDescriptor(descs,
          pilTrnGetKeyword("OptDistY", i, j))) == NULL) {
        cpl_msg_error(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("OptDistY", i, j));
        cpl_free(pafName);
        return NULL;
      }

      writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisY", i, j),
                          desc->descValue->d);

    }
  }

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("OptXrms"),
                           &rms, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisXrms"), rms);

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("OptYrms"),
                           &rms, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFOptDisYrms"), rms);


 /*
  *  Write curvature model coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("CurvatureOrd"), 
                        &ord, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvModOrd"), ord);

  if (readIntDescriptor(descs, pilTrnGetKeyword("CurvatureOrdX"), 
                        &xOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvModXord"), xOrd);

  if (readIntDescriptor(descs, pilTrnGetKeyword("CurvatureOrdY"), 
                        &yOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvModYord"), yOrd);

  for (i = 0; i <= ord; i++) {
    for (j = 0; j <= xOrd; j++) {
      for (k = 0; k <= yOrd; k++) {

        if ((desc = findDescriptor(descs, 
            pilTrnGetKeyword("Curvature", i, j, k))) == NULL) {
          cpl_msg_error(modName, "Cannot read descriptor %s", 
                      pilTrnGetKeyword("Curvature", i, j, k));
          cpl_free(pafName);
          return NULL;
        }

        writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvMod", i, j, k), 
                            desc->descValue->d);

      }
    }
  }

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("CurvatureXrms"),
                           &rms, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvXrms"), rms);

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("CurvatureYrms"),
                           &rms, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCrvYrms"), rms);


 /*
  *  Write contamination model coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("ZeroOrdX"),
                        &xOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFZeroXord"), xOrd);

  for (i = 0; i <= xOrd; i++) {
    for (j = 0; j <= xOrd; j++) {

      if ((desc = findDescriptor(descs,
          pilTrnGetKeyword("ZeroX", i, j))) == NULL) {
        cpl_msg_error(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("ZeroX", i, j));
        cpl_free(pafName);
        return NULL;
      }

      writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFZeroX", i, j),
                          desc->descValue->d);

    }
  }

  if (readIntDescriptor(descs, pilTrnGetKeyword("ZeroOrdY"),
                        &yOrd, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFZeroYord"), yOrd);

  for (i = 0; i <= yOrd; i++) {
    for (j = 0; j <= yOrd; j++) {

      if ((desc = findDescriptor(descs,
          pilTrnGetKeyword("ZeroY", i, j))) == NULL) {
        cpl_msg_error(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("ZeroY", i, j));
        cpl_free(pafName);
        return NULL;
      }

      writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFZeroY", i, j),
                          desc->descValue->d);

    }
  }

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("ZeroXrms"),
                           &rms, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFZeroXrms"), rms);

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("ZeroYrms"),
                           &rms, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFZeroYrms"), rms);

  if (readDoubleDescriptor(descs, pilTrnGetKeyword("ZeroOrderWidth"),
                           &pWidth, NULL) == VM_FALSE) {
    cpl_free(pafName);
    return NULL;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFZeroWidth"), pWidth);



  fclose(fp);

  return pafName;

}
