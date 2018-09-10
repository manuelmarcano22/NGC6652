/* $Id: vmimgresampling.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <assert.h>

#include <vimoswcs.h>

#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmimage.h"
#include "vmimagearray.h"
#include "vmwcsutils.h"
#include "vmimgresampling.h"


#define BI_LINEAR  1
#define BI_CUBIC   2


/*
 * @brief
 *   Swap two double values.
 *
 * @return Nothing.
 * 
 * @param a  Pointer to first value
 * @param b  Pointer to second value
 * 
 * The function exchanges the double values stored in @em a and @em b
 * in place.
 *
 * @author P. Montegriffo
 */

/*
 * This one is no longer necessary, after fixing the getFieldOfView()
 *

static void swap(double *a, double *b)
{
  double dmp = *a;

  *a = *b;
  *b = dmp;

   return;
}

*********/
 

/*
 * @brief
 *   Compute an interpolated value.
 *
 * @return The interpolated value.
 * 
 * @param image  Image pixel buffer.
 * @param nx     Number of image pixels along x.
 * @param ny     Number of image pixels along y.
 * @param ifunc  Selection flag for the interpolation function.
 * @param x      x coordinate of point to interpolate
 * @param y      y coordinate of point to interpolate
 * 
 * The function omputes an interpolated value at a position given by @em x
 * and @em y using a bi-linear or a bi-cubic interpolation function. The
 * interpolation function is selected by the flag @em ifunc. Passing 2 for
 * @em ifunc selects the bi-cubic interpolation. Any other value selects
 * the bi-linear interpolation function.
 * 
 * @note
 *   Coefficents and formulae for the bi-cubic interpolation are taken from
 *   'Numerical Recipes', II Ed., p. 119f.
 *
 * @author P. Montegriffo
 */

static double interpolate(float *image, int nx, int ny, int ifunc,
                          double x, double y)
{

  int  i, j, ii, kk;
  int  ind[4];

  double val;
  double v[16], t, u;
  double xx, c[16];

  static int wt[16][16] = {
    {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
    {-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
    {2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
    {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
    {0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
    {0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
    {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
    {9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
    {-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
    {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
    {-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
    {4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}
  };


  val = 0.0;                /*** Value for undefined pixels ***/

  i = (int) floor (x);
  j = (int) floor (y);

  if (i <= 0 || j <= 0 || i >= nx - 1 || j >= ny - 1 )
    return (0.0);

  /* Indexes for 4 closest neighbors positions */

  ind[0] = nx * j + i;
  ind[1] = nx * j + i + 1;
  ind[2] = nx * (j + 1) + i + 1;
  ind[3] = nx * (j + 1) + i;


  /* Values for 4 closest neighbors pixels */

  for (ii = 0; ii < 4; ii++) {
    v[ii] = *(image + ind[ii]);

    /* FIXME:
     *   What's the use of this print statment, especially since the 
     *   hardcoded values for x and y do not correspond to the size
     *   of a VIMOS image. Did I miss something. (RP)
     *
     *   Temporarily disabled, since not present in the original version. (RP)
     */
    /*
      if ((y >= 2015.) && (y <= 2016.) && (x >= 2041.0)) {
        printf("ii= %d, ind[ii]=%d, v[ii]=%f\n", ii, ind[ii], 
               *(image + ind[ii]));
      }
    */
  }

  t = x - (double) i;
  u = y - (double) j;

  switch (ifunc) {
      case BI_CUBIC:
        if ( i <= 1 || j <= 1 || i >= nx - 2 || j >= ny - 2 )
          return (0.0);

        /* Compute the gradients along x,y and the cross derivative */
     
        for (ii = 0; ii < 4; ii++) {
          v[ii + 4] = (*(image + ind[ii] + 1) -
                       *(image +ind[ii] - 1)) / 2.0;

          v[ii + 8] = (*(image + ind[ii] + nx) -
                       *(image + ind[ii] - nx)) / 2.0;

          v[ii + 12] = (*(image + ind[ii] + nx + 1) -
                        *(image + ind[ii] - nx + 1) - 
                        *(image+ind[ii]+nx-1) +
                        *(image+ind[ii]-nx-1)) / 4.0;
        }
     
        /* Get coefficents for bicubic interpolation */
     
        for (ii = 0; ii < 16; ii++) {
          xx = 0.0;
          for (kk=0; kk<16; kk++)
            xx += wt[ii][kk]*v[kk];
          c[ii] = xx;
        }         

        /* Compute interpolated value */
     
        for (ii = 3; ii >= 0; ii--)
          val = t * val + ((c[ii * 4 + 3] * u + c[ii * 4 + 2]) * u +
                           c[ii * 4 + 1]) * u + c[ii * 4];
     
        break;

      default:
        val = (double)( (1.0 - t) * (1.0 - u) * v[0] 
                       +       t  * (1.0 - u) * v[1] 
                       +       t  *        u  * v[2] 
                       +(1.0 - t) *        u  * v[3] );
        
        /* FIXME:
         *   What's the use of this print statment, especially since the 
         *   hardcoded values for x and y do not correspond to the size
         *   of a VIMOS image. Did I miss something. (RP)
         *
         *   Temporarily disabled, since not present in the original
         *   version. (RP)
         */
        /*
          if (y >= 2015. && y <= 2016. && x >= 2041.0) {
            printf("bi_linear: i=%d, j=%d, t=%f, u=%f, val=%f, v[0]=%f, "
                   "v[1]=%f, v[2]=%f, v[3]=%f\n", i, j, t, u, val, v[0],
                   v[1], v[2], v[3]);
           }
        */
        break;
  }

  return (val);

}


/**
 * @brief
 *   Determine the total field of view of a set of images.
 *
 * @return The function returns a pointer to the WCS structure of the
 *   reference image if no error occurred, otherwise @c NULL is returned.
 *
 * @param images  Set of images.
 * @param ra1     Start right ascension of the field of view.
 * @param dec1    Start declination  of the field of view.
 * @param ra2     End right ascension of the field of view.
 * @param dec2    End declination  of the field of view.
 *
 * The function determines the union of the sky regions covered by the
 * images in the input array @em images. The first image in the input
 * array is, arbitrarily choosen, as reference, and is returned if
 * the field of view determination was successful.
 *
 * @author P. Montegriffo, R. Palsa
 */

/* Remove, because of wrongly computed FOV: the pixel having
 * simultaneously the max observed RA and Dec is not necessarily
 * at the corner of the union of all images (think of the case
 * when the field is rotated); the same is true for the pixel
 * having simultaneously the min observed RA and Dec.
 * 

static struct WorldCoor *getFieldOfView(VimosImageArray *images,
                                        double *ra1, double *dec1,
                                        double *ra2, double *dec2)
{

  int i;

  double wra1, wra2, wdec1, wdec2;

  struct  WorldCoor *wcs_r = 0L;
  struct  WorldCoor *wcs_w = 0L;

  VimosImage *image = (VimosImage *)imageArrayGet(images, 0);


   *
   * Read reference WCS from the first (arbitrarily choosen) image
   * in the list.
   * 

  if (!(wcs_r = rdimage(image->descs)))
    return 0L;

  vimoswcsrange(wcs_r, ra1, ra2, dec1, dec2);

  for (i = 1; i < imageArraySize(images); i++) {
    image = (VimosImage *)imageArrayGet(images, i);

    if (!(wcs_w = rdimage(image->descs))) {
      vimoswcsfree(wcs_r);
      return 0L;
    }

    vimoswcsrange (wcs_w, &wra1, &wra2, &wdec1, &wdec2);

    *ra1 = wra1 < *ra1 ? wra1 : *ra1;
    *ra2 = wra2 > *ra2 ? wra2 : *ra2;

    *dec1 = wdec1 < *dec1 ? wdec1 : *dec1;
    *dec2 = wdec2 > *dec2 ? wdec2 : *dec2;

    vimoswcsfree(wcs_w);
  }

  return wcs_r;

}

**********/

/*
 * This is the new version to compute the FOV. That is, this is
 * the one that really finds it.
 */

static struct WorldCoor *getFieldOfViewNew(VimosImageArray *images,
                                        double *xpix1, double *ypix1,
                                        double *xpix2, double *ypix2)
{

  int offscl;
  int i;

  double ra1, ra2, dec1, dec2;
  double xmin, ymin, xmax, ymax;
  double x, y;

  struct  WorldCoor *wcs_r = 0L;
  struct  WorldCoor *wcs_w = 0L;


  VimosImage *image = (VimosImage *)imageArrayGet(images, 0);


  /*
   * Read reference WCS from the first (arbitrarily choosen) image
   * in the list.
   */

  if (!(wcs_r = rdimage(image->descs)))
    return 0L;

  xmin = 1.0;
  ymin = 1.0;
  xmax = wcs_r->nxpix;
  ymax = wcs_r->nypix;

  for (i = 1; i < imageArraySize(images); i++) {
    image = (VimosImage *)imageArrayGet(images, i);

    if (!(wcs_w = rdimage(image->descs))) {
      vimoswcsfree(wcs_r);
      return 0L;
    }

    pix2vimoswcs(wcs_w, 1., 1., &ra1, &dec1);
    pix2vimoswcs(wcs_w, wcs_w->nxpix, wcs_w->nypix, &ra2, &dec2);

    vimoswcs2pix(wcs_r, ra1, dec1, &x, &y, &offscl);
    if (x < xmin)
      xmin = x;
    if (y < ymin)
      ymin = y;

    vimoswcs2pix(wcs_r, ra2, dec2, &x, &y, &offscl);
    if (x > xmax)
      xmax = x;
    if (y > ymax)
      ymax = y;

    vimoswcsfree(wcs_w);
  }

  *xpix1 = xmin;
  *ypix1 = ymin;
  *xpix2 = xmax;
  *ypix2 = ymax;

  return wcs_r;

}


/**
 * @brief
 *   Setup a new world coordinate system.
 *
 * @return The function returns a reference to the newly created WCS
 *   structure if no error occurred, otherwise the return value is
 *   @c NULL.
 *
 * @param wcs   Reference world coordinate system.
 * @param ra1   Start right ascension of the field of view.
 * @param dec1  Start declination  of the field of view.
 * @param ra2   End right ascension of the field of view.
 * @param dec2  End declination  of the field of view.
 *
 * The functions computes the pixel position of corners of the field of
 * view given by @em ra1, @em dec1 and @em ra2, @em dec2 with respect
 * to the coordinate system defined by @em wcs. The pixel coordinates
 * of the reference pixel defined by @em wcs are adjusted accordingly.
 *
 * @author P. Montegriffo
 */

/*
 * This function is outdated, after fixing the getFieldOfView().
 *

static struct WorldCoor *setupWcs(struct WorldCoor *wcs, double ra1,
                                  double dec1, double ra2, double dec2)
{

   int offscl;

   double xpix1, ypix1, xpix2, ypix2;
   double crpix1, crpix2;

   struct WorldCoor *wcs_r = 0L;
  

    * Initialize returning structure * 

    * FIXME:
    *  The original code sets the CD matrix to NULL, so that wcskinit()
    *  ignores it and uses the values from CDELTi and CROTAi instead
    *  to setup the WCS structure. But note that the original code
    *  used a hardcoded 0.0 for the rotation angle.
    *
    *  Check if this was done on purpose. There was no hint about that
    *  in the original code. (RP)
    *
    *  Original code follows:
    *
    *  double * cd = NULL;
    *  double   crota = 0.0;
    *
    *  .
    *  .
    *  .
    *
    *  wcs_r = wcskinit((int)wcs->nxpix, (int)wcs->nypix,
    *                   "RA---TAN", "DEC--TAN",
    *                   wcs->crpix[0], wcs->crpix[1],
    *                   wcs->crval[0], wcs->crval[1], cd,
    *                   wcs->cdelt[0], wcs->cdelt[1], crota,
    *                   (int)wcs->equinox, wcs->epoch);
    *
    * 

   wcs_r = wcskinit((int)wcs->nxpix, (int)wcs->nypix,
                    "RA---TAN", "DEC--TAN",
                    wcs->crpix[0], wcs->crpix[1],
                    wcs->crval[0], wcs->crval[1], wcs->cd,
                    0., 0., 0., (int)wcs->equinox, wcs->epoch);

   if (!wcs_r)
     return 0L;


    * Compute corner pixels coordinates in the new framework * 

   wcs2pix(wcs_r, ra1, dec1, &xpix1, &ypix1, &offscl);
   wcs2pix(wcs_r, ra2, dec2, &xpix2, &ypix2, &offscl);

   if (xpix1 > xpix2)
     swap(&xpix1, &xpix2);

   if (ypix1 > ypix2)
     swap(&ypix1, &ypix2);


    * Set new image dimensions and reference pixel * 

    * FIXME:
    *   The original code had the following 2 line just before
    *   the 2 floor() calls below.
    *
    *   wcs_r->nxpix = nint(xpix2) - nint(xpix1) + 1.0;
    *   wcs_r->nypix = nint(ypix2) - nint(ypix1) + 1.0;
    *
    *   They had no effect since their result was overwritten by the
    *   two calls to floor(). Anyway the function nint() is not available
    *   on every system, it belongs to the sunmath library available
    *   on SUNs. So if this alternative should be used it can be emulated by
    *
    *   wcs_r->nxpix = (int)(xpix2 + 0.5) - (int)(xpix1 + 0.5) + 1.0;
    *   wcs_r->nypix = (int)(ypix2 + 0.5) - (int)(ypix1 + 0.5) + 1.0;
    * 

   wcs_r->nxpix = floor(xpix2 - xpix1 + 1.0);
   wcs_r->nypix = floor(ypix2 - ypix1 + 1.0);

   crpix1 = wcs_r->crpix[0] - floor(xpix1);
   crpix2 = wcs_r->crpix[1] - floor(ypix1);

   if (wcsreset(wcs_r, crpix1, crpix2, wcs_r->crval[0], wcs_r->crval[1],
                wcs_r->cdelt[0], wcs_r->cdelt[1], 0., wcs_r->cd)) {
     wcsfree(wcs_r);
     return 0L;
   }

   return wcs_r;

}

****************/

/*
 * This is the new function replacing the above one - note how
 * simple it is...
 */

static struct WorldCoor *setupWcsNew(struct WorldCoor *wcs, double xpix1,
                                  double ypix1, double xpix2, double ypix2)
{

   double crpix1, crpix2;

   wcs->nxpix = floor(xpix2 - xpix1 + 1.0);
   wcs->nypix = floor(ypix2 - ypix1 + 1.0);

   crpix1 = wcs->crpix[0] - floor(xpix1);
   crpix2 = wcs->crpix[1] - floor(ypix1);

   if (vimoswcsreset(wcs, crpix1, crpix2, wcs->crval[0], wcs->crval[1],
                wcs->cdelt[0], wcs->cdelt[1], 0., wcs->cd)) {
     return 0L;
   }

   return wcs;

}

/**
 * @brief
 *   Update an image keyword list with a new world coordinate system.
 *
 * @return The function returns 0 if the keyword list was updated
 *   successfully, otherwise 1 is returned.
 *
 * @param image  Image object.
 * @param wcs    WCS information.
 *
 * The function updates the world coordinate system related keywords
 * @li CRPIX1,
 * @li CRPIX2,
 * @li CRVAL1,
 * @li CRVAL2,
 * and the CD matrix coefficients of the image @em image, i.e. it's keyword
 * list, with the new values taken from @em wcs. If the keyword do not
 * exist they are created.
 *
 * @author P. Sartoretti
 */

static int updateWcsKeywords(VimosImage *image, struct WorldCoor *wcs)
{

  register int  i, j, k;


  if (writeDoubleDescriptor(&image->descs,
                            pilTrnGetKeyword("Crpix", 1), wcs->xrefpix,
                            pilTrnGetComment("Crpix")) == VM_FALSE) {
    return 1;
  }

  if (writeDoubleDescriptor(&image->descs,
                            pilTrnGetKeyword("Crpix", 2), wcs->yrefpix,
                            pilTrnGetComment("Crpix")) == VM_FALSE) {
    return 1;
  }

  if (writeDoubleDescriptor(&image->descs,
                            pilTrnGetKeyword("Crval", 1), wcs->xref,
                            pilTrnGetComment("Crval")) == VM_FALSE) {
    return 1;
  }

  if (writeDoubleDescriptor(&image->descs,
                            pilTrnGetKeyword("Crval", 2), wcs->yref,
                            pilTrnGetComment("Crval")) == VM_FALSE) {
    return 1;
  }


  k = 0;

  for(i = 1; i <= 2; i++) {
    for(j = 1; j <= 2; j++) {
      if (writeDoubleDescriptor(&image->descs,
                                pilTrnGetKeyword("CD", i, j), wcs->cd[k],
                                pilTrnGetComment("CD")) == VM_FALSE) {
	return 1;
      }

      k++;
    }
  }
      
  return 1;

}


/**
 * @brief
 *   Resample a set of images to a common coordinate grid.
 *
 * @return The function returns TBD.
 *
 * @param set     Set of images.
 * @param method  Interpolation method.
 *
 * TBD
 *
 * @author P. Montegriffo, R. Palsa
 */

VimosImageArray *VmImResampleImages(VimosImageArray *set,
                                    ResamplingMethod method)
{

  const char *fctid = "VmImResampleImages";


  int i;
  int nx, ny;

  size_t k;
  size_t npixel;

/*  double ra1, ra2, dec1, dec2;  */

  double xpix1, ypix1, xpix2, ypix2;

  VimosImage *srcImage, *dstImage;

  struct WorldCoor *wcs_r, *wcs;


  /*
   * Get the full field of view, i.e. the ranges of right ascension
   * and declination, covered by the union of the fields of view of
   * all images in the input set.
   */

  cpl_msg_info(fctid, "Computing sky area covered by image set ...");

/*
 * The original code was incorrectly defining the new image corners
 * from the ranges in RA and Dec... Instead, the new image corners
 * are simply the min-max range of all image corners, computed in
 * the WCS of the reference frame!
 *
 *  if (!(wcs_r = getFieldOfView(set, &ra1, &dec1, &ra2, &dec2))) {
 *    cpl_msg_error(fctid, "Sky area computation failed!");
 *    return 0L;
 *  }
 *
 *  cpl_msg_info(fctid, "R.A. covered range: %10.6f --> %10.6f", ra1, ra2);
 *  cpl_msg_info(fctid, "Dec. covered range: %10.6f --> %10.6f", dec1, dec2);
 */

  if (!(wcs_r = getFieldOfViewNew(set, &xpix1, &ypix1, &xpix2, &ypix2))) {
    cpl_msg_error(fctid, "Sky area computation failed!");
    return 0L;
  }

  /*
   * Setup new WCS for the full field of view
   */

/*
 * See note above about incorrectly computed FOV.
 *
 *  if (!(wcs = setupWcs(wcs_r, ra1, dec1, ra2, dec2))) {
 *    cpl_msg_error(fctid, "Field of view coordinate system creation failed!");
 *
 *    vimoswcsfree(wcs_r);
 *    return 0L;
 *  }
 */

  if (!(wcs = setupWcsNew(wcs_r, xpix1, ypix1, xpix2, ypix2))) {
    cpl_msg_error(fctid, "Field of view coordinate system creation failed!");

    vimoswcsfree(wcs_r);
    return 0L;
  }
  else {
    cpl_msg_info(fctid, "New image dimensions: NAXIS1=%d, NAXIS2=%d",
                (int)wcs->nxpix, (int)wcs->nypix);
    cpl_msg_debug(fctid, "New reference pixel:  CRPIX1=%10.3f, CRPIX2=%10.3f",
                wcs->crpix[0], wcs->crpix[1]);
  }

  /*
   * No need to free memory, it's the input wcs_r the it's output to
   * from setupWcsNew(), i.e., wcs_r = wcs (same pointer).
   * 
   *
   *  vimoswcsfree(wcs_r);
   *
   *****/


  /*
   * Remap imput images to the new coordinate grid
   */

  nx = (int)wcs->nxpix;
  ny = (int)wcs->nypix;
  
  npixel = nx * ny;


  for (i = 0; i < imageArraySize(set); i++) {

    cpl_msg_info(fctid, "Resampling image %d.", i);


    /*
     * Get source image from the input set and retrieve its world
     * coordinate system information.
     */

    srcImage = (VimosImage *)imageArrayGet(set, i);

    if (!(wcs_r = rdimage(srcImage->descs))) {
      cpl_msg_error(fctid, "Retrieving coordinate system information from "
                  "image %d failed!", i);

      vimoswcsfree(wcs);
      return 0L;
    }


    if (!(dstImage = newImageAndAlloc(nx, ny))) {
      cpl_msg_error(fctid, "Not enough memory!");

      vimoswcsfree(wcs_r);
      vimoswcsfree(wcs);

      return 0L;
    }


    /*
     * Copy the image header to the resampled image and update 
     * the WCS related keywords.
     */

    if (copyAllDescriptors(srcImage->descs, &dstImage->descs) == VM_FALSE) {
      cpl_msg_error(fctid, "Copying header keywords failed!");
      
      deleteImage(dstImage);

      vimoswcsfree(wcs_r);
      vimoswcsfree(wcs);

      return 0L;
    }

    if (!updateWcsKeywords(dstImage, wcs)) {
      cpl_msg_error(fctid, "Updating world coordinate system information "
                  "failed!");
      deleteImage(dstImage);

      vimoswcsfree(wcs_r);
      vimoswcsfree(wcs);

      return 0L;
    }


    /*
     * Loop over pixels. For each pixel position of the resampled image get
     * its position on the sky and use this position to compute its pixel
     * value, by interpolation, from the input image pixels.
     */

    /* FIXME:
     *  For performance reasons this loop over all pixels of the result
     *  image should be done inside the resampling function. Calling
     *  interpolate() for each pixel is not very efficient. (RP)
     */

    for (k = 0; k < npixel; k++) {

      int offscl;

      double xpix, ypix, x, y;
      double ra, dec;
      

      /*
       * Compute world coordinates for pixel center
       */

      xpix = (double)(k % nx) + 1.0;
      ypix = (double)(k / nx) + 1.0;

      pix2vimoswcs(wcs, xpix, ypix, &ra, &dec);

      /*
       * Compute pixel coordinates of the source image corresponding to
       * the just computed ra and dec.
       */

      vimoswcs2pix(wcs_r, ra, dec, &x, &y, &offscl);

      x = x - 1.0 ;
      y = y - 1.0 ;

/*
if (k % 64 == 0) {
printf("XXXX %d %d %f %f\n", (int)xpix - 1, (int)ypix - 1, x, y);
}
*/

      /*
       * Compute interpolated value
       */

      /* FIXME:
       *   Instead of DOUBLENULLVALUE (which comes from cfitsio) a simple
       *   0. does the same job. To be checked! (RP)
       */

      if (!offscl)
        dstImage->data[k] = interpolate(srcImage->data, wcs_r->nxpix,
                                        wcs_r->nypix, method, x, y);
      else
        dstImage->data[k] = DOUBLENULLVALUE;

    }


    /*
     * Replace the source image in the input set with the resampled
     * result. The source is deallocated.
     */

    deleteImage((VimosImage *)imageArrayRemove(set, i));
    imageArraySet(set, i, dstImage);

    vimoswcsfree(wcs_r);
  }


  /*
   * Cleanup
   */

  vimoswcsfree(wcs);
    
  return set;

}
