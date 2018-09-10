/* $Id: vmwcsutils.c,v 1.4 2013-03-25 11:43:04 cgarcia Exp $
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

#include <pilmacros.h>
#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmmath.h"
#include "vmimgutils.h"
#include "vmwcsutils.h"
#include "cpl.h"


/*-----------------------------------------------------------------------------
 * Function  : wcstopix
 * In        : Number of elements, a_star structure, wcs structure
 * Out       : Projected a_star->ximage & a_star->yimage coordinates
 * Purpouse  : Project sky coordinates on pixel coordinates 
 * Note      : Uses WCS Tool's 'wcs2pix' subroutine
 *---------------------------------------------------------------------------*/

void
wcstopix(int n_ele, VimosTable * a_star, struct WorldCoor * wcs)
{
   int i;
   VimosColumn *raCol, *decCol, *xpixCol, *ypixCol;
   VimosColumn  *goffCol, *magCol; 
   char modName[] = "wcstopix";

   if((raCol = findColInTab(a_star, "RA")) == NULL) {
     cpl_msg_error(modName, "Column RA not found in Astrometric table");
     return;
   } 
   
  if(( decCol = findColInTab(a_star, "DEC")) == NULL) {
   cpl_msg_error(modName, "Column DEC not found in Astrometric table");
     return;
  } 

  if((magCol = findColInTab(a_star, "MAG")) == NULL)
    cpl_msg_warning(modName, "Column MAG not found in Astrometric Table");

   /****** The first time the function is called the astrometric table 
           only contains the column ID, RA, DEC, and MAG and does not contain 
           the columns X_IMAGE, Y_IMAGE GOFF so, if they are not yet
	   there  
	   create them *******/


   if((xpixCol = findColInTab(a_star, "X_IMAGE"))==NULL) {
     xpixCol = newDoubleColumn(n_ele, "X_IMAGE");
     tblAppendColumn(a_star, xpixCol);
   }
   if((ypixCol = findColInTab(a_star, "Y_IMAGE")) == NULL) {
     ypixCol = newDoubleColumn(n_ele, "Y_IMAGE");   
     tblAppendColumn(a_star, ypixCol);
   }
   if((goffCol = findColInTab(a_star, "GOFF")) == NULL) {
     goffCol = newIntColumn(n_ele, "GOFF");
     tblAppendColumn(a_star, goffCol);
   }
     

   for(i=0; i<n_ele; i++) {
     vimoswcs2pix (wcs, raCol->colValue->dArray[i], decCol->colValue->dArray[i],
              &(xpixCol->colValue->dArray[i]), 
	      &(ypixCol->colValue->dArray[i]), 
	      &(goffCol->colValue->iArray[i]));
   }

   return ;
}

/*-----------------------------------------------------------------------------
 * Function  : pixtowcs
 * In        : Number of elements, a_star structure, wcs structure
 * Out       : Projected a_star->ximage & a_star->yimage coordinates
 * Purpouse  : Project pixel coordinates on sky coordinates
 * Note      : Use WCS Tool's 'pix2wcs' subroutine
 *---------------------------------------------------------------------------*/

void
pixtowcs(int n_ele, VimosTable * o_star, struct WorldCoor * wcs)
{
   int i;
   VimosColumn *xpixCol, *ypixCol, *xwldCol, *ywldCol;

   xpixCol = findColInTab(o_star, "X_IMAGE");
   ypixCol = findColInTab(o_star, "Y_IMAGE");
   xwldCol = findColInTab(o_star, "X_WORLD");
   ywldCol = findColInTab(o_star, "Y_WORLD");

   for(i=0; i<n_ele; i++) {
     xwldCol->colValue->dArray[i] = 0.0;
     ywldCol->colValue->dArray[i] = 0.0;
     pix2vimoswcs (wcs, xpixCol->colValue->dArray[i], ypixCol->colValue->dArray[i],
              &(xwldCol->colValue->dArray[i]), 
	      &(ywldCol->colValue->dArray[i]));
   }
   return ;
}

/**
 * @memo 
 *   Fit the matching observed and catalog stars to find new plate center. 
 *   
 * @return VM_TRUE/VM_FALSE. The input wcs is updated
 * 
 * @param wcs         wcs structure to update 
 * @param o_star      star match table
 * @param nmatch      number of matching stars
 * @doc 
 *   Fit the matched image and catalog stars to find the new plate center.
 *   Uses FitMatch from WCSTools (minimization with NR's AMOEBA). 
 *
 * @author P. Montegriffo (i/o modified by P. Sartoretti)
 *         modified by Bianca Garilli
 */

/*-----------------------------------------------------------------------------
 * Function  : VmFitMatch
 * In        : wcs structure, o_star structure, a_star structure, index table 
 *             for matching stars and number of matches
 * Out       : updated wcs structure
 * Purpouse  : fit the matched image and catalog stars 
 *             to find the new plate center.
 * Note      : Uses FitMatch from WCSTools (minimization with NR's AMOEBA)
 *---------------------------------------------------------------------------*/

VimosBool vimosFitMatch(struct WorldCoor *wcs,VimosTable *o_star, int nmatch)
{
  int      i,j,nmatch1,iterate,nfit;
  double * in_ximage=0;
  double * in_yimage=0;
  double * in_xworld=0;
  double * in_yworld=0;
  double *xe = 0;
  double *ye = 0;
  double *resid = 0;
  double tmp;
  double mx, my, xsig, ysig, rsig, siglim;
  double xsum = 0.0;
  double ysum = 0.0;
  double rsum = 0.0;
  double xmean, ymean, rmean, xsumsq, ysumsq, diff;
  double dmatch ,dmatch1 ;

  VimosColumn      *inXwldCol, *inYwldCol;
  VimosColumn      *inXimaCol, *inYimaCol;
  char   modName[] = "vimosFitMatch";
  double rmsTol = 0.05;   /* Max residuals rms in arcsec (should this be a parameter?)*/

  /* Allocate memory for working arrays */
  
   if(!(in_ximage = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ||
      !(in_yimage = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ||
      !(in_xworld = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ||
      !(in_yworld = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ){
     cpl_msg_error(modName, "Could not alloc memory for finding plate solution");
     return VM_FALSE;
   }
   
   /* Set image and reference matching stars coordinates */
   if(!(inXimaCol = findColInTab(o_star, "X_IMAGE"))) {
     cpl_msg_error(modName, "Astrometric Table: Column with X-pixel coord "
		 "not found");
     return VM_FALSE;
   }
   if(!(inYimaCol = findColInTab(o_star, "Y_IMAGE"))) {
     cpl_msg_error(modName, "Astrometric Table: Column with Y-pixel coord "
		 "not found");
     return VM_FALSE;
   }
   if(!(inXwldCol = findColInTab(o_star, "X_WORLD"))) {
     cpl_msg_error(modName, "Astrometric Table: Column with RA coord "
		 "not found");
     return VM_FALSE;
   }
   if(!(inYwldCol = findColInTab(o_star, "Y_WORLD"))) {
     cpl_msg_error(modName, "Astrometric Table: Column with RA coord "
		 "not found");
     return VM_FALSE;
   }

   
   
   for(i=0; i<nmatch; i++) {
     in_ximage[i] = inXimaCol->colValue->dArray[i];
     in_yimage[i] = inYimaCol->colValue->dArray[i];
     in_xworld[i] = inXwldCol->colValue->dArray[i]; 
     in_yworld[i] = inYwldCol->colValue->dArray[i]; 
   }
   

   /*
    * Fit matching stars
    * ---> Fit only x,y shift and rotation 
    *
    * do not refine over esiduals: we want more control over 
    * this step, thus we do it in here instead of in wcslib
    */

   setnfit(-125);
   nfit = 3;

   resid = (double *)cpl_malloc(nmatch * sizeof(double));
   xe = (double *)cpl_malloc(nmatch * sizeof(double));
   ye = (double *)cpl_malloc(nmatch * sizeof(double));
   dmatch = (double)nmatch;
   dmatch1 = (double)(nmatch - 1);

   iterate = 0;

   /* do 1 iteration only (more is useless) */
   while (iterate<3) {
     /*
       FitMatch (nmatch, in_ximage, in_yimage, in_xworld, in_yworld,
                 out_ximage, out_yimage, 0.0, 0.0, wcs, 1);
     */

     FitMatch (nmatch, in_ximage, in_yimage, in_xworld, in_yworld, wcs, 1);
     
     /* The found solution seems good, but no check whatsoever is made
	(the routine always gives a result, even if unreasonable) Some kind of check on the 
	residuals should be done */
     
     /* Compute residuals at each star location */
     for (i = 0; i < nmatch; i++) {
       pix2vimoswcs (wcs, in_ximage[i], in_yimage[i], &mx, &my);
       xe[i] = (mx - in_xworld[i]) * 3600.0;
       ye[i] = (my - in_yworld[i]) * 3600.0;
       resid[i] = sqrt (xe[i]*xe[i] + ye[i]*ye[i]);
       cpl_msg_debug (modName,"%3d (%12.8f,%12.8f) -> %12.8f %12.8f %6.3f %6.3f %6.3f\n",
		i, in_xworld[i], in_yworld[i], mx, my, xe[i], ye[i], resid[i]);
       xsum = xsum + xe[i];
       ysum = ysum + ye[i];
       rsum = rsum + resid[i];
     }
     
     /* Compute means and standard deviations */
     xmean = xsum / dmatch;
     ymean = ysum / dmatch;
     rmean = rsum / dmatch;
     xsumsq = 0.0;
     ysumsq = 0.0;
     for (i = 0; i < nmatch; i++) {
       diff = xe[i] - xmean;
       xsumsq = xsumsq + (diff * diff);
       diff = ye[i] - ymean;
       ysumsq = ysumsq + (diff * diff);
     }
     xsig = sqrt (xsumsq / dmatch1);
     ysig = sqrt (ysumsq / dmatch1);
     rsig = sqrt ((xsumsq + ysumsq)/ dmatch1);
     cpl_msg_debug (modName,"Mean x: %12.8f/%12.8f y: %12.8f/%12.8f r: %12.8f/%12.8f\n",
	      xmean, xsig, ymean, ysig, rmean, rsig);     
     if (xsig < rmsTol || ysig < rmsTol ) {
       break;
     }
     
     cpl_msg_warning (modName,"Residuals for WCS fit exeed the limit, "
		    "discarding and iterating.");
     siglim = 2.0 * rsig;

     /* sort by increasing total residual */
     for (i = 0; i < nmatch-1; i++) {
       for (j = i+1; j < nmatch; j++) {
	 if (resid[j] < resid[i]) {
	   tmp = in_ximage[i]; in_ximage[i] = in_ximage[j]; in_ximage[j] = tmp;
	   tmp = in_yimage[i]; in_yimage[i] = in_yimage[j]; in_yimage[j] = tmp;
	   tmp = in_xworld[i]; in_xworld[i] = in_xworld[j]; in_xworld[j] = tmp;
	   tmp = in_yworld[i]; in_yworld[i] = in_yworld[j]; in_yworld[j] = tmp;
	   tmp = resid[i]; resid[i] = resid[j]; resid[j] = tmp;
	 }
       }
     }
     
     nmatch1 = nfit+1;
     /* Cut off points at residual of two sigma */
     for (i = nmatch-1; i >= nmatch1 ; i--) {
       if (resid[i] <= siglim ) {
	 nmatch1 = i+1;
	 break;
       }
     }
     nmatch=nmatch1;
     iterate++;
     
     
   } /* end of while */

   cpl_msg_info(modName, "Final RMS of fitted shift and rotation: (x, y) = (%f, %f)",
	      xsig, ysig);

   if (xsig > rmsTol || ysig > rmsTol ) {
     cpl_msg_error (modName,"Could not reach a reasonable fit.");
     return VM_FALSE;
   }
   

   /* Free memory */
   
   if(in_ximage) cpl_free(in_ximage);
   if(in_yimage) cpl_free(in_yimage);
   if(in_xworld) cpl_free(in_xworld);
   if(in_yworld) cpl_free(in_yworld);
   
   return VM_TRUE;
}

/**
 * @memo 
 *   Fit the matching observed and catalog stars with a polynomial
 *    to find plate solution. 
 *   
 * @return VM_TRUE/VM_FALSE. The input wcs is updated
 * 
 * @param wcs         wcs structure to update 
 * @param o_star      star table
 * @param a_star      astrometric Table
 * @param nmatch      number of matching stars
 * @doc
 *   Fit the matched observed and catalog stars with a SAO polynomial:
 *   it is simply an interface to FitPlate wcslib function
 *
 * @author B.Garilli
 */



VimosBool vimosFitPlate(struct WorldCoor * wcs,
           VimosTable * o_star,
           VimosTable * a_star,
           int      nmatch,
           int      ncoeff,
           double * chisq)
{
  
  int      i;
  double * x_1=0;
  double * y_1=0;
  double * x_2=0;
  double * y_2=0;
  double * sig=0;
  double * xfit_coeff;
  double * yfit_coeff;
  int    * ifit;
  double   chi2;
  char     modName[] = "vimosFitPlate";
  VimosColumn *oXimaCol, *oYimaCol, *aXwldCol, *aYwldCol;
  double   xsp, ysp, dx, dy;


  
 
   if(!(aXwldCol = findColInTab(a_star, "RA"))) {
     cpl_msg_error(modName, "Astrometric Table: Column with RA coord "
		 "not found");
     return VM_FALSE;
   }
   if(!(aYwldCol = findColInTab(a_star, "DEC"))) {
     cpl_msg_error(modName, "Astrometric Table: Column with RA coord "
		 "not found");
     return VM_FALSE;
   }
   
   if(!(oXimaCol = findColInTab(o_star, "X_IMAGE"))) {
     cpl_msg_error(modName, "Star Table: Column with X-pixel coord "
		 "not found");
     return VM_FALSE;
   }
   if(!(oYimaCol = findColInTab(o_star, "Y_IMAGE"))) {
     cpl_msg_error(modName, "Star Table: Column with Y-pixel coord "
		 "not found");
     return VM_FALSE;
   }
   
 
   
   if(!(x_1 = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ||
      !(y_1 = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ||
      !(sig = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ||
      !(x_2 = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ||
      !(y_2 = (double *)cpl_calloc((size_t)nmatch , sizeof (double))) ||
      !(ifit = (int *)cpl_calloc ((size_t)VMMAXPAR , sizeof(int))) ||
      !(xfit_coeff = (double *)cpl_calloc((size_t)VMMAXPAR , sizeof(double))) ||
      !(yfit_coeff = (double *)cpl_calloc((size_t)VMMAXPAR , sizeof(double))) ){
      cpl_msg_error(modName,"Could not alloc memory for finding plate solution");
      return VM_FALSE;
   }
   
   /* could not really make the wcs functions running silent!! */

   for(i=0; i<nmatch; i++) {
     x_1[i] = (oXimaCol->colValue->dArray[i]);
     y_1[i] = (oYimaCol->colValue->dArray[i]);
     x_2[i] = (aXwldCol->colValue->dArray[i]);
     y_2[i] = (aYwldCol->colValue->dArray[i]);
     sig[i] = 0.5;
   }
   /* warning: in FitPlate comment there is a mistake: the
      correct calling sequence is: 
      wcs structure, image coordinates, wordls coordinates, etc */

   {

     /*
      * This patch is added by C.Izzo to adapt coordinates
      * spreading around the 360 degrees line. Without this
      * the fit fails. Note that original tables remain 
      * unmodified.
      */

     int gap360 = 0;

     for(i=1; i<nmatch; i++) {
       if (fabs(x_2[i] - x_2[i-1]) > 250.) {
         gap360 = 1;
         break;
       }
     }
     if (gap360) {
       for(i=0; i<nmatch; i++) {
         if (x_2[i] < 250.) {
           x_2[i] += 360.;
         }
       }
     }
   }

   if(!(FitPlate(wcs, x_1, y_1, x_2, y_2, nmatch, ncoeff,0))) {
     chi2 = 0; 
     for (i = 0; i < nmatch; i++) {
       pix2vimoswcs (wcs, x_1[i], y_1[i], &xsp, &ysp);
       dx =3600.0 * (xsp - x_2[i]);
       dy = 3600.0 * (ysp - y_2[i]);
       chi2 += dx*dx + dy*dy;
     }
     *chisq = (chi2)/(double)MAX(nmatch-ncoeff,1);
     } else {
       cpl_msg_error(modName, "Cannot fit wcs plate model");
       return VM_FALSE;
     }
 
   /* Free memory */
   
   if(x_1) cpl_free(x_1);
   if(y_1) cpl_free(y_1);
   if(sig) cpl_free(sig);
   if(x_2) cpl_free(x_2);
   if(y_2) cpl_free(y_2);
   if(ifit) cpl_free(ifit);
   if(xfit_coeff) cpl_free(xfit_coeff);
   if(yfit_coeff) cpl_free(yfit_coeff);
   
   return VM_TRUE;
}


/*
 * @memo 
 *   Compute temperature scale factor.
 *
 * @param scale     Location where to store the temperature scale factor
 * @param dscList   Descriptor list
 * @param flag      Flag to enable/disable temperature check
 * @param tolerance Temperature tolerance value
 *
 * @return The function returns 0 on success and 1 otherwise.
 * 
 * Computes the temperature scale factor from the beam temperature. The 
 * beam temerature is cross-checked with the ambient temperature. In case
 * the beam temperature deviates more than @em tolerance from the ambient
 * temperature the ambient temperature is used for the computation.
 *
 * The computed temperature scale factor is stored in the variable pointed
 * to by the @em scale argument.
 */    

static int
vmGetTemperatureScale(double *scale, VimosDescriptor *dscList,
                      unsigned int flag, double tolerance)
{

    const char *id = "vmGetTemperatureScale";

    register const char *key;

    int q;

    double tbeam;         /* Current instrument temperature */
    double treference;    /* Instrument temperature at matrix creation time */
    double tcoeff = 6.0e-4;


    key = pilTrnGetKeyword("Quadrant");
    if (readIntDescriptor(dscList, key, &q, NULL) == VM_FALSE) {
        cpl_msg_error(id, "Missing keyword `%s'", key);
        return 1;
    }

    key = pilTrnGetKeyword("CcdSkyTemp");
    if (readDoubleDescriptor(dscList, key, &treference, NULL) == VM_FALSE) {
        cpl_msg_error(id, "Missing keyword `%s'", key);
        return 1;
    }

    if (flag) {

        /*
         * Get "safe" beam temperature.
         */

        if (getBeamTemperature(dscList, &tbeam, tolerance, q))
          return 1;

    }
    else {

        key = pilTrnGetKeyword("BeamTemperature", q);
        if (readDoubleDescriptor(dscList, key, &tbeam, NULL) == VM_FALSE) {
            cpl_msg_error(id, "Missing keyword `%s'", key);
            return 1;
        }

    }

    /*
     * Compute temperature scale factor.
     * Note: This is copied from vmmcs!
     */

    *scale = 1. + tcoeff * (tbeam - treference);

    return 0;


}


/*
 * @memo 
 *   Read the CCD to Sky matrix from a descriptor list
 *   
 * @param descs    image descriptors
 * @param x_coeff  Ccd Sky To x coeffs (returned)
 * @param y_coeff  Ccd To  Sky  y coeffs (returned)
 *
 * @return VM_TRUE/VM_FALSE
 * 
 * Reads the CCD to Sky transformation matrix from a descriptor list.
 *
 * @author B.Garilli
 */    

static VimosBool
getCcdSky(VimosDescriptor *descs, double *x_coeff, double *y_coeff)
{ 

    const char *id = "getCcdSky";

    const char *key;
    char valcoef[80] = "0.";

    int quad, i, k, j;
    int ncoeffx, ncoeffy;



    key = pilTrnGetKeyword("Quadrant");
    if ((readIntDescriptor(descs, key, &quad, NULL)) == VM_FALSE) {
        cpl_msg_error(id, "Missing keyword `%s'", key);
        return VM_FALSE;
    }


    /*
     * Read the coefficients of the CCD to Sky transformation
     */

    key = pilTrnGetKeyword("CcdSkyXord");
    if (readIntDescriptor(descs, key, &ncoeffx, NULL) == VM_FALSE) {
        cpl_msg_error(id, "Missing keyword `%s'", key);
        return VM_FALSE;
    }

    key = pilTrnGetKeyword("CcdSkyYord");
    if (readIntDescriptor(descs, key, &ncoeffy, NULL) == VM_FALSE) {
        cpl_msg_error(id, "Missing keyword %s", key);
        return VM_FALSE;
    }
 
    k=0;
    for (i = 0; i <= ncoeffx; i++) {
        for (j = 0; j <= ncoeffx; j++, k++) {
            key = pilTrnGetKeyword("CcdSkyX", i, j);
            if (readStringDescriptor(descs, key, valcoef, NULL) == VM_FALSE) {
                cpl_msg_error(id, "Missing keyword `%s'", key);
                return VM_FALSE;
            }
            x_coeff[k] = atof(valcoef);
        }
    }

    k=0;
    for (i = 0; i <= ncoeffy; i++) {
        for (j = 0; j <= ncoeffy; j++, k++) {
            key = pilTrnGetKeyword("CcdSkyY", i, j);
            if (readStringDescriptor(descs, key, valcoef, NULL) == VM_FALSE) {
                cpl_msg_error(id, "Missing keyword `%s'", key);
                return VM_FALSE;
            }
            y_coeff[k] = atof(valcoef);
        }
    }

    return VM_TRUE;

}


/*-----------------------------------------------------------------------------
 * Function  : rdimage
 * In        : VimosImage
 * Out       : wcs structure
 *
 * Purpouse  : Compute wcs structure from image header. If flag is set 
 *             write also CO matrix in wcs structure
 *---------------------------------------------------------------------------*/

struct WorldCoor *rdimage(VimosDescriptor *descs)
{
   struct WorldCoor * wcs;
   int    i,j,k,naxes,nxpix,nypix;
   double cd[4],equinox,epoch;
   char   comment[80];
   char   modName[]  = "rdimage";
   char       *descName, *ctype1, *ctype2;


   /* Evaluate image data format */
 
   if ((descs == NULL)) {
     cpl_msg_error(modName, "Null input image");
     return NULL;
   }
   

   wcs = (struct WorldCoor *)cpl_calloc(1, sizeof(struct WorldCoor));

   /* Set WCSLIB flags so that structures will be reinitialized */
   wcs->cel.flag = 0;
   wcs->lin.flag = 0;
   wcs->vimoswcsl.flag = 0;

   /* Initialize to no plate fit */
   wcs->ncoeff1 = 0;
   wcs->ncoeff2 = 0;

   /* Initialize to no CD matrix */
   cd[0] = 0.0;
   wcs->rotmat = 0;
   wcs->rot = 0.0;
   
   
   /* Header parameters independent of projection */
   
   naxes = 2;
   wcs->naxes = naxes;
   wcs->lin.naxis = naxes;
   wcs->nxpix = 0;
   
   
   /* Read image WCS header keywords and if missing 
      substitute with corresponding filter table values */
   
   if(!(readIntDescriptor(descs, pilTrnGetKeyword("Naxis",1),
                          &nxpix, comment))) {
     cpl_msg_error(modName,"Descriptor NAXIS not found");
     return NULL;
   }
   if(!(readIntDescriptor(descs, pilTrnGetKeyword("Naxis",2),
                          &nypix, comment))) {
     cpl_msg_error(modName,"Descriptor NAXIS not found");
     return NULL;
   }

   wcs->nxpix = (double) nxpix;
   wcs->nypix = (double) nypix;

   descName = (char *)pilKeyTranslate("Ctype",1);
   ctype1 = (char *)cpl_malloc(9*sizeof(char)); //ctype[4][9] in vimoswcs.h
   if(!(readStringDescriptor(descs, descName,
                             ctype1, comment))) {
    cpl_msg_error(modName,"Descriptor %s not found", descName);
    cpl_free(descName);

    return NULL;
   }
   cpl_free(descName);
   descName = (char *)pilKeyTranslate("Ctype",2);
   ctype2 = (char *)cpl_malloc(9*sizeof(char)); //ctype[4][9] in vimoswcs.h
   if(!(readStringDescriptor(descs, descName,
                             ctype2, comment))) {
    cpl_msg_error(modName,"Descriptor %s not found", descName);
    cpl_free(descName);

    return NULL;
   }
   cpl_free(descName);


   strcpy (wcs->ctype[0], ctype1);
   strcpy (wcs->ctype[1], ctype2);
   
   /* Set projection type in WCS data structure */
   if (vimoswcstype (wcs, ctype1, ctype2)) {
     vimoswcsfree (wcs);
     return (NULL);
   }
   
   /* Reference pixel coordinates and WCS value */
   
   wcs->crpix[0] = 1.0;
   wcs->crpix[1] = 1.0;
   
   if (!(readDoubleDescriptor(descs, pilTrnGetKeyword("Crpix",1),
                              &(wcs->crpix[0]), comment))){
     cpl_msg_error(modName, "Cannot find CRPIX descriptors");
     return NULL;
   }
   if (!(readDoubleDescriptor(descs, pilTrnGetKeyword("Crpix",2),
                              &(wcs->crpix[1]), comment))){
     cpl_msg_error(modName, "Cannot find CRPIX descriptors");
     return NULL;
   }

    
   wcs->xrefpix = wcs->crpix[0];
   wcs->yrefpix = wcs->crpix[1];
   wcs->crval[0] = 0.0;
   wcs->crval[1] = 0.0;

   if (!(readDoubleDescriptor(descs, pilTrnGetKeyword("Crval",1),
                              &(wcs->crval[0]), comment))){
     cpl_msg_error(modName, "Cannot find world coordinates of "
                 "telescope pointing");
     return NULL;
   }
   if (!(readDoubleDescriptor(descs, pilTrnGetKeyword("Crval",2),
                              &(wcs->crval[1]), comment))){
     cpl_msg_error(modName, "Cannot find world coordinates of "
                 "telescope pointing");
     return NULL;
   }

   wcs->xref = wcs->crval[0];
   wcs->yref = wcs->crval[1];
   if (wcs->coorflip) {
     wcs->cel.ref[0] = wcs->crval[1];
     wcs->cel.ref[1] = wcs->crval[0];
   }else {
     wcs->cel.ref[0] = wcs->crval[0];
     wcs->cel.ref[1] = wcs->crval[1];
   }
   wcs->longpole = 999.0;
   wcs->cel.ref[2] = wcs->longpole;
   wcs->latpole = 999.0;
   wcs->cel.ref[3] = wcs->latpole;
   wcs->lin.crpix = wcs->crpix;
   wcs->lin.cdelt = wcs->cdelt;
   wcs->lin.pc = wcs->pc;

   wcs->prj.r0 = 0.0;


   /* Get CD matrix */
   k = 0;
   for(i=0; i<2; i++) {     
     for(j=0; j<2; j++) {
       if(readDoubleDescriptor(descs,pilTrnGetKeyword("CD",i+1,j+1),
                               &cd[k], comment)) 
         k++;
       else {     
         cpl_msg_error(modName, "Cannot find CD matrix");
         return NULL;
       }
     }
   }
   
   vimoswcscdset (wcs, cd);


   /* Coordinate reference frame, equinox, and epoch */


   if(!readDoubleDescriptor(descs, pilTrnGetKeyword("Equinox"), &equinox,
                            comment)) {
     cpl_msg_error(modName, "Cannot find EQUINOX in image header");
     return NULL;
   }
   readDoubleDescriptor(descs, "EPOCH", &epoch, comment);

    wcs->equinox =  equinox;
    if (equinox > 1980)
        strcpy (wcs->radecsys,"FK5");
    else
        strcpy (wcs->radecsys,"FK4");
    if (epoch > 0)
        wcs->epoch = epoch;
    else
        wcs->epoch = 0.0;


   wcs->vimoswcson = 1;

   wcs->lin.crpix = wcs->crpix;
   wcs->lin.cdelt = wcs->cdelt;
   wcs->lin.pc = wcs->pc;
   if (strlen (wcs->radecsys) == 0 || wcs->prjcode == VIMOSWCS_LIN)
        strcpy (wcs->radecsys, "LINEAR");
   wcs->sysvimoswcs = vimoswcscsys (wcs->radecsys);

   if (wcs->sysvimoswcs == VIMOSWCS_B1950)
        strcpy (wcs->radecout, "FK4");
    else if (wcs->sysvimoswcs == VIMOSWCS_J2000)
        strcpy (wcs->radecout, "FK5");
    else
        strcpy (wcs->radecout, wcs->radecsys);
    wcs->sysout = vimoswcscsys (wcs->radecout);
    wcs->eqout = wcs->equinox;
    strcpy (wcs->radecin, wcs->radecsys);
    wcs->sysin = vimoswcscsys (wcs->radecin);
    wcs->eqin = wcs->equinox;
    wcs->printsys = 1;
    wcs->tabsys = 0;
    wcs->linmode = 0;
   


    cpl_msg_debug(modName,"NAXIS1 = %10.f; NAXIS2 = %10.f "
                "CRVAL1 = %10f; CRVAL2 = %10f; CRPIX1 = %10.3f; "
                "CRPIX2 = %10.3f; CDELT1 = %10.4g; CDELT2 = %10.4g; "
                "CROTA = %f; EPOCH = %f; EQUINOX = %f; RADECSYS = %s", 
                wcs->nxpix,wcs->nypix, wcs->crval[0], wcs->crval[1], 
                wcs->crpix[0], wcs->crpix[1],wcs->cdelt[0], wcs->cdelt[1],
                wcs->rot, wcs->epoch,wcs->equinox,wcs->radecsys);
   
   return (wcs);
}


/**
 * @memo 
 *   Correct a Sky to CCD transformation matrix for temperature effects
 *   
 * @param descs      Image descriptor list
 * @param ostar      Table with grid of positions computed without
 *                   plate solution
 * @param flag       Flag to enable/disable temperature check
 * @param tolerance  Tolerance value for temperature check.
 *
 * @return VM_TRUE/VM_FALSE
 * 
 * The function applies corrections for instrument distortions and
 * temperature variation effects to a list of X and Y image coordinates
 * found in the table @em ostar. The table @em ostar is updated with the
 * corrected positions.
 *
 * If @em flag is different from 0, the beam temperature used for the
 * computation of the temperature corrections is checked against the
 * ambient temperature. The beam temperature is used only if it differs
 * not more than @em tolerance from the ambient temperature, otherwise
 * the ambient temperature is used for the computation of the corrections.
 *
 * Note: The validation of the beam temperature assumes that the ambient
 *       temperature value is correct if the keyword entry is found in
 *       the descriptor list @em descs.
 */    

VimosBool
computeVirtualPixels(VimosDescriptor *descs, VimosTable *ostar,
                     unsigned int flag, double tolerance)
{ 

  const char *id = "computeVirtualPixels";

  register const char *key;

  int i, j, k, n, nrows;
  int nx, ny;
  int xCenter, yCenter;

  double *x_coeff, *y_coeff;
  double x, y;
  double scale = 1.;

  VimosColumn *xpixCol, *ypixCol;


  nrows = colGetSize(ostar->cols);
  if (nrows == 0) {
      cpl_msg_error(id, "Source list is empty!");
      return VM_FALSE;
  }


  /*
   * Compute scale factor of temperature variations
   */

  if (vmGetTemperatureScale(&scale, descs, flag, tolerance)) {
      cpl_msg_error(id, "Cannot compute temperature corrections!");
      return VM_FALSE;
  }


  /*
   * Read CCD to sky matrix coefficients
   */

  key = pilTrnGetKeyword("CcdSkyXord");
  if (readIntDescriptor(descs, key, &nx, NULL) == VM_FALSE) {
      cpl_msg_error(id, "Missing keyword `%s'", key);
      return VM_FALSE;
  }

  key = pilTrnGetKeyword("CcdSkyYord");
  if (readIntDescriptor(descs, key, &ny, NULL) == VM_FALSE) {
      cpl_msg_error(id, "Missing keyword %s", key);
      return VM_FALSE;
  } 

  x_coeff = (double *)cpl_calloc(ipow(nx + 1, 2), sizeof(double));
  y_coeff = (double *)cpl_calloc(ipow(ny + 1, 2), sizeof(double));

  if (getCcdSky(descs, x_coeff, y_coeff) == VM_FALSE) {
      cpl_msg_error(id, "Cannot read CCD to sky transformation matrix!");

      cpl_free(x_coeff);
      cpl_free(y_coeff);

      return VM_FALSE;
  }


  /*
   * We compute "virtual" pixel coordinates, i.e. we correct pixel
   * coordinates for distortions and temperature variations.
   *
   * Note: Corrections are relative to CCD central point computed
   *       from the NAXIS keyword!
   */

  key = pilTrnGetKeyword("Naxis", 1);
  if (readIntDescriptor(descs, key, &xCenter, NULL) == VM_FALSE) {
      cpl_msg_error(id, "Missing keyword `%s'", key);

      cpl_free(x_coeff);
      cpl_free(y_coeff);

      return VM_FALSE;
  }

  key = pilTrnGetKeyword("Naxis", 2);
  if (readIntDescriptor(descs, key, &yCenter, NULL) == VM_FALSE) {
      cpl_msg_error(id, "Missing keyword `%s'", key);

      cpl_free(x_coeff);
      cpl_free(y_coeff);

      return VM_FALSE;
  }

  xCenter = xCenter / 2;
  yCenter = yCenter / 2;

  xpixCol = findColInTab(ostar, "X_IMAGE");
  ypixCol = findColInTab(ostar, "Y_IMAGE");

  for (n = 0; n < nrows; n++) {
      x=xpixCol->colValue->dArray[n];
      y=ypixCol->colValue->dArray[n];
      xpixCol->colValue->dArray[n] = -xCenter;
      ypixCol->colValue->dArray[n] = -yCenter;

      k=0;
      for (i = 0; i <= nx; i++) {
          for (j = 0; j <= nx; j++,k++) {
              xpixCol->colValue->dArray[n] += x_coeff[k] * ipow(x,j) *
                  ipow(y,i);  
          }
      }
      xpixCol->colValue->dArray[n] = xpixCol->colValue->dArray[n] * scale +
          xCenter;

      k=0;
      for (i = 0; i <= ny; i++) {
          for (j = 0; j <= ny; j++,k++) {
              ypixCol->colValue->dArray[n] += y_coeff[k] * ipow(x,j) *
                  ipow(y,i);  
          }
      }
      ypixCol->colValue->dArray[n] = ypixCol->colValue->dArray[n] * scale +
          yCenter;

  }


  cpl_free(x_coeff);
  cpl_free(y_coeff);

  return VM_TRUE;

}


/**
 * @memo 
 *   Update image header 
 *   
 * @return VM_TRUE/VM_FALSE
 * 
 * @param image     image to update
 * @param wcs       wcs structure
 * @param rms       rms array
 *
 * @doc 
 *   Update image header with the new projection coefficients.
 *
 * @author P. Montegriffo (i/o modified by P.Sartoretti)
 */     
/*-----------------------------------------------------------------------------
 * Function  : upheader
 * In        : image & wcs structures
 * Out       : Updated image header 
 *
 * Purpouse  : Updates image header 
 *---------------------------------------------------------------------------*/

int upheader (VimosImage *image, struct WorldCoor * wcs, double rms[])
{
  int    i,j, k=0;
  char modName[] = "upheader";
  

 

  if(!(writeDoubleDescriptor(&(image->descs), 
                             pilTrnGetKeyword("PixelScale"), 
                             (fabs(wcs->cdelt[0]) *3600.0),
                             "pixel scale of telescope in CCD plane"))) {
    cpl_msg_error(modName,"Cannot write pixelscale");
    return VM_FALSE;
  }
  
  if(!(writeDoubleDescriptor(&(image->descs),pilTrnGetKeyword("Crpix",1),
       (wcs->xrefpix), "X position of telescope optical axis in pixels"))) {
    cpl_msg_error(modName, "Cannot write CRPIX in image header ");
    return VM_FALSE;
  }
  if(!(writeDoubleDescriptor(&(image->descs),pilTrnGetKeyword("Crpix",2),
       (wcs->yrefpix), ""))) {
    cpl_msg_error(modName, "Cannot write CRPIX2 in image header");
    return VM_FALSE;
  }
  if(!(writeDoubleDescriptor(&(image->descs),pilTrnGetKeyword("Crval",1),
       wcs->xref, ""))) {
    cpl_msg_error(modName, "Cannot update CRVAL1 in image header");
    return VM_FALSE;
  }
  if(!(writeDoubleDescriptor(&(image->descs),pilTrnGetKeyword("Crval",2),
       wcs->yref, ""))) {
    cpl_msg_error(modName, "Cannot update CRVAL2 in image table");
    return VM_FALSE;
  }
  
  for (i=1; i<=2; i++) {
    for(j=1; j<=2; j++) {
      if(writeDoubleDescriptor(&(image->descs), 
                               pilTrnGetKeyword("CD",i,j),
                               wcs->cd[k], "FITS CD transformation matrix")){
        k++;
      }
      else {
        cpl_msg_error(modName, "Cannot write CD matrix in image");
        return VM_FALSE;
      }
    }
  }

  for(i=0; i<wcs->ncoeff1; i++) {
    if(!(writeDoubleDescriptor(&(image->descs), 
                          pilTrnGetKeyword("CO1",i+1), 
                        wcs->x_coeff[i], "X transformation coefficients"))) {
      cpl_msg_error(modName, "Cannot write X transformation coefficients "
                  "in header");
      return VM_FALSE;
    }
  }
  for(i=0; i<wcs->ncoeff2; i++) {
    if(!(writeDoubleDescriptor(&(image->descs), 
                               pilTrnGetKeyword("CO2",i+1), 
                      wcs->y_coeff[i], "Y transformation coefficients"))) {
      cpl_msg_error(modName, "Cannot write ccd-sky Y transformation "
                  "coefficients in image header");
      return VM_FALSE;
    }
  }
 
  if (!(writeDoubleDescriptor(&(image->descs), 
                              pilTrnGetKeyword("InvCO1rms"), rms[0], 
                              ""))) {
    cpl_msg_error(modName, "Cannot write ccd-sky X RMS in header");
    return VM_FALSE;
  } 
  if (!(writeDoubleDescriptor(&(image->descs), 
                              pilTrnGetKeyword("InvCO2rms"), rms[1], 
                              ""))) {
    cpl_msg_error(modName, "Cannot write ccd-sky  Y RMS in header");
    return VM_FALSE;
  }
  if (!(writeDoubleDescriptor(&(image->descs), 
                              pilTrnGetKeyword("CO1rms"), rms[2], 
                              ""))) {
    cpl_msg_error(modName, "Cannot write sky-ccd X RMS in header");
    return VM_FALSE;
  } 
  if (!(writeDoubleDescriptor(&(image->descs), 
                              pilTrnGetKeyword("CO2rms"), rms[3], 
                              ""))) {
    cpl_msg_error(modName, "Cannot write ccd-sky  Y RMS in header");
    return VM_FALSE;
  }

   return VM_TRUE;
}
