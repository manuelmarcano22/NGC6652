/* $Id: vmmosflat.c,v 1.3 2013-08-22 16:58:58 cgarcia Exp $
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
 * $Date: 2013-08-22 16:58:58 $
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmgrismtable.h"
#include "vmextractiontable.h"
#include "vmadf.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmmosflat.h"
#include "cpl.h"


#define EXTRA              20.
#define REBIN_LOWER_THRESH 1.e-4


/**
 * @name vmmosflat MOS Flat Fields
 *
 * The module provides the functions for creating and applying a
 * spectroscopic master spectral flat field.
 */

/**@{*/

/**
 * @memo
 *   Normalize a spectral Master Flat Field.
 *
 * @return Pair of images, the first the normalized flat field, the
 *         second a misterious "count" image that nobody apparently
 *         know about its meaning.
 *
 * @param imageData       Master spectral Flat to be normalized.
 * @param extractionTable Input Extraction Table
 * @param polyDegX        Degree of polynomial fitting any spectral
 *                        flat field trend along the X (spatial)
 *                        direction.
 * @param polyDegY        Degree of polynomial fitting any spectral
 *                        flat field trend along the Y (dispersion)
 *                        direction.
 */

VimosImage *
VmSpNormFF(VimosImage *imageData, VimosExtractionTable *extractionTable,
           int normMethod, int polyDegX, int polyDegY, FilterMethod method,
           int xSize, int ySize)
{
  VimosImage *ffData;

  if (normMethod == 1)
    ffData = VmSpNormPoly(imageData, extractionTable, polyDegX, polyDegY);
  else 
    ffData = VmSpNormSmooth(imageData, extractionTable, xSize, ySize, method);

  return ffData;
}


/**
 * @memo
 *   Normalize a spectral Master Flat Field.
 *
 * @return Normalized flat field.
 *
 * @param imageData       Master spectral Flat to be normalized.
 * @param extractionTable Input Extraction Table
 * @param polyDegX        Degree of polynomial fitting any spectral
 *                        flat field trend along the X (spatial)
 *                        direction.
 * @param polyDegY        Degree of polynomial fitting any spectral
 *                        flat field trend along the Y (dispersion)
 *                        direction.
 */

VimosImage *
VmSpNormPoly(VimosImage *imageData, VimosExtractionTable *extractionTable,
             int polyDegX, int polyDegY)
{
  int i, j, k;
  int xOut, yOut;
  int numRows;
  int numPoints;
  VimosUlong32 index;
  int nPixBelow, nPixAbove;
  int imageXlen, imageYlen;
  int ycount=0;
  
  double *coeffsX = NULL;
  double *coeffsY = NULL;
  double ffVal;
  double ffValX;
  double ffValY;
  double ffSum;
  double datVal;
  double xp;
  double xOutF;
  double frac;
    
  VimosImage  *ffData;
  VimosImage  *countData;
  
  VimosDpoint *bufferX = NULL;
  VimosDpoint *bufferY = NULL;
 
  VimosExtractionSlit *slit;
  
  char modName[]   ="VmSpNormPoly";



  cpl_msg_debug(modName, "Normalize Flat Field");


  imageXlen = imageData->xlen;
  imageYlen = imageData->ylen;

  countData = newImageAndAlloc(imageXlen, imageYlen);
  ffData = newImageAndAlloc(imageXlen, imageYlen);
  
  for (xOut = 0; xOut < imageXlen; xOut++) {
    for (yOut = 0; yOut < imageYlen; yOut++) {
      ffData->data[xOut+imageXlen*yOut] = 0.0;
      countData->data[xOut+imageXlen*yOut] = 0.0;
    }
  }
  
  slit = extractionTable->slits;

  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN LO", 
                    &nPixBelow, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN HI",
                    &nPixAbove, NULL);
  numPoints = nPixBelow+nPixAbove+1;
  
  bufferY = newDpoint(numPoints);
  
  while (slit) {
    
    numRows = slit->numRows;

    for (i = 0 ; i < numPoints; i++) {
      bufferY[i].x = i;
      bufferY[i].y = 0.0;
    }
    
    deleteDpoint(bufferX);
    bufferX = newDpoint(numRows);
    for (i = 1; i < numRows-1; i++) {
      bufferX[i-1].x = slit->ccdX->data[i];
      bufferX[i-1].y = 0.0;
    }
    ffSum = 0.0;
    
    for (i = 2; i < numRows-2; i++) {
      ycount = 0;
      for (j = -nPixBelow+1; j < nPixAbove; j++) { 

        yOut = slit->ccdY->data[i] + j;
        xOutF = slit->ccdX->data[i]+computeDistModel1D(slit->crvPol[i],yOut);
        xOut = xOutF;
        frac = xOutF-xOut;
        if (xOut >= 0 && xOut+1 < imageXlen && yOut >=0 && yOut < imageYlen) {
          index = xOut+yOut*imageXlen;
          datVal = ( (1.0-frac)*imageData->data[index] +
                       frac*imageData->data[index+1] );
          bufferX[i-2].y += datVal;
          bufferY[ycount].y += datVal;
          ycount++;
          ffSum += datVal;
        }
      }
    }
    if (coeffsX != NULL) {
      cpl_free(coeffsX);
      coeffsX = NULL;
    }
    if (coeffsY != NULL) {
      cpl_free(coeffsY);
      coeffsY = NULL;
    }
    
    coeffsX = fit1DPoly(polyDegX, bufferX, numRows-4, NULL);
    coeffsY = fit1DPoly(polyDegY, bufferY, ycount, NULL);

    for (i = 1; i < numRows-1; i++) {
      ycount = 0;
      for (j = -nPixBelow+1; j < nPixAbove; j++) {
        
        yOut = slit->ccdY->data[i] + j;
        xOutF = computeDistModel1D(slit->crvPol[i], yOut);
        xOutF = slit->ccdX->data[i]+xOutF;
        xOut = xOutF;

        if ( (xOut >= 0) && (xOut < imageXlen) &&
             (yOut >= 0) && (yOut < imageYlen) ) {
          frac = xOutF-xOut;
          ffValX = coeffsX[0];
          for (k = 1; k <= polyDegX; k++) {
            xp = ipow(slit->ccdX->data[i], k) ;
            ffValX += coeffsX[k] * xp ; 
          }
          
          ffValY = coeffsY[0];
          for (k = 1; k <= polyDegY; k++) {
            xp = ipow(1.0*(ycount+1), k) ;
            ffValY += coeffsY[k] * xp ; 
          }
          ycount++;
          ffVal = ffValX*ffValY/ffSum;

          index = xOut+yOut*imageXlen;
          ffData->data[index] += (1.0-frac)*ffVal;
          ffData->data[index+1] += frac*ffVal;
          countData->data[index] += 1.0-frac;
          countData->data[index+1] += frac;
        }
      }

      /* TBD CHECK ON FLAT_TOL and write new table */
    }
        
    slit = slit->next;
    
  }
  deleteDpoint(bufferX);
  deleteDpoint(bufferY);
  
  for (xOut = 0; xOut < imageXlen; xOut++) {
    for (yOut = 0; yOut < imageYlen; yOut++) {
      index = xOut+yOut*imageXlen;
      if (countData->data[index] > 0.0) {
        ffData->data[index] /= countData->data[index];
        ffData->data[index] = imageData->data[index] / ffData->data[index];
      }
      else {
        ffData->data[index] = 1.0;
      }
    }
  }
  
  copyAllDescriptors(imageData->descs, &(ffData->descs));

  deleteImage(countData);
  
  return ffData;
}


/**
 * @memo
 *   Normalize a spectral Master Flat Field.
 *
 * @return One image containing the first the normalized flat field.
 *
 * @param imageData       Master spectral Flat to be normalized.
 * @param extractionTable Input Extraction Table
 * @param xSize           size in x of the smoothing filter box
 * @param ySize           size in y of the smoothing filter box
 * @param method          smoothing filter method
 *
 */

VimosImage *
VmSpNormSmooth(VimosImage *imageData, VimosExtractionTable *extractionTable,
               int xSize, int ySize, FilterMethod method)
{
  int i, j;
  int xOut, yOut;
  int numRows;
  int numPoints;
  VimosUlong32 index;
  int nPixBelow, nPixAbove;
  int imageXlen, imageYlen;

  double xOutF;
  double frac;
  double smoothedValue;
    
  VimosImage *extractedSpectrum;
  VimosImage *smoothedSpectrum;
  VimosImage *countData;
  VimosImage *ffData;
 
  VimosExtractionSlit *slit;
  
  char modName[]   ="VmSpNormSmooth";



  cpl_msg_info(modName, "Normalize Flat Field");


  imageXlen = imageData->xlen;
  imageYlen = imageData->ylen;

  countData = newImageAndAlloc(imageXlen, imageYlen);
  ffData = newImageAndAlloc(imageXlen, imageYlen);
  
  for (xOut = 0; xOut < imageXlen; xOut++) {
    for (yOut = 0; yOut < imageYlen; yOut++) {
      countData->data[xOut+imageXlen*yOut] = 0.0;
      ffData->data[xOut+imageXlen*yOut] = 0.0;
    }
  }
  
  slit = extractionTable->slits;

  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN LO", 
                    &nPixBelow, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN HI",
                    &nPixAbove, NULL);
  numPoints = nPixBelow+nPixAbove+1;

  
  while (slit) {
    
    numRows = slit->numRows;

    if (numRows < 2) {
      slit = slit->next;
      continue;
    }

    extractedSpectrum = newImageAndAlloc(numRows, nPixBelow+nPixAbove+1);

    for (i = 0; i < numRows; i++) {
      for (j = -nPixBelow; j < nPixAbove+1; j++) { 

        yOut = slit->ccdY->data[i] + j;
        xOutF = slit->ccdX->data[i]+computeDistModel1D(slit->crvPol[i],yOut);
        xOut = xOutF;
        frac = xOutF-xOut;
        if (xOut >= 0 && xOut+1 < imageXlen && yOut >=0 && yOut < imageYlen) {
          index = xOut+yOut*imageXlen;

          /* "Warp" spectrum to a rectangular area */
          extractedSpectrum->data[i+(j+nPixBelow)*numRows] = 
              ((1.0-frac) * imageData->data[index] + 
               frac*imageData->data[index+1] );
        }
      }
    }

    /* smooth it */
    smoothedSpectrum = VmFrFilter(extractedSpectrum, 1, ySize, method, 0);

    /* "Warp" back to output image */
    for (i = 0; i < numRows; i++) {
      for (j = -nPixBelow; j < nPixAbove+1; j++) { 

        yOut = slit->ccdY->data[i] + j;
        xOutF = slit->ccdX->data[i]+computeDistModel1D(slit->crvPol[i],yOut);
        xOut = xOutF;
        frac = xOutF-xOut;
        if (xOut >= 0 && xOut+1 < imageXlen && yOut >=0 && yOut < imageYlen) {
          index = xOut+yOut*imageXlen;
          smoothedValue = smoothedSpectrum->data[i+(j+nPixBelow)*numRows];
          if ((1.0-frac)*smoothedValue > REBIN_LOWER_THRESH) {
            ffData->data[index] += (1.0-frac)*smoothedValue;
            countData->data[index] += 1.0-frac;
          }                      

          if (frac*smoothedValue > REBIN_LOWER_THRESH) {
            ffData->data[index+1] += frac*smoothedValue;
            countData->data[index+1] += frac;
          }
        }
      }
    }

    deleteImage(extractedSpectrum);
    deleteImage(smoothedSpectrum);

    slit = slit->next;
    
  }

  /* Correct weights and normalize */  
  for (xOut = 0; xOut < imageXlen; xOut++) {
    for (yOut = 0; yOut < imageYlen; yOut++) {
      index = xOut+yOut*imageXlen;
      if (countData->data[index] > 0.0) {
        ffData->data[index] /= countData->data[index];
        ffData->data[index] = imageData->data[index] / ffData->data[index];
      }
      else {
        ffData->data[index] = 1.0;
      }
    }
  }
  
  copyAllDescriptors(imageData->descs, &(ffData->descs));
  
  deleteImage(countData);

  return ffData;
}


VimosImage **
VmSpStackFF(VimosImage **flatImages, int numFlat,
            VimosExtractionTable *extrTable, int fuzz)
{
  char          modName[] = "VmSpStack";
  VimosTable   *adf;
  VimosAdfType adfType;
  VimosImage   *inputImage;
  VimosImage   **outputImages=NULL;
  VimosImage   *zeroOrderIma=NULL;
  VimosImage   *firstOrderIma;
  VimosDpoint  *minY;
  VimosDpoint  *maxY;
  VimosDpoint  *zeroY;
  VimosExtractionSlit  *slit;
  int           doZero;
  int           i, j;
  int           imaXlen;
  int           imaYlen;
  int           xOut;
  int           yOut;
  int           nSlits;
  int           initialize = 1;
  int           zWidth = 10;
  char          comment[80];


  cpl_msg_info(modName, "Stacking %d flat fields", numFlat);

  outputImages = (VimosImage **) cpl_malloc(2*sizeof(VimosImage *));

  if (readIntDescriptor(extrTable->descs, pilTrnGetKeyword("ZeroOrderFlag"), 
                        &doZero, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot find descriptor %s", 
                pilTrnGetKeyword("ZeroOrderFlag"));
    return NULL;
  }
  if (doZero && numFlat > 1)
  {
    for (i=0; i<numFlat; i++) 
    {
      inputImage = flatImages[i];
 
      imaXlen = inputImage->xlen;
      imaYlen = inputImage->ylen;

      /* create empty ADF */
      adfType = VM_ADF_TYPE_UDF;
      adf = newADF();
      /* read ADF from image */
      readADF(adf, inputImage);

      if (!strcmp(adf->name, VM_ADF_IMA)) {
        adfType = VM_ADF_TYPE_IMA;
      }
      if (!strcmp(adf->name, VM_ADF_MOS)) {
        adfType = VM_ADF_TYPE_MOS;
      }
      if (!strcmp(adf->name, VM_ADF_IFU)) {
        adfType = VM_ADF_TYPE_IFU;
      }
      if (adfType == VM_ADF_TYPE_UDF) {
        return NULL;
      }      

      if (adfType != VM_ADF_TYPE_MOS) {
	cpl_msg_error(modName, "Trying to stack exposures other than MOS...");
	return NULL;
      }
      else {
	/* the first time around, create output image */
	if (initialize) 
	{
	  zeroOrderIma = newImageAndAlloc(imaXlen,imaYlen);
	  copyAllDescriptors(adf->descs, &(zeroOrderIma)->descs);
	  firstOrderIma = newImageAndAlloc(imaXlen,imaYlen);
	  copyAllDescriptors(adf->descs, &(firstOrderIma)->descs);
	  initialize = 0;
	}

	/* determine frame region that is effectively exposed */
	determineExposedMosArea(adf,extrTable,&nSlits,&minY,&maxY,&zeroY);

	/* write in the output image the appropriate section of input image */
	for (j = 0; j < nSlits; j++) 
	{
	  for (xOut = floor(minY->x) - fuzz; 
	       xOut <= ceil(maxY->x) + fuzz; xOut++) 
	  {
	    for (yOut=floor(minY->y)-fuzz; yOut<=ceil(maxY->y)+fuzz; yOut++) 
	    {
	      if (xOut >= 0 && xOut < imaXlen && yOut >=0 && yOut < imaYlen) {
		firstOrderIma->data[xOut+imaXlen*yOut] = 
		  inputImage->data[xOut+imaXlen*yOut];
	      }
	    }
	    for(yOut=floor(zeroY->y-zWidth);yOut<=ceil(zeroY->y+zWidth); yOut++)
	    {
	      if (xOut >= 0 && xOut < imaXlen && yOut >=0 && yOut < imaYlen) {
		zeroOrderIma->data[xOut+imaXlen*yOut] = 
		  inputImage->data[xOut+imaXlen*yOut];
	      }
	    }
	  }
	  minY = minY->prev;
	  maxY = maxY->prev;
	  zeroY = zeroY->prev;
	}
	deleteADF(adf);
      }
    }

   /*
    * Check whether some slits are left in the extraction table. The
    * remaining slits are the ones that are never illuminated by any
    * of the identified shutter positions. I believe that this should
    * never happen, because in that case the next routines would try
    * to use such slits, looking for them in the stacked images to
    * try fitting the optical distorsion and the spectral curvature
    * models. (C.Izzo)
    */

    if (extrTable->slits) {
      slit = extrTable->slits;
      j = 0;
      while (slit) {
        j++;
        slit = slit->next;
      }
      cpl_msg_error(modName, "%d unexposed slits found!", j);
      deleteImage(firstOrderIma);
      deleteImage(zeroOrderIma);
      cpl_free(outputImages);
      return NULL;
    }

    outputImages[0] = firstOrderIma;
    outputImages[1] = zeroOrderIma;
  }
  else
  {
    outputImages[0] = duplicateImage(flatImages[0]);
    copyAllDescriptors(flatImages[0]->descs, &(outputImages[0])->descs);    
    outputImages[1] = NULL;
  }

  return outputImages;
}


VimosImage *
VmSpApplyFF(VimosImage *inputImage, VimosImage *flatImage,
            VimosExtractionTable *extrTable)
{
  VimosImage   *normImage;
  VimosUlong32 numPoints;
  VimosUlong32 index;
  VimosExtractionSlit *slit;
  VimosDpoint         *nearPoints = NULL;

  float   xZero,yZero;
  float   zeroWidth;
  double  dataVal;
  double *curve = NULL;
  int     imageXlen, imageYlen;
  int     numRows;
  int     doZero;
  int     replaceWidth;
  int     xOut, yOut;
  int     i, j, k;
  int     interpOrder = 2;
  char    modName[] = "VmSpApplyFF";

  cpl_msg_debug(modName,"Applying Flat Field");

  imageXlen = inputImage->xlen;
  imageYlen = inputImage->ylen;
  numPoints = imageXlen*imageYlen;
   
  if (readIntDescriptor(extrTable->descs, pilTrnGetKeyword("ZeroOrderFlag"), 
                        &doZero, NULL) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot find descriptor %s", 
                pilTrnGetKeyword("ZeroOrderFlag"));
    return NULL;
  }

  /* interpolation to remove the zero order contamination (if necessary) */
  if (doZero)
  {
    if (readFloatDescriptor(extrTable->descs,
                            pilTrnGetKeyword("ZeroOrderWidth"),
                            &zeroWidth, NULL) == VM_FALSE)
    {
      cpl_msg_error(modName, "Cannot find descriptor %s", 
		  pilTrnGetKeyword("ZeroOrderWidth"));
      return NULL;
    }

    replaceWidth = ceil(zeroWidth) + 1;
    slit = extrTable->slits;
    nearPoints = newDpoint(2 * replaceWidth);

    /* loop over all  slits */
    while (slit) {
      /* spatial length of slit */
      numRows = slit->numRows;
      for (i = 0; i < numRows; i++ ) 
      {
	yZero = slit->ccdY->data[i] + slit->zeroY->data[i];
	xZero = slit->ccdX->data[i] + slit->zeroX->data[i];
        if (yZero<0 || yZero>=imageYlen || xZero<0 || xZero>=imageXlen) {
          continue;
        }

	/* build up the array of points below and above the zero order
	   contamination, needed to carry out the interpolation */
	xOut = xZero;
	/* the pixels below the contaminated area */
	for (j = 0; j < replaceWidth; j++)
	{
	  yOut = (int) yZero - (int) (1.5 * replaceWidth) + j;
	  index = xOut+yOut*imageXlen;
	  dataVal = inputImage->data[index];
	  /* and store in buffer */
	  nearPoints[j].x = yOut;
	  nearPoints[j].y = dataVal;
	}
	/* and the pixels above the contaminated area */
	for (j = 0; j < replaceWidth; j++)
	{
	  yOut = (int) yZero + ceil(0.5 * replaceWidth) + j;
	  index = xOut+yOut*imageXlen;
	  dataVal = inputImage->data[index];
	  /* and store in buffer */
	  nearPoints[j+replaceWidth].x = yOut;
	  nearPoints[j+replaceWidth].y = dataVal;
	}
	/* do a poly fit to the data */
	curve = fit1DPoly(interpOrder, nearPoints, 2*replaceWidth, NULL);
	if (curve == NULL) {
	  return NULL;
	}
	/* and compute the interpolated values in the contaminated area */
	for (j = 0; j < replaceWidth; j++)
	{
	  yOut = (int) yZero - (int)(0.5 * replaceWidth) + j;
	  index = xOut+yOut*imageXlen;
	  dataVal = 0.;
	  for (k = 0; k <= interpOrder; k++) {
	    dataVal += ipow((double)yOut,k) * curve[k];
	  }
	  inputImage->data[index] = dataVal;
	}
      }
      /* go to next slit */
      slit = slit->next;    
    }
  }

  normImage = newImageAndAlloc(imageXlen, imageYlen);
  
  for (index = 0; index < numPoints; index++) {
    if (flatImage->data[index] != 0.) 
    {
      normImage->data[index] = inputImage->data[index] / 
	  flatImage->data[index];
    }
    else
    {
      normImage->data[index] = inputImage->data[index];
    }
  } 

  copyAllDescriptors(inputImage->descs, &(normImage->descs));
 
  return normImage;
}
VimosBool determineExposedMosArea(VimosTable *adf, 
          VimosExtractionTable *extTable, int *nSlits,
          VimosDpoint **minY, VimosDpoint **maxY, VimosDpoint **zeroY)
{
  const char modName[] = "determineExposedMosArea";
  int        i,j,k;
  int        nPixBelow;
  int        nPixAbove;
  int        crvOrder;
  float      xMin;
  float      xMax;
  float      yMin;
  float      yMax;
  float      yZero;
  float      delXmin = 0.;
  float      delXmax = 0.;
  float      delXlo;
  float      delXhi;
  char       comment[80];

  VimosExtractionSlit *exSlit;
  VimosExtractionSlit *expoSlit;
  VimosBool firstTime;

  /* validate input */
  if (extTable == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return (VM_FALSE);
  }

  exSlit = extTable->slits;

  expoSlit = determineExposedSlits(adf, &(exSlit), nSlits);
  if (expoSlit == NULL) {
    cpl_msg_error(modName, "Function determineExposedSlits failure");
    return(VM_FALSE);
  }

  extTable->slits = exSlit;
  *minY = newDpoint(*nSlits);
  if (*minY == NULL) {
    cpl_msg_error(modName, "Function newDpoint failure");
    return(VM_FALSE);
  }
  *maxY = newDpoint(*nSlits);
  if (*maxY == NULL) {
    cpl_msg_error(modName, "Function newDpoint failure");
    return(VM_FALSE);
  }
  *zeroY = newDpoint(*nSlits);
  if (*zeroY == NULL) {
    cpl_msg_error(modName, "Function newDpoint failure");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(extTable->descs, "ESO PRO SPECT LLEN LO", &nPixBelow,
        comment)) {
    cpl_msg_error(modName, "Cannot read descriptor %s", "ESO PRO SPECT LLEN LO");
    return(VM_FALSE);
  }
  if (!readIntDescriptor(extTable->descs, "ESO PRO SPECT LLEN HI", &nPixAbove,
        comment)) {
    cpl_msg_error(modName, "Cannot read descriptor %s", "ESO PRO SPECT LLEN HI");
    return(VM_FALSE);
  }

  firstTime = VM_TRUE;

  while (expoSlit) {

   /*
    *  Smallest rectangle containing the slit (that may be oblique).
    *  Note that this method doesn't necessarily work for curved
    *  slits where the max or min Y may not correspond to the slit
    *  ends. For the same reason, in this case yZero would not be 
    *  the Y coordinate of the slit center. Nevertheless, this is
    *  quite tolerable (C.Izzo)
    */

    xMin = expoSlit->ccdX->data[0];
    xMax = expoSlit->ccdX->data[(expoSlit->numRows)-1];
    yMin = expoSlit->ccdY->data[0];
    yMax = expoSlit->ccdY->data[(expoSlit->numRows)-1];
    yZero = 0.5 * (yMin + yMax);

    crvOrder = expoSlit->crvPol[0]->order;

    for (i = 0; i <= (nPixAbove+nPixBelow); i++) {
      delXlo = 0.;
      delXhi = 0.;
      k = i - nPixBelow;

     /*
      *  Note: in computing the X position due to the curvature, the
      *  offset parameter is not used  (C.Izzo)
      */

      for (j = 0; j <= crvOrder; j++) {
        delXlo += ipow((double) k, j) * 
                  expoSlit->crvPol[0]->coefs[j];
        delXhi += ipow((double) k, j) * 
                  expoSlit->crvPol[(expoSlit->numRows)-1]->coefs[j];
      }
    
      if (i) {
        if (delXlo < delXmin) delXmin = delXlo;
        if (delXhi > delXmax) delXmax = delXhi;
      }
      else {
       /* Initialization */
        delXmin = delXlo;
        delXmax = delXhi;
      }
    }

   /*
    *  Smallest rectangle containg the whole spectrum:
    */

    yMax = MAX(yMin, yMax) + nPixAbove;
    yMin = MIN(yMin, yMax) - nPixBelow;
    xMax += delXmax;
    xMin += delXmin;

   /*
    *  Here yZero is upgraded to be the middle point Y coordinate
    *  of the slit for the first order, plus the middle point Y
    *  diff coordinate of the slit for the zeroth order. In total,
    *  is the Y coordinate of the middle point of the slit for the
    *  zeroth order. Note that it could be set like that since the 
    *  beginning, as it wasn't used in the meantime. (C.Izzo)
    */
    yZero += 0.5 * (expoSlit->zeroY->data[0] + 
             expoSlit->zeroY->data[(expoSlit->numRows)-1]);

    if (firstTime == VM_TRUE) {
      firstTime = VM_FALSE;
    }
    else {
      *minY = (*minY)->next;
      *maxY = (*maxY)->next;
      *zeroY = (*zeroY)->next;
    }
    (*minY)->x = (double) xMin;
    (*minY)->y = (double) yMin;
    (*maxY)->x = (double) xMax;
    (*maxY)->y = (double) yMax;
    (*zeroY)->y = (double) yZero;

   /*
    *  Note that in this way the pointer is not set back to the
    *  beginning of the list. The last values of *minY, *maxY
    *  and *zeroY are returned, being passed by address. This 
    *  must be taken into account in the caller (C.Izzo)
    */

    expoSlit = expoSlit->next;
  }

  return (VM_TRUE);
}

VimosExtractionSlit *determineExposedSlits(VimosTable *adf, 
                     VimosExtractionSlit **exSlit, int *nSlits)
{
  const char           modName[] = "determineExposedSlits";
  char                 comment[80];
  int                  quadNum;
  float                yMin;
  float                yMax;
  float                refLow;
  float                refHigh;
  float                slYmin;
  float                slYmax;
  VimosExtractionSlit *lastSlit;
  VimosExtractionSlit *lastLostSlit;
  VimosExtractionSlit *lostSlit;
  VimosExtractionSlit *returnSlit;
  VimosExtractionSlit *t1Slit;
  VimosExtractionSlit *t2Slit;
  VimosBool            firstTime, foundNoIllum, firstNoIllum;


  /* 
   *  find out the minimum and maximum y mask coordinate in the 
   *  illuminated area 
   */

  if (readIntDescriptor(adf->descs, pilTrnGetKeyword("Quadrant"), 
                         &quadNum, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("Quadrant"));
    return NULL;
  }

  if (readFloatDescriptor(adf->descs, pilTrnGetKeyword("MshuPosH", quadNum), 
                           &yMax, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("MshuPosH", quadNum));
    return NULL;
  }

  if (readFloatDescriptor(adf->descs, pilTrnGetKeyword("MshuPosL", quadNum), 
                          &yMin, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("MshuPosL", quadNum));
    return NULL;
  }

  if (readFloatDescriptor(adf->descs, pilTrnGetKeyword("MshuRefH", quadNum), 
                           &refHigh, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("MshuRefH", quadNum));
    return NULL;
  }

  if (readFloatDescriptor(adf->descs, pilTrnGetKeyword("MshuRefL", quadNum), 
                          &refLow, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("MshuRefL", quadNum));
    return NULL;
  }

  /* correct the position of the shutters as given by the ICS for the 
     positioning offset between VMMCS and ICS */
  yMin = yMin - refLow;
  yMax = refHigh - yMax;

  *nSlits = 0;
  lastSlit = NULL;
  lostSlit = NULL;
  t1Slit = NULL;
  firstTime = VM_TRUE;

  if (!slitMinMaxY(*exSlit, &slYmin, &slYmax)) {
    cpl_msg_error(modName, "Function slitMinMaxY failure");
    return(NULL);
  }

  /* FIXME:
   *  Note that this check may lead to "lost slits" for
   *  any one of the shutter positions: a slit at the limit
   *  of the windows may have a minimum inside a window,
   *  and a maximum inside another. This slit would be
   *  rejected in both cases (C.Izzo)
   */

  /* determination of the new entry point in the slits list */
  /* skip those slits located outside the shutter window */
  while (slYmin < yMin || slYmin > yMax || slYmax < yMin || slYmax > yMax) {
    if (firstTime) {
     /* FIXME:
      *  As soon as the "while" is entered (and this can be only
      *  if the very first slit is outside the shutter window),
      *  the list of rejected lists (lostSlit) is initialized
      *  to the input slit list. The temporary pointer t1Slit
      *  is positioned to the next element, which is also checked
      *  and we stay in the "while" as long as slits are outside 
      *  the shutter window. (C.Izzo)
      */
      lostSlit = *exSlit;
 /* t1Slit = newExtractionSlit(); This is unnecessary (memory leak) (C.Izzo) */
      t1Slit = (*exSlit)->next;
      firstTime = VM_FALSE;
    }
    else {
      t1Slit = t1Slit->next;
    }
    if (t1Slit == NULL) {
     /* FIXME:
      *  I believe that this error message is inappropriate:
      *  we reach this point as soon as all slits are rejected,
      *  that is, no slits are found within the shutter window.
      *  This is an acceptable situation. A NULL is correctly
      *  returned, but it is the higher level function that 
      *  should decide how to deal with this (legal) case. (C.Izzo)
      */
      /* cpl_msg_error(modName, "NULL next pointer in the structure"); */
      cpl_msg_error(modName, "No slits are found within this shutter window");
      return (NULL);
    }
    if (!slitMinMaxY(t1Slit, &slYmin, &slYmax)) {
      cpl_msg_error(modName, "Function slitMinMaxY failure");
      return(NULL);
    }
  }

  /*
   *  Here t1Slit points to the first slit found within the shutter
   *  window. The only case in which it might be t1Slit = NULL is 
   *  that the "while" was never entered (i.e., the first slit is
   *  inside the shutter window) (C.Izzo)
   */

  if (t1Slit) {
    lastSlit = t1Slit;
    returnSlit = t1Slit;
    lastLostSlit = t1Slit->prev;
    lastLostSlit->next = NULL; /* To set lastLostSlit->next = NULL is */
    t1Slit->prev->next = NULL; /* equivalent to set t1Slit->prev->next = NULL,
                                  therefore one of the 2 lines is superfluous.
                                  (C.Izzo) */

   /* NOTE: The *exSlit list is not duplicated when we set t1Slit
      to point to one element of the slit list. Therefore, to set 
      t1Slit->prev->next = NULL is equivalent to cut the original
      list at this point. We have then a list of rejected slits 
      starting at *exSlit, and ending with the first accepted slit 
      (excluded).  (C.Izzo) */
  }
  else {
   /*
    *  The "while" was never entered (i.e., the first slit is
    *  inside the shutter window) (C.Izzo)
    */
    lastSlit = *exSlit;
    returnSlit = *exSlit;
   /* FIXME: 
    * lastLostSlit in this case is NOT initialized. It contains
    * garbage. This might be a bug to solve (C.Izzo) 
    */
  }    

 /*
  *  Also returnSlit points currently to the first illuminated slit.
  *  To complete the separation of the list of illuminated and not
  *  illuminated slits, its previous element must be set to NULL.
  *  (C.Izzo)
  */
  returnSlit->prev = NULL;
  firstTime = VM_TRUE;
  foundNoIllum = VM_FALSE;

  /* determination of the remaining exposed slits */
  if (lastSlit->next != NULL) {   /* Does another slit exists, after the */
    do {                          /* last illuminated one?  (C.Izzo)     */
 /* t2Slit = newExtractionSlit(); This is unnecessary (memory leak) (C.Izzo)
      if (t2Slit == NULL) {
        cpl_msg_error(modName, "Function newExtractionSlit failure");
        return(NULL);
      }  ***/
      t2Slit = lastSlit->next;   /* Points to next slit to check  (C.Izzo) */
      if (!slitMinMaxY(t2Slit, &slYmin, &slYmax)) {
        cpl_msg_error(modName, "Function slitMinMaxY failure");
        return(NULL);
      }
      firstNoIllum = VM_TRUE;

  /* skip those slits located outside the shutter window */
      while (slYmin < yMin || slYmin > yMax || slYmax < yMin || slYmax > yMax) {
        foundNoIllum = VM_TRUE;
       /*
        *  lostSlit is NULL only if no slit was lost at the
        *  beginning of the list  (C.Izzo) 
        */
        if (!lostSlit) {
         /*
          *  If it is so, then we have found the beginning of the list
          *  of lost slits. Then, cut it away from the original list  
          *  (C.Izzo)
          */
          lostSlit = t2Slit;
          lostSlit->prev = NULL;
          firstNoIllum = VM_FALSE;
        }
        else if (firstNoIllum) {
         /*
          *  Instead, if this is just another rejected slit, and it is
          *  the first one to be found in this round, just append it 
          *  to the list of lost slits. This was surely initialized,
          *  because we get here only if some slits were rejected
          *  at the beginning of the list (C.Izzo)
          */
          lastLostSlit->next = t2Slit;
          t2Slit->prev = lastLostSlit;
          firstNoIllum = VM_FALSE;
        }

        t2Slit = t2Slit->next;
        if (t2Slit == NULL) {
          /* close linked list */
          lastSlit->next = NULL;
         /*
          *  *exSlit, passed by address, is also returned: it contains
          *  now just the rejected slits (the only necessary for further 
          *  checks whith different shutter positions). But note: now
          *  the original list is destroyed (C.Izzo).
          */
          *exSlit = lostSlit;
          return returnSlit;
        }
        else {
          if (!slitMinMaxY(t2Slit, &slYmin, &slYmax)) {
            cpl_msg_error(modName, "Function slitMinMaxY failure");
            return(NULL);
          }
        }
      }

      if (foundNoIllum) {
       /*
        *  If the last "while" was entered, upgrade the position
        *  of the last lost slit. (C.Izzo)
        */
        lastLostSlit = t2Slit->prev;
        lastLostSlit->next = NULL;
        t2Slit->prev->next = NULL; /* Superfluous, see note above (C.Izzo) */
        foundNoIllum = VM_FALSE;
      }

      t2Slit->prev = lastSlit;
      if (firstTime) {
       /*
        *  Disconnect the still-to-be-checked slit list from the
        *  original list, and connect it to the list of rejected
        *  (C.Izzo)
        */
        t2Slit->prev = returnSlit;
        returnSlit->next = t2Slit;
       /*
        *  Initialization of the returned number of illuminated
        *  slits. Every time we get here, we terminated an interval
        *  of rejected slits. The first round was at the beginning
        *  of the program, there we have found the first illuminated
        *  slit. If we get here we have found another illuminated
        *  slit. That makes 2  (C.Izzo).
        */
        (*nSlits) = 2;
        firstTime = VM_FALSE;
      }
      else {
       /*
        *  Every time we get here, we terminated an interval
        *  of rejected slits. One illuminated slit (the one that
        *  interrupted the scan) is counted (C.Izzo)
        */
        lastSlit->next = t2Slit;
        (*nSlits)++;
      }
      lastSlit = t2Slit;

    } while (t2Slit->next != NULL);
  }
  else {
    returnSlit->next = NULL;
    (*nSlits) = 1;
  }
  *exSlit = lostSlit;
  return returnSlit;
}

VimosBool slitMinMaxY(VimosExtractionSlit *exSlit, float *slYmin, float *slYmax)
{
  float min, max;
  int i;

  min = max = exSlit->maskY->data[0];

  for (i = 1; i < exSlit->numRows; i++) {
    if (exSlit->maskY->data[i] < min) min = exSlit->maskY->data[i];
    if (exSlit->maskY->data[i] > max) max = exSlit->maskY->data[i];
  }
  
  *slYmax = max;
  *slYmin = min;
  return (VM_TRUE);
}
/**@}*/

