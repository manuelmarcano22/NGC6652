/* $Id: vmmosmodels.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include "vmfit.h"
#include "vmmath.h"
#include "vmmosutils.h"
#include "vmifuutils.h"
#include "vmmosmodels.h"
#include "cpl.h"


/**
 * @brief vmmosmodels MOS Optical Models
 *
 * The module provides functions to model the instrument optical distortions
 * and the spectral curvature.
 */

/**@{*/

int
VmSpDerCurves(VimosImage **imageData, VimosExtractionTable *extractionTable,
              int fuzz)
{
  char                 modName[] = "VmSpDerCurves";
  int                  i, j, k;
  int                  outOfImage = 0;
  int                  xOut, yOut;
  int                  numRows;
  int                  numPoints;
  int                  numFitPointsLower, numFitPointsUpper;
  VimosUlong32         index;
  int                  nPixBelow, nPixAbove;
  int                  imageXlen, imageYlen;
  int                  crvOrder; 

  double               datVal;
  double               xOutF,yOutF;
  double               frac;
  double               fitVariance;
  double              *curve = NULL;
  float                xOffset;
  float                xStart, xEnd;
  float                xZero,yZero;
  double               zeroUpper,zeroLower;
  float                zeroWidth;
  int                  extraWidth = 10;
  
  int                  doZero;
  int mycount = 0;
  VimosImage          *firstOrderIma;
  VimosImage          *zeroOrderIma;
  VimosDpoint         *lower = NULL;
  VimosDpoint         *upper = NULL;
  VimosFloatArray     *tmpArray;
  VimosExtractionSlit *slit;
  
  /*ALEX */
  char insMode[80];

  cpl_msg_info(modName, "Deriving Curvature Polynomials");

  /*ALEX: determine numSlit depending on MOS or IFU: if IFU, numSlit equal to
    number of good fibers (to prevent zeros in medianWirth on offsetYArray */
  readStringDescriptor(imageData[0]->descs,
                       pilTrnGetKeyword("InstrumentMode"), 
		       insMode, NULL);

  firstOrderIma = newImage(imageData[0]->xlen,imageData[0]->ylen,
                           imageData[0]->data);
  firstOrderIma->descs = imageData[0]->descs;
  imageXlen = firstOrderIma->xlen;
  imageYlen = firstOrderIma->ylen;
  xStart = 0.;

  if (readIntDescriptor(extractionTable->descs,
                        pilTrnGetKeyword("ZeroOrderFlag"), 
                        &doZero, NULL) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot find descriptor %s",
                pilTrnGetKeyword("ZeroOrderFlag"));
    return EXIT_FAILURE;
  }

  if (doZero) {
    zeroOrderIma = newImage(imageData[1]->xlen,imageData[1]->ylen,
                            imageData[1]->data);
    zeroOrderIma->descs = imageData[1]->descs;
  }

  /* read the order of the Curvature Polynomials in the header of the
     Extraction Table */
  readIntDescriptor(extractionTable->descs, "ESO PRO CRV POL ORD",
                     &crvOrder, NULL);

  /* read length of spectrum in pixels in wavelength direction */
  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN LO", 
                    &nPixBelow, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN HI",
                    &nPixAbove, NULL);
  /* length of spectrum inpixels in the wavelength direction */
  numPoints = nPixBelow + nPixAbove + 1;
  
  /* create buffers to store slit edges as function of wavelength */
  upper = newDpoint(numPoints);
  lower = newDpoint(numPoints);
  
  /* first slit in table */
  slit = extractionTable->slits;
  /* loop over all  slits */
  while (slit) {

    /* ALEX: when IFU, do nothing for dead fibers */
    if (slit->IFUfibTrans >= 0.0)
      {
	/* spatial length of slit */
	numRows = slit->numRows;
    
	/* create temp buffer */
	tmpArray = newFloatArray(numRows + 2*fuzz);  
	numFitPointsLower = 0;
	numFitPointsUpper = 0;

	/* loop over 'wavelength' */
	for (j = -nPixBelow; j <= nPixAbove; j++) {
	  numPoints = 0;

	  /* extract a cut across the slit, with a bit of extra space */
	  for (i = -fuzz; i < numRows+fuzz; i++) {
	    if (i < 0)
	      {
		k = 0;
		yOut = slit->ccdY->data[k] + j;
		xOutF = slit->ccdX->data[k] + 
		  computeDistModel1D(slit->crvPol[k], yOut) + i;
	      }
	    else if (i >= numRows)
	      {
		k = numRows - 1;
		yOut = slit->ccdY->data[k] + j;
		xOutF = slit->ccdX->data[k] + 
		  computeDistModel1D(slit->crvPol[k], yOut) + i - (numRows-1);
	      }
	    else
	      {
		k = i;
		yOut = slit->ccdY->data[k] + j;
		xOutF = slit->ccdX->data[k] + 
		  computeDistModel1D(slit->crvPol[k], yOut);
	      }
	    /* use curvature polynomial of edge to tracte these  data */
	    xOut = xOutF;

            /*
             * 100 = safety marging from X image borders.
             */

	    if (yOut<0 || yOut >= imageYlen || 
                xOut < 100 || xOut >= imageXlen - 100) {
	      outOfImage = 1; 
	      continue;
	    } 
	    outOfImage = 0;
	    numPoints++;

	    /* quick interpolation of data */
	    frac = xOutF-xOut;
	    index = xOut+yOut*imageXlen;
	    datVal = ( (1.0-frac)*firstOrderIma->data[index] +
		       frac*firstOrderIma->data[index+1] );
	    /* and store in buffer */
	    tmpArray->data[i+fuzz] = datVal;

	    if (i == 0)
	      xStart = slit->ccdX->data[0] + 
		computeDistModel1D(slit->crvPol[0],yOut);
	    if (i == numRows-1)
	      xEnd = slit->ccdX->data[numRows-1] +
		computeDistModel1D(slit->crvPol[numRows-1],yOut);
	  }
      
	  if (numPoints < numRows) {
	    continue;
	  }
	  else {
	    if (!strncmp(insMode,"IFU",3)) {
	      if (!(findIfuBorders(tmpArray, &(upper[numFitPointsUpper].y),
				   &(lower[numFitPointsLower].y)))) {

/*		cpl_msg_error(modName, 
		"function findIfuBorders returned error"); */
/*    		cpl_msg_warning(modName, 
		"%d %d %f %d \n",slit->IFUslitNo,slit->IFUfibNo,xStart,yOut);*/

		mycount++;

	      } else {
		upper[numFitPointsUpper].y += xStart-fuzz;
		lower[numFitPointsLower].y += xStart-fuzz;
		/* write pixel index as 'X' coordinate in buffer */
		upper[numFitPointsUpper].x =  j;
		lower[numFitPointsLower].x =  j;
		
		numFitPointsUpper++;
		numFitPointsLower++;
	      }
	    } else {  /* MOS */
	      findSpectrumBorders(tmpArray, &(upper[numFitPointsUpper].y),
				  &(lower[numFitPointsLower].y), fuzz); 
	      /* add the starting offset to the 'Y' coordinate in buffer */
	      /* NOTE: we add 1 to the lower edge and subtract 1 to the upper
		 edge to reduce the problems with handling of edge effects */

              if (upper[numFitPointsUpper].y > 0.0) {
	        upper[numFitPointsUpper].y += xStart-fuzz;
	        upper[numFitPointsUpper].x =  j;
	        numFitPointsUpper++;
              }

              if (lower[numFitPointsLower].y > 0.0) {
	        lower[numFitPointsLower].y += xStart-fuzz;
	        lower[numFitPointsLower].x =  j;
	        numFitPointsLower++;
              }
	    }
	  }
	}

	if (numFitPointsLower == 0 || numFitPointsUpper == 0) {
	  outOfImage = 0; 
	  slit = slit->next;
	  continue;
	}

        if (fabs(numFitPointsLower - numFitPointsUpper) / numFitPointsUpper 
            > 0.25 ) {
          outOfImage = 0;
          slit = slit->next;
          continue;
        }

	/* delete array with coefficients, if it exists */
	if (curve!= NULL) {
	  cpl_free(curve);
	  curve = NULL;
	}
	/* delete the existing curvature polynomials in slit, in case the order
	   is not the same as we are fitting */
	for (i = 0; i < numRows; i++ ) {
	  deleteDistModel1D(slit->crvPol[i]);
	  /* and create of the right order */
	  slit->crvPol[i] = newDistModel1D(crvOrder);
	}
	/* fit polynomial to lower edge of slit */
	curve = fit1DPoly(crvOrder, lower, numFitPointsLower, &fitVariance);
	if (curve == NULL) {
	  printf("VmSpDerCurves: Invalid fit\n");
	  /* no fit, abort */
	  return EXIT_FAILURE;
	}
	/* copy coefficients to curvature polynomial */
	for (i = 0; i <= crvOrder; i++) {
	  slit->crvPol[0]->coefs[i] = curve[i];
	}

	slit->crvPolRms->data[0] = sqrt(fitVariance);

	/* set the right offset for the curvature polynomial */
	slit->crvPol[0]->offset = slit->ccdY->data[0];
	/* compute where lower edge for central wavelength is according to 
	   model */

	xOffset =  slit->crvPol[0]->coefs[0] - slit->ccdX->data[0];
	/* and correct models for this offset, so we can fit both a global
	   curvature model AND an optical distortion model (at least in X) to
	   the results of this function */    
	slit->crvPol[0]->coefs[0] = 0.;
	for (i = 0; i < numRows; i++)
	  slit->ccdX->data[i] += (int) xOffset;
    
	/* free buffer */
	cpl_free(curve);
	curve = NULL;
    
	/* same story for upper edge of slit */
	curve = fit1DPoly(crvOrder, upper, numFitPointsUpper, &fitVariance);
	if (curve == NULL) {
	  return EXIT_FAILURE;
	}
	for (i = 0; i <= crvOrder; i++) {
	  slit->crvPol[numRows-1]->coefs[i] = curve[i];
	}
	slit->crvPolRms->data[numRows-1] = sqrt(fitVariance);
	slit->crvPol[numRows-1]->offset = slit->ccdY->data[numRows-1];
	slit->crvPol[numRows-1]->coefs[0] = 0.;
  
	/* destruct temp */
	deleteFloatArray(tmpArray);

	/* This is a new pice of code for the fitting of the Zero order
	   spectrum position */

	if (doZero) {
	  if (readFloatDescriptor(extractionTable->descs,
				  pilTrnGetKeyword("ZeroOrderWidth"), 
				  &zeroWidth, NULL) == VM_FALSE)
	    {
	      cpl_msg_error(modName, "Cannot find descriptor %s", 
			  pilTrnGetKeyword("ZeroOrderWidth"));
	      return EXIT_FAILURE;
	    }

	  yZero = 0.5 * (slit->ccdY->data[0]+slit->zeroY->data[0] + 
			 slit->ccdY->data[numRows-1] + 
			 slit->zeroY->data[numRows-1]);
	  xZero = 0.5 * (slit->ccdX->data[0]+slit->zeroX->data[0] + 
			 slit->ccdX->data[numRows-1] + 
			 slit->zeroX->data[numRows-1]);
	  tmpArray = newFloatArray((int) zeroWidth + 2 * extraWidth);
	  xOut = xZero;
	  for (i = 0; i < (int) zeroWidth + 2 * extraWidth; i++)
	    {
	      yOutF = yZero - 0.5 * (float) (zeroWidth + 2 * extraWidth) + 
		(float) i;
	      yOut = yOutF;
	      if (yOut < 0 || yOut >= imageYlen) {
		outOfImage = 1; 
		break;
	      }
	      frac = yOutF - yOut;
	      index = xOut+yOut*imageXlen;
	      datVal = ( (1.0-frac)*zeroOrderIma->data[index] +
			 frac*zeroOrderIma->data[index+1] );
	      /* and store in buffer */
	      tmpArray->data[i] = datVal;
	    }
	  if (outOfImage)
	    outOfImage = 0;
	  else 
	    findSpectrumBorders(tmpArray, &zeroUpper, &zeroLower, extraWidth); 
    
	  deleteFloatArray(tmpArray);
	  tmpArray = newFloatArray(numRows + 2*fuzz);  

	  yZero = 0.5 * (zeroLower + zeroUpper);
	  yOut = yZero;
	  for (i = -fuzz; i < numRows+fuzz; i++) {
	    if (i < 0)
	      k = 0;
	    else if (i >= numRows)
	      k = numRows - 1;
	    else
	      k = i;
	    xOutF = slit->ccdX->data[k] + slit->zeroX->data[k] + i;
	    xOut = xOutF;
	    /* quick interpolation of data */
	    frac = xOutF-xOut;
	    index = xOut+yOut*imageXlen;
	    datVal = ( (1.0-frac)*zeroOrderIma->data[index] +
		       frac*zeroOrderIma->data[index+1] );
	    /* and store in buffer */
	    tmpArray->data[i+fuzz] = datVal;
	    if (i == 0)
	      xStart = slit->ccdX->data[0] + slit->zeroX->data[0];
	    if (i == numRows-1)
	      xEnd= slit->ccdX->data[numRows-1] + slit->zeroX->data[numRows-1];
	  }
	  findSpectrumBorders(tmpArray, &zeroUpper, &zeroLower, fuzz); 
	  zeroLower += xStart-fuzz;
	  zeroUpper += xStart-fuzz;

	  slit->zeroX->data[0] = zeroLower - slit->ccdX->data[0];
	  slit->zeroY->data[0] = yZero - slit->ccdY->data[0];
	  slit->zeroX->data[numRows-1]=zeroUpper - slit->ccdX->data[numRows-1];
	  slit->zeroY->data[numRows-1] = yZero - slit->ccdY->data[numRows-1];
	}
	/* end of the new part dealing with zero order stuff */

      } /* ALEX: end if on dead fibers */

    /* go to next slit */
    slit = slit->next;    
  }
  
  /* clean up... */
  if (curve != NULL) {
    cpl_free(curve);
  }
  deleteDpoint(upper);
  deleteDpoint(lower);
  
  return EXIT_SUCCESS;
}


/**
 * @memo
 *   Fit a Curvature Distorsion Model
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 * 
 * @param extractionTable Input Extraction Table
 * @param grismTable      Input (and perhaps output) Grism Table
 * @param grismFlag       Flag to indicate whether to update or not
 *                        the Grism Table with the new Curvature Model
 *
 * @doc
 *   Just fit a polynomial in two variables to the coefficients
 *   of the polynomials fitting the edges of spectra found in a 
 *   spectral Flat Field
 *
 * @author M. Scodeggio
 */

int
VmSpCurveModel(VimosExtractionTable *extractionTable, VimosTable *grismTable,
               int grismFlag)
{
  int i; 
  int numSlits;
  int slitDone;
  int polOrder, xOrder, yOrder;
  double rmsFit;
  
  VimosExtractionSlit *slit;
  VimosPixel *coefArray;
  VimosDistModelFull *crvMod;
  
  char modName[]   ="VmSpCurveModel";

  
  cpl_msg_debug(modName, "Compute curvature model");

  numSlits = 0;
  slit = extractionTable->slits;

  while (slit) {
    numSlits++;
    slit = slit->next;
  }
  
  coefArray = newPixel(2*numSlits);
  readIntDescriptor(extractionTable->descs, "ESO PRO CRV POL ORD", 
                    &polOrder, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO CRV MOD XORD", 
                    &xOrder, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO CRV MOD YORD", 
                    &yOrder, NULL);
    
  /* note that for the moment we keep the order in X equeal to the one in Y!!*/
  crvMod = newDistModelFull(polOrder, xOrder, xOrder);
  
  for (i = 0; i <= polOrder; i++) {
    slit = extractionTable->slits;
    slitDone = 0;
    while (slit) {
      coefArray[slitDone].x = slit->maskX->data[0];
      coefArray[slitDone].y = slit->maskY->data[0];
      coefArray[slitDone].i = slit->crvPol[0]->coefs[i];
      slitDone++;
      coefArray[slitDone].x = slit->maskX->data[slit->numRows-1];
      coefArray[slitDone].y = slit->maskY->data[slit->numRows-1];
      coefArray[slitDone].i = slit->crvPol[slit->numRows-1]->coefs[i];
      slitDone++;
      slit = slit->next;
    }
    /* ugly!!!! */
    deleteDistModel2D(crvMod->coefs[i]);
      
    if (fitDistModel2D(coefArray, 2*numSlits, xOrder, 0.0, 0.0,
                       &(crvMod->coefs[i]), &rmsFit) == VM_FALSE)
      return EXIT_FAILURE;
  }

  writeCurvatureModel(&(extractionTable->descs), crvMod);

  if (grismFlag)
    writeCurvatureModel(&(grismTable->descs), crvMod);
  
  deletePixel(coefArray);
  deleteDistModelFull(crvMod);
  
  return EXIT_SUCCESS;
}
/**
 * @memo
 *   Fit an Optical Distorsion Model
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param extractionTable Input Extraction Table
 * @param grismTable      Input (and perhaps output) Grism Table
 * @param grismFlag       Flag to indicate whether to update or
 *                        not the Grism Table with the new Optical 
 *                        Distortion Model.
 * @doc
 *   [MISSING]
 *
 * @author M. Scodeggio
 */

int VmSpOptModel
(VimosExtractionTable *extractionTable, VimosTable *grismTable, int grismFlag)
{
  int numSlits;
  int slitDone;
  int order;
  double rmsFit;
  
  VimosExtractionSlit *slit;
  VimosPixel *coefArray;
  VimosDistModel2D *optModX;
  VimosDistModel2D *optModY;
  int               doZero;

  char modName[]   ="VmSpOptModel";

  cpl_msg_debug(modName, "Fit curvature model");

  numSlits = 0;
  slit = extractionTable->slits;

  while (slit) {
    numSlits++;
    slit = slit->next;
  }
  
  coefArray = newPixel(2*numSlits);
  
  slit = extractionTable->slits;
  slitDone = 0;
  while (slit) {
    coefArray[slitDone].x = slit->maskX->data[0];
    coefArray[slitDone].y = slit->maskY->data[0];
    coefArray[slitDone].i = slit->ccdX->data[0];
    slitDone++;
    coefArray[slitDone].x = slit->maskX->data[slit->numRows-1];
    coefArray[slitDone].y = slit->maskY->data[slit->numRows-1];
    coefArray[slitDone].i = slit->ccdX->data[slit->numRows-1];
    slitDone++;
    slit = slit->next;
  }

  readIntDescriptor(extractionTable->descs, "ESO PRO OPT DIS XORD", 
                    &order, NULL);

  if (fitDistModel2D(coefArray, 2*numSlits, order, 
                     0.0, 0.0, &optModX, &rmsFit) == VM_FALSE)
    return EXIT_FAILURE;

  slit = extractionTable->slits;
  slitDone = 0;
  while (slit) {
    coefArray[slitDone].x = slit->maskX->data[0];
    coefArray[slitDone].y = slit->maskY->data[0];
    coefArray[slitDone].i = slit->ccdY->data[0];
    slitDone++;
    coefArray[slitDone].x = slit->maskX->data[slit->numRows-1];
    coefArray[slitDone].y = slit->maskY->data[slit->numRows-1];
    coefArray[slitDone].i = slit->ccdY->data[slit->numRows-1];
    slitDone++;
    slit = slit->next;
  }

  readIntDescriptor(extractionTable->descs, "ESO PRO OPT DIS YORD", 
                    &order, NULL);
 
  if (fitDistModel2D(coefArray, 2*numSlits, order, 
                     0.0, 0.0, &optModY, &rmsFit) == VM_FALSE)
    return EXIT_FAILURE;

  writeOptDistModel(&(extractionTable->descs), optModX, optModY);

  if (grismFlag) {
    writeOptDistModel(&(grismTable->descs), optModX, optModY);
  }
  
  if (readIntDescriptor(extractionTable->descs, 
                        pilTrnGetKeyword("ZeroOrderFlag"),
                        &doZero, NULL) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot find descriptor %s",
                pilTrnGetKeyword("ZeroOrderFlag"));
    return EXIT_FAILURE;
  }

  if (doZero)
  {
    slit = extractionTable->slits;
    slitDone = 0;
    while (slit) {
      coefArray[slitDone].x = slit->maskX->data[0];
      coefArray[slitDone].y = slit->maskY->data[0];
      coefArray[slitDone].i = slit->zeroX->data[0];
      slitDone++;
      coefArray[slitDone].x = slit->maskX->data[slit->numRows-1];
      coefArray[slitDone].y = slit->maskY->data[slit->numRows-1];
      coefArray[slitDone].i = slit->zeroX->data[slit->numRows-1];
      slitDone++;
      slit = slit->next;
    }

    readIntDescriptor(extractionTable->descs, "ESO PRO ZERO XORD", 
		      &order, NULL);

    if (fitDistModel2D(coefArray, 2*numSlits, order, 
                       0.0, 0.0, &optModX, &rmsFit) == VM_FALSE)
      return EXIT_FAILURE;

    slit = extractionTable->slits;
    slitDone = 0;
    while (slit) {
      coefArray[slitDone].x = slit->maskX->data[0];
      coefArray[slitDone].y = slit->maskY->data[0];
      coefArray[slitDone].i = slit->zeroY->data[0];
      slitDone++;
      coefArray[slitDone].x = slit->maskX->data[slit->numRows-1];
      coefArray[slitDone].y = slit->maskY->data[slit->numRows-1];
      coefArray[slitDone].i = slit->zeroY->data[slit->numRows-1];
      slitDone++;
      slit = slit->next;
    }

    readIntDescriptor(extractionTable->descs, "ESO PRO ZERO YORD", 
		      &order, NULL);
 
    if (fitDistModel2D(coefArray, 2*numSlits, order, 
                       0.0, 0.0, &optModY, &rmsFit) == VM_FALSE)
      return EXIT_FAILURE;
    
    writeContaminationModel(&(extractionTable->descs), optModX, optModY);

    if (grismFlag) {
      writeContaminationModel(&(grismTable->descs), optModX, optModY);
    }
  }

  deletePixel(coefArray);
  deleteDistModel2D(optModX);
  deleteDistModel2D(optModY);
  
  return EXIT_SUCCESS;
}
/**@}*/
