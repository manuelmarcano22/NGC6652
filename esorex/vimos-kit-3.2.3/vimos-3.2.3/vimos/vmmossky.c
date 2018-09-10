/* $Id: vmmossky.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <math.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmmatrix.h"
#include "vmwindowtable.h"
#include "vmextractiontable.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmmossky.h"
#include "cpl.h"


#define DO_MEDIAN 1
#define DO_FIT    2


VimosImage **
VmSpSkyFra(VimosImage *inputImage, VimosExtractionTable *extractionTable, 
           VimosWindowTable *windowTable, SkyMethod skyMode, int polyDeg,
           int limMargin, int nSkyMin, int maxIter, float nSigBelow,
           float nSigAbove)
{
  int i, j, k;
  
  int yOut;
  int xLow,xHigh;
  int startX,endX;
  int numGoodRows;
  int numRows;
  VimosUlong32 numPoints;
  unsigned int index;
  int nPixBelow, nPixAbove;
  int imageXlen, imageYlen;
  int nReject,nIter;
  
  float skyVal;
  float xp, var, resid;
  
  double rmsVal,rms;
  double *skyCoefs=NULL;
  
  VimosFloatArray *datVal=NULL;
  VimosImage      *outData;
  VimosImage      *skyData;
  VimosImage     **outputImages;

  VimosDpoint     *datPoint=NULL;
  
  VimosExtractionSlit *extSlit;
  VimosWindowSlit *winSlit;
  
  char modName[]   ="VmSpSkyFra";

  cpl_msg_debug(modName, "Subtracting sky in short slits");

  imageXlen = inputImage->xlen;
  imageYlen = inputImage->ylen;

  outData   = newImageAndAlloc(imageXlen, imageYlen);
  /*  openNewMidasImage(outDataName, outData); */
  
  skyData   = newImageAndAlloc(imageXlen, imageYlen);
  /*  openNewMidasImage(outSkyName, skyData); */
    
  numPoints = imageXlen*imageYlen;
  
  for (index = 0; index < numPoints; index++) {
    outData->data[index] = inputImage->data[index];
    skyData->data[index] = 0.0;
  }
  
  extSlit = extractionTable->slits;
  winSlit = windowTable->slits;

  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN LO", 
                    &nPixBelow, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN HI",
                    &nPixAbove, NULL);

  numPoints = nPixBelow+nPixAbove+1;
  
  while ( (extSlit) && (winSlit) ) {

    if (extSlit->slitNo != winSlit->slitNo) {
      cpl_msg_error(modName, "Extraction and Window tables out of sync !!");
      return NULL;
    }
    
    numRows = extSlit->numRows;

    /* we just copy  long slits */
    if (winSlit->specLong == VM_FALSE) {
      switch (skyMode) {
      default: 
      case SKY_MEDIAN : 
	{
	  datVal = newFloatArray(numRows-limMargin*2);
	  break;
	}
      case SKY_FIT : 
	{
	  datPoint = newDpoint(numRows-limMargin*2);
	  break;
	}
      }
    
      for (j = -nPixBelow; j <= nPixAbove; j++) {
	/* Y coordinate of row. This is constant for the spectrum since 
	   this is a short slit */
	yOut = extSlit->ccdY->data[0] + j;
	/* X coordinate of lower edge of spectrum */
	xLow = extSlit->ccdX->data[0] + 
	  computeDistModel1D(extSlit->crvPol[0], yOut);
	/* X coordinate of higher edge of spectrum */
	xHigh = extSlit->ccdX->data[numRows-1] + 
	  computeDistModel1D(extSlit->crvPol[numRows-1], yOut);

	/* loop over rows to get data for this wavelength */
	if (yOut >= 0 && yOut < imageYlen) {

	  /* determine slit pixels inside the CCD */
	  startX = limMargin;
	  endX = numRows-limMargin;
	  for (i = xLow+limMargin; i < xHigh-limMargin; i++) {
	    if(i<0) startX++;
	    if(i>=imageXlen) endX--;
	  }
	
	  numGoodRows = 0;      
	  for (i = startX; i < endX; i++) {
	    if (windowObjectInRow(winSlit, i) == VM_FALSE) {
	      switch (skyMode) {
	      default: 
	      case SKY_MEDIAN : 
		{
		  datVal->data[numGoodRows] = 
		    inputImage->data[xLow + i + yOut*imageXlen];
		  break;
		}
	      case SKY_FIT : 
		{
		  datPoint[numGoodRows].x = i;
		  datPoint[numGoodRows].y = 
		    inputImage->data[xLow + i + yOut*imageXlen];
      /* printf("%d %f\n",i,inputImage->data[xLow + i + yOut*imageXlen]); */
		  break;
		}
	      }
	      numGoodRows++;
	    }
	  }

	  /* at least 10 pixel to compute sky */
	  if(numGoodRows<10) continue;
      
	  nIter=0;
	  do {
	    nReject=0;
	    switch (skyMode) {
	    default: 
	    case SKY_MEDIAN : 
	      {
		skyVal = medianWirth(datVal->data, numGoodRows);
		var = 0.0;
		for (i = 0; i < numGoodRows; i++) {
		  xp = ipow(datVal->data[i] - skyVal, 2);
		  var += xp;
		}
		var /= (float) numGoodRows;
		rmsVal = (double) sqrt(var);
		for (i = 0; i < numGoodRows; i++) {
		  resid = datVal->data[i] - skyVal;
		  if (resid > nSigAbove*rmsVal || resid < -nSigBelow * rmsVal)
		  {
		    for (k = i; k < numGoodRows-1; k++) 
		      datVal->data[k] = datVal->data[k+1];
		    numGoodRows--;
		    nReject++;
		  }
		}
		nIter++;
		break;
	      }
	    case SKY_FIT : 
	      {
		if (skyCoefs != NULL) skyCoefs=NULL;
		rmsVal=0.0;
		skyCoefs = fit1DPoly(polyDeg, datPoint, numGoodRows, &rmsVal);
		rms = sqrt(rmsVal);	      
		for (i = 0; i < numGoodRows; i++) {
		  skyVal = skyCoefs[0];
		  for (k = 1; k <= polyDeg; k++) {
		    xp = ipow((float)i, k) ;
		    skyVal += skyCoefs[k] * xp ; 
		  }
		  resid = datPoint[i].y - skyVal;
		  if (resid > nSigAbove * (float) rms || resid < -nSigBelow * 
		      (float) rms)
		  {
		    for (k = i; k < numGoodRows-1; k++) {
		      datPoint[k].x = datPoint[k+1].x;
		      datPoint[k].y = datPoint[k+1].y;
		    }
		    numGoodRows--;
		    nReject++;
		  }
		}
		nIter++;
		break;
	      }
	    }
	  } while (nReject > 0 && nIter < maxIter && numGoodRows >= nSkyMin);
	
	  switch (skyMode) {
	  default: 
	  case SKY_MEDIAN : 
	    {
	      for (i = 0; i < numRows; i++) {
		outData->data[xLow + i + yOut*imageXlen] = 
		  inputImage->data[xLow + i + yOut*imageXlen] - skyVal;
		skyData->data[xLow + i + yOut*imageXlen] = skyVal;
	      }
	      break;
	    }
	  case SKY_FIT : 
	    {
	      for (i = 0; i < numRows; i++) {
		skyVal = skyCoefs[0];
		for (k = 1; k <= polyDeg; k++) {
		  xp = ipow((float)i, k) ;
		  skyVal += skyCoefs[k] * xp ; 
		}
		outData->data[xLow + i + yOut*imageXlen] = 
		  inputImage->data[xLow + i + yOut*imageXlen] - skyVal;
		skyData->data[xLow + i + yOut*imageXlen] = skyVal;
	      }
	      cpl_free(skyCoefs);
	      break;
	    }
	  }
	  
	}
      }

      switch (skyMode) {
      default: 
      case SKY_MEDIAN : 
	{
	  deleteFloatArray(datVal);
	  break;
	}
      case SKY_FIT : 
	{
	  deleteDpoint(datPoint);
	  break;
	}
      }
    }

    extSlit = extSlit->next;
    winSlit = winSlit->next;
  }

  copyAllDescriptors(inputImage->descs, &(outData->descs));
  copyAllDescriptors(inputImage->descs, &(skyData->descs));
  outputImages = (VimosImage **) cpl_malloc(2*sizeof(VimosImage *));
  outputImages[0] = outData;
  outputImages[1] = skyData;
  
  return outputImages;  
}


VimosImage **
VmSpSkyExt(VimosImage **outSpEx2D, VimosWindowTable *windowTable,
			int skyMode, int polyDeg, int limMargin, int nSkyMin,
			int maxIter, float nSigBelow, float nSigAbove)
{
  int i, j, k;
  int x, y;
  int numRows;
  int numGoodRows;
  VimosUlong32 numPoints;
  VimosUlong32 index;
  int imageXlen, imageYlen;
  int nReject, nIter;  
  
  float skyVal;
  float resid;
  
  double xp, var;
  double rmsVal;
  double *skyCoefs = 0;
  
  VimosFloatArray *datVal=NULL;
  VimosImage *imageData;
  VimosImage *outData;
  VimosImage *skyInData;
  VimosImage *skyData;
  VimosDpoint *datPoint=NULL;
  VimosImage **outputImages;
  
  VimosWindowSlit *winSlit;
  char modName[]   ="VmSpSkyExt";

  cpl_msg_debug(modName, "Subtracting Sky 2D");

  imageData = outSpEx2D[0];
  /*  imageData = openOldFitsFile(imageName, 1, 1); */
  
  imageXlen = imageData->xlen;
  imageYlen = imageData->ylen;

  
  skyInData = outSpEx2D[1];
  
  outData = newImageAndAlloc(imageXlen, imageYlen);
  skyData   = newImageAndAlloc(imageXlen, imageYlen);
    
  numPoints = imageXlen*imageYlen;
  
  for (index = 0; index < numPoints; index++) {
      outData->data[index] = imageData->data[index];
      skyData->data[index] = 0.0;
  }
  
  winSlit = windowTable->slits;
  
  while  (winSlit)  {
    numRows = winSlit->specEnd-winSlit->specStart+1;

    /* we just copy  short slits */
    if (winSlit->specLong == VM_FALSE) {
      for (x = 0; x < imageXlen; x++) {
        for (y = winSlit->specStart; y <= winSlit->specEnd; y++) {
          outData->data[x + y*imageXlen] = imageData->data[x + y*imageXlen];
          skyData->data[x + y*imageXlen] = skyInData->data[x + y*imageXlen];
        }
      }
      winSlit = winSlit->next;
      continue;
    }
    switch (skyMode) {
    default: 
    case DO_MEDIAN : 
      {
        datVal = newFloatArray(numRows-2*limMargin);
        break;
      }
    case DO_FIT : 
      {
        datPoint = newDpoint(numRows-2*limMargin);
        break;
      }
    }
    
    
    for (j = 0; j < imageXlen; j++) {
      y = j;
      /* loop over rows to get data for this wavelength */

      numGoodRows = 0;      
      for (i = limMargin; i < numRows-limMargin; i++) {
        if (windowObjectInRow(winSlit, i) == VM_FALSE) {
          switch (skyMode) {
          default: 
          case DO_MEDIAN : 
            {
              datVal->data[numGoodRows] = 
                imageData->data[y + (i+winSlit->specStart)*imageXlen];
              break;
            }
          case DO_FIT : 
            {
              datPoint[numGoodRows].x = i;
              datPoint[numGoodRows].y =                
                imageData->data[y + (i+winSlit->specStart)*imageXlen];
              break;
            }
          }
          numGoodRows++;
        }
      }
      
      nIter = 0;
      do {
        nReject = 0;
	switch (skyMode) {
	default: 
	case SKY_MEDIAN : 
	  {
	    skyVal = medianWirth(datVal->data, numGoodRows);
	    var = 0.0;
	    for (i = 0; i < numGoodRows; i++) {
	      xp = ipow(datVal->data[i] - skyVal, 2);
	      var += xp;
	    }
	    var /= (double) numGoodRows;
	    rmsVal = sqrt(var);
	    for (i = 0; i < numGoodRows; i++) {
	      resid = datVal->data[i] - skyVal;
	      if (resid > nSigAbove * rmsVal || resid < -nSigBelow * rmsVal)
		{
		  for (k = i; k < numGoodRows-1; k++)
		    datVal->data[k] = datVal->data[k+1];
		  numGoodRows--;
		  nReject++;
		}
	    }
            nIter++;
	    break;
	  }
	case SKY_FIT : 
	  {
	    if (skyCoefs != NULL) 
              cpl_free(skyCoefs);
            rmsVal = 0.0;
	    skyCoefs = fit1DPoly(polyDeg, datPoint, numGoodRows, &rmsVal);
            rmsVal = sqrt(rmsVal);
	    for (i = 0; i < numGoodRows; i++) {
	      skyVal = skyCoefs[0];
	      for (k = 1; k <= polyDeg; k++) {
		xp = ipow((float)i, k) ;
		skyVal += skyCoefs[k] * xp ; 
	      }
	      resid = datPoint[i].y - skyVal;
	      if (resid > nSigAbove * (float) rmsVal || resid < -nSigBelow * 
		                                        (float) rmsVal)
		{
		  for (k = i; k < numGoodRows-1; k++) {
		    datPoint[k].x = datPoint[k+1].x;
		    datPoint[k].y = datPoint[k+1].y;
		  }
		  numGoodRows--;
		  nReject++;
		}
	    }
            nIter++;
	    break;
	  }
	}
      } while (nReject > 0 && nIter < maxIter && numGoodRows >= nSkyMin);
	
      switch (skyMode) {
      default: 
      case DO_MEDIAN : 
        {
          for (i = 0; i < numRows; i++) {
            outData->data[y + (i+winSlit->specStart)*imageXlen] = 
              imageData->data[y + (i+winSlit->specStart)*imageXlen] - skyVal;
            skyData->data[y + (i+winSlit->specStart)*imageXlen] = skyVal;
          }
          break;
        }
      case DO_FIT : 
        {
          for (i = 0; i < numRows; i++) {
            skyVal = skyCoefs[0];
            for (k = 1; k <= polyDeg; k++) {
              xp = ipow((float)i, k) ;
              skyVal += skyCoefs[k] * xp ; 
            }

            if (y >= 0 && y < imageXlen) {
              if (i+winSlit->specStart >= 0 &&
                  i+winSlit->specStart < imageYlen) {
                outData->data[y + (i+winSlit->specStart)*imageXlen] = 
                imageData->data[y + (i+winSlit->specStart)*imageXlen] - skyVal;
                skyData->data[y + (i+winSlit->specStart)*imageXlen] = skyVal;
              }
              else
                cpl_msg_debug(modName, "y out of range %d %d %d %d\n", 
                            i+winSlit->specStart, i, winSlit->specStart, 
                            imageYlen);
            }
            else
              cpl_msg_debug(modName, "x out of range %d", y);
          }
          cpl_free(skyCoefs);
          skyCoefs = 0;
          break;
        }
      }
    }

    winSlit = winSlit->next;
    switch (skyMode) {
    default: 
    case DO_MEDIAN : 
      {
        deleteFloatArray(datVal);
        break;
      }
    case DO_FIT : 
      {
        deleteDpoint(datPoint);
        break;
      }
    }
  }

  copyAllDescriptors(imageData->descs, &(outData->descs));
  copyAllDescriptors(skyInData->descs, &(skyData->descs));
  
  outputImages = (VimosImage **) cpl_malloc(2*sizeof(VimosImage *));
  outputImages[0] = outData;
  outputImages[1] = skyData;

  return outputImages;
}

