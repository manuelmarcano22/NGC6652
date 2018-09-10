/* $Id: vmmosextraction.c,v 1.5 2013-08-07 18:30:07 cgarcia Exp $
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
 * $Date: 2013-08-07 18:30:07 $
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>
#include <pilerrno.h>

#include "vmimage.h"
#include "vmmatrix.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmifutable.h"
#include "vmgrismtable.h"
#include "vmextractiontable.h"
#include "vmwindowtable.h"
#include "vmobjecttable.h"
#include "vmdistmodels.h"
#include "vmadf.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmmosextraction.h"
#include "cpl.h"


#define OFFSET_TOLERANCE (2) /* FIXME Hardcoded! */
#define MIN_OFFSET_OK (1)    /* FIXME Hardcoded! */



/**
 * @brief vmmosextraction MOS Spectrum Extraction
 *
 * The module provides functions to detect object spectra and extract them.
 */

/**@{*/

/* Returns difference in seconds between two dates in ISO standard format */

static int diffdate(char *date1, char *date2) 
{
    int   jdate1, jdate2;
    int   year1, year2;
    int   month1, month2;
    int   day1, day2;
    int   hour1, hour2;
    int   min1, min2;
    int   sec1, sec2;
    int   diff;
    int   f;

    if (date1 == NULL || date2 == NULL) {
        cpl_msg_error("diffdate", "NULL input dates");
        return 0;
    }
    if (date1[4] != '-' ||
        date1[7] != '-' ||
        date1[10] != 'T' ||
        date1[13] != ':' ||
        date1[16] != ':') {
        cpl_msg_error("diffdate", "Date is not in ISO standard format (%s)",
                      date1);
        return 0;
    }
    if (date2[4] != '-' ||
        date2[7] != '-' ||
        date2[10] != 'T' ||
        date2[13] != ':' ||
        date2[16] != ':') {
        cpl_msg_error("diffdate", "Date is not in ISO standard format (%s)",
                      date2);
        return 0;
    }

    if (!strncmp(date1, date2, 19)) {
        diff = 0;
    }
    else {
        sscanf(date1, "%d-%d-%dT%d:%d:%d",
               &year1, &month1, &day1, &hour1, &min1, &sec1);
        sscanf(date2, "%d-%d-%dT%d:%d:%d",
               &year2, &month2, &day2, &hour2, &min2, &sec2);
        if (month1 < 3)
            f = -1;
        else
            f = 0;
        jdate1 = floor((1461 * (f + 4800 + year1))/4)
               + floor(((month1 - 2 - (f * 12)) * 367) / 12)
               - floor(3 * floor((year1 + 4900 + f) / 100) / 4)
               + day1 - 32075;
        jdate2 = floor((1461 * (f + 4800 + year2))/4)
               + floor(((month2 - 2 - (f * 12)) * 367) / 12)
               - floor(3 * floor((year2 + 4900 + f) / 100) / 4)
               + day2 - 32075;
        diff = (jdate1 - jdate2) * 24 * 3600;
        diff += (hour1 - hour2) * 3600 + (min1 - min2) * 60 + (sec1 - sec2);

    }

    return diff;

}

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosWindowObject *findObjectsInProfile(VimosFloatArray *profile, 
                                float detLevel, int numLev, float objFrac, 
                                float limFrac, int limMargin, 
                                int minObjectSize, int minCompositeSize);
   
   Description:
   Find objects in a 1D profile 
   
   First, the function thersholds the profile. For each contiguous region of
   pixels (longer than 4 pixels) it applied a watershed algorithm to see if
   more than one object are in that region. For each object found, a Window
   Object is created with the following parameters filled: 
   ->objStart
   ->objEnd
   ->objPos
   ->posDef = VM_FALSE
   ->pardef = VM_FALSE
   ->objX = 0.0
   ->objY = 0.0

   Input:
   VimosFloatArray *profile
   Float Array containg the spatial profile to scan.
   
   float detLevel
   Detection Level in rms units

   int numLev
   Number of levels to use in the watershed algorithm (should be something like 
   32, see VimosMath.h for doc).

   float objFrac
   Flux fraction to use in watershed (see VimosMath.h for doc)

   float limFrac
   not yet used

   int limMargin
   not yet used ********* excluded here and used in upper levels: VmSpDetObj, 
     VmSpSkyFra, VmSpSkyExt***********

   int minObjectSize
   minimum size for considering the detection an object

   int minCompositeSize
   minimum size for trying to deblend objects. If a candidate is found of size 
   smaller than minCompositeSize, it is taken as one object

   Return Value:
   pointer to newly created Window Object(s). NULL if no objects were found. 
   
   Updates:
   21 Apr 99: Created (TAO)

-------------------------------------------------------------------------------- 
*/

static VimosWindowObject *
findObjectsInProfile(VimosFloatArray *profile, float detLevel, int numLev,
                     float objFrac, float limFrac, int minObjectSize,
                     int minCompositeSize)
{
  /* find contiguous ranges above threshold
     apply watershed for each contiguous range
     extend limits
     create windowobjects and return them
  */
  int i, j, index;
  int numPoints;
  int objectNum;
  int numObjects=0;
  int totalNumObjects;
  int specStart;
  int pos;
  (void)limFrac;
    
  /* should become a parameter, but is not in DRS document */
  int smoothWidth = 1;
  float threshold;
  float sum,weights,norm;
  
  float meanInt, rmsInt, junk1, junk2;

  VimosFloatArray *fBuffer;
  VimosIntArray   *iBuffer;
  
  VimosWindowObject *tmpObject;
  VimosWindowObject *lastObject;
  char               modName[] = "findObjectsInProfile";

  /* this is not enough data */
  if (profile->len < 5) {
    cpl_msg_error(modName,"There is not enough data");
    return(NULL);
  }
 
  /* create internally used buffers */
  iBuffer = newIntArray(profile->len);
  if (iBuffer == NULL) {
    cpl_msg_error(modName,"function newIntArray failure");
    return(NULL);
  }
  fBuffer = newFloatArray(profile->len);
  if (fBuffer == NULL) {
    cpl_msg_error(modName,"The function newFloatArray has returned NULL");
    return(NULL);
  }
  /* initilse counter and pointer */
  totalNumObjects = 0;
  lastObject = NULL;
  
  /* determine mean and variance of the data using the biweight estimator */
  xbiwt(profile->data, profile->len, &meanInt, &rmsInt, &junk1, &junk2);
  /* then apply a 2-sigma cut only to positive outliers */
  numPoints=0; 
  for (i = 0; i < profile->len; i++)
    {
      if (profile->data[i] <= (1. * rmsInt + meanInt)) 
        {
          fBuffer->data[numPoints] = profile->data[i];
          numPoints++;
        }
    }
  
  /* and redetermine the mean and variance, used to define the threshold */
  xbiwt(fBuffer->data, numPoints, &meanInt, &rmsInt, &junk1, &junk2);

  if (rmsInt <= 0.0000001) {
    deleteIntArray(iBuffer);
    deleteFloatArray(fBuffer);
    return(NULL);
  }

  threshold = meanInt + detLevel * rmsInt;
  /*printf("meanInt =%f, rmsInt =%f, thres = %f\n", 
    meanInt, rmsInt, threshold); */
  /* start at the beginning  (where else...) */
  i = 0;
  /* loop though data */
  do {
    /* skip points that are below thrshold */
    while ( (i < profile->len) && (profile->data[i] < threshold) ) {
      
      i++;
    }
    if(i >= (profile->len-minObjectSize)) break;

    /* if we are not at the end of the profile, scan for points above the
       threshold */
    else {
      /* initialise counter */
      numPoints = 0;
      /* store index of beginning for later use */
      specStart = i;
      /* scan as long as data is above threshold */
      while ( (i < profile->len) && (profile->data[i] >= threshold) ) {
        /* copy data in buffer */
        /*printf("i =%d, profile =%f\n", i, profile->data[i]); */
        fBuffer->data[numPoints] = profile->data[i];
        /* increment counter */
        numPoints++;
        /* increment index */
        i++;
      }
      
      /**************
        Here add check on objects: exclude profiles that do not have a 
        "object shape" (maximum is at the borders)
      ***********/

      if (numPoints > minCompositeSize) {
        /* now we have a contiguous range of pixels above the threshold, to
           which we apply the watershed to see if this range contains more than
           one object */
        numObjects = waterShed(fBuffer->data, numPoints, numLev, objFrac,
                               smoothWidth, iBuffer->data);
        /* here I set numObjects = 1, when the watershed returns 0
         * This means that objects candidates are only one: cannot be
         * separated, it is a candidate extended object.
	 * The same is true when numObjects = 1
         */
        if (numObjects == 0 || numObjects == 1) {
          numObjects = 1;
          for (index = 0; index < numPoints; index++) {
            iBuffer->data[index] = 1;
          }
        }
         
      } else if (numPoints > minObjectSize) {
        /* range too short for deblending, only one object */
        numObjects = 1;
        for (index = 0; index < numPoints; index++) {
          iBuffer->data[index] = 1;
        }  
      } else {
        /* if only two pixels are above the threshold, we don't think 
           there is an object */
        continue;
      }
      /* extract the info of each object and create WindowObjects for them */
      for (objectNum = 1; objectNum <= numObjects; objectNum++) {
        tmpObject = newWindowObject();
        if (tmpObject == NULL) {
          cpl_msg_error(modName,"The function newWindowObject has returned NULL");
          return(NULL);
        }
        index = 0;
        /* One could write the loop as
           while (iBuffer->data[index++] != objectNum) ;
           but I don't like this kind of thing..
           I do not check on the range of index here, because I know at least
           one pixel in iBuffer has the right value
        */ 
        
        while (iBuffer->data[index] != objectNum) {
          index++;
        }
        tmpObject->objStart = specStart + index;


        /* Next, we scan for the end of this object. Here I have to check on
           the range of index, namely to trap the case that
           iBuffer->data[numPoints-1] == objectNum */
       
        while ( (index < numPoints) && (iBuffer->data[index] == objectNum) ) {
          index++;
        }
        tmpObject->objEnd = specStart + index -1;
        
        /* printf("specStart=%d, objstart=%d, objend=%d\n", specStart,
                tmpObject->objStart, tmpObject->objEnd); */

        /* record sequential number for this slit */
        tmpObject->objNo = totalNumObjects + objectNum;


	/* compute the object position as the baricenter of the data 
	   and store the object profile */

	tmpObject->objProfile=newFloatArray(tmpObject->objEnd-tmpObject->objStart+1);
	for (j=tmpObject->objStart, sum=0., weights=0., norm=0.; 
	     j<=tmpObject->objEnd; j++) {
	  tmpObject->objProfile->data[j-tmpObject->objStart] = profile->data[j] - threshold;
	  weights += profile->data[j];
	  sum     += j*profile->data[j];
	  norm    += tmpObject->objProfile->data[j-tmpObject->objStart];
	}
	tmpObject->objPos = sum/weights;


	/* P_FIXME normalize the object profile */
	for (j=0; j<tmpObject->objProfile->len; j++) {
	  tmpObject->objProfile->data[j] /= norm;
	}
	
	  
	/* object width as the intensity of the 2 nearest pixel */
	
	pos=(int)(tmpObject->objPos);
	tmpObject->objWidth = profile->data[pos] + profile->data[pos+1];  



        /* mask positions (tmpObject->objX and ->objY) should be filled
           outside this routine, since I don't have the information here. Set
           them to zero, just in case... */
        tmpObject->objX = tmpObject->objY = 0.0;
        /* no other parameter are defined */
        tmpObject->posDef = tmpObject->parDef = VM_FALSE;
        
        /* link objects */
        if (lastObject == NULL) {
          /* this is the first object of the linked list */
          lastObject = tmpObject;
        }
        else {
          /* append to existing linked list */
          lastObject->next = tmpObject;
          tmpObject->prev = lastObject;
          /* and copy address of last in list */
          lastObject = tmpObject;
        }
      }
      
    }
    /* keep on counting */
    totalNumObjects += numObjects;
      
  }

  /* loop until we scanned the whole slit  */
  while  (i < profile->len);

  /* looks like we found something in this slit */
  if (totalNumObjects > 0) {
    /* rewind linked list */
    while (lastObject->prev != NULL) {
      lastObject = lastObject->prev;
    }
  }

  /* delete internal buffers */
  deleteIntArray(iBuffer);
  deleteFloatArray(fBuffer);
  
  /* return what we found */
  return(lastObject);
  
    
}


/* private function for specExtract1D. Does the actual optimal  extraction */

static VimosBool
horneExtract(VimosImage *inData, VimosImage *skyData, float ron, float conad, 
             VimosFloatArray *outSpec, VimosFloatArray *outSkySpec, 
             int horne, int ncomb)
{
  int i, j;
  int specLen;
  int numRows;
  int index;
  int iter;
  int maxIter   = 2;         /* Not less than 2 !!! */
  int smoothBox = 21;        /* Not less than 5 !!! */
  
  double sumWeight, sum, sumSky, variance, weight;
  float *profile;
  float *buffer;

  double value;

  
  specLen = inData->xlen;
  numRows = inData->ylen;
  
  /*
   * Initial spectrum estimate 
   */

  for (i = 0; i < specLen; i++)
    for (j = 0, outSpec->data[i] = 0.0; j < numRows; j++)
      if (inData->data[i + j * specLen] > 0.0)
        outSpec->data[i] += inData->data[i + j * specLen];


  if (horne) {

    profile = (float *)cpl_calloc(specLen * numRows, sizeof(float));
    buffer  = (float *)cpl_calloc(specLen, sizeof(float));

    for (iter = 0; iter < maxIter; iter++) {

      /*  
       * Spatial profile
       */

      for (i = 0; i < specLen; i++) {
        for (j = 0; j < numRows; j++) {
          index = i + j * specLen;
          if (inData->data[index] > 0.0 && outSpec->data[i] > MIN_DIVISOR) 
            profile[index] = inData->data[index] / outSpec->data[i];
          else
            profile[index] = 0.0;
        }
      }

      for (j = 0; j < numRows; j++) {
  
        /*
         * Smooth each row in the dispersion direction, and enforce positivity
         */
  
        for (i = 0; i < specLen - smoothBox; i++) {
          value = medianPixelvalue(profile + i + j * specLen, smoothBox);
          if (value < 0)
            value = 0.0;
          buffer[i + smoothBox / 2] = value;
        }
  
        /*
         * Replace the end portions (i.e., not median filtered) with a mean
         */
  
        value = computeAverageFloat(profile + j * specLen, smoothBox / 2);
        if (value < 0)
            value = 0.0;
        for (i = 0; i < smoothBox / 2; i++)
          buffer[i] = value;
  
        value = 
        computeAverageFloat(profile + specLen - smoothBox / 2 + j * specLen,
                            smoothBox / 2);
  
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
        if (value > MIN_DIVISOR)
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
        for (j = 0; j < numRows; j++) {
          index = i + j * specLen;
          if (inData->data[index] > 0.0) {
            variance = (ron / conad) * (ron / conad) 
                     + fabs(outSpec->data[i] * profile[index] 
                            + skyData->data[index])
                     / conad;
            variance /= ncomb;  /* If input dataset is sum of ncomb images */
            value = inData->data[index] - outSpec->data[i] * profile[index];
            if (fabs(value) / sqrt(variance) < 5.0) {
              weight = profile[index] / variance;
              sum += weight * inData->data[index];
              sumSky += weight * skyData->data[index];
              sumWeight += weight * profile[index];
            }
          }
        }
        if (sumWeight > MIN_DIVISOR) {
          outSpec->data[i] = sum / sumWeight;
          outSkySpec->data[i] = sumSky / sumWeight;
        }
        else {
          outSpec->data[i] = 0.0;
          outSkySpec->data[i] = 0.0;
        }
      }
    }
    cpl_free(profile);
    cpl_free(buffer);
  }
  else {

  /*
   * Add sky estimation for the simple aperture extraction.
   */

    for (i = 0; i < specLen; i++)
      for (j = 0, outSkySpec->data[i] = 0.0; j < numRows; j++)
        if (skyData->data[i + j * specLen] > 0.0)
          outSkySpec->data[i] += skyData->data[i + j * specLen];

  }
  
  return VM_TRUE;

}


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool specExtract2D(VimosImage *imageIn, VimosExtractionTable *exTable, 
                           VimosImage *imageOut, VimosSpecSampleOption
                           sampleOption, char *kernelType);
   
   Description:
   2D extracts the spectra defined in the Extraction Table from Frame format
   to Stacked Spectra Format.
   Code is based on the function warp_image from eclipse 3.1.1.
   
   Input: 
   VimosImage   *imageIn
   Pointer to the image containg the spectra in Frame Format

   VimosExtractionTable *exTable
   Pointer of Extraction Table containing the information about the slits

   VimosSpecSampleOption sampleOption
   Sample linear in wavelength (VM_SP_LIN_LAMBDA) or in log wavelength
   (VM_SPEC_LOG_LAMBDA) 

   char *kernelType
   Not used

   Input/Output: 
   VimosImage   *imageOut
   Pointer to the output image. The output image should exist (nothing is
   allocated in this function) and should be large enough to contain the
   spectra. 

   
   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   14 Jun 00: Return VimosBool instead of void (Maura)
   10 Feb 99: Created (TAO)

------------------------------------------------------------------------------ 
*/

static VimosBool
specExtract2D(VimosImage *imageIn, VimosExtractionTable *exTable, 
              VimosImage *imageOut, VimosSpecSampleOption sampleOption, 
              char *kernelType/* not used*/)
{
    const char   modName[] = "specExtract2D";

    int          i, j, k;
    double       cur;
    double       neighbors[16];
    double       rsc[8];
    double       sumrs;
    double       x, y;
    double       stepFactor = 0.;
    
    VimosUlong32 px, py;
    VimosUlong32 pos;
    int          tabx, taby;
    double      *kernel;
    VimosLong32 *leaps;

    int          inpXlen;
    int          outXlen;
    int          specLen;
    int          numRows;
    int          numSpec;
    float        wLenStart;
    float        wLenEnd;
    float        wLenInc;
    float        lambda;
    double       dValue;
    char         comment[80];

    VimosExtractionSlit  *slit;
    
    (void)kernelType; /* To avoid compiler warning */

    pilErrno = 0;

    /* error handling: test entries */
    if (imageIn == NULL) {
      cpl_msg_error(modName, "NULL input image");
      return(VM_FALSE);
    }
    
    /* setup arrays for the interpolation*/
    inpXlen = imageIn->xlen;

    leaps = (VimosLong32 *) cpl_malloc(16*sizeof(VimosLong32));

    /* check if space was allocated */
    if (leaps == NULL) {
      cpl_msg_error(modName, "Allocation Error");
      return(VM_FALSE);
    } 
    
    /* setup the interpolation kernel and the vector with the leaps */
    if (!setupInterpolation(&kernel, &leaps, inpXlen)) {
      cpl_msg_error(modName, "Function setupInterpolation failure");
      return(VM_FALSE);
    }

    /* get required length of spectrum and wavelength step */

    if (readDoubleDescriptor(exTable->descs, pilTrnGetKeyword("WlenStart"), 
                             &dValue, comment) == VM_FALSE) {
      cpl_msg_error(modName, "Descriptor %s not found", 
                  pilTrnGetKeyword("WlenStart"));
      return VM_FALSE;
    }
    wLenStart = (float) dValue;

    if (readDoubleDescriptor(exTable->descs, pilTrnGetKeyword("WlenEnd"), 
                             &dValue, comment) == VM_FALSE) {
      cpl_msg_error(modName, "Descriptor %s not found", 
                  pilTrnGetKeyword("WlenEnd"));
      return VM_FALSE;
    }
    wLenEnd = (float) dValue;

    if (readDoubleDescriptor(exTable->descs, pilTrnGetKeyword("WlenInc"), 
                             &dValue, comment) == VM_FALSE) {
      cpl_msg_error(modName, "Descriptor %s not found", 
                  pilTrnGetKeyword("WlenInc"));
      return VM_FALSE;
    }
    wLenInc = (float) dValue;
   
    /* compute required length */
    specLen = (wLenEnd-wLenStart)/wLenInc + 1;

    switch (sampleOption) {
      default: 
      case VM_SP_LIN_LAMBDA : {
        break;
      }
      case VM_SP_LOG_LAMBDA : {
        stepFactor = log10(wLenEnd/wLenStart)/(1.0*(specLen-1));
        stepFactor = pow(10.0, stepFactor);
        break;
      }
    }
        
    /* length X-axis of output image */
    outXlen = imageOut->xlen;
    /* if output image too small : abort */
    if (specLen < outXlen) {
      cpl_msg_error(modName, "Output image too small in X");
      return(VM_FALSE);
    }
    
    /* check Y-size of output image */
    slit = exTable->slits;
    numRows = 0;
    /* loop through slits and count number of output rows needed */
    while (slit) {
      numRows += slit->numRows;
      slit = slit->next;
   }
    
    /* if image too small: abort */
    if (numRows < imageOut->ylen) {
      cpl_msg_error(modName, "Output image too small in Y");
      return(VM_FALSE);
    }
    numSpec = 0;
    
    /* reset to first slit */
    slit = exTable->slits;

    /* while there are slits defined...  */
    while (slit) {
      /* get number of spectra to do for this slit */
      numRows = slit->numRows;
      /* loop of row within slit to update locations of spectra in image */
      for (i = 0; i < numRows; i++) {
        slit->y->data[i] = numSpec+i;
      }
      
      lambda = wLenStart;
      
      /* loop over wavelength */
      for (j = 0; j < specLen; j++) {
        /* loop of row within slit */
        for (i = 0; i < numRows; i++) {

          /* check that this row has a wavelength solution */
          if (slit->invDisQuality->data[i] == 0)
            continue;

          /* Compute the  source for this pixel 
             first by inverting wavelength to Y-pixel  */
          y = computeDistModel1D(slit->invDis[i], lambda);
          if (pilErrno) {
            cpl_msg_error(modName, "Function computeDistModel1D failure (y)");
            return(VM_FALSE);
          }
         
          /* adding position of slit on CCD */
          y += slit->ccdY->data[i];
          
          /* use curvature polynomial to find X-pixel */
          x = computeDistModel1D(slit->crvPol[i], y);
          if (pilErrno) {
            cpl_msg_error(modName, "Function computeDistModel1D failure (x)");
            return(VM_FALSE);
          }
          x += slit->ccdX->data[i];
         
          /* Which is the closest integer  neighbour?    */
          px = (VimosUlong32) x;
          py = (VimosUlong32) y;
          
          /* if too close to border: set to zero */
          if ((px < 1) ||
              (px > (VimosUlong32) (imageIn->xlen-3)) ||
              (py < 1) ||
              (py > (VimosUlong32) (imageIn->ylen-3))) {
            imageOut->data[j+(numSpec+i)*outXlen] = 0.0;
          } else {
            /* Now feed the positions for the closest 16 neighbours  */
            pos = px + py * inpXlen;
            for (k = 0; k < 16; k++) {
              neighbors[k] =  
                (double)(imageIn->data[(VimosUlong32)((VimosLong32)pos
                                                      +leaps[k])]);
            }
            
            /* Which tabulated value index shall we use for the kernel? */
            tabx = (x - (double)px) * (double)(TABSPERPIX); 
            taby = (y - (double)py) * (double)(TABSPERPIX); 
            
            /* get resampling coefficients from kernel array  */
            /* rsc[0..3] in x, rsc[4..7] in y   */
            
            rsc[0] = kernel[TABSPERPIX + tabx];
            rsc[1] = kernel[tabx];
            rsc[2] = kernel[TABSPERPIX - tabx];
            rsc[3] = kernel[2 * TABSPERPIX - tabx];
            rsc[4] = kernel[TABSPERPIX + taby];
            rsc[5] = kernel[taby];
            rsc[6] = kernel[TABSPERPIX - taby];
            rsc[7] = kernel[2 * TABSPERPIX - taby];
            
            /* compute normalization */
            sumrs = (rsc[0]+rsc[1]+rsc[2]+rsc[3]) *
              (rsc[4]+rsc[5]+rsc[6]+rsc[7]);
            
            
            /* Compute interpolated pixel */
            cur = rsc[4] * (rsc[0]*neighbors[0]   +
                            rsc[1]*neighbors[1]   +
                            rsc[2]*neighbors[2]   +
                            rsc[3]*neighbors[3])  +
                  rsc[5] * (rsc[0]*neighbors[4]   +
                            rsc[1]*neighbors[5]   +
                            rsc[2]*neighbors[6]   +
                            rsc[3]*neighbors[7])  +
                  rsc[6] * (rsc[0]*neighbors[8]   +
                            rsc[1]*neighbors[9]   +
                            rsc[2]*neighbors[10]  +
                            rsc[3]*neighbors[11]) +
                  rsc[7] * (rsc[0]*neighbors[12]  +
                            rsc[1]*neighbors[13]  +
                            rsc[2]*neighbors[14]  +
                            rsc[3]*neighbors[15]); 
            
            /* write output value */
            imageOut->data[j+(numSpec+i)*outXlen] = (float)(cur/sumrs);
          }       
        }
        switch (sampleOption) {
          default: 
          case VM_SP_LIN_LAMBDA : {
            lambda = wLenStart + (j + 1) * wLenInc;
            break;
          }
          case VM_SP_LOG_LAMBDA : {
            lambda *= stepFactor;
            break;
          }
        }
      }
      /* increment counter of number of rows done */
      numSpec += numRows;

      /* get pointer to next slit */
      slit = slit->next;
    }
   
    /* free the space of the kernel */
    cpl_free(kernel);
    return(VM_TRUE);
}


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool specExtract1D(VimosImage *imageIn, VimosImage *skyImage, 
                           VimosWindowTable *winTable, VimosTable *ccdTable, 
                           VimosImage *outImage, VimosObjectTable *objTable)

   
   Description:

   1D extracts the spectra defined in the Window Table from Stacked Spectra
   Format to Packed  Spectra Format using optimal extraction as described by
   Keith Horne in PASP, 98, 609 (1986).
   The code is based on the function horne from the ESO MOS context, but was
   reworked to adapt to the Virmos data structure (and to make the code
   readable....) 
   
   Input: 
   VimosImage   *imageIn
   Pointer to the image containg the spectra in Stacked Spetcra Format

   VimosImage   *skyImage
   Pointer to the image containg the sky spectra in Stacked Spetcra Format

   VimosWindowTable *winTable
   Pointer to Window  Table containing the information about the slits and
   objects 

   VimosTable *ccdTable
   Pointer to CCD Table (for Read-Out Noise)

   Input/Output: 
   VimosImage   *outImage
   Pointer to the output image. The output image should exist (nothing is
   allocated in this function) and should be large enough to contain the
   spectra. 

   Output:
   VimosObjectTable *objTable
   Address of pointyer to Object Table of output image. Object Table is 
   NOT created by specExtract1D
   
   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   15 Jun 00: Return VimosBool instead of void (Maura)  
   11 Feb 99: Created (TAO)

------------------------------------------------------------------------------ 
*/

static VimosBool
specExtract1D(VimosImage *imageIn, VimosImage *skyImage,
              VimosWindowTable *winTable,
              VimosImage *outImage, VimosImage *outSkyImage,
              VimosObjectTable *objTable, int horne, int ncomb)
{
  int i,j;
  
  int specLen;
  int numRows, rowsBelow, rowsAbove, rowStart, rowEnd;
  int objNumber;
  double dValue;
  char comment[80];
  char modName[] = "specExtract1D";  

  VimosUlong32 index1, index2;
    
  float readOutNoise;
  float adu2Electron;

  VimosImage *inData;
  VimosImage *skyData;
  VimosFloatArray *outSpec;
  VimosFloatArray *outSkySpec;
  
  VimosWindowSlit *slit;
  VimosWindowObject *winObj;
  VimosObjectObject *lastObj=NULL;
  VimosObjectObject *tmpObj;

  /* get RON & CONAD from input data image */
  if (readDoubleDescriptor(imageIn->descs, pilTrnGetKeyword("ReadNoise", 1), 
                           &dValue, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not present in CCD Table", 
                pilTrnGetKeyword("ReadNoise", 1));
    return VM_FALSE;
  }
  readOutNoise = (float) dValue;
  if (readDoubleDescriptor(imageIn->descs,
                           pilTrnGetKeyword("Adu2Electron", 1), 
                           &dValue, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not present in CCD Table", 
                pilTrnGetKeyword("Adu2Electron", 1));
    return VM_FALSE;
  }
  adu2Electron = (float) dValue;

  /* explicit initialization */
  outSkySpec = outSpec = NULL;
  skyData = inData = NULL;

  /* length of spectra in wavelength */  
  specLen = imageIn->xlen;

  /* first slit in Window Table */  
  slit = winTable->slits;
  
  /* counter for all objects in image */  
  objNumber = 0;

  /* loop over all slits in Window Table */
  while (slit) {
    /* get first object for this slit */
    winObj = slit->objs;
    /* loop over all objects in this slits */
    while (winObj){
      /* cleanup previous allocations, if they exist */
      deleteImageAndAlloc(inData);
      deleteImageAndAlloc(skyData);
      deleteFloatArray(outSpec);
      deleteFloatArray(outSkySpec);
      

      /* 
       * Number of spatial rows for this object.
       * Extra rows are added to allow Horne to work
       */

      numRows = winObj->objEnd - winObj->objStart + 1;  

      rowsBelow = rowsAbove = numRows / 2;  /* Attempt to double range */

      /* ... then ensure we don't go beyond slit margins */

      if (winObj->objStart < rowsBelow)
        rowsBelow = winObj->objStart;

      rowStart = slit->specStart + winObj->objStart - rowsBelow;

      if (slit->specEnd - slit->specStart - winObj->objEnd < rowsAbove)
        rowsAbove = slit->specEnd - slit->specStart - winObj->objEnd;

      rowEnd = slit->specStart + winObj->objEnd + rowsAbove;

      numRows += rowsBelow + rowsAbove;

      /* allocate work arrays */
      inData   = newImageAndAlloc(specLen, numRows);
      if (inData == NULL) {
	cpl_msg_error(modName,"function newImageAndAlloc failure");
	return(VM_FALSE);
      }
      skyData  = newImageAndAlloc(specLen, numRows);
      if (skyData == NULL) {
	cpl_msg_error(modName,"function newImageAndAlloc failure");
	return(VM_FALSE);
      }

      outSpec  = newFloatArray(specLen);
      if (outSpec == NULL) {
	cpl_msg_error(modName,"function newFloatArray failure");
	return(VM_FALSE);
      }

      outSkySpec  = newFloatArray(specLen);
      if (outSkySpec == NULL) {
        cpl_msg_error(modName,"function newFloatArray failure");
        return(VM_FALSE);
      }


      /* ALEX: when IFU, do nothing for dead fibers (IFUfibTrans = -1) */
      if (slit->IFUfibTrans != -1.)
	{
	 /* extract the 2D array of data for this object and store in
	    work arrays */
	 for (i = 0; i < specLen; i++ ) {
	   for (j = rowStart; j < rowEnd; j++) {
	     index1 = i + (j - (slit->specStart + winObj->objStart 
                                                - rowsBelow)) * specLen; 
	     index2 = j;
	     index2 *=specLen;
	     index2 +=i;
             if (index2 < imageIn->xlen * imageIn->ylen) {
	       inData->data[index1] = imageIn->data[index2];
	       skyData->data[index1] = skyImage->data[index2];
             }
	   }
	 }

	 /* do the actual extraction on the work arrays */
	 if (!horneExtract(inData, skyData, readOutNoise, 
			   adu2Electron, outSpec, outSkySpec, horne, ncomb)) {
	   cpl_msg_error(modName,"function horneExtract failure");
	   return(VM_FALSE);
	 }

	 /* copy data to output frame */
	 /* should we also give a sky spectrum??? */
	 for (i = 0; i < specLen; i++ ) {
	   outImage->data[i + objNumber*specLen] = outSpec->data[i];
	   outSkyImage->data[i + objNumber*specLen] = outSkySpec->data[i];
	 }
	}
      else if (slit->IFUfibTrans == -1.)
	{

	 for (i = 0; i < specLen; i++ )
	   outImage->data[i + objNumber*specLen] = 0.;

	} /* ALEX: end if dead fibers */


      /* update Object Table */
      tmpObj = newObjectObject();
      if (tmpObj == NULL) {
	cpl_msg_error(modName,"function newObjectObject failure");
	return(VM_FALSE);
      }
      tmpObj->slitNo = slit->slitNo;
      tmpObj->rowNum = objNumber;
      tmpObj->objX = winObj->objX;
      tmpObj->objY = winObj->objY;
      tmpObj->objRA = winObj->objRA;
      tmpObj->objDec = winObj->objDec;

      /* 
       * ALEX: add IFUslitNo and IFUfibNo parameters of tmpObj:
       * copy them from winSlit. When MOS, they will be zero by default; 
       * when IFU they will be properly set 
       */

      tmpObj->IFUslitNo = slit->IFUslitNo;
      tmpObj->IFUfibNo = slit->IFUfibNo;

      /* NOTE: should we include IFUfibTrans parameter 
	 in object Table as well ? (AZ) */

      /* add object to Object Table */
      if (objTable->objs == NULL) {
	/* first object to add */
	objTable->objs = tmpObj;
      }
      else {
	/* not first object */
	/* just update links */
	lastObj->next = tmpObj;
	tmpObj->prev = lastObj;
      }
      /* store pointer to last object to be used in next iter */
      lastObj = tmpObj;

      /* get pointer to next object in this slit */
      winObj = winObj->next;
      /* increment counter of objects done */
      objNumber++;

    }
    /* get pointer to next slit in Window Table */

    slit = slit->next;
  }
  
  /* cleanup local allocations */
  deleteImageAndAlloc(inData);
  deleteImageAndAlloc(skyData);
  deleteFloatArray(outSpec);
  deleteFloatArray(outSkySpec);
      

  return(VM_TRUE);
}


/**
 * @memo
 *   Find objects in a 2D extracted MOS frame.
 *
 * @return window Table.
 *
 * @param imageData2D            
 * @param extractionTable
 * @param numLev
 * @param limMargin
 * @param detLevel
 * @param objFrac
 * @param limFrac
 * @param minObjectSize
 * @param minCompositeSize
 * @param slitTolerance
 * @param specFrac
 *
 * @doc  
 *   This functions is used in the data reduction recipes VmMosObsSingle
 *   and VmMosObsSpectra to detect object in spectra. It is used on the
 *   MOS image in frame format, for searching the high signal
 *   object, then it is used on the sky subtracted 2D extracted format  
 *   iteration of VmSpDetObj that does a second search on data after the sky
 *   has been subtracted to find low signal objects.
 * 
 * @author 
 *   Tim Osterlo (modified P.Sartoretti to work also with 2D extracted image)
 */

VimosWindowTable *
VmSpDetObj(VimosImage *imageData, VimosExtractionTable *extractionTable,
           int numLev, int limMargin, float detLevel, float objFrac, 
           float limFrac, int minObjectSize, int minCompositeSize,
           float slitTolerance, float specFrac)
{
/*   int imageXlen, imageYlen; */
  int outXlen, outYlen;
  int pixStart, pixEnd;
  int pixStartSpec, pixEndSpec;
  int specLo, specHi;
  int x, y;
  int lastPix;
  int numRows;
  int yStart;
  int objPos;
  int flag = 0;

  float wLenStart;
  float wLenEnd;
  float wLenInc;
  float wLenCen;
  float avgPos;
  float avgRng;
  float tmpFloat;
  
  VimosImage *outData;
  VimosExtractionSlit  *extractionSlit;
  VimosWindowTable    *windowTable;
  VimosWindowSlit     *windowSlit;
  VimosWindowSlit     *lastWindowSlit;
  VimosWindowObject   *windowObject;
  
  VimosFloatArray *profile;
  
  char modName[]      = "VmSpDetObj";


  cpl_msg_debug(modName, "Finding Objects");

/*  imageXlen = imageData->xlen;
  imageYlen = imageData->ylen; */

  
  /* get required length of spectrum and wavelength step */
  readFloatDescriptor(extractionTable->descs, "ESO PRO WLEN START", 
                      &wLenStart, NULL);
  readFloatDescriptor(extractionTable->descs, "ESO PRO WLEN END", 
                      &wLenEnd, NULL);
  readFloatDescriptor(extractionTable->descs, "ESO PRO WLEN INC", 
                      &wLenInc, NULL);
  readFloatDescriptor(extractionTable->descs, "ESO PRO WLEN CEN", 
                      &wLenCen, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN LO", 
                      &specLo, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN HI", 
                      &specHi, NULL);
  readFloatDescriptor(extractionTable->descs, "ESO PRO AVG POS", 
                      &avgPos, NULL);
  readFloatDescriptor(extractionTable->descs, "ESO PRO AVG RNG", 
                      &avgRng, NULL);
  
  /* checks (likely redundant) on params */
  if (wLenStart > wLenEnd) {
    tmpFloat  = wLenStart;
    wLenStart = wLenEnd;
    wLenEnd   = tmpFloat;
  }  
  wLenInc = fabs(wLenInc);
  
  /* compute required length */
  outXlen = (wLenEnd-wLenStart)/wLenInc + 1;
  
  /* check Y-size of output image */
  extractionSlit = extractionTable->slits;
  outYlen = 0;
  /* loop through slits and count number of output rows needed */
  while (extractionSlit) {
    outYlen += extractionSlit->numRows;
    extractionSlit =extractionSlit->next;
  }

  /*
   * Here check on image dimensions to see if input image is 
   * already 2D extracted 
   */
  
  if ((imageData->xlen != outXlen) && (imageData->ylen != outYlen)) 
    flag = 1;

  if (flag) {
    outData = newImageAndAlloc(outXlen, outYlen);

    specExtract2D(imageData, extractionTable, outData, VM_SP_LIN_LAMBDA,
		  "default");
    /*  createFitsImage("ima2Dextr.fits", outData, "cat"); */
  }
  else {
    outData = imageData;
  }
  
  
  /* first pixel to use for averaging */
  pixStart = (avgPos - avgRng/2.0 - wLenStart)/wLenInc;
  /* checks on index (probably redundant) */
  

 /*
  *  If the image is already 2D extracted check if object region 
  *  is within the useful spectrum (this is because only the
  *  useful spectrum is sky-subtracted)
  *  Instead, if the image is in frame format (i.e. not sky subtracted)
  *  we can still choose a region for looking for object that will
  *  then be excluded from the useful spectrum, therefore the check is done
  *  on image dimensions 
  */

 
  pixStartSpec = (wLenCen - wLenStart)/wLenInc - specLo;
  if (pixStartSpec < 0) pixStartSpec =0;
  pixEndSpec = (wLenCen -wLenStart)/wLenInc + specHi;
  if (pixEndSpec > outXlen) pixEndSpec = outXlen;
  
  if (flag) {
    pixStartSpec = 0;
    pixEndSpec = outXlen;
  }
    
  if (pixStart < pixStartSpec) {
    pixStart = pixStartSpec;
  }
  if (pixStart >= pixEndSpec) {
    pixStart = pixEndSpec;
  }
  
  /* last  pixel to use for averaging */
  pixEnd = (avgPos + avgRng/2.0 - wLenStart)/wLenInc;
  /* checks on index (probably redundant) */
  if (pixEnd < pixStartSpec) {
    pixEnd = pixStartSpec;
  }
  if (pixEnd >= pixEndSpec) {
    pixEnd = pixEndSpec-1;
  }
  
/*    printf("pixstart = %d, pixend =%d, pixstartspec =%d, pixendspec=%d, Xlen=%d\n", pixStart, pixEnd, pixStartSpec, pixEndSpec, outXlen); */
  /* create output table */
  windowTable = newWindowTable();
  /* copt descriptors of Extraction Table to Window Table */
  copyExtTab2WinTab(extractionTable, windowTable);
  /* initialise pointer to tail of Window Slits to NULL */
  lastWindowSlit = NULL;

  /* just to be sure... */
  profile = NULL;
  /* get first slit in Extraction Table */
  extractionSlit = extractionTable->slits;
  
  while (extractionSlit) {
    yStart = extractionSlit->y->data[0]+limMargin;
/*      printf("yStart =%d\n", yStart); */
    numRows = extractionSlit->numRows-2*limMargin;
    profile = newFloatArray(numRows);
    for (y = 0; y < numRows; y++) {
      profile->data[y] = 0.0;
    }
    
    for (y = yStart; y < yStart+numRows; y++) {
      /* 
         We sort the spectral range we want to average, then we add only the
         first 80% of the data. This eleminates cosmix???
         Is this clever or very stupid??? 
         Note that outData in corrupted now (sorting is inplace)
      */
      sort(pixEnd-pixStart+1, outData->data+outXlen*y+pixStart);
      /**************
      FIXME: This way of removing cosmics, done in the original version by 
      T.O. does NOT work. To find objects we need to remove cosmic rays
      properly. I did it for testing in the recipe VmMosObsSingle using 
      VmImCosmics and it worked. However the cosmic cleaning must be 
      done in this function, and not in the recipe. 

      lastPix = pixStart 
      ************/

      lastPix = pixStart + specFrac*(pixEnd-pixStart); 
      for (x = pixStart; x <= lastPix; x++) {

        profile->data[y-yStart] += outData->data[x + outXlen*y];
      }
      /* average */
      
      profile->data[y-yStart] /= 1.0*(lastPix-pixStart+1);
      /* if(extractionSlit->slitNo == 12) printf("%d profile=%f\n", y-yStart,
	 profile->data[y-yStart]);*/
    }

/*      printf("EXTRSLIT no =%d\n", extractionSlit->slitNo); */
    windowObject = findObjectsInProfile(profile, detLevel, numLev, objFrac,
	   limFrac, minObjectSize, minCompositeSize);
    
    /* create windowslit */
    windowSlit = newWindowSlit();
    /* put objects in slit */
    if (windowObject != NULL)
      windowSlit->objs = windowObject;
    /* copy the sequential number of the slit */
    windowSlit->slitNo = extractionSlit->slitNo;
    /* where spectrum is in data file */
    windowSlit->specStart = yStart-limMargin;
    windowSlit->specEnd = yStart + numRows +limMargin -1;
    /* calculate whether slit is long or short */
    windowSlit->specLong = slitLongOrShort(extractionSlit, slitTolerance);
    
    /* fill how many objects in this slit */
    if (windowObject != NULL)
      windowSlit->numObj = numObjectsInWindowObject(windowObject);
    else {
      windowSlit->numObj = 0;
      cpl_msg_info(modName, "No object found in slit %d\n", 
		 extractionSlit->slitNo);
    }

    while (windowObject) { 
      objPos = windowObject->objPos;
      windowObject->objStart += limMargin;
      windowObject->objEnd += limMargin;
      windowObject->objPos += limMargin;
      windowObject->objX = extractionSlit->maskX->data[objPos];
      windowObject->objY = extractionSlit->maskY->data[objPos];
      windowObject = windowObject->next; 
    }
      
    /* link Window Slit to list */
    if (lastWindowSlit == NULL) {
      /* first slit we do */
      windowTable->slits = windowSlit;
    }
    else {
      /* not first, add to tail */
      lastWindowSlit->next = windowSlit;
      windowSlit->prev = lastWindowSlit;
    }
    /* store tail of linked list */
    lastWindowSlit = windowSlit;
    
    /* do next slit */
    extractionSlit = extractionSlit->next;
    deleteFloatArray(profile); 
  }
  if (flag) deleteImage(outData);

  /*****Test: write window table *******/
  extractionSlit = extractionTable->slits;
  windowSlit = windowTable->slits;
  while ( (extractionSlit) && (windowSlit) ) {
    /*  printf("EXTslitno =%d, WINslitnp =%d, numObj=%d\n", extractionSlit->slitNo,  */
    /*    windowSlit->slitNo,windowSlit->numObj ); */
    extractionSlit=extractionSlit->next;
    windowSlit = windowSlit->next;
  }
/*    nullImage = newImage(0, 0, (float *) NULL); */
/*    openNewFitsImage("windowTable.fits", nullImage); */
/*    writeFitsWindowTable(windowTable, nullImage->fptr); */
/*    closeFitsImage(nullImage,1); */
  
  return windowTable;
}


VimosImage **
VmSpEx2D(VimosImage **outSpSkyFra, VimosExtractionTable *extractionTable, 
         VimosSpecSampleOption sampleOption)
{
/*  int imageXlen, imageYlen; */
  int outXlen, outYlen;
  float wLenStart;
  float wLenEnd;
  float wLenInc;
  float tmpFloat;
  
  VimosImage  *imageData;
  VimosImage  *outData;
  VimosImage  *skyData;
  VimosImage  *outSkyData;
  VimosImage **outputImages;  
  VimosExtractionSlit *slit;

  char modName[] = "VmSpEx2D";

  cpl_msg_debug (modName, "2D extract spectra");
  
  imageData = outSpSkyFra[0];
/*  imageXlen = imageData->xlen;
  imageYlen = imageData->ylen; */

  skyData = outSpSkyFra[1];
  
  sampleOption = VM_SP_LIN_LAMBDA;
  
  slit = extractionTable->slits;

  /* get required length of spectrum and wavelength step */
  readFloatDescriptor(extractionTable->descs, pilTrnGetKeyword("WlenStart"), 
                      &wLenStart, NULL);
  readFloatDescriptor(extractionTable->descs, pilTrnGetKeyword("WlenEnd"), 
                      &wLenEnd, NULL);
  readFloatDescriptor(extractionTable->descs, pilTrnGetKeyword("WlenInc"), 
                      &wLenInc, NULL);
  /* checks (likely redundant) on params */
  if (wLenStart > wLenEnd) {
    tmpFloat  = wLenStart;
    wLenStart = wLenEnd;
    wLenEnd   = tmpFloat;
  }
  /* you never know... */
  wLenInc = fabs(wLenInc);  
  /* compute required length */
  outXlen = (wLenEnd-wLenStart)/wLenInc + 1;
  
  /* determine Y-size of output image */
  slit = extractionTable->slits;
  outYlen = 0;
  /* loop through slits and count number of output rows needed */
  while (slit) {
    outYlen += slit->numRows;
    slit = slit->next;
  }
  outData = newImageAndAlloc(outXlen, outYlen);
  outSkyData = newImageAndAlloc(outXlen, outYlen);
  /*  openNewMidasImage(outDataName, outData); */

  specExtract2D(imageData, extractionTable, outData, sampleOption, "default");
  specExtract2D(skyData, extractionTable, outSkyData, sampleOption, "default");


  copyAllDescriptors(imageData->descs, &(outData->descs)); 
  writeIntDescriptor(&(outData->descs),
		     pilTrnGetKeyword("Naxis",1) , outXlen, "");
  writeIntDescriptor(&(outData->descs),
		     pilTrnGetKeyword("Naxis",2) , outYlen, ""); 
  writeDoubleDescriptor(&(outData->descs), pilTrnGetKeyword("Crval",1), 
		       (double)(wLenStart), ""); 
  writeDoubleDescriptor(&(outData->descs), pilTrnGetKeyword("Crval",2),
                        1., "");
  writeDoubleDescriptor(&(outData->descs), pilTrnGetKeyword("Crpix",1),
                        1., "");
  writeDoubleDescriptor(&(outData->descs), pilTrnGetKeyword("Crpix",2),
                        1., "");
  writeDoubleDescriptor(&(outData->descs), pilTrnGetKeyword("Cdelt",1),
		       (double)(wLenInc), "");
  writeFloatDescriptor(&(outData->descs), pilTrnGetKeyword("Cdelt",2), 1., "");
  writeStringDescriptor(&(outData->descs), pilTrnGetKeyword("Ctype",1), 
			 "LAMBDA", "");
  writeStringDescriptor(&(outData->descs), pilTrnGetKeyword("Ctype",2), 
			"PIXEL", "");
  writeStringDescriptor(&(outData->descs), "ESO PRO VMTYPE", "STACKED", "");
 
  copyAllDescriptors(imageData->descs, &(outSkyData->descs));
  writeIntDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Naxis",1), 
		     outXlen, "");
  writeIntDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Naxis",2), 
		     outYlen, "");
  writeDoubleDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Crval",1), 
		       (double)(wLenStart), "");
  writeDoubleDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Crval",2),
                        1.,"");
  writeDoubleDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Crpix",1),
                        1.,"");
  writeDoubleDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Crpix",2),
                        1.,"");
  writeDoubleDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Cdelt",1), 
		       (double)(wLenInc), "");
  writeDoubleDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Cdelt",2),
                        1.,"");
  writeStringDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Ctype",1), 
			 "LAMBDA", "");
  writeStringDescriptor(&(outSkyData->descs), pilTrnGetKeyword("Ctype",2), 
			"PIXEL", "");
  writeStringDescriptor(&(outSkyData->descs), "ESO PRO VMTYPE", "STACKED", "");

  outputImages = (VimosImage **) cpl_malloc(2*sizeof(VimosImage *));
  outputImages[0] = outData;
  outputImages[1] = outSkyData;

  return outputImages;
}


VimosImage **
VmSpEx1D(VimosImage **outSpSkyEx, VimosWindowTable *windowTable,
         VimosObjectTable *objectTable, int horne, int ncomb)
{
  int imageXlen, imageYlen; 
  int outXlen, outYlen;
  
  VimosImage  *imageData;
  VimosImage  *skyImage;
  VimosImage  *outData;
  VimosImage  *outSkyData;
  VimosImage **outputs;

  char modName[] = "VmSpEx1D";

  cpl_msg_info (modName, "1D extract Spectra");

  imageData = outSpSkyEx[0];
/*   imageXlen = imageData->xlen; */
  imageYlen = imageData->ylen;
     
  /* same for sky spectra */
  skyImage = outSpSkyEx[1];

  /* output X size is the same as input */
  outXlen = imageXlen;
  /* Y size is numberof objects defined in Window Table */
  outYlen = numObjsInWindowTable(windowTable);

  /* create output image */
  outData = newImageAndAlloc(outXlen, outYlen);
  outSkyData = newImageAndAlloc(outXlen, outYlen);
  
  /* copy all descriptors from Window Table to Object Table*/
  copyWinTab2ObjTab(windowTable, objectTable);

  /* do extraction */
  specExtract1D(imageData, skyImage, windowTable, outData, 
		outSkyData, objectTable, horne, ncomb);

  copyAllDescriptors(imageData->descs, &(outData->descs));
  writeIntDescriptor(&(outData->descs), "NAXIS1", outXlen, "");
  writeIntDescriptor(&(outData->descs), "NAXIS2", outYlen, "");
  writeStringDescriptor(&(outData->descs), "CTYPE2", "OBJECT", "");
  writeStringDescriptor(&(outData->descs), "ESO PRO VMTYPE", "PACKED", "");

  copyAllDescriptors(imageData->descs, &(outSkyData->descs));
  writeIntDescriptor(&(outSkyData->descs), "NAXIS1", outXlen, "");
  writeIntDescriptor(&(outSkyData->descs), "NAXIS2", outYlen, "");
  writeStringDescriptor(&(outSkyData->descs), "CTYPE2", "SKY", "");
  writeStringDescriptor(&(outSkyData->descs), "ESO PRO VMTYPE", "PACKED", "");

  outputs = (VimosImage **) cpl_malloc(2*sizeof(VimosImage *));
  outputs[0] = outData;
  outputs[1] = outSkyData;
  
  return outputs;
 
}


VimosImage *
VmSpStack2D(VimosImage **images2D, VimosWindowTable **winTables, 
            VimosExtractionTable *extTable, VimosIntArray *twoDMap, 
            int frameCount, CombMethod combMethod, 
            CombParameters *combParameter, VimosFloatArray *offSets,
            int offSetsFlag, int mapFlag, int rescaleFlag)
{
  char        task[]               = "VmSpStack2D";


  int                   wSlitRefNo       = 0;
  int                   foundOffSet      = 0;
  int                   foundOffSetOk    = 0;
  int                   i,j,ok;
  int                   slitXlen;
  int                   newSlitYlen,slitYlen;
  int                   xlen2D,ylen2D;
  int                   startY; 
  int                   ref,medxlen;
  int                   nRow32low,nRow32high;

  float                 **singleOffSets,**singleOffSetsOk;
/*  float                 minOffSet,maxOffSet; */
  float                 maxWidth,objPos;
  float                 min,refOffSet,refExpectedOffSets;
  float                 max,newMax,yadd;
  float                 levelMedian, rescaleFactor;

  double                refValue,dValue,pixScale;
  double                coordOffSet;
  

  VimosImage           *slitCombImage    = NULL;
  VimosImage           *outImage         = NULL;

  VimosFloatArray      *expectedOffSets  = NULL;
  VimosFloatArray      *imageMedian      = NULL;
  VimosIntArray        *tmpY             = NULL;

  VimosImage           **slitImages       = NULL;
  VimosImage           **shiftSlitImages  = NULL;

  VimosWindowSlit      *wSlitRef         = NULL;

  VimosExtractionSlit  *extrSlit         = NULL;

  VimosWindowObject    *object           = NULL;

  VimosImage           *tmpImage         = NULL;

  /* Error handling */
  
  if (images2D == NULL) {
    cpl_msg_debug(task, "NULL input images");
    return(NULL);
  }

  if (winTables == NULL) {
    cpl_msg_debug(task, "NULL input tables");
    return(NULL);
  }


  cpl_msg_debug(task, "Stacking 2D images");
 

  /* Reference slit number */
  wSlitRef=winTables[0]->slits;
  wSlitRefNo=0;
  while(wSlitRef) {
    wSlitRefNo++;
    wSlitRef=wSlitRef->next;
  }



  /*
   *  Offset determination (optional)                  
   */

  if(offSetsFlag) {

    /* Determine expected offsets */
    expectedOffSets = newFloatArray(frameCount);
    readDoubleDescriptor(images2D[0]->descs, pilTrnGetKeyword("PixelScale"),
                         &pixScale, NULL);
    readDoubleDescriptor(images2D[0]->descs, pilTrnGetKeyword("Alpha"),
                         &refValue, NULL); 
    for(i=1; i<frameCount; i++) {    
      readDoubleDescriptor(images2D[i]->descs, pilTrnGetKeyword("Alpha"),
                           &dValue, NULL); 
      coordOffSet = refValue - dValue; 
      expectedOffSets->data[i] = -1 * (coordOffSet * 3600 /  pixScale);  
    }
    
    
    /* Allocation */
    singleOffSets = (float **) cpl_calloc((frameCount-1),sizeof(float *));
    singleOffSetsOk = (float **) cpl_calloc((frameCount-1),sizeof(float *));
    for(i=0; i<frameCount-1; i++) {
      singleOffSets[i] = (float *) cpl_calloc(wSlitRefNo+1,sizeof(float));
      singleOffSetsOk[i] = (float *) cpl_calloc(wSlitRefNo+1,sizeof(float));
    }
    
    
    wSlitRef=winTables[0]->slits;
    foundOffSet=0;
    
    while(wSlitRef) {
      
      /* if an object doesn't exists in the reference slit go to next object in all slits */
      if(!wSlitRef->objs){
	if(wSlitRef->next) {
	  wSlitRef=wSlitRef->next;
	  for(i=1; i<frameCount; i++) winTables[i]->slits=winTables[i]->slits->next;
	  continue;
	} else break; /* last slit */
      }
      
      /* check in the object found in the reference slit exists in all slits */
      ok=1;
      for(i=1; i<frameCount; i++) if(!winTables[i]->slits->objs) ok=0;
      if(ok) foundOffSet++;
      
      /* compute offsets (if ok) and go to next slit in all tables */
      for(i=1; i<frameCount; i++) {
	if(ok) {
	  /* find the most brillant object */
	  object = winTables[i]->slits->objs;
	  maxWidth = object->objWidth;
	  objPos = object->objPos;
	  object = object->next;
	  while(object) {
	    if(object->objWidth > maxWidth){
	      maxWidth = object->objWidth;
	      objPos = object->objPos;
	    }
	    object = object->next;
	  }
	  /* determine shift */
	  singleOffSets[i-1][foundOffSet-1] = wSlitRef->objs->objPos - objPos;
	}
	if(winTables[i]->slits->next) winTables[i]->slits=winTables[i]->slits->next;
      }
      wSlitRef=wSlitRef->next;
    }

    /* rewind windowSlits */
    for(i=1; i<frameCount; i++) {
      while(winTables[i]->slits->prev) winTables[i]->slits=winTables[i]->slits->prev;
    }
    
    
    /* compute offset median rejecting bad values  */
    
    offSets->data[0] = 0.00;
    for(i=0; i<frameCount-1; i++) {
      foundOffSetOk=0;
/*      minOffSet = expectedOffSets->data[i+1] - OFFSET_TOLERANCE;
      maxOffSet = expectedOffSets->data[i+1] + OFFSET_TOLERANCE; */
      
      /* rejecting */
      for(j=0; j<foundOffSet; j++) {
/* if(singleOffSets[i][j]>= minOffSet && singleOffSets[i][j]<= maxOffSet){ */
	  singleOffSetsOk[i][foundOffSetOk] = singleOffSets[i][j];
	  foundOffSetOk++;
/* } */
      }

      cpl_msg_debug(task,"Image %d: rejected %d offsets on %d", i+1, 
		 foundOffSet - foundOffSetOk, foundOffSet);

      if(foundOffSetOk>=MIN_OFFSET_OK){ 
      	offSets->data[i+1] = median(singleOffSetsOk[i],foundOffSetOk);
      } else {
      	cpl_msg_error(task,"Image %d: too few available offsets", i+1);
      	return NULL;
      }

    }


    /* determine the image to be used as reference in shifting */
    min = offSets->data[0];
    ref = 0;
    for(i=0; i<frameCount; i++) {
      if (offSets->data[i] < min) {
	min = offSets->data[i];
	ref = i;
      }
    }

    /* rescale offsets  */
    refOffSet = offSets->data[ref];
    refExpectedOffSets = expectedOffSets->data[ref];
    for(i=0; i<frameCount; i++) {
      offSets->data[i] -= refOffSet;
      expectedOffSets->data[i] -= refExpectedOffSets;
    } 

    cpl_msg_debug(task,"Image %d is the reference image",ref+1);
    for(i=0; i<frameCount; i++) {
      if(i==ref) continue;
      else {
	cpl_msg_debug(task,"Image %d shifts:",i+1);
	cpl_msg_debug(task,"Expected: %f Computed: %f Diff: %f",
		   expectedOffSets->data[i],offSets->data[i], 
                   expectedOffSets->data[i] - offSets->data[i]);
      }
    }

  } else {
    /* the reference image is the one with offSet=0 */
    for(i=0; i<frameCount; i++) {
      if(offSets->data[i]==0.00) {
	ref=i;
	break;
      }
    }
  }


  /*
   *  Images shifting and combining
   */
  

  wSlitRef=winTables[ref]->slits; /* shift using as reference the ref image */
  if(extTable) extrSlit=extTable->slits;

  startY = 0;

  /* determine output 2D image dimension */
  max = offSets->data[0];
  for(i=0; i<frameCount; i++) { 
    newMax = offSets->data[i];
    if(newMax > max) max=newMax;
  }
  yadd = (int)(max) + 1;
  xlen2D = slitXlen = images2D[0]->xlen;
  ylen2D = images2D[0]->ylen + wSlitRefNo*yadd;
  medxlen = xlen2D/2;


  /* allocations */
  if(!(outImage=newImageAndAlloc(xlen2D,ylen2D))){
    cpl_msg_error(task, "Failure in creating output 2D image");
    return(NULL);
  }

  copyAllDescriptors(images2D[0]->descs, &(outImage->descs));

/*  outImage->descs=images2D[0]->descs;  ?!?!?!? */

  /* change naxis2 value */
  writeIntDescriptor(&(outImage->descs),"NAXIS2",ylen2D,"# pixels/axis");


  if (!(slitImages = (VimosImage **)cpl_calloc(frameCount, sizeof(VimosImage *)))) {
    cpl_msg_error(task, "Failure creating list of temporary slit images");
    return NULL;
  }

  if (!(shiftSlitImages = (VimosImage **)cpl_calloc(frameCount, sizeof(VimosImage *)))) {
    cpl_msg_error(task, "Failure creating list of temporary shifted slit images");
    return NULL;
  }

  
  while(wSlitRef) {
    
    /* extract slit images */
    for(i=0; i<frameCount; i++) {
      slitYlen = wSlitRef->specEnd - wSlitRef->specStart + 1;
      tmpImage = newImage(slitXlen,slitYlen,extractFloatImage(images2D[i]->data,
                 xlen2D,ylen2D,0,wSlitRef->specStart,slitXlen,slitYlen));
      slitImages[i] = tmpImage;    
    }    
    
    /* images shifting */
    newSlitYlen = slitYlen + yadd;
    for(i=0; i<frameCount; i++) {
      shiftSlitImages[i] =  imageShift(slitImages[i],0.,offSets->data[i],slitXlen,newSlitYlen,-32000.);
    }      

    /* median image level for rescaling */
    if (rescaleFlag) {
      imageMedian = newFloatArray(frameCount);

      for(i=0; i<frameCount; i++) {
	nRow32low = ceil(offSets->data[i])+2;
	nRow32high = ceil(yadd - offSets->data[i])+2;
	imageMedian->data[i] = medianPixelvalue((shiftSlitImages[i])->data+nRow32low*slitXlen+1, (newSlitYlen-nRow32low-nRow32high)*slitXlen);
      }      
      levelMedian = medianPixelvalue(imageMedian->data, frameCount);
      for(i=0; i<frameCount; i++) {
	rescaleFactor = imageMedian->data[i] / levelMedian;
	for (j=0; j < slitXlen*newSlitYlen; j++) {
	  if (shiftSlitImages[i]->data[j]>-31000.)
	    shiftSlitImages[i]->data[j] /= rescaleFactor;
	}
      }
    }

    /* slit images combining  */
    if(!(slitCombImage = frComb32000(shiftSlitImages, frameCount, combMethod, 
                         combParameter, 1))){
      cpl_msg_error(task, "Failure combining images for slit %d",wSlitRef->slitNo);
      return(NULL);
    }

    /* 2D map */
    if(mapFlag) {
      for(j=startY; j<startY+newSlitYlen; j++){
	for(i=0; i<frameCount; i++) {
	  if(shiftSlitImages[i]->data[medxlen + xlen2D * (j-startY)] != -32000) 
	    twoDMap->data[j]++;
	}
      }
    }
      
    /* insert slitCombImage in output image  */
    if(VM_FALSE == insertFloatImage(outImage->data, xlen2D, ylen2D, 
				    0, startY, slitXlen, newSlitYlen, slitCombImage->data)){
      cpl_msg_error(task, "Failure inserting combined image of slit %d",wSlitRef->slitNo);
      return(NULL);
    }

    /* update the extraction table (optional) */
    if (offSetsFlag && extrSlit) {

      tmpY = newIntArray(newSlitYlen);
      for (j = 0; j < newSlitYlen; j++) 
        tmpY->data[j] = j + startY + 1;
      deleteIntArray(extrSlit->y);
      extrSlit->y = tmpY;

      /*
       * Added, to avoid a segfault when the destructor is called...
       * Models are not needed anymore anyway.
       */

      if (extrSlit->crvPol && extrSlit->invDis) {
        for (j = 0; j < extrSlit->numRows; j++) {
          deleteDistModel1D(extrSlit->crvPol[j]);
          deleteDistModel1D(extrSlit->invDis[j]);
        }
        cpl_free(extrSlit->crvPol);
        extrSlit->crvPol = NULL;
        cpl_free(extrSlit->invDis);
        extrSlit->invDis = NULL;
      }

      extrSlit->numRows = newSlitYlen;

    }

    /* go to next slit */
    startY+= newSlitYlen;
    wSlitRef=wSlitRef->next;
    if(extrSlit) extrSlit=extrSlit->next;
  }
    
  return outImage; 

}


/**
 * @memo
 *   Construct an Extraction Table from ADF and Grism (or IFU) table.
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param inputImage      Image having ADF in header.
 * @param grismTable      Input Grism Table
 * @param ifuTable        Input IFU Table
 * @param extrTable       Input ExtractionTable (optional)
 *
 * @doc
 *   Just fit a polynomial in two variables to the coefficients
 *   of the polynomials fitting the edges of spectra found in a
 *   spectral Flat Field
 *
 * @author M. Scodeggio
 */

VimosExtractionTable *
VmSpExTab(VimosImage *inputImage, VimosTable *grismTable,
          VimosIfuTable *ifuTable, VimosExtractionTable *extrTable)
{
  char                  modName[] = "VmSpExTab";
  VimosTable           *adf;
  VimosAdfType          adfType;
  VimosExtractionTable *extractionTable=NULL;
  
  cpl_msg_info(modName, "Building Extraction Table from ADF");

  adfType = VM_ADF_TYPE_UDF;

  adf = newADF();

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
  
 /* 
  *  Create Extraction Table 
  */
  if (extrTable)
  {
    extractionTable = extrTable;
    deleteExtractionSlit(extractionTable->slits);
  }
  else
  {
    extractionTable = newExtractionTable(); 

    /*
     * Copy grism and filter name to the extraction table.
     */

    vimosDscCopy(&extractionTable->descs, inputImage->descs,
                 "^ESO INS (FILT|GRIS)[1-4] (NAME|ID)",
                 pilTrnGetKeyword("Table"));

    /*
     * Get just a few descriptors from the Grism Table.
     */

    if ((copyGrsTab2ExtTab(grismTable, extractionTable) == VM_FALSE)) {
      cpl_msg_error(modName, 
                  "Failure in copying Grism Table to Extraction Table");
      return NULL;
    }

    /*
     * The distorsion models are read from the image header
     * (curiously the ADF Table contains just a copy of the
     * input image header).
     */
    
    if ((copyAdf2ExtTab(adf, extractionTable) == VM_FALSE)) {
      cpl_msg_error(modName, 
                  "Failure in copying ADF info to Extraction Table");
      return NULL;
    }
  }

  /* 
   *  Use ADF to compute extraction related stuff 
   */

  if ((computeExtractionTable(adf, ifuTable, extractionTable) == VM_FALSE))
  {
    cpl_msg_error(modName, "Failure in computing Extraction Table");
    return NULL;
  }

  return extractionTable;
}

/*
 * Images containing 2D-extracted spectra and belonging to the same
 * mask are compared, the smaller images are expanded to the largest
 * one, adding the appropriate number of extra rows to each spectrum.
 * (C.Izzo, May 16, 2005).
 */

int VmSpMatch(VimosImage **image, int n, VimosExtractionTable *extTable,
              VimosTable **winList, int w)
{

  char modName[] = "VmSpMatch";

  int gap       = 0;            /* True if in a gap between spectra        */
  int row       = 0;            /* Counter of spectra rows                 */
  int maxy      = 0;            /* Y size of largest input image           */
  int maxima    = 0;            /* Sequence number of largest input image  */
  int countSlit = 0;            /* Slit counter                            */
  int numRows;                  /* Number of rows in a slit spectrum       */
  int delta;                    /* Difference in Y size for each spectrum  */
  int found;
  int winLength;
  int xlen      = image[0]->xlen;
  int numSlit   = numSlitsInExtTable(extTable);
  int i, j, k, m, mm;

  int    *maxSlitSize;
  int    *slitSize;
  float  *data;
  int    *slitEnd;
  int    *slitNumber;
  int    *gapPos;

  VimosIntArray       *y;
  VimosExtractionSlit *slit;
  VimosTable          *winTable;

 /* char recipeStart[80], refRecipeStart[80]; */
  char creationDate[80], refCreationDate[80];


  /*
   * Check that all images have the same X size
   */

  for (i = 1; i < n; i++) {
    if (xlen != image[i]->xlen) {
      cpl_msg_error(modName, "Incompatible images: different sizes in x");
      return EXIT_FAILURE;
    }
  }


  /*
   * Find among input images the ones that have no -32000 gaps,
   * and insert such gaps according to the ADF content.
   */

  for (i = 0; i < n; i++) {
    gap = 0;
    for (j = 0, k = xlen / 2; j < image[i]->ylen; j++, k += xlen) {
      if (image[i]->data[k] < -30000) {
        gap = 1;
        break;
      }
    }

    if (gap)
      continue;

    /*
     * This is an image without -32000 gaps. This means that we need
     * a window table to know where the slit spectra are. The right
     * window table is the one having the same ESO PRO REC1 START
     * descriptor - i.e., the one produced in the same recipe run
     * that produced the 2D-extracted image.
     */

/*** Removed part based on same ESO PRO REC1 START instead of same DATE

    found = 0;
    if (readStringDescriptor(image[i]->descs, 
                             pilTrnGetKeyword("RecipeStart", 1),
                             refRecipeStart, NULL) == VM_FALSE) {
      cpl_msg_error(modName, "Missing descriptor %s in image header!",
                    pilTrnGetKeyword("RecipeStart", 1));
      return EXIT_FAILURE;
    }

    for (j = 0; j < w; j++) {
      if (readStringDescriptor(winList[j]->descs, 
                               pilTrnGetKeyword("RecipeStart", 1),
                               recipeStart, NULL) == VM_FALSE) {
        cpl_msg_error(modName, "Missing descriptor %s in window table header!",
                      pilTrnGetKeyword("RecipeStart", 1));
        return EXIT_FAILURE;
      }
      if (!strncmp(recipeStart, refRecipeStart, 23)) {
        found = 1;
        winTable = winList[j];
        break;
      }
    }

    if (!found) {
      cpl_msg_error(modName, "A matching window table was not found for "
                    "the spectral image produced by the vmmosobsstare "
                    "recipe at %s! Please specify the corresponding "
                    "window table.", refRecipeStart);
      return EXIT_FAILURE;
    }

End of removed part ***/

/*** New part based on very close (few seconds) creation times */

    found = 0;
    if (readStringDescriptor(image[i]->descs, "DATE",
                             refCreationDate, NULL) == VM_FALSE) {
      cpl_msg_error(modName, "Missing keyword DATE in image header!");
      return EXIT_FAILURE;
    }

    for (j = 0; j < w; j++) {
      if (readStringDescriptor(winList[j]->descs, "DATE",
                               creationDate, NULL) == VM_FALSE) {
        cpl_msg_error(modName, "Missing keyword DATE in window table header!");
        return EXIT_FAILURE;
      }
      if (abs(diffdate(creationDate, refCreationDate)) < 5) {
        found = 1;
        winTable = winList[j];
        break;
      }
    }

    if (!found) {
      cpl_msg_error(modName, "A matching window table was not found for "
                    "the spectral image produced by the vmmosobsstare "
                    "recipe at (about) DATE = %s! Please specify the "
                    "corresponding window table.", refCreationDate);
      return EXIT_FAILURE;
    }

/*** End of new part ***/

    winLength = tblGetSize(winTable, "SLIT");
    slitNumber = tblGetIntData(winTable, "SLIT");
    slitEnd = tblGetIntData(winTable, "SPEC_END");

    gapPos = (int *)cpl_calloc(numSlit, sizeof(int));

    countSlit = 0;
    for (j = 0; j < winLength; j++) {
      if (j)
        if (slitEnd[j] == slitEnd[j-1])
          continue;
      gapPos[countSlit] = slitEnd[j];
      countSlit++;
    }


    /*
     * Inserting the gaps!
     */

    data = (float *)cpl_calloc(image[i]->xlen * (image[i]->ylen + 3 * numSlit),
                               sizeof(float));

    row = 0;
    for (m = 0; m < 2; m++) {    /* Pad with -32000 first 2 rows */
      for (mm = 0; mm < xlen; mm++)
        data[row*xlen + mm] = -32000;
      row++; 
    }

    countSlit = 0;
    for (j = 0; j < image[i]->ylen; j++) {
      for (mm = 0; mm < xlen; mm++)
        data[row*xlen + mm] = image[i]->data[j*xlen + mm];
      row++;
      if (j == gapPos[countSlit]) {
        countSlit++;
        if (countSlit < numSlit) {
          for (m = 0; m < 3; m++) {    /* Pad with -32000 three rows */
            for (mm = 0; mm < xlen; mm++)
              data[row*xlen + mm] = -32000; 
            row++; 
          }
        }
      }
    }

    for (mm = 0; mm < xlen; mm++)      /* Pad last row with -32000 */
      data[row*xlen + mm] = -32000;

    cpl_free(gapPos);

    cpl_free(image[i]->data);
    image[i]->data = data;
    image[i]->ylen += 3 * numSlit;

/***
createFitsImage("staccata.fits", image[i], "test");
***/


  }


  /*
   * Find the largest input image
   */

  for (i = 0; i < n; i++) {
    if (maxy < image[i]->ylen) {
      maxy = image[i]->ylen;
      maxima = i;
    }
  }


  /*
   * Modify the extraction table to describe the slits in the largest
   * image, not counting the -32000 gaps (that will be suppressed).
   */

  maxSlitSize = (int *)cpl_calloc(numSlit, sizeof(int));
  slit = extTable->slits;
  countSlit = 0;
  row = 0;
  gap = 1;

  for (i = 0, j = xlen / 2; i < maxy; i++, j += xlen) {
    if (image[maxima]->data[j] < -30000) {
      if (!gap) {
        gap = 1;

        if (!slit) {
          cpl_msg_error(modName, "Found in image more slits than in ADF!");
          cpl_free(maxSlitSize);
          return EXIT_FAILURE;
        }

        y = newIntArray(numRows);
        for (k = 0; k < numRows; k++)
          y->data[k] = k + row + 1;

        deleteIntArray(slit->y);
        slit->y = y;

        if (slit->crvPol && slit->invDis) {
          for (k = 0; k < slit->numRows; k++) {
            deleteDistModel1D(slit->crvPol[k]);
            deleteDistModel1D(slit->invDis[k]);
          }
          cpl_free(slit->crvPol);
          slit->crvPol = NULL;
          cpl_free(slit->invDis);
          slit->invDis = NULL;
        }

        slit->numRows = numRows;
        slit = slit->next;
        maxSlitSize[countSlit] = numRows;
        countSlit++;
        row += numRows;
      }
    }
    else {
      if (gap) {
        gap = 0;
        numRows = 0;
      }
      numRows++;
    }
  }

  if (!gap) {

    if (!slit) {
      cpl_msg_error(modName, "Found in image one more slits than in ADF!");
      cpl_free(maxSlitSize);
      return EXIT_FAILURE;
    }

    y = newIntArray(numRows);
    for (k = 0; k < numRows; k++)
      y->data[k] = k + row + 1;

    deleteIntArray(slit->y);
    slit->y = y;

    if (slit->crvPol && slit->invDis) {
      for (k = 0; k < slit->numRows; k++) {
        deleteDistModel1D(slit->crvPol[k]);
        deleteDistModel1D(slit->invDis[k]);
      }
      cpl_free(slit->crvPol);
      slit->crvPol = NULL;
      cpl_free(slit->invDis);
      slit->invDis = NULL;
    }

    slit->numRows = numRows;
    slit = slit->next;
    maxSlitSize[countSlit] = numRows;
    countSlit++;

  }

  if (slit) {
    cpl_msg_error(modName, "Found in image less slits than in ADF!");
    cpl_free(maxSlitSize);
    return EXIT_FAILURE;
  }


  /*
   * Map all images into new data buffers, where the -32000 gaps
   * are eliminated, and the spectra are expanded to the size of
   * the largest image.
   */

  for (i = 0; i < n; i++) {

    /*
     * Allocate the output image: this will have the same size of
     * the largest input image, minus the -32000 gaps. For each slit
     * there is a 3 pixels wide gap.
     */

    data = (float *)cpl_calloc(image[i]->xlen * (maxy - 3 * numSlit), 
                               sizeof(float));
    slitSize = (int *)cpl_calloc(numSlit, sizeof(int));

    countSlit = 0;
    gap = 1;

    for (j = 0, k = xlen / 2; j < image[i]->ylen; j++, k += xlen) {
      if (image[i]->data[k] < -30000) {
        if (!gap) {
          gap = 1;
          slitSize[countSlit] = numRows;
          countSlit++;
        }
      }
      else {
        if (gap) {
          gap = 0;
          numRows = 0;
        }
        numRows++;
      }
    }

    countSlit = 0;
    row = 0;
    gap = 1;

    for (j = 0, k = xlen / 2; j < image[i]->ylen; j++, k += xlen) {
      if (image[i]->data[k] < -30000)
        gap = 1;
      else {
        if (gap) {
          gap = 0;
          delta = maxSlitSize[countSlit] - slitSize[countSlit];
          countSlit++;
          for (m = 0; m < delta; m++) {    /* Pad with zeroes */
            for (mm = 0; mm < xlen; mm++)
              data[row*xlen + mm] = 0;
            row++;
          }
        }
        if (delta < 0) {
          /*
           * It may rarely happen that even if one image is larger than
           * the other, a particular slit is larger in the smaller image.
           * In this case, we suppress some data lines to recover self-
           * consistency. This case may happen only with images having
           * equal or very close sizes.
           */
          delta++;
          continue;
        }
        for (mm = 0; mm < xlen; mm++)
          data[row*xlen + mm] = image[i]->data[j*xlen + mm];
        row++;
      }
    }

    cpl_free(image[i]->data);
    image[i]->data = data;
    image[i]->ylen = maxy - 3 * numSlit;

  }

  return EXIT_SUCCESS;

}

/**@}*/
