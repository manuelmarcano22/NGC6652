/* $Id: vmifufibers.c,v 1.3 2013-07-11 11:27:43 cgarcia Exp $
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
 * $Date: 2013-07-11 11:27:43 $
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

#include <pilmessages.h>
#include <piltranslator.h>

#include "vmimage.h"
#include "vmimageset.h"
#include "vmmatrix.h"
#include "vmtable.h"
#include "vmifutable.h"
#include "vmgrismtable.h"
#include "vmextractiontable.h"
#include "vmobjecttable.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmifufibers.h"


#define PSEUDOSLIT 400
#define DO_LINE 1
#define DO_CONTINUUM 2
#define NTERMS 6
#define NPOINTS 13
#define HALFPOINTS ((NPOINTS -1)/2)

#define THE_SKY_LINE 5577.1

/* conversion factor between FWHM and SIGMA for a gaussian function */
#define FWHM_TO_SIGMA 2.35482

/* gaussians are computed in the range [-NSIGMA*sigma ; +NSIGMA*sigma] */
#define NSIGMA 10.0

/* parameter used to set the number of adjacent fibers (per side, except
   beginning and end of the pseudoslit, see code) that we consider are 
   crosstalk-contributing to the flux of a given fiber */
#define HOWMANYFIBS 2

/* noise (of the theoretical profiles) * LIMIT defines the threshold
   for cosmic rejection */
#define LIMIT 5.0


int
VmIfuComputePSF(VimosImage *imageData, VimosIfuTable *ifuTable, 
                VimosObjectTable *objectTable, char *toBeFitted, 
                int nIntervals, int quadNo)
{
  int i,j,k,l,ii, uu;
  int numSkyLines, skyPixInc, skyPix, skyLowPix;
  int Fit, numQuant, stepInt, totalGrouped;

  int nGoodFibs, nDeadFibs, nTotGood, nTotBad;
  int specLen;

  float *coeffs, skyLam, psfy, step;
  float XRefPixVal, YRefPixVal, XRefPix, YRefPix, lambdaInc, skyPixIncF;

  char modName[] = "VmIfuComputePSF";
  char comment[80];

  VimosIntArray *frequency, *numGrouped;
  VimosFloatArray *skyLines, *lineSigma, *psfYLimit, *anX, *anY;
  VimosFloatArray *objectPsfY, *sortObjectPsfY, *tmpObjectPsfY;

  VimosIfuQuad *ifuQuads;
  VimosIfuSlit *theIfuSlits;
  VimosIfuFiber *theIfuFibers;
  VimosObjectObject *objects;

  pilMsgInfo (modName, "Computing PSF groups");

  puts("WARNING! LAMBDA of SKY LINES in ANGSTROMS!");

  /* initialize this option to fit only one line */
  Fit = 0;

  if (!strncmp(toBeFitted,"ONE",3))
    {
     Fit = 0;
    }
  else if (!strncmp(toBeFitted,"ALL",3))
    {
     Fit = 1;
    }
  else
    {
     pilMsgError(modName, "Unable to set fit parameter");
     return EXIT_FAILURE;
    }

  switch(Fit)
    {
    default:
    case 0 :
      {
       numSkyLines = 1;
       skyLines = newFloatArray(1);
       skyLines->data[0] = THE_SKY_LINE;
       /* define the array with the sigma of each sky line, to be computed */
       lineSigma = newFloatArray(numSkyLines);
       break;
      }
    case 1 :
      {
       readSkyLines(objectTable->descs, &numSkyLines, &skyLines);
       /* define the array with the sigma of each sky line, to be computed */
       lineSigma = newFloatArray(numSkyLines);
       break;
      }
    }

  specLen = imageData->xlen;

  /* define how many quantiles must be computed for grouping fibers in */
  /* each pseudoslit */
  /* if > 5 intervals (4 quantiles) default to 5 intervals */
  if (nIntervals > 5) 
    {
     nIntervals = 5;
     pilMsgInfo(modName,"N. of intervals greater than 5, DEFAULTING to 5");
    }

  /* write keyword in IFU Table, will be used by VmIfuSky */
  writeIntDescriptor(&ifuTable->descs, "ESO PRO SKYGROUP",
                         nIntervals, "");

  /* By definition, (n-1) variate-values divide the total frequency in n */
  /* intervals */
  numQuant = nIntervals - 1;
  frequency = newIntArray(numQuant);

  /* initialize array with values corresponding to the selected quantiles */
  /* and set first value to zero: this is used for first interval */
  psfYLimit = newFloatArray(nIntervals);
  psfYLimit->data[0] = 0.;

  /* array for the total number of fibers grouped in each PSFY interval */
  numGrouped = newIntArray(nIntervals);

  /* setting frequency step. REMEMBER: frequency range from 1 to 400, this */ 
  /* is taken into account later when evaluating the quantile (when defining */
  /* the rank of the kth-smallest) */
  for (i=1; i<numQuant; i++)
    {
     frequency->data[i] = frequency->data[i-1] + stepInt;
    }

  /* read descriptors from the Object Table */
  /* FROM DRS DOCUMENT: */
  /* Value of reference pixel in X (packed-spectra format): nm */
  /*
  readFloatDescriptor(objectTable->descs, "CRVAL1", &XRefPixVal, comment);
  readFloatDescriptor(objectTable->descs, "CRVAL2", &YRefPixVal, comment);
  readFloatDescriptor(objectTable->descs, "CRPIX1", &XRefPix, comment);
  readFloatDescriptor(objectTable->descs, "CRPIX2", &YRefPix, comment);
  */
  /*  FROM DRS DOCUMENT: */
  /* WARNING: Wavelength step (packed-spectra format): Angstroms */
  /*
  readFloatDescriptor(objectTable->descs, "CDELT1", &lambdaInc, comment);
  */

  /* READ THESE DESCRIPTORS FROM IMAGE DATA, NOT OBJECT TABLE */
  /* FROM DRS DOCUMENT:
     Value of reference pixel in X (packed-spectra format): A */
  readFloatDescriptor(imageData->descs, "CRVAL1", &XRefPixVal, comment);
  readFloatDescriptor(imageData->descs, "CRVAL2", &YRefPixVal, comment);
  readFloatDescriptor(imageData->descs, "CRPIX1", &XRefPix, comment);
  readFloatDescriptor(objectTable->descs, "CRPIX2", &YRefPix, comment);

  /*  FROM DRS DOCUMENT: */
  /* WARNING: Wavelength step (packed-spectra format): Angstroms */
  readFloatDescriptor(imageData->descs, "CDELT1", &lambdaInc, comment);

  /* NOTE: VALUES IN HEADER ARE ALL IN ANGSTROMS (24/04/02)*/
  /* WARNING! TRANSFORM lambdaInc from angstroms to nm for consistency */
  /*  lambdaInc = lambdaInc/10.;*/


  /* initialize arrays that will be used to fit sky lines */
  anX=newFloatArray(NPOINTS);
  anY=newFloatArray(NPOINTS);
  coeffs=floatVector(1,NTERMS);

  /* initialize counters for total number of good fibers and dead fibers on */
  /* this image */
  nTotGood = 0;
  nTotBad = 0;

  /* start loop to compute fluxes of sky lines for all the 1600 spectra */
  /* loop on the IFU Table */
  ifuQuads = ifuTable->quads;
  while (ifuQuads)
    {
      /* take the pseudoslits in the quadrant this image refers to */
     if (ifuQuads->quadNo == quadNo)
       {
        theIfuSlits = ifuQuads->ifuSlits;
        while(theIfuSlits)
          {
	   /* counters for dead fibers and used fibers in this slit */
	   nDeadFibs = 0;
	   nGoodFibs = 0;

	   /* count how many good fibers (not dead ones) in this slit */
	   theIfuFibers = theIfuSlits->fibers;
	   while(theIfuFibers)
             {
	      if (theIfuFibers->fiberTrans != -1.) nGoodFibs +=1;
	      theIfuFibers = theIfuFibers->next;
	     }
	   /* initialize working arrays */
	   objectPsfY = newFloatArray(nGoodFibs);
	   sortObjectPsfY = newFloatArray(nGoodFibs);
	   tmpObjectPsfY = newFloatArray(nGoodFibs);

	   pilMsgInfo(modName,"%d fibers used for PSF grouping on this image",
	     nGoodFibs);


	   /* total frequency is no. of good fibers in the slit */
	   /* step in frequency for the quantiles */
	   step = ( (float)nGoodFibs ) / ( (float)nIntervals );

	   /* take integer part */
	   stepInt = (int)step ;

	   /* if step it is not integer, set frequency of first quantile */
	   /* to stepInt+1 and use stepInt to set the frequencies of other */
	   /* quantiles */
	   if ( (step-((float)stepInt)) != 0.0 ) frequency->data[0]=stepInt+1;
	   if ( (step-((float)stepInt)) == 0.0 ) frequency->data[0] = stepInt;


	   /* rewind */
	   theIfuFibers = theIfuSlits->fibers;

	   while(theIfuFibers)
             {
	      /* if it is not a dead fiber */
	      if ( (theIfuFibers->fiberTrans) != -1.) 
		{
		 /* for each fiber, start a loop on the objects */
		 objects = objectTable->objs;

		 /* initialize counter of objects in this IFU image */
		 l = 0;

		 while (objects)
		   {
		    /* if current fiber correspond to this object spectrum */ 
		    if (objects->IFUslitNo == theIfuSlits->ifuSlitNo &&
			objects->IFUfibNo == theIfuFibers->fibNo )
		      {
			/* cleanup the lineSigma array for this object */
			for (k=1; k<=numSkyLines; k++) 
			  {
			   lineSigma->data[k-1] = 0.0;
			  }

			/* loop on sky lines for this spectrum */
			for (k=1; k<=numSkyLines; k++)
			  {
			   /* initialize arrays for fitting sky lines */
			   for (uu=0; uu<NPOINTS; uu++)
			     {
			      anX->data[uu] = 0.0;
			      anY->data[uu] = 0.0;
			     }
			   for (uu=1; uu<=NTERMS; uu++) coeffs[uu]=0.0;

			   /*NOTE: HEADER VALUES ARE IN ANGSTROMS (24/04/02)*/
			   /* WARNING! LAMBDA of SKY LINES is in ANGSTROMS: */
			   /* convert to nm */
			   /*	   skyLam = skyLines->data[k-1]/10.;*/
			   skyLam = skyLines->data[k-1];

			   /* determine pixel position of this sky line */
			   skyPixIncF = (ABS(skyLam - XRefPixVal)) / lambdaInc;

			   /* cast to integer: truncation */
			   skyPixInc = (int)(skyPixIncF);

			   /* NOT THIS: rounding to next integer, i.e. 9.0=9 */
			   /* and 9.1 = 10 */
			   /* if ((skyPixIncF/((float)skyPixInc))>0.) 
			      skyPixInc+=1;
			   */

			   /* position of sky line in pixels */
			   if ((skyLam-XRefPixVal) >=0.) 
			     skyPix=XRefPix+skyPixInc;
			   if ((skyLam-XRefPixVal) <0.)  
			     skyPix=XRefPix-skyPixInc;

			   /* lower pixel for working arrays */
			   skyLowPix = skyPix - HALFPOINTS;

			   /* create working arrays */
			   for (j=0; j<NPOINTS; j++)
			     {
			      anX->data[j] = skyLam - (HALFPOINTS+j)*lambdaInc;

			      /*
				anY->data[j]=imageData->data[objects->rowNum] +
				skyLowPix + j;
			      */

			      anY->data[j] = imageData->data[skyLowPix + j +
						    (objects->rowNum)*specLen];
			     }

			   fit1DGauss(anX,anY,coeffs,NTERMS);

			   if (coeffs[1]<0.0 || coeffs[2]<0.0 || coeffs[3]<0.0)
			     {
			      pilMsgInfo(modName, "Slit: %3d, Fib: %3d - Negative value: coeffs[1]=%10.6f coeffs[2]=%10.6f, coeffs[3]=%10.6f\n", theIfuSlits->ifuSlitNo, theIfuFibers->fibNo, coeffs[1], coeffs[2], coeffs[3]);
			     }


			   /* WARNING!!! VERIFY THIS! we go from 1 to nterms*/
			   lineSigma->data[k-1] = coeffs[3];

			  } /* end loop on sky lines */

			/* evaluate mean sigma for this object: use median */
			/* from medianWirth */
			if (numSkyLines != 1)
			  {
			   objectPsfY->data[l]=medianWirth(lineSigma->data, 
							   numSkyLines);
			  }
			else
			  {
			   objectPsfY->data[l] = lineSigma->data[0];
			   if (lineSigma->len != 1)
			     pilMsgError(modName,"Wrong number of elements!");
			  }

			/* update IFU table for this fiber: write sigmaY */
			theIfuFibers->sigmaY = objectPsfY->data[l];

		      } /* end if this object correspond to current fiber */

		    objects = objects->next;

		    /* increment object counter */
		    l++;

		   } /* end loop on objects in this IFU image */

		 /* just to be sure ... */
		 if (l != (PSEUDOSLIT*4))
		   {
		    pilMsgError(modName,"WRONG no. of spectra in this image");
		    pilMsgInfo(modName,"SLit: %4d FIB: %4d OBJ: %4d",
			       theIfuSlits->ifuSlitNo,
			       theIfuFibers->fibNo, objects->rowNum);
		   }

		} /* end if it is not a dead fiber */
	      else if ( (theIfuFibers->fiberTrans) == -1.)
		{
		 nDeadFibs += 1;
		}

	      theIfuFibers = theIfuFibers->next;
	     } /* end loop on fibers for this pseudoslit */

	   /* start to compute PSFY group for this pseudoslit */

	   /* for the moment copy the original, to test... */
	   for (ii=0; ii<nGoodFibs; ii++) 
	     sortObjectPsfY->data[ii] = objectPsfY->data[ii];

	   /* sort in ascending order */
	   sort(nGoodFibs, sortObjectPsfY->data);

	   /* this is because kth_smallest modifies the input array */
	   for (ii=0; ii<nGoodFibs; ii++) 
	     tmpObjectPsfY->data[ii] = sortObjectPsfY->data[ii];

	   /* define the values corresponding to the selected quantiles */
	   /* first element of psfYLimit already initialized to zero, here
	    set the others */

           /* NOTE: the values in frequency->data[ii] are defined by means */ 
	   /* of nGoodFibs. The rank in the tmpObjectPsfY array follows the */
	   /* usual array indexing and goes from 0 to (nGoodFibs-1) => */
	   /* rank set to (frequency->data[ii]-1) */
	   for (ii=0; ii<numQuant; ii++)
	     psfYLimit->data[ii+1] = kthSmallest(tmpObjectPsfY->data, 
					nGoodFibs, (frequency->data[ii]-1));

	   /* loop on fibers to define the psfY group and write the */
	   /* Sky_Group parameter in the IFU Table */
	   theIfuFibers = theIfuSlits->fibers;
	   while(theIfuFibers)
	     {
	      if (theIfuFibers->fiberTrans != -1.)
		{
		 psfy = theIfuFibers->sigmaY;
		 /* check if psfY belong to a given interval BUT last one */
		 for (ii=0; ii<nIntervals-1; ii++)
		   {
		    if ( (psfy > psfYLimit->data[ii]) && 
			 (psfy <= psfYLimit->data[ii+1]) )
		      {
		       theIfuFibers->sigmaYGroup = ii+1;
		       numGrouped->data[ii] += 1;
		      }
		   
		   }
		 /* last one interval */
		 if (psfy > psfYLimit->data[nIntervals-1])
		   {
		    theIfuFibers->sigmaYGroup = nIntervals;
		    numGrouped->data[nIntervals-1] += 1;
		   }
		}
	      else if (theIfuFibers->fiberTrans == -1.)
		{
		 theIfuFibers->sigmaYGroup = -1;
		}

	      theIfuFibers = theIfuFibers->next;   
	     } /* end loop on fibers for this pseudoslit */

	   for (ii=0; ii<nIntervals-1; ii++) 
	     totalGrouped += numGrouped->data[ii];

	   if (totalGrouped != nGoodFibs) 
	     {
	      pilMsgError(modName,"WRONG grouping: n. of grouped fibs != n. "
                          "of good fibs in pseudoslit %d",
			  theIfuSlits->ifuSlitNo);
	      return EXIT_FAILURE;
	     }

	   pilMsgInfo(modName,"Done slit no. %3d ",theIfuSlits->ifuSlitNo);

	   if ( (totalGrouped + nDeadFibs) != PSEUDOSLIT )
	     {
	      pilMsgError(modName,
			  "ERROR: %d good fibs and %d dead fibs for slit %d",
			  totalGrouped,nDeadFibs,theIfuSlits->ifuSlitNo);
	      return EXIT_FAILURE;
	     }

	   /* cleanup previous allocations for this slit*/
	   deleteFloatArray(objectPsfY);
	   deleteFloatArray(sortObjectPsfY);
	   deleteFloatArray(tmpObjectPsfY);

	   /* increment counter of total n. of good fibers on this image */
	   nTotGood += totalGrouped;
	   nTotBad += nDeadFibs;

	   theIfuSlits = theIfuSlits->next;
	  } /* end loop on pseudoslits for this quadrant */

       } /* end if this is the right quadrant */

     ifuQuads = ifuQuads->next;
    } /* end loop on quadrants in the IFU table */

  if ((nTotGood+nTotBad) != (4*PSEUDOSLIT))
    {
     pilMsgError(modName,"Wrong number of fibers considered in this image");
     return EXIT_FAILURE;
    }

  deleteFloatArray(anX);
  deleteFloatArray(anY);
  freeFloatVector(coeffs,1,NTERMS);

  deleteFloatArray(skyLines);
  deleteFloatArray(lineSigma);
  deleteFloatArray(psfYLimit);

  deleteObjectObject(objects);
  deleteIntArray(frequency);
  deleteIntArray(numGrouped);

  deleteIfuQuad(ifuQuads);
  deleteIfuSlit(theIfuSlits);
  deleteIfuFiber(theIfuFibers);

  return EXIT_SUCCESS;
}


VimosImage *
VmIfuCrossTalk(VimosImage *imageData, VimosIfuTable *ifuTable, 
               VimosExtractionTable *extractionTable, 
               float wLenStart, float wLenEnd, int xLen, int yLen, 
               int quadNum)
{
  int numFibers, numRows, fibNum;
  int firstX, whereZero, lowerPix, upperPix, countCosmic, left, right;
  int jj, ii, i, j, k, l, ll, m, n, peakCentX, start;

  VimosUlong32 peak, pminY, pmaxY;

  float sigma, z, noise, diff, tmpGauss;

  double x, y, maxY, minY;

  char modName[] = "VmIfuCrossTalk";

  VimosIntArray *startXPeak, *startYPeak, *endXPeak, *endYPeak, *numRowsArray;

  VimosFloatArray *moduleProfile, *medModuleProfile, *medModuleProfile2;
  VimosFloatArray *crossCut, *normCrossCut, *crossCut2, *normCrossCut2;
  VimosFloatArray *tmpNormCrossCut, *tmpNormCrossCut2, *xFWHMArray, *startX;
  VimosFloatArray *takeThis, *medianPeaks2, *medianPeaks, *peakX;

  VimosImage *peakImage, *startXImage, *crossTalkMap;

  VimosExtractionSlit *extSlits;
  VimosIfuQuad *ifuQuads;
  VimosIfuSlit *theIfuSlits;
  VimosIfuFiber *theIfuFibers;

  VimosImage *outData;

  pilMsgInfo (modName, "Start computing cross-talk correction...");

  /* create output image and initialize to zero, then overwrite on it
     the cross-dispersion cuts corrected for crosstalk */
  outData = newImageAndAlloc(xLen, yLen);
  for (j=0; j<yLen; j++)
     {
      for (k=0; k<xLen; k++)
         {
           outData->data[k + j*xLen] = 0.0;
         }
     }

  /* define two VimosImage (Y size equal to the input, as many columns as the 
     fibers in this pseudoslit) images containing one the position of the first
     X pixel and the other the position of the peak of each fiber along Y */
  startXImage = newImageAndAlloc(PSEUDOSLIT, yLen);
  peakImage = newImageAndAlloc(PSEUDOSLIT, yLen);

  /* define the crosstalk map image. For each cross dispersion cut a, a 
     crosstalk map is computed. Each row of the crosstalk map contains the 
     crosstalk contribution of one fiber of the cross-dispersion cut */
  crossTalkMap = newImageAndAlloc(xLen, PSEUDOSLIT);

  /* define the array for the cross dispersion cut and for cosmic-cleaned 
     cross dispersion cut*/
  crossCut = newFloatArray(imageData->xlen);
  crossCut2 = newFloatArray(xLen);

  /* define normalized crossCut to compute medians, and crossCut2 to 
     re-compute medians  */
  normCrossCut = newFloatArray(xLen);
  normCrossCut2 = newFloatArray(xLen);

  /* define temporary files for medianWirth: copies of normCrossCut and 
     normCrossCut2 */
  tmpNormCrossCut = newFloatArray(xLen);
  tmpNormCrossCut2 = newFloatArray(xLen);

  /* working array to store rows of peakXImage (during cross-cut loop) */
  peakX = newFloatArray(PSEUDOSLIT);

  /* define array containing the peak fluxes obtained by the median process */
  medianPeaks = newFloatArray(PSEUDOSLIT);

  /* define array containing the peak fluxes obtained by the median process:
     this is for when re-computing medians */
  medianPeaks2 = newFloatArray(PSEUDOSLIT);

  /* define theoretical gaussian profile for the pseudoslit */
  moduleProfile = newFloatArray(xLen);

  /* define a second theoretical gaussian profile for the pseudoslit, will have
     each gaussian scaled to the fiber peak flux */
  medModuleProfile = newFloatArray(xLen);

  /* define the theoretical gaussian profile for the pseudoslit, with each 
     gaussian scaled to the fiber peak flux: determined by second pass of 
     median process */
  medModuleProfile2 = newFloatArray(xLen);

  /* starting X position of each fiber along a cross dispersion cut */
  startX = newFloatArray(PSEUDOSLIT);

  /* work array to be used for crosstalk correction */
  takeThis = newFloatArray(xLen);

  /* define "left" and "right" as the number of adjacent fibers to consider on 
     the two sides of the fiber which is undergoing crosstalk correction */
  left =  HOWMANYFIBS;
  right = HOWMANYFIBS;

  /* for each IFU pseudoslit, check on the image the maximum/minimum Y pixel
     for all fibers, to determine which is the Y coordinate where to start with
     the cross dispersion cuts. This is done using wLenStart, wLenEnd, and the
     inverse dispersion solution, ccdX[0] and ifuPeakX in extraction table */

  ifuQuads = ifuTable->quads;

  while(ifuQuads)
    {
     /* take the right quadrant in the IFU table */
     if (ifuQuads->quadNo == quadNum)
       {
        /* if this is the right quadrant, take pseudoslits */
        theIfuSlits = ifuQuads->ifuSlits;
        while(theIfuSlits)
          {
           theIfuFibers = theIfuSlits->fibers;
	   /* initialize fiber counter */
	   numFibers = 0;

	   /* initialize arrays for x,y lower and upper limits of peak
	      position for each fiber (i.e. at wLenStart and wLenEnd) */
	   startXPeak = newIntArray(PSEUDOSLIT);
	   startYPeak = newIntArray(PSEUDOSLIT);
	   endXPeak = newIntArray(PSEUDOSLIT);
	   endYPeak = newIntArray(PSEUDOSLIT);

	   /* initialize the arrays with fwhm and numRows for all fibers in 
	      the pseudoslit */
           xFWHMArray = newFloatArray(PSEUDOSLIT);
	   numRowsArray = newIntArray(PSEUDOSLIT);

	   /* initialize min/max in Y for all the fibers in that pseudoslit */
	   maxY = 0.;
	   minY = (double)yLen;

	   /* initialize startXImage and peakImage to -1.0 for the pseudoslit*/
	   for (jj=0; jj<yLen; jj++)
	      {
	       for (ll=0; ll<PSEUDOSLIT; ll++)
		  {
		   startXImage->data[ll +jj*PSEUDOSLIT] = -1.0;
		   peakImage->data[ll +jj*PSEUDOSLIT] = -1.0;
		  }
	      }

           while(theIfuFibers)
             {
              /* set the value of FWHM for the fibers */
              xFWHMArray->data[numFibers] = theIfuFibers->fiberPwidth;

              /* take slits from the extraction table */
              extSlits = extractionTable->slits;

              while(extSlits)
                {
                 /* if the extraction slit correspond to the IFU fiber... */
		 if ( (extSlits->IFUslitNo == theIfuSlits->ifuSlitNo) &&
                      (extSlits->IFUfibNo == theIfuFibers->fibNo) )
                   {
		    /* take the fiber sequential number in the pseudoslit:
		       needed to write startXImage */
		    /* theIfuFibers->fibNo goes from 1 to 400, fibNum goes from
		       0 to 399 */
		    fibNum = theIfuFibers->fibNo - 1;

		    /* set the value of numRows for the fiber */
		    numRows = extSlits->numRows;
		    numRowsArray->data[fibNum] = numRows;

		    /* we need the value of the fiber peak X pixel WITH
		       RESPECT to firstX, to get the right invDis and crvPol */

		    /* this is the first pixel in the fiber */
		    firstX = (int)(extSlits->ccdX->data[0]);

		    /* this is the fiber peak X pixel WITH RESPECT to firstX
		       cast to integer (truncation) */
		    /* NOTE: this casting is needed to get the right pixel
		       in the slit (in the range from 0 to numRows-1) */

		    /* CHECK THIS !!!!! */
		    peakCentX = (int)(extSlits->IFUfibPeakX) - firstX;

		    /* PER PROVARE A COMPILARE: */
		    /*peakCentX = 0.;*/

		    y = 0.;
		    x = 0.;

		    /* compute peak Y position at wLenStart */
		    y = extSlits->ccdY->data[peakCentX] + 
		        computeDistModel1D(extSlits->invDis[peakCentX], 
					   wLenStart);

		    /* check minimum */
                    if (y < minY) minY = y;

		    /* use curvature polinomial to get peak x */
		    x = extSlits->ccdX->data[peakCentX] + 
		        computeDistModel1D(extSlits->crvPol[peakCentX], y);

		    /* cast to integer (truncation) */
		    startXPeak->data[fibNum] = x;
		    startYPeak->data[fibNum] = y;

		    y = 0.;
		    x = 0.;

		    /* compute peak Y position at wLenEnd */
		    y = extSlits->ccdY->data[peakCentX] + 
		        computeDistModel1D(extSlits->invDis[peakCentX], 
					   wLenEnd);

		    /* check maximum */
		    if (y > maxY) maxY = y;
		    /* use curvature polinomial to get peak x */
		    x = extSlits->ccdX->data[peakCentX] + 
		        computeDistModel1D(extSlits->crvPol[peakCentX],y);

		    /* cast to integer (truncation) */
		    endXPeak->data[fibNum] = x;
		    endYPeak->data[fibNum] = y;

		    /* check if endYPeak > yLen: if so, set it to last pixel 
		       (yLen-1). For endXPeak: being the  fiber peak , it's 
		       always less than xLen ? */
		    if (endYPeak->data[fibNum] >= yLen)
		      endYPeak->data[fibNum] = (yLen - 1);

		    /* increment fiber counter */
		    numFibers += 1;

                   } /* end if the extraction slit = the IFU fiber */

		 extSlits = extSlits->next;

                } /* end loop on extraction slits */

	      /* ...just to check ...*/
	      if ( (maxY > (double)yLen) || (minY < 0.) ) 
		pilMsgError(modName,
			    "%3d: ERROR! Wrong max. Y pixel for pseudoslit",
			    theIfuSlits->ifuSlitNo);

              theIfuFibers = theIfuFibers->next;

             } /* end loop on fibers for this IFU pseudoslit */

	   if ((numFibers) != PSEUDOSLIT)
	     pilMsgError(modName,
			 "ERROR! total of: %5d fibers for the pseudoslit: %3d",
			 numFibers, theIfuSlits->ifuSlitNo);

           /* take the integer min/max in Y (truncation) */
           pminY = (VimosUlong32) minY;
           pmaxY = (VimosUlong32) maxY;

	   /* if pmaxY > image size, set it to image size */
	   if (pmaxY >= (VimosUlong32)yLen) pmaxY = (yLen - 1);


	   /* now restart loop on fibers and extraction slits to compute,
	      between pminY and pmaxY for each fiber: 1) peakCentX along 
	      dispersion and write it in peakImage 2) firstX and write it in 
	      startXImage (needed for evaluating median fluxes) */

	   /* initialize fiber counter */
	   numFibers = 0;

           theIfuFibers = theIfuSlits->fibers;
           while(theIfuFibers)
             {
              /* take slits from the extraction table */
              extSlits = extractionTable->slits;

              while(extSlits)
                {
                 /* if the extraction slit correspond to the IFU fiber... */
		 if ( (extSlits->IFUslitNo == theIfuSlits->ifuSlitNo) &&
                      (extSlits->IFUfibNo == theIfuFibers->fibNo) )
		   {
		    /* this is the first pixel in the fiber : NOT NEEDED HERE*/
		    firstX = (int)(extSlits->ccdX->data[0]);

		    fibNum = theIfuFibers->fibNo - 1;

		    /* CHECK THE LIMITS FOR THIS LOOP */
		    for (jj=pminY; (VimosUlong32)jj<=pmaxY; jj++)
		       {
			peakImage->data[fibNum + jj*PSEUDOSLIT] = 
			    extSlits->ccdX->data[peakCentX] + 
			    computeDistModel1D(extSlits->crvPol[peakCentX],jj);
			
			startXImage->data[fibNum + jj*PSEUDOSLIT] = 
			    extSlits->ccdX->data[0] + 
			    computeDistModel1D(extSlits->crvPol[0],jj);

			/* check that, for some fibers, "tracing" betweeen 
			   min-max does not give pixels out of the CCD size */
			if (peakImage->data[fibNum + jj*PSEUDOSLIT] < 0)
			  peakImage->data[fibNum + jj*PSEUDOSLIT] = 0;

			if (peakImage->data[fibNum + jj*PSEUDOSLIT] 
			    >= xLen)  peakImage->data[fibNum + 
						 jj*PSEUDOSLIT] = xLen-1;

			if (startXImage->data[fibNum + jj*PSEUDOSLIT] < 0)
			  startXImage->data[fibNum + jj*PSEUDOSLIT]  = 0;

			if (startXImage->data[fibNum + jj*PSEUDOSLIT] 
			    >= xLen) startXImage->data[fibNum + 
						 jj*PSEUDOSLIT] = xLen-1;
			
		       }

		    /* increment fiber counter */
		    numFibers += 1;
		   }

		 extSlits = extSlits->next;

                } /* end loop on extraction slits */

              theIfuFibers = theIfuFibers->next;

             } /* end loop on fibers for this IFU pseudoslit */
	   if ((numFibers) != PSEUDOSLIT)
	     pilMsgError(modName,
			 "ERROR! total no. of: %5d fibers for pseudoslit: %3d",
			 numFibers, theIfuSlits->ifuSlitNo);

	   /* some checks, maybe redundant */
	   for (m=0; m<PSEUDOSLIT; m++)
	     {
	      for (jj=pminY; (VimosUlong32)jj<=pmaxY; jj++)
	         {
		  if (peakImage->data[m + jj*PSEUDOSLIT] == -1.0) 
		       puts("ERROR in peakImage!");
		  if (startXImage->data[m + jj*PSEUDOSLIT] == -1.0) 
		       puts("ERROR in startXImage!");
		 }

	     } /* end check on limits for each fiber */

	   /**********************************************************/
	   /* now we can start with the actual cross-dispersion cuts */
	   /* for this pseudoslit                                    */
	   /**********************************************************/

	   /* TO BE CHECKED: we have (pmaxY-pminY+1) cross disp. cuts to be
	      done for this pseudoslit, so set loop from =pminY to <=pmaxY */

	   /* counter for rows in peakImage */
	   i = 0;
	   for (j=pminY; (VimosUlong32)j<=pmaxY; j++)
	      {

	       /* some initializations, just to be safe..*/

	       /* CHECK THAT THESE NUMBERS DO NOT PRODUCE MISTAKES */
	       for (k=0; k<xLen; k++)
		  {
		   crossCut->data[k] = 0.0;
		   crossCut2->data[k] = 0.0;
		   normCrossCut->data[k] = 0.0;
		   normCrossCut2->data[k] = 0.0;
		   moduleProfile->data[k] = 0.0;
		   medModuleProfile->data[k] = 0.0;
		   medModuleProfile2->data[k] = 0.0;
		  }
	       
	       for (l=0; l<PSEUDOSLIT; l++)
		  {
		   startX->data[l] = 0.0;
		   peakX->data[l] = 0.0;
		   medianPeaks->data[l] = 0.0;
		   medianPeaks2->data[l] = 0.0;
		  }

	       /* now read crossCut values from input image */
               for (k=0; k<xLen; k++)
		 crossCut->data[k] = imageData->data[k + j*xLen];

	       /* get fiber peaks and first Xpixels for cross-dispersion cut */
	       for (l=0; l<PSEUDOSLIT; l++)
		 {
		  startX->data[l] = startXImage->data[l + j*PSEUDOSLIT];
		  peakX->data[l] = peakImage->data[l + j*PSEUDOSLIT];
		 }

	       /* create module profile */
	       for (l=0; l<PSEUDOSLIT; l++)
		  {
		   /* take peak, cast to integer (truncation) */
		   peak = (VimosUlong32) (peakX->data[l]);

		   /* this is the relation between FWHM and sigma */
		   sigma = xFWHMArray->data[l] / FWHM_TO_SIGMA ;

		   /* for each fiber, compute the gaussian in the range 
		      (-10sigma, +10sigma) around the peak X pixel */
		   lowerPix = peak - ( (VimosUlong32)( sigma*NSIGMA) );
		   if (lowerPix < 0) lowerPix = 0;

		   upperPix = peak + ( (VimosUlong32)( sigma*NSIGMA) );
		   if (upperPix >= xLen) upperPix = xLen -1;

		   for (k=lowerPix; k<upperPix; k++)
		      {
		       z = (pow (((float)k - (float)peak), 2.) ) 
			 / (2. * (pow(sigma,2.)));
		       tmpGauss = exp(-z);
		       moduleProfile->data[k] += tmpGauss;
		      }
		  }

	       /* normalize crossCut array to the theoretical module profile */
	       for (k=0; k<xLen; k++)
		 normCrossCut->data[k] = crossCut->data[k] / 
		                         moduleProfile->data[k];

	       /* for the i-th fiber on the normalized crosscut, take its
		  numRows[i] pixels and do the median to estimate peak flux */
	       for (l=0; l<PSEUDOSLIT; l++)
		  {
		   /* cast to integer (truncation) */
		   start = (VimosUlong32) (startX->data[l]);

		   /* in principle, this re-initialization for each fiber is 
		      not needed but, since medianWirth modifies the input 
		      array, do it just to be safe*/
		   for(k=0; k<xLen; k++) tmpNormCrossCut->data[k]=
					   normCrossCut->data[k];

		   /* for each fiber, apply medianWirth to numRows pixels in 
		      the tmpcrosscut array, starting from the first pixel 
		      belonging to the fiber */
		   /* CHECK THIS!!! */
		   medianPeaks->data[l] = medianWirth(
					      &tmpNormCrossCut->data[start],
					      numRowsArray->data[l]);
		  }

	       /* now recompute the theoretical module profile, this time 
		  scaling each gaussian to its peak flux */
	       for (l=0; l<PSEUDOSLIT; l++)
		  {
		   /* take peak, cast to integer (truncation) */
		   peak = (VimosUlong32) (peakX->data[l]);

		   /* this is the relation between FWHM and sigma */
		   sigma = xFWHMArray->data[l] / FWHM_TO_SIGMA ;

		   /* for each fiber, compute the gaussian in the range 
		      (-10sigma, +10sigma) around the peak X pixel */

		   lowerPix = peak - ( (VimosUlong32)( sigma*NSIGMA) );
		   if (lowerPix < 0) lowerPix = 0;

		   upperPix = peak + ( (VimosUlong32)( sigma*NSIGMA) );
		   if (upperPix >= xLen) upperPix = xLen -1;
		   
		   for (k=lowerPix; k<upperPix; k++)
		      {
		       z = ( pow(((float)k - (float)peak), 2.) ) 
			 / (2. * (pow(sigma,2.)));
		       tmpGauss = medianPeaks->data[l] * exp(-z);
		       medModuleProfile->data[k] += tmpGauss;
		      }
		  }

	       /* compare medModuleProfile with crossCut to find out if
		  there are any cosmics... In case, substitute pixel value
	          in crossCut with the pixel value in medModuleProfile */

	       countCosmic = 0;

	       for (k=0; k<xLen; k++)
		  {
		   noise = sqrt(medModuleProfile->data[k]);
		   diff = abs(crossCut->data[k] - medModuleProfile->data[k]);
		   if (diff > (noise * LIMIT))
		     {
		      crossCut2->data[k] = medModuleProfile->data[k];
		      countCosmic ++; 
		     }
		  }

	       /******************************************************/
	       /* done cosmic cleaning */
	       /******************************************************/

	       /* now on the "new" crossCut array make a second step to
		  evaluate the medians and recompute medModuleProfile */

	       /* first normalize crossCut2 to theoretical module profile */
	       for (k=0; k<xLen; k++)
		 normCrossCut2->data[k] = crossCut2->data[k] / 
		                          moduleProfile->data[k];

	       /* recompute medians for each fiber */
	       for (l=0; l<PSEUDOSLIT; l++)
		  {
		   /* cast to integer (truncation) */
		   start = (VimosUlong32) (startX->data[l]);

		   /* in principle, this re-initialization for each fiber is 
		      not needed but, since medianWirth modifies the input 
		      array, do it just to be safe*/
		   for(k=0; k<xLen; k++) tmpNormCrossCut2->data[k]=
					   normCrossCut2->data[k];

		   /* for each fiber, apply medianWirth to numRows pixels in 
		      the tmpcrosscut array, starting from the first pixel 
		      belonging to the fiber */
		   /* CHECK THIS!!! */
		   medianPeaks2->data[l] = medianWirth(
					      &tmpNormCrossCut2->data[start],
					      numRowsArray->data[l]);

		  } /* end loop to recompute medians of peak fluxes */

	       /* last, recompute the theoretical module profile, scaling
		  each gaussian to the new peak flux */
	       /* THIS PROFILE WILL NOT BE USED IN THE COMPUTATIONS BELOW,
		  BUT IT CAN BE USEFUL DURING TEST PHASE (WILL DISAPPEAR
		  LATER?) */

	       for (l=0; l<PSEUDOSLIT; l++)
		  {
		   /* take peak, cast to integer (truncation) */
		   peak = (VimosUlong32) (peakX->data[l]);

		   /* this is the relation between FWHM and sigma */
		   sigma = xFWHMArray->data[l] / FWHM_TO_SIGMA ;

		   /* for each fiber, compute the gaussian in the range 
		      (-10sigma, +10sigma) around the peak X pixel */

		   lowerPix = peak - ( (VimosUlong32)( sigma*NSIGMA) );
		   if (lowerPix < 0) lowerPix = 0;

		   upperPix = peak + ( (VimosUlong32)( sigma*NSIGMA) );
		   if (upperPix >= xLen) upperPix = xLen -1;

		   for (k=lowerPix; k<upperPix; k++)
		      {
		       z = ( pow(((float)k - (float)peak), 2.) ) 
			 / (2. * (pow(sigma,2.)));
		       tmpGauss = medianPeaks2->data[l] * exp(-z);
		       medModuleProfile2->data[k] += tmpGauss;
		      }

		  } /* end loop to compute final theoretical profile */


	       /* initialize the crosstalk map */
	       for (l=0; l<PSEUDOSLIT; l++)
		  {
		   for (k=0; k<xLen; k++)
		     crossTalkMap->data[k + l*xLen] = 0.0;
		  }

	       /* start writing the crosstalk map */
	       for (l=0; l<PSEUDOSLIT; l++)
		  {
		   /* take peak, cast to integer (truncation) */
		   /* WARNING!!!!! CHECK IF PEAKX IS OK!!!! */
		   peak = (VimosUlong32) (peakX->data[l]);

		   /* this is the relation between FWHM and sigma */
		   sigma = xFWHMArray->data[l] / FWHM_TO_SIGMA ;

		   /* for each fiber, compute the gaussian in the range 
		      (-10sigma, +10sigma) around the peak X pixel */

		   lowerPix = peak - ( (VimosUlong32)( sigma*NSIGMA) );
		   if (lowerPix < 0) lowerPix = 0;

		   upperPix = peak + ( (VimosUlong32)( sigma*NSIGMA) );
		   if (upperPix >= xLen) upperPix = xLen -1;
		   
		   for (k=lowerPix; k<upperPix; k++)
		      {
		       z = ( pow(((float)k - (float)peak), 2.) ) 
			 / (2. * (pow(sigma,2.)));

		       crossTalkMap->data[k+l*xLen] = 
		         medianPeaks2->data[l] * exp(-z);
		     
		      }
		  } /* end loop to write the crosstalk map for this cut */


	       /************************/
	       /* crosstalk correction */
	       /************************/

	       /* take each fiber in the cross-dispersion cut */
	       for (m=0; m<PSEUDOSLIT; m++)
		  {
		   /* this is to count how many times it happens that the 
		      crosstalk subtraction gives a value <0 of the data */
		   whereZero = 0;

		   /* "left" and "right" are the number of adjacent fibers to 
		      consider on the two sides of the fiber which is
		      undergoing crosstalk correction */

		   /* WARNING!!! CHECK THIS PARAMETRIZATION! */
		   /* limit for first fibers... */
		   if (m < HOWMANYFIBS) left = m;
		   /* limit for last fibers... */
		   if ( (m + HOWMANYFIBS) >= PSEUDOSLIT ) 
		        right = (PSEUDOSLIT - m -1);

		   /* take the pixels of the fiber in the crosscut */
		   for (n=(m-left); n<=(m+right) && n!=m; n++)
		      {
		       /* initialize working array for the crosstalk
			  contribution of this fiber */
		       for (k=0; k<xLen; k++) takeThis->data[k] = 0.0;

		       for (k=0; k<xLen; k++)
			 takeThis->data[k] = crossTalkMap->data[k + n*xLen];

		       /* subtract crosstalk given by takeThis from the cross
			  dispersion pixels corresponding to the m-th fiber */

		       /* take numRows for the fiber to be corrected */
		       numRows = numRowsArray->data[m];

		       /* WARNING!!! CHECK startX!!!! */
		       /* cast to integer (truncation) */
		       /* WARNING!!! CHECK THE m INDEX (PREVIOUSLY I WROTE l)*/
		       start = (VimosUlong32) (startX->data[m]);

		       /* actual correction */
		       for (ii=start; ii<=(start+numRows -1); ii++)
		          {
			   crossCut2->data[ii] -= takeThis->data[ii];

			   /* if less than zero, set to zero. This can happen
			      due to possonian noise */
			   if (crossCut2->data[ii] < 0.0) 
			     {
			      crossCut2->data[ii] = 0.0;
			      whereZero += 1;
			     }
		          }

		      } /* end loop on +/- 2 fibers on each side */

		  } /* done crosstalk corr. for this cross-dispersion cut */

	       /* now, overwrite the "new" (corrected) cross-dispersion cut
	          in the output image */

	       for (ii=0; ii<xLen; ii++)
	          {
	           outData->data[ii + j*xLen] = crossCut2->data[ii];
	          }

	       /* counter for rows in peakXImage: incremented */
	       i++;

	      } /* end loop on cross-dispersion cuts for this IFU pseudoslit */

	   /* just to check */
	   if ((VimosUlong32)i != (pmaxY - pminY +1))
	     {
	      pilMsgError(modName,
		  "ERROR! %5d cross-dispersion cuts done for pseudoslit: %5ld",
			  i, theIfuSlits->ifuSlitNo);
	      pilMsgError(modName,
		      "       while the evaluated (pmaxY - pminY +1) is: %5ld",
			  (pmaxY - pminY +1));
	      return NULL;
	     }
	   else if  ((VimosUlong32)i == (pmaxY - pminY +1))
	     {
	      pilMsgInfo(modName,
			 "%5ld cross dispersion cuts done for pseudoslit %3d", 
			 i, theIfuSlits->ifuSlitNo);
	     }
	   else
	     {
	      pilMsgError(modName,"Unable to compute no. of cross-disp. cuts");
	      return NULL;
	     }
	   deleteIntArray(startXPeak);
	   deleteIntArray(startYPeak);
	   deleteIntArray(endXPeak);
	   deleteIntArray(endYPeak);

           deleteFloatArray(xFWHMArray);
	   deleteIntArray(numRowsArray);

           theIfuSlits = theIfuSlits->next;
          } /* end loop on IFU pseudoslits */

       } /* end "if the quadrant is the right one" */

     ifuQuads = ifuQuads->next;
    } /* end loop on IFU quadrants */

  pilMsgInfo(modName, "Done crosstalk correction");



  /* cleaning up */
  deleteFloatArray(moduleProfile);
  deleteFloatArray(medModuleProfile);
  deleteFloatArray(medModuleProfile2);
  deleteFloatArray(takeThis);
  deleteFloatArray(crossCut);
  deleteFloatArray(normCrossCut);
  deleteFloatArray(crossCut2);
  deleteFloatArray(normCrossCut2);
  deleteFloatArray(startX);
  deleteFloatArray(medianPeaks);
  deleteFloatArray(medianPeaks2);

  deleteImage(peakImage);
  deleteImage(startXImage);
  deleteImage(crossTalkMap);

  return outData;

}


int
VmIfuGetTransmission(VimosImageSet *imageSet, 
                     VimosIfuTable *theMainIfuTable, char *toBeFitted)
{
  int i,j,k;

  int refL, refM, refFibNo, refSlitNo, refQuadNo, quad;
  float refTrans, lineFlux, refLineFlux;

  VimosFloatArray *anX, *anY;
  float *coeffs, skyLam;

  float XRefPixVal, YRefPixVal, XRefPix, YRefPix, lambdaInc, skyPixIncF;
  int numSkyLines, skyPixInc;
  int Fit=0;
  int quadNo, ifuSlit, ifuFib;
  int skyPix, skyLowPix;
  VimosFloatArray *skyLines;

  VimosIfuQuad *ifuQuads;
  VimosIfuSlit *theIfuSlits;
  VimosIfuFiber *theIfuFibers;
  VimosSingleImage *ifuIma;
  VimosObjectObject *ifuImaObj;
  VimosObjectTable *ifuImaObjTab;

  puts("WARNING! DEFINE THE REFERENCE PIXEL!");
  puts("WARNING! LAMBDA of SKY LINES in ANGSTROMS!");

  
  /* read L,M for the reference fiber */
  /* IS IT OK TO READ THEM FROM IFU TABLE DESCRIPTORS? */
  readIntDescriptor(theMainIfuTable->descs, "ESO PRO REF L", &refL, "");
  readIntDescriptor(theMainIfuTable->descs, "ESO PRO REF M", &refM, "");

  /* read sky lines from objectTable of first image */
  if (toBeFitted == "ONE") Fit = 0;
  if (toBeFitted == "ALL") Fit = 1;

  switch(Fit)
    {
    default:
    case 0 :
      {
       numSkyLines = 1;
       skyLines = newFloatArray(1);
       skyLines->data[0] = THE_SKY_LINE;
       break;
      }
    case 1 :
      {
       readSkyLines(imageSet->images->objectTable->descs, &numSkyLines, 
                     &skyLines);
       break;
      }
    }

  /* now get transmission for reference fiber: take the Standard IFU Table */
  /* and look for ref fiber transmission, quadrant, slit and fiber number */

  ifuQuads = theMainIfuTable->quads;
  i = 0;
  refTrans = 0.;
  refFibNo = 0;
  refSlitNo = 0;
  refQuadNo = 0;

  while(ifuQuads)
    {
      theIfuSlits = ifuQuads->ifuSlits;
      while(theIfuSlits)
        {
         theIfuFibers = theIfuSlits->fibers;
         while(theIfuFibers)
           {
            if ( (theIfuFibers->fiberL == refL) && 
                 (theIfuFibers->fiberM == refM) )
              {
               refTrans = theIfuFibers->fiberTrans;
	       refFibNo = theIfuFibers->fibNo;
	       refSlitNo = theIfuSlits->ifuSlitNo;
	       refQuadNo = ifuQuads->quadNo;
               i++;
              }
            theIfuFibers = theIfuFibers->next;
           }

        theIfuSlits = theIfuSlits->next;

       }
     ifuQuads = ifuQuads->next;
    }

  /* check that we find just one reference fiber */
  if (i != 1)
    {
     puts("error in selection of reference fiber");
     return EXIT_FAILURE;
    }

  /* initialize arrays */
  anX=newFloatArray(NPOINTS);
  anY=newFloatArray(NPOINTS);
  coeffs=floatVector(1,NTERMS);

  /* now first call to evalLineFlux for the reference fiber and store */
  /* the reference Line Flux */
  ifuIma = imageSet->images;
  while (ifuIma)
    {
     /* read the quadrant number for this image */
     readIntDescriptor(ifuIma->theImage->descs, "ESO QUAD", &quad, "");

     if (quad == refQuadNo)
       {
	ifuImaObjTab = ifuIma->objectTable;

        /* read descriptors from each Object Table */
	/* FROM DRS DOCUMENT:
	   Value of reference pixel in X (packed-spectra format): nm */
        readFloatDescriptor(ifuImaObjTab->descs, "CRVAL1", &XRefPixVal, "");
        readFloatDescriptor(ifuImaObjTab->descs, "CRVAL2", &YRefPixVal, "");
        readFloatDescriptor(ifuImaObjTab->descs, "CRPIX1", &XRefPix, "");
        readFloatDescriptor(ifuImaObjTab->descs, "CRPIX2", &YRefPix, "");
	/*  FROM DRS DOCUMENT:
	    WARNING: Wavelength step (packed-spectra format): Angstroms */
        readFloatDescriptor(ifuImaObjTab->descs, "CDELT1", &lambdaInc, "");

	/* WARNING! TRANSFORM lambdaInc from angstroms to nm for consistency */
	lambdaInc = lambdaInc/10.;

        ifuImaObj = ifuImaObjTab->objs;
        while (ifuImaObj)
          {
           if ( (ifuImaObj->IFUslitNo == refSlitNo) &&
                (ifuImaObj->IFUfibNo == refFibNo) )
             {
	      refLineFlux = 0.;

	      /* loop on sky lines */
       	      for (k=1; k<=numSkyLines; k++)
       	         {
       	          /* WARNING! EXPECTING LAMBDA of SKY LINES in ANGSTROMS */
                  skyLam = skyLines->data[k-1]/10.;

	          /* determine pixel position of this sky line */
                  skyPixIncF = (ABS(skyLam - XRefPixVal))/
                                   lambdaInc;

		  /* cast to integer: truncation */
                  skyPixInc = (int)(skyPixIncF);

		  /* rounding to next integer, i.e. 9.0 = 9 and 9.1 = 10 */
		  /* NO MORE USED:*/
		  /* if ( (skyPixIncF/((float)skyPixInc))>0.)  skyPixInc+=1; */

		  /* position of sky line in pixels */
                  if ((skyLam-XRefPixVal) >= 0.) skyPix = XRefPix + skyPixInc;
                  if ((skyLam-XRefPixVal) < 0.)  skyPix = XRefPix - skyPixInc;

		  /* lower pixel for working arrays */
                  skyLowPix = skyPix - (HALFPOINTS);

                  /* create working arrays */
                  for (j=0; j<NPOINTS; j++)
		    {
		     anX->data[j] = skyLam - (HALFPOINTS+j)*lambdaInc;
                     anY->data[j] = ifuIma->theImage->data[ifuImaObj->rowNum]+
		                    skyLowPix + j;
                    }

                  refLineFlux += evalLineFlux(anX,anY,coeffs,NTERMS);
                 }
	      /* ended loop on sky lines, evaluate mean flux */
              refLineFlux /= (float)numSkyLines;

             }

	    ifuImaObj = ifuImaObj->next;

          }
      }

    ifuIma = ifuIma->next;

   }

  /* now start loop to compute relative fiber trans for all 6400 fibers */
  ifuIma = imageSet->images;

  while (ifuIma)
    {
     ifuImaObjTab = ifuIma->objectTable;

     /* read quadrant number */
     readIntDescriptor(ifuImaObjTab->descs, "ESO PRO QUAD", &quadNo, "");

     ifuImaObj = ifuImaObjTab->objs;     

     while (ifuImaObj)
       {
        ifuSlit = ifuImaObj->IFUslitNo;
        ifuFib = ifuImaObj->IFUfibNo;

	/* loop on sky lines */
       	for (k=1; k<=numSkyLines; k++)
       	   {
            /* cleanup previous allocations, if they exist */
            deleteFloatArray(anX);
            deleteFloatArray(anY);
            freeFloatVector(coeffs,1,NTERMS);

            /* initialize arrays */
            anX=newFloatArray(NPOINTS);
            anY=newFloatArray(NPOINTS);
            coeffs=floatVector(1,NTERMS);

       	    /* WARNING! EXPECTING LAMBDA of SKY LINES in ANGSTROMS */
            skyLam = skyLines->data[k-1]/10.;

	    /* determine pixel position of this sky line */
            skyPixIncF = (ABS(skyLam - XRefPixVal))/
                                   lambdaInc;
            skyPixInc = (int)(skyPixIncF);

	    /* rounding to next integer, i.e. 9.0 = 9 and 9.1 = 10 */
            if ( (skyPixIncF/((float)skyPixInc)) > 0.)  skyPixInc += 1;

	    /* position of sky line in pixels */
            if ((skyLam-XRefPixVal) >= 0.) skyPix = XRefPix + skyPixInc;
            if ((skyLam-XRefPixVal) < 0.)  skyPix = XRefPix - skyPixInc;


	    /* lower pixel for working arrays */
            skyLowPix = skyPix - HALFPOINTS;

            /* create working arrays */
            for (j=0; j<NPOINTS; j++)
	       {
		anX->data[j] = skyLam - (HALFPOINTS + j)*lambdaInc;
                anY->data[j] = ifuIma->theImage->data[ifuImaObj->rowNum] +
		               skyLowPix + j;
               }

            lineFlux += evalLineFlux(anX,anY,coeffs,NTERMS);
           }
	/* ended loop on sky lines, evaluate mean flux */
        lineFlux /= (float)numSkyLines;

        /* loop on the Standard IFU Table to update fiber transmission */
        ifuQuads = theMainIfuTable->quads;
        while (ifuQuads)
          {
	   if (ifuQuads->quadNo == quadNo)
             {
              theIfuSlits = ifuQuads->ifuSlits;
              while(theIfuSlits)
                {
       	         if (theIfuSlits->ifuSlitNo == ifuSlit)
                   {
                    theIfuFibers = theIfuSlits->fibers;
                    while(theIfuFibers)
                      {
       		       if (theIfuFibers->fibNo == ifuFib)
                         {
       		          theIfuFibers->fiberTrans = (refTrans * lineFlux)/
       			                              refLineFlux;
                         }
                       theIfuFibers = theIfuFibers->next;
                      }
       	           }
                 theIfuSlits = theIfuSlits->next;
                }
             }
           ifuQuads = ifuQuads->next;
          }

        ifuImaObj = ifuImaObj->next;
       }

     ifuIma = ifuIma->next;
    }

  
  deleteFloatArray(anX);
  deleteFloatArray(anY);
  deleteFloatArray(skyLines);
  freeFloatVector(coeffs,1,NTERMS);

  deleteObjectTable(ifuImaObjTab);
  return EXIT_SUCCESS;
}


VimosImage *
VmIfuApplyTransmission(VimosImage *imageData, VimosIfuTable *ifuTable, 
                       VimosObjectTable *objectTable, int quadNum,
                       int imageXlen, int imageYlen)
{
  int objNumber;
  int i, k;
  int outXlen, outYlen;
  int refL, refM;

  int nTotGood, nTotDead;

  float refTrans, scal;

  char modName[] = "VmIfuApplyTransmission";
  char comment[80];

  VimosImage *outData;

  VimosIfuQuad *ifuQuads;
  VimosIfuSlit *theIfuSlits;
  VimosIfuFiber *theIfuFibers;
  VimosObjectObject *objects;

  pilMsgInfo (modName, "Apply Relative Transmission Correction");

  /* WARNING: FIBERS TRANSMISSION IS CURRENTLY SCALED TO REFTRANS */
  /* (TRANSMISSION OF THE REFERENCE FIBER): DO WE WANT TO SCALE ALL */
  /* SPECTRA SO TO ALWAYS HAVE TRANSMISSION = 1 INSTEAD OF */
  /* TRANSMISSION = REFTRANSMISSION ? */

  puts("WARNING - fiber transm. scaled to transm. of reference fiber: do");
  puts("          we want to scale to transm. = 1 for all spectra?");

  /* output X size is the same as input */
  outXlen = imageXlen;
  /* Y size is the same as input */
  outYlen = imageYlen;

  /* create output image */
  outData = newImageAndAlloc(outXlen, outYlen);

  /* get the reference fiber L,M coordinates from IFU Table descriptors */
  readIntDescriptor(ifuTable->descs, "ESO PRO REF L", &refL, comment);
  readIntDescriptor(ifuTable->descs, "ESO PRO REF M", &refM, comment);

  /* loop on IFU Table to find the transmission for the reference fiber */
  i = 0;
  refTrans = 0.;
  ifuQuads = ifuTable->quads;
  while(ifuQuads)
    {
     theIfuSlits = ifuQuads->ifuSlits;
     while(theIfuSlits)
       {
	theIfuFibers = theIfuSlits->fibers;
	while(theIfuFibers)
	  {
	   if ( (theIfuFibers->fiberL == refL) && 
		(theIfuFibers->fiberM == refM) )
	     {
              refTrans = theIfuFibers->fiberTrans;
	      i++;
	     }
	   theIfuFibers = theIfuFibers->next;
	  }

        theIfuSlits = theIfuSlits->next;

       }
     ifuQuads = ifuQuads->next;
    }

  /* check that we find just one reference fiber */
  if (i != 1)
    {
     pilMsgError(modName,"Error in selection of reference fiber");
     return NULL;
    }

  /* one never knows... */
  if (refTrans == -1.)
    {
     pilMsgError(modName,"Reference fiber is a dead fiber");
     return NULL;
    }

  /* for each object in object table, look for the transmission in the IFU */
  /* Table and scale the spectrum accordingly to the reference one */

  objects = objectTable->objs;

  /* object counter */
  objNumber = 0;

  /* counter for good and dead fibers on the image */
  nTotGood = 0;
  nTotDead = 0;

  /* re-load the quadrants and select only the right one */
  ifuQuads = ifuTable->quads;
  while(ifuQuads)
    {
     if (ifuQuads->quadNo == quadNum)
       {
        /* if this is the right quadrant, take objects */
        while (objects)
          {
	    /* look at the slits in IFU Table */
           theIfuSlits = ifuQuads->ifuSlits;
           while(theIfuSlits)
             {
	      /* if the IFU slit is the right one the look at the fibers */
              if (objects->IFUslitNo == theIfuSlits->ifuSlitNo)
                {
                 theIfuFibers = theIfuSlits->fibers;
                 while(theIfuFibers)             
                   {
		    /* do not act on dead fibers */
		    if (theIfuFibers->fiberTrans != -1.)
		      {
		       /* if you match the object with the fiber */
		       if (objects->IFUfibNo ==  theIfuFibers->fibNo) 
			 {
			  nTotGood += 1;
			  /* scale transmission to refTrans */
		          scal = refTrans / (theIfuFibers->fiberTrans);

			  /* loop to scale the spectrum */
			  objNumber = objects->rowNum;

			  for (k=0; k<imageXlen; k++)
			    {
			     /* set outData to imageData*scal */
			     outData->data[k + objNumber*imageXlen] = 
			       imageData->data[k + objNumber*imageXlen]*scal;
			    }

			 }
		      } /* end if it is not a dead fiber */
		    else if ( (theIfuFibers->fiberTrans == -1.) &&
			      (objects->IFUfibNo ==  theIfuFibers->fibNo))
		      nTotDead += 1;

		    theIfuFibers = theIfuFibers->next;
                   }
                }

              theIfuSlits = theIfuSlits->next;
             }

           objects = objects->next;
          }
       }

     ifuQuads = ifuQuads->next;
    }

  /* check */
  pilMsgInfo(modName,"nTotGood + nTotDead, %d",(nTotGood + nTotDead) );
  if ( (nTotGood + nTotDead) != 1600)
    {
     pilMsgError(modName,"Wrong number of good + dead fibers");
     return NULL;
    }
  /* debug */
  pilMsgInfo(modName,"N good: %d, N. dead: %d", nTotGood, nTotDead);

  deleteIfuQuad(ifuQuads);
  deleteIfuSlit(theIfuSlits);
  deleteIfuFiber(theIfuFibers);
  deleteObjectObject(objects);

  copyAllDescriptors(imageData->descs, &(outData->descs));
  return outData;
  
}
