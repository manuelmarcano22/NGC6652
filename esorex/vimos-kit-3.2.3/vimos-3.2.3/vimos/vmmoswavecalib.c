/* $Id: vmmoswavecalib.c,v 1.2 2013/03/25 11:43:04 cgarcia Exp $
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
 * $Date: 2013/03/25 11:43:04 $
 * $Revision: 1.2 $
 * $Name:  $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmimage.h"
#include "vmmatrix.h"
#include "vmtable.h"
#include "vmgrismtable.h"
#include "vmextractiontable.h"
#include "vmlinecatalog.h"
#include "vmdistmodels.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmmosutils.h"
#include "vmmoswavecalib.h"
#include "cpl.h"


/**
 * @name vmmoswavecalib MOS Wavelength Calibration
 *
 * @doc
 *   The module collects MOS wavelenght calibration related functions.
 */

/**@{*/

/**
 * @memo
 *   Fit an Inverse Dispersion Relation Model.
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param inputImage      Arc Lamp exposure
 * @param extractionTable Extraction Table
 * @param lineCat         Line Catalog
 * @param arcWidth        Window used for determining the peak of a
 *                        centered arc line (in CCD pixels).
 * @param refineIDS       flag to signal the need to refine the first guess
 *                        IDS from the Extraction Table
 *
 * @doc 
 *   Identify lines in Arc Spectrum, fit dispersion relation for 
 *   every row in te Extraction Table. Calibrate X positioning of slit
 *
 * @author M. Scodeggio
 */

int VmSpDerDispOld(VimosImage *inputImage, VimosExtractionTable *extractionTable, 
                VimosTable *lineCat, int arcWidth, int refineIDS)
{
  const char           modName[] = "VmSpDerDisp";
  int                  i, j, k;
  int                  xOut, yOut;
  int                  numRows, totNumRows;
  int                  numPoints;
  int                  foundLines;
  int                  nPixBelow, nPixAbove;
  int                  imageXlen, imageYlen;
  int                  numLines;
  int                  numLinesDone;
  int                  disOrder;

  int                  widthBuffer = 2;
  
  int                  yExpect;
  float                yExpectFloat;
  int                  validFit, invalidFit, totInvalidFit;
  int                  slitCount = 0;

  float                lambdaCen;
  float                pos;
  /*** float                offset;  ***/
  
  int                  tolMax;
  /***  float                fitTol=2.0; ***/
  float                lambda2pix;
  float                mm2pix;
  float                pixWidth;
  float                lambdaWidth;
  float                nSigma=3.0;
  int                  windowWidth, zeros;
  int                  bigWindowWidth = 30;
  VimosFloatArray     *searchSpectrum;
  
  double              *coeffs = NULL;
  double               datVal;
  double               xOutF;
  double               frac;
  double               fitVariance;

  VimosDpoint         *buffer = NULL;
  VimosDistModel1D    *tunedModel;
  VimosExtractionSlit *slit;
  VimosFloatArray     *spectrum;
  VimosFloatArray     *tmpSpectrum;
  /*** VimosFloatArray     *offsetArray;  ***/
  /*** VimosDistModel2D    *optModX, *optModY;  ***/
  VimosColumn         *wLen;


  char comment[80];

  if (NULL == (searchSpectrum = newFloatArray(bigWindowWidth + 1))) {
    return EXIT_FAILURE;
  }

  cpl_msg_info (modName, "Finding Inverse Dispersion relation Matrix");
  /*
   *  Note: all code related to offsets of the Y coordinate returned
   *  from model are commented with ***'s before and after, to be
   *  easily restored. There is no reason to use offsets, as they
   *  are not saved to header, and the central wavelength has nothing
   *  special with respect to any other wavelength.
   */

  imageXlen = inputImage->xlen;
  imageYlen = inputImage->ylen;

  slit = extractionTable->slits;

  readIntDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("DispersionOrd"), &disOrder, comment);
  readIntDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("NumPixBelow"), &nPixBelow, comment);
  readIntDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("NumPixAbove"), &nPixAbove, comment);
  readFloatDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("WlenCen"), &lambdaCen, comment);
  readFloatDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("OptDistY", 0, 1), &mm2pix, comment);
  readFloatDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("Dispersion", 1, 0, 0), &lambda2pix, comment);

  numPoints = nPixBelow+nPixAbove+1;

  spectrum = newFloatArray(numPoints);
  tmpSpectrum = newFloatArray(numPoints);
  
  numLines = lineCat->cols->len;
  wLen = findColInTab(lineCat, "WLEN");
  
  buffer = newDpoint(numLines);
  
  /* get the number of spectra rows in this Extraction Table */
  totNumRows = numRowsInExtSlits(slit);
  /* create buffer to hold offset for all rows and all skylines */
  /*** offsetArray = newFloatArray(totNumRows); ***/

  numLinesDone = 0;
  totInvalidFit = 0;

  while (slit) {

    /* ALEX: when IFU, do nothing if dead fiber */
    if (slit->IFUfibTrans != -1.)
      {

       numRows = slit->numRows;

       pixWidth = mm2pix * slit->width + widthBuffer;
       lambdaWidth = pixWidth / lambda2pix;
       if (arcWidth < 0) windowWidth = 3*pixWidth;
       else              windowWidth = arcWidth;
       
       /*** Begin duplicated part ***/
       i = numRows/2;
       if (refineIDS)
	 {
	   /*
	    *  Extract spectrum from central row (far from possible 
	    * contaminations)
	    *  for "tuning" the first guess.
	    */
      
	   slit->invDisQuality->data[i] = 1;
	   /* extract spectrum for this row */
	   for (j = -nPixBelow; j <= nPixAbove; j++) {
	     /* compute Y-CCD-pixel of spectrum */
	     yOut = slit->ccdY->data[i] + j;
	     /* compute X-pixel of spectrum */
	     xOutF = slit->ccdX->data[i] + computeDistModel1D(slit->crvPol[i],
							      yOut);
	     xOut = xOutF;  /* Make it integer (truncation) */
	     if (xOut >= 0 && xOut+1 < imageXlen) {
	       if (yOut >= 0 && yOut < imageYlen) {
		 /* simple linear interpolation */
		 frac = xOutF-xOut;
		 datVal= ( (1.0-frac)*inputImage->data[xOut + yOut*imageXlen] +
			    frac*inputImage->data[xOut+1+yOut*imageXlen] );
	       }
	       else {
		 /*
		  *  The spectrum is dispersed beyond the CCD frame.
		  *  Fill with zeroes...
		  */
		 datVal = 0.;
	       }
	       /* store in buffer */
	       spectrum->data[nPixBelow+j] = datVal;
	     }
	     else {
	       /*
		*  Mark as invalid the part of the slit extending beyond the
		*  CCD frame.
		*/
	       slit->invDisQuality->data[i] = 0;
	       break;
	     }
	   }

	   /*
	    *  NEW: pre-tuning of first guess.   (C.Izzo)
	    */
      
	   if (slit->invDisQuality->data[i]) {
	     tunedModel = findMaxMatchIndex(slit->invDis[i], lineCat, 
					    lambdaWidth, spectrum, nPixBelow);
	   }
	   else {
	     tunedModel = newDistModel1D(slit->invDis[i]->order);
	     tunedModel->offset = slit->invDis[i]->offset;
	     for(j=0; j<=slit->invDis[i]->order; j++)
	       tunedModel->coefs[j] = slit->invDis[i]->coefs[j];
	   }
	 }
       /***  End of duplicated part ***/
       else {
	 tunedModel = newDistModel1D(slit->invDis[i]->order);
	 tunedModel->offset = slit->invDis[i]->offset;
	 for(j=0; j<=slit->invDis[i]->order; j++)
	   tunedModel->coefs[j] = slit->invDis[i]->coefs[j];
       }

       invalidFit = 0;
       slitCount++;
       
       for (i = 0; i < numRows; i++) 
	 {
	   /* TBD: checks on convergence!!! */
	   
	   /* extract spectrum for this row */
	   for (j = -nPixBelow; j <= nPixAbove; j++) {
	     /* compute Y-CCD-pixel of spectrum */  
	     yOut = slit->ccdY->data[i] + j;
	     /* compute X-pixel of spectrum */
	     xOutF = slit->ccdX->data[i] + computeDistModel1D(slit->crvPol[i], 
							      yOut);
	     xOut = xOutF;  /* Make it integer (truncation) */
	     if (xOut >= 0 && xOut+1 < imageXlen) {
	       if (yOut >= 0 && yOut < imageYlen) {
		 /* simple linear interpolation */
		 frac = xOutF-xOut;
		 datVal= ( (1.0-frac)*inputImage->data[xOut + yOut*imageXlen] +
			    frac*inputImage->data[xOut+1+yOut*imageXlen] );
	       }
	       else {
		 /*
		  *  Fill with zeroes parts of the spectrum dispersed beyond
		  *  the CCD frame.
		  */
		 datVal = 0.;
	       }
	       /* store in buffer */
	       spectrum->data[nPixBelow+j] = datVal;
	     }
	     else {
	       /*
		*  Mark as invalid the part of the slit extending beyond the
		*  CCD frame.
		*/
	       slit->invDisQuality->data[i] = 0;
	       break;
	     }
	   }
	   /*
	    *  Go immediately to next slit row if current row is outside CCD
	    */
	   if (slit->invDisQuality->data[i] == 0) continue;

	   /* 
	    *  Loop over lines in Line Catalog
	    */
	   foundLines = 0;
	   for (j = 0 ; j < numLines; j++) {
	     /*
	      * Compute expected offset from slit position of location of 
	      * arc line j 
	      */
	     
	     yExpectFloat = computeDistModel1D(tunedModel, 
					       wLen->colValue->fArray[j]);
	     yExpect = yExpectFloat;

	     if (refineIDS) {
	     /*
	      *  Check whether yExpect and the search interval are included
	      *  in the extracted spectral range: if not, skip current line
	      */
	       if ((yExpect+nPixBelow-bigWindowWidth/2 < 0) ||
		  (yExpect+nPixBelow+bigWindowWidth/2 > numPoints)) continue;
	       /*
		*  Extract a wide portion of spectrum around the expected
		*  position, were to look for the closest peak.
		*/

	       for (k = 0; k <= bigWindowWidth; k++) {
		 searchSpectrum->data[k] = 
		   log10(MAX(spectrum->data[yExpect+nPixBelow+1-
					   bigWindowWidth/2+k], 1.));
	       }
	       
	       /*
		*  Upgrade the expected position
		*/
	  
	       yExpect += findClosestPeak(searchSpectrum->data, bigWindowWidth)
		 - bigWindowWidth/2;
	     }

	     /* 
	      *  Extract a smaller window, centered on this raw position
	      *  of the closest peak, and find peak's accurate position.
	      */
	     if ((yExpect+nPixBelow-windowWidth/2 < 0) ||
		 (yExpect+nPixBelow+windowWidth/2 > numPoints)) continue;

	     for (k = 0; k <= windowWidth; k++) {
	       tmpSpectrum->data[k] = 
		 spectrum->data[yExpect+nPixBelow+1-windowWidth/2 + k];
	     }

             /*
              * Added to check that we are not close to the edge of a 
              * truncated spectrum
              */

             zeros = 0;
             for (k = 0; k <= windowWidth; k++)
               if (fabs(tmpSpectrum->data[k]) < 0.000000001)
                 zeros++;

             if (zeros > 0)
               continue;          /* Line not found */

	     if (VM_TRUE == findPeak1D(tmpSpectrum->data, windowWidth+1, 
				       &pos, 2)) {
	       pos += yExpect-windowWidth/2;
	     }
	     else {
	       continue;          /* Line not found */
	     }
        

	     /* 
	      *  Add position of arc line in buffer 
	      */
	     buffer[foundLines].x = wLen->colValue->fArray[j] - lambdaCen;
	     buffer[foundLines].y = pos;
	     foundLines++;
	   }

	   /* 
	    *  Free old array of coefficients, if exists 
	    */
	   if (coeffs != NULL) {
	     cpl_free(coeffs);
	     coeffs = NULL;
	   }

	   /* 
	    *  Fit new inverse dispersion relation 
	    */

	   do {
	     tolMax = 0;
	     validFit = 1;
	     coeffs=fit1DPoly(slit->invDis[i]->order,buffer,foundLines,
			      &fitVariance);
	     /*	printf("%d %d %d %f \n", slit->slitNo, i, foundLines, 
		sqrt(rmsError));*/
	     if (coeffs == NULL) {
/*
	       cpl_msg_warning(modName, "One invalid fit");
*/
	       slit->invDisQuality->data[i] = validFit = 0;
               if (foundLines >= 2*slit->invDis[i]->order) {
                 invalidFit++;
               }
	       break;
	     }

	     /* 
	      *  Copy coefficients in inverse dispersion relation for this row
	      */
	     for (k = 0; k <= slit->invDis[i]->order; k++) {
	       slit->invDis[i]->coefs[k] = coeffs[k];
	     }
	     slit->invDisRms->data[i] = sqrt(fitVariance);
	     slit->invDis[i]->offset = 0.0;
	     /* cull bad lines */
	     for (k = 0; k < foundLines; k++) {
	       yExpectFloat = computeDistModel1D(slit->invDis[i], buffer[k].x);
	       yExpect = yExpectFloat;
	       if( fabs(yExpectFloat-buffer[k].y) > nSigma*sqrt(fitVariance) ){
		 tolMax = 1;
		 if (k != foundLines-1) {
		   for (j = k; j < foundLines-1; j++) {
		     buffer[j].x = buffer[j+1].x;
		     buffer[j].y = buffer[j+1].y;
		   }
		   k--;
		 }
		 foundLines--;
	       }
	     }
	     /* until all lines are ok */
	   } while ( (foundLines > slit->invDis[i]->order+2) && (tolMax > 0) );

	   if (foundLines < 2*slit->invDis[i]->order) {
/*
	     cpl_msg_warning(modName, 
			   "Too few lines for a valid fit: fit rejected!");
*/
	     slit->invDisQuality->data[i] = validFit = 0;
             invalidFit++;
	   }


	   /* store wavelength shift on CCD */
	   slit->invDis[i]->offset = lambdaCen;
	   
	   if (validFit) {
	     /*** offsetArray->data[numLinesDone] =  
		  computeDistModel1D(slit->invDis[i], lambdaCen);  ***/
	     numLinesDone++;
	   }
	 }

       if (invalidFit) {
         totInvalidFit += invalidFit;
         cpl_msg_warning(modName, "Slit %d: %d of %d rows IDS fit failure.", 
                       slitCount, invalidFit, numRows);
       }

       if (refineIDS)
	 deleteDistModel1D(tunedModel);  /***  Added  ***/

      } /* ALEX: end if dead fiber */

    slit = slit->next;
  }

  if (totInvalidFit == totNumRows) {
    cpl_msg_error(modName, "The IDS fit failed in all slit rows!");
    return EXIT_FAILURE;
  }
  else {
    if (totInvalidFit) {
      cpl_msg_warning(modName, "All slits: %d of %d rows IDS fit failure.", 
                    totInvalidFit, totNumRows);
    }
  }

  /* determine median offset for CCD */
  /*** offset = medianWirth(offsetArray->data, numLinesDone); ***/

  /* Test: remove the offset!!!  (C.Izzo)  */
  /*** offset = 0.;  ***/

  /* get parameters of Optical Distortion Model from Extraction Table */
  /*** readOptDistModel(extractionTable->descs, &optModX, &optModY);  ***/
  /* add offset */
  /*** optModY->coefs[0][0] += offset;  ***/
  /* and write coefficients back to Table */
  /*** writeOptDistModel(&(extractionTable->descs), optModX, optModY);  ***/

  /* loop over slits to update slit positions */
  /*** slit = extractionTable->slits; ***/

  /*** while (slit) {

       * ALEX: the update of crvPol and invDis is done ALSO for dead fibers *
       * (no IF on IFUfibTrans) as in VmSpCalShifts *
    numRows = slit->numRows; ***/
    /* loop over all rows of a slit */
    /*** for (i = 0; i < numRows; i++) { ***/
      /* add offset to slit position */
 /*** slit->ccdY->data[i] += offset;
      slit->crvPol[i]->offset += offset;
      slit->invDis[i]->coefs[0] -= offset;
    }

    slit = slit->next;
  }

  deleteFloatArray(searchSpectrum);

  writeIntDescriptor(&(extractionTable->descs),
                     pilTrnGetKeyword("NumPixBelow"), nPixBelow + offset, "");
  writeIntDescriptor(&(extractionTable->descs),
                     pilTrnGetKeyword("NumPixAbove"), nPixAbove - offset, "");
  ***/

  return EXIT_SUCCESS;
}

/**
 * @memo
 *   Fit an Inverse Dispersion Relation Model.
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param inputImage      Arc Lamp exposure
 * @param extractionTable Extraction Table
 * @param lineCat         Line Catalog
 * @param arcWidth        Window used for determining the peak of a
 *                        centered arc line (in CCD pixels).
 * @param refineIDS       flag to signal the need to refine the first guess
 *                        IDS from the Extraction Table
 *
 * @doc 
 *   Identify lines in Arc Spectrum, fit dispersion relation for 
 *   every row in te Extraction Table. Calibrate X positioning of slit
 *
 * @author C. Izzo
 */

int VmSpDerDisp(VimosImage *inputImage, VimosExtractionTable *extractionTable, 
                VimosTable *lineCat, int arcWidth, float level)
{
  const char           modName[] = "VmSpDerDisp";
  int                  i, j, k;
  int                  xOut, yOut;
  int                  numRows, totNumRows;
  int                  numPoints;
  int                  foundLines;
  int                  nPixBelow, nPixAbove;
  int                  imageXlen, imageYlen;
  int                  numLines;
  int                  disOrder;

  int                  widthBuffer = 0;
  
  int                  validFit, invalidFit, totInvalidFit;
  int                  slitCount = 0;

  float                lambdaCen, lambdaChecklo, lambdaCheckhi;
  
  float                lambda2pix;
  float                mm2pix;
  float                pixWidth;
/*  float                lambdaWidth; */
  float                pixelCen, pixelLo, pixelHi;
  float               *check[3];
/*  int                  windowWidth; */
  
  double              *coeffs = NULL;
  double               datVal;
  double               xOutF;
  double               frac;
  double               fitVariance;

  int                  npeaks;
  double              *peaks;
  double              *lines;
  double             **output;
  double               min_disp;
  double               max_disp;
  double               min, max, span;
  double               tolerance = 0.05;

  VimosDpoint         *buffer = NULL;
  VimosExtractionSlit *slit;
  VimosFloatArray     *spectrum;
  VimosColumn         *wLen;


  char comment[80];

  cpl_msg_info (modName, "Finding Inverse Dispersion relation Matrix");

  /*
   *  Note: all code related to offsets of the Y coordinate returned
   *  from model are commented with ***'s before and after, to be
   *  easily restored. There is no reason to use offsets, as they
   *  are not saved to header, and the central wavelength has nothing
   *  special with respect to any other wavelength.
   */

  imageXlen = inputImage->xlen;
  imageYlen = inputImage->ylen;

  slit = extractionTable->slits;

  readIntDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("DispersionOrd"), &disOrder, comment);
  readIntDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("NumPixBelow"), &nPixBelow, comment);
  readIntDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("NumPixAbove"), &nPixAbove, comment);
  readFloatDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("WlenCen"), &lambdaCen, comment);
  readFloatDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("OptDistY", 0, 1), &mm2pix, comment);
  readFloatDescriptor(extractionTable->descs, 
              pilTrnGetKeyword("Dispersion", 1, 0, 0), &lambda2pix, comment);

  max_disp = min_disp = 1 / lambda2pix;
  max_disp += max_disp / 5.5;
  min_disp -= min_disp / 5.5;

  numPoints = nPixBelow+nPixAbove+1;

  spectrum = newFloatArray(numPoints);
  
  numLines = lineCat->cols->len;
  wLen = findColInTab(lineCat, "WLEN");
  lines = cpl_malloc(numLines * sizeof(double));

  for (j = 0; j < numLines; j++)
    lines[j] = wLen->colValue->fArray[j];

  buffer = newDpoint(numLines);
  
  /* get the number of spectra rows in this Extraction Table */
  totNumRows = numRowsInExtSlits(slit);

  totInvalidFit = 0;

  while (slit) {

       numRows = slit->numRows;

       pixWidth = mm2pix * slit->width + widthBuffer;
  /*     lambdaWidth = pixWidth / lambda2pix; 
       if (arcWidth < 0) windowWidth = 3*pixWidth;
       else              windowWidth = arcWidth; */

       min = max = slit->ccdY->data[0];
       for (i = 1; i < numRows; i++) {
         if (min > slit->ccdY->data[i])
           min = slit->ccdY->data[i];
         if (max < slit->ccdY->data[i])
           max = slit->ccdY->data[i];
       }

       span = (max - min) * mm2pix + 3;
       
       invalidFit = 0;
       slitCount++;
       
       for (i = 0; i < numRows; i++) 
	 {
	   
           slit->invDisQuality->data[i] = 1;

	   /* extract spectrum for this row */
	   for (j = -nPixBelow; j <= nPixAbove; j++) {
	     /* compute Y-CCD-pixel of spectrum */  
	     yOut = slit->ccdY->data[i] + j;
	     /* compute X-pixel of spectrum */
	     xOutF = slit->ccdX->data[i] + computeDistModel1D(slit->crvPol[i], 
							      yOut);
	     xOut = xOutF;  /* Make it integer (truncation) */
	     if (xOut >= 0 && xOut+1 < imageXlen) {
	       if (yOut >= 0 && yOut < imageYlen) {
		 /* simple linear interpolation */
		 frac = xOutF-xOut;
		 datVal= ( (1.0-frac)*inputImage->data[xOut + yOut*imageXlen] +
			    frac*inputImage->data[xOut+1+yOut*imageXlen] );
	       }
	       else {
		 /*
		  *  Fill with zeroes parts of the spectrum dispersed beyond
		  *  the CCD frame.
		  */
		 datVal = 0.;
	       }
	       /* store in buffer */
	       spectrum->data[nPixBelow+j] = datVal;
	     }
	     else {
	       /*
		*  Mark as invalid the part of the slit extending beyond the
		*  CCD frame.
		*/
	       slit->invDisQuality->data[i] = 0;
	       break;
	     }
	   }
	   /*
	    *  Go immediately to next slit row if current row is outside CCD
	    */
	   if (slit->invDisQuality->data[i] == 0) continue;

           peaks = collectPeaks(spectrum->data, numPoints, 
                                level, pixWidth, &npeaks);

           if (peaks) {

             output = identPeaks(peaks, npeaks, lines, numLines,
                      min_disp, max_disp, tolerance, &foundLines);

             if (output) {
               for (j = 0; j < foundLines; j++) {
	         buffer[j].x = output[1][j] - lambdaCen;
	         buffer[j].y = output[0][j] - nPixBelow - 1;
	       }
               cpl_free(output[0]);
               cpl_free(output[1]);
               cpl_free(output);
             }

             cpl_free(peaks);

           }
           else
             foundLines = 0;

	   validFit = 0;
           if (foundLines >= 1.5*slit->invDis[i]->order) {

	     /* 
	      *  Free old array of coefficients, if exists 
	      */

	     if (coeffs != NULL) {
	       cpl_free(coeffs);
	       coeffs = NULL;
	     }

	     /* 
	      *  Fit new inverse dispersion relation 
	      */

	     coeffs = fit1DPoly(slit->invDis[i]->order, buffer, foundLines,
                                &fitVariance);

	     if (coeffs != NULL && fitVariance < 1.0) {

               /*
                *  Copy coefficients of dispersion relation for this row
                */

               for (k = 0; k <= slit->invDis[i]->order; k++)
                 slit->invDis[i]->coefs[k] = coeffs[k];

               slit->invDisRms->data[i] = sqrt(fitVariance);
               validFit = 1;
             }
           }

           if (!validFit) {
	     slit->invDisQuality->data[i] = 0;
             invalidFit++;
           }

	   /* store wavelength shift on CCD */
	   slit->invDis[i]->offset = lambdaCen;
	   
	 }

       /*
        * Check if all fits are belonging to the same slit (it may not
        * be the case, due to misalignments).
        */

       check[0] = cpl_malloc(numRows * sizeof(float));
       check[1] = cpl_malloc(numRows * sizeof(float));
       check[2] = cpl_malloc(numRows * sizeof(float));

       lambdaChecklo = lambdaCen - (nPixBelow / lambda2pix) / 2;
       lambdaCheckhi = lambdaCen + (nPixAbove / lambda2pix) / 2;

       for (j = 0, k = 0; k < numRows; k++) {
         if (slit->invDisQuality->data[k] != 0) {
           check[0][j] = computeDistModel1D(slit->invDis[k], lambdaChecklo);
           check[1][j] = computeDistModel1D(slit->invDis[k], lambdaCen);
           check[2][j] = computeDistModel1D(slit->invDis[k], lambdaCheckhi);
           ++j;
         }
       }

       if (j > 0) {
         pixelLo = medianPixelvalue(check[0], j);
         pixelCen = medianPixelvalue(check[1], j);
         pixelHi = medianPixelvalue(check[2], j);

         if (fabs(pixelCen) > 30) {
           for (k = 0; k < numRows; k++) {
             if (slit->invDisQuality->data[k] != 0) {
               slit->invDisQuality->data[k] = 0;
               ++invalidFit;
             }
           }
         }

         for (j = 0, k = 0; k < numRows; k++) {
           if (slit->invDisQuality->data[k] != 0) {
             if (fabs(check[0][j] - pixelLo) > span
               || fabs(check[1][j] - pixelCen) > span
               || fabs(check[2][j] - pixelHi) > span) {
               slit->invDisQuality->data[k] = 0;
               ++invalidFit;
             }
             ++j;
           }
         }
       }

       cpl_free(check[0]);
       cpl_free(check[1]);
       cpl_free(check[2]);

       if (invalidFit) {
         totInvalidFit += invalidFit;
         cpl_msg_warning(modName, "Slit %d: %d of %d rows IDS fit failure.", 
                       slitCount, invalidFit, numRows);
       }

    slit = slit->next;
  }

  deleteDpoint(buffer);
  deleteFloatArray(spectrum);
  cpl_free(lines);

  if (totInvalidFit == totNumRows) {
    cpl_msg_error(modName, "The IDS fit failed in all slit rows!");
    return EXIT_FAILURE;
  }
  else {
    if (totInvalidFit) {
      cpl_msg_warning(modName, "All slits: %d of %d rows IDS fit failure.", 
                    totInvalidFit, totNumRows);
    }
  }

  return EXIT_SUCCESS;
}


/**
 * @memo
 *   Fit a Inverse Dispersion Relation to wavelength calibrations
 *   obtained for each row of each slit spectrum.
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
 *   of the polynomials fitting the wavelength/pixel relations
 *   for each spectra found in an Arc Lamp Exposure.
 *
 * @author M. Scodeggio
 */

int VmSpDispMatrix(VimosExtractionTable *extractionTable,
                   VimosTable *grismTable, int grismFlag)
{
  int     i, j; 
  int     numSlits;
  int     slitDone;
  int     polOrder, xOrder, yOrder;
  int     numRows, numGoodRows;
  double  rmsFit;
  double  minSlity, maxSlity;
  double *coeff;
  double *c = NULL;
  double  medianCoeff;
  int     xord = 2;

  VimosExtractionSlit *slit;
  VimosPixel *coefArray;
  VimosDpoint *list;
  VimosDistModelFull *idsMat;
  
  char comment[80];
  
  numSlits = 0;
  slit = extractionTable->slits;

  while (slit) {
    numSlits++;
    slit = slit->next;
  }
  
  coefArray = newPixel(2*numSlits);
  list = newDpoint(2*numSlits);
  readIntDescriptor(extractionTable->descs, pilTrnGetKeyword("DispersionOrd"),
                    &polOrder, comment);
  readIntDescriptor(extractionTable->descs, pilTrnGetKeyword("DispersionOrdX"),
                    &xOrder, comment);
  readIntDescriptor(extractionTable->descs, pilTrnGetKeyword("DispersionOrdY"),
                    &yOrder, comment);
    
  /* note that for the moment we keep the order in X equal to the one in Y!! */
  idsMat = newDistModelFull(polOrder, xOrder, xOrder);

  slit = extractionTable->slits;
  minSlity = maxSlity = slit->maskY->data[0];
  while (slit) {
    if (slit->maskY->data[0] > maxSlity)
      maxSlity = slit->maskY->data[0];
    if (slit->maskY->data[0] < minSlity)
      minSlity = slit->maskY->data[0];
    slit = slit->next;
  }
  
  for (i = 0; i <= polOrder; i++) {

    slit = extractionTable->slits;
    slitDone = 0;
    while (slit) {

      numRows = slit->numRows;
      coeff = (double *) cpl_malloc(numRows*sizeof(double));

      /* 
       *  Skip columns where the solution is not available 
       */

      for (numGoodRows = 0, j = 0; j < numRows; j++) {
        if (slit->invDisQuality->data[j] != 0) {
          coeff[numGoodRows] = slit->invDis[j]->coefs[i];
          numGoodRows++;
        }
      }

      if (numGoodRows) {
        medianCoeff = medianDouble(coeff, numGoodRows);

        if (maxSlity - minSlity > 1.0) {
          coefArray[slitDone].x = slit->maskX->data[numRows/2];
          coefArray[slitDone].y = slit->maskY->data[numRows/2];
          coefArray[slitDone].i = medianCoeff;
        }
        else {
          list[slitDone].x = slit->maskX->data[numRows/2];
          list[slitDone].y = medianCoeff;
        }
        slitDone++;

        cpl_free(coeff);
      }
      slit = slit->next;

    }

    if (maxSlity - minSlity > 1.0) {

      /* ugly!!!! */
      deleteDistModel2D(idsMat->coefs[i]);
  
      /*
       *  Note: the offsets are always set to zero.
       */

      if (fitDistModel2D(coefArray, slitDone, xOrder, 0.0, 0.0,
                     &(idsMat->coefs[i]), &rmsFit) == VM_FALSE)
        return EXIT_FAILURE;
    }
    else {
      c = fit1DPoly(xord, list, slitDone, NULL);
      if (c) {
        for (j = 0; j <= xord; j++)
          idsMat->coefs[i]->coefs[j][0] = c[j];
        free(c);
        c = NULL;
      }
      else {
        deleteDpoint(list);
        deletePixel(coefArray);
        deleteDistModelFull(idsMat);
        return EXIT_FAILURE;
      }
    }
  }

  writeInvDispMatrix(&(extractionTable->descs), idsMat);

  if (grismFlag)
    writeInvDispMatrix(&(grismTable->descs), idsMat);
    
  deleteDpoint(list);
  deletePixel(coefArray);
  deleteDistModelFull(idsMat);

  return EXIT_SUCCESS;
}


int VmSpCalShifts(VimosImage *inputImage, VimosTable *grismTable, 
		  VimosExtractionTable *extractionTable, int grismFlag, 
		  int lineWidth, int fuzz)
{
  int i, j, k;
  int jMin, jMax;
  int xOut, yOut;
  int numRows;
  int numPoints;
  int numSlits;
  int bigWindowWidth = 30;
  int windowWidth = 20;
  VimosUlong32 index;
  
  int nPixBelow, nPixAbove;
  int imageXlen, imageYlen;
  int numLinesDone;
  int numSkyLines;
  int yExpect;
  int outOfImage = 0;
  
  float pos;
  float offSet;
  float yExpectFloat;
  float xStart, xEnd;
  VimosFloatArray     *searchSpectrum;

  double upVal, loVal;
  double datVal;
  double xOutF;
  double frac;
  
  VimosExtractionSlit *slit;
  VimosFloatArray      *spectrum;
  VimosFloatArray      *tmpSpectrum;
  VimosFloatArray      *skyLines;
  VimosFloatArray      *offsetArray;
  VimosFloatArray      *offsetYArray;
  VimosFloatArray      *tmpArray;
  VimosDistModel2D     *optModX, *optModY;
  
 
  char modName[] = "VmSpCalShifts";

  /*ALEX */
  char insMode[80];

  cpl_msg_debug(modName,"Calibrating shifts");

  if (NULL == (searchSpectrum = newFloatArray(bigWindowWidth + 1))) {
    return EXIT_FAILURE;
  }

 
  imageXlen = inputImage->xlen;
  imageYlen = inputImage->ylen;
  slit = extractionTable->slits;
  
  if (lineWidth > 0) windowWidth = lineWidth;

  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN LO", 
                    &nPixBelow, NULL);
  readIntDescriptor(extractionTable->descs, "ESO PRO SPECT LLEN HI",
                    &nPixAbove, NULL); 
  readSkyLines(extractionTable->descs, &numSkyLines, &skyLines);
  
  /* length in pixels of an extracted spetral row */
  numPoints = nPixBelow+nPixAbove+1;


  /*ALEX: determine numSlit depending on MOS or IFU: if IFU, numSlit equal to
    number of good fibers (to prevent zeros in medianWirth on offsetYArray */
  readStringDescriptor(inputImage->descs, pilTrnGetKeyword("InstrumentMode"), 
		       insMode, NULL);
  if (!strncmp(insMode,"MOS",3))
    {
      /* get the number of slits in this Extraction Table */
      numSlits = numSlitsInExtTable(extractionTable);
    }
  else if (!strncmp(insMode,"IFU",3))
    {
      numSlits = 0;
      while (slit) {
	if ( (slit->IFUfibTrans) >= 0.0)
	  {
	    numSlits++;
	  }
	slit = slit->next;
      }
    }
  else
    {
      cpl_msg_error(modName,"Unable to determine Instrument Mode");
      return EXIT_FAILURE;
    }
  /* ALEX end */

  /* create buffers to hold extracted spectral row */
  spectrum = newFloatArray(numPoints);
  tmpSpectrum = newFloatArray(numPoints);

  offsetYArray = newFloatArray(numSlits);

  /* ALEX: restart */
  slit = extractionTable->slits;

  /* loop over all slits: determine offset in Y0  */
  while (slit) {

    /* ALEX: when IFU, do the following except for dead fibers */
    if ( (slit->IFUfibTrans) >= 0.0) 
      {
       /* get number of rows for this slit */
	numRows = slit->numRows;
	/* create buffer to hold offset for all rows and all skylines */
	offsetArray = newFloatArray(numRows*numSkyLines);
	/* done nothing sofar */  
	numLinesDone = 0;

	/* loop over all rows of this slit */
	for (i = 0; i < numRows; i++) {
	  /* check that this row has a wavelength solution */
	  if (slit->invDisQuality->data[i] == 0)
	    continue;
	  /* extract spectrum for this row  */
	  for (j = -nPixBelow; j <= nPixAbove; j++) {
	    /* compute Y-pixel of spectrum */
	    yOut = slit->ccdY->data[i] + j;
	    /* compute X-pixel of spectrum */
	    xOutF = slit->ccdX->data[i]+computeDistModel1D(slit->crvPol[i],yOut);
	    xOut = xOutF;
	    if (yOut >= 0 && yOut < imageYlen
             && xOut > 1 && xOut < imageXlen - 1) {
	      /* simple linear interpolation */
	      frac = xOutF-xOut;
	      datVal = ( (1.0-frac)*inputImage->data[xOut + yOut*imageXlen] +
			 frac*inputImage->data[xOut+1+yOut*imageXlen] );
	    }
	    else {
	      /*
	       *  The spectrum is dispersed beyond the CCD frame.
	       *  Fill with zeroes...
	       */
	      datVal = 0.;
	    }
	    /* store in buffer */
	    spectrum->data[nPixBelow+j] = datVal;
	    
	  }
	  
	  /* loop over skylines */
	  for (j = 0 ; j < numSkyLines; j++) {
	    /* compute expected offset from slit position of location of sky */
	    /* line j */
	    yExpectFloat = computeDistModel1D(slit->invDis[i], skyLines->data[j]);
	    yExpect = yExpectFloat;
	    
	    /*
	     *  Check whether yExpect and the search interval are included
	     *  in the extracted spectral range: if not, skip current line
	     */
	    if ((yExpect+nPixBelow-windowWidth/2 < 0) ||
		(yExpect+nPixBelow+windowWidth/2 > numPoints)) continue;
	    
	    /* 
	     *  Extract a window, centered on this raw position
	     *  of the closest peak, and find peak's accurate position.
	     */
	    for (k = 0; k <= windowWidth; k++) {
	      tmpSpectrum->data[k] = 
		spectrum->data[yExpect+nPixBelow+1-windowWidth/2 + k];
	    }
	    if(VM_TRUE == findPeak1D(tmpSpectrum->data,windowWidth+1,&pos,2)) {
	      pos += yExpect-windowWidth/2;
	    }
	    else {
	      continue;          /* Line not found */
	    }
	    
	    /* add position of sky line in buffer */
	    offsetArray->data[numLinesDone] = pos-yExpectFloat;
	    numLinesDone++;        
	  }            
	}
	/* determine median offset in Y for Slit */

        offSet = medianWirth(offsetArray->data, numLinesDone);

	/* loop over slit to update model positions */
	/* ALEX: this updating is done for good fibers only */

	for (i = 0; i < numRows; i++) {
	  slit->ccdY->data[i] += offSet;
	  slit->crvPol[i]->offset += offSet;
	}
        offsetYArray->data[slit->slitNo - 1] = offSet;

      } /* ALEX: end if is a good fiber */

    /* ALEX: if is a dead fiber set offset to 0 */
    else if ( (slit->IFUfibTrans) < 0.0)
      {
	offSet = 0.;
      }

    /* process next slit */    
    slit = slit->next;
  }

  offSet =  medianWirth(offsetYArray->data, numSlits);
  deleteFloatArray(offsetYArray);

  /* get parameters of Optical Distortion Model from Extraction Table */
  readOptDistModel(extractionTable->descs, &optModX, &optModY);
  /* add offset */
  optModY->coefs[0][0] += offSet;

  writeIntDescriptor(&(extractionTable->descs), "ESO PRO SPECT LLEN LO",
		     nPixBelow + offSet, "");
  writeIntDescriptor(&(extractionTable->descs), "ESO PRO SPECT LLEN HI",
		     nPixAbove - offSet, "");

  /* now find offset in X */

  /* create buffer to hold offset for all slits  and all skylines */
  deleteFloatArray(offsetArray);
  offsetArray = newFloatArray(numSlits*2*numSkyLines);
  /* reset counter */
  numLinesDone = 0;
  
  /* first slit in table */
  slit = extractionTable->slits;
  /* loop over all  slits */
  while (slit) {

    /* spatial length of slit */
    numRows = slit->numRows;
    jMin = 0;
    jMax = numRows-1;


    while ((slit->invDisQuality->data[jMin] == 0) && jMin < (numRows-1))
      jMin++;
    while ((slit->invDisQuality->data[jMax] == 0) && jMax > 0)
      jMax--;

    /* create temp buffer */
    tmpArray = newFloatArray(numRows + 2*fuzz);


    /* ALEX: when IFU, do the following except for dead fibers */
    if ( (slit->IFUfibTrans) >= 0.0) 
      {

	/* loop over 'wavelength' */
	for (j = 0; j < numSkyLines; j++) {
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
	    if(yOut < 0 || yOut >=imageYlen || xOut < 0 || xOut >= imageXlen) {
	      outOfImage = 1; 
	      continue;
	    } 
	    outOfImage = 0;
	    numPoints++;
	    
	    /* quick interpolation of data */
	    frac = xOutF-xOut;
	    index = xOut+yOut*imageXlen;
	    datVal = ( (1.0-frac)*inputImage->data[index] +
		       frac*inputImage->data[index+1] );
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
	    findSpectrumBorders(tmpArray, &upVal, &loVal, fuzz); 
	    /* add position of sky line in buffer */
            if (upVal > 0.0 && loVal > 0.0) {
	      offsetArray->data[numLinesDone] = loVal - fuzz;
	      numLinesDone++;
	      offsetArray->data[numLinesDone] = upVal - fuzz - numRows;
	      numLinesDone++;
            }
	  }
	}
	deleteFloatArray(tmpArray);

      } /* ALEX: end if on dead fibers */

    slit = slit->next;
  }
  /* determine median offset in X for CCD */
  offSet = medianWirth(offsetArray->data, numLinesDone);

  /* add offset */
  optModX->coefs[0][0] += offSet;
  /* and write coefficients back to Table */
  writeOptDistModel(&(extractionTable->descs), optModX, optModY);

  /* loop over slits to update slit positions */

  /* ALEX: this updating is done for good fibers only */

  slit = extractionTable->slits;
  while (slit) {
    if ( (slit->IFUfibTrans) >= 0.0) 
      {
	numRows = slit->numRows;
	/* loop over all rows of a slit */
	for (i = 0; i < numRows; i++) {
	  /* add offset to slit position */
	  slit->ccdX->data[i] += offSet;
	}
      }
    slit = slit->next;
  }
 
  if (grismFlag) {
    
    writeOptDistModel(&(grismTable->descs), optModX, optModY);
    
  }
 
  deleteFloatArray(spectrum);
  /*  deleteFloatArray(tmpSpectrum);*/
  deleteFloatArray(offsetArray);
  deleteFloatArray(skyLines);
  deleteDistModel2D(optModX);
  deleteDistModel2D(optModY);
 
  return EXIT_SUCCESS;
}


/**
 * @memo 
 *   Build the list of wavelength intervals where the spectrum to be
 *   calibrated will be integrated.
 *
 * @return List of wavelength intervals.
 * 
 * @param lineCat Line catalogue.
 * @param width   Slit width (in the same wavelength units used in 
 *                line catalog).
 *   
 * @doc
 *   The wavelength intervals where a spectrum will be integrated
 *   depend on the slit width and on the wavelength of the lines
 *   listed in the catalogue. When intervals centered on different 
 *   lines overlap, they are merged into a single interval. The 
 *   lines in catalog must be listed with increasing wavelength.
 *
 * @author C.Izzo
 */

VimosDpoint *getWavIntervals(VimosTable *lineCat, float width)
{
  const char   modName[]    = "getWavIntervals";
  VimosColumn *wLen         = findColInTab(lineCat, "WLEN");
  VimosDpoint *wavIntervals = NULL;
  VimosDpoint *interval;
  double      *start;
  double      *end;
  int          numLines;
  int          numIntervals = 0;
  int          i            = 0;

  if (wLen) {
    numLines = lineCat->cols->len;
    start    = (double *) cpl_malloc(numLines*sizeof(double));
    end      = (double *) cpl_malloc(numLines*sizeof(double));
   /*
    *  First interval
    */
    start[0] = wLen->colValue->fArray[i] - width/2;
    end[0]   = wLen->colValue->fArray[i] + width/2;
   /*
    *  Note: current interval just get longer if next interval
    *  begins inside it.
    */
    for (i=1; i<numLines; i++) {
      if ((wLen->colValue->fArray[i] - wLen->colValue->fArray[i-1]) > width) {
        numIntervals++;
        start[numIntervals] = wLen->colValue->fArray[i] - width/2;
      }
      end[numIntervals] = wLen->colValue->fArray[i] + width/2;
    }
    numIntervals++;
    cpl_msg_debug(modName, "%d integration intervals found:", numIntervals);
    interval = wavIntervals = newDpoint(numIntervals);
    for (i=0; i<numIntervals; i++) {
      interval->x = start[i];
      interval->y = end[i];
      cpl_msg_debug(modName, "from %f to %f", interval->x, interval->y);
      interval = interval->next;
    }
    cpl_free(start);
    cpl_free(end);
  }
  return(wavIntervals);
}

/**
 * @memo 
 *   Delete a list of wavelength intervals.
 *
 * @return nothing
 * 
 * @param wavIntervals List of wavelength intervals.
 *   
 * @doc
 *   This is just a wrapper for deleteDpoint()
 *
 * @author C.Izzo
 */

void forgetWavIntervals(VimosDpoint *wavIntervals)
{
  deleteDpoint(wavIntervals);
}

/**
 * @memo 
 *   Compute the Match Index between a spectrum and a list
 *   of wavelength intervals.
 *
 * @return Match Index.
 * 
 * @param lambdaToCcd  Wavelength to CCD pixel conversion to be used.
 * @param wavIntervals Intervals of integration.
 * @param spectrum     Spectrum to be matched.
 * @param offset       Position of central wavelength in spectrum.
 *   
 * @doc
 *   The list of lambda intervals is converted into a list of
 *   CCD pixel intervals according to the given lambdaToCcd 
 *   transform. Next, the integral of the input spectrum is 
 *   computed on the given intervals: the value of this integral 
 *   is the Match Index. It is suggested to replace the spectrum
 *   with its own logarithm to avoid too strong a dependency of 
 *   the Match Index on the brightest lines (this is not done
 *   within this function, for obvious reasons of efficiency).
 *
 * @author C.Izzo
 */

double computeMatchIndex(VimosDistModel1D *lambdaToCcd, 
       VimosDpoint *wavIntervals, VimosFloatArray *spectrum, int offset) 
{
  VimosDpoint     *pixIntervals;
  VimosDpoint     *dWav;
  VimosDpoint     *dPix;
  double           matchIndex     = -1.;
  int              spectrumLength = spectrum->len;
  int              numWav         = 0;
  int              numPix         = 0;
  int              i, j, low, high;

 /*
  *  Conversion of wavelength intervals into CCD pixels intervals.
  */
  for (dWav=wavIntervals; dWav; dWav=dWav->next) 
    numWav++;

  dPix = pixIntervals = newDpoint(numWav);

  for (dWav=wavIntervals; dWav; dWav=dWav->next) {

    dPix->x = offset + computeDistModel1D(lambdaToCcd, dWav->x); /* Start */
    dPix->y = offset + computeDistModel1D(lambdaToCcd, dWav->y); /* End   */

    if (dPix->x < 0.) {     /* Interval start falls out of spectrum... */
      if (dPix->y > 0.) {   /* ... but interval end is in.             */
        dPix->x = 0.;
        dPix = dPix->next;
        numPix++;
      }
    }
    else if (dPix->x < spectrumLength) {  /* Interval start is in...   */
      if (dPix->y > spectrumLength) {     /* ... but interval end not. */
        dPix->y = spectrumLength;
      }
      dPix = dPix->next;
      numPix++;
    }
  }
 /*
  *  Compute match index
  */
  matchIndex = 0.;
  for (i=0, dPix=pixIntervals; i<numPix; i++, dPix=dPix->next) {
    low  = (int) (dPix->x + 0.5);    /* Nearest integer     */
    high = (int) (dPix->y + 1.5);    /* Nearest integer + 1 */
    for (j=low; j<high; j++) matchIndex += spectrum->data[j];
  }
  deleteDpoint(pixIntervals);
  return(matchIndex);
}

/**
 * @memo 
 *   Reduce intensity differences in a spectrum.
 *
 * @return New spectrum.
 * 
 * @param spectrum     Spectrum to process.
 *   
 * @doc
 *   Currently the "equalization" of intensities is performed
 *   by simply make the logarithm of the spectrum itself. Negative
 *   numbers are avoided.
 *
 * @author C.Izzo
 */

VimosFloatArray *equalizeSpectrum(VimosFloatArray *spectrum) 
{
  VimosFloatArray *logSpectrum;
  int              spectrumLength = spectrum->len;
  int              i;

  if ((logSpectrum = newFloatArray(spectrumLength))) {
   /*
    *  Logarithm of spectrum (avoiding negative numbers)
    */
    for (i=0; i<spectrumLength; i++)
      logSpectrum->data[i] = log10(MAX(spectrum->data[i], 1.));
  }
  return(logSpectrum);
}

/**
 * @memo 
 *   Search the maximum of the Match Index in terms of the
 *   lambdaToCcd transform coefficients.
 *
 * @return Lambda to CCD transform corresponding to the max Match
 *         Index found.
 * 
 * @param lambdaToCcd  First guess for the wavelength calibration
 * @param lineCat      Line catalog.
 * @param width        Slit width (in the same wavelength units used
 *                     in line catalog).
 * @param spectrum     Spectrum to be matched.
 * @param offset       Position of central wavelength in spectrum.
 *   
 * @doc
 *   The Match Index is computed for a grid of values of the lambda
 *   to CCD transform coefficients. Each coefficient will be varied 
 *   on a range inversely proportional to its corresponding power
 *   of a reference lambda in the polynomial transform (variation 
 *   on the CCD is depending linearly on the variations of the 
 *   coefficients, modulated by powers of lambda). The chosen value 
 *   for lambda is its highest possible value. The max Match Index 
 *   for all the values in the grid is found, and the corresponding 
 *   wavelength to CCD pixel transform is returned. The (local) minimum 
 *   is not searched (i.e., no fit is performed): the match should be
 *   just enough to ensure line identification, and proceed with the 
 *   usual peak search algorithm.
 *
 * @author C.Izzo
 */

VimosDistModel1D *findMaxMatchIndex(VimosDistModel1D *lambdaToCcd,
                  VimosTable *lineCat, float width, 
                  VimosFloatArray *spectrum, int offset)
{
  VimosFloatArray  *eqSpectrum;
  VimosDpoint      *wavIntervals;
  VimosDistModel1D *grid;
  VimosDistModel1D *start;
  VimosDistModel1D *best           = NULL;
  int               nCoeff;
  int              *count;
  double           *coeffRange;
  double           *coeffStep;
  int               nStep          = 50;
  int               pixRange       = 50;
  double            matchIndex     = 0.;
  double            bestMatchIndex = 0.;
  int               i              = 0;
  int               j              = 0;
  double            lambda;


  if (lambdaToCcd) {
    nCoeff = lambdaToCcd->order + 1;
    if ((eqSpectrum = equalizeSpectrum(spectrum))) {
      if ((wavIntervals = getWavIntervals(lineCat, width))) {
       /*
        *  Get a typical value of lambda in the polynomial, 
        *  necessary to evaluate the variation range to assign 
        *  to each coefficient. Since in the model lambda enters 
        *  as an offset from the central wavelength, half of the 
        *  spectral range might be a reasonable choice (currently 
        *  hardcoded to a reasonable value for FORS data).
        */
    
        lambda = 2000.;
    
       /*
        *  Get max variation for each polynomial coefficient.
        */
    
        if ((coeffRange = cpl_calloc(nCoeff, sizeof(double)))) {
          coeffRange[0] = pixRange;

          if ((coeffStep = cpl_calloc(nCoeff, sizeof(double)))) {
            coeffStep[0] = ((float) pixRange)/nStep;

            for (i = 1; i < nCoeff; i++) {
              coeffRange[i] = coeffRange[i-1]/lambda;
              coeffStep[i] = coeffRange[i]/nStep;
            }

            if ((grid = newDistModel1D(lambdaToCcd->order))) {
              for (i = 0; i < nCoeff; i++) {
                grid->coefs[i] = lambdaToCcd->coefs[i];
              }
              grid->offset = lambdaToCcd->offset;

              if ((start = newDistModel1D(lambdaToCcd->order))) {
                for (i = 0; i < nCoeff; i++) {
                  start->coefs[i] = grid->coefs[i] - coeffRange[i]/2.;
                }
                start->offset = grid->offset;
    
                if ((count = cpl_calloc(nCoeff, sizeof(int)))) {
                  for (i = 1; i < nCoeff; i++) count[i] = 1;
                 /*
                  *  Find best Match Index in a n-dimensional grid
                  *  of coefficents values.
                  */
                  if ((best = newDistModel1D(lambdaToCcd->order))) {
                    for (i = 0; i < nCoeff; i++) {
                      best->coefs[i] = lambdaToCcd->coefs[i];
                    }
                    best->offset = lambdaToCcd->offset;

                    grid->coefs[0] = start->coefs[0];
                    grid->coefs[1] = start->coefs[1];
                   /*
                    *  The navel of the function!
                    */
                    i = 0;
                    while (i < 2) {
                      if (count[i] < nStep) {
                        count[i]++;
                        grid->coefs[i] += coeffStep[i];
                        i = 0;
                        matchIndex = 
                        computeMatchIndex(grid, wavIntervals, 
                                          eqSpectrum, offset);
                        if (matchIndex > bestMatchIndex) {
                          bestMatchIndex = matchIndex;
                          for (j = 0; j < nCoeff; j++) {
                            best->coefs[j] = grid->coefs[j];
                          }
                        }
                      }
                      else {
                        count[i] = 1;
                        grid->coefs[i] = start->coefs[i];    /* reset */
                        i++;
                      }
                    }
                  }
                  cpl_free(count);
                }
                deleteDistModel1D(start);
              }
              deleteDistModel1D(grid);
            }
            cpl_free(coeffStep);
          }
          cpl_free(coeffRange);
        }
        forgetWavIntervals(wavIntervals);
      }
      deleteFloatArray(eqSpectrum);
    }
  }
  return(best);
}

/**
 * @memo 
 *   Find approximate position of the peak (if any) closest to the input
 *   histogram center.
 *
 * @return Approximate position of the closest-to-center peak. In case
 *         of failure, a negative number is returned.
 * 
 * @param data  Input histogram.
 * @param size  Number of bins in input histogram.
 *   
 * @doc
 *   This routine was added to increase robustness in the correct 
 *   identification of lines in an arc spectrum. The input is the 
 *   arc spectrum portion where a given line is expected: in case 
 *   more than one peak were included in the interval, the closest 
 *   one to the expected position (i.e., the center) is chosen. 
 *   Starting from the interval center, the spectrum is searched 
 *   for the closest-to-center flux excess, and then the rough 
 *   position of this excess is returned. This routine will allow 
 *   to use a wider window around the expected arc line position, 
 *   increasing the probability of including the searched line 
 *   - and limiting the effect of the equally increased probability 
 *   of including also other nearby lines. For this function to be 
 *   useful, therefore, it should be called after the first-guess 
 *   inverse dispersion relation was pre-tuned (using findMaxMatchIndex).
 *
 * @see findMaxMatchIndex()
 *
 * @author C.Izzo
 */

int findClosestPeak(float *data, int size) 
{
  float  max, min, level;
  float  percent = 1./4.;
  int    pos, posRight, posLeft, distRight, distLeft;
  int    i;

  if (data == NULL || size <= 10) return -1;
  
  pos = posRight = posLeft = size/2;
  distRight = distLeft = 0;

 /*
  *  Find min and max
  */
  min = max = data[0];
  for (i=1; i<size; i++) {
    if (data[i] > max) max = data[i];
    if (data[i] < min) min = data[i];
  }

 /*
  *  If the max equals the min we have a flat input, therefore
  *  no peak is found. Conventionally, return the position of the
  *  interval centre.
  */

  if (max-min < MIN_DIVISOR) return pos;

 /*
  *  Discrimination level: when pixels with values above this
  *  level are found, we are (probably) within a peak.
  */
  level = percent * max + (1 - percent) * min;

  if (data[pos] < level) {
    for (; (posRight < size) && (data[posRight] <= level); posRight++);
    distRight = posRight - size/2;
    for (; (posLeft  >= 0  ) && (data[posLeft] <= level); posLeft--);
    distLeft = size/2 - posLeft;

    if (distRight < distLeft) {
      posLeft = posRight;
      for (; (posRight < size) && (data[posRight] > level); posRight++);
    }
    else {
      posRight = posLeft;
      for (; (posLeft  >= 0  ) && (data[posLeft]  > level); posLeft--);
    }
  }
  else if (data[pos] > level) {
    for (; (posRight < size) && (data[posRight] > level); posRight++);
    for (; (posLeft  >= 0  ) && (data[posLeft]  > level); posLeft--);
  }

  pos = (posLeft + posRight) / 2;

  return(pos);

}


double **identPeaks(double *peak, int npeaks, double *line, int nlines,
                    double min_disp, double max_disp, double tolerance,
                    int *identified)
{

  int      i, j, k, l;
  int      nlint, npint;
  double   lratio, pratio;
  double   lo_start, lo_end, hi_start, hi_end, denom;
  double   disp, variation, prev_variation;
  int      max, maxpos, minl, mink;
  int      ambiguous;
  int      npeaks_lo, npeaks_hi;
  int     *peak_lo    = cpl_malloc(npeaks * sizeof(int));
  int     *peak_hi    = cpl_malloc(npeaks * sizeof(int));
  int    **ident      = cpl_malloc(npeaks * sizeof(int *));
  int     *nident     = cpl_calloc(npeaks, sizeof(int));
  int     *lident     = cpl_calloc(nlines, sizeof(int));

  double **output     = cpl_calloc(2, sizeof(double *));
  double  *xpos       = cpl_calloc(npeaks, sizeof(double));
  double  *lambda     = cpl_calloc(npeaks, sizeof(double));
  int     *ilambda    = cpl_calloc(npeaks, sizeof(int));
  double  *tmp_xpos    = cpl_calloc(npeaks, sizeof(double));
  double  *tmp_lambda  = cpl_calloc(npeaks, sizeof(double));
  double  *tmp_ilambda = cpl_calloc(npeaks, sizeof(double));
  int     *flag       = cpl_calloc(npeaks, sizeof(int));
  int      n = 0;
  int      nn;
  int      nseq       = 0;
  int      gap;
  int     *seq_length = cpl_calloc(npeaks, sizeof(int));
  int      found;


  *identified = 0;
  output[0] = xpos;
  output[1] = lambda;

  for (i = 0; i < npeaks; i++)
    ident[i] = cpl_malloc(3 * npeaks * sizeof(int));


  /*
   * This is just the number of intervals - one less than the number
   * of points (catalog wavelengths, or detected peaks). I do this
   * for performance reasons.
   */

  nlint = nlines;
  npint = npeaks;

  --nlint;
  --npint;


  /*
   * Here the big loops on catalog lines begins.
   */

  for (i = 1; i < nlint; i++) {


    /*
     * For each catalog line I take the previous and the next one, and 
     * compute the ratio of their distances from the central wavelengths.
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

      prev_variation = 1000;
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
           * is marked in the matrix ident[][], where the first index
           * correspond to a peak sequence number, and the second
           * index is the counter of the identifications made during
           * this whole process. The array of counters is nident[].
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

    if (nident[i] > 0) {


      /*
       * Reset the histogram of wavelengths assigned to the i-th peak.
       */

      for (j = 0; j < nlines; j++)
        lident[j] = 0;


      /*
       * Collect all the nident[i] wavelengths assigned to the i-th peak.
       */

      for (j = 0; j < nident[i]; j++)
        ++lident[ident[i][j]];


      /*
       * What wavelength was most frequently assigned to the i-th peak?
       */

      max = 0;
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
   * Check on identified peaks. Contaminations might be present and
   * should be excluded. Contamination from multiplexed spectra
   * consists of correctly identified lines that really belong to other
   * spectra. Generic light contamination and line misidentification
   * should have been almost all removed in the previous steps, but
   * it may still be present. The identified peaks are sorted into
   * separated self-consistent sequences. The longest of those sequences
   * is the one that is returned.
   */

  if (n > 1) {
    nn = 0;
    nseq = 0;
    for (k = 0; k < n; k++) {
      if (flag[k] == 0) {
        flag[k] = 1;
        xpos[nn] = tmp_xpos[k];
        lambda[nn] = tmp_lambda[k];
        ilambda[nn] = tmp_ilambda[k];
        ++seq_length[nseq];
        ++nn;
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
        ++nseq;
        k = 0;
      }
    }

    max = 0;
    for (i = 0; i < nseq; i++) {
      if (seq_length[i] > max) {
        max = seq_length[i];
        maxpos = i;
      }
    }

    nn = 0;
    for (i = 0; i < maxpos; i++)
      nn += seq_length[i];
  
    n = max;
    for (i = 0; i < n; i++, nn++) {
      xpos[i] = xpos[nn];
      lambda[i] = lambda[nn];
      ilambda[i] = ilambda[nn];
    }


    /* 
     * Are some wavelengths missing?
     */

    for (i = 1; i < n; i++) {
      gap = ilambda[i] - ilambda[i-1];
      if (gap > 1) {
        j = 1;
        disp = (lambda[i] - lambda[i-1]) / (xpos[i] - xpos[i-1]);
        hi_start = xpos[i-1] + (line[ilambda[i-1] + j] - lambda[i-1]) / disp;
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
      }
    }

    /*
     * Try to extrapolate forward
     */

    found = 1;
    while (ilambda[n-1] < nlines - 1 && found) {
      hi_start = xpos[n-1] + (line[ilambda[n-1]+1] - lambda[n-1]) / max_disp;
      hi_end   = xpos[n-1] + (line[ilambda[n-1]+1] - lambda[n-1]) / min_disp;
      found = 0;
      for (i = 0; i < npeaks; i++) {
        if (peak[i] > hi_start && peak[i] < hi_end) {
          xpos[n] = peak[i];
          lambda[n] = line[ilambda[n-1]+1];
          ilambda[n] = ilambda[n-1]+1;
          ++n;
          found = 1;
          break;
        }
      }
    }


    /*
     * Try to extrapolate backward
     */

    found = 1;
    while (ilambda[0] > 0 && found) {
      hi_start = xpos[0] + (line[ilambda[0]-1] - lambda[0]) / min_disp;
      hi_end   = xpos[0] + (line[ilambda[0]-1] - lambda[0]) / max_disp;
      found = 0;
      for (i = 0; i < npeaks; i++) {
        if (peak[i] > hi_start && peak[i] < hi_end) {
          for (j = n; j > 0; j--) {
            xpos[j] = xpos[j-1];
            lambda[j] = lambda[j-1];
            ilambda[j] = ilambda[j-1];
          }
          xpos[0] = peak[i];
          lambda[0] = line[ilambda[0]-1];
          ilambda[0] = ilambda[0]-1;
          ++n;
          found = 1;
          break;
        }
      } 
    }
  }


  /*
   * At this point all peaks are processed. Free memory, and return 
   * the result.
   */

/************************************************
  for (i = 0; i < npeaks; i++) {
    printf("Peak %d:\n   ", i);
    for (j = 0; j < nident[i]; j++)
      printf("%.2f, ", line[ident[i][j]]);
    printf("\n");
  }

  printf("\n");

  for (i = 0; i < n; i++)
    printf("%.2f, %.2f\n", xpos[i], lambda[i]);
**************************************************/

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

  *identified = n;

  if (n == 0) {
    cpl_free(output[0]);
    cpl_free(output[1]);
    cpl_free(output);
    return NULL;
  }

  return output;

}

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

double *collectPeaks(float *row, int length, 
                     float level, float exp_width, int *npeak)
{
  int     i, j;
  int     nint   = length - 1;
  int     n      = 0;
  int     width  = 2 * ceil(exp_width / 2) + 1;
  int     start  = width / 2;
  int     end    = length - width / 2;
  int     step;
  int     minbox = 21;
  float   min;
  float  *smo;
  float  *flat;
  double *peak   = cpl_calloc(length/2, sizeof(double));


  /*
   * If lines have a flat top (as in the case of broad slit), smooth
   * before determining the max.
   */

  if (width > 3) {
    smo = cpl_calloc(length, sizeof(float));
    start = width / 2;
    end = length - width / 2;
    for (i = 0; i < start; i++)
      smo[i] = row[i];
    for (i = start; i < end; i++) {
      for (j = i - start; j <= i + start; j++)
        smo[i] += row[j];
      smo[i] /= width;
    }
    for (i = end; i < length; i++)
      smo[i] = row[i];
  }
  else {
    smo = row;
  }


  /*
   * Pass a min filter to subtract the continuum
   */

  flat = cpl_calloc(length, sizeof(float));
  start = minbox / 2;
  end = length - minbox / 2;

  for (i = start; i < end; i++) {
    min = smo[i-start];
    for (j = i - start + 1; j <= i + start; j++)
      if (min > smo[j])
        min = smo[j];
    flat[i] = min;
  }

  if (width > 3) {
    cpl_free(smo);
  }

  for (i = 0; i < start; i++)
    flat[i] = row[i] - flat[start];
  for (i = start; i < end; i++)
    flat[i] = row[i] - flat[i];
  for (i = end; i < length; i++)
    flat[i] = row[i] - flat[end - 1];

  if (width > 20) 
    step = width / 2;
  else 
    step = 1;

  /*
   * Collect all relative maxima along row, that are higher than the
   * specified level.
   */

  for (i = step; i < nint - step + 1; i += step) {
    if (flat[i] > level) {
      if (flat[i] >= flat[i-step] && flat[i] > flat[i+step]) {
        if (flat[i-step] != 0.0 && flat[i+step] != 0.0) {
          peak[n] = i + step * values_to_dx(flat[i-step], flat[i], flat[i+step]);
          ++n;
        }
      }
    }
  }

  *npeak = n;

  cpl_free(flat);

  if (n == 0) {
    cpl_free(peak);
    return NULL;
  }

  return peak;

}


double *collectPeaks_double(double *row, int length, 
                     float level, float exp_width, int *npeak)
{
  int     i, j;
  int     nint   = length - 1;
  int     n      = 0;
  int     width  = 2 * ceil(exp_width / 2) + 1;
  int     start  = width / 2;
  int     end    = length - width / 2;
  int     step;
  int     minbox = 21;
  float   min;
  double *smo;
  float  *flat;
  double *peak   = cpl_calloc(length/2, sizeof(double));


  /*
   * If lines have a flat top (as in the case of broad slit), smooth
   * before determining the max.
   */

  if (width > 3) {
    smo = cpl_calloc(length, sizeof(float));
    start = width / 2;
    end = length - width / 2;
    for (i = 0; i < start; i++)
      smo[i] = row[i];
    for (i = start; i < end; i++) {
      for (j = i - start; j <= i + start; j++)
        smo[i] += row[j];
      smo[i] /= width;
    }
    for (i = end; i < length; i++)
      smo[i] = row[i];
  }
  else {
    smo = row;
  }


  /*
   * Pass a min filter to subtract the continuum
   */

  flat = cpl_calloc(length, sizeof(float));
  start = minbox / 2;
  end = length - minbox / 2;

  for (i = start; i < end; i++) {
    min = smo[i-start];
    for (j = i - start + 1; j <= i + start; j++)
      if (min > smo[j])
        min = smo[j];
    flat[i] = min;
  }

  if (width > 3) {
    cpl_free(smo);
  }

  for (i = 0; i < start; i++)
    flat[i] = row[i] - flat[start];
  for (i = start; i < end; i++)
    flat[i] = row[i] - flat[i];
  for (i = end; i < length; i++)
    flat[i] = row[i] - flat[end - 1];

  if (width > 20) 
    step = width / 2;
  else 
    step = 1;

  /*
   * Collect all relative maxima along row, that are higher than the
   * specified level.
   */

  for (i = step; i < nint - step + 1; i += step) {
    if (flat[i] > level) {
      if (flat[i] >= flat[i-step] && flat[i] > flat[i+step]) {
        if (flat[i-step] != 0.0 && flat[i+step] != 0.0) {
          peak[n] = i + step * values_to_dx(flat[i-step], flat[i], flat[i+step]);
          ++n;
        }
      }
    }
  }

  *npeak = n;

  cpl_free(flat);

  if (n == 0) {
    cpl_free(peak);
    return NULL;
  }

  return peak;

}


/**@}*/
