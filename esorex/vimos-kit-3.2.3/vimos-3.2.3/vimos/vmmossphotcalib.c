/* $Id: vmmossphotcalib.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <string.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmimage.h"
#include "vmmatrix.h"
#include "vmtable.h"
#include "vmgrismtable.h"
#include "vmextractiontable.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmdistmodels.h"
#include "vmmossphotcalib.h"
#include "vmimgutils.h"
#include "cpl.h"


static int
polint2(float xa[], float ya[], int nDat, float x, float *y, float *dy)
{
  char modName[] = "polint2";
  int i,m,ns=0;
  float den,dif,dift,ho,hp,w;
  float *c,*d;

  dif=fabs(x-xa[0]);
  c = (float *) cpl_calloc(nDat, sizeof(float));
  d = (float *) cpl_calloc(nDat, sizeof(float));
  for (i=0; i<nDat; i++) 
  {
    if ( (dift=fabs(x-xa[i])) < dif) 
    {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=0; m<nDat-1; m++) 
  {
    for (i=0; i<nDat-m-1; i++) 
    {
      ho=xa[i]-x;
      hp=xa[i+m+1]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0)
      {
	cpl_msg_error(modName,"Error: two identical x values");
	return EXIT_FAILURE;
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (nDat-m) ? c[ns+1] : d[ns--]));
  }
  cpl_free(d);
  cpl_free(c);
  return EXIT_SUCCESS;
}


static double *
readCalSphotModel(VimosDescriptor *desc, int modelDeg)
{
  const char modName[] = "readCalSphotModel";
  char       comment[80];
  double     dValue;
  double    *coeffs;
  int        i;

  coeffs = cpl_malloc((modelDeg+1)*sizeof(double)) ;
  for (i = 0; i <= modelDeg; i++)
  {
    if (readDoubleDescriptor(desc, pilTrnGetKeyword("SphotModel", i), 
			     &dValue, comment) == VM_FALSE) {
      cpl_free(coeffs);
      coeffs = NULL;
      cpl_msg_error(modName, "Cannot read descriptor %s", 
		  pilTrnGetKeyword("SphotModel", i)); 
    }
    else {
      coeffs[i] = dValue;
    }
  }    
  return coeffs;
}


static VimosBool
writeCalSphotModel(double *modelCoeffs, int modelDeg, VimosDescriptor *desc)
{
  const char modName[] = "writeCalSphotModel";
  int        i;

  if (writeIntDescriptor(&desc, pilTrnGetKeyword("SphotOrder"), 
                         modelDeg, "") == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found", 
                pilTrnGetKeyword("SphotOrder"));
    return VM_FALSE;
  }

  for (i = 0; i <= modelDeg; i++)
  {
    if (writeDoubleDescriptor(&desc, pilTrnGetKeyword("SphotModel", i), 
			      modelCoeffs[i], "") == VM_FALSE) {
      cpl_msg_error(modName, "Descriptor %s not found", 
		  pilTrnGetKeyword("SphotModel", i));
      return VM_FALSE;
    }
  }
  return VM_TRUE;
}

static VimosFloatArray *
resample1D(VimosFloatArray *input, float startLambda, float delLambda,
           float firstLambda, float lastLambda, float newDelLambda)
{
  char     modName[] = "resample1D";
  int      i, j;
  int      nPixIn;
  int      nPixOut;
  int      pix;
  float    lambda;
  float    flux;
  float    outFlux;
  float    sigma;
  
  VimosFloatArray *outSpectrum, *inputLambda;

  nPixIn = input->len;

  if (firstLambda < startLambda)
  {
    cpl_msg_error(modName, 
		"Lower output wavelength smaller than lower input one");
    return NULL;
  }
  if (lastLambda > (startLambda + (float) nPixIn * delLambda))
  {
    cpl_msg_error(modName, 
		"Upper output wavelength larger than upper input one");
    return NULL;
  }
    
  nPixOut = (int) ((lastLambda - firstLambda) / newDelLambda);
  outSpectrum = newFloatArray(nPixOut);
  inputLambda = newFloatArray(nPixIn);
  
  for (i=0; i<nPixIn; i++)
  {  
    inputLambda->data[i] = startLambda + delLambda * (float) i;
    input->data[i] /= delLambda;
  }  

  /* interpolation of the input spectrum, on a wavelength scale one/fourth 
     of the required output one */
  for (i=0; i<nPixOut; i++)
  {
    lambda = firstLambda + newDelLambda * (float) i;
    pix = (int) ((lambda - startLambda) / delLambda);
    outFlux = 0.;

    if (pix < 2)
    {
      /* special case for the first two pixels */
      for (j=0; j<4; j++)
      {
	polint2(inputLambda->data, input->data, 3, lambda, &flux, &sigma);
	outFlux += flux * newDelLambda / 4.;
	lambda += (float) j * newDelLambda / 4.;
      }
    }
    else if (pix > (nPixIn - 2))
    {
      /* special case for the last two pixels */
      for (j=0; j<4; j++)
      {
	polint2(&inputLambda->data[nPixIn-4], &input->data[nPixIn-4], 3, 
	       lambda, &flux, &sigma);
	outFlux += flux * newDelLambda / 4.;
	lambda += (float) j * newDelLambda / 4.;
      }
    }
    else
    {
      /* general case */
      for (j=0; j<4; j++)
      {
	polint2(&inputLambda->data[pix-2], &input->data[pix-2], 4, lambda, 
	       &flux, &sigma);
	outFlux += flux * newDelLambda / 4.;
	lambda += (float) j * newDelLambda / 4.;
      }
    }
    outSpectrum->data[i] = outFlux;
  }

  return outSpectrum;
}


int
VmSpCalPhot(VimosImage *inputImage, VimosTable *sphotStdTable, int stdRow,
            int sensFuncDeg)
{
  char     modName[] = "VmSpCalPhot";
  int      constDLambda=1;
  int      i,j;
  int      nStdPoints;
  int      delLambdaNew;
  int      nConstDLambda;
  int      stdDLambda;
  int      inpSpecLen;
  int      firstFlux;
  int      lastFlux;
  int      nInt;
  int      nBins;
  int      pix;
  float    specDLambda;
  float    specStartLambda;
  float    firstLambda;
  float    lastLambda;
  float    newDLambda;
  float    lambda;
  float    magStd;
  float    expTime;
  float    ratio;
  float    specCounts;
  float    stdFlux;
  double   dValue;
  double   fitRms;
  double  *sensFuncCoeffs = NULL;
  char     comment[80];

  float    fnuZero = 3.68e-20;  /*absolute flux zero point, from IRAF */
  float    vLight = 2.9979246e18; /* speed of light in Angstrom/sec */
  
  VimosColumn  *lambdaCol, *magCol, *delLambdaCol;
  VimosFloatArray *inputSpectrum, *resampSpectrum;
  VimosDpoint *sensFunc = NULL;

  /* read in the input spectrum */
  readDoubleDescriptor(inputImage->descs, pilTrnGetKeyword("ExposureTime"), 
		       &dValue, comment);
  expTime = (float) dValue;
  readDoubleDescriptor(inputImage->descs, pilTrnGetKeyword("Cdelt",1), 
		       &dValue, comment);
  specDLambda = (float) dValue;
  readDoubleDescriptor(inputImage->descs, pilTrnGetKeyword("Crval",1), 
		       &dValue, comment);
  specStartLambda = (float) dValue;
  readIntDescriptor(inputImage->descs,pilTrnGetKeyword("Naxis",1),&inpSpecLen,
		    comment);

  inputSpectrum = newFloatArray(inpSpecLen);

  for (i=0; i<inpSpecLen; i++)
    inputSpectrum->data[i] = inputImage->data[i+stdRow*inpSpecLen];

  /* read in the spectrophotometric standard table columns */
  lambdaCol = findColInTab(sphotStdTable, "LAMBDA");
  magCol = findColInTab(sphotStdTable, "MAG");
  delLambdaCol = findColInTab(sphotStdTable, "DELTA_LAMBDA");
  nStdPoints = delLambdaCol->len;

  /* find out if the wavelength sampling of the spectrophotometric standard 
     is constant, or if we must work over fractions of the table separately */
  for (i=1; i<nStdPoints; i++)
  {
    stdDLambda = delLambdaCol->colValue->fArray[0];
    if (delLambdaCol->colValue->fArray[i]!=stdDLambda)
    {
      constDLambda=0;
      delLambdaNew=i;
      nConstDLambda=i;
    }
  }

  /* constant sampling, so we can work with the whole spectrum at once */
  if (constDLambda)
  {
    /* determine resampling factor for the input spectrum */
    ratio = (float) stdDLambda / specDLambda;
    if (ratio >= 2.)
    {
      newDLambda = stdDLambda / (int) ratio;
      nInt = (int) ratio;
    }
    else if (ratio < 2. && ratio >= 1.5)
    {    
      newDLambda = stdDLambda / 2.;
      nInt = 2;
    }
    else if (ratio < 1.5 && ratio >= 0.8)
    {    
      newDLambda = (float) stdDLambda;
      nInt = 1;
    }
    else
    {
      /* if the sampling of the standard star if at a higher resolution than
	 our spectrum, we cannot use the table */
      cpl_msg_error(modName, "Mismatch between standard Table and spectrum");
      return EXIT_FAILURE;
    }
  
    /* find first pixel with flux > 0 in spectrum, and first wavelength 
       matching the table ones */
    i=0;
    while (inputSpectrum->data[i] < 1)
      i++;
    firstFlux = i;
    firstLambda = specStartLambda + specDLambda * (float) firstFlux;
    i=0;
    while (lambdaCol->colValue->fArray[i] < firstLambda)
      i++;
    firstLambda = lambdaCol->colValue->fArray[i];
    
    /* find last pixel with flux > 0 in spectrum, and last wavelength 
       matching the table ones */
    i=inpSpecLen-1;
    while (inputSpectrum->data[i] < 1)
      i--;
    lastFlux = i;
    lastLambda = specStartLambda + specDLambda * (float) lastFlux;
    i=0;
    while (lambdaCol->colValue->fArray[i] < lastLambda)
      i++;
    lastLambda = lambdaCol->colValue->fArray[i-1];
    
    /* resample the input spectrum */
    resampSpectrum = resample1D(inputSpectrum,specStartLambda,specDLambda,
				firstLambda,lastLambda,newDLambda);

    /* integrate the flux over the same table wavelength bins, compute
       the calibration flux, and the aperture sensitivity function */
    nBins = (int)((lastLambda - firstLambda) / stdDLambda);    
    sensFunc = newDpoint(nBins);

    for (i=0; i<nBins; i++)
    {
      specCounts = 0.;
      lambda = firstLambda + (float) (stdDLambda * i);
      for (j=0; j<nInt; j++)
      {
	pix = (int) ((lambda - firstLambda) / newDLambda);
	specCounts += resampSpectrum->data[pix+j];
      }

      /* normalize the counts to the exposure time and passband width */
      specCounts /= (expTime * stdDLambda);

      /* compute the calibration flux from the tabulated mag */
      pix = (int) ((lambda - lambdaCol->colValue->fArray[0]) / stdDLambda);
      magStd = magCol->colValue->fArray[pix];
      stdFlux = fnuZero * pow(10., -0.4*magStd) * vLight / SQR(lambda);
      /* derive the sensitivity function */
      sensFunc[i].x = (double) lambda;
      sensFunc[i].y = (double) (2.5 * log10(specCounts / stdFlux));
    }
  }
  else
  {
    /* we must resample the spectrum in chunks of constant stdDLambda */
  }

  /* fit the sensitivity function with a polynomial, and write it
     in the spectrophotometric table */
  sensFuncCoeffs = fit1DPoly(sensFuncDeg, sensFunc, nBins, &fitRms);
  writeCalSphotModel(sensFuncCoeffs, sensFuncDeg, sphotStdTable->descs);

  return EXIT_SUCCESS;
}


VimosImage *
VmSpApplyPhotOld(VimosImage *inputImage, VimosTable *sphotStdTable)
{
  char     comment[80];
  int      imageXlen, imageYlen;
  int      sensFuncDeg;
  int      index;
  int      i, j;
  float    specDLambda;
  float    specStartLambda;
  float    lambda;
  float    expTime;
  float    normCounts;
  float    normFlux;
  float    sensValue;
  double   dValue;
  double   sensFunc;
  double  *sensFuncCoeffs = NULL;

  VimosImage  *outputImage;
  imageXlen = inputImage->xlen;
  imageYlen = inputImage->ylen;

  /* read in spectral parameters */
  readDoubleDescriptor(inputImage->descs, pilTrnGetKeyword("ExposureTime"), 
		       &dValue, comment);
  expTime = (float) dValue;
  readDoubleDescriptor(inputImage->descs, pilTrnGetKeyword("Cdelt",1), 
		       &dValue, comment);
  specDLambda = (float) dValue;
  readDoubleDescriptor(inputImage->descs, pilTrnGetKeyword("Crval",1), 
		       &dValue, comment);
  specStartLambda = (float) dValue;

  /* create output image */
  outputImage = newImageAndAlloc(imageXlen, imageYlen);
  copyAllDescriptors(inputImage->descs, &(outputImage)->descs);

  /* read in the order of the polynomial fit to the sensitivity function */
  readIntDescriptor(sphotStdTable->descs, pilTrnGetKeyword("SphotOrder"),
		    &sensFuncDeg, comment);
  /* read in the coefficients of the polynomial */
  sensFuncCoeffs = readCalSphotModel(sphotStdTable->descs, sensFuncDeg);

  /* loop over lambda (i.e. the image X axis) */
  for (i = 0; i < imageXlen; i++)
  {
    lambda = specStartLambda + (float) (specDLambda * i);

    /* re-construct the value of the sensitivity function */
    sensFunc = sensFuncCoeffs[0];
    for (j = 1; j <= sensFuncDeg; j++)
      sensFunc += sensFuncCoeffs[j] * ipow(lambda, j);
    sensFunc /= 2.5;
    sensValue = pow(10., (float) sensFunc);

    /* loop over the variuos spectra (i.e. the image Y axis) */
    for (j = 0; j < imageYlen; j++)
    {
      index = i + j * imageXlen;
      /* normalize the counts */
      normCounts = inputImage->data[index] / (expTime * specDLambda);
      /* compute the flux */
      normFlux = normCounts / sensValue;
      /* write it to the output image */
      outputImage->data[index] = normFlux;
    }
  }

  return outputImage;
}


VimosImage *
VmSpApplyPhot(VimosImage *image, VimosTable *sphotTable, VimosTable *atmTable)
{
  const char modName[] = "VmSpApplyPhot";

  double     *data;
  int         tlength;
  int         xlength = image->xlen;
  int         ylength = image->ylen;
  int         index   = 0;
  int         i, j;
  double      step, start, gain, time, airmass;
  VimosImage *extinction  = NULL;
  VimosImage *outputImage = NULL;


  if (!atmTable && !sphotTable) {
    cpl_msg_error(modName, 
                  "Missing both atmospheric extinction and response curves!");
    return NULL;
  }

  if (sphotTable) {
    data = tblGetDoubleData(sphotTable, "RESPONSE");

    if (data == NULL) {
      cpl_msg_error(modName, "Missing RESPONSE column in input table");
      return NULL;
    }

    tlength = tblGetSize(sphotTable, "RESPONSE");

    if (xlength != tlength) {
      cpl_msg_error(modName, "Input table and input image are incompatible");
      return NULL;
    }
  }

  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("Cdelt", 1), 
                       &step, NULL) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found", 
                  pilTrnGetKeyword("Cdelt", 1));
    return NULL;
  }

  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("Crval", 1), 
                       &start, NULL) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found", 
                  pilTrnGetKeyword("Crval", 1));
    return NULL;
  }

/***
  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("SummedExposureTime"),
                       &time, NULL) == VM_FALSE) {
***/
    if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("ExposureTime"),
                         &time, NULL) == VM_FALSE) {
      cpl_msg_error(modName, "Neither descriptor %s nor descriptor "
                    "%s were found", pilTrnGetKeyword("SummedExposureTime"),
                    pilTrnGetKeyword("ExposureTime"));
      return NULL;
    }
/***
  }
***/

  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("Adu2Electron", 1),
                       &gain, NULL) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found", 
                  pilTrnGetKeyword("Adu2Electron", 1));
    return NULL;
  }

  if (atmTable) {

    /*
     * Airmass priorities... The airmass is read from ESO PRO AIRMASS,
     * to give precedence to the software airmass of a product that
     * was already airmass corrected. If ESO PRO AIRMASS is missing,
     * the airmass is computed from Ra, Dec, and time. If even this
     * fails, one last attempt to get the airmass from the standard
     * FITS keyword AIRMASS is made.
     */

    if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("AirMass"),
                             &airmass, NULL) == VM_FALSE) {
      if (VmComputeAirmass(image, &airmass) == EXIT_FAILURE) {
        if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("Airmass"),
                             &airmass, NULL) == VM_FALSE) {
          cpl_msg_error(modName, "Descriptor %s not found", 
                        pilTrnGetKeyword("Airmass"));
          return NULL;
        }
      }
    }

    cpl_msg_info(modName, "Mean airmass: %.4f", airmass);

    /*
     * Map the atmospheric extinction factors to the same lambda sampling
     * of the extracted spectrum, and convert to actual flux loss.
     */

    extinction = newImageAndAlloc(xlength, 1);
    mapTable(extinction, start, step, atmTable, "WAVE", "EXTINCTION");
    constArithLocal(extinction, 0.4 * airmass, VM_OPER_MUL);

    for (i = 0; i < xlength; i++)
      if (extinction->data[i] > 0.0)
        extinction->data[i] = pow(10., extinction->data[i]);

  }

  /*
   * Create an output image from the reduced spectra, correcting
   * them to airmass = 0, and for the instrument response.
   */

  outputImage = newImageAndAlloc(xlength, ylength);
  copyAllDescriptors(image->descs, &(outputImage)->descs);

  if (atmTable && sphotTable) {
    for (i = 0; i < ylength; i++) {
      for (j = 0; j < xlength; j++) {
        outputImage->data[index] = image->data[index] 
                                 * extinction->data[j] 
                                 * data[j];
        index++;
      }
    }
  }
  else if (atmTable) {
    for (i = 0; i < ylength; i++) {
      for (j = 0; j < xlength; j++) {
        outputImage->data[index] = image->data[index] * extinction->data[j];
        index++;
      }
    }
  }
  else if (sphotTable) {
    for (i = 0; i < ylength; i++) {
      for (j = 0; j < xlength; j++) {
        outputImage->data[index] = image->data[index] * data[j];
        index++;
      }
    }
  }

  deleteImage(extinction);

  if (sphotTable)
    constArithLocal(outputImage, gain / time / step, VM_OPER_MUL);

  return outputImage;

}
