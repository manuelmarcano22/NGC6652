/* $Id: vmimgutils.c,v 1.3 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Revision: 1.3 $
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
#include <pilastroutils.h>
#include <piltranslator.h>
#include <pildate.h>
#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmdetector.h"
#include "vmimgutils.h"
#include "vmutils.h"
#include "cpl.h"


#define COMMENT_LENGTH    (80)


/**
 * @name vmimgutils Image processing utilities
 *
 * The module collects utility functions for generic image operations.
 */

/**@{*/


/**
 * @memo 
 *   Trim away overscans from an image.
 *
 * @return EXIT_SUCCESS/EXIT_FAILURE
 * 
 * @param image   Image with overscan regions
 *   
 * @doc
 *   This function is used trim away the overscans regions from an image.
 *
 * @author C. Izzo, P. Sartoretti
 */

int
trimOverscans(VimosImage *image) 
{
  char            modName[] = "trimOverscans";
  float          *trimdata;
  int             startPortsX, startPortsY, endPortsX, endPortsY;
  int             sizePortsX, sizePortsY, nports;
  VimosPort      *ports, *currport;
  double          crpix1, crpix2;

  cpl_msg_debug(modName, "Trimming Overscans");

  if ((ports = getPorts(image, &nports)) == NULL) {
    cpl_msg_error(modName,"Cannot read Pre/OverScans from input image");
    return EXIT_FAILURE;
  }
  
  for (currport = ports; currport; currport = currport->next) {

    if (currport->prScan->nX + currport->ovScan->nX) {
      break;
    }

   /*
    *  No overscans, no need to trim.
    */

    cpl_msg_debug(modName, "No overscans, no need to trim.");
    deletePortList(ports);
    return EXIT_SUCCESS;
  }

 /*
  * Get the total area covered by the entire readout window:
  */

  getTotalReadoutWindow(ports, 
      &startPortsX, &startPortsY, &sizePortsX, &sizePortsY);

  deletePortList(ports);

  endPortsX = sizePortsX + startPortsX;
  endPortsY = sizePortsY + startPortsY;
  
  cpl_msg_debug(modName,
            "Extract image excluding overscans: start=(%d,%d) end=(%d,%d)",
            startPortsX, startPortsY, endPortsX, endPortsY);

 /*
  * Extract image excluding overscans 
  */

  trimdata = extractFloatImage(image->data, image->xlen, image->ylen,
             startPortsX, startPortsY, sizePortsX, sizePortsY);
  
  cpl_free(image->data);
  image->data = trimdata;
  image->xlen = sizePortsX;
  image->ylen = sizePortsY;
  

 /*
  * Update header
  */

  writeIntDescriptor(&(image->descs), pilTrnGetKeyword("Naxis", 1),
                     sizePortsX, pilTrnGetComment("Naxis"));
  writeIntDescriptor(&(image->descs), pilTrnGetKeyword("Naxis", 2),
                     sizePortsY, pilTrnGetComment("Naxis"));
  
  if (VM_TRUE == readDoubleDescriptor(image->descs,
                                      pilTrnGetKeyword("Crpix", 1),
                                      &crpix1, NULL)) { 
    writeDoubleDescriptor(&(image->descs), pilTrnGetKeyword("Crpix", 1),
                          (double) (crpix1-startPortsX),
                          pilTrnGetComment("Crpix"));
  }
  else {
    cpl_msg_error(modName, "Cannot read descriptor %s",
                pilTrnGetKeyword("Crpix", 1));
    return EXIT_FAILURE;
  }
    
  if (VM_TRUE == readDoubleDescriptor(image->descs,
                                      pilTrnGetKeyword("Crpix", 2),
                                      &crpix2, NULL)) { 
    writeDoubleDescriptor(&(image->descs), pilTrnGetKeyword("Crpix", 2),
                          (double) (crpix2-startPortsY),
                          pilTrnGetComment("Crpix"));
  }
  else {
    cpl_msg_error(modName, "Cannot read descriptor %s",
                pilTrnGetKeyword("Crpix", 2));
    return EXIT_FAILURE;
  }
                      
  if ((deleteSetOfDescriptors(&(image->descs), "*OVSC*")) == 0) {
      cpl_msg_warning(modName, "Cannot delete overscan keywords: not found");
    }
  if ((deleteSetOfDescriptors(&(image->descs), "*PRSC*")) == 0) {
    cpl_msg_warning(modName, "Cannot delete prescan keywords: not found");
  } 
  
  return EXIT_SUCCESS;
}
    

/**
 * @memo 
 *   Create overscans on a master bias image.
 *
 * @return Master bias with overscans, or NULL in case of failure.
 * 
 * @param mbias    Master bias without overscans
 * @param refImage Reference image
 *   
 * @doc
 *   Allocate and return a new master bias frame idetincal
 *   to the input one, but with overscans added. The overscan
 *   regions are inserted at the positions where they appear
 *   on the input reference image.
 *
 * @author P.Sartoretti and C.Izzo
 */

VimosImage *
growOverscans(VimosImage *mbias, VimosImage *refImage) 
{
  char         modName[] = "growOverscans";
  VimosImage  *ima_out = NULL;
  int          startPortsX, startPortsY;
  int          sizePortsX, sizePortsY, nports;
  VimosPort   *ports;
  int          deltaX, deltaY, sizeOveX, sizeOveY;
  float       *preScanX, *preScanY, *oveScanX, *oveScanY;
  
  cpl_msg_debug(modName, "Regrow overscans");

  if ((mbias == NULL) || (refImage == NULL)) {
    cpl_msg_error(modName, "Null input images");
    return NULL;
  }

  deltaX  = refImage->xlen - mbias->xlen;
  deltaY  = refImage->ylen - mbias->ylen;

  if ((deltaX == 0) && (deltaY == 0)) {
    cpl_msg_debug(modName, "Input and Reference image have the same"
                  " dimensions: return input master bias");
    return mbias;
  }
  else {
    if (deltaX && deltaY) {
      cpl_msg_error(modName, "Cannot grow overscans both in X and Y directions");
      return NULL;
    }
    ima_out = newImageAndAlloc(refImage->xlen, refImage->ylen);
    if ((ports = getPorts(refImage, &nports)) == NULL) {
      cpl_msg_error(modName,
        "Cannot read overscan keywords from reference image");
      return NULL;
    }

   /*
    * Get the total area covered by total readout window:
    */

    getTotalReadoutWindow(ports,
        &startPortsX, &startPortsY, &sizePortsX, &sizePortsY);

    if (VM_FALSE == (insertFloatImage(ima_out->data, ima_out->xlen,
                     ima_out->ylen, startPortsX, startPortsY, sizePortsX, 
                     sizePortsY, mbias->data))) {
      cpl_msg_error(modName, "Failure in insertFloatImage");
      return NULL;
    }

    if (deltaX) {

     /*
      * This means scan direction along X, then it should be 
      * startPortY=0, sizePortY=0, and deltaY=0
      */

      if (startPortsX) {
        preScanX = extractFloatImage(mbias->data, mbias->xlen, 
                   mbias->ylen, 0, 0, startPortsX, mbias->ylen);
        if ((insertFloatImage(ima_out->data, ima_out->xlen, ima_out->ylen,
                   0, 0, startPortsX, ima_out->ylen, preScanX)) == VM_FALSE) {
          cpl_msg_error(modName, "Cannot create preScan");
          return NULL;
        }
        cpl_free(preScanX);
      }

      sizeOveX = ima_out->xlen - (sizePortsX + startPortsX);
        
      if (sizeOveX) {
        oveScanX = extractFloatImage(mbias->data, mbias->xlen, mbias->ylen, 
                   mbias->xlen-sizeOveX, 0, sizeOveX, mbias->ylen);
        if ((insertFloatImage(ima_out->data, ima_out->xlen, ima_out->ylen,
                              ima_out->xlen-sizeOveX, 0, sizeOveX,
                              ima_out->ylen, oveScanX)) == VM_FALSE) {
          cpl_msg_error(modName, 
            "Cannot create overScan : Failure in insertFloatImage");
          return NULL;
        }
        cpl_free(oveScanX);
      }
    }
    else {

     /*
      * This is the case of readOut in Y, therefore startPortX=0 and
      * sizePortX=0, and deltaX=0
      */

      if (startPortsY) {
        preScanY = extractFloatImage(mbias->data, mbias->xlen, 
                   mbias->ylen, 0, 0, mbias->xlen, startPortsY);
        if ((insertFloatImage(ima_out->data, ima_out->xlen, ima_out->ylen,
                   0, 0, ima_out->xlen, startPortsY, preScanY)) == VM_FALSE) {
          cpl_msg_error(modName, "Cannot create preScan : Failure in "
                      "insertFloatImage");
          return NULL;
        }
        cpl_free(preScanY);
      }

      sizeOveY = ima_out->ylen - (startPortsY+sizePortsY);

      if (sizeOveY) {
        oveScanY = extractFloatImage(mbias->data, mbias->xlen, mbias->ylen, 
                   0, mbias->ylen-sizeOveY, mbias->xlen, sizeOveY);
        if ((insertFloatImage(ima_out->data, ima_out->xlen, ima_out->ylen,
                              0, ima_out->ylen-sizeOveY, ima_out->xlen,
                              sizeOveY, oveScanY)) == VM_FALSE) {
          cpl_msg_error(modName, 
            "Cannot create overScan : Failure in insertFloatImage");
          return NULL;
        }
        cpl_free(oveScanY);
      }
    }
    deletePortList(ports);
  }
  return ima_out;
}


/**
 * @memo 
 *   Transform mask coordinates (mm) into CCD coordinates (pixels)
 *
 * @return List of CCD coordinates
 * 
 * @param masklist      list of mask coordinates 
 * @param noslits       number of points in masklist
 * @param maskimaDesc   mask image descriptors
 *
 * @doc
 *   Read Mask To CCD transformation coefficients in the image header
 *   and transform mask coordinates of masklist into pixels
 *
 * @author P. Sartoretti 
 */    
VimosPixel *
MaskToCcd(VimosPixel *masklist, int noslit, VimosDescriptor *maskimaDescs)
{
  char        modName[] = "MaskToCcd";

  VimosPixel *ccdlist;
  char        valcoefX[80], valcoefY[80];
  char        comment[COMMENT_LENGTH];
  int         xord, yord, i, j, n;

 /*
  * Set flag = 1 to
  * read the coefficients/compute the poly that are read by the MPS. 
  * It is not a standard 3 degree polynomial, but it is an incomplete 
  * 6 deg polynomial containing terms until the 3 degree in x AND in y
  * (for ex. terms x^3y^3, but not terms like x^6, y^4...)
  * REMEMBER to do the same in VmImCalOpt
  */

  int         flag = 1;


 /* 
  * Check input 
  */
 
  cpl_msg_debug(modName,"transorming mm to pixels");

  if (masklist == NULL) {
    cpl_msg_error(modName, "No list of mask coordinates to convert"); 
    return NULL;
  }
  if  (maskimaDescs == NULL) {
    cpl_msg_error(modName, 
      "No image descriptors: cannot read coefficients for conversion");
    return NULL;
  }

 /* 
  * Allocate output list 
  */
  
  ccdlist = newPixel(noslit);
 
  if (readIntDescriptor(maskimaDescs, pilTrnGetKeyword("MaskCcdXord"), 
                        &xord, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Integer descriptor %s not found", 
                pilTrnGetKeyword("MaskCcdXord"));
    return NULL;
  }
    
  if (readIntDescriptor(maskimaDescs, pilTrnGetKeyword("MaskCcdYord"), &yord,
                        comment) == VM_FALSE) {
    cpl_msg_error(modName, "Integer descriptor %s not found", 
                pilTrnGetKeyword("MaskCcdYord"));
    return NULL;
  }

  for (n = 0; n < noslit; n++) {

    if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("MaskCcdX0"), 
                             valcoefX, comment) == VM_TRUE) {
      ccdlist[n].x += atof(valcoefX);
    }
    else { 
      cpl_msg_error(modName, "String descriptor %s not found", 
                  pilTrnGetKeyword("MaskCcdX0"));
      return NULL;
    }
    if (readStringDescriptor(maskimaDescs,  pilTrnGetKeyword("MaskCcdXX"), 
                             valcoefX, comment) == VM_TRUE) {
      ccdlist[n].x += atof(valcoefX) * masklist[n].x;
    }
    else { 
      cpl_msg_error(modName, "String descriptor %s not found", 
                  pilTrnGetKeyword("MaskCcdXX"));
      return NULL;
    }
    if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("MaskCcdXY"), 
                             valcoefX, comment) == VM_TRUE) {
      ccdlist[n].x += atof(valcoefX) * masklist[n].y;
    }
    else { 
      cpl_msg_error(modName, "String descriptor %s not found", 
                  pilTrnGetKeyword("MaskCcdXY"));
      return NULL;
    }
    if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("MaskCcdY0"), 
                             valcoefY, comment) == VM_TRUE) {
      ccdlist[n].y += atof(valcoefY);
    }
    else { 
      cpl_msg_error(modName, "String descriptor %s not found", 
                  pilTrnGetKeyword("MaskCcdY0"));
      return NULL;
    }
    if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("MaskCcdYY"), 
                             valcoefY, comment) == VM_TRUE) {
      ccdlist[n].y += atof(valcoefY) * masklist[n].y;
    }
    else { 
      cpl_msg_error(modName, "String descriptor %s not found", 
                  pilTrnGetKeyword("MaskCcdYY"));
      return NULL;
    }
    if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("MaskCcdYX"), 
                             valcoefY, comment)) {
      ccdlist[n].y += atof(valcoefY) * masklist[n].x;
    }
    else { 
      cpl_msg_error(modName, "String descriptor %s not found", 
                  pilTrnGetKeyword("MaskCcdYX"));
      return NULL;
    }
    if (flag) {
      for (i = 0; i <= xord; i++) {
        for (j = 0; j <= xord; j++) {
          if (readStringDescriptor(maskimaDescs, 
                                   pilTrnGetKeyword("MaskCcdX", i, j),
                                   valcoefX, comment) == VM_TRUE) {
            ccdlist[n].x += atof(valcoefX) * ipow(masklist[n].x, j) * 
                            ipow(masklist[n].y, i);
          }
          else { 
            cpl_msg_error(modName, "String descriptor %s not found", 
                        pilTrnGetKeyword("MaskCcdX", i, j));
            return NULL;
          }
        }
      }
      for (i = 0; i <= yord; i++) {
        for (j = 0; j <= yord; j++) { 
          if (readStringDescriptor(maskimaDescs,
                                   pilTrnGetKeyword("MaskCcdY", i, j),
                                   valcoefY, comment) == VM_TRUE) {
            ccdlist[n].y += atof(valcoefY) * ipow(masklist[n].x, j) * 
                            ipow(masklist[n].y, i);
          }
          else {
            cpl_msg_error(modName, "String descriptor %s not found", 
                        pilTrnGetKeyword("MaskCcdY", i, j));
            return NULL;
          }
        }
      }
    }
    else {
      for (i = 0; i <= xord; i++) {
        for (j = 0; j <= xord-i; j++) {
          if (readStringDescriptor(maskimaDescs, 
                                   pilTrnGetKeyword("MaskCcdX", i, j),
                                   valcoefX, comment) == VM_TRUE) {
            ccdlist[n].x += atof(valcoefX) * ipow(masklist[n].x, j) * 
                            ipow(masklist[n].y, i);
          }
          else { 
            cpl_msg_error(modName, "String descriptor %s not found", 
                        pilTrnGetKeyword("MaskCcdX", i, j));
            return NULL;
          }
        }
      }
      for (i = 0; i <= yord; i++) {
        for (j = 0; j <= yord-i; j++) { 
          if (readStringDescriptor(maskimaDescs,
                                   pilTrnGetKeyword("MaskCcdY", i, j),
                                   valcoefY, comment) == VM_TRUE) {
            ccdlist[n].y += atof(valcoefY) * ipow(masklist[n].x, j) * 
                            ipow(masklist[n].y, i);
          }
          else {
            cpl_msg_error(modName, "String descriptor %s not found",
                        pilTrnGetKeyword("MaskCcdY", i, j));
            return NULL;
          }
        }
      }
    }
  }
  return ccdlist;
}


/**
 * @memo 
 *   Transform CCD coordinates (pixels) into mask coordinates (mm)
 *
 * @return List of mask coordinates
 * 
 * @param ccdlist      list of ccd coordinates 
 * @param npoints      number of points in ccdlist
 * @param maskimaDescs mask image descriptors
 *
 * @doc
 *   Read CCD to Mask transformation coefficients in the image header
 *   and transform CCD coordinates of ccdList into mask coordinates
 *
 * @author P. Sartoretti  modified B.Garilli (now reading is in another s/r)
 */   

VimosPixel *
CcdToMask(VimosPixel *ccdlist, int npoints, VimosDescriptor *maskimaDescs)
{

  VimosPixel *masklist;
  char        valcoefX[80],valcoefY[80];
  char        comment[COMMENT_LENGTH];
  char        modName[] = "CcdToMask";
  int         xord,yord,i,j,k,n;

  double *x_coeff;
  double *y_coeff;
  double scale;

 /*
  * Set flag = 1 to
  * read the coefficients/compute the poly that are read by the MPS. 
  * It is not a standard 3 degree polynomial, but it is an incomplete 
  * 6 deg polynomial containing terms until the 3 degree in x AND in y
  * (for ex. terms x^3y^3, but not terms like x^6, y^4...)
  */

  int         flag = 1;


 /* 
  * Check input 
  */
 
  cpl_msg_debug(modName,"transforming pixels to mm");

  if (ccdlist == NULL) {
    cpl_msg_error(modName, "No list of CCD coordinates to convert"); 
    return NULL;
  }

  if (maskimaDescs == NULL) {
    cpl_msg_error(modName, 
      "No image descriptors: can not read coefficients for conversion");
    return NULL;
  }


 /* 
  * Allocate output list 
  */

  masklist = newPixel(npoints);


  if (readIntDescriptor(maskimaDescs, pilTrnGetKeyword("CcdMaskXord"), &xord,
                        comment) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found", 
                pilTrnGetKeyword("CcdMaskXord"));
    return NULL;
  }
    
  if (readIntDescriptor(maskimaDescs, pilTrnGetKeyword("CcdMaskYord"), &yord,
                        comment) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found", 
                pilTrnGetKeyword("CcdMaskYord"));
    return NULL;
  }

  n = (xord + 1) * (xord + 1);
  n += 3;
  x_coeff = (double*) cpl_calloc(n,sizeof(double));
  n = (yord + 1) * (yord + 1);
  n += 3;
  y_coeff = (double*) cpl_calloc(n, sizeof(double));




  if (flag) {
    if (!readMaskCcd(maskimaDescs, x_coeff, y_coeff, &scale)) {
      cpl_msg_error(modName, 
                  "Could not read coefficients for conversion");
      return NULL;
    }
    for (n = 0; n < npoints; n++) {
      masklist[n].x = x_coeff[0] + x_coeff[1] * ccdlist[n].x + 
        x_coeff[2]*ccdlist[n].y;
      masklist[n].y = y_coeff[0] + y_coeff[1] * ccdlist[n].y + 
        y_coeff[2]*ccdlist[n].x;
      
      
      k = 3;
      for (i = 0; i <= xord; i++) {
        for (j = 0; j <= xord; j++,k++) {
          masklist[n].x += x_coeff[k] * ipow(ccdlist[n].x, j) *
                           ipow(ccdlist[n].y, i);  
        }
      }
      masklist[n].x = masklist[n].x * scale;
      k = 3;
      for (i = 0; i <= yord; i++) {
        for (j = 0; j <= yord; j++, k++) {
          masklist[n].y += y_coeff[k] * ipow(ccdlist[n].x,j) *
                           ipow(ccdlist[n].y, i);
        }
      }
      masklist[n].y= masklist[n].y * scale;
    }

  } else {
    /* This part is unchanged (i.e. for each point re-reads coeffs)
       As it is not used, I did not bother to change it BG */
    for (n = 0; n < npoints; n++) {
      if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("CcdMaskX0"), 
                               valcoefX, comment) == VM_TRUE) {
        masklist[n].x += atof(valcoefX);
      }
      if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("CcdMaskXX"), 
                               valcoefX, comment) == VM_TRUE) {
        masklist[n].x += atof(valcoefX) * ccdlist[n].x;
      }
      if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("CcdMaskXY"), 
                               valcoefX, comment) == VM_TRUE) {
        masklist[n].x += atof(valcoefX) * ccdlist[n].y;
      }
      if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("CcdMaskY0"), 
                               valcoefY, comment) == VM_TRUE) {
        masklist[n].y += atof(valcoefY);
      }
      if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("CcdMaskYY"), 
                               valcoefY, comment) == VM_TRUE) {
        masklist[n].y += atof(valcoefY) * ccdlist[n].y;
      }
      if (readStringDescriptor(maskimaDescs, pilTrnGetKeyword("CcdMaskYX"), 
                               valcoefY, comment) == VM_TRUE) {
        masklist[n].y += atof(valcoefY) * ccdlist[n].x;
      }
      
      for (i = 0; i <= xord; i++) {
        for (j = 0; j <= xord - i; j++) {
          if (readStringDescriptor(maskimaDescs, 
                                   pilTrnGetKeyword("CcdMaskX", i, j),
                                   valcoefX, comment) == VM_TRUE) {
            masklist[n].x += atof(valcoefX) * ipow(ccdlist[n].x, j) * 
                             ipow(ccdlist[n].y, i);
          }
          else { 
            cpl_msg_error(modName, "Descriptor %s not found", 
                        pilTrnGetKeyword("CcdMaskX", i, j));
            return NULL;
          }
        }
      }

      for (i = 0; i <= yord; i++) {
        for (j = 0; j <= yord - i; j++) {
          if (readStringDescriptor(maskimaDescs, 
                                   pilTrnGetKeyword("CcdMaskY", i, j),
                                   valcoefY, comment) == VM_TRUE) {
            masklist[n].y += atof(valcoefY) * ipow(ccdlist[n].x, j) * 
                             ipow(ccdlist[n].y, i);
          }
          else {
            cpl_msg_error(modName, "Descriptor %s not found", 
                        pilTrnGetKeyword("CcdMaskY", i, j));
            return NULL;
          }
        }
      }
    }
  }

  return masklist;

}


VimosBool
readMaskCcd(VimosDescriptor *desc, double *x_coeff, double *y_coeff,
            double *scale)

{
  const char modName[] = "readMaskCcd";

  int      quad, i, k, j;
  int      xord,yord;

  char     comment[COMMENT_LENGTH];
  char       valcoef[80] = "0.";
  double   tCoeff = 6.0e-4;
  double   focuTempValue,tRef;
 /*
  *  Note: the optical distorsion model is made of two 2D distorsion
  *  models, one for the output X coordinate, and the other for the
  *  output Y coordinate. Both models have an equal "order" for the
  *  input x and y. Note also that when an Optical Distorsion Model
  *  is read, the offsets are always left to zero (as in the header
  *  no corresponding keyword is present). This implies that offsets
  *  are _not_ foreseen for the DRS optical distorsion models.  (C.Izzo)
  */

  if (readIntDescriptor(desc, pilTrnGetKeyword("CcdMaskXord"), &xord, 
                        comment) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("CcdMaskXord"));
    return VM_FALSE;
  }
  
  if (readIntDescriptor(desc, pilTrnGetKeyword("CcdMaskYord"), &yord,
                        comment) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found", 
                pilTrnGetKeyword("CcdMaskYord"));
    return VM_FALSE;
  }



    if (readStringDescriptor(desc, pilTrnGetKeyword("CcdMaskX0"), 
                             valcoef, comment) == VM_TRUE) {
      x_coeff[0] = atof(valcoef);
    }
    if (readStringDescriptor(desc, pilTrnGetKeyword("CcdMaskXX"), 
                             valcoef, comment) == VM_TRUE) {
      x_coeff[1] = atof(valcoef);
    }
    if (readStringDescriptor(desc, pilTrnGetKeyword("CcdMaskXY"), 
                             valcoef, comment) == VM_TRUE) {
      x_coeff[2] = atof(valcoef);
    }
    if (readStringDescriptor(desc, pilTrnGetKeyword("CcdMaskY0"), 
                             valcoef, comment) == VM_TRUE) {
      y_coeff[0] = atof(valcoef);
    }
    if (readStringDescriptor(desc, pilTrnGetKeyword("CcdMaskYY"), 
                             valcoef, comment) == VM_TRUE) {
      y_coeff[1] = atof(valcoef) ;
    }
    if (readStringDescriptor(desc, pilTrnGetKeyword("CcdMaskYX"), 
                             valcoef, comment) == VM_TRUE) {
      y_coeff[2] = atof(valcoef);
    }


 k=3;
 for(i=0; i <= xord; i++) {
   for(j=0; j <= xord; j++,k++) {
     if (!readStringDescriptor(desc, pilTrnGetKeyword("CcdMaskX",i,j), valcoef, comment)) {
       cpl_msg_warning(modName, "X Coefficient %d %d of the CCD-SKY transformation not found ", i,j);
       return VM_FALSE;
     }
     x_coeff[k] = atof(valcoef);
   }
 }
 k=3;
 for(i=0; i <= yord; i++) {
   for(j=0; j <= yord; j++,k++) {
     if (!readStringDescriptor(desc, pilTrnGetKeyword("CcdMaskY",i,j), valcoef, comment)) {
       cpl_msg_warning(modName, "Y Coefficient %d %d of the CCD-SKY transformation not found ", i,j);
       return VM_FALSE;
     }
     y_coeff[k] = atof(valcoef);
   }
 }

 if ((readIntDescriptor(desc, pilTrnGetKeyword("Quadrant"), &quad,
                          comment)) == VM_FALSE) {
     cpl_msg_error(modName, "Cannot read %s", pilTrnGetKeyword("Quadrant"));
     return VM_FALSE;
 }

 if (!readDoubleDescriptor(desc,pilTrnGetKeyword("CcdMaskTemp"),
                           &tRef, comment)) {
     cpl_msg_warning(modName,"Cannot find descriptor %s",
                   pilTrnGetKeyword("CcdSkyTemp"));
     return VM_FALSE;
 }

 if (!readDoubleDescriptor(desc, pilTrnGetKeyword("BeamTemperature", quad),
                           &focuTempValue, comment)) {
   cpl_msg_warning(modName,"Cannot find descriptor %s",
                 pilTrnGetKeyword("BeamTemperature", quad));
   return VM_FALSE;
 }

 /* Compute scale factor: this is copied from vmmcs */
 *scale = 1.+tCoeff*(focuTempValue-tRef);

  
  return VM_TRUE;
}


/**
 * @memo 
 *   Compute the mean gain factor of all image readout ports.
 *
 * @return Mean gain factor (e-/ADU), a negative number in case
 *         of failure.
 * 
 * @param image Input image
 *
 * @doc
 *   Search in the image header the gain factor for each port,
 *   and return the average value.
 *
 * @author C.Izzo
 */    

double
getMeanGainFactor(VimosImage *image)
{

    char       comment[COMMENT_LENGTH];

    int        i, nports;

    double     gain;
    double     sumGain = 0.0;

    VimosPort *ports;

    ports = getPorts(image, &nports);
    if (ports) {
        sumGain = 0.0;
        for (i = 1; i <= nports; i++) {
            if (readDoubleDescriptor(image->descs, 
                                     pilTrnGetKeyword("Adu2Electron", i), 
                                     &gain, comment) == VM_FALSE) {
                gain = -1.0;
                return gain;
            }
            sumGain += gain;
        }
        gain = sumGain / nports;
  
        deletePortList(ports);
    }

    return gain;

}


/**
 * @memo
 *   Create the PAF files containing the CCD to Mask transformation
 *   for X and Y coordinates for a given quadrant.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param descs      List of descriptor containing the CCD to Mask
 *                   transformation coefficients
 * @param namePAF    PAF file root name
 * @param pafName    Returned created PAF files names
 *
 * @doc
 *   Both PAF files are created with a name depending on the coordinate
 *   transformation and on the quadrant; a common header is supplied
 *   and the polynomial transformation coefficients are written.
 *
 * @author C. Izzo
 */

int
createMaskCcdPAF(VimosDescriptor *descs, const char *namePAF, char **pafName)
{
  const char       modName[] = "createMaskCcdPAF";
  char            *pafName_noext;
  VimosDescriptor *desc;
  char             comment[COMMENT_LENGTH];
  FILE            *fp;
  int              xOrd, yOrd;
  int              i, j, len, quad;

  cpl_msg_debug(modName,"write CCD2Mask into PAF file");

  readIntDescriptor(descs, pilTrnGetKeyword("Quadrant"), &quad, comment);

  len = sizeof(char) * (strlen(namePAF)+7);
  *pafName = (char *)cpl_malloc(len);
  sprintf(*pafName, "%s_%d.cmf", namePAF, quad);

  fp = fopen(*pafName, "w");

  if (!fp) 
    return EXIT_FAILURE;


 /*
  *  Write header
  */

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafHeaderStart"),NULL);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafType"),"Configuration");
  pafName_noext = (char *)cpl_malloc(len - 4);
  sprintf(pafName_noext, "%s_%d", namePAF, quad);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafId"), pafName_noext);
  cpl_free(pafName_noext);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafName"), *pafName);
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
    return EXIT_FAILURE;
  }

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskDate"), 
                      desc->descValue->s);
 /*
  * Temp
  */

  desc = findDescriptor(descs, pilTrnGetKeyword("BeamTemperature", quad));

  if (!desc) {
    cpl_free(pafName);
    return EXIT_FAILURE;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskTemp"), 
                      desc->descValue->d);


 /*
  *  Write transformation Mask to CCD coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("MaskCcdXord"), &xOrd, comment)
                        == VM_FALSE) {
    cpl_free(pafName);
    return EXIT_FAILURE;
  }

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdXord"), xOrd);

  if (readIntDescriptor(descs, pilTrnGetKeyword("MaskCcdYord"), &yOrd, comment)
                        == VM_FALSE) {
    cpl_free(pafName);
    return EXIT_FAILURE;
  } 

  writeIntPAFEntry(fp, (char *) pilTrnGetKeyword("PAFMaskCcdYord"), yOrd);

  if ((desc = findDescriptor(descs,pilTrnGetKeyword("MaskCcdX0")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdX0"), 
                        desc->descValue->s);
  }
  
  if ((desc = findDescriptor(descs,pilTrnGetKeyword("MaskCcdXX")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdXX"), 
                        desc->descValue->s);
  }
  
  if ((desc = findDescriptor(descs, pilTrnGetKeyword("MaskCcdXY")))) {
    writeStringPAFEntry(fp, (char *) pilTrnGetKeyword("PAFMaskCcdXY"), 
                        desc->descValue->s);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("MaskCcdY0")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdY0"), 
                        desc->descValue->s);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("MaskCcdYY")))) {
    writeStringPAFEntry(fp, (char *) pilTrnGetKeyword("PAFMaskCcdYY"), 
                        desc->descValue->s);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("MaskCcdYX")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdYX"), 
                        desc->descValue->s);
  }

  for (i = 0; i <= xOrd; i++)  {
    for (j = 0; j <= xOrd; j++)  {
      if ((desc = findDescriptor(descs, pilTrnGetKeyword("MaskCcdX",i,j)))) {
        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdX",i,j), 
                            desc->descValue->s);
      }
      else {
        cpl_msg_error(modName, 
                    "Missing Mask-to-CCD transformation coefficient %s", 
                    desc->descName);
      }
    }
  }

  for (i = 0; i <= yOrd; i++)  {
    for (j = 0; j <= yOrd; j++)  {
      if ((desc = findDescriptor(descs, pilTrnGetKeyword("MaskCcdY",i,j)))) {
        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdY",i,j), 
                            desc->descValue->s);
      }
      else {
        cpl_msg_error(modName, 
                    "Missing Mask-to-CCD transformation coefficient %s", 
                    desc->descName);
      }
    }
  }

  if ((desc = findDescriptor(descs,pilTrnGetKeyword("MaskCcdXrms")))) {
    writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdXrms"), 
                        desc->descValue->d);
  }

  if ((desc = findDescriptor(descs,pilTrnGetKeyword("MaskCcdYrms")))) {
    writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFMaskCcdYrms"), 
                        desc->descValue->d);
  }


 /*
  *  Write transformation CCD to Mask coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("CcdMaskXord"), &xOrd, comment)
                        == VM_FALSE) {
    return EXIT_FAILURE;
  } 
  else {
    writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskXord"), xOrd);
  }

  if (readIntDescriptor(descs, pilTrnGetKeyword("CcdMaskYord"), &yOrd, comment)
                        == VM_FALSE) {
    return EXIT_FAILURE;
  } 
  else {
    writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskYord"), yOrd);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskX0")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskX0"), 
                        desc->descValue->s);
  }
  
  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskXX")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskXX"), 
                        desc->descValue->s);
  }
  
  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskXY")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskXY"), 
                        desc->descValue->s);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskY0")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskY0"), 
                        desc->descValue->s);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskYY")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskYY"), 
                        desc->descValue->s);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskYX")))) {
    writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskYX"), 
                        desc->descValue->s);
  }


  for (i = 0; i <= xOrd; i++)  {
    for (j = 0; j <= xOrd; j++)  {
      if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskX",i,j)))) {
        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskX",i,j), 
                            desc->descValue->s);
      }
      else {
        cpl_msg_error(modName, 
                    "Missing Mask-to-CCD transformation coefficient %s", 
                    desc->descName);
      }
    }
  }

  for (i = 0; i <= yOrd; i++)  {
    for (j = 0; j <= yOrd; j++)  {
      if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskY",i,j)))) {
        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskY",i,j), 
                            desc->descValue->s);
      }
      else {
        cpl_msg_error(modName, 
                    "Missing Mask-to-CCD transformation coefficient %s", 
                    desc->descName);
      }
    }
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskXrms")))) {
    writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskXrms"), 
                        desc->descValue->d);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdMaskYrms")))) {
    writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdMaskYrms"), 
                        desc->descValue->d);
  }

  fclose(fp);
  return EXIT_SUCCESS;
}


/**
 * @memo
 *   Create the PAF files containing the CCD to Sky transformation
 *   for X and Y coordinates for a given quadrant.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param descs      List of descriptor containing the CCD to Mask
 *                   transformation coefficients
 * @param namePAF    PAF file root name
 * @param pafName    Returned created PAF files names
 *
 * @doc
 *   Both PAF files are created with a name depending on the coordinate
 *   transformation and on the quadrant; a common header is supplied
 *   and the polynomial transformation coefficients are written.
 *
 * @author B.Garilli
 */

int
createCcdSkyPAF(VimosDescriptor *descs, char *namePAF, char **pafName)
{
  const char       modName[] = "createCcdSkyPAF";
  char            *pafName_noext;
  VimosDescriptor *desc;
  char             comment[COMMENT_LENGTH];
  FILE            *fp;
  int              xOrd, yOrd;
  int              i, j, len, quad;


  cpl_msg_debug(modName, "write CCD2Sky into PAF file");

  readIntDescriptor(descs, pilTrnGetKeyword("Quadrant"), &quad, comment);

  len = sizeof(char) * (strlen(namePAF)+7);
  *pafName = (char *)cpl_malloc(len);
  sprintf(*pafName, "%s_%d.cmf", namePAF, quad);

  fp = fopen(*pafName, "w");

  if (!fp) 
    return EXIT_FAILURE;


 /*
  *  Write header 
  */

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafHeaderStart"),NULL);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafType"),"Configuration");
  pafName_noext = (char *)cpl_malloc(len - 4);
  sprintf(pafName_noext, "%s_%d", namePAF, quad);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafId"), pafName_noext);
  cpl_free(pafName_noext);
  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PafName"), *pafName);
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
    return EXIT_FAILURE;
  }

  writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdSkyDate"), 
                      desc->descValue->s);


 /*
  * Temp %%%
  */

  desc = findDescriptor(descs, pilTrnGetKeyword("BeamTemperature", quad));

  if (!desc) {
    cpl_free(pafName);
    return EXIT_FAILURE;
  }

  writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdSkyTemp"),
                      desc->descValue->d);


 /*
  *  Write transformation CCD to Sky coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("CcdSkyXord"), &xOrd, comment)
                        == VM_FALSE) {
    cpl_free(pafName);
    return EXIT_FAILURE;
  } 
  
  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdSkyXord"), xOrd);

  if (readIntDescriptor(descs, pilTrnGetKeyword("CcdSkyYord"), &yOrd, comment)
                        == VM_FALSE) {
    cpl_free(pafName);
    return EXIT_FAILURE;
  } 

  writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdSkyYord"), yOrd);


  for (i = 0; i <= xOrd; i++)  {
    for (j = 0; j <= xOrd; j++)  {
      if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdSkyX",i,j)))) {
        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdSkyX",i,j), 
                            desc->descValue->s);
      }
      else {
        cpl_msg_error(modName, 
                    "Missing Sky-to-CCD transformation coefficient %s", 
                    desc->descName);
      }
    }
  }

  for (i = 0; i <= yOrd; i++)  {
    for (j = 0; j <= yOrd; j++)  {
      if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdSkyY",i,j)))) {
        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdSkyY",i,j), 
                            desc->descValue->s);
      }
      else {
        cpl_msg_error(modName, 
                    "Missing Sky-to-CCD transformation coefficient %s", 
                    desc->descName);
      }
    }
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdSkyXrms")))) {
    writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdSkyXrms"), 
                        desc->descValue->d);
  }

  if ((desc = findDescriptor(descs, pilTrnGetKeyword("CcdSkyYrms")))) {
    writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFCcdSkyYrms"), 
                        desc->descValue->d);
  }

 /*
  *  Write transformation Sky to CCD coefficients
  */

  if (readIntDescriptor(descs, pilTrnGetKeyword("SkyCcdXord"), &xOrd, comment)
                        == VM_FALSE) {
    return EXIT_FAILURE;
  } 
  else {
    writeIntPAFEntry(fp, (char *)pilTrnGetKeyword("PAFSkyCcdXord"), xOrd);
  }

  if (readIntDescriptor(descs, pilTrnGetKeyword("SkyCcdYord"), &yOrd, comment)
                        == VM_FALSE) {
    return EXIT_FAILURE;
  } 
  else {
    writeIntPAFEntry(fp, (char *) pilTrnGetKeyword("PAFSkyCcdYord"), yOrd);
  }

 

  for (i = 0; i <= xOrd; i++)  {
    for (j = 0; j <= xOrd; j++)  {
      if ((desc = findDescriptor(descs, pilTrnGetKeyword("SkyCcdX",i,j)))) {
        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFSkyCcdX",i,j), 
                            desc->descValue->s);
      }
      else {
        cpl_msg_error(modName, 
                    "Missing Sky-to-CCD transformation coefficient %s", 
                    desc->descName);
      }
    }
  }

  for (i = 0; i <= yOrd; i++)  {
    for (j = 0; j <= yOrd; j++)  {
      if ((desc = findDescriptor(descs, pilTrnGetKeyword("SkyCcdY",i,j)))) {
        writeStringPAFEntry(fp, (char *)pilTrnGetKeyword("PAFSkyCcdY",i,j), 
                            desc->descValue->s);
      }
      else {
        cpl_msg_error(modName, 
                    "Missing Sky-to-CCD transformation coefficient %s", 
                    desc->descName);
      }
    }
  }

  if ((desc = findDescriptor(descs,pilTrnGetKeyword("SkyCcdXrms")))) {
    writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFSkyCcdXrms"), 
                        desc->descValue->d);
  }

  if ((desc = findDescriptor(descs,pilTrnGetKeyword("SkyCcdYrms")))) {
    writeDoublePAFEntry(fp, (char *)pilTrnGetKeyword("PAFSkyCcdYrms"), 
                        desc->descValue->d);
  }



  fclose(fp);
  return EXIT_SUCCESS;
}


/**
 * @memo
 *   Compute the average readout noise in adu for an image.
 *
 * @return The value of the computed readout noise is returned. In case
 *   of any error the return value is set to -1.
 *
 * @param image  Pointer to an existing image.
 *
 * @doc
 *   The function first determines the used readout ports from
 *   the header keywords. For each readout port the readout noise
 *   is determined from the overscan areas. Finally the average
 *   readout noise is determined from the individual values.
 *
 *   Note that the computed average readout noise is measured in
 *   units of adu!
 *
 * @see estimateImageRon
 *
 * @author R. Palsa
 */

double
computeAverageRon(VimosImage *image)
{

  int portCount = 0;

  double meanRon;

  VimosFloatArray *ronList;

  VimosPort *portList = NULL;



  /*
   * Get the list of readout ports
   */

  if (!(portList = getPorts(image, &portCount)))
    return -1;


  /*
   * Estimate the readout noise for each individual port and calculate
   * the average readout noise.
   */

  if (!(ronList = estimateImageRon(image, portList))) {
    deletePortList(portList);
    return -1;
  }

  meanRon = computeAverageFloat(ronList->data, ronList->len);


  /*
   * Cleanup
   */

  deletePortList(portList);
  deleteFloatArray(ronList);

  return meanRon;

}


/**
 * @memo
 *   Read the average readout noise in ADU from an image keyword header.
 *
 * @return The value of the mean readout noise read from header is returned. 
 *   In case of any error the return value is set to -1.
 *
 * @param image  Pointer to an existing image.
 *
 * @doc
 *   The function first determines the used readout ports from
 *   the header keywords. For each readout port the readout noise
 *   is determined from the corresponding header entry. Finally 
 *   the average readout noise is determined from the individual 
 *   values.
 *
 *   Note that the computed average readout noise is measured in
 *   units of ADU!
 *
 * @see getImageRon
 *
 * @author C. Izzo
 */

double
getAverageRon(VimosImage *image)
{

  double meanRon;
  int i;
  VimosFloatArray *ronList;


  /*
   * Estimate the readout noise for each individual port and calculate
   * the average readout noise.
   */

  if (!(ronList = getImageRon(image)))
    return -1;

  for (i=0; i<ronList->len; i++) printf("*** %f ***\n", ronList->data[i]);

  meanRon = computeAverageFloat(ronList->data, ronList->len);


  /*
   * Cleanup
   */

  deleteFloatArray(ronList);

  return meanRon;

}


/**
 * @memo
 *   Evaluate the average noise in ADU for a raw data image.
 *
 * @return The value of the computed noise is returned. 
 *   In case of any error the return value is set to -1.
 *
 * @param image   Pointer to an existing image.
 * @param meanRon Average image RON, in ADU. If negative, is computed
 *                within the function.
 * @param gain    Gain factor (e-/ADU).
 *
 * @doc
 *   The function first determines the used readout ports from
 *   the header keywords. The mean bias level and readout noise, 
 *   both in ADU, are determined within the overscan regions. 
 *   The rest of the image is extracted, subtracted the mean bias 
 *   level, its flux converted into electrons, square-rooted, and 
 *   finally averaged and converted back to ADU. Finally, the total 
 *   image noise is computed as the geometric mean of this noise 
 *   and RON. This is an estimate of the image typical noise, 
 *   outside the prescan/overscan regions.
 *
 *   Note that the input image is assumed to be a raw, not reduced
 *   image, with intensities in ADU and still having overscan regions. 
 *
 * @see estimateImageRon
 *
 * @author C. Izzo
 */

double
evaluateAverageNoise(VimosImage *image, float meanRon, float gain)
{

  const char  modName[] = "evaluateAverageNoise";
  int         portCount = 0;
  int         i;
  int         startPortsX, startPortsY, endPortsX, endPortsY;
  int         sizePortsX, sizePortsY, sizePorts;
  float      *trimdata;
  VimosPort  *ports = NULL;
  VimosPort  *currport = NULL;
  VimosImage *imageMinusBias = NULL;
  double      meanShotNoise, meanNoise, meanLevel;


  if ((ports = getPorts(image, &portCount)) == NULL) {
    cpl_msg_debug(modName, "Failure in reading port structure");
    return -1.;
  }
  
  for (currport = ports; currport; currport = currport->next) {

    if (currport->prScan->nX + currport->ovScan->nX)
      break;

    cpl_msg_debug(modName, "No prescans/overscans found.");
    deletePortList(ports);
    return -1.;
  }

  if (meanRon < 0.) {

   /*
    * First, determine the average RON
    */

    meanRon = computeAverageRon(image);

    if (meanRon < 0.) {

     /*
      * Cannot compute RON, probably because overscans are missing.
      * Then read it from the header.
      */

      cpl_msg_warning(modName, "RON level is read from header instead of "
                    "being computed");

      meanRon = getAverageRon(image);

      if (meanRon < 0.) {
        cpl_msg_debug(modName, "Failure determining RON.");
        deletePortList(ports);
        return -1.;
      }
    }
  }

 /*
  * Second, apply a rough bias subtraction
  */

  if (!(imageMinusBias = duplicateImage(image))) {
    deletePortList(ports);
    return -1.;
  }

  if (VM_FALSE == subtractOverscan(imageMinusBias->data, 
                  imageMinusBias->xlen, imageMinusBias->ylen, ports)) {
    cpl_msg_debug(modName, "Failure in subtracting mean bias.");
    deletePortList(ports);
    deleteImage(imageMinusBias);
    return -1.0;
  }

  sizePorts = getTotalReadoutWindow(ports, 
              &startPortsX, &startPortsY, &sizePortsX, &sizePortsY);

  deletePortList(ports);

  endPortsX = sizePortsX + startPortsX;
  endPortsY = sizePortsY + startPortsY;

  cpl_msg_debug(modName,
              "Extract image excluding overscans: start=(%d,%d) end=(%d,%d)",
              startPortsX, startPortsY, endPortsX, endPortsY);

  trimdata = extractFloatImage(imageMinusBias->data, image->xlen, 
                               image->ylen, startPortsX, startPortsY, 
                               sizePortsX, sizePortsY);

  deleteImage(imageMinusBias);

  if (!trimdata) {
    cpl_msg_debug(modName, "Failure in extracting image");
    return -1.0;
  }

 /*
  * Here trimdata values are scrambled - it's irrelevant for the
  * algorithm
  */

  meanLevel = computeAverageFloat(trimdata, sizePorts);

/* printf("*** Mean level: %f\n", meanLevel); */

  /*
   * Convert flux to electron counts, derive one sigma
   * noise as square root of flux, and convert back to ADU.
   */

  for (i = 0; i < sizePorts; i++) {
    if (trimdata[i] > 0.5) {
      trimdata[i] = sqrt(trimdata[i]*gain)/gain;
    }
    else {
      trimdata[i] = 1.0;
    }
  }

  meanShotNoise = computeAverageFloat(trimdata, sizePorts);

/* printf("*** Mean shot noise: %f\n", meanShotNoise); */

  cpl_free(trimdata);

  meanNoise = sqrt(meanShotNoise * meanShotNoise + meanRon * meanRon);

/* printf("*** Mean noise: %f\n", meanNoise); */

  return meanNoise;

}


/**
 * @memo
 *   Compute average airmass from an image.
 *
 * @return The function returns #EXIT_SUCCESS# if no error occured,
 *   otherwise the return value is #EXIT_FAILURE#.
 *
 * @param image    Image of the an observation.
 * @param airmass  Computed airmass.
 *
 * @doc
 *   The function gets the rightascension and declination of telescope
 *   pointing from the image \textbf{image}. Note that currently the
 *   offset of the field center of the individual quadrants from the
 *   pointing direction is not taken into account. Also the latitude of
 *   the telescope, the siderial time in seconds elapsed since siderial
 *   midnight and the observation's exposure time. If any of this
 *   parameters cannot be retrieved from the image the function returns
 *   an error.
 *
 *   Using these data the function computes the approximated, average 
 *   airmass of the observation and returns this value through the
 *   parameter \textbf{airmass}. In case of an error \textbf{aimass}
 *   is set to -1.
 *
 * @author R. Palsa
 */

int
VmComputeAirmass(VimosImage *image, double *airmass)
{

  const char fctid[] = "VmComputeAirmass";


  char comment[COMMENT_LENGTH];

  double alpha, delta, latitude;
  double siderialTime, exposureTime;
  double a;


  *airmass = -1;


  /*
   * Get telescope pointing, telescope latitude, siderial time, and
   * the exposure time from the image header.
   */

  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("Alpha"), &alpha,
                           comment) != VM_TRUE ||
      readDoubleDescriptor(image->descs, pilTrnGetKeyword("Delta"), &delta,
                           comment) != VM_TRUE) {
    cpl_msg_error(fctid, "Cannot get telescope pointing!");
    return EXIT_FAILURE;
  }

  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("Latitude"),
                           &latitude, comment) != VM_TRUE) {
    cpl_msg_error(fctid, "Cannot get telescope latitude!");
    return EXIT_FAILURE;
  }

  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("SiderialTime"),
                           &siderialTime, comment) != VM_TRUE) {
    cpl_msg_error(fctid, "Cannot get siderial time at observation start!");
    return EXIT_FAILURE;
  }

  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("ExposureTime"),
                           &exposureTime, comment) != VM_TRUE) {
    cpl_msg_error(fctid, "Cannot get exposure time of observation!");
    return EXIT_FAILURE;
  }


  /*
   * Compute average airmass.
   */

  if ((a = pilAstroComputeAirmass(alpha, delta, siderialTime,
                                  exposureTime, latitude)) < 0.) {
    cpl_msg_error(fctid, "Airmass computation failed!");
    return EXIT_FAILURE;
  }

  *airmass = a;

  return EXIT_SUCCESS;

}


/*
 * This function return a "reasonable" focus temperature for a given
 * quadrant, obtained by comparing it with the other focus temperatures
 * and the ambient temperatures. This is to add redundancy to the system 
 * in case a temperature sensor fails. 
 */

int
getBeamTemperature(VimosDescriptor *desc, double *temperature, 
                   double tolerance, int quadrant)
{

  const char id[] = "getBeamTemperature";

  double ta;        /* Ambient temperature         */
  double tb;        /* Other beams temperatures    */
  double tm;        /* Mean of good temperatures   */
  int    count;     /* Number of good temperatures */
  int    i;


  if (readDoubleDescriptor(desc, pilTrnGetKeyword("AmbientTemperature"),
                           &ta, NULL) == VM_FALSE) {
    cpl_msg_warning(id, "Cannot find descriptor %s",
                  pilTrnGetKeyword("AmbientTemperature"));
    return 1;
  }

  if (readDoubleDescriptor(desc, pilTrnGetKeyword("BeamTemperature", quadrant),
                           temperature, NULL) == VM_FALSE) {
    cpl_msg_warning(id, "Cannot find descriptor %s",
                  pilTrnGetKeyword("BeamTemperature", quadrant));
    return 1;
  }

  if (fabs(*temperature - ta) < tolerance)
    return 0;


  /*
   * If we reach here, beam temperature is beyond tolerance.
   * We select the temperatures from other beams that are within
   * tolerance, and replace the beam temperature with their mean value.
   */

  tm = 0.0;
  count = 0;
  for (i = 1; i <= 4; i++) {

    if (i == quadrant)
      continue;

    if (readDoubleDescriptor(desc, pilTrnGetKeyword("BeamTemperature", i),
                             &tb, NULL) == VM_FALSE) {
      /*
       * These are data that do not have in one header all four beams
       * temperatures. The old method is applied (for backward compatibility):
       */

      cpl_msg_warning(id, "Beam temperature (%f) out of range! "
                    "Using ambient temperature (%f) instead!",
                    *temperature, ta);
      *temperature = ta;

      return 0;

    }
    else {
      if (fabs(tb - ta) < tolerance) {
        tm += tb;
        count++;
      }
    }

  }

  if (count) {
    cpl_msg_warning(id, "Beam temperature (%f) out of range! "
                  "Using estimate from other beam temperatures (%f) instead!",
                  *temperature, tm / count);
    *temperature = tm / count;
  }
  else {

    /*
     * All temperatures are beyond tolerance. Could they possibly be
     * all wrong? Here we have a philosophical choice: either they
     * are correct, then we use the current beam temperature, or they
     * are wrong, then we use the ambient temperature. Currently the
     * latter is chosen. If the former were chosen, this else block
     * should be simply removed.
     */

    cpl_msg_warning(id, "Beam temperature (%f) out of range! "
                  "Using ambient temperature (%f) instead!",
                  *temperature, ta);
    *temperature = ta;
  }

  return 0;

}
/**@}*/
