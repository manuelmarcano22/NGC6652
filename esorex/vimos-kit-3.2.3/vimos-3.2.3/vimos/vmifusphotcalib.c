/* $Id: vmifusphotcalib.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <pilmessages.h>

#include "vmimage.h"
#include "vmimageset.h"
#include "vmmatrix.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmccdtable.h"
#include "vmifutable.h"
#include "vmobjecttable.h"
#include "vmifutypes.h"
#include "vmmossphotcalib.h"
#include "vmifusphotcalib.h"


#define NUMFIBS 6400
#define IFU_HEAD_SIZE 80


int
VmIfuCalPhot(VimosImageSet *imageSet, int fitOrd, float int_frac)
{
  int i, k, right, got, quadNum, ngood;
  int referenceL, referenceM, specLen, objNumber = 0;
  int fibL[NUMFIBS], fibM[NUMFIBS], fibLS[NUMFIBS];
  int fibMS[NUMFIBS], toSort[NUMFIBS];
  int refSlitNo = 0, refFibNo = 0, refQuad = 0, stdRow = 0, xLen, yLen;
  int jj, ll, ok;

  float referenceFlux, refTotIntensity = 0, wLenInc;
  float threshold, fraction;
  float fibIntFlux[NUMFIBS], fibIntFluxS[NUMFIBS];

  double dValue;
  char comment[80];
  char modName[] = "VmIfuCalPhot";

  VimosFloatArray *spectrum;
  VimosImage *ifuImage, *ifuMask, *stdImage = NULL;

  VimosSingleImage *oneImage;

  VimosObjectTable *objTable;
  VimosObjectObject *imaObjs;

  VimosTable *stdTable = NULL;
  VimosIfuTable *ifuTable;
  VimosIfuQuad *theQuads;
  VimosIfuSlit *theIfuSlits = NULL;
  VimosIfuFiber *theIfuFibers;

  VimosPixelRegion *pixRegion;
  VimosPixelList *thePixList;
  VimosPixelData *thePixData;

  pilMsgInfo (modName, "Start computing spectro-photometric calibration");

  puts("TBD: numbering for L,M from 0 or 1 (below) - CHECK IT");

  oneImage = imageSet->images;

  /* this is the subsctript for fibM, fibL, fibIntFlux arrays */
  got = 0;

  while (oneImage)
  {
    /* take the object table and IFU table of this image */
    objTable = oneImage->objectTable;
    ifuTable = oneImage->ifuTable;

    /* length of spectra in wavelength */  
    specLen = oneImage->theImage->xlen;
    
    /* Allocate spectrum */
    spectrum = newFloatArray(specLen);

    /* CHECK wLenInc DOUBLE OR FLOAT! */
    /* read wlen inc keyword: expected in nm (DRS Document) */
    readDoubleDescriptor(objTable->descs, "ESO PRO WLEN INC", 
                         &dValue, comment);
    wLenInc = (float) dValue;

    /* counter for all objects in image */  
    objNumber = 0;

    /* read which quadrant this image refers */
    readIntDescriptor(objTable->descs, "ESO PRO QUAD", &quadNum, comment);
    theQuads = ifuTable->quads;

    /* take from ifuTable only those slits in the quadrant corresponding */
    /* to the image in use */
    while (theQuads)
    {
      if (theQuads->quadNo == quadNum)
      {
        theIfuSlits = theQuads->ifuSlits;
      }
      theQuads = theQuads->next;
    }

    /* take objects in the image object table */
    imaObjs = objTable->objs;
    while (imaObjs)
    {
      /* loop on IFU slits for this quadrant to find L,M of each object */
      while (theIfuSlits)
      {
        /* if this is the "right" slit go on to take fibers */
        if (theIfuSlits->ifuSlitNo == imaObjs->IFUslitNo)
        {
          theIfuFibers = theIfuSlits->fibers;

          /* loop on fibers of this slit until find the one */
          /* corresponding to the object */
          while (theIfuFibers)
          {
            if (theIfuFibers->fibNo == imaObjs->IFUfibNo)
            {
              /* take this spectrum, to integrate */
              for (i = 0; i < specLen; i++ )
              {
                spectrum->data[i] = 
                    oneImage->theImage->data[i + ((imaObjs->rowNum)*
                        specLen)];
              }
              fibIntFlux[got] = integrateSpectrum(spectrum, wLenInc);

              fibL[got] = theIfuFibers->fiberL;
              fibM[got] = theIfuFibers->fiberM;
              got++;
            }
            theIfuFibers = theIfuFibers->next;
          }
        }

        theIfuSlits = theIfuSlits->next;
      }

      imaObjs = imaObjs->next;

      /* increment counter of objects done */
      objNumber++;

    }

    oneImage = oneImage->next; 
    deleteFloatArray(spectrum);
  }

  if (got != NUMFIBS) printf("WRONG NUMBER of DETECTED OBJECTS: %d",got);

  /* now retrieve array of subscripts to sort in flux */
  Indexx(NUMFIBS, fibIntFlux, toSort);

  /* actual sorting. This mantains correspondence between flux array and */
  /* L,M coordinates of fibers                                           */
  for (i=0; i<NUMFIBS; i++)
  {
    fibIntFluxS[i] = fibIntFlux[toSort[i]];
    fibLS[i] = fibL[toSort[i]];
    fibMS[i] = fibM[toSort[i]];
  }

  /* WARNING! CHECK THE SORTING */
  /* take the flux reference fiber, i.e. the one with the max. flux */
  referenceFlux = fibIntFluxS[0];
  referenceL = fibLS[0];
  referenceM = fibMS[0];

  ifuImage = newImageAndAlloc(IFU_HEAD_SIZE,IFU_HEAD_SIZE);

  /* NOTE: pixel (0,0) correspond to (L,M) = (1,1) */
  /* put the fiber (fibL, fibM) = (1,1) at the pixel (0,0) */
  /* CHECK IF THIS WORKS CORRECTLY */
  puts("Check renumbering of fibers -> pixels");
  for (k=0; k<NUMFIBS; k++)
  {
    ifuImage->data[(fibL[k]-1) + (fibM[k]-1)*IFU_HEAD_SIZE] = fibIntFlux[k];
  }

  /* now we can use the image for thresholding & looking for adjacent pixels*/
  ifuMask = newImageAndAlloc(IFU_HEAD_SIZE,IFU_HEAD_SIZE);

  /* set the threshold value */
  threshold = referenceFlux * int_frac;

  ngood = 0;
  thresholdImage(ifuImage, threshold, ifuMask, &ngood);

  /* now detect adjacent pixels, i.e. regions */
  pixRegion = newPixelRegion();

  pixRegion = findRegionsOnPixelMap(ifuImage, ifuMask, pixRegion,
      ngood);

  /* if more than one blob, look for the one with reference fiber */
  if (pixRegion->numRegions > 1)
  {
    puts("VmIfuCalPhot: number of regions above threshold > 1");
  }

  thePixList = pixRegion->pixelList;

  right = 0;
  while (thePixList)
  {
    thePixData = thePixList->pixelData;
    while (thePixData)
    {
      if ( (thePixData->x == referenceL) &&
          (thePixData->y == referenceM) ) right++;
      thePixData = thePixData->next;
    }
    if (right > 1) puts("something wrong in detecting reference region");
    if (right == 1) refTotIntensity = thePixList->totIntensity;

    thePixList = thePixList->next;
  }

  if(right == 0) puts("something wrong in detecting reference region");

  /* loop on IFU Table to find slit and fiber number of reference spectrum */
  oneImage = imageSet->images;
  ifuTable = oneImage->ifuTable;
  theQuads = ifuTable->quads;
  while (theQuads)
  {
    theIfuSlits = theQuads->ifuSlits;
    while (theIfuSlits)
    {
      theIfuFibers = theIfuSlits->fibers;
      while (theIfuFibers)
      {
        if ((theIfuFibers->fiberL == referenceL) &&
            (theIfuFibers->fiberM == referenceM))
        {
          refFibNo = theIfuFibers->fibNo;
          refSlitNo = theIfuSlits->ifuSlitNo;
          refQuad = theQuads->quadNo;
        }
        theIfuFibers = theIfuFibers->next;
      }
      theIfuSlits = theIfuSlits->next;
    }
    theQuads = theQuads->next;
  }

  /* determine the fraction factor: fraction of light lost by the ref. fibre*/
  fraction = refTotIntensity/referenceFlux;

  /* now loop on imageSet and objectTables to find the image with reference 
     spectrumand relative spectro-phot table, find row number for the 
     reference spectrum, correct flux loss for reference spectrum, and evaluate
     spectro-photometric calibration */

  ok = 0;
  oneImage = imageSet->images;
  while (oneImage)
  {
    /* take the object table and IFU table of this image */
    objTable = oneImage->objectTable;

    /* read which quadrant this image refers to */
    readIntDescriptor(objTable->descs, "ESO PRO QUAD", &quadNum, "");

    /* if this image refers to the quadrant with reference spectrum, go on */
    if (quadNum == refQuad)
    {
      /* take objects in the image object table */
      imaObjs = objTable->objs;
      while (imaObjs)
      {
        if ((imaObjs->IFUslitNo == refSlitNo) &&
            (imaObjs->IFUfibNo == refFibNo))
        {
          if (ok > 1) 
          {
            pilMsgError(modName, ": More than one image for reference spectrum. Exiting");
            return EXIT_FAILURE;
          }

          /* take row number of the std spectrum in the image */
          stdRow = imaObjs->rowNum;

          xLen = oneImage->theImage->xlen;
          yLen = oneImage->theImage->ylen;

          /* correct the ref. fibre spectrum for flux loss
		 before deriving response function as for MOS */

          /* CHECK if THIS IS OK! */
          puts("Check flux correction!");
          for (ll=0; ll<xLen; ll++)
            oneImage->theImage->data[ll + stdRow*yLen] *= fraction;

          stdImage = newImageAndAlloc(xLen, yLen);
          /* take image for spectro-photometric calibration */
          for (jj=0; jj<yLen; jj++)
          {
            for (ll=0; ll<xLen; ll++)
            {
              stdImage->data[ll +jj*yLen] = 
                  oneImage->theImage->data[ll +jj*yLen];
            }
          }

          /* take spectro-photometric table of this image */
          stdTable = oneImage->sphotStdTable;

          ok++;
        } /* end if this object is the reference spectrum */

        imaObjs = imaObjs->next;
      } /* end loop on objects for this image */

    } /* end if the image comes from the quadrant with reference spectrum */

    oneImage = oneImage->next;
  }

  /* now call VmSpCalPhot for spectro-photometric calibration */
  if (EXIT_FAILURE == VmSpCalPhot(stdImage, stdTable, stdRow, fitOrd))
  {
    pilMsgError(modName," Failure in deriving spectro-photom. calibration");
    return EXIT_FAILURE;
  }

  /* NOTE: ifu table is NOT updated with refL and refM*/


  deleteImage(ifuImage);
  deleteImage(ifuMask);
  deleteImage(stdImage);

  deleteObjectTable(objTable);
  deleteObjectObject(imaObjs);
  deleteTable(stdTable);
  deleteIfuTable(ifuTable);

  deletePixelRegion(pixRegion);

  return EXIT_SUCCESS;
}
