/* $Id: vmifuimage.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <piltranslator.h>
#include <pilmessages.h>
#include <pilutils.h>

#include "vmimage.h"
#include "vmimageset.h"
#include "vmcube.h"
#include "vmtable.h"
#include "vmifutable.h"
#include "vmobjecttable.h"
#include "vmmath.h"
#include "vmifuimage.h"


#define XYSIZE 80


VimosImage *
VmIfu2DImage(VimosImageSet *imageSet, float wStart, float wEnd,
             float wLenStart, float wLenEnd, float wLenInc, int specLen)
{
  int objNum;
  int indexImage;
  int pixStart, pixEnd, pixRange, j;
  int quadrant, objL, objM, i, tmpSpecLen;
  float theFlux;

  char comment[80];
  char modName[] = "VmIfu2DImage";

  VimosFloatArray *tmpSpec, *spectrum = NULL;

  VimosObjectTable *objTable;
  VimosObjectObject *imaObjs;
  VimosIfuTable *ifuTable;
  VimosIfuQuad *theQuads;
  VimosIfuSlit *theIfuSlits;
  VimosSingleImage *oneImage;

  VimosImage *outImage;

  pilMsgInfo (modName, "Start computing 2D reconstructed Image");

  /* if full lambda range, set start and end to header values */
  if (wStart == 0.0)
    {
      wStart = wLenStart;
      wEnd = wLenEnd;
      tmpSpecLen = (wEnd - wStart) / wLenInc + 1;
    }

  /* first initialize to 0 the whole 2D image, so that we will not have 
     problems if less than 4 images are provided */

  outImage = newImageAndAlloc(XYSIZE, XYSIZE);
  indexImage = XYSIZE * XYSIZE;

  for (i=0; i<indexImage; i++) outImage->data[indexImage] = 0.0;

  /* start loop on images in the Image Set */
  oneImage = imageSet->images;
  while (oneImage)
    {
      /* take the object table and IFU table of this image */
      objTable = oneImage->objectTable;
      ifuTable = oneImage->ifuTable;

      /* take objects in the image object table */
      imaObjs = objTable->objs;
      objNum = 0;

      /* read which quadrant this image refers */
      readIntDescriptor(objTable->descs, pilTrnGetKeyword("Quadrant"),
			&quadrant, comment);

      theQuads = ifuTable->quads;
      
      while (theQuads)
	{
	  /* take from ifuTable only the slits in the quad */
	  /* corresponding to the image in use */
	  if (theQuads->quadNo == quadrant)
	    {
	      theIfuSlits = theQuads->ifuSlits;
	    }
	  theQuads = theQuads->next;
	}
      
      while (imaObjs)
	{
	  
	  tmpSpec = selectFiberForObject(theIfuSlits, imaObjs, 
					 (oneImage->theImage->data), specLen,
					 objNum, &objL, &objM);
	  
	  
	  /* clean previous allocations if they exist */
	  deleteFloatArray(spectrum);
	  /* here take just the range to be integrated */
	  /* if wStart and wEnd = 0, they are already set to */
	  /* wLenStart and wLenEnd */
	  pixStart = (wStart - wLenStart) / wLenInc;
	  pixEnd = (wEnd - wLenStart) / wLenInc;
	  pixRange = pixEnd - pixStart +1;
	  spectrum = newFloatArray(pixRange);

	  i = 0; 
	  j=0;
	  if (wStart == 0.0 && wEnd == 0.0)
	    {
	      for (i=pixStart; i<=pixEnd; i++)
		{
		  spectrum->data[j] = tmpSpec->data[i];
		  j++;
		}
	    }
	  else
	    {
	      spectrum=tmpSpec;
	    }
	  
	  /* compute the integrated flux */
	  theFlux = integrateSpectrum(spectrum, wLenInc);

	  /* CHECK IF THIS IS OK! */
	  /* assume fiber coordinates (L,M) start from 1, not 0 */
	  outImage->data[(objL-1) +XYSIZE*(objM-1)] = theFlux;
	  
	  /* take next object */
	  imaObjs = imaObjs->next;
	  objNum++;
	  
	} /* close loop on objects in this image */

      oneImage = oneImage->next; 
      
    } /* done the image Set */

  /* WARNING!!!! THIS FOR THE MOMENT!!!!!!!! */
  oneImage = imageSet->images;
  copyAllDescriptors(oneImage->theImage->descs, &(outImage->descs));
  
  return(outImage);

}


VimosCube *
VmIfu3DCube(VimosImageSet *imageSet, int specLen)
{

  int objNum;
  int j;
  int quadrant, objL, objM,i;
  VimosUlong32 ilong, indexCube;
  char comment[80];
  char modName[] = "VmIfu3DCube";

  VimosFloatArray *tmpSpec;

  VimosObjectTable *objTable;
  VimosObjectObject *imaObjs;
  VimosIfuTable *ifuTable;
  VimosIfuQuad *theQuads;
  VimosIfuSlit *theIfuSlits;
  VimosCube *outCube;
  VimosSingleImage *oneImage;

  pilMsgInfo (modName, "Start computing 3D Cube");

  /* first initialize to 0 the whole cube, so that we will not have problems
     if less than 4 images are provided */

  outCube = newCubeAndAlloc(XYSIZE, XYSIZE, specLen);
  indexCube = (XYSIZE * XYSIZE * specLen);
  for (ilong=0; ilong<indexCube; ilong++) outCube->data[indexCube] = 0.0;

  /* start loop on images in the Image Set */
  oneImage = imageSet->images;
  while (oneImage)
    {
      /* take the object table and IFU table of this image */
      objTable = oneImage->objectTable;
      ifuTable = oneImage->ifuTable;
      
      /* take objects in the image object table */
      imaObjs = objTable->objs;
      objNum = 0;
      
      /* read which quadrant this image refers */
      readIntDescriptor(objTable->descs, pilTrnGetKeyword("Quadrant"), 
			&quadrant, comment);
      
      theQuads = ifuTable->quads;
      
      while (theQuads)
	{
	  /* take from ifuTable only the slits in the quad */
	  /* corresponding to the image in use */
	  if (theQuads->quadNo == quadrant)
	    {
	      theIfuSlits = theQuads->ifuSlits;
	    }
	  theQuads = theQuads->next;
	}

      while (imaObjs)
	{
	  
	  tmpSpec = selectFiberForObject(theIfuSlits, imaObjs, 
					 (oneImage->theImage->data), specLen,
					 objNum, &objL, &objM);
	  
	  /* write the cube */
	  /* CHECK IF THIS IS OK! */
	  /* assume fiber coordinates (L,M) start from 1, not 0 */
	  i=0;
	  for (j=0; j<specLen; j++)
	    {
	      outCube->data[j*((objL-1) +XYSIZE*(objM-1))] = 
		tmpSpec->data[i];
	      i++;
	    }
	  
	  /* take next object */
	  imaObjs = imaObjs->next;
	  objNum++;
	  
	} /* close loop on objects in this image */

      oneImage = oneImage->next; 
      
    } /* done the image Set */
  
  return(outCube);

}
