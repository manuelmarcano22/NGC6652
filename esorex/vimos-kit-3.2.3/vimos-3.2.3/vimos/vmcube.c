/* $Id: vmcube.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <pilmemory.h>

#include "vmtypes.h"
#include "vmimage.h"
#include "vmcube.h"
#include "cpl.h"


/*
  Create a VimosCube, not allocate space for the data
*/

VimosCube *newCube(int xlen, int ylen, int zlen, float *data) 
{
  VimosCube *theCube;
  
  theCube = (VimosCube *) cpl_malloc(sizeof(VimosCube));

  /* check if space was allocated */
  if (theCube == NULL) 
    {
       abort();
    }
    
  /* a VimosCube is only an interface to a Fits Image. */
  /* The pointer should point to the data of the Fits Image */
 
  theCube->data = data;  
  theCube->xlen = xlen;
  theCube->ylen = ylen;
  theCube->zlen = zlen;
  theCube->descs = NULL;
  theCube->fptr = NULL;
  
  return(theCube);
}


/*
  Create a VimosCube, allocate space for the data
*/

VimosCube *newCubeAndAlloc(int xlen, int ylen, int zlen)
{
  VimosCube *tmpCube;
  
  tmpCube = newCube(xlen, ylen, zlen, NULL);
    
  tmpCube->data = (float *) cpl_malloc(xlen*ylen*zlen*sizeof(float));
  if (tmpCube->data == NULL)
    {
     /* cleanup */
     deleteCube(tmpCube);
     abort();
    }
  
  return(tmpCube);
  
}


/*
  Delete a VimosCube 
*/

void deleteCube(VimosCube *Cube) 
{
  if (Cube != NULL) {
    deleteAllDescriptors(Cube->descs);
    if (Cube->data != NULL) {
      cpl_free(Cube->data);
    }
    cpl_free(Cube);
  }
}

void deleteCubeAndAlloc(VimosCube *Cube)
{
  if (Cube != NULL) {
    if (Cube->data != NULL) {
      cpl_free(Cube->data);
      Cube->data = NULL;
    }
  }
  deleteCube(Cube);
}



/*
  Find (L,M) coordinates of the fiber corresponding to an Object
  in an objectTable
*/

VimosFloatArray *selectFiberForObject(VimosIfuSlit *ifuSlits, 
                                      VimosObjectObject *object,
                                      float *theData, int specLen, 
                                      int objNum,int *theL, int *theM)
{
  int i;

  VimosFloatArray *theSpec = NULL;
  VimosIfuFiber *theIfuFibers;


  /* loop on IFU slits for this quadrant to find L,M of each object */
  while (ifuSlits)
       {
        /* if this is the "right" slit take fibers */
        if (ifuSlits->ifuSlitNo == object->IFUslitNo)
          {
           theIfuFibers = ifuSlits->fibers;

           /* loop on fibers of this slit until find the one */
           /* corresponding to the object */
           while (theIfuFibers)
             {
              if (theIfuFibers->fibNo == object->IFUfibNo)
                {
               	 /* clean previous allocation */
	         deleteFloatArray(theSpec);

		 /* take (L,M) coordinates for the fiber */
		 *theL = theIfuFibers->fiberL;
		 *theM = theIfuFibers->fiberM;

		 /* create working array for this object spectrum */
		 theSpec = newFloatArray(specLen);
                 for (i = 0; i < specLen; i++ )
                    {
	             theSpec->data[i] = theData[i + objNum*specLen];
	            }

		 /* here compute flux or whatever else */
                }
              theIfuFibers = theIfuFibers->next;
             }
          }

        ifuSlits = ifuSlits->next;

       } /* close loop on IFU slits for this object */

  return(theSpec);
}


/*
  Open a new Fits VimosCube
*/

VimosBool  openNewFitsCube(char *cubeName, VimosCube *cube)
{
  int naxis;
  long npix[3];
  int status, bitpix;

  npix[0] = cube->xlen;
  npix[1] = cube->ylen;
  npix[2] = cube->zlen;

  naxis = 3;
  bitpix = -32;
  
  status = 0;

  fits_create_file (&cube->fptr, cubeName, &status);
  fits_create_img(cube->fptr,  bitpix, naxis, npix, &status); 
  
  if (status)
    {
     return(VM_FALSE);
    }
  else
    {
     return(VM_TRUE);
    }
    
}


/*
  Close a Fits VimosCube 
*/

VimosBool closeFitsCube(VimosCube *cube, int flag) 
{
  int status;
  long fpixel, npixels;

  status = 0;

  npixels  = cube->xlen * cube->ylen * cube->zlen;         
  fpixel   = 1;

  if (flag != 0 )
    {
     if (fits_write_img(cube->fptr, TFLOAT, fpixel, npixels, cube->data,
                        &status) )
     return(VM_FALSE);      
    }

  status = 0;
  if ( fits_close_file(cube->fptr, &status) ) return(VM_FALSE);

  return(VM_TRUE);
}
