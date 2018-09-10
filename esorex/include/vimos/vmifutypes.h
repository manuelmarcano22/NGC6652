/* $Id: vmifutypes.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_IFUTYPES_H
#define VM_IFUTYPES_H

#include <pilmacros.h>

#include <vmimage.h>


PIL_BEGIN_DECLS

/* THESE FUNCTIONS & C. ARE STILL TO BE TESTED */

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  These are functions for thresholding and detecting adjacent
  pixels on ifu images (Eclipse- like)

------------------------------------------------------------------------------
*/


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  structure VimosPixelData

  Description:
  Linked list structure that contains coords and intensity of a pixel in 
  a region. 
  Created for thresholding an IFU image and looking for adjacent
  "IFU pixels" (VmIfuCalPhot function).
  Original from Eclipse: pixel_position

  Layout:
  int x;
  int y;
  float value;
  VimosPixelData *prev;
  VimosPixelData *next;

  Updates:
  15 Feb 00: Created (AZ)

-------------------------------------------------------------------------------
*/

typedef struct _VIMOS_PIXEL_DATA_ {
  int x;
  int y;
  float value;
  struct _VIMOS_PIXEL_DATA_ *prev;
  struct _VIMOS_PIXEL_DATA_ *next;
} VimosPixelData;



/*
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  structure VimosPixelList

  Description:
  Linked list structure that stores pixels from a region. 
  Created to be used when thresholding an IFU image and looking for adjacent
  "IFU pixels" (VmIfuCalPhot).
  Original from Eclipse: pixel_list

  Layout:
  int numPixInRegion
  float totIntensity
  VimosPixelData *pixelData
  VimosPixelList *prev;
  VimosPixelList *next;

  Updates:
  15 Feb 00: Created (AZ)

-------------------------------------------------------------------------------
*/

typedef struct _VIMOS_PIXEL_LIST_ {
  int numPixInRegion;
  float totIntensity;
  VimosPixelData *pixelData;
  struct _VIMOS_PIXEL_LIST_ *prev;
  struct _VIMOS_PIXEL_LIST_ *next;
} VimosPixelList;



/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  structure VimosPixelRegion

  Description:
  Linked list structure that contains regions of pixels.
  Created to be used when thresholding an IFU image and looking for adjacent
  "IFU pixels".
  More than one region detected above threshold is allowed: just in case two 
  stars fall in the same image, we don't want to lose the right (reference 
  fibre) one.
  Needed by VmIfuCalPhot.

  Layout:
  int numRegions
  VimosPixelList *pixelList
 
  Updates:
  15 Feb 00: Created (AZ)

-------------------------------------------------------------------------------
*/

typedef struct _VIMOS_PIXEL_REGION_ {
  int numRegions;
  VimosPixelList *pixelList;
  struct _VIMOS_PIXEL_REGION_ *prev;
  struct _VIMOS_PIXEL_REGION_ *next;
} VimosPixelRegion;


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  structure VimosPixelAcc

  Description:
  The following structure contains all information necessary for
  recursive search of consistent zones of pixels.  see the functions
  below to have an idea of how it is used. floodfill_from_pixel() is
  recursive get_pixelrank_if_valid() is the checking function is for
  recursive evals.
  Original from Eclipse: pixel_accumulator

  Layout:
  float      accIntensity            accumulated pixelvalues

  int        PixelsInCurrentBlock    number of pixels accumulated so far in 
				     the current consistent block of pixels

  int        imWidth                 image width, necessary to compute 
				     positions in the image 

  int        CurrentRank             rank of the current pixel being examined
				     mostly used to optimize speed 

  int        Total1Pixels            total number of pixels set to 1 in the 
				     currently examined pixel_map 

  int        RemainingPixels         number of pixels which still have to be 
				     examined 

  VimosBool    *ValidPixels          flag per pixel set to 1 in the pixel_map
				     a flag set to TRUE means the pixel has not
				     been associated to a zone which has been 
				     found. If it is set to FALSE, the pixel
				     has already been associated with a white 
				     zone 

  VimosUlong32    *pos               array contains all the positions of the 
				     white pixels in the input image. Positions
				     are not referenced by a couple (x,y) but 
				     by a single number POS = x + ImWidth*y 

  VimosImage *refImage               input image 

  Updates:
  15 Feb 00: Created (AZ)

-------------------------------------------------------------------------------
*/

typedef struct _VIMOS_PIXEL_ACC_
{
  float          accIntensity;
  int            PixelsInCurrentBlock;
  int            imWidth;
  int            CurrentRank;
  int            Total1Pixels;
  int            RemainingPixels;
  VimosBool     *ValidPixels;
  VimosUlong32  *pos;
  VimosImage    *refImage;
} VimosPixelAcc;


VimosPixelData *newPixelData();
void deletePixelData(VimosPixelData *aPixelData);

VimosPixelList *newPixelList();
void deletePixelList(VimosPixelList *aPixelList);

VimosPixelRegion *newPixelRegion();
void deletePixelRegion(VimosPixelRegion *aPixelRegion);






/*
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void  thresholdImage(VimosImage *ifuImage, float threshold, 
                     VimosImage *ifuMask, int *ngood)

   Description:
   Write a mask image with only pixels above threshold value.
   Needed for IFU.

   ORIGINAL FROM ECLIPSE: thresh_image_to_pixelmap

   Input: 
   VimosImage ifuImage
   Starting image to be thresholded

   float threshold
   Threshold value

   int ngood
   Number of pixels on the mask above the threshold value
   
   Output:
   VimosImage ifuMask
   Image with "bad" pixels (below threshold) set to "ZERO"

   Return Value:
   void
   
   Updates:
   14 Feb 00: Created (AZ)
-------------------------------------------------------------------------------
*/
void  thresholdImage(VimosImage *ifuImage, float threshold, 
                     VimosImage *ifuMask, int *ngood);



/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosPixelRegion *findRegionsOnPixelMap(VimosImage *refImage, 
                                          VimosImage *ifuMask,
                                          VimosPixelRegion *pixRegion,
                                          int ngoodpix)
  Description:
  Find regions of "connected" pixels and store them in a PixelRegion structure.
 
  ORIGINAL FROM ECLIPSE: find_centers_on_pixelmap()

  Input:
  VimosImage *refImage

  VimosImage *ifuMask
  Image with "bad" pixels (below threshold) set to "ZERO"

  VimosPixelRegion *pixRegion

  int ngoodpix

  Updates:
  15 Feb 00: Created (AZ)

-------------------------------------------------------------------------------
*/

VimosPixelRegion *findRegionsOnPixelMap(VimosImage *refImage, 
                                        VimosImage *ifuMask,
                                        VimosPixelRegion *pixRegion,
                                        int ngoodpix);

PIL_END_DECLS

#endif /* VM_IFUTYPES_H */
