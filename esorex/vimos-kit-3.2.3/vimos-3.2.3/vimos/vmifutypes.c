/* $Id: vmifutypes.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include "vmimage.h"
#include "vmifutypes.h"
#include "cpl.h"


/* THESE FUNCTIONS & C. ARE STILL TO BE TESTED */

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  These are functions for thresholding and detecting adjacent
  pixels on ifu images (Eclipse- like)

------------------------------------------------------------------------------
*/

/*
  Returns a pointer to a new Pixel Data
*/

VimosPixelData *newPixelData()
{
 VimosPixelData *newPData;
  /* allocate memory for Pixel Data */
  newPData = (VimosPixelData *) cpl_malloc(sizeof(VimosPixelData));

  /* error occured, exit */
  if (newPData == NULL)
    {
     abort();
    }

  /* fill up fields with default values */
  /* setting values to -1. to avoid always finding pixel (0,0) */
  /* IS IT OK? */
  newPData->x     = -1;
  newPData->y     = -1;
  newPData->value = -1.0;

  newPData->prev = newPData->next = NULL;
  /* return to caller */
  return(newPData);

}

/*
  Deletes all Pixel Data contained in the list
*/

void deletePixelData(VimosPixelData *aPixelData)
{
 VimosPixelData  *tmpPData;
 VimosPixelData  *nextPData;

  /* store start of the pixeldata list */
  tmpPData = aPixelData;

  /* traverse  list */
  while (tmpPData)
    {
     /* get address of next pixel data */
     nextPData = tmpPData->next;
     /* free current object */
     cpl_free(tmpPData);
     /* next one to process */
     tmpPData = nextPData;
    }
}


/*
  Returns a pointer to a new Pixel List
*/

VimosPixelList *newPixelList()
{
 VimosPixelList *newPList;
  /* allocate memory for Pixel List */
  newPList = (VimosPixelList *) cpl_malloc(sizeof(VimosPixelList));

  /* error occured, exit */
  if (newPList == NULL)
    {
     abort();
    }

  /* fill up fields with default values */
  newPList->numPixInRegion = 0;
  newPList->totIntensity   = 0.;
  newPList->pixelData      = NULL;

  newPList->prev = newPList->next = NULL;
  /* return to caller */
  return(newPList);

}

/*
  Deletes all Pixel List contained in the list
*/

void deletePixelList(VimosPixelList *aPixelList)
{
 VimosPixelList  *tmpPList;
 VimosPixelList  *nextPList;

  /* store start of the list */
  tmpPList = aPixelList;

  /* traverse  list */
  while (tmpPList)
    {
     /* delete the fibers in this slit */
     deletePixelData(tmpPList->pixelData);
     /* get address of next pixel list */
     nextPList = tmpPList->next;
     /* free current object */
     cpl_free(tmpPList);
     /* next one to process */
     tmpPList = nextPList;
    }
}



/*
  Returns a pointer to a new Pixel Region
*/

VimosPixelRegion *newPixelRegion()
{
  VimosPixelRegion *newPRegion;
  /* allocate memory for PixelRegion */
  newPRegion = (VimosPixelRegion *) cpl_malloc(sizeof(VimosPixelRegion));
  
  /* error occured, exit */
  if (newPRegion == NULL)
    {
     abort();
    }

  /* fill up fields with default values */
  newPRegion->numRegions = 0;
  newPRegion->pixelList = NULL;

  newPRegion->prev = newPRegion->next = NULL;
  /* return to caller */
  return(newPRegion);
  
}

/*
  Deletes all Pixel Regions contained in the list
*/

void deletePixelRegion(VimosPixelRegion *aPixelRegion)
{

  VimosPixelRegion *tmpPRegion, *nextPRegion;
  if (aPixelRegion == NULL)
    {
     return;
    }

  /* store start of the list */
  tmpPRegion = aPixelRegion;

  /* traverse  list */
  while (tmpPRegion)
    {

     deletePixelList(tmpPRegion->pixelList);
     /* get address of next pixel region */
     nextPRegion = tmpPRegion->next;
     /* free current object */
     cpl_free(tmpPRegion);
     /* next one to process */
     tmpPRegion = nextPRegion;
    }

}




/*
  Write a mask image with only pixels above threshold value.
*/

void  thresholdImage(VimosImage *ifuImage, float threshold, 
                     VimosImage *ifuMask, int *ngood)
{

#define ZERO  0.0
#define ONE  1.0

  int i;
  int imaSize;

  imaSize = ifuImage->xlen * ifuImage->ylen;
  *ngood = imaSize;

  /* loop on all input image pixels */
  for (i=0; i<imaSize; i++)
     {
      if (ifuImage->data[i] < threshold)
        {
         ifuMask->data[i] = ZERO;
         (*ngood)--;
        }
      else
        {
         ifuMask->data[i] = ONE;
        }
     }

#undef ZERO
#undef ONE

}

static int getPixelrankIfValid(VimosPixelAcc *pixAcc, int rank, int neighbor);

static int getPixelrankIfValid(VimosPixelAcc *pixAcc, int rank, int neighbor)
{
  /* Neighboring pixels are stored as: */
  /*      2                            */
  /*    1  3                           */
  /*      4                            */

  VimosUlong32 newpos;
  int i;
 
  /* This switch gets the neighbor position */
  switch (neighbor)
    {
     /* Left neighbor */
     case 1:
     newpos = pixAcc->pos[rank] - 1;
     break;

     /* Up neighbor */
     case 2:
     newpos = pixAcc->pos[rank] + pixAcc->imWidth;
     break;
 
     /* Right neighbor */
     case 3:
     newpos = pixAcc->pos[rank] + 1;
     break;
 
     /* Down neighbor */
     case 4:
     newpos = pixAcc->pos[rank] - pixAcc->imWidth;
     break;

     /* Invalid neighbor */
     default:
     puts("internal error: no such neighbor in getPixelrankIfValid");
     return 0;
    }

  /* Find out if the neighbor position is in a list, and if it is, */
  /* find if it has bot been tagged in another region yet          */

  for (i=pixAcc->CurrentRank ; i<pixAcc->Total1Pixels ; i++)
     {
      if ((pixAcc->pos[i] == newpos) && (pixAcc->ValidPixels[i] == TRUE))
          /* Found neighbor in the list: return its rank */
          return(i) ;
     }

  /* the neighbor is not a valid pixel: return 0 */
  return 0 ;

}


static void floodfillFromPixel(VimosPixelAcc *pixAcc, int rank,
                     VimosPixelList *pixList, VimosPixelData *lastPixData);


static void floodfillFromPixel(VimosPixelAcc *pixAcc, int rank,
                      VimosPixelList *pixList, VimosPixelData *lastPixData)
{
  int vrank;
  int posx, posy;
  VimosPixelData *pixData;

  /* At least the given pixel is white, accumulate its X and Y */
  /* coordinates and tag it as taken (FALSE)                   */
  posx = pixAcc->pos[rank] % pixAcc->imWidth;
  posy = pixAcc->pos[rank] / pixAcc->imWidth;

  pixAcc->accIntensity +=
                pixAcc->refImage->data[posx+posy*pixAcc->imWidth];
  pixAcc->ValidPixels[rank] = FALSE;
  pixAcc->RemainingPixels --;
  pixAcc->PixelsInCurrentBlock++;

  /* store coordinates and intensity of each pixel in this region in a */
  /* pixel list                                                        */
  pixData = newPixelData();
  pixData->x = posx;
  pixData->y = posy;
  pixData->value = pixAcc->refImage->data[posx+posy*pixAcc->imWidth];

  /* start linking the pixel to a list */
  /* link pixel to list */
  if (lastPixData == NULL) 
    {
     /* first pixel we do */
     pixList->pixelData = pixData;
    }
  else
    {
     /* not first, add to tail */
     lastPixData->next = pixData;
     pixData->prev = lastPixData;
    }
  /* store tail of linked list */
  lastPixData = pixData;

  /* Go recursive with 4 neighbors */
  if ((vrank = getPixelrankIfValid(pixAcc, rank, 1)) != 0)
      floodfillFromPixel(pixAcc, vrank, pixList, lastPixData);
 
  if ((vrank = getPixelrankIfValid(pixAcc, rank, 2)) != 0)
      floodfillFromPixel(pixAcc, vrank, pixList, lastPixData);
 
  if ((vrank = getPixelrankIfValid(pixAcc, rank, 3)) != 0)
      floodfillFromPixel(pixAcc, vrank, pixList, lastPixData);
 
  if ((vrank = getPixelrankIfValid(pixAcc, rank, 4)) != 0)
      floodfillFromPixel(pixAcc, vrank, pixList, lastPixData);

}

/*
  Find regions of "connected" pixels and store them in a PixelRegion structure.
*/

VimosPixelRegion *findRegionsOnPixelMap (VimosImage *refImage, 
                                         VimosImage *ifuMask, 
                                         VimosPixelRegion *pixRegion,
                                         int ngoodpix)
{

#define ZERO  0.0
#define ONE  1.0

  int i;
  int foundRegions;
  int numMapPix, nbpix;

  VimosPixelAcc pixAcc;
  VimosPixelList *lastPixList, *pixList;
  VimosPixelData *lastPixData;

  /* Initialize number of found white blobs to zero */
  foundRegions = 0;

  /* get total number of pixels in the ifu mask */
  numMapPix = ifuMask->xlen * ifuMask->ylen;

  /* if more than 60% of the input pixels are GOOD, do nothing ... WHY????*/
  if (ngoodpix > (int)(0.6 * (double)numMapPix))
    {
     puts("findRegionsOnPixelMap: more than 60% of good pixel. Exiting");
     abort();
    }

  /* Allocate memory: we will store a pixel position and a boolean flag */
  /* per valid pixel in the input binary map.                           */
 
  pixAcc.pos = (VimosUlong32*)cpl_calloc(ngoodpix, sizeof(VimosUlong32));
  pixAcc.ValidPixels = (VimosBool*)cpl_malloc(ngoodpix * sizeof(VimosBool));
 
  /* Set all flags to TRUE: no pixel has been placed in any region yet */
  for (i=0 ; i<ngoodpix ; i++) pixAcc.ValidPixels[i] = TRUE;

  /* Store all valid pixel positions in an array */
  nbpix = 0;

  for (i=0 ; i<numMapPix ; i++)
     {
      if (ifuMask->data[i] == ONE) pixAcc.pos[nbpix++] = i;
     }

  /* Initialize all before entering recursivity */
 
  pixAcc.CurrentRank = 0;
  pixAcc.imWidth = ifuMask->xlen;
  pixAcc.Total1Pixels = ngoodpix;
  pixAcc.RemainingPixels = ngoodpix;

  /*  pixAcc.accXpos = 0.0; */
  /*  pixAcc.accYpos = 0.0; */

  pixAcc.accIntensity = 0.;
  pixAcc.PixelsInCurrentBlock = 0;
  pixAcc.refImage = refImage ;

  i=0;
  lastPixList = NULL;

  while (pixAcc.RemainingPixels)
    {
     /* Find first valid pixel in the list */
     while ( (pixAcc.ValidPixels[i] == FALSE) && 
             (i<pixAcc.Total1Pixels))
       {
        i++;
       }
 
     /* If we reached the table size, no more pixels to compute */
     if (i == pixAcc.Total1Pixels) break ;

     pixAcc.CurrentRank = i;

     lastPixData = NULL;
     pixList = newPixelList();

     /* Go recursive on all white blobs */
     floodfillFromPixel(&pixAcc, i, pixList, lastPixData);

     /* One more found blob */
     foundRegions++;

     pixList->totIntensity = pixAcc.accIntensity;
     pixList->numPixInRegion = pixAcc.PixelsInCurrentBlock;

     /* link pixelList to list */
     if (lastPixList == NULL) 
       {
        /* first pixelList we do */
        pixRegion->pixelList = pixList;
       }
     else
       {
        /* not first, add to tail */
        lastPixList->next = pixList;
        pixList->prev = lastPixList;
       }
     /* store tail of linked list */
     lastPixList = pixList;
 
     /* Reset accumulators and loop */
     pixAcc.accIntensity = 0.;
     pixAcc.PixelsInCurrentBlock = 0;
    }

  pixRegion->numRegions = foundRegions;

  cpl_free(pixAcc.ValidPixels) ;
  cpl_free(pixAcc.pos) ;

#undef ZERO
#undef ONE

  return(pixRegion);
}
