/* $Id: vmimageset.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <pilmemory.h>

#include "vmimage.h"
#include "vmccdtable.h"
#include "vmifutable.h"
#include "vmextractiontable.h"
#include "vmwindowtable.h"
#include "vmobjecttable.h"
#include "vmimageset.h"
#include "cpl.h"

#define MAXFILES 50


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosSingleImage *newSingleImage()

 
  Description: 
  Returns a pointer to a new VimosSingleImage structure.
 
  Input:
  void   
 
  Return Value (succes):
  Pointer to a newly allocated VimosSingleImage structure.

  Updates:
  03 Feb 00: Created (AZ)

---------------------------------------------------------------------
*/

VimosSingleImage *
newSingleImage()
{
 VimosSingleImage *theSingleImage;
  /* allocate memory for one image */
  theSingleImage = (VimosSingleImage *) cpl_malloc(sizeof(VimosSingleImage));

  /* error occured, exit */
  if (theSingleImage == NULL) {
    abort();
  }

  /* fill up fields with default values */

  theSingleImage->theImage = NULL; 
  theSingleImage->windowTable = NULL;
  theSingleImage->ifuTable = NULL;
  theSingleImage->objectTable = NULL;
  theSingleImage->extractionTable = NULL;
  theSingleImage->ccdTable = NULL;
  theSingleImage->sphotStdTable = NULL;

  theSingleImage->imageName = "";

  theSingleImage->next = theSingleImage->prev = NULL;

  return(theSingleImage);
}


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void deleteSingleImage(VimosSingleImage *oneImage)

  Description:
  Deletes all VimosSingleImages contained in the linked list VimosSingleImage.

  Input:
  VimosSingleImage *oneImage
  Pointer to VimosSingleImage list to be deleted

  Return Value:
  void
   
  Updates:
  03 Feb 00: Created (AZ)

---------------------------------------------------------------------
*/

void
deleteSingleImage(VimosSingleImage *oneImage)
{
  VimosSingleImage *tmpSingleImage;
  VimosSingleImage *nextSingleImage;

  /* store start of the list */
  tmpSingleImage = oneImage;

  /* traverse  list */
  while (tmpSingleImage)
    {
     /* delete the image */
     deleteImage(tmpSingleImage->theImage);
     deleteCcdTable(tmpSingleImage->ccdTable);

     /*WARNING!!!! COMMENTED, WE NEED TO DEFINE THIS FUNCTION*/
     /* deleteSPhotTable(theSingleImage->sphotStdTable); */

     if (tmpSingleImage->objectTable) 
       deleteObjectTable(tmpSingleImage->objectTable);
     if (tmpSingleImage->ifuTable)
       deleteIfuTable(tmpSingleImage->ifuTable);
     if (tmpSingleImage->windowTable) 
       deleteWindowTable(tmpSingleImage->windowTable);
     if (tmpSingleImage->extractionTable)
       deleteExtractionTable(tmpSingleImage->extractionTable);

     /* get address of next SingleImage */
     nextSingleImage = tmpSingleImage->next;
     /* free current SingleImage */
     cpl_free(tmpSingleImage);
     /* next one to process */
     tmpSingleImage = nextSingleImage;

     /* image name? */
    }
}


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosImageSet *newImageSet()

 
  Description: 
  Returns a pointer to a new VimosImageSet structure.
 
  Input:
  void   
 
  Return Value (succes):
  Pointer to a newly allocated VimosImageSet structure.

  Updates:
  03 Feb 00: Created (AZ)

---------------------------------------------------------------------
*/

VimosImageSet *
newImageSet()
{
  VimosImageSet *theImageSet;
  /* allocate new VimosImageSet */
  theImageSet = (VimosImageSet *) cpl_malloc(sizeof(VimosImageSet));

  /* if error: exit */
  if (theImageSet == NULL) {
    abort();
  }

  theImageSet->images = NULL;
  theImageSet->nImas = 0;

  return(theImageSet);
}


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  void deleteImageSet(VimosImageSet *imageSet)

  Description:
  Deletes the linked list VimosImageSet.

  Input:
  VimosImageSet *imageSet
  Pointer to VimosImageSet list to be deleted

  Return Value:
  void
   
  Updates:
  03 Feb 00: Created (AZ)

---------------------------------------------------------------------
*/

void
deleteImageSet(VimosImageSet *theImageSet)
{
 if (theImageSet == NULL) {
    return;
  }

 deleteSingleImage(theImageSet->images);

}


/*
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   void storeImagesAndTables(VimosImageSet *imageSet, char *imageList, 
                             char extname[],
                             int ccdT, int sphT, int ifuT, int extT, 
			     int winT, int objT);
  
   Description:
   read an input list of images and creates an imageSet.
   For each image, this function can store in the imageSet the IFU, 
   Extraction, Window and Object Table, depending on the parameters 
   ifuT, extT, winT, objT.
 
   Input: 
   VimosImageSet imageSet
   Pointer to an imageSet

   char imageList
   Input file containing the list of image names to be read and stored

   int ccdT, int sphT, int ifuT, int extT, int winT, int objT
   For each image to be stored, these parameters say if the CCD, 
   spectrophotometric, IFU, Extraction, Window and Object Tables must be 
   stored as well. 
   Set to 1 to store the table. Set to 0 to skip the table.

   Return Value:
   void
   
   Updates:
   01 Mar 00: Created (AZ)

-------------------------------------------------------------------------------
*/

void
storeImagesAndTables(VimosImageSet *imageSet, char *imageList, char extname[],
                     int ccdT, int sphT, int ifuT, int extT, int winT,
                     int objT)
{
  char inImaName[80];
  int nImaInSet, numImages, nChar;

  VimosImage *tmpImage;
  VimosSingleImage *lastOneImage, *oneImage;
  VimosTable *tmpCcdTable;
  VimosIfuTable *tmpIfuTab;
  VimosExtractionTable *tmpExtTab;
  VimosWindowTable *tmpWinTab;
  VimosObjectTable *tmpObjTab;

  FILE *ifp = NULL;

  
  if ((ifp = fopen(imageList,"r")) == NULL)
    {
     printf("Can't open file %s. \n", imageList);
     abort();
    }

  lastOneImage = NULL;
  nImaInSet = 0;

  while(fgets(inImaName,MAXFILES,ifp) != NULL) numImages++;
 
  if(ifp) rewind(ifp);

  while(fgets(inImaName,MAXFILES,ifp) != NULL)
  {
    nChar = strlen(inImaName);
    if(inImaName[nChar-1] == '\n') inImaName[nChar-1] = '\0';
    if(inImaName[nChar-2] == ' ') inImaName[nChar-2] = '\0';

    oneImage = newSingleImage();
    tmpImage = openOldFitsFile(inImaName, 0, 1);
    /*    oneImage->theImage = openOldFitsFile(inImaName, 0, 1); */
    oneImage->theImage = openFitsImageExtension(tmpImage->fptr, extname);
    oneImage->imageName = inImaName;

    /* read tables according to the selection */

    if (ccdT == 1)
      {
       tmpCcdTable = newCcdTable();
       readFitsCcdTable(tmpCcdTable, oneImage->theImage->fptr);
       oneImage->ccdTable = tmpCcdTable;
      }
    /* WARNING!!! WE NEED TO SET THESE READ FUNCTIONS FOR SPhotTable!!! */
    /*
    if (sphT == 1)
      {
       tmpSphotStdTable = newSPhotTable();
       readFitsSPhotTable(tmpSphotStdTable, oneImage->theImage->fptr);
       oneImage->sphotStdTable = tmpSphotStdTable;
      }
    */
    if (ifuT == 1)
      {
       tmpIfuTab = newIfuTable();
       readFitsIfuTable(tmpIfuTab, oneImage->theImage->fptr);
       oneImage->ifuTable = tmpIfuTab;
      }
    if (extT == 1)
      {
       tmpExtTab = newExtractionTable();
       readFitsExtractionTable(tmpExtTab, oneImage->theImage->fptr);
       oneImage->extractionTable = tmpExtTab;
      }
    if (winT == 1)
      {
       tmpWinTab = newWindowTable();
       readFitsWindowTable(tmpWinTab, oneImage->theImage->fptr);
       oneImage->windowTable = tmpWinTab;
      }
    if (objT == 1)
      {
       tmpObjTab = newObjectTable();
       readFitsObjectTable(tmpObjTab, oneImage->theImage->fptr);
       oneImage->objectTable = tmpObjTab;
      }

    /* link oneImage to list */
    if (lastOneImage == NULL)
      {
       /* first image we do */
       imageSet->images = oneImage;
      }
    else
      {
       /* not first, add to tail */
       lastOneImage->next = oneImage;
       oneImage->prev = lastOneImage;
      }
    /* store tail of linked list */
    lastOneImage = oneImage;

    nImaInSet++;

  }

  if (ifp) fclose(ifp);

  deleteIfuTable(tmpIfuTab);
  deleteCcdTable(tmpCcdTable);

  /* WARNING!!! WE NEED TO SET THESE READ FUNCTIONS FOR SPhotTable!!! */
  /*
  deleteSPhotTable(tmpSphotStdTable);
  */
  deleteExtractionTable(tmpExtTab);
  deleteWindowTable(tmpWinTab);
  deleteObjectTable(tmpObjTab);
  deleteSingleImage(lastOneImage);
  deleteImage(tmpImage);

  return;
}


/*
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  void closeFitsImageSet(VimosImageSet *imageSet, int flag);
  
   Description:
   Closes a set if IFU images
 
   Input: 
   VimosImageSet imageSet
   Pointer to an imageSet

   int flag

   Return Value:
   void
   
   Updates:
   01 Mar 00: Created (AZ)

-------------------------------------------------------------------------------
*/

void
closeFitsImageSet(VimosImageSet *imageSet, int flag)
{
  VimosSingleImage *allImages;

  allImages = imageSet->images;

    while (allImages)
      {
       closeFitsImage(allImages->theImage, flag);
       allImages = allImages->next;
      }

    return;
}
