/* $Id: vmccdtable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <string.h>

#include <fitsio.h>

#include <pilmemory.h>
#include <pilstrutils.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmtable.h"
#include "vmmath.h"
#include "vmccdtable.h"
#include "cpl.h"


#define MAX_BAD_PIXEL_FRACTION (0.15)


static VimosTable *
createCcdTable()
{

    VimosTable      *ccdTable;
    VimosColumn     *currentColumn;

    if ((ccdTable = newCcdTable())) {

        /*
         * Columns initialization
         */

        ccdTable->cols = newColumn();           /* Column X */
        currentColumn = ccdTable->cols;
        strcpy(currentColumn->colName,"X");
        currentColumn->colType = VM_INT;

        currentColumn->next = newColumn();      /* Column Y */
        currentColumn = currentColumn->next;
        strcpy(currentColumn->colName,"Y");
        currentColumn->colType = VM_INT;

        ccdTable->numColumns = 2;
    }

    return(ccdTable);   /* ccdTable == NULL in case of failures above */

}


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosTable *newCcdTable()

   Description:
   Allocate a new Ccd Table. This is a VimosTable with name VM_CCD (which
   has value "CCD").
   Only a Descriptor "TABLE" is present, with value "CCD".
   
   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosTable

   Return Value (error):
   function calls exit(-1)
   
   Updates:
   11 Feb 99: Created (TAO)

-------------------------------------------------------------------------------
*/

VimosTable *
newCcdTable()
{

  VimosTable *newTab;
  
  /* allocate new VimosTable */
  newTab = newTable();

  /* check if space was allocated */
  if (newTab == NULL) {
    cpl_msg_error("newCcdTable","The function newTable has returned NULL");
    return(NULL);
  }
  
  /* copy "CCD" into name of table */
  strcpy(newTab->name, VM_CCD);
  newTab->descs = newStringDescriptor(pilTrnGetKeyword("Table"), VM_CCD,
                                      pilTrnGetComment("Table"));

  if (newTab->descs == NULL) {
    cpl_free(newTab);
    cpl_msg_error("newCcdTable","The function newStringDescriptor has "
                "returned NULL");
    return(NULL);
  }

  /* should create all the descriptors and columns here ?*/
  
  return(newTab);
  
}


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteCcdTable(VimosTable *cTable)
   
   Description:
   Delete a Ccd Table. This is just an esthetic wrapper for
   deleteTable(cTable)
   
   Input: 
   VimosTable *cTable
   Pointer of Ccd Table to be deleted
   
   Return Value:
   void
   
   Updates:
   11 Feb 99: Created (TAO)

-------------------------------------------------------------------------------
*/

void
deleteCcdTable(VimosTable *ccdTable)
{
  /* just a wrapper for deleteTabel() */
  deleteTable(ccdTable);
}

VimosTable *badPixelImage2CcdTable(VimosImage *badPixelImage)
{
  int              i,npix;
  int              x,y;
  int              nbadpix = 0;

  VimosTable      *ccdTable;
  VimosColumn     *currentColumn;


  if ((ccdTable = createCcdTable())) {
      const char *hint = cpl_strdup(pilTrnGetKeyword("Table"));

      /*
       * The table is created from the bad pixel map, which
       * is assumed to already have a correct and complete
       * list of header keywords.
       */

      vimosDscCopy(&ccdTable->descs, badPixelImage->descs,
                   "[A-Z].*", hint);

      cpl_free((void *)hint);

    npix = badPixelImage->xlen * badPixelImage->ylen;

    for (i = 0; i < npix; i++) nbadpix += badPixelImage->data[i];

    currentColumn = ccdTable->cols;
    currentColumn->colValue->iArray = (int *)cpl_malloc(nbadpix * sizeof(int));
    currentColumn->len = nbadpix;

    currentColumn = currentColumn->next;
    currentColumn->colValue->iArray = (int *)cpl_malloc(nbadpix * sizeof(int));
    currentColumn->len = nbadpix;

    for (x = 0; (x < badPixelImage->xlen) && nbadpix; x++) {
      for (y = 0; (y < badPixelImage->ylen) && nbadpix; y++) {
        if (badPixelImage->data[x + y * badPixelImage->xlen] > .5) {
          nbadpix--;
          currentColumn = ccdTable->cols;
          currentColumn->colValue->iArray[nbadpix] = x + 1;
          currentColumn = currentColumn->next;
          currentColumn->colValue->iArray[nbadpix] = y + 1;
        }
      }
    }
  }
  return(ccdTable);   /* ccdTable == NULL in case of failures above */
}


int cleanBadPixels(VimosImage *image, VimosTable *ccdTable, int correctLevel)
{
  char         modName[] = "cleanBadPixels";

  int         *isBadPix;
  int          i, j, k, d;
  int          xlen, ylen, totPix;
  int          nBadPixels = 0;
  int          sign, foundFirst;
  VimosColumn *currColumn;
  VimosColumn *colX;
  VimosColumn *colY;
  VimosColumn *colLevel;
  int         *xValue = NULL;
  int         *yValue = NULL;
  int         *level  = NULL;
  float        save = 0.;
  double       sumd;
  int          cx, cy;
  int          nPairs;
  float        estimate[4];
  int          sx[] = {0, 1, 1, 1};
  int          sy[] = {1,-1, 0, 1};
  int          searchHorizon = 100;
  int          percent = MAX_BAD_PIXEL_FRACTION * 100;

  char         instMode[80];
  int          imageMode, ccdTableMode;


  if (image == NULL) {
    cpl_msg_error(modName, "Null input image: aborting bad pixels cleaning");
    return(EXIT_FAILURE);
  }
  if (ccdTable == NULL) {
    cpl_msg_warning(modName, "No CCD Table: no bad pixel cleaning");
    return(EXIT_SUCCESS);
  }

  readStringDescriptor(image->descs, pilTrnGetKeyword("InstrumentMode"),
                       instMode, NULL);
  imageMode = (instMode[1] == 'M' ? 1 : 0);  /* 0 = MOS, 1 = IMG */

  /*
   * Check that the given CCD_TABLE is the right one, i.e., the one with
   * the right instrument mode. If the instrument mode is missing, this
   * is a pure cosmic ray events table, that is assumed correct anyway.
   */

  if (VM_TRUE == readStringDescriptor(ccdTable->descs, 
                 pilTrnGetKeyword("InstrumentMode"), instMode, NULL)) {
    ccdTableMode = (instMode[1] == 'M' ? 1 : 0);  /* 0 = MOS, 1 = IMG */
  }
  else
    ccdTableMode = imageMode;

  if (ccdTableMode != imageMode) {
    cpl_msg_error(modName, "The CCD table doesn't match the image to be cleaned");
    return(EXIT_FAILURE);
  }

 /* 
  *  Read bad pixels listed in CCD table 
  */

  if (ccdTable->cols == NULL)  {
    nBadPixels = 0;
  }
  else {
    currColumn = ccdTable->cols;
    nBadPixels = currColumn->len;
  }

  if (nBadPixels) {
    xlen = image->xlen;
    ylen = image->ylen;
    totPix = xlen * ylen;
    if (((float) nBadPixels) / ((float) totPix) < MAX_BAD_PIXEL_FRACTION) {
      isBadPix = (int *) cpl_malloc(totPix * sizeof(int));
      for (i = 0; i < totPix; i++) isBadPix[i] = 0;
    }
    else {
      cpl_msg_warning(modName, 
      "Too many bad pixels (> %d%%): skip bad pixel correction", percent);
      return(EXIT_FAILURE);
    }
  }
  else {
    cpl_msg_info(modName, "No bad pixels to correct");
    return(EXIT_SUCCESS);
  }
  if ((colX = findColumn(currColumn,"X")) 
                                    && (colY = findColumn(currColumn,"Y"))) {
    xValue = colX->colValue->iArray;
    yValue = colY->colValue->iArray;
    for (i = 0; i < nBadPixels; i++) {
        isBadPix[xValue[i]-1 + (yValue[i]-1) * xlen] = 1;
    }
  }
  else {
    cpl_msg_error(modName, "Missing coordinates columns in CCD table");
    cpl_free(isBadPix);
    return(EXIT_FAILURE);
  }

  if (correctLevel) {
    if ((colLevel = findColumn(currColumn, "level"))) {
      level = colLevel->colValue->iArray;
    }
    else {
      cpl_msg_error(modName, "Missing level column in CCD Table");
      cpl_free(isBadPix);
      return(EXIT_FAILURE);
    }
  }

  for (i = 0; i < nBadPixels; i++) { 
    if (correctLevel && correctLevel != level[i]) continue;
   /*
    *  Search for the closest good pixel along the 4 fundamental directions
    *  (in both senses):
    *                            \ | /
    *                             \|/
    *                           --- ---
    *                             /|\
    *                            / | \
    *
    *  Then collect pairs of values to interpolate linearly.
    */
    nPairs = 0;
    for (j = 0; j < 4; j++) {
      estimate[nPairs] = 0.;  /* Pairs interpolation results are stored here */
      sumd = 0.;
      foundFirst = 0;
      for (k = 0; k < 2; k++) {
        sign = 2 * k - 1;
        d = 0;
        cx = xValue[i]-1;
        cy = yValue[i]-1;
        do {
          cx += sign * sx[j];
          cy += sign * sy[j];
          if (cx < 0 || cx >= xlen || cy < 0 || cy >= ylen) break;
          d++;
        } while (isBadPix[cx + cy * xlen] && d < searchHorizon);
        if (cx >= 0 && cx < xlen && cy >= 0 && cy < ylen &&
            d < searchHorizon) {
         /*
          *  In this block is cripted the linear interpolation...
          */
          save = image->data[cx + cy * xlen];
          estimate[nPairs] += save / d;
          sumd += 1. / (double) d;
          if (k) {
            estimate[nPairs] /= sumd;
            nPairs++;
          }
          else {
            foundFirst = 1;
          }
        }
        else {
         /* 
          * Image borders were crossed, incomplete pair of values
          */
          if (k) {
            if (foundFirst) {
              estimate[nPairs] = save;
              nPairs++;
            }
          }
        }
      }
    }
   /*
    * Replace pixel value of the input image, corresponding to
    * the current bad pixel, with the median of the estimates
    * resulted from the 4 linear interpolations.
    */
    if (nPairs > 2) {
      image->data[(xValue[i]-1) + (yValue[i]-1) * xlen] = 
                                               medianWirth(estimate, nPairs);
    }
    else if (nPairs == 2) {
      image->data[(xValue[i]-1) + (yValue[i]-1) * xlen] = 
                                            (estimate[0] + estimate[1]) / 2.;
    }
    else if (nPairs == 1) {
      image->data[(xValue[i]-1) + (yValue[i]-1) * xlen] = estimate[0];
    }
    else {
      cpl_msg_warning(modName, "Cannot correct bad pixel %d,%d\n",
                                                       xValue[i], yValue[i]);
    }
  }
  cpl_free(isBadPix);
  return(EXIT_SUCCESS);
}

int shiftCcdTableCoords(VimosTable *ccdTable, int vertical, int shift)
{
  VimosColumn *coord;
  int          numBad;
  int          i;
  int          exitStatus = EXIT_FAILURE;

  if (ccdTable) {
    if (ccdTable->numColumns > 0) {
      numBad = ccdTable->cols->len;
      coord = ccdTable->cols;
      if (vertical) coord = coord->next;
      for (i = 0; i < numBad; i++) coord->colValue->iArray[i] += shift;
      exitStatus = EXIT_SUCCESS;
    }
  }
  return(exitStatus);
}


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool readFitsCcdTable(VimosTable *ccdTable, fitsfile *fptr)

   Description:
   Read a FITS Table into a CCD Table. This requires the FITS file to 
   be already open, and does not close it at the end.
   
   Input:
   VimosCCDTable *ccdTable  
   Pointer to a CCD Table to copy the FITS table into.

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool (true if the read is succesful, false otherwise)
   
   Updates:
   15 Dec 99: Created (BG)
   02 Feb 00: Changed the handling of the table, from a standalone FITS file,
              to an extension of the primary FITS file (MS)

-------------------------------------------------------------------------------
*/

VimosBool
readFitsCcdTable(VimosTable *ccdTable, fitsfile *fptr)
{

    int status;

    /* validate input */
    if (ccdTable == NULL) {
        cpl_msg_error("readFitsCcdTable","NULL input table");
        return(VM_FALSE);
    }
    if (fptr == NULL) {
        cpl_msg_error("readFitsCcdTable","NULL pointer to fitsfile");
        return (VM_FALSE);
    }

    if (strcmp(ccdTable->name, VM_CCD)) {
        cpl_msg_error("readFitsCcdTable","Invalid input table");
        return(VM_FALSE);
    }
  
    status = 0;

    /* open Table */
    if (fits_movnam_hdu(fptr, BINARY_TBL, "CCD", 0, &status) ){
        cpl_msg_error("readFitsCcdTable","The function fits_movnam_hdu has "
                    "returned an error (code %d)", status);
        return(VM_FALSE);
    }

    ccdTable->fptr = fptr;

    if (!readFitsTable(ccdTable, ccdTable->fptr)) {
        cpl_msg_info("readFitsCcdTable", "Error in reading the FITS file");
        return VM_FALSE;
    }

    return(VM_TRUE);
}  


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool writeFitsCcdTable(VimosTable *ccdTable, fitsfile *fptr)

   Description:
   Write a CCD  Table into a FITS Table. This requires the FITS file
   to be already open, and does not close it at the end. If a CCD
   Table extension is already present in the file, it is first removed, and
   then the new one is written into a new extension.  
   
   Input:
   VimosTable *ccdTable  
   Pointer to the CCD Table to write to a FITS Table

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool (true if the write is succesful, false otherwise)
   
   Updates:
   15 Dec 99: Created (BG)
   02 Feb 00: Changed the handling of the table, from a standalone FITS file,
              to an extension of the primary FITS file (MS)

-------------------------------------------------------------------------------
*/

VimosBool
writeFitsCcdTable(VimosTable *ccdTable, fitsfile *fptr)
{
  int status;
 
  /* validate input */
  if (ccdTable == NULL) {
    cpl_msg_error("writeFitsCcdTable","NULL input table");
    return(VM_FALSE);
  }

  /* check table name, should be "CCD" */
  if ( strcmp(ccdTable->name, VM_CCD) ) {
    cpl_msg_error("writeFitsCcdTable","Invalid input table");
    return(VM_FALSE);
  }

  status = 0;
  ccdTable->fptr = fptr;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "CCD", 0, &status) ) {
    if(fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_error("writeFitsCcdTable","The function fits_delete_hdu has "
                  "returned an error (code %d)", status);
      return(VM_FALSE);
    }
  } else {
    status = 0;
  }

  /* create 1st HDU for BINARY TABLE */
  if (fits_create_tbl (ccdTable->fptr,BINARY_TBL,0,0,NULL, NULL,NULL,
		       "CCD", &status)) {
    cpl_msg_error("writeFitsCcdTable","The function fits_create_tbl has "
                "returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "CCD", 0, &status) ) {
    cpl_msg_error("writeFitsCcdTable","The function fits_movnam_hdu has "
                "returned an error (code %d)", status);
    return (VM_FALSE);
  }

  /* write the table to the Fits file */
  if (!writeDescsToFitsTable(ccdTable->descs, ccdTable->fptr)) {
    cpl_msg_error("writeFitsCcdTable","The function writeDescsToFitsTable "
                "has returned an error");
    return(VM_FALSE);
  }

  return(VM_TRUE);
}
