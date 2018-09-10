/* $Id: vmimgastrometry.c,v 1.6 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilstrutils.h>
#include <piltranslator.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmtablearray.h"
#include "vmastrometrictable.h"
#include "vmstarmatchtable.h"
#include "vmwcsutils.h"
#include "vmimgextraction.h"
#include "vmimgastrometry.h"
#include "cpl.h"


/**
 * @defgroup vmimgastrometry Astrometry Utilities
 *
 * The module  provides the utility functions to determine
 * the correct astrometry (CD and CO matrices) for a VIMOS image.
 */

/**@{*/

/**
 * @brief
 *   Fit a CD matrix (pointing and rotation)
 *
 * @return @c VM_TRUE or @c VM_FALSE
 * 
 * @param calTables  Set of star match tables from different quadrants.
 * @param refImage   Reference image.
 *
 * The function merges the star match tables from the input set @em set
 * into a single star match table having the origin of the coordinate
 * system a pixel (0, 0) (given by the CRPIXi keywords). The new CD
 * matrix is computed fitting pixels versus right ascension and
 * declination.
 *
 * The reference image @em refImage is used to obtain the image size
 * in pixels of the images from which the input star match tables were
 * created. It is assumed that all images have the same size.
 *
 * The headers of each star match table in the input set is updated with
 * the computed CD matrix.
 *
 * @note
 *   To work properly the pixel coordinates of the input star match tables
 *   must be corrected for temperature effects.
 */

VimosBool
VmAstroComputeCD(VimosTableArray *calTables, const VimosImage *refImage)
{

  const char *modName = "VmAstroComputeCD";

  char *comment = "";

  int j;
  int nxpix, nypix;
  
  VimosTable  *stmcTable;

  struct WorldCoor *wcs = 0L;


  /*
   * Reset the input star match tables to same world coordinate system
   * and merge them into a single new table
   */

  
  if (!(stmcTable = shiftStarMatch(tblArrayGetData(calTables)))) {
    cpl_msg_error(modName,"Merging star match tables failed!");
    return VM_FALSE;
  }
    

  /*
   * Read common world coordinate system from the merged star match table
   * and update the WCS structure with the image size in pixels from
   * the header of the reference image.
   */
  
  if (readIntDescriptor(refImage->descs, pilTrnGetKeyword("Naxis", 1),
                        &nxpix, 0L) == VM_FALSE)
      return VM_FALSE;

  if (readIntDescriptor(refImage->descs, pilTrnGetKeyword("Naxis", 2),
                        &nypix, 0L) == VM_FALSE)
      return VM_FALSE;

  writeIntDescriptor(&(stmcTable->descs), "NAXIS", 2, comment);
  writeIntDescriptor(&(stmcTable->descs), pilTrnGetKeyword("Naxis", 1),
                     nxpix, comment);
  writeIntDescriptor(&(stmcTable->descs), pilTrnGetKeyword("Naxis", 2),
                     nypix, comment);

  if (!(wcs = rdimage(stmcTable->descs)))
    return VM_FALSE;

  removeDescriptor(&(stmcTable->descs), pilTrnGetKeyword("Naxis", 1));
  removeDescriptor(&(stmcTable->descs), pilTrnGetKeyword("Naxis", 2));

  wcs->nxpix = nxpix;
  wcs->nypix = nypix;


  /*
   * Fit CD matrix on the merged table.
   */

  if (!vimosFitMatch(wcs, stmcTable, stmcTable->cols->len))
    return VM_FALSE;
    
  deleteTable(stmcTable);


  /*
   * Update star match tables with new world coordinate system (but
   * CRPIX is the original one).
   */

  for (j = 0; j < 4; j++) {
    int k = 0;
    int ki, kj;

    VimosTable *table = (VimosTable *)tblArrayGet(calTables, j);

    for (ki = 1; ki <= 2; ki++) {
      writeDoubleDescriptor(&table->descs, pilTrnGetKeyword("Crval", ki),
                             wcs->crval[ki-1], pilTrnGetComment("Crval"));
      writeStringDescriptor(&table->descs, pilTrnGetKeyword("Ctype", ki),
                            wcs->ctype[ki-1], pilTrnGetComment("Ctype"));

      for(kj = 1; kj <= 2; kj++) {
        writeDoubleDescriptor (&table->descs, pilTrnGetKeyword("CD", ki, kj),
                               wcs->cd[k], pilTrnGetComment("CD"));
        k++;
      }
    }
  }

  vimoswcsfree(wcs);

  return VM_TRUE;

}


/**
 * @brief
 *   Fit a CO matrix (scale and distortions)
 *
 * @param calImage   The image for which CO has to be computed
 * @param tolrad     Search radius.
 * @param sigClip    Threshold for sigma clipping (?)
 * @param stmcTable  Star match table.
 * @param tflag      Flag to enable/disable temperature checks
 * @param tdevmax    Tolerance in C for beam temperature
 *
 * @return @c VM_TRUE or @c VM_FALSE
 * 
 * The function reads the CD matrix from the image @em calImage,
 * or, if the a star match table was given, i.e. @em stmcTable
 * is different from @c NULL, it is taken from @em stmcTable.
 * 
 * An artificial star table with a grid of 10 times 10 points is
 * created and the CCD to sky coordinate transformation is applied
 * to them, computing right ascension and declination. The right
 * ascension and declination are used to build an astrometric table,
 * and from both tables, star table and astrometric table, a star match
 * table is created and the CO matrix is fitted. The result is an image
 * @em calImage having a CO matrix which comprises wcs and CcdToSky.
 *
 * For details on the temperature check see the documentation for
 * @b computeVirtualPixels
 *
 * @see computeVirtualPixels
 */

VimosBool
VmAstroComputeCO(VimosImage *calImage, float tolrad, float sigClip,
                 VimosTable *stmcTable, unsigned int tflag, double tdevmax)
{
  const char *modName = "VmAstroComputeCO";


  char        comment[80], filterName[80];
  int         quad;
  int         naxis1, naxis2, nrow1;
  int         count;
  int         i,j,kj;
  double      rms[4];

  double rmsTol = 1.E-2;     /* tolerance for fitting in arcseconds */

  /*
   * Star Match Parameters: in this case magnitudes are fake, so anything
   * is OK
   */

  double tolMag1 = 30 , tolMag2 = 0.1;
  int minStar = 1;

  VimosTable *astroTable;
  VimosTable *newStmcTable;
  VimosTable *new2StmcTable;

  VimosColumn *RaCol, *xWorldCol, *NumCol, *IdCol, *xpixCol, *ypixCol;
  VimosColumn *DecCol, *yWorldCol;
  VimosColumn *xvirtCol, *yvirtCol;

  struct WorldCoor *wcs = 0;


  /*
   * We get Sky To Ccd coeffs from image
   */

  readIntDescriptor(calImage->descs, pilTrnGetKeyword("Naxis",1), &naxis1,
                    comment);
  readIntDescriptor(calImage->descs, pilTrnGetKeyword("Naxis",2), &naxis2,
                    comment);

  readIntDescriptor(calImage->descs, pilTrnGetKeyword("Quadrant"), &quad,
                    comment);
  readStringDescriptor(calImage->descs, pilTrnGetKeyword("FilterName", quad),
                       filterName, comment);

  /*
   * Now select the first character of the filter name that should indicate 
   * the wavelength band: 
   */
  
  filterName[1] = '\0';


  /*
   * In case of 4 quadrants, take wcs From StarMatch. The four quadrant
   * case is recognized by passing a non-NULL pointer for stmcTable.
   */

  if (stmcTable) {

    /*
     * As it is a table, we don't have naxis k/w..., an we need them for wcs!
     * thus we take them from first image, and just add to table, temporarily
     */

    writeIntDescriptor(&(stmcTable->descs), pilTrnGetKeyword("Naxis", 1),
                       naxis1, "");
    writeIntDescriptor(&(stmcTable->descs), pilTrnGetKeyword("Naxis", 2),
                       naxis2, "");

    if (!(wcs = rdimage(stmcTable->descs))) {
      return VM_FALSE;
    } 

    removeDescriptor(&(stmcTable->descs), pilTrnGetKeyword("Naxis", 1));
    removeDescriptor(&(stmcTable->descs), pilTrnGetKeyword("Naxis", 2));
    
    /*
     * and put it into image, because ImStarMatch re-reads it from image
     */

    upheader(calImage, wcs,  rms);

  }
  else {

    /*
     * In the case of 1 Quadrant, it is the same as from image
     */

    if (!(wcs = rdimage(calImage->descs))) {
      return VM_FALSE;
    } 
  }


  /*
   * What we want in StarMatch Table now is a grid of countxcount points 
   * we don't care whatever there was before
   */

  count = 10;

  newStmcTable = resetStarMatchTable(count, naxis1, naxis2);
  if (!newStmcTable) {
    cpl_msg_error(modName,"Could not make new StarMatch Table");
    return VM_FALSE;
  }
       
  if (computeVirtualPixels(calImage->descs, newStmcTable, tflag,
                           tdevmax) == VM_FALSE) {
    cpl_msg_error(modName, "Unable to apply Sky To CCd matrix");
    return VM_FALSE;
  }

  /*
   * Finally, use columns X-Image and Y_Image from StarMatchTable, convert
   * them to wcs using Sky To Ccd and NEW CD
   */

  nrow1 = newStmcTable->cols->len;
  pixtowcs(nrow1, newStmcTable, wcs);   


  /*
   * We now have a starMatch which in reality is a sort of AstroTable
   * as it contains xpixvirtual Ypixvirtual and Ra DEC. These are the real
   * RA DEC I put it in an AstroTable structure
   */

  astroTable = resetAstroTable(count,filterName);

  if (!astroTable) {
    cpl_msg_error(modName,"Could not make new astro Table");
    return VM_FALSE;
  }

  IdCol = findColInTab(astroTable, "ID");
  RaCol = findColInTab(astroTable, "RA");
  DecCol = findColInTab(astroTable, "DEC");

  NumCol = findColInTab(newStmcTable, "NUMBER");
  xWorldCol = findColInTab(newStmcTable, "X_WORLD");
  yWorldCol = findColInTab(newStmcTable, "Y_WORLD");


  /* delete next 4 before flight */

      xvirtCol = findColInTab(astroTable, "X_IMAGE");
      yvirtCol = findColInTab(astroTable, "Y_IMAGE");
      xpixCol = findColInTab(newStmcTable, "X_IMAGE");
      ypixCol = findColInTab(newStmcTable, "Y_IMAGE");


  for (i = 0; i < nrow1; i++) {
    RaCol->colValue->dArray[i] = xWorldCol->colValue->dArray[i];
    DecCol->colValue->dArray[i] = yWorldCol->colValue->dArray[i];
    sprintf(IdCol->colValue->sArray[i], "%d", NumCol->colValue->iArray[i]);
    
    /* delete next 2 before flight */
        xvirtCol->colValue->dArray[i] = xpixCol->colValue->dArray[i];
        yvirtCol->colValue->dArray[i] = ypixCol->colValue->dArray[i];
  }


  /*
   * And now reset the X-IMAGE Y_IMAGE columns of stmcTable
   * to the original grid
   */

  xpixCol = findColInTab(newStmcTable, "X_IMAGE");
  kj=0;
  for(i=0; i<count; i++) { 
    for(j=0; j<count; j++,kj++) {
      xpixCol->colValue->dArray[kj] = (j+1)*naxis1/(count + 1);
    }
  }

  ypixCol = findColInTab(newStmcTable, "Y_IMAGE");
  kj=0;
  for(i=0; i<count; i++) { 
    for(j=0; j<count; j++,kj++) {
      ypixCol->colValue->dArray[kj] = (i+1)*naxis2/(count + 1);
    }
  }


  /*
   * Thus we have a StarMatch Table with Xpix Ypix, and an astroTable with
   * RA DEC  we should now start from Xpix Ypix 
   * and fit a CO  to RA DEC, keeping CRVAL and CRPIX fixed (last param=0)
   */

  if ((new2StmcTable = VmImBuildStarMatchTable(calImage, newStmcTable, 
                                                astroTable,
                                      minStar, tolrad, tolMag1, tolMag2,
                                      sigClip)) == NULL) {
    cpl_msg_error(modName, "Failure in making the Star Match Table");
    return VM_FALSE;
  }
         

  if (fitCO(calImage, astroTable, new2StmcTable, minStar, tolrad, tolMag1,
            tolMag2, sigClip, rmsTol) != VM_TRUE) {
    cpl_msg_error(modName, "Failure in making the Star Match Table");
    return VM_FALSE;
  }

  deleteTable(newStmcTable);
  deleteTable(new2StmcTable);
  deleteTable(astroTable);
  
  return VM_TRUE;

}


/**
 * @brief
 *   Build Astrometric Table
 *
 * @return VimosTable
 * 
 * @param count       N. of entries
 * @param filterName  Image filter
 *
 * Build up an astrometric Table (empty) of count*count rows
 * for filter FilterName (used by VmImComputeCO)
 *
 * @author B.Garilli
 */

VimosTable *resetAstroTable (int count, char *filterName) 
{

  int         totcount, i;
  char        colName[6];

  VimosColumn *currentColumn;
  VimosTable *astroTable;

  totcount = count*count;

  astroTable = newAstrometricTable(); 
  astroTable->numColumns = 7;

  /*  only number, X_image and Y_image are actually needed */  

  astroTable->cols = newStringColumn(totcount,"ID");
  astroTable->cols->len= totcount;
  currentColumn = astroTable->cols;
  for(i=0; i<totcount; i++) 
    currentColumn->colValue->sArray[i] = cpl_strdup("    ");
        
  currentColumn->next = newDoubleColumn(totcount, "RA");
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++) 
    currentColumn->colValue->dArray[i] = 0.; 
      
  currentColumn->next = newDoubleColumn(totcount, "DEC");
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++)
    currentColumn->colValue->dArray[i] = 0.;

  sprintf(colName,"MAG_%s",filterName);
  currentColumn->next = newDoubleColumn(totcount, colName);
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++)
    currentColumn->colValue->dArray[i] = 0.;

  currentColumn->next = newDoubleColumn(totcount, "X_IMAGE");
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++)
    currentColumn->colValue->dArray[i] = 0.;

  currentColumn->next = newDoubleColumn(totcount, "Y_IMAGE");
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++)
    currentColumn->colValue->dArray[i] = 0.;
  

  currentColumn->next = newIntColumn(totcount, "GOFF");
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++)
    currentColumn->colValue->iArray[i] = 0;
      
  currentColumn->next = NULL;

  return(astroTable);

}


/**
 * @brief
 *   Build StarMatchTable Table
 *
 * @return VimosTable
 * 
 * @param count       N. of entries
 *
 * Build up an StarTable Table of count*count rows
 * Fill it up with an equally spaced grid of count*count points
 * (used by VmImComputeCO)
 *
 * @author B.Garilli
 */

VimosTable *resetStarMatchTable (int  count, int naxis1, int naxis2) 
{

  int          totcount, i, j, kj;
  VimosColumn *currentColumn;
  VimosTable  *starMatchTable;

  starMatchTable = newStarMatchTableEmpty(); 
 
  /*  only number, X_image and Y_image are actually needed */  
  totcount = count*count;

  starMatchTable->numColumns = 6;

  starMatchTable->cols = newIntColumn(totcount, "NUMBER");
  starMatchTable->cols->len = totcount;

  currentColumn = starMatchTable->cols;
  for(i=0; i<totcount; i++) 
    currentColumn->colValue->iArray[i] = i+1;
        
  currentColumn->next = newDoubleColumn(totcount, "MAG");
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++) 
    currentColumn->colValue->dArray[i] = 0.;
      
  currentColumn->next = newDoubleColumn(totcount, "X_IMAGE");
  currentColumn = currentColumn->next;
  kj=0;
  for(i=0; i<count; i++) {
    for(j=0; j<count; j++,kj++) {
      currentColumn->colValue->dArray[kj] = (j+1)*naxis1/(count + 1);
    }
  }
      
  currentColumn->next = newDoubleColumn(totcount, "Y_IMAGE");
  currentColumn = currentColumn->next;
  kj=0;
  for(i=0; i<count; i++) { 
    for(j=0; j<count; j++,kj++) {
      currentColumn->colValue->dArray[kj] = (i+1)*naxis2/(count + 1);
    }
  }

      
  currentColumn->next = newDoubleColumn(totcount, "X_WORLD");
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++) 
    currentColumn->colValue->dArray[i] = 0.; 
      
  currentColumn->next = newDoubleColumn(totcount, "Y_WORLD");
  currentColumn = currentColumn->next;
  for(i=0; i<totcount; i++)
    currentColumn->colValue->dArray[i] = 0.;
  currentColumn->next = NULL;

  return(starMatchTable);

}


/**
 * @brief
 *   Build StarMatchTable Table
 *
 * @return VimosTable
 * 
 * @param **starMatchTable       Array of 4 StarMatchTables (1 per quadrant)
 *
 * put together the 4 tables, different for each quadrant
 * we must change the image coordinates, and then change the CRPIX1 and 2
 * values to 0,0 
 *
 * for positive X quads (quad 1 and 4) it is:
 * new X = oldX-crpix1-0.5
 * for neg. X quads (2 and 3)
 * new X = -crpix1+oldX+0.5
 * for pos. Y quads (1 and 2)
 * new Y = oldY-crpix2-0.5
 * for neg. Y quads (3 and 4)
 * new Y = -crpix2+oldY+0.5
 *
 * @author B.Garilli
 */

VimosTable *shiftStarMatch (VimosTable **starMatchTable)
{ 
  char   modName[] = "shiftStarMatch";
  char comment[80] = "";
  char smydesc[80];

  int quad, i,j,k, nrow, nrowTot=0;
  double crpix1, crpix2, mydesc;

  VimosTable *megaStarMatchTable;
  VimosColumn *outXimaCol, *outYimaCol, *outXwldCol, *outYwldCol;
  VimosColumn *outNameCol, *outMagCol;
  VimosColumn *inXimaCol, *inYimaCol, *inXwldCol, *inYwldCol;
  VimosColumn *inNameCol, *inMagCol;


  megaStarMatchTable = newStarMatchTableEmpty(); 
 
  /* I need a number of descriptors (essentially the wcs realted ones) */

  readDoubleDescriptor(starMatchTable[0]->descs, pilTrnGetKeyword("Equinox"),
                       &mydesc, comment);
  writeDoubleDescriptor(&(megaStarMatchTable->descs),
                        pilTrnGetKeyword("Equinox"), mydesc, comment);
  readStringDescriptor(starMatchTable[0]->descs, pilTrnGetKeyword("Radecsys"),
                       smydesc, comment);
  writeStringDescriptor(&(megaStarMatchTable->descs),
                        pilTrnGetKeyword("Radecsys"), smydesc, comment);
  for (i=1; i<=2; i++) {
    readDoubleDescriptor(starMatchTable[0]->descs, pilTrnGetKeyword("Crval",i),
                         &mydesc, comment);
    writeDoubleDescriptor(&(megaStarMatchTable->descs),
                          pilTrnGetKeyword("Crval",i), mydesc, comment);
    readStringDescriptor(starMatchTable[0]->descs, pilTrnGetKeyword("Ctype",i),
                         smydesc, comment);
    writeStringDescriptor(&(megaStarMatchTable->descs),
                          pilTrnGetKeyword("Ctype",i),smydesc, comment);
    for(j=1; j<=2; j++,k++) {
      readDoubleDescriptor(starMatchTable[0]->descs, 
                           pilTrnGetKeyword("CD", i, j), &mydesc, comment);
      writeDoubleDescriptor(&(megaStarMatchTable->descs),
                            pilTrnGetKeyword("CD", i, j), mydesc, comment);
    }
  }

  for (i=0; i<4;i++) {
    nrowTot += starMatchTable[i]->cols->len;
  }

  megaStarMatchTable->numColumns = 6;
  megaStarMatchTable->cols = newIntColumn(nrowTot, "NUMBER");
  outNameCol = megaStarMatchTable->cols;
  megaStarMatchTable->cols->next = newDoubleColumn(nrowTot, "MAG");
  outMagCol = megaStarMatchTable->cols->next;
  megaStarMatchTable->cols->next->next = newDoubleColumn(nrowTot, "X_IMAGE");
  outXimaCol = megaStarMatchTable->cols->next->next;
  megaStarMatchTable->cols->next->next->next = newDoubleColumn(nrowTot,
                                                               "Y_IMAGE");
  outYimaCol = megaStarMatchTable->cols->next->next->next;
  megaStarMatchTable->cols->next->next->next->next = newDoubleColumn(nrowTot, "X_WORLD");
  outXwldCol = megaStarMatchTable->cols->next->next->next->next ;
  megaStarMatchTable->cols->next->next->next->next->next = newDoubleColumn(nrowTot, "Y_WORLD");
  outYwldCol = megaStarMatchTable->cols->next->next->next->next->next;

  k=0;
  for (i=0; i<4;i++) {
    nrow = starMatchTable[i]->cols->len; 

    readIntDescriptor(starMatchTable[i]->descs, pilTrnGetKeyword("Quadrant"),
                      &quad, comment);
    readDoubleDescriptor(starMatchTable[i]->descs, pilTrnGetKeyword("Crpix",1),
                         &crpix1, comment);
    readDoubleDescriptor(starMatchTable[i]->descs, pilTrnGetKeyword("Crpix",2),
                         &crpix2, comment);

    if(!(inXimaCol = findColInTab(starMatchTable[i], "X_IMAGE"))) {
      cpl_msg_error(modName, "Star Table: Column with X-pixel coord "
                  "not found");
      return NULL;
    }
    if(!(inYimaCol = findColInTab(starMatchTable[i], "Y_IMAGE"))) {
      cpl_msg_error(modName, "Star Table: Column with Y-pixel coord "
                  "not found");
      return NULL;
    }
    if(!(inMagCol = findColInTab(starMatchTable[i], "MAG"))) {
      cpl_msg_error(modName, "Star Table: Column with Y-pixel coord "
                  "not found");
      return NULL;
    }
    if(!(inNameCol = findColInTab(starMatchTable[i], "NUMBER"))) {
      cpl_msg_error(modName, "Star Table: Column with Y-pixel coord "
                  "not found");
      return NULL;
    }
    if(!(inXwldCol = findColInTab(starMatchTable[i], "RA"))) {
      cpl_msg_error(modName, "Star Table: Column with Y-pixel coord "
                  "not found");
      return NULL;
    }
    if(!(inYwldCol = findColInTab(starMatchTable[i], "DEC"))) {
      cpl_msg_error(modName, "Star Table: Column with Y-pixel coord "
                  "not found");
      return NULL;
    }
    
    for (j=0;j<nrow;j++) {
      if (quad == 1 || quad == 4 ) {
        outXimaCol->colValue->dArray[k] = inXimaCol->colValue->dArray[j] -
          crpix1 -0.0;
      } else {
        outXimaCol->colValue->dArray[k] = -crpix1 + 
          inXimaCol->colValue->dArray[j]+0.0;
      }
      if (quad == 1 || quad == 2 ) {
        outYimaCol->colValue->dArray[k] = inYimaCol->colValue->dArray[j] -
          crpix2 -0.0;
      } else {
        outYimaCol->colValue->dArray[k] = -crpix2 +
          inYimaCol->colValue->dArray[j]+0.0;
      }
      outXwldCol->colValue->dArray[k] = inXwldCol->colValue->dArray[j];
      outYwldCol->colValue->dArray[k] = inYwldCol->colValue->dArray[j];
      outMagCol->colValue->dArray[k] = inMagCol->colValue->dArray[j];
      outNameCol->colValue->iArray[k] = inNameCol->colValue->iArray[j];
        
      k += 1;
    }
    
    crpix1 = 0.;
    crpix2=0.;
    
    writeDoubleDescriptor(&(megaStarMatchTable->descs), 
                          pilTrnGetKeyword("Crpix",1), crpix1, comment);
    writeDoubleDescriptor(&(megaStarMatchTable->descs), 
                          pilTrnGetKeyword("Crpix",2), crpix2, comment);
    
  }

  return(megaStarMatchTable);

}


VimosBool fitCO(VimosImage *image, VimosTable *astroTable,
                VimosTable *starMatch,int minstar, double tolrad,
                double tolMag1, double tolMag2, float sigClip,double rmsTol)
{

  struct WorldCoor    * wcs = 0;     /* WCS structure */

  int      i,nastro,nstarMatch;
  int      vmncoef = 13;           /* Number of fitting coefficents */
  char     modName[] = "fitCO";

  int    * imatch = NULL;     /* index for matching stars */
  int      nmatch ;        /* Number of matches returned */
  double   chisq = 0;         /* Reduced chi-2 of plate model */
  double   rms[4] ;            /* RMS error in X and Y */



  /* Initialize the WCS structure from image descriptors  
   */
    
  for (i=0;i<=3;i++) rms[i] = 0;
  
  if  (!(wcs = rdimage(image->descs))) {       
    return VM_FALSE;
  } 
 
  /*  Find matching pairs and transformation coeff. */
   
  cpl_msg_info(modName,"Begin  to fit WCS on image");
   



  /* BG: this call to searchmatch is useless, as the input to
     fitCO is an already matched table, thus the number of matched stars
     is the length of the table */

/*    nstarMatch = starMatch->cols->len; */
/*    nastro = astroTable->cols->len; */
/*    nmatch = nstarMatch; */
/*    imatch = searchmatch(wcs, starMatch, nstarMatch, astroTable, nastro, */
/*                     &nmatch, tolrad, tolMag1, tolMag2, sigClip,minstar); */
        

  nstarMatch = starMatch->cols->len;
  nmatch = nstarMatch;

  cpl_msg_info(modName, "Fitting %d matching stars with a %d-coefficients "
             "polynomial", nmatch, vmncoef);
  
  if ((vimosFitPlate(wcs, starMatch, astroTable, nmatch, vmncoef, &chisq) == 
      VM_FALSE)) { 
    cpl_free(imatch);
    return VM_FALSE;
  }
   
    
  /* Project the reference stars into pixels on image plane */
  nastro = astroTable->cols->len;
    
  wcstopix(nastro,astroTable,wcs);
   
  /* 
   * Search for matching stars. First convert search radius from
   * arcsec to pixel.
   */
    
  cpl_msg_info(modName,"Searching for matching stars");

/*    imatch = searchmatch(wcs, starMatch, nstarMatch, astroTable, nastro, */
/*                     &nmatch, tolrad, tolMag1, tolMag2, sigClip,minstar); */

  tolrad /= fabs(3600. * wcs->cdelt[0]);

  imatch = VmSearchMatches(starMatch, astroTable, tolrad, 
                            tolMag1, tolMag2, sigClip,
                            minstar, &nmatch);

  if(!imatch || nmatch < minstar) {
    cpl_msg_warning(modName, "Insufficent number of matching stars: %d found",
                  nmatch);
    cpl_free(imatch);
    return VM_FALSE;
  }      
    
  /* Quality control : this is commented as now the fit is done
     only on a 10x10 grid of points: therefore the quality control
     is already done by calcres (add a check there, if wished..) */
    
    
  cpl_msg_info(modName,"Number of matching pairs is %d",nmatch);
    
  /* Project the reference stars into pixels and observed into WCS */
  /* and  compute RMS */
    
  wcstopix(nastro,astroTable,wcs);
  pixtowcs(nstarMatch,starMatch,wcs);
  calcres(starMatch, astroTable, imatch, nmatch, rms);
    
  cpl_msg_info(modName,"Computed RMS from model fit in X and Y:"
             "CCD->Sky: X_RMS error = %g (arcsec); "
             "Y_RMS error = %g (arcsec), Sky->CCD: XRMS error = %g (pixels)" 
             "Y_RMS error =%g (pixels)", rms[2],rms[3], rms[0],rms[1]);
  
  if (rms[2] > rmsTol || rms[3] > rmsTol) {
    cpl_msg_warning(modName, "CCD to Sky RMS is greater than expected: %g,%g "
                  "against %g", rms[2],rms[3],rmsTol);
  }
    
  /* update the image header */
  if(!(upheader(image, wcs, rms))) {
    cpl_msg_error(modName,"Image header cannot be be updated");
    return VM_FALSE;
  }
    
  /* Free Memory */
    
  if(wcs)     vimoswcsfree(wcs);
    
  return VM_TRUE;
}


/**
 * @brief 
 *   Compute the RMS error of the CCD <-> Sky transformation. 
 *   
 * @return VM_TRUE/VM_FALSE
 * 
 * @param wcs      wcs structure
 * @param o_star   star table
 * @param a_star   astrometric table
 * @param imatch   table of matching indexes
 * @param nmatch   number of matches
 * @param rms      pointer to the returned array
 *
 * Compute the RMS error of the CCD <-> Sky transformation.
 *
 * @author P. Montegriffo (i/o modified by P.Sartoretti)
 */

/*-----------------------------------------------------------------------------
 * Function  : calcres
 * In        : wcs structure, o_star structure, a_star structure, index table 
 *             for matching stars and number of matches
 * Out       : vector rms with mean dx,dy, x_rms, y_rms, covariance dxy
 * Purpouse  : Compute the RMS error in X and Y after that the transformation
 *             is applied; values are converted in [arcsec]
 * Note      : 
 *---------------------------------------------------------------------------*/

int calcres(VimosTable *o_star, VimosTable *a_star, int imatch[], int nmatch,
            double * rms)
{
  int    i;
  double dx, dy, dmatch, dxwld=0., dywld=0.;
  double dxsum = 0.;
  double dysum = 0.;
  double dywldsum = 0.;
  double dxwldsum = 0.;
  char   modName[] = "calcres";
  VimosColumn    *aXimaCol, *aYimaCol, *aXwldCol, *aYwldCol;
  VimosColumn    *oXimaCol, *oYimaCol, *oXwldCol, *oYwldCol;

  /****
       I have changed this function in order to estimate the RMS both of
       the direct (Sky->CCD) and the inverse transformation (CCD->Sky).
       Write all these values in rms[], and comment out what it was written
       before in rms[] (and never used). 
  ****/ 

  if(!(aXimaCol = findColInTab(a_star, "X_IMAGE"))) {
    cpl_msg_error(modName, "Astrometric Table: Column with X-pixel coord "
                "not found");
    return VM_FALSE;
  }
  if(!(aYimaCol = findColInTab(a_star, "Y_IMAGE"))) {
    cpl_msg_error(modName, "Astrometric Table: Column with Y-pixel coord "
                "not found");
    return VM_FALSE;
  }
  if(!(aXwldCol = findColInTab(a_star, "RA"))) {
    cpl_msg_error(modName, "Astrometric Table: Column with RA coord "
                "not found");
    return VM_FALSE;
  }
  if(!(aYwldCol = findColInTab(a_star, "DEC"))) {
    cpl_msg_error(modName, "Astrometric Table: Column with RA coord "
                "not found");
    return VM_FALSE;
  }
  if(!(oXimaCol = findColInTab(o_star, "X_IMAGE"))) {
    cpl_msg_error(modName, "Star Table: Column with X-pixel coord "
                "not found");
    return VM_FALSE;
  }
  if(!(oYimaCol = findColInTab(o_star, "Y_IMAGE"))) {
    cpl_msg_error(modName, "Star Table: Column with Y-pixel coord "
                "not found");
    return VM_FALSE;
  }
  if(!(oXwldCol = findColInTab(o_star, "X_WORLD"))) {
    cpl_msg_error(modName, "Star Table: Column with X-world coord "
                "not found");
    return VM_FALSE;
  }
  if(!(oYwldCol = findColInTab(o_star, "Y_WORLD"))) {
    cpl_msg_error(modName, "Star Table: Column with Y-world coord "
                "not found");
    return VM_FALSE;
  }
    
  for (i = 0; i < nmatch; i++) {
    dxwld = (aXwldCol->colValue->dArray[imatch[2*i+1]]) - 
      (oXwldCol->colValue->dArray[imatch[2*i]]);
    dywld = (aYwldCol->colValue->dArray[imatch[2*i+1]]) - 
      (oYwldCol->colValue->dArray[imatch[2*i]]); 
    dx = (aXimaCol->colValue->dArray[imatch[2*i+1]]) - 
      (oXimaCol->colValue->dArray[imatch[2*i]]);
    dy = (aYimaCol->colValue->dArray[imatch[2*i+1]]) - 
      (oYimaCol->colValue->dArray[imatch[2*i]]);

    /*
     * Next two lines added by C.Izzo, RMS computation was wrong 
     * in case coordinates are crossing the 360 degrees line.
     */

    if (fabs(fabs(dxwld) - 360.) < 0.1)
      dxwld = fabs(dxwld) - 360.;
    
    dxwld = fabs(dxwld)* 3600.;
    dywld = fabs(dywld) * 3600.;
    dx = fabs(dx);
    dy = fabs(dy);
    dxsum = dxsum + dx;
    dysum = dysum + dy;
    dxwldsum = dxwldsum + dxwld;
    dywldsum = dywldsum + dywld;
  }
  dmatch = (double) nmatch;
  *rms = dxsum / dmatch;
  *(rms+1) = dysum / dmatch;
  *(rms+2) = dxwldsum / dmatch;
  *(rms+3) = dywldsum / dmatch;
      
  return VM_TRUE;
}
/**@}*/
