/* $Id: vmextractiontable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <math.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>
#include <pilerrno.h>

#include "vmtable.h"
#include "vmadfifutable.h"
#include "vmdistmodels.h"
#include "vmextractiontable.h"
#include "cpl.h"


/*
  Construct a new Extraction Slit.
*/
VimosExtractionSlit *newExtractionSlit()
{
  VimosExtractionSlit *newSlit;
  
  /* allocate memory */
  newSlit = (VimosExtractionSlit *) cpl_malloc(sizeof(VimosExtractionSlit));

  /* check if space was allocated */
  if (newSlit == NULL) {
    cpl_msg_error("newExtractionSlit","Allocation Error");
    return(NULL);
  }

  /* set all field to 0 */
  newSlit->slitNo  = 0;
  newSlit->numRows = 0;
  newSlit->IFUslitNo = 0;
  newSlit->IFUfibNo = 0;
  newSlit->IFUfibPeakX = 0.;
  newSlit->IFUfibTrans = 0.;
  newSlit->width = 0.;
  newSlit->y       = NULL;
  newSlit->ccdX    = NULL;
  newSlit->ccdY    = NULL;
  newSlit->maskX   = NULL;
  newSlit->maskY   = NULL;
  newSlit->numSpec = NULL;
  newSlit->crvPol  = NULL;
  newSlit->crvPolRms = NULL;
  newSlit->invDis  = NULL;
  newSlit->invDisRms = NULL;
  newSlit->invDisQuality = NULL;
  newSlit->zeroX   = NULL;
  newSlit->zeroY   = NULL;
  newSlit->prev    = NULL;
  newSlit->next    = NULL;
  

  /* return address of new Extraction Slit */
  return(newSlit);
  
}

/*
  Delete an Extraction Slit.
*/
void deleteExtractionSlit(VimosExtractionSlit *slit)
{
  VimosExtractionSlit *tmpSlit;
  VimosExtractionSlit *nxtSlit;
  int i;
  
  /* check input */
  if (slit == NULL) {
    return;
  }
  
  
  /* loop through all slits in the slit */
  tmpSlit = slit;
  while (tmpSlit) {
    /* get oiter to next slit */
    nxtSlit = tmpSlit->next;
  
    /* free members */
    deleteIntArray(tmpSlit->y);
    deleteFloatArray(tmpSlit->ccdX);
    deleteFloatArray(tmpSlit->ccdY);
    deleteFloatArray(tmpSlit->maskX);
    deleteFloatArray(tmpSlit->maskY);
    deleteIntArray(tmpSlit->numSpec);
    deleteFloatArray(tmpSlit->zeroX);
    deleteFloatArray(tmpSlit->zeroY);
    deleteFloatArray(tmpSlit->crvPolRms);
    deleteFloatArray(tmpSlit->invDisRms);
    if (tmpSlit->crvPol && tmpSlit->invDis) {
      for (i = 0; i < tmpSlit->numRows; i++) {
        deleteDistModel1D(tmpSlit->crvPol[i]);
        deleteDistModel1D(tmpSlit->invDis[i]);
      }
      cpl_free(tmpSlit->crvPol);
      cpl_free(tmpSlit->invDis);
    }
    cpl_free(tmpSlit);
    
    /* address of next slit */
    tmpSlit = nxtSlit;
  }
  
}



/*
  Construct a new Extraction Table.
*/
VimosExtractionTable *newExtractionTable()
{
  const char            modName[] = "newExtractionTable";
  VimosExtractionTable *newTab;
  
  /* allocate memory */
  newTab = (VimosExtractionTable *) cpl_malloc(sizeof(VimosExtractionTable));

  /* check if space was allocated */
  if (newTab == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  
  /* copy "EXT" into name of table */
  strcpy(newTab->name, VM_EXT);
  
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_EXT, "");
  if (newTab->descs == NULL) {
    /* cleanup */
    cpl_free(newTab);
    cpl_msg_error(modName, "Function newStringDescriptor failure");
    return(NULL);
  }
  newTab->slits = NULL;
  newTab->fptr = NULL;
  
  /* return address of new Extraction Table */
  return(newTab);
  
}

/*
  Delete a Extraction Table. This is just an esthetic wrapper for 
  deleteTable(extTable)
*/
void deleteExtractionTable(VimosExtractionTable *extTable)
{
  VimosDescriptor *tmpDesc;
  VimosDescriptor *nxtDesc;
  
  if (extTable == NULL) {
    return;
  }
  
  deleteExtractionSlit(extTable->slits);

  tmpDesc = extTable->descs;  
  while (tmpDesc) {
    nxtDesc = tmpDesc->next;
    deleteDescriptor(tmpDesc);
    tmpDesc = nxtDesc;
  }

}

VimosBool copyGrsTab2ExtTab(VimosTable *grsTable, 
                            VimosExtractionTable *extTable)
{
  const char  modName[] = "copyGrsTab2ExtTab";

 /* 
  *  Some descriptors should not be copied!! e.g. TABLE 
  */  
  if (!copyAllDescriptors(grsTable->descs, &(extTable->descs))) {
    cpl_msg_error(modName, "Function copyAllDescriptors failure");
    return(VM_FALSE);
  }

  if ((writeStringDescriptor(&(extTable->descs), pilTrnGetKeyword("Table"), 
                            VM_EXT, "")) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot write descriptor %s", 
                pilTrnGetKeyword("Table"));
    return VM_FALSE;
  }

  if (!writeStringDescriptor(&(extTable->descs), "EXTNAME", "EXR", "")) {
    cpl_msg_error(modName, "Function writeStringDescriptor failure");
    return(VM_FALSE);
  }

  return VM_TRUE;
}


VimosBool copyAdf2ExtTab(VimosTable *adf, VimosExtractionTable *extTable)
{
  const char  modName[] = "copyAdf2ExtTab";
  VimosDescriptor  *firstDesc;
  VimosDescriptor  *copyDesc;
  double            descValue;

  firstDesc = findDescriptor(adf->descs, pilTrnGetKeyword("Quadrant"));
  copyDesc = copyOfDescriptor(firstDesc);
  if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
    cpl_msg_debug(modName, "The function addDesc2Table has returned an error");
    return VM_FALSE;
  } 
  
  firstDesc = findDescriptor(adf->descs, "ESO PRO IDS*");
  while (firstDesc) {
    if (strstr(firstDesc->descName, "DAYTIM") || 
	strstr(firstDesc->descName, "ORD") || 
	strstr(firstDesc->descName, "RMS") || 
	strstr(firstDesc->descName, "TEMP")) {
      copyDesc = copyOfDescriptor(firstDesc);
      if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
	cpl_msg_debug(modName, "The function addDesc2Table has returned an error");
	return(VM_FALSE);
      } 
    }
    else {
      descValue = atof(firstDesc->descValue->s);
      copyDesc = newDoubleDescriptor(firstDesc->descName, descValue, 
				     firstDesc->descComment);
      if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
	cpl_msg_debug(modName, "The function addDesc2Table has returned an error");
	return(VM_FALSE);
      } 
    }
    firstDesc = firstDesc->next;
    firstDesc = findDescriptor(firstDesc, "ESO PRO IDS*");
  }

  firstDesc = findDescriptor(adf->descs, "ESO PRO OPT*");
  while (firstDesc) {
    if (strstr(firstDesc->descName, "DAYTIM") || 
	strstr(firstDesc->descName, "ORD") || 
	strstr(firstDesc->descName, "RMS") || 
	strstr(firstDesc->descName, "TEMP")) {
      copyDesc = copyOfDescriptor(firstDesc);
      if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
	cpl_msg_debug(modName, "The function addDesc2Table has returned an error");
	return(VM_FALSE);
      } 
    }
    else {
      descValue = atof(firstDesc->descValue->s);
      copyDesc = newDoubleDescriptor(firstDesc->descName, descValue, 
				     firstDesc->descComment);
      if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
	cpl_msg_debug(modName, "The function addDesc2Table has returned an error");
	return(VM_FALSE);
      } 
    }
    firstDesc = firstDesc->next;
    firstDesc = findDescriptor(firstDesc, "ESO PRO OPT*");
  }

  firstDesc = findDescriptor(adf->descs, "ESO PRO CRV*");
  while (firstDesc) {
    if (strstr(firstDesc->descName, "ORD") || 
	strstr(firstDesc->descName, "RMS")) {
      copyDesc = copyOfDescriptor(firstDesc);
      if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
	cpl_msg_debug(modName, 
                    "The function addDesc2Table has returned an error");
	return(VM_FALSE);
      } 
    }
    else {
      descValue = atof(firstDesc->descValue->s);
      copyDesc = newDoubleDescriptor(firstDesc->descName, descValue, 
				     firstDesc->descComment);
      if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
	cpl_msg_debug(modName, "The function addDesc2Table has returned an error");
	return(VM_FALSE);
      } 
    }
    firstDesc = firstDesc->next;
    firstDesc = findDescriptor(firstDesc, "ESO PRO CRV*");
  }

  firstDesc = findDescriptor(adf->descs, "ESO PRO ZERO*");
  while (firstDesc) {
    if (strstr(firstDesc->descName, "ORD") || 
	strstr(firstDesc->descName, "RMS") ||
	strstr(firstDesc->descName, "WIDTH")) {
      copyDesc = copyOfDescriptor(firstDesc);
      if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
	cpl_msg_debug(modName, "The function addDesc2Table has returned an error");
	return(VM_FALSE);
      } 
    }
    else {
      descValue = atof(firstDesc->descValue->s);
      copyDesc = newDoubleDescriptor(firstDesc->descName, descValue, 
				     firstDesc->descComment);
      if (!addDesc2Desc(copyDesc, &(extTable->descs))) {
	cpl_msg_debug(modName, "The function addDesc2Table has returned an error");
	return(VM_FALSE);
      } 
    }
    firstDesc = firstDesc->next;
    firstDesc = findDescriptor(firstDesc, "ESO PRO ZERO*");
  }

  return(VM_TRUE);
}  


/*
 * Using the positions of the slits as defined in the ADF, and using the
 * distortion models that are available in the Extraction Table (e.g. the
 * parameters of these models been copied from a Grism Table using
 * copyGrsTab2ExtTab), the columns of the Extraction Table are computed. 
 * Such a computed Extraction Table is used in the on-linepipeline to reduce
 * the data, while in the off-line pipeline it can be used as a first
 * approximation for the locations of the spectra on the CCD.
 * A few descriptors are copied from the ADF to the Extraction Table.
 */

/* 
   23 Oct 01: Added IFUfibPeakX parameter in the call to 
              calcSlitLocationsOnCCD (AZ)
*/

VimosBool computeExtractionTable(VimosTable *adf, VimosIfuTable *ifuTab,
                                 VimosExtractionTable *extTable)
{
  const char           modName[] = "computeExtractionTable";
  int                  i;
  int                  count = 0;
  int                  quadNum;
  int                  yCoord;
  int                  zeroCor;
  char                 comment[72];
  char                *descName;
  char                *obsType = "UNDEFINED";
  float                lambda0;  
  double               dValue;
  VimosBool            rdOK;
  VimosIfuMode         ifuMode;
  VimosDescriptor     *tDesc;
  VimosAdfSlitHolder  *slitHolder=NULL;
  VimosAdfSlitHolder  *tSlitHolder;
  VimosExtractionSlit *exSlit;
  VimosExtractionSlit *tSlit;
  VimosDistModel2D    *optModX;
  VimosDistModel2D    *optModY;
  VimosDistModel2D    *contModX;
  VimosDistModel2D    *contModY;
  VimosDistModelFull  *crvMod;
  VimosDistModelFull  *invDispMat;

 /* 
  *  Check table name, should be MOS- or IFU ADF 
  */
  if ( (strcmp(adf->name, VM_ADF_MOS) ) &&
       (strcmp(adf->name, VM_ADF_IFU) ) ) {
    cpl_msg_error(modName, "Invalid input table");
    return(VM_FALSE);
  }

 /* 
  *  Check table name, should be Extraction Table 
  */
  if ( strcmp(extTable->name, VM_EXT) ) {
    cpl_msg_error(modName, "Invalid input extraction table");
    return(VM_FALSE);
  }
  cpl_msg_debug(modName, "Computing extraction Table");

  ifuMode = VM_IFU_MODE_UNDEF;

 /* 
  *  Check if distortion model is present in Extraction Table 
  */

  readIntDescriptor(extTable->descs, pilTrnGetKeyword("Quadrant"), 
                      &quadNum, comment);

 /* 
  *  Set type of data 
  */

  if ( !strcmp(adf->name, VM_ADF_MOS) ) obsType = "MOS";
  if ( !strcmp(adf->name, VM_ADF_IFU) ) obsType = "IFU";

  if (writeStringDescriptor(&(extTable->descs), pilTrnGetKeyword("ObsType"), 
                            obsType, "") == VM_FALSE) {
    cpl_msg_error(modName, "Cannot write keyword %s to Extraction Table", 
                pilTrnGetKeyword("ObsType"));
  }

  if (obsType == "MOS") {
    /*
     *  Copy Mask indentification from ADF
     */

    descName = (char *) pilKeyTranslate("MaskId", quadNum);
    rdOK = copyFromHeaderToHeader(adf->descs, descName,
                                  &(extTable->descs), NULL);
    if (rdOK == VM_TRUE) {
      descName = (char *) pilKeyTranslate("AdfId", quadNum);
      rdOK = copyFromHeaderToHeader(adf->descs,descName,
                                    &(extTable->descs), NULL);
    }

    if (rdOK == VM_FALSE) {
      cpl_msg_error(modName, "Failure in copying %s from ADF to "
                  "Extraction Table", descName);
      return(VM_FALSE);
    }
  }

  if (obsType == "IFU") {
    if (!writeStringDescriptor(&(extTable->descs), pilTrnGetKeyword("MaskId", quadNum),
			       "IFU Mask", "")) {
      cpl_msg_error(modName, "Failure in writing keyword %s to Extraction Table",
		  pilTrnGetKeyword("MaskId", quadNum));
      return(VM_FALSE);
    }
    if (!writeStringDescriptor(&(extTable->descs), pilTrnGetKeyword("AdfId", quadNum),
			       "IFU Mask", "")) {
      cpl_msg_error(modName, "Failure in writing keyword %s to Extraction Table",
		  pilTrnGetKeyword("AdfId", quadNum));
      return(VM_FALSE);
    }

    rdOK = copyFromHeaderToHeader(adf->descs,
                                  pilTrnGetKeyword("IfuMode"),
                                  &(extTable->descs), NULL);
    if (rdOK == VM_FALSE) {
      cpl_msg_error(modName, 
                  "Failure in copying %s from ADF to Extraction Table", 
                  pilTrnGetKeyword("IfuMode"));
      return VM_FALSE;
    }
    /* ALEX!!! define descName, otherwise still set to ADF ID !!!*/
    descName = (char *) pilKeyTranslate("IfuMode");

    tDesc = findDescriptor(adf->descs, descName);

    if ( !strncmp(tDesc->descValue->s, "ON", 2) ) {
      ifuMode = VM_IFU_MODE_SMALL;
    }
    else if ( !strncmp(tDesc->descValue->s, "OFF", 3) ) {
      ifuMode = VM_IFU_MODE_LARGE;
    }
  }

  rdOK = readDoubleDescriptor(extTable->descs, pilTrnGetKeyword("WlenCen"), 
                              &dValue, comment);
  if (rdOK == VM_FALSE) {
    cpl_msg_error(modName, "Keyword %s not found",  pilTrnGetKeyword("WlenCen"));
    return VM_FALSE;
  }
  lambda0 = (float) dValue;

  if ( !strcmp(adf->name, VM_ADF_MOS) ) {
    slitHolder = extractSlitsFromADF(adf);
    if (!slitHolder) {
      cpl_msg_error(modName, "Function extractSlitsFromADF failure");
      return VM_FALSE;
    }
  } 
  else if ( !strcmp(adf->name, VM_ADF_IFU) ) {
    slitHolder = extractSlitsFromIFU(adf, ifuTab, ifuMode);
    if (!slitHolder) {
      cpl_msg_error(modName, "Function extractSlitsFromIFU failure");
      return VM_FALSE;
    }
  }
  else {
    cpl_msg_error(modName, "Unrecognized instrument mode");
  }

  if (!readOptDistModel(extTable->descs, &optModX, &optModY)) {
    cpl_msg_error(modName, "Function readOptDistModel failure");
    return VM_FALSE;
  }

  if (!readCurvatureModel(extTable->descs, &crvMod)) {
    cpl_msg_error(modName,"Function readCurvatureModel failure");
    return VM_FALSE;
  }

  if (!readInvDispMatrix(extTable->descs, &invDispMat)) {
    cpl_msg_error(modName, "Function readInvDispMatrix failure");
    return VM_FALSE;
  }

  if (readIntDescriptor(extTable->descs, pilTrnGetKeyword("ZeroOrderFlag"), 
                         &zeroCor, comment) == VM_FALSE) {
    cpl_msg_error(modName, "Keyword %s not found", 
                pilTrnGetKeyword("ZeroOrderFlag"));
    return VM_FALSE;
  }

  if (zeroCor) {
    if (!readContaminationModel(extTable->descs, &contModX, &contModY)) {
      cpl_msg_error(modName, 
                  "Function readContaminationModel returned an error");
      return VM_FALSE;
    }
  }
  else {
    contModX = NULL;
    contModY = NULL;
  }

  tSlit = NULL;
  if (slitHolder == NULL) return VM_FALSE;
  else tSlitHolder = slitHolder;

  yCoord = 0;
  
  while (tSlitHolder) {
    count++;
    exSlit = newExtractionSlit();
    if (exSlit == NULL) {
      cpl_msg_error(modName, "Function newExtractionSlit failure");
      return VM_FALSE;
    }

    if (!calcSlitLocationsOnCCD(tSlitHolder->slit, tSlitHolder->slitType, 
                           optModX, optModY, crvMod, invDispMat,  
                           contModX, contModY,
                           &(exSlit->ccdX), &(exSlit->ccdY), &(exSlit->maskX),
			   &(exSlit->maskY), &(exSlit->crvPol), 
			   &(exSlit->crvPolRms), &(exSlit->invDis), 
			   &(exSlit->invDisRms), lambda0, &(exSlit->numRows), 
                           &(exSlit->IFUslitNo), &(exSlit->IFUfibNo),
			   &(exSlit->IFUfibPeakX), &(exSlit->IFUfibTrans),
                           &(exSlit->zeroX), &(exSlit->zeroY), 
                           &(exSlit->invDisQuality))) {
      cpl_msg_error(modName, "Function calcSlitLocationsOnCCD failure");
      return VM_FALSE;
    }

    if (exSlit->numRows > 0) {
      if (tSlit == NULL) {
        extTable->slits = exSlit;
      }
      else {
        tSlit->next = exSlit;
        exSlit->prev = tSlit;
      }
      tSlit = exSlit;
      tSlit->slitNo = tSlitHolder->slitNo;
      exSlit->y = newIntArray(exSlit->numRows);
      if (exSlit->y == NULL) {
        deleteExtractionSlit(exSlit);
        cpl_msg_error(modName, "Function newIntArray failure");
        return(VM_FALSE);
      }

      exSlit->numSpec = newIntArray(exSlit->numRows);
      if (exSlit->numSpec == NULL) {
        deleteExtractionSlit(exSlit);
        cpl_msg_error(modName, "Function newIntArray failure");
        return(VM_FALSE);
      }

      for (i=0; i < tSlit->numRows; i++) {
        tSlit->y->data[i] = yCoord;
        exSlit->numSpec->data[i] = 1;
        yCoord++;
      }

      if (readFloatDescriptor(adf->descs, pilTrnGetKeyword("SlitDimY", count), 
                              &(exSlit->width), comment) == VM_FALSE) {
        cpl_msg_debug(modName, "Missing slit width - set to 0.5 mm");
        exSlit->width = .5;
      }
    }
    else {
      deleteExtractionSlit(exSlit);
    }
    tSlitHolder = tSlitHolder->next;
      
  }
  return(VM_TRUE);
}

/*
   Return number of spectra rows in a Extraction Slit. Since an Extraction
   Slit is a linked list of Extraction Slits, more than one slit can be
   contained in the list. The sum of all rows in all slits is returned.
*/
int numRowsInExtSlits(VimosExtractionSlit *slit)
{
  const char modName[] = "numRowsInExtSlits";
  int numRows;

  if (slit == NULL)
  {
    pilErrno = 1;
    cpl_msg_error(modName, "NULL imput slit");
    return 0;
  }

  /* initialize */
  numRows = 0;
  
  /* traverse list and count */
  while (slit) {
    numRows += slit->numRows;
    slit = slit->next;
  }
  
  return(numRows);
}


int numSlitsInExtTable(VimosExtractionTable *exTab)
{
  const char modName[] = "numSlitsInExtTable";
  int numSlits;
  VimosExtractionSlit *slit;
  
  if (exTab == NULL)
  {
    pilErrno = 1;
    cpl_msg_error(modName, "NULL input extraction Table");
    return 0;
  }

  /* initialize */
  numSlits = 0;
  slit = exTab->slits;
  
  /* traverse list and count */
  while (slit) {
    numSlits++;
    slit = slit->next;
  }
  
  return(numSlits);
  
}



/*
 * Slit closest to mask center is found on given Extraction Table.
 * Shutter positions are ignored - slit might be not illuminated.
 * Mask center is assumed to have coordinates (0,0).
 */

VimosExtractionSlit *slitClosestToCenter(VimosExtractionTable *extractionTable)
{
  const char modName[] = "slitClosestToCenter";

  double               distance, newdist;
  int                  crow;
  VimosExtractionSlit *slit;
  VimosExtractionSlit *closestSlit;

  if (extractionTable == NULL) {
    pilErrno = 1;
    cpl_msg_error(modName, "NULL input extraction Table");
    return 0;
  }

  closestSlit = slit = extractionTable->slits;
  crow = slit->numRows / 2;
  distance = slit->maskX->data[crow] * slit->maskX->data[crow]
           + slit->maskY->data[crow] * slit->maskY->data[crow];
  slit = slit->next;

  while (slit) {
    crow = slit->numRows / 2;
    newdist = slit->maskX->data[crow] * slit->maskX->data[crow]
            + slit->maskY->data[crow] * slit->maskY->data[crow];
    if (newdist < distance) {
      distance = newdist;
      closestSlit = slit;
    }
    slit = slit->next;
  }

  return(closestSlit);

}


VimosBool slitLongOrShort(VimosExtractionSlit *slit, float tolerance)
{
  int i;
  float yStart;
  float maxDrift;
  float drift;
  
  /* pathetic case */
  if (slit->numRows < 2) {
    return(VM_FALSE);
  }
  
  yStart = slit->ccdY->data[0];
  maxDrift = 0.0;
  
  for (i = 1; i < slit->numRows; i++) {
    drift = fabs(slit->ccdY->data[i] - yStart);
    if (drift > maxDrift) {
      maxDrift = drift;
    }
  }
  
  if (maxDrift > tolerance) {
    return(VM_TRUE);
  }
  else {
    return(VM_FALSE);
  }
  
}

VimosBool readFitsExtractionTable(VimosExtractionTable *extTable, fitsfile
				  *fptr)
{
  const char modName[] = "readFitsExtractionTable";
  int    i, j, k;
  int    nCols;
  int    nRows;
  int    null;
  int   *crvPolCol;
  int   *invDisCol;
  int    slitCol, yCol, ccdXCol, ccdYCol, maskXCol, maskYCol, specNoCol;
  int    zeroXCol, zeroYCol, disQualCol,crvPolRmsCol,invDisRmsCol;
  int    slitLen;
  int    newSlitNo;
  int    slitNo;
  int    IFUslitNoCol;
  int    IFUfibNoCol;
  /*ALEX: added IFUfibPeakX and IFUfibTrans columns */
  int    IFUfibPeakXCol;
  int    IFUfibTransCol;

  float  lambda0;
  double dValue;
  int    crvOrder, invOrder;

  char   colName[80];
  char   comment[80];
  int    nfound, ii;
  char **ttype;
  int    status = 0;

  VimosExtractionSlit *exSlit;
  VimosExtractionSlit *lastSlit;

  /* validate input */
  if (extTable == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  /* validate input */
  if (strcmp(extTable->name, VM_EXT)) {
    cpl_msg_error(modName, "Invalid input table");
    return(VM_FALSE);
  }
  
  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, "EXR", 0, &status)) {
    cpl_msg_error(modName, "fits_movnam_hdu returned error code %d", status);
    return(VM_FALSE);
  }

  extTable->fptr = fptr;
  
  /* read Table */
  if (!readDescsFromFitsTable(&(extTable->descs), extTable->fptr)) {
    cpl_msg_error(modName, "readDescsFromFitsTable returned an error");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(extTable->descs, "ESO PRO CRV POL ORD", &crvOrder,
			comment)) {
    cpl_msg_error(modName, "Cannot read descriptor %s", "ESO PRO CRV POL ORD");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(extTable->descs, "ESO PRO IDS REL ORD", &invOrder,
		        comment)) {
    cpl_msg_error(modName, "Cannot read descriptor %s", "ESO PRO IDS REL ORD");
    return(VM_FALSE);
  }

  if (!readDoubleDescriptor(extTable->descs, "ESO PRO WLEN CEN", &dValue, 
			   comment)) {
    cpl_msg_error(modName, "Cannot read descriptor %s", "ESO PRO WLEN CEN");
    return(VM_FALSE);
  }

  lambda0 = (float) dValue;
  readIntDescriptor(extTable->descs, "TFIELDS", &nCols, comment);
  readIntDescriptor(extTable->descs, "NAXIS2", &nRows, comment);

  if (!readIntDescriptor(extTable->descs, "TFIELDS", &nCols, comment)) {
    cpl_msg_error(modName, "Cannot read descriptor %s", "TFIELDS");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(extTable->descs, "NAXIS2", &nRows, comment)) {
    cpl_msg_error(modName,"Cannot read descriptor %s", "NAXIS2");
    return(VM_FALSE);
  }

  /* allocate space for the column labels */
  ttype = (char **) cpl_malloc(nCols*sizeof(char *));
  for (ii = 0; ii < nCols; ii++) {
    ttype[ii] = (char *) cpl_malloc(FLEN_VALUE*sizeof(char));

    /* check if space was allocated */
    if (ttype[ii] == NULL) {
      cpl_msg_error(modName, "Allocation Error");
      return(VM_FALSE);
    }
  }
      /* read the column names from the TTYPEn keywords */
  if (fits_read_keys_str(extTable->fptr, "TTYPE", 1, nCols, ttype, &nfound, 
			&status)) {
    cpl_msg_error(modName, "fits_read_keys_str returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"SLIT",&slitCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"Y",&yCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"CCD_X",&ccdXCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"CCD_Y",&ccdYCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"MASK_X",&maskXCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"MASK_Y",&maskYCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"SPEC_NO",&specNoCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"ZERO_X",&zeroXCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"ZERO_Y",&zeroYCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"CRVPOL_RMS",&crvPolRmsCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"INVDIS_RMS",&invDisRmsCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }

    /* NOTE FOR IFU: these added to ease IFU reduction.  */
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"IFUSLIT_NO",&IFUslitNoCol,
                     &status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"IFUFIB_NO",&IFUfibNoCol,
		     &status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }

  /*ALEX: read IFUfibPeakX and IFUfibTrans columns */
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"IFUFIBPEAKX",&IFUfibPeakXCol,
		     &status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }

  if (fits_get_colnum(extTable->fptr,CASEINSEN,"IFUFIBTRANS",&IFUfibTransCol,
		     &status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }

  crvPolCol = (int *) cpl_malloc((crvOrder+1)*sizeof(int));
  /* check if space was allocated */
  if (crvPolCol == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(VM_FALSE);
  }

  for (i = 0; i <= crvOrder; i++) {
    sprintf(colName,"CRV_POL_%d",i);
    if (fits_get_colnum(extTable->fptr,CASEINSEN,colName,&(crvPolCol[i]),
		       &status)) {
      cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
      return(VM_FALSE);
    }
  }

  invDisCol = (int *) cpl_malloc((invOrder+1)*sizeof(int));
  /* check if space was allocated */
  if (invDisCol == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(VM_FALSE);
  } 

  for (i = 0; i <= invOrder; i++) {
    sprintf(colName,"INV_DIS_%d",i);
    if (fits_get_colnum(extTable->fptr,CASEINSEN,colName,&(invDisCol[i]),
		       &status)) {
      cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
      return(VM_FALSE);
    }
  }
  if (fits_get_colnum(extTable->fptr,CASEINSEN,"DIS_QUAL",&disQualCol,&status)) {
    cpl_msg_error(modName, "fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }

   if (fits_read_col_int(extTable->fptr,slitCol,1,1,1,null,&newSlitNo,&null,
			&status)) {
     cpl_msg_error(modName, "fits_read_col_int returned error code %d", status);
     return(VM_FALSE);
   } 

  lastSlit = NULL;

  i = 0;
  while (i < nRows) {
    exSlit = newExtractionSlit();
    if (exSlit == NULL) {
      cpl_msg_error(modName, "The function newExtractionSlit has returned NULL");
      return(VM_FALSE);
    }

    slitLen = 0;
    
    slitNo = newSlitNo;
    while ( (i+slitLen < nRows) && (newSlitNo == slitNo) ) {
      if (fits_read_col_int(extTable->fptr,slitCol,i+slitLen+1,1,1,null,
			   &newSlitNo,&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_int returned error %d", status);
        return(VM_FALSE);
      }

      /*ALEX: add this "if" also here, otherwise it misses the first
       fiber and duplcates the 1600th one*/
      if (newSlitNo == slitNo) {
	if (fits_read_col_int(extTable->fptr,IFUslitNoCol,i+slitLen+1,1,1,null,
			      &(exSlit->IFUslitNo),&null,&status)) {
	  cpl_msg_error(modName, "fits_read_col_int returned error %d", status);
	  return(VM_FALSE);
	} 
	if (fits_read_col_int(extTable->fptr,IFUfibNoCol,i+slitLen+1,1,1,null,
			      &(exSlit->IFUfibNo),&null,&status)) {
	  cpl_msg_error(modName, "fits_read_col_int returned error %d", status);
	  return(VM_FALSE);
	} 
	
	/*ALEX: read IFUfibPeakX and IFUfibTrans columns */
	if (fits_read_col_flt(extTable->fptr,IFUfibPeakXCol,i+slitLen+1,1,1,null,
			      &(exSlit->IFUfibPeakX),&null,&status)) {
	  cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
	  return(VM_FALSE);
	}  

	if (fits_read_col_flt(extTable->fptr,IFUfibTransCol,i+slitLen+1,1,1,null,
			      &(exSlit->IFUfibTrans),&null,&status)) {
	  cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
	  return(VM_FALSE);
	}  
      }/*ALEX: end "if"*/

      slitLen++;
    }

    if ( i+slitLen < nRows) {
      slitLen--;
    }
    exSlit->slitNo = slitNo;
    exSlit->numRows = slitLen;
    exSlit->y = newIntArray(slitLen); 
    if (exSlit->y == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newIntArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->ccdX = newFloatArray(slitLen);
    if (exSlit->ccdX == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newFloatArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->ccdY = newFloatArray(slitLen);
    if (exSlit->ccdY == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newFloatArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->maskX = newFloatArray(slitLen);
    if (exSlit->maskX == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newFloatArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->maskY = newFloatArray(slitLen);
    if (exSlit->maskY == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newFloatArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->numSpec = newIntArray(slitLen);
    if (exSlit->numSpec == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newIntArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->zeroX = newFloatArray(slitLen);
    if (exSlit->zeroX == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newFloatArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->zeroY = newFloatArray(slitLen);
    if (exSlit->zeroY == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newFloatArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->crvPol = 
      (VimosDistModel1D **)cpl_malloc(slitLen*sizeof(VimosDistModel1D *));
    /* check if space was allocated */
    if (exSlit->crvPol == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "Allocation Error");
      return(VM_FALSE);
    }
    exSlit->crvPolRms = newFloatArray(slitLen);
    if (exSlit->crvPolRms == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newFloatArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->invDis = 
      (VimosDistModel1D **)cpl_malloc(slitLen*sizeof(VimosDistModel1D *));
    if (exSlit->invDis == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "Allocation Error");
      return(VM_FALSE);
    }
    exSlit->invDisRms = newFloatArray(slitLen);
    if (exSlit->invDisRms == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newFloatArray has return NULL");
      return(VM_FALSE);
    }   
    exSlit->invDisQuality = newIntArray(slitLen);
    if (exSlit->invDisQuality == NULL) {
      deleteExtractionSlit(exSlit);
      cpl_msg_error(modName, "The function newIntArray has return NULL");
      return(VM_FALSE);
    }   

    for (j = 0; j < slitLen; j++) {
      if (fits_read_col_int(extTable->fptr,yCol,i+j+1,1,1,null,
			    &(exSlit->y->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_int returned error %d", status);
        return(VM_FALSE);
      }  
      if (fits_read_col_flt(extTable->fptr,ccdXCol,i+j+1,1,1,null,
			    &(exSlit->ccdX->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
        return(VM_FALSE);
      }  
      if (fits_read_col_flt(extTable->fptr,ccdYCol,i+j+1,1,1,null,
			    &(exSlit->ccdY->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
        return(VM_FALSE);
      }   
      if (fits_read_col_flt(extTable->fptr,maskXCol,i+j+1,1,1,null,
			    &(exSlit->maskX->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
        return(VM_FALSE);
      }   
      if (fits_read_col_flt(extTable->fptr,maskYCol,i+j+1,1,1,null,
			    &(exSlit->maskY->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
        return(VM_FALSE);
      }  
      if (fits_read_col_int(extTable->fptr,specNoCol,i+j+1,1,1,null,
			    &(exSlit->numSpec->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_int returned error %d", status);
        return(VM_FALSE);
      }  
      if (fits_read_col_flt(extTable->fptr,zeroXCol,i+j+1,1,1,null,
			    &(exSlit->zeroX->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
        return(VM_FALSE);
      }   
      if (fits_read_col_flt(extTable->fptr,zeroYCol,i+j+1,1,1,null,
			    &(exSlit->zeroY->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
        return(VM_FALSE);
      }  
      exSlit->crvPol[j] = newDistModel1D(crvOrder);  
      if (exSlit->crvPol[j] == NULL) {
	deleteExtractionSlit(exSlit);
	cpl_msg_error(modName, "The function newDistModel1D has returned NULL");
	return(VM_FALSE);
      }

      for (k = 0; k <= crvOrder; k++) {
	if (fits_read_col_dbl(extTable->fptr,crvPolCol[k],i+j+1,1,1,null,
			      &(exSlit->crvPol[j]->coefs[k]),&null,&status)) {
	  cpl_msg_error(modName, "fits_read_col_dbl returned error %d", status);
          return(VM_FALSE);
	}   
      }
      exSlit->crvPol[j]->offset = exSlit->ccdY->data[j];
      if (fits_read_col_flt(extTable->fptr,crvPolRmsCol,i+j+1,1,1,null,
			    &(exSlit->crvPolRms->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
        return(VM_FALSE);
      }  
      
      exSlit->invDis[j] = newDistModel1D(invOrder);
      if (exSlit->invDis[j] == NULL) {
	deleteExtractionSlit(exSlit);
	cpl_msg_error(modName, "The function newDistModel1D has returned NULL");
	return(VM_FALSE);
      }
      for (k = 0; k <= invOrder; k++) {
	if (fits_read_col_dbl(extTable->fptr,invDisCol[k],i+j+1,1,1,null,
			      &(exSlit->invDis[j]->coefs[k]),&null,&status)) { 
	  cpl_msg_error(modName, "fits_read_col_dbl returned error %d", status);
          return(VM_FALSE);
	}    
	
	exSlit->invDis[j]->offset = lambda0;
      }
      if (fits_read_col_flt(extTable->fptr,invDisRmsCol,i+j+1,1,1,null,
			    &(exSlit->invDisRms->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_flt returned error %d", status);
        return(VM_FALSE);
      }  
      if (fits_read_col_int(extTable->fptr,disQualCol,i+j+1,1,1,null,
			    &(exSlit->invDisQuality->data[j]),&null,&status)) {
        cpl_msg_error(modName, "fits_read_col_int returned error %d", status);
        return(VM_FALSE);
      } 

    } 
    i += slitLen;
    if (lastSlit == NULL) {
      extTable->slits = exSlit;
    }
    else {
      lastSlit->next = exSlit;
      exSlit->prev = lastSlit;
    }
    lastSlit = exSlit;
    
  }
   
  cpl_free(crvPolCol);
  cpl_free(invDisCol);
  for (ii = 0; ii < 3; ii++){     /* free the memory for the column labels */
    cpl_free( ttype[ii] );
  }

  return (VM_TRUE);
}


VimosBool writeFitsExtractionTable(VimosExtractionTable *extTable, fitsfile
				   *fptr)
{
  int i, j, k;
  int nRows;
  int rowNum;
  int slitLen;
  
  int crvOrder, invOrder;

  char colName[80];
  int status, nbytes;
  char *ttype[84], *tform[84], comment[80];

  

  VimosExtractionSlit *exSlit;
  
  /* validate input */
  if (extTable == NULL) {
    cpl_msg_error("writeFitsExtractionTable","NULL input table");
    return (VM_FALSE);
  }

  /* validate input */
  if ( strcmp(extTable->name, VM_EXT) ) {
    cpl_msg_error("writeFitsExtractionTable","Invalid input table");
    return (VM_FALSE);
  }
  
  extTable->fptr = fptr;
  if (!readIntDescriptor(extTable->descs, "ESO PRO CRV POL ORD", &crvOrder,
			 comment)) {
    cpl_msg_error("writeFitsExtractionTable","The function readIntDescriptor has returned an error");
    return(VM_FALSE);
  }
  if (!readIntDescriptor(extTable->descs, "ESO PRO IDS REL ORD", &invOrder,
			 comment)) {
    cpl_msg_error("writeFitsExtractionTable","The function readIntDescriptor has returned an error");
    return(VM_FALSE);
  }

  exSlit = extTable->slits;

  nRows = numRowsInExtSlits(exSlit);

  status = 0;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "EXR", 0, &status) ) {
    if(fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_error("writeFitsExtractionTable","The function fits_delete_hdu has returned an error (code %d)", status);
      return (VM_FALSE);
    }
  } else {
    status = 0;
  }

  for (i = 0; i <= 15+crvOrder+1+invOrder+1; i++){  
    /* allocate space for the column labels */
    ttype[i] = (char *) cpl_malloc(FLEN_VALUE);  /* max label length = 69 */
    /* check if space was allocated */
    if (ttype[i] == NULL) {
      cpl_msg_error("writeFitsExtractionTable","Allocation Error");
      return(VM_FALSE);
    }
    tform[i] = (char *) cpl_malloc(FLEN_VALUE);
    /* check if space was allocated */
    if (tform[i] == NULL) {
      cpl_msg_error("writeFitsExtractionTable","Allocation Error");
      return(VM_FALSE);
    }
  }

  ttype[0] = "SLIT";
  tform[0] = "1J";
  ttype[1]="IFUSLIT_NO";
  tform[1] = "1J";
  ttype[2]="IFUFIB_NO"; 
  tform[2] = "1J";
  ttype[3] = "Y";
  tform[3] = "1J";
  ttype[4] = "CCD_X";
  tform[4] = "1E";
  ttype[5]= "CCD_Y";
  tform[5] = "1E";
  ttype[6]= "MASK_X";
  tform[6] = "1E";
  ttype[7]="MASK_Y";
  tform[7] = "1E";
  ttype[8]="SPEC_NO";
  tform[8] = "1J";
  ttype[9]="ZERO_X";
  tform[9] = "1E";
  ttype[10]="ZERO_Y";
  tform[10] = "1E";

  for (i = 0; i <= crvOrder; i++) {
    sprintf(colName,"CRV_POL_%d",i);
    sprintf(ttype[10+i+1], "%s", colName);
    tform[10+i+1] = "1D";
  }
  ttype[11+crvOrder+1]="CRVPOL_RMS";
  tform[11+crvOrder+1] = "1E";

  for (i = 0; i <= invOrder; i++) {
    sprintf(colName,"INV_DIS_%d",i);
    sprintf(ttype[11+crvOrder+1+i+1], "%s", colName); 
    tform[11+crvOrder+1+i+1] = "1D";
  }
  ttype[12+crvOrder+1+invOrder+1]="INVDIS_RMS";
  tform[12+crvOrder+1+invOrder+1] = "1E";

  ttype[13+crvOrder+1+invOrder+1] = "DIS_QUAL";
  tform[13+crvOrder+1+invOrder+1] = "1J";

  /*ALEX PROVA */
  ttype[14+crvOrder+1+invOrder+1] = "IFUFIBPEAKX";
  tform[14+crvOrder+1+invOrder+1] = "1E";
  ttype[15+crvOrder+1+invOrder+1] = "IFUFIBTRANS";
  tform[15+crvOrder+1+invOrder+1] = "1E";

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl(fptr, BINARY_TBL,0, 16+crvOrder+1+invOrder+1, ttype, 
		      tform,NULL, "EXR", &status)) { 
    cpl_msg_error("writeFitsExtractionTable","The function fits_create_tbl has returned an error (code %d)", status);
    return (VM_FALSE); 
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "EXR", 0, &status)) {
    cpl_msg_error("writeFitsExtractionTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  /* write the table descriptors to the Fits file */
  
  if (fits_read_key (extTable->fptr,TINT,"NAXIS1",&nbytes,NULL,&status)) {
    cpl_msg_error("writeFitsExtractionTable","The function fits_read_key has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(extTable->descs), "NAXIS1", nbytes, "")) {
    cpl_msg_error("writeFitsExtractionTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(extTable->descs), "NAXIS2", nRows, "")) {
    cpl_msg_error("writeFitsExtractionTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }

  /*ALEX: added +2 columns: 11+crvOrder+1+invOrder+1+1+2 */
  if (!writeIntDescriptor(&(extTable->descs), "TFIELDS",
			  16+crvOrder+1+invOrder+1,"")) {
    cpl_msg_error("writeFitsExtractionTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeDescsToFitsTable(extTable->descs, extTable->fptr)) {
    cpl_msg_error("writeFitsExtractionTable","The function writeDescsToFitsTable has returned an error");
    return (VM_FALSE);
  }

  exSlit = extTable->slits;
  rowNum = 1;
  while (exSlit) {
    slitLen = exSlit->numRows;
    
    for (j = 0; j < slitLen; j++) {
      if (fits_write_col_int(extTable->fptr,1,rowNum,1,1,&(exSlit->slitNo),
			    &status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_int has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_int(extTable->fptr,2,rowNum,1,1,&(exSlit->IFUslitNo),
			     &status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_int has returned an error (code %d)", status);
        return(VM_FALSE); 
      }
      if (fits_write_col_int(extTable->fptr,3,rowNum,1,1,&(exSlit->IFUfibNo),
			     &status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_int has returned an error (code %d)", status);
        return(VM_FALSE); 
      }
      if (fits_write_col_int(extTable->fptr,4,rowNum,1,1,&(exSlit->y->data[j]),
			    &status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_flt(extTable->fptr,5,rowNum,1,1,
			    &(exSlit->ccdX->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_flt(extTable->fptr,6,rowNum,1,1,
			    &(exSlit->ccdY->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_flt(extTable->fptr,7,rowNum,1,1,
			    &(exSlit->maskX->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_flt(extTable->fptr,8,rowNum,1,1,
			    &(exSlit->maskY->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_int(extTable->fptr,9,rowNum,1,1,
			    &(exSlit->numSpec->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_flt(extTable->fptr,10,rowNum,1,1,
			    &(exSlit->zeroX->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_flt(extTable->fptr,11,rowNum,1,1,
			    &(exSlit->zeroY->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
  
      for (k = 0; k <= crvOrder; k++) {
	if (fits_write_col_dbl(extTable->fptr,11+k+1,rowNum,1,1,
			      &(exSlit->crvPol[j]->coefs[k]),&status)) {
          cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_dbl has returned an error (code %d)", 
                       status);
	  return(VM_FALSE);
	} 
      }
      if (fits_write_col_flt(extTable->fptr,12+crvOrder+1,rowNum,1,1,
			     &(exSlit->crvPolRms->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }

      for (k = 0; k <= invOrder; k++) {
	if (fits_write_col_dbl(extTable->fptr,12+k+1+crvOrder+1,rowNum,1,1,
			      &(exSlit->invDis[j]->coefs[k]),&status)) {
          cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_dbl has returned an error (code %d)", 
                       status);
	  return(VM_FALSE);
	}
      }
      if (fits_write_col_flt(extTable->fptr,13+crvOrder+1+invOrder+1,
			     rowNum,1,1,&(exSlit->invDisRms->data[j]),
			     &status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }
      if (fits_write_col_int(extTable->fptr,14+crvOrder+1+invOrder+1,
		      rowNum,1,1,&(exSlit->invDisQuality->data[j]),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }

      /*ALEX: added IFUfibPeakX and IFUfibTrans columns */
      if (fits_write_col_flt(extTable->fptr,15+crvOrder+1+invOrder+1,
			     rowNum,1,1,
			    &(exSlit->IFUfibPeakX),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }

      if (fits_write_col_flt(extTable->fptr,16+crvOrder+1+invOrder+1,
			     rowNum,1,1,
			    &(exSlit->IFUfibTrans),&status)) {
        cpl_msg_error("writeFitsExtractionTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return(VM_FALSE);
      }


      rowNum++;
    }
    
    exSlit = exSlit->next;
  }
   
  return(VM_TRUE);
}
