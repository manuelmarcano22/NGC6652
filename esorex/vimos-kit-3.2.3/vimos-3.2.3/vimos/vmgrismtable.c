/* $Id: vmgrismtable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmtable.h"
#include "vmmatrix.h"
#include "vmgrismtable.h"
#include "cpl.h"


/*
  Allocate a new Grism Table. This is a VimosTable with name VM_GRS (which
  should have value "GRS").
  No descriptors or columns are allocated (but this may change in a later
  version).
*/
VimosTable *newGrismTable()
{
  VimosTable *newTab;
  
  /* allocate new VimosTable */
  newTab = newTable();

  /* check if space was allocated */
  if (newTab == NULL) {
    cpl_msg_error("newGrismTable","The function newTable has returned NULL");
    return(NULL);
  }
  
  /* copy "GRS" into name of table */
  strcpy(newTab->name, VM_GRS);
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_GRS, "");
  if (newTab->descs == NULL) {
    /* cleanup */
    cpl_free(newTab);
    cpl_msg_error("newGrismTable","The function newStringDescriptor has returned NULL");
    return(NULL);
  }
  
  /* should create all the descriptors and columns here ?*/
  
  /* return address of new Grism Table */
  return(newTab);
  
}

/*
  Delete a Grism Table. This is just an esthetic wrapper for 
  deleteTable(grsTable)
*/
void deleteGrismTable(VimosTable *grsTable)
{
  /* just a wrapper for deleteTabel() */
  deleteTable(grsTable);
}


/* read the info on sky lines from the Descriptors of a table */
VimosBool readSkyLines(VimosDescriptor *desc, int *numSkyLines, 
                       VimosFloatArray **skyLines)
{
  int i;
  VimosBool rdOK;
  double dValue;
  char  descName[80];
  char comment[80];

  if (desc == NULL) {
    *numSkyLines = 0;
    *skyLines = NULL;
    cpl_msg_error("readSkyLines","NULL input pointer");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(desc, "ESO PRO SKY NO", numSkyLines, comment)) {
    cpl_msg_error("readSkyLines","The function readIntDescriptor has returned an error");
    return(VM_FALSE);
  }

  *skyLines = newFloatArray(*numSkyLines);
  if (*skyLines == NULL) {
    cpl_msg_error("readSkyLines","The function newFloatArray has returned an error");
    return(VM_FALSE);
  }
  
  for (i = 0; i < *numSkyLines; i++) {
    sprintf(descName,"ESO PRO SKY WLEN%d",i+1);
    rdOK = readDoubleDescriptor(desc,descName,&dValue,comment);
    (*skyLines)->data[i] = (float) dValue;
    if (!rdOK) {
      deleteFloatArray(*skyLines);
      *skyLines = NULL;
      numSkyLines = 0;
      cpl_msg_error("readSkyLines","The function readDoubleDescriptor has returned an error");
      return(VM_FALSE);
    }
  }

  return(VM_TRUE);
}


VimosBool readFitsGrismTable(VimosTable *grsTable, fitsfile *fptr)
{
  int status;
  
/* validate input */
  if (grsTable == NULL) {
    cpl_msg_error("readFitsGrismTable","NULL input table");
    return(VM_FALSE);
  }
  
  status = 0;

  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, "GRS", 0, &status)) {
    cpl_msg_error("readFitsGrismTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return(VM_FALSE);
  }

  grsTable->fptr = fptr;
  /* read Fits Table */
  if (!readDescsFromFitsTable(&(grsTable->descs), grsTable->fptr)) {
    cpl_msg_error("readFitsGrismTable","The function readDescsFromFitsTable has returned an error");
    return(VM_FALSE);
  }

  return(VM_TRUE);
}



VimosBool writeFitsGrismTable(VimosTable *grsTable, fitsfile *fptr)
{
  int status;

  /* validate input */
  if (grsTable == NULL) {
    cpl_msg_error("writeFitsGrismTable","NULL input table");
    return(VM_FALSE);
  }
  
  /* check table name, should be "GRS" */
  if ( strcmp(grsTable->name, VM_GRS) ) {
    cpl_msg_error("writeFitsGrismTable","Invalid input table");
    return(VM_FALSE);
  }

  status = 0;
  grsTable->fptr = fptr;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "GRS", 0, &status)) {
    if (fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_error("writeFitsGrismTable","The function fits_delete_hdu has returned an error (code %d)", status);
      return(VM_FALSE);
    }
  } else {
    status = 0;
  }

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl(fptr,BINARY_TBL,0,0,NULL,NULL,NULL,"GRS",&status)) {
    cpl_msg_error("writeFitsGrismTable","The function fits_create_tbl has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "GRS", 0, &status)) {
    cpl_msg_error("writeFitsGrismTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return(VM_FALSE);
  }

  /* write the table to the Fits file */
  if (!writeDescsToFitsTable(grsTable->descs, grsTable->fptr)) {
    cpl_msg_error("writeFitsGrismTable","The function writeDescsToFitsTable has returned an error");
    return(VM_FALSE);
  }

  return(VM_TRUE);
}
