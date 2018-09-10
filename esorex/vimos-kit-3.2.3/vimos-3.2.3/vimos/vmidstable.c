/* $Id: vmidstable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <piltranslator.h>

#include "vmtable.h"
#include "vmmatrix.h"
#include "vmidstable.h"
#include "cpl.h"


/*
  Allocate a new Ids Table. This is a VimosTable with name VM_IDS (which
  should have value "IDS").
  No descriptors or columns are allocated (but this may change in a later
  version).
*/

VimosTable *newIdsTable()
{
  char modName[] = "newIdsTable";
  VimosTable *newTab;
  
  /* allocate new VimosTable */
  newTab = newTable();

  /* check if space was allocated */
  if (newTab == NULL) {
    cpl_msg_error(modName, "The function newTable has returned NULL");
    return(NULL);
  }
  
  /* copy "IDS" into name of table */
  strcpy(newTab->name, VM_IDS);
  newTab->descs = newStringDescriptor(pilTrnGetKeyword("Table"),
                                      VM_IDS,
                                      pilTrnGetComment("Table"));
  if (newTab->descs == NULL) {
    /* cleanup */
    cpl_free(newTab);
    cpl_msg_error(modName, "The function newStringDescriptor has returned NULL");
    return(NULL);
  }
  
  /* should create all the descriptors and columns here ?*/

  return(newTab);
  
}

/*
  Delete a IDS Table. This is just an esthetic wrapper for 
  deleteTable(idsTable)
*/
void deleteIdsTable(VimosTable *idsTable)
{
  /* just a wrapper for deleteTabel() */
  deleteTable(idsTable);
}

VimosBool readFitsIdsTable(VimosTable *idsTable, fitsfile *fptr)
{

  char modName[] = "readFitsIdsTable";
  int status;
  
/* validate input */
  if (idsTable == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  status = 0;

  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, "IDS", 0, &status)) {
    cpl_msg_error(modName, "The function fits_movnam_hdu has returned an error (code %d)", status);
    return(VM_FALSE);
  }

  idsTable->fptr = fptr;
  /* read Fits Table */
  if (!readDescsFromFitsTable(&(idsTable->descs), idsTable->fptr)) {
    cpl_msg_error(modName, "The function readDescsFromFitsTable has returned an error");
    return(VM_FALSE);
  }

  return(VM_TRUE);
}



VimosBool writeFitsIdsTable(VimosTable *idsTable, fitsfile *fptr)
{

  char modName[] = "writeFitsIdsTable";
  int status;

  /* validate input */
  if (idsTable == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  /* check table name, should be "IDS" */
  if ( strcmp(idsTable->name, VM_IDS) ) {
    cpl_msg_error(modName, "Invalid input table");
    return(VM_FALSE);
  }

  status = 0;
  idsTable->fptr = fptr;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "IDS", 0, &status)) {
    if (fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_error(modName, "The function fits_delete_hdu has returned an error (code %d)", status);
      return(VM_FALSE);
    }
  } else {
    status = 0;
  }

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl(fptr,BINARY_TBL,0,0,NULL,NULL,NULL,"IDS",&status)) {
    cpl_msg_error(modName, "The function fits_create_tbl has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "IDS", 0, &status)) {
    cpl_msg_error(modName, "The function fits_movnam_hdu has returned an error (code %d)", status);
    return(VM_FALSE);
  }

  /* write the table to the Fits file */
  if (!writeDescsToFitsTable(idsTable->descs, idsTable->fptr)) {
    cpl_msg_error(modName, "The function writeDescsToFitsTable has returned an error");
    return(VM_FALSE);
  }

  return(VM_TRUE);
}
