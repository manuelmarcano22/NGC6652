/* $Id: vmsphottable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include "vmsphottable.h"
#include "cpl.h"


/*
  Allocate a new SpectroPhotometric Table. This is a VimosTable with 
  name VM_SPHOT (which should have value "SPH").
  No descriptors or columns are allocated (but this may change in a later
  version).
*/
VimosTable *newSphotTable()
{
  VimosTable *newTab;
  
  /* allocate new VimosTable */
  newTab = newTable();

  /* check if new VimosTable was allocated */
  if (newTab == NULL) {
    cpl_msg_error("newSphotTable","The function newTable has returned NULL");
    return(NULL);
  }
  
  /* copy "SPHOT" into name of table */
  strcpy(newTab->name, VM_SPHOT);
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_SPHOT, "");
  if (newTab->descs == NULL) {
    /* cleanup */
    cpl_free(newTab);
    cpl_msg_error("newSphotTable","The function newStringDescriptor has returned NULL");
    return(NULL);
  }

  /* should create all the descriptors and columns here ? */
  
  /* should create all the descriptors and columns here ?*/
  
  /* return address of new Line Catalog */
  return(newTab);
  
}

/*
  Delete a SpectroPhotometric Table. This is just an esthetic wrapper for 
  deleteTable(sTab)
*/
void deleteSphotTable(VimosTable *sTab)
{
  /* just a wrapper for deleteTable() */
  deleteTable(sTab);
}

VimosBool readFitsSphotTable(VimosTable *sphotTable, fitsfile *fptr)
{
  int    i;
  int    nCols;
  int    nRows;
  int    null;
  int    lambdaColNo, magColNo, deltaLColNo;
  int    status, nfound, ii;
  char **ttype;
  char   comment[80];

  VimosColumn *lambdaCol;
  VimosColumn *magCol;
  VimosColumn *deltaLCol;

  char         modName[] = "readFitsSphotTable";

  status = 0;

  /* validate input */
  if (sphotTable == NULL) {
    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  /* validate input */
  if ( strcmp(sphotTable->name, VM_SPHOT) ) {
    cpl_msg_debug(modName, "Invalid input table");
    return(VM_FALSE);
  }
  
  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, "SPH", 0, &status)) {
    cpl_msg_debug(modName, 
    "The function fits_movnam_hdu returned error code %d", status);
    return(VM_FALSE);
  }
  
  sphotTable->fptr = fptr;
  
  /* read Table */
  if (!readDescsFromFitsTable(&(sphotTable->descs), sphotTable->fptr)) {
    cpl_msg_debug(modName, "Function readDescsFromFitsTable returned an error");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(sphotTable->descs, "TFIELDS", &nCols, comment)) {
    cpl_msg_debug(modName, "The function readIntDescriptor returned an error");
    return(VM_FALSE);
  }
  sphotTable->numColumns=nCols;

  if (!readIntDescriptor(sphotTable->descs, "NAXIS2", &nRows, comment)) {
    cpl_msg_debug(modName, "The function readIntDescriptor returned an error");
    return(VM_FALSE);
  }

  /* allocate space for the column labels */
  ttype = (char **) cpl_malloc(nCols*sizeof(char *));
  for (ii = 0; ii < nCols; ii++) {
    ttype[ii] = (char *) cpl_malloc(FLEN_VALUE*sizeof(char));
    /* check if space was allocated */
    if (ttype[ii] == NULL) {
      cpl_msg_debug(modName, "Allocation Error");
      return(VM_FALSE);
    }
  }

  /* read the column names from the TTYPEn keywords */
  if (fits_read_keys_str(sphotTable->fptr, "TTYPE", 1, nCols, ttype, &nfound, 
			 &status)) {
    cpl_msg_debug(modName, 
    "Function fits_read_keys_str returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(sphotTable->fptr,CASEINSEN,"LAMBDA",&lambdaColNo,&status)) {
    cpl_msg_debug(modName, 
    "The function fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(sphotTable->fptr,CASEINSEN,"MAG",&magColNo,&status)) {
    cpl_msg_debug(modName, 
    "The function fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(sphotTable->fptr,CASEINSEN,"DELTA_LAMBDA",&deltaLColNo,&status)) {
    cpl_msg_debug(modName, 
    "The function fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  
  /* create column for LAMBDA */
  lambdaCol = newFloatColumn(nRows, "LAMBDA");
  if (lambdaCol == NULL) {
    cpl_msg_debug(modName, "The function newFloatColumn has returned NULL");
    return(VM_FALSE);
  }

  /* put column in Table */
  sphotTable->cols = lambdaCol;
  
  /* create column for MAG */
  magCol = newFloatColumn(nRows, "MAG");
  if (magCol == NULL) {
    cpl_msg_debug(modName, "The function newStringColumn has returned NULL");
    return(VM_FALSE);
  }

  /* link column into table */
  lambdaCol->next = magCol;
  magCol->prev = lambdaCol;
  
  /* create column for FLUX */
  deltaLCol = newFloatColumn(nRows, "DELTA_LAMBDA");
  if (deltaLCol == NULL) {
    cpl_msg_debug(modName, "The function newFloatColumn has returned NULL");
    return(VM_FALSE);
  }
  /* link column into table */
  magCol->next = deltaLCol;
  deltaLCol->prev = magCol;
  
  for (i = 1; i <= nRows; i++) {
      if (fits_read_col_flt(sphotTable->fptr,lambdaColNo,i,1,1,null,
			    &(lambdaCol->colValue->fArray[i-1]),&null,&status)){
	cpl_msg_debug(modName, 
        "Function fits_read_col_flt returned error code %d", status);
	return(VM_FALSE);
      } 
      
      if (fits_read_col_flt(sphotTable->fptr,magColNo,i,1,1,null,
			    &(magCol->colValue->fArray[i-1]),&null,&status)) {
	cpl_msg_debug(modName, 
		    "Function fits_read_col_flt returned error code %d", status);
	return(VM_FALSE);
      }   

      if (fits_read_col_flt(sphotTable->fptr,deltaLColNo,i,1,1,null,
			    &(deltaLCol->colValue->fArray[i-1]),&null,&status)) {
	cpl_msg_debug(modName, 
		    "Function fits_read_col_flt returned error code %d", status);
	return(VM_FALSE);
      } 
  }  

  for (i = 0; i < 2; i++){     /* free the memory for the column labels */
    cpl_free( ttype[i] );
  }
  
  return (VM_TRUE);  
}



VimosBool writeFitsSphotTable(VimosTable *sphotTable, fitsfile *fptr)
{
  int nRows;
  int i;
  int status,nbytes;
  char *ttype[84], *tform[84];

  VimosColumn *lambdaCol;
  VimosColumn *magCol;
  VimosColumn *deltaLCol;

  char         modName[] = "writeFitsSphotTable";
  
  /* validate input */
  if (sphotTable == NULL) {
    cpl_msg_debug(modName, "NULL input table");
    return (VM_FALSE);
  }
  if (strcmp(sphotTable->name, VM_SPHOT) ) {
    cpl_msg_debug(modName, "Invalid input table");
    return (VM_FALSE);
  }

  nRows = sphotTable->cols->len;
  sphotTable->fptr = fptr;
  status = 0;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "SPH", 0, &status) ) {
    if(fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_debug(modName, 
      "Function fits_delete_hdu returned error code %d", status);
      return (VM_FALSE);
    }
  } else {
    status = 0;
  }

  for (i = 0; i <= 2; i++){      /* allocate space for the column labels */
    ttype[i] = (char *) cpl_malloc(FLEN_VALUE);  /* max label length = 69 */
    /* check if space was allocated */
    if (ttype[i] == NULL) {
      cpl_msg_debug(modName, "Allocation error");
      return(VM_FALSE);
    }
    tform[i] = (char *) cpl_malloc(FLEN_VALUE);
    /* check if space was allocated */
    if (tform[i] == NULL) {
      cpl_msg_debug(modName, "Allocation error");
      return(VM_FALSE);
    }
  }

  ttype[0] = "LAMBDA";
  tform[0] = "1E";
  ttype[1]="MAG";
  tform[1] = "1E";
  ttype[2]="DELTA_LAMBDA"; 
  tform[2] = "1E";

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl(fptr, BINARY_TBL,0, 3, ttype, tform,NULL, 
		      "SPH", &status)) {
    cpl_msg_debug(modName, 
    "Function fits_create_tbl returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "SPH", 0, &status)) {
    cpl_msg_debug(modName, 
    "Function fits_movnam_hdu returned error code %d", status);
    return(VM_FALSE);
  }

  /* write the table descriptors to the Fits file */
  
  if (fits_read_key (sphotTable->fptr,TINT,"NAXIS1",&nbytes,NULL,&status)) {
    cpl_msg_debug(modName, 
    "Function fits_read_key returned error code %d", status);
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(sphotTable->descs), "NAXIS1", nbytes, "")) {
    cpl_msg_debug(modName, "Function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(sphotTable->descs), "NAXIS2", nRows, "")) {
    cpl_msg_debug(modName, "Function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(sphotTable->descs),"TFIELDS",3, "")) {
    cpl_msg_debug(modName, "Function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  } 
  if (!writeDescsToFitsTable(sphotTable->descs, sphotTable->fptr)) {
    cpl_msg_debug(modName, "Function writeDescsToFitsTable returned an error");
    return(VM_FALSE);
  }

  lambdaCol = sphotTable->cols;
  magCol = lambdaCol->next;
  deltaLCol = magCol->next;
  
  for (i = 1; i <= nRows; i++) {
    if (fits_write_col_flt(sphotTable->fptr,1,i,1,1,
		&(lambdaCol->colValue->fArray[i-1]),&status)) {
      cpl_msg_debug(modName, 
      "Function fits_write_col_flt returned error code %d", status);
      return(VM_FALSE);
    }
    if (fits_write_col_flt(sphotTable->fptr,2,i,1,1,
		&(magCol->colValue->fArray[i-1]),&status)) {
      cpl_msg_debug(modName, 
      "Function fits_write_col_flt returned error code %d", status);
      return(VM_FALSE);
    }
    if (fits_write_col_flt(sphotTable->fptr,3,i,1,1,
		&(deltaLCol->colValue->fArray[i-1]),&status)) {
      cpl_msg_debug(modName, 
      "Function fits_write_col_flt returned error code %d", status);
      return(VM_FALSE);
    }
  }  
  return(VM_TRUE);
}
