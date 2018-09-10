/* $Id: vmlinecatalog.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include "vmlinecatalog.h"
#include "cpl.h"


/*
  Allocate a new line Catalog. This is a VimosTable with name VM_LIN (which
  should have value "IN").
  No descriptors or columns are allocated (but this may change in a later
  version).
*/
VimosTable *newLineCatalog()
{
  VimosTable *newTab;
  
  /* allocate new VimosTable */
  newTab = newTable();

  /* check if new VimosTable was allocated */
  if (newTab == NULL) {
    cpl_msg_error("newLineCatalog","The function newTable has returned NULL");
    return(NULL);
  }
  
  /* copy "LIN" into name of table */
  strcpy(newTab->name, VM_LIN);
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_LIN, "");
  if (newTab->descs == NULL) {
    /* cleanup */
    cpl_free(newTab);
    cpl_msg_error("newLineCatalog","The function newStringDescriptor has returned NULL");
    return(NULL);
  }

  /* should create all the descriptors and columns here ? */
  
  /* should create all the descriptors and columns here ?*/
  
  /* return address of new Line Catalog */
  return(newTab);
  
}

/*
  Delete a Line Catalog. This is just an esthetic wrapper for 
  deleteTable(lCat)
*/
void deleteLineCatalog(VimosTable *lCat)
{
  /* just a wrapper for deleteTabel() */
  deleteTable(lCat);
}

VimosBool readFitsLineCatalog(VimosTable *linCat, fitsfile *fptr)
{
  int    i;
  int    nCols;
  int    nRows;
  int    null;
  int    wlenColNo, nameColNo, fluxColNo, commColNo;
  int    status, nfound, ii;
  long   lenName, lenComm, width;
  int    typecode;
  char **ttype;
  char   comment[80];
  char   nstr[] = {'\0'};

  VimosColumn *wlenCol;
  VimosColumn *nameCol;
  VimosColumn *fluxCol;
  VimosColumn *commCol;

  char         modName[] = "readFitsLineCatalog";

  status = 0;

  /* validate input */
  if (linCat == NULL) {
    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  /* validate input */
  if ( strcmp(linCat->name, VM_LIN) ) {
    cpl_msg_debug(modName, "Invalid input table");
    return(VM_FALSE);
  }
  
  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, "LIN", 0, &status)) {
    cpl_msg_debug(modName, 
    "The function fits_movnam_hdu returned error code %d", status);
    return(VM_FALSE);
  }
  
  linCat->fptr = fptr;
  
  /* read Table */
  if (!readDescsFromFitsTable(&(linCat->descs), linCat->fptr)) {
    cpl_msg_debug(modName, "Function readDescsFromFitsTable returned an error");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(linCat->descs, "TFIELDS", &nCols, comment)) {
    cpl_msg_debug(modName, "The function readIntDescriptor returned an error");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(linCat->descs, "NAXIS2", &nRows, comment)) {
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
  if (fits_read_keys_str(linCat->fptr, "TTYPE", 1, nCols, ttype, &nfound, 
			 &status)) {
    cpl_msg_debug(modName, 
    "Function fits_read_keys_str returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(linCat->fptr,CASEINSEN,"WLEN",&wlenColNo,&status)) {
    cpl_msg_debug(modName, 
    "The function fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(linCat->fptr,CASEINSEN,"NAME",&nameColNo,&status)) {
    cpl_msg_debug(modName, 
    "The function fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_coltype(linCat->fptr, nameColNo, 
                       &typecode, &lenName, &width, &status)) {
    cpl_msg_debug(modName,
    "The function fits_get_coltype returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(linCat->fptr,CASEINSEN,"FLUX",&fluxColNo,&status)) {
    cpl_msg_debug(modName, 
    "The function fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(linCat->fptr,CASEINSEN,"COMMENT",&commColNo,&status)) {
    cpl_msg_debug(modName, 
    "Function fits_get_colnum returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_get_coltype(linCat->fptr, commColNo, 
                       &typecode, &lenComm, &width, &status)) {
    cpl_msg_debug(modName,
    "The function fits_get_coltype returned error code %d", status);
    return(VM_FALSE);
  }
  
  /* create column for wavelengths */
  wlenCol = newFloatColumn(nRows, "WLEN");
  if (wlenCol == NULL) {
    cpl_msg_debug(modName, "The function newFloatColumn has returned NULL");
    return(VM_FALSE);
  }

  /* put column in Table */
  linCat->cols = wlenCol;
  
  /* create column for NAME */
  nameCol = newStringColumn(nRows, "NAME");
  if (nameCol == NULL) {
    cpl_msg_debug(modName, "The function newStringColumn has returned NULL");
    return(VM_FALSE);
  }

  /* link column into table */
  wlenCol->next = nameCol;
  nameCol->prev = wlenCol;
  
  /* create column for FLUX */
  fluxCol = newFloatColumn(nRows, "FLUX");
  if (fluxCol == NULL) {
    cpl_msg_debug(modName, "The function newFloatColumn has returned NULL");
    return(VM_FALSE);
  }
  /* link column into table */
  nameCol->next = fluxCol;
  fluxCol->prev = nameCol;
  
  /* create column for COMMENT */
  commCol = newStringColumn(nRows, "COMMENT");
  if (commCol == NULL) {
    cpl_msg_debug(modName, "The function newStringColumn has returned NULL");
    return(VM_FALSE);
  }
  /* link column into table */
  fluxCol->next = commCol;
  commCol->prev = fluxCol;

  for (i = 1; i <= nRows; i++) {
      if (fits_read_col_flt(linCat->fptr,wlenColNo,i,1,1,null,
			    &(wlenCol->colValue->fArray[i-1]),&null,&status)){
	cpl_msg_debug(modName, 
        "Function fits_read_col_flt returned error code %d", status);
	return(VM_FALSE);
      } 
      
      nameCol->colValue->sArray[i-1] = 
                           (char *) cpl_malloc(((int)lenName + 1)*sizeof(char));

      if (fits_read_col_str(linCat->fptr,nameColNo,i,1,1,nstr,
			    &(nameCol->colValue->sArray[i-1]),&null,&status)){
	cpl_msg_debug(modName, 
        "Function fits_read_col_str returned error code %d", status);
	return(VM_FALSE);
      }  

    /* just to be sure ... */
    nameCol->colValue->sArray[i-1][lenName] = '\0';   
    if (fits_read_col_flt(linCat->fptr,fluxColNo,i,1,1,null,
			  &(fluxCol->colValue->fArray[i-1]),&null,&status)) {
      cpl_msg_debug(modName, 
      "Function fits_read_col_flt returned error code %d", status);
      return(VM_FALSE);
    }   

    commCol->colValue->sArray[i-1] = 
                         (char *) cpl_malloc(((int)lenComm + 1)*sizeof(char));

    if (fits_read_col_str(linCat->fptr,commColNo,i,1,1,nstr,
			  &(commCol->colValue->sArray[i-1]),&null,&status)) {
      cpl_msg_debug(modName, 
      "The function fits_read_col_str returned error code %d", status);
      return(VM_FALSE);
    }    
 
    commCol->colValue->sArray[i-1][lenComm] = '\0';
  }
   
  for (i = 0; i < 3; i++){     /* free the memory for the column labels */
    cpl_free( ttype[i] );
  }

  return (VM_TRUE);  
}



VimosBool writeFitsLineCatalog(VimosTable *linCat, fitsfile *fptr)
{
  int nRows;
  int i;
  int status,nbytes;
  char *ttype[84], *tform[84];

  VimosColumn *wlenCol;
  VimosColumn *nameCol;
  VimosColumn *fluxCol;
  VimosColumn *commCol;

  char         modName[] = "writeFitsLineCatalog";
  
  /* validate input */
  if (linCat == NULL) {
    cpl_msg_debug(modName, "NULL input table");
    return (VM_FALSE);
  }
  if ( strcmp(linCat->name, VM_LIN) ) {
    cpl_msg_debug(modName, "Invalid input table");
    return (VM_FALSE);
  }

  nRows = linCat->cols->len;
  linCat->fptr = fptr;
  status = 0;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "LIN", 0, &status) ) {
    if(fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_debug(modName, 
      "Function fits_delete_hdu returned error code %d", status);
      return (VM_FALSE);
    }
  } else {
    status = 0;
  }

  for (i = 0; i <= 3; i++){      /* allocate space for the column labels */
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

  ttype[0] = "WLEN";
  tform[0] = "1E";
  ttype[1]="NAME";
  tform[1] = "9A";
  ttype[2]="FLUX"; 
  tform[2] = "1E";
  ttype[3] = "COMMENT";
  tform[3] = "20A";

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl(fptr, BINARY_TBL,0, 4, ttype, tform,NULL, 
		      "LIN", &status)) {
    cpl_msg_debug(modName, 
    "Function fits_create_tbl returned error code %d", status);
    return(VM_FALSE);
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "LIN", 0, &status)) {
    cpl_msg_debug(modName, 
    "Function fits_movnam_hdu returned error code %d", status);
    return(VM_FALSE);
  }

  /* write the table descriptors to the Fits file */
  
  if (fits_read_key (linCat->fptr,TINT,"NAXIS1",&nbytes,NULL,&status)) {
    cpl_msg_debug(modName, 
    "Function fits_read_key returned error code %d", status);
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(linCat->descs), "NAXIS1", nbytes, "")) {
    cpl_msg_debug(modName, "Function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(linCat->descs), "NAXIS2", nRows, "")) {
    cpl_msg_debug(modName, "Function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(linCat->descs),"TFIELDS",4, "")) {
    cpl_msg_debug(modName, "Function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  } 
  if (!writeDescsToFitsTable(linCat->descs, linCat->fptr)) {
    cpl_msg_debug(modName, "Function writeDescsToFitsTable returned an error");
    return(VM_FALSE);
  }

  wlenCol = linCat->cols;
  nameCol = wlenCol->next;
  fluxCol = nameCol->next;
  commCol = fluxCol->next;
  
  for (i = 1; i <= nRows; i++) {
    if (fits_write_col_flt(linCat->fptr,1,i,1,1,
		&(wlenCol->colValue->fArray[i-1]),&status)) {
      cpl_msg_debug(modName, 
      "Function fits_write_col_flt returned error code %d", status);
      return(VM_FALSE);
    }
    if (fits_write_col_str(linCat->fptr,2,i,1,1,
		&(nameCol->colValue->sArray[i-1]),&status)) {
      cpl_msg_debug(modName, 
      "Function fits_write_col_str returned error code %d", status);
      return(VM_FALSE);
    }
    if (fits_write_col_flt(linCat->fptr,3,i,1,1,
		&(fluxCol->colValue->fArray[i-1]),&status)) {
      cpl_msg_debug(modName, 
      "Function fits_write_col_flt returned error code %d", status);
      return(VM_FALSE);
    }
    if(fits_write_col_str(linCat->fptr,4,i,1,1,
		&(commCol->colValue->sArray[i-1]),&status)) {
      cpl_msg_debug(modName, 
      "Function fits_write_col_str returned error code %d", status);
      return(VM_FALSE);
    }
  }  
  return(VM_TRUE);
}
