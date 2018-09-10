/* $Id: vmobjecttable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <string.h>

#include <fitsio.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilerrno.h>

#include "vmtable.h"
#include "vmwindowtable.h"
#include "vmobjecttable.h"
#include "cpl.h"


/*
 * Returns a pointer to a new Object Object
 */

VimosObjectObject *newObjectObject(void)
{
  VimosObjectObject *newObject;
  
  /* allocate memory for Object Object */
  newObject = (VimosObjectObject *) cpl_malloc(sizeof(VimosObjectObject));
  
  /* check if space was allocated */
  if (newObject == NULL) {
    cpl_msg_error("newObjectObject","Allocation Error");
    return(NULL);
  }
  
  /* fill up fields with default values */
  newObject->slitNo   = 0;
  newObject->IFUslitNo = 0;
  newObject->IFUfibNo = 0;
  newObject->rowNum   = 0;
  newObject->objNo    = 0;
  newObject->objX     = 0.0;
  newObject->objY     = 0.0;
  newObject->objRA    = 0.0;
  newObject->objDec   = 0.0;

  newObject->prev = newObject->next = NULL;
  
  /* return to caller */
  return(newObject);
  
}

/*
Deletes all Object Objects contained in the list
*/
void deleteObjectObject(VimosObjectObject *obj)
{
  VimosObjectObject    *tmpObj;
  VimosObjectObject    *nextObj;
  
  /* store start of the list */
  tmpObj = obj;
  
  /* traverse  list */
  while (tmpObj) {
    /* get address of next object */
    nextObj = tmpObj->next;
    /* free current object */
    cpl_free(tmpObj);
    /* next one to process */
    tmpObj = nextObj;
  }
}


  
  
/*
  Allocate a new Object Table. 
*/
VimosObjectTable *newObjectTable()
{
  VimosObjectTable *newTab;
  
  /* allocate new VimosTable */
  newTab = (VimosObjectTable *) cpl_malloc(sizeof(VimosObjectTable));

  /* check if space was allocated */
  if (newTab == NULL) {
    cpl_msg_error("newObjectTable","Allocation Error");
    return(NULL);
  }
  
  /* copy "OBJ" into name of table and initialize fields*/
  strcpy(newTab->name, VM_OBJ);
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_OBJ, "");
  if (newTab->descs == NULL) {
    /* cleanup */
    cpl_free(newTab);
    cpl_msg_error("newObjectTable","The function newStringDescriptor has returned NULL");
    return(NULL);
  }

  newTab->objs = NULL;
  newTab->fptr = NULL;
  
  /* return address of new Object Table */
  return(newTab);
  
}

/*

*/
void deleteObjectTable(VimosObjectTable *oTable)
{
  VimosDescriptor *tmpDesc;
  VimosDescriptor *nxtDesc;
  
  if (oTable == NULL) {
    return;
  }
  
  deleteObjectObject(oTable->objs);

  tmpDesc = oTable->descs;  
  while (tmpDesc) {
    nxtDesc = tmpDesc->next;
    deleteDescriptor(tmpDesc);
    tmpDesc = nxtDesc;
  }

}





/* 
   Copy ALL descriptors (including those that may not be relevant) of a
   Window Table into Object Table. The input tables are checked that they
   are a Window Table and a Object Table.
*/
VimosBool copyWinTab2ObjTab(VimosWindowTable *winTable, 
                            VimosObjectTable *objTable)
{
  /* is this a WindowTable ?*/
  if (strcmp(winTable->name, VM_WIN)) {
    cpl_msg_error("copyWinTab2ObjTab","There is no Window Table");
    return(VM_FALSE);
  }
  /* is this a Object Table ?*/
  if (strcmp(objTable->name, VM_OBJ)) {
    cpl_msg_error("copyWinTab2ObjTab","There is no Object Table");
    return(VM_FALSE);
  }

  /* copy all descriptors from Window Table to Object Table*/
  if (!copyAllDescriptors(winTable->descs, &(objTable->descs))) {
    cpl_msg_error("copyWinTab2ObjTab","The function copyAllDescriptors has returned an error");
    return(VM_FALSE);
  }
  /* overwrite Descriptor TABLE with "OBJ"*/
  if (!writeStringDescriptor(&(objTable->descs), "ESO PRO TABLE", VM_OBJ, "")){
    cpl_msg_error("copyWinTab2ObjTab","The function writeStringDescriptor has returned an error");
    return(VM_FALSE);
  }
  if (!writeStringDescriptor(&(objTable->descs), "EXTNAME", "OBJ", "")) {
    cpl_msg_error("copyWinTab2ObjTab","The function writeStringDescriptor has returned an error");
    return(VM_FALSE);
  }

  return(VM_TRUE);
}




/* Return number of objects in a VimosObjectObject. Since a VimosObjectObject
  is a linked list of ObjectObjects, more than one objects can be contained in
  the list.  
*/
int numObjectsInObjectObject(VimosObjectObject *oObject)
{
  int numObjects;
  
  /* redundant, but nevertheless...*/
  if (oObject == NULL) {
    cpl_msg_error("numObjectsInObjectObject","NULL input pointer");
    pilErrno = 1;
    return 0;
  }
  
  /* initialize */
  numObjects = 0;
  
  /* traverse list and count */
  while (oObject) {
    numObjects++;
    oObject = oObject->next;
  }
  
  /* return number of objects */
  return(numObjects);
}


VimosBool readFitsObjectTable(VimosObjectTable *objTable, fitsfile *fptr)
{
  int    i;
  int    nCols;
  int    nRows;
  int    dummyInt;
  int    null;
  int    slitCol, IFUslitNoCol, IFUfibNoCol, yCol;
  int    objNoCol, objXCol, objYCol, objRACol, objDecCol;

  int    status = 0;
  int    ii;
  char **ttype;
  char    comment[80];

  VimosObjectObject *object;
  VimosObjectObject *lastObject;
  
  /* validate input */
  if (objTable == NULL) {
    cpl_msg_error("readFitsObjectTable","NULL input table");
    return(VM_FALSE);
  }
  /* validate input */
  if ( strcmp(objTable->name, VM_OBJ) ) {
    cpl_msg_error("readFitsObjectTable","Invalid input table");
    return(VM_FALSE);
  }
  
  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, "OBJ", 0, &status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  
  objTable->fptr = fptr;
  
  /* read Table */
  if (!readDescsFromFitsTable(&(objTable->descs), objTable->fptr)) {
    cpl_msg_error("readFitsObjectTable","The function readDescsFromFitsTable has returned an error");
    return(VM_FALSE);
  }
  if (!readIntDescriptor(objTable->descs, "TFIELDS", &nCols, comment)) {
    cpl_msg_error("readFitsObjectTable","The function readIntDescriptor has returned an error");
    return(VM_FALSE);
  }
  if (!readIntDescriptor(objTable->descs, "NAXIS2", &nRows, comment)) {
    cpl_msg_error("readFitsObjectTable","The function readIntDescriptor has returned an error");
    return(VM_FALSE);
  }

 /* allocate space for the column labels */
  ttype = (char **) cpl_malloc(nCols*sizeof(char *));
  for (ii = 0; ii < nCols; ii++) {
    ttype[ii] = (char *) cpl_malloc(FLEN_VALUE*sizeof(char));
    /* check if space was allocated */
    if (ttype[ii] == NULL) {
      cpl_msg_error("readFitsObjectTable","Allocation Error");
      return(VM_FALSE);
    }
  }
  /* read the column names from the TTYPEn keywords */
  if (fits_read_keys_str(objTable->fptr, "TTYPE", 1, nCols, ttype, &dummyInt, 
			 &status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_read_keys_str has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"SLIT",&slitCol,&status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"IFUSLIT_NO",&IFUslitNoCol,
		      &status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"IFUFIBER_NO",&IFUfibNoCol,
		      &status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"Y",&yCol,&status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"OBJ_NO",&objNoCol,&status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"OBJ_X",&objXCol,&status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"OBJ_Y",&objYCol,&status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"OBJ_RA",&objRACol,&status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_get_colnum(objTable->fptr,CASEINSEN,"OBJ_DEC",&objDecCol,&status)) {
    cpl_msg_error("readFitsObjectTable","The function fits_get_colnum has returned an error (code %d)", status);
    return(VM_FALSE);
  }

  lastObject = NULL;

  i = 1;
  while (i <= nRows) {
    object = newObjectObject();
    if (object == NULL) {
      cpl_msg_error("readFitsObjectTable","The function newObjectObject has returned NULL");
      return(VM_FALSE);
    }
    if (fits_read_col_int(objTable->fptr,slitCol,i,1,1,null,&(object->slitNo),
			  &null,&status)) {
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    } 
    if (fits_read_col_int(objTable->fptr,IFUslitNoCol,i,1,1,null,
			  &(object->IFUslitNo),&null,&status)) {
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    }  
    if (fits_read_col_int(objTable->fptr,IFUfibNoCol,i,1,1,null,
			  &(object->IFUfibNo),&null,&status)) {
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    } 
    if (fits_read_col_int(objTable->fptr,yCol,i,1,1,null,&(object->rowNum),
			  &null,&status)) {
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    }  
    if (fits_read_col_int(objTable->fptr,objNoCol,i,1,1,null,&(object->objNo),
			  &null,&status)) {
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    }  
    if (fits_read_col_flt(objTable->fptr,objXCol,i,1,1,null,&(object->objX),
			  &null,&status)) {
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_flt has returned an error (code %d)", status);
      return(VM_FALSE);
    }  
    if (fits_read_col_flt(objTable->fptr,objYCol,i,1,1,null,&(object->objY),
			  &null,&status)) {
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_flt has returned an error (code %d)", status);
      return(VM_FALSE);
    }  
    if (fits_read_col_dbl(objTable->fptr,objRACol,i,1,1,null,&(object->objRA),
			  &null,&status)) {
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_dbl has returned an error (code %d)", status);
      return(VM_FALSE);
    }   
    if (fits_read_col_dbl(objTable->fptr,objDecCol,i,1,1,null,
			  &(object->objDec),&null,&status)) { 
      cpl_msg_error("readFitsObjectTable","The function fits_read_col_dbl has returned an error (code %d)", status);
      return(VM_FALSE);
    }  
    i++;
    
    if (lastObject == NULL) {
      objTable->objs = object;
    }
    else {
      lastObject->next = object;
      object->prev = lastObject;
    }
    lastObject = object;
        
  }
   
  for (ii = 0; ii < 3; ii++){     /* free the memory for the column labels */
    cpl_free( ttype[ii] );
  }

  return (VM_TRUE);
}



VimosBool writeFitsObjectTable(VimosObjectTable *objTable, fitsfile *fptr)
{
  int nRows;
  int rowNum;
  int i, status, nbytes;
  char *ttype[84], *tform[84];

  VimosObjectObject *object;
  
  
  /* validate input */
  if (objTable == NULL) {
    cpl_msg_error("writeFitsObjectTable","NULL input table");
    return(VM_FALSE);
  }
  /* validate input */
  if ( strcmp(objTable->name, VM_OBJ) ) {
    cpl_msg_error("writeFitsObjectTable","Invalid input table");
    return(VM_FALSE);
  }
  
  objTable->fptr = fptr;

  nRows = numObjectsInObjectObject(objTable->objs);
  if (nRows == 0) {
    cpl_msg_error("writeFitsObjectTable","The function numObjectsInObjectObject has returned an error");
    return(VM_FALSE);
  }
 
  status = 0;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "OBJ", 0, &status) ) {
    if (fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_delete_hdu has returned an error (code %d)", status);
      return(VM_FALSE);
    }
  } else {
    status = 0;
  }
 
  for (i = 0; i <= 8; i++){      /* allocate space for the column labels */
    ttype[i] = (char *) cpl_malloc(FLEN_VALUE);  /* max label length = 69 */
    /* check if space was allocated */
    if (ttype[i] == NULL) {
      cpl_msg_error("writeFitsObjectTable","Allocation Error");
      return(VM_FALSE);
    }
    tform[i] = (char *) cpl_malloc(FLEN_VALUE);
    /* check if space was allocated */
    if (tform[i] == NULL) {
      cpl_msg_error("writeFitsObjectTable","Allocation Error");
      return(VM_FALSE);
    }
  }

  ttype[0] = "SLIT";
  tform[0] = "1J";
  ttype[1] = "IFUSLIT_NO";
  tform[1] = "1J";
  ttype[2] = "IFUFIBER_NO";
  tform[2] = "1J";
  ttype[3] = "Y";
  tform[3] = "1J";
  ttype[4] = "OBJ_NO"; 
  tform[4] = "1J";
  ttype[5] = "OBJ_X";
  tform[5] = "1E";
  ttype[6] = "OBJ_Y";
  tform[6] = "1E";
  ttype[7] = "OBJ_RA";
  tform[7] = "1D";
  ttype[8] = "OBJ_DEC";
  tform[8] = "1D";

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl(fptr, BINARY_TBL,0, 9, ttype, tform,NULL, 
		      "OBJ", &status)) {
    cpl_msg_error("writeFitsObjectTable","The function fits_create_tbl has returned an error (code %d)", status);
    return(VM_FALSE);
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "OBJ", 0, &status)) {
    cpl_msg_error("writeFitsObjectTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return(VM_FALSE);
  }

  /* write the table descriptors to the Fits file */
  
  if (fits_read_key (objTable->fptr,TINT,"NAXIS1",&nbytes,NULL,&status)) {
    cpl_msg_error("writeFitsObjectTable","The function fits_read_key has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(objTable->descs), "NAXIS1", nbytes, "")) {
    cpl_msg_error("writeFitsObjectTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(objTable->descs), "NAXIS2", nRows, "")) {
    cpl_msg_error("writeFitsObjectTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(objTable->descs),"TFIELDS",9, "")) {
    cpl_msg_error("writeFitsObjectTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeDescsToFitsTable(objTable->descs, objTable->fptr)) {
    cpl_msg_error("writeFitsObjectTable","The function writeDescsToFitsTable has returned an error");
    return(VM_FALSE);
  }


  object = objTable->objs;
  rowNum = 1;
  while (object) {
    if (fits_write_col_int(objTable->fptr,1,rowNum,1,1,&(object->slitNo),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    }
    if (fits_write_col_int(objTable->fptr,2,rowNum,1,1,&(object->IFUslitNo),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    }
    if (fits_write_col_int(objTable->fptr,3,rowNum,1,1,&(object->IFUfibNo),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    }
    if (fits_write_col_int(objTable->fptr,4,rowNum,1,1,&(object->rowNum),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    }
    if (fits_write_col_int(objTable->fptr,5,rowNum,1,1,&(object->objNo),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_int has returned an error (code %d)", status);
      return(VM_FALSE);
    }
    if (fits_write_col_flt(objTable->fptr,6,rowNum,1,1,&(object->objX),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_flt has returned an error (code %d)", status);
      return(VM_FALSE);
    }
    if (fits_write_col_flt(objTable->fptr,7,rowNum,1,1,&(object->objY),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_flt has returned an error (code %d)", status);
      return(VM_FALSE);  
    } 
    if (fits_write_col_dbl(objTable->fptr,8,rowNum,1,1,&(object->objRA),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_dbl has returned an error (code %d)", status);
      return(VM_FALSE);   
    }
    if (fits_write_col_dbl(objTable->fptr,9,rowNum,1,1,&(object->objDec),
			   &status)) {
      cpl_msg_error("writeFitsObjectTable","The function fits_write_col_dbl has returned an error (code %d)", status);
      return(VM_FALSE);
    }

    rowNum++;
    object = object->next;
  }
  
  return(VM_TRUE);  
}
