/* $Id: vmwindowtable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include "cpl.h"

/*
 * Returns a pointer to a new Window Object
 */

VimosWindowObject *newWindowObject(void)
{
  VimosWindowObject *newObject;
  
  /* allocate memory for Window Object */
  newObject = (VimosWindowObject *) cpl_malloc(sizeof(VimosWindowObject));
  
  /* check if space was allocated */
  if (newObject == NULL) {
    cpl_msg_error("newWindowObject","Allocation Error");
    return(NULL);
  }
  
  /* fill up fields with default values */
  newObject->objStart = 0;
  newObject->objEnd   = 0;
  newObject->objNo    = 0;
  newObject->objPos   = 0.0;
  newObject->objWidth = 0.0;
  newObject->objX     = 0.0;
  newObject->objY     = 0.0;
  newObject->objProfile = NULL;
  newObject->posDef   = VM_FALSE;
  newObject->objRA    = 0.0;
  newObject->objDec   = 0.0;
  newObject->parDef   = VM_FALSE;
/*    newObject->objPeak  = 0.0;
      newObject->skyLevel = 0.0; */

  newObject->prev = newObject->next = NULL;
  
  /* return to caller */
  return(newObject);
  
}

/*
Deletes all Window Objects contained in the list
*/
void deleteWindowObject(VimosWindowObject *obj)
{
  VimosWindowObject    *tmpObj;
  VimosWindowObject    *nextObj;
  
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
  Returns a pointer to a new Window Slit
*/

VimosWindowSlit *newWindowSlit()
{
  VimosWindowSlit *newSlit;
  
  /* allocate memory for WindowSlit */
  newSlit = (VimosWindowSlit *) cpl_malloc(sizeof(VimosWindowSlit));
  
  /* check if space was allocated */
  if (newSlit == NULL) {
    cpl_msg_error("newWindowSlit","Allocation Error");
    return(NULL);
  }
  
  /* fill up fields with default values */
  newSlit->slitNo = 0;
  newSlit->IFUslitNo = 0;
  newSlit->IFUfibNo = 0;
  newSlit->IFUfibTrans = 0.;
  newSlit->specLong = VM_FALSE;
  newSlit->specStart = 0;
  newSlit->specEnd = 0;
  newSlit->prev = newSlit->next = NULL;
  newSlit->numObj = 0;
  newSlit->objs = NULL;
  
  /* return to caller */
  return(newSlit);
  
}

/*
Deletes all Window Slits contained in the list
*/
void deleteWindowSlit(VimosWindowSlit *wSlit)
{
  VimosWindowSlit    *tmpSlit;
  VimosWindowSlit    *nextSlit;
  
  /* store start of the list */
  tmpSlit = wSlit;
  
  /* traverse  list */
  while (tmpSlit) {
    /* delete the objects in this slit */
    deleteWindowObject(tmpSlit->objs);
    /* get address of next slit */
    nextSlit = tmpSlit->next;
    /* free currebt slit */
    cpl_free(tmpSlit);
    /* next one to process */
    tmpSlit = nextSlit;
  }
}

  
  
/*
  Allocate a new Window Table. 
*/
VimosWindowTable *newWindowTable()
{
  VimosWindowTable *newTab;
  
  /* allocate new VimosTable */
  newTab = (VimosWindowTable *) cpl_malloc(sizeof(VimosWindowTable));

  /* check if space was allocated */
  if (newTab == NULL) {
    cpl_msg_error("newWindowTable","Allocation Error");
    return(NULL);
  }
  
  /* copy "WIN" into name of table and initialize fields*/
  strcpy(newTab->name, VM_WIN);
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_WIN, "");
  if (newTab->descs == NULL) {
    /* cleanup */
    cpl_free(newTab);
    cpl_msg_error("newWindowTable","The function newStringDescriptor has returned NULL");
    return(NULL);
  }

  newTab->slits = NULL;
  newTab->fptr = NULL;
  
  /* return address of new Window Table */
  return(newTab);
  
}

/*

*/
void deleteWindowTable(VimosWindowTable *wTable)
{
  VimosDescriptor *tmpDesc;
  VimosDescriptor *nxtDesc;
  
  if (wTable == NULL) {
    return;
  }
  
  deleteWindowSlit(wTable->slits);

  tmpDesc = wTable->descs;  
  while (tmpDesc) {
    nxtDesc = tmpDesc->next;
    deleteDescriptor(tmpDesc);
    tmpDesc = nxtDesc;
  }

}





/* 
   Copy ALL descriptors (including those that may not be relevant) of an
   Extraction Table into Window Table. The input tables are checked that they
   are an Extraction Table and a Window Table.
*/
VimosBool copyExtTab2WinTab(VimosExtractionTable *extTable, 
                            VimosWindowTable *winTable)
{
  /* is this a Window Table ?*/
  if (strcmp(winTable->name, VM_WIN)) {
    cpl_msg_error("copyExtTab2WinTab","There is no Window Table");
    return(VM_FALSE);
  }
  /* is this a Extraction  Table ?*/
  if (strcmp(extTable->name, VM_EXT)) {
    cpl_msg_error("copyExtTab2WinTab","There is no Extraction Table");
    return(VM_FALSE);
  }
    
  /* some descriptors should not be copied!! e.g. TABLE */  
  if (!copyAllDescriptors(extTable->descs, &(winTable->descs))) {
    cpl_msg_error("copyExtTab2WinTab","The function copyAllDescriptors has returned an error");
    return(VM_FALSE);
  }
  if (!writeStringDescriptor(&(winTable->descs), "ESO PRO TABLE", VM_WIN, "")){
    cpl_msg_error("copyExtTab2WinTab","The function writeStringDescriptor has returned an error");
    return(VM_FALSE);
  }
  if (!writeStringDescriptor(&(winTable->descs), "EXTNAME", "WIN", "")) {
    cpl_msg_error("copyExtTab2WinTab","The function writeStringDescriptor has returned an error");
    return(VM_FALSE);
  }
  return(VM_TRUE);
}


/* 
   Return number of slits in a WindowTable 
*/
int numSlitsInWindowTable(VimosWindowTable *wTable)
{

  /* check input */
  if (wTable == NULL) {
    cpl_msg_error("numSlitsInWindowTable","NULL input table");
    pilErrno = 1;
    return(0);
  }
        
  /* is this a Window Table ?*/
  if (strcmp(wTable->name, VM_WIN)) {
    cpl_msg_error("numSlitsInWindowTable","The table in input is not a Window Table");
    pilErrno = 1;
    return(0);
  }
  
  /* delegate call */
  return(numSlitsInWindowSlit(wTable->slits));
      
}


/* Return number of objects in the Window Table */
int numObjsInWindowTable(VimosWindowTable *wTable)
{
  /* validate input */
  if (strcmp(wTable->name, VM_WIN)) {
    cpl_msg_error("numObjsInWindowTable","The table in input is not a Window Table");
    pilErrno = 1;
    return(0);
  }

  /* delegate call */
  return(numObjectsInWindowSlit(wTable->slits));
  
}



/*
  Return number of objects in a WindowObject. Since a WindowObject is a
  linked list of WindowObjects, more than one objects can be contained in the
  list.
*/
int numObjectsInWindowObject(VimosWindowObject *wObject)
{
  int numObjects;
  
  if (wObject == NULL) {
    cpl_msg_error("numObjectsInWindowObject","There is no Window Object");
    pilErrno = 1;
    return(0);
  }

  /* initialize */
  numObjects = 0;
  
  /* traverse list and count */
  while (wObject) {
    numObjects++;
    wObject = wObject->next;
  }
  
  /* return number of objects */
  return(numObjects);
}



/*
  Return number of objects in a WindowSlit. Since a WindowSlit is a linked list
  of WindowSlits, more than one slit can be contained in the list. The sum 
  of all objects in all slits is returned.
  
*/
int numObjectsInWindowSlit(VimosWindowSlit *wSlit)
{
  int numObjects;
  
  if (wSlit == NULL) {
    cpl_msg_error("numObjectsInWindowSlit","There is no Window Slit");
    pilErrno = 1;
    return(0);
  }

  /* initialize */
  numObjects = 0;
  
  /* traverse list and count */
  while (wSlit) {
    numObjects += wSlit->numObj;
    wSlit = wSlit->next;
  }
  
  return(numObjects);
  
}

/*
  Return number of slits where no object was detected in a WindowSlit. 
  Since a WindowSlit is a linked list of WindowSlits, more than one slit 
  can be contained in the list. The sum of such "empty"  slits 
  is returned.
  
*/
int numEmptySlitsInWindowSlit(VimosWindowSlit *wSlit)
{
  int numSlits;
  
  if (wSlit == NULL) {
    cpl_msg_error("numObjectsInWindowSlit","There is no Window Slit");
    pilErrno = 1;
    return(0);
  }

  /* initialize */
  numSlits = 0;
  
  /* traverse list and count */
  while (wSlit) {
    if (wSlit->numObj == 0)
      numSlits += 1;
    wSlit = wSlit->next;
  }
  
  return(numSlits);
  
}


/* Return number of slits in a Window Slit. Since a Window Slit is a linked
  list of Window Slits, more than one slit can be contained in the list.  */
int numSlitsInWindowSlit(VimosWindowSlit *wSlit)
{
  int numSlits;

  if (wSlit == NULL) {
    cpl_msg_error("numSlitsInWindowSlit","There is no Window Slit");
    pilErrno = 1;
    return(0);
  }

  /* initialize */  
  numSlits = 0;
  
  /* traverse list and count */
  while (wSlit) {
    numSlits++;
    wSlit = wSlit->next;
  }
  
  return(numSlits);
  
}



VimosBool windowObjectInRow(VimosWindowSlit *slit, int row)
{
  VimosWindowObject *obj;
  
  if (slit == NULL) {
    cpl_msg_error("windowObjectInRow","There is no Window Slit");
    return(VM_FALSE);
  }
  
  obj = slit->objs;
  
  while (obj) {
    if ( (row >= obj->objStart) && (row <= obj->objEnd)) {
      return(VM_TRUE);
    }
    obj = obj->next;
  }
  
  return(VM_FALSE);
    
}

/* add offSer to  objects in refWinTable and put them in winTable */
VimosBool shiftWindowObjects(VimosWindowTable *refWinTable, VimosWindowTable *winTable,  float offSet)
{

  int                   objNo            = 0;
  int                   wSlitWidth       = 0;

  float                 objPos           = 0.0;

  VimosWindowObject    *tmpObj;

  VimosWindowSlit      *wSlit,*refWSlit;




  wSlit = winTable->slits;
  refWSlit = refWinTable->slits;
  
  while(wSlit){
    wSlitWidth = wSlit->specEnd - wSlit->specStart; 
    wSlit->objs = NULL;
    
    /*  shift all objects */
    objNo=1;
    while(refWSlit->objs) {
      
      objPos = refWSlit->objs->objPos - offSet;
      
      if(objPos>0 && objPos < wSlitWidth) {
	tmpObj = newWindowObject();
	tmpObj->objNo = objNo;
	tmpObj->objPos = objPos;
	tmpObj->objStart = refWSlit->objs->objStart - offSet > 0 ?
	  refWSlit->objs->objStart - offSet : 0;
	tmpObj->objEnd = refWSlit->objs->objEnd - offSet < wSlitWidth ?
	  refWSlit->objs->objEnd - offSet : wSlitWidth;  
	tmpObj->objX = refWSlit->objs->objX - offSet;
	tmpObj->objY = refWSlit->objs->objY - offSet;
	if(!wSlit->objs) wSlit->objs=tmpObj;  
	else {
	  wSlit->objs->next = tmpObj;
	  tmpObj->prev = wSlit->objs;
	  wSlit->objs = tmpObj;
	}
	objNo++;
      }
      
      if(refWSlit->objs->next) refWSlit->objs = refWSlit->objs->next;
      else break;
    }
    
    /* rewind refWSlit->obj */
    if(refWSlit->objs)  while(refWSlit->objs->prev) refWSlit->objs = refWSlit->objs->prev;
    
    /* go to next slit */
    if(wSlit->next) {
      /* rewind current slit objects*/
      if(wSlit->objs) while(wSlit->objs->prev) wSlit->objs = wSlit->objs->prev;
      wSlit = wSlit->next;
      refWSlit=refWSlit->next;
    } else break;
  }
  
  /* rewind refWSlit */
  while(refWSlit->prev) refWSlit = refWSlit->prev;
  
  return(VM_TRUE);
  
  
}

VimosBool readFitsWindowTable(VimosWindowTable *winTable, fitsfile *fptr)
{
  int    i, ii;
  int    nCols;
  int    nRows;
  int    dummyInt;
  int    null;
  int    slitCol, IFUslitNoCol, IFUfibNoCol, specLenCol; 
  int    specStartCol, specEndCol;
  int    objNoCol, objStartCol, objEndCol;
  int    objPosCol, objXCol, objYCol, objRACol, objDecCol;
  int    objWidCol;
  /*  objPeakCol, skyLevCol; */

  int IFUfibTransCol;

  int    objNum;
  int    newSlitNo;
  int    specLong;

  int    status = 0;
  char **ttype;
  char   comment[80];

  VimosWindowSlit *winSlit;
  VimosWindowSlit *lastSlit;
  VimosWindowObject *winObj;
  VimosWindowObject *lastObject;
  
  /* validate input */
  if (winTable == NULL) {
    cpl_msg_error("readFitsWindowTable","NULL input table");
    return (VM_FALSE);
  }
  /* validate input */
  if ( strcmp(winTable->name, VM_WIN) ) {
    cpl_msg_error("readFitsWindowTable","Invalid input table");
    return (VM_FALSE);
  }
  
  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, "WIN", 0, &status) ){
    cpl_msg_error("readFitsWindowTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  
  winTable->fptr = fptr;
  
  /* read Table */
  if (!readDescsFromFitsTable(&(winTable->descs), winTable->fptr)) {
    cpl_msg_error("readFitsWindowTable","The function readDescsFromFitsTable has returned an error");
    return (VM_FALSE);
  }
  if (!readIntDescriptor(winTable->descs, "TFIELDS", &nCols, comment)) {
    cpl_msg_error("readFitsWindowTable","The function readIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!readIntDescriptor(winTable->descs, "NAXIS2", &nRows, comment)) {
    cpl_msg_error("readFitsWindowTable","The function readIntDescriptor has returned an error");
    return (VM_FALSE);
  }

  /* allocate space for the column labels */
  ttype = (char **) cpl_malloc(nCols*sizeof(char *));
  for (i=0; i < nCols; i++) {
    ttype[i] = (char *) cpl_malloc(FLEN_VALUE*sizeof(char));
    /* check if space was allocated */
    if (ttype[i] == NULL) {
      cpl_msg_error("readFitsWindowTable","Allocation Error");
      return (VM_FALSE);
    }
  }
      /* read the column names from the TTYPEn keywords */
  if (fits_read_keys_str(winTable->fptr, "TTYPE", 1, nCols, ttype, &dummyInt, 
			 &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_read_keys_str has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  if (fits_get_colnum(winTable->fptr,CASEINSEN,"SLIT",&slitCol,&status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"IFUSLIT_NO",&IFUslitNoCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"IFUFIB_NO",&IFUfibNoCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"SPEC_LEN", &specLenCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"SPEC_START",&specStartCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"SPEC_END", &specEndCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_START", &objStartCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_END", &objEndCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_NO", &objNoCol,&status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_POS", &objPosCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_WIDTH", &objWidCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }

/*  fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_PEAK", &objPeakCol,&status);
  fits_get_colnum(winTable->fptr,CASEINSEN,"SKY_LEV", &skyLevCol,&status); */


  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_X", &objXCol,&status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_Y", &objYCol,&status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_RA", &objRACol,&status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"OBJ_DEC", &objDecCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  /* ALEX: added IFUfibTrans */
  if (fits_get_colnum(winTable->fptr,CASEINSEN,"IFUFIB_TRANS",&IFUfibTransCol,
		      &status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }


  lastSlit = NULL;
  if (fits_read_col_int(winTable->fptr,slitCol,1,1,1,null,&newSlitNo,
			&null,&status)) {
    cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
    return (VM_FALSE);
  } 

  i = 1;
  while (i <= nRows) {
    winSlit = newWindowSlit();
    if (winSlit == NULL) {
      cpl_msg_error("readFitsWindowTable","The function newWindowSlit has returned NULL");
      return (VM_FALSE);
    }
    lastObject = NULL;
    objNum = 0;
    winSlit->slitNo = newSlitNo;

    while ( (i+objNum <= nRows) && (newSlitNo == winSlit->slitNo) ) {
      winObj = newWindowObject();
      if (winObj == NULL) {
	cpl_msg_error("readFitsWindowTable","The function newWindowObject has returned NULL");
        return (VM_FALSE);
      }
      if (fits_read_col_int(winTable->fptr,specLenCol,i+objNum,1,1,null,
			    &specLong,&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      } 

      if (specLong) {
        winSlit->specLong = VM_TRUE;
      }
      else {
        winSlit->specLong = VM_FALSE;
      }

      /* some redundant reading.. */   
      if (fits_read_col_int(winTable->fptr,IFUslitNoCol,i+objNum,1,1,null,
			    &(winSlit->IFUslitNo),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      }
      if (fits_read_col_int(winTable->fptr,IFUfibNoCol,i+objNum,1,1,null,
			    &(winSlit->IFUfibNo),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      }

      /* ALEX: added IFUfibTrans */
      if (fits_read_col_flt(winTable->fptr,IFUfibTransCol,i+objNum,1,1,null,
			    &(winSlit->IFUfibTrans),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_flt has returned an error (code %d)", status);
        return (VM_FALSE);
      }


      if (fits_read_col_int(winTable->fptr,specStartCol,i+objNum,1,1,null,
			    &(winSlit->specStart),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      } 
      if (fits_read_col_int(winTable->fptr,specEndCol,i+objNum,1,1,null,
			    &(winSlit->specEnd),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      }
      if (fits_read_col_int(winTable->fptr,objStartCol,i+objNum,1,1,null,
			    &(winObj->objStart),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      } 
      if (fits_read_col_int(winTable->fptr,objEndCol,i+objNum,1,1,null,
			    &(winObj->objEnd),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      } 
      if (fits_read_col_int(winTable->fptr,objNoCol,i+objNum,1,1,null,
			    &(winObj->objNo),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      } 
      if (fits_read_col_flt(winTable->fptr,objPosCol,i+objNum,1,1,null,
			    &(winObj->objPos),&null,&status)) {
	cpl_msg_error("readFitsWindowTable","The function fits_read_col_flt has returned an error (code %d)", status);
        return (VM_FALSE);
      }
      if (fits_read_col_flt(winTable->fptr,objWidCol,i+objNum,1,1,null,
			    &(winObj->objWidth),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_dbl has returned an error (code %d)", status);
        return (VM_FALSE);
      } 
 

/*       if (null) {
         winObj->parDef = VM_FALSE;
         winObj->objWidth = 0.0;
         winObj->objPeak = 0.0;
         winObj->skyLevel = 0.0;
       } else {
         winObj->parDef = VM_TRUE;
	 fits_read_col_flt(winTable->fptr,objPeakCol,i+objNum,1,1,null,
		 	  &(winObj->objPeak),&null,&status); 
	fits_read _col_flt(winTable->fptr,skyLevCol,i+objNum,1,1,null,
			   &(winObj->skyLevel),&null,&status); 
       } */


      if (fits_read_col_flt(winTable->fptr,objXCol,i+objNum,1,1,null,
			    &(winObj->objX),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_dbl has returned an error (code %d)", status);
        return (VM_FALSE);
      } 
      if (fits_read_col_flt(winTable->fptr,objYCol,i+objNum,1,1,null,
			    &(winObj->objY),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_dbl has returned an error (code %d)", status);
        return (VM_FALSE);
      } 
      if (fits_read_col_dbl(winTable->fptr,objRACol,i+objNum,1,1,null,
			    &(winObj->objRA),&null,&status)) {
        cpl_msg_error("readFitsWindowTable","The function fits_read_col_dbl has returned an error (code %d)", status);
        return (VM_FALSE);
      } 

      if (null) {
        winObj->posDef = VM_FALSE;
        winObj->objRA = 0.0;
        winObj->objDec = 0.0;
      } else {
        winObj->posDef = VM_TRUE;
	if (fits_read_col_dbl(winTable->fptr,objDecCol,i+objNum,1,1,null,
			      &(winObj->objDec),&null,&status)) {
          cpl_msg_error("readFitsWindowTable","The function fits_read_col_dbl has returned an error (code %d)", status);
          return (VM_FALSE);
        } 
      }


         
      if (lastObject == NULL) {
        lastObject = winObj;
        winSlit->objs = winObj;
      }
      else {
        lastObject->next = winObj;
        winObj->prev = lastObject;
      }
      lastObject = winObj;
      
      objNum++;
      if ((i + objNum) <= nRows)
      {
	if (fits_read_col_int(winTable->fptr,slitCol,i+objNum,1,1,null,
			      &newSlitNo,&null,&status)) {
	  cpl_msg_error("readFitsWindowTable","The function fits_read_col_int has returned an error (code %d)", status);
	  return (VM_FALSE);
	} 
      }
    }
    winSlit->numObj = objNum;
    
    i += objNum;
    
    if (lastSlit == NULL) {
      winTable->slits = winSlit;
    }
    else {
      lastSlit->next = winSlit;
      winSlit->prev = lastSlit;
    }
    lastSlit = winSlit;
        
  }

  for (ii = 0; ii < 3; ii++){     /* free the memory for the column labels */
    cpl_free( ttype[ii] );
  }

  return (VM_TRUE);
}


/*

  Updates:
  03 Sept 02: added read/write of IFUfibTrans column (AZ)

 */

VimosBool writeFitsWindowTable(VimosWindowTable *winTable, fitsfile *fptr)
{
  int nRows;
  int rowNum;
  int specLong;
  
  int objNum;
  int i;
  int status, nbytes;
  char *ttype[84], *tform[84];
  
  VimosWindowSlit *winSlit;
  VimosWindowObject *winObject;
  
  /* validate input */
  if (winTable == NULL) {
    cpl_msg_error("writeFitsWindowTable","NULL input table");
    return (VM_FALSE);
  }
  if ( strcmp(winTable->name, VM_WIN) ) {
    cpl_msg_error("writeFitsWindowTable","Invalid input table");
    return (VM_FALSE);
  }
  

  nRows = numObjectsInWindowSlit(winTable->slits);
  nRows += numEmptySlitsInWindowSlit(winTable->slits);

  winTable->fptr = fptr;
  status = 0;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "WIN", 0, &status)) {
    if (fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_error("writeFitsWindowTable","The function fits_delete_hdu has returned an error (code %d)", status);
      return (VM_FALSE);
    }
  } else {
    status = 0;
  }

  /* ALEX: added IFUfibTrans => i<=15 instead of 14 */
  
  for (i = 0; i <= 15; i++){      /* allocate space for the column labels */
    ttype[i] = (char *) cpl_malloc(FLEN_VALUE);  /* max label length = 69 */
    /* check if space was allocated */
    if (ttype[i] == NULL) {
      cpl_msg_error("writeFitsWindowTable","Allocation Error");
      return (VM_FALSE);
    }
    tform[i] = (char *) cpl_malloc(FLEN_VALUE);
    /* check if space was allocated */
    if (tform[i] == NULL) {
      cpl_msg_error("writeFitsWindowTable","Allocation Error");
      return (VM_FALSE);
    }
  }
  
  ttype[0] = "SLIT";
  tform[0] = "1J";
  ttype[1] = "IFUSLIT_NO";
  tform[1] = "1J";
  ttype[2] = "IFUFIB_NO";
  tform[2] = "1J";
  ttype[3]="SPEC_LEN";
  tform[3] = "1J";
  ttype[4]="SPEC_START"; 
  tform[4] = "1J";
  ttype[5] = "SPEC_END";
  tform[5] = "1J";
  ttype[6] = "OBJ_START";
  tform[6] = "1J";
  ttype[7]= "OBJ_END";
  tform[7] = "1J";
  ttype[8]= "OBJ_NO";
  tform[8] = "1J";
  ttype[9]="OBJ_POS";
  tform[9] = "1E";
  ttype[10]="OBJ_WIDTH";
  tform[10] = "1E";
  ttype[11]="OBJ_X";
  tform[11] = "1E";
  ttype[12]="OBJ_Y";
  tform[12] = "1E";
  ttype[13]="OBJ_RA";
  tform[13] = "1D";
  ttype[14]="OBJ_DEC";
  tform[14] = "1D";
  /* ALEX: added IFUfibTrans */
  ttype[15]="IFUFIB_TRANS";
  tform[15] = "1E";

  /* append a new empty binary table onto the FITS file */

  /* ALEX: added IFUfibTrans => 16 instead of 15 */
  /* orig.: fits_create_tbl(fptr, BINARY_TBL,0, 15, ttype, tform,NULL, .... */

  if (fits_create_tbl(fptr, BINARY_TBL,0, 16, ttype, tform,NULL, 
		      "WIN", &status)) {
    cpl_msg_error("writeFitsWindowTable","The function fits_create_tbl has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "WIN", 0, &status)) {
    cpl_msg_error("writeFitsWindowTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  /* write the table descriptors to the Fits file */
  
  if (fits_read_key (winTable->fptr,TINT,"NAXIS1",&nbytes,NULL,&status)) {
    cpl_msg_error("writeFitsWindowTable","The function fits_read_key has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(winTable->descs), "NAXIS1", nbytes, "")) {
    cpl_msg_error("writeFitsWindowTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeIntDescriptor(&(winTable->descs), "NAXIS2", nRows, "")) {
    cpl_msg_error("writeFitsWindowTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }

  /* ALEX: added IFUfibTrans => 16 instead of 15 */
  if (!writeIntDescriptor(&(winTable->descs),"TFIELDS",16, "")) {
    cpl_msg_error("writeFitsWindowTable","The function writeIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!writeDescsToFitsTable(winTable->descs, winTable->fptr)) {
    cpl_msg_error("writeFitsWindowTable","The function writeDescsToFitsTable has returned an error");
    return (VM_FALSE);
  }

  winSlit = winTable->slits;
  rowNum = 1;
  while (winSlit) {
    winObject = winSlit->objs;
    objNum = 1;
    if (winObject == NULL) 
      winObject = newWindowObject();
    
    while (winObject) {
      if (fits_write_col_int(winTable->fptr,1,rowNum,1,1,&(winSlit->slitNo),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE);   
      }
      if (fits_write_col_int(winTable->fptr,2,rowNum,1,1,&(winSlit->IFUslitNo),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE); 
      }  
      if (fits_write_col_int(winTable->fptr,3,rowNum,1,1,&(winSlit->IFUfibNo),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE);   
      }

      if (winSlit->specLong == VM_TRUE) {
        specLong = 1;
      }
      else {
        specLong = 0;
      }
      if (fits_write_col_int(winTable->fptr,4,rowNum,1,1,&specLong,&status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE);   
      }
      if (fits_write_col_int(winTable->fptr,5,rowNum,1,1,&(winSlit->specStart),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE);   
      }
      if (fits_write_col_int(winTable->fptr,6,rowNum,1,1,&(winSlit->specEnd),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE);     
      }   
      if (fits_write_col_int(winTable->fptr,7,rowNum,1,1,
			     &(winObject->objStart),&status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE);
      }        
      if (fits_write_col_int(winTable->fptr,8,rowNum,1,1,&(winObject->objEnd),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE);       
      }
      if (fits_write_col_int(winTable->fptr,9,rowNum,1,1,&(objNum),&status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_int has returned an error (code %d)", status);
        return (VM_FALSE);  
      } 
      if (fits_write_col_flt(winTable->fptr,10,rowNum,1,1,&(winObject->objPos),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return (VM_FALSE); 
      }

/****  Remove objWidth

      if (fits_write_col_flt(winTable->fptr,11,rowNum,1,1,&(winObject->objWidth),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_flt has returned an error (code %d)", status);
        return (VM_FALSE); 
      }  

****/

      if (fits_write_col_flt(winTable->fptr,12,rowNum,1,1,&(winObject->objX),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_dbl has returned an error (code %d)", status);
        return (VM_FALSE);   
      }
      if (fits_write_col_flt(winTable->fptr,13,rowNum,1,1,&(winObject->objY),
			     &status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_dbl has returned an error (code %d)", status);
        return (VM_FALSE);   
      }

      if (winObject->posDef) {
	if (fits_write_col_dbl(winTable->fptr,14,rowNum,1,1,
			       &(winObject->objRA),&status)) {
          cpl_msg_error("writeFitsWindowTable","The function fits_write_col_dbl has returned an error (code %d)", status);
	  return (VM_FALSE);   
	}
	if (fits_write_col_dbl(winTable->fptr,15,rowNum,1,1,
			       &(winObject->objDec),&status)) {
          cpl_msg_error("writeFitsWindowTable","The function fits_write_col_dbl has returned an error (code %d)", status);
	  return (VM_FALSE);   
	}
      }
      else { 
        /* should write a null value */
	if (fits_write_col_dbl(winTable->fptr,14,rowNum,1,1,
			       &(winObject->objRA),&status)) {
          cpl_msg_error("writeFitsWindowTable","The function fits_write_col_dbl has returned an error (code %d)", status);
	  return (VM_FALSE);
	}   
	if (fits_write_col_dbl(winTable->fptr,15,rowNum,1,1,
			       &(winObject->objDec),&status)) {
          cpl_msg_error("writeFitsWindowTable","The function fits_write_col_dbl has returned an error (code %d)", status);
	  return (VM_FALSE);
	}
      }

      /* ALEX: added IFUfibTrans */
      if (fits_write_col_flt(winTable->fptr,16,rowNum,1,1,
			     &(winSlit->IFUfibTrans),&status)) {
        cpl_msg_error("writeFitsWindowTable","The function fits_write_col_flt has returned an error (code %d)",status);
        return (VM_FALSE);   
      }
      
      rowNum++;
      objNum++;
      winObject = winObject->next;
    }

    winSlit = winSlit->next;
  }
  
  return(VM_TRUE);
}
