/* $Id: vmobjecttable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_OBJECTTABLE_H
#define VM_OBJECTTABLE_H

#include <fitsio.h>

#include <pilmacros.h>

#include <vmwindowtable.h>


PIL_BEGIN_DECLS

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   defines of name strings for the Object Tables

   Description:
   The following strings are defined to use in the Oject Tables for
   identifying the name of the table and other things

   Values:
   VM_OBJ           string giving the name of the Object Table
   more to follow.......
   
   Updates:
   11 Feb 99: Created (TAO)

--------------------------------------------------------------------------------
*/
#define VM_OBJ           "OBJ"
/*  more to follow ... */

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   stucture VimosObjectObject
   
   Description: 
   Linked list structure that contains the information of one object
   in an Object Table. An Object Table is simply a listof VimosObjectObjects
   plus descriptors. 

   Layout:
     int               slitNo;
     int               IFUslitNo;
     int               IFUfibNo;
     int               rowNum;
     int               objNo;
     float             objX;
     float             objY;
     double            objRA;
     double            objDec;
     VimosObjectObject *prev;
     VimosObjectObject *next;

   Updates:
   11 Feb 99: Created (TAO)
   03 Feb 00: Added IFU slit and fiber bumber information (MS)

--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_OBJECT_OBJ_
{
  int                       slitNo;
  int                       IFUslitNo;
  int                       IFUfibNo;
  int                       rowNum;
  int                       objNo;
  float                     objX;
  float                     objY;
  double                    objRA;
  double                    objDec;
  /**/  
  struct _VIMOS_OBJECT_OBJ_ *prev;
  struct _VIMOS_OBJECT_OBJ_ *next;
} VimosObjectObject;


typedef struct _VIMOS_OBJECT_TABLE_ 
{
  char            name[VM_DESC_LENGTH];
  VimosDescriptor *descs;
  VimosObjectObject *objs;
  fitsfile             *fptr;                  
} VimosObjectTable;

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosWindowObject  *newObjectObject()
   
   Description: 
   Returns a pointer to a new VimosObjectObject structure.
   
   Input:
   void   
   
   Return Value (succes):
   Pointer to a newly allocated VimosobjectObject structure. All fields are
   initialized to 0

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   11 Feb 99: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosObjectObject  *newObjectObject();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteObjectObject(VimosObjectObject *oObject)
   
   Description: 

   Deletes all VimosObjectObejcts contained in the linked list oObject
   
   Input:
   VimosObjectObject *oObject
   Pointer to VimosObjectObject list to be deleted
   
   Return Value:
   void
   
   Updates:
   11 Feb 99: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteObjectObject(VimosObjectObject *oObject);





/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosObjectTable *newObjectTable()

   Description:
   Allocate a new Object Table. All fields are initialized to 0, except the
   name field is set to VM_OBJ
   
   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosObjectTable

   Return Value (error):
   NULL
      
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   11 Feb 99: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosObjectTable *newObjectTable();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteObjectTable(VimosObjectTable *oTable)
   
   Description:
   Delete a Object Table. 
   
   Input: 
   VimosObjectTable *oTable
   Pointer of Object Table to be deleted
   
   Return Value:
   void
   
   Updates:
   11 feb 99: Created (TAO)

-------------------------------------------------------------------------------- 
*/
void deleteObjectTable(VimosObjectTable *oTable);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool copyWinTab2ObjTab(VimosWindowTable *winTable, 
                               VimosObjectTable *objTable)

   Description:
   Copy ALL descriptors (including those that may not be relevant) of an
   Window Table into Object Table. 
   
   Input:
   VimosWindowTable *winTable  
   Pointer to the Window Table to copy the descriptors from

   VimosObjectTable *objTable  
   Pointer to the Object Table to copy the descriptors to
   
   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE

   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   11 Feb 99: Created (TAO)

--------------------------------------------------------------------------------
*/

VimosBool copyWinTab2ObjTab(VimosWindowTable *winTable, 
                            VimosObjectTable *objTable);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   int numObjectsInObjectObject(VimosObjectObject *obj)
   
   Description:
   Return number of object in a VimosObjectObject. 
   
   Input:
   VimosObjectObject *obj  
   Pointer to the Object . 
   
   Return Value:
   Number of objects in the VimosObjectObject 
   
   Updates:
   11 Feb 99: Created (TAO)
   
--------------------------------------------------------------------------------
*/
int numObjectsInObjectObject(VimosObjectObject *obj);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool readFitsObjectTable(VimosObjectTable *objTable, fitsfile *fptr)

   Description:
   Read a FITS Table into an Object Table. This requires the FITS file to 
   be already open, and does not close it at the end.
   
   Input:
   VimosObjectTable *objTable  
   Pointer to the Object Table to copy the FITS table into.

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool (true if the read is succesful, false otherwise)
   
   Updates:
   15 Dec 99: Created (BG)
   02 Feb 00: Changed the handling of the table, from a standalone FITS file,
              to an extension of the primary FITS file (MS)

--------------------------------------------------------------------------------
*/
VimosBool readFitsObjectTable(VimosObjectTable *objTable, fitsfile *fptr);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool writeFitsObjectTable(VimosObjectTable *winTable, fitsfile *fptr)

   Description:
   Write a Object Table into a FITS Table. This requires the FITS file
   to be already open, and does not close it at the end. If an Extraction
   Table extension is already present in the file, it is first removed, and
   then the new one is written into a new extension.  
   
   Input:
   VimosObjectTable *objTable  
   Pointer to the Object Table to write to a FITS Table

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool (true if the write is succesful, false otherwise)
   
   Updates:
   15 Dec 99: Created (BG)
   02 Feb 00: Changed the handling of the table, from a standalone FITS file,
              to an extension of the primary FITS file (MS)

--------------------------------------------------------------------------------
*/
VimosBool writeFitsObjectTable(VimosObjectTable *objTable, fitsfile *fptr);

PIL_END_DECLS

#endif /* VM_OBJECTTABLE_H */
