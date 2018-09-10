/* $Id: vmwindowtable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_WINDOWTABLE_H
#define VM_WINDOWTABLE_H

#include <fitsio.h>

#include <pilmacros.h>

#include <vmextractiontable.h>


PIL_BEGIN_DECLS

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   defines of name strings for the Window Tables

   Description:
   The following strings are defined to use in the Window Tables for
   identifying the name of the table and of the columns

   Values:
   VM_WIN           string giving the name of the Extraction Table
   VM_WIN_SLIT      name of slit ID column
   VM_WIN_SPC_LEN   name of column of spectrum length
   VM_WIN_OBJ_NO    name of column of object number
   more to follow.......
   
   Updates:
   03 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
#define VM_WIN           "WIN"
/* Column names */
#define VM_WIN_SLIT      "SLIT"
#define VM_WIN_SPC_LEN   "SPC LEN" 
#define VM_WIN_OBJ_NO    "OBJ NO"

/*  more to follow ... */

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   stucture VimosWindowObject
   
   Description: 
   Linked list structure that contains the information of objects in one slit, 
   as needed by the Window Table. VimosWindowObjects are pointed to in a 
   VimosWindowSlit

   Layout:
     int               objStart;
     int               objEnd;
     int               objNo;
     float             objPos;
     float             objWidth;
     double            objX;
     double            objY;
     VimosFloatArray   *objProfile;
     Bool              posDef;
     double            objRA;
     double            objDec;
     Bool              parDef;
     VimosWindowObject *prev;
     VimosWindowObject *next;

   Updates:
   09 Nov 98: Created (TAO)
   03 Feb 00: Removed objWidth, objPeak, skyLevel, as they are not used (MS)
   12 Aug 02: Addes objWidth  and objProfile (PF)

--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_WINDOW_OBJ_
{
  int                       objStart;
  int                       objEnd;
  int                       objNo;
  float                     objPos;
  float                     objWidth;
  float                     objX;
  float                     objY;
  VimosFloatArray           *objProfile;
  VimosBool                 posDef;
  double                    objRA;
  double                    objDec;
  VimosBool                 parDef;
  /*  float                     objPeak;
      float                     skyLevel; */
  /**/  
  struct _VIMOS_WINDOW_OBJ_ *prev;
  struct _VIMOS_WINDOW_OBJ_ *next;
} VimosWindowObject;

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   structure VimosWindowSlit
   
   Description: 
   Linked list structure that contains the information of the slits defined in 
   a Vimos spectroscopic observation.

   Layout:
      int              slitNo;
      int              IFUslitNo;
      int              IFUfibNo;
      float            IFUfibTrans;
      VimosBool        specLong;
      int              specStart;
      int              specEnd;
      VimosWindowSlit  *prev;
      VimosWindowSlit  *next;
      int               numObj;
      VimosWindowObject *objs;

   Updates:
   09 Nov 98: Created (TAO)
   03 Feb 00: Added IFU slit and fiber number information (MS)
   03 Sept 02: Added IFU fib. trans column, defaulted to zero for MOS, will
               be updated by VmIfuWindowTable for IFU (AZ)

--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_WINDOW_SLIT_
{
  int        slitNo;
  int        IFUslitNo;
  int        IFUfibNo;
  float      IFUfibTrans;
  VimosBool  specLong;
  int        specStart;
  int        specEnd;
  struct     _VIMOS_WINDOW_SLIT_ *prev;
  struct     _VIMOS_WINDOW_SLIT_ *next;
  /* Number of objects in this slit */
  int        numObj;
  VimosWindowObject *objs;
} VimosWindowSlit;



typedef struct _VIMOS_WINDOW_TABLE_ 
{
  char            name[VM_DESC_LENGTH];
  VimosDescriptor *descs;
  VimosWindowSlit *slits;
  fitsfile             *fptr;                  
} VimosWindowTable;

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosWindowObject  *newWindowObject()
   
   Description: 
   Returns a pointer to a new WindowObject structure.
   
   Input:
   void   
   
   Return Value (succes):
   Pointer to a newly allocated WindowObject structure. All fields are
   initialized to 0

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   09 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosWindowObject  *newWindowObject();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteWindowObject(VimosWindowObject *wObject)
   
   Description: 

   Deletes all VimosWindowObejcts contained in the linked list wObject
   
   Input:
   VimosWindowObject *wObject
   Pointer to VimosWindowObject list to be deleted
   
   Return Value:
   void
   
   Updates:
   09 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteWindowObject(VimosWindowObject *wObj);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosWindowSlit  *newWindowSlit()
   
   Description: 
   Returns a pointer to a new VimosWindowSlit structure.
   
   Input:
   void   
   
   Return Value (success):
   Pointer to a newly allocated VimosWindowSlit structure. All fields are
   initialized to 0

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   09 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosWindowSlit  *newWindowSlit();





/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteWindowSlit(VimosWindowSlit *wSlit)
   
   Description: 

   Deletes all VimosWindowSlits (and VimosWindowObjects contained in the 
   VimosWindowSlits) contained in the linked list wSlit
   
   Input:
   VimosWindowSlit *wSlit
   Pointer to VimosWindowSlit list to be deleted
   
   Return Value:
   void
   
   Updates:
   09 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteWindowSlit(VimosWindowSlit *wSlit);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosWindowTable *newWindowTable()

   Description:
   Allocate a new Window Table. All fields are initialized to 0, except the
   name field is set to VM_WIN
   
   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosWindowTable

   Return Value (error):
   NULL
      
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   20 Nov 98: Adapted for new declaration of VimosWindowTable
   06 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosWindowTable *newWindowTable();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteWindowTable(VimosWindowTable *wTable)
   
   Description:
   Delete a Window Table. 
   
   Input: 
   VimosWindowTable *wTable
   Pointer of Window Table to be deleted
   
   Return Value:
   void
   
   Updates:
   09 Nov 98: Created (TAO)
   20 Nov 98: Adapted for new declaration of VimosWindowTable

-------------------------------------------------------------------------------- 
*/
void deleteWindowTable(VimosWindowTable *wTable);





/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool copyExtTab2WinTab(VimosExtractionTable *extTable, 
                               VimosWindowTable *winTable)

   Description:
   Copy ALL descriptors (including those that may not be relevant) of an
   Extraction Table into Window Table. 
   
   Input:
   VimosTable *extTable  
   Pointer to the Extraction Table to copy the descriptors from

   VimosWindowTable *winTable  
   Pointer to the Window Table to copy the descriptors to
   
   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   06 Nov 98: Created (TAO)
   20 Nov 98: Adapted for new declaration of VimosWindowTable

--------------------------------------------------------------------------------
*/
VimosBool copyExtTab2WinTab(VimosExtractionTable *extTable, 
                            VimosWindowTable *winTable);




/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   int numSlitsInWindowTable(VimosWindowTable *wTable)

   Description:
   Return number of slits in a Window Table. 
   
   Input:
   VimosWindowTable *wtable  
   Pointer to the Window Table. 
   
   Return Value:
   Number of slits in the Window Table
   
   Updates:
   06 Nov 98: Created (TAO)
   20 Nov 98: Adapted for new declaration of VimosWindowTable

--------------------------------------------------------------------------------
*/
int numSlitsInWindowTable(VimosWindowTable *wTable);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   int numObjsInWindowTable(VimosWindowTable *wTable)
   
   Description:
   Return number of objects in a Window Table. 
   
   Input:
   VimosWindowTable *wtable  
   Pointer to the Window Table. 
   
   Return Value:
   Number of objects in the Window Table
   
   Updates:
   09 Nov 98: Created (TAO)
   20 Nov 98: Adapted for new declaration of VimosWindowTable
   
--------------------------------------------------------------------------------
*/
int numObjsInWindowTable(VimosWindowTable *wTable);




/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   int numObjectsInObject(VimosWindowObject *wObject)

   Description:
   Return number of objects in a Window Object. Since a Window Object is a
   linked list of Window Objects, more than one objects can be contained in the
   list.
   
   Input:
   VimosWindowObject *wObject
   Pointer to the Window Object for which the number of objects has to be 
   counted.
   
   Return Value:
   Number of objects in the Window Object
   
   Updates:
   09 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- */
int numObjectsInWindowObject(VimosWindowObject *wObject);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   int numObjectInWindowSlit(VimosWindowSlit *wSlit)

   Description:

   Return number of objects in a Window Slit. Since a Window Slit is a 
   linked list of Window Slits, more than one slit can be contained in the
   list. The sum of all objects in all slits is returned.
   
   Input:
   VimosWindowSlit *wSlit
   Pointer to the Window Slit for which the number of objects has to be 
   counted.
   
   Return Value:
   Number of objects in the Window Slit
   
   Updates:
   09 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
int numObjectsInWindowSlit(VimosWindowSlit *wSlit);





/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   int numSlitsInWindowSlit(VimosWindowSlit *wSlit)

   Description:
   Return number of slits in a Window Slit. Since a Window Slit is a linked 
   list of Window Slits, more than one slit can be contained in the list.
   
   Input:
   VimosWindowSlit *wSlit
   Pointer to the Window Slit for which the number of slits has to be counted.
   
   Return Value:
   Number of slits in the Window Slit
   
   Updates:
   09 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
int numEmptySlitsInWindowSlit(VimosWindowSlit *wSlit);

int numSlitsInWindowSlit(VimosWindowSlit *wSlit);

VimosBool windowObjectInRow(VimosWindowSlit *slit, int row);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

VimosBool shiftWindowObjects(VimosWindowTable *refWinTable, 
                             VimosWindowTable *winTable,  float offSet);

   Description:
   Add an offSet to objects in refWinTable and put them in winTable 
   
   Input:
   VimosWindowTable *refWinTable -
      Pointer to the Window Table to be used as reference
   VimosWindowTable *winTable -
      Pointer to the Window Table in which the 'shifted' window objects
      are put 
   float offSet -
      offset to be added

   Return Value:
   VM_TRUE
   
   Updates:
   12 Aug 02 Created (PF)

--------------------------------------------------------------------------------
*/

VimosBool shiftWindowObjects(VimosWindowTable *refWinTable, VimosWindowTable *winTable,  float offSet);

/*
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool readFitsWindowTable(VimosWindowTable *winTable, fitsfile *fptr)

   Description:
   Read a FITS Table into a Window Table. This requires the FITS file to 
   be already open, and does not close it at the end.
   
   Input:
   VimosWindowTable *winTable  
   Pointer to the Window Table to copy the FITS table into.

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
VimosBool readFitsWindowTable(VimosWindowTable *winTable, fitsfile *fptr);



/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool writeFitsWindowTable(VimosWindowTable *winTable, fitsfile *fptr)

   Description:
   Write a Window Table into a FITS Table. This requires the FITS file
   to be already open, and does not close it at the end. If an Extraction
   Table extension is already present in the file, it is first removed, and
   then the new one is written into a new extension.  
   
   Input:
   VimosWindowTable *winTable  
   Pointer to the Window Table to write to a FITS Table

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
VimosBool writeFitsWindowTable(VimosWindowTable *winTable, fitsfile *fptr);

PIL_END_DECLS

#endif /* VM_WINDOWTABLE_H */
