/* $Id: vmifutable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_IFUTABLE_H
#define VM_IFUTABLE_H

#include <pilmacros.h>

#include <vmtable.h>


PIL_BEGIN_DECLS
/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   defines of name strings for the IFU Table
 
   Description:
   The following strings are defined to use in the IFU Table for
   identifying the name of the table and of the columns
 
   Values:  
   VM_IFU          string giving the name of the IFU Table
   VM_IFU_x        name of column for fiber x coordinate on mask
   VM_IFU_y        name of column for fiber y coordinate on mask
   VM_IFU_L        name of column for fiber L coordinate on microlens
   VM_IFU_M        name of column for fiber M coordinate on microlens
   VM_IFU_PWIDTH   name of column for FWHM of spectrum on CCD (spatial
                   direction) in frame format
   VM_IFU_QUAD     name of column for vimos quadrant of mask this fiber is on
   VM_IFU_ROW      name of column for slit of fiber in mask plane
   VM_IFU_FIB      name of sequential number of fiber in each slit
   VM_IFU_TRANS    name of column for relative transmission of fiber with
                   respect to reference fiber
   VM_IFU_SIGMAY   name of column for sigma of spectrum on CCD in dispersion
                   direction (Y)
   VM_IFU_SIGMAYGROUP   name of column for grouping spectra on each slit 
                        according to their values of sigmaY (needed for sky
			spectrum determination)


   Updates:
   22 Nov 99: Created (AZ)
 
--------------------------------------------------------------------------------
*/

#define VM_IFU      "IFU"
/* Column names */
#define VM_IFU_x        "x"
#define VM_IFU_y        "y"
#define VM_IFU_L        "L"
#define VM_IFU_M        "M"
#define VM_IFU_PWIDTH   "PWIDTH"
#define VM_IFU_QUAD     "QUAD"
#define VM_IFU_ROW      "ROW"
#define VM_IFU_FIB      "FIB"
#define VM_IFU_TRANS    "TRANS"
#define VM_IFU_SIGMAY      "SIGMAY"
#define VM_IFU_SIGMAYGROUP "SIGMAYGROUP"

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  stucture VimosIfuFiber

  Description:
  Linked list structure that contains the information of the fibers
  defined for each IFU slit on each mask.

  Layout:
     int   fibNo
     int   fiberL
     int   fiberM
     float fiberx
     float fibery
     float fiberTrans
     float fiberPwidth
     float sigmaY
     int   sigmaYGroup
     VimosIfuFiber *prev
     VimosIfuFiber *next

   Updates:
   22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

typedef struct _VIMOS_IFU_FIBER_
{
  int   fibNo;
  int   fiberL;
  int   fiberM;
  float fiberx;
  float fibery;
  float fiberTrans;
  float fiberPwidth;
  float sigmaY;
  int   sigmaYGroup;
  struct _VIMOS_IFU_FIBER_ *prev;
  struct _VIMOS_IFU_FIBER_ *next;
} VimosIfuFiber;


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  stucture VimosIfuSlit

  Description:
  Linked list structure that contains the information of the slits defined
  for each IFU quadrant.

  Layout:
     int   ifuSlitNo
     VimosIfuFiber *fibers
     VimosIfuSlit  *prev
     VimosIfuSlit  *next

   Updates:
   22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

typedef struct _VIMOS_IFU_SLIT_
{
  int   ifuSlitNo;
  VimosIfuFiber *fibers;
  struct _VIMOS_IFU_SLIT_  *prev;
  struct _VIMOS_IFU_SLIT_  *next;
} VimosIfuSlit;



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  structure VimosIfuQuad

  Description:
  Linked list structure that contains the information of the
  vimos quadrants.

  Layout:
     int           quadNo
     VimosIfuSlit  *ifuSlits
     int           numIfuSlits
     VimosIfuQuad  *prev
     VimosIfuQuad  *next

  Updates:
  22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

typedef struct _VIMOS_IFU_QUAD_
{
  int           quadNo;
  VimosIfuSlit  *ifuSlits;
  int numIfuSlits;
  struct _VIMOS_IFU_QUAD_  *prev;
  struct _VIMOS_IFU_QUAD_  *next;
} VimosIfuQuad;


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  structure VimosIfuTable  

  Description:
  Here the true structure of a vimos IFU Table is defined.

  Layout:
     
  char            name
  VimosDescriptor *descs
  VimosIfuQuad    *quads
  int             numIfuQuads
  int             numIfufibs
  fitsfile       *fptr

  Updates:
  22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

typedef struct _VIMOS_IFU_TABLE_
{
  char            name[VM_DESC_LENGTH];
  VimosDescriptor *descs;
  VimosIfuQuad    *quads;
  int             numIfuQuads;
  int             numIfuFibs;
  fitsfile             *fptr;                  
} VimosIfuTable;






/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  VimosIfuFiber *newIfuFiber()

  Description: 
  Returns a pointer to a new IfuFiber structure.

  Input:
  void   

   Return Value (succes):
   Pointer to a newly allocated IfuFiber structure. All fields are
   initialized to 0

   Return Value (error):
   NULL

   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

VimosIfuFiber *newIfuFiber();

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
  void deleteIfuFiber(VimosIfuFiber *ifuFiber)

  Description:
  Deletes all VimosIfuFibers contained in the linked list ifuFiber.

  Input:
  VimosIfuFiber *ifuFiber
  Pointer to VimosIfuFiber list to be deleted

  Return Value:
  void
   
  Updates:
  22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

 void deleteIfuFiber(VimosIfuFiber *ifuFiber);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosIfuSlit  *newIfuSlit()

   Description: 
   Returns a pointer to a new VimosIfuSlit structure.
   
   Input:
   void   
   
   Return Value (success):
   Pointer to a newly allocated VimosIfuSlit structure. All fields are
   initialized to 0
 
   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

  VimosIfuSlit  *newIfuSlit();

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteIfuSlit(VimosIfuSlit *ifuSlit);

   Description: 
   Deletes all VimosIfuSlits (and VimosIfuFibers contained in the 
   VimosIfuSlit) contained in the linked list ifuSlit
   
   Input:
   VimosIfuSlit *ifuSlit
   Pointer to VimosIfuSlit list to be deleted

   Return Value:
   void
   
   Updates:
   22 Nov 99: Created (AZ)   

--------------------------------------------------------------------------------
*/

void deleteIfuSlit(VimosIfuSlit *ifuSlit);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosIfuQuad  *newIfuQuad()

   Description: 
   Returns a pointer to a new VimosIfuQuad structure.
   
   Input:
   void   
   
   Return Value (success):
   Pointer to a newly allocated VimosIfuQuad structure. All fields are
   initialized to 0
 
   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

  VimosIfuQuad  *newIfuQuad();

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteIfuQuad(VimosIfuQuad *ifuQuad);

   Description: 
   Deletes all VimosIfuQuads (and VimosIfuSlits, VimosIfuFibers contained
   in the VimosIfuQuad) contained in the linked list ifuQuad
   
   Input:
   VimosIfuQuad *ifuQuad
   Pointer to VimosIfuQuad list to be deleted

   Return Value:
   void
   
   Updates:
   22 Nov 99: Created (AZ)   

--------------------------------------------------------------------------------
*/

void deleteIfuQuad(VimosIfuQuad *ifuQuad);

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   VimosIfuTable *newIfuTable()

   Description:
   Allocate a new IFU Table. All fields are initialized to 0, except the
   name field is set to VM_IFU

   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosIfuTable

   Return Value (error):
   NULL
      
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

  VimosIfuTable *newIfuTable();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   void deleteIfuTable(VimosIfuTable *ifuTable)

   Description:
   Delete an IFU Table.


   Input: 
   VimosIfuTable *ifuTable
   Pointer of IFU Table to be deleted
   
   Return Value:
   void
   
   Updates:
   22 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

   void deleteIfuTable(VimosIfuTable *ifuTable);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   VimosIfuSlit *computeIfuslit(int l_start, int m_start, int l_step, 
                               int m_step, int module_step_m, 
                               float x_start, float x_step, float y,
                               float x_module_step)


   Description:
   Computes (L,M) and (X,Y) coordinates of each fiber in a given slit

   Input: 
//   VimosIfuSlit *ifuSlit
   Pointer of IFU Slit
   
   Return Value (succes):
   VimosIfuSlit *ifuSlit

   Return Value (error):
   NULL 
   
   Updates:
   13 Jun 00: Return NULL on error (Maura)
   23 Nov 99: Created (AZ)

--------------------------------------------------------------------------------
*/

VimosIfuSlit *computeIfuSlit(int l_start, int m_start, int l_step, 
                             int m_step, int module_step_m,
                             float x_start, float x_step, float y,
                             float x_module_step);


VimosBool writeTable(VimosIfuTable *anIfuTable);

/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   VimosBool readFitsIfuTable(VimosIfuTable *ifuTable, fitsfile *fptr)

   Description:
   Read a FITS Table into an IFU Table. This requires the FITS file to 
   be already open, and does not close it at the end.

   Input:
   VimosIfuTable *ifuTable
   Pointer to the IFU Table to copy the FITS Table into.

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool

   Updates:
   11 Jan 00: Created (AZ)
   02 Feb 00: Changed the handling of the table, from a standalone FITS file,
              to an extension of the primary FITS file (MS)

-------------------------------------------------------------------------------
*/

VimosBool readFitsIfuTable(VimosIfuTable *ifuTable, fitsfile *fptr);


/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   VimosBool writeFitsIfuTable(VimosIfuTable *ifuTable, fitsfile *fptr)

   Description:
   Write an IFU Table into a FITS Table. This requires the FITS file
   to be already open, and does not close it at the end. If an IFU
   Table extension is already present in the file, it is first removed, and
   then the new one is written into a new extension.

   Input:
   VimosIfuTable *ifuTable
   Pointer to the IFU Table to write to a FITS Table.

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool

   Updates:
   11 Jan 00: Created (AZ)
   02 Feb 00: Changed the handling of the table, from a standalone FITS file,
              to an extension of the primary FITS file (MS)

-------------------------------------------------------------------------------
*/

VimosBool writeFitsIfuTable(VimosIfuTable *ifuTable, fitsfile *fptr);

PIL_END_DECLS

#endif /* VM_IFUTABLE_H */
