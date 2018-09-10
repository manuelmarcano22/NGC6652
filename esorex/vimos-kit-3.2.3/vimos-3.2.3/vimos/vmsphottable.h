/* $Id: vmsphottable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_SPHOTTABLE_H
#define VM_SPHOTTABLE_H

#include <pilmacros.h>

#include <vmtable.h>


PIL_BEGIN_DECLS

#define    VM_SPHOT     "SPH"

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosTable *newSphotTable()

   Description:
   Allocate a new SpectroPhotometric Table. This is a VimosTable with 
   name VM_SPHOT (which has value "SPH").
   Only a Descriptor "TABLE" is present, with value "SPH".
   
   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosTable

   Return Value (error):
   function calls exit(-1)
   
   Updates:
   23 Jul 02: Created (PF)

--------------------------------------------------------------------------------
*/
VimosTable  *newSphotTable(/* void */);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteSphotTable(VimosTable *sTab)
   
   Description:
   Delete a SpectroPhotometric Table. This is just an esthetic wrapper for
   deleteTable(sTab)
   
   Input: 
   VimosTable *sTab
   Pointer to SpectroPhotometric Table to be deleted
   
   Return Value:
   void
   
   Updates:
   22 Jul 02: Created (PF)

-------------------------------------------------------------------------------- */
void        deleteSphotTable(VimosTable *sTab);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool readFitsSphotTable(VimosTable *sphotTable, fitsfile *fptr)

   Description:
   Read a FITS Table into a SpectroPhotometric Table. This requires 
   the FITS file to be already open, and does not close it at the end.
   
   Input:
   VimosTable *sphotTable 
   Pointer to the SpectroPhotometric Table to copy the FITS table into.

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool (true if the read is succesful, false otherwise)
   
   Updates:
   23 Jul 02: Created (PF)

--------------------------------------------------------------------------------
*/
VimosBool readFitsSphotTable(VimosTable *sphotTable, fitsfile *fptr);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool writeFitsSphotTable(VimosTable *sphotTable, fitsfile *fptr)

   Description:
   Write a SpectroPhotometric Table into a FITS Table. This requires the FITS file
   to be already open, and does not close it at the end.
   
   Input:
   VimosTable *sphotTable  
   Pointer to the SpectroPhotometric Table to write to a FITS Table

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool (true if the write is succesful, false otherwise)
   
   Updates:
   23 Jul 02: Created (PF)

--------------------------------------------------------------------------------
*/
VimosBool writeFitsSphotTable(VimosTable *sphotTable, fitsfile *fptr);

PIL_END_DECLS

#endif /* VM_SPHOTTABLE_H */
