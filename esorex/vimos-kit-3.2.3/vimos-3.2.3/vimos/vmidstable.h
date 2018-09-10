/* $Id: vmidstable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_IDSTABLE_H
#define VM_IDSTABLE_H

#include <pilmacros.h>

#include <vmtable.h>


PIL_BEGIN_DECLS

#define    VM_IDS     "IDS"

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosTable *newIdsTable()

   Description:
   Allocate a new IDS Table. This is a VimosTable with name VM_IDS (which
   has value "IDS").
   Only a Descriptor "TABLE" is present, with value "IDS".
   
   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosTable

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: Call exit(-1) instead of return(NULL) on error (TAO)
   03 Nov 98: Created (TAO)
   10 Nov 98: Changed to fit new table structure (TAO)

--------------------------------------------------------------------------------
*/
VimosTable  *newIdsTable(void);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteWindowTable(VimosTable *gTable)
   
   Description:
   Delete a IDS Table. This is just an esthetic wrapper for
   deleteTable(gTable)
   
   Input: 
   VimosTable *gTable
   Pointer of IDS Table to be deleted
   
   Return Value:
   void
   
   Updates:
   10 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- */
void        deleteIdsTable(VimosTable *gTable);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool readFitsIdsTable(VimosTable *idsTable, fitsfile *fptr)
   Description:
   Read a FITS Table into a IDS Table. This function is just an esthetic
   wrapper for readMidasTable, except that it checks that the type of the
   VimosTable is VM_IDS (=="IDS"). It requires the FITS file to 
   be already open, and does not close it at the end.
   
   Input:
   VimosTable *idsTable  
   Pointer to the IDS Table to copy the FITS table into.

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
VimosBool readFitsIdsTable(VimosTable *idsTable, fitsfile *fptr);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool writeFitsIdsTable(VimosTable *idsTable, fitsfile *fptr)

   Description:
   Write a IDS Table into a FITS Table. This function is just an esthetic
   wrapper for writeFitsTable, except that it checks that the type of the
   VimosTable is VM_IDS (=="IDS"). It requires the FITS file
   to be already open, and does not close it at the end. If a IDS
   Table extension is already present in the file, it is first removed, and
   then the new one is written into a new extension.

   Input:
   VimosTable *idsTable  
   Pointer to the IDS Table to write to a FITS Table

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
VimosBool writeFitsIdsTable(VimosTable *idsTable, fitsfile *fptr);

PIL_END_DECLS

#endif /* VM_IDSTABLE_H */
