/* $Id: vmlinecatalog.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_LINECATALOG_H
#define VM_LINECATALOG_H

#include <pilmacros.h>

#include <vmtable.h>


PIL_BEGIN_DECLS

#define    VM_LIN     "LIN"

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosTable *newLineCatalog()

   Description:
   Allocate a new Line Catalog. This is a VimosTable with name VM_LIN (which
   has value "LIN").
   Only a Descriptor "TABLE" is present, with value "LIN".
   
   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosTable

   Return Value (error):
   function calls exit(-1)
   
   Updates:
   04 Feb 99: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosTable  *newLineCatalog(void);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteLineCatalog(VimosTable *lCat)
   
   Description:
   Delete a Line Catalog. This is just an esthetic wrapper for
   deleteTable(lCat)
   
   Input: 
   VimosTable *lCat
   Pointer to Line Catalog to be deleted
   
   Return Value:
   void
   
   Updates:
   04 Feb 99: Created (TAO)

-------------------------------------------------------------------------------- */
void        deleteLineCatalog(VimosTable *lCat);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool readFitsLineCatalog(VimosTable *linCat, fitsfile *fptr)

   Description:
   Read a FITS Table into a Line Catalog. This requires the FITS file to 
   be already open, and does not close it at the end.
   
   Input:
   VimosTable *linCat  
   Pointer to the Line Catalog to copy the FITS table into.

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
VimosBool readFitsLineCatalog(VimosTable *linCat, fitsfile *fptr);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool writeFitsLineCatalog(VimosTable *linCat, fitsfile *fptr)

   Description:
   Write a Line Catalog into a FITS Table. This requires the FITS file
   to be already open, and does not close it at the end. If an Extraction
   Table extension is already present in the file, it is first removed, and
   then the new one is written into a new extension.  
   
   Input:
   VimosTable *linCat  
   Pointer to the Line Catalog to write to a FITS Table

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
VimosBool writeFitsLineCatalog(VimosTable *linCat, fitsfile *fptr);

PIL_END_DECLS

#endif /* VM_LINECATALOG_H */
