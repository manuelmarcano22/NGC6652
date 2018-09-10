/* $Id: vmextractiontable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_EXTRACTIONTABLE_H
#define VM_EXTRACTIONTABLE_H

#include <fitsio.h>

#include <pilmacros.h>

#include <vmmatrix.h>
#include <vmifutable.h>
#include <vmdistmodels.h>


PIL_BEGIN_DECLS

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   defines of name strings for the Extraction Tables

   Description:
   The following strings are defined to use in the Extraction Tables for
   identifying the name of the table and of the columns

   Values:
   VM_EXT      string giving the name of the Extraction Table
   
   Updates:
   03 Nov 98: Created (TAO)
   02 Feb 00: Modified the value of the string to "EXR", to avoid conflicts
              with the FITS extensions mechanism (ext and ext...ension) (MS).
--------------------------------------------------------------------------------
*/
#define    VM_EXT     "EXR"



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   Structure VimosExtractionSlit

   Description:
   Structure that more or less corresponds to a row in an Extraction Table

   Layout:
     int                  slitNo
     int                  numRows
     int                  IFUslitNo
     int                  IFUfibNo
     float                IFUfibPeakX
     float                IFUfibTrans
     float                width   (Added for instrument independency (C.Izzo))
     float                *y
     float                *ccdX
     float                *ccdY
     float                *maskX
     float                *maskY
     int                  numSpec
     VimosDistModel1D     **crvPol;
     VimosDistModel1D     **invDis;
     VimosIntArray        invDisQuality
     float                *zeroX
     float                *zeroY
     VimosExtractionSlit  *prev
     VimosExtractionSlit  *next

   Updates:
   25 Jul 02: Added IFUfibTrans field (AZ)
   23 May 02: Added the rms values for the curvature polynomial and
              inverse disperion relation fit (MS)
   23 Oct 01: Added IFUfibPeakX field, used for IFU data reduction (AZ)
   24 Nov 98: Replaced coefficients of curvature polynomial and inverse
              dispersion relation by VimosDistModel1Ds
   20 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_EXTRACTION_SLIT_ 
{
  int                            slitNo;
  int                            numRows;
  int                            IFUslitNo;
  int                            IFUfibNo;
  float                          IFUfibPeakX;
  float                          IFUfibTrans;
  float                          width;
  VimosIntArray                  *y;
  VimosFloatArray                *ccdX;
  VimosFloatArray                *ccdY;
  VimosFloatArray                *maskX;
  VimosFloatArray                *maskY;
  VimosIntArray                  *numSpec;
  VimosDistModel1D               **crvPol;
  VimosFloatArray                *crvPolRms;
  VimosDistModel1D               **invDis;
  VimosFloatArray                *invDisRms;
  VimosIntArray                  *invDisQuality;
  VimosFloatArray                *zeroX;
  VimosFloatArray                *zeroY;
  struct _VIMOS_EXTRACTION_SLIT_ *prev;
  struct _VIMOS_EXTRACTION_SLIT_ *next;
} VimosExtractionSlit;



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   Structure VimosExtractionTable

   Description:
   Structure that more or less corresponds to an Extraction Table

   Layout:
     char                 name[VM_DESC_LENGTH]  always VM_EXT
     VimosDescriptor      *descs
     VimosExtractionSlit  *slits
     fitsfile             *fptr
 

   Updates:
   20 Nov 98: Created (TAO)
   29 Dec 99: Added the pointer to a FITS file, to allow fitsio based
              input/output
--------------------------------------------------------------------------------
*/

typedef struct _VIMOS_EXTRACTION_TABLE_ 
{
  char                 name[VM_DESC_LENGTH];
  VimosDescriptor      *descs;
  VimosExtractionSlit  *slits;
  fitsfile             *fptr;                  
} VimosExtractionTable;



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosExtractionSLit *newExtractionSLit()

   Description:
   Allocate a new Extraction Slit. 
   
   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosExtractionSlit

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: Call exit on error instead of return NULL (TAO)
   23 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosExtractionSlit  *newExtractionSlit();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteExtractionSlit(VimosExtractionSlit *slit)
   
   Description:
   Delete recursively all Extraction Slit in *slit. 
   
   Input: 
   VimosExtarctionSlit *slit
   Pointer of Extraction Slit to be deleted
   
   Return Value:
   void
   
   Updates:
   23 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- 
*/
void deleteExtractionSlit(VimosExtractionSlit *slit);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosExtractionTable *newExtractionTable()

   Description:
   Allocate a new Extraction Table. The field "name" is set to VM_EXT (==
   "EXR"), while also a Descriptor is "TABLE" created with value VM_EXT
   
   Input:
   void
   
   Return Value (success):
   Pointer to the new VimosExtractionTable

   Return Value (error):
   NULL

   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   02 Feb 00: Modified the name from EXT to EXR, to avoid conflicts with the
              FITS extension mechanism (MS)
   07 Dec 98: Added TABLE descriptor
   24 Nov 98: Calls exit(-1) instead of return(NULL) on error (TAO)
   20 Nov 98: Changed to fit new table structure, again (TAO)
   10 Nov 98: Changed to fit new table structure (TAO)
   03 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosExtractionTable  *newExtractionTable();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteExtractionTable(VimosExtractionTable *eTable)
   
   Description:
   Delete a Extraction Table.
   
   Input: 
   VimosExtractionTable *eTable
   Pointer of Extraction Table to be deleted
   
   Return Value:
   void
   
   Updates:
   10 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- 
*/
void deleteExtractionTable(VimosExtractionTable *eTable);

VimosBool copyGrsTab2ExtTab(VimosTable *grsTable, 
                            VimosExtractionTable *extTable);

VimosBool copyAdf2ExtTab(VimosTable *adf, VimosExtractionTable *extTable);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool computeExtractionTable(VimosTable *adf, VimosIfuTable *ifuTab,
                                    VimosTable *extTable)

   Description:
   Using the positions of the slits as defined in the ADf, and using the
   distortion models that are available in the Extraction Table (e.g. the
   parameters of these models been copied from a Grism Table using
   copyGrsTab2ExtTab), the columns of the Extraction Table are computed. 
   Such a computed Extraction Table is used in the on-linepipeline to reduce
   the data, while in the off-line pipeline it can be used as a first
   approximation for the locations of the spectra on the CCD.
   A few descriptors are copied from the ADF to the Extraction Table.
   
   Input:
   VimosTable *adf  
   Pointer to the ADF that has the slit locations in the mask
   
   VimosIfuTable *ifuTab
   IFU Table where to get the IFU slits information

   VimosTable *extTable
   Pointer to the Extraction Table to be computed. The optical distortion
   models have to be present in the Extraction Table

   Return Value (success):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   This is returned if there is no distortion model present in the
   Extraction Table. 
   
   Updates:
   13 Dec 99: Modified the handling of IFU Observations, because the
              slits information has to be collected from the IFU Table,
	      instead of the ADF.  (MS)
   24 Nov 98: Added Curvature Model and Invers Dispersion Relation (TAO)
   18 Nov 98: Included copying of few descriptors from ADF, so checks on
              Quadrant etc are removed. (TAO)
   12 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- 
*/
VimosBool computeExtractionTable(VimosTable *adf, VimosIfuTable *ifuTab,
				 VimosExtractionTable *extTable);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   int numRowsInExtSlit(VimosExtractionSlit *slit)

   Description:

   Return number of spectra rows in a Extraction Slit. Since an Extraction
   Slit is a linked list of Extraction Slits, more than one slit can be
   contained in the list. The sum of all rows in all slits is returned.
   
   Input:
   VimosExtractionSlit *slit
   Pointer to the Extraction Slit for which the number of rows has to be counted.
   
   Return Value:
   Number of objects in the Extraction Slit
   
   Updates:
   09 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- 
*/
int numRowsInExtSlits(VimosExtractionSlit *slit);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   int numSlitsInExtTable(VimosExtractionTable *exTab)

   Description:

   Return number of slits in an Extraction Table.
   
   Input:
   VimosExtractionTable *exTab
   Pointer to the Extraction Table to count
   
   Return Value:
   Number of objects in the Extraction Table
   
   Updates:
   10 May 99: Created (TAO)

-------------------------------------------------------------------------------- 
*/
int numSlitsInExtTable(VimosExtractionTable *exTab);

VimosExtractionSlit *slitClosestToCenter(VimosExtractionTable *);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   VimosBool  slitLongOrShort(VimosExtractionSlit *slit, float tolerance)
 
   Description:
 
   Determine whether the slit is to be considered long or short for the sky
   subtraction. For long slits, the sky is subtracted from the 2D extracted
   spectra, for short spectra it is done from the spectra in Frame format
   (i.e. the un-extracted spectra).
   
   Input:
   VimosExtractionSlit *slit 
   Pointer to the Extraction Slit to check

   float tolerance 
   Tolerance for drift of slit on the CCD.  If the position of the slit on the
   CCD is straight within tolerance pixels, the slit is classified as short,
   otherwise as long.
   
   Return Value: VM_TRUE if slit is long, VM_FALSE if short
   
   Updates:
   22 Apr 99: Created (TAO)

-------------------------------------------------------------------------------- 
*/
VimosBool slitLongOrShort(VimosExtractionSlit *slit, float tolerance);

/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool readFitsExtractionTable(VimosExtractionTable *extTable, fitsfile
                                     *fptr)

   Description:
   Read a FITS Table into a Extraction Table. This requires the FITS file to 
   be already open, and does not close it at the end.
   
   Input:
   VimosExtractionTable *extTable  
   Pointer to the Extraction Table to copy the FITS table into.

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool (true if read is succesful, false otherwise)
   
   Updates:
   15 Dec 99: Created (BG)
   02 Feb 00: Changed the handling of the table, from a standalone FITS file,
              to an extension of the primary FITS file (MS)

--------------------------------------------------------------------------------
*/
VimosBool readFitsExtractionTable(VimosExtractionTable *extTable, fitsfile
				  *fptr);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool writeFitsExtractionTable(VimosExtractionTable *extTable, fitsfile
                                      *fptr)

   Description:
   Write an Extraction Table into a FITS Table. This requires the FITS file
   to be already open, and does not close it at the end. If an Extraction
   Table extension is already present in the file, it is first removed, and
   then the new one is written into a new extension.
   
   Input:
   VimosExtractionTable *extTable  
   Pointer to the Extraction Table to write to a FITS Table

   fitsfile *fptr
   Pointer to the FITS file where the table is to be found

   Return Value:
   VimosBool (true if write is succesful, false otherwise)
   
   Updates:
   15 Dec 99: Created (BG)
   02 Feb 00: Changed the handling of the table, from a standalone FITS file,
              to an extension of the primary FITS file (MS)

--------------------------------------------------------------------------------
*/
VimosBool writeFitsExtractionTable(VimosExtractionTable *extTable, fitsfile
				   *fptr);

PIL_END_DECLS

#endif /* VM_EXTRACTIONTABLE_H */
