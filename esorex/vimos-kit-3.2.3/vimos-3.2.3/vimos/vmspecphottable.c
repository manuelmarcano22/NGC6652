/* $Id: vmspecphottable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmtable.h"
#include "vmspecphottable.h"


/*
 * Column definitions for a spectrophotometric table.
 */

static const char *spec_phot_table_columns[] = {
  "WAVE",                /* Angstrom                                    */
  "STD_FLUX",            /* Catalog flux in erg/cm/cm/s/Angstrom        */
  "OBS_FLUX",            /* Extracted spectrum, in electrons/s/Angstrom */
  "RAW_EFFICIENCY",      /* Ratio between input and detected photons    */
  "EFFICIENCY",          /* Smoothed version of the above               */
  "RAW_RESPONSE",        /* Ratio between "STD_FLUX" and "OBS_FLUX"     */
  "RESPONSE"             /* Smoothed version of the above               */
};

static size_t nc = sizeof(spec_phot_table_columns) / sizeof(const char *);

/**
 * @name vmspecphottable Spectrophotometric Table
 *
 * @doc
 *   The module provides constructors for the creation of spectrophotometric 
 *   table objects. These constructors are specialized versions derived from
 *   the generic table constructor. No specialized version of a destructor
 *   is provided, since the generic table destructor is sufficient.
 */
 
/**@{*/

/**
 * @memo
 *   Create an empty spectrophotometric table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a #NULL# pointer is returned.
 *
 * @doc
 *   The function uses the generic table constructor to create a generic
 *   table object. The table name is set to #VM_SPHOT# and this value is
 *   added to the list of descriptors with the descriptor name referenced
 *   by the alias 'Table'.
 *
 * @author R. Palsa, C. Izzo
 */

VimosTable *newSpecPhotTableEmpty()
{

  VimosTable *table = newTable();

  if (table) {

    /*
     * Set the table tag
     */

    strcpy(table->name, VM_SPHOT);
    table->descs = newStringDescriptor(pilTrnGetKeyword("Table"), VM_SPHOT,
				       "Type of table");
  }

  return table;

}


/**
 * @memo
 *   Create a spectrophotometric table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a #NULL# pointer is returned.
 *
 * @doc
 *   The function creates an empty spectrophotometric table using
 *   \textbf{newSpecPhotTableEmpty()} and adds the columns
 *   \begin{itemize}
 *     \item #WAVE#,
 *     \item #STD_FLUX#,
 *     \item #OBS_FLUX#,
 *     \item #RAW_EFFICIENCY#,
 *     \item #EFFICIENCY#,
 *     \item #RAW_RESPONSE#,
 *     \item #RESPONSE#,
 *   \end{itemize}
 *   which define a spectrophotometric table, to this object. The columns 
 *   are created with \textbf{numRows} rows.
 *
 * @see newSpecPhotTableEmpty
 *
 * @author R. Palsa
 */

VimosTable *newSpecPhotTable(size_t numRows)
{
 
  register size_t i;
 
  VimosTable *table = newSpecPhotTableEmpty();
 
  if (table) {
 
    /*
     * Create the spectrophotomoetric table columns.
     */
 
    for (i = 0; i < nc; i++) {
      if ((tblAppendColumn(table, newDoubleColumn(numRows, 
		           spec_phot_table_columns[i]))) == EXIT_FAILURE) {
        deleteTable(table);
        return 0;
      }
    }
    
  }
  
  return table;
 
}


/** 
 * @memo 
 *  Read Spectrophotometric Table from FITS file (the file has to be open)
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param sphTable        Spectrophotometric Table
 * @param fptr            Pointer to FITS file 
 *
 * @doc 
 *   Read Spectrophotometric Table from FITS file (the file has to be 
 *   already open)
 *
 * @author R. Palsa, C. Izzo
 */

VimosBool readFitsSpecPhotTable(VimosTable *sphTable, fitsfile *fptr)
{

  int status = 0;
  char     modName[] = "readFitsSpecPhotTable";
 

  /* 
   * Validate input 
   */

  if (sphTable == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return VM_FALSE;
  }

  if (fptr == NULL) {
    cpl_msg_error(modName, "NULL pointer to FITS file");
    return VM_FALSE;
  }
  
  if (strcmp(sphTable->name, VM_SPHOT) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }


  /* 
   * Open Table 
   */

  if (fits_movnam_hdu(fptr, BINARY_TBL, VM_SPHOT, 0, &status)) {
    cpl_msg_error(modName, "The function fits_movnam_hdu has returned an "
		  "error (code %d)", status);
    return(VM_FALSE);
  }

  sphTable->fptr = fptr;

   
  /* 
   * Read Table 
   */

  if (!readFitsTable(sphTable, sphTable->fptr)) {
    cpl_msg_info(modName, "Error in reading the FITS file");
    return VM_FALSE;
  }

  if (!checkSpecPhotTable(sphTable)) {
    cpl_msg_error(modName, "Invalid spectrophotometric table");
    return VM_FALSE;
  }

  return VM_TRUE;
}


/** 
 * @memo 
 *  Check a Spectrophotometric Table 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param sphTable    Spectrophotometric Table
 *
 * @doc 
 *   Check if all necessary descriptors and columns of the Spectrophotometric
 *   Table are present.
 *
 * @author P. Sartoretti
 */

VimosBool checkSpecPhotTable(VimosTable *sphTable) 
{

  char     modName[]  = "checkSpecPhotTable";
  size_t   i;

  
  if (sphTable == NULL) {
    cpl_msg_error(modName, "Null input table");
    return VM_FALSE;
  }

  if (strcmp(sphTable->name, VM_SPHOT)) {
    cpl_msg_error(modName, "Invalid input table");
    return (VM_FALSE);
  }
 
  for (i = 0; i < nc; i++) {
    if (findColInTab(sphTable, spec_phot_table_columns[i]) == NULL) { 
      cpl_msg_error(modName, "Column %s not found", spec_phot_table_columns[i]);
      return VM_FALSE;
    }
  }

  return VM_TRUE;

}


VimosBool specPhotTableHeader(VimosTable *table, VimosDescriptor *descs)
{

  char modName[]  = "specPhotTableHeader";
  int  quadrant;


  if (table == NULL) {
    cpl_msg_error(modName, "Null input table");
    return VM_FALSE;
  }

  if (descs == NULL) {
    cpl_msg_error(modName, "Null input descriptors");
    return VM_FALSE;
  }

  if (strcmp(table->name, VM_SPHOT)) {
    cpl_msg_error(modName, "Invalid input table");
    return (VM_FALSE);
  }

  if (readIntDescriptor(descs, pilTrnGetKeyword("Quadrant"),
                        &quadrant, NULL) != VM_TRUE)
    return 0;

  if (vimosDscCopy(&(table->descs), descs, ".*-OBS$", 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs, pilTrnGetKeyword("Instrument"), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs, "^ESO OBS (DID|ID|PROG ID)", 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs, pilTrnGetKeyword("INS.DID"), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   pilTrnGetKeyword("InstrumentMode"), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   pilTrnGetKeyword("FilterId", quadrant), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   pilTrnGetKeyword("FilterName", quadrant), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs, pilTrnGetKeyword("DET.DID"), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   pilTrnGetKeyword("Adu2Electron", 1), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   pilTrnGetKeyword("ReadNoise", 1), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   pilTrnGetKeyword("Electron2Adu", 1), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   pilTrnGetKeyword("GrismId", quadrant), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   pilTrnGetKeyword("GrismName", quadrant), 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs,
                   "^ESO DET READ (CLOCK|SPEED|MODE)", 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs, "^ESO OCS (CON QUAD|DID)", 0))
      return 0;

  if (vimosDscCopy(&(table->descs), descs, pilTrnGetKeyword("Airmass"), 0))
      return 0;

  return VM_TRUE;

}


/** 
 * @memo 
 *  Write Spectrophotometric Table to FITS file 
 *
 * @return VM_TRUE or VM_FALSE.
 * 
 * @param filename   Name of the FITS file 
 * @param sphTable   Spectrophotometric Table
 *
 * @doc 
 *   Write Spectrophotometric Table to FITS file
 *
 * @author P. Sartoretti
 */

VimosBool writeFitsSpecPhotTable(char *filename, VimosTable *sphTable)
{

  char modName[] = "writeFitsSpecPhotTable";

  
  if (sphTable == NULL) {
    cpl_msg_error(modName, "Null input Table");
    return VM_FALSE;
  }
  
  if (strcmp(sphTable->name, VM_SPHOT)) {
    cpl_msg_error(modName, "Invalid input table");
    return VM_FALSE;
  }

  /* 
   * Check if descriptors and columns are there 
   */

  if (!checkSpecPhotTable(sphTable)) {
    cpl_msg_error(modName, "Check on table failed: incomplete table");
    return VM_FALSE;
  }

  if (!createFitsTable(filename, sphTable, VM_SPHOT)) {
    cpl_msg_error(modName, "Error in writing FITS table");
    return VM_FALSE;
  }

  return VM_TRUE;

}

/**@}*/
