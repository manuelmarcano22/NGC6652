/* $Id: vmstdfluxtable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include "vmstdfluxtable.h"


/*
 * Column definitions for a standard star flux table.
 */

static const char *std_flux_table_columns[] = {
  "WAVE",                /* Angstrom                              */
  "FLUX",                /* Catalog flux in erg/cm/cm/s/Angstrom  */
  "BIN"                  /* Extracted spectrum, in Angstrom       */
};

static size_t nc = sizeof(std_flux_table_columns) / sizeof(const char *);

/**
 * @name vmstdfluxtable Standard Flux Table
 *
 * @doc
 *   The module provides constructors for the creation of standard flux 
 *   table objects. These constructors are specialized versions derived from
 *   the generic table constructor. No specialized version of a destructor
 *   is provided, since the generic table destructor is sufficient.
 */
 
/**@{*/

/**
 * @memo
 *   Create an empty standard star flux table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a #NULL# pointer is returned.
 *
 * @doc
 *   The function uses the generic table constructor to create a generic
 *   table object. The table name is set to #VM_SFLUX# and this value is
 *   added to the list of descriptors with the descriptor name referenced
 *   by the alias 'Table'.
 *
 * @author R. Palsa, C. Izzo
 */

VimosTable *newStdFluxTableEmpty()
{

  VimosTable *table = newTable();

  if (table) {

    /*
     * Set the table tag
     */

    strcpy(table->name, VM_SFLUX);
    table->descs = newStringDescriptor(pilTrnGetKeyword("Table"), VM_SFLUX,
				       "Type of table");
  }

  return table;

}


/**
 * @memo
 *   Create a standard star flux table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a #NULL# pointer is returned.
 *
 * @doc
 *   The function creates an empty standard star flux table using
 *   \textbf{newStdFluxTableEmpty()} and adds the columns
 *   \begin{itemize}
 *     \item #WAVE#,
 *     \item #FLUX#,
 *     \item #BIN#,
 *   \end{itemize}
 *   which define a standard star flux table, to this object. The columns 
 *   are created with \textbf{numRows} rows.
 *
 * @see newStdFluxTableEmpty
 *
 * @author R. Palsa
 */

VimosTable *newStdFluxTable(size_t numRows)
{
 
  register size_t i;
 
  VimosTable *table = newStdFluxTableEmpty();
 
  if (table) {
 
    /*
     * Create the standard star flux table columns.
     */
 
    for (i = 0; i < nc; i++) {
      if ((tblAppendColumn(table, newDoubleColumn(numRows, 
		           std_flux_table_columns[i]))) == EXIT_FAILURE) {
        deleteTable(table);
        return 0;
      }
    }
    
  }
  
  return table;
 
}


/** 
 * @memo 
 *  Read Standard Star Flux Table from FITS file (the file has to be open)
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param stdTable        Standard Star Flux Table
 * @param fptr            Pointer to FITS file 
 *
 * @doc 
 *   Read Standard Star Flux Table from FITS file (the file has to be 
 *   already open)
 *
 * @author R. Palsa, C. Izzo
 */

VimosBool readFitsStdFluxTable(VimosTable *stdTable, fitsfile *fptr)
{

  int status = 0;
  char     modName[] = "readFitsStdFluxTable";
 

  /* 
   * Validate input 
   */

  if (stdTable == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return VM_FALSE;
  }

  if (fptr == NULL) {
    cpl_msg_error(modName, "NULL pointer to FITS file");
    return VM_FALSE;
  }
  
  if (strcmp(stdTable->name, VM_SFLUX) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }


  /* 
   * Open Table 
   */

  if (fits_movnam_hdu(fptr, BINARY_TBL, VM_SFLUX, 0, &status)) {
    cpl_msg_error(modName, "The function fits_movnam_hdu has returned an "
		  "error (code %d)", status);
    return(VM_FALSE);
  }

  stdTable->fptr = fptr;

   
  /* 
   * Read Table 
   */

  if (!readFitsTable(stdTable, stdTable->fptr)) {
    cpl_msg_info(modName, "Error in reading the FITS file");
    return VM_FALSE;
  }

  if (!checkStdFluxTable(stdTable)) {
    cpl_msg_error(modName, "Invalid standard star flux table");
    return VM_FALSE;
  }

  return VM_TRUE;
}


/** 
 * @memo 
 *  Check a Standard Star Flux Table 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param stdTable    Standard Star Flux Table
 *
 * @doc 
 *   Check if all necessary descriptors and columns of the Standard Star Flux
 *   Table are present.
 *
 * @author P. Sartoretti
 */

VimosBool checkStdFluxTable(VimosTable *stdTable) 
{

  char     modName[]  = "checkStdFluxTable";
  size_t   i;

  
  if (stdTable == NULL) {
    cpl_msg_error(modName, "Null input table");
    return VM_FALSE;
  }

  if (strcmp(stdTable->name, VM_SFLUX)) {
    cpl_msg_error(modName, "Invalid input table");
    return (VM_FALSE);
  }
 
  for (i = 0; i < nc; i++) {
    if (findColInTab(stdTable, std_flux_table_columns[i]) == NULL) { 
      cpl_msg_error(modName, "Column %s not found", std_flux_table_columns[i]);
      return VM_FALSE;
    }
  }

  return VM_TRUE;

}


/** 
 * @memo 
 *  Write Standard Star Flux Table to FITS file 
 *
 * @return VM_TRUE or VM_FALSE.
 * 
 * @param filename   Name of the FITS file 
 * @param stdTable   Standard Star Flux Table
 *
 * @doc 
 *   Write Standard Star Flux Table to FITS file
 *
 * @author P. Sartoretti
 */

VimosBool writeFitsStdFluxTable(char *filename, VimosTable *stdTable)
{

  char modName[] = "writeFitsStdFluxTable";

  
  if (stdTable == NULL) {
    cpl_msg_error(modName, "Null input Table");
    return VM_FALSE;
  }
  
  if (strcmp(stdTable->name, VM_SFLUX)) {
    cpl_msg_error(modName, "Invalid input table");
    return VM_FALSE;
  }

  /* 
   * Check if descriptors and columns are there 
   */

  if (!checkStdFluxTable(stdTable)) {
    cpl_msg_error(modName, "Check on table failed: incomplete table");
    return VM_FALSE;
  }

  if (!createFitsTable(filename, stdTable, VM_SFLUX)) {
    cpl_msg_error(modName, "Error in writing FITS table");
    return VM_FALSE;
  }

  return VM_TRUE;

}

/**@}*/
