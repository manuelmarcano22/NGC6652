/* $Id: vmextincttable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include "vmextincttable.h"


/*
 * Column definitions for a spectral atmospheric extinction table.
 */

static const char *extinct_table_columns[] = {
  "WAVE",                /* Angstrom                                    */
  "EXTINCTION"           /* Atmospheric extinction in mag/airmass       */
};

static size_t nc = sizeof(extinct_table_columns) / sizeof(const char *);

/**
 * @name vmextinctttable Spectral Atmospheric Extinction Table
 *
 * @doc
 *   The module provides constructors for the creation of Spectral Atmospheric 
 *   Extinction table objects. These constructors are specialized versions 
 *   derived from the generic table constructor. No specialized version of a
 *   destructor is provided, since the generic table destructor is sufficient.
 */
 
/**@{*/

/**
 * @memo
 *   Create an empty spectral atmospheric extinction table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a #NULL# pointer is returned.
 *
 * @doc
 *   The function uses the generic table constructor to create a generic
 *   table object. The table name is set to #VM_ATMEXT# and this value is
 *   added to the list of descriptors with the descriptor name referenced
 *   by the alias 'Table'.
 *
 * @author R. Palsa, C. Izzo
 */

VimosTable *newExtinctTableEmpty()
{

  VimosTable *table = newTable();

  if (table) {

    /*
     * Set the table tag
     */

    strcpy(table->name, VM_ATMEXT);
    table->descs = newStringDescriptor(pilTrnGetKeyword("Table"), VM_ATMEXT,
				       "Type of table");
  }

  return table;

}


/**
 * @memo
 *   Create a spectral atmospheric extinction table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a #NULL# pointer is returned.
 *
 * @doc
 *   The function creates an empty spectral atmospheric extinction table using
 *   \textbf{newExtinctTableEmpty()} and adds the columns
 *   \begin{itemize}
 *     \item #WAVE#,
 *     \item #EXTINCTION#,
 *   \end{itemize}
 *   which define a spectral atmospheric extinction table, to this object. 
 *   The columns are created with \textbf{numRows} rows.
 *
 *
 * @see newExtinctTableEmpty
 *
 * @author R. Palsa
 */

VimosTable *newExtinctTable(size_t numRows)
{
 
  register size_t i;
 
  VimosTable *table = newExtinctTableEmpty();
 
  if (table) {
 
    /*
     * Create the spectral atmospheric extinction table columns.
     */
 
    for (i = 0; i < nc; i++) {
      if ((tblAppendColumn(table, newDoubleColumn(numRows, 
		           extinct_table_columns[i]))) == EXIT_FAILURE) {
        deleteTable(table);
        return 0;
      }
    }
    
  }
  
  return table;
 
}


/** 
 * @memo 
 *  Read spectral atmospheric extinction Table from FITS file
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param extTable        Spectral Atmospheric Extinction Table
 * @param fptr            Pointer to FITS file 
 *
 * @doc 
 *   Read spectral atmospheric extinction Table from FITS file 
 *
 * @author R. Palsa, C. Izzo
 */

VimosBool readFitsExtinctTable(VimosTable *extTable, fitsfile *fptr)
{

  int status = 0;
  char     modName[] = "readFitsExtinctTable";
 

  /* 
   * Validate input 
   */

  if (extTable == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return VM_FALSE;
  }

  if (fptr == NULL) {
    cpl_msg_error(modName, "NULL pointer to FITS file");
    return VM_FALSE;
  }
  
  if (strcmp(extTable->name, VM_ATMEXT) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }


  /* 
   * Open Table 
   */

  if (fits_movnam_hdu(fptr, BINARY_TBL, VM_ATMEXT, 0, &status)) {
    cpl_msg_error(modName, "The function fits_movnam_hdu has returned an "
		  "error (code %d)", status);
    return(VM_FALSE);
  }

  extTable->fptr = fptr;

   
  /* 
   * Read Table 
   */

  if (!readFitsTable(extTable, extTable->fptr)) {
    cpl_msg_info(modName, "Error in reading the FITS file");
    return VM_FALSE;
  }

  if (!checkExtinctTable(extTable)) {
    cpl_msg_error(modName, "Invalid spectral atmospheric extinction table");
    return VM_FALSE;
  }

  return VM_TRUE;
}


/** 
 * @memo 
 *  Check a spectral atmospheric extinction Table 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param extTable    Spectral atmospheric extinction Table
 *
 * @doc 
 *   Check if all necessary descriptors and columns of the spectral 
 *   atmospheric extinction Table are present.
 *
 * @author P. Sartoretti
 */

VimosBool checkExtinctTable(VimosTable *extTable) 
{

  char     modName[]  = "checkExtinctTable";
  size_t   i;

  
  if (extTable == NULL) {
    cpl_msg_error(modName, "Null input table");
    return VM_FALSE;
  }

  if (strcmp(extTable->name, VM_ATMEXT)) {
    cpl_msg_error(modName, "Invalid input table");
    return (VM_FALSE);
  }
 
  for (i = 0; i < nc; i++) {
    if (findColInTab(extTable, extinct_table_columns[i]) == NULL) { 
      cpl_msg_error(modName, "Column %s not found", extinct_table_columns[i]);
      return VM_FALSE;
    }
  }

  return VM_TRUE;

}


/** 
 * @memo 
 *  Write spectral atmospheric extinction Table to FITS file 
 *
 * @return VM_TRUE or VM_FALSE.
 * 
 * @param filename   Name of the FITS file 
 * @param extTable   Spectral atmospheric extinction Table
 *
 * @doc 
 *   Write spectral atmospheric extinction Table to FITS file
 *
 * @author P. Sartoretti
 */

VimosBool writeFitsExtinctTable(char *filename, VimosTable *extTable)
{

  char modName[] = "writeFitsExtinctTable";

  
  if (extTable == NULL) {
    cpl_msg_error(modName, "Null input Table");
    return VM_FALSE;
  }
  
  if (strcmp(extTable->name, VM_ATMEXT)) {
    cpl_msg_error(modName, "Invalid input table");
    return VM_FALSE;
  }

  /* 
   * Check if descriptors and columns are there 
   */

  if (!checkExtinctTable(extTable)) {
    cpl_msg_error(modName, "Check on table failed: incomplete table");
    return VM_FALSE;
  }

  if (!createFitsTable(filename, extTable, VM_ATMEXT)) {
    cpl_msg_error(modName, "Error in writing FITS table");
    return VM_FALSE;
  }

  return VM_TRUE;

}

/**@}*/
