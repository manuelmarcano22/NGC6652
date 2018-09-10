/* $Id: vmstartable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include "vmstartable.h"


/*
 * Column definitions for a star table
 */

static const char *star_table_columns[] = {
  "NUMBER",
  "X_IMAGE",
  "Y_IMAGE",
  "X_WORLD",
  "Y_WORLD",
  "MAG"
};

static size_t nc = sizeof(star_table_columns) / sizeof(const char *);


/**
 * @name vmstartable Star Table
 *
 * @doc
 *   The module provides constructors for the creation of star table
 *   objects. These constructors are specialized versions derived from
 *   the generic table constructor. No specialized version of a destructor
 *   is provided, since the generic table destructor is sufficient.
 */

/**@{*/

/**
 * @memo
 *   Create an empty star table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a NULL pointer is returned.
 *
 * @doc
 *   The function uses the generic table constructor to create a generic
 *   table object. The table name is set to #VM_STAR# and this value is
 *   added to the list of descriptors with the descriptor name referenced
 *   by the alias 'Table'.
 *
 * @author P. Sartoretti, R. Palsa
 */

VimosTable *newStarTableEmpty(void)
{

  VimosTable *table = newTable();

  if (table) {

    /*
     * Set the table tag
     */

    strcpy(table->name, VM_STAR);
    table->descs = newStringDescriptor(pilTrnGetKeyword("Table"), VM_STAR,
					"Type of table");
  }

  return table;

}


/**
 * @memo
 *   Create a star table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a NULL pointer is returned.
 *
 * @doc
 *   The function creates an empty star table using 
 *   \textbf{newStarTableEmpty()} and adds the columns
 *   \begin{itemize}
 *     \item #NUMBER#,
 *     \item #X_IMAGE#,
 *     \item #Y_IMAGE#,
 *     \item #X_WORLD#,
 *     \item #Y_WORLD# and
 *     \item #MAG#,
 *   \end{itemize}
 *   which define a star table, to this object. The columns are created
 *   with \textbf{numRows} rows.
 *
 * @see newStarTableEmpty
 *
 * @author R. Palsa
 */

VimosTable *newStarTable(size_t numRows)
{
  
  register size_t i;

  VimosTable *table = newStarTableEmpty();


  if (table) {

    /*
     * Create the star table columns with the appropriate size.
     */


    if ((tblAppendColumn(table, newIntColumn(numRows, 
			 star_table_columns[0]))) == EXIT_FAILURE) {
      deleteTable(table);
      return 0;
    }

    for (i = 1; i < nc; i++) {
      if ((tblAppendColumn(table, newDoubleColumn(numRows, 
			 star_table_columns[i]))) == EXIT_FAILURE) {
	deleteTable(table);
	return 0;
      }
    }

    table->numColumns = nc;
  }

  return table;

}

/** 
 * @memo 
 *  Read Star Table from fitsfile (the file has to be open)
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param starTable        photometric Table
 * @param fptr             pointer to fitsfile 
 *
 * @doc 
 *   Read Star Table from fitsfile (the file has to be already open)
 *
 * @author P. Sartoretti
 */


VimosBool readFitsStarTable(VimosTable *starTable, fitsfile *fptr)
{

 int status = 0;
 char     modName[] = "readFitsStarTable";


  /* validate input */
  if (starTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return (VM_FALSE);
  }

  if (fptr == NULL) {
    cpl_msg_error(modName,"NULL pointer to fitsfile");
    return (VM_FALSE);
  }
  
  if ( strcmp(starTable->name, VM_STAR) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }


/* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, VM_STAR, 0, &status)) {
    cpl_msg_error(modName,"The function fits_movnam_hdu has returned an  "
		"error (code %d)", status);
    return(VM_FALSE);
  }

  starTable->fptr = fptr;

   

 /* read Table */

  if (!readFitsTable(starTable, starTable->fptr)) {
    cpl_msg_info(modName, "Error in reading the FITS file");
    return VM_FALSE;
  }
  if (!checkStarTable(starTable)) {
    cpl_msg_info(modName, "Star Table is not complete");
    return VM_FALSE;
  }
       
  return VM_TRUE;
}
/** 
 * @memo 
 *  Check the Star Table 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param starTable        star Table to be checked 
 *
 * @doc 
 *   Check if all necessary descriptors and columns of the star Table 
 *   are present
 *
 * @author P. Sartoretti
 */
VimosBool checkStarTable(VimosTable *starTable)
{
  char     modName[] = "checkStarTable";


  if (starTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return (VM_FALSE);
  }

  if ( strcmp(starTable->name, VM_STAR) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }
  /* check if all the necessary descriptors and columns are there */
 
  if (!findDescInTab(starTable, pilTrnGetKeyword("AirMass"))) {
    cpl_msg_error(modName, "Descriptor AirMass not found");
    return VM_FALSE;
  }
  if (!findDescInTab(starTable, pilTrnGetKeyword("MagZero"))) {
    cpl_msg_error(modName, "Descriptor MagZero not found");
    return VM_FALSE;
  }
  
  if (findColInTab(starTable, "NUMBER") == NULL) { 
    cpl_msg_error(modName, "Column NUMBER ot found");
    return VM_FALSE;
  }
  if (findColInTab(starTable, "MAG") == NULL) { 
    cpl_msg_error(modName, "Column MAG ot found");
    return VM_FALSE;
  }
  if (findColInTab(starTable, "X_IMAGE") == NULL) { 
    cpl_msg_error(modName, "Column X_IMAGE ot found");
    return VM_FALSE;
  }
  if (findColInTab(starTable, "Y_IMAGE") == NULL) { 
    cpl_msg_error(modName, "Column Y_IMAGE ot found");
    return VM_FALSE;
  }
  if (findColInTab(starTable, "X_WORLD") == NULL) { 
    cpl_msg_error(modName, "Column X_WORLD ot found");
    return VM_FALSE;
  }
  if (findColInTab(starTable, "Y_WORLD") == NULL) { 
    cpl_msg_error(modName, "Column Y_WORLD ot found");
    return VM_FALSE;
  }
  return VM_TRUE;
}
/** 
 * @memo 
 *  Write Star Table to fitsfile 
 *
 * @return VM_TRUE or VM_FALSE.
 * 
 * @param filename        name of the fitsfile 
 * @param starTable        star Table
 *
 *
 * @doc 
 *   Write Star Table to fitsfile 
 *
 * @author P. Sartoretti
 */

VimosBool writeFitsStarTable(char *filename, VimosTable *starTable) 
{
 char modName[]  = "writeFitsStarTable";

 /* validate input */
  if (starTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return(VM_FALSE);
  }
  if ( strcmp(starTable->name, VM_STAR) ) {
    cpl_msg_error(modName,"Invalid input table");
    return(VM_FALSE);
  }
  if (!checkStarTable(starTable)) {
    cpl_msg_info(modName, "Star Table is not complete");
    return VM_FALSE;
  } 
   if (!createFitsTable(filename, starTable, VM_STAR)) {
    cpl_msg_error(modName, "Error in writing fits table");
    return VM_FALSE;
  }
  return VM_TRUE;
}
/**@}*/
