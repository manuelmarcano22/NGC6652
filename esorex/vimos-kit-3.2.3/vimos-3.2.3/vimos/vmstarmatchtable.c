/* $Id: vmstarmatchtable.c,v 1.4 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Revision: 1.4 $
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
#include "vmstarmatchtable.h"


/*
 * Column definitions for a star match table
 */

static const char *star_match_table_columns[] = {
  "NUMBER",
  "ID",
  "X_IMAGE",
  "Y_IMAGE",
  "X_WORLD",
  "Y_WORLD",
  "MAG",
  "RA",
  "DEC",
  "MAG_U",
  "MAG_B",
  "MAG_V",
  "MAG_R",
  "MAG_I",
  "MAG_z"
};

static size_t nc = sizeof(star_match_table_columns) / sizeof(const char *);

/**
 * @name vmstarmatchtable Star Match Table
 *
 * @doc
 *   The module provides constructors for the creation of star match table
 *   objects. These constructors are specialized versions derived from
 *   the generic table constructor. No specialized version of a destructor
 *   is provided, since the generic table destructor is sufficient.
 */
 
/**@{*/

/**
 * @memo
 *   Create an empty star match table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a #NULL# pointer is returned.
 *
 * @doc
 *   The function uses the generic table constructor to create a generic
 *   table object. The table name is set to #VM_MATCH# and this value is
 *   added to the list of descriptors with the descriptor name referenced
 *   by the alias 'Table'.
 *
 * @author P. Sartoretti, R. Palsa
 */

VimosTable *newStarMatchTableEmpty()
{

  VimosTable *table = newTable();

  if (table) {

    /*
     * Set the table tag
     */

    strcpy(table->name, VM_MATCH);
    table->descs = newStringDescriptor(pilTrnGetKeyword("Table"), VM_MATCH,
				       "Type of table");
  }

  return table;

}


/**
 * @memo
 *   Create a star match table object.
 *
 * @return The function returns a pointer to the created table object if
 *   no error occured, otherwise a #NULL# pointer is returned.
 *
 * @doc
 *   The function creates an empty star match table using
 *   \textbf{newStarMatchTableEmpty()} and adds the columns
 *   \begin{itemize}
 *     \item #NUMBER#,
 *     \item #ID#,
 *     \item #X_IMAGE#,
 *     \item #Y_IMAGE#,
 *     \item #X_WORLD#,
 *     \item #Y_WORLD#,
 *     \item #MAG#,
 *     \item #RA#,
 *     \item #DEC#,
 *     \item #MAG_U#,
 *     \item #MAG_B#,
 *     \item #MAG_V#,
 *     \item #MAG_R#,
 *     \item #MAG_I# and
 *     \item #MAG_z#,
 *   \end{itemize}
 *   which define a star match table, to this object. The columns are created
 *   with \textbf{numRows} rows.
 *
 *
 * @see newStarMatchTableEmpty
 *
 * @author R. Palsa
 */

VimosTable *newStarMatchTable(size_t numRows)
{
 
  register size_t i;
 
  VimosTable *table = newStarMatchTableEmpty();
 
  if (table) {
 
    /*
     * Create the star match table columns with the appropriate size.
     */
 
    if ((tblAppendColumn(table, newIntColumn(numRows, 
		       star_match_table_columns[0]))) == EXIT_FAILURE) {
      deleteTable(table);
      return 0;
    }
    
    
    if ((tblAppendColumn(table, newStringColumn(numRows, 
		     star_match_table_columns[1]))) == EXIT_FAILURE) {
      deleteTable(table);
      return 0;
    }
    

    for (i = 2; i < nc; i++) {
      if ((tblAppendColumn(table, newDoubleColumn(numRows,  
	       star_match_table_columns[i]))) == EXIT_FAILURE) {
        deleteTable(table);
        return 0;
      }
    } 
  }
  
  return table;
 
}


/** 
 * @memo 
 *  Read Star Match Table from fitsfile (the file has to be open)
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param stmTable        photometric Table
 * @param fptr             pointer to fitsfile 
 *
 * @doc 
 *   Read Star Match Table from fitsfile (the file has to be already open)
 *
 * @author P. Sartoretti
 */

VimosBool readFitsStarMatchTable(VimosTable *stmTable, fitsfile *fptr)
{

 int status = 0;
 char     modName[] = "readFitsStarMatchTable";
 

  /* validate input */
  if (stmTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return (VM_FALSE);
  }

  if (fptr == NULL) {
    cpl_msg_error(modName,"NULL pointer to fitsfile");
    return (VM_FALSE);
  }
  
  if ( strcmp(stmTable->name, VM_MATCH) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }


/* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, VM_MATCH, 0, &status)) {
    cpl_msg_error(modName,"The function fits_movnam_hdu has returned an  "
		"error (code %d)", status);
    return(VM_FALSE);
  }

  stmTable->fptr = fptr;

   
 /* read Table */

  if(!readFitsTable(stmTable, stmTable->fptr)) {
    cpl_msg_info(modName, "Error in reading the FITS file");
    return VM_FALSE;
  }

  if (stmTable->numColumns == 0) {
    cpl_msg_warning(modName, "Empty input star match table");
    return VM_TRUE;
  }

  if(!checkStarMatchTable(stmTable)) {
    cpl_msg_error(modName, "check on table failed: invalid table");
    return VM_FALSE;
  }
  return VM_TRUE;
}

/** 
 * @memo 
 *  Check the Star Match Table 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param stmTable        star Table 
 *
 * @doc 
 *   Check if all necessary descriptors and columns of the star match Table 
 *   are present
 *
 * @author P. Sartoretti
 */
VimosBool checkStarMatchTable(VimosTable *stmTable) 
{
  char     modName[]  = "checkStarMatchTable";
  int      i, j ;
  
  if(stmTable == NULL) {
    cpl_msg_error(modName,"Null inputTable");
    return VM_FALSE;
  }
  if ( strcmp(stmTable->name, VM_MATCH) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }
 

  /* if(!findDescInTab(stmTable, pilTrnGetKeyword("AirMass"))) {
   cpl_msg_error(modName, "Descriptor AirMass not found");
   return VM_FALSE;
   } */
 
  if(!findDescInTab(stmTable, pilTrnGetKeyword("MagZero"))) {
    cpl_msg_error(modName, "Descriptor MagZero not found");
    return VM_FALSE;
  }
  /*BG: TabPixelScale (which is CDELT) check commented out: we have
    normally CD matrix, NOT CDELT */
  /*  if (!findDescInTab(stmTable, pilTrnGetKeyword("TabPixelScale"))) {
    cpl_msg_error(modName, "Descriptor PixelScale not found");
    return VM_FALSE;
    }*/

  /* BG wcs always written with its own appropriate k/w (i.e. CRPIX,
     and  NOT ESO PRO CRPIX, etc etc */
  for(i=1; i<=2; i++) {
    if(!findDescInTab(stmTable, pilTrnGetKeyword("Crpix",i))) {
      cpl_msg_error(modName, "Descriptor CRPIX not found");
      return VM_FALSE;
    }
    if(!findDescInTab(stmTable, pilTrnGetKeyword("Crval",i))) {
     cpl_msg_error(modName, "Descriptor CRVAL not found");
     return VM_FALSE;
   }
  }
  
 
  if(findDescInTab(stmTable, pilTrnGetKeyword("CD",1,1))) {
    for(i=1; i<=2; i++) {  
      for(j=1; j<=2; j++) {
	if(!(findDescInTab(stmTable, pilTrnGetKeyword("CD",i,j)))) {
	  cpl_msg_error(modName, "CD Matrix Incomplete");
	  return VM_FALSE;
	}
      }
    }
  }
  
  else {
    cpl_msg_error(modName, "Translation matrix descriptors not found");
    return VM_FALSE;
  }
  


  if(!(findDescInTab(stmTable, pilTrnGetKeyword("Equinox")))) {
    cpl_msg_error(modName, "Descriptor EQUINOX not found");
    return VM_FALSE;
  }
								 
  /*if(findColInTab(stmTable, "NUMBER") == NULL) { 
    cpl_msg_error(modName, "Column NUMBER ot found");
    return VM_FALSE;
  }
  */
  if(findColInTab(stmTable, "ID") == NULL) { 
    cpl_msg_error(modName, "Column ID ot found");
    return VM_FALSE;
  }
  if(findColInTab(stmTable, "MAG") == NULL) { 
    cpl_msg_error(modName, "Column MAG ot found");
    return VM_FALSE;
  }
  if(findColInTab(stmTable, "X_IMAGE") == NULL) { 
    cpl_msg_error(modName, "Column X_IMAGE ot found");
    return VM_FALSE;
  }
  if(findColInTab(stmTable, "Y_IMAGE") == NULL) { 
    cpl_msg_error(modName, "Column Y_IMAGE ot found");
    return VM_FALSE;
  }
  if(findColInTab(stmTable, "X_WORLD") == NULL) { 
    cpl_msg_error(modName, "Column X_WORLD ot found");
    return VM_FALSE;
  }
  if(findColInTab(stmTable, "Y_WORLD") == NULL) { 
    cpl_msg_error(modName, "Column Y_WORLD ot found");
    return VM_FALSE;
  }  
  if(findColInTab(stmTable, "RA") == NULL) { 
    cpl_msg_error(modName, "Column RA ot found");
    return VM_FALSE;
  }
  if(findColInTab(stmTable, "DEC") == NULL) { 
    cpl_msg_error(modName, "Column DEC ot found");
    return VM_FALSE;
  }

  return VM_TRUE;
}

/** 
 * @memo 
 *  Write Star Match Table to fitsfile 
 *
 * @return VM_TRUE or VM_FALSE.
 * 
 * @param filename        name of the fitsfile 
 * @param stmTable        star match Table
 *
 *
 * @doc 
 *   Write Star Match Table to fitsfile 
 *
 * @author P. Sartoretti
 */

VimosBool writeFitsStarMatchTable(char *filename, VimosTable *stmTable)
{
  char         modName[] = "writeFitsStarMatchTable";
  
  if(stmTable == NULL) {
    cpl_msg_error(modName, "Null input Table");
    return VM_FALSE;
  }
  
  if (strcmp(stmTable->name, VM_MATCH) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }

  if ((stmTable->cols) && stmTable->cols->len > 0) {
    /* check if descriptors and columns are there */
    if(!checkStarMatchTable(stmTable)) {
      cpl_msg_error(modName, "check on table failed: incomplete table");
      return VM_FALSE;
    }
  }
  if(!createFitsTable(filename, stmTable, VM_MATCH)) {
    cpl_msg_error(modName, "Error in writing fits table");
    return VM_FALSE;
  }
  return VM_TRUE;
}
/**@}*/
