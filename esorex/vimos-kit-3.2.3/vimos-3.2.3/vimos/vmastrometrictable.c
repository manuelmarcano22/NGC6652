/* $Id: vmastrometrictable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>

#include <fitsio.h>

#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmtable.h"
#include "vmmath.h"
#include "vmastrometrictable.h"


/**
 * @name vmastrometrictable Astrometric Table
 *
 * The module astrometricTableFits includes the functions to treat 
 * fits astrometric tables.
 */

/**@{*/

/*
 *  Allocate a new AST Table. This is a VimosTable with name VM_AST (which
 * should have value "AST").
 */

VimosTable *newAstrometricTable()
{
  VimosTable *newTab;

  /* allocate new VimosTable */
  newTab = newTable();

  /* if error, exit(-1) */
  if (newTab == NULL) {
    return NULL;
  }

  /* copy "STAR" into name of table */
  strcpy(newTab->name, VM_AST);
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_AST, "");

  /* return address of new Astrometric Table */
  return(newTab);

}

/*
  Delete an Astrometric Table. This is just an esthetic wrapper for
  deleteTable(astroTable)
*/
void deleteAstrometricTable(VimosTable *astroTable)
{
  deleteTable(astroTable);
}


/** 
 * @memo 
 *  Read astrometric Table from fitsfile (the file has to be open)
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param astTable        astrometric Table
 * @param fptr            pointer to fitsfile 
 *
 * @doc 
 *   Read astrometric Table from fitsfile (the file has to be already open)
 *
 * @author P. Sartoretti
 */


VimosBool readFitsAstrometricTable(VimosTable *astTable, fitsfile *fptr)
{

 int status = 0;
 char     modName[] = "readFitsAstrometricTable";

  /* validate input */
  if (astTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return (VM_FALSE);
  }

  if (fptr == NULL) {
    cpl_msg_error(modName,"NULL pointer to fitsfile");
    return (VM_FALSE);
  }
  
  if ( strcmp(astTable->name, VM_AST) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }


/* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, VM_AST, 0, &status)) {
    cpl_msg_error(modName,"The function fits_movnam_hdu has returned an  "
		"error (code %d)", status);
    return(VM_FALSE);
  }

  astTable->fptr = fptr;

   
 /* read Table */

  if(!readFitsTable(astTable, astTable->fptr)) {
    cpl_msg_error(modName, "Error in reading the FITS file");
    return VM_FALSE;
  }
  
  /* check if all the descriptors and columns (with rigth name) are 
   present*/
  if(!checkAstrometricTable(astTable)) {
    cpl_msg_error(modName, "Astrometric table is incomplete");
    return VM_FALSE;
  }
  return VM_TRUE;
}


/** 
 * @memo 
 *  Write astrometric Table to fitsfile 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param astTable        astrometric Table
 * @param filename        name of the fitsfile 
 *
 * @doc 
 *   Write astrometric Table to fitsfile 
 *
 * @author P. Sartoretti
 */

VimosBool writeFitsAstrometricTable(char *filename, VimosTable *astTable)
{
  char modName[]  = "writeFitsAstrometricTable";

 /* 
  * validate input 
  */
  if (astTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return(VM_FALSE);
  }
  
  /*
   * check table 
   */
  if(!checkAstrometricTable(astTable)) {
    cpl_msg_error(modName, "Astrometric table is incomplete");
    return VM_FALSE;
  }
  if(!createFitsTable(filename, astTable, VM_AST)) {
    cpl_msg_error(modName, "Error in writing fits table");
    return VM_FALSE;
  }
  return VM_TRUE;
}


/** 
 * @memo 
 *  Check the astrometric Table 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param astTable        astrometric Table 
 *
 * @doc 
 *   Check if all keywords and columns of the astrometric Table are there
 *
 * @author P. Sartoretti
 */
  
VimosBool checkAstrometricTable(VimosTable *astTable)
{
  char modName[]   = "checkAstrometricTable";
  VimosColumn *column;
  int count = 0;

  if(astTable == NULL) {
    cpl_msg_info(modName, "Null Input Table");
    return VM_FALSE;
  }

  if (strcmp(astTable->name, VM_AST)) {
    cpl_msg_error(modName, "Invalid input table");
    return(VM_FALSE);
  }

  if ((column = findColInTab(astTable, "ID"))) { 
    if (column->colType != VM_STRING) {
      cpl_msg_error(modName, "Column ID has wrong type - should be VM_STRING");
      return VM_FALSE;
    }
  }
  else {
    cpl_msg_error(modName, "Column ID not found");
    return VM_FALSE;
  }

  if ((column = findColInTab(astTable, "RA"))) {
    if (column->colType != VM_DOUBLE) {
      cpl_msg_error(modName, "Column RA has wrong type - should be VM_DOUBLE");
      return VM_FALSE;
    }
  } 
  else { 
    cpl_msg_error(modName, "Column RA not found");
    return VM_FALSE;
  }
  
  if ((column = findColInTab(astTable, "DEC"))) { 
    if (column->colType != VM_DOUBLE) {
      cpl_msg_error(modName, "Column DEC has wrong type - should be VM_DOUBLE");
      return VM_FALSE; 
    } 
  }  
  else { 
    cpl_msg_error(modName, "Column DEC not found");
    return VM_FALSE;
  }
 
  if ((column = findColInTab(astTable, "MAG_U"))) { 
    if (column->colType != VM_DOUBLE) {
      cpl_msg_error(modName, "Column MAG_U wrong type - should be VM_DOUBLE");
      return VM_FALSE; 
    } 
  }  
  else { 
    cpl_msg_warning(modName, "Column MAG_U not found");
    count++;
  }

  if ((column = findColInTab(astTable, "MAG_B"))) { 
    if (column->colType != VM_DOUBLE) {
      cpl_msg_error(modName, "Column MAG_B wrong type - should be VM_DOUBLE");
      return VM_FALSE; 
    } 
  }  
  else { 
    cpl_msg_warning(modName, "Column MAG_B not found");
    count++;
  }

  if ((column = findColInTab(astTable, "MAG_V"))) { 
    if (column->colType != VM_DOUBLE) {
      cpl_msg_error(modName, "Column MAG_V wrong type - should be VM_DOUBLE");
      return VM_FALSE; 
    } 
  }  
  else { 
    cpl_msg_warning(modName, "Column MAG_V not found");
    count++;
  }

  if ((column = findColInTab(astTable, "MAG_R"))) { 
    if (column->colType != VM_DOUBLE) {
      cpl_msg_error(modName, "Column MAG_R wrong type - should be VM_DOUBLE");
      return VM_FALSE; 
    } 
  }  
  else { 
    cpl_msg_warning(modName, "Column MAG_R not found");
    count++;
  }

  if ((column = findColInTab(astTable, "MAG_I"))) { 
    if (column->colType != VM_DOUBLE) {
      cpl_msg_error(modName, "Column MAG_I wrong type - should be VM_DOUBLE");
      return VM_FALSE; 
    } 
  }  
  else { 
    cpl_msg_warning(modName, "Column MAG_I not found");
    count++;
  }

  if ((column = findColInTab(astTable, "MAG_z"))) { 
    if (column->colType != VM_DOUBLE) {
      cpl_msg_error(modName, "Column MAG_z wrong type - should be VM_DOUBLE");
      return VM_FALSE; 
    } 
  }  
  else { 
    cpl_msg_warning(modName, "Column MAG_z not found");
    count++;
  }

  if (count == 6) { 
    cpl_msg_error(modName, "No magnitude column found");
    return VM_FALSE;
  }

  return VM_TRUE;
}
/**@}*/
