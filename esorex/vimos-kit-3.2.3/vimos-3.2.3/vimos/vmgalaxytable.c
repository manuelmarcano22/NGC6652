/* $Id: vmgalaxytable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include "vmmath.h"
#include "vmgalaxytable.h"


/**
 * @name vmgalaxytable Galaxy Table
 *
 * @doc
 *   The module  provides the functions to treat galaxy tables.
 */

/**@{*/

/*
 *  Allocate a new GALAXY Table. This is a VimosTable with name VM_GAL (which
 * should have value "GAL").
 */

VimosTable *newGalaxyTable()
{
  VimosTable *newTab;

  /* allocate new VimosTable */
  newTab = newTable();

  /* if error, exit(-1) */
  if (newTab == NULL) {
    return NULL;
  }

  /* copy "STAR" into name of table */
  strcpy(newTab->name, VM_GAL);
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_GAL, "");

  /* return address of new Ccd Table */
  return(newTab);

}


/** 
 * @memo 
 *  Read galaxy Table from fitsfile (the file has to be open)
 *
 * @return VM_TRUE VM_FALSE. 
 *
 * @param galTable        galaxy Table
 * @param fptr            pointer to fitsfile 
 *
 * @doc 
 *   Read galaxy Table from fitsfile (the file has to be already open)
 *
 * @author P. Sartoretti
 */



VimosBool readFitsGalaxyTable(VimosTable *galTable, fitsfile *fptr)
{

 int status = 0;
 char     modName[] = "readFitsGalaxyTable";

  /* validate input */
  if (galTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return (VM_FALSE);
  }

  if (fptr == NULL) {
    cpl_msg_error(modName,"NULL pointer to fitsfile");
    return (VM_FALSE);
  }
  
  if ( strcmp(galTable->name, VM_GAL) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }
  
  
  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL,
		      (char *) pilTrnGetCategory("GalaxyTable"), 0, &status)) {
    cpl_msg_error(modName,"The function fits_movnam_hdu has returned an  "
		"error (code %d)", status);
    return(VM_FALSE);
  }
  
  galTable->fptr = fptr;
  
  /* read Table */
  
  if(!readFitsTable(galTable, galTable->fptr)) {
    cpl_msg_info(modName, "Error in reading the FITS file");
    return VM_FALSE;
  }

  /* check table */
  if(!checkGalaxyTable(galTable)) {
    cpl_msg_error(modName, "Incomplete table");
    return VM_FALSE;
  }
  
  return VM_TRUE;
}

/** 
 * @memo 
 *  Check the galaxy Table 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param galTable        galaxy Table 
 *
 * @doc 
 *   Check if all necessary descriptors and columns of the galaxy Table are 
 *   present
 *
 * @author P. Sartoretti
 */

VimosBool checkGalaxyTable(VimosTable *galTable)
{
  char modName[]  = "checkGalaxyTable";
  
  if (galTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return(VM_FALSE);
  }
  if ( strcmp(galTable->name, VM_GAL) ) {
    cpl_msg_error(modName,"Invalid input table");
    return(VM_FALSE);
  }

  if(findColInTab(galTable, "X_IMAGE") == NULL) { 
    cpl_msg_error(modName, "Column X_IMAGE ot found");
    return VM_FALSE;
  }
  if(findColInTab(galTable, "Y_IMAGE") == NULL) { 
    cpl_msg_error(modName, "Column Y_IMAGE ot found");
    return VM_FALSE;
  }
  if(findColInTab(galTable, "X_WORLD") == NULL) { 
    cpl_msg_error(modName, "Column X_WORLD ot found");
    return VM_FALSE;
  }
  if(findColInTab(galTable, "Y_WORLD") == NULL) { 
    cpl_msg_error(modName, "Column Y_WORLD ot found");
    return VM_FALSE;
  }
  return VM_TRUE;
}
 
/** 
 * @memo 
 *  Write galaxy Table to fitsfile 
 *
 * @return VM_TRUE or VM_FALSE.
 * 
 * @param filename        name of the fitsfile 
 * @param galTable        galaxy Table
 *
 *
 * @doc 
 *   Write galaxy Table to fitsfile 
 *
 * @author P. Sartoretti
 */

VimosBool writeFitsGalaxyTable(char *filename, VimosTable *galTable)
{
  char modName[]  = "writeFitsGalaxyTable";

 /* validate input */
  if (galTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return(VM_FALSE);
  }
  
  if(!checkGalaxyTable(galTable)) {
    cpl_msg_error(modName, "Invalid input Table");
    return VM_FALSE;
  }
  
  /* write Table */
  if(!createFitsTable(filename, galTable, pilTrnGetCategory("GalaxyTable"))) {
    cpl_msg_error(modName, "Error in writing fits table");
    return VM_FALSE;
  }
  return VM_TRUE;
}
/**@}*/
