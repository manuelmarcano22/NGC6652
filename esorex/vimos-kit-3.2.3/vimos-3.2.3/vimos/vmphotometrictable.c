/* $Id: vmphotometrictable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <piltranslator.h>

#include "vmtable.h"
#include "vmmath.h"
#include "vmphotometrictable.h"


/**
 * @name vmphotometrictable Photometric Table
 *
 * @doc
 *   The module provides the functions to treat photometric tables.
 */

/**@{*/

/*
 *  Allocate a new Photometric Table. This is a VimosTable with name 
 *  VM_IPC (which should have value "IPC").
 */

VimosTable *newPhotometricTable()
{
  VimosTable *newTab;

  /* allocate new VimosTable */
  newTab = newTable();

  if (newTab == NULL) {
    return NULL;
  }

  /* copy "IPC" into name of table */
  strcpy(newTab->name, VM_IPC);
  newTab->descs = newStringDescriptor("ESO PRO TABLE", VM_IPC, "");

  /* return address of new Photometric Table */
  return(newTab);

}

/** 
 * @memo 
 *  Read photometric Table from fitsfile (the file has to be open)
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE. 
 *
 * @param photTable        photometric Table
 * @param fptr             pointer to fitsfile 
 *
 * @doc 
 *   Read photometric Table from fitsfile (the file has to be already open)
 *
 * @author P. Sartoretti
 */


VimosBool readFitsPhotometricTable(VimosTable *photTable, fitsfile *fptr)
{
  int status = 0;
  char     modName[] = "readFitsPhotometricTable";


  /* validate input */
  if (photTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return (VM_FALSE);
  }

  if (fptr == NULL) {
    cpl_msg_error(modName,"NULL pointer to fitsfile");
    return (VM_FALSE);
  }

  if (strcmp(photTable->name, VM_IPC) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }

/* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, VM_IPC, 0, &status)) {
    cpl_msg_error(modName,"The function fits_movnam_hdu has returned an  "
                "error (code %d)", status);
    return(VM_FALSE);
  }

  photTable->fptr = fptr;


 /* read Table */

  if (!readFitsTable(photTable, photTable->fptr)) {
    cpl_msg_info(modName, "Error in reading the FITS file");
    return VM_FALSE;
  }
  if (!checkPhotometricTable(photTable)) {
    cpl_msg_info(modName, "Photometric Table is not complete");
    return VM_FALSE;
  }

  return VM_TRUE;
}
/** 
 * @memo 
 *  Check the photometric Table 
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE. 
 *
 * @param photTable        photometric Table 
 *
 * @doc 
 *   Check if all necessary descriptors and columns of the photometric Table 
 *   are present
 *
 * @author P. Sartoretti
 */

VimosBool checkPhotometricTable(VimosTable *photTable)
{
  char     modName[] = "checkPhotometricTable";


  if (photTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return (VM_FALSE);
  }

  if ( strcmp(photTable->name, VM_IPC) ) {
    cpl_msg_error(modName,"Invalid input table");
    return (VM_FALSE);
  }
  /* check if all the necessary descriptors are there */

  if (!findDescInTab(photTable, pilTrnGetKeyword("MagZero"))) {
    cpl_msg_error(modName, "Descriptor MagZero not found");
    return VM_FALSE;
  }
   if (!findDescInTab(photTable, pilTrnGetKeyword("Extinction"))) {
    cpl_msg_error(modName, "Descriptor Extinction not found");
    return VM_FALSE;
  }
 return VM_TRUE;
}


/** 
 * @memo 
 *  Write photometric Table to fitsfile 
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE.
 * 
 * @param filename        name of the fitsfile 
 * @param photTable       photometric Table
 *
 *
 * @doc 
 *   Write photometric Table to fitsfile 
 *
 * @author P. Sartoretti
 */
VimosBool writeFitsPhotometricTable(char *filename, VimosTable *photTable)
{
 char modName[]  = "writeFitsPhotometricTable";

 /* validate input */
  if (photTable == NULL) {
    cpl_msg_error(modName,"NULL input table");
    return(VM_FALSE);
  }
  if ( strcmp(photTable->name, VM_IPC) ) {
    cpl_msg_error(modName,"Invalid input table");
    return(VM_FALSE);
  }
  if (!checkPhotometricTable(photTable)) {
    cpl_msg_info(modName, "Photometric Table is not complete");
    return VM_FALSE;
  }
   if (!createFitsTable(filename, photTable, VM_IPC)) {
    cpl_msg_error(modName, "Error in writing fits table");
    return VM_FALSE;
  }
  return VM_TRUE;
}
/**@}*/
