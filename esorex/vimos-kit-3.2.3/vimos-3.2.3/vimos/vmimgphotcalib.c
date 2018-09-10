/* $Id: vmimgphotcalib.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmtable.h"
#include "vmimage.h"


/**
 * @brief vmimgphotcalib Imaging Photometric Calibration
 *
 * The module provides the function for applying a photometric
 * calibration to an image.
 */

/**@{*/

/**
 * @memo 
 *   Apply photometric calibration to an image.
 *
 * @return VM_TRUE in case of success, otherwise VM_FALSE
 *
 * @param ima       Input image
 * @param photTable Photometric table
 *
 * @doc
 *   Copy ALL Photometric Calibration Table in the image header 
 *
 * @author P. Sartoretti Modified B.Garilli
 */

int
VmImApplyPhot(VimosImage *ima, VimosTable *photTable) 
{
  char modName[] = "VmImApplyPhot";

  if ((ima == NULL) || (photTable == NULL)) {
    cpl_msg_error(modName, "Null input");
    return VM_FALSE;
  }

  /* NOTE by BG: Not ONLY the zero point is needed, but also the other
   parameters in phot. table (extincion, colour term, etc.) with their rms*/

 if (!(copyFromHeaderToHeader(photTable->descs, pilTrnGetKeyword("MagZero"), 
			     &(ima->descs), NULL))) {
   cpl_msg_error(modName, "Missing descriptor %s", pilTrnGetKeyword("MagZero"));
   return VM_FALSE;
 }
 copyFromHeaderToHeader(photTable->descs, "ESO PRO MAGZERO RMS",
                        &(ima->descs), NULL); 

 if (!(copyFromHeaderToHeader(photTable->descs,
                              pilTrnGetKeyword("Extinction"), 
                              &(ima->descs), NULL))) {
   cpl_msg_warning(modName, "Missing descriptor %s",
                 pilTrnGetKeyword("Extinction"));
 }
 copyFromHeaderToHeader(photTable->descs, "ESO PRO EXTINC RMS",
                        &(ima->descs), NULL); 

 if (!(copyFromHeaderToHeader(photTable->descs, pilTrnGetKeyword("Colour"), 
                              &(ima->descs), NULL))) {
   cpl_msg_warning(modName, "Missing descriptor %s", pilTrnGetKeyword("Colour"));
 }
 if (!(copyFromHeaderToHeader(photTable->descs, pilTrnGetKeyword("ColorTerm"), 
                              &(ima->descs), NULL))) {
   cpl_msg_warning(modName, "Missing descriptor %s",
                 pilTrnGetKeyword("ColorTerm"));
 }
 copyFromHeaderToHeader(photTable->descs, "ESO PRO COLTERM RMS",
                        &(ima->descs), NULL); 



 return VM_TRUE;
}
/**@}*/
