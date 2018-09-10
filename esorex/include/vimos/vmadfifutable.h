/* $Id: vmadfifutable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_ADFIFUTABLE_H
#define VM_ADFIFUTABLE_H

#include <pilmacros.h>

#include <vmtable.h>
#include <vmadf.h>
#include <vmifutable.h>


PIL_BEGIN_DECLS

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   VimosAdfSlitHolder *extractSlitsFromIFU(VimosTable *adf, 
                                 VimosIfuTable *ifuTab, VimosIfuMode ifuMode)
 
   Description:
   
   Input:
   VimosTable *adf
   Aperture Definition File, where to get info on the VIMOS Quadrant
 
   VimosIfuTable *ifuTab
   IFU Table where to get slits information
 
   VimosIfuMode ifuMode
   Distinguishes between low and high resolution IFU mode.
   
   Return Value (success):
 
   Return Value (error):
   NULL
   
   Updates:
   13 Dec 99: Created (MS)
 
--------------------------------------------------------------------------------
 
*/
VimosAdfSlitHolder *extractSlitsFromIFU(VimosTable *adf, 
			      VimosIfuTable *ifuTab, VimosIfuMode ifuMode);

PIL_END_DECLS
 
#endif /* VM_ADFIFUTABLE_H */
