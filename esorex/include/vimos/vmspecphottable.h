/* $Id: vmspecphottable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
 
#ifndef VM_SPECPHOTTABLE_H
#define VM_SPECPHOTTABLE_H

#include <pilmacros.h>

#include <fitsio.h>

#include <vmtable.h>
 

PIL_BEGIN_DECLS

#define VM_SPHOT "SPH"

VimosTable *newSpecPhotTableEmpty(void);
VimosTable *newSpecPhotTable(size_t);

VimosBool readFitsSpecPhotTable(VimosTable *, fitsfile *);
VimosBool specPhotTableHeader(VimosTable *, VimosDescriptor *);
VimosBool writeFitsSpecPhotTable(char *, VimosTable *);

VimosBool checkSpecPhotTable(VimosTable *);

PIL_END_DECLS
 
#endif /* VM_SPECPHOTTABLE_H */
