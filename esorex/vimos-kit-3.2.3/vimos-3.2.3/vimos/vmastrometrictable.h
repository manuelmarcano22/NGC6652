/* $Id: vmastrometrictable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_ASTROMETRICTABLE_H
#define VM_ASTROMETRICTABLE_H

#include <fitsio.h>

#include <pilmacros.h>

#include <vmtable.h>


PIL_BEGIN_DECLS

#define VM_AST "AST"

VimosTable *newAstrometricTable(void);
void deleteAstrometricTable(VimosTable *);

VimosBool readFitsAstrometricTable(VimosTable *astTable, fitsfile *fptr);
VimosBool writeFitsAstrometricTable(char *filename, VimosTable *astTable);
VimosBool checkAstrometricTable(VimosTable *astTable);

PIL_END_DECLS

#endif /* VM_ASTROMETRICTABLE_H */
