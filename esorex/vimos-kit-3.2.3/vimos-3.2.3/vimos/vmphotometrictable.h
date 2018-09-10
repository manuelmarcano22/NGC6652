/* $Id: vmphotometrictable.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_PHOTOMETRICTABLE_H
#define VM_PHOTOMETRICTABLE_H

#include <stdio.h>
#include <math.h>

#include <fitsio.h>

#include <pilmacros.h>

#include <vmtable.h>
#include <vmimage.h>


PIL_BEGIN_DECLS

#define VM_IPC "IPC"


VimosTable *newPhotometricTable();
VimosBool readFitsPhotometricTable(VimosTable *, fitsfile *);
VimosBool writeFitsPhotometricTable(char *, VimosTable *);
VimosBool checkPhotometricTable(VimosTable *);

PIL_END_DECLS

#endif /* VM_PHOTOMETRICTABLE_H */
