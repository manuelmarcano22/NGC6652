/* $Id: vmimgastrometry.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_IMGASTROMETRY_H
#define VM_IMGASTROMETRY_H

#include <pilmacros.h>

#include <vmimage.h>
#include <vmtable.h>
#include <vmtablearray.h>


PIL_BEGIN_DECLS

VimosBool computeCO(VimosImage *, const char *, int, float, float);
VimosBool computeCD(char **, unsigned int);

VimosBool VmAstroComputeCD(VimosTableArray *, const VimosImage *);
VimosBool VmAstroComputeCO(VimosImage *, float , float, VimosTable *,
			   unsigned int, double);

VimosTable *resetStarMatchTable(int, int, int);
VimosTable *resetAstroTable(int, char *);
VimosTable *shiftStarMatch(VimosTable **);
VimosBool fitCO(VimosImage *, VimosTable *, VimosTable *, int, double,
		double, double, float, double);

int calcres(VimosTable *, VimosTable *, int *, int, double *);

PIL_END_DECLS

#endif /* VM_IMGASTROMETRY_H */