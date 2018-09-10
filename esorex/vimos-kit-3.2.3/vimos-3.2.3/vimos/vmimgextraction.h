/* $Id: vmimgextraction.h,v 1.3 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VM_IMGEXTRACTION_H
#define VM_IMGEXTRACTION_H

#include <pilmacros.h>

#include <vmimage.h>
#include <vmtable.h>
#include <vmwcsutils.h>

 
PIL_BEGIN_DECLS

int *VmSearchMatches(VimosTable *, VimosTable *, double, double, double,
		     float, int, int *);
VimosTable *VmImStarCat(VimosImage *, char *, double, int);
VimosTable *VmImDetectObjects(VimosImage *, VimosImage *);
VimosTable *VmImBuildGalaxyTable(VimosTable *, VimosImage *);
VimosTable *VmImBuildStarTable(VimosTable *, float, float);
VimosTable *VmImBuildStarMatchTable(VimosImage *, VimosTable *, VimosTable *,
				    int, double, double, double, float);
VimosTable *VmImBuildStarMatchTable_skyccd(VimosImage *, VimosTable *, VimosTable *,
				    int, double, double, double, float);

PIL_END_DECLS
 
#endif /* VM_IMGEXTRACTION_H */
