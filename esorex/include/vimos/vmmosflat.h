/* $Id: vmmosflat.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_MOSFLAT_H
#define VM_MOSFLAT_H

#include <pilmacros.h>

#include <vmimage.h>
#include <vmtable.h>
#include <vmextractiontable.h> 


PIL_BEGIN_DECLS

VimosImage *VmSpNormFF(VimosImage *, VimosExtractionTable *, int, int, int,
                       FilterMethod, int, int);
VimosImage *VmSpNormPoly(VimosImage *, VimosExtractionTable *, int, int);
VimosImage *VmSpNormSmooth(VimosImage *, VimosExtractionTable *, int, int,
                           FilterMethod);

VimosImage **VmSpStackFF(VimosImage **, int, VimosExtractionTable *, int);

VimosImage *VmSpApplyFF(VimosImage *, VimosImage *, VimosExtractionTable *);


VimosBool determineExposedMosArea(VimosTable *,VimosExtractionTable *,
                                  int *, VimosDpoint **, VimosDpoint **,
                                  VimosDpoint **zeroY);
VimosExtractionSlit *determineExposedSlits(VimosTable *,
                                           VimosExtractionSlit **,
                                           int *nSlits);
VimosBool slitMinMaxY(VimosExtractionSlit *, float *, float *);

PIL_END_DECLS

#endif /* VM_MOSFLAT_H */
