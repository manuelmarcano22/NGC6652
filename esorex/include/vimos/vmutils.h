/* $Id: vmutils.h,v 1.6 2013-08-23 10:23:02 cgarcia Exp $
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
 * $Date: 2013-08-23 10:23:02 $
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VM_UTILS_H
#define VM_UTILS_H

#include <stdio.h>

#include <pilmacros.h>
#include <pilframeset.h>

#include <vmimage.h>
#include <cpl_propertylist.h>
#include <cxstring.h>

PIL_BEGIN_DECLS

void vimos_print_banner(void);
const char * vimos_get_version(void);
const char * vimos_get_license(void);

int writeDoublePAFEntry(FILE *, char *name, double value);
int writeIntPAFEntry(FILE *, char *name, int value);
int writeStringPAFEntry(FILE *, char *name, char *value);

int UpdateProductDescriptors(VimosImage *ima_in, const char *category);

int applyListSelection(VimosImage **, float *, int, double, double,
		       unsigned int);

int remapFloatsLikeImages(VimosImage **original, VimosImage **sorted,
                          float *array, int imageCount);
int remapDoublesLikeImages(VimosImage **original, VimosImage **sorted,
                           double *array, int imageCount);

void sortN (int ncol, float **ra, int sortCol, int fromRow, int forRows);

cxint vm_plist_update(cpl_propertylist *, cpl_propertylist *, const cxchar *);

int vm_dfs_setup_product_header(PilFrame *, const char *, PilSetOfFrames *);

int getArcLampTimes(VimosImage *, double *);

PIL_END_DECLS

#endif /* VM_UTILS_H */

