/* $Id: vmwcsutils.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_WCSUTILS_H
#define VM_WCSUTILS_H

#include <pilmacros.h>

#include <vimoswcs.h>

#include <vmtable.h>
#include <vmimage.h>


PIL_BEGIN_DECLS

/* FIXME:
 *   Functions NOT part of the libwcs API, i.e. not declared there, are used
 *   here. This a dirty hack and not acceptable.
 *   At least we do a complete declaration here to be able to visually compare
 *   the actual definitions and the one we expect. The correct way would have
 *   been to reimplement these functions ! (RP)
 */

extern void setnfit(int nfit);
extern void setresid_refine(int refine);

extern int FitPlate(struct WorldCoor *wcs, double *x, double *y, double *x1,
                    double *y1, int np, int ncoeff0, int debug);
extern int FitMatch(int nmatch, double *sbx, double *sby, double *gbra,
                    double *gbdec, struct WorldCoor *wcs, int debug);


/*
 * Maximum number of coefficents for distortion function
 */

#define COMMENT_LENGTH 80
#define VMMAXPAR 20

/*--------------------------------------*
 * Module project.c
 * changes from version 1.5 (BG)
 * VmFitMatch, VmFitPlate, VmPlateFunc moved here from VmImSTarMatch
 *--------------------------------------*/

void wcstopix(int n_ele, VimosTable *a_star, struct WorldCoor *wcs);
void pixtowcs(int n_ele, VimosTable *o_star, struct WorldCoor *wcs);

VimosBool vimosFitMatch(struct WorldCoor *wcs, VimosTable *o_star, int nmatch);
VimosBool vimosFitPlate(struct WorldCoor *wcs, VimosTable *o_star,
                        VimosTable *a_star, int nmatch, int ncoeff,
                        double * chisq);

/*--------------------------------------*
 * Module io_imahead.c
 * changes from version 1.5 (BG)
 * upheader moved here from VmImSTarMatch
 * computeVirtualPixels getCcdSky added 
 * rdfits and rdfits2 taken out (not used/needed)
 * rdimage deeply changed
 *--------------------------------------*/

struct WorldCoor * rdimage(VimosDescriptor *descs);

VimosBool computeVirtualPixels(VimosDescriptor *descs,VimosTable *ostar,
                               unsigned int flag, double tolerance);

int upheader (VimosImage * image,  struct WorldCoor * wcs, double rms[]);

PIL_END_DECLS

#endif /* VM_WCSUTILS_H */
