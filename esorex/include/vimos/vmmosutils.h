/* $Id: vmmosutils.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_MOSUTILS_H
#define VM_MOSUTILS_H

#include <pilmacros.h>

#include <cpl_table.h>

#include <vmimage.h>
#include <vmmatrix.h>
#include <vmtable.h>
#include <vmextractiontable.h>


PIL_BEGIN_DECLS

int modelWavcal(VimosExtractionTable *, int);

int extractSpecLevel(VimosImage *, VimosExtractionSlit *, double wave, 
                     int dy, double *flux);

int extractSpecFlux(VimosImage *, VimosExtractionSlit *, double wave,
                    int dy, double *flux, double *error);

int spectralResolution(VimosImage *, float, double *, double *, int);

int testLineSaturation(VimosImage *, VimosTable *);

double distortionsRms(VimosImage *, VimosTable *, double);
double distortionsRms_CPL(VimosImage *, cpl_table *, double);

int getGrism(VimosImage *);
int getGrismAgain(VimosTable *);

int
findCentralPosition(VimosImage *, VimosDescriptor *, double X, double Y, 
                    double slitLength, float width, VimosTable *,
                    double *deltaX, double *deltaY);

int alignWavePattern(VimosImage *, double X, double Y, double slitLength,
                     double *deltaX, double *deltaY);


void findSpectrumBorders(VimosFloatArray *, double *, double *, int);

int *sortByShutterPosition(VimosImage **, int imageCount, int *shutPosCount);

char *createSpectralDistModelsPAF(VimosDescriptor *, char *namePAF);
char *createSpectralDistPAF(VimosDescriptor *descs, char *namePAF);
char *createIdsPAF(VimosDescriptor *, char *namePAF);

PIL_END_DECLS

#endif /* VM_MOSUTILS_H */
