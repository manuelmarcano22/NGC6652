/* $Id: vmmoswavecalib.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_MOSWAVECALIB_H
#define VM_MOSWAVECALIB_H

#include <pilmacros.h>

#include <vmmatrix.h>
#include <vmimage.h>
#include <vmtable.h>
#include <vmlinecatalog.h>
#include <vmdistmodels.h>


PIL_BEGIN_DECLS

enum _EXTR_METHOD_
{
  EXTR_UNDEF = 0,
  EXTR_LOCAL,
  EXTR_GLOBAL,
  EXTR_INPUT
};

typedef enum _EXTR_METHOD_ ExtrMethod;


int VmSpDerDispOld(VimosImage *, VimosExtractionTable *, VimosTable *, int, int);
int VmSpDerDisp(VimosImage *, VimosExtractionTable *, VimosTable *, int, float);
int VmSpDispMatrix(VimosExtractionTable *, VimosTable *, int);
int VmSpCalShifts(VimosImage *, VimosTable *, VimosExtractionTable *, int, 
		  int, int);

VimosDpoint *getWavIntervals(VimosTable *, float);
void forgetWavIntervals(VimosDpoint *);
double computeMatchIndex(VimosDistModel1D *, VimosDpoint *, VimosFloatArray *,
                         int);
VimosFloatArray *equalizeSpectrum(VimosFloatArray *);

VimosDistModel1D *findMaxMatchIndex(VimosDistModel1D *, VimosTable *, float,
                                    VimosFloatArray *, int);

int findClosestPeak(float *, int);

double **identPeaks(double *, int, double *, int, 
                    double, double, double, int *);

double *collectPeaks(float *, int, float, float, int *);
double *collectPeaks_double(double *, int, float, float, int *);


PIL_END_DECLS

#endif /* VM_MOSWAVECALIB_H */
