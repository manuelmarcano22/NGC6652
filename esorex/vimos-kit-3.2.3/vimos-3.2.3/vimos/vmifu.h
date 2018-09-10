/* $Id: vmifu.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_IFU_H
#define VM_IFU_H

#include <pilmacros.h>
#include <cpl_memory.h> 
#include <cpl_image.h>
#include <cpl_table.h>

PIL_BEGIN_DECLS

void flux_constant();

cpl_image *cpl_image_general_median_filter(cpl_image *, int, int, int);
cpl_image *cpl_image_vertical_median_filter(cpl_image *, 
                                            int, int, int, int, int);

cpl_image *removeBias(cpl_image *, cpl_image *);
cpl_image *removeBiasLevel(cpl_image *image);
int fiberPeak(cpl_image *, int, float *, float *);
cpl_table *ifuIdentify(cpl_image *, int);
int ifuIdentifyUpgrade(cpl_image *, int, float *, cpl_table *, int, int);
cpl_table **ifuTrace(cpl_image *, int, int, int, int, cpl_table *);
cpl_table *ifuTraceDetected(cpl_image *, int, int, int, int, cpl_table *);
cpl_table **ifuFit(cpl_table *, int, float, int);
cpl_image *ifuBack(cpl_image *, int, cpl_table **, int *, int *, int);
cpl_table *ifuGap(cpl_image *, cpl_table *, int, int, int);
int ifuSignal(cpl_table *, int, int);
cpl_table *ifuProfile(cpl_image *, cpl_table *, cpl_table *, cpl_table *);
cpl_table *ifuGauss(cpl_table *, int, int);
cpl_table *rebinProfile(cpl_table *, int, int, double, double);
cpl_table *ifuDetect(cpl_image *, int, float);
cpl_table *ifuFitDetected(cpl_table *, int, float, int);
cpl_table *ifuMatch(cpl_table *, cpl_table *, int, double *, double *);
int ifuFillTracings(cpl_table *coeff, cpl_table *model);
cpl_table *ifuComputeTraces(cpl_table *, int, int, int);
cpl_table *ifuAlign(cpl_table *, cpl_table *, double, double);
cpl_table *ifuSimpleExtraction(cpl_image *, cpl_table *);
cpl_table *ifuExtraction(cpl_image *, cpl_table *);
cpl_table *ifuVerySimpleExtraction(cpl_image *, cpl_table *);
cpl_table *ifuTransmission(cpl_image *, int, int, double *, double *);
int ifuApplyTransmission(cpl_image *image, cpl_table *table);
double *ifuIntegrateSpectra(cpl_table *spectra, int, int);
int ifuImage(cpl_image *, double *, int, int);
int ifuRange(int, double *, double *, double *);
int ifuRangeTransmission(int, double *, double *);
int ifuExtractionParameters(int, int, int, int, int *, int *, int *, int *);
double *ifuFirstIds(int grism, int quadrant, int slit,
                    int *order, double *lambda);
cpl_table *ifuComputeIds(cpl_table *, cpl_table *, double *, int, double, 
                         int, double);
double *ifuComputeIdsBlind(cpl_table *, cpl_table *, double, int, double,
                           double);
int ifuResampleSpectra(cpl_image *, cpl_table *, cpl_table *,
                       int, double, double, double);
double ifuAlignSkylines(cpl_table *, cpl_table *, double, int);
int findCentralFiber(cpl_table *, int);
cpl_image *ifuSubtractSky(cpl_image *);
cpl_image *ifuSumSpectrum(cpl_image *);
int extractIfuFlux(cpl_image *, double, double, double, double *, double *);
int ifuSetZeroLevel(cpl_image *);

PIL_END_DECLS

#endif   /* VM_IFU_H */
