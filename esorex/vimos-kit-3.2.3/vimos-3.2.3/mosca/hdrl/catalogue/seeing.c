/* $Id: 4
 *
 * This file is part of the CASU Pipeline utilities
 * Copyright (C) 2015 European Southern Observatory
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: jim $
 * $Date: 2015/08/12 11:16:55 $
 * $Revision: 1.3 $
 * $Name:  $
 */

#include <stdio.h>
#include <math.h>

#include "imcore.h"
#include "util.h"
#include "floatmath.h"

static void sortit (float [], intptr_t);

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out the median seeing
  
    \par Name:
        imcore_seeing
    \par Purpose:
        Work out the median seeing
    \par Description:
        The areal profiles for an array of objects is examined. The point
        where the areal profile falls to half its peak value is found for
        each object and the final seeing estimate is the median of these 
        results
    \par Language:
        C
    \param ap
        The current ap structure
    \param nrows
        The number rows in the object catalogue
    \param ellipt
        The array of ellipticities from the object catalogue
    \param pkht
        The array of peak heights from the object catalogue
    \param areal
        The array of areal profiles from the object catalogue
    \param work
        A work array (should probably allocate this local at some stage)
    \param fwhm
        The output FWHM estimate
    \returns
        Nothing
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern void imcore_seeing(ap_t *ap, intptr_t nrows, float *ellipt, float *pkht, 
                          float **areal, float *work, float *fwhm) {
    intptr_t i,ii,iaper;
    float aper,delaper,area,logf5t,logf2,arg;

    /* Convenience variables */

    logf5t = logf(0.5/ap->thresh);
    logf2 = logf(2.0);

    /* Do the seeing calculation */
/* CONSTANTS: 0.5, 2.0, 3, 0.25, 2.0 */
    ii = 0;
    for (i = 0; i < nrows; i++) {
        if (ellipt[i] < 0.2 && pkht[i] < 30000.0 && pkht[i] > 10.0*ap->thresh) {
            aper = (logf5t + logf(pkht[i]))/logf2 + 1.0;
            iaper = (intptr_t)aper;
            delaper = aper - iaper;
            if (iaper > 0 && iaper < NAREAL && areal[1][i] > 0.0) {
                area = (1.0-delaper)*areal[iaper-1][i] + 
                    delaper*areal[iaper][i];
                work[ii++] = CPL_MATH_2_SQRTPI*sqrtf(area);
            }
        }
    }

    /* Sort the resulting array and choose a location that allows for
       contamination by galaxies */

    if (ii >= 3) {    
        sortit(work,ii);
        *fwhm = work[ii/3 - 1];

        /* Allow for finite pixel size */

        arg = 0.25*CPL_MATH_PI*powf(*fwhm,2.0) - 1;
        *fwhm = 2.0*sqrt(MAX(0.0,arg/CPL_MATH_PI));
    } else 
        *fwhm = 0.0;

}

/**@}*/

static void sortit (float ia[], intptr_t n) {
    intptr_t i, j, ii, jj, ifin;
    float it;
 /* CONSTANTS: 4, 2, 3, */
    jj = 4;
    while (jj < n) 
        jj = 2 * jj;
    jj = MIN(n,(3 * jj)/4 - 1);
    while (jj > 1) {
        jj = jj/2;
        ifin = n - jj;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + jj;
            if (ia[i] <= ia[j]) 
                continue;
            it = ia[j];
            do {
                ia[j] = ia[i];
                j = i;
                i = i - jj;
                if (i < 0) 
                    break;
            } while (ia[i] > it);
            ia[j] = it;
        }
    }
    return;
}

/*

$Log: seeing.c,v $
Revision 1.3  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.2  2015/08/07 13:06:54  jim
Fixed copyright to ESO

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.2  2014/04/09 09:09:51  jim
Detabbed

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
