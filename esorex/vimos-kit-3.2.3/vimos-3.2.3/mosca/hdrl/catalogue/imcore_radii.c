/* $Id: imcore_radii.c,v 1.2 2015/08/07 13:06:54 jim Exp $
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
 * $Date: 2015/08/07 13:06:54 $
 * $Revision: 1.2 $
 * $Name:  $
 */

#include <stdio.h>
#include <math.h>

#include "imcore.h"
#include "imcore_radii.h"
#include "floatmath.h"
#include "util.h"
#include "ap.h"

static float fraction (float x, float y, float r_out);

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out the half-light radius for an object
  
    \par Name:
        imcore_halfrad
    \par Purpose:
        Work out the half-light radius for an object
    \par Description:
        Given the the array of core apertures and core fluxes, work out
        the half light radius
    \par Language:
        C
    \param rcores
        The list of aperture radii used
    \param cflux
        The list of fluxes through the aperture radii
    \param halflight
        An estimate of half the light of the object
    \param peak
        The peak flux of the object
    \param naper
        The number of radii used
    \returns
        The half-light radius of the object
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern float imcore_halflight(float rcores[], float cflux[], float halflight,
                              float peak, intptr_t naper) {
    float delr,halfrad;
    intptr_t i,gotone;

    /* Work out the half-light value from either isophotal flux or the 
       flux at Rcore. The find out roughly where the curve of growth
       exceeds this */

    gotone = 0;
    for (i = 0; i < naper; i++) {
        if (cflux[i] > halflight) {
            gotone = 1;
            break;
        }
    }
    if (! gotone) 
        i = naper - 1;

    /* Now work out what the radius of half light is */

    if (i == 0) {
        delr = (cflux[i] - halflight)/MAX(1.0,cflux[i]-peak);
        halfrad = rcores[0]*(1.0 - delr) + delr*sqrt(1.0/CPL_MATH_PI);
    } else {
        delr = (cflux[i] - halflight)/MAX(1.0,(cflux[i] - cflux[i-1]));
        halfrad = rcores[i-1]*delr + rcores[i]*(1.0-delr);
    }
    return(halfrad);
}   

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out the exponential radius for an object
  
    \par Name:
        imcore_exprad
    \par Purpose:
        Work out the exponential radius for an object
    \par Description:
        Given the detection threshold, peak flux and lowest areal 
        profile, work out the exponential radius for an object
    \par Language:
        C
    \param thresh
        The detection threshold
    \param peak
        The peak flux of the object
    \param areal0
        The lowest level areal profile for the object
    \param rcores
        The list of aperture radii used
    \param naper
        The number of radii used
    \returns
        The expoential radius of the object
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern float imcore_exprad(float thresh, float peak, float areal0, 
                           float rcores[], intptr_t naper) {
    float pk,r_t,rad;

    /* Work out the radius... */
/* CONSTANTS: 1.5, 5 */
    pk = MAX(1.5*thresh,peak);
    r_t = sqrtf(areal0/CPL_MATH_PI);
    rad = 5.0*r_t/logf(pk/thresh);
    rad = MAX(r_t,MIN(5.0*r_t,MIN(rad,rcores[naper-1])));
    return(rad);
}
    
/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out the Kron radius for an object
  
    \par Name:
        imcore_kronrad
    \par Purpose:
        Work out the Kron radius for an object
    \par Description:
        Given the lowest areal profile and the circular aperture fluxes
        already done, calculate the Kron radius
    \par Language:
        C
    \param areal0
        The lowest level areal profile for the object
    \param rcores
        The list of aperture radii used
    \param cflux
        The aperture fluxes for each radius
    \param naper
        The number of radii used
    \returns
        The Kron radius of the object
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern float imcore_kronrad(float areal0, float rcores[], float cflux[], 
                            intptr_t naper) {
    intptr_t i,imax;
    float r_t,rad,wt,sum;

    /* Work out the radius... */
/* CONSTANTS: 0.5, 7, 5 */
    r_t = sqrtf(areal0/CPL_MATH_PI);
    rad = 0.5*rcores[0]*cflux[0];
    sum = cflux[0];
    imax = MIN(naper,7);
    for (i = 1; i < imax; i++) {
        wt = MAX(0.0,cflux[i]-cflux[i-1]);
        rad += 0.5*(rcores[i] + rcores[i-1])*wt;
        sum += wt;
    }
    rad /= sum;
    rad = MAX(r_t,MIN(5.0*r_t,MIN(2.0*rad,rcores[naper-1])));    
    return(rad);
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out the Petrosian
  
    \par Name:
        imcore_petrad
    \par Purpose:
        Work out the Petrosian radius for an object
    \par Description:
        Given the lowest areal profile and the circular aperture fluxes
        already done, calculate the Petrosian radius
    \par Language:
        C
    \param areal0
        The lowest level areal profile for the object
    \param rcores
        The list of aperture radii used
    \param cflux
        The aperture fluxes for each radius
    \param naper
        The number of radii used
    \returns
        The Petrosian radius of the object
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern float imcore_petrad (float areal0, float rcores[], float cflux[], 
                            intptr_t naper) {
    intptr_t j;
    float eta,r_t,etaold,r1,r2,r3,r4,r5,r_petr;

    /* Work out petrosian radius */

    r_t = sqrtf(areal0/CPL_MATH_PI);
    eta = 1.0;
    etaold = eta;
    j = 1;
    while (eta > 0.2 && j < naper) {
        etaold = eta;
        r1 = rcores[j]*rcores[j]/(rcores[j-1]*rcores[j-1]) - 1.0;
        r2 = cflux[j]/cflux[j-1] - 1.0;
        eta = r2/r1;
        j++;
    }
    if (j == naper) {
        r_petr = rcores[naper-1];
    } else {
        r1 = rcores[j]*rcores[j];
        r2 = rcores[j-1]*rcores[j-1];
        r3 = rcores[j-2]*rcores[j-2];
        r4 = (etaold - 0.2)/(etaold - eta);
        r5 = (0.2 - eta)/(etaold - eta);
        r_petr = r4*sqrt(0.5*(r1 + r2)) + r5*sqrt(0.5*(r2 + r3));
    }
    r_petr = MAX(r_t,MIN(5.0*r_t,MIN(2.0*r_petr,rcores[naper-1])));
    return(r_petr);
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out the fluxes for special radii
  
    \par Name:
        imcore_flux
    \par Purpose:
        Work out the fluxes for special radii
    \par Description:
        The fluxes for the 'special' radii (Kron etc) are worked out
        by an interpolation of the pre-existing aperture photometry to the
        new radius
    \par Language:
        C
    \param ap
        The current ap structure
    \param parm
        The parameters for each object already detected
    \param nbit
        The number of detected objects in the current Plessey structure.
    \param apers
        The radii of the standard apertures
    \param fluxes
        The fluxes computed through the standard apertures
    \param nr
        The number of special apertures
    \param rcores
        The radii the special apertures
    \param rfluxes
        The fluxes computed through the special apertures.
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

void imcore_flux(ap_t *ap, float parm[IMNUM][NPAR], intptr_t nbit, float apers[], 
                 float fluxes[], intptr_t nr, float rcores[], float rfluxes[]) {
    float *map,t,xj,yj,sumiso,sumcf,delr;
    unsigned char *mflag,mf;
    intptr_t nx,ny;
    intptr_t xmin,xmax,ymin,ymax,ix1,ix2,iy1,iy2,i,j,kk,n;

    /* Set up some local variables */
 
    map = ap->indata;
    mflag = ap->mflag;
    nx = ap->lsiz;
    ny = ap->csiz;

    /* Section for nbit == 1 */

    if (nbit == 1) {

        /* Generate image-blend outer boundaries */

        xmin = parm[0][1] - apers[0] - 0.5;
        xmax = parm[0][1] + apers[0] + 0.5;
        ymin = parm[0][2] - apers[0] - 0.5;
        ymax = parm[0][2] + apers[0] + 0.5;
        ix1 = MAX(0,(intptr_t)xmin-1);
        ix2 = MIN(nx-1,(intptr_t)xmax);
        iy1 = MAX(0,(intptr_t)ymin-1);
        iy2 = MIN(ny-1,(intptr_t)ymax);

        /* Now go through pixel region and add up the contributions inside
           the aperture */

        fluxes[0] = 0.0;
        for(j = iy1; j <= iy2; j++) {
            kk = j*nx;
            for(i = ix1; i <= ix2; i++) {
                mf = mflag[kk+i];
                if (mf == MF_CLEANPIX || mf == MF_OBJPIX || 
                    mf == MF_SATURATED) {
                    t = map[kk+i];   
                    xj = (float)i - parm[0][1] + 1.0;
                    yj = (float)j - parm[0][2] + 1.0;
                    fluxes[0] += fraction(xj,yj,apers[0])*t;
                } 
            }
        }
        if (fluxes[0] <= 0) 
            fluxes[0] = parm[0][0];
        
    /* Section for blended images */
    
    } else {
    
        /* Interpolate circular aperture fluxes */
    
        sumiso = 0.0;
        sumcf = 0.0;
        for (j = 0; j < nbit; j++) {
            sumiso += parm[j][0];
            n = 1;
            while (rcores[n] < apers[j] && n < nr-1)
                n++;
            delr = (rcores[n] - apers[j])/(rcores[n] - rcores[n-1]);
            fluxes[j] = rfluxes[j*nr+n]*(1.0 - delr) + rfluxes[j*nr+n-1]*delr;
            sumcf += fluxes[j];
        }

        /* Constrain the result so that the ratios are the same as for the
           isophotal fluxes */

        for (j = 0; j < nbit; j++) {
            fluxes[j] = sumcf*parm[j][0]/MAX(1.0,sumiso);
            if (fluxes[j] < 0.0)
                fluxes[j] = parm[j][0];
        }
    }
}

/**@}*/

/* returns fraction of pixel bounded by 0 -  r_out
 * x,y coordinates relative to centre
 * Uses linear approximation ok if pixel located >>1 away from centre */

static float fraction (float x, float y, float r_out) {
    float r,t,x_a,x_b,frac,tanao2,cosa,tanp2a,sqrt2o2;

    r = sqrtf(x*x + y*y);
    sqrt2o2 = 0.5*CPL_MATH_SQRT2;

    /* is it worth bothering? */

    if(r > r_out+sqrt2o2) 
        return(0.0);

    /* is it trivially all in? */

    if(r < r_out-sqrt2o2) 
        return(1.0);

    /* bugger - have to do some work then ... ok first ...
     * use 8-fold symmetry to convert to 0-45 degree range */

    x = fabsf(x);
    y = fabsf(y);
    if(y > x) {
        t = x;
        x = y;
        y = t;
    }

    /* If the angles are too close to cardinal points, then fudge something */

    if (x > 0.0 && y > 0.0) {
        tanao2 = 0.5*y/x;
        tanp2a = x/y;
        cosa = x/sqrt(x*x + y*y);
    } else {
        tanao2 = 0.00005;
        tanp2a = 10000.0;
        cosa = 1.0;
    }
/* CONSTANTS: 0.5, 0.00005, 10000 */
    /* only outer radius - compute linear intersections top and bot of pixel */

    x_a = x - tanao2 + (r_out - r)/cosa;
    if(x_a < x+0.5) {

        /* intersects */

        x_b = x + tanao2 + (r_out - r)/cosa;

        /* three cases to consider */

        if(x_a < x-0.5)
            frac = 0.5*MAX(0.0,x_b-(x-0.5))*MAX(0.0,x_b-(x-0.5))*tanp2a;
        else {
            if(x_b > x+0.5)
                frac = 1.0 - 0.5*(x+0.5-x_a)*(x+0.5-x_a)*tanp2a;
            else
                frac = 0.5-(x-x_a)+0.5*(x_b-x_a);
        }
    } else  /* missed entirely */
        frac = 1.0;

    return(frac);
}

/*

$Log: imcore_radii.c,v $
Revision 1.2  2015/08/07 13:06:54  jim
Fixed copyright to ESO

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.3  2015/03/12 09:16:51  jim
Modified to remove some compiler moans

Revision 1.2  2014/04/09 09:09:51  jim
Detabbed

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
