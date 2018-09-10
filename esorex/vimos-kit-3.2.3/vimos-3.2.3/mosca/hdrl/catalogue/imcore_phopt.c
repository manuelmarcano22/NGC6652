/* $Id: imcore_phopt.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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
#include <stdlib.h>
#include <math.h>

#include <cpl.h>

#include "floatmath.h"
#include "util.h"
#include "imcore.h"
#include "ap.h"

/* Function Prototypes */

static void dchole (double a[IMNUM+1][IMNUM+1], double b[IMNUM+1], intptr_t n);
static double fraction (double x, double y, double r_out);

/* does multiple profile fitting to determine intensities */

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do multiple profile fitting
  
    \par Name:
        imcore_phopt
    \par Purpose:
        Do multiple profile fitting
    \par Description:
        Given a Plessey array and some parameters determined from a moments
        analysis for each of the objects detected in the array, this routine
        does multiple profile fitting for the given aperture set
    \par Language:
        C
    \param ap
        The current ap structure
    \param parm
        The input/output object parameters
    \param nbit
        The number of objects detected in the current Plessey structure
    \param naper
        The number of apertures 
    \param apertures
        Array of aperture radii
    \param cflux
        Array of aperture fluxes
    \param badpix
        Array saying how many bad pixels were included in the data for
        each object at each radius
    \param nrcore
        The index of the apertures array that defines where the radius = Rcore
    \param avconf
        TODO
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

extern void imcore_phopt(ap_t *ap, float parm[IMNUM][NPAR], intptr_t nbit, 
                         intptr_t naper, float apertures[], float cflux[], 
                         float badpix[], intptr_t nrcore, float avconf[]) {
    double aa[IMNUM+1][IMNUM+1],bb[IMNUM+1];
    float d,arg,*map,rcirc,ff;
    float cn,parrad,xmin,xmax,xi,yi,ymin,ymax;
    float t,xj,yj,cnsq,tj,xk,yk,tk;
    intptr_t i,ii,j,kk,ix1,ix2,iy1,iy2,nx,ny,k,iaper;
    unsigned char *mflag,mf;
    float *conf;

    /* Set up some local variables */

    map = ap->indata;
    conf = ap->confdata;
    mflag = ap->mflag;
    nx = ap->lsiz;
    ny = ap->csiz;   

    /* Loop for each of the apertures */

    for (iaper = 0; iaper < naper; iaper++) {
        rcirc = apertures[iaper];
        parrad = rcirc + 0.5;
        cn = 1.0/(CPL_MATH_PI*rcirc*rcirc);  /* profile normalising constant */
        cnsq = cn*cn;
    
        /* set up covariance matrix - analytic special case for cores */

        for(i = 0; i < nbit; i++) {
            aa[i][i] = cn;                 /* overlaps totally area=pi*r**2 */
            if(nbit > 1) {
                xi = parm[i][1];
                yi = parm[i][2];
                for(j = i+1; j < nbit; j++) {
                    d = sqrtf((xi-parm[j][1])*(xi-parm[j][1]) 
                        + (yi-parm[j][2])*(yi-parm[j][2]));
                    if(d >= 2.0*rcirc) {
                        aa[j][i] = 0.0;
                        aa[i][j] = aa[j][i];
                    } else {
                        arg = d/(2.0*rcirc);
                        aa[j][i] = cnsq*2.0*rcirc*rcirc*
                            (acosf(arg)-arg*(sqrtf(1.0-arg*arg)));
                        aa[i][j] = aa[j][i];
                    }
                }
            }
        }

        /* clear accumulators */

        for(i = 0; i < nbit; i++) 
            bb[i] = 0.0;

        /* generate image-blend outer boundaries */

        xmin = 1.0e6;
        xmax = -1.0e6;
        ymin = 1.0e6;
        ymax = -1.0e6;
        for(i = 0; i < nbit; i++) {
            xi = parm[i][1];
            yi = parm[i][2];
            xmin = MIN(xmin, xi);
            xmax = MAX(xmax, xi);
            ymin = MIN(ymin, yi);
            ymax = MAX(ymax, yi);
        }
        ix1 = MAX(0,(intptr_t)(xmin-parrad)-1);
        ix2 = MIN(nx-1,(intptr_t)(xmax+parrad));
        iy1 = MAX(0,(intptr_t)(ymin-parrad)-1);
        iy2 = MIN(ny-1,(intptr_t)(ymax+parrad));

        /* now go through pixel region */

        for(ii = iy1; ii <= iy2; ii++) {
            kk = ii*nx;
            for(i = ix1; i <= ix2; i++) {
                mf = mflag[kk+i];
                if (mf == MF_ZEROCONF || mf == MF_STUPID_VALUE) {
                    for (j = 0; j < nbit; j++) {
                        xj = i - parm[j][1] + 1.0;
                        yj = ii - parm[j][2] + 1.0;
                        tj = fraction(xj,yj,rcirc);
                        aa[j][j] -= tj*tj*cnsq;
                        for (k = j + 1; k < nbit; k++) {
                            xk = i - parm[k][1] + 1.0;
                            yk = ii - parm[k][2] + 1.0;
                            tk = fraction(xk,yk,rcirc);
                            aa[k][j] -= tk*tj*cnsq;
                            aa[j][k] = aa[k][j];
                        }
                        if (iaper == nrcore)
                            badpix[j] += tj;
                    }
                } else if (mf == MF_CLEANPIX || mf == MF_OBJPIX ||
                           mf == MF_SATURATED) {
                    t = map[kk+i];   
                    for(j = 0; j < nbit; j++) {
                        xj = i - parm[j][1] + 1.0;
                        yj = ii - parm[j][2] + 1.0;
                        ff = fraction(xj,yj,rcirc);
                        bb[j] += ff*t;
                        if (iaper == nrcore)
                            avconf[j] += ff*(float)conf[kk+i];
                    }
                } 
            }
        }

        /* Trivial solution for single object */

        if (nbit == 1) {
            cflux[iaper] = bb[0];

        /* solve for profile intensities */

        } else {
            for (i = 0; i < nbit; i++) 
                aa[i][i] = MAX(aa[i][i],cnsq);
            dchole(aa,bb,nbit);
            for(i = 0; i < nbit; i++) 
                cflux[i*naper+iaper] = cn*bb[i];
        }
    }
}

/**@}*/

/* CHOLEsky decomposition of +ve definite symmetric matrix to solve Ax = b */

static void dchole (double a[IMNUM+1][IMNUM+1], double b[IMNUM+1], intptr_t n) {
  double sum, l[IMNUM+1][IMNUM+1], y[IMNUM+1];
  double aveigv, offset;
  intptr_t i, j, k;

restart:
    l[0][0] = sqrt(a[0][0]);

    for(k = 1; k < n; k++) {
        for(j = 0; j <= k-1; j++) {
            sum = a[j][k];
            if(j != 0) 
                for(i = 0; i <= j-1; i++) 
                    sum -= l[i][k]*l[i][j];
            l[j][k] = sum/l[j][j];
        }
        sum = a[k][k];
        for(i = 0; i <= k-1; i++) 
            sum -= l[i][k]*l[i][k];
        if(sum <= 0.0) {
/*          fprintf(stderr, "dchole: warning: matrix ill-conditioned\n"); */
            aveigv = a[0][0];
            for(i = 1; i < n; i++) 
                aveigv += a[i][i];
            /* max eigenvalue < trace */
            offset = 0.1*aveigv/((double) n);
            for(i = 0; i < n; i++) 
                a[i][i] += offset;
/*          fprintf(stderr, "dchole: Offset added to diagonal = %f\n", offset); */
            goto restart;
        }
        l[k][k] = sqrt(sum);
    }

    /* solve Ly = b */

    y[0] = b[0]/l[0][0];
    for(i = 1; i < n; i++) {
        sum = b[i];
        for(k = 0; k <= i-1; k++) 
            sum -= l[k][i]*y[k];
        y[i] = sum/l[i][i];
    }

    /* solve L(T)x = y */

    b[n-1] = y[n-1]/l[n-1][n-1];
    for(i = n-2; i >= 0; i--) {
        sum = y[i];
        for(k = i+1; k < n; k++) 
            sum -= l[i][k]*b[k];
        b[i] = sum/l[i][i];
    }
}

/* returns fraction of pixel bounded by 0 -  r_out
 * x,y coordinates relative to centre
 * Uses linear approximation ok if pixel located >>1 away from centre */

static double fraction (double x, double y, double r_out) {
    double r,r2,t,x_a,x_b,frac,tanao2,cosa,tanp2a,sqrt2o2;

    r2 = (x*x + y*y);
    sqrt2o2 = 0.5*CPL_MATH_SQRT2;

    /* is it worth bothering? */

    if(r2 > (r_out+sqrt2o2) * (r_out+sqrt2o2))
        return(0.0);

    /* is it trivially all in? */

    if(r_out >= sqrt2o2 && r2 < (r_out-sqrt2o2) * (r_out-sqrt2o2))
        return(1.0);

    /* bugger - have to do some work then ... ok first ...
     * use 8-fold symmetry to convert to 0-45 degree range */

    x = fabs(x);
    y = fabs(y);
    if(y > x) {
        t = x;
        x = y;
        y = t;
    }

    r = sqrt(r2);

    /* If the angles are too close to cardinal points, then fudge something */

    if (x > 0.0 && y > 0.0) {
        tanao2 = 0.5*y/x;
        tanp2a = x/y;
        cosa = x/r;
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

$Log: imcore_phopt.c,v $
Revision 1.3  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.2  2015/08/07 13:06:54  jim
Fixed copyright to ESO

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.3  2015/01/09 11:41:29  jim
Added average confidence calculation

Revision 1.2  2014/04/09 09:09:51  jim
Detabbed

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
