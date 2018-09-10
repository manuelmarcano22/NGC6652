/* $Id: imcore_extend.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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
#include <string.h>

#include <cpl.h>

#include "imcore.h"
#include "floatmath.h"
#include "util.h"

#define NACC 10
#define NCOEF 4

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do aperture integration
  
    \par Name:
        imcore_extend
    \par Purpose:
        Do aperture integration
    \par Description:
        The integrated flux of an object is calculated using matched
        ellipses
    \par Language:
        C
    \param ap
        The current ap structure
    \param xniso
        The isophotal flux
    \param xbar
        The X position of the object
    \param ybar
        The Y position of the object
    \param sxx
        Second moment in X
    \param syy
        Second moment in Y
    \param sxy
        Second moment cross term
    \param areal0
        The first areal profile 
    \param tmax
        The peak flux of the object
    \param ttotal
        The output total integrated flux
    \retval CASU_OK
        If all went OK. Currently this is the only value.
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_extend(ap_t *ap, float xniso, float xbar, float ybar, 
                         float sxx, float sxy, float syy, float areal0, 
                         float tmax, float *ttotal) {
    float srr,ecc,xx,ctheta,stheta,a,b,stretch,rad,sfac,climsq,clim;
    float pt1,pt2,pt3,c,pa,pb,pc,arg1,xliml,xlimu,y,t,x,xnew,ynew,ellrad;
    float xmax,xlim1,xlim2,xcord[NACC],xdat[NACC],polycf[NCOEF],rt1,rt2,t1;
    float xlimit,theta,accum[NACC],*map,skysig,thresh;
    intptr_t jmin,jmax,imax,imin,jj,kk,ii,iupd,j,i,ir,nx,ny;
    unsigned char *mflag;

    /* Initialise a few things */

    map = ap->indata;
    nx = ap->lsiz;
    ny = ap->csiz;
    skysig = ap->sigma;
    thresh = ap->thresh;
    mflag = ap->mflag;

    /* Calculate the eccentricity and position angle of the object */
/* CONSTANTS: 0.5, 4, 16, 0.9, 5, 3, 2 */
    srr = MAX(0.5,sxx+syy);
    ecc = sqrtf((syy-sxx)*(syy-sxx) + 4.0*sxy*sxy)/srr;
    ecc = MIN(0.9,ecc);
    xx = 0.5*(1.0 + ecc)*srr - sxx;
    if (sxy == 0)
        theta = 0.0;
    else
        if (xx == 0.0) 
            theta = CPL_MATH_PI_2;
        else 
            theta = atanf(sxy/xx);
    ctheta = cos(theta);
    stheta = sin(theta);

    /* Eccentricity modified by noise effect.  NB: 50 == 16*pi */

    ecc = sqrtf(MAX((syy-sxx)*(syy-sxx) 
                    - 16.0*CPL_MATH_PI*skysig*srr*srr*srr/(xniso*xniso)
                    + 4.0*sxy*sxy,0.0))/srr;
    ecc = MIN(0.9,ecc);
    
    /* Set initial aperture to be isophotal area */

    a = sqrtf(srr*(1.0 + ecc));
    b = sqrtf(srr*(1.0 - ecc));
    stretch = sqrt(areal0/(CPL_MATH_PI*a*b));

    /* Number of isophotal radii to extend */

    rad = MAX(1.1,(tmax - skysig)/thresh);
    sfac = MIN(5.0,MAX(2.0,3.0/sqrtf(logf(rad))));
    a *= sfac*stretch;
    b *= sfac*stretch;

    /* Clear accumulator */

    memset(accum,0,NACC*sizeof(float));
    
    /* Generate image boundaries. First for y */

    climsq = (a*ctheta)*(a*ctheta) + (b*stheta)*(b*stheta);
    climsq = MAX(1.0,climsq);
    clim = sqrtf(climsq);
    pt1 = sinf(2.0*theta)*(b*b-a*a);
    pt2 = (b*ctheta)*(b*ctheta) + (a*stheta)*(a*stheta);
    pt3 = (a*b)*(a*b);
    jmin = MAX(1,(intptr_t)(ybar - clim));
    jmax = MIN(ny,(intptr_t)(ybar + clim + 1.0));
    for (jj = jmin; jj <= jmax; jj++) {

        /* Now for x */

        kk = (jj-1)*nx;
        c  = (float)jj - ybar;
        pa = climsq;
        pb = pt1*c;
        pc = pt2*c*c - pt3;
        arg1 = pb*pb - 4.0*pa*pc;
        arg1 = sqrtf(MAX(arg1,0.0));
        xliml = (-pb - arg1)/(2.0*pa);
        xlimu = (-pb + arg1)/(2.0*pa);
        imin = MAX(1,(intptr_t)(xbar + xliml));
        imax = MIN(nx,(intptr_t)(xbar + xlimu + 1.0));
        y = c;
        for(ii = imin; ii <= imax; ii++) {
            if (mflag[kk+ii-1] == MF_CLEANPIX || mflag[kk+ii-1] == MF_OBJPIX ||
                mflag[kk+ii-1] == MF_SATURATED) {
                t = map[kk+ii-1];
                x = (float)ii - xbar;

                /* Accumulate elliptical isophotal areas */

                xnew = x*ctheta - y*stheta;
                ynew = x*stheta + y*ctheta;
                ellrad = 2.0*sqrtf((ynew/a)*(ynew/a) + (xnew/b)*(xnew/b));
                iupd = ((intptr_t)((2.0-ellrad)*(float)NACC)) + 1;
                iupd = MAX(1,iupd);
                iupd = MIN(NACC,iupd);
                for(j = 1; j <= iupd; j++) 
                    accum[NACC-j] += t;
            }
        }
    }

    /* Now find limiting intensity */

    if (xniso < 0.0) 
        for(i = 0; i < NACC; i++) 
            accum[i] = -accum[i];
    imcore_median(accum,NACC,3);
    xmax = 0.0;
    xlim1 = -1.0;
    xlim2 = -1.0;
    for(i = 0; i < NACC; i++) {
        xcord[i] = i+1;
        xmax = MAX(xmax,accum[i]);
        xdat[i] = accum[i];
    }
    imcore_polynm(xdat,xcord,NACC,polycf,NCOEF,0);
    pa = polycf[1];
    pb = polycf[2]*2.0;
    pc = polycf[3]*3.0;
    arg1 = sqrtf(MAX(0.0,pb*pb - 4.0*pa*pc));
    if (pc != 0.0) {
        rt1 = (-pb + arg1)/(2.0*pc);
        rt2 = (-pb - arg1)/(2.0*pc);
        if(rt1 < (float)NACC && rt1 > 1.0) {
            ir = (intptr_t)rt1;
            t1 = rt1 - (float)ir;
            xlim1 = (1.0 - t1)*accum[ir-1] + t1*accum[ir];
        }
        if(rt2 < (float)NACC && rt2 > 1.0) {
            ir = (intptr_t)rt2;
            t1 = rt2 - ir;
            xlim2 = (1.0 - t1)*accum[ir-1] + t1*accum[ir];
        }
    }
    xlimit = MAX(xlim1,xlim2);
    if(xlimit < 0.0) 
        xlimit = xmax;

    /* Update total intensity */

    if(xniso < 0.0) 
        xlimit = -xlimit;
    *ttotal = xlimit;
    return(CASU_OK);
}

/**@}*/

/*

$Log: imcore_extend.c,v $
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
