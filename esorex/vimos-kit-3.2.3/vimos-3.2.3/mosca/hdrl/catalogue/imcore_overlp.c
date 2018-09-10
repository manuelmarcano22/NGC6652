/* $Id: imcore_overlp.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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

#include "ap.h"
#include "util.h"
#include "imcore.h"
#include "floatmath.h"

#define NITER 6

static void moments_thr(ap_t *, float [], intptr_t []);
static void sort_on_zsm_rev(intptr_t, plstruct *);
static void update_ov(float [], float, float, float, float);
static void check_term(ap_t *, intptr_t *, float [IMNUM][NPAR+1],
		       intptr_t [IMNUM][2], 
                       intptr_t *);

static float oldthr;
static float curthr;
static float nexthr;
static float lasthr;
static float xbar_start;
static float ybar_start;

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Deblend overlapping images
  
    \par Name:
        imcore_overlp
    \par Purpose:
        Deblend overlapping images
    \par Description:
        An array of pixels that are believed to be part of a single large
        object are anaylsed with successively higher threhsolds to see if
        they resolve into multiple objects
    \par Language:
        C
    \param ap
        The current input ap structure
    \param parm
        The parameter array for the deblended objects
    \param nbit
        The output number of objects found in the deblended object
    \param xbar
        The X position of the input object
    \param ybar
        The Y position of the input object
    \param total
        The total flux of the input object
    \param npix
        The number of pixels in the original object
    \param tmax
        The peak flux of the original object
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

extern void imcore_overlp(ap_t *ap, float parm[IMNUM][NPAR], intptr_t *nbit, 
                          float xbar, float ybar, float total, intptr_t npix, 
                          float tmax) {
    plstruct *pl;
    intptr_t npl,ipix,ipixo2,npl2,nbitprev,nobj,toomany,i,isnew,k,kk,j,ibitx[IMNUM];
    intptr_t ibity[IMNUM],iwas,iupdate[IMNUM],npl3,iter,lastone,ic,jj;
    intptr_t ii,conv,ipks[IMNUM][2];
    float fconst,offset,tmul,smul,xintmn,itmaxlim,algthr,radmax,xb,yb,radius2;
    float results[IMNUM][NPAR+1],distmax,dx,dy,parmnew[IMNUM][NPAR],sumint;
    float xlevol,radold,slope,xx,xlevel,radius,xdat[NAREAL+1],xcor[NAREAL+1];
    float dlbydr,wt,dist,xeff,polycf[3],ttt,radthr,delb,deli,ratio;
    float bitx[IMNUM],bitl[IMNUM],sxx,syy;
    ap_t ap2;

    /* Initialise a few variables */

    pl = ap->plarray;
    npl = ap->npl_pix;
    ipix = ap->ipnop;
    oldthr = ap->thresh;
    fconst = ap->fconst;
    offset = ap->areal_offset;
    xbar_start = xbar;
    ybar_start = ybar;

    /* Initialise some constants that you might need later */
/* CONSTANTS: 1.2589678, 2.5, 2,0.9 */
    tmul = 1.2589678;             /* 1/4 mag deblending contour increment */
    smul = 2.5;                   /* starting contour increment */
    ipixo2 = MAX(2,(ipix + 1)/2); /* ipix is the minimum image pixel size */
    xintmn = oldthr*ipixo2;       /* minimum intensity for fragments */
    itmaxlim = 0.9*tmax;          /* upper limit deblending 90% of peak */
    lasthr = itmaxlim;
    algthr = logf(oldthr);        /* Convenient bit of shorthand */
    radmax = sqrtf(((float)npix)/CPL_MATH_PI); /* Isophotal size of input data array */

    /* Get a maximum of IDBLIM points brighter than the new detection threshold
       by reverse sorting the input array. If there are still more than IDBLIM
       above the threshold, then revise the thresold until there aren't. Then
       use the variable npl2 to stop the rest of the routine accessing any of
       the fainter pixels in the list. This is to restrict processing time
       for large extended objects */
 
    curthr = smul*oldthr;
    sort_on_zsm_rev(npl,pl);
    while (1) {
        npl2 = 0;
        while (pl[npl2].zsm > curthr && npl2 < npl-1)
            npl2++;    
        if (npl2 > IDBLIM) 
            curthr += oldthr;
        else
            break;
    }

    /* If there are fewer pixels above the new threshold than the minimum 
       specified in the input parameters, then you have no reason to be here */

    if (npl2 < ipix) {
        *nbit = 1;
        return;
    }

    /* Get a new ap structure */

    ap2.lsiz = ap->lsiz;
    ap2.csiz = ap->csiz;
    ap2.multiply = 1;
    ap2.ipnop = ipixo2;
    ap2.areal_offset = offset;
    ap2.fconst = fconst;
    ap2.mflag = cpl_calloc((ap2.lsiz)*(ap2.csiz),sizeof(unsigned char*));

    /* Main analysis loop at new thresholds */

    *nbit = 0;
    nbitprev = 0;
    imcore_apinit(&ap2);
    while (1) {
        nexthr = MAX(curthr+oldthr,curthr*tmul);
        
        /* Locate objects in this cluster */

        ap2.thresh = curthr;
        imcore_apclust(&ap2,npl2,pl);
        check_term(&ap2,&nobj,results,ipks,&toomany);
        imcore_apreinit(&ap2);
        if (nobj == 0)
            break;

        /* For each image check the number of points above the next threshold
           and flag. Do a moments analysis of each object */

        for (i = 0; i < nobj; i++) {

            /* Ok, has this object already been detected?  If so, then
               load the new results into the parmnew array. We'll check
               whether it's worthwhile updating the master parameters
               list later */

            isnew = 1;
            xb = results[i][1];
            yb = results[i][2];
            sxx = MAX(1.0,results[i][4]);
            syy = MAX(1.0,results[i][6]);
            for (k = 0; k < nbitprev; k++) {
                dx = xb - parm[k][1];
                dy = yb - parm[k][2];
                radius2 = dx*dx/sxx + dy*dy/syy;
                if ((ibitx[k] == ipks[i][0] && ibity[k] == ipks[i][1]) || 
                    radius2 < 1.0) {
                    isnew = 0;
                    for (kk = 0; kk < NPAR; kk++)
                        parmnew[k][kk] = results[i][kk];
                    break;
                }
            }

            /* If this is a new one and it's above the minimum threshold 
               then check to make sure you don't already have too many. 
               If you do, then flag this and break out to the next iteration.
               If there is room for another, then store the moments analysis
               profile */

            if (isnew && results[i][0] > xintmn) {
                if (*nbit >= IMNUM) {
                    *nbit = IMNUM;
                    toomany = 1;
                    break;
                }
                ibitx[*nbit] = ipks[i][0];
                ibity[*nbit] = ipks[i][1];
                for (kk = 0; kk < NPAR; kk++)
                    parm[*nbit][kk] = results[i][kk];
                (*nbit)++;
            }
        } /* End of object loop */

        /* If too many objects were found, then skip the next bit...otherwise
           go through and update parameters if necessary. This block of
           code is a bit of a place holder waiting for something better to
           be worked out...*/

        if (! toomany) {
            if (*nbit > nbitprev && nbitprev > 0) {
                for (i = 0; i < nbitprev; i++) 
                    iupdate[i] = 0;
                for (j = nbitprev; j < *nbit; j++) {
                    distmax = 0.0;
                    iwas = 0;
                    for (i = 0; i < nbitprev; i++) {
                        if (parmnew[i][0] > 0.0) {
                            dx = parmnew[i][1] - parm[i][1];
                            dy = parmnew[i][2] - parm[i][2];
                            radius2 = dx*dx + dy*dy;
                            if (radius2 > distmax) {
                                iwas = i;
                                distmax = radius2;
                            }
                        }
                    }
                    iupdate[iwas] = 1;
                }
                for (i = 0; i < nbitprev; i++)
                    if (iupdate[i] == 1 && parmnew[i][0] > 0.0) 
                        for (j = 0; j < NPAR; j++) 
                            parm[i][j] = parmnew[i][j];
            }

            /* Reset the update flag and prepare for next iteration*/

            for (i = 0; i <= *nbit; i++)
                parmnew[i][0] = -1.0;
            nbitprev = *nbit;
        }

        /* Where do we cut in the list now? */


        npl3 = 0;
        while (pl[npl3].zsm > nexthr && npl3 < npl2-1)
            npl3++;    
        npl2 = npl3;

        /* Now, do we need to move onto the next threshold? */

        if (npl2 == 0 || toomany || nexthr >= itmaxlim) 
            break;

        /* If so, then reset some variables and continue */

        curthr = nexthr;

    } /* End of main analysis loop */

    /* Free up some workspace */

    cpl_free(ap2.mflag);
    imcore_apclose(&ap2);

    /* If there is only one then we can exit now */

    if (*nbit == 1)
        return;

    /* Find out which images terminated properly and remove those that didn't */
    j = -1;
    for (k = 0; k < *nbit; k++) {
        /* Commented this out as checking for terminations seems to miss some
           if the total flux above the threshold for an object is negative */

/*         if (ibitl[k] == 1 && parm[k][0] > xintmn) { */ 
        if (parm[k][0] > xintmn) {
            j++;
            if (j != k)
                for (i = 0; i < NPAR; i++)
                    parm[j][i] = parm[k][i];
        }
    }
    *nbit = j + 1;
    for (j = 0; j < *nbit; j++) {
        bitx[j] = 0.0;
        bitl[j] = 0.0;
    }

    /* For each image find true areal profile levels and iterate to find 
       local continuum */
    /* CONSTANTS: 0.5, 5, 2,0.01,50.,3,7 */
    iter = 0;
    sumint = 0.0;
    lastone = 0;
    while (iter < NITER) {
        iter++;
        
        /* Loop for each of the objects and create a level vs radius array */

        for (k = 0; k < *nbit; k++) {
            if (parm[k][0] < 0.0)
                continue;
            xlevol = logf(parm[k][7] + parm[k][3] - bitl[k]); /* Pk + newthresh - cont */
            xlevel = xlevol;
            radold = 0.0;
            radius = radold;
            slope = 1.0;
            ic = 0;
            for (i = 1; i <= NAREAL; i++) {
                jj = NPAR - i;
                ii = NAREAL - i;
                xx = (float)ii + offset;
                if (parm[k][jj] > 0.5) {
                    if (ii == 0) 
                        xlevel = logf(parm[k][3] - bitl[k] + 0.5);
                    else 
                        xlevel = logf(powf(2.0,xx) - oldthr + parm[k][3] -
                                      bitl[k] - 0.5);
                    radius = sqrt(parm[k][jj]/CPL_MATH_PI);
                    xdat[ic] = xlevel;
                    xcor[ic++] = radius;
                    dlbydr = (xlevol - xlevel)/MAX(0.01,radius-radold);
                    wt = MIN(1.0,MAX((radius-radold)*5.0,0.1));
                    slope = (1.0 - 0.5*wt)*slope + 0.5*wt*MIN(5.0,dlbydr);
                    radold = radius;
                    xlevol = xlevel;
                }
            }

            /* If this is not the last iteration then work out the effect
               on the local continuum from each image */

            if (! lastone) {
                for (i = 0; i < *nbit; i++) {
                    if (parm[i][0] >= 0.0 && i != k) {
                        dx = parm[k][1] - parm[i][1];
                        dy = parm[k][2] - parm[i][2];
                        dist = sqrtf(dx*dx + dy*dy);
                        xeff = xlevel - MAX(0.0,MIN(50.0,slope*(dist-radius)));
                        bitx[i] += expf(xeff);
                    }
                }

            /* If this is the last iteration loop, then update the parameters
               before exiting*/
            
            } else {
                if (ic > 2) {
                    imcore_polynm(xdat,xcor,ic,polycf,3,0);
                    ttt = polycf[1] + 2.0*polycf[2]*radius;
                } else
                    ttt = 0.0;
                slope = MAX(0.1,MAX(-ttt,slope));
                radthr = radius + (xlevel - algthr)/slope;
                if (radthr > radmax) {
                    slope = 1.0;
                    radthr = radmax;
                }
                
                /* Pixel area */
               
                delb = parm[k][8]*(parm[k][3] - bitl[k]);
                parm[k][8] = CPL_MATH_PI*radthr*radthr;

                /* Peak height */

                parm[k][7] += (parm[k][3] - bitl[k]);
                
                /* Intensity */

                deli = 2.0*CPL_MATH_PI*((parm[k][3] - bitl[k])*(1.0 + slope*radius) -
                                 oldthr*(1.0 + slope*radthr))/(slope*slope);
                parm[k][0] += delb + MAX(0.0,deli);
                for (i = 0; i < 7; i++)
                    parm[k][i+9] = -1.0;
                if (parm[k][0] > xintmn)
                    sumint += parm[k][0];
            }
        }

        /* If this is not the last iteration then check and see how the
           continuum estimates are converging. If they appear to be converging
           then let the next iteration be the last one. */

        if (! lastone) {
            conv = 1;
            for (i = 0; i < *nbit; i++) {
                if (parm[i][0] >= 0.0) {
                    if (abs(bitx[i] - bitl[i]) > 3.0)
                        conv = 0;
                    bitl[i] = bitx[i];
                    bitx[i] = 0;
                    bitl[i] = MIN(bitl[i],NINT(parm[i][3]-oldthr));
                }
            }
            lastone = (conv || (iter == NITER-1));
        } else {
            break;
        }
    }

    /* Find the scaling if needed */

    if (sumint == 0.0) {
        *nbit = 1;
        return;
    } else 
        ratio = total/sumint;
    for (i = 0; i < *nbit; i++)
        parm[i][0] = ratio*parm[i][0];
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        sort_on_zsm_rev
    \par Purpose:
        Sort a Plessey array
    \par Description:
        The entries in a Plessey array are sorted from brightest to faintest
        based on the smoothed flux.
    \par Language:
        C
    \param npts
        The number of pixels in the Plessey array
    \param pts
        The Plessey array
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/
        
static void sort_on_zsm_rev(intptr_t npts, plstruct *pts) {
    intptr_t i,j,ii,jj,ifin;
    plstruct tmp;
/* CONSTANTS: 4,2,4 */
    jj = 4;
    while (jj < npts)
        jj = 2*jj;
    jj = MIN(npts,(3*jj)/4-1);
    while (jj > 1) {
        jj = jj/2;
        ifin = npts - jj;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + jj;
            if (pts[i].zsm > pts[j].zsm)
                continue;
            tmp = pts[j];
            do {
                pts[j] = pts[i];
                j = i;
                i = i - jj;
                if (i < 0) 
                    break;
            } while (pts[i].zsm <= tmp.zsm);
            pts[j] = tmp;
        }
    }
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        moments_thr
    \par Purpose:
        Work out moments for an object with the current threhsold.
    \par Description:
        The moments an object in the current Plessey array are calculated
        for pixels above the current threshold
    \par Language:
        C
    \param ap
        The current ap structure
    \param results
        The moments results
    \param ipk
        The x,y coordinates of the peak pixel
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void moments_thr(ap_t *ap, float results[NPAR+1], intptr_t ipk[2]) {
    intptr_t i,np,nnext;
    float x,y,xoff,yoff,xsum,ysum,xsumsq,ysumsq,tsum,xysum,t,tmax;
    float xbar,ybar,sxx,syy,sxy,fconst,offset,xsum_w,ysum_w,wsum,w;
    plstruct *plarray;

    /* Copy some stuff to local variables */

    fconst = ap->fconst;
    offset = ap->areal_offset;
    plarray = ap->plarray;
    np = ap->npl_pix;

    /* Initialise a few things */

    xoff = xbar_start;
    yoff = ybar_start;
    xsum = 0.0;
    ysum = 0.0;
    xsum_w = 0.0;
    ysum_w = 0.0;
    wsum = 0.0;
    xsumsq = 0.0;
    ysumsq = 0.0;
    tsum = 0.0;
    xysum = 0.0;
    tmax = plarray[0].z - curthr;
    ipk[0] = plarray[0].x;
    ipk[1] = plarray[0].y;
    for (i = 8; i < NPAR; i++)
        results[i] = 0.0;

    /* Do a moments analysis on an object */
    nnext = 0;
    for (i = 0; i < np; i++) {
        x = (float)plarray[i].x - xoff;
        y = (float)plarray[i].y - yoff;
        t = plarray[i].z - curthr;
        w = plarray[i].zsm - curthr;
        if (w > nexthr)
            nnext++;
        xsum += t*x;
        ysum += t*y;
        tsum += t;
        xsum_w += w*t*x;
        ysum_w += w*t*y;
        wsum += w*t;
        xsumsq += (x*x)*t;
        ysumsq += (y*y)*t;
        xysum += x*y*t;
        update_ov(results+8,t,oldthr,fconst,offset);
        if (t > tmax) {
            ipk[0] = plarray[i].x;
            ipk[1] = plarray[i].y;
            tmax = t;
        }
    }

    /* Check that the total intensity is enough and if it is, then do
       the final results. Use negative total counts to signal an error */

    if (tsum > 0.0) {
        results[0] = tsum;
    } else {
        results[0] = -1.0;
        tsum = 1.0;
    }
    xbar = xsum/tsum;
    ybar = ysum/tsum;
    sxx = MAX(0.0,(xsumsq/tsum-xbar*xbar));
    syy = MAX(0.0,(ysumsq/tsum-ybar*ybar));
    sxy = xysum/tsum - xbar*ybar;
    wsum = MAX(1.0,wsum);
    xbar = xsum_w/wsum;
    ybar = ysum_w/wsum;
    xbar += xoff;
    ybar += yoff;
    xbar = MAX(1.0,MIN(xbar,ap->lsiz));
    ybar = MAX(1.0,MIN(ybar,ap->csiz));

    /* Store the results now */

    results[1] = xbar;
    results[2] = ybar;
    results[3] = curthr;
    results[4] = sxx;
    results[5] = sxy;
    results[6] = syy;
    results[7] = tmax;
    results[NPAR] = ((nnext > ap->ipnop && nexthr < lasthr) ? 0 : 1);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        update_ov
    \par Purpose:
        Update areal profiles for a given threshold
    \par Description:
        The areal profiles are updated for a given flux and threshold
    \par Language:
        C
    \param iap
        The areal profile array
    \param t
        The input flux
    \param thresh
        The current detection threshold
    \param fconst
        Scaling parameter (log_2(e))
    \param offset
        The offset between areal levels
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/
        
static void update_ov(float iap[NAREAL], float t, float thresh, float fconst, 
                      float offset) {
    intptr_t nup,i;

    /* Get out of here if the intensity is too small */

    if (t <= 0.0) 
        return;

    /* Otherwise update the relevant profile counts */

    nup = MAX(1,MIN(NAREAL,(intptr_t)(logf(t+thresh)*fconst-offset)+1));
    for (i = 0; i < nup; i++)
        iap[i] += 1.0;
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        check_term
    \par Purpose:
        Check for terminations in ap structure
    \par Description:
        The ap structure is analysed to see if there are any objects that
        can be terminated
    \par Language:
        C
    \param ap
        The current ap structure
    \param nobj
        The output number of objects detected
    \param parm
        The parameter array for all objects
    \param peaks
        The location of the peaks for each of the objects
    \param toomany
        If set, then too many objects have been detected
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void check_term(ap_t *ap, intptr_t *nobj, float parm[IMNUM][NPAR+1], 
                       intptr_t peaks[IMNUM][2], intptr_t *toomany) {
    intptr_t ip,i,ipks[2];
    float momresults[NPAR+1];

    /* Search through all possible parents */

    *nobj = 0;
    *toomany = 0;
    for (ip = 1; ip <= ap->maxip; ip++) {
        if (ap->parent[ip].pnop != -1) {
/*             if (ap->parent[ip].pnop == ap->parent[ip].growing) { */
  
                /* That's a termination: */
 
                if ((ap->parent[ip].pnop >= ap->ipnop &&
                     ap->parent[ip].touch == 0)) {
                    imcore_extract_data(ap,ip);
                    moments_thr(ap,momresults,ipks);
                    if (momresults[0] > 0.0) {
                        if (*nobj == IMNUM-1) {
                            *toomany = 1;
                            break;
                        }
                        for (i = 0; i <= NPAR; i++) 
                            parm[*nobj][i] = momresults[i];
                        for (i = 0; i < 2; i++)
                            peaks[*nobj][i] = ipks[i];
                        (*nobj)++;
                    }
                }
                imcore_restack(ap,ip);
/*          } else { */

/*              /\* This parent still active: *\/ */
 
/*              ap->parent[ip].growing = ap->parent[ip].pnop; */
/*          } */
        }
    }
}

/**@}*/

/*

$Log: imcore_overlp.c,v $
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
