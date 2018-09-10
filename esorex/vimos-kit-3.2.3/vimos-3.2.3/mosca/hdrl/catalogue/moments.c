/* $Id: moments.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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

#include "imcore.h"
#include "util.h"

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do moments analysis on an object
  
    \par Name:
        imcore_moments
    \par Purpose:
        Do moments analysis on an object in a Plessey array
    \par Description:
        The Plessey array is given from an ap structure for a single object.
        This routine does a basic zeroth, first and second moments analysis.
    \par Language:
        C
    \param ap
        The current ap structure
    \param results
        The output array with the moments results.
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

extern void imcore_moments(ap_t *ap, float results[]) {
    intptr_t i,np;
    float x,y,xoff,yoff,xsum,ysum,xsumsq,ysumsq,tsum,xysum,t,tmax;
    float xbar,ybar,sxx,syy,sxy,xintmin,w,wsum,xsum_w,ysum_w;
    plstruct *plarray;

    /* Initialise a few things */

    xintmin = ap->xintmin;
    plarray = ap->plarray;
    np = ap->npl_pix;
    xoff = (float)plarray[0].x;
    yoff = (float)plarray[0].y;
    xsum = 0.0;
    ysum = 0.0;
    xsum_w = 0.0;
    ysum_w = 0.0;
    wsum = 0.0;
    xsumsq = 0.0;
    ysumsq = 0.0;
    tsum = 0.0;
    xysum = 0.0;
    tmax = plarray[0].z;
    /* Do a moments analysis on an object */

    for (i = 0; i < np; i++) {
        x = (float)plarray[i].x - xoff;
        y = (float)plarray[i].y - yoff;
        t = plarray[i].z;
        w = plarray[i].zsm; 
        if (t < 0.0)
            continue;
        xsum += t*x;
        ysum += t*y;
        tsum += t;
        xsum_w += w*t*x;
        ysum_w += w*t*y;
        wsum += w*t;
        tmax = MAX(tmax,plarray[i].z);
        xsumsq += (x*x)*t;
        ysumsq += (y*y)*t;
        xysum += x*y*t;
    }

    /* Check that the total intensity is enough and if it is, then do
       the final results */

    if (tsum >= xintmin) {
        xbar = xsum/tsum;
        ybar = ysum/tsum;
        sxx = MAX(0.0,(xsumsq/tsum-xbar*xbar));
        syy = MAX(0.0,(ysumsq/tsum-ybar*ybar));
        sxy = xysum/tsum - xbar*ybar;
        xbar = xsum_w/wsum;
        ybar = ysum_w/wsum;
        xbar += xoff;
        ybar += yoff;
        xbar = MAX(1.0,MIN(xbar,ap->lsiz));
        ybar = MAX(1.0,MIN(ybar,ap->csiz));
        results[0] = 1.0;
        results[1] = xbar;
        results[2] = ybar;
        results[3] = tsum;
        results[4] = sxx;
        results[5] = sxy;
        results[6] = syy;
        results[7] = tmax;
    } else {
        results[0] = -1.0;
    }
}

/**@}*/

/*

$Log: moments.c,v $
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
