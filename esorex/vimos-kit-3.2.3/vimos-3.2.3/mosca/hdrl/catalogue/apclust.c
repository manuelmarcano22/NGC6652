/* $Id: apclust.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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
#include <cpl.h>

#include "imcore.h"
#include "util.h"

static void minmax_xy(intptr_t, plstruct *, intptr_t *, intptr_t *, intptr_t *, intptr_t *);

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Detect multiple objects from a given Plessey array
  
    \par Name:
        imcore_apclust
    \par Purpose:
        Detect multiple peaks in a Plessey array and write to a new ap structure
    \par Description:
        The Plessey array is given from an ap structure for a single object.
        Given a second ap structure with a revised threshold this routine
        will attempt to detect multiple objects within that Plessey array.
    \par Language:
        C
    \param ap
        The new input ap structure
    \param np
        The number of pixels within the input Plessey array
    \param plstr
        The Plessey array from the original structure with the lower
        detection threshold.
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

void imcore_apclust(ap_t *ap, intptr_t np, plstruct *plstr) {

    intptr_t i,i1,loop,ik;
    intptr_t is;    /* parent name for image in this slice */
    intptr_t ip;    /* parent name for image on last line */
    intptr_t ib;    /* data block name */
    intptr_t k,j,ix1,ix2,iy1,iy2,nwork,nx,ny,kk;
    float i2compare,icompare;
    short int *work;
 
    /* A couple of useful things */

    i2compare = ap->thresh;
    icompare = i2compare * ap->multiply;

    /* Get the min and max positions. Create a raster with the IDs of the
       pixels in the pixel list (algorithm prefers data to be in a raster) */

    minmax_xy(np,plstr,&ix1,&ix2,&iy1,&iy2);
    nx = ix2 - ix1 + 1;
    ny = iy2 - iy1 + 1;
    nwork = nx*ny;
    work = cpl_malloc(nwork*sizeof(*work));
    for (i = 0; i < nwork; i++)
        work[i] = -1;
    for (k = 0; k < np; k++) {
        i = plstr[k].x - 1;
        j = plstr[k].y - 1;
        kk = (j - iy1)*nx + i - ix1;
        work[kk] = k;
    }

    /* Now do the job */
            
    for (j = iy1; j <= iy2; j++) {
        for (i = ix1; i <= ix2; i++) {
            kk = (j - iy1)*nx + i - ix1;
            k = work[kk];
            if (k < 0) {
                ap->lastline[i + 1] = 0;
            } else {
                if (plstr[k].zsm > icompare) {

                    /* Pixel is above threshold, find which parent it belongs to. */

                    is = ap->lastline[i];       /* Parent last pixel this line */
                    ip = ap->lastline[i + 1];   /* Guess belongs to above line */
                    if (ip == 0) {

                        /* New parent, or, horizontal slice: */

                        if (is == 0) {

                            /* Ah - new parent. */
/* CONSTANTS: 3/4, 3/8 */
                            if (ap->ipstack > ap->maxpa*3/4) {
                                for (ik = 0; ik < ap->maxpa*3/8; ik++)
                                    imcore_apfu(ap);
                            }
                            ip = ap->pstack[ap->ipstack++];
                            ap->parent[ip].first = ap->bstack[ap->ibstack];
                            ap->parent[ip].pnop = 0;
                            ap->parent[ip].pnbp = 0;
                            ap->parent[ip].growing = 0;
                            if (j == 0)

                                /* It touches first line: */

                                ap->parent[ip].touch = 1;
                            else
                                ap->parent[ip].touch = 0;

                            /* For hunt thru list for terminates: */

                            if (ip > ap->maxip)
                                ap->maxip = ip;
                        } else {

                            /* Slice with no vertical join: */

                            ip = is;
                        }
                    } else if ((ip > 0 && is > 0) && (ip != is)) {

                        /* merge: Join linked lists: */

                        ap->blink[ap->parent[ip].last] = ap->parent[is].first;

                        /* Copy `last block': */

                        ap->parent[ip].last = ap->parent[is].last;
                        ap->parent[ip].pnop += ap->parent[is].pnop;
                        ap->parent[ip].pnbp += ap->parent[is].pnbp;

                        /* Fix `lastline' correlator array: */

                        ib = ap->parent[is].first;
                        loop = 1;
                        while (loop) {
                            i1 = ap->plessey[ib].x;
                            if (ap->lastline[i1 + 1] == is)
                                ap->lastline[i1 + 1] = ip;
                            if (ap->parent[is].last == ib)
                                loop = 0;
                            else
                                ib = ap->blink[ib];
                        }

                        /* Mark parent inactive: */

                        ap->parent[is].pnop = -1;
                        ap->parent[is].pnbp = -1;

                        /* return name to stack: */

                        ap->pstack[--ap->ipstack] = is;
                    }

                    /* Add in pixel to linked list: */

                    ib = ap->bstack[ap->ibstack++];

                    /* Patch forward link into last data block: */

                    if (ap->parent[ip].pnop > 0)
                        ap->blink[ap->parent[ip].last] = ib;

                    /* Remember last block in chain: */

                    ap->parent[ip].last = ib;

                    /* Store the data: */

                    ap->plessey[ib].x = i;
                    ap->plessey[ib].y = j;
                    ap->plessey[ib].z = plstr[k].z;
                    ap->plessey[ib].zsm = plstr[k].zsm;

                    /* increment active count: */

                    ap->parent[ip].pnop++;

                    /* remember which parent this pixel was for next line: */

                    ap->lastline[i + 1] = ip;

                } else {

                    /* Pixel was below threshold, mark lastline: */

                    ap->lastline[i + 1] = 0;
                }
            }
        }
    }

    /* Check for images touching left & right edges:
       OR the touch flag with 2 for left, 4 for right: */
/* CONSTANTS 2, 4 */
    if(ap->lastline[1] > 0 )
        ap->parent[ap->lastline[1]].touch |= 2;
    if(ap->lastline[ap->lsiz] > 0)
        ap->parent[ap->lastline[ap->lsiz]].touch |= 4;
    cpl_free(work);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        minmax_xy
    \par Purpose:
        Work out the min and max x,y positions within a Plessey array
    \par Description:
        All the pixels within the given Plessey array are searched and the
        min and max positions are returned. These values have 1 subtracted 
        off so they can be used as array indices.
    \par Language:
        C
    \param np
        The number of pixels in the Plessey array
    \param plstr
        The Plessey array
    \param ix1, ix2, iy1, iy2
        The min and max x,y coordinates.
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void minmax_xy(intptr_t np, plstruct *plstr, intptr_t *ix1, intptr_t *ix2,
                      intptr_t *iy1,
                      intptr_t *iy2) {
    intptr_t i;

    /* Get the minmax of the positions of the pixels in a plstruct.  Take
       1 away from each position so that it runs from 0 rather than 1 */

    *ix1 = plstr[0].x - 1;
    *ix2 = plstr[0].x - 1;
    *iy1 = plstr[0].y - 1;
    *iy2 = plstr[0].y - 1; 
    for (i = 1; i < np; i++) {
        *ix1 = MIN(*ix1,plstr[i].x - 1);
        *ix2 = MAX(*ix2,plstr[i].x - 1);
        *iy1 = MIN(*iy1,plstr[i].y - 1);
        *iy2 = MAX(*iy2,plstr[i].y - 1);
    }
}

/**@}*/

/*

$Log: apclust.c,v $
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
