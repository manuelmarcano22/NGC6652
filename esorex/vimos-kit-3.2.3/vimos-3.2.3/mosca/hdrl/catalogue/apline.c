/* $Id: apline.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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
 * $Date: 
 * $Revision: 1.3 $
 * $Name:  $
 */

#include <stdio.h>
#include "imcore.h"
#include "util.h"

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Detect objects on a line of data
  
    \par Name:
        imcore_apline
    \par Purpose:
        Detect objects on a line of data
    \par Description:
        Pixels above the defined threshold are detected on a line of data. 
        A search is done in the ap structure to see if the pixel forms part
        of an object that is already known or whether this is part of a new
        object. The information is stored in the ap structure linked list.
    \par Language:
        C
    \param ap
        The input ap structure
    \param dat
        The line of data
    \param conf
        The confidence map for the line of data
    \param smoothed
        A smoothed version of the data line.
    \param smoothedc
        A smoothed version of the confidence map line.
    \param j
        A number that tells you which row this is in the image
    \param bpm
        A bad pixel mask for that line
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

void imcore_apline(ap_t *ap, float dat[], float conf[], float smoothed[], 
                   float smoothedc[], intptr_t j, unsigned char *bpm) {
    intptr_t i,i1,loop,nn;
    intptr_t is;    /* parent name for image in this slice */
    intptr_t ip;    /* parent name for image on last line */
    intptr_t ib;    /* data block name */
    float i2compare,icompare;
    unsigned char *mflag;

    i2compare = ap->thresh;
    icompare = i2compare * ap->multiply;
    mflag = ap->mflag;

    for (i = 0; i < ap->lsiz; i++) {
        if (smoothedc[i] > icompare && conf[i] != 0) {

            /* Pixel is above threshold, find which parent it belongs to. */

            is = ap->lastline[i];       /* Parent last pixel this line */
            ip = ap->lastline[i + 1];   /* Guess belongs to above line */
            if (ip == 0) {

                /* New parent, or, horizontal slice: */

                if (is == 0) {

                    /* Ah - new parent. */

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
            ap->plessey[ib].z = dat[i];
            nn = j*ap->lsiz + i;
            if (mflag[nn] != MF_SATURATED)
                ap->plessey[ib].zsm = MIN(ap->saturation,smoothed[i]);
            else
                ap->plessey[ib].zsm = ap->saturation;
            mflag[nn] = MF_POSSIBLEOBJ;

            /* increment active count: */

            ap->parent[ip].pnop++;
            if (bpm != NULL) 
                ap->parent[ip].pnbp += bpm[i];

            /* remember which parent this pixel was for next line: */

            ap->lastline[i + 1] = ip;

        } else {

            /* Pixel was below threshold, mark lastline: */

            ap->lastline[i + 1] = 0;
        }
    }

    /* Check for images touching left & right edges:
       OR the touch flag with 2 for left, 4 for right: */
/* CONSTANTS: 2,4 */
    if(ap->lastline[1] > 0 )
        ap->parent[ap->lastline[1]].touch |= 2;
    if(ap->lastline[ap->lsiz] > 0)
        ap->parent[ap->lastline[ap->lsiz]].touch |= 4;
}

/**@}*/

/*

$Log: apline.c,v $
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
