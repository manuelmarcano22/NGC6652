/* $Id: apinit.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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


#include <stdlib.h>
#include <cpl.h>
#include "imcore.h"

#define freespace(_p) if (_p != NULL) {cpl_free(_p); _p = NULL;}

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Initialise the ap structure
  
    \par Name:
        imcore_apinit
    \par Purpose:
        Initialise the ap structure given some pre-existing information
    \par Description:
        The ap structure is initialised. In order for this to be done
        properly the value of ap->lsiz (the length of the image rows)
        must be set before calling this routine.
    \par Language:
        C
    \param ap
        The input ap structure
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

void imcore_apinit(ap_t *ap) {
    intptr_t maxpa;

    intptr_t i;
/* CONSTANTS 2, 8, */
    maxpa = ap->lsiz / 2;               /* max possible parents */
    ap->lastline = cpl_calloc(ap->lsiz + 1, sizeof(short int));
    ap->maxip = 0;
    ap->maxpa = maxpa;
    ap->pstack = cpl_malloc(maxpa*sizeof(*ap->pstack));
    ap->parent = cpl_malloc(maxpa*sizeof(*(ap->parent)));
    for(i = 0; i < maxpa; i++) {
        ap->pstack[i] = i;
        ap->parent[i].pnop = -1;        /* mark all parents inactive */
        ap->parent[i].pnbp = -1;        /* mark all parents inactive */
    }
    ap->ipstack = 1;
    ap->maxbl = MAXBL;
    ap->bstack = cpl_malloc(ap->maxbl*sizeof(*ap->bstack));
    ap->blink = cpl_malloc(ap->maxbl*sizeof(*ap->blink));
    ap->plessey = cpl_malloc(ap->maxbl*sizeof(*ap->plessey));
    for (i = 0; i < MAXBL; i++)
        ap->bstack[i] = i;
    ap->ibstack = 2;    /* block 1 will get overwritten; don't use it */
    ap->nimages = 0;

    /* set up exponential areal-profile levels: */

    ap->areal[0] = 1;
    for (i = 1; i < 8; i++)
        ap->areal[i] = ap->areal[i-1]*2;

    /* allocate some space for a processing array */

    ap->npl = ap->lsiz;
    ap->npl_pix = 0;
    ap->plarray = cpl_malloc(ap->npl*sizeof(plstruct));

    /* set these to null values as you may not need the background structure */

    ap->backmap.nby = -1;
    ap->backmap.bvals = NULL;

    /* Initialise some info about the input images */

    ap->indata = NULL;
    ap->confdata = NULL;
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Re-initialise the ap structure
  
    \par Name:
        imcore_apreinit
    \par Purpose:
        Re-initialise the ap structre
    \par Description:
        The ap structure is reinitialised to the state it was in before
        we started detecting objects in it. All information about detected
        objects is erased.
    \par Language:
        C
    \param ap
        The input ap structure
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

void imcore_apreinit(ap_t *ap) {
    intptr_t i;
/* CONSTANTS 2 */
    for (i = 0; i < ap->lsiz+1; i++)
        ap->lastline[i] = 0;
    ap->maxip = 0;
    for(i = 0; i < ap->maxpa; i++) {
        ap->pstack[i] = i;
        ap->parent[i].pnop = -1;        /* mark all parents inactive */
        ap->parent[i].pnbp = -1;        /* mark all parents inactive */
    }
    ap->ipstack = 1;
    ap->ibstack = 2;    /* block 1 will get overwritten; don't use it */
    ap->nimages = 0;
    ap->npl_pix = 0;

}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Close ap structure
  
    \par Name:
        imcore_apclose
    \par Purpose:
        Close the ap structre
    \par Description:
        The ap structure has all of its allocated memory freed.
    \par Language:
        C
    \param ap
        The input ap structure
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

void imcore_apclose(ap_t *ap) {
    intptr_t i;
    freespace(ap->lastline);
    freespace(ap->pstack);
    freespace(ap->parent);
    freespace(ap->bstack);
    freespace(ap->blink);
    freespace(ap->plessey);
    freespace(ap->plarray);
    if (ap->backmap.bvals != NULL) {
        for (i = 0; i < ap->backmap.nby; i++)
            freespace(ap->backmap.bvals[i]);
        freespace(ap->backmap.bvals);
    }
}

/*

$Log: apinit.c,v $
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
