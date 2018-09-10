/* $Id: terminate.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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
#include "imcore.h"
#include "util.h"

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Free information for an object from the ap structure
  
    \par Name:
        imcore_restack
    \par Purpose:
        Free information for an object from the ap structure
    \par Description:
        The starting address for an object in the ap structure is given.
        Information relating to that object is erased and the space made
        available
    \par Language:
        C
    \param ap
        The current ap structure
    \param ip
        The parent number for the object
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

extern void imcore_restack(ap_t *ap, intptr_t ip) {
    intptr_t i,ib,nn,np;
    unsigned char *mflag;

    /* Reset the mflag */

    np = ap->parent[ip].pnop;
    ib = ap->parent[ip].first;
    mflag = ap->mflag;
    for (i = 0; i < np; i++) {
        nn = ap->plessey[ib].y*ap->lsiz + ap->plessey[ib].x;
        mflag[nn] = MF_POSSIBLEOBJ;
        ib = ap->blink[ib];
    }

    /* Stash all blocks back in a burst: */
    
    ib = ap->parent[ip].first;
    for(i = ap->ibstack - ap->parent[ip].pnop; i < ap->ibstack-1;  i++) {
        ap->bstack[i] = ib;
        ib = ap->blink[ib];
    }

    /* and the last one: */

    ap->bstack[ap->ibstack-1] = ib;
    ap->ibstack -= ap->parent[ip].pnop;

    /* Put parent name back on stack: */

    ap->pstack[--ap->ipstack] = ip;

    /* Mark that parent inactive: */

    ap->parent[ip].pnop = -1;
    ap->parent[ip].pnbp = -1;
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Check for objects that have terminated
  
    \par Name:
        imcore_terminate
    \par Purpose:
        Check for objects that have terminated
    \par Description:
        The parents in the current ap structure are examined to see which
        have not grown since the last pass. Any that have not grown are
        sent to the processing routine.
    \par Language:
        C
    \param ap
        The current ap structure
    \param cattype
        The type of catalogue to be produced
    \param gain
        The header keyword with the gain in e-/ADU
    \param nobjects
        Number of detected objects
    \param tab
        Output catalogue table
    \param res
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

extern void imcore_terminate(ap_t *ap, int cattype, float gain,
                             intptr_t *nobjects,
                             cpl_table *tab, hdrl_imcore_result * res) {
    intptr_t ip;
    int status;

    /* Search through all possible parents!  */

    for (ip = 1; ip <= ap->maxip; ip++) {
        if(ap->parent[ip].pnop != -1) {
            if(ap->parent[ip].pnop == ap->parent[ip].growing) {

                /* That's a termination: */

                if((ap->parent[ip].pnop >= ap->ipnop &&
                    ap->parent[ip].touch == 0) &&
                    (ap->parent[ip].pnbp < (ap->parent[ip].pnop)/2)) {
                    imcore_extract_data(ap,ip);
                    
                    /* Call the processing routine */

                    status = imcore_process_results(ap,cattype,gain,nobjects,
                                                    tab, res);
                    if (status != CASU_OK) {
                        imcore_restack(ap,ip);
                        continue;
                    }
                }
                imcore_restack(ap,ip);
            } else {

                /* This parent still active: */

                ap->parent[ip].growing = ap->parent[ip].pnop;
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Get rid of the largest contributor in an ap structure
  
    \par Name:
        imcore_apfu
    \par Purpose:
        Get rid of the largest contributor in an ap structure
    \par Description:
        The parents in the current ap structure are examined to see which
        has the largest number of pixels. That parent is junked.
    \par Language:
        C
    \param ap
        The current ap structure
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

extern void imcore_apfu(ap_t *ap) {
    intptr_t ip, big, ipbig;

    /* Search through all possible parents and just junk the biggest
       one to free space:  */

    big = 0;
    ipbig = 0;
    for (ip = 1; ip <= ap->maxip; ip++) {
        if(ap->parent[ip].pnop != -1) {
            if(ap->parent[ip].pnop > big) {
                big = ap->parent[ip].pnop;
                ipbig = ip;
            }
        }
    }
    if(big > 0) {
        imcore_restack(ap, ipbig);

        /* clearout lastline references to this parent: */

        for (ip = 0; ip <= ap->lsiz; ip++)
            if(ap->lastline[ip] == ipbig) ap->lastline[ip] = 0;
    }
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Put data into the Plessey array for an object
  
    \par Name:
        imcore_extract_data
    \par Purpose:
        Put data into the Plessey array for an object
    \par Description:
        The information for the object from a given parent is extracted from
        the link list in the ap structure and put into the Plessey array in
        preparation for analysis.
    \par Language:
        C
    \param ap
        The current ap structure
    \param ip
        The parent in question
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

extern void imcore_extract_data(ap_t *ap, intptr_t ip) {
    intptr_t ib,i,np,nn;
    unsigned char *mflag;

    /* Check the size of the workspace and see if it's big enough. If it
       isn't then increase the size until it is */

    np = ap->parent[ip].pnop;
    if (ap->npl < np) {
        ap->plarray = cpl_realloc(ap->plarray,np*sizeof(plstruct));
        ap->npl = np;
    }

    /* Pull the info out now */

    ib = ap->parent[ip].first;
    ap->npl_pix = np;
    mflag = ap->mflag;
    for (i = 0; i < np; i++) {
        ap->plarray[i].x = ap->plessey[ib].x + 1;
        ap->plarray[i].y = ap->plessey[ib].y + 1;
        ap->plarray[i].z = ap->plessey[ib].z;
        ap->plarray[i].zsm = ap->plessey[ib].zsm;
        nn = ap->plessey[ib].y*ap->lsiz + ap->plessey[ib].x;
        mflag[nn] = MF_OBJPIX;
        ib = ap->blink[ib];
    }
}

/**@}*/

/*

$Log: terminate.c,v $
Revision 1.3  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.2  2015/08/07 13:06:54  jim
Fixed copyright to ESO

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.3  2015/01/09 11:42:36  jim
Fixed routines to remove globals

Revision 1.2  2014/04/09 09:09:51  jim
Detabbed

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
