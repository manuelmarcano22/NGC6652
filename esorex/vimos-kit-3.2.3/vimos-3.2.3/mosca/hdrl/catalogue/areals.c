/* $Id: areals.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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

#include "imcore.h"
#include "util.h"
#include "floatmath.h"

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out the areal profiles for an object
  
    \par Name:
        imcore_areals
    \par Purpose:
        Work out the areal profiles for an object
    \par Description:
        The pixel list for an object is used to define the areal profiles
        for that object and a given detection threshold.
    \par Language:
        C
    \param ap
        The input ap structure
    \param iareal
        The output areal profile array
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

extern void imcore_areals(ap_t *ap, intptr_t iareal[NAREAL]) {
    intptr_t i,nup,j,np;
    float t,thresh,fconst,offset;
    plstruct *plarray;

    /* Initialise some stuff */

    np = ap->npl_pix;
    plarray = ap->plarray;
    thresh = ap->thresh;
    fconst = ap->fconst;
    offset = ap->areal_offset;

    /* Zero the areal profile array */

    (void)memset(iareal,0,NAREAL*sizeof(intptr_t));

    /* Loop through the array and form the areal profiles */

    for (i = 0; i < np; i++) {
        t = plarray[i].z;
        if (t <= thresh) 
            continue;
        nup = MIN(NAREAL,(intptr_t)(logf(t)*fconst - offset)+1);
        nup = MAX(1,nup);
        for (j = 0; j < nup; j++)
            iareal[j]++;
    }
}

/**@}*/
                
/*

$Log: areals.c,v $
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
