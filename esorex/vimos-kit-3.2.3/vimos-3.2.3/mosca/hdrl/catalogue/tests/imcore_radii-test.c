/* $Id: imcore_radii-test.c,v 1.1 2015/10/15 11:27:22 jim Exp $
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
 * $Date: 2015/10/15 11:27:22 $
 * $Revision: 1.1 $
 * $Name:  $
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpl_init.h>
#include <cpl_test.h>
//#include <casu_tfits.h>
//#include <casu_utils.h>
//#include <casu_mods.h>
#include "../imcore.h"
#include "../imcore_radii.h"
#include "../casu_utilfunctions.h"

#define NRADII 13

int main(void) {
    float radii[] = {2.5,3.53553,5.0,7.07107,10,14,20,25,30,35,40,50,60};
    float fluxes[] = {13670.3,19834.4,23923.2,25124.0,25332.3,25488.9,
                      25648.7,25842.8,25950.9,25893.9,25982,1,25297.6,
                      24919.1};
    float halfrad = 2.35;
    float kronrad = 6.18;
    float petrrad = 12.45;
    float halflight,rad;
    float peak = 1007.07;
    float areal = 120.0;

    /* Initialise */

    cpl_test_init(PACKAGE_BUGREPORT,CPL_MSG_WARNING);

    /* Test halflight */

    halflight = fluxes[4]/2.0;
    rad = imcore_halflight(radii,fluxes,halflight,peak,NRADII);
    cpl_test_rel(rad,halfrad,0.01);

    /* Test Kron */

    rad = imcore_kronrad(areal,radii,fluxes,NRADII);
    cpl_test_rel(rad,kronrad,0.01);

    /* Test Petrosian */

    rad = imcore_petrad(areal,radii,fluxes,NRADII);
    cpl_test_rel(rad,petrrad,0.01);

    /* Get out of here */


    return(cpl_test_end(0));

}

/*

$Log: imcore_radii-test.c,v $
Revision 1.1  2015/10/15 11:27:22  jim
new


*/
