/*
 * This file is part of the CASU Pipeline utilities
 * Copyright (C) 2015,2016 European Southern Observatory
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

#include <stdio.h>
#include <stdlib.h>

#include <cpl_init.h>
#include <cpl_test.h>
//#include <casu_mods.h>
#include "../imcore.h"

#include "../imcore.h"


int main(void) {
    ap_t ap;
    float ell[] = {0.009,0.050,0.025,0.033,0.010,0.012,0.014,0.005,0.011,0.011};
    float pk[] = {1007.065,101.293,204.195,493.250,559.111,609.799,642.603,
                  698.117,740.227,797.222};
    int iareal[][10] = {{120,53,73,104,97,109,104,107,110,110},
                        {100,39,60,75,80,81,86,90,94,95},
                        {75,21,42,63,66,64,68,69,69,70},
                        {63,2,21,45,45,45,45,47,56,56},
                        {45,0,4,25,29,29,37,37,37,37},
                        {25,0,0,9,9,13,14,21,21,21},
                        {9,0,0,0,0,0,0,1,1},
                        {0,0,0,0,0,0,0,0}};

    float *work,fwhm,**areals;
    int i,j;

    /* Initialise */

    cpl_test_init(PACKAGE_BUGREPORT,CPL_MSG_WARNING);

    /* Set up apm structure */

    ap.lsiz = 2048;
    ap.csiz = 2048;
    ap.xtnum = 0;
    ap.thresh = 11.0936;

    /* Areals */

    areals = cpl_malloc(8*sizeof(float *));
    for (i = 0; i < 8; i++) {
        areals[i] = cpl_malloc(10*sizeof(float));
        for (j = 0; j < 10; j++)
            areals[i][j] = (float)iareal[i][j];
    }

    /* Run the test */

    work = cpl_malloc(10*sizeof(float));
    imcore_seeing(&ap, 10, ell, pk, areals, work, &fwhm);
    cpl_free(work);
    for (i = 0; i < 8; i++)
        cpl_free(areals[i]);
    cpl_free(areals);
    cpl_test_rel(fwhm, 4.50384, 0.01);

    /* Get out of here */

    return(cpl_test_end(0));

}
