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

#include "../imcore.h"
#include <cpl.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    ap_t ap;
    float skymed,skysig;
    int retval,i,j;

    /* Initialise */

    cpl_test_init(PACKAGE_BUGREPORT,CPL_MSG_WARNING);

    /* Create an input apm structure */
    
    ap.lsiz = 2048;
    ap.csiz = 2048;
    ap.xtnum = 0;
    ap.inframe = cpl_image_new(2048,2048,CPL_TYPE_FLOAT);
    ap.conframe = cpl_image_new(2048,2048,CPL_TYPE_FLOAT);
    imcore_apinit(&ap);
    ap.indata = cpl_image_get_data_float(ap.inframe);
    ap.confdata = cpl_image_get_data_float(ap.conframe);
    ap.mflag = cpl_calloc(2048*2048,sizeof(unsigned char));

    /* Create a background */

    cpl_image_fill_noise_uniform(ap.inframe,-10.0,10.0);
    cpl_image_add_scalar(ap.inframe,5000.0);
    cpl_image_fill_noise_uniform(ap.conframe,99,101);

    /* Get the background value */

    retval = imcore_backstats(&ap,-100.0,&skymed,&skysig);
    cpl_test_eq(retval,CASU_OK);
    cpl_test_rel(skymed,5000.0,0.01);
    cpl_test_rel(skysig, 20 / sqrt(12), 0.1);

    /* Create a background map */

    hdrl_imcore_result * res = cpl_calloc(sizeof(*res), 1);
    res->background = cpl_image_new(ap.lsiz, ap.csiz, CPL_TYPE_FLOAT);
    retval = imcore_background(&ap,64,-100.0, 1, res);

    cpl_test_eq(retval,CASU_OK);
    for (j = 0; j < ap.backmap.nby; j++) 
        for (i = 0; i < ap.backmap.nbx; i++) 
            cpl_test_rel((ap.backmap.bvals)[i][j],5000.0,0.01);
    cpl_test_rel(cpl_image_get_median(ap.inframe),5000.0,0.1);
    imcore_backest(&ap,1000.0,1000.0,&skymed,&skysig);
    cpl_test_rel(skymed,5000.0,0.01);
    cpl_test_lt(0.,skysig);

    //cpl_table_delete(res->catalogue);
    //cpl_image_delete(res->segmentation_map);
    cpl_image_delete(res->background);
    cpl_free(res); // TODO check results
    /* Get out of here */

    imcore_apclose(&ap);
    cpl_free(ap.mflag);
    cpl_image_delete(ap.inframe);
    cpl_image_delete(ap.conframe);
    return(cpl_test_end(0));

}
