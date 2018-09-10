/* $Id: hdrl_imagelist_basic-test.c,v 1.6 2013-10-02 12:49:29 yjung Exp $
 *
 * This file is part of the HDRL
 * Copyright (C) 2013 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 * $Author: yjung $
 * $Date: 2013-10-02 12:49:29 $
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                    Includes
 -----------------------------------------------------------------------------*/

#include "hdrl_imagelist.h"
#include "hdrl_test.h"

#include <cpl.h>

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef ARRAY_LEN
#define ARRAY_LEN(a) sizeof((a))/sizeof((a)[0])
#endif
#define    IMAGESZ     265
#define    IMAGENB     10

/*----------------------------------------------------------------------------*/
/**
 * @defgroup hdrl_imagelist_basic_test   
            Testing of hdrl_imagelist_basic module
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 @brief   Unit tests of hdrl_image_basic
 **/
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    hdrl_imagelist  *   himlist ;
    hdrl_image      *   himg ;
    hdrl_image      *   himg2 ;
    cpl_image       *   contrib ;
    cpl_image       *   contrib2 ;
    cpl_size            nimages = IMAGENB ;

    /* Create an image list */
    himlist = hdrl_imagelist_new() ;
    for (cpl_size i = 0 ; i < nimages ; i++) {
        /*
        What about having an

        hdrl_image_load_image(hdrl_image * hdrlimage
                           const char * filename,
                           cpl_type     im_type,
                           cpl_size     pnum,
                           cpl_size     xtnum);

        hdrl_image_load_error(hdrl_image * hdrlimage
                           const char * filename,
                           cpl_type     im_type,
                           cpl_size     pnum,
                           cpl_size     xtnum);

        with the same API like cpl?
        two functions are better in the case the error is in another fitsfile
        */

        himg = hdrl_image_new(IMAGESZ, IMAGESZ);
        if (i == nimages / 2)
            hdrl_image_add_scalar(himg, (hdrl_value){1000, 100});
        else
            hdrl_image_add_scalar(himg, (hdrl_value){i, 1});
        hdrl_imagelist_set(himlist, himg, i);
    }
   
    // Mean collapse
    hdrl_imagelist_collapse_mean(himlist, &himg, &contrib);
    hdrl_imagelist_collapse(himlist, HDRL_COLLAPSE_MEAN, &himg2, &contrib2);
    hdrl_test_image_abs(himg, himg2, 0);
    cpl_test_image_abs(contrib, contrib2, 0);
    cpl_image_delete(contrib) ;
    hdrl_image_delete(himg) ;
    cpl_image_delete(contrib2) ;
    hdrl_image_delete(himg2) ;

    // Weighted Mean collapse
    hdrl_imagelist_collapse_weighted_mean(himlist, &himg, &contrib);
    hdrl_imagelist_collapse(himlist, HDRL_COLLAPSE_WEIGHTED_MEAN, &himg2, &contrib2);
    hdrl_test_image_abs(himg, himg2, 0);
    cpl_test_image_abs(contrib, contrib2, 0);
    cpl_image_delete(contrib) ;
    hdrl_image_delete(himg) ;
    cpl_image_delete(contrib2) ;
    hdrl_image_delete(himg2) ;

    // Median collapse
    hdrl_imagelist_collapse_median(himlist, &himg, &contrib);
    hdrl_imagelist_collapse(himlist, HDRL_COLLAPSE_MEDIAN, &himg2, &contrib2);
    hdrl_test_image_abs(himg, himg2, 0);
    cpl_test_image_abs(contrib, contrib2, 0);
    cpl_image_delete(contrib) ;
    hdrl_image_delete(himg) ;
    cpl_image_delete(contrib2) ;
    hdrl_image_delete(himg2) ;

    // Sigclip collapse
    cpl_image * rej_low, * rej_high;
    hdrl_imagelist_collapse_sigclip(himlist, 1.0, 3.0, 10, &himg, &contrib,
                                    NULL, NULL);
    cpl_image_delete(contrib) ;
    hdrl_image_delete(himg) ;
    hdrl_imagelist_collapse_sigclip(himlist, 1.0, 3.0, 10, &himg, &contrib,
                                    &rej_low, NULL);
    cpl_image_delete(contrib) ;
    cpl_image_delete(rej_low) ;
    hdrl_image_delete(himg) ;
    hdrl_imagelist_collapse_sigclip(himlist, 1.0, 3.0, 10, &himg, &contrib,
                                    NULL, &rej_high);
    cpl_image_delete(contrib) ;
    cpl_image_delete(rej_high) ;
    hdrl_image_delete(himg) ;
    hdrl_imagelist_collapse_sigclip(himlist, 1.0, 3.0, 10, &himg, &contrib,
                                    &rej_low, &rej_high);
    cpl_image_delete(contrib) ;
    cpl_image_delete(rej_low) ;
    cpl_image_delete(rej_high) ;
    hdrl_image_delete(himg) ;

    hdrl_imagelist_delete(himlist) ;
    return cpl_test_end(0);
}

