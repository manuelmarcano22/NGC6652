/* $Id: casu_imcore-test.c,v 1.3 2015/10/15 11:18:32 jim Exp $
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
 * $Date: 2015/10/15 11:18:32 $
 * $Revision: 1.3 $
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
#include "../casu_utilfunctions.h"
#include "../../hdrl_catalogue.h"

#include <math.h>

#define NTEST 10
#define NCOLS 63

int main(void) {
    intptr_t retval,i,nrow,ncol;
    int status, nl;

    cpl_image *im,*bkg,*cnf;
    cpl_table *tab;
    double sigma=2.0, norm2, tot[NTEST], sky=500.0, diff;
    double xpos[] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,
                     1000.0};
    double ypos[] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,
                     1000.0};
    double norm[] = {1000.0,100.0,200.0,500.0,550.0,600.0,650.0,700.0,
                     750.0,800.0};
    casu_fits *inf,*inconf;
    cpl_propertylist *pl;

    /* Initialise */

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Check inherited status */

    status = CASU_FATAL;
    hdrl_imcore_result * res = cpl_calloc(sizeof(*res), 1);
    retval = hdrl_casu_imcore(NULL, NULL, NULL, 5, 1.5, 1, 3.5, 1, 64, 6, 1.0, 2.0,
                         HDRL_SATURATION_INIT, res, &status);
             /*casu_imcore(casu_fits *infile, casu_fits *conf, const cpl_wcs * wcs,
                       int ipix,
                       float threshold, int icrowd, float rcore, int nbsize,
                       int cattype, float filtfwhm, float gainloc,
                       float saturation, hdrl_imcore_result * res, int *status) */
                    
    cpl_test_eq(status,CASU_FATAL);
    cpl_test_eq(status,retval);
    cpl_test_null(res->catalogue);

    cpl_image_delete(res->segmentation_map); // TODO check results
    cpl_image_delete(res->background);
    casu_tfits_delete(res->catalogue);
    cpl_free(res);
    

    /* Generate a field with some stars and a confidence map */

    bkg = cpl_image_new(1024,1024,CPL_TYPE_FLOAT);
    im = cpl_image_new(1024,1024,CPL_TYPE_FLOAT);
    cnf = cpl_image_new(1024,1024,CPL_TYPE_FLOAT);
    norm2 = 2.0*CPL_MATH_PI*sigma*sigma;
    cpl_image_fill_noise_uniform(bkg,-10.0,10.0);
    cpl_image_add_scalar(bkg,sky);
    cpl_image_fill_noise_uniform(cnf,99.9,100.1);
    for (i = 0; i < NTEST; i++) {
        cpl_image_fill_gaussian(im,xpos[i],ypos[i],norm[i]*norm2,sigma,sigma);
        tot[i] = cpl_image_get_flux(im);
        cpl_image_add(bkg,im);
    }
    pl = cpl_propertylist_new();
    inf = casu_fits_wrap(bkg,NULL,NULL,pl);
    inconf = casu_fits_wrap(cnf,NULL,NULL,NULL);
    cpl_image_delete(im);
    cpl_propertylist_delete(pl);

    /* Give it a WCS */

    pl = casu_fits_get_ehu(inf);
    cpl_propertylist_update_string(pl,"CTYPE1","RA---TAN");
    cpl_propertylist_update_string(pl,"CTYPE2","DEC--TAN");
    cpl_propertylist_update_double(pl,"CRVAL1",30.0);
    cpl_propertylist_update_double(pl,"CRVAL2",12.0);
    cpl_propertylist_update_double(pl,"CRPIX1",512.0);
    cpl_propertylist_update_double(pl,"CRPIX2",512.0);
    cpl_propertylist_update_double(pl,"CD1_1",-1.0/3600);
    cpl_propertylist_update_double(pl,"CD1_2",0.0);
    cpl_propertylist_update_double(pl,"CD2_1",0.0);
    cpl_propertylist_update_double(pl,"CD2_2",1.0/3600);
    cpl_propertylist_update_int(pl, "NAXIS1", 1024);
    cpl_propertylist_update_int(pl, "NAXIS2", 1024);
    cpl_wcs *wcs = cpl_wcs_new_from_propertylist(pl);

    /* Run imcore */

    status = CASU_OK;
    res = cpl_calloc(sizeof(*res), 1);
    retval = hdrl_casu_imcore(inf, inconf, wcs, 5, 1.5, 0, 5.0, 1, 64, 6, 3.0, 
                         1.0, HDRL_SATURATION_INIT, res, &status);
    cpl_test_eq(status,CASU_OK);
    cpl_test_eq(status,retval);
    cpl_test_nonnull(res->catalogue);
    casu_fits_delete(inf);
    casu_fits_delete(inconf);
    
    cpl_wcs_delete(wcs);
    cpl_image_delete(res->segmentation_map);
    cpl_image_delete(res->background);
    
    /* Check the results. Start by checking the number of rows and columns. 
       Sort the table by X */
    tab = casu_tfits_get_table(res->catalogue);
    cpl_test_nonnull(tab);
    ncol = cpl_table_get_ncol(tab);
    cpl_test_eq(ncol, NCOLS);
    nrow = cpl_table_get_nrow(tab);
    cpl_test_eq(nrow, NTEST);
    pl = cpl_propertylist_new();
    cpl_propertylist_append_bool(pl, "X_coordinate", 0);
    cpl_table_sort(tab, pl);
    cpl_propertylist_delete(pl);

    /* Test the column content of the table */
    for (i = 0; i < NTEST; i++) {
        cpl_test_abs(xpos[i],  cpl_table_get_float(tab, "X_coordinate", (cpl_size)i, &nl),  0.2);
        cpl_test_abs(ypos[i],  cpl_table_get_float(tab, "Y_coordinate", (cpl_size)i, &nl),  0.2);
        diff = fabs( cpl_table_get_float(tab, "Aper_flux_5", (cpl_size)i, &nl) - tot[i]);
        diff /= cpl_table_get_float(tab, "Aper_flux_5_err", (cpl_size)i, &nl);
        cpl_test_lt(diff, 1.6);
        cpl_test_eq( cpl_table_get_float(tab, "Classification", (cpl_size)i, &nl),  -1.0);
    }

    /* Compare some header info */
    pl = casu_tfits_get_ehu(res->catalogue);
    cpl_test_rel( cpl_propertylist_get_float(pl, "ESO QC IMAGE_SIZE")/2.355, sigma, 0.1);

    /* Get out of here */
    casu_tfits_delete(res->catalogue);
    cpl_free(res);
    return( cpl_test_end(0) );
}

/*

$Log: casu_imcore-test.c,v $
Revision 1.3  2015/10/15 11:18:32  jim
More comprehensive

Revision 1.2  2015/08/07 13:06:54  jim
Fixed copyright to ESO

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.1  2015/01/09 11:39:55  jim
new entry


*/
