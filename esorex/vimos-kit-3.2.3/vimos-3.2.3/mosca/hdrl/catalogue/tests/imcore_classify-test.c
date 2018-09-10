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


#include "../casu_utilfunctions.h"
#include "../imcore.h"
#include <hdrl.h>
#include <cpl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NTEST 10
#define NCOLS 63

int main(void) {
    int retval,i,nrow,ncol,nl;
    cpl_image *im,*bkg,*cnf;
    cpl_table *tab;
    double sigma=2.0,norm2,tot[NTEST],sky=500.0,diff,val;
    double xpos[] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,
                     1000.0};
    double ypos[] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,
                     1000.0};
    double norm[] = {1000.0,100.0,200.0,500.0,550.0,600.0,650.0,700.0,
                     750.0,800.0};
    casu_fits *inf,*inconf;
    cpl_propertylist *pl;

    /* Initialise */

    cpl_test_init(PACKAGE_BUGREPORT,CPL_MSG_WARNING);

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

    /* Run imcore */

    hdrl_imcore_result * res = cpl_calloc(sizeof(*res), 1);
    retval = imcore_conf(inf, inconf, 5, 1.5, 0, 5.0, 1, 64, 6, 3.0, 1.0,
                    HDRL_SATURATION_INIT, res);
    cpl_test_eq(retval,CASU_OK);
    cpl_test_nonnull(res->catalogue);

    cpl_image_delete(res->segmentation_map); // TODO check results
    cpl_image_delete(res->background);

    /* Check the results. Start by checking the number of rows and columns.
       Sort the table by X */

    tab = casu_tfits_get_table(res->catalogue);
    cpl_test_nonnull(tab);
    ncol = cpl_table_get_ncol(tab);
    cpl_test_eq(ncol, NCOLS);
    nrow = cpl_table_get_nrow(tab);
    cpl_test_eq(nrow,NTEST);
    pl = cpl_propertylist_new();
    cpl_propertylist_append_bool(pl,"X_coordinate",0);
    cpl_table_sort(tab,pl);
    cpl_propertylist_delete(pl);

    /* Test the column content of the table */

    for (i = 0; i < NTEST; i++) {
        cpl_test_abs(xpos[i],cpl_table_get_float(tab,"X_coordinate",(cpl_size)i,
                                                 &nl),0.2);
        cpl_test_abs(ypos[i],cpl_table_get_float(tab,"Y_coordinate",(cpl_size)i,
                                                 &nl),0.2);
        diff = fabs(cpl_table_get_float(tab,"Aper_flux_5",(cpl_size)i,&nl) -
                    tot[i]);
        diff /= cpl_table_get_float(tab,"Aper_flux_5_err",(cpl_size)i,&nl);
        cpl_test_lt(diff,1.6);
    }

    /* Run classify and test the values of the classification */

    retval = imcore_classify(res->catalogue,NULL,5,6);
    cpl_test_eq(retval,CASU_OK);
    pl = casu_tfits_get_ehu(res->catalogue);
    cpl_test_rel(cpl_propertylist_get_float(pl,"ESO QC IMAGE_SIZE"),4.47,0.02);
    cpl_test_eq(cpl_propertylist_get_bool(pl,"ESO DRS CLASSIFD"),1);
    cpl_test_rel( cpl_propertylist_get_float(pl, "APCOR3"),  0.132, 0.01);
    for (i = 0; i < NTEST; i++) {
        val = cpl_table_get(casu_tfits_get_table(res->catalogue),"Classification",
                            (cpl_size)i,NULL);
        cpl_test_rel(val,-1.0,0.001);
    }

    /* Get out of here */

    casu_fits_delete(inf);
    casu_fits_delete(inconf);
    casu_tfits_delete(res->catalogue);
    cpl_free(res);
    return(cpl_test_end(0));

}
