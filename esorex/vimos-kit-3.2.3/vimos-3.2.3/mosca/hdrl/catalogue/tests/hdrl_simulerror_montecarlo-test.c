/*
 * This file is part of the HDRL
 * Copyright (C) 2016 European Southern Observatory
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
 *
 *
 * hdrl_simulerror_montecarlo-test.c
 *
 *  Created on: Jul 7, 2016
 *      Author: agabasch
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "hdrl_fringe.h"
#include "hdrl_prototyping.h"
#include "hdrl_random.h"
#include <cpl.h>
#include <math.h>

#include "../imcore.h"
#include "../casu_utilfunctions.h"
#include "../../hdrl_catalogue.h"


#define COMP_TOL_REL 1.0/3.0
#define COMP_TOL_ABS 1e-2
#define IMG_XSIZE 120
#define IMG_YSIZE 180




/*----------------------------------------------------------------------------*/
/**
  @defgroup hdrl_simulerror_montecarlo  Testing of the Catalogue error module
 */
/*----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/
void test_basic(void);






void test_basic(void){
    int status;
    cpl_image *im, *bkg, *cnf;
    cpl_table *tab, *tabfinal = NULL;

    double xpos = 80.0;
    double ypos = 100.0;
    double norm = 3000.0;
    double norm2,  sigma=2.0, sky=500.0;
    casu_fits *inf, *inconf;
    cpl_propertylist *pl;
    hdrl_imcore_result * res = NULL;
    res = cpl_calloc(sizeof(*res), 1);

    /* Generate a field with some stars and a confidence map */
    bkg    = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    im     = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    cnf    = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    /* Generate confidence uniform confidence map */
    cpl_image_add_scalar(cnf, 100.);
    /* Generate sky and a gaussian source without errors */
    norm2 = 2.0*CPL_MATH_PI*sigma*sigma;

    cpl_image_add_scalar(bkg, sky);
    cpl_image_fill_gaussian(im, xpos, ypos, norm*norm2, sigma,sigma);
    cpl_image_add(bkg, im);

    pl = cpl_propertylist_new();
    inf    = casu_fits_wrap(im , NULL, NULL, pl);
    inconf = casu_fits_wrap(cnf , NULL, NULL, NULL);

    /* Give it a WCS */
    cpl_propertylist_update_string(pl,"CTYPE1", "RA---TAN");
    cpl_propertylist_update_string(pl,"CTYPE2", "DEC--TAN");
    cpl_propertylist_update_double(pl,"CRVAL1", 30.0);
    cpl_propertylist_update_double(pl,"CRVAL2", 12.0);
    cpl_propertylist_update_double(pl,"CRPIX1", 512.0);
    cpl_propertylist_update_double(pl,"CRPIX2", 512.0);
    cpl_propertylist_update_double(pl,"CD1_1", -1.0/3600);
    cpl_propertylist_update_double(pl,"CD1_2", 0.0);
    cpl_propertylist_update_double(pl,"CD2_1", 0.0);
    cpl_propertylist_update_double(pl,"CD2_2", 1.0/3600);
    cpl_propertylist_update_int(pl, "NAXIS1", IMG_XSIZE);
    cpl_propertylist_update_int(pl, "NAXIS2", IMG_YSIZE);
    cpl_wcs *wcs = cpl_wcs_new_from_propertylist(pl);
    cpl_propertylist_delete(pl);

    /*Perform the montecarlo simulation*/
    for(int iterate = 0; iterate < 100; iterate++) {
        /* Add realistic poisonnian noise */
        hdrl_random_state * rng = hdrl_random_state_new(1, NULL);
        int size = cpl_image_get_size_x(bkg) * cpl_image_get_size_y(bkg);
        float * pim = cpl_image_get_data_float(im);
        float * pbkg= cpl_image_get_data_float(bkg);

        for(int i=0;i<size;i++) {
            pim[i] = (float)hdrl_random_poisson(rng, (double)pbkg[i]);
        }
        hdrl_random_state_delete(rng);

        //cpl_image_save(bkg, "image_noerror.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
        //cpl_image_save(im, "image_error.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);

        /* Run imcore */
        status = CASU_OK;
        hdrl_casu_imcore(inf, inconf, wcs, 5, 2.5, 0, 3.0, 1, 32,
                         HDRL_CATALOGUE_ALL, 3.0, 1.0, HDRL_SATURATION_INIT,
                         res, &status);

        tab = casu_tfits_get_table(res->catalogue);

        if(iterate == 0) {
            tabfinal = cpl_table_duplicate(tab);
        }
        else {
            cpl_table_insert(tabfinal, tab, iterate);
        }

        cpl_image_delete(res->segmentation_map);
        cpl_image_delete(res->background);
        casu_tfits_delete(res->catalogue);
    }

    /* Do the checks */

    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "X_coordinate_err"), cpl_table_get_column_stdev(tabfinal, "X_coordinate"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Y_coordinate_err"), cpl_table_get_column_stdev(tabfinal, "Y_coordinate"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Peak_height_err"), cpl_table_get_column_stdev(tabfinal, "Peak_height"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_1_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_1"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_2_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_2"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_3_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_3"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_4_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_4"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_5_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_5"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_6_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_6"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_7_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_7"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_8_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_8"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_9_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_9"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_10_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_10"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_11_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_11"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_12_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_12"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Aper_flux_13_err"), cpl_table_get_column_stdev(tabfinal, "Aper_flux_13"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Petr_flux_err"), cpl_table_get_column_stdev(tabfinal, "Petr_flux"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Kron_flux_err"), cpl_table_get_column_stdev(tabfinal, "Kron_flux"), COMP_TOL_REL);
    cpl_test_rel(cpl_table_get_column_mean(tabfinal, "Half_flux_err"), cpl_table_get_column_stdev(tabfinal, "Half_flux"), COMP_TOL_REL);

    //cpl_table_save(tabfinal,    NULL, NULL, "tabfinal.fits",    CPL_IO_CREATE);

    //Cleanup

    casu_fits_delete(inf);
    casu_fits_delete(inconf);
    cpl_wcs_delete(wcs);
    cpl_table_delete(tabfinal);

    cpl_image_delete(bkg);
    cpl_free(res);




}

/*----------------------------------------------------------------------------*/
/**
 @brief   Unit tests of hdrl_catalogue module
 **/
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    test_basic();
    return cpl_test_end(0);

}

