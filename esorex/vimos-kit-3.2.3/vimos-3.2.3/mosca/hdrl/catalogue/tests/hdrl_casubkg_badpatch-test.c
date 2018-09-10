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
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpl_init.h>
#include <cpl_test.h>

#include "../imcore.h"
#include "../casu_utilfunctions.h"
#include "../../hdrl_catalogue.h"

#include <math.h>

#define NTEST 10
#define COMP_TOLERANCE 1e-2
#define CORNER_OFFSET 10
#define CORNER_REL_TOL 1e-2
#define X_OS_P1 60
#define X_OS_P2 95
#define STAR 7
#define PATCH_SIZE 100 //even number (patch is a square)
#define IMG_XSIZE 1100
#define IMG_YSIZE 1500
#define CELL_SIZE 32
#define MIN_RAMP 10.0
#define MAX_RAMP 100.0
#define APER_FLUX_NUM "Aper_flux_3" //string

/*----------------------------------------------------------------------------*/
/**
  @brief Check effect of bad pixel square patches in various conditions
  @note Magic numerical constants inherited from 'casu_imcore-test.c'
  @return cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrl_casubkg_badpatch_compute(void) {
    intptr_t retval, i, j, x_os = 0, y_os = 0; //nrow, ncol;
    int status, nl, rej;

    cpl_image *im, *bkg, *cnf, *cnf_p1, *cnf_p2, *aux;
    cpl_table *tab, *tab_p1, *tab_p2;
    
    double xpos[] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0, 1000.0};
    double ypos[] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0, 1000.0};
    double norm[] = {1000.0,100.0,200.0,500.0,550.0,600.0,650.0,700.0, 750.0, 800.0};
    double norm2,  sigma=2.0, sky=500.0; //diff, tot[NTEST];
    casu_fits *inf, *inconf, *inconf_p1, *inconf_p2;
    cpl_propertylist *pl;

    /* Check inherited status */
    status = CASU_FATAL;
    hdrl_imcore_result * res    = cpl_calloc(sizeof(*res), 1);
    hdrl_imcore_result * res_p1 = cpl_calloc(sizeof(*res), 1);
    hdrl_imcore_result * res_p2 = cpl_calloc(sizeof(*res), 1);
    retval = hdrl_casu_imcore(NULL, NULL, NULL, 5, 1.5, 1, 3.5, 1, CELL_SIZE, 
                              6, 1.0, 2.0,
							  HDRL_SATURATION_INIT, res, &status);
                              
    cpl_test_eq(status, CASU_FATAL);
    cpl_test_eq(status, retval);
    cpl_test_null(res->catalogue);

    cpl_image_delete(res->segmentation_map);
    cpl_image_delete(res->background);
    casu_tfits_delete(res->catalogue);
    cpl_free(res);

    /* Generate a field with some stars and a confidence map */
    bkg = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    im  = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    cnf = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    aux = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    norm2 = 2.0*CPL_MATH_PI*sigma*sigma;
    cpl_image_fill_noise_uniform(bkg, -10.0, 10.0);
    cpl_image_add_scalar(bkg, sky);
    cpl_image_fill_noise_uniform(cnf, 99.9, 100.1);
    cnf_p1 = cpl_image_duplicate(cnf);
    cnf_p2 = cpl_image_duplicate(cnf);
    for (i = 0; i < NTEST; i++) {
        cpl_image_fill_gaussian(im, xpos[i], ypos[i], norm[i]*norm2, sigma,sigma);
        //~ tot[i] = cpl_image_get_flux(im);
        cpl_image_add(bkg, im);
    }
    cpl_image_delete(im);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   //debug
        cpl_image_save(bkg, "bkg.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    
    //Add gradient
    double dlt_ramp = MAX_RAMP - MIN_RAMP;
    for(i = 1; i <= IMG_XSIZE; i++)
        for(j = 1; j <= IMG_YSIZE; j++)
            cpl_image_set(aux, i, j, MIN_RAMP+ dlt_ramp*(double)i/IMG_XSIZE);
    cpl_image_add(bkg, aux);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   //debug
        cpl_image_save(bkg, "bkg_grad.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    
    cpl_image_delete(aux);
    
    //NEW - add bad pixel square "patch"
    //Generation of confidence map "Patch Close"
    x_os = X_OS_P1;
    for(i = 0; i <= PATCH_SIZE; i++)
        for(j = 0; j <= PATCH_SIZE; j++)
            cpl_image_set(cnf_p1, xpos[STAR]-(PATCH_SIZE/2)+i +x_os, 
                                  ypos[STAR]-(PATCH_SIZE/2)+j +y_os, 0);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   //debug
        cpl_image_save(cnf_p1, "cnf_p1.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    
    //Generation of confidence map "Patch Away"                           
    x_os = X_OS_P2;
    for(i = 0; i <= PATCH_SIZE; i++)
        for(j = 0; j <= PATCH_SIZE; j++)
            cpl_image_set(cnf_p2, xpos[STAR]-(PATCH_SIZE/2)+i +x_os, 
                                  ypos[STAR]-(PATCH_SIZE/2)+j +y_os, 0);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   //debug
        cpl_image_save(cnf_p2, "cnf_p2.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    
    pl = cpl_propertylist_new();
    inf       = casu_fits_wrap(bkg   , NULL, NULL, pl);
    inconf    = casu_fits_wrap(cnf   , NULL, NULL, NULL);
    inconf_p1 = casu_fits_wrap(cnf_p1, NULL, NULL, NULL);
    inconf_p2 = casu_fits_wrap(cnf_p2, NULL, NULL, NULL);
    cpl_propertylist_delete(pl);
    
    /* Give it a WCS */
    pl = casu_fits_get_ehu(inf);
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

    /* Run imcore */
    status = CASU_OK;
    res = cpl_calloc(sizeof(*res), 1);
    retval = hdrl_casu_imcore(inf, inconf, wcs, 5, 1.5, 0, 5.0, 1, CELL_SIZE, 
                              HDRL_CATALOGUE_ALL, 3.0, 
                              1.0, HDRL_SATURATION_INIT, res, &status);
    cpl_test_eq(status, CASU_OK);
    cpl_test_eq(status, retval);
    cpl_test_nonnull(res->catalogue);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   {   //debug
        cpl_image_save(res->segmentation_map, "rsegmap.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
        cpl_image_save(res->background, "rbkg.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    }
    casu_fits_delete(inconf);
    
    /* Run imcore for PC*/
    status = CASU_OK;
    retval = hdrl_casu_imcore(inf, inconf_p1, wcs, 5, 1.5, 0, 5.0, 1, CELL_SIZE, 
                              HDRL_CATALOGUE_ALL, 3.0, 
                              1.0, HDRL_SATURATION_INIT, res_p1, &status);
    cpl_test_eq(status, CASU_OK);
    cpl_test_eq(status, retval);
    cpl_test_nonnull(res_p1->catalogue);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   {   //debug
        cpl_image_save(res_p1->segmentation_map, "rsegmap_p1.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
        cpl_image_save(res_p1->background, "rbkg_p1.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    }
    casu_fits_delete(inconf_p1);

    /* Run imcore for PA*/
    status = CASU_OK;
    retval = hdrl_casu_imcore(inf, inconf_p2, wcs, 5, 1.5, 0, 5.0, 1, CELL_SIZE, 
                              HDRL_CATALOGUE_ALL, 3.0, 
                              1.0, HDRL_SATURATION_INIT, res_p2, &status);
    cpl_test_eq(status, CASU_OK);
    cpl_test_eq(status, retval);
    cpl_test_nonnull(res_p2->catalogue);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   {   //debug
        cpl_image_save(res_p2->segmentation_map, "rsegmap_p2.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
        cpl_image_save(res_p2->background, "rbkg_p2.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    }
    casu_fits_delete(inconf_p2);
    casu_fits_delete(inf);
    cpl_wcs_delete(wcs);

    //Tests
    cpl_test_image_rel(res->background, res_p1->background, COMP_TOLERANCE);
    cpl_test_image_rel(res->background, res_p2->background, COMP_TOLERANCE);
    
    cpl_size corner_x = cpl_image_get_size_x(res->background) -CORNER_OFFSET;
    cpl_size corner_y = cpl_image_get_size_y(res->background) -CORNER_OFFSET;
    cpl_test_abs(cpl_image_get(res   ->background, corner_x, corner_y, &rej),
                 cpl_image_get(res_p1->background, corner_x, corner_y, &rej),
                 CORNER_REL_TOL);
    cpl_test_abs(cpl_image_get(res   ->background, corner_x, corner_y, &rej),
                 cpl_image_get(res_p2->background, corner_x, corner_y, &rej),
                 CORNER_REL_TOL);

    cpl_test_image_rel(res->segmentation_map, res_p1->segmentation_map, COMP_TOLERANCE);
    cpl_test_image_rel(res->segmentation_map, res_p2->segmentation_map, COMP_TOLERANCE);
    
    
    /* Check the results. Start by checking the number of rows and columns. 
    Sort the table by X */
    tab    = casu_tfits_get_table(res   ->catalogue);
    tab_p1 = casu_tfits_get_table(res_p1->catalogue);
    tab_p2 = casu_tfits_get_table(res_p2->catalogue);
    cpl_test_nonnull(tab);
    cpl_test_nonnull(tab_p1);
    cpl_test_nonnull(tab_p2);
    
    for (i = 0; i < NTEST; i++) {
        cpl_test_rel( cpl_table_get_float(tab,    APER_FLUX_NUM, (cpl_size)i, &nl),
                      cpl_table_get_float(tab_p1, APER_FLUX_NUM, (cpl_size)i, &nl),
                      COMP_TOLERANCE );
        cpl_test_rel( cpl_table_get_float(tab,    APER_FLUX_NUM, (cpl_size)i, &nl),
                      cpl_table_get_float(tab_p2, APER_FLUX_NUM, (cpl_size)i, &nl),
                      COMP_TOLERANCE );
    }
    
    //~ ncol = cpl_table_get_ncol(tab);
    //~ cpl_test_eq(ncol, 80);
    //~ nrow = cpl_table_get_nrow(tab);
    //~ cpl_test_eq(nrow, NTEST);
    //~ pl = cpl_propertylist_new();
    //~ cpl_propertylist_append_bool(pl, "X_coordinate", 0);
    //~ cpl_table_sort(tab, pl);
    //~ cpl_propertylist_delete(pl);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   {   //debug
        cpl_table_save(tab,    NULL, NULL, "rcat.fits",    CPL_IO_CREATE);
        cpl_table_save(tab_p1, NULL, NULL, "rcat_p1.fits", CPL_IO_CREATE);
        cpl_table_save(tab_p2, NULL, NULL, "rcat_p2.fits", CPL_IO_CREATE);
    }

    /* Test the column content of the table */
    //~ for (i = 0; i < NTEST; i++) {
        //~ cpl_test_abs(xpos[i],  cpl_table_get_float(tab, "X_coordinate", (cpl_size)i, &nl),  0.2);
        //~ cpl_test_abs(ypos[i],  cpl_table_get_float(tab, "Y_coordinate", (cpl_size)i, &nl),  0.2);
        //~ diff = fabs( cpl_table_get_float(tab, "Aper_flux_5", (cpl_size)i, &nl) - tot[i]);
        //~ diff /= cpl_table_get_float(tab, "Aper_flux_5_err", (cpl_size)i, &nl);
        //~ cpl_test_lt(diff, 1.5);
        //~ cpl_test_eq( cpl_table_get_float(tab, "Classification", (cpl_size)i, &nl),  -1.0);
    //~ }
    

    //Cleanup
    cpl_image_delete(res->segmentation_map);
    cpl_image_delete(res->background);
    casu_tfits_delete(res->catalogue);
    cpl_free(res);
    
    cpl_image_delete(res_p1->segmentation_map);
    cpl_image_delete(res_p1->background);
    casu_tfits_delete(res_p1->catalogue);
    cpl_free(res_p1);
    
    cpl_image_delete(res_p2->segmentation_map);
    cpl_image_delete(res_p2->background);
    casu_tfits_delete(res_p2->catalogue);
    cpl_free(res_p2);

    return cpl_error_get_code();
    
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Unit tests of bad patches in background
 **/
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);
    hdrl_casubkg_badpatch_compute();
    cpl_test_error(CPL_ERROR_NONE);
    return cpl_test_end(0);
}

