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
#include <math.h>

#include <cpl_init.h>
#include <cpl_test.h>

#include "../imcore.h"
#include "../casu_utilfunctions.h"
#include "../../hdrl_catalogue.h"

#define COMP_TOL_REL 1e-4
#define COMP_TOL_ABS 1e-2
#define IMG_XSIZE 110*2
#define IMG_YSIZE 150*2
#define CELL_SIZE 32
#define APER_FLUX_NUM "Aper_flux_3" //string
#define NUM_FACTORS 5 
#define FACTORS {2.0, 10.0, 100.0, 1.0e4, 1.0e6}
#define NUM_BIASES 5
#define BIASES {2.0, 100.0, 5000.0, 5.0e4, 5.0e5}

/*----------------------------------------------------------------------------*/
/**
  @brief Check effect of multiplication on a given "Aper_flux",
  "Isophotal_flux", "FWHM", "Kron_radius", "Classification", "X_coordinate", 
  "X_coordinate_err", "Y_coordinate" and "Y_coordinate_err".
  @note Magic numerical constants inherited from 'casu_imcore-test.c'
  @return cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrl_casumul_compute(double factor) {
    intptr_t retval;
    int status, nl;
    cpl_image *im, *bkg, *cnf, *bkg2;
    cpl_table *tab, *tab_p2;
    
    double xpos = 100.0;
    double ypos = 100.0;
    double norm = 5000.0;
    double norm2,  sigma=2.0, sky=500.0;
    casu_fits *inf, *inf2, *inconf;
    cpl_propertylist *pl;

    /* Check inherited status */
    status = CASU_FATAL;
    hdrl_imcore_result * res    = cpl_calloc(sizeof(*res), 1);
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
    bkg    = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    im     = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    cnf    = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    norm2 = 2.0*CPL_MATH_PI*sigma*sigma;
    cpl_image_fill_noise_uniform(bkg, -10.0, 10.0);
    cpl_image_add_scalar(bkg, sky);
    cpl_image_fill_noise_uniform(cnf, 99.9, 100.1);

    cpl_image_fill_gaussian(im, xpos, ypos, norm*norm2, sigma,sigma);
    cpl_image_add(bkg, im);

    cpl_image_delete(im);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   //debug
        cpl_image_save(bkg, "bkg.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);

    bkg2 = cpl_image_duplicate(bkg);
    
    //MUL                       
    cpl_image_multiply_scalar(bkg2, factor);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   //debug
        cpl_image_save(bkg2, "bkg2.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    
    pl = cpl_propertylist_new();
    inf    = casu_fits_wrap(bkg , NULL, NULL, pl);
    inf2   = casu_fits_wrap(bkg2, NULL, NULL, pl);
    inconf = casu_fits_wrap(cnf , NULL, NULL, NULL);
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
        cpl_image_save(res->segmentation_map, "rsegmap.fits", CPL_TYPE_FLOAT,
                        NULL, CPL_IO_CREATE);
        cpl_image_save(res->background, "rbkg.fits", CPL_TYPE_FLOAT, NULL,
                       CPL_IO_CREATE);
    }
    casu_fits_delete(inf);
    
    /* Run imcore for Mul*/
    status = CASU_OK;
    retval = hdrl_casu_imcore(inf2, inconf, wcs, 5, 1.5, 0, 5.0, 1, CELL_SIZE, 
                              HDRL_CATALOGUE_ALL, 3.0, 
                              1.0, HDRL_SATURATION_INIT, res_p2, &status);
    cpl_test_eq(status, CASU_OK);
    cpl_test_eq(status, retval);
    cpl_test_nonnull(res_p2->catalogue);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   {   //debug
        cpl_image_save(res_p2->segmentation_map, "rsegmap_p2.fits",
                        CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
        cpl_image_save(res_p2->background, "rbkg_p2.fits", CPL_TYPE_FLOAT, NULL,
                       CPL_IO_CREATE);
    }
    casu_fits_delete(inf2);
    casu_fits_delete(inconf);
    cpl_wcs_delete(wcs);

    //Tests
    cpl_test_image_rel(res->segmentation_map, res_p2->segmentation_map,
                       COMP_TOL_REL);
    
    tab    = casu_tfits_get_table(res   ->catalogue);
    tab_p2 = casu_tfits_get_table(res_p2->catalogue);
    cpl_test_nonnull(tab);
    cpl_test_nonnull(tab_p2);
    
    cpl_test_rel( cpl_table_get_float(tab,    APER_FLUX_NUM, (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, APER_FLUX_NUM, (cpl_size)0, &nl) / factor,
                  COMP_TOL_REL );
    cpl_test_rel( cpl_table_get_float(tab,    "Isophotal_flux", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Isophotal_flux", (cpl_size)0, &nl) / factor,
                  COMP_TOL_REL );
    cpl_test_rel( cpl_table_get_float(tab,    "FWHM", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "FWHM", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_rel( cpl_table_get_float(tab,    "Kron_radius", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Kron_radius", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_eq ( cpl_table_get_float(tab,    "Classification", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Classification", (cpl_size)0, &nl) );
                  
    cpl_test_rel( cpl_table_get_float(tab,    "X_coordinate", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "X_coordinate", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_abs( cpl_table_get_float(tab,    "X_coordinate_err", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "X_coordinate_err", (cpl_size)0, &nl),
                  COMP_TOL_ABS );
    cpl_test_rel( cpl_table_get_float(tab,    "Y_coordinate", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Y_coordinate", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_abs( cpl_table_get_float(tab,    "Y_coordinate_err", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Y_coordinate_err", (cpl_size)0, &nl),
                  COMP_TOL_ABS );
                      
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   {   //debug
        cpl_table_save(tab,    NULL, NULL, "rcat.fits",    CPL_IO_CREATE);
        cpl_table_save(tab_p2, NULL, NULL, "rcat_p2.fits", CPL_IO_CREATE);
    }


    //Cleanup
    cpl_image_delete(res->segmentation_map);
    cpl_image_delete(res->background);
    casu_tfits_delete(res->catalogue);
    cpl_free(res);

    cpl_image_delete(res_p2->segmentation_map);
    cpl_image_delete(res_p2->background);
    casu_tfits_delete(res_p2->catalogue);
    cpl_free(res_p2);

    return cpl_error_get_code();
    
}

/*----------------------------------------------------------------------------*/
/**
  @brief Check effect of addition on a given "Aper_flux",
  "Isophotal_flux", "FWHM", "Kron_radius", "Classification", "X_coordinate", 
  "X_coordinate_err", "Y_coordinate" and "Y_coordinate_err".
  @note Magic numerical constants inherited from 'casu_imcore-test.c'
  @return cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrl_casuadd_compute(double bias) {
    intptr_t retval;
    int status, nl;
    cpl_image *im, *bkg, *cnf, *bkg2;
    cpl_table *tab, *tab_p2;
    
    double xpos = 100.0;
    double ypos = 100.0;
    double norm = 1000.0;
    double norm2,  sigma=2.0, sky=500.0;
    casu_fits *inf, *inf2, *inconf;
    cpl_propertylist *pl;

    /* Check inherited status */
    status = CASU_FATAL;
    hdrl_imcore_result * res    = cpl_calloc(sizeof(*res), 1);
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
    bkg    = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    im     = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    cnf    = cpl_image_new(IMG_XSIZE, IMG_YSIZE, CPL_TYPE_FLOAT);
    norm2 = 2.0*CPL_MATH_PI*sigma*sigma;
    cpl_image_fill_noise_uniform(bkg, -10.0, 10.0);
    cpl_image_add_scalar(bkg, sky);
    cpl_image_fill_noise_uniform(cnf, 99.9, 100.1);

    cpl_image_fill_gaussian(im, xpos, ypos, norm*norm2, sigma,sigma);
    cpl_image_add(bkg, im);

    cpl_image_delete(im);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   //debug
        cpl_image_save(bkg, "bkg.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);

    bkg2 = cpl_image_duplicate(bkg);
    
    //ADD                       
    cpl_image_add_scalar(bkg2, bias);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   //debug
        cpl_image_save(bkg2, "bkg2.fits", CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
    
    pl = cpl_propertylist_new();
    inf    = casu_fits_wrap(bkg , NULL, NULL, pl);
    inf2   = casu_fits_wrap(bkg2, NULL, NULL, pl);
    inconf = casu_fits_wrap(cnf , NULL, NULL, NULL);
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
        cpl_image_save(res->segmentation_map, "rsegmap.fits", CPL_TYPE_FLOAT,
                        NULL, CPL_IO_CREATE);
        cpl_image_save(res->background, "rbkg.fits", CPL_TYPE_FLOAT, NULL,
                       CPL_IO_CREATE);
    }
    casu_fits_delete(inf);
    
    /* Run imcore for ADD*/
    status = CASU_OK;
    retval = hdrl_casu_imcore(inf2, inconf, wcs, 5, 1.5, 0, 5.0, 1, CELL_SIZE, 
                              HDRL_CATALOGUE_ALL, 3.0, 
                              1.0, HDRL_SATURATION_INIT, res_p2, &status);
    cpl_test_eq(status, CASU_OK);
    cpl_test_eq(status, retval);
    cpl_test_nonnull(res_p2->catalogue);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   {   //debug
        cpl_image_save(res_p2->segmentation_map, "rsegmap_p2.fits",
                        CPL_TYPE_FLOAT, NULL, CPL_IO_CREATE);
        cpl_image_save(res_p2->background, "rbkg_p2.fits", CPL_TYPE_FLOAT,
                       NULL, CPL_IO_CREATE);
    }
    casu_fits_delete(inf2);
    casu_fits_delete(inconf);
    cpl_wcs_delete(wcs);

    //Tests
    cpl_test_image_rel(res->segmentation_map, res_p2->segmentation_map,
                       COMP_TOL_REL);
    
    tab    = casu_tfits_get_table(res   ->catalogue);
    tab_p2 = casu_tfits_get_table(res_p2->catalogue);
    cpl_test_nonnull(tab);
    cpl_test_nonnull(tab_p2);

    cpl_test_rel( cpl_table_get_float(tab,    APER_FLUX_NUM, (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, APER_FLUX_NUM, (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_rel( cpl_table_get_float(tab,    "Isophotal_flux", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Isophotal_flux", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_rel( cpl_table_get_float(tab,    "FWHM", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "FWHM", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_rel( cpl_table_get_float(tab,    "Kron_radius", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Kron_radius", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_eq ( cpl_table_get_float(tab,    "Classification", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Classification", (cpl_size)0, &nl) );
                  
    cpl_test_rel( cpl_table_get_float(tab,    "X_coordinate", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "X_coordinate", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_abs( cpl_table_get_float(tab,    "X_coordinate_err", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "X_coordinate_err", (cpl_size)0, &nl),
                  COMP_TOL_ABS );
    cpl_test_rel( cpl_table_get_float(tab,    "Y_coordinate", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Y_coordinate", (cpl_size)0, &nl),
                  COMP_TOL_REL );
    cpl_test_abs( cpl_table_get_float(tab,    "Y_coordinate_err", (cpl_size)0, &nl),
                  cpl_table_get_float(tab_p2, "Y_coordinate_err", (cpl_size)0, &nl),
                  COMP_TOL_ABS );

    if (cpl_msg_get_level() == CPL_MSG_DEBUG)   {   //debug
        cpl_table_save(tab,    NULL, NULL, "rcat.fits",    CPL_IO_CREATE);
        cpl_table_save(tab_p2, NULL, NULL, "rcat_p2.fits", CPL_IO_CREATE);
    }


    //Cleanup
    cpl_image_delete(res->segmentation_map);
    cpl_image_delete(res->background);
    casu_tfits_delete(res->catalogue);
    cpl_free(res);
    
    cpl_image_delete(res_p2->segmentation_map);
    cpl_image_delete(res_p2->background);
    casu_tfits_delete(res_p2->catalogue);
    cpl_free(res_p2);

    return cpl_error_get_code();
    
}


/*----------------------------------------------------------------------------*/
/**
  @brief   Unit tests of addition or multiplication of some scalar value
 **/
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    double bias[] = BIASES;
    for(int i = 0; i < NUM_BIASES; i++)
        hdrl_casuadd_compute(bias[i]);

    double factor[] = FACTORS;
    for(int i = 0; i < NUM_FACTORS; i++)  
        hdrl_casumul_compute(factor[i]);
    
    cpl_test_error(CPL_ERROR_NONE);
    return cpl_test_end(0);
}
