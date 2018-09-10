/*
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                    Includes
 -----------------------------------------------------------------------------*/

#include "hdrl_bpm_2d.h"
#include "hdrl_image.h"

#include <cpl.h>
#include <math.h>

/*----------------------------------------------------------------------------*/
/**
  @defgroup hdrl_bpm_2d_test    Testing of the HDRL bpm_2d
 */
/*----------------------------------------------------------------------------*/

void test_parlist(void)
{
    /* parameter parsing smoketest */
    hdrl_parameter * hpar;
    hdrl_parameter * fil_def =
        hdrl_bpm_2d_parameter_create_filtersmooth(4, 5, 6,
                                                  CPL_FILTER_MEDIAN,
                                                  CPL_BORDER_NOP,
                                                  7, 9);
    hdrl_parameter * leg_def =
        hdrl_bpm_2d_parameter_create_legendresmooth(4, 5, 6,
                                                    20, 21,
                                                    11, 12, 2, 10);
    cpl_parameterlist * pos = hdrl_bpm_2d_parameter_create_parlist(
                "RECIPE", "bpm", "FILTER", fil_def, leg_def);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(cpl_parameterlist_get_size(pos), 17);

    hpar = hdrl_bpm_2d_parameter_parse_parlist(pos, "RECIPE.invalid");
    cpl_test_null(hpar);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);

    hpar = hdrl_bpm_2d_parameter_parse_parlist(pos, "RECIPE.bpm");
    cpl_parameterlist_delete(pos);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(hdrl_bpm_2d_parameter_get_filter(hpar), CPL_FILTER_MEDIAN);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_border(hpar), CPL_BORDER_NOP);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_kappa_low(hpar), 4);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_kappa_high(hpar), 5);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_maxiter(hpar), 6);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_smooth_x(hpar), 7);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_smooth_y(hpar), 9);
    /* TODO wrong par error? */
    cpl_test_eq(hdrl_bpm_2d_parameter_get_steps_x(hpar), 0);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_steps_y(hpar), 0);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_filter_size_x(hpar), 0);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_filter_size_y(hpar), 0);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_order_x(hpar), 0);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_order_y(hpar), 0);
    hdrl_parameter_delete(hpar);

    pos = hdrl_bpm_2d_parameter_create_parlist(
                "RECIPE", "bpm", "LEGENDRE", fil_def, leg_def);
    cpl_test_error(CPL_ERROR_NONE);

    hpar = hdrl_bpm_2d_parameter_parse_parlist(pos, "RECIPE.bpm");
    cpl_parameterlist_delete(pos);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(hdrl_bpm_2d_parameter_get_kappa_low(hpar), 4);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_kappa_high(hpar), 5);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_maxiter(hpar), 6);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_steps_x(hpar), 20);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_steps_y(hpar), 21);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_filter_size_x(hpar), 11);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_filter_size_y(hpar), 12);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_order_x(hpar), 2);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_order_y(hpar), 10);
    /* TODO wrong par error? */
    cpl_test_eq(hdrl_bpm_2d_parameter_get_filter(hpar), CPL_FILTER_MEDIAN);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_border(hpar), CPL_BORDER_FILTER);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_smooth_x(hpar), 0);
    cpl_test_eq(hdrl_bpm_2d_parameter_get_smooth_y(hpar), 0);
    hdrl_parameter_delete(hpar);

    hdrl_parameter_delete(leg_def);
    hdrl_parameter_delete(fil_def);
}

/*----------------------------------------------------------------------------*/
/**
  @brief Check hdrl_bpm_2d_compute() in various conditions
  @return cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrl_bpm_2d_test_compute(void)
{
    cpl_image       *   data = NULL;
    cpl_image       *   errors = NULL;
    cpl_mask        *   mask_out = NULL;
    hdrl_parameter  *   bpm_param = NULL;

    /* Create BPM parameters */
    bpm_param = hdrl_bpm_2d_parameter_create_filtersmooth(3., 3., 2, 
            CPL_FILTER_MEDIAN, CPL_BORDER_FILTER, 3, 3) ;

    { /* Test sigmaclipped mean */
        /* gauss mean 100 sigma 3.5 and 2 outliers */
        double values[] = {92, 93, 94, 94, 95, 95, 96, 96, 96, 97,
                        97, 97, 97, 98, 98, 98, 98, 99, 99, 99,
                        99, 100, 100, 100, 100, 100, 101, 101, 101, 101,
                        102, 102, 102, 102, 103, 103, 103, 103, 104, 104,
                        104, 105, 105, 106, 106, 107, 108, 500, 600 };

        data = cpl_image_wrap(7, 7, CPL_TYPE_DOUBLE, values);
        errors = cpl_image_new(7, 7, CPL_TYPE_DOUBLE);
        cpl_image_add_scalar(errors, 1);

        cpl_image_set(errors, 7, 7, 100000.);
        cpl_image_set(errors, 6, 7, 10000.);

        hdrl_image * sigimage = hdrl_image_create(data, errors);

        mask_out = hdrl_bpm_2d_compute(sigimage, bpm_param) ;
        cpl_mask_delete(mask_out);
        
        hdrl_parameter_delete(bpm_param) ;
        bpm_param = hdrl_bpm_2d_parameter_create_legendresmooth(3., 3.,
                2,  20, 20, 11, 11, 3, 3) ;

        mask_out = hdrl_bpm_2d_compute(sigimage, bpm_param) ;

        //cpl_mask_save(mask_out, "bpm1.fits", NULL, CPL_IO_CREATE);

        cpl_mask_delete(mask_out);
        cpl_image_unwrap(data);
        cpl_image_delete(errors);
        hdrl_image_delete(sigimage);
    }

    {
        cpl_mask * data_bpm = cpl_mask_new(200, 300);
        data = cpl_image_new(200, 300, CPL_TYPE_FLOAT);

        cpl_image_fill_noise_uniform(data, 90, 110);
        cpl_image_set(data, 50, 50, 300.);
        cpl_image_set(data, 100, 100, 300.);
        cpl_image_set(data, 150, 150, 300.);
        cpl_image_set(data, 110, 260, 300.);

        cpl_mask_set(data_bpm, 120,120, CPL_BINARY_1);
        cpl_mask_set(data_bpm, 120,121, CPL_BINARY_1);
        cpl_mask_set(data_bpm, 120,122, CPL_BINARY_1);
        cpl_mask_set(data_bpm, 121,120, CPL_BINARY_1);
        cpl_mask_set(data_bpm, 121,121, CPL_BINARY_1);
        cpl_mask_set(data_bpm, 121,122, CPL_BINARY_1);
        cpl_mask_set(data_bpm, 122,120, CPL_BINARY_1);
        cpl_mask_set(data_bpm, 122,121, CPL_BINARY_1);
        cpl_mask_set(data_bpm, 122,122, CPL_BINARY_1);
        /*set one outlier on a bad pixel*/
        cpl_image_set(data, 122, 122, 300.);

        cpl_image_reject_from_mask(data, data_bpm);

        errors=cpl_image_power_create(data, 0.5);
        hdrl_image * image = hdrl_image_create(data, errors);

        /* Note that for CPL_FILTER_STDEV one gets more false positives!! */
        hdrl_parameter_delete(bpm_param) ;
        bpm_param = hdrl_bpm_2d_parameter_create_filtersmooth(3., 3., 5, 
                CPL_FILTER_MEDIAN, CPL_BORDER_FILTER, 3, 3) ;
        mask_out = hdrl_bpm_2d_compute(image, bpm_param) ;
        cpl_test_eq(cpl_mask_get(mask_out, 50, 50),   CPL_BINARY_1);
        cpl_test_eq(cpl_mask_get(mask_out, 100, 100), CPL_BINARY_1);
        cpl_test_eq(cpl_mask_get(mask_out, 150, 150), CPL_BINARY_1);
        cpl_test_eq(cpl_mask_get(mask_out, 110, 260), CPL_BINARY_1);
        cpl_test_eq(cpl_mask_get(mask_out, 122, 122), CPL_BINARY_0);
        cpl_mask_delete(mask_out);

        hdrl_parameter_delete(bpm_param) ;
        bpm_param = hdrl_bpm_2d_parameter_create_legendresmooth(3., 3.,
                5, 20, 20, 11, 11, 3, 3) ;
        mask_out = hdrl_bpm_2d_compute(image, bpm_param) ;
        cpl_test_eq(cpl_mask_get(mask_out, 50, 50),   CPL_BINARY_1);
        cpl_test_eq(cpl_mask_get(mask_out, 100, 100), CPL_BINARY_1);
        cpl_test_eq(cpl_mask_get(mask_out, 150, 150), CPL_BINARY_1);
        cpl_test_eq(cpl_mask_get(mask_out, 110, 260), CPL_BINARY_1);
        cpl_test_eq(cpl_mask_get(mask_out, 122, 122), CPL_BINARY_0);

        cpl_mask_delete(data_bpm);
        cpl_mask_delete(mask_out);
        cpl_image_delete(data);
        cpl_image_delete(errors);
        hdrl_image_delete(image);
    }
    hdrl_parameter_delete(bpm_param) ;
    return cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Unit tests of BPM module
 **/
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);
    hdrl_bpm_2d_test_compute();
    test_parlist();
    cpl_test_error(CPL_ERROR_NONE);
    return cpl_test_end(0);
}
