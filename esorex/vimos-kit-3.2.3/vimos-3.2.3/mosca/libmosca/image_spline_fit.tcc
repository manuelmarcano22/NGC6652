/* $Id: image_spline_fit.tcc,v 1.7 2013-09-09 13:14:32 cgarcia Exp $
 *
 * This file is part of the MOSCA library
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-09-09 13:14:32 $
 * $Revision: 1.7 $
 * $Name: not supported by cvs2svn $
 */


#ifndef IMAGE_SPLINE_FIT_TCC
#define IMAGE_SPLINE_FIT_TCC

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include "mosca_image.h"
#include "image_spline_fit.h"

cpl_vector * fit_cubic_bspline(cpl_vector * values, int nknots, double threshold);


/**
 * @brief
 *   TODO 
 *
 * @param TODO
 * @return TODO
 *
 */
template<typename T>
void mosca::image_cubicspline_1d_fit(mosca::image& image, 
                                     int nknots, double flux_threshold,
                                     mosca::axis fitting_axis)
{
    cpl_size ny = image.size_y();
    cpl_size nx = image.size_x();
    mosca::axis image_axis;
    int collapse_dir;

    /* If the fitting direction is given in terms of dispersion or spatial
       axis, then decide whether it is X or Y. If not, it is already X or Y */
    image_axis = image.axis_to_image(fitting_axis);
    
    /* Lines are contiguous pixels along the fitting direction.
       i.e., if image_axis == mosca::X_AXIS the lines are rows and there are
       ny lines.
       Each line contains nfit pixels. For image_axis == mosca::X_AXIS this
       means nfit = nx. */
    if(image_axis == mosca::X_AXIS)
        collapse_dir = 0;
    else
        collapse_dir = 1;
        
    
    //We collapse in the opposite direction of the fitting
    cpl_image * collapse_median = 
          cpl_image_collapse_median_create(image.get_cpl_image(),
                                           collapse_dir, 0, 0);
    
    //We get the vector of the collapsed image
    cpl_vector * collapse_line;
    if(image_axis == mosca::X_AXIS)
        collapse_line = cpl_vector_new_from_image_row(collapse_median, 1);
    else
        collapse_line = cpl_vector_new_from_image_column(collapse_median, 1);

    //We fit a spline to the line-image
    cpl_vector * spline_fit = fit_cubic_bspline(collapse_line, nknots, flux_threshold);

    //Place in each column or row the spline fit
    T *      p_image      = image.get_data<T>();
    double * p_spline_fit = cpl_vector_get_data(spline_fit);
    for (cpl_size j = 0; j< ny; ++j)
    {
        for (cpl_size i = 0; i< nx; ++i, ++p_image)
        {
            if(image_axis == mosca::X_AXIS)
                *p_image = p_spline_fit[i];
            else
                *p_image =  p_spline_fit[j];
        }
    }

    cpl_image_delete(collapse_median);    
    cpl_vector_delete(collapse_line);
    cpl_vector_delete(spline_fit);
}

//TODO: Make sure that this is not exported
cpl_vector * fit_cubic_bspline(cpl_vector * values, int nknots, double threshold)
{
    //Allocate result
    cpl_vector * result;
    cpl_size     npix = cpl_vector_get_size(values);
    result = cpl_vector_new(npix);

    //This is valid for cubic splines, if not nbreak = ncoeffs + 2 - k, 
    //where k is the degree of the spline
    int ncoeffs = nknots + 2;
    
    /* Get the threshold in terms of the maximum */
    double max_value = cpl_vector_get_max(values);
    double thres_value = threshold * max_value;
    
    /* Create a "mask" of pixels not to use */
    cpl_array * values_mask = cpl_array_new(npix, CPL_TYPE_INT);
    int nfit = 0;
    for (cpl_size i = 0; i < npix; ++i)
        if(cpl_vector_get(values, i) < thres_value)
            cpl_array_set_int(values_mask, i, 0);
        else
        {
            cpl_array_set_int(values_mask, i, 1);
            nfit++;
        }

    /* allocate a cubic bspline workspace (k = 4) */
    gsl_bspline_workspace *bspl_wspc;
    gsl_vector *B;
    gsl_matrix * X;
    bspl_wspc = gsl_bspline_alloc(4, nknots);
    B = gsl_vector_alloc(ncoeffs);
    X = gsl_matrix_alloc(nfit, ncoeffs);
    
    /* allocate objects for the fitting */
    gsl_matrix *cov;
    gsl_vector * y_fit;
    gsl_vector * weigth;
    gsl_multifit_linear_workspace * mfit_wspc;
    gsl_vector *spl_coeffs;
    double chisq;
    y_fit = gsl_vector_alloc(nfit);
    weigth = gsl_vector_alloc(nfit);
    mfit_wspc = gsl_multifit_linear_alloc(nfit, ncoeffs);
    spl_coeffs = gsl_vector_alloc(ncoeffs);
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    
    /* use uniform breakpoints on [0, npix], which is the range of fitting */
    gsl_bspline_knots_uniform(0.0, (double)npix, bspl_wspc);
    

    /* construct the fit matrix X */
    for (int i = 0, ifit = 0; i < npix; ++i)
    {
        int null;
        if(cpl_array_get(values_mask, (cpl_size)i, &null) == 1)
        {
            double xi = i;
            double yi = cpl_vector_get(values, i);

            /* Fill the vector to fit */
            gsl_vector_set(y_fit, ifit, yi);
            gsl_vector_set(weigth, ifit, 1.);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, B, bspl_wspc);

            /* fill in row i of X */
            for (int j = 0; j < ncoeffs; ++j)
            {
                double Bj = gsl_vector_get(B, j);
                gsl_matrix_set(X, ifit, j, Bj);
            }
            ifit++;
        }
    } 
    
    /* do the fit */
    gsl_multifit_wlinear(X, weigth, y_fit, spl_coeffs, cov, &chisq, mfit_wspc);
    
    /* output the fit */

    for(int i = 0; i < npix; i ++)
    {
        double yi, yerr;
        gsl_bspline_eval((double)i, B, bspl_wspc);
        gsl_multifit_linear_est(B, spl_coeffs, cov, &yi, &yerr);
        cpl_vector_set(result, i, yi);
    }
    
    return result;
}

/**
 * Smooth an image fitting a polynomial in spatial direction 
 * @param image       Image to be smoothed
 * @param polyorder   the order of the polynomial to fit
 * 
 * This function will smooth a MOS image in the dispersion direction 
 * (here assumed to be the X axis) fitting each row with a polynomial.
 * TODO: merge this with image_smooth_fit_1d_pol_spa() and place it in MOSCA.
 * TODO: It is named image_smooth_fit_1d_pol but it is not yet generalized
 * to any direction, only to spatial.   
 */
template<typename T>
void mosca::image_pol_1d_fit(mosca::image& image, int polyorder, 
                             double threshold)
{

    //We collapse in the oposite direction of the fitting
    cpl_image * collapse_median = 
          cpl_image_collapse_median_create(image.get_cpl_image(),
                                           1, 0, 0);

    //We get the vector of the collapsed image
    cpl_vector * collapse_line;
    collapse_line = cpl_vector_new_from_image_column(collapse_median, 1);

    
    //We fit a polynomial to the line-image
    double max_value = cpl_vector_get_max(collapse_line);
    double thres_value = threshold * max_value;
    
    cpl_vector * filtered_values = cpl_vector_duplicate(collapse_line);
    cpl_vector * xpos = cpl_vector_duplicate(collapse_line);
    int i_size_fit = 0;
    size_t npix = cpl_vector_get_size(collapse_line);
    for (cpl_size i = 0; i < npix; ++i)
    {
        if(cpl_vector_get(collapse_line, i) > thres_value )
        {
            cpl_vector_set(filtered_values, i_size_fit, cpl_vector_get(collapse_line, i));
            cpl_vector_set(xpos, i_size_fit, (double)i);
            i_size_fit++;
        }
    }
    cpl_vector * values_final = 
            cpl_vector_extract(filtered_values, 0, i_size_fit - 1, 1);
    cpl_vector * xpos_final = 
            cpl_vector_extract(xpos, 0, i_size_fit - 1, 1);

    cpl_polynomial * pol_fit = cpl_polynomial_fit_1d_create
            (xpos_final, values_final, polyorder, NULL);
    
    cpl_size ny = image.size_y();
    cpl_size nx = image.size_x();
    if (pol_fit)
    {
        //Place in each column or row the polynomial fit
        T *      p_image      = image.get_data<T>();
        for (cpl_size j = 0; j< ny; ++j)
        {
            for (cpl_size i = 0; i< nx; ++i, ++p_image)
            {
                *p_image = cpl_polynomial_eval_1d(pol_fit, j, NULL);
            }
        }
    }
    else {
        cpl_msg_warning(cpl_func, "Invalid flat field dispersion fit ");
    }
    
    cpl_vector_delete(filtered_values);
    cpl_vector_delete(xpos);
    cpl_vector_delete(values_final);
    cpl_vector_delete(xpos_final);
    cpl_vector_delete(collapse_line);
    cpl_polynomial_delete(pol_fit);
    cpl_image_delete(collapse_median);
}

#endif
