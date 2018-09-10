/* $Id: vector_utils.tcc,v 1.1 2013-07-24 07:44:56 cgarcia Exp $
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
 * $Date: 2013-07-24 07:44:56 $
 * $Revision: 1.1 $
 * $Name: not supported by cvs2svn $
 */


#ifndef VECTOR_UTILS_TCC
#define VECTOR_UTILS_TCC

#include <stdexcept>
#include <algorithm>
#include "vector_utils.h"
#include "gsl/gsl_multifit.h"
#include "cpl_vector.h"

/**
 * @brief
 *   Collapse a list of images in a container with median 
 *
 * @return The mean image
 *
 */
template<typename T>
void mosca::vector_divide(std::vector<T>& input,
                          std::vector<T>& input_err,
                          const std::vector<int>& dividend)
{
    if((input.size() != input_err.size()) || (input.size() != dividend.size()))
        throw std::invalid_argument("Vector sizes do not match");

    for(size_t i = 0 ; i < input.size(); ++i)
    {
        input[i] = input[i] / dividend[i];
        input_err[i] = input_err[i] / dividend[i];
    }    

}


/* TODO: It is wrong, it uses only the forward elements.
 * Use some cpl stuff here 
 */
template<typename T>
void mosca::vector_smooth(std::vector<T>& input,
                          std::vector<T>& input_err,
                          size_t smooth_size)
{
    if(input.size() != input_err.size())
        throw std::invalid_argument("Vector sizes do not match");

    mosca::vector_smooth(input, smooth_size);
    mosca::vector_smooth(input_err, smooth_size);
}

template<typename T>
void mosca::vector_smooth(std::vector<T>& input,
                          size_t smooth_size)
{
    std::vector<bool> mask(input.size(), true);
    mosca::vector_smooth(input, mask, smooth_size);
}

/* TODO: It is wrong, it uses only the forward elements.
 * Use some cpl stuff here 
 */
template<typename T>
void mosca::vector_smooth(std::vector<T>& input,
                          const std::vector<bool>& mask,
                          size_t smooth_size)
{
    if(smooth_size >= input.size())
        throw std::invalid_argument("Smooth size too large");

    if(mask.size() != input.size())
        throw std::invalid_argument("Vector and mask size mismatch");

    cpl_vector * line = cpl_vector_new(input.size());
    cpl_size n_unmasked = 0;
    for(size_t i = 0; i<input.size(); i++)
    {
        if(mask[i])
        {
            cpl_vector_set(line, n_unmasked, input[i]);
            n_unmasked++;
        }
    }
    cpl_vector_set_size(line, n_unmasked);
    
    //We smooth the line-image
    cpl_vector * smooth_line = 
            cpl_vector_filter_median_create(line, smooth_size);

    int i_unmasked = 0;
    for(size_t i = 0; i<input.size(); i++)
    {
        if(mask[i])
        {
            input[i] = T(cpl_vector_get(smooth_line, i_unmasked));
            i_unmasked++;
        }
    }
    
    //Cleanup
    cpl_vector_delete(smooth_line);
    cpl_vector_delete(line);
}

template<typename T>
void mosca::vector_cubicspline::fit(std::vector<T>& yval, 
                                    size_t& nknots)
{
    std::vector<T> xval;
    for (size_t i = 0; i < yval.size(); ++i)
        xval.push_back(T(i));
    m_x_min_knot = 0;
    m_x_max_knot = yval.size() - 1;
    std::vector<bool> mask(yval.size(), true);
    fit(xval, yval, mask, nknots);
}

template<typename T>
void mosca::vector_cubicspline::fit(std::vector<T>& yval, 
                                    const std::vector<bool>& mask, 
                                    size_t& nknots)
{
    std::vector<T> xval;
    for (size_t i = 0; i < yval.size(); ++i)
        xval.push_back(T(i));
    m_x_min_knot = 0;
    m_x_max_knot = yval.size() - 1;
    fit(xval, yval, mask, nknots);
}

template<typename T>
void mosca::vector_cubicspline::fit(const std::vector<T>& xval, std::vector<T>& yval, 
                                    const std::vector<bool>& mask, 
                                    size_t& nknots,
                                    double x_min_knot, double x_max_knot)
{

    if(xval.size() != yval.size())
        throw std::invalid_argument("xval and yval sizes do not match");

    if(nknots <= 1)
        throw std::invalid_argument("number of knots must be at least 2");
    
    //Allocate result
    size_t     nval = yval.size();

    //This is valid for cubic splines, if not nbreak = ncoeffs + 2 - k, 
    //where k is the degree of the spline
    int ncoeffs = nknots + 2;

    //Get the nknots range
    if(x_min_knot == x_max_knot)
    {
        m_x_min_knot = *std::min_element(xval.begin(), xval.end());
        m_x_max_knot = *std::max_element(xval.begin(), xval.end());
    }
    else
    {
        m_x_min_knot = x_min_knot;
        m_x_max_knot = x_max_knot;
    }
    
    /* Create a "mask" of pixels not to use by combining the input mask */
    std::vector<bool> input_mask(mask);
    for (size_t i = 0; i < nval; ++i)
        if(xval[i] < m_x_min_knot || xval[i] > m_x_max_knot)
            input_mask[i] = false;
    int nfit = std::count(mask.begin(), mask.end(), true);

    /* Throw if the fit is going to fail */
    if(nfit < ncoeffs)
    {
        ncoeffs = nfit;
        nknots = ncoeffs - 2;
    }
    if(nfit < 3)
        throw std::length_error("Number of fitting points too small");
    
    /* Deallocate fit if not empty */
    if(m_bspline_workspace != NULL)
        m_clear_fit();
    
    /* allocate a cubic bspline workspace (k = 4) */
    gsl_matrix * X;
    m_bspline_workspace = gsl_bspline_alloc(4, nknots);
    m_basis = gsl_vector_alloc(ncoeffs);
    X = gsl_matrix_alloc(nfit, ncoeffs);
    
    /* allocate objects for the fitting */
    gsl_vector * y_fit;
    gsl_vector * weigth;
    gsl_multifit_linear_workspace * mfit_wspc;
    double chisq;
    y_fit = gsl_vector_alloc(nfit);
    weigth = gsl_vector_alloc(nfit);
    mfit_wspc = gsl_multifit_linear_alloc(nfit, ncoeffs);
    m_spline_coeffs = gsl_vector_alloc(ncoeffs);
    m_covar = gsl_matrix_alloc(ncoeffs, ncoeffs);
    
    /* use uniform breakpoints on x_min, x_max, which is the range of fitting */
    gsl_bspline_knots_uniform(m_x_min_knot, m_x_max_knot,
                              m_bspline_workspace);
    

    /* construct the fit matrix X */
    for (size_t i = 0, ifit = 0; i < nval; ++i)
    {
        double xi = xval[i];
        double yi = yval[i];

        if(input_mask[i])
        {
            /* Fill the vector to fit */
            gsl_vector_set(y_fit, ifit, yi);
            gsl_vector_set(weigth, ifit, 1.);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(xi, m_basis, m_bspline_workspace);

            /* fill in row i of X */
            for (int j = 0; j < ncoeffs; ++j)
            {
                double Bj = gsl_vector_get(m_basis, j);
                gsl_matrix_set(X, ifit, j, Bj);
            }
            ifit++;
        }
    } 
    
    /* do the fit */
    gsl_multifit_wlinear(X, weigth, y_fit, m_spline_coeffs, m_covar, &chisq, mfit_wspc);
    
    /* output the fit */

    for(size_t i = 0; i < nval; i++)
    {
        double yi, yerr;
        if(xval[i] >= m_x_min_knot && xval[i] <= m_x_max_knot)
        {
            gsl_bspline_eval(xval[i], m_basis, m_bspline_workspace);
            gsl_multifit_linear_est(m_basis, m_spline_coeffs, m_covar, &yi, &yerr);
            yval[i] = T(yi);
        }
        else 
            yval[i] = 0;
   }

    gsl_vector_free(y_fit);
    gsl_vector_free(weigth);
    gsl_multifit_linear_free(mfit_wspc);
}

template<typename T>
void mosca::vector_polynomial::fit(std::vector<T>& yval, size_t& polyorder)
{
    std::vector<T> xval;
    for (size_t i = 0; i < yval.size(); ++i)
        xval.push_back(T(i));
    std::vector<bool> mask(yval.size(), true);
    fit(xval, yval, mask, polyorder);
}

template<typename T>
void mosca::vector_polynomial::fit(std::vector<T>& yval, 
                                   const std::vector<bool>& mask, 
                                   size_t& polyorder)
{
    std::vector<T> xval;
    for (size_t i = 0; i < yval.size(); ++i)
        xval.push_back(T(i));
    fit(xval, yval, mask, polyorder);
}

/**
 * Smooth a vector fitting a polynomial 
 * @param input       Vector to be smoothed
 * @param polyorder   the order of the polynomial to fit
 * 
 * This function will smooth a MOS image in the dispersion direction 
 * (here assumed to be the X axis) fitting each row with a polynomial.
 * TODO: merge this with image_smooth_fit_1d_pol_spa() and place it in MOSCA.
 * TODO: It is named image_smooth_fit_1d_pol but it is not yet generalized
 * to any direction, only to spatial.   
 */
template<typename T>
void mosca::vector_polynomial::fit(const std::vector<T>& xval, std::vector<T>& yval, 
                                   const std::vector<bool>& mask, 
                                   size_t& polyorder)
{
    if(xval.size() != yval.size() ||
       xval.size() != mask.size())
        throw std::invalid_argument("xval, yval and mask sizes do not match");
    
    //Allocate result
    size_t     nval = yval.size();

    int nfit = std::count(mask.begin(), mask.end(), true);

    cpl_vector * filtered_values = cpl_vector_new(nfit);
    cpl_vector * xpos = cpl_vector_new(nfit);
    for (size_t i = 0, ifit = 0; i < nval; ++i)
    {
        if(mask[i])
        {
            cpl_vector_set(filtered_values, (cpl_size)ifit, yval[i]);
            cpl_vector_set(xpos, ifit, xval[i]);
            ifit++;
        }
    }
 
    /* Throw if the fit is going to fail */
    if(cpl_vector_get_size(xpos) < polyorder + 1)
        polyorder = cpl_vector_get_size(xpos) - 1;
    if(cpl_vector_get_size(xpos) < 1)
        throw std::length_error("Number of fitting points too small");

    /* Deallocate fit if not empty */
    if(m_pol_fit != NULL)
        m_clear_fit();

    /* Allocate fit */
    m_pol_fit = cpl_polynomial_fit_1d_create
            (xpos, filtered_values, polyorder, NULL);
    
    if (m_pol_fit)
    {
        for(size_t i = 0; i < nval; i++)
            yval[i] = T(cpl_polynomial_eval_1d(m_pol_fit, xval[i], NULL));
    }
    else //TODO: Should this throw an exception instead?
        std::fill(yval.begin(), yval.end(), T()); 
    
    cpl_vector_delete(filtered_values);
    cpl_vector_delete(xpos);
}

#endif
