/* $Id: vector_utils.h,v 1.1 2013-07-24 07:44:56 cgarcia Exp $
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

#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <vector>
#include "cpl_polynomial.h"
#include "gsl/gsl_bspline.h"

namespace mosca
{
template<typename T>
void vector_divide(std::vector<T>& input, 
                   std::vector<T>& input_err,
                   const std::vector<int>& dividend);

template<typename T>
void vector_smooth(std::vector<T>& input, 
                   std::vector<T>& input_err,
                   size_t smooth_size);

template<typename T>
void vector_smooth(std::vector<T>& input, 
                   size_t smooth_size);

//TODO: Document that smooth_size is actually half the size
template<typename T>
void vector_smooth(std::vector<T>& input, 
                   const std::vector<bool>& mask,
                   size_t smooth_size);

//TODO: Generalize this to cubicspline_fit. Fit X vs Y. Currently it is a particular
//case in which the X points are asumed to equispaced.
class vector_cubicspline
{
public:
    
    vector_cubicspline();
    
    ~vector_cubicspline();

    /* This will fit a spline fit to the points in X,Y space given by:
     * X: a regular spaced vector from 0 to input.size() -1
     * Y: the values in parameter input
     * This won't work if the input values represent a non-linear trend.
     * It modifies the yval vector with the fit.
     */
    template<typename T>
    void fit(std::vector<T>& yval, size_t& nknots);
    
    /* This will fit a spline fit to the points in X,Y space given by:
     * X: a regular spaced vector from 0 to input.size() -1
     * Y: the values in parameter input
     * This won't work if the input values represent a non-linear trend.
     * It modifies the yval vector with the fit.
     */
    template<typename T>
    void fit(std::vector<T>& yval, const std::vector<bool>& mask, size_t& nknots);
    
    /* This will fit a spline fit to the points in X,Y space given by
     * xval, yval. 
     * It modifies the yval vector with the fit.
     * The nknots are set uniformly in the x_min_knot, x_max_knot range.
     * If they are the same, the minimum and maximum values of xval are used.
     * The fit cannot be evaluated outside the knots range. 
     * If xval has values outside the knots range, yval is set to 0 
     * The mask specifies which vector values are valid (mask == 1)
     */
    template<typename T>
    void fit(const std::vector<T>& xval, std::vector<T>& yval,
             const std::vector<bool>& mask, 
             size_t& nknots,
             double x_min_knot = 0, double x_max_knot = 0);

    /* This evaluates the fitting at point xi. Take into account that if
     * the fit version with just an yval vector is used, then this xi is
     * in the scale of 0 to input.size() -1
     */
    double eval(double xi) const;
    
private:
    
    void m_clear_fit();
    
    gsl_bspline_workspace * m_bspline_workspace;

    gsl_matrix * m_covar;

    gsl_vector *m_spline_coeffs;
    
    gsl_vector *m_basis;
    
    double m_x_min_knot;

    double m_x_max_knot;

};

class vector_polynomial
{
public:

    vector_polynomial();
    
    ~vector_polynomial();
    
    /* This will fit a polynomial fit to the points in X,Y space given by:
     * X: a regular spaced vector from 0 to yval.size() -1
     * Y: the values in parameter input
     * This won't work if the input values represent a non-linear trend
     * It modifies the yval vector with the fit.
     */
    template<typename T>
    void fit(std::vector<T>& yval, size_t& polyorder);

    /* This will fit a polynomial fit to the points in X,Y space given by:
     * X: a regular spaced vector from 0 to yval.size() -1
     * Y: the values in parameter input
     * This won't work if the input values represent a non-linear trend
     * It modifies the yval vector with the fit.
     * The points where mask is not true are excluded.
     */
    template<typename T>
    void fit(std::vector<T>& yval, 
             const std::vector<bool>& mask, size_t& polyorder);

    /* This will fit a polynomial fit to the points in X,Y space given by
     * xval, yval. 
     * It modifies the yval vector with the fit.
     * TODO: Document that polyorder is modified
     */
    template<typename T>
    void fit(const std::vector<T>& xval, std::vector<T>& yval, 
             const std::vector<bool>& mask, 
             size_t& polyorder);

    /* This evaluates the fitting at point xi. Take into account that if
     * the fit version with just an yval vector is used, then this xi is
     * in the scale of 0 to input.size() -1
     */
    double eval(double xi) const;
    
private:
    
    void m_clear_fit();

    cpl_polynomial * m_pol_fit;

};
}


#include "vector_utils.tcc"

#endif
