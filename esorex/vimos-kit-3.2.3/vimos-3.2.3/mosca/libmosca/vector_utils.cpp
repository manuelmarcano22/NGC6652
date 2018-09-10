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


#include "vector_utils.h"
#include "gsl/gsl_multifit.h"

mosca::vector_cubicspline::vector_cubicspline() :
m_bspline_workspace(NULL), m_covar(NULL), m_spline_coeffs(NULL), m_basis(NULL),
m_x_min_knot(0.), m_x_max_knot(0.)
{
}

mosca::vector_cubicspline::~vector_cubicspline() 
{
    if(m_bspline_workspace != NULL)
        m_clear_fit();
}

void mosca::vector_cubicspline::m_clear_fit()
{
    gsl_bspline_free(m_bspline_workspace);
    gsl_matrix_free(m_covar);
    gsl_vector_free(m_spline_coeffs);
    gsl_vector_free(m_basis);
}

double mosca::vector_cubicspline::eval(double xi) const
{
    double valfit = 0, valfiterr = 0;
    if(xi > m_x_max_knot || xi < m_x_min_knot)
        throw std::domain_error("evaluating spline outside its domain");
    if(m_bspline_workspace != NULL)
    {
        gsl_bspline_eval(xi, m_basis, m_bspline_workspace);
        gsl_multifit_linear_est(m_basis, m_spline_coeffs, m_covar, 
                                &valfit, &valfiterr);
    }
    return valfit;
}

mosca::vector_polynomial::vector_polynomial() :
m_pol_fit(NULL)
{
    
}

mosca::vector_polynomial::~vector_polynomial() 
{
    if(m_pol_fit != NULL)
        m_clear_fit();
}

void mosca::vector_polynomial::m_clear_fit()
{
    cpl_polynomial_delete(m_pol_fit);
}

double mosca::vector_polynomial::eval(double xi) const
{
    double valfit = 0;
    if(m_pol_fit != NULL)
        valfit = cpl_polynomial_eval_1d(m_pol_fit, xi, NULL);
    return valfit;
}
