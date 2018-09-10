/* $Id: detected_slit.cpp,v 1.5 2013-08-21 15:13:55 cgarcia Exp $
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
 * $Date: 2013-08-21 15:13:55 $
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#include <stdexcept>
#include <iostream>
#include <cmath>
#include "detected_slit.h"

mosca::detected_slit::detected_slit(int slit_ident, 
                                    double disp_bottom, double spa_bottom, 
                                    double disp_top,    double spa_top,
                                    int position_spatial_corrected,
                                    int length_spatial_corrected,
                                    const std::vector<double>& trace_bottom_pol_coeffs,                  
                                    const std::vector<double>& trace_top_pol_coeffs) :
 m_slit_id(slit_ident), m_disp_bottom(disp_bottom),
 m_spa_bottom(spa_bottom), m_disp_top(disp_top), m_spa_top(spa_top),
 m_position_spatial_corrected(position_spatial_corrected),
 m_length_spatial_corrected(length_spatial_corrected),
 m_trace_bottom_pol_coeffs(trace_bottom_pol_coeffs),
 m_trace_top_pol_coeffs(trace_top_pol_coeffs)
{
    
    m_trace_bottom_pol = cpl_polynomial_new(1);
    m_trace_top_pol =  cpl_polynomial_new(1);
    cpl_size idx;
    idx = 0;
    std::vector<double>::reverse_iterator coeff;
    for(idx = trace_bottom_pol_coeffs.size() - 1, coeff = m_trace_bottom_pol_coeffs.rbegin();
            coeff != m_trace_bottom_pol_coeffs.rend(); ++coeff, --idx)
        cpl_polynomial_set_coeff(m_trace_bottom_pol, &idx, *coeff);

    for(idx = trace_top_pol_coeffs.size() -1, coeff = m_trace_top_pol_coeffs.rbegin();
            coeff != m_trace_top_pol_coeffs.rend(); ++coeff, --idx)
        cpl_polynomial_set_coeff(m_trace_top_pol, &idx, *coeff);
}

mosca::detected_slit::detected_slit()
{
    m_slit_id     = 0;
    m_disp_bottom = 0;
    m_spa_bottom  = 0;
    m_disp_top    = 0;
    m_spa_top     = 0;
    m_position_spatial_corrected = 0;
    m_length_spatial_corrected   = 0;
    m_trace_bottom_pol = NULL;
    m_trace_top_pol = NULL;

}

mosca::detected_slit::detected_slit(const detected_slit& rhs)
{
    m_slit_id     = rhs.m_slit_id;
    m_disp_bottom = rhs.m_disp_bottom;
    m_spa_bottom  = rhs.m_spa_bottom;
    m_disp_top    = rhs.m_disp_top;
    m_spa_top     = rhs.m_spa_top;
    m_position_spatial_corrected = rhs.m_position_spatial_corrected;
    m_length_spatial_corrected = rhs.m_length_spatial_corrected;
    m_trace_bottom_pol = cpl_polynomial_duplicate(rhs.m_trace_bottom_pol);
    m_trace_top_pol = cpl_polynomial_duplicate(rhs.m_trace_top_pol);

}

mosca::detected_slit::~detected_slit()
{
    if(m_trace_bottom_pol != NULL)
        cpl_polynomial_delete(m_trace_bottom_pol);
    if(m_trace_top_pol != NULL)
        cpl_polynomial_delete(m_trace_top_pol);
}

void mosca::detected_slit::get_extent(double & disp_bottom, double & spa_bottom, 
                                      double & disp_top,    double & spa_top)
{
    disp_bottom = m_disp_bottom;
    spa_bottom  = m_spa_bottom;
    disp_top    = m_disp_top;
    spa_top     = m_spa_top;
}

void mosca::detected_slit::get_extent_pix(int& disp_bottom, int& spa_bottom, 
                                          int& disp_top,    int& spa_top) const
{
    if(m_disp_bottom < m_disp_top)
    {
        disp_bottom = (int)std::floor(m_disp_bottom);
        disp_top    = (int)std::ceil(m_disp_top);
    }
    else
    {
        disp_bottom = (int)std::ceil(m_disp_bottom);
        disp_top    = (int)std::floor(m_disp_top);

    }
    if(m_spa_bottom < m_spa_top)
    {
        spa_bottom  = (int)std::floor(m_spa_bottom);
        spa_top     = (int)std::ceil(m_spa_top);
    }
    else
    {
        spa_bottom  = (int)std::ceil(m_spa_bottom);
        spa_top     = (int)std::floor(m_spa_top);
    }
}

//TODO:Starting in 0 or 1?
bool mosca::detected_slit::within_trace(double dispersion_pos,
                                        double spatial_pos) const
{
    
    double trace_bottom = cpl_polynomial_eval_1d(m_trace_bottom_pol, 
                                                 dispersion_pos, NULL);
    double trace_top    = cpl_polynomial_eval_1d(m_trace_top_pol, 
                                                 dispersion_pos, NULL);
    
    if(spatial_pos >= trace_bottom && spatial_pos <= trace_top)
        return true;

    return false;
}

double mosca::detected_slit::spatial_correct(double dispersion_pos,
                                             double spatial_pos) const
{
    double trace_bottom = cpl_polynomial_eval_1d(m_trace_bottom_pol, 
                                                 dispersion_pos, NULL);
    double trace_top    = cpl_polynomial_eval_1d(m_trace_top_pol, 
                                                 dispersion_pos, NULL);
    //The -1 factor is there because originally the length of the spatial
    //corrected slit has a +1 added in function mos_spatial_calibration 
    double factor = (trace_top-trace_bottom)/(m_length_spatial_corrected - 1);
    
    double spa_corrected_from_top = (spatial_pos - trace_top) / factor;
    
    double spa_corrected = m_position_spatial_corrected + 
            m_length_spatial_corrected + spa_corrected_from_top;
    
    
    return spa_corrected;
}


int mosca::detected_slit::slit_id() const
{
    return m_slit_id;
}

int mosca::detected_slit::get_position_spatial_corrected() const
{
    return m_position_spatial_corrected;
}

int mosca::detected_slit::get_length_spatial_corrected() const
{
    return m_length_spatial_corrected;
}
