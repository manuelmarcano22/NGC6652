/* $Id: wavelength_calibration.cpp,v 1.5 2013/08/07 15:47:01 cgarcia Exp $
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
 * $Date: 2013/08/07 15:47:01 $
 * $Revision: 1.5 $
 * $Name:  $
 */

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <cpl_table.h>
#include <cpl_msg.h>
#include <cpl_polynomial.h>
#include "wavelength_calibration.h"
#include "statistics.h"

mosca::wavelength_calibration::wavelength_calibration()
{   
}

mosca::wavelength_calibration::wavelength_calibration
(const std::string& fits_disp_coeff, double refwave) :
m_refwave(refwave)
{
    /* We read from extension 1 of slit location table */
    cpl_table* idscoeff = cpl_table_load(fits_disp_coeff.c_str(), 1, 1);

    from_idscoeff(idscoeff, refwave); 
   
    cpl_table_delete(idscoeff);
}

mosca::wavelength_calibration::wavelength_calibration
(const cpl_table * idscoeff, double refwave) :
m_refwave(refwave)
{
    from_idscoeff(idscoeff, refwave);
}

mosca::wavelength_calibration::wavelength_calibration
(const wavelength_calibration& rhs)
{
    m_refwave = rhs.m_refwave;

    std::vector<cpl_polynomial *>::const_iterator poly;
    for(poly = rhs.m_wave_coeff.begin(); poly != rhs.m_wave_coeff.end(); ++poly)
    {
        if(*poly != NULL)
            m_wave_coeff.push_back(cpl_polynomial_duplicate(*poly));
        else 
            m_wave_coeff.push_back(NULL);
    }
}

//TODO: Create a copy operator (copy constructor is not trivial)

void mosca::wavelength_calibration::from_idscoeff
(const cpl_table * idscoeff, double refwave)
{
    cpl_size ncol = cpl_table_get_ncol(idscoeff) - 2;

    for(cpl_size irow = 0 ; irow < cpl_table_get_nrow(idscoeff); irow++)
    {
        std::vector<double> pol_coeff;
        int null = 0;
        for(cpl_size idx_coeff = 0; idx_coeff < ncol; idx_coeff++)
        {
            std::ostringstream colname;
            colname<<std::left<<"c"<<idx_coeff;
            if(cpl_table_has_column(idscoeff, colname.str().c_str()))
            {
                pol_coeff.push_back(cpl_table_get_double
                        (idscoeff, colname.str().c_str(), irow, &null));
                if(null)
                    break;
            }
        }
        cpl_polynomial * poly = NULL;
        if(!null)
        {
            poly = cpl_polynomial_new(1);
            cpl_size idx;
            std::vector<double>::reverse_iterator coeff;
            for( idx = pol_coeff.size() - 1, coeff = pol_coeff.rbegin(); coeff != pol_coeff.rend();
                    ++coeff, --idx)
                cpl_polynomial_set_coeff(poly, &idx, *coeff);
        }
        m_wave_coeff.push_back(poly);

        m_nlines.push_back(cpl_table_get_int
                               (idscoeff, "nlines", irow, &null));
    }
    
    m_refwave = refwave;
}

mosca::wavelength_calibration::~wavelength_calibration()
{
    std::vector<cpl_polynomial *>::iterator poly;
    for(poly = m_wave_coeff.begin(); poly != m_wave_coeff.end(); ++poly)
    {
        if(*poly != NULL)
            cpl_polynomial_delete(*poly);
    }
}

double mosca::wavelength_calibration::get_pixel(double spatial_corrected_pos,
                                                double wavelength) const
{
    double pixel = -1;
    size_t row = (size_t)(spatial_corrected_pos);
    if(row >= m_wave_coeff.size())
        return pixel;
    
    cpl_polynomial * poly = m_wave_coeff[row];
    if(poly == NULL)
        return pixel;

    pixel = cpl_polynomial_eval_1d(poly, wavelength - m_refwave, NULL);
    return pixel;
}

//TODO:Starting at 0 or 1?
double mosca::wavelength_calibration::get_wave(double spatial_corrected_pos,
                                               double dispersion_pos) const
{
    double wavelength = -1; //Denotes and error. TODO
    /* TODO: I think this is wrong, the reference wavelength probably has to be added 
     * Already done, but needs to be checked */
    size_t row = (size_t)(spatial_corrected_pos);
    if(row >= m_wave_coeff.size())
        return wavelength;
    
    cpl_polynomial * poly = m_wave_coeff[row];
    if(poly == NULL)
        return wavelength;

    cpl_polynomial *  inv_poly = cpl_polynomial_duplicate(poly);

    cpl_size   zero_coeff = 0;
    double coeff = cpl_polynomial_get_coeff(inv_poly, &zero_coeff);
    cpl_polynomial_set_coeff(inv_poly, &zero_coeff, coeff - dispersion_pos);

    wavelength = 0;
    cpl_polynomial_solve_1d(inv_poly, wavelength, &wavelength, 1);
    if(cpl_error_get_code() == CPL_ERROR_DIVISION_BY_ZERO ||
       cpl_error_get_code() == CPL_ERROR_CONTINUE)
    {
        cpl_error_reset();
        cpl_polynomial_delete(inv_poly);
        return -1;
    }

    cpl_polynomial_delete(inv_poly);
    return wavelength + m_refwave;
}


void mosca::wavelength_calibration::min_max_wave(double& min_wave,
                                                 double& max_wave,
                                                 int size_dispersion,
                                                 int min_spa_row,
                                                 int max_spa_row) const
{
     
    std::vector<double> wave_pix_0;
    std::vector<double> wave_pix_size;
    /* TODO: Check this limits */
    for(int irow = min_spa_row ; irow < max_spa_row; ++irow)
    {
        cpl_polynomial * poly = m_wave_coeff[irow];
        if(poly != NULL)
        {
            double wave_0 = cpl_polynomial_eval_1d(poly, 0., NULL);
            double wave_size = cpl_polynomial_eval_1d(poly, (double)size_dispersion, NULL);
            wave_pix_0.push_back(wave_0);
            wave_pix_size.push_back(wave_size);
        }
    }
    if(wave_pix_0.size() != 0)
    {
        min_wave = *std::min_element(wave_pix_0.begin(), wave_pix_0.end());
        max_wave = *std::max_element(wave_pix_size.begin(), wave_pix_size.end());
    }
    else 
        min_wave = max_wave = 0;
}

double mosca::wavelength_calibration::mean_dispersion(double start_wave, 
                                                      double end_wave,
                                                      int min_spa_row,
                                                      int max_spa_row) const
{
    std::vector<double> pix_startwave;
    std::vector<double> pix_endwave;

    for(int irow = min_spa_row ; irow < max_spa_row; ++irow)
    {
        if(m_nlines[irow] != 0) //For some reason there are wavelength calibrations that are not NULL in the coeffs but have just zero lines.
        {
            double pix_start = get_pixel((double)irow, start_wave);
            double pix_end = get_pixel((double)irow, end_wave);
            if(pix_start != -1 )
                pix_startwave.push_back(pix_start);
            if(pix_end != -1)
                pix_endwave.push_back(pix_end);
        }
    }
    
    if(!pix_startwave.empty() && !pix_endwave.empty())
    {
        double mean_pix_start= mosca::mean(pix_startwave.begin(), pix_startwave.end());
        double mean_pix_end = mosca::mean(pix_endwave.begin(), pix_endwave.end());
        double pixels_span = std::fabs(mean_pix_end - mean_pix_start);
        
        double mean_disp = (end_wave - start_wave) / pixels_span;

        return mean_disp;
    }
    else 
        return 0;
}

bool mosca::wavelength_calibration::has_valid_cal
(double spatial_corrected_pos) const
{
    size_t row = (size_t)(spatial_corrected_pos);
    if(row >= m_wave_coeff.size())
        return false;
    cpl_polynomial * poly = m_wave_coeff[row];
    if(poly == NULL)
        return false;
    
    return true;
}

bool mosca::wavelength_calibration::is_monotonical
(size_t spa_row, double start_wave, double end_wave, double dispersion) const
{
    if(spa_row >= m_wave_coeff.size())
        return false;
    cpl_polynomial * poly = m_wave_coeff[spa_row];
    if(poly == NULL)
        return false;
    
    for(double wave = start_wave; wave <= end_wave; wave += dispersion)
    {
        double grad;
        cpl_polynomial_eval_1d(poly, wave - m_refwave, &grad);
        
        if(grad < 0)
            return false;
    }

    return true;
}

double mosca::wavelength_calibration::get_refwave() const
{
    return m_refwave;
}
