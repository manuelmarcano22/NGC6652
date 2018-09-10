/*
 * This file is part of the FORS Data Reduction Pipeline
 * Copyright (C) 2002-2010 European Southern Observatory
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

/*
 * spectrum.cpp
 *
 *  Created on: 2014 3 28
 *      Author: cgarcia
 */

#include <spectrum.h>
#include <cmath>
#include <stdexcept>

mosca::spectrum::spectrum(const cpl_image * spec, double start_wave,
                          double step_wave) :
 m_gsl_acc(NULL), m_gsl_interp(NULL)
{
    if(cpl_image_get_size_y(spec) != 1)
        throw std::invalid_argument("Only images with NY=1 supported");
    
    cpl_image * spec_double = cpl_image_cast(spec, CPL_TYPE_DOUBLE);
    cpl_size spec_size = cpl_image_get_size_x(spec);
    m_flux.insert(m_flux.end(), cpl_image_get_data_double(spec_double), 
                  cpl_image_get_data_double(spec_double) + spec_size);
    
    for(size_t i = 0; i< m_flux.size(); ++i)
        m_wave.push_back(start_wave + i * step_wave);

    cpl_image_delete(spec_double);
}

mosca::spectrum::spectrum(std::vector<double>& flux_vec, 
                          std::vector<double>& wave_vec) :
  m_flux(flux_vec), m_wave(wave_vec), m_gsl_acc(NULL), m_gsl_interp(NULL)
{
    if(flux_vec.size() != wave_vec.size())
        throw std::invalid_argument("Vectors do not have the same size");
}

mosca::spectrum::spectrum(const spectrum& other) :
m_flux(other.m_flux), m_wave(other.m_wave), m_gsl_acc(NULL), m_gsl_interp(NULL)
{
}

//TODO: Create a copy operator (destructor is not trivial)

mosca::spectrum::spectrum() :
m_flux(), m_wave(), m_gsl_acc(NULL), m_gsl_interp(NULL)
{
}

mosca::spectrum::~spectrum()
{
    if(m_gsl_interp != NULL)
    {
        gsl_interp_free (m_gsl_interp);
        gsl_interp_accel_free (m_gsl_acc);
    }
}

std::vector<double> mosca::spectrum::flux() const
{
    return m_flux;
}

std::vector<double> mosca::spectrum::wave() const
{
    return m_wave;
}

double mosca::spectrum::integrate
(double start_wave, double end_wave, bool ignore_neg_flux, 
 float ignore_threshold) const
{
    //Get only the positive flux if requested (and not done before)
    if(ignore_neg_flux && m_wave_nonzero.size() == 0)
        m_create_filtered_flux();
        
    //Get the pointers
    cpl_size npts;
    const double * wave_p;
    const double * flux_p;
    if(ignore_neg_flux)
    {
        wave_p = &(m_wave_nonzero[0]);
        flux_p = &(m_flux_nonzero[0]);
        npts = m_flux_nonzero.size();
    }
    else
    {
        wave_p = &(m_wave[0]);
        flux_p = &(m_flux[0]);
        npts = m_flux.size();
    }

    //Get the integration limits
    double start_int = std::max(start_wave, wave_p[0]);
    double end_int   = std::min(end_wave, wave_p[npts-1]);

    if(start_int >= end_int)
        return 0;

    //If the final integration range is less than ignore_threshold then ignore it
    if(ignore_neg_flux)
        if((end_int - start_int)/ (end_wave - start_wave) < ignore_threshold)
            return 0;
                    
    if(m_gsl_interp == NULL)
    {
        m_gsl_acc
            = gsl_interp_accel_alloc ();
        m_gsl_interp
            = gsl_interp_alloc (gsl_interp_linear, npts);
        gsl_interp_init (m_gsl_interp, wave_p, flux_p, npts);
    }
    
    double integrated_flux = 
            gsl_interp_eval_integ(m_gsl_interp, wave_p, flux_p, 
                                  start_int, end_int, m_gsl_acc);
    
    //Correct for the integration interval
    if(ignore_neg_flux)
        integrated_flux *= (end_wave - start_wave) / (end_int - start_int);
    return integrated_flux;
}

void mosca::spectrum::m_create_filtered_flux() const
{
    m_wave_nonzero.resize(m_wave.size());
    m_flux_nonzero.resize(m_wave.size());
    size_t i_nonzero = 0;
    for(size_t i_bin = 0; i_bin < m_wave.size(); ++i_bin)
    {
        if(m_flux[i_bin] > 0)
        {
            m_wave_nonzero[i_nonzero] = m_wave[i_bin];
            m_flux_nonzero[i_nonzero] = m_flux[i_bin];
            i_nonzero++;
        }
    }
    m_wave_nonzero.resize(i_nonzero);
    m_flux_nonzero.resize(i_nonzero);
}

mosca::spectrum mosca::spectrum::rebin
(double start_wave, double end_wave, double step_wave)
{
    //create regular step in wavelength
    std::vector<double> wave_rebinned;
    size_t nwave = std::floor((end_wave - start_wave) / step_wave);
    
    for(size_t i = 0; i< nwave; ++i)
        wave_rebinned.push_back(start_wave + i * step_wave);

    //TODO: This has to be fixed
    std::vector<double> flux_rebinned(m_flux); 
    
    spectrum spec_rebinned(flux_rebinned, wave_rebinned);
    return spec_rebinned;
}
