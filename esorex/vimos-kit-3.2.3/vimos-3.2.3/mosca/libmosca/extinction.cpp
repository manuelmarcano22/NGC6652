/*
 * This file is part of the MOSCA library
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
 * extinction.cpp
 *
 *  Created on: 2014 4 2
 *      Author: cgarcia
 */

#include <cmath>
#include <extinction.h>


mosca::extinction::extinction(const cpl_table * ext) :
 m_gsl_acc(NULL), m_gsl_interp(NULL)
{
    m_extinction = cpl_table_duplicate(ext);
    cpl_table_cast_column(m_extinction, "WAVE", "WAVE_D", CPL_TYPE_DOUBLE);
    cpl_table_cast_column(m_extinction, "EXTINCTION", 
                          "EXTINCTION_D", CPL_TYPE_DOUBLE);
}

mosca::extinction::extinction(const mosca::extinction& other) :
 m_gsl_acc(NULL), m_gsl_interp(NULL)
{
    m_extinction = cpl_table_duplicate(other.m_extinction);
}

mosca::extinction::extinction() :
m_extinction(NULL), m_gsl_acc(NULL), m_gsl_interp(NULL)
{
}

//TODO: Create a copy operator (copy destructor is not trivial)

mosca::extinction::~extinction()
{
    if(m_extinction != NULL)
        cpl_table_delete(m_extinction);
    if(m_gsl_interp != NULL)
    {
        gsl_interp_free (m_gsl_interp);
        gsl_interp_accel_free (m_gsl_acc);
    }
}

mosca::spectrum mosca::extinction::correct_spectrum
(const spectrum& spec, double airmass)
{
    std::vector<double> flux = spec.flux(); 
    std::vector<double> wave = spec.wave();
    
    for(size_t i = 0; i < flux.size(); ++i)
    {
        double wavelength = wave[i];
        double interp_ext = eval_at_wave(wavelength);
        double ext_factor = std::pow(10., interp_ext*0.4*airmass);
        double corrected_flux = flux[i]*ext_factor;
        flux[i] = corrected_flux;
    }
    
    mosca::spectrum spec_corrected(flux, wave);
    
    return spec_corrected;
}

double mosca::extinction::eval_at_wave(double wavelength)
{
    if(m_extinction == NULL)
        return 0.;

    cpl_size npts = cpl_table_get_nrow(m_extinction);
    double * wave = 
        cpl_table_get_data_double(m_extinction, "WAVE_D");
    double * ext = 
        cpl_table_get_data_double(m_extinction, "EXTINCTION_D");
    
    if(wavelength > wave[npts-1])
        return ext[npts-1];
    if(wavelength < wave[0])
        return ext[0];
                
    if(m_gsl_interp == NULL)
    {
        m_gsl_acc    = gsl_interp_accel_alloc ();
        m_gsl_interp = gsl_interp_alloc (gsl_interp_linear, npts);
        gsl_interp_init (m_gsl_interp,wave, ext, npts);
    }
    
    double ext_at_wave =
        gsl_interp_eval (m_gsl_interp, wave, ext, wavelength, m_gsl_acc);
    return ext_at_wave;
}

