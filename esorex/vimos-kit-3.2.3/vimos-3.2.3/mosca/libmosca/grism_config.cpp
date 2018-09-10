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
 * fiera_ccd.cpp
 *
 *  Created on: 2013 11 25
 *      Author: cgarcia
 */

#include <cpl.h>
#include <stdexcept>
#include "grism_config.h"

namespace mosca 
{

grism_config::grism_config(double nominal_disp,
                           double start_wavelength, 
                           double end_wavelength,
                           double wavelength_ref) :
 m_nominal_dispersion(nominal_disp), 
 m_start_wave(start_wavelength), 
 m_end_wave(end_wavelength),
 m_wave_ref(wavelength_ref)
{
}

grism_config::grism_config():
m_nominal_dispersion(0.), m_start_wave(0.), m_end_wave(0.), m_wave_ref(0.)
{
}

grism_config::~grism_config()
{
}

double grism_config::nominal_dispersion() const
{
    return m_nominal_dispersion;
}

double grism_config::start_wave() const
{
    return m_start_wave;
}

double grism_config::end_wave() const
{
    return m_end_wave;
}

double grism_config::wave_ref() const
{
    return m_wave_ref;
}

} /* namespace mosca */
