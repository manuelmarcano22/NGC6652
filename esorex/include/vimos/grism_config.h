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

#ifndef GRISM_CONFIG_H
#define GRISM_CONFIG_H

namespace mosca 
{

class grism_config
{
public:
    
    grism_config(double nominal_dispersion,
                 double start_wave, 
                 double end_wave,
                 double wave_ref);
    
    grism_config();
    
    ~grism_config();
    
    double nominal_dispersion() const;

    double start_wave() const;
    
    double end_wave() const;
    
    double wave_ref() const;
    
protected:
    
    double m_nominal_dispersion;
    
    double m_start_wave;

    double m_end_wave;

    double m_wave_ref;
};


} /* namespace mosca */
#endif /* GRISM_CONFIG_H */
