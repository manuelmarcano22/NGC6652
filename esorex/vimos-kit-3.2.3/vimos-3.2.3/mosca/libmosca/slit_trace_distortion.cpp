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
#include <cmath>
#include "slit_trace_distortion.h"

namespace mosca 
{

slit_trace_distortion::slit_trace_distortion(cpl_table * polytraces) 
{
    m_polytraces = cpl_table_duplicate(polytraces);
}

slit_trace_distortion::slit_trace_distortion() :
    m_polytraces(NULL)  
{
}

slit_trace_distortion::slit_trace_distortion(slit_trace_distortion& rhs)
{
    m_polytraces = cpl_table_duplicate(rhs.m_polytraces);
}

slit_trace_distortion::~slit_trace_distortion()
{
    if(m_polytraces  != NULL)
        cpl_table_delete(m_polytraces );
}

//TODO: This is not get verified
cpl_image * slit_trace_distortion::calibrate_spatial
(cpl_image * image, cpl_table *slits, double reference, 
 double start_wavelength, double end_wavelength, double dispersion)
{
    if (image == NULL) 
        return NULL;

    if (dispersion <= 0.0) 
        return NULL;

    if (end_wavelength - start_wavelength < dispersion) 
        return NULL;
    
    
    cpl_image * spatial_calibrated;
    spatial_calibrated = m_calibrate_spatial
            (image, slits, m_polytraces, reference, 
             start_wavelength, end_wavelength, dispersion);

    return spatial_calibrated;
}


bool slit_trace_distortion::to_distorted(double spa_coord_undistorted, 
                                         double disp_coord,
                                         double &spa_coord_distorted,
                                         cpl_table *slits)
{
    return m_to_distorted(spa_coord_undistorted, disp_coord, spa_coord_distorted,
            slits, m_polytraces);
}

bool slit_trace_distortion::to_undistorted(double spa_coord_distorted, 
                                           double disp_coord,
                                           double &spa_coord_undistorted,
                                           cpl_table *slits)
{
    return m_to_undistorted(spa_coord_undistorted, disp_coord, spa_coord_distorted,
            slits, m_polytraces);
}

} /* namespace mosca */
