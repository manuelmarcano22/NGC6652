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

#ifndef GLOBAL_DISTORTION_H
#define GLOBAL_DISTORTION_H

#include "spatial_distortion.h"

namespace mosca 
{

class global_distortion : public spatial_distortion
{
public:
    
    global_distortion(cpl_table * global_dist);
    
    global_distortion();
    
    global_distortion(global_distortion& rhs);
    
    ~global_distortion();
    
    cpl_image * calibrate_spatial(cpl_image * image, cpl_table *slits, 
                                  double reference, 
                                  double start_wavelength, double end_wavelength, 
                                  double dispersion);
    
    bool to_distorted(double spa_coord_undistorted, 
                      double disp_coord,
                      double &spa_coord_distorted,
                      cpl_table *slits);

    bool to_undistorted(double spa_coord_distorted, 
                        double disp_coord,
                        double &spa_coord_undistorted,
                        cpl_table *slits);

private:
    
    cpl_table * m_create_curv_coeff_table(cpl_table *slits);
    
    cpl_polynomial * m_read_polynomial_row(cpl_size row);
    
    cpl_table * m_global_dist;
};


} /* namespace mosca */
#endif /* GLOBAL_DISTORTION_H */
