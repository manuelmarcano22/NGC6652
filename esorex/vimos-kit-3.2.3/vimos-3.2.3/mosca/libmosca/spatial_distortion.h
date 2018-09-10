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

#ifndef SPATIAL_DISTORTION_H
#define SPATIAL_DISTORTION_H

#include <utility>

namespace mosca 
{

class spatial_distortion
{
public:
    
    spatial_distortion();
    
    virtual ~spatial_distortion();
    
    //TODO: Do this with mosca::calibrated_slit 
    virtual cpl_image * calibrate_spatial(cpl_image * image, cpl_table *slits, 
                                          double reference, 
                                          double start_wavelength, double end_wavelength, 
                                          double dispersion) = 0;
    
    //TODO: Merge this with the spatial_correct functionality within detected_slit
    virtual bool to_distorted(double spa_coord_undistorted, 
                              double disp_coord,
                              double &spa_coord_distorted,
                              cpl_table *slits) = 0;

    virtual bool to_undistorted(double spa_coord_distorted, 
                                double disp_coord,
                                double &spa_coord_undistorted,
                                cpl_table *slits) = 0;

protected:
    
    cpl_image * m_calibrate_spatial
           (cpl_image * image, cpl_table * slits, cpl_table * polytraces, 
            double reference, double start_wavelength, double end_wavelength, 
            double dispersion);
    
    bool m_to_distorted
    (double spa_coord_undistorted, 
     double disp_coord,
     double &spa_coord_distorted,
     cpl_table *slits, cpl_table * polytraces);
    
    bool m_to_undistorted
    (double spa_coord_distorted, 
     double disp_coord,
     double &spa_coord_undistorted,
     cpl_table *slits, cpl_table * polytraces);

    bool m_get_curv_polynomials(cpl_table * polytraces, 
                                cpl_table * slits,
                                cpl_size i_slit, 
                                cpl_polynomial * polytop, 
                                cpl_polynomial * polybot);
};


} /* namespace mosca */
#endif /* SPATIAL_DISTORTION_H */
