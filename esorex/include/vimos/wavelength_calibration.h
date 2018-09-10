/* $Id: wavelength_calibration.h,v 1.2 2013/07/24 14:53:01 cgarcia Exp $
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
 * $Date: 2013/07/24 14:53:01 $
 * $Revision: 1.2 $
 * $Name:  $
 */

#ifndef WAVELENGTH_CALIBRATION_H
#define WAVELENGTH_CALIBRATION_H

#include <vector>
#include "cpl_polynomial.h"
#include "cpl_table.h"

namespace mosca
{


/**
 * This class the wavelength calibration for all the rows in a rectified
 * image (so the rows correspond to rectified spatial values). 
 */
class wavelength_calibration
{
public:
    wavelength_calibration();
    
    wavelength_calibration(const std::string& fits_disp_coeff, double refwave);

    wavelength_calibration(const cpl_table * idscoeff, double refwave);

    wavelength_calibration(const wavelength_calibration& rhs);

    ~wavelength_calibration();
  
    void from_idscoeff(const cpl_table * idscoeff, double refwave);
    
    double get_pixel(double spatial_corrected_pos, double wavelength) const;

    double get_wave(double spatial_corrected_pos, double dispersion_pos) const;
    
    void min_max_wave(double& min_wave, double &max_wave, int size_dispersion,
                      int min_spa_row,
                      int max_spa_row) const;

    double mean_dispersion(double start_wave, 
                           double end_wave,
                           int min_spa_row,
                           int max_spa_row) const;
    
    bool has_valid_cal(double spatial_corrected_pos) const;

    bool is_monotonical(size_t spa_row, double start_wave, double end_wave, 
                        double dispersion) const;

    double get_refwave() const;

    
private:
    
    std::vector<cpl_polynomial *> m_wave_coeff;

    std::vector<int> m_nlines;
    
    double m_refwave;
};

}
#endif
