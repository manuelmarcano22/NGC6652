/* 
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

#ifndef VIMOS_CALIBRATED_SLITS_H
#define VIMOS_CALIBRATED_SLITS_H

#include <vector>
#include <string>
#include "calibrated_slit.h"
#include "wavelength_calibration.h"
#include "grism_config.h"
#include "mosca_image.h"
#include "vimos_detected_slits.h"

namespace vimos
{

class calibrated_slits : public std::vector<mosca::calibrated_slit>
{
public:
    calibrated_slits(const vimos::detected_slits& det_slits,
                     const mosca::wavelength_calibration& wave_calib,
                     const mosca::grism_config& grism_cfg,
                     size_t ima_nx, size_t ima_ny);
    
private:
    
};

cpl_mask ** get_all_slits_valid_masks(const vimos::calibrated_slits& slits,
                                      mosca::axis disp_axis);

}

#endif
