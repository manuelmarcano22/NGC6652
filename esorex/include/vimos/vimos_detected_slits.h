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

#ifndef VIMOS_DETECTED_SLITS_H
#define VIMOS_DETECTED_SLITS_H

#include <vector>
#include <string>
#include "detected_slit.h"
#include "wavelength_calibration.h"
#include "grism_config.h"
#include "mosca_image.h"

namespace vimos
{

typedef std::vector<mosca::detected_slit> detected_slits;

detected_slits detected_slits_load_fits(const std::string& fitsfile_slit_loc,
                                        const std::string& fitsfile_curv_coeff,
                                        int image_size_disp);

detected_slits detected_slits_from_tables(cpl_table * slit_loc,
                                          cpl_table * curv_coeff,
                                          int image_size_disp);

}

#endif
