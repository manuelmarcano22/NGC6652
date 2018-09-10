/* 
 * This file is part of the VIMOS pipeline
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


#include "cpl.h"
#include "vimos_calibrated_slits.h"


vimos::calibrated_slits::calibrated_slits
(const vimos::detected_slits& det_slits,
 const mosca::wavelength_calibration& wave_calib,
 const mosca::grism_config& grism_cfg,
 size_t ima_nx, size_t ima_ny) :
 std::vector<mosca::calibrated_slit>()
{

    size_t n_slits = det_slits.size();
    for(size_t i_slit = 0; i_slit < n_slits; i_slit++)
    {
        mosca::calibrated_slit calib_slit(det_slits[i_slit], wave_calib,
                                          grism_cfg, ima_nx, ima_ny);
        push_back(calib_slit);
    }

}


cpl_mask ** vimos::get_all_slits_valid_masks(const vimos::calibrated_slits& slits,
                                            mosca::axis disp_axis)
{

    /* We work on a slit per slit basis */
    size_t n_slits = slits.size();

    /* Get the masks of each of the slits */
    cpl_mask ** slit_masks = (cpl_mask **)cpl_malloc(n_slits * sizeof(cpl_mask*));
    for(size_t i_slit = 0; i_slit < n_slits; i_slit++)
    {
        slit_masks[i_slit] = slits[i_slit].get_mask_valid(disp_axis);
    }
    
    return slit_masks;
}
