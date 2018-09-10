/* $Id: detected_slit.h,v 1.4 2013-10-24 08:52:49 cgarcia Exp $
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
 * $Date: 2013-10-24 08:52:49 $
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifndef CALIBRATED_SLIT_H
#define CALIBRATED_SLIT_H

#include <vector>
#include "cpl_mask.h"
#include "detected_slit.h"
#include "wavelength_calibration.h"
#include "grism_config.h"
#include "mosca_image.h"

namespace mosca
{

/**
 * This class represents a slit that has been wavelength calibrated
 */
class calibrated_slit : public detected_slit
{
public:
    /* 
     * Constructor from a detected slit plus a wavelength calibration
     * @param det_slit    The detected list
     * @param wave_calib  The wavelength calibration
     * TODO: Decide whether all this should go in the constructor
     * (maybe grism_cfg, ima_nx, ima_ny should just be arguments of 
     * get_valid_region) 
     */
    calibrated_slit(const detected_slit& det_slit,
                    const wavelength_calibration& wave_calib,
                    const grism_config& grism_cfg,
                    size_t ima_nx, size_t ima_ny);
    
    calibrated_slit();

    calibrated_slit(const calibrated_slit& rhs);

    virtual ~calibrated_slit();
    
    cpl_mask * get_mask_valid(axis disp_axis) const;

    bool has_valid_wavecal() const;

private:

    //The wavelength calibration
    wavelength_calibration m_wave_calib;
    
    //The grism configuration 
    grism_config m_grism_config;    
    
    //Size of the full image in X
    size_t m_ima_nx;
    
    //Size of the full image in Y
    size_t m_ima_ny;
    
};

}

#endif
