/* $Id: detected_slit.cpp,v 1.5 2013-08-21 15:13:55 cgarcia Exp $
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
 * $Date: 2013-08-21 15:13:55 $
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#include <stdexcept>
#include <iostream>
#include <cmath>
#include "calibrated_slit.h"

mosca::calibrated_slit::calibrated_slit
(const detected_slit& det_slit,
 const wavelength_calibration& wave_calib, 
 const grism_config& grism_cfg,
 size_t ima_nx, size_t ima_ny) :
 detected_slit(det_slit), m_wave_calib(wave_calib), m_grism_config(grism_cfg),
 m_ima_nx(ima_nx), m_ima_ny(ima_ny)
{
}

mosca::calibrated_slit::calibrated_slit(const calibrated_slit& rhs) :
detected_slit(rhs), m_wave_calib(rhs.m_wave_calib), 
m_grism_config(rhs.m_grism_config),
m_ima_nx(rhs.m_ima_nx), m_ima_ny(rhs.m_ima_ny)
{
    
}

mosca::calibrated_slit::~calibrated_slit()
{
}

cpl_mask * mosca::calibrated_slit::get_mask_valid(mosca::axis disp_axis) const
{
    
    cpl_mask * slit_mask = cpl_mask_new(m_ima_nx, m_ima_ny);
    
    double start_wave = m_grism_config.start_wave();
    double end_wave = m_grism_config.end_wave();

    int disp_bottom, disp_top, spa_bottom, spa_top;
    get_extent_pix(disp_bottom, spa_bottom, disp_top, spa_top);
    //TODO:: Probably these two lines has to be done at the time of slit creation, to make sure that only slits within the image limits can be created.
    if(spa_bottom < 1)
        spa_bottom = 1;
    if(disp_axis == mosca::X_AXIS)
    {
        if(spa_top > (int)m_ima_ny)
             spa_top = m_ima_ny;
    }
    else
    {
        if(spa_top > (int)m_ima_nx)
             spa_top = m_ima_nx;
    }

    for(cpl_size i_disp = disp_bottom; i_disp <= disp_top; i_disp++)
    {
        for(cpl_size i_spa = spa_bottom; i_spa <= spa_top; i_spa++)
        {
            if(within_trace((double)(i_disp), (double)(i_spa)))
            {

                double spatial_corrected = spatial_correct
                        ((double)(i_disp),
                                (double)(i_spa));
                double wavelength = m_wave_calib.get_wave
                        (spatial_corrected, (double)i_disp);
                if(wavelength > start_wave && wavelength < end_wave)
                {
                    if(disp_axis == mosca::X_AXIS)
                        cpl_mask_set(slit_mask, i_disp, i_spa, CPL_BINARY_1);
                    else 
                        cpl_mask_set(slit_mask, i_spa, i_disp, CPL_BINARY_1);
                }
            }
        }
    }
    
    return slit_mask;
}

bool mosca::calibrated_slit::has_valid_wavecal() const
{
    int disp_bottom, disp_top, spa_bottom, spa_top;
    get_extent_pix(disp_bottom, spa_bottom, disp_top, spa_top);
    if(spa_bottom < 1)
        spa_bottom = 1;
    for(cpl_size i_disp = disp_bottom; i_disp <= disp_top; i_disp++)
    {
        for(cpl_size i_spa = spa_bottom; i_spa <= spa_top; i_spa++)
        {
            double spatial_corrected = spatial_correct
                    ((double)(i_disp),
                     (double)(i_spa));
            if(m_wave_calib.has_valid_cal(spatial_corrected))
                return true;
        }
    }
    return false;
}
