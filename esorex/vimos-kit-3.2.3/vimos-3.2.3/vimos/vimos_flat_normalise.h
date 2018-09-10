/* $Id: vimos_flat_normalise.h,v 1.3 2013-09-09 12:14:12 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-09-09 12:14:12 $
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VIMOS_FLAT_NORMALISE_H
#define VIMOS_FLAT_NORMALISE_H

#include <cpl.h>
#include "mosca_image.h"
#include "wavelength_calibration.h"
#include "vimos_detected_slits.h"
#include "calibrated_slit.h"    

namespace vimos
{

class flat_normaliser
{
    
public:

    flat_normaliser();
    ~flat_normaliser();
    
    int mos_normalise(mosca::image& flat, 
                      const mosca::wavelength_calibration& wave_cal,
                      cpl_image *spatial, 
                      const std::vector<mosca::calibrated_slit>& calib_slits,
                      cpl_table *slits, cpl_table *polytraces, 
                      double blue, double red, 
                      double dispersion,
                      int spa_smooth_radius, int disp_smooth_radius, 
                      int spa_fit_polyorder, int disp_fit_nknots,
                      double fit_threshold);
    
    const mosca::image& get_normalisation_image() const;
    
    const std::vector<std::vector<float> >& get_wave_profiles() const;
    
    std::vector<float> get_wave_profiles_norm(double mflat_exptime,
                                              const std::vector<float>& slit_widths,
                                              const std::vector<float>& slit_lengths) const;
    
    cpl_image * get_wave_profiles_im() const;

    cpl_image * get_wave_profiles_im_mapped(const vimos::detected_slits& det_slits, 
                                            const mosca::wavelength_calibration& wave_cal,
                                            double firstLambda, 
                                            double lastLambda, 
                                            double dispersion) const;

    static int get_middle_slit_valid_calib
       (const mosca::wavelength_calibration& wave_cal, 
        int slit_end_pos, int slit_begin_pos);

private:
    
    mosca::image m_normalisation_image;
    
    std::vector<std::vector<float> > m_wave_profiles;
    
    std::vector<float> m_wave_profiles_norm;
    
};

}

#endif   /* VIMOS_FLAT_NORMALISE_H */
