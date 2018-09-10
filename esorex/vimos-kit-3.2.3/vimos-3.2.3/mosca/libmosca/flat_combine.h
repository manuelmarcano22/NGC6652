/* $Id: flat_combine.h,v 1.10 2013-08-13 13:12:32 cgarcia Exp $
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
 * $Date: 2013-08-13 13:12:32 $
 * $Revision: 1.10 $
 * $Name: not supported by cvs2svn $
 */

#ifndef FLAT_COMBINE_H
#define FLAT_COMBINE_H

#include <memory>
#include "calibrated_slit.h"
#include "wavelength_calibration.h"
#include "mosca_image.h"
#include "reduce_method.h"

namespace mosca
{
template<typename T, typename Iter, typename CombineMethod>
std::auto_ptr<mosca::image> 
flat_combine(Iter image_start, Iter image_end, 
             const std::vector<mosca::calibrated_slit>& slits,
             const mosca::wavelength_calibration& wave_cal,
             const mosca::grism_config& grism_cfg,
             size_t smooth_size,
             CombineMethod comb_method= mosca::reduce_mean());

template<typename T, typename CombineMethod>
std::auto_ptr<mosca::image> 
flat_combine(std::vector<mosca::image>& image_list, 
             const std::vector<mosca::calibrated_slit>& slits,
             const mosca::wavelength_calibration& wave_cal,
             const mosca::grism_config& grism_cfg,
             size_t smooth_size, CombineMethod com_method);
    
}


#include "flat_combine.tcc"

#endif
