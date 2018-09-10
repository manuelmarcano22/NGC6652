/*
 * This file is part of the FORS Data Reduction Pipeline
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

/*
 * image_normalisation.h
 *
 *  Created on: 2014 3 28
 *      Author: cgarcia
 */

#ifndef IMAGE_NORMALISATION_H_
#define IMAGE_NORMALISATION_H_

#include <vector>
#include "mosca_image.h"

namespace mosca {

template<typename T>
mosca::image image_normalise
(mosca::image& image,
 mosca::image& slit_image_weight,
 int spa_smooth_radius, int disp_smooth_radius, 
 int spa_fit_polyorder, int disp_fit_nknots, double fit_threshold,
 std::vector<T>& slit_spa_profile, std::vector<T>& slit_disp_profile);

class no_flux_exception :  public std::exception
{
      virtual const char* what() const throw();
};


} /* namespace mosca */


#include "image_normalisation.tcc" 

#endif /* IMAGE_NORMALISATION_H_ */
