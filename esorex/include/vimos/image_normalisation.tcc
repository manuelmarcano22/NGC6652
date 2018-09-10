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
 * image_normalisation.cpp
 *
 *  Created on: 2014 3 28
 *      Author: cgarcia
 */

#include <cpl.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <numeric>
#include <exception>
#include "image_normalisation.h"
#include "vector_utils.h"

const char* mosca::no_flux_exception::what() const throw()
{
    const char* ret =
            "The sum of all the flux contributions for the provided slit "
            "is zero, making normalisation not possible";
    return ret;
}

template<typename  T>
mosca::image mosca::image_normalise
(mosca::image& slit_image,
 mosca::image& slit_image_weight,
 int spa_smooth_radius, int disp_smooth_radius, 
 int spa_fit_polyorder, int disp_fit_nknots, double fit_threshold,
 std::vector<T>& slit_spa_norm_profile, std::vector<T>& slit_disp_norm_profile)
{
    //TODO: Check dispersion axis are the same
    if(slit_image.size_x() != slit_image_weight.size_x() ||
       slit_image.size_y() != slit_image_weight.size_y())
        throw std::invalid_argument("image and weight sizes do not match");

    //TODO: Use size_disp rather than size_x in the relevant places
    mosca::image slit_weighted = slit_image;
    std::transform (slit_image.get_data<T>(),
                    slit_image.get_data<T>() + slit_image.size_x() * slit_image.size_y(),
                    slit_image_weight.get_data<T>(),
                    slit_weighted.get_data<T>(), std::multiplies<T>());

    //Collapsing the data to get the profiles in each direction
    std::vector<T> slit_spa_profile = 
            slit_weighted.collapse<T>(mosca::DISPERSION_AXIS);
    std::vector<T> slit_disp_profile = 
            slit_weighted.collapse<T>(mosca::SPATIAL_AXIS);

    //Collapsing the weights 
    std::vector<T> weight_spa_profile = 
            slit_image_weight.collapse<T>(mosca::DISPERSION_AXIS);
    std::vector<T> weight_disp_profile = 
            slit_image_weight.collapse<T>(mosca::SPATIAL_AXIS);

    //Getting the profiles properly weighted
    std::vector<T> slit_spa_profile_w;
    std::transform (slit_spa_profile.begin(), slit_spa_profile.end(),
                    weight_spa_profile.begin(), std::back_inserter(slit_spa_profile_w),
                    std::divides<T>());
    std::vector<T> slit_disp_profile_w;
    std::transform (slit_disp_profile.begin(), slit_disp_profile.end(),
                    weight_disp_profile.begin(), std::back_inserter(slit_disp_profile_w),
                    std::divides<T>());

    T * p_ima = slit_weighted.get_data<T>();
    T total_flux = 
         std::accumulate(p_ima, p_ima + slit_image.size_x() * slit_image.size_y(), T(0));
    T * p_weight = slit_image_weight.get_data<T>();
    T total_weight = 
         std::accumulate(p_weight, p_weight + slit_image.size_x() * slit_image.size_y(), T(0));

    
    if(total_flux == T(0) || total_weight == T(0))
	{
        slit_spa_norm_profile = slit_spa_profile;
        slit_disp_norm_profile = slit_disp_profile;
    	return slit_image;
	}
    //If we are doing any fitting/smoothing in that direction, 
    //initialise it to the current profile, if not initialise it to a constant
    if (spa_smooth_radius > 0 || spa_fit_polyorder > 0)
        slit_spa_norm_profile = slit_spa_profile_w;
    else 
        slit_spa_norm_profile = std::vector<T>(slit_spa_profile_w.size(), 
                T(total_flux / total_weight));
    
    if (disp_smooth_radius > 0 || disp_fit_nknots > 0)
        slit_disp_norm_profile = slit_disp_profile_w;
    else 
        slit_disp_norm_profile = std::vector<T>(slit_disp_profile_w.size(), 
                T(total_flux / total_weight));

    if (spa_smooth_radius > 0)
    {
        std::vector<bool> mask;
        std::transform(weight_spa_profile.begin(), weight_spa_profile.end(), 
                       std::back_inserter(mask), std::bind1st(std::not_equal_to<T>(), T(0)));
        mosca::vector_smooth<T>(slit_spa_norm_profile, mask, spa_smooth_radius);
    }

    if (spa_fit_polyorder > 0)
    {
        std::vector<bool> mask;
        const double max_el = *std::max_element(slit_spa_norm_profile.begin(), slit_spa_norm_profile.end());
        const double th_this = fit_threshold * max_el;
        std::transform(slit_spa_norm_profile.begin(), slit_spa_norm_profile.end(),
                       std::back_inserter(mask), std::bind2nd(std::greater_equal<T>(), th_this));

        size_t used_spa_fit_polyorder = spa_fit_polyorder;  
        mosca::vector_polynomial polfit;
        polfit.fit<T>(slit_spa_norm_profile, mask, used_spa_fit_polyorder);
    }
    
    if (disp_smooth_radius > 0)
    {
        std::vector<bool> mask;
        std::transform(weight_disp_profile.begin(), weight_disp_profile.end(), 
                       std::back_inserter(mask), std::bind1st(std::not_equal_to<T>(), T(0)));
        mosca::vector_smooth<T>(slit_disp_norm_profile, mask, disp_smooth_radius);
    }

    if (disp_fit_nknots > 0)
    {
        std::vector<bool> mask;
        std::transform(weight_disp_profile.begin(), weight_disp_profile.end(), 
                       std::back_inserter(mask), std::bind1st(std::not_equal_to<T>(), T(0)));
        size_t used_disp_fit_nknots = disp_fit_nknots;  
        mosca::vector_cubicspline splfit;
        splfit.fit<T>(slit_disp_norm_profile, mask, 
                      used_disp_fit_nknots);
    }
    
    cpl_size nx = slit_image.size_x();
    cpl_size ny = slit_image.size_y();
    mosca::image result(nx, ny, mosca::type_trait<T>::cpl_eq_type,
                        slit_image.dispersion_axis());
    T * p_res = result.get_data<T>();
    p_weight = slit_image_weight.get_data<T>();
    for (cpl_size j = 0; j< ny; ++j)
    {
        for (cpl_size i = 0; i< nx; ++i, ++p_res, ++p_weight)
        {
            if(*p_weight != 0)
            {
                if(slit_image.dispersion_axis() == mosca::X_AXIS)
                    *p_res = slit_spa_norm_profile[j] * slit_disp_norm_profile[i] / 
                    total_flux * total_weight;
                else
                    *p_res = slit_spa_norm_profile[i] * slit_disp_norm_profile[j] / 
                    total_flux * total_weight;
            }
            else
                *p_res = 1;
                
        }
    }

    return result;
}
