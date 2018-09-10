/* $Id: image_smooth.cpp,v 1.7 2013/09/06 13:11:02 cgarcia Exp $
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
 * $Date: 2013/09/06 13:11:02 $
 * $Revision: 1.7 $
 * $Name:  $
 */


#ifndef IMAGE_SMOOTH_CPP
#define IMAGE_SMOOTH_CPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "mosca_image.h"
#include "image_smooth.h"


/**
 * @brief
 *   Smooth an image with a 
 *
 * @param spectrum  A 1D emission line spectrum
 * @param level     Significance level
 * @param exp_width Expected lines FWHM (in pixels)
 *
 * @return List of peaks candidates positions
 * TODO: share most of the code with image_spline_fit and image_pol_fit
 *
 */
template<typename T>
void mosca::image_smooth_1d_median(mosca::image& image,
                                   int half_width, mosca::axis smooth_axis)
{
    if(image.get_cpl_image() == NULL)
        return;
    
    cpl_size ny = image.size_y();
    cpl_size nx = image.size_x();
    cpl_size nsmooth;
    mosca::axis image_axis;
    int collapse_dir;

    /* If the smooth direction is given in terms of dispersion or spatial
       axis, then decide whether it is X or Y. If not, it is already X or Y */
    image_axis = image.axis_to_image(smooth_axis);
    
    /* Lines are contiguos pixels along the smoothing direction.
       i.e., if image_axis == mosca::X_AXIS the lines are rows and there are
       ny lines.
       Each line contains nsmooth pixels. For image_axis == mosca::X_AXIS this
       means nsmooth = nx. */
    if(image_axis == mosca::X_AXIS)
    {
        nsmooth = nx;
        collapse_dir = 0;
    }
    else
    {
        nsmooth = ny;
        collapse_dir = 1;
    }

    if(2*half_width > nsmooth)
        throw std::out_of_range("2*half_width is larger than the image size");
    
    //We collapse in the opposite direction of the fitting
    cpl_image * collapse_median = 
          cpl_image_collapse_median_create(image.get_cpl_image(),
                                           collapse_dir, 0, 0);
    
    //We get the vector of the collapsed image
    cpl_vector * collapse_line;
    if(image_axis == mosca::X_AXIS)
        collapse_line = cpl_vector_new_from_image_row(collapse_median, 1);
    else
        collapse_line = cpl_vector_new_from_image_column(collapse_median, 1);

    //We smooth the line-image
    cpl_vector * smooth_line = cpl_vector_filter_median_create(collapse_line,
                                                               half_width);

    //Place in each column or row the smoothed line
    T *      p_image      = image.get_data<T>();
    double * p_smoothed   = cpl_vector_get_data(smooth_line);
    for (cpl_size j = 0; j< ny; ++j)
    {
        for (cpl_size i = 0; i< nx; ++i, ++p_image)
        {
            if(image_axis == mosca::X_AXIS)
                *p_image = p_smoothed[i];
            else
                *p_image =  p_smoothed[j];
        }
    }

    cpl_image_delete(collapse_median);    
    cpl_vector_delete(collapse_line);
    cpl_vector_delete(smooth_line);
}
        
#endif
