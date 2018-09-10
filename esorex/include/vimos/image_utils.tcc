/* $Id: image_utils.tcc,v 1.8 2013-10-16 11:30:00 cgarcia Exp $
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
 * $Date: 2013-10-16 11:30:00 $
 * $Revision: 1.8 $
 * $Name: not supported by cvs2svn $
 */


#ifndef IMAGE_UTILS_TCC
#define IMAGE_UTILS_TCC

#include <iterator>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "mosca_image.h"
#include "image_utils.h"
#include "reduce_method.h"
#include "hdrl.h"


/**
 * @brief
 *   Collapse a list of images in a container with a reduction method
 *   (mean, median, etc...) 
 *
 * @return The mean image
 *
 */
template<typename Iter, typename ReduceMethod>
mosca::image mosca::imagelist_reduce
(Iter image_start, Iter image_end, ReduceMethod reduce_method)
{

    //cpl_errorstate      error_prevstate = cpl_errorstate_get();

    /* Maybe all this can be done with iterators in hdrl */
    hdrl_imagelist * im_list = hdrl_imagelist_new();
    mosca::axis disp_axis = image_start->dispersion_axis();
    cpl_size idx;
    Iter it;

    for(idx = 0, it = image_start; it != image_end; ++it, ++idx)
    {
        mosca::image& im = *it;
        if(im.dispersion_axis() != disp_axis)
            throw std::invalid_argument("Dispersion axes are not the same");
        hdrl_image * tmp = hdrl_image_create(im.get_cpl_image(),
                                             im.get_cpl_image_err());
        
        hdrl_imagelist_set(im_list, tmp, idx);
    }

    //TODO: This depend on the dispersion axis
    hdrl_image * mean;
    cpl_image * contrib;
    hdrl_parameter * reduce_method_hdrl = reduce_method.hdrl_reduce();
    if(hdrl_imagelist_collapse(im_list, reduce_method_hdrl, 
                              &mean, &contrib) != CPL_ERROR_NONE)
    {
        cpl_msg_error(cpl_func,"%s", cpl_error_get_message());
        cpl_msg_error(__func__,"Could not collapse the images");
    }

    hdrl_imagelist_delete(im_list);
    hdrl_parameter_delete(reduce_method_hdrl);
    
    cpl_image * mean_img =  cpl_image_duplicate(hdrl_image_get_image(mean));
    cpl_image * mean_err =  cpl_image_duplicate(hdrl_image_get_error(mean));
    
    hdrl_image_delete(mean);
    cpl_image_delete(contrib);
    
    return mosca::image(mean_img, mean_err, true, disp_axis);
}

template<typename Container>
Container operator/ (Container& image_list, 
                     mosca::image& dividend)
{
    cpl_errorstate      error_prevstate = cpl_errorstate_get();

    Container out_list;
    mosca::axis disp_axis = image_list.begin()->dispersion_axis();
    
    cpl_image * dividend_im = dividend.get_cpl_image();
    cpl_image * dividend_err = dividend.get_cpl_image_err();
    
    cpl_image * dividend_sq = cpl_image_power_create(dividend_im, 2);
    cpl_image * dividend_sq_sq = cpl_image_power_create(dividend_sq, 2);
    cpl_image * dividend_err_sq = cpl_image_power_create(dividend_err, 2);
    
    for(typename Container::iterator it = image_list.begin(); 
            it != image_list.end(); ++it)
    {
        mosca::image& target = *it;
        if(target.dispersion_axis() != disp_axis)
            throw std::invalid_argument("Dispersion axes are not the same");

        cpl_image * target_im = target.get_cpl_image();
        cpl_image * target_err = target.get_cpl_image_err();
        
        cpl_image * target_sq = cpl_image_power_create(target_im, 2);
        cpl_image * target_err_sq = cpl_image_power_create(target_err, 2);

        cpl_image * out = cpl_image_divide_create(target_im, dividend_im);
        
        /* Error propagation: z = y/x, 
         * err_z = sqrt(err_x^2/y^2 + x^2 * err_y^2/y^4)
         */
        cpl_image * out_err = cpl_image_divide_create(target_err_sq, dividend_sq);
        cpl_image * tmp2 = cpl_image_divide_create(target_sq, dividend_sq_sq);
        cpl_image * tmp3 = cpl_image_multiply_create(tmp2, dividend_err_sq);
        cpl_image_add(out_err, tmp3);
        cpl_image_power(out_err, 0.5);
        
        mosca::image image_divided(out, out_err, true, disp_axis);
        
        (*std::inserter(out_list, out_list.end())) = image_divided;

        /* Book keeping */
        cpl_image_delete(target_sq);
        cpl_image_delete(target_err_sq);
        cpl_image_delete(tmp2);
        cpl_image_delete(tmp3);
    }
    cpl_image_delete(dividend_sq);
    cpl_image_delete(dividend_sq_sq);
    cpl_image_delete(dividend_err_sq);
    
    if(!cpl_errorstate_is_equal(error_prevstate))
    {
        cpl_msg_error(cpl_func,"%s", cpl_error_get_message());
        cpl_msg_error(__func__,"Could not divide images");
    }

    return out_list;
}

#endif
