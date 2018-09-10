/* $Id: mosca_image.tcc,v 1.6 2013-07-24 07:40:54 cgarcia Exp $
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
 * $Date: 2013-07-24 07:40:54 $
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifndef MOSCA_IMAGE_TCC
#define MOSCA_IMAGE_TCC

#include <stdexcept>
#include "mosca_image.h"
#include "type_traits.h"

cpl_size mosca::image::size_x() const
{
    return cpl_image_get_size_x(m_image);
}

cpl_size mosca::image::size_y() const
{
    return cpl_image_get_size_y(m_image);
}

cpl_size mosca::image::npix() const
{
    return size_x() * size_y();
}

cpl_size mosca::image::size_spatial() const
{
    if(m_spatial_axis == mosca::X_AXIS)
        return cpl_image_get_size_x(m_image);
    else
        return cpl_image_get_size_y(m_image);
}

cpl_size mosca::image::size_dispersion() const
{
    if(m_dispersion_axis == mosca::Y_AXIS)
        return cpl_image_get_size_y(m_image);
    else
        return cpl_image_get_size_x(m_image);
}

mosca::axis mosca::image::dispersion_axis() const
{
    return m_dispersion_axis;
}

mosca::axis mosca::image::spatial_axis() const
{
    return m_spatial_axis;
}

template<typename T>
T* mosca::image::get_data()
{
    if(mosca::type_trait<T>::cpl_eq_type != cpl_image_get_type(m_image))
        throw std::invalid_argument("type requested does not match image data type");
    T* data = NULL;
    if(m_image != NULL)
        data = static_cast<T*>(cpl_image_get_data(m_image));
    return data;
}

template<typename T>
const T* mosca::image::get_data() const
{
    if(mosca::type_trait<T>::cpl_eq_type != cpl_image_get_type(m_image))
        throw std::invalid_argument("type requested does not match image data type");
    const T* data = NULL;
    if(m_image != NULL)
        data = static_cast<const T*>(cpl_image_get_data_const(m_image));
    return data;
}

template<typename T>
T* mosca::image::get_data_err()
{
    if(mosca::type_trait<T>::cpl_eq_type != cpl_image_get_type(m_image_err))
        throw std::invalid_argument("type requested does not match image data type");
    T* err_data= NULL;
    if(m_image_err != NULL)
        err_data = static_cast<T*>(cpl_image_get_data(m_image_err));
    return err_data;
}

template<typename T>
const T* mosca::image::get_data_err() const
{
    if(mosca::type_trait<T>::cpl_eq_type != cpl_image_get_type(m_image_err))
        throw std::invalid_argument("type requested does not match image data type");
    const T* err_data = NULL;
    if(m_image_err != NULL)
        err_data = static_cast<const T*>(cpl_image_get_data_const(m_image_err));
    return err_data;
}

template<typename T>
std::vector<T> mosca::image::collapse(mosca::axis collapse_axis) const
{
    mosca::axis image_axis;
    int collapse_dir;

    /* If the collapse direction is given in terms of dispersion or spatial
       axis, then decide whether it is X or Y. If not, it is already X or Y */
    image_axis = axis_to_image(collapse_axis);

    if(image_axis == mosca::X_AXIS)
        collapse_dir = 1;
    else
        collapse_dir = 0;

    //TODO: Do this with STL directly, since this is just a sum right now.
    //Addtionally, no copy will be needed.
    image collapse_median (cpl_image_collapse_create(m_image,
                                collapse_dir), true);

    cpl_size size = collapse_median.size_x() * collapse_median.size_y();
    std::vector<T> collapsed_vector(size);

    T * p_image =  collapse_median.get_data<T>();
    for(size_t i = 0; i < size; i++)
        collapsed_vector[i] = p_image[i];

    return collapsed_vector;
}


#endif
