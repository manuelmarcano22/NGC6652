/* $Id: mosca_image.h,v 1.13 2013-09-06 13:10:41 cgarcia Exp $
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
 * $Date: 2013-09-06 13:10:41 $
 * $Revision: 1.13 $
 * $Name: not supported by cvs2svn $
 */

#ifndef MOSCA_IMAGE_H
#define MOSCA_IMAGE_H


#include <vector>
#include "cpl_image.h"

namespace mosca
{

enum axis
{
    X_AXIS,
    Y_AXIS,
    DISPERSION_AXIS,
    SPATIAL_AXIS
};

/**
 * This class represents a 2d image. It uses cpl_image as an internal
 * image representation
 */
class image
{
public:
    image(cpl_size nx, cpl_size ny, cpl_type type = CPL_TYPE_DOUBLE,
          axis dispersion_axis = X_AXIS);
    
    image(cpl_image * image, bool take_over = false, 
          axis dispersion_axis = X_AXIS);
    
    image(cpl_image * image, cpl_image * image_err, bool take_over = false, 
          axis dispersion_axis = X_AXIS);
    
    image();

    image(const image& other);
    
    image& operator= (const image& other);
    
    virtual ~image();
    
    inline cpl_size size_x() const;

    inline cpl_size size_y() const;

    inline cpl_size npix() const;

    inline cpl_size size_spatial() const;

    inline cpl_size size_dispersion() const;

    inline axis dispersion_axis() const;
    
    inline axis spatial_axis() const;

    //Return an axis in terms of image axis (X or Y)
    axis axis_to_image(axis an_axis) const;
    
    template<typename T>
    std::vector<T> collapse(axis collapse_axis) const;
    
    template<typename T>
    inline T* get_data();

    template<typename T>
    inline const T* get_data() const;

    template<typename T>
    inline T* get_data_err();

    template<typename T>
    inline const T* get_data_err() const;

    image trim(cpl_size disp_bottom, cpl_size spa_bottom,
               cpl_size disp_top,    cpl_size spa_top) const;

    cpl_image * get_cpl_image();

    cpl_image * get_cpl_image_err();

    const cpl_image * get_cpl_image() const;

    const cpl_image * get_cpl_image_err() const;

private:
    
    void set_axis(axis dispersion_axis);
    
    axis m_dispersion_axis;

    axis m_spatial_axis;

    bool m_take_over;
    
    cpl_image *  m_image;
    
    cpl_image *  m_image_err;
};

}

#endif

#include "mosca_image.tcc"
