/* $Id: mosca_image.cpp,v 1.17 2013-09-06 13:10:23 cgarcia Exp $
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
 * $Date: 2013-09-06 13:10:23 $
 * $Revision: 1.17 $
 * $Name: not supported by cvs2svn $
 */

#include <stdexcept>
#include <iostream>
#include "mosca_image.h"
#include "cpl.h"

mosca::image::image(cpl_size nx, cpl_size ny, cpl_type type, 
                    mosca::axis disp_axis) :
  m_image(NULL), m_image_err(NULL)
{
    m_image     = cpl_image_new(nx, ny, type);
    m_image_err = cpl_image_new(nx, ny, type);
    m_take_over = true;
    set_axis(disp_axis);
}
    
mosca::image::image()
{
    m_image = NULL;
    m_image_err = NULL;
    m_take_over = false;
    m_spatial_axis = mosca::X_AXIS;
    m_dispersion_axis = mosca::Y_AXIS;
}
    
mosca::image::image(cpl_image * ima, bool take_over, 
                    mosca::axis disp_axis) :
  m_image(NULL), m_image_err(NULL)
{
    m_image = ima;
    m_image_err = NULL;
    m_take_over = take_over;
    set_axis(disp_axis);
}

mosca::image::image(cpl_image * ima, cpl_image * ima_err, 
                    bool take_over, mosca::axis disp_axis) :
  m_image(NULL), m_image_err(NULL)
{
    if(m_image_err != NULL)
        if((cpl_image_get_size_x(ima) != cpl_image_get_size_x(ima_err)) ||
                (cpl_image_get_size_y(ima) != cpl_image_get_size_y(ima_err)))
            throw std::invalid_argument("Data and error should have the same size");
    m_image = ima;
    m_image_err = ima_err;
    m_take_over = take_over;
    set_axis(disp_axis);
}

mosca::image::image(const image& rhs) :
  m_image(NULL), m_image_err(NULL)
{
    if(rhs.m_image != NULL)
        m_image = cpl_image_duplicate(rhs.m_image);
    if(rhs.m_image_err != NULL)
        m_image_err = cpl_image_duplicate(rhs.m_image_err);
    m_take_over = true;
    set_axis(rhs.m_dispersion_axis);    
}

mosca::image::~image()
{
    if(m_take_over) 
    {
        if(m_image != NULL)
            cpl_image_delete(m_image);
        if(m_image_err != NULL)
            cpl_image_delete(m_image_err);
    }
}

void mosca::image::set_axis(mosca::axis disp_axis)
{
    m_dispersion_axis = disp_axis;
    if(disp_axis == mosca::X_AXIS)
        m_spatial_axis    = mosca::Y_AXIS;
    else if(disp_axis == mosca::Y_AXIS)
        m_spatial_axis    = mosca::X_AXIS;
    else
        throw std::invalid_argument("Invalid dispersion axis. "
                                    "Use X_AXIS or Y_AXIS");
}

mosca::image& mosca::image::operator= (const mosca::image& rhs)
{
    if(rhs.m_image != NULL)
        m_image = cpl_image_duplicate(rhs.m_image);
    if(rhs.m_image_err != NULL)
        m_image_err = cpl_image_duplicate(rhs.m_image_err);
    m_take_over = true;
    set_axis(rhs.m_dispersion_axis);    
    return *this;
}

//starting with 0 or 1?
mosca::image mosca::image::trim(cpl_size disp_bottom, cpl_size spa_bottom,
                                cpl_size disp_top,    cpl_size spa_top) const
{
    cpl_image * trimmed_image = NULL;
    cpl_image * trimmed_image_err = NULL;
    if(m_dispersion_axis == mosca::X_AXIS)
    {
        trimmed_image = cpl_image_extract
            (m_image, disp_bottom, spa_bottom, disp_top, spa_top);
        if(cpl_image_count_rejected(m_image) != 0)
        {
            cpl_mask * trimmed_mask = cpl_mask_extract(cpl_image_get_bpm(m_image), 
                    disp_bottom, spa_bottom, disp_top, spa_top);
            cpl_image_set_bpm(trimmed_image, trimmed_mask);
        }
        if(m_image_err != NULL)
            trimmed_image_err = cpl_image_extract
                (m_image_err, disp_bottom, spa_bottom, disp_top, spa_top);
    }
    else
    {
        trimmed_image = cpl_image_extract
            (m_image, spa_bottom, disp_bottom, spa_top, disp_top);
        if(cpl_image_count_rejected(m_image) != 0)
        {
            cpl_mask * trimmed_mask = cpl_mask_extract(cpl_image_get_bpm(m_image), 
                    spa_bottom, disp_bottom, spa_top, disp_top);
            cpl_image_set_bpm(trimmed_image, trimmed_mask);
        }
        if(m_image_err != NULL)
            trimmed_image_err = cpl_image_extract
                (m_image_err, spa_bottom, disp_bottom, spa_top, disp_top);
    }
    mosca::image result(trimmed_image, trimmed_image_err, 
                        true, m_dispersion_axis); 
    return result;
}

mosca::axis mosca::image::axis_to_image(mosca::axis an_axis) const
{
    mosca::axis image_axis;
    /* If an axis is given in terms of dispersion or spatial
       axis, then decide whether it is X or Y. If not, it is already X or Y */ 
    if(an_axis == mosca::DISPERSION_AXIS)
        image_axis = dispersion_axis();
    else if(an_axis == mosca::SPATIAL_AXIS)
        image_axis = spatial_axis();
    else
        image_axis = an_axis;
    
    return image_axis;
}


cpl_image* mosca::image::get_cpl_image()
{
    return m_image;
}

cpl_image* mosca::image::get_cpl_image_err()
{
    return m_image_err;
}

const cpl_image* mosca::image::get_cpl_image() const
{
    return m_image;
}

const cpl_image* mosca::image::get_cpl_image_err() const
{
    return m_image_err;
}
