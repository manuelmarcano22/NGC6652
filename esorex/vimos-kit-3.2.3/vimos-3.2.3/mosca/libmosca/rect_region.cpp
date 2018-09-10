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
 * rect_region.cpp
 *
 *  Created on: 2013 11 25
 *      Author: cgarcia
 */

#include <stdexcept>
#include <algorithm>
#include "rect_region.h"


namespace mosca {

//Which convention: starting from 0 or from 1?
rect_region::rect_region(int llx_pix, int lly_pix, int urx_pix, int ury_pix) :
        m_llx(llx_pix), m_lly(lly_pix), m_urx(urx_pix), m_ury(ury_pix),
        m_hdrl_param(NULL),
        m_is_empty(false)
{
    if(urx_pix < llx_pix || ury_pix < lly_pix)
        throw std::invalid_argument("Upper right coordinates smaller than "
                "lower left coordinates");
}

rect_region::rect_region() :
  m_llx(0), m_lly(0), m_urx(0), m_ury(0),
  m_hdrl_param(NULL),
  m_is_empty(true)
{
}

rect_region::~rect_region()
{
    if (m_hdrl_param != NULL)
        hdrl_parameter_destroy(m_hdrl_param);
}

int rect_region::llx() const
{
    return m_llx;
}

int rect_region::lly() const
{
    return m_lly;
}

int rect_region::urx() const
{
    return m_urx;
}

int rect_region::ury() const
{
    return m_ury;
}

int rect_region::length_x() const
{
    if(is_empty())
        return 0;
    return m_urx - m_llx;
}

int rect_region::length_y() const
{
    if(is_empty())
        return 0;
    return m_ury - m_lly;
}

int rect_region::pix_inc_x() const
{
    if(is_empty())
        return 0;
    return m_urx - m_llx + 1;
}
int rect_region::pix_inc_y() const
{
    if(is_empty())
        return 0;
    return m_ury - m_lly + 1;
}

int rect_region::area_geom() const
{
    if(is_empty())
        return 0;
    return (length_x() * length_y());
}

int rect_region::area_pix_inc() const
{
    if(is_empty())
        return 0;
    return ((length_x() + 1)* (length_y() +1));
}

void rect_region::set_llx(int llx_pix)
{
    m_llx = llx_pix;
}

void rect_region::set_lly(int lly_pix)
{
    m_lly = lly_pix;
}

void rect_region::set_urx(int urx_pix)
{
    m_urx = urx_pix;
}

void rect_region::set_ury(int ury_pix)
{
    m_ury = ury_pix;
}

bool rect_region::is_empty() const
{
    return m_is_empty;
}


//The returned hdrl_parameter doesn't have to be deallocated
hdrl_parameter * rect_region::hdrl_param()
{
    if(is_empty())
        return NULL;
    else
        if (m_hdrl_param == NULL)
            m_hdrl_param =  hdrl_rect_region_parameter_create
            (m_llx, m_lly, m_urx, m_ury);
    return m_hdrl_param; 
}

//Convert from 0 starting to 1 starting convention
rect_region rect_region::coord_0to1(void) const
{
    if(is_empty())
        return rect_region();
    return rect_region(m_llx + 1, m_lly + 1, m_urx + 1, m_ury + 1);
}

//Convert from 1 starting to 0 starting convention
rect_region rect_region::coord_1to0(void) const
{
    if(is_empty())
        return rect_region();
    return rect_region(m_llx - 1, m_lly - 1, m_urx - 1, m_ury - 1);
}

//Return the minimun rectangular region that encloses the other two
rect_region rect_region_minenclose
(const rect_region & reg1, const rect_region & reg2)
{
    std::vector<rect_region> regions;
    regions.push_back(reg1);
    regions.push_back(reg2);
    return rect_region_minenclose(regions);
}

//Return the minimun rectangular region that encloses the other three
rect_region rect_region_minenclose
(const rect_region & reg1, const rect_region & reg2, const rect_region & reg3)
{
    std::vector<rect_region> regions;
    regions.push_back(reg1);
    regions.push_back(reg2);
    regions.push_back(reg3);
    return rect_region_minenclose(regions);
}

//Return the minimun rectangular region that encloses  all other ones
rect_region rect_region_minenclose(const std::vector<rect_region> & regions)
{
    std::vector<int> llx;
    std::vector<int> lly;
    std::vector<int> urx;
    std::vector<int> ury;
    for(size_t ireg = 0; ireg < regions.size(); ireg++)
    {
        if(regions[ireg].is_empty())
            std::invalid_argument("Input regions cannot be empty");
        llx.push_back(regions[ireg].llx());
        lly.push_back(regions[ireg].lly());
        urx.push_back(regions[ireg].urx());
        ury.push_back(regions[ireg].ury());
    }
    int min_llx = *std::min_element(llx.begin(), llx.end());
    int min_lly = *std::min_element(lly.begin(), lly.end());
    int max_urx = *std::max_element(urx.begin(), urx.end());
    int max_ury = *std::max_element(ury.begin(), ury.end());

    return rect_region(min_llx, min_lly, max_urx, max_ury);
}

bool operator == (const mosca::rect_region& lhs, const mosca::rect_region& rhs) 
{
    return (lhs.m_is_empty == rhs.m_is_empty &&
            lhs.m_llx      == rhs.m_llx &&
            lhs.m_lly      == rhs.m_lly &&
            lhs.m_urx      == rhs.m_urx &&
            lhs.m_ury      == rhs.m_ury);
}

bool operator != (const mosca::rect_region& lhs, const mosca::rect_region& rhs) 
{
  return ( !(lhs == rhs) );
}

} /* namespace mosca */

