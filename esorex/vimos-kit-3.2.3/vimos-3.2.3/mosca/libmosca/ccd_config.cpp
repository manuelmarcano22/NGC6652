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
 * detector_setup.cpp
 *
 *  Created on: 2013 11 25
 *      Author: cgarcia
 */

#include <stdexcept>
#include <algorithm>
#include <vector>
#include "ccd_config.h"

namespace mosca {

ccd_config::ccd_config(const std::vector<port_config>& port_configs,
                       double pix_size,
                       size_t bin_factor_x, 
                       size_t bin_factor_y) :
        m_port_configs(port_configs), m_pixel_size(pix_size),
        m_binning_factor_x(bin_factor_x), m_binning_factor_y(bin_factor_y)
{
}

ccd_config::ccd_config() 
{
}

ccd_config::~ccd_config()
{
}

double ccd_config::nominal_gain(size_t port) const
{
    check_port(port);        
    return m_port_configs[port].nominal_gain;
}

double ccd_config::nominal_ron(size_t port) const
{
    check_port(port);
    return m_port_configs[port].nominal_ron;
}

double ccd_config::computed_gain(size_t port) const
{
    check_port(port);        
    return m_port_configs[port].computed_gain;
}

double ccd_config::computed_ron(size_t port) const
{
    check_port(port);
    return m_port_configs[port].computed_ron;
}

void ccd_config::set_computed_gain(size_t port, double gain)
{
    check_port(port);
    m_port_configs[port].computed_gain = gain;
}

void ccd_config::set_computed_ron(size_t port, double ron)
{
    check_port(port);
    m_port_configs[port].computed_ron = ron;
}

size_t ccd_config::nports() const
{
    return m_port_configs.size();
}

size_t ccd_config::valid_npix_x(size_t port) const
{
    check_port(port);
    return m_port_configs[port].validpix_region.pix_inc_x(); 
}

size_t ccd_config::valid_npix_y(size_t port)  const
{
    check_port(port);
    return m_port_configs[port].validpix_region.pix_inc_y(); 
}

const rect_region& ccd_config::overscan_region(size_t port) const
{
    check_port(port);
    return m_port_configs[port].overscan_region;
}

const rect_region& ccd_config::prescan_region(size_t port) const
{
    check_port(port);
    return m_port_configs[port].prescan_region;
}

const rect_region& ccd_config::validpix_region(size_t port) const
{
    check_port(port);
    return m_port_configs[port].validpix_region;
}

rect_region ccd_config::port_region(size_t port) const
{
    return rect_region_minenclose(m_port_configs[port].prescan_region,
            m_port_configs[port].overscan_region,
            m_port_configs[port].validpix_region);
}


size_t ccd_config::binning_factor_x() const
{
    return m_binning_factor_x;
}

size_t ccd_config::binning_factor_y() const
{
    return m_binning_factor_y;
}

double ccd_config::pixel_size() const
{
    return m_pixel_size;
}

void ccd_config::check_port(size_t port) const
{
    if(port > nports() -1)
        std::invalid_argument("port does not exist");
}

rect_region ccd_config::whole_image_region() const
{
    std::vector<rect_region> regions;
    for(size_t iport = 0; iport < nports(); iport++)
    {
        if(!m_port_configs[iport].prescan_region.is_empty())
            regions.push_back(m_port_configs[iport].prescan_region);
        if(!m_port_configs[iport].overscan_region.is_empty())
            regions.push_back(m_port_configs[iport].overscan_region);
        regions.push_back(m_port_configs[iport].validpix_region);
    }

    return rect_region_minenclose(regions);
}

rect_region ccd_config::whole_valid_region() const
{
    std::vector<rect_region> regions;
    for(size_t iport = 0; iport < nports(); iport++)
        regions.push_back(m_port_configs[iport].validpix_region);

    return rect_region_minenclose(regions);
}

int ccd_config::whole_image_npix_x() const
{
    return whole_image_region().pix_inc_x();
}

int ccd_config::whole_image_npix_y() const
{
    return whole_image_region().pix_inc_y();
}
bool operator == (const mosca::ccd_config::port_config& lhs,
                  const mosca::ccd_config::port_config& rhs) 
{
    return (lhs.nominal_gain    == rhs.nominal_gain    &&
            lhs.nominal_ron     == rhs.nominal_ron     &&
            lhs.overscan_region == rhs.overscan_region &&
            lhs.prescan_region  == rhs.prescan_region  &&
            lhs.validpix_region == rhs.validpix_region);
}

bool operator == (const mosca::ccd_config& lhs, const mosca::ccd_config& rhs) 
{
    return (lhs.m_port_configs     == rhs.m_port_configs      &&
            lhs.m_pixel_size       == rhs.m_pixel_size        &&
            lhs.m_binning_factor_x == rhs.m_binning_factor_x  &&
            lhs.m_binning_factor_y == rhs.m_binning_factor_y);
}

bool operator != (const mosca::ccd_config& lhs, const mosca::ccd_config& rhs) 
{
  return ( !(lhs == rhs) );
}


} /* namespace mosca */

