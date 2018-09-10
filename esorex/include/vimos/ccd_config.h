/*
 * This file is part of the MOSCA library
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
 * ccd_config.h
 *
 *  Created on: 2013 11 25
 *      Author: cgarcia
 */

#ifndef CCD_CONFIG_H
#define CCD_CONFIG_H

#include <vector>
#include "rect_region.h"

namespace mosca 
{

class ccd_config
{
public:
    
    struct port_config
    {
        double nominal_gain;
        
        double nominal_ron;

        double computed_gain;
        
        double computed_ron;

        rect_region overscan_region;
        
        rect_region prescan_region;
        
        rect_region validpix_region;
        
    };
    
    ccd_config(const std::vector<port_config>& port_configs,
               double pixel_size, 
               size_t binning_factor_x, 
               size_t binning_factor_y);
    
    ccd_config();
    
    //Channel starts from 0
    double nominal_gain(size_t port) const;
    
    double nominal_ron(size_t port) const;
    
    //Channel starts from 0
    double computed_gain(size_t port) const;
    
    double computed_ron(size_t port) const;
    
    void set_computed_gain(size_t port, double gain);
    
    void set_computed_ron(size_t port, double ron);
    
    size_t nports() const;
    
    size_t valid_npix_x(size_t port) const;

    size_t valid_npix_y(size_t port) const;
    
    const rect_region& overscan_region(size_t port) const;

    const rect_region& prescan_region(size_t port) const;

    const rect_region& validpix_region(size_t port) const;
    
    rect_region port_region(size_t port) const;
    
    rect_region whole_image_region() const;

    rect_region whole_valid_region() const;
    
    int whole_image_npix_x() const;
        
    int whole_image_npix_y() const;
        
    size_t binning_factor_x() const;

    size_t binning_factor_y() const;

    double pixel_size() const;
    
    virtual ~ccd_config();
    
protected:
    
    void check_port(size_t port) const;
    
    std::vector<port_config> m_port_configs;

    double m_pixel_size;
    
    size_t m_binning_factor_x;

    size_t m_binning_factor_y;
    
    friend bool operator ==
      (const ccd_config& lhs, const ccd_config& rhs);
    
    friend bool operator ==
      (const port_config& lhs, const port_config& rhs);
    
};


bool operator == (const mosca::ccd_config& lhs, const mosca::ccd_config& rhs);

bool operator != (const mosca::ccd_config& lhs, const mosca::ccd_config& rhs);



} /* namespace mosca */
#endif /* CCD_CONFIG_H */
