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
 * fiera_ccd.cpp
 *
 *  Created on: 2013 11 25
 *      Author: cgarcia
 */

#include <cpl.h>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include "fiera_config.h"

namespace mosca 
{

fiera_config::fiera_config(const cpl_propertylist * header)
{
    if(header == NULL)
        std::invalid_argument("empty header");
    
    char * keyname;
    size_t n_ports = 1;

    if (cpl_propertylist_has(header, "ESO DET OUTPUTS"))
        n_ports = cpl_propertylist_get_int(header, "ESO DET OUTPUTS");

    if(n_ports !=1 && n_ports !=4)
        std::invalid_argument("Setup nor supported: Fiera detectors with "
                "only 1 or 4 ports are supported");

    //Binning 
    m_binning_factor_x = cpl_propertylist_get_int(header, "ESO DET WIN1 BINX"); 
    m_binning_factor_y = cpl_propertylist_get_int(header, "ESO DET WIN1 BINY");
    
    //Pixel size
    m_pixel_size = cpl_propertylist_get_double(header, "ESO DET CHIP1 PSZX");

    //Chip name
    m_chip_id = 
         std::string(cpl_propertylist_get_string(header, "ESO DET CHIP1 ID"));

    //Chip dimensions
    //int det_nx = 
    //        cpl_propertylist_get_int(header, "ESO DET CHIP1 NX");
    //int det_ny = 
    //        cpl_propertylist_get_int(header, "ESO DET CHIP1 NY");

    //Image dimensions
    int ima_nx = 
            cpl_propertylist_get_int(header, "NAXIS1");
    int ima_ny = 
            cpl_propertylist_get_int(header, "NAXIS2");

    //Relative positions of the ports readouts with respect to the whole chip
    std::vector<int> port_xpos(n_ports);
    std::vector<int> port_ypos(n_ports);
    for(size_t i = 0; i < n_ports; i++) 
    {
        keyname = cpl_sprintf("ESO DET OUT%zu X", i+1);
        port_xpos[i] = 
            cpl_propertylist_get_int(header, keyname);
        cpl_free(keyname);
        keyname = cpl_sprintf("ESO DET OUT%zu Y", i+1);
        port_ypos[i] = 
            cpl_propertylist_get_int(header, keyname);
        cpl_free(keyname);
    }
    //This is the left lower corner of the final image with respect to the
    //whole chip
    int min_port_xpos = *std::min_element(port_xpos.begin(), port_xpos.end());
    int min_port_ypos = *std::min_element(port_ypos.begin(), port_ypos.end());
    
    //Loop on the ports to get their configuration
    for(size_t i = 0; i < n_ports; i++)
    {
        keyname = cpl_sprintf("ESO DET OUT%zu GAIN", i+1);
        double nom_gain = 
            cpl_propertylist_get_double(header, keyname);
        cpl_free(keyname);
        keyname = cpl_sprintf("ESO DET OUT%zu RON", i+1);
        double nom_ron = 
            cpl_propertylist_get_double(header, keyname);
        cpl_free(keyname);
        keyname = cpl_sprintf("ESO DET OUT%zu NX", i+1);
        
        int port_nx = 
            cpl_propertylist_get_int(header, keyname);
        cpl_free(keyname);
        keyname = cpl_sprintf("ESO DET OUT%zu NY", i+1);
        int port_ny = 
            cpl_propertylist_get_int(header, keyname);
        cpl_free(keyname);
        keyname = cpl_sprintf("ESO DET OUT%zu PRSCX", i+1);
        int pre_x= 
            cpl_propertylist_get_int(header, keyname);
        cpl_free(keyname);
        keyname = cpl_sprintf("ESO DET OUT%zu PRSCY", i+1);
        int pre_y = 0;
        if (cpl_propertylist_has(header, keyname))
            pre_y = cpl_propertylist_get_int(header, keyname);
        cpl_free(keyname);
        keyname = cpl_sprintf("ESO DET OUT%zu OVSCX", i+1);
        int ovs_x= 
            cpl_propertylist_get_int(header, keyname);
        cpl_free(keyname);
        keyname = cpl_sprintf("ESO DET OUT%zu OVSCY", i+1);
        int ovs_y = 0; 
        if (cpl_propertylist_has(header, keyname))
            ovs_y = cpl_propertylist_get_int(header, keyname);
        cpl_free(keyname);
        
        if(pre_x != 0 && pre_y !=0)
            std::invalid_argument("Two prescan regions are not supported");
        if(ovs_x != 0 && ovs_y !=0)
            std::invalid_argument("Two overscan regions are not supported");
        
        //Get the port offset with respect to the final image
        int port_off_x;
        if(port_xpos[i] == 1)
            port_off_x = port_xpos[i] - min_port_xpos;
        else 
            port_off_x = ima_nx - port_nx - pre_x - ovs_x;

        int port_off_y;
        if(port_ypos[i] == 1)
            port_off_y = port_ypos[i] - min_port_ypos;
        else 
            port_off_y = ima_ny - port_ny - pre_y - ovs_y;
        
        //Computing Prescan region
        rect_region prescan_reg;
        if(pre_x != 0) 
        {
            //Prescan is before any other real pixel
            if(port_xpos[i] == 1)
                prescan_reg = rect_region(port_off_x + 0, 
                                          port_off_y + 0, 
                                          port_off_x + pre_x - 1, 
                                          port_off_y + port_ny - 1);
            //Prescan is after other real pixels because the port is read
            //opposite to the normal direction
            else
                prescan_reg = rect_region(port_off_x + port_nx + ovs_x, 
                                          port_off_y + 0, 
                                          port_off_x + port_nx + pre_x + ovs_x - 1, 
                                          port_off_y + port_ny - 1);
        }
        
        if(pre_y != 0)
        {
            if(port_ypos[i] == 1)
                prescan_reg = rect_region(port_off_x + 0, 
                                          port_off_y + 0, 
                                          port_off_x + port_nx - 1, 
                                          port_off_y + pre_y - 1);
            else
                prescan_reg = rect_region(port_off_x + 0,
                                          port_off_y + port_ny + ovs_y, 
                                          port_off_x + port_nx - 1, 
                                          port_off_y + port_ny + pre_y + ovs_y - 1); 
        }

        rect_region overscan_reg;
        if(ovs_x != 0)
        {
            //Overscan is at the end of any real pixel. The image might
            //also contain prescan, so this has to be taken into account
            if(port_xpos[i] == 1)
                overscan_reg = rect_region(port_off_x + port_nx + pre_x, 
                                           port_off_y + 0, 
                                           port_off_x + port_nx + pre_x + ovs_x - 1, 
                                           port_off_y + port_ny - 1);
            else
                overscan_reg = rect_region(port_off_x + 0, 
                                           port_off_y + 0, 
                                           port_off_x + pre_x - 1, 
                                           port_off_y + port_ny - 1);
        }
        
        if(ovs_y != 0)
        {
            if(port_ypos[i] == 1)
                overscan_reg = rect_region(port_off_x + 0,
                                           port_off_y + port_ny + pre_y, 
                                           port_off_x + port_nx - 1,
                                           port_off_y + port_ny + pre_y + ovs_y - 1);
            else 
                overscan_reg = rect_region(port_off_x + 0, 
                                           port_off_y + 0, 
                                           port_off_x + port_nx - 1,
                                           port_off_y + pre_y - 1);
        }

        rect_region validpix_reg;
        if(port_ypos[i] == 1)
            if(port_xpos[i] == 1)
                validpix_reg = rect_region(port_off_x + pre_x, 
                                           port_off_y + pre_y, 
                                           port_off_x + port_nx + pre_x - 1, 
                                           port_off_y + port_ny + pre_y - 1); 
            else
                validpix_reg = rect_region(port_off_x + ovs_x, 
                                           port_off_y + pre_y, 
                                           port_off_x + port_nx + ovs_x - 1, 
                                           port_off_y + port_ny + pre_y - 1); 
        else
            if(port_xpos[i] == 1)
                validpix_reg = rect_region(port_off_x + pre_x, 
                                           port_off_y + ovs_y, 
                                           port_off_x + port_nx + pre_x - 1, 
                                           port_off_y + port_ny + ovs_y - 1); 
            else
                validpix_reg = rect_region(port_off_x + ovs_x, 
                                           port_off_y + ovs_y, 
                                           port_off_x + port_nx + ovs_x - 1, 
                                           port_off_y + port_ny + ovs_y - 1); 
                
        //The convention for the regions is to start from 0 for the first pixel.
        port_config config;
        config.nominal_gain    = nom_gain;
        config.nominal_ron     = nom_ron;
        config.computed_gain   = 0.;
        config.computed_ron    = 0.;
        config.overscan_region = overscan_reg;
        config.prescan_region  = prescan_reg;
        config.validpix_region = validpix_reg;
        
        if(validpix_reg.is_empty())
            std::invalid_argument("Port with empty valid pixel area");
        
        //Add to the port configurations
        m_port_configs.push_back(config); 
    }
    
}

fiera_config::fiera_config()
{
}

fiera_config::~fiera_config()
{
}

} /* namespace mosca */
