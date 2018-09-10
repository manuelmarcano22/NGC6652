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
 * vimos_detmodel.cc
 *
 *  Created on: 2014 7 31
 *      Author: cgarcia
 */


#include <stdexcept>
#include <sstream>
#include "fiera_config.h"
#include "vimos_detmodel.h"

//The detector model is the usual for CCDs: a readout noise component
//and a photon noise component scaled by the gain. 
//The readout noise is taken from the overscan region of each readout port
//The gain is taken as the nominal gain from the detector.
//If ccd_config is created with fiera_config(header), the nominal gain
//is taken from the header as ESO DET OUT? GAIN, for instance, also per 
//readout port.
cpl_image * vimos_image_variance_from_detmodel(cpl_image * image,
                                               const mosca::ccd_config& ccd_config)
{

    //Check that the overscan has not been trimmed:
    if(cpl_image_get_size_x(image) != ccd_config.whole_image_npix_x() ||
       cpl_image_get_size_y(image) != ccd_config.whole_image_npix_y())
        throw std::invalid_argument("Pre/Overscan has already been trimmed. "
                                    "Cannot compute variance for detector");
    
    //Images which contain the model parameters per port
    cpl_size im_nx = cpl_image_get_size_x(image);
    cpl_size im_ny = cpl_image_get_size_y(image);
    cpl_image * ron_sq_im  = cpl_image_new(im_nx, im_ny, CPL_TYPE_FLOAT);
    cpl_image * gain_im = cpl_image_new(im_nx, im_ny, CPL_TYPE_FLOAT);
    cpl_image * os_im   = cpl_image_new(im_nx, im_ny, CPL_TYPE_FLOAT);
    
    //Loop on each port:
    for(size_t iport = 0; iport < ccd_config.nports(); iport++)
    {
        mosca::rect_region os_reg = 
                ccd_config.overscan_region(iport).coord_0to1();
        if(os_reg.is_empty())
        {
            cpl_image_delete(ron_sq_im);
            cpl_image_delete(gain_im);
            cpl_image_delete(os_im);
            throw std::invalid_argument("Overscan area is empty. Cannot compute "
                    "detector noise model");
        }
        
        double ron  = ccd_config.computed_ron(iport);
        
        double gain = ccd_config.nominal_gain(iport);
        
        double os_level = cpl_image_get_median_window
          (image, os_reg.llx(), os_reg.lly(), os_reg.urx(), os_reg.ury());

        mosca::rect_region port_reg = 
                        ccd_config.port_region(iport).coord_0to1();

        //TODO: Use cpl_image_fill_window when it exists.
        for(cpl_size ix = port_reg.llx(); ix <= port_reg.urx(); ++ix)
            for(cpl_size iy = port_reg.lly(); iy <= port_reg.ury(); ++iy)
            {
                cpl_image_set(ron_sq_im, ix, iy, ron * ron);
                cpl_image_set(gain_im, ix, iy, gain);
                cpl_image_set(os_im, ix, iy, os_level);
            }
    }
    
    cpl_image * ima_sub_os = cpl_image_subtract_create(image, os_im);
    /* flux_photoel = flux_adu / gain
     * poisson_variance_photel = flux_photoel  -> Poisson law
     * poisson_variance_adu    = poisson_variance_photel * gain * gain
     * poisson_variance_adu    = flux_adu * gain */
    cpl_image * poisson_var_adu = cpl_image_multiply_create(ima_sub_os, gain_im);

    cpl_image * variance = cpl_image_add_create(poisson_var_adu, ron_sq_im);
    
    cpl_image_delete(ron_sq_im);
    cpl_image_delete(gain_im);
    cpl_image_delete(os_im);
    cpl_image_delete(ima_sub_os);
    cpl_image_delete(poisson_var_adu);
    
    return variance;
}

cpl_image * vimos_image_variance_from_detmodel
(cpl_image * image, cpl_propertylist * header, cpl_propertylist * mbias_header)
{
    mosca::fiera_config config_ccd(header);
    
    if(mbias_header == NULL)
        return NULL;

    //Get number of ports
    size_t nports = config_ccd.nports();

    //Loop on the ports
    for(size_t iport = 0; iport < nports; iport++)
    {
        std::ostringstream key_stream;
        key_stream<<"ESO QC DET OUT"<<iport+1<<" RON";
        double computed_ron =
                cpl_propertylist_get_double(mbias_header,
                                            key_stream.str().c_str());
        config_ccd.set_computed_ron(iport, computed_ron);
    }

    return vimos_image_variance_from_detmodel(image, config_ccd);
}
