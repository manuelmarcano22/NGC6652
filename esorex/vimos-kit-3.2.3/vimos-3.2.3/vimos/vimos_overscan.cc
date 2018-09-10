/* $Id: moses.c,v 1.116 2013/10/15 09:27:38 cgarcia Exp $
 *
 * This file is part of the MOSES library
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013/10/15 09:27:38 $
 * $Revision: 1.116 $
 * $Name:  $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdexcept>
#include <string>
#include <cpl.h>
#include <hdrl.h>
#include "vimos_overscan.h"
#include "ccd_config.h"
#include "statistics.h"


/**
 * @brief
 *   Subtract the overscan from a CCD exposure
 *
 * @param image       Image containing the data to correct
 * @param header      Header of the image
 *
 * @return A newly allocated overscan subtracted image
 *
 */
mosca::image  vimos_preoverscan::subtract_prescan(mosca::image& image, 
                                   const mosca::ccd_config& ccd_config)
{
    int box_hsize = HDRL_OVERSCAN_FULL_BOX;

    //Get number of ports
    size_t nports = ccd_config.nports();
    
    //Create a copy where to store the results
    cpl_image * image_err = cpl_image_duplicate(image.get_cpl_image_err());
    cpl_mask * old_bpm = cpl_image_set_bpm(image_err, 
                            cpl_mask_duplicate(cpl_image_get_bpm(image.get_cpl_image())));
    cpl_mask_delete(old_bpm);
    hdrl_image * target_image = hdrl_image_create(image.get_cpl_image(), image_err);
    cpl_image_delete(image_err);

    //Loop on the ports
    for(size_t iport = 0; iport < nports; iport++)
    {
        hdrl_parameter * collapse_method =  
                hdrl_collapse_median_parameter_create();
        mosca::rect_region ps_reg = ccd_config.prescan_region(iport).coord_0to1();
        hdrl_parameter * prescan_region = ps_reg.hdrl_param(); //deleted by ps_reg

        hdrl_direction correction_direction;
        if(ccd_config.prescan_region(iport).length_x() > 
           ccd_config.prescan_region(iport).length_y())
            correction_direction = HDRL_Y_AXIS;
        else 
            correction_direction = HDRL_X_AXIS;
        
        double ron = ccd_config.computed_ron(iport);
        hdrl_parameter * prescan_params =
                hdrl_overscan_parameter_create(correction_direction,
                                               ron,
                                               box_hsize,
                                               collapse_method,
                                               prescan_region);
                
        hdrl_overscan_compute_result * os_computation =
                hdrl_overscan_compute(image.get_cpl_image(), prescan_params);

        hdrl_image * os_image = 
                hdrl_overscan_compute_result_get_correction(os_computation);
        m_median_correction = 
                cpl_image_get_median(hdrl_image_get_image(os_image));
        
        mosca::rect_region por_reg = ccd_config.validpix_region(iport).coord_0to1(); 
        hdrl_parameter * port_region = por_reg.hdrl_param(); //deleted by por_reg

        hdrl_overscan_correct_result * os_correction =
            hdrl_overscan_correct(target_image, port_region, os_computation);

        hdrl_image * os_corrected_ima = 
            hdrl_overscan_correct_result_get_corrected(os_correction);
        
        hdrl_image * port_image = hdrl_image_extract(os_corrected_ima,
                       ccd_config.validpix_region(iport).coord_0to1().llx(),
                       ccd_config.validpix_region(iport).coord_0to1().lly(),
                       ccd_config.validpix_region(iport).coord_0to1().urx(),
                       ccd_config.validpix_region(iport).coord_0to1().ury());

        hdrl_image_copy(target_image, port_image, 
                        ccd_config.validpix_region(iport).coord_0to1().llx(), 
                        ccd_config.validpix_region(iport).coord_0to1().lly());
        
        hdrl_overscan_compute_result_delete(os_computation);
        hdrl_overscan_correct_result_delete(os_correction);
        hdrl_image_delete(port_image);
        hdrl_parameter_delete(prescan_params);
        
    }
    
    //The rest of the VIMOS pipeline works with float images, so we cast it
    mosca::image os_subtracted(
       cpl_image_cast(hdrl_image_get_image_const(target_image), CPL_TYPE_FLOAT), 
       cpl_image_cast(hdrl_image_get_error_const(target_image), CPL_TYPE_FLOAT),
       true, mosca::Y_AXIS);
    
    hdrl_image_delete(target_image);
    return os_subtracted;
}

//TODO: This is by far not efficient. Too many copies are happening
std::vector<mosca::image> vimos_preoverscan::subtract_prescan
(std::vector<mosca::image>& ima_list, const mosca::ccd_config& ccd_config)
{
    int nima = ima_list.size();
    std::vector<mosca::image> ima_list_os_sub;

    std::vector<double> mean_corr;
    for(int ima = 0; ima < nima; ++ima)
    {
        ima_list_os_sub.push_back(subtract_prescan(ima_list[ima], ccd_config));
        mean_corr.push_back(get_median_correction());
    }
    m_median_correction = mosca::mean(mean_corr.begin(), mean_corr.end());
    return ima_list_os_sub;
}

/**
 * @brief
 *   Subtract the overscan from a CCD exposure
 *
 * @param image       Image containing the data to correct
 * @param header      Header of the image
 *
 * @return A newly allocated overscan subtracted image
 * TODO: Refactor this with subtract_prescan
 */
mosca::image  vimos_preoverscan::subtract_overscan(mosca::image& image, 
                                   const mosca::ccd_config& ccd_config)
{
    int box_hsize = HDRL_OVERSCAN_FULL_BOX;

    //Get number of ports
    size_t nports = ccd_config.nports();
    
    //Create a copy where to store the results
    cpl_image * image_err = cpl_image_duplicate(image.get_cpl_image_err());
    cpl_mask * old_bpm = cpl_image_set_bpm(image_err, 
                            cpl_mask_duplicate(cpl_image_get_bpm(image.get_cpl_image())));
    cpl_mask_delete(old_bpm);
    hdrl_image * target_image = hdrl_image_create(image.get_cpl_image(), image_err);
    cpl_image_delete(image_err);

    //Loop on the ports
    for(size_t iport = 0; iport < nports; iport++)
    {
        hdrl_parameter * collapse_method =  
                hdrl_collapse_median_parameter_create();
        mosca::rect_region ps_reg = ccd_config.overscan_region(iport).coord_0to1();
        hdrl_parameter * overscan_region = ps_reg.hdrl_param(); //deleted by ps_reg

        hdrl_direction correction_direction;
        if(ccd_config.overscan_region(iport).length_x() > 
           ccd_config.overscan_region(iport).length_y())
            correction_direction = HDRL_Y_AXIS;
        else 
            correction_direction = HDRL_X_AXIS;
        
        double ron = ccd_config.computed_ron(iport);
        hdrl_parameter * overscan_params =
                hdrl_overscan_parameter_create(correction_direction,
                                               ron,
                                               box_hsize,
                                               collapse_method,
                                               overscan_region);
                
        hdrl_overscan_compute_result * os_computation =
                hdrl_overscan_compute(image.get_cpl_image(), overscan_params);

        hdrl_image * os_image = 
                hdrl_overscan_compute_result_get_correction(os_computation);
        m_median_correction = 
                cpl_image_get_median(hdrl_image_get_image(os_image));
        
        mosca::rect_region por_reg = ccd_config.validpix_region(iport).coord_0to1(); 
        hdrl_parameter * port_region = por_reg.hdrl_param(); //deleted by por_reg

        hdrl_overscan_correct_result * os_correction =
            hdrl_overscan_correct(target_image, port_region, os_computation);

        hdrl_image * os_corrected_ima = 
            hdrl_overscan_correct_result_get_corrected(os_correction);
        
        hdrl_image * port_image = hdrl_image_extract(os_corrected_ima,
                       ccd_config.validpix_region(iport).coord_0to1().llx(),
                       ccd_config.validpix_region(iport).coord_0to1().lly(),
                       ccd_config.validpix_region(iport).coord_0to1().urx(),
                       ccd_config.validpix_region(iport).coord_0to1().ury());

        hdrl_image_copy(target_image, port_image, 
                        ccd_config.validpix_region(iport).coord_0to1().llx(), 
                        ccd_config.validpix_region(iport).coord_0to1().lly());
        
        hdrl_overscan_compute_result_delete(os_computation);
        hdrl_overscan_correct_result_delete(os_correction);
        hdrl_image_delete(port_image);
        hdrl_parameter_delete(overscan_params);
        
    }
    
    //The rest of the VIMOS pipeline works with float images, so we cast it
    mosca::image os_subtracted(
       cpl_image_cast(hdrl_image_get_image_const(target_image), CPL_TYPE_FLOAT), 
       cpl_image_cast(hdrl_image_get_error_const(target_image), CPL_TYPE_FLOAT),
       true, mosca::Y_AXIS);
    
    hdrl_image_delete(target_image);
    return os_subtracted;
}

//TODO: This is by far not efficient. Too many copies are happening
std::vector<mosca::image> vimos_preoverscan::subtract_overscan
(std::vector<mosca::image>& ima_list, const mosca::ccd_config& ccd_config)
{
    int nima = ima_list.size();
    std::vector<mosca::image> ima_list_os_sub;

    std::vector<double> mean_corr;
    for(int ima = 0; ima < nima; ++ima)
    {
        ima_list_os_sub.push_back(subtract_overscan(ima_list[ima], ccd_config));
        mean_corr.push_back(get_median_correction());
    }
    m_median_correction = mosca::mean(mean_corr.begin(), mean_corr.end());
    return ima_list_os_sub;
}

mosca::image vimos_preoverscan::trimm_preoverscan(mosca::image& image, 
                                           const mosca::ccd_config& ccd_config)
{
    mosca::rect_region crop_region = ccd_config.whole_valid_region();
    mosca::rect_region crop_region_1 = crop_region.coord_0to1(); 
    if(crop_region_1.is_empty())
        throw std::invalid_argument("Region to crop is empty");
    mosca::image image_trimmed = 
            image.trim(crop_region_1.lly(), crop_region_1.llx(),
                       crop_region_1.ury(), crop_region_1.urx());
    return image_trimmed;
}

//TODO: This is by far not efficient. Too many copies are happening
std::vector<mosca::image> vimos_preoverscan::trimm_preoverscan
(std::vector<mosca::image>& ima_list, const mosca::ccd_config& ccd_config)
{
    int nima = ima_list.size();
    std::vector<mosca::image> ima_list_trimmed;

    for(int ima = 0; ima < nima; ++ima)
        ima_list_trimmed.push_back(trimm_preoverscan(ima_list[ima], ccd_config));
    return ima_list_trimmed;
}

void vimos_preoverscan::fix_wcs_trimm(cpl_propertylist * header)
{
    mosca::fiera_config ccd_config(header);
    mosca::rect_region crop_region = ccd_config.whole_valid_region();
    mosca::rect_region crop_region_1 = crop_region.coord_0to1(); 
    if(crop_region_1.is_empty())
        throw std::invalid_argument("Cannot fix WCS from overscan trimming");
    cpl_propertylist_update_double(header, "CRPIX1",
        cpl_propertylist_get_double(header, "CRPIX1") - crop_region_1.llx() + 1);
    cpl_propertylist_update_double(header, "CRPIX2",
        cpl_propertylist_get_double(header, "CRPIX2") - crop_region_1.lly() + 1);
}

double vimos_preoverscan::get_median_correction() const
{
    return m_median_correction;
}

cpl_image * vimos_subtract_prescan(cpl_image * image, cpl_image * image_var, 
                                   cpl_propertylist * header)
{
    mosca::fiera_config config_ccd(header);
    
    cpl_image * image_err = cpl_image_power_create(image_var, 0.5); 
    mosca::image target(image, image_err, false, mosca::Y_AXIS);
    
    vimos_preoverscan scan_corr;
    
    mosca::image target_corr = 
            scan_corr.subtract_prescan(target, config_ccd);
    
    cpl_image_delete(image_err);
    return cpl_image_duplicate(target_corr.get_cpl_image()); 
}

cpl_image * vimos_subtract_overscan(cpl_image * image, cpl_image * image_var, 
                                   cpl_propertylist * header)
{
    mosca::fiera_config config_ccd(header);
    
    cpl_image * image_err = cpl_image_power_create(image_var, 0.5); 
    mosca::image target(image, image_err, false, mosca::Y_AXIS);
    
    vimos_preoverscan scan_corr;
    
    mosca::image target_corr = 
            scan_corr.subtract_overscan(target, config_ccd);
    
    cpl_image_delete(image_err);
    return cpl_image_duplicate(target_corr.get_cpl_image()); 
}

cpl_image * vimos_trimm_preoverscan(cpl_image * image, 
                                    cpl_propertylist * header)
{
    mosca::fiera_config config_ccd(header);
    
    mosca::image target(image, false, mosca::Y_AXIS);
    
    vimos_preoverscan scan_corr;
    
    mosca::image target_corr = 
            scan_corr.trimm_preoverscan(target, config_ccd);
    
    return cpl_image_duplicate(target_corr.get_cpl_image()); 
}



//void fors_trimm_fill_info(cpl_propertylist * header,
//                          const mosca::ccd_config& ccd_config)
//{
//    mosca::rect_region crop_region = ccd_config.whole_valid_region();
//    mosca::rect_region crop_region_1 = crop_region.coord_0to1(); 
//
//    cpl_propertylist_append_int(header, "ESO QC TRIMM LLX",
//                                crop_region_1.llx());
//    cpl_propertylist_append_int(header, "ESO QC TRIMM LLY",
//                                crop_region_1.lly());
//    cpl_propertylist_append_int(header, "ESO QC TRIMM URX",
//                                crop_region_1.urx());
//    cpl_propertylist_append_int(header, "ESO QC TRIMM URY",
//                                crop_region_1.ury());
//}
