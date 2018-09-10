/* 
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

#ifndef FLAT_COMBINE_TCC
#define FLAT_COMBINE_TCC

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include "flat_combine.h"
#include "mosca_image.h"
#include "reduce_method.h"
#include "image_utils.h"
#include "vector_utils.h"


/**
 * @brief
 *   Get a master spectroscopy flat 
 *
 * @param image_start An iterator with the first flat image
 * @param image_end   An iterator with the end of the list of flat images.
 *
 * @tparam T    The type of data of the images: float, double, etc...
 * @tparam Iter The type of iterator. If dereferenced it should return a 
 *              mosca::image object
 * TODO: Better interpolation of each value in the SED image
 */
template<typename T, typename Iter, typename CombineMethod>
std::auto_ptr<mosca::image> 
mosca::flat_combine(Iter flat_start, Iter flat_end, 
                    const std::vector<mosca::calibrated_slit>& slits,
                    const mosca::wavelength_calibration& wave_cal,
                    const mosca::grism_config& grism_cfg,
                    size_t smooth_size, CombineMethod comb_method)
{

    /* TODO: Check sizes of the flats are the same */
    
    /* TODO: Shouldn't be hard-coded DOUBLE */
    cpl_image * master_flat_im = cpl_image_new(flat_start->size_dispersion(), flat_start->size_spatial(), CPL_TYPE_DOUBLE); 
    cpl_image * master_flat_err = cpl_image_new(flat_start->size_dispersion(), flat_start->size_spatial(), CPL_TYPE_DOUBLE); 
    
    /* We work on a slit per slit basis */
    for(std::vector<mosca::calibrated_slit>::const_iterator slit_it = slits.begin();
        slit_it != slits.end() ; ++slit_it)
    {
        
        /* Create the slit mask */
        cpl_mask * slit_mask_whole =
                slit_it->get_mask_valid(flat_start->dispersion_axis());

        /* If this slit doesn't contain any belonging pixel */
        if(cpl_mask_count(slit_mask_whole) == 0)
            continue;
        
        /* TODO: This limits are wrong. A larger margin has to be applied. 
         * Together with the mask, there shouldn't be a problem if a security 
         * margin to the disp_top and spa_top is added.  
         */
        int disp_bottom, spa_bottom, disp_top, spa_top;
        slit_it->get_extent_pix(disp_bottom, spa_bottom,
                                disp_top,    spa_top);
        if(spa_bottom < 1)
            spa_bottom = 1; //TODO:: Probably this has to be done at the time of slit creation, to make sure that only slits within the image limits can be created.
        if(spa_top > flat_start->size_spatial())
            spa_top = flat_start->size_spatial(); //TODO: The same as above

        /* Trim the mask */
        cpl_mask * trimmed_mask = cpl_mask_extract(slit_mask_whole, 
                disp_bottom, spa_bottom, disp_top, spa_top);
        cpl_image * slit_mask = cpl_image_new_from_mask(trimmed_mask);
        int * slit_mask_pix = cpl_image_get_data_int(slit_mask);
        
        double start_wave = grism_cfg.start_wave();
        double end_wave = grism_cfg.end_wave();
        double mean_dispersion = wave_cal.mean_dispersion
               (start_wave, end_wave,
                slit_it->get_position_spatial_corrected(), 
                slit_it->get_position_spatial_corrected() + 
                     slit_it->get_length_spatial_corrected());
        size_t n_pix_sed = 
               (size_t)std::ceil((end_wave - start_wave) / mean_dispersion) + 2;
        
        std::vector<mosca::image> slit_flats;
        
        for(Iter flat_it = flat_start; flat_it != flat_end ; ++flat_it)
        {
            /* Trim the image to the slit area only */
            /* TODO: This assumes dispersion axis in X?? */
            mosca::image slit_flat = flat_it->trim(disp_bottom, spa_bottom,
                                                   disp_top, spa_top);
            //Multiply by the slit mask. Only valid pixels are not 0
            cpl_image_multiply(slit_flat.get_cpl_image(), slit_mask);
            cpl_image_multiply(slit_flat.get_cpl_image_err(), slit_mask);
            slit_flats.push_back(slit_flat);
        }
        
        /* Get the average flat. We don't care so much if areas
         * which don't belong to this slit are present */
        mosca::image average_slit_flat = 
                mosca::imagelist_reduce(slit_flats.begin(), 
                                        slit_flats.end(), 
                                        mosca::reduce_mean());

        /* Divide each individual flat by the average flat. 
         * This removes the pixel to pixel variation, therefore only differences
         * in large scale variations in wavelength are retained */
        std::vector<mosca::image> pix2pix_slit_flats = 
                slit_flats / average_slit_flat;

        /* All the seds of all the flats */
        std::vector<std::vector<double> > seds;
        
        /* Get the spectral energy distribution for each individual flat */ 
        for(std::vector<mosca::image>::iterator flat_it =  pix2pix_slit_flats.begin(); 
            flat_it != pix2pix_slit_flats.end(); ++flat_it)
        {

            mosca::image& slit_flat = *flat_it;
            T * flat_data = slit_flat.get_data<T>();
            T * flat_err = slit_flat.get_data_err<T>();
            cpl_mask * flat_mask = cpl_image_get_bpm(slit_flat.get_cpl_image());

            /* Compute the SED for this flat */
            std::vector<double> sed(n_pix_sed, 0.);
            std::vector<double> sed_err(n_pix_sed, 0.);
            std::vector<int> nsum(n_pix_sed, 0);
            for(cpl_size i = 0; i< flat_it->size_spatial(); ++i)
            {
                for(cpl_size j = 0; j< flat_it->size_dispersion(); ++j)
                {
                    /* When a pixel in column j is corrected from slit distortion, 
                     * it can happen that it corresponds to a column below or beyond
                     * the spectrum limits in the image. This margin_distor will
                     * ensure that +- that margin will be taken into account.
                     */
                    if(slit_mask_pix[j + slit_flat.size_dispersion() * i])
                    {

                        double spatial_corrected = slit_it->spatial_correct
                                ((double)(j + disp_bottom),
                                 (double)(i + spa_bottom));
                        double wavelength = wave_cal.get_wave
                                (spatial_corrected, (double)(j + disp_bottom));
                         
                        /* TODO: j, i depends on the spectral axis (vertical or horizontal) */
                        if(cpl_mask_get(flat_mask, j + 1, i + 1) == CPL_BINARY_0)
                        {
                            size_t idx_sed = (size_t)std::ceil((wavelength - start_wave) / mean_dispersion + 1.0);
                            if(idx_sed > 0 && idx_sed < n_pix_sed - 1)
                            {
                                sed[idx_sed] += flat_data[j + slit_flat.size_dispersion() * i];
                                sed_err[idx_sed] += flat_err[j + slit_flat.size_dispersion() * i];
                                nsum[idx_sed] += 1;
                            }
                        }
                    }
                }
            }
            for(int i_sed=0 ; i_sed< n_pix_sed; ++i_sed) 
                if(nsum[i_sed] == 0)
                    nsum[i_sed] = 1;
            mosca::vector_divide(sed, sed_err, nsum);

            /* Smooth the SED */
            if(smooth_size > 1 && smooth_size < sed.size())
                mosca::vector_smooth(sed, sed_err, smooth_size);
            
            /* Add this SED to the list of SEDs */
            seds.push_back(sed);
        }

        /* Divide the original flats by its SED */
        std::vector<mosca::image>::iterator flat_it;
        std::vector<std::vector<double> >::iterator sed_it;
        for(sed_it = seds.begin(), flat_it =  slit_flats.begin(); 
                flat_it != slit_flats.end(); ++flat_it, ++sed_it)
        {
            mosca::image& slit_flat = *flat_it;
            T * flat_data = slit_flat.get_data<T>();
            T * flat_err = slit_flat.get_data_err<T>();
            for(cpl_size j = 0; j< flat_it->size_dispersion(); ++j)
            {
                for(cpl_size i = 0; i< flat_it->size_spatial(); ++i)
                {
                    if(slit_mask_pix[j + slit_flat.size_dispersion() * i])
                    {
                        double spatial_corrected = slit_it->spatial_correct
                                ((double)(j + disp_bottom),
                                        (double)(i + spa_bottom));
                        double wavelength = wave_cal.get_wave
                                (spatial_corrected, (double)(j + disp_bottom));
                        size_t idx_sed = (size_t)std::ceil((wavelength - start_wave) / mean_dispersion + 1.0);
                        if(idx_sed > 0 && idx_sed < n_pix_sed - 1)
                        {
                            double sed_val = (*sed_it)[idx_sed]; 
                            if( sed_val != 0)
                            {
                                flat_data[j + slit_flat.size_dispersion() * i] /= sed_val;
                                /* TODO: Use the error in sed_err */
                                flat_err[j + slit_flat.size_dispersion() * i] /= sed_val;
                            }
                        }
                    }
                }
            }
        }

        /* Now we can stack the flats in the "usual" way. */
        mosca::image stacked_slit_flat_no_sed = 
                mosca::imagelist_reduce(slit_flats.begin(), 
                                        slit_flats.end(), comb_method);

        /* We compute the average SED */
        /* TODO: Error propagation here */
        /* TODO: there are nans */
        /* TODO: It is always 1, strange... */
        std::vector<double> avg_sed(n_pix_sed, 0.);
        for(size_t ipix = 0; ipix < n_pix_sed; ++ipix)
        {
            double sum = 0;
            for(size_t ised = 0; ised < seds.size(); ++ised)
                sum += seds[ised][ipix];
            avg_sed[ipix] = sum / seds.size();
        }

        /* Now we multiply the master slit flat by the average SED */
        double * average_slit_flat_no_sed_im = stacked_slit_flat_no_sed.get_data<double>();
        double * average_slit_flat_no_sed_err = stacked_slit_flat_no_sed.get_data_err<double>();
        for(cpl_size j = 0; j< stacked_slit_flat_no_sed.size_dispersion(); ++j)
        {
            for(cpl_size i = 0; i< stacked_slit_flat_no_sed.size_spatial(); ++i)
            {
                if(slit_mask_pix[j + stacked_slit_flat_no_sed.size_dispersion() * i])
                {
                    double spatial_corrected = slit_it->spatial_correct
                            ((double)(j + disp_bottom),
                             (double)(i + spa_bottom));
                    double wavelength = wave_cal.get_wave
                            (spatial_corrected, (double)(j + disp_bottom));
                    size_t idx_sed = (size_t)std::ceil((wavelength - start_wave) / mean_dispersion + 1.0);
                    if(idx_sed > 0 && idx_sed < n_pix_sed - 1)
                    {
                        average_slit_flat_no_sed_im[j + stacked_slit_flat_no_sed.size_dispersion() * i] *= avg_sed[idx_sed];
                        /* TODO: Use the error in sed_err */
                        average_slit_flat_no_sed_err[j + stacked_slit_flat_no_sed.size_dispersion() * i] *= avg_sed[idx_sed];
                    }
                }
            }
        }
        
        
        /* The master slit flat is placed in the master flat */
        cpl_image_copy(master_flat_im, stacked_slit_flat_no_sed.get_cpl_image(),
                       1, spa_bottom);
        cpl_image_copy(master_flat_err, stacked_slit_flat_no_sed.get_cpl_image_err(),
                       1, spa_bottom);
        
        /* Clean up */
        cpl_mask_delete(slit_mask_whole);
        cpl_mask_delete(trimmed_mask);
        cpl_image_delete(slit_mask); 
    }
    
    std::auto_ptr<mosca::image> 
        master_flat(new mosca::image(master_flat_im, master_flat_err, true));
    
    return master_flat;
}

/**
 * @brief
 *   Get a master spectroscopy flat 
 *
 * @param image_list A vector with all the flats.
 * @param image_end   An iterator with the end of the list of flat images.
 *
 * @tparam T    The type of data of the images: float, double, etc...
 * @tparam Iter The type of iterator. If dereferenced it should return a 
 *              mosca::image object
 *
 */
template<typename T, typename CombineMethod>
std::auto_ptr<mosca::image> mosca::flat_combine
(std::vector<mosca::image>& image_list, 
 const std::vector<mosca::calibrated_slit>& slits,
 const mosca::wavelength_calibration& wave_cal,
 const mosca::grism_config& grism_cfg,
 size_t smooth_size, CombineMethod comb_method)
{
    typedef std::vector<mosca::image>::iterator iter_type;
    return mosca::flat_combine<T, iter_type, CombineMethod >
        (image_list.begin(), image_list.end(), 
                slits, wave_cal, grism_cfg, smooth_size, comb_method);
}


#endif
