/* $Id: vimos_flat_normalise.cc,v 1.9 2013-10-24 16:44:34 cgarcia Exp $
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
 * $Date: 2013-10-24 16:44:34 $
 * $Revision: 1.9 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>
#include <functional>
#include "vimos_flat_normalise.h"
#include "moses.h"
#include "image_smooth.h"
#include "vector_utils.h"
#include "image_spline_fit.h"
#include "image_normalisation.h"

#define STRETCH_FACTOR   (1.20)


vimos::flat_normaliser::flat_normaliser() : m_normalisation_image()
{
}

vimos::flat_normaliser::~flat_normaliser()
{
}
    
/**
 * @brief 
 *   Normalise a flat field exposure
 *  
 * @param flat        Image containing the original flat field spectra
 * @param wave_cal    The wavelength calibration
 * @param spatial     Spatial calibration image
 * @param slits       Table with slits positions
 * @param polytraces  Coefficients of spectral curvature polynomials
 * @param blue        Start lambda to process
 * @param red         End lambda to process
 * @param dispersion  Mean spectral dispersion
 * @param spa_smooth_radius   Number of pixels for smoothing kernel in spatial axis
 * @param disp_smooth_radius  Number of pixels for smoothing kernel in dispersion axis
 * @param spa_fit_polyorder   Order for polynomial fit along spatial axis
 * @param disp_fit_knots      Number of knots in the spline fitting along dispersion axis
 * @param fit_threshold       values below this will be igonred in the fits 
 *  
 * @return The smoothed flat field exposure used for normalisation
 *
 * TODO: rewrite
 * The input @em flat frame should be already bias subtracted, and should
 * be oriented so that the dispersion direction is horizontal with @em blue
 * on the left and @em red on the right. The flat field spectra are spatially 
 * rectified, heavily smoothed, and then mapped back on the CCD. The original 
 * @em flat image is divided IN PLACE by its smoothed counterpart, which is 
 * also returned. If the polynomial @em polyorder is set to a negative number 
 * the smoothing consists of a linear fit along the spatial direction 
 * (excluding 3+3 pixels at the spectral edges), and by a median filtering 
 * along the dispersion direction using a window with the specified 
 * @em sradius; alternatively, if @em polyorder is not negative, the smoothing
 * will consist of a polynomial fitting of the illumination profile along
 * the dispersion direction, performed independently for each row of the
 * spatially remapped spectra.
 * TODO: Use vimos::detected_slits rather than cpl_table* slits
 */ 
int vimos::flat_normaliser::mos_normalise
(mosca::image& flat, const mosca::wavelength_calibration& wave_cal,
 cpl_image *spatial, 
 const std::vector<mosca::calibrated_slit>& calib_slits,
 cpl_table *slits, cpl_table *polytraces, 
 double blue, double red, 
 double dispersion,
 int spa_smooth_radius, int disp_smooth_radius, 
 int spa_fit_polyorder, int disp_fit_nknots,
 double fit_threshold)
{
    const char     *func = "mos_mosflat_normalise";

    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_image      *rectified;
    cpl_image      *smo_flat;
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;

    int            *slit_id;
    float          *sdata;
    float          *xdata;
    float          *wdata;
    double          vtop, vbot, value;
    double          top, bot;
    double          coeff;
    double          ytop, ybot;
    double          ypos;
    double          fvalue;
    int             ivalue;
    int             npseudo;

    int             pixel_above, pixel_below, refpixel, start_pixel, end_pixel;
    int             nx, ny;
    int             xlow, ylow, xhig, yhig;
    int             nslits;
    int            *position;
    int            *length;
    int             missing_top, missing_bot;
    int             order;
    int             null;
    int             i, j;
    cpl_size        k;

    
    /* For debug puposes only: cpl_image      *smo_rectified; */

    if (flat.get_cpl_image() == NULL || slits == NULL || polytraces == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 1;
    }
 
    if (dispersion <= 0.0) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return 1;
    }

    if (red - blue < dispersion) {
        cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
        return 1;
    }
    
    cpl_image * flat_im = flat.get_cpl_image();
    cpl_image * flat_err = flat.get_cpl_image_err();
    
    double reference = wave_cal.get_refwave();
    rectified = mos_spatial_calibration(flat_im, slits, polytraces, reference,
                                        blue, red, dispersion, 0, NULL);
    
    nx = cpl_image_get_size_x(rectified);
    ny = cpl_image_get_size_y(rectified);

    smo_flat = cpl_image_new(cpl_image_get_size_x(spatial), 
                             cpl_image_get_size_y(spatial), CPL_TYPE_FLOAT);
    wdata = cpl_image_get_data_float(smo_flat);

    nslits   = cpl_table_get_nrow(slits);
    order    = cpl_table_get_ncol(polytraces) - 2;
    position = cpl_table_get_data_int(slits, "position");
    length   = cpl_table_get_data_int(slits, "length");
    slit_id  = cpl_table_get_data_int(slits, "slit_id");

    /*
     * The spatial resampling is performed for a certain number of
     * pixels above and below the position of the reference wavelength:
     */
    
    pixel_above = (int)(STRETCH_FACTOR * (red - reference) / dispersion);
    pixel_below = (int)(STRETCH_FACTOR * (reference - blue) / dispersion);

    xlow = 1;
    xhig = nx;
    m_wave_profiles.clear();
    
    for (i = 0; i < nslits; i++) {

        if (length[i] == 0)
        {
            std::vector<float> empty_sed(nx, 0.);
            m_wave_profiles.push_back(empty_sed);
            m_wave_profiles_norm.push_back(1.);
            continue;
        }

        /*
         * We DON'T write:
         *
         * ylow = position[i];
         * yhig = ylow + length[i];
         *
         * because the cpl_image pixels are counted from 1, and because in 
         * cpl_image_extract() the coordinates of the last pixel are inclusive.
         */

        ylow = position[i] + 1;
        yhig = ylow + length[i] - 1;

        cpl_image * exslit_orig =
                cpl_image_extract(rectified, xlow, ylow, xhig, yhig);
        mosca::image im_exslit_orig(exslit_orig, true);

        //Get a mask of valid pixels of the slit in the rectified space
        std::vector<mosca::calibrated_slit>::const_iterator slit_it; 
        for(std::vector<mosca::calibrated_slit>::const_iterator it = calib_slits.begin();
            it != calib_slits.end(); ++it) //Life would be easier with a lambda...
        {
            if(it->slit_id() == slit_id[i])
            {
                slit_it = it; 
                break;
            }
        }
        cpl_mask * slit_mask_whole =
                slit_it->get_mask_valid(flat.dispersion_axis());
        cpl_image * slit_mask_im = cpl_image_new_from_mask(slit_mask_whole);
        cpl_image * slit_mask_im_d =  cpl_image_cast(slit_mask_im, CPL_TYPE_FLOAT);
        cpl_image * slit_mask_rectified = 
                mos_spatial_calibration(slit_mask_im_d, slits, polytraces, 
                                        reference, blue, red, dispersion, 
                                        0, NULL);
        //We don't consider fractional pixels yet
        cpl_image_threshold(slit_mask_rectified, 0.75, 1.25, 0., 1.);
        cpl_image * exslit_mask_rectified =
                cpl_image_extract(slit_mask_rectified, xlow, ylow, xhig, yhig);
        mosca::image slit_mask_im_mos(exslit_mask_rectified, true);
        cpl_image_delete(slit_mask_rectified);
        cpl_image_delete(slit_mask_im);
        cpl_image_delete(slit_mask_im_d);

        int  final_spa_smooth_radius = spa_smooth_radius;
        if (im_exslit_orig.size_spatial() / 2 < spa_smooth_radius)
        {
            final_spa_smooth_radius = im_exslit_orig.size_spatial() / 2;
            cpl_msg_warning(cpl_func, "Slit too narrow for requested "
                            "smoothing radius %d. Using %d", 
                            spa_smooth_radius, final_spa_smooth_radius);
        }

        
        std::vector<float> slit_spa_norm_profile;
        std::vector<float> slit_disp_norm_profile;
        mosca::image normslit = mosca::image_normalise(im_exslit_orig,
                                      slit_mask_im_mos,
                                      final_spa_smooth_radius, disp_smooth_radius,
                                      spa_fit_polyorder, disp_fit_nknots, 
                                      fit_threshold, slit_spa_norm_profile,
                                      slit_disp_norm_profile);

        //Get a position around the middle of the slit with a valid wave calib
        int middle_slit = get_middle_slit_valid_calib(wave_cal, yhig, ylow);
        
        //Get the pixel of the waveref at the middle of the slit
        double pix_waveref = wave_cal.get_pixel(middle_slit,
                                                wave_cal.get_refwave());
        int pix_waveref_left = (int)std::floor(pix_waveref);
        int pix_waveref_rigth = (int)std::ceil(pix_waveref);
        //Normalise by the value at waveref
        double prof_val_waveref = 1;
        if(pix_waveref_left>=0 && pix_waveref_rigth < flat.size_dispersion())
            prof_val_waveref = (slit_disp_norm_profile[pix_waveref_left] +
                                slit_disp_norm_profile[pix_waveref_rigth]) / 2.;
        for(size_t i_prof = 0; i_prof < slit_disp_norm_profile.size(); i_prof++)
            slit_disp_norm_profile[i_prof] /= prof_val_waveref;
        
        m_wave_profiles.push_back(slit_disp_norm_profile);
        m_wave_profiles_norm.push_back(prof_val_waveref);
        
        cpl_image * exslit = normslit.get_cpl_image();
        
        /*
         * Recover from the table of spectral curvature coefficients
         * the curvature polynomials.
         */

        refpixel = (int)(cpl_table_get_double(slits, "xtop", i, NULL));

        start_pixel = refpixel - pixel_below;
        if (start_pixel < 0)
            start_pixel = 0;

        end_pixel = refpixel + pixel_above;
        if (end_pixel > nx)
            end_pixel = nx;

        missing_top = 0;
        polytop = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i, &null);
            if (null) {
                cpl_polynomial_delete(polytop);
                missing_top = 1;
                break;
            }
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        missing_bot = 0;
        polybot = cpl_polynomial_new(1);
        for (k = 0; k <= order; k++) {
            coeff = cpl_table_get_double(polytraces, clab[k], 2*i+1, &null);
            if (null) {
                cpl_polynomial_delete(polybot);
                missing_bot = 1;
                break;
            }
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }

        if (missing_top && missing_bot) {
            cpl_msg_debug(func, "Slit %d was not traced: no extraction!",
                          slit_id[i]);
            continue;
        }

        /*
         * In case just one of the two edges was not traced, the other
         * edge curvature model is duplicated and shifted to the other
         * end of the slit: better than nothing!
         */

        if (missing_top) {
            cpl_msg_debug(func, "Upper edge of slit %d was not traced: "
                          "the spectral curvature of the lower edge "
                          "is used instead.", slit_id[i]);
            polytop = cpl_polynomial_duplicate(polybot);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polybot, &k);
            coeff += ytop - ybot;
            cpl_polynomial_set_coeff(polytop, &k, coeff);
        }

        if (missing_bot) {
            cpl_msg_debug(func, "Lower edge of slit %d was not traced: "
                          "the spectral curvature of the upper edge "
                          "is used instead.", slit_id[i]);
            polybot = cpl_polynomial_duplicate(polytop);
            ytop = cpl_table_get_double(slits, "ytop", i, NULL);
            ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
            k = 0;
            coeff = cpl_polynomial_get_coeff(polytop, &k);
            coeff -= ytop - ybot;
            cpl_polynomial_set_coeff(polybot, &k, coeff);
        }


        /*
         * Now map smoothed image to CCD.
         * Note that the npseudo value related to this slit is equal
         * to the number of spatial pseudo-pixels decreased by 1
         * (compare with function mos_spatial_calibration()).
         */

        nx = cpl_image_get_size_x(flat_im);
        ny = cpl_image_get_size_y(flat_im);

        sdata = cpl_image_get_data_float(spatial);
        cpl_binary * slit_mask_data = cpl_mask_get_data(slit_mask_whole);
        xdata = cpl_image_get_data_float(exslit);
        npseudo = cpl_image_get_size_y(exslit) - 1;

        /*
         * Write interpolated smoothed values to CCD image
         */
        int disp_bottom, spa_bottom, disp_top, spa_top;
        slit_it->get_extent_pix(disp_bottom, spa_bottom, disp_top, spa_top);
        for (j = start_pixel; j < end_pixel; j++) 
        {
            for (int yint = spa_bottom; yint <= spa_top; yint++) 
            {
                if (slit_mask_data[j + nx*yint] == CPL_BINARY_0)
                    continue;

                /*
                 * The line:
                 *     value = sdata[j + nx*yint];
                 * should be equivalent to:
                 *     value = npseudo*(top-yint)/(top-bot);
                 */


                value = sdata[j + nx*yint];   /* Spatial coordinate on rectified space */
                ivalue = std::floor(value);        /* Nearest spatial pixels:   */
                fvalue = value - ivalue;      /* ivalue and ivalue+1       */
                if (ivalue < npseudo && ivalue >= 0) {
                    vtop = xdata[j + nx*(npseudo-ivalue)];
                    vbot = xdata[j + nx*(npseudo-ivalue-1)];
                    wdata[j + nx*yint] = vtop*(1-fvalue) + vbot*fvalue;
                }
                else if (ivalue == npseudo)
                    wdata[j + nx*yint] = xdata[j]; //j+nx*(npseudo-npseudo)
                
            }
        }
        cpl_polynomial_delete(polytop);
        cpl_polynomial_delete(polybot);
        cpl_mask_delete(slit_mask_whole);
    }
    
    //TODO: What happens with inter-slit flux conservation. In principle the overall level
    //difference across slits should be conserved.

    cpl_image_delete(rectified);

    cpl_image_divide(flat_im, smo_flat);
    if(flat_err != NULL)
        cpl_image_divide(flat_err, smo_flat);

    m_normalisation_image = mosca::image(smo_flat, true);
    
    return 0;
}


//TODO: This is probably not the best place for this. It is a static function in
//order to be used directly by vimos_science_map_disp_profile
int vimos::flat_normaliser::get_middle_slit_valid_calib
(const mosca::wavelength_calibration& wave_cal, 
 int slit_end_pos, int slit_begin_pos)
{
    int slit_mean_pos = slit_begin_pos + (slit_end_pos - slit_begin_pos) / 2;
    int slit_pos_good_wavecal = -1;

    //Try to find a slit position with valid wavelength calibration
    //First from the middle to the top
    for(int i_pos = slit_mean_pos; i_pos <= slit_end_pos; i_pos++)
        if(wave_cal.has_valid_cal(i_pos))
        {
            slit_pos_good_wavecal = i_pos;
            break;
        }
    //If not successful, try from the middle to the bottom
    if(slit_pos_good_wavecal == -1)
        for(int i_pos = slit_mean_pos; i_pos >= slit_begin_pos; i_pos--)
            if(wave_cal.has_valid_cal(i_pos))
            {
                slit_pos_good_wavecal = i_pos;
                break;
            }
    if(slit_pos_good_wavecal == -1)
        throw std::runtime_error("Slit doesn't have any good wavelength calibration");
    return slit_pos_good_wavecal;    
}

const mosca::image& vimos::flat_normaliser::get_normalisation_image() const
{
    return m_normalisation_image;
}
    
const std::vector<std::vector<float> >& vimos::flat_normaliser::get_wave_profiles() const
{
    return m_wave_profiles;
}

std::vector<float> vimos::flat_normaliser::get_wave_profiles_norm
(double mflat_exptime,
 const std::vector<float>& slit_widths,
 const std::vector<float>& slit_lengths) const
{
    if(m_wave_profiles_norm.size() != slit_widths.size() ||
       m_wave_profiles_norm.size() != slit_lengths.size())
        throw std::invalid_argument("Vector sizes do not match");
    std::vector<float> wave_profiles_norm_scaled;
    
    for(size_t i = 0; i < m_wave_profiles_norm.size(); ++i)
    {
        float scale = mflat_exptime * slit_widths[i] * slit_lengths[i];
        if(scale == 0 )
            scale = 1;  //In some cases the length has been detected as 0.
        wave_profiles_norm_scaled.push_back(m_wave_profiles_norm[i] / scale);
    }

    return wave_profiles_norm_scaled;
}

cpl_image * vimos::flat_normaliser::get_wave_profiles_im() const
{
    cpl_image * wave_profiles_im = 
            cpl_image_new(m_wave_profiles.front().size(), m_wave_profiles.size(), 
                          CPL_TYPE_FLOAT);
    float * wave_profiles_data = cpl_image_get_data_float(wave_profiles_im);
    for(size_t i_slit = 0 ; i_slit < m_wave_profiles.size(); ++i_slit)
    {
        wave_profiles_data = 
            std::copy(m_wave_profiles[i_slit].begin(), 
                      m_wave_profiles[i_slit].end(), wave_profiles_data);
    }
    
    return wave_profiles_im;
}

cpl_image * vimos::flat_normaliser::get_wave_profiles_im_mapped
(const vimos::detected_slits& det_slits,
 const mosca::wavelength_calibration& wave_cal,
 double firstLambda, double lastLambda, double dispersion) const
{
    int nl = (lastLambda - firstLambda) / dispersion;

    cpl_image * mapped_flat_sed = 
            cpl_image_new(nl, m_wave_profiles.size(), CPL_TYPE_FLOAT);
    
    for(size_t i_slit = 0; i_slit < det_slits.size(); ++i_slit)
    {
        int slit_begin_pos = 
                det_slits[i_slit].get_position_spatial_corrected();
        int slit_end_pos = slit_begin_pos +  
                det_slits[i_slit].get_length_spatial_corrected();
        
        if(slit_begin_pos != -1)
        {
            int slit_pos_good_wavecal = vimos::flat_normaliser::get_middle_slit_valid_calib
                    (wave_cal, slit_end_pos, slit_begin_pos);
            for(cpl_size i_wave = 0; i_wave < nl; ++i_wave)
            {
                double wave = firstLambda + i_wave * dispersion;
                double pixel = wave_cal.get_pixel(slit_pos_good_wavecal, wave);
                int i_pix = std::ceil(pixel+0.5);
                if(i_pix>= 0 && i_pix <m_wave_profiles.front().size())
                    cpl_image_set(mapped_flat_sed, i_wave+1, i_slit+1,
                                  m_wave_profiles[i_slit][i_pix]);
            }
        }
    }
    
    return mapped_flat_sed;
}
