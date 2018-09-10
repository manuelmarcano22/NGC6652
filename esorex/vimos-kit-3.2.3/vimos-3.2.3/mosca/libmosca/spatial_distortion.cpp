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
#include <cmath>
#include "spatial_distortion.h"

namespace mosca 
{


spatial_distortion::spatial_distortion()
{
}

spatial_distortion::~spatial_distortion()
{
}

/**
 * 
 * @param image
 * @param slits
 * @param polytraces
 * @param reference
 * @param start_wavelength
 * @param end_wavelength
 * @param dispersion
 * @return
 * 
 *  It assumes flux conservation
 */
#define STRETCH_FACTOR   (1.20)
cpl_image * spatial_distortion::m_calibrate_spatial
(cpl_image * image, cpl_table * slits, cpl_table * polytraces, 
 double reference, double start_wavelength, double end_wavelength, 
 double dispersion)
{
    /* Max order is 5 */
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;
    cpl_image     **exslit;
    cpl_image      *resampled;
    float          *sdata;
    float          *xdata;
    double          vtop, vbot, value;
    double          top, bot;
    double          ypos, yfra;
    double          factor;
    int             yint, ysize;
    int             nslits;
    int             npseudo;
    cpl_size        nx, ny;
    //int             order;
    int             i, j; 
    cpl_size        k;

    int       pixel_above, pixel_below, refpixel, start_pixel, end_pixel;



    nx = cpl_image_get_size_x(image);
    ny = cpl_image_get_size_y(image);
    sdata = (float*)cpl_image_get_data(image);


    nslits   = cpl_table_get_nrow(slits);
    //order    = cpl_table_get_ncol(polytraces) - 2;

    /*
     * The spatial resampling is performed for a certain number of 
     * pixels above and below the position of the reference wavelength:
     */

    pixel_above = STRETCH_FACTOR * (end_wavelength - reference) / dispersion;
    pixel_below = STRETCH_FACTOR * (reference - start_wavelength) / dispersion;

    exslit = (cpl_image**)cpl_calloc(nslits, sizeof(cpl_image *));

    for (i = 0; i < nslits; i++) 
    {

        /*
         * Note that the x coordinate of the reference pixels on the CCD
         * is taken arbitrarily at the top end of each slit. This wouldn't
         * be entirely correct in case of curved slits, or in presence of
         * heavy distortions: in such cases the spatial resampling is
         * really performed across a wide range of wavelengths. But
         * the lag between top and bottom spectral curvature models 
         * would introduce even in such cases negligible effects on
         * the spectral spatial resampling.
         */

        refpixel = cpl_table_get_double(slits, "xtop", i, NULL);

        start_pixel = refpixel - pixel_below;
        if (start_pixel < 0)
            start_pixel = 0;

        end_pixel = refpixel + pixel_above;
        if (end_pixel > nx)
            end_pixel = nx;

        /*
         * Recover from the table of spectral curvature coefficients
         * the curvature polynomials.
         */
        polytop = cpl_polynomial_new(1);
        polybot = cpl_polynomial_new(1);
        if(!m_get_curv_polynomials(polytraces, slits, i, polytop, polybot))
            return NULL;

        /*
         * Allocate image for current extracted slit
         */

        //We base the size of the corrected slit on the detected size of the slit
        double ytop = cpl_table_get_double(slits, "ytop", i, NULL);
        double ybot = cpl_table_get_double(slits, "ybottom", i, NULL);
        npseudo = std::ceil(ytop-ybot);

        if (npseudo < 1) {
            cpl_polynomial_delete(polytop);
            cpl_polynomial_delete(polybot);
            continue;
        }

        exslit[i] = cpl_image_new(nx, npseudo+1, CPL_TYPE_FLOAT);
        xdata = (float*)cpl_image_get_data(exslit[i]);

        /*
         * Write interpolated values to slit image.
         */

        for (j = start_pixel; j < end_pixel; j++)
        {
            top = cpl_polynomial_eval_1d(polytop, j, NULL);
            bot = cpl_polynomial_eval_1d(polybot, j, NULL);
            factor = (top-bot)/npseudo;
            for (k = 0; k <= npseudo; k++) {
                ypos = top - k*factor;
                yint = std::floor(ypos);
                yfra = ypos - yint;
                if (yint >= 0 && yint < ny-1) {
                    vtop = sdata[j + nx*yint];
                    vbot = sdata[j + nx*(yint+1)];
                    //This means that the top and bottom traces are crossing,
                    //which is physically impossible, so let's set it to 0. 
                    if(factor <= 0 )  
                        value = 0;
                    else if(vtop == FLT_MAX || vbot == FLT_MAX)
                        value = FLT_MAX;
                    else
                    {
                        value = vtop*(1-yfra) + vbot*yfra;
                        value *= factor;
                    }
                    xdata[j + nx*(npseudo-k)] = value;
                }
            }
        }
        cpl_polynomial_delete(polytop);
        cpl_polynomial_delete(polybot);
    }

    /*
     * Now all the slits images are copied to a single image
     */

    ysize = 0;
    for (i = 0; i < nslits; i++)
        if (exslit[i])
            ysize += cpl_image_get_size_y(exslit[i]);

    resampled = cpl_image_new(nx, ysize, CPL_TYPE_FLOAT);

    yint = -1;
    for (i = 0; i < nslits; i++) {
        if (exslit[i]) {
            yint += cpl_image_get_size_y(exslit[i]);
            cpl_image_copy(resampled, exslit[i], 1, ysize - yint);
            cpl_image_delete(exslit[i]);
        }
    }

    cpl_free(exslit);

    return resampled;
}

bool spatial_distortion::m_to_undistorted
(double spa_coord_distorted, 
 double disp_coord,
 double &spa_coord_undistorted,
 cpl_table *slits, cpl_table * polytraces)
{
    /* Max order is 5 */
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;
    int             npseudo;
    int             i_slit = -1; 


    /*
     * The spatial resampling is performed for a certain number of 
     * pixels above and below the position of the reference wavelength:
     */

    int slit_start;
    for (int i = 0; i < cpl_table_get_nrow(slits); i++) 
    {
        double ytop = cpl_table_get_double(slits, "ytop", i_slit, NULL);
        double ybot = cpl_table_get_double(slits, "ybottom", i_slit, NULL);
        slit_start = cpl_table_get_int(slits, "position", i, NULL);
        if(spa_coord_distorted > ybot && spa_coord_distorted < ytop)
        {
            i_slit = i;
            break;
        }
    }
    if (i_slit == -1)
        return false;

    /*
     * Recover from the table of spectral curvature coefficients
     * the curvature polynomials.
     */
    polytop = cpl_polynomial_new(1);
    polybot = cpl_polynomial_new(1);
    if(!m_get_curv_polynomials(polytraces, slits, i_slit, polytop, polybot))
        return false;

    //We base the size of the corrected slit on the detected size of the slit
    double ytop = cpl_table_get_double(slits, "ytop", i_slit, NULL);
    double ybot = cpl_table_get_double(slits, "ybottom", i_slit, NULL);
    npseudo = std::ceil(ytop-ybot);

    if (npseudo < 1) {
        return false;
    }


    double spa_slit_top = 
            cpl_polynomial_eval_1d(polytop, disp_coord, NULL);
    double spa_slit_bot = 
            cpl_polynomial_eval_1d(polybot, disp_coord, NULL);
    double spa_factor = (spa_slit_top-spa_slit_bot)/npseudo;

    double ydiff = spa_coord_distorted - spa_slit_bot;
    spa_coord_undistorted = ydiff / spa_factor + slit_start;
    
    //Cleanup
    cpl_polynomial_delete(polytop);
    cpl_polynomial_delete(polybot);

    return true;
}

bool spatial_distortion::m_to_distorted
(double spa_coord_undistorted, 
 double disp_coord,
 double &spa_coord_distorted,
 cpl_table *slits, cpl_table * polytraces)
{
    cpl_polynomial *polytop;
    cpl_polynomial *polybot;
    int             npseudo;
    cpl_size        i_slit = -1; 


    /*
     * The spatial resampling is performed for a certain number of 
     * pixels above and below the position of the reference wavelength:
     */

    int slit_start = 0;
    for (int i = 0; i < cpl_table_get_nrow(slits); i++) 
    {
        slit_start = cpl_table_get_int(slits, "position", i, NULL);
        if(std::floor(spa_coord_undistorted) >= slit_start)
        {
            i_slit = i;
            break;
        }
    }
   
    /*
     * Recover from the table of spectral curvature coefficients
     * the curvature polynomials.
     */
    polytop = cpl_polynomial_new(1);
    polybot = cpl_polynomial_new(1);
    if(!m_get_curv_polynomials(polytraces, slits, i_slit, polytop, polybot))
        return false;    

    //We base the size of the corrected slit on the detected size of the slit
    double ytop = cpl_table_get_double(slits, "ytop", i_slit, NULL);
    double ybot = cpl_table_get_double(slits, "ybottom", i_slit, NULL);
    npseudo = std::ceil(ytop-ybot);

    if (npseudo < 1) {
        return false;
    }

    double spa_slit_top = 
            cpl_polynomial_eval_1d(polytop, disp_coord, NULL);
    double spa_slit_bot = 
            cpl_polynomial_eval_1d(polybot, disp_coord, NULL);
    double spa_factor = (spa_slit_top-spa_slit_bot)/npseudo;


    double ydiff = spa_coord_undistorted - slit_start;
    spa_coord_distorted = ydiff * spa_factor + spa_slit_bot;
    
    //Cleanup
    cpl_polynomial_delete(polytop);
    cpl_polynomial_delete(polybot);

    return true;
}

bool spatial_distortion::m_get_curv_polynomials(cpl_table * polytraces, 
                                                cpl_table * slits,
                                                cpl_size i_slit, 
                                                cpl_polynomial * polytop, 
                                                cpl_polynomial * polybot)
{
    const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
    int             missing_top, missing_bot;
    cpl_size        k;
    int             null;
    double          coeff;
    int             order;

    order    = cpl_table_get_ncol(polytraces) - 2;

    missing_top = 0;
    for (k = 0; k <= order; k++) {
        coeff = cpl_table_get_double(polytraces, clab[k], 2*i_slit, &null);
        if (null) {
            cpl_polynomial_delete(polytop);
            missing_top = 1;
            break;
        }
        cpl_polynomial_set_coeff(polytop, &k, coeff);
    }

    missing_bot = 0;
    for (k = 0; k <= order; k++) {
        coeff = cpl_table_get_double(polytraces, clab[k], 2*i_slit+1, &null);
        if (null) {
            cpl_polynomial_delete(polybot);
            missing_bot = 1;
            break;
        }
        cpl_polynomial_set_coeff(polybot, &k, coeff);
    }

    if (missing_top && missing_bot) 
        return false;


    /*
     * In case just one of the two edges was not traced, the other
     * edge curvature model is duplicated and shifted to the other
     * end of the slit: better than nothing!
     * TODO: This probably should be done in the tracing class, rather than here
     */

    if (missing_top) {
        polytop = cpl_polynomial_duplicate(polybot);
        double ytop = cpl_table_get_double(slits, "ytop", i_slit, NULL);
        double ybot = cpl_table_get_double(slits, "ybottom", i_slit, NULL);
        k = 0;
        coeff = cpl_polynomial_get_coeff(polybot, &k);
        coeff += ytop - ybot;
        cpl_polynomial_set_coeff(polytop, &k, coeff);
    }

    if (missing_bot) {
        polybot = cpl_polynomial_duplicate(polytop);
        double ytop = cpl_table_get_double(slits, "ytop", i_slit, NULL);
        double ybot = cpl_table_get_double(slits, "ybottom", i_slit, NULL);
        k = 0;
        coeff = cpl_polynomial_get_coeff(polytop, &k);
        coeff -= ytop - ybot;
        cpl_polynomial_set_coeff(polybot, &k, coeff);
    }
    
    return true;
}

} /* namespace mosca */
