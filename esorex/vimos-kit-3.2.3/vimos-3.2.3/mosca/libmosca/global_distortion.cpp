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
#include "global_distortion.h"

namespace mosca 
{

global_distortion::global_distortion(cpl_table * global_dist) 
{
    if(cpl_table_get_nrow(global_dist) != 13)
        throw std::invalid_argument("Global distortion table must contain 13 rows");
    m_global_dist = cpl_table_duplicate(global_dist);
}

global_distortion::global_distortion() :
        m_global_dist(NULL)
{
}

global_distortion::global_distortion(global_distortion& rhs)
{
    m_global_dist = cpl_table_duplicate(rhs.m_global_dist);
}

global_distortion::~global_distortion()
{
    if(m_global_dist != NULL)
        cpl_table_delete(m_global_dist);
}

cpl_image * global_distortion::calibrate_spatial
(cpl_image * image, cpl_table *slits, double reference, 
 double start_wavelength, double end_wavelength, double dispersion)
{
    if (image == NULL) 
        return NULL;

    if (dispersion <= 0.0) 
        return NULL;

    if (end_wavelength - start_wavelength < dispersion) 
        return NULL;
    cpl_table * polytraces = m_create_curv_coeff_table(slits);
    
    
    cpl_image * spatial_calibrated;
    spatial_calibrated = m_calibrate_spatial
            (image, slits, polytraces, reference, 
             start_wavelength, end_wavelength, dispersion);

    cpl_table_delete(polytraces);

    return spatial_calibrated;
}


bool global_distortion::to_distorted(double spa_coord_undistorted, 
                                         double disp_coord,
                                         double &spa_coord_distorted,
                                         cpl_table *slits)
{
    cpl_table * polytraces = m_create_curv_coeff_table(slits);
    bool status = m_to_distorted(spa_coord_undistorted, disp_coord, 
            spa_coord_distorted,
            slits, polytraces);
    cpl_table_delete(polytraces);
    return status;
}

bool global_distortion::to_undistorted(double spa_coord_distorted, 
                                           double disp_coord,
                                           double &spa_coord_undistorted,
                                           cpl_table *slits)
{
    cpl_table * polytraces = m_create_curv_coeff_table(slits);
    bool status = m_to_undistorted(spa_coord_undistorted, disp_coord, 
            spa_coord_distorted,
            slits, polytraces);
    cpl_table_delete(polytraces);
    return status;
}


/**
 * @brief
 *   Build the curvature coefficients table from a global distortions table
 *
 * @param global      Global distortions table
 * @param maskslits   Table with slits positions on mask
 * @param slits       Table with slits positions on CCD
 *
 * @return Curvature coefficients table
 *
 * The output curvature coefficients table has the same structure of the
 * output of the function @c mos_poly_trace(). The column "slit_id" is 
 * obtained from the "slit_id" column of the input @em maskslits table.
 * The coefficients columns are obtained as
 * @code
 *               c0 = poly7(mx, my)
 *               c1 = poly8(mx, my)
 *               c2 = poly9(mx, my)
 * @endcode
 * where polyX is the polynomial obtained from row X
 * of the input global distortions table, and (mx, my) are the 
 * coordinates of the slits ends listed in the input @em maskslits 
 * table. The slits that are completely outside the CCD are excluded
 * from the table.
 */
cpl_table * global_distortion::m_create_curv_coeff_table(cpl_table *slits)
{
    const char     *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};
                                                 /* Max order is 5 */

    cpl_polynomial *crv[3];
    cpl_vector     *point;
    cpl_table      *polytraces;
    double         *dpoint;
    double         *xtop;
    double         *ytop;
    double         *xbottom;
    double         *ybottom;
    int            *slit_id;
    int            *valid_id;
    int             nslits, nvalid;
    int             found;
    int             i, j, k;



    nslits  = cpl_table_get_nrow(slits);
    slit_id = cpl_table_get_data_int(slits, "slit_id");
    xtop    = cpl_table_get_data_double(slits, "xtop");
    ytop    = cpl_table_get_data_double(slits, "ytop");
    xbottom = cpl_table_get_data_double(slits, "xbottom");
    ybottom = cpl_table_get_data_double(slits, "ybottom");

    polytraces = cpl_table_new(2*nslits);
    cpl_table_new_column(polytraces, "slit_id", CPL_TYPE_INT);
    for (i = 0; i < 3; i++)
        cpl_table_new_column(polytraces, clab[i], CPL_TYPE_DOUBLE);

    crv[0] = m_read_polynomial_row(10);
    crv[1] = m_read_polynomial_row(11);
    crv[2] = m_read_polynomial_row(12);

    point = cpl_vector_new(2);
    dpoint = cpl_vector_get_data(point);

    for (i = 0; i < nslits; i++) {
        for (j = 0; j < 2; j++) {  /* For top and bottom trace of each slit */

            cpl_table_set_int(polytraces, "slit_id", 2*i+j, slit_id[i]);

            if (j) {
                dpoint[0] = xbottom[i];
                dpoint[1] = ybottom[i];                
            }
            else {
                dpoint[0] = xtop[i];
                dpoint[1] = ytop[i];                
            }

            for (k = 0; k < 3; k++)
                if (crv[j])
                {
                    cpl_table_set_double(polytraces, clab[k], 2*i+j,
                                         cpl_polynomial_eval(crv[k], point));
                }
        }
    }

    cpl_vector_delete(point);
    for (j = 0; j < 3; j++)
        cpl_polynomial_delete(crv[j]);

    /*
     * Eliminate slits which are _entirely_ outside the CCD
     */
 
    nvalid  = cpl_table_get_nrow(slits);
    valid_id = cpl_table_get_data_int(slits, "slit_id");
    cpl_table_unselect_all(polytraces);
    for (i = 0; i < nslits; i++) {
        found = 0;
        for (j = 0; j < nvalid; j++) {
            if (slit_id[i] == valid_id[j]) {
                found = 1;
                break;
            }
        }
        if (!found) {
            cpl_table_select_row(polytraces, 2*i);
            cpl_table_select_row(polytraces, 2*i + 1);
        }
    }
    cpl_table_erase_selected(polytraces);
 
    nslits = cpl_table_get_nrow(polytraces);

    return polytraces;
}

/*
 * The following two static functions are used to read and write from the 
 * global distortion table the different model components. Conventionally
 * the table consists of 6 columns and 10 rows. Each row is just ordered 
 * storage for model coefficients, and these functions guarantee that the
 * coefficients are read in and written out correctly, independent on their
 * physical meaning. The first 6 table rows are a description of the IDS
 * coefficients, followed by a row containing only the used reference 
 * wavelength. The remaining 3 are a description of the spectral curvature.
 * The first row is a description of coefficient c0, the second of coefficient
 * c1, etc., of the IDS. The 8th row is a description of coefficient c0,
 * the 9th of coefficient c1, etc., of the spectral curvature. All are
 * bivariate polynomialx on x,y mask coordinates. If the input table
 * to the write routine is NULL, it is allocated and initialised. Also
 * the input polynomial could be NULL, and nothing would be written to 
 * the table. If both pointers are NULL the function is basically a
 * constructor of the global distortion table.
 */
#define MAX_COLNAME      (80)
cpl_polynomial *global_distortion::m_read_polynomial_row(cpl_size row)
{
    cpl_polynomial *poly = NULL;
    cpl_size        p[2];
    cpl_size        degree = 2;
    int             null;
    double          coeff;

    char   name[MAX_COLNAME];


    for (p[0] = 0; p[0] <= degree; p[0]++) {
        for (p[1] = 0; p[1] <= degree - p[0]; p[1]++) {
            snprintf(name, MAX_COLNAME, "a%" CPL_SIZE_FORMAT"%" CPL_SIZE_FORMAT"", p[0], p[1]);
            coeff = cpl_table_get_double(m_global_dist, name, row, &null);
            if (null)
                continue;
            if (poly == NULL)
                poly = cpl_polynomial_new(2);
            cpl_polynomial_set_coeff(poly, p, coeff);
        }
    }

    return poly;
}

} /* namespace mosca */
