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


#include <cpl.h>
#include <iostream>
#include <stdexcept>
#include "vimos_detected_slits.h"
#include "calibrated_slit.h"

std::vector<mosca::detected_slit>  
vimos::detected_slits_load_fits(const std::string& fitsfile_slit_loc,
                               const std::string& fitsfile_curv_coeff,
                               int image_size_disp)
{
    /* We read from extension 1 of slit location table */
    cpl_table* slit_loc_table = 
        cpl_table_load(fitsfile_slit_loc.c_str(), 1, 0);

    /* We read from extension 1 of curvature coefficients table */
    cpl_table* curv_coeff_table = 
        cpl_table_load(fitsfile_curv_coeff.c_str(), 1, 0);

    return vimos::detected_slits_from_tables(slit_loc_table, curv_coeff_table, 
                                            image_size_disp);
}

std::vector<mosca::detected_slit>  
vimos::detected_slits_from_tables(cpl_table * slit_loc,
                                 cpl_table * curv_coeff,
                                 int image_size_disp)
{
    std::vector<mosca::detected_slit> slits;
    
    if(2*cpl_table_get_nrow(slit_loc) != cpl_table_get_nrow(curv_coeff))
        throw std::invalid_argument("Slit and curv coeff Tables do not match");
    
    for(cpl_size irow = 0 ; irow < cpl_table_get_nrow(slit_loc); irow++)
    {
        int null;
        int slit_id = cpl_table_get_int(slit_loc, 
                                        "slit_id", irow, &null); 
        double disp_bottom = cpl_table_get_double(slit_loc,
                                               "xbottom", irow, &null); 
        disp_bottom = 1;
        double disp_top    = cpl_table_get_double(slit_loc,
                                               "xtop", irow, &null);
        disp_top = image_size_disp;
        double spa_bottom  = cpl_table_get_double(slit_loc, 
                                               "ybottom", irow, &null); 
        double spa_top     = cpl_table_get_double(slit_loc, 
                                               "ytop", irow, &null);
        
        int slit_id_check = cpl_table_get_int(curv_coeff, 
                                              "slit_id", 2*irow, &null);
        
        int position_spatial_corrected  = cpl_table_get_int(slit_loc, 
                                               "position", irow, &null); 
        int length_spatial_corrected     = cpl_table_get_int(slit_loc, 
                                               "length", irow, &null);
        
        if(slit_id != slit_id_check)
            throw std::runtime_error("Slit identification doesn't match");

        cpl_size ndegree = cpl_table_get_ncol(curv_coeff) - 1;

        std::vector<double> bottom_trace_coeff;
        std::vector<double> top_trace_coeff;
        for(cpl_size idx_coeff = 0; idx_coeff < ndegree; idx_coeff++)
        {
            char * colname = cpl_sprintf("c%"CPL_SIZE_FORMAT,idx_coeff);
            top_trace_coeff.push_back(cpl_table_get_double
                    (curv_coeff, colname, 2*irow, NULL));
            bottom_trace_coeff.push_back(cpl_table_get_double
                    (curv_coeff, colname, 2*irow + 1, NULL));
            cpl_free(colname);
        }
        
        
        slits.push_back(mosca::detected_slit(slit_id, disp_bottom, spa_bottom,
                                             disp_top, spa_top,
                                             position_spatial_corrected,
                                             length_spatial_corrected,
                                             bottom_trace_coeff, 
                                             top_trace_coeff));
    }
    
    return slits;
}

