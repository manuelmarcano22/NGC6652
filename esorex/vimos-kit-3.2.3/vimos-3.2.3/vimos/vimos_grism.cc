/*
 * This file is part of the VIMOS Data Reduction Pipeline
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdexcept>
#include "cpl.h"
#include "vimos_grism.h"

std::auto_ptr<mosca::grism_config> vimos_grism_config_from_table
(cpl_table * grism_table)
{
    std::auto_ptr<mosca::grism_config> grism_cfg;
    
    if (!cpl_table_has_column(grism_table, "dispersion") ||
        !cpl_table_has_column(grism_table, "reference") ||
        !cpl_table_has_column(grism_table, "startwavelength") ||
        !cpl_table_has_column(grism_table, "endwavelength") ) 
         throw std::invalid_argument("Table doesn't not contain "
                 "a grism configuration");
    
    
    if (cpl_table_get_column_type(grism_table, "dispersion") != CPL_TYPE_DOUBLE ||
        cpl_table_get_column_type(grism_table, "reference") != CPL_TYPE_DOUBLE ||
        cpl_table_get_column_type(grism_table, "startwavelength") != CPL_TYPE_DOUBLE ||
        cpl_table_get_column_type(grism_table, "endwavelength") != CPL_TYPE_DOUBLE) 
        throw std::invalid_argument("Unexpected type for GRISM_TABLE. "
                "Expected double");

    double nominal_dispersion = 
            cpl_table_get_double(grism_table,  "dispersion", 0, NULL);
    double wave_ref = 
            cpl_table_get_double(grism_table,  "reference", 0, NULL);
    double startwavelength = 
            cpl_table_get_double(grism_table,  "startwavelength", 0, NULL);
    double endwavelength = 
            cpl_table_get_double(grism_table,  "endwavelength", 0, NULL);

    grism_cfg.reset
     (new mosca::grism_config(nominal_dispersion,
                              startwavelength, endwavelength, wave_ref));
        
    return grism_cfg;
}
