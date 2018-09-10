/*
 * This file is part of the HDRL
 * Copyright (C) 2016 European Southern Observatory
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

#ifndef HDRL_CATALOGUE_H
#define HDRL_CATALOGUE_H

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/
#include "hdrl_types.h"
#include "hdrl_image.h"
#include <cpl.h>
#include <float.h>

CPL_BEGIN_DECLS

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define HDRL_SATURATION_INIT INFINITY

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/


typedef struct {
    cpl_table * catalogue;
    cpl_image * segmentation_map;
    cpl_image * background;
    cpl_propertylist * qclist;
} hdrl_catalogue_result;

typedef enum {
    HDRL_CATALOGUE_BKG = 1,
    HDRL_CATALOGUE_SEGMAP = 2,
    HDRL_CATALOGUE_CAT_COMPLETE = 4,
    HDRL_CATALOGUE_ALL = 7
} hdrl_catalogue_options;

hdrl_catalogue_result *
hdrl_catalogue_compute(const cpl_image * image, const cpl_image * confidence_map,
                       const cpl_wcs * wcs, hdrl_parameter * param);
void hdrl_catalogue_result_delete(hdrl_catalogue_result * result);

hdrl_parameter * hdrl_catalogue_parameter_create(int obj_min_pixels,
                       double obj_threshold, cpl_boolean obj_deblending,
                       double obj_core_radius,
                       cpl_boolean bkg_estimate, int bkg_mesh_size,
                       double bkg_smooth_size, double det_eff_gain,
                       double det_saturation,
                       hdrl_catalogue_options resulttype);
cpl_error_code hdrl_catalogue_parameter_set_option(hdrl_parameter * par,
                                                   hdrl_catalogue_options opt);
cpl_boolean hdrl_catalogue_parameter_check(const hdrl_parameter * self);
hdrl_parameter * hdrl_catalogue_parameter_parse_parlist(
        const cpl_parameterlist *   parlist,
        const char              *   prefix);
cpl_parameterlist * hdrl_catalogue_parameter_create_parlist(
        const char      *   base_context,
        const char      *   prefix,
        hdrl_parameter  *   par_);


/*-----------------------------------------------------------------------------
             Private declarations - must not be used outside of hdrl
 -----------------------------------------------------------------------------*/

#ifdef HDRL_USE_PRIVATE

#endif

CPL_END_DECLS

#endif
