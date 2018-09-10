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

/*
 * vimos_response.h
 *
 *  Created on: 2014 4 2
 *      Author: cgarcia
 */

#ifndef VIMOS_RESPONSE_H_
#define VIMOS_RESPONSE_H_

#include "cpl.h"
#include "vimos_detected_slits.h"

CPL_BEGIN_DECLS

cpl_table *vimos_compute_response
(cpl_image *spectra, cpl_image * mapped_flat_sed,
 cpl_propertylist * flat_sed_header, cpl_table *objects,
 double startwave, double dispersion, double gain,
 double exptime, cpl_table *ext_table, double airmass, cpl_table *flux_table,
 const std::vector<double>& ignored_waves,
 const std::vector<std::pair<double, double> >& ignored_wave_ranges,
 int nknots, int degree, cpl_table *& response_interp,
 double& flat_sed_norm_factor,
 const vimos::detected_slits& det_slits);


int vimos_science_correct_flat_sed
(cpl_image *spectra,  cpl_table *objects,
 cpl_image * mapped_flat_sed,
 cpl_propertylist * flat_sed_header,
 cpl_propertylist * specphot_header,
 const vimos::detected_slits& det_slits);


int vimos_science_correct_flat_sed_mapped
(cpl_image *spectra,  cpl_table *objects,
 cpl_image * mapped_flat_sed,
 cpl_propertylist * flat_sed_header,
 cpl_propertylist * specphot_header,
 const vimos::detected_slits& det_slits);


CPL_END_DECLS

#endif /* VIMOS_RESPONSE_H_ */
