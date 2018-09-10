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
 * vimos_detmodel.h
 *
 *  Created on: 2014 7 31
 *      Author: cgarcia
 */

#include "cpl.h"


#ifndef VIMOS_DETMODEL_H_
#define VIMOS_DETMODEL_H_

#ifdef __cplusplus

#include "ccd_config.h"

cpl_image * vimos_image_variance_from_detmodel
(cpl_image * image, const mosca::ccd_config& ccd_config);

#endif

//TODO: Remove this if all the recipes can call C++ code.
CPL_BEGIN_DECLS

cpl_image * vimos_image_variance_from_detmodel
(cpl_image * image, cpl_propertylist * header, cpl_propertylist * mbias_header);

CPL_END_DECLS

#endif /* VIMOS_DETMODEL_H_ */
