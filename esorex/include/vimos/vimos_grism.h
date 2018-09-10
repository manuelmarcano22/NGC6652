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

#ifndef VIMOS_GRISM_H_
#define VIMOS_GRISM_H_

#include <memory>
#include "cpl.h"
#include "grism_config.h"

std::auto_ptr<mosca::grism_config> vimos_grism_config_from_table
(cpl_table * grism_table);

#endif /* VIMOS_GRISM_H_ */
