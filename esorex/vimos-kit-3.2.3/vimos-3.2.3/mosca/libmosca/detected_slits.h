/* $Id: detected_slits.h,v 1.3 2013-07-24 11:20:19 cgarcia Exp $
 *
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

/*
 * $Author: cgarcia $
 * $Date: 2013-07-24 11:20:19 $
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifndef DETECTED_SLITS_H
#define DETECTED_SLITS_H

#include <vector>
#include <string>
#include "detected_slit.h"

namespace mosca
{

typedef std::vector<detected_slit> detected_slits;

detected_slits detected_slits_load_fits(const std::string& fitsfile_slit_loc,
                                        const std::string& fitsfile_curv_coeff,
                                        int image_disp_size);

}

#endif
