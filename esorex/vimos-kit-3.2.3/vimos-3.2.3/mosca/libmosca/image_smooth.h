/* $Id: image_smooth.h,v 1.6 2013-05-17 10:05:16 cgarcia Exp $
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
 * $Date: 2013-05-17 10:05:16 $
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifndef IMAGE_SMOOTH_H
#define IMAGE_SMOOTH_H

#include "mosca_image.h"

namespace mosca
{
template<typename T>
void image_smooth_1d_median(mosca::image& image, 
                            int half_width, mosca::axis smooth_axis);
}


#include "image_smooth.cpp"

#endif
