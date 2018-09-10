/* $Id: image_spline_fit.h,v 1.5 2013-09-06 17:23:47 cgarcia Exp $
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
 * $Date: 2013-09-06 17:23:47 $
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#ifndef IMAGE_SPLINE_FIT_H
#define IMAGE_SPLINE_FIT_H

#include "mosca_image.h"

namespace mosca
{
template<typename T>
void image_cubicspline_1d_fit(mosca::image& image, 
                              int nknots, double flux_threshold,
                              mosca::axis fitting_axis);

//TODO: Move it somewhere else
template<typename T>
void image_pol_1d_fit(mosca::image& image, int polyorder, 
                      double threshold);
}


#include "image_spline_fit.tcc"

#endif
