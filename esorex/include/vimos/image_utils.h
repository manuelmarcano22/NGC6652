/* $Id: image_utils.h,v 1.4 2013-08-13 12:55:26 cgarcia Exp $
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
 * $Date: 2013-08-13 12:55:26 $
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifndef IMAGE_UTILS_H
#define IMAGE_UTILS_H

#include "mosca_image.h"

namespace mosca
{
template<typename Iter, typename ReduceMethod>
mosca::image imagelist_reduce
(Iter image_start, Iter image_end, 
 ReduceMethod reduce_method = mosca::reduce_mean());
}

template<typename Container>
Container operator/ (Container& image_list, 
                     mosca::image& dividend);


#include "image_utils.tcc"

#endif
