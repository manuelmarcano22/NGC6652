/* $Id: VimosUtils.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2012 European Southern Observatory
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
 * $Date: 2013-03-25 11:43:04 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VIMOSUTILS_H
#define VIMOSUTILS_H

#include <string>
#include "pilframeset.h"

/**
 * This function transforms a double that representa a sexagesimal number
 * like hhmmss or ggmmss.sss into a double like hh.hhhhh or gg.gggg. Negative 
 * values are properly converted.
 * @param xxmmss The input number in "sexagesimal format"
 * @return The result as xx.xxxx
 */
double vimos_utils_sexagesimal_to_double(double xxmmss);


/**
 * This function checks that all the frames have the same value of a given
 * keyword and returns the value of such a keyword. So far it works
 * with T=double,int. 
 * @param frames    The input frames
 * @param tag       The input frames will be filtered to contain only these tag
 * @param keyname   The name of the keyword.
 * @param keyvalue  The keyword value (returned)
 * @return true if all the headers have the same value of the keyword
 */
template<typename T> 
bool vimos_check_equal_keys
(PilSetOfFrames * frames, 
 const std::string& tag, 
 const std::string&  keyname,
 T& keyvalue);

#endif

#include "vimos_utils.tcc"
