/* $Id: VimosUtils.cc,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <cmath>
#include "vimos_utils.h"

double vimos_utils_sexagesimal_to_double(double xxmmss)
{
    double double_number;

    double sign = std::fabs(xxmmss) / xxmmss;
    double xx = std::fabs((int)(xxmmss/10000));
    double mm = (std::fabs(xxmmss) - xx * 10000) / 100;
    double ss = (std::fabs(xxmmss) - xx * 10000 - mm * 100);

    double_number = sign * (xx + mm / 60 + ss / 3600);
    return double_number;
}
