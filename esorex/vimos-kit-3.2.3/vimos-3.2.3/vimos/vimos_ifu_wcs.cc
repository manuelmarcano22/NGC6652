/* $Id: VimosIfuWCS.cc,v 1.3 2013/03/25 11:43:04 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2004 European Southern Observatory
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
 * $Date: 2013/03/25 11:43:04 $
 * $Revision: 1.3 $
 * $Name:  $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <cmath>
#include "cpl_math_const.h"
#include "cpl_msg.h"
#include "vimos_ifu_wcs.h"

two_d_linear_wcs vimos_ifu_get_2d_wcs_from_pointing
(double raVimos,  double decVimos, 
 double posAngle,
 cpl_size nxFibers, cpl_size nyFibers, double fiberSize, 
 double epoch, double equinox)
{
    /*   VIMOS IFU offset from the field center: */
    /*      667" along x (West  for PA=0, dec=0) */
    /*      1.7" along -y (South for PA=0, dec=0) */
    double a = 667.8; //This should come from ESO OCS CON WCS IFUCENX
    double b = -1.7; //This should come from ESO OCS CON WCS IFUCENY
    /*  Define IFU distance from the field center and angle */
    double di = std::sqrt(a * a + b * b);
    double oa = std::atan(b / a) * CPL_MATH_DEG_RAD; //degrees
    double ofs = (90. - di / 3600.) * CPL_MATH_RAD_DEG; //radians


    /*   first correct the position angle of the offset by adding the PA on sky of the OB */
    double paifu = -90. - oa + posAngle; //degrees
    /*   then transform the RA, DEC and corrected PA (paifu) in radiants: */
    double sra = raVimos * CPL_MATH_RAD_DEG; //radians
    double sdec = decVimos * CPL_MATH_RAD_DEG; //radians
    double ppa = paifu * CPL_MATH_RAD_DEG;//radians
    /*   transformation equations */
    double decifu = asin(cos(ppa) * cos(ofs) * cos(sdec) + sin(ofs) * sin(sdec));
    double ha = asin(-cos(ofs) * sin(ppa) / cos(decifu));
    double raifu = sra - ha;

    two_d_linear_wcs wcs = two_d_linear_wcs
            (raifu * CPL_MATH_DEG_RAD, decifu * CPL_MATH_DEG_RAD, fiberSize, 
            nxFibers ,nyFibers,
            posAngle,
            false, epoch, equinox);
    
    return wcs;
}
