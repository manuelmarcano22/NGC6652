/* $Id: TwoDLinearWCS.cc,v 1.5 2013-08-07 15:28:51 cgarcia Exp $
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
 * $Date: 2013-08-07 15:28:51 $
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#include <cmath>
#include "two_d_linear_wcs.h"


#include "vimoswcs.h"
#include "vimoswcscat.h"


two_d_linear_wcs::two_d_linear_wcs()
{
    std::string  proj = "TAN";
    wcs_ = vimoswcsxinit(0.,0.,0.,0.,0.,0,0,0.,0,0.,
                         const_cast<char *>(proj.c_str()));
}

two_d_linear_wcs::two_d_linear_wcs
(const double& centralRA, const double& centralDEC, const double& plateScale,
 const cpl_size& dimX, const cpl_size& dimY, 
 const double& rotAngle, bool flip,
 const double& wcs_epoch, const int wcs_equinox)
{
    std::string  proj = "TAN";
    wcs_ = vimoswcsxinit(centralRA, centralDEC, plateScale, dimX/2., dimY/2., 
            dimX, dimY, rotAngle, wcs_epoch, wcs_equinox,
            const_cast<char *>(proj.c_str()));
    if(!flip) vimoswcsdeltset(wcs_, wcs_->cdelt[0],wcs_->cdelt[1],rotAngle);
    else      vimoswcsdeltset(wcs_,-wcs_->cdelt[0],wcs_->cdelt[1],rotAngle);
}

two_d_linear_wcs::~two_d_linear_wcs()
{
    if(wcs_ != NULL)
        vimoswcsfree(wcs_);
}

void two_d_linear_wcs::to_world
(const double& xPixel, const double& yPixel, double& ra, double& dec) const
{
    pix2vimoswcs(wcs_, xPixel, yPixel, &ra, &dec);
    ra /= 15.;
}

bool two_d_linear_wcs::to_pixel
(const double& ra, const double& dec, double& xPixel, double& yPixel) const
{
    int off;
    vimoswcs2pix(wcs_, ra * 15., dec, &xPixel, &yPixel, &off);
    if(off) return false;
    return true;
}

double two_d_linear_wcs::epoch() const
{
    return (wcs_->epoch);
}

double two_d_linear_wcs::equinox() const
{
    return (wcs_->equinox);
}

double two_d_linear_wcs::plate_scale() const
{
    double xPixSize=std::fabs(((*wcs_).xinc)*3600.);
    double yPixSize=std::fabs(((*wcs_).yinc)*3600.);
    return ((xPixSize + yPixSize) / 2.);
}

double two_d_linear_wcs::crpix1() const
{
    return wcs_->crpix[0];
}

double two_d_linear_wcs::crpix2() const
{
    return wcs_->crpix[1];
}

double two_d_linear_wcs::crval1() const
{
    return wcs_->crval[0];
}

double two_d_linear_wcs::crval2() const
{
    return wcs_->crval[1];
}

double two_d_linear_wcs::cdelt1() const
{
    return wcs_->cdelt[0];
}

double two_d_linear_wcs::cdelt2() const
{
    return wcs_->cdelt[1];
}

std::string two_d_linear_wcs::ctype1() const
{
    return wcs_->ctype[0];
}

std::string two_d_linear_wcs::ctype2() const
{
    return wcs_->ctype[1];
}

std::string two_d_linear_wcs::cunit1() const
{
    return wcs_->units[0];
}

std::string two_d_linear_wcs::cunit2() const
{
    return wcs_->units[1];
}

double * two_d_linear_wcs::cd_matrix() const
{
    return wcs_->cd;
}

double two_d_linear_wcs::orientation() const
{
    return wcs_->rot;
}
