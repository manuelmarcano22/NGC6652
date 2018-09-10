/* $Id: TwoDLinearWCS.h,v 1.3 2013/03/25 11:43:04 cgarcia Exp $
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
 * $Date: 2013/03/25 11:43:04 $
 * $Revision: 1.3 $
 * $Name:  $
 */

#ifndef TWODLINEARWCS_H
#define TWODLINEARWCS_H

#include <string>
#include "vimoswcs.h"
#include "cpl_type.h"

/**
 * This class represents a linear astrometric solution of an image.
 * Conversion from/to image coordinates to/from world coordinates are provided.
 */
class two_d_linear_wcs
{

public:

    /**
     * Default constructor. 
     */
    two_d_linear_wcs();

    /**
     * Constructor from a simple linear solution
     * @param centralRA RA of the center of the image(degrees)
     * @param centralDEC Dec of the center of the image(degrees)
     * @param PixelSize  Pixel size in arc seconds (i.e. plate scale)
     * @param dimX     Dimension in X direction (size of image)
     * @param dimY     Dimension in Y direction (size of image)
     * @param rotAngle Rotation angle wrt the north (clockwise positive)(degrees)
     * @param flip     There is a flip in the image
     * @param epoch    Epoch of the solution
     * @param equinox  Equinox to which the coordinates are refered.
     */
    two_d_linear_wcs(const double& centralRA, const double& centralDEC, 
                  const double& plateScale, const cpl_size& dimX, 
                  const cpl_size& dimY, 
                  const double& rotAngle = 0, bool flip = false,
                  const double& epoch = 2000.,const int equinox = 2000);

    /**
     * Destructor
     */
    ~two_d_linear_wcs();

    /**
     * Converts pixel positions to astrometric coordinates.
     * @param xPixel X coordinate in the image (pixel) 
     * @param yPixel Pixel in the image (pixel) 
     * @param ra     RA of the point (hours)
     * @param dec    Declination of the point (degrees)
     */
    void to_world
    (const double& xPixel, const double& yPixel, 
     double& ra, double& dec) const;

    /**
     * Converts astrometric coordinates to pixel positions.
     * @param ra     RA of the point (hours)
     * @param dec    Declination of the point (degrees)
     * @param xPixel X coordinate in the image (pixel) 
     * @param yPixel Pixel in the image (pixel)
     * @returns true if it success 
     */
    bool to_pixel
    (const double& ra, const double& dec, 
     double& xPixel, double& yPixel) const;

    /**
     * Returns the epoch of the astrometric solution
     */
    double epoch() const;

    /**
     * Returns the equinox of the astrometric solution
     */
    double equinox() const;

    /**
     * Returns the platescale
     */
    double plate_scale() const;

    /**
     * Returns the orientation in degrees of the solution (N through E).
     */
    double orientation() const;

    /**
     * Returns the CRPIX1 keyword of WCS FITS standard: 
     * X coordinate of reference pixel.
     */
    double crpix1() const;

    /**
     * Returns the CRPIX2 keyword of WCS FITS standard: 
     * Y coordinate of reference pixel.
     */
    double crpix2() const;

    /**
     * Returns the CDELT1 keyword of WCS FITS standard: 
     * Increment of WCS coord per X unit
     */
    double cdelt1() const;

    /**
     * Returns the CDELT2 keyword of WCS FITS standard: 
     * Increment of WCS coord per Y unit
     */
    double cdelt2() const;

    /**
     * Returns the CRVAL1 keyword of WCS FITS standard:
     * Value of X WCS coordinate for reference pixel.
     */
    double crval1() const;

    /**
     * Returns the CRVAL2 keyword of WCS FITS standard:
     * Value of Y WCS coordinate for reference pixel.
     */
    double crval2() const;

    /**
     * Returns the type of first WCS coordinate 
     */
    std::string ctype1() const;

    /**
     * Returns the type of second WCS coordinate 
     */
    std::string ctype2() const;

    /**
     * Returns the units of first WCS coordinate 
     */
    std::string cunit1() const;

    /**
     * Returns the units of second WCS coordinate 
     */
    std::string cunit2() const;

    /**
     * Returns the CD?_? matrix from WCS FITS standard
     */
    double * cd_matrix() const;

private:

    /**
     * Implementation of WCS using WCSLIB
     */
    struct WorldCoor *wcs_;

};

#endif
