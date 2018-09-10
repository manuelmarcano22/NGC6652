/* $Id: pilastroutils.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
 *
 * This file is part of the VIMOS pipeline library
 * Copyright (C) 2000-2004 European Southern Observatory
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: cizzo $
 * $Date: 2008-10-21 09:10:13 $
 * $Revision: 1.1.1.1 $
 * $Name: not supported by cvs2svn $
 */

#define _XOPEN_SOURCE
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pilmessages.h"
#include "pilastroutils.h"

#ifndef DEG_TO_RAD
#  define DEG_TO_RAD  ((M_PI) / 180.)
#endif

#ifndef SEC_TO_DEG
#  define SEC_TO_DEG  (15. / 3600.)
#endif

#ifndef SEC_TO_RAD
#  define SEC_TO_RAD  ((SEC_TO_DEG) * (DEG_TO_RAD))
#endif


/**
 * @defgroup pilAstroUtils pilAstroUtils
 *
 * The module @b pilAstroUtils provides a set of astronomy related utility
 * functions.
 */

/**@{*/

/*
 * @brief
 *   Compute the zenith distance for an observation.
 *
 * @return The function returns sec(z) if no error occurred, otherwise
 *   the return value is 0.
 *
 * @param hourAngle  Hour angle in radians.
 * @param delta      Declination in radians.
 * @param latitude   Latitude of the observatory in radians.
 *
 * The function computes the secans of the zenith distance for an
 * observation taken at an angle @em hourAngle from the meridian,
 * which can take values in the range extending from $-\pi$ to $\pi$,
 * and the declination @em delta with possible values between
 * $-0.5\pi$ and $0.5\pi$. The latitude @em latitude of the
 * observing site may take values in the range $0$ to $2\pi$.
 */

 static double pilAstroGetZDistance(double hourAngle, double delta,
					  double latitude)
{

  double p0 = sin(latitude) * sin(delta);
  double p1 = cos(latitude) * cos(delta);
  double z = p0 + cos(hourAngle) * p1;


  if (z < MIN_DIVISOR)
    return 0.;

  return 1. / z;
  
}


/*
 * @brief
 *   Compute approximated airmass value.
 *
 * @return The function returns the airmass.
 *
 * @param secZ  Secans of the zenith distance.
 *
 * The function uses the approximation given by Young and Irvine
 * (Young A. T., Irvine W. M., 1967, Astron. J. 72, 945) to compute
 * the airmass for a given sec(z) @em secZ. This approximation
 * takes into account atmosphere refraction and curvature, but is in
 * principle only valid at sea level.
 */

 static double pilAstroGetAirmassYoung(double secZ)
{

  return secZ * (1. - 0.0012 * (pow(secZ, 2) - 1.));

}


/**
 * @brief
 *   Compute the airmass of an observation.
 *
 * @return The function returns the computed average airmass. In case an
 *   error occurred the returned airmass is set to -1.
 *
 * @param alpha         Rightascension in degrees.
 * @param delta         Declination in degrees.
 * @param siderialTime  Local sidereal time in seconds elapsed since
 *                      siderial midnight.
 * @param exposureTime  Integration time in seconds.
 * @param latitude      Latitude of the observatory site in degrees.
 *
 * The function calculates the average airmass for the line of sight
 * given by the rightascension @em alpha and the declination
 * @em delta. The latitude @em latitude in degrees of the
 * observatory site and the local siderial time @em siderialTime
 * at observation start has to be given, as well as the duration of the
 * observation, i.e. the exposure time @em exposureTime.
 * 
 * The airmass is computed using the approximations of Young and Irvine
 * (Young A. T., Irvine W. M., 1967, Astron. J. 72, 945) and Stetson
 * (Stetson P., 1987, PASP 99, 191).
 */

double pilAstroComputeAirmass(double alpha, double delta, double siderialTime,
			      double exposureTime, double latitude)  
{

  const char fctid[] = "pilAstroComputeAirmass";


  /* Weights for Stetson's formula */

  const double weights[] = {1. / 6., 2. / 3., 1. / 6.};
  const int nweights = sizeof(weights) / sizeof(double);

  /*
   * Accuracy limit for airmass approximation (cf. Young A. T., Irvine W. M.,
   * 1967, Astron. J. 72, 945).
   */

  const double airmassUpperLimit = 4.;


  int i;

  double z;
  double hourAngle, airmass;


  /*
   * Compute hour angle of the observation in degrees.
   */

  hourAngle = siderialTime * SEC_TO_DEG - alpha;


  /*
   * Range adjustments. Angle between line of sight and the meridian
   * is needed.
   */

  if (hourAngle < -180.)
    hourAngle += 360.;

  if (hourAngle > 180.)
    hourAngle -= 360.;


  /* 
   * Convert angles from degrees to radians
   */
  
  delta *= DEG_TO_RAD;
  latitude *= DEG_TO_RAD;
  hourAngle *= DEG_TO_RAD;


  /*
   * Calculate airmass of the observation using the approximation given
   * by Young and Irvine (Young A. T., Irvine W. M., 1967, Astron. J. 72,
   * 945) for the individual airmass values. For finite exposure times
   * these airmass values are averaged using the weights given by Stetson
   * (Stetson P., 1987, PASP 99, 191)
   */

  z = pilAstroGetZDistance(hourAngle, delta, latitude);

  if (fabs(z) < MIN_DIVISOR) {
    pilMsgDebug(fctid, "Airmass computation failed. Object is below the "
		"horizon.");
    return -1.;
  }
    
  airmass = pilAstroGetAirmassYoung(z);

  if (exposureTime > 0.) {
    double timeStep = exposureTime / (nweights - 1) * SEC_TO_RAD;
    airmass *= weights[0];

    for (i = 1; i < nweights; i++) {
      z = pilAstroGetZDistance(hourAngle + i * timeStep, delta, latitude);

      if (fabs(z) < MIN_DIVISOR) {
	pilMsgDebug(fctid, "Airmass computation failed. Object is below the "
		    "horizon.");
	return -1.;
      }

      airmass += weights[i] * pilAstroGetAirmassYoung(z);
    }
  }


  if (airmass > airmassUpperLimit) 
    pilMsgWarning(fctid, "Airmass larger than %d", airmassUpperLimit);
 
  return airmass;

}
/**@}*/
