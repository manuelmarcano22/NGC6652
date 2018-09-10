/*
 * This file is part of the CASU Pipeline utilities
 * Copyright (C) 2015,2016 European Southern Observatory
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

#include "imcore.h"
#include "floatmath.h"
#include "util.h"
#include <stdlib.h>

/* least-squares fit of order m polynomial to n data points */
/*---------------------------------------------------------------------------*/
/**
 * TODO MISSING DOXYGEN
 */


/*---------------------------------------------------------------------------*/
/** TODO: complete documentation
    \ingroup cataloguemodules
    \brief Work out the median seeing

    \par Name:
        imcore_polynm
    \par Purpose:
        Determine polynomial coefficients
    \par Description:
        Determine polynomial coefficients
    \par Language:
        C
    \param xdat
        x data points
    \param xcor
        TODO
    \param n
        Number of data

    \param polycf
        polynomial coefficients

    \param m
        number of coefficients

    \param ilim
        TODO
    \returns
        Nothing

    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/


extern void imcore_polynm (float xdat[], float xcor[], intptr_t n,
                           float polycf[], 
                           intptr_t m, intptr_t ilim) {
  double a[25][25], b[25], temp;
  intptr_t i, j, k;

/*   if(n < m)  bomboutx(1, "polynm: too few data points"); */
/*   if(m > 25) bomboutx(1, "polynm: order of polynomial too large"); */
/* CONSTANTS: 25, */
  /* clear arrays */
  for(i = 0; i < 25; i++) {
    b[i] = 0.0;
    for(j = 0; j < 25; j++) a[i][j] = 0.0;
  }

  /* cumulate sums */
  for(i = 0; i < n; i++) {
    for(k = 0; k < m; k++) {
      temp = 1.0;
      if(k+ilim != 0)temp = pow(xcor[i], (float) (k+ilim));
      b[k] += xdat[i]*temp;

      for(j = 0; j <= k; j++) {
        temp = 1.0;
        if(k+j+2*ilim != 0)temp = pow(xcor[i], (float) (k+j+2*ilim));
        a[j][k] += temp;
      }
    }
  }

  for(k = 1; k < m; k++) {
    for(j = 0; j < k; j++) a[k][j] = a[j][k];
  }

  /* solve linear equations */
  imcore_solve(a, b, m);

  for(i = 0; i < m; i++) polycf[i] = b[i];
}
