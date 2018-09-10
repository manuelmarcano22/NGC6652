/*

$Id: solve.c,v 1.3 2015/09/22 15:09:20 jim Exp $

* This file is part of the CASU Pipeline utilities
* Copyright (C) 2015 European Southern Observatory
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
#include <stdio.h>
#include <stdlib.h>

#include "imcore.h"
#include "floatmath.h"
#include "util.h"

/* gauss elimination to solve ax=b */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Use Gauss-Jordan elimination to solve ax=b

    \par Name:
        imcore_solve
    \par Purpose:
        Use Gauss-Jordan elimination to solve ax=b
    \par Description:
        Use Gauss-Jordan elimination to solve ax=b
        Standard algorithm
    \par Language:
        C
    \param a
        Input matrix
    \param b
        Input vector
    \param m
        rank of matrix

    \returns
        Nothing

    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern void imcore_solve (double a[25][25], double b[25], intptr_t m) {
  double temp, big, pivot, rmax;
  intptr_t i, iu, j, k, l = 0, jl, ib, ir;

  iu = m-1;
  for(i = 0; i < iu; i++) {
    big = 0.0;

    /* find largest remaining term in ith column for pivot */
    for(k = i; k < m; k++) {
      rmax = fabs(a[i][k]);
      if(rmax > big) {
        big = rmax;
        l = k;
      }
    }

    /* check for non-zero term */
    if(big == 0.0) {
      for(ib = 0; ib < m; ib++) b[ib] = 0.0;
/*        fprintf(stderr, "solve: Zero determinant\n"); */
      return;
    }

    if(i != l) {
      /* switch rows */
      for(j = 0; j < m; j++) {
        temp    = a[j][i];
        a[j][i] = a[j][l];
        a[j][l] = temp;
      }
      temp = b[i];
      b[i] = b[l];
      b[l] = temp;
    }

    /* pivotal reduction */
    pivot = a[i][i];
    jl = i+1;

    for(j = jl; j < m; j++) {
      temp = a[i][j]/pivot;
      b[j] -= temp*b[i];
      for(k = i; k < m; k++) a[k][j] -= temp*a[k][i];
    }
  }

  /* back substitution for solution */
  for(i = 0; i < m; i++) {
    ir = m-1-i;
    if(a[ir][ir] != 0.0) {
      temp = b[ir];
      if(ir != m-1) {
        for(j = 1; j <= i; j++) {
          k = m-j;
          temp -= a[k][ir]*b[k];
        }
      }
      b[ir] = temp/a[ir][ir];
    }
    else
      b[ir] = 0.0;
  }
}

/*

$Log: solve.c,v $
Revision 1.3  2015/09/22 15:09:20  jim
Fixed guards and comments

Revision 1.2  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.2  2014/04/09 09:09:51  jim
Detabbed

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
