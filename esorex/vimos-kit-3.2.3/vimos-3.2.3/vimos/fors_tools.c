/* $Id: fors_tools.c,v 1.1 2010-06-16 11:07:49 cizzo Exp $
 *
 * This file is part of the FORS Library
 * Copyright (C) 2002-2006 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
 * $Author: cizzo $
 * $Date: 2010-06-16 11:07:49 $
 * $Revision: 1.1 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fors_tools.h>

#include <cpl.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>

/*----------------------------------------------------------------------------*/
/**
 * @defgroup fors_tools       High level functions
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/**
 * @brief    Compute average airmass
 * @param    header               header to read from
 * @return   average airmass
 */
double
fors_get_airmass(const cpl_propertylist *header)
{
  double airmass_start, airmass_end;
  airmass_start = cpl_propertylist_get_double(header, "ESO TEL AIRM START");
  if (cpl_error_get_code()) {
      cpl_msg_error(cpl_func, "Could not read ESO TEL AIRM START from header");
      return -1; 
  }
  
  airmass_end = cpl_propertylist_get_double(header, "ESO TEL AIRM END");
  if (cpl_error_get_code()) {
      cpl_msg_error(cpl_func, "Could not read ESO TEL AIRM END from header");
      return -1; 
  }
  
  return 0.5 * (airmass_start + airmass_end);
}

/*----------------------------------------------------------------------------*/
/** 
 * @brief Same as cpl_tools_get_kth_double
 */
/*----------------------------------------------------------------------------*/
double fors_tools_get_kth_double(
    double  *   a,
    int         n,
    int         k)
{
    double x ;
    int    i, j, l, m ;

    cpl_ensure(a, CPL_ERROR_NULL_INPUT, 0.00) ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                double temp = a[i];
                a[i] = a[j];
                a[j] = temp;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}

/**@}*/
