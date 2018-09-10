/* $Id: pildate.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <string.h>
#include <time.h>

#define TIME_ISO8601_LENGTH (20)

/**
 * @defgroup pilDate pilDate
 *
 * The module @b pilDate provides functions related to time and date.
 */

/**@{*/

/**
 * @brief
 *   Get current date and time in ISO8601 format.
 *
 * @return Pointer to a statically allocated string in the function
 *   (no need to free it). In case of failure, returns a @c NULL.
 */

char *
pilDateGetISO8601(void)
{
  static char timeISO8601[TIME_ISO8601_LENGTH];
  time_t      seconds = time((time_t *)0);

  if (strftime(timeISO8601, TIME_ISO8601_LENGTH,
	       "%Y-%m-%dT%T", localtime(&seconds))) {
    return timeISO8601;
  }
  return NULL;
}
/**@}*/
