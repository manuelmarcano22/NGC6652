/* $Id: pilstrutils.c,v 1.2 2013-04-23 14:26:36 cgarcia Exp $
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-04-23 14:26:36 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <assert.h>

#include "pilmemory.h"


/**
 * @defgroup pilstrutils String Utility Functions
 *
 * The module implements various string-related utility functions suitable
 * for creating, searching and modifying C strings.
 *
 * @par Synopsis:
 * @code
 *   #include <pilstrutils.h>
 * @endcode
 */

/**@{*/

/**
 * @brief
 *   Duplicate a string.
 *
 * @param string  String to be duplicated.
 *
 * @return Newly allocated copy of the original string.
 *
 * Duplicates the input string @em string. The newly allocated copy returned
 * to the caller can be deallocated using @b pil_free().
 */

char *
pil_strdup(const char *string)
{

#ifdef HAVE_STRDUP

    if (!string) {
        return NULL;
    }
    else {
        /* Cast needed because of xmemory */
        return strdup((char *)string);
    }

#else

    char *s;
    size_t sz;

    if (string) {
        sz = strlen(string) + 1;
        s = pil_malloc(sz * sizeof(char));
        memcpy(s, string, sz);
    }
    else
        s = NULL;

    return s;

#endif

}


/**
 * @brief
 *   Duplicate the first n charactes of a string.
 *
 * @param string  Source string
 * @param n       Maximum number of characters to be duplicated.
 *
 * @return Newly allocated copy of the first @em n characters of @em string.
 *
 * Duplicates the first @em n characters of the source string @em string,
 * returning the copied characters in newly allocated string of the size
 * @em n + 1. The returned string is always null terminated. If the length
 * of @em string is less than @em n the returned string is padded with nulls.
 * The newly allocated string can be deallocated using @b pil_free(). 
 */

 char *
pil_strndup(const char *string, size_t n)
{

    char *s;

    if (string) {
        s = pil_calloc(n + 1, sizeof(char));
        memcpy(s, string, n);
        s[n] = '\0';
    }
    else
        s = NULL;

    return s;

}
/**@}*/
