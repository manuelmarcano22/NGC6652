/* $Id: pilmacros.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
 * $Author: cizzo $
 * $Date: 2008-10-21 09:10:13 $
 * $Revision: 1.1.1.1 $
 * $Name: not supported by cvs2svn $
 */

/*
 * This file MUST not include an other pipeline library header file!
 */

#ifndef _PIL_MACROS_H
#define _PIL_MACROS_H


/*
 * An alias for __extension__
 */

#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 8)
#  define PIL_GNUC_EXTENSION __extension__
#else
#  define PIL_GNUC_EXTENSION
#endif

/*
 * Wrap the gcc __PRETTY_FUNCTION__ and __FUNCTION__ variables with macros.
 */

#if defined (__GNUC__) && (__GNUC__ < 3)
#  define PIL_GNUC_FUNCTION         __FUNCTION__
#  define PIL_GNUC_PRETTY_FUNCTION  __PRETTY_FUNCTION__
#else /* !__GNUC__ */
#  define PIL_GNUC_FUNCTION         ""
#  define PIL_GNUC_PRETTY_FUNCTION  ""
#endif /* !__GNUC__ */

#define PIL_STRINGIFY(macro)         PIL_STRINGIFY_ARG(macro)
#define PIL_STRINGIFY_ARG(contents)  #contents

/*
 * String identifier for the current code position
 */

#if defined (__GNUC__) && (__GNUC__ < 3)
#  define PIL_CODE_POS  __FILE__ ":" PIL_STRINGIFY(__LINE__) ":" __PRETTY_FUNCTION__ "()"
#else
#  define PIL_CODE_POS  __FILE__ ":" PIL_STRINGIFY(__LINE__)
#endif

/*
 * C code guard
 */

#undef PIL_BEGIN_DECLS
#undef PIL_END_DECLS

#ifdef __cplusplus
#  define PIL_BEGIN_DECLS extern "C" {
#  define PIL_END_DECLS }
#else
#  define PIL_BEGIN_DECLS     /* empty */
#  define PIL_END_DECLS       /* empty */
#endif

/*
 * Commonly used macros. If they are provided by the system it is assumed
 * that the provided definition is correct.
 */

#ifndef NULL
#  ifdef __cplusplus
#    define NULL        (0L)
#  else /* !__cplusplus */
#    define NULL        ((void*) 0)
#  endif /* !__cplusplus */
#endif

#ifndef FALSE
#  define FALSE  (0)
#endif

#ifndef TRUE
#  define TRUE  (!FALSE)
#endif

#ifndef MAX
#  define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#  define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#  define ABS(a)     (((a) < 0) ? -(a) : (a))
#endif

#ifndef SQR
#  define SQR(a)     ((a) * (a))
#endif

/*
 * Number of elements in an array
 */

#define PIL_N_ELEMENTS(array)  (sizeof (array) / sizeof ((array)[0]))

/*
 * The maximum length, in bytes of an input line. It is set to the
 * POSIX2 value.
 */

#ifndef PIL_LINE_LENGTH_MAX
#  define PIL_LINE_LENGTH_MAX  2048
#endif

/*
 * Define the maximum number of characters a filename including the
 * path may have, excluding the trailing zero.
 */

#ifndef PIL_PATHNAME_MAX
#  define PIL_PATHNAME_MAX  4095
#endif

#endif /* _PIL_MACROS_H */
