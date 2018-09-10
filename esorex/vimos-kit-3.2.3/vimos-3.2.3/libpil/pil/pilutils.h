/* $Id: pilutils.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifndef _PIL_UTILS_H
#define _PIL_UTILS_H

#include <stdlib.h>

#include <pilmacros.h>


PIL_BEGIN_DECLS

/*
 * String manipulating functions
 */

extern char *strlower(char *);
extern char *strupper(char *);
extern char *strtrim(char *, unsigned int);
extern char *strskip(char *, int (*)(int));
char *strsplit(char *, unsigned int, unsigned int);


/*
 * String checking functions
 */

extern int strempty(const char *, const char *);


/*
 * String selection functions
 */

extern int strselect(const char *, const char **, size_t);

PIL_END_DECLS

#endif /* _PIL_UTILS_H */
