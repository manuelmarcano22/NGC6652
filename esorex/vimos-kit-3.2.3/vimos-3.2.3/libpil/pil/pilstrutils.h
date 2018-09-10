/* $Id: pilstrutils.h,v 1.2 2013-04-23 14:24:10 cgarcia Exp $
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
 * $Date: 2013-04-23 14:24:10 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef _PIL_STRUTILS_H
#define _PIL_STRUTILS_H

#include <stdlib.h>

#include <pilmacros.h>

 
PIL_BEGIN_DECLS

/*
 * Utilities returning a newly allocated string.
 */

char *pil_strdup(const char *);
extern char *pil_strndup(const char *, size_t);

PIL_END_DECLS

#endif /* _PIL_STRUTILS_H */
