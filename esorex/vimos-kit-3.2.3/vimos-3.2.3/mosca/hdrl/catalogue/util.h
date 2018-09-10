/*

$Id: util.h,v 1.2 2015/09/22 15:09:20 jim Exp $

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

#ifndef UTIL_H
#define UTIL_H

/* --Macros-- */

#undef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#undef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define ARRAYLEN(a) (sizeof(a)/sizeof((a)[0]))

#define NINT(a) ((intptr_t) ((a)+((a) < 0 ? -0.5 : 0.5)))

#define ARRAYELEM2(a, s, x, y) ((a)[(x)*(s) + (y)])

#endif
/*

$Log: util.h,v $
Revision 1.2  2015/09/22 15:09:20  jim
Fixed guards and comments

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
