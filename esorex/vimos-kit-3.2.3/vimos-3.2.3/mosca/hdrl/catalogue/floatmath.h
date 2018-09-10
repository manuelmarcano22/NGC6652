/*

$Id: floatmath.h,v 1.2 2015/09/22 15:09:20 jim Exp $

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

#ifndef FLOATMATH_H
#define FLOATMATH_H

#include <math.h>

/* Emulate 'float' versions of the math.h functions we use on systems which
 * don't have them */

/* sqrtf */
#define sqrtf(a) ((float) sqrt((double) (a)))

/* fabsf */
#define fabsf(a) ((float) fabs((double) (a)))

/* logf */
#define logf(a) ((float) log((double) (a)))

/* log10f */
#define log10f(a) ((float) log10((double) (a)))

/* expf */
#define expf(a) ((float) exp((double) (a)))

/* sinf */
#define sinf(a) ((float) sin((double) (a)))

/* cosf */
#define cosf(a) ((float) cos((double) (a)))

/* tanf */
#define tanf(a) ((float) tan((double) (a)))

/* asinf */
#define asinf(a) ((float) asin((double) (a)))

/* acosf */
#define acosf(a) ((float) acos((double) (a)))

/* atanf */
#define atanf(a) ((float) atan((double) (a)))

/* atan2f */
#define atan2f(x, y) ((float) atan2((double) (x), (double) (y)))

/* powf */
#define powf(x, y) ((float) pow((double) (x), (double) (y)))

#endif  /* __FLOATMATH_H__ */

/*

$Log: floatmath.h,v $
Revision 1.2  2015/09/22 15:09:20  jim
Fixed guards and comments

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
