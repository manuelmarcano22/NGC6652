/*

$Id: imcore_radii.h,v 1.3 2015/09/22 15:09:20 jim Exp $

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

#ifndef IMCORE_RADII_H
#define IMCORE_RADII_H

#include "ap.h"

extern float imcore_halflight(float [], float [], float, float, intptr_t);
extern float imcore_exprad(float, float, float, float [], intptr_t);
extern float imcore_kronrad(float, float [], float [], intptr_t);
extern float imcore_petrad(float, float [], float [], intptr_t);
extern void imcore_flux(ap_t *, float [IMNUM][NPAR], intptr_t, float [], 
                             float [], intptr_t, float [], float []);

#endif

/*

$Log: imcore_radii.h,v $
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
