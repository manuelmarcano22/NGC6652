/*

$Id: classify.h,v 1.3 2015/09/22 15:09:20 jim Exp $

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

#ifndef CLASSIFY_H
#define CLASSIFY_H

/* Required information */

#include <cpl.h>
//#include "../casu_tfits.h"
#include "casu_utilfunctions.h"

/* Prototype */

extern int imcore_classify(casu_tfits *, cpl_propertylist *, float, int);

#endif

/*

$Log: classify.h,v $
Revision 1.3  2015/09/22 15:09:20  jim
Fixed guards and comments

Revision 1.2  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
