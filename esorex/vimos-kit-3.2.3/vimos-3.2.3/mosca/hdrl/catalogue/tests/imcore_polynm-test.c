/* $Id: imcore_polynm-test.c,v 1.1 2015/10/15 16:18:25 jim Exp $
 *
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

/*
 * $Author: jim $
 * $Date: 2015/10/15 16:18:25 $
 * $Revision: 1.1 $
 * $Name:  $
 */

#include <stdio.h>
#include <stdlib.h>

#include <cpl_init.h>
#include <cpl_test.h>
//#include <casu_tfits.h>
//#include <casu_utils.h>
//#include <casu_mods.h>
#include "../casu_utilfunctions.h"

#include "../imcore.h"
#include "../imcore_radii.h"

#define NRADII 13

int main(void) {
    float x[] = {1.0,3.0,5.0,-10.0};
    float y[] = {-0.5,-27.5,-86.5,-424.0};
    float coefs[3] = {1.0,2.5,-4.0};
    float poly[3];

    /* Initialise */

    cpl_test_init(PACKAGE_BUGREPORT,CPL_MSG_WARNING);

    /* Do the test */

    imcore_polynm(y,x,4,poly,3,0);
    cpl_test_rel(poly[0],coefs[0],0.01);
    cpl_test_rel(poly[1],coefs[1],0.01);
    cpl_test_rel(poly[2],coefs[2],0.01);

    /* Get out of here */

    return(cpl_test_end(0));

}

/*

$Log: imcore_polynm-test.c,v $
Revision 1.1  2015/10/15 16:18:25  jim
New



*/
