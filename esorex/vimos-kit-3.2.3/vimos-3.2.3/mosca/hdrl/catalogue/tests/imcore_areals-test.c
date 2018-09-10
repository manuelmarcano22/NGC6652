/* $Id: imcore_areals-test.c,v 1.1 2015/10/15 16:18:25 jim Exp $
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
//#include <casu_mods.h>
#include "../casu_utilfunctions.h"
#include "../imcore.h"

#define NT 117

int main(void) {
    int x[] = {398,399,400,397,398,399,400,401,402,403,396,397,398,399,400,
               401,402,403,404,395,396,397,398,399,400,401,402,403,404,405,
               395,396,397,398,399,400,401,402,403,404,405,395,396,397,398,
               399,400,401,402,403,404,405,395,396,397,398,399,400,401,402,
               403,404,405,406,395,396,397,398,399,400,401,402,403,404,405,
               395,396,397,398,399,400,401,402,403,404,405,395,396,397,398,
               399,400,401,402,403,404,405,396,397,398,399,400,401,402,403,
               404,397,398,399,400,401,402,403,398,399,400,401};
    int y[] = {394,394,394,395,395,395,395,395,395,395,396,396,396,396,396,
               396,396,396,396,397,397,397,397,397,397,397,397,397,397,397,
               398,398,398,398,398,398,398,398,398,398,398,399,399,399,399,
               399,399,399,399,399,399,399,400,400,400,400,400,400,400,400,
               400,400,400,400,401,401,401,401,401,401,401,401,401,401,401,
               402,402,402,402,402,402,402,402,402,402,402,403,403,403,403,
               403,403,403,403,403,403,403,404,404,404,404,404,404,404,404,
               404,405,405,405,405,405,405,405,406,406,406,406};
    float z[] = {8.87152,12.515,7.69699,10.8527,22.2509,21.7368,13.0388,
                 12.1853,17.1976,7.43948,15.2245,29.1964,37.9117,57.9371,
                 71.5542,57.1288,34.7726,15.5934,11.5374,15.995,21.3606,
                 60.4006,103.46,147.55,168.274,147.476,98.9157,51.7186,20.188,
                 3.04248,5.77832,49.3103,98.2057,187.557,268.353,310.638,
                 274.295,183.969,94.6933,47.9889,20.245,26.3758,59.1781,
                 152.389,275.916,395.107,450.251,397.53,272.322,147.053,54.767,
                 11.8971,13.3888,73.3689,165.899,298.455,449.707,493.25,441.585,
                 299.31,157.474,70.1224,15.5313,8.76074,20.7188,54.5798,141.249,
                 264.87,382.736,435.452,393.871,268.175,138.485,65.9307,28.7812,
                 19.379,36.6449,93.5458,186.823,270.95,305.093,260.879,183.683,
                 100.676,32.6281,16.6497,5.94965,17.8105,57.256,106.32,145.264,
                 164.271,137.093,88.9384,60.7841,31.8582,10.0435,4.69162,
                 15.2187,32.5385,61.0381,74.5399,67.3727,43.3964,25.0956,
                 16.7595,-0.37323,21.3832,19.2497,18.5883,9.37448,19.6048,
                 11.5006,13.0159,14.5852,13.66,-1.04889};
    ap_t ap;
    intptr_t i,iareal[NAREAL];

    /* Initialise */

    cpl_test_init(PACKAGE_BUGREPORT,CPL_MSG_WARNING);

    /* Set up apm structure */

    ap.lsiz = 2048;
    ap.csiz = 2048;
    ap.xtnum = 0;
    ap.inframe = NULL;
    ap.conframe = NULL;
    imcore_apinit(&ap);
    ap.npl_pix = NT;
    ap.plarray = cpl_realloc(ap.plarray,NT*sizeof(plstruct));
    for (i = 0; i < NT; i++) {
        ap.plarray[i].x = x[i];
        ap.plarray[i].y = y[i];
        ap.plarray[i].z = z[i];
        ap.plarray[i].zsm = z[i];
    }
    ap.xintmin = 0.0;
    ap.areal_offset = 3.47165;
    ap.thresh = 11.0936;
    ap.fconst = 1.4427;

    /* Run the test */

    imcore_areals(&ap,iareal);
    cpl_test_eq(iareal[0],104);
    cpl_test_eq(iareal[1],75);
    cpl_test_eq(iareal[2],63);
    cpl_test_eq(iareal[3],45);
    cpl_test_eq(iareal[4],25);
    cpl_test_eq(iareal[5],9);
    cpl_test_eq(iareal[6],0);
    cpl_test_eq(iareal[7],0);
        
    /* Get out of here */

    imcore_apclose(&ap);
    return(cpl_test_end(0));

}

/*

$Log: imcore_areals-test.c,v $
Revision 1.1  2015/10/15 16:18:25  jim
New


*/
