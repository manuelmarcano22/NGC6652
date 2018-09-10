/*

$Id: ap.h,v 1.2 2015/09/22 15:09:20 jim Exp $

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

#ifndef AP_H
#define AP_H

#include <cpl.h>

#define MAXBL       250000
#define NAREAL      8

typedef struct {
    intptr_t x;
    intptr_t y;
    float z;
    float zsm;
    intptr_t iobj;
} plstruct;

typedef struct {
    short int areal[NAREAL]; /* height above thresh of areal-prof cuts */
    intptr_t lsiz;           /* size of a line */
    intptr_t csiz;           /* size of a column */
    intptr_t maxip;          /* max no. of parents ever used. */
    intptr_t maxbl;          /* size of pixel-storage block stack */
    intptr_t maxpa;          /* size of parent-stack. */
    intptr_t ipnop;          /* parent-number-of-pixels, min size of image */
    intptr_t nimages;        /* count of images */
    intptr_t ipstack;        /* parent-name stack pointer */
    intptr_t ibstack;        /* pixel-block name stack pointer */
    float thresh;       /* threshold for image detection */
    float background;   /* background value */
    float sigma;        /* median background sigma */
    intptr_t multiply;       /* smoothing multiplication */
    float xintmin;      /* minimum intensity for consideration */
    intptr_t mulpix;         /* minimum size for considering multiple images */
    float areal_offset; /* offset in areal profile levels */
    float fconst;       /* Normalisation constant for areal profiles */
    float saturation;   /* saturation level from background analysis */
    intptr_t icrowd;         /* true if deblending routine is to be used */

    intptr_t *blink;         /* block-link array */
    intptr_t *bstack;        /* stack of pixel names */
    struct {            /* Image control block array */
        intptr_t first;      /* link to first data block */
        intptr_t last;       /* current last block   */
        intptr_t pnop;       /* Parent no. pixels (-1 = inactive) */
        intptr_t growing;
        intptr_t touch;      /* 0 = does not touch an edge */
        intptr_t pnbp;       /* Parent no of bad pixels */
    } *parent;

    short int *pstack;  /* stack of parent names */
    plstruct *plessey;  /* x,y,i storage array */
    short int *lastline;/* Parents on last line */

    cpl_image *inframe;  /* Pointer to original image data frame */
    cpl_image *conframe; /* Pointer to original confidence map frame */
    intptr_t xtnum;           /* Extension number to use */
    float *indata;       /* Pointer to original image data */
    float *confdata;       /* Pointer to original confidence map data */
    unsigned char *mflag; /* Pointer to mflag array for tracking merges */
    cpl_mask *opmask;   /* Object pixel mask */
    float rcore;        /* Core radius for aperture photometry */
    float filtfwhm;     /* FWHM of smoothing kernel in detection algorithm */
    plstruct *plarray;  /* Plessey structure workspace for passing data to 
                           various processing routines */
    intptr_t npl;            /* Size of the above */
    intptr_t npl_pix;        /* Number of pixels in the above structure */
    float fwhm;          /* Value of the seeing */
    
    struct {
        intptr_t nbx;        /* X dimension of background map */
        intptr_t nby;        /* Y dimension of background map */
        intptr_t nbsize;     /* Size of a side of background map cell */
        float **bvals;  /* Pointer to background map */
    } backmap;
} ap_t;

typedef struct {
    float x;            /* x position                           */
    float y;            /* y position                           */
    float total;        /* total integrated intensity           */
    intptr_t area;           /* image area in pixels                 */
    float peak;         /* peak image intensity above sky       */
    float xx;           /* 2nd moment x                         */
    float xy;           /* 2nd moment cross term                */
    float yy;           /* 2nd moment y                         */
    float ecc;          /* Eccentricity                         */
    intptr_t areal[NAREAL];  /* areal profile of image               */
} apmCat_t;

#endif

/*

$Log: ap.h,v $
Revision 1.2  2015/09/22 15:09:20  jim
Fixed guards and comments

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.4  2014/12/11 12:23:34  jim
new version

Revision 1.3  2014/04/09 09:09:51  jim
Detabbed

Revision 1.2  2014/03/26 15:25:19  jim
Modified for floating point confidence maps

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/

