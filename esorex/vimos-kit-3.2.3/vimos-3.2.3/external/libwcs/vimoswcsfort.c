/*** File saoimage/wcslib/wcsf77.c
 *** June 2, 2000
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics

 * Module:	wcsfort.c (World Coordinate Systems)
 * Purpose:	Fortran interface to C WCS subroutines in wcsinit.c and wcs.c
 * Subroutine:	wcsinit_ (hstring, wcs) sets WCS structure from an image header
 * Subroutine:	wcsxinit_ (cra,cdec,secpix,xrpix,yrpix,nxpix,nypix,rotate,
		    equinox,epoch,proj)	sets WCS structure from arguments
 * Subroutine:	wcskinit_ (nxpix,nypix,ctype1,ctype2,crpix1,crpix2,crval1,
 *			   crval2, cd,cdelt1,cdelt2,crota,equinox,epoch)
 *		    sets a WCS structure from keyword-based arguments
 * Subroutine:	wcssize_ (wcs, cra, cdec, dra, ddec)
 *		    returns the image center and size in WCS units
 * Subroutine:	wcsdist_ (x1,y1,x2,y2,diff)
 *		    compute angular distance between ra/dec or lat/long
 * Subroutine:	wcsoutinit (wcs,coor) sets up the output coordinate system
 * Subroutine:	getwcsout_ (iwcs, coorsys)
 *		Return WCS output coordinate system used by pix2wcs
 * Subroutine:	wcsininit (wcs,coor) sets up the input coordinate system
 * Subroutine:	getwcsin_ (iwcs, coorsys)
 *		Return WCS output coordinate system used by pix2wcs
 * Subroutine:	getwcsim_ (iwcs, coorsys) Return WCS coordinate system of image 
 * Subroutine:	setwcslin_ (iwcs, mode)
 *		Set output string mode for LINEAR coordinates
 * Subroutine:	setwcsdeg_ (iwcs, mode)
 *		Set output string mode as decimal degrees if not zero
 * Subroutine:	pix2wcs_ (wcs,xpix,ypix,xpos,ypos)
 *		    pixel coordinates -> sky coordinates
 * Subroutine:	pix2wcst_ (wcs,xpix,ypix,wcstring,lstr)
 *		    pixels -> sky coordinate string
 * Subroutine:	wcs2pix_ (wcs,xpos,ypos,xpix,ypix)
 *		    sky coordinates -> pixel coordinates

 * Copyright:   2000 Smithsonian Astrophysical Observatory
 *              You may do anything you like with this file except remove
 *              this copyright.  The Smithsonian Astrophysical Observatory
 *              makes no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.
 */

#include <string.h>		/* strstr, NULL */
#include <stdio.h>		/* stderr */
#include <math.h>		/* stderr */
#include "vimoswcs.h"
#ifndef VMS
#include <stdlib.h>
#endif

static struct WorldCoor **pvimoswcs;
static int nvimoswcs = 0;

/* set up a WCS structure from a FITS or IRAF image header */

/* Call WCSINIT (HSTRING, IWCS)
 * where HSTRING is a Character string and IWCS is a returned integer
 * identifying the data structure which is created */

void
vimoswcsinit_ (hstring, ivimoswcs, nc)

char	*hstring;	/* Character string containing FITS header information
			   in the format <keyword>= <value> {/ <comment>} */
int	*ivimoswcs;		/* Pointer to wcs structure (returned) */
int	nc;		/* Number of characters in hstring (supplied by Fortran */
{
    int id;
    struct WorldCoor *tvimoswcs, *vimoswcsinit();

    tvimoswcs = vimoswcsinit (hstring);

    /* Set index to -1 if no WCS is found in header */
    if (tvimoswcs == NULL)
	*ivimoswcs = -1;

    /* Otherwise, use first available index into vector of pointers */
    else if (nvimoswcs == 0) {
	pvimoswcs = (struct WorldCoor **) calloc (10, sizeof (void *));
	nvimoswcs = 10;
	}
    else {
	for (id = 0; id < nvimoswcs; id++) {
	    if (pvimoswcs[id] == NULL) {
		*ivimoswcs = id;
		pvimoswcs[id] = tvimoswcs;
		break;
		}
	    }
	}
    return;
}

/* set up a WCS structure from subroutine arguments*/

void
vimoswcsxinit_ (ivimoswcs,cra,cdec,secpix,xrpix,yrpix,nxpix,nypix,rotate,equinox,
	   epoch,proj,np)

int	*ivimoswcs;		/* Pointer to wcs structure (returned) */
double	*cra;		/* Center right ascension in degrees */
double	*cdec;		/* Center declination in degrees */
double	*secpix;	/* Number of arcseconds per pixel */
double	*xrpix;		/* Reference pixel X coordinate */
double	*yrpix;		/* Reference pixel X coordinate */
int	*nxpix;		/* Number of pixels along x-axis */
int	*nypix;		/* Number of pixels along y-axis */
double	*rotate;	/* Rotation angle (clockwise positive) in degrees */
int	*equinox;	/* Equinox of coordinates, 1950 and 2000 supported */
double	*epoch;		/* Epoch of coordinates, used for FK4/FK5 conversion
			 * no effect if 0 */
char	*proj;		/* Projection */
int	np;		/* Length of projection (supplied by Fortran) */

{
    int id;
    struct WorldCoor *tvimoswcs;
    
    tvimoswcs = vimoswcsxinit (*cra,*cdec,*secpix,*xrpix,*yrpix,*nxpix,*nypix,
		     *rotate,*equinox,*epoch,proj);

    /* Set index to -1 if no WCS is found in header */
    if (tvimoswcs == NULL)
	*ivimoswcs = -1;
    else {
	for (id = 0; id < 10; id++) {
	    if (pvimoswcs[id] == NULL) {
		*ivimoswcs = id;
		pvimoswcs[id] = tvimoswcs;
		break;
		}
	    }
	}
    return;
}

/* set up a WCS structure from subroutine arguments matching FITS keywords*/

void
vimoswcskinit_ (ivimoswcs,naxis1,naxis2,ctype1,ctype2,crpix1,crpix2,crval1,crval2,cd,
	   cdelt1,cdelt2,crota,equinox,epoch, nc1,nc2)

int	*ivimoswcs;		/* Pointer to wcs structure (returned) */
int	*naxis1;	/* Number of pixels along x-axis */
int	*naxis2;	/* Number of pixels along y-axis */
char	*ctype1;	/* FITS WCS projection for axis 1 */
char	*ctype2;	/* FITS WCS projection for axis 2 */
double	*crpix1, *crpix2; /* Reference pixel coordinates */
double	*crval1, *crval2; /* Coordinates at reference pixel in degrees */
double	*cd;		/* Rotation matrix, used if all 4 elements nonzero */
double	*cdelt1, *cdelt2; /* scale in degrees/pixel, ignored if cd is not NULL */
double	*crota;		/* Rotation angle in degrees, ignored if cd is set */
int	*equinox;	/* Equinox of coordinates, 1950 and 2000 supported */
double	*epoch;		/* Epoch of coordinates, used for FK4/FK5 conversion
			 * no effect if 0 */
int	nc1, nc2;	/* Lengths of CTYPEs (supplied by Fortran) */

{
    int id;
    struct WorldCoor *tvimoswcs;

    tvimoswcs = vimoswcskinit (*naxis1,*naxis2,ctype1,ctype2,*crpix1,*crpix2,*crval1,
		     *crval2, cd, *cdelt1, *cdelt2, *crota, *equinox, *epoch);

    /* Set index to -1 if no WCS is found in header */
    if (tvimoswcs == NULL)
	*ivimoswcs = -1;
    else {
	for (id = 0; id < nvimoswcs; id++) {
	    if (pvimoswcs[id] == NULL) {
		*ivimoswcs = id;
		pvimoswcs[id] = tvimoswcs;
		break;
		}
	    }
	}
    return;
}


void
vimoswcsclose_ (ivimoswcs)

int	*ivimoswcs;		/* Index to WCS pointer array */
{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];
    if (vimoswcs != NULL) {
	vimoswcsfree (vimoswcs);
	pvimoswcs[*ivimoswcs] = NULL;
	}
    return;
}


/* Return RA and Dec of image center, plus size in RA and Dec */

void
vimoswcssize_ (ivimoswcs, cra, cdec, dra, ddec, radecsys)

int	*ivimoswcs;		/* Index to WCS pointer array */
double	*cra;		/* Right ascension of image center (rad) (returned) */
double	*cdec;		/* Declination of image center radg) (returned) */
double	*dra;		/* Half-width in radians (returned) */
double	*ddec;		/* Half-height in radians (returned) */
char	*radecsys;	/* Equinox (returned) */

{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */
    double width, height;

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    /* If there is no WCS defined, return all zeroes */
    if (vimoswcs == NULL) {
	(void)fprintf (stderr,"No WCS info available\n");
	*cra = 0.0;
	*cdec = 0.0;
	*dra = 0.0;
	*ddec = 0.0;
	}

    else {
	vimoswcsfull (vimoswcs, cra, cdec, &width, &height);
	*dra = 0.5 * width;
	*ddec = 0.5 * height;

	/* Get coordinate system from structure */
	strcpy (radecsys, vimoswcs->radecsys);
	}
    return;
}


/* Compute distance in degrees between two sky coordinates */

void vimoswcsdist_ (x1,y1,x2,y2,diff)

double	*x1,*y1;	/* (RA,Dec) or (Long,Lat) in degrees */
double	*x2,*y2;	/* (RA,Dec) or (Long,Lat) in degrees */
double	*diff;		/* Distance in degrees */

{
    double xr1, xr2, yr1, yr2;
    double pos1[3], pos2[3], w, cosb;
    int i;

    /* Convert two vectors to direction cosines */
    xr1 = degrad (*x1);
    yr1 = degrad (*y1);
    cosb = cos (yr1);
    pos1[0] = cos (xr1) * cosb;
    pos1[1] = sin (xr1) * cosb;
    pos1[2] = sin (yr1);

    xr2 = degrad (*x2);
    yr2 = degrad (*y2);
    cosb = cos (yr2);
    pos2[0] = cos (xr2) * cosb;
    pos2[1] = sin (xr2) * cosb;
    pos2[2] = sin (yr2);

    /* Modulus squared of half the difference vector */
    w = 0.0;
    for (i = 0; i < 3; i++) {
	w = w + (pos1[i] - pos2[i]) * (pos1[i] - pos2[i]);
	}
    w = w / 4.0;
    if (w > 1.0) w = 1.0;

    /* Angle beween the vectors */
    *diff = 2.0 * atan2 (sqrt (w), sqrt (1.0 - w));
    *diff = raddeg (*diff);
    return;
}


/* Initialize WCS output coordinate system used by pix2vimoswcs */

void
vimoswcsoutinit_ (ivimoswcs, coorsys, nc)

int	*ivimoswcs;		/* Index to WCS structure */
char	*coorsys;	/* Output coordinate system */
int	nc;		/* Length of coorsys (set by Fortran) */

{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	vimoswcsoutinit (vimoswcs, coorsys);

    return;
}


/* Return WCS output coordinate system used by pix2vimoswcs */

void
getvimoswcsout_ (ivimoswcs, coorsys, nc)

int	*ivimoswcs;		/* Index to WCS structure */
char	*coorsys;	/* Output coordinate system (returned) */
int	nc;		/* Length of coorsys (set by Fortran) */

{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	strncpy (coorsys, getvimoswcsout (vimoswcs), nc);

    return;
}


/* Initialize WCS input coordinate system used by vimoswcs2pix */

void
vimoswcsininit_ (ivimoswcs, coorsys, nc)

int	*ivimoswcs;		/* Index to WCS structure */
char	*coorsys;	/* Input coordinate system */
int	nc;		/* Length of coorsys (Set by Fortran) */

{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	vimoswcsininit (vimoswcs, coorsys);

    return;
}


/* Return WCS input coordinate system used by vimoswcs2pix */

void
getvimoswcsin_ (ivimoswcs, coorsys, nc)

int	*ivimoswcs;		/* Index to WCS structure */
char	*coorsys;	/* Input coordinate system (returned) */
int	nc;		/* Length of coorsys (Set by Fortran) */

{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	strncpy (coorsys, getvimoswcsin (vimoswcs), nc);

    return;
}


/* Return WCS coordinate system of image */

void
getvimoswcsim_ (ivimoswcs, coorsys, nc)

int	*ivimoswcs;		/* Index to WCS structure */
char	*coorsys;	/* Image coordinate system (returned) */
int	nc;		/* Length of coorsys (Set by Fortran) */

{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	strncpy (coorsys, getradecsys (vimoswcs), nc);

    return;
}


/* Set WCS output in degrees or hh:mm:ss dd:mm:ss */

void
setvimoswcsdeg_ (ivimoswcs, mode)

int	*ivimoswcs;		/* Index to WCS structure */
int	*mode;		/* mode = 0: h:m:s d:m:s
			   mode = 1: fractional degrees */
{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	vimoswcs->degout = *mode;

    return;
}


/* Set output string mode for LINEAR coordinates */

void
setvimoswcslin_ (ivimoswcs, mode)

int	*ivimoswcs;		/* World coordinate system structure */
int	*mode;		/* mode = 0: x y linear
			   mode = 1: x units x units
			   mode = 2: x y linear units */
{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	vimoswcs->linmode = *mode;

    return;
}


/* Convert pixels to sky coordinate string */

void
pix2vimoswcst_ (ivimoswcs,xpix,ypix,vimoswcstring,lstr)

int	*ivimoswcs;		/* World coordinate system structure */
double	*xpix, *ypix;	/* x and y image coordinates in pixels */
char	*vimoswcstring;	/* World coordinate string (returned) */
int	lstr;		/* Length of world coordinate string (set by Fortran) */

{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	pix2vimoswcst (vimoswcs, *xpix, *ypix, vimoswcstring, lstr);

    return;
}


/* Convert pixel coordinates to World Coordinates */

void
pix2vimoswcs_ (ivimoswcs, xpix, ypix, xpos, ypos)

int	*ivimoswcs;		/* World coordinate system structure */
double	*xpix, *ypix;	/* x and y image coordinates in pixels */
double	*xpos, *ypos;	/* RA and Dec in radians (returned) */
{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	pix2vimoswcs (vimoswcs, *xpix, *ypix, xpos, ypos);

    return;
}


/* Convert World Coordinates to pixel coordinates */

void
vimoswcs2pix_ (ivimoswcs, xpos, ypos, xpix, ypix, offscl)

int	*ivimoswcs;		/* World coordinate system structure */
double	*xpos, *ypos;	/* World coordinates in degrees */
double	*xpix, *ypix;	/* Image coordinates in pixels */
int	*offscl;
{
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */

    if (*ivimoswcs >= 0 && *ivimoswcs < nvimoswcs)
	vimoswcs = pvimoswcs[*ivimoswcs];
    else
	vimoswcs = pvimoswcs[0];

    if (vimoswcs != NULL)
	vimoswcs2pix (vimoswcs, *xpos, *ypos, xpix, ypix, offscl);

    return;
}
/* Jan  4 1996	new program
 * Jan 12 1996	Add WCSSET to set WCS without an image
 *
 * Jun 15 1998	rename vimoswcsf77.c to vimoswcsfort.c and move to libvimoswcs
 * Jun 15 1998	rename WCSSET WCSXINIT and update arguments
 * Jul  7 1998	Change setlinmode to setvimoswcslin_; add setvimoswcsdeg_
 *
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 21 1999	Add pix2wcst_(); drop unused variables after lint
 * Dec 10 1999	Add error handling for iwcs; document all subroutines
 *
 * Jun  2 2000	Fix WCS structure pointers
 */
