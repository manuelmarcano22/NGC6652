/*** File saoimage/wcslib/platepos.c
 *** February 29, 2000
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics

 * Module:	platepos.c (Plate solution WCS conversion
 * Purpose:	Compute WCS from plate fit
 * Subroutine:	platepos() converts from pixel location to RA,Dec 
 * Subroutine:	platepix() converts from RA,Dec to pixel location   

    These functions are based on the astrmcal.c portion of GETIMAGE by
    J. Doggett and the documentation distributed with the Digital Sky Survey.

*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "vimoswcs.h"

int
platepos (xpix, ypix, vimoswcs, xpos, ypos)

/* Routine to determine accurate position for pixel coordinates */
/* returns 0 if successful otherwise 1 = angle too large for projection; */
/* based on amdpos() from getimage */

/* Input: */
double	xpix;		/* x pixel number  (RA or long without rotation) */
double	ypix;		/* y pixel number  (dec or lat without rotation) */
struct WorldCoor *vimoswcs;	/* WCS parameter structure */

/* Output: */
double	*xpos;		/* Right ascension or longitude in degrees */
double	*ypos;		/* Declination or latitude in degrees */

{
    double x, y, x2, y2, x3, y3, r2;
    double xi, xir, eta, etar, raoff, ra, dec, ra0, dec0;
    double twopi = 6.28318530717959;
    double ctan, ccos;
    int ncoeff1 = vimoswcs->ncoeff1;
    int ncoeff2 = vimoswcs->ncoeff2;

    /*  Ignore magnitude and color terms 
    double mag = 0.0;
    double color = 0.0; */

    /* Convert from pixels to millimeters */
    x = xpix - vimoswcs->crpix[0];
    y = ypix - vimoswcs->crpix[1];
    x2 = x * x;
    y2 = y * y;
    x3 = x * x2;
    y3 = y * y2;
    r2 = x2 + y2;

    /*  Compute xi,eta coordinates in degrees from x,y and plate model */
    xi =  vimoswcs->x_coeff[ 0]		+ vimoswcs->x_coeff[ 1]*x +
	  vimoswcs->x_coeff[ 2]*y	+ vimoswcs->x_coeff[ 3]*x2 +
	  vimoswcs->x_coeff[ 4]*y2	+ vimoswcs->x_coeff[ 5]*x*y;

    if (ncoeff1 > 6)
	  xi = xi + vimoswcs->x_coeff[ 6]*x3	+ vimoswcs->x_coeff[ 7]*y3;

    if (ncoeff1 > 8) {
	xi = xi + vimoswcs->x_coeff[ 8]*x2*y	+ vimoswcs->x_coeff[ 9]*x*y2 +
		  vimoswcs->x_coeff[10]*(r2)	+ vimoswcs->x_coeff[11]*x*r2 +
		  vimoswcs->x_coeff[12]*y*r2;
	}

    eta = vimoswcs->y_coeff[ 0]		+ vimoswcs->y_coeff[ 1]*x +
	  vimoswcs->y_coeff[ 2]*y	+ vimoswcs->y_coeff[ 3]*x2 +
	  vimoswcs->y_coeff[ 4]*y2	+ vimoswcs->y_coeff[ 5]*x*y;

    if (ncoeff2 > 6)
	eta = eta + vimoswcs->y_coeff[ 6]*x3	+ vimoswcs->y_coeff[ 7]*y3;

    if (ncoeff2 > 8) {
	eta = eta + vimoswcs->y_coeff[ 8]*x2*y + vimoswcs->y_coeff[ 9]*y2*x +
		    vimoswcs->y_coeff[10]*r2   + vimoswcs->y_coeff[11]*x*r2 +
		    vimoswcs->y_coeff[12]*y*r2;
	}

    /* Convert to radians */
    xir = degrad (xi);
    etar = degrad (eta);

    /* Convert to RA and Dec */
    ra0 = degrad (vimoswcs->crval[0]);
    dec0 = degrad (vimoswcs->crval[1]);
    ctan = tan (dec0);
    ccos = cos (dec0);
    raoff = atan2 (xir / ccos, 1.0 - etar * ctan);
    ra = raoff + ra0;
    if (ra < 0.0) ra = ra + twopi;
    *xpos = raddeg (ra);

    dec = atan (cos (raoff) / ((1.0 - (etar * ctan)) / (etar + ctan)));
    *ypos = raddeg (dec);
    return 0;
}


int
platepix (xpos, ypos, vimoswcs, xpix, ypix)

/* Routine to determine pixel coordinates for sky position */
/* returns 0 if successful otherwise 1 = angle too large for projection; */
/* based on amdinv() from getimage */

/* Input: */
double	xpos;		/* Right ascension or longitude in degrees */
double	ypos;		/* Declination or latitude in degrees */
struct WorldCoor *vimoswcs;	/* WCS parameter structure */

/* Output: */
double	*xpix;		/* x pixel number  (RA or long without rotation) */
double	*ypix;		/* y pixel number  (dec or lat without rotation) */

{
    double xi,eta,x,y,xy,x2,y2,x2y,y2x,x3,y3,r2,dx,dy;
    double tdec,ctan,ccos,traoff, craoff, etar, xir;
    double f,fx,fy,g,gx,gy;
    double ra0, dec0, ra, dec;
    double tolerance = 0.0000005;
    int    max_iterations = 50;
    int    i;
    int	ncoeff1 = vimoswcs->ncoeff1;
    int	ncoeff2 = vimoswcs->ncoeff2;

    /* Convert RA and Dec in radians to standard coordinates on a plate */
    ra = degrad (xpos);
    dec = degrad (ypos);
    tdec = tan (dec);
    ra0 = degrad (vimoswcs->crval[0]);
    dec0 = degrad (vimoswcs->crval[1]);
    ctan = tan (dec0);
    ccos = cos (dec0);
    traoff = tan (ra - ra0);
    craoff = cos (ra - ra0);
    etar = (1.0 - ctan * craoff / tdec) / (ctan + (craoff / tdec));
    xir = traoff * ccos * (1.0 - (etar * ctan));
    xi = raddeg (xir);
    eta = raddeg (etar);

    /* Set initial value for x,y */
    x = xi * vimoswcs->dc[0] + eta * vimoswcs->dc[1];
    y = xi * vimoswcs->dc[2] + eta * vimoswcs->dc[3];

    /* if (vimoswcs->x_coeff[1] == 0.0)
	x = xi - vimoswcs->x_coeff[0];
    else
	x = (xi - vimoswcs->x_coeff[0]) / vimoswcs->x_coeff[1];
    if (vimoswcs->y_coeff[2] == 0.0)
	y = eta - vimoswcs->y_coeff[0];
    else
	y = (eta - vimoswcs->y_coeff[0]) / vimoswcs->y_coeff[2]; */

    /* Iterate by Newton's method */
    for (i = 0; i < max_iterations; i++) {

	/* X plate model */
	xy = x * y;
	x2 = x * x;
	y2 = y * y;
	x3 = x2 * x;
	y3 = y2 * y;
	x2y = x2 * y;
	y2x = y2 * x;
	r2 = x2 + y2;

	f = vimoswcs->x_coeff[0]	+ vimoswcs->x_coeff[1]*x +
	    vimoswcs->x_coeff[2]*y	+ vimoswcs->x_coeff[3]*x2 +
	    vimoswcs->x_coeff[4]*y2	+ vimoswcs->x_coeff[5]*xy;

	/*  Derivative of X model wrt x */
	fx = vimoswcs->x_coeff[1]	+ vimoswcs->x_coeff[3]*2.0*x +
	     vimoswcs->x_coeff[5]*y;

	/* Derivative of X model wrt y */
	fy = vimoswcs->x_coeff[2]	+ vimoswcs->x_coeff[4]*2.0*y +
	     vimoswcs->x_coeff[5]*x;

	if (ncoeff1 > 6) {
	    f = f + vimoswcs->x_coeff[6]*x3	+ vimoswcs->x_coeff[7]*y3;
	    fx = fx + vimoswcs->x_coeff[6]*3.0*x2;
	    fy = fy + vimoswcs->x_coeff[7]*3.0*y2;
	    }

	if (ncoeff1 > 8) {
	    f = f +
		vimoswcs->x_coeff[8]*x2y	+ vimoswcs->x_coeff[9]*y2x +
		vimoswcs->x_coeff[10]*r2 + vimoswcs->x_coeff[11]*x*r2 +
		vimoswcs->x_coeff[12]*y*r2;

	    fx = fx +	vimoswcs->x_coeff[8]*2.0*xy + 
			vimoswcs->x_coeff[9]*y2 +
	 		vimoswcs->x_coeff[10]*2.0*x +
			vimoswcs->x_coeff[11]*(3.0*x2+y2) +
			vimoswcs->x_coeff[12]*2.0*xy;

	    fy = fy +	vimoswcs->x_coeff[8]*x2 +
			vimoswcs->x_coeff[9]*2.0*xy +
			vimoswcs->x_coeff[10]*2.0*y +
			vimoswcs->x_coeff[11]*2.0*xy +
			vimoswcs->x_coeff[12]*(3.0*y2+x2);
	    }

	/* Y plate model */
	g = vimoswcs->y_coeff[0]	+ vimoswcs->y_coeff[1]*x +
	    vimoswcs->y_coeff[2]*y	+ vimoswcs->y_coeff[3]*x2 +
	    vimoswcs->y_coeff[4]*y2	+ vimoswcs->y_coeff[5]*xy;

	/* Derivative of Y model wrt x */
	gx = vimoswcs->y_coeff[1]	+ vimoswcs->y_coeff[3]*2.0*x +
	     vimoswcs->y_coeff[5]*y;

	/* Derivative of Y model wrt y */
	gy = vimoswcs->y_coeff[2]	+ vimoswcs->y_coeff[4]*2.0*y +
	     vimoswcs->y_coeff[5]*x;

	if (ncoeff2 > 6) {
	    g = g + vimoswcs->y_coeff[6]*x3	+ vimoswcs->y_coeff[7]*y3;
	    gx = gx + vimoswcs->y_coeff[6]*3.0*x2;
	    gy = gy + vimoswcs->y_coeff[7]*3.0*y2;
	    }

	if (ncoeff2 > 8) {
	    g = g +
		vimoswcs->y_coeff[8]*x2y	+ vimoswcs->y_coeff[9]*y2x +
		vimoswcs->y_coeff[10]*r2	+ vimoswcs->y_coeff[11]*x*r2 +
		vimoswcs->y_coeff[12]*y*r2;

	    gx = gx +	vimoswcs->y_coeff[8]*2.0*xy + 
			vimoswcs->y_coeff[9]*y2 +
	 		vimoswcs->y_coeff[10]*2.0*x +
			vimoswcs->y_coeff[11]*(3.0*x2+y2) +
			vimoswcs->y_coeff[12]*2.0*xy;

	    gy = gy +	vimoswcs->y_coeff[8]*x2 +
			vimoswcs->y_coeff[9]*2.0*xy +
			vimoswcs->y_coeff[10]*2.0*y +
			vimoswcs->y_coeff[11]*2.0*xy +
			vimoswcs->y_coeff[12]*(3.0*y2+x2);
	    }

	f = f - xi;
	g = g - eta;
	dx = ((-f * gy) + (g * fy)) / ((fx * gy) - (fy * gx));
	dy = ((-g * fx) + (f * gx)) / ((fx * gy) - (fy * gx));
	x = x + dx;
	y = y + dy;
	if ((fabs(dx) < tolerance) && (fabs(dy) < tolerance)) break;
	}

    /* Convert from plate pixels to image pixels */
    *xpix = x + vimoswcs->crpix[0];
    *ypix = y + vimoswcs->crpix[1];

    /* If position is off of the image, return offscale code */
    if (*xpix < 0.5 || *xpix > vimoswcs->nxpix+0.5)
	return -1;
    if (*ypix < 0.5 || *ypix > vimoswcs->nypix+0.5)
	return -1;

    return 0;
}


/* Set plate fit coefficients in structure from arguments */
int
SetPlate (vimoswcs, ncoeff1, ncoeff2, coeff)

struct WorldCoor *vimoswcs;  /* World coordinate system structure */
int	ncoeff1;		/* Number of coefficients for x */
int	ncoeff2;		/* Number of coefficients for y */
double	*coeff;		/* Plate fit coefficients */

{
    int i;

    if (novimoswcs (vimoswcs) || (ncoeff1 < 1 && ncoeff2 < 1))
	return 1;

    vimoswcs->ncoeff1 = ncoeff1;
    vimoswcs->ncoeff2 = ncoeff2;
    vimoswcs->prjcode = VIMOSWCS_PLT;

    for (i = 0; i < 20; i++) {
	if (i < ncoeff1)
	    vimoswcs->x_coeff[i] = coeff[i];
	else
	    vimoswcs->x_coeff[i] = 0.0;
	}

    for (i = 0; i < 20; i++) {
	if (i < ncoeff2)
	    vimoswcs->y_coeff[i] = coeff[ncoeff1+i];
	else
	    vimoswcs->y_coeff[i] = 0.0;
	}
    return 0;
}


/* Return plate fit coefficients from structure in arguments */
int
GetPlate (vimoswcs, ncoeff1, ncoeff2, coeff)

struct WorldCoor *vimoswcs;  /* World coordinate system structure */
int	*ncoeff1;	/* Number of coefficients for x */
int	*ncoeff2;	/* Number of coefficients for y) */
double	*coeff;		/* Plate fit coefficients */

{
    int i;

    if (novimoswcs (vimoswcs))
	return 1;

    *ncoeff1 = vimoswcs->ncoeff1;
    *ncoeff2 = vimoswcs->ncoeff2;

    for (i = 0; i < *ncoeff1; i++)
	coeff[i] = vimoswcs->x_coeff[i];

    for (i = 0; i < *ncoeff2; i++)
	coeff[*ncoeff1+i] = vimoswcs->y_coeff[i];

    return 0;
}


/* Set FITS header plate fit coefficients from structure */
void
SetFITSPlate (header, vimoswcs)

char    *header;        /* Image FITS header */
struct WorldCoor *vimoswcs;  /* WCS structure */

{
    char keyword[16];
    int i;

    for (i = 0; i < vimoswcs->ncoeff1; i++) {
	sprintf (keyword,"CO1_%d",i+1);
	hputnr8 (header, keyword, -15, vimoswcs->x_coeff[i]);
	}
    for (i = 0; i < vimoswcs->ncoeff2; i++) {
	sprintf (keyword,"CO2_%d",i+1);
	hputnr8 (header, keyword, -15, vimoswcs->y_coeff[i]);
	}
    return;
}

/* Mar 27 1998	New subroutines for direct image pixel <-> sky polynomials
 * Apr 10 1998	Make terms identical for both x and y polynomials
 * Apr 10 1998	Allow different numbers of coefficients for x and y
 * Apr 16 1998	Drom NCOEFF header parameter
 * Apr 28 1998  Change projection flags to WCS_*
 * Sep 10 1998	Check for xc1 and yc2 divide by zero after Allen Harris, SAO
 *
 * Oct 21 1999	Drop unused variables after lint
 *
 * Feb 29 2000	Use inverse CD matrix to get initial X,Y in platepix()
 *		as suggested by Paolo Montegriffo from Bologna Ast. Obs.
 */