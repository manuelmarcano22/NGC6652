/*** File libwcs/wcs.c
 *** March 22, 2001
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics

 * Module:	vimoswcs.c (World Coordinate Systems)
 * Purpose:	Convert FITS WCS to pixels and vice versa:
 * Subroutine:	vimoswcsxinit (cra,cdec,secpix,xrpix,yrpix,nxpix,nypix,rotate,equinox,epoch,proj)
 *		sets a WCS structure from arguments
 * Subroutine:	vimoswcskinit (nxpix,nypix,ctype1,ctype2,crpix1,crpix2,crval1,crval2,
		cd,cdelt1,cdelt2,crota,equinox,epoch)
 *		sets a WCS structure from keyword-based arguments
 * Subroutine:	vimoswcsreset (wcs,crpix1,crpix2,crval1,crval2,cdelt1,cdelt2,crota,cd, equinox)
 *		resets an existing WCS structure from arguments
 * Subroutine:	vimoswcsdeltset (wcs,cdelt1,cdelt2,crota) sets rotation and scaling
 * Subroutine:	vimoswcscdset (wcs, cd) sets rotation and scaling from a CD matrix
 * Subroutine:	vimoswcspcset (wcs,cdelt1,cdelt2,pc) sets rotation and scaling
 * Subroutine:	vimoswcseqset (wcs, equinox) resets an existing WCS structure to new equinox
 * Subroutine:	isvimoswcs(wcs) returns 1 if WCS structure is filled, else 0
 * Subroutine:	novimoswcs(wcs) returns 0 if WCS structure is filled, else 1
 * Subroutine:	vimoswcscent (wcs) prints the image center and size in WCS units
 * Subroutine:	wcssize (wcs, cra, cdec, dra, ddec) returns image center and size
 * Subroutine:	wcsfull (wcs, cra, cdec, width, height) returns image center and size
 * Subroutine:	wcsrange (wcs, ra1, ra2, dec1, dec2) returns image coordinate limits

 * Subroutine:	wcsshift (wcs,cra,cdec) resets the center of a WCS structure
 * Subroutine:	wcsdist (x1,y1,x2,y2) compute angular distance between ra/dec or lat/long
 * Subroutine:	wcsdiff (x1,y1,x2,y2) compute angular distance between ra/dec or lat/long
 * Subroutine:	wcscominit (wcs,command) sets up a command format for execution by wcscom
 * Subroutine:	wcsoutinit (wcs,coor) sets up the coordinate system used by pix2wcs
 * Subroutine:	getwcsout (wcs) returns current output coordinate system used by pix2wcs
 * Subroutine:	wcsininit (wcs,coor) sets up the coordinate system used by wcs2pix
 * Subroutine:	getwcsin (wcs) returns current input coordinate system used by wcs2pix
 * Subroutine:	setwcsdeg(wcs, new) sets WCS output in degrees or hh:mm:ss
 * Subroutine:	getradecsys(wcs) returns current coordinate system type
 * Subroutine:	wcscom (wcs,file,x,y,wcstr) executes a command using the current world coordinates
 * Subroutine:	setwcslin (wcs, mode) sets output string mode for LINEAR
 * Subroutine:	pix2wcst (wcs,xpix,ypix,wcstring,lstr) pixels -> sky coordinate string
 * Subroutine:	pix2wcs (wcs,xpix,ypix,xpos,ypos) pixel coordinates -> sky coordinates
 * Subroutine:	wcsc2pix (wcs,xpos,ypos,coorsys,xpix,ypix,offscl) sky coordinates -> pixel coordinates
 * Subroutine:	wcs2pix (wcs,xpos,ypos,xpix,ypix,offscl) sky coordinates -> pixel coordinates
 * Subroutine:  wcszin (izpix) sets third dimension for pix2wcs() and pix2wcst()
 * Subroutine:  wcszout (wcs) returns third dimension from wcs2pix()
 * Subroutine:	setwcsfile (filename)  Set file name for error messages 
 * Subroutine:	setwcserr (errmsg)  Set error message 
 * Subroutine:	wcserr()  Print error message 
 * Subroutine:	setdefwcs (wcsproj)  Set flag to choose AIPS or WCSLIB WCS subroutines 
 * Subroutine:	getdefwcs()  Get flag to switch between AIPS and WCSLIB subroutines 
 * Subroutine:	savewcscoor (wcscoor)
 * Subroutine:	getwcscoor()  Return preset output default coordinate system 
 * Subroutine:	savewcscom (i, wcscom)  Save specified WCS command 
 * Subroutine:	setwcscom (wcs)  Initialize WCS commands 
 * Subroutine:	getwcscom (i)  Return specified WCS command 
 * Subroutine:	wcsfree (wcs)  Free storage used by WCS structure
 * Subroutine:	freewcscom (wcs)  Free storage used by WCS commands 

 * Copyright:   2001 Smithsonian Astrophysical Observatory
 *              You may do anything you like with this file except remove
 *              this copyright.  The Smithsonian Astrophysical Observatory
 *              makes no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.

 */

#include <string.h>		/* strstr, NULL */
#include <stdio.h>		/* stderr */
#include <math.h>
#include "vimoswcs.h"
#include "fitshead.h"
#ifndef VMS
#include <stdlib.h>
#endif

static char vimoswcserrmsg[80];
static char vimoswcsfile[256]={""};
static void vimoswcslibrot();
void vimoswcsrotset();
static int vimoswcsproj0 = 0;
static int izpix = 0;
static double zpix = 0.0;

void
vimoswcsfree (
struct WorldCoor *vimoswcs)	/* WCS structure */
{
    if (novimoswcs (vimoswcs)) {

	/* Free WCS structure if allocated but not filled */
	if (vimoswcs)
	    free (vimoswcs);

	return;
	}

    freevimoswcscom (vimoswcs);
    if (vimoswcs->lin.imgpix != NULL)
	free (vimoswcs->lin.imgpix);
    if (vimoswcs->lin.piximg != NULL)
	free (vimoswcs->lin.piximg);
    free (vimoswcs);
    return;
}

/* Set up a WCS structure from subroutine arguments */

struct WorldCoor *
vimoswcsxinit (
double	cra,	/* Center right ascension in degrees */
double	cdec,	/* Center declination in degrees */
double	secpix,	/* Number of arcseconds per pixel */
double	xrpix,	/* Reference pixel X coordinate */
double	yrpix,	/* Reference pixel X coordinate */
int	nxpix,	/* Number of pixels along x-axis */
int	nypix,	/* Number of pixels along y-axis */
double	rotate,	/* Rotation angle (clockwise positive) in degrees */
int	equinox, /* Equinox of coordinates, 1950 and 2000 supported */
double	epoch,	/* Epoch of coordinates, used for FK4/FK5 conversion
		 * no effect if 0 */
char	*proj)	/* Projection */

{
    struct WorldCoor *vimoswcs;
    double cdelt1, cdelt2;

    vimoswcs = (struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

    /* Set WCSLIB flags so that structures will be reinitialized */
    vimoswcs->cel.flag = 0;
    vimoswcs->lin.flag = 0;
    vimoswcs->vimoswcsl.flag = 0;

    /* Image dimensions */
    vimoswcs->naxes = 2;
    vimoswcs->lin.naxis = 2;
    vimoswcs->nxpix = nxpix;
    vimoswcs->nypix = nypix;

    vimoswcs->vimoswcsproj = vimoswcsproj0;

    vimoswcs->crpix[0] = xrpix;
    vimoswcs->crpix[1] = yrpix;
    vimoswcs->xrefpix = vimoswcs->crpix[0];
    vimoswcs->yrefpix = vimoswcs->crpix[1];
    vimoswcs->lin.crpix = vimoswcs->crpix;

    vimoswcs->crval[0] = cra;
    vimoswcs->crval[1] = cdec;
    vimoswcs->xref = vimoswcs->crval[0];
    vimoswcs->yref = vimoswcs->crval[1];
    vimoswcs->cel.ref[0] = vimoswcs->crval[0];
    vimoswcs->cel.ref[1] = vimoswcs->crval[1];
    vimoswcs->cel.ref[2] = 999.0;

    strcpy (vimoswcs->c1type,"RA");
    strcpy (vimoswcs->c2type,"DEC");

/* Allan Brighton: 28.4.98: for backward compat., remove leading "--" */
    while (proj && *proj == '-')
	proj++;
    strcpy (vimoswcs->ptype,proj);
    strcpy (vimoswcs->ctype[0],"RA---");
    strcpy (vimoswcs->ctype[1],"DEC--");
    strcat (vimoswcs->ctype[0],proj);
    strcat (vimoswcs->ctype[1],proj);

    if (vimoswcstype (vimoswcs, vimoswcs->ctype[0], vimoswcs->ctype[1])) {
	vimoswcsfree (vimoswcs);
	return (NULL);
	}
    
    /* Approximate world coordinate system from a known plate scale */
    cdelt1 = -secpix / 3600.0;
    cdelt2 = secpix / 3600.0;
    vimoswcsdeltset (vimoswcs, cdelt1, cdelt2, rotate);
    vimoswcs->lin.cdelt = vimoswcs->cdelt;
    vimoswcs->lin.pc = vimoswcs->pc;

    /* Coordinate reference frame and equinox */
    vimoswcs->equinox =  (double) equinox;
    if (equinox > 1980)
	strcpy (vimoswcs->radecsys,"FK5");
    else
	strcpy (vimoswcs->radecsys,"FK4");
    if (epoch > 0)
	vimoswcs->epoch = epoch;
    else
	vimoswcs->epoch = 0.0;
    vimoswcs->vimoswcson = 1;

    vimoswcs->sysvimoswcs = vimoswcscsys (vimoswcs->radecsys);
    vimoswcsoutinit (vimoswcs, vimoswcs->radecsys);
    vimoswcsininit (vimoswcs, vimoswcs->radecsys);
    vimoswcs->eqout = 0.0;
    vimoswcs->printsys = 1;
    vimoswcs->tabsys = 0;

    /* Initialize special WCS commands */
    setvimoswcscom (vimoswcs);

    return (vimoswcs);
}


/* Set up a WCS structure from subroutine arguments based on FITS keywords */

struct WorldCoor *
vimoswcskinit (
int	naxis1,		/* Number of pixels along x-axis */
int	naxis2,		/* Number of pixels along y-axis */
char	*ctype1,	/* FITS WCS projection for axis 1 */
char	*ctype2,	/* FITS WCS projection for axis 2 */
double	crpix1, 
double crpix2,	/* Reference pixel coordinates */
double	crval1, 
double crval2,	/* Coordinates at reference pixel in degrees */
double	*cd,		/* Rotation matrix, used if not NULL */
double	cdelt1,
double cdelt2,	/* scale in degrees/pixel, ignored if cd is not NULL */
double	crota,		/* Rotation angle in degrees, ignored if cd is not NULL */
int	equinox, /* Equinox of coordinates, 1950 and 2000 supported */
double	epoch)	/* Epoch of coordinates, used for FK4/FK5 conversion
		 * no effect if 0 */
{
    struct WorldCoor *vimoswcs;

    vimoswcs = (struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

    /* Set WCSLIB flags so that structures will be reinitialized */
    vimoswcs->cel.flag = 0;
    vimoswcs->lin.flag = 0;
    vimoswcs->vimoswcsl.flag = 0;

    /* Image dimensions */
    vimoswcs->naxes = 2;
    vimoswcs->lin.naxis = 2;
    vimoswcs->nxpix = naxis1;
    vimoswcs->nypix = naxis2;

    vimoswcs->vimoswcsproj = vimoswcsproj0;

    vimoswcs->crpix[0] = crpix1;
    vimoswcs->crpix[1] = crpix2;
    vimoswcs->xrefpix = vimoswcs->crpix[0];
    vimoswcs->yrefpix = vimoswcs->crpix[1];
    vimoswcs->lin.crpix = vimoswcs->crpix;

    if (vimoswcstype (vimoswcs, ctype1, ctype2)) {
	vimoswcsfree (vimoswcs);
	return (NULL);
	}
    if (vimoswcs->latbase == 90)
	crval2 = 90.0 - crval2;
    else if (vimoswcs->latbase == -90)
	crval2 = crval2 - 90.0;

    vimoswcs->crval[0] = crval1;
    vimoswcs->crval[1] = crval2;
    vimoswcs->xref = vimoswcs->crval[0];
    vimoswcs->yref = vimoswcs->crval[1];
    vimoswcs->cel.ref[0] = vimoswcs->crval[0];
    vimoswcs->cel.ref[1] = vimoswcs->crval[1];
    vimoswcs->cel.ref[2] = 999.0;

    if (cd != NULL)
	vimoswcscdset (vimoswcs, cd);

    else if (cdelt1 != 0.0)
	vimoswcsdeltset (vimoswcs, cdelt1, cdelt2, crota);

    else {
	vimoswcsdeltset (vimoswcs, 1.0, 1.0, crota);
	setvimoswcserr ("WCSRESET: setting CDELT to 1");
	}
    vimoswcs->lin.cdelt = vimoswcs->cdelt;
    vimoswcs->lin.pc = vimoswcs->pc;

    /* Coordinate reference frame and equinox */
    vimoswcs->equinox =  (double) equinox;
    if (equinox > 1980)
	strcpy (vimoswcs->radecsys,"FK5");
    else
	strcpy (vimoswcs->radecsys,"FK4");
    if (epoch > 0)
	vimoswcs->epoch = epoch;
    else
	vimoswcs->epoch = 0.0;
    vimoswcs->vimoswcson = 1;

    strcpy (vimoswcs->radecout, vimoswcs->radecsys);
    vimoswcs->sysvimoswcs = vimoswcscsys (vimoswcs->radecsys);
    vimoswcsoutinit (vimoswcs, vimoswcs->radecsys);
    vimoswcsininit (vimoswcs, vimoswcs->radecsys);
    vimoswcs->eqout = 0.0;
    vimoswcs->printsys = 1;
    vimoswcs->tabsys = 0;

    /* Initialize special WCS commands */
    setvimoswcscom (vimoswcs);

    return (vimoswcs);
}


/* Set projection in WCS structure from FITS keyword values */

int
vimoswcstype (

struct WorldCoor *vimoswcs,	/* World coordinate system structure */
char	*ctype1,	/* FITS WCS projection for axis 1 */
char	*ctype2)	/* FITS WCS projection for axis 2 */

{
    int i, iproj;
    int nctype = 30;
    char ctypes[30][4];

    strcpy (ctypes[0], "LIN");
    strcpy (ctypes[1], "AZP");
    strcpy (ctypes[2], "TAN");
    strcpy (ctypes[3], "SIN");
    strcpy (ctypes[4], "STG");
    strcpy (ctypes[5], "ARC");
    strcpy (ctypes[6], "ZPN");
    strcpy (ctypes[7], "ZEA");
    strcpy (ctypes[8], "AIR");
    strcpy (ctypes[9], "CYP");
    strcpy (ctypes[10], "CAR");
    strcpy (ctypes[11], "MER");
    strcpy (ctypes[12], "CEA");
    strcpy (ctypes[13], "COP");
    strcpy (ctypes[14], "COD");
    strcpy (ctypes[15], "COE");
    strcpy (ctypes[16], "COO");
    strcpy (ctypes[17], "BON");
    strcpy (ctypes[18], "PCO");
    strcpy (ctypes[19], "GLS");
    strcpy (ctypes[20], "PAR");
    strcpy (ctypes[21], "AIT");
    strcpy (ctypes[22], "MOL");
    strcpy (ctypes[23], "CSC");
    strcpy (ctypes[24], "QSC");
    strcpy (ctypes[25], "TSC");
    strcpy (ctypes[26], "NCP");
    strcpy (ctypes[27], "DSS");
    strcpy (ctypes[28], "PLT");
    strcpy (ctypes[29], "TNX");

    if (!strncmp (ctype1, "LONG",4))
	strncpy (ctype1, "XLON",4);

    strcpy (vimoswcs->ctype[0], ctype1);
    strcpy (vimoswcs->c1type, ctype1);
    strcpy (vimoswcs->ptype, ctype1);

    /* Linear coordinates */
    if (!strncmp (ctype1,"LINEAR",6))
	vimoswcs->prjcode = VIMOSWCS_LIN;

    /* Pixel coordinates */
    else if (!strncmp (ctype1,"PIXEL",6))
	vimoswcs->prjcode = VIMOSWCS_PIX;

    /* Set up right ascension, declination, latitude, or longitude */
    else if (ctype1[0] == 'R' || ctype1[0] == 'D' ||
	     ctype1[0] == 'A' || ctype1[1] == 'L') {
	vimoswcs->c1type[0] = ctype1[0];
	vimoswcs->c1type[1] = ctype1[1];
	if (ctype1[2] == '-') {
	    vimoswcs->c1type[2] = 0;
	    iproj = 3;
	    }
	else {
	    vimoswcs->c1type[2] = ctype1[2];
	    iproj = 4;
	    if (ctype1[3] == '-') {
		vimoswcs->c1type[3] = 0;
		}
	    else {
		vimoswcs->c1type[3] = ctype1[3];
		vimoswcs->c1type[4] = 0;
		}
	    }
	if (ctype1[iproj] == '-') iproj = iproj + 1;
	if (ctype1[iproj] == '-') iproj = iproj + 1;
	if (ctype1[iproj] == '-') iproj = iproj + 1;
	if (ctype1[iproj] == '-') iproj = iproj + 1;
	vimoswcs->ptype[0] = ctype1[iproj];
	vimoswcs->ptype[1] = ctype1[iproj+1];
	vimoswcs->ptype[2] = ctype1[iproj+2];
	vimoswcs->ptype[3] = 0;
	sprintf (vimoswcs->ctype[0],"%-4s%4s",vimoswcs->c1type,vimoswcs->ptype);
	for (i = 0; i < 8; i++)
	    if (vimoswcs->ctype[0][i] == ' ') vimoswcs->ctype[0][i] = '-';

	/*  Find projection type  */
	vimoswcs->prjcode = 0;  /* default type is linear */
	for (i = 1; i < nctype; i++) {
	    if (!strncmp(vimoswcs->ptype, ctypes[i], 3))
		vimoswcs->prjcode = i;
	    }

	/* Handle "obsolete" NCP projection */
	if (vimoswcs->prjcode == VIMOSWCS_NCP) {
	    if (vimoswcs->vimoswcsproj == VIMOSWCS_BEST)
		vimoswcs->vimoswcsproj = VIMOSWCS_OLD;
	    else if (vimoswcs->vimoswcsproj == VIMOSWCS_ALT)
		vimoswcs->vimoswcsproj = VIMOSWCS_NEW;
	    }

	/* Work around bug in WCSLIB handling of CAR projection */
	else if (vimoswcs->prjcode == VIMOSWCS_CAR) {
	    if (vimoswcs->vimoswcsproj == VIMOSWCS_BEST)
		vimoswcs->vimoswcsproj = VIMOSWCS_OLD;
	    else if (vimoswcs->vimoswcsproj == VIMOSWCS_ALT)
		vimoswcs->vimoswcsproj = VIMOSWCS_NEW;
	    }

	/* Work around bug in WCSLIB handling of COE projection */
	else if (vimoswcs->prjcode == VIMOSWCS_COE) {
	    if (vimoswcs->vimoswcsproj == VIMOSWCS_BEST)
		vimoswcs->vimoswcsproj = VIMOSWCS_OLD;
	    else if (vimoswcs->vimoswcsproj == VIMOSWCS_ALT)
		vimoswcs->vimoswcsproj = VIMOSWCS_NEW;
	    }

	else if (vimoswcs->vimoswcsproj == VIMOSWCS_BEST)
	    vimoswcs->vimoswcsproj = VIMOSWCS_NEW;

	else if (vimoswcs->vimoswcsproj == VIMOSWCS_ALT)
	    vimoswcs->vimoswcsproj = VIMOSWCS_OLD;

	if (vimoswcs->vimoswcsproj == VIMOSWCS_OLD && (
	    vimoswcs->prjcode != VIMOSWCS_STG && vimoswcs->prjcode != VIMOSWCS_AIT &&
	    vimoswcs->prjcode != VIMOSWCS_MER && vimoswcs->prjcode != VIMOSWCS_GLS &&
	    vimoswcs->prjcode != VIMOSWCS_ARC && vimoswcs->prjcode != VIMOSWCS_TAN &&
	    vimoswcs->prjcode != VIMOSWCS_TNX && vimoswcs->prjcode != VIMOSWCS_SIN &&
	    vimoswcs->prjcode != VIMOSWCS_PIX && vimoswcs->prjcode != VIMOSWCS_LIN &&
	    vimoswcs->prjcode != VIMOSWCS_CAR && vimoswcs->prjcode != VIMOSWCS_COE &&
	    vimoswcs->prjcode != VIMOSWCS_NCP))
	    vimoswcs->vimoswcsproj = VIMOSWCS_NEW;

	/* Handle NOAO corrected TNX as uncorrected TAN if oldvimoswcs is set */
	if (vimoswcs->vimoswcsproj == VIMOSWCS_OLD && vimoswcs->prjcode == VIMOSWCS_TNX) {
	    vimoswcs->ctype[0][6] = 'A';
	    vimoswcs->ctype[0][7] = 'N';
	    vimoswcs->prjcode = VIMOSWCS_TAN;
	    }
	}

    /* If not sky coordinates, assume linear */
    else {
	vimoswcs->prjcode = VIMOSWCS_LIN;
	return (0);
	}

    /* Second coordinate type */
    if (!strncmp (ctype2, "NPOL",4)) {
	ctype2[0] = ctype1[0];
	strncpy (ctype2+1, "LAT",3);
	vimoswcs->latbase = 90;
	strcpy (vimoswcs->radecsys,"NPOLE");
	vimoswcs->sysvimoswcs = VIMOSWCS_NPOLE;
	}
    else if (!strncmp (ctype2, "SPA-",4)) {
	ctype2[0] = ctype1[0];
	strncpy (ctype2+1, "LAT",3);
	vimoswcs->latbase = -90;
	strcpy (vimoswcs->radecsys,"SPA");
	vimoswcs->sysvimoswcs = VIMOSWCS_SPA;
	}
    else
	vimoswcs->latbase = 0;
    strcpy (vimoswcs->ctype[1], ctype2);
    strcpy (vimoswcs->c2type, ctype2);

    /* Linear coordinates */
    if (!strncmp (ctype2,"LINEAR",6))
	vimoswcs->prjcode = VIMOSWCS_LIN;

    /* Pixel coordinates */
    else if (!strncmp (ctype2,"PIXEL",6))
	vimoswcs->prjcode = VIMOSWCS_PIX;

    /* Set up right ascension, declination, latitude, or longitude */
    else if (ctype2[0] == 'R' || ctype2[0] == 'D' ||
	     ctype2[0] == 'A' || ctype2[1] == 'L') {
	vimoswcs->c2type[0] = ctype2[0];
	vimoswcs->c2type[1] = ctype2[1];
	if (ctype2[2] == '-') {
	    vimoswcs->c2type[2] = 0;
	    iproj = 3;
	    }
	else {
	    vimoswcs->c2type[2] = ctype2[2];
	    iproj = 4;
	    if (ctype2[3] == '-') {
		vimoswcs->c2type[3] = 0;
		}
	    else {
		vimoswcs->c2type[3] = ctype2[3];
		vimoswcs->c2type[4] = 0;
		}
	    }
	if (ctype2[iproj] == '-') iproj = iproj + 1;
	if (ctype2[iproj] == '-') iproj = iproj + 1;
	if (ctype2[iproj] == '-') iproj = iproj + 1;
	if (ctype2[iproj] == '-') iproj = iproj + 1;
	vimoswcs->ptype[0] = ctype2[iproj];
	vimoswcs->ptype[1] = ctype2[iproj+1];
	vimoswcs->ptype[2] = ctype2[iproj+2];
	vimoswcs->ptype[3] = 0;

	if (!strncmp (ctype1, "DEC", 3) ||
	    !strncmp (ctype1+1, "LAT", 3))
	    vimoswcs->coorflip = 1;
	else
	    vimoswcs->coorflip = 0;
	if (ctype2[1] == 'L' || ctype2[0] == 'A') {
	    vimoswcs->degout = 1;
	    vimoswcs->ndec = 5;
	    }
	else {
	    vimoswcs->degout = 0;
	    vimoswcs->ndec = 3;
	    }
	sprintf (vimoswcs->ctype[1],"%-4s%4s",vimoswcs->c2type,vimoswcs->ptype);
	for (i = 0; i < 8; i++)
	    if (vimoswcs->ctype[1][i] == ' ') vimoswcs->ctype[1][i] = '-';
	}

    /* If not sky coordinates, assume linear */
    else {
	vimoswcs->prjcode = VIMOSWCS_LIN;
	}

    return (0);
}


int
vimoswcsreset (

struct WorldCoor *vimoswcs,		/* World coordinate system data structure */
double crpix1,
double crpix2,		/* Reference pixel coordinates */
double crval1,
double crval2,		/* Coordinates at reference pixel in degrees */
double cdelt1, 
double cdelt2,		/* scale in degrees/pixel, ignored if cd is not NULL */
double crota,			/* Rotation angle in degrees, ignored if cd is not NULL */
double *cd)			/* Rotation matrix, used if not NULL */
{

    if (novimoswcs (vimoswcs))
	return (-1);

    /* Set WCSLIB flags so that structures will be reinitialized */
    vimoswcs->cel.flag = 0;
    vimoswcs->lin.flag = 0;
    vimoswcs->vimoswcsl.flag = 0;

    /* Reference pixel coordinates and WCS value */
    vimoswcs->crpix[0] = crpix1;
    vimoswcs->crpix[1] = crpix2;
    vimoswcs->xrefpix = vimoswcs->crpix[0];
    vimoswcs->yrefpix = vimoswcs->crpix[1];
    vimoswcs->lin.crpix = vimoswcs->crpix;

    vimoswcs->crval[0] = crval1;
    vimoswcs->crval[1] = crval2;
    vimoswcs->xref = vimoswcs->crval[0];
    vimoswcs->yref = vimoswcs->crval[1];
    if (vimoswcs->coorflip) {
	vimoswcs->cel.ref[1] = vimoswcs->crval[0];
	vimoswcs->cel.ref[0] = vimoswcs->crval[1];
	}
    else {
	vimoswcs->cel.ref[0] = vimoswcs->crval[0];
	vimoswcs->cel.ref[1] = vimoswcs->crval[1];
	}
    /* Keep ref[2] and ref[3] from input */

    /* Initialize to no plate fit */
    vimoswcs->ncoeff1 = 0;
    vimoswcs->ncoeff2 = 0;

    if (cd != NULL)
	vimoswcscdset (vimoswcs, cd);

    else if (cdelt1 != 0.0)
	vimoswcsdeltset (vimoswcs, cdelt1, cdelt2, crota);

    else {
	vimoswcs->xinc = 1.0;
	vimoswcs->yinc = 1.0;
	setvimoswcserr ("WCSRESET: setting CDELT to 1");
	}

    /* Coordinate reference frame, equinox, and epoch */
    if (!strncmp (vimoswcs->ptype,"LINEAR",6) ||
	!strncmp (vimoswcs->ptype,"PIXEL",5))
	vimoswcs->degout = -1;

    vimoswcs->vimoswcson = 1;
    return (0);
}

void
vimoswcseqset (

struct WorldCoor *vimoswcs,		/* World coordinate system data structure */
double equinox)			/* Desired equinox as fractional year */
{
    extern void fk425e(), fk524e();

    if (novimoswcs (vimoswcs))
	return;

    /* Leave WCS alone if already at desired equinox */
    if (vimoswcs->equinox == equinox)
	return;

    /* Convert center from B1950 (FK4) to J2000 (FK5) */
    if (equinox == 2000.0 && vimoswcs->equinox == 1950.0) {
	if (vimoswcs->coorflip) { 
	    fk425e (&vimoswcs->crval[1], &vimoswcs->crval[0], vimoswcs->epoch);
	    vimoswcs->cel.ref[1] = vimoswcs->crval[0];
	    vimoswcs->cel.ref[0] = vimoswcs->crval[1];
	    }
	else {
	    fk425e (&vimoswcs->crval[0], &vimoswcs->crval[1], vimoswcs->epoch);
	    vimoswcs->cel.ref[0] = vimoswcs->crval[0];
	    vimoswcs->cel.ref[1] = vimoswcs->crval[1];
	    }
	vimoswcs->xref = vimoswcs->crval[0];
	vimoswcs->yref = vimoswcs->crval[1];
	vimoswcs->equinox = 2000.0;
	strcpy (vimoswcs->radecsys, "FK5");
	vimoswcs->sysvimoswcs = VIMOSWCS_J2000;
	vimoswcs->cel.flag = 0;
	vimoswcs->vimoswcsl.flag = 0;
	}

    /* Convert center from J2000 (FK5) to B1950 (FK4) */
    else if (equinox == 1950.0 && vimoswcs->equinox == 2000.0) {
	if (vimoswcs->coorflip) { 
	    fk524e (&vimoswcs->crval[1], &vimoswcs->crval[0], vimoswcs->epoch);
	    vimoswcs->cel.ref[1] = vimoswcs->crval[0];
	    vimoswcs->cel.ref[0] = vimoswcs->crval[1];
	    }
	else {
	    fk524e (&vimoswcs->crval[0], &vimoswcs->crval[1], vimoswcs->epoch);
	    vimoswcs->cel.ref[0] = vimoswcs->crval[0];
	    vimoswcs->cel.ref[1] = vimoswcs->crval[1];
	    }
	vimoswcs->xref = vimoswcs->crval[0];
	vimoswcs->yref = vimoswcs->crval[1];
	vimoswcs->equinox = 1950.0;
	strcpy (vimoswcs->radecsys, "FK4");
	vimoswcs->sysvimoswcs = VIMOSWCS_B1950;
	vimoswcs->cel.flag = 0;
	vimoswcs->vimoswcsl.flag = 0;
	}
    vimoswcsoutinit (vimoswcs, vimoswcs->radecsys);
    vimoswcsininit (vimoswcs, vimoswcs->radecsys);
    return;
}


/* Set scale and rotation in WCS structure */

void
vimoswcscdset (
struct WorldCoor *vimoswcs,	/* World coordinate system structure */
double *cd)			/* CD matrix, ignored if NULL */
{
    extern int vimosmatinv();
    double tcd;

    if (cd == NULL)
	return;

    vimoswcs->rotmat = 1;
    vimoswcs->cd[0] = cd[0];
    vimoswcs->cd[1] = cd[1];
    vimoswcs->cd[2] = cd[2];
    vimoswcs->cd[3] = cd[3];
    (void) vimosmatinv (2, vimoswcs->cd, vimoswcs->dc);

    /* Compute scale */
    vimoswcs->xinc = sqrt (cd[0]*cd[0] + cd[2]*cd[2]);
    vimoswcs->yinc = sqrt (cd[1]*cd[1] + cd[3]*cd[3]);

    /* Deal with x=Dec/y=RA case */
    if (vimoswcs->coorflip) {
	tcd = cd[1];
	cd[1] = -cd[2];
	cd[2] = -tcd;
	}
    vimoswcslibrot (vimoswcs);
    vimoswcs->vimoswcson = 1;

    /* Compute image rotation */
    vimoswcsrotset (vimoswcs);

    vimoswcs->cdelt[0] = vimoswcs->xinc;
    vimoswcs->cdelt[1] = vimoswcs->yinc;

    return;
}


/* Set scale and rotation in WCS structure from axis scale and rotation */

void
vimoswcsdeltset (

struct WorldCoor *vimoswcs,	/* World coordinate system structure */
double cdelt1,		/* degrees/pixel in first axis (or both axes) */
double cdelt2,		/* degrees/pixel in second axis if nonzero */
double crota)		/* Rotation counterclockwise in degrees */
{
    extern int vimosmatinv();
    double *pci;
    double crot, srot;
    int i, j;

    vimoswcs->cdelt[0] = cdelt1;
    if (cdelt2 != 0.0)
	vimoswcs->cdelt[1] = cdelt2;
    else
	vimoswcs->cdelt[1] = cdelt1;
    vimoswcs->xinc = vimoswcs->cdelt[0];
    vimoswcs->yinc = vimoswcs->cdelt[1];
    pci = vimoswcs->pc;
    for (i = 0; i < vimoswcs->lin.naxis; i++) {
	for (j = 0; j < vimoswcs->lin.naxis; j++) {
	    if (i ==j)
		*pci = 1.0;
	    else
		*pci = 0.0;
	    pci++;
	    }
	}
    vimoswcs->rotmat = 0;

    /* If image is reversed, value of CROTA is flipped, too */
    vimoswcs->rot = crota;
    crot = cos (degrad(vimoswcs->rot));
    if (cdelt1 * cdelt2 > 0)
	srot = sin (-degrad(vimoswcs->rot));
    else
	srot = sin (degrad(vimoswcs->rot));

    /* Set CD matrix */
    vimoswcs->cd[0] = vimoswcs->cdelt[0] * crot;
    if (vimoswcs->cdelt[0] < 0)
	vimoswcs->cd[1] = -fabs (vimoswcs->cdelt[1]) * srot;
    else
	vimoswcs->cd[1] = fabs (vimoswcs->cdelt[1]) * srot;
    if (vimoswcs->cdelt[1] < 0)
	vimoswcs->cd[2] = fabs (vimoswcs->cdelt[0]) * srot;
    else
	vimoswcs->cd[2] = -fabs (vimoswcs->cdelt[0]) * srot;
    vimoswcs->cd[3] = vimoswcs->cdelt[1] * crot;
    (void) vimosmatinv (2, vimoswcs->cd, vimoswcs->dc);

    /* Set rotation matrix */
    vimoswcslibrot (vimoswcs);

    /* Set image rotation and mirroring */
    if (vimoswcs->coorflip) {
	if (vimoswcs->cdelt[0] < 0 && vimoswcs->cdelt[1] > 0) {
	    vimoswcs->imflip = 1;
	    vimoswcs->imrot = vimoswcs->rot - 90.0;
	    if (vimoswcs->imrot < -180.0) vimoswcs->imrot = vimoswcs->imrot + 360.0;
	    vimoswcs->pa_north = vimoswcs->rot;
	    vimoswcs->pa_east = vimoswcs->rot - 90.0;
	    if (vimoswcs->pa_east < -180.0) vimoswcs->pa_east = vimoswcs->pa_east + 360.0;
	    }
	else if (vimoswcs->cdelt[0] > 0 && vimoswcs->cdelt[1] < 0) {
	    vimoswcs->imflip = 1;
	    vimoswcs->imrot = vimoswcs->rot + 90.0;
	    if (vimoswcs->imrot > 180.0) vimoswcs->imrot = vimoswcs->imrot - 360.0;
	    vimoswcs->pa_north = vimoswcs->rot;
	    vimoswcs->pa_east = vimoswcs->rot - 90.0;
	    if (vimoswcs->pa_east < -180.0) vimoswcs->pa_east = vimoswcs->pa_east + 360.0;
	    }
	else if (vimoswcs->cdelt[0] > 0 && vimoswcs->cdelt[1] > 0) {
	    vimoswcs->imflip = 0;
	    vimoswcs->imrot = vimoswcs->rot + 90.0;
	    if (vimoswcs->imrot > 180.0) vimoswcs->imrot = vimoswcs->imrot - 360.0;
	    vimoswcs->pa_north = vimoswcs->imrot;
	    vimoswcs->pa_east = vimoswcs->rot + 90.0;
	    if (vimoswcs->pa_east > 180.0) vimoswcs->pa_east = vimoswcs->pa_east - 360.0;
	    }
	else if (vimoswcs->cdelt[0] < 0 && vimoswcs->cdelt[1] < 0) {
	    vimoswcs->imflip = 0;
	    vimoswcs->imrot = vimoswcs->rot - 90.0;
	    if (vimoswcs->imrot < -180.0) vimoswcs->imrot = vimoswcs->imrot + 360.0;
	    vimoswcs->pa_north = vimoswcs->imrot;
	    vimoswcs->pa_east = vimoswcs->rot + 90.0;
	    if (vimoswcs->pa_east > 180.0) vimoswcs->pa_east = vimoswcs->pa_east - 360.0;
	    }
	}
    else {
	if (vimoswcs->cdelt[0] < 0 && vimoswcs->cdelt[1] > 0) {
	    vimoswcs->imflip = 0;
	    vimoswcs->imrot = vimoswcs->rot;
	    vimoswcs->pa_north = vimoswcs->rot + 90.0;
	    if (vimoswcs->pa_north > 180.0) vimoswcs->pa_north = vimoswcs->pa_north - 360.0;
	    vimoswcs->pa_east = vimoswcs->rot + 180.0;
	    if (vimoswcs->pa_east > 180.0) vimoswcs->pa_east = vimoswcs->pa_east - 360.0;
	    }
	else if (vimoswcs->cdelt[0] > 0 && vimoswcs->cdelt[1] < 0) {
	    vimoswcs->imflip = 0;
	    vimoswcs->imrot = vimoswcs->rot + 180.0;
	    if (vimoswcs->imrot > 180.0) vimoswcs->imrot = vimoswcs->imrot - 360.0;
	    vimoswcs->pa_north = vimoswcs->imrot + 90.0;
	    if (vimoswcs->pa_north > 180.0) vimoswcs->pa_north = vimoswcs->pa_north - 360.0;
	    vimoswcs->pa_east = vimoswcs->imrot + 180.0;
	    if (vimoswcs->pa_east > 180.0) vimoswcs->pa_east = vimoswcs->pa_east - 360.0;
	    }
	else if (vimoswcs->cdelt[0] > 0 && vimoswcs->cdelt[1] > 0) {
	    vimoswcs->imflip = 1;
	    vimoswcs->imrot = vimoswcs->rot;
	    vimoswcs->pa_north = vimoswcs->imrot + 90.0;
	    if (vimoswcs->pa_north > 180.0) vimoswcs->pa_north = vimoswcs->pa_north - 360.0;
	    vimoswcs->pa_east = vimoswcs->rot;
	    }
	else if (vimoswcs->cdelt[0] < 0 && vimoswcs->cdelt[1] < 0) {
	    vimoswcs->imflip = 1;
	    vimoswcs->imrot = vimoswcs->rot + 180.0;
	    if (vimoswcs->imrot > 180.0) vimoswcs->imrot = vimoswcs->imrot - 360.0;
	    vimoswcs->pa_north = vimoswcs->imrot + 90.0;
	    if (vimoswcs->pa_north > 180.0) vimoswcs->pa_north = vimoswcs->pa_north - 360.0;
	    vimoswcs->pa_east = vimoswcs->rot + 90.0;
	    if (vimoswcs->pa_east > 180.0) vimoswcs->pa_east = vimoswcs->pa_east - 360.0;
	    }
	}

    return;
}


/* Set scale and rotation in WCS structure */

void
vimoswcspcset (

struct WorldCoor *vimoswcs,	/* World coordinate system structure */
double cdelt1,		/* degrees/pixel in first axis (or both axes) */
double cdelt2,		/* degrees/pixel in second axis if nonzero */
double *pc)		/* Rotation matrix, ignored if NULL */
{
    extern int vimosmatinv();
    double *pci, *pc0i;
    int i, j, naxes;

    if (pc == NULL)
	return;

    naxes = vimoswcs->naxes;
    vimoswcs->cdelt[0] = cdelt1;
    if (cdelt2 != 0.0)
	vimoswcs->cdelt[1] = cdelt2;
    else
	vimoswcs->cdelt[1] = cdelt1;
    vimoswcs->xinc = vimoswcs->cdelt[0];
    vimoswcs->yinc = vimoswcs->cdelt[1];

    /* Set rotation matrix */
    pci = vimoswcs->pc;
    pc0i = pc;
    for (i = 0; i < naxes; i++) {
	for (j = 0; j < naxes; j++) {
	    *pci = *pc0i;
	    pci++;
	    pc0i++;
	    }
	}

    /* Set CD matrix */
    if (naxes < 3) {
	vimoswcs->cd[0] = pc[0] * vimoswcs->cdelt[0];
	vimoswcs->cd[1] = pc[1] * vimoswcs->cdelt[1];
	vimoswcs->cd[2] = pc[2] * vimoswcs->cdelt[0];
	vimoswcs->cd[3] = pc[3] * vimoswcs->cdelt[1];
	}
    else if (naxes == 3) {
	vimoswcs->cd[0] = pc[0] * vimoswcs->cdelt[0];
	vimoswcs->cd[1] = pc[1] * vimoswcs->cdelt[1];
	vimoswcs->cd[2] = pc[3] * vimoswcs->cdelt[0];
	vimoswcs->cd[3] = pc[4] * vimoswcs->cdelt[1];
	}
    else if (naxes == 4) {
	vimoswcs->cd[0] = pc[0] * vimoswcs->cdelt[0];
	vimoswcs->cd[1] = pc[1] * vimoswcs->cdelt[1];
	vimoswcs->cd[2] = pc[4] * vimoswcs->cdelt[0];
	vimoswcs->cd[3] = pc[5] * vimoswcs->cdelt[1];
	}
    (void) vimosmatinv (naxes, vimoswcs->cd, vimoswcs->dc);
    vimoswcs->rotmat = 1;

    (void)vimoslinset (&vimoswcs->lin);
    vimoswcs->vimoswcson = 1;

    vimoswcsrotset (vimoswcs);

    return;
}


/* Set up rotation matrix for WCSLIB projection subroutines */

static void
vimoswcslibrot (

struct WorldCoor *vimoswcs)	/* World coordinate system structure */

{
    int i, mem, naxes;

    naxes = vimoswcs->naxes;
    mem = naxes * naxes * sizeof(double);
    if (vimoswcs->lin.piximg == NULL)
	vimoswcs->lin.piximg = (double*)malloc(mem);
    if (vimoswcs->lin.piximg != NULL) {
	if (vimoswcs->lin.imgpix == NULL)
	    vimoswcs->lin.imgpix = (double*)malloc(mem);
	if (vimoswcs->lin.imgpix != NULL) {
	    vimoswcs->lin.flag = LINSET;
	    if (naxes == 2) {
		for (i = 0; i < 4; i++) {
		    vimoswcs->lin.piximg[i] = vimoswcs->cd[i];
		    }
		}
	    else if (naxes == 3) {
		for (i = 0; i < 9; i++)
		    vimoswcs->lin.piximg[i] = 0.0;
		vimoswcs->lin.piximg[0] = vimoswcs->cd[0];
		vimoswcs->lin.piximg[1] = vimoswcs->cd[1];
		vimoswcs->lin.piximg[3] = vimoswcs->cd[2];
		vimoswcs->lin.piximg[4] = vimoswcs->cd[3];
		vimoswcs->lin.piximg[8] = 1.0;
		}
	    else if (naxes == 4) {
		for (i = 0; i < 16; i++)
		    vimoswcs->lin.piximg[i] = 0.0;
		vimoswcs->lin.piximg[0] = vimoswcs->cd[0];
		vimoswcs->lin.piximg[1] = vimoswcs->cd[1];
		vimoswcs->lin.piximg[4] = vimoswcs->cd[2];
		vimoswcs->lin.piximg[5] = vimoswcs->cd[3];
		vimoswcs->lin.piximg[10] = 1.0;
		vimoswcs->lin.piximg[15] = 1.0;
		}
	    (void) vimosmatinv (naxes, vimoswcs->lin.piximg, vimoswcs->lin.imgpix);
	    vimoswcs->lin.crpix = vimoswcs->crpix;
	    vimoswcs->lin.cdelt = vimoswcs->cdelt;
	    vimoswcs->lin.pc = vimoswcs->pc;
	    vimoswcs->lin.flag = LINSET;
	    }
	}
    return;
}


/* Compute image rotation */

void
vimoswcsrotset (

struct WorldCoor *vimoswcs)	/* World coordinate system structure */
{
    int off;
    double cra, cdec, xc, xn, xe, yc, ynn, ye;

    /* If image is one-dimensional, leave rotation angle alone */
    if (vimoswcs->nxpix < 1.5 || vimoswcs->nypix < 1.5) {
	vimoswcs->imrot = vimoswcs->rot;
	vimoswcs->pa_north = vimoswcs->rot + 90.0;
	vimoswcs->pa_east = vimoswcs->rot + 180.0;
	return;
	}


    /* Do not try anything if image is LINEAR (not Cartesian projection) */
    if (vimoswcs->sysvimoswcs == VIMOSWCS_LINEAR)
	return;

    vimoswcs->xinc = fabs (vimoswcs->xinc);
    vimoswcs->yinc = fabs (vimoswcs->yinc);

    /* Compute position angles of North and East in image */
    xc = vimoswcs->xrefpix;
    yc = vimoswcs->yrefpix;
    pix2vimoswcs (vimoswcs, xc, yc, &cra, &cdec);
    if (vimoswcs->coorflip) {
	vimoswcs2pix (vimoswcs, cra+vimoswcs->yinc, cdec, &xe, &ye, &off);
	vimoswcs2pix (vimoswcs, cra, cdec+vimoswcs->xinc, &xn, &ynn, &off);
	}
    else {
	vimoswcs2pix (vimoswcs, cra+vimoswcs->xinc, cdec, &xe, &ye, &off);
	vimoswcs2pix (vimoswcs, cra, cdec+vimoswcs->yinc, &xn, &ynn, &off);
	}
    vimoswcs->pa_north = raddeg (atan2 (ynn-yc, xn-xc));
    if (vimoswcs->pa_north < -90.0)
	vimoswcs->pa_north = vimoswcs->pa_north + 360.0;
    vimoswcs->pa_east = raddeg (atan2 (ye-yc, xe-xc));
    if (vimoswcs->pa_east < -90.0)
	vimoswcs->pa_east = vimoswcs->pa_east + 360.0;

    /* Compute image rotation angle from North */
    if (vimoswcs->pa_north < -90.0)
	vimoswcs->imrot = 270.0 + vimoswcs->pa_north;
    else
	vimoswcs->imrot = vimoswcs->pa_north - 90.0;

    /* Compute CROTA */
    if (vimoswcs->coorflip) {
	vimoswcs->rot = vimoswcs->imrot + 90.0;
	if (vimoswcs->rot > 180)
	    vimoswcs->rot = vimoswcs->rot - 360.0;
	}
    else
	vimoswcs->rot = vimoswcs->imrot;

    /* Set image mirror flag based on axis orientation */
    vimoswcs->imflip = 0;
    if (vimoswcs->pa_east - vimoswcs->pa_north < -80.0 &&
	vimoswcs->pa_east - vimoswcs->pa_north > -100.0)
	vimoswcs->imflip = 1;
    if (vimoswcs->pa_east - vimoswcs->pa_north < 280.0 &&
	vimoswcs->pa_east - vimoswcs->pa_north > 260.0)
	vimoswcs->imflip = 1;
    if (vimoswcs->pa_north - vimoswcs->pa_east > 80.0 &&
	vimoswcs->pa_north - vimoswcs->pa_east < 100.0)
	vimoswcs->imflip = 1;
    if (vimoswcs->coorflip) {
	if (vimoswcs->imflip)
	    vimoswcs->yinc = -vimoswcs->yinc;
	}
    else {
	if (!vimoswcs->imflip)
	    vimoswcs->xinc = -vimoswcs->xinc;
	}

    return;
}


/* Return 1 if WCS structure is filled, else 0 */

int
isvimoswcs (

struct WorldCoor *vimoswcs)		/* World coordinate system structure */

{
    if (vimoswcs == NULL)
	return (0);
    else
	return (vimoswcs->vimoswcson);
}


/* Return 0 if WCS structure is filled, else 1 */

int
novimoswcs (

struct WorldCoor *vimoswcs)		/* World coordinate system structure */

{
    if (vimoswcs == NULL)
	return (1);
    else
	return (!vimoswcs->vimoswcson);
}


/* Reset the center of a WCS structure */

void
vimoswcsshift (

struct WorldCoor *vimoswcs,	/* World coordinate system structure */
double	rra,		/* Reference pixel right ascension in degrees */
double	rdec,		/* Reference pixel declination in degrees */
char	*coorsys)	/* FK4 or FK5 coordinates (1950 or 2000) */

{
    if (novimoswcs (vimoswcs))
	return;

/* Approximate world coordinate system from a known plate scale */
    vimoswcs->crval[0] = rra;
    vimoswcs->crval[1] = rdec;
    vimoswcs->xref = vimoswcs->crval[0];
    vimoswcs->yref = vimoswcs->crval[1];


/* Coordinate reference frame */
    strcpy (vimoswcs->radecsys,coorsys);
    vimoswcs->sysvimoswcs = vimoswcscsys (coorsys);
    if (vimoswcs->sysvimoswcs == VIMOSWCS_B1950)
	vimoswcs->equinox = 1950.0;
    else
	vimoswcs->equinox = 2000.0;

    return;
}

/* Print position of WCS center, if WCS is set */

void
vimoswcscent (

struct WorldCoor *vimoswcs)		/* World coordinate system structure */

{
    double	xpix,ypix, xpos1, xpos2, ypos1, ypos2;
    char vimoswcstring[32];
    double width, height, secpix, secpixh, secpixw;
    int lstr = 32;

    if (novimoswcs (vimoswcs))
	(void)fprintf (stderr,"No WCS information available\n");
    else {
	if (vimoswcs->prjcode == VIMOSWCS_DSS)
	    (void)fprintf (stderr,"WCS plate center  %s\n", vimoswcs->center);
	xpix = 0.5 * vimoswcs->nxpix;
	ypix = 0.5 * vimoswcs->nypix;
	(void) pix2vimoswcst (vimoswcs,xpix,ypix,vimoswcstring, lstr);
	(void)fprintf (stderr,"WCS center %s %s %s %s at pixel (%.2f,%.2f)\n",
		     vimoswcs->ctype[0],vimoswcs->ctype[1],vimoswcstring,vimoswcs->ptype,xpix,ypix);

	/* Image width */
	(void) pix2vimoswcs (vimoswcs,1.0,ypix,&xpos1,&ypos1);
	(void) pix2vimoswcs (vimoswcs,vimoswcs->nxpix,ypix,&xpos2,&ypos2);
	if (vimoswcs->sysvimoswcs == VIMOSWCS_LINEAR) {
	    width = xpos2 - xpos1;
	    if (width < 100.0)
	    (void)fprintf (stderr, "WCS width = %.5f %s ",width, vimoswcs->units[0]);
	    else
	    (void)fprintf (stderr, "WCS width = %.3f %s ",width, vimoswcs->units[0]);
	    }
	else {
	    width = vimoswcsdist (xpos1,ypos1,xpos2,ypos2);
	    if (width < 1/60.0)
		(void)fprintf (stderr, "WCS width = %.2f arcsec ",width*3600.0);
	    else if (width < 1.0)
		(void)fprintf (stderr, "WCS width = %.2f arcmin ",width*60.0);
	    else
		(void)fprintf (stderr, "WCS width = %.3f degrees ",width);
	    }
	secpixw = width / (vimoswcs->nxpix - 1.0);

	/* Image height */
	(void) pix2vimoswcs (vimoswcs,xpix,1.0,&xpos1,&ypos1);
	(void) pix2vimoswcs (vimoswcs,xpix,vimoswcs->nypix,&xpos2,&ypos2);
	if (vimoswcs->sysvimoswcs == VIMOSWCS_LINEAR) {
	    height = ypos2 - ypos1;
	    if (height < 100.0)
	    (void)fprintf (stderr, " height = %.5f %s ",height, vimoswcs->units[1]);
	    else
	    (void)fprintf (stderr, " height = %.3f %s ",height, vimoswcs->units[1]);
	    }
	else {
	    height = vimoswcsdist (xpos1,ypos1,xpos2,ypos2);
	    if (height < 1/60.0)
		(void) fprintf (stderr, " height = %.2f arcsec",height*3600.0);
	    else if (height < 1.0)
		(void) fprintf (stderr, " height = %.2f arcmin",height*60.0);
	    else
		(void) fprintf (stderr, " height = %.3f degrees",height);
	    }
	secpixh = height / (vimoswcs->nypix - 1.0);

	/* Image scale */
	if (vimoswcs->sysvimoswcs == VIMOSWCS_LINEAR) {
	    (void) fprintf (stderr,"\n");
	    (void) fprintf (stderr,"WCS  %.5f %s/pixel, %.5f %s/pixel\n",
			    vimoswcs->xinc,vimoswcs->units[0],vimoswcs->yinc,vimoswcs->units[1]);
	    }
	else {
	    if (vimoswcs->xinc != 0.0 && vimoswcs->yinc != 0.0)
		secpix = (fabs(vimoswcs->xinc) + fabs(vimoswcs->yinc)) * 0.5 * 3600.0;
	    else if (secpixh > 0.0 && secpixw > 0.0)
		secpix = (secpixw + secpixh) * 0.5 * 3600.0;
	    else if (vimoswcs->xinc != 0.0 || vimoswcs->yinc != 0.0)
		secpix = (fabs(vimoswcs->xinc) + fabs(vimoswcs->yinc)) * 3600.0;
	    else
		secpix = (secpixw + secpixh) * 3600.0;
	    if (secpix < 100.0)
		(void) fprintf (stderr, "  %.3f arcsec/pixel\n",secpix);
	    else if (secpix < 3600.0)
		(void) fprintf (stderr, "  %.3f arcmin/pixel\n",secpix/60.0);
	    else
		(void) fprintf (stderr, "  %.3f degrees/pixel\n",secpix/3600.0);
	    }
	}
    return;
}

/* Return RA and Dec of image center, plus size in RA and Dec */

void
vimoswcssize (

struct WorldCoor *vimoswcs,	/* World coordinate system structure */
double	*cra,		/* Right ascension of image center (deg) (returned) */
double	*cdec,		/* Declination of image center (deg) (returned) */
double	*dra,		/* Half-width in right ascension (deg) (returned) */
double	*ddec)		/* Half-width in declination (deg) (returned) */

{
    double width, height;

    /* Find right ascension and declination of coordinates */
    if (isvimoswcs(vimoswcs)) {
	vimoswcsfull (vimoswcs, cra, cdec, &width, &height);
	*dra = 0.5 * width / cos (degrad (*cdec));
	*ddec = 0.5 * height;
	}
    else {
	*cra = 0.0;
	*cdec = 0.0;
	*dra = 0.0;
	*ddec = 0.0;
	}
    return;
}


/* Return RA and Dec of image center, plus size in degrees */

void
vimoswcsfull (vimoswcs, cra, cdec, width, height)

struct WorldCoor *vimoswcs;	/* World coordinate system structure */
double	*cra;		/* Right ascension of image center (deg) (returned) */
double	*cdec;		/* Declination of image center (deg) (returned) */
double	*width;		/* Width in degrees (returned) */
double	*height;	/* Height in degrees (returned) */

{
    double xpix, ypix, xpos1, xpos2, ypos1, ypos2;
    double xcent, ycent;

    /* Find right ascension and declination of coordinates */
    if (isvimoswcs(vimoswcs)) {
	xpix = 0.5 * vimoswcs->nxpix;
	ypix = 0.5 * vimoswcs->nypix;
	(void) pix2vimoswcs (vimoswcs,xpix,ypix,&xcent, &ycent);
	*cra = xcent;
	*cdec = ycent;

	/* Compute image width in degrees */
	(void) pix2vimoswcs (vimoswcs,1.0,ypix,&xpos1,&ypos1);
	(void) pix2vimoswcs (vimoswcs,vimoswcs->nxpix,ypix,&xpos2,&ypos2);
	if (strncmp (vimoswcs->ptype,"LINEAR",6) &&
	    strncmp (vimoswcs->ptype,"PIXEL",5))
	    *width = vimoswcsdist (xpos1,ypos1,xpos2,ypos2);
	else
	    *width = sqrt (((ypos2-ypos1) * (ypos2-ypos1)) +
		     ((xpos2-xpos1) * (xpos2-xpos1)));

	/* Compute image height in degrees */
	(void) pix2vimoswcs (vimoswcs,xpix,1.0,&xpos1,&ypos1);
	(void) pix2vimoswcs (vimoswcs,xpix,vimoswcs->nypix,&xpos2,&ypos2);
	if (strncmp (vimoswcs->ptype,"LINEAR",6) &&
	    strncmp (vimoswcs->ptype,"PIXEL",5))
	    *height = vimoswcsdist (xpos1,ypos1,xpos2,ypos2);
	else
	    *height = sqrt (((ypos2-ypos1) * (ypos2-ypos1)) +
		      ((xpos2-xpos1) * (xpos2-xpos1)));
	}

    else {
	*cra = 0.0;
	*cdec = 0.0;
	*width = 0.0;
	*height = 0.0;
	}

    return;
}


/* Return minimum and maximum RA and Dec of image in degrees */

void
vimoswcsrange (vimoswcs, ra1, ra2, dec1, dec2)

struct WorldCoor *vimoswcs;	/* World coordinate system structure */
double	*ra1;		/* Minimum right ascension of image (deg) (returned) */
double	*ra2;		/* Maximum right ascension of image (deg) (returned) */
double	*dec1;		/* Minimum declination of image (deg) (returned) */
double	*dec2;		/* Maximum declination of image (deg) (returned) */

{
    double xpos1, xpos2, xpos3, xpos4, ypos1, ypos2, ypos3, ypos4;

    if (isvimoswcs(vimoswcs)) {

	/* Compute image corner coordinates in degrees */
	(void) pix2vimoswcs (vimoswcs,1.0,1.0,&xpos1,&ypos1);
	(void) pix2vimoswcs (vimoswcs,1.0,vimoswcs->nypix,&xpos2,&ypos2);
	(void) pix2vimoswcs (vimoswcs,vimoswcs->nxpix,1.0,&xpos3,&ypos3);
	(void) pix2vimoswcs (vimoswcs,vimoswcs->nxpix,vimoswcs->nypix,&xpos4,&ypos4);

	/* Find minimum right ascension or longitude */
	*ra1 = xpos1;
	if (xpos2 < *ra1) *ra1 = xpos2;
	if (xpos3 < *ra1) *ra1 = xpos3;
	if (xpos4 < *ra1) *ra1 = xpos4;

	/* Find maximum right ascension or longitude */
	*ra2 = xpos1;
	if (xpos2 > *ra2) *ra2 = xpos2;
	if (xpos3 > *ra2) *ra2 = xpos3;
	if (xpos4 > *ra2) *ra2 = xpos4;

	/* Find minimum declination or latitude */
	*dec1 = ypos1;
	if (ypos2 < *dec1) *dec1 = ypos2;
	if (ypos3 < *dec1) *dec1 = ypos3;
	if (ypos4 < *dec1) *dec1 = ypos4;

	/* Find maximum declination or latitude */
	*dec2 = ypos1;
	if (ypos2 > *dec2) *dec2 = ypos2;
	if (ypos3 > *dec2) *dec2 = ypos3;
	if (ypos4 > *dec2) *dec2 = ypos4;
	}

    else {
	*ra1 = 0.0;
	*ra2 = 0.0;
	*dec1 = 0.0;
	*dec2 = 0.0;
	}

    return;
}


/* Compute distance in degrees between two sky coordinates */

double
vimoswcsdist (x1,y1,x2,y2)

double	x1,y1;	/* (RA,Dec) or (Long,Lat) in degrees */
double	x2,y2;	/* (RA,Dec) or (Long,Lat) in degrees */

{
	double xr1, xr2, yr1, yr2;
	double pos1[3], pos2[3], w, diff, cosb;
	int i;

	/* Convert two vectors to direction cosines */
	xr1 = degrad (x1);
	yr1 = degrad (y1);
	cosb = cos (yr1);
	pos1[0] = cos (xr1) * cosb;
	pos1[1] = sin (xr1) * cosb;
	pos1[2] = sin (yr1);

	xr2 = degrad (x2);
	yr2 = degrad (y2);
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
	diff = 2.0 * atan2 (sqrt (w), sqrt (1.0 - w));
	diff = raddeg (diff);
	return (diff);
}


/* Compute distance in degrees between two sky coordinates  away from pole */

double
vimoswcsdiff (x1,y1,x2,y2)

double	x1,y1;	/* (RA,Dec) or (Long,Lat) in degrees */
double	x2,y2;	/* (RA,Dec) or (Long,Lat) in degrees */

{
    double xdiff, ydiff, ycos, diff;

    ycos = cos (degrad ((y2 + y1) / 2.0));
    xdiff = x2 - x1;
    if (xdiff > 180.0)
	xdiff = xdiff - 360.0;
    if (xdiff < -180.0)
	xdiff = xdiff + 360.0;
    xdiff = xdiff / ycos;
    ydiff = (y2 - y1);
    diff = sqrt ((xdiff * xdiff) + (ydiff * ydiff));
    return (diff);
}


/* Initialize catalog search command set by -vimoswcscom */

void
vimoswcscominit (vimoswcs, i, command)

struct WorldCoor *vimoswcs;		/* World coordinate system structure */
int	i;			/* Number of command (0-9) to initialize */
char	*command;		/* command with %s where coordinates will go */

{
    int lcom,icom;

    if (isvimoswcs(vimoswcs)) {
	lcom = strlen (command);
	if (lcom > 0) {
	    if (vimoswcs->command_format[i] != NULL)
		free (vimoswcs->command_format[i]);
	    vimoswcs->command_format[i] = (char *) calloc (lcom+2, 1);
	    if (vimoswcs->command_format[i] == NULL)
		return;
	    for (icom = 0; icom < lcom; icom++) {
		if (command[icom] == '_')
		    vimoswcs->command_format[i][icom] = ' ';
		else
		    vimoswcs->command_format[i][icom] = command[icom];
		}
	    vimoswcs->command_format[i][lcom] = 0;
	    }
	}
    return;
}


/* Execute Unix command with world coordinates (from x,y) and/or filename */

void
vimoswcscom ( vimoswcs, i, filename, xfile, yfile, vimoswcstring )

struct WorldCoor *vimoswcs;		/* World coordinate system structure */
int	i;			/* Number of command (0-9) to execute */
char	*filename;		/* Image file name */
double	xfile,yfile;		/* Image pixel coordinates for WCS command */
char	*vimoswcstring;		/* WCS String from pix2vimoswcst() */
{
    char command[120];
    char comform[120];
    char xystring[32];
    char *fileform, *posform, *imform;
    int ier;

    if (novimoswcs (vimoswcs)) {
	(void)fprintf(stderr,"WCSCOM: no WCS\n");
	return;
	}

    if (vimoswcs->command_format[i] != NULL)
	strcpy (comform, vimoswcs->command_format[i]);
    else
	strcpy (comform, "sgsc -ah %s");

    if (comform[0] > 0) {

	/* Create and execute search command */
	fileform = strsrch (comform,"%f");
	imform = strsrch (comform,"%x");
	posform = strsrch (comform,"%s");
	if (imform != NULL) {
	    *(imform+1) = 's';
	    (void)sprintf (xystring, "%.2f %.2f", xfile, yfile);
	    if (fileform != NULL) {
		*(fileform+1) = 's';
		if (posform == NULL) {
		    if (imform < fileform)
			(void)sprintf(command, comform, xystring, filename);
		    else
			(void)sprintf(command, comform, filename, xystring);
		    }
		else if (fileform < posform) {
		    if (imform < fileform)
			(void)sprintf(command, comform, xystring, filename,
				      vimoswcstring);
		    else if (imform < posform)
			(void)sprintf(command, comform, filename, xystring,
				      vimoswcstring);
		    else
			(void)sprintf(command, comform, filename, vimoswcstring,
				      xystring);
		    }
		else
		    if (imform < posform)
			(void)sprintf(command, comform, xystring, vimoswcstring,
				      filename);
		    else if (imform < fileform)
			(void)sprintf(command, comform, vimoswcstring, xystring,
				      filename);
		    else
			(void)sprintf(command, comform, vimoswcstring, filename,
				      xystring);
		}
	    else if (posform == NULL)
		(void)sprintf(command, comform, xystring);
	    else if (imform < posform)
		(void)sprintf(command, comform, xystring, vimoswcstring);
	    else
		(void)sprintf(command, comform, vimoswcstring, xystring);
	    }
	else if (fileform != NULL) {
	    *(fileform+1) = 's';
	    if (posform == NULL)
		(void)sprintf(command, comform, filename);
	    else if (fileform < posform)
		(void)sprintf(command, comform, filename, vimoswcstring);
	    else
		(void)sprintf(command, comform, vimoswcstring, filename);
	    }
	else
	    (void)sprintf(command, comform, vimoswcstring);
	ier = system (command);
	if (ier)
	    (void)fprintf(stderr,"WCSCOM: %s failed %d\n",command,ier);
	}
    return;
}

/* Initialize WCS output coordinate system for use by PIX2VIMOSWCS() */

void
vimoswcsoutinit (vimoswcs, coorsys)

struct WorldCoor *vimoswcs;	/* World coordinate system structure */
char	*coorsys;	/* Input world coordinate system:
			   FK4, FK5, B1950, J2000, GALACTIC, ECLIPTIC
			   fk4, fk5, b1950, j2000, galactic, ecliptic */
{
    int sysout, i;

    if (novimoswcs (vimoswcs))
	return;

    /* If argument is null, set to image system and equinox */
    if (coorsys == NULL || strlen (coorsys) < 1 ||
	!strcmp(coorsys,"IMSYS") || !strcmp(coorsys,"imsys")) {
	sysout = vimoswcs->sysvimoswcs;
	strcpy (vimoswcs->radecout, vimoswcs->radecsys);
	vimoswcs->eqout = vimoswcs->equinox;
	if (sysout == VIMOSWCS_B1950) {
	    if (vimoswcs->eqout != 1950.0) {
		vimoswcs->radecout[0] = 'B';
		sprintf (vimoswcs->radecout+1,"%.4f", vimoswcs->equinox);
		i = strlen(vimoswcs->radecout) - 1;
		if (vimoswcs->radecout[i] == '0')
		    vimoswcs->radecout[i] = (char)0;
		i = strlen(vimoswcs->radecout) - 1;
		if (vimoswcs->radecout[i] == '0')
		    vimoswcs->radecout[i] = (char)0;
		i = strlen(vimoswcs->radecout) - 1;
		if (vimoswcs->radecout[i] == '0')
		    vimoswcs->radecout[i] = (char)0;
		}
	    else
		strcpy (vimoswcs->radecout, "B1950");
	    }
	else if (sysout == VIMOSWCS_J2000) {
	    if (vimoswcs->eqout != 2000.0) {
		vimoswcs->radecout[0] = 'J';
		sprintf (vimoswcs->radecout+1,"%.4f", vimoswcs->equinox);
		i = strlen(vimoswcs->radecout) - 1;
		if (vimoswcs->radecout[i] == '0')
		    vimoswcs->radecout[i] = (char)0;
		i = strlen(vimoswcs->radecout) - 1;
		if (vimoswcs->radecout[i] == '0')
		    vimoswcs->radecout[i] = (char)0;
		i = strlen(vimoswcs->radecout) - 1;
		if (vimoswcs->radecout[i] == '0')
		    vimoswcs->radecout[i] = (char)0;
		}
	    else
		strcpy (vimoswcs->radecout, "J2000");
	    }
	}

    /* Ignore new coordinate system if it is not supported */
    else {
	if ((sysout = vimoswcscsys (coorsys)) < 0)
	return;

	/* Do not try to convert linear or alt-az coordinates */
	if (sysout != vimoswcs->sysvimoswcs &&
	    (vimoswcs->sysvimoswcs == VIMOSWCS_LINEAR || vimoswcs->sysvimoswcs == VIMOSWCS_ALTAZ))
	    return;

	strcpy (vimoswcs->radecout, coorsys);
	vimoswcs->eqout = vimoswcsceq (coorsys);
	}

    vimoswcs->sysout = sysout;
    if (vimoswcs->vimoswcson) {

	/* Set output in degrees flag and number of decimal places */
	if (vimoswcs->sysout == VIMOSWCS_GALACTIC || vimoswcs->sysout == VIMOSWCS_ECLIPTIC ||
	    vimoswcs->sysout == VIMOSWCS_PLANET) {
	    vimoswcs->degout = 1;
	    vimoswcs->ndec = 5;
	    }
	else if (vimoswcs->sysout == VIMOSWCS_ALTAZ) {
	    vimoswcs->degout = 1;
	    vimoswcs->ndec = 5;
	    }
	else if (vimoswcs->sysout == VIMOSWCS_NPOLE || vimoswcs->sysout == VIMOSWCS_SPA) {
	    vimoswcs->degout = 1;
	    vimoswcs->ndec = 5;
	    }
	else {
	    vimoswcs->degout = 0;
	    vimoswcs->ndec = 3;
	    }
	}
    return;
}


/* Return current value of WCS output coordinate system set by -vimoswcsout */
char *
getvimoswcsout(vimoswcs)
struct	WorldCoor *vimoswcs; /* World coordinate system structure */
{
    if (novimoswcs (vimoswcs))
	return (NULL);
    else
	return(vimoswcs->radecout);
}


/* Initialize WCS input coordinate system for use by VIMOSWCS2PIX() */

void
vimoswcsininit (vimoswcs, coorsys)

struct WorldCoor *vimoswcs;	/* World coordinate system structure */
char	*coorsys;	/* Input world coordinate system:
			   FK4, FK5, B1950, J2000, GALACTIC, ECLIPTIC
			   fk4, fk5, b1950, j2000, galactic, ecliptic */
{
    int sysin, i;

    if (novimoswcs (vimoswcs))
	return;

    /* If argument is null, set to image system and equinox */
    if (coorsys == NULL || strlen (coorsys) < 1) {
	vimoswcs->sysin = vimoswcs->sysvimoswcs;
	strcpy (vimoswcs->radecin, vimoswcs->radecsys);
	vimoswcs->eqin = vimoswcs->equinox;
	if (vimoswcs->sysin == VIMOSWCS_B1950) {
	    if (vimoswcs->eqin != 1950.0) {
		vimoswcs->radecin[0] = 'B';
		sprintf (vimoswcs->radecin+1,"%.4f", vimoswcs->equinox);
		i = strlen(vimoswcs->radecin) - 1;
		if (vimoswcs->radecin[i] == '0')
		    vimoswcs->radecin[i] = (char)0;
		i = strlen(vimoswcs->radecin) - 1;
		if (vimoswcs->radecin[i] == '0')
		    vimoswcs->radecin[i] = (char)0;
		i = strlen(vimoswcs->radecin) - 1;
		if (vimoswcs->radecin[i] == '0')
		    vimoswcs->radecin[i] = (char)0;
		}
	    else
		strcpy (vimoswcs->radecin, "B1950");
	    }
	else if (vimoswcs->sysin == VIMOSWCS_J2000) {
	    if (vimoswcs->eqin != 2000.0) {
		vimoswcs->radecin[0] = 'J';
		sprintf (vimoswcs->radecin+1,"%.4f", vimoswcs->equinox);
		i = strlen(vimoswcs->radecin) - 1;
		if (vimoswcs->radecin[i] == '0')
		    vimoswcs->radecin[i] = (char)0;
		i = strlen(vimoswcs->radecin) - 1;
		if (vimoswcs->radecin[i] == '0')
		    vimoswcs->radecin[i] = (char)0;
		i = strlen(vimoswcs->radecin) - 1;
		if (vimoswcs->radecin[i] == '0')
		    vimoswcs->radecin[i] = (char)0;
		}
	    else
		strcpy (vimoswcs->radecin, "J2000");
	    }
	}

    /* Ignore new coordinate system if it is not supported */
    if ((sysin = vimoswcscsys (coorsys)) < 0)
	return;

    vimoswcs->sysin = sysin;
    vimoswcs->eqin = vimoswcsceq (coorsys);
    strcpy (vimoswcs->radecin, coorsys);
    return;
}


/* Return current value of WCS input coordinate system set by vimoswcsininit */
char *
getvimoswcsin (vimoswcs)
struct	WorldCoor *vimoswcs; /* World coordinate system structure */
{
    if (novimoswcs (vimoswcs))
	return (NULL);
    else
	return (vimoswcs->radecin);
}


/* Set WCS output in degrees or hh:mm:ss dd:mm:ss, returning old flag value */
int
setvimoswcsdeg(vimoswcs, new)
struct	WorldCoor *vimoswcs;	/* World coordinate system structure */
int	new;		/* 1: degrees, 0: h:m:s d:m:s */
{
    int old;

    if (novimoswcs (vimoswcs))
	return (0);
    old = vimoswcs->degout;
    vimoswcs->degout = new;
    if (new == 1 && old == 0 && vimoswcs->ndec == 3)
	vimoswcs->ndec = 6;
    if (new == 0 && old == 1 && vimoswcs->ndec == 5)
	vimoswcs->ndec = 3;
    return(old);
}


/* Set number of decimal places in pix2vimoswcst output string */
int
vimoswcsndec (vimoswcs, ndec)
struct	WorldCoor *vimoswcs;	/* World coordinate system structure */
int	ndec;		/* Number of decimal places in output string */
			/* If < 0, return current unchanged value */
{
    if (novimoswcs (vimoswcs))
	return (0);
    else if (ndec >= 0)
	vimoswcs->ndec = ndec;
    return (vimoswcs->ndec);
}



/* Return current value of coordinate system */
char *
getradecsys(vimoswcs)
struct     WorldCoor *vimoswcs; /* World coordinate system structure */
{
    if (novimoswcs (vimoswcs))
	return (NULL);
    else
	return (vimoswcs->radecsys);
}


/* Set output string mode for LINEAR coordinates */

void
setvimoswcslin (vimoswcs, mode)
struct	WorldCoor *vimoswcs; /* World coordinate system structure */
int	mode;		/* mode = 0: x y linear
			   mode = 1: x units x units
			   mode = 2: x y linear units */
{
    if (isvimoswcs (vimoswcs))
	vimoswcs->linmode = mode;
    return;
}


/* Convert pixel coordinates to World Coordinate string */

int
pix2vimoswcst (vimoswcs, xpix, ypix, vimoswcstring, lstr)

struct	WorldCoor *vimoswcs;	/* World coordinate system structure */
double	xpix,ypix;	/* Image coordinates in pixels */
char	*vimoswcstring;	/* World coordinate string (returned) */
int	lstr;		/* Length of world coordinate string (returned) */
{
	double	xpos,ypos;
	char	rastr[32], decstr[32];
	int	minlength, lunits, lstring;

	if (novimoswcs (vimoswcs)) {
	    if (lstr > 0)
		vimoswcstring[0] = 0;
	    return(0);
	    }

	pix2vimoswcs (vimoswcs,xpix,ypix,&xpos,&ypos);

	/* Keep ra/longitude within range
	if (xpos < 0.0)
	    xpos = xpos + 360.0;

	else if (xpos > 360.0)
	    xpos = xpos - 360.0; */

	/* If point is off scale, set string accordingly */
	if (vimoswcs->offscl) {
	    (void)sprintf (vimoswcstring,"Off map");
	    return (1);
	    }

	/* Print coordinates in degrees */
	else if (vimoswcs->degout == 1) {
	    minlength = 9 + (2 * vimoswcs->ndec);
	    if (lstr > minlength) {
		deg2str (rastr, 32, xpos, vimoswcs->ndec);
		deg2str (decstr, 32, ypos, vimoswcs->ndec);
		if (vimoswcs->tabsys)
		    (void)sprintf (vimoswcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (vimoswcstring,"%s %s", rastr, decstr);
		lstr = lstr - minlength;
		}
	    else {
		if (vimoswcs->tabsys)
		    strncpy (vimoswcstring,"*********	**********",lstr);
		else
		    strncpy (vimoswcstring,"*******************",lstr);
		lstr = 0;
		}
	    }

	/* print coordinates in sexagesimal notation */
	else if (vimoswcs->degout == 0) {
	    minlength = 18 + (2 * vimoswcs->ndec);
	    if (lstr > minlength) {
		if (vimoswcs->sysout == VIMOSWCS_J2000 || vimoswcs->sysout == VIMOSWCS_B1950) {
		    ra2str (rastr, 32, xpos, vimoswcs->ndec);
		    dec2str (decstr, 32, ypos, vimoswcs->ndec-1);
		    }
		else {
		    dec2str (rastr, 32, xpos, vimoswcs->ndec);
		    dec2str (decstr, 32, ypos, vimoswcs->ndec);
		    }
		if (vimoswcs->tabsys)
		    (void)sprintf (vimoswcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (vimoswcstring,"%s %s", rastr, decstr);
	        lstr = lstr - minlength;
		}
	    else {
		if (vimoswcs->tabsys)
		    strncpy (vimoswcstring,"*************	*************",lstr);
		else
		    strncpy (vimoswcstring,"**************************",lstr);
		lstr = 0;
		}
	    }

	/* Label galactic coordinates */
	if (vimoswcs->sysout == VIMOSWCS_GALACTIC) {
	    if (lstr > 9 && vimoswcs->printsys)
		strcat (vimoswcstring," galactic");
	    }

	/* Label ecliptic coordinates */
	else if (vimoswcs->sysout == VIMOSWCS_ECLIPTIC) {
	    if (lstr > 9 && vimoswcs->printsys) {
		if (vimoswcs->tabsys)
		    strcat (vimoswcstring,"	ecliptic");
		else
		    strcat (vimoswcstring," ecliptic");
		}
	    }

	/* Label planet coordinates */
	else if (vimoswcs->sysout == VIMOSWCS_PLANET) {
	    if (lstr > 9 && vimoswcs->printsys) {
		if (vimoswcs->tabsys)
		    strcat (vimoswcstring,"	planet");
		else
		    strcat (vimoswcstring," planet");
		}
	    }

	/* Label alt-az coordinates */
	else if (vimoswcs->sysout == VIMOSWCS_ALTAZ) {
	    if (lstr > 7 && vimoswcs->printsys) {
		if (vimoswcs->tabsys)
		    strcat (vimoswcstring,"	alt-az");
		else
		    strcat (vimoswcstring," alt-az");
		}
	    }

	/* Label north pole angle coordinates */
	else if (vimoswcs->sysout == VIMOSWCS_NPOLE) {
	    if (lstr > 7 && vimoswcs->printsys) {
		if (vimoswcs->tabsys)
		    strcat (vimoswcstring,"	long-npa");
		else
		    strcat (vimoswcstring," long-npa");
		}
	    }

	/* Label south pole angle coordinates */
	else if (vimoswcs->sysout == VIMOSWCS_SPA) {
	    if (lstr > 7 && vimoswcs->printsys) {
		if (vimoswcs->tabsys)
		    strcat (vimoswcstring,"	long-spa");
		else
		    strcat (vimoswcstring," long-spa");
		}
	    }

	/* Label equatorial coordinates */
	else if (vimoswcs->sysout==VIMOSWCS_B1950 || vimoswcs->sysout==VIMOSWCS_J2000) {
	    if (lstr > strlen(vimoswcs->radecout)+1 && vimoswcs->printsys) {
		if (vimoswcs->tabsys)
		    strcat (vimoswcstring,"	");
		else
		    strcat (vimoswcstring," ");
		strcat (vimoswcstring, vimoswcs->radecout);
		}
	    }

	/* Output linear coordinates */
	else {
	    num2str (rastr, xpos, 0, vimoswcs->ndec);
	    num2str (decstr, ypos, 0, vimoswcs->ndec);
	    lstring = strlen (rastr) + strlen (decstr) + 1;
	    lunits = strlen (vimoswcs->units[0]) + strlen (vimoswcs->units[1]) + 2;
	    if (vimoswcs->sysvimoswcs == VIMOSWCS_LINEAR && vimoswcs->linmode == 1) {
		if (lstr > lstring + lunits) {
		    if (strlen (vimoswcs->units[0]) > 0) {
			strcat (rastr, " ");
			strcat (rastr, vimoswcs->units[0]);
			}
		    if (strlen (vimoswcs->units[1]) > 0) {
			strcat (decstr, " ");
			strcat (decstr, vimoswcs->units[1]);
			}
		    lstring = lstring + lunits;
		    }
		}
	    if (lstr > lstring) {
		if (vimoswcs->tabsys)
		    (void)sprintf (vimoswcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (vimoswcstring,"%s %s", rastr, decstr);
		}
	    else {
		if (vimoswcs->tabsys)
		    strncpy (vimoswcstring,"**********	*********",lstr);
		else
		    strncpy (vimoswcstring,"*******************",lstr);
		}
	    if (vimoswcs->sysvimoswcs == VIMOSWCS_LINEAR && vimoswcs->linmode != 1 &&
		lstr > lstring + 7)
		strcat (vimoswcstring, " linear");
	    if (vimoswcs->sysvimoswcs == VIMOSWCS_LINEAR && vimoswcs->linmode == 2 &&
		lstr > lstring + lunits + 7) {
		if (strlen (vimoswcs->units[0]) > 0) {
		    strcat (vimoswcstring, " ");
		    strcat (vimoswcstring, vimoswcs->units[0]);
		    }
		if (strlen (vimoswcs->units[1]) > 0) {
		    strcat (vimoswcstring, " ");
		    strcat (vimoswcstring, vimoswcs->units[1]);
		    }
		    
		}
	    }
	return (1);
}


/* Convert pixel coordinates to World Coordinates */

void
pix2vimoswcs (

struct WorldCoor *vimoswcs,		/* World coordinate system structure */
double	xpix,
double ypix,	/* x and y image coordinates in pixels */
double	*xpos,
double *ypos)	/* RA and Dec in degrees (returned) */
{
    double	xpi, ypi, xp, yp;
    double	eqin, eqout;
    int vimoswcspos();
    extern int dsspos();
    extern int platepos();
    extern int worldpos();
    extern int tnxpos();
    extern void vimoswcscon();

    if (novimoswcs (vimoswcs))
	return;
    vimoswcs->xpix = xpix;
    vimoswcs->ypix = ypix;
    vimoswcs->zpix = zpix;
    vimoswcs->offscl = 0;

    /* If this WCS is converted from another WCS rather than pixels, convert now */
    if (vimoswcs->vimoswcs != NULL) {
	pix2vimoswcs (vimoswcs->vimoswcs, xpix, ypix, &xpi, &ypi);
	}
    else {
	xpi = xpix;
	ypi = ypix;
	}

    /* Convert image coordinates to sky coordinates */

    /* Use Digitized Sky Survey plate fit */
    if (vimoswcs->prjcode == VIMOSWCS_DSS) {
	if (dsspos (xpi, ypi, vimoswcs, &xp, &yp))
	    vimoswcs->offscl = 1;
	}

    /* Use SAO plate fit */
    else if (vimoswcs->prjcode == VIMOSWCS_PLT) {
	if (platepos (xpi, ypi, vimoswcs, &xp, &yp))
	    vimoswcs->offscl = 1;
	}

    /* Use NOAO IRAF corrected plane tangent projection */
    else if (vimoswcs->prjcode == VIMOSWCS_TNX) {
	if (tnxpos (xpi, ypi, vimoswcs, &xp, &yp))
	    vimoswcs->offscl = 1;
	}

    /* Use Classic AIPS projections */
    else if (vimoswcs->vimoswcsproj == VIMOSWCS_OLD || vimoswcs->prjcode <= 0) {
	if (worldpos (xpi, ypi, vimoswcs, &xp, &yp))
	    vimoswcs->offscl = 1;
	}

    /* Use Mark Calabretta's WCSLIB projections */
    else if (vimoswcspos (xpi, ypi, vimoswcs, &xp, &yp))
	vimoswcs->offscl = 1;
	    	

    /* Do not change coordinates if offscale */
    if (vimoswcs->offscl) {
	*xpos = 0.0;
	*ypos = 0.0;
	}
    else {

	/* Convert coordinates to output system, if not LINEAR */
        if (vimoswcs->prjcode > 0) {

	    /* Convert coordinates to desired output system */
	    eqin = vimoswcs->equinox;
	    eqout = vimoswcs->eqout;
	    vimoswcscon (vimoswcs->sysvimoswcs,vimoswcs->sysout,eqin,eqout,&xp,&yp,vimoswcs->epoch);
	    }
	if (vimoswcs->latbase == 90)
	    yp = 90.0 - yp;
	else if (vimoswcs->latbase == -90)
	    yp = yp - 90.0;
	vimoswcs->xpos = xp;
	vimoswcs->ypos = yp;
	*xpos = xp;
	*ypos = yp;
	}
    return;
}


/* Convert World Coordinates to pixel coordinates */

void
vimoswcs2pix (vimoswcs, xpos, ypos, xpix, ypix, offscl)

struct WorldCoor *vimoswcs;	/* World coordinate system structure */
double	xpos,ypos;	/* World coordinates in degrees */
double	*xpix,*ypix;	/* Image coordinates in pixels */
int	*offscl;	/* 0 if within bounds, else off scale */
{
    vimoswcsc2pix (vimoswcs, xpos, ypos, vimoswcs->radecin, xpix, ypix, offscl);
    return;
}

/* Convert World Coordinates to pixel coordinates */

void
vimoswcsc2pix (vimoswcs, xpos, ypos, coorsys, xpix, ypix, offscl)

struct WorldCoor *vimoswcs;	/* World coordinate system structure */
double	xpos,ypos;	/* World coordinates in degrees */
char	*coorsys;	/* Input world coordinate system:
			   FK4, FK5, B1950, J2000, GALACTIC, ECLIPTIC
			   fk4, fk5, b1950, j2000, galactic, ecliptic
			   * If NULL, use image WCS */
double	*xpix,*ypix;	/* Image coordinates in pixels */
int	*offscl;	/* 0 if within bounds, else off scale */
{
    double xp, yp, xpi, ypi;
    double eqin, eqout;
    int sysin;
    int vimoswcspix();
    extern int dsspix();
    extern int platepix();
    extern int worldpix();
    extern int tnxpix();
    extern void vimoswcscon();

    if (novimoswcs (vimoswcs))
	return;

    *offscl = 0;
    xp = xpos;
    yp = ypos;
    if (vimoswcs->latbase == 90)
	yp = 90.0 - yp;
    else if (vimoswcs->latbase == -90)
	yp = yp - 90.0;
    if (coorsys == NULL) {
	sysin = vimoswcs->sysvimoswcs;
	eqin = vimoswcs->equinox;
	}
    else {
	sysin = vimoswcscsys (coorsys);
	eqin = vimoswcsceq (coorsys);
	}
    vimoswcs->zpix = 1.0;

    /* Convert coordinates to same system as image */
    eqout = vimoswcs->equinox;
    vimoswcscon (sysin, vimoswcs->sysvimoswcs, eqin, eqout, &xp, &yp, vimoswcs->epoch);

    /* Convert sky coordinates to image coordinates */

    /* Use Digitized Sky Survey plate fit */
    if (vimoswcs->prjcode == VIMOSWCS_DSS) {
	if (dsspix (xp, yp, vimoswcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use SAO polynomial plate fit */
    else if (vimoswcs->prjcode == VIMOSWCS_PLT) {
	if (platepix (xp, yp, vimoswcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use NOAO IRAF corrected plane tangent projection */
    else if (vimoswcs->prjcode == VIMOSWCS_TNX) {
	if (tnxpix (xp, yp, vimoswcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use Classic AIPS projections */
    else if (vimoswcs->vimoswcsproj == VIMOSWCS_OLD || vimoswcs->prjcode <= 0) {
	if (worldpix (xp, yp, vimoswcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use Mark Calabretta's WCSLIB projections */
    else if (vimoswcspix (xp, yp, vimoswcs, &xpi, &ypi)) {
	*offscl = 1;
	}

    /* If this WCS is converted from another WCS rather than pixels, convert now */
    if (vimoswcs->vimoswcs != NULL) {
	vimoswcsc2pix (vimoswcs->vimoswcs, xpi, ypi, NULL, xpix, ypix, offscl);
	}
    else {
	*xpix = xpi;
	*ypix = ypi;

	/* Set off-scale flag to 2 if off image but within bounds of projection */
	if (!*offscl) {
	    if (xpi < 0.5 || ypi < 0.5)
		*offscl = 2;
	    else if (xpi > vimoswcs->nxpix + 0.5 || ypi > vimoswcs->nypix + 0.5)
		*offscl = 2;
	    }
	}

    vimoswcs->offscl = *offscl;
    vimoswcs->xpos = xpos;
    vimoswcs->ypos = ypos;
    vimoswcs->xpix = *xpix;
    vimoswcs->ypix = *ypix;
    return;
}


int
vimoswcspos (xpix, ypix, vimoswcs, xpos, ypos)

/* Input: */
double  xpix;          /* x pixel number  (RA or long without rotation) */
double  ypix;          /* y pixel number  (dec or lat without rotation) */
struct WorldCoor *vimoswcs;  /* WCS parameter structure */

/* Output: */
double  *xpos;           /* x (RA) coordinate (deg) */
double  *ypos;           /* y (dec) coordinate (deg) */
{
    int offscl;
    int i;
    int vimoswcsrev();
    double vimoswcscrd[4], imgcrd[4], pixcrd[4];
    double phi, theta;
    
    *xpos = 0.0;
    *ypos = 0.0;

    pixcrd[0] = xpix;
    pixcrd[1] = ypix;
    if (vimoswcs->prjcode == VIMOSWCS_CSC || vimoswcs->prjcode == VIMOSWCS_QSC ||
	vimoswcs->prjcode == VIMOSWCS_TSC)
	pixcrd[2] = (double) (izpix + 1);
    else
	pixcrd[2] = zpix;
    pixcrd[3] = 1.0;
    for (i = 0; i < 4; i++)
	imgcrd[i] = 0.0;
    offscl = vimoswcsrev (vimoswcs->ctype, &vimoswcs->vimoswcsl, pixcrd, &vimoswcs->lin, imgcrd,
		    &vimoswcs->prj, &phi, &theta, vimoswcs->crval, &vimoswcs->cel, vimoswcscrd);
    if (offscl == 0) {
	*xpos = vimoswcscrd[vimoswcs->vimoswcsl.lng];
	*ypos = vimoswcscrd[vimoswcs->vimoswcsl.lat];
	}
    return (offscl);
}

int
vimoswcspix (xpos, ypos, vimoswcs, xpix, ypix)

/* Input: */
double  xpos;           /* x (RA) coordinate (deg) */
double  ypos;           /* y (dec) coordinate (deg) */
struct WorldCoor *vimoswcs;  /* WCS parameter structure */

/* Output: */
double  *xpix;          /* x pixel number  (RA or long without rotation) */
double  *ypix;          /* y pixel number  (dec or lat without rotation) */

{
    int offscl;
    int vimoswcsfwd();
    double vimoswcscrd[4], imgcrd[4], pixcrd[4];
    double phi, theta;

    *xpix = 0.0;
    *ypix = 0.0;
    if (vimoswcs->vimoswcsl.flag != VIMOSWCSSET) {
	if (vimoswcsset (vimoswcs->lin.naxis, vimoswcs->ctype, &vimoswcs->vimoswcsl) )
	    return (1);
	}

    /* Set input for WCSLIB subroutines */
    vimoswcscrd[0] = 0.0;
    vimoswcscrd[1] = 0.0;
    vimoswcscrd[2] = 0.0;
    vimoswcscrd[3] = 0.0;
    vimoswcscrd[vimoswcs->vimoswcsl.lng] = xpos;
    vimoswcscrd[vimoswcs->vimoswcsl.lat] = ypos;

    /* Initialize output for WCSLIB subroutines */
    pixcrd[0] = 0.0;
    pixcrd[1] = 0.0;
    pixcrd[2] = 1.0;
    pixcrd[3] = 1.0;
    imgcrd[0] = 0.0;
    imgcrd[1] = 0.0;
    imgcrd[2] = 1.0;
    imgcrd[3] = 1.0;

    /* Invoke WCSLIB subroutines for coordinate conversion */
    offscl = vimoswcsfwd (vimoswcs->ctype, &vimoswcs->vimoswcsl, vimoswcscrd, vimoswcs->crval, &vimoswcs->cel,
		     &phi, &theta, &vimoswcs->prj, imgcrd, &vimoswcs->lin, pixcrd);

    if (!offscl) {
	*xpix = pixcrd[0];
	*ypix = pixcrd[1];
	if (vimoswcs->prjcode == VIMOSWCS_CSC || vimoswcs->prjcode == VIMOSWCS_QSC ||
	    vimoswcs->prjcode == VIMOSWCS_TSC)
	    vimoswcs->zpix = pixcrd[2] - 1.0;
	else
	    vimoswcs->zpix = pixcrd[2];
	}
    return (offscl);
}


/* Set third dimension for cube projections */

int
vimoswcszin (izpix0)

int	izpix0;		/* coordinate in third dimension
			   (if < 0, return current value without changing it */
{
    if (izpix0 > -1) {
	izpix = izpix0;
	zpix = (double) izpix0;
	}
    return (izpix);
}


/* Return third dimension for cube projections */

int
vimoswcszout (vimoswcs)

struct WorldCoor *vimoswcs;  /* WCS parameter structure */
{
    return ((int) vimoswcs->zpix);
}

/* Set file name for error messages */
void
setvimoswcsfile (
char *filename)	/* FITS or IRAF file with WCS */
{   if (strlen (filename) < 256)
	strcpy (vimoswcsfile, filename);
    else
	strncpy (vimoswcsfile, filename, 255);
    return; }

/* Set error message */
void
setvimoswcserr (
char *errmsg)	/* Error mesage  < 80 char */
{ strcpy (vimoswcserrmsg, errmsg); return; }

/* Print error message */
void
vimoswcserr()
{   if (strlen (vimoswcsfile) > 0)
	fprintf (stderr, "%s in file %s\n",vimoswcserrmsg, vimoswcsfile);
    else
	fprintf (stderr, "%s\n",vimoswcserrmsg);
    return; }


/* Flag to use AIPS WCS subroutines instead of WCSLIB */
void
setdefvimoswcs (
int wp)
{ vimoswcsproj0 = wp; return; }

int
getdefvimoswcs ()
{ return (vimoswcsproj0); }

/* Save output default coordinate system */
static char vimoswcscoor0[16];

void
savevimoswcscoor (
char *vimoswcscoor)
{ strcpy (vimoswcscoor0, vimoswcscoor); return; }

/* Return preset output default coordinate system */
char * getvimoswcscoor ()
{ return (vimoswcscoor0); }


/* Save default commands */
static char *vimoswcscom0[10];

void
savevimoswcscom (
int i,
char *vimoswcscomm)
{
    int lcom;
    if (i < 0) i = 0;
    else if (i > 9) i = 9;
    lcom = strlen (vimoswcscomm) + 2;
    vimoswcscom0[i] = (char *) calloc (lcom, 1);
    if (vimoswcscom0[i] != NULL)
	strcpy (vimoswcscom0[i], vimoswcscomm);
    return;
}

void
setvimoswcscom (
struct WorldCoor *vimoswcs)  /* WCS parameter structure */
{
    char envar[16];
    int i;
    char *str;
    if (novimoswcs(vimoswcs))
	return;
    for (i = 0; i < 10; i++) {
	if (i == 0)
	    strcpy (envar, "WCS_COMMAND");
	else
	    sprintf (envar, "WCS_COMMAND%d", i);
	if (vimoswcscom0[i] != NULL)
	    vimoswcscominit (vimoswcs, i, vimoswcscom0[i]);
	else if ((str = getenv (envar)) != NULL)
	    vimoswcscominit (vimoswcs, i, str);
	else if (i == 1)
	    vimoswcscominit (vimoswcs, i, "suac -ah %s");	/* F1= Search USNO A Catalog */
	else if (i == 2)
	    vimoswcscominit (vimoswcs, i, "sgsc -ah %s");	/* F2= Search HST GSC */
	else if (i == 3)
	    vimoswcscominit (vimoswcs, i, "sact -ah %s"); /* F3= Search USNO ACT Catalog */
	else if (i == 4)
	    vimoswcscominit (vimoswcs, i, "sppm -ah %s");	/* F4= Search PPM Catalog */
	else if (i == 5)
	    vimoswcscominit (vimoswcs, i, "ssao -ah %s");	/* F5= Search SAO Catalog */
	else
	    vimoswcs->command_format[i] = NULL;
	}
    return;
}

char *
getvimoswcscom (
int i)
{ return (vimoswcscom0[i]); }


void
freevimoswcscom (
struct WorldCoor *vimoswcs)  /* WCS parameter structure */
{
    int i;
    for (i = 0; i < 10; i++) {
	if (vimoswcscom0[i] != NULL) {
	    free (vimoswcscom0[i]);
	    vimoswcscom0[i] = NULL;
	    }
	}
    if (isvimoswcs (vimoswcs)) {
	for (i = 0; i < 10; i++) {
	    if (vimoswcs->command_format[i] != NULL) {
		free (vimoswcs->command_format[i]);
		}
	    }
	}
    return;
}


/* Oct 28 1994	new program
 * Dec 21 1994	Implement CD rotation matrix
 * Dec 22 1994	Allow RA and DEC to be either x,y or y,x
 *
 * Mar  6 1995	Add Digital Sky Survey plate fit
 * May  2 1995	Add prototype of PIX2WCST to WCSCOM
 * May 25 1995	Print leading zero for hours and degrees
 * Jun 21 1995	Add WCS2PIX to get pixels from WCS
 * Jun 21 1995	Read plate scale from FITS header for plate solution
 * Jul  6 1995	Pass WCS structure as argument; malloc it in WCSINIT
 * Jul  6 1995	Check string lengths in PIX2WCST
 * Aug 16 1995	Add galactic coordinate conversion to PIX2WCST
 * Aug 17 1995	Return 0 from isvimoswcs if vimoswcs structure is not yet set
 * Sep  8 1995	Do not include malloc.h if VMS
 * Sep  8 1995	Check for legal WCS before trying anything
 * Sep  8 1995	Do not try to set WCS if missing key keywords
 * Oct 18 1995	Add WCSCENT and WCSDIST to print center and size of image
 * Nov  6 1995	Include stdlib.h instead of malloc.h
 * Dec  6 1995	Fix format statement in PIX2WCST
 * Dec 19 1995	Change MALLOC to CALLOC to initialize array to zeroes
 * Dec 19 1995	Explicitly initialize rotation matrix and yinc
 * Dec 22 1995	If SECPIX is set, use approximate WCS
 * Dec 22 1995	Always print coordinate system
 *
 * Jan 12 1996	Use plane-tangent, not linear, projection if SECPIX is set
 * Jan 12 1996  Add WCSSET to set WCS without an image
 * Feb 15 1996	Replace all calls to HGETC with HGETS
 * Feb 20 1996	Add tab table output from PIX2WCST
 * Apr  2 1996	Convert all equinoxes to B1950 or J2000
 * Apr 26 1996	Get and use image epoch for accurate FK4/FK5 conversions
 * May 16 1996	Clean up internal documentation
 * May 17 1996	Return width in right ascension degrees, not sky degrees
 * May 24 1996	Remove extraneous print command from WCSSIZE
 * May 28 1996	Add NOWCS and WCSSHIFT subroutines
 * Jun 11 1996	Drop unused variables after running lint
 * Jun 12 1996	Set equinox as well as system in WCSSHIFT
 * Jun 14 1996	Make DSS keyword searches more robust
 * Jul  1 1996	Allow for SECPIX1 and SECPIX2 keywords
 * Jul  2 1996	Test for CTYPE1 instead of CRVAL1
 * Jul  5 1996	Declare all subroutines in vimoswcs.h
 * Jul 19 1996	Add subroutine WCSFULL to return real image size
 * Aug 12 1996	Allow systemless coordinates which cannot be converted
 * Aug 15 1996	Allow LINEAR WCS to pass numbers through transparently
 * Aug 15 1996	Add WCSERR to print error message under calling program control
 * Aug 16 1996	Add latitude and longitude as image coordinate types
 * Aug 26 1996	Fix arguments to HLENGTH in WCSNINIT
 * Aug 28 1996	Explicitly set OFFSCL in WCS2PIX if coordinates outside image
 * Sep  3 1996	Return computed pixel values even if they are offscale
 * Sep  6 1996	Allow filename to be passed by WCSCOM
 * Oct  8 1996	Default to 2000 for EQUINOX and EPOCH and FK5 for RADECSYS
 * Oct  8 1996	If EPOCH is 0 and EQUINOX is not set, default to 1950 and FK4
 * Oct 15 1996  Add comparison when testing an assignment
 * Oct 16 1996  Allow PIXEL CTYPE which means WCS is same as image coordinates
 * Oct 21 1996	Add WCS_COMMAND environment variable
 * Oct 25 1996	Add image scale to WCSCENT
 * Oct 30 1996	Fix bugs in WCS2PIX
 * Oct 31 1996	Fix CD matrix rotation angle computation
 * Oct 31 1996	Use inline degree <-> radian conversion functions
 * Nov  1 1996	Add option to change number of decimal places in PIX2WCST
 * Nov  5 1996	Set wcs->crot to 1 if rotation matrix is used
 * Dec  2 1996	Add altitide/azimuth coordinates
 * Dec 13 1996	Fix search format setting from environment
 *
 * Jan 22 1997	Add ifdef for Eric Mandel (SAOtng)
 * Feb  5 1997	Add wcsout for Eric Mandel
 * Mar 20 1997	Drop unused variable STR in WCSCOM
 * May 21 1997	Do not make pixel coordinates mod 360 in PIX2WCST
 * May 22 1997	Add PIXEL prjcode = -1;
 * Jul 11 1997	Get center pixel x and y from header even if no WCS
 * Aug  7 1997	Add NOAO PIXSCALi keywords for default WCS
 * Oct 15 1997	Do not reset reference pixel in WCSSHIFT
 * Oct 20 1997	Set chip rotation
 * Oct 24 1997	Keep longitudes between 0 and 360, not -180 and +180
 * Nov  5 1997	Do no compute crot and srot; they are now computed in worldpos
 * Dec 16 1997	Set rotation and axis increments from CD matrix
 *
 * Jan  6 1998	Deal with J2000 and B1950 as EQUINOX values (from ST)
 * Jan  7 1998	Read INSTRUME and DETECTOR header keywords
 * Jan  7 1998	Fix tab-separated output
 * Jan  9 1998	Precess coordinates if either FITS projection or *DSS plate*
 * Jan 16 1998	Change PTYPE to not include initial hyphen
 * Jan 16 1998	Change WCSSET to WCSXINIT to avoid conflict with Calabretta
 * Jan 23 1998	Change PCODE to PRJCODE to avoid conflict with Calabretta
 * Jan 27 1998	Add LATPOLE and LONGPOLE for Calabretta projections
 * Feb  5 1998	Make cd and dc into vectors; use matinv() to invert cd
 * Feb  5 1998	In wcsoutinit(), check that corsys is a valid pointer
 * Feb 18 1998	Fix bugs for Calabretta projections
 * Feb 19 1998	Add wcs structure access subroutines from Eric Mandel
 * Feb 19 1998	Add wcsreset() to make sure derived values are reset
 * Feb 19 1998	Always set oldwcs to 1 if NCP projection
 * Feb 19 1998	Add subroutine to set oldwcs default
 * Feb 20 1998	Initialize projection types one at a time for SunOS C
 * Feb 23 1998	Add TNX projection from NOAO; treat it as TAN
 * Feb 23 1998	Compute size based on max and min coordinates, not sides
 * Feb 26 1998	Add code to set center pixel if part of detector array
 * Mar  6 1998	Write 8-character values to RADECSYS
 * Mar  9 1998	Add naxis to WCS structure
 * Mar 16 1998	Use separate subroutine for IRAF TNX projection
 * Mar 20 1998	Set PC matrix if more than two axes and it's not in header
 * Mar 20 1998	Reset lin flag in WCSRESET if CDELTn
 * Mar 20 1998	Set CD matrix with CDELTs and CROTA in wcsinit and wcsreset
 * Mar 20 1998	Allow initialization of rotation angle alone
 * Mar 23 1998	Use dsspos() and dsspix() instead of platepos() and platepix()
 * Mar 24 1998	Set up PLT/PLATE plate polynomial fit using platepos() and platepix()
 * Mar 25 1998	Read plate fit coefficients from header in getwcscoeffs()
 * Mar 27 1998	Check for FITS WCS before DSS WCS
 * Mar 27 1998	Compute scale from edges if xinc and yinc not set in wcscent()
 * Apr  6 1998	Change plate coefficient keywords from PLTij to COi_j
 * Apr  6 1998	Read all coefficients in line instead of with subroutine
 * Apr  7 1998	Change amd_i_coeff to i_coeff
 * Apr  8 1998	Add wcseqset to change equinox after wcs has been set
 * Apr 10 1998	Use separate counters for x and y polynomial coefficients
 * Apr 13 1998	Use CD/CDELT+CDROTA if oldwcs is set
 * Apr 14 1998	Use codes instead of strings for various coordinate systems
 * Apr 14 1998	Separate input coordinate conversion from output conversion
 * Apr 14 1998	Use wcscon() for most coordinate conversion
 * Apr 17 1998	Always compute cdelt[]
 * Apr 17 1998	Deal with reversed axis more correctly
 * Apr 17 1998	Compute rotation angle and approximate CDELTn for polynomial
 * Apr 23 1998	Deprecate xref/yref in favor of crval[]
 * Apr 23 1998	Deprecate xrefpix/yrefpix in favor of crpix[]
 * Apr 23 1998	Add LINEAR to coordinate system types
 * Apr 23 1998	Always use AIPS subroutines for LINEAR or PIXEL
 * Apr 24 1998	Format linear coordinates better
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * Apr 28 1998  Change projection flags to WCS_*
 * Apr 28 1998	Add subroutine wcsc2pix for coordinates to pixels with system
 * Apr 28 1998	Add setlinmode() to set output string mode for LINEAR coordinates
 * Apr 30 1998	Fix bug by setting degree flag for lat and long in wcsinit()
 * Apr 30 1998	Allow leading "-"s in projecting in wcsxinit()
 * May  1 1998	Assign new output coordinate system only if legitimate system
 * May  1 1998	Do not allow oldwcs=1 unless allowed projection
 * May  4 1998	Fix bug in units reading for LINEAR coordinates
 * May  6 1998	Initialize to no CD matrix
 * May  6 1998	Use TAN instead of TNX if oldwcs
 * May 12 1998	Set 3rd and 4th coordinates in wcspos()
 * May 12 1998	Return *xpos and *ypos = 0 in pix2wcs() if offscale
 * May 12 1998	Declare undeclared external subroutines after lint
 * May 13 1998	Add equinox conversion to specified output equinox
 * May 13 1998	Set output or input system to image with null argument
 * May 15 1998	Return reference pixel, cdelts, and rotation for DSS
 * May 20 1998	Fix bad bug so setlinmode() is no-op if wcs not set
 * May 20 1998	Fix bug so getwcsout() returns null pointer if wcs not set
 * May 27 1998	Change WCS_LPR back to WCS_LIN; allow CAR in oldwcs
 * May 28 1998	Go back to old WCSFULL, computing height and width from center
 * May 29 1998	Add wcskinit() to initialize WCS from arguments
 * May 29 1998	Add wcstype() to set projection from arguments
 * May 29 1998	Add wcscdset(), and wcsdeltset() to set scale from arguments
 * Jun  1 1998	In wcstype(), reconstruct ctype for WCS structure
 * Jun 11 1998	Split off header-dependent subroutines to wcsinit.c
 * Jun 18 1998	Add wcspcset() for PC matrix initialization
 * Jun 24 1998	Add string lengths to ra2str(), dec2str, and deg2str() calls
 * Jun 25 1998	Use AIPS software for CAR projection
 * Jun 25 1998	Add wcsndec to set number of decimal places in output string
 * Jul  6 1998	Add wcszin() and wcszout() to use third dimension of images
 * Jul  7 1998	Change setlinmode() to setwcslin(); setdegout() to setwcsdeg()
 * Jul 10 1998	Initialize matrices correctly for naxis > 2 in wcs<>set()
 * Jul 16 1998	Initialize coordinates to be returned in wcspos()
 * Jul 17 1998	Link lin structure arrays to wcs structure arrays
 * Jul 20 1998	In wcscdset() compute sign of xinc and yinc from CD1_1, CD 2_2
 * Jul 20 1998	In wcscdset() compute sign of rotation based on CD1_1, CD 1_2
 * Jul 22 1998	Add wcslibrot() to compute lin() rotation matrix
 * Jul 30 1998	Set wcs->naxes and lin.naxis in wcsxinit() and wcskinit()
 * Aug  5 1998	Use old WCS subroutines to deal with COE projection (for ESO)
 * Aug 14 1998	Add option to print image coordinates with wcscom()
 * Aug 14 1998	Add multiple command options to wcscom*()
 * Aug 31 1998	Declare undeclared arguments to wcspcset()
 * Sep  3 1998	Set CD rotation and cdelts from sky axis position angles
 * Sep 16 1998	Add option to use North Polar Angle instead of Latitude
 * Sep 29 1998	Initialize additional WCS commands from the environment
 * Oct 14 1998	Fix bug in wcssize() which didn't divide dra by cos(dec)
 * Nov 12 1998	Fix sign of CROTA when either axis is reflected
 * Dec  2 1998	Fix non-arcsecond scale factors in wcscent()
 * Dec  2 1998	Add PLANET coordinate system to pix2wcst()

 * Jan 20 1999	Free lin.imgpix and lin.piximg in wcsfree()
 * Feb 22 1999	Fix bug setting latitude reference value of latbase != 0
 * Feb 22 1999	Fix bug so that quad cube faces are 0-5, not 1-6
 * Mar 16 1999	Always initialize all 4 imgcrds and pixcrds in wcspix()
 * Mar 16 1999	Always return (0,0) from wcs2pix if offscale
 * Apr  7 1999	Add code to put file name in error messages
 * Apr  7 1999	Document utility subroutines at end of file
 * May  6 1999	Fix bug printing height of LINEAR image
 * Jun 16 1999	Add wcsrange() to return image RA and Dec limits
 * Jul  8 1999	Always use FK5 and FK4 instead of J2000 and B1950 in RADECSYS
 * Aug 16 1999	Print dd:mm:ss dd:mm:ss if not J2000 or B1950 output
 * Aug 20 1999	Add WCS string argument to wcscom(); don't compute it there
 * Aug 20 1999	Change F3 WCS command default from Tycho to ACT
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 21 1999	Drop declarations of unused functions after lint
 * Oct 25 1999	Drop out of wcsfree() if wcs is null pointer
 * Nov 17 1999	Fix bug which caused software to miss NCP projection
 *
 * Jan 24 2000	Default to AIPS for NCP, CAR, and COE proj.; if -z use WCSLIB
 * Feb 24 2000	If coorsys is null in wcsc2pix, wcs->radecin is assumed
 * May 10 2000	In wcstype(), default to WCS_LIN, not error (after Bill Joye)
 * Jun 22 2000	In wcsrotset(), leave rotation angle alone in 1-d image
 * Jul  3 2000	Initialize wcscrd[] to zero in wcspix()
 *
 * Feb 20 2001	Add recursion to wcs2pix() and pix2wcs() for dependent WCS's
 * Mar 20 2001	Add braces to avoid ambiguity in if/else groupings
 * Mar 22 2001	Free WCS structure in wcsfree even if it is not filled
 */
