/*** File libwcs/imgetwcs.c
 *** January 11, 2001
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** (remotely based on UIowa code)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fitshead.h"
#include "vimoswcs.h"
#include "lvimoswcs.h"

/* Get the C* WCS fields in  a FITS header based on a reference catalog
 * do it by finding stars in the image and in the reference catalog and
 * finding the rotation and offsets which result in a best-fit.
 * verbose generates extra info on stderr.
 * try using deeper reference star catalog searches if there is trouble.
 * return 1 if all ok, else 0
 */

/* These parameters can be set on the command line */
static double secpix0 = PSCALE;		/* Set image scale--override header */
static double secpix2 = PSCALE;		/* Set image scale 2--override header */
static double *cd0 = NULL;		/* Set CD matrix--override header */
static double rot0 = 361.0;		/* Initial image rotation */
static int comsys = VIMOSWCS_J2000;		/* Command line center coordinte system */
static int wp0 = 0;			/* Initial width of image */
static int hp0 = 0;			/* Initial height of image */
static double ra0 = -99.0;		/* Initial center RA in degrees */
static double dec0 = -99.0;		/* Initial center Dec in degrees */
static double xref0 = -99999.0;		/* Reference pixel X coordinate */
static double yref0 = -99999.0;		/* Reference pixel Y coordinate */
static int ptype0 = -1;			/* Projection type to fit */
static int  nctype = 28;		/* Number of possible projections */
static char ctypes[28][4];		/* 3-letter codes for projections */
static int usecdelt = 0;		/* Use CDELT if 1, else CD matrix */

/* Set a nominal world coordinate system from image header info.
 * If the image center is not FK5 (J2000) equinox, convert it
 * Return a WCS structure if OK, else return NULL
 */

struct WorldCoor *
GetFITSVIMOSWCS (filename, header, verbose, cra, cdec, dra, ddec, secpix, wp, hp,
	    sysout, eqout)

char	*filename;	/* FITS or IRAF file name */
char	*header;	/* Image FITS header */
int	verbose;	/* Extra printing if =1 */
double	*cra;		/* Center right ascension in degrees (returned) */
double	*cdec;		/* Center declination in degrees (returned) */
double	*dra;		/* Right ascension half-width in degrees (returned) */
double	*ddec;		/* Declination half-width in degrees (returned) */
double	*secpix;	/* Arcseconds per pixel (returned) */
int	*wp;		/* Image width in pixels (returned) */
int	*hp;		/* Image height in pixels (returned) */
int	*sysout;	/* Coordinate system to return (0=image, returned) */
double	*eqout;		/* Equinox to return (0=image, returned) */
{
    int nax;
    double eq1, xref, yref, degpix, ra1, dec1;
    struct WorldCoor *vimoswcs;
    char rstr[64], dstr[64], temp[16], cstr[16];

    /* Set image dimensions */
    nax = 0;
    if (hp0 > 0 || wp0 > 0) {
	*hp = hp0;
	*wp = wp0;
	if (*hp > 0 && *wp > 0)
	    nax = 2;
	else
	    nax = 1;
	hputi4 (header, "NAXIS", nax);
	hputi4 (header, "NAXIS1", *wp);
	hputi4 (header, "NAXIS2", *hp);
	}
    else if (hgeti4 (header,"NAXIS",&nax) < 1)
	return (NULL);
    else {
	if (hgeti4 (header,"NAXIS1",wp) < 1)
	    return (NULL);
	else {
	    if (hgeti4 (header,"NAXIS2",hp) < 1)
		return (NULL);
	    }
	}

    /* Set plate center from command line, if it is there */
    if (ra0 > -99.0 && dec0 > -99.0) {
	hputnr8 (header, "CRVAL1" ,8,ra0);
	hputnr8 (header, "CRVAL2" ,8,dec0);
	hputra (header, "RA", ra0);
	hputdec (header, "DEC", dec0);
	if (comsys == VIMOSWCS_B1950) {
	    hputi4 (header, "EPOCH", 1950);
	    hputi4 (header, "EQUINOX", 1950);
	    hputs (header, "RADECSYS", "FK4");
	    }
	else {
	    hputi4 (header, "EPOCH", 2000);
	    hputi4 (header, "EQUINOX", 2000);
	    if (comsys == VIMOSWCS_GALACTIC)
		hputs (header, "RADECSYS", "GALACTIC");
	    else if (comsys == VIMOSWCS_ECLIPTIC)
		hputs (header, "RADECSYS", "ECLIPTIC");
	    else
		hputs (header, "RADECSYS", "FK5");
	    }
	if (hgetr8 (header, "SECPIX", secpix)) {
	    degpix = *secpix / 3600.0;
	    hputnr8 (header, "CDELT1", 8, -degpix);
	    hputnr8 (header, "CDELT2", 8, degpix);
	    hdel (header, "CD1_1");
	    hdel (header, "CD1_2");
	    hdel (header, "CD2_1");
	    hdel (header, "CD2_2");
	    } 
	}
    if (ptype0 > -1 && ptype0 < nctype) {
	strcpy (temp,"RA---");
	strcat (temp, ctypes[ptype0]);
	hputc (header, "CTYPE1", temp);
	strcpy (temp,"DEC--");
	strcat (temp, ctypes[ptype0]);
	hputc (header, "CTYPE2", temp);
	}

    /* Set reference pixel from command line, if it is there */
    if (xref0 > -99999.0 && yref0 > -99999.0) {
	hputr8 (header, "CRPIX1", xref0);
	hputr8 (header, "CRPIX2", yref0);
	}
    else if (hgetr8 (header, "CRPIX1", &xref) < 1) {
	xref = (double) *wp / 2.0;
	yref = (double) *hp / 2.0;
	hputnr8 (header, "CRPIX1", 3, xref);
	hputnr8 (header, "CRPIX2", 3, yref);
	}

    /* Set plate scale from command line, if it is there */
    if (secpix0 != 0.0 || cd0 != NULL) {
        if (secpix2 != 0.0) {
	    *secpix = 0.5 * (secpix0 + secpix2);
	    hputnr8 (header, "SECPIX1", 5, secpix0);
	    hputnr8 (header, "SECPIX2", 5, secpix2);
	    degpix = -secpix0 / 3600.0;
	    hputnr8 (header, "CDELT1", 8, degpix);
	    degpix = secpix2 / 3600.0;
	    hputnr8 (header, "CDELT2", 8, degpix);
	    hdel (header, "CD1_1");
	    hdel (header, "CD1_2");
	    hdel (header, "CD2_1");
	    hdel (header, "CD2_2");
	    }
	else if (secpix0 != 0.0) {
	    *secpix = secpix0;
	    hputnr8 (header, "SECPIX", 5, *secpix);
	    degpix = *secpix / 3600.0;
	    hputnr8 (header, "CDELT1", 8, -degpix);
	    hputnr8 (header, "CDELT2", 8, degpix);
	    hdel (header, "CD1_1");
	    hdel (header, "CD1_2");
	    hdel (header, "CD2_1");
	    hdel (header, "CD2_2");
	    }
	else {
	    hputr8 (header, "CD1_1", cd0[0]);
	    hputr8 (header, "CD1_2", cd0[1]);
	    hputr8 (header, "CD2_1", cd0[2]);
	    hputr8 (header, "CD2_2", cd0[3]);
	    hdel (header, "CDELT1");
	    hdel (header, "CDELT2");
	    hdel (header, "CROTA1");
	    hdel (header, "CROTA2");
	    }
	if (!ksearch (header,"CRVAL1")) {
	    hgetra (header, "RA", &ra0);
	    hgetdec (header, "DEC", &dec0);
	    hputnr8 (header, "CRVAL1", 8, ra0);
	    hputnr8 (header, "CRVAL2", 8, dec0);
	    }
	if (!ksearch (header,"CRPIX1")) {
	    xref = (double) *wp / 2.0;
	    yref = (double) *hp / 2.0;
	    hputnr8 (header, "CRPIX1", 3, xref);
	    hputnr8 (header, "CRPIX2", 3, yref);
	    }
	if (!ksearch (header,"CTYPE1")) {
	    if (comsys == VIMOSWCS_GALACTIC) {
		hputc (header, "CTYPE1", "GLON-TAN");
		hputc (header, "CTYPE2", "GLAT-TAN");
		}
	    else {
		hputc (header, "CTYPE1", "RA---TAN");
		hputc (header, "CTYPE2", "DEC--TAN");
		}
	    }
	}

    /* Set rotation angle from command line, if it is there */
    if (rot0 < 361.0) {
	hputnr8 (header, "CROTA1", 5, rot0);
	hputnr8 (header, "CROTA2", 5, rot0);
	}

    /* Initialize WCS structure from FITS header */
    vimoswcs = vimoswcsinit (header);

    /* If incomplete WCS in header, drop out */
    if (novimoswcs (vimoswcs)) {
	setvimoswcsfile (filename);
	vimoswcserr();
	if (verbose)
	    fprintf (stderr,"Insufficient information for initial WCS\n");
	return (NULL);
	}

    /* Set flag to get appropriate equinox for catalog search */
    if (!*sysout)
	*sysout = vimoswcs->sysvimoswcs;
    if (*eqout == 0.0)
	*eqout = vimoswcs->equinox;
    eq1 = vimoswcs->equinox;
    if (vimoswcs->coorflip) {
	ra1 = vimoswcs->crval[1];
	dec1 = vimoswcs->crval[0];
	}
    else {
	ra1 = vimoswcs->crval[0];
	dec1 = vimoswcs->crval[1];
	}

    /* Print reference pixel position and value */
    if (verbose && (eq1 != *eqout || vimoswcs->sysvimoswcs != *sysout)) {
	if (vimoswcs->degout) {
	    deg2str (rstr, 32, ra1, 6);
            deg2str (dstr, 32, dec1, 6);
	    }
	else {
	    ra2str (rstr, 32, ra1, 3);
            dec2str (dstr, 32, dec1, 2);
	    }
	vimoswcscstr (cstr, vimoswcs->sysvimoswcs, vimoswcs->equinox, vimoswcs->epoch);
	fprintf (stderr,"Reference pixel (%.2f,%.2f) %s %s %s\n",
		 vimoswcs->xrefpix, vimoswcs->yrefpix, rstr, dstr, cstr);
	}

    /* Get center and size for catalog searching */
    vimoswcssize (vimoswcs, cra, cdec, dra, ddec);

    /* Convert center to desired output coordinate system */
    vimoswcscon (vimoswcs->sysvimoswcs, *sysout, vimoswcs->equinox, *eqout, cra, cdec, vimoswcs->epoch);

    /* Compute plate scale to return if it was not set on the command line */
    if (secpix0 <= 0.0)
	*secpix = 3600.0 * 2.0 * *ddec / (double) *hp;

    /* Set reference pixel to center of image if it has not been set */
    if (vimoswcs->xref == 0.0 && vimoswcs->yref == 0.0) {
	vimoswcs->xref = *cra;
	vimoswcs->yref = *cdec;
	ra1 = *cra;
	dec1 = *cdec;
	if (vimoswcs->xrefpix == 0.0 && vimoswcs->yrefpix == 0.0) {
	    vimoswcs->xrefpix = (double) vimoswcs->nxpix * 0.5;
	    vimoswcs->yrefpix = (double) vimoswcs->nypix * 0.5;
	    }
	vimoswcs->xinc = *dra * 2.0 / (double) vimoswcs->nxpix;
	vimoswcs->yinc = *ddec * 2.0 / (double) vimoswcs->nypix;
	/* hchange (header,"PLTRAH","PLT0RAH");
	vimoswcs->plate_fit = 0; */
	}

    /* Convert center to desired coordinate system */
    else if (vimoswcs->sysvimoswcs != *sysout && vimoswcs->equinox != *eqout) {
	vimoswcscon (vimoswcs->sysvimoswcs, *sysout, vimoswcs->equinox, *eqout, &ra1, &dec1, vimoswcs->epoch);
	if (vimoswcs->coorflip) {
	    vimoswcs->yref = ra1;
	    vimoswcs->xref = dec1;
	    }
	else {
	    vimoswcs->xref = ra1;
	    vimoswcs->yref = dec1;
	    }
	}
    vimoswcs->crval[0] = vimoswcs->xref;
    vimoswcs->crval[1] = vimoswcs->yref;
    vimoswcs->cel.ref[0] = vimoswcs->crval[0];
    vimoswcs->cel.ref[1] = vimoswcs->crval[1];
    vimoswcs->cel.flag = 0;
    vimoswcs->vimoswcsl.flag = 0;
    vimoswcs->equinox = *eqout;
    vimoswcs->sysvimoswcs = *sysout;
    vimoswcs->sysout = *sysout;
    vimoswcs->eqout = *eqout;
    vimoswcs->sysin = *sysout;
    vimoswcs->eqin = *eqout;
    vimoswcscstr (cstr,*sysout,*eqout,vimoswcs->epoch);
    strcpy (vimoswcs->radecsys, cstr);
    strcpy (vimoswcs->radecout, cstr);
    strcpy (vimoswcs->radecin, cstr);
    vimoswcsininit (vimoswcs, vimoswcs->radecsys);
    vimoswcsoutinit (vimoswcs, vimoswcs->radecsys);

    if (usecdelt) {
	hputnr8 (header, "CDELT1", 9, vimoswcs->xinc);
	if (vimoswcs->naxes > 1) {
	    hputnr8 (header, "CDELT2", 9, vimoswcs->yinc);
	    hputnr8 (header, "CROTA2", 9, vimoswcs->rot);
	    }
	hdel (header, "CD1_1");
	hdel (header, "CD1_2");
	hdel (header, "CD2_1");
	hdel (header, "CD2_2");
	}
    else {
	hputnr8 (header, "CD1_1", 9, vimoswcs->cd[0]);
	if (vimoswcs->naxes > 1) {
	    hputnr8 (header, "CD1_2", 9, vimoswcs->cd[1]);
	    hputnr8 (header, "CD2_1", 9, vimoswcs->cd[2]);
	    hputnr8 (header, "CD2_2", 9, vimoswcs->cd[3]);
	    }
	}

    /* Print reference pixel position and value */
    if (verbose) {
	if (vimoswcs->degout) {
	    deg2str (rstr, 32, ra1, 6);
            deg2str (dstr, 32, dec1, 6);
	    }
	else {
	    ra2str (rstr, 32, ra1, 3);
            dec2str (dstr, 32, dec1, 2);
	    }
	vimoswcscstr (cstr,*sysout,*eqout,vimoswcs->epoch);
	fprintf (stderr,"Reference pixel (%.2f,%.2f) %s %s %s\n",
		vimoswcs->xrefpix, vimoswcs->yrefpix, rstr, dstr, cstr);
	}

    /* Image size for catalog search */
    if (verbose) {
	if (vimoswcs->degout) {
	    deg2str (rstr, 32, *cra, 6);
            deg2str (dstr, 32, *cdec, 6);
	    }
	else {
	    ra2str (rstr, 32, *cra, 3);
	    dec2str (dstr, 32, *cdec, 2);
	    }
	vimoswcscstr (cstr, *sysout, *eqout, vimoswcs->epoch);
	fprintf (stderr,"Search at %s %s %s", rstr, dstr, cstr);
	fprintf (stderr," +- %s %s\n", rstr, dstr);
	fprintf (stderr,"Image width=%d height=%d, %g arcsec/pixel\n",
				*wp, *hp, *secpix);
	}

    return (vimoswcs);
}

void
setcdelt()			/* Set flag to use CDELTn, not CD matrix */
{usecdelt = 1; return;}

void
setnpix (	/* Set image size */
int nx,
int ny)
{ wp0 = nx; hp0 = ny; return; }

void
setrot (
double rot)
{ rot0 = rot; return; }

void
setsecpix (		/* Set first axis arcseconds per pixel */
double secpix)
{ secpix0 = secpix; return; }

void
setsecpix2 (		/* Set second axis arcseconds per pixel */
double secpix)
{ secpix2 = secpix; return; }

void
setcd (		/* Set initial CD matrix */
double *cd)
{ int i;
  if (cd0 != NULL) free (cd0);
  cd0 = (double *) calloc (4, sizeof (double));
  for (i = 0; i < 4; i++) cd0[i] = cd[i];
  return; }

void
setsys (		/* Set WCS coordinates as FK4 */
int comsys0)
{ comsys = comsys0; return; }

void
setcenter (	/* Set center sky coordinates in strings */
char *rastr,
char *decstr)
{
    ra0 = str2ra (rastr);
    dec0 = str2dec (decstr);
    return;
}

void
setdcenter (		/* Set center sky coordinates in degrees */
double ra, 
double dec)
{
    ra0 = ra;
    dec0 = dec;
    return;
}

void
setrefpix (
double x, 
double y)
{ xref0 = x; yref0 = y; return; }

void
setproj (
char*	ptype)
{
    int i;

    /* Set up array of projection types */
    strcpy (ctypes[0], "DSS");
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
    strcpy (ctypes[27], "TNX");

    ptype0 = -1;
    for (i = 0; i < nctype; i++) {
	if (!strcmp(ptype, ctypes[i]))
	    ptype0 = i;
	}
    return;
}


/* Feb 29 1996	New program
 * May 23 1996	Use pre-existing WCS for center, if it is present
 * May 29 1996	Simplify program by always using WCS structure
 * Jun 12 1996	Be more careful with nominal WCS setting
 * Jul  3 1996	Set epoch from old equinox if not already set
 * Jul 19 1996	Set image center in WCS if DSS WCS
 * Aug  5 1996	Check for SECPIX1 as well as SECPIX
 * Aug  7 1996	Save specified number of decimal places in header parameters
 * Aug  7 1996	Rename old center parameters
 * Aug 26 1996	Decrease default pixel tolerance from 20 to 10
 * Sep  1 1996	Set plate scale default in lvimoswcs.h
 * Sep  3 1996	Fix bug to set plate scale from command line
 * Oct 15 1996	Break off from imsetvimoswcs.c
 * Oct 16 1996	Clean up center setting so eqref is used
 * Oct 17 1996	Do not print error messages unless verbose is set
 * Oct 30 1996	Keep equinox from image if EQREF is zero
 * Nov  1 1996	Declare undeclared subroutines; remove unused variables
 * Nov  4 1996	Add reference pixel and projection to vimoswcsset() call
 * Nov 14 1996	Add GetLimits() to deal with search limits around the poles
 * Nov 15 1996	Drop GetLimits(); code moved to individual catalog routines
 * Dec 10 1996	Fix precession and make equinox double

 * Feb 19 1997	If eqout is 0, use equinox of image
 * Feb 24 1997	Always convert center to output equinox (bug fix)
 * Mar 20 1997	Declare EQ2 double instead of int, fixing a bug
 * Jul 11 1997	Allow external (command line) setting of reference pixel coords
 * Sep 26 1997	Set both equinox and epoch if input center coordinates
 * Nov  3 1997	Separate WCS reference pixel from search center
 * Dec  8 1997	Set CDELTn using SECPIX if it is in the header
 *
 * Jan  6 1998	Do not print anything unless verbose is set
 * Jan 29 1998	Use flag to allow AIPS classic WCS subroutines
 * Mar  1 1998	Set x and y plate scales from command line if there are two
 * Mar  2 1998	Do not reset plate solution switch
 * Mar  6 1998	Add option to reset center sky coordinates in degrees
 * Mar  6 1998	Add option to set projection type
 * Mar 18 1998	Initialize reference pixel if CDELTn's being set
 * Mar 20 1998	Only set CTYPEn if CDELTn or CRVALn are set
 * Apr 10 1998	Add option to set search area to circle or box
 * Apr 20 1998	Move GetArea() to scat.c
 * Apr 24 1998	Always convert image reference coordinate to catalog equinox
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * Jun  1 1998	Print error message if WCS cannot be initialized
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jun 25 1998	Leave use of AIPS wcs to vimoswcs.c file
 * Sep 17 1998	Add sysout to argument list and use scscon() for conversion
 * Sep 25 1998	Make sysout==0 indicate output in image coordinate system
 * Oct 28 1998	Set coordinate system properly to sysout/eqout in GetFITSVIMOSWCS()
 *
 * Apr  7 1999	Add filename argument to GetFITSVIMOSWCS
 * Apr 29 1999	Add option to set image size
 * Jun  2 1999	Fix sign of CDELT1 if secpix2 and secpix0 are set
 * Jul  7 1999	Fix conversion of center coordinates to refsys
 * Jul  9 1999	Fix bug which reset command-line-set reference pixel coordinate
 * Oct 21 1999	Fix declarations after lint
 * Nov  1 1999	Add option to write CD matrix
 * Nov  1 1999	If CDELTn set from command line delete previous header CD matrix
 * Nov 12 1999	Add galactic coordinates as command line option
 * Nov 16 1999	Set radecsys correctly for command line galactic
 *
 * Feb 15 2000	Add option to override the header CD matrix (like CDELTs)
 * Feb 29 2000	Fix bug, converting reference pixel WCS coordinates everywhere
 * Mar 27 2000	Drop unused subroutine setradius()
 * May 24 2000	Print degrees in debugging messages if output format
 * Sep 14 2000	Set xinc and yinc correctly if center pixel in header
 *
 * Jan 11 2001	All printing to stderr
 */
