/*** File libwcs/wcsinit.c
 *** October 11, 2000
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics

 * Module:	wcsinit.c (World Coordinate Systems)
 * Purpose:	Convert FITS WCS to pixels and vice versa:
 * Subroutine:	wcsinit (hstring) sets a WCS structure from an image header
 * Subroutine:	wcsninit (hstring,lh) sets a WCS structure from an image header
 * Subroutine:	wcseq (hstring, wcs) set radecsys and equinox from image header

 * Copyright:   2000 Smithsonian Astrophysical Observatory
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

static void vimoswcseq();
void vimoswcsrotset();
char vimoswcserrmsg[80];
extern int dsspos();

/* set up a WCS structure from a FITS image header lhstring bytes long */

struct WorldCoor *
vimoswcsninit (hstring, lhstring)

char	*hstring;	/* character string containing FITS header information
		   	in the format <keyword>= <value> [/ <comment>] */
int	lhstring;	/* Length of FITS header in bytes */
{
    hlength (hstring, lhstring);
    return (vimoswcsinit (hstring));
}

/* set up a WCS structure from a FITS image header */

struct WorldCoor *
vimoswcsinit (hstring)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> [/ <comment>] */
{
    struct WorldCoor *vimoswcs;
    char ctype1[32], ctype2[32];
    char *hcoeff;		/* pointer to first coeff's in header */
    char decsign;
    double rah,ram,ras, dsign,decd,decm,decs;
    double dec_deg,ra_hours, secpix, ra0, ra1, dec0, dec1;
    double cdelt1, cdelt2, cd[4], pc[16];
    char keyword[16];
    int ieq, i, naxes, cd11p, cd12p, cd21p, cd22p;
    /* int ix1, ix2, iy1, iy2, idx1, idx2, idy1, idy2;
    double dxrefpix, dyrefpix;
    char temp[48];
    char *ic; */
    double mjd;
    double rot;
    extern int tnxinit();
    extern int platepos();

    vimoswcs = (struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

    /* Set WCSLIB flags so that structures will be reinitialized */
    vimoswcs->cel.flag = 0;
    vimoswcs->lin.flag = 0;
    vimoswcs->vimoswcsl.flag = 0;

    /* Initialize to no plate fit */
    vimoswcs->ncoeff1 = 0;
    vimoswcs->ncoeff2 = 0;

    /* Initialize to no CD matrix */
    cdelt1 = 0.0;
    cdelt2 = 0.0;
    cd[0] = 0.0;
    cd[1] = 0.0;
    cd[2] = 0.0;
    cd[3] = 0.0;
    pc[0] = 0.0;
    vimoswcs->rotmat = 0;
    vimoswcs->rot = 0.0;

    /* Header parameters independent of projection */
    naxes = 0;
    hgeti4 (hstring, "NAXIS", &naxes);
    if (naxes < 1) {
	setvimoswcserr ("WCSINIT: No NAXIS keyword");
	vimoswcsfree (vimoswcs);
	return (NULL);
	}
    vimoswcs->naxes = naxes;
    vimoswcs->lin.naxis = naxes;
    vimoswcs->nxpix = 0;
    hgetr8 (hstring, "NAXIS1", &vimoswcs->nxpix);
    if (vimoswcs->nxpix < 1) {
	setvimoswcserr ("WCSINIT: No NAXIS1 keyword");
	vimoswcsfree (vimoswcs);
	return (NULL);
	}
    vimoswcs->nypix = 0;
    hgetr8 (hstring, "NAXIS2", &vimoswcs->nypix);
    hgets (hstring, "INSTRUME", 16, vimoswcs->instrument);
    hgeti4 (hstring, "DETECTOR", &vimoswcs->detector);
    vimoswcs->vimoswcsproj = getdefvimoswcs();
    for (i = 0; i < 16; i++) vimoswcs->pc[i] = 0.0;
    for (i = 0; i < naxes; i++) vimoswcs->pc[(i*naxes)+i] = 1.0;

    /* World coordinate system reference coordinate information */
    if (hgets (hstring,"CTYPE1", 16, ctype1)) {

	/* Read second coordinate type */
	if (!hgets (hstring,"CTYPE2", 16, ctype2)) {
	    setvimoswcserr ("WCSINIT: No CTYPE2 -> no WCS");
	    vimoswcsfree (vimoswcs);
	    return (NULL);
	    }

	/* Read third and fourth coordinate types, if present */
	strcpy (vimoswcs->ctype[0], ctype1);
	strcpy (vimoswcs->ctype[1], ctype2);
	strcpy (vimoswcs->ctype[2], "");
	hgets (hstring, "CTYPE3", 9, vimoswcs->ctype[2]);
	strcpy (vimoswcs->ctype[3], "");
	hgets (hstring, "CTYPE4", 9, vimoswcs->ctype[3]);

	/* Set projection type in WCS data structure */
	if (vimoswcstype (vimoswcs, ctype1, ctype2)) {
	    vimoswcsfree (vimoswcs);
	    return (NULL);
	    }

	/* Get units, if present, for linear coordinates */
	if (vimoswcs->prjcode == VIMOSWCS_LIN) {
	    if (!hgets (hstring, "CUNIT1", 16, vimoswcs->units[0])) {
		if (!mgets (hstring, "WAT1", "units", 16, vimoswcs->units[0])) {
		    vimoswcs->units[0][0] = 0;
		    }
		}
	    if (!hgets (hstring, "CUNIT2", 16, vimoswcs->units[1])) {
		if (!mgets (hstring, "WAT2", "units", 16, vimoswcs->units[1])) {
		    vimoswcs->units[1][0] = 0;
		    }
		}
	    }

	/* Reference pixel coordinates and WCS value */
	vimoswcs->crpix[0] = 1.0;
	hgetr8 (hstring,"CRPIX1",&vimoswcs->crpix[0]);
	vimoswcs->crpix[1] = 1.0;
	hgetr8 (hstring,"CRPIX2",&vimoswcs->crpix[1]);
	vimoswcs->xrefpix = vimoswcs->crpix[0];
	vimoswcs->yrefpix = vimoswcs->crpix[1];
	vimoswcs->crval[0] = 0.0;
	hgetr8 (hstring,"CRVAL1",&vimoswcs->crval[0]);
	vimoswcs->crval[1] = 0.0;
	hgetr8 (hstring,"CRVAL2",&vimoswcs->crval[1]);
	if (vimoswcs->sysvimoswcs == VIMOSWCS_NPOLE)
	    vimoswcs->crval[1] = 90.0 - vimoswcs->crval[1];
	if (vimoswcs->sysvimoswcs == VIMOSWCS_SPA)
	    vimoswcs->crval[1] = vimoswcs->crval[1] - 90.0;
	vimoswcs->xref = vimoswcs->crval[0];
	vimoswcs->yref = vimoswcs->crval[1];
	if (vimoswcs->coorflip) {
	    vimoswcs->cel.ref[0] = vimoswcs->crval[1];
	    vimoswcs->cel.ref[1] = vimoswcs->crval[0];
	    }
	else {
	    vimoswcs->cel.ref[0] = vimoswcs->crval[0];
	    vimoswcs->cel.ref[1] = vimoswcs->crval[1];
	    }
	vimoswcs->longpole = 999.0;
	hgetr8 (hstring,"LONGPOLE",&vimoswcs->longpole);
	vimoswcs->cel.ref[2] = vimoswcs->longpole;
	vimoswcs->latpole = 999.0;
	hgetr8 (hstring,"LATPOLE",&vimoswcs->latpole);
	vimoswcs->cel.ref[3] = vimoswcs->latpole;
	vimoswcs->lin.crpix = vimoswcs->crpix;
	vimoswcs->lin.cdelt = vimoswcs->cdelt;
	vimoswcs->lin.pc = vimoswcs->pc;

	/* Projection constants */
	vimoswcs->prj.r0 = 0.0;
	hgetr8 (hstring, "PROJR0", &vimoswcs->prj.r0);
	for (i = 1; i < 10; i++) {
	    vimoswcs->prj.p[i] = 0.0;
	    sprintf (keyword,"PROJP%d",i);
	    hgetr8 (hstring, keyword, &vimoswcs->prj.p[i]);
	    }

	/* Use polynomial fit instead of projection, if present */
	vimoswcs->ncoeff1 = 0;
	vimoswcs->ncoeff2 = 0;
	cd11p = hgetr8 (hstring,"CD1_1",&cd[0]);
	cd12p = hgetr8 (hstring,"CD1_2",&cd[1]);
	cd21p = hgetr8 (hstring,"CD2_1",&cd[2]);
	cd22p = hgetr8 (hstring,"CD2_2",&cd[3]);
	if (vimoswcs->vimoswcsproj != VIMOSWCS_OLD && (hcoeff = ksearch (hstring,"CO1_1")) != NULL) {
	    vimoswcs->prjcode = VIMOSWCS_PLT;
	    (void)strcpy (vimoswcs->ptype, "PLATE");
	    for (i = 0; i < 20; i++) {
		sprintf (keyword,"CO1_%d", i+1);
		vimoswcs->x_coeff[i] = 0.0;
		if (hgetr8 (hcoeff, keyword, &vimoswcs->x_coeff[i]))
		    vimoswcs->ncoeff1 = i + 1;
		}
	    hcoeff = ksearch (hstring,"CO2_1");
	    for (i = 0; i < 20; i++) {
		sprintf (keyword,"CO2_%d",i+1);
		vimoswcs->y_coeff[i] = 0.0;
		if (hgetr8 (hcoeff, keyword, &vimoswcs->y_coeff[i]))
		    vimoswcs->ncoeff2 = i + 1;
		}

	    /* Compute a nominal scale factor */
	    platepos (vimoswcs->crpix[0], vimoswcs->crpix[1], vimoswcs, &ra0, &dec0);
	    platepos (vimoswcs->crpix[0], vimoswcs->crpix[1]+1.0, vimoswcs, &ra1, &dec1);
	    vimoswcs->yinc = dec1 - dec0;
	    vimoswcs->xinc = -vimoswcs->yinc;

	    /* Compute image rotation angle */
	    vimoswcs->vimoswcson = 1;
	    vimoswcsrotset (vimoswcs);
	    rot = degrad (vimoswcs->rot);

	    /* Compute scale at reference pixel */
	    platepos (vimoswcs->crpix[0], vimoswcs->crpix[1], vimoswcs, &ra0, &dec0);
	    platepos (vimoswcs->crpix[0]+cos(rot),
		      vimoswcs->crpix[1]+sin(rot), vimoswcs, &ra1, &dec1);
	    vimoswcs->cdelt[0] = -vimoswcsdist (ra0, dec0, ra1, dec1);
	    vimoswcs->xinc = vimoswcs->cdelt[0];
	    platepos (vimoswcs->crpix[0]+sin(rot),
		      vimoswcs->crpix[1]+cos(rot), vimoswcs, &ra1, &dec1);
	    vimoswcs->cdelt[1] = vimoswcsdist (ra0, dec0, ra1, dec1);
	    vimoswcs->yinc = vimoswcs->cdelt[1];

	    /* Set CD matrix from header */
	    vimoswcs->cd[0] = cd[0];
	    vimoswcs->cd[1] = cd[1];
	    vimoswcs->cd[2] = cd[2];
	    vimoswcs->cd[3] = cd[3];
	    (void) vimosmatinv (2, vimoswcs->cd, vimoswcs->dc);
	    }

	/* Else use CD matrix, if and elements are present */
	else if (cd11p || cd12p || cd21p || cd22p) {
	    vimoswcs->rotmat = 1;
	    vimoswcscdset (vimoswcs, cd);
	    }

	/* Else get scaling from CDELT1 and CDELT2 */
	else if (hgetr8 (hstring,"CDELT1",&cdelt1) != 0) {
	    hgetr8 (hstring,"CDELT2",&cdelt2);

	    /* If CDELT1 or CDELT2 is 0 or missing */
	    if (cdelt1 == 0.0 || (vimoswcs->nypix > 1 && cdelt2 == 0.0)) {
		if (ksearch (hstring,"SECPIX") != NULL ||
		    ksearch (hstring,"PIXSCALE") != NULL ||
		    ksearch (hstring,"PIXSCAL1") != NULL ||
		    ksearch (hstring,"SECPIX1") != NULL) {
		    secpix = 0.0;
		    hgetr8 (hstring,"SECPIX",&secpix);
		    if (secpix == 0.0)
			hgetr8 (hstring,"PIXSCALE",&secpix);
		    if (secpix == 0.0) {
			hgetr8 (hstring,"SECPIX1",&secpix);
			if (secpix != 0.0) {
			    if (cdelt1 == 0.0)
				cdelt1 = -secpix / 3600.0;
			    if (cdelt2 == 0.0) {
				hgetr8 (hstring,"SECPIX2",&secpix);
				cdelt2 = secpix / 3600.0;
				}
			    }
			else {
			    if (cdelt1 == 0.0) {
				hgetr8 (hstring,"PIXSCAL1",&secpix);
				cdelt1 = -secpix / 3600.0;
				}
			    if (cdelt2 == 0.0) {
				hgetr8 (hstring,"PIXSCAL2",&secpix);
				cdelt2 = secpix / 3600.0;
				}
			    }
			}
		    else {
			if (cdelt1 == 0.0)
			    cdelt1 = -secpix / 3600.0;
			if (cdelt2 == 0.0)
			    cdelt2 = secpix / 3600.0;
			}
		    }
		}
	    if (cdelt2 == 0.0 && vimoswcs->nypix > 1)
		cdelt2 = -cdelt1;
	    vimoswcs->cdelt[2] = 1.0;
	    vimoswcs->cdelt[3] = 1.0;

	    /* Use rotation matrix, if present */
	    for (i = 0; i < 16; i++)
		vimoswcs->pc[i] = 0.0;
	    if (hgetr8 (hstring,"PC001001",&pc[0]) != 0) {
		hgetr8 (hstring,"PC001002",&pc[1]);
		if (naxes < 3) {
		    hgetr8 (hstring,"PC002001",&pc[2]);
		    pc[3] = vimoswcs->pc[0];
		    hgetr8 (hstring,"PC002002",&pc[3]);
		    }
		if (naxes == 3) {
		    hgetr8 (hstring,"PC001003",&pc[2]);
		    hgetr8 (hstring,"PC002001",&pc[3]);
		    pc[4] = vimoswcs->pc[0];
		    hgetr8 (hstring,"PC002002",&pc[4]);
		    hgetr8 (hstring,"PC002003",&pc[5]);
		    hgetr8 (hstring,"PC003001",&pc[6]);
		    hgetr8 (hstring,"PC003002",&pc[7]);
		    pc[8] = 1.0;
		    hgetr8 (hstring,"PC003003",&pc[8]);
		    }
		if (naxes > 3) {
		    hgetr8 (hstring,"PC001003",&pc[2]);
		    hgetr8 (hstring,"PC001004",&pc[3]);
		    hgetr8 (hstring,"PC002001",&pc[4]);
		    pc[5] = vimoswcs->pc[0];
		    hgetr8 (hstring,"PC002002",&pc[5]);
		    hgetr8 (hstring,"PC002003",&pc[6]);
		    hgetr8 (hstring,"PC002004",&pc[7]);
		    hgetr8 (hstring,"PC003001",&pc[8]);
		    hgetr8 (hstring,"PC003002",&pc[9]);
		    pc[10] = 1.0;
		    hgetr8 (hstring,"PC003003",&pc[10]);
		    hgetr8 (hstring,"PC003004",&pc[11]);
		    hgetr8 (hstring,"PC004001",&pc[12]);
		    hgetr8 (hstring,"PC004002",&pc[13]);
		    hgetr8 (hstring,"PC004003",&pc[14]);
		    pc[15] = 1.0;
		    hgetr8 (hstring,"PC004004",&pc[15]);
		    }
		vimoswcspcset (vimoswcs, cdelt1, cdelt2, pc);
		}

	    /* Otherwise, use CROTAn */
	    else {
		rot = 0.0;
		hgetr8 (hstring,"CROTA2",&rot);
		if (rot == 0.)
		    hgetr8 (hstring,"CROTA1",&rot);
		vimoswcsdeltset (vimoswcs, cdelt1, cdelt2, rot);
		}
	    }

	/* If no scaling is present, set to 1 per pixel, no rotation */
	else {
	    vimoswcs->xinc = 1.0;
	    vimoswcs->yinc = 1.0;
	    vimoswcs->cdelt[0] = 1.0;
	    vimoswcs->cdelt[1] = 1.0;
	    vimoswcs->rot = 0.0;
	    vimoswcs->rotmat = 0;
	    setvimoswcserr ("WCSINIT: setting CDELT to 1");
	    }

	/* Initialize TNX, defaulting to TAN if there is a problem */
	if (vimoswcs->prjcode == VIMOSWCS_TNX) {
	    if (tnxinit (hstring, vimoswcs)) {
		vimoswcs->ctype[1][6] = 'A';
		vimoswcs->ctype[1][7] = 'N';
		vimoswcs->prjcode = VIMOSWCS_TAN;
		}
	    }

	/* Coordinate reference frame, equinox, and epoch */
	if (strncmp (vimoswcs->ptype,"LINEAR",6) &&
	    strncmp (vimoswcs->ptype,"PIXEL",5))
	    vimoswcseq (hstring,vimoswcs);
	else {
	    vimoswcs->degout = -1;
	    vimoswcs->ndec = 5;
	    }

	vimoswcs->vimoswcson = 1;
	}

    /* Plate solution coefficients */
    else if (ksearch (hstring,"PLTRAH") != NULL) {
	vimoswcs->prjcode = VIMOSWCS_DSS;
	hcoeff = ksearch (hstring,"PLTRAH");
	hgetr8 (hcoeff,"PLTRAH",&rah);
	hgetr8 (hcoeff,"PLTRAM",&ram);
	hgetr8 (hcoeff,"PLTRAS",&ras);
	ra_hours = rah + (ram / (double)60.0) + (ras / (double)3600.0);
	vimoswcs->plate_ra = hrrad (ra_hours);
	decsign = '+';
	hgets (hcoeff,"PLTDECSN", 1, &decsign);
	if (decsign == '-')
	    dsign = -1.;
	else
	    dsign = 1.;
	hgetr8 (hcoeff,"PLTDECD",&decd);
	hgetr8 (hcoeff,"PLTDECM",&decm);
	hgetr8 (hcoeff,"PLTDECS",&decs);
	dec_deg = dsign * (decd+(decm/(double)60.0)+(decs/(double)3600.0));
	vimoswcs->plate_dec = degrad (dec_deg);
	hgetr8 (hstring,"EQUINOX",&vimoswcs->equinox);
	hgeti4 (hstring,"EQUINOX",&ieq);
	if (ieq == 1950)
	    strcpy (vimoswcs->radecsys,"FK4");
	else
	    strcpy (vimoswcs->radecsys,"FK5");
	vimoswcs->epoch = vimoswcs->equinox;
	hgetr8 (hstring,"EPOCH",&vimoswcs->epoch);
	(void)sprintf (vimoswcs->center,"%2.0f:%2.0f:%5.3f %c%2.0f:%2.0f:%5.3f %s",
		       rah,ram,ras,decsign,decd,decm,decs,vimoswcs->radecsys);
	hgetr8 (hstring,"PLTSCALE",&vimoswcs->plate_scale);
	hgetr8 (hstring,"XPIXELSZ",&vimoswcs->x_pixel_size);
	hgetr8 (hstring,"YPIXELSZ",&vimoswcs->y_pixel_size);
	hgetr8 (hstring,"CNPIX1",&vimoswcs->x_pixel_offset);
	hgetr8 (hstring,"CNPIX2",&vimoswcs->y_pixel_offset);
	hcoeff = ksearch (hstring,"PPO1");
	for (i = 0; i < 6; i++) {
	    sprintf (keyword,"PPO%d", i+1);
	    vimoswcs->ppo_coeff[i] = 0.0;
	    hgetr8 (hcoeff,keyword,&vimoswcs->ppo_coeff[i]);
	    }
	hcoeff = ksearch (hstring,"AMDX1");
	for (i = 0; i < 20; i++) {
	    sprintf (keyword,"AMDX%d", i+1);
	    vimoswcs->x_coeff[i] = 0.0;
	    hgetr8 (hcoeff, keyword, &vimoswcs->x_coeff[i]);
	    }
	hcoeff = ksearch (hstring,"AMDY1");
	for (i = 0; i < 20; i++) {
	    sprintf (keyword,"AMDY%d",i+1);
	    vimoswcs->y_coeff[i] = 0.0;
	    hgetr8 (hcoeff, keyword, &vimoswcs->y_coeff[i]);
	    }
	vimoswcs->vimoswcson = 1;
	(void)strcpy (vimoswcs->c1type, "RA");
	(void)strcpy (vimoswcs->c2type, "DEC");
	(void)strcpy (vimoswcs->ptype, "DSS");
	vimoswcs->degout = 0;
	vimoswcs->ndec = 3;

	/* Compute a nominal reference pixel at the image center */
	strcpy (vimoswcs->ctype[0], "RA---DSS");
	strcpy (vimoswcs->ctype[1], "DEC--DSS");
	vimoswcs->crpix[0] = 0.5 * vimoswcs->nxpix;
	vimoswcs->crpix[1] = 0.5 * vimoswcs->nypix;
	vimoswcs->xrefpix = vimoswcs->crpix[0];
	vimoswcs->yrefpix = vimoswcs->crpix[1];
	dsspos (vimoswcs->crpix[0], vimoswcs->crpix[1], vimoswcs, &ra0, &dec0);
	vimoswcs->crval[0] = ra0;
	vimoswcs->crval[1] = dec0;
	vimoswcs->xref = vimoswcs->crval[0];
	vimoswcs->yref = vimoswcs->crval[1];

	/* Compute a nominal scale factor */
	dsspos (vimoswcs->crpix[0], vimoswcs->crpix[1]+1.0, vimoswcs, &ra1, &dec1);
	vimoswcs->yinc = dec1 - dec0;
	vimoswcs->xinc = -vimoswcs->yinc;

	/* Compute image rotation angle */
	vimoswcs->vimoswcson = 1;
	vimoswcsrotset (vimoswcs);
	rot = degrad (vimoswcs->rot);

	/* Compute image scale at center */
	dsspos (vimoswcs->crpix[0]+cos(rot),
		vimoswcs->crpix[1]+sin(rot), vimoswcs, &ra1, &dec1);
	vimoswcs->cdelt[0] = -vimoswcsdist (ra0, dec0, ra1, dec1);
	dsspos (vimoswcs->crpix[0]+sin(rot),
		vimoswcs->crpix[1]+cos(rot), vimoswcs, &ra1, &dec1);
	vimoswcs->cdelt[1] = vimoswcsdist (ra0, dec0, ra1, dec1);

	/* Set all other image scale parameters */
	vimoswcsdeltset (vimoswcs, vimoswcs->cdelt[0], vimoswcs->cdelt[1], vimoswcs->rot);
	}

    /* Approximate world coordinate system if plate scale is known */
    else if (ksearch (hstring,"SECPIX") != NULL ||
	     ksearch (hstring,"PIXSCALE") != NULL ||
	     ksearch (hstring,"PIXSCAL1") != NULL ||
	     ksearch (hstring,"SECPIX1") != NULL) {
	secpix = 0.0;
	hgetr8 (hstring,"SECPIX",&secpix);
	if (secpix == 0.0)
	    hgetr8 (hstring,"PIXSCALE",&secpix);
	if (secpix == 0.0) {
	    hgetr8 (hstring,"SECPIX1",&secpix);
	    if (secpix != 0.0) {
		cdelt1 = -secpix / 3600.0;
		hgetr8 (hstring,"SECPIX2",&secpix);
		cdelt2 = secpix / 3600.0;
		}
	    else {
		hgetr8 (hstring,"PIXSCAL1",&secpix);
		cdelt1 = -secpix / 3600.0;
		hgetr8 (hstring,"PIXSCAL2",&secpix);
		cdelt2 = secpix / 3600.0;
		}
	    }
	else {
	    cdelt2 = secpix / 3600.0;
	    cdelt1 = -cdelt2;
	    }

	/* Get rotation angle from the header, if it's there */
	rot = 0.0;
	hgetr8 (hstring,"CROTA1", &rot);
	if (vimoswcs->rot == 0.)
	    hgetr8 (hstring,"CROTA2", &rot);

	/* Set CD and PC matrices */
	vimoswcsdeltset (vimoswcs, cdelt1, cdelt2, rot);

	/* By default, set reference pixel to center of image */
	vimoswcs->crpix[0] = vimoswcs->nxpix * 0.5;
	vimoswcs->crpix[1] = vimoswcs->nypix * 0.5;

	/* Get reference pixel from the header, if it's there */
	if (ksearch (hstring,"CRPIX1") != NULL) {
	    hgetr8 (hstring,"CRPIX1",&vimoswcs->crpix[0]);
	    hgetr8 (hstring,"CRPIX2",&vimoswcs->crpix[1]);
	    }

	/* Use center of detector array as reference pixel
	else if (ksearch (hstring,"DETSIZE") != NULL ||
		 ksearch (hstring,"DETSEC") != NULL) {
	    hgets (hstring, "DETSIZE", 32, temp);
	    ic = strchr (temp, ':');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ',');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ':');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ']');
	    if (ic != NULL)
		*ic = (char) 0;
	    sscanf (temp, "%d %d %d %d", &idx1, &idx2, &idy1, &idy2);
	    dxrefpix = 0.5 * (double) (idx1 + idx2 - 1);
	    dyrefpix = 0.5 * (double) (idy1 + idy2 - 1);
	    hgets (hstring, "DETSEC", 32, temp);
	    ic = strchr (temp, ':');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ',');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ':');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ']');
	    if (ic != NULL)
		*ic = (char) 0;
	    sscanf (temp, "%d %d %d %d", &ix1, &ix2, &iy1, &iy2);
	    vimoswcs->crpix[0] = dxrefpix - (double) (ix1 - 1);
	    vimoswcs->crpix[1] = dyrefpix - (double) (iy1 - 1);
	    } */
	vimoswcs->xrefpix = vimoswcs->crpix[0];
	vimoswcs->yrefpix = vimoswcs->crpix[1];

	vimoswcs->crval[0] = 0.0;
	if (!hgetra (hstring,"RA",&vimoswcs->crval[0])) {
	    setvimoswcserr ("WCSINIT: No RA with SECPIX, no WCS");
	    vimoswcsfree (vimoswcs);
	    return (NULL);
	    }
	vimoswcs->crval[1] = 0.0;
	if (!hgetdec (hstring,"DEC",&vimoswcs->crval[1])) {
	    setvimoswcserr ("WCSINIT No DEC with SECPIX, no WCS");
	    vimoswcsfree (vimoswcs);
	    return (NULL);
	    }
	vimoswcs->xref = vimoswcs->crval[0];
	vimoswcs->yref = vimoswcs->crval[1];
	vimoswcs->coorflip = 0;

	vimoswcs->cel.ref[0] = vimoswcs->crval[0];
	vimoswcs->cel.ref[1] = vimoswcs->crval[1];
	vimoswcs->cel.ref[2] = 999.0;
	hgetr8 (hstring,"LONGPOLE",&vimoswcs->cel.ref[2]);
	vimoswcs->cel.ref[3] = 999.0;
	hgetr8 (hstring,"LATPOLE",&vimoswcs->cel.ref[3]);
	    
	(void) vimoswcstype (vimoswcs, "RA---TAN", "DEC--TAN");
	vimoswcs->coorflip = 0;
	vimoswcs->degout = 0;
	vimoswcs->ndec = 3;

	/* Coordinate reference frame and equinox */
	vimoswcseq (hstring,vimoswcs);

	/* Epoch of image (from observation date, if possible) */
	if (hgetr8 (hstring, "MJD-OBS", &mjd))
	    vimoswcs->epoch = 1900.0 + (mjd - 15019.81352) / 365.242198781;
	else if (!hgetdate (hstring,"DATE-OBS",&vimoswcs->epoch)) {
	    if (!hgetdate (hstring,"DATE",&vimoswcs->epoch)) {
		if (!hgetr8 (hstring,"EPOCH",&vimoswcs->epoch))
		    vimoswcs->epoch = vimoswcs->equinox;
		}
	    }
	vimoswcs->vimoswcson = 1;
	}

    else {
	setvimoswcserr ("WCSINIT: No image scale");
	vimoswcsfree (vimoswcs);
	return (NULL);
	}

    vimoswcs->lin.crpix = vimoswcs->crpix;
    vimoswcs->lin.cdelt = vimoswcs->cdelt;
    vimoswcs->lin.pc = vimoswcs->pc;
    if (strlen (vimoswcs->radecsys) == 0 || vimoswcs->prjcode == VIMOSWCS_LIN)
	strcpy (vimoswcs->radecsys, "LINEAR");
    vimoswcs->sysvimoswcs = vimoswcscsys (vimoswcs->radecsys);

    if (vimoswcs->sysvimoswcs == VIMOSWCS_B1950)
	strcpy (vimoswcs->radecout, "FK4");
    else if (vimoswcs->sysvimoswcs == VIMOSWCS_J2000)
	strcpy (vimoswcs->radecout, "FK5");
    else
	strcpy (vimoswcs->radecout, vimoswcs->radecsys);
    vimoswcs->sysout = vimoswcscsys (vimoswcs->radecout);
    vimoswcs->eqout = vimoswcs->equinox;
    strcpy (vimoswcs->radecin, vimoswcs->radecsys);
    vimoswcs->sysin = vimoswcscsys (vimoswcs->radecin);
    vimoswcs->eqin = vimoswcs->equinox;
    vimoswcs->printsys = 1;
    vimoswcs->tabsys = 0;
    vimoswcs->linmode = 0;

    /* Initialize special WCS commands */
    setvimoswcscom (vimoswcs);

    return (vimoswcs);
}


static void
vimoswcseq (hstring, vimoswcs)

char	*hstring;
struct WorldCoor *vimoswcs;
{
    int ieq = 0;
    int eqhead = 0;
    char systring[32], eqstring[32];

    /* Set equinox from EQUINOX, EPOCH, or RADECSYS; default to 2000 */
    systring[0] = 0;
    eqstring[0] = 0;
    hgets (hstring,"EQUINOX",16,eqstring);
    if (eqstring[0] == 'J') {
	vimoswcs->equinox = atof (eqstring+1);
	ieq = atoi (eqstring+1);
	strcpy (systring, "FK5");
	}
    else if (eqstring[0] == 'B') {
	vimoswcs->equinox = atof (eqstring+1);
	ieq = atoi (eqstring+1);
	strcpy (systring, "FK4");
	}
    else if (hgeti4 (hstring,"EQUINOX",&ieq)) {
	hgetr8 (hstring,"EQUINOX",&vimoswcs->equinox);
	eqhead = 1;
	}

    else if (hgeti4 (hstring,"EPOCH",&ieq)) {
	if (ieq == 0) {
	    ieq = 1950;
	    vimoswcs->equinox = 1950.0;
	    }
	else {
            hgetr8 (hstring,"EPOCH",&vimoswcs->equinox);
	    eqhead = 1;
	    }
	}

    else if (hgets (hstring,"RADECSYS", 16, systring)) {
	if (!strncmp (systring,"FK4",3)) {
	    vimoswcs->equinox = 1950.0;
	    ieq = 1950;
	    }
	else if (!strncmp (systring,"ICRS",4)) {
	    vimoswcs->equinox = 2000.0;
	    ieq = 2000;
	    }
	else if (!strncmp (systring,"FK5",3)) {
	    vimoswcs->equinox = 2000.0;
	    ieq = 2000;
	    }
	else if (!strncmp (systring,"GAL",3)) {
	    vimoswcs->equinox = 2000.0;
	    ieq = 2000;
	    }
	else if (!strncmp (systring,"ECL",3)) {
	    vimoswcs->equinox = 2000.0;
	    ieq = 2000;
	    }
	}

    if (ieq == 0) {
	vimoswcs->equinox = 2000.0;
	ieq = 2000;
	if (vimoswcs->c1type[0] == 'R' || vimoswcs->c1type[0] == 'D')
	    strcpy (systring,"FK5");
	}

    /* Epoch of image (from observation date, if possible) */
    if (!hgetdate (hstring,"DATE-OBS",&vimoswcs->epoch)) {
	if (!hgetdate (hstring,"DATE",&vimoswcs->epoch)) {
	    if (!hgetr8 (hstring,"EPOCH",&vimoswcs->epoch))
		vimoswcs->epoch = vimoswcs->equinox;
	    }
	}
    if (vimoswcs->epoch == 0.0)
	vimoswcs->epoch = vimoswcs->equinox;

    /* Set coordinate system from keyword, if it is present */
    if (systring[0] == (char) 0)
	 hgets (hstring,"RADECSYS", 16, systring);
    if (systring[0] != (char) 0) {
	strcpy (vimoswcs->radecsys,systring);
	if (!eqhead) {
	    if (!strncmp (vimoswcs->radecsys,"FK4",3))
		vimoswcs->equinox = 1950.0;
	    else if (!strncmp (vimoswcs->radecsys,"FK5",3))
		vimoswcs->equinox = 2000.0;
	    else if (!strncmp (vimoswcs->radecsys,"ICRS",4))
		vimoswcs->equinox = 2000.0;
	    else if (!strncmp (vimoswcs->radecsys,"GAL",3) && ieq == 0)
		vimoswcs->equinox = 2000.0;
	    }
	}

    /* Set galactic coordinates if GLON or GLAT are in C1TYPE */
    else if (vimoswcs->c1type[0] == 'G')
	strcpy (vimoswcs->radecsys,"GALACTIC");
    else if (vimoswcs->c1type[0] == 'E')
	strcpy (vimoswcs->radecsys,"ECLIPTIC");
    else if (vimoswcs->c1type[0] == 'S')
	strcpy (vimoswcs->radecsys,"SGALACTC");
    else if (vimoswcs->c1type[0] == 'H')
	strcpy (vimoswcs->radecsys,"HELIOECL");
    else if (vimoswcs->c1type[0] == 'A')
	strcpy (vimoswcs->radecsys,"ALTAZ");
    else if (vimoswcs->c1type[0] == 'L')
	strcpy (vimoswcs->radecsys,"LINEAR");

    /* Otherwise set coordinate system from equinox */
    /* Systemless coordinates cannot be translated using b, j, or g commands */
    else if (vimoswcs->sysvimoswcs != VIMOSWCS_NPOLE) {
	if (ieq > 1980)
	    strcpy (vimoswcs->radecsys,"FK5");
	else
	    strcpy (vimoswcs->radecsys,"FK4");
	}
    vimoswcs->sysvimoswcs = vimoswcscsys (vimoswcs->radecsys);

    return;
}

/* Jun 11 1998	Split off header-dependent WCS initialization from other subs
 * Jun 15 1998	Fix major bug in vimoswcsinit() when synthesizing WCS from header
 * Jun 18 1998	Fix bug in CD initialization; split PC initialization off
 * Jun 18 1998	Split PC initialization off into subroutine vimoswcspcset()
 * Jun 24 1998	Set equinox from RADECSYS only if EQUINOX and EPOCH not present
 * Jul  6 1998  Read third and fourth axis CTYPEs
 * Jul  7 1998  Initialize eqin and eqout to equinox,
 * Jul  9 1998	Initialize rotation matrices correctly
 * Jul 13 1998	Initialize rotation, scale for polynomial and DSS projections
 * Aug  6 1998	Fix CROTA computation for DSS projection
 * Sep  4 1998	Fix CROTA, CDELT computation for DSS and polynomial projections
 * Sep 14 1998	If DATE-OBS not found, check for DATE
 * Sep 14 1998	If B or J present in EQUINOX, use that info to set system
 * Sep 29 1998  Initialize additional WCS commands from the environment
 * Sep 29 1998	Fix bug which read DATE as number rather than formatted date
 * Dec  2 1998	Read projection constants from header (bug fix)
 *
 * Feb  9 1999	Set rotation angle correctly when using DSS projection
 * Feb 19 1999	Fill in CDELTs from scale keyword if absent or zero
 * Feb 19 1999	Add PIXSCALE as possible default arcseconds per pixel
 * Apr  7 1999	Add error checking for NAXIS and NAXIS1 keywords
 * Apr  7 1999	Do not set systring if epoch is 0 and not RA/Dec
 * Jul  8 1999	In RADECSYS, use FK5 and FK4 instead of J2000 and B1950
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 21 1999	Remove unused variables; declare dsspos() after lint
 *
 * Jan 24 2000	Set CD matrix from header even if using polynomial
 * Jan 27 2000	Fix MJD to epoch conversion for when MJD-OBS is the only date
 * Jan 28 2000	Set CD matrix for DSS projection, too
 * Jan 28 2000	Use wcsproj instead of oldwcs
 * Feb 29 2000	Compute inverse CD matrix even if polynomial solution
 * May  9 2000	Add PROJR0 keyword for WCSLIB projections
 * Sep 29 2000	Use CDi_j matrix if any elements are present
 * Oct 11 2000	Fix bug initializing CD matrix
 */
