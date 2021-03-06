/*** File libwcs/imrotate.c
 *** January 18, 2001
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitshead.h"
#include "imio.h"

static void RotVIMOSWCSFITS();

/* Rotate an image by 90, 180, or 270 degrees, with an optional
 * reflection across the vertical axis.
 * verbose generates extra info on stdout.
 * return NULL if successful or rotated image.
 */

char *
RotFITS (pathname, header, image, rotate, mirror, bitpix2, verbose)

char	*pathname;	/* Name of file which is being changed */
char	*header;	/* FITS header */
char	*image;		/* Image pixels */
int	rotate;		/* Angle to by which to rotate image (90, 180, 270) */
int	mirror;		/* 1 to reflect image around vertical axis */
int	bitpix2;	/* Number of bits per pixel in output image */
int	verbose;

{
    int bitpix1, ny, nx, nax;
    int x1, y1, x2, y2, nbytes;
    char *rotimage;
    char history[128];
    char *filename;

    if (rotate == 1)
	rotate = 90;
    else if (rotate == 2)
	rotate = 180;
    else if (rotate == 3)
	rotate = 270;
    else if (rotate < 0)
	rotate = rotate + 360;

    filename = strrchr (pathname,'/');
    if (filename)
	filename = filename + 1;
    else
	filename = pathname;

    /* Get image size */
    nax = 0;
    if (hgeti4 (header,"NAXIS",&nax) < 1)
	return (NULL);
    else {
	if (hgeti4 (header,"NAXIS1",&nx) < 1)
	    return (NULL);
	else {
	    if (hgeti4 (header,"NAXIS2",&ny) < 1)
		return (NULL);
	    }
	}
    bitpix1 = 16;
    hgeti4 (header,"BITPIX", &bitpix1);
    if (bitpix2 == 0)
	bitpix2 = bitpix1;

    /* Delete WCS fields in header */
    if (rotate != 0 || mirror)
	RotVIMOSWCSFITS (header, rotate, mirror, verbose);

    /* Compute size of image in bytes */
    switch (bitpix2) {
	case 16:
	    nbytes = nx * ny * 2;
	    break;
	case 32:
	    nbytes = nx * ny * 4;
	    break;
	case -16:
	    nbytes = nx * ny * 2;
	    break;
	case -32:
	    nbytes = nx * ny * 4;
	    break;
	case -64:
	    nbytes = nx * ny * 8;
	    break;
	default:
	    return (NULL);
	}

    /* Allocate buffer for rotated image */
    rotimage = (char *) malloc (nbytes);
    if (rotimage == NULL)
	return (NULL);

    if (bitpix1 != bitpix2) {
	sprintf (history,"Copy of image %s bits per pixel %d -> %d",
		filename, bitpix1, bitpix2);
	hputc (header,"HISTORY",history);
	if (verbose)
	    fprintf (stderr,"%s\n",history);
	}

    /* Mirror image without rotation */
    if (rotate < 45.0 && rotate > -45.0) {
	if (mirror) {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = nx - x1 - 1;
		    y2 = y1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,nx,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s reflected",filename);
	    hputc (header,"HISTORY",history);
	    }
	else {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = x1;
		    y2 = y1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,nx,x2,y2);
		    }
		}
	    }
	}

    /* Rotate by 90 degrees */
    else if (rotate >= 45 && rotate < 135) {
	if (mirror) {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = ny - y1 - 1;
		    y2 = nx - x1 - 1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s reflected, rotated 90 degrees",
		     filename);
            hputc (header,"HISTORY",history);
	    }
	else {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = ny - y1 - 1;
		    y2 = x1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s rotated 90 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	hputi4 (header,"NAXIS1",ny);
	hputi4 (header,"NAXIS2",nx);
	}

    /* Rotate by 180 degrees */
    else if (rotate >= 135 && rotate < 225) {
	if (mirror) {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = x1;
		    y2 = ny - y1 - 1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,nx,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s reflected, rotated 180 degrees",
		     filename);
            hputc (header,"HISTORY",history);
	    }
	else {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = nx - x1 - 1;
		    y2 = ny - y1 - 1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,nx,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s rotated 180 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	}

    /* Rotate by 270 degrees */
    else if (rotate >= 225 && rotate < 315) {
	if (mirror) {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = y1;
		    y2 = x1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s reflected, rotated 270 degrees",
		     filename);
            hputc (header,"HISTORY",history);
	    }
	else {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = y1;
		    y2 = nx - x1 - 1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s rotated 270 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	hputi4 (header,"NAXIS1",ny);
	hputi4 (header,"NAXIS2",nx);
	}

    /* If rotating by more than 315 degrees, assume top-bottom reflection */
    else if (rotate >= 315 && mirror) {
	for (y1 = 0; y1 < ny; y1++) {
	    for (x1 = 0; x1 < nx; x1++) {
		x2 = y1;
		y2 = x1;
		movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		}
	    }
	sprintf (history,"Copy of image %s reflected top to bottom",filename);
        hputc (header,"HISTORY",history);
	}
    
    if (verbose)
	fprintf (stderr,"%s\n",history);

    return (rotimage);
}


/* rotate all the C* fields.
 * return 0 if at least one such field is found, else -1.  */

static void
RotVIMOSWCSFITS (header, angle, mirror, verbose)

char	*header;	/* FITS header */
int	angle;		/* Angle to be rotated (0, 90, 180, 270) */
int	mirror;		/* 1 if mirrored left to right, else 0 */
int	verbose;	/* Print progress if 1 */

{
    static char flds[15][8];
    char ctype1[16], ctype2[16];
    double ctemp1, ctemp2, ctemp3, ctemp4, naxis1, naxis2;
    int i, n, ndec1, ndec2, ndec3, ndec4;

    strcpy (flds[0], "CTYPE1");
    strcpy (flds[1], "CTYPE2");
    strcpy (flds[2], "CRVAL1");
    strcpy (flds[3], "CRVAL2");
    strcpy (flds[4], "CDELT1");
    strcpy (flds[5], "CDELT2");
    strcpy (flds[6], "CRPIX1");
    strcpy (flds[7], "CRPIX2");
    strcpy (flds[8], "CROTA1");
    strcpy (flds[9], "CROTA2");
    strcpy (flds[10], "IMWCS");
    strcpy (flds[11], "CD1_1");
    strcpy (flds[12], "CD1_2");
    strcpy (flds[13], "CD2_1");
    strcpy (flds[14], "CD2_2");

    n = 0;
    hgetr8 (header, "NAXIS1", &naxis1);
    hgetr8 (header, "NAXIS2", &naxis2);

    /* Find out if there any WCS keywords in this header */
    for (i = 0; i < sizeof(flds)/sizeof(flds[0]); i++) {
	if (ksearch (header, flds[i]) != NULL) {
	    n++;
	    if (verbose)
		fprintf (stderr,"%s: found\n", flds[i]);
	    }
	}

    /* Return if no WCS keywords to change */
    if (n == 0) {
	if (verbose)
	    fprintf (stderr,"RotVIMOSWCSFITS: No WCS in header\n");
	return;
	}

    /* Reset CTYPEn and CRVALn if axes have been exchanged */
    if (angle == 90 || angle == 270) {
	if (hgets (header, "CTYPE1", 16, ctype1) &&
	    hgets (header, "CTYPE2", 16, ctype2)) {
	    hputs (header, "CTYPE1", ctype2);
	    hputs (header, "CTYPE2", ctype1);
	    }
	if (hgetr8 (header, "CRVAL1", &ctemp1) &&
	    hgetr8 (header, "CRVAL2", &ctemp2)) { 
	    hgetndec (header, "CRVAL1", &ndec1);
	    hgetndec (header, "CRVAL2", &ndec2);
	    hputnr8 (header, "CRVAL1", ndec2, ctemp2);
	    hputnr8 (header, "CRVAL2", ndec1, ctemp1);
	    }
	if (hgets (header, "CUNIT1", 16, ctype1) &&
	    hgets (header, "CUNIT2", 16, ctype2)) {
	    hputs (header, "CUNIT1", ctype2);
	    hputs (header, "CUNIT2", ctype1);
	    }
	}

    /* Negate rotation angle if mirrored */
    if (mirror) {
	if (hgetr8 (header, "CROTA1", &ctemp1)) {
	    hgetndec (header, "CROTA1", &ndec1);
	    hputnr8 (header, "CROTA1", ndec1, -ctemp1);
	    }
	if (hgetr8 (header, "CROTA2", &ctemp2)) {
	    hgetndec (header, "CROTA2", &ndec2);
	    hputnr8 (header, "CROTA2", ndec2, -ctemp2);
	    }
	if (hgetr8 (header, "LTM1_1", &ctemp1)) {
	    hgetndec (header, "LTM1_1", &ndec1);
	    hputnr8 (header, "LTM1_1", ndec1, -ctemp1);
	    }
	if (hgetr8 (header, "CD1_1", &ctemp1))
	    hputr8 (header, "CD1_1", -ctemp1);
	if (hgetr8 (header, "CD1_2", &ctemp1))
	    hputr8 (header, "CD1_2", -ctemp1);
	if (hgetr8 (header, "CD2_1", &ctemp1))
	    hputr8 (header, "CD2_1", -ctemp1);
	}

    /* Unbin CRPIX and CD matrix */
    if (hgetr8 (header, "LTM1_1", &ctemp1)) {
	if (ctemp1 != 1.0) {
	    if (hgetr8 (header, "LTM2_2", &ctemp2)) {
		if (ctemp1 == ctemp2) {
		    double ltv1 = 0.0;
		    double ltv2 = 0.0;
		    if (hgetr8 (header, "LTV1", &ltv1))
			hdel (header, "LTV1");
		    if (hgetr8 (header, "LTV2", &ltv1))
			hdel (header, "LTV2");
		    if (hgetr8 (header, "CRPIX1", &ctemp3))
			hputr8 (header, "CRPIX1", (ctemp3-ltv1)/ctemp1);
		    if (hgetr8 (header, "CRPIX2", &ctemp3))
			hputr8 (header, "CRPIX2", (ctemp3-ltv2)/ctemp1);
		    if (hgetr8 (header, "CD1_1", &ctemp3))
			hputr8 (header, "CD1_1", ctemp3/ctemp1);
		    if (hgetr8 (header, "CD1_2", &ctemp3))
			hputr8 (header, "CD1_2", ctemp3/ctemp1);
		    if (hgetr8 (header, "CD2_1", &ctemp3))
			hputr8 (header, "CD2_1", ctemp3/ctemp1);
		    if (hgetr8 (header, "CD2_2", &ctemp3))
			hputr8 (header, "CD2_2", ctemp3/ctemp1);
		    hdel (header, "LTM1_1");
		    hdel (header, "LTM2_2");
		    }
		}
	    }
	}

    /* Reset CRPIXn */
    if (hgetr8 (header, "CRPIX1", &ctemp1) &&
	hgetr8 (header, "CRPIX2", &ctemp2)) { 
	hgetndec (header, "CRPIX1", &ndec1);
	hgetndec (header, "CRPIX2", &ndec2);
	if (mirror) {
	    if (angle == 0)
		hputnr8 (header, "CRPIX1", ndec1, naxis1-ctemp1);
	    else if (angle == 90) {
		hputnr8 (header, "CRPIX1", ndec2, naxis2-ctemp2);
		hputnr8 (header, "CRPIX2", ndec1, naxis1-ctemp1);
		}
	    else if (angle == 180) {
		hputnr8 (header, "CRPIX1", ndec1, ctemp1);
		hputnr8 (header, "CRPIX2", ndec2, naxis2-ctemp2);
		}
	    else if (angle == 270) {
		hputnr8 (header, "CRPIX1", ndec2, ctemp2);
		hputnr8 (header, "CRPIX2", ndec1, ctemp1);
		}
	    }
	else {
	    if (angle == 90) {
		hputnr8 (header, "CRPIX1", ndec2, naxis2-ctemp2);
		hputnr8 (header, "CRPIX2", ndec1, ctemp1);
		}
	    else if (angle == 180) {
		hputnr8 (header, "CRPIX1", ndec1, naxis1-ctemp1);
		hputnr8 (header, "CRPIX2", ndec2, naxis2-ctemp2);
		}
	    else if (angle == 270) {
		hputnr8 (header, "CRPIX1", ndec2, ctemp2);
		hputnr8 (header, "CRPIX2", ndec1, naxis1-ctemp1);
		}
	    }
	}

    /* Reset CDELTn (degrees per pixel) */
    if (hgetr8 (header, "CDELT1", &ctemp1) &&
	hgetr8 (header, "CDELT2", &ctemp2)) { 
	hgetndec (header, "CDELT1", &ndec1);
	hgetndec (header, "CDELT2", &ndec2);
	if (mirror) {
	    if (angle == 0)
		hputnr8 (header, "CDELT1", ndec1, -ctemp1);
	    else if (angle == 90) {
		hputnr8 (header, "CDELT1", ndec2, -ctemp2);
		hputnr8 (header, "CDELT2", ndec1, -ctemp1);
		}
	    else if (angle == 180) {
		hputnr8 (header, "CDELT1", ndec1, ctemp1);
		hputnr8 (header, "CDELT2", ndec2, -ctemp2);
		}
	    else if (angle == 270) {
		hputnr8 (header, "CDELT1", ndec2, ctemp2);
		hputnr8 (header, "CDELT2", ndec1, ctemp1);
		}
	    }
	else {
	    if (angle == 90) {
		hputnr8 (header, "CDELT1", ndec2, -ctemp2);
		hputnr8 (header, "CDELT2", ndec1, ctemp1);
		}
	    else if (angle == 180) {
		hputnr8 (header, "CDELT1", ndec1, -ctemp1);
		hputnr8 (header, "CDELT2", ndec2, -ctemp2);
		}
	    else if (angle == 270) {
		hputnr8 (header, "CDELT1", ndec2, ctemp2);
		hputnr8 (header, "CDELT2", ndec1, -ctemp1);
		}
	    }
	}

    /* Reset CD matrix, if present */
    ctemp1 = 0.0;
    ctemp2 = 0.0;
    ctemp3 = 0.0;
    ctemp4 = 0.0;
    if (hgetr8 (header, "CD1_1", &ctemp1)) {
	hgetr8 (header, "CD1_2", &ctemp2);
	hgetr8 (header, "CD2_1", &ctemp3);
	hgetr8 (header, "CD2_2", &ctemp4);
	hgetndec (header, "CD1_1", &ndec1);
	hgetndec (header, "CD1_2", &ndec2);
	hgetndec (header, "CD2_1", &ndec3);
	hgetndec (header, "CD2_2", &ndec4);
	if (mirror) {
	    if (angle == 0) {
		hputnr8 (header, "CD1_2", ndec2, -ctemp2);
		hputnr8 (header, "CD2_1", ndec3, -ctemp3);
		}
	    else if (angle == 90) {
		hputnr8 (header, "CD1_1", ndec4, -ctemp4);
		hputnr8 (header, "CD1_2", ndec3, -ctemp3);
		hputnr8 (header, "CD2_1", ndec2, -ctemp2);
		hputnr8 (header, "CD2_2", ndec1, -ctemp1);
		}
	    else if (angle == 180) {
		hputnr8 (header, "CD1_1", ndec1, ctemp1);
		hputnr8 (header, "CD1_2", ndec2, ctemp2);
		hputnr8 (header, "CD2_1", ndec3, -ctemp3);
		hputnr8 (header, "CD2_2", ndec4, -ctemp4);
		}
	    else if (angle == 270) {
		hputnr8 (header, "CD1_1", ndec4, ctemp4);
		hputnr8 (header, "CD1_2", ndec3, ctemp3);
		hputnr8 (header, "CD2_1", ndec2, ctemp2);
		hputnr8 (header, "CD2_2", ndec1, ctemp1);
		}
	    }
	else {
	    if (angle == 90) {
		hputnr8 (header, "CD1_1", ndec4, -ctemp4);
		hputnr8 (header, "CD1_2", ndec3, -ctemp3);
		hputnr8 (header, "CD2_1", ndec2, ctemp2);
		hputnr8 (header, "CD2_2", ndec1, ctemp1);
		}
	    else if (angle == 180) {
		hputnr8 (header, "CD1_1", ndec1, -ctemp1);
		hputnr8 (header, "CD1_2", ndec2, -ctemp2);
		hputnr8 (header, "CD2_1", ndec3, -ctemp3);
		hputnr8 (header, "CD2_2", ndec4, -ctemp4);
		}
	    else if (angle == 270) {
		hputnr8 (header, "CD1_1", ndec4, ctemp4);
		hputnr8 (header, "CD1_2", ndec3, ctemp3);
		hputnr8 (header, "CD2_1", ndec2, -ctemp2);
		hputnr8 (header, "CD2_2", ndec1, -ctemp1);
		}
	    }
	}

    /* Delete any polynomial solution */
    /* (These could maybe be switched, but I don't want to work them out yet */
    if (ksearch (header, "CO1_1")) {
	int i;
	char keyword[16];

	for (i = 1; i < 13; i++) {
	    sprintf (keyword,"CO1_%d", i);
	    hdel (header, keyword);
	    }
	for (i = 1; i < 13; i++) {
	    sprintf (keyword,"CO2_%d", i);
	    hdel (header, keyword);
	    }
	}

    return;
}

/* May 29 1996	Change name from rotFITS to RotFITS
 * Jun  4 1996	Fix bug when handling assymetrical images
 * Jun  5 1996	Print filename, not pathname, in history
 * Jun 10 1996	Remove unused variables after running lint
 * Jun 13 1996	Replace image with rotated image
 * Jun 18 1996	Fix formatting bug in history
 *
 * Jul 11 1997	If rotation is 360, flip top bottom if mirror flat is set
 *
 * Feb 23 1998	Do not delete WCS if image not rotated or mirrored
 * May 26 1998	Rotate WCS instead of deleting it
 * May 27 1998	Include imio.h

 * Jun  8 1999	Return new image pointer instead of flag; do not free old image
 * Jun  9 1999	Make history buffer 128 instead of 72 to avoid overflows
 * Jun 10 1999	Drop image0; use image
 * Oct 21 1999	Fix hputnr8() calls after lint
 *
 * Jan 11 2001	Print all messages to stderr
 * Jan 17 2001	Reset coordinate direction if image is mirrored
 * Jan 18 2001	Reset WCS scale if image is binned
 */
