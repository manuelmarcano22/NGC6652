/*** File libwcs/fitswcs.c
 *** March 8, 2001
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics

 * Module:      fitswcs.c (FITS file WCS reading and deleting)
 * Purpose:     Read and delete FITS image world coordinate system keywords
 * Subroutine:  GetVIMOSWCSFITS (filename)
 *		Open a FITS or IRAF image file and returns its WCS structure
 * Subroutine:  GetFITShead (filename)
 *		Open a FITS or IRAF image file and returns a FITS header
 * Subroutine:  DelVIMOSWCSFITS (filename, verbose)
 *		Delete all standard WCS keywords from a FITS header
 * Subroutine:	PrintVIMOSWCS (header, verbose)
 *		Check the WCS fields and print any that are found if verbose.
 * Subroutine:	SetFITSVIMOSWCS (header, wcs)
 *		Set FITS WCS keywords from WCS data structure
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fitsfile.h"
#include "vimoswcs.h"


struct WorldCoor *
GetVIMOSWCSFITS (filename)

char *filename;	/* FITS or IRAF file filename */

{
    char *header;		/* FITS header */
    struct WorldCoor *vimoswcs;	/* World coordinate system structure */
    char *GetFITShead();
    char *cwcs;			/* Multiple wcs string (name or character) */

    /* Read the FITS or IRAF image file header */
    header = GetFITShead (filename);
    if (header == NULL)
	return (NULL);

    /* Set the world coordinate system from the image header */
    cwcs = strchr (filename, '%');
    if (cwcs != NULL)
	cwcs++;
    vimoswcs = vimoswcsinitn (header, cwcs);
    if (vimoswcs == NULL) {
	setvimoswcsfile (filename);
	vimoswcserr ();
	}
    free (header);

    return (vimoswcs);
}

char *
GetFITShead (filename)

char *filename;	/* FITS or IRAF file filename */

{
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    char *irafheader;		/* IRAF image header */
    int nbiraf, nbfits;

    /* Open IRAF image if .imh extension is present */
    if (isiraf (filename)) {
	if ((irafheader = irafrhead (filename, &nbiraf)) != NULL) {
	    if ((header = iraf2fits (filename, irafheader, nbiraf, &lhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s\n",filename);
		free (irafheader);
		irafheader = NULL;
		return (NULL);
		}
	    free (irafheader);
	    irafheader = NULL;
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    return (NULL);
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	if ((header = fitsrhead (filename, &lhead, &nbfits)) == NULL) {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return (NULL);
	    }
	}

    return (header);
}


/* delete all the C* fields.
 * return 0 if at least one such field is found, else -1.  */

int
DelVIMOSWCSFITS (header, verbose)

char *header;
int verbose;

{
    static char flds[15][8];
    int i;
    int n, nfields;
    double eq;
    char rastr[16],decstr[16];

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
    nfields = 15;

    for (i = 0; i < nfields; i++) {
	if (hdel (header, flds[i])) {
	    n++;
	    if (verbose)
		fprintf (stderr,"%s: deleted\n", flds[i]);
	    }
	}
    if (verbose && n == 0)
	fprintf (stderr,"DelVIMOSWCSFITS: No WCS in header\n");

    /* Delete RA DEC EPOCH, replacing with saved values, if present */
    if (ksearch (header,"WRA")) {
	hdel (header, "RA");
	n++;
	hchange (header, "WRA","RA");
	if (ksearch (header,"WDEC")) {
	    hdel (header, "DEC");
	    n++;
	    hchange (header, "WDEC", "DEC");
	    }
	if (ksearch (header,"WEPOCH")) {
	    hdel (header, "EPOCH");
	    n++;
	    hchange (header, "WEPOCH", "EPOCH");
	    }
	if (ksearch (header,"WEQUINOX")) {
	    hdel (header, "EQUINOX");
	    n++;
	    hchange (header, "WEQUINOX", "EQUINOX");
	    }
	if (ksearch (header, "EPOCH")) {
	    hdel (header, "EQUINOX");
	    n++;
	    if (verbose)
		fprintf (stderr,"EQUINOX deleted\n");
	    }
	hdel (header, "RADECSYS");
	n++;
	if (verbose)
	    fprintf (stderr,"RADECSYS deleted\n");
	hdel (header, "SECPIX1");
	n++;
	if (verbose)
	    fprintf (stderr,"SECPIX1 deleted\n");
	hdel (header, "SECPIX2");
	n++;
	if (verbose)
	    fprintf (stderr,"SECPIX2 deleted\n");
	if (verbose) {
	    hgets (header,"RA", 16, rastr);
	    hgets (header,"DEC", 16, decstr);
	    eq = 0.0;
	    hgetr8 (header,"EPOCH",&eq);
	    if (eq == 0.0)
		hgetr8 (header,"EQUINOX",&eq);
	    fprintf (stderr,"DelVIMOSWCS: Center reset to %s %s %.1f\n", rastr,decstr, eq);
	    }
	}
    else if (ksearch (header, "EPOCH") && !ksearch (header, "PLTRAH")) {
	if (hdel (header,"EQUINOX")) {
	    if (verbose)
		fprintf (stderr,"EQUINOX: deleted\n");
	    n++;
	    }
	else if (verbose)
	    fprintf (stderr,"DelVIMOSWCS: EPOCH, but not EQUINOX found\n");
	}

    /* Delete SAO polynomial, if present */
    if (ksearch (header, "CO1_1")) {
	int i;
	char keyword[16];

	for (i = 1; i < 13; i++) {
	    sprintf (keyword,"CO1_%d", i);
	    hdel (header, keyword);
	    if (verbose)
		fprintf (stderr,"%s deleted\n", keyword);
	    n++;
	    }
	for (i = 1; i < 13; i++) {
	    sprintf (keyword,"CO2_%d", i);
	    hdel (header, keyword);
	    if (verbose)
		fprintf (stderr,"%s deleted\n", keyword);
	    n++;
	    }
	}

    /* Delete rotation matrix, if present */
    if (ksearch (header, "CO1_1")) {
	int i, j;
	char keyword[16];
	for (i = 1; i < 6; i++) {
	    for (j = 1; i < 6; i++) {
		sprintf (keyword,"PC%03d%03d", i, j);
		hdel (header, keyword);
		if (verbose)
		    fprintf (stderr,"%s deleted\n", keyword);
		n++;
		}
	    }
	}

    if (n > 0 && verbose)
	fprintf (stderr,"%d keywords deleted\n", n);

    return (n);
}


/* check the WCS fields and print any that are found if verbose.
 * return 0 if all are found, else -1.
 */
int
PrintVIMOSWCS (header, verbose)
char	*header;	/* FITS header */
int	verbose;	/* 1 to print WCS header keyword values */

{
    char str[80];
    double v;
    int n, i;
    char keyword[16];

    n = 0;

    if (hgets (header,"IMWCS",80,str)) {
	if (verbose) fprintf (stderr,"IMWCS = %s\n", str);
	n++;
	}
    if (hgets (header,"CTYPE1",16,str)) {
	if (verbose) fprintf (stderr,"CTYPE1 = %s\n", str);
	n++;
	}
    if (hgetr8 (header, "CRVAL1", &v)) {
	if (verbose) fprintf (stderr,"CRVAL1 = %.8f\n", v);
	n++;
	}
    if (hgetr8 (header, "CRPIX1", &v)) {
	if (verbose) fprintf (stderr,"CRPIX1 = %.8f\n", v);
	n++;
	}

    if (hgets (header,"CTYPE2",16,str)) {
	if (verbose) fprintf (stderr,"CTYPE2 = %s\n", str);
	n++;
	}
    if (hgetr8 (header, "CRVAL2", &v)) {
	if (verbose) fprintf (stderr,"CRVAL2 = %.8f\n", v);
	n++;
	}
    if (hgetr8 (header, "CRPIX2", &v)) {
	if (verbose) fprintf (stderr,"CRPIX2 = %.8f\n", v);
	n++;
	}

    /* Polynomial plate fit */
    if (hgetr8 (header, "CO1_1", &v)) {
	if (verbose) fprintf (stderr,"CO1_1 = %.8g\n", v);
	for (i = 1; i < 20; i++) {
	    sprintf (keyword,"CO1_%d",i+1);
	    if (hgetr8 (header, keyword, &v)) {
		if (verbose) fprintf (stderr,"%s = %.8g\n", keyword, v);
		n++;
		}
	    }
	}
    if (hgetr8 (header, "CO2_1", &v)) {
	if (verbose) fprintf (stderr,"CO2_1 = %.8g\n", v);
	for (i = 1; i < 20; i++) {
	    sprintf (keyword,"CO2_%d",i+1);
	    if (hgetr8 (header, keyword, &v)) {
		if (verbose) fprintf (stderr,"%s = %.8g\n", keyword, v);
		n++;
		}
	    }
	}

    /* Plate scale and rotation from CD matrix */
    if (hgetr8 (header, "CD1_1", &v)) {
	if (verbose) fprintf (stderr,"CD1_1 = %.8g\n", v);
	n++;
	if (hgetr8 (header, "CD1_2", &v)) {
	    if (verbose) fprintf (stderr,"CD1_2 = %.gf\n", v);
	    n++;
	    }
	if (hgetr8 (header, "CD2_1", &v)) {
	    if (verbose) fprintf (stderr,"CD2_1 = %.gf\n", v);
	    n++;
	    }
	if (hgetr8 (header, "CD2_2", &v)) {
	    if (verbose) fprintf (stderr,"CD2_2 = %.gf\n", v);
	    n++;
	    }
	}

    /* Plate scale and rotation from CDELTn and CROTAn */
    else {
	if (hgetr8 (header, "CDELT1", &v)) {
	    if (verbose) fprintf (stderr,"CDELT1 = %.8f\n", v);
	    n++;
	    }
	if (hgetr8 (header, "CROTA1", &v)) {
	    if (verbose) fprintf (stderr,"CROTA1 = %.3f\n", v);
	    n++;
	    }
	if (hgetr8 (header, "CDELT2", &v)) {
	    if (verbose) fprintf (stderr,"CDELT2 = %.8f\n", v);
	    n++;
	    }
	if (hgetr8 (header, "CROTA2", &v)) {
	    if (verbose) fprintf (stderr,"CROTA2 = %.3f\n", v);
	    n++;
	    }
	}

    return (n > 8 ? 0 : -1);
}

static char vimoswcsproj[8]="TAN";		/* WCS projection name */
void
setvimoswcsproj (type)
char *type;
{ strcpy (vimoswcsproj, type); return; }


/* Set FITS C* fields, assuming RA/DEC refers to the reference pixel, CRPIX1/CRPIX2 */

void
SetFITSVIMOSWCS (header, vimoswcs)

char	*header;	/* Image FITS header */
struct WorldCoor *vimoswcs;	/* WCS structure */

{
    double ep;
    char wcstemp[16];

    /* Rename old center coordinates */
    if (!ksearch (header,"WRA") && ksearch (header,"RA"))
	hchange (header,"RA","WRA");
    if (!ksearch (header,"WDEC") && ksearch (header,"DEC"))
	hchange (header,"DEC","WDEC");

    if (!ksearch (header,"WEQUINOX") && ksearch (header,"EQUINOX"))
	hchange (header, "EQUINOX", "WEQUINOX");

    /* Only change EPOCH if it is used instead of EQUINOX */
    else if (!ksearch (header,"WEPOCH") && ksearch (header,"EPOCH"))
	hchange (header, "EPOCH", "WEPOCH");


    /* Set new center coordinates */
    hputra (header,"RA",vimoswcs->xref);
    hputdec (header,"DEC",vimoswcs->yref);
    hputr8 (header, "EQUINOX", vimoswcs->equinox);
    if (hgetr8 (header, "WEPOCH", &ep))
	hputr8 (header, "EPOCH", vimoswcs->equinox);
    else if (!hgetr8 (header, "EPOCH", &ep))
	hputr8 (header, "EPOCH", vimoswcs->equinox);

    if (vimoswcs->radecsys[0] == 'B' || vimoswcs->radecsys[0] == 'b')
	hputs (header, "RADECSYS", "FK4");
    else if (vimoswcs->radecsys[0] == 'J' || vimoswcs->radecsys[0] == 'j')
	hputs (header, "RADECSYS", "FK5");
    else
	hputs (header, "RADECSYS", vimoswcs->radecsys);

    /* Set standard FITS WCS keywords */
    strcpy (wcstemp, "RA---");
    strcat (wcstemp, vimoswcsproj);
    hputs  (header, "CTYPE1", wcstemp);
    strcpy (wcstemp, "DEC--");
    strcat (wcstemp, vimoswcsproj);
    hputs  (header, "CTYPE2", wcstemp);

    /* Reference pixel in WCS and image coordinates */
    hputnr8 (header, "CRVAL1", 9, vimoswcs->xref);
    hputnr8 (header, "CRVAL2", 9, vimoswcs->yref);
    hputnr8 (header, "CRPIX1", 4, vimoswcs->xrefpix);
    hputnr8 (header, "CRPIX2", 4, vimoswcs->yrefpix);

    /* CD matrix (proposed FITS standard) */
    if (vimoswcs->rotmat) {
	hputnr8 (header, "CD1_1", 9, vimoswcs->cd[0]);
	hputnr8 (header, "CD1_2", 9, vimoswcs->cd[1]);
	hputnr8 (header, "CD2_1", 9, vimoswcs->cd[2]);
	hputnr8 (header, "CD2_2", 9, vimoswcs->cd[3]);
	hdel (header, "CDELT1");
	hdel (header, "CDELT2");
	hdel (header, "CROTA1");
	hdel (header, "CROTA2");
	}

    /* Scale and rotation (old FITS standard) */
    else {
	hputnr8 (header, "CDELT1", 9, vimoswcs->xinc);
	hputnr8 (header, "CDELT2", 9, vimoswcs->yinc);
	hputnr8 (header, "CROTA1", 3, vimoswcs->rot);
	hputnr8 (header, "CROTA2", 3, vimoswcs->rot);
	hdel (header, "CD1_1");
	hdel (header, "CD1_2");
	hdel (header, "CD2_1");
	hdel (header, "CD2_2");
	}

    /* Plate scale at reference pixel */
    if (-vimoswcs->xinc != vimoswcs->yinc) {
	if (ksearch (header,"SECPIX"))
	    hdel (header,"SECPIX");
	hputnr8 (header, "SECPIX1", 4, -vimoswcs->xinc*3600.0);
	hputnr8 (header, "SECPIX2", 4, vimoswcs->yinc*3600.0);
	}
    else {
	if (ksearch (header,"SECPIX1"))
	    hdel (header,"SECPIX1");
	if (ksearch (header,"SECPIX2"))
	    hdel (header,"SECPIX2");
	hputnr8 (header, "SECPIX", 4, vimoswcs->yinc*3600.0);
	}

    /* Plate fit coefficients, if present */
    if (vimoswcs->ncoeff1 > 0) {
	char keyword[16];
	int i;
	for (i = 0; i < vimoswcs->ncoeff1; i++) {
	    sprintf (keyword, "CO1_%d",i+1);
	    hputr8 (header, keyword, vimoswcs->x_coeff[i]);
	    }
	}
    if (vimoswcs->ncoeff2 > 0) {
	char keyword[16];
	int i;
	for (i = 0; i < vimoswcs->ncoeff2; i++) {
	    sprintf (keyword, "CO2_%d",i+1);
	    hputr8 (header, keyword, vimoswcs->y_coeff[i]);
	    }
	}

    return;
}
/* May 29 1996	Change name from delVIMOSWCSFITS to DelVIMOSWCSFITS
 * May 31 1996	Print single message if no WCS is found in header
 * May 31 1996	Use stream I/O instead of standard I/O
 * Jun 10 1996	Combine imgetvimoswcs.c and imdelvimoswcs.c into fitsvimoswcs.c
 * Jun 17 1996	Delete IMWCS record, too
 * Jul 16 1996	Update arguments for header-reading subroutines
 * Aug  6 1996  Fixed small defects after lint
 * Aug  8 1996  Restore old image center after deleting WCS
 * Aug 26 1996	Fix subroutine arguments after lint
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Feb 21 1997  Add GetFITShead() subroutine and use it
 * Mar 20 1997	Remove unused variables
 * Nov  6 1997	Add PrintVIMOSWCS() from IMWCS
 *
 * Jan  7 1998	Return NULL WCS structure if no FITS header can be read
 * Feb 18 1998	Move SetFITSVIMOSWCS() here from imsetvimoswcs.c
 * Feb 24 1998	Delete CD matrix in DelVIMOSWCS()
 * Mar 20 1998	Write CD matrix in SetFITSVIMOSWCS()
 * Mar 27 1998	Add plate constants in SetFITSVIMOSWCS()
 * Mar 27 1998	Delete plate constants in DelFITSVIMOSWCS()
 * Apr  6 1998	Change coefficient keywords from PLTij to COi_j
 * Apr  7 1998	Change amd_i_coeff to i_coeff
 * Apr 10 1998	Write out polynomial coefficients correctly
 * Apr 13 1998	Print polynomial coefficients, if in header
 * Apr 16 1998	Drop NCOEFF header parameter
 * Apr 17 1998	Do not write W* keywords if they are already there
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun  1 1998	Print error message if WCS cannot be initialized
 * Jun 11 1998	Change WCSTYPE to WCSPROJ to avoid conflict
 * Jul 23 1998	In DelVIMOSWCS, delete specific number of fields
 * Jul 24 1998	Make irafheader char instead of int
 * Jul 27 1998	Set irafheader pointer to NULL after use
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Oct  5 1998	Use isiraf() to determine file type
 * Oct 28 1998	Delete EQUINOX, RADECSYS, SECPIX1, SECPIX2 from imvimoswcs
 *
 * Apr  7 1999	Add file name to error message if WCS error
 * Jul  8 1999	Write RADECSYS as FK5 or FK4 instead of J2000 or B1950
 * Jul 21 1999	Add SECPIX plate scale output to SetFITSVIMOSWCS()
 * Oct 21 1999	Fix declarations after lint
 *
 * Mar 23 2000	Fix bug in IRAF header error message
 *
 * Jan 11 2001	Print all messages to stderr
 * Jan 31 2001	Add code to extract WCS name or character from filename
 * Mar  8 2001	Change WCS character separator from : to % in FITS filenames
 */
