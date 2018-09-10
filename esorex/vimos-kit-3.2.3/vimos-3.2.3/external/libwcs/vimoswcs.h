/* libwcs/wcs.h
 * March 20, 2001
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics */

#ifndef _vimoswcs_h_
#define _vimoswcs_h_

#include "vimoswcslib.h"
#include "fitshead.h"

struct WorldCoor {
  double	xref;		/* X reference coordinate value (deg) */
  double	yref;		/* Y reference coordinate value (deg) */
  double	xrefpix;	/* X reference pixel */
  double	yrefpix;	/* Y reference pixel */
  double	xinc;		/* X coordinate increment (deg) */
  double	yinc;		/* Y coordinate increment (deg) */
  double	rot;		/* rotation around axis (deg) (N through E) */
  double	cd[4];		/* rotation matrix */
  double	dc[4];		/* inverse rotation matrix */
  double	equinox;	/* Equinox of coordinates default to 1950.0 */
  double	epoch;		/* Epoch of coordinates default to equinox */
  double	nxpix;		/* Number of pixels in X-dimension of image */
  double	nypix;		/* Number of pixels in Y-dimension of image */
  double	plate_ra;	/* Right ascension of plate center */
  double	plate_dec;	/* Declination of plate center */
  double	plate_scale;	/* Plate scale in arcsec/mm */
  double	x_pixel_offset;	/* X pixel offset of image lower right */
  double	y_pixel_offset;	/* Y pixel offset of image lower right */
  double	x_pixel_size;	/* X pixel_size */
  double	y_pixel_size;	/* Y pixel_size */
  double	ppo_coeff[6];	/* pixel to plate coefficients for DSS */
  double	x_coeff[20];	/* X coefficients for plate model */
  double	y_coeff[20];	/* Y coefficients for plate model */
  double	xpix;		/* X (RA) coordinate (pixels) */
  double	ypix;		/* Y (dec) coordinate (pixels) */
  double	zpix;		/* Z (face) coordinate (pixels) */
  double	xpos;		/* X (RA) coordinate (deg) */
  double	ypos;		/* Y (dec) coordinate (deg) */
  double	crpix[4];	/* Values of CRPIXn keywords */
  double	crval[4];	/* Values of CRVALn keywords */
  double	cdelt[4];	/* Values of CDELTn keywords */
  double	pc[16];		/* Values of PCiiijjj keywords */
  double	projp[10];	/* Constants for various projections */
  double	longpole;	/* Longitude of North Pole in degrees */
  double	latpole;	/* Latitude of North Pole in degrees */
  double	rodeg;		/* Radius of the projection generating sphere */
  double	imrot;		/* Rotation angle of north pole */
  double	pa_north;	/* Position angle of north (0=horizontal) */
  double	pa_east;	/* Position angle of east (0=horizontal) */
  int		imflip;		/* If not 0, image is reflected around axis */
  int		prjcode;	/* projection code (-1-32) */
  int		latbase;	/* Latitude base 90 (NPA), 0 (LAT), -90 (SPA) */
  int		ncoeff1;	/* Number of x-axis plate fit coefficients */
  int		ncoeff2;	/* Number of y-axis plate fit coefficients */
  int		changesys;	/* 1 for FK4->FK5, 2 for FK5->FK4 */
  				/* 3 for FK4->galactic, 4 for FK5->galactic */
  int		printsys;	/* 1 to print coordinate system, else 0 */
  int		ndec;		/* Number of decimal places in PIX2WCST */
  int		degout;		/* 1 to always print degrees in PIX2WCST */
  int		tabsys;		/* 1 to put tab between RA & Dec, else 0 */
  int		rotmat;		/* 0 if CDELT, CROTA; 1 if CD */
  int		coorflip;	/* 0 if x=RA, y=Dec; 1 if x=Dec, y=RA */
  int		offscl;		/* 0 if OK, 1 if offscale */
  int		vimoswcson;		/* 1 if WCS is set, else 0 */
  int		naxes;		/* Number of axes in image */
  int		vimoswcsproj;	/* WCS_OLD: AIPS worldpos() and worldpix()
				   WCS_NEW: Mark Calabretta's WCSLIB subroutines
				   WCS_BEST: WCSLIB for all but CAR,COE,NCP
				   WCS_ALT:  AIPS for all but CAR,COE,NCP */
  int		linmode;	/* 0=system only, 1=units, 2=system+units */
  int		detector;	/* Instrument detector number */
  char		instrument[32];	/* Instrument name */
  char		ctype[4][9];	/* Values of CTYPEn keywords */
  char		c1type[9];	/*  1st coordinate type code:
					RA--, GLON, ELON */
  char		c2type[9];	/*  2nd coordinate type code:
					DEC-, GLAT, ELAT */
  char		ptype[9];	/*  projection type code:
				    SIN, TAN, ARC, NCP, GLS, MER, AIT, etc */
  char		units[4][32];	/* Units if LINEAR */
  char		radecsys[32];	/* Reference frame: FK4, FK4-NO-E, FK5, GAPPT*/
  char		radecout[32];	/* Output reference frame: FK4,FK5,GAL,ECL */
  char		radecin[32];	/* Input reference frame: FK4,FK5,GAL,ECL */
  double	eqin;		/* Input equinox (match sysin if 0.0) */
  double	eqout;		/* Output equinox (match sysout if 0.0) */
  int		sysin;		/* Input coordinate system code */
  int		sysvimoswcs;		/* WCS coordinate system code */
  int		sysout;		/* Output coordinate system code */
				/* WCS_B1950, WCS_J2000, WCS_GALACTIC,
				 * WCS_ECLIPTIC, WCS_LINEAR, WCS_ALTAZ  */
  char		center[32];	/* Center coordinates (with frame) */
  struct vimoswcsprm vimoswcsl;		/* WCSLIB main projection parameters */
  struct linprm lin;		/* WCSLIB image/pixel conversion parameters */
  struct celprm cel;		/* WCSLIB projection type */
  struct prjprm prj;		/* WCSLIB projection parameters */
  struct IRAFsurface *lngcor;	/* RA/longitude correction structure */
  struct IRAFsurface *latcor;	/* Dec/latitude correction structure */
  char *command_format[10];	/* WCS command formats */
				/* where %s is replaced by WCS coordinates */
				/* where %f is replaced by the image filename */
				/* where %x is replaced by image coordinates */
  double	ltm[4];		/* Image rotation matrix */
  double	ltv[2];		/* Image offset */
  int		idpix[2];	/* First pixel to use in image (x, y) */
  int		ndpix[2];	/* Number of pixels to use in image (x, y) */
  struct WorldCoor *vimoswcs;	/* WCS upon which this WCS depends */
};

/* Projections (1-26 are WCSLIB) */
#define VIMOSWCS_PIX -1	/* Pixel WCS */
#define VIMOSWCS_LIN  0	/* Linear projection */
#define VIMOSWCS_AZP  1	/* Zenithal/Azimuthal Perspective */
#define VIMOSWCS_TAN  2	/* Gnomonic = Tangent Plane */
#define VIMOSWCS_SIN  3	/* Orthographic/synthesis */
#define VIMOSWCS_STG  4	/* Stereographic */
#define VIMOSWCS_ARC  5	/* Zenithal/azimuthal equidistant */
#define VIMOSWCS_ZPN  6	/* Zenithal/azimuthal PolyNomial */
#define VIMOSWCS_ZEA  7	/* Zenithal/azimuthal Equal Area */
#define VIMOSWCS_AIR  8	/* Airy */
#define VIMOSWCS_CYP  9	/* CYlindrical Perspective */
#define VIMOSWCS_CAR 10	/* Cartesian */
#define VIMOSWCS_MER 11	/* Mercator */
#define VIMOSWCS_CEA 12	/* Cylindrical Equal Area */
#define VIMOSWCS_CPS 13	/* Conic PerSpective (COP) */
#define VIMOSWCS_COD 14	/* COnic equiDistant */
#define VIMOSWCS_COE 15	/* COnic Equal area */
#define VIMOSWCS_COO 16	/* COnic Orthomorphic */
#define VIMOSWCS_BON 17	/* Bonne */
#define VIMOSWCS_PCO 18	/* Polyconic */
#define VIMOSWCS_GLS 19	/* Sanson-Flamsteed (GLobal Sinusoidal) */
#define VIMOSWCS_PAR 20	/* Parabolic */
#define VIMOSWCS_AIT 21	/* Hammer-Aitoff */
#define VIMOSWCS_MOL 22	/* Mollweide */
#define VIMOSWCS_CSC 23	/* COBE quadrilateralized Spherical Cube */
#define VIMOSWCS_QSC 24	/* Quadrilateralized Spherical Cube */
#define VIMOSWCS_TSC 25	/* Tangential Spherical Cube */
#define VIMOSWCS_NCP 26	/* Special case of SIN */
#define VIMOSWCS_DSS 27	/* Digitized Sky Survey plate solution */
#define VIMOSWCS_PLT 28	/* Plate fit polynomials (SAO) */
#define VIMOSWCS_TNX 29	/* Gnomonic = Tangent Plane (NOAO with corrections) */

/* Coordinate systems */
#define VIMOSWCS_J2000	1	/* J2000(FK5) right ascension and declination */
#define VIMOSWCS_B1950	2	/* B1950(FK4) right ascension and declination */
#define VIMOSWCS_GALACTIC	3	/* Galactic longitude and latitude */
#define VIMOSWCS_ECLIPTIC	4	/* Ecliptic longitude and latitude */
#define VIMOSWCS_ALTAZ	5	/* Azimuth and altitude/elevation */
#define VIMOSWCS_LINEAR	6	/* Linear with optional units */
#define VIMOSWCS_NPOLE	7	/* Longitude and north polar angle */
#define VIMOSWCS_SPA		8	/* Longitude and south polar angle */
#define VIMOSWCS_PLANET	9	/* Longitude and latitude on planet */
#define VIMOSWCS_XY		10	/* X-Y Cartesian coordinates */

/* Method to use */
#define VIMOSWCS_BEST	0	/* Use best WCS projections */
#define VIMOSWCS_ALT		1	/* Use not best WCS projections */
#define VIMOSWCS_OLD		2	/* Use AIPS WCS projections */
#define VIMOSWCS_NEW		3	/* Use WCSLIB 2.5 WCS projections */

#ifndef PI
#define PI	3.141592653589793238462643
#endif

/* Conversions among hours of RA, degrees and radians. */
#define degrad(x)	((x)*PI/180.)
#define raddeg(x)	((x)*180./PI)
#define hrdeg(x)	((x)*15.)
#define deghr(x)	((x)/15.)
#define hrrad(x)	degrad(hrdeg(x))
#define radhr(x)	deghr(raddeg(x))

/* TNX surface fitting structure and flags */
struct IRAFsurface {
  double xrange;	/* 2. / (xmax - xmin), polynomials */
  double xmaxmin;	/* - (xmax + xmin) / 2., polynomials */
  double yrange;	/* 2. / (ymax - ymin), polynomials */
  double ymaxmin;	/* - (ymax + ymin) / 2., polynomials */
  int	 type;		/* type of curve to be fitted */
  int    xorder;	/* order of the fit in x */
  int    yorder;	/* order of the fit in y */
  int    xterms;	/* cross terms for polynomials */
  int    ncoeff;	/* total number of coefficients */
  double *coeff;	/* pointer to coefficient vector */
  double *xbasis;	/* pointer to basis functions (all x) */
  double *ybasis;	/* pointer to basis functions (all y) */
};

/* TNX permitted types of surfaces */
#define  TNX_CHEBYSHEV    1
#define  TNX_LEGENDRE     2
#define  TNX_POLYNOMIAL   3

/* TNX cross-terms flags */
#define	TNX_XNONE	0	/* no x-terms (old no) */
#define	TNX_XFULL	1	/* full x-terms (new yes) */
#define	TNX_XHALF	2	/* half x-terms (new) */

#ifdef __cplusplus /* allan: 28.4.98: added C++ prototypes */
extern "C" {
#endif

    /* WCS subroutines in wcs.c */
    struct WorldCoor *vimoswcsinit (char* hstring);
    struct WorldCoor *vimoswcsninit (
	char* hstring,	/* FITS header */
	int len);		/* Length of FITS header */
    struct WorldCoor *vimoswcsinitn (
	char* hstring,	/* FITS header */
	char* vimoswcsname);	/* WCS name */
    struct WorldCoor *vimoswcsninitn (
	char* hstring,	/* FITS header */
	int len,		/* Length of FITS header */
	char* vimoswcsname);	/* WCS name */
    struct WorldCoor *vimoswcsinitc (
	char* hstring,	/* FITS header */
	char vimoswcschar);	/* WCS character (A-Z) */
    struct WorldCoor *vimoswcsninitc (
	char* hstring,	/* FITS header */
	int len,		/* Length of FITS header */
	char vimoswcschar);	/* WCS character (A-Z) */
    void vimoswcsfree (
	struct WorldCoor *vimoswcs);	/* World coordinate system structure */

    int isvimoswcs(			/* Returns 1 if vimoswcs structure set, else 0 */
	struct WorldCoor *vimoswcs);	/* World coordinate system structure */
    int novimoswcs(			/* Returns 0 if vimoswcs structure set, else 1 */
	struct WorldCoor *vimoswcs);	/* World coordinate system structure */

    int pix2vimoswcst (
        struct WorldCoor *vimoswcs,  /* World coordinate system structure */
        double xpix, 
        double ypix,            /* Image coordinates in pixels */
        char   *vimoswcstring,       /* World coordinate string (returned) */
        int    lstr             /* Length of world coordinate string (returned) */
        );

    void pix2vimoswcs (
        struct WorldCoor *vimoswcs,  /* World coordinate system structure */
        double xpix, 
        double ypix,            /* Image coordinates in pixels */
        double *xpos, 
        double *ypos            /* RA and Dec in degrees (returned) */
        );

    void vimoswcsc2pix (
        struct WorldCoor *vimoswcs,  /* World coordinate system structure */
        double xpos,
        double ypos,            /* World coordinates in degrees */
	char *coorsys,		/* Coordinate system (B1950, J2000, etc) */
        double *xpix,
        double *ypix,           /* Image coordinates in pixels */
        int     *offscl);

    void vimoswcs2pix (
        struct WorldCoor *vimoswcs,  /* World coordinate system structure */
        double xpos,
        double ypos,            /* World coordinates in degrees */
        double *xpix,
        double *ypix,           /* Image coordinates in pixels */
        int     *offscl);

    double vimoswcsdist(		/* Compute angular distance between 2 sky positions */
	double ra0,		/* World coordinates in degrees */
	double dec0,
	double ra1,		/* World coordinates in degrees */
	double dec1);

    double vimoswcsdiff(		/* Compute angular distance between 2 sky positions */
	double ra0,		/* World coordinates in degrees */
	double dec0,
	double ra1,		/* World coordinates in degrees */
	double dec1);

    struct WorldCoor* vimoswcsxinit(
        double  cra,    /* Center right ascension in degrees */
        double  cdec,   /* Center declination in degrees */
        double  secpix, /* Number of arcseconds per pixel */
        double  xrpix,  /* Reference pixel X coordinate */
        double  yrpix,  /* Reference pixel X coordinate */
        int     nxpix,  /* Number of pixels along x-axis */
        int     nypix,  /* Number of pixels along y-axis */
        double  rotate, /* Rotation angle (clockwise positive) in degrees */
        int     equinox, /* Equinox of coordinates, 1950 and 2000 supported */
        double  epoch,  /* Epoch of coordinates, used for FK4/FK5 conversion
                         * no effect if 0 */
        char    *proj); /* Projection */

    struct WorldCoor* vimoswcskinit( /* set up WCS structure from keyword values */
	int     naxis1,		/* Number of pixels along x-axis */
	int     naxis2,		/* Number of pixels along y-axis */
	char    *ctype1,	/* FITS WCS projection for axis 1 */
	char    *ctype2,	/* FITS WCS projection for axis 2 */
	double  crpix1,		/* Reference pixel coordinates */
	double  crpix2,		/* Reference pixel coordinates */
	double  crval1,		/* Coordinate at reference pixel in degrees */
	double  crval2,		/* Coordinate at reference pixel in degrees */
	double  *cd,            /* Rotation matrix, used if not NULL */
	double  cdelt1,		/* scale in degrees/pixel, if cd is NULL */
	double  cdelt2,		/* scale in degrees/pixel, if cd is NULL */
	double  crota,          /* Rotation angle in degrees, if cd is NULL */
	int     equinox, /* Equinox of coordinates, 1950 and 2000 supported */
	double  epoch);  /* Epoch of coordinates, for FK4/FK5 conversion */

    void vimoswcsshift(		/* Change center of WCS */
        struct WorldCoor *vimoswcs,  /* World coordinate system structure */
        double  cra,            /* New center right ascension in degrees */
        double  cdec,           /* New center declination in degrees */
        char    *coorsys);      /* FK4 or FK5 coordinates (1950 or 2000) */

    void vimoswcsfull(
        struct WorldCoor *vimoswcs,  /* World coordinate system structure */
        double  *cra,           /* Right ascension of image center (deg) (returned) */
        double  *cdec,          /* Declination of image center (deg) (returned) */
        double  *width,         /* Width in degrees (returned) */
        double  *height);       /* Height in degrees (returned) */

    void vimoswcsrange(
        struct WorldCoor *vimoswcs,  /* World coordinate system structure */
        double  *ra1,           /* Min. right ascension of image (deg) (returned) */
        double  *ra2,           /* Max. right ascension of image (deg) (returned) */
        double  *dec1,          /* Min. declination of image (deg) (returned) */
        double  *dec2);         /* Max. declination of image (deg) (returned) */

    void setvimoswcserr(		/* Set WCS error message for later printing */
	char *errmsg);		/* Error mesage < 80 char */
    void vimoswcserr();		/* Print WCS error message to stderr */

    void setdefvimoswcs(		/* Set flag to use AIPS WCS instead of WCSLIB */
	int oldvimoswcs);		/* 1 for AIPS WCS subroutines, else WCSLIB */
    int getdefvimoswcs();		/* Return flag for AIPS WCS set by setdefvimoswcs */

    char *getradecsys(		/* Return name of image coordinate system */
        struct WorldCoor *vimoswcs);	/* World coordinate system structure */
	
    void vimoswcsoutinit(		/* Set output coordinate system for pix2vimoswcs */
        struct WorldCoor *vimoswcs,	/* World coordinate system structure */
	char *coorsys);		/* Coordinate system (B1950, J2000, etc) */

    char *getvimoswcsout(		/* Return current output coordinate system */
        struct WorldCoor *vimoswcs);	/* World coordinate system structure */

    void vimoswcsininit(		/* Set input coordinate system for vimoswcs2pix */
        struct WorldCoor *vimoswcs,	/* World coordinate system structure */
	char *coorsys);		/* Coordinate system (B1950, J2000, etc) */

    char *getvimoswcsin(		/* Return current input coordinate system */
        struct WorldCoor *vimoswcs);	/* World coordinate system structure */

    int setvimoswcsdeg(		/* Set WCS coordinate output format */
        struct WorldCoor *vimoswcs,	/* World coordinate system structure */
	int degout);		/* 1= degrees, 0= hh:mm:ss dd:mm:ss */

    int vimoswcsndec(		/* Set or get number of output decimal places */
        struct WorldCoor *vimoswcs,	/* World coordinate system structure */
	int ndec);		/* Number of decimal places in output string
				   if < 0, return current ndec unchanged */

    void setvimoswcslin(		/* Set pix2vimoswcst() mode for LINEAR coordinates */
        struct WorldCoor *vimoswcs,	/* World coordinate system structure */
	int mode);		/* 0: x y linear, 1: x units x units
				   2: x y linear units */

    int vimoswcszin(
	int izpix);		/* Set coordinate in third dimension (face) */

    int vimoswcszout (		/* Return coordinate in third dimension */
        struct WorldCoor *vimoswcs);	/* World coordinate system structure */

    void savevimoswcscoor(		/* Save output coordinate system */
	char *vimoswcscoor);		/* coordinate system (J2000, B1950, galactic) */
    char *getvimoswcscoor();		/* Return output coordinate system */
    void savevimoswcscom(		/* Save WCS shell command */
            int i,
            char *vimoswcscom);		/* Shell command using output WCS string */
    char *getvimoswcscom();		/* Return WCS shell command */
    void setvimoswcsfile(		/* Set filename for WCS error message */
            char *filename);	/* FITS or IRAF file name */


    /* Coordinate conversion subroutines in vimoswcscon.c */
    void vimoswcsconv(	/* Convert between coordinate systems and equinoxes */
	int sys1,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
	int sys2,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
	double eq1,	/* Input equinox (default of sys1 if 0.0) */
	double eq2,	/* Output equinox (default of sys2 if 0.0) */
	double ep1,	/* Input Besselian epoch in years */
	double ep2,	/* Output Besselian epoch in years */
	double *dtheta,	/* Longitude or right ascension in degrees
			   Input in sys1, returned in sys2 */
	double *dphi,	/* Latitude or declination in degrees
			   Input in sys1, returned in sys2 */
	double *ptheta,	/* Longitude or right ascension proper motion in deg/year
			   Input in sys1, returned in sys2 */
	double *pphi,	/* Latitude or declination proper motion in deg/year */
	double *px,	/* Parallax in arcseconds */
	double *rv);	/* Radial velocity in km/sec */
    void vimoswcsconp(	/* Convert between coordinate systems and equinoxes */
	int sys1,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
	int sys2,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
	double eq1,	/* Input equinox (default of sys1 if 0.0) */
	double eq2,	/* Output equinox (default of sys2 if 0.0) */
	double ep1,	/* Input Besselian epoch in years */
	double ep2,	/* Output Besselian epoch in years */
	double *dtheta,	/* Longitude or right ascension in degrees
			   Input in sys1, returned in sys2 */
	double *dphi,	/* Latitude or declination in degrees
			   Input in sys1, returned in sys2 */
	double *ptheta,	/* Longitude or right ascension proper motion in degrees/year
			   Input in sys1, returned in sys2 */
	double *pphi);	/* Latitude or declination proper motion in degrees/year
			   Input in sys1, returned in sys2 */
    void vimoswcscon(	/* Convert between coordinate systems and equinoxes */
	int sys1,	/* Input coordinate system (J2000, B1950, ECLIPTIC, GALACTIC */
	int sys2,	/* Output coordinate system (J2000, B1950, ECLIPTIC, G ALACTIC */
	double eq1,	/* Input equinox (default of sys1 if 0.0) */
	double eq2,	/* Output equinox (default of sys2 if 0.0) */
	double *dtheta,	/* Longitude or right ascension in degrees
			   Input in sys1, returned in sys2 */
	double *dphi,	/* Latitude or declination in degrees
			   Input in sys1, returned in sys2 */
	double epoch);	/* Besselian epoch in years */

    int vimoswcscsys(	/* Return code for coordinate system in string */
	char *coorsys);	 /* Coordinate system (B1950, J2000, etc) */

    double vimoswcsceq (	/* Set equinox from string (return 0.0 if not obvious) */
	char *vimoswcstring);  /* Coordinate system (B1950, J2000, etc) */

    void vimoswcscstr (	/* Set coordinate system type string from system and equinox */
	char   *cstr,	 /* Coordinate system string (returned) */
	int    sysvimoswcs,	/* Coordinate system code */
	double equinox,	/* Equinox of coordinate system */
	double epoch);	/* Epoch of coordinate system */

    void vimoswcsdeltset ( /* set scaling and rotation from CDELTs and CROTA2 */
    struct WorldCoor *vimoswcs, /* World coordinate system structure */
    double cdelt1, /* degrees/pixel in first axis (or both axes) */
    double cdelt2, /* degrees/pixel in second axis if nonzero */
    double crota); /* Rotation counterclockwise in degrees */

    void    vimoswcserr();
    int getdefvimoswcs ();
    int vimoswcstype (
            struct WorldCoor *vimoswcs, /* World coordinate system structure */
            char    *ctype1,    /* FITS WCS projection for axis 1 */
            char    *ctype2);    /* FITS WCS projection for axis 2 */
    int    vimoswcsreset (
            struct WorldCoor *vimoswcs,     /* World coordinate system data structure */
            double crpix1,
            double crpix2,      /* Reference pixel coordinates */
            double crval1,
            double crval2,      /* Coordinates at reference pixel in degrees */
            double cdelt1, 
            double cdelt2,      /* scale in degrees/pixel, ignored if cd is not NULL */
            double crota,           /* Rotation angle in degrees, ignored if cd is not NULL */
            double *cd);         /* Rotation matrix, used if not NULL */

    void vimoswcseqset (
    struct WorldCoor *vimoswcs,     /* World coordinate system data structure */
    double equinox);         /* Desired equinox as fractional year */

    void vimoswcscdset (
    struct WorldCoor *vimoswcs, /* World coordinate system structure */
    double *cd);         /* CD matrix, ignored if NULL */


#ifdef __cplusplus /* allan: 28.4.98: added C++ prototypes */
};
#endif

#ifndef __cplusplus
/* WCS subroutines in vimoswcs.c */
int vimoswcstype();		/* Set projection type from header CTYPEs */
void vimoswcscdset();	/* Set scaling and rotation from CD matrix */
void vimoswcspcset();	/* set scaling and rotation from CDELTs and PC matrix */
void vimoswcscent();		/* Print the image center and size in WCS units */
void vimoswcssize();		/* Return RA and Dec of image center, size in RA and Dec */
void vimoswcscominit();	/* Initialize catalog search command set by -vimoswcscom */
void vimoswcscom();		/* Execute catalog search command set by -vimoswcscom */
int vimoswcsreset();		/* Change WCS using arguments */
void vimoswcseqset();	/* Change equinox of reference pixel coordinates in WCS */
void setvimoswcscom();	/* Set WCS shell commands from stored values */
void freevimoswcscom();	/* Free memory used to store WCS shell commands */

#endif

/* Oct 26 1994	New file
 * Dec 21 1994	Add rotation matrix
 * Dec 22 1994	Add flag for coordinate reversal

 * Mar  6 1995	Add parameters for Digital Sky Survey plate fit
 * Jun  8 1995	Add parameters for coordinate system change
 * Jun 21 1995	Add parameter for plate scale
 * Jul  6 1995	Add parameter to note whether WCS is set
 * Aug  8 1995	Add parameter to note whether to print coordinate system
 * Oct 16 1995	Add parameters to save image dimensions and center coordinates

 * Feb 15 1996	Add coordinate conversion functions
 * Feb 20 1996	Add flag for tab tables
 * Apr 26 1996	Add epoch of positions (actual date of image)
 * Jul  5 1996	Add subroutine declarations
 * Jul 19 1996	Add WCSFULL declaration
 * Aug  5 1996	Add WCSNINIT to initialize WCS for non-terminated header
 * Oct 31 1996	Add DCnn inverse rotation matrix
 * Nov  1 1996	Add NDEC number of decimal places in output
 *
 * May 22 1997	Change range of pcode from 1-8 to -1-8 for linear transform
 * Sep 12 1997	Add chip rotation MROT, XMPIX, YMPIX
 *
 * Jan  7 1998	Add INSTRUME and DETECTOR for HST metric correction
 * Jan 16 1998	Add Mark Calabretta's WCSLIB data structures
 * Jan 16 1998	Add LONGPOLE, LATPOLE, and PROJP constants for Calabretta
 * Jan 22 1998	Add ctype[], crpix[], crval[], and cdelt[] for Calabretta
 * Jan 23 1998	Change wcsset() to wcsxinit() and pcode to prjcode
 * Jan 23 1998	Define projection type flags
 * Jan 26 1998	Remove chip rotation
 * Jan 26 1998	Add chip correction polynomial
 * Feb  3 1998	Add number of coefficients for residual fit
 * Feb  5 1998	Make cd and dc matrices vectors, not individual elements
 * Feb 19 1998	Add projection names
 * Feb 23 1998	Add TNX projection from NOAO
 * Mar  3 1998	Add NOAO plate fit and residual fit
 * Mar 12 1998	Add variables for TNX correction surface
 * Mar 23 1998	Add PLT plate fit polynomial projection; reassign DSS
 * Mar 23 1998	Drop plate_fit flag from structure
 * Mar 25 1998	Add npcoeff to wcs structure for new plate fit WCS
 * Apr  7 1998	Change amd_i_coeff to i_coeff
 * Apr  8 1998	Add wcseqset() and wcsreset() subroutine declarations
 * Apr 10 1998	Rearrange order of nonstandard WCS types
 * Apr 13 1998	Add setdefwcs() subroutine declaration
 * Apr 14 1998	Add coordinate systems and wcscoor()
 * Apr 24 1998	Add units
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * Apr 28 1998	Change projection flags to WCS_*
 * Apr 28 1998	Add wcsc2pix()
 * May  7 1998	Add C++ declarations
 * May 13 1998	Add eqin and eqout for conversions to and from equinoxes
 * May 14 1998	Add declarations for coordinate conversion subroutines
 * May 27 1998	Add blsearch()
 * May 27 1998	Change linear projection back to WCS_LIN from WCS_LPR
 * May 27 1998	Move hget.c and hput.c C++ declarations to fitshead.h
 * May 27 1998	Include fitshead.h
 * May 29 1998	Add wcskinit()
 * Jun  1 1998	Add wcserr()
 * Jun 11 1998	Add initialization support subroutines
 * Jun 18 1998	Add wcspcset()
 * Jun 25 1998	Add wcsndec()
 * Jul  6 1998	Add wcszin() and wcszout() to use third dimension of images
 * Jul  7 1998	Change setdegout() to setwcsdeg(); setlinmode() to setwcslin()
 * Jul 17 1998	Add savewcscoor(), getwcscoor(), savewcscom(), and getwcscom()
 * Aug 14 1998	Add freewcscom(), setwcscom(), and multiple WCS commands
 * Sep  3 1998	Add pa_north, pa_east, imrot and imflip to wcs structure
 * Sep 14 1998	Add latbase for AXAF North Polar angle (NPOL not LAT-)
 * Sep 16 1998	Make WCS_system start at 1; add NPOLE
 * Sep 17 1998	Add wcscstr()
 * Sep 21 1998	Add wcsconp() to convert proper motions, too.
 * Dec  2 1998	Add WCS type for planet surface

 * Jan 20 1999	Add declaration of wcsfree()
 * Jun 16 1999	Add declaration of wcsrange()
 * Oct 21 1999	Add declaration of setwcsfile()
 *
 * Jan 28 2000	Add flags for choice of WCS projection subroutines
 * Jun 26 2000	Add XY coordinate system
 * Nov  2 2000	Add wcsconv() to convert coordinates when parallax or rv known
 *
 * Jan 17 2001	Add idpix and ndpix for trim section, ltm for readout rotation
 * Jan 31 2001	Add wcsinitn(), wcsninitn(), wcsinitc(), and wcsninitc()
 * Feb 20 2001	Add wcs->wcs to main data structure
 * Mar 20 2001	Close unclosed comment in wcsconv() argument list
 */
#endif
