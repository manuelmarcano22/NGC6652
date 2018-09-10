/*
 				types.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN, (IAP, Leiden observatory & ESO)
*
*	Contents:	global type definitions.
*
*	Last modify:	25/08/99
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

typedef		float		PIXTYPE;
typedef enum {PSF_NONE, PSF_CHI, PSF_PROTO, PSF_RESIDUALS, PSF_RAWDATA,
		PSF_SAMPLES, PSF_SNAPSHOTS, PSF_WEIGHTS, PSF_PCPROTO}
	checkenum;

/*------------------------------- preferences -------------------------------*/
typedef struct
  {
  char		prefs_name[MAXCHAR];		/* prefs filename */
  char		psf_name[MAXCHAR];		/* PSF filename */
  int		retisize[2], nretisize;		/* Retina size */
  double	minsn;				/* Minimum S/N for patterns */
  double	maxelong;			/* Maximum A/B for patterns */
  double	maxvar;				/* Maximum FWHM variability */
  double	fwhmrange[2];			/* Allowed FWHM range */
  int		nfwhmrange;	       		/* nb of params */
  double	prof_accuracy;			/* Required PSF accuracy */
  double	psf_step;			/* Oversampling (pixels) */
  int		nsuper;				/* nb of supersampled pixels */
  int		autoselect_flag;		/* Auto. select FWHMs ? */
  int		recenter_flag;			/* Recenter PSF-candidates? */
  checkenum	check_type[MAXCHECK];		/* check-image types */
  int		ncheck_type;			/* nb of params */
  char		*(check_name[MAXCHECK]);	/* check-image names */
  int		ncheck_name;			/* nb of params */
  char		*(context_name[MAXCONTEXT]);	/* Names of context-keys */
  int		ncontext_name;			/* nb of params */
  int		context_group[MAXCONTEXT];	/* Context group */
  int		ncontext_group;			/* nb of params */
  int		group_deg[MAXCONTEXT];		/* Degree for each group */
  int		ngroup_deg;			/* nb of params */
  int		badpix_flag;			/* Filter bad pixels? */
  int		badpix_nmax;			/* Max number of bad pixels */
  int		pc_flag;			/* Include PCs? */
  char		pc_name[MAXCHAR];		/* PC filename */
  int		pc_npc;				/* Max. number of PCs */
  enum {QUIET, NORM, FULL}	verbose_type;	/* How much it displays info */
  }	prefstruct;

