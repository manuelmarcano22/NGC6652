 /*
 				prefs.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Keywords for the configuration file.
*
*	Last modify:	25/08/99
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define	MAXLIST		32		/* max. nb of list members */

/* NOTES:
One must have:	MAXLIST >= 1 (preferably >= 16!)
*/

/*-------------------------------- initialization ---------------------------*/
int	idummy;

pkeystruct key[] =
 {
  {"BADPIXEL_FILTER", P_BOOL, &prefs.badpix_flag},
  {"BADPIXEL_NMAX", P_INT, &prefs.badpix_nmax, 0,100000000},
  {"CHECKIMAGE_NAME", P_STRINGLIST, prefs.check_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncheck_name},
  {"CHECKIMAGE_TYPE", P_KEYLIST, prefs.check_type, 0,0, 0.0,0.0,
   {"NONE", "CHI", "PROTOTYPES", "RESIDUALS", "RAWDATA", "SAMPLES",
	"SNAPSHOTS", "WEIGHTS", "PC_CONVOLVED"},
   0, MAXCHECK, &prefs.ncheck_type},
  {"CONTEXT_KEYS", P_STRINGLIST, prefs.context_name, 0,0,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_name},
  {"CONTEXT_GROUPS", P_INTLIST, prefs.context_group, 1,MAXCONTEXT,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_group},
  {"GROUP_DEGREES", P_INTLIST, prefs.group_deg, 1,32,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ngroup_deg},
  {"PC_INCLUDE", P_BOOL, &prefs.pc_flag},
  {"PC_NAME", P_STRING, prefs.pc_name},
  {"PC_NPC",  P_INT, &prefs.pc_npc, 0,1000000},
  {"PSF_ACCURACY", P_FLOAT, &prefs.prof_accuracy, 0,0, 0.0,1.0},
  {"PSF_AUTOSELECT", P_BOOL, &prefs.autoselect_flag},
  {"PSF_NSUPER", P_INT, &prefs.nsuper,0,1000000000},
  {"PSF_FWHMRANGE", P_FLOATLIST, prefs.fwhmrange, 0,0, 0.0,1e3, {""},
     2,2, &prefs.nfwhmrange},
  {"PSF_MAXELONG", P_FLOAT, &prefs.maxelong, 0,0, 1.0, BIG},
  {"PSF_MINSN", P_FLOAT, &prefs.minsn, 0,0, 1e-6,1e15},
  {"PSF_NAME", P_STRING, prefs.psf_name},
  {"PSF_RECENTER", P_BOOL, &prefs.recenter_flag},
  {"PSF_SAMPLING", P_FLOAT, &prefs.psf_step, 1,1024, 0.0,1.0e3},
  {"PSF_SIZE", P_INTLIST, prefs.retisize, 1,1024, 0.0,0.0, {""},
     1,2, &prefs.nretisize},
  {"PSF_VARIABILITY", P_FLOAT, &prefs.maxvar, 0,0, 0.0, BIG},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL","FULL",""}},
  {""}
 };

char			keylist[sizeof(key)/sizeof(pkeystruct)][16];
static const char	notokstr[] = {" \t=,;\n\r\""};

char *(default_prefs[]) =
 {
  "BADPIXEL_FILTER	N",
  "BADPIXEL_NMAX	0",
  "PC_INCLUDE		N",
  "PC_NAME		default.pc",
  "PC_NPC		0",
  "PSF_ACCURACY		1e-2",
  "PSF_AUTOSELECT	Y",
  "PSF_FWHMRANGE	0.5,10.0",
  "PSF_MAXELONG		2.0",
  "PSF_MINSN 		20",
  "PSF_NSUPER		64",
  "PSF_RECENTER		N",
  "PSF_SAMPLING		0.0",
  "PSF_VARIABILITY     	0.5",
  "VERBOSE_TYPE        NORMAL",
  ""
 };
