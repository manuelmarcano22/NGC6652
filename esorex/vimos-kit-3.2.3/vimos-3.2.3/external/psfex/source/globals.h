 /*
 				globals.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SOMfit
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	global declarations.
*
*	Last modify:	01/03/99
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------- miscellaneous variables ---------------------------*/
prefstruct		prefs;

/*------------------------------- functions ---------------------------------*/
extern void		error(int, char *, char *),
			makesom(char **, int),
			readprefs(char *filename, char **argkey,
				char **argval, int narg),
			swapbytes(void *, int, int),
			useprefs(void),
			warning(char *, char *);

extern int		copyvignet_center(float *, int, int, float *,
				int, int, float, float, float, float, float);


extern float		hmedian(float *ra, int n);
