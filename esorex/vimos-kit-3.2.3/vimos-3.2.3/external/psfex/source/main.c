 /*
 				main.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	parsing and main loop.
*
*	Last modify:	02/06/98
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"

#define		SYNTAX \
"psfex input_Catalog(s) [-c <config_file>] [-<keyword> <value>]\n\n"

/********************************** main ************************************/

main(int argc, char *argv[])
  {
   int		a,ab,na, narg;
   char		**argkey, **argval;

  if (argc<2)
    {
    fprintf(OUTPUT, "\n		%s  Version %s\n", BANNER, VERSION);
    fprintf(OUTPUT, "\nFor information, please contact: %s\n", COPYRIGHT);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }

  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/* Default parameters */
  ab=na=0;
  narg = 0;
  strcpy(prefs.prefs_name, "default.psfex");
  prefs.verbose_type = NORM;

  for (a=1; a<argc; a++)
    {
    if (argv[a][0] == '-' && a<(argc-1))
      {
      if (strlen(argv[a])<3)
        switch((int)tolower((int)argv[a][1]))
          {
          case 'c':
            strcpy(prefs.prefs_name, argv[++a]);
            break;
          default:
            error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);
          }
      else
        {
        argkey[narg] = &argv[a][1];
        argval[narg++] = argv[++a];
        }       
      }
    else
      {
      for(ab = a; (a<argc) && (*argv[a]!='-'); a++);
      na = (a--) - ab;
      }
    }

  if (!na)
    error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);

  readprefs(prefs.prefs_name, argkey, argval, narg);
  useprefs();
  free(argkey);
  free(argval);


  makeit(argv+ab, na);

  NFPRINTF(OUTPUT, "All done");
  NPRINTF(OUTPUT, "\n");

  exit(EXIT_SUCCESS);
  }
