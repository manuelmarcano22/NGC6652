 /*
 				makeit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	Learning and testing.
*
*	Last modify:	13/10/99
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/


#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fitscat.h"
#include	"sample.h"
#include	"psf.h"
#include	"vignet.h"
#include	"check.h"

/********************************** makeit ***********************************/
/*
*/
void	makeit(char **incatnames, int ncat)

  {
   setstruct		*set;
   psfstruct		*psf;
   pcstruct		*pc, *pcc, *pco;
   static char		str[MAXCHAR];
   float		psfstep;
   int			i;

/* Load all the samples */
  NFPRINTF(OUTPUT,"Loading samples...");
  set = load_samples(incatnames, ncat);

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "%d samples loaded\n", set->nsample);

  if (!set->nsample)
    error(EXIT_FAILURE, "*Error*: No appropriate source found!!","");

  psfstep = (float)(prefs.psf_step? prefs.psf_step
				: (set->fwhm/2.35)*(1.0-1.0/INTERPFAC));

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "PSF sampled each %.2f pixel%s (%s)\n",
			psfstep,
			psfstep>=2.0?"s":"",
			prefs.psf_step?"manual":"automatic");

/* Init the PSF */
  NFPRINTF(OUTPUT,"Initializing PSF modules...");
  psf = psf_init(prefs.context_name, prefs.context_group, prefs.ncontext_name,
		prefs.group_deg, prefs.ngroup_deg,
		set->retisize[0], set->retisize[1],
		psfstep,
		set->nsample);

/* Make the basic PSF-model (1st pass) */
  NFPRINTF(OUTPUT,"Modeling the PSF.");
  psf_make(psf, set);

/* Remove bad PSF candidates */
  psf_clean(psf, set, 1);
  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "%d samples accepted\n", set->nsample);

/* Make the basic PSF-model (2nd pass) */
  NFPRINTF(OUTPUT,"Modeling the PSF.");
  psf_make(psf, set);

/* Remove bad PSF candidates */
  psf_clean(psf, set, 1);
  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "%d samples accepted\n", set->nsample);

/* Refine the PSF-model */
  psf_refine(psf, set, prefs.nsuper);

/* Just check the Chi2 */
  psf_clean(psf, set, 0);
  NFPRINTF(OUTPUT, "");

/* Load the PCs */
  if (prefs.pc_flag)
    {
    NFPRINTF(OUTPUT,"Including principal components...");
    pc = pc_load(prefs.pc_name);
    pcc = pc_convolve(pc, psf);
    pco = pc_orthogon(pc, pcc, psf->pixstep);
    }
  else
    pc = pcc = pco = NULL;

/* Save result */
  NFPRINTF(OUTPUT,"Saving the PSF description...");

  psf_save(psf, pco, pc, prefs.psf_name);

/* Save "Check-images" */
  for (i=0; i<prefs.ncheck_type; i++)
    if (prefs.check_type[i])
      {
      sprintf(str, "Saving CHECK-image #%d...", i+1);
      NFPRINTF(OUTPUT, str);
      psf_writecheck(psf, pco, set, prefs.check_name[i], prefs.check_type[i]);
      }

/* Free memory */
  end_set(set);
  psf_end(psf);
  if (prefs.pc_flag)
    {
    pc_end(pc);
    pc_end(pcc);
    pc_end(pco);
    }

  return;
  }

