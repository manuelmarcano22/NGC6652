  /*
 				basis.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	Basis functions for the PSF model.
*
*	Last modify:	26/02/99
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"poly.h"

/****** basis_alard ***********************************************************
PROTO   double *basis_alard(int w, int h, int deg, int ngauss, double fwhm,
                            int *nbasis)
PURPOSE Compute Alard & Lupton's (1998) image basis.
INPUT   Width of vignets,
        Height of vignets,
        Max. degree of polynom,
        Number of gaussians.
        Pointer to an int.
OUTPUT  2D*ncoeff*n2 array. The number of basis images is returned in nbasis.
NOTES   Uses functions in poly.c.
AUTHOR  E. Bertin (IAP)
VERSION 26/02/99
 ***/
double	*basis_alard(int w, int h, int deg, int ngauss, double fwhm,
			int *nbasis)
  {
   polystruct	*poly;
   double	x[2], *basis,*basist,*basist0, *pbasis,*pbasist,
		sigma,alpha,expo,dval;
   int		c,g,n, ix,iy, hw,hh, ncoeff,npix;

/* Prepare a POLYnom structure (function of x and y) */
  g = 0;	/* Only 1 group, with label 0 */
  poly = poly_init(&g, 2, &deg, 1);

  pbasis = poly->basis;
  ncoeff = poly->ncoeff;
  npix = w*h;
  *nbasis = ncoeff*ngauss;
/* Allocate memory to store the basis functions */
  QMALLOC(basis, double, npix**nbasis);

  alpha = 1/(2.0*(fwhm/2.35)*(fwhm/2.35));
  hw = w/2;
  hh = h/2;

  basist0 = basis;
  for (iy=-hh; iy<hh; iy++)
    for (ix=-hw; ix<hw; ix++)
      {
      basist = basist0++;
      poly_func(poly, x);
      for (n=0; n<ngauss; n++)
        {
        if ((expo=alpha/(n+1.0)*(ix*ix+iy*iy))<70.0)
        dval = exp(-expo);
        for (pbasist=pbasis, c=ncoeff; c--; basist+=npix)
          *basist = dval**(pbasist++);
        }
      }

  poly_end(poly);

  return basis;
  }

