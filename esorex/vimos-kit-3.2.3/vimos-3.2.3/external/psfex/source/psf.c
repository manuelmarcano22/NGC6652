  /*
 				psf.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	Stuff related to building the PSF.
*
*	Last modify:	13/10/99
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
#include	"fitscat.h"
#include	"sample.h"
#include	"poly.h"
#include	"psf.h"
#include	"vignet.h"

/********************************* psf_clean *********************************/
/*
Filter PSF candidates.
*/
#define	EPS	(1e-4)  /* a small number */
void	psf_clean(psfstruct *psf, setstruct *set, int clean_flag)
  {
   samplestruct	*sample;
   double	chi2,chimean,chivar,chisig,chisig1,chival, locut,hicut;
   float	*chi, *chit,*chit2,
		chimed, chi2max;
   int		i, n, nsample;

/* First compute residuals for each sample (chi^2) */
  NFPRINTF(OUTPUT,"Computing residuals...");
  psf_makeresi(psf, set, prefs.recenter_flag);

/* Store the chi's (sqrt(chi2) pdf close to gaussian) */
  NFPRINTF(OUTPUT,"Computing Chi2 statistics...");
  nsample = set->nsample;
  QMALLOC(chi, float, nsample);
  chit = chi;
  for (sample=set->sample, n=nsample; n--; sample++)
    *(chit++) = (float)sqrt(sample->chi2);
/* Produce k-sigma-clipped statistiscs */
  locut = -BIG;
  hicut = BIG;
  chisig = BIG;
  chisig1 = 1.0;
  for (i=clean_flag?100:1; i-- && chisig>=0.1 && fabs(chisig/chisig1-1.0)>EPS;)
    {
    chisig1 = chisig;
    chimed = hmedian(chi, nsample);
    chimean = chivar = 0.0;
    chit2 = chit = chi;
    for (n=nsample; n--;)
      {
      chival = *(chit++);
      if (chival>locut && chival<hicut)
        {
        chimean += (*(chit2++) = chival);
        chivar += chival*chival;
        }
      else
        nsample--;
      }

    chimean /= (double)nsample;
    chisig = sqrt((chivar-chimean*chimean*nsample)/(nsample-(nsample>1?1:0)));
    locut = chimed - 3.0*chisig;
    hicut = chimed + 3.0*chisig;
    }

  free(chi);

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "<Chi2/dof> = %.3f\n",chivar/(nsample-(nsample>1?1:0)));

/* Clip outliers */
  if (clean_flag)
    {
    NFPRINTF(OUTPUT,"Filtering PSF-candidates...");
    chi2max = (float)hicut;
    chi2max *= chi2max;
    nsample=set->nsample;
    for (sample=set->sample, n=0; n<nsample;)
    if ((sample++)->chi2>chi2max)
      {
      sample=remove_sample(set, n);
      nsample--;
      }
    else
      n++;
    }

  return;
  }


/****** psf_init *************************************************************
PROTO   psfstruct *psf_init(int *dim, int ndim)
PURPOSE Allocate and initialize a PSF structure.
INPUT   1D array of degrees of the polynom,
        array of char pointers to the context names,
        number of dimensions.
OUTPUT  psfstruct pointer.
NOTES   The maximum degrees and number of dimensions allowed are set in poly.h.
AUTHOR  E. Bertin (IAP)
VERSION 04/03/99
 ***/
psfstruct	*psf_init(char **names, int *group, int ndim,
			int *dim, int ngroup,
			int wpsf, int hpsf, float psfstep, int nsample)
  {
   psfstruct	*psf;
   static char	str[MAXCHAR];
   char		**names2, **names2t;
   int		*group2, *dim2,
		d, ndim2,ngroup2, npix;

/* Allocate memory for the PSF structure itself */
  QCALLOC(psf, psfstruct, 1);
  psf->dim = PSF_NMASKDIM;	/* This is constant */
  QMALLOC(psf->size, int, psf->dim);

/* The polynom */
  names2 = NULL;
  group2 = dim2 = NULL;
  if (ndim2=ndim)
    {
    QMEMCPY(names, names2, char *, ndim);
    QMEMCPY(group, group2, int, ndim);
    }
  if (ngroup2=ngroup)
    QMEMCPY(dim, dim2, int, ngroup);
  psf->poly = poly_init(group2, ndim2, dim2, ngroup2);

/* Compute the maximum advised number of degrees of freedom */
  while ((int)(psf->poly->ncoeff/(psfstep*psfstep*PSF_FREEDFACTOR)+0.499)
	   >nsample)
    {
    poly_end(psf->poly);
    if (ngroup2)
      {
/*---- If still too many degrees of freedom, try to lower degrees */
      d=ngroup2%10;
      sprintf(str, "%d%s", ngroup2, d==1?"st":(d==2?"nd":(d==3?"rd":"th")));
      if (!(--dim2[ngroup2-1]))
        {
/*---- If degree is 0, just remove all the group components */
        for (d=0; d<ndim2; d++)
          if (group2[d]==ngroup2 && d!=(--ndim2))
            {
            names2[d]=names2[ndim2];
            group2[d]=group2[ndim2];
            }
        ngroup2--;
        warning(str, " context group removed (not enough samples)");
        }
      else
        warning(str, " context group-degree lowered (not enough samples)");
      psf->poly = poly_init(group2, ndim2, dim2, ngroup2);
      }
    else
      error(EXIT_FAILURE, "*Error*: Not enough sources!!","");
    }

  psf->pixstep = psfstep;
  psf->npix = psf->size[0] = wpsf;
  psf->npix *= (psf->size[1] = hpsf);
  psf->npix *= (psf->size[2] = psf->poly->ncoeff);
  QMALLOC(psf->comp, float, psf->npix);
  npix = psf->size[0]*psf->size[1];
  QMALLOC(psf->loc, float, npix);
  QMALLOC(psf->resi, float, npix);

/* Context arrays */
  if (ndim2)
    {
    QMALLOC(psf->contextoffset, double, ndim2);
    QMALLOC(psf->contextscale, double, ndim2);
    QMALLOC(psf->contextname, char *, ndim2);
    for (names2t=names2, d=0; d<ndim2; d++)
      {
      QMALLOC(psf->contextname[d], char, 80);
      strcpy(psf->contextname[d], *(names2t++));
      }
    }


/* Free temporary arrays */
  if (ndim)
    {
    free(names2);
    free(group2);
    free(dim2);
    }

 return psf;
  }


/****** psf_end **************************************************************
PROTO   void psf_end(psfstruct *psf)
PURPOSE Free a PSF structure and everything it contains.
INPUT   psfstruct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 04/07/98
 ***/
void	psf_end(psfstruct *psf)
  {
   int	d, ndim;

  ndim = psf->poly->ndim;
  for (d=0; d<ndim; d++)
    free(psf->contextname[d]);
  free(psf->contextname);
  poly_end(psf->poly);
  free(psf->comp);
  free(psf->loc);
  free(psf->resi);
  free(psf->size);
  free(psf);

  return;
  }


/******************************* psf_build **********************************/
/*
Build the local PSF (function of "coordinates").
*/
void	psf_build(psfstruct *psf, double *pos)
  {
   double	*basis;
   float	*ppc, *pl, fac;
   int		n,p, npix;

  npix = psf->size[0]*psf->size[1];
/* Reset the Local PSF mask */
  memset(psf->loc, 0, npix*sizeof(float));

  poly_func(psf->poly, pos);
  basis = psf->poly->basis;

  ppc = psf->comp;
/* Sum each component */
  for (n = (psf->dim>2?psf->size[2]:1); n--;)
    {
    pl = psf->loc;
    fac = (float)*(basis++);
    for (p=npix; p--;)
      *(pl++) +=  fac**(ppc++);
    }

  return;
  }


/****************************** psf_makeresi *********************************/
/*
Compute PSF residuals (chi2).
*/
void	psf_makeresi(psfstruct *psf, setstruct *set, int centflag)
  {
   samplestruct		*sample;
   static double	pos[MAXCONTEXT], amat[9], bmat[3];
   double		*dresi, *dresit, *amatt,
			*cvigx,*cvigxt, *cvigy,*cvigyt,
			nm1, chi2, dx,dy, ddx,ddy, dval,dvalx,dvaly,dwval,
			radmin2,radmax2, hcw,hch, yb, mx2,my2,mxy;
   float		*vigresi, *vig, *vigw, *fresi,*fresit,
			*cbasis,*cbasist, *cdata,*cdatat, *cvigw,*cvigwt,
			norm, fval, vigstep;
   int			i,j,n,x,y, ndim,npix,nsample, cw,ch,ncpix, okflag;

  vigstep = 1/psf->pixstep;
  nsample = set->nsample;
  npix = set->vigsize[0]*set->vigsize[1];
  ndim = psf->poly->ndim;
  QCALLOC(dresi, double, npix);

  if (centflag)
    {
/*-- Compute Centering sub-vignet size (containing most of the signal) */
    cw=ch=(int)(2*set->fwhm+1.0);
    if (cw>set->vigsize[0])
      cw=set->vigsize[0];
    if (ch>set->vigsize[1])
      ch=set->vigsize[1];
/*-- Allocate memory for the sub-vignet */
    ncpix = cw*ch;
    QMALLOC(cdata, float, ncpix);
    QMALLOC(cbasis, float, ncpix);
    QMALLOC(cvigw, float, ncpix);
    QMALLOC(cvigx, double, ncpix);
    QMALLOC(cvigy, double, ncpix);
/*-- Initialize gradient image */
    hcw = (double)(cw/2);
    hch = (double)(ch/2);
    cvigxt = cvigx;
    cvigyt = cvigy;
    for (y=0; y<ch; y++)
      {
      yb = y-hch;
      for (x=0; x<cw; x++)
        {
        *(cvigxt++) = x-hcw;
        *(cvigyt++) = yb;
        }
      }
 /*-- Set convergence boundaries */
    radmin2 = PSF_MINSHIFT*PSF_MINSHIFT;
    radmax2 = PSF_MAXSHIFT*PSF_MAXSHIFT;
    okflag = 0;
    }

/* Compute the chi2 */
  for (sample=set->sample, n=nsample; n--; sample++)
    {
/*-- Build the local PSF */
    for (i=0; i<ndim; i++)
      pos[i] = (sample->context[i]-set->contextoffset[i])
		/set->contextscale[i];
    psf_build(psf, pos);

/*-- Delta-x and Delta-y in vignet-pixel units */
    dx = sample->dx;
    dy = sample->dy;

    if (centflag)
      {
/*---- Copy the data into the sub-vignet */
      vignet_copy(sample->vig, set->vigsize[0], set->vigsize[1],
		cdata, cw,ch, 0,0, VIGNET_CPY);
/*---- Weight the data */
      vignet_copy(sample->vigweight, set->vigsize[0], set->vigsize[1],
		cvigw, cw,ch, 0,0, VIGNET_CPY);

      for (cdatat=cdata, cvigwt=cvigw, i=ncpix; i--;)
        *(cdatat++) *= *(cvigwt++);

      for (j=0; j<PSF_NITER; j++)
        {
/*------ Map the PSF model at the current position */
        vignet_resample(psf->loc, psf->size[0], psf->size[1],
		cbasis, cw,ch, -dx*vigstep, -dy*vigstep, vigstep, 1.0);

/*------ Build the a and b matrices */
        memset(amat, 0, 9*sizeof(double));
        *bmat = bmat[1] = bmat[2] = mx2=my2=mxy = 0.0;
        for (cvigxt=cvigx,cvigyt=cvigy,cvigwt=cvigw,
		cbasist=cbasis,cdatat=cdata, i=ncpix; i--;)
          {
          dval = (double)*(cbasist++);
          *bmat += (dwval = dval*(double)*(cdatat++));
          bmat[1] += dwval*(dvalx = *(cvigxt++) - dx);
          bmat[2] += dwval*(dvaly = *(cvigyt++) - dy);
          mx2 += dval*dvalx*dvalx;
          my2 += dval*dvaly*dvaly;
          mxy += dval*dvalx*dvaly;
          amatt=amat;
          *(amatt++) += (dval *= dval*(double)*(cvigwt++));
          *(amatt++) += dval*dvalx;
          *(amatt++) += dval*dvaly;
          *(++amatt) += dval*dvalx*dvalx;
          *(++amatt) += dval*dvalx*dvaly;
          *(amatt+3) += dval*dvaly*dvaly;
          }

/*------ Solve the system */
        cholsolve(amat,bmat, 3);

/*------ Convert to a shift */
        dx += 0.5*(ddx = (bmat[1]*mx2 + bmat[2]*mxy) / bmat[0]); 
        dy += 0.5*(ddy = (bmat[2]*my2 + bmat[1]*mxy) / bmat[0]); 
/*------ Exit if it converges or diverges */
        if (ddx*ddx+ddy*ddy < radmin2)
          {
          okflag = 1;
          break;
	  }
        else if (dx*dx+dy*dy > radmax2)
          break;
        }
      if (okflag)
        {
        sample->dx = dx;
        sample->dy = dy;
        }
      }


/*-- Map the PSF model at the current position */
    vignet_resample(psf->loc, psf->size[0], psf->size[1],
	sample->vigresi, set->vigsize[0], set->vigsize[1],
	-dx*vigstep, -dy*vigstep, vigstep, 1.0);
/*-- Subtract the PSF model and compute Chi2 */
    norm = sample->norm;
    chi2 = 0.0;
    dresit = dresi;
    for (vigresi=sample->vigresi, vig=sample->vig, vigw=sample->vigweight,
		i=npix; i--; vigresi++)
      {
      *vigresi = fval = *(vig++)-*vigresi*norm;
      chi2 += (double)(fval *= fval**(vigw++));
      *(dresit++) += fval;
      }

    sample->chi2 = (nm1 = (double)(npix - 1)) > 0.0? chi2/nm1 : chi2;
    }

/* Normalize and convert to floats the Residual array */
  QMALLOC(fresi, float, npix); 
  nm1 = nsample > 1?  (double)(nsample - 1): 1.0;
  for (dresit=dresi,fresit=fresi, i=npix; i--;)
      *(fresit++) = sqrt(*(dresit++)/nm1);

/*-- Map the residuals to PSF coordinates */
  vignet_resample(fresi, set->vigsize[0], set->vigsize[1],
	psf->resi, psf->size[0], psf->size[1], 0.0,0.0, psf->pixstep, 1.0);

/* Free memory */
  free(dresi);
  free(fresi);
  if (centflag)
    {
    free(cvigx);
    free(cvigy);
    free(cvigw);
    free(cbasis);
    free(cdata);
    }

  return;
  }


/******************************* psf_refine **********************************/
/*
Refine PSF by resolving "aliased" pixels.
*/
void	psf_refine(psfstruct *psf, setstruct *set, int npsf)
  {
   polystruct		*poly;
   samplestruct		*sample;
   static double	pos[MAXCONTEXT];
   static char		str[MAXCHAR];
   double		*pix, *desmat,*desmatt,*desmatt2, *desmat0,*desmat02,
			*bmat,*bmatt, *basis,*basist, *basist2,
			*sigvig,*sigvigt, *alphamat,*alphamatt,
			*betamat,*betamatt, *coeffmat,*coeffmatt,
			dx,dy, dval, norm;
   float		*vig,*vigt,*vigt2, *wvig, *diracpsf,*diracpsft,
			*diracvig,*diracvigt, *ppix,*ppixt,
			*psforder,*psfordert,
			psfthresh, vigstep, val;
   int			*psfmask,*psfmaskt, *desindex,*desindext,*desindext2,
			*desindex0,*desindex02;
   int			i,j,jo,k,l,c,n, npix,nvpix, ndata,ncoeff,nsample,
			ncontext, nunknown, matoffset, dindex;

/* Exit if no pixel is to be "refined"  */
  if (!npsf)
    return;

/* First Select the brightest pixels */
  NFPRINTF(OUTPUT,"Selecting pixels...");
  npix = psf->size[0]*psf->size[1];
  if (npsf>npix)
    npsf=npix;
  QMEMCPY(psf->comp, psforder, float, npix);
  for (psfordert=psforder, i=npix; i--; psfordert++)
    *psfordert = fabs(*psfordert);
  hmedian(psforder, npix);
  psfthresh = psforder[npix-npsf];
  free(psforder);

/* Mark pixels which have to be reexamined */
  QCALLOC(psfmask, int, npix);
  npsf = 0;
  for (psfmaskt=psfmask, ppix=psf->comp, i=npix; i--; psfmaskt++)
    if (fabs(*(ppix++))>=psfthresh)
      {
      npsf++;
      *psfmaskt = 1;
      }

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "%d PSF pixels retained for super-resolution\n", npsf);

  nvpix = set->vigsize[0]*set->vigsize[1];
  poly = psf->poly;
  ncontext = set->ncontext;
  ncoeff = poly->ncoeff;
  nsample = set->nsample;
  nunknown = ncoeff*npsf;
  vigstep = 1/psf->pixstep;
/* Prepare a PSF mask that will contain Dirac peaks only... */
  QCALLOC(diracpsf, float, npix);
/* ... and a vignet that will contain interpolating coefficients */
  QCALLOC(diracvig, float, nvpix);

  NFPRINTF(OUTPUT,"Processing samples...");
/* Size of the compressed design matrix along the "data" axis */
  ndata = (1+(int)(INTERPW*psf->pixstep))*(1+(int)(INTERPH*psf->pixstep)) + 1;
  matoffset =nunknown-ncoeff;		/* Offset between matrix coeffs */
/* Set-up the (compressed) design matrix and data vector */
  QCALLOC(desmat, double, npsf*ndata);
  QCALLOC(desindex, int, npsf*ndata);
  QMALLOC(bmat, double, nvpix);
/* ... a matrix containing the context coefficient submatrix... */
  QMALLOC(coeffmat, double, ncoeff*ncoeff);
/* ... a vignet that will contain the current vignet residuals... */
  QMALLOC(vig, float, nvpix);
/* ... a vignet that will contain the current 1/sigma map... */
  QMALLOC(sigvig, double, nvpix);
/* ... and allocate some more for storing the normal equations */
  QCALLOC(alphamat, double, nunknown*nunknown);
  QCALLOC(betamat, double, nunknown);

/* Go through each sample */
  for (sample=set->sample, n=0; n<nsample ; n++, sample++)
    {
    sprintf(str, "Processing sample #%d", n+1);
    NFPRINTF(OUTPUT, str);
/*-- Delta-x and Delta-y in PSF-pixel units */
    dx = -sample->dx*vigstep;
    dy = -sample->dy*vigstep;

/*-- Build the local PSF */
    for (i=0; i<ncontext; i++)
      pos[i] = (sample->context[i]-set->contextoffset[i])
		/set->contextscale[i];
    psf_build(psf, pos);

/*-- Build the current context coefficient sub-matrix */
    basis = poly->basis;
    for (basist=basis, coeffmatt=coeffmat, l=ncoeff; l--;)
      for (dval=*(basist++), basist2=basis, i=ncoeff; i--;)
        *(coeffmatt++) = dval**(basist2++);

/*-- Map the PSF model at the current position */
    vignet_resample(psf->loc, psf->size[0], psf->size[1],
		vig, set->vigsize[0], set->vigsize[1], dx, dy, vigstep, 1.0);
/*-- Subtract the PSF model */
    norm = (double)sample->norm;
    for (vigt=vig, vigt2=sample->vig, i=nvpix; i--; vigt++)
      *vigt = (float)(*(vigt2++) - *vigt*norm);
/*-- Precompute the 1/sigma-map for the current sample */
    for (sigvigt=sigvig, wvig=sample->vigweight, i=nvpix; i--;)
      *(sigvigt++) = sqrt(*(wvig++)) ;

/*-- Go through each relevant PSF pixel */
    desmatt = desmat;
    desindext = desindex;
    for (psfmaskt=psfmask, diracpsft=diracpsf, i=npix; i--; diracpsft++)
      if (*(psfmaskt++))
        {
/*------ Put a Dirac peak at the current PSF pixel */
        *diracpsft = (float)norm;
        vignet_resample(diracpsf, psf->size[0], psf->size[1],
		diracvig, set->vigsize[0],set->vigsize[1], dx,dy, vigstep,1.0);

/*------ Retrieve coefficient for each relevant data pixel */
        for (diracvigt=diracvig, sigvigt=sigvig,
		desmatt2=desmatt, desindext2=desindext, j=jo=0; j++<nvpix;)
          if (fabs(dval = *(diracvigt++) * *(sigvigt++)) > (1/BIG))
            {
            *(desmatt2++) = dval;
            *(desindext2++) = (j-jo);
            jo = j;
            }

        *desindext2 = 0;

        desindext += ndata;
        desmatt += ndata;

/*------ Set the current PSF pixel back to zero */
        *diracpsft = 0.0;
        }

/*-- Fill the b matrix with data points */
    for (vigt=vig, sigvigt=sigvig, bmatt=bmat, j=nvpix; j--;)
      *(bmatt++) = *(vigt++) * *(sigvigt++);

/*-- Compute the matrix of normal equations */
    betamatt = betamat;
    for (desmat0=desmat, desindex0=desindex, k=0; k<npsf;
		desmat0+=ndata, desindex0+=ndata, k++)
      {
      for (desmat02=desmat0, desindex02=desindex0, j=k; j<npsf;
		desmat02+=ndata, desindex02+=ndata, j++)
        {
        dval = 0.0;
        desmatt=desmat0;
        desmatt2=desmat02;
        desindext=desindex0;
        desindext2=desindex02;
        dindex=*desindext-*desindext2;
        while (*desindext && *desindext2)
          {
          while (*desindext && dindex<0)
            {
            dindex+=*(++desindext);
            desmatt++;
            }
          while (*desindext2 && dindex>0)
            {
            dindex-=*(++desindext2);
            desmatt2++;
            }
          while (*desindext && !dindex)
            {
            dval += *(desmatt++)**(desmatt2++);
            dindex = *(++desindext)-*(++desindext2);
            }
          }
        if (fabs(dval) > (1/BIG))
          {
          alphamatt = alphamat+(j+k*npsf*ncoeff)*ncoeff;
          for (coeffmatt=coeffmat, l=ncoeff; l--; alphamatt+=matoffset)
            for (i=ncoeff; i--;)
              *(alphamatt++) += dval**(coeffmatt++);
          }
        }
      dval = 0.0;
      desmatt=desmat0;
      desindext=desindex0;
      bmatt=bmat-1;
      while (*desindext)
        dval += *(desmatt++)**(bmatt+=*(desindext++));
      for (basist=basis,i=ncoeff; i--;)
        *(betamatt++) += dval**(basist++);
      }
    }

/* Free some memory... */
  free(coeffmat);
  free(desmat);
  free(desindex);
  free(bmat);
  free(diracpsf);
  free(diracvig);
  free(vig);
  free(sigvig);

  NFPRINTF(OUTPUT,"Solving the system...");
  cholsolve(alphamat,betamat,nunknown);

  NFPRINTF(OUTPUT,"Updating the PSF...");
  for (ppix=psf->comp,betamatt=betamat, psfmaskt=psfmask, i=npix; i--; ppix++)
    if (*(psfmaskt++))
      for (ppixt=ppix, c=ncoeff; c--; ppixt+=npix)
        *ppixt += *(betamatt++);

/* Free all */
  free(alphamat);
  free(betamat);
  free(psfmask);

  return;
  }


/********************************* psf_make **********************************/
/*
Make the PSF.
*/
void	psf_make(psfstruct *psf, setstruct *set)
  {
   polystruct	*poly;
   samplestruct	*sample;
   double	*pstack,*wstack, *basis, *pix,*wpix, *coeff, *pos, *post;
   float	*comp;
   int		i,c,n, ncoeff,npix,nsample;

  poly = psf->poly;

/* First copy the offset and scaling information from the set structure */
  for (i=0; i<psf->poly->ndim; i++)
    {
    psf->contextoffset[i] = set->contextoffset[i];
    psf->contextscale[i] = set->contextscale[i];
    }

  ncoeff = poly->ncoeff;
  npix = psf->size[0]*psf->size[1];
  nsample = set->nsample;
  QMALLOC(pstack, double, nsample);
  QMALLOC(wstack, double, nsample);
  QMALLOC(basis, double, poly->ncoeff*nsample);
  QMALLOC(pos, double, poly->ndim?(nsample*poly->ndim):1);

  for (sample=set->sample, post=pos, n=nsample; n--; sample++)
    {
    update_retina(set, sample, psf->pixstep);
    for (i=0; i<poly->ndim; i++)
      *(post++) = (sample->context[i]-set->contextoffset[i])
		/set->contextscale[i];
    }

/* Make a polynomial fit to each pixel */
  for (i=0; i<npix; i++)
    {
/*-- Stack ith pixel from each PSF candidate */
    for (sample=set->sample, pix=pstack,wpix=wstack, n=nsample; n--; sample++)
      {
      *(pix++) = (double)*(sample->retina+i);
      *(wpix++) = (double)*(sample->retiweight+i);
      }

/*-- Polynomial fitting */
   poly_fit(poly, i?NULL:pos, pstack, wstack, nsample, basis);

/*-- Store as a PSF component */
    for (coeff=poly->coeff, comp=psf->comp+i,  c=ncoeff; c--; comp+=npix)
      *comp = *(coeff++);
    }

  free(pstack);
  free(wstack);
  free(basis);
  free(pos);

  return;
  }


/********************************* psf_save **********************************/
/*
Save the PSF data as a FITS file.
*/
void	psf_save(psfstruct *psf, pcstruct *pcc, pcstruct *pc, char *filename)
  {
   catstruct	*cat;
   tabstruct	*tab;
   keystruct	*key;
   codestruct	*code;
   float	*w;
   char		*head, str[80], str2[80];
   int		i, npc, temp;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  cat = new_cat(1);
  init_cat(cat);
  tab = new_tab("PSF_DATA");
  add_tab(tab, cat, 0);
/* Let's allocate more than strictly necessary to be sure... */
  QREALLOC(tab->headbuf, char, (3+tab->headnblock)*FBSIZE);
  head = tab->headbuf;
/* ... and blank the extra space */
  memset(head+FBSIZE*tab->headnblock, ' ', 3*FBSIZE);

  fitsadd(head, "POLNAXIS", "Number of context parameters");
  fitswrite(head, "POLNAXIS", &psf->poly->ndim, H_INT, T_LONG);
  for (i=0; i<psf->poly->ndim; i++)
    {
    sprintf(str, "POLGRP%1d", i+1);
    fitsadd(head, str, "Polynom group for this context parameter");
    temp = psf->poly->group[i]+1;
    fitswrite(head, str, &temp, H_INT, T_LONG);
    sprintf(str, "POLNAME%1d", i+1);
    fitsadd(head, str, "Name of this context parameter");
    fitswrite(head, str, psf->contextname[i], H_STRING, T_STRING);
    sprintf(str, "POLZERO%1d", i+1);
    fitsadd(head, str, "Offset value for this context parameter");
    fitswrite(head, str, &psf->contextoffset[i], H_EXPO, T_DOUBLE);
    sprintf(str, "POLSCAL%1d", i+1);
    fitsadd(head, str, "Scale value for this context parameter");
    fitswrite(head, str, &psf->contextscale[i], H_EXPO, T_DOUBLE);
    }

  fitsadd(head, "POLNGRP", "Number of context groups");
  fitswrite(head, "POLNGRP", &psf->poly->ngroup, H_INT, T_LONG);
  for (i=0; i<psf->poly->ngroup; i++)
    {
    sprintf(str, "POLDEG%1d", i+1);
    fitsadd(head, str, "Polynom degree for this context group");
    fitswrite(head, str, &psf->poly->degree[i], H_INT, T_LONG);
    }

/* Add and write important scalars as FITS keywords */
  fitsadd(head, "PSF_SAMP", "Sampling step of the PSF data");
  fitswrite(head, "PSF_SAMP", &psf->pixstep, H_FLOAT, T_FLOAT);
  fitsadd(head, "PSFNAXIS", "Dimensionality of the PSF data");
  fitswrite(head, "PSFNAXIS", &psf->dim, H_INT, T_LONG);
  for (i=0; i<psf->dim; i++)
    {
    sprintf(str, "PSFAXIS%1d", i+1);
    fitsadd(head, str, "Number of element along this axis");
    fitswrite(head, str, &psf->size[i], H_INT, T_LONG);
    }

/* Create and fill the arrays */
  key = new_key("PSF_MASK");
  key->naxis = psf->dim;
  QMALLOC(key->naxisn, int, key->naxis);
  for (i=0; i<psf->dim; i++)
    key->naxisn[i] = psf->size[i];
  strcat(key->comment, "Tabulated PSF data");
  key->htype = H_FLOAT;
  key->ttype = T_FLOAT;
  key->nbytes = psf->npix*t_size[T_FLOAT];
  key->nobj = 1;
  key->ptr = psf->comp;
  add_key(key, tab, 0);

/* Find the useful FITS header size */
  tab->headnblock = fitsfind(head, "END     ")/(FBSIZE/80)+1;

/* Save the convolved PC components if available */
  if (pc)
    {
    tab = new_tab("PC_DATA");
    add_tab(tab, cat, 0);
/*-- Let's allocate more than strictly necessary to be sure... */
    QREALLOC(tab->headbuf, char, (3+tab->headnblock)*FBSIZE);
    head = tab->headbuf;
/*-- ... and blank the extra space */
    memset(head+FBSIZE*tab->headnblock, ' ', 3*FBSIZE);

/*-- Add and write important scalars as FITS keywords */
    fitsadd(head, "PCNAXIS", "Dimensionality of the PC data");
    fitswrite(head, "PCNAXIS", &pcc->dim, H_INT, T_LONG);
    for (i=0; i<pcc->dim; i++)
      {
      sprintf(str, "PCAXIS%1d", i+1);
      fitsadd(head, str, "Number of element along this axis");
      fitswrite(head, str, &pcc->size[i], H_INT, T_LONG);
      }

/*-- Create and fill the arrays */
    key = new_key("PC_CONVMASK");
    key->naxis = pcc->dim;
    QMALLOC(key->naxisn, int, key->naxis);
    for (i=0; i<pcc->dim; i++)
      key->naxisn[i] = pcc->size[i];
    strcat(key->comment, "Convolved data vector, pixel-by-pixel");
    key->htype = H_FLOAT;
    key->ttype = T_FLOAT;
    key->nbytes = pcc->npix*t_size[T_FLOAT];
    key->nobj = 1;
    key->ptr = pcc->comp;
    add_key(key, tab, 0);

    npc = pc->size[pc->dim-1];
    key = new_key("PC_MASK");
    key->naxis = pc->dim;
    QMALLOC(key->naxisn, int, key->naxis);
    for (i=0; i<pc->dim; i++)
      key->naxisn[i] = pc->size[i];
    strcat(key->comment, "Original data vector, pixel-by-pixel");
    key->htype = H_FLOAT;
    key->ttype = T_FLOAT;
    key->nbytes = pc->npix*t_size[T_FLOAT];
    key->nobj = 1;
    key->ptr = pc->comp;
    add_key(key, tab, 0);

    key = new_key("PC_MX2");
    key->naxis = 1;
    QMALLOC(key->naxisn, int, key->naxis);
    *key->naxisn = npc;
    strcat(key->comment, "Variance along x for each component");
    key->htype = H_EXPO;
    key->ttype = T_DOUBLE;
    key->nbytes = npc*t_size[T_DOUBLE];
    key->nobj = 1;
    key->ptr = pc->mx2;
    add_key(key, tab, 0);

    key = new_key("PC_MY2");
    key->naxis = 1;
    QMALLOC(key->naxisn, int, key->naxis);
    *key->naxisn = npc;
    strcat(key->comment, "Variance along y for each component");
    key->htype = H_EXPO;
    key->ttype = T_DOUBLE;
    key->nbytes = npc*t_size[T_DOUBLE];
    key->nobj = 1;
    key->ptr = pc->my2;
    add_key(key, tab, 0);

    key = new_key("PC_MXY");
    key->naxis = 1;
    QMALLOC(key->naxisn, int, key->naxis);
    *key->naxisn = npc;
    strcat(key->comment, "Covariance for each component");
    key->htype = H_EXPO;
    key->ttype = T_DOUBLE;
    key->nbytes = npc*t_size[T_DOUBLE];
    key->nobj = 1;
    key->ptr = pc->mxy;
    add_key(key, tab, 0);

    key = new_key("PC_FLUX");
    key->naxis = 1;
    QMALLOC(key->naxisn, int, key->naxis);
    *key->naxisn = npc;
    strcat(key->comment, "Flux in each component");
    key->htype = H_EXPO;
    key->ttype = T_DOUBLE;
    key->nbytes = npc*t_size[T_DOUBLE];
    key->nobj = 1;
    key->ptr = pc->flux;
    add_key(key, tab, 0);

    key = new_key("PC_BRATIO");
    key->naxis = 1;
    QMALLOC(key->naxisn, int, key->naxis);
    *key->naxisn = npc;
    strcat(key->comment, "B/T for each component");
    key->htype = H_EXPO;
    key->ttype = T_DOUBLE;
    key->nbytes = npc*t_size[T_DOUBLE];
    key->nobj = 1;
    key->ptr = pc->bt;
    add_key(key, tab, 0);

  if (pcc->code)
    {
    code = pcc->code;
    fitsadd(head, "NCODE", "Total number of Codebook vectors");
    fitswrite(head, "NCODE", &code->ncode, H_INT, T_LONG);
    fitsadd(head, "NCODEPAR", "Total number of Codebook vectors");
    fitswrite(head, "NCODEPAR", &code->nparam, H_INT, T_LONG);
    key = new_key("CODE_PC");
    key->naxis = 2;
    QMALLOC(key->naxisn, int, key->naxis);
    key->naxisn[0] = npc;
    key->naxisn[1] = code->ncode;
    strcat(key->comment, "Codebook vectors");
    key->htype = H_EXPO;
    key->ttype = T_FLOAT;
    key->nbytes = key->naxisn[0]*key->naxisn[1]*t_size[key->ttype];
    key->nobj = 1;
    key->ptr = code->pc;
    add_key(key, tab, 0);
    for (i=0; i<code->nparam; i++)
      {
      sprintf(str, "CODE_P%d", i+1);
      key = new_key(str);
      key->naxis = 2;
      QMALLOC(key->naxisn, int, key->naxis);
      key->naxisn[0] = 1;
      key->naxisn[1] = code->ncode;
      sprintf(str, "Codebook parameter #%d", i+1);
      strcat(key->comment, str);
      key->htype = H_EXPO;
      key->ttype = T_FLOAT;
      key->nbytes = key->naxisn[0]*key->naxisn[1]*t_size[key->ttype];
      key->nobj = 1;
      key->ptr = code->param[i];
      add_key(key, tab, 0);
      sprintf(str, "CODE_M%d", i+1);
      sprintf(str2, "Step modulus for Codebook parameter #%d", i+1);
      fitsadd(head, str, str2);
      fitswrite(head, str, &code->parammod[i], H_INT, T_LONG);
      }
    }

/*-- Find the useful FITS header size */
    tab->headnblock = fitsfind(head, "END     ")/(FBSIZE/80)+1;
    }

/* Then, just save everything and free memory */
  save_cat(cat, filename);

/* But don't touch my arrays!! */
  blank_keys(tab->prevtab);
  blank_keys(tab);
  free_cat(cat, 1);

  return;
  }


