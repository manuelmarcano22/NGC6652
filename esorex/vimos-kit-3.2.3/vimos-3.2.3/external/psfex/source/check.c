  /*
 				check.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	Production of check-images for the PSF.
*
*	Last modify:	14/03/99
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
#include	"check.h"


/******************************* psf_writecheck ******************************/
/*
Write a FITS image for check.
*/
void	psf_writecheck(psfstruct *psf, pcstruct *pc, setstruct *set,
		char *filename, checkenum checktype)
  {
   catstruct		*cat;
   tabstruct		*tab;
   samplestruct		*sample;
   char			*head;
   static double	dpos[POLY_MAXDIM], *dpost;
   double		*dpix0, *dpix, dstep,dstart;
   float		*pix,*pix0, *fpix, val;
   int			i,j,x,y, w,h,n, npc,nt, nw,nh, step;

/* Create the new cat (well it is not a "cat", but simply a FITS table */
  cat = new_cat(1);
  init_cat(cat);
  tab = cat->tab;
  tab->naxis = 2;	/* This is an image */
  QMALLOC(tab->naxisn, int, tab->naxis);
  fitsremove(tab->headbuf,"HISTORY ");
  fitsremove(tab->headbuf,"EXTEND  ");
  head = tab->headbuf;
  switch(checktype)
    {
    case PSF_CHI:
/*---- sqrt(chi2) map in PSF pixel-space */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      nw = (int)(sqrt((double)set->nsample*set->nvig)/set->vigsize[0]+1);
      nw = 1;
      nh = 1;
      w = psf->size[0];
      h = psf->size[1];
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      pix = pix0;
      fpix = psf->resi;
      for (i=w*h; i--;)
        *(pix++) = *(fpix++);
      break;

    case PSF_PROTO:
/*---- PSF data for all components are arranged as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      npc = psf->size[2];
      nw = npc<10? npc:10;
      nh = (npc-1)/nw + 1;
      w = psf->size[0];
      h = psf->size[1];
      step = (nw-1)*w;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
/*---- Normalize the components in the image corners: pos=(0.5,0.5,..) */
      for (dpost=dpos, i=psf->poly->ndim; i--;)
        *(dpost++) = 0.5;
      poly_func(psf->poly, dpos);
      dpost = psf->poly->basis;
      fpix = psf->comp;
      for (n=0; n<npc; n++)
        {
        val = *(dpost++);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++)*val;
        }
      break;

    case PSF_RESIDUALS:
/*---- Residual vectors for all samples are arranged as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      nw = (int)(sqrt((double)set->nsample*set->nvig)/set->vigsize[0]+1);
      nw = ((nw-1)/10+1)*10;
      nh = set->nsample/(nw+1) + 1;
      w = set->vigsize[0];
      h = set->vigdim>1? set->vigsize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      sample = set->sample;
      for (n=0; n<set->nsample; n++)
        {
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = (sample++)->vigresi;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++);
        }
      break;

    case PSF_RAWDATA:
/*----  View original samples as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      nw = (int)(sqrt((double)set->nsample*set->nvig)/set->vigsize[0]+1);
      nw = ((nw-1)/10+1)*10;
      nh = set->nsample/(nw+1) + 1;
      w = set->vigsize[0];
      h = set->vigdim>1? set->vigsize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      sample = set->sample;
      for (n=0; n<set->nsample; n++)
        {
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = (sample++)->vig;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++);
        }
      break;

    case PSF_SAMPLES:
/*----  View all training samples as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      nw = (int)(sqrt((double)set->nsample*set->nreti)/set->retisize[0]+1);
      nw = ((nw-1)/10+1)*10;
      nh = set->nsample/(nw+1) + 1;
      w = set->retisize[0];
      h = set->retidim>1? set->retisize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      sample = set->sample;
      for (n=0; n<set->nsample; n++)
        {
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = (sample++)->retina;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++);
        }
      break;

    case PSF_SNAPSHOTS:
/*----  View reconstructed PSFs as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      npc = psf->poly->ndim;
      nw = PSF_SNAPWIDTH;
      for (nt=PSF_SNAPWIDTH*PSF_SNAPWIDTH, i=npc-2; (i--)>0;)
        nt *= PSF_SNAPWIDTH;
      nh = nt/nw;
      w = set->retisize[0];
      h = set->retidim>1? set->retisize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      dstep = 1.0/PSF_SNAPWIDTH;
      dstart = (1.0-dstep)/2.0;
      memset(dpos, 0, POLY_MAXDIM*sizeof(double));
      for (i=0; i<npc; i++)
        dpos[i] = -dstart;
      for (n=0; n<nt; n++)
        {
        psf_build(psf, dpos);
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = psf->loc;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++);
        for (i=0; i<npc; i++)
          if (dpos[i]<dstart-0.01)
            {
            dpos[i] += dstep;
            break;
            }
          else
            dpos[i] = -dstart;
        }
      break;

    case PSF_WEIGHTS:
/*----  View all training sample weights as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      nw = (int)(sqrt((double)set->nsample*set->nreti)/set->retisize[0]+1);
      nw = ((nw-1)/10+1)*10;
      nh = set->nsample/(nw+1) + 1;
      w = set->retisize[0];
      h = set->retidim>1? set->retisize[1] : 1;
      tab->naxisn[0] = nw*w;
      tab->naxisn[1] = nh*h;
      step = (nw-1)*w;
      tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
      QCALLOC(pix0, float, tab->tabsize);
      tab->bodybuf = (char *)pix0; 
      sample = set->sample;
      for (n=0; n<set->nsample; n++)
        {
        pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
        fpix = (sample++)->retiweight;
        for (y=h; y--; pix += step)
          for (x=w; x--;)
            *(pix++) = *(fpix++);
        }
      break;

    case PSF_PCPROTO:
/*---- convolved PC data for all components are arranged as small vignets */
      tab->bitpix =  BP_FLOAT;
      tab->bytepix = t_size[T_FLOAT];
      if (pc)
        {
        npc = pc->size[2]*pc->size[3];
        nw = pc->size[2];
        nh = pc->size[3];
        w = pc->size[0];
        h = pc->size[1];
        step = (nw-1)*w;
        tab->naxisn[0] = nw*w;
        tab->naxisn[1] = nh*h;
        tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
        QCALLOC(pix0, float, tab->tabsize);
        tab->bodybuf = (char *)pix0; 
        fpix = pc->comp;
        for (n=0; n<npc; n++)
          {
          pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w;
          for (y=h; y--; pix += step)
            for (x=w; x--;)
              *(pix++) = *(fpix++);
          }
        }
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Yet unavailable CHECKIMAGE type",
	"");
    }

/* Then, just save everything and free memory */
  save_cat(cat, filename);
  free_cat(cat, 1);

  return;
  }

