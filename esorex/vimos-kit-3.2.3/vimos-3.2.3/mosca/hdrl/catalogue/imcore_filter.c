/* $Id: imcore_filter.c,v 1.3 2015/08/12 11:16:55 jim Exp $
 *
 * This file is part of the CASU Pipeline utilities
 * Copyright (C) 2015 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: jim $
 * $Date: 2015/08/12 11:16:55 $
 * $Revision: 1.3 $
 * $Name:  $
 */

#include <stdlib.h>

#include "util.h"
#include "imcore.h"

/* Function prototypes */

static void sortm (float ia[], intptr_t ib[], intptr_t n);
static void quicksort (float x[], intptr_t point[], intptr_t l, intptr_t nfilt);
static void filt1d(float [], intptr_t, intptr_t);
static void padext(float [], intptr_t);
static void hanning(float [], intptr_t);

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do bilinear median and linear filtering on background values
  
    \par Name:
        imcore_bfilt
    \par Purpose:
        Do bilinear median and linear filtering on background values
    \par Description:
        A map is smoothed using sliding median and mean filters.
    \par Language:
        C
    \param xbuf
        The input map to be smoothed
    \param nx
        The X dimension of the map
    \param ny
        The Y dimension of the map
    \returns
        Nothing
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/
 
extern void imcore_bfilt (float **xbuf, intptr_t nx, intptr_t ny) {
   float *ybuf, *save;
   intptr_t mfilt = 5, j, k;
 
   /*CONSTANTS 5,1000 */
/* Allocate temporary storage */
   ybuf = (float *) cpl_malloc(MAX(nx,ny) * sizeof(float));
   save = (float *) cpl_malloc((nx+1) * ny * sizeof(float));
/*    if(!ybuf || !save) bombout(1, "malloc"); */
 
/* median filter across */
   for(k = 0; k < ny; k++) {
     for(j = 0; j < nx; j++) {
       save[(nx+1)*k+j] = xbuf[k][j];
       ybuf[j] = xbuf[k][j];
     }
     filt1d(ybuf, nx, mfilt);
     for(j = 0; j < nx; j++) xbuf[k][j] = ybuf[j];
   }
 
/* and now down */
   for(k = 0; k < nx; k++) {
     for(j = 0; j < ny; j++) ybuf[j] = xbuf[j][k];
     filt1d(ybuf, ny, mfilt);
     for(j = 0; j < ny; j++)
/* make sure median filtered values are not large than original */
       if(save[(nx+1)*j+k] > -1000.0)
       xbuf[j][k] = MIN(save[(nx+1)*j+k], ybuf[j]);
   }

/* now repeat with linear filters across */
   for(k = 0; k < ny; k++) {
     for(j = 0; j < nx; j++) ybuf[j] = xbuf[k][j];
     hanning(ybuf, nx);
     for(j = 0; j < nx; j++) xbuf[k][j] = ybuf[j];
   }

/* and now down */
   for(k = 0; k < nx; k++) {
     for(j = 0; j < ny; j++) ybuf[j] = xbuf[j][k];
     hanning(ybuf, ny);
     for(j = 0; j < ny; j++) xbuf[j][k] = ybuf[j];
   }

/* Free temporary storage */
   cpl_free((void *) ybuf);
   cpl_free((void *) save);
}/* --------------------------------------- ------------------------------ */

/**@}*/
 
/* does median filtering allowing for unmeasured entries */
 
static void filt1d (float ybuf[], intptr_t mpt, intptr_t mfilt) {
   float *wbuf;
   intptr_t i, irc;
/* Allocate temporary storage */
   wbuf = (float *) cpl_malloc(mpt * sizeof(float));
/*    if(!wbuf) bombout(1, "malloc"); */
 /* CONSTANTS: 1000 */
   irc = 0;
   for(i = 0; i < mpt; i++){
     if(ybuf[i] > -1000.0){
       wbuf[irc] = ybuf[i];
       irc++;
     }
   }
   if(irc == 0) {
     cpl_free((void *) wbuf);
     return;
   }
   imcore_median(wbuf, irc, mfilt);
   irc = 0;
   for(i = 0; i < mpt; i++){
     if(ybuf[i] > -1000.0){
       ybuf[i] = wbuf[irc];
       irc++;
     }
   }
   padext(ybuf, mpt);
/* Free temporary storage */
   cpl_free((void *) wbuf);
}


/* pads out array with missing points and linearly extrapolates the ends */
 
static void padext (float x[], intptr_t n) {
   intptr_t i, j, ilow, ihih=0, ic;
   float xlow, xhih, slope, t1 ,t2;
/* elements <= -1000.0 are treated as missing */

   /* CONSTANTS: -1000 */
   i = 0;
   while(i < n && x[i] <= -1000.0) i++;
   ilow = i;
   for(i = ilow+1; i < n; i++){
     if(x[i] <= -1000.0) {
       ic = 1;
       if (i < n - 1) {
           while(x[i+ic] <= -1000.0) {
               ic++;
               if (i+ic >= n-1)
                   break;
           }
       }
       if(i+ic < n-1){
         xlow = x[i-1];
         xhih = x[i+ic];
         for(j = 0; j < ic; j++){
           t2 = ((float) j+1)/((float) ic+1);
           t1 = 1.0 - t2;
           x[i+j] = t1*xlow+t2*xhih;
         }
       }
     } else {
       ihih = i;
     }
   }
/* linear extrapolation of ends */
   if(ilow > 0){
     if (ilow < n - 1)   
       slope = x[ilow+1]-x[ilow];
     else 
       slope = 0.0;
     for(i = 0; i < ilow; i++) x[i] = x[ilow]-slope*(ilow-i);
   }
   if(ihih < n-1) {
     if (ihih > 0)
       slope = x[ihih]-x[ihih-1];
     else
       slope = 0.0;
     for(i = ihih+1; i < n; i++) x[i] = x[ihih]+slope*(i-ihih);
   }
}

/* performs linear filtering on array xbuf */

static void hanning (float xbuf[], intptr_t npt) {
  float *ybuf;
  float sum = 0.0, xmns, xmnf;
  intptr_t nfilt = 3, i, il, ilow, nelem;

  /* CONSTANTS: 2,3,4,0.25 */
  if(npt <= nfilt)
    return;

  /* set first and last edges equal */
  il   = nfilt/2;
  ilow = MAX(3,nfilt/4);
  ilow = (ilow/2)*2 + 1;

  for(i = 0; i < ilow; i++)
    sum += xbuf[i];

  xmns = sum/((float) ilow);

  sum=0.0;
  for(i = 0; i < ilow; i++)
    sum += xbuf[npt-1-i];

  xmnf = sum/((float) ilow);

  /* allocate ybuf array */
  nelem = npt + nfilt;  /* Max. number of elements req'd */

  ybuf = (float *) cpl_malloc(nelem * sizeof(float));
/*   if(!ybuf) */
/*     bombout(1, "malloc"); */

  /* reflect edges before filtering */
  for(i = 0; i < il; i++) {
    ybuf[i] = 2.0 * xmns - xbuf[il+ilow-1-i];
    ybuf[npt+i+il] = 2.0 * xmnf - xbuf[npt-i-ilow-1];
  }

  for(i = 0; i < npt; i++)
    ybuf[i+il] = xbuf[i];

  /* do linear filtering on rest */
  for(i = 0; i < npt; i++)
    xbuf[i] = 0.25 * (ybuf[i] + 2.0 * ybuf[i+1] + ybuf[i+2]);  /* 1-2-1 Hanning weighting */

  cpl_free((void *) ybuf);
}

/* performs median filtering on array xbuf */
/*---------------------------------------------------------------------------*/
/** TODO MISSING DOXYGEN
 *
 */


/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief compute median

    \par Name:
        imcore_median
    \par Purpose:
        Compute median
    \par Description:
        Compute median
    \par Language:
        C
    \param xbuf
        x buffer
    \param npt
        Number of points
    \param nfilt
        Size of median filter
    \retval void

    \author
        Jim Lewis, CASU
 */

extern void imcore_median (float xbuf[], intptr_t npt, intptr_t nfilt) {
  float *ybuf, *array;
  float xmns, xmnf;
  intptr_t *point;
  intptr_t nfo2p1, i, il, ilow, j, jl, jh, nelem, l=0;

  if((nfilt/2)*2 == nfilt) nfilt++;
  if(npt <= nfilt) return;
  nfo2p1 = nfilt/2;
/* CONSTANTS: 2,3,4 */
  /* allocate ybuf, array, point */
  nelem = npt + nfilt;  /* Max. number of elements req'd */
  ybuf = (float *) cpl_malloc(nelem * sizeof(float));
/*   if(!ybuf) */
/*     bombout(1, "malloc"); */
  array = (float *) cpl_malloc(nfilt * sizeof(float));
  point = (intptr_t *) cpl_malloc(nfilt * sizeof(intptr_t));
/*   if(!array || !point) */
/*     bombout(1, "malloc"); */

  /* set first and last edges equal */
  il   = nfilt/2;
  ilow = MAX(3, nfilt/4);
  ilow = (ilow/2)*2 + 1;

  for(i = 0; i < ilow; i++) array[i] = xbuf[i];
  sortm(array, point, ilow);
  xmns = array[ilow/2];

  for(i = 0; i < ilow; i++) array[i] = xbuf[npt-1-i];
  sortm(array, point, ilow);
  xmnf = array[ilow/2];

  /* reflect edges before filtering */
  for(i = 0; i < il; i++) {
    ybuf[i] = 2.0 * xmns - xbuf[il+ilow-1-i];
    ybuf[npt+i+il] = 2.0 * xmnf - xbuf[npt-i-ilow-1];
  }
  for(i = 0; i < npt; i++) ybuf[i+il] = xbuf[i];

  /* do median filtering on rest */
  for(i = 0; i < nfilt; i++) {
    array[i] = ybuf[i];
    point[i] = i+1;
  }

  sortm(array, point, nfilt);

  xbuf[0] = array[nfo2p1];
  jl = nfilt;
  jh = nfilt+npt-1;
  for(j = jl; j < jh; j++) {

    for(i = 0; i < nfilt; i++) {
      if(point[i] != 1) {
        point[i]--;
        continue;
      }
      point[i] = nfilt;
      array[i] = ybuf[j];
      l = i;
    }
    quicksort(array, point, l, nfilt);
    xbuf[j-jl+1] = array[nfo2p1];
  }

  /* Free temporary arrays */
  cpl_free((void *) point);
  cpl_free((void *) array);
  cpl_free((void *) ybuf);
}

static void sortm (float ia[], intptr_t ib[], intptr_t n) {
  intptr_t i,j, ii, jj, ifin, iu;
  float it;
/* CONSTATS: 2,3,4 */
  jj = 2;
  while(jj < n) jj = 2 * jj;
  jj = MIN(n,(3 * jj)/4 - 1);
  while(jj > 1) {
    jj = jj/2;
    ifin = n - jj;
    for(ii = 0; ii < ifin; ii++) {
      i = ii;
      j = i + jj;
      if(ia[i] <= ia[j]) continue;
      it = ia[j];
      iu = ib[j];
      do {
        ia[j] = ia[i];
        ib[j] = ib[i];
        j = i;
        i = i - jj;
        if (i < 0) break;
      } while(ia[i] > it);
      ia[j] = it;
      ib[j] = iu;
    }
  }
}

static void quicksort (float x[], intptr_t point[], intptr_t l, intptr_t nfilt) {
  float test, temp;
  intptr_t i, it, j, npt, ii;

  test = x[l];
  j = nfilt;
  for(i = 0; i < nfilt; i++) {
    if(i != l && test <= x[i]) {
      j = i;
      break;
    }
  }
  if(j - 1 == l) return;
  
  if(j < l) {
    temp = x[l];
    it = point[l];
    npt = l - j;
    for(i = 0; i < npt; i++) {
      ii = l - i - 1;
      x[ii+1] = x[ii];
      point[ii+1] = point[ii];
    }
    x[j] = temp;
    point[j] = it;
  }
  else if(j > l) {
    temp = x[l];
    it = point[l];
    j--;
    npt = j - l;
    if(npt != 0) {
      for(i = 0; i < npt; i++) {
        ii = l + i + 1;
        x[ii-1] = x[ii];
        point[ii-1] = point[ii];
      }
    }
    x[j] = temp;
    point[j] = it;
  }
}


/*

$Log: imcore_filter.c,v $
Revision 1.3  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.2  2015/08/07 13:06:54  jim
Fixed copyright to ESO

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.3  2014/04/09 09:09:51  jim
Detabbed

Revision 1.2  2013/08/28 05:43:22  jim
Fixed padext for the situation where all the row is bad

Revision 1.1.1.1  2013-08-27 12:07:48  jim
Imported


*/
