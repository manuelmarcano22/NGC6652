/*
 * This file is part of the CASU Pipeline utilities
 * Copyright (C) 2015,2016 European Southern Observatory
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cpl.h>

#include "floatmath.h"
#include "util.h"
#include "imcore.h"
#include "casu_utilfunctions.h"
#include "hdrl_utils.h"

static void sortit (float [], intptr_t);

cpl_image *
hdrl_sigclipfilter_image_grid(const cpl_image * ima, cpl_matrix * x, cpl_matrix * y,
                             cpl_size  filtersize_x, cpl_size filtersize_y);

/**@{*/


/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Model and create background map

    \par Name:
        imcore_background
    \par Purpose:
        Model and create background map
    \par Description:
        The image data array is split into cells. In each cell a robust
        background estimate is obtained. The cell raster is gently smoothed
        and then used to create a full background map with a bi-linear
        interpolation scheme.
    \par Language:
        C
    \param ap
        The current ap structure
    \param nbsize
        The size of the cells in pixels
    \param nullval
        A null value used to flag bad pixels
    \param bkg_subtr
        switch: if bkg_subtr !=0 fill background map
    \param res
        TODO
    \retval CASU_OK
        If all went well. This is currently the only value.
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_background(ap_t *ap, intptr_t nbsize, float nullval, int bkg_subtr,
                             hdrl_imcore_result * res) {
    float fracx,fracy,avsky,fnbsize,dely,delx;
    float t1,t2,dsky,*map,**bvals,*work;
    intptr_t ifracx,ifracy,nbsizx,nbsizy,nbx,nby;
    intptr_t nbsizo2,kk,k,iby,ibyp1,ibx,ibxp1;
    unsigned char *mflag;
    intptr_t nx,ny;

    /* Set up some variables */

    map = ap->indata;
    mflag = ap->mflag;
    nx = ap->lsiz;
    ny = ap->csiz;

    /* check to see if nbsize is close to exact divisor */

    nbsize = MIN(MIN(nx, ny), nbsize);
    fracx = ((float)nx)/((float)nbsize);
    fracy = ((float)ny)/((float)nbsize);
    ifracx = (intptr_t)(fracx + 0.1);
    ifracy = (intptr_t)(fracy + 0.1);
    nbsizx = nx/ifracx;
    nbsizy = ny/ifracy;
    nbsize = MAX(NINT(0.9*nbsize), MIN(nbsize, MIN(nbsizx,nbsizy)));
    nbsize = MIN(nx,MIN(ny,nbsize)); /* trap for small maps */

    /* Divide the map into partitions */

    nbx = nx/nbsize;
    nby = ny/nbsize;

    /* Same for background values array */

    bvals = cpl_malloc(nby*sizeof(float *));
    for (intptr_t l = 0; l < nby; l++)
        bvals[l] = cpl_malloc(nbx*sizeof(float));

    /* Store some of this away for use later */

    ap->backmap.nbx = nbx;
    ap->backmap.nby = nby;
    ap->backmap.nbsize = nbsize;
    ap->backmap.bvals = bvals;

    /* create cpl image with mask from raw data */
    cpl_image * image = cpl_image_wrap_float((cpl_size)nx, (cpl_size)ny, map);
    cpl_mask * image_mask = cpl_image_get_bpm(image);
    cpl_binary * image_mask_data = cpl_mask_get_data(image_mask);

    /* Add additional bad pixels if present*/
    for (intptr_t i = 0; i < nx * ny ; i++){
        if (map[i] == nullval || mflag[i] == MF_ZEROCONF ||
            mflag[i] == MF_STUPID_VALUE || mflag[i] == MF_SATURATED ) {
            image_mask_data[i] = CPL_BINARY_1;
        }
    }

    const cpl_size steps_x = nbx;
    const cpl_size steps_y = nby;

    int filter_size_x = (int)(nbsize/2);
    int filter_size_y = (int)(nbsize/2);

    cpl_size sx = CX_MAX(nx / steps_x, 1);
    cpl_size sy = CX_MAX(ny / steps_y, 1);

    /* sigclip stepped grid */
    cpl_matrix * x = hdrl_matrix_linspace(sx / 2, nx, sx);
    cpl_matrix * y = hdrl_matrix_linspace(sy / 2, ny, sy);
    cpl_image * imgtmp_mod =
        hdrl_sigclipfilter_image_grid(image, x, y,
                                      filter_size_x, filter_size_y);
    /* interpolate remaining bad pixels */
    cpl_detector_interpolate_rejected(imgtmp_mod);

    cpl_matrix_delete(x);
    cpl_matrix_delete(y);

    for (intptr_t l = 0; l < nby; l++) {
        for (intptr_t j = 0; j < nbx; j++) {
            int rej;
            bvals[l][j] = cpl_image_get(imgtmp_mod, j + 1, l + 1, &rej);
        }
    }

    cpl_image_delete(imgtmp_mod);
    cpl_image_unwrap(image);

    /* filter raw background values */

    imcore_bfilt(bvals,nbx,nby);

    /* compute average sky level */

    work = cpl_malloc(nbx*nby*sizeof(*work));
    k = 0;
    for(intptr_t l = 0; l < nby; l++)
        for(intptr_t j = 0; j < nbx; j++)
            work[k++] = bvals[l][j];
    sortit(work,k);
    avsky = work[(k)/2];
    freespace(work);

    /* ok now correct map for background variations and put avsky back on */

    nbsizo2 = nbsize/2;
    fnbsize = 1.0/((float)nbsize);
    for (k = 0; k < ny; k++) {
        kk = k*nx;

        /* Nearest background pixel vertically */

        iby = (k + 1 + nbsizo2)/nbsize;
        ibyp1 = iby + 1;
        iby = MIN(nby,MAX(1,iby));
        ibyp1 = MIN(nby,ibyp1);
        dely = (k + 1 - nbsize*iby + nbsizo2)*fnbsize;

        for (intptr_t j = 0; j < nx; j++) {
            if (map[kk+j] == nullval)
                continue;

            /* nearest background pixel across */

            ibx = (j + 1 + nbsizo2)/nbsize;
            ibxp1 = ibx + 1;
            ibx = MIN(nbx,MAX(1,ibx));
            ibxp1 = MIN(nbx,ibxp1);
            delx = (j + 1 - nbsize*ibx + nbsizo2)*fnbsize;

            /* bilinear interpolation to find background */

            t1 = (1.0 - dely)*bvals[iby-1][ibx-1] + dely*bvals[ibyp1-1][ibx-1];
            t2 = (1.0 - dely)*bvals[iby-1][ibxp1-1] + dely*bvals[ibyp1-1][ibxp1-1];
            dsky = avsky - (1.0 - delx)*t1 - delx*t2;
            if (bkg_subtr) {
                map[kk+j] += dsky;
                /* Fill the background map */
                if (res->background) {
                    cpl_image_set(res->background,j+1,k+1,avsky - dsky);
                }
            }
        }
    }

    return(CASU_OK);
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out robust background estimate over a whole input image

    \par Name:
        imcore_backstats
    \par Purpose:
        Work out a robust background estimate over a whole input image
    \par Description:
        The image is analysed to work out a robust estimate of the background
        median, and sigma
    \par Language:
        C
    \param ap
        The current ap structure
    \param nullval
        A null value used to flag bad pixels
    \param skymed
        Output sky median
    \param skysig
        Output sky noise
    \retval CASU_OK
        If all went well.
    \retval CASU_WARN
        If there aren't enough good values to do the calculation
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_backstats(ap_t *ap, float nullval,
                            float *skymed, float *skysig)
{
    intptr_t nx, ny;
    float *map;
    unsigned char *mflag;

    /* Get some info from the ap structure */
    map = ap->indata;
    nx = ap->lsiz;
    ny = ap->csiz;
    mflag = ap->mflag;

    cpl_image * ima_wrp = cpl_image_wrap(nx, ny, CPL_TYPE_FLOAT, map);
    cpl_mask * image_mask = cpl_image_get_bpm(ima_wrp);
    cpl_binary * image_mask_data = cpl_mask_get_data(image_mask);

    /* Add additional bad pixels if present*/
    for (intptr_t i = 0; i < nx * ny ; i++){
        if (map[i] == nullval || mflag[i] == MF_ZEROCONF ||
            mflag[i] == MF_STUPID_VALUE || mflag[i] == MF_SATURATED ) {
            image_mask_data[i] = CPL_BINARY_1;
        }
    }

    /* do kappa-sigma clipping to thresh out basic outliers */

    int niter = 30;
    cpl_size rej_new;
    for ( int i = 0; i < niter; i++ ) {
        double mad, kappa_low = 2.5, kappa_high = 2.5;

        double median = cpl_image_get_mad(ima_wrp, &mad);
        double stdev = mad * CPL_MATH_STD_MAD;

        double lo_cut = median - kappa_low * stdev;
        double hi_cut = median + kappa_high * stdev;

        cpl_size rej_orig = cpl_image_count_rejected(ima_wrp);

        if (lo_cut < hi_cut) {
            cpl_mask_threshold_image(image_mask, ima_wrp, lo_cut, hi_cut,
                            CPL_BINARY_0);
        }
        rej_new = cpl_image_count_rejected(ima_wrp);

        if (rej_orig == rej_new) break;
    }

    if (rej_new ==  nx * ny) {
        *skymed = 0.0;
        *skysig = 0.0;
        cpl_image_unwrap(ima_wrp);
        return(CASU_WARN);
    }
    else {

        /* Set the final answer. Here the normal mean and standard deviation is used
         * as all outliers should be masked and thus the mean is more accurate */

        *skymed = cpl_image_get_mean(ima_wrp);
        *skysig = cpl_image_get_stdev(ima_wrp);
    }

    cpl_msg_debug(cpl_func,"Background: Clipped mean: %g Clipped stdev: %g",
                  *skymed, *skysig);

    cpl_image_unwrap(ima_wrp);
    return(CASU_OK);
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Work out estimated sky for a pixel position

    \par Name:
        imcore_backest
    \par Purpose:
        Work out estimated sky for a pixel position
    \par Description:
        Given the coarse background grid, calculate the background at a
        given image pixel position by doing a bi-linear interpolation of
        it's position within the grid.
    \par Language:
        C
    \param ap
        The current ap structure
    \param x
        The X position in question
    \param y
        The Y position in question
    \param skylev
        Output sky level at x,y
    \param skyrms
        Output sky noise at x,y
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

extern void imcore_backest(ap_t *ap, float x, float y, float *skylev,
                           float *skyrms) {
    intptr_t i,j,nbx,nby,nbsize,nbsizo2,iby,ibyp1,ibx,ibxp1;
    float **bvals,fnbsize,dely,delx,t1,t2;

    /* Define some local variables */

    nbx = ap->backmap.nbx;
    nby = ap->backmap.nby;
    nbsize = ap->backmap.nbsize;
    bvals = ap->backmap.bvals;

    /* Get closest pixel to the input location */

    i = NINT(x);
    j = NINT(y);

    /* Now, work out where in the map to do the interpolation */
/* CONSTANTS: 2, 0.25 */
    nbsizo2 = nbsize/2;
    fnbsize = 1.0/((float)nbsize);
    iby = (j + nbsizo2)/nbsize;
    ibyp1 = iby + 1;
    iby = MIN(nby,MAX(1,iby));
    ibyp1 = MIN(nby,ibyp1);
    dely = (j  - nbsize*iby + nbsizo2)*fnbsize;
    ibx = (i + nbsizo2)/nbsize;
    ibxp1 = ibx + 1;
    ibx = MIN(nbx,MAX(1,ibx));
    ibxp1 = MIN(nbx,ibxp1);
    delx = (i - nbsize*ibx + nbsizo2)*fnbsize;

    /* Now do a linear interpolation to find the background. Calculate MAD of
       the four adjacent background cells as an estimate of the RMS */

    t1 = (1.0 - dely)*bvals[iby-1][ibx-1] + dely*bvals[ibyp1-1][ibx-1];
    t2 = (1.0 - dely)*bvals[iby-1][ibxp1-1] + dely*bvals[ibyp1-1][ibxp1-1];
    *skylev = (1.0 - delx)*t1 + delx*t2;
    *skyrms = 0.25*(fabsf(bvals[iby-1][ibx-1] - *skylev) +
                    fabsf(bvals[ibyp1-1][ibx-1] - *skylev) +
                    fabsf(bvals[iby-1][ibxp1-1] - *skylev) +
                    fabsf(bvals[ibyp1-1][ibxp1-1] - *skylev));
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Analyse histogram to work out median and sigma

    \par Name:
        imcore_medsig
    \par Purpose:
        Analyse histogram to work out median and sigma
    \par Description:
        Given a histogram work out the median and sigma of a distribution.
        A starting point in the histogram can be given, which allows you
        to ignore bins at the lower end. A target value defines at what point
        we stop summing the histogram.
    \par Language:
        C
    \param shist
        The input histogram
    \param nh
        The number of bins in the histogram
    \param ist
        The first bin position that we will look at
    \param itarg
        The target value for the summation.
    \param med
        Output median
    \param sig
        Output sigma from the first quartile.
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

extern void imcore_medsig(intptr_t *shist, intptr_t nh, intptr_t ist, intptr_t itarg,
                          float *med, float *sig) {
    intptr_t isum, medata;
    float ffrac,sigmed;

    /* median */
/* CONSTANTS: 3,4,0.5,1.48, CPL_MATH_STD_MAD */
    isum = 0;
    medata = ist;
    while (isum <= (itarg+1)/2 && (medata-MINHISTVAL) < nh) {
        medata++;
        isum += shist[medata-MINHISTVAL];
    }
    if (shist[medata-MINHISTVAL] == 0) {
        ffrac = 0.0;
    } else {
        ffrac = (float)(isum - (itarg+1)/2)/(float)shist[medata-MINHISTVAL];
    }
    *med = (float)medata - ffrac + 0.5;

    /* sigma */

    isum = 0;
    medata = ist;
    while (isum <= (itarg+3)/4 && (medata-MINHISTVAL) < nh) {
        medata++;
        isum += shist[medata-MINHISTVAL];
    }
    if (shist[medata-MINHISTVAL] == 0) {
        ffrac = 0.0;
    } else {
        ffrac = (float)(isum - (itarg+3)/4)/(float)shist[medata-MINHISTVAL];
    }
    sigmed = (float)medata - ffrac + 0.5;
    *sig = CPL_MATH_STD_MAD*(*med - sigmed);
    *sig = MAX(0.5,*sig);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        sortit
    \par Purpose:
        Sort a single array into ascending order
    \par Description:
        A single array is sorted into ascending order
    \par Language:
        C
    \param ia
        The input data array
    \param n
        The number of values in the data array
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void sortit (float ia[], intptr_t n) {
    intptr_t i, j, ii, jj, ifin;
    float it;

    /* CONSTANTS: 2,3,4 */
    jj = 4;
    while (jj < n)
        jj = 2 * jj;
    jj = MIN(n,(3 * jj)/4 - 1);
    while (jj > 1) {
        jj = jj/2;
        ifin = n - jj;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + jj;
            if (ia[i] <= ia[j])
                continue;
            it = ia[j];
            do {
                ia[j] = ia[i];
                j = i;
                i = i - jj;
                if (i < 0)
                    break;
            } while (ia[i] > it);
            ia[j] = it;
        }
    }
    return;
}


/* ---------------------------------------------------------------------------*/
/**
 * @brief filter image on a grid
 *
 * @param ima           image to filter
 * @param x             row vector of coordinates to filter on
 * @param y             row vector of coordinates to filter on
 * @param filtersize_x  size of the median filter
 * @param filtersize_y  size of the median filter
 * @return filtered image of size of the grid
 *
 */
/* ---------------------------------------------------------------------------*/

cpl_image *
hdrl_sigclipfilter_image_grid(const cpl_image * ima, cpl_matrix * x, cpl_matrix * y,
                              cpl_size  filtersize_x, cpl_size filtersize_y)
{
    cpl_error_ensure(ima != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                    "NULL input image");
    cpl_error_ensure(filtersize_x > 0 && filtersize_y > 0 ,
                     CPL_ERROR_INCOMPATIBLE_INPUT, return NULL,
                                     "All function parameters must be greater then Zero");

    const cpl_size nx = cpl_image_get_size_x(ima);
    const cpl_size ny = cpl_image_get_size_y(ima);
    const cpl_size steps_x = cpl_matrix_get_nrow(x);
    const cpl_size steps_y = cpl_matrix_get_nrow(y);

    cpl_image * ima_local = cpl_image_new(steps_x, steps_y, CPL_TYPE_DOUBLE);
    cpl_image_get_bpm(ima_local);

#pragma omp parallel for
    for (cpl_size iy = 0; iy < steps_y; iy++) {
        cpl_size middlep_y = cpl_matrix_get(y, iy, 0);
        for (cpl_size ix = 0; ix < steps_x; ix++) {
            cpl_image * ima_cut;
            hdrl_image * ima_cut_hdrl;
            cpl_size middlep_x = cpl_matrix_get(x, ix, 0);

            cpl_size lowerlimit_x = CX_MAX(middlep_x - filtersize_x, 1);
            cpl_size lowerlimit_y = CX_MAX(middlep_y - filtersize_y, 1);
            cpl_size upperlimit_x = CX_MIN(middlep_x + filtersize_x, nx);
            cpl_size upperlimit_y = CX_MIN(middlep_y + filtersize_y, ny);

            ima_cut = cpl_image_extract(ima, lowerlimit_x, lowerlimit_y,
                                        upperlimit_x, upperlimit_y);
            ima_cut_hdrl = hdrl_image_create(ima_cut, NULL);

            hdrl_value median = hdrl_image_get_sigclip_mean(ima_cut_hdrl, 2.5 ,2.5 ,3);

            cpl_image_set(ima_local, ix + 1, iy + 1, median.data);
            /* reject also if less than one 4th bad pixels as in original
             * catalogue code */
            if (isnan(median.data) ||
                cpl_image_count_rejected(ima_cut) >= 0.25 * 2 * (filtersize_x *
                                                                 filtersize_y)) {
                cpl_image_reject(ima_local, ix + 1, iy + 1);

            }
            cpl_image_delete(ima_cut);
            hdrl_image_delete(ima_cut_hdrl);
        }
    }
    return ima_local;
}

/**@}*/
