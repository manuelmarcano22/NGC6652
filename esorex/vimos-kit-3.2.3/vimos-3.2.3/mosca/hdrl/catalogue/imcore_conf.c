/* $Id: imcore_conf.c,v 1.4 2015/08/12 11:16:55 jim Exp $
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
 * $Revision: 1.4 $
 * $Name:  $
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <cpl.h>
//#include "../casu_fits.h"
//#include "../casu_wcsutils.h"
#include "casu_utilfunctions.h"

#include "ap.h"
#include "util.h"
#include "imcore.h"
#include "floatmath.h"
#include "imcore_version.h"

#define FATAL_ERR(_a) {freetable(tab); cpl_msg_error(fctid,_a); tidy(); return(CASU_FATAL);}

#define NW 5

static float *g_smoothed = NULL;
static float *g_smoothedc = NULL;
static unsigned char *g_mflag = NULL;
static float *g_indata = NULL;
static float *g_confdata = NULL;
static float *g_confsqrt = NULL;
static ap_t g_ap;
static intptr_t g_freeconf = 0;

static float g_weights[NW*NW];
static long g_nx;
static long g_ny;

static void crweights(float);
static void convolve(intptr_t);
static void tidy(void);

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do source extraction
  
    \par Name:
        imcore_conf
    \par Purpose:
        Do source extraction
    \par Description:
        The main driving routine for the imcore source extraction program.
    \par Language:
        C
    \param infile
        The input image
    \param conf
        The input confidence map
    \param ipix
        The minimum allowable size of an object
    \param threshold
        The detection threshold in sigma above sky
    \param icrowd
        If set then the deblending software will be used
    \param rcore
        The core radius in pixels
    \param nbsize
        The smoothing box size for background map estimation
    \param cattype
        The type of catalogue to be produced
    \param filtfwhm
        The FWHM of the smoothing kernel in the detection algorithm
    \param gain
        The header keyword with the gain in e-/ADU
    \param outcat
        The output table of object
    \param res
        TODO
    \retval CASU_OK 
        if everything is ok
    \retval CASU_WARN,CASU_FATAL
        errors in the called routines
    \par QC headers:
        The following values will go into the table extension propertylist
        - \b SATURATION
            Saturation level in ADU
        - \b MEAN_SKY
            Mean sky brightness in ADU
        - \b SKY_NOISE
            Pixel noise at sky level in ADU
    \par DRS headers:
        The following values will go into the image extension propertylist
        - \b SKYLEVEL
            Mean sky brightness in ADU
        - \b SKYNOISE
            Pixel noise at sky level in ADU
        The following values will go into the table extension propertylist
        - \b THRESHOL
            The detection threshold in ADU
        - \b MINPIX 
            The minimum number of pixels per image
        - \b CROWDED
            Flag for crowded field analysis
        - \b RCORE
            The core radius for default profile fit in pixels
        - \b FILTFWHM
            The FWHM of the smoothing kernel in the detection algorithm
        - \b SEEING
            The average FWHM of stellar objects in pixels
        - \b XCOL
            The column containing the X position
        - \b YCOL
            The column containing the Y position
        - \b NXOUT
            The X dimension of the original image array
        - \b NYOUT
            The Y dimension of the original image array
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_conf(casu_fits *infile, casu_fits *conf, intptr_t ipix, 
                       float threshold, intptr_t icrowd, float rcore,
                       int bkg_subtr, intptr_t nbsize,
                       int cattype, float filtfwhm, float gain,
                       float saturation, hdrl_imcore_result * res) {

    intptr_t i,retval,mulpix,j,nw2,nobjects,imcore_xcol,imcore_ycol;
    int status;
    float fconst,nullval,skymed,skysig,thresh,xintmin,offset;
    float *current,*currentc;
    long npix,nxc,nyc,npts;
    cpl_image *map,*cmap;
    cpl_propertylist *plist,*extra;
    cpl_table *tab;
    char card[64];
    const char *fctid = "imcore_conf";

    /* Initialise output */

    res->catalogue = NULL;

    /* Useful constants */

    fconst = CPL_MATH_LOG2E;
    nullval = 0.0;
    nobjects = 0;

    /* Open input image */
/* CONSTANTS: 100, 1.5, 2, 8 , 3/4, 3/8*/
    tab = NULL;
    map = casu_fits_get_image(infile);
    if ((g_indata = cpl_image_get_data_float(map)) == NULL) 
        FATAL_ERR("Error getting image data");
    g_nx = (long)cpl_image_get_size_x(map);
    g_ny = (long)cpl_image_get_size_y(map);
    npts = g_nx*g_ny;

    /* Open the associated confidence map, if it exists */

    if (conf != NULL) {
        cmap = casu_fits_get_image(conf);
        if ((g_confdata = cpl_image_get_data(cmap)) == NULL)
            FATAL_ERR("Error getting confidence map data");
        nxc = (long)cpl_image_get_size_x(cmap);
        nyc = (long)cpl_image_get_size_y(cmap);
        if ((g_nx != nxc) || (g_ny != nyc))
            FATAL_ERR("Input image and confidence dimensions don't match");
        g_freeconf = 0;
    } else {
        g_confdata = cpl_malloc(npts*sizeof(*g_confdata));
        for (i = 0; i < npts; i++) 
            g_confdata[i] = 100;
        g_freeconf = 1;
        cmap = NULL;
    }
    
    /* Get g_mflag array for flagging saturated pixels */

    npix = g_nx*g_ny;
    g_mflag = cpl_calloc(npix,sizeof(*g_mflag));

    /* Open the g_ap structure and define some stuff in it */

    g_ap.lsiz = g_nx;
    g_ap.csiz = g_ny;
    g_ap.inframe = map;
    g_ap.conframe = cmap;
    g_ap.xtnum = casu_fits_get_nexten(infile);
    imcore_apinit(&g_ap);
    g_ap.indata = g_indata;
    g_ap.confdata = g_confdata;
    g_ap.multiply = 1;
    g_ap.ipnop = ipix;
    g_ap.mflag = g_mflag;
    g_ap.rcore = rcore;
    g_ap.filtfwhm = filtfwhm;
    g_ap.icrowd = icrowd;
    g_ap.fconst = fconst;

    /* Open the output catalogue FITS table */

    imcore_tabinit(&g_ap,&imcore_xcol,&imcore_ycol,cattype,&tab, res);

    //cpl_image_save(res->segmentation_map,"segmap.fits",CPL_TYPE_INT,NULL,CPL_IO_CREATE);
    /* Set up the data flags */

    for (i = 0; i < npix; i++) 
        if (g_confdata[i] == 0)
            g_mflag[i] = MF_ZEROCONF;
        else if (g_indata[i] < STUPID_VALUE)
            g_mflag[i] = MF_STUPID_VALUE;
        else 
            g_mflag[i] = MF_CLEANPIX;

    /* Flag up regions where the value is above the saturation level*/

    for (i = 0; i < npix ; i++)
        if (g_mflag[i] == MF_CLEANPIX && g_indata[i] > saturation)
            g_mflag[i] = MF_SATURATED;

    /* Compute the background variation and remove it from the data*/

    retval = imcore_background(&g_ap,nbsize,nullval,bkg_subtr,res);
    if (retval != CASU_OK) 
        FATAL_ERR("Error calculating background");


    /* Compute background statistics */
    imcore_backstats(&g_ap,nullval,&skymed,&skysig);

    if (retval != CASU_OK) 
        FATAL_ERR("Error calculating background stats");

    /* Take mean sky level out of data. */

    if (bkg_subtr) {
        for (i = 0; i < g_nx*g_ny; i++) {
            g_indata[i] -= skymed;
        }
    }
    /* Work out isophotal detection threshold levels */

    thresh = threshold*skysig;
    if (!bkg_subtr && thresh < skymed) {
        FATAL_ERR("Bad background corrected input. Background estimation "
                  "disabled but image median larger than threshold * sigma.");
    }
    
    /* Minimum intensity for consideration */

    xintmin = 1.5*thresh*((float)ipix);

    /* Minimum size for considering multiple images */

    mulpix = MAX(8,2*ipix);

    /* Actual areal profile levels: T, 2T, 4T, 8T,...but written wrt T
       i.e. threshold as a power of 2 */

    offset = logf(thresh)*fconst;

    /* Get a bit of workspace for buffers */

    g_smoothed = cpl_malloc(g_nx*sizeof(*g_smoothed));
    g_smoothedc = cpl_malloc(g_nx*sizeof(*g_smoothedc));

    /* Define a few things more things in g_ap structure */

    g_ap.mulpix = mulpix;
    g_ap.areal_offset = offset;
    g_ap.thresh = thresh;
    g_ap.xintmin = xintmin;

    g_ap.sigma = skysig;
    if (bkg_subtr) {
        g_ap.background = skymed;
        g_ap.saturation = (float)saturation-skymed;
    }
    else {
        g_ap.background = 0.;
        g_ap.saturation = (float)saturation;
    }

    /* Set the weights */

    crweights(filtfwhm);
    nw2 = NW/2;
    /* sqrt confdata, need as many rows as the convolution is wide */
    g_confsqrt = cpl_malloc(g_nx * NW * sizeof(*g_confsqrt));
    for (j = 0; j < NW; j++) {
        for (i = 0; i < g_nx; i++)
            g_confsqrt[j * g_nx + i] =
                sqrt(0.01*(float)g_confdata[j * g_nx + i]);
    }

    /* Right, now for the extraction loop.  Begin by defining a group of
       three rows of data and confidence */

    for (j = nw2; j < g_ny-nw2; j++) {
        current = g_indata + j*g_nx;
        if (j != nw2) {
            /* rotate buffer, could be a more efficient structure, but this
             * sufficent for now */
            memmove(g_confsqrt, g_confsqrt + g_nx,
                    g_nx * (NW - 1) * sizeof(*g_confsqrt));
            /* fill last row of buffer */
            for (i = 0; i < g_nx; i++) {
                g_confsqrt[(NW - 1) * g_nx + i] =
                    sqrt(0.01*(float)g_confdata[(j + nw2) * g_nx + i]);
            }
        }
        /* current row is the center of the buffer */
        currentc = g_confsqrt + nw2 * g_nx;
        convolve(j);
   
        /* Do the detection now */

        imcore_apline(&g_ap,current,currentc,g_smoothed,g_smoothedc,j,NULL);

        /* Make sure we're not overruning the stacks */

        if (g_ap.ibstack > (g_ap.maxbl - g_ap.lsiz))
            imcore_apfu(&g_ap);
        if (g_ap.ipstack > (g_ap.maxpa*3/4))
            for (i = 0; i < g_ap.maxpa*3/8; i++)
                imcore_apfu(&g_ap);

        /* See if there are any images to terminate */

        if (g_ap.ipstack > 1)
            imcore_terminate(&g_ap,cattype,gain,&nobjects,tab,res);
    }

    /* Post process. First truncate the cpl_table to the correct size and then
       work out an estimate of the seeing */

    cpl_table_set_size(tab,nobjects);
    retval = imcore_do_seeing(&g_ap,cattype,nobjects,tab);
    if (retval != CASU_OK)
        FATAL_ERR("Error working out seeing");
    imcore_tabclose(&g_ap,cattype);

    /* Create a property list with extra parameters. First QC parameters */

    extra = cpl_propertylist_duplicate(casu_fits_get_ehu(infile));
    cpl_propertylist_update_float(extra,"ESO QC SATURATION",g_ap.saturation);
    cpl_propertylist_set_comment(extra,"ESO QC SATURATION",
                                 "[adu] Saturation level");
    cpl_propertylist_update_float(extra,"ESO QC MEAN_SKY",g_ap.background);
    cpl_propertylist_set_comment(extra,"ESO QC MEAN_SKY",
                                 "[adu] Median sky brightness");
    cpl_propertylist_update_float(extra,"ESO QC SKY_NOISE",g_ap.sigma);
    cpl_propertylist_set_comment(extra,"ESO QC SKY_NOISE",
                                 "[adu] Pixel noise at sky level");

    /* Now DRS parameters */

    cpl_propertylist_update_float(extra,"ESO DRS THRESHOL",g_ap.thresh);
    cpl_propertylist_set_comment(extra,"ESO DRS THRESHOL",
                                 "[adu] Isophotal analysis threshold");
    cpl_propertylist_update_int(extra,"ESO DRS MINPIX",g_ap.ipnop);
    cpl_propertylist_set_comment(extra,"ESO DRS MINPIX",
                                 "[pixels] Minimum size for images");
    cpl_propertylist_update_int(extra,"ESO DRS CROWDED",g_ap.icrowd);
    cpl_propertylist_set_comment(extra,"ESO DRS CROWDED",
                                 "Crowded field analysis flag");
    cpl_propertylist_update_float(extra,"ESO DRS RCORE",g_ap.rcore);
    cpl_propertylist_set_comment(extra,"ESO DRS RCORE",
                                 "[pixels] Core radius for default profile fit");
    cpl_propertylist_update_float(extra,"ESO DRS SEEING",g_ap.fwhm);
    cpl_propertylist_set_comment(extra,"ESO DRS SEEING",
                                 "[pixels] Average FWHM");
    cpl_propertylist_update_float(extra,"ESO DRS FILTFWHM",g_ap.filtfwhm);
    cpl_propertylist_set_comment(extra,"ESO DRS FILTFWHM",
                                 "[pixels] FWHM of smoothing kernel");
    cpl_propertylist_update_int(extra,"ESO DRS XCOL",imcore_xcol);
    cpl_propertylist_set_comment(extra,"ESO DRS XCOL","Column for X position");
    cpl_propertylist_update_int(extra,"ESO DRS YCOL",imcore_ycol);
    cpl_propertylist_set_comment(extra,"ESO DRS YCOL","Column for Y position");
    cpl_propertylist_update_int(extra,"ESO DRS NXOUT",g_nx);
    cpl_propertylist_set_comment(extra,"ESO DRS NXOUT",
                                 "X Dimension of input image");
    cpl_propertylist_update_int(extra,"ESO DRS NYOUT",g_ny);
    cpl_propertylist_set_comment(extra,"ESO DRS NYOUT",
                                 "Y Dimension of input image");
    snprintf(card,64,"IMCORE version: %s",imcore_version);
    cpl_propertylist_append_string(extra,"HISTORY",card);
    
    /* Now wrap all this stuff up and send it back */

    plist = cpl_propertylist_duplicate(casu_fits_get_phu(infile));
    status = CASU_OK;
    (void)casu_tabwcs(extra,imcore_xcol,imcore_ycol,&status);
    res->catalogue = casu_tfits_wrap(tab,NULL,plist,extra);

    /* Tidy and exit */  

    tidy();
    return(CASU_OK);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        crweights
    \par Purpose:
        Create convolution kernel g_weights
    \par Description:
        The convolution kernel for a Gaussian with a given FWHM is created
    \par Language:
        C
    \param filtfwhm
        The FWHM of the Gaussian kernel
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void crweights(float filtfwhm) {
    intptr_t i,j,nw2,n;
    double gsigsq,di,dj;
    float renorm;

    /* Get the kernel size */

    nw2 = NW/2;
    
    /* Set the normalisation constants */
/* 2, 2.35 one should use CPL_MATH_FWHM_SIG */
    gsigsq = 1.0/(2.0*pow(MAX(1.0,(double)filtfwhm)/2.35,2.0));
    renorm = 0.0;

    /* Now work out the weights */

    n = -1;
    for (i = -nw2; i <= nw2; i++) {
        di = (double)i;
        di *= gsigsq*di;
        for (j = -nw2; j <= nw2; j++) {
            dj = (double)j;
            dj *= gsigsq*dj;
            n++;
            g_weights[n] = (float)exp(-(di + dj));
            renorm += g_weights[n];
        }
    }

    /* Now normalise the weights */

    n = -1;
    for (i = -nw2; i <= nw2; i++) {
        for (j = -nw2; j <= nw2; j++) {
            n++;
            g_weights[n] /= renorm;
        }
    }
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        convolve
    \par Purpose:
        Smooth the original data and confidence
    \par Description:
        Smooth a line of original data and confidence by convolving with a
        Gaussian kernel.
    \par Language:
        C
    \param ir
        The row number of the input image to smooth
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void convolve(intptr_t ir) {
    intptr_t i,nw2,ix,jx,jy,n;
    float *idata,*cdata;

    /* Zero the summations */

    for (i = 0; i < g_nx; i++) {
        g_smoothed[i] = 0.0;
        g_smoothedc[i] = 0.0;
    }

    /* Now big is the smoothing kernel? */

    nw2 = NW/2;

    /* Now loop for each column */

    for (ix = nw2; ix < g_nx-nw2; ix++) {
        n = -1;
        float sum = 0., sumc = 0.;
        for (jy = ir-nw2; jy <= ir+nw2; jy++) {
            idata = g_indata + jy*g_nx;
            /* g_confsqrt [0..NW] rows */
            cdata = g_confsqrt + (jy - ir + nw2) * g_nx;
            for (jx = ix-nw2; jx <= ix+nw2; jx++) {
                n++;
                sum += g_weights[n] * idata[jx];
                sumc += g_weights[n] * idata[jx] * cdata[jx];
            }
        }
        g_smoothed[ix] = sum;
        g_smoothedc[ix] = sumc;
    }
}

static void tidy(void) {

    if (g_freeconf) 
        freespace(g_confdata);
    freespace(g_confsqrt);
    freespace(g_smoothed);
    freespace(g_smoothedc);
    freespace(g_mflag);
    imcore_apclose(&g_ap);
}

/**@}*/

/* 

$Log: imcore_conf.c,v $
Revision 1.4  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.3  2015/08/07 13:06:54  jim
Fixed copyright to ESO

Revision 1.2  2015/08/06 05:34:02  jim
Fixes to get rid of compiler moans

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.5  2015/01/09 11:42:36  jim
Fixed routines to remove globals

Revision 1.4  2014/12/11 12:23:34  jim
new version

Revision 1.3  2014/04/09 09:09:51  jim
Detabbed

Revision 1.2  2014/03/26 15:25:19  jim
Modified for floating point confidence maps

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
