/* $Id: classify.c,v 1.3 2015/08/12 11:16:55 jim Exp $
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>

#include <cpl.h>
//#include <../casu_utils.h>
#include "casu_utilfunctions.h"

#include "classify.h"

#define MAXHIST   111
#define STEP      0.05
#define NSAMPLE   150
#define MAXLOOP   5
#define BLIMDEF   15.0;
#define FLIMDEF   11.0;
#define CMINDEF   7.5
#define CMAXDEF   15.0
#define NAREAL    8
#define PI        3.14159265358979323846
#define DEGRAD 57.2957795130823229
/* CONSTANT: use rather CPL_MATH_DEG_RAD,CPL_MATH_PI */
/* the term 2.5*log10(value) comes from the formula to compute a magnitude */
#define COREMAG(A,B,C) 2.5*log10((double)(max(C,A-B)))
#define MAX(A,B) (A > B ? A : B)
#define NCOLS 63

/* Make the data arrays and header values global */

static long g_nrows;
static float g_thresh, g_skylevel, g_skynoise, g_rcore;

/* Derived values */

static intptr_t g_poor;
static float g_sigell, g_fitell, g_elllim, g_sigellf, g_fitellf, g_sigpa, g_fitpa;
static float g_blim, g_flim, g_cmin, g_cmax;
static float g_fit1, g_fit2, g_fit3, g_fit4, g_fit5, g_fit6, g_fit7;
static float g_fit_final, g_sigma_final;
static float *g_lower1, *g_lower2, *g_lower3, *g_upper1, *g_upper2, *g_upper3, *g_uppere;
static float g_avsig1, g_avsig2, g_avsig3, g_wt1, g_wt2, g_wt3;

/* Classification values */

static intptr_t g_nstar, g_ngal, g_njunk, g_ncmp;

/* Values for the data quality and aperture corrections */

static float g_avsat, g_corlim, g_cormin, g_apcpkht;
static float g_apcor,  g_apcor1, g_apcor2, g_apcor3, g_apcor4, g_apcor5;
static float g_apcor6, g_apcor7;

/* Data arrays */

static float *g_workspace = NULL;
static cpl_table *g_catcopy = NULL;
static float *g_areal[NAREAL];
static float *g_core_flux,  *g_core1_flux, *g_core2_flux, *g_core3_flux;
static float *g_core4_flux, *g_core5_flux, *g_core6_flux;
static float *g_peak_height, *g_peak_mag, *g_ellipticity, *g_iso_flux;
static float *g_total_flux, *g_cls, *g_sig, *g_xpos, *g_ypos, *g_pa, *g_skylev;


/* Column definitions */

#define NCOL32 14
#define NCOLFULL 15

static const char *g_colsfull[NCOLFULL] = {"Aper_flux_3","Aper_flux_1","Aper_flux_4",
                                     "Aper_flux_5","Aper_flux_6","Peak_height",
                                     "Ellipticity","Isophotal_flux",
                                     "Isophotal_flux","Aper_flux_7",
                                     "X_coordinate","Y_coordinate",
                                     "Position_angle","Sky_level",
                                     "Aper_flux_2"};

static intptr_t g_ncols;
static float g_xmin;
static float g_xmax;
static float g_ymin;
static float g_ymax;
static float g_pixlim;

#define FRAMECUT 0.05

/* Subroutine prototypes */

static void anhist(float *, intptr_t, float *, float *);
static void boundaries(float *, float *, float *, float, float, float, float,
                       intptr_t, float, float, float *, float *, float *,
		       float *);
static void boundpk(float *, float *, float, float, float *, float *,
                    float *, float *);
static void classify_run(void);
static void classstats(float *, float *, intptr_t, float, float *, float *);
static void classstats_ap0(float *, float *);
static void classstats_ap67(float *, float *, float *, float *);
static void classstats_el(void);
static void classstats_pa(void);
static void classstats_ellf(float);
static void classstats_final(void);
static void medstat(float *, intptr_t, float *, float *);
static void sort1(float *, intptr_t);
static void sort2(float *, float *, intptr_t);

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do star/galaxy classification
  
    \par Name:
        imcore_classify
    \par Purpose:
        Do star/galaxy classification from a imcore catalogue
    \par Description:
        The information in an imcore catalogue is scanned and each object
        is classified based on a number of shape criteria.
    \par Language:
        C
    \param catalogue
        The input imcore catalogue
    \param plist
        The propertylist from the primary header of the input image
    \param minsize
        The minimum size in pixels of objects to be used in the analysis
    \param cattype
        The catalogue type
    \retval CASU_OK
        If everything is OK
    \retval CASU_FATAL
        If the input catalogue is unrecognised.
    \par QC headers:
        The following QC parameters are read from the catalogue extension
        header propertylist
        - \b MEAN_SKY
            The mean sky found by imcore
        - \b SKY_NOISE
            The sky noise found by imcore
        The following QC parameters are written to the catalogue extension
        header propertylist
        - \b IMAGE_SIZE
            The average FWHM of stellar objects in the catalogue
        - \b ELLIPTICITY
            The average stellar ellipticity
        - \b POSANG
            The average position angle on the image
        - \b APERTURE_CORR
            The aperture correction for an aperture with a radius of Rcore.
        - \b NOISE_OBJ
            The number of noise objects detected on the image
    \par DRS headers:
        The following DRS parameters are read from the catalogue extension
        header propertylist
        - \b THRESHOL
            The detection threshold used by imcore
        - \b RCORE
            The core radius used by imcore
        - \b SEEING
            The averaged seeing found by imcore
        - \b NXOUT
            The number of pixels in a row of the original image
        - \b NYOUT
            The number of pixels in a column of the original image
        The following DRS parameters are written to the catalogue extension
        header propertylist
        - \b CLASSIFD 
            Set if the catalogue has been classified
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_classify(casu_tfits *catalogue, cpl_propertylist *plist, 
                           float minsize, int cattype) {
    float fwhm,*work,moff;
    float pkht,ell,core,ap,delap,area,junk,arg;
    char *cols[ MAX(NCOL32, NCOLFULL) ], colname[32];
    cpl_propertylist *extra;
    cpl_table *cat;
    intptr_t i,n,iap,nxout,nyout;

    /* Get some DQC info from the extra propertylist generated by imcore */

    extra = casu_tfits_get_ehu(catalogue);
    g_thresh = cpl_propertylist_get_float(extra,"ESO DRS THRESHOL");
    g_skylevel = cpl_propertylist_get_float(extra,"ESO QC MEAN_SKY");
    g_skynoise = cpl_propertylist_get_float(extra,"ESO QC SKY_NOISE");
    g_rcore = cpl_propertylist_get_float(extra,"ESO DRS RCORE");
    fwhm = cpl_propertylist_get_float(extra,"ESO DRS SEEING");
    nxout = cpl_propertylist_get_int(extra,"ESO DRS NXOUT");
    nyout = cpl_propertylist_get_int(extra,"ESO DRS NYOUT");
    g_xmin = FRAMECUT*(float)nxout;
    g_xmax = (1.0 - FRAMECUT)*(float)nxout;
    g_ymin = FRAMECUT*(float)nyout;
    g_ymax = (1.0 - FRAMECUT)*(float)nyout;
    g_pixlim = minsize;

    /* Get the number of columns and decide which column labels to use */

    cat = casu_tfits_get_table(catalogue);
    g_ncols = cpl_table_get_ncol(cat);
    for (i = 0; i < NCOLFULL; i++)
        cols[i] = (char *)g_colsfull[i];

    /* Make a copy of the table as you are going to muck about with the
       column values. Get the column data */

    g_catcopy = cpl_table_duplicate(cat);
    g_nrows = cpl_table_get_nrow(cat);
    g_core_flux = cpl_table_get_data_float(g_catcopy,cols[0]);
    g_core1_flux = cpl_table_get_data_float(g_catcopy,cols[1]);
    g_core2_flux = cpl_table_get_data_float(g_catcopy,cols[2]);
    g_core3_flux = cpl_table_get_data_float(g_catcopy,cols[3]);
    g_core4_flux = cpl_table_get_data_float(g_catcopy,cols[4]);
    g_peak_height = cpl_table_get_data_float(g_catcopy,cols[5]);
    g_ellipticity = cpl_table_get_data_float(g_catcopy,cols[6]);
    g_iso_flux = cpl_table_get_data_float(g_catcopy,cols[7]);
    g_total_flux = cpl_table_get_data_float(g_catcopy,cols[8]);
    g_core5_flux = cpl_table_get_data_float(g_catcopy,cols[9]);
    g_xpos = cpl_table_get_data_float(g_catcopy,cols[10]);
    g_ypos = cpl_table_get_data_float(g_catcopy,cols[11]);
    g_pa = cpl_table_get_data_float(g_catcopy,cols[12]);
    g_skylev = cpl_table_get_data_float(g_catcopy,cols[13]);
    g_core6_flux = cpl_table_get_data_float(g_catcopy,cols[14]);
    g_cls = cpl_table_get_data_float(cat,"Classification");
    g_sig = cpl_table_get_data_float(cat,"Statistic");

    /* Get some workspace */

    g_workspace = cpl_malloc(2*g_nrows*sizeof(float));
    g_peak_mag = g_workspace;
    work = g_workspace + g_nrows;
    
    /* Convert fluxes to "magnitudes" */

    /* CONSTANTS: 0.6 , 32, 5, sqrt(2), log(2), 2.0, 0.25 */
    for (i = 0; i < g_nrows; i++) {
        g_core_flux[i] = COREMAG(g_core_flux[i],0.0,1.0);
        g_core1_flux[i] = COREMAG(g_core1_flux[i],0.0,1.0);
        g_core2_flux[i] = COREMAG(g_core2_flux[i],0.0,1.0);
        g_core3_flux[i] = COREMAG(g_core3_flux[i],0.0,1.0);
        g_core4_flux[i] = COREMAG(g_core4_flux[i],0.0,1.0);
        g_core5_flux[i] = COREMAG(g_core5_flux[i],0.0,1.0);
        moff = 1.0/(1.0 - pow((g_thresh/MAX(g_peak_height[i],g_thresh)),0.6));
        g_iso_flux[i] = COREMAG(moff*g_iso_flux[i],0.0,1.0);
        g_peak_mag[i] = COREMAG(g_peak_height[i],g_skynoise,0.1);
    }
    if (g_core6_flux != NULL) 
        for (i = 0; i < g_nrows; i++)
            g_core6_flux[i] = COREMAG(g_core6_flux[i],0.0,1.0);
    if (g_ncols == 32)
        for (i = 0; i < g_nrows; i++)
            g_total_flux[i] = COREMAG(g_total_flux[i],0.0,1.0);

    /*  Now get the g_areal profile information. You'll need this in a sec */

    for (i = 0; i < NAREAL; i++) {
        sprintf(colname,"Areal_%zd_profile",i+1);
        g_areal[i] = cpl_table_get_data_float(g_catcopy,colname);
    }

    /* What is the seeing like? */

    g_poor = 0;
    if (fwhm > max(5.0,g_rcore*sqrt(2.0)))
        g_poor = 1;

    /* Ok, now call the routine that does all the work */

    classify_run();

    /* Right, now get a better estimate of the seeing */

    n = 0;
    for (i = 0; i < g_nrows; i++) {
        pkht = g_peak_height[i];
        ell = g_ellipticity[i];
        core = g_core_flux[i];
        if (g_cls[i] == -1.0 && ell < g_elllim && core < g_corlim && 
            pkht > 10.0*g_thresh) { 
            ap = log(0.5*pkht/g_thresh)/log(2.0) + 1.0;
            iap = (intptr_t)ap;
            delap = ap - (float)iap;
            if (iap > 0 && iap < NAREAL && g_areal[1][i] > 0.0) {
                area = g_areal[iap-1][i]*(1.0 - delap) + g_areal[iap][i]*delap;
                work[n++] = 2.0*sqrt(area/PI);
            }
        }
    }
    if (n > 2) { 
        medstat(work,n,&fwhm,&junk);
       
        /* Allow for finite pixel size */

        arg = 0.25*PI*fwhm*fwhm - 1;
        fwhm = 2.0*sqrt(max(0.0,arg/PI));
       
    } else
        fwhm = -1.0;

    /* Tidy up a bit */

    freespace(g_workspace);
    freetable(g_catcopy);

    /* Write header results into extra property list. First the QC */

    cpl_propertylist_update_float(extra,"ESO QC IMAGE_SIZE",fwhm);
    cpl_propertylist_set_comment(extra,"ESO QC IMAGE_SIZE",
                                 "[pixels] Average FWHM of stellar objects");
    cpl_propertylist_update_float(extra,"ESO QC ELLIPTICITY",g_fitell);
    cpl_propertylist_set_comment(extra,"ESO QC ELLIPTICITY",
                                 "Average stellar ellipticity (1-b/a)");
    cpl_propertylist_update_float(extra,"ESO QC POSANG",g_fitpa);
    cpl_propertylist_set_comment(extra,"ESO QC POSANG",
                                 "[degrees] Median position angle");
    cpl_propertylist_update_float(extra,"ESO QC APERTURE_CORR",g_apcor3);
    cpl_propertylist_set_comment(extra,"ESO QC APERTURE_CORR",
                                 "Stellar ap-corr 1x core flux");
    cpl_propertylist_update_int(extra,"ESO QC NOISE_OBJ",g_njunk);
    cpl_propertylist_set_comment(extra,"ESO QC NOISE_OBJ",
                                 "Number of noise objects");
    cpl_propertylist_update_float(extra,"ESO QC SATURATION",g_avsat);
    
    /* Now some helpful DRS keywords */

    cpl_propertylist_update_bool(extra,"ESO DRS CLASSIFD",1);
    cpl_propertylist_set_comment(extra,"ESO DRS CLASSIFD",
                                 "Catalogue has been classified");

    /* Now the aperture correction keywords */
    cpl_propertylist_update_float(extra,"APCORPK",g_apcpkht);
    cpl_propertylist_set_comment(extra,"APCORPK","Stellar aperture correction - peak height");
    cpl_propertylist_update_float(extra,"APCOR1",g_apcor1);
    cpl_propertylist_set_comment(extra,"APCOR1","Stellar aperture correction - 1/2x core flux");
    cpl_propertylist_update_float(extra,"APCOR2",g_apcor2);
    cpl_propertylist_set_comment(extra,"APCOR2","Stellar aperture correction - core/sqrt(2) flux");
    cpl_propertylist_update_float(extra,"APCOR3",g_apcor3);
    cpl_propertylist_set_comment(extra,"APCOR3","Stellar aperture correction - 1x core flux");
    cpl_propertylist_update_float(extra,"APCOR4",g_apcor4);
    cpl_propertylist_set_comment(extra,"APCOR4","Stellar aperture correction - sqrt(2)x core flux");
    cpl_propertylist_update_float(extra,"APCOR5",g_apcor5);
    cpl_propertylist_set_comment(extra,"APCOR5","Stellar aperture correction - 2x core flux");
    cpl_propertylist_update_float(extra,"APCOR6",g_apcor6);
    cpl_propertylist_set_comment(extra,"APCOR6","Stellar aperture correction - 2*sqrt(2)x core flux");
    cpl_propertylist_update_float(extra,"APCOR7",g_apcor7);
    cpl_propertylist_set_comment(extra,"APCOR7","Stellar aperture correction - 4x core flux");

    /* Write header information to help GAIA */

    cpl_propertylist_update_string(extra,"SYMBOL1","{Ellipticity Position_angle Areal_1_profile Classification} {el");
    cpl_propertylist_update_string(extra,"SYMBOL2","lipse blue (1.0-$Ellipticity) $Position_angle+90 {} $Classific");
    cpl_propertylist_update_string(extra,"SYMBOL3","ation==1} {sqrt($Areal_1_profile*(1.0-$Ellipticity)/3.142)} : {");
    cpl_propertylist_update_string(extra,"SYMBOL4","Ellipticity Position_angle Areal_1_profile Classification} {el");
    cpl_propertylist_update_string(extra,"SYMBOL5","lipse red (1.0-$Ellipticity) $Position_angle+90 {} $Classific");
    cpl_propertylist_update_string(extra,"SYMBOL6","ation==-1} {sqrt($Areal_1_profile*(1.0-$Ellipticity)/3.142)} :");
    cpl_propertylist_update_string(extra,"SYMBOL7","{Ellipticity Position_angle Areal_1_profile Classification} {el");
    cpl_propertylist_update_string(extra,"SYMBOL8","lipse green (1.0-$Ellipticity) $Position_angle+90 {} $Classifi");
    cpl_propertylist_update_string(extra,"SYMBOL9","cation==0} {sqrt($Areal_1_profile*(1.0-$Ellipticity)/3.142)}");

    /* Get out of here */

    return(CASU_OK);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        anhist
    \par Purpose:
        Analyse a histogram to give a median and sigma
    \par Description:
        The entries in a data array are histogrammed. The histogram is 
        analsyed to work out a median and a sigma. A certain amount of 
        smoothing is done and extra searches are also done to make sure 
        that they median position estimate is the best it can be.
    \par Language:
        C
    \param data
        The input data array
    \param n
        The size of the input data array
    \param medval
        The output median value
    \param sigma
        The output sigma value 
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void anhist(float *data, intptr_t n, float *medval, float *sigma) {
    intptr_t i,*histo,ilev,imax,ismax;
    float *sval,hmax,smax,hlim,ratio;

    /* CONSTANTS: -10, 100, 1.48, 0.5, 3.0, sqrt(2.) */
    /* Get some workspace for the histogram */

    histo = cpl_calloc(MAXHIST,sizeof(intptr_t));
    sval = cpl_calloc(MAXHIST,sizeof(float));

    /* Sort data into the histogram */

    for (i = 0; i < n; i++) {
        ilev = casu_nint(data[i]/STEP);
        if (ilev >= -10 && ilev <= 100) {
            ilev += 10;
            histo[ilev] += 1;
        }
    }

    /* Now find the maximum of the histogram and its position... */

    hmax = 0.0;
    imax = 0;
    for (i = 0; i < MAXHIST; i++) {
        if (histo[i] > hmax) {
            hmax = (float)histo[i];
            imax = i;
        }
    }

    /* Trap for hmax == 0 */

    if (hmax == 0.0) {
        if (n >= 10) {
            *medval = data[(n+1)/2-1];
            *sigma = CPL_MATH_STD_MAD*0.5*(data[(3*n+3)/4-1] - data[(n+3)/4-1]);
        } else {
            *medval = 0.0;
            *sigma = 1.0;
        }
        freespace(histo);
        freespace(sval);
        return;
    }

    /* Now do three point running average to see if there are other local
       maxima */

    smax = 0.0;
    ismax = 0;
    for (i = 1; i < MAXHIST-1; i++) {
        sval[i] = (histo[i-1] + histo[i] + histo[i+1])/3.0;
        if (sval[i] > smax) {
            smax = sval[i];
            ismax = i;
        }
    }
    if (ismax < imax) {
        imax = ismax;
        hmax = (float)histo[imax];
    }

    /* Now check for lower local maxima */

    for (i = imax-1; i > 0; i--) {
        if (sval[i] >= sval[i+1] && sval[i] >= sval[i-1]) {
            if (sval[i] > 0.5*smax)
                ismax = i;
        }
    }
    if (ismax < imax) {
        imax = ismax;
        hmax = (float)histo[imax];
    }
    
    /* Now work out where the peak is */
    
    *medval = min((float)(imax-10)*STEP,data[(n+1)/2-1]);
    hlim = casu_nint(0.5*hmax);
    i = 1;
    while (histo[imax-i] > hlim && imax-i > 1)
        i++;
    ratio = hmax/max(1.0,(float)histo[imax-i]);
    *sigma = (float)i*STEP/(sqrt(2.0)*max(1.0,log(ratio)));
    *sigma = max(*sigma,0.5*STEP);

    /* Tidy and exit */
    
    freespace(histo);
    freespace(sval);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        boundaries
    \par Purpose:
        Work out boundaries of the stellar locus 
    \par Description:
        A number of flux estimates are given along with comparison statistics
        between those fluxes.
    \par Language:
        C
    \param core1
        The flux array for the primary aperture (in magnitudes)
    \param core2
        The flux array for the secondary aperture in the case of good seeing
    \param core3
        The flux array for the secondary aperture in the case of poor seeing
    \param medval1
        The median of magnitude difference of core1 and core2
    \param sigma1
        The sigma of magnitude difference of core1 and core2
    \param medval2
        The median of magnitude difference of core1 and core3
    \param sigma2
        The sigma of magnitude difference of core1 and core3
    \param small
        This parameter changes sign in internal computation
    \param area1
        The area of the larger of the apertures for core1,core2 comparison
    \param area2
        The area of the larger of the apertures for core1,core3 comparison
    \param wt
        An output weight for the estimate based on the scatter in the
        magnitude differences
    \param avsig
        An output average magnitude difference for the apertures used.
    \param lower
        An array delimiting the lower boundary
    \param upper
        An array delimiting the upper boundary
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void boundaries(float *core1, float *core2, float *core3, float medval1,
                       float sigma1, float medval2, float sigma2, intptr_t small, 
                       float area1, float area2, float *wt, float *avsig, 
                       float *lower, float *upper) {
    intptr_t i,n;
    float c1,c2,dc,*work,xnoise,xmag,xflux,ratio,asign,junk;

    /* Get a workspace */

    work = cpl_malloc(g_nrows*sizeof(float));

    /* Initialise the lower boundary */

    lower[0] = g_cmin;
    lower[1] = g_cmax;
    asign = ((small == 1) ? -1.0 : 1.0);
    
    /* Now collect the data */
    /* CONSTANTS: 3.0, 2.5, 0.25, 4 */
    n = 0;
    for (i = 0; i < g_nrows; i++) {
        c1 = core1[i];
        if (! g_poor) {
            c2 = core2[i];
            dc = asign*(c2 - c1);
            if (dc > medval1 - 3.0*sigma1 && c1 < g_blim - 3.0)
                work[n++] = dc - medval1;
        } else {
            c2 = core3[i];
            dc = c2 - c1;
            if (dc > medval2 - 3.0*sigma2 && c1 < g_blim - 3.0)
                work[n++] = dc - medval2;
        }
    }
 
    /* Find the median */

    medstat(work,n,avsig,&junk);
    freespace(work);

    /* Work out sigma levels for both types of seeing */

    if (! g_poor) {
        *wt = min(5.0,max(1.0,*avsig/sigma1));
        xnoise = sqrt(area1)*g_skynoise;
    } else {
        *wt = min(2.5,max(1.0,*avsig/sigma2));
        xnoise = sqrt(area2)*g_skynoise;
    }

    /* Now work out the boundaries */
    /* The term pow(10.0,(double)(0.4*xmag)); comes from the magitude formula */
    for (i = 0; i < NSAMPLE; i++) {
        xmag = 5.0 + (float)(i+1)*0.1;
        xflux = pow(10.0,(double)(0.4*xmag));
        ratio = COREMAG(1.0+xnoise/xflux,0.0,0.0);
        if (! g_poor) {
            lower[i] = medval1 - 3.0*sqrt(sigma1*sigma1 + ratio*ratio);
            upper[i] = medval1 + 3.0*sqrt(sigma1*sigma1 + 0.5*ratio*ratio);
        } else {
            lower[i] = medval2 - 3.0*sqrt(sigma2*sigma2 + ratio*ratio);
            upper[i] = medval2 + 3.0*sqrt(sigma2*sigma2 + 0.5*ratio*ratio);
        }
    }
    upper[0] = ((g_poor == 0) ? medval1 : medval2);
    upper[1] = upper[0];
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        boundpk
    \par Purpose:
        Work out boundaries of the stellar locus using peak flux
    \par Description:
        The boundaries of the stellar locus are located using the core
        flux and the peak flux.
    \par Language:
        C
    \param core
        The flux array for the primary aperture (in magnitudes)
    \param pkht
        The array of magnitudes based on the peak height
    \param medval
        The median of magnitude difference of core and pkht
    \param sigma
        The sigma of magnitude difference of core and pkht
    \param wt
        An output weight for the estimate based on the scatter in the
        magnitude differences
    \param avsig
        An output average magnitude difference for the two estimates
    \param lower
        An array delimiting the lower boundary
    \param upper
        An array delimiting the upper boundary
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void boundpk(float *core, float *pkht, float medval, float sigma, 
                    float *wt, float *avsig, float *lower, float *upper) {
    intptr_t i,n;
    float c,p,*work,xnoise,xmag,pmag,xflux,pflux,ratio,junk;

    /* Get the space for the boundry lines and a workspace */

    work = cpl_malloc(g_nrows*sizeof(float));

    /* Collect the data */

    n = 0;
    for (i = 0; i < g_nrows; i++) {
        c = core[i];
        p = pkht[i];
        if (c - p > medval - 3.0*sigma && c < g_blim - 3.0)
            work[n++] = c - p - medval;
    }

    /* Find the median */

    medstat(work,n,avsig,&junk);
    freespace(work);
    *wt = min(5.0,max(1.0,*avsig/sigma));

    /* Now work out boundaries */
    /* CONSTANTS: 3, 5, 0.4, */
    /* The term pow(10.0,(double)(0.4*xmag)) comes from the magnitude formula */
    /* The term 2.5*log10 comes from the magnitude formula */
    xnoise = sqrt(PI*g_rcore*g_rcore)*g_skynoise;
    for (i = 0; i < NSAMPLE; i++) {
        xmag = 5.0 + (float)(i+1)*0.1;
        pmag = xmag - medval;
        xflux = pow(10.0,(double)(0.4*xmag));
        pflux = pow(10.0,(double)(0.4*pmag));
        ratio = 2.5*log10((double)(1.0+max(xnoise/xflux,g_skynoise/pflux)));
        lower[i] = medval - 3.0*sqrt(sigma*sigma + ratio*ratio);
        upper[i] = medval + 3.0*sqrt(sigma*sigma + 0.5*ratio*ratio);
    }
    upper[0] = medval;
    upper[1] = upper[0];
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        classify_run
    \par Purpose:
        Main driver routine
    \par Description:
        This is the main driver routine for classify. It calls all the
        statistical routines and the boundary finding routines. It works 
        out the aperture corrections and does the final classification for
        all objects in the catalogue.
    \par Language:
        C
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void classify_run(void) {
    float fluxlim,ell,pkht,core,sig1,sig2,sig3,denom,w1,w2,w3;
    float core_small,core_large,core_midd,statistic,statcut,sigtot;
    float fit0,sigma0,xnoise,xmag,ratio,xflux,ratell,ratscl,ellbound;
    float *lower,*upper,sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7;
    float *work,avsatnew,junk;
    intptr_t i,iarg,ii;

    /* Update faint limit to cope with short exposures */

    g_blim = BLIMDEF;
    g_flim = FLIMDEF;
    /* the following formula comes from the magnitude computation */
    fluxlim = 2.5*log10((double)(5.0*sqrt(PI*g_rcore*g_rcore)*g_skynoise));
    g_flim = min(g_flim,max(6.0,fluxlim+3.0));
    g_corlim = min(g_blim,max(12.5,fluxlim+5.0));
    g_cormin = min(g_blim,max(12.5,fluxlim+5.0));
    /* CONSTANTS: 2.5, 5, PI, 6.0, 12.5, 0.5, 20, 0.2 */
    /* Work out min and max core flux */

    g_cmin = CMINDEF;
    g_cmax = CMAXDEF;
    for (i = 0; i < g_nrows; i++) {
        xflux = g_core_flux[i];
        g_cmin = min(g_cmin,xflux);
        g_cmax = max(g_cmax,xflux);
    }
    g_cmin = max(fluxlim-0.5,g_cmin);
    g_cmax += 0.1;
    g_cmax = min(g_cmax,20.0);

    /* Work out g_ellipticity stats for likely stellar objects */

    classstats_el();

    /* Ok, get the classification statistics for each of the tests.  First
       the core flux vs 1/2*core flux */

    classstats(g_core_flux,g_core1_flux,1,0.2,&g_fit1,&sigma1);

    /* Core flux vs 2*core flux */

    classstats(g_core_flux,g_core3_flux,0,0.1,&g_fit2,&sigma2);

    /* Core flux vs sqrt(2)*core flux */

    classstats(g_core_flux,g_core2_flux,0,0.0,&g_fit4,&sigma4);

    /* Core flux vs 2*sqrt(2)*core flux */

    classstats(g_core_flux,g_core4_flux,0,0.1,&g_fit5,&sigma5);

    /* Core flux vs Peak height */

    classstats(g_core_flux,g_peak_mag,1,0.2,&g_fit3,&sigma3);

    /* Faint end g_ellipticity */

    classstats_ellf(fluxlim);

    /* Work out position angle stats for likely stellar objects */

    classstats_pa();

    /* Get workspace for the boundary arrays */
   
    g_lower1 = cpl_malloc(NSAMPLE*sizeof(float));
    g_lower2 = cpl_malloc(NSAMPLE*sizeof(float));
    g_lower3 = cpl_malloc(NSAMPLE*sizeof(float));
    g_upper1 = cpl_malloc(NSAMPLE*sizeof(float));
    g_upper2 = cpl_malloc(NSAMPLE*sizeof(float));
    g_upper3 = cpl_malloc(NSAMPLE*sizeof(float));

    /* Right, work out the boundaries for the classification tests 
       First core vs sqrt(2)*core or core vs 0.5*core depending upon
       the seeing */

    boundaries(g_core_flux,g_core1_flux,g_core2_flux,g_fit1,sigma1,g_fit4,sigma4,
               1,PI*g_rcore*g_rcore,2.0*PI*g_rcore*g_rcore,&g_wt1,&g_avsig1,g_lower1,
               g_upper1);

    /* Now core vs 2*core or core vs 2*sqrt(2)*core */

    boundaries(g_core_flux,g_core3_flux,g_core4_flux,g_fit2,sigma2,g_fit5,sigma5,
               0,4.0*PI*g_rcore*g_rcore,8.0*PI*g_rcore*g_rcore,&g_wt2,&g_avsig2,g_lower2,
               g_upper2);

    /* Now core vs peak height */

    boundpk(g_core_flux,g_peak_mag,g_fit3,sigma3,&g_wt3,&g_avsig3,g_lower3,g_upper3); 

     
    /* Do final classification statistics and find the saturation limit */

    classstats_final();

    /* Define final boundaries */
    /* The term pow(10.0,(0.4*fluxlim+1.5)); comes from the magitude formula */
    lower = cpl_malloc(NSAMPLE*sizeof(float));
    upper = cpl_malloc(NSAMPLE*sizeof(float));
    g_uppere = cpl_malloc(NSAMPLE*sizeof(float));
    xnoise = sqrt(PI*g_rcore*g_rcore)*g_skynoise;
    ratell = xnoise/pow(10.0,0.4*(fluxlim+1.5));
    ratell = COREMAG(1.0+ratell,0.0,0.0);
    ratscl = (pow((g_fitellf + 2.0*g_sigellf - g_fitell),2.0) - 4.0*g_sigell*g_sigell)/(4.0*ratell*ratell);
    ratscl = max(0.25,min(10.0,ratscl));
    for (i = 0; i < NSAMPLE; i++) {
        xmag = 5.0 + 0.1*(float)(i+1);
        xflux = pow(10.0,0.4*xmag);
        ratio = 2.5*log10(1.0+xnoise/xflux);
        lower[i] = g_fit_final - 5.0*sqrt(g_sigma_final*g_sigma_final + ratio*ratio);
        upper[i] = g_fit_final + sqrt(9.0*g_sigma_final*g_sigma_final + 0.0*ratio*ratio);
        g_uppere[i] = g_fitell + 2.0*sqrt(g_sigell*g_sigell + ratscl*ratio*ratio);
        g_uppere[i] = min(0.5,g_uppere[i]);
    }
    g_elllim = min(0.5,max(0.2,g_fitell+2.0*g_sigell));
    fluxlim = 2.5*log10((double)(2.5*sqrt(PI*g_rcore*g_rcore)*g_skynoise));

    /* Ok, final classification loop now... */
    /* CONSTANTS: 0.4, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 9.0, 0.5, -1, -2, 10000s */
    g_nstar = 0;
    g_ngal = 0;
    g_njunk = 0;
    g_ncmp = 0;
    for (i = 0; i < g_nrows; i++) {
        ell = g_ellipticity[i];
        pkht = g_peak_mag[i];
        core = g_core_flux[i];
        iarg = casu_nint(10.0*(core - 5.0));
        iarg = max(1,min(NSAMPLE,iarg)) - 1;
        if (! g_poor) {
            sig1 = max(0.01,(g_fit1 - g_lower1[iarg])/3.0);
            sig2 = max(0.01,(g_fit2 - g_lower2[iarg])/3.0);
        } else {
            sig1 = max(0.01,(g_fit4 - g_lower1[iarg])/3.0);
            sig2 = max(0.01,(g_fit5 - g_lower2[iarg])/3.0);
        }
        sig3 = max(0.01,(g_fit3 - g_lower3[iarg])/3.0);
        denom = (g_wt1/sig1 + g_wt2/sig2 + g_wt3/sig3);
        w1 = (g_wt1/sig1)/denom;
        w2 = (g_wt2/sig2)/denom;
        w3 = (g_wt3/sig3)/denom;
        if (! g_poor) {
            core_small = g_core1_flux[i];
            core_large = g_core3_flux[i];
            statistic = (core - core_small - g_fit1)*w1 + 
                (max(-3.0*sig2,core_large - core - g_fit2))*w2 + 
                (core - pkht - g_fit3)*w3;
        } else {
            core_midd = g_core2_flux[i];
            core_large = g_core4_flux[i];
            statistic = (core_midd - core - g_fit4)*w1 +
                (max(-3.0*sig2,core_large - core - g_fit5))*w2 + 
                (core - pkht - g_fit3)*w3;
        }
        g_cls[i] = -1.0;
        statcut = upper[iarg] + 3.0*g_sigma_final*(exp(max(0.0,core-g_corlim+1.0)) - 1.0);
        if (statistic  >= statcut) 
            g_cls[i] = 1.0;
        else if (statistic <= lower[iarg])
            g_cls[i] = 0.0;

        /* Save distance from the stellar locus */

        sigtot = (g_fit_final - lower[iarg])/5.0;
        g_sig[i] = (statistic - g_fit_final)/sigtot;

        /* Right, now here are lots of overrides for special circumstances */
        /* Too spikey? -> junk */

        if (core - pkht - g_fit3 < -4.0*sig3) 
            g_cls[i] = 0.0;

        /* Elliptical star? -> compact */

        ellbound = max(g_elllim,g_uppere[iarg]);
        if (ell > ellbound && g_cls[i] == -1.0 && core < g_flim && g_sig[i] > -2.0)
            g_cls[i] = -2.0;

        /* Saturated? -> star */

        if (core > g_corlim && statistic >= lower[iarg])
            g_cls[i] = -1.0;

        /* Too elliptical? -> junk */

        if (ell > 0.9 && core < g_corlim)
            g_cls[i] = 0.0;

        /* Too faint? -> junk */

        if (core < fluxlim)
            g_cls[i] = 0.0;

        /* Now count how many you have of each */

        if (g_cls[i] == -1.0)
            g_nstar++;
        else if (g_cls[i] == 1.0)
            g_ngal++;
        else if (g_cls[i] == -2.0)
            g_ncmp++;
        else
            g_njunk++;
    }

    /* Do stats to get the aperture corrections */

    if (g_ncols == NCOLS) {
        classstats_ap67(g_core5_flux,g_core3_flux,&g_fit6,&sigma6);
        classstats_ap67(g_core_flux,g_core6_flux,&g_fit7,&sigma7);
        g_fit6 += g_fit2;
    }
    classstats_ap0(&fit0,&sigma0);
    if (g_ncols == NCOLS) 
        fit0 = max(g_fit6,fit0);
    else
        fit0 = max(g_fit5,fit0);
    g_apcpkht = fit0 + g_fit3; /* pkht */
    switch (g_ncols) {
    case 32:
        g_apcor1 = fit0 + g_fit1;  /* 0.5*core */
        g_apcor = fit0;          /* core */
        g_apcor2 = fit0 - g_fit4;  /* sqrt(2)*core */
        g_apcor3 = fit0 - g_fit2;  /* 2*core */
        g_apcor4 = fit0 - g_fit5;  /* 2*sqrt(2)*core */
        g_apcor5 = 0.0;          /* 4*core */
        break;
    case NCOLS:
        g_apcor1 = fit0 + g_fit1;      /* 0.5*core */
        g_apcor2 = fit0 + g_fit7;      /* 1/sqrt(2) * core */
        g_apcor3 = fit0;             /* core */
        g_apcor4 = fit0 - g_fit4;      /* core * sqrt(2) */
        g_apcor5 = fit0 - g_fit2;      /* 2*core */
        g_apcor6 = fit0 - g_fit5;      /* 2*sqrt(2)*core */
        g_apcor7 = fit0 - g_fit6;      /* 4*core */
        break;
    }

    /* Now do a better job on the saturation */

    ii = 0;
    work = cpl_malloc(g_nrows*sizeof(float));
    for (i = 0; i < g_nrows; i++) {
        ell = g_ellipticity[i];
        core = g_core_flux[i];
        pkht = max(g_thresh,g_peak_height[i]) + g_skylev[i];
        if (((ell < g_elllim && core > g_flim && g_cls[i] == -1 && g_sig[i] >= 5.0 &&
             g_areal[0][i] >= g_pixlim) || pkht >= 0.9*g_avsat) && g_xpos[i] >= g_xmin && 
             g_xpos[i] <= g_xmax && g_ypos[i] >= g_ymin && g_ypos[i] <= g_ymax) {
            work[ii++] = pkht;
        }
    }
    if (ii > 0) {
        medstat(work,ii,&avsatnew,&junk);
        avsatnew = max(10000.0+g_skylevel,avsatnew);
    } else {
        avsatnew = 10000.0 + g_skylevel;
    }
    g_avsat = avsatnew;
    freespace(work);    

    /* Ok, now get rid of some workspace */
        
    freespace(g_lower1);
    freespace(g_lower2);
    freespace(g_lower3);
    freespace(g_upper1);
    freespace(g_upper2);
    freespace(g_upper3);
    freespace(lower);
    freespace(upper);
    freespace(g_uppere);

}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        classstats
    \par Purpose:
        Work out the median difference between two magnitude estimates
    \par Description:
        The difference between the magnitudes of two different apertures
        is calculated for all objects in a list. 
    \par Language:
        C
    \param core1
        The flux array for the first aperture (in magnitudes)
    \param core2
        The flux array for the second aperture (in magnitudes)
    \param small
        Set if the second aperture is smaller than the first
    \param cutlev
        An upper limit to the allowed value of the magnitude difference
    \param medval
        The output median magnitude difference
    \param sigma
        The output sigma in magnitude difference
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void classstats(float *core1, float *core2, intptr_t small, float cutlev,
                       float *medval, float *sigma) {

    intptr_t i,iloop,n;
    float *work,*dc,sigmaold,amult;

    /* Initialise the output values to something stupid */
/* CONSTANTS: 1.e6, 3.0,0.01 */
    *medval = 0.0;
    *sigma = 1.0e6;
    amult = (small == 1 ? -1.0 : 1.0);

    /* Get some workspace */

    work = cpl_malloc(g_nrows*sizeof(float));
    dc = cpl_malloc(g_nrows*sizeof(float));

    /* Work out differences */

    for (i = 0; i < g_nrows; i++)
        dc[i] = amult*(core2[i] - core1[i]);

    /* Do an iteration loop */

    for (iloop = 0; iloop < MAXLOOP; iloop++) {
        sigmaold = *sigma;
        n = 0;

        /* Ok, gather up all the stats */

        for (i = 0; i < g_nrows; i++) {
            
            /* Clipping criteria */

            if (g_ellipticity[i] < g_elllim && core1[i] < g_blim && core1[i] > g_flim &&
                fabs(dc[i] - *medval) < 3.0*(*sigma) &&
                g_xpos[i] >= g_xmin && g_xpos[i] <= g_xmax && g_ypos[i] >= g_ymin &&
                g_ypos[i] <= g_ymax && g_areal[0][i] >= g_pixlim) {
                if (iloop > 0 || (iloop == 0 && dc[i] >= cutlev)) 
                    work[n++] = dc[i];
            }
        }

        /* Sort the work array and find the median and sigma */

        if (n > 0) {
            sort1(work,n);
            if (iloop == 0) {
                anhist(work,n,medval,sigma);
            } else {
                medstat(work,n,medval,sigma);
                *sigma = min(sigmaold,*sigma);
            }
        } else {
            *medval = 0.0; 
            *sigma = 0.01;
        }

        /* Just in case... */

        *sigma = max(*sigma,0.01);
    }

    /* Tidy and exit */

    cpl_free(work);
    cpl_free(dc);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        classstats_el
    \par Purpose:
        Work out the median ellipticity of the sample
    \par Description:
        The median ellipticity of the sample is calculated interatively
    \par Language:
        C
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void classstats_el(void) {
    intptr_t iloop,n,i;
    float *work;

    /* Initialise the mean and sigma to something stupid */
/* CONSTANTS: 1e6,0.5, 2.0 ,0.25,0.05,0.2 */
    g_sigell = 1.0e6;
    g_fitell = 0.0;

    /* Get some workspace */

    work = cpl_malloc(g_nrows*sizeof(float));

    /* Do iteration loop */

    for (iloop = 0; iloop < MAXLOOP; iloop++) {
        n = 0;
        for (i = 0; i < g_nrows; i++) {
            if (g_ellipticity[i] < 0.5 && g_core_flux[i] < g_blim && 
                g_core_flux[i] > g_flim && 
                fabs(g_ellipticity[i] - g_fitell) < 2.0*g_sigell &&
                g_xpos[i] >= g_xmin && g_xpos[i] <= g_xmax && g_ypos[i] >= g_ymin &&
                g_ypos[i] <= g_ymax && g_areal[0][i] >= g_pixlim)
                work[n++] = g_ellipticity[i];
        }
        if (n > 2)
            medstat(work,n,&g_fitell,&g_sigell);
        else {
            g_fitell = 0.25;
            g_sigell = 0.05;
        }
    }
    g_elllim = min(0.5,max(0.2,g_fitell+2.0*g_sigell));

    /* Get out of here */

    freespace(work);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        classstats_pa
    \par Purpose:
        Work out the median position angle of the sample
    \par Description:
        The median position angle of the sample is calculated interatively
    \par Language:
        C
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void classstats_pa(void) {
    intptr_t iloop,n,i;
    float *work;

    /* Initialise the mean and sigma to something stupid */
/* CONSTANTS: 1e6, 2, 0.05 */
    g_sigpa = 1.0e6;
    g_fitpa = 0.0;

    /* Get some workspace */

    work = cpl_malloc(g_nrows*sizeof(float));

    /* Do iteration loop */

    for (iloop = 0; iloop < MAXLOOP; iloop++) {
        n = 0;
        for (i = 0; i < g_nrows; i++) {
            if (g_core_flux[i] < g_blim && g_core_flux[i] > g_flim && 
                fabs(g_pa[i] - g_fitpa) < 2.0*g_sigpa &&
                g_xpos[i] >= g_xmin && g_xpos[i] <= g_xmax && g_ypos[i] >= g_ymin &&
                g_ypos[i] <= g_ymax && g_areal[0][i] >= g_pixlim)
                work[n++] = g_pa[i];
        }
        if (n > 2)
            medstat(work,n,&g_fitpa,&g_sigpa);
        else {
            g_fitpa = 0.0;
            g_sigpa = 0.05;
        }
    }

    /* Get out of here */

    freespace(work);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        classstats_ellf
    \par Purpose:
        Work out the median ellipticity for faint objects
    \par Description:
        The median ellipticity is calculated for objects that are fainter
        than a given flux limit.
    \par Language:
        C
    \param fluxlim
        The flux limit (in magnitudes)
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void classstats_ellf(float fluxlim) {
    intptr_t iloop,n,i;
    float *work;

    /* Initialise the mean and sigma to something stupid */
/* CONSTANTS: 1e6,0.75,0.25,0.05,2.0 */
    g_sigellf = 1.0e6;
    g_fitellf = 0.0;

    /* Get some workspace */

    work = cpl_malloc(g_nrows*sizeof(float));

    /* Do iteration loop */

    for (iloop = 0; iloop < MAXLOOP; iloop++) {
        n = 0;
        for (i = 0; i < g_nrows; i++) {
            if (g_ellipticity[i] < 0.75 && g_core_flux[i] > fluxlim+1.0 && 
                g_core_flux[i] < fluxlim+2.0 &&
                fabs(g_ellipticity[i] - g_fitellf) < 2.0*g_sigellf)
                work[n++] = g_ellipticity[i];
        }
        if (n > 2)
            medstat(work,n,&g_fitellf,&g_sigellf);
        else {
            g_fitellf = 0.25;
            g_sigellf = 0.05;
        }
    }

    /* Get out of here */

    freespace(work);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        classstats_ap0
    \par Purpose:
        Work out the median magnitude difference for large apertures
    \par Description:
        The median magnitude difference for isophotal fluxes and a large
        aperture is calculated.
    \par Language:
        C
    \param medval
        The output median magnitude difference
    \param sigma
        The output sigma of the magnitude difference
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void classstats_ap0(float *medval, float *sigma) {

    intptr_t i,iloop,n;
    float *work,*dc,c2,sigmanew;

    /* Initialise the output values to something stupid */
/* CONSTANTS: 1e6,0.2,2.0,0.5,3.0,5.0,1.48, 0.025,0.01,0.25,s CPL_MATH_STD_MAD*/
    *medval = 0.0;
    *sigma = 1.0e6;
    g_elllim = min(0.5,max(0.2,g_fitell+2.0*g_sigell));

    /* Get some workspace */

    work = cpl_malloc(g_nrows*sizeof(float));
    dc = cpl_malloc(g_nrows*sizeof(float));

    /* Work out differences */

    for (i = 0; i < g_nrows; i++) {
        c2 = max(0.0,max(g_iso_flux[i],g_core5_flux[i]));
        dc[i] = c2 - g_core_flux[i];
    }

    /* Do an iteration loop */

    for (iloop = 0; iloop < MAXLOOP; iloop++) {
        n = 0;

        /* Ok, gather up all the stats */

        for (i = 0; i < g_nrows; i++) {
            
            /* Clipping criteria */

            if (g_ellipticity[i] < g_elllim && g_core_flux[i] < g_blim && 
                g_core_flux[i] > g_flim && 
                fabs(dc[i] - *medval) < 3.0*(*sigma) && 
                g_cls[i] == -1.0 && g_sig[i] < 5.0 &&
                g_xpos[i] >= g_xmin && g_xpos[i] <= g_xmax && g_ypos[i] >= g_ymin &&
                g_ypos[i] <= g_ymax && g_areal[0][i] >= g_pixlim) 
                if (iloop > 0 || (iloop == 0 && dc[i] >= 0.0)) {
                    work[n++] = dc[i];
                }
        }

        /* Sort the work array and find the median and sigma */
        /* TODO AMO: we should use CPL_MATH_STD_MAD instead of 1.48 */
        if (n > 0) {
            sort1(work,n);
            if (iloop == 0) {
                anhist(work,n,medval,sigma);
                *sigma = CPL_MATH_STD_MAD*(*medval - work[(intptr_t)(0.25*(float)(n+3))-1]);
                *sigma = max(0.025,*sigma);
            } else {
                medstat(work,n,medval,&sigmanew);
                *sigma = min(*sigma,sigmanew);
                *sigma = max(0.01,*sigma);
            }
        } else {
            *medval = 0.0;
            *sigma = 0.01;
        }

        /* Just in case... */

        *sigma = max(*sigma,0.01);
    }

    /* Tidy and exit */

    freespace(work);
    freespace(dc);
}


/*---------------------------------------------------------------------------*/
/**
    \par Name:
        classstats_ap67
    \par Purpose:
        TODO
    \par Description:
        TODO
    \par Language:
        C
    \param mag1
        TODO
    \param mag2
        TODO
    \param medval
        TODO
    \param sigma
        The output sigma of the magnitude difference
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */

static void classstats_ap67(float *mag1, float *mag2, float *medval, 
                            float *sigma) {

    intptr_t i,iloop,n;
    float *work,*dc,sigmanew;

    /* Initialise the output values to something stupid */

    *medval = 0.0;
    *sigma = 1.0e6;
    g_elllim = min(0.5,max(0.2,g_fitell+2.0*g_sigell));
    /* CONSTANTS: 1e6,0.2,2.0,3.0,5.0,1.48,0.025,0.25,0.01, CPL_MATH_STD_MAD */
    /* Get some workspace */

    work = cpl_malloc(g_nrows*sizeof(float));
    dc = cpl_malloc(g_nrows*sizeof(float));

    /* Work out differences */

    for (i = 0; i < g_nrows; i++) 
        dc[i] = mag1[i] - mag2[i];

    /* Do an iteration loop */

    for (iloop = 0; iloop < MAXLOOP; iloop++) {
        n = 0;

        /* Ok, gather up all the stats */

        for (i = 0; i < g_nrows; i++) {
            
            /* Clipping criteria */

            if (g_ellipticity[i] < g_elllim && g_core_flux[i] < g_blim && 
                g_core_flux[i] > g_flim && 
                fabs(dc[i] - *medval) < 3.0*(*sigma) && 
                g_cls[i] == -1.0 && g_sig[i] < 5.0 &&
                g_xpos[i] >= g_xmin && g_xpos[i] <= g_xmax && g_ypos[i] >= g_ymin &&
                g_ypos[i] <= g_ymax && g_areal[0][i] >= g_pixlim) {
                if (iloop > 0 || (iloop == 0 && dc[i] >= 0.0)) {
                    work[n++] = dc[i];
                }
            }
        }

        /* Sort the work array and find the median and sigma */

        if (n > 0) {
            sort1(work,n);
            if (iloop == 0) {
                anhist(work,n,medval,sigma);
                *sigma = CPL_MATH_STD_MAD*(*medval - work[(intptr_t)(0.25*(float)(n+3))-1]);
                *sigma = max(0.025,*sigma);
            } else {
                medstat(work,n,medval,&sigmanew);
                *sigma = min(*sigma,sigmanew);
                *sigma = max(0.01,*sigma);
            }
        } else {
            *medval = 0.0;
            *sigma = 0.01;
        }

        /* Just in case... */

        *sigma = max(*sigma,0.01);
    }

    /* Tidy and exit */

    cpl_free(work);
    cpl_free(dc);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        classstats_final
    \par Purpose:
        Routine to define the median classification statistic
    \par Description:
        The result of the previous boundary and statistical routines are
        combined to define the median classification statistic and its 
        sigma. The curve is also investigated to see where saturation occurs
    \par Language:
        C
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void classstats_final(void) {
    intptr_t n,i,iloop,iarg,ii,iend,ncls,kk,k;
    float *work,ell,core,sig1,sig2,sig3,denom,w1,w2,w3,core_small;
    float core_large,*statistic,core_midd,pkht,xcor,cfit,csig;
    float *work1,junk,corlim1,corval1,corlim2,corval2,sigmaold;

    /* Initialise */

    g_sigma_final = 1.0e6;
    g_fit_final = 0.0;
    ncls = 0;
    /* CONSTANTS 1e6, 5, 3, 2, 12.5, 3.0, 10, 0.01, 0.25, 0.5, 10000 */
    /* Get some workspace */

    work = cpl_malloc(g_nrows*sizeof(float));
    work1 = cpl_malloc(g_nrows*sizeof(float));
    statistic = cpl_malloc(g_nrows*sizeof(float));

    /* Calculate the statistic now */

    for (i = 0; i < g_nrows; i++) {
        ell = g_ellipticity[i];
        pkht = g_peak_mag[i];
        core = g_core_flux[i];
        iarg = casu_nint(10.0*(core - 5.0));
        iarg = max(1,min(NSAMPLE,iarg)) - 1;
        if (! g_poor) {
            sig1 = max(0.01,(g_fit1 - g_lower1[iarg])/3.0);
            sig2 = max(0.01,(g_fit2 - g_lower2[iarg])/3.0);
        } else {
            sig1 = max(0.01,(g_fit4 - g_lower1[iarg])/3.0);
            sig2 = max(0.01,(g_fit5 - g_lower2[iarg])/3.0);
        }
        sig3 = max(0.01,(g_fit3 - g_lower3[iarg])/3.0);
        denom = (g_wt1/sig1 + g_wt2/sig2 + g_wt3/sig3);
        w1 = (g_wt1/sig1)/denom;
        w2 = (g_wt2/sig2)/denom;
        w3 = (g_wt3/sig3)/denom;
        if (! g_poor) {
            core_small = g_core1_flux[i];
            core_large = g_core3_flux[i];
            statistic[i] = (core - core_small - g_fit1)*w1 + 
                (core_large - core - g_fit2)*w2 + (core - pkht - g_fit3)*w3;
        } else {
            core_midd = g_core2_flux[i];
            core_large = g_core4_flux[i];
            statistic[i] = (core_midd - core - g_fit4)*w1 +
                (core_large - core - g_fit5)*w2 + (core - pkht - g_fit3)*w3;
        }
    }

    /* Iteration loop.  Use only lower g_ellipticity images and relevant
       peak height range */

    for (iloop = 0; iloop < MAXLOOP; iloop++) {
        sigmaold = g_sigma_final;
        n = 0;
        for (i = 0; i < g_nrows ; i++) {

            ell = g_ellipticity[i];
            core = g_core_flux[i];
            if (ell < g_elllim && core < g_blim && core > g_flim && 
                fabs((double)(statistic[i] - g_fit_final)) < 3.0*g_sigma_final &&
                g_areal[0][i] >= g_pixlim)
                work[n++] = statistic[i];

            /* This information is to be used later to find the curvature of 
               saturated region */

            if (core > g_corlim && iloop == MAXLOOP-2) {
                g_cls[ncls] = statistic[i];
                g_sig[ncls++] = core;
            }
        }

        /* Median defines general fit */

        if (n > 2) {
            sort1(work,n);
            if (iloop == 0 && n > 10) {
                anhist(work,n,&g_fit_final,&g_sigma_final);
            } else {
                medstat(work,n,&g_fit_final,&g_sigma_final);
            }
            g_sigma_final = max(0.01,min(sigmaold,g_sigma_final));
        } else {
            g_fit_final = 0.0;
            g_sigma_final = 0.01;
        }
    }

    /* Now work out the curvature in the saturated region */

    sort2(g_sig,g_cls,ncls);
    ii = 0;
    xcor = 12.5;
    iend = 0;
    i = -1;
    corlim1 = 0.0;
    corlim2 = 0.0;
    corval1 = 0.0;
    corval2 = 0.0;
    while (iend == 0 && i < ncls-1) {
        i++;
        if (g_sig[i] > xcor+0.25 && ii >= 3) {
            medstat(work,ii,&cfit,&csig);
            for (iloop = 0; iloop < 3; iloop++) {
                kk = 0;
                for (k = 0; k < ii; k++) {
                    if (work[k] <= cfit + 3.0*csig)
                        work1[kk++] = work[k];
                }
                medstat(work1,kk,&cfit,&junk);
            }
            if (cfit <= g_fit_final + 3.0*g_sigma_final) {
                corlim1 = xcor;
                corval1 = cfit;
            } else {
                corlim2 = xcor;
                corval2 = cfit;
                iend = 1;
            }
        } else {
            work[ii++] = g_cls[i];
        }
    }

    /* Estimate where core measure and statistic become unreliable */

    if (iend == 1) 
        g_corlim = corlim2 - 0.5*(corval2 - g_fit_final - 3.0*g_sigma_final)/(corval2 - corval1);
    else 
        g_corlim = corlim1;
    g_corlim = max(g_cormin,g_corlim);
    kk = 0;
    for (i = 0; i < g_nrows; i++) {
        core = g_core_flux[i];
        if (core >= g_corlim)
            work[kk++] = g_peak_height[i] + g_skylevel;
    }
    if (kk > 0) {
        medstat(work,kk,&g_avsat,&junk);
        g_avsat = max(10000.0+g_skylevel,g_avsat);
    } else {
        g_avsat = 10000.0 + g_skylevel;
    }

    /* Tidy and exit */
  
    freespace(work);
    freespace(work1);
    freespace(statistic);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        medstat
    \par Purpose:
        Work out the median and sigma of an array
    \par Description:
        The median of an array is calculated by sorting and choosing the 
        middle value. The sigma is estimated by halving the space between
        the first and third quartiles and scaling by the appropriate factor
    \par Language:
        C
    \param array
        The input data array
    \param n
        The number of values in the data array
    \param medval
        The output median value
    \param sigval
        The output sigma value
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void medstat(float *array, intptr_t n, float *medval, float *sigval) {
    intptr_t lev1,lev2,lev3;

    /* Sort the array first, then choose the median.  The sigma is defined
       as half the distance between the two quartile points multiplied by
       the appropriate scaling factor (1.48) */
/* CONSTANTS: 3, 4, 1.48, 0.5. We should use CPL_MATH_STD_MAD instead of 1.48 */
    if (n == 0) {
        *medval = 0.0;
        *sigval = 0.0;
        return;
    }
    sort1(array,n);
    lev1 = (n + 1)/2;
    lev2 = (3*n + 3)/4;
    lev3 = (n + 3)/4;
    *medval = array[lev1-1];
    *sigval = CPL_MATH_STD_MAD*0.5*(array[lev2-1] - array[lev3-1]);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        sort1
    \par Purpose:
        Sort a single array into ascending order
    \par Description:
        A single array is sorted into ascending order
    \par Language:
        C
    \param a
        The input data array
    \param n
        The number of values in the data array
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void sort1(float *a, intptr_t n) {
    intptr_t iii,ii,i,ifin,j;
    float b;
/* CONSTANTS: 4, 2, 3 */
    iii = 4;
    while (iii < n)
        iii *= 2;
    iii = min(n,(3*iii)/4 - 1);

    while (iii > 1) {
        iii /= 2;
        ifin = n - iii;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + iii;
            if (a[i] > a[j]) {
                b = a[j];
                while (1) {
                    a[j] = a[i];
                    j = i;
                    i = i - iii;
                    if (i < 0 || a[i] <= b) 
                        break;
                }
                a[j] = b;
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        sort2
    \par Purpose:
        Sort a single array into ascending order and cosort a second
    \par Description:
        A single array is sorted into ascending order and a second
        array is cosorted
    \par Language:
        C
    \param a1
        The first input data array
    \param a2
        The second input data array
    \param n
        The number of values in the data array
    \returns
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

static void sort2(float *a1, float *a2, intptr_t n) {
    intptr_t iii,ii,i,ifin,j;
    float b1,b2;
    /* CONSTANTS: 4, 2, 3 */
    iii = 4;
    while (iii < n)
        iii *= 2;
    iii = min(n,(3*iii)/4 - 1);

    while (iii > 1) {
        iii /= 2;
        ifin = n - iii;
        for (ii = 0; ii < ifin; ii++) {
            i = ii;
            j = i + iii;
            if (a1[i] > a1[j]) {
                b1 = a1[j];
                b2 = a2[j];
                while (1) {
                    a1[j] = a1[i];
                    a2[j] = a2[i];
                    j = i;
                    i = i - iii;
                    if (i < 0 || a1[i] <= b1) 
                        break;
                }
                a1[j] = b1;
                a2[j] = b2;
            }
        }
    }
}

/**@}*/

/*

$Log: classify.c,v $
Revision 1.3  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.2  2015/08/07 13:06:54  jim
Fixed copyright to ESO

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.4  2015/03/03 10:48:11  jim
Fixed some memory leaks

Revision 1.3  2014/04/09 11:08:21  jim
Get rid of a couple of compiler moans

Revision 1.2  2014/04/09 09:09:51  jim
Detabbed

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
