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
#include <math.h>
#include "imcore.h"
#include "imcore_radii.h"
#include "util.h"
#include "floatmath.h"
 
#define COL_NUMBER      1
#define COL_FLUXISO     2
#define COL_X           3
#define COL_XERR        4
#define COL_Y           5
#define COL_YERR        6
#define COL_SIGMA       7
#define COL_ELLIPT      8
#define COL_PA          9
#define COL_AREAL1     10
#define COL_AREAL2     11
#define COL_AREAL3     12
#define COL_AREAL4     13
#define COL_AREAL5     14
#define COL_AREAL6     15
#define COL_AREAL7     16
#define COL_AREAL8     17
#define COL_PEAKHEIGHT 18
#define COL_PKHTERR    19
#define COL_APFLUX1    20
#define COL_APFLUX1ERR 21
#define COL_APFLUX2    22
#define COL_APFLUX2ERR 23
#define COL_APFLUX3    24
#define COL_APFLUX3ERR 25
#define COL_APFLUX4    26
#define COL_APFLUX4ERR 27
#define COL_APFLUX5    28
#define COL_APFLUX5ERR 29
#define COL_APFLUX6    30
#define COL_APFLUX6ERR 31
#define COL_APFLUX7    32
#define COL_APFLUX7ERR 33
#define COL_APFLUX8    34
#define COL_APFLUX8ERR 35
#define COL_APFLUX9    36
#define COL_APFLUX9ERR 37
#define COL_APFLUX10    38
#define COL_APFLUX10ERR 39
#define COL_APFLUX11    40
#define COL_APFLUX11ERR 41
#define COL_APFLUX12    42
#define COL_APFLUX12ERR 43
#define COL_APFLUX13    44
#define COL_APFLUX13ERR 45
#define COL_PETRAD      46
#define COL_KRONRAD     47
#define COL_HALFRAD     48
#define COL_PETFLUX     49
#define COL_PETFLUXERR  50
#define COL_KRONFLUX    51
#define COL_KRONFLUXERR 52
#define COL_HALFFLUX    53
#define COL_HALFFLUXERR 54
#define COL_ERRFLAG     55
#define COL_SKYLEVEL    56
#define COL_SKYSIGMA    57
#define COL_AVCONF      58
#define COL_RA          59
#define COL_DEC         60
#define COL_CLASS       61
#define COL_STAT        62
#define COL_FWHM        63

/* Number of columns in the table */
 
#define NCOLS 63

static const char *ttype[NCOLS]={"Sequence_number","Isophotal_flux",
                                 "X_coordinate","X_coordinate_err",
                                 "Y_coordinate","Y_coordinate_err",
                                 "Gaussian_sigma","Ellipticity","Position_angle",
                                 "Areal_1_profile","Areal_2_profile","Areal_3_profile",
                                 "Areal_4_profile","Areal_5_profile","Areal_6_profile",
                                 "Areal_7_profile","Areal_8_profile",
                                 "Peak_height","Peak_height_err",
                                 "Aper_flux_1","Aper_flux_1_err",
                                 "Aper_flux_2","Aper_flux_2_err",
                                 "Aper_flux_3","Aper_flux_3_err",
                                 "Aper_flux_4","Aper_flux_4_err",
                                 "Aper_flux_5","Aper_flux_5_err",
                                 "Aper_flux_6","Aper_flux_6_err",
                                 "Aper_flux_7","Aper_flux_7_err",
                                 "Aper_flux_8","Aper_flux_8_err",
                                 "Aper_flux_9","Aper_flux_9_err",
                                 "Aper_flux_10","Aper_flux_10_err",
                                 "Aper_flux_11","Aper_flux_11_err",
                                 "Aper_flux_12","Aper_flux_12_err",
                                 "Aper_flux_13","Aper_flux_13_err",
                                 "Petr_radius","Kron_radius","Half_radius",
                                 "Petr_flux","Petr_flux_err",
                                 "Kron_flux","Kron_flux_err","Half_flux","Half_flux_err",
                                 "Error_bit_flag","Sky_level","Sky_rms",
                                 "Av_conf",
                                 "RA","DEC","Classification","Statistic",
                                 "FWHM"};

static const char *tunit[NCOLS]={"Number","ADU",
                                 "Pixels","Pixels",
                                 "Pixels","Pixels",
                                 "Pixels","Number","Degrees",
                                 "Pixels","Pixels","Pixels",
                                 "Pixels","Pixels","Pixels",
                                 "Pixels","Pixels",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "ADU","ADU",
                                 "Pixels","Pixels","Pixels",
                                 "ADU","ADU",
                                 "ADU","ADU","ADU","ADU",
                                 "Number","ADU","ADU","Number",
                                 "Degrees","Degrees","Flag","N-sigma",
                                 "Pixels"};

static cpl_type g_tform[NCOLS]={CPL_TYPE_INT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,
                              CPL_TYPE_DOUBLE,CPL_TYPE_DOUBLE,CPL_TYPE_FLOAT,
                              CPL_TYPE_FLOAT,CPL_TYPE_FLOAT};

#define NRADS 13
static float g_rmults[] = {0.5, 1.0/CPL_MATH_SQRT2, 1.0, CPL_MATH_SQRT2,
                           2.0, 2.0*CPL_MATH_SQRT2, 
                           4.0, 5.0, 6.0, 7.0, 8.0, 10, 12.0};
static intptr_t g_nrcore = 2;
static intptr_t g_n2rcore = 4;
static float g_apertures[NRADS];

static intptr_t g_areal_cols[NAREAL] = {COL_AREAL1, COL_AREAL2, COL_AREAL3,
                                   COL_AREAL4, COL_AREAL5, COL_AREAL6,
                                   COL_AREAL7, COL_AREAL8};

/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Initialise type 6 catalogue
  
    \par Name:
        imcore_tabinit_6
    \par Purpose:
        Initialise type 6 catalogue
    \par Description:
        Type 6 catalogue is initialised and the xy columns are identified.
    \par Language:
        C
    \param imcore_xcol
        TODO
    \param imcore_ycol
        TODO
    \param tab
        TODO
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
 
extern void imcore_tabinit_6(intptr_t *imcore_xcol, intptr_t *imcore_ycol, 
                             cpl_table **tab) {

    /* Call the generic routine to open a new output table */

    imcore_tabinit_gen(NCOLS,ttype,tunit,g_tform,tab);

    /* Define RA and Dec columns */

    *imcore_xcol = COL_X;
    *imcore_ycol = COL_Y;
        
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do seeing estimate for type 6 catalogue
  
    \par Name:
        imcore_do_seeing_6
    \par Purpose:
        Do seeing estimate for type 6 catalogue
    \par Description:
        Areal profiles in a type 6 catalogue are analysed and a seeing
        estimate is extracted
    \par Language:
        C
    \param ap
        The current ap structure
    \param nobjects
        TODO
    \param tab
        TODO
    \retval CASU_OK
        If all is well. Currently this is the only return value
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_do_seeing_6(ap_t *ap, intptr_t nobjects, cpl_table *tab) {
    intptr_t retval,i;
    char *areal_colnames[NAREAL];

    /* Sort out the areal profile column names */

    for (i = 0; i < NAREAL; i++) 
        areal_colnames[i] = (char *)ttype[g_areal_cols[i]-1];
 
    /* Just call the generic seeing routine */
 
    retval = imcore_do_seeing_gen(ap,ttype[COL_ELLIPT-1],ttype[COL_PEAKHEIGHT-1],
                                  areal_colnames,nobjects,tab);

    /* Get out of here */
 
    return(retval);
}
        
/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Process results for type 6 catalogue
  
    \par Name:
        imcore_process_results_6
    \par Purpose:
        Create the results for objects in a type 6 catalogue
    \par Description:
        The pixel processing is done for all the parameters wanted for
        a type 6 catalogue
    \par Language:
        C
    \param ap
        The current ap structure
    \param gain
        The header keyword with the gain in e-/ADU
    \param nobjects
        TODO
    \param tab
        TODO
   \param res
        TODO
    \retval CASU_OK
        If all is well. 
    \retval CASU_FATAL
        If peak flux < 0
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_process_results_6(ap_t *ap, float gain, intptr_t *nobjects,
                                    cpl_table *tab, hdrl_imcore_result * res) {
    float momresults[8],ttotal,parmall[IMNUM][NPAR],cflux[NRADS*IMNUM];
    float sxx,syy,srr,sxy,ecc,temp,xx,theta,radeg,ell,iso_flux;
    float apflux1,apflux2,apflux3,apflux4,apflux5,yy,sigma,peak,areal1,apflux6;
    float apflux7,apflux8,apflux9,apflux10,apflux11,apflux12,apflux13,zero;
    float areal2,areal3,areal4,areal5,areal6,areal7,areal8;
    float skylev,skyrms,half_rad[IMNUM],half_flux[IMNUM],kron_flux[IMNUM];
    float petr_flux[IMNUM],kron_rad[IMNUM],petr_rad[IMNUM],badpix[IMNUM];
    float theta_ra,skyvar[IMNUM],cc,dd,sigsq,xxe,yye,peake,kron_fluxe;
    float apflux1e,apflux2e,apflux3e,apflux4e,apflux5e,apflux6e,apflux7e;
    float apflux8e,apflux9e,apflux10e,apflux11e,apflux12e,apflux13e,half_fluxe;
    float petr_fluxe,avconf[IMNUM],rcore_area;
    float fwhm;
    intptr_t iareal[NAREAL],nbit,i,k,nr,mbit,j;
    intptr_t nrows;

    /* Do a basic moments analysis and work out the areal profiles*/

    imcore_moments(ap,momresults);
    if (momresults[0] < 0)
        return(CASU_FATAL);
    imcore_areals(ap,iareal);

    /* See if this object makes the cut in terms of its size.  If not, then
       just return with good status */

    if (iareal[0] < ap->ipnop || momresults[3] < ap->xintmin)
        return(CASU_OK);

    /* Work out the total flux */

    imcore_extend(ap,momresults[3],momresults[1],momresults[2],
                  momresults[4],momresults[5],momresults[6],
                  (float)iareal[0],momresults[7],&ttotal);

    /* Try and deblend the images if it is requested and justified */

    if (iareal[0] >= ap->mulpix && ap->icrowd)
        imcore_overlp(ap,parmall,&nbit,momresults[1],momresults[2],
                      momresults[3],iareal[0],momresults[7]);
    else
        nbit = 1;
    if (nbit == 1) {
        parmall[0][0] = momresults[3];
        parmall[0][1] = momresults[1];
        parmall[0][2] = momresults[2];
        parmall[0][3] = ap->thresh;
        for (i = 4; i < 8; i++)
            parmall[0][i] = momresults[i];
        for (i = 0; i < NAREAL; i++)
            parmall[0][i+8] = (float)iareal[i];
    } else {
        mbit = 0;
        for (i = 0; i < nbit; i++) {
            if (parmall[i][1] > 1.0 && parmall[i][1] < ap->lsiz &&
                parmall[i][2] > 1.0 && parmall[i][2] < ap->csiz) {
                for (j = 0; j < NPAR; j++) 
                    parmall[mbit][j] = parmall[i][j];
                mbit++;
            }
        }
        nbit = mbit;
        if (nbit == 0)
            return(CASU_OK);
    } 

    /* Create a list of g_apertures */

    for (i = 0; i < NRADS; i++) {
        g_apertures[i] = g_rmults[i]*(ap->rcore);
        skyvar[i] = CPL_MATH_PI*g_apertures[i]*g_apertures[i];
    }
    rcore_area = CPL_MATH_PI*pow(ap->rcore,2.0);

    /* Initialise the badpix accumulator */

    for (i = 0; i < nbit; i++) {
        badpix[i] = 0.0;
        avconf[i] = 0.0;
    }

    /* Get the core fluxes in all g_apertures */

    imcore_phopt(ap,parmall,nbit,NRADS,g_apertures,cflux,badpix,g_nrcore,avconf);
    for (i = 0; i < nbit; i++)
        avconf[i] /= rcore_area;

    /* Get half-light radius for all images */

    for (k = 0; k < nbit; k++) {
        half_flux[k] = 0.5*(MAX(parmall[k][0],cflux[k*NRADS+g_n2rcore]));
        half_rad[k] = imcore_halflight(g_apertures,cflux+k*NRADS,half_flux[k],
                                       parmall[k][7],NRADS);
    }

    /* Get Kron radius for all images and get the flux */
/* CONSTANTS: 2.0, 8, 180, 1e-4,0.99,0.5,90, */
    for (k = 0; k < nbit; k++) {
        areal1 = parmall[k][8];
        kron_rad[k] = imcore_kronrad(areal1,g_apertures,cflux+k*NRADS,NRADS);
    }
    imcore_flux(ap,parmall,nbit,kron_rad,kron_flux,NRADS,g_apertures,cflux);

    /* Get Petrosian radius for all images and get the flux */

    for (k = 0; k < nbit; k++) {
        areal1 = parmall[k][8];
        petr_rad[k] = imcore_petrad(areal1,g_apertures,cflux+k*NRADS,NRADS);
    }
    imcore_flux(ap,parmall,nbit,petr_rad,petr_flux,NRADS,g_apertures,cflux); 

    /* Massage the results and write them to the fits table */

    sigsq = powf(ap->sigma,2.0);
    radeg = 180.0/CPL_MATH_PI;
    for (k = 0; k < nbit; k++) {
        sxx = parmall[k][4];
        sxy = parmall[k][5];
        syy = parmall[k][6];
        if(sxy > 0.0)
          sxy = MAX(1.0e-4,MIN(sxy,sqrtf(sxx*syy)));
        else
          sxy = MIN(-1.0e-4,MAX(sxy,-sqrtf(sxx*syy)));
        srr = MAX(0.5,sxx+syy);
        ecc = sqrtf((syy-sxx)*(syy-sxx)+4.0*sxy*sxy)/srr;
        temp = MAX((1.0-ecc)/(1.0+ecc),0.0);
        ell = 1.0 - sqrtf(temp);
        ell = MIN(0.99,MAX(0.0,ell));
        xx = 0.5*(1.0+ecc)*srr-sxx;
        if(xx == 0.0)
            theta = 0.0;
        else
            theta = 90.0 - radeg*atanf(sxy/xx);
        theta_ra = theta/radeg;
        cc = (1.0 + ecc)*pow(cos(theta_ra),2.0) + (1.0 - ecc)*pow(sin(theta_ra),2.0);
        dd = (1.0 + ecc)*pow(sin(theta_ra),2.0) + (1.0 - ecc)*pow(cos(theta_ra),2.0);

        /* Create a list of values */

        nrows = cpl_table_get_nrow(tab);
        (*nobjects)++;
        if (*nobjects > nrows) 
            (void)cpl_table_set_size(tab,nrows+INITROWS);
        nr = *nobjects - 1;
        iso_flux = parmall[k][0];
        apflux1 = cflux[k*NRADS + 0];
        apflux2 = cflux[k*NRADS + 1];
        apflux3 = cflux[k*NRADS + 2];
        apflux4 = cflux[k*NRADS + 3];
        apflux5 = cflux[k*NRADS + 4];
        apflux6 = cflux[k*NRADS + 5];
        apflux7 = cflux[k*NRADS + 6];
        apflux8 = cflux[k*NRADS + 7];
        apflux9 = cflux[k*NRADS + 8];
        apflux10 = cflux[k*NRADS + 9];
        apflux11 = cflux[k*NRADS + 10];
        apflux12 = cflux[k*NRADS + 11];
        apflux13 = cflux[k*NRADS + 12];
        peak = parmall[k][7];
        xx = parmall[k][1];
        xxe = sqrtf((2.0*sigsq/(CPL_MATH_PI*peak*peak)) + cc/(2.0*CPL_MATH_PI*gain*peak) + 
                    0.0001);
        yy = parmall[k][2];
        yye = sqrtf((2.0*sigsq/(CPL_MATH_PI*peak*peak)) + dd/(2.0*CPL_MATH_PI*gain*peak) + 
                    0.0001);
        sigma = sqrt(srr);
        fwhm = sqrt( sigma * sigma / 2.) * CPL_MATH_FWHM_SIG;
        /* heuristic correction of moment based fwhm obtained via simulated 2d gaussians,
         * with gaussians 4.3 corrects slighltly better but it is not very
         * significant and 4.0 is the same factor sextractor uses */
        fwhm -= 1. / (4. * fwhm);
        areal1 = parmall[k][8];
        areal2 = parmall[k][9];
        areal3 = parmall[k][10];
        areal4 = parmall[k][11];
        areal5 = parmall[k][12];
        areal6 = parmall[k][13];
        areal7 = parmall[k][14];
        if (nbit > 1 && k == 0) 
            areal8 = 0.0;
        else
            areal8 = parmall[k][15];
        imcore_backest(ap,xx,yy,&skylev,&skyrms);
        peake = sqrtf(peak/gain + sigsq + skyrms*skyrms);
        kron_fluxe = sqrt(kron_flux[k]/gain + 
                          (sigsq + skyrms*skyrms)*CPL_MATH_PI*powf(kron_rad[k],2.0));
        half_fluxe = sqrt(MAX(half_flux[k],0.0)/gain + 
                         (sigsq + skyrms*skyrms)*CPL_MATH_PI*powf(half_rad[k],2.0));
        petr_fluxe = sqrt(petr_flux[k]/gain + 
                          (sigsq + skyrms*skyrms)*CPL_MATH_PI*powf(petr_rad[k],2.0));
        apflux1e = sqrt(MAX(0.0,apflux1/gain) + skyvar[0]*(sigsq + skyrms*skyrms));
        apflux2e = sqrt(MAX(0.0,apflux2/gain) + skyvar[1]*(sigsq + skyrms*skyrms));
        apflux3e = sqrt(MAX(0.0,apflux3/gain) + skyvar[2]*(sigsq + skyrms*skyrms));
        apflux4e = sqrt(MAX(0.0,apflux4/gain) + skyvar[3]*(sigsq + skyrms*skyrms));
        apflux5e = sqrt(MAX(0.0,apflux5/gain) + skyvar[4]*(sigsq + skyrms*skyrms));
        apflux6e = sqrt(MAX(0.0,apflux6/gain) + skyvar[5]*(sigsq + skyrms*skyrms));
        apflux7e = sqrt(MAX(0.0,apflux7/gain) + skyvar[6]*(sigsq + skyrms*skyrms));
        apflux8e = sqrt(MAX(0.0,apflux8/gain) + skyvar[7]*(sigsq + skyrms*skyrms));
        apflux9e = sqrt(MAX(0.0,apflux9/gain) + skyvar[8]*(sigsq + skyrms*skyrms));
        apflux10e = sqrt(MAX(0.0,apflux10/gain) + skyvar[9]*(sigsq + skyrms*skyrms));
        apflux11e = sqrt(MAX(0.0,apflux11/gain) + skyvar[10]*(sigsq + skyrms*skyrms));
        apflux12e = sqrt(MAX(0.0,apflux12/gain) + skyvar[11]*(sigsq + skyrms*skyrms));
        apflux13e = sqrt(MAX(0.0,apflux13/gain) + skyvar[12]*(sigsq + skyrms*skyrms));

        /* Store away the results for this object */

        zero = 0.0;
        cpl_table_set_int(tab,ttype[COL_NUMBER-1],nr,*nobjects);
        cpl_table_set_float(tab,ttype[COL_FLUXISO-1],nr,iso_flux);
        cpl_table_set_float(tab,ttype[COL_X-1],nr,xx);
        cpl_table_set_float(tab,ttype[COL_XERR-1],nr,xxe);
        cpl_table_set_float(tab,ttype[COL_Y-1],nr,yy);
        cpl_table_set_float(tab,ttype[COL_YERR-1],nr,yye);
        cpl_table_set_float(tab,ttype[COL_SIGMA-1],nr,sigma);
        cpl_table_set_float(tab,ttype[COL_ELLIPT-1],nr,ell);
        cpl_table_set_float(tab,ttype[COL_PA-1],nr,theta);
        cpl_table_set_float(tab,ttype[COL_AREAL1-1],nr,areal1);
        cpl_table_set_float(tab,ttype[COL_AREAL2-1],nr,areal2);
        cpl_table_set_float(tab,ttype[COL_AREAL3-1],nr,areal3);
        cpl_table_set_float(tab,ttype[COL_AREAL4-1],nr,areal4);
        cpl_table_set_float(tab,ttype[COL_AREAL5-1],nr,areal5);
        cpl_table_set_float(tab,ttype[COL_AREAL6-1],nr,areal6);
        cpl_table_set_float(tab,ttype[COL_AREAL7-1],nr,areal7);
        cpl_table_set_float(tab,ttype[COL_AREAL8-1],nr,areal8);
        cpl_table_set_float(tab,ttype[COL_PEAKHEIGHT-1],nr,peak);
        cpl_table_set_float(tab,ttype[COL_PKHTERR-1],nr,peake);
        cpl_table_set_float(tab,ttype[COL_APFLUX1-1],nr,apflux1);
        cpl_table_set_float(tab,ttype[COL_APFLUX1ERR-1],nr,apflux1e);
        cpl_table_set_float(tab,ttype[COL_APFLUX2-1],nr,apflux2);
        cpl_table_set_float(tab,ttype[COL_APFLUX2ERR-1],nr,apflux2e);
        cpl_table_set_float(tab,ttype[COL_APFLUX3-1],nr,apflux3);
        cpl_table_set_float(tab,ttype[COL_APFLUX3ERR-1],nr,apflux3e);
        cpl_table_set_float(tab,ttype[COL_APFLUX4-1],nr,apflux4);
        cpl_table_set_float(tab,ttype[COL_APFLUX4ERR-1],nr,apflux4e);
        cpl_table_set_float(tab,ttype[COL_APFLUX5-1],nr,apflux5);
        cpl_table_set_float(tab,ttype[COL_APFLUX5ERR-1],nr,apflux5e);
        cpl_table_set_float(tab,ttype[COL_APFLUX6-1],nr,apflux6);
        cpl_table_set_float(tab,ttype[COL_APFLUX6ERR-1],nr,apflux6e);
        cpl_table_set_float(tab,ttype[COL_APFLUX7-1],nr,apflux7);
        cpl_table_set_float(tab,ttype[COL_APFLUX7ERR-1],nr,apflux7e);
        cpl_table_set_float(tab,ttype[COL_APFLUX8-1],nr,apflux8);
        cpl_table_set_float(tab,ttype[COL_APFLUX8ERR-1],nr,apflux8e);
        cpl_table_set_float(tab,ttype[COL_APFLUX9-1],nr,apflux9);
        cpl_table_set_float(tab,ttype[COL_APFLUX9ERR-1],nr,apflux9e);
        cpl_table_set_float(tab,ttype[COL_APFLUX10-1],nr,apflux10);
        cpl_table_set_float(tab,ttype[COL_APFLUX10ERR-1],nr,apflux10e);
        cpl_table_set_float(tab,ttype[COL_APFLUX11-1],nr,apflux11);
        cpl_table_set_float(tab,ttype[COL_APFLUX11ERR-1],nr,apflux11e);
        cpl_table_set_float(tab,ttype[COL_APFLUX12-1],nr,apflux12);
        cpl_table_set_float(tab,ttype[COL_APFLUX12ERR-1],nr,apflux12e);
        cpl_table_set_float(tab,ttype[COL_APFLUX13-1],nr,apflux13);
        cpl_table_set_float(tab,ttype[COL_APFLUX13ERR-1],nr,apflux13e);
        cpl_table_set_float(tab,ttype[COL_PETRAD-1],nr,0.5*petr_rad[k]);
        cpl_table_set_float(tab,ttype[COL_KRONRAD-1],nr,0.5*kron_rad[k]);
        cpl_table_set_float(tab,ttype[COL_HALFRAD-1],nr,half_rad[k]);
        cpl_table_set_float(tab,ttype[COL_PETFLUX-1],nr,petr_flux[k]);
        cpl_table_set_float(tab,ttype[COL_PETFLUXERR-1],nr,petr_fluxe);
        cpl_table_set_float(tab,ttype[COL_KRONFLUX-1],nr,kron_flux[k]);
        cpl_table_set_float(tab,ttype[COL_KRONFLUXERR-1],nr,kron_fluxe);
        cpl_table_set_float(tab,ttype[COL_HALFFLUX-1],nr,half_flux[k]);
        cpl_table_set_float(tab,ttype[COL_HALFFLUXERR-1],nr,half_fluxe);
        cpl_table_set_float(tab,ttype[COL_ERRFLAG-1],nr,badpix[k]);
        cpl_table_set_float(tab,ttype[COL_SKYLEVEL-1],nr,skylev);
        cpl_table_set_float(tab,ttype[COL_SKYSIGMA-1],nr,skyrms);
        cpl_table_set_float(tab,ttype[COL_AVCONF-1],nr,avconf[k]);
        cpl_table_set_float(tab,ttype[COL_FWHM-1],nr, fwhm);
        /* Store away some dummy values to avoid problems later on */
        cpl_table_set_double(tab,ttype[COL_RA-1],nr,zero);
        cpl_table_set_double(tab,ttype[COL_DEC-1],nr,zero);
        cpl_table_set_float(tab,ttype[COL_CLASS-1],nr,100.0);
        cpl_table_set_float(tab,ttype[COL_STAT-1],nr,zero);
    }

    /*Now that everything is okay - fill in the segmenataion map*/
    cpl_msg_debug(cpl_func, "nobjects: %zd",*nobjects);
    if (res->segmentation_map) {
        for (intptr_t index = 0; index < ap->npl_pix; index++) {
            cpl_image_set(res->segmentation_map, ap->plarray[index].x,
                          ap->plarray[index].y,*nobjects);
        }
    }

    /* Get outta here */

    return(CASU_OK);
}

/**@}*/
