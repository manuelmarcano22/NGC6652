/*

$Id: imcore.h,v 1.3 2015/09/22 15:09:20 jim Exp $

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

/* Required information */

#ifndef IMCORE_H
#define IMCORE_H

#include "ap.h"
#include "classify.h"
#include "casu_utilfunctions.h"
#include "hdrl.h"
#include <cpl.h>

/* imcore specific macros */

#define MINHISTVAL -1000 /* Minimum value to histogram for background stats */
#define MAXHISTVAL 65535 /* Maximum value to histogram for background stats */
#define MAXHIST (MAXHISTVAL-MINHISTVAL+1) /* maximum size of histogram array */
#define HIST_ELEM(a, i) ((a)[(i) - MINHISTVAL]);

/* Catalogue generation parameters */

#define MINSATURATE 20000 /* Minimum background saturation level */
#define IMNUM 200         /* Maximum number of images to be deblended */
#define NPAR 16           /* Number of parameters in a basic results array */
#define NAREAL 8          /* Number of areal profiles */
#define IDBLIM 10000      /* Maximum number of pixels to use in deblending */
#define INITROWS 2048     /* Allocation size for rows in the output table */

#define STUPID_VALUE -1000 /* Minimum value of a pixel */

/* MFLAG values used for tracking the quality of individual pixels */

#define MF_CLEANPIX     0
#define MF_OBJPIX       1
#define MF_SATURATED    2
#define MF_ZEROCONF     3
#define MF_STUPID_VALUE 4
#define MF_3SIG         5
#define MF_POSSIBLEOBJ  6

/* Flags for types of catalogues to be created */

#define CAT_INTWFC      1   /* Original 32 column catalogue used for INT WFC */
#define CAT_WFCAM       2   /* Catalogue proposed for WFCAM */
#define CAT_BASIC       3   /* Very basic with positions, moments and areals */
#define CAT_OBJMASK     4   /* No table, just an object mask */
#define CAT_VIRCAM      6   /* 80 column VIRCAM table */

/* Tidy-up macros */

#define freespace(_p) if (_p != NULL) {cpl_free(_p); _p = NULL;}
#define closefile(_p) if (_p != NULL) {fclose(_p); _p = NULL;}

/* External Function Prototypes. First, main processing routines */

extern int imcore_conf(casu_fits *infile, casu_fits *conf, intptr_t ipix,
                       float threshold, intptr_t icrowd, float rcore,
                       int bkg_subtr, intptr_t nbsize,
                       int cattype, float filtfwhm, float gain,
                       float saturation,
                       hdrl_imcore_result * res);
extern int imcore_opm(casu_fits *infile, casu_fits *conf, intptr_t ipix, 
                      float threshold, intptr_t nbsize, float filtfwhm, 
                      intptr_t niter);
extern int imcore_background(ap_t *ap, intptr_t nbsize, float nullval, int bkg_subtr,
                             hdrl_imcore_result * res);
extern int imcore_backstats(ap_t *ap, float nullval,
                            float *skymed, float *skysig);
extern void imcore_backest(ap_t *ap, float x, float y, float *skylev, 
                           float *skyrms);
extern void imcore_medsig(intptr_t *hist, intptr_t nh, intptr_t ist, intptr_t itarg, 
                          float *med, float *sig);
extern int imcore_extend(ap_t *, float, float, float, float, float, 
                         float, float, float, float *);
extern void imcore_overlp(ap_t *ap, float [IMNUM][NPAR], intptr_t *, float, 
                          float, float, intptr_t, float);
extern void imcore_phopt(ap_t *, float [IMNUM][NPAR], intptr_t, intptr_t, float [], 
                         float [], float [], intptr_t, float []);
extern void imcore_seeing(ap_t *, intptr_t, float *, float *, float **, 
                          float *, float *);

/* Filter routines */

extern void imcore_bfilt(float **, intptr_t, intptr_t);
extern void imcore_median(float[], intptr_t, intptr_t);

/* Polynomial solution */

extern void imcore_polynm(float [], float [], intptr_t, float [], intptr_t, intptr_t);
extern void imcore_solve (double a[25][25], double b[25], intptr_t m);

/* Routines that generate the catalogues */

extern void imcore_tabinit(ap_t *, intptr_t *, intptr_t *,
                           hdrl_catalogue_options, cpl_table **,
                           hdrl_imcore_result * res);
extern int imcore_tabclose(ap_t *, int);
extern int imcore_do_seeing(ap_t *, int, intptr_t, cpl_table *);
extern int imcore_process_results(ap_t *, int, float, intptr_t *,
				   cpl_table *, hdrl_imcore_result * );
extern void imcore_tabinit_6(intptr_t *, intptr_t *, cpl_table **);
extern int imcore_do_seeing_6(ap_t *, intptr_t, cpl_table *);
extern int imcore_process_results_6(ap_t *, float, intptr_t *, cpl_table *,
				     hdrl_imcore_result *);
extern void imcore_tabinit_gen(intptr_t, const char *[], const char *[], 
                               cpl_type[], cpl_table **);
extern int imcore_do_seeing_gen(ap_t *, const char *, const char *, 
                                char *[], intptr_t, cpl_table *);

/* AP routines */

extern void imcore_apclose(ap_t *);
extern void imcore_apfu(ap_t *);
extern void imcore_apinit(ap_t *);
extern void imcore_apreinit(ap_t *);
extern void imcore_apline(ap_t *, float [], float [], float [], float [], 
                          intptr_t, unsigned char *);
extern void imcore_apclust(ap_t *, intptr_t, plstruct *);
extern void imcore_moments (ap_t *, float []);
extern void imcore_areals(ap_t *, intptr_t [NAREAL]);
extern void imcore_restack(ap_t *, intptr_t);
extern void imcore_terminate(ap_t *, int, float, intptr_t *, cpl_table *,
                             hdrl_imcore_result * res);
extern void imcore_extract_data(ap_t *, intptr_t);

#endif

/*

$Log: imcore.h,v $
Revision 1.3  2015/09/22 15:09:20  jim
Fixed guards and comments

Revision 1.2  2015/08/12 11:16:55  jim
Modified procedure names to protect namespace

Revision 1.1.1.1  2015/06/12 10:44:32  jim
Initial import

Revision 1.3  2015/01/09 11:42:36  jim
Fixed routines to remove globals

Revision 1.2  2014/04/09 09:09:51  jim
Detabbed

Revision 1.1.1.1  2013/08/27 12:07:48  jim
Imported


*/
