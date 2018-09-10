/* $Id: vimos_calib_mult.c,v 1.19 2013-10-22 16:57:09 cgarcia Exp $
 *
 * This file is part of the FORS Data Reduction Pipeline
 * Copyright (C) 2010 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-10-22 16:57:09 $
 * $Revision: 1.19 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vimos_calib_mult.h>

#include <memory>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <cpl.h>
#include <moses.h>
#include <vimos_dfs.h>

#include "vimos_grism.h"
#include "vimos_lines.h"
#include "vimos_detected_slits.h"
#include "vimos_calibrated_slits.h"
#include "vimos_overscan.h"
#include "vimos_detmodel.h"
#include "vimos_flat_normalise.h"
#include "flat_combine.h"

std::auto_ptr<mosca::image> vimos_calmul_flat_mos_create_master_flat
(vimos::calibrated_slits& calibrated_slits, 
 const mosca::wavelength_calibration& wave_cal,
 const mosca::grism_config& grism_cfg,
 cpl_image *master_bias, cpl_image * master_bias_err,
 cpl_propertylist * master_bias_header,
 cpl_frameset * frameset,
 const char * flat_tag);

int vimos_calmul_flats_save
(std::auto_ptr<mosca::image>& master_flat_d, 
 std::auto_ptr<mosca::image>& norm_flat,
 cpl_frameset * frameset, const char * flat_tag, 
 const char * master_screen_flat_tag, const char * master_norm_flat_tag, 
 cpl_parameterlist * parlist, cpl_propertylist * qc_list, int i);

cpl_propertylist * vimos_calmul_flat_qc(mosca::image& master_flat, 
                                      cpl_table * slits,
                                      int nx, int ny, int nflats,
                                      mosca::grism_config& grism_cfg,
                                      double alltime, double xwidth,
                                      double ywidth, double gain, 
                                      double focu_scale);

#define vimos_calmul_exit(message)            \
{                                             \
if ((const char *)message != NULL) cpl_msg_error(recipe, message);  \
cpl_free(pipefile);                           \
cpl_free(fiterror);                           \
cpl_free(fitlines);                           \
cpl_free(wcoeff);                             \
if(globals != NULL)                           \
  for (i = 0; i < nglobal; i++) cpl_table_delete(globals[i]); \
cpl_free(globals);                            \
cpl_image_delete(master_bias);                \
cpl_image_delete(coordinate);                 \
cpl_image_delete(checkwave);                  \
cpl_image_delete(flat);                       \
cpl_image_delete(master_flat);                \
cpl_image_delete(rainbow);                    \
cpl_image_delete(rectified);                  \
cpl_image_delete(residual);                   \
cpl_image_delete(smo_flat);                   \
cpl_image_delete(spatial);                    \
cpl_image_delete(spectra);                    \
cpl_image_delete(wavemap);                    \
cpl_image_delete(delta);                      \
cpl_mask_delete(refmask);                     \
cpl_propertylist_delete(header);              \
cpl_propertylist_delete(save_header);         \
cpl_propertylist_delete(qclist);              \
cpl_table_delete(grism_table);                \
cpl_table_delete(idscoeff);                   \
cpl_table_delete(restable);                   \
cpl_table_delete(maskslits);                  \
cpl_table_delete(traces);                     \
cpl_table_delete(polytraces);                 \
cpl_table_delete(slits);                      \
cpl_table_delete(subslits);                   \
cpl_table_delete(restab);                     \
cpl_table_delete(global);                     \
cpl_vector_delete(lines);                     \
cpl_msg_indent_less();                        \
return -1;                                    \
}

#define vimos_calmul_exit_memcheck(message)     \
{                                               \
if ((const char *)message != NULL) cpl_msg_info(recipe, message); \
printf("free instrume (%p)\n", instrume);       \
cpl_free(instrume);                             \
printf("free pipefile (%p)\n", pipefile);       \
cpl_free(pipefile);                             \
printf("free fiterror (%p)\n", fiterror);       \
cpl_free(fiterror);                             \
printf("free fitlines (%p)\n", fitlines);       \
cpl_free(fitlines);                             \
printf("free bias (%p)\n", bias);               \
cpl_image_delete(bias);                         \
printf("free master_bias (%p)\n", master_bias); \
cpl_image_delete(master_bias);                  \
printf("free coordinate (%p)\n", coordinate);   \
cpl_image_delete(coordinate);                   \
printf("free checkwave (%p)\n", checkwave);     \
cpl_image_delete(checkwave);                    \
printf("free flat (%p)\n", flat);               \
cpl_image_delete(flat);                         \
printf("free master_flat (%p)\n", master_flat); \
cpl_image_delete(master_flat);                  \
printf("free rainbow (%p)\n", rainbow);         \
cpl_image_delete(rainbow);                      \
printf("free rectified (%p)\n", rectified);     \
cpl_image_delete(rectified);                    \
printf("free residual (%p)\n", residual);       \
cpl_image_delete(residual);                     \
printf("free smo_flat (%p)\n", smo_flat);       \
cpl_image_delete(smo_flat);                     \
printf("free spatial (%p)\n", spatial);         \
cpl_image_delete(spatial);                      \
printf("free spectra (%p)\n", spectra);         \
cpl_image_delete(spectra);                      \
printf("free wavemap (%p)\n", wavemap);         \
cpl_image_delete(wavemap);                      \
printf("free delta (%p)\n", delta);             \
cpl_image_delete(delta);                        \
printf("free refmask (%p)\n", refmask);         \
cpl_mask_delete(refmask);                       \
printf("free header (%p)\n", header);           \
cpl_propertylist_delete(header);                \
printf("free save_header (%p)\n", save_header); \
cpl_propertylist_delete(save_header);           \
printf("free qclist (%p)\n", qclist);           \
cpl_propertylist_delete(qclist);                \
printf("free grism_table (%p)\n", grism_table); \
cpl_table_delete(grism_table);                  \
printf("free idscoeff (%p)\n", idscoeff);       \
cpl_table_delete(idscoeff);                     \
printf("free restable (%p)\n", restable);       \
cpl_table_delete(restable);                     \
printf("free maskslits (%p)\n", maskslits);     \
cpl_table_delete(maskslits);                    \
printf("free traces (%p)\n", traces);           \
cpl_table_delete(traces);                       \
printf("free polytraces (%p)\n", polytraces);   \
cpl_table_delete(polytraces);                   \
printf("free slits (%p)\n", slits);             \
cpl_table_delete(slits);                        \
printf("free subslits (%p)\n", subslits);       \
cpl_table_delete(subslits);                     \
printf("free restab (%p)\n", restab);           \
cpl_table_delete(restab);                       \
printf("free global (%p)\n", global);           \
cpl_table_delete(global);                       \
printf("free lines (%p)\n", lines);             \
cpl_vector_delete(lines);                       \
cpl_msg_indent_less();                          \
return 0;                                       \
}

/**
 * @addtogroup vimos_calib_mult
 */

/**@{*/

/**
 * @brief    Interpret the command line options and execute the data processing
 *
 * @param    parlist     The parameters list
 * @param    frameset    The set-of-frames
 *
 * @return   0 if everything is ok
 */

int vimos_calib_mult(cpl_frameset *frameset, cpl_parameterlist *parlist,
                     cpl_table *allmaskslits)
{
    const char *recipe = "vmmoscalib";

    /*
     * Input parameters
     */

    double      dispersion;
    double      peakdetection;
    int         wdegree;
    int         wradius;
    double      wreject;
    int         wmodelss;
    int         wmodemos;
    const char *ignore_lines;
    const char *used_linesets;
    int         cdegree;
    int         cmode;
    double      startwavelength;
    double      endwavelength;
    double      reference;
    int         slit_ident;
    int         spa_polyorder;
    int         disp_nknots;
    int         sradius;
    int         dradius;
    float       fit_threshold;
    float       disp_tolerance;
    float       ratio_tolerance;

    /*
     * CPL objects
     */

    cpl_image        *master_bias = NULL;
    cpl_image        *master_bias_err = NULL;
    cpl_image        *flat        = NULL;
    cpl_image        *master_flat = NULL;
    cpl_image        *smo_flat    = NULL;
    cpl_image        *spectra     = NULL;
    cpl_image        *wavemap     = NULL;
    cpl_image        *delta       = NULL;
    cpl_image        *residual    = NULL;
    cpl_image        *checkwave   = NULL;
    cpl_image        *rectified   = NULL;
    cpl_image        *dummy       = NULL;
    cpl_image        *refimage    = NULL;
    cpl_image        *coordinate  = NULL;
    cpl_image        *rainbow     = NULL;
    cpl_image        *spatial     = NULL;

    cpl_mask         *refmask     = NULL;

    cpl_table        *grism_table = NULL;
    cpl_table        *idscoeff    = NULL;
    cpl_table        *restable    = NULL;
    cpl_table        *slits       = NULL;
    cpl_table        *subslits    = NULL;
    cpl_table        *positions   = NULL;
    cpl_table        *maskslits   = NULL;
    cpl_table        *traces      = NULL;
    cpl_table        *polytraces  = NULL;
    cpl_table        *restab      = NULL;
    cpl_table        *global      = NULL;
    cpl_table       **globals     = NULL;
    int               nglobal     = 0;

    cpl_vector       *lines       = NULL;

    cpl_propertylist *header      = NULL;
    cpl_propertylist *save_header = NULL;
    cpl_propertylist *qclist      = NULL;

    /*
     * Auxiliary variables
     */

    /* cpl_table  *idscoeff_lss = NULL; */
    char        version[80];
    const char *arc_tag;
    const char *flat_tag;
    const char *master_screen_flat_tag;
    const char *master_norm_flat_tag;
    const char *reduced_lamp_tag;
    const char *disp_residuals_tag;
    const char *disp_coeff_tag;
    const char *wavelength_map_tag;
    const char *spectra_detection_tag;
    const char *spectral_resolution_tag;
    const char *slit_map_tag;
    const char *curv_traces_tag;
    const char *curv_coeff_tag;
    const char *spatial_map_tag;
    const char *slit_location_tag;
    const char *global_distortion_tag = "GLOBAL_DISTORTION_TABLE";
    const char *disp_residuals_table_tag;
    const char *detected_lines_tag;
    const char *delta_image_tag;
    const char *flat_sed_tag = NULL;
    const char *key_gris_name;
    const char *key_gris_id;
    const char *key_filt_name;
    const char *key_filt_id;
    const char *key_mask_id;
    const char *keyname;
    int         quadrant;
    int         mos;
    int         nslits;
    int         ngroups;
    /* double     *xpos; 
    double      mxpos; */
    double      mean_rms;
    double      alltime, arctime;
    int         nflats;
    int         nlines;
    double     *line;
    double     *fiterror = NULL;
    int        *fitlines = NULL;
    int         nx, ny;
    double      gain;
    int         ccd_xsize, ccd_ysize;
    int         cslit, cslit_id;
    double      xwidth, ywidth = 0;
    const int   nchunks = 5; // This is ~CCD y size / mos_region_size (moses.c)
    int         failures;
    int         rotate = 1;
    int         rotate_back = -1;
    int         flat_was_not_saved = 1;
    int         first = 1;
    int         i, j;

    const char *instrume = NULL;
    char       *pipefile = NULL;
    char       *grism;

    double     *wcoeff   = NULL;
    double      scale = 0;
    double      slit_width = 0;
    float       lambdaHe = 0;
    float       lambdaNe = 0;
    float       lambdaAr = 0;
    float       lambdaRed = 0;
    float       lambdaYel = 0;
    float       lambdaBlu = 0;
    double      ar_flux     = 0.0;
    double      ar_flux_err = 0.0;
    double      ne_flux     = 0.0;
    double      ne_flux_err = 0.0;
    double      he_flux     = 0.0;
    double      he_flux_err = 0.0;
    double      r_resol     = 0.0;
    double      r_resol_err = 0.0;
    double      y_resol     = 0.0;
    double      y_resol_err = 0.0;
    double      b_resol     = 0.0;
    double      b_resol_err = 0.0;
    double      qc_mean_rms = 0.0;

    snprintf(version, 80, "%s-%s", PACKAGE, PACKAGE_VERSION);

    cpl_msg_set_indentation(2);


    /* 
     * Get configuration parameters
     */

    cpl_msg_info(recipe, "Recipe %s configuration parameters:", recipe);
    cpl_msg_indent_more();

    if (cpl_frameset_count_tags(frameset, "CONFIG_TABLE") > 1)
        vimos_calmul_exit("Too many in input: CONFIG_TABLE");

    grism_table = dfs_load_table(frameset, "CONFIG_TABLE", 1);

    dispersion = dfs_get_parameter_double(parlist, 
                    "vimos.vmmoscalib.dispersion", grism_table);

    if (dispersion <= 0.0)
        vimos_calmul_exit("Invalid spectral dispersion value");

    peakdetection = dfs_get_parameter_double(parlist, 
                    "vimos.vmmoscalib.peakdetection", grism_table);
    if (peakdetection <= 0.0)
        vimos_calmul_exit("Invalid peak detection level");

    wdegree = dfs_get_parameter_int(parlist, 
                    "vimos.vmmoscalib.wdegree", grism_table);

    if (wdegree < 1)
        vimos_calmul_exit("Invalid polynomial degree");

    if (wdegree > 5)
        vimos_calmul_exit("Max allowed polynomial degree is 5");

    wradius = dfs_get_parameter_int(parlist, "vimos.vmmoscalib.wradius", NULL);

    if (wradius < 0)
        vimos_calmul_exit("Invalid search radius");

    wreject = dfs_get_parameter_double(parlist, 
                                       "vimos.vmmoscalib.wreject", NULL);

    if (wreject <= 0.0)
        vimos_calmul_exit("Invalid rejection threshold");

    wmodelss = dfs_get_parameter_int(parlist, 
                                     "vimos.vmmoscalib.wmodelss", NULL);

    if (wmodelss < 0 || wmodelss > 2)
        vimos_calmul_exit("Invalid wavelength solution interpolation mode");

    wmodemos = dfs_get_parameter_int(parlist, 
                                     "vimos.vmmoscalib.wmodemos", NULL);

    if (wmodemos < 0 || wmodemos > 2)
        vimos_calmul_exit("Invalid wavelength solution interpolation mode");

    ignore_lines= dfs_get_parameter_string(parlist, 
                   "vimos.vmmoscalib.ignore_lines", NULL);

    used_linesets= dfs_get_parameter_string(parlist, 
                   "vimos.vmmoscalib.used_linesets", NULL);

    cdegree = dfs_get_parameter_int(parlist, 
                    "vimos.vmmoscalib.cdegree", grism_table);

    if (cdegree < 1)
        vimos_calmul_exit("Invalid polynomial degree");

    if (cdegree > 5)
        vimos_calmul_exit("Max allowed polynomial degree is 5");

    cmode = dfs_get_parameter_int(parlist, "vimos.vmmoscalib.cmode", NULL);

    if (cmode < 0 || cmode > 2)
        vimos_calmul_exit("Invalid curvature solution interpolation mode");

    startwavelength = dfs_get_parameter_double(parlist, 
                    "vimos.vmmoscalib.startwavelength", grism_table);
    if (startwavelength > 1.0)
        if (startwavelength < 3000.0 || startwavelength > 13000.0)
            vimos_calmul_exit("Invalid wavelength");

    endwavelength = dfs_get_parameter_double(parlist, 
                    "vimos.vmmoscalib.endwavelength", grism_table);
    if (endwavelength > 1.0) {
        if (endwavelength < 3000.0 || endwavelength > 13000.0)
            vimos_calmul_exit("Invalid wavelength");
        if (startwavelength < 1.0)
            vimos_calmul_exit("Invalid wavelength interval");
    }

    if (startwavelength > 1.0)
        if (endwavelength - startwavelength <= 0.0)
            vimos_calmul_exit("Invalid wavelength interval");

    reference = dfs_get_parameter_double(parlist,
                "vimos.vmmoscalib.reference", grism_table);

    if (reference < startwavelength || reference > endwavelength)
        vimos_calmul_exit("Invalid reference wavelength");

    slit_ident = dfs_get_parameter_bool(parlist, 
                    "vimos.vmmoscalib.slit_ident", NULL);

    if (slit_ident == 0) {
        cpl_msg_warning(recipe, "Slit identification is mandatory in case of "
                        "spectral multiplexing. Setting --slit_ident=true...");
        slit_ident = 1;
    }

    spa_polyorder = dfs_get_parameter_int(parlist, "vimos.vmmoscalib.s_degree", NULL);

    disp_nknots = dfs_get_parameter_int(parlist, "vimos.vmmoscalib.d_nknots", NULL);

    sradius = dfs_get_parameter_int(parlist, "vimos.vmmoscalib.sradius", NULL);

    dradius = dfs_get_parameter_int(parlist, "vimos.vmmoscalib.dradius", NULL);

    fit_threshold = dfs_get_parameter_double(parlist, 
            "vimos.vmmoscalib.fit_threshold", NULL);

    ratio_tolerance = dfs_get_parameter_double(parlist, 
            "vimos.vmmoscalib.line_ident_tol", NULL);
    
    disp_tolerance = 0.2; 
            
    if (cpl_error_get_code())
        vimos_calmul_exit("Failure getting the configuration parameters");


    /* 
     * Check input set-of-frames
     */

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Check input set-of-frames:");
    cpl_msg_indent_more();

    {
        cpl_frameset *subframeset = cpl_frameset_duplicate(frameset);
        cpl_frameset_erase(subframeset, "CONFIG_TABLE");
        cpl_frameset_erase(subframeset, "LINE_CATALOG");

        if (!dfs_equal_keyword(subframeset, "ESO OCS CON QUAD"))
            vimos_calmul_exit("Input frames are not from the same quadrant");

        cpl_frameset_delete(subframeset);
    }

    mos = cpl_frameset_count_tags(frameset, "MOS_ARC_SPECTRUM");

    if (mos == 0)
        vimos_calmul_exit("Missing input arc lamp frame");

    if (mos > 1)
        vimos_calmul_exit("Just one input arc lamp frame is allowed"); 

    arc_tag                  = "MOS_ARC_SPECTRUM";
    flat_tag                 = "MOS_SCREEN_FLAT";
    master_screen_flat_tag   = "MOS_COMBINED_SCREEN_FLAT";
    master_norm_flat_tag     = "MOS_MASTER_SCREEN_FLAT";
    reduced_lamp_tag         = "MOS_ARC_SPECTRUM_EXTRACTED";
    disp_residuals_tag       = "MOS_DISP_RESIDUALS";
    disp_coeff_tag           = "MOS_DISP_COEFF";
    wavelength_map_tag       = "MOS_WAVELENGTH_MAP";
    spectra_detection_tag    = "MOS_SPECTRA_DETECTION";
    spectral_resolution_tag  = "MOS_SPECTRAL_RESOLUTION";
    slit_map_tag             = "MOS_SLIT_MAP";
    curv_traces_tag          = "MOS_CURV_TRACES";
    curv_coeff_tag           = "MOS_CURV_COEFF";
    spatial_map_tag          = "MOS_SPATIAL_MAP";
    slit_location_tag        = "MOS_SLIT_LOCATION";
    disp_residuals_table_tag = "MOS_DISP_RESIDUALS_TABLE";
    detected_lines_tag       = "MOS_DETECTED_LINES";
    delta_image_tag          = "MOS_DELTA_IMAGE";
    flat_sed_tag             = "MOS_FLAT_SED";

    if (cpl_frameset_count_tags(frameset, "MASTER_BIAS") == 0)
        vimos_calmul_exit("Missing required input: MASTER_BIAS");

    if (cpl_frameset_count_tags(frameset, "MASTER_BIAS") > 1)
        vimos_calmul_exit("Too many in input: MASTER_BIAS");

    if (cpl_frameset_count_tags(frameset, "LINE_CATALOG") == 0)
        vimos_calmul_exit("Missing required input: LINE_CATALOG");

    if (cpl_frameset_count_tags(frameset, "LINE_CATALOG") > 1)
        vimos_calmul_exit("Too many in input: LINE_CATALOG");

    nflats = cpl_frameset_count_tags(frameset, flat_tag);

    if (nflats < 1) {
        cpl_msg_error(recipe, "Missing required input: %s", flat_tag);
        vimos_calmul_exit(NULL);
    }

    cpl_msg_indent_less();

    /*
     * Load the master bias
     */
    cpl_msg_info(recipe, "Load master bias...");
    master_bias = dfs_load_image(frameset, "MASTER_BIAS", 
            CPL_TYPE_FLOAT, 0, 1);
    master_bias_err = dfs_load_image(frameset, "MASTER_BIAS", 
            CPL_TYPE_FLOAT, 1, 1);
    cpl_propertylist * master_bias_header = 
            dfs_load_header(frameset, "MASTER_BIAS", 0);
    if (master_bias == NULL || master_bias_err == NULL || master_bias_header == NULL)
        vimos_calmul_exit("Cannot load master bias");

    /*
     * Load flats
     */
    if (nflats > 1)
        cpl_msg_info(recipe, "Load %d flat field frames and sum them...",
                     nflats);
    else
        cpl_msg_info(recipe, "Load flat field exposure...");

    cpl_msg_indent_more();

    header = dfs_load_header(frameset, flat_tag, 0);

    if (header == NULL)
        vimos_calmul_exit("Cannot load flat field frame header");

    alltime = cpl_propertylist_get_double(header, "EXPTIME");

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_calmul_exit("Missing keyword EXPTIME in flat field frame header");

    cpl_propertylist_delete(header);

    for (i = 1; i < nflats; i++) {

        header = dfs_load_header(frameset, NULL, 0);

        if (header == NULL)
            vimos_calmul_exit("Cannot load flat field frame header");

        alltime += cpl_propertylist_get_double(header, "EXPTIME");

        if (cpl_error_get_code() != CPL_ERROR_NONE)
            vimos_calmul_exit("Missing keyword EXPTIME in flat field "
                            "frame header");

        cpl_propertylist_delete(header);

    }

    /*
     * Compute tracing master flat
     */
    cpl_frameset * frameset_copy = cpl_frameset_duplicate(frameset);
    cpl_image * master_flat_raw = 
            dfs_load_image(frameset, flat_tag, CPL_TYPE_FLOAT, 0, 0);
    cpl_propertylist * master_flat_header = 
            dfs_load_header(frameset_copy, flat_tag, 0);
    //Subtract prescan
    cpl_image * master_flat_var = 
        vimos_image_variance_from_detmodel(master_flat_raw, 
                master_flat_header, master_bias_header);
    cpl_image * master_flat_notrim = 
         vimos_subtract_overscan(master_flat_raw, master_flat_var, master_flat_header);
    //Trimm flat
    master_flat = vimos_trimm_preoverscan(master_flat_notrim, master_flat_header);
    //Subtract master bias
    cpl_image_subtract(master_flat, master_bias);

    if (master_flat == NULL)
        vimos_calmul_exit("Cannot load flat field");

    for (i = 1; i < nflats; i++) {
        cpl_image * flat_raw =
                dfs_load_image(frameset, NULL, CPL_TYPE_FLOAT, 0, 0);
        //Subtract prescan
        cpl_propertylist * flat_header = 
                dfs_load_header(frameset_copy, NULL, 0);
        cpl_image * flat_var = 
            vimos_image_variance_from_detmodel(flat_raw, 
                    flat_header, master_bias_header);
        cpl_image * flat_notrim = 
             vimos_subtract_overscan(flat_raw, flat_var, flat_header);
        //Trimm flat
        flat = vimos_trimm_preoverscan(flat_notrim, flat_header);
        //Subtract master bias
        cpl_image_subtract(flat, master_bias);
        if (flat) {
            cpl_image_add(master_flat, flat);
            cpl_image_delete(flat); flat = NULL;
        }
        else
            vimos_calmul_exit("Cannot load flat field");
    }

    /*
     * Get some info from arc lamp header
     */

    header = dfs_load_header(frameset, arc_tag, 0);

    if (header == NULL)
        vimos_calmul_exit("Cannot load arc lamp header");

    instrume = cpl_propertylist_get_string(header, "INSTRUME");
    if (instrume == NULL)
        vimos_calmul_exit("Missing keyword INSTRUME in arc lamp header");
    instrume = cpl_strdup(instrume);

    arctime = cpl_propertylist_get_double(header, "EXPTIME");

    quadrant = cpl_propertylist_get_int(header, "ESO OCS CON QUAD");

    switch (quadrant) {
    case 1:
        key_gris_name = "ESO INS GRIS1 NAME";
        key_gris_id = "ESO INS GRIS1 ID";
        key_filt_name = "ESO INS FILT1 NAME";
        key_filt_id = "ESO INS FILT1 ID";
        key_mask_id = "ESO INS MASK1 ID";
        break;
    case 2:
        key_gris_name = "ESO INS GRIS2 NAME";
        key_gris_id = "ESO INS GRIS2 ID";
        key_filt_name = "ESO INS FILT2 NAME";
        key_filt_id = "ESO INS FILT2 ID";
        key_mask_id = "ESO INS MASK2 ID";
        break;
    case 3:
        key_gris_name = "ESO INS GRIS3 NAME";
        key_gris_id = "ESO INS GRIS3 ID";
        key_filt_name = "ESO INS FILT3 NAME";
        key_filt_id = "ESO INS FILT3 ID";
        key_mask_id = "ESO INS MASK3 ID";
        break;
    case 4:
        key_gris_name = "ESO INS GRIS4 NAME";
        key_gris_id = "ESO INS GRIS4 ID";
        key_filt_name = "ESO INS FILT4 NAME";
        key_filt_id = "ESO INS FILT4 ID";
        key_mask_id = "ESO INS MASK4 ID";
        break;
    }

    grism = cpl_strdup(cpl_propertylist_get_string(header, key_gris_name));

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_calmul_exit("Missing keyword ESO INS GRISn NAME in arc lamp "
                         "frame header");

    cpl_msg_info(recipe, "The grism is: %s", grism);
    cpl_msg_info(recipe, "The spectral multiplexing factor is: %d",
                 1 + (int)cpl_table_get_column_max(allmaskslits, "multiplex"));

/*
    if (!dfs_equal_keyword(frameset, key_gris_id))
        vimos_calmul_exit("Input frames are not from the same grism");

    if (!dfs_equal_keyword(frameset, key_filt_id))
        vimos_calmul_exit("Input frames are not from the same filter");
*/

    gain = cpl_propertylist_get_double(header, "ESO DET OUT1 CONAD");

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_calmul_exit("Missing keyword ESO DET OUT1 CONAD in arc lamp "
                        "frame header");

    cpl_msg_info(recipe, "The gain factor is: %.2f e-/ADU", gain);

    if (wmodemos > 1 && cpl_table_get_column_max(allmaskslits, "curved")) {
        wmodemos = 1;
        cpl_msg_warning(recipe, "There are curved slits on this mask, and "
                        "global distortion solution is not yet supported "
                        "in this case. Setting --wmodemos=1...");
    }

    /*
     * Get the ID of the slit closest to center.
     */

    cslit = mos_slit_closest_to_center(allmaskslits, 0, 0);
    cslit_id = cpl_table_get_int(allmaskslits, "slit_id", cslit, NULL);

    /* Leave the header on for the next step... */

    /*
     * Loading arclamp
     */
    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Load arc lamp exposure...");
    cpl_msg_indent_more();

    cpl_image * spectra_raw = 
            dfs_load_image(frameset, arc_tag, CPL_TYPE_FLOAT, 0, 0);

    if (spectra_raw == NULL)
        vimos_calmul_exit("Cannot load arc lamp exposure");

    //Subtract prescan
    cpl_propertylist * arc_header = 
            dfs_load_header(frameset_copy, arc_tag, 0);
    cpl_image * spectra_var = 
        vimos_image_variance_from_detmodel(spectra_raw, 
                arc_header, master_bias_header);
    cpl_msg_info(recipe, "Subtract pre/overscan..");
    cpl_image * spectra_notrimm =
            vimos_subtract_overscan(spectra_raw, spectra_var, arc_header);
    //Trimm arc
    spectra = vimos_trimm_preoverscan(spectra_notrimm, arc_header);

    cpl_msg_info(recipe, "Remove the master bias...");

    if (cpl_image_subtract(spectra, master_bias))
        vimos_calmul_exit("Cannot remove bias from arc lamp exposure");

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Load input line catalog...");
    cpl_msg_indent_more();

    /*
     * Get the reference lines
     */
    
    lines = vimos_lines_filter(frameset, ignore_lines, used_linesets);
    if(lines == NULL)
        vimos_calmul_exit("Cannot get reference lines");
    nlines = cpl_vector_get_size(lines);

    /*
     * Rotate frames horizontally with red to the right
     */

    cpl_image_turn(spectra, rotate);
    cpl_image_turn(master_flat, rotate);

    ccd_xsize = nx = cpl_image_get_size_x(spectra);      // added...
    ccd_ysize = ny = cpl_image_get_size_y(spectra);

    /*
     * Detecting spectra on the CCD
     */

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Detecting spectra on CCD...");
    cpl_msg_indent_more();

    if (mos_saturation_process(spectra))
        vimos_calmul_exit("Cannot process saturation");

    if (mos_subtract_background(spectra))
        vimos_calmul_exit("Cannot subtract the background");

    failures = 0;

    for (i = 0; i < nchunks; i++) {

        mos_set_multiplex(i);

        refmask = cpl_mask_new(nx, ny);

        checkwave = mos_wavelength_calibration_raw(spectra, lines, dispersion, 
                                                   peakdetection, wradius, 
                                                   wdegree, wreject, reference,
                                                   &startwavelength, 
                                                   &endwavelength,
                                                   NULL, NULL, NULL, NULL, 
                                                   NULL, NULL, refmask, NULL,
                                                   disp_tolerance,
                                                   ratio_tolerance);

        if (checkwave == NULL) {
            failures++;
            continue;
        }

        {
            header = cpl_propertylist_new();
            cpl_propertylist_update_double(header, "CRPIX1", 1.0);
            cpl_propertylist_update_double(header, "CRPIX2", 1.0);
            cpl_propertylist_update_double(header, "CRVAL1",
                                           startwavelength + dispersion/2);
            cpl_propertylist_update_double(header, "CRVAL2", 1.0);
            cpl_propertylist_update_double(header, "CD1_1", dispersion);
            cpl_propertylist_update_double(header, "CD1_2", 0.0);
            cpl_propertylist_update_double(header, "CD2_1", 0.0);
            cpl_propertylist_update_double(header, "CD2_2", 1.0);
            cpl_propertylist_update_string(header, "CTYPE1", "LINEAR");
            cpl_propertylist_update_string(header, "CTYPE2", "PIXEL");

            if (first) {
                first = 0;
                if (dfs_save_image_null(frameset, NULL, parlist, 
                                        spectra_detection_tag,
                                        recipe, version)) {
                    vimos_calmul_exit(NULL);
                }
            }
            if (dfs_save_image_ext(checkwave, spectra_detection_tag, header)) {
                vimos_calmul_exit(NULL);
            }
            cpl_propertylist_delete(header); header = NULL;
        }

        cpl_image_delete(checkwave); checkwave = NULL;

        cpl_msg_debug(recipe, "Locate slits at reference wavelength on CCD...");
        subslits = mos_locate_spectra(refmask);

        if (!subslits) {
            cpl_error_reset();
            cpl_mask_delete(refmask); refmask = NULL;
            failures++;
            continue;
        }

        if (refimage == NULL) {
            refimage = cpl_image_new_from_mask(refmask);
        }
        else {
            dummy = cpl_image_new_from_mask(refmask);
            cpl_image_add(refimage, dummy);
            cpl_image_delete(dummy);
        }

        cpl_mask_delete(refmask); refmask = NULL;

        if (slits == NULL) {
            slits = cpl_table_duplicate(subslits);
        }
        else {
            cpl_table_insert(slits, subslits, cpl_table_get_nrow(slits));
            cpl_table_delete(subslits); subslits = NULL;
        }
    }

    mos_set_multiplex(1);  // Active

    if (failures == nchunks) {

        /*
         * Spectra weren't found anywhere, this is a failure
         */

        vimos_calmul_exit("Wavelength calibration failure.");
    }

    {
        save_header = dfs_load_header(frameset, arc_tag, 0);
        cpl_image_turn(refimage, rotate_back);
        if (dfs_save_image(frameset, refimage, slit_map_tag, NULL,
                           parlist, recipe, version))
            vimos_calmul_exit(NULL);
        cpl_propertylist_delete(save_header); save_header = NULL;
    }

    cpl_image_delete(refimage); refimage = NULL;

    if (slit_ident) {

        /*
         * Slit identification: this recipe may continue even
         * in case of failed identification (i.e., the position table is 
         * not produced, but an error is not set). In case of failure,
         * the spectra would be still extracted, even if they would not
         * be associated to slits on the mask.
         * 
         * The reason for making the slit identification an user option 
         * (via the parameter slit_ident) is to offer the possibility 
         * to avoid identifications that are only apparently successful, 
         * as it would happen in the case of an incorrect slit description 
         * in the data header.
         */

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Slit identification...");
        cpl_msg_indent_more();

        mos_rotate_slits(allmaskslits, -rotate, 0, 0);
        positions = mos_identify_slits(slits, allmaskslits, NULL);

        if (positions) {
            cpl_table_delete(slits);
            slits = positions;

            /*
             * Eliminate slits which are _entirely_ outside the CCD
             */

            cpl_table_and_selected_double(slits, 
                                          "ybottom", CPL_GREATER_THAN, ny-1);
            cpl_table_or_selected_double(slits, 
                                          "ytop", CPL_LESS_THAN, 0);
            cpl_table_erase_selected(slits);

            nslits = cpl_table_get_nrow(slits);

            if (nslits == 0)
                vimos_calmul_exit("No slits found on the CCD");

            cpl_msg_info(recipe, "%d slits are entirely or partially "
                         "contained in CCD", nslits);

        }
        else {
            vimos_calmul_exit("Failure of slit identification");
        }
    }

    mos_set_multiplex(-1);  // Disabled

    /*
     * Now loop on all groups of un-multiplexed spectra.
     */

    ngroups = 1 + cpl_table_get_column_max(allmaskslits, "group");

    globals = (cpl_table **)cpl_calloc(ngroups, sizeof(cpl_table *));

    wcoeff = (double *)cpl_calloc(wdegree + 1, sizeof(double));

    for (i = 0; i < ngroups; i++) {

        cpl_table_select_all(slits);
        cpl_table_and_selected_int(slits, "group", CPL_EQUAL_TO, i);
        subslits = cpl_table_extract_selected(slits);
        
        if (subslits) {
            cpl_propertylist *sort_col = cpl_propertylist_new();
            cpl_propertylist_append_bool(sort_col, "ytop", 1);
            cpl_table_sort(subslits, sort_col);
            cpl_propertylist_delete(sort_col);
        }

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Processing spectra from multiplex group %d:", i);
#ifdef CPL_SIZE_FORMAT
        cpl_msg_info(recipe, 
              "(%" CPL_SIZE_FORMAT " spectra out of %" CPL_SIZE_FORMAT ")", 
              cpl_table_get_nrow(subslits), cpl_table_get_nrow(slits));
#else
        cpl_msg_info(recipe, "(%d spectra out of %d)", 
                     cpl_table_get_nrow(subslits), cpl_table_get_nrow(slits));
#endif
        cpl_msg_indent_more();


        /*
         * Determination of spectral curvature
         */

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Determining spectral curvature...");
        cpl_msg_indent_more();

        cpl_msg_info(recipe, "Tracing master flat field spectra edges...");

        traces = mos_trace_flat(master_flat, subslits, reference, 
                                startwavelength, endwavelength, dispersion);

        if (!traces) {
            vimos_calmul_exit("Tracing failure");
        }

        cpl_msg_info(recipe, "Fitting flat field spectra edges...");
        polytraces = mos_poly_trace(subslits, traces, cdegree);

        if (!polytraces) {
            vimos_calmul_exit("Trace fitting failure");
        }

        if (cmode) {
            cpl_msg_info(recipe, "Computing global spectral curvature...");
            mos_global_trace(subslits, polytraces, cmode);
        }

        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist, curv_traces_tag,
                                   recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_table_ext(traces, curv_traces_tag, NULL)) {
            vimos_calmul_exit(NULL);
        }

        cpl_table_delete(traces); traces = NULL;

        coordinate = cpl_image_new(ccd_xsize, ccd_ysize, CPL_TYPE_FLOAT);
        spatial = mos_spatial_calibration(spectra, subslits, polytraces, 
                                          reference, startwavelength, 
                                          endwavelength, dispersion, 0, 
                                          coordinate);
        if(!spatial)
            vimos_calmul_exit(" Cannot correct spatial distortion. "
                              "No slits could be traced.");
            
        /*
         * QC parameters for flat
         */
        //TODO: Remove this? This was done for the previous method of the master flat  
        if ((cpl_table_and_selected_int(subslits,
                                  "slit_id", CPL_EQUAL_TO, cslit_id) == 1)) {

            double     flux, flux_err;


            cpl_table_select_all(subslits);

            cslit = mos_slit_closest_to_center(subslits, ccd_xsize, ccd_ysize);

            /*
             * Refresh base property list, because previous saving 
             * modified it - e.g., the keyword ARCFILE is missing
             */

            cpl_propertylist_delete(save_header);
            save_header = dfs_load_header(frameset, flat_tag, 0);
            cpl_propertylist_update_int(save_header, 
                                        "ESO PRO DATANCOM", nflats);

            cpl_propertylist_update_string(save_header, "ESO QC DID",
                                           "1.1");        
            cpl_propertylist_set_comment(save_header, "ESO QC DID",
                                         "QC1 dictionary");

            cpl_propertylist_update_double(save_header, 
                                           "ESO PRO WLEN CEN", reference);
            cpl_propertylist_update_double(save_header,
                                           "ESO PRO WLEN INC", dispersion);
            cpl_propertylist_update_double(save_header,
                                       "ESO PRO WLEN START", startwavelength);
            cpl_propertylist_update_double(save_header,
                                           "ESO PRO WLEN END", endwavelength);

            scale = cpl_propertylist_get_double(save_header, 
                                                "ESO TEL FOCU SCALE");

            if (cpl_error_get_code()) {
                cpl_error_reset();
                scale = 1.718;
                cpl_msg_warning(recipe, "Cannot read keyword TEL FOCU SCALE "
                                "(defaulted to %f arcsec/mm)", scale);
            }

            /*
             * QC1 parameters
             */

            keyname = "ESO QC MOS SLIT WIDTH";

            if (cpl_table_has_column(subslits, "ywidth"))
                ywidth = cpl_table_get(subslits, "ywidth", cslit, NULL);

            slit_width = scale * ywidth;

            cpl_propertylist_update_double(save_header, keyname, slit_width);
            cpl_propertylist_set_comment(save_header, keyname,
                                       "Width of slit closest to center (arcsec)");                            

            mos_extract_flux(master_flat, subslits, xwidth, ywidth, 
                             2, gain, &flux, &flux_err);

            flux_err /= alltime; // The master is simply the sum of all flats
            flux     /= alltime;

            cpl_msg_info(recipe, 
                         "Flux at wavelength %.2f: %.2f +/- %.2f ADU/mm^2/s\n",
                         reference, flux, flux_err);

            keyname = "ESO QC MOS FLAT FLUX";

            cpl_propertylist_update_double(save_header, keyname, flux);
            cpl_propertylist_set_comment(save_header, keyname,
                                       "Flux at reference wavelength (ADU/mm^2/s)");

            keyname = "ESO QC MOS FLAT FLUXERR";

            cpl_propertylist_update_double(save_header, keyname, flux_err);
            cpl_propertylist_set_comment(save_header, keyname,
                             "Error on flux at reference wavelength (ADU/mm^2/s)");

            cpl_image_turn(master_flat, rotate_back);
//            if (dfs_save_image(frameset, master_flat, master_screen_flat_tag,
//                               save_header, parlist, recipe, version))
//                vimos_calmul_exit(NULL);

            flat_was_not_saved = 0;
            cpl_image_turn(master_flat, rotate);
//            cpl_image_delete(master_flat); master_flat = NULL;
            cpl_propertylist_delete(save_header); save_header = NULL;

        }

        cpl_table_select_all(subslits);


        /*
         * Final wavelength calibration of spectra having their curvature
         * removed
         */

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Perform final wavelength calibration...");
        cpl_msg_indent_more();

        nx = cpl_image_get_size_x(spatial);
        ny = cpl_image_get_size_y(spatial);

        idscoeff = cpl_table_new(ny);
        restable = cpl_table_new(nlines);
        rainbow = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
        residual = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
        fiterror = (double *)cpl_calloc(ny, sizeof(double));
        fitlines = (int *)cpl_calloc(ny, sizeof(int));

        //Table with positions of the detected lines used for wavelength calibration
        cpl_table * detected_lines = cpl_table_new(1);
        
        rectified = mos_wavelength_calibration_final(spatial, subslits, lines, 
                                                 dispersion, peakdetection, 
                                                 wradius, wdegree, wreject,
                                                 reference, &startwavelength, 
                                                 &endwavelength, fitlines, 
                                                 fiterror, idscoeff, rainbow, 
                                                 residual, restable, 
                                                 detected_lines,
                                                 disp_tolerance,
                                                 ratio_tolerance);

        if (rectified == NULL)
            vimos_calmul_exit("Wavelength calibration failure.");

        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist, disp_residuals_table_tag,
                                   recipe, version)) {
                vimos_calmul_exit(NULL);
            }
            if (dfs_save_image_null(frameset, NULL, parlist, detected_lines_tag,
                                   recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_table_ext(restable, disp_residuals_table_tag, NULL)) 
            vimos_calmul_exit(NULL);

        if (dfs_save_table_ext(detected_lines, detected_lines_tag, NULL))
            vimos_calmul_exit(NULL);

        cpl_table_delete(detected_lines); detected_lines = NULL;

        cpl_table_delete(restable); restable = NULL;

        cpl_table_wrap_double(idscoeff, fiterror, "error"); fiterror = NULL;
        cpl_table_set_column_unit(idscoeff, "error", "pixel");
        cpl_table_wrap_int(idscoeff, fitlines, "nlines"); fitlines = NULL;

        for (j = 0; j < ny; j++)
            if (!cpl_table_is_valid(idscoeff, "c0", j))
                cpl_table_set_invalid(idscoeff, "error", j);

        if (wmodemos > 0) {
            mos_interpolate_wavecalib_slit(idscoeff, subslits, 1, wmodemos - 1);

            cpl_image_delete(rectified);

            rectified = mos_wavelength_calibration(spatial, reference,
                                               startwavelength, endwavelength,
                                               dispersion, idscoeff, 0);

            cpl_image_delete(rainbow);
            rainbow = mos_map_idscoeff(idscoeff, nx, reference,
                                       startwavelength, endwavelength);
        }

        cpl_image_delete(spatial); spatial = NULL;

        delta = mos_map_pixel(idscoeff, reference, startwavelength,
                              endwavelength, dispersion, 2);

        header = cpl_propertylist_new();
        cpl_propertylist_update_double(header, "CRPIX1", 1.0);
        cpl_propertylist_update_double(header, "CRPIX2", 1.0);
        cpl_propertylist_update_double(header, "CRVAL1",
                                       startwavelength + dispersion/2);
        cpl_propertylist_update_double(header, "CRVAL2", 1.0);
        /* cpl_propertylist_update_double(header, "CDELT1", dispersion);
        cpl_propertylist_update_double(header, "CDELT2", 1.0); */
        cpl_propertylist_update_double(header, "CD1_1", dispersion);
        cpl_propertylist_update_double(header, "CD1_2", 0.0);
        cpl_propertylist_update_double(header, "CD2_1", 0.0);
        cpl_propertylist_update_double(header, "CD2_2", 1.0);
        cpl_propertylist_update_string(header, "CTYPE1", "LINEAR");
        cpl_propertylist_update_string(header, "CTYPE2", "PIXEL");

        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist, delta_image_tag,
                                    recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_image_ext(delta, delta_image_tag, header)) {
            vimos_calmul_exit(NULL);
        }

        cpl_image_delete(delta); delta = NULL;
        cpl_propertylist_delete(header); header = NULL;

        mean_rms = mos_distortions_rms(rectified, lines, startwavelength, 
                                       dispersion, 6, 0);

        cpl_msg_info(recipe, "Mean residual: %f pixel", mean_rms);

        mean_rms = cpl_table_get_column_mean(idscoeff, "error");

        cpl_msg_info(recipe, "Mean model accuracy: %f pixel (%f A)", 
                     mean_rms, mean_rms * dispersion);

        restab = mos_resolution_table(rectified, startwavelength, dispersion, 
                                      60000, lines);

        if (restab) {
            cpl_msg_info(recipe, "Mean spectral resolution: %.2f", 
                       cpl_table_get_column_mean(restab, "resolution"));
            cpl_msg_info(recipe, 
                   "Mean reference lines FWHM: %.2f +/- %.2f pixel",
                   cpl_table_get_column_mean(restab, "fwhm") / dispersion,
                   cpl_table_get_column_mean(restab, "fwhm_rms") / dispersion);

            if (i == 0) {
                if (dfs_save_image_null(frameset, NULL, parlist, 
                                        spectral_resolution_tag,
                                        recipe, version)) {
                    vimos_calmul_exit(NULL);
                }
            }

            if (dfs_save_table_ext(restab, spectral_resolution_tag, NULL)) {
                vimos_calmul_exit(NULL);
            }

            cpl_propertylist_delete(qclist); qclist = NULL;

        }
        else
            vimos_calmul_exit("Cannot compute the spectral resolution table");

        cpl_propertylist * disp_coeff_header = cpl_propertylist_new();
        cpl_propertylist_update_double(disp_coeff_header,
                                       "ESO PRO WLEN CEN", reference);
        cpl_propertylist_update_double(disp_coeff_header,
                                       "ESO PRO WLEN INC", dispersion);
        cpl_propertylist_update_double(disp_coeff_header,
                                       "ESO PRO WLEN START", startwavelength);
        cpl_propertylist_update_double(disp_coeff_header,
                                       "ESO PRO WLEN END", endwavelength);
        if (i == 0) {
            if (dfs_save_image_null(frameset, disp_coeff_header, parlist,
                                    disp_coeff_tag,
                                    recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_table_ext(idscoeff, disp_coeff_tag, NULL)) {
            vimos_calmul_exit(NULL);
        }
        cpl_propertylist_delete(disp_coeff_header);

        //Get the mosca wave calib
        //TODO: Use the same format as in FORS. 
        //Better: make a wavelength_calibration constructor with polynomials
        //and a VIMOS function that reads the VIMOS specific idscoeff format.
        cpl_table *idscoeff_lite = cpl_table_duplicate(idscoeff);
        if(cpl_table_has_column(idscoeff_lite, "x"))
            cpl_table_erase_column(idscoeff_lite, "x");
        if(cpl_table_has_column(idscoeff_lite, "y"))
            cpl_table_erase_column(idscoeff_lite, "y");
        mosca::wavelength_calibration wave_cal(idscoeff_lite, reference);
        
        /* Check that the wavelength solution is monotonically increasing */
        for(size_t spa_row = 0 ; spa_row < (size_t)ny; spa_row++)
            if(wave_cal.has_valid_cal((double)spa_row))
                if(!wave_cal.is_monotonical(spa_row, startwavelength,
                                            endwavelength, 
                                            dispersion))
                {
                    std::stringstream error_msg;
                    error_msg <<"The wavelength solution at row "<<spa_row<<
                         " does not increase monotonically, "
                         "which is physically impossible. Try with new parameters.";
                    cpl_msg_warning(cpl_func,"%s",error_msg.str().c_str());
                }

        header = cpl_propertylist_new();
        cpl_propertylist_update_double(header, "CRPIX1", 1.0);
        cpl_propertylist_update_double(header, "CRPIX2", 1.0);
        cpl_propertylist_update_double(header, "CRVAL1", 
                                       startwavelength + dispersion/2);
        cpl_propertylist_update_double(header, "CRVAL2", 1.0);
        /* cpl_propertylist_update_double(header, "CDELT1", dispersion);
        cpl_propertylist_update_double(header, "CDELT2", 1.0); */
        cpl_propertylist_update_double(header, "CD1_1", dispersion);
        cpl_propertylist_update_double(header, "CD1_2", 0.0);
        cpl_propertylist_update_double(header, "CD2_1", 0.0);
        cpl_propertylist_update_double(header, "CD2_2", 1.0);
        cpl_propertylist_update_string(header, "CTYPE1", "LINEAR");
        cpl_propertylist_update_string(header, "CTYPE2", "PIXEL");
        cpl_propertylist_update_int(header, "ESO PRO DATANCOM", 1);

        {

            double flux, flux_err, resol, resol_err;
            int    selected;

            if (grism[0] == 'L') {
                if (grism[3] == 'r') {      /* LR_red    */
                    lambdaHe  = 7065.19;
                    lambdaNe  = 0.0;
                    lambdaAr  = 7723.80;
                    lambdaRed = 9122.97;
                    lambdaYel = 7635.11;
                    lambdaBlu = 5875.62;
                }
                if (grism[3] == 'b') {      /* LR_blue   */
                    lambdaHe  = 5015.68;
                    lambdaNe  = 6598.96;
                    lambdaAr  = 0.0;
                    lambdaRed = 6598.95;
                    lambdaYel = 5015.68;
                    lambdaBlu = 3888.65;
                }
            }

            if (grism[0] == 'M') {          /* MR        */
                lambdaHe  = 7065.19;
                lambdaNe  = 7032.41;
                lambdaAr  = 7723.80;
                lambdaRed = 8264.521;
                lambdaYel = 6678.200;
                lambdaBlu = 5015.675;
            }

            if (grism[0] == 'H') {
                if (grism[3] == 'r') {      /* HR_red    */
                    lambdaHe  = 7065.19;
                    lambdaNe  = 7032.41;
                    lambdaAr  = 7723.80;
                    lambdaRed = 9122.966;
                    lambdaYel = 7948.175;
                    lambdaBlu = 6929.468;
                }
                if (grism[3] == 'o') {      /* HR_orange */
                    lambdaHe  = 7065.19;
                    lambdaNe  = 7032.41;
                    lambdaAr  = 7723.80;
                    lambdaRed = 7948.175;
                    lambdaYel = 6929.468;
                    lambdaBlu = 5875.618;
                }
                if (grism[3] == 'b') {      /* HR_blue   */
                    lambdaHe  = 5015.68;
                    lambdaNe  = 5944.83;
                    lambdaAr  = 0.0;
                    lambdaRed = 6598.953;
                    lambdaYel = 5875.618;
                    lambdaBlu = 5015.675;
                }
            }

            cpl_table_select_all(subslits);

            /*
             * QC1 group header
             */

            if (cpl_table_and_selected_int(subslits,
                                  "slit_id", CPL_EQUAL_TO, cslit_id) == 1) {

                cslit = mos_slit_closest_to_center(subslits, 
                                                   ccd_xsize, ccd_ysize);

                scale = cpl_propertylist_get_double(arc_header, 
                                                    "ESO TEL FOCU SCALE");

                if (cpl_error_get_code()) {
                    cpl_error_reset();
                    scale = 1.718;
                    cpl_msg_warning(recipe, 
                                    "Cannot read keyword TEL FOCU SCALE "
                                    "(defaulted to %f arcsec/mm)", scale);
                }

                /*
                 * QC1 parameters
                 */

                if (cpl_table_has_column(subslits, "ywidth"))
                    ywidth = cpl_table_get(subslits, "ywidth", cslit, NULL);

                slit_width = scale * ywidth;

                if (lambdaHe > 1.) {
                    mos_extract_flux_mapped(rectified, subslits, xwidth, 
                                            ywidth, lambdaHe, startwavelength, 
                                            dispersion, 4, gain, &flux, 
                                            &flux_err);

                    flux     /= arctime;
                    flux_err /= arctime;

                    cpl_msg_info(recipe, 
                                 "Flux of He %.2f: %.2f +/- %.2f ADU/mm^2/s",
                                 lambdaHe, flux, flux_err);

                    he_flux = flux;
                    he_flux_err = flux_err;
                }
                else
                    cpl_msg_warning(recipe, "No He lines in %s spectral range: "
                                    "corresponding QC1 parameters are not "
                                    "computed.", grism);

                if (lambdaNe > 1.) {
                    mos_extract_flux_mapped(rectified, subslits, xwidth, 
                                            ywidth, lambdaNe, 
                                            startwavelength, dispersion, 
                                            4, gain, &flux, &flux_err);

                    flux     /= arctime;
                    flux_err /= arctime;

                    cpl_msg_info(recipe, "Flux of Ne %.2f: %.2f +/- %.2f "
                                 "ADU/mm^2/s", lambdaNe, flux, flux_err);

                    ne_flux = flux;
                    ne_flux_err = flux_err;
                }
                else
                    cpl_msg_warning(recipe, 
                                    "No Ne lines in %s spectral range: "
                                    "corresponding QC1 parameters are not "
                                    "computed.", grism);

                if (lambdaAr > 1.) {
                    mos_extract_flux_mapped(rectified, subslits, xwidth, 
                                            ywidth, lambdaAr, 
                                            startwavelength, dispersion, 
                                            4, gain, &flux, &flux_err);

                    flux     /= arctime;
                    flux_err /= arctime;

                    cpl_msg_info(recipe, 
                                 "Flux of Ar %.2f: %.2f +/- %.2f ADU/mm^2/s",
                                 lambdaAr, flux, flux_err);

                    ar_flux = flux;
                    ar_flux_err = flux_err;
                }
                else
                    cpl_msg_warning(recipe, "No Ar lines in %s spectral range: "
                                    "corresponding QC1 parameters are not "
                                    "computed.", grism);

                /*
                 * IDS coefficients
                 */

                for (j = 0; j <= wdegree; j++) {
                    char  *label = cpl_sprintf("c%d", j);
                    double mcoeff;

                    mcoeff = 0.0;    // Zero by definition when j == 0

                    if (j) {
                        if (mos_median_in_slit(idscoeff, subslits, 
                                               cslit, label, &mcoeff)) {
                            cpl_free(label);
                            break;
                        }
                    }

                    cpl_free(label);

                    wcoeff[j] = mcoeff;
                }
            }

            if (restab) {

                /*
                 * About spectral resolution:
                 */

                cpl_table_and_selected_double(restab, "wavelength", 
                                              CPL_GREATER_THAN, 
                                              lambdaRed - 1.0);
                selected =
                cpl_table_and_selected_double(restab, "wavelength", 
                                              CPL_LESS_THAN, lambdaRed + 1.0);

                if (selected == 1) {
                    cpl_table *one_line = cpl_table_extract_selected(restab);

                    resol = cpl_table_get_double(one_line, 
                                                 "resolution", 0, NULL);
                    resol_err = cpl_table_get_double(one_line, 
                                                 "resolution_rms", 0, NULL);

                    cpl_table_delete(one_line);
                }
                else {
                    resol = 0.0;
                    resol_err = 0.0;
                }
    
                cpl_table_select_all(restab);
    
                cpl_msg_info(recipe, 
                             "Spectral resolution at %.2f: %.2f +/- %.2f",
                             lambdaRed, resol, resol_err);

                r_resol     += resol;
                r_resol_err += resol_err;
    
                cpl_table_and_selected_double(restab, "wavelength",
                                              CPL_GREATER_THAN, 
                                              lambdaYel - 1.0);
                selected =
                cpl_table_and_selected_double(restab, "wavelength",
                                              CPL_LESS_THAN, lambdaYel + 1.0);

                if (selected == 1) {
                    cpl_table *one_line = cpl_table_extract_selected(restab);
    
                    resol = cpl_table_get_double(one_line, 
                                                 "resolution", 0, NULL);
                    resol_err = cpl_table_get_double(one_line, 
                                                 "resolution_rms", 0, NULL);

                    cpl_table_delete(one_line);
                }
                else {
                    resol = 0.0;
                    resol_err = 0.0;
                }

                cpl_table_select_all(restab);
    
                cpl_msg_info(recipe, 
                             "Spectral resolution at %.2f: %.2f +/- %.2f",
                             lambdaYel, resol, resol_err);

                y_resol     += resol;
                y_resol_err += resol_err;

                cpl_table_and_selected_double(restab, "wavelength",
                                              CPL_GREATER_THAN, 
                                              lambdaBlu - 1.0);
                selected =
                cpl_table_and_selected_double(restab, "wavelength",
                                              CPL_LESS_THAN, lambdaBlu + 1.0);

                if (selected == 1) {
                    cpl_table *one_line = cpl_table_extract_selected(restab);
    
                    resol = cpl_table_get_double(one_line, 
                                                 "resolution", 0, NULL);
                    resol_err = cpl_table_get_double(one_line, 
                                                 "resolution_rms", 0, NULL);
    
                    cpl_table_delete(one_line);
                }
                else {
                    resol = 0.0;
                    resol_err = 0.0;
                }

                cpl_table_select_all(restab);

                cpl_msg_info(recipe, 
                             "Spectral resolution at %.2f: %.2f +/- %.2f",
                             lambdaBlu, resol, resol_err);

                b_resol     += resol;
                b_resol_err += resol_err;

                qc_mean_rms += mean_rms;
            }
        }

        cpl_msg_info(recipe, "Computing global distortions model");

        cpl_table_select_all(allmaskslits);
        cpl_table_and_selected_int(allmaskslits, "group", CPL_EQUAL_TO, i);
        maskslits = cpl_table_extract_selected(allmaskslits);

        global = mos_global_distortion(subslits, maskslits, idscoeff,
                                       polytraces, reference);

        if (global && 0) {
            cpl_table *stest;
            cpl_table *ctest;
            cpl_table *dtest;
            cpl_image *itest;

            stest = mos_build_slit_location(global, maskslits, ccd_ysize);

            ctest = mos_build_curv_coeff(global, maskslits, stest);
            if (dfs_save_table(frameset, ctest, "CURVS", NULL,
                               parlist, recipe, version))
                vimos_calmul_exit(NULL);

            itest = mos_spatial_calibration(spectra, stest, ctest,
                                            reference, startwavelength,
                                            endwavelength, dispersion,
                                            0, NULL);
            cpl_table_delete(ctest); ctest = NULL;
            cpl_image_delete(itest); itest = NULL;
            if (dfs_save_table(frameset, stest, "SLITS", NULL,
                               parlist, recipe, version))
                vimos_calmul_exit(NULL);

            dtest = mos_build_disp_coeff(global, stest);
            if (dfs_save_table(frameset, dtest, "DISPS", NULL,
                               parlist, recipe, version))
                vimos_calmul_exit(NULL);

            cpl_table_delete(dtest); dtest = NULL;
            cpl_table_delete(stest); stest = NULL;
        }

        cpl_table_delete(maskslits); maskslits = NULL;

        if (global) {

/* Saving one by one, in different extensions. Now is
 * commented out, only the average is saved (at the end):

            if (i == 0) {
                if (dfs_save_image_null(frameset, parlist, 
                                        global_distortion_tag,
                                        recipe, version)) {
                    vimos_calmul_exit(NULL);
                }
            }

            if (dfs_save_table_ext(global, global_distortion_tag, NULL))
                vimos_calmul_exit(NULL);

            cpl_table_delete(global); global = NULL;
*/

            globals[nglobal] = global;
            global = NULL;
            nglobal++;
        }

        cpl_table_delete(restab); restab = NULL;
        cpl_table_select_all(subslits);

        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist, reduced_lamp_tag,
                                    recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_image_ext(rectified, reduced_lamp_tag, header)) {
            vimos_calmul_exit(NULL);
        }

        cpl_image_delete(rectified); rectified = NULL;
        cpl_propertylist_delete(header); header = NULL;

        {
            if (i == 0) {
                if (dfs_save_image_null(frameset, NULL, parlist, disp_residuals_tag,
                                        recipe, version)) {
                    vimos_calmul_exit(NULL);
                }
            }

            if (dfs_save_image_ext(residual, disp_residuals_tag, NULL)) {
                vimos_calmul_exit(NULL);
            }

            cpl_image_delete(residual); residual = NULL;
        }


        wavemap = mos_map_wavelengths(coordinate, rainbow, subslits, 
                                      polytraces, reference, 
                                      startwavelength, endwavelength, 
                                      dispersion);


        /* Read grism configuration */
        std::auto_ptr<mosca::grism_config> grism_cfg = 
                vimos_grism_config_from_table(grism_table);

        /* Get the detected slit locations */
        //TODO: This is done on rotated images. Change the things to work 
        //on mosca::images with Y_AXIS
        int size_spec = cpl_image_get_size_x(spectra);
        vimos::detected_slits det_slits = 
            vimos::detected_slits_from_tables(subslits, polytraces, size_spec);
        
        /* Get the calibrated slits */
        vimos::calibrated_slits calib_slits(det_slits, wave_cal, *grism_cfg,
                                           (size_t)nx, (size_t)ccd_ysize); 
        for(std::vector<mosca::calibrated_slit>::const_iterator 
                slit_it = calib_slits.begin();
            slit_it != calib_slits.end() ; slit_it++)
        {
            if(!slit_it->has_valid_wavecal())
                cpl_msg_warning(cpl_func, "Slit %d does not contain valid "
                        "wavelength calibration. Skipping it for master flat", 
                        slit_it->slit_id());
        }

        /* Compute master flat.
         * TODO: master flat has already been computed above using the old method
         * Here we use the new method and is the one saved. The other is not yet
         * deleted in case it is used for something else.
         */
        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Perform flat field combination...");
        
        std::auto_ptr<mosca::image> master_flat_d;
        std::auto_ptr<mosca::fiera_config> ccd_config;
        master_flat_d = vimos_calmul_flat_mos_create_master_flat(calib_slits,
                wave_cal, *grism_cfg, master_bias, master_bias_err,
                master_bias_header, frameset, flat_tag);
        if(master_flat_d.get() == 0 || cpl_error_get_code() != CPL_ERROR_NONE)
            vimos_calmul_exit("Cannot combine flat frames");

        /*
         * Perform Flat field normalisation
         */
        
        cpl_msg_info(recipe, "Perform flat field normalisation...");
        
        std::auto_ptr<mosca::image> norm_flat;
        norm_flat.reset(new mosca::image(cpl_image_cast(master_flat_d->get_cpl_image(), 
                                                        CPL_TYPE_FLOAT),
                                         cpl_image_cast(master_flat_d->get_cpl_image_err(), 
                                                        CPL_TYPE_FLOAT), true));

        vimos::flat_normaliser normaliser;
        normaliser.mos_normalise(*norm_flat, wave_cal,
                                 coordinate, calib_slits, subslits, polytraces,
                                 startwavelength, endwavelength,
                                 dispersion, sradius, dradius,
                                 spa_polyorder, disp_nknots, fit_threshold);
        if( cpl_error_get_code() != CPL_ERROR_NONE)
            vimos_calmul_exit("Cannot normalize master flat");

        //Saving the flats
        cpl_propertylist * qc_list = vimos_calmul_flat_qc(*master_flat_d,
                 subslits, nx, ny, nflats, *grism_cfg, alltime, xwidth, ywidth, 
                 gain, scale);
        //cpl_propertylist * qc_list = cpl_propertylist_new();
        vimos_calmul_flats_save(master_flat_d,  norm_flat, frameset, flat_tag, 
                               master_screen_flat_tag, master_norm_flat_tag, 
                               parlist, qc_list, i);

        // Get the spectral shape of the slits and save it
        cpl_image *wave_profiles = normaliser.get_wave_profiles_im_mapped
                (det_slits, wave_cal, startwavelength, endwavelength, dispersion);
        // Get the normalisation factors used 
        std::vector<float> sed_norm;
        std::vector<float> slit_widths;
        std::vector<float> slit_lengths;
        cpl_size i_slit = 0;
        for(vimos::detected_slits::iterator slit_it = det_slits.begin();
            slit_it != det_slits.end(); slit_it++, i_slit++)
        {
            int null;
            slit_lengths.push_back(slit_it->get_length_spatial_corrected());
            if(slit_ident) //If slit identification failed, the pipeline fails.
                slit_widths.push_back(cpl_table_get_double(subslits, "ywidth", i_slit, &null));
        }
        sed_norm = normaliser.get_wave_profiles_norm(alltime, 
                                                     slit_widths, slit_lengths);
        cpl_propertylist * sed_header = cpl_propertylist_new();
        for(size_t ised = 0 ; ised < sed_norm.size();ised++)
        {
            std::ostringstream norm_key;
            norm_key<< "ESO QC FLAT SED_"<<calib_slits[ised].slit_id()<<" NORM ";
            cpl_propertylist_append_float(sed_header, norm_key.str().c_str(),
                                          sed_norm[ised]);
        }
        cpl_propertylist_append_bool(sed_header, "ESO QC FLAT SED CORR_SLITWID",
                                     slit_ident);
        cpl_propertylist_update_double(sed_header, "CRPIX1", 1.0);
        cpl_propertylist_update_double(sed_header, "CRPIX2", 1.0);
        cpl_propertylist_update_double(sed_header, "CRVAL1", 
                                       startwavelength + dispersion/2);
        cpl_propertylist_update_double(sed_header, "CRVAL2", 1.0);
        cpl_propertylist_update_double(sed_header, "CD1_1", dispersion);
        cpl_propertylist_update_double(sed_header, "CD1_2", 0.0);
        cpl_propertylist_update_double(sed_header, "CD2_1", 0.0);
        cpl_propertylist_update_double(sed_header, "CD2_2", 1.0);
        cpl_propertylist_update_string(sed_header, "CTYPE1", "LINEAR");
        cpl_propertylist_update_string(sed_header, "CTYPE2", "SLIT");
        //Saving the flat sed
        std::ostringstream prof_filename_oss;
        prof_filename_oss << flat_sed_tag << ".fits";
        std::string prof_filename = prof_filename_oss.str();
        std::transform(prof_filename.begin(), prof_filename.end(), prof_filename.begin(), ::tolower);
        const cpl_frame * ref_flat_frame = cpl_frameset_find_const(frameset, flat_tag); 
        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist, flat_sed_tag,
                                    recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }
        if (dfs_save_image_ext(wave_profiles, flat_sed_tag,
                              sed_header))
            vimos_calmul_exit(NULL);

        
        /* Saving the rest of the calibrations */
        cpl_image_delete(rainbow); rainbow = NULL;

        save_header = dfs_load_header(frameset, arc_tag, 0);
        vimos_preoverscan ps_scan;
        ps_scan.fix_wcs_trimm(save_header);

        cpl_image_turn(wavemap, rotate_back);

        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist, wavelength_map_tag,
                                    recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_image_ext(wavemap, wavelength_map_tag, NULL)) {
            vimos_calmul_exit(NULL);
        }

        cpl_image_delete(wavemap); wavemap = NULL;

        cpl_image_turn(coordinate, rotate_back);

        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist, spatial_map_tag,
                                    recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_image_ext(coordinate, spatial_map_tag, NULL)) {
            vimos_calmul_exit(NULL);
        }

        cpl_propertylist_delete(save_header); save_header = NULL;

        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist,
                                    curv_coeff_tag,
                                    recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_table_ext(polytraces, curv_coeff_tag, NULL)) {
            vimos_calmul_exit(NULL);
        }
        
        mos_rotate_slits(subslits, rotate, ccd_ysize, ccd_xsize);

        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist,
                                    slit_location_tag,
                                    recipe, version)) {
                vimos_calmul_exit(NULL);
            }
        }

        if (dfs_save_table_ext(subslits, slit_location_tag, NULL)) {
            vimos_calmul_exit(NULL);
        }

        if (cpl_error_get_code()) {
            cpl_msg_error(cpl_error_get_where(), "error: %s", cpl_error_get_message());
            vimos_calmul_exit(NULL);
        }
        
        //Cleanup
        cpl_table_delete(subslits); subslits = NULL;
        cpl_table_delete(polytraces); polytraces = NULL;
        cpl_image_delete(coordinate); coordinate = NULL;
        
    }

    if (flat_was_not_saved) {
        save_header = dfs_load_header(frameset, flat_tag, 0);
        cpl_propertylist_update_int(save_header, "ESO PRO DATANCOM", nflats);

        cpl_image_turn(master_flat, rotate_back);
//        if (dfs_save_image(frameset, master_flat, master_screen_flat_tag,
//                           save_header, parlist, recipe, version))
//            vimos_calmul_exit(NULL);

        cpl_image_delete(master_flat); master_flat = NULL;
        cpl_propertylist_delete(save_header); save_header = NULL;
    }

    /*
     * Average and save all good global distortion tables
     */

    global = mos_average_global_distortion(globals, nglobal, 8.4, 0.2);

    for (i = 0; i < nglobal; i++)
         cpl_table_delete(globals[i]);
    cpl_free(globals);

    if (global) {
        if (dfs_save_table(frameset, global, global_distortion_tag, NULL,
                           parlist, recipe, version))
                    vimos_calmul_exit(NULL);
    }
    else {
        cpl_msg_warning(recipe, 
                        "The global distortion table cannot be computed");
    }


    /*
     * This tail is added to fulfill a request from the QC team:
     * The QC parameters should refer to the whole observation,
     * and not to the multiplex groups taken separately. Therefore
     * all parameters should be moved to the primary array.
     * Unfortunately the QC can only be established when the
     * loop on the multiplex groups is completed. I prefer to
     * reopen the reduced arc lamp product, edit its headers,
     * and rewrite everything to disk. It's a overhead? Maybe.
     * But the maintenance is easier if I handle all this conventional
     * thing here at the end. It will be easier to eliminate it
     * in future, if necessary.
     */

    {
        cpl_frame *frame;
        const char *name    = "mos_arc_spectrum_extracted.fits";
        const char *tmpname = "TMP_mos_arc_spectrum_extracted.fits";
        int         status  = rename(name, tmpname);

        if (status) {
            vimos_calmul_exit("Cannot rename product for QC handling.");
        }

        header = dfs_load_header(frameset, arc_tag, 0);
        cpl_propertylist_update_double(header, "CRPIX1", 1.0);
        cpl_propertylist_update_double(header, "CRPIX2", 1.0);
        cpl_propertylist_update_double(header, "CRVAL1",
                                       startwavelength + dispersion/2);
        cpl_propertylist_update_double(header, "CRVAL2", 1.0);
        cpl_propertylist_update_double(header, "CD1_1", dispersion);
        cpl_propertylist_update_double(header, "CD1_2", 0.0);
        cpl_propertylist_update_double(header, "CD2_1", 0.0);
        cpl_propertylist_update_double(header, "CD2_2", 1.0);
        cpl_propertylist_update_string(header, "CTYPE1", "LINEAR");
        cpl_propertylist_update_string(header, "CTYPE2", "PIXEL");
        cpl_propertylist_update_int(header, "ESO PRO DATANCOM", 1);


        cpl_propertylist_update_string(header, "ESO QC DID",
                                       "1.1");        
        cpl_propertylist_set_comment(header, "ESO QC DID",
                                     "QC1 dictionary");

        cpl_propertylist_update_double(header,
                                       "ESO PRO WLEN CEN", reference);
        cpl_propertylist_update_double(header,
                                       "ESO PRO WLEN INC", dispersion);
        cpl_propertylist_update_double(header,
                                       "ESO PRO WLEN START", startwavelength);
        cpl_propertylist_update_double(header,
                                       "ESO PRO WLEN END", endwavelength);

        /*
         * QC1 parameters
         */

        keyname = "ESO QC MOS SLIT WIDTH";

        cpl_propertylist_update_double(header, keyname, slit_width);
        cpl_propertylist_set_comment(header, keyname,
                                   "Width of slit closest to center (arcsec)");                            


        if (lambdaHe > 1.) {
            keyname = "ESO QC MOS HE LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaHe);
            cpl_propertylist_set_comment(header, keyname,
                         "He arc lamp line for flux determination (Angstrom)");                            

            keyname = "ESO QC MOS HE FLUX";

            cpl_propertylist_update_double(header, keyname, he_flux);
            cpl_propertylist_set_comment(header, keyname,
                         "Flux at chosen He wavelength (ADU/mm^2/s)");                            

            keyname = "ESO QC MOS HE FLUXERR";

            cpl_propertylist_update_double(header, keyname, he_flux_err);
            cpl_propertylist_set_comment(header, keyname,
                         "Error on flux at chosen He wavelength (ADU/mm^2/s)");                            
        }

        if (lambdaNe > 1.) {
            keyname = "ESO QC MOS NE LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaNe);
            cpl_propertylist_set_comment(header, keyname,
                         "Ne arc lamp line for flux determination (Angstrom)");

            keyname = "ESO QC MOS NE FLUX";

            cpl_propertylist_update_double(header, keyname, ne_flux);
            cpl_propertylist_set_comment(header, keyname,
                                  "Flux at chosen Ne wavelength (ADU/mm^2/s)");

            keyname = "ESO QC MOS NE FLUXERR";

            cpl_propertylist_update_double(header, keyname, ne_flux_err);
            cpl_propertylist_set_comment(header, keyname,
                          "Error on flux at chosen Ne wavelength (ADU/mm^2/s)");

        }

        if (lambdaAr > 1.) {
            keyname = "ESO QC MOS AR LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaAr);
            cpl_propertylist_set_comment(header, keyname,
                         "Ar arc lamp line for flux determination (Angstrom)");

            keyname = "ESO QC MOS AR FLUX";

            cpl_propertylist_update_double(header, keyname, ar_flux);
            cpl_propertylist_set_comment(header, keyname,
                                  "Flux at chosen Ar wavelength (ADU/mm^2/s)");

            keyname = "ESO QC MOS AR FLUXERR";

            cpl_propertylist_update_double(header, keyname, ar_flux_err);
            cpl_propertylist_set_comment(header, keyname,
                          "Error on flux at chosen Ar wavelength (ADU/mm^2/s)");
        }

        for (j = 0; j <= wdegree; j++) {
            /* char *label = cpl_sprintf("c%d", j); */
            char *unit; 
            char *comment;

            keyname = cpl_sprintf("ESO QC MOS WAVECAL COEFF%d", j);

            switch (j) {
            case 0:
                unit = cpl_strdup("pixel");
                break;
            case 1:
                unit = cpl_strdup("pixel/Angstrom");
                break;
            default:
                unit = cpl_sprintf("pixel/Angstrom^%d", j);
                break;
            }

            comment = cpl_sprintf("Median coefficient %d of IDS (%s)", j, unit);

            cpl_propertylist_update_double(header, keyname, wcoeff[j]);
            cpl_propertylist_set_comment(header, keyname, comment);

            cpl_free(comment);
            cpl_free(unit);
        }

        cpl_free(wcoeff); wcoeff = NULL;

        /*
         * These parameters are now useless, I set them to zero.
         */

        keyname = "ESO QC MOS REFWAVE MEAN";

        cpl_propertylist_update_double(header, keyname, 0.0);

        keyname = "ESO QC MOS REFWAVE RMS";

        cpl_propertylist_update_double(header, keyname, 0.0);

        /*
         * About spectral resolution:
         */

        keyname = "ESO QC MOS RESOLUTION1 LAMBDA";

        cpl_propertylist_update_double(header, keyname, lambdaRed);
        cpl_propertylist_set_comment(header, keyname, 
               "Line used in spectral resolution determination (Angstrom)");

        keyname = "ESO QC MOS RESOLUTION1";

        r_resol /= ngroups;

        cpl_propertylist_update_double(header, keyname, r_resol);
        cpl_propertylist_set_comment(header, keyname, 
                         "Mean spectral resolution at red end of spectrum");

        keyname = "ESO QC MOS RESOLUTION1 RMS";

        r_resol_err /= ngroups * sqrt(ngroups);

        cpl_propertylist_update_double(header, keyname, r_resol_err);
        cpl_propertylist_set_comment(header, keyname, 
                         "Error on mean spectral resolution");

        keyname = "ESO QC MOS RESOLUTION2 LAMBDA";

        cpl_propertylist_update_double(header, keyname, lambdaYel);
        cpl_propertylist_set_comment(header, keyname, 
              "Line used in spectral resolution determination (Angstrom)");

        keyname = "ESO QC MOS RESOLUTION2";

        y_resol /= ngroups;
        cpl_propertylist_update_double(header, keyname, y_resol);
        cpl_propertylist_set_comment(header, keyname, 
                         "Mean spectral resolution at center of spectrum");

        keyname = "ESO QC MOS RESOLUTION2 RMS";

        y_resol_err /= ngroups * sqrt(ngroups);

        cpl_propertylist_update_double(header, keyname, y_resol_err);
        cpl_propertylist_set_comment(header, keyname, 
                         "Error on mean spectral resolution");

        keyname = "ESO QC MOS RESOLUTION3 LAMBDA";

        cpl_propertylist_update_double(header, keyname, lambdaBlu);
        cpl_propertylist_set_comment(header, keyname, 
               "Line used in spectral resolution determination (Angstrom)");

        keyname = "ESO QC MOS RESOLUTION3";

        b_resol /= ngroups;
        cpl_propertylist_update_double(header, keyname, b_resol);
        cpl_propertylist_set_comment(header, keyname, 
               "Mean spectral resolution at blue end of spectrum");

        keyname = "ESO QC MOS RESOLUTION3 RMS";

        b_resol_err /= ngroups * sqrt(ngroups);

        cpl_propertylist_update_double(header, keyname, b_resol_err);
        cpl_propertylist_set_comment(header, keyname, 
                                     "Error on mean spectral resolution");

        qc_mean_rms /= ngroups;

        keyname = "ESO QC MOS IDS RMS";

        cpl_propertylist_update_double(header, keyname, qc_mean_rms);
        cpl_propertylist_set_comment(header, keyname, 
                           "Mean accuracy of dispersion solution (pixel)");

        frame = cpl_frameset_find(frameset, reduced_lamp_tag);

        cpl_dfs_setup_product_header(header, frame, frameset, parlist,
                                     recipe, version, "PRO-1.15", NULL);

        cpl_propertylist_erase_regexp(header,
        "^ESO DPR |^ARCFILE$|^ORIGFILE$|^ESO PRO CRV |^ESO PRO IDS |^ESO PRO ZERO |^ESO PRO OPT |^ESO PRO CCD |^ESO PRO SKY ", 0);

        cpl_propertylist_save(header, name, CPL_IO_CREATE);
        cpl_propertylist_delete(header);

        for (i = 1; i <= ngroups; i++) {
            cpl_image *image = cpl_image_load(tmpname, CPL_TYPE_FLOAT, 0, i);
            header = cpl_propertylist_load(tmpname, i);
            cpl_image_save(image, name, CPL_BPP_IEEE_FLOAT,
                           header, CPL_IO_EXTEND);
            cpl_image_delete(image);
            cpl_propertylist_delete(header);
        }
        system("rm TMP_mos_arc_spectrum_extracted.fits");
    }

    cpl_propertylist_delete(arc_header);
    cpl_table_delete(slits); slits = NULL;
    cpl_vector_delete(lines); lines = NULL;
    cpl_free(grism); grism = NULL;
    instrume = NULL;
    cpl_image_delete(spectra); spectra = NULL;

    return 0;
}

std::auto_ptr<mosca::image> vimos_calmul_flat_mos_create_master_flat
(vimos::calibrated_slits& calibrated_slits, 
 const mosca::wavelength_calibration& wave_cal,
 const mosca::grism_config& grism_cfg,
 cpl_image *master_bias, cpl_image * master_bias_err,
 cpl_propertylist * master_bias_header,
 cpl_frameset * frameset,
 const char * flat_tag)
{
    const char     * recipe_name = "fors_calib";
    cpl_errorstate   error_prevstate = cpl_errorstate_get();
    std::auto_ptr<mosca::image> master_flat;

    cpl_msg_indent_more();

    /* Get the flat frames */
    cpl_frameset * flatframes = vimos_frameset_extract(frameset, flat_tag);
    size_t nflats = cpl_frameset_get_size(flatframes);

    //Get the variance of the master bias
    cpl_image * master_bias_var = cpl_image_power_create(master_bias_err, 2);

    /* Reading individual raw flats */
    //TODO: This has copy overhead. Substitute with shared_ptr
    std::vector<mosca::image> basiccal_flats;
    for (size_t i_flat = 0; i_flat < nflats; i_flat++)
    {
        cpl_frame * flatframe = cpl_frameset_get_position(flatframes, i_flat);
        cpl_image * flat_raw = cpl_image_load(cpl_frame_get_filename(flatframe),
                CPL_TYPE_FLOAT, 0, 0);
        cpl_propertylist * flat_header = 
                cpl_propertylist_load(cpl_frame_get_filename(flatframe), 0);
        
        if (!flat_raw)
            return master_flat;
        
        /* Create variances map */
        cpl_image * flat_var = 
            vimos_image_variance_from_detmodel(flat_raw, 
                    flat_header, master_bias_header);

        if(!cpl_errorstate_is_equal(error_prevstate))
            return master_flat;

        /* Subtract overscan */
        cpl_image * flat_notrim = 
             vimos_subtract_overscan(flat_raw, flat_var, flat_header);
        if(!cpl_errorstate_is_equal(error_prevstate))
            return master_flat;

        /* Trimm pre/overscan */
        cpl_image * flat_trimmed = vimos_trimm_preoverscan(flat_notrim, flat_header);
        cpl_image * flat_trimmed_var = vimos_trimm_preoverscan(flat_var, flat_header);

        cpl_image_delete(flat_raw);
        cpl_image_delete(flat_var);
        cpl_image_delete(flat_notrim);
        if(!cpl_errorstate_is_equal(error_prevstate))
            return master_flat;

        /* Subtract master bias */

        cpl_image_subtract(flat_trimmed, master_bias);
        cpl_image_add(flat_trimmed_var, master_bias_var);
        if(!cpl_errorstate_is_equal(error_prevstate))
            return master_flat;

        /* Transforming into mosca::image, which takes ownership. The images
         * are rotated, since all he calibrations, wave calib, slits, ...
         * refer to rotated images */
        cpl_image * flat_trimmed_err = flat_trimmed_var;
        cpl_image_power(flat_trimmed_err, 0.5);
        cpl_image_set_bpm(flat_trimmed_err, 
                          cpl_mask_duplicate(cpl_image_get_bpm(flat_trimmed)));
        cpl_image_turn(flat_trimmed, 1);
        cpl_image_turn(flat_trimmed_err, 1);

        mosca::image new_flat(flat_trimmed, flat_trimmed_err, true, mosca::X_AXIS);
        basiccal_flats.push_back(new_flat);
        //Only the structure is freed, the images are taken over by new_flat
        cpl_propertylist_delete(flat_header);
    }
    cpl_image_delete(master_bias_var);

    if(!cpl_errorstate_is_equal(error_prevstate))
    {
        cpl_msg_error(recipe_name, "Could not read the flats");
        return master_flat;
    }   

    /* Computing master flat */
    cpl_msg_info(cpl_func, "Computing master flat");
        
    int smooth_size = 10; //TODO: Hardcoded value!!  
    mosca::reduce_mean reduce_method;
    master_flat = mosca::flat_combine<float, mosca::reduce_mean>
            (basiccal_flats, calibrated_slits, wave_cal, grism_cfg, smooth_size, reduce_method);
    //We multiply to get the same master flat as the simple summed.
    cpl_image_multiply_scalar(master_flat->get_cpl_image(), nflats);
    cpl_image_multiply_scalar(master_flat->get_cpl_image_err(), nflats);

    //Cleanup
    cpl_frameset_delete(flatframes);

    cpl_msg_indent_less();
    return master_flat;
}

int vimos_calmul_flats_save
(std::auto_ptr<mosca::image>& master_flat_d, 
 std::auto_ptr<mosca::image>& norm_flat,
 cpl_frameset * frameset, const char * flat_tag, 
 const char * master_screen_flat_tag, const char * master_norm_flat_tag, 
 cpl_parameterlist * parlist, cpl_propertylist * qc_list, int i)
{
    const char *recipe_name = "vmmoscalib";
    char        version[80];
    snprintf(version, 80, "%s-%s", PACKAGE, PACKAGE_VERSION);
    
    cpl_msg_indent_more();

    /* Saving regular flat */
    cpl_image_turn(master_flat_d->get_cpl_image(), -1);
    cpl_image_turn(master_flat_d->get_cpl_image_err(), -1);
    if (i == 0) {
        if (dfs_save_image_null(frameset, NULL, parlist, master_screen_flat_tag,
                recipe_name, version)) {
            return(-1);
        }
    }
    dfs_save_image_ext(master_flat_d->get_cpl_image(), 
                   master_screen_flat_tag, qc_list);
    dfs_save_image_ext(master_flat_d->get_cpl_image_err(), 
                       master_screen_flat_tag, qc_list);

    if(cpl_error_get_code() != CPL_ERROR_NONE)
        return -1;

    /* Saving normalised flats */
    if(norm_flat.get() != NULL)
    {
        cpl_image_turn(norm_flat->get_cpl_image(), -1);
        cpl_image_turn(norm_flat->get_cpl_image_err(), -1);
        if (i == 0) {
            if (dfs_save_image_null(frameset, NULL, parlist, master_norm_flat_tag,
                    recipe_name, version)) {
                return(-1);
            }
        }
        dfs_save_image_ext(norm_flat->get_cpl_image(), 
                       master_norm_flat_tag, qc_list);
        dfs_save_image_ext(norm_flat->get_cpl_image_err(), 
                           master_norm_flat_tag, NULL);
        if(cpl_error_get_code() != CPL_ERROR_NONE)
            return -1;
    }


    cpl_msg_indent_less();

    return 0;
}


cpl_propertylist * vimos_calmul_flat_qc(mosca::image& master_flat, 
                                      cpl_table * slits,
                                      int nx, int ny, int nflats,
                                      mosca::grism_config& grism_cfg,
                                      double alltime, double xwidth,
                                      double ywidth, double gain, 
                                      double focu_scale)
{
    double     slit_width;
    double     flux, flux_err;
    cpl_propertylist * qc_list = cpl_propertylist_new();
    const char *recipe_name = "vmmoscalib";
    const char *keyname;
    int         cslit;


    cslit = mos_slit_closest_to_center(slits, nx, ny);

    cpl_propertylist_update_string(qc_list, "ESO QC DID",
                                   "1.1");        
    cpl_propertylist_set_comment(qc_list, "ESO QC DID",
                                 "QC1 dictionary");
            
    cpl_propertylist_update_int(qc_list, "ESO PRO DATANCOM", nflats);


    cpl_propertylist_update_double(qc_list, 
                                   "ESO PRO WLEN CEN", grism_cfg.wave_ref());
    cpl_propertylist_update_double(qc_list, 
                                   "ESO PRO WLEN INC", grism_cfg.nominal_dispersion());
    cpl_propertylist_update_double(qc_list, 
                                   "ESO PRO WLEN START", grism_cfg.start_wave());
    cpl_propertylist_update_double(qc_list,
                                   "ESO PRO WLEN END", grism_cfg.end_wave());


    /*
     * QC1 parameters
     */

    keyname = "ESO QC MOS SLIT WIDTH";

    if (cpl_table_has_column(slits, "ywidth"))
        ywidth = cpl_table_get(slits, "ywidth", cslit, NULL);

    slit_width = focu_scale * ywidth;

    cpl_propertylist_update_double(qc_list, keyname, slit_width);
    cpl_propertylist_set_comment(qc_list, keyname,
                               "Width of slit closest to center (arcsec)");                            

    cpl_image * master_flat_f = cpl_image_cast(master_flat.get_cpl_image(),
            CPL_TYPE_FLOAT);
    //We divide first to make sure that the saturation check
    //holds (it should be valid only for individual frames)
    cpl_image_divide_scalar(master_flat_f, nflats);
    mos_extract_flux(master_flat_f, slits, xwidth, ywidth, 
                     2, gain, &flux, &flux_err);

    flux_err /= alltime / nflats; // The master is simply the sum of all flats
    flux     /= alltime / nflats;

    cpl_msg_info(recipe_name, 
                 "Flux at wavelength %.2f: %.2f +/- %.2f ADU/mm^2/s\n",
                 grism_cfg.wave_ref(), flux, flux_err);

    keyname = "ESO QC MOS FLAT FLUX";

    cpl_propertylist_update_double(qc_list, keyname, flux);
    cpl_propertylist_set_comment(qc_list, keyname,
                                 "Flux at reference wavelength (ADU/mm^2/s)");

    keyname = "ESO QC MOS FLAT FLUXERR";

    cpl_propertylist_update_double(qc_list, keyname, flux_err);
    cpl_propertylist_set_comment(qc_list, keyname,
                     "Error on flux at reference wavelength (ADU/mm^2/s)");

    cpl_image_delete(master_flat_f);
    return qc_list;
}

/**@}*/
