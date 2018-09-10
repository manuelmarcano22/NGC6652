/* $Id: vimos_calib_impl.c,v 1.22 2013-10-22 16:57:09 cgarcia Exp $
 *
 * This file is part of the VIMOS Data Reduction Pipeline
 * Copyright (C) 2006 European Southern Observatory
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
 * $Revision: 1.22 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <memory>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <cpl.h>
#include <moses.h>
#include <vimos_dfs.h>
#include <vimos_calib_impl.h>

#include "vimos_grism.h"
#include "vimos_lines.h"
#include "vimos_detected_slits.h"
#include "vimos_calibrated_slits.h"
#include "vimos_overscan.h"
#include "vimos_detmodel.h"
#include "vimos_flat_normalise.h"
#include "flat_combine.h"


std::auto_ptr<mosca::image> vimos_calib_flat_mos_create_master_flat
(vimos::calibrated_slits& calibrated_slits, 
 const mosca::wavelength_calibration& wave_cal,
 const mosca::grism_config& grism_cfg,
 cpl_image *master_bias, cpl_image * master_bias_err,
 cpl_propertylist * master_bias_header,
 cpl_frameset * frameset,
 const char * flat_tag);

int vimos_calib_flats_save
(std::auto_ptr<mosca::image>& master_flat_d, 
 std::auto_ptr<mosca::image>& norm_flat,
 cpl_frameset * frameset, const char * flat_tag, 
 const char * master_screen_flat_tag, const char * master_norm_flat_tag, 
 cpl_parameterlist * parlist, cpl_propertylist * qc_list);

cpl_propertylist * vimos_calib_flat_qc(mosca::image& master_flat, 
                                      cpl_table * slits,
                                      int nx, int ny, int nflats,
                                      mosca::grism_config& grism_cfg,
                                      double alltime, double xwidth,
                                      double ywidth, double gain, 
                                      double focu_scale);


#define vimos_calib_exit(message)             \
{                                             \
if ((const char *)message != NULL) cpl_msg_error(recipe, message);  \
cpl_free(fiterror);                           \
cpl_free(fitlines);                           \
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
cpl_table_delete(restab);                     \
cpl_table_delete(global);                     \
cpl_vector_delete(lines);                     \
cpl_msg_indent_less();                        \
return -1;                                    \
}

#define vimos_calib_exit_memcheck(message)       \
{                                               \
if (message) cpl_msg_info(recipe, message);     \
printf("free instrume (%p)\n", instrume);       \
cpl_free(instrume);                             \
printf("free fiterror (%p)\n", fiterror);       \
cpl_free(fiterror);                             \
printf("free fitlines (%p)\n", fitlines);       \
cpl_free(fitlines);                             \
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
 * @addtogroup vimos_calib_impl
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

int vimos_calib_impl(cpl_frameset *frameset, cpl_parameterlist *parlist)
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
    int         compute_global;
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
    cpl_table        *positions   = NULL;
    cpl_table        *maskslits   = NULL;
    cpl_table        *traces      = NULL;
    cpl_table        *polytraces  = NULL;
    cpl_table        *restab      = NULL;
    cpl_table        *global      = NULL;

    cpl_vector       *lines       = NULL;

    cpl_propertylist *arc_header      = NULL;
    cpl_propertylist *header      = NULL;
    cpl_propertylist *save_header = NULL;
    cpl_propertylist *qclist      = NULL;

    /*
     * Auxiliary variables
     */

    cpl_table  *idscoeff_lss = NULL;
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
    const char *key_gris_name = NULL;
    const char *key_gris_id;
    const char *key_filt_name;
    const char *keyname;
    int         quadrant;
    int         mos;
    int         treat_as_lss = 0;
    int         nslits;
    double     *xpos;
    double      mxpos;
    double      mean_rms;
    double      alltime, arctime;
    int         nflats;
    int         nlines;
    double     *fiterror = NULL;
    int        *fitlines = NULL;
    int         nx, ny;
    double      gain;
    int         ccd_xsize, ccd_ysize;
    int         cslit;
    double      xwidth, ywidth;
    int         rotate = 1;
    int         rotate_back = -1;
    int         i;
    double      scale;

    const char *instrume = NULL;
    char       *grism;


    snprintf(version, 80, "%s-%s", PACKAGE, PACKAGE_VERSION);

    cpl_msg_set_indentation(2);

    /* 
     * Get configuration parameters
     */

    cpl_msg_info(recipe, "Recipe %s configuration parameters:", recipe);
    cpl_msg_indent_more();

    if (cpl_frameset_count_tags(frameset, "CONFIG_TABLE") > 1)
        vimos_calib_exit("Too many in input: CONFIG_TABLE");

    grism_table = dfs_load_table(frameset, "CONFIG_TABLE", 1);

    dispersion = dfs_get_parameter_double(parlist, 
                    "vimos.vmmoscalib.dispersion", grism_table);

    if (dispersion <= 0.0)
        vimos_calib_exit("Invalid spectral dispersion value");

    peakdetection = dfs_get_parameter_double(parlist, 
                    "vimos.vmmoscalib.peakdetection", grism_table);
    if (peakdetection <= 0.0)
        vimos_calib_exit("Invalid peak detection level");

    wdegree = dfs_get_parameter_int(parlist, 
                    "vimos.vmmoscalib.wdegree", grism_table);

    if (wdegree < 1)
        vimos_calib_exit("Invalid polynomial degree");

    if (wdegree > 5)
        vimos_calib_exit("Max allowed polynomial degree is 5");

    wradius = dfs_get_parameter_int(parlist, "vimos.vmmoscalib.wradius", NULL);

    if (wradius < 0)
        vimos_calib_exit("Invalid search radius");

    wreject = dfs_get_parameter_double(parlist, 
                                       "vimos.vmmoscalib.wreject", NULL);

    if (wreject <= 0.0)
        vimos_calib_exit("Invalid rejection threshold");

    wmodelss = dfs_get_parameter_int(parlist, 
                                     "vimos.vmmoscalib.wmodelss", NULL);

    if (wmodelss < 0 || wmodelss > 2)
        vimos_calib_exit("Invalid wavelength solution interpolation mode");

    wmodemos = dfs_get_parameter_int(parlist, 
                                     "vimos.vmmoscalib.wmodemos", NULL);

    if (wmodemos < 0 || wmodemos > 2)
        vimos_calib_exit("Invalid wavelength solution interpolation mode");

    ignore_lines= dfs_get_parameter_string(parlist, 
                   "vimos.vmmoscalib.ignore_lines", NULL);

    used_linesets= dfs_get_parameter_string(parlist, 
                   "vimos.vmmoscalib.used_linesets", NULL);

    cdegree = dfs_get_parameter_int(parlist, 
                    "vimos.vmmoscalib.cdegree", grism_table);

    if (cdegree < 1)
        vimos_calib_exit("Invalid polynomial degree");

    if (cdegree > 5)
        vimos_calib_exit("Max allowed polynomial degree is 5");

    cmode = dfs_get_parameter_int(parlist, "vimos.vmmoscalib.cmode", NULL);

    if (cmode < 0 || cmode > 2)
        vimos_calib_exit("Invalid curvature solution interpolation mode");

    startwavelength = dfs_get_parameter_double(parlist, 
                    "vimos.vmmoscalib.startwavelength", grism_table);
    if (startwavelength > 1.0)
        if (startwavelength < 3000.0 || startwavelength > 13000.0)
            vimos_calib_exit("Invalid wavelength");

    endwavelength = dfs_get_parameter_double(parlist, 
                    "vimos.vmmoscalib.endwavelength", grism_table);
    if (endwavelength > 1.0) {
        if (endwavelength < 3000.0 || endwavelength > 13000.0)
            vimos_calib_exit("Invalid wavelength");
        if (startwavelength < 1.0)
            vimos_calib_exit("Invalid wavelength interval");
    }

    if (startwavelength > 1.0)
        if (endwavelength - startwavelength <= 0.0)
            vimos_calib_exit("Invalid wavelength interval");

    reference = dfs_get_parameter_double(parlist,
                "vimos.vmmoscalib.reference", grism_table);

    if (reference < startwavelength || reference > endwavelength)
        vimos_calib_exit("Invalid reference wavelength");

    slit_ident = dfs_get_parameter_bool(parlist, 
                    "vimos.vmmoscalib.slit_ident", NULL);

    compute_global = slit_ident;
    
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
        vimos_calib_exit("Failure getting the configuration parameters");


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
            vimos_calib_exit("Input frames are not from the same quadrant");

        cpl_frameset_delete(subframeset);
    }

    mos = cpl_frameset_count_tags(frameset, "MOS_ARC_SPECTRUM");

    if (mos == 0)
        vimos_calib_exit("Missing input arc lamp frame");

    if (mos > 1)
        vimos_calib_exit("Just one input arc lamp frame is allowed"); 

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
        vimos_calib_exit("Missing required input: MASTER_BIAS");

    if (cpl_frameset_count_tags(frameset, "MASTER_BIAS") > 1)
        vimos_calib_exit("Too many in input: MASTER_BIAS");

    if (cpl_frameset_count_tags(frameset, "LINE_CATALOG") == 0)
        vimos_calib_exit("Missing required input: LINE_CATALOG");

    if (cpl_frameset_count_tags(frameset, "LINE_CATALOG") > 1)
        vimos_calib_exit("Too many in input: LINE_CATALOG");

    nflats = cpl_frameset_count_tags(frameset, flat_tag);

    if (nflats < 1) {
        cpl_msg_error(recipe, "Missing required input: %s", flat_tag);
        vimos_calib_exit(NULL);
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
        vimos_calib_exit("Cannot load master bias");

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
        vimos_calib_exit("Cannot load flat field frame header");

    alltime = cpl_propertylist_get_double(header, "EXPTIME");

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_calib_exit("Missing keyword EXPTIME in flat field frame header");

    cpl_propertylist_delete(header);

    for (i = 1; i < nflats; i++) {

        header = dfs_load_header(frameset, NULL, 0);

        if (header == NULL)
            vimos_calib_exit("Cannot load flat field frame header");

        alltime += cpl_propertylist_get_double(header, "EXPTIME");

        if (cpl_error_get_code() != CPL_ERROR_NONE)
            vimos_calib_exit("Missing keyword EXPTIME in flat field "
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
    if(cpl_error_get_code() != CPL_ERROR_NONE)
        throw std::invalid_argument("Could not get RON from master bias"
            " (missing QC DET OUT? RON keywords)");
    cpl_image * master_flat_notrim = 
         vimos_subtract_overscan(master_flat_raw, master_flat_var, master_flat_header);
    //Trimm flat
    master_flat = vimos_trimm_preoverscan(master_flat_notrim, master_flat_header);
    //Subtract master bias
    cpl_image_subtract(master_flat, master_bias);

    if (master_flat == NULL)
        vimos_calib_exit("Cannot load flat field");

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
            vimos_calib_exit("Cannot load flat field");
    }

    /*
     * Get some info from arc lamp header
     */

    arc_header = dfs_load_header(frameset, arc_tag, 0);

    if (arc_header == NULL)
        vimos_calib_exit("Cannot load arc lamp header");

    instrume = cpl_propertylist_get_string(arc_header, "INSTRUME");
    if (instrume == NULL)
        vimos_calib_exit("Missing keyword INSTRUME in arc lamp header");
    instrume = cpl_strdup(instrume);

    arctime = cpl_propertylist_get_double(arc_header, "EXPTIME");

    quadrant = cpl_propertylist_get_int(arc_header, "ESO OCS CON QUAD");

    switch (quadrant) {
    case 1:
        key_gris_name = "ESO INS GRIS1 NAME";
        key_gris_id = "ESO INS GRIS1 ID";
        key_filt_name = "ESO INS FILT1 NAME";
        break;
    case 2:
        key_gris_name = "ESO INS GRIS2 NAME";
        key_gris_id = "ESO INS GRIS2 ID";
        key_filt_name = "ESO INS FILT2 NAME";
        break;
    case 3:
        key_gris_name = "ESO INS GRIS3 NAME";
        key_gris_id = "ESO INS GRIS3 ID";
        key_filt_name = "ESO INS FILT3 NAME";
        break;
    case 4:
        key_gris_name = "ESO INS GRIS4 NAME";
        key_gris_id = "ESO INS GRIS4 ID";
        key_filt_name = "ESO INS FILT4 NAME";
        break;
    }

    grism = cpl_strdup(cpl_propertylist_get_string(arc_header, key_gris_name));

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_calib_exit("Missing keyword ESO INS GRISn NAME in arc lamp "
                         "frame header");

    cpl_msg_info(recipe, "The grism is: %s", grism);

/*
    if (!dfs_equal_keyword(frameset, key_gris_id))
        vimos_calib_exit("Input frames are not from the same grism");

    if (!dfs_equal_keyword(frameset, key_filt_id))
        vimos_calib_exit("Input frames are not from the same filter");
*/

    gain = cpl_propertylist_get_double(arc_header, "ESO DET OUT1 CONAD");

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_calib_exit("Missing keyword ESO DET OUT1 CONAD in arc lamp "
                        "frame header");

    cpl_msg_info(recipe, "The gain factor is: %.2f e-/ADU", gain);

    cpl_msg_info(recipe, "Produce mask slit position table...");

    maskslits = mos_load_slits_vimos(arc_header, 0);

    if (wmodemos > 1 && cpl_table_get_column_max(maskslits, "curved")) {
        wmodemos = 1;
        cpl_msg_warning(recipe, "There are curved slits on this mask, and "
                        "global distortion solution is not yet supported "
                        "in this case. Setting --wmodemos=1...");
    }

    /*
     * Get the (default) width of the slit closest to center,
     * in case slit identification is off (or failes).
     */

    cslit = mos_slit_closest_to_center(maskslits, 0, 0);

    xwidth = cpl_table_get(maskslits, "xwidth", cslit, NULL);
    ywidth = cpl_table_get(maskslits, "ywidth", cslit, NULL);

    /*
     * Check if all slits have the same X offset: in such case, 
     * treat the observation as a long-slit one!
     */

    mxpos = cpl_table_get_column_median(maskslits, "ytop");
    xpos = cpl_table_get_data_double(maskslits, "ytop");
    nslits = cpl_table_get_nrow(maskslits);

    treat_as_lss = 1;
    for (i = 0; i < nslits; i++) {
        if (fabs(mxpos-xpos[i]) > 0.3) {
            treat_as_lss = 0;
            break;
        }
        cpl_propertylist * ref_slits = cpl_propertylist_new();
        cpl_propertylist_copy_property_regexp(ref_slits, arc_header, "ESO INS REF", 0);
        bool reference_slits_present = (cpl_propertylist_get_size(ref_slits) > 0); 
        if(nslits == 1 && reference_slits_present)
            treat_as_lss = 0;
        cpl_propertylist_delete(ref_slits);
    }

    if (treat_as_lss) {
        cpl_msg_warning(recipe, "All MOS slits have the same offset: %.2f mm\n"
                        "The long-slit data reduction strategy is applied!", 
                        mxpos);
        cpl_table_delete(maskslits); maskslits = NULL;
    }

    if (slit_ident == 0) {
        cpl_table_delete(maskslits); maskslits = NULL;
    }


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
        vimos_calib_exit("Cannot load arc lamp exposure");

    /*
     * Subtract prescan
     */
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
        vimos_calib_exit("Cannot remove bias from arc lamp exposure");

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Load input line catalog...");
    cpl_msg_indent_more();

    /*
     * Get the reference lines
     */
    
    lines = vimos_lines_filter(frameset, ignore_lines, used_linesets);
    if(lines == NULL)
        vimos_calib_exit("Cannot get reference lines");
    nlines = cpl_vector_get_size(lines);

    /*
     * Rotate frames horizontally with red to the right
     */

    cpl_image_turn(spectra, rotate);
    cpl_image_turn(master_flat, rotate);

    ccd_xsize = nx = cpl_image_get_size_x(spectra);      // added...
    ccd_ysize = ny = cpl_image_get_size_y(spectra);

    if (treat_as_lss) {

        /*
         * In the case of LSS data, find first a "one slit"
         * solution. This will be later on split into many-slits
         * solutions. This is done for greater accuracy.
         */

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Perform wavelength calibration...");
        cpl_msg_indent_more();

        wavemap = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
        residual = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);

        fiterror = (double*)cpl_calloc(ny, sizeof(double));
        fitlines = (int *)cpl_calloc(ny, sizeof(int));
        idscoeff = cpl_table_new(ny);
        refmask  = cpl_mask_new(nx, ny);

	if (mos_saturation_process(spectra))
	    vimos_calib_exit("Cannot process saturation");

	if (mos_subtract_background(spectra))
	    vimos_calib_exit("Cannot subtract the background");

        rectified = mos_wavelength_calibration_raw(spectra, lines, dispersion,
                                                   peakdetection, wradius,
                                                   wdegree, wreject, reference,
                                                   &startwavelength, 
                                                   &endwavelength, fitlines, 
                                                   fiterror, idscoeff, wavemap,
                                                   residual, NULL, 
                                                   refmask, NULL,
                                                   disp_tolerance,
                                                   ratio_tolerance);

        if (rectified == NULL)
            vimos_calib_exit("Wavelength calibration failure.");

        if (!cpl_table_has_valid(idscoeff, "c0"))
            vimos_calib_exit("Wavelength calibration failure.");

        if (wmodelss) {
            cpl_binary *mdata;
            float      *wdata;
            int         npix = nx * ny;

            cpl_image_delete(rectified); rectified = NULL;
            cpl_image_delete(wavemap); wavemap = NULL;
            mos_interpolate_wavecalib(idscoeff, wavemap, wmodelss, 2);
            wavemap = mos_map_idscoeff(idscoeff, nx, reference,
                                       startwavelength, endwavelength);
            rectified = mos_wavelength_calibration(spectra, reference,
                                                   startwavelength, 
                                                   endwavelength, dispersion, 
                                                   idscoeff, 0);
            mdata = cpl_mask_get_data(refmask);
            wdata = cpl_image_get_data_float(wavemap);

            for (i = 0; i < npix; i++) {
                mdata[i] = CPL_BINARY_0;
            }
            for (i = 0; i < npix; i++) {
                if (fabs(wdata[i] - reference) < 2*dispersion) {
                    mdata[i] = CPL_BINARY_1;
                }
            }
        }

        /*
         * This is necessary to move on to a many-slits solution
         */

        mos_refmask_find_gaps(refmask, master_flat, 5000.);


        cpl_table_wrap_double(idscoeff, fiterror, "error"); fiterror = NULL;
        cpl_table_set_column_unit(idscoeff, "error", "pixel");
        cpl_table_wrap_int(idscoeff, fitlines, "nlines"); fitlines = NULL;

        for (i = 0; i < ny; i++)
            if (!cpl_table_is_valid(idscoeff, "c0", i))
                cpl_table_set_invalid(idscoeff, "error", i);

        idscoeff_lss = idscoeff;

    }
    else {

        /*
         * Here the generic MOS calibration is carried out.
         */

        /*
         * Detecting spectra on the CCD
         */

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Detecting spectra on CCD...");
        cpl_msg_indent_more();

        ccd_xsize = nx = cpl_image_get_size_x(spectra);
        ccd_ysize = ny = cpl_image_get_size_y(spectra);

        refmask = cpl_mask_new(nx, ny);

        if (mos_saturation_process(spectra))
	    vimos_calib_exit("Cannot process saturation");

        if (mos_subtract_background(spectra))
	    vimos_calib_exit("Cannot subtract the background");

        // Begin code to smooth image

        if (0) {
            cpl_mask *kernel = cpl_mask_new(3, 7);
            cpl_mask_not(kernel);
            dummy = cpl_image_duplicate(spectra);
            cpl_image_filter_mask(dummy, spectra, kernel, CPL_FILTER_MEDIAN,
                                  CPL_BORDER_FILTER);
            cpl_image_delete(spectra);
            spectra = dummy;
            cpl_mask_delete(kernel);
        }

        // End code to smooth image

        checkwave = mos_wavelength_calibration_raw(spectra, lines, dispersion, 
                                                   peakdetection, wradius, 
                                                   wdegree, wreject, reference,
                                                   &startwavelength, 
                                                   &endwavelength,
                                                   NULL, NULL, NULL, NULL, 
                                                   NULL, NULL, refmask, NULL,
                                                   disp_tolerance,
                                                   ratio_tolerance);

        if (checkwave == NULL)
            vimos_calib_exit("Wavelength calibration failure.");

        /*
         * Save check image to disk
         */

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

        if (dfs_save_image(frameset, checkwave, spectra_detection_tag, 
                           header, parlist, recipe, version))
            vimos_calib_exit(NULL);

        cpl_image_delete(checkwave); checkwave = NULL;
        cpl_propertylist_delete(header); header = NULL;

    }

    cpl_msg_info(recipe, "Locate slits at reference wavelength on CCD...");
    slits = mos_locate_spectra(refmask);

    if (!slits) {
        cpl_msg_error(cpl_error_get_where(), "error: %s", cpl_error_get_message());
        vimos_calib_exit("No slits could be detected!");
    }

    refimage = cpl_image_new_from_mask(refmask);
    cpl_mask_delete(refmask); refmask = NULL;

    {
        save_header = dfs_load_header(frameset, arc_tag, 0);
        cpl_image_turn(refimage, rotate_back);
        if (dfs_save_image(frameset, refimage, slit_map_tag, NULL,
                           parlist, recipe, version))
            vimos_calib_exit(NULL);
        cpl_propertylist_delete(save_header); save_header = NULL;
    }

    cpl_image_delete(refimage); refimage = NULL;

    if (slit_ident && !treat_as_lss) {

        /*
         * Attempt slit identification: this recipe may continue even
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
        cpl_msg_info(recipe, "Attempt slit identification (optional)...");
        cpl_msg_indent_more();

        mos_rotate_slits(maskslits, -rotate, 0, 0);
        positions = mos_identify_slits(slits, maskslits, NULL);

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
                vimos_calib_exit("No slits found on the CCD");

            cpl_msg_info(recipe, "%d slits are entirely or partially "
                         "contained in CCD", nslits);

        }
        else {
            compute_global = 0;
            cpl_msg_info(recipe, "Global distortion model cannot be computed");
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                vimos_calib_exit(NULL);
            }
        }
    }
    else if(slit_ident && treat_as_lss)
    {
        /* In this case we perform a slit identification based purely on 
         * the slits X positions rather than pattern matching.
         */
        cpl_msg_info(cpl_func, "Performing LSS slit identification");
        
        /* First we need to get all the slits, including the reference ones */
        cpl_table * maskslits_inclref = mos_load_slits_vimos(arc_header, 1);

        mos_rotate_slits(maskslits_inclref, -rotate, 0, 0);

        positions = mos_identify_slits_linear(slits, maskslits_inclref);

        if (positions) {
            cpl_table_delete(slits);
            slits = positions;

            /*
             * Very ugly kludge for standard star masks which have a duplicated
             * slit id 80 in SLIT2 ID and SLIT6 ID 
             * This is a problem with the mask generation.
             * TODO: Remove this if all the data in the archive is fixed eventually.
             */
            for(int i_slit = 0; i_slit < cpl_table_get_nrow(slits); i_slit++)
            {
                int null;
                if(cpl_table_get_int(slits, "slit_id", i_slit, &null) == 80)
                {
                    cpl_table_set_int(slits, "slit_id", i_slit, 0);
                    break;
                }
            }

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
                vimos_calib_exit("No slits found on the CCD");

            cpl_msg_info(recipe, "%d slits are entirely or partially "
                         "contained in CCD", nslits);

        }
    }


    /*
     * Determination of spectral curvature
     */

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Determining spectral curvature...");
    cpl_msg_indent_more();

    cpl_msg_info(recipe, "Tracing master flat field spectra edges...");
    traces = mos_trace_flat(master_flat, slits, reference, 
                            startwavelength, endwavelength, dispersion);

    if (!traces)
        vimos_calib_exit("Tracing failure");

    cpl_msg_info(recipe, "Fitting flat field spectra edges...");
    polytraces = mos_poly_trace(slits, traces, cdegree);

    if (!polytraces)
        vimos_calib_exit("Trace fitting failure");

    if (cmode) {
        cpl_msg_info(recipe, "Computing global spectral curvature model...");
        mos_global_trace(slits, polytraces, cmode);
    }

    if (dfs_save_table(frameset, traces, curv_traces_tag, NULL, parlist,
                       recipe, version))
        vimos_calib_exit(NULL);

    cpl_table_delete(traces); traces = NULL;

    coordinate = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
    spatial = mos_spatial_calibration(spectra, slits, polytraces, reference, 
                                      startwavelength, endwavelength, 
                                      dispersion, 0, coordinate);

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

    rectified = mos_wavelength_calibration_final(spatial, slits, lines, 
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
        vimos_calib_exit("Wavelength calibration failure.");

    if (dfs_save_table(frameset, restable, disp_residuals_table_tag, NULL,
                       parlist, recipe, version))
        vimos_calib_exit(NULL);
    cpl_table_delete(restable); restable = NULL;

    if(dfs_save_table(frameset, detected_lines, detected_lines_tag, NULL,
                        parlist, recipe, version))
        vimos_calib_exit(NULL);
    cpl_table_delete(detected_lines); detected_lines = NULL;

    cpl_table_wrap_double(idscoeff, fiterror, "error"); fiterror = NULL;
    cpl_table_set_column_unit(idscoeff, "error", "pixel");
    cpl_table_wrap_int(idscoeff, fitlines, "nlines"); fitlines = NULL;

    if (treat_as_lss) {

        //const char *clab[6] = {"c0", "c1", "c2", "c3", "c4", "c5"};

/*
        int    *position = cpl_table_get_data_int   (slits, "position");
        int    *length   = cpl_table_get_data_int   (slits, "length");
        double *ytop     = cpl_table_get_data_double(slits, "ytop");
        double *ybottom  = cpl_table_get_data_double(slits, "ybottom");
        int     j, k, dr, irow;
*/

        /*
         * Replace the LSS solutions to the poorer MOS solutions
         */

        mos_interpolate_wavecalib(idscoeff, NULL, wmodelss, 2);

/*
        for (j = 0, i = ny - 1; i >= 0; i--) {
            if (i < position[j]) {
                ++j;
            }
            dr = position[j] + length[j] - i - 1;
            irow = floor(ytop[j] - dr*(ytop[j] - ybottom[j])/length[j] + 0.5);
            for (k = 0; k <= wdegree; k++) {
                cpl_table_set_double(idscoeff, clab[k], i, 
                cpl_table_get_double(idscoeff_lss, clab[k], irow, NULL));
            }
            cpl_table_set_double(idscoeff, "error", i,
            cpl_table_get_double(idscoeff_lss, "error", irow, NULL));
            cpl_table_set_int(idscoeff, "nlines", i,
            cpl_table_get_int(idscoeff_lss, "nlines", irow, NULL));
        }
*/

        cpl_table_delete(idscoeff_lss);

        cpl_image_delete(rectified);

        rectified = mos_wavelength_calibration(spatial, reference,
                                               startwavelength, endwavelength,
                                               dispersion, idscoeff, 0);
    }
    else {
        for (i = 0; i < ny; i++)
            if (!cpl_table_is_valid(idscoeff, "c0", i))
                cpl_table_set_invalid(idscoeff, "error", i);

        if (wmodemos > 0) {
            mos_interpolate_wavecalib_slit(idscoeff, slits, 1, wmodemos - 1);

            cpl_image_delete(rectified);

            rectified = mos_wavelength_calibration(spatial, reference,
                                               startwavelength, endwavelength,
                                               dispersion, idscoeff, 0);

            cpl_image_delete(rainbow);
            rainbow = mos_map_idscoeff(idscoeff, nx, reference,
                                       startwavelength, endwavelength);
        }
    }

    cpl_image_delete(spatial); spatial = NULL;

    delta = mos_map_pixel(idscoeff, reference, startwavelength,
                          endwavelength, dispersion, 2);

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

    if (dfs_save_image(frameset, delta, delta_image_tag,
                       header, parlist, recipe, version))
        vimos_calib_exit(NULL);

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
        cpl_msg_info(recipe, "Mean reference lines FWHM: %.2f +/- %.2f pixel",
                   cpl_table_get_column_mean(restab, "fwhm") / dispersion,
                   cpl_table_get_column_mean(restab, "fwhm_rms") / dispersion);


        if (dfs_save_table(frameset, restab, spectral_resolution_tag, NULL,
                           parlist, recipe, version))
            vimos_calib_exit(NULL);

    }
    else
        vimos_calib_exit("Cannot compute the spectral resolution table");

    cpl_vector_delete(lines); lines = NULL;

    /*
     * Global distortion models
     */

    if (compute_global && !treat_as_lss) {

        cpl_msg_info(recipe, "Computing global distortions model");
        global = mos_global_distortion(slits, maskslits, idscoeff, 
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
                vimos_calib_exit(NULL);

            itest = mos_spatial_calibration(spectra, stest, ctest, 
                                            reference, startwavelength, 
                                            endwavelength, dispersion, 
                                            0, NULL);
            cpl_table_delete(ctest); ctest = NULL;
            cpl_image_delete(itest); itest = NULL;
            if (dfs_save_table(frameset, stest, "SLITS", NULL,
                               parlist, recipe, version))
                vimos_calib_exit(NULL);

            dtest = mos_build_disp_coeff(global, stest);
            if (dfs_save_table(frameset, dtest, "DISPS", NULL,
                               parlist, recipe, version))
                vimos_calib_exit(NULL);

            cpl_table_delete(dtest); dtest = NULL;
            cpl_table_delete(stest); stest = NULL;
        }

        if (global) {
            if (dfs_save_table(frameset, global, global_distortion_tag, NULL,
                               parlist, recipe, version))
                vimos_calib_exit(NULL);
            cpl_table_delete(global); global = NULL;
        }

//        cpl_image_delete(spectra); spectra = NULL;
        cpl_table_delete(maskslits); maskslits = NULL;
    }

    header = dfs_load_header(frameset, arc_tag, 0);
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
    scale = cpl_propertylist_get_double(header, "ESO TEL FOCU SCALE");

    if (cpl_error_get_code()) {
        cpl_error_reset();
        scale = 1.718;
        cpl_msg_warning(recipe, "Cannot read keyword TEL FOCU SCALE "
                        "(defaulted to %f arcsec/mm)", scale);
    }


    {

        double slit_width;
        float  lambdaHe = 0;
        float  lambdaNe = 0;
        float  lambdaAr = 0;
        float  lambdaRed = 0;
        float  lambdaYel = 0;
        float  lambdaBlu = 0;
        double flux, flux_err, resol, resol_err;
        int    selected;


        cslit = mos_slit_closest_to_center(slits, nx, ny);

        /*
         * QC1 group header
         */

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

        if (cpl_table_has_column(slits, "ywidth"))
            ywidth = cpl_table_get(slits, "ywidth", cslit, NULL);

        slit_width = scale * ywidth;

        cpl_propertylist_update_double(header, keyname, slit_width);
        cpl_propertylist_set_comment(header, keyname,
                                   "Width of slit closest to center (arcsec)");                            

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

        if (lambdaHe > 1.) {
            mos_extract_flux_mapped(rectified, slits, xwidth, ywidth, 
                                    lambdaHe, startwavelength, dispersion, 
                                    4, gain, &flux, &flux_err);

            flux     /= arctime;
            flux_err /= arctime;

            cpl_msg_info(recipe, "Flux of He %.2f: %.2f +/- %.2f ADU/mm^2/s",
                         lambdaHe, flux, flux_err);

            keyname = "ESO QC MOS HE LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaHe);
            cpl_propertylist_set_comment(header, keyname,
                         "He arc lamp line for flux determination (Angstrom)");                            

            keyname = "ESO QC MOS HE FLUX";

            cpl_propertylist_update_double(header, keyname, flux);
            cpl_propertylist_set_comment(header, keyname,
                         "Flux at chosen He wavelength (ADU/mm^2/s)");                            

            keyname = "ESO QC MOS HE FLUXERR";

            cpl_propertylist_update_double(header, keyname, flux_err);
            cpl_propertylist_set_comment(header, keyname,
                         "Error on flux at chosen He wavelength (ADU/mm^2/s)");                            

        }
        else
            cpl_msg_warning(recipe, 
                            "No He lines in %s spectral range: corresponding "
                            "QC1 parameters are not computed.", grism);

        if (lambdaNe > 1.) {
            mos_extract_flux_mapped(rectified, slits, xwidth, ywidth,
                                    lambdaNe, 
                                    startwavelength, dispersion, 
                                    4, gain, &flux, &flux_err);

            flux     /= arctime;
            flux_err /= arctime;

            cpl_msg_info(recipe, "Flux of Ne %.2f: %.2f +/- %.2f ADU/mm^2/s",
                         lambdaNe, flux, flux_err);

            keyname = "ESO QC MOS NE LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaNe);
            cpl_propertylist_set_comment(header, keyname,
                         "Ne arc lamp line for flux determination (Angstrom)");

            keyname = "ESO QC MOS NE FLUX";

            cpl_propertylist_update_double(header, keyname, flux);
            cpl_propertylist_set_comment(header, keyname,
                                  "Flux at chosen Ne wavelength (ADU/mm^2/s)");

            keyname = "ESO QC MOS NE FLUXERR";

            cpl_propertylist_update_double(header, keyname, flux_err);
            cpl_propertylist_set_comment(header, keyname,
                          "Error on flux at chosen Ne wavelength (ADU/mm^2/s)");

        }
        else
            cpl_msg_warning(recipe, 
                            "No Ne lines in %s spectral range: corresponding "
                            "QC1 parameters are not computed.", grism);

        if (lambdaAr > 1.) {
            mos_extract_flux_mapped(rectified, slits, xwidth, ywidth,
                                    lambdaAr, 
                                    startwavelength, dispersion, 
                                    4, gain, &flux, &flux_err);
//            mos_extract_flux(spectra, slits, 3, gain, &flux, &flux_err);

            flux     /= arctime;
            flux_err /= arctime;

            cpl_msg_info(recipe, "Flux of Ar %.2f: %.2f +/- %.2f ADU/mm^2/s",
                         lambdaAr, flux, flux_err);

            keyname = "ESO QC MOS AR LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaAr);
            cpl_propertylist_set_comment(header, keyname,
                         "Ar arc lamp line for flux determination (Angstrom)");

            keyname = "ESO QC MOS AR FLUX";

            cpl_propertylist_update_double(header, keyname, flux);
            cpl_propertylist_set_comment(header, keyname,
                                  "Flux at chosen Ar wavelength (ADU/mm^2/s)");

            keyname = "ESO QC MOS AR FLUXERR";

            cpl_propertylist_update_double(header, keyname, flux_err);
            cpl_propertylist_set_comment(header, keyname,
                          "Error on flux at chosen Ar wavelength (ADU/mm^2/s)");
        }
        else
            cpl_msg_warning(recipe, 
                            "No Ar lines in %s spectral range: corresponding "
                            "QC1 parameters are not computed.", grism);

        /*
         * IDS coefficients
         */

        for (i = 0; i <= wdegree; i++) {
            char  *label = cpl_sprintf("c%d", i);
            char  *unit;
            char  *comment;
            double mcoeff;


            mcoeff = 0.0;    // Zero by definition when i == 0
            if (i) {
                if (mos_median_in_slit(idscoeff, slits, 
                                       cslit, label, &mcoeff)) {
                    cpl_free(label);
                    break;
                }
            }

            keyname = cpl_sprintf("ESO QC MOS WAVECAL COEFF%d", i);

            switch (i) {
            case 0:
                unit = cpl_strdup("pixel");
                break;
            case 1:
                unit = cpl_strdup("pixel/Angstrom");
                break;
            default:
                unit = cpl_sprintf("pixel/Angstrom^%d", i);
                break;
            }

            comment = cpl_sprintf("Median coefficient %d of IDS (%s)", i, unit);

            cpl_propertylist_update_double(header, keyname, mcoeff);
            cpl_propertylist_set_comment(header, keyname, comment);

            cpl_free(comment);
            cpl_free(unit);
            cpl_free(label);
        }

        /*
         * These parameters are now useless, I set them to zero.
         */

        keyname = "ESO QC MOS REFWAVE MEAN";

        cpl_propertylist_update_double(header, keyname, 0.0);

        keyname = "ESO QC MOS REFWAVE RMS";

        cpl_propertylist_update_double(header, keyname, 0.0);

        if (restab) {

            /*
             * About spectral resolution:
             */

            keyname = "ESO QC MOS RESOLUTION1 LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaRed);
            cpl_propertylist_set_comment(header, keyname, 
                   "Line used in spectral resolution determination (Angstrom)");

            keyname = "ESO QC MOS RESOLUTION1";

            cpl_table_and_selected_double(restab, "wavelength", 
                                          CPL_GREATER_THAN, lambdaRed - 1.0);
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
    
            cpl_msg_info(recipe, "Spectral resolution at %.2f: %.2f +/- %.2f",
                         lambdaRed, resol, resol_err);
    
            cpl_propertylist_update_double(header, keyname, resol);
            cpl_propertylist_set_comment(header, keyname, 
                             "Mean spectral resolution at red end of spectrum");

            keyname = "ESO QC MOS RESOLUTION1 RMS";

            cpl_propertylist_update_double(header, keyname, resol_err);
            cpl_propertylist_set_comment(header, keyname, 
                             "Error on mean spectral resolution");

            keyname = "ESO QC MOS RESOLUTION2 LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaYel);
            cpl_propertylist_set_comment(header, keyname, 
                  "Line used in spectral resolution determination (Angstrom)");

            keyname = "ESO QC MOS RESOLUTION2";

            cpl_table_and_selected_double(restab, "wavelength",
                                          CPL_GREATER_THAN, lambdaYel - 1.0);
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
    
            cpl_msg_info(recipe, "Spectral resolution at %.2f: %.2f +/- %.2f",
                         lambdaYel, resol, resol_err);

            cpl_propertylist_update_double(header, keyname, resol);
            cpl_propertylist_set_comment(header, keyname, 
                             "Mean spectral resolution at center of spectrum");

            keyname = "ESO QC MOS RESOLUTION2 RMS";

            cpl_propertylist_update_double(header, keyname, resol_err);
            cpl_propertylist_set_comment(header, keyname, 
                             "Error on mean spectral resolution");

            keyname = "ESO QC MOS RESOLUTION3 LAMBDA";

            cpl_propertylist_update_double(header, keyname, lambdaBlu);
            cpl_propertylist_set_comment(header, keyname, 
                   "Line used in spectral resolution determination (Angstrom)");

            keyname = "ESO QC MOS RESOLUTION3";

            cpl_table_and_selected_double(restab, "wavelength",
                                          CPL_GREATER_THAN, lambdaBlu - 1.0);
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

            cpl_msg_info(recipe, "Spectral resolution at %.2f: %.2f +/- %.2f",
                         lambdaBlu, resol, resol_err);

            cpl_propertylist_update_double(header, keyname, resol);
            cpl_propertylist_set_comment(header, keyname, 
                   "Mean spectral resolution at blue end of spectrum");

            keyname = "ESO QC MOS RESOLUTION3 RMS";

            cpl_propertylist_update_double(header, keyname, resol_err);
            cpl_propertylist_set_comment(header, keyname, 
                                         "Error on mean spectral resolution");

            keyname = "ESO QC MOS IDS RMS";

            cpl_propertylist_update_double(header, keyname, mean_rms);
            cpl_propertylist_set_comment(header, keyname, 
                               "Mean accuracy of dispersion solution (pixel)");
        }
    }

    /* Read grism configuration */
    std::auto_ptr<mosca::grism_config> grism_cfg = 
            vimos_grism_config_from_table(grism_table);
    cpl_table_delete(grism_table); grism_table = NULL;

    /* Get the detected slit locations */
    //TODO: This is done on rotated images. Change the things to work 
    //on mosca::images with Y_AXIS
    int size_spec = cpl_image_get_size_x(spectra);
    vimos::detected_slits det_slits = 
        vimos::detected_slits_from_tables(slits, polytraces, size_spec);
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
     * Here we use the new method and is the one saved as a product. The other 
     * is not yet deleted in case it is used for something else.
     */
    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Perform flat field combination...");
    
    std::auto_ptr<mosca::image> master_flat_d;
    std::auto_ptr<mosca::fiera_config> ccd_config;
    master_flat_d = vimos_calib_flat_mos_create_master_flat(calib_slits,
            wave_cal, *grism_cfg, master_bias, master_bias_err,
            master_bias_header, frameset, flat_tag);
    if(master_flat_d.get() == 0 || cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_calib_exit("Cannot combine flat frames");

    /*
     * Perform Flat field normalisation
     */
    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Perform flat field normalisation...");
    cpl_msg_indent_more();

    std::auto_ptr<mosca::image> norm_flat;
    norm_flat.reset(new mosca::image(cpl_image_cast(master_flat_d->get_cpl_image(), 
                                                    CPL_TYPE_FLOAT),
                                     cpl_image_cast(master_flat_d->get_cpl_image_err(), 
                                                    CPL_TYPE_FLOAT), true));

    vimos::flat_normaliser normaliser;
    normaliser.mos_normalise(*norm_flat, wave_cal,
                             coordinate, calib_slits, slits, polytraces,
                             startwavelength, endwavelength,
                             dispersion, sradius, dradius,
                             spa_polyorder, disp_nknots, fit_threshold);
    if( cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_calib_exit("Cannot normalize master flat");
    
    //Saving the flats
    cpl_propertylist * qc_list = vimos_calib_flat_qc(*master_flat_d,
             slits, nx, ny, nflats, *grism_cfg, alltime, xwidth, ywidth, 
             gain, scale);
    vimos_calib_flats_save(master_flat_d,  norm_flat, frameset, flat_tag, 
                           master_screen_flat_tag, master_norm_flat_tag, 
                           parlist, qc_list);


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
            slit_widths.push_back(cpl_table_get_double(slits, "ywidth", i_slit, &null));
        else
            slit_widths.push_back(1);

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
    cpl_propertylist_append_string(sed_header, CPL_DFS_PRO_TYPE, "REDUCED");
    cpl_propertylist_append_string(sed_header, CPL_DFS_PRO_CATG, flat_sed_tag);
    const cpl_frame * ref_flat_frame = cpl_frameset_find_const(frameset, flat_tag); 
    cpl_dfs_save_image(frameset, NULL, parlist, frameset, 
                       ref_flat_frame, wave_profiles, CPL_BPP_IEEE_FLOAT,
                       recipe, sed_header, NULL, PACKAGE "/" PACKAGE_VERSION,  
                       prof_filename.c_str());

    cpl_image_delete(master_flat); master_flat = NULL;

    /* Final saving of products */
    
    cpl_free(grism); grism = NULL;
    instrume = NULL;
    cpl_table_delete(restab); restab = NULL;
    cpl_image_delete(spectra); spectra = NULL;

    if (dfs_save_image(frameset, rectified, reduced_lamp_tag, header,
                       parlist, recipe, version))
        vimos_calib_exit(NULL);

    cpl_image_delete(rectified); rectified = NULL;

    {
        save_header = dfs_load_header(frameset, arc_tag, 0);

        cpl_propertylist_update_double(save_header, "CRPIX2", 1.0);
        cpl_propertylist_update_double(save_header, "CRVAL2", 1.0);
        /* cpl_propertylist_update_double(save_header, "CDELT2", 1.0); */
        cpl_propertylist_update_double(save_header, "CD1_1", 1.0);
        cpl_propertylist_update_double(save_header, "CD1_2", 0.0);
        cpl_propertylist_update_double(save_header, "CD2_1", 0.0);
        cpl_propertylist_update_double(save_header, "CD2_2", 1.0);
        cpl_propertylist_update_string(save_header, "CTYPE1", "LINEAR");
        cpl_propertylist_update_string(save_header, "CTYPE2", "PIXEL");

        if (dfs_save_image(frameset, residual, disp_residuals_tag, save_header,
                           parlist, recipe, version))
            vimos_calib_exit(NULL);

        cpl_image_delete(residual); residual = NULL;
        cpl_propertylist_delete(save_header); save_header = NULL;
    }

    if (!treat_as_lss) {

        /*
         * Wavemap was already produced if treat_as_lss = true
         */

        wavemap = mos_map_wavelengths(coordinate, rainbow, slits, 
                                      polytraces, reference, 
                                      startwavelength, endwavelength, 
                                      dispersion);
    }

    cpl_image_delete(rainbow); rainbow = NULL;

    save_header = dfs_load_header(frameset, arc_tag, 0);
    vimos_preoverscan ps_scan;
    ps_scan.fix_wcs_trimm(save_header);

    cpl_image_turn(wavemap, rotate_back);
    if (dfs_save_image(frameset, wavemap, wavelength_map_tag, save_header,
                       parlist, recipe, version))
        vimos_calib_exit(NULL);

    cpl_image_delete(wavemap); wavemap = NULL;

    cpl_image_turn(coordinate, rotate_back);
    if (dfs_save_image(frameset, coordinate, spatial_map_tag, save_header,
                       parlist, recipe, version))
        vimos_calib_exit(NULL);

    cpl_image_delete(coordinate); coordinate = NULL;
    cpl_propertylist_delete(save_header); save_header = NULL;

    if (dfs_save_table(frameset, polytraces, curv_coeff_tag, NULL,
                       parlist, recipe, version))
        vimos_calib_exit(NULL);

    cpl_table_delete(polytraces); polytraces = NULL;

    if (dfs_save_table(frameset, idscoeff, disp_coeff_tag, header,
                       parlist, recipe, version))
        vimos_calib_exit(NULL);

    mos_rotate_slits(slits, rotate, ccd_ysize, ccd_xsize);
    if (dfs_save_table(frameset, slits, slit_location_tag, NULL,
                       parlist, recipe, version))
        vimos_calib_exit(NULL);

    cpl_table_delete(slits); slits = NULL;

    if (cpl_error_get_code()) {
        cpl_msg_error(cpl_error_get_where(), "error: %s", cpl_error_get_message());
        vimos_calib_exit(NULL);
    }
    cpl_propertylist_delete(header); header = NULL;
    cpl_table_delete(idscoeff); idscoeff = NULL;

    return 0;
}

std::auto_ptr<mosca::image> vimos_calib_flat_mos_create_master_flat
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

int vimos_calib_flats_save
(std::auto_ptr<mosca::image>& master_flat_d, 
 std::auto_ptr<mosca::image>& norm_flat,
 cpl_frameset * frameset, const char * flat_tag, 
 const char * master_screen_flat_tag, const char * master_norm_flat_tag, 
 cpl_parameterlist * parlist, cpl_propertylist * qc_list)
{
    cpl_propertylist * save_header;
    const char *recipe_name = "vmmoscalib";
    char        version[80];
    snprintf(version, 80, "%s-%s", PACKAGE, PACKAGE_VERSION);
    
    cpl_msg_indent_more();

    size_t nflats = cpl_frameset_count_tags(frameset, flat_tag);
    save_header = dfs_load_header(frameset, flat_tag, 0);
    cpl_propertylist_update_int(save_header, "ESO PRO DATANCOM", nflats);
    cpl_propertylist_append(save_header, qc_list);

    /* Saving regular flat */
    cpl_image_turn(master_flat_d->get_cpl_image(), -1);
    cpl_image_turn(master_flat_d->get_cpl_image_err(), -1);
    dfs_save_image(frameset, master_flat_d->get_cpl_image(), 
                   master_screen_flat_tag, save_header, parlist, 
                   recipe_name, version);
    dfs_save_image_ext(master_flat_d->get_cpl_image_err(), 
                       master_screen_flat_tag, NULL);

    if(cpl_error_get_code() != CPL_ERROR_NONE)
    {
        cpl_propertylist_delete(save_header);
        return -1;
    }

    /* Saving normalised flats */
    if(norm_flat.get() != NULL)
    {
        cpl_image_turn(norm_flat->get_cpl_image(), -1);
        cpl_image_turn(norm_flat->get_cpl_image_err(), -1);
        dfs_save_image(frameset, norm_flat->get_cpl_image(), 
                       master_norm_flat_tag, save_header, parlist, 
                       recipe_name, version);
        dfs_save_image_ext(norm_flat->get_cpl_image_err(), 
                           master_norm_flat_tag, NULL);
        if(cpl_error_get_code() != CPL_ERROR_NONE)
        {
            cpl_propertylist_delete(save_header);
            return -1;
        }
    }


    cpl_propertylist_delete(save_header);

    cpl_msg_indent_less();

    return 0;
}


cpl_propertylist * vimos_calib_flat_qc(mosca::image& master_flat, 
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
