/* $Id: vimos_science_impl.c,v 1.24 2013-04-24 14:05:53 cgarcia Exp $
 *
 * This file is part of the VIMOS Data Reduction Pipeline
 * Copyright (C) 2006-2010 European Southern Observatory
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
 * $Date: 2013-04-24 14:05:53 $
 * $Revision: 1.24 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vimos_science_impl.h>

#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>
#include <math.h>
#include <cpl.h>
#include <moses.h>
#include <fors_tools.h>
#include <vimos_dfs.h>

#include "vimos_overscan.h"
#include "vimos_detmodel.h"
#include "vimos_response.h"
#include "vimos_detected_slits.h"
#include "vimos_flat_normalise.h"

static int vimos_science_response_fill_ignore
(const cpl_table * flux_table, const cpl_table * telluric_table,
 const std::string& grism_name,
 const std::string& resp_ignore_mode, const std::string& resp_ignore_lines,
 std::vector<double>& resp_ignore_lines_list,
 std::vector<std::pair<double, double> >& resp_ignore_ranges_list);

static bool vimos_science_photcal_apply_flat_corr
(cpl_frameset * frameset, const char * specphot_tag, 
 const char * master_specphot_tag, const char * flat_sed_tag,
 int* fluxcal, bool response_apply_flat_corr, int standard);

static bool vimos_science_response_apply_flat_corr
(cpl_frameset * frameset, const char * flat_sed_tag,
 std::string& resp_use_flat_sed, cpl_table * grism_table, int standard);

static void vimos_shift_stdstar(cpl_image * stdstar_spectra, 
                                cpl_image * stdstar_spectra_error,
                                double dispersion, 
                                double resp_shift);

#define vimos_science_exit(message)           \
{                                             \
if ((const char *)message != NULL) cpl_msg_error(recipe, message);  \
cpl_free(exptime);                            \
cpl_image_delete(dummy);                      \
cpl_image_delete(mapped);                     \
cpl_image_delete(mapped_sky);                 \
cpl_image_delete(mapped_cleaned);             \
cpl_image_delete(skylocalmap);                \
cpl_image_delete(skymap);                     \
cpl_image_delete(smapped);                    \
cpl_table_delete(offsets);                    \
cpl_table_delete(sky);                        \
if(master_bias != NULL) cpl_image_delete(master_bias);              \
if(master_bias_err != NULL) cpl_image_delete(master_bias_err);      \
cpl_image_delete(spectra);                    \
cpl_image_delete(coordinate);                 \
cpl_image_delete(norm_flat);                  \
cpl_image_delete(rainbow);                    \
cpl_image_delete(rectified);                  \
cpl_image_delete(wavemap);                    \
cpl_propertylist_delete(header);              \
cpl_propertylist_delete(save_header);         \
cpl_table_delete(grism_table);                \
cpl_table_delete(idscoeff);                   \
cpl_table_delete(maskslits);                  \
cpl_table_delete(overscans);                  \
cpl_table_delete(polytraces);                 \
cpl_table_delete(slits);                      \
cpl_table_delete(wavelengths);                \
cpl_vector_delete(lines);                     \
cpl_msg_indent_less();                        \
return -1;                                    \
}

/**
 * @addtogroup vimos_science_impl
 */

/**@{*/

int vimos_science_impl(cpl_frameset *frameset, cpl_parameterlist *parlist, 
                       int section)
{
    const char *recipe = "vmmosscience";


    /*
     * Input parameters
     */

    double      dispersion;
    int         skyalign;
    const char *wcolumn = "WLEN";
    const char *resp_ignore_mode;
    const char *resp_ignore_lines;
    double      resp_shift;
    double      startwavelength;
    double      endwavelength;
    double      reference;
    int         flatfield;
    int         skyglobal;
    int         skylocal;
    int         skymedian;
    int         cosmics;
    int         slit_margin;
    int         ext_radius;
    int         cont_radius;
    int         ext_mode;
    double      detection;
    int         time_normalise;
    int         anyframe;
    int         resp_fit_nknots;
    int         resp_fit_degree;
    int         photometry = 0;

    /*
     * CPL objects
     */

    cpl_imagelist    *all_science;
    cpl_imagelist    *all_science_err;
    cpl_image       **images;

    cpl_image        *master_bias    = NULL;
    cpl_image        *master_bias_err= NULL;
    cpl_image        *norm_flat      = NULL;
    cpl_image        *norm_flat_err  = NULL;
    cpl_image        *spectra        = NULL;
    cpl_image        *spectra_err    = NULL;
    cpl_image        *rectified      = NULL;
    cpl_image        *coordinate     = NULL;
    cpl_image        *rainbow        = NULL;
    cpl_image        *mapped         = NULL;
    cpl_image        *mapped_var     = NULL;
    cpl_image        *mapped_sky     = NULL;
    cpl_image        *mapped_sky_var = NULL;
    cpl_image        *mapped_cleaned = NULL;
    cpl_image        *science_ima_nosky =NULL;
    cpl_image        *science_ima_nosky_var =NULL;
    cpl_image        *smapped        = NULL;
    cpl_image        *smapped_var    = NULL;
    cpl_image        *smapped_nosky  = NULL;
    cpl_image        *smapped_nosky_var  = NULL;
    cpl_image        *wavemap        = NULL;
    cpl_image        *skymap         = NULL;
    cpl_image        *skylocalmap    = NULL;
    cpl_image        *dummy          = NULL;

    cpl_table        *grism_table    = NULL;
    cpl_table        *overscans      = NULL;
    cpl_table        *wavelengths    = NULL;
    cpl_table        *idscoeff       = NULL;
    cpl_table        *slits          = NULL;
    cpl_table        *maskslits      = NULL;
    cpl_table        *polytraces     = NULL;
    cpl_table        *offsets        = NULL;
    cpl_table        *sky            = NULL;
    cpl_table        *photcal        = NULL;
    cpl_table        *response_interp = NULL;

    cpl_vector       *lines          = NULL;

    cpl_propertylist *header         = NULL;
    cpl_propertylist *wcal_header    = NULL;
    cpl_propertylist *save_header    = NULL;
    cpl_propertylist *qclist         = NULL;

    /*
     * Auxiliary variables
     */

    char        version[80];
    const char *instrume = NULL;
    const char *science_tag = NULL;
    const char *master_norm_flat_tag = NULL;
    const char *disp_coeff_tag = NULL;
    const char *disp_coeff_sky_tag = NULL;
    const char *wavelength_map_sky_tag = NULL;
    const char *curv_coeff_tag = NULL;
    const char *slit_location_tag = NULL;
    const char *reduced_science_tag = NULL;
    const char *reduced_sky_tag = NULL;
    const char *reduced_flux_science_tag = NULL;
    const char *reduced_error_tag = NULL;
    const char *reduced_flux_error_tag = NULL;
    const char *mapped_science_tag = NULL;
    const char *mapped_flux_science_tag = NULL;
    const char *unmapped_science_tag = NULL;
    const char *unmapped_science_err_tag = NULL;
    const char *mapped_science_sky_tag = NULL;
    const char *mapped_sky_tag = NULL;
    const char *unmapped_sky_tag = NULL;
    const char *global_sky_spectrum_tag = NULL;
    const char *flat_sed_tag = NULL;
    const char *object_table_tag = NULL;
    const char *skylines_offsets_tag = NULL;
    const char *specphot_tag = NULL;
    const char *specphot_in_tag = NULL;
    const char *telluric_contamination_tag = NULL;
    int         mos;
    int         nexp, expno = -1;
    int         treat_as_lss = 0;
    int         have_phot = 0;
    int         nslits;
    int         nscience;
    int         quadrant;
    int         curved;
    double     *exptime = NULL;
    double      airmass = -1;
    double      alltime;
    double      mean_rms;
    int         nlines;
    double     *line;
    int         nx, ny;
    double      gain;
    double      ron;
    int         standard;
    int         highres;
    int         rotate = 1;
    int         rotate_back = -1;
    int         i;
    double      wstart;
    double      wstep;
    int         wcount;


    snprintf(version, 80, "%s-%s", PACKAGE, PACKAGE_VERSION);

    cpl_msg_set_indentation(2);


    /* 
     * Get configuration parameters
     */

    cpl_msg_info(recipe, "Recipe %s configuration parameters:", recipe);
    cpl_msg_indent_more();

    if (cpl_frameset_count_tags(frameset, "CONFIG_TABLE") > 1)
        vimos_science_exit("Too many in input: CONFIG_TABLE");

    grism_table = dfs_load_table(frameset, "CONFIG_TABLE", 1);

    skyalign = dfs_get_parameter_int(parlist, 
                    "vimos.vmmosscience.skyalign", NULL);

    if (skyalign > 2)
        vimos_science_exit("Max polynomial degree for sky alignment is 2");

    flatfield = dfs_get_parameter_bool(parlist, "vimos.vmmosscience.flatfield", 
                                       NULL);

    skyglobal = dfs_get_parameter_bool(parlist, "vimos.vmmosscience.skyglobal", 
                                       NULL);
    skylocal  = dfs_get_parameter_bool(parlist, "vimos.vmmosscience.skylocal", 
                                       NULL);
    skymedian = dfs_get_parameter_bool(parlist, "vimos.vmmosscience.skymedian", 
                                       NULL);

    if (skylocal && skyglobal)
        vimos_science_exit("Cannot do both local and global sky subtraction");

    if (skylocal && skymedian)
        vimos_science_exit("Cannot do sky subtraction both on extracted "
                           "and non-extracted spectra");

    cosmics = dfs_get_parameter_bool(parlist, 
                                     "vimos.vmmosscience.cosmics", NULL);

    if (cosmics)
        if (!(skyglobal || skylocal))
            vimos_science_exit("Cosmic rays correction requires "
                              "either skylocal=true or skyglobal=true");

    slit_margin = dfs_get_parameter_int(parlist, 
                                        "vimos.vmmosscience.slit_margin",
                                        NULL);
    if (slit_margin < 0)
        vimos_science_exit("Value must be zero or positive");

    ext_radius = dfs_get_parameter_int(parlist, 
                                       "vimos.vmmosscience.ext_radius",
                                       NULL);
    if (ext_radius < 0)
        vimos_science_exit("Value must be zero or positive");

    cont_radius = dfs_get_parameter_int(parlist, 
                                        "vimos.vmmosscience.cont_radius",
                                        NULL);
    if (cont_radius < 0)
        vimos_science_exit("Value must be zero or positive");

    ext_mode = dfs_get_parameter_int(parlist, "vimos.vmmosscience.ext_mode",
                                       NULL);
    if (ext_mode < 0 || ext_mode > 1)
        vimos_science_exit("Invalid object extraction mode");

    resp_fit_nknots = dfs_get_parameter_int(parlist, 
                                           "vimos.vmmosscience.resp_fit_nknots",
                                            NULL);
    if(resp_fit_nknots == -2)  //-2 means read from grism table
        resp_fit_nknots = dfs_get_parameter_int(parlist,
                                           "vimos.vmmosscience.resp_fit_nknots",
                                            grism_table);

    if (resp_fit_nknots >= 0 &&  resp_fit_nknots < 2)
        vimos_science_exit("Invalid instrument response spline knots");

    resp_fit_degree = dfs_get_parameter_int(parlist, 
                                           "vimos.vmmosscience.resp_fit_degree",
                                            NULL);
    if(resp_fit_degree == -2)
        resp_fit_degree = dfs_get_parameter_int(parlist,
                                           "vimos.vmmosscience.resp_fit_degree",
                                            grism_table);

    if (resp_fit_degree >= 0 && resp_fit_degree < 1)
        vimos_science_exit("Invalid instrument response polynomial degree");
    
    if (resp_fit_degree > 0 && resp_fit_nknots > 0 )
        vimos_science_exit("Only spline or polynomial fitting of response allowed, but not both");

    resp_ignore_mode = dfs_get_parameter_string(parlist, 
                    "vimos.vmmosscience.resp_ignore_mode", NULL);

    resp_ignore_lines = dfs_get_parameter_string(parlist, 
                    "vimos.vmmosscience.resp_ignore_points", NULL);

    resp_shift = dfs_get_parameter_double(parlist, 
                    "vimos.vmmosscience.resp_shift", NULL);

    std::string resp_use_flat_sed = dfs_get_parameter_string(parlist, 
                    "vimos.vmmosscience.resp_use_flat_sed", NULL);
    std::transform(resp_use_flat_sed.begin(), resp_use_flat_sed.end(),
                   resp_use_flat_sed.begin(), ::tolower);

    detection = dfs_get_parameter_double(parlist, 
                    "vimos.vmmosscience.detection", NULL);

    if (detection <= 0.0)
        vimos_science_exit("Invalid detection threshold");

    time_normalise = dfs_get_parameter_bool(parlist, 
                             "vimos.vmmosscience.time_normalise", NULL);

    anyframe = dfs_get_parameter_bool(parlist, "vimos.vmmosscience.anyframe", 
                                      NULL);

    if (cpl_error_get_code())
        vimos_science_exit("Failure getting the configuration parameters");

    /* 
     * Check input set-of-frames
     */

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Check input set-of-frames:");
    cpl_msg_indent_more();

    {
        cpl_frameset *subframeset = cpl_frameset_duplicate(frameset);
        cpl_frameset_erase(subframeset, "EXTINCT_TABLE");
        cpl_frameset_erase(subframeset, "STD_FLUX_TABLE");
        cpl_frameset_erase(subframeset, "CONFIG_TABLE");

        if (!dfs_equal_keyword(subframeset, "ESO OCS CON QUAD"))
            vimos_science_exit("Input frames are not from the same quadrant");

        cpl_frameset_delete(subframeset);
    }

    mos = cpl_frameset_count_tags(frameset, "MOS_SCIENCE");
    standard = 0;

    if (mos == 0) {
        mos = cpl_frameset_count_tags(frameset, "MOS_STANDARD");
        standard = 1;
    }

    if (mos == 0)
        vimos_science_exit("Missing input scientific frame");

    nscience = mos;

    if (standard) {
        science_tag              = "MOS_STANDARD";
        reduced_science_tag      = "MOS_STANDARD_REDUCED";
        reduced_flux_science_tag = "MOS_STANDARD_FLUX_REDUCED";
        unmapped_science_tag     = "MOS_UNMAPPED_STANDARD";
        unmapped_science_err_tag = "MOS_UNMAPPED_STANDARD_ERR";
        mapped_science_tag       = "MOS_STANDARD_EXTRACTED";
        mapped_flux_science_tag  = "MOS_STANDARD_FLUX_EXTRACTED";
        mapped_science_sky_tag   = "MOS_STANDARD_SKY_EXTRACTED";
        mapped_sky_tag           = "MOS_STANDARD_SKY";
        reduced_sky_tag          = "MOS_STD_SKY_REDUCED";
        reduced_error_tag        = "MOS_STD_ERROR_REDUCED";
        reduced_flux_error_tag   = "MOS_STD_ERROR_FLUX_REDUCED";
        disp_coeff_sky_tag       = "MOS_STD_DISP_COEFF_SKY";
        wavelength_map_sky_tag   = "MOS_STD_WAVELENGTH_MAP_SKY";
        skylines_offsets_tag     = "MOS_STD_SKYLINES_OFFSETS_SLIT";
        unmapped_sky_tag         = "MOS_STD_UNMAPPED_SKY";
        global_sky_spectrum_tag  = "MOS_STD_GLOBAL_SKY_SPECTRUM";
        object_table_tag         = "OBJECT_STD_TABLE";
    }
    else {
        science_tag              = "MOS_SCIENCE";
        reduced_science_tag      = "MOS_SCIENCE_REDUCED";
        reduced_flux_science_tag = "MOS_SCIENCE_FLUX_REDUCED";
        unmapped_science_tag     = "MOS_UNMAPPED_SCIENCE";
        unmapped_science_err_tag = "MOS_UNMAPPED_SCIENCE_ERR";
        mapped_science_tag       = "MOS_SCIENCE_EXTRACTED";
        mapped_flux_science_tag  = "MOS_SCIENCE_FLUX_EXTRACTED";
        mapped_science_sky_tag   = "MOS_SCIENCE_SKY_EXTRACTED";
        mapped_sky_tag           = "MOS_SCIENCE_SKY";
        reduced_sky_tag          = "MOS_SCI_SKY_REDUCED";
        reduced_error_tag        = "MOS_SCI_ERROR_REDUCED";
        reduced_flux_error_tag   = "MOS_SCI_ERROR_FLUX_REDUCED";
        disp_coeff_sky_tag       = "MOS_SCI_DISP_COEFF_SKY";
        wavelength_map_sky_tag   = "MOS_SCI_WAVELENGTH_MAP_SKY";
        skylines_offsets_tag     = "MOS_SCI_SKYLINES_OFFSETS_SLIT";
        unmapped_sky_tag         = "MOS_SCI_UNMAPPED_SKY";
        global_sky_spectrum_tag  = "MOS_SCI_GLOBAL_SKY_SPECTRUM";
        object_table_tag         = "OBJECT_SCI_TABLE";
    }

    flat_sed_tag               = "MOS_FLAT_SED";
    master_norm_flat_tag       = "MOS_MASTER_SCREEN_FLAT";
    disp_coeff_tag             = "MOS_DISP_COEFF";
    curv_coeff_tag             = "MOS_CURV_COEFF";
    slit_location_tag          = "MOS_SLIT_LOCATION";
    specphot_tag               = "MOS_SPECPHOT_TABLE";
    specphot_in_tag            = "MOS_SPECPHOT_TABLE";
    telluric_contamination_tag = "TELLURIC_CONTAMINATION";

    if (cpl_frameset_count_tags(frameset, "MASTER_BIAS") == 0)
        vimos_science_exit("Missing required input: MASTER_BIAS");

    if (cpl_frameset_count_tags(frameset, "MASTER_BIAS") > 1)
        vimos_science_exit("Too many in input: MASTER_BIAS");

    if (skyalign >= 0)
        if (cpl_frameset_count_tags(frameset, "SKY_LINE_CATALOG") > 1)
            vimos_science_exit("Too many in input: SKY_LINE_CATALOG");

    if (cpl_frameset_count_tags(frameset, telluric_contamination_tag) > 1)
        vimos_science_exit("Too many in input: TELLURIC_CONTAMINATION");

    if (cpl_frameset_count_tags(frameset, flat_sed_tag) > 1)
        vimos_science_exit("Too many in input: FLAT_SED_*");

    if (cpl_frameset_count_tags(frameset, disp_coeff_tag) == 0) {
        cpl_msg_error(recipe, "Missing required input: %s", disp_coeff_tag);
        vimos_science_exit(NULL);
    }

    if (cpl_frameset_count_tags(frameset, disp_coeff_tag) > 1) {
        cpl_msg_error(recipe, "Too many in input: %s", disp_coeff_tag);
        vimos_science_exit(NULL);
    }

    if (cpl_frameset_count_tags(frameset, master_norm_flat_tag) > 1) {
        if (flatfield) {
            cpl_msg_error(recipe, "Too many in input: %s", 
                          master_norm_flat_tag);
            vimos_science_exit(NULL);
        }
        else {
            cpl_msg_warning(recipe, "%s in input are ignored, "
                            "since flat field correction was not requested", 
                            master_norm_flat_tag);
        }
    }

    if (cpl_frameset_count_tags(frameset, master_norm_flat_tag) == 1) {
        if (!flatfield) {
            cpl_msg_warning(recipe, "%s in input is ignored, "
                            "since flat field correction was not requested", 
                            master_norm_flat_tag);
        }
    }

    if (cpl_frameset_count_tags(frameset, master_norm_flat_tag) == 0) {
        if (flatfield) {
            cpl_msg_error(recipe, "Flat field correction was requested, "
                          "but no %s are found in input",
                          master_norm_flat_tag);
            vimos_science_exit("Not flat field found");
        }
    }

    if (standard) {

        if (cpl_frameset_count_tags(frameset, "EXTINCT_TABLE") == 0) {
            cpl_msg_warning(recipe, "An EXTINCT_TABLE was not found in input: "
                            "instrument response curve will not be produced.");
            standard = 0;
        }

        if (cpl_frameset_count_tags(frameset, "EXTINCT_TABLE") > 1)
            vimos_science_exit("Too many in input: EXTINCT_TABLE");

        if (cpl_frameset_count_tags(frameset, "STD_FLUX_TABLE") == 0) {
            cpl_msg_warning(recipe, "A STD_FLUX_TABLE was not found in input: "
                            "instrument response curve will not be produced.");
            standard = 0;
        }

        if (cpl_frameset_count_tags(frameset, "STD_FLUX_TABLE") > 1)
            vimos_science_exit("Too many in input: STD_FLUX_TABLE");

        if (!dfs_equal_keyword(frameset, "ESO OBS TARG NAME")) {
            cpl_msg_warning(recipe, "The target name of observation does not "
                            "match the standard star catalog: "
                            "instrument response curve will not be produced.");
            standard = 0;
        }
        //If we still "have" a standard, perform also photometry
        if(standard)
           photometry = 1;
    }

    have_phot = cpl_frameset_count_tags(frameset, specphot_in_tag);
    if (have_phot == 0) {
        specphot_in_tag = "MOS_MASTER_RESPONSE";
    }
    have_phot = cpl_frameset_count_tags(frameset, specphot_in_tag);

    if (!standard) {
        if (have_phot == 0) 
            photometry = 0;
        else if (have_phot > 1) 
        {
            vimos_science_exit("Too many in input: SPECPHOT_TABLE");
        }
        else 
            photometry = 1;
    }

    if (photometry) {

        if (cpl_frameset_count_tags(frameset, "EXTINCT_TABLE") == 0) {
            vimos_science_exit("An EXTINCT_TABLE was not found in input: "
                    "the requested photometric calibrated spectra "
                    "cannot be produced.");
        }

        if (cpl_frameset_count_tags(frameset, "EXTINCT_TABLE") > 1)
            vimos_science_exit("Too many in input: EXTINCT_TABLE");

    }

    cpl_msg_indent_less();

    //Getting wavelengths parameters from the headers
    wcal_header = dfs_load_header(frameset, disp_coeff_tag, 0);

    dispersion = cpl_propertylist_get_double(wcal_header, "ESO PRO WLEN INC"); 

    reference = cpl_propertylist_get_double(wcal_header, "ESO PRO WLEN CEN");

    startwavelength = cpl_propertylist_get_double(wcal_header, "ESO PRO WLEN START");

    endwavelength = cpl_propertylist_get_double(wcal_header, "ESO PRO WLEN END");

    if (cpl_error_get_code() != CPL_ERROR_NONE) 
    {
        cpl_msg_error(cpl_func, "Cannot get wavelength calibration parameters "
                "from dispersion coefficients table");
        return -1;
    }
    
    cpl_msg_info(cpl_func, " Using start wavelength from calibrations: %f", startwavelength);
    cpl_msg_info(cpl_func, " Using end wavelength from calibrations: %f", endwavelength);
    cpl_msg_info(cpl_func, " Using reference wavelength from calibrations: %f", reference);
    cpl_msg_info(cpl_func, " Using dispersion from calibrations: %f", dispersion);

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
        vimos_science_exit("Cannot load master bias. Maybe error extension is missing");
    cpl_image * master_bias_var = cpl_image_power_create(master_bias_err, 2);
    cpl_image_delete(master_bias_err); master_bias_err = NULL;

    /*
     * Loading input science
     */

    exptime = (double *)cpl_calloc(nscience, sizeof(double));

    //TODO: it looks like the code after the if is never executed, since the 
    //sof always comes with one single science. The case of more science is 
    //treated specially by the vimos_science_mult() function in vimos_science.cc
    //which in turn calls this function but then only for each individual frame.
    //The problem is that there is then a lot of code duplication.
    //Note that the average in vimos_science.cc takes into account the
    //shifts, whereas the code below simply performs a median.
    //Also, the average here is of the raw frames, whereas in vimos_science.cc
    //it makes an average of the unmmaped frames, which are already
    //half-processed.
    if (nscience > 1) {

        cpl_msg_info(recipe, "Load %d scientific frames and median them...",
                     nscience);
        cpl_msg_indent_more();

        all_science = cpl_imagelist_new();
        all_science_err = cpl_imagelist_new();

        header = dfs_load_header(frameset, science_tag, 0);

        if (header == NULL)
            vimos_science_exit("Cannot load scientific frame header");

        alltime = exptime[0] = cpl_propertylist_get_double(header, "EXPTIME");

        if (cpl_error_get_code() != CPL_ERROR_NONE)
            vimos_science_exit("Missing keyword EXPTIME in scientific "
                              "frame header");

        if (standard) {
            airmass = fors_get_airmass(header);
            if (airmass < 0.0)
                vimos_science_exit("Missing airmass information in "
                                   "scientific frame header");
        }

        instrume = cpl_propertylist_get_string(header, "INSTRUME");
        if (instrume == NULL)
            vimos_science_exit("Missing keyword INSTRUME "
                               "in scientific frame header");

        cpl_propertylist_delete(header); header = NULL;

        cpl_msg_info(recipe, "Scientific frame 1 exposure time: %.2f s", 
                     exptime[0]);

        for (i = 1; i < nscience; i++) {

            header = dfs_load_header(frameset, NULL, 0);

            if (header == NULL)
                vimos_science_exit("Cannot load scientific frame header");
    
            exptime[i] = cpl_propertylist_get_double(header, "EXPTIME");

            alltime += exptime[i];
    
            if (cpl_error_get_code() != CPL_ERROR_NONE)
                vimos_science_exit("Missing keyword EXPTIME in scientific "
                                  "frame header");
    
            cpl_propertylist_delete(header); header = NULL;

            cpl_msg_info(recipe, "Scientific frame %d exposure time: %.2f s", 
                         i+1, exptime[i]);
        }

        cpl_image * spectra_raw =
                dfs_load_image(frameset, science_tag, CPL_TYPE_FLOAT, 0, 0);

        cpl_propertylist * science_header = 
                dfs_load_header(frameset, science_tag, 0);
        
        /* Create variances map */
        cpl_image * science_var = 
            vimos_image_variance_from_detmodel(spectra_raw, 
                    science_header, master_bias_header);

        //Subtract prescan
        cpl_image * science_notrim = 
             vimos_subtract_overscan(spectra_raw, science_var, science_header);

        //Trimm science
        spectra = vimos_trimm_preoverscan(science_notrim, science_header);
        cpl_image * spectra_var = vimos_trimm_preoverscan(science_var, science_header);
        
        //Subtract master bias
        cpl_image_subtract(spectra, master_bias);
        cpl_image_add(spectra_var, master_bias_var);

        //Getting the science error
        spectra_err = spectra_var;
        cpl_image_power(spectra_err, 0.5);

        if (spectra == NULL)
            vimos_science_exit("Cannot load scientific frame");

        cpl_image_divide_scalar(spectra, exptime[0]);
        cpl_image_divide_scalar(spectra_err, exptime[0]);
        cpl_imagelist_set(all_science, spectra, 0);
        cpl_imagelist_set(all_science_err, spectra_err, 0);
        cpl_image_delete(spectra_raw); spectra_raw = NULL;
        cpl_image_delete(science_var); science_var = NULL;
        cpl_image_delete(science_notrim); science_notrim = NULL;

        for (i = 1; i < nscience; i++) {

            spectra_raw =
                    dfs_load_image(frameset, science_tag, CPL_TYPE_FLOAT, 0, 0);
            science_header = 
                    dfs_load_header(frameset, science_tag, 0);
            
            /* Create variances map */
            science_var = 
                vimos_image_variance_from_detmodel(spectra_raw, 
                        science_header, master_bias_header);
            
            //Subtract prescan
            science_notrim = 
                 vimos_subtract_overscan(spectra_raw, science_var, science_header);
            
            //Trimm science
            spectra = vimos_trimm_preoverscan(science_notrim, science_header);
            spectra_var = vimos_trimm_preoverscan(science_var, science_header);

            //Subtract master bias
            cpl_image_subtract(spectra, master_bias);
            cpl_image_add(spectra_var, master_bias_var);

            //Getting the science error
            spectra_err = spectra_var;
            cpl_image_power(spectra_err, 0.5);

            if (spectra) {
                cpl_image_divide_scalar(spectra, exptime[i]);
                cpl_image_divide_scalar(spectra_err, exptime[0]);
                cpl_imagelist_set(all_science, spectra, i);
                cpl_imagelist_set(all_science_err, spectra_err, i);
            }
            else
                vimos_science_exit("Cannot load scientific frame");

            cpl_image_delete(spectra_raw);
            cpl_image_delete(science_var);
            cpl_image_delete(science_notrim);
        }

        spectra = cpl_imagelist_collapse_median_create(all_science);
        spectra_err = cpl_imagelist_collapse_median_create(all_science_err);
        cpl_image_multiply_scalar(spectra, alltime);
        cpl_image_multiply_scalar(spectra_err, alltime);

        cpl_imagelist_delete(all_science);
        cpl_imagelist_delete(all_science_err);
    }
    else {
        cpl_msg_info(recipe, "Load scientific exposure...");
        cpl_msg_indent_more();

        header = dfs_load_header(frameset, science_tag, 0);

        if (header == NULL)
            vimos_science_exit("Cannot load scientific frame header");

        alltime = exptime[0] = cpl_propertylist_get_double(header, "EXPTIME");

        if (cpl_error_get_code() != CPL_ERROR_NONE)
            vimos_science_exit("Missing keyword EXPTIME in scientific "
                              "frame header");

        if (standard || photometry) {
            airmass = fors_get_airmass(header);
            if (airmass < 0.0)
                vimos_science_exit("Missing airmass information in "
                                   "scientific frame header");
        }

        instrume = cpl_propertylist_get_string(header, "INSTRUME");
        if (instrume == NULL)
            vimos_science_exit("Missing keyword INSTRUME "
                               "in scientific frame header");

        cpl_propertylist_delete(header); header = NULL;

        cpl_msg_info(recipe, "Scientific frame exposure time: %.2f s", 
                     exptime[0]);

        cpl_image * spectra_raw =
                dfs_load_image(frameset, science_tag, CPL_TYPE_FLOAT, 0, 0);
        cpl_propertylist * science_header = 
                dfs_load_header(frameset, science_tag, 0);

        /* Create variances map */
        cpl_image * science_var = 
            vimos_image_variance_from_detmodel(spectra_raw, 
                    science_header, master_bias_header);
        
        //Subtract prescan
        cpl_image * science_notrim = 
             vimos_subtract_overscan(spectra_raw, science_var, science_header);

        //Trimm science
        spectra = vimos_trimm_preoverscan(science_notrim, science_header);
        cpl_image * spectra_var = vimos_trimm_preoverscan(science_var, science_header);

        //Subtract master bias
        cpl_image_subtract(spectra, master_bias);
        cpl_image_add(spectra_var, master_bias_var);

        //Getting the science error
        spectra_err = spectra_var;
        cpl_image_power(spectra_err, 0.5);
        
        cpl_image_delete(spectra_raw);
        cpl_image_delete(science_var);
        cpl_image_delete(science_notrim);
    }
    cpl_image_delete(master_bias_var); master_bias_var = NULL;
    cpl_image_delete(master_bias); master_bias = NULL;

    if (spectra == NULL)
        vimos_science_exit("Cannot load scientific frame");

    cpl_free(exptime); exptime = NULL;

    cpl_msg_indent_less();


    /*
     * Get some info from header
     */

    header = dfs_load_header(frameset, science_tag, 0);

    if (header == NULL)
        vimos_science_exit("Cannot load scientific frame header");

    quadrant = cpl_propertylist_get_int(header, "ESO OCS CON QUAD");

    std::ostringstream grism_name_key;
    grism_name_key<< "ESO INS GRIS"<<quadrant<<" NAME";
    std::string grism_name = 
            cpl_propertylist_get_string(header, grism_name_key.str().c_str());

    if (standard) {
        if (!anyframe) {

           /*
            *  As a default this recipe wouldn't reduce a frame belonging 
            *  to a quadrant that is not expected to contain a standard 
            *  star. The current template for getting data for photometric 
            *  calibration consists of a sequence of 4 exposures, where 
            *  the same standard star is moved through all quadrants. 
            *  The first exposure has the standard star in quadrant 1, 
            *  the second exposure in quadrant 2, and so on. So a standard 
            *  star will be reduced only if it is part of a template of 
            *  4 exposures (i.e., TPL NEXP == 4), and if its sequence 
            *  number in the template is equal to its quadrant number 
            *  (i.e., TPL EXPNO == OCS CON QUAD).
            */

            nexp = cpl_propertylist_get_int(header, "ESO TPL NEXP");

            if (cpl_error_get_code() != CPL_ERROR_NONE)
                vimos_science_exit("Missing keyword ESO TPL NEXP in "
                                    "scientific frame header");

            if (nexp == 4) {
                expno = cpl_propertylist_get_int(header, "ESO TPL EXPNO");

            if (cpl_error_get_code() != CPL_ERROR_NONE)
                vimos_science_exit("Missing keyword ESO TPL EXPNO in "
                                   "scientific frame header");
            }

            if (quadrant != expno) {
                cpl_msg_warning(recipe, "The MOS_STANDARD frame is not "
                                "expected to contain a standard star: "
                                "instrument response curve will not be "
                                "produced. Set --anyframe=true to skip "
                                "this check.");
                standard = 0;
            }
        }
    }

/*
    if (!dfs_equal_keyword(frameset, key_gris_id))
        vimos_science_exit("Input frames are not from the same grism");

    if (!dfs_equal_keyword(frameset, key_filt_id))
        vimos_science_exit("Input frames are not from the same filter");
*/

    gain = cpl_propertylist_get_double(header, "ESO DET OUT1 CONAD");

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_science_exit("Missing keyword ESO DET OUT1 CONAD in scientific "
                          "frame header");

    cpl_msg_info(recipe, "The gain factor is: %.2f e-/ADU", gain);

    ron = cpl_propertylist_get_double(header, "ESO DET OUT1 RON");

    if (cpl_error_get_code() != CPL_ERROR_NONE)
        vimos_science_exit("Missing keyword ESO DET OUT1 RON in scientific "
                          "frame header");

    ron /= gain;     /* Convert from electrons to ADU */

    cpl_msg_info(recipe, "The read-out-noise is: %.2f ADU", ron);

    /*
     * Check for curved slits
     */

    maskslits = mos_load_slits_vimos(header, 0);
    curved = cpl_table_get_column_max(maskslits, "curved");
    cpl_table_delete(maskslits); maskslits = NULL;

    cpl_msg_info(recipe, "There are curved slits");

    /* 
     * In VIMOS data are "never" LSS. This is why this part is
     * commented out.
     */

//    maskslits = mos_load_slits_vimos(header);
//
//    /*
//     * Check if all slits have the same X offset: in such case,
//     * treat the observation as a long-slit one!
//     */
//
//    mxpos = cpl_table_get_column_median(maskslits, "ytop");
//    xpos = cpl_table_get_data_double(maskslits, "ytop");
//    nslits = cpl_table_get_nrow(maskslits);
//     
//    treat_as_lss = 1;
//    for (i = 0; i < nslits; i++) { 
//        if (fabs(mxpos-xpos[i]) > 0.01) {
//            treat_as_lss = 0;
//            break;
//        }
//    }
//
//    treat_as_lss = 0;  // FIXME Prevent LSS treatment for the moment!!!
//
//    cpl_table_delete(maskslits); maskslits = NULL;
//
//    if (treat_as_lss) {
//        cpl_msg_warning(recipe, "All MOS slits have the same offset: %.2f\n"
//                        "The LSS data reduction strategy is applied!",
//                        mxpos);
//        skylines_offsets_tag   = "MOS_SKYLINES_OFFSETS_LONG";
//    }

    if (treat_as_lss) {
        if (skylocal) {
            if (cosmics)
                vimos_science_exit("Cosmic rays correction for long-slit-like "
                                  "data requires --skyglobal=true");
            skymedian = skylocal;
            skylocal = 0;
        }
    }
    else {
        if (cpl_frameset_count_tags(frameset, curv_coeff_tag) == 0) {
            cpl_msg_error(recipe, "Missing required input: %s", curv_coeff_tag);
            vimos_science_exit(NULL);
        }

        if (cpl_frameset_count_tags(frameset, curv_coeff_tag) > 1) {
            cpl_msg_error(recipe, "Too many in input: %s", curv_coeff_tag);
            vimos_science_exit(NULL);
        }

        if (cpl_frameset_count_tags(frameset, slit_location_tag) == 0) {
            cpl_msg_error(recipe, "Missing required input: %s", 
                          slit_location_tag);
            vimos_science_exit(NULL);
        }
        
        if (cpl_frameset_count_tags(frameset, slit_location_tag) > 1) {
            cpl_msg_error(recipe, "Too many in input: %s", slit_location_tag);
            vimos_science_exit(NULL);
        }
    }

    /* Leave the header on for the next step... */


    /*
     * Rotate frames horizontally with red to the right
     */

    cpl_image_turn(spectra, rotate);
    cpl_image_turn(spectra_err, rotate);

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Load normalised flat field (if present)...");
    cpl_msg_indent_more();

    if (flatfield) {

        cpl_frame  *frame = cpl_frameset_find(frameset, master_norm_flat_tag);
        const char *filename = cpl_frame_get_filename(frame);

        if (cpl_fits_count_extensions(filename) > 1) {
            //Flats have their error extensions after each section
            norm_flat = dfs_load_image(frameset, master_norm_flat_tag, 
                                       CPL_TYPE_FLOAT, 1+2*(section-1), 1);
            norm_flat_err = dfs_load_image(frameset, master_norm_flat_tag, 
                                       CPL_TYPE_FLOAT, 2*section, 1);
        }
        else {
            norm_flat = dfs_load_image(frameset, master_norm_flat_tag, 
                                       CPL_TYPE_FLOAT, 0, 1);
            norm_flat_err = dfs_load_image(frameset, master_norm_flat_tag, 
                                       CPL_TYPE_FLOAT, 1, 1);
        }

        if (norm_flat && norm_flat_err) {
            cpl_image_turn(norm_flat, rotate);
            cpl_image_turn(norm_flat_err, rotate);
            cpl_msg_info(recipe, "Apply flat field correction...");
            if (cpl_image_divide(spectra, norm_flat) != CPL_ERROR_NONE) {
                cpl_msg_error(recipe, "Failure of flat field correction: %s",
                              cpl_error_get_message());
                vimos_science_exit("Cannot perform flat fielding");
            }
            cpl_image * dupl_flat = cpl_image_duplicate(norm_flat);

            cpl_image * spectra_var = cpl_image_power_create(spectra_err, 2);
            cpl_image * norm_flat_var = cpl_image_power_create(norm_flat_err, 2);
            cpl_image_multiply(norm_flat_var, spectra);
            cpl_image_multiply(norm_flat_var, spectra);

            /* Now  dupl->variance = sigma2^2 * data1^2 / data2^2 */

            cpl_image_add(spectra_var, norm_flat_var);

            /* Now  left->variance = sigma1^2 + sigma2^2 * data1^2 / data2^2 */

            cpl_image_divide(spectra_var, dupl_flat);
            cpl_image_divide(spectra_var, dupl_flat);
            
            
            cpl_image_delete(spectra_err);
            spectra_err = cpl_image_power_create(spectra_var, 0.5);

            
            cpl_image_delete(norm_flat); norm_flat = NULL;
            cpl_image_delete(norm_flat_var); norm_flat_var = NULL;
            cpl_image_delete(dupl_flat); dupl_flat = NULL;
            cpl_image_delete(spectra_var); spectra_var = NULL;
        }
        else {
            cpl_msg_error(recipe, "Cannot load input %s for flat field "
                          "correction", master_norm_flat_tag);
            vimos_science_exit(NULL);
        }

    }


    if (skyalign >= 0) {
        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Load input sky line catalog...");
        cpl_msg_indent_more();

        wavelengths = dfs_load_table(frameset, "SKY_LINE_CATALOG", 1);

        if (wavelengths) {

            /*
             * Cast the wavelengths into a (double precision) CPL vector
             */

            nlines = cpl_table_get_nrow(wavelengths);

            if (nlines == 0)
                vimos_science_exit("Empty input sky line catalog");

            if (cpl_table_has_column(wavelengths, wcolumn) != 1) {
                cpl_msg_error(recipe, "Missing column %s in input line "
                              "catalog table", wcolumn);
                vimos_science_exit(NULL);
            }

            line = (double* )cpl_malloc(nlines * sizeof(double));
    
            for (i = 0; i < nlines; i++)
                line[i] = cpl_table_get(wavelengths, wcolumn, i, NULL);

            cpl_table_delete(wavelengths); wavelengths = NULL;

            lines = cpl_vector_wrap(nlines, line);
        }
        else {
            cpl_msg_info(recipe, "No sky line catalog found in input - fine!");
        }
    }


    /*
     * Load the spectral curvature table, or provide a dummy one
     * in case of LSS-like data (single slit with flat curvature)
     */

    if (!treat_as_lss) {
        polytraces = dfs_load_table(frameset, curv_coeff_tag, section);
        if (polytraces == NULL)
            vimos_science_exit("Cannot load spectral curvature table");
    }


    /*
     * Load the slit location table, or provide a dummy one in case
     * of LSS-like data (single slit spanning whole image)
     */

    if (treat_as_lss) {
        slits = cpl_table_new(1);
        cpl_table_new_column(slits, "slit_id", CPL_TYPE_INT);
        cpl_table_set_int(slits, "slit_id", 0, 1);
        cpl_table_new_column(slits, "position", CPL_TYPE_INT);
        cpl_table_set_int(slits, "position", 0, 0);
        cpl_table_new_column(slits, "length", CPL_TYPE_INT);
        cpl_table_set_int(slits, "length", 0, ny);
    }
    else {
        slits = dfs_load_table(frameset, slit_location_tag, section);
        if (slits == NULL)
            vimos_science_exit("Cannot load slits location table");
        mos_rotate_slits(slits, -rotate, nx, ny);
    }


    /*
     * Load the wavelength calibration table
     */

    idscoeff = dfs_load_table(frameset, disp_coeff_tag, section);
    if (idscoeff == NULL)
        vimos_science_exit("Cannot load wavelength calibration table");

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Processing scientific spectra...");
    cpl_msg_indent_more();

    /*
     * This one will also generate the spatial map from the spectral 
     * curvature table (in the case of multislit data)
     */

    if (treat_as_lss) {
        smapped = cpl_image_duplicate(spectra);
        smapped_var = cpl_image_power_create(spectra_err, 2);
    }
    else {
        coordinate = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);

        smapped = mos_spatial_calibration(spectra, slits, polytraces, reference,
                                          startwavelength, endwavelength,
                                          dispersion, 1, coordinate);
        cpl_image * spectra_var = cpl_image_power_create(spectra_err, 2);
        smapped_var = mos_spatial_calibration(spectra_var, slits, polytraces, reference,
                                          startwavelength, endwavelength,
                                          dispersion, 1, coordinate);
        cpl_image_delete(spectra_var);
    }


    /*
     * Generate a rectified wavelength map from the wavelength calibration 
     * table
     */

    rainbow = mos_map_idscoeff(idscoeff, nx, reference, startwavelength, 
                               endwavelength);

    if (dispersion > 1.0)
        highres = 0;
    else
        highres = 1;

    if (skyalign >= 0) {
        if (skyalign) {
            cpl_msg_info(recipe, "Align wavelength solution to reference "
            "skylines applying %d order residual fit...", skyalign);
        }
        else {
            cpl_msg_info(recipe, "Align wavelength solution to reference "
            "skylines applying median offset...");
        }

        if (treat_as_lss) {
            offsets = mos_wavelength_align_lss(smapped, reference, 
                                               startwavelength, endwavelength, 
                                               idscoeff, lines, highres, 
                                               skyalign, rainbow, 4);
        }
        else {
            offsets = mos_wavelength_align(smapped, slits, reference, 
                                           startwavelength, endwavelength, 
                                           idscoeff, lines, highres, skyalign, 
                                           rainbow, 4);
        }

        if (offsets) {
            if (standard)
                cpl_msg_warning(recipe, "Alignment of the wavelength solution "
                                "to reference sky lines may be unreliable in "
                                "this case!");

            if (dfs_save_table(frameset, offsets, skylines_offsets_tag, NULL, 
                               parlist, recipe, version))
                vimos_science_exit(NULL);

            cpl_table_delete(offsets); offsets = NULL;
        }
        else {
            cpl_msg_warning(recipe, "Alignment of the wavelength solution "
                            "to reference sky lines could not be done!");
            skyalign = -1;
        }

    }

    if (treat_as_lss) {
        wavemap = rainbow;
        rainbow = NULL;
    }
    else {
        wavemap = mos_map_wavelengths(coordinate, rainbow, slits, 
                                      polytraces, reference, 
                                      startwavelength, endwavelength,
                                      dispersion);
    }

    cpl_image_delete(rainbow); rainbow = NULL;
    cpl_image_delete(coordinate); coordinate = NULL;

    /*
     * Here the wavelength calibrated slit spectra are created. This frame
     * contains sky_science.
     */

    mapped_sky = mos_wavelength_calibration(smapped, reference,
                                            startwavelength, endwavelength,
                                            dispersion, idscoeff, 1);

    mapped_sky_var = mos_wavelength_calibration(smapped_var, reference,
                                                startwavelength, endwavelength,
                                                dispersion, idscoeff, 1);

    cpl_msg_indent_less();
    cpl_msg_info(recipe, "Check applied wavelength against skylines...");
    cpl_msg_indent_more();

    mean_rms = mos_distortions_rms(mapped_sky, lines, startwavelength,
                                   dispersion, 6, highres);

    cpl_vector_delete(lines); lines = NULL;

    cpl_msg_info(recipe, "Mean residual: %f", mean_rms);

    mean_rms = cpl_table_get_column_mean(idscoeff, "error");

    cpl_msg_info(recipe, "Mean model accuracy: %f pixel (%f A)",
                 mean_rms, mean_rms * dispersion);

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

    if (time_normalise) {
        dummy = cpl_image_divide_scalar_create(mapped_sky, alltime);
        if (dfs_save_image(frameset, dummy, mapped_science_sky_tag, header, 
                           parlist, recipe, version))
            vimos_science_exit(NULL);
        cpl_image * mapped_sky_err =cpl_image_power_create(mapped_sky_var, 0.5);
        cpl_image_divide_scalar(mapped_sky_err, alltime);
        if (dfs_save_image_ext(mapped_sky_err, mapped_science_sky_tag, header))
            vimos_science_exit("Cannot save error extension");
        cpl_image_delete(mapped_sky_err);
        cpl_image_delete(dummy); dummy = NULL;
    }
    else {
        if (dfs_save_image(frameset, mapped_sky, mapped_science_sky_tag, 
                           header, parlist, recipe, version))
            vimos_science_exit(NULL);
        cpl_image * mapped_sky_err =cpl_image_power_create(mapped_sky_var, 0.5);
        if (dfs_save_image_ext(mapped_sky_err, mapped_science_sky_tag, header))
            vimos_science_exit("Cannot save error extension");
        cpl_image_delete(mapped_sky_err);
    }

    if (skyglobal || skylocal) {

        cpl_msg_indent_less();

        if (skyglobal) {
            cpl_msg_info(recipe, "Global sky determination...");
            cpl_msg_indent_more();
            skymap = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
            sky = mos_sky_map_super(spectra, wavemap, dispersion, 
                                    2.0, 50, skymap);
            if (sky) {
                science_ima_nosky = cpl_image_subtract_create(spectra, skymap);
                //TODO: This assumes that the sky computation has no error
                science_ima_nosky_var = cpl_image_power_create(spectra_err, 2);
            }
            else {
                science_ima_nosky = cpl_image_duplicate(spectra);
                //TODO: This assumes that the sky computation has no error
                science_ima_nosky_var = cpl_image_power_create(spectra_err, 2);
                cpl_image_delete(skymap); skymap = NULL;
            }
        }
        else {
            cpl_msg_info(recipe, "Local sky determination...");
            cpl_msg_indent_more();

            if (curved) {
                cpl_image *oneslitwave;

                skymap = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
                nslits = cpl_table_get_nrow(slits);

                for (i = 0; i < nslits; i++) {
                    oneslitwave = cpl_image_duplicate(wavemap);
                    if (mos_slit_wavemap(oneslitwave, i, slits, polytraces, 
                            reference, startwavelength, endwavelength,
                            dispersion) == CPL_ERROR_DATA_NOT_FOUND) {
                        cpl_image_delete(oneslitwave);
                        cpl_error_reset();
                        continue;
                    }

                    sky = mos_sky_map_super(spectra, oneslitwave, dispersion,
                                            1.0, 10, skymap);
                    if (sky == NULL) {
                        cpl_image_delete(oneslitwave);
                        cpl_error_reset();
                        continue;
                    }

                    cpl_image_delete(oneslitwave);
                    cpl_table_delete(sky); sky = NULL;
                }
                //TODO: This assumes that the sky computation has no error
                science_ima_nosky = cpl_image_subtract_create(spectra,
                                                              skymap);
                science_ima_nosky_var = cpl_image_power_create(spectra_err, 2);
            }
            else {
                skymap = mos_subtract_sky(spectra, slits, polytraces, 
                                          reference, startwavelength, 
                                          endwavelength, dispersion);
				science_ima_nosky = cpl_image_duplicate(spectra);
                //TODO: This assumes that the sky computation has no error
                science_ima_nosky_var = cpl_image_power_create(spectra_err, 2);
            }
        }

        if (skymap) {
            if (skyglobal) {
                if (time_normalise)
                    cpl_table_divide_scalar(sky, "sky", alltime);
                if (dfs_save_table(frameset, sky, global_sky_spectrum_tag, 
                                   NULL, parlist, recipe, version))
                    vimos_science_exit(NULL);
    
                cpl_table_delete(sky); sky = NULL;
            }

            save_header = dfs_load_header(frameset, science_tag, 0);

            if (time_normalise)
                cpl_image_divide_scalar(skymap, alltime);
            cpl_image_turn(skymap, rotate_back);
            if (dfs_save_image(frameset, skymap, unmapped_sky_tag,
                               save_header, parlist, recipe, version))
                vimos_science_exit(NULL);

            cpl_image_delete(skymap); skymap = NULL;

            cpl_image_turn(science_ima_nosky, rotate_back);

            if (cosmics) {
                cpl_msg_info(recipe, "Removing cosmic rays...");
                mos_clean_cosmics(science_ima_nosky, gain, -1., -1.);
            }

            if (dfs_save_image(frameset, science_ima_nosky, unmapped_science_tag,
                               save_header, parlist, recipe, version))
                vimos_science_exit(NULL);
            cpl_image * science_ima_nosky_err = 
                    cpl_image_power_create(science_ima_nosky_var, 0.5);
            cpl_image_turn(science_ima_nosky_err, rotate_back);
            if (dfs_save_image(frameset, science_ima_nosky_err, unmapped_science_err_tag,
                               save_header, parlist, recipe, version))
                vimos_science_exit(NULL);
            cpl_image_turn(science_ima_nosky, rotate);
            cpl_image_delete(science_ima_nosky_err); science_ima_nosky_err = NULL;

            cpl_propertylist_delete(save_header); save_header = NULL;
/* Cosmics was here */
            /*
             * The spatially rectified image, that contained the sky,
             * is replaced by a sky-subtracted spatially rectified image:
             */


            if (treat_as_lss) {
                smapped_nosky = cpl_image_duplicate(science_ima_nosky);
                smapped_nosky_var = cpl_image_duplicate(science_ima_nosky_var);
            }
            else {
                smapped_nosky = mos_spatial_calibration(science_ima_nosky, slits, polytraces, 
                                                  reference, startwavelength, 
                                                  endwavelength, dispersion, 
                                                  1, NULL);
                smapped_nosky_var = mos_spatial_calibration(science_ima_nosky_var, slits, polytraces, 
                                                  reference, startwavelength, 
                                                  endwavelength, dispersion, 1, NULL);
            }
        }
        else {
            cpl_msg_warning(recipe, "Sky subtraction failure");
            if (cosmics)
                cpl_msg_warning(recipe, "Cosmic rays removal not performed!");
            cosmics = skylocal = skyglobal = 0;
        }
    }
    else
    {
        smapped_nosky = cpl_image_duplicate(smapped);
        smapped_nosky_var = cpl_image_duplicate(smapped_var);
    }
    cpl_image_delete(science_ima_nosky); science_ima_nosky=NULL;
    cpl_image_delete(science_ima_nosky_var); science_ima_nosky_var=NULL;

    //cpl_image_delete(spectra); spectra = NULL;

    if (skyalign >= 0) {
        save_header = dfs_load_header(frameset, science_tag, 0);
        cpl_image_turn(wavemap, rotate_back);
        if (dfs_save_image(frameset, wavemap, wavelength_map_sky_tag,
                           save_header, parlist, recipe, version))
            vimos_science_exit(NULL);
        cpl_propertylist_delete(save_header); save_header = NULL;
    }

    cpl_image_delete(wavemap); wavemap = NULL;

    mapped = mos_wavelength_calibration(smapped_nosky, reference,
                                        startwavelength, endwavelength,
                                        dispersion, idscoeff, 1);

    mapped_var = mos_wavelength_calibration(smapped_nosky_var, reference,
                                        startwavelength, endwavelength,
                                        dispersion, idscoeff, 1);

    cpl_image_delete(smapped); smapped = NULL;
    cpl_image_delete(smapped_var); smapped_var = NULL;
    cpl_image_delete(smapped_nosky); smapped_nosky = NULL;
    cpl_image_delete(smapped_nosky_var); smapped_nosky_var = NULL;

    if (skymedian) {
        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Local sky determination...");
        cpl_msg_indent_more();
    
        skylocalmap = mos_sky_local_old(mapped, slits);       
        cpl_image_subtract(mapped, skylocalmap);
        cpl_image_delete(skylocalmap); skylocalmap = NULL;
    }

    //Determine whether the response will contain the flat correction
    bool response_apply_flat_corr = 
            vimos_science_response_apply_flat_corr
                (frameset, flat_sed_tag, resp_use_flat_sed, 
                        grism_table, standard);

    //Determine wheter the flux calibration will use the flat correction
    bool photcal_apply_flat_corr = vimos_science_photcal_apply_flat_corr
            (frameset, specphot_tag, specphot_in_tag,
             flat_sed_tag, &photometry, 
             response_apply_flat_corr, standard);

    //TODO: Place this in a better place. The whole recipe needs refactoring...
    cpl_image        * mapped_flat_sed = NULL;
    cpl_propertylist * flat_sed_header = NULL;
    vimos::detected_slits det_slits = 
            vimos::detected_slits_from_tables(slits, polytraces, nx);
    mosca::wavelength_calibration wave_cal(idscoeff, reference);
    if(photcal_apply_flat_corr || response_apply_flat_corr)
    {
        const cpl_frame * flat_sed_frame = 
                cpl_frameset_find_const(frameset, flat_sed_tag);
        const char * filename = cpl_frame_get_filename(flat_sed_frame);
        if (cpl_fits_count_extensions(filename) > 0)
        {
            mapped_flat_sed = dfs_load_image(frameset, flat_sed_tag,  
                                             CPL_TYPE_FLOAT, section, 1);
            flat_sed_header = dfs_load_header(frameset, flat_sed_tag, section);
        }
        else
        {
            mapped_flat_sed = dfs_load_image(frameset, flat_sed_tag, CPL_TYPE_FLOAT, 0, 1);
            flat_sed_header = dfs_load_header(frameset, flat_sed_tag, 0);
        }
    }
    
    if (skyglobal || skymedian || skylocal) {

        skylocalmap = cpl_image_subtract_create(mapped_sky, mapped);

        if (time_normalise) {
            dummy = cpl_image_divide_scalar_create(skylocalmap, alltime);
            if (dfs_save_image(frameset, dummy, mapped_sky_tag, header,
                               parlist, recipe, version))
                vimos_science_exit(NULL);
            cpl_image_delete(dummy); dummy = NULL;
        }
        else {
            if (dfs_save_image(frameset, skylocalmap, mapped_sky_tag, header,
                               parlist, recipe, version))
                vimos_science_exit(NULL);
        }

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Object detection...");
        cpl_msg_indent_more();

        if (cosmics || nscience > 1) {
            dummy = mos_detect_objects(mapped, slits, slit_margin, ext_radius, 
                                       cont_radius, detection);
        }
        else {
            mapped_cleaned = cpl_image_duplicate(mapped);
            mos_clean_cosmics(mapped_cleaned, gain, -1., -1.);
            dummy = mos_detect_objects(mapped_cleaned, slits, slit_margin, 
                                       ext_radius, cont_radius, detection);

            cpl_image_delete(mapped_cleaned); mapped_cleaned = NULL;
        }

        cpl_image_delete(dummy); dummy = NULL;

        mos_rotate_slits(slits, rotate, ny, nx);
        if (dfs_save_table(frameset, slits, object_table_tag, NULL, parlist, 
                           recipe, version))
            vimos_science_exit(NULL);

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Object extraction...");
        cpl_msg_indent_more();

        images = mos_extract_objects(mapped, mapped_var, skylocalmap, slits, 
                                     ext_mode, ron, gain, 1);

        cpl_image_delete(skylocalmap); skylocalmap = NULL;

        if (images) {
            if (standard) {
              if (resp_fit_degree < 0 && resp_fit_nknots < 0 )
                  vimos_science_exit("Either spline or polynomial fitting "
                                     "of response has to be specified");

                cpl_table *ext_table  = NULL;
                cpl_table *flux_table = NULL;
                cpl_table *telluric_table = NULL;
                double flat_sed_norm_factor;

                ext_table = dfs_load_table(frameset, "EXTINCT_TABLE", 1);
                flux_table = dfs_load_table(frameset, "STD_FLUX_TABLE", 1);

                if(cpl_frameset_count_tags(frameset, telluric_contamination_tag))
                    telluric_table = dfs_load_table(frameset, "TELLURIC_CONTAMINATION", 1);

                std::vector<double> resp_ignore_lines_list;
                std::vector<std::pair<double, double> > resp_ignore_ranges_list;
                
                //Fill the points and lines to be ignored by response
                if(vimos_science_response_fill_ignore(flux_table,
                                                     telluric_table,
                                                     grism_name,
                                                     resp_ignore_mode,
                                                     resp_ignore_lines,
                                                     resp_ignore_lines_list,
                                                     resp_ignore_ranges_list))
                    vimos_science_exit("Cannot parse the response ignored lines");

                //Shift the standard star spectrum if needed
                if(std::abs(resp_shift) > 0)
                    vimos_shift_stdstar(images[0], images[1],
                                        dispersion, resp_shift);
                
                //Compute the response
                photcal = vimos_compute_response(images[0], 
                                                 mapped_flat_sed,
                                                 flat_sed_header,
                                                 slits,
                                                 startwavelength,
                                                 dispersion, gain,
                                                 alltime, ext_table,
                                                 airmass, flux_table,
                                                 resp_ignore_lines_list,
                                                 resp_ignore_ranges_list,
                                                 resp_fit_nknots,
                                                 resp_fit_degree,
                                                 response_interp,
                                                 flat_sed_norm_factor,
                                                 det_slits);

                cpl_table_delete(ext_table);
                cpl_table_delete(flux_table);
                if(cpl_frameset_count_tags(frameset, telluric_contamination_tag))
                    cpl_table_delete(telluric_table);

                if (photcal) {

                    float *data;
                    char   keyname[40];

                    qclist = dfs_load_header(frameset, science_tag, 0);

                    if (qclist == NULL)
                        vimos_science_exit("Cannot reload scientific "
                                "frame header");

                    /*
                     * QC1 parameters
                     */

                    cpl_propertylist_update_string(qclist, "ESO QC DID",
                                                   "2.0");                        
                    cpl_propertylist_set_comment(qclist, "ESO QC DID",
                                                 "QC1 dictionary");

                    wstart = 3700.;
                    wstep  = 400.;
                    wcount = 15;

                    dummy = cpl_image_new(wcount, 1, CPL_TYPE_FLOAT);
                    data = cpl_image_get_data_float(dummy);
                    map_table(dummy, wstart, wstep, photcal,
                              "WAVE", "EFFICIENCY");

                    for (i = 0; i < wcount; i++) {
                        sprintf(keyname, "ESO QC MOS EFFICIENCY%d LAMBDA",
                                i + 1);

                        cpl_propertylist_update_double(qclist, keyname,
                                                       wstart + wstep * i);
                        cpl_propertylist_set_comment(qclist, keyname,
                                                     "Wavelength of efficiency evaluation (Angstrom)");

                        sprintf(keyname, "ESO QC MOS EFFICIENCY%d", i + 1);

                        cpl_propertylist_update_double(qclist, keyname,
                                                       data[i]);
                        cpl_propertylist_set_comment(qclist, keyname,
                                                     "Efficiency (e-/photom)");                            
                    }

                    if(mapped_flat_sed != NULL)
                    {
                        cpl_propertylist_update_int(qclist, 
                            "ESO QC RESP FLAT SED_CORR", 1);
                        cpl_propertylist_set_comment(qclist, 
                            "ESO QC RESP FLAT SED_CORR",
                            "Response corrected from flat dispersion profile");

                        cpl_propertylist_update_double(qclist, 
                            "ESO QC RESP FLAT SED_NORM", flat_sed_norm_factor);
                        cpl_propertylist_set_comment(qclist, 
                            "ESO QC RESP FLAT SED_NORM",
                            "Normalisation factor applied to flat sed");
                        
                        cpl_propertylist_update_bool(qclist, 
                            "ESO QC RESP FLAT SED CORR_SLITWID",
                            cpl_propertylist_get_bool(flat_sed_header, 
                                    "ESO QC FLAT SED CORR_SLITWID"));
                        cpl_propertylist_set_comment(qclist, 
                            "ESO QC RESP FLAT SED CORR_SLITWID",
                            "Whether the normalisation factor includes slit width");
                    }
                    else
                    {
                        cpl_propertylist_update_int(qclist, 
                            "ESO QC RESP FLAT SED_CORR", 0);
                        cpl_propertylist_set_comment(qclist, 
                            "ESO QC RESP FLAT SED_CORR",
                            "Response corrected from flat dispersion profile");
                    }
                    
                    cpl_image_delete(dummy); dummy = NULL;


                    dfs_save_table(frameset, photcal, specphot_tag, qclist,
                                       parlist, recipe, version);
                    //TODO: Create a dfs_save_tables to save more than one table in a FITS
                    cpl_table_save(response_interp, NULL, NULL,
                                   "mos_specphot_table.fits", CPL_IO_EXTEND);

                    if(cpl_error_get_code() != CPL_ERROR_NONE)
                        vimos_science_exit(NULL);

                    cpl_propertylist_delete(qclist); qclist = NULL;

                    if (have_phot) {
                        cpl_table_delete(photcal);
                    }
                }
            }

            if (photometry) {
                cpl_image *calibrated;
                cpl_table *ext_table;
                cpl_propertylist * specphot_header;

                ext_table = dfs_load_table(frameset, "EXTINCT_TABLE", 1);

                if (cpl_frameset_count_tags(frameset, specphot_tag) == 0) 
                {
                    response_interp = dfs_load_table(frameset, 
                            specphot_in_tag, 2);
                    specphot_header = dfs_load_header(frameset, specphot_in_tag, 0);
                }
                else 
                {
                    response_interp = dfs_load_table(frameset, specphot_tag, 2);
                    specphot_header = dfs_load_header(frameset, specphot_tag, 0);
                }

                cpl_image * science_images = cpl_image_duplicate(images[0]);

                if(photcal_apply_flat_corr)
                {
                    cpl_msg_info(cpl_func, "Applying flat SED correction");
                    vimos_science_correct_flat_sed(science_images, slits, 
                            mapped_flat_sed, flat_sed_header, specphot_header,
                            det_slits);
                }
                
                calibrated = mos_apply_photometry(science_images, response_interp,
                                                  ext_table, startwavelength,
                                                  dispersion, gain, alltime,
                                                  airmass);
                cpl_table_delete(ext_table);

                cpl_propertylist_update_string(header, "BUNIT",
                                   "10^(-16) erg/(cm^2 s Angstrom)");

                if (dfs_save_image(frameset, calibrated,
                                   reduced_flux_science_tag, header,
                                   parlist, recipe, version)) {
                    cpl_image_delete(calibrated);
                    vimos_science_exit(NULL);
                }

                cpl_image_delete(calibrated); calibrated = NULL;
                cpl_image_delete(science_images); science_images = NULL;
            }

            if (time_normalise) {
                cpl_propertylist_update_string(header, "BUNIT", "ADU/s");
                cpl_image_divide_scalar(images[0], alltime);
            }
            else {
                cpl_propertylist_update_string(header, "BUNIT", "ADU");
            }

            if (dfs_save_image(frameset, images[0], reduced_science_tag, header,
                               parlist, recipe, version))
                vimos_science_exit(NULL);

            if (time_normalise)
                cpl_image_divide_scalar(images[1], alltime);

            if (dfs_save_image(frameset, images[1], reduced_sky_tag, header,
                               parlist, recipe, version))
                vimos_science_exit(NULL);

            if (photometry) {
                cpl_image *calibrated;
                cpl_table *ext_table;
                cpl_propertylist * specphot_header;

                ext_table = dfs_load_table(frameset, "EXTINCT_TABLE", 1);

                if (cpl_frameset_count_tags(frameset, specphot_tag) == 0) 
                {
                    response_interp = dfs_load_table(frameset, 
                            specphot_in_tag, 2);
                    specphot_header = dfs_load_header(frameset, specphot_in_tag, 0);
                }
                else 
                {
                    response_interp = dfs_load_table(frameset, specphot_tag, 2);
                    specphot_header = dfs_load_header(frameset, specphot_tag, 0);
                }

                cpl_image * science_images = cpl_image_duplicate(images[0]);
                cpl_image * science_images_err = cpl_image_duplicate(images[2]);
                if(photcal_apply_flat_corr)
                {
                    //We assume that the mapped flat profile has no associated error
                    vimos_science_correct_flat_sed(science_images, slits, 
                            mapped_flat_sed, flat_sed_header, specphot_header,
                            det_slits);
                    vimos_science_correct_flat_sed(science_images_err, slits, 
                            mapped_flat_sed, flat_sed_header, specphot_header,
                            det_slits);
                }

                calibrated = mos_propagate_photometry_error(science_images, 
                                                  science_images_err, response_interp,
                                                  ext_table, startwavelength,
                                                  dispersion, gain, alltime,
                                                  airmass);

                cpl_table_delete(ext_table);

                cpl_propertylist_update_string(header, "BUNIT",
                                   "10^(-16) erg/(cm^2 s Angstrom)");

                if (dfs_save_image(frameset, calibrated,
                                   reduced_flux_error_tag, header,
                                   parlist, recipe, version)) {
                    cpl_image_delete(calibrated);
                    vimos_science_exit(NULL);
                }

                cpl_image_delete(calibrated); calibrated = NULL;
            }
    
            if (time_normalise) {
                cpl_propertylist_update_string(header, "BUNIT", "ADU/s");
                cpl_image_divide_scalar(images[2], alltime);
            }
            else {
                cpl_propertylist_update_string(header, "BUNIT", "ADU");
            }

            if (dfs_save_image(frameset, images[2], reduced_error_tag, header,
                               parlist, recipe, version))
                vimos_science_exit(NULL);

            cpl_image_delete(images[0]);
            cpl_image_delete(images[1]);
            cpl_image_delete(images[2]);

            cpl_free(images);
        }
        else {
            cpl_msg_warning(recipe, "No objects found: the products "
                            "%s, %s, and %s are empty.", 
                            reduced_science_tag, reduced_sky_tag, 
                            reduced_error_tag);

            if (dfs_save_image_null(frameset, NULL, parlist, reduced_science_tag,
                                    recipe, version))
                vimos_science_exit(NULL);

            if (dfs_save_image_null(frameset, NULL, parlist, reduced_sky_tag, 
                                    recipe, version))
                vimos_science_exit(NULL);

            if (dfs_save_image_null(frameset, NULL, parlist, reduced_error_tag,
                                    recipe, version))
                vimos_science_exit(NULL);
            
            if(standard) //No photometric solution has been created.
                photometry = 0;
        }
    }


    if (skyalign >= 0) {
        if (dfs_save_table(frameset, idscoeff, disp_coeff_sky_tag, NULL, 
                           parlist, recipe, version))
            vimos_science_exit(NULL);
    }

    cpl_table_delete(idscoeff); idscoeff = NULL;

    if (skyglobal || skymedian || skylocal) {

        if (photometry) {
            cpl_image *calibrated;
            cpl_table *ext_table;
            cpl_propertylist * specphot_header;

            ext_table = dfs_load_table(frameset, "EXTINCT_TABLE", 1);

            if (cpl_frameset_count_tags(frameset, specphot_tag) == 0)
            {
                response_interp = dfs_load_table(frameset, 
                        specphot_in_tag, 2);
                specphot_header = dfs_load_header(frameset, specphot_in_tag, 0);
            }
            else 
            {
                response_interp = dfs_load_table(frameset, specphot_tag, 2);
                specphot_header = dfs_load_header(frameset, specphot_tag, 0);
            }

            cpl_image * mapped_images = cpl_image_duplicate(mapped);

            if(photcal_apply_flat_corr)
            {
                cpl_msg_info(cpl_func, "Applying flat SED correction");
                vimos_science_correct_flat_sed_mapped(mapped_images, slits, 
                        mapped_flat_sed, flat_sed_header, specphot_header,
                        det_slits);
            }
           
            calibrated = mos_apply_photometry(mapped_images, response_interp,
                                              ext_table, startwavelength,
                                              dispersion, gain, alltime,
                                              airmass);

            cpl_table_delete(ext_table);

            cpl_propertylist_update_string(header, "BUNIT",
                                   "10^(-16) erg/(cm^2 s Angstrom)");

            if (dfs_save_image(frameset, calibrated,
                               mapped_flux_science_tag, header,
                               parlist, recipe, version)) {
                cpl_image_delete(calibrated);
                vimos_science_exit("Cannot save 2D mapped flux calibrated science");
            }

            cpl_image_delete(mapped_images); calibrated = NULL;
            cpl_image_delete(calibrated); calibrated = NULL;
        }

        if (time_normalise) {
            cpl_propertylist_update_string(header, "BUNIT", "ADU/s");
            cpl_image_divide_scalar(mapped, alltime);
        }
        else {
            cpl_propertylist_update_string(header, "BUNIT", "ADU");
        }

        if (dfs_save_image(frameset, mapped, mapped_science_tag, header, 
                           parlist, recipe, version))
            vimos_science_exit("Cannot save 2D mapped science");

        cpl_image * mapped_err = cpl_image_power_create(mapped_var, 0.5);
        if (time_normalise) 
            cpl_image_divide_scalar(mapped_err, alltime);
        if (dfs_save_image_ext(mapped_err, mapped_science_tag, header))
            vimos_science_exit("Cannot save error extension");
        cpl_image_delete(mapped_err);
    }

    cpl_image_delete(spectra); spectra = NULL;
    cpl_image_delete(spectra_err); spectra_err = NULL;
    cpl_image_delete(mapped); mapped = NULL;
    cpl_image_delete(mapped_var); mapped_var = NULL;
    cpl_image_delete(mapped_sky); mapped_sky = NULL;
    cpl_image_delete(mapped_sky_var); mapped_sky_var = NULL;
    cpl_table_delete(slits); slits = NULL;

    cpl_propertylist_delete(header); header = NULL;
    cpl_table_delete(grism_table); grism_table = NULL;
    cpl_table_delete(polytraces); polytraces = NULL;
    if(photcal_apply_flat_corr || response_apply_flat_corr)
        cpl_image_delete(mapped_flat_sed);


    if (cpl_error_get_code()) {
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        vimos_science_exit(NULL);
    }
    else 
        return 0;
    return 0;
}

static int vimos_science_response_fill_ignore
(const cpl_table * flux_table, const cpl_table * telluric_table,
 const std::string& grism_name,
 const std::string& resp_ignore_mode, const std::string& resp_ignore_lines,
 std::vector<double>& resp_ignore_lines_list,
 std::vector<std::pair<double, double> >& resp_ignore_ranges_list)
{
    //Reading response mode
    bool mask_stellar_absorption = false;
    bool mask_telluric = false;
    bool mask_commandline = false;
    std::string mode(resp_ignore_mode);
    while(mode.length() > 0)
    {
        //Parsing ignore_lines (values are separated by comma)
        std::string::size_type found = mode.find(',');
        std::string mode_str;
        if(found != std::string::npos)
        {
            mode_str = mode.substr(0, found);
            mode = mode.substr(found+1);
        }
        else
        {
            mode_str = mode;
            mode = "";
        }
        if(mode_str == "stellar_absorption")
            mask_stellar_absorption = true;
        if(mode_str == "telluric")
            mask_telluric = true;
        if(mode_str == "command_line")
            mask_commandline = true;
    }
    
    //Adding lines from the standard star table
    if(mask_stellar_absorption)
    {
        if(cpl_table_has_column(flux_table, "STLLR_ABSORP"))
        {
            cpl_size stdtable_size = cpl_table_get_nrow(flux_table);
            for(cpl_size irow = 0; irow < stdtable_size; ++irow)
                if(cpl_table_get_int(flux_table, "STLLR_ABSORP", irow, NULL))
                    resp_ignore_lines_list.push_back
                    (cpl_table_get_float(flux_table, "WAVE", irow, NULL));
        }
        else
            cpl_msg_warning(cpl_func, " Column STLLR_ABSORP not found in std "
                    "star table. Value 'stellar_absorption' in 'resp_ignore_mode' is ignored.");
    }
    
    //Adding regions from the telluric contamination table
    if(mask_telluric && telluric_table != NULL)
    {
        if(cpl_table_has_column(telluric_table, grism_name.c_str()))
        {
            cpl_size telltable_size = cpl_table_get_nrow(telluric_table);
            for(cpl_size irow = 0; irow < telltable_size; ++irow)
                if(cpl_table_get_int(telluric_table, grism_name.c_str(), irow, NULL))
                {
                    double wave_start = 
                            cpl_table_get_float(telluric_table, "START_WAVE", irow, NULL);
                    double wave_end = 
                            cpl_table_get_float(telluric_table, "END_WAVE", irow, NULL);
                    resp_ignore_ranges_list.push_back(std::make_pair(wave_start, wave_end));
                }
        }
    }
    
    //Adding lines and ranges from the command line
    if(mask_commandline)
    {
        std::string lines(resp_ignore_lines);
        while(lines.length() > 0)
        {
            std::string::size_type found = lines.find(',');
            std::string line_str;
            if(found != std::string::npos)
            {
                line_str = lines.substr(0, found);
                lines = lines.substr(found+1);
            }
            else
            {
                line_str = lines;
                lines = "";
            }
            //We have a line or an interval. Let's check which of them
            std::string::size_type found_interval = line_str.find('-');
            if(found_interval != std::string::npos)     //It is an interval
            {
                double wave_start;
                std::istringstream iss1(line_str.substr(0, found_interval));
                //We cannot use simpy iss1 >> std::ws && iss1.eof() because in
                //some STL implementations iss >> std::ws gives a failure
                //if there are no more characters (the standard is not clear)
                //See http://llvm.org/bugs/show_bug.cgi?id=19497
                if ( !(iss1 >> wave_start) || !(iss1.eof() || (iss1 >> std::ws 
&& iss1.eof())) )
                {
                    cpl_msg_error(cpl_func, "Cannot interpret number in resp_ignore_lines");
                    return 1;
                }
                double wave_end;
                std::istringstream iss2(line_str.substr(found_interval + 1));
                if ( !(iss2 >> wave_end) || !(iss2.eof() || (iss2 >> std::ws && iss2.eof())) )
                {
                    cpl_msg_error(cpl_func, "Cannot interpret number in resp_ignore_lines");
                    return 1;
                }
                resp_ignore_ranges_list.push_back(std::make_pair(wave_start, wave_end));
            
            }
            else   //It is a single line
            {
                double wave;
                std::istringstream iss(line_str);
                if ( !(iss >> wave) || !(iss.eof() || (iss >> std::ws && iss.eof())) )
                {
                    cpl_msg_error(cpl_func, "Cannot interpret number in resp_ignore_lines");
                    return 1;
                }
                resp_ignore_lines_list.push_back(wave);
            }
        }
    }
    
    std::sort(resp_ignore_lines_list.begin(), resp_ignore_lines_list.end());  
    std::sort(resp_ignore_ranges_list.begin(), resp_ignore_ranges_list.end());  
    
    cpl_msg_indent_more();
    std::string all_lines;
    if(resp_ignore_lines_list.size() != 0)
    {
        std::ostringstream oss;
        for(size_t i = 0 ; i < resp_ignore_lines_list.size(); i++)
            oss<<resp_ignore_lines_list[i]<<" ";
        cpl_msg_info(cpl_func, "Total list of lines to ignore in response: %s", 
                     oss.str().c_str());
    }
    if(resp_ignore_ranges_list.size() != 0)
    {
        std::ostringstream oss;
        for(size_t i = 0 ; i < resp_ignore_ranges_list.size(); i++)
            oss<<resp_ignore_ranges_list[i].first<<"-"<<resp_ignore_ranges_list[i].second<<" ";
        cpl_msg_info(cpl_func, "Total list of ranges to ignore in response: %s", 
                     oss.str().c_str());
    }
    
    return 0;
}

#define MAX_COLNAME      (80)

static bool vimos_science_response_apply_flat_corr
(cpl_frameset * frameset, const char * flat_sed_tag,
 std::string& resp_use_flat_sed, cpl_table * grism_table, int standard)
{
    if(standard)
    {
        bool have_flat_sed = false;
        if(cpl_frameset_count_tags(frameset, flat_sed_tag) > 0)
            have_flat_sed= true;

        bool requested_in_grism_table = false;
        int null;
        if(cpl_table_get_int(grism_table, "RESP_USE_FLAT_SED", 0, &null))
            requested_in_grism_table = true;

        if(have_flat_sed)
        {
            if(resp_use_flat_sed == "false" || (resp_use_flat_sed == "grism_table" &&
                    !requested_in_grism_table))
            {
                cpl_msg_warning(cpl_func, "Flat SED is part of the input but "
                        "no correction has been requested");
                return false;
            }
            return true;
        }
        else
        {
            if(resp_use_flat_sed == "true" || (resp_use_flat_sed == "grism_table" &&
                    requested_in_grism_table))
                throw std::invalid_argument("Flat SED correction requested "
                        "but MOS_FLAT_SED it is not part of input.");
            return false;
        }
    }
    else
        return false;
    return false;

}

static bool vimos_science_photcal_apply_flat_corr
(cpl_frameset * frameset, const char * specphot_tag, 
 const char * master_specphot_tag, const char * flat_sed_tag,
 int* fluxcal, bool response_apply_flat_corr, int standard)
{
    
    //We are correcting a stdstar. Use what it was used to compute the response
    if(standard)
        return response_apply_flat_corr;
    else //We have a science
    {
        if(*fluxcal) //Logic only if flux calibration is going to be performed
        {
            bool have_flat_sed = false;
            if(cpl_frameset_count_tags(frameset, flat_sed_tag) > 0)
                have_flat_sed= true;
            
            //Check if input response has been flat-sed-corrected 
            bool resp_has_flat_corr = false;
            cpl_table * response_interp;
            if (cpl_frameset_count_tags(frameset, master_specphot_tag) != 0) {
                response_interp = dfs_load_table(frameset, 
                        master_specphot_tag, 2);
                if(cpl_table_has_column(response_interp, "RESPONSE_FFSED"))
                    resp_has_flat_corr = true;
            }
            else if (cpl_frameset_count_tags(frameset, specphot_tag) != 0){
                response_interp = dfs_load_table(frameset, specphot_tag, 2);
                if(cpl_table_has_column(response_interp, "RESPONSE_FFSED"))
                    resp_has_flat_corr = true;
            }

            if(resp_has_flat_corr && !have_flat_sed)
            {
                cpl_msg_warning(cpl_func, "The response is corrected with the"
                   " FLAT_SED but there is no FLAT_SED provided for the "
                   "science frame. Therefore the data are not flux-calibrated");
                *fluxcal = 0;
                return false;
            }
            
            if(!resp_has_flat_corr && have_flat_sed)
            {
                cpl_msg_warning(cpl_func, "FLAT_SED in science sof ignored as"
                        " no FLAT_SED correction was applied to the response");
                return false;
            }
            
            if(!resp_has_flat_corr && !have_flat_sed)
                return false;
            
            if(resp_has_flat_corr && have_flat_sed)
                return true;
        }
        else        //No flux calibration 
            return false;
    }
    return false;
}

static void vimos_shift_stdstar(cpl_image * stdstar_spectra, 
                                cpl_image * stdstar_spectra_error,
                                double dispersion, 
                                double resp_shift)
{
    cpl_size shift_pix = (cpl_size)std::floor(resp_shift / dispersion + 0.5);
    if(std::abs(shift_pix) >= cpl_image_get_size_x(stdstar_spectra))
        throw std::invalid_argument("stdstar shift is larger than spectrum size");
    cpl_image_shift(stdstar_spectra, shift_pix, 0);
    cpl_image_shift(stdstar_spectra_error, shift_pix, 0);
}
/**@}*/
