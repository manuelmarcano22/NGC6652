/* $Id: vimos_science.c,v 1.36 2013-08-23 10:17:39 cgarcia Exp $
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
 * $Date: 2013-08-23 10:17:39 $
 * $Revision: 1.36 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vimos_science_impl.h>
#include <vmutils.h>
#include <cpl.h>
#include <moses.h>
#include <fors_tools.h>
#include <vimos_dfs.h>
#include <vimos_response.h>
#include "wavelength_calibration.h"
#include "hdrl.h"

static int vimos_science_create(cpl_plugin *);
static int vimos_science_exec(cpl_plugin *);
static int vimos_science_destroy(cpl_plugin *);
static int vimos_science(cpl_parameterlist *, cpl_frameset *, int);
static int vimos_science_mult(cpl_parameterlist *, cpl_frameset *);

static bool vimos_science_photcal_apply_flat_corr_multiframe
(cpl_frameset * frameset, const char * specphot_tag, 
 const char * master_specphot_tag, const char * flat_sed_tag,
 int* fluxcal);

static char vimos_science_description[] =
"This recipe is used to reduce scientific spectra using the extraction\n"
"mask and the products created by the recipe vimos_calib. The spectra are\n"
"bias subtracted, flat fielded (if a normalised flat field is specified)\n"
"and remapped eliminating the optical distortions. The wavelength calibration\n"
"can be optionally upgraded using a number of sky lines: if no sky lines\n"
"catalog of wavelengths is specified, an internal one is used instead.\n"
"If the alignment to the sky lines is performed, the input dispersion\n"
"coefficients table is upgraded and saved to disk, and a new CCD wavelengths\n"
"map is created. A configuration table (typically depending on the grism in\n"
"use) may also be specified: this table contains a default recipe parameter\n" 
"setting to control the way spectra are extracted for a specific instrument\n"
"mode, as it is used for automatic run of the pipeline on Paranal and in\n" 
"Garching. If this table is specified, it will modify the default recipe\n" 
"parameter setting, with the exception of those parameters which have been\n" 
"explicitly modifyed on the command line. If a configuration table is not\n"
"specified, the input recipe parameters values will always be read from the\n"
"command line, or from an esorex configuration file if present, or from their\n"
"generic default values (which are rarely meaningful).\n" 
"Either a scientific or a standard star exposure can be specified in input.\n"
"Only in case of a standard star exposure input, the atmospheric extinction\n"
"table and a table with the physical fluxes of the observed standard star\n"
"must be specified in input, and a spectro-photometric table is created in\n"
"output. This table can then be input again to this recipe, always with an\n"
"atmospheric extinction table, and if a photometric calibration is requested\n"
"then flux calibrated spectra (in units of erg/cm/cm/s/Angstrom) are also\n"
"written in output.\n\n"
"Input files:\n\n"
"  DO category:                Type:       Explanation:         Required:\n"
"  MOS_SCIENCE                 Raw         Scientific exposure     Y\n"
"  or MOS_STANDARD             Raw         Standard star exposure  Y\n"
"\n"
"  MASTER_BIAS                 Calib       Master bias             Y\n"
"  SKY_LINE_CATALOG            Calib       Sky lines catalog       .\n"
"  MOS_MASTER_SCREEN_FLAT      Calib       Normalised flat field   .\n"
"  MOS_DISP_COEFF              Calib       Inverse dispersion      Y\n"
"  MOS_CURV_COEFF              Calib       Spectral curvature      Y\n"
"  MOS_SLIT_LOCATION           Calib       Slits positions table   Y\n"
"  CONFIG_TABLE                Calib       Configuration table     .\n"
"\n"
"  In case MOS_STANDARD is specified in input,\n"
"\n"
"  EXTINCT_TABLE               Calib       Atmospheric extinction  Y\n"
"  STD_FLUX_TABLE              Calib       Standard star flux      Y\n"
"  TELLURIC_CONTAMINATION      Calib       Telluric regions list    .\n"
"\n"
"\n"
"  In case a photometric calibration is requested for scientific\n"
"  data, the following inputs are mandatory:\n"
"\n"
"  EXTINCT_TABLE              Calib       Atmospheric extinction  Y\n"
"  MOS_SPECPHOT_TABLE         Calib       Response curves         Y\n"
"\n"
"  If requested for standard star data, the SPECPHOT_TABLE can be dropped:\n"
"  in this case the correction is applied using the SPECPHOT_TABLE produced\n"
"  in the same run.\n\n"
"Output files (acronym _STD_ replaces _SCI_ and _STANDARD_ replaces _SCIENCE_\n"
"in case input is MOS_STANDARD rather than MOS_SCIENCE):\n\n"
"  DO category:                Data type:  Explanation:\n"
"  MOS_SCIENCE_REDUCED           Image  Extracted scientific spectra\n"
"  MOS_SCI_SKY_REDUCED           Image  Extracted sky spectra\n"
"  MOS_SCI_ERROR_REDUCED         Image  Errors on extracted spectra\n"
"  MOS_UNMAPPED_SCIENCE          Image  Sky subtracted scientific spectra\n"
"  MOS_SCIENCE_EXTRACTED         Image  Rectified scientific spectra\n"
"  MOS_SCIENCE_SKY_EXTRACTED     Image  Rectified science spectra with sky\n"
"  MOS_SCIENCE_SKY               Image  Rectified sky spectra\n"
"  MOS_SCI_UNMAPPED_SKY          Image  Sky on CCD\n"
"  MOS_SCI_GLOBAL_SKY_SPECTRUM   Table  Global sky spectrum\n"
"  OBJECT_SCI_TABLE              Table  Positions of detected objects\n"
"\n"
"  Only if fringing correction is requested (dithered exposures):\n"
"  MOS_SCI_FRINGES               Image  Fringe map\n"
"\n"
"  Only if the sky-alignment of the wavelength solution is requested:\n"
"  MOS_SCI_SKYLINES_OFFSETS_SLIT Table  Sky lines offsets\n"
"  MOS_SCI_DISP_COEFF_SKY        Table  Upgraded dispersion coefficients\n"
"  MOS_SCI_WAVELENGTH_MAP_SKY    Image  Upgraded wavelength map\n"
"\n"
"  Only if a MOS_STANDARD is specified in input:\n"
"  MOS_SPECPHOT_TABLE            Table  Efficiency and response curves\n"
"\n"
"  Only if MOS_SPECPHOT_TABLE or MOS_MASTER_RESPONSE are specified in input:\n"
"  MOS_SCIENCE_FLUX_REDUCED      Image  Flux calibrated scientific spectra\n"
"  MOS_SCI_ERROR_FLUX_REDUCED    Image  Errors on flux calibrated spectra\n"
"  MOS_SCIENCE_FLUX_EXTRACTED    Image  Flux calibrated slit spectra\n\n";

/**
 * @brief    Build the list of available plugins, for this module. 
 *
 * @param    list    The plugin list
 *
 * @return   0 if everything is ok, -1 otherwise
 *
 * Create the recipe instance and make it available to the application 
 * using the interface. This function is exported.
 */

int cpl_plugin_get_info(cpl_pluginlist *list)
{
    cpl_recipe *recipe = (cpl_recipe*)cpl_calloc(1, sizeof *recipe);
    cpl_plugin *plugin = &recipe->interface;

    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmmosscience",
                    "Extraction of scientific spectra",
                    vimos_science_description,
                    "Carlo Izzo",
                    PACKAGE_BUGREPORT,
                    vimos_get_license(),
                    vimos_science_create,
                    vimos_science_exec,
                    vimos_science_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;
}


/**
 * @brief    Setup the recipe options    
 *
 * @param    plugin  The plugin
 *
 * @return   0 if everything is ok
 *
 * Defining the command-line/configuration parameters for the recipe.
 */

static int vimos_science_create(cpl_plugin *plugin)
{
    cpl_recipe    *recipe;
    cpl_parameter *p;


    /* 
     * Check that the plugin is part of a valid recipe 
     */

    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE) 
        recipe = (cpl_recipe *)plugin;
    else 
        return -1;

    /* 
     * Create the parameters list in the cpl_recipe object 
     */

    recipe->parameters = cpl_parameterlist_new(); 

    /*
     * Sky lines alignment
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.skyalign",
                                CPL_TYPE_INT,
                                "Polynomial order for sky lines alignment, "
                                "or -1 to avoid alignment",
                                "vimos.vmmosscience",
                                -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "skyalign");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Apply flat field
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.flatfield",
                                CPL_TYPE_BOOL,
                                "Apply flat field",
                                "vimos.vmmosscience",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "flatfield");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Global sky subtraction
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.skyglobal",
                                CPL_TYPE_BOOL,
                                "Subtract global sky spectrum from CCD",
                                "vimos.vmmosscience",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "skyglobal");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Local sky subtraction on extracted spectra
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.skymedian",
                                CPL_TYPE_BOOL,
                                "Sky subtraction from extracted slit spectra",
                                "vimos.vmmosscience",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "skymedian");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Local sky subtraction on CCD spectra
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.skylocal",
                                CPL_TYPE_BOOL,
                                "Sky subtraction from CCD slit spectra",
                                "vimos.vmmosscience",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "skylocal");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Cosmic rays removal
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.cosmics",
                                CPL_TYPE_BOOL,
                                "Eliminate cosmic rays hits (only if global "
                                "or local sky subtraction is also requested)",
                                "vimos.vmmosscience",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "cosmics");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Slit margin
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.slit_margin",
                                CPL_TYPE_INT,
                                "Number of pixels to exclude at each slit "
                                "in object detection and extraction",
                                "vimos.vmmosscience",
                                3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "slit_margin");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Extraction radius
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.ext_radius",
                                CPL_TYPE_INT,
                                "Maximum extraction radius for detected "
                                "objects (pixel)",
                                "vimos.vmmosscience",
                                6);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_radius");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Contamination radius
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.cont_radius",
                                CPL_TYPE_INT,
                                "Minimum distance at which two objects "
                                "of equal luminosity do not contaminate "
                                "each other (pixel)",
                                "vimos.vmmosscience",
                                0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "cont_radius");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Object extraction method
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.ext_mode",
                                CPL_TYPE_INT,
                                "Object extraction method: 0 = aperture, "
                                "1 = Horne optimal extraction",
                                "vimos.vmmosscience",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_mode");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Object detection threshold
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.detection",
                                CPL_TYPE_DOUBLE,
                                "Object detection threshold (ADU)",
                                "vimos.vmmosscience",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detection");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Number of nknots for the spline modeling the instrument response.
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.resp_fit_nknots",
                                CPL_TYPE_INT,
                                "Number of knots in spline fitting of the "
                                "instrument response. "
                                "(-1: No fitting. -2: Read from grism table)",
                                "vimos.vmmosscience",
                                -2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "resp_fit_nknots");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Order of polynomial modeling the instrument response.
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.resp_fit_degree",
                                CPL_TYPE_INT,
                                "Degree of polynomial in fitting of the "
                                "instrument response. "
                                "(-1: No fitting. -2: Read from grism table)",
                                "vimos.vmmosscience",
                                -2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "resp_fit_degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Mode for the wavelengths to ignore during response computation
     */
    
    p = cpl_parameter_new_value("vimos.vmmosscience.resp_ignore_mode",
                                CPL_TYPE_STRING,
                                "Types of lines/regions to ignore in response. "
                                "Valid ones are 'stellar_absorption', "
                                "'telluric' and 'command_line' (from parameter resp_ignore_lines)",
                                "vimos.vmmosscience",
                                "stellar_absorption,telluric,command_line");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "resp_ignore_mode");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * List of wavelengths (discrete and ranges) to ignore during response computation
     */
    
    p = cpl_parameter_new_value("vimos.vmmosscience.resp_ignore_points",
                                CPL_TYPE_STRING,
                                "Extra lines/regions to ignore in response. "
                                "Use a comma separated list of values. A range "
                                "can be specified like 4500.0-4600.0",
                                "vimos.vmmosscience",
                                "");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "resp_ignore_points");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Whether to use the flat sed to normalise spectra before applying
     * photometric calibration
     */
    
    p = cpl_parameter_new_value("vimos.vmmosscience.resp_use_flat_sed",
                                CPL_TYPE_STRING,
                                "Use the flat SED to normalise the observed spectra. "
                                "Value are true, false, grism_table.",
                                "vimos.vmmosscience",
                                "grism_table");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "resp_use_flat_sed");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Parameter to shift extracted standard star before computing the response 
     */
    
    p = cpl_parameter_new_value("vimos.vmmosscience.resp_shift",
                                CPL_TYPE_DOUBLE,
                                "The extracted standard star will be shifted "
                                "these many angstroms before using it to "
                                "compute the response. This is useful for "
                                "observed std stars not centered in the slits. "
                                "Positive values will shift the spectrum to "
                                "the red. Shift is given in Angstroms but no "
                                "fraction of pixels will be shifted.",
                                "vimos.vmmosscience",
                                0.);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "resp_shift");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Normalise output by exposure time
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.time_normalise",
                                CPL_TYPE_BOOL,
                                "Normalise output spectra by the exposure time",
                                "vimos.vmmosscience",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "time_normalise");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Look for a standard star in any frame classified as MOS_STANDARD
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.anyframe",
                                CPL_TYPE_BOOL,
                                "Look for a standard star in any frame "
                                "classified as MOS_STANDARD",
                                "vimos.vmmosscience",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "anyframe");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    /*
     * Type of alignment of dithered frames
     */

    p = cpl_parameter_new_enum("vimos.vmmosscience.alignment",
                               CPL_TYPE_STRING,
                               "Type of alignment of dithered frames, "
                               "either to the nearest neighbour pixel "
                               "or to fractions of pixel",
                               "vimos.vmmosscience",
                               "integer", 2,
                               "integer", "float");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "alignment");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    /*
     * Stacking method of dithered frames
     */

    p = cpl_parameter_new_enum("vimos.vmmosscience.stack_method",
                               CPL_TYPE_STRING,
                               "Frames combination method",
                               "vimos.vmmosscience",
                               "average", 4,
                               "average", "median", "minmax", "ksigma");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "stack_method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.vmmosscience.minrejection",
                                CPL_TYPE_INT,
                                "Number of lowest values to be rejected",
                                "vimos.vmmosscience",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "minrejection");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.vmmosscience.maxrejection",
                                CPL_TYPE_INT,
                                "Number of highest values to be rejected",
                                "vimos.vmmosscience",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "maxrejection");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.vmmosscience.klow",
                                CPL_TYPE_DOUBLE,
                                "Low threshold in ksigma method",
                                "vimos.vmmosscience",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "klow");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.vmmosscience.khigh",
                                CPL_TYPE_DOUBLE,
                                "High threshold in ksigma method",
                                "vimos.vmmosscience",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "khigh");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.vmmosscience.kiter",
                                CPL_TYPE_INT,
                                "Max number of iterations in ksigma method",
                                "vimos.vmmosscience",
                                999);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "kiter");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    /*
     * Dithering
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.dither",
                                CPL_TYPE_BOOL,
                                "Align dithered frames before stacking"
                                "(for multiple input frames)",
                                "vimos.vmmosscience",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "dither");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    /*
     * Compute offsets or read them from header?
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.dither.compute",
                                CPL_TYPE_BOOL,
                                "Compute offsets of dithered images from "
                                "detected objects (true), or read offsets "
                                "from header (false)",
                                "vimos.vmmosscience",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "compute");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    /*
     * Fringing correction
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.fringing",
                                CPL_TYPE_BOOL,
                                "Apply fringing correction "
                                "(only for dithered observations)",
                                "vimos.vmmosscience",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "fringing");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    /*
     * Min required offset for fringing correction
     */

    p = cpl_parameter_new_value("vimos.vmmosscience.fringing.offset",
                                CPL_TYPE_DOUBLE,
                                "Minimum required offset between exposures "
                                "for applying the sky fringing correction.",
                                "vimos.vmmosscience",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "offset");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    return 0;
}


/**
 * @brief    Execute the plugin instance given by the interface
 *
 * @param    plugin  the plugin
 *
 * @return   0 if everything is ok
 */

static int vimos_science_exec(cpl_plugin *plugin)
{
    cpl_recipe *recipe;
    int             status = 1;
   
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else 
        return -1;

    /* Issue a banner */
    vimos_print_banner();
    
    try
    {
        status = vimos_science_mult(recipe->parameters, recipe->frames); 
    }
    catch(std::exception& ex)
    {
        cpl_msg_error(cpl_func, "Recipe error: %s", ex.what());
    }
    catch(...)
    {
        cpl_msg_error(cpl_func, "An uncaught error during recipe execution");
    }

    return status;
}


/**
 * @brief    Destroy what has been created by the 'create' function
 *
 * @param    plugin  The plugin
 *
 * @return   0 if everything is ok
 */

static int vimos_science_destroy(cpl_plugin *plugin)
{
    cpl_recipe *recipe;
    
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE) 
        recipe = (cpl_recipe *)plugin;
    else 
        return -1;

    cpl_parameterlist_delete(recipe->parameters); 

    return 0;
}


/**
 * @brief    Interpret the command line options and execute the data processing
 *
 * @param    parlist     The parameters list
 * @param    frameset    The set-of-frames
 *
 * @return   0 if everything is ok
 */

static int vimos_science_mult(cpl_parameterlist *parlist, 
                              cpl_frameset *frameset)
{
    const char       *recipe            = "vmmosscience";
    char              version[80];
    const char       *slit_location_tag = "MOS_SLIT_LOCATION";
    const char       *curv_coeff_tag    = "MOS_CURV_COEFF";
    const char       *disp_coeff_tag    = "MOS_DISP_COEFF";
    const char       *global_dist_tag   = "GLOBAL_DISTORTION_TABLE";
    const char       *science_tag       = "MOS_SCIENCE";
    const char       *name              = NULL;
    cpl_frame        *frame             = NULL;
    cpl_table        *slits             = NULL;
    cpl_table        *subslits          = NULL;
    cpl_table        *maskslits         = NULL;
    cpl_table        *polytraces        = NULL;
    cpl_table        *idscoeff          = NULL;
    cpl_table        *global            = NULL;
    cpl_propertylist *header            = NULL;
    cpl_propertylist *wcal_header            = NULL;
    cpl_frame_type    type;
    cpl_frame_group   group;
    int               skyalign;
    int               i, multiplex, ngroups;
    int               status;
    int               mos;
    int               have_slit;
    int               have_curv;
    int               have_disp;
    int               error;
    const int         nx = 4096;        // Hard coded CCD size (rotated).
    const int         ny = 2048;


    snprintf(version, 80, "%s-%s", PACKAGE, PACKAGE_VERSION);

    have_slit = cpl_frameset_count_tags(frameset, slit_location_tag);
    have_curv = cpl_frameset_count_tags(frameset, curv_coeff_tag);
    have_disp = cpl_frameset_count_tags(frameset, disp_coeff_tag);

    if (have_slit == 0 && have_curv == 0 && have_disp == 0) {

        /*
         * Preliminary part, added to make use of a global distortion table
         * for creating the slit, the curv, and the disp tables, in case
         * they are (all) missing.
         */

        double      dispersion;
        double      startwavelength;
        double      endwavelength;
        double      refwave;

        wcal_header = dfs_load_header(frameset, disp_coeff_tag, 0);

        dispersion = cpl_propertylist_get_double(wcal_header, "ESO PRO WLEN INC"); 
                
        refwave = cpl_propertylist_get_double(wcal_header, "ESO PRO WLEN CEN");
        
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
        cpl_msg_info(cpl_func, " Using reference wavelength from calibrations: %f", refwave);
        cpl_msg_info(cpl_func, " Using dispersion from calibrations: %f", dispersion);
        
        if (cpl_frameset_count_tags(frameset, global_dist_tag) > 1) {
            cpl_msg_error(recipe, "Too many in input: %s", global_dist_tag);
            return -1;
        }

        if (cpl_frameset_count_tags(frameset, global_dist_tag) == 0) {
            cpl_msg_error(recipe, "Missing required input: %s, %s, and %s. "
                          "As a possible alternative, a %s may be "
                          "specified.", slit_location_tag, curv_coeff_tag,
                          disp_coeff_tag, global_dist_tag);
            return -1;
        }

        global = dfs_load_table(frameset, global_dist_tag, 1);

        if (global == NULL) {
            cpl_msg_error(recipe, "Cannot load global distortion table");
            return -1;
        }

        mos = cpl_frameset_count_tags(frameset, science_tag);

        if (mos == 0) {
            science_tag = "MOS_STANDARD";
            mos = cpl_frameset_count_tags(frameset, science_tag);
        }

        if (mos == 0) {
            cpl_msg_error(recipe, "Missing input scientific frame");
            return -1;
        }

        header = dfs_load_header(frameset, science_tag, 0);

        if (header == NULL) {
            cpl_msg_error(recipe, "Cannot load scientific frame header");
            return -1;
        }

        maskslits = mos_load_slits_vimos(header, 0);
        cpl_propertylist_delete(header); header = NULL;

        mos_assign_multiplex_group(maskslits);
        ngroups = 1 + cpl_table_get_column_max(maskslits, "group");

        mos_rotate_slits(maskslits, -1, 0, 0);

        for (i = 0; i < ngroups; i++) {
            cpl_image *spectra;
            cpl_image *dummy;

            cpl_table_select_all(maskslits);
            cpl_table_and_selected_int(maskslits, "group", CPL_EQUAL_TO, i);
            subslits = cpl_table_extract_selected(maskslits);
            slits = mos_build_slit_location(global, subslits, ny);
            polytraces = mos_build_curv_coeff(global, subslits, slits);

            /*
             * This is just to add "length" and "position" to the
             * slits table. A real overkill...
             */

            spectra = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
            dummy = mos_spatial_calibration(spectra, slits, polytraces, 
                                            refwave, startwavelength, 
                                            endwavelength, dispersion, 
                                            0, NULL);
            cpl_image_delete(dummy);
            cpl_image_delete(spectra);

            idscoeff = mos_build_disp_coeff(global, slits);
            cpl_table_delete(subslits);

            /*
             * Saving results, and mark them in the frameset as inputs,
             * not outputs
             */

            mos_rotate_slits(slits, 1, ny, nx);

            error = 0;

            if (i == 0) {
                error += dfs_save_image_null(frameset, NULL, parlist, 
                                             curv_coeff_tag,
                                             recipe, version);
                error += dfs_save_image_null(frameset, NULL, parlist, 
                                             disp_coeff_tag,
                                             recipe, version);
                error += dfs_save_image_null(frameset, NULL, parlist, 
                                             slit_location_tag,
                                             recipe, version);

                if (error) {
                    cpl_table_delete(maskslits);
                    cpl_table_delete(slits);
                    cpl_table_delete(idscoeff);
                    cpl_table_delete(polytraces);
                    return -1;
                }
            }

            error += dfs_save_table_ext(polytraces, curv_coeff_tag, NULL);
            error += dfs_save_table_ext(idscoeff, disp_coeff_tag, NULL);
            error += dfs_save_table_ext(slits, slit_location_tag, NULL);

            cpl_table_delete(slits);
            cpl_table_delete(idscoeff);
            cpl_table_delete(polytraces);

            if (error) {
                return -1;
            }
        }

        cpl_table_delete(global);

        frame = cpl_frameset_find(frameset, curv_coeff_tag);
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        frame = cpl_frameset_find(frameset, disp_coeff_tag);
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        frame = cpl_frameset_find(frameset, slit_location_tag);
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
    }

    if (cpl_frameset_count_tags(frameset, slit_location_tag) > 1) {
        cpl_msg_error(recipe, "Too many in input: %s", slit_location_tag);
        return -1;
    }

    frame = cpl_frameset_find(frameset, slit_location_tag);

    if (frame == NULL) {
        cpl_msg_error(recipe, "Missing required input: %s", slit_location_tag);
        return -1;
    }

    multiplex = cpl_fits_count_extensions(cpl_frame_get_filename(frame));

    if (multiplex == 1) {

        /*
         * No spectral multiplexing: classical handling
         */

        return vimos_science(parlist, frameset, 1);
    }

    if (0) {   // protect the rest
        cpl_msg_error(recipe, "The spectral multiplexing is %d, "
                      "I cannot handle such data yet!", multiplex);
        return -1;
    }


    /*
     * Parameters to be checked
     */

    skyalign = dfs_get_parameter_int(parlist,
                                     "vimos.vmmosscience.skyalign", NULL);

    if (skyalign > 2) {
        cpl_msg_error(recipe, "Max polynomial degree for sky alignment is 2");
        return -1;
    }


    /*
     * Check input
     */

    if (cpl_frameset_count_tags(frameset, curv_coeff_tag) > 1) {
        cpl_msg_error(recipe, "Too many in input: %s", curv_coeff_tag);
        return -1;
    }

    frame = cpl_frameset_find(frameset, curv_coeff_tag);

    if (frame == NULL) {
        cpl_msg_error(recipe, "Missing required input: %s", curv_coeff_tag);
        return -1;
    }

    if (multiplex != cpl_fits_count_extensions(cpl_frame_get_filename(frame))) {
        cpl_msg_error(recipe, "Expected %d data sections in input %s: " 
                      "%"CPL_SIZE_FORMAT" found.", multiplex, curv_coeff_tag,
                      cpl_fits_count_extensions(cpl_frame_get_filename(frame)));
        return -1;
    }

    if (cpl_frameset_count_tags(frameset, disp_coeff_tag) > 1) {
        cpl_msg_error(recipe, "Too many in input: %s", disp_coeff_tag);
        return -1;
    }

    frame = cpl_frameset_find(frameset, disp_coeff_tag);

    if (frame == NULL) {
        cpl_msg_error(recipe, "Missing required input: %s", disp_coeff_tag);
        return -1;
    }

    name = cpl_frame_get_filename(frame);

    if (multiplex != cpl_fits_count_extensions(name)) {
        cpl_msg_error(recipe, "Expected %d data sections in input %s: " 
                      "%"CPL_SIZE_FORMAT" found.", multiplex, disp_coeff_tag,
                      cpl_fits_count_extensions(name));
        return -1;
    }


    /*
     * Repeat for each extension
     */

    for (i = 0; i < multiplex; i++) {
        if (vimos_science(parlist, frameset, i + 1)) {
            cpl_msg_error(recipe, "Failure processing multiplex group %d", i);
            return -1;
        }

        /*
         * Search for all products
         */

        frame = cpl_frameset_get_first(frameset);
        while (frame) {
            char *tmpname;

            name  = cpl_frame_get_filename(frame);
            type  = cpl_frame_get_type(frame);
            group = cpl_frame_get_group(frame);

            if (group != CPL_FRAME_GROUP_PRODUCT) {
                frame = cpl_frameset_get_next(frameset);
                continue;
            }

            tmpname = cpl_sprintf("TMP_%s", name);

            if (i == 0) {

                if (type == CPL_FRAME_TYPE_IMAGE) {

                    /*
                     * If the product is an image, we need to
                     * recreate it with an empty primary array.
                     */

                    cpl_image *
                        image = cpl_image_load(name, CPL_TYPE_FLOAT, 0, 0);
                    cpl_propertylist *
                        this_header = cpl_propertylist_load(name, 0);
                    cpl_image * image_err = NULL;
                    cpl_propertylist * header_err = NULL;
                    //If it has error extension, read the extension
                    if (cpl_fits_count_extensions(name) > 0) {
                        image_err = cpl_image_load(name, CPL_TYPE_FLOAT, 0, 1);
                        header_err = cpl_propertylist_load(name, 1);
                    }
                    cpl_propertylist_save(this_header, name, CPL_IO_CREATE);
                    cpl_propertylist_erase_regexp(this_header,
                                       "^ESO DPR |^ARCFILE$|^ORIGFILE$", 0);
                    cpl_image_save(image, name, CPL_BPP_IEEE_FLOAT,
                                   this_header, CPL_IO_EXTEND);
                    cpl_image_delete(image);
                    cpl_propertylist_delete(this_header);
                    //If it has error extension, copy that one as well
                    if (image_err != NULL) {
                        cpl_image_save(image_err, name, CPL_BPP_IEEE_FLOAT,
		  	   header_err, CPL_IO_EXTEND);
                        cpl_image_delete(image_err);
                        cpl_propertylist_delete(header_err);
                    }
                }
                else if (type == CPL_FRAME_TYPE_TABLE) {

                    /*
                     * If the product is a table, there is
                     * nothing to be done: the table was created
                     * by the called recipe.
                     */

                     ;

                }
                else {

                    /*
                     * If the product is a header, it must be
                     * duplicated at the first extension.
                     */

                    cpl_propertylist *this_header = cpl_propertylist_load(name, 0);
                    cpl_propertylist_save(this_header, name, CPL_IO_EXTEND);
                    cpl_propertylist_delete(this_header);
                }


                /*
                 * Move products to product frameset.
                 */

                status = rename(name, tmpname);

                if (status) {
                    cpl_msg_error(recipe, "Cannot rename product (%s to %s).",
                                  name, tmpname);
                    return -1;
                }
            }
            else {
                if (type == CPL_FRAME_TYPE_IMAGE) {
                    cpl_image *
                    image = cpl_image_load(name, CPL_TYPE_FLOAT, 0, 0);
                    cpl_propertylist *
                    this_header = cpl_propertylist_load(name, 0);
                    cpl_propertylist_erase_regexp(this_header,
                                       "^ESO DPR |^ARCFILE$|^ORIGFILE$", 0);
                    cpl_image_save(image, tmpname, CPL_BPP_IEEE_FLOAT,
                                   this_header, CPL_IO_EXTEND);
                    if (cpl_error_get_code() != CPL_ERROR_NONE) {
                        cpl_error_reset();
                        cpl_msg_warning(recipe, 
                                    "Some fringes maps cannot be produced.");
                        cpl_propertylist_save(this_header, tmpname, CPL_IO_CREATE);
                        cpl_image_save(image, tmpname, CPL_BPP_IEEE_FLOAT,
                                       this_header, CPL_IO_EXTEND);
                    }
                    cpl_image_delete(image);
                    cpl_propertylist_delete(this_header);
                    //If it has error extension, copy that one as well
                    if (cpl_fits_count_extensions(name) > 0) {
                        cpl_image * image_err = 
                            cpl_image_load(name, CPL_TYPE_FLOAT, 0, 1);
                        cpl_propertylist * header_err =
                            cpl_propertylist_load(name, 1);
                        cpl_image_save(image_err, tmpname, CPL_BPP_IEEE_FLOAT,
		  	   header_err, CPL_IO_EXTEND);
                        cpl_image_delete(image_err);
                        cpl_propertylist_delete(header_err);
                    }
                }
                else if (type == CPL_FRAME_TYPE_TABLE) {
                    cpl_table *table = cpl_table_load(name, 1, 1);
                    cpl_propertylist *header = cpl_propertylist_load(name, 0);
                    cpl_table_save(table, NULL, header, tmpname, CPL_IO_EXTEND);
                    cpl_table_delete(table);
                    cpl_propertylist_delete(header);
                }
                else {

                    /*
                     * If the product is a header, it must be
                     * duplicated at the first extension.
                     */

                    cpl_propertylist *this_header = cpl_propertylist_load(name, 0);
                    cpl_propertylist_save(this_header, name, CPL_IO_EXTEND);
                    cpl_propertylist_delete(this_header);
                }
            }

            cpl_free(tmpname);
            if (i + 1 < multiplex) {

                /*
                 * Delete the products from the frameset, with
                 * the exception of the last round.
                 */

                cpl_frameset_erase_frame(frameset, frame);
                frame = cpl_frameset_get_first(frameset);
            }
            else {
                frame = cpl_frameset_get_next(frameset);
            }
        }
    }

    /*
     * When all is done, rename the products back to their proper name
     */

    for(cpl_size iframe = 0; iframe < cpl_frameset_get_size(frameset); ++iframe)
    {
        char *tmpname;
        frame = cpl_frameset_get_position(frameset, iframe);

        name  = cpl_frame_get_filename(frame); 
        group = cpl_frame_get_group(frame);

        if (group != CPL_FRAME_GROUP_PRODUCT) {
            continue;
        }

        tmpname = cpl_sprintf("TMP_%s", name); 
        status   = rename(tmpname, name);
        cpl_free(tmpname);

        if (status) {
            cpl_msg_error(recipe, "Cannot rename product (%s to %s).",
                          tmpname, name);
            return -1;
        }

    }

    return 0;
}

static int vimos_science(cpl_parameterlist *parlist, cpl_frameset *frameset,
                         int section)
{
    const char *recipe      = "vmmosscience";
    const char *science_tag = "MOS_SCIENCE";
    int         mos         = cpl_frameset_count_tags(frameset, science_tag);

    if (mos > 1) {
        char              version[80];
        const char       *curv_coeff_tag       = "MOS_CURV_COEFF";
        const char       *slit_location_tag    = "MOS_SLIT_LOCATION";
        const char       *disp_coeff_tag       = "MOS_DISP_COEFF";
        const char       *disp_coeff_sky_tag   = "MOS_SCI_DISP_COEFF_SKY";
        const char       *unmapped_science_tag = "MOS_UNMAPPED_SCIENCE";
        const char       *unmapped_science_err_tag = "MOS_UNMAPPED_SCIENCE_ERR";
        const char       *mapped_science_tag   = "MOS_SCIENCE_EXTRACTED";
        const char       *mapped_science_sky_tag = "MOS_SCIENCE_SKY_EXTRACTED";
        const char    *mapped_flux_science_tag = "MOS_SCIENCE_FLUX_EXTRACTED";
        const char       *object_table_tag     = "OBJECT_SCI_TABLE";
        const char       *unmapped_sky_tag     = "MOS_SCI_UNMAPPED_SKY";
        const char       *mapped_sky_tag       = "MOS_SCIENCE_SKY";
        const char       *reduced_science_tag  = "MOS_SCIENCE_REDUCED";
        const char   *reduced_flux_science_tag = "MOS_SCIENCE_FLUX_REDUCED";
        const char       *reduced_sky_tag      = "MOS_SCI_SKY_REDUCED";
        const char       *reduced_error_tag    = "MOS_SCI_ERROR_REDUCED";
        const char     *reduced_flux_error_tag = "MOS_SCI_ERROR_FLUX_REDUCED";
        const char       *fringes_tag          = "MOS_SCI_FRINGES";
        const char       *specphot_tag         = "MOS_SPECPHOT_TABLE";
        const char       *flat_sed_tag         = "MOS_FLAT_SED";
        const char  *skylines_offsets_tag = "MOS_SCI_SKYLINES_OFFSETS_SLIT";
        const char  *wavelength_map_sky_tag   = "MOS_SCI_WAVELENGTH_MAP_SKY";

        const char       *stack_method;
        const char       *alignment;
        double            refwave;
        double            startwavelength;
        double            endwavelength;
        double            dispersion;
        int               fringing;
        int               find_off;
        int               dither;
        int               cosmics;
        int               slit_margin;
        int               ext_radius;
        int               cont_radius;
        int               ext_mode;
        double            detection;
        int               photometry;
        int               skylocal, skyglobal;
        int               rotate      = 1;
        /* int               rotate_back = -1;*/
        cpl_frameset     *work;
        cpl_frame       **mos_science;
        cpl_frame        *frame;
        cpl_parameter    *param_time;
        cpl_table        *reference   = NULL;
        cpl_table        *objects     = NULL;
        cpl_image       **images      = NULL;
        cpl_image        *image       = NULL;
        cpl_image        *image_err   = NULL;
        cpl_image        *sky_image   = NULL;
        cpl_imagelist    *imagelist   = NULL;
        cpl_imagelist    *imagelist_err = NULL;
        cpl_image        *fringes     = NULL;
        cpl_image        *stacked     = NULL;
        cpl_image        *stacked_err = NULL;
        cpl_image        *mapped      = NULL;
        cpl_image        *mapped_var  = NULL;
        cpl_image        *smapped     = NULL;
        cpl_image        *smapped_nosky = NULL;
        cpl_image        *smapped_nosky_var = NULL;
        cpl_image        *sky_stacked = NULL;
        cpl_image        *sci_sky_mapped  = NULL;
        cpl_image        *sky_mapped  = NULL;
        cpl_image        *sky_smapped = NULL;
        cpl_table        *polytraces  = NULL;
        cpl_table        *idscoeff    = NULL;
        cpl_table        *slits       = NULL;
        cpl_table        *grism_table = NULL;
        cpl_table        *offsets     = NULL;
        cpl_propertylist *wcal_header      = NULL;
        cpl_propertylist *header      = NULL;
        cpl_propertylist *sort_col    = NULL;
        char             *name;
        double            offset, min_fring_offset, min_offset;
        double            a_offset, d_offset, ref_alpha, ref_delta, cos_delta;
        double            pix_scale, value;
        double            gain;
        double            ron;
        double            alltime;
        double            airmass;
        int               time_normal;
        int               min_reject;
        int               max_reject;
        double            klow;
        double            khigh;
        int               kiter;
        int               int_alignment = 0;
        int               status;
        int               nx, ny;
        int               i;
        int               done = 0;
        int               skyalign;


        snprintf(version, 80, "%s-%s", PACKAGE, PACKAGE_VERSION);

        grism_table = dfs_load_table(frameset, "CONFIG_TABLE", 1);

        cosmics = dfs_get_parameter_bool(parlist,
                                         "vimos.vmmosscience.cosmics", NULL);

        alignment = dfs_get_parameter_string(parlist,
                    "vimos.vmmosscience.alignment", NULL);

        if (strcmp(alignment, "integer") == 0) {
            int_alignment = 1;
        }

        stack_method = dfs_get_parameter_string(parlist,
                       "vimos.vmmosscience.stack_method", NULL);

        if (strcmp(stack_method, "minmax") == 0) {
            min_reject = dfs_get_parameter_int(parlist,
                         "vimos.vmmosscience.minrejection", NULL);
            if (min_reject < 0) {
                cpl_msg_error(recipe, "Invalid number of lower rejections");
                return -1;
            }
    
            max_reject = dfs_get_parameter_int(parlist,
                         "vimos.vmmosscience.maxrejection", NULL);

            if (max_reject < 0) {
                cpl_msg_error(recipe, "Invalid number of upper rejections");
                return -1;
            }
        }
    
        if (strcmp(stack_method, "ksigma") == 0) {
            klow  = dfs_get_parameter_double(parlist,
                                             "vimos.vmmosscience.klow", NULL);
            if (klow < 0.1) {
                cpl_msg_error(recipe, "Invalid lower K-sigma");
                return -1;
            }
    
            khigh = dfs_get_parameter_double(parlist,
                                             "vimos.vmmosscience.khigh", NULL);
            if (khigh < 0.1) {
                cpl_msg_error(recipe, "Invalid lower K-sigma");
                return -1;
            }
    
            kiter = dfs_get_parameter_int(parlist,
                                          "vimos.vmmosscience.kiter", NULL);
            if (kiter < 1) {
                cpl_msg_error(recipe, "Invalid number of iterations");
                return -1;
            }
        }

        slit_margin = dfs_get_parameter_int(parlist,
                                            "vimos.vmmosscience.slit_margin",
                                            NULL);
        ext_radius = dfs_get_parameter_int(parlist,
                                           "vimos.vmmosscience.ext_radius",
                                           NULL);
        cont_radius = dfs_get_parameter_int(parlist,
                                            "vimos.vmmosscience.cont_radius",
                                            NULL);
        ext_mode = dfs_get_parameter_int(parlist,
                                         "vimos.vmmosscience.ext_mode",
                                         NULL);

        detection = dfs_get_parameter_double(parlist,
                                             "vimos.vmmosscience.detection", 
                                             NULL);

        dither = dfs_get_parameter_bool(parlist, 
                                          "vimos.vmmosscience.dither", NULL);

        fringing = dfs_get_parameter_bool(parlist, 
                                          "vimos.vmmosscience.fringing", NULL);

        if (fringing && !dither) {
            cpl_msg_warning(recipe, "Fringing correction cannot be "
                            "applied when --dither=false.");
            fringing = 0;
        }

        min_fring_offset = dfs_get_parameter_double(parlist, 
                                 "vimos.vmmosscience.fringing.offset", NULL);

        find_off = dfs_get_parameter_bool(parlist, 
                                 "vimos.vmmosscience.dither.compute", NULL);

        skyglobal = dfs_get_parameter_bool(parlist, 
                                 "vimos.vmmosscience.skyglobal", NULL);
        skylocal  = dfs_get_parameter_bool(parlist, 
                                 "vimos.vmmosscience.skylocal", NULL);
        //skymedian = dfs_get_parameter_bool(parlist, 
        //                          "vimos.vmmosscience.skymedian", NULL);

        if (!skyglobal && !skylocal) {
            cpl_msg_error(recipe, "In the case of dithered observations "
                          "the sky subtraction must be performed before "
                          "the alignment of the CCD frames: either -skylocal "
                          "or -skyglobal must be set.");
            return -1;
        }

        //Getting wavelengths parameters form the headers
        wcal_header = dfs_load_header(frameset, disp_coeff_tag, 0);

        dispersion = cpl_propertylist_get_double(wcal_header, "ESO PRO WLEN INC"); 
                
        refwave = cpl_propertylist_get_double(wcal_header, "ESO PRO WLEN CEN");
        
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
        cpl_msg_info(cpl_func, " Using reference wavelength from calibrations: %f", refwave);
        cpl_msg_info(cpl_func, " Using dispersion from calibrations: %f", dispersion);

        cpl_table_delete(grism_table); grism_table = NULL;

        frame = cpl_frameset_find(frameset, science_tag);

        while (frame) {
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_RAW);
            frame = cpl_frameset_find(frameset, NULL);
        }

        frame = cpl_frameset_find(frameset, "MASTER_BIAS");
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        frame = cpl_frameset_find(frameset, "SKY_LINE_CATALOG");
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        frame = cpl_frameset_find(frameset, "MOS_MASTER_SCREEN_FLAT");
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        frame = cpl_frameset_find(frameset, "MOS_DISP_COEFF");
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        frame = cpl_frameset_find(frameset, "MOS_CURV_COEFF");
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        frame = cpl_frameset_find(frameset, "MOS_SLIT_LOCATION");
        if (frame)
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);


        /*
         * Check if a photometric table exists and if so activate the
         * photometric correction.
         */

        int have_phot;

        have_phot = cpl_frameset_count_tags(frameset, specphot_tag);
        if (have_phot == 0) {
            specphot_tag = "MOS_MASTER_RESPONSE";
        }
        have_phot = cpl_frameset_count_tags(frameset, specphot_tag);

        if (have_phot == 0) 
            photometry = 0;
        else if (have_phot > 1) 
        {
            cpl_msg_error(recipe, "Too many in input: %s", specphot_tag);
            return -1;
        }
        else 
        {
            cpl_msg_info(recipe, "Doing spectrophotometry using file %s", 
                         specphot_tag);
            photometry = 1;
        }

        if (photometry) {
            if (cpl_frameset_count_tags(frameset, "EXTINCT_TABLE") == 0) {
                cpl_msg_error(recipe,
                        "An EXTINCT_TABLE was not found in input: "
                        "the requested photometric calibrated "
                        "spectra cannot be produced.");
                return -1;
            }

            if (cpl_frameset_count_tags(frameset, "EXTINCT_TABLE") > 1) {
                cpl_msg_error(recipe, "Too many in input: EXTINCT_TABLE");
                return -1;
            }


        }


        /*
         * Disable time normalisation for single runs
         */

        param_time = cpl_parameterlist_find(parlist, 
                                       "vimos.vmmosscience.time_normalise");
        time_normal = cpl_parameter_get_bool(param_time);
        cpl_parameter_set_bool(param_time, 0);


        /*
         * Collect all scientific frames from input SOF
         */

        mos_science = (cpl_frame **)cpl_calloc(mos, sizeof(cpl_frame*));

        mos_science[0] = cpl_frameset_find(frameset, science_tag);
        for (i = 1; i < mos; i++) {
            mos_science[i] = cpl_frameset_find(frameset, NULL);
        }


        /*
         * Get total exposure time
         */

        alltime = 0.0;

        for (i = 0; i < mos; i++) {
            double time;

            if (i)
                header = dfs_load_header(frameset, NULL, 0);
            else
                header = dfs_load_header(frameset, science_tag, 0);

            if (header == NULL) {
                cpl_msg_error(recipe, "Cannot load scientific frame header");
                return -1;
            }

            alltime += time = cpl_propertylist_get_double(header, "EXPTIME");

            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_msg_error(recipe, "Missing keyword EXPTIME in scientific "
                              "frame header");
                cpl_propertylist_delete(header);
                return -1;
            }

            cpl_propertylist_delete(header);

            cpl_msg_info(recipe, "Scientific frame %d exposure time: %.2f s",
                         i+1, time);
        }


        /*
         * Process scientific frames one at a time
         */

        work = cpl_frameset_duplicate(frameset);

        for (i = 0; i < mos; i++) {
            cpl_frameset_erase(work, science_tag);
            /* Spectrophotometric disable on individual frames */
            cpl_frameset_erase(work, specphot_tag); 
            cpl_frameset_insert(work, cpl_frame_duplicate(mos_science[i]));
            if (vimos_science_impl(work, parlist, section)) {
                cpl_frameset_delete(work);
                cpl_msg_error(recipe, 
                              "Failure processing science exposure %d (%s)", i, 
                              cpl_frame_get_filename(mos_science[i]));
                return -1;
            }


            /*
             * Cleanup work SOF from products 
             * (not really necessary, but clean)
             */

            cpl_frameset_erase(work, skylines_offsets_tag);
            cpl_frameset_erase(work, unmapped_sky_tag);
            cpl_frameset_erase(work, unmapped_science_tag);
            cpl_frameset_erase(work, unmapped_science_err_tag);
            cpl_frameset_erase(work, wavelength_map_sky_tag);
            cpl_frameset_erase(work, mapped_sky_tag);
            cpl_frameset_erase(work, object_table_tag);
            cpl_frameset_erase(work, reduced_science_tag);
            cpl_frameset_erase(work, reduced_sky_tag);
            cpl_frameset_erase(work, reduced_error_tag);
            cpl_frameset_erase(work, disp_coeff_sky_tag);
            cpl_frameset_erase(work, mapped_science_sky_tag);
            cpl_frameset_erase(work, mapped_science_tag);

            if (photometry) {
                cpl_frameset_erase(work, reduced_flux_error_tag);
                cpl_frameset_erase(work, reduced_flux_science_tag);
                cpl_frameset_erase(work, mapped_flux_science_tag);
            }


            /*
             * Of all products, only keep the object tables, the
             * unmapped sky, and the unmapped science, unmapped science error.
             */

            done = 0;
            name = cpl_sprintf("object_table_%d.fits", i);
            status = rename("object_sci_table.fits", name);
            if (status) {
                done = 1;
                cpl_msg_error(recipe, "Cannot rename product "
                              "(object_sci_table.fits to %s).", name);
            }
            cpl_free(name);
            name = cpl_sprintf("mos_unmapped_sky_%d.fits", i);
            status |= rename("mos_sci_unmapped_sky.fits", name);
            if (status && !done) {
                done = 1;
                cpl_msg_error(recipe, "Cannot rename product "
                              "(mos_sci_unmapped_sky.fits to %s).", name);
            }
            cpl_free(name);
            name = cpl_sprintf("mos_unmapped_science_%d.fits", i);
            status |= rename("mos_unmapped_science.fits", name);
            if (status && !done) {
                done = 1;
                cpl_msg_error(recipe, "Cannot rename product "
                              "(mos_unmapped_science.fits to %s).", name);
            }
            cpl_free(name);
            name = cpl_sprintf("mos_unmapped_science_err_%d.fits", i);
            status |= rename("mos_unmapped_science_err.fits", name);
            if (status && !done) {
                done = 1;
                cpl_msg_error(recipe, "Cannot rename product "
                              "(mos_unmapped_science_err.fits to %s).", name);
            }

            if (status) {
                cpl_free(mos_science);
                cpl_frameset_delete(work);
                cpl_free(name);
                return -1;
            }

            cpl_free(name);
        }

        cpl_free(mos_science);
        cpl_frameset_delete(work);

        if (fringing) {

            /*
             * Produce fringing map by median stacking
             */

            imagelist = cpl_imagelist_new();

            for (i = 0; i < mos; i++) {
                name = cpl_sprintf("mos_unmapped_science_%d.fits", i);
                image = cpl_image_load(name, CPL_TYPE_FLOAT, 0, 0);
                cpl_imagelist_set(imagelist, image,
                                  cpl_imagelist_get_size(imagelist));
            }

            fringes = cpl_imagelist_collapse_median_create(imagelist);
            cpl_imagelist_delete(imagelist);
        }


        /*
         * Now stack all frames and extract final objects
         */

        reference = cpl_table_load("object_table_0.fits", 1, 1);
        stacked = cpl_image_load("mos_unmapped_science_0.fits", 
                                 CPL_TYPE_FLOAT, 0, 0);
        stacked_err = cpl_image_load("mos_unmapped_science_err_0.fits", 
                                 CPL_TYPE_FLOAT, 0, 0);
        sky_stacked = cpl_image_load("mos_unmapped_sky_0.fits", 
                                     CPL_TYPE_FLOAT, 0, 0);

        min_offset = 0.0;

        if (dither) {

            /*
             * Compute minimum offset
             */

            offsets = cpl_table_new(mos);
            cpl_table_new_column(offsets, "offset", CPL_TYPE_DOUBLE);
            cpl_table_set_double(offsets, "offset", 0, 0.0);

            if (!find_off) {
                header = cpl_propertylist_load("object_table_0.fits", 0);
                pix_scale = cpl_propertylist_get_double(header, 
                                                        "ESO INS PIXSCALE");
                ref_alpha = cpl_propertylist_get_double(header, "RA");
                ref_delta = cpl_propertylist_get_double(header, "DEC");
                cos_delta = fabs(cos(CPL_MATH_PI * ref_delta / 180.0));
                cpl_propertylist_delete(header);
            }

            for (i = 1; i < mos; i++) {
                name = cpl_sprintf("object_table_%d.fits", i);
                if (find_off) {
                    objects = cpl_table_load(name, 1, 1);
                    if (mos_compute_offset(reference, objects, &offset)) {
                        cpl_msg_error(recipe, "Cannot compute offset.");
                        cpl_table_delete(objects);
                        cpl_table_delete(offsets);
                        cpl_image_delete(stacked);
                        cpl_image_delete(stacked_err);
                        cpl_image_delete(sky_stacked);
                        return -1;
                    }
                    cpl_table_delete(objects);
                }
                else {
                    header = cpl_propertylist_load(name, 0);
                    if (cos_delta > 0.0001) {
                        value = cpl_propertylist_get_double(header, "RA");
                        a_offset = 3600.0 
                               * ((ref_alpha - value) / cos_delta) / pix_scale;
                    }
                    else {
                        a_offset = 0.0;
                    }
                    value = cpl_propertylist_get_double(header, "DEC");
                    d_offset = 3600.0 * (ref_delta - value) / pix_scale;
                    if (fabs(a_offset) > fabs(d_offset)) {
                        offset = -a_offset;
                    }
                    else {
                        offset = d_offset;
                    }
                    cpl_propertylist_delete(header);
                }
                cpl_free(name);

                cpl_table_set_double(offsets, "offset", i, offset);
    
                cpl_msg_info(recipe, "Frame %d offset relative "
                             "to frame 1: %.2f pix\n", i + 1, offset); 
    
                if (int_alignment) {
                    offset = floor(offset + 0.5);
                    cpl_msg_info(recipe, "Nearest neighbour offset applied is "
                                 "%.0f pix\n", offset); 
                }
            }

            sort_col = cpl_propertylist_new();
            cpl_propertylist_append_bool(sort_col, "offset", 0);
            cpl_table_sort(offsets, sort_col);
            cpl_propertylist_delete(sort_col);

            min_offset = cpl_table_get_double(offsets, "offset", 1, NULL)
                       - cpl_table_get_double(offsets, "offset", 0, NULL);

            for (i = 2; i < mos; i++) {
                double 
                delta = cpl_table_get_double(offsets, "offset", i    , NULL)
                      - cpl_table_get_double(offsets, "offset", i - 1, NULL);
                if (min_offset > delta)
                    min_offset = delta;
            }

            cpl_table_delete(offsets);

            cpl_msg_info(recipe, "Minimum observed offset between frames is "
                         "%.2f pixels.", min_offset);
        }

        if (fringing) {
            if (min_offset < min_fring_offset) {
                cpl_msg_warning(recipe, "At least a %.2f pixel offset is "
                                "required for fringing correction. The "
                                "sky fringing correction is now DISABLED: "
                                "run the recipe again with a lower threshold "
                                "if you want to try the fringing correction "
                                "anyway.", min_fring_offset);
                fringing = 0;
            }
        }

        if (fringing) {
            cpl_image_subtract(stacked, fringes);
            cpl_image_add(sky_stacked, fringes);
        }

// Added just for eliminating wcs like in other images
        cpl_image_save(stacked, "mos_unmapped_science_0.fits", 
                       CPL_BPP_IEEE_FLOAT, NULL, CPL_IO_DEFAULT);
// end of debug line

        for (i = 1; i < mos; i++) {
            if (dither) {
                name = cpl_sprintf("object_table_%d.fits", i);
                if (find_off) {
                    objects = cpl_table_load(name, 1, 1);
                    mos_compute_offset(reference, objects, &offset);
                    cpl_table_delete(objects);
                }
                else {
                    header = cpl_propertylist_load(name, 0);
                    if (cos_delta > 0.0001) {
                        value = cpl_propertylist_get_double(header, "RA");
                        a_offset = 3600.0
                               * ((ref_alpha - value) / cos_delta) / pix_scale;
                    }
                    else {
                        a_offset = 0.0;
                    }
                    value = cpl_propertylist_get_double(header, "DEC");
                    d_offset = 3600.0 * (ref_delta - value) / pix_scale;
                    if (fabs(a_offset) > fabs(d_offset)) {
                        offset = -a_offset;
                    }
                    else {
                        offset = d_offset;
                    }
                    cpl_propertylist_delete(header);
                }
                cpl_free(name);
            }
            else {
                offset = 0.0;
            }

            name = cpl_sprintf("mos_unmapped_science_%d.fits", i);
            image = cpl_image_load(name, CPL_TYPE_FLOAT, 0, 0);

            name = cpl_sprintf("mos_unmapped_science_err_%d.fits", i);
            image_err = cpl_image_load(name, CPL_TYPE_FLOAT, 0, 0);
            
            if (fringing) {
                cpl_image_subtract(image, fringes);
            }

            mos_image_shift(image, offset, 0.0);
            cpl_image_save(image, name, CPL_BPP_IEEE_FLOAT, 
                           NULL, CPL_IO_DEFAULT);
            cpl_free(name);

            /*
             * Here stacking sky frames, they are always averaged
             */

            name = cpl_sprintf("mos_unmapped_sky_%d.fits", i);
            sky_image = cpl_image_load(name, CPL_TYPE_FLOAT, 0, 0);

            if (fringing) {
                cpl_image_add(sky_image, fringes);
            }

            mos_image_shift(sky_image, offset, 0.0);
            cpl_free(name);
            cpl_image_add(sky_stacked, sky_image);
            cpl_image_delete(sky_image);

            if (i == 1) {
                imagelist = cpl_imagelist_new();
                cpl_imagelist_set(imagelist, stacked,
                                      cpl_imagelist_get_size(imagelist));
                imagelist_err = cpl_imagelist_new();
                cpl_imagelist_set(imagelist_err, stacked_err,
                                      cpl_imagelist_get_size(imagelist_err));
            }
            cpl_imagelist_set(imagelist, image,
                             cpl_imagelist_get_size(imagelist));
            cpl_imagelist_set(imagelist_err, image_err,
                                  cpl_imagelist_get_size(imagelist_err));
        }
        cpl_table_delete(reference);

        cpl_image_divide_scalar(sky_stacked, mos);


        /*
         * Do the image stacking
         */
        hdrl_parameter * stackmethod_par = NULL;
        
        if (strcmp(stack_method, "average") == 0) {
            stackmethod_par = hdrl_collapse_mean_parameter_create();
        }

        if (strcmp(stack_method, "median") == 0) {
            stackmethod_par = hdrl_collapse_median_parameter_create();
        }

        if (strcmp(stack_method, "minmax") == 0) {
            stackmethod_par = hdrl_collapse_minmax_parameter_create(
                    min_reject, max_reject);
        }

        if (strcmp(stack_method, "ksigma") == 0) {
            stackmethod_par = hdrl_collapse_sigclip_parameter_create
              (klow, khigh, kiter);
        }

        
        //Transform to HDRL
        hdrl_imagelist * im_list = hdrl_imagelist_new();
        for(size_t idx = 0; idx < cpl_imagelist_get_size(imagelist); idx++)
        {
            cpl_image * im = cpl_imagelist_get(imagelist, idx);
            cpl_image * im_err = cpl_imagelist_get(imagelist_err, idx);
            hdrl_image * tmp = hdrl_image_create(im, im_err);
            
            hdrl_imagelist_set(im_list, tmp, idx);
        }
        
        //Do the stacking
        cpl_image * contrib;
        hdrl_image * stacked_hdrl = NULL;
        if(hdrl_imagelist_collapse(im_list, stackmethod_par ,&stacked_hdrl, 
                                &contrib) != CPL_ERROR_NONE)
        {   
            cpl_msg_error(recipe, "Cannot stack science images");
            cpl_image_delete(sky_stacked);
            return -1;
        }
        cpl_image_delete(contrib);
        hdrl_imagelist_delete(im_list);
        cpl_imagelist_delete(imagelist);
        cpl_imagelist_delete(imagelist_err);
        
        stacked = cpl_image_cast(hdrl_image_get_image(stacked_hdrl), CPL_TYPE_FLOAT);
        stacked_err = cpl_image_cast(hdrl_image_get_error(stacked_hdrl), CPL_TYPE_FLOAT);
        hdrl_image_delete(stacked_hdrl);
        
        
        /*
         * Restore main configuration
         */

        cpl_parameter_set_bool(param_time, time_normal);

        header = dfs_load_header(frameset, science_tag, 0);

        if (photometry) {
            airmass = fors_get_airmass(header);
            if (airmass < 0.0) {
                cpl_msg_error(recipe, "Missing airmass information in "
                              "scientific frame header");
                cpl_image_delete(fringes);
                cpl_image_delete(stacked);
                cpl_image_delete(stacked_err);
                cpl_image_delete(sky_stacked);
                cpl_propertylist_delete(header);
                return -1;
            }
        }

        gain = cpl_propertylist_get_double(header, "ESO DET OUT1 CONAD");
        ron = cpl_propertylist_get_double(header, "ESO DET OUT1 RON");
        ron /= gain;     /* Convert from electrons to ADU */

        cpl_propertylist_update_double(header, "ESO PRO EXPTTOT", alltime);

        if (fringing) {
            if (dfs_save_image(frameset, fringes, fringes_tag,
                               header, parlist, recipe, version)) {
                cpl_image_delete(fringes);
                cpl_image_delete(stacked);
                cpl_image_delete(stacked_err);
                cpl_image_delete(sky_stacked);
                cpl_propertylist_delete(header);
                return -1;
            }
            cpl_image_delete(fringes);
        }

        if (dfs_save_image(frameset, stacked, unmapped_science_tag,
                           header, parlist, recipe, version)) {
            cpl_image_delete(stacked);
            cpl_image_delete(stacked_err);
            cpl_image_delete(sky_stacked);
            cpl_propertylist_delete(header);
            return -1;
        }

        if (dfs_save_image(frameset, sky_stacked, unmapped_sky_tag,
                           header, parlist, recipe, version)) {
            cpl_image_delete(stacked);
            cpl_image_delete(stacked_err);
            cpl_image_delete(sky_stacked);
            cpl_propertylist_delete(header);
            return -1;
        }

        cpl_image_turn(stacked, rotate);
        cpl_image_turn(stacked_err, rotate);
        cpl_image_turn(sky_stacked, rotate);

        nx = cpl_image_get_size_x(stacked);
        ny = cpl_image_get_size_y(stacked);

        polytraces = dfs_load_table(frameset, curv_coeff_tag, section);

        slits = dfs_load_table(frameset, slit_location_tag, section);
        mos_rotate_slits(slits, -rotate, nx, ny);

// FIXME: Qui potremmo tirare su una tabella allineata al cielo, se esiste.
        idscoeff = dfs_load_table(frameset, disp_coeff_tag, section);

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Processing stacked scientific spectra...");
        cpl_msg_indent_more();

        smapped_nosky = mos_spatial_calibration(stacked, slits, polytraces, 
                                          refwave, startwavelength, 
                                          endwavelength, dispersion, 
                                          1, NULL);
        cpl_image * stacked_var = cpl_image_power_create(stacked_err, 2);
        smapped_nosky_var = mos_spatial_calibration(stacked_var, slits, polytraces, 
                                          refwave, startwavelength, 
                                          endwavelength, dispersion, 1, NULL);

        cpl_image_delete(stacked);
        cpl_image_delete(stacked_err);
        cpl_image_delete(stacked_var);

        sky_smapped = mos_spatial_calibration(sky_stacked, slits, polytraces, 
                                              refwave, startwavelength, 
                                              endwavelength, dispersion, 
                                              1, NULL);

        cpl_image_delete(sky_stacked);

        mapped = mos_wavelength_calibration(smapped_nosky, refwave,
                                            startwavelength, endwavelength,
                                            dispersion, idscoeff, 1);

        mapped_var = mos_wavelength_calibration(smapped_nosky_var, refwave,
                                            startwavelength, endwavelength,
                                            dispersion, idscoeff, 1);

        sky_mapped = mos_wavelength_calibration(sky_smapped, refwave,
                                                startwavelength, endwavelength,
                                                dispersion, idscoeff, 1);
        cpl_image_delete(smapped_nosky);
        cpl_image_delete(smapped_nosky_var);

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

        if (time_normal) {
            cpl_image *dummy = cpl_image_divide_scalar_create(mapped, 
                                                              alltime / mos);
            if (dfs_save_image(frameset, dummy, mapped_science_tag, 
                               header, parlist, recipe, version)) {
                cpl_image_delete(dummy);
                cpl_propertylist_delete(header);
                cpl_image_delete(mapped);
                cpl_image_delete(sky_mapped);
                cpl_image_delete(sci_sky_mapped);
                cpl_table_delete(slits);
                return -1;
            }
            cpl_image * mapped_err = cpl_image_power_create(mapped_var, 0.5);
            cpl_image_divide_scalar(mapped_err, alltime / mos);
            if (dfs_save_image_ext(mapped_err, mapped_science_tag, header))
            {
                cpl_image_delete(dummy);
                cpl_propertylist_delete(header);
                cpl_image_delete(mapped);
                cpl_image_delete(sky_mapped);
                cpl_image_delete(sci_sky_mapped);
                cpl_table_delete(slits);
                cpl_image_delete(mapped_err);
                return -1;
            }

            cpl_image_delete(mapped_err);
            cpl_image_delete(dummy);
        }
        else {
            if (dfs_save_image(frameset, mapped, mapped_science_tag,
                               header, parlist, recipe, version)) {
                cpl_propertylist_delete(header);
                cpl_image_delete(mapped);
                cpl_image_delete(sky_mapped);
                cpl_image_delete(sci_sky_mapped);
                cpl_table_delete(slits);
                return -1;
            }
            cpl_image * mapped_err = cpl_image_power_create(mapped_var, 0.5);
            if (dfs_save_image_ext(mapped_err, mapped_science_tag, header))
            {
                cpl_propertylist_delete(header);
                cpl_image_delete(mapped);
                cpl_image_delete(sky_mapped);
                cpl_image_delete(sci_sky_mapped);
                cpl_table_delete(slits);
                cpl_image_delete(mapped_err);
                return -1;
            }
            cpl_image_delete(mapped_err);
        }

        sci_sky_mapped = cpl_image_add_create(mapped, sky_mapped);

        if (time_normal) {
            cpl_image *dummy = cpl_image_divide_scalar_create(sci_sky_mapped,
                                                              alltime / mos);
            if (dfs_save_image(frameset, dummy, mapped_science_sky_tag,
                               header, parlist, recipe, version)) {
                cpl_image_delete(dummy);
                cpl_propertylist_delete(header);
                cpl_image_delete(mapped);
                cpl_image_delete(sky_mapped);
                cpl_image_delete(sci_sky_mapped);
                cpl_table_delete(slits);
                return -1;
            }
            //It is regarded as if the sky subtraction doesn't have any error
            cpl_image * mapped_sky_err = cpl_image_power_create(mapped_var, 0.5);
            cpl_image_divide_scalar(mapped_sky_err, alltime / mos);
            if (dfs_save_image_ext(mapped_sky_err, mapped_science_sky_tag, header))
            {
                cpl_image_delete(dummy);
                cpl_propertylist_delete(header);
                cpl_image_delete(mapped);
                cpl_image_delete(sky_mapped);
                cpl_image_delete(sci_sky_mapped);
                cpl_table_delete(slits);
                cpl_image_delete(mapped_sky_err);
                return -1;
            }

            cpl_image_delete(mapped_sky_err);
            cpl_image_delete(dummy);
        }
        else {
            if (dfs_save_image(frameset, sci_sky_mapped, mapped_science_sky_tag,
                               header, parlist, recipe, version)) {
                cpl_propertylist_delete(header);
                cpl_image_delete(mapped);
                cpl_image_delete(sky_mapped);
                cpl_image_delete(sci_sky_mapped);
                cpl_table_delete(slits);
                return -1;
            }
            //It is regarded as if the sky subtraction doesn't have any error
            cpl_image * mapped_sky_err = cpl_image_power_create(mapped_var, 0.5);
            if (dfs_save_image_ext(mapped_sky_err, mapped_science_sky_tag, header))
            {
                cpl_propertylist_delete(header);
                cpl_image_delete(mapped);
                cpl_image_delete(sky_mapped);
                cpl_image_delete(sci_sky_mapped);
                cpl_table_delete(slits);
                cpl_image_delete(mapped_sky_err);
                return -1;
            }
            cpl_image_delete(mapped_sky_err);
        }

        cpl_image_delete(sci_sky_mapped);

        //Determine wheter the flux calibration will use the flat correction
        bool photcal_apply_flat_corr = vimos_science_photcal_apply_flat_corr_multiframe
                (frameset, specphot_tag, specphot_tag,
                 flat_sed_tag, &photometry);

        //TODO: Place this in a better place. The whole recipe needs refactoring...
        cpl_image        * mapped_flat_sed = NULL;
        cpl_propertylist * flat_sed_header = NULL;
        vimos::detected_slits det_slits = 
                vimos::detected_slits_from_tables(slits, polytraces, nx);
        mosca::wavelength_calibration wave_cal(idscoeff, refwave);
        if(photcal_apply_flat_corr)
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

        if (photometry) {
            cpl_image *calibrated;
            cpl_table *ext_table;
            cpl_table *response_interp;
            cpl_propertylist * specphot_header;

            ext_table = dfs_load_table(frameset, "EXTINCT_TABLE", 1);

            response_interp = dfs_load_table(frameset, specphot_tag, 2);

            specphot_header = dfs_load_header(frameset, specphot_tag, 0);

            cpl_image * mapped_images = cpl_image_duplicate(mapped);

            if(photcal_apply_flat_corr)
            {
                cpl_msg_info(cpl_func, "Applying flat SED correction");
                if(vimos_science_correct_flat_sed_mapped(mapped_images, slits,
                        mapped_flat_sed, flat_sed_header, specphot_header,
                        det_slits))
                {
                    return -1;
                }
            }

            calibrated = mos_apply_photometry(mapped_images, response_interp,
                                              ext_table, startwavelength,
                                              dispersion, gain, alltime / mos,
                                              airmass);
            cpl_table_delete(ext_table);
            cpl_table_delete(response_interp);

            cpl_propertylist_update_string(header, "BUNIT",
                                   "10^(-16) erg/(cm^2 s Angstrom)");

            if (dfs_save_image(frameset, calibrated,
                               mapped_flux_science_tag, header,
                               parlist, recipe, version)) {
                cpl_image_delete(calibrated);
                cpl_propertylist_delete(header);
                return -1;
            }

            cpl_image_delete(calibrated);
            cpl_image_delete(mapped_images);
        }

        cpl_propertylist_update_string(header, "BUNIT", "ADU");

        if (dfs_save_image(frameset, sky_mapped, mapped_sky_tag,
                           header, parlist, recipe, version)) {
            cpl_propertylist_delete(header);
            cpl_image_delete(mapped);
            cpl_image_delete(sky_mapped);
            cpl_table_delete(slits);
            return -1;
        }

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Final object detection...");
        cpl_msg_indent_more();

        if (cosmics || strcmp(stack_method, "average")) {
            image = mos_detect_objects(mapped, slits, slit_margin, ext_radius,
                                       cont_radius, detection);
        }
        else {
            cpl_image *mapped_cleaned = cpl_image_duplicate(mapped);
            mos_clean_cosmics(mapped_cleaned, gain, -1., -1.);
            image = mos_detect_objects(mapped_cleaned, slits, slit_margin,
                                       ext_radius, cont_radius, detection);

            cpl_image_delete(mapped_cleaned);
        }

        cpl_image_delete(image);

        mos_rotate_slits(slits, rotate, ny, nx);
        if (dfs_save_table(frameset, slits, object_table_tag, NULL, parlist,
                           recipe, version)) {
            cpl_propertylist_delete(header);
            cpl_image_delete(mapped);
            cpl_image_delete(sky_mapped);
            cpl_table_delete(slits);
            return -1;
        }

        cpl_msg_indent_less();
        cpl_msg_info(recipe, "Final object extraction...");
        cpl_msg_indent_more();

        images = mos_extract_objects(mapped, mapped_var, sky_mapped, slits,
                                     ext_mode, ron, gain, mos);

        cpl_image_delete(mapped);
        cpl_image_delete(sky_mapped);

        if (images) {

            if (photometry) {
                cpl_image *calibrated;
                cpl_table *ext_table;
                cpl_table *response_interp;
                cpl_propertylist * specphot_header;

                ext_table = dfs_load_table(frameset, "EXTINCT_TABLE", 1);

                response_interp = dfs_load_table(frameset, specphot_tag, 2);

                specphot_header = dfs_load_header(frameset, specphot_tag, 0);

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
                                                  dispersion, gain, 
                                                  alltime / mos, airmass);

                cpl_table_delete(ext_table);
                cpl_table_delete(response_interp);

                cpl_propertylist_update_string(header, "BUNIT",
                                   "10^(-16) erg/(cm^2 s Angstrom)");

                if (dfs_save_image(frameset, calibrated,
                                   reduced_flux_science_tag, header,
                                   parlist, recipe, version)) {
                    cpl_image_delete(calibrated);
                    cpl_image_delete(images[0]);
                    cpl_image_delete(images[1]);
                    cpl_image_delete(images[2]);
                    cpl_free(images);
                    cpl_propertylist_delete(header);
                    return -1;
                }

                cpl_image_delete(calibrated); calibrated = NULL;
                cpl_image_delete(science_images); science_images = NULL;
            }

            if (time_normal) {
                cpl_propertylist_update_string(header, "BUNIT", "ADU/s");
                cpl_image_divide_scalar(images[0], alltime / mos);
            }
            else {
                cpl_propertylist_update_string(header, "BUNIT", "ADU");
            }

            if (dfs_save_image(frameset, images[0], reduced_science_tag, 
                               header, parlist, recipe, version)) {
                cpl_image_delete(images[0]);
                cpl_image_delete(images[1]);
                cpl_image_delete(images[2]);
                cpl_free(images);
                cpl_propertylist_delete(header);
                return -1;
            }


            if (time_normal) {
                cpl_image_divide_scalar(images[1], alltime / mos);
            }

            if (dfs_save_image(frameset, images[1], reduced_sky_tag, header,
                               parlist, recipe, version)) {
                cpl_image_delete(images[1]);
                cpl_image_delete(images[2]);
                cpl_free(images);
                cpl_propertylist_delete(header);
                return -1;
            }

            cpl_image_delete(images[1]);

            if (photometry) {
                cpl_image *calibrated;
                cpl_table *ext_table;
                cpl_table *response_interp;
                cpl_propertylist * specphot_header;

                ext_table = dfs_load_table(frameset, "EXTINCT_TABLE", 1);

                response_interp = dfs_load_table(frameset, specphot_tag, 2);
                specphot_header = dfs_load_header(frameset, specphot_tag, 0);

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
                                                  dispersion, gain, alltime / mos,
                                                  airmass);

                cpl_table_delete(ext_table);
                cpl_table_delete(response_interp);

                cpl_propertylist_update_string(header, "BUNIT",
                                   "10^(-16) erg/(cm^2 s Angstrom)");

                if (dfs_save_image(frameset, calibrated,
                                   reduced_flux_error_tag, header,
                                   parlist, recipe, version)) {
                    cpl_image_delete(calibrated);
                    cpl_image_delete(images[0]);
                    cpl_image_delete(images[2]);
                    cpl_free(images);
                    cpl_propertylist_delete(header);
                    return -1;
                }
                cpl_image_delete(calibrated);
            }

            if (time_normal) {
                cpl_propertylist_update_string(header, "BUNIT", "ADU/s");
                cpl_image_divide_scalar(images[2], alltime / mos);
            }
            else {
                cpl_propertylist_update_string(header, "BUNIT", "ADU");
            }

            if (dfs_save_image(frameset, images[2], reduced_error_tag, header,
                               parlist, recipe, version)) {
                cpl_image_delete(images[2]);
                cpl_free(images);
                cpl_propertylist_delete(header);
                return -1;
            }

            cpl_image_delete(images[2]);
            cpl_image_delete(images[0]);

            cpl_free(images);
            cpl_table_delete(slits);
            cpl_table_delete(polytraces);
            cpl_table_delete(idscoeff);
        }
        else {
            cpl_msg_warning(recipe, "No objects found: the products "
                            "%s, %s, and %s will be empty.",
                            reduced_science_tag, reduced_sky_tag,
                            reduced_error_tag);

            if (dfs_save_image_null(frameset, NULL, parlist, reduced_science_tag, 
                                    recipe, version))
                return -1;

            if (dfs_save_image_null(frameset, NULL, parlist, reduced_sky_tag, 
                                    recipe, version))
                return -1;

            if (dfs_save_image_null(frameset, NULL, parlist, reduced_error_tag,
                                    recipe, version))
                return -1;
        }

        cpl_propertylist_delete(header);

        skyalign = dfs_get_parameter_int(parlist,
                                         "vimos.vmmosscience.skyalign", NULL);
        if(skyalign >= 0)
          system("rm object_table_*.fits "
                 "mos_unmapped_science_*.fits mos_unmapped_sky_*.fits "
                 "mos_sci_skylines_offsets_slit.fits "
                 "mos_sci_wavelength_map_sky.fits mos_sci_disp_coeff_sky.fits");
        else
          system("rm object_table_*.fits "
                 "mos_unmapped_science_*.fits mos_unmapped_sky_*.fits ");

        return 0;
    }

    return vimos_science_impl(frameset, parlist, section);
}

static bool vimos_science_photcal_apply_flat_corr_multiframe
(cpl_frameset * frameset, const char * specphot_tag, 
 const char * master_specphot_tag, const char * flat_sed_tag,
 int* fluxcal)
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
    return false;
}
