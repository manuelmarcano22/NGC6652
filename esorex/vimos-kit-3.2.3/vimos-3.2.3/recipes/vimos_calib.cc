/* $Id: vimos_calib.c,v 1.16 2013-10-22 16:55:57 cgarcia Exp $
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
 * $Date: 2013-10-22 16:55:57 $
 * $Revision: 1.16 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <vimos_calib_impl.h>
#include <vimos_calib_mult.h>
#include <cpl.h>
#include <moses.h>
#include <vimos_dfs.h>
#include <vmutils.h>
#include <exception>

#define vimos_multi_exit(message)                    \
{                                                    \
if (message) cpl_msg_error("vmmoscalib", message);  \
cpl_propertylist_delete(header);                     \
cpl_table_delete(maskslits);                         \
cpl_table_delete(slits);                             \
cpl_msg_indent_less();                               \
return -1;                                           \
}

static int vimos_calib_create(cpl_plugin *);
static int vimos_calib_exec(cpl_plugin *);
static int vimos_calib_destroy(cpl_plugin *);
static int vimos_calib(cpl_parameterlist *, cpl_frameset *);

static char vimos_calib_description[] =
"This recipe is used to identify reference lines on MOS arc lamp\n"
"exposures, and trace the spectral edges on the corresponding flat field\n"
"exposures. This information is used to determine the spectral extraction\n"
"mask to be applied in the scientific data reduction, performed with the\n"
"recipe vimos_science. The input arc lamp and flat field exposures are\n"
"assumed to be obtained quasi-simultaneously, so that they would be\n" 
"described by exactly the same instrument distortions.\n"
"A line catalog must be specified, containing the wavelengths of the\n"
"reference arc lamp lines used for the wavelength calibration. A grism\n"
"table (typically depending on the instrument mode, and in particular on\n"
"the grism used) may also be specified: this table contains a default\n"
"recipe parameter setting to control the way spectra are extracted for\n"
"a specific instrument mode, as it is used for automatic run of the\n"
"pipeline on Paranal and in Garching. If this table is specified, it\n"
"will modify the default recipe parameter setting, with the exception of\n"
"those parameters which have been explicitly modifyed on the command line.\n"
"If a grism table is not specified, the input recipe parameters values\n"
"will always be read from the command line, or from an esorex configuration\n"
"file if present, or from their generic default values (that are rarely\n"
"meaningful). Finally a master bias frame must be input to this recipe.\n" 
"In the table below the MOS_CURV_COEFF, MOS_CURV_TRACES, MOS_SPATIAL_MAP\n"
"MOS_ARC_SPECTRUM_EXTRACTED, MOS_SPECTRA_DETECTION, MOS_SLIT_MAP, and\n" 
"MOS_SLIT_LOCATION, are never created in case of long-slit-like data.\n" 
"The products MOS_SPECTRA_DETECTION, MOS_SLIT_MAP, and MOS_DISP_RESIDUALS,\n" 
"are just created if the --check parameter is set to true. The product\n"
"GLOBAL_DISTORTION_TABLE is just created if more than 12 separate spectra\n"
"are found in the CCD.\n\n"
"Input files:\n\n"
"  DO category:               Type:       Explanation:         Required:\n"
"  MOS_SCREEN_FLAT            Raw         Flat field exposures    Y\n"
"  MOS_ARC_SPECTRUM           Raw         Arc lamp exposure       Y\n"
"  MASTER_BIAS                Calib       Bias frame              Y\n"
"  LINE_CATALOG               Calib       Line catalog            Y\n"
"  CONFIG_TABLE               Calib       Grism table             .\n\n"
"Output files:\n\n"
"  DO category:               Data type:  Explanation:\n"
"  MOS_COMBINED_SCREEN_FLAT   FITS image  Combined (sum) flat field\n"
"  MOS_MASTER_SCREEN_FLAT     FITS image  Normalised flat field\n"
"  MOS_ARC_SPECTRUM_EXTRACTED FITS image  Wavelength calibrated arc spectrum\n"
"  MOS_DISP_COEFF             FITS table  Inverse dispersion coefficients\n"
"  MOS_DISP_RESIDUALS         FITS image  Residuals in wavelength calibration\n"
"  MOS_DISP_RESIDUALS_TABLE   FITS table  Residuals in wavelength calibration\n"
"  MOS_DELTA_IMAGE            FITS image  Offset vs linear wavelength calib\n"
"  MOS_WAVELENGTH_MAP         FITS image  Wavelength for each pixel on CCD\n"
"  MOS_SPECTRA_DETECTION      FITS image  Check for preliminary detection\n"
"  MOS_SLIT_MAP               FITS image  Map of central wavelength on CCD\n"
"  MOS_CURV_TRACES            FITS table  Spectral curvature traces\n"
"  MOS_CURV_COEFF             FITS table  Spectral curvature coefficients\n"
"  MOS_SPATIAL_MAP            FITS image  Spatial position along slit on CCD\n"
"  MOS_SPECTRAL_RESOLUTION    FITS table  Resolution at reference arc lines\n"
"  MOS_SLIT_LOCATION          FITS table  Slits on product frames and CCD\n"
"  GLOBAL_DISTORTION_TABLE    FITS table  Global distortions table\n\n";

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
    cpl_recipe *recipe = static_cast<cpl_recipe *>(cpl_calloc(1, sizeof *recipe ));
    cpl_plugin *plugin = &recipe->interface;

    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmmoscalib",
                    "Determination of the extraction mask",
                    vimos_calib_description,
                    "Carlo Izzo",
                    PACKAGE_BUGREPORT,
    "This file is currently part of the VIMOS Instrument Pipeline\n"
    "Copyright (C) 2002-2006 European Southern Observatory\n\n"
    "This program is free software; you can redistribute it and/or modify\n"
    "it under the terms of the GNU General Public License as published by\n"
    "the Free Software Foundation; either version 2 of the License, or\n"
    "(at your option) any later version.\n\n"
    "This program is distributed in the hope that it will be useful,\n"
    "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
    "GNU General Public License for more details.\n\n"
    "You should have received a copy of the GNU General Public License\n"
    "along with this program; if not, write to the Free Software Foundation,\n"
    "Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA\n",
                    vimos_calib_create,
                    vimos_calib_exec,
                    vimos_calib_destroy);

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

static int vimos_calib_create(cpl_plugin *plugin)
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
     * Dispersion
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.dispersion",
                                CPL_TYPE_DOUBLE,
                                "Expected spectral dispersion (Angstrom/pixel)",
                                "vimos.vmmoscalib",
                                0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "dispersion");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Peak detection level
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.peakdetection",
                                CPL_TYPE_DOUBLE,
                                "Initial peak detection threshold (ADU)",
                                "vimos.vmmoscalib",
                                0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "peakdetection");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* 
     * Degree of wavelength calibration polynomial
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.wdegree",
                                CPL_TYPE_INT,
                                "Degree of wavelength calibration polynomial",
                                "vimos.vmmoscalib",
                                0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wdegree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Reference lines search radius
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.wradius",
                                CPL_TYPE_INT,
                                "Search radius if iterating pattern-matching "
                                "with first-guess method (pixel)",
                                "vimos.vmmoscalib",
                                4);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wradius");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Rejection threshold in dispersion relation polynomial fitting
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.wreject",
                                CPL_TYPE_DOUBLE,
                                "Rejection threshold in dispersion "
                                "relation fit (pixel)",
                                "vimos.vmmoscalib",
                                1.);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wreject");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Wavelength solution interpolation (for LSS data)
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.wmodelss",
                                CPL_TYPE_INT,
                                "Interpolation mode of wavelength solution "
                                "(0 = no interpolation, 1 = fill gaps, "
                                "2 = global model)",
                                "vimos.vmmoscalib",
                                2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wmodelss");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    /*
     * Wavelength solution interpolation (for MOS data)
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.wmodemos",
                                CPL_TYPE_INT,
                                "Interpolation mode of wavelength solution "
                                "(0 = no interpolation, 1 = local (slit) "
                                "solution, 2 = global model)",
                                "vimos.vmmoscalib",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wmodemos");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    /*
     * Catalog lines to ignore in wavelength calibration
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.ignore_lines",
                                CPL_TYPE_STRING,
                                "Catalog lines nearest to wavelengths in this "
                                "list will be ignored for wavelength calibration",
                                "vimos.vmmoscalib",
                                "");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ignore_lines");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Linesets to use
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.used_linesets",
                                CPL_TYPE_STRING,
                                "Linesets to use. Valid are 'standard' and"
                                "'extended' (see column LINE_SET in the "
                                "line catalogue)",
                                "vimos.vmmoscalib",
                                "standard");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "used_linesets");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Degree of spectral curvature polynomial
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.cdegree",
                                CPL_TYPE_INT,
                                "Degree of spectral curvature polynomial",
                                "vimos.vmmoscalib",
                                0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "cdegree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Curvature solution interpolation (for MOS-like data)
     */
 
    p = cpl_parameter_new_value("vimos.vmmoscalib.cmode",
                                CPL_TYPE_INT,
                                "Interpolation mode of curvature solution "
                                "applicable to MOS-like data (0 = no "
                                "interpolation, 1 = fill gaps, 2 = global "
                                "model)",
                                "vimos.vmmoscalib",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "cmode");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Start wavelength for spectral extraction
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.startwavelength",
                                CPL_TYPE_DOUBLE,
                                "Start wavelength in spectral extraction",
                                "vimos.vmmoscalib",
                                0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "startwavelength");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * End wavelength for spectral extraction
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.endwavelength",
                                CPL_TYPE_DOUBLE,
                                "End wavelength in spectral extraction",
                                "vimos.vmmoscalib",
                                0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "endwavelength");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Reference wavelength for wavelength calibration
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.reference",
                                CPL_TYPE_DOUBLE,
                                "Reference wavelength for slit map",
                                "vimos.vmmoscalib",
                                0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "reference");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Try slit identification
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.slit_ident",
                                CPL_TYPE_BOOL,
                                "Attempt slit identification."
                                "For multiplexing data slit identification "
                                "is always performed. Switching it off will "
                                "cause systematic errors in the absolute flux "
                                "calibration if the FLAT_SED correction is "
                                "used (see pipeline manual for details)",
                                "vimos.vmmoscalib",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "slit_ident");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Polynomial degree for flat field polynomial fitting along spatial direction 
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.s_degree",
                                CPL_TYPE_INT,
                                "Polynomial degree for the flat field fitting "
                                "along spatial direction",
                                "vimos.vmmoscalib",
                                -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "s_degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Smooth box radius for flat field along spatial direction
     * (if s_degree < 0)
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.sradius",
                                CPL_TYPE_INT,
                                "Smooth box radius for flat field along "
                                "spatial direction (used if s_degree < 0)",
                                "vimos.vmmoscalib",
                                -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "sradius");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Number of knots in flat field fitting splines along dispersion direction 
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.d_nknots",
                                CPL_TYPE_INT,
                                "Number of knots in flat field fitting "
                                "splines along dispersion direction",
                                "vimos.vmmoscalib",
                                -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "d_nknots");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Smooth box radius for flat field along dispersion direction
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.dradius",
                                CPL_TYPE_INT,
                                "Smooth box radius for flat field along "
                                "dispersion direction (if d_knots < 0)",
                                "vimos.vmmoscalib",
                                10);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "dradius");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Threshold percentage for flat spline fitting with respect to the maximum
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.fit_threshold",
                                CPL_TYPE_DOUBLE,
                                "Threshold percentage for flat spline fitting "
                                "with respect to the maximum",
                                "vimos.vmmoscalib",
                                0.01);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "fit_threshold");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
     * Tolerance used for ratios during line identification 
     */

    p = cpl_parameter_new_value("vimos.vmmoscalib.line_ident_tol",
                                CPL_TYPE_DOUBLE,
                                "Tolerance for the ratio of detected lines "
                                "vs reference lines. This is used during "
                                "for arc line identification.",
                                "vimos.vmmoscalib",
                                0.05);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "line_ident_tol");
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

static int vimos_calib_exec(cpl_plugin *plugin)
{
    cpl_recipe  *   recipe ;
    int             status = 1;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin ;
    else return -1 ;

    /* Issue a banner */
    vimos_print_banner();

    try
    {
        status = vimos_calib(recipe->parameters, recipe->frames);
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

static int vimos_calib_destroy(cpl_plugin *plugin)
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

static int vimos_calib(cpl_parameterlist *parlist, cpl_frameset *frameset)
{
    const char *recipe = "vmmoscalib";

    cpl_table        *maskslits = NULL;
    cpl_table        *slits     = NULL;
    cpl_propertylist *header    = NULL;
    int               mos;
    int               multiplex;
    int               status;


    /*
     * Check whether there is spectral multiplexing here
     */

    mos = cpl_frameset_count_tags(frameset, "MOS_ARC_SPECTRUM");

    if (mos == 0)
        vimos_multi_exit("Missing input arc lamp frame");

    if (mos > 1)
        vimos_multi_exit("Just one input arc lamp frame is allowed");

    header = dfs_load_header(frameset, "MOS_ARC_SPECTRUM", 0);

    if (header == NULL)
        vimos_multi_exit("Cannot load arc lamp header");

    maskslits = mos_load_slits_vimos(header, 0);
    cpl_propertylist_delete(header); header = NULL;
    multiplex = mos_check_multiplex(maskslits);
//    cpl_table_save(maskslits, NULL, NULL, "groups.fits", CPL_IO_DEFAULT);

    if (cpl_error_get_code())
        vimos_multi_exit("Error retrieving slit information from header");
    
    if (multiplex == 1) {

        /*
         * No spectral multiplexing: classical handling
         */

        cpl_table_delete(maskslits); maskslits = NULL;

        return vimos_calib_impl(frameset, parlist);
    }

    /*
     * Here we deal with spectral multiplexing
     */

    if (0) {   // protect the rest
        cpl_msg_info(recipe, "Spectral multiplexing: %d", multiplex);
        vimos_multi_exit("Cannot handle such data");
    }

    mos_assign_multiplex_group(maskslits);

//    cpl_table_save(maskslits, NULL, NULL, "groups.fits", CPL_IO_DEFAULT);

    status = vimos_calib_mult(frameset, parlist, maskslits);

    cpl_table_delete(maskslits); maskslits = NULL;

    return status;
}
