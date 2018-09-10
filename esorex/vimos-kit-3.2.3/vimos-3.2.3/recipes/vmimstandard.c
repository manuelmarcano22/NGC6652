/* $Id: vmimstandard.c,v 1.13 2012-01-26 16:07:06 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2005 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2012-01-26 16:07:06 $
 * $Revision: 1.13 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <string.h>

#include <cxmemory.h>
#include <cxstring.h>

#include <cpl.h>

#include <pilmemory.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilrecipe.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilqc.h>
#include <pilfits.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmgalaxytable.h"
#include "vmstarmatchtable.h"
#include "vmastrometrictable.h"
#include "vmmath.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgastrometry.h"
#include "vmimgextraction.h"
#include "vmimgphotcalib.h"
#include "vmimgutils.h"
#include "vmqcutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


static int vmimstandard_create(cpl_plugin *);
static int vmimstandard_exec(cpl_plugin *);
static int vmimstandard_destroy(cpl_plugin *);

static cxint vmimstandard(PilSetOfFrames *);

static char vmimstandard_description[] =
"This recipe is used to determine the instrumental magnitude of the\n"
"stars matching the entries of a photometric catalog. The same procedure\n"
"applied by the recipe vmimobsstare is used to reduce the standard field\n"
"image. SExtractor v2.1.6 is used for the source detection task.\n\n"
"Input files:\n\n"
"  DO category:              Type:       Explanation:         Required:\n"
"  IMG_STANDARD              Raw         Standard stars image    Y\n"
"  MASTER_BIAS               Calib       Master bias             Y\n"
"  MASTER_DARK               Calib       Master dark             .\n"
"  IMG_MASTER_SKY_FLAT       Calib       Master sky flat         Y\n"
"  PHOTOMETRIC_CATALOG       Calib       Photometric catalog     Y\n"
"  PHOTOMETRIC_TABLE         Calib       Photometric table       Y\n"
"  CCD_TABLE                 Calib       Bad pixel table         .\n\n"
"Output files:\n\n"
"  DO category:              Data type:  Explanation:\n"
"  IMG_STAR_MATCH_TABLE      FITS table  List of found standard stars\n"
"  IMG_STANDARD_REDUCED      FITS image  Reduced scientific exposure\n"
"  IMG_GALAXY_TABLE          FITS table  List of detected objects\n\n"
"The galaxy table is the output of SExtractor converted into FITS\n"
"format. The star match table is the list of identified standard stars,\n"
"with their positions on sky and CCD and their instrumental and catalog\n"
"magnitudes.\n\n"
"For more details, please refer to the VIMOS Pipeline User's Guide.";


/**
 * Definition of supported recipe engine modes.
 */

enum EngineMode {
  MODE_UNDEF = 0,
  MODE_SCIENCE,
  MODE_PHOTOMETRY
};


/**
 * Definition of the label strings for all methods the recipe function
 * supports for removing the bias.
 */

static const char *biasMethodNames[] = {
  "Master",
  "Zmaster"
};

static const BiasMethod biasMethods[] = {
  BIAS_MASTER,
  BIAS_ZMASTER
};

static unsigned int nBiasMethods = sizeof(biasMethods) / sizeof(BiasMethod);


/**
 * Build table of contents, i.e. the list of available plugins, for
 * this module. This function is exported.
 */

int cpl_plugin_get_info(cpl_pluginlist *list)
{

    cpl_recipe *recipe = cpl_calloc(1, sizeof(*recipe));
    cpl_plugin *plugin = &recipe->interface;


    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmimstandard",
                    "Reduce an imaging standard star exposure.",
                    vmimstandard_description,
                    "ESO VIMOS Pipeline Team and VIMOS Consortium",
                    PACKAGE_BUGREPORT,
                    vimos_get_license(),
                    vmimstandard_create,
                    vmimstandard_exec,
                    vmimstandard_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}


/**
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static int vmimstandard_create(cpl_plugin *plugin)
{
    cpl_recipe    *recipe;
    cpl_parameter *p;

    cx_string *path = cx_string_new();
    cxint status = 0;

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
    if (recipe->parameters == NULL)
        return 1;

    /*
     * Fill the parameter list
     */


    p = cpl_parameter_new_enum("vimos.Parameters.bias.removing.method",
                                CPL_TYPE_STRING,
                                "Bias removal method.",
                                "vimos.Parameters",
                                "Zmaster", 2, "Zmaster", "Master");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "BiasMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "BiasMethod");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.badpixel.clean",
                                CPL_TYPE_BOOL,
                                "Bad pixel correction on reduced science "
                                "image.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.clean",
                                CPL_TYPE_BOOL,
                                "Cosmic ray removal from input science image.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanCosmic");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanCosmic");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.threshold",
                                CPL_TYPE_DOUBLE,
                                "Threshold for cosmic ray candidate "
                                "selection.",
                                "vimos.Parameters",
                                4.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CosmicsThreshold");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CosmicsThreshold");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.ratio",
                                CPL_TYPE_DOUBLE,
                                "Ratio for discrimination of objects and "
                                "cosmics.",
                                "vimos.Parameters",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CosmicsRatio");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CosmicsRatio");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.magnitude.limit",
                                CPL_TYPE_DOUBLE,
                                "Magnitude upper limit for selection of "
                                "detected objects.",
                                "vimos.Parameters",
                                100.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MagLimit");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MagLimit");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.stellarity",
                                CPL_TYPE_DOUBLE,
                                "Stellarity index used to select stars "
                                "in the image.",
                                "vimos.Parameters",
                                0.01);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "StarIndex");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "StarIndex");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.stars",
                                CPL_TYPE_INT,
                                "Minimum number of stars required.",
                                "vimos.Parameters",
                                3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MinStars");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MinStars");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.radius",
                                CPL_TYPE_DOUBLE,
                                "Aperture used for object identification.",
                                "vimos.Parameters",
                                5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SearchRadius");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SearchRadius");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.magnitude.initial",
                                CPL_TYPE_DOUBLE,
                                "Magnitude tolerance for initial selection "
                                "of matching stars.",
                                "vimos.Parameters",
                                2.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MagInitial");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MagInitial");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.magnitude.final",
                                CPL_TYPE_DOUBLE,
                                "Magnitude tolerance for final selection of "
                                "matching stars.",
                                "vimos.Parameters",
                                1.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MagFinal");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MagFinal");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.ksigma",
                                CPL_TYPE_DOUBLE,
                                "Sigma clipping factor used when matching "
                                "stars.",
                                "vimos.Parameters",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "KSigmaClip");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "KSigmaClip");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.temperature.check",
                                CPL_TYPE_BOOL,
                                "Check beam temperature when updating "
                                "the world coordinate system.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "TemperatureCheck");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "TemperatureCheck");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.temperature.tolerance",
                                CPL_TYPE_DOUBLE,
                                "Maximum allowed difference between beam "
                                "and ambient temperature.",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "TemperatureTolerance");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "TemperatureTolerance");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.catalog.remote",
                                CPL_TYPE_BOOL,
                                "Enable online access to a catalog server. "
                                "(Not yet implemented!)",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "RemoteCatalog");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "RemoteCatalog");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.quality.enable",
                                CPL_TYPE_BOOL,
                                "Compute QC1 parameters",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ComputeQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ComputeQC");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.comatrix.radius",
                                CPL_TYPE_DOUBLE,
                                "Aperture used for object identification "
                                "during computation of the CO matrix",
                                "vimos.Parameters",
                                6.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "COMatrixSearchRadius");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "COMatrixSearchRadius");
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.Parameters.reduceall",
                                CPL_TYPE_BOOL,
                                "Reduce any frame.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ReduceAnyFrame");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ReduceAnyFrame");
    cpl_parameterlist_append(recipe->parameters, p);

#else

    p = cpl_parameter_new_value("vimos.Parameters.reduceall",
                                CPL_TYPE_BOOL,
                                "Reduce any frame.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ReduceAnyFrame");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ReduceAnyFrame");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

    p = cpl_parameter_new_value("vimos.SExtractor.AnalysisThresh",
                                CPL_TYPE_DOUBLE,
                                "Surface brightness threshold for "
                                "FWHM computation.",
                                "vimos.SExtractor",
                                2.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.AnalysisThresh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.AnalysisThresh");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.AssocData",
                                CPL_TYPE_STRING,
                                "Column indices for columns to copy "
                                "to the catalog output.",
                                "vimos.SExtractor",
                                "2,3,4");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.AssocData");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.AssocData");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.AssocName",
                                CPL_TYPE_STRING,
                                "Name of the ASSOC file.",
                                "vimos.SExtractor",
                                "sky.list");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.AssocName");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.AssocName");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.AssocParams",
                                CPL_TYPE_STRING,
                                "Column indices in the ASSOC file to use as "
                                "coordinates and weight for cross-matching.",
                                "vimos.SExtractor",
                                "2,3,4");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.AssocParams");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.AssocParams");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.AssocRadius",
                                CPL_TYPE_DOUBLE,
                                "Search radius for ASSOC.",
                                "vimos.SExtractor",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.AssocRadius");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.AssocRadius");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.AssocSelecType",
                                CPL_TYPE_STRING,
                                "Output selector for detected sources.",
                                "vimos.SExtractor",
                                "MATCHED");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.AssocSelecType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.AssocSelecType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.AssocType",
                                CPL_TYPE_STRING,
                                "Method for cross-matching.",
                                "vimos.SExtractor",
                                "MAG_SUM");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.AssocType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.AssocType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.BackFilterSize",
                                CPL_TYPE_INT,
                                "Size (in background meshes) of the "
                                "background filtering mask.",
                                "vimos.SExtractor",
                                3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.BackFilterSize");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.BackFilterSize");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.BackPhotoThick",
                                CPL_TYPE_INT,
                                "Thickness of the background LOCAL annulus.",
                                "vimos.SExtractor",
                                24);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.BackPhotoThick");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.BackPhotoThick");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.BackPhotoType",
                                CPL_TYPE_STRING,
                                "Select background for magnitude computation.",
                                "vimos.SExtractor",
                                "GLOBAL");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.BackPhotoType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.BackPhotoType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.BackSize",
                                CPL_TYPE_INT,
                                "Size of a background mesh.",
                                "vimos.SExtractor",
                                64);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.BackSize");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.BackSize");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.BackValue",
                                CPL_TYPE_DOUBLE,
                                "Constant to subtract from the images "
                                "for MANUAL background subtraction.",
                                "vimos.SExtractor",
                                0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.BackValue");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.BackValue");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.CatalogType",
                                CPL_TYPE_STRING,
                                "Select output catalog format (only "
                                "ASCII_HEAD is supported).",
                                "vimos.SExtractor",
                                "ASCII_HEAD");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.CatalogType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.CatalogType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.CheckImageName",
                                CPL_TYPE_STRING,
                                "Check image file name.",
                                "vimos.SExtractor",
                                "check.fits");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.CheckImageName");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.CheckImageName");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.CheckImageType",
                                CPL_TYPE_STRING,
                                "Select information to put in "
                                "the `check image'.",
                                "vimos.SExtractor",
                                "NONE");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.CheckImageType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.CheckImageType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.Clean",
                                CPL_TYPE_STRING,
                                "If `Y', clean catalog before writing.",
                                "vimos.SExtractor",
                                "Y");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.Clean");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.Clean");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.CleanParam",
                                CPL_TYPE_DOUBLE,
                                "Cleaning efficiency.",
                                "vimos.SExtractor",
                                1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.CleanParam");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.CleanParam");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.DeblendMinCont",
                                CPL_TYPE_DOUBLE,
                                "Minimum contrast for deblending.",
                                "vimos.SExtractor",
                                0.005);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.DeblendMinCont");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.DeblendMinCont");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.DeblendNthresh",
                                CPL_TYPE_INT,
                                "Number of deblending sub-thresholds.",
                                "vimos.SExtractor",
                                32);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.DeblendNthresh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.DeblendNthresh");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.DetectMinArea",
                                CPL_TYPE_INT,
                                "Minimum number of pixels above threshold "
                                "triggering detection.",
                                "vimos.SExtractor",
                                5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.DetectMinArea");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.DetectMinArea");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.DetectThresh",
                                CPL_TYPE_DOUBLE,
                                "Detection threshold (relative to "
                                "background RMS).",
                                "vimos.SExtractor",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.DetectThresh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.DetectThresh");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.DetectType",
                                CPL_TYPE_STRING,
                                "Device type the image originates from.",
                                "vimos.SExtractor",
                                "CCD");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.DetectType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.DetectType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.Filter",
                                CPL_TYPE_STRING,
                                "If `Y', filter data before extraction.",
                                "vimos.SExtractor",
                                "N");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.Filter");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.Filter");
    cpl_parameterlist_append(recipe->parameters, p);

    cx_string_set(path, VIMOS_SEXTRACTOR_CONFIG);
    cx_string_append(path, "/vimos.conv");

    p = cpl_parameter_new_value("vimos.SExtractor.FilterName",
                                CPL_TYPE_STRING,
                                "Name of the file containing the filter "
                                "definition.",
                                "vimos.SExtractor",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.FilterName");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.FilterName");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.FilterThresh",
                                CPL_TYPE_STRING,
                                "Lower, upper threshold (in background "
                                "sigmas) for filtering (retina-filtering "
                                "only).",
                                "vimos.SExtractor",
                                "");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.FilterThresh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.FilterThresh");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.FlagImage",
                                CPL_TYPE_STRING,
                                "Flag image file name.",
                                "vimos.SExtractor",
                                "flag.fits");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.FlagImage");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.FlagImage");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.FlagType",
                                CPL_TYPE_STRING,
                                "Flag combination method.",
                                "vimos.SExtractor",
                                "OR");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.FlagType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.FlagType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.InterpMaxXlag",
                                CPL_TYPE_INT,
                                "Maximum X gap allowed in interpolation.",
                                "vimos.SExtractor",
                                16);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.InterpMaxXlag");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.InterpMaxXlag");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.InterpMaxYlag",
                                CPL_TYPE_INT,
                                "Maximum Y gap allowed in interpolation.",
                                "vimos.SExtractor",
                                16);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.InterpMaxYlag");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.InterpMaxYlag");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.InterpType",
                                CPL_TYPE_STRING,
                                "Interpolation method.",
                                "vimos.SExtractor",
                                "ALL");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.InterpType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.InterpType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.MagGamma",
                                CPL_TYPE_DOUBLE,
                                "Gamma of emulsion (only used in PHOTO mode).",
                                "vimos.SExtractor",
                                4.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.MagGamma");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.MagGamma");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.MagZeropoint",
                                CPL_TYPE_DOUBLE,
                                "Zero-point offset to apply to magnitudes.",
                                "vimos.SExtractor",
                                0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.MagZeropoint");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.MagZeropoint");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.MaskType",
                                CPL_TYPE_STRING,
                                "Masking of neighbours for photometry.",
                                "vimos.SExtractor",
                                "CORRECT");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.MaskType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.MaskType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.MemoryBufSize",
                                CPL_TYPE_INT,
                                "Number of scan-lines in the image buffer.",
                                "vimos.SExtractor",
                                1024);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.MemoryBufSize");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.MemoryBufSize");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.MemoryObjStack",
                                CPL_TYPE_INT,
                                "Maximum number of objects the object "
                                "stack can contain.",
                                "vimos.SExtractor",
                                2000);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.MemoryObjStack");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.MemoryObjStack");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.MemoryPixStack",
                                CPL_TYPE_INT,
                                "Maximum number of pixels the pixel "
                                "stack can contain.",
                                "vimos.SExtractor",
                                100000);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.MemoryPixStack");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.MemoryPixStack");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.PhotApertures",
                                CPL_TYPE_DOUBLE,
                                "Aperture diameters used for MAG_APER.",
                                "vimos.SExtractor",
                                5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.PhotApertures");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.PhotApertures");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.PhotAutoParams",
                                CPL_TYPE_STRING,
                                "MAG_AUTO controls: 1st order moment "
                                "scaling parameter, minimum Rmin.",
                                "vimos.SExtractor",
                                "2.5,3.5");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.PhotAutoParams");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.PhotAutoParams");
    cpl_parameterlist_append(recipe->parameters, p);

    cx_string_set(path, VIMOS_SEXTRACTOR_PATH);
    cx_string_append(path, "/sex");

    p = cpl_parameter_new_value("vimos.SExtractor.SExtractor",
                                CPL_TYPE_STRING,
                                "Path to the sextractor executable.",
                                "vimos.SExtractor",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.SExtractor");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.SExtractor");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.SaturLevel",
                                CPL_TYPE_DOUBLE,
                                "Pixel values above are considered saturated.",
                                "vimos.SExtractor",
                                60000.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.SaturLevel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.SaturLevel");
    cpl_parameterlist_append(recipe->parameters, p);

    cx_string_set(path, VIMOS_SEXTRACTOR_CONFIG);
    cx_string_append(path, "/vimos.nnw");

    p = cpl_parameter_new_value("vimos.SExtractor.StarNnwName",
                                CPL_TYPE_STRING,
                                "Name of the neuronal network weights file.",
                                "vimos.SExtractor",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.StarNnwName");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.StarNnwName");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.Timeout",
                                CPL_TYPE_DOUBLE,
                                "Time after which sextractor will be aborted.",
                                "vimos.SExtractor",
                                300.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.Timeout");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.Timeout");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.VerboseType",
                                CPL_TYPE_STRING,
                                "Selects verbosity level.",
                                "vimos.SExtractor",
                                "QUIET");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.VerboseType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.VerboseType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.WeightType",
                                CPL_TYPE_STRING,
                                "Selects the weighting scheme.",
                                "vimos.SExtractor",
                                "NONE");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.WeightType");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.WeightType");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.SExtractor.Window",
                                CPL_TYPE_STRING,
                                "Frame window used in SExtractor detection.",
                                "vimos.SExtractor",
                                "100,100,1948,2340");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.Window");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.Window");
    cpl_parameterlist_append(recipe->parameters, p);

    cx_string_delete(path);


    /*
     * Initialize the VIMOS recipe subsystems (configuration data base,
     * alias tables, messaging facilities) from the current CPL setup.
     */

    status = vmCplRecipeStart(cpl_plugin_get_name(plugin), VERSION);

    if (status) {
        return 1;
    }

    return 0;

}


/*
 * Execute the plugin instance given by the interface.
 */

static cxint
vmimstandard_exec(cpl_plugin *plugin)
{
    cpl_recipe *recipe;
    
    cxint status = 0;

    PilSetOfFrames *sof = NULL;


    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE) 
        recipe = (cpl_recipe *)plugin;
    else 
        return -1;

    if (recipe->parameters == NULL || recipe->frames == NULL) {
        return 1;
    }

    /* Issue a banner */
    vimos_print_banner();

    /*
     * Convert recipe inputs
     */

    sof = newPilSetOfFrames();

    if (sof == NULL) {
        return 1;
    }

    status = vmCplFramesetExport(recipe->frames, sof);

    if (status) {
        deletePilSetOfFrames(sof);
        return 1;
    }

    status = pilRecValidateSet(sof);

    if (!status) {
        deletePilSetOfFrames(sof);
        return 1;
    }

    status = vmCplParlistExport(recipe->parameters);

    if (status) {
        deletePilSetOfFrames(sof);
        return 1;
    }


    /*
     * Execute the data reduction task.
     */

    vmCplRecipeTimerStart(NULL);

    if (vmimstandard(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmimstandard");
        
        if (status == 0) {

            /*
             * Update recipe interface with the product frames found in 'sof'
             * and destroy it.
             */

            status = vmCplFramesetImport(recipe->frames, sof);
        }

    }
    else
        status = 1;

    vmCplRecipeTimerStop(NULL);


    /*
     * Release locally acquired resources
     */

    deletePilSetOfFrames(sof);
    
    return status == 0 ? 0 : 1;

}


static cxint
vmimstandard_destroy(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;


    /*
     * Stop the VIMOS recipe subsystems first.
     */

    vmCplRecipeStop();


    /*
     * We just destroy what was created during the plugin initialization
     * phase, i.e. the parameter list. The frame set is managed by the 
     * application which called us, so we must not touch it.
     */

    if (recipe->parameters != NULL) {
        cpl_parameterlist_delete(recipe->parameters);
    }

    return 0;

}


/** 
 * @brief
 *   Processing of a science or photometric standard star field observation
 *   from a single quadrant.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise @c EXIT_FAILURE is returned.
 *
 * @param sof   Set of frames containing the references to a raw imaging
 *              observation, a master bias and a master flat field. 
 *              Depending on the options used auxiliary frames must be part
 *              of the set of frames (see below).
 *
 * This function provides the the basic reduction steps for imaging
 * observations. This includes at least the removal of the bias,
 * and the flat field correction. Optionally, a cleaning of bad pixels
 * and or cosmic ray hits is preformed on request. The engine can
 * run in two different modes, depending which of the two different
 * image types the engine accepts as input, is actually found in the
 * input set of frames. The two modes are processing a science
 * observation or a photometric standard star field.
 * 
 * Depending on the mode additional or different data reduction steps are
 * executed. In case a science observation is processed, the photometric
 * calibration is applied. This step is skipped for photometric standard
 * star observations.
 * 
 * For science and standard star fields objects are detected and a source
 * list is created, but only for standard star fields the stars are selected
 * from this source list and associated to entries in a reference catalog,
 * which must either be part of the input set of frames or is retrieved
 * online from a catalog server (not yet implemented!). The list of 
 * identifications can be used in the post processor engine
 * @b VmImCalPhot().
 * 
 * Finally, for science and standard star fields the CD and the CCD to
 * sky transformation matrix are used to create the CO matrix which is
 * written to the image header.
 * 
 * On successful termination the engine creates a reduced image and a list
 * of detected sources for science and standard star field observations. For
 * standard star field observations a list containing the identified stars
 * is created in addition.
 * 
 * Control options and additional parameters are read from the recipe
 * configuration database. The recipe function accepts the following
 * task parameters:
 * @li BiasMethod
 * @li CleanBadPixel
 * @li CleanCosmic
 * @li CosmicsRatio
 * @li CosmicsThreshold
 * @li KSigmaClip
 * @li MagInitial
 * @li MagFinal
 * @li MagLimit
 * @li MinStars
 * @li RemoteCatalog
 * @li SearchRadius
 * @li COMatrixSearchRadius
 * @li StarIndex
 * If any of these task parameters is not set in the recipe configuration
 * database the recipe function uses the builtin defaults for these
 * task parameters.
 * 
 * To detect sources in the image the function calls SExtractor. For this
 * purpose SExtractor can be configured in the same way as the engine 
 * itself. For a description of the SExtractor parameters refer to the 
 * SExtractor manual.
 *
 * @author P. Sartoretti, C. Izzo, B. Garilli, R. Palsa
 */   

static cxint
vmimstandard(PilSetOfFrames *sof)
{

  /* FIXME:
   * This code is a copy of the code in vmimobssingle.c, so any change
   * should be applied identically in both modules. This situation will
   * disappear as soon as a full implementation within the CPL scheme
   * will be made.
   */

  const char  task[] = "vmimstandard";

  char  *biasMethodTag = 0;
  char  *rawTag = 0;
  char  *reducedTag = 0;
  char  *productTag = 0;
  char   productName[PATHNAME_MAX + 1];
  char   arcfile[VM_DESC_LENGTH];
  char   filter[VM_DESC_LENGTH];
  char   tplid[VM_DESC_LENGTH];
  float  seeingStart, seeingEnd;
  double beamTemp, ambiTemp;
  double extinction;
  int    quadrant;
  int    nexp, expno;

  unsigned int cleanBadPixel, cleanCosmic, computeQC, remoteCatalog;
  unsigned int reduceAnyFrame;
  unsigned int needsCatalog = 0;
  unsigned int useDark = 1;
  unsigned int tempCheck;

  int biasMethodEntry;
  int minStars;
  int scienceCount = 0;
  int photometryCount = 0;

  float starIndex, magLimit;
  float thresholdCosmics, ratioCosmics;

  double searchRadius;
  double coMatrixSearchRadius;
  double magTolerance1, magTolerance2;
  double tempTolerance;
  double sigmaClip;

  PilFrame *ccdFrame, *biasFrame, *darkFrame;
  PilFrame *flatFrame, *ipcFrame, *astFrame;
  PilFrame *rawFrame;
  PilFrame *productFrame;

  VimosImage *rawImage;
  VimosImage *biasImage, *darkImage, *flatImage;
  VimosImage *resizedBias = 0;
  VimosImage *rawImageFF = 0;
  VimosImage *astFile;

  VimosTable *ccdTable = 0;
  VimosTable *ipcTable = 0;
  VimosTable *astTable = 0;
  VimosTable *galTable = 0;
  VimosTable *stmcTable = 0;
  VimosTable *starTable = 0;

  BiasMethod biasMethod = BIAS_UNDEF;

  enum EngineMode mode = MODE_UNDEF;

  int error;

  int i;

/*
  FILE *fp;
*/


  /*
   * Get task parameters from the recipe database
   */

  /*
   * Retrieve the bias removal method
   */

  biasMethodTag = (char *)pilDfsDbGetString("Parameters", "BiasMethod");

  if ((biasMethodEntry = strselect(biasMethodTag, biasMethodNames,
                                  nBiasMethods)) < 0) {
    cpl_msg_error(task, "%s: Invalid bias removal method!", biasMethodTag);
    return EXIT_FAILURE;
  }
  else 
    biasMethod = biasMethods[biasMethodEntry];


  /*
   * Check if bad pixel correction is required.
   */

  cleanBadPixel = pilDfsDbGetBool("Parameters", "CleanBadPixel", 1);


  /*
   * Check if cosmic ray hit correction is required.
   */

  cleanCosmic = pilDfsDbGetBool("Parameters", "CleanCosmic", 1);
  thresholdCosmics = pilDfsDbGetFloat("Parameters", "CosmicsThreshold", 4.0);
  ratioCosmics = pilDfsDbGetFloat("Parameters", "CosmicsRatio", 2.0);


  /*
   * Retrieve the setup for detection of stars in the input
   * science observation.
   */

  magLimit = pilDfsDbGetFloat("Parameters", "MagLimit", 100.0);

  starIndex = pilDfsDbGetFloat("Parameters", "StarIndex", 0.5);

  if (starIndex < 0. || starIndex > 1.0) {
    cpl_msg_error(task, "Stellarity index is out of range!");
    return EXIT_FAILURE;
  }

  if ((minStars = pilDfsDbGetInt("Parameters", "MinStars", 3)) < 1) {
    cpl_msg_error(task, "Minimum number of required stars is out of range!");
    return EXIT_FAILURE;
  }


  /*
   * Retrieve the setup for CO matrix computation.
   */

  coMatrixSearchRadius = 
          pilDfsDbGetDouble("Parameters", "COMatrixSearchRadius", 6.);

  if (coMatrixSearchRadius < 0.) {
    cpl_msg_error(task, "CO Matrix Search Radius is out of range!");
    return EXIT_FAILURE;
  }


  /*
   * Retrieve the setup for matching the detected stars with an
   * astrometric catalog.
   */

  searchRadius = pilDfsDbGetDouble("Parameters", "SearchRadius", 2.);

  if (searchRadius < 0.) {
    cpl_msg_error(task, "Search Radius is out of range!");
    return EXIT_FAILURE;
  }

  magTolerance1 = pilDfsDbGetDouble("Parameters", "MagInitial", 2.5);
  magTolerance2 = pilDfsDbGetDouble("Parameters", "MagFinal", 1.5);

  sigmaClip = pilDfsDbGetDouble("Parameters", "KSigmaClip", 2.0);

  if (sigmaClip < 0.) {
    cpl_msg_error(task, "Sigma clipping factor is out of range!");
    return EXIT_FAILURE;
  }


  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   * Check if only frames with TPL EXPNO = quadrant number should be reduced.
   * This is effective just in photometric mode, and only if the frame
   * is part of a template sequence.
   */

  reduceAnyFrame = pilDfsDbGetBool("Parameters", "ReduceAnyFrame", 0);


  /*
   * Check whether temperature checks should be used when
   * updating the world coordinate system with a CO matrix
   * and get the temperature tolerance to apply.
   */

  tempCheck = pilDfsDbGetBool("Parameters", "TemperatureCheck", 1);
  tempTolerance = pilDfsDbGetDouble("Parameters", "TemperatureTolerance", 3.);


  /*
   * Check if the astrometric catalog should be downloaded from
   * a remote server, or if the astrometric catalog is expected to
   * be present in the input set of frames.
   */

  /* FIXME:
   * This is not yet supported. Using catlib requires this recipe to
   * become a c++ program! (RP)
   */

  remoteCatalog = pilDfsDbGetBool("Parameters", "RemoteCatalog", 0);


  /*
   * Make sure that all required input frames are present. Exactly
   * one science or standard star field observation from a single
   * quadrant is allowed.
   */

  photometryCount = pilSofFrameCount(sof, pilTrnGetCategory("ImgStandard"));
  scienceCount = pilSofFrameCount(sof, pilTrnGetCategory("ImgScience"));

  if (photometryCount > 0) {
    if (photometryCount > 1) {
      cpl_msg_error(task, "More than one standard star field observation "
                  "found in input!");
      return EXIT_FAILURE;
    }
    else {
      cpl_msg_info(task, "Running in photometry mode ...");
      mode = MODE_PHOTOMETRY;

     /*
      * This is for disabling the magnitude check while building
      * the star match table. It could be restored just if a
      * preliminary photometric calibration is applied also to 
      * the input standard star image before submitting it to 
      * SExtractor. It is questionable whether a preliminary 
      * assumption on the zeropoint should be used in the
      * construction of the star match table, that will be
      * then used to compute the real zeropoint.
      */

      magTolerance1 = -1.0;
      magTolerance2 = -1.0;

    }
  }
  else {
    if (scienceCount > 0) {
      if (scienceCount > 1) {
        cpl_msg_error(task, "More than one science field observation "
                    "found in input!");
        return EXIT_FAILURE;
      }
      else {
        cpl_msg_info(task, "Running in science mode ...");
        mode = MODE_SCIENCE;
      }
    }
    else {
      cpl_msg_error(task, "No science or standard star field observation "
                  "found in input!");
      return EXIT_FAILURE;
    }
  }


  switch (mode) {
      case MODE_PHOTOMETRY:
        rawTag = (char *)pilTrnGetCategory("ImgStandard");
        reducedTag = (char *)pilTrnGetCategory("redImgStandard");
        needsCatalog = 1;

        break;
        
      default:
        rawTag = (char *)pilTrnGetCategory("ImgScience");
        reducedTag = (char *)pilTrnGetCategory("redImgScience");

        break;
  }


  /*
   * Check for the master bias
   */

  if (!(biasFrame = pilSofLookup(sof, pilTrnGetCategory("MasterBias")))) {
    cpl_msg_error(task, "No master bias frame in input!");

    return EXIT_FAILURE;
  }

  pilFrmSetType(biasFrame, PIL_FRAME_TYPE_CALIB);

  /*
   * Check for the master dark. If no master dark frame is found
   * the dark subtraction is skipped.
   */

  if (!(darkFrame = pilSofLookup(sof, pilTrnGetCategory("MasterDark")))) {
    cpl_msg_warning(task, "No master dark frame in input. Skipping dark "
                  "subtraction!");
    useDark = 0;
  }
  else
    pilFrmSetType(darkFrame, PIL_FRAME_TYPE_CALIB);


  /*
   * Check for the master flat field frame
   */

  if (!(flatFrame = pilSofLookup(sof,
                                 pilTrnGetCategory("ImgMasterSkyFlat")))) {
    cpl_msg_error(task, "No master flat field frame in input!");
    return EXIT_FAILURE;
  }

  pilFrmSetType(flatFrame, PIL_FRAME_TYPE_CALIB);

  /*
   * If the bad pixel correction is requested a bad pixel table is
   * required in the input set of frames. If no bad pixel table can
   * be found this is an error. If cosmic ray cleaning is enabled
   * the bad pixel table is used, if it is present in the input set
   * to avoid bad pixels in the surrounding when correcting for cosmics,
   * but a missing bad pixel table will not cause an error.
   */

  ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable"));
  if (ccdFrame)
    pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

  if (cleanBadPixel || cleanCosmic) {
    if (!ccdFrame) {
      if (cleanBadPixel) {
        cpl_msg_error(task, "Bad pixel cleaning requires a CCD table in "
                    "input!");
        return EXIT_FAILURE;
      }
    }
    else {
      cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));
    }
  }


  /*
   * Check for the photometric calibration data if the engine is running
   * in science mode.
   */

  ipcFrame = pilSofLookup(sof, pilTrnGetCategory("PhotometricTable"));
  if (ipcFrame)
    pilFrmSetType(ipcFrame, PIL_FRAME_TYPE_CALIB);

  if (mode == MODE_SCIENCE || (mode == MODE_PHOTOMETRY && computeQC)) {
    if (!ipcFrame) {
      if (mode == MODE_SCIENCE)
        cpl_msg_warning(task, "No photometric calibration data in input: "
                      "no photometric calibration will be applied to the "
                      "reduced data.");
      else
        cpl_msg_error(task, "No photometric calibration data in input!");

      return EXIT_FAILURE;
    }
  }


  /*
   * Load the observation data.
   */

  if ((rawFrame = pilSofLookup(sof, rawTag))) {
    if ((rawImage = openOldFitsFile(pilFrmGetName(rawFrame), 1, 0)) == NULL) {
      cpl_msg_error(task, "Cannot load observation %s!",
                  pilFrmGetName(rawFrame));

      return EXIT_FAILURE;
    }
    else
      closeFitsImage(rawImage, 0);

    pilFrmSetType(rawFrame, PIL_FRAME_TYPE_RAW);
  }

  if ((mode == MODE_PHOTOMETRY) && (!reduceAnyFrame)) {

    /*
     *  In this case the recipe wouldn't reduce a frame belonging to
     *  a quadrant that is not expected to contain any standard star.
     *  This option is intended just for speeding up the online pipeline
     *  processing. The current template for getting data for photometric
     *  calibration consists of a sequence of four exposures, where the
     *  same standard stars field is moved through all quadrants. The
     *  first exposure has the standard field in quadrant 1, the second
     *  exposure in quadrant 2, and so on. So a standard star field 
     *  image will be reduced only if it is part of a template of 
     *  four exposures (i.e., TPL NEXP == 4), and if its sequence 
     *  number in the template is equal to its quadrant number
     *  (i.e., TPL EXPNO == OCS CON QUAD).
     */ 

    if (readIntDescriptor(rawImage->descs, pilTrnGetKeyword("TplExposures"),
                          &nexp, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Integer descriptor %s not found",
                  pilTrnGetKeyword("TplExposures"));

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }

    if (nexp == 4) {
      if (readIntDescriptor(rawImage->descs, 
                            pilTrnGetKeyword("TplExposureNumber"),
                            &expno, NULL) == VM_FALSE) {
        cpl_msg_error(task, "Integer descriptor %s not found",
                    pilTrnGetKeyword("TplExposureNumber"));

        deleteImage(rawImage);

        return EXIT_FAILURE;
      }

      if (readIntDescriptor(rawImage->descs, pilTrnGetKeyword("Quadrant"),
                            &quadrant, NULL) == VM_FALSE) {
        cpl_msg_error(task, "Integer descriptor %s not found",
                    pilTrnGetKeyword("Quadrant"));

        deleteImage(rawImage);
  
        return EXIT_FAILURE;
      }

      if (quadrant != expno) {
        cpl_msg_warning(task, "Frame %s is not reduced, because it is not "
                      "expected to contain any standard stars. This recipe "
                      "run will be aborted!", pilFrmGetName(rawFrame));
        deleteImage(rawImage);
        return EXIT_FAILURE;
      }

    }
  }

  if (computeQC) {

   /*
    * Get the ARCFILE name, and other data characteristics that should
    * be included in QC1 log, from header of input raw frame.
    */
 
    if (readStringDescriptor(rawImage->descs, pilTrnGetKeyword("ArchiveFile"),
                             arcfile, NULL) == VM_FALSE) {
      cpl_msg_error(task, "String descriptor %s not found",
                  pilTrnGetKeyword("ArchiveFile"));

      deleteImage(rawImage);

      return EXIT_FAILURE;

    }


   /*
    * ESO OCS CON QUAD
    */

    if (readIntDescriptor(rawImage->descs, pilTrnGetKeyword("Quadrant"),
                          &quadrant, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Integer descriptor %s not found",
                  pilTrnGetKeyword("Quadrant"));

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }


   /*
    * ESO INS FILT# NAME
    */

    if (readStringDescriptor(rawImage->descs,
                             pilTrnGetKeyword("FilterName", quadrant),
                             filter, NULL) == VM_FALSE) {
      cpl_msg_error(task, "String descriptor %s not found",
                    pilTrnGetKeyword("FilterName", quadrant));

      deleteImage(rawImage);

      return EXIT_FAILURE;

    }


   /*
    * TPL ID
    */

    if (readStringDescriptor(rawImage->descs, pilTrnGetKeyword("TplId"),
                             tplid, NULL) == VM_FALSE) {
      cpl_msg_error(task, "String descriptor %s not found",
                  pilTrnGetKeyword("TplId"));

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }


   /*
    * ESO TEL AMBI FWHM START
    */

    if (readFloatDescriptor(rawImage->descs, pilTrnGetKeyword("SeeingStart"),
                            &seeingStart, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Float descriptor %s not found",
                  pilTrnGetKeyword("SeeingStart"));

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }


   /*
    * ESO TEL AMBI FWHM START
    */

    if (readFloatDescriptor(rawImage->descs, pilTrnGetKeyword("SeeingEnd"),
                            &seeingEnd, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Float descriptor %s not found",
                  pilTrnGetKeyword("SeeingEnd"));

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }


   /*
    * ESO INS FOCU# TEMP
    */

    if (readDoubleDescriptor(rawImage->descs, 
                             pilTrnGetKeyword("BeamTemperature", quadrant),
                             &beamTemp, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Double descriptor %s not found",
                  pilTrnGetKeyword("BeamTemperature", quadrant));

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }


   /*
    * ESO TEL AMBI TEMP
    */

    if (readDoubleDescriptor(rawImage->descs, 
                             pilTrnGetKeyword("AmbientTemperature"),
                             &ambiTemp, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Double descriptor %s not found",
                  pilTrnGetKeyword("AmbientTemperature"));

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }

  }
   

  /*
   * If running in photometry mode and the photometric catalog is not
   * retrieved from a server it must be present in the input set of frames.
   * This is done before the actual processing starts to avoid wasting
   * cpu time in case the connection to a catalog server is not possible
   * due to network failure.
   */

  error = 1;

  astFrame = pilSofLookup(sof, pilTrnGetCategory("PhotometricCatalog"));
  if (astFrame)
    pilFrmSetType(astFrame, PIL_FRAME_TYPE_CALIB);

  if (needsCatalog) {
    if (!remoteCatalog) {
      if (astFrame) {
        if ((astFile = openOldFitsFile(pilFrmGetName(astFrame), 0, 0))) {
          if ((astTable = newAstrometricTable())) {
            if (readFitsAstrometricTable(astTable, astFile->fptr) == VM_TRUE) {
              closeFitsImage(astFile, 0);
              error = 0;
            }
            else {
              cpl_msg_error(task, "Failure reading photometric standard star "
                          "catalog");
              deleteTable(astTable);
            }
          }
          else
            cpl_msg_error(task, "Not enough memory");
        }
        else
          cpl_msg_error(task, "Failure opening astrometric table");
      }
      else
        cpl_msg_error(task, "No input astrometric table found");

      if (error) {
        deleteImage(rawImage);
        return EXIT_FAILURE;
      }

    }
    else {

      /* FIXME:
       *   Call VmImStarCat() or something similar here. This needs
       *   access to the observation's FITS header. So this step
       *   must be done after loading the observation data.
       */

      cpl_msg_error(task, "Remote access to catalog server is not supported!");
      deleteImage(rawImage);

      return EXIT_FAILURE;
    }
  }


  /*
   * Remove the bias from the observation. The master bias does
   * not have overscan areas anymore so they are faked by enlarging
   * the master bias using the observation as reference.
   */

  cpl_msg_info(task, "Restoring overscan regions in master bias ...");

  if (!(biasImage = openOldFitsFile(pilFrmGetName(biasFrame), 1, 0))) {
    cpl_msg_error(task, "Cannot load master bias %s!", pilFrmGetName(biasFrame));

    if (astTable)
      deleteTable(astTable);

    deleteImage(rawImage);

    return EXIT_FAILURE;
  }
  else 
    closeFitsImage(biasImage, 0);

  if (!(resizedBias = growOverscans(biasImage, rawImage))) {
    cpl_msg_error(task, "Restoring overscan regions failed!");

    if (astTable)
      deleteTable(astTable);

    deleteImage(biasImage);
    deleteImage(rawImage);

    return EXIT_FAILURE;
  }

  if (resizedBias != biasImage) {
    deleteImage(biasImage);
    biasImage = resizedBias;
  }
    
  /*
   * Now go for the bias subtraction
   * 
   * NOTE: This includes the overscan subtraction if biasMethod = BIAS_ZMASTER,
   * which is the default. Since the overscan region in the master bias has 
   * been created with growOverscans, there is some extra noise added here,
   * but this can probably be ignored. A better approach is to use the method
   * used in the MOS pipeline: subtract overscan for each individual frame,
   * then subtract bias.
   */

  if (VmSubBias(rawImage, biasImage, biasMethod) == EXIT_FAILURE) {
    cpl_msg_error(task, "Bias correction failed!");
    
    if (astTable)
      deleteTable(astTable);

    deleteImage(biasImage);
    deleteImage(rawImage);

    return EXIT_FAILURE;
  }
  else 
    deleteImage(biasImage);


  /*
   * Dark correction. It is skipped if it was not requested.
   */

  /* FIXME:
   * Dark subtraction should also be skipped if the average dark value
   * is smaller than a given limit. But should this be evaluated at this
   * level? (RP)
   */

  if (useDark) {
    cpl_msg_info(task, "Correcting for dark current ...");

    if (!(darkImage = openOldFitsFile(pilFrmGetName(darkFrame), 1, 0))) {
      cpl_msg_error(task, "Cannot load master dark %s!",
                  pilFrmGetName(darkFrame));

      if (astTable)
        deleteTable(astTable);

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }
    else 
      closeFitsImage(darkImage, 0);

    if (VmSubDark(rawImage, darkImage) == EXIT_FAILURE) {
      cpl_msg_error(task, "Dark subtraction failed!");

      if (astTable)
        deleteTable(astTable);

      deleteImage(darkImage);
      deleteImage(rawImage);

      return EXIT_FAILURE;
    }
    else 
      deleteImage(darkImage);
  }


  /*
   * Flat field correction.
   */

  cpl_msg_info(task, "Performing flat field correction ...");

  if (!(flatImage = openOldFitsFile(pilFrmGetName(flatFrame), 1, 0))) {
    cpl_msg_error(task, "Cannot load master flat field %s!",
                pilFrmGetName(flatFrame));

    if (astTable)
      deleteTable(astTable);

    deleteImage(rawImage);

    return EXIT_FAILURE;
  }
  else 
    closeFitsImage(flatImage, 0);

  if (!(rawImageFF = imageArith(rawImage, flatImage, VM_OPER_DIV))) {
    cpl_msg_error(task, "Flat field correction failed!");

    if (astTable)
      deleteTable(astTable);

    deleteImage(flatImage);
    deleteImage(rawImage);

    return EXIT_FAILURE;
  }
  else {
    copyAllDescriptors(rawImage->descs, &(rawImageFF->descs));

    deleteImage(flatImage);
    deleteImage(rawImage);
  }
  

  /*
   * Load bad pixel data in case bad pixel correction and/or cosmic ray
   * cleaning should be done.
   */

  if (cleanBadPixel || (cleanCosmic && ccdFrame)) {
    if (!(ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
      cpl_msg_error(task, "Cannot load bad pixel data from %s!",
                  pilFrmGetName(ccdFrame));
      return EXIT_FAILURE;
    }
    else {
      closeFitsTable(ccdTable, 0);
    }
  }

 
  /*
   *  Bad pixel cleaning
   */
  
  if (cleanBadPixel) {
    cpl_msg_info(task, "Cleaning bad pixels ...");

    if (cleanBadPixels(rawImageFF, ccdTable, 0) == EXIT_FAILURE) {
      cpl_msg_error(task, "Bad pixel cleaning failed!");

      if (ccdTable)
        deleteTable(ccdTable);

      if (astTable)
        deleteTable(astTable);

      deleteImage(rawImageFF);

      return EXIT_FAILURE;
    }
  }

  
  /*
   * Correct for cosmic ray hits.
   */

  if (cleanCosmic) {

    /* FIXME:
     *  I would prefer the image mode here, assuming there is more 
     *  background than objects, but this requires cleaning up the 
     *  methods for histograms. Let's postpone it a little. (RP)
     */

    float bkgEstimate = imageMedian(rawImageFF);
    float gain = getMeanGainFactor(rawImageFF);
    float ron = computeAverageRon(rawImageFF);

    if (gain < MIN_DIVISOR) {
      cpl_msg_error(task, "Missing or invalid gain factor encountered in %s!",
                  pilFrmGetName(rawFrame));

      if (ccdTable)
        deleteTable(ccdTable);

      if (astTable)
        deleteTable(astTable);

      deleteImage(rawImageFF);

      return EXIT_FAILURE;
    }
    else 
      cpl_msg_info(task, "Removing cosmic ray hits ...");

    if (VmCosmicClean(rawImageFF, ccdTable, 0, bkgEstimate, gain, ron,
                      thresholdCosmics, ratioCosmics) == EXIT_FAILURE) {
      cpl_msg_error(task, "Removal of cosmic ray hits failed!");

      if (ccdTable)
        deleteTable(ccdTable);

      if (astTable)
        deleteTable(astTable);

      deleteImage(rawImageFF);

      return EXIT_FAILURE;
    }
  }

  if (ccdTable)
    deleteTable(ccdTable);


  /* 
   * If running in science mode the photometric calibration is applied.
   * This just updates header keywords. The image data is not altered in
   * any way. In case an observation of a photometric standard star field
   * is processed, this step has to be skipped.
   */
  
  if (mode == MODE_SCIENCE || (mode == MODE_PHOTOMETRY && computeQC)) {

    if (ipcFrame) {
      if (!(ipcTable = openOldFitsTable(pilFrmGetName(ipcFrame), 0))) {
        cpl_msg_error(task, "Cannot load photometric calibration data from %s!",
                    pilFrmGetName(ipcFrame));

        if (astTable)
          deleteTable(astTable);

        return EXIT_FAILURE;
      }
      else 
        closeFitsTable(ipcTable, 0);
    }

  }

  if (mode == MODE_SCIENCE && ipcTable) {

    cpl_msg_info(task, "Performing photometric calibration ...");

    if (!(VmImApplyPhot(rawImageFF, ipcTable))) {
      cpl_msg_error(task, "Photometric calibration failed!");

      if (astTable)
        deleteTable(astTable);

      if (ipcTable)
        deleteTable(ipcTable);

      deleteImage(rawImageFF);

      return EXIT_FAILURE;
    }

    deleteTable(ipcTable);
    ipcTable = 0;

  }


  /*
   * Detect objects in the image. The function uses SExtrator for doing the
   * actual source detection and generates the input SExtractor requires
   * from the image passed to it and the recipe database. The result ASCII
   * source list created by SExtractor is converted into a galaxy table
   * having the sky to CCD transformation matrix and the correction for
   * temperature effects applied.
   */

  cpl_msg_info(task, "Detecting objects ...");

  /* FIXME:
   * Currently the SExtractor output is restricted to what is really
   * making up a galaxy table. Is a configurable output required at
   * this point? (RP)
   */

  if (!(galTable = VmImDetectObjects(rawImageFF, 0))) {
    cpl_msg_error(task, "Object detection failed!");

    if (astTable)
      deleteTable(astTable);

    if (ipcTable)
      deleteTable(ipcTable);

    deleteImage(rawImageFF);

    return EXIT_FAILURE;
  }

  if (computeQC) {

    VimosColumn *stellarity  = findColumn(galTable->cols, "CLASS_STAR");
    VimosColumn *flags  = findColumn(galTable->cols, "FLAGS");
    int   ngals    = galTable->cols->len; 
    int   nstars   = 0;                   /* Stars in galaxy table */
    int   minstars = 10;                  /* Min stars required for QC */


    for (i = 0; i < ngals; i++)
      if (stellarity->colValue->dArray[i] > starIndex 
          && flags->colValue->iArray[i] == 0)
        nstars++;

    if (nstars < minstars) {
      cpl_msg_warning(task, "Less than %d stars in Galaxy Table (%d): "
                    "image quality is not computed.", minstars, nstars);
    }
    else {
      cpl_msg_info(task, "Computing QC1 parameters (image quality)...");

      if (pilQcGroupStart() == EXIT_SUCCESS) {

        VimosColumn *fwhm_world  = findColumn(galTable->cols, "FWHM_WORLD");
        VimosColumn *ellipticity = findColumn(galTable->cols, "ELLIPTICITY");
        VimosColumn *theta       = findColumn(galTable->cols, "THETA_IMAGE");
        VimosColumn *magnitude   = findColumn(galTable->cols, "MAG_BEST");

        char         filterEntry[VM_DESC_LENGTH];
        char         temperEntry[VM_DESC_LENGTH];
        double       average, newAverage, tolerance, delta;
        double       nSigma = 1;
        float        value, minvalue;
        float       *buffer = cpl_malloc(nstars * sizeof(float));
        float       *profile = NULL;
        float       *triprofile = NULL;
        float        meanbin = 1.0; /* Mean no. of objects per theta bin */
        float        mbinsz;
        float        magbeg, magend, maglim, maglimerr;
        int         *ibuffer;
        int          ivalue;
        int          minstar;
        int          binsize, nbin;
        int          nIter, excluded, j;
        int          pos, maxpos;
        int          maxIter = 4;


        pilQcWriteString("PRO.CATG", pilTrnGetCategory("GalaxyTable"), 
                         "Product category");

        if (pilQcWriteString("ARCFILE", arcfile, "Archive File Name")
                                                           == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot write ARCFILE to QC1 logging! (Ignored)");

        if (pilQcWriteString("TPL.ID", tplid, "Template signature ID")
                                                         == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot write TPL.ID to QC1 logging! (Ignored)");

        if (pilQcWriteInt("OCS.CON.QUAD", quadrant, NULL, "Quadrant")
                                                             == EXIT_FAILURE)
          cpl_msg_error(task,
                      "Cannot write OCS.CON.QUAD to QC1 logging! (Ignored)");

        sprintf(filterEntry, "INS.FILT%d.NAME", quadrant);
        if (pilQcWriteString(filterEntry, filter, "Filter name")
                                                         == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot write %s to QC1 logging! (Ignored)",
                      filterEntry);

        if (pilQcWriteDouble("TEL.AMBI.FWHM.START", seeingStart, "arcsec",
                             "ASM seeing at observation start") == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot write TEL.AMBI.FWHM.START "
                      "to QC1 logging! (Ignored)");

        if (pilQcWriteDouble("TEL.AMBI.FWHM.END", seeingEnd, "arcsec",
                             "ASM seeing at observation end") == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot write TEL.AMBI.FWHM.END "
                      "to QC1 logging! (Ignored)");

        sprintf(temperEntry, "INS.FOCU%d.TEMP", quadrant);
        if (pilQcWriteDouble(temperEntry, beamTemp, "C",
                             "Beam temperature") == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot write %s to QC1 logging! (Ignored)", 
                      temperEntry);

        if (pilQcWriteDouble("TEL.AMBI.TEMP", ambiTemp, "C",
                             "ASM ambient temperature") == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot write TEL.AMBI.TEMP "
                      "to QC1 logging! (Ignored)");

        qcCopyValue(rawImageFF->descs, pilTrnGetKeyword("RotatorAngleStart"),
                    NULL, "Adaptor initial rotator angle");

        qcCopyValue(rawImageFF->descs, pilTrnGetKeyword("RotatorAngleEnd"),
                    NULL, "Adaptor final rotator angle");

        /* QC.IMAGE.QUALITY  and  QC.IMAGE.QUALITY.ERROR */

        for (j = 0, i = 0; i < ngals; i++) {
          if (stellarity->colValue->dArray[i] > starIndex
              && flags->colValue->iArray[i] == 0) {
            buffer[j] = (float)fwhm_world->colValue->dArray[i];
            j++;
          }
        }

        /*
         *  Let the first guess for the average be a median: more robust.
         */

        average = medianPixelvalue(buffer, nstars);

        excluded = 1;              /* Just to let it enter the while block */
        nIter = 0;

        while (excluded > 0 && nIter < maxIter) {

          /* 
           * First calculate the nSigma level
           */

          tolerance = 0.0;
          for (i = 0; i < j; i++)
            tolerance += fabs(buffer[i] - average);

          tolerance /= j;
          tolerance *= MEANDEV_TO_SIGMA * nSigma;
    
          cpl_msg_debug(task, "Average FWHM: %f", average);
          cpl_msg_debug(task, "Tolerance   : %f", tolerance);


          /*
           * Then we exclude values departing from average more
           * than NSIGMA sigmas.
           */

          excluded = 0;
          newAverage = 0;
          for (i = 0; i < j; i++) {
            delta = fabs(buffer[i] - average);
            if (delta < tolerance) {
              buffer[i - excluded] = buffer[i];
              newAverage = (((double) i) / (i + 1)) * newAverage
                         + buffer[i] / (i + 1);
            }
            else {
              excluded++;
            }
          }
          average = newAverage;
          j -= excluded;
          cpl_msg_debug(task,"Excluded       : %d\n",excluded);
          nIter++;
        }

        qcWriteValueDouble(galTable->descs, average, 
                           "QC.IMAGE.QUALITY", "arcsec", 
                           "Mean FWHM of stars");

        tolerance /= nSigma;
        /* tolerance /= j; */

        qcWriteValueDouble(galTable->descs, tolerance, 
                           "QC.IMAGE.QUALITY.ERROR", "arcsec", 
                           "Error on mean FWHM of stars");


        /* QC.STAR.COUNT */

        qcWriteValueInt(galTable->descs, nstars, "QC.STAR.COUNT", NULL,
                        "Number of stars used");
        
        /* QC.STELLARITY.MEAN */

        average = computeAverageDouble(stellarity->colValue->dArray, ngals);

        qcWriteValueDouble(galTable->descs, average, "QC.STELLARITY.MEAN", 
                           NULL, "Mean stellarity of all objects");
        
        /* QC.STELLARITY.RMS */

        tolerance = 0.0;
        for (i = 0; i < ngals; i++)
          tolerance += fabs(stellarity->colValue->dArray[i] - average);
 
        tolerance /= ngals;
        tolerance *= MEANDEV_TO_SIGMA;

        qcWriteValueDouble(galTable->descs, tolerance, "QC.STELLARITY.RMS", 
                           NULL, "RMS of stellarity of all objects");

        /* QC.STAR.STELLARITY.MEAN */

        for (j = 0, i = 0; i < ngals; i++) {
          if (stellarity->colValue->dArray[i] > starIndex
              && flags->colValue->iArray[i] == 0) {
            buffer[j] = (float)stellarity->colValue->dArray[i];
            j++;
          }
        }

        average = computeAverageFloat(buffer, nstars);

        qcWriteValueDouble(galTable->descs, average, "QC.STAR.STELLARITY.MEAN", 
                           NULL, "Mean stellarity of all stars");
        
        /* QC.STAR.STELLARITY.RMS */

        tolerance = 0.0;
        for (i = 0; i < nstars; i++)
          tolerance += fabs(buffer[i] - average);
 
        tolerance /= nstars;
        tolerance *= MEANDEV_TO_SIGMA;

        qcWriteValueDouble(galTable->descs, tolerance, 
                           "QC.STAR.STELLARITY.RMS", NULL, 
                           "RMS of stellarity of all stars");

        /* QC.STAR.ELLIPTICITY.MEAN */

        for (j = 0, i = 0; i < ngals; i++) {
          if (stellarity->colValue->dArray[i] > starIndex
              && flags->colValue->iArray[i] == 0) {
            buffer[j] = (float)ellipticity->colValue->dArray[i];
            j++;
          }
        }

        average = computeAverageFloat(buffer, nstars);

        qcWriteValueDouble(galTable->descs, average, 
                           "QC.STAR.ELLIPTICITY.MEAN", NULL, 
                           "Mean ellipticity of all stars");

        /* QC.STAR.ELLIPTICITY.RMS */

        tolerance = 0.0;
        for (i = 0; i < nstars; i++)
          tolerance += fabs(buffer[i] - average);

        tolerance /= nstars;
        tolerance *= MEANDEV_TO_SIGMA;

        qcWriteValueDouble(galTable->descs, tolerance, 
                           "QC.STAR.ELLIPTICITY.RMS", NULL, 
                           "RMS of ellipticity of all stars");

        cpl_free(buffer);

        /* QC.STAR.ORIENTATION.MEAN */
        /* QC.STAR.ORIENTATION.RMS  */

        binsize = 180 * meanbin / nstars;
        if (binsize > 10)
          binsize = 10;
        if (binsize == 7 || binsize == 8)
          binsize = 9;
        if (binsize == 0)
          binsize = 1;

        nbin = 180 / binsize;

        profile = (float *)cpl_calloc(nbin, sizeof(float));
        triprofile = (float *)cpl_calloc(3 * nbin, sizeof(float));

        for (i = 0; i < ngals; i++) {
          if (stellarity->colValue->dArray[i] > starIndex
              && flags->colValue->iArray[i] == 0) {
            value = theta->colValue->dArray[i] + 90.;
            j = floor(value / binsize);
            if (j == nbin)
              j = 0;
            if (j >= 0 && j < nbin)
              profile[j]++;
            else
              cpl_msg_error(task, 
                          "THIS MESSAGE SHOULD NEVER BE PRINTED (%d)!", j);
          }
        }


        /*
         *  Find position of max value in histogram
         */

        maxpos = 0;
        value = profile[maxpos];
        minvalue = profile[maxpos];
        for (i = 0; i < nbin; i++) {
          if (profile[i] > value) {
            value = profile[i];
            maxpos = i;
          }
          if (profile[i] < minvalue)
            minvalue = profile[i];
        }

        /*
         * Build shifted profile
         */

        for (j = 0; j < 3; j++)
          for (i = 0; i < nbin; i++)
            triprofile[i + j * nbin] = profile[i];

        for (i = 0, j = nbin / 2 + maxpos; i < nbin; i++, j++)
          profile[i] = triprofile[j];

        if (findPeak1D(profile, nbin, &value, 3) == VM_FALSE) {
          average = 0.0;
          tolerance = 90.0;
        }
        else {
          average = value + nbin / 2 + maxpos;
          if (average > nbin)
            average -= nbin;

          tolerance = 0.0;
          for (i = 0; i < nbin; i++) {
            tolerance += (i - average)*(i - average)*(profile[i] - minvalue);
            value += (profile[i] - minvalue);
          }

          tolerance /= value;
          tolerance = sqrt(tolerance / value);

          tolerance *= binsize;
          average *= binsize;
          average -= 90.0;

        }

        qcWriteValueDouble(galTable->descs, average,
                           "QC.STAR.ORIENTATION.MEAN", "degree",
                           "Mean orientation of all stars ellipticity");

        qcWriteValueDouble(galTable->descs, tolerance,
                           "QC.STAR.ORIENTATION.RMS", "degree",
                           "RMS of orientation of all stars ellipticity");

        /* QC.MAGLIM */
        /* QC.MAGLIM.ERROR */

        magbeg = 5.0;    /* Histogram range */
        magend = 27.0;
        minstar = 20;    /* Min number of stars required in max bin */

        maglim = 0.0;    /* Value assigned for undefined limiting magnitude */
        maglimerr = 0.0; /* Value assigned for undefined error */

        for (j = 1; j < 10; j++) {
          mbinsz = j / 10.0;
          nbin = (magend - magbeg) / mbinsz;
          ibuffer = cpl_calloc(nbin, sizeof(int));
          for (i = 0; i < ngals; i++) {
            if (stellarity->colValue->dArray[i] > starIndex
                && flags->colValue->iArray[i] == 0) {
              pos = floor((magnitude->colValue->dArray[i] - magbeg) / mbinsz);
              if (pos >= 0 && pos < nbin)
                ibuffer[pos]++;
            }
          }

          maxpos = 0;
          ivalue = ibuffer[maxpos];
          for (i = 0; i < nbin; i++) {
            if (ibuffer[i] > ivalue) {
              ivalue = ibuffer[i];
              maxpos = i;
            }
          }

          cpl_free(ibuffer);

          if (ivalue > minstar) {
            maglim = magbeg + mbinsz * maxpos + mbinsz / 2.0;
            maglimerr = mbinsz;
            break;
          }

        }

        qcWriteValueDouble(galTable->descs, maglim, "QC.MAGLIM", 
                           "mag", "Limiting magnitude");

        qcWriteValueDouble(galTable->descs, maglimerr, "QC.MAGLIM.ERROR", 
                           "mag", "Error on limiting magnitude");

/* %%% */

        if (pilQcGroupEnd() == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

      }
      else
        cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

    }

  }


  /*
   * Select all stars among the detected objects and associate them 
   * to the entries in the astrometric table. The list of identified
   * stars is needed by the to post processors deriving either the
   * image distortions or the photometric calibration.
   */

  if (mode == MODE_PHOTOMETRY) {

    cpl_msg_info(task, "Selecting stars from the list of detected "
               "objects ...");

    if (!(starTable = VmImBuildStarTable(galTable, starIndex, magLimit))) {
      cpl_msg_error(task, "Selection of stars failed!");

      if (astTable)
        deleteTable(astTable);

      if (ipcTable)
        deleteTable(ipcTable);

      deleteTable(galTable);
      deleteImage(rawImageFF);
    
      return EXIT_FAILURE;
    }
    
    if ((starTable->cols) && starTable->cols->len > 0) {

      cpl_msg_info(task, "Associating detected stars with reference "
                   "catalog ...");
    
      if (computeVirtualPixels(rawImageFF->descs, starTable,
                               tempCheck, tempTolerance) != VM_TRUE) {
        cpl_msg_error(task, "Distortion correction of the list of "
                    "stars failed!");
        return EXIT_FAILURE;
      }
    }
    else {
     cpl_msg_warning(task, "Empty Star Table - no star could be selected!");
    }

    stmcTable = VmImBuildStarMatchTable(rawImageFF, starTable, astTable,
                                        minStars, searchRadius,
                                        magTolerance1, magTolerance2,
                                        sigmaClip);
    
    if (!stmcTable) {
      cpl_msg_error(task, "Association failed! Cannot build Star Match "
                  "table.");

      if (astTable)
        deleteTable(astTable);

      if (ipcTable)
        deleteTable(ipcTable);

      deleteTable(galTable);
      deleteTable(starTable);
      deleteImage(rawImageFF);
    
      return EXIT_FAILURE;
    }
 
    if (astTable)
      deleteTable(astTable);

    deleteTable(starTable); 

    if (stmcTable->cols == 0 || stmcTable->cols->len <= 0) {
      mode = MODE_SCIENCE;   /* MODE_PHOTOMETRY no longer applicable */
    }

  }


  if (computeQC) {

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      /*
       *  Currently the central 1600x1800 pixels are used for
       *  sky level estimation.
       */  
 
      int     winSizeX = 1600;
      int     winSizeY = 1800;
      int     winStartX = (rawImageFF->xlen - winSizeX) / 2;
      int     winStartY = (rawImageFF->ylen - winSizeY) / 2;
      int     j, sx, sy, nsec, npix, nvalues;
      double  meanValue, rmsValue, expTime;
      float  *region;
      float  *subregion;
      float  *buffer;
      float  *lvalues;
      char    filterEntry[VM_DESC_LENGTH];


      pilQcWriteString("PRO.CATG", reducedTag, "Product category");

      if (pilQcWriteString("ARCFILE", arcfile, "Archive File Name")
                                                         == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot write ARCFILE to QC1 logging! (Ignored)");

      if (pilQcWriteString("TPL.ID", tplid, "Template signature ID")
                                                       == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot write TPL.ID to QC1 logging! (Ignored)");

      if (pilQcWriteInt("OCS.CON.QUAD", quadrant, NULL, "Quadrant")
                                                           == EXIT_FAILURE)
        cpl_msg_error(task,
                    "Cannot write OCS.CON.QUAD to QC1 logging! (Ignored)");

      sprintf(filterEntry, "INS.FILT%d.NAME", quadrant);
      if (pilQcWriteString(filterEntry, filter, "Filter name")
                                                       == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot write %s to QC1 logging! (Ignored)", 
                    filterEntry);

      if (mode == MODE_PHOTOMETRY) {

        double airmass, magzero, rmszerop;
        double i_zeropoint, i_zero_err, f_extinction = 0, f_ext_err = 0;
        float  gain = getMeanGainFactor(rawImageFF);
        float *fmagValues;
        float  sigma, mean, dummy, dummy1;
        int    numStars, numGoodStars;
        char   magFilColumnName[10];
        VimosColumn *magFilColumn;
        VimosColumn *magBestColumn;

        cpl_msg_info(task, "Computing QC1 parameters (frame zeropoint)...");
  

       /*
        * Add more information from the photometric table and the
        * star match table.
        */

//        qcCopyValue(ipcTable->descs, pilTrnGetKeyword("Extinction"),
//                    "mag", "Extinction coefficient");

//        qcCopyValue(stmcTable->descs, pilTrnGetKeyword("AirMass"),
//                    NULL, "Airmass");

       /* FIXME: currently hardcoded, magzero should be read from 
        * descriptor pilTrnGetKeyword("MagZero") of stmcTable->descs .
        * Anyway magzero is always zero in header...
        */

        magzero = 0.0;

       /*
        * Zeropoint of Star Match Table
        */

        if (readDoubleDescriptor(ipcTable->descs,
                                 pilTrnGetKeyword("Extinction"),
                                 &extinction, NULL) == VM_FALSE) {

          cpl_msg_error(task, "Descriptor %s not found",
                        pilTrnGetKeyword("Extinction"));

          if (stmcTable)
            deleteTable(stmcTable);

          deleteTable(ipcTable);
          deleteTable(galTable);
          deleteImage(rawImageFF);

          return EXIT_FAILURE;

        }

        qcWriteValueDouble(rawImageFF->descs, extinction, "PRO.EXTINCT",
                           "mag", "Extinction coefficient");

        if (readDoubleDescriptor(ipcTable->descs,
                                 pilTrnGetKeyword("MagZero"),
                                 &i_zeropoint, NULL) == VM_FALSE) {

          cpl_msg_error(task, "Descriptor %s not found",
                        pilTrnGetKeyword("MagZero"));

          if (stmcTable)
            deleteTable(stmcTable);

          deleteTable(ipcTable);
          deleteTable(galTable);
          deleteImage(rawImageFF);

          return EXIT_FAILURE;

        }

        if (readDoubleDescriptor(ipcTable->descs,
                                 pilTrnGetKeyword("MagZeroRms"),
                                 &i_zero_err, NULL) == VM_FALSE) {

          cpl_msg_error(task, "Descriptor %s not found",
                        pilTrnGetKeyword("MagZeroRms"));

          if (stmcTable)
            deleteTable(stmcTable);

          deleteTable(ipcTable);
          deleteTable(galTable);
          deleteImage(rawImageFF);

          return EXIT_FAILURE;

        }

        if (readDoubleDescriptor(stmcTable->descs, pilTrnGetKeyword("AirMass"),
                                 &airmass, NULL) == VM_FALSE) {
          cpl_msg_error(task, "Double descriptor %s not found",
                      pilTrnGetKeyword("AirMass"));

          if (stmcTable)
            deleteTable(stmcTable);

          deleteTable(ipcTable);
          deleteTable(galTable);
          deleteImage(rawImageFF);

          return EXIT_FAILURE;

        }

        qcWriteValueDouble(rawImageFF->descs, airmass, "PRO.AIRMASS",
                           NULL, "Airmass");

        sprintf(magFilColumnName, "%s_%s", "MAG", filter);

        magFilColumn = findColumn(stmcTable->cols, magFilColumnName);
        magBestColumn = findColumn(stmcTable->cols, "MAG");

        numStars = stmcTable->cols->len;

        fmagValues = (float *)cpl_malloc(numStars * sizeof(float));

        numGoodStars = 0;

        for (i = 0; i < numStars; i++) {

          if (magFilColumn->colValue->dArray[i] < 50) {
            fmagValues[numGoodStars] = magzero
                          + extinction * airmass
                          + magFilColumn->colValue->dArray[i]
                          - magBestColumn->colValue->dArray[i];
            numGoodStars++;
          }

        }

        if (numGoodStars > 1) {
          xbiwt(fmagValues, numGoodStars, &mean, &sigma, &dummy, &dummy1);
          magzero = (double)mean ;
          rmszerop = (double)sigma ;
        } 
        else if (numGoodStars == 1) {
          magzero = (double)fmagValues[0];
          rmszerop = 0.0;
        }
        else {
          magzero = 0.0;
          rmszerop = 0.0;
        }

        qcWriteValueInt(rawImageFF->descs, numGoodStars, "QC.ZEROPOINT.NSTARS",
                        NULL, "Stars used in zeropoint computation");


       /*
        * Write to reduced image the computed zeropoint - NOT gain
        * corrected!!! This is for astronomers!
        */

        if (writeDoubleDescriptor(&(rawImageFF->descs), 
                                  pilTrnGetKeyword("MagZero"), magzero, 
                                  "Frame zeropoint") == VM_FALSE)
          cpl_msg_error(task, "Cannot write ESO QC ZEROPOINT "
                      "to reduced image header! (Ignored)");

        if (numGoodStars > 0) {

          /*
           * Apply gain correction for QC1 parameters: zeropoint is
           * now in "electrons" - this is for quality controllers!
           */

          magzero += 2.5 * log10(gain);

          f_extinction = extinction + (i_zeropoint - magzero)/airmass;
          f_ext_err = sqrt(i_zero_err*i_zero_err + rmszerop*rmszerop) / airmass;

        }

        qcWriteValueDouble(rawImageFF->descs, magzero, "QC.ZEROPOINT",
                           NULL, "Gain corrected frame zeropoint");

        qcWriteValueDouble(rawImageFF->descs, rmszerop, "QC.ZEROPOINT.ERROR",
                           NULL, "Frame zeropoint error");

        qcWriteValueDouble(rawImageFF->descs, f_extinction, "QC.EXTINCTION",
                           NULL, "Atmospheric extinction");

        qcWriteValueDouble(rawImageFF->descs, f_ext_err, "QC.EXTINCTION.ERROR",
                           NULL, "Error on atmospheric extinction");

        cpl_msg_info(task, " ");

        if (numGoodStars > 0) {
          cpl_msg_info(task, "Frame zeropoint: %5.2f +/- %5.2f (filter %s)", 
                       magzero, rmszerop, filter);
          cpl_msg_info(task, "Atm. extinction: %5.2f +/- %5.2f (filter %s)", 
                       f_extinction, f_ext_err, filter);
          cpl_msg_info(task, "Number of stars used: %d", numGoodStars);
        }
        else
          cpl_msg_info(task, "No good stars found, no zeropoint computation "
                     "(filter %s)", filter);

      }

      /*
       * Sky background estimation. 
       */

      region = extractFloatImage(rawImageFF->data, rawImageFF->xlen, 
                                 rawImageFF->ylen, winStartX, winStartY, 
                                 winSizeX, winSizeY);

      if (!region) {
        cpl_msg_error(task, "Memory allocation!");
        deleteImage(rawImageFF);
        return EXIT_FAILURE;
      }

      nsec = 10;  /* Divide region sides by this number to define subregions */

      npix = (winSizeX / nsec) * (winSizeY / nsec);
      buffer = cpl_calloc(nsec * nsec, sizeof(float));
      j = 0;

      for (sx = 0; sx < winSizeX; sx += winSizeX / nsec) {
        for (sy = 0; sy < winSizeY; sy += winSizeY / nsec) {
          subregion = extractFloatImage(region, winSizeX, winSizeY, sx, sy,
                                        winSizeX / nsec, winSizeY / nsec);
          buffer[j] = medianPixelvalue(subregion, npix);
          j++;
          cpl_free(subregion);
        }
      }

      nvalues = 10;  /* Number of lowest values to consider */
      lvalues = cpl_calloc(nvalues, sizeof(float));

      for (i = 0; i < nvalues; i++)
        lvalues[i] = kthSmallest(buffer, j, i);

      cpl_free(buffer);

      meanValue = computeAverageFloat(lvalues, nvalues);

      rmsValue = 0.0;
      for (i = 0; i < nvalues; i++)
        rmsValue += fabs(lvalues[i] - meanValue); 

      rmsValue /= nvalues;
      rmsValue *= MEANDEV_TO_SIGMA;

      readDoubleDescriptor(rawImageFF->descs,
                           pilTrnGetKeyword("ExposureTime"), &expTime, NULL);

      meanValue /= expTime;
      rmsValue /= expTime;

      qcWriteValueDouble(rawImageFF->descs, meanValue, "QC.SKYBACK",
                         "ADU/s", "Sky level");

      qcWriteValueDouble(rawImageFF->descs, rmsValue, "QC.SKYBACK.ERROR",
                         "ADU/s", "Sky level error");

      cpl_msg_info(task, " ");

      if (pilQcGroupEnd() == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

    }
    else
      cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

    deleteTable(ipcTable);

  }  /* End of QC1 computation */

  /*
   * Write the CO matrix to the header of the reduced image. It is
   * generated from the CD matrix and the sky to CCD transformation
   * found in the image header, i.e. no star match table needs to
   * be passed.
   */

  /* FIXME:
   *   The value 10.0 is passed as sigma clipping factor which is
   *   arbitrarily choosen. BG just said it should be a high value
   *   since the tables are created artificially. But in this case it
   *   would be better to remove it completely from the call and
   *   handle this internally choosing a value appropriate for the
   *   create tables. (RP)
   */

  cpl_msg_info(task, "Updating image world coordinate system ...");

  if (VmAstroComputeCO(rawImageFF, coMatrixSearchRadius, 10., 0L,
                       tempCheck, tempTolerance) != VM_TRUE) {
    cpl_msg_error(task, "CO matrix computation failed!");

    if (stmcTable)
      deleteTable(stmcTable);

    deleteTable(galTable);
    deleteImage(rawImageFF);
    
    return EXIT_FAILURE;
  }


  /*
   * Prepare the products. This includes creation of the temporary product
   * file name from its product category, update of the product header,
   * and the creation of the temporary product file in the current working
   * directory. Finally each product has to be registered with the
   * appropriate attributes in the set of frames.
   */

  /* FIXME:
   * For the moment also keywords which are not task specific
   * are handled here, since this is the last possibility to access
   * the linked list of keywords without reopening the file.
   * This will change in future!
   */


  /*
   * Reduced image
   */

  productTag = reducedTag;

  vmstrlower(strcpy(productName, productTag));
  strcat(productName, ".fits");


  insertDoubleDescriptor(&rawImageFF->descs,
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(rawImageFF),
                         pilTrnGetComment("DataMin"),
                         "ESO*", 1);

  insertDoubleDescriptor(&rawImageFF->descs,
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(rawImageFF),
                         pilTrnGetComment("DataMax"),
                         "ESO*", 1);

  insertDoubleDescriptor(&rawImageFF->descs,
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(rawImageFF), 
                         pilTrnGetComment("DataMedian"),
                         "ESO PRO*", 1);

  insertDoubleDescriptor(&rawImageFF->descs, 
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(rawImageFF),
                         pilTrnGetComment("DataStdDeviation"),
                         "ESO PRO*", 1);

  insertDoubleDescriptor(&rawImageFF->descs,
                         pilTrnGetKeyword("DataMean"),
                         imageMean(rawImageFF),
                         pilTrnGetComment("DataMean"),
                         "ESO PRO*", 1);

  deleteSetOfDescriptors(&rawImageFF->descs, "ESO DPR*");
  //These keywords shouldn't be present in the input raw, but we nevertheless
  //remove them in case they are present (See PIPE-5422)
  deleteSetOfDescriptors(&(rawImageFF->descs), "PCOUNT");  
  deleteSetOfDescriptors(&(rawImageFF->descs), "GCOUNT");  

  if (createFitsImage(productName, rawImageFF, productTag) != VM_TRUE) {
    cpl_msg_error(task, "Cannot create local product file %s!", productName);
 
    if (stmcTable)
      deleteTable(stmcTable);

    deleteTable(galTable);
    deleteImage(rawImageFF);
 
    return EXIT_FAILURE;
  }
  else {
    productFrame = newPilFrame(productName, productTag);

    pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_IMAGE);

    if (mode != MODE_PHOTOMETRY) {
      pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    }
    else {
      pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    }

    pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);
 
    pilSofInsert(sof, productFrame);
  }

  deleteImage(rawImageFF);


  /*
   * Source list
   */

  if (galTable) {
    productTag = (char *)pilTrnGetCategory("GalaxyTable");

    vmstrlower(strcpy(productName, productTag));
    /* strcat(productName, ".TFITS"); */
    strcat(productName, ".fits");
  
    if (writeFitsGalaxyTable(productName, galTable) != VM_TRUE) {
      cpl_msg_error(task, "Cannot create local product file %s!", productName);
 
      if (stmcTable)
        deleteTable(stmcTable);

      deleteTable(galTable);
 
      return EXIT_FAILURE;
    }
    else {
      pilFitsHdrCopy(productName, 0, NULL, ".*-OBS$", 1);
      pilFitsHdrCopy(productName, 0, NULL, "^ESO .*", 1);
      for (i = 0; i < 13; i++)
          pilFitsHdrDeleteKeys(productName, "ESO QC*", 1);

      productFrame = newPilFrame(productName, productTag);

      pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_TABLE);
      pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);
 
      pilSofInsert(sof, productFrame);
    }

    deleteTable(galTable);
  }


  /*
   * List of identified stars. This table is needed for the photometry
   * post processor.
   */

  if (stmcTable) {
    productTag = (char *)pilTrnGetCategory("StarMatchTable");

    vmstrlower(strcpy(productName, productTag));
    /* strcat(productName, ".TFITS"); */
    strcat(productName, ".fits");
  
    if (writeFitsStarMatchTable(productName, stmcTable) != VM_TRUE) {
      cpl_msg_error(task, "Cannot create local product file %s!", productName);
 
      deleteTable(stmcTable);
 
      return EXIT_FAILURE;
    }
    else {
      pilFitsHdrCopy(productName, 0, NULL, ".*-OBS$", 1);
      pilFitsHdrCopy(productName, 0, NULL, "^ESO .*", 1);

      productFrame = newPilFrame(productName, productTag);

      pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_TABLE);

      if (mode == MODE_PHOTOMETRY) {
        pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_PRIMARY);
      }
      else {
        pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      }

      pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);
 
      pilSofInsert(sof, productFrame);
    }

    deleteTable(stmcTable);
  }

  return EXIT_SUCCESS;

}


