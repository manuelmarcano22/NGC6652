/* $Id: vmimobsjitter.c,v 1.14 2013-08-07 16:45:33 cgarcia Exp $
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
 * $Date: 2013-08-07 16:45:33 $
 * $Revision: 1.14 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <string.h>

#include <cxmemory.h>
#include <cxstring.h>

#include <cpl_recipe.h>
#include <cpl_plugininfo.h>
#include <cpl_parameterlist.h>
#include <cpl_frameset.h>

#include <pilmemory.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilrecipe.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilfits.h>

#include "vmimage.h"
#include "vmimagearray.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmgalaxytable.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgastrometry.h"
#include "vmimgresampling.h"
#include "vmimgextraction.h"
#include "vmimgphotcalib.h"
#include "vmimgutils.h"
#include "vmqcutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"

static int vmimobsjitter_create(cpl_plugin *);
static int vmimobsjitter_exec(cpl_plugin *);
static int vmimobsjitter_destroy(cpl_plugin *);

static cxint vmimobsjitter(PilSetOfFrames *);


static char vmimobsjitter_description[] =
"This recipe is used to apply basic reduction steps to a sequence of\n"
"exposures made in direct imaging mode, and to combine them in a single\n"
"image. Each input image is processed in the same way as by recipe\n"
"vmimobsstare, therefore what characterises the vmimobsjitter is just\n"
"the final combination of the input frames.\n\n"
"Input files:\n\n"
"  DO category:              Type:       Explanation:         Required:\n"
"  IMG_SCIENCE               Raw         Science exposure        Y\n"
"  MASTER_BIAS               Calib       Master bias             Y\n"
"  MASTER_DARK               Calib       Master dark             .\n"
"  IMG_MASTER_SKY_FLAT       Calib       Master sky flat         Y\n"
"  PHOT_COEFF_TABLE          Calib       Photometric table       .\n"
"  CCD_TABLE                 Calib       Bad pixel table         .\n\n"
"Output files:\n\n"
"  DO category:              Data type:  Explanation:\n"
"  IMG_SCIENCE_REDUCED       FITS image  Reduced scientific exposure\n"
"  IMG_GALAXY_TABLE          FITS table  List of detected objects\n"
"  IMG_FRINGES               FITS image  Fringes map \n\n"
"If a photometric table is specified, the magnitude zeropoint, the\n"
"atmospheric extinction coefficient, and the colour term are copied\n"
"from the photometric table to the header of the reduced image. The\n" 
"galaxy table is the output of SExtractor converted into FITS format.\n" 
"A CCD table must be specified in input only if a bad pixel cleaning\n"
"is requested.\n\n"
"For more details, please refer to the VIMOS Pipeline User's Guide.";

        
/*
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


/*
 * Definition of the label strings for all methods the recipe function
 * supports for pixel value interpolation.
 */

static const char *resamplingMethodNames[] = {
        "BiLinear",
        "BiCubic"
};

static const ResamplingMethod resamplingMethods[] = {
        RESAMPLING_BILINEAR,
        RESAMPLING_BICUBIC
};

static unsigned int nResamplingMethods = sizeof(resamplingMethods) /
                                         sizeof(ResamplingMethod);


/*
 * Definition of the label strings for all methods the recipe function
 * supports for combining frames and their associated method code.
 */

static const char *stackMethodNames[] = {
        "Auto",
        "Ksigma",
        "MinMax",
        "Median",
        "Average"
};

static const CombMethod stackMethods[] = {
        COMB_AUTO,
        COMB_KSIGMA,
        COMB_REJECT,
        COMB_MEDIAN,
        COMB_AVERAGE
};

static unsigned int nStackMethods = sizeof(stackMethods) / sizeof(CombMethod);



/*
 * Build table of contents, i.e. the list of available plugins, for
 * this module. This function is exported.
 */

cxint
cpl_plugin_get_info(cpl_pluginlist *list)
{

    cpl_recipe *recipe = cx_calloc(1, sizeof *recipe);
    cpl_plugin *plugin = &recipe->interface;


    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmimobsjitter",
                    "Reduce and combine a set of jittered imaging exposures.",
                    vmimobsjitter_description,
                    "ESO VIMOS Pipeline Team and VIMOS Consortium",
                    PACKAGE_BUGREPORT,
                    vimos_get_license(),
                    vmimobsjitter_create,
                    vmimobsjitter_exec,
                    vmimobsjitter_destroy);

    cpl_pluginlist_append(list, plugin);

    return 0;

}


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmimobsjitter_create(cpl_plugin *plugin)
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
     * Create the list of options we accept and hook it into the
     * recipe interface.
     */

    recipe->parameters = cpl_parameterlist_new();
    if (recipe->parameters == NULL) {
        return 1;
    }


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
                                "images before combination.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.quality.enable",
                                CPL_TYPE_BOOL,
                                "Compute QC1 parameters",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ComputeQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ComputeQC");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.interpolation",
                               CPL_TYPE_STRING,
                               "Pixel value interpolation method.",
                               "vimos.Parameters",
                               "BiLinear", 2, "BiLinear", "BiCubic");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "Resampling");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "Resampling");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.fringing",
                                CPL_TYPE_BOOL,
                                "Apply fringing corrections.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "FringingCorr");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "FringingCorr");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.stacking.method",
                               CPL_TYPE_STRING,
                               "Frame combination method",
                               "vimos.Parameters",
                               "Median", 5, "Average", "Median", "MinMax",
                               "Ksigma", "Auto");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "StackMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "StackMethod");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.stacking.ksigma.low",
                                CPL_TYPE_DOUBLE,
                                "Low threshold for K-sigma clipping method.",
                                "vimos.Parameters",
                                5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "KSigmaLow");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "KSigmaLow");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.stacking.ksigma.high",
                                CPL_TYPE_DOUBLE,
                                "High threshold for K-sigma clipping method.",
                                "vimos.Parameters",
                                5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "KSigmaHigh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "KSigmaHigh");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.stacking.minmax.minimum",
                                CPL_TYPE_INT,
                                "Number of lowest rejected values for "
                                "rejection method.",
                                "vimos.Parameters",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MinRejection");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MinRejection");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.stacking.minmax.maximum",
                                CPL_TYPE_INT,
                                "Number of highest rejected values for "
                                "rejection method.",
                                "vimos.Parameters",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MaxRejection");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MaxRejection");
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

    p = cpl_parameter_new_value("vimos.Parameters.matching.ksigma",
                                CPL_TYPE_DOUBLE,
                                "Sigma clipping factor used for CO "
                                "matrix computation.",
                                "vimos.Parameters",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "KSigmaClip");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "KSigmaClip");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

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


    p = cpl_parameter_new_value("vimos.Parameters.matching.stellarity",
                                CPL_TYPE_DOUBLE,
                                "Stellarity index used to select stars "
                                "in the image.",
                                "vimos.Parameters",
                                0.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "StarIndex");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "StarIndex");
    cpl_parameterlist_append(recipe->parameters, p);


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
vmimobsjitter_exec(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;

    cxint status = 0;

    PilSetOfFrames *sof = NULL;


    if (recipe->parameters == NULL || recipe->frames == NULL) {
        return 1;
    }


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

    if (vmimobsjitter(sof) == EXIT_SUCCESS) {

        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmimobsjitter");

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
vmimobsjitter_destroy(cpl_plugin *plugin)
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
 *   Reduce a sequence of imaging exposure taken in jitter mode.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occured,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param sof   Set of frames containing the references to a set of raw
 *              imaging observation taken in shift and stare mode, a master
 *              bias and a master flat field. Depending on the options used
 *              auxiliary frames must be part of the set of frames (see
 *              below).
 *
 * TBD
 *
 * @author P. Sartoretti, R. Palsa
 */   

static cxint
vmimobsjitter(PilSetOfFrames *sof)
{

    const char task[] = "vmimobsjitter";

    char *biasMethodTag = 0L;
    char *stackMethodTag = 0L;
    char *sampleMethodTag = 0L;
    char *productTag = 0L;
    char productName[PATHNAME_MAX + 1];

    size_t  frameCount;
    size_t  minFrames;

    unsigned int cleanBadPixel, computeQC;
    unsigned int fringingCorr;
    unsigned int useDark = 1;
    unsigned int tempCheck;

    int i, nFrames;
    int biasMethodEntry;
    int stackMethodEntry;
    int sampleMethodEntry;

    double coMatrixSearchRadius;
    double sigmaClip;
    double tempTolerance;
    double expTime, sumExpTime;

    float  starIndex;

    PilFrame *ccdFrame, *biasFrame, *darkFrame;
    PilFrame *flatFrame, *ipcFrame;
    PilFrame *rawFrame;
    PilFrame *productFrame;

    VimosImage *biasImage, *darkImage, *flatImage;
    VimosImage *resizedBias = 0L;
    VimosImage *rawImageFF = 0L;
    VimosImage *rawImage;
    VimosImage **set;
    VimosImage *mosaic;
    VimosImage *imaFringes = 0L;

    VimosImageArray *images;

    VimosTable *ccdTable = 0L;
    VimosTable *ipcTable = 0L;
    VimosTable *galTable = 0L;

    BiasMethod biasMethod = BIAS_UNDEF;

    ResamplingMethod sampleMethod = RESAMPLING_UNDEF;

    CombMethod stackMethod = COMB_UNDEF;
    CombParameters combParameter;

  
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
   * Check if the fringing should be corrected.
   */

    fringingCorr = pilDfsDbGetBool("Parameters", "FringingCorr", 0);


    /*
     * Retrieve the interpolation type used for image resampling.
     */

    sampleMethodTag = (char *)pilDfsDbGetString("Parameters", "Resampling");

    if ((sampleMethodEntry = strselect(sampleMethodTag, resamplingMethodNames,
            nResamplingMethods)) < 0) {
        cpl_msg_error(task, "%s: Invalid resampling method!", sampleMethodTag);
        return EXIT_FAILURE;
    }
    else
        sampleMethod = resamplingMethods[sampleMethodEntry];


  /*
   * Retrieve the frame combination method and all method
   * dependent parameters.
   */

    stackMethodTag = (char *)pilDfsDbGetString("Parameters", "StackMethod");

    if ((stackMethodEntry = strselect(stackMethodTag, stackMethodNames,
            nStackMethods)) < 0) {
        cpl_msg_error(task, "%s: Invalid frame combination method.", stackMethodTag);
        return EXIT_FAILURE;
    }
    else 
        stackMethod = stackMethods[stackMethodEntry];

    switch (stackMethod) {
    case COMB_KSIGMA:
        minFrames = MIN_FRAMES_KSIGMA;
        combParameter.kSigmaLow = pilDfsDbGetDouble("Parameters",
                "KSigmaLow", 5.0);
        combParameter.kSigmaHigh = pilDfsDbGetDouble("Parameters",
                "KSigmaHigh", 5.0);
        break;

    case COMB_REJECT:
        minFrames = MIN_FRAMES_REJECT;
        combParameter.minRejection = pilDfsDbGetInt("Parameters",
                "MinRejection", 1);
        combParameter.maxRejection = pilDfsDbGetInt("Parameters",
                "MaxRejection", 1);
        break;

    case COMB_MEDIAN:
        minFrames = MIN_FRAMES_MEDIAN;
        break;

    case COMB_AVERAGE:
        minFrames = MIN_FRAMES_AVERAGE;
        break;

    default:
        cpl_msg_warning(task, "Invalid stacking method. Using default "
                "method 'Average'!");
        stackMethod = COMB_AVERAGE;
        minFrames = MIN_FRAMES_AVERAGE;
        break;
    }


  /*
   * Retrieve the setup for CO matrix computation.
   */

    coMatrixSearchRadius =
            pilDfsDbGetDouble("Parameters", "COMatrixSearchRadius", 6.);

    if (coMatrixSearchRadius < 0.) {
        cpl_msg_error(task, "Search Radius is out of range!");
        return EXIT_FAILURE;
    }

    sigmaClip = pilDfsDbGetDouble("Parameters", "KSigmaClip", 2.0);

    if (sigmaClip < 0.) {
        cpl_msg_error(task, "Sigma clipping factor is out of range!");
        return EXIT_FAILURE;
    }


  /*
   * Check whether temperature checks should be used when
   * updating the world coordinate system with a CO matrix
   * and get the temperature tolerance to apply.
   */

    tempCheck = pilDfsDbGetBool("Parameters", "TemperatureCheck", 1);
    tempTolerance = pilDfsDbGetDouble("Parameters", "TemperatureTolerance", 3.);


    /*
     * Check if enough science observations are available in the input
     * set of frames. Should be more than one.
     */

    frameCount = pilSofFrameCount(sof, pilTrnGetCategory("ImgScience"));

    if (frameCount < minFrames) {
        cpl_msg_error(task, "Too few (%zd) science field observations in input "
                "for combination method '%s'!", frameCount, stackMethodTag);
        return EXIT_FAILURE;
    }


    /*
     * Check if QC1 parameters should be computed
     */

    computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);

    starIndex = pilDfsDbGetFloat("Parameters", "StarIndex", 0.5);


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
     * be found this is an error.
     */

    ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable"));
    if (ccdFrame)
        pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

    if (cleanBadPixel) {
        if (!ccdFrame) {
            cpl_msg_error(task, "Bad pixel cleaning requires a CCD table in input!");
            return EXIT_FAILURE;
        }
        else {
            cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));
        }
    }


    /*
     * Check for the photometric calibration data.
     */

    ipcFrame = pilSofLookup(sof, pilTrnGetCategory("PhotometricCoeffTable"));
    if (!ipcFrame)
      //Backwards-compatible. This was the tag used before Nov 2012. Probably it could be removed in the future
      ipcFrame = pilSofLookup(sof, pilTrnGetCategory("PhotometricTable"));
      if (!ipcFrame)
        cpl_msg_warning(task, "No photometric calibration data in input: "
                  "no photometric calibration will be applied to the "
                  "reduced data.");
    if (ipcFrame)
        pilFrmSetType(ipcFrame, PIL_FRAME_TYPE_CALIB);


    /* FIXME:
     *   The current implementation loads and unloads the calibration images
     *   for each science frame to be processed. An alternative implementation
     *   could load all calibration frames in the beginning and keep them
     *   in memory for the whole preprocessing phase, but this implies a
     *   higher demand on memory. Which solution is best in terms of performance
     *   and memory consumption has to be checked.
     */

  /*
   * Load bad pixel data in case bad pixel correction cleaning
   * should be done.
   */

    if (cleanBadPixel) {
        if (!(ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
            cpl_msg_error(task, "Cannot load bad pixel data from %s!",
                    pilFrmGetName(ccdFrame));
            return EXIT_FAILURE;
        }
        else
            closeFitsTable(ccdTable, 0);
    }


    /*
     * Load the photometric calibration table.
     */

    if (ipcFrame) {
        if (!(ipcTable = openOldFitsTable(pilFrmGetName(ipcFrame), 0))) {
            cpl_msg_error(task, "Cannot load photometric calibration data from %s!",
                    pilFrmGetName(ipcFrame));

            if (ccdTable)
                deleteTable(ccdTable);

            return EXIT_FAILURE;
        }
        else 
            closeFitsTable(ipcTable, 0);
    }


    /*
     * Create a container to collect the preprocessed images.
     */

    if (!(images = newImageArray(frameCount))) {
        cpl_msg_error(task, "Not enough memory");
        return EXIT_FAILURE;
    }


    /*
     * Loop over observation and apply basic preprocessing steps.
     */

    i = 0;
    rawFrame = pilSofLookupNext(sof, pilTrnGetCategory("ImgScience"));

    while (rawFrame) {

        cpl_msg_info(task, "Preprocessing observation %d ...", i);


        /*
         * Load the observation data.
         */

        if ((rawImage = openOldFitsFile(pilFrmGetName(rawFrame), 1, 0)) == 0L) {
            cpl_msg_error(task, "Cannot load observation %s!",
                    pilFrmGetName(rawFrame));

            if (ccdTable)
                deleteTable(ccdTable);

            destroyImageArray(images);

            return EXIT_FAILURE;
        }
        else
            closeFitsImage(rawImage, 0);

        pilFrmSetType(rawFrame, PIL_FRAME_TYPE_RAW);


        /*
         * Remove the bias from the observation. The master bias does
         * not have overscan areas anymore so they are faked by enlarging
         * the master bias using the observation as reference.
         */

        cpl_msg_info(task, "Restoring overscan regions in master bias ...");

        if (!(biasImage = openOldFitsFile(pilFrmGetName(biasFrame), 1, 0))) {
            cpl_msg_error(task, "Cannot load master bias %s!",
                    pilFrmGetName(biasFrame));

            if (ccdTable)
                deleteTable(ccdTable);

            deleteImage(rawImage);
            destroyImageArray(images);

            return EXIT_FAILURE;
        }
        else 
            closeFitsImage(biasImage, 0);


        if (!(resizedBias = growOverscans(biasImage, rawImage))) {
            cpl_msg_error(task, "Restoring overscan regions failed!");

            if (ccdTable)
                deleteTable(ccdTable);

            deleteImage(biasImage);
            deleteImage(rawImage);
            destroyImageArray(images);

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

            if (ccdTable)
                deleteTable(ccdTable);

            deleteImage(biasImage);
            deleteImage(rawImage);
            destroyImageArray(images);

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

                if (ccdTable)
                    deleteTable(ccdTable);

                deleteImage(rawImage);
                destroyImageArray(images);

                return EXIT_FAILURE;
            }
            else 
                closeFitsImage(darkImage, 0);

            if (VmSubDark(rawImage, darkImage) == EXIT_FAILURE) {
                cpl_msg_error(task, "Dark subtraction failed!");

                if (ccdTable)
                    deleteTable(ccdTable);

                deleteImage(darkImage);
                deleteImage(rawImage);
                destroyImageArray(images);

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

            if (ccdTable)
                deleteTable(ccdTable);

            deleteImage(rawImage);
            destroyImageArray(images);

            return EXIT_FAILURE;
        }
        else 
            closeFitsImage(flatImage, 0);

        if (!(rawImageFF = imageArith(rawImage, flatImage, VM_OPER_DIV))) {
            cpl_msg_error(task, "Flat field correction failed!");

            if (ccdTable)
                deleteTable(ccdTable);

            deleteImage(flatImage);
            deleteImage(rawImage);
            destroyImageArray(images);

            return EXIT_FAILURE;
        }
        else {
            copyAllDescriptors(rawImage->descs, &(rawImageFF->descs));

            deleteImage(flatImage);
            deleteImage(rawImage);

            rawImage = rawImageFF;
        }


        /*
         *  Bad pixel cleaning
         */

        if (cleanBadPixel) {
            cpl_msg_info(task, "Cleaning bad pixels ...");

            if (cleanBadPixels(rawImage, ccdTable, 0) == EXIT_FAILURE) {
                cpl_msg_error(task, "Bad pixel cleaning failed!");

                if (ccdTable)
                    deleteTable(ccdTable);

                deleteImage(rawImage);
                destroyImageArray(images);

                return EXIT_FAILURE;
            }
        }


        if (ipcTable) {

            /* 
             * Apply the photometric calibration. This just updates header
             * keywords. The image data is not altered in any way.
             */

            cpl_msg_info(task, "Performing photometric calibration ...");

            if (!(VmImApplyPhot(rawImage, ipcTable))) {
                cpl_msg_error(task, "Photometric calibration failed!");

                if (ccdTable)
                    deleteTable(ccdTable);

                if (ipcTable)
                    deleteTable(ipcTable);

                deleteImage(rawImage);
                destroyImageArray(images);

                return EXIT_FAILURE;
            }
        }


      /*
       * Save preprocessed image in the image list.
       */

        imageArraySet(images, i, rawImage);

        rawFrame = pilSofLookupNext(sof, 0L);
        i++;

    }

    if (ccdTable)
        deleteTable(ccdTable);

    if (ipcTable)
        deleteTable(ipcTable);


    /*
     * Apply fringing correction if requested
     */

    if (fringingCorr) {

        cpl_msg_info(task, "Applying fringing correction...");

        set = imageArrayGetData(images);
        nFrames = imageArraySize(images);

        if (nFrames == 2)
            imaFringes = frCombMinMaxReject(set, 0, 1, 2);
        else
            imaFringes = frCombMedian(set, nFrames, 0);

        copyAllDescriptors(set[0]->descs, &(imaFringes)->descs);

        for (i = 0; i < nFrames; i++)
            imageArithLocal(set[i], imaFringes, VM_OPER_SUB);

    }


    /*
     * Create common coordinate grid and combine the individual images
     * of the jitter sequence.
     */

    if (!VmImResampleImages(images, sampleMethod)) {
        cpl_msg_error(task, "Resampling images to common coordinate grid failed!");

        deleteImage(imaFringes);
        destroyImageArray(images);

        return EXIT_FAILURE;
    }

    set = imageArrayGetData(images);
    nFrames = imageArraySize(images);

    if (!(mosaic = frComb(set, nFrames, stackMethod, &combParameter, 0))) {
        cpl_msg_error(task, "Combination of frames failed!");

        deleteImage(imaFringes);
        destroyImageArray(images);

        return EXIT_FAILURE;
    }

    if (copyAllDescriptors(set[0]->descs, &(mosaic->descs)) == VM_FALSE) {
        cpl_msg_error(task, "Cannot copy keywords to combined image!");

        deleteImage(imaFringes);
        deleteImage(mosaic);
        destroyImageArray(images);

        return EXIT_FAILURE;
    }

    sumExpTime = 0.0;
    for (i = 0; i < nFrames; i++) {
        readDoubleDescriptor(set[i]->descs, pilTrnGetKeyword("ExposureTime"), 
                &expTime, NULL);
        sumExpTime += expTime;
    }


    /*
     * Detect objects in the image. The function uses SExtrator
     * for doing the actual source detection and generates the input
     * SExtractor requires from the image passed to it and the
     * recipe database. The result ASCII source list created by
     * SExtractor is converted into a galaxy table having the sky to
     * CCD transformation matrix and the correction for temperature
     * effects applied.
     */

    cpl_msg_info(task, "Detecting objects ...");

    /* FIXME:
     * Currently the SExtractor output is restricted to what is really
     * making up a galaxy table. Is a configurable output required at
     * this point? (RP)
     */

    if (!(galTable = VmImDetectObjects(mosaic, 0))) {
        cpl_msg_error(task, "Object detection failed!");

        deleteImage(imaFringes);
        deleteImage(mosaic);
        destroyImageArray(images);

        return EXIT_FAILURE;
    }

    if (computeQC) {

        double *stellarity  = tblGetDoubleData(galTable, "CLASS_STAR");
        int    *flags       = tblGetIntData(galTable, "FLAGS");
        int     ngals       = galTable->cols->len;
        int     nstars      = 0;            /* Stars in galaxy table */
        int     minstars    = 10;           /* Min stars required for QC */


        for (i = 0; i < ngals; i++)
            if (stellarity[i] > starIndex && flags[i] == 0)
                nstars++;

        if (nstars < minstars) {
            cpl_msg_warning(task, "Less than %d stars in Galaxy Table (%d): "
                    "image quality is not computed.", minstars, nstars);
        }
        else {

            double *fwhm_world  = tblGetDoubleData(galTable, "FWHM_WORLD");
            double *ellipticity = tblGetDoubleData(galTable, "ELLIPTICITY");
            double *theta       = tblGetDoubleData(galTable, "THETA_IMAGE");
            double *magnitude   = tblGetDoubleData(galTable, "MAG_BEST");
            double *xImage      = tblGetDoubleData(galTable, "X_IMAGE");
            double *yImage      = tblGetDoubleData(galTable, "Y_IMAGE");

            char    descName[VM_DESC_LENGTH];

            double  average, newAverage, tolerance, delta;
            double  nSigma = 1;
            float   value, minvalue;
            float  *buffer = cpl_malloc(ngals * sizeof(float));
            float  *profile = NULL;
            float  *triprofile = NULL;
            float   meanbin = 1.0; /* Mean no. of objects per theta bin */
            float   mbinsz;
            float   magbeg, magend, maglim, maglimerr;
            int    *ibuffer;
            int     ivalue;
            int     minstar;
            int     binsize, nbin;
            int     nIter, excluded, j;
            int     pos, maxpos;
            int     maxIter = 4;
            int     xStart, xEnd, xStep, xpos;
            int     yStart, yEnd, yStep, ypos;
            int     cell;
            int     subngals;

            cpl_msg_info(task, "Computing QC1 parameters (image quality)...");


            /* QC.IMAGE.QUALITY  and  QC.IMAGE.QUALITY.ERROR */

            for (j = 0, i = 0; i < ngals; i++) {
                if (stellarity[i] > starIndex && flags[i] == 0) {
                    buffer[j] = (float)fwhm_world[i];
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

            writeDoubleDescriptor(&(galTable->descs), "ESO QC IMAGE QUALITY",
                                  average, "Mean FWHM of stars [arcsec]");

            tolerance /= nSigma;
            /* tolerance /= j; */

            writeDoubleDescriptor(&galTable->descs, "ESO QC IMAGE QUALITY ERROR",
                                  tolerance, "Error on mean FWHM of stars [arcsec]");

            /* QC.STAR.COUNT */

            writeIntDescriptor(&galTable->descs, "ESO QC STAR COUNT",
                               nstars, "Number of stars used");

            /* QC.STELLARITY.MEAN */

            average = computeAverageDouble(stellarity, ngals);

            writeDoubleDescriptor(&galTable->descs, "ESO QC STELLARITY MEAN",
                                  average, "Mean stellarity of all objects");

            /* QC.STELLARITY.RMS */

            tolerance = 0.0;
            for (i = 0; i < ngals; i++)
                tolerance += fabs(stellarity[i] - average);

            tolerance /= ngals;
            tolerance *= MEANDEV_TO_SIGMA;

            writeDoubleDescriptor(&galTable->descs, "ESO QC STELLARITY RMS",
                                  tolerance, "RMS of stellarity of all objects");

            /* QC.STAR.STELLARITY.MEAN */

            for (j = 0, i = 0; i < ngals; i++) {
                if (stellarity[i] > starIndex && flags[i] == 0) {
                    buffer[j] = (float)stellarity[i];
                    j++;
                }
            }

            average = computeAverageFloat(buffer, nstars);

            writeDoubleDescriptor(&galTable->descs, "ESO QC STAR STELLARITY MEAN",
                                  average, "Mean stellarity of all stars");

            /* QC.STAR.STELLARITY.RMS */

            tolerance = 0.0;
            for (i = 0; i < nstars; i++)
                tolerance += fabs(buffer[i] - average);

            tolerance /= nstars;
            tolerance *= MEANDEV_TO_SIGMA;

            writeDoubleDescriptor(&galTable->descs, "ESO QC STAR STELLARITY RMS",
                                  tolerance, "RMS of stellarity of all stars");

            /* QC.STAR.ELLIPTICITY.MEAN */

            for (j = 0, i = 0; i < ngals; i++) {
                if (stellarity[i] > starIndex && flags[i] == 0) {
                    buffer[j] = (float)ellipticity[i];
                    j++;
                }
            }

            average = computeAverageFloat(buffer, nstars);

            writeDoubleDescriptor(&galTable->descs, "ESO QC STAR ELLIPTICITY MEAN",
                                  average, "Mean ellipticity of all stars");

            /* QC.STAR.ELLIPTICITY.RMS */

            tolerance = 0.0;
            for (i = 0; i < nstars; i++)
                tolerance += fabs(buffer[i] - average);

            tolerance /= nstars;
            tolerance *= MEANDEV_TO_SIGMA;

            writeDoubleDescriptor(&galTable->descs, "ESO QC STAR ELLIPTICITY RMS",
                                  tolerance, "RMS of ellipticity of all stars");


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
                if (stellarity[i] > starIndex && flags[i] == 0) {
                    value = theta[i] + 90.;
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

            writeDoubleDescriptor(&galTable->descs, "ESO QC STAR ORIENTATION MEAN",
                                  average, 
                                  "Mean orientation of all stars ellipticity [degree]");

            writeDoubleDescriptor(&galTable->descs, "ESO QC STAR ORIENTATION RMS",
                                  tolerance,
                                  "RMS of orientation of all stars ellipticity [degree]");


            /*
             * Here the image quality is evaluated on 9 rectangular cells
             * of the image.
             */

            yStart = 80;
            yStep = 760;
            cell = 1;
            for (ypos = 0; ypos < 3; ypos++, yStart += yStep) {
                yEnd = yStart + yStep;
                xStart = 73;
                xStep = 634;
                for (xpos = 0; xpos < 3; xpos++, xStart += xStep, cell++) {
                    xEnd = xStart + xStep;

                    nstars = 0;
                    for (i = 0; i < ngals; i++)
                        if (xImage[i] >= xStart && xImage[i] < xEnd)
                            if (yImage[i] >= yStart && yImage[i] < yEnd)
                                if (stellarity[i] > starIndex && flags[i] == 0)
                                    nstars++;

                    /* QC.STAR.COUNT.CELLn */

                    sprintf(descName, "ESO QC STAR COUNT CELL%d", cell);
                    writeIntDescriptor(&galTable->descs, descName,
                                       nstars, "Number of stars used");


                    if (nstars < minstars) {
                        cpl_msg_warning(task, "Less than %d stars in cell %d (%d): "
                                "no sub-image quality evaluation here.",
                                minstars, cell, nstars);
                        continue;
                    }


                    /* QC.IMAGE.QUALn  and  QC.IMAGE.QUALn.ERROR */

                    for (j = 0, i = 0; i < ngals; i++) {
                        if (xImage[i] >= xStart && xImage[i] < xEnd) {
                            if (yImage[i] >= yStart && yImage[i] < yEnd) {
                                if (stellarity[i] > starIndex && flags[i] == 0) {
                                    buffer[j] = (float)fwhm_world[i];
                                    j++;
                                }
                            }
                        }
                    }


                    /*
                     *  Let the first guess for the average be a median: more robust.
                     */

                    average = medianPixelvalue(buffer, nstars);

                    excluded = 1;       /* Just to let it enter the while block */
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

                        cpl_msg_debug(task, "Average FWHM in cell %d: %f",
                                      cell, average);
                        cpl_msg_debug(task, "Tolerance in cell %d   : %f",
                                      cell, tolerance);


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
                        cpl_msg_debug(task, "Excluded in cell %d    : %d", 
                                      cell, excluded);
                        nIter++;
                    }

                    sprintf(descName, "ESO QC IMAGE QUAL%d", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          average, "Mean FWHM of stars [arcsec]");

                    tolerance /= nSigma;

                    sprintf(descName, "ESO QC IMAGE QUAL%d ERROR", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          tolerance, "Error on mean FWHM of stars [arcsec]");


                    /* QC.STELLARn.MEAN */

                    for (subngals = 0, i = 0; i < ngals; i++) {
                        if (xImage[i] >= xStart && xImage[i] < xEnd) {
                            if (yImage[i] >= yStart && yImage[i] < yEnd) {
                                buffer[subngals] = (float)stellarity[i];
                                subngals++;
                            }
                        }
                    }

                    average = computeAverageFloat(buffer, subngals);

                    sprintf(descName, "ESO QC STELLAR%d MEAN", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          average, "Mean stellarity of all objects");


                    /* QC.STELLARn.RMS */

                    tolerance = 0.0;
                    for (i = 0; i < subngals; i++)
                        tolerance += fabs(buffer[i] - average);

                    tolerance /= subngals;
                    tolerance *= MEANDEV_TO_SIGMA;

                    sprintf(descName, "ESO QC STELLAR%d RMS", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          tolerance, "RMS of stellarity of all objects");

                    /* QC.STAR.STELLARn.MEAN */

                    for (j = 0, i = 0; i < ngals; i++) {
                        if (xImage[i] >= xStart && xImage[i] < xEnd) {
                            if (yImage[i] >= yStart && yImage[i] < yEnd) {
                                if (stellarity[i] > starIndex && flags[i] == 0) {
                                    buffer[j] = (float)stellarity[i];
                                    j++;
                                }
                            }
                        }
                    }

                    average = computeAverageFloat(buffer, nstars);

                    sprintf(descName, "ESO QC STAR STELLAR%d MEAN", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          average, "Mean stellarity of all stars");

                    /* QC.STAR.STELLARn.RMS */

                    tolerance = 0.0;
                    for (i = 0; i < nstars; i++)
                        tolerance += fabs(buffer[i] - average);

                    tolerance /= nstars;
                    tolerance *= MEANDEV_TO_SIGMA;

                    sprintf(descName, "ESO QC STAR STELLAR%d RMS", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          tolerance, "RMS of stellarity of all stars");

                    /* QC.STAR.ELLIPTICn.MEAN */

                    for (j = 0, i = 0; i < ngals; i++) {
                        if (xImage[i] >= xStart && xImage[i] < xEnd) {
                            if (yImage[i] >= yStart && yImage[i] < yEnd) {
                                if (stellarity[i] > starIndex && flags[i] == 0) {
                                    buffer[j] = (float)ellipticity[i];
                                    j++;
                                }
                            }
                        }
                    }

                    average = computeAverageFloat(buffer, nstars);

                    sprintf(descName, "ESO QC STAR ELLIPTIC%d MEAN", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          average, "Mean ellipticity of all stars");

                    /* QC.STAR.ELLIPTIC%d.RMS */

                    tolerance = 0.0;
                    for (i = 0; i < nstars; i++)
                        tolerance += fabs(buffer[i] - average);

                    tolerance /= nstars;
                    tolerance *= MEANDEV_TO_SIGMA;

                    sprintf(descName, "ESO QC STAR ELLIPTIC%d RMS", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          tolerance, "RMS of ellipticity of all stars");



                    /* QC.STAR.ORIENTATn.MEAN */
                    /* QC.STAR.ORIENTATn.RMS  */

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
                        if (xImage[i] >= xStart && xImage[i] < xEnd) {
                            if (yImage[i] >= yStart && yImage[i] < yEnd) {
                                if (stellarity[i] > starIndex && flags[i] == 0) {
                                    value = theta[i] + 90.;
                                    j = floor(value / binsize);
                                    if (j == nbin)
                                        j = 0;
                                    if (j >= 0 && j < nbin)
                                        profile[j]++;
                                    else
                                        cpl_msg_error(task, "YOU SHOULD NOT SEE THIS (%d)!", j);
                                }
                            }
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
                            tolerance += (i-average)*(i-average)*(profile[i]-minvalue);
                            value += (profile[i] - minvalue);
                        }

                        tolerance /= value;
                        tolerance = sqrt(tolerance / value);

                        tolerance *= binsize;
                        average *= binsize;
                        average -= 90.0;

                    }

                    sprintf(descName, "ESO QC STAR ORIENTAT%d MEAN", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          average, "Mean orientation of all stars ellipticity [degree]");

                    sprintf(descName, "ESO QC STAR ORIENTAT%d RMS", cell);
                    writeDoubleDescriptor(&galTable->descs, descName,
                                          tolerance, "RMS of orientation of all stars ellipticity [degree]");

                }
            }


            cpl_free(buffer);


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
                    if (stellarity[i] > starIndex && flags[i] == 0) {
                        pos = floor((magnitude[i] - magbeg) / mbinsz);
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

            writeDoubleDescriptor(&galTable->descs, "ESO QC MAGLIM",
                                  maglim, "Limiting magnitude [mag]");

            writeDoubleDescriptor(&galTable->descs, "ESO QC MAGLIM ERROR",
                                  maglimerr, "Error on limiting magnitude [mag]");

            /* %%% */


        }


        /*
         * Sky background evaluation
         */

        /*
         *  Currently the central 1600x1800 pixels are used for
         *  sky level estimation.
         */

        int     winSizeX = 1600;
        int     winSizeY = 1800;
        int     winStartX;
        int     winStartY;
        int     j, sx, sy, nsec, npix, nvalues;
        double  meanValue, rmsValue;
        float  *region;
        float  *subregion;
        float  *buffer;
        float  *lvalues;


        if (fringingCorr) {
            winStartX = (imaFringes->xlen - winSizeX) / 2;
            winStartY = (imaFringes->ylen - winSizeY) / 2;
        }
        else {
            winStartX = (mosaic->xlen - winSizeX) / 2;
            winStartY = (mosaic->ylen - winSizeY) / 2;
        }

        if (fringingCorr) {
            region = extractFloatImage(imaFringes->data, imaFringes->xlen,
                    imaFringes->ylen, winStartX, winStartY,
                    winSizeX, winSizeY);
        }
        else {
            region = extractFloatImage(mosaic->data, mosaic->xlen,
                    mosaic->ylen, winStartX, winStartY,
                    winSizeX, winSizeY);
        }

        if (!region) {
            cpl_msg_error(task, "Memory allocation!");
            deleteImage(imaFringes);
            deleteImage(mosaic);
            deleteTable(galTable);
            destroyImageArray(images);
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
        cpl_free(region);

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

        meanValue /= expTime;
        rmsValue /= expTime;

        writeDoubleDescriptor(&mosaic->descs, "ESO QC SKYBACK",
                              meanValue, "Sky level [ADU/s]");

        writeDoubleDescriptor(&mosaic->descs, "ESO QC SKYBACK ERROR",
                              rmsValue, "Sky level error [ADU/s]");

        cpl_msg_info(task, " ");


    }

    destroyImageArray(images);

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

    if (VmAstroComputeCO(mosaic, coMatrixSearchRadius, 10., 0,
            tempCheck, tempTolerance) != VM_TRUE) {
        cpl_msg_error(task, "CO matrix computation failed!");

        deleteImage(imaFringes);
        deleteTable(galTable);
        deleteImage(mosaic);

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

    productTag = (char *)pilTrnGetCategory("redImgScience");

    vmstrlower(strcpy(productName, productTag));
    strcat(productName, ".fits");


    insertDoubleDescriptor(&mosaic->descs,
                           pilTrnGetKeyword("ExposureTime"),
                           expTime,
                           pilTrnGetComment("ExposureTime"),
                           "ESO PRO*", 1);

    insertDoubleDescriptor(&mosaic->descs,
                           pilTrnGetKeyword("SummedExposureTime"),
                           sumExpTime,
                           pilTrnGetComment("SummedExposureTime"),
                           "ESO PRO*", 1);

    insertIntDescriptor(&mosaic->descs,
                        pilTrnGetKeyword("NFramesCombined"),
                        nFrames,
                        pilTrnGetComment("NFramesCombined"),
                        "ESO PRO*", 1);

    insertDoubleDescriptor(&mosaic->descs,
                           pilTrnGetKeyword("DataMin"),
                           imageMinimum(mosaic),
                           pilTrnGetComment("DataMin"),
                           "ESO*", 1);

    insertDoubleDescriptor(&mosaic->descs,
                           pilTrnGetKeyword("DataMax"),
                           imageMaximum(mosaic),
                           pilTrnGetComment("DataMax"),
                           "ESO*", 1);

    insertDoubleDescriptor(&mosaic->descs,
                           pilTrnGetKeyword("DataMedian"),
                           imageMedian(mosaic),
                           pilTrnGetComment("DataMedian"),
                           "ESO PRO*", 1);

    insertDoubleDescriptor(&mosaic->descs,
                           pilTrnGetKeyword("DataStdDeviation"),
                           imageSigma(mosaic), 
                           pilTrnGetComment("DataStdDeviation"),
                           "ESO PRO*", 1);

    insertDoubleDescriptor(&mosaic->descs,
                           pilTrnGetKeyword("DataMean"),
                           imageMean(mosaic),
                           pilTrnGetComment("DataMean"),
                           "ESO PRO*", 1);

    /*
     * Create the product file on disk, set the product attributes and
     * update the set of frames.
     */

    if (createFitsImage(productName, mosaic, productTag) != VM_TRUE) {
        cpl_msg_error(task, "Cannot create local product file %s!", productName);

        deleteImage(imaFringes);
        deleteTable(galTable);
        deleteImage(mosaic);

        return EXIT_FAILURE;
    }
    else {
        productFrame = newPilFrame(productName, productTag);

        pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);
        pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_IMAGE);
        pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_PRIMARY);
        pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);

        pilSofInsert(sof, productFrame);
    }

    deleteImage(mosaic);


    /*
     * Source list
     */

    productTag = (char *)pilTrnGetCategory("GalaxyTable");

    vmstrlower(strcpy(productName, productTag));
    /* strcat(productName, ".TFITS"); */
    strcat(productName, ".fits");

    if (writeFitsGalaxyTable(productName, galTable) != VM_TRUE) {
        cpl_msg_error(task, "Cannot create local product file %s!", productName);

        deleteImage(imaFringes);
        deleteTable(galTable);

        return EXIT_FAILURE;
    }
    else {
        pilFitsHdrCopy(productName, 0, NULL, ".*-OBS$", 1);
        pilFitsHdrCopy(productName, 0, NULL, "^ESO .*", 1);

        productFrame = newPilFrame(productName, productTag);

        pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);
        pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_TABLE);
        pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_SECONDARY);
        pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);

        pilSofInsert(sof, productFrame);
    }

    deleteTable(galTable);

    /* FIXME:
     * For the moment also keywords which are not task specific
     * are handled here, since this is the last possibility to access
     * the linked list of keywords without reopening the file.
     * This will change in future!
     */

    if (fringingCorr) {

        productTag = (char *)pilTrnGetCategory("ImgFringes");

        vmstrlower(strcpy(productName, productTag));
        strcat(productName, ".fits");


        insertDoubleDescriptor(&imaFringes->descs,
                               pilTrnGetKeyword("ExposureTime"),
                               expTime,
                               pilTrnGetComment("ExposureTime"),
                               "ESO PRO*", 1);

        insertDoubleDescriptor(&imaFringes->descs,
                               pilTrnGetKeyword("SummedExposureTime"),
                               sumExpTime,
                               pilTrnGetComment("SummedExposureTime"),
                               "ESO PRO*", 1);

        insertIntDescriptor(&imaFringes->descs,
                            pilTrnGetKeyword("NFramesCombined"),
                            nFrames,
                            pilTrnGetComment("NFramesCombined"),
                            "ESO PRO*", 1);

        insertDoubleDescriptor(&imaFringes->descs,
                               pilTrnGetKeyword("DataMin"),
                               imageMinimum(imaFringes),
                               pilTrnGetComment("DataMin"),
                               "ESO*", 1);

        insertDoubleDescriptor(&imaFringes->descs,
                               pilTrnGetKeyword("DataMax"),
                               imageMaximum(imaFringes),
                               pilTrnGetComment("DataMax"),
                               "ESO*", 1);

        insertDoubleDescriptor(&imaFringes->descs,
                               pilTrnGetKeyword("DataMedian"),
                               imageMedian(imaFringes),
                               pilTrnGetComment("DataMedian"),
                               "ESO PRO*", 1);

        insertDoubleDescriptor(&imaFringes->descs,
                               pilTrnGetKeyword("DataStdDeviation"),
                               imageSigma(imaFringes), 
                               pilTrnGetComment("DataStdDeviation"),
                               "ESO PRO*", 1);

        insertDoubleDescriptor(&imaFringes->descs,
                               pilTrnGetKeyword("DataMean"),
                               imageMean(imaFringes),
                               pilTrnGetComment("DataMean"),
                               "ESO PRO*", 1);

        /*
         * Create the product file on disk, set the product attributes and
         * update the set of frames.
         */

        if (createFitsImage(productName, imaFringes, productTag) != VM_TRUE) {
            cpl_msg_error(task, "Cannot create local product file %s!", productName);

            deleteImage(imaFringes);

            return EXIT_FAILURE;
        }
        else {
            productFrame = newPilFrame(productName, productTag);

            pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);
            pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_IMAGE);
            pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_PRIMARY);
            pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);

            pilSofInsert(sof, productFrame);
        }

        deleteImage(imaFringes);

    }

    return EXIT_SUCCESS;

}


