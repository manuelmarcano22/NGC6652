/* $Id: vmskyccd.c,v 1.10 2012-11-08 17:56:57 cgarcia Exp $
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
 * $Date: 2012-11-08 17:56:57 $
 * $Revision: 1.10 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <string.h>

#include <cxmemory.h>
#include <cxmessages.h>
#include <cxstring.h>

#include <cpl_recipe.h>
#include <cpl_plugininfo.h>
#include <cpl_parameterlist.h>
#include <cpl_frameset.h>

#include <pilmemory.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilrecipe.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmgalaxytable.h"
#include "vmstartable.h"
#include "vmstarmatchtable.h"
#include "vmastrometrictable.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgextraction.h"
#include "vmimgastrometry.h"
#include "vmimgphotcalib.h"
#include "vmimgutils.h"
#include "vmcpl.h"


static cxint vmskyccd(PilSetOfFrames *);


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
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it availble to the application using the interface.
 */

static cxint
vmskyccd_create(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;

    cpl_parameter *p;

    cx_string *path = cx_string_new();

    cxint status = 0;


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

    p = cpl_parameter_new_value("vimos.Parameters.zeropoint",
                                CPL_TYPE_BOOL,
                                "Apply zeropoint from input photometric "
                                "table before running SExtractor.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyZeropoint");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyZeropoint");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.bias.removing.method",
                                CPL_TYPE_STRING,
                                "Method used for removing the bias.",
                                "vimos.Parameters",
                                "Zmaster", 2, "Zmaster", "Master");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "BiasMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "BiasMethod");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.badpixel.clean",
                                CPL_TYPE_BOOL,
                                "Bad pixel correction on input science "
                                "image.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.clean",
                                CPL_TYPE_BOOL,
                                "Cosmic ray removal from the input science "
                                "image.",
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
                                "Stellarity index used to select stars in "
                                "the image.",
                                "vimos.Parameters",
                                0.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "StarIndex");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "StarIndex");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.radius",
                                CPL_TYPE_DOUBLE,
                                "Aperture used for object identification.",
                                "vimos.Parameters",
                                2.0);
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
                                "Check beam temperature when updating the "
                                "world coordinate system.",
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
vmskyccd_exec(cpl_plugin *plugin)
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

    if (vmskyccd(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmskyccd");
        
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
vmskyccd_destroy(cpl_plugin *plugin)
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
 * @memo 
 *   Upgrade Filter Table with CCD to Sky transformation.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE.
 *
 * @param sof names of input files:  Astrometric Image,
 *            corresponding Astrometric Table, Master Bias, Master Flat Field
 *            optional master Dark, optionl CCD Table. 
 *
 * @doc 
 *   This recipe first performs the basic operation steps
 *   on an astrometric field (see VmImObsSingle for details)
 *
 *   Objects are detected by Sextractor and a source list
 *   is created, the stars are selected
 *   from this source list and associated to entries in a reference catalog,
 *   which must either be part of the input set of frames or is retrieved
 *   online from a catalog server (not yet implemented!). 
 *
 *   The list of identifications is used to derive the Sky To CCD
 *   matrix coefficients, i.e. the coefficients of the matrix which
 *   corrects for optical distortions (the name is misleading,
 *   and has historical reasons. A more correct naming wouldbe
 *   RealToVirtual matrix). The algortihm is is the same as in
 *   \textbf{VmImCalOpt()}.
 *   Also the inverse matrix is computed, even if it is presently never used,
 *   and the results are written in the appropriate PAF file
 *
 *   On successful termination the engine creates a reduced image and 
 *   the PAF file, to be sent to Instrument Workstation
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
 *   \begin{itemize}
 *     \item BiasMethod
 *     \item CleanBadPixel
 *     \item CleanCosmic
 *     \item CosmicsRatio
 *     \item CosmicsThreshold
 *     \item KSigmaClip
 *     \item MagInitial
 *     \item MagFinal
 *     \item MagLimit
 *     \item RemoteCatalog
 *     \item SearchRadius
 *     \item StarIndex
 *   \end{itemize}
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the builtin defaults for these
 *   task parameters.
 *
 *   To detect sources in the image the function calls SExtractor. For this
 *   purpose SExtractor can be configured in the same way as the engine 
 *   itself, i.e. through the recipe configuration file. The only restriction
 *   is that, due to the sheer number of parameters the SExtractor parameters
 *   are not offered on the command line, but have to be placed in the
 *   configuration group \textbf{SExtractor} of the recipe configuration
 *   file.
 *
 * @author  B.Garilli
 */

static cxint 
vmskyccd(PilSetOfFrames *sof)
{

  const char  task[] = "vmskyccd";

  const char  parameter[]        = "Parameters";
  char       *rawTag = 0;
  char       *reducedTag = 0;
  char       *productTag = 0;
  char       *starMatchProductTag = 0;
  char        productName[PATHNAME_MAX + 1];
  char        starMatchProductName[PATHNAME_MAX + 1];
  char        namePAF[]          = "IMG_sky2ccd";
  char        pipeNamePAF[]      = "IMG_sky2ccd_";
  char       *pafFileName;
  char        descStringVal[80];


  int         numStars, allStars;

  PilFrame *ccdFrame, *biasFrame, *darkFrame;
  PilFrame *flatFrame, *astFrame;
  PilFrame *rawFrame;
  PilFrame *productFrame;
  PilFrame *starMatchProductFrame;
  PilFrame *ipcFrame = NULL;

  VimosImage *rawImage;
  VimosImage *biasImage, *darkImage, *flatImage;
  VimosImage *resizedBias = 0;
  VimosImage *rawImageBS = 0, *rawImageDK = 0, *rawImageFF = 0;
  VimosImage *astFile;

  VimosTable  *starMatchTable;
  VimosTable *ccdTable = 0;
  VimosTable *astTable = 0;
  VimosTable *galTable = 0;
  VimosTable *ipcTable = 0;
  VimosTable *starTable;

  VimosPixel *surfaceX  = NULL;
  VimosPixel *surfaceY  = NULL;
  VimosPixel *newPos    = NULL;

  char *biasMethodTag = 0;
  BiasMethod biasMethod = BIAS_UNDEF;
  int biasMethodEntry;

  int frameCount = 0;

  unsigned int applyZeropoint, cleanBadPixel, cleanCosmic, remoteCatalog;
  unsigned int useDark = 1;
  unsigned int tempCheck;

  float thresholdCosmics, ratioCosmics;

  float starIndex, magLimit;

  double      searchRadius;
  double      magTolerance1, magTolerance2;
  double      tempTolerance;
  double      sigmaClip;

  double      xrms = 0.;
  double      yrms = 0.;
  double      deltaXrms, deltaYrms, oldXrms, oldYrms;
  double     *coeffX;
  double     *coeffY;
  float       diffX, diffY;
  char       *coeffControlString;
  int        *exclude;
  int         fitOrderX, fitOrderY;
  int         noFit, noReject;
  int         noFound;
  int         noIter;
  int         noIterations;
  int         minStarNo;
  int         nextCycle = 1; /* flag for the rejection cycle */
  int         i, j, k, n, error;
  int         dummy;     /* Dummy int argument for fitSurfacePolynomial() */

  struct WorldCoor    * wcs = 0;     /* WCS structure */
  double     *ximage, *yimage, *ximage_orig, *yimage_orig;
  double     *xVirtual, *yVirtual;
  double     *ximage_tmp, *yimage_tmp;
  int        *star, *star_orig;
  int        *star_tmp;

  int debug = 1;

  /*
   * Get task parameters from the recipe database
   */

  applyZeropoint = pilDfsDbGetBool("Parameters", "ApplyZeropoint", 1);

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


  /*
   * Retrieve the setup for matching the detected stars with an
   * astrometric catalog.
   */

  searchRadius = pilDfsDbGetDouble("Parameters", "SearchRadius", 2.);

  if (searchRadius < 0.) {
    cpl_msg_error(task, "Search Radius is out of range!");
    return EXIT_FAILURE;
  }

  if (applyZeropoint) {
    magTolerance1 = pilDfsDbGetDouble("Parameters", "MagInitial", 2.5);
    magTolerance2 = pilDfsDbGetDouble("Parameters", "MagFinal", 1.5);
  }
  else {

   /*
    * In this way the comparison of SExtractor magnitudes with catalog 
    * magnitudes is avoided altogether.
    */

    magTolerance1 = -1.;
    magTolerance2 = -1.;
  }

  sigmaClip = pilDfsDbGetDouble("Parameters", "KSigmaClip", 2.0);

  if (sigmaClip < 0.) {
    cpl_msg_error(task, "Sigma clipping factor is out of range!");
    return EXIT_FAILURE;
  }

 /*
  * Get degrees of the bivariate polynomial used to fit the mask-CCD
  * transformation. They should be set to 3, 3 to match the MPS request.
  */ 

  fitOrderX = pilDfsDbGetInt(parameter, "PolyOrderX", 3);
  fitOrderY = pilDfsDbGetInt(parameter, "PolyOrderY", 3);

  if (fitOrderX < 2 || fitOrderY < 2) {
    cpl_msg_error(task, "Invalid choice of polynomial order for mask to "
                "CCD transformation: order should be at least 2");
    return EXIT_FAILURE;
  }

 /*
  * Get max number of model fit iterations.
  */

  noIterations = 10;//pilDfsDbGetInt(parameter, "Iterations", 1);


  /*
   * Check whether temperature checks should be used when
   * updating the world coordinate system with a CO matrix
   * and get the temperature tolerance to apply.
   */

  tempCheck = pilDfsDbGetBool(parameter, "TemperatureCheck", 1);
  tempTolerance = pilDfsDbGetDouble(parameter, "TemperatureTolerance", 3.);


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
   * quadrant is processed.
   */

  frameCount = pilSofFrameCount(sof, pilTrnGetCategory("ImgAstrometry"));

  if ( frameCount != 1) {
    cpl_msg_error(task, "No or more than one science or standard star "
                "field observation found in input!");
    return EXIT_FAILURE;
  }
  else {
    rawTag = (char *)pilTrnGetCategory("ImgAstrometry");
    reducedTag = (char *)pilTrnGetCategory("redImgAstrometry");
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
   * If the photometric calibration was requested, the photometric table
   * must be in the input set of frames.
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


  if (applyZeropoint) {

    error = 1;

    if (ipcFrame) {
      if ((ipcTable = openOldFitsTable(pilFrmGetName(ipcFrame), 0))) {
        closeFitsTable(ipcTable, 0);
        error = 0;
      }
      else
        cpl_msg_error(task, "Cannot load photometric calibration data from %s!",
                    pilFrmGetName(ipcFrame));
    }
    else
      cpl_msg_error(task, "No photometric calibration data in input!");

    if (error)
      return EXIT_FAILURE;

  }


  /*
   * If the astrometric catalog is not retrieved from a server it must 
   * be present in the input set of frames. This is done before the 
   * actual processing starts to avoid wasting cpu time in case the 
   * connection to a catalog server is not possible due to network 
   * failure.
   */

  astFrame = pilSofLookup(sof, pilTrnGetCategory("AstrometricTable"));
  if (astFrame)
    pilFrmSetType(astFrame, PIL_FRAME_TYPE_CALIB);

  if (!remoteCatalog) {

    if (astFrame) {
      if ((astFile = openOldFitsFile(pilFrmGetName(astFrame), 0, 0))) {
        if ((astTable = newAstrometricTable())) {
          if (readFitsAstrometricTable(astTable, astFile->fptr) == VM_TRUE) {
            closeFitsImage(astFile, 0);
            error = 0;
          }
          else {
            cpl_msg_error(task, "Cannot read astrometric catalog from %s",
                        pilFrmGetName(astFrame));
            deleteTable(astTable);
          }
        }
        else
          cpl_msg_error(task, "Not enough memory");
      }
      else
        cpl_msg_error(task, 
                    "Cannot load astrometric reference catalog data from %s!",
                    pilFrmGetName(astFrame));
    }
    else
      cpl_msg_error(task, "No input astrometric catalog found");

    if (error) {
      if (ipcTable)
        deleteTable(ipcTable);
      return EXIT_FAILURE;
    }

  } 
  else {

    /* FIXME:
     *   Call VmImStarCat() or something similar here
     */
      
    cpl_msg_error(task, "Remote access to catalog server is not supported!");

    if (ipcTable)
      deleteTable(ipcTable);

    return EXIT_FAILURE;
  }

  /*
   * Load the observation data.
   */

  rawFrame = pilSofLookup(sof, rawTag);
  if (rawFrame)
    pilFrmSetType(rawFrame, PIL_FRAME_TYPE_RAW);

  if (debug == 1) {
    if ((rawFrame = pilSofLookup(sof, rawTag))) {
      if ((rawImage = 
           openOldFitsFile(pilFrmGetName(rawFrame), 1, 0)) == NULL) {
        cpl_msg_error(task, "Cannot load observation %s!",
                    pilFrmGetName(rawFrame));

        if (ipcTable)
          deleteTable(ipcTable);

        if (astTable)
          deleteTable(astTable);

        return EXIT_FAILURE;
      }
    }
   

   /*
    * Remove the bias from the observation. The master bias does
    * not have overscan areas anymore so they are faked by enlarging
    * the master bias using the observation as reference.
    */

    cpl_msg_info(task, "Restoring overscan regions in master bias...");

    if (!(biasImage = openOldFitsFile(pilFrmGetName(biasFrame), 1, 0))) {
      cpl_msg_error(task, "Cannot load master bias %s!", 
                  pilFrmGetName(biasFrame));

      if (ipcTable)
        deleteTable(ipcTable);

      if (astTable)
        deleteTable(astTable);

      deleteImage(rawImage);

      return EXIT_FAILURE;
    }
    else 
      closeFitsImage(biasImage, 0);

    if (!(resizedBias = growOverscans(biasImage, rawImage))) {
      cpl_msg_error(task, "Restoring overscan regions failed!");

      if (ipcTable)
        deleteTable(ipcTable);

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

      if (ipcTable)
        deleteTable(ipcTable);
    
      if (astTable)
        deleteTable(astTable);

      deleteImage(biasImage);
      deleteImage(rawImage);

      return EXIT_FAILURE;
    }
    else {
      deleteImage(biasImage);
      rawImageBS = rawImage;
    }


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

        if (ipcTable)
          deleteTable(ipcTable);
  
        if (astTable)
          deleteTable(astTable);

        deleteImage(rawImageBS);
  
        return EXIT_FAILURE;
      }
      else 
        closeFitsImage(darkImage, 0);

      if (VmSubDark(rawImageBS, darkImage) == EXIT_FAILURE) {
        cpl_msg_error(task, "Dark subtraction failed!");

        if (ipcTable)
          deleteTable(ipcTable);

        if (astTable)
          deleteTable(astTable);

        deleteImage(darkImage);
        deleteImage(rawImageBS);

        return EXIT_FAILURE;
      }
      else {
        deleteImage(darkImage);

        /* FIXME: Needed because VmSubDark does not return a new image. Why? */
        rawImageDK = rawImageBS;
      }
    }
    else {
      rawImageDK = rawImageBS;
    }

   /*
    * Flat field correction.
    */

    cpl_msg_info(task, "Performing flat field correction ...");

    if (!(flatImage = openOldFitsFile(pilFrmGetName(flatFrame), 1, 0))) {
      cpl_msg_error(task, "Cannot load master flat field %s!",
                  pilFrmGetName(flatFrame));

      if (ipcTable)
        deleteTable(ipcTable);

      if (astTable)
        deleteTable(astTable);

      deleteImage(rawImageDK);

      return EXIT_FAILURE;
    }
    else
      closeFitsImage(flatImage, 0);

    if (!(rawImageFF = imageArith(rawImageDK, flatImage, VM_OPER_DIV))) {
      cpl_msg_error(task, "Flat field correction failed!");

      if (ipcTable)
        deleteTable(ipcTable);

      if (astTable)
        deleteTable(astTable);

      deleteImage(flatImage);
      deleteImage(rawImageDK);

      return EXIT_FAILURE;
    }
    else {
      copyAllDescriptors(rawImageDK->descs, &(rawImageFF->descs));

      deleteImage(flatImage);
      deleteImage(rawImageDK);
    }
  

   /*
    * Load bad pixel data in case bad pixel correction and/or cosmic ray
    * cleaning should be done.
    */

    if (cleanBadPixel || (cleanCosmic && ccdFrame)) {
      if (!(ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
        cpl_msg_error(task, "Cannot load bad pixel data from %s!",
                    pilFrmGetName(ccdFrame));

        if (ipcTable)
          deleteTable(ipcTable);
  
        if (astTable)
          deleteTable(astTable);

        deleteImage(flatImage);
        deleteImage(rawImageDK);

        return EXIT_FAILURE;
      }
      else
        closeFitsTable(ccdTable, 0);
    }

 
    /*
     *  Bad pixel cleaning
     */
  
    if (cleanBadPixel) {
      cpl_msg_info(task, "Cleaning bad pixels ...");

      if (cleanBadPixels(rawImageFF, ccdTable, 0) == EXIT_FAILURE) {
        cpl_msg_error(task, "Bad pixel cleaning failed!");

        if (ipcTable)
          deleteTable(ipcTable);

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

        if (ipcTable)
          deleteTable(ipcTable);

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

        if (ipcTable)
          deleteTable(ipcTable);

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

    createFitsImage("SkyCcdBias.fits", rawImageFF, productTag) ;
  } 
  else {
    rawImageFF = openOldFitsFile("SkyCcdBias.fits", 1, 0);
  }

 /*
  * Before running SExtractor, apply the photometric calibration 
  * if requested.
  */

  if (applyZeropoint) {

    cpl_msg_info(task, "Performing photometric calibration ...");

    if (!(VmImApplyPhot(rawImageFF, ipcTable))) {
      cpl_msg_error(task, "Photometric calibration failed!");

      deleteTable(ipcTable);

      if (astTable)
        deleteTable(astTable);

      deleteImage(rawImageFF);

      return EXIT_FAILURE;

    }

    deleteTable(ipcTable);

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

  if (!(galTable = VmImDetectObjects(rawImageFF, 0))) {
    cpl_msg_error(task, "Object detection failed!");

    if (astTable)
      deleteTable(astTable);

    deleteImage(rawImageFF);

    return EXIT_FAILURE;
  }


 /*
  * Select all stars among the detected objects and associate them 
  * to the entries in the astrometric table. 
  */

  cpl_msg_info(task, "Selecting stars from the list of detected objects ...");

  if (!(starTable = VmImBuildStarTable(galTable, starIndex, magLimit))) {
    cpl_msg_error(task, "Selection of stars failed!");
    
    if (astTable)
      deleteTable(astTable);

    deleteTable(galTable);
    deleteImage(rawImageFF);
    return EXIT_FAILURE;
  }

  deleteTable(galTable);
  /* writeFitsStarTable("starTable.TFITS", starTable) ; */
  writeFitsStarTable("starTable.fits", starTable) ;

 /* 
  * Before matching, I must convert into virtual pixels,
  * with whichever coeffs I have in header. But I also 
  * need to keep the original coordinates, these I save
  * in temporary vectors.
  */

  allStars = starTable->cols->len;

  star_orig = (int *)cpl_calloc(allStars, sizeof(int));
  ximage_orig = (double *)cpl_calloc(allStars, sizeof(double));
  yimage_orig = (double *)cpl_calloc(allStars, sizeof(double));


  star_tmp = tblGetIntData(starTable, "NUMBER");
  ximage_tmp = tblGetDoubleData(starTable, "X_IMAGE");
  yimage_tmp = tblGetDoubleData(starTable, "Y_IMAGE");

  for (i=0; i< allStars ; i++) {
    star_orig[i] = star_tmp[i];
    ximage_orig[i] = ximage_tmp[i];
    yimage_orig[i] = yimage_tmp[i];
  }

  cpl_msg_debug(task, "Applying SkyToCcd matrix, correcting for temperature "
              "effects...");
  
  if (computeVirtualPixels(rawImageFF->descs, starTable, tempCheck,
                           tempTolerance) == VM_FALSE) {
    cpl_msg_error(task, "Cannot apply SkyToCcd matrix!");
    return 0;
  }
  

  cpl_msg_info(task, "Associating detected stars with reference catalog ...");

  minStarNo = MAX((fitOrderX + 1) * (fitOrderX + 1),
                  (fitOrderY + 1) * (fitOrderY + 1));
  /*PDB redirected to alternative VmImBuildStarMatchTable_skyccd function
    here that contains changed indexing necessary for skyccd to work properly*/
  starMatchTable = VmImBuildStarMatchTable_skyccd(rawImageFF, starTable, astTable,
                                           minStarNo, searchRadius,
                                           magTolerance1, magTolerance2,
                                           sigmaClip);

  if (!starMatchTable) {
    cpl_msg_error(task, "Association failed! Cannot build Star Match table.");
    if (astTable)
      deleteTable(astTable);
    deleteTable(starTable);
    deleteImage(rawImageFF);
    return EXIT_FAILURE;
  }
  
  if (astTable)
    deleteTable(astTable);
  
  deleteTable(starTable);
  
  /* Save the star match product */
  starMatchProductTag = (char *)pilTrnGetCategory("ImgStarMatchAstrometry");
  vmstrlower(strcpy(starMatchProductName, starMatchProductTag));
  strcat(starMatchProductName, ".fits");
  if (createFitsTable(starMatchProductName, starMatchTable, 
                      starMatchProductTag) != VM_TRUE) {
    cpl_msg_error(task, "Cannot create local product file %s!",
                  starMatchProductName);

    deleteTable(starMatchTable);
    deleteImage(rawImageFF);

    return EXIT_FAILURE;
  }
  else {
    starMatchProductFrame = newPilFrame(starMatchProductName, 
                                        starMatchProductTag);

    pilFrmSetType(starMatchProductFrame , PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(starMatchProductFrame , PIL_FRAME_FORMAT_TABLE);
    pilFrmSetProductLevel(starMatchProductFrame , PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(starMatchProductFrame , PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, starMatchProductFrame );
  }
  
  if (starMatchTable->cols) {
      numStars = starMatchTable->cols->len;
  }
  else {
      numStars = 0;
      goto skip;
  }
  error = 1;
  
  cpl_msg_info(task, "Determining the conversion matrix");


 /*    
  * Procedure should be as follows: from starMatch Table we have 
  * xVirtual, yVirtual (i.e. corrected for temp effects) on image 
  * and RA DEC we must recover the original X and Y, then we should 
  * transform the RA DEC into Xvirtual, Yvirtual by using solely the 
  * WCS present in header of starMatch. Then compute a transformation 
  * from pixels to pixels as in MaskToCcd.
  */

  star = tblGetIntData(starMatchTable, "NUMBER");

 /*  
  * we transform the RA DEC into Xvirtual  Yvirtual
  * WARNING: this is done "in-place" thus X-IMAGE Y_IMAGE
  * are not waht they were before (I think...)
  */

  cpl_msg_debug(task, "Projecting astrometric stars into image plane!");

  wcstopix(starMatchTable->cols->len, starMatchTable, wcs);
  
  xVirtual = tblGetDoubleData(starMatchTable, "X_IMAGE");
  yVirtual = tblGetDoubleData(starMatchTable, "Y_IMAGE");
  //xyVirtual *were* measured pixel coords corrected for Temp and distortion
  // (i.e. where the centroid would be measured to have been if there were
  //no distortion...)
  // With the change that I made to vmimgextraction.c xyVirtual are now the
  // astrometric sky co-ords converted to pix using *only* the standard WCS
  // mapping (no distortion or temperature). This is what is needed for the
  // surface fit below.

 /* 
  * We recover the original sex coordinates 
  */

  ximage = (double *)cpl_calloc(numStars, sizeof(double));
  yimage = (double *)cpl_calloc(numStars, sizeof(double));

  for (i = 0; i < numStars; i++) {
    for (j = 0; j < allStars; j++) {
      if (star[i] == star_orig[j]) {
        ximage[i] = ximage_orig[j];
        yimage[i] = yimage_orig[j];
        break;
      }
    }
  }


 /* 
  * Compute transformation coefficients 
  * i.e. from real to Virtual Pixels, so called CCD To Sky
  */

  if (numStars < minStarNo) {
    cpl_msg_error(task, "Fit order is too high: "
                "%d stars, %d X degree polynomial %d Y degree polynomial", 
                numStars, fitOrderX, fitOrderY); 
    deleteTable(starMatchTable);
    return EXIT_FAILURE;
  }
  

 /*
  * Prepare pixel lists to fit 
  */

  noFound = numStars;
  surfaceX = newPixel(numStars);
  surfaceY = newPixel(numStars);

  for (i = 0; i < noFound; i++) {
    surfaceX[i].x = ximage[i];
    surfaceX[i].y = yimage[i];
    surfaceX[i].i = xVirtual[i];
    surfaceY[i].x = ximage[i];
    surfaceY[i].y = yimage[i];
    surfaceY[i].i = yVirtual[i];
  }

 /* 
  * we fit with a rejection on the points
  * that deviate from the fit  
  */

  noIter = 0;
  noFit = noFound;
  oldXrms = 99.99;
  oldYrms = 99.99;
  exclude = intVector(0, noFound);

  for (i = 0; i < noFound; i++)
    exclude[i] = 0;

  coeffControlString = createVimosCtrlStr(fitOrderX, fitOrderY); 

  do {

    coeffX = fitSurfacePolynomial(surfaceX, noFit, coeffControlString, 
                                  fitOrderX * 2, &dummy, &xrms);
    coeffY = fitSurfacePolynomial(surfaceY, noFit, coeffControlString, 
                                    fitOrderY * 2, &dummy, &yrms);

   /* 
    * Do we need to iterate some more ?? 
    */

    if (noIter < noIterations) {

     /* 
      * Measure differences between new and old positions, and reject
      * points that are 4-sigma outliers
      */

      deletePixel(newPos);
      /*PDB changed noFit to noFound below*/
      newPos = newPixel(noFound);
      for (n = 0; n < noFound; n++) {
	/*PDB added !exclude[n] condition below*/
        if (!exclude[n]) {
	  //printf("%lf %lf\n",ximage[n],yimage[n]);
	  k=0;
	  for (i = 0; i <= fitOrderX; i++) {
	    k=i;
	    for (j = 0; j <= fitOrderX; j++) {
	      newPos[n].x += coeffX[k] * ipow(ximage[n],j) * ipow(yimage[n],i); 
	      //printf("%d %lf %lf %lf %d %lf %d\n",noIter, newPos[n].x,coeffX[k],ximage[n],j,yimage[n],i);
	      k = k + fitOrderX + 1;
	    }
	  }
	  k=0;
	  for (i = 0; i <= fitOrderY; i++) {
	    k=i;
	    for (j = 0; j <= fitOrderY; j++) { 
	      newPos[n].y += coeffY[k] * ipow(ximage[n],j) * ipow(yimage[n],i); 
	      //printf("%d %lf %lf %lf %d %lf %d\n",noIter, newPos[n].y,coeffY[k],ximage[n],j,yimage[n],i);
	      k = k + fitOrderY + 1;
	    }
	  }
	}
      }

      for (i = 0; i < noFound; i++) {
	/* PDB added exclude[i] condition below */
        if (!exclude[i]) {
	  diffX = fabs(newPos[i].x - xVirtual[i]);
	  diffY = fabs(newPos[i].y - yVirtual[i]);
	  if (diffX > 4. * sqrt(xrms) || diffY > 4. * sqrt(yrms))
	    exclude[i] = 1;
	}
      }
      j = 0;
      for (i = 0; i < noFound; i++) {
        if (!exclude[i]) {
          surfaceX[j].x = ximage[i];
          surfaceX[j].y = yimage[i];
          surfaceX[j].i = xVirtual[i];
          surfaceY[j].x = ximage[i];
          surfaceY[j].y = yimage[i];
          surfaceY[j].i = yVirtual[i];
          j++;
        }
      }

      deltaXrms = oldXrms - sqrt(xrms);
      deltaYrms = oldYrms - sqrt(yrms);
      oldXrms = sqrt(xrms);
      oldYrms = sqrt(yrms);
      noReject = noFit - j;
      noFit = j;
      noIter++;

      if (noReject == 0) {

       /* 
        * If no points are rejected we don't need any more iterations 
        */
        nextCycle = 0;
      }
      else if (deltaXrms < (0.1 * oldXrms) && deltaYrms < (0.1 * oldYrms)) {

       /* 
        * If the rms doesn't decrease any more, we can stop with a last
        * iteration, to make sure we use the latest rejections.
        */
        noIter = noIterations;
      }
    }
    else {

     /* 
      * Otherwise we just stop the iterations cycle 
      */
      nextCycle = 0;
    }
  } while (nextCycle);

  cpl_msg_info (task, "Number of fitted points %d", noFit);

  
  if ( (xrms >= 1)  || (yrms >= 1) ) {
    cpl_msg_debug(task, "The error in the object position is > 1 pixel: "
               "The object coordinates on the CCD are determined with "
               "accuracy XRMS = %f, YRMS = %f", sqrt(xrms), sqrt(yrms));
  }
  


  writeDoubleDescriptor(&rawImageFF->descs, pilTrnGetKeyword("CcdSkyXrms"),
                        sqrt(xrms), ""); 
  writeDoubleDescriptor(&rawImageFF->descs, pilTrnGetKeyword("CcdSkyYrms"),
                        sqrt(yrms), "");
  
  writeIntDescriptor(&(rawImageFF->descs), pilTrnGetKeyword("CcdSkyXord"),
                     fitOrderX,  "");
  writeIntDescriptor(&(rawImageFF->descs), pilTrnGetKeyword("CcdSkyYord"),
                     fitOrderY,"");
   
 
  k = 0;
  for (i = 0; i <= fitOrderX; i++) {
    k = i;
    for (j = 0; j <= fitOrderX; j++) {
      sprintf(descStringVal, "%.8G", coeffX[k]);
      writeStringDescriptor((&(rawImageFF->descs)),
                       pilTrnGetKeyword("CcdSkyX", i, j), descStringVal, "");
      k = k + fitOrderX + 1;
    }
  }
  
  k = 0;
  for (i = 0; i <= fitOrderY; i++) {
    k = i;
    for (j = 0; j <= fitOrderY; j++) {
      sprintf(descStringVal, "%.8G", coeffY[k]);
      writeStringDescriptor((&(rawImageFF->descs)),
            pilTrnGetKeyword("CcdSkyY", i, j), descStringVal, "");
      k = k + fitOrderX + 1;
    }
  } 
 
 /*
  * Compute a transformation from Virtual pixels to Real pixels
  * so called SkyToCcd.
  * IMPORTANT NOTE (BG)  Sky To CCD matrix is NEVER used. It is 
  * just computed and written to PAF. We might want to drop the 
  * writing at a certain point.
  */

  j = 0;

  for (i = 0; i < noFound; i++) {
    if (!exclude[i]) {
      surfaceX[j].x = xVirtual[i];
      surfaceX[j].y = yVirtual[i];
      surfaceX[j].i = ximage[i];
      surfaceY[j].x = xVirtual[i];
      surfaceY[j].y = yVirtual[i];
      surfaceY[j].i = yimage[i];
      j++;
    }
  }

  coeffX = fitSurfacePolynomial(surfaceX, noFit, coeffControlString, 
                                fitOrderX * 2, &dummy, &xrms);
  coeffY = fitSurfacePolynomial(surfaceY, noFit, coeffControlString, 
                                fitOrderY * 2, &dummy, &yrms);
 

  writeDoubleDescriptor(&rawImageFF->descs, pilTrnGetKeyword("SkyCcdXrms"),
                        sqrt(xrms), ""); 
  writeDoubleDescriptor(&rawImageFF->descs, pilTrnGetKeyword("SkyCcdYrms"),
                        sqrt(yrms), "");
  
  writeIntDescriptor(&(rawImageFF->descs), pilTrnGetKeyword("SkyCcdXord"),
                     fitOrderX,  "");
  writeIntDescriptor(&(rawImageFF->descs), pilTrnGetKeyword("SkyCcdYord"),
                     fitOrderY,"");
 
  k = 0;
  for (i = 0; i <= fitOrderX; i++) {
    k = i;
    for (j = 0; j <= fitOrderX; j++) {
      sprintf(descStringVal, "%.8G", coeffX[k]);
      writeStringDescriptor((&(rawImageFF->descs)),
                       pilTrnGetKeyword("SkyCcdX", i, j), descStringVal, "");
      k = k + fitOrderX + 1;
    }
  }
  
  k = 0;
  for (i = 0; i <= fitOrderY; i++) {
    k = i;
    for (j = 0; j <= fitOrderY; j++) {
      sprintf(descStringVal, "%.8G", coeffY[k]);
      writeStringDescriptor((&(rawImageFF->descs)),
            pilTrnGetKeyword("SkyCcdY", i, j), descStringVal, "");
      k = k + fitOrderY + 1;
    }
  } 
 
  cpl_free(coeffControlString);

 /*
  * Setting up the output PAF file  
  */

  createCcdSkyPAF(rawImageFF->descs, namePAF, &pafFileName);

  if (createCcdSkyPAF(rawImageFF->descs, pipeNamePAF, &pafFileName) 
      == EXIT_SUCCESS) {
    productFrame = newPilFrame(pafFileName, pilTrnGetCategory("PAFCategory"));

    pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);  
    pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_PAF);
    pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, productFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", pafFileName);

    deleteImage(rawImageFF);
    deleteTable(starMatchTable);

    return EXIT_FAILURE;
  }

skip:

 /*
  * Reduced image
  */

  productTag = reducedTag;
  vmstrlower(strcpy(productName, productTag));
  strcat(productName, ".fits");
  insertDoubleDescriptor(&rawImageFF->descs, pilTrnGetKeyword("DataMin"),
                         imageMinimum(rawImageFF), pilTrnGetComment("DataMin"),
                         "ESO*", 1);
  insertDoubleDescriptor(&rawImageFF->descs, pilTrnGetKeyword("DataMax"),
                         imageMaximum(rawImageFF), pilTrnGetComment("DataMax"),
                         "ESO*", 1);
  writeDoubleDescriptor(&rawImageFF->descs, pilTrnGetKeyword("DataMean"),
                        imageMean(rawImageFF), pilTrnGetComment("DataMean"));
  writeDoubleDescriptor(&rawImageFF->descs,
                        pilTrnGetKeyword("DataStdDeviation"), 
                        imageSigma(rawImageFF), 
                        pilTrnGetComment("DataStdDeviation"));
  writeDoubleDescriptor(&rawImageFF->descs, pilTrnGetKeyword("DataMedian"),
                        imageMedian(rawImageFF),
                        pilTrnGetComment("DataMedian"));
  writeIntDescriptor(&rawImageFF->descs, pilTrnGetKeyword("MatchNstars"),
                     numStars,
                     pilTrnGetComment("MatchNstars"));
  writeIntDescriptor(&rawImageFF->descs, pilTrnGetKeyword("MatrixFittedPoints"),
                     noFit,
                     pilTrnGetComment("MatrixFittedPoints"));
  insertStringDescriptor(&rawImageFF->descs, pilTrnGetKeyword("DoCategory"),
                         productTag, pilTrnGetKeyword("DoCategory"),"ESO*", 1);
  deleteSetOfDescriptors(&rawImageFF->descs, "ESO DPR*");
  //These keywords shouldn't be present in the input raw, but we nevertheless
  //remove them in case they are present (See PIPE-5422)
  deleteSetOfDescriptors(&(rawImageFF->descs), "PCOUNT");  
  deleteSetOfDescriptors(&(rawImageFF->descs), "GCOUNT");  
    
  if (createFitsImage(productName, rawImageFF, productTag) != VM_TRUE) {
    cpl_msg_error(task, "Cannot create local product file %s!", productName);

    deleteTable(starMatchTable);
    deleteImage(rawImageFF);

    return EXIT_FAILURE;
  } 
  else {
    productFrame = newPilFrame(productName, productTag);

    pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, productFrame);
  }

  deleteImage(rawImageFF);
  deleteTable(starMatchTable);

  return EXIT_SUCCESS;

}


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
                    "vmskyccd",
    "Determine instrument distortions from an astrometric field exposure.",
    "This recipe is used to determine the CCD to Sky distortion and its\n"
    "inverse. SExtractor v2.1.6 is used for detecting in the field of\n" 
    "view the objects matching the entries of an astrometric catalog. The\n" 
    "deviations from the theoretical relation between CCD and celestial\n"
    "coordinates (WCS) are then determined and modelled.\n\n"
    "Input files:\n\n"
    "  DO category:              Type:       Explanation:         Required:\n"
    "  IMG_ASTROMETRY            Raw         Astrometric field       Y\n"
    "  MASTER_BIAS               Calib       Master bias             Y\n"
    "  MASTER_DARK               Calib       Master dark             .\n"
    "  IMG_MASTER_SKY_FLAT       Calib       Master flat field       Y\n"
    "  ASTROMETRIC_TABLE         Calib       Astrometric catalog     Y\n"
    "  PHOT_COEFF_TABLE          Calib       Photometric table       .\n"
    "  CCD_TABLE                 Calib       Bad pixel table         .\n\n"
    "The only product of this recipe is a PAF file, copied (or moved) to\n"
    "the product directory, that is identical to the produced Instrument\n"
    "WorkStation configuration file IMG_sky2ccd_Q.cmf (where Q indicates\n"
    "the VIMOS quadrant number) that is created in the same directory where\n"
    "the recipe is launched. A CCD table must be specified in input only\n"
    "if a bad pixel cleaning is requested.\n\n"
    "For more details, please refer to the VIMOS Pipeline User's Guide.",

    "ESO VIMOS Pipeline Team and VIMOS Consortium",

    PACKAGE_BUGREPORT,

    "This file is part of the VIMOS Instrument Pipeline\n"
    "Copyright (C) 2002-2005 European Southern Observatory\n\n"
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
    "Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA\n",

                    vmskyccd_create,
                    vmskyccd_exec,
                    vmskyccd_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
