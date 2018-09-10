/* $Id: vmmosobsjitter.c,v 1.5 2012-10-30 10:36:05 cgarcia Exp $
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
 * $Date: 2012-10-30 10:36:05 $
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <string.h>
#include <math.h>

#include <cxmemory.h>

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
#include <pilutils.h>
#include <pilfits.h>

#include "vmimage.h"
#include "vmimagearray.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmgrismtable.h"
#include "vmextractiontable.h"
#include "vmobjecttable.h"
#include "vmwindowtable.h"
#include "vmextincttable.h"
#include "vmspecphottable.h"
#include "vmfit.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmmosflat.h"
#include "vmmoswavecalib.h"
#include "vmmossky.h"
#include "vmmosfringes.h"
#include "vmmosextraction.h"
#include "vmmossphotcalib.h"
#include "vmmosutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"

static int vmmosobsjitter_create(cpl_plugin *);
static int vmmosobsjitter_exec(cpl_plugin *);
static int vmmosobsjitter_destroy(cpl_plugin *);
static cxint vmmosobsjitter(PilSetOfFrames *);


static char vmmosobsjitter_description[] =
"This recipe is used to apply basic reduction steps to a sequence of\n"
"exposures made in MOS mode, to combine them in a single image, to\n"
"locate objects, and to optimally extract their spectra. Each input\n"
"image is processed in the same way as by recipe vmmosobsstare,\n"
"therefore what mainly characterises the vmmosobsjitter recipe is\n"
"the task of combining the input frames.\n\n"
"Input files:\n\n"
"  DO category:             Type:       Explanation:          Required:\n"
"  MOS_SCIENCE              Raw         Science exposure          Y\n"
"  MASTER_BIAS              Calib       Master bias               Y\n"
"  MASTER_DARK              Calib       Master dark               .\n"
"  MOS_MASTER_SCREEN_FLAT   Calib       Normalised flat field     .\n"
"  EXTRACT_TABLE            Calib       Extraction table          .\n"
"  GRISM_TABLE              Calib       Grism table               Y\n"
"  MOS_FRINGES_SKY          Calib       Sky+fringes map ('Raw')   .\n"
"  MOS_FRINGES              Calib       Fringe map ('Resampled')  .\n"
"  EXTINCT_TABLE            Calib       Atmospheric extinction    .\n"
"  MOS_SPECPHOT_TABLE       Calib       Response curve            .\n"
"  CCD_TABLE                Calib       Bad pixel table           .\n\n"
"Output files:\n\n"
"  DO category:             Data type:  Explanation:\n"
"  MOS_SCIENCE_REDUCED      FITS image  Extracted objects spectra\n"
"  MOS_SCIENCE_FLUX_REDUCED FITS image  Flux calibrated objects spectra\n"
"  MOS_SCIENCE_EXTRACTED    FITS image  Sky subtracted slit spectra\n"
"  MOS_SCIENCE_SKY          FITS image  Sky slit spectra\n"
"  MOS_SKY_REDUCED          FITS image  Extracted sky spectra\n"
"  MOS_FRINGES_SKY          FITS image  Sky+fringes map (method 'Raw')\n"
"  MOS_FRINGES              FITS image  Fringe map (method 'Resampled')\n"
"  OBJECT_TABLE             FITS table  Objects spectra identification\n"
"  WINDOW_TABLE             FITS table  Objects positions in slit\n\n"
"A flat field correction is applied only if a normalised master flat\n"
"field (produced by the recipe vmspflat) is specified.\n\n" 
"The extraction table is the product of the local spectral distortions\n"
"modelling performed by the recipe vmspcaldisp. If an extraction table\n"
"is not specified, then the global distortion models read from the\n" 
"science frame header are used.\n\n"
"The grism table contains necessary information to control the way\n" 
"spectra are extracted, starting from the reference wavelength (header\n" 
"entry PRO WLEN CEN), on a specific range of pixels above and below\n"
"its position on the CCD (header entries PRO SPECT LLEN LO and PRO\n" 
"SPECT LLEN HI). Other parameters, used in the extraction of the\n"
"science slit spectra, are the start and the end wavelength of the\n" 
"image of the extracted slit spectra (header entries PRO WLEN START\n" 
"and PRO WLEN END), and the step of the sampling along the dispersion\n" 
"direction (header entry PRO WLEN INC). Finally, the wavelengths of the\n" 
"sky lines used in the alignment of the spectral distortion models,\n" 
"necessary to keep into account the possible coordinates shifts\n" 
"introduced by a variation of the instrument flexures between the\n" 
"science and the calibration exposures, are listed in the header\n" 
"keywords PRO SKY WLENi, with i ranging from 1 to the number specified\n" 
"in the keyword PRO SKY NO.\n\n"
"A CCD table must be specified in input only if a bad pixel cleaning is\n"
"requested.\n\n"
"The slit spectra are remapped with the instrument distortions removed\n"
"and at a fixed wavelength step. A sky value is estimated for each\n"
"wavelength and then subtracted from the data. The result is stored\n"
"in the MOS_SCIENCE_EXTRACTED image, while the image MOS_SCIENCE_SKY\n"
"contains the subtracted sky model. The 1D extracted spectra are stored\n" 
"in the MOS_SCIENCE_REDUCED image, while the corresponding sky spectra\n"
"extracted with the same method are stored in the MOS_SKY_REDUCED image.\n"
"The positions of the extracted slit spectra and of the detected\n"
"objects that they may contain are listed in the window table.\n\n"
"If a spectro-photometric table (produced by the recipe vmmosstandard)\n"
"is specified together with an atmospheric extinction table, and a flux\n"
"calibration is requested, then a MOS_SCIENCE_FLUX_REDUCED image is also\n"
"created. This image is identical to the MOS_SCIENCE_REDUCED, but the\n"
"spectra it contains are flux calibrated, and expressed in units of\n"
"erg/cm/cm/s/Angstrom.\n\n"
"Sky fringes may be subtracted from the data using two different\n" 
"methods. The 'Raw' method will subtract the median sky+fringes\n"
"pattern from each input image, and then fit away the residuals\n"
"possibly introduced by sky variations between exposures. The\n" 
"'Resampled' method will fit the sky first, and then subtract the\n"
"median image of the residuals (stored in the product MOS_FRINGES).\n"
"For more details, please refer to the VIMOS Pipeline User's Guide.";     


#define MAX_COMMENT_LENGTH (80)

/*
 * Definition of the label strings for all methods the recipe function
 * supports for combining frames and bias removal, with their associated
 * method code.
 */

static const char *combMethodNames[] = {
  "Auto",
  "Ksigma",
  "MinMax",
  "Median",
  "Average"
};

static const CombMethod combMethods[] = {
  COMB_AUTO,
  COMB_KSIGMA,
  COMB_REJECT,
  COMB_MEDIAN,
  COMB_AVERAGE
};

static unsigned int nCombMethods = sizeof(combMethods) / sizeof(CombMethod);

static const char *biasMethodNames[] = {
  "Master",
  "Zmaster"
};

static const BiasMethod biasMethods[] = {
  BIAS_MASTER,
  BIAS_ZMASTER
};

static unsigned int nBiasMethods = sizeof(biasMethods) / sizeof(BiasMethod);


static const char *skyMethodNames[] = {
  "Median",
  "Fit"
};

static const SkyMethod skyMethods[] = {
  SKY_MEDIAN,
  SKY_FIT
};

static unsigned int nSkyMethods = sizeof(skyMethods) / sizeof(SkyMethod);


typedef VimosSpecSampleOption SamplingMethod;   /* Just an alias */

#ifdef ONLINE_MODE

static const char *samplingMethodNames[] = {
  "Linear",
  "Log"
};

static const SamplingMethod samplingMethods[] = {
  VM_SP_LIN_LAMBDA,
  VM_SP_LOG_LAMBDA
};

static unsigned int nSamplingMethods
                         = sizeof(samplingMethods) / sizeof(SamplingMethod);

#endif



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
                    "vmmosobsjitter",
                    "Reduce and combine a set of jittered MOS observations.",                    
                    vmmosobsjitter_description,
                    "ESO VIMOS Pipeline Team and VIMOS Consortium",
                    PACKAGE_BUGREPORT,
                    vimos_get_license(),
                    vmmosobsjitter_create,
                    vmmosobsjitter_exec,
                    vmmosobsjitter_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}

/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmmosobsjitter_create(cpl_plugin *plugin)
{

    cpl_recipe    *recipe;
    cpl_parameter *p;

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


    p = cpl_parameter_new_enum("vimos.Parameters.stacking.method",
                               CPL_TYPE_STRING,
                               "Frames combination method",
                               "vimos.Parameters",
                               "Average", 5, "Average", "Median", "MinMax",
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


    p = cpl_parameter_new_enum("vimos.Parameters.sky.method",
                                CPL_TYPE_STRING,
                                "Sky level determination method.",
                                "vimos.Parameters",
                                "Median", 2, "Fit", "Median");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SkyMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SkyMethod");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.sky.order",
                                CPL_TYPE_INT,
                                "Degree of polynomial used when the "
                                "SkyMethod is set to Fit.",
                                "vimos.Parameters",
                                2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PolyOrder");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PolyOrder");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.sky.ksigma.low",
                                CPL_TYPE_DOUBLE,
                                "Low threshold for K-sigma rejection "
                                "in sky fitting.",
                                "vimos.Parameters",
                                1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SkyKSigmaLow");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SkyKSigmaLow");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.sky.ksigma.high",
                                CPL_TYPE_DOUBLE,
                                "High threshold for K-sigma rejection "
                                "in sky fitting.",
                                "vimos.Parameters",
                                1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SkyKSigmaHigh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SkyKSigmaHigh");
    cpl_parameterlist_append(recipe->parameters, p);


#ifdef ONLINE_MODE

    p = cpl_parameter_new_enum("vimos.Parameters.sampling",
                                CPL_TYPE_STRING,
                                "Spectrum sampling in wavelength "
                                "during extraction.",
                                "vimos.Parameters",
                                "Linear", 2, "Linear", "Log");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SamplingMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SamplingMethod");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

    p = cpl_parameter_new_value("vimos.Parameters.extraction.fuzz",
                                CPL_TYPE_INT,
                                "Extra pixels from expected position "
                                "of spectrum edge in extraction.",
                                "vimos.Parameters",
                                5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "Fuzz");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "Fuzz");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.detection.exclude",
                                CPL_TYPE_INT,
                                "Number of excluded pixels at slit ends in "
                                "object search or in sky level determination.",
                                "vimos.Parameters",
                                2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SlitMargin");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SlitMargin");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.slit.tolerance",
                                CPL_TYPE_DOUBLE,
                                "Tolerance for drift of slit on the CCD, "
                                "used to tell long slits from short ones.",
                                "vimos.Parameters",
                                0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SlitTolerance");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SlitTolerance");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

    p = cpl_parameter_new_value("vimos.Parameters.sky.linewidth",
                                CPL_TYPE_INT,
                                "Size of spectrum to extract around "
                                "any skyline.",
                                "vimos.Parameters",
                                16);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LineWidth");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LineWidth");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.detection.sigma",
                                CPL_TYPE_DOUBLE,
                                "Object detection level in units of sigma.",
                                "vimos.Parameters",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "DetectionLevel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "DetectionLevel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.detection.levels",
                                CPL_TYPE_INT,
                                "Number of levels in the watershed method "
                                "in object detection.",
                                "vimos.Parameters",
                                32);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "WatershedLevels");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "WatershedLevels");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.detection.fraction",
                                CPL_TYPE_DOUBLE,
                                "Flux fraction to use in watershed.",
                                "vimos.Parameters",
                                0.01);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "WatershedFraction");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "WatershedFraction");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.detection.minsize",
                                CPL_TYPE_INT,
                                "Minimal size for an object candidate to "
                                "be considered an object.",
                                "vimos.Parameters",
                                2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MinObjectSize");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MinObjectSize");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.detection.maxsize",
                                CPL_TYPE_INT,
                                "Maximal size for an object candidate for "
                                "not trying deblend into sub-objects.",
                                "vimos.Parameters",
                                7);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MaxObjectSize");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MaxObjectSize");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.badpixel.clean",
                                CPL_TYPE_BOOL,
                                "Bad pixel correction on MOS science "
                                "exposure.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flux.calibration",
                                CPL_TYPE_BOOL,
                                "Extracted spectra are flux calibrated.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CalibrateFlux");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CalibrateFlux");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.fringing",
                                CPL_TYPE_BOOL,
                                "Apply fringing corrections.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "FringingCorr");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "FringingCorr");
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_enum("vimos.Parameters.fringing.removing.method",
                                CPL_TYPE_STRING,
                                "Sky fringes removal method.",
                                "vimos.Parameters",
                                "Raw", 2, "Raw", "Resampled");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "FringingMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "FringingMethod");
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.Parameters.fringing.offset",
                                CPL_TYPE_DOUBLE,
                                "Minimum required offset between exposures "
                                "for applying the sky fringing correction.",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "FringingOffset");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "FringingOffset");
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.Parameters.extraction.optimal",
                                CPL_TYPE_BOOL,
                                "Use 1D Horne extraction",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "HorneExtraction");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "HorneExtraction");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.slit.model",
                                CPL_TYPE_BOOL,
                                "Model wavelength solution within each slit.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ModelSlit");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ModelSlit");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.slit.order",
                                CPL_TYPE_INT,
                                "Order of polynomial for wavelength solution "
                                "modeling within each slit.",
                                "vimos.Parameters",
                                0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ModelSlitOrder");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ModelSlitOrder");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.sky.align",
                                CPL_TYPE_BOOL,
                                "Use sky lines to refine the wavelength "
                                "calibration",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "UseSkylines");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "UseSkylines");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.save.intermediate",
                                CPL_TYPE_BOOL,
                                "Save intermediate reduction steps",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SaveIntermediate");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SaveIntermediate");
    cpl_parameterlist_append(recipe->parameters, p);


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
vmmosobsjitter_exec(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;

    cxint status = 0;

    PilSetOfFrames *sof = NULL;


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

    if (vmmosobsjitter(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmmosobsjitter");
        
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
vmmosobsjitter_destroy(cpl_plugin *plugin)
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
 *   Combine spectral data taken in shift-and-stare mode.
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param sof  Set of frames with references to at least two MOS 
 *             exposures, a master bias, a grism table, and a CCD
 *             table. A spectral master flat field, a master dark, 
 *             a spectrophotometric table, and an extraction table, 
 *             may be optionally given. If an extraction table is 
 *             specified, all the spectral distorsion model are 
 *             read from its header, otherwise they are read from 
 *             the header of the science exposures (first guesses). 
 *             If bad pixel cleaning is requested, a CCD table must 
 *             be given.
 *
 * @doc 
 *   Reduce a sequence of at least two MOS exposures taken with the same 
 *   mask. Reduce each exposure as in VmMosObsSingle, but all the 
 *   extracted 2D spectra are stacked before applying the optimal
 *   extraction of the detected objects and creating the object table.
 *   The spectrophotometric calibration on 1D and 2D spectra may be 
 *   finally applied.
 *
 *   \begin{itemize}
 *
 *     \item BiasMethod:         Method for bias removal from MOS
 *                               frame. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Master:             Prescan and overscan regions are trimmed
 *                                 away from the MOS image after the master
 *                                 bias removal.
 *
 *       \item Zmaster:            After master bias removal the residual
 *                                 signal found in the MOS image overscan
 *                                 regions is modelled and subtracted.
 *                                 Next, prescan and overscan regions
 *                                 are trimmed away.
 *     \end{itemize}
 *
 *     \item SkyMethod:          Method used for finding the sky level
 *                               to subtract at each wavelength. Legal
 *                               settings are:
 *
 *     \begin{itemize}
 *
 *       \item Median:             Median of sky level.
 *
 *       \item Fit:                Polynomial fitting of sky.
 *
 *     \end{itemize}
 *
 *     \item PolyOrder:          Degree of the polynomial used when the
 *                               SkyMethod is set to Fit.
 *
 *     \item Sampling:           How the spectrum is sampled in wavelength
 *                               during extraction. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Linear:             Linear sampling
 *
 *       \item Log:                Logarithmic sampling
 *
 *     \end{itemize}
 *
 *     \item Fuzz:               Extra number of pixels to be used in
 *                               spectra extraction.
 *
 *     \item SlitMargin:         Number of pixels at slit ends to be
 *                               excluded from object search, or from
 *                               sky level determination.
 *
 *     \item SlitTolerance:      Tolerance for drift of slit on the CCD.
 *                               If the position of the slit on the CCD
 *                               is straight within tolerance pixels,
 *                               the slit is classified as short,
 *                               otherwise as long.
 *
 *     \item LineWidth:          Size of spectrum extracted around any
 *                               skyline.
 *
 *     \item DetectionLevel:     Object detection level in units of sigma.
 *
 *     \item WatershedLevels:    Number of levels in the watershed
 *                               method used to see if more than one
 *                               object is present in a contiguous
 *                               region of pixels above the threshold.
 *
 *     \item WatershedFraction   Flux fraction to use in watershed.
 *
 *     \item MinObjectSize:      Minimal size for an object candidate to
 *                               be considered an object.
 *
 *     \item MaxObjectSize:      Maximal size for an object candidate for
 *                               NOT trying to deblend it into sub-objects.
 *                               If a candidate of size smaller than
 *                               minObjectSize is found, it is taken as
 *                               one object.
 *
 *     \item CleanBadPixel:      Bad pixel correction on arc lamp frame.
 *                               If this option is turned on, the recipe
 *                               expects a CCD\_TABLE in the input set
 *                               of frames.
 *
 *     \item FringingCorr:       Fringing correction on input frames.
 *
 *     \item FringingMethod:     Fringing correction method. Legal 
 *                               settings are:
 *
 *     \begin{itemize}
 *
 *       \item Raw:                The median image obtained from the bias
 *                                 subtracted input images is subtracted
 *                                 from each input image. In this way both
 *                                 sky level and sky fringes are eliminated
 *                                 from the data. Possible residuals due to
 *                                 sky and/or fringes variability are then
 *                                 removed from the data according to the
 *                                 setting of the parameter SkyMethod.
 *
 *       \item Resampled:          The median image obtained from the images 
 *                                 containing the rectified and sky subtracted 
 *                                 spectra is subtracted from each of the 
 *                                 images.
 *     \end{itemize}
 *
 *
 *     \item HorneExtract:       1D extraction using Horne algorithm.
 *                               Alternatively a simple sum is used.
 *
 *     \item ModelSlit           Modeling of wavelength solution within
 *                               each slit (local global model).
 *
 *     \item ModelSlitOrder      Order of polynomial for modeling of
 *                               wavelength solution within each slit.
 *
 *   \end{itemize}
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 *
 * @see VmMosObsSingle()
 * 
 * @author P.Franzetti, P. Sartoretti, C. Izzo, R. Palsa
 */   

static cxint
vmmosobsjitter(PilSetOfFrames *sof)
{
  const char   task[]               = "vmmosobsjitter";

  const char  parameter[]           = "Parameters";

  const char *mosCategory          = pilTrnGetCategory("MosScience");
  const char *mosExtractedCategory = pilTrnGetCategory("MosScienceExtracted");
  const char *mosSkyCategory       = pilTrnGetCategory("MosScienceSky");
  const char *mosReducedCategory   = pilTrnGetCategory("MosScienceReduced");
  const char *mosFluxReducedCategory
                                   = pilTrnGetCategory("MosScienceFluxReduced");
  const char *mosSkyReducedCategory= pilTrnGetCategory("MosSkyReduced");
  const char *mosFringesSkyCategory= pilTrnGetCategory("MosFringesSky");
  const char *mosFringesCategory   = pilTrnGetCategory("MosFringes");
  const char *flatCategory         = pilTrnGetCategory("MosMasterScreenFlat");
  const char *sphotTableCategory   = pilTrnGetCategory("MosSphotTable");
  const char *biasCategory         = pilTrnGetCategory("MasterBias");
  const char *darkCategory         = pilTrnGetCategory("MasterDark");
  const char *ccdTableCategory     = pilTrnGetCategory("CcdTable");
  const char *grismTableCategory   = pilTrnGetCategory("GrismTable");
  const char *objectTableCategory  = pilTrnGetCategory("ObjectTable");
  const char *windowTableCategory  = pilTrnGetCategory("WindowTable");
  const char *extrTableCategory    = pilTrnGetCategory("ExtractTable");
  const char *atmTableCategory     = pilTrnGetCategory("ExtinctTable");

/*+*/
  const char  *lineCatalogCategory = pilTrnGetCategory("LineCatalog");


  char                  *biasMethodTag = NULL;
  BiasMethod             biasMethod = BIAS_UNDEF;
  int                    biasMethodEntry;

  char                  *combMethodTag = NULL;
  int                    combMethodEntry;

  CombMethod             combMethod = COMB_MEDIAN;
  CombParameters         combParameter;

  char                  *skyMethodTag = NULL;
  SkyMethod              skyMethod = SKY_UNDEF;
  int                    skyMethodEntry;
 
#ifdef ONLINE_MODE

  char                  *samplingMethodTag = NULL;
  int                    samplingMethodEntry;

#endif

  SamplingMethod         samplingMethod = VM_SP_LIN_LAMBDA;

  char                   output2DName[PATHNAME_MAX + 1];
  char                   output2DSkyName[PATHNAME_MAX + 1];
  char                   output1DName[PATHNAME_MAX + 1];
  char                   outputFringesName[PATHNAME_MAX + 1];
  char                   objectTableName[PATHNAME_MAX + 1];
  char                   windowTableName[PATHNAME_MAX + 1];

  char                   tmpName[PATHNAME_MAX + 1];

  VimosBool              updateOK = VM_TRUE;

  unsigned int           cleanBadPixel;
  unsigned int           calibrateFlux;
  unsigned int           fringingCorr;
  char                  *fringMethodTag = NULL;
  int                    fringMethod = 0;
  float                  skySigmaLow, skySigmaHigh;
  float                  minFringOffset;
  unsigned int           horneExtraction;
  unsigned int           modelSlit;
  unsigned int           useSkylines;
  unsigned int           saveIntermediate;
  unsigned int           error;

  int                    i, j;

  int                    order               = 0;
  int                    polyDeg             = 0;
  int                    fuzz                = 0;
  int                    slitMargin          = 0;
  int                    lineWidth           = 0;
  int                    numLevels           = 0;
  int                    minObjectSize       = 0;
  int                    minCompositeSize    = 0;
  int                    mosCount            = 0;
  int                    friCount            = 0;
  int                    minFrames           = 0;

  float                  slitTolerance       = 0.0;
  float                  detLevel            = 0.0;
  float                  objFrac             = 0.0;
  float                  specFrac            = 0.0;

  double                 expTime, sumExpTime; 
  double                 refAlpha, refDelta, cosDelta, dValue, pixScale;
  float                  minOffset;
  float                 *aOffset;
  float                 *dOffset;
  float                 *offset;

  double         cdelt, crval;

  PilFrame              *biasFrame, *darkFrame, *ccdFrame, *grismFrame;
  PilFrame              *mosFrame, *flatFrame, *extrFrame, *sphotFrame;
  PilFrame              *fringFrame, *outputFrame;
  PilFrame              *atmFrame;

/*+*/
  PilFrame              *lineCatFrame;
  VimosImage            *lineCatFile;
  VimosTable            *lineCat = NULL;
  double                 rmsValue, meanFwhm = 0.0;

  VimosImage            *grismFile           = NULL;
  VimosImage            *sphotFile           = NULL;
  VimosImage            *atmFile             = NULL;
  VimosImage            *extrFile            = NULL; 
  VimosImage            *objectFile          = NULL;
  VimosImage            *windowFile          = NULL;
  VimosImage            *biasImage           = NULL;
  VimosImage            *darkImage           = NULL;
  VimosImage            *flatImage           = NULL;
  VimosImage           **imaSpEx1D           = NULL;
  VimosImage            *imaSpEx1DCal        = NULL;
  VimosImage           **mosList             = NULL;
  VimosImage           **mosList2D           = NULL;
  VimosImage           **skyList2D           = NULL;
  VimosImage           **outSpSkyFra         = NULL;
  VimosImage           **outSpEx2D           = NULL;
  VimosImage           **outSpSkyEx          = NULL;
  VimosImage           **outSpSkyExStack     = NULL;
  VimosImage           **tempo               = NULL;
  VimosImage            *aritempo[2];
  VimosImage            *copySpStack2D;
  VimosImage            *tmpImage;
  VimosImage            *imaFringes          = NULL;
  VimosImage            *imaFringes2D        = NULL;
  VimosTable            *ccdTable            = NULL;
  VimosTable            *grismTable          = NULL;
  VimosTable            *sphotTable          = NULL;
  VimosTable            *atmTable            = NULL;

  VimosFloatArray       *offSets             = NULL;
  VimosIntArray         *twoDMap             = NULL;

  VimosExtractionTable  *combExtractionTable = NULL;
  VimosExtractionTable **extTablesList       = NULL;
  VimosWindowTable      *combWindowTable     = NULL;
  VimosWindowTable     **winTablesList       = NULL;
  VimosObjectTable      *objectTable         = NULL;


 /* FIXME:
  * This parameter is still unused:
  */

  float          limFrac = 0.0;


 /*
  * Get task parameters from the recipe database
  */

 /*
  * Determine the frame bias removal method.
  */

  biasMethodTag = (char *)pilDfsDbGetString(parameter, "BiasMethod");

  if ((biasMethodEntry =
        strselect(biasMethodTag, biasMethodNames, nBiasMethods)) < 0) {
    cpl_msg_error(task, "%s: Invalid bias removal method.", biasMethodTag);
    return EXIT_FAILURE;
  }

  biasMethod = biasMethods[biasMethodEntry];

 /*
  * Determine the frame stacking method and all method dependent
  * parameters.
  */

  combMethodTag = (char *)pilDfsDbGetString(parameter, "StackMethod");

  if ((combMethodEntry =
        strselect(combMethodTag, combMethodNames, nCombMethods)) < 0) {
    cpl_msg_error(task, "%s: Invalid frame combination method.", combMethodTag);
    return EXIT_FAILURE;
  }

  combMethod = combMethods[combMethodEntry];

  switch (combMethod) {
    
  case COMB_KSIGMA:

    minFrames = MIN_FRAMES_KSIGMA;
    combParameter.kSigmaLow =
                  pilDfsDbGetDouble(parameter, "KSigmaLow", 5.0);
    combParameter.kSigmaHigh =
                  pilDfsDbGetDouble(parameter, "KSigmaHigh", 5.0);
    break;

  case COMB_REJECT:

    combParameter.minRejection =
                  pilDfsDbGetDouble(parameter, "MinRejection", 0.0);
    combParameter.maxRejection =
                  pilDfsDbGetDouble(parameter, "MaxRejection", 0.0);

    minFrames = combParameter.minRejection + combParameter.maxRejection + 1;
    break;

  case COMB_MEDIAN:
    minFrames = MIN_FRAMES_MEDIAN;
    break;

  case COMB_AVERAGE:
    minFrames = MIN_FRAMES_AVERAGE;
    break;

  default:
    cpl_msg_warning(task, "Invalid stacking method. Using default "
                  "method: Median");
    combMethod = COMB_MEDIAN;
    minFrames = MIN_FRAMES_MEDIAN;
    break;
  }


  /*
   * Get the sky level determination method.
   */

  skyMethodTag = (char *)pilDfsDbGetString(parameter, "SkyMethod");

  if ((skyMethodEntry =
        strselect(skyMethodTag, skyMethodNames, nSkyMethods)) < 0) {
    cpl_msg_error(task,
                "%s: Invalid sky level determination method.", skyMethodTag);
    return EXIT_FAILURE;
  }

  skyMethod = skyMethods[skyMethodEntry];

  if (skyMethod == SKY_FIT) {
    polyDeg = pilDfsDbGetInt(parameter, "PolyOrder", 0);
    if (polyDeg < 0) {
      cpl_msg_error(task, "Invalid choice of polynomial for sky level "
                  "modeling: degree should not be less than zero");
      return EXIT_FAILURE;
    }

    skySigmaLow = pilDfsDbGetDouble(parameter, "SkyKSigmaLow", 3.0);
    skySigmaHigh = pilDfsDbGetDouble(parameter, "SkyKSigmaHigh", 3.0);

    if (skySigmaLow < 0.001 || skySigmaHigh < 0.001) {
      cpl_msg_error(task, "Invalid choice of K-sigma rejection for sky level "
                  "modeling: number of sigmas should be positive.");
      return EXIT_FAILURE;
    }

  }


  /*
   * Determine wavelength sampling mode in spectral extraction
   */

#ifdef ONLINE_MODE

  samplingMethodTag = (char *)pilDfsDbGetString(parameter, "SamplingMethod");

  if ((samplingMethodEntry = strselect(samplingMethodTag,
                             samplingMethodNames, nSamplingMethods)) < 0) {
    cpl_msg_error(task,
                "%s: Invalid wavelength sampling mode", samplingMethodTag);
    return EXIT_FAILURE;
  }

  samplingMethod = samplingMethods[samplingMethodEntry];

#endif


  /*
   * Get "fuzz" parameter (allowing extra pixels around spectra to
   * extract, when searching for its edges).
   */

  fuzz = pilDfsDbGetInt(parameter, "Fuzz", 5);

  if (fuzz < 0) {
    cpl_msg_error(task, "Fuzz parameter must be positive");
    return EXIT_FAILURE;
  }


  /*
   * Get number of pixels at slit ends to be excluded from object search
   * and from sky level determination.
   */

  slitMargin = pilDfsDbGetInt(parameter, "SlitMargin", 2);

  if (fuzz < 0) {
    cpl_msg_error(task, "Fuzz parameter must be positive");
    return EXIT_FAILURE;
  }


  /*
   * Get tolerance for drift of slit on the CCD.
   */

  slitTolerance = pilDfsDbGetDouble(parameter, "SlitTolerance", 0.0);

  if (slitTolerance < 0) {
    cpl_msg_error(task, "SlitTolerance parameter must be positive");
    return EXIT_FAILURE;
  }

  /*
   * To ensure zero-tolerance (i.e., all slits are classified as "long",
   * no matter what), we set internalli a negative tolerance. 
   */

  if (slitTolerance < 0.00001)
    slitTolerance = -1.0;


  /*
   * Get size of spectrum extracted around any skyline.
   */

  lineWidth = pilDfsDbGetInt(parameter, "LineWidth", 16);

  if (lineWidth < 0) {
    cpl_msg_error(task, "LineWidth parameter must be positive");
    return EXIT_FAILURE;
  }


  /*
   * Get object detection level in units of sigma.
   */

  detLevel = pilDfsDbGetDouble(parameter, "DetectionLevel", 2.0);

  if (detLevel < 0.0) {
    cpl_msg_error(task, "DetectionLevel parameter must be positive");
    return EXIT_FAILURE;
  }


  /*
   * Get number of levels to use in watershed method.
   */

  numLevels = pilDfsDbGetInt(parameter, "WatershedLevels", 32);

  if (numLevels < 2) {
    cpl_msg_error(task, "WatershedLevels parameter must be at least 2");
    return EXIT_FAILURE;
  }


  /*
   * Get flux fraction to use in watershed.
   */

  objFrac = pilDfsDbGetDouble(parameter, "WatershedFraction", 0.01);

  if (objFrac < 0.0) {
    cpl_msg_error(task, "WatershedFraction parameter must be positive");
    return EXIT_FAILURE;
  }

  /*
   * Get spectrum fraction to be collapsed
   */

  specFrac = pilDfsDbGetDouble(parameter, "SpectrumFraction", 0.8);

  if (specFrac < 0 || specFrac >1) {
    cpl_msg_error(task, "SpectrumFraction parameter must be >0 and <1");
    return EXIT_FAILURE;
  }


  /*
   * Get minimal size for an object candidate to be considered an object.
   */

  minObjectSize = pilDfsDbGetInt(parameter, "MinObjectSize", 2);

  if (minObjectSize < 1) {
    cpl_msg_error(task, "MinObjectSize parameter must be positive");
    return EXIT_FAILURE;
  }


  /*
   * Get max size for an object candidate for NOT trying to deblend it
   * into sub-objects.
   */

  minCompositeSize = pilDfsDbGetInt(parameter, "MaxObjectSize", 7);

  if (minCompositeSize < minObjectSize) {
    cpl_msg_error(task,
                "MaxObjectSize parameter must be greater than MinObjectSize");
    return EXIT_FAILURE;
  }


  /*
   * Check if the bad pixels should be corrected.
   */

  cleanBadPixel = pilDfsDbGetBool(parameter, "CleanBadPixel", 0);


  /*
   * Check if the spectro-photometric calibration should be applied.
   */

  calibrateFlux = pilDfsDbGetBool(parameter, "CalibrateFlux", 0);


  /*
   * Check if the fringing should be corrected, and with what method.
   */

  fringingCorr = pilDfsDbGetBool(parameter, "FringingCorr", 0);

  if (fringingCorr) {
    fringMethodTag = (char *)pilDfsDbGetString(parameter, "FringingMethod");
    if (strcmp(fringMethodTag, "Raw") == 0)
      fringMethod = 0;
    if (strcmp(fringMethodTag, "Resampled") == 0)
      fringMethod = 1;

    minFringOffset = pilDfsDbGetDouble(parameter, "FringingOffset", 4.0);
  }


  /*
   * Check what kind of 1D extraction should be applied
   */

  horneExtraction = pilDfsDbGetBool(parameter, "HorneExtraction", 1);


  /*
   * Optional modeling of wavelength solution within each slit
   */

  modelSlit = pilDfsDbGetBool(parameter, "ModelSlit", 0);


  /*
   * Order of polynomial for modeling of wavelength solution within each slit.
   */

  order = pilDfsDbGetInt(parameter, "ModelSlitOrder", 0);


  /*
   * Check if the wavelength calibration should be refined on skylines.
   */

  useSkylines = pilDfsDbGetBool(parameter, "UseSkylines", 1);


  /*
   * Check if the intermediate data reduction products should be saved
   * to disk.
   */

  saveIntermediate = pilDfsDbGetBool(parameter, "SaveIntermediate", 0);


  /*
   * Load input frames
   */

  cpl_msg_info(task, "Loading input frames...");

  sphotFrame = pilSofLookup(sof, sphotTableCategory);
  if (sphotFrame)
    pilFrmSetType(sphotFrame, PIL_FRAME_TYPE_CALIB);

  atmFrame = pilSofLookup(sof, atmTableCategory);
  if (atmFrame)
    pilFrmSetType(atmFrame, PIL_FRAME_TYPE_CALIB);

  if (calibrateFlux) {

    /*
     * Check that a SPH table is present in input frames.
     */

    if (!sphotFrame) {
      cpl_msg_error(task, "Missing SPH Table: "
                    "spectrophotometric calibration cannot be applied.");
      return EXIT_FAILURE;
    }


    /*
     * Check that an atmospheric extinction table is present in input frames.
     */

    if (!atmFrame) {
      cpl_msg_error(task, "Missing atmospheric extinction table: "
                    "spectrophotometric calibration cannot be applied.");
      return EXIT_FAILURE;
    }
  }


  /*
   * Get raw MOS frames
   */

  mosCount = pilSofFrameCount(sof, mosCategory);
  if (mosCount < minFrames) {
    cpl_msg_error(task, "At least %d raw MOS exposures should be input for "
                "specified combine method '%s'", minFrames, combMethodTag);
    return EXIT_FAILURE;
  }

  if ((mosList = (VimosImage **)cpl_calloc(mosCount, sizeof(VimosImage *)))) {

    mosFrame = pilSofLookupNext(sof, mosCategory);

    for (i = 0; i < mosCount; i++) {
      if ((mosList[i] = openOldFitsFile(pilFrmGetName(mosFrame), 1, 0))) {
        pilFrmSetType(mosFrame, PIL_FRAME_TYPE_RAW);
        closeFitsImage(mosList[i], 0);
      }
      else {
        cpl_msg_error(task, "Failure opening MOS science image %d", i + 1);
        for (j = 0; j < i; j++)
          deleteImage(mosList[j]);
        cpl_free(mosList);
        return EXIT_FAILURE;
      }
      mosFrame = pilSofLookupNext(sof, NULL);
    }
  }
  else {
    cpl_msg_error(task, "Failure creating list of input MOS exposures");
    return EXIT_FAILURE;
  }


  if (fringingCorr) {

    /*
     * Find jittering offsets, to see if they are compatible with a
     * fringing correction.
     */

    readDoubleDescriptor(mosList[0]->descs, pilTrnGetKeyword("PixelScale"),
                         &pixScale, NULL);
    readDoubleDescriptor(mosList[0]->descs, pilTrnGetKeyword("Alpha"),
                         &refAlpha, NULL);
    readDoubleDescriptor(mosList[0]->descs, pilTrnGetKeyword("Delta"),
                         &refDelta, NULL);

    cosDelta = fabs(cos(PI_NUMB * refDelta / 180.));

    aOffset = malloc(mosCount * sizeof(float));
    dOffset = malloc(mosCount * sizeof(float));

    aOffset[0] = 0.0;
    dOffset[0] = 0.0;

    for (i = 1; i < mosCount; i++) {

      if (cosDelta > 0.0001) {
        readDoubleDescriptor(mosList[i]->descs, pilTrnGetKeyword("Alpha"),
                             &dValue, NULL);
        aOffset[i] = 3600 * ((refAlpha - dValue) / cosDelta) / pixScale;
      }
      else {
        aOffset[i] = 0.0;
      }

      readDoubleDescriptor(mosList[i]->descs, pilTrnGetKeyword("Delta"),
                           &dValue, NULL);
      dOffset[i] = 3600 * (refDelta - dValue) / pixScale;

    }

    sort(mosCount, aOffset);
    sort(mosCount, dOffset);

    if (aOffset[mosCount-1] - aOffset[0] > dOffset[mosCount-1] - dOffset[0]) {
      if (cosDelta < 0.0001)
        cpl_msg_warning(task, "Close to the pole - unsafe offsets in RA!");
      offset = aOffset;
      cpl_free(dOffset);
    }
    else {
      offset = dOffset;
      cpl_free(aOffset);
    }

    for (i = 1; i < mosCount; i++)
      offset[i - 1] = offset[i] - offset[i - 1];

    minOffset = offset[0];
    for (i = 1; i < mosCount - 1; i++)
      if (minOffset > offset[i])
        minOffset = offset[i];

    cpl_free(offset);

    if (minOffset < minFringOffset) {
      cpl_msg_warning(task, "Minimum observed offset between frames is %.2f "
                    "pixel, while at least %.2f is required for fringing "
                    "correction. The sky fringing correction is now DISABLED: "
                    "run the recipe again with a lower threshold if you "
                    "want to try the fringing correction anyway.", 
                    minOffset, minFringOffset);
      fringingCorr = 0;
    }
    else {
      cpl_msg_info(task, 
                 "Minimum observed offset between frames is %.2f pixels.",
                 minOffset);
    }

  }


  /*
   * Exposure times
   */

  sumExpTime = 0.0;
  for (i = 0; i < mosCount; i++) {
    readDoubleDescriptor(mosList[i]->descs, pilTrnGetKeyword("ExposureTime"),
                         &expTime, NULL);
    sumExpTime += expTime;
  }


  /*
   * Get the master bias frame
   */

  error = 1;

  if ((biasFrame = pilSofLookup(sof, biasCategory))) {
    pilFrmSetType(biasFrame, PIL_FRAME_TYPE_CALIB);
    if ((biasImage = openOldFitsFile(pilFrmGetName(biasFrame), 1, 0))) {
      closeFitsImage(biasImage, 0);
      error = 0;
    }
    else
      cpl_msg_error(task, "Failure opening master bias frame");
  }
  else
    cpl_msg_error(task, "No master bias found in input");

  if (error) {
    for (i = 0; i < mosCount; i++)
      deleteImage(mosList[i]);
    cpl_free(mosList);
    return EXIT_FAILURE;
  }


  /*
   * Recreate bias overscans using as a reference the first MOS frame
   * in the list.
   */

  if ((tmpImage = growOverscans(biasImage, mosList[0]))) {
    if (biasImage != tmpImage) {
      deleteImage(biasImage);
      biasImage = tmpImage;
    }
  }
  else {
    cpl_msg_error(task, "Failure in growing overscans in master bias");
    for (i = 0; i < mosCount; i++)
      deleteImage(mosList[i]);
    cpl_free(mosList);
    deleteImage(biasImage);
    return EXIT_FAILURE;
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

  cpl_msg_info(task, "Bias removal...");

  for (i = 0; i < mosCount; i++) {
    if (VmSubBias(mosList[i], biasImage, biasMethod) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot remove bias from raw MOS image %d", i + 1);
      for (i = 0; i < mosCount; i++)
        deleteImage(mosList[i]);
      cpl_free(mosList);
      deleteImage(biasImage);
      return EXIT_FAILURE;
    }
  }

  deleteImage(biasImage);

  if (saveIntermediate) {
    for (i = 0; i < mosCount; i++) {
      sprintf(tmpName, "after_bias_subtraction%d.fits", i + 1);
      createFitsImage(tmpName, mosList[i], "test");
    }
  }


  /*
   * Get the (optional) master dark frame
   */

  if ((darkFrame = pilSofLookup(sof, darkCategory))) {
    pilFrmSetType(darkFrame, PIL_FRAME_TYPE_CALIB);
    if ((darkImage = openOldFitsFile(pilFrmGetName(darkFrame), 1, 0))) {
      closeFitsImage(darkImage, 0);
    }
    else {
      cpl_msg_error(task, "Failure opening master dark frame");
      for (i = 0; i < mosCount; i++)
        deleteImage(mosList[i]);
      cpl_free(mosList);
      return EXIT_FAILURE;
    }
  }
  else
    cpl_msg_warning(task, "No master dark in input, dark subtraction "
                  "will not be performed");


  /*
   * Now go for the (optional) dark subtraction
   */

  if (darkImage) {
    cpl_msg_info(task, "Dark subtraction...");

    for (i = 0; i < mosCount; i++) {
      if (VmSubDark(mosList[i], darkImage) == EXIT_FAILURE) {
        cpl_msg_error(task, "Cannot subtract dark from MOS image %d", i + 1);
        for (i = 0; i < mosCount; i++)
          deleteImage(mosList[i]);
        cpl_free(mosList);
        deleteImage(darkImage);
        return EXIT_FAILURE;
      }
    }
    deleteImage(darkImage);
  }


  /*
   * Get the (optional) master flat field frame.
   */

  if ((flatFrame = pilSofLookup(sof, flatCategory))) {
    pilFrmSetType(flatFrame, PIL_FRAME_TYPE_CALIB);
    if ((flatImage = openOldFitsFile(pilFrmGetName(flatFrame), 1, 0)))
      closeFitsImage(flatImage, 0);
    else {
      cpl_msg_error(task, "Failure opening master MOS flat field frame");
      for (i = 0; i < mosCount; i++)
        deleteImage(mosList[i]);
      cpl_free(mosList);
      return EXIT_FAILURE;
    }
  }


  /*
   * The grism table enters the show:
   */

  error = 1;

  if ((grismFrame = pilSofLookup(sof, grismTableCategory))) {
    pilFrmSetType(grismFrame, PIL_FRAME_TYPE_CALIB);
    if ((grismFile = openOldFitsFile(pilFrmGetName(grismFrame), 0, 0))) {
      if ((grismTable = newGrismTable())) {
        if (readFitsGrismTable(grismTable, grismFile->fptr) == VM_TRUE) {
          closeFitsImage(grismFile, 0);
          error = 0;
        }
        else {
          cpl_msg_error(task, "Failure reading grism table");
          deleteTable(grismTable);
        }
      }
      else
        cpl_msg_error(task, "Not enough memory");

      deleteImage(grismFile);

    }
    else
      cpl_msg_error(task, "Failure in opening grism table");
  }
  else
    cpl_msg_error(task, "No input grism table found");

  if (error) {
    deleteImage(flatImage);
    for (i = 0; i < mosCount; i++)
      deleteImage(mosList[i]);
    cpl_free(mosList);
    return EXIT_FAILURE;
  }


  /*
   * Setup the list of (initially) identical extraction tables.
   */

  if ((extTablesList = (VimosExtractionTable **)cpl_calloc(mosCount, 
                        sizeof(VimosExtractionTable *)))) {

    if ((extrFrame = pilSofLookup(sof, extrTableCategory))) {

      pilFrmSetType(extrFrame, PIL_FRAME_TYPE_CALIB);

      /*
       * There is an extraction table in input. Read it in as many
       * times as many input science frames we have. Note that
       * defining a copy contructor would probably make this code
       * more efficient.
       */

      if ((extrFile = openOldFitsFile(pilFrmGetName(extrFrame), 0, 0))) {

        for (i = 0; i < mosCount; i++) {
          if ((extTablesList[i] = newExtractionTable())) {
            if (readFitsExtractionTable(extTablesList[i], extrFile->fptr)
                                                               == VM_TRUE) {
             /*
              * Replace in extraction table the fit of the wavelength
              * solution coefficients trend along each slit.
              */

              if (modelSlit)
                modelWavcal(extTablesList[i], order);

            }
            else {
              cpl_msg_error(task, 
                          "Failure reading the extraction table (%d)", i);
              i++;
              error = 1;
              break;
            }
          }
          else {
            cpl_msg_error(task, "Cannot create new extraction table (%d)", i);
            error = 1;
            break;
          }
        }

        if (!error) {
          if ((combExtractionTable = newExtractionTable())) {
            if (readFitsExtractionTable(combExtractionTable, extrFile->fptr)
                                                              == VM_FALSE) {
              cpl_msg_error(task, "Failure reading the extraction table");
              error = 1;
            }
          }
          else {
            cpl_msg_error(task, "Cannot create new extraction table");
            error = 1;
          }
        }

        closeFitsImage(extrFile, 0);
        deleteImage(extrFile);

        if (error) {
          for (j = 0; j < i; j++)
            deleteExtractionTable(extTablesList[j]);
          cpl_free(extTablesList);
          deleteExtractionTable(combExtractionTable);
        }

      }
      else {
        cpl_msg_error(task, "Failure opening the extraction table");
        error = 1;
      }
    }
    else {

      /*
       * If no extraction table was found in the input Set-Of-Frames,
       * derive the extraction table from global models. The global
       * distorsion models coefficients are read from the science frames.
       * One extraction table is created for each exposure: they are
       * initially all identical, and then they differentiate at the
       * CalShift correction.
       */

      for (i = 0; i < mosCount; i++) {
        if (!(extTablesList[i] = VmSpExTab(mosList[0], 
                                           grismTable, NULL, NULL))) { 
          cpl_msg_error(task, "Cannot compute extraction table %d", i + 1);
          error = 1;
          break;
        }
      }

      if (!error) {
        if (!(combExtractionTable = VmSpExTab(mosList[0], 
                                              grismTable, NULL, NULL))) { 
          cpl_msg_error(task, "Cannot compute the combined extraction table.");
          error = 1;
        }
      }

      if (error) {
        for (j = 0; j < i; j++)
          deleteExtractionTable(extTablesList[j]);
        cpl_free(extTablesList);
        deleteExtractionTable(combExtractionTable);
      }

    }
  }
  else {
    cpl_msg_error(task, "Failure creating extraction tables list");
    error = 1;
  }

  if (error) {
    deleteTable(grismTable);
    deleteImage(flatImage);
    for (j = 0; j < mosCount; j++)
      deleteImage(mosList[j]);
    cpl_free(mosList);
    return EXIT_FAILURE;
  }


  /*
   * Apply flat field (if present).
   */

  if (flatImage) {

    for (i = 0; i < mosCount; i++) {
      if ((tmpImage = VmSpApplyFF(mosList[i], flatImage, extTablesList[i]))) {
        deleteImage(mosList[i]);
        mosList[i] = tmpImage;
      }
      else {
        cpl_msg_error(task, "Failure flat-fielding MOS exposure %d", i + 1);
        deleteImage(flatImage);
        deleteTable(grismTable);
        for (i = 0; i < mosCount; i++)
          deleteImage(mosList[i]);
        cpl_free(mosList);
        for (i = 0; i < mosCount; i++)
          deleteExtractionTable(extTablesList[i]);
        cpl_free(extTablesList);
        return EXIT_FAILURE;
      }
    }
  
    deleteImage(flatImage);

    if (saveIntermediate) {
      for (i = 0; i < mosCount; i++) {
        sprintf(tmpName, "after_flat_fielding%d.fits", i + 1);
        createFitsImage(tmpName, mosList[i], "test");
      }
    }

  } 

  fringFrame = pilSofLookup(sof, mosFringesSkyCategory);
  if (!fringFrame)
    fringFrame = pilSofLookup(sof, mosFringesCategory);
  if (fringFrame)
    pilFrmSetType(fringFrame, PIL_FRAME_TYPE_CALIB);

  if (fringingCorr) {

    /*
     * Get an (optional) image of fringes + sky, product of previous runs
     * with this recipe. If this image is found, then the fringing correction 
     * will not be computed, and this image of fringes + sky will be used 
     * instead. The category mosFringesSkyCategory is found if the image 
     * was the product of the fringing correction "Raw", while the category 
     * mosFringesCategory is found if the image was the product of the
     * fringing correction "Resampled". The chosen fringing correction
     * method must correspond to the correct image category.
     */

    friCount = pilSofFrameCount(sof, mosFringesSkyCategory)
             + pilSofFrameCount(sof, mosFringesCategory);

    if (friCount > 1) {
      cpl_msg_error(task, 
                  "Just one fringes image may be given in input, %d found",
                  friCount);
      error = 1;
    }
    else {
      if ((fringFrame = pilSofLookup(sof, mosFringesSkyCategory))) {
        if (fringMethod != 0) {
          cpl_msg_error(task, "The input fringes image, %s, is compatible with "
                      "the fringing correction method 'Raw'. Please select "
                      "the 'Raw' method.", pilFrmGetName(fringFrame));
          error = 1;
        }
      }
      else if ((fringFrame = pilSofLookup(sof, mosFringesCategory))) {
        if (fringMethod != 1) {
          cpl_msg_error(task, "The input fringes image, %s, is compatible with "
                      "the fringing correction method 'Resampled'. Please "
                      "select the 'Resampled' method.", 
                      pilFrmGetName(fringFrame));
          error = 1;
        }
      }
      
      if (fringFrame != NULL && error == 0) {
        if ((tmpImage = openOldFitsFile(pilFrmGetName(fringFrame), 1, 0))) {
          closeFitsImage(tmpImage, 0);

          if (fringMethod == 0)
            imaFringes = tmpImage;
          else
            imaFringes2D = tmpImage;

        }
        else {
          cpl_msg_error(task, "Failure opening fringes frame");
          error = 1;
        }
      }
    }

    if (error) {
      deleteTable(grismTable);
      for (i = 0; i < mosCount; i++)
        deleteImage(mosList[i]);
      cpl_free(mosList);
      for (i = 0; i < mosCount; i++)
        deleteExtractionTable(extTablesList[i]);
      cpl_free(extTablesList);
      return EXIT_FAILURE;
    }

  }


  /*
   * The "Raw" fringing correction consists on subtracting the median
   * image from all input images. Flat fielding should be avoided in this
   * case, because the flat field contains fringes too. Alternatively, 
   * a flat field normalized in such a way that the fringes were modeled 
   * and removed together with the continuum should be provided: this 
   * might be obtained by applying a very small smoothing-box in recipe 
   * vmspflat.
   */

  if (fringingCorr == 1 && fringMethod == 0) {
    cpl_msg_info(task, "Applying fringing correction method 'Raw':");
    if (imaFringes)
      cpl_msg_info(task, "Using input image of fringes + sky...");
    else {
      cpl_msg_info(task, "Creating own image of fringes + sky...");

      if (mosCount == 2)
        imaFringes = frCombMinMaxReject(mosList, 0, 1, 2);
      else
        imaFringes = frCombMedian(mosList, mosCount, 0);

      copyAllDescriptors(mosList[0]->descs, &(imaFringes)->descs);
    }
/*
 * FIXED - This part is commented out: this code has been moved after
 * the alignment to the sky lines. The sky should not be subtracted 
 * before  that!
 * 
 *  for (i = 0; i < mosCount; i++)
 *    imageArithLocal(mosList[i], imaFringes, VM_OPER_SUB);
 *
 *  if (saveIntermediate) {
 *    for (i = 0; i < mosCount; i++) {
 *      sprintf(tmpName, "after_raw_fringing_and_sky_removal%d.fits", i + 1);
 *      createFitsImage(tmpName, mosList[i], "test");
 *    }
 *  }
 */
  }


  /*
   * A bad pixel table is required in the input set of frames for
   * image cleaning.
   */

  ccdFrame = pilSofLookup(sof, ccdTableCategory);
  if (ccdFrame)
    pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

  if (cleanBadPixel) {

    error = 1;

    if (ccdFrame) {
      cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));
      if ((ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
        closeFitsTable(ccdTable, 0);
        error = 0;
      }
      else
        cpl_msg_error(task, "Failure in opening CCD table");
    }
    else
      cpl_msg_error(task, "No CCD table in input: cannot clean bad pixels.");

    if (error) {
      deleteImage(imaFringes);
      deleteTable(grismTable);
      for (i = 0; i < mosCount; i++)
        deleteImage(mosList[i]);
      cpl_free(mosList);
      for (i = 0; i < mosCount; i++)
        deleteExtractionTable(extTablesList[i]);
      cpl_free(extTablesList);
      return EXIT_FAILURE;
    }


    /*
     *  Image cleaning
     */

    cpl_msg_info(task, "Cleaning bad pixels in MOS science exposures...");

    for (i = 0; i < mosCount; i++) {
      if (EXIT_FAILURE == cleanBadPixels(mosList[i], ccdTable, 0)) {
        cpl_msg_error(task, "Cannot clean MOS science exposure %d", i + 1);
        deleteImage(imaFringes);
        deleteTable(grismTable);
        deleteTable(ccdTable);
        for (i = 0; i < mosCount; i++)
          deleteImage(mosList[i]);
        cpl_free(mosList);
        for (i = 0; i < mosCount; i++)
          deleteExtractionTable(extTablesList[i]);
        cpl_free(extTablesList);
        return EXIT_FAILURE;
      }
    }

    deleteTable(ccdTable);

    if (saveIntermediate) {
      for (i = 0; i < mosCount; i++) {
        sprintf(tmpName, "after_bad_pixel_cleaning%d.fits", i + 1);
        createFitsImage(tmpName, mosList[i], "test");
      }
    }
  }


  if (useSkylines) {

    /*
     * Shift calibration: here the extraction tables get different
     * from each other.
     */
  
    cpl_msg_info(task, "Calibrating shifts in MOS science exposures...");
  
    for (i = 0; i < mosCount; i++) {
      if (VmSpCalShifts(mosList[i], grismTable, extTablesList[i], 
                        0, lineWidth, fuzz) == EXIT_FAILURE) {
        cpl_msg_error(task, "Failure correcting offset in MOS image %d", i + 1);
        for (i = 0; i < mosCount; i++)
          deleteImage(mosList[i]);
        cpl_free(mosList);
        for (i = 0; i < mosCount; i++)
          deleteExtractionTable(extTablesList[i]);
        cpl_free(extTablesList);
        deleteTable(grismTable);
        deleteImage(imaFringes);
        return EXIT_FAILURE;
      }
    }
  }

  deleteTable(grismTable);

  if (fringingCorr == 1 && fringMethod == 0) {
    for (i = 0; i < mosCount; i++)
      imageArithLocal(mosList[i], imaFringes, VM_OPER_SUB);

    if (saveIntermediate) {
      for (i = 0; i < mosCount; i++) {
        sprintf(tmpName, "after_raw_fringing_and_sky_removal%d.fits", i + 1);
        createFitsImage(tmpName, mosList[i], "test");
      }
    }
  }


  /*
   * Object detection in raw MOS exposures. Note that this is much more
   * efficient if the fringing correction method "Raw" is applied, because
   * at this stage the sky fringes are already removed.
   */

  cpl_msg_info(task, "Detecting objects in 2D images...");

  if ((winTablesList = (VimosWindowTable **)cpl_calloc(mosCount, 
                        sizeof(VimosWindowTable *)))) {

    for (i = 0; i < mosCount; i++) {

      if (!(winTablesList[i] = VmSpDetObj(mosList[i], extTablesList[i], 
                                          numLevels, slitMargin, detLevel, 
                                          objFrac, limFrac, minObjectSize, 
                                          minCompositeSize, slitTolerance,
                                          specFrac))) {
        cpl_msg_error(task, "Failure deriving window table for MOS exposure %d",
                      i + 1);
        for (j = 0; j < i; j++)
          deleteWindowTable(winTablesList[j]);
        cpl_free(winTablesList);
        error = 1;
        break;
      }

    }
  }
  else {
    cpl_msg_error(task, "Failure creating window tables list");
    error = 1;
  }

  if (error) {
    for (j = 0; j < mosCount; j++)
      deleteImage(mosList[j]);
    cpl_free(mosList);
    for (j = 0; j < mosCount; j++)
      deleteExtractionTable(extTablesList[j]);
    cpl_free(extTablesList);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  /*
   *  Sky subtraction for MOS exposures. If the "Raw" fringing correction
   *  was applied, this step removes the residuals of the sky+fringes 
   *  subtraction made at the beginning of the recipe. Such residuals
   *  would be the result of sky variability during the jittered exposure.
   */

  if (fringingCorr == 1 && fringMethod == 0)
     cpl_msg_info(task, 
     "Removing the residuals of sky+fringes subtraction from MOS exposures...");
  else
     cpl_msg_info(task, "Sky subtraction from MOS exposures...");

  error = 1;

  if ((mosList2D = (VimosImage **)cpl_calloc(mosCount, sizeof(VimosImage *)))) {
    if ((skyList2D = (VimosImage **)cpl_calloc(mosCount, sizeof(VimosImage *))))
      error = 0;
    else
      cpl_msg_error(task, "Failure creating list of 2D sky images");
  }
  else
    cpl_msg_error(task, "Failure creating list of 2D extracted images");

  if (error) {
    cpl_free(mosList2D);
    cpl_free(skyList2D);
    for (j = 0; j < mosCount; j++)
      deleteWindowTable(winTablesList[j]);
    cpl_free(winTablesList);
    for (j = 0; j < mosCount; j++)
      deleteImage(mosList[j]);
    cpl_free(mosList);
    for (j = 0; j < mosCount; j++)
      deleteExtractionTable(extTablesList[j]);
    cpl_free(extTablesList);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }


  for (i = 0; i < mosCount; i++) {

    /*  
     *  Sky-sub in short slits 
     */

    if (!(outSpSkyFra = VmSpSkyFra(mosList[i], extTablesList[i], 
                                   winTablesList[i], skyMethod, 
                                   polyDeg, slitMargin, 10, 3, 
                                   skySigmaLow, skySigmaHigh))) {

      cpl_msg_error(task, "Failure in sky subtraction in short slits"
                  "for MOS science exposure %d", i + 1);
      for (j = 0; j < mosCount; j++)
        deleteImage(mosList[j]);
      cpl_free(mosList);
      for (j = 0; j < mosCount; j++)
        deleteExtractionTable(extTablesList[j]);
      cpl_free(extTablesList);
      deleteImage(imaFringes);
      return EXIT_FAILURE;

    }
    

    if (i == 0 && fringingCorr == 1 && fringMethod == 0) {

      /*
       * If the "Raw" fringing correction was applied, we need to sum
       * the sky+fringes map to the residual map. In order to do this, 
       * we must rectify also the sky+fringes map.
       * This is a quick trick for rectifying also the sky+fringes map.
       */

      aritempo[0] = outSpSkyFra[0];
      aritempo[1] = imaFringes;

      tempo = VmSpEx2D(aritempo, extTablesList[i], samplingMethod);
      
      imaFringes2D = tempo[1];
      deleteImage(tempo[0]);
      cpl_free(tempo);

    }


    /* 
     *  2D extraction 
     */

    outSpEx2D = VmSpEx2D(outSpSkyFra, extTablesList[i], samplingMethod);

    deleteImage(outSpSkyFra[0]);
    deleteImage(outSpSkyFra[1]);
    cpl_free(outSpSkyFra);
    
    if (!outSpEx2D) {
      cpl_msg_error(task, "Failure in 2D extraction of spectra "
                  "from MOS science exposure %d", i + 1);
      for (j = 0; j < mosCount; j++)
        deleteImage(mosList[j]);
      cpl_free(mosList);
      for (j = 0; j < mosCount; j++)
        deleteExtractionTable(extTablesList[j]);
      cpl_free(extTablesList);
      for (j = 0; j < mosCount; j++)
        deleteWindowTable(winTablesList[j]);
      cpl_free(winTablesList);
      deleteImage(imaFringes2D);
      deleteImage(imaFringes);
      return EXIT_FAILURE;
    }

    if (saveIntermediate) {
      sprintf(tmpName, "after_2D_extraction%d.fits", i + 1);
      createFitsImage(tmpName, outSpEx2D[0], "test");
    }

    
    /*  
     *  Sky-sub in long slits 
     */

    if (!(outSpSkyEx = VmSpSkyExt(outSpEx2D, winTablesList[i], 
                                  skyMethod, polyDeg, slitMargin, 
                                  10, 3, skySigmaLow, skySigmaHigh))) {
      cpl_msg_error(task, "Failure in sky subtraction in long slits");
      deleteImage(outSpEx2D[0]);
      deleteImage(outSpEx2D[1]);
      cpl_free(outSpEx2D);
      for (j = 0; j < mosCount; j++)
        deleteImage(mosList[j]);
      cpl_free(mosList);
      for (j = 0; j < mosCount; j++)
        deleteExtractionTable(extTablesList[j]);
      cpl_free(extTablesList);
      for (j = 0; j < mosCount; j++)
        deleteWindowTable(winTablesList[j]);
      cpl_free(winTablesList);
      deleteImage(imaFringes2D);
      deleteImage(imaFringes);
      return EXIT_FAILURE;
    }

    if (saveIntermediate) { 
      sprintf(tmpName, "after_2D_sky_subtraction%d.fits", i + 1);
      createFitsImage(tmpName, outSpSkyEx[0], "test");
      sprintf(tmpName, "after_2D_sky_subtraction_SKY%d.fits", i + 1);
      createFitsImage(tmpName, outSpSkyEx[1], "test");
    }

    deleteImage(outSpEx2D[0]);
    deleteImage(outSpEx2D[1]);
    cpl_free(outSpEx2D);
    deleteImage(mosList[i]);
    mosList2D[i] = outSpSkyEx[0];
    skyList2D[i] = outSpSkyEx[1];
  } 


  /*
   * In the "Resampled" fringing correction method, the fringes map is
   * obtained _after_ the sky was modeled and removed. Of course this is
   * risky, because it is based on the assumption that the sky modeling
   * method would be capable of telling sky from fringes. However, it
   * is sometimes preferable to give up the construction of a true
   * fringes map: fringes are fitted with the sky, and the fringes
   * map would just contain the residuals of such fit.
   */

  if (fringingCorr == 1) {
    if (fringMethod == 1) {
      cpl_msg_info(task, "Applying fringing correction 'Resampled':");
      if (imaFringes2D)
        cpl_msg_info(task, "Using input image of fringes...");
      else {
        cpl_msg_info(task, "Creating own image of fringes...");

        if (mosCount == 2)
          imaFringes2D = frCombMinMaxReject(mosList2D, 0, 1, 2);
        else
          imaFringes2D = frCombMedian(mosList2D, mosCount, 0);

        copyAllDescriptors(mosList2D[0]->descs, &(imaFringes2D)->descs);
      }

      for (i = 0; i < mosCount; i++)
        imageArithLocal(mosList2D[i], imaFringes2D, VM_OPER_SUB);

      if (saveIntermediate) {
        for (i = 0; i < mosCount; i++) {
          sprintf(tmpName, "after_resampled_fringing_removal%d.fits", i + 1);
          createFitsImage(tmpName, mosList2D[i], "test");
        }
      }
    }

    /*
     * Note that the subtracted fringes must be added to the sky,
     * no matter what kind of fringing correction was applied.
     */

    for (i = 0; i < mosCount; i++)
      imageArithLocal(skyList2D[i], imaFringes2D, VM_OPER_ADD);

  }


  /*
   *  Fringing corrections (Consortium). It might be enabled for comparing
   *  the results.
   */

/***
  if (fringingCorr) {
    cpl_msg_info(task, "Applying fringing corrections 3...");

    if (VmSpFringCorr(mosList2D, winTablesList, mosCount, 0, 0) 
                                                      == EXIT_FAILURE) {
      cpl_msg_error(task, "Failure applying fringing corrections");
      for (i = 0; i < mosCount; i++) 
        deleteExtractionTable(extTablesList[i]);
      cpl_free(extTablesList);
      for (j = 0; j < mosCount; j++)
        deleteImage(mosList2D[j]);
      cpl_free(mosList2D);
      for (j = 0; j < mosCount; j++)
        deleteImage(skyList2D[j]);
      cpl_free(skyList2D);
      for (j = 0; j < mosCount; j++)
        deleteWindowTable(winTablesList[j]);
      cpl_free(winTablesList);
      deleteImage(imaFringes2D);
      deleteImage(imaFringes);
      return EXIT_FAILURE;
    }
  }
***/

      
  /*
   *  Stack 2D images. The stacking is made keeping into account the
   *  jittering, of course. Note that also the sky+fringes images are
   *  stacked together in the same way, in order to have an image of 
   *  the total sky+fringes. This of course will destroy the neat and
   *  clean look of the fringing map, but it is a faithful representation
   *  of all that has been subtracted from the stacked data. Note that
   *  this is necessary for the estimation of the variance of the signal
   *  in the optimal extraction algorithm.
   */

  offSets = newFloatArray(mosCount);

  error = 0;
  if (!(outSpSkyExStack = (VimosImage **)cpl_calloc(2, sizeof(VimosImage *)))) {
    cpl_msg_error(task, "Failure in creating 2D images array");
    error = 1;
  } 
  else {
    cpl_msg_info(task, "Stacking 2D images...");
    if (!(outSpSkyExStack[0] = VmSpStack2D(mosList2D, winTablesList, 
                                           extTablesList[0], NULL, 
                                           mosCount, combMethod, 
                                           &combParameter, offSets, 
                                           1, 0, 0))) {
      cpl_msg_error(task, "Failure in creating combined 2D image");
      error = 1;
    } 
    else {
      if (!(twoDMap = newIntArray(outSpSkyExStack[0]->ylen))) {
        cpl_msg_error(task, "Failure in creating 2D image map");
        error = 1;
      } 
      else {
        cpl_msg_info(task, "Stacking 2D sky images...");
        if (!(outSpSkyExStack[1] = VmSpStack2D(skyList2D, winTablesList, 
                                               extTablesList[0], twoDMap, 
                                               mosCount, combMethod,
                                               &combParameter, offSets, 
                                               0, 1, 0))) {
          cpl_msg_error(task, "Failure in creating combined 2D sky image");
          error = 1;
        }
      }
    }
  }

  for (j = 0; j < mosCount; j++)
    deleteImage(mosList2D[j]);
  cpl_free(mosList2D);
  for (j = 0; j < mosCount; j++)
    deleteImage(skyList2D[j]);
  cpl_free(skyList2D);
  for (j = 0; j < mosCount; j++)
    deleteWindowTable(winTablesList[j]);
  cpl_free(winTablesList);

  if (error) {
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }


  /*
   * Determine sky alignment
   */

 /*
  *  Line Catalog
  */

  if ((lineCatFrame = pilSofLookup(sof, lineCatalogCategory))) {
    pilFrmSetType(lineCatFrame, PIL_FRAME_TYPE_CALIB);
    if ((lineCatFile = openOldFitsFile(pilFrmGetName(lineCatFrame), 0, 0))) {
      if ((lineCat = newLineCatalog())) {
        if (readFitsLineCatalog(lineCat, lineCatFile->fptr) == VM_TRUE) {
          closeFitsImage(lineCatFile, 0);
          readDoubleDescriptor(outSpSkyExStack[1]->descs, pilTrnGetKeyword("Cdelt", 1), &meanFwhm, 0);
          meanFwhm *= 4;
          rmsValue = distortionsRms(outSpSkyExStack[1], lineCat, meanFwhm);
          cpl_msg_info(task, "IDS RMS after alignment: %.3f pixel", rmsValue);
        }
        else {
          cpl_msg_error(task, "Failure reading line catalog");
          deleteTable(lineCat);
        }
      }
      else
        cpl_msg_error(task, "Not enough memory");
    }
    else
      cpl_msg_error(task, "Failure opening line catalog");
  }


  /* 
   *  Detection on sky subtracted image 
   */

  copySpStack2D = duplicateImage(outSpSkyExStack[0]);
  copyAllDescriptors(outSpSkyExStack[0]->descs, &(copySpStack2D)->descs);

  specFrac = 0.9;

  if (!(combWindowTable = VmSpDetObj(copySpStack2D, extTablesList[0], 
                                     numLevels, slitMargin, detLevel, 
                                     objFrac, limFrac, minObjectSize, 
                                     minCompositeSize, slitTolerance,
                                     specFrac))) {
    cpl_msg_error(task, "Failure deriving the window table for 2D comb image");
    deleteImage(copySpStack2D);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  vimosDscCopy(&combWindowTable->descs, copySpStack2D->descs,
               pilTrnGetKeyword("PROG.ID"), NULL);

  vimosDscCopy(&combWindowTable->descs, copySpStack2D->descs,
               pilTrnGetKeyword("OBS.ID"), NULL);

  for (j = 0; j < mosCount; j++)
    deleteExtractionTable(extTablesList[j]);
  cpl_free(extTablesList);

  deleteImage(copySpStack2D);


  /*
   *  1D extraction
   */

  if (numObjsInWindowTable(combWindowTable) < 1) {
    cpl_msg_error(task, "No objects found!");
    deleteWindowTable(combWindowTable);
    deleteImage(outSpSkyExStack[0]);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  error = 1;
  cpl_msg_info(task, "Extracting 1D spectra...");
  if ((objectTable = newObjectTable())) {
    if ((imaSpEx1D = VmSpEx1D(outSpSkyExStack, combWindowTable, 
                              objectTable, horneExtraction, mosCount))) { 
      error = 0;
    }
    else {
      cpl_msg_error(task, "Failure in 1D extraction of spectra");
      deleteObjectTable(objectTable);
    }
  }
  else
    cpl_msg_error(task, "Failure creating the object table");

  vimosDscCopy(&objectTable->descs, combWindowTable->descs,
               pilTrnGetKeyword("PROG.ID"), NULL);

  vimosDscCopy(&objectTable->descs, combWindowTable->descs,
               pilTrnGetKeyword("OBS.ID"), NULL);


  if (error) {
    deleteImage(outSpSkyExStack[0]);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }


  /* 
   * Spectrophotometric calibration to 1D spectrum
   */ 

  if (calibrateFlux) {

    cpl_msg_info(task, "Applying spectro-photometric calibration ...");
  

    /* 
     * Spectrophotometric table opening 
     */
  
    error = 1;
    if ((sphotFile = openOldFitsFile(pilFrmGetName(sphotFrame), 0, 0))) {
      if ((sphotTable = newSpecPhotTableEmpty())) {
        if (readFitsSpecPhotTable(sphotTable, sphotFile->fptr) == VM_TRUE) {
          closeFitsImage(sphotFile, 0);
          error = 0;
        }
        else {
          cpl_msg_error(task, "Failure reading spectro-photometric table");
          deleteTable(sphotTable);
        }
      }
      else
        cpl_msg_error(task, "Not enough memory");

      deleteImage(sphotFile);

    }
    else 
      cpl_msg_error(task, "Failure in opening spectro-photometric table");
  
    if (error) {
      deleteImage(imaSpEx1D[0]);
      deleteImage(imaSpEx1D[1]);
      cpl_free(imaSpEx1D);
      deleteImage(outSpSkyExStack[0]);
      deleteImage(outSpSkyExStack[1]);
      cpl_free(outSpSkyExStack);
      deleteObjectTable(objectTable);
      deleteWindowTable(combWindowTable);
      deleteImage(imaFringes2D);
      deleteImage(imaFringes);
      return EXIT_FAILURE;
    }


    /*
     * Load the atmospheric extinction table
     */

    error = 1;
    if ((atmFile = openOldFitsFile(pilFrmGetName(atmFrame), 0, 0))) {
      if ((atmTable = newExtinctTableEmpty())) {
        if (readFitsExtinctTable(atmTable, atmFile->fptr) == VM_TRUE) {
          closeFitsImage(atmFile, 0);
          error = 0;
        }
        else {
          cpl_msg_error(task, "Failure reading atmospheric extinction table");
          deleteTable(atmTable);
        }
      }
      else
        cpl_msg_error(task, "Not enough memory");

      deleteImage(atmFile);

    }
    else
      cpl_msg_error(task, "Failure in opening atmospheric extinction table");

    if (error) {
      deleteTable(sphotTable);
      deleteImage(imaSpEx1D[0]);
      deleteImage(imaSpEx1D[1]);
      cpl_free(imaSpEx1D);
      deleteImage(outSpSkyExStack[0]);
      deleteImage(outSpSkyExStack[1]);
      cpl_free(outSpSkyExStack);
      deleteObjectTable(objectTable);
      deleteWindowTable(combWindowTable);
      deleteImage(imaFringes2D);
      deleteImage(imaFringes);
      return EXIT_FAILURE;
    }

    imaSpEx1DCal = VmSpApplyPhot(imaSpEx1D[0], sphotTable, atmTable);

    if (imaSpEx1DCal == NULL) {
      cpl_msg_warning(task, "Spectro-photometric calibration failure: "
                      "the corresponding product, %s, will not be created", 
                      mosFluxReducedCategory);
      calibrateFlux = 0;
    }

    deleteTable(sphotTable);
    deleteTable(atmTable);

  }


  /*
   * Update the products header
   *
   * Note that for the moment also keywords which are not task specific
   * are handled here, since this is the last possibility to access
   * the linked list of keywords without reopening the file.
   * This may change in future!
   */

  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("ExposureTime"),
                         expTime,
                         pilTrnGetComment("ExposureTime"),
                         "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("SummedExposureTime"),
                         sumExpTime,
                         pilTrnGetComment("SummedExposureTime"),
                         "ESO*", 1);

  updateOK = updateOK && insertIntDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("NFramesCombined"),
                         mosCount,
                         pilTrnGetComment("NFramesCombined"),
                         "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(imaSpEx1D[0]),
                         pilTrnGetComment("DataMin"), 
                         "ESO*", 1);
  
  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(imaSpEx1D[0]),
                         pilTrnGetComment("DataMax"), 
                         "ESO*", 1);
  
  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(imaSpEx1D[0]),
                         pilTrnGetComment("DataMean"));
  
  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(imaSpEx1D[0]),
                         pilTrnGetComment("DataStdDeviation"));
  
  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(imaSpEx1D[0]),
                         pilTrnGetComment("DataMedian")); 
  
/*
  updateOK = updateOK && insertStringDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("DoCategory"), 
                         mosReducedCategory,
                         pilTrnGetComment("DoCategory"), 
                         "ESO*", 1);
*/
 
  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product 1 header");
    deleteImage(imaSpEx1DCal);
    deleteImage(imaSpEx1D[0]);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[0]);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }
  
  deleteSetOfDescriptors(&(imaSpEx1D[0]->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(imaSpEx1D[0]->descs), "ESO DPR*");
  deleteSetOfDescriptors(&(imaSpEx1D[0]->descs), "CD1_*");
  deleteSetOfDescriptors(&(imaSpEx1D[0]->descs), "CD2_*");
                                         /* to use IRAF splot task */
  
  
  removeDescriptor(&(imaSpEx1D[0]->descs), 
                   pilTrnGetKeyword("TplExposureNumber"));
  
  vmstrlower(strcpy(output1DName, mosReducedCategory));
  strcat(output1DName, ".fits");

/* Align WCS to convention used with IFU */

  readDoubleDescriptor(imaSpEx1D[0]->descs,
                       pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
  readDoubleDescriptor(imaSpEx1D[0]->descs,
                       pilTrnGetKeyword("Crval", 1), &crval, 0);
  crval -= cdelt / 2;
  writeDoubleDescriptor(&(imaSpEx1D[0]->descs),
                        pilTrnGetKeyword("Crval", 1), crval,
                        pilTrnGetComment("Crval"));
/* */

  if (createFitsImage(output1DName, imaSpEx1D[0], mosReducedCategory)) {
    outputFrame = newPilFrame(output1DName, mosReducedCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", output1DName);
    deleteImage(imaSpEx1DCal);
    deleteImage(imaSpEx1D[0]);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[0]);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }
  
  deleteImage(imaSpEx1D[0]);

  if (calibrateFlux) {

    updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("ExposureTime"),
                           1.0,
                           pilTrnGetComment("ExposureTime"),
                           "ESO*", 1);

    updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("SummedExposureTime"),
                           sumExpTime,
                           pilTrnGetComment("SummedExposureTime"),
                           "ESO*", 1);

    updateOK = updateOK && insertIntDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("NFramesCombined"),
                           mosCount,
                           pilTrnGetComment("NFramesCombined"),
                           "ESO*", 1);

    updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("DataMin"),
                           imageMinimum(imaSpEx1DCal),
                           pilTrnGetComment("DataMin"),
                           "ESO*", 1);

    updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("DataMax"),
                           imageMaximum(imaSpEx1DCal),
                           pilTrnGetComment("DataMax"),
                           "ESO*", 1);

    updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("DataMean"),
                           imageMean(imaSpEx1DCal),
                           pilTrnGetComment("DataMean"));

    updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("DataStdDeviation"),
                           imageSigma(imaSpEx1DCal),
                           pilTrnGetComment("DataStdDeviation"));

    updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("DataMedian"),
                           imageMedian(imaSpEx1DCal),
                           pilTrnGetComment("DataMedian"));

    updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1DCal->descs),
                           pilTrnGetKeyword("AirMass"), 0.0,
                           pilTrnGetComment("AirMass"));

    if (!updateOK) {
      cpl_msg_error(task, "Failure updating product 1 header");
      deleteImage(imaSpEx1DCal);
      deleteImage(imaSpEx1D[1]);
      cpl_free(imaSpEx1D);
      deleteImage(outSpSkyExStack[0]);
      deleteImage(outSpSkyExStack[1]);
      cpl_free(outSpSkyExStack);
      deleteObjectTable(objectTable);
      deleteWindowTable(combWindowTable);
      deleteImage(imaFringes2D);
      deleteImage(imaFringes);
      return EXIT_FAILURE;
    } 

    deleteSetOfDescriptors(&(imaSpEx1DCal->descs), "ESO ADA*");
    deleteSetOfDescriptors(&(imaSpEx1DCal->descs), "ESO DPR*");
    deleteSetOfDescriptors(&(imaSpEx1DCal->descs), "CD1_*");
    deleteSetOfDescriptors(&(imaSpEx1DCal->descs), "CD2_*");
                                           /* to use IRAF splot task */

    removeDescriptor(&(imaSpEx1DCal->descs),
                     pilTrnGetKeyword("TplExposureNumber"));

    vmstrlower(strcpy(output1DName, mosFluxReducedCategory));
    strcat(output1DName, ".fits");

/* Align WCS to convention used with IFU */
  
  readDoubleDescriptor(imaSpEx1DCal->descs,
                       pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
  readDoubleDescriptor(imaSpEx1DCal->descs,
                       pilTrnGetKeyword("Crval", 1), &crval, 0);
  crval -= cdelt / 2;
  writeDoubleDescriptor(&(imaSpEx1DCal->descs),
                        pilTrnGetKeyword("Crval", 1), crval,
                        pilTrnGetComment("Crval"));
/* */

    if (createFitsImage(output1DName, imaSpEx1DCal, mosFluxReducedCategory)) {
      outputFrame = newPilFrame(output1DName, mosFluxReducedCategory);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
    }
    else {
      cpl_msg_error(task, "Cannot create local product file %s!", output1DName);
      deleteImage(imaSpEx1DCal);
      deleteImage(imaSpEx1D[1]);
      cpl_free(imaSpEx1D);
      deleteImage(outSpSkyExStack[0]);
      deleteImage(outSpSkyExStack[1]);
      cpl_free(outSpSkyExStack);
      deleteObjectTable(objectTable);
      deleteWindowTable(combWindowTable);
      deleteImage(imaFringes2D);
      deleteImage(imaFringes);
      return EXIT_FAILURE;
    }

    deleteImage(imaSpEx1DCal);

  }


  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("ExposureTime"),
                         expTime,
                         pilTrnGetComment("ExposureTime"),
                         "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("SummedExposureTime"),
                         sumExpTime,
                         pilTrnGetComment("SummedExposureTime"),
                         "ESO*", 1);

  updateOK = updateOK && insertIntDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("NFramesCombined"),
                         mosCount,
                         pilTrnGetComment("NFramesCombined"),
                         "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(outSpSkyExStack[0]),
                         pilTrnGetComment("DataMin"), 
                         "ESO*", 1);
  
  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(outSpSkyExStack[0]),
                         pilTrnGetComment("DataMax"), 
                         "ESO*", 1);
  
  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(outSpSkyExStack[0]),
                         pilTrnGetComment("DataMean"));
  
  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(outSpSkyExStack[0]),
                         pilTrnGetComment("DataStdDeviation"));
  
  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(outSpSkyExStack[0]),
                         pilTrnGetComment("DataMedian"));
  
/*
  updateOK = updateOK && insertStringDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("DoCategory"), 
                         mosExtractedCategory,
                         pilTrnGetComment("DoCategory"), 
                         "ESO*", 1);
*/
  
  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product 2 header");
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[0]);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  deleteSetOfDescriptors(&(outSpSkyExStack[0]->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(outSpSkyExStack[0]->descs), "ESO DPR*");

  deleteSetOfDescriptors(&(outSpSkyExStack[0]->descs), "CD1_*");  
  deleteSetOfDescriptors(&(outSpSkyExStack[0]->descs), "CD2_*");
                                            /* to use IRAF splot task */

  removeDescriptor(&(outSpSkyExStack[0]->descs),
                   pilTrnGetKeyword("TplExposureNumber"));

  vmstrlower(strcpy(output2DName, mosExtractedCategory));
  strcat(output2DName, ".fits");

/* Align WCS to convention used with IFU */

  readDoubleDescriptor(outSpSkyExStack[0]->descs,
                       pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
  readDoubleDescriptor(outSpSkyExStack[0]->descs,
                       pilTrnGetKeyword("Crval", 1), &crval, 0);
  crval -= cdelt / 2;
  writeDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                        pilTrnGetKeyword("Crval", 1), crval,
                        pilTrnGetComment("Crval"));
/* */

  if (createFitsImage(output2DName, 
                      outSpSkyExStack[0], mosExtractedCategory)) {
    outputFrame = newPilFrame(output2DName, mosExtractedCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", output2DName);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[0]);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  deleteImage(outSpSkyExStack[0]);


  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("ExposureTime"),
                         expTime,
                         pilTrnGetComment("ExposureTime"),
                         "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("SummedExposureTime"),
                         sumExpTime,
                         pilTrnGetComment("SummedExposureTime"),
                         "ESO*", 1);

  updateOK = updateOK && insertIntDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("NFramesCombined"),
                         mosCount,
                         pilTrnGetComment("NFramesCombined"),
                         "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(outSpSkyExStack[1]),
                         pilTrnGetComment("DataMin"), 
                         "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(outSpSkyExStack[1]),
                         pilTrnGetComment("DataMax"), 
                         "ESO*", 1);

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(outSpSkyExStack[1]),
                         pilTrnGetComment("DataMean"));

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(outSpSkyExStack[1]),
                         pilTrnGetComment("DataStdDeviation"));

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(outSpSkyExStack[1]),
                         pilTrnGetComment("DataMedian"));

/*
  updateOK = updateOK && insertStringDescriptor(&(outSpSkyExStack[1]->descs),
                         pilTrnGetKeyword("DoCategory"), 
                         mosSkyCategory,
                         pilTrnGetComment("DoCategory"), 
                         "ESO*", 1);
*/

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product 3 header");
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  deleteSetOfDescriptors(&(outSpSkyExStack[1]->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(outSpSkyExStack[1]->descs), "ESO DPR*");

  deleteSetOfDescriptors(&(outSpSkyExStack[1]->descs), "CD1_*");
  deleteSetOfDescriptors(&(outSpSkyExStack[1]->descs), "CD2_*");
                                          /* to use IRAF splot task */

  removeDescriptor(&(outSpSkyExStack[1]->descs),
                   pilTrnGetKeyword("TplExposureNumber"));

  vmstrlower(strcpy(output2DSkyName, mosSkyCategory));
  strcat(output2DSkyName, ".fits");

/* Align WCS to convention used with IFU */

  readDoubleDescriptor(outSpSkyExStack[1]->descs,
                       pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
  readDoubleDescriptor(outSpSkyExStack[1]->descs,
                       pilTrnGetKeyword("Crval", 1), &crval, 0);
  crval -= cdelt / 2;
  writeDoubleDescriptor(&(outSpSkyExStack[1]->descs),
                        pilTrnGetKeyword("Crval", 1), crval,
                        pilTrnGetComment("Crval"));
/* */

  if (createFitsImage(output2DSkyName, outSpSkyExStack[1], mosSkyCategory)) {
    outputFrame = newPilFrame(output2DSkyName, mosSkyCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task,
                "Cannot create local product file %s!", output2DSkyName);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  deleteImage(outSpSkyExStack[1]);
  cpl_free(outSpSkyExStack);

  writeStringDescriptor(&(objectTable->descs), pilTrnGetKeyword("DoCategory"),
                        objectTableCategory, "Category of pipeline product");

  deleteSetOfDescriptors(&(objectTable->descs), "ESO QC*");

  vmstrlower(strcpy(objectTableName, objectTableCategory));
  /* strcat(objectTableName, ".TFITS"); */
  strcat(objectTableName, ".fits");

  error = 1;

  if ((objectFile = newImage(0, 0, NULL))) {
    if (VM_TRUE == openNewFitsImage(objectTableName, objectFile)) {
      if (VM_TRUE == writeFitsObjectTable(objectTable, objectFile->fptr)) {
        closeFitsImage(objectFile, 0);
        pilFitsHdrCopy(objectTableName, 0, NULL, ".*-OBS$", 1);
        pilFitsHdrCopy(objectTableName, 0, NULL, "^ESO .*", 1);
        outputFrame = newPilFrame(objectTableName, objectTableCategory);

        pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
        pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
        pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
        pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

        pilSofInsert(sof, outputFrame);
        error = 0;
      }
    }
  }
  else {
    cpl_msg_error(task,
                "Cannot create local product file %s!", objectTableName);
  }

  deleteImage(objectFile);
  deleteObjectTable(objectTable);

  if (error) {
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  writeStringDescriptor(&(combWindowTable->descs), 
                        pilTrnGetKeyword("DoCategory"), windowTableCategory, 
                        "Category of pipeline product");

  deleteSetOfDescriptors(&(combWindowTable->descs), "ESO QC*");

  vmstrlower(strcpy(windowTableName, windowTableCategory));
  /* strcat(windowTableName, ".TFITS"); */
  strcat(windowTableName, ".fits");

  error = 1;

  if ((windowFile = newImage(0, 0, NULL))) {
    if (VM_TRUE == openNewFitsImage(windowTableName, windowFile)) {
      if (VM_TRUE == writeFitsWindowTable(combWindowTable, windowFile->fptr)) {
        closeFitsImage(windowFile, 0);
        pilFitsHdrCopy(windowTableName, 0, NULL, ".*-OBS$", 1);
        pilFitsHdrCopy(windowTableName, 0, NULL, "^ESO .*", 1);
        outputFrame = newPilFrame(windowTableName, windowTableCategory);

        pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
        pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
        pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
        pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

        pilSofInsert(sof, outputFrame);
        error = 0;
      }
    }
  }
  else {
    cpl_msg_error(task,
                "Cannot create local product file %s!", windowTableName);
  }

  deleteImage(windowFile);
  deleteWindowTable(combWindowTable);

  if (error) {
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("ExposureTime"),
                         expTime,
                         pilTrnGetComment("ExposureTime"),
                         "ESO*", 1);
  
  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("SummedExposureTime"),
                         sumExpTime,
                         pilTrnGetComment("SummedExposureTime"),
                         "ESO*", 1);
  
  updateOK = updateOK && insertIntDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("NFramesCombined"),
                         mosCount,
                         pilTrnGetComment("NFramesCombined"),
                         "ESO*", 1);
  
  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(imaSpEx1D[1]),
                         pilTrnGetComment("DataMin"),
                         "ESO*", 1);
  
  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(imaSpEx1D[1]),
                         pilTrnGetComment("DataMax"),
                         "ESO*", 1);
  
  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(imaSpEx1D[1]),
                         pilTrnGetComment("DataMean"));

  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(imaSpEx1D[1]),
                         pilTrnGetComment("DataStdDeviation"));

  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(imaSpEx1D[1]),
                         pilTrnGetComment("DataMedian"));

/*
  updateOK = updateOK && insertStringDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DoCategory"),
                         mosSkyReducedCategory,
                         pilTrnGetComment("DoCategory"),
                         "ESO*", 1);
*/

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product 5 header");
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  deleteSetOfDescriptors(&(imaSpEx1D[1]->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(imaSpEx1D[1]->descs), "ESO DPR*");
  deleteSetOfDescriptors(&(imaSpEx1D[1]->descs), "CD1_*");
  deleteSetOfDescriptors(&(imaSpEx1D[1]->descs), "CD2_*");
                                         /* to use IRAF splot task */


  removeDescriptor(&(imaSpEx1D[1]->descs),
                   pilTrnGetKeyword("TplExposureNumber"));

  vmstrlower(strcpy(output1DName, mosSkyReducedCategory));
  strcat(output1DName, ".fits");

/* Align WCS to convention used with IFU */

  readDoubleDescriptor(imaSpEx1D[1]->descs,
                       pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
  readDoubleDescriptor(imaSpEx1D[1]->descs,
                       pilTrnGetKeyword("Crval", 1), &crval, 0);
  crval -= cdelt / 2;
  writeDoubleDescriptor(&(imaSpEx1D[1]->descs),
                        pilTrnGetKeyword("Crval", 1), crval,
                        pilTrnGetComment("Crval"));
/* */

  if (createFitsImage(output1DName, imaSpEx1D[1], mosSkyReducedCategory)) {
    outputFrame = newPilFrame(output1DName, mosSkyReducedCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", output1DName);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[0]);
    deleteImage(outSpSkyExStack[1]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    deleteImage(imaFringes2D);
    deleteImage(imaFringes);
    return EXIT_FAILURE;
  }

  if (fringingCorr == 1 && friCount == 0) {

    if (fringMethod == 0) {
      tmpImage = imaFringes;
      mosFringesCategory = mosFringesSkyCategory;
    }
    else
      tmpImage = imaFringes2D;

    updateOK = updateOK && insertDoubleDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("ExposureTime"),
                           expTime,
                           pilTrnGetComment("ExposureTime"),
                           "ESO*", 1);

    updateOK = updateOK && insertDoubleDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("SummedExposureTime"),
                           sumExpTime,
                           pilTrnGetComment("SummedExposureTime"),
                           "ESO*", 1);

    updateOK = updateOK && insertIntDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("NFramesCombined"),
                           mosCount,
                           pilTrnGetComment("NFramesCombined"),
                           "ESO*", 1);

    updateOK = updateOK && insertDoubleDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("DataMin"),
                           imageMinimum(tmpImage),
                           pilTrnGetComment("DataMin"),
                           "ESO*", 1);

    updateOK = updateOK && insertDoubleDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("DataMax"),
                           imageMaximum(tmpImage),
                           pilTrnGetComment("DataMax"),
                           "ESO*", 1);

    updateOK = updateOK && writeDoubleDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("DataMean"),
                           imageMean(tmpImage),
                           pilTrnGetComment("DataMean"));

    updateOK = updateOK && writeDoubleDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("DataStdDeviation"),
                           imageSigma(tmpImage),
                           pilTrnGetComment("DataStdDeviation"));
  
    updateOK = updateOK && writeDoubleDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("DataMedian"),
                           imageMedian(tmpImage),
                           pilTrnGetComment("DataMedian"));

    updateOK = updateOK && insertStringDescriptor(&(tmpImage->descs),
                           pilTrnGetKeyword("DoCategory"),
                           mosFringesCategory,
                           pilTrnGetComment("DoCategory"),
                           "ESO*", 1);

    if (!updateOK) {
      cpl_msg_error(task, "Failure updating product 4 header");
      deleteImage(tmpImage);
      deleteImage(imaSpEx1D[1]);
      cpl_free(imaSpEx1D);
      return EXIT_FAILURE;
    }

    deleteSetOfDescriptors(&(tmpImage->descs), "ESO ADA*");
    deleteSetOfDescriptors(&(tmpImage->descs), "ESO DPR*");
  
    deleteSetOfDescriptors(&(tmpImage->descs), "CD1_*");
    deleteSetOfDescriptors(&(tmpImage->descs), "CD2_*");
                                            /* to use IRAF splot task */

    removeDescriptor(&(tmpImage->descs),
                     pilTrnGetKeyword("TplExposureNumber"));

    vmstrlower(strcpy(outputFringesName, mosFringesCategory));
    strcat(outputFringesName, ".fits");

    if (fringMethod == 1) {
      /* Align WCS to convention used with IFU */

      readDoubleDescriptor(imaSpEx1D[1]->descs,
                           pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
      readDoubleDescriptor(imaSpEx1D[1]->descs,
                           pilTrnGetKeyword("Crval", 1), &crval, 0);
      crval -= cdelt / 2;
      writeDoubleDescriptor(&(imaSpEx1D[1]->descs),
                            pilTrnGetKeyword("Crval", 1), crval,
                            pilTrnGetComment("Crval"));
      /* */
    }

    if (createFitsImage(outputFringesName, tmpImage, mosFringesCategory)) {
      outputFrame = newPilFrame(outputFringesName, mosFringesCategory);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);
  
      pilSofInsert(sof, outputFrame);
    }
    else {
      cpl_msg_error(task,
                  "Cannot create local product file %s!", outputFringesName);
      deleteImage(tmpImage);
      deleteImage(imaSpEx1D[1]);
      cpl_free(imaSpEx1D);
      return EXIT_FAILURE;
    }

    deleteImage(tmpImage);

  }

  deleteImage(imaSpEx1D[1]);
  cpl_free(imaSpEx1D);

  return EXIT_SUCCESS;

}


