/* $Id: vmmosstandard.c,v 1.4 2011-03-14 15:28:58 cizzo Exp $
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
 * $Author: cizzo $
 * $Date: 2011-03-14 15:28:58 $
 * $Revision: 1.4 $
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
#include <piltranslator.h>
#include <pilrecipe.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilutils.h>
#include <pilfits.h>
#include <pilqc.h>
#include <vmqcutils.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmgrismtable.h"
#include "vmextractiontable.h"
#include "vmwindowtable.h"
#include "vmobjecttable.h"
#include "vmextincttable.h"
#include "vmspecphottable.h"
/* #include "vmsphottable.h" */
#include "vmstdfluxtable.h"
#include "vmfit.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmmosflat.h"
#include "vmmoswavecalib.h"
#include "vmmossky.h"
#include "vmmosextraction.h"
/* #include "vmmossphotcalib.h"  */
#include "vmmosutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


#define MAX_COMMENT_LENGTH (80)


static cxint vmmosstandard(PilSetOfFrames *);


/*
 * Definition of the label strings for all methods the recipe function
 * supports for combining frames and bias removal, with their associated
 * method code. 
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
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it availble to the application using the interface.
 */

static cxint
vmmosstandard_create(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;

    cpl_parameter *p;

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

    p = cpl_parameter_new_enum("vimos.Parameters.bias.removing.method",
                                CPL_TYPE_STRING,
                                "Bias removal method.",
                                "vimos.Parameters",
                                "Zmaster", 2, "Zmaster", "Master");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "BiasMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "BiasMethod");
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
                                0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PolyOrder");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PolyOrder");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.sky.ksigma.low",
                                CPL_TYPE_DOUBLE,
                                "Low threshold for K-sigma rejection "
                                "in sky fitting.",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SkyKSigmaLow");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SkyKSigmaLow");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.sky.ksigma.high",
                                CPL_TYPE_DOUBLE,
                                "High threshold for K-sigma rejection "
                                "in sky fitting.",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SkyKSigmaHigh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SkyKSigmaHigh");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_enum("vimos.Parameters.sampling",
                                CPL_TYPE_STRING,
                                "Spectrum sampling in wavelength during "
                                "extraction.",
                                "vimos.Parameters",
                                "Linear", 2, "Linear", "Log");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SamplingMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SamplingMethod");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

    p = cpl_parameter_new_value("vimos.Parameters.extraction.fuzz",
                                CPL_TYPE_INT,
                                "Extra pixels from expected position of "
                                "spectrum edge in extraction.",
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
                                "Size of spectrum to extract around any "
                                "skyline.",
                                "vimos.Parameters",
                                16);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LineWidth");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LineWidth");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.detection.sigma",
                                CPL_TYPE_DOUBLE,
                                "Object detection level in units of sigma.",
                                "vimos.Parameters",
                                10.0);
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
                                "Minimal size for an object candidate to be "
                                "considered an object.",
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
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.extraction.optimal",
                                CPL_TYPE_BOOL,
                                "Use 1D Horne extraction.",
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


    p = cpl_parameter_new_value("vimos.Parameters.response.order",
                                CPL_TYPE_INT,
                                "Order of the polynomial used to smooth "
                                "the instrument response curve.",
                                "vimos.Parameters",
                                5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ResponseOrder");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ResponseOrder");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.reduceall",
                                CPL_TYPE_BOOL,
                                "Reduce any frame.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ReduceAnyFrame");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ReduceAnyFrame");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.quality.enable",
                                CPL_TYPE_BOOL,
                                "Compute QC1 parameters",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ComputeQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ComputeQC");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.slit.select",
                                CPL_TYPE_INT,
                                "Slit containing the standard star.",
                                "vimos.Parameters",
                                0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SelectSlit");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SelectSlit");
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
vmmosstandard_exec(cpl_plugin *plugin)
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

    if (vmmosstandard(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmmosstandard");
        
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
vmmosstandard_destroy(cpl_plugin *plugin)
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
 *   Reduce a single MOS standard star exposure, and produce the
 *   response and the efficiency curves.
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param sof  Set of frames with references to a MOS standard star exposure, 
 *             a master bias, a grism table, and catalog spectral flux from 
 *             the same observed standard star. A spectral master flat field, 
 *             a CCD table, an extraction table, and a master dark, may be 
 *             given. If an extraction table is specified, all the spectral 
 *             distorsion models are read from its header, otherwise they 
 *             are read from the header of the science exposure (first 
 *             guesses). If bad pixel cleaning is requested, a CCD table 
 *             must be given.
 *
 * @doc 
 *   Subtract bias and dark and apply flat field (if given); calibrate 
 *   shift of the slit positions using skylines listed in the Extraction 
 *   Table; detect objects in frame; one and not more than one object
 *   must be detected. Slits are classified in long and short ones; 
 *   subtract sky from short slits; 2D extraction and wavelength 
 *   calibration of object spectrum and sky spectra; subtract sky 
 *   from long slits; extract 1D spectra of the object and derive 
 *   the response and the efficiency curves.
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
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
 *     \item CleanBadPixel:      Bad pixel correction on input MOS exposure.
 *                               If this option is turned on, the recipe
 *                               expects a CCD\_TABLE in the input set of
 *                               frames.
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
 *    \end{itemize}
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 * 
 * @author P. Sartoretti, C. Izzo, R. Palsa
 */   

static cxint 
vmmosstandard(PilSetOfFrames *sof)
{

  const char  task[] = "vmmosstandard";

  const char  parameter[]          = "Parameters";

  const char *mosCategory          = pilTrnGetCategory("MosStandard");
  const char *mosExtractedCategory = pilTrnGetCategory("MosStandardExtracted");
  const char *mosReducedCategory   = pilTrnGetCategory("MosStandardReduced");
  const char *mosSkyCategory       = pilTrnGetCategory("MosStandardSky");
  const char *mosSkyReducedCategory= pilTrnGetCategory("MosStandardSkyReduced");
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
  const char *fluxTableCategory    = pilTrnGetCategory("StdFluxTable");

  char          *biasMethodTag = NULL;
  BiasMethod     biasMethod = BIAS_UNDEF;
  int            biasMethodEntry;

  char          *skyMethodTag = NULL;
  SkyMethod      skyMethod = SKY_UNDEF;
  int            skyMethodEntry;

#ifdef ONLINE_MODE

  char          *samplingMethodTag = NULL;
  int            samplingMethodEntry;

#endif

  SamplingMethod samplingMethod = VM_SP_LIN_LAMBDA;

  char           output2DName[PATHNAME_MAX + 1];
  char           output2DSkyName[PATHNAME_MAX + 1];
  char           output1DName[PATHNAME_MAX + 1];
  char           objectTableName[PATHNAME_MAX + 1];
  char           windowTableName[PATHNAME_MAX + 1];
  char           sphotTableName[PATHNAME_MAX + 1];
  char           parName[30]; 

  VimosBool      updateOK = VM_TRUE;

  unsigned int   computeQC;
  unsigned int   cleanBadPixel;
  unsigned int   reduceAnyFrame;
  unsigned int   horneExtraction;
  unsigned int   modelSlit;
  unsigned int   useSkylines;
  unsigned int   error;

  int            order            = 0;
  int            responseOrder    = 0;
  int            selectSlit       = 0;
  int            polyDeg          = 0;
  int            fuzz             = 0;
  int            slitMargin       = 0;
  int            lineWidth        = 0;
  int            numLevels        = 0;
  int            minObjectSize    = 0;
  int            minCompositeSize = 0;

  float          slitTolerance    = 0.0;
  float          detLevel         = 0.0;
  float          objFrac          = 0.0;
  float          specFrac         = 0.0;

  double         time, gain, airmass;
  double         start, step;
  double         wstart, wstep;
  int            wcount;

  float          skySigmaLow, skySigmaHigh;

  PilFrame      *biasFrame, *darkFrame, *ccdFrame, *grismFrame, *mosFrame;
  PilFrame      *flatFrame, *extrFrame, *outputFrame;
  PilFrame      *atmFrame, *fluxFrame;

  VimosImage           *grismFile;
  VimosImage           *extrFile;
  VimosImage           *windowFile = NULL;
  VimosImage           *objectFile = NULL;
  VimosImage           *biasImage = NULL;
  VimosImage           *darkImage = NULL;
  VimosImage           *flatImage = NULL;
  VimosImage           *atmFile = NULL;
  VimosImage           *fluxFile = NULL;
  VimosImage           *mosImage  = NULL;
  VimosImage          **imaSpEx1D = NULL;
  VimosImage           *imaSpEx1DCal = NULL;
  VimosImage          **outSpSkyFra = NULL;
  VimosImage          **outSpEx2D = NULL;
  VimosImage          **outSpSkyEx = NULL;
  VimosImage           *tmpImage;
  VimosTable           *ccdTable = NULL;
  VimosTable           *grismTable = NULL;
  VimosTable           *atmTable = NULL;
  VimosTable           *fluxTable = NULL;
  VimosTable           *sphotTable = NULL;
  VimosExtractionTable *extractionTable = NULL;
  VimosWindowTable     *windowTable, *dummyWindow;
  VimosWindowSlit      *windowSlit;
  VimosObjectTable     *objectTable;

  VimosImage           *spectrum;
  VimosImage           *extinction;
  VimosImage           *smo_response;
  VimosImage           *response;
  VimosImage           *smo_efficiency;
  VimosImage           *efficiency;
  VimosImage           *stdFlux;
  VimosImage           *convStdFlux;


  /* FIXME:
   * This parameter is still unused:
   */

  int     i, count, preCount, pos, prePos, quadrant, nexp, expno;
  int     newStart, newEnd;
  float   limFrac = 0.0;
  double  lambda;
  double *data;

  double         cdelt, crval;


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
   * Determine the sky level determination method.
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
    polyDeg = pilDfsDbGetInt(parameter, "PolyOrder", 2);
    if (polyDeg < 0) {
      cpl_msg_error(task, "Invalid choice of polynomial for sky "
                  "level modeling: degree should be greater than zero");
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
   * no matter what), we set internally a negative tolerance.
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
   * Get order of polynomial used to smooth the response function
   */

  responseOrder = pilDfsDbGetInt(parameter, "ResponseOrder", 5);

  if (responseOrder < 2) {
    cpl_msg_error(task, "Order of the polynomial fitting the instrument "
                  "response must be at least 2");
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

  cleanBadPixel = pilDfsDbGetBool(parameter, "CleanBadPixel", 1);


  /*
   * Check what kind of 1D extraction should be applied
   */

  horneExtraction = pilDfsDbGetBool(parameter, "HorneExtraction", 0);


  /*
   * Optional modeling of wavelength solution within each slit
   */

  modelSlit = pilDfsDbGetBool(parameter, "ModelSlit", 0);


  /*
   * Order of polynomial for modeling of wavelength solution within each slit.
   */

  order = pilDfsDbGetInt(parameter, "ModelSlitOrder", 0);


  /*
   * Order of polynomial for modeling of wavelength solution within each slit.
   */

  selectSlit = pilDfsDbGetInt(parameter, "SelectSlit", 0);

  if (order < 0) {
    cpl_msg_error(task, "Invalid order for wavelength solution modeling.");
    return EXIT_FAILURE;
  }


 /*
  * Check if the wavelength calibration should be refined on skylines.
  */

  useSkylines = pilDfsDbGetBool(parameter, "UseSkylines", 1);


  /*
   * Check if only frames with TPL EXPNO = quadrant number should be reduced.
   * This is effective only if the frame is part of a template sequence.
   */

  reduceAnyFrame = pilDfsDbGetBool("Parameters", "ReduceAnyFrame", 0);


  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   * Load input frames
   */

  cpl_msg_info(task, "Loading input frames...");


  /*
   * Get raw MOS frame
   */

  if ((mosFrame = pilSofLookup(sof, mosCategory))) {
    if ((mosImage = openOldFitsFile(pilFrmGetName(mosFrame), 1, 0))) {
      pilFrmSetType(mosFrame, PIL_FRAME_TYPE_RAW);
      closeFitsImage(mosImage, 0);
    }
    else {
      cpl_msg_error(task, "Failure opening MOS science image");
      return EXIT_FAILURE;
    }
  }
  else {
    cpl_msg_error(task, "No MOS science frame found in input");
    return EXIT_FAILURE;
  }


  if (readIntDescriptor(mosImage->descs, pilTrnGetKeyword("Quadrant"),
                        &quadrant, NULL) == VM_FALSE) {
    cpl_msg_error(task, "Integer descriptor %s not found",
                  pilTrnGetKeyword("Quadrant"));

    deleteImage(mosImage);
    return EXIT_FAILURE;
  }


  if (!reduceAnyFrame) {

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

    if (readIntDescriptor(mosImage->descs,
                          pilTrnGetKeyword("TplExposures"),
                          &nexp, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Integer descriptor %s not found",
                    pilTrnGetKeyword("TplExposures"));

      deleteImage(mosImage);
      return EXIT_FAILURE;
    }

    if (nexp == 4) {
      if (readIntDescriptor(mosImage->descs,
                            pilTrnGetKeyword("TplExposureNumber"),
                            &expno, NULL) == VM_FALSE) {
        cpl_msg_error(task, "Integer descriptor %s not found",
                    pilTrnGetKeyword("TplExposureNumber"));

        deleteImage(mosImage);
        return EXIT_FAILURE;
      }

      if (readIntDescriptor(mosImage->descs,
                            pilTrnGetKeyword("Quadrant"),
                            &quadrant, NULL) == VM_FALSE) {
        cpl_msg_error(task, "Integer descriptor %s not found",
                      pilTrnGetKeyword("Quadrant"));

        deleteImage(mosImage);
        return EXIT_FAILURE;
      }

      if (quadrant != expno) {
        cpl_msg_warning(task, "Frame %s is not reduced, because it is not "
                        "expected to contain any standard stars. This recipe "
                        "run will be aborted!", pilFrmGetName(mosFrame));
        deleteImage(mosImage);
        return EXIT_FAILURE;
      }
    }
  }


  /*
   * Get exposure time, gain, airmass
   */

  if (readDoubleDescriptor(mosImage->descs, pilTrnGetKeyword("ExposureTime"),
                           &time, NULL) == VM_TRUE) {
    if (time < MIN_DIVISOR) {
      cpl_msg_error(task, "Zero or negative exposure time in input.");
      deleteImage(mosImage);
      return EXIT_FAILURE;
    }
  }
  else {
    cpl_msg_error(task, "Cannot read %s from input exposure.",
                  pilTrnGetKeyword("ExposureTime"));
    deleteImage(mosImage);
    return EXIT_FAILURE;
  }

  if (readDoubleDescriptor(mosImage->descs, pilTrnGetKeyword("Airmass"),
                           &airmass, NULL) == VM_FALSE) {
    if (readDoubleDescriptor(mosImage->descs, "ESO TEL AIRM START",
                             &airmass, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Cannot read the airmass from input exposure.");
      deleteImage(mosImage);
      return EXIT_FAILURE;
    }
  }

  if (airmass < 0.0) {
    cpl_msg_error(task, "Negative airmass in input.");
    deleteImage(mosImage);
    return EXIT_FAILURE;
  }

  gain = getMeanGainFactor(mosImage);
  if (gain < MIN_DIVISOR) {
    cpl_msg_error(task, "Wrong or no gain factor in input.");
    deleteImage(mosImage);
    return EXIT_FAILURE;
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
    deleteImage(mosImage);
    return EXIT_FAILURE;
  }


  /*
   * Recreate bias overscans using as a reference the MOS frame
   */

  if ((tmpImage = growOverscans(biasImage, mosImage))) {
    if (biasImage != tmpImage) {
      deleteImage(biasImage);
      biasImage = tmpImage;
    }
  }
  else {
    cpl_msg_error(task, "Failure in growing overscans in master bias");
    deleteImage(mosImage);
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

  if (VmSubBias(mosImage, biasImage, biasMethod) == EXIT_FAILURE) {
    cpl_msg_error(task, "Cannot remove bias from raw MOS image");
    deleteImage(mosImage);
    deleteImage(biasImage);
    return EXIT_FAILURE;
  }

  deleteImage(biasImage);


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
      deleteImage(mosImage);
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

    if (VmSubDark(mosImage, darkImage) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot subtract dark from MOS image");
      deleteImage(mosImage);
      deleteImage(darkImage);
      return EXIT_FAILURE;
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
      deleteImage(mosImage);
      return EXIT_FAILURE;
    }
  }


  /*
   * Get the (optional) extraction table.
   */

  if ((extrFrame = pilSofLookup(sof, extrTableCategory))) {

    pilFrmSetType(extrFrame, PIL_FRAME_TYPE_CALIB);

    error = 1;

    if ((extrFile = openOldFitsFile(pilFrmGetName(extrFrame), 0, 0))) {
      if ((extractionTable = newExtractionTable())) {
        if (readFitsExtractionTable(extractionTable, extrFile->fptr) 
                                                               == VM_TRUE) {
          closeFitsImage(extrFile, 0);
          error = 0;
        }
        else {
          cpl_msg_error(task, "Failure reading the extraction table");
          deleteExtractionTable(extractionTable);
        }
      }
      else
        cpl_msg_error(task, "Not enough memory");

      deleteImage(extrFile);

    }
    else
      cpl_msg_error(task, "Failure opening the extraction table");
  }

  if (error) {
    deleteImage(mosImage);
    deleteImage(flatImage);
    return EXIT_FAILURE;
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
    deleteExtractionTable(extractionTable);
    deleteImage(flatImage);
    deleteImage(mosImage);
    return EXIT_FAILURE;
  }


  if (!extractionTable) {

    /* 
     * If no extraction table was found in the input Set-Of-Frames, 
     * derive the extraction table from global models. The global
     * distorsion models coefficients are read from the science frame.
     */

    if (!(extractionTable = VmSpExTab(mosImage, grismTable, NULL, NULL))) {
      cpl_msg_error(task, "Cannot compute the extraction table");
      deleteExtractionTable(extractionTable);
      deleteImage(flatImage);
      deleteImage(mosImage);
      deleteTable(grismTable);
      return EXIT_FAILURE;
    }
  }
  else {

    /*
     * Replace in extraction table the fit of the wavelength solution
     * coefficients trend along each slit.
     */

    if (modelSlit)
      modelWavcal(extractionTable, order);

  }

 
  /*
   * Apply flat field (if present).
   */

  if (flatImage) {
    if ((tmpImage = VmSpApplyFF(mosImage, flatImage, extractionTable))) {
      deleteImage(mosImage);
      mosImage = tmpImage;
    }
    else {
      cpl_msg_error(task, "Failure applying flat field to MOS science exposure");
      deleteExtractionTable(extractionTable);
      deleteImage(flatImage);
      deleteImage(mosImage);
      deleteTable(grismTable);
      return EXIT_FAILURE;
    }

    deleteImage(flatImage);
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
      deleteExtractionTable(extractionTable);
      deleteImage(mosImage);
      deleteTable(grismTable);
      return EXIT_FAILURE;
    }

    /*
     *  Image cleaning
     */

    cpl_msg_info(task, "Cleaning bad pixels in MOS science exposure");

    if (EXIT_FAILURE == cleanBadPixels(mosImage, ccdTable, 0)) {
      cpl_msg_error(task, "Cannot clean MOS science exposure");
      deleteExtractionTable(extractionTable);
      deleteImage(mosImage);
      deleteTable(grismTable);
      deleteTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  deleteTable(ccdTable);

  if (useSkylines) {

    cpl_msg_info(task, "Shift wavelength calibration using skylines...");

    if ((VmSpCalShifts(mosImage, grismTable, 
                       extractionTable, 0, lineWidth, fuzz)) == EXIT_FAILURE) {
      cpl_msg_error(task, "Failure correcting offset of slit positions");
      deleteExtractionTable(extractionTable);
      deleteImage(mosImage);
      deleteTable(grismTable);
      return EXIT_FAILURE;
    }

  }

  deleteTable(grismTable);

  if (!(windowTable = VmSpDetObj(mosImage, extractionTable, numLevels, 
                                 slitMargin, detLevel, objFrac, limFrac, 
                                 minObjectSize, minCompositeSize, 
                                 slitTolerance,specFrac))) {
    cpl_msg_error(task, "Failure deriving the window table");
    deleteExtractionTable(extractionTable);
    deleteImage(mosImage);
    return EXIT_FAILURE;
  }

  vimosDscCopy(&windowTable->descs, mosImage->descs,
               pilTrnGetKeyword("PROG.ID"), NULL);

  vimosDscCopy(&windowTable->descs, mosImage->descs,
               pilTrnGetKeyword("OBS.ID"), NULL);
 
  if (!(outSpSkyFra = VmSpSkyFra(mosImage, extractionTable, windowTable, 
                           skyMethod, polyDeg, slitMargin, 10, 3, 
                           skySigmaLow, skySigmaHigh))) {
    cpl_msg_error(task, "Failure in sky subtraction in short slits");
    deleteImage(mosImage);
    deleteExtractionTable(extractionTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }
 

  /*
   * Debug:  *+
     createFitsImage("outSkyFraD.fits", outSpSkyFra[0], "cat");
     createFitsImage("outSkyFraS.fits", outSpSkyFra[1], "cat");
    +* */
  
  outSpEx2D = VmSpEx2D(outSpSkyFra, extractionTable, samplingMethod); 

  deleteImage(outSpSkyFra[0]);
  deleteImage(outSpSkyFra[1]);
  cpl_free(outSpSkyFra);

  if (!outSpEx2D) {
    cpl_msg_error(task, "Failure in 2D extraction of spectra");
    deleteExtractionTable(extractionTable);
    deleteWindowTable(windowTable);
    deleteImage(mosImage);
    return EXIT_FAILURE;
  }

  /* 
   * Debug: *+
     printf("X = %d, Y=%d\n", outSpEx2D[0]->xlen, outSpEx2D[0]->ylen);
     createFitsImage("outsp2dD.fits", outSpEx2D[0], "cat"); 
     createFitsImage("outsp2dS.fits", outSpEx2D[1], "cat"); 
    +* */

  if (!(outSpSkyEx = VmSpSkyExt(outSpEx2D, windowTable, skyMethod, 
                                polyDeg, slitMargin, 10, 3, 
                                skySigmaLow, skySigmaHigh))) {
    cpl_msg_error(task, "Failure in sky subtraction in long slits");
    deleteImage(mosImage);
    deleteImage(outSpEx2D[0]);
    deleteImage(outSpEx2D[1]);
    cpl_free(outSpEx2D);
    deleteExtractionTable(extractionTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }

  deleteImage(outSpEx2D[1]);

  /*
   * Debug: *+
     createFitsImage("outspskyexD.fits", outSpSkyEx[0], "cat"); 
     createFitsImage("outspskyexS.fits", outSpSkyEx[1], "cat"); 
   +* */
 
  /*
   * Here iterate the object finding but this time using the sky-subtracted
   * data. 
   */

  /* FIXME:
   * Question: why iterate, if the windowTable is not modified?
   */

  if (!(dummyWindow = VmSpDetObj(outSpEx2D[0], extractionTable, numLevels, 
                                 slitMargin, detLevel, objFrac, limFrac, 
                                 minObjectSize, minCompositeSize, 
                                 slitTolerance,specFrac))) {
    cpl_msg_error(task, "Iteration on object detection failed");
    deleteImage(mosImage);
    deleteImage(outSpEx2D[0]);
    cpl_free(outSpEx2D);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteExtractionTable(extractionTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }

  deleteWindowTable(dummyWindow);
  deleteExtractionTable(extractionTable);
  deleteImage(outSpEx2D[0]);
  cpl_free(outSpEx2D);

  if (numObjsInWindowTable(windowTable) < 1) {
    deleteImage(mosImage);
    cpl_msg_error(task, "No objects found!");
    deleteWindowTable(windowTable);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    return EXIT_FAILURE;
  }

  if (numObjsInWindowTable(windowTable) > 1) {
    if (selectSlit == 0) {
      cpl_msg_error(task, "More than one object found! Use parameter "
                    "SelectSlit to indicate the slit with the standard star.");
      deleteImage(mosImage);
      deleteWindowTable(windowTable);
      deleteImage(outSpSkyEx[0]);
      deleteImage(outSpSkyEx[1]);
      cpl_free(outSpSkyEx);
      return EXIT_FAILURE;
    }
    else {

      /*
       * Destroy all objects but the one in the selected slit.
       */

      windowSlit = windowTable->slits;
      while (windowSlit) {
        if (windowSlit->slitNo == selectSlit) {
          if (windowSlit->numObj == 1) {
            cpl_msg_info(task, "Selecting just slit %d", selectSlit);
          }
          else {
            if (windowSlit->numObj == 0)
              cpl_msg_error(task, "No object found in selected slit!");
            if (windowSlit->numObj > 1)
              cpl_msg_error(task, "More than one object found in selected "
                            "slit! Perhaps the standard star spectrum is "
                            "contaminated. Modifying the DetectionLevel "
                            "parameter may help...");
            deleteImage(mosImage);
            deleteWindowTable(windowTable);
            deleteImage(outSpSkyEx[0]);
            deleteImage(outSpSkyEx[1]);
            cpl_free(outSpSkyEx);
            return EXIT_FAILURE;
          }
        }
        else {
          if (windowSlit->numObj) {
            deleteWindowObject(windowSlit->objs);
            windowSlit->objs = NULL;
            windowSlit->numObj = 0;
          }
        }
        windowSlit = windowSlit->next;
      }
    }
  }

  error = 1;

  if ((objectTable = newObjectTable())) {
    if ((imaSpEx1D = VmSpEx1D(outSpSkyEx, windowTable, 
                              objectTable, horneExtraction, 1))) {
      error = 0;
    }
    else {
      cpl_msg_error(task, "Failure in 1D extraction of spectra");
      deleteObjectTable(objectTable);
    }
  }
  else
    cpl_msg_error(task, "Failure creating the object table");

  vimosDscCopy(&objectTable->descs, windowTable->descs,
               pilTrnGetKeyword("PROG.ID"), NULL);

  vimosDscCopy(&objectTable->descs, windowTable->descs,
               pilTrnGetKeyword("OBS.ID"), NULL);


  if (error) {
    deleteImage(mosImage);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }

  imaSpEx1DCal = imaSpEx1D[0];


  /*
   * Load the atmospheric extinction table
   */
   
  error = 1;
  
  if ((atmFrame = pilSofLookup(sof, atmTableCategory))) {
    pilFrmSetType(atmFrame, PIL_FRAME_TYPE_CALIB);
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
  }
  else
    cpl_msg_error(task, "No input atmospheric extinction table found");
    
  if (error) {
    deleteImage(mosImage);
    deleteImage(imaSpEx1DCal);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }


  /*
   * Load the standard star flux table
   */

  error = 1;

  if ((fluxFrame = pilSofLookup(sof, fluxTableCategory))) {
    pilFrmSetType(fluxFrame, PIL_FRAME_TYPE_CALIB);
    if ((fluxFile = openOldFitsFile(pilFrmGetName(fluxFrame), 0, 0))) {
      if ((fluxTable = newStdFluxTableEmpty())) {
        if (readFitsStdFluxTable(fluxTable, fluxFile->fptr) == VM_TRUE) {
          closeFitsImage(fluxFile, 0);
          error = 0;
        }
        else {
          cpl_msg_error(task, "Failure reading standard star flux table");
          deleteTable(fluxTable);
        }
      }
      else
        cpl_msg_error(task, "Not enough memory");

      deleteImage(fluxFile);

    }
    else
      cpl_msg_error(task, "Failure in opening standard star flux table");
  }
  else
    cpl_msg_error(task, "No input standard star flux table found");

  if (error) {
    deleteImage(mosImage);
    deleteImage(imaSpEx1DCal);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteTable(atmTable);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }

  /*
   * Computation of the efficiency curve and of the response curve
   * begins here.
   */

  readDoubleDescriptor(imaSpEx1DCal->descs,
                       pilTrnGetKeyword("Cdelt", 1), &step, NULL);

  readDoubleDescriptor(imaSpEx1DCal->descs,
                       pilTrnGetKeyword("Crval", 1), &start, NULL);

  /*
   * Convert extracted spectrum in electrons per second per Angstrom.
   */

  spectrum = duplicateImage(imaSpEx1DCal);
  constArithLocal(spectrum, gain / time / step, VM_OPER_MUL);

  /*
   * Map the atmospheric extinction factors to the same lambda sampling
   * of the extracted spectrum, and convert to actual flux loss.
   */

  extinction = duplicateImage(imaSpEx1DCal);
  mapTable(extinction, start, step, atmTable, "WAVE", "EXTINCTION");
  deleteTable(atmTable);

  count = 0;
  for (i = 0; i < extinction->xlen; i++) {
    if (extinction->data[i] > 0.0) {
      if (count == 0) {
        pos = i;
      }
      count++;
    }
    else {
      if (count) {
        break;
      }
    }
  }

  constArithLocal(extinction, 0.4 * airmass, VM_OPER_MUL);

  for (i = 0; i < extinction->xlen; i++)
    if (extinction->data[i] > 0.0)
      extinction->data[i] = pow(10., extinction->data[i]);

  /*
   * Correct the scientific spectrum to airmass 0
   */

  imageArithLocal(spectrum, extinction, VM_OPER_MUL);
  deleteImage(extinction);

  /*
   * Map the standard star catalog flux to the same lambda sampling
   * of the extracted spectrum.
   */

  stdFlux = duplicateImage(imaSpEx1DCal);
  mapTable(stdFlux, start, step, fluxTable, "WAVE", "FLUX");
  deleteTable(fluxTable);

  preCount = count;
  prePos = pos;
  count = 0;
  pos = 0;
  for (i = 0; i < stdFlux->xlen; i++) {
    if (stdFlux->data[i] > 0.0) {
      if (count == 0) {
        pos = i;
      }
      count++;
    }
    else {
      if (count) {
        break;
      }
    }
  }

  /*
   * Intersection with previous selection
   */

  newStart = prePos > pos ? prePos : pos;
  newEnd = (prePos + preCount) < (pos + count) ? 
           (prePos + preCount) : (pos + count);
  pos = newStart;
  count = newEnd - newStart;

  /*
   * Convert the flux to photons (per second per Angstrom).
   * stdFlux is in units of erg / cm^2 / s / Angstrom. This
   * must be multiplied by the efficient area of the telescope,
   * 5.18E+5 cm^2, and divided by hv (v = frequency). With 
   * hc = 1.98E-8 erg*Angstrom one obtains the following:
   */

  convStdFlux = duplicateImage(imaSpEx1DCal);

  for (i = 0; i < convStdFlux->xlen; i++) {
    lambda = start + step * i;
    convStdFlux->data[i] = 0.0026 * lambda * stdFlux->data[i];
  }

  efficiency = duplicateImage(imaSpEx1DCal);

  for (i = 0; i < efficiency->xlen; i++) {
    if (convStdFlux->data[i] > 0.0)
      efficiency->data[i] = spectrum->data[i] / convStdFlux->data[i];
    else
      efficiency->data[i] = 0.0;
  }

  deleteImage(convStdFlux);

  preCount = count;
  prePos = pos;
  count = 0;
  pos = 0;
  for (i = 0; i < efficiency->xlen; i++) {
    if (efficiency->data[i] > 0.01) {
      if (count == 0) {
        pos = i;
      }
      count++;
    }
    else {
      if (count > 300) {
        break;
      }
    }
  }

  /*
   * Intersection with previous selection
   */

  newStart = prePos > pos ? prePos : pos;
  newEnd = (prePos + preCount) < (pos + count) ?
           (prePos + preCount) : (pos + count);
  pos = newStart + 2;
  count = newEnd - newStart - 2;

  if (count <= 0) {
    cpl_msg_warning(task, "No calibration is possible with given input data");
    pos = prePos;
    count = preCount;
/*
    deleteImage(efficiency);
    deleteImage(stdFlux);
    deleteImage(spectrum);
    deleteImage(mosImage);
    deleteImage(imaSpEx1DCal);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
*/
  }

  tmpImage = newImageAndAlloc(count, 1);

  for (i = 0; i < count; i++)
    tmpImage->data[i] = efficiency->data[pos + i];

  polySmooth(tmpImage, responseOrder, 50);

  smo_efficiency = duplicateImage(efficiency);

  for (i = 0; i < count; i++)
    smo_efficiency->data[pos + i] = tmpImage->data[i];

  deleteImage(tmpImage);

  response = duplicateImage(imaSpEx1DCal);

  for (i = 0; i < response->xlen; i++) {
    if (efficiency->data[i] > 0.01 && stdFlux->data[i] != 0.0)
      response->data[i] = spectrum->data[i] / stdFlux->data[i];
    else
      response->data[i] = 0.0;
  }

/*   response = imageArith(spectrum, stdFlux, VM_OPER_DIV); */

  tmpImage = newImageAndAlloc(count, 1);

  for (i = 0; i < count; i++)
    tmpImage->data[i] = response->data[pos + i];

  polySmooth(tmpImage, responseOrder, 50);

  smo_response = duplicateImage(response);

  for (i = 0; i < count; i++)
    smo_response->data[pos + i] = tmpImage->data[i];

  deleteImage(tmpImage);

  for (i = 0; i < response->xlen; i++) {
    if (efficiency->data[i] > 0.01) {
      response->data[i] = 1 / response->data[i];
      smo_response->data[i] = 1 / smo_response->data[i];
    }
    else {
      response->data[i] = 0;
      smo_response->data[i] = 0;
    }
  }

  /*
   * Assemble the product spectrophotometric table.
   */

  sphotTable = newSpecPhotTable(response->xlen);

  data = tblGetDoubleData(sphotTable, "WAVE");
  for (i = 0; i < response->xlen; i++)
    data[i] = start + step * i;

  data = tblGetDoubleData(sphotTable, "STD_FLUX");
  for (i = 0; i < response->xlen; i++)
    data[i] = stdFlux->data[i];

  deleteImage(stdFlux);

  data = tblGetDoubleData(sphotTable, "OBS_FLUX");
  for (i = 0; i < response->xlen; i++)
    data[i] = spectrum->data[i];

  deleteImage(spectrum);

  data = tblGetDoubleData(sphotTable, "RAW_EFFICIENCY");
  for (i = 0; i < response->xlen; i++)
    data[i] = efficiency->data[i];

  deleteImage(efficiency);

  data = tblGetDoubleData(sphotTable, "EFFICIENCY");
  for (i = 0; i < response->xlen; i++)
    data[i] = smo_efficiency->data[i];

  deleteImage(smo_efficiency);

  data = tblGetDoubleData(sphotTable, "RAW_RESPONSE");
  for (i = 0; i < response->xlen; i++)
    data[i] = response->data[i];

  deleteImage(response);

  data = tblGetDoubleData(sphotTable, "RESPONSE");
  for (i = 0; i < smo_response->xlen; i++)
    data[i] = smo_response->data[i];

  deleteImage(smo_response);

  specPhotTableHeader(sphotTable, mosImage->descs);

  vimosDscCopy(&sphotTable->descs, mosImage->descs,
               pilTrnGetKeyword("TargetName"), NULL);

  vimosDscCopy(&sphotTable->descs, mosImage->descs,
               pilTrnGetKeyword("MaskId", quadrant), NULL);

  if (computeQC) {

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      cpl_msg_info(task, "Computing QC1 parameters...");

      /*
       * Write PRO.CATG, ARCFILE, TPL ID, and OCS CON QUAD, to QC1 group.
       */

      pilQcWriteString("PRO.CATG", sphotTableCategory, "Product category");

      qcCopyValue(mosImage->descs, pilTrnGetKeyword("ArchiveFile"),
                  NULL, "Archive File Name");

      qcCopyValue(mosImage->descs, pilTrnGetKeyword("TplId"),
                  NULL, "Template signature ID");

      qcCopyValue(mosImage->descs, pilTrnGetKeyword("Quadrant"),
                  NULL, "Quadrant");

      readIntDescriptor(mosImage->descs, pilTrnGetKeyword("Quadrant"),
                        &quadrant, NULL);

      qcCopyValue(mosImage->descs, pilTrnGetKeyword("FilterName", quadrant),
                  NULL, "Filter name");

      qcCopyValue(mosImage->descs, pilTrnGetKeyword("GrismName", quadrant),
                  NULL, "Grism name");

      qcCopyValue(mosImage->descs, pilTrnGetKeyword("MaskId", quadrant),
                  NULL, "Mask identification");

      wstart = 3700.;
      wstep = 400.;
      wcount = 15;
      tmpImage = newImageAndAlloc(wcount, 1);
      mapTableDouble(tmpImage, wstart, wstep, sphotTable, "WAVE", "EFFICIENCY");

      for (i = 0; i < wcount; i++) {
        sprintf(parName, "QC.MOS.EFFICIENCY%d.LAMBDA", i + 1);
        qcWriteValueDouble(sphotTable->descs, wstart + wstep * i, parName, 
                           "Angstrom", "Wavelength of efficiency evaluation");
  
        sprintf(parName, "QC.MOS.EFFICIENCY%d", i + 1);
        qcWriteValueDouble(sphotTable->descs, tmpImage->data[i], parName, 
                           "e-/photon", "Efficiency");
      }

      deleteImage(tmpImage);

      if (pilQcGroupEnd() == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

    }
    else
      cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

  } /* End of QC1 computation. */

  deleteImage(mosImage);


  /*
   * Update the products header
   *
   * Note that for the moment also keywords which are not task specific
   * are handled here, since this is the last possibility to access
   * the linked list of keywords without reopening the file.
   * This may change in future!
   */

  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1DCal->descs),
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(imaSpEx1DCal),
                         "Minimum pixel value", "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1DCal->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(imaSpEx1DCal),
                         "Maximum pixel value", "ESO*", 1);

  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1DCal->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(imaSpEx1DCal),
                         "Mean pixel value");

  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1DCal->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(imaSpEx1DCal),
                         "Standard deviation of pixel values");

  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1DCal->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(imaSpEx1DCal),
                         "Median pixel value");

  updateOK = updateOK && insertStringDescriptor(&(imaSpEx1DCal->descs),
                         pilTrnGetKeyword("DoCategory"), mosReducedCategory,
                         "Category of pipeline product", "ESO*", 1);

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    deleteTable(sphotTable);
    deleteImage(imaSpEx1DCal);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  } 

  deleteSetOfDescriptors(&(imaSpEx1DCal->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(imaSpEx1DCal->descs), "ESO DPR*");

  deleteSetOfDescriptors(&(imaSpEx1DCal->descs), "CD1_*");  
  deleteSetOfDescriptors(&(imaSpEx1DCal->descs), "CD2_*");  
                                              /* to use IRAF splot task */

  removeDescriptor(&(imaSpEx1DCal->descs), 
                   pilTrnGetKeyword("TplExposureNumber"));

  vmstrlower(strcpy(output1DName, mosReducedCategory));
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

  if (createFitsImage(output1DName, imaSpEx1DCal, mosReducedCategory)) {
    outputFrame = newPilFrame(output1DName, mosReducedCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", output1DName);
    deleteTable(sphotTable);
    deleteImage(imaSpEx1DCal);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }
 
  deleteImage(imaSpEx1DCal);


  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyEx[0]->descs),
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(outSpSkyEx[0]),
                         "Minimum pixel value", "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyEx[0]->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(outSpSkyEx[0]),
                         "Maximum pixel value", "ESO*", 1);

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyEx[0]->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(outSpSkyEx[0]),
                         "Mean pixel value");

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyEx[0]->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(outSpSkyEx[0]),
                         "Standard deviation of pixel values");

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyEx[0]->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(outSpSkyEx[0]),
                         "Median pixel value");


  updateOK = updateOK && insertStringDescriptor(&(outSpSkyEx[0]->descs),
                         pilTrnGetKeyword("DoCategory"), mosExtractedCategory,
                         "Category of pipeline product", "ESO*", 1);

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    deleteTable(sphotTable);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }

  deleteSetOfDescriptors(&(outSpSkyEx[0]->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(outSpSkyEx[0]->descs), "ESO DPR*");

  deleteSetOfDescriptors(&(outSpSkyEx[0]->descs), "CD1_*");
  deleteSetOfDescriptors(&(outSpSkyEx[0]->descs), "CD2_*");
                                           /* to use IRAF splot task */

  removeDescriptor(&(outSpSkyEx[0]->descs), 
                   pilTrnGetKeyword("TplExposureNumber"));

  vmstrlower(strcpy(output2DName, mosExtractedCategory));
  strcat(output2DName, ".fits");

/* Align WCS to convention used with IFU */

  readDoubleDescriptor(outSpSkyEx[0]->descs,
                       pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
  readDoubleDescriptor(outSpSkyEx[0]->descs,
                       pilTrnGetKeyword("Crval", 1), &crval, 0);
  crval -= cdelt / 2;
  writeDoubleDescriptor(&(outSpSkyEx[0]->descs),
                        pilTrnGetKeyword("Crval", 1), crval,
                        pilTrnGetComment("Crval"));
/* */

  if (createFitsImage(output2DName, outSpSkyEx[0], mosExtractedCategory)) {
    outputFrame = newPilFrame(output2DName, mosExtractedCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", output2DName);
    deleteTable(sphotTable);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[0]);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }
 
  deleteImage(outSpSkyEx[0]);


  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyEx[1]->descs),
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(outSpSkyEx[1]),
                         "Minimum pixel value", "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyEx[1]->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(outSpSkyEx[1]),
                         "Maximum pixel value", "ESO*", 1);

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyEx[1]->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(outSpSkyEx[1]),
                         "Mean pixel value");

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyEx[1]->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(outSpSkyEx[1]),
                         "Standard deviation of pixel values");

  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyEx[1]->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(outSpSkyEx[1]),
                         "Median pixel value");


  updateOK = updateOK && insertStringDescriptor(&(outSpSkyEx[1]->descs),
                         pilTrnGetKeyword("DoCategory"), mosSkyCategory,
                         "Category of pipeline product", "ESO*", 1);

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    deleteTable(sphotTable);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }

  deleteSetOfDescriptors(&(outSpSkyEx[1]->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(outSpSkyEx[1]->descs), "ESO DPR*");
  deleteSetOfDescriptors(&(outSpSkyEx[1]->descs), "CD1_*");
  deleteSetOfDescriptors(&(outSpSkyEx[1]->descs), "CD2_*");
                                          /* to use IRAF splot task */

  removeDescriptor(&(outSpSkyEx[1]->descs), 
                   pilTrnGetKeyword("TplExposureNumber"));

  vmstrlower(strcpy(output2DSkyName, mosSkyCategory));
  strcat(output2DSkyName, ".fits");

/* Align WCS to convention used with IFU */

  readDoubleDescriptor(outSpSkyEx[1]->descs,
                       pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
  readDoubleDescriptor(outSpSkyEx[1]->descs,
                       pilTrnGetKeyword("Crval", 1), &crval, 0);
  crval -= cdelt / 2;
  writeDoubleDescriptor(&(outSpSkyEx[1]->descs),
                        pilTrnGetKeyword("Crval", 1), crval,
                        pilTrnGetComment("Crval"));
/* */

  if (createFitsImage(output2DSkyName, outSpSkyEx[1], mosSkyCategory)) {
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
    deleteTable(sphotTable);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyEx[1]);
    cpl_free(outSpSkyEx);
    deleteObjectTable(objectTable);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }

  deleteImage(outSpSkyEx[1]);

  cpl_free(outSpSkyEx);

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
    deleteTable(sphotTable);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    deleteWindowTable(windowTable);
    return EXIT_FAILURE;
  }

  writeStringDescriptor(&(windowTable->descs), pilTrnGetKeyword("DoCategory"),
                        windowTableCategory, "Category of pipeline product");

  deleteSetOfDescriptors(&(windowTable->descs), "ESO QC*");

  vmstrlower(strcpy(windowTableName, windowTableCategory));
  /* strcat(windowTableName, ".TFITS"); */
  strcat(windowTableName, ".fits");

  error = 1;

  if ((windowFile = newImage(0, 0, NULL))) {
    if (VM_TRUE == openNewFitsImage(windowTableName, windowFile)) {
      if (VM_TRUE == writeFitsWindowTable(windowTable, windowFile->fptr)) {
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
  deleteWindowTable(windowTable);

  if (error) {
    deleteTable(sphotTable);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    return EXIT_FAILURE;
  }

  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataMin"), 
                         imageMinimum(imaSpEx1D[1]),
                         "Minimum pixel value", "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(imaSpEx1D[1]),
                         "Maximum pixel value", "ESO*", 1);

  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(imaSpEx1D[1]),
                         "Mean pixel value");

  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(imaSpEx1D[1]),
                         "Standard deviation of pixel values");

  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(imaSpEx1D[1]),
                         "Median pixel value");

  updateOK = updateOK && insertStringDescriptor(&(imaSpEx1D[1]->descs),
                         pilTrnGetKeyword("DoCategory"), mosSkyReducedCategory,
                         "Category of pipeline product", "ESO*", 1);

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    deleteTable(sphotTable);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
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
    deleteTable(sphotTable);
    deleteImage(imaSpEx1D[1]);
    cpl_free(imaSpEx1D);
    return EXIT_FAILURE;
  }

  deleteImage(imaSpEx1D[1]);
  cpl_free(imaSpEx1D);

  vmstrlower(strcpy(sphotTableName, sphotTableCategory));
  /* strcat(sphotTableName, ".TFITS"); */
  strcat(sphotTableName, ".fits");

  if (writeFitsSpecPhotTable(sphotTableName, sphotTable) != VM_TRUE) {
    cpl_msg_error(task, "Cannot create local product file %s!", sphotTableName);

    deleteTable(sphotTable);

    return EXIT_FAILURE; 
  }
  else {
    pilFitsHdrCopy(sphotTableName, 0, NULL, ".*-OBS$", 1); 
    pilFitsHdrCopy(sphotTableName, 0, NULL, "^ESO .*", 1);

    outputFrame = newPilFrame(sphotTableName, sphotTableCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);
    pilSofInsert(sof, outputFrame);
  }

  deleteTable(sphotTable);

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
                    "vmmosstandard",
    "Reduce a MOS standard star exposure.",
    "This recipe is used to calibrate and extract a standard star spectrum\n"
    "obtained in MOS mode, and compare it with a catalog spectrum to obtain\n"
    "the efficiency and the response curves.\n\n"
    "Input files:\n\n"
    "  DO category:               Type:      Explanation:         Required:\n"
    "  MOS_STANDARD               Raw        Science exposure        Y\n"
    "  MASTER_BIAS                Calib      Master bias             Y\n"
    "  EXTINCT_TABLE              Calib      Atmospheric extinction  Y\n"
    "  STD_FLUX_TABLE             Calib      Standard star fluxes    Y\n"
    "  MASTER_DARK                Calib      Master dark             .\n"
    "  MOS_MASTER_SCREEN_FLAT     Calib      Normalised flat field   .\n"
    "  EXTRACT_TABLE              Calib      Extraction table        .\n"
    "  GRISM_TABLE                Calib      Grism table             Y\n"
    "  CCD_TABLE                  Calib      Bad pixel table         .\n\n"
    "Output files:\n\n"
    "  DO category:               Data type: Explanation:\n"
    "  MOS_STANDARD_REDUCED       FITS image Extracted spectrum\n"
    "  MOS_STANDARD_EXTRACTED     FITS image Sky subtracted slit spectra\n"
    "  MOS_STANDARD_SKY           FITS image Sky slit spectra\n"
    "  MOS_STANDARD_SKY_EXTRACTED FITS image Extracted sky spectrum\n"
    "  OBJECT_TABLE               FITS table Objects spectra identification\n"
    "  WINDOW_TABLE               FITS table Objects positions in slit\n"
    "  MOS_SPECPHOT_TABLE         FITS table Efficiency and response curves\n\n"
    "The standard star spectrum is calibrated and extracted as described\n"
    "for the recipe vmmosobsstare. The extracted spectrum is compared to\n"
    "catalog fluxes of the same star to obtain the efficiency curve and\n"
    "the response curve to use for flux calibration. For this comparison\n"
    "the extracted spectrum is converted into electrons/sec/Angstrom. The\n"
    "catalog fluxes, as well as the atmospheric extinction factors, are\n"
    "remapped to the same wavelength bins of the extracted spectrum. After\n"
    "removing the effect of the atmospheric extinction (i.e., correcting\n"
    "the extracted spectrum to airmass zero), the response curve is derived\n"
    "as the ratio between the catalog spectrum and the calibrated spectrum.\n"
    "The efficiency curve is obtained by dividing the extracted spectrum\n" 
    "in electrons/sec/Angstrom by the the catalog fluxes converted from\n"
    "erg/cm/cm/s/Angstrom to photons/s/Angstrom (at 100% efficiency each\n" 
    "photon would correspond to one electron).\n\n"
    "For more details, please refer to the VIMOS Pipeline User's Guide.",

    "ESO VIMOS Pipeline Team",

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

                    vmmosstandard_create,
                    vmmosstandard_exec,
                    vmmosstandard_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
