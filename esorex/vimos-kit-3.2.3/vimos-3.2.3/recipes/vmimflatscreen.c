/* $Id: vmimflatscreen.c,v 1.4 2013-07-11 11:49:24 cgarcia Exp $
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
 * $Date: 2013-07-11 11:49:24 $
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
#include <pilerrno.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilrecipe.h>
#include <pilqc.h>
#include <pilutils.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmqcutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


#define MAX_COMMENT_LENGTH  (80)
#define SATURATION_LEVEL    (65535)
#define QC1_NSTEP           (10)

static cxint vmimflatscreen(PilSetOfFrames *);


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

static const char *filterMethodNames[] = {
  "Auto",
  "Median",
  "Average"
};

static const FilterMethod filterMethods[] = {
  FILTER_AUTO,
  FILTER_MEDIAN,
  FILTER_AVERAGE
};

static unsigned int nFilterMethods = sizeof(filterMethods)
                                   / sizeof(FilterMethod);

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
 * recipe and make it available to the application using the interface.
 */

static cxint
vmimflatscreen_create(cpl_plugin *plugin)
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

    p = cpl_parameter_new_value("vimos.Parameters.smooth.boxsize",
                                CPL_TYPE_INT,
                                "Size of running box used for smoothing",
                                "vimos.Parameters",
                                7);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SmoothBoxSize");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SmoothBoxSize");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.smooth.method",
                               CPL_TYPE_STRING,
                               "Smooth method for image trends removal",
                               "vimos.Parameters",
                               "Median", 3, "Median", "Average", "Auto");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SmoothMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SmoothMethod");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.stacking.singleframes",
                                CPL_TYPE_BOOL,
                                "Frame combination method is ignored.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "AllowSingleFrames");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "AllowSingleFrames");
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
                                "Low threshold for K-sigma clipping method",
                                "vimos.Parameters",
                                5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "KSigmaLow");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "KSigmaLow");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.stacking.ksigma.high",
                                CPL_TYPE_DOUBLE,
                                "High threshold for K-sigma clipping method",
                                "vimos.Parameters",
                                5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "KSigmaHigh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "KSigmaHigh");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.minmax.minimum",
                                CPL_TYPE_INT,
                                "Number of lowest rejected values for "
                                "rejection method",
                                "vimos.Parameters",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MinRejection");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MinRejection");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.minmax.maximum",
                                CPL_TYPE_INT,
                                "Number of highest rejected values for "
                                "rejection method",
                                "vimos.Parameters",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MaxRejection");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MaxRejection");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.flat.frames.validate",
                                CPL_TYPE_BOOL,
                                "Consistency check on input frames",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ValidateFrames");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ValidateFrames");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.tolerance.low",
                                CPL_TYPE_DOUBLE,
                                "Minimum level for frame acceptance",
                                "vimos.Parameters",
                                1000.);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LevelToleranceLow");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LevelToleranceLow");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.tolerance.high",
                                CPL_TYPE_DOUBLE,
                                "Maximum level for frame acceptance",
                                "vimos.Parameters",
                                60000.);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LevelToleranceHigh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LevelToleranceHigh");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.tolerance.pattern",
                                CPL_TYPE_DOUBLE,
                                "Threshold for flux pattern consistency",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PatternTolerance");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PatternTolerance");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

    p = cpl_parameter_new_enum("vimos.Parameters.bias.removing.method",
                               CPL_TYPE_STRING,
                               "Bias removal method",
                               "vimos.Parameters",
                               "Zmaster", 2, "Zmaster", "Master");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "BiasMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "BiasMethod");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.badpixel.clean",
                                CPL_TYPE_BOOL,
                                "Bad pixel correction on master flat",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.clean",
                                CPL_TYPE_BOOL,
                                "Cosmic ray events cleaning from raw flats",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanCosmic");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanCosmic");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.threshold",
                                CPL_TYPE_DOUBLE,
                                "Threshold to identify cosmic rays "
                                "candidates.",
                                "vimos.Parameters",
                                4.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CosmicThreshold");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CosmicThreshold");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.ratio",
                                CPL_TYPE_DOUBLE,
                                "Min ratio between peak and nearby pixels "
                                "for a cosmic ray.",
                                "vimos.Parameters",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CosmicRatio");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CosmicRatio");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.quality.enable",
                                CPL_TYPE_BOOL,
                                "Compute QC1 parameters.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ComputeQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ComputeQC");
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
vmimflatscreen_exec(cpl_plugin *plugin)
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

    if (vmimflatscreen(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmimflatscreen");
        
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
vmimflatscreen_destroy(cpl_plugin *plugin)
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
 *   Create master screen flat field from a set of raw screen flat fields.
 * 
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param sof  Set of frames containing the references to raw screen
 *             flat field frames, a master bias, an optional master 
 *             dark, and an optional CCD table.
 *
 * @doc 
 *   The recipe function creates a master screen flat field frame from 
 *   a set of raw screen flat field frames, after bias removal and 
 *   optional dark subtraction. The raw frames, the master bias frame,
 *   and the optional master dark frame, are passed via the set of 
 *   frames \textbf{sof}. A CCD table may also be specified, to be 
 *   used in bad pixel correction and/or quality check of the output 
 *   master screen flat field. If requested, input frames consistency 
 *   is checked, and cosmic ray events are removed from each input frame 
 *   before combination. The result is divided by its smoothed image,
 *   to eliminate large scale trends and to normalize it. On successful 
 *   termination the created master screen flat field is added to the 
 *   set of frames.
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
 *
 *   \begin{itemize}
 *
 *     \item SmoothMethod        Method for smoothing the combined input 
 *                               flat fields, for determining (and then 
 *                               removing) the large scale trends related 
 *                               either to the CCD intrinsic response 
 *                               variation, or to a non uniform illumination 
 *                               of the screen. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Median:             Median filter
 *
 *       \item Average:            Average filter
 *
 *     \end{itemize}
 *
 *     \item SmoothBoxSize:      Size of the running box used for smoothing.
 *
 *     \item StackMethod:        Method for raw screen flat field frames 
 *                               combination. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Average:            Average of frames
 *
 *       \item Median:             Median stacking
 *
 *       \item MinMax:             Minimum-Maximum rejection. If this
 *                                 option is chosen, the upper and lower
 *                                 percentiles of pixel values to be
 *                                 rejected before averaging may be 
 *                                 specified using the parameters 
 *                                 MinRejection and MaxRejection.
 *                                 The default is to reject just the 
 *                                 highest and the lowest pixel value.
 *
 *       \item Ksigma:             Kappa-Sigma clipping. If this option
 *                                 is chosen, the upper and lower rejection
 *                                 thresholds, given in units of sigma, may
 *                                 be specified using the parameters KSigmaLow
 *                                 and KSigmaHigh. Pixel values deviating
 *                                 from the median more than the specified
 *                                 thresholds are rejected before final
 *                                 averaging.
 *
 *       \item Auto:               Automatic method (not yet available, it
 *                                 falls back to Average).
 *
 *     \end{itemize}
 *
 *     \item BiasMethod:         Method for bias removal from raw screen flat
 *                               field frames. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Master:             Prescan and overscan regions are trimmed
 *                                 away from raw screen flat fields after 
 *                                 master bias subtraction.
 *
 *       \item Zmaster:            After master bias subtraction the residual
 *                                 signal found in each raw screen flat field
 *                                 overscan regions is modelled and subtracted 
 *                                 from the image. Next, prescan and overscan 
 *                                 regions are trimmed away.
 *
 *     \end{itemize}
 *
 *     \item KSigmaLow:          Lower rejection threshold, given in units
 *                               of sigma, used when parameter StackMethod
 *                               is set to Ksigma.
 *
 *     \item KSigmaHigh:         Upper rejection threshold, given in units
 *                               of sigma, used when parameter StackMethod
 *                               is set to Ksigma.
 *
 *     \item MinRejection:       Lower percent of sorted pixel values to
 *                               be rejected, when parameter StackMethod
 *                               is set to MinMax. If set to 0.0, just
 *                               one pixel value, the lowest, is rejected.
 *
 *     \item MaxRejection:       Upper percent of sorted pixel values to
 *                               be rejected, when parameter StackMethod
 *                               is set to MinMax. If set to 0.0, just
 *                               one pixel value, the highest, is rejected.
 *
 *     \item ValidateFrames:     Rejection of raw flat field frames having
 *                               too low or too high level, or having too
 *                               much different intensity pattern.
 *
 *     \item LevelToleranceLow:  Minimum mean level for an input frame 
 *                               to be accepted, in ADU. This level is 
 *                               meant to avoid frames having a poor
 *                               signal. It can be disabled by setting
 *                               it to a negative value.
 *
 *     \item LevelToleranceHigh: Maximum mean level for an input frame 
 *                               to be accepted, in ADU. This level is
 *                               meant to avoid frames too close to the
 *                               saturation level (65535), where detector
 *                               response linearity tends to be lost. It 
 *                               can be disabled by setting it greater 
 *                               than the saturation level.
 *
 *     \item PatternTolerance:   Threshold controlling frames rejection
 *                               based on pattern differences, in units of
 *                               sigmas.
 *
 *     \item CleanBadPixel:      Bad pixel correction on output master screen
 *                               flat field. If this option is turned on, the 
 *                               recipe expects to find a CCD\_TABLE in the
 *                               input set of frames.
 *
 *     \item CleanCosmic:        Cosmic ray events removal from input raw 
 *                               screen flat field frames. If a CCD\_TABLE 
 *                               is found among the input set of frames, 
 *                               bad pixels will not be used in computing 
 *                               the values to replace the cosmic rays events.
 *
 *     \item CosmicThreshold:    Threshold for the detection of cosmic rays,
 *                               having the same meaning of the parameter ns
 *                               of the MIDAS command FILTER/COSMIC: it is
 *                               the number of theoretical noise sigmas above
 *                               smoothed level that must be reached to make
 *                               a pixel a cosmic ray event location candidate.
 *                               This parameter is effective when CleanCosmic
 *                               is set to "true".
 *
 *     \item CosmicRatio:        Critical ratio for discrimination of objects
 *                               and cosmic rays, having the same meaning
 *                               of the parameter rc of the MIDAS command
 *                               FILTER/COSMIC: it is the ratio between the
 *                               peak of a cosmic ray event candidate, and the
 *                               average of the 8 nearby pixels that must be
 *                               reached to identify the candidate as a real
 *                               cosmic ray event.
 *
 *     \item ComputeQC:          Compute QC1 parameters and write them to
 *                               the QC1 PAF file and to the primary product
 *                               header.
 *
 *   \end{itemize}
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 *
 *   The master screen flat field contains just information on the small 
 *   scale variations (i.e. less than SmoothBoxSize pixel) of the CCD 
 *   response, and it should not be used in the reduction of scientific 
 *   data.
 *   
 * @author C. Izzo, P. Sartoretti, R. Palsa
 */   

static cxint 
vmimflatscreen(PilSetOfFrames *sof)
{

  const char     task[] = "vmimflatscreen";

  const char     parameter[] = "Parameters";

  const char    *flatTag  = pilTrnGetCategory("ImgScreenFlat");
  const char    *mFlatTag = pilTrnGetCategory("ImgMasterScreenFlat");
  const char    *cFlatTag = pilTrnGetCategory("ImgCombScreenFlat");

  char          *combMethodTag = NULL;
  char          *biasMethodTag = NULL;
  char          *filterMethodTag = NULL;
  char           mFlatName[PATHNAME_MAX + 1];
  char           cFlatName[PATHNAME_MAX + 1]; 
  char           comment[MAX_COMMENT_LENGTH];

  VimosBool      updateOK = VM_TRUE;

  size_t         flatCount, goodFrames, minFrames;

  unsigned int   i, j;
  unsigned int   cleanBadPixel, cleanCosmic, validateFrames;
  unsigned int   computeQC, singleFrames;
  unsigned int   error;

  int            biasMethodEntry;
  int            combMethodEntry;
  int            filterMethodEntry;

  int            quadrant;

  float         *flatLevel = NULL;
  float         *flatNoise = NULL;
  float         *flatRon = NULL;
  float         *flatGain = NULL;
  float         *exposureTime = NULL;

  float          expTolerance;
  float          relExpDiff;

  double         time, gain;
  double         sumExposures = 0.0;
  double         rescaleFactor;
  double         lowerThreshold, upperThreshold;
  double         patternTolerance;
  double         cosmicThreshold = 0.0;
  double         cosmicRatio = 0.0;

  PilFrame      *currFrame, *ccdFrame, *biasFrame, *darkFrame, *outputFrame;

  BiasMethod     biasMethod = BIAS_UNDEF;

  CombMethod     combMethod = COMB_UNDEF;
  CombParameters combParameter;

  CombMethod     filterMethod = FILTER_UNDEF;
  int            smoothBoxSize;

  VimosImage    *mFlat = NULL;
  VimosImage   **flatList, **original;
  VimosImage    *biasImage = NULL;
  VimosImage    *darkImage = NULL;
  VimosImage    *biasOver;
  VimosImage    *mFlatNorm, *mFlatFilt, *mFlatClean;

  VimosTable    *ccdTable = NULL;


  /*
   * Get task parameters from the recipe database
   */

  /*
   * Determine the frame bias removal method first.
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
    minFrames = MIN_FRAMES_REJECT;
    combParameter.minRejection = pilDfsDbGetInt(parameter, "MinRejection", 1);
    combParameter.maxRejection = pilDfsDbGetInt(parameter, "MaxRejection", 1);
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
    combMethod = COMB_AVERAGE;
    minFrames = MIN_FRAMES_AVERAGE;
    break;
  }


  /*
   * Check if a single frame in input should be tolerated. This
   * means that if a single frame is found in input, then the
   * stacking method will be ignored.
   */

  singleFrames = pilDfsDbGetBool("Parameters", "AllowSingleFrames", 0);

  /*
   * Determine the frame smoothing method and all method dependent
   * parameters.
   */

  filterMethodTag = (char *)pilDfsDbGetString(parameter, "SmoothMethod");

  if ((filterMethodEntry =
        strselect(filterMethodTag, filterMethodNames, nFilterMethods)) < 0) {
    cpl_msg_error(task, "%s: Invalid frame smoothing method.", filterMethodTag);
    return EXIT_FAILURE;
  }

  filterMethod = filterMethods[filterMethodEntry];

  smoothBoxSize = pilDfsDbGetInt(parameter, "SmoothBoxSize", 11);

  if (smoothBoxSize < 3) {
    cpl_msg_error(task, "Invalid smooth box size, %d: it should be at least 3", 
                smoothBoxSize);
    return EXIT_FAILURE;
  }

  /*
   * Check if raw frames with non-nominal exposure level should be ignored
   * by the image stacking task. Default is to remove non-nominal frames.
   */

  validateFrames = pilDfsDbGetBool(parameter, "ValidateFrames", 0);

  if (validateFrames) {
    patternTolerance = pilDfsDbGetDouble(parameter, "PatternTolerance", 3.0);
    lowerThreshold = pilDfsDbGetDouble(parameter, 
                                       "LevelToleranceLow", 1000.);
    upperThreshold = pilDfsDbGetDouble(parameter, 
                                       "LevelToleranceHigh", 50000.);

    if (patternTolerance < 1.0) {
      cpl_msg_error(task, "Invalid tolerance, it should be at least one sigma.");
      return EXIT_FAILURE;
    }

  }


  /*
   * Check if the bad pixels should be corrected.
   */

  cleanBadPixel = pilDfsDbGetBool(parameter, "CleanBadPixel", 1);


  /*
   * Check if the cosmic rays events should be corrected.
   */

  cleanCosmic = pilDfsDbGetBool(parameter, "CleanCosmic", 1);

  if (cleanCosmic) {
    cosmicThreshold = pilDfsDbGetDouble(parameter, "CosmicThreshold", 4.0);
    cosmicRatio = pilDfsDbGetDouble(parameter, "CosmicRatio", 2.0);
    if (cosmicThreshold < 1.0 || cosmicRatio < 1.0) {
      cpl_msg_error(task, "Invalid cosmic ray filtering parameters: both "
                  "CosmicThreshold and CosmicRatio must be greater than one");
      return EXIT_FAILURE;
    }
  }


  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   * Make sure that there are enough raw data frames in the
   * input set.
   */

  flatCount = (int)pilSofFrameCount(sof, flatTag);

  if (flatCount < minFrames) {
    if (flatCount == 1 && singleFrames) {
      combMethod = COMB_UNDEF;
    }
    else {
      cpl_msg_error(task, "Not enough raw screen flat field frames in "
                  "input for stacking method '%s'!", combMethodTag);
      return EXIT_FAILURE;
    }
  }

  if (validateFrames && flatCount < 2) {
    validateFrames = 0;
    cpl_msg_warning(task, "Too few screen flat field frames (%zd) in input. "
                  "Skipping frame validation task!", flatCount);
  }

  cpl_msg_info(task, "Loading input frames...");

  /*
   * If bad pixel correction is enabled, a bad pixel table is required 
   * in the input set of frames. If no bad pixel table is present this 
   * is an error. If cleaning of cosmic rays is enabled, the bad pixel 
   * table will be used, if present, to avoid bad pixels while smoothing 
   * away cosmic ray events. In this case, a missing bad pixel table is 
   * not an error.
   */

  if ((ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable")))) {
    pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);
    cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));
    if ((ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
      closeFitsTable(ccdTable, 0);
    }
    else {
      cpl_msg_error(task, "Failure in opening CCD table");
      return EXIT_FAILURE;
    }
  }
  else {
    if (cleanBadPixel) {
      cpl_msg_error(task, "No CCD table in input: cannot clean bad pixels.");
      return EXIT_FAILURE;
    }
    cpl_msg_warning(task, "No CCD table in input: flat field normalization "
                  "may introduce ghosts from bad pixels.");
  }


 /*
  * Get the master bias frame
  */

  error = 1;

  if ((biasFrame = pilSofLookup(sof, pilTrnGetCategory("MasterBias")))) {
    pilFrmSetType(biasFrame, PIL_FRAME_TYPE_CALIB);
    if ((biasImage = openOldFitsFile(pilFrmGetName(biasFrame), 1, 0))) {
      closeFitsImage(biasImage, 0);
      error = 0;
    }
    else 
      cpl_msg_error(task, "Failure opening master bias frame");
  }
  else 
    cpl_msg_error(task, "No master bias in input");

  if (error) {
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Get the (optional) master dark frame
  */

  if ((darkFrame = pilSofLookup(sof, pilTrnGetCategory("MasterDark")))) {
    pilFrmSetType(darkFrame, PIL_FRAME_TYPE_CALIB);
    if ((darkImage = openOldFitsFile(pilFrmGetName(darkFrame), 1, 0))) {
      closeFitsImage(darkImage, 0);
    }
    else {
      cpl_msg_error(task, "Failure opening master dark frame");
      deleteImage(biasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }
  else 
    cpl_msg_warning(task, "No master dark in input, dark subtraction "
                  "will not be performed");


 /*
  * Load the raw screen flat field frames.
  */

  error = 1;

  flatList = (VimosImage **)cpl_calloc(flatCount, sizeof(VimosImage *));
  if (flatList) {

    error = 0;
    currFrame = pilSofLookupNext(sof, flatTag);

    for (i = 0; i < flatCount; i++) {
      if ((flatList[i] = openOldFitsFile(pilFrmGetName(currFrame), 1, 0))) {
        pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
        closeFitsImage(flatList[i], 0);
      }
      else {
        error = 1;
        cpl_msg_error(task, 
                    "Failure opening screen flat field frame %d", i + 1);
        for (j = 0; j < i; j++) 
          deleteImage(flatList[j]);
        cpl_free(flatList);
        break;
      }
      currFrame = pilSofLookupNext(sof, NULL);
    }
  }
  else
    cpl_msg_error(task, 
                "Failure in creating list of input raw screen flat fields");

  if (error) {
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteCcdTable(ccdTable); 
    return EXIT_FAILURE;
  }

 /*
  * Recreate bias overscans using as reference the first input flat field
  */

  if ((biasOver = growOverscans(biasImage, flatList[0]))) {
    if (biasImage != biasOver) {
      deleteImage(biasImage);
      biasImage = biasOver;
    }
  }
  else {
    cpl_msg_error(task, "Failure in growing overscans in master bias");
    for (i = 0; i < flatCount; i++) 
      deleteImage(flatList[i]);
    cpl_free(flatList);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Get exposure times
  */

  error = 1;

  if ((exposureTime = (float *)cpl_calloc(flatCount, sizeof(float)))) {
    error = 0;
    for (i = 0; i < flatCount; i++) {
      if ((readDoubleDescriptor(flatList[i]->descs,
                                pilTrnGetKeyword("ExposureTime"), 
                                &time, comment)) == VM_TRUE) {
        if (time < MIN_DIVISOR) {
          error = 1;
          cpl_msg_error(task, 
            "Zero or negative exposure time for flat field %d", i + 1);
          cpl_free(exposureTime);
          break;
        }
        else {
          exposureTime[i] = time;
        }
      }
      else {
        error = 1;
        cpl_msg_error(task, "Cannot read %s from flat field frame %d", 
                    pilTrnGetKeyword("ExposureTime"), i + 1);
        cpl_free(exposureTime);
        break;
      }
    }
  }
  else 
    cpl_msg_error(task, "Problems with memory allocation");

  if (error) {
    for (j = 0; j < flatCount; j++)
      deleteImage(flatList[j]);
    cpl_free(flatList);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Get gain factors (e-/ADU)
  */

  error = 1;

  if ((flatGain = (float *)cpl_calloc(flatCount, sizeof(float)))) {
    error = 0;
    for (i = 0; i < flatCount; i++) {
      gain = getMeanGainFactor(flatList[i]);
      if (gain > MIN_DIVISOR) {
        flatGain[i] = gain;
      }
      else {
        error = 1;
        cpl_msg_error(task, "Wrong or no gain factor for flat field %d", i + 1);
        cpl_free(flatGain);
        break;
      }
    }
  }
  else 
    cpl_msg_error(task, "Problems with memory allocation");

  if (error) {
    for (j = 0; j < flatCount; j++)
      deleteImage(flatList[j]);
    cpl_free(flatList);
    cpl_free(exposureTime);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Noise contributions for each input screen flat field frame must be 
  * estimated before any processing (as bias subtraction, normalization, 
  * etc.) is applied.
  */

  if ((flatRon = (float *)cpl_calloc(flatCount, sizeof(float)))) {
    for (i = 0; i < flatCount; i++) {
      flatRon[i] = computeAverageRon(flatList[i]);
    }
  }
  else {
    cpl_msg_error(task, "Problems with memory allocation");

    for (i = 0; i < flatCount; i++)
      deleteImage(flatList[i]);
    cpl_free(flatList);
    cpl_free(exposureTime);
    cpl_free(flatGain);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  if ((flatNoise = (float *)cpl_calloc(flatCount, sizeof(float)))) {
    for (i = 0; i < flatCount; i++) {
      flatNoise[i] = evaluateAverageNoise(flatList[i], 
                     flatRon[i], flatGain[i]);
    }
  }
  else {
    cpl_msg_error(task, "Problems with memory allocation");

    for (i = 0; i < flatCount; i++)
      deleteImage(flatList[i]);
    cpl_free(flatList);
    cpl_free(exposureTime);
    cpl_free(flatGain);
    cpl_free(flatRon);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteCcdTable(ccdTable);
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

  for (i = 0; i < flatCount; i++) {
    if (VmSubBias(flatList[i], biasImage, biasMethod) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot remove bias from input flat field %d", i + 1);
      for (i = 0; i < flatCount; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);
      cpl_free(exposureTime);
      cpl_free(flatGain);
      cpl_free(flatRon);
      cpl_free(flatNoise);
      deleteImage(biasImage);
      deleteImage(darkImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  deleteImage(biasImage);


 /*
  * Now go for the (optional) dark subtraction
  */

  if (darkImage) {
    cpl_msg_info(task, "Dark subtraction...");

    for (i = 0; i < flatCount; i++) {
      if (VmSubDark(flatList[i], darkImage) == EXIT_FAILURE) {
        cpl_msg_error(task, 
                    "Cannot subtract dark from input flat field %d", i + 1);
        for (i = 0; i < flatCount; i++)
          deleteImage(flatList[i]);
        cpl_free(flatList);
        cpl_free(exposureTime);
        cpl_free(flatGain);
        cpl_free(flatRon);
        cpl_free(flatNoise);
        deleteImage(darkImage);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }
    }
    deleteImage(darkImage);
  }


 /*
  * Now flat fields levels can be computed
  */

  if ((flatLevel = (float *)cpl_calloc(flatCount, sizeof(float)))) {
    for (i = 0; i < flatCount; i++) {
      flatLevel[i] = imageMedian(flatList[i]);
      cpl_msg_info(task, "Level of flat field %-d is %10.4f +/- %-.4f ADU",
                                      i + 1, flatLevel[i], flatNoise[i]);
    }
  }
  else {
    cpl_msg_error(task, "Problems with memory allocation");
    for (i = 0; i < flatCount; i++)
      deleteImage(flatList[i]);
    cpl_free(flatList);
    cpl_free(flatGain);
    cpl_free(flatRon);
    cpl_free(flatNoise);
    cpl_free(exposureTime);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * When RON and flat fields levels are known, but before any image
  * rescaling is applied, cosmic ray events might be removed from all 
  * input images.
  */

  if (cleanCosmic) {
    cpl_msg_info(task, "Cleaning cosmic ray events...");
    for (i = 0; i < flatCount; i++) {
      if (VmCosmicClean(flatList[i], ccdTable, 0, flatLevel[i], flatGain[i],
                        flatRon[i], cosmicThreshold, cosmicRatio) 
                        == EXIT_FAILURE) {
        cpl_msg_error(task, "Cannot clean cosmic rays from raw screen "
                    "flat field frame %d", i + 1);
        for (i = 0; i < flatCount; i++)
          deleteImage(flatList[i]);
        cpl_free(flatList);
        cpl_free(flatNoise);
        cpl_free(flatGain);
        cpl_free(flatRon);
        cpl_free(flatLevel);
        cpl_free(exposureTime);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }
      cpl_msg_info(task, "  Flat field %d of %zd done", i + 1, flatCount);
    }
  }

  cpl_free(flatRon);
  cpl_free(flatGain);

  if (validateFrames) {

   /*
    * If requested, exclude frames with too low and too high median
    * flux, and with illumination pattern not consistent with the set.
    */

   /*
    * Remove from the input list all flat field frames which do not have
    * flux within a reasonable range. The lowest and highest values are
    * user-defined.
    */

    cpl_msg_info(task, "Selecting good flat field frames ...");

    original = (VimosImage **)cpl_calloc(flatCount, sizeof(VimosImage *));
    if (original) {
      for (i = 0; i < flatCount; i++)
        original[i] = flatList[i];     /* Hold original image sorting */
    }

    if (lowerThreshold > MIN_DIVISOR || upperThreshold < SATURATION_LEVEL) {
      cpl_msg_info(task, "Valid range for mean screen flat field flux: "
                 "%4.0f to %5.0f ADU", lowerThreshold, upperThreshold);
  
     /*
      * Good images are moved to the beginning of the flatList:
      */

      goodFrames = applyListSelection(flatList, flatLevel, flatCount,
                                      lowerThreshold, upperThreshold, 1);
    }
    else 
      goodFrames = flatCount;


    if (goodFrames < flatCount) {

     /*
      * If some images were filtered, then associated array of values
      * must be resorted accordingly
      */

      error = 1;
      if (remapFloatsLikeImages(original, flatList, flatNoise, flatCount)
                                                          == EXIT_SUCCESS) {
        if (remapFloatsLikeImages(original, flatList, exposureTime, flatCount)
                                                          == EXIT_SUCCESS) {
          error = 0;
        }
      }

      if (error) {
        cpl_msg_error(task, "Problems remapping arrays after frame selection");
        for (i = 0; i < flatCount; i++)
          deleteImage(flatList[i]);
        cpl_free(flatList);
        cpl_free(original);
        cpl_free(flatNoise);
        cpl_free(flatLevel);
        cpl_free(exposureTime);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }

    }


   /*
    * Check if the remaining input frames have the same intensity
    * distribution, i.e., if the input flat field frames show a similar
    * illumination pattern. This makes sense only if there is more
    * than one image left.
    */

    if (goodFrames > 1) {
      cpl_msg_info(task, "Checking for consistent signal distribution ...");

     /*
      * All flat fields values are normalized to a common level
      * (arbitrarily, the one of the first flat field image).
      */
    
      for (i = 1; i < goodFrames; i++) {
        rescaleFactor = flatLevel[0] / flatLevel[i];
        constArithLocal(flatList[i], rescaleFactor, VM_OPER_MUL);
        flatNoise[i] *= rescaleFactor;
      }

      for (i = 0; i < flatCount; i++)
        original[i] = flatList[i];     /* Hold original image sorting */

      goodFrames = qcSelectConsistentImages(flatList, flatNoise, goodFrames,
                                            patternTolerance);

      if (goodFrames == 0 && pilErrno == P_EGENERIC) {
	cpl_msg_error(task, "Selection of consistent images failed!");
	return EXIT_FAILURE;
      }

      if (remapFloatsLikeImages(original, flatList, exposureTime, flatCount)
                                                          == EXIT_FAILURE) {
        cpl_msg_error(task, "Problems remapping arrays after frame selection");
        for (i = 0; i < flatCount; i++)
          deleteImage(flatList[i]);
        cpl_free(flatList);
        cpl_free(original);
        cpl_free(flatNoise);
        cpl_free(flatLevel);
        cpl_free(exposureTime);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }
    }

    if (goodFrames) {
      cpl_msg_info(task, "%zd out of %zd consistent screen flat field frames "
                 "selected from input set", goodFrames, flatCount);
    }

    cpl_free(original);

  }
  else
    goodFrames = flatCount;

 /*
  * Destroy inconsistent frames from the rearranged list of flat field
  * frames, free arrays.
  */

  for (i = goodFrames; i < flatCount; i++)
    deleteImage(flatList[i]);
  cpl_free(flatNoise);

  if (!goodFrames) {
    cpl_msg_error(task, "No screen flat field frames with consistent signal "
                "distribution found in input. All frames excluded!");
    cpl_free(flatList);
    cpl_free(exposureTime);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Stack the surviving flat field frames, if requested, and if possible.
  */

  if (combMethod == COMB_UNDEF) {
    mFlat = duplicateImage(flatList[0]);
  }
  else if (goodFrames >= minFrames) {
    cpl_msg_info(task, "Combining %zd frames with method '%s'", goodFrames,
               combMethodNames[combMethodEntry]);

    if (combMethod == COMB_AVERAGE) {

     /* 
      * Average is weighted on the original levels - i.e., eliminate
      * the rescaling that was done in case validateFrames was "true".
      */ 

      if (validateFrames) {
        for (i = 1; i < goodFrames; i++) {
          rescaleFactor = flatLevel[i] / flatLevel[0];
          constArithLocal(flatList[i], rescaleFactor, VM_OPER_MUL);
        }
      }
      mFlat = frCombAverage(flatList, goodFrames);

    }
    else {
      mFlat = frComb(flatList, goodFrames, combMethod, &combParameter, 0);
    }

  }
  else {
    cpl_msg_error(task, "Not enough good frames for combine method "
                "'%s'!", combMethodNames[combMethodEntry]);

    for (i = 0; i < goodFrames; i++)
      deleteImage(flatList[i]);
    cpl_free(flatList);
    deleteCcdTable(ccdTable);
    cpl_free(exposureTime);
    cpl_free(flatLevel);
    return EXIT_FAILURE;
  }

  cpl_free(flatLevel);

  if (!mFlat) {
    cpl_msg_error(task, "Stacking of raw screen flat field frames failure!");
    deleteCcdTable(ccdTable);
    cpl_free(exposureTime);
    for (i = 0; i < goodFrames; i++)
      deleteImage(flatList[i]);
    cpl_free(flatList);
    return EXIT_FAILURE;
  }


  /*
   * Create master screen flat field header
   */

  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               pilTrnGetKeyword("MjdObs"), NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               pilTrnGetKeyword("DateObs"), NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               pilTrnGetKeyword("ArchiveFile"), NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^C[A-Z]*[1,2]", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO OBS (DID|PROG ID)", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO OBS ID", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO TPL [.DINPSV.]", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO INS (DID|MODE)", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO INS FILT[1-4] .*", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO INS LAMP[1-5] .*", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO DET [A-Z].*", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO OCS (DID|CON QUAD)", NULL);

  if (computeQC) {

    if (goodFrames > 1) {

      expTolerance = 0.04;  /* 4% variation to consider 2 exposures equal */
      relExpDiff = fabs(exposureTime[0] - exposureTime[1]) / exposureTime[0];
      
      if (relExpDiff < expTolerance) {

        /*
         *  Currently the central 1600x1800 pixels are used.
         */

        int          winSizeX = 1600;
        int          winSizeY = 1800;
        int          winStepX = 1600 / QC1_NSTEP;
        int          winStepY = 1800 / QC1_NSTEP;
        int          winStartX = (flatList[0]->xlen - winSizeX) / 2;
        int          winStartY = (flatList[0]->ylen - winSizeY) / 2;
        unsigned int npix = winSizeX * winSizeY;
        unsigned int ngains = (winSizeX * winSizeY) / (winStepX * winStepY);

        double       value, meanValue, phnValue, fpnValue, structValue, sigma;

        float       *sample1 = extractFloatImage(flatList[0]->data,
                               flatList[0]->xlen, flatList[0]->ylen,
                               winStartX, winStartY, winSizeX, winSizeY);
        float       *sample2 = extractFloatImage(flatList[1]->data,
                               flatList[1]->xlen, flatList[1]->ylen,
                               winStartX, winStartY, winSizeX, winSizeY);
        float       *diff    = cpl_malloc(winSizeX * winSizeY * sizeof(float));
        float       *convFactors;
        float       *box;
        float       *diffBox;

        int          x, y;


        cpl_msg_info(task, "Computing QC1 parameters...");

        if (!(sample1 && sample2 && diff)) {
          cpl_msg_error(task, "Memory allocation!");
          cpl_free(exposureTime);
          for (i = 0; i < goodFrames; i++)
            deleteImage(flatList[i]);
          cpl_free(flatList);
          deleteCcdTable(ccdTable);
          return EXIT_FAILURE;
        }

        if (readIntDescriptor(flatList[0]->descs, pilTrnGetKeyword("Quadrant"),
                              &quadrant, NULL) == VM_FALSE) {
          cpl_msg_error(task, "Integer descriptor %s not found",
                      pilTrnGetKeyword("Quadrant"));
          cpl_free(sample1);
          cpl_free(sample2);
          cpl_free(diff);
          cpl_free(exposureTime);
          for (i = 0; i < goodFrames; i++)
            deleteImage(flatList[i]);
          cpl_free(flatList);
          deleteCcdTable(ccdTable);
          return EXIT_FAILURE;
        }

        if (pilQcGroupStart() == EXIT_SUCCESS) {

          /*
           * Write ARCFILE, TPL ID, INS FILT# NAME, and OCS CON QUAD,
           * to QC1 group.
           */

          pilQcWriteString("PRO.CATG", mFlatTag, "Product category");

          qcCopyValue(flatList[0]->descs, pilTrnGetKeyword("ArchiveFile"),
                      NULL, "Archive File Name");
    
          qcCopyValue(flatList[0]->descs, pilTrnGetKeyword("TplId"),
                      NULL, "Template signature ID");

          qcCopyValue(flatList[0]->descs, pilTrnGetKeyword("Quadrant"),
                      NULL, "Quadrant");

          qcCopyValue(flatList[0]->descs,
                      pilTrnGetKeyword("FilterName", quadrant),
                      NULL, "Filter name");

          for (i = 0; i < npix; i++)
            diff[i] = sample1[i] - sample2[i];

          cpl_free(sample2);

          convFactors = cpl_malloc(ngains * sizeof(float));

          if (!convFactors) {
            cpl_msg_error(task, "Memory allocation!");
            cpl_free(sample1);
            cpl_free(diff);
            cpl_free(exposureTime);
            for (i = 0; i < goodFrames; i++)
              deleteImage(flatList[i]);
            cpl_free(flatList);
            deleteCcdTable(ccdTable);
            return EXIT_FAILURE;
          }

          for (i = 0, x = 0; x < winSizeX; x += winStepX) {
            for (y = 0; y < winSizeY; y += winStepY, i++) {
              box = extractFloatImage(sample1, winSizeX, winSizeY, x, y, 
                                      winStepX, winStepY);
              diffBox = extractFloatImage(diff, winSizeX, winSizeY, x, y, 
                                          winStepX, winStepY);

              meanValue = medianPixelvalue(box, winStepX * winStepY);
              cpl_free(box);
              value = computeVarianceFloat2D(diffBox, winStepX, winStepY);
              value /= 2;
              convFactors[i] = meanValue / value;

              cpl_free(diffBox);
            }
          }

          meanValue = computeAverageFloat(convFactors, ngains);

          value = 0.0;
          for (i = 0; i < ngains; i++)
            value += (convFactors[i] - meanValue) 
                   * (convFactors[i] - meanValue);

          value /= ngains;
          sigma = sqrt(value);

          meanValue = medianPixelvalue(convFactors, ngains);

          cpl_free(convFactors);

          qcWriteValueDouble(mFlat->descs, meanValue, "QC.CONAD", "e-/ADU", 
                             "Conversion factor from electrons to ADU");

          qcWriteValueDouble(mFlat->descs, sigma, "QC.CONAD.RMS", "e-/ADU", 
                             "Error conversion factor from electrons to ADU");

          value = computeVarianceFloat2D(diff, winSizeX, winSizeY);
          value /= 2;
          phnValue = sqrt(value);

          qcWriteValueDouble(mFlat->descs, phnValue, "QC.FLAT.PHN", "ADU",
                             "Photon noise");

          sample2 = extractFloatImage(flatList[1]->data,
                                      flatList[1]->xlen, flatList[1]->ylen,
                                      winStartX + 10, winStartY + 10,
                                      winSizeX, winSizeY);

          if (!sample2) {
            cpl_msg_error(task, "Memory allocation!");
            cpl_free(exposureTime);
            for (i = 0; i < goodFrames; i++)
              deleteImage(flatList[i]);
            cpl_free(flatList);
            deleteCcdTable(ccdTable);
            return EXIT_FAILURE;
          }

          for (i = 0; i < npix; i++)
            diff[i] = sample1[i] - sample2[i];

          cpl_free(sample2);

          value = computeVarianceFloat2D(diff, winSizeX, winSizeY);
          value /= 2;
          if (value > phnValue * phnValue)
            fpnValue = sqrt(value - phnValue * phnValue);
          else
            fpnValue = 0.0;

          qcWriteValueDouble(mFlat->descs, fpnValue, "QC.FLAT.FPN", "ADU",
                             "Fixed pattern noise");

          cpl_free(diff);

          meanValue = computeAverageFloat(sample1, npix);

          value = 0.0;
          for (i = 0; i < npix; i++)
            value += (sample1[i] - meanValue) * (sample1[i] - meanValue);

          value /= npix;
          if (value > fpnValue + phnValue)
            structValue = sqrt(value - fpnValue - phnValue);
          else
            structValue = 0.0;

          qcWriteValueDouble(mFlat->descs, structValue, "QC.FLAT.STRUCT", 
                             "ADU", "Flat field structure");

          value = medianPixelvalue(sample1, npix) / exposureTime[0];

          qcWriteValueDouble(mFlat->descs, value, "QC.FLAT.EFFICIENCY", 
                             "ADU/s", "Signal per unit exposure");

          cpl_free(sample1);

          /*
           * QC.FLAT.MASTER.MEDIAN 
           */

          winStartX = (mFlat->xlen - winSizeX) / 2;
          winStartY = (mFlat->ylen - winSizeY) / 2;

          sample1 = extractFloatImage(mFlat->data, mFlat->xlen, mFlat->ylen,
                                      winStartX, winStartY,
                                      winSizeX, winSizeY);

          if (!(sample1 && sample2 && diff)) {
            cpl_msg_error(task, "Memory allocation!");
            cpl_free(exposureTime);
            for (i = 0; i < goodFrames; i++)
              deleteImage(flatList[i]);
            cpl_free(flatList);
            deleteCcdTable(ccdTable);
            return EXIT_FAILURE;
          }

          meanValue = computeAverageFloat(sample1, npix);
          value = 0.0;
          for (i = 0; i < npix; i++)
            value += (sample1[i] - meanValue) * (sample1[i] - meanValue);

          value /= npix;

          qcWriteValueDouble(mFlat->descs, sqrt(value), "QC.FLAT.MASTER.RMS",
                             "ADU", "RMS of not normalized master screen flat");

          value = medianPixelvalue(sample1, npix);

          qcWriteValueDouble(mFlat->descs, value, "QC.FLAT.MASTER.MEDIAN",
                             "ADU", 
                             "Median of not normalized master screen flat");

          if (pilQcGroupEnd() == EXIT_FAILURE)
            cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

        }
        else
          cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

      }
      else
        cpl_msg_warning(task, 
                      "Exposure times of first two input flat fields differ: "
                      "no QC1 computation can be carried out.");
    }
    else
      cpl_msg_warning(task, 
                    "Less than two input frames survived the verification "
                    "test: no QC1 computation can be carried out.");

  }

 /*
  * Compute total exposure time of surviving frames
  */

  for (i = 0; i < goodFrames; i++)
    sumExposures += exposureTime[i];

  cpl_free(exposureTime);

  for (i = 0; i < goodFrames; i++)
    deleteImage(flatList[i]);
  cpl_free(flatList);

  cpl_msg_info(task, "Master screen flat field trend removal...");

  if (ccdTable) {
    mFlatClean = duplicateImage(mFlat);
    copyAllDescriptors(mFlat->descs, &(mFlatClean->descs));
    if (cleanBadPixels(mFlatClean, ccdTable, 0) == EXIT_FAILURE) {
      cpl_msg_error(task, "Problems modeling screen flat field continuum");
      deleteImage(mFlat);
      deleteImage(mFlatClean);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
    mFlatFilt = VmFrFilter(mFlatClean, smoothBoxSize, smoothBoxSize, 
                           filterMethod, 0);
    deleteImage(mFlatClean);
  }
  else
    mFlatFilt = VmFrFilter(mFlat, smoothBoxSize, smoothBoxSize, 
                           filterMethod, 0);

  mFlatNorm = imageArith(mFlat, mFlatFilt, VM_OPER_DIV);

  deleteImage(mFlatFilt);   /* Destroy the no longer needed smoothed flat */
  mFlatFilt = mFlatNorm;    /* Now mFlatFilt holds the flattened flat */

  cpl_msg_info(task, "Master screen flat field normalization...");

  if ((mFlatNorm = VmImNorm(mFlatFilt, MEDIAN)) == NULL) {
    cpl_msg_error(task, "Problems generating the master screen flat field");
    deleteImage(mFlat);
    deleteImage(mFlatFilt);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  deleteImage(mFlatFilt);


  copyAllDescriptors(mFlat->descs, &(mFlatNorm->descs));

  if (cleanBadPixel) {
    error = 1;
    cpl_msg_info(task, "Cleaning bad pixels from result frames...");
    if (cleanBadPixels(mFlat, ccdTable, 0) == EXIT_SUCCESS) {
      if (cleanBadPixels(mFlatNorm, ccdTable, 0) == EXIT_SUCCESS) {
        error = 0;
      }
    }

    if (error) {
      cpl_msg_error(task, "Failure cleaning bad pixels");
      deleteImage(mFlat);
      deleteImage(mFlatNorm);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  deleteCcdTable(ccdTable);

 /*
  * Update the products header
  */

  /*
   * Output the normalized master screen flat field.
   */

  vmstrlower(strcpy(mFlatName, mFlatTag));
  strcat(mFlatName, ".fits");

  /* FIXME: Is this really applicable */
  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlatNorm->descs), 
                                     pilTrnGetKeyword("ExposureTime"),
                                     1., 
                                     pilTrnGetComment("ExposureTime"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlatNorm->descs), 
                                     pilTrnGetKeyword("DataMin"), 
                                     imageMinimum(mFlatNorm), 
                                     pilTrnGetComment("DataMin"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlatNorm->descs), 
                                     pilTrnGetKeyword("DataMax"), 
                                     imageMaximum(mFlatNorm), 
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlatNorm->descs),
                                     pilTrnGetKeyword("SummedExposureTime"), 
                                     sumExposures, 
                                     pilTrnGetComment("SummedExposureTime"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertIntDescriptor(&(mFlatNorm->descs), 
                                  pilTrnGetKeyword("NFramesCombined"),
                                  goodFrames,
                                  pilTrnGetComment("NFramesCombined"),
                                  "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlatNorm->descs), 
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(mFlatNorm), 
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlatNorm->descs), 
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(mFlatNorm), 
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlatNorm->descs), 
                                     pilTrnGetKeyword("DataMean"), 
                                     imageMean(mFlatNorm), 
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    deleteImage(mFlat);
    deleteImage(mFlatNorm);
    return EXIT_FAILURE;
  }


  if (createFitsImage(mFlatName, mFlatNorm, mFlatTag) != VM_TRUE) {
    cpl_msg_error(task, "Cannot create local product file %s!", mFlatName);

    deleteImage(mFlatNorm);
    deleteImage(mFlat);

    return EXIT_FAILURE;
  }
  else {
    outputFrame = newPilFrame(mFlatName, mFlatTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }

  deleteImage(mFlatNorm);


 /*
  * Output the combination of the raw flat fields
  */

  vmstrlower(strcpy(cFlatName, cFlatTag));
  strcat(cFlatName, ".fits");

  /* FIXME: Is this really applicable */
  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("ExposureTime"),
                                     1., 
                                     pilTrnGetComment("ExposureTime"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataMin"), 
                                     imageMinimum(mFlat), 
                                     pilTrnGetComment("DataMin"), 
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataMax"), 
                                     imageMaximum(mFlat),
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs),
                                     pilTrnGetKeyword("SummedExposureTime"), 
                                     sumExposures, 
                                     pilTrnGetComment("SummedExposureTime"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertIntDescriptor(&(mFlat->descs), 
                                  pilTrnGetKeyword("NFramesCombined"),
                                  goodFrames,
                                  pilTrnGetComment("NFramesCombined"),
                                  "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(mFlat), 
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(mFlat), 
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataMean"), 
                                     imageMean(mFlat), 
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));

  if (createFitsImage(cFlatName, mFlat, cFlatTag) != VM_TRUE) {
    cpl_msg_error(task, "Cannot create local product file %s!",
                cFlatName);
    deleteImage(mFlat);

    return EXIT_FAILURE;
  }
  else {
    outputFrame = newPilFrame(cFlatName, cFlatTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }

  deleteImage(mFlat);

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
                    "vmimflatscreen",
    "Create a master screen flat field from set of raw screen flat field "
    "frames.",
    "This recipe is used to create a master screen flat field from a set\n"
    "of raw screen flat fields. The master screen flat field is not used\n"
    "directly in the flat field correction of scientific data, but it is\n"
    "optionally used just in the creation of a master sky flat field by\n"
    "the recipe vmimflatsky.\n\n"
    "Input files:\n\n"
    "  DO category:              Type:       Explanation:         Required:\n"
    "  IMG_SCREEN_FLAT           Raw         Screen flat exposure    Y\n"
    "  MASTER_BIAS               Calib       Master bias             Y\n"
    "  MASTER_DARK               Calib       Master dark             .\n"
    "  CCD_TABLE                 Calib       Bad pixel table         .\n\n"
    "Output files:\n\n"
    "  DO category:              Data type:  Explanation:\n"
    "  IMG_MASTER_SCREEN_FLAT    FITS image  Master screen flat field\n"
    "  IMG_COMBINED_SCREEN_FLAT  FITS image  Combined screen flat field\n\n"
    "The primary product is the normalised master screen flat field.\n"
    "A secondary product is the combined screen flat field, that is the\n"
    "result of the combination of all inputs but without any normalisation\n"
    "applied, and is just used for data quality control. A CCD table must\n"
    "be specified in input only if a bad pixel cleaning is requested.\n\n"
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

                    vmimflatscreen_create,
                    vmimflatscreen_exec,
                    vmimflatscreen_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
