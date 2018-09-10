/* $Id: vmimflatsky.c,v 1.5 2013-08-07 16:45:33 cgarcia Exp $
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
#define SATURATION_LEVEL (65535.)


static cxint vmimflatsky(PilSetOfFrames *);


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
vmimflatsky_create(cpl_plugin *plugin)
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
                               "Median", 5, "Average", "Median", "MinMax",
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
                                50000.);
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
                               "Bias removal method.",
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

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.quality.apply",
                                CPL_TYPE_BOOL,
                                "Quality control of the master flat.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyQC");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

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
vmimflatsky_exec(cpl_plugin *plugin)
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

    if (vmimflatsky(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmimflatsky");
        
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
vmimflatsky_destroy(cpl_plugin *plugin)
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
 *   Create master flat field from a set of raw sky flat fields.
 * 
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param sof  Set of frames containing the references to raw sky
 *             flat field frames, an optional master screen flat
 *             field, a master bias, an optional master dark, and 
 *             an optional CCD table.
 *
 * @doc 
 *   The recipe function creates a master flat field frame from 
 *   a set of raw sky flat field frames, after bias removal and 
 *   optional dark subtraction. Optionally, an input master screen 
 *   flat field may be used to better model the small scale variations
 *   of the instrument response. A CCD table may also be specified, 
 *   to be used in bad pixel correction and/or quality check of the 
 *   output master screen flat field. All the input frames are passed 
 *   via the set of frames \textbf{sof}. If requested, input sky flat 
 *   fields consistency is checked, and cosmic ray events are removed 
 *   before combination. If an input master screen flat field was
 *   given, the combined sky flat fields would be smoothed, and then
 *   multiplied by the master screen flat field; otherwise, it would
 *   be just normalized to its mean value. On successful termination 
 *   the created master flat field is added to the set of frames.
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
 *
 *   \begin{itemize}
 *
 *     \item SmoothMethod        Method for smoothing the combined input 
 *                               flat fields, for determining the large 
 *                               scale trends related to the CCD intrinsic 
 *                               response variation. This is used only if a
 *                               master screen flat field is also supplied
 *                               in input. Legal settings are:
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
 *     \item StackMethod:        Method for raw sky flat field frames 
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
 *       \item Ksigma:             Kappa-Sigma clipping. If this option is
 *                                 chosen, the upper and lower rejection
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
 *     \item BiasMethod:         Method for bias removal from raw sky flat
 *                               field frames. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Master:             Prescan and overscan regions are trimmed
 *                                 away from raw sky flat fields after 
 *                                 master bias subtraction.
 *
 *       \item Zmaster:            After master bias subtraction the residual
 *                                 signal found in each raw sky flat field
 *                                 overscan regions is modelled and subtracted 
 *                                 from the image. Next, prescan and overscan 
 *                                 regions are trimmed away.
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
 *                               based on pattern differences, in units
 *                               of sigmas.
 *
 *     \item CleanBadPixel:      Bad pixel correction on output master sky
 *                               flat field. If this option is turned on, 
 *                               the recipe expects to find a CCD\_TABLE 
 *                               in the input set of frames.
 *
 *     \item CleanCosmic:        Cosmic ray events removal from input raw 
 *                               sky flat field frames. If a CCD\_TABLE is
 *                               found among the input set of frames, bad
 *                               pixels will not be used in computing the
 *                               values to replace the cosmic rays events.
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
 *     \item ApplyQC:            Apply a simple quality control check to
 *                               the created master flat field frame
 *                               (CURRENTLY NOT AVAILABLE).
 *
 *   \end{itemize}
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 *
 * @author C. Izzo, P. Sartoretti, R. Palsa
 */   

static cxint 
vmimflatsky(PilSetOfFrames *sof)
{

  const char     task[] = "vmimflatsky";

  const char     parameter[] = "Parameters";

  const char    *flatTag  = pilTrnGetCategory("ImgSkyFlat");
  const char    *mflatTag = pilTrnGetCategory("ImgMasterSkyFlat");

  char          *combMethodTag = NULL;
  char          *biasMethodTag = NULL;
  char          *filterMethodTag = NULL;
  char           masterFlatName[PATHNAME_MAX + 1];
  char           comment[MAX_COMMENT_LENGTH];

  VimosBool      updateOK = VM_TRUE;

  size_t         flatCount, goodFrames, minFrames;

  unsigned int   i, j;
  unsigned int   cleanBadPixel, cleanCosmic, validateFrames, applyQc;
  unsigned int   computeQC, singleFrames;
  unsigned int   error;

  int            biasMethodEntry;
  int            combMethodEntry;
  int            filterMethodEntry;

  float         *flatLevel = NULL;
  float         *flatNoise = NULL;
  float         *flatRon = NULL;
  float         *flatGain = NULL;
  float         *exposureTime = NULL;

  double         time, gain;
  double         sumExposures = 0.0;
  double         rescaleFactor;
  double         lowerThreshold, upperThreshold;
  double         patternTolerance;
  double         cosmicThreshold = 0.0;
  double         cosmicRatio = 0.0;

  PilFrame      *currFrame, *ccdFrame, *biasFrame, *darkFrame, *screenFrame;
  PilFrame      *outputFrame;

  BiasMethod     biasMethod = BIAS_UNDEF;

  CombMethod     combMethod = COMB_UNDEF;
  CombParameters combParameter;

  CombMethod     filterMethod = FILTER_UNDEF;
  int            smoothBoxSize;
  int            quadrant;

  VimosImage    *mFlat = NULL;
  VimosImage   **flatList, **original;
  VimosImage    *biasImage = NULL;
  VimosImage    *darkImage = NULL;
  VimosImage    *screenImage = NULL;
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
   * Check if the flat field quality control should be applied to the
   * final product (TO BE IMPLEMENTED).
   */

  applyQc = pilDfsDbGetBool(parameter, "ApplyQC", 0);


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
      cpl_msg_error(task, "Not enough raw sky flat field frames in "
                  "input for stacking method '%s'!", combMethodTag);
      return EXIT_FAILURE;
    }
  }

  if (validateFrames && flatCount < 2) {
    validateFrames = 0;
    cpl_msg_warning(task, "Too few sky flat field frames (%zd) in input. "
                  "Skipping frame validation task!", flatCount);
  }

  cpl_msg_info(task, "Loading input frames...");


 /*
  * Get the (optional) master screen flat field frame
  */

  if ((screenFrame = pilSofLookup(sof,
				  pilTrnGetCategory("ImgMasterScreenFlat")))) {
    pilFrmSetType(screenFrame, PIL_FRAME_TYPE_CALIB);
    if ((screenImage = openOldFitsFile(pilFrmGetName(screenFrame), 1, 0))) {
      closeFitsImage(screenImage, 0);
    }
    else {
      cpl_msg_error(task, "Failure opening master screen flat field frame");
      return EXIT_FAILURE;
    }
  }
  else 
    cpl_msg_warning(task, "No master screen flat field in input, input sky "
                  "flat fields will be just stacked");


  /*
   * If bad pixel correction is enabled, a bad pixel table is required 
   * in the input set of frames. If no bad pixel table is present this 
   * is an error. If cleaning of cosmic rays is enabled, the bad pixel 
   * table will be used, if present, to avoid bad pixels while smoothing 
   * away cosmic ray events. In this case, a missing bad pixel table is 
   * not an error.
   */

  ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable"));
  if (ccdFrame)
    pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

  if (cleanBadPixel || cleanCosmic || screenImage) {
    if (ccdFrame) {
      cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));
      if ((ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
        closeFitsTable(ccdTable, 0);
      }
      else {
        cpl_msg_error(task, "Failure in opening CCD table");
        deleteImage(screenImage);
        return EXIT_FAILURE;
      }
    }
    else {
      if (cleanBadPixel) {
        cpl_msg_error(task, "No CCD table in input: cannot clean bad pixels.");
        return EXIT_FAILURE;
      }
      cpl_msg_warning(task, "No CCD table in input: flat field continuum "
                    "may include ghosts from bad pixels.");
    }
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
    deleteImage(screenImage);
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
      deleteImage(screenImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }
  else 
    cpl_msg_warning(task, "No master dark in input, dark subtraction "
                  "will not be performed");


 /*
  * Load the raw sky flat field frames.
  */

  error = 1;

  if ((flatList = (VimosImage **)cpl_calloc(flatCount, sizeof(VimosImage *)))) {

    error = 0;
    currFrame = pilSofLookupNext(sof, flatTag);

    for (i = 0; i < flatCount; i++) {
      if ((flatList[i] = openOldFitsFile(pilFrmGetName(currFrame), 1, 0))) {
        pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
        closeFitsImage(flatList[i], 0);
      }
      else {
        error = 1;
        cpl_msg_error(task, "Failure opening sky flat field frame %d", i + 1);
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
                "Failure in creating list of input raw sky flat fields");

  if (error) {
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteImage(screenImage);
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
    deleteImage(screenImage);
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
    deleteImage(screenImage);
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
    deleteImage(screenImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Noise contributions for each input sky flat field frame must be 
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
    deleteImage(screenImage);
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
    deleteImage(screenImage);
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
      deleteImage(screenImage);
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
        deleteImage(screenImage);
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
    deleteImage(screenImage);
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
        cpl_msg_error(task, "Cannot clean cosmic rays from raw sky "
                    "flat field frame %d", i + 1);
        for (i = 0; i < flatCount; i++)
          deleteImage(flatList[i]);
        cpl_free(flatList);
        cpl_free(flatNoise);
        cpl_free(flatGain);
        cpl_free(flatRon);
        cpl_free(flatLevel);
        cpl_free(exposureTime);
        deleteImage(screenImage);
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
      cpl_msg_info(task, "Valid range for mean sky flat field flux: "
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
        deleteImage(screenImage);
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
        deleteImage(screenImage);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }
    }

    if (goodFrames) {
      cpl_msg_info(task, "%zd out of %zd consistent sky flat field frames "
                 "selected from input set", goodFrames, flatCount);
    }

    cpl_free(original);

  }
  else {
    goodFrames = flatCount;

    /*
     * All flat fields values are normalized to a common level
     * (arbitrarily, the one of the first flat field image).
     */

    for (i = 1; i < goodFrames; i++) {
      rescaleFactor = flatLevel[0] / flatLevel[i];
      constArithLocal(flatList[i], rescaleFactor, VM_OPER_MUL);
    }
  }

 /*
  * Destroy inconsistent frames from the rearranged list of flat field
  * frames, free arrays.
  */

  for (i = goodFrames; i < flatCount; i++)
    deleteImage(flatList[i]);
  cpl_free(flatNoise);

  if (!goodFrames) {
    cpl_msg_error(task, "No sky flat field frames with consistent signal "
                "distribution found in input. All frames excluded!");
    cpl_free(flatList);
    cpl_free(exposureTime);
    deleteImage(screenImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Stack the surviving flat field frames, if requested, and if possible.
  */

  if (combMethod == COMB_UNDEF) {
    mFlat = duplicateImage(flatList[0]);
  }

  if (goodFrames >= minFrames) {
    cpl_msg_info(task, "Combining %zd frames with method '%s'", goodFrames,
               combMethodNames[combMethodEntry]);

    if (combMethod == COMB_AVERAGE) {

     /* 
      * Average is weighted on the original levels - i.e., eliminate
      * the rescaling that was done in case validateFrames was "true".
      */ 

/*
      if (validateFrames) {
*/
        for (i = 1; i < flatCount; i++) {
          rescaleFactor = flatLevel[i] / flatLevel[0];
          constArithLocal(flatList[i], rescaleFactor, VM_OPER_MUL);
        }
/*
      }
*/
      mFlat = frCombAverage(flatList, goodFrames);

    }
    else {
      mFlat = frComb(flatList, goodFrames, combMethod, &combParameter, 0);
    }
  }
  else {
    cpl_msg_error(task, "Not enough good frames for combine method "
                    "'%s'!", combMethodNames[combMethodEntry]);
    cpl_free(exposureTime);
    for (i = 0; i < goodFrames; i++)
      deleteImage(flatList[i]);
    cpl_free(flatList);
    deleteImage(screenImage);
    deleteCcdTable(ccdTable);
    cpl_free(flatLevel);
    return EXIT_FAILURE;
  }

  cpl_free(flatLevel);

  if (!mFlat) {
    cpl_msg_error(task, "Stacking of raw sky flat field frames failure!");
    deleteImage(screenImage);
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

/*  copyAllDescriptors(flatList[0]->descs, &(mFlat->descs));   */


  if (computeQC) {

    /*
     *  Currently the central 1600x1800 pixels are used
     */

    int          winSizeX = 1600;
    int          winSizeY = 1800;
    int          winStartX = (flatList[0]->xlen - winSizeX) / 2;
    int          winStartY = (flatList[0]->ylen - winSizeY) / 2;

    double       meanValue;

    float       *sample1;


    cpl_msg_info(task, "Computing QC1 parameters...");

    if (readIntDescriptor(flatList[0]->descs, pilTrnGetKeyword("Quadrant"),
                          &quadrant, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Integer descriptor %s not found",
                  pilTrnGetKeyword("Quadrant"));
      cpl_free(exposureTime);

      for (i = 0; i < goodFrames; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);

      deleteImage(screenImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }

    for (i = 0; i < goodFrames; i++) {

      if (pilQcGroupStart() == EXIT_SUCCESS) {

        sample1 = extractFloatImage(flatList[i]->data,
                                    flatList[i]->xlen, flatList[i]->ylen,
                                    winStartX, winStartY, winSizeX, winSizeY);

        if (!sample1) {
          cpl_msg_error(task, "Memory allocation!");

          cpl_free(exposureTime);

          for (i = 0; i < goodFrames; i++)
            deleteImage(flatList[i]);
          cpl_free(flatList);

          deleteImage(screenImage);
          deleteCcdTable(ccdTable);
          return EXIT_FAILURE;
        }

        /*
         * Write ARCFILE, TPL ID, INS FILT# NAME, and OCS CON QUAD, 
         * to QC1 group.
         */

        pilQcWriteString("PRO.CATG", mflatTag, "Product category");
        pilQcWriteString("DPR.CATG", flatTag, "DO category");

        qcCopyValue(flatList[i]->descs, pilTrnGetKeyword("ArchiveFile"),
                    NULL, "Archive File Name");

        qcCopyValue(flatList[i]->descs, pilTrnGetKeyword("TplId"),
                    NULL, "Template signature ID");

        qcCopyValue(flatList[i]->descs, pilTrnGetKeyword("Quadrant"),
                    NULL, "Quadrant");

        qcCopyValue(flatList[i]->descs, 
                    pilTrnGetKeyword("FilterName", quadrant), 
                    NULL, "Filter name");


        /*
         * Compute QC parameters.
         */

        /* QC.SKY.FLAT.FLUX */

        meanValue = computeAverageFloat(sample1, winSizeX * winSizeY);

        meanValue /= exposureTime[i];

        if (pilQcWriteDouble("QC.SKY.FLAT.FLUX", meanValue, 
                             "ADU/s", "Mean sky flux") == EXIT_FAILURE)
          cpl_msg_error(task, "Could not copy value to QC1 PAF! (Ignored)");

        cpl_free(sample1);

        if (pilQcGroupEnd() == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

      }
      else
        cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

    }

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

  if (screenImage) {
    cpl_msg_info(task, "Smooth the combined flat field...");

   /*
    * Flat field correction:
    */

    imageArithLocal(mFlat, screenImage, VM_OPER_DIV);

   /*
    * Smoothing:
    */

    if (ccdTable) {
      mFlatClean = duplicateImage(mFlat);
      copyAllDescriptors(mFlat->descs, &(mFlatClean->descs));
      if (cleanBadPixels(mFlatClean, ccdTable, 0) == EXIT_FAILURE) {
        cpl_msg_error(task, "Problems modeling flat field continuum");
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

    if (mFlatFilt)
      copyAllDescriptors(mFlat->descs, &(mFlatFilt->descs));

    deleteImage(mFlat);
    mFlat = mFlatFilt;

    cpl_msg_info(task, "Apply to the smoothed flat field the master "
               "screen flat field pattern...");

    imageArithLocal(mFlat, screenImage, VM_OPER_MUL);

    deleteImage(screenImage);
  }

  cpl_msg_info(task, "Master flat field normalization...");

  if ((mFlatNorm = VmImNorm(mFlat, MEDIAN)) == NULL) {
    cpl_msg_error(task, "Problems generating the master flat field");
    deleteImage(mFlat);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  copyAllDescriptors(mFlat->descs, &(mFlatNorm->descs));

  deleteImage(mFlat);

  if (cleanBadPixel) {
    cpl_msg_info(task, "Cleaning bad pixels from result frame...");
    if (cleanBadPixels(mFlatNorm, ccdTable, 0) == EXIT_FAILURE) {
      cpl_msg_error(task, "Failure cleaning bad pixels");
      deleteImage(mFlatNorm);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  if (computeQC) {

    /*
     *  Currently the central 1600x1800 pixels are used
     */

    int          winSizeX = 1600;
    int          winSizeY = 1800;
    int          winStartX = (mFlatNorm->xlen - winSizeX) / 2;
    int          winStartY = (mFlatNorm->ylen - winSizeY) / 2;
    unsigned int npix = winSizeX * winSizeY;

    double       value, meanValue, structValue;

    float       *sample1 = extractFloatImage(mFlatNorm->data,
                           mFlatNorm->xlen, mFlatNorm->ylen,
                           winStartX, winStartY, winSizeX, winSizeY);
    float       *sample2 = extractFloatImage(mFlatNorm->data,
                           mFlatNorm->xlen, mFlatNorm->ylen,
                           winStartX + 10, winStartY + 10, winSizeX, winSizeY);
    float       *diff    = cpl_malloc(winSizeX * winSizeY * sizeof(float));


    if (!(sample1 && sample2 && diff)) {
      cpl_msg_error(task, "Memory allocation!");
      deleteImage(mFlatNorm);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      /*
       * Write ARCFILE, TPL ID, INS FILT# NAME, and OCS CON QUAD, 
       * to QC1 group.
       */

      pilQcWriteString("PRO.CATG", mflatTag, "Product category");

      qcCopyValue(mFlatNorm->descs, pilTrnGetKeyword("ArchiveFile"),
                  NULL, "Archive File Name");

      qcCopyValue(mFlatNorm->descs, pilTrnGetKeyword("TplId"),
                  NULL, "Template signature ID");

      qcCopyValue(mFlatNorm->descs, pilTrnGetKeyword("Quadrant"),
                  NULL, "Quadrant");

      qcCopyValue(mFlatNorm->descs,
                  pilTrnGetKeyword("FilterName", quadrant),
                  NULL, "Filter name");

      meanValue = computeAverageFloat(sample1, npix);

      structValue = 0.0;
      for (i = 0; i < npix; i++)
        structValue += (sample1[i] - meanValue) * (sample1[i] - meanValue);

      structValue /= npix;

      qcWriteValueDouble(mFlatNorm->descs, sqrt(structValue), 
                         "QC.SKY.FLAT.RMS", NULL, 
                         "Structure of master sky flat field");

      for (i = 0; i < npix; i++)
        diff[i] = sample1[i] - sample2[i];

      cpl_free(sample1);
      cpl_free(sample2);

      value = computeVarianceFloat2D(diff, winSizeX, winSizeY);
      value /= 2.;
 
      structValue -= value;

      if (structValue > 0.0)
        structValue = sqrt(structValue);
      else
        structValue = 0.0;

      qcWriteValueDouble(mFlatNorm->descs, structValue, "QC.SKY.FLAT.STRUCT", 
                         NULL, "Structure of master sky flat field");

      cpl_free(diff);

      if (pilQcGroupEnd() == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");
    }
    else
      cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

  }


 /*
  * Apply the quality control check to the created master flat field.
  */

  if (applyQc) {
    cpl_msg_info(task, "Quality check on created master flat field - "
               "NOT YET IMPLEMENTED");

/*** TEMPORARILY COMMENTED OUT ***
    if (qcCheckFlatLevel(mFlat, ccdTable, ...??) == EXIT_FAILURE) {
      cpl_msg_error(task, "Quality check on master sky flat field failed!");

      deleteCcdTable(ccdTable);
      deleteImage(mFlat);

      return EXIT_FAILURE;
    }
***/

  }

  deleteCcdTable(ccdTable);

 /*
  * Update the products header
  */

 /*
  * Output the normalized master flat field:
  */

  /* FIXME: Is this really applicable */
  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlatNorm->descs), 
                                     pilTrnGetKeyword("ExposureTime"), 1., 
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
    deleteImage(mFlatNorm);
    return EXIT_FAILURE;
  }


  /*
   * Create the product file on disk, set the product attributes and
   * update the set of frames.
   */

  vmstrlower(strcpy(masterFlatName, mflatTag));
  strcat(masterFlatName, ".fits");

  if (createFitsImage(masterFlatName, mFlatNorm, mflatTag) == VM_FALSE) {
      cpl_msg_error(task, 
                  "Cannot create local product file %s!", masterFlatName);
      deleteImage(mFlatNorm);
      return EXIT_FAILURE;
  }
  else {
    outputFrame = newPilFrame(masterFlatName, mflatTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }


  /*
   * Cleanup.
   */

  deleteImage(mFlatNorm);

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
                    "vmimflatsky",
    "Create a master sky flat field from set of raw sky flat field frames.",
    "This recipe is used to create a master sky flat field from a set of\n"
    "raw sky flat fields. The master sky flat field is the dataset used\n"
    "for the flat field correction of scientific data.\n\n"
    "Input files:\n\n"
    "  DO category:              Type:       Explanation:         Required:\n"
    "  IMG_SKY_FLAT              Raw         Sky flat exposure       Y\n"
    "  MASTER_BIAS               Calib       Master bias             Y\n"
    "  MASTER_DARK               Calib       Master dark             .\n"
    "  IMG_MASTER_SCREEN_FLAT    Calib       Master screen flat      .\n"
    "  CCD_TABLE                 Calib       Bad pixel table         .\n\n"
    "Output files:\n\n"
    "  DO category:              Data type:  Explanation:\n"
    "  IMG_MASTER_SKY_FLAT       FITS image  Master sky flat field\n\n"
    "If a master screen flat is specified, the input sky flat exposures\n"
    "are flat-fielded with it, and then just used to determine the large\n"
    "scale trends of the instrument response; the result is then multiplied\n"
    "by the master screen flat. A CCD table must be specified in input only\n"
    "if a bad pixel cleaning is requested.\n\n"
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

                    vmimflatsky_create,
                    vmimflatsky_exec,
                    vmimflatsky_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
