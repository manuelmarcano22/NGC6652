/* $Id: vmdark.c,v 1.6 2013-08-07 16:41:46 cgarcia Exp $
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
 * $Date: 2013-08-07 16:41:46 $
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <string.h>
#include <math.h>

#include <cxmemory.h>

#include <cpl.h>

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


#define SOME_PERCENT (0.01)
#define MAX_COMMENT_LENGTH (80)


static cxint vmdark(PilSetOfFrames *);


/*
 * Definition of the label strings for all methods the recipe function
 * supports for combining frames and their associated method code.
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

static const char *biasMethodNames[] = {
  "Master",
  "Zmaster"
};

static const BiasMethod biasMethods[] = {
  BIAS_MASTER,
  BIAS_ZMASTER
};


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmdark_create(cpl_plugin *plugin)
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
                               "Frame combination method.",
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

    p = cpl_parameter_new_value("vimos.Parameters.dark.frames.validate",
                                CPL_TYPE_BOOL,
                                "Consistency check on raw dark frames.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ValidateFrames");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ValidateFrames");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.dark.tolerance.level",
                                CPL_TYPE_DOUBLE,
                                "Threshold for level consistency.",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LevelTolerance");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LevelTolerance");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.dark.tolerance.pattern",
                                CPL_TYPE_DOUBLE,
                                "Threshold for flux pattern consistency.",
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
                                "Bad pixel correction on master dark.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.clean",
                                CPL_TYPE_BOOL,
                                "Cosmic ray events cleaning from raw darks.",
                                "vimos.Parameters",
                                TRUE);
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
                                "Compute QC1 parameters",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ComputeQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ComputeQC");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.quality.apply",
                                CPL_TYPE_BOOL,
                                "Quality control of the master dark.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyQC");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.quality.maxoffset",
                                CPL_TYPE_DOUBLE,
                                "Maximum allowed deviation from nominal "
                                "dark level",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MaxDeviation");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MaxDeviation");
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
vmdark_exec(cpl_plugin *plugin)
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

    if (vmdark(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmdark");
        
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
vmdark_destroy(cpl_plugin *plugin)
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
 *   Create master dark from a set of raw dark frames.
 * 
 * @return EXIT_SUCCESS or EXIT_FAILURE.   
 *
 * @param sof  Set of frames containing the references to raw dark
 *             frames, a master bias, and (optionally) a CCD Table. 
 *  
 * @doc
 *   The recipe function creates a master dark frame by stacking a set
 *   of raw dark frames, after bias removal. The raw frames and the
 *   master bias frame are passed via the set of frames \textbf{sof}.
 *   A CCD table may also be specified, to be used in bad pixel 
 *   correction and/or quality check of the output master dark. 
 *   If requested, input frames consistency is checked, and cosmic 
 *   ray events are removed from each input frame before combination.
 *   On successful termination the created master dark is added to 
 *   the set of frames.
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
 *
 *   \begin{itemize}
 *
 *     \item StackMethod:      Method for raw dark frames combination.
 *                             Legal settings are:
 *
 *       \begin{itemize}
 *
 *         \item Average:        Average of frames
 *
 *         \item Median:         Median stacking
 *
 *         \item MinMax:         Minimum-Maximum rejection. If this 
 *                               option is chosen, the upper and lower 
 *                               percentiles of pixel values to be 
 *                               rejected before averaging may be 
 *                               specified using the parameters 
 *                               MinRejection and MaxRejection. 
 *                               The default is to reject just the 
 *                               highest and the lowest pixel value.
 *
 *         \item Ksigma:         Kappa-Sigma clipping. If this option 
 *                               is chosen, the upper and lower rejection 
 *                               thresholds, given in units of sigma, may 
 *                               be specified using the parameters KSigmaLow 
 *                               and KSigmaHigh. Pixel values deviating
 *                               from the median more than the specified
 *                               thresholds are rejected before final
 *                               averaging.
 *
 *         \item Auto:           Automatic method (not yet available, it 
 *                               falls back to Average).
 *
 *     \end{itemize}
 *
 *     \item BiasMethod:       Method for bias removal from raw dark frames.
 *                             Legal settings are:
 *
 *       \begin{itemize}
 *
 *         \item Master:         Prescan and overscan regions are trimmed
 *                               away from raw darks after master bias 
 *                               subtraction.
 *                             
 *         \item Zmaster:        After master bias subtraction the residual
 *                               signal found in each raw dark overscan 
 *                               regions is modelled and subtracted from 
 *                               the image. Next, prescan and overscan 
 *                               regions are trimmed away.
 *
 *     \end{itemize}
 *                             
 *     \item KSigmaLow:        Lower rejection threshold, given in units 
 *                             of sigma, used when parameter StackMethod 
 *                             is set to Ksigma.
 *
 *     \item KSigmaHigh:       Upper rejection threshold, given in units 
 *                             of sigma, used when parameter StackMethod 
 *                             is set to Ksigma.
 *
 *     \item MinRejection:     Lower percent of sorted pixel values to 
 *                             be rejected, when parameter StackMethod 
 *                             is set to MinMax. If set to 0.0, just
 *                             one pixel value, the lowest, is rejected.
 *
 *     \item MaxRejection:     Upper percent of sorted pixel values to
 *                             be rejected, when parameter StackMethod 
 *                             is set to MinMax. If set to 0.0, just
 *                             one pixel value, the highest, is rejected.
 *
 *     \item ValidateFrames:   Rejection of raw dark frames departing
 *                             from the mean dark level, measured across 
 *                             all the frames, more than a given threshold.
 *
 *     \item LevelTolerance:   Threshold controlling frames rejection
 *                             based on their mean flux, in units of sigmas.
 *
 *     \item PatternTolerance: Threshold controlling frames rejection
 *                             based on pattern differences, in units of 
 *                             sigmas.
 *
 *     \item CleanBadPixel:    Bad pixel correction on output master dark.
 *                             If this option is turned on, the recipe 
 *                             expects to find a CCD\_TABLE in the input 
 *                             set of frames.
 *
 *     \item CleanCosmic:      Cosmic ray events removal from input raw
 *                             dark frames. If a CCD\_TABLE is found among
 *                             the input set of frames, bad pixels will not
 *                             be used in computing the values to replace
 *                             the cosmic rays events.
 *
 *     \item CosmicThreshold:  Threshold for the detection of cosmic rays,
 *                             having the same meaning of the parameter ns
 *                             of the MIDAS command FILTER/COSMIC: it is
 *                             the number of theoretical noise sigmas above 
 *                             smoothed level that must be reached to make 
 *                             a pixel a cosmic ray event location candidate.
 *                             This parameter is effective when CleanCosmic
 *                             is set to "true".
 *
 *     \item CosmicRatio:      Critical ratio for discrimination of objects 
 *                             and cosmic rays, having the same meaning 
 *                             of the parameter rc of the MIDAS command 
 *                             FILTER/COSMIC: it is the ratio between the
 *                             peak of a cosmic ray event candidate, and the
 *                             average of the 8 nearby pixels that must be
 *                             reached to identify the candidate as a real
 *                             cosmic ray event.
 *
 *     \item ApplyQC:          Apply a simple quality control check to the 
 *                             created master dark frame. The test compares 
 *                             the dark level of the master dark with the 
 *                             nominal dark constant. If the offset from 
 *                             the nominal level is larger than MaxDeviation 
 *                             sigma a warning is issued, but the function 
 *                             does not fail. The nominal dark constant is 
 *                             taken from a CCD\_TABLE. If this option is 
 *                             turned on the recipe expects to find a 
 *                             CCD\_TABLE in input. 
 *
 *     \item MaxDeviation:     Maximum allowed deviation of result master 
 *                             dark from nominal dark level, given in 
 *                             terms of sigmas.
 *
 *   \end{itemize}
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 *
 * @author C. Izzo, R. Palsa, P. Sartoretti
 */

static cxint 
vmdark(PilSetOfFrames *sof)
{

  const char     task[] = "vmdark";

  const char     parameter[] = "Parameters";

  const char    *darkTag  = pilTrnGetCategory("Dark");
  const char    *mdarkTag = pilTrnGetCategory("MasterDark");

  char          *combMethodTag = NULL;
  char          *biasMethodTag = NULL;
  char           masterDarkName[PATHNAME_MAX + 1];
  char           comment[MAX_COMMENT_LENGTH];

  VimosBool      updateOK = VM_TRUE;

  size_t         darkCount, goodFrames, minFrames;

  unsigned int   i, j;
  unsigned int   cleanBadPixel, cleanCosmic, validateFrames;
  unsigned int   computeQC, applyQc;
  unsigned int   singleFrames;
  unsigned int   error;

  int            biasMethodEntry;
  int            combMethodEntry;

  float         *darkLevel = NULL;
  float         *darkNoise = NULL;
  float         *darkRon = NULL;
  float         *darkGain = NULL;
  float         *exposureTime = NULL;
  float          maxExposure, minExposure;

  double         time, gain;
  double         sumExposures = 0.0;
  double         medianLevel, maxNoise;
  double         levelTolerance = 0.0;
  double         levelThreshold = 0.0;
  double         lowerThreshold = 0.0;
  double         upperThreshold = 0.0;
  double         patternTolerance = 0.0;
  double         cosmicThreshold = 0.0;
  double         cosmicRatio = 0.0;
  double         maxDeviation;

  PilFrame      *currFrame, *ccdFrame, *mBiasFrame, *outputFrame;

  BiasMethod     biasMethod = BIAS_UNDEF;
  CombMethod     combMethod = COMB_UNDEF;
  CombParameters combParameter;

  VimosImage    *mDark = NULL;
  VimosImage   **darkList, **original, *mBiasImage, *mBiasOver;

  VimosTable    *ccdTable = NULL;


  /*
   * Get task parameters from the recipe database
   */

  /*
   * Determine the frame bias removal method first.
   */

  biasMethodTag = (char *)pilDfsDbGetString(parameter, "BiasMethod");

  if ((biasMethodEntry = strselect(biasMethodTag, biasMethodNames, 
                                   PIL_N_ELEMENTS(biasMethods))) < 0) {
    cpl_msg_error(task, "%s: Invalid bias removal method.", biasMethodTag);
    return EXIT_FAILURE;
  }

  biasMethod = biasMethods[biasMethodEntry];


  /*
   * Determine the frame stacking method and all
   * method dependent parameters.
   */

  combMethodTag = (char *)pilDfsDbGetString(parameter, "StackMethod");

  if ((combMethodEntry = strselect(combMethodTag, combMethodNames, 
                                   PIL_N_ELEMENTS(combMethods))) < 0) {
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
   * Check if raw frames with non-nominal exposure level should be ignored
   * by the image stacking task. Default is to remove non-nominal frames.
   */

  validateFrames = pilDfsDbGetBool(parameter, "ValidateFrames", 0);

  if (validateFrames) {
    levelTolerance = pilDfsDbGetDouble(parameter, "LevelTolerance", 3.0);
    patternTolerance = pilDfsDbGetDouble(parameter, "PatternTolerance", 3.0);

    if (levelTolerance < 1.0 || patternTolerance < 1.0) {
      cpl_msg_error(task, "Invalid tolerance, it should be at least one sigma");
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
   * Check if the dark quality control should be applied to the
   * final product.
   */

  applyQc = pilDfsDbGetBool(parameter, "ApplyQC", 0);
  maxDeviation = pilDfsDbGetDouble(parameter, "MaxDeviation", 3.0);


  /*
   * Make sure that there are enough raw data frames in the
   * input set.
   */

  darkCount = (int)pilSofFrameCount(sof, darkTag);

  if (darkCount < minFrames) {
    if (darkCount == 1 && singleFrames) {
      combMethod = COMB_UNDEF;
    }
    else {
      cpl_msg_error(task, "Not enough raw dark frames in input for stacking "
                  "method '%s'!", combMethodTag);
      return EXIT_FAILURE;
    }
  }

  if (validateFrames && darkCount < 2) {
    validateFrames = 0;
    cpl_msg_warning(task, "Too few dark frames (%zd) in input. Skipping "
                  "frame validation task!", darkCount);
  }

  cpl_msg_info(task, "Loading input frames...");

  /*
   * If bad pixel correction or the quality control check are enabled,
   * a bad pixel table is required in the input set of frames. If no 
   * bad pixel table is present this is an error. If cleaning of cosmic 
   * rays is enabled, the bad pixel table will be used, if present, to
   * avoid bad pixels while smoothing away cosmic ray events. In this
   * case, a missing bad pixel table is not an error.
   */

  ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable"));
  if (ccdFrame)
      pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

  if (cleanBadPixel || cleanCosmic || applyQc) {
    if (ccdFrame) {
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
      if (cleanBadPixel || applyQc) {
        cpl_msg_error(task, "No CCD table in input:");
        if (cleanBadPixel) cpl_msg_error(task, "Cannot clean bad pixels.");
        if (applyQc) cpl_msg_error(task, "Cannot quality check product.");
        return EXIT_FAILURE;
      }
    }
  }


 /*
  * Get the master bias frame
  */

  error = 1;

  if ((mBiasFrame = pilSofLookup(sof, pilTrnGetCategory("MasterBias")))) {
    pilFrmSetType(mBiasFrame, PIL_FRAME_TYPE_CALIB);
    if ((mBiasImage = openOldFitsFile(pilFrmGetName(mBiasFrame), 1, 0))) {
      closeFitsImage(mBiasImage, 0);
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
  * Load the raw dark frames.
  */

  error = 1;

  if ((darkList = (VimosImage **)cpl_calloc(darkCount, sizeof(VimosImage *)))) {

    error = 0;
    currFrame = pilSofLookupNext(sof, darkTag);

    for (i = 0; i < darkCount; i++) {
      if ((darkList[i] = openOldFitsFile(pilFrmGetName(currFrame), 1, 0))) {
        pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
        closeFitsImage(darkList[i], 0);
      }
      else {
        error = 1;
        cpl_msg_error(task, "Failure opening dark frame %d", i + 1);
        for (j = 0; j < i; j++) 
          deleteImage(darkList[j]);
        cpl_free(darkList);
        break;
      }
      currFrame = pilSofLookupNext(sof, NULL);
    }
  }
  else
    cpl_msg_error(task, "Failure in creating list of input raw darks");

  if (error) {
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable); 
    return EXIT_FAILURE;
  }

 /*
  * Recreate bias overscans using as reference the first input dark
  */

  if ((mBiasOver = growOverscans(mBiasImage, darkList[0]))) {
    if (mBiasImage != mBiasOver) {
      deleteImage(mBiasImage);
      mBiasImage = mBiasOver;
    }
  }
  else {
    cpl_msg_error(task, "Failure in growing overscans in master bias");
    for (i = 0; i < darkCount; i++) 
      deleteImage(darkList[i]);
    cpl_free(darkList);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Get exposure times
  */

  error = 1;

  if ((exposureTime = (float *)cpl_calloc(darkCount, sizeof(float)))) {
    error = 0;
    for (i = 0; i < darkCount; i++) {
      if ((readDoubleDescriptor(darkList[i]->descs, 
                                pilTrnGetKeyword("ExposureTime"), 
                                &time, comment)) == VM_TRUE) {
        if (time < MIN_DIVISOR) {
          error = 1;
          cpl_msg_error(task, 
            "Zero or negative exposure time for dark %d", i + 1);
          cpl_free(exposureTime);
          break;
        }
        else {
          exposureTime[i] = time;
        }
      }
      else {
        error = 1;
        cpl_msg_error(task, "Cannot read %s from dark frame %d", 
                    pilTrnGetKeyword("ExposureTime"), i + 1);
        cpl_free(exposureTime);
        break;
      }
    }
  }
  else 
    cpl_msg_error(task, "Problems with memory allocation");

  if (error) {
    for (j = 0; j < darkCount; j++)
      deleteImage(darkList[j]);
    cpl_free(darkList);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Get gain factors (e-/ADU)
  */

  error = 1;

  if ((darkGain = (float *)cpl_calloc(darkCount, sizeof(float)))) {
    error = 0;
    for (i = 0; i < darkCount; i++) {
      gain = getMeanGainFactor(darkList[i]);
      if (gain > MIN_DIVISOR) {
        darkGain[i] = gain;
      }
      else {
        error = 1;
        cpl_msg_error(task, "Wrong or no gain factor for dark %d", i + 1);
        cpl_free(darkGain);
        break;
      }
    }
  }
  else 
    cpl_msg_error(task, "Problems with memory allocation");

  if (error) {
    for (j = 0; j < darkCount; j++)
      deleteImage(darkList[j]);
    cpl_free(darkList);
    cpl_free(exposureTime);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Noise contributions for each input dark frame must be estimated
  * before any processing (as bias subtraction, normalization, etc.) 
  * is applied.
  */

  if ((darkRon = (float *)cpl_calloc(darkCount, sizeof(float)))) {
    for (i = 0; i < darkCount; i++) {
      darkRon[i] = computeAverageRon(darkList[i]);
    }
  }
  else {
    cpl_msg_error(task, "Problems with memory allocation");

    for (i = 0; i < darkCount; i++)
      deleteImage(darkList[i]);
    cpl_free(darkList);
    cpl_free(exposureTime);
    cpl_free(darkGain);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  if ((darkNoise = (float *)cpl_calloc(darkCount, sizeof(float)))) {
    for (i = 0; i < darkCount; i++) {
      darkNoise[i] = evaluateAverageNoise(darkList[i], 
                     darkRon[i], darkGain[i]);
    }
  }
  else {
    cpl_msg_error(task, "Problems with memory allocation");

    for (i = 0; i < darkCount; i++)
      deleteImage(darkList[i]);
    cpl_free(darkList);
    cpl_free(exposureTime);
    cpl_free(darkGain);
    cpl_free(darkRon);
    deleteImage(mBiasImage);
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

  for (i = 0; i < darkCount; i++) {
    if (VmSubBias(darkList[i], mBiasImage, biasMethod) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot remove bias from input dark %d", i + 1);
      for (i = 0; i < darkCount; i++)
        deleteImage(darkList[i]);
      cpl_free(darkList);
      cpl_free(exposureTime);
      cpl_free(darkGain);
      cpl_free(darkRon);
      cpl_free(darkNoise);
      deleteImage(mBiasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  deleteImage(mBiasImage);


 /*
  * Now dark levels can be computed
  */

  if ((darkLevel = (float *)cpl_calloc(darkCount, sizeof(float)))) {
    for (i = 0; i < darkCount; i++) {
      darkLevel[i] = imageMedian(darkList[i]);
      cpl_msg_info(task, "Level of dark %-d is %10.4f +/- %-.4f ADU",
                                      i + 1, darkLevel[i], darkNoise[i]);
    }
  }
  else {
    cpl_msg_error(task, "Problems with memory allocation");
    for (i = 0; i < darkCount; i++)
      deleteImage(darkList[i]);
    cpl_free(darkList);
    cpl_free(darkGain);
    cpl_free(darkRon);
    cpl_free(darkNoise);
    cpl_free(exposureTime);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * When RON and dark levels are known, but before any image
  * normalization is applied, cosmic ray events might be removed 
  * from all input images.
  */

  if (cleanCosmic) {
    cpl_msg_info(task, "Cleaning cosmic ray events...");
    for (i = 0; i < darkCount; i++) {
      if (VmCosmicClean(darkList[i], ccdTable, 0, darkLevel[i], darkGain[i],
                        darkRon[i], cosmicThreshold, cosmicRatio) 
                        == EXIT_FAILURE) {
        cpl_msg_error(task, 
                    "Cannot clean cosmic rays from raw dark frame %d", i + 1);
        for (i = 0; i < darkCount; i++)
          deleteImage(darkList[i]);
        cpl_free(darkList);
        cpl_free(darkNoise);
        cpl_free(darkGain);
        cpl_free(darkRon);
        cpl_free(darkLevel);
        cpl_free(exposureTime);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }
      cpl_msg_info(task, "  Dark %d of %zd done", i + 1, darkCount);
    }
  }

  cpl_free(darkRon);
  cpl_free(darkGain);


 /*
  * All dark values are normalized to the unit of time.
  * Zero or negative exposures have been already taken care of,
  * when reading them from the header.
  */

  for (i = 0; i < darkCount; i++) {
    constArithLocal(darkList[i], exposureTime[i], VM_OPER_DIV);
    darkLevel[i] /= exposureTime[i];
    darkNoise[i] /= exposureTime[i];
    cpl_msg_info(task, "Flux of dark %-d is %10.4f +/- %-.4f ADU/hour", i + 1, 
               3600 * darkLevel[i], 3600 * darkNoise[i]);
  }

  if (validateFrames) {

   /*
    * If requested, exclude frames with median flux not consistent 
    * with the set, and with illumination pattern not consistent 
    * with the set.
    */

   /*
    * Determine the typical dark flux, and associate to it
    * the max noise level found. This needs to be refined
    * in future, modifying the applyListSelection function
    * to keep into account different noise levels in the
    * frames to select.
    */

    if (darkCount < 3) {
      medianLevel = computeAverageFloat(darkLevel, darkCount);
    }
    else {
      medianLevel = medianPixelvalue(darkLevel, darkCount);
    }

    maxNoise = darkNoise[0];
    for (i = 1; i < darkCount; i++) {
      if (maxNoise < darkNoise[i]) 
        maxNoise = darkNoise[i];
    }

   /*
    * Remove from the input list all dark frames which do not have
    * flux comparable with the typical flux. The selection criteria
    * is bound to the max noise found among all images.
    */

    levelThreshold = levelTolerance * maxNoise;

    cpl_msg_info(task, "Selecting frames with consistent dark level ...");
    cpl_msg_info(task, 
               "Valid range for mean dark flux: %10.4f +/- %-.4f ADU/hour",
               3600 * medianLevel, 3600 * levelThreshold);

    lowerThreshold = medianLevel - levelThreshold;
    upperThreshold = medianLevel + levelThreshold;

    original = (VimosImage **)cpl_calloc(darkCount, sizeof(VimosImage *));
    if (original) {
      for (i = 0; i < darkCount; i++)
        original[i] = darkList[i];     /* Hold original image sorting */
    }

   /*
    * Good images are moved to the beginning of the darkList:
    */

    goodFrames = applyListSelection(darkList, darkLevel, darkCount,
                                    lowerThreshold, upperThreshold, 1);

    if (goodFrames < darkCount) {

     /*
      * If some images were filtered, then associated array of values
      * must be resorted accordingly
      */

      error = 1;
      if (remapFloatsLikeImages(original, darkList, darkNoise, darkCount) 
                                                          == EXIT_SUCCESS) {
        if (remapFloatsLikeImages(original, darkList, exposureTime, darkCount) 
                                                          == EXIT_SUCCESS) {
          error = 0;
        }
      }

      if (error) {
        cpl_msg_error(task, "Problems remapping arrays after frame selection");
        for (i = 0; i < darkCount; i++)
          deleteImage(darkList[i]);
        cpl_free(darkList);
        cpl_free(original);
        cpl_free(darkNoise);
        cpl_free(darkLevel);
        cpl_free(exposureTime);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }

    }


   /*
    * Check if the remaining input frames have the same intensity
    * distribution, i.e., if the input dark frames show a similar
    * illumination pattern. This makes sense only if there is more
    * than one image left.
    */

    if (goodFrames > 1) {
      cpl_msg_info(task, "Checking for consistent signal distribution ...");

      for (i = 0; i < darkCount; i++)
        original[i] = darkList[i];     /* Hold original image sorting */

      goodFrames = qcSelectConsistentImages(darkList, darkNoise, goodFrames,
                                            patternTolerance);

      if (goodFrames == 0 && pilErrno == P_EGENERIC) {
	cpl_msg_error(task, "Selection of consistent images failed!");
	return EXIT_FAILURE;
      }

      if (remapFloatsLikeImages(original, darkList, exposureTime, darkCount) 
                                                          == EXIT_FAILURE) {
        cpl_msg_error(task, "Problems remapping arrays after frame selection");
        for (i = 0; i < darkCount; i++)
          deleteImage(darkList[i]);
        cpl_free(darkList);
        cpl_free(original);
        cpl_free(darkNoise);
        cpl_free(darkLevel);
        cpl_free(exposureTime);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }
    }

    if (goodFrames) {
      cpl_msg_info(task, 
           "%zd out of %zd consistent dark frames selected from input set", 
           goodFrames, darkCount);
    }

    cpl_free(original);

  }
  else 
    goodFrames = darkCount;

 /*
  * Destroy inconsistent frames from the rearranged list of dark
  * frames (it has no effect if goodFrames == goodFrames), free arrays.
  */

  for (i = goodFrames; i < darkCount; i++)
    deleteImage(darkList[i]);
  cpl_free(darkNoise);
  cpl_free(darkLevel);

  if (!goodFrames) {
    cpl_msg_error(task, "No dark frames with consistent signal "
                "distribution found in input. All frames excluded!");
    cpl_free(darkList);
    cpl_free(exposureTime);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Stack the surviving dark frames, if requested, and if possible.
  */

  if (combMethod == COMB_UNDEF) {
    mDark = duplicateImage(darkList[0]);
  }
  else if (goodFrames >= minFrames) {
    cpl_msg_info(task, "Combining %zd frames with method '%s'", goodFrames,
               combMethodNames[combMethodEntry]);

    if (combMethod == COMB_AVERAGE) {

     /* 
      * Average is weighted on the exposure times - i.e., eliminate
      * normalization if exposure times are different.
      */

      minExposure = maxExposure = exposureTime[0];
      for (i = 1; i < goodFrames; i++) {
        if (exposureTime[i] < minExposure)
          minExposure = exposureTime[i];
        if (exposureTime[i] > maxExposure)
          maxExposure = exposureTime[i];
      }

      if ((maxExposure - minExposure) / minExposure > SOME_PERCENT) {

       /*
        * Exposure times are different - eliminate normalization
        */

        for (i = 0; i < goodFrames; i++) {
          constArithLocal(darkList[i], exposureTime[i], VM_OPER_MUL);
        }
      }

      mDark = frCombAverage(darkList, goodFrames);

      if ((maxExposure - minExposure) / minExposure > SOME_PERCENT) {

       /*
        * Normalize result by the average exposure time
        */

        constArithLocal(mDark, computeAverageFloat(exposureTime, goodFrames),
                                                                VM_OPER_DIV);
      }
    }
    else {
      mDark = frComb(darkList, goodFrames, combMethod, &combParameter, 0);
    }
  }
  else {
    cpl_msg_error(task, "Not enough good frames for requested stacking "
                "method '%s'!", combMethodNames[combMethodEntry]);
    for (i = 0; i < goodFrames; i++)
      deleteImage(darkList[i]);
    cpl_free(darkList);
    cpl_free(exposureTime);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Compute total exposure time of surviving frames
  */

  for (i = 0; i < goodFrames; i++)
    sumExposures += exposureTime[i];

  cpl_free(exposureTime);

  if (!mDark) {
    cpl_msg_error(task, 
       "Stacking of raw dark frames failed! No master dark was created!");
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


  /*
   * Create the master bias header
   */

  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               pilTrnGetKeyword("MjdObs"), NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               pilTrnGetKeyword("DateObs"), NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               pilTrnGetKeyword("ArchiveFile"), NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               "^C[A-Z]*[1,2]", NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               "^ESO OBS (DID|PROG ID)", NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               "^ESO OBS ID", NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               "^ESO TPL [.DINPSV.]", NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               "^ESO INS (DID|MODE)", NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               "^ESO DET [A-Z].*", NULL);
  vimosDscCopy(&mDark->descs, darkList[0]->descs,
               "^ESO OCS (DID|CON QUAD)", NULL);


  if (cleanBadPixel) {
    cpl_msg_info(task, "Cleaning bad pixels from result frame");
    if (EXIT_FAILURE == cleanBadPixels(mDark, ccdTable, 0)) {
      cpl_msg_error(task, "Cannot clean bad pixels");
      deleteImage(mDark);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

 /*
  * Apply the quality control check to the created master dark.
  */

  if (applyQc) {
    cpl_msg_info(task, "Quality check on created master dark...");
    if (qcCheckDarkLevel(mDark, ccdTable, maxDeviation, 1, 1)
        == EXIT_FAILURE) {
      cpl_msg_error(task, "Quality check on master dark failed!");

      deleteCcdTable(ccdTable);
      deleteImage(mDark);

      return EXIT_FAILURE;
    }
  }

  deleteCcdTable(ccdTable);

  if (computeQC) {

    cpl_msg_info(task, "Computing QC1 parameters...");

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      /*
       *  Currently the central 1600x1800 pixels are used
       */

      int          winSizeX = 1600;
      int          winSizeY = 1800;
      int          winStartX = (mDark->xlen - winSizeX) / 2;
      int          winStartY = (mDark->ylen - winSizeY) / 2;
      unsigned int npix = winSizeX * winSizeY;

      double       value, medianValue, meanValue, rmsValue;

      float       *sample = extractFloatImage(mDark->data,
                            mDark->xlen, mDark->ylen,
                            winStartX, winStartY, winSizeX, winSizeY);

      if (!sample) {
        cpl_msg_error(task, "Memory allocation!");
        deleteImage(mDark);
        return EXIT_FAILURE;
      }


      /*
       * Write PRO.CATG, ARCFILE, TPL ID, and OCS CON QUAD, to QC1 group.
       */

      pilQcWriteString("PRO.CATG", mdarkTag, "Product category");

      qcCopyValue(darkList[0]->descs, pilTrnGetKeyword("ArchiveFile"),
                  NULL, "Archive File Name");

      qcCopyValue(darkList[0]->descs, pilTrnGetKeyword("TplId"),
                  NULL, "Template signature ID");

      qcCopyValue(darkList[0]->descs, pilTrnGetKeyword("Quadrant"),
                  NULL, "Quadrant");


      /* QC.DARK.MASTER.MEAN */

      meanValue = computeAverageFloat(sample, npix);

      qcWriteValueDouble(mDark->descs, meanValue, "QC.DARK.MASTER.MEAN", 
                         "ADU/s", "Mean master dark level");

      /* QC.DARK.MASTER.RMS */

      value = 0.0;
      for (i = 0; i < npix; i++)
        value += (sample[i] - meanValue) * (sample[i] - meanValue);

      value /= npix;
      rmsValue = sqrt(value);

      qcWriteValueDouble(mDark->descs, rmsValue, "QC.DARK.MASTER.RMS", 
                         "ADU/s", "RMS of master dark");

      /* QC.DARK.MASTER.MEDIAN */

      medianValue = medianPixelvalue(sample, npix);

      qcWriteValueDouble(mDark->descs, medianValue, "QC.DARK.MASTER.MEDIAN", 
                         "ADU/s", "Median master dark level");

      cpl_free(sample);

      /* QC.DARK.CURRENT */

      gain = getMeanGainFactor(mDark);

      value = 3600 * gain * medianValue;

      qcWriteValueDouble(mDark->descs, value, "QC.DARK.CURRENT",
                         "e-/pix/h", "Dark current");

      /* QC.DARK.CURRENT.RMS */

      value = 3600 * gain * rmsValue;

      qcWriteValueDouble(mDark->descs, value, "QC.DARK.CURRENT.RMS",
                         "e-/pix/h", "Dark current RMS");

      if (pilQcGroupEnd() == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");
    }
    else
      cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

  } /* End of QC1 computation. */


  for (i = 0; i < goodFrames; i++)
    deleteImage(darkList[i]);
  cpl_free(darkList);


 /*
  * Update the product header
  */

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mDark->descs), 
                                     pilTrnGetKeyword("ExposureTime"),
                                     1.,
                                     pilTrnGetComment("ExposureTime"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mDark->descs), 
                                     pilTrnGetKeyword("DataMin"),
                                     imageMinimum(mDark), 
                                     pilTrnGetComment("DataMin"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mDark->descs), 
                                     pilTrnGetKeyword("DataMax"),
                                     imageMaximum(mDark), 
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mDark->descs), 
                                     pilTrnGetKeyword("SummedExposureTime"),
                                     sumExposures, 
                                     pilTrnGetComment("SummedExposureTime"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertIntDescriptor(&(mDark->descs), 
                                  pilTrnGetKeyword("NFramesCombined"),
                                  goodFrames,
                                  pilTrnGetComment("NFramesCombined"),
                                  "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mDark->descs), 
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(mDark), 
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mDark->descs), 
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(mDark), 
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mDark->descs), 
                                     pilTrnGetKeyword("DataMean"), 
                                     imageMean(mDark), 
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    deleteImage(mDark);
    return EXIT_FAILURE;
  }


  /*
   * Create the product file on disk, set the product attributes and
   * update the set of frames.
   */

  vmstrlower(strcpy(masterDarkName, mdarkTag));
  strcat(masterDarkName, ".fits");

  if (createFitsImage(masterDarkName, mDark, mdarkTag)) {
    outputFrame = newPilFrame(masterDarkName, mdarkTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", masterDarkName);
    deleteImage(mDark);
    return EXIT_FAILURE;
  }


  /*
   * Cleanup.
   */

  deleteImage(mDark);

  return EXIT_SUCCESS;
}


/*
 * Build table of contents, i.e. the list of available plugins, for
 * this module. This function is exported.
 */

cxint
cpl_plugin_get_info(cpl_pluginlist *list)
{

    cpl_recipe *recipe = cpl_calloc(1, sizeof *recipe);
    cpl_plugin *plugin = &recipe->interface;


    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmdark",

    "Create a master dark from set of raw dark frames.",

    "This recipe is used to create a master dark frame from a set of raw\n"
    "dark frames.\n\n"
    "Input files:\n\n"
    "  DO category:  Type:       Explanation:     Required:\n"
    "  DARK          Raw         Dark exposure       Y\n"
    "  MASTER_BIAS   Calib       Master bias         Y\n"
    "  CCD_TABLE     Calib       Bad pixel table     .\n\n"
    "Output files:\n\n"
    "  DO category:  Data type:  Explanation:\n"
    "  MASTER_DARK   FITS image  Master dark\n\n"
    "A CCD table must be specified only if a bad pixel cleaning is\n"
    "requested.\n\n"
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

                    vmdark_create,
                    vmdark_exec,
                    vmdark_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
