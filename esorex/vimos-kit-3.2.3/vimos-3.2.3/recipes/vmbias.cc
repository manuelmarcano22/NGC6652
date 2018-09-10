/* $Id: vmbias.c,v 1.8 2013-07-11 11:44:09 cgarcia Exp $
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-07-11 11:44:09 $
 * $Revision: 1.8 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <vector>
#include <sstream>
#include <cmath>
#include <string.h>
#include <math.h>
#include <memory>

#include <cpl.h>

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
#include "vimos_detmodel.h"
#include "vimos_overscan.h"

#include "statistics.h"
#include "fiera_config.h"

static cxint vmbias(PilSetOfFrames *sof, cpl_frameset *frameset,
                    cpl_parameterlist * parlist);

void vimos_bias_compute_ron
(const std::vector<mosca::image>& raw_biases, mosca::ccd_config& bias_ccd_config);

void vimos_bias_compute_ron_single
(const cpl_image * raw_bias, mosca::ccd_config& bias_ccd_config);

std::auto_ptr<mosca::image> vimos_bias_stack
(std::vector<mosca::image>& images, 
 CombMethod method, CombParameters combParameter);

void vimos_bias_write_qc(cpl_propertylist * qc_param, 
                         const mosca::image& master_bias,
                         const std::vector<mosca::image>& biases,
                         const mosca::image& first_raw_bias,
                         vimos_preoverscan& scan_corr,
                         size_t goodFrames, mosca::ccd_config& bias_ccd_config);

/*
 * Definition of the label strings for all methods the recipe function
 * supports for combining frames and their associated method code.
 */

static const cxchar *methodNames[] = {
  "Auto",
  "Ksigma",
  "MinMax",
  "Median",
  "Average"
};

static const CombMethod methods[] = {
  COMB_AUTO,
  COMB_KSIGMA,
  COMB_REJECT,
  COMB_MEDIAN,
  COMB_AVERAGE
};


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it availble to the application using the interface.
 */

static cxint
vmbias_create(cpl_plugin *plugin)
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
                               "Frame combination method",
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


    p = cpl_parameter_new_value("vimos.Parameters.stacking.minmax.minimum",
                                CPL_TYPE_INT,
                                "Number of lowest rejected values for "
                                "rejection method",
                                "vimos.Parameters",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MinRejection");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MinRejection");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.stacking.minmax.maximum",
                                CPL_TYPE_INT,
                                "Number of highest rejected values for "
                                "rejection method",
                                "vimos.Parameters",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MaxRejection");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MaxRejection");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.bias.frames.validate",
                                CPL_TYPE_BOOL,
                                "Consistency check on raw bias frames",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ValidateFrames");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ValidateFrames");
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("vimos.Parameters.bias.tolerance.level",
                                CPL_TYPE_DOUBLE,
                                "Threshold for level consistency",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LevelTolerance");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LevelTolerance");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.bias.tolerance.pattern",
                                CPL_TYPE_DOUBLE,
                                "Threshold for flux pattern consistency",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PatternTolerance");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PatternTolerance");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

    p = cpl_parameter_new_value("vimos.Parameters.bias.overscan.remove",
                                CPL_TYPE_BOOL,
                                "Remove overscan regions from master bias",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "RemoveOverscan");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "RemoveOverscan");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.bias.badpixel.clean",
                                CPL_TYPE_BOOL,
                                "Bad pixel correction on master bias",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.bias.cosmics.clean",
                                CPL_TYPE_BOOL,
                                "Cosmic ray removal from each raw bias",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanCosmic");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanCosmic");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.bias.quality.enable",
                                CPL_TYPE_BOOL,
                                "Compute QC1 parameters",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ComputeQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ComputeQC");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.bias.quality.apply",
                                CPL_TYPE_BOOL,
                                "Quality control of master bias",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyQC");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.bias.quality.maxoffset",
                                CPL_TYPE_DOUBLE,
                                "Maximum allowed deviation from nominal "
                                "bias level",
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
vmbias_exec(cpl_plugin *plugin)
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

    if (vmbias(sof, recipe->frames, recipe->parameters) == EXIT_SUCCESS) {

        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmbias");
        
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
vmbias_destroy(cpl_plugin *plugin)
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


/* 
 * @brief
 *   Create a master bias from set of raw bias frames.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the function returns @c EXIT_FAILURE.
 *
 * @param sof   Set of frames containing the references to
 *              raw bias images and (optionally) a CCD Table.
 *
 * The recipe function creates a master bias frame from a set of
 * raw bias frames. The raw frames are passed via the set of
 * frames @b sof. On successful termination the created master
 * bias is added to this set of frames.
 * 
 * Control options and additional parameters are read from the recipe
 * configuration database. The recipe function accepts the following
 * task parameters:
 * @li StackMethod
 * @li KSigmaLow
 * @li KSigmaHigh
 * @li MinRejection
 * @li MaxRejection
 * @li ValidateFrames
 * @li LevelTolerance
 * @li PatternTolerance
 * @li RemoveOverscan
 * @li CorrectBadPixel
 * @li CorrectCosmic
 * @li ApplyQC
 * @li MaxDeviation
 * If any of these task parameters is not set in the recipe configuration
 * database the recipe function uses the builtin defaults for these
 * task parameters.
 * 
 * The master bias is created by combining the raw bias frames using
 * one of the following combination methods:
 * @li Average of frames
 * @li Median stacking
 * @li Minimum-Maximum rejection
 * @li Kappa-Sigma clipping
 * @li Automatic (not yet available)
 * The default combination method is 'Automatic'. This method determines
 * the number of available raw bias frames and selects the method
 * giving the best possible result (currently this is not implemented,
 * it falls back to average).
 * 
 * The recipe function can be configured to ignore input frames showing
 * a bias level and/or image pattern that is inconsistent with the
 * majority of the input frames. This consistency check is controlled
 * by the task parameter @b ValidateFrames. The default is to perform
 * a validation of the input frames. Bias frames are ignored if their bias
 * level is not within @b LevelTolerance sigma of the median bias level
 * measured across all raw bias frames. The variance used for this bias
 * selection is computed as the average standard deviation from the median
 * bias level. A bias frame is also excluded if the image pattern is
 * different compared to the other bias frames. This is measured by comparing
 * the deviation from 0 of the median pixel value of the difference image
 * with @b PatternTolerance times the expected noise.
 * 
 * The function allows to remove or keep the pre- and overscan regions
 * in the created master bias frame. The removal of the pre- and overscan
 * regions is controlled by the task parameter @b RemoveOverscan.
 * The default is to remove the pre- and overscan areas.
 * 
 * Optionally the input raw bias frames may be corrected for cosmic
 * ray events. The cosmic ray removal is controlled by the task 
 * parameter @b CorrectCosmic. In case a CCD_TABLE is also given, bad
 * pixels will not be used in computing the interpolated values to replace
 * the cosmic rays events with. The default is not to correct for cosmic
 * ray events.
 * 
 * Optionally the created master bias frame may be corrected for bad
 * pixels. If this option is turned on, the recipe expects to find
 * a CCD_TABLE in the input set of frames. The bad pixel correction
 * is controlled by the task parameter @b CorrectBadPixel. The
 * default is not to correct for bad pixels.
 * 
 * If the task parameter @b ApplyQC is set to @b true the recipe
 * function applies a simple quality control task to the created
 * master bias frame. The check compares the bias level of the master
 * bias with the nominal bias value taken from the header of a reference
 * master bias which is looked up in the set of frames @em sof. If the
 * offset from the nominal level is larger than @b MaxDeviation sigma a
 * warning is issued, but the function does not fail.
 * If this option is turned on the recipe expects to find a master bias in
 * input, or it will return an error.
 */

static cxint
vmbias(PilSetOfFrames *sof, cpl_frameset *frameset,
       cpl_parameterlist * parlist)
{

    const char task[] = "vmbias";
    char        version[80];

    const char *biasTag = pilTrnGetCategory("Bias");
    const char *mbiasTag = pilTrnGetCategory("MasterBias");

    const char *methodTag = 0;
    char mbiasName[PATHNAME_MAX + 1];

    size_t  biasCount, goodFrames;
    size_t  minFrames;

    unsigned int i, j;
    unsigned int validateFrames, removeOverscan;
    unsigned int cleanBadPixel, cleanCosmic, computeQC, applyQc;
    unsigned int singleFrames;

    int methodEntry;

    float *biasLevel = 0;
    float *biasNoise = 0;

    double medianLevel, meanNoise;
    double levelTolerance, levelThreshold, lowerThreshold, upperThreshold;
    double patternTolerance;
    double maxDeviation;

    PilFrame *currFrame, *ccdFrame, *rbiasFrame, *productFrame;

    CombMethod method = COMB_UNDEF;
    CombParameters combParameter;

    VimosImage *mBias = 0;
    VimosImage *rBias = 0;
    VimosImage **biasList, **refList;

    VimosTable *ccdTable = 0;


    snprintf(version, 80, "%s-%s", PACKAGE, PACKAGE_VERSION);


    //Set frame group
    vimos_dfs_set_groups(frameset);
    
    /*
     * Get task parameters from the recipe database
     */

    /*
     * Determine the frame stacking method first and all
     * method dependent parameters.
     */
  
    methodTag = pilDfsDbGetString("Parameters", "StackMethod");
    methodEntry = strselect(methodTag, methodNames, PIL_N_ELEMENTS(methods));

    if (methodEntry < 0) {
        cpl_msg_error(task, "%s: Invalid frame combination method.", methodTag);
        return EXIT_FAILURE;
    }
    else {
        method = methods[methodEntry];
    }

    switch (method) {
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
        method = COMB_AVERAGE;
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

    validateFrames = pilDfsDbGetBool("Parameters", "ValidateFrames", 0);


    /*
     * Check if the overscan areas should be removed after stacking the
     * frames.
     */

    removeOverscan = pilDfsDbGetBool("Parameters", "RemoveOverscan", 1);


    /*
     * Check if the bad pixels should be corrected.
     */

    cleanBadPixel = pilDfsDbGetBool("Parameters", "CleanBadPixel", 0);


    /*
     * Check if the cosmic rays events should be corrected.
     */

    cleanCosmic = pilDfsDbGetBool("Parameters", "CleanCosmic", 0);


    /*
     * Check if QC1 parameters should be computed
     */

    computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


    /*
     * Check if the bias quality control should be applied to the
     * final product.
     */

    applyQc = pilDfsDbGetBool("Parameters", "ApplyQC", 0);
    maxDeviation = pilDfsDbGetDouble("Parameters", "MaxDeviation", 3.0);


    /* 
     * Make sure that there are enough raw data frames in the
     * input set.
     */

    biasCount = (int)pilSofFrameCount(sof, biasTag);

    if (biasCount < minFrames) {
        if (biasCount == 1 && singleFrames) {
            method = COMB_UNDEF;
            methodTag = "Undefined";
            if (computeQC) {
                cpl_msg_warning(task, "QC1 parameters are not computed: at "
                              "least two bias exposures would be required.");
                computeQC = 0;
            }
        }
        else {
            cpl_msg_error(task, "Not enough raw bias frames in input for "
                        "stacking method '%s'!", methodTag);
            return EXIT_FAILURE;
        }
    }


    /*
     * If the bad pixel correction is enabled, a bad pixel table is required
     * in the input set of frames. If no bad pixel table is present this is an
     * error. If cleaning of cosmic rays is enabled, the bad pixel table will
     * be used, if present, to avoid bad pixels while smoothing away cosmic
     * ray events. In this case, a missing bad pixel table is not an error.
     */


    ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable"));
    if (ccdFrame) 
        pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

    if (cleanBadPixel || cleanCosmic) {
        if (!ccdFrame) {
            if (cleanBadPixel) {
                cpl_msg_error(task, "Bad pixel cleaning requires a CCD "
                            "table in input!");
                return EXIT_FAILURE;
            }
        }
        else {
            cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));

            if (!(ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
                cpl_msg_error(task, "Cannot load CCD table %s!",
                            pilFrmGetName(ccdFrame));
                return EXIT_FAILURE;
            }
            else 
                closeFitsTable(ccdTable, 0);
        }
    }


    /*
     * If the quality control check is requested a reference master bias is
     * expected in the input set of frames. If it is not found it is an error.
     */

    rbiasFrame = pilSofLookup(sof, pilTrnGetCategory("MasterBias"));
    if (rbiasFrame)
        pilFrmSetType(rbiasFrame, PIL_FRAME_TYPE_CALIB);

    if (applyQc) {

        if (rbiasFrame == NULL) {
            cpl_msg_error(task, "Product frame quality control requires "
                        "a master bias in input!");
            return EXIT_FAILURE;
        }
        else {
            cpl_msg_debug(task, "Reference bias is %s",
                        pilFrmGetName(rbiasFrame));

            if (!(rBias = openOldFitsFile(pilFrmGetName(rbiasFrame), 0, 0))) {
                cpl_msg_error(task, "Cannot load reference bias %s!",
                            pilFrmGetName(rbiasFrame));
                return EXIT_FAILURE;
            }
            else {
                if (loadFitsHeader(rBias) == VM_FALSE) {
                    cpl_msg_error(task, "Cannot load reference bias header");

                    closeFitsImage(rBias, 0);
                    deleteImage(rBias);

                    if (ccdTable)
                        deleteTable(ccdTable);
                }
                else
                    closeFitsImage(rBias, 0);
            }
        }
    }


    /*
     * Load the raw bias frames.
     */

    biasList = (VimosImage **)cpl_calloc(biasCount, sizeof(VimosImage *));
    if (biasList == NULL) {
        cpl_msg_error(task, "Not enough memory!");
        deleteTable(ccdTable);

        return EXIT_FAILURE;
    }

    currFrame = pilSofLookupNext(sof, biasTag);

    for (i = 0; i < biasCount; i++) {

        biasList[i] = openOldFitsFile(pilFrmGetName(currFrame), 1, 0);

        if (biasList[i] == NULL) {
            cpl_msg_error(task, "Cannot load bias frame %d", i + 1);

            for (j = 0; j < i; j++)
                deleteImage(biasList[j]);
            cpl_free(biasList);

            deleteTable(ccdTable);
            return EXIT_FAILURE;
        }
        else {
            pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
            closeFitsImage(biasList[i], 0);
            currFrame = pilSofLookupNext(sof, 0);
        }
    }



    /*
     * Do a preselection of the raw bias frames based on their average
     * median intensity level. A bias frames whose level deviates too
     * much from the average level are removed from the input and are
     * not treated in all following processing steps.
     */

    if (validateFrames && biasCount < 2) {
        validateFrames = 0;
        cpl_msg_warning(task, "Too few bias frames (%zd) in input. Skipping "
                      "frame selection task!", biasCount);
    }

    /*
     * Initially assume that all frames are valid. This is needed to
     * make sure that goodFrames is setup correctly if the tests
     * are skipped.
     */

    goodFrames = biasCount;

    if (validateFrames || cleanCosmic) {

        if (!(biasLevel = (float *)cpl_calloc(biasCount, sizeof(float)))) {
            cpl_msg_error(task, "Not enought memory!");

            for (i = 0; i < biasCount; i++)
                deleteImage(biasList[i]);

            cpl_free(biasList);

            deleteTable(ccdTable);
            return EXIT_FAILURE;
        }

        if (!(biasNoise = (float *)cpl_calloc(biasCount, sizeof(float)))) {
            cpl_msg_error(task, "Not enough memory!");

            for (i = 0; i < biasCount; i++) {
                deleteImage(biasList[i]);
            }

            cpl_free(biasList);
            cpl_free(biasLevel);

            deleteTable(ccdTable);
            return EXIT_FAILURE;
        }


        /*
         * Bias mean level and read out noise must be evaluated if either 
         * cosmic rays cleaning or frame validation was requested.
         */

        cpl_msg_info(task, "Computing bias level(s) ...");

        for (i = 0; i < biasCount; i++) {
            biasLevel[i] = imageMean(biasList[i]);
            biasNoise[i] = computeAverageRon(biasList[i]);

            if (biasNoise[i] < 0.) {
                cpl_msg_warning(task, "Get bias %d RON from keyword header",
                              i + 1);
                biasNoise[i] = getAverageRon(biasList[i]);

                if (biasNoise[i] < 0.) {

                    cpl_msg_error(task, "Cannot compute bias RON!");

                    for (i = 0; i < biasCount; i++) {
                        deleteImage(biasList[i]);
                    }

                    cpl_free(biasList);
                    cpl_free(biasLevel);

                    deleteTable(ccdTable);

                    return EXIT_FAILURE;
                }

            }

            cpl_msg_info(task, "Level of bias %-d is %10.4f +/- "
                       "%-.4f ADU", i + 1, biasLevel[i], biasNoise[i]);
        }


        /*
         * Remove cosmic ray hits from all input frames
         */

        if (cleanCosmic) {
            cpl_msg_info(task, "Cleaning cosmic ray events...");

            for (i = 0; i < biasCount; i++) {
                if (VmCosmicClean(biasList[i], ccdTable, 0, biasLevel[i], 1.,
                                  biasNoise[i], -1., -1.) == EXIT_FAILURE) {
                    cpl_msg_error(task, "Cannot clean cosmic ray hits from "
                                "raw bias frame %d", i + 1);

                    for (i = 0; i < biasCount; i++)
                        deleteImage(biasList[i]);

                    cpl_free(biasList);
                    cpl_free(biasNoise);
                    cpl_free(biasLevel);

                    deleteTable(ccdTable);
                    return EXIT_FAILURE;
                }
                cpl_msg_info(task, "  Bias %d of %zd done", i + 1, biasCount);
            }
        }


        /*
         * Select valid frames only for the frame combination.
         */

        if (validateFrames) {

            cpl_msg_info(task, "Checking bias levels for consistency ...");

            if (biasCount < 3) {
                medianLevel = computeAverageFloat(biasLevel, biasCount);
            }
            else {
                medianLevel = medianPixelvalue(biasLevel, biasCount);
            }

            meanNoise = computeAverageFloat(biasNoise, biasCount);
            cpl_msg_info(task, "Mean bias noise: %-.4f adu", meanNoise);


            /*  
             * Remove all bias frames from the input list which do not have
             * a bias level comparable to the median bias level of all bias
             * frames in the list. The selection criteria is bound to the 
             * mean noise of the frames.
             */

            levelTolerance = pilDfsDbGetDouble("Parameters", "LevelTolerance",
                                               3.0);
            levelThreshold = levelTolerance * meanNoise;

            cpl_msg_info(task, "Selecting frames with consistent bias "
                       "level ...");
            cpl_msg_info(task, "Valid range for bias levels: %-10.4f +/- "
                       "%-.4f adu", medianLevel, levelThreshold);

            /*
             * Make sure that bias frames having a level smaller than the
             * mean noise are excluded.
             */

            upperThreshold = medianLevel + levelThreshold;
            lowerThreshold = medianLevel - levelThreshold;

            if (lowerThreshold < 0) {
                lowerThreshold = meanNoise;
                cpl_msg_warning(task, "Lower limit bias level would be "
                              "below 0. Lower limit reset to %-.4f adu",
                              lowerThreshold);
            }

            /*
             * The original sorting order of the input bias list must be kept.
             * It is needed later on for changing the sorting order of the
             * value arrays associated to the list, like the level and noise
             * array.
             */

            refList = (VimosImage **)cpl_calloc(biasCount,
                                                sizeof(VimosImage *));
            if (!refList) {
                cpl_msg_error(task, "Not enough memory!");

                for (i = 0; i < biasCount; i++)
                    deleteImage(biasList[i]);

                cpl_free(biasList);
                cpl_free(biasNoise);
                cpl_free(biasLevel);

                deleteTable(ccdTable);
                return EXIT_FAILURE;
            }
            else 
                for (i = 0; i < biasCount; i++)
                    refList[i] = biasList[i];


            /*
             * Good images are moved to the beginning of the biasList:
             */

            goodFrames = applyListSelection(biasList, biasLevel, biasCount, 
                                            lowerThreshold, upperThreshold, 1);

            if (goodFrames < biasCount) {
                if (remapFloatsLikeImages(refList, biasList, biasNoise,
                                          biasCount) == EXIT_FAILURE) {
                    cpl_msg_error(task, "Bias frame selection failed!");

                    for (i = 0; i < biasCount; i++)
                        deleteImage(biasList[i]);

                    cpl_free(refList);
                    cpl_free(biasList);
                    cpl_free(biasNoise);
                    cpl_free(biasLevel);

                    deleteTable(ccdTable);
                    return EXIT_FAILURE;
                }
            }

            cpl_free(refList);

            cpl_msg_info(task, "%zd raw bias frames selected.", goodFrames);


            /*
             * Check if the remaining input frames have the same intensity
             * distribution, i.e., if the input bias frames show a similar
             * illumination pattern. This makes sense only if there is more
             * than one image left.
             */

            if (goodFrames > 1) {
                cpl_msg_info(task, "Checking for consistent signal "
                           "distribution ...");

                patternTolerance = pilDfsDbGetDouble("Parameters",
                                                     "PatternTolerance", 3.0);
                goodFrames = qcSelectConsistentImages(biasList, biasNoise,
                                                      goodFrames,
                                                      patternTolerance);

                if (goodFrames == 0 && pilErrno == P_EGENERIC) {
                    cpl_msg_error(task, "Selection of consistent images "
                                "failed!");
                    return EXIT_FAILURE;
                }

                cpl_msg_info(task, "%zd out of %zd consistent bias frames "
                           "selected from input set", goodFrames, biasCount);
            }
        }

        cpl_free(biasLevel);
        cpl_free(biasNoise);


        /*
         * Erase inconsistent frames from the rearranged list of bias
         * frames.
         */

        for (i = goodFrames; i < biasCount; i++)
            deleteImage(biasList[i]);

        if (!goodFrames) {
            cpl_msg_error(task, "No bias frames with consistent signal "
                        "distribution found in input. All frames excluded!");

            cpl_free(biasList);

            deleteTable(ccdTable);
            return EXIT_FAILURE;
        }
    }


    /* Get instrument setting */
    cpl_propertylist * bias_header = 
            cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position(frameset, 0)), 0);
    mosca::fiera_config bias_ccd_config(bias_header);
    if(cpl_error_get_code())
    {
        cpl_msg_error(cpl_func, "Could not get instrument setting");
        return EXIT_FAILURE;
    }
 
    /* Cast the images into mosca::image and create the variance maps*/
    std::vector<mosca::image> bias_images;
    for(size_t i_bias = 0; i_bias < goodFrames; ++i_bias)
    {
        cpl_image * this_bias =
                cpl_image_wrap(biasList[i_bias]->xlen, biasList[i_bias]->ylen, CPL_TYPE_FLOAT, 
                               biasList[i_bias]->data);
        
        /* Compute the RON and update the ccd config just for this image */
        vimos_bias_compute_ron_single(this_bias, bias_ccd_config);
        
        cpl_image * this_bias_var = 
                vimos_image_variance_from_detmodel(this_bias, bias_ccd_config);
        
        //We convert to error
        cpl_image_power(this_bias_var, 0.5);

        //The temporary mosca::image doesn't take over ownership, since
        //the vector push_back will do a deep copy of it.
        bias_images.push_back(mosca::image(this_bias, this_bias_var,
                                           false, mosca::Y_AXIS));
        cpl_image_unwrap(this_bias);
        cpl_image_delete(this_bias_var);
        deleteImage(biasList[i_bias]);
    }
    cpl_free(biasList);

    /* Compute the RON and update the ccd config */
    vimos_bias_compute_ron(bias_images, bias_ccd_config);

    /* 
     * Subtract overscan from the raw biases
     */
    vimos_preoverscan scan_corr;
    std::vector<mosca::image> bias_ovs = 
            scan_corr.subtract_overscan(bias_images, bias_ccd_config);
    if(cpl_error_get_code())
    {
        cpl_msg_error(cpl_func, "Cannot subtract pre/overscan");
        return EXIT_FAILURE;
    }
    //Leave the first raw bias, since it is used to compute some QC values.
    bias_images.erase(bias_images.begin()+1, bias_images.end());
   
    /* Trimm pre/overscan */
    std::vector<mosca::image> bias_trimmed;
    if (removeOverscan) {
        bias_trimmed = 
                scan_corr.trimm_preoverscan(bias_ovs, bias_ccd_config);
        if(cpl_error_get_code())
        {
            cpl_msg_error(cpl_func, "Cannot trimm pre/overscan");
            return EXIT_FAILURE;
        }
    }
    else
        bias_trimmed = bias_ovs;
    bias_ovs.clear();

    /* Stack */
    cpl_msg_info(task, "Combining %zd frames with method '%s'", goodFrames,
                methodNames[methodEntry]);
    std::auto_ptr<mosca::image> master_bias = vimos_bias_stack(bias_trimmed, method, combParameter);
    if(cpl_error_get_code() || master_bias.get() == NULL)
    {
        cpl_msg_error(cpl_func, "Cannot stack biases");
        return EXIT_FAILURE;
    }
    
    /*
     * Bad pixel cleaning of the result
     */
    mBias = newImageAndAlloc(master_bias->size_x(), master_bias->size_y());
    mosca::image master_bias_f(cpl_image_cast(master_bias->get_cpl_image(), CPL_TYPE_FLOAT), 
                               cpl_image_cast(master_bias->get_cpl_image_err(), CPL_TYPE_FLOAT),
                               true);
    //Dirty casting into VimosImage
    float * data_save = mBias->data; 
    mBias->data = master_bias_f.get_data<float>();
    
    if (cleanBadPixel) {
        cpl_msg_info(task, "Cleaning bad pixels on result frame ...");

        if (cleanBadPixels(mBias, ccdTable, 0) == EXIT_FAILURE) {
            cpl_msg_error(task, "Bad pixel cleaning failed!");
            mBias->data = data_save;

            deleteTable(ccdTable);

            return EXIT_FAILURE;
        }
    }
    mBias->data = data_save;
    deleteImage(mBias);

    /* QC */
    cpl_propertylist * qc_param = cpl_propertylist_new();
    vimos_bias_write_qc(qc_param, *master_bias, bias_trimmed, bias_images[0],
                        scan_corr, goodFrames, bias_ccd_config);

    /* Save product */
    //Append the QC param to the raw bias header 
    cpl_propertylist_append(bias_header, qc_param);
    //Update the WCS by removing overscan
    cpl_propertylist_set_double(bias_header, "CRPIX1", 
           cpl_propertylist_get_double(bias_header, "CRPIX1") - 
           bias_ccd_config.validpix_region(0).llx());
    cpl_propertylist_set_double(bias_header, "CRPIX2", 
           cpl_propertylist_get_double(bias_header, "CRPIX2") - 
           bias_ccd_config.validpix_region(0).lly());
    //Finally save it
    if (dfs_save_image(frameset, master_bias_f.get_cpl_image(), mbiasTag, 
                       bias_header, parlist, task, version)) {
        cpl_msg_error(cpl_func, "Cannot save master bias");
        return EXIT_FAILURE;
    }
    cpl_propertylist_delete(bias_header);
    //And save the error extension 
    cpl_propertylist * errbias_hdr = cpl_propertylist_new();
    cpl_propertylist_append_string(errbias_hdr, "EXTNAME", "ERRORMAP");
    if (dfs_save_image_ext(master_bias_f.get_cpl_image_err(), mbiasTag, errbias_hdr)) {
        cpl_msg_error(cpl_func, "Cannot save master bias error");
        return EXIT_FAILURE;
    }
    cpl_propertylist_delete(errbias_hdr);
 

    /*
     * Apply the quality control check to the created master bias
     * For the master bias creation the QC check will only issue
     * a warning but it will not fail. For the check a reference
     * master bias is required.
     */

    if (applyQc) {
        int state;

        state = qcCheckBiasLevel(mBias, rBias, maxDeviation, 1, 1);

        if (state == EXIT_FAILURE) {
            cpl_msg_error(task, "Quality control of master bias failed!");

            deleteTable(ccdTable);

            return EXIT_FAILURE;
        }
    }

    /*
     * Cleanup.
     */

    deleteTable(ccdTable);
 
    return EXIT_SUCCESS;

}


/*
 * Build table of contents, i.e. the list of available plugins, for
 * this module. This function is exported.
 */

cxint
cpl_plugin_get_info(cpl_pluginlist *list)
{

    cpl_recipe *recipe = (cpl_recipe*)cpl_calloc(1, sizeof *recipe);
    cpl_plugin *plugin = &recipe->interface;


    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmbias",

    "Create a master bias from set of raw bias frames.",

    "This recipe is used to create a master bias frame from a set of raw\n"
    "bias frames.\n\n" 
    "Input files:\n\n"
    "  DO category:  Type:       Explanation:     Required:\n"
    "  BIAS          Raw         Bias exposure       Y\n"
    "  CCD_TABLE     Calib       Bad pixel table     .\n\n"
    "Output files:\n\n"
    "  DO category:  Data type:  Explanation:\n"
    "  MASTER_BIAS   FITS image  Master bias\n\n"
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
    "Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA\n",

                    vmbias_create,
                    vmbias_exec,
                    vmbias_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}

void vimos_bias_compute_ron
(const std::vector<mosca::image>& raw_biases, mosca::ccd_config& bias_ccd_config)
{
    for(size_t iport = 0; iport< bias_ccd_config.nports(); ++iport)
    {
        mosca::rect_region os_region = 
                bias_ccd_config.overscan_region(iport).coord_0to1();
        std::vector<double> rons;
        for(size_t ibias = 0; ibias< raw_biases.size(); ++ibias)
        {
            const mosca::image& bias = raw_biases[ibias];
            //These coordinates are in spectral, spatial coordinates
            mosca::image bias_os_im= bias.trim(os_region.lly(),
                    os_region.llx(), os_region.ury(), os_region.urx());
            float * data = bias_os_im.get_data<float>();
            rons.push_back(mosca::robust_variance(data, 
                                                  data + bias_os_im.npix()));
        }
        double ron = std::sqrt(mosca::mean(rons.begin(), rons.end()));
        bias_ccd_config.set_computed_ron(iport, ron);
    }
}

void vimos_bias_compute_ron_single
(const cpl_image * raw_bias, mosca::ccd_config& bias_ccd_config)
{
    for(size_t iport = 0; iport< bias_ccd_config.nports(); ++iport)
    {
        mosca::rect_region os_region = 
                bias_ccd_config.overscan_region(iport).coord_0to1();
        std::vector<double> rons;
        cpl_image * bias_os_im= cpl_image_extract
                (raw_bias, os_region.llx(), os_region.lly(), os_region.urx(), os_region.ury());
        float * data = cpl_image_get_data_float(bias_os_im);
        double ron = std::sqrt(mosca::robust_variance(data, 
                                 data + cpl_image_get_size_x(bias_os_im) *
                                 cpl_image_get_size_y(bias_os_im)));
        bias_ccd_config.set_computed_ron(iport, ron);
        cpl_image_delete(bias_os_im);
    }
}


std::auto_ptr<mosca::image> vimos_bias_stack
(std::vector<mosca::image>& images, 
 CombMethod method, CombParameters combParameter)
{
    std::auto_ptr<mosca::image> master_bias;
    cpl_image * contrib;

    if(images.size() == 0)
        return master_bias;
    
    hdrl_parameter * stackmethod_par = NULL;
    hdrl_imagelist * images_hdrl = NULL;
    hdrl_image * master_bias_hdrl = NULL;

    switch (method) {
    case COMB_AVERAGE :
        stackmethod_par = hdrl_collapse_mean_parameter_create();
        break;
//    case stack_method::WMEAN :
//        stackmethod_par = hdrl_collapse_weighted_mean_parameter_create();
//        break;
    case COMB_MEDIAN : 
        stackmethod_par = hdrl_collapse_median_parameter_create();
        break;
    case COMB_REJECT :
        stackmethod_par = hdrl_collapse_minmax_parameter_create(
                combParameter.minRejection, combParameter.maxRejection);
        break;
    case COMB_KSIGMA : 
        stackmethod_par = hdrl_collapse_sigclip_parameter_create
          (combParameter.kSigmaLow, combParameter.kSigmaHigh, 2);
        break;
    default:
        cpl_msg_error(cpl_func, "Unknown stack method");
        return master_bias;
        break;
    }

    //Transform to HDRL
    hdrl_imagelist * im_list = hdrl_imagelist_new();
    for(size_t idx = 0; idx < images.size(); idx++)
    {
        const mosca::image& im = images[idx];
        hdrl_image * tmp = hdrl_image_create(im.get_cpl_image(),
                                             im.get_cpl_image_err());
                
        hdrl_imagelist_set(im_list, tmp, idx);
    }

    //Do the stacking
    hdrl_imagelist_collapse(im_list, stackmethod_par ,&master_bias_hdrl, 
                            &contrib);
    cpl_image_delete(contrib);
    hdrl_imagelist_delete(im_list);

    //Transform back to mosca. 
    cpl_image * img =  cpl_image_duplicate(hdrl_image_get_image(master_bias_hdrl));
    cpl_image * err =  cpl_image_duplicate(hdrl_image_get_error(master_bias_hdrl));
    
    //The dispersion axis is irrelevant in this case, but we put the real one
    master_bias.reset(new mosca::image(img, err, true, mosca::Y_AXIS));
    
    if(stackmethod_par != NULL)
        hdrl_parameter_destroy(stackmethod_par);
    if(images_hdrl != NULL)
        hdrl_imagelist_delete(images_hdrl);
    if(master_bias_hdrl != NULL)
        hdrl_image_delete(master_bias_hdrl);

    return master_bias;
}

void vimos_bias_write_qc(cpl_propertylist * qc_param, 
                         const mosca::image& master_bias,
                         const std::vector<mosca::image>& biases,
                         const mosca::image& first_raw_bias,
                         vimos_preoverscan& scan_corr,
                         size_t goodFrames, mosca::ccd_config& bias_ccd_config)
{
    cpl_msg_info(cpl_func, "Computing QC1 parameters...");

    /*
     *  Currently the central 1600x1800 pixels are used
     */

    int winSizeX = 1600;
    int winSizeY = 1800;
    int winStartX = (biases[0].size_x() - winSizeX) / 2;
    int winStartY = (biases[0].size_y() - winSizeY) / 2;
    int good;
    size_t i;

    size_t npix = winSizeX * winSizeY;

    double value, meanValue, ronValue, fpnValue, rmsValue, structValue;
    double threshold;

    float *sample1 = NULL;
    float *sample2 = NULL;
    float *sample3 = NULL;
    double *sample4 = NULL;
    float *sample5 = NULL;
    double *diff = NULL;
    
    cpl_image * sample1_im = cpl_image_extract(biases[0].get_cpl_image(),
            winStartX, winStartY, 
            winStartX + winSizeX, winStartY + winSizeY);
    cpl_image * sample2_im = cpl_image_extract(biases[1].get_cpl_image(),
            winStartX, winStartY, 
            winStartX + winSizeX, winStartY + winSizeY);
    cpl_image * sample5_im = cpl_image_extract(first_raw_bias.get_cpl_image(),
            winStartX, winStartY, 
            winStartX + winSizeX, winStartY + winSizeY);
    
    sample1 = cpl_image_get_data_float(sample1_im);
    sample2 = cpl_image_get_data_float(sample2_im);
    sample5 = cpl_image_get_data_float(sample5_im);

    diff = (double *)cpl_malloc(winSizeX * winSizeY * sizeof(double));

    /*
     *  Compute images difference
     */

    for (i = 0; i < npix; i++)
        diff[i] = sample1[i] - sample2[i];

    /*
     * Compute QC parameters.
     */

    /* QC.BIAS.MEAN */

    meanValue = computeAverageFloat(sample1, winSizeX * winSizeY);

    double meanValueRaw = computeAverageFloat(sample5, winSizeX * winSizeY);

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS MEAN",
                                   meanValueRaw);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS MEAN",
                          "Mean bias level [ADU]");

    /* Old ronValue computed here. It is used for other QC params */

    value = computeVarianceDouble2D(diff, winSizeX, winSizeY);
    ronValue = sqrt(value / 2);

    /* QC.BIAS.FPN */

    cpl_image * sample3_im = cpl_image_extract(biases[1].get_cpl_image(),
            winStartX+10, winStartY+10, 
            winStartX + 10 + winSizeX, winStartY + 10 + winSizeY);
    sample3= cpl_image_get_data_float(sample3_im);

    /*
     *  Compute difference of first bias with shifted second bias.
     */

    for (i = 0; i < npix; i++)
        diff[i] = sample1[i] - sample3[i];


    value = computeVarianceDouble2D(diff, winSizeX, winSizeY);
    if (value / 2 > ronValue * ronValue)
        fpnValue = sqrt(value / 2 - ronValue * ronValue);
    else
        fpnValue = 0.0;

    cpl_propertylist_append_double(qc_param,  "ESO QC BIAS FPN",
                          fpnValue);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS FPN",
                                 "Bias fixed pattern noise [ADU]");

    cpl_free(diff);

    value = 0.0;
    for (i = 0; i < npix; i++)
        value += (sample1[i] - meanValue) * (sample1[i] - meanValue);

    value /= npix;
    rmsValue = sqrt(value);

    cpl_propertylist_append_double(qc_param,  "ESO QC BIAS RMS",
                          rmsValue);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS RMS",
                                 "RMS of bias [ADU]");

    if (value > fpnValue * fpnValue + ronValue * ronValue)
        structValue = sqrt(value - fpnValue * fpnValue -
                ronValue * ronValue);
    else
        structValue = 0;

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS STRUCT",
                                   structValue);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS STRUCT",
                                 "Bias structure [ADU]");

    /* QC.BIAS.MEDIAN */

    value = medianPixelvalue(sample5, winSizeX * winSizeY);

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS MEDIAN",
                                   value);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS MEDIAN",
                                 "Median bias level [ADU]");

    /* QC.BIAS.MASTER.MEAN */

    /*
     * For the master bias the sample window must be recomputed - here 
     * overscans are missing...
     */

    winStartX = (master_bias.size_x()- winSizeX) / 2;
    winStartY = (master_bias.size_y()- winSizeY) / 2;

    cpl_image * sample4_im = cpl_image_extract(master_bias.get_cpl_image(),
            winStartX, winStartY, 
            winStartX + winSizeX, winStartY + winSizeY);
    sample4 = cpl_image_get_data_double(sample4_im);

    meanValue = computeAverageDouble(sample4, winSizeX * winSizeY);

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS MASTER MEAN",
                          meanValue);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS MASTER MEAN",
                                 "Mean master bias level [ADU]");

    /* QC.BIAS.MASTER.RMS */

    meanValue = computeAverageDouble(sample4, winSizeX * winSizeY);

    value = 0.0;
    for (i = 0; i < npix; i++)
        value += (sample4[i] - meanValue) * (sample4[i] - meanValue);

    value /= npix;
    rmsValue = sqrt(value);

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS MASTER RMS",
                                   rmsValue);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS MASTER RMS",
                                 "RMS of master bias [ADU]");

    /* QC.BIAS.MASTER.NOISE  &  QC.BIAS.MASTER.FPN */

    /*
     * This is not yet the fixed pattern noise: the ron contribution
     * must still be evaluated and subtracted.
     */

    fpnValue = computeVarianceDouble2D(sample4, winSizeX, winSizeY);

    /*
     * The expected noise is QC.RON / sqrt(goodFrames). From this
     * we set a tolerance against outlyers equal to 3 times this
     * value. The noise is evaluated from the master central region, 
     * excluding outlyers the (supposedly belonging to the fixed 
     * pattern noise). Along the way, the fixed pattern noise is 
     * also evaluated.
     */

    threshold = 3 * ronValue / sqrt((double)goodFrames);

    /*
     * Here the RON contribution is evaluated:
     */

    meanValue = computeAverageDouble(sample4, winSizeX * winSizeY);

    good = 0;
    value = 0.0;
    for (i = 0; i < npix; i++) {
        if (fabs(sample4[i] - meanValue) < threshold) {
            value += (sample4[i] - meanValue) *
                    (sample4[i] - meanValue);
            good++;
        }
    }

    /*
     * Threshold was clearly underestimated, skip it
     */

    if (good == 0) {
        value = 0.0;
        for (i = 0; i < npix; i++)
            value += (sample4[i] - meanValue) *
            (sample4[i] - meanValue);
        good = npix;
    }

    value /= good;  /* Mean variance excluding outlyers */

    if (fpnValue > value)
        fpnValue = sqrt(fpnValue - value);
    else
        fpnValue = 0.0;

    ronValue = sqrt(value);

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS MASTER NOISE",
                                   ronValue);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS MASTER NOISE",
                                 "Noise of master bias [ADU]");

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS MASTER FPN",
                                   fpnValue);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS MASTER NOISE",
                                 "Fixed pattern noise of master bias [ADU]");

    /* QC.BIAS.MASTER.STRUCT */

    meanValue = computeAverageDouble(sample4, winSizeX * winSizeY);

    value = 0.0;
    for (i = 0; i < npix; i++)
        value += (sample4[i] - meanValue) * (sample4[i] - meanValue);

    value /= npix;

    if (rmsValue * rmsValue > fpnValue * fpnValue +
            ronValue * ronValue) {
        structValue = sqrt(rmsValue * rmsValue -
                fpnValue * fpnValue - ronValue * ronValue);
    }
    else {
        structValue = 0.0;
    }

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS MASTER STRUCT",
                                   structValue);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS MASTER STRUCT",
                                 "Structure of master bias [ADU]");

    /* QC.BIAS.MASTER.MEDIAN */

    value = medianPixelvalueDouble(sample4, winSizeX * winSizeY);

    cpl_propertylist_append_double(qc_param, "ESO QC BIAS MASTER MEDIAN",
                                   value);
    cpl_propertylist_set_comment(qc_param, "ESO QC BIAS MASTER MEDIAN",
                                 "Median master bias level [ADU]");

    //Save computed readout noise
    double avg_ron = 0;
    for(size_t iport = 0; iport< bias_ccd_config.nports(); ++iport)
    {
        std::ostringstream key_stream;
        key_stream<<"ESO QC DET OUT"<<iport+1<<" RON";
        cpl_propertylist_append_double(qc_param, key_stream.str().c_str(),
                                       bias_ccd_config.computed_ron(iport));
        cpl_propertylist_set_comment(qc_param, key_stream.str().c_str(),
                                     "Readout noise per port from overscan [ADU]");
        avg_ron += bias_ccd_config.computed_ron(iport);
    }
    avg_ron /= bias_ccd_config.nports();

    /* QC.RON */
    cpl_propertylist_append_double(qc_param,  "ESO QC RON",
                          avg_ron);
    cpl_propertylist_set_comment(qc_param, "ESO QC RON", 
                                 "Average read out noise over all the ports [ADU]");


    cpl_image_delete(sample1_im);
    cpl_image_delete(sample2_im);
    cpl_image_delete(sample3_im);
    cpl_image_delete(sample4_im);    
}
