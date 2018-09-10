/* $Id: vmspflat.c,v 1.4 2011-03-14 15:28:58 cizzo Exp $
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
#include <cxmessages.h>

#include <cpl_recipe.h>
#include <cpl_plugininfo.h>
#include <cpl_parameterlist.h>
#include <cpl_frameset.h>

#include <pilmemory.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilutils.h>
#include <pilrecipe.h>
#include <pilqc.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmifutable.h"
#include "vmgrismtable.h"
#include "vmextractiontable.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmmosflat.h"
#include "vmmosmodels.h"
#include "vmmosextraction.h"
#include "vmmosutils.h"
#include "vmifuutils.h"
#include "vmqcutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


static cxint vmspflat(PilSetOfFrames *);


#define MAX_COMMENT_LENGTH  (80)
#define SATURATION_LEVEL (65535.)


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

typedef enum FlatMethod {
  FLAT_UNDEF,
  FLAT_AUTO,
  FLAT_MEDIAN,
  FLAT_AVERAGE,
  FLAT_POLY
} FlatMethod;

static const char *flatMethodNames[] = {
  "Auto",
  "Median",
  "Average",
  "Polynomial"
};

static const FlatMethod flatMethods[] = {
  FLAT_AUTO,
  FLAT_MEDIAN,
  FLAT_AVERAGE,
  FLAT_POLY
};

static unsigned int nFlatMethods = sizeof(flatMethods) / sizeof(FlatMethod);

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
vmspflat_create(cpl_plugin *plugin)
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
                                "Frames combination method.",
                                "vimos.Parameters",
                                "Average", 5, "Average", "Auto", "Ksigma", 
                                "MinMax", "Median");
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


    p = cpl_parameter_new_value("vimos.Parameters.extraction.fuzz",
                                CPL_TYPE_INT,
                                "Extra pixels from expected position of "
                                "spectrum edge in spectral extraction.",
                                "vimos.Parameters",
                                5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "Fuzz");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "Fuzz");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.flat.removing.trend",
                                CPL_TYPE_STRING,
                                "Flat field trends removal method.",
                                "vimos.Parameters",
                                "Median", 4, "Polynomial", "Median",
                                "Average", "Auto");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "FlatMethod");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "FlatMethod");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.removing.orderx",
                                CPL_TYPE_INT,
                                "Degree of polynomial used for fitting "
                                "the flat spectrum along X.",
                                "vimos.Parameters",
                                3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PolyOrderX");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PolyOrderX");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.removing.ordery",
                                CPL_TYPE_INT,
                                "Degree of polynomial used for fitting "
                                "the flat spectrum along Y.",
                                "vimos.Parameters",
                                8);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PolyOrderY");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PolyOrderY");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.removing.boxsize",
                                CPL_TYPE_INT,
                                "Size of running box used for smoothing "
                                "along the dispersion direction (currently "
                                "unused).",
                                "vimos.Parameters",
                                11);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SmoothBoxSize");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SmoothBoxSize");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.flat.frames.validate",
                                CPL_TYPE_BOOL,
                                "Consistency check on input frames.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ValidateFrames");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ValidateFrames");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.tolerance.low",
                                CPL_TYPE_DOUBLE,
                                "Minimum level for frame acceptance.",
                                "vimos.Parameters",
                                1000.);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LevelToleranceLow");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LevelToleranceLow");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.tolerance.high",
                                CPL_TYPE_DOUBLE,
                                "Maximum level for frame acceptance.",
                                "vimos.Parameters",
                                60000.);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LevelToleranceHigh");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LevelToleranceHigh");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.tolerance.pattern",
                                CPL_TYPE_DOUBLE,
                                "Threshold for flux pattern consistency.",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PatternTolerance");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PatternTolerance");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

    p = cpl_parameter_new_value("vimos.Parameters.badpixel.clean",
                                CPL_TYPE_BOOL,
                                "Bad pixel correction on master flat.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.clean",
                                CPL_TYPE_BOOL,
                                "Cosmic ray events cleaning from raw flats.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanCosmic");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanCosmic");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.quality.apply",
                                CPL_TYPE_BOOL,
                                "Quality control of the master flat.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyQC");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

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
vmspflat_exec(cpl_plugin *plugin)
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

    if (vmspflat(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmspflat");
        
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
vmspflat_destroy(cpl_plugin *plugin)
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
 *   Create master spectral flat field from a sequence of raw flat fields
 *   exposure taken with different shutter positions. Determine optical
 *   distorsion and spectral curvature models, and write them to the
 *   output flat field header.
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param sof  Pointer to a list of input frames containing
 *             the reference to flat field images, a master
 *             bias, and optional master dark, CCD Table,
 *             Grism Table, and (in case) an IFU Table.
 *
 * @doc
 *   The recipe function combines a set of spectroscopic flat fields
 *   into a normalized master spectral flat field, after bias removal
 *   and optional dark subtraction. A CCD table may also be specified,
 *   to be used in bad pixel correction on the output master spectral
 *   flat field. An IFU table must be specified in case of an IFU
 *   observation. All the input frames are passed via the set of frames
 *   \textbf{sof}. If requested, cosmic ray events are removed from
 *   each raw flat field before combination. A spectral extraction
 *   table, based on the information contained in the input grism
 *   table and in the input image header, is generated: if the input
 *   grism table carries also the coefficients of an optical distorsion
 *   model and of a spectral curvature model, they supersede those
 *   found in the image header. The input frames having the same
 *   shutter setting are stacked, and each spectrum from each stacked
 *   frame is extracted and inserted in the result master spectral flat
 *   field. The large scale trends are then removed from the result,
 *   which is finally normalized. On successful termination the created
 *   master flat field is added to the set of frames.
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
 *
 *   \begin{itemize}
 *
 *     \item FlatMethod:         Method used in modeling the 2D spectrum
 *                               of each slit, for determining the large
 *                               scale trends related to the CCD intrinsic
 *                               response variation with changing wavelength
 *                               and position. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Median:             Median filter
 *
 *       \item Average:            Average filter
 *
 *       \item Polynomial:         Polynomial fitting of spectrum
 *
 *     \end{itemize}
 *
 *     \item SmoothBoxSize:      Size of the running box used for smoothing,
 *                               expressed as number of pixels along the
 *                               dispersion direction (ignored if FlatMethod
 *                               is set to Polynomial).
 *
 *     \item PolyOrderX:         Degree of polynomial used for fitting
 *                               the flat spectrum in the X direction.
 *                               Used if FlatMethod is set to Polynomial.
 *
 *     \item PolyOrderY:         Degree of polynomial used for fitting
 *                               the flat spectrum in the Y direction.
 *                               Used if FlatMethod is set to Polynomial.
 *
 *     \item StackMethod:        Method for raw spectral flat field frames
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
 *     \item BiasMethod:         Method for bias removal from raw spectral
 *                               flat field frames. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Master:             Prescan and overscan regions are trimmed
 *                                 away from raw spectral flat fields after
 *                                 master bias subtraction.
 *
 *       \item Zmaster:            After master bias subtraction the residual
 *                                 signal found in each raw spectral flat
 *                                 field overscan regions is modelled and
 *                                 subtracted from the image. Next, prescan
 *                                 and overscan regions are trimmed away.
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
 *                               much different intensity pattern. Images
 *                               consistency is checked only against images
 *                               with the same shutter position.
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
 *     \item Fuzz:               Extra number of pixels to be used in
 *                               spectra extraction, when searching for
 *                               a spectrum edge on the CCD.
 *
 *     \item CleanBadPixel:      Bad pixel correction on output master
 *                               spectral flat field. If this option is
 *                               turned on, the recipe expects a CCD\_TABLE
 *                               in the input set of frames.
 *
 *     \item CleanCosmic:        Cosmic ray events removal from input raw
 *                               spectral flat field frames. If a CCD\_TABLE
 *                               is found in the input set of frames, bad
 *                               pixels will not be used in computing the
 *                               values to replace the cosmic rays events.
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
 * @author C. Izzo
 */

static cxint 
vmspflat(PilSetOfFrames *sof)
{

  const char  task[] = "vmspflat";

  const char  parameter[]        = "Parameters";

  const char *mosFlatCategory    = pilTrnGetCategory("MosScreenFlat");
  const char *mosMasterCategory  = pilTrnGetCategory("MosMasterScreenFlat");
  const char *mosCombCategory    = pilTrnGetCategory("MosCombScreenFlat");
  const char *masterBiasCategory = pilTrnGetCategory("MasterBias");
  const char *masterDarkCategory = pilTrnGetCategory("MasterDark");
  const char *ccdTableCategory   = pilTrnGetCategory("CcdTable");
  const char *grismTableCategory = pilTrnGetCategory("GrismTable");
  const char *ifuFlatCategory    = pilTrnGetCategory("IfuScreenFlat");
  const char *ifuMasterCategory  = pilTrnGetCategory("IfuMasterScreenFlat");
  const char *ifuTableCategory   = pilTrnGetCategory("IfuTable");
  const char *zeroOrderCategory  = pilTrnGetCategory("MosZeroOrder");
  const char *pafFileCategory    = pilTrnGetCategory("PAFCategory");

  char       *flatCategory       = NULL;
  char       *masterCategory     = NULL;

  char       *combMethodTag      = NULL;
  char       *biasMethodTag      = NULL;
  char       *flatMethodTag      = NULL;
  char       *pafFileName        = NULL;
  char        namePAF[]          = "vmCrvOpt";
  char        masterFlatName[PATHNAME_MAX + 1];
  char        zeroOrderName[PATHNAME_MAX + 1];
  char        comment[MAX_COMMENT_LENGTH];

  VimosBool   updateOK           = VM_TRUE;

  int         flatCount, mosFlatCount, ifuFlatCount;
  int         minFrames;
  int         goodFrames;

  int         i, j, k;
  int         cleanBadPixel, cleanCosmic, validateFrames;
  int         singleFrames, computeQC;
  int         error;

  int         biasMethodEntry;
  int         combMethodEntry;
  int         flatMethodEntry;

  BiasMethod     biasMethod = BIAS_UNDEF;

  CombMethod     combMethod = COMB_UNDEF;
  CombParameters combParameter;

  CombMethod     flatMethod = FLAT_UNDEF;
  int            smoothBoxSize;
  int            polyDegX, polyDegY;
  int            normMethod;
  FilterMethod   filterMethod = FILTER_UNDEF;


  int            fuzz;

  int            grismFlag = 0;

  int            shutter;
  int            shutPosCount;
  int           *flatCountShutPos;
  int            nCombFrames;

  float         *flatLevel = NULL;
  float         *flatRon = NULL;
  float         *flatGain = NULL;

  float         *extrData;
  float          distance, newDistance;
  float          shutPosMed;
  float          shutLow;
  float          shutHigh;
  float          lambda;
  
  int            startX, startY, windowDimen;
  int            quadrant;

  double         focalScale;
  double         flux, flux_err;
  double         time, gain;
  double         sumExposures = 0.0;
  double         lowerThreshold, upperThreshold;
  double         patternTolerance;

  int            crvOrdX, crvOrdY, fitCrv;
  int            numSlits;
  double         x, y, minSlitx, maxSlitx, minSlity, maxSlity;

  PilFrame      *biasFrame, *darkFrame, *ccdFrame, *grismFrame, *ifuFrame;
  PilFrame      *currFrame, *outputFrame;


  VimosImage   **flatList;
  VimosImage   **flatSubList;
  VimosImage   **combList;
  VimosImage    *masterFlat;
  VimosImage   **stackImages;
  VimosImage    *biasImage = NULL;
  VimosImage    *biasOver;
  VimosImage    *darkImage = NULL;
  VimosImage    *grismFile;
  VimosImage    *ifuFile = NULL;
  VimosTable    *ccdTable = NULL;
  VimosTable    *grismTable = NULL;
  VimosIfuTable *ifuTable = NULL;

  VimosExtractionTable *extractionTable = NULL;
  VimosExtractionTable *saveExtract;
  VimosDistModel2D     *optModX, *optModY;
  VimosDistModelFull   *crvMod;

  VimosExtractionSlit  *goodSlit;
  VimosExtractionSlit  *slit;

  int                   flag;


/*  Uncomment for debug purposes *
  VimosImage    *extraction_table_file;
 * */
 

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
   * Check if a single frame for each shutter position should be tolerated. 
   * In such a case, instead of signalling an error the stacking method will 
   * be ignored.
   */

  singleFrames = pilDfsDbGetBool("Parameters", "AllowSingleFrames", 0);

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
   * Determine the flat spectrum modeling method with its dependent
   * parameters.
   */

  flatMethodTag = (char *)pilDfsDbGetString(parameter, "FlatMethod");

  if ((flatMethodEntry = strselect(flatMethodTag,
                                   flatMethodNames, nFlatMethods)) < 0) {
    cpl_msg_error(task, "%s: Invalid flat modeling method.", flatMethodTag);
    return EXIT_FAILURE;
  }

  flatMethod = flatMethods[flatMethodEntry];

  if (flatMethod == FLAT_POLY) {
    normMethod = 1;
    polyDegX = pilDfsDbGetInt(parameter, "PolyOrderX", 3);
    polyDegY = pilDfsDbGetInt(parameter, "PolyOrderY", 8);
    if (polyDegX < 2 || polyDegY < 2) {
      cpl_msg_error(task, "Invalid choice of polynomial for spectral "
                  "flat modeling: degree should be at least 2");
      return EXIT_FAILURE;
    }
  }
  else {
    normMethod = 0;

    switch(flatMethod) {
    case FLAT_AUTO   :
    case FLAT_MEDIAN : filterMethod = FILTER_MEDIAN; break;
    case FLAT_AVERAGE: filterMethod = FILTER_AVERAGE; break;
    default : cpl_msg_error(task, "Invalid flat field normalization method");
              return EXIT_FAILURE;
    }

    smoothBoxSize = pilDfsDbGetInt(parameter, "SmoothBoxSize", 11);

    if (smoothBoxSize < 3) {
      cpl_msg_error(task,
                  "Invalid smooth box size, %d: it should be at least 3",
                  smoothBoxSize);
      return EXIT_FAILURE;
    }
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
                                       "LevelToleranceHigh", 60000.);

    if (patternTolerance < 1.0) {
      cpl_msg_error(task,
                  "Invalid tolerance, it should be at least one sigma.");
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

  cleanCosmic = pilDfsDbGetBool(parameter, "CleanCosmic", 0);


  /*
   * Check if the flat field quality control should be applied to the
   * final product (TO BE IMPLEMENTED).
   */

  //applyQc = pilDfsDbGetBool(parameter, "ApplyQC", 0);


  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   * Get number of raw data frames in the input set.
   */

  mosFlatCount = pilSofFrameCount(sof, mosFlatCategory);
  ifuFlatCount = pilSofFrameCount(sof, ifuFlatCategory);

  flatCount = mosFlatCount + ifuFlatCount;

  if (flatCount == 0) {
    cpl_msg_error(task, "No input raw spectral flat field frames found!");
    return EXIT_FAILURE;
  }

  if (mosFlatCount && ifuFlatCount) {
    cpl_msg_error(task,"Both MOS and IFU flat fields found in input");
    return EXIT_FAILURE;
  }

  if (mosFlatCount) {
    flatCount = mosFlatCount;
    flatCategory = (char *)mosFlatCategory;
    masterCategory = (char *)mosMasterCategory;
  }

  if (ifuFlatCount) {
    flatCount = ifuFlatCount;
    flatCategory = (char *)ifuFlatCategory;
    masterCategory = (char *)ifuMasterCategory;
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

  ccdFrame = pilSofLookup(sof, ccdTableCategory);
  if (ccdFrame)
    pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

  if (cleanBadPixel || cleanCosmic) {
    if (ccdFrame) {
      cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));
      if ((ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
        closeFitsTable(ccdTable, 0);
      }
      else {
        cpl_msg_error(task, "Failure opening CCD table");
        return EXIT_FAILURE;
      }
    }
    else {
      if (cleanBadPixel) {
        cpl_msg_error(task, "No CCD table in input: cannot clean bad pixels.");
        return EXIT_FAILURE;
      }
    }
  }


 /*
  * The grism table is required.
  */

  error = 1;

  if ((grismFrame = pilSofLookup(sof, grismTableCategory))) {
    pilFrmSetType(grismFrame, PIL_FRAME_TYPE_CALIB);
    if ((grismFile = openOldFitsFile(pilFrmGetName(grismFrame), 0, 0))) {
      if ((grismTable = newGrismTable())) {
        if (readFitsGrismTable(grismTable, grismFile->fptr) == VM_TRUE) {
          error = 0;
          closeFitsImage(grismFile, 0);
        }
        else {
          cpl_msg_error(task, "Failure reading grism table");
          deleteTable(grismTable);
        }
      }
      else
        cpl_msg_error(task, "Not enough memory");
    }
    else
      cpl_msg_error(task, "Failure opening grism table");
  }
  else
    cpl_msg_error(task,"No input Grism Table found");

  if (error) {
    deleteTable(ccdTable);
    return EXIT_FAILURE;
  }

 /*
  * The IFU table should be present for IFU data.
  */

  if (ifuFlatCount) {

    error = 1;

    if ((ifuFrame = pilSofLookup(sof, ifuTableCategory))) {
      if ((ifuFile = openOldFitsFile(pilFrmGetName(ifuFrame), 0, 0))) {
        if ((ifuTable = newIfuTable())) {
          if (readFitsIfuTable(ifuTable, ifuFile->fptr) == VM_TRUE) {
            error = 0;
            closeFitsImage(ifuFile, 0);
          }
          else {
            cpl_msg_error(task, "Failure reading IFU Table");
            deleteIfuTable(ifuTable);
          }
        }
        else
          cpl_msg_error(task, "Not enough memory");
      }
      else
        cpl_msg_error(task, "Failure opening IFU table");
    }
    else
      cpl_msg_error(task,"No input IFU Table found");

    if (error) {
      deleteTable(ccdTable);
      deleteTable(grismTable);
      return EXIT_FAILURE;
    }
  }


 /*
  * Get the master bias frame
  */

  error = 1;

  if ((biasFrame = pilSofLookup(sof, masterBiasCategory))) {
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
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    return EXIT_FAILURE;
  }


 /*
  * Get the (optional) master dark frame
  */

  if ((darkFrame = pilSofLookup(sof, masterDarkCategory))) {
    pilFrmSetType(darkFrame, PIL_FRAME_TYPE_CALIB);
    if ((darkImage = openOldFitsFile(pilFrmGetName(darkFrame), 1, 0))) {
      closeFitsImage(darkImage, 0);
    }
    else {
      cpl_msg_error(task, "Failure opening master dark frame");
      deleteImage(biasImage);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }
  }
  else
    cpl_msg_warning(task, "No master dark in input, dark subtraction "
                  "will not be performed");


 /*
  * Load just the headers of the raw spectral flat field frames.
  */

  error = 1;

  if ((flatList = (VimosImage **)cpl_calloc(flatCount, sizeof(VimosImage*)))) {

    error = 0;
    currFrame = pilSofLookupNext(sof, flatCategory);

    for (i = 0; i < flatCount; i++) {
      flatList[i] = openOldFitsFile(pilFrmGetName(currFrame), 0, 0);
      pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
      if ((loadFitsHeader(flatList[i])) == VM_FALSE) {
        error = 1;
        cpl_msg_error(task,
                    "Failure opening spectral flat field frame %d", i + 1);
        for (j = 0; j < i; j++) {
          closeFitsImage(flatList[j], 0);
          deleteImage(flatList[j]);
        }
        if (flatList[i])
          closeFitsImage(flatList[i], 0);
        cpl_free(flatList);
        break;
      }
      currFrame = pilSofLookupNext(sof, NULL);
    }
  }
  else
    cpl_msg_error(task,
                "Failure creating list of input raw spectral flat fields");

  if (error) {
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    return EXIT_FAILURE;
  }


 /*
  * Determine the quadrant from the first image on list
  */

  if (readIntDescriptor(flatList[0]->descs, pilTrnGetKeyword("Quadrant"),
                        &quadrant, comment) == VM_FALSE) {
    cpl_msg_error(task, "Cannot read descriptor %s",
                pilTrnGetKeyword("Quadrant"));
    for (i = 0; i < flatCount; i++) {
      closeFitsImage(flatList[i], 0);
      deleteImage(flatList[i]);
    }
    cpl_free(flatList);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    return EXIT_FAILURE;
  }

  if (readDoubleDescriptor(flatList[0]->descs, pilTrnGetKeyword("FocalScale"),
                        &focalScale, NULL) == VM_FALSE) {
    focalScale = 1.718;
    cpl_msg_warning(task, "Cannot read descriptor %s (defaulted to %f)",
                  pilTrnGetKeyword("FocalScale"), focalScale);
  }


  /* ALEX: if IFU, we want to first of all set to zero the first guess of
     curvature model coefficients in the header of the first image in the list,
     otherwise the extraction table will be very bad */

  if (ifuFlatCount) {
    if ((ifuDeleteCrvMod(flatList[0])) == VM_FALSE) {
      cpl_msg_error(task,"Failure setting to zero curvature model");
      return EXIT_FAILURE;
    }
  }

 /*
  * Create the extraction table from the header of any of the flat fields,
  * from the grism table, and if present from the IFU table.
  */

  if (!(extractionTable = VmSpExTab(flatList[0], grismTable, ifuTable,
                                    NULL))) {
    cpl_msg_error(task, "Cannot derive the extraction table");
    for (i = 0; i < flatCount; i++) {
      closeFitsImage(flatList[i], 0);
      deleteImage(flatList[i]);
    }
    cpl_free(flatList);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    return EXIT_FAILURE;
  }


 /*
  * Recreate bias overscans using as reference any of input flat fields.
  */

  if ((biasOver = growOverscans(biasImage, flatList[0]))) {
    if (biasImage != biasOver) {
      deleteImage(biasImage);
      biasImage = biasOver;
    }
  }
  else {
    cpl_msg_error(task, "Failure growing overscans in master bias");
    for (i = 0; i < flatCount; i++) {
      closeFitsImage(flatList[i], 0);
      deleteImage(flatList[i]);
    }
    cpl_free(flatList);
    deleteExtractionTable(extractionTable);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    return EXIT_FAILURE;
  }


 /*
  * The list of input images must now be sorted according to shutter
  * position. Resorting is done in place. The number of different
  * shutter positions found is returned in shutPosCount, and the
  * number of images corresponding to each shutter position is given
  * in the integer array flatCountShutPos[shutPosCount].
  */

  flatCountShutPos = sortByShutterPosition(flatList, flatCount, &shutPosCount);

 /*
  * Check if with so many frames per shutter position it is possible
  * to validate and stack frames using the specified method.
  */

  for (shutter = 0; shutter < shutPosCount; shutter++) {

    if (flatCountShutPos[shutter] < minFrames) {
      if (flatCountShutPos[shutter] == 1 && singleFrames) {
        combMethod = COMB_UNDEF;
      }
      else {
        cpl_msg_error(task, "Not enough raw spectral flat field frames in "
                    "input for stacking method '%s'!", combMethodTag);

        for (j = 0; j < flatCount; j++) {
          closeFitsImage(flatList[j], 0);
          deleteImage(flatList[j]);
        }
        cpl_free(flatList);
        cpl_free(flatCountShutPos);
        deleteExtractionTable(extractionTable);
        deleteImage(biasImage);
        deleteImage(darkImage);
        deleteTable(ccdTable);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        return EXIT_FAILURE;
      }
    }

    if (validateFrames && flatCountShutPos[shutter] < 2) {
      validateFrames = 0;
      cpl_msg_warning(task, "Too few spectral flat field frames (%d) in input. "
                    "Skipping frame selection task!", flatCount);
    }

  }


 /*
  * Prepare storage for stacked frames
  */

  combList = (VimosImage **)cpl_calloc(shutPosCount, sizeof(VimosImage *));

  if (!combList) {
    cpl_msg_error(task, "Not enough memory");
    for (i = 0; i < flatCount; i++) {
      closeFitsImage(flatList[i], 0);
      deleteImage(flatList[i]);
    }
    cpl_free(flatList);
    cpl_free(flatCountShutPos);
    deleteExtractionTable(extractionTable);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    return EXIT_FAILURE;
  }

 /*
  * Here begins the mega-loop on the different shutter positions.
  * Frames corresponding to different shutter positions are
  * processed separately, for a more constrained memory usage,
  * and then stacked.
  */

  flatSubList = flatList;

  for (shutter = 0; shutter < shutPosCount; shutter++) {

   /*
    * Just load into memory the images for a given shutter position.
    */

    for (i = 0;  i < flatCountShutPos[shutter]; i++) {
      if ((loadFitsData(flatSubList[i])) == VM_FALSE) {
        cpl_msg_error(task, "Failure loading spectral flat field frame %d/%d "
                    "of shutter position %d/%d", i + 1, 
                    flatCountShutPos[shutter], shutter + 1, shutPosCount);

        for (j = 0, k = 0; k < shutter; k++)
          j += flatCountShutPos[k];
        for (k = j; k < flatCount; k++) {
          closeFitsImage(flatList[k], 0);
          deleteImage(flatList[k]);
        }
        cpl_free(flatList);

        for (i = 0; i < shutPosCount; i++)
          deleteImage(combList[i]);
        cpl_free(combList);
        cpl_free(flatCountShutPos);
        deleteExtractionTable(extractionTable);
        deleteImage(biasImage);
        deleteImage(darkImage);
        deleteTable(ccdTable);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        return EXIT_FAILURE;
      }
    }


   /*
    * Get gain factors (e-/ADU)
    */

    error = 1;

    if ((flatGain = 
        (float *)cpl_calloc(flatCountShutPos[shutter], sizeof(float)))) {
      error = 0;
      for (i = 0; i < flatCountShutPos[shutter]; i++) {
        gain = getMeanGainFactor(flatSubList[i]);
        if (gain > MIN_DIVISOR) {
          flatGain[i] = gain;
        }
        else {
          error = 1;
          cpl_msg_error(task, "Wrong or no gain factor for flat field %d/%d "
                      "of shutter position %d/%d", i + 1, 
                      flatCountShutPos[shutter], shutter + 1, shutPosCount);
          cpl_free(flatGain);
          break;
        }
      }
    }
    else
      cpl_msg_error(task, "Not enough memory");
  
    if (error) {
      for (j = 0, k = 0; k < shutter; k++)
        j += flatCountShutPos[k];
      for (k = j; k < flatCount; k++) {
        closeFitsImage(flatList[k], 0);
        deleteImage(flatList[k]);
      }
      cpl_free(flatList);
      for (i = 0; i < shutPosCount; i++)
        deleteImage(combList[i]);
      cpl_free(combList);
      cpl_free(flatCountShutPos);
      deleteExtractionTable(extractionTable);
      deleteImage(biasImage);
      deleteImage(darkImage);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }


   /*
    * Noise contributions for each input spectral flat field frame must be
    * estimated before any processing (as bias subtraction, normalization,
    * etc.) is applied.
    */

    if ((flatRon = 
        (float *)cpl_calloc(flatCountShutPos[shutter], sizeof(float)))) {
      for (i = 0; i < flatCountShutPos[shutter]; i++) {
        flatRon[i] = computeAverageRon(flatSubList[i]);
      }
    }
    else {
      cpl_msg_error(task, "Not enough memory");

      for (j = 0, k = 0; k < shutter; k++)
        j += flatCountShutPos[k];
      for (k = j; k < flatCount; k++) {
        closeFitsImage(flatList[k], 0);
        deleteImage(flatList[k]);
      }
      cpl_free(flatList);
      for (i = 0; i < shutPosCount; i++)
        deleteImage(combList[i]);
      cpl_free(combList);
      cpl_free(flatCountShutPos);
      cpl_free(flatGain);
      deleteExtractionTable(extractionTable);
      deleteImage(biasImage);
      deleteImage(darkImage);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
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

    cpl_msg_info(task, "Bias removal for shutter position %d/%d frames...",
               shutter + 1, shutPosCount);

    for (i = 0; i < flatCountShutPos[shutter]; i++) {
      if (VmSubBias(flatSubList[i], biasImage, biasMethod) == EXIT_FAILURE) {
        cpl_msg_error(task, "Cannot remove bias from input flat field %d/%d "
                    "of shutter position %d/%d", i + 1,
                    flatCountShutPos[shutter], shutter + 1, shutPosCount);
        for (j = 0, k = 0; k < shutter; k++)
          j += flatCountShutPos[k];
        for (k = j; k < flatCount; k++) {
          closeFitsImage(flatList[k], 0);
          deleteImage(flatList[k]);
        }
        cpl_free(flatList);
        for (i = 0; i < shutPosCount; i++)
          deleteImage(combList[i]);
        cpl_free(combList);
        cpl_free(flatCountShutPos);
        cpl_free(flatGain);
        cpl_free(flatRon);
        deleteExtractionTable(extractionTable);
        deleteImage(biasImage);
        deleteImage(darkImage);
        deleteTable(ccdTable);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        return EXIT_FAILURE;
      }
    }


   /*
    * Now go for the (optional) dark subtraction
    */

    if (darkImage) {
      cpl_msg_info(task, 
                 "Dark subtraction for shutter position %d/%d frames...",
                 shutter + 1, shutPosCount);
  
      for (i = 0; i < flatCountShutPos[shutter]; i++) {
        if (VmSubDark(flatSubList[i], darkImage) == EXIT_FAILURE) {
          cpl_msg_error(task,
                      "Cannot subtract dark from input flat field %d/%d "
                      "of shutter position %d/%d", i + 1,
                      flatCountShutPos[shutter], shutter + 1, shutPosCount);
          for (j = 0, k = 0; k < shutter; k++)
            j += flatCountShutPos[k];
          for (k = j; k < flatCount; k++) {
            closeFitsImage(flatList[k], 0);
            deleteImage(flatList[k]);
          }
          cpl_free(flatList);
          for (i = 0; i < shutPosCount; i++)
            deleteImage(combList[i]);
          cpl_free(combList);
          cpl_free(flatCountShutPos);
          cpl_free(flatGain);
          cpl_free(flatRon);
          deleteExtractionTable(extractionTable);
          deleteImage(darkImage);
          deleteTable(ccdTable);
          deleteTable(grismTable);
          deleteIfuTable(ifuTable);
          return EXIT_FAILURE;
        }
      }
    }


   /*
    * Flat fields levels are computed on a window centered on the 
    * slit closest to the middle of the illuminated part of the
    * image.
    */

    if ((flatLevel = 
        (float *)cpl_calloc(flatCountShutPos[shutter], sizeof(float)))) {

     /* FIXME:
      * Here the "flag" is used to switch between determining the
      * illumination levels by averaging the whole image (i.e., 
      * including a lot of background), or determining the levels
      * averaging just a window centered on a spectrum. The problem
      * with the latter method is that the little window has a high
      * chance to miss its target. The reason of that is that the 
      * distorsion models first guesses, currently contained in the 
      * extraction table, are not so precise to position the little 
      * window accurately enough. As a result, instead of the flat 
      * fields levels one gets just the level of the background.
      *
      * UPGRADE: now the first guesses seem to be accurate enough
      * to allow centering well enough the central slit. It is true
      * that some background may get in, but not so much as it can
      * happen with flag = 0! The method with flag=2 has been added: 
      * similarly to what is done with the computation of the QC1
      * parameters, the slit closer to center is taken, and the
      * mean count level around the reference wavelength is computed.
      * (C.Izzo, September 2005).
      */

      flag = 2;
      if (flag == 1) {
        for (i = 0; i < flatCountShutPos[shutter]; i++) {
          flatLevel[i] = imageMean(flatSubList[i]);

          cpl_msg_info(task, "Level of flat field %d/%d of shutter position "
                     "%d/%d is %10.4f ADU", i + 1, flatCountShutPos[shutter], 
                     shutter + 1, shutPosCount, flatLevel[i]);
        }
      }
      else if (flag == 0) {

        readFloatDescriptor(flatSubList[0]->descs,
                            pilTrnGetKeyword("MshuPosL", quadrant),
                            &shutLow, comment);
  
        readFloatDescriptor(flatSubList[0]->descs,
                            pilTrnGetKeyword("MshuPosH", quadrant),
                            &shutHigh, comment);

        shutPosMed = (shutHigh - shutLow) / 2 + shutLow;

        slit = extractionTable->slits;
        distance = fabs(slit->maskY->data[0] - shutPosMed);
        slit = slit->next;

        while (slit) {
          newDistance = fabs(slit->maskY->data[0] - shutPosMed);
          if (distance > newDistance) {
            distance = newDistance;
            goodSlit = slit;
          }
          slit = slit->next;
        }

       /*
        * Determine extraction window dimension.
        */

        windowDimen = (int)(goodSlit->numRows * 0.8);

        startX = goodSlit->ccdX->data[(int)(goodSlit->numRows/2)] 
               - windowDimen / 2;

        startY = goodSlit->ccdY->data[0] - windowDimen / 2;

        for (i = 0; i < flatCountShutPos[shutter]; i++) {

         /*
          * Extract sub-image centered on the spectrum and compute mean
          */

          extrData = extractFloatImage(flatSubList[i]->data, 
                                       flatSubList[i]->xlen, 
                                       flatSubList[i]->ylen, 
                                       startX, startY, 
                                       windowDimen, windowDimen);

          if (!extrData) {
            for (j = 0, k = 0; k < shutter; k++)
              j += flatCountShutPos[k];
            for (k = j; k < flatCount; k++) {
              closeFitsImage(flatList[k], 0);
              deleteImage(flatList[k]);
            }
            cpl_free(flatList);
            for (i = 0; i < shutPosCount; i++)
              deleteImage(combList[i]);
            cpl_free(combList);
            cpl_free(flatCountShutPos);
            cpl_free(flatGain);
            cpl_free(flatRon);
            deleteExtractionTable(extractionTable);
            deleteImage(darkImage);
            deleteTable(ccdTable);
            deleteTable(grismTable);
            deleteIfuTable(ifuTable);
            return EXIT_FAILURE;
          }

          flatLevel[i] = (float)computeAverageFloat(extrData, 
                                                    windowDimen * windowDimen);

          cpl_free(extrData);

          cpl_msg_info(task, "Level of flat field %d/%d of shutter position "
                     "%d/%d is %10.4f ADU", i + 1, flatCountShutPos[shutter], 
                     shutter + 1, shutPosCount, flatLevel[i]);
        }
      }  /*** End of flag == 0 ***/
      else if (flag == 2) {
        slit = slitClosestToCenter(extractionTable);
        readFloatDescriptor(extractionTable->descs,
                            pilTrnGetKeyword("WlenCen"), &lambda, NULL);
        for (i = 0; i < flatCountShutPos[shutter]; i++) {
          flux = 0.0;
          extractSpecLevel(flatSubList[i], slit, lambda, 2, &flux);
          flatLevel[i] = flux;

          cpl_msg_info(task, "Level of flat field %d/%d of shutter position "
                     "%d/%d is %10.4f ADU", i + 1, flatCountShutPos[shutter], 
                     shutter + 1, shutPosCount, flatLevel[i]);
        }
      }
    }
    else {
      cpl_msg_error(task, "Not enough memory");
      for (j = 0, k = 0; k < shutter; k++)
        j += flatCountShutPos[k];
      for (k = j; k < flatCount; k++) {
        closeFitsImage(flatList[k], 0);
        deleteImage(flatList[k]);
      }
      cpl_free(flatList);
      for (i = 0; i < shutPosCount; i++)
        deleteImage(combList[i]);
      cpl_free(combList);
      cpl_free(flatCountShutPos);
      cpl_free(flatGain);
      cpl_free(flatRon);
      deleteExtractionTable(extractionTable);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }


   /*
    * When RON and flat fields levels are known, but before any image
    * rescaling is applied, cosmic ray events might be removed from all
    * input images.
    */
  
    if (cleanCosmic) {
      cpl_msg_info(task, "Cleaning cosmic ray events...");
      for (i = 0; i < flatCountShutPos[shutter]; i++) {
        if (VmCosmicClean(flatSubList[i], ccdTable, 0, flatLevel[i], 
                          flatGain[i], flatRon[i], -1., -1.) 
                          == EXIT_FAILURE) {
          cpl_msg_error(task, "Cannot clean cosmic rays from raw sky "
                      "flat field frame %d/%d of shutter position %d/%d", 
                      i + 1, flatCountShutPos[shutter], shutter + 1,
                      shutPosCount);
          for (j = 0, k = 0; k < shutter; k++)
            j += flatCountShutPos[k];
          for (k = j; k < flatCount; k++) {
            closeFitsImage(flatList[k], 0);
            deleteImage(flatList[k]);
          }
          cpl_free(flatList);
          for (i = 0; i < shutPosCount; i++)
            deleteImage(combList[i]);
          cpl_free(combList);
          cpl_free(flatCountShutPos);
          cpl_free(flatGain);
          cpl_free(flatRon);
          cpl_free(flatLevel);
          deleteExtractionTable(extractionTable);
          deleteTable(ccdTable);
          deleteTable(grismTable);
          deleteIfuTable(ifuTable);
          return EXIT_FAILURE;
        }
        cpl_msg_info(task, "  Flat field %d/%d of shutter position %d/%d done", 
                   i + 1, flatCountShutPos[shutter], shutter + 1,
                   shutPosCount);
      }
    }

    cpl_free(flatRon);
    cpl_free(flatGain);


    if (validateFrames) {
  
      if (lowerThreshold > MIN_DIVISOR || upperThreshold < SATURATION_LEVEL) {

        cpl_msg_info(task, "Selecting good flat field frames for shutter "
                   "position %d/%d...", shutter + 1, shutPosCount);
  
        cpl_msg_info(task, "Valid range for mean spectral flat field flux: "
                   "%4.0f to %5.0f ADU", lowerThreshold, upperThreshold);

       /*
        * Remove from the input list all flat field frames which do not have
        * flux within a reasonable range. The lowest and highest values are
        * user-defined.
        */

        cpl_msg_info(task, "Checking levels of frames with shutter position %d", 
                   shutter);

       /*
        * Good images are moved to the beginning of the sublist:
        */

        goodFrames = applyListSelection(flatSubList, flatLevel, 
                                        flatCountShutPos[shutter],
                                        lowerThreshold, upperThreshold, 1);

       /*
        * Free the data buffer of rejected frames to spare memory.
        */

/*
        for (i = goodFrames; i < flatCountShutPos[shutter]; i++) {
          cpl_free(flatSubList[i]->data);
          flatSubList[i]->data = NULL;
        }
*/

      }
      else
        goodFrames = flatCountShutPos[shutter];
    }
    else
      goodFrames = flatCountShutPos[shutter];


    cpl_free(flatLevel);

   /*
    * Stack the surviving flat field frames if requested and if possible,
    * and cleanup the raw flat field frames afterwards.
    */

    if (combMethod == COMB_UNDEF) {
      combList[shutter] = duplicateImage(flatSubList[0]);
    }
    else if (goodFrames >= minFrames) {
      cpl_msg_info(task, "Combining %d frames with shutter position %d "
                       "using method '%s'", goodFrames, shutter + 1,
                       combMethodNames[combMethodEntry]);
      combList[shutter] = frComb(flatSubList, goodFrames, combMethod, 
                                 &combParameter,0);
    }
    else {
      cpl_msg_error(task, "Not enough good frames for combine method "
                  "'%s'!", combMethodNames[combMethodEntry]);
      for (j = 0, k = 0; k < shutter; k++)
        j += flatCountShutPos[k];
      for (k = j; k < flatCount; k++) {
        closeFitsImage(flatList[k], 0);
        deleteImage(flatList[k]);
      }
      cpl_free(flatList);
      for (i = 0; i < shutPosCount; i++)
        deleteImage(combList[i]);
      cpl_free(combList);
      cpl_free(flatCountShutPos);
      deleteExtractionTable(extractionTable);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }

    if (!combList[shutter]) {
      cpl_msg_error(task, "Stacking of raw spectral flat field frames "
                  "having shutter position %d, failed!", shutter);
      for (j = 0, k = 0; k < shutter; k++)
        j += flatCountShutPos[k];
      for (k = j; k < flatCount; k++) {
        closeFitsImage(flatList[k], 0);
        deleteImage(flatList[k]);
      }
      cpl_free(combList);
      for (i = 0; i < flatCount; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);
      cpl_free(flatCountShutPos);
      deleteExtractionTable(extractionTable);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }


   /*
    * To the combined frame copy always the header of the first image
    * with that shutter position.
    */

    copyAllDescriptors(flatSubList[0]->descs, &(combList[shutter]->descs));

    if (shutter == 0) {

     /*
      * Get exposure times
      */

     /* FIXME:   (but how?)
      * The total exposure time is here computed as the total exposure 
      * time of all selected frames having shutter position 0. It is 
      * clear that this is wrong in case all shutter positions are not 
      * 1) equally exposed, and 2) equally rejected by the validation 
      * algorithm!
      */

      error = 0;
      for (i = 0; i < goodFrames; i++) {
        if (readDoubleDescriptor(flatSubList[i]->descs,
                                 pilTrnGetKeyword("ExposureTime"),
                                 &time, comment) == VM_TRUE) {
          if (time < MIN_DIVISOR) {
            error = 1;
            cpl_msg_error(task, "Zero or negative exposure time in flat field");
            break;
          }
          else {
            sumExposures += time;
          }
        }
        else {
          error = 1;
          cpl_msg_error(task, "Cannot read %s from flat field frame.",
                      pilTrnGetKeyword("ExposureTime"));
          break;
        }
      }
    
      if (error) {
        for (i = 0; i < shutPosCount; i++)
          deleteImage(combList[i]);
        cpl_free(combList);
        for (j = 0, k = 0; k < shutter; k++)
          j += flatCountShutPos[k];
        for (k = j; k < flatCount; k++) {
          closeFitsImage(flatList[k], 0);
          deleteImage(flatList[k]);
        }
        cpl_free(flatList);
        cpl_free(flatCountShutPos);
        deleteExtractionTable(extractionTable);
        deleteTable(ccdTable);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        return EXIT_FAILURE;
      }

      nCombFrames = goodFrames;

    }


   /*
    * Free all images with the current shutter position.
    */

    for (i = 0; i < flatCountShutPos[shutter]; i++) {
      closeFitsImage(flatSubList[i], 0);
      deleteImage(flatSubList[i]);
    }

    flatSubList += flatCountShutPos[shutter];

  }  /* End of loop on shutter positions */


  deleteImage(biasImage);
  deleteImage(darkImage);
  cpl_free(flatList);
  cpl_free(flatCountShutPos);

  if (!(stackImages = VmSpStackFF(combList,
                                  shutPosCount, extractionTable, fuzz))) {
    cpl_msg_error(task, "Failure stacking flat fields");
    for (i = 0; i < shutPosCount; i++)
      deleteImage(combList[i]);
    cpl_free(combList);
    deleteExtractionTable(extractionTable);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    return EXIT_FAILURE;
  }

  for (i = 0; i < shutPosCount; i++)
    deleteImage(combList[i]);
  cpl_free(combList);


 /*
  * Rebuilding the extraction table, since it was eaten up by VmSpStackFF()
  */

  deleteExtractionTable(extractionTable);
  if (!(extractionTable = VmSpExTab(stackImages[0], grismTable, ifuTable,  
                                    NULL))) {
    cpl_msg_error(task, "Failure rebuilding the extraction table");
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    deleteImage(stackImages[0]);
    deleteImage(stackImages[1]);
    cpl_free(stackImages);
    return EXIT_FAILURE;
  }

  if (cleanBadPixel) {

   /*
    * Bad pixel cleaning
    */

    error = 0;

    cpl_msg_info(task, "Cleaning bad pixels from result frames...");

    if (cleanBadPixels(stackImages[0], ccdTable, 0) == EXIT_FAILURE) {
      cpl_msg_error(task,
                  "Failure cleaning bad pixels in stacked flat field frame");
      error = 1;
    }

    if (stackImages[1]) {
      if (cleanBadPixels(stackImages[1], ccdTable, 0) == EXIT_FAILURE) {
        cpl_msg_error(task, "Failure cleaning bad pixels in stacked zero "
                    "order flat field frame");
        error = 1;
      }
    }

    if (error) {
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      deleteExtractionTable(extractionTable);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }
  }

  deleteTable(ccdTable);


  /* if IFU mode, now re-shift extraction table to take into account 
     flexures before going to VmSpDerCurves & c.*/

  if (ifuFlatCount) {
    if (EXIT_FAILURE == ifuExtrShift(stackImages,extractionTable)) {
      cpl_msg_error(task, "Failure upgrading spectral curvature "
                  "polynomials on stacked flat field");
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      deleteExtractionTable(extractionTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }
  }

  /* uncomment to output extraction table after offsetting */
  /*
  extraction_table_file = newImage(0, 0, NULL);
  openNewFitsImage("extraction_off.fits", extraction_table_file);
  writeFitsExtractionTable(extractionTable, extraction_table_file->fptr);
  closeFitsImage(extraction_table_file, 0);
  */


  /*
   * Determine whether it is possible to fit a new curvature and 
   * optical distortion model
   */

  slit = extractionTable->slits;

  minSlitx = maxSlitx = slit->maskX->data[0];
  minSlity = maxSlity = slit->maskY->data[0];

  numSlits = 0;
  while (slit) {
    x = slit->maskX->data[0];
    y = slit->maskY->data[0];
    if (x > maxSlitx)
      maxSlitx = x;
    if (x < minSlitx)
      minSlitx = x;
    if (y > maxSlity)
      maxSlity = y;
    if (y < minSlity)
      minSlity = y;
    numSlits++;
    slit = slit->next;
  }

  readIntDescriptor(extractionTable->descs,
                    pilTrnGetKeyword("CurvatureOrdX"), &crvOrdX, NULL);

  readIntDescriptor(extractionTable->descs,
                    pilTrnGetKeyword("CurvatureOrdY"), &crvOrdY, NULL);

  /*
   * Conditions required for trying to improve the default distortion
   * models: slits should be at least two times the number of free
   * parameters, and their coordinates should span for at least 50 mm
   * on the mask.
   */

  fitCrv = 0;

  if (maxSlitx - minSlitx > 50.0 && maxSlity - minSlity > 50.0)
    if (numSlits > (crvOrdX + 1) * (crvOrdY + 1))
      fitCrv = 1;

  if (fitCrv) {

    /*
     * Light median filter to eliminate hot columns: they may disturb
     * determination of spectral curvature
     */

    /***

    filtIma = VmFrMedFil(stackImages[0], 3, 3, 0);

    if (filtIma) {
      copyAllDescriptors(stackImages[0]->descs, &(filtIma->descs));
      deleteImage(stackImages[0]);
      stackImages[0] = filtIma;
    }

    if (stackImages[1]) {
      filtIma = VmFrMedFil(stackImages[1], 3, 3, 0);
      if (filtIma) {
        copyAllDescriptors(stackImages[1]->descs, &(filtIma->descs));
        deleteImage(stackImages[1]);
        stackImages[1] = filtIma;
      }
    }

    ***/

    if (EXIT_FAILURE == VmSpDerCurves(stackImages, extractionTable, fuzz)) {
      cpl_msg_error(task, "Failure upgrading spectral curvature "
                  "polynomials on stacked flat field");
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      deleteExtractionTable(extractionTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }

    if (EXIT_FAILURE == VmSpCurveModel(extractionTable, 
                                       grismTable, grismFlag)) {
      cpl_msg_error(task, "Failure upgrading spectral curvature model "
                  "for stacked flat field");
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      deleteExtractionTable(extractionTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }

    if (EXIT_FAILURE == VmSpOptModel(extractionTable, grismTable, grismFlag)) {
      cpl_msg_error(task, "Failure upgrading optical distorsion model "
                  "for stacked flat field");
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      deleteExtractionTable(extractionTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }


    /* uncomment to output extraction table after writing models in header */
    /*
    extraction_table_file = newImage(0, 0, NULL);
    openNewFitsImage("extraction_optMod.fits", extraction_table_file);
    writeFitsExtractionTable(extractionTable, extraction_table_file->fptr);
    closeFitsImage(extraction_table_file, 0);
    */
  

   /*
    * Upgrading the extraction table, according to the new found models.
    */

    saveExtract = extractionTable;

/* It was unrecoverable:

    if (!(extractionTable = VmSpExTab(stackImages[0], grismTable, ifuTable,
                                      extractionTable))) {
      deleteExtractionTable(saveExtract);
      cpl_msg_error(task, "Failure upgrading the extraction table");
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      return EXIT_FAILURE;
    }

... and now it's like that: */

    if (!(extractionTable = VmSpExTab(stackImages[0], grismTable, ifuTable,
                                      extractionTable))) {
      deleteExtractionTable(saveExtract);
      cpl_msg_warning(task, "Failure computing new distortion models: "
                    "using default models.");
      if (!(extractionTable = VmSpExTab(stackImages[0], grismTable, 
                                        ifuTable, NULL))) {
        cpl_msg_error(task, "Failure upgrading the extraction table");
        deleteImage(stackImages[0]);
        deleteImage(stackImages[1]);
        cpl_free(stackImages);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        return EXIT_FAILURE;
      }
    }


  }
  else
    cpl_msg_warning(task, "Not enough data for distortion models computation: "
                  "using default models.");

  deleteTable(grismTable);
  deleteIfuTable(ifuTable);

  /* BG changed do not normalize in IFU case */
  if (ifuFlatCount) {

    masterFlat=duplicateImage(stackImages[0]);
    copyAllDescriptors(stackImages[0]->descs, &(masterFlat->descs));

  }
  else {
    if (!(masterFlat = VmSpNormFF(stackImages[0], extractionTable, normMethod,
				  polyDegX, polyDegY, filterMethod, 
                                  smoothBoxSize, smoothBoxSize))) {
      cpl_msg_error(task, "Failure normalizing stacked flat field");
      deleteExtractionTable(extractionTable);
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      return EXIT_FAILURE;
    }
  }


  if (computeQC) {

    cpl_msg_info(task, "Computing QC1 parameters...");

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      /*
       * Write PRO.CATG, ARCFILE, TPL ID, and OCS CON QUAD, to QC1 group.
       */

      pilQcWriteString("PRO.CATG", mosMasterCategory, "Product category");

      qcCopyValue(stackImages[0]->descs, 
                  pilTrnGetKeyword("ArchiveFile"), NULL, "Archive File Name");

      qcCopyValue(stackImages[0]->descs, 
                  pilTrnGetKeyword("TplId"), NULL, "Template signature ID");

      qcCopyValue(stackImages[0]->descs, 
                  pilTrnGetKeyword("Quadrant"), NULL, "Quadrant");

      qcCopyValue(stackImages[0]->descs, 
                  pilTrnGetKeyword("FilterName", quadrant), NULL, "Filter");

      qcCopyValue(stackImages[0]->descs, 
                  pilTrnGetKeyword("GrismName", quadrant), NULL, "Grism");

      qcCopyValue(stackImages[0]->descs, 
                  pilTrnGetKeyword("MaskId", quadrant), NULL, "Mask");

      qcCopyValue(extractionTable->descs, pilTrnGetKeyword("WlenCen"),
                  "Angstrom", "Reference wavelength");


      /*
       * Flux at reference wavelength in slit closest to center.
       */

      slit = slitClosestToCenter(extractionTable);

      qcWriteValueDouble(extractionTable->descs, slit->width * focalScale,
                         "QC.MOS.SLIT.WIDTH", "arcsec",
                         "Width of slit closest to center");

      readFloatDescriptor(extractionTable->descs, 
                          pilTrnGetKeyword("WlenCen"), &lambda, NULL);

      extractSpecFlux(stackImages[0], slit, lambda, 2, &flux, &flux_err);

      flux_err /= time;     /* Assuming all frames have the same exposure */
      flux     /= time;

      cpl_msg_info(task, "Flux at wavelength %.2f: %.2f +/- %.2f ADU/mm^2/s\n",
                 lambda, flux, flux_err);

      qcWriteValueDouble(masterFlat->descs, flux, "QC.MOS.FLAT.FLUX",
                         "ADU/mm^2/s", "Flux at reference wavelength");

      qcWriteValueDouble(masterFlat->descs, flux_err, "QC.MOS.FLAT.FLUXERR",
                         "ADU/mm^2/s", 
                         "Error on flux at reference wavelength");

      if (pilQcGroupEnd() == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

    }
    else
      cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

  } /* End of QC1 computation. */


  /*
   * Copy the curvature and the optical distorsion models to the master
   * flat field header.
   */

  if (readOptDistModel(extractionTable->descs, &optModX, &optModY) 
                                                             == VM_FALSE) {
    cpl_msg_error(task, "Failure copying the optical distorsion model "
                "from Extraction Table");
    deleteExtractionTable(extractionTable);
    deleteImage(masterFlat);
    deleteImage(stackImages[0]);
    deleteImage(stackImages[1]);
    cpl_free(stackImages);
    return EXIT_FAILURE;
  }


 /* FIXME:
  * Replace the next call with a call to writeOptDistModel() as soon as
  * distorsion models coefficients in image headers are no longer character
  * strings!!!
  */

  if (writeOptDistModelString(&(masterFlat->descs), optModX, optModY) 
                                                             == VM_FALSE) {
    cpl_msg_error(task, "Failure copying the optical distorsion model "
                "to Master Flat header");
    deleteExtractionTable(extractionTable);
    deleteDistModel2D(optModX);
    deleteDistModel2D(optModY);
    deleteImage(masterFlat);
    deleteImage(stackImages[0]);
    deleteImage(stackImages[1]);
    cpl_free(stackImages);
    return EXIT_FAILURE;
  }

  deleteDistModel2D(optModX);
  deleteDistModel2D(optModY);

  if (readCurvatureModel(extractionTable->descs, &crvMod) == VM_FALSE) {
    cpl_msg_error(task, "Failure copying the curvature model "
                "from Extraction Table");
    deleteExtractionTable(extractionTable);
    deleteImage(masterFlat);
    deleteImage(stackImages[0]);
    deleteImage(stackImages[1]);
    cpl_free(stackImages);
    return EXIT_FAILURE;
  }


 /* FIXME:
  * Replace the next call with a call to writeCurvatureModel() as soon as
  * distorsion models coefficients in image headers are no longer character
  * strings!!!
  */

  if (writeCurvatureModelString(&(masterFlat->descs), crvMod) == VM_FALSE) {
    cpl_msg_error(task, "Failure copying the curvature model "
                "to Master Flat header");
    deleteExtractionTable(extractionTable);
    deleteDistModelFull(crvMod);
    deleteImage(masterFlat);
    deleteImage(stackImages[0]);
    deleteImage(stackImages[1]);
    cpl_free(stackImages);
    return EXIT_FAILURE;
  }

  deleteDistModelFull(crvMod);

 /*
  *  Uncomment this to output the final extraction table to disk.
  */
  /*
  extraction_table_file = newImage(0, 0, NULL);
  openNewFitsImage("extraction.fits", extraction_table_file);
  writeFitsExtractionTable(extractionTable, extraction_table_file->fptr);
  closeFitsImage(extraction_table_file, 0);
  */
  /*
  */

  deleteExtractionTable(extractionTable);


 /*
  * Update the products header
  *
  * Note that for the moment also keywords which are not task specific
  * are handled here, since this is the last possibility to access
  * the linked list of keywords without reopening the file.
  * This may change in future!
  */

 /*
  * Normalized master flat field first:
  */

  updateOK = updateOK && insertDoubleDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(masterFlat),
                         "Minimum pixel value", "ESO*", 1);

  updateOK = updateOK && insertDoubleDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(masterFlat),
                         "Maximum pixel value", "ESO*", 1);

  updateOK = updateOK && writeDoubleDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("DataMean"),
                         imageMean(masterFlat),
                         "Mean pixel value");

  updateOK = updateOK && writeDoubleDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(masterFlat),
                         "Standard deviation of pixel values");

  updateOK = updateOK && writeDoubleDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(masterFlat),
                         "Median pixel value");

  updateOK = updateOK && writeIntDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("NFramesCombined"),
                         nCombFrames, "Number of frames combined");

  updateOK = updateOK && writeDoubleDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("ExposureTime"),
                         sumExposures, "Exposure time");

  updateOK = updateOK && writeDoubleDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("SummedExposureTime"),
                         sumExposures,
                         "Total exposure time of combined frames");

  updateOK = updateOK && insertStringDescriptor(&(masterFlat->descs),
                         pilTrnGetKeyword("DoCategory"), masterCategory,
                         "Category of pipeline product", "ESO*", 1);

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    deleteImage(masterFlat);
    deleteImage(stackImages[0]);
    deleteImage(stackImages[1]);
    cpl_free(stackImages);
    return EXIT_FAILURE;
  }

  deleteSetOfDescriptors(&(masterFlat->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(masterFlat->descs), "ESO DPR*");

  removeDescriptor(&(masterFlat->descs),
                                     pilTrnGetKeyword("TplExposureNumber"));

  if (ifuFlatCount == 0) {

   /*
    * NOT-normalized master flat field first:
    */

    updateOK = updateOK && insertDoubleDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("DataMin"),
                           imageMinimum(stackImages[0]),
                           "Minimum pixel value", "ESO*", 1);
  
    updateOK = updateOK && insertDoubleDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("DataMax"),
                           imageMaximum(stackImages[0]),
                           "Maximum pixel value", "ESO*", 1);
  
    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("DataMean"),
                           imageMean(stackImages[0]),
                           "Mean pixel value");
  
    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("DataStdDeviation"),
                           imageSigma(stackImages[0]),
                           "Standard deviation of pixel values");
  
    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("DataMedian"),
                           imageMedian(stackImages[0]),
                           "Median pixel value");

    updateOK = updateOK && writeIntDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("NFramesCombined"),
                           nCombFrames, "Number of frames combined");
  
    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("ExposureTime"),
                           sumExposures, "Exposure time");
  
    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("SummedExposureTime"),
                           sumExposures,
                           "Total exposure time of combined frames");
  
    updateOK = updateOK && insertStringDescriptor(&(stackImages[0]->descs),
                           pilTrnGetKeyword("DoCategory"), mosCombCategory,
                           "Category of pipeline product", "ESO*", 1);
  
    if (!updateOK) {
      cpl_msg_error(task, "Failure updating product header");
      deleteImage(stackImages[0]);
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      return EXIT_FAILURE;
    }

    deleteSetOfDescriptors(&(stackImages[0]->descs), "ESO ADA*");
    deleteSetOfDescriptors(&(stackImages[0]->descs), "ESO DPR*");

    removeDescriptor(&(stackImages[0]->descs),
                                       pilTrnGetKeyword("TplExposureNumber"));
  }

  if (stackImages[1]) {

   /*
    * Zero order frame:
    */

    updateOK = updateOK && insertDoubleDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("DataMin"),
                           imageMinimum(stackImages[1]),
                           "Minimum pixel value", "ESO*", 1);

    updateOK = updateOK && insertDoubleDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("DataMax"),
                           imageMaximum(stackImages[1]),
                           "Maximum pixel value", "ESO*", 1);

    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("DataMean"),
                           imageMean(stackImages[1]),
                           "Mean pixel value");

    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("DataStdDeviation"),
                           imageSigma(stackImages[1]),
                           "Standard deviation of pixel values");

    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("DataMedian"),
                           imageMedian(stackImages[1]),
                           "Median pixel value");

    updateOK = updateOK && writeIntDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("NFramesCombined"),
                           nCombFrames, "Number of frames combined");

    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("ExposureTime"), sumExposures,
                           "Exposure time");

    updateOK = updateOK && writeDoubleDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("SummedExposureTime"),
                           sumExposures,
                           "Total exposure time of combined frames");

    updateOK = updateOK && insertStringDescriptor(&(stackImages[1]->descs),
                           pilTrnGetKeyword("DoCategory"), zeroOrderCategory,
                           "Category of pipeline product", "ESO*", 1);

    if (!updateOK) {
      cpl_msg_error(task, "Failure updating product header");
      deleteImage(masterFlat);
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      return EXIT_FAILURE;
    }

    deleteSetOfDescriptors(&(stackImages[1]->descs), "ESO ADA*");
    deleteSetOfDescriptors(&(stackImages[1]->descs), "ESO DPR*");

    removeDescriptor(&(stackImages[1]->descs),
                                     pilTrnGetKeyword("TplExposureNumber"));

  }

 /*
  * Create the product files on disk, set the product attributes and
  * update the set of frames.
  */

 /*
  * Output the primary product: the spectral master flat field.
  */

  vmstrlower(strcpy(masterFlatName, masterCategory));
  strcat(masterFlatName, ".fits");

  if (createFitsImage(masterFlatName, masterFlat, masterCategory)) {
    outputFrame = newPilFrame(masterFlatName, masterCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", masterFlatName);
    deleteImage(masterFlat);
    deleteImage(stackImages[0]);
    deleteImage(stackImages[1]);
    cpl_free(stackImages);
    return EXIT_FAILURE;
  }


  if (ifuFlatCount == 0) {

   /*
    * Output the not-normalized spectral master flat field.
    * This is new...
    */

    vmstrlower(strcpy(masterFlatName, mosCombCategory));
    strcat(masterFlatName, ".fits");

    if (createFitsImage(masterFlatName, stackImages[0], mosCombCategory)) {
      outputFrame = newPilFrame(masterFlatName, mosCombCategory);
  
      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);
  
      pilSofInsert(sof, outputFrame);
    }
    else {
      cpl_msg_error(task, "Cannot create local product file %s!",
                  masterFlatName);
      deleteImage(masterFlat);
      deleteImage(stackImages[0]);
      deleteImage(stackImages[1]);
      cpl_free(stackImages);
      return EXIT_FAILURE;
    }
  }


  if (stackImages[1]) {

   /*
    * This is the zero order image
    */

    vmstrlower(strcpy(zeroOrderName, zeroOrderCategory));
    strcat(zeroOrderName, ".fits");

    if (createFitsImage(zeroOrderName, stackImages[1], zeroOrderCategory)) {
      outputFrame = newPilFrame(zeroOrderName, zeroOrderCategory);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
    }
  }
  else
    cpl_msg_info(task, "No zero order image created.");

  deleteImage(stackImages[0]);
  deleteImage(stackImages[1]);
  cpl_free(stackImages);

  if ((pafFileName = createSpectralDistModelsPAF(masterFlat->descs, namePAF))
                                                         != NULL) {
    outputFrame = newPilFrame(pafFileName, pafFileCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_PAF);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
    cpl_free(pafFileName);
  }
  else
    cpl_msg_error(task, "Cannot create local product file %s!", namePAF);

  deleteImage(masterFlat);

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
                    "vmspflat",
    "Build master flat field from a set of raw MOS flat fields.",
    "This recipe is used to create the MOS normalised master flat field\n" 
    "from a set of MOS flat field exposures. A _not_ normalised master\n" 
    "flat field is also required by the recipe vmspcaldisp to determine\n" 
    "the spectral curvature model. The recipe vmspflat does compute its\n" 
    "own curvature model, but this is generally incompatible with the\n" 
    "wavelength calibration and the optical distortion models computed\n"
    "by vmspcaldisp. For this reason this curvature model is just used\n" 
    "internally in the process of flat field normalisation, where the\n"
    "compatibility with the Y component of the optical distortion model\n" 
    "is irrelevant.\n\n"
    "Input files:\n\n"
    "  DO category:               Type:       Explanation:         Required:\n"
    "  MOS_SCREEN_FLAT            Raw         Flat field exposures    Y\n"
    "  MASTER_BIAS                Calib       Master bias             Y\n"
    "  MASTER_DARK                Calib       Master dark             .\n"
    "  GRISM_TABLE                Calib       Grism table             Y\n"
    "  CCD_TABLE                  Calib       Bad pixel table         .\n\n"
    "Output files:\n\n"
    "  DO category:               Data type:  Explanation:\n"
    "  MOS_MASTER_SCREEN_FLAT     FITS image  Master flat field\n"
    "  MOS_COMBINED_SCREEN_FLAT   FITS image  Combined flat field\n"
    "  (none)                     PAF         Distortion models\n\n"
    "At least one raw flat field exposure should be present in the input\n" 
    "SOF. The acquisition of input flat fields may be done using different\n" 
    "mask shutter settings (to avoid contamination between different\n" 
    "spectral orders in LR grisms used with multiplexed masks. The vmspflat\n" 
    "recipe will properly combine all the input frames according to a\n" 
    "specified method.\n\n"
    "The bad pixel table needs to be specified only if the cleaning of bad\n"
    "pixels is requested.\n\n"
    "The grism table contains necessary information to control the way\n" 
    "spectra are extracted, and the determination of the spectral distortion\n"
    "models. The vmspflat recipe gets from the grism table the wavelength\n" 
    "that should be used as reference (header entry PRO WLEN CEN), and the\n" 
    "spectrum extension in CCD pixels above and below the position of the\n" 
    "reference wavelength (header entries PRO SPECT LLEN LO and PRO SPECT\n" 
    "LLEN HI).\n\n"
    "The PAF file just contains the curvature and the optical distortion\n" 
    "models computed in the process of spectra normalisation. Such models\n" 
    "are available just for debug purposes, and are not to be used in\n" 
    "further data reduction steps.",

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

                    vmspflat_create,
                    vmspflat_exec,
                    vmspflat_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
