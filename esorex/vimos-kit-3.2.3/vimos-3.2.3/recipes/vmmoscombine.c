/* $Id: vmmoscombine.c,v 1.3 2011-03-14 15:28:58 cizzo Exp $
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
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <fitsio.h>

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


#define MAX_COMMENT_LENGTH (80)

/*
 * Definition of the label strings for all methods the recipe function
 * supports for combining frames.
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

static cxint vmmoscombine(PilSetOfFrames *);


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmmoscombine_create(cpl_plugin *plugin)
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


    p = cpl_parameter_new_value("vimos.Parameters.flux.calibration",
                                CPL_TYPE_BOOL,
                                "Extracted spectra are flux calibrated.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CalibrateFlux");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CalibrateFlux");
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
vmmoscombine_exec(cpl_plugin *plugin)
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

    if (vmmoscombine(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmmoscombine");
        
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
vmmoscombine_destroy(cpl_plugin *plugin)
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
 *   Combine reduced and 2D-extracted spectra from different OBs.
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param sof  Set of frames with references to at least two frames
 *             in 2D-stacked format (MOS_SCIENCE_EXTRACTED) created by
 *             either vmmosobsstare or vmmosobsjitter, and an atmospheric 
 *             extinction table. A spectro-photometric table may be
 *             added if a relative flux correction is requested.
 *
 * @doc 
 *   Reduce a sequence of at least two frames in 2D-stacked format 
 *   created by either vmmosobsstare or vmmosobsjitter from data 
 *   obtained with the same mask. It cannot be assumed that the
 *   2D-extracted slit spectra have correspondingly the same sizes
 *   in all images, as they might be produced by different runs
 *   of vmmosobsjitter on data sequences with different jittering
 *   span. For this reasons all images are resized to the same Y size
 *   of the largest input image, inserting the appropriate number
 *   of rows (padded with zeroes) at each slit spectrum. The object
 *   detection and alignment algorithms are exactly the same as
 *   those applied by the recipe vmmosobsjitter. Only the sky 
 *   subtraction is not repeated. The spectra are corrected to
 *   airmass zero before being combined.
 *
 *   \begin{itemize}
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
 *   \end{itemize}
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 *
 * @author C.Izzo
 */   

static cxint
vmmoscombine(PilSetOfFrames *sof)
{
  const char  task[]               = "vmmoscombine";

  const char  parameter[]          = "Parameters";

  const char *mosExtractedCategory = pilTrnGetCategory("MosScienceExtracted");
  const char *mosReducedCategory   = pilTrnGetCategory("MosScienceReduced");
  const char *mosFluxReducedCategory
                                   = pilTrnGetCategory("MosScienceFluxReduced");
  const char *sphotTableCategory   = pilTrnGetCategory("MosSphotTable");
  const char *grismTableCategory   = pilTrnGetCategory("GrismTable");
  const char *objectTableCategory  = pilTrnGetCategory("ObjectTable");
  const char *windowTableCategory  = pilTrnGetCategory("WindowTable");
  const char *atmTableCategory     = pilTrnGetCategory("ExtinctTable");


  char                  *combMethodTag = NULL;
  int                    combMethodEntry;

  CombMethod             combMethod = COMB_MEDIAN;
  CombParameters         combParameter;

  char                   output2DName[PATHNAME_MAX + 1];
  char                   output1DName[PATHNAME_MAX + 1];
  char                   objectTableName[PATHNAME_MAX + 1];
  char                   windowTableName[PATHNAME_MAX + 1];


  VimosBool              updateOK = VM_TRUE;

  unsigned int           calibrateFlux;
  unsigned int           error;

  int                    i, j;

  int                    slitMargin          = 2;
  int                    numLevels           = 0;
  int                    minObjectSize       = 0;
  int                    minCompositeSize    = 0;
  int                    mosCount            = 0;
  int                    winCount            = 0;
  int                    minFrames           = 0;
  int                    quadrant;
  int                    firstMaskId, maskId;

  float                  detLevel            = 0.0;
  float                  objFrac             = 0.0;
  float                  specFrac            = 0.0;

  double                 sumExpTime; 
  double                *appExpTime;
  double                *truExpTime;

  PilFrame              *grismFrame;
  PilFrame              *winFrame;
  PilFrame              *mosFrame, *sphotFrame;
  PilFrame              *outputFrame;
  PilFrame              *atmFrame;

  VimosImage            *winFile             = NULL;
  VimosImage            *grismFile           = NULL;
  VimosImage            *sphotFile           = NULL;
  VimosImage            *atmFile             = NULL;
  VimosImage            *objectFile          = NULL;
  VimosImage            *windowFile          = NULL;
  VimosImage           **imaSpEx1D           = NULL;
  VimosImage            *imaSpEx1DCal        = NULL;
  VimosImage           **mosList2D           = NULL;
  VimosImage           **outSpSkyExStack     = NULL;
  VimosImage            *copySpStack2D;
  VimosImage            *copyMos;
  VimosImage            *tmpImage;
  VimosTable            *grismTable          = NULL;
  VimosTable            *sphotTable          = NULL;
  VimosTable            *atmTable            = NULL;
  VimosTable           **winList             = NULL;

  VimosFloatArray       *offSets             = NULL;

  VimosExtractionTable  *extTable            = NULL;
  VimosWindowTable      *combWindowTable;
  VimosWindowTable     **winTablesList;
  VimosObjectTable      *objectTable;


 /* FIXME:
  * This parameter is still unused:
  */

  float          limFrac = 0.0;


 /*
  * Get task parameters from the recipe database
  */

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
                  "method: Average");
    combMethod = COMB_MEDIAN;
    minFrames = MIN_FRAMES_AVERAGE;
    break;
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

  specFrac = pilDfsDbGetDouble(parameter, "SpectrumFraction", 0.9);

  if (specFrac < 0 || specFrac > 1) {
    cpl_msg_error(task, "SpectrumFraction parameter must be > 0 and < 1");
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
   * Check if the spectro-photometric calibration should be applied.
   */

  calibrateFlux = pilDfsDbGetBool(parameter, "CalibrateFlux", 0);


  /*
   * Load input frames
   */

  cpl_msg_info(task, "Loading input frames...");


  if (calibrateFlux) {

    /*
     * Check that a SPH table is present in input frames.
     */

    sphotFrame = pilSofLookup(sof, sphotTableCategory);
    if (!sphotFrame) {
      cpl_msg_error(task, "Missing SPH Table: "
                    "spectrophotometric calibration cannot be applied.");
      return EXIT_FAILURE;
    }
    else
      pilFrmSetType(sphotFrame, PIL_FRAME_TYPE_CALIB);
  }


  /*
   * Check that an atmospheric extinction table is present in input frames.
   * It must be present in any case, because different observations
   * may be obtained at significantly different airmasses.
   */

  atmFrame = pilSofLookup(sof, atmTableCategory);
  if (!atmFrame) {
    cpl_msg_error(task, "Missing atmospheric extinction table: "
                  "input frames cannot be combined.");
    return EXIT_FAILURE;
  }
  else
    pilFrmSetType(atmFrame, PIL_FRAME_TYPE_CALIB);


  /*
   * Get 2D-extracted MOS frames
   */

  mosCount = pilSofFrameCount(sof, mosExtractedCategory);
  if (mosCount < minFrames) {
    cpl_msg_error(task, "At least %d %s exposures should be input for "
                  "specified combine method '%s'", minFrames, 
                  mosExtractedCategory, combMethodTag);
    return EXIT_FAILURE;
  }

  if ((mosList2D = (VimosImage **)cpl_calloc(mosCount, sizeof(VimosImage *)))) {

    mosFrame = pilSofLookupNext(sof, mosExtractedCategory);

    for (i = 0; i < mosCount; i++) {
      if ((mosList2D[i] = openOldFitsFile(pilFrmGetName(mosFrame), 1, 0))) {
        closeFitsImage(mosList2D[i], 0);
        pilFrmSetType(mosFrame, PIL_FRAME_TYPE_CALIB);
      }
      else {
        cpl_msg_error(task, "Failure opening %s image %d", mosExtractedCategory,
                      i + 1);
        for (j = 0; j < i; j++)
          deleteImage(mosList2D[j]);
        cpl_free(mosList2D);
        return EXIT_FAILURE;
      }
      mosFrame = pilSofLookupNext(sof, NULL);
    }
  }
  else {
    cpl_msg_error(task, "Failure creating list of input MOS exposures");
    return EXIT_FAILURE;
  }


  if (readIntDescriptor(mosList2D[0]->descs, pilTrnGetKeyword("Quadrant"),
                        &quadrant, NULL) == VM_FALSE) {
    cpl_msg_error(task, "Cannot read descriptor %s",
                  pilTrnGetKeyword("Quadrant"));
    for (i = 0; i < mosCount; i++)
      deleteImage(mosList2D[i]);
    cpl_free(mosList2D);
    return EXIT_FAILURE;
  }

  for (i = 0; i < mosCount; i++) {
    if (readIntDescriptor(mosList2D[i]->descs, 
                          pilTrnGetKeyword("MaskId", quadrant), &maskId, NULL)
        == VM_TRUE) {
      if (i) {
        if (firstMaskId == maskId) {
          continue;
        }
        cpl_msg_error(task, "Input mages do not come from the same mask");
      }
      else {
        firstMaskId = maskId;
        continue;
      }
    }
    else
      cpl_msg_error(task, "Missing descriptor %s in input images", 
                    pilTrnGetKeyword("MaskId", quadrant));
    for (i = 0; i < mosCount; i++)
      deleteImage(mosList2D[i]);
    cpl_free(mosList2D);
    return EXIT_FAILURE;
  }


  /*
   * Exposure times. If SummedExposureTime is not present, then
   * it's a product of vmmosobsstare, not vmmosobsjitter.
   * A distinction is made between different exposure times found
   * in header:
   *
   *   appExpTime[i] = Exposure time assigned to each frame to be
   *                   used for flux measurements. This is read
   *                   from the decriptor EXPTIME aka "ExposureTime".
   *   truExpTime[i] = Real exposure time for each input frame, to
   *                   be used for error computation purposes. This
   *                   is read from the descriptor ESO PRO EXPTTOT
   *                   aka "SummedExposureTime".
   *   sumExpTime    = The sum of all the truExpTime[i].
   *
   * If the "SummedExposureTime" is not present, appExpTime and 
   * truExpTime are coincident.
   */

  appExpTime = (double *)cpl_calloc(mosCount, sizeof(double));
  truExpTime = (double *)cpl_calloc(mosCount, sizeof(double));

  sumExpTime = 0.0;
  for (i = 0; i < mosCount; i++) {
    readDoubleDescriptor(mosList2D[i]->descs,
                         pilTrnGetKeyword("ExposureTime"),
                         appExpTime + i, NULL);
    if (readDoubleDescriptor(mosList2D[i]->descs, 
                             pilTrnGetKeyword("SummedExposureTime"),
                             truExpTime + i, NULL) == VM_FALSE) {
      truExpTime[i] = appExpTime[i];
    }
    sumExpTime += truExpTime[i];
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
    for (i = 0; i < mosCount; i++)
      deleteImage(mosList2D[i]);
    cpl_free(mosList2D);
    cpl_free(appExpTime);
    cpl_free(truExpTime);
    return EXIT_FAILURE;
  }


  /*
   * Create a working extraction table
   */

  if (!(extTable = VmSpExTab(mosList2D[0], grismTable, NULL, NULL))) { 
    cpl_msg_error(task, "Cannot compute extraction table");
    deleteTable(grismTable);
    for (j = 0; j < mosCount; j++)
      deleteImage(mosList2D[j]);
    cpl_free(mosList2D);
    cpl_free(appExpTime);
    cpl_free(truExpTime);
    return EXIT_FAILURE;
  }

  deleteTable(grismTable);


  /*
   * Load any window table specified in input
   */

  winCount = pilSofFrameCount(sof, windowTableCategory);

  if (winCount) {

    cpl_msg_info(task, "Loading %d input window tables...", winCount);

    winList = (VimosTable **)cpl_calloc(winCount, sizeof(VimosTable *));

    winFrame = pilSofLookupNext(sof, windowTableCategory);

    for (i = 0; i < winCount; i++) {
      if ((winFile = openOldFitsFile(pilFrmGetName(winFrame), 0, 0))) {
        int status = 0;
        pilFrmSetType(winFrame, PIL_FRAME_TYPE_CALIB);
        if (readDescsFromFitsImage(&(winFile->descs), winFile) == VM_FALSE) {
          cpl_msg_error(task, "readDescsFromFitsImage returned an error");
          return 0;
        }
        fits_movnam_hdu(winFile->fptr, BINARY_TBL, "WIN", 0, &status);
        winList[i] = newTable();
        if (readFitsTable(winList[i], winFile->fptr) == VM_TRUE) {
          closeFitsImage(winFile, 0);
          copyFromHeaderToHeader(winFile->descs, "DATE",
                                 &(winList[i]->descs), NULL);
        }
        else {
          cpl_msg_error(task, "Failure reading window table");
          for (j = 0; j < i; j++)
            deleteTable(winList[j]);
          cpl_free(winList);
          for (j = 0; j < mosCount; j++)
            deleteImage(mosList2D[j]);
          cpl_free(mosList2D);
          cpl_free(appExpTime);
          cpl_free(truExpTime);
          return EXIT_FAILURE;
        }
      }
      else
        cpl_msg_error(task, "Failure in opening window table");

      winFrame = pilSofLookupNext(sof, NULL);
    }
  }


  /*
   * Make all input images compatible with each others, and upgrade
   * the extraction table accordingly.
   */

  if (VmSpMatch(mosList2D, mosCount, extTable, 
                winList, winCount) == EXIT_FAILURE) {
    cpl_msg_error(task, "Cannot match input spectra");
    deleteExtractionTable(extTable);
    for (j = 0; j < mosCount; j++)
      deleteImage(mosList2D[j]);
    cpl_free(mosList2D);
    cpl_free(appExpTime);
    cpl_free(truExpTime);
    return EXIT_FAILURE;
  }


  cpl_msg_info(task, "Detecting objects in 2D-extracted spectra...");

  error = 0;

  if ((winTablesList = (VimosWindowTable **)cpl_calloc(mosCount, 
                        sizeof(VimosWindowTable *)))) {

    for (i = 0; i < mosCount; i++) {

      copyMos = duplicateImage(mosList2D[i]);
      copyAllDescriptors(mosList2D[i]->descs, &(copyMos)->descs);

      if (!(winTablesList[i] = VmSpDetObj(copyMos, extTable,
                                          numLevels, slitMargin, detLevel,
                                          objFrac, limFrac, minObjectSize,
                                          minCompositeSize, 0.0, specFrac))) {
        cpl_msg_error(task, "Failure deriving the window table for image %d",
                      i + 1);
        deleteImage(copyMos);
        for (j = 0; j < i; j++)
          deleteWindowTable(winTablesList[j]);
        cpl_free(winTablesList);
        error = 1;
        break;
      }
      deleteImage(copyMos);
    }
  }
  else {
    cpl_msg_error(task, "Failure creating window tables list");
    error = 1;
  }


  if (error) {
    deleteExtractionTable(extTable);
    for (j = 0; j < mosCount; j++)
      deleteImage(mosList2D[j]);
    cpl_free(mosList2D);
    cpl_free(appExpTime);
    cpl_free(truExpTime);
    return EXIT_FAILURE;
  }


  /*
   * Load the atmospheric extinction table and correct all spectra to
   * airmass zero.
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
    deleteExtractionTable(extTable);
    for (j = 0; j < mosCount; j++)
      deleteImage(mosList2D[j]);
    cpl_free(mosList2D);
    cpl_free(appExpTime);
    cpl_free(truExpTime);
    return EXIT_FAILURE;
  }

  for (i = 0; i < mosCount; i++) {
    tmpImage = VmSpApplyPhot(mosList2D[i], NULL, atmTable);
    if (!tmpImage) {
      deleteExtractionTable(extTable);
      deleteTable(atmTable);
      for (j = 0; j < mosCount; j++)
        deleteImage(mosList2D[j]);
      cpl_free(mosList2D);
      cpl_free(appExpTime);
      cpl_free(truExpTime);
      return EXIT_FAILURE;
    }
    deleteImage(mosList2D[i]);
    mosList2D[i] = tmpImage;
  }

  deleteTable(atmTable);


  /*
   *  Stack 2D images. The stacking is made keeping into account the
   *  different positions of the objects along the pseudo-slits. 
   *  The sky here is just an image padded with zeroes, since the
   *  sky was already removed by previous processing.
   *
   *  According to the stacking method applied, images are normalised
   *  in different ways. 
   *
   *  If the stacking method is "Average", the images are rescaled 
   *  according to their true exposure times, i.e. each image is 
   *  divided by its apparent exposure time and multiplied by its 
   *  true exposure time. In this way the "Average" stacking method 
   *  will be the weighted average of the input frames, with the 
   *  true exposure times used as weights, with the product image 
   *  having an apparent exposure time equal to the mean of the true 
   *  exposure times. In symbols:
   *
   *    F_i                       = counts in one second.
   *    T_i                       = true exposure time of frame i.
   *    F_i * T_i                 = counts in frame i after rescaling.
   *    N                         = number of frames to average.
   *    SUM(T_i * F_i) / N        = output of "Average" stacking method.
   *    SUM(T_i * F_i) / SUM(T_i) = average counts weighted with T_i values,
   *                                i.e., total counts per one second. 
   *    SUM(T_i) / N              = mean exposure time, i.e., ratio of 
   *                                the two previous items, QED.
   *
   *  If the stacking method is different from "Average", the images 
   *  are rescaled according to the mean true exposure time of all
   *  exposures. In this way the product image will have also in this
   *  case an apparent exposure time equal to the mean of the true 
   *  exposure times.
   */

  if (combMethod == COMB_AVERAGE) {
    for (i = 0; i < mosCount; i++) {
      constArithLocal(mosList2D[i], 
                      truExpTime[i]/appExpTime[i], VM_OPER_MUL);
    }
  }
  else {
    for (i = 0; i < mosCount; i++) {
      constArithLocal(mosList2D[i], 
                      sumExpTime/mosCount/appExpTime[i], VM_OPER_MUL);
    }
  }

  offSets = newFloatArray(mosCount);

  if (!(outSpSkyExStack = (VimosImage **)cpl_calloc(2, sizeof(VimosImage *)))) {
    cpl_msg_error(task, "Failure in creating 2D images array");
    deleteFloatArray(offSets);
    deleteExtractionTable(extTable);
    for (j = 0; j < mosCount; j++)
      deleteImage(mosList2D[j]);
    cpl_free(mosList2D);
    cpl_free(appExpTime);
    cpl_free(truExpTime);
    return EXIT_FAILURE;
  } 
  else {
    cpl_msg_info(task, "Stacking 2D images...");
    if (!(outSpSkyExStack[0] = VmSpStack2D(mosList2D, winTablesList, 
                                           extTable, NULL, 
                                           mosCount, combMethod, 
                                           &combParameter, offSets, 
                                           1, 0, 0))) {
      cpl_msg_error(task, "Failure in creating combined 2D image");
      deleteFloatArray(offSets);
      deleteExtractionTable(extTable);
      for (j = 0; j < mosCount; j++)
        deleteImage(mosList2D[j]);
      cpl_free(mosList2D);
      cpl_free(outSpSkyExStack);
      return EXIT_FAILURE;
    } 

    outSpSkyExStack[1] = newImageAndAlloc(outSpSkyExStack[0]->xlen, 
                                          outSpSkyExStack[0]->ylen);
  }

  deleteFloatArray(offSets);

  for (j = 0; j < mosCount; j++)
    deleteImage(mosList2D[j]);
  cpl_free(mosList2D);

  for (j = 0; j < mosCount; j++)
    deleteWindowTable(winTablesList[j]);
  cpl_free(winTablesList);


  /* 
   *  Detection on total image.
   */

  copySpStack2D = duplicateImage(outSpSkyExStack[0]);
  copyAllDescriptors(outSpSkyExStack[0]->descs, &(copySpStack2D)->descs);

  if (!(combWindowTable = VmSpDetObj(copySpStack2D, extTable, 
                                     numLevels, slitMargin, detLevel, 
                                     objFrac, limFrac, minObjectSize, 
                                     minCompositeSize, 0.0, specFrac))) {
    cpl_msg_error(task, "Failure deriving the window table for combined image");
    deleteImage(copySpStack2D);
    deleteImage(outSpSkyExStack[0]);
    deleteImage(outSpSkyExStack[1]);
    deleteExtractionTable(extTable);
    cpl_free(outSpSkyExStack);
    return EXIT_FAILURE;
  }

  vimosDscCopy(&combWindowTable->descs, copySpStack2D->descs,
               pilTrnGetKeyword("PROG.ID"), NULL);

  vimosDscCopy(&combWindowTable->descs, copySpStack2D->descs,
               pilTrnGetKeyword("OBS.ID"), NULL);

  deleteExtractionTable(extTable);
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
    return EXIT_FAILURE;
  }

  error = 1;
  cpl_msg_info(task, "Extracting 1D spectra...");
  if ((objectTable = newObjectTable())) {
    if ((imaSpEx1D = VmSpEx1D(outSpSkyExStack, combWindowTable, 
                              objectTable, 0, mosCount))) { 
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
    return EXIT_FAILURE;
  }

  deleteImage(imaSpEx1D[1]);
  deleteImage(outSpSkyExStack[1]);


  /* 
   * Spectrophotometric calibration to 1D spectra
   */ 

  insertDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("ExposureTime"),
                         sumExpTime/mosCount,
                         pilTrnGetComment("ExposureTime"),
                         "ESO*", 1);

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
      cpl_free(imaSpEx1D);
      deleteImage(outSpSkyExStack[0]);
      cpl_free(outSpSkyExStack);
      deleteObjectTable(objectTable);
      deleteWindowTable(combWindowTable);
      return EXIT_FAILURE;
    }

    imaSpEx1DCal = VmSpApplyPhot(imaSpEx1D[0], sphotTable, NULL);

    if (imaSpEx1DCal == NULL) {
      cpl_msg_warning(task, "Spectro-photometric calibration failure: "
                      "the corresponding product, %s, will not be created", 
                      mosFluxReducedCategory);
      calibrateFlux = 0;
    }

    deleteTable(sphotTable);

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
  
  updateOK = updateOK && writeDoubleDescriptor(&(imaSpEx1D[0]->descs),
                         pilTrnGetKeyword("AirMass"), 0.0,
                         pilTrnGetComment("AirMass"));
  
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
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[0]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
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
    cpl_free(imaSpEx1D);
    deleteImage(outSpSkyExStack[0]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    return EXIT_FAILURE;
  }
  
  deleteImage(imaSpEx1D[0]);
  cpl_free(imaSpEx1D);

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
      deleteImage(outSpSkyExStack[0]);
      cpl_free(outSpSkyExStack);
      deleteObjectTable(objectTable);
      deleteWindowTable(combWindowTable);
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
      deleteImage(outSpSkyExStack[0]);
      cpl_free(outSpSkyExStack);
      deleteObjectTable(objectTable);
      deleteWindowTable(combWindowTable);
      return EXIT_FAILURE;
    }

    deleteImage(imaSpEx1DCal);

  }


  updateOK = updateOK && insertDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("ExposureTime"),
                         sumExpTime/mosCount,
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
  
  updateOK = updateOK && writeDoubleDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("AirMass"), 0.0,
                         pilTrnGetComment("AirMass"));
  
/*
  updateOK = updateOK && insertStringDescriptor(&(outSpSkyExStack[0]->descs),
                         pilTrnGetKeyword("DoCategory"), 
                         mosExtractedCategory,
                         pilTrnGetComment("DoCategory"), 
                         "ESO*", 1);
*/
  
  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product 2 header");
    deleteImage(outSpSkyExStack[0]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
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
    deleteImage(outSpSkyExStack[0]);
    cpl_free(outSpSkyExStack);
    deleteObjectTable(objectTable);
    deleteWindowTable(combWindowTable);
    return EXIT_FAILURE;
  }

  deleteImage(outSpSkyExStack[0]);
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
    deleteWindowTable(combWindowTable);
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

  if (error)
    return EXIT_FAILURE;

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
                    "vmmoscombine",
    "Combine reduced MOS observations from different OBs.",
    "This recipe is used to sum the contributes from a sequence of\n"
    "2D-extracted spectral frames generated by the recipes vmmosobsstare,\n"
    "vmmosobsjitter, and vmmoscombine itself. The only restriction is\n"
    "that all frames must be obtained with the same mask. Sky subtraction,\n" 
    "fringing correction, flat fielding, etc., are not applied, since they\n" 
    "had been applied in the previous reduction of each input dataset.\n" 
    "Each input image is corrected to airmass zero before the stacking,\n"
    "and for this reason an atmospheric extinction table must always be\n"
    "specified. Optionally a correction for the instrument response can\n" 
    "also be applied.\n\n"
    "Input files:\n\n"
    "  DO category:             Type:       Explanation:          Required:\n"
    "  MOS_SCIENCE_EXTRACTED    Product     Combined slit spectra     Y\n"
    "  GRISM_TABLE              Calib       Grism table               Y\n"
    "  EXTINCT_TABLE            Calib       Atmospheric extinction    Y\n"
    "  MOS_SPECPHOT_TABLE       Calib       Response curve            .\n\n"
    "Output files:\n\n"
    "  DO category:             Data type:  Explanation:\n"
    "  MOS_SCIENCE_REDUCED      FITS image  Extracted objects spectra\n"
    "  MOS_SCIENCE_FLUX_REDUCED FITS image  Flux calibrated objects spectra\n"
    "  MOS_SCIENCE_EXTRACTED    FITS image  Combined slit spectra\n"
    "  OBJECT_TABLE             FITS table  Objects spectra identification\n"
    "  WINDOW_TABLE             FITS table  Objects positions in slit\n\n"
    "The positions of the extracted slit spectra and of the detected\n"
    "objects that they may contain are listed in the window table.\n\n"
    "If a spectro-photometric table (produced by the recipe vmmosstandard)\n"
    "is specified and a flux calibration is requested, then a\n" 
    "MOS_SCIENCE_FLUX_REDUCED image is also created. This image is identical\n" 
    "to the MOS_SCIENCE_REDUCED, but the spectra it contains are flux\n" 
    "calibrated, and expressed in units of erg/cm/cm/s/Angstrom.\n\n"
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

                    vmmoscombine_create,
                    vmmoscombine_exec,
                    vmmoscombine_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
