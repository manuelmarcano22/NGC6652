/* $Id: vmspcaldisp.c,v 1.5 2013-08-07 16:45:33 cgarcia Exp $
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
#include "vmlinecatalog.h"
#include "vmextractiontable.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmmosflat.h"
#include "vmmosmodels.h"
#include "vmmoswavecalib.h"
#include "vmmosextraction.h"
#include "vmmosutils.h"
#include "vmqcutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


static cxint vmspcaldisp(PilSetOfFrames *);


#define MAX_COMMENT_LENGTH  (80)
#define NSIGMA (3)


/*
 * Definition of the label strings for all methods the recipe function
 * supports for bias removal, with their associated method code.
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
 * supports for 2D extraction of the arc lamp frame, after the IDS has
 * been computed, with their associated method code.
 */

static const char *extrMethodNames[] = {
  "Local",
  "Global",
  "Input"
};

static const ExtrMethod extrMethods[] = {
  EXTR_LOCAL,
  EXTR_GLOBAL,
  EXTR_INPUT
};

static unsigned int nExtrMethods = sizeof(extrMethods) / sizeof(ExtrMethod);


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it availble to the application using the interface.
 */

static cxint
vmspcaldisp_create(cpl_plugin *plugin)
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


    p = cpl_parameter_new_value("vimos.Parameters.extraction.fuzz",
                                CPL_TYPE_INT,
                                "Extra pixels from expected position of "
                                "spectrum edge in spectral extraction.",
                                "vimos.Parameters",
                                10);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "Fuzz");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "Fuzz");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.lamp.frames.validate",
                                CPL_TYPE_BOOL,
                                "Exposure check on input arc lamp frames.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ValidateFrames");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ValidateFrames");
    cpl_parameterlist_append(recipe->parameters, p);

#endif


    p = cpl_parameter_new_value("vimos.Parameters.ids.refine",
                                CPL_TYPE_BOOL,
                                "Refine Inverse Dispersion Solution.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "RefineIDS");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "RefineIDS");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.extraction.window",
                                CPL_TYPE_INT,
                                "Size of search window around expected "
                                "arc line positions.",
                                "vimos.Parameters",
                                10);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ExtractionWindow");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ExtractionWindow");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flat.apply",
                                CPL_TYPE_BOOL,
                                "Flat field correction for input arc "
                                "lamp frames.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyFlatField");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyFlatField");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.badpixel.clean",
                                CPL_TYPE_BOOL,
                                "Bad pixel correction on arc lamp image.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.clean",
                                CPL_TYPE_BOOL,
                                "Cosmic ray events cleaning in arc lamp "
                                "image.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanCosmic");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanCosmic");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.line.ident.method",
                               CPL_TYPE_STRING,
                               "Arc line identification method.",
                               "vimos.Parameters",
                               "FirstGuess", 2, "FirstGuess", "Blind");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LineIdent");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LineIdent");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.line.ident.level",
                                CPL_TYPE_DOUBLE,
                                "Threshold for peak detection.",
                                "vimos.Parameters",
                                500.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LineIdentLevel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LineIdentLevel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.extraction",
                                CPL_TYPE_STRING,
                                "Arc lamp 2D extraction according to "
                                "computed Local or Global IDS.",
                                "vimos.Parameters",
                                "Local", 3, "Local", "Global", "Input");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ArcExtraction");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ArcExtraction");
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
vmspcaldisp_exec(cpl_plugin *plugin)
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

    if (vmspcaldisp(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmspcaldisp");
        
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
vmspcaldisp_destroy(cpl_plugin *plugin)
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
 *   Upgrade the inverse dispersion solution identifying lines of known
 *   wavelength on an arc lamp frame.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param sof  Pointer to a list of input frames containing references
 *             to arc lamp frames (just one for each mask shutter position), 
 *             master bias, optional master dark, optional master spectral 
 *             flat field, line catalog, CCD table, grism table, and (in 
 *             case) IFU table. 
 *
 * @doc
 *   This recipe derives an improved wavelength calibration identifying 
 *   on the stacked arc lamp frame the lines whose wavelengths are listed 
 *   in the input catalog. The optical distortion model components are
 *   also computed, and if an input master flat field is given also the
 *   spectral curvature model is improved. Bias is preliminarly removed 
 *   from the raw arc lamp frames, dark is optionally subtracted, flat 
 *   fielding, cosmic rays removal and bad pixel correction are an option. 
 *   An IFU table must be specified in case of an IFU observation. All 
 *   the input frames are passed via the set of frames \textbf{sof}. 
 *   On successful termination the global and local solutions are written 
 *   to a new extraction table, and all the global distortion models are
 *   written to a PAF file.
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
 *
 *   \begin{itemize}
 *
 *     \item BiasMethod:         Method for bias removal from arc lamp
 *                               frame. Legal settings are:
 *
 *     \begin{itemize}
 *
 *       \item Master:             Prescan and overscan regions are trimmed
 *                                 away from the arc lamp image after the
 *                                 master bias removal.
 *
 *       \item Zmaster:            After master bias removal the residual
 *                                 signal found in the arc lamp image
 *                                 overscan regions is modelled and
 *                                 subtracted. Next, prescan and overscan 
 *                                 regions are trimmed away.
 *     \end{itemize}
 *
 *     \item Fuzz:               Extra number of pixels to be used in
 *                               spectra extraction, when searching for
 *                               a spectrum edge on the CCD. This value
 *                               is ignored if an input master spectral
 *                               flat field is not found in input.
 *
 *     \item ValidateFrames:     Exposure check on raw arc lamp frames.
 *
 *     \item RefineIDS           Flag to signal the need to compute a
 *                               refined IDS with respect to what is
 *                               stored in the Extraction Table
 *
 *     \item ExtractionWindow:   Size of extraction window centered on 
 *                               the expected position of any arc lamp 
 *                               line.
 *
 *     \item CleanBadPixel:      Bad pixel correction on arc lamp frame.
 *                               If this option is turned on, the recipe 
 *                               expects a CCD\_TABLE in the input set 
 *                               of frames.
 *
 *     \item CleanCosmic:        Cosmic ray events removal from input arc
 *                               lamp frame. If a CCD\_TABLE is found in 
 *                               the input set of frames, bad pixels will 
 *                               not be used in computing the values to 
 *                               replace the cosmic rays events.
 *
 *     \item ApplyFlatField:     Flat field correction on arc lamp frame.
 *                               If this option is turned on, the recipe
 *                               expects a MOS\_MASTER\_SCREEN\_FLAT in
 *                               the input set of frames.
 *
 *   \end{itemize}
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 *
 * @author C. Izzo
 */   

static cxint 
vmspcaldisp(PilSetOfFrames *sof)
{

  const char   task[] = "vmspcaldisp";

  const char   parameter[]         = "Parameters";

  const char  *mosMasterCategory   = pilTrnGetCategory("MosMasterScreenFlat");
  const char  *mosCombCategory     = pilTrnGetCategory("MosCombScreenFlat");
  const char  *ifuMasterCategory   = pilTrnGetCategory("IfuMasterScreenFlat");
  const char  *masterBiasCategory  = pilTrnGetCategory("MasterBias");
  const char  *masterDarkCategory  = pilTrnGetCategory("MasterDark");
  const char  *ccdTableCategory    = pilTrnGetCategory("CcdTable");
  const char  *grismTableCategory  = pilTrnGetCategory("GrismTable");
  const char  *extrTableCategory   = pilTrnGetCategory("ExtractTable");
  const char  *ifuTableCategory    = pilTrnGetCategory("IfuTable");
  const char  *lineCatalogCategory = pilTrnGetCategory("LineCatalog");
  const char  *mosArcLampCategory  = pilTrnGetCategory("MosArcSpectrum");
  const char  *ifuArcLampCategory  = pilTrnGetCategory("IfuArcSpectrum");
  const char  *pafFileCategory     = pilTrnGetCategory("PAFCategory");
  const char  *arc2DCatg       = pilTrnGetCategory("MosArcSpectrumExtracted");

  char          *arcLampCategory, *masterCategory;

  char          *biasMethodTag = NULL;
  BiasMethod     biasMethod    = BIAS_UNDEF;
  int            biasMethodEntry;

  char          *extrMethodTag = NULL;
  ExtrMethod     extrMethod    = EXTR_UNDEF;
  int            extrMethodEntry;

  char          *lineIdentificationTag = NULL;
  unsigned int   lineIdentification;

  char          *pafFileName   = NULL;
  char           namePAF[]     = "MOS_wavecal";
  char           pipeNamePAF[] = "MOS_wavecal_";

  char           outputExtrName[PATHNAME_MAX + 1];
  char           extrArcName[PATHNAME_MAX + 1];
  char           comment[MAX_COMMENT_LENGTH];
  char           parName[30];
  char           unit[20];
  char           dateObs[25];
  char           grismName[10];

  size_t         numFlat;
  size_t         numComb;
  size_t         numArcMos = 0;
  size_t         numArcIfu = 0;
  size_t         lampCount;

  unsigned int   h, i, j, k;
  unsigned int   validateFrames, cleanBadPixel, cleanCosmic, applyFlat;
  unsigned int   computeQC;
  unsigned int   error;
  unsigned int   windowDimen, numRows;
  unsigned int   refineIDS;
  unsigned int   modelSlit;

  int            order = 0;
  int            fitOpt, fitCrv, fitIds, maxOrdOpt;
  int            fuzz;
  int            arcWidth;
  int            X, prX;
  int            Y, prY;
  int            goodSlitX, goodSlitY;
  int            numSlits;
  int            startX, startY, xmax;
  int            quadrant;

  float         *lampRon = NULL;
  float         *lampGain = NULL;

  float          posH, posL, mShuMedPos, linePosFloat;
  float          newDistance, distance;
  float         *extrLine;
  float         *lineShape;
  float         *buffer;
  float          threshold, xmaxFloat;
  float          lineFlux, lineFluxRem, tmpPos;
  float          lambdaHe, lambdaNe, lambdaAr;
  float          lambdaRed, lambdaYel, lambdaBlu;
  float          lambdaCen, pixelCen;
  float          level;

  double         focalScale;
  double         temperature, time;
  double         sumExposures = 0.0;
  double         meanValue, medianValue, rmsValue;
  double         lflux, lflux_err;
  double         meanFwhm = 0.0;
  double         cdelt, crval;
  int            npix, idsOrd;
  int            count;

  int            idsOrdX, idsOrdY;
  int            crvOrdX, crvOrdY;
  int            optOrdX, optOrdY;
  double         optRmsX, optRmsY;
  double         x, y, minSlitx, maxSlitx, minSlity, maxSlity;
  double         slitLength;
  double         deltaX, deltaY;

  PilFrame      *masterFlatFrame, *grismFrame, *ifuFrame, *outputFrame;
  PilFrame      *arcLampFrame, *lineCatFrame, *biasFrame, *darkFrame;
  PilFrame      *ccdFrame, *combFlatFrame;
  PilFrame      *extrFrame = NULL;

  VimosImage    *masterFlatImage = NULL;
  VimosImage    *combFlatImage = NULL;
  VimosImage    *biasImage, *biasOver;
  VimosImage    *darkImage = NULL;
  VimosImage    *arcLampImage;
  VimosImage    *grismFile, *lineCatFile;
  VimosImage    *ifuFile = NULL;
  VimosImage    *extrFile = NULL;
  VimosImage   **lampList = NULL;
  VimosImage   **stackImages = NULL;
  VimosImage   **extr2D = NULL;
  VimosImage    *extracted2D = NULL;
  VimosImage   **imabuf = NULL;
  VimosImage    *tmpImage = NULL;
  VimosTable    *lineCat = NULL;
  VimosTable    *ccdTable = NULL;
  VimosTable    *grismTable = NULL;
  VimosIfuTable *ifuTable = NULL;
  VimosPixel    *coefArrayX = NULL;
  VimosPixel    *coefArrayY = NULL;

  VimosExtractionTable *extractionTable = NULL;
  VimosExtractionTable *useExtractTable;
  VimosExtractionSlit  *slit;
  VimosExtractionSlit  *goodSlit;
  VimosColumn          *wLen;
  VimosColumn          *flux;

  VimosDistModel2D     *optModX = NULL;
  VimosDistModel2D     *optModY = NULL;
  VimosDistModelFull   *crvMod;
  VimosDistModelFull   *idsMat;


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
  * Determine the frame 2D extraction method.
  */

  extrMethodTag = (char *)pilDfsDbGetString(parameter, "ArcExtraction");

  if ((extrMethodEntry =
        strselect(extrMethodTag, extrMethodNames, nExtrMethods)) < 0) {
    cpl_msg_error(task, "%s: Invalid 2d extraction method.", extrMethodTag);
    return EXIT_FAILURE;
  }

  extrMethod = extrMethods[extrMethodEntry];


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
   * Get "refineIDS" parameter (flag to signal the need to refine the IDS
   * with respect to what is sored in the Extraction Table)
   */

  refineIDS = pilDfsDbGetBool(parameter, "RefineIDS", 0);

 /*
  * Get "arcWidth" parameter (search window around expected arc line
  * position).
  */

  arcWidth = pilDfsDbGetInt(parameter, "ExtractionWindow", -1);


 /*
  * Check if the bad pixels should be corrected.
  */

  cleanBadPixel = pilDfsDbGetBool(parameter, "CleanBadPixel", 0);


 /*
  * Check if raw frames with non-nominal exposure level should be ignored
  * by the task.
  */

  validateFrames = pilDfsDbGetBool(parameter, "ValidateFrames", 0);


 /*
  * Check if the cosmic rays events should be corrected.
  */

  cleanCosmic = pilDfsDbGetBool(parameter, "CleanCosmic", 0);


 /*
  * Check if the flat fielding correction should be applied.
  */

  applyFlat = pilDfsDbGetBool(parameter, "ApplyFlatField", 0);


  /*
   * Optional modeling of wavelength solution within each slit
   */

  modelSlit = pilDfsDbGetBool(parameter, "ModelSlit", 0);


  /*
   * Order of polynomial for modeling of wavelength solution within each slit.
   */

  order = pilDfsDbGetInt(parameter, "ModelSlitOrder", 0);

  if (order < 0) {
    cpl_msg_error(task, "Invalid order for wavelength solution modeling.");
    return EXIT_FAILURE;
  }


  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool(parameter, "ComputeQC", 1);


  /*
   * Get line identification method
   */

  lineIdentificationTag = (char *)pilDfsDbGetString(parameter, "LineIdent");

  if (strcmp(lineIdentificationTag, "FirstGuess") == 0)
    lineIdentification = 0;
  else {
    lineIdentification = 1;
    level = pilDfsDbGetDouble(parameter, "LineIdentLevel", 500.);
  }


  /*
   * Load input frames
   */

  cpl_msg_info(task, "Loading input frames...");


 /*
  * If bad pixel correction is enabled, a bad pixel table is required
  * in the input set of frames. If no bad pixel table is present this
  * is an error. If cleaning of cosmic rays is enabled, the bad pixel
  * table will be used, if present, to avoid bad pixels while smoothing
  * away cosmic ray events. In this case, a missing bad pixel table is
  * not an error.
  */

  error = 0;

  ccdFrame = pilSofLookup(sof, ccdTableCategory);
  if (ccdFrame)
    pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

  if (cleanBadPixel || cleanCosmic) {
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
    else {
      if (cleanBadPixel)
        cpl_msg_error(task, "No CCD table in input: cannot clean bad pixels.");
      else
        error = 0;
    }
  }

  if (error)
    return EXIT_FAILURE;

 /*
  * The grism table is required.
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
    }
    else
      cpl_msg_error(task, "Failure in opening grism table");
  }
  else
    cpl_msg_error(task,"No input grism table found");

  if (error) {
    deleteTable(ccdTable);
    return EXIT_FAILURE;
  }


 /*
  * Load the raw arc lamp exposures (either IFU or MOS)
  */

  numArcMos = pilSofFrameCount(sof, mosArcLampCategory);
  numArcIfu = pilSofFrameCount(sof, ifuArcLampCategory);

  if ((numArcMos + numArcIfu) == 0) {
    cpl_msg_info(task, "No arc lamp frame found in input!");
    deleteTable(ccdTable);
    deleteTable(grismTable);
    return EXIT_FAILURE;
  }

  if (numArcMos && numArcIfu) {
    cpl_msg_error(task, "Both IFU and MOS arc lamp frames found in input!");
    deleteTable(ccdTable);
    deleteTable(grismTable);
    return EXIT_FAILURE;
  }

  if (numArcMos) {
    lampCount = numArcMos;
    arcLampCategory = (char *)mosArcLampCategory;
  }

  if (numArcIfu) {
    lampCount = numArcIfu;
    arcLampCategory = (char *)ifuArcLampCategory;
  }

  error = 1;

  if ((lampList = (VimosImage **)cpl_calloc(lampCount, sizeof(VimosImage *)))) {

    error = 0;
    arcLampFrame = pilSofLookupNext(sof, arcLampCategory);

    for (i = 0; i < lampCount; i++) {
      if ((lampList[i] = openOldFitsFile(pilFrmGetName(arcLampFrame), 1, 0))) {
        pilFrmSetType(arcLampFrame, PIL_FRAME_TYPE_RAW);
        closeFitsImage(lampList[i], 0);

        readDoubleDescriptor(lampList[i]->descs,
                             pilTrnGetKeyword("ExposureTime"), &time, NULL);
        sumExposures += time;
      }
      else {
        error = 1;
        cpl_msg_error(task, "Failure opening arc lamp frame %d", i + 1);
        for (j = 0; j < i; j++)
          deleteImage(lampList[j]);
        cpl_free(lampList);
        break;
      }
      arcLampFrame = pilSofLookupNext(sof, NULL);
    }
  }
  else
    cpl_msg_error(task, "Failure creating list of input arc lamp exposures");

  if (error) {
    deleteTable(ccdTable);
    deleteTable(grismTable);
    return EXIT_FAILURE;
  }

  if (readDoubleDescriptor(lampList[0]->descs, pilTrnGetKeyword("FocalScale"),
                           &focalScale, NULL) == VM_FALSE) {
    focalScale = 1.718;
    cpl_msg_warning(task, "Cannot read descriptor %s (defaulted to %f)",
                  pilTrnGetKeyword("FocalScale"), focalScale);
  }

  readIntDescriptor(lampList[0]->descs, pilTrnGetKeyword("Quadrant"), 
                    &quadrant, NULL);

  readDoubleDescriptor(lampList[0]->descs,
                       pilTrnGetKeyword("BeamTemperature", quadrant),
                       &temperature, NULL);

  readStringDescriptor(lampList[0]->descs, pilTrnGetKeyword("DateObs"), 
                       dateObs, NULL);

  readStringDescriptor(lampList[0]->descs, 
                       pilTrnGetKeyword("GrismName", quadrant), 
                       grismName, NULL);


 /*
  * The IFU table should be present for IFU data.
  */

  ifuFrame = pilSofLookup(sof, ifuTableCategory);
  if (ifuFrame)
    pilFrmSetType(ifuFrame, PIL_FRAME_TYPE_CALIB);

  if (numArcIfu) {

    error = 1;

    if (ifuFrame) {
      if ((ifuFile = openOldFitsFile(pilFrmGetName(ifuFrame), 0, 0))) {
        if ((ifuTable = newIfuTable())) {
          if (readFitsIfuTable(ifuTable, ifuFile->fptr) == VM_TRUE) {
            closeFitsImage(ifuFile, 0);
            error = 0;
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
      cpl_msg_error(task, "No input IFU table found");

    if (error) {
      deleteTable(ccdTable);
      deleteTable(grismTable);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
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
    for (i = 0; i < lampCount; i++)
      deleteImage(lampList[i]);
    cpl_free(lampList);
    return EXIT_FAILURE;
  }


 /*
  * Recreate bias overscans using as a reference the arc lamp frame
  */

  if ((biasOver = growOverscans(biasImage, lampList[0]))) {
    if (biasImage != biasOver) {
      deleteImage(biasImage);
      biasImage = biasOver;
    }
  }
  else {
    cpl_msg_error(task, "Failure in growing overscans in master bias");
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    deleteImage(biasImage);
    for (i = 0; i < lampCount; i++)
      deleteImage(lampList[i]);
    cpl_free(lampList);
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
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteImage(biasImage);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      return EXIT_FAILURE;
    }
  }
  else
    cpl_msg_warning(task, "No master dark in input, dark subtraction "
                  "will not be performed");


 /*
  *  Master Spectral Flat Field (mandatory only if applyFlat == true).
  */

  if (numArcMos)
    masterCategory = (char *)mosMasterCategory;
  else
    masterCategory = (char *)ifuMasterCategory;

  numFlat = pilSofFrameCount(sof, masterCategory);

  if (numFlat > 1) {
    cpl_msg_error(task,
                "More than one master spectral flat field found in input");
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    deleteImage(biasImage);
    deleteImage(darkImage);
    for (i = 0; i < lampCount; i++)
      deleteImage(lampList[i]);
    cpl_free(lampList);
    return EXIT_FAILURE;
  }


  error = 1;

  if (numFlat) {
    masterFlatFrame = pilSofLookup(sof, masterCategory);
    pilFrmSetType(masterFlatFrame, PIL_FRAME_TYPE_CALIB);
    if ((masterFlatImage = openOldFitsFile(pilFrmGetName(masterFlatFrame), 
                                           1, 0))) {
      closeFitsImage(masterFlatImage, 0);
      error = 0;
    }
    else
      cpl_msg_error(task, "Failure opening master spectral flat image %s", 
                  pilFrmGetName(masterFlatFrame));
  }
  else 
    if (applyFlat)
      cpl_msg_error(task, "Missing master spectral flat image");
    else
      error = 0;

  if (error) {
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    deleteImage(biasImage);
    deleteImage(darkImage);
    for (i = 0; i < lampCount; i++)
      deleteImage(lampList[i]);
    cpl_free(lampList);
    return EXIT_FAILURE;
  }


  if (numArcMos) {
    numComb = pilSofFrameCount(sof, mosCombCategory);

    if (numComb > 1) {
      cpl_msg_error(task, "More than one combined spectral flat field "
                  "found in input");
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteImage(biasImage);
      deleteImage(darkImage);
      deleteImage(masterFlatImage);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      return EXIT_FAILURE;
    }


    error = 1;

    if (numComb) {
      combFlatFrame = pilSofLookup(sof, mosCombCategory);
      pilFrmSetType(combFlatFrame, PIL_FRAME_TYPE_CALIB);
      if ((combFlatImage = openOldFitsFile(pilFrmGetName(combFlatFrame),
                                           1, 0))) {
        closeFitsImage(combFlatImage, 0);
        error = 0;
      }
      else
        cpl_msg_error(task, "Failure opening combined spectral flat image %s",
                    pilFrmGetName(combFlatFrame));
    }

  }


 /*
  *  Line Catalog
  */

  error = 1;

  if ((lineCatFrame = pilSofLookup(sof, lineCatalogCategory))) {
    pilFrmSetType(lineCatFrame, PIL_FRAME_TYPE_CALIB);
    if ((lineCatFile = openOldFitsFile(pilFrmGetName(lineCatFrame), 0, 0))) {
      if ((lineCat = newLineCatalog())) {
        if (readFitsLineCatalog(lineCat, lineCatFile->fptr) == VM_TRUE) {
          closeFitsImage(lineCatFile, 0);
          error = 0;
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
  else
    cpl_msg_error(task,"No input line catalog found");

  if (error) {
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteImage(masterFlatImage);
    deleteImage(combFlatImage);
    for (i = 0; i < lampCount; i++)
      deleteImage(lampList[i]);
    cpl_free(lampList);
    return EXIT_FAILURE;
  }


 /*
  * Noise contributions for each input arc lamp frame must be estimated 
  * before any processing (as bias subtraction, normalization, etc.).
  */

  if ((lampRon = (float *)cpl_calloc(lampCount, sizeof(float)))) {
    for (i = 0; i < lampCount; i++) {
      lampRon[i] = computeAverageRon(lampList[i]);
    }
  }
  else {
    cpl_msg_error(task, "Not enough memory");
    deleteTable(lineCat);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    deleteImage(biasImage);
    deleteImage(darkImage);
    deleteImage(masterFlatImage);
    deleteImage(combFlatImage);
    for (i = 0; i < lampCount; i++)
      deleteImage(lampList[i]);
    cpl_free(lampList);
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

  for (i = 0; i < lampCount; i++) {
    if (VmSubBias(lampList[i], biasImage, biasMethod) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot remove bias from arc lamp frame %d", i + 1);
      deleteTable(lineCat);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteImage(biasImage);
      deleteImage(darkImage);
      deleteImage(masterFlatImage);
      deleteImage(combFlatImage);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      cpl_free(lampRon);
      return EXIT_FAILURE;
    }
  }

  deleteImage(biasImage);


 /*
  * Now go for the (optional) dark subtraction
  */

  if (darkImage) {
    cpl_msg_info(task, "Dark subtraction...");

    for (i = 0; i < lampCount; i++) {
      if (VmSubDark(lampList[i], darkImage) == EXIT_FAILURE) {
        cpl_msg_error(task, 
                    "Cannot subtract dark from arc lamp frame %d", i + 1);
        deleteTable(lineCat);
        deleteTable(ccdTable);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        deleteImage(darkImage);
        deleteImage(masterFlatImage);
        deleteImage(combFlatImage);
        for (i = 0; i < lampCount; i++)
          deleteImage(lampList[i]);
        cpl_free(lampList);
        cpl_free(lampRon);
        return EXIT_FAILURE;
      }
    }
    deleteImage(darkImage);
  }


 /*
  *  Optional image cleaning
  */

  if (cleanCosmic) {

   /*
    * Compute gain for arc lamp frames for cosmic rays cleaning
    * (RON is already computed, and median levels are computed 
    * within VmCosmicClean() ).
    */

    error = 1;
    if ((lampGain = (float *)cpl_calloc(lampCount, sizeof(float)))) {
      error = 0;
      for (i = 0; i < lampCount; i++) {
        lampGain[i] = getMeanGainFactor(lampList[i]);
        if (lampGain[i] < MIN_DIVISOR) {
          error = 1;
          cpl_msg_error(task, 
                      "Wrong or no gain factor for arc lamp frame %d", i + 1);
          cpl_free(lampGain);
          break;
        }
      }
    }
    else
      cpl_msg_error(task, "Not enough memory");

    if (error) {
      deleteTable(lineCat);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteImage(masterFlatImage);
      deleteImage(combFlatImage);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      cpl_free(lampRon);
      return EXIT_FAILURE;
    }


   /*
    * Cosmic rays cleaning
    */

    cpl_msg_info(task, "Cleaning cosmic rays in arc lamp frames.");

    for (i = 0; i < lampCount; i++) {

      if (EXIT_FAILURE == VmCosmicClean(lampList[i], ccdTable, 0, -1., 
                                        lampGain[i], lampRon[i], -1., -1.)) {
        deleteTable(lineCat);
        deleteTable(ccdTable);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        deleteImage(masterFlatImage);
        deleteImage(combFlatImage);
        for (i = 0; i < lampCount; i++)
          deleteImage(lampList[i]);
        cpl_free(lampList);
        cpl_free(lampGain);
        cpl_free(lampRon);
        return EXIT_FAILURE;
      }
      cpl_msg_info(task, "  Arc lamp %d of %zd done", i + 1, lampCount);
    }

    cpl_free(lampGain);

  }

  if ((extrFrame = pilSofLookup(sof, extrTableCategory))) {

    cpl_msg_info(task, "An Extraction Table was found in input: using as "
               "first-guesses the spectral distortion models found there, "
               "instead of the distortion models from the input image.");

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
    deleteTable(lineCat);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    deleteImage(masterFlatImage);
    deleteImage(combFlatImage);
    for (i = 0; i < lampCount; i++)
      deleteImage(lampList[i]);
    cpl_free(lampList);
    cpl_free(lampGain);
    cpl_free(lampRon);
    return EXIT_FAILURE;
  }

  if (masterFlatImage) {
   
    if (extractionTable) {
      readOptDistModel(extractionTable->descs, &optModX, &optModY);
      writeOptDistModelString(&(masterFlatImage->descs), optModX, optModY);
      readCurvatureModel(extractionTable->descs, &crvMod);
      writeCurvatureModelString(&(masterFlatImage->descs), crvMod);
      readInvDispMatrix(extractionTable->descs, &idsMat);
      writeInvDispMatrixString(&(masterFlatImage->descs), idsMat);
      deleteExtractionTable(extractionTable);
      deleteDistModel2D(optModX);
      optModX = NULL;
      deleteDistModel2D(optModY);
      optModY = NULL;
      deleteDistModelFull(crvMod);
      deleteDistModelFull(idsMat);
    }


   /*
    *  Get the Extraction Table from the Master Flat Field.
    */

    if (!(extractionTable = VmSpExTab(masterFlatImage, grismTable,
                                      ifuTable, NULL))) {
      cpl_msg_error(task, "Problems in deriving the extraction table");
      deleteTable(lineCat);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteImage(masterFlatImage);
      deleteImage(combFlatImage);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      cpl_free(lampRon);
      return EXIT_FAILURE;
    }

    if (applyFlat) {
      for (i = 0; i < lampCount; i++)
        if ((tmpImage = VmSpApplyFF(lampList[i], masterFlatImage, 
                                    extractionTable))) {
        deleteImage(lampList[i]);
        lampList[i] = tmpImage;
      }
      else {
        cpl_msg_error(task, "Failure applying flat field to arc lamp frame %d", 
                    i + 1);
        deleteTable(lineCat);
        deleteTable(ccdTable);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        deleteImage(masterFlatImage);
        deleteImage(combFlatImage);
        deleteExtractionTable(extractionTable);
        for (i = 0; i < lampCount; i++)
          deleteImage(lampList[i]);
        cpl_free(lampList);
        cpl_free(lampRon);
        return EXIT_FAILURE;
      }
    }
  }
  else {

    if (extractionTable) {
      readOptDistModel(extractionTable->descs, &optModX, &optModY);
      writeOptDistModelString(&(lampList[0]->descs), optModX, optModY);
      readCurvatureModel(extractionTable->descs, &crvMod);
      writeCurvatureModelString(&(lampList[0]->descs), crvMod);
      readInvDispMatrix(extractionTable->descs, &idsMat);
      writeInvDispMatrixString(&(lampList[0]->descs), idsMat);
      deleteExtractionTable(extractionTable);
      deleteDistModel2D(optModX);
      optModX = NULL;
      deleteDistModel2D(optModY);
      optModY = NULL;
      deleteDistModelFull(crvMod);
      deleteDistModelFull(idsMat);
    }

   /*
    *  Get the Extraction Table from the first arc lamp frame.
    */

    if (!(extractionTable = VmSpExTab(lampList[0], grismTable, 
                                      ifuTable, NULL))) {
      cpl_msg_error(task, "Problems in deriving the extraction table");
      deleteTable(lineCat);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      cpl_free(lampRon);
      return EXIT_FAILURE;
    }
  }
  

 /*
  * Exposure check
  */

  if (validateFrames) {

    cpl_msg_info(task, "Performing exposure check...");


   /* 
    * Read parameters that MUST be the same for all exposures from 
    * the first file
    */

    readIntDescriptor(lampList[0]->descs, pilTrnGetKeyword("SeqPortStartX", 1),
                      &X, NULL);
    readIntDescriptor(lampList[0]->descs, pilTrnGetKeyword("SeqPortPrscX", 1),
                      &prX, NULL);
    readIntDescriptor(lampList[0]->descs, pilTrnGetKeyword("SeqPortStartY", 1),
                      &Y, NULL);
    readIntDescriptor(lampList[0]->descs, pilTrnGetKeyword("SeqPortPrscY", 1),
                      &prY, NULL);

    for (i = 0; i < lampCount; i++) {

     /* 
      * Compute threshold for line detection.
      */

      threshold = NSIGMA * lampRon[i];

     /*  
      * Choose the closest slit to mask shutter center.
      */

      readFloatDescriptor(lampList[i]->descs, 
                          pilTrnGetKeyword("MshuPosL", quadrant), &posL, NULL);
      readFloatDescriptor(lampList[i]->descs, 
                          pilTrnGetKeyword("MshuPosH", quadrant), &posH, NULL);

      mShuMedPos = posL + fabs(posH - posL) / 2;

      slit = extractionTable->slits;

      /* ALEX: skip dead fibers if first, when IFU */
      while (slit->IFUfibTrans < 0.0)
        {
          slit = slit->next;
        }

      distance = fabs(slit->maskY->data[0] - mShuMedPos);
      slit = slit->next;

      while (slit) {

        /* ALEX: skip dead fibers when IFU */
        if (slit->IFUfibTrans >= 0.0)
          {
           newDistance = fabs(slit->maskY->data[0] - mShuMedPos);
           if (distance > newDistance) {
             distance = newDistance;
             goodSlit = slit;
           }
          }

        slit = slit->next;
      }

      windowDimen = (unsigned int)(goodSlit->numRows * 0.5);

      goodSlitX = goodSlit->ccdX->data[goodSlit->numRows / 2] + X + prX - 1;

      goodSlitY = goodSlit->ccdY->data[0] + Y + prY - 1;

      wLen = findColInTab(lineCat, "WLEN");
      flux = findColInTab(lineCat, "FLUX");
      lineFluxRem = 0.0;
      linePosFloat = 0.;

     /*
      * Choose the line with the highest flux.
      */

      for (h = 0; h <= (unsigned int)wLen->len; h++) {
        lineFlux = flux->colValue->fArray[h];
        if (lineFlux > lineFluxRem) {
          lineFluxRem = lineFlux;
          tmpPos = computeDistModel1D(goodSlit->invDis[windowDimen],
                                      wLen->colValue->fArray[h]);

          if ((goodSlitY + (int) tmpPos - windowDimen / 2) > 0 &&
              (goodSlitY + (int) tmpPos + windowDimen / 2) < 4096)
            linePosFloat = tmpPos;
        }
      }

      startX = goodSlitX - windowDimen / 2;
      startY = goodSlitY + (int) linePosFloat - windowDimen / 2;

     /*  
      * Compute lamp level around line position 
      */

      lineShape = (float *)cpl_malloc((windowDimen + 1) * sizeof(float));

      for (j = 0; j <= windowDimen; j++) {
        extrLine = extractFloatImage(lampList[i]->data, lampList[i]->xlen,
                                     lampList[i]->ylen, startX, startY + j, 
                                     windowDimen, 1);
        lineShape[j] = medianPixelvalue(extrLine, windowDimen);
        cpl_free(extrLine);
      }

     /* 
      * Find maximum 
      */

      error = 0;
      if (VM_TRUE == findPeak1D(lineShape, windowDimen, &xmaxFloat, 2)) {

       /* 
        * Compare levels 
        */

        xmax = xmaxFloat;
        if (lineShape[xmax - 1] < threshold || lineShape[xmax] < threshold
            || lineShape[xmax + 1] < threshold) {
          cpl_msg_error(task, "Arc lamp %d is underexposed.", i + 1);
          error = 1;
        }
      } 
      else {
        cpl_msg_error(task, 
                    "Failure in finding check arc line in lamp %d.", i + 1);
        error = 1;
      }

      if (error) {
        deleteTable(lineCat);
        deleteTable(ccdTable);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        deleteExtractionTable(extractionTable);
        deleteImage(masterFlatImage);
        deleteImage(combFlatImage);
        for (i = 0; i < lampCount; i++)
          deleteImage(lampList[i]);
        cpl_free(lampList);
        cpl_free(lampRon);
        cpl_free(lineShape);
        return EXIT_FAILURE;
      }
      cpl_free(lineShape);
    }
  }                              /* End of validation */

  cpl_free(lampRon);


 /*
  * Stack arc lamp images
  */

  if (!(stackImages = VmSpStackFF(lampList,
                                  lampCount, extractionTable, fuzz))) {
    cpl_msg_error(task, "Failure stacking arc lamp exposures");
    deleteTable(lineCat);
    deleteTable(ccdTable);
    deleteTable(grismTable);
    deleteIfuTable(ifuTable);
    deleteExtractionTable(extractionTable);
    deleteImage(masterFlatImage);
    deleteImage(combFlatImage);
    for (i = 0; i < lampCount; i++)
      deleteImage(lampList[i]);
    cpl_free(lampList);
    return EXIT_FAILURE;
  }

  arcLampImage = stackImages[0];
  deleteImage(stackImages[1]);
  cpl_free(stackImages);

 /*
  * Bad pixel cleaning
  */

  if (cleanBadPixel) {
    cpl_msg_info(task, "Cleaning bad pixels from stacked arc lamp frame...");

    if (cleanBadPixels(arcLampImage, ccdTable, 0) == EXIT_FAILURE) {
      cpl_msg_error(task,
                  "Failure cleaning bad pixels from stacked arc lamp frame.");
      deleteTable(lineCat);
      deleteTable(ccdTable);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteExtractionTable(extractionTable);
      deleteImage(masterFlatImage);
      deleteImage(combFlatImage);
      deleteImage(arcLampImage);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      return EXIT_FAILURE;
    }

    deleteTable(ccdTable);

  }

 /*
  * Rebuilding the extraction table, since it was eaten up by VmSpStackFF()
  */

  deleteExtractionTable(extractionTable);

  if (masterFlatImage) {
    if (!(extractionTable = VmSpExTab(masterFlatImage, grismTable, ifuTable,
                                      NULL))) {
      cpl_msg_error(task, "Failure rebuilding the extraction table");
      deleteTable(lineCat);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteImage(arcLampImage); 
      deleteImage(masterFlatImage);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      return EXIT_FAILURE;
    }
  }
  else {
    if (!(extractionTable = VmSpExTab(arcLampImage, grismTable, ifuTable,
                                      NULL))) {
      cpl_msg_error(task, "Failure rebuilding the extraction table");
      deleteTable(lineCat);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteImage(arcLampImage); 
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      return EXIT_FAILURE;
    }
  }

  deleteImage(masterFlatImage);

  if (extrMethod != EXTR_INPUT) {

    /*
     * Computation of the Optical Distortion Model
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

    readIntDescriptor(extractionTable->descs, pilTrnGetKeyword("OptDistOrdX"),
                      &optOrdX, NULL);

    readIntDescriptor(extractionTable->descs, pilTrnGetKeyword("OptDistOrdY"), 
                      &optOrdY, NULL);

    maxOrdOpt = optOrdX > optOrdY ? optOrdX : optOrdY;

    readIntDescriptor(extractionTable->descs, 
                      pilTrnGetKeyword("CurvatureOrdX"), &crvOrdX, NULL);

    readIntDescriptor(extractionTable->descs, 
                      pilTrnGetKeyword("CurvatureOrdY"), &crvOrdY, NULL);

    readIntDescriptor(extractionTable->descs, 
                      pilTrnGetKeyword("DispersionOrdX"), &idsOrdX, NULL);

    readIntDescriptor(extractionTable->descs, 
                      pilTrnGetKeyword("DispersionOrdY"), &idsOrdY, NULL);

    /*
     * Conditions required for trying to improve the default distortion
     * models: slits should be at least two times the number of free
     * parameters, and their coordinates should span for at least 50 mm
     * on the mask.
     */

    fitOpt = 0;
    fitCrv = 0;
    fitIds = 1;   /* Fit always */

    if (maxSlitx - minSlitx > 50.0 && maxSlity - minSlity > 50.0) {
      if (numSlits > 2 * (maxOrdOpt + 1) * (maxOrdOpt + 1))
        fitOpt = 1;
      if (numSlits > (crvOrdX + 1) * (crvOrdY + 1))
        fitCrv = 1;
      if (numSlits > 2 * (idsOrdX + 1) * (idsOrdY + 1))
        fitIds = 1;
    }

/*
    if (lineIdentification == 1)
      fitOpt = 0;
*/

    if (fitOpt) {

      cpl_msg_info(task, 
                 "Computing optical distortion model based on %d slits...", 
                 numSlits);

      coefArrayX = newPixel(numSlits);
      coefArrayY = newPixel(numSlits);

      slit = extractionTable->slits;

      numSlits = 0;
      while (slit) {
        x = slit->ccdX->data[0];
        y = slit->ccdY->data[0];
        slitLength = slit->ccdX->data[slit->numRows - 1] - slit->ccdX->data[0];
        cpl_msg_debug(task, "Slit %d:", numSlits);
        cpl_msg_debug(task, "Old X: %.3f", x);
        cpl_msg_debug(task, "Old Y: %.3f", y);

/* */
    if (alignWavePattern(arcLampImage, x, y, slitLength, &deltaX, &deltaY) 
                                                    == EXIT_SUCCESS) {   
/* */

/* +/
    if (findCentralPosition(arcLampImage, extractionTable->descs,
                    x, y, slitLength, slit->width, lineCat,
                    &deltaX, &deltaY) == 0) {
/+ */
          cpl_msg_debug(task, "New X: %.3f = %.3f + %.3f", x + deltaX, x, deltaX);
          cpl_msg_debug(task, "New Y: %.3f = %.3f + %.3f", y + deltaY, y, deltaY);
          coefArrayX[numSlits].x = slit->maskX->data[0];
          coefArrayX[numSlits].y = slit->maskY->data[0];
          coefArrayX[numSlits].i = x + deltaX;
          coefArrayY[numSlits].x = slit->maskX->data[0];
          coefArrayY[numSlits].y = slit->maskY->data[0];
          coefArrayY[numSlits].i = y + deltaY;
          numSlits++;
        }
        else {
          cpl_msg_debug(task, "Out of CCD!");
        }
        slit = slit->next;
      }
  
      if (numSlits > 2 * (maxOrdOpt + 1) * (maxOrdOpt + 1)) {

        fitDistModel2D(coefArrayX, numSlits, optOrdX, 
                       0.0, 0.0, &optModX, &optRmsX);
        fitDistModel2D(coefArrayY, numSlits, optOrdY, 
                       0.0, 0.0, &optModY, &optRmsY);

        cpl_msg_info(task,
                   "X component of optical distortion model RMS = %.3f pixel", 
                   optRmsX);
        cpl_msg_info(task, 
                   "Y component of optical distortion model RMS = %.3f pixel", 
                   optRmsY);

        /*
         * The fit is accepted only if the RMS of the residuals is better
         * than 1 pixel.
         */

        if (optRmsX < 1.0 && optRmsY < 1.0) {

          writeOptDistModel(&(extractionTable->descs), optModX, optModY);
  
          writeDoubleDescriptor(&(extractionTable->descs), 
                                pilTrnGetKeyword("OptXrms"),
                                optRmsX, pilTrnGetComment("OptXrms"));

          writeDoubleDescriptor(&(extractionTable->descs), 
                                pilTrnGetKeyword("OptYrms"),
                                optRmsY, pilTrnGetComment("OptYrms"));
  
          /*
           * Upgrading the extraction table, according to the new OPT.
           */

          if (!(extractionTable = VmSpExTab(arcLampImage, grismTable, ifuTable,
                                            extractionTable))) {
            cpl_msg_error(task, "Failure upgrading the extraction table with "
                        "new optical distortion model");
            deleteTable(lineCat);
            deleteTable(grismTable);
            deleteIfuTable(ifuTable);
            deleteExtractionTable(extractionTable);
            deleteImage(combFlatImage);
            deleteImage(arcLampImage);
            for (i = 0; i < lampCount; i++)
              deleteImage(lampList[i]);
            cpl_free(lampList);
            return EXIT_FAILURE;
          }
        }
        else {
          cpl_msg_warning(task, 
                        "Bad fit in optical distortion model computation: "
                        "preferring default optical distortion model.");
        }
      }
      else {
        cpl_msg_warning(task, "Not enough slits survived: "
                      "using default optical distortion model.");
      }

      deletePixel(coefArrayX);
      deletePixel(coefArrayY);
      if (optModX != NULL)
        deleteDistModel2D(optModX);
      if (optModY != NULL)
        deleteDistModel2D(optModY);

    }
    else
      cpl_msg_warning(task, "Not enough data for optical distortion model "
                    "computation: using default optical distortion model.");

    if (combFlatImage) {

      if (fitCrv) {

        /*
         * Computation of the Spectral Curvature Model
         */

        stackImages = (VimosImage **) cpl_calloc(2, sizeof(VimosImage *));
        stackImages[0] = combFlatImage;

        if (EXIT_SUCCESS == VmSpDerCurves(stackImages, extractionTable, fuzz)) {
          if (EXIT_FAILURE == VmSpCurveModel(extractionTable, grismTable, 0)) {
            cpl_msg_warning(task, "Failure computing the global spectral "
                          "curvature model from master flat field: using "
                          "the first guess solution");
          }
        }
        else
          cpl_msg_warning(task, "Failure computing the local spectral curvature "
                        "models from master flat field: using the first guess "
                        "solution");

        deleteImage(combFlatImage);
        deleteImage(stackImages[1]);
        cpl_free(stackImages);

        /*
         * Upgrading the extraction table, according to the new CRV.
         */
  
        if (!(extractionTable = VmSpExTab(arcLampImage, grismTable, ifuTable,
                                          extractionTable))) {
          cpl_msg_error(task, "Failure upgrading the extraction table with "
                      "new spectral curvature model");
          deleteTable(lineCat);
          deleteTable(grismTable);
          deleteIfuTable(ifuTable);
          deleteExtractionTable(extractionTable);
          deleteImage(arcLampImage);
          for (i = 0; i < lampCount; i++)
            deleteImage(lampList[i]);
          cpl_free(lampList);
          return EXIT_FAILURE;
        }
      }
      else
        cpl_msg_warning(task, "Not enough data for curvature model computation: "
                      "using default optical distortion model.");
  
    }

    /*
     * Computation of the Inverse Dispersion Solution
     */

#ifdef ONLINE_MODE

      /*FIXME:
       * Always apply refineIDS for HR grisms when realtime processing.
       * Always apply ExractionWindow=10 for LR grisms when realtime processing.
       * THIS IS REALLY DIRTY!
       */

    if (grismName[0] == 'H') {
      cpl_msg_warning(task, "Enforcing RefineIDS=true for online "
                      "processing of HR data!");
      refineIDS = 1;
    }

    if (grismName[0] == 'L') {
      cpl_msg_warning(task, "Enforcing ExtractionWindow=5 for online "
                      "processing of LR data!");
      arcWidth = 5;
    }

#endif

    if (lineIdentification == 0) {

      if (EXIT_FAILURE == VmSpDerDispOld(arcLampImage, extractionTable, 
                                      lineCat, arcWidth, refineIDS)) {
        cpl_msg_error(task, "Failure in deriving wavelength dispersion relation");
        deleteTable(lineCat);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        deleteExtractionTable(extractionTable);
        deleteImage(arcLampImage); 
        for (i = 0; i < lampCount; i++)
          deleteImage(lampList[i]);
        cpl_free(lampList);
        return EXIT_FAILURE;
      }
    }
    else {
      if (EXIT_FAILURE == VmSpDerDisp(arcLampImage, extractionTable,
                                      lineCat, arcWidth, level)) {
        cpl_msg_error(task, "Failure in deriving wavelength dispersion relation");
        deleteTable(lineCat);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        deleteExtractionTable(extractionTable);
        deleteImage(arcLampImage);
        for (i = 0; i < lampCount; i++)
          deleteImage(lampList[i]);
        cpl_free(lampList);
        return EXIT_FAILURE;
      }
    }


    /*
     * Refine optical distortion model if there are enough slits.
     * Currently disabled, it's not really necessary.
     */

    fitOpt = 0;
    if (fitOpt) {
      readFloatDescriptor(extractionTable->descs, pilTrnGetKeyword("WlenCen"), 
                          &lambdaCen, NULL);

      slit = extractionTable->slits;

      coefArrayY = newPixel(numSlits);

      numSlits = 0;
      while(slit) {
        numRows = slit->numRows;
        buffer = cpl_malloc(numRows * sizeof(float));
        for (j = 0, k = 0; k < numRows; k++) {
          if (slit->invDisQuality->data[k] != 0) {
            buffer[j] = computeDistModel1D(slit->invDis[k], lambdaCen);
            ++j;
          }
        }

        if (j > 0) {
          pixelCen = medianPixelvalue(buffer, j);

          coefArrayY[numSlits].x = slit->maskX->data[0];
          coefArrayY[numSlits].y = slit->maskY->data[0];
          coefArrayY[numSlits].i = slit->ccdY->data[0]  /* + pixelCen */  ;

printf("x = %.2f  y = %.2f  i = %.2f + %.2f\n", coefArrayY[numSlits].x, coefArrayY[numSlits].y, coefArrayY[numSlits].i, pixelCen);

          ++numSlits;

        }
        cpl_free(buffer);
        slit = slit->next;
      }

      fitDistModel2D(coefArrayY, numSlits, optOrdY,
                     0.0, 0.0, &optModY, &optRmsY);

      optRmsY = sqrt(optRmsY);

      cpl_msg_info(task, "Y component of optical distortion model "
                 "RMS = %.3f pixel", optRmsY);

      /*
       * The fit is accepted only if the RMS of the residuals is better
       * than 1 pixel.
       */

      if (optRmsY < 1.0) {

        writeOptDistModel(&(extractionTable->descs), NULL, optModY);

        writeDoubleDescriptor(&(extractionTable->descs),
                              pilTrnGetKeyword("OptYrms"),
                              optRmsY, pilTrnGetComment("OptYrms"));

        /*
         * Upgrading the extraction table, according to the new OPT.
         */

        if (!(extractionTable = VmSpExTab(arcLampImage, grismTable, ifuTable,
                                          extractionTable))) {
          cpl_msg_error(task, "Failure upgrading the extraction table with "
                      "new optical distortion model");
          deleteTable(lineCat);
          deleteTable(grismTable);
          deleteIfuTable(ifuTable);
          deleteExtractionTable(extractionTable);
          deleteImage(arcLampImage);
          for (i = 0; i < lampCount; i++)
            deleteImage(lampList[i]);
          cpl_free(lampList);
          return EXIT_FAILURE;
        }

        /*
         * REcomputation of the Inverse Dispersion Solution
         */
    
        if (EXIT_FAILURE == VmSpDerDisp(arcLampImage, extractionTable,
                                        lineCat, arcWidth, level)) {
          cpl_msg_error(task, "Failure in deriving wavelength dispersion");
          deleteTable(lineCat);
          deleteTable(grismTable);
          deleteIfuTable(ifuTable);
          deleteExtractionTable(extractionTable);
          deleteImage(arcLampImage);
          for (i = 0; i < lampCount; i++)
            deleteImage(lampList[i]);
          cpl_free(lampList);
          return EXIT_FAILURE;
        }

      }
      else {
        cpl_msg_warning(task, 
                      "The optical distortion models could not be improved.");
      }
      deletePixel(coefArrayY);
      deleteDistModel2D(optModY);
    }

    /*
     * Compute a new global IDS only if there are enough slits, and
     * only if they span enough portions of the mask (same criterion
     * as for the OPT).
     */

    if (fitIds) {
  
      if (EXIT_FAILURE == VmSpDispMatrix(extractionTable, grismTable, 0)) {
        cpl_msg_error(task, 
                    "Failure in computing Inverse Dispersion Relation Matrix");
        deleteTable(lineCat);
        deleteTable(grismTable);
        deleteIfuTable(ifuTable);
        deleteExtractionTable(extractionTable);
        deleteImage(arcLampImage); 
        for (i = 0; i < lampCount; i++)
          deleteImage(lampList[i]);
        cpl_free(lampList);
        return EXIT_FAILURE;
      }

    }
    else
      cpl_msg_warning(task, "Not enough data for global dispersion relation "
                    "computation: using default model.");

    if (modelSlit)
      modelWavcal(extractionTable, order);

  }

  imabuf = cpl_calloc(2, sizeof(VimosImage *));
  imabuf[0] = arcLampImage;
  imabuf[1] = newImageAndAlloc(arcLampImage->xlen, arcLampImage->ylen);

  if (extrMethod == EXTR_GLOBAL) {

    readOptDistModel(extractionTable->descs, &optModX, &optModY);
    writeOptDistModelString(&(arcLampImage->descs), optModX, optModY);
    readCurvatureModel(extractionTable->descs, &crvMod);
    writeCurvatureModelString(&(arcLampImage->descs), crvMod);
    readInvDispMatrix(extractionTable->descs, &idsMat);
    writeInvDispMatrixString(&(arcLampImage->descs), idsMat);
    deleteDistModel2D(optModX);
    deleteDistModel2D(optModY);
    deleteDistModelFull(crvMod);
    deleteDistModelFull(idsMat);

    if (!(useExtractTable = VmSpExTab(arcLampImage, grismTable,
                                      ifuTable, NULL))) {
      cpl_msg_error(task, "Problems in deriving the global extraction table");
      deleteTable(lineCat);
      deleteTable(grismTable);
      deleteIfuTable(ifuTable);
      deleteImage(arcLampImage);
      deleteExtractionTable(extractionTable);
      for (i = 0; i < lampCount; i++)
        deleteImage(lampList[i]);
      cpl_free(lampList);
      return EXIT_FAILURE;
    }

    deleteExtractionTable(extractionTable);
    extractionTable = useExtractTable;

  }

  extr2D = VmSpEx2D(imabuf, extractionTable, VM_SP_LIN_LAMBDA);
  extracted2D = extr2D[0];

  deleteImage(extr2D[1]); 
  cpl_free(extr2D);

  deleteImage(imabuf[1]);
  cpl_free(imabuf);

/* * 
  createFitsImage("extr2D.fits", extracted2D, "cat");
 * */

  /* %%% */

  if (computeQC) {

    cpl_msg_info(task, "Computing QC1 parameters...");

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      /*
       * Write PRO.CATG, ARCFILE, TPL ID, and OCS CON QUAD, to QC1 group.
       */

      pilQcWriteString("PRO.CATG", extrTableCategory, "Product category");

      qcCopyValue(lampList[0]->descs, pilTrnGetKeyword("ArchiveFile"),
                  NULL, "Archive File Name");

      qcCopyValue(lampList[0]->descs, pilTrnGetKeyword("TplId"),
                  NULL, "Template signature ID");

      qcCopyValue(lampList[0]->descs, pilTrnGetKeyword("Quadrant"),
                  NULL, "Quadrant");

      qcCopyValue(lampList[0]->descs, pilTrnGetKeyword("FilterName", quadrant),
                  NULL, "Filter name");

      qcCopyValue(lampList[0]->descs, pilTrnGetKeyword("GrismName", quadrant),
                  NULL, "Grism name");

      qcCopyValue(lampList[0]->descs, pilTrnGetKeyword("MaskId", quadrant),
                  NULL, "Mask identification");

      qcCopyValue(extractionTable->descs, pilTrnGetKeyword("WlenCen"),
                  NULL, "Reference wavelength");

      if (grismName[0] == 'L') {
        if (grismName[3] == 'r') {      /* LR_red    */
          lambdaHe = 7065.19;
          lambdaNe = 0.0;
          lambdaAr = 7723.80;
        }
        if (grismName[3] == 'b') {      /* LR_blue   */
          lambdaHe = 5015.68;
          lambdaNe = 6598.96;
          lambdaAr = 0.0;
        }
      }

      if (grismName[0] == 'M') {        /* MR        */
          lambdaHe = 7065.19;
          lambdaNe = 7032.41;
          lambdaAr = 7723.80;
      }

      if (grismName[0] == 'H') {
        if (grismName[3] == 'r') {      /* HR_red    */
          lambdaHe = 7065.19;
          lambdaNe = 7032.41;
          lambdaAr = 7723.80;
        }
        if (grismName[3] == 'o') {      /* HR_orange */
          lambdaHe = 7065.19;
          lambdaNe = 7032.41;
          lambdaAr = 7723.80;
        }
        if (grismName[3] == 'b') {      /* HR_blue   */
          lambdaHe = 5015.68;
          lambdaNe = 5944.83;
          lambdaAr = 0.0;
        }
      }

      slit = slitClosestToCenter(extractionTable);

      qcWriteValueDouble(extractionTable->descs, slit->width * focalScale,
                         "QC.MOS.SLIT.WIDTH", "arcsec",
                         "Width of slit closest to center");

      if (lambdaHe > 1.) {
        extractSpecFlux(arcLampImage, slit, lambdaHe, 3, &lflux, &lflux_err);

        lflux     /= time;   /* Assuming all exposure times are the same */
        lflux_err /= time;

        cpl_msg_info(task, "Flux of He %.2f: %.2f +/- %.2f ADU/mm^2/s", 
                   lambdaHe, lflux, lflux_err);

        qcWriteValueDouble(extractionTable->descs, lambdaHe,
                           "QC.MOS.HE.LAMBDA", "Angstrom",
                           "He arc lamp line for flux determination");
  
        qcWriteValueDouble(extractionTable->descs, lflux, 
                           "QC.MOS.HE.FLUX", "ADU/mm^2/s", 
                           "Flux at chosen He wavelength");
  
        qcWriteValueDouble(extractionTable->descs, lflux_err, 
                           "QC.MOS.HE.FLUXERR", "ADU/mm^2/s",
                           "Error on flux at chosen He wavelength");
      }
      else
        cpl_msg_warning(task, "No He lines in %s spectral range: corresponding "
                   "QC1 parameters are not computed.", grismName);

      if (lambdaNe > 1.) {
        extractSpecFlux(arcLampImage, slit, lambdaNe, 3, &lflux, &lflux_err);

        lflux_err /= time;   /* Assuming all exposure times are the same */
        lflux     /= time;

        cpl_msg_info(task, "Flux of Ne %.2f: %.2f +/- %.2f ADU/mm^2/s", 
                   lambdaNe, lflux, lflux_err);

        qcWriteValueDouble(extractionTable->descs, lambdaNe,
                           "QC.MOS.NE.LAMBDA", "Angstrom",
                           "Ne arc lamp line for flux determination");
  
        qcWriteValueDouble(extractionTable->descs, lflux, 
                           "QC.MOS.NE.FLUX", "ADU/mm^2/s", 
                           "Flux at chosen Ne wavelength");
  
        qcWriteValueDouble(extractionTable->descs, lflux_err, 
                           "QC.MOS.NE.FLUXERR", "ADU/mm^2/s",
                           "Error on flux at chosen Ne wavelength");
      }
      else
        cpl_msg_warning(task, "No Ne lines in %s spectral range: corresponding "
                   "QC1 parameters are not computed.", grismName);

      if (lambdaAr > 1.) {
        extractSpecFlux(arcLampImage, slit, lambdaAr, 3, &lflux, &lflux_err);

        lflux_err /= time;   /* Assuming all exposure times are the same */
        lflux     /= time;

        cpl_msg_info(task, "Flux of Ar %.2f: %.2f +/- %.2f ADU/mm^2/s", 
                   lambdaAr, lflux, lflux_err);

        qcWriteValueDouble(extractionTable->descs, lambdaAr,
                           "QC.MOS.AR.LAMBDA", "Angstrom",
                           "Ar arc lamp line for flux determination");
  
        qcWriteValueDouble(extractionTable->descs, lflux, 
                           "QC.MOS.AR.FLUX", "ADU/mm^2/s", 
                           "Flux at reference wavelength");
  
        qcWriteValueDouble(extractionTable->descs, lflux_err, 
                           "QC.MOS.AR.FLUXERR", "ADU/mm^2/s",
                           "Error on flux at reference wavelength");
      }
      else
        cpl_msg_warning(task, "No Ar lines in %s spectral range: corresponding "
                   "QC1 parameters are not computed.", grismName);


      /*
       * QC MOS WAVECAL COEFFi  (continue working with slit closest to
       *                         center obtained above).
       */

      readIntDescriptor(extractionTable->descs, 
                        pilTrnGetKeyword("DispersionOrd"), &idsOrd, NULL);
      
      numRows = slit->numRows;

      buffer = (float *)cpl_malloc(numRows * sizeof(float));

      for (i = 0; i <= (unsigned int)idsOrd; i++) {
        npix = 0;
        for (j = 0; j < numRows; j++) {
          if (slit->invDisQuality->data[j]) {
            buffer[npix] = slit->invDis[j]->coefs[i];
            npix++;
          }
        }

        if (npix == 0) {
          medianValue = 0.0;
          cpl_msg_warning(task, "Slit-closest-to-center was completely lost! "
                        "QC1 parameters QC.MOS.WAVECAL.COEFFi are not "
                        "computed.");
          break;
        }
        else
          medianValue = medianPixelvalue(buffer, npix);

        sprintf(parName, "QC.MOS.WAVECAL.COEFF%d", i);

        switch (i) {
        case 0: 
          sprintf(unit, "pixel"); break;
        case 1: 
          sprintf(unit, "pixel/Angstrom"); break;
        default:
          sprintf(unit, "pixel/Angstrom^%d", i);
        }

        sprintf(comment, "Median coefficient %d of IDS", i);

        qcWriteValueDouble(extractionTable->descs, medianValue,
                           parName, unit, comment);

      }

      cpl_free(buffer);


      /*
       * QC MOS REFWAVE MEAN / RMS
       */

      slit = extractionTable->slits;

      numRows = 0;
      while (slit) {
        numRows += slit->numRows;
        slit = slit->next;
      }

      buffer = (float *)cpl_malloc(numRows * sizeof(float));

      slit = extractionTable->slits;

      npix = 0;
      while (slit) {
        numRows = slit->numRows;
        for (j = 0; j < numRows; j++) {
          if (slit->invDisQuality->data[j]) {
            buffer[npix] = slit->invDis[j]->coefs[0];
            npix++;
          }
        }
        slit = slit->next;
      }

      meanValue = computeAverageFloat(buffer, npix);

      qcWriteValueDouble(extractionTable->descs, meanValue,
                         "QC.MOS.REFWAVE.MEAN", "pixel",
                         "MEAN of CCD positions of reference wavelength");

      rmsValue = 0.0;
      for (i = 0; i < (unsigned int)npix; i++)
        rmsValue += fabs(buffer[i] - meanValue);

      rmsValue /= npix;
      rmsValue *= MEANDEV_TO_SIGMA;

      qcWriteValueDouble(extractionTable->descs, rmsValue,
                         "QC.MOS.REFWAVE.RMS", "pixel",
                         "RMS of CCD positions of reference wavelength");

      cpl_free(buffer);


      /*
       *  QC MOS RESOLUTIONi  and  QC MOS RESOLUTIONi RMS
       */


      if (grismName[0] == 'L') {
        if (grismName[3] == 'r') {      /* LR_red    */
          lambdaRed = 9122.97;
          lambdaYel = 7635.11;
          lambdaBlu = 5875.62;
        }
        if (grismName[3] == 'b') {      /* LR_blue   */
          lambdaRed = 6598.95;
          lambdaYel = 5015.68;
          lambdaBlu = 3888.65;
        }
      }

      if (grismName[0] == 'M') {        /* MR        */
          lambdaRed = 8264.521;
          lambdaYel = 6678.200;
          lambdaBlu = 5015.675;
      }

      if (grismName[0] == 'H') {
        if (grismName[3] == 'r') {      /* HR_red    */
          lambdaRed = 9122.966;
          lambdaYel = 7948.175;
          lambdaBlu = 6929.468;
        }
        if (grismName[3] == 'o') {      /* HR_orange */
          lambdaRed = 7948.175;
          lambdaYel = 6929.468;
          lambdaBlu = 5875.618;
        }
        if (grismName[3] == 'b') {      /* HR_blue   */
          lambdaRed = 6598.953;
          lambdaYel = 5875.618;
          lambdaBlu = 5015.675;
        }
      }

      count = 0;
      if (spectralResolution(extracted2D, lambdaRed, 
                             &meanValue, &rmsValue, 60000) == EXIT_SUCCESS) {
        count++;
        meanFwhm += lambdaRed / meanValue;
      }

      cpl_msg_info(task, "Spectral resolution at %.2f: %.2f +/- %.2f", 
                 lambdaRed, meanValue, rmsValue);

      qcWriteValueDouble(extractionTable->descs, lambdaRed,
                         "QC.MOS.RESOLUTION1.LAMBDA", "Angstrom",
                         "Line used in spectral resolution determination");

      qcWriteValueDouble(extractionTable->descs, meanValue,
                         "QC.MOS.RESOLUTION1", NULL,
                         "Mean spectral resolution at red end of spectrum");

      qcWriteValueDouble(extractionTable->descs, rmsValue,
                         "QC.MOS.RESOLUTION1.RMS", NULL,
                         "RMS of spectral resolution at red end of spectrum");


      if (spectralResolution(extracted2D, lambdaYel,
                             &meanValue, &rmsValue, 60000) == EXIT_SUCCESS) {
        count++;
        meanFwhm += lambdaYel / meanValue;
      }

      cpl_msg_info(task, "Spectral resolution at %.2f: %.2f +/- %.2f", 
                 lambdaYel, meanValue, rmsValue);

      qcWriteValueDouble(extractionTable->descs, lambdaYel,
                         "QC.MOS.RESOLUTION2.LAMBDA", "Angstrom",
                         "Line used in spectral resolution determination");

      qcWriteValueDouble(extractionTable->descs, meanValue,
                         "QC.MOS.RESOLUTION2", NULL,
                         "Mean spectral resolution at center of spectrum");

      qcWriteValueDouble(extractionTable->descs, rmsValue,
                         "QC.MOS.RESOLUTION2.RMS", NULL,
                         "RMS of spectral resolution at center of spectrum");

      if (spectralResolution(extracted2D, lambdaBlu,
                             &meanValue, &rmsValue, 60000) == EXIT_SUCCESS) {
        count++;
        meanFwhm += lambdaBlu / meanValue;
      }

      cpl_msg_info(task, "Spectral resolution at %.2f: %.2f +/- %.2f", 
                 lambdaBlu, meanValue, rmsValue);

      qcWriteValueDouble(extractionTable->descs, lambdaBlu,
                         "QC.MOS.RESOLUTION3.LAMBDA", "Angstrom",
                         "Line used in spectral resolution determination");

      qcWriteValueDouble(extractionTable->descs, meanValue,
                         "QC.MOS.RESOLUTION3", NULL,
                         "Mean spectral resolution at blue end of spectrum");

      qcWriteValueDouble(extractionTable->descs, rmsValue,
                         "QC.MOS.RESOLUTION3.RMS", NULL,
                         "RMS of spectral resolution at blue end of spectrum");

      if (count)
        meanFwhm /= 3.;          /* In Angstrom */
      else {
        readDoubleDescriptor(extracted2D->descs, 
                             pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
        meanFwhm = 4 * cdelt;
      }


      /*
       * Determination of the distortion models RMS. This is computed
       * by a peak detection on all the rows of the 2D extracted image,
       * for each line from the catalog. The window where the peak is
       * searched is twice the FWHM of the lines (derived from the
       * spectral resolution at the center of the spectrum). 
       */

      rmsValue = distortionsRms(extracted2D, lineCat, meanFwhm);

      if (extrMethod == EXTR_GLOBAL) {
        cpl_msg_info(task, "Global IDS RMS: %.3f pixel", rmsValue);
        qcWriteValueDouble(extractionTable->descs, rmsValue,
                           "QC.MOS.IDS.RMS", "pixel", "Global IDS rms");
        writeDoubleDescriptor(&(extractionTable->descs), 
                              pilTrnGetKeyword("IdsYrms"),
                              rmsValue, pilTrnGetComment("IdsYrms"));
      }
      else if (extrMethod == EXTR_INPUT) {
        cpl_msg_info(task, "First guess IDS rms: %.3f pixel", rmsValue);
        qcWriteValueDouble(extractionTable->descs, rmsValue,
                           "QC.MOS.IDS.RMS", "pixel", "Global IDS rms");
        writeDoubleDescriptor(&(extractionTable->descs), 
                              pilTrnGetKeyword("IdsYrms"),
                              rmsValue, pilTrnGetComment("IdsYrms"));
      }
      else {
        cpl_msg_info(task, "Local IDS's mean rms: %.3f pixel", rmsValue);
        qcWriteValueDouble(extractionTable->descs, rmsValue,
                           "QC.MOS.IDS.RMS", "pixel", "Local IDS's mean rms");
      }

      if (pilQcGroupEnd() == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

    }
    else
      cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

  } /* End of QC1 computation. */


  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("MjdObs"), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("Adu2Electron", 1), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("Electron2Adu", 1), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("ReadNoise", 1), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("READ.CLOCK"), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("READ.MODE"), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("READ.SPEED"), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("PROG.ID"), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("OBS.ID"), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("WINi.BINX", 1), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("WINi.BINY", 1), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("SeqWindowSizeX", 1), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("SeqWindowSizeY", 1), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("FilterName", quadrant), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("FilterId", quadrant), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("GrismName", quadrant), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("GrismId", quadrant), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               pilTrnGetKeyword("MaskId", quadrant), NULL);
  vimosDscCopy(&extractionTable->descs, lampList[0]->descs,
               "^ESO INS LAMP[1-5] .*", NULL);


  for (i = 0; i < lampCount; i++)
    deleteImage(lampList[i]);
  cpl_free(lampList);

  deleteImage(arcLampImage);

  deleteTable(grismTable);


  writeStringDescriptor(&(extractionTable->descs), "DATE-OBS", dateObs, "");

  writeStringDescriptor(&(extractionTable->descs), 
                        pilTrnGetKeyword("DoCategory"), extrTableCategory, "");

  writeIntDescriptor(&(extractionTable->descs), pilTrnGetKeyword("Quadrant"), 
                     quadrant, "");

  writeDoubleDescriptor(&(extractionTable->descs), 
                        pilTrnGetKeyword("BeamTemperature", quadrant),
                        temperature, "");


  if (numArcMos)
    writeStringDescriptor(&(extractionTable->descs),
                          pilTrnGetKeyword("InstrumentMode"), "MOS", "");
  else
    writeStringDescriptor(&(extractionTable->descs),
                          pilTrnGetKeyword("InstrumentMode"), "IFU", "");

  if (fitOpt)
    writeIntDescriptor(&(extractionTable->descs), pilTrnGetKeyword("IdsFlag"), 
                       1, "IDS based on fitted distorsion models");
  else
    writeIntDescriptor(&(extractionTable->descs), pilTrnGetKeyword("IdsFlag"), 
                       0, "IDS based on first guess distorsion models");


  cpl_msg_info(task, "Testing arc lamp lines saturation...");

  if (testLineSaturation(extracted2D, lineCat)) {
    cpl_msg_warning(task, "Saturated lines found in arc lamp spectra!");

/* IN CASE ABORT IS PREFERRED, UNCOMMENT THIS:
 *
 *  cpl_msg_error(task, "Saturated lines found in arc lamp spectra!");
 *  deleteExtractionTable(extractionTable);
 *  deleteIfuTable(ifuTable);
 *  deleteImage(extracted2D);
 *  deleteTable(lineCat);
 *  return EXIT_FAILURE;
 */
  }

  deleteTable(lineCat);

 /*
  * Create the product files on disk, set the product attributes and
  * update the set of frames.
  */

 /*
  * Output the primary product: the Extraction Table
  */

  vmstrlower(strcpy(outputExtrName, extrTableCategory));
  /* strcat(outputExtrName, ".TFITS"); */
  strcat(outputExtrName, ".fits");

  error = 1;

  if ((extrFile = newImage(0, 0, NULL))) {
    if (openNewFitsImage(outputExtrName, extrFile) == VM_TRUE) {
      if (writeFitsExtractionTable(extractionTable, extrFile->fptr) 
          == VM_TRUE) {
        error = 0;
        closeFitsImage(extrFile, 0);
        copyToPrimary(outputExtrName, pilTrnGetKeyword("InstrumentMode"));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("Quadrant"));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("MjdObs"));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("READ.CLOCK"));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("PROG.ID"));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("OBS.ID"));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("Adu2Electron", 1));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("Electron2Adu", 1));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("WINi.BINX", 1));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("WINi.BINY", 1));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("SeqWindowSizeX", 1));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("SeqWindowSizeY", 1));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("FilterName", quadrant));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("FilterId", quadrant));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("GrismName", quadrant));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("GrismId", quadrant));
        copyToPrimary(outputExtrName, pilTrnGetKeyword("MaskId", quadrant));
        copyToPrimary(outputExtrName, "ESO INS LAMP1 NAME");
        copyToPrimary(outputExtrName, "ESO INS LAMP2 NAME");
        copyToPrimary(outputExtrName, "ESO INS LAMP3 NAME");
        copyToPrimary(outputExtrName, "ESO INS LAMP4 NAME");
        copyToPrimary(outputExtrName, "ESO INS LAMP5 NAME");
        copyToPrimary(outputExtrName, "ESO INS LAMP1 STATE");
        copyToPrimary(outputExtrName, "ESO INS LAMP2 STATE");
        copyToPrimary(outputExtrName, "ESO INS LAMP3 STATE");
        copyToPrimary(outputExtrName, "ESO INS LAMP4 STATE");
        copyToPrimary(outputExtrName, "ESO INS LAMP5 STATE");

        outputFrame = newPilFrame(outputExtrName, extrTableCategory);

        pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
        pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
        pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
        pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);
        pilSofInsert(sof, outputFrame);
      }
    }
  }

  if (error) {
    cpl_msg_error(task, "Cannot create local product file %s!", outputExtrName);
    deleteExtractionTable(extractionTable);
    deleteIfuTable(ifuTable);
    deleteImage(extracted2D);
    return EXIT_FAILURE;
  }


  createSpectralDistPAF(extractionTable->descs, namePAF);

  if ((pafFileName = createSpectralDistPAF(extractionTable->descs, 
                                           pipeNamePAF)) != NULL) {
    outputFrame = newPilFrame(pafFileName, pafFileCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_PAF);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
    cpl_free(pafFileName);
  }
  else {
    deleteExtractionTable(extractionTable);
    deleteImage(extracted2D);
    cpl_msg_error(task, "Cannot create local product file %s!", namePAF);
    return EXIT_FAILURE;
  }

  deleteExtractionTable(extractionTable);

  vmstrlower(strcpy(extrArcName, arc2DCatg));
  strcat(extrArcName, ".fits");


/* */
  readDoubleDescriptor(extracted2D->descs,
                       pilTrnGetKeyword("Cdelt", 1), &cdelt, 0);
  readDoubleDescriptor(extracted2D->descs,
                       pilTrnGetKeyword("Crval", 1), &crval, 0);
  crval -= cdelt / 2;
  writeDoubleDescriptor(&(extracted2D->descs),
                        pilTrnGetKeyword("Crval", 1), crval, 
                        pilTrnGetComment("Crval"));
/* */

  if (createFitsImage(extrArcName, extracted2D, arc2DCatg)) {
    outputFrame = newPilFrame(extrArcName, arc2DCatg);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
    deleteImage(extracted2D);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", extrArcName);
    deleteImage(extracted2D);
    return EXIT_FAILURE;
  }


 /*
  *  Cleanup
  */

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
                    "vmspcaldisp",
    "Determine spectral distortion models from flat field and arc lamp "
    "exposure.",
    "This recipe is used to determine all the spectral distortions and\n" 
    "transformations from an arc lamp exposure and a master flat field\n" 
    "produced by the recipe vmspflat.\n\n"
    "Input files:\n\n"
    "  DO category:               Type:       Explanation:         Required:\n"
    "  MOS_ARC_SPECTRUM           Raw         Arc lamp exposure       Y\n"
    "  MASTER_BIAS                Calib       Master bias             Y\n"
    "  MASTER_DARK                Calib       Master dark             .\n"
    "  MOS_MASTER_SCREEN_FLAT     Calib       Normalised flat field   .\n"
    "  MOS_COMBINED_SCREEN_FLAT   Calib       Combined flat field     .\n"
    "  GRISM_TABLE                Calib       Grism table             Y\n"
    "  LINE_CATALOG               Calib       Line catalog            Y\n"
    "  EXTRACT_TABLE              Calib       Extraction table        .\n"
    "  CCD_TABLE                  Calib       Bad pixel table         .\n\n"
    "Output files:\n\n"
    "  DO category:               Data type:  Explanation:\n"
    "  EXTRACT_TABLE              FITS table  Extraction table\n"
    "  MOS_ARC_SPECTRUM_EXTRACTED FITS image  Sky subtracted slit spectra\n"
    "  (none)                     PAF         Distortion models\n\n"
    "At least one raw arc lamp exposure should be present in the input SOF.\n"
    "The normalised and the combined flat fields are the products of the\n"
    "recipe vmspflat run on flat field data obtained with the same mask.\n"
    "Neither of them is required for running vmspcaldisp, but if a\n"
    "combined flat field is not given then no spectral curvature model\n" 
    "can be computed, and a first-guess is used in its place. A normalised\n" 
    "master flat field needs to be specified only if a flat field correction\n" 
    "is requested.\n\n"
    "The bad pixel table needs to be specified only if the cleaning of bad\n"
    "pixels is requested.\n\n"
    "An extraction table generated from previous runs might be input to the\n"
    "recipe, in order to iterate the modeling of the spectral distortions.\n\n"
    "The grism table contains necessary information to control the way\n" 
    "spectra are extracted and the determination of the distortion models.\n" 
    "The vmspcaldisp recipe gets from the grism table the wavelength that\n" 
    "should be used as reference (header entry PRO WLEN CEN), and the\n" 
    "spectrum extension in CCD pixels above and below the position of the\n" 
    "reference wavelength (header entries PRO SPECT LLEN LO and PRO SPECT\n" 
    "LLEN HI). Other parameters, used in the construction of the extracted\n"
    "arc lamp slit spectra, are the start and the end wavelength of the\n" 
    "image containing the extracted spectra (header entries PRO WLEN START\n" 
    "and PRO WLEN END), and the step of the sampling along the dispersion\n"
    "direction (header entry PRO WLEN INC).\n\n"
    "The primary recipe product is the extraction table, that contains\n"
    "information about the local modelling of spectral distortions.\n" 
    "A secondary product of the vmspcaldisp recipe is the image of\n"
    "extracted slit spectra, that allows a visual check of the distortion\n"
    "models quality. A last product is the PAF file carrying all the\n" 
    "information related to the spectral distortions. This PAF file is\n" 
    "copied (or moved) to the product directory, and it is identical to\n" 
    "the IWS configuration file MOS_wavecal_G_Q.cmf where Q indicates the\n" 
    "VIMOS quadrant number, and G the grism name) that is created in the\n"
    "directory where vmspcaldisp is launched.",

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

                    vmspcaldisp_create,
                    vmspcaldisp_exec,
                    vmspcaldisp_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
