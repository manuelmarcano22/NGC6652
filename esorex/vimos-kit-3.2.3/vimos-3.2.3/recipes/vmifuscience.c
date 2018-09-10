/* $Id: vmifuscience.c,v 1.11 2013-08-07 16:42:52 cgarcia Exp $
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
 * $Date: 2013-08-07 16:42:52 $
 * $Revision: 1.11 $
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
#include <cpl_propertylist.h>
#include <cpl_frameset.h>

#include <cpl_memory.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <cpl_imagelist.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilrecipe.h>
#include <pilqc.h>
#include <pilutils.h>

#include "vmimage.h"
#include "vmccdtable.h"
#include "vmmosutils.h"
#include "vmifu.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmqcutils.h"
#include "vmmossphotcalib.h"
#include "vmextincttable.h"
#include "vmspecphottable.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


static cxint vmifuscience(PilSetOfFrames *);


/*
 * Definition of the label strings for all methods the recipe function
 * supports and their associated method code.
 */

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
vmifuscience_create(cpl_plugin *plugin)
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


    p = cpl_parameter_new_value("vimos.Parameters.badpixel.clean",
                                CPL_TYPE_BOOL,
                                "Bad pixel correction on master flat field.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.apply.transmission",
                                CPL_TYPE_BOOL,
                                "Apply transmission correction to extracted "
                                "scientific spectra",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyTransmission");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyTransmission");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.flux.calibration",
                                CPL_TYPE_BOOL,
                                "Extracted spectra are flux calibrated.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CalibrateFlux");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CalibrateFlux");
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


     /* Added by Peter Weilbacher: */

    p = cpl_parameter_new_value("vimos.Parameters.sky.individual",
                                CPL_TYPE_BOOL,
                                "Use sky lines to refine the wavelength "
                                "calibration individually for each spectrum",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "UseSkyIndividual");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "UseSkyIndividual");
    cpl_parameterlist_append(recipe->parameters, p);

    /* End of addition by Peter Weilbacher */


    p = cpl_parameter_new_value("vimos.Parameters.quality.enable",
                                CPL_TYPE_BOOL,
                                "Compute QC1 parameters",
                                "vimos.Parameters",
                                FALSE);
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
vmifuscience_exec(cpl_plugin *plugin)
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

    if (vmifuscience(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmifuscience");
        
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
vmifuscience_destroy(cpl_plugin *plugin)
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
 *   Extract and calibrate IFU science spectra.
 * 
 * @return EXIT_SUCCESS or EXIT_FAILURE.   
 *
 * @param sof  Set of frames containing one IFU scientific exposure,
 *             a tracing solution, a wavelength calibration, a transmission 
 *             correction, a master bias, and (optionally) a CCD Table. 
 *  
 * @doc
 *   This recipe would subtract the master bias from the science exposure, 
 *   and then trace the brightest scientific spectra. In case of successful
 *   tracing of at least one scientific spectrum, the input tracing solution
 *   is aligned to the (partial) scientific tracing; otherwise, the alignment
 *   is not performed. On the basis of such (perhaps aligned) tracing 
 *   solution, the science spectra would be extracted, resampled at a 
 *   constant wavelength step, and then transmission corrected. Finally, 
 *   all the extracted and calibrated spectra would be integrated, and 
 *   the obtained values would be used to compose a spatial image of the 
 *   field-of-view. A CCD table may also be specified, to be used in bad 
 *   pixel correction.
 *
 * @author C. Izzo, R. Palsa
 */

static cxint 
vmifuscience(PilSetOfFrames *sof)
{

  const char     task[]           = "vmifuscience";

  const char     parameter[]      = "Parameters";

  const char  *outIdsTag          = "IFU_SKY_IDS";
  const char  *idsTag             = pilTrnGetCategory("IfuIds");
  const char  *traceTag           = pilTrnGetCategory("IfuTrace");
  const char  *outTraceTag        = "IFU_SKY_TRACE";
  const char  *transTag           = pilTrnGetCategory("IfuTransmission");
  const char  *scienceTag         = pilTrnGetCategory("IfuScience");
  const char  *reducedTag         = pilTrnGetCategory("IfuScienceReduced");
  const char  *fluxReducedTag     = pilTrnGetCategory("IfuScienceFluxReduced");
  const char  *sphotTableCategory = pilTrnGetCategory("IfuSphotTable");
  const char  *atmTableCategory   = pilTrnGetCategory("ExtinctTable");
  const char  *fovTag             = pilTrnGetCategory("IfuFov");

  char        *biasMethodTag      = NULL;

  char        *idsName;
  char         outIdsName[PATHNAME_MAX + 1];
  char        *traceName;
  char         outTraceName[PATHNAME_MAX + 1];
  char        *transName;
  char         reducedName[PATHNAME_MAX + 1];
  char         fovName[PATHNAME_MAX + 1];

  char         grismName[20];
  char         filterName[20];
  char         ifuShutter[10];

  VimosBool    updateOK = VM_TRUE;

  size_t       scienceCount, traceCount, idsCount, transCount;
  size_t       detectedCount;
  float        minSignal, minSigLimit = 10.0; /* Added by Peter Weilbacher */

  int          i, j;
  unsigned int cleanBadPixel;
  unsigned int calibrateFlux;
  unsigned int useSkylines;
  unsigned int useSkyIndividual; /* Added by Peter Weilbacher */
  unsigned int applyTrans;
  unsigned int computeQC;
  unsigned int error;

  int          biasMethodEntry;

  float       *data;

  PilFrame    *ccdFrame, *mBiasFrame;
  PilFrame    *currFrame, *outputFrame;
  PilFrame    *sphotFrame;
  PilFrame    *atmFrame;

  BiasMethod   biasMethod = BIAS_UNDEF;

  VimosImage  *mBiasImage, *mBiasOver;
  VimosImage  *scienceImage = NULL;
  VimosImage  *sphotFile = NULL;
  VimosImage  *atmFile = NULL;
  VimosImage  *sciExtracted = NULL;
  VimosImage  *sciFluxExtracted = NULL;
  VimosImage  *sciFov = NULL;

  VimosTable  *ccdTable = NULL;
  VimosTable  *sphotTable = NULL;
  VimosTable  *atmTable = NULL;

  int          quadrant;
  int          grism;
  int          xlen, ylen;
  int          extension;
  int          row;             /* Reference image row                      */
  int          above;           /* Pixels to trace above reference row      */
  int          below;           /* Pixels to trace below reference row      */
  int          zero;            /* Rough position of zero order contam      */
  float        tolerance = 0.3; /* Max deviation from fit for rejection     */
  int          maxReject;       /* Max rejections to consider a fiber dead  */
  int          shortRange;      /* Radius of short trace                    */
  int          firstSlit, lastSlit, slit;
  int          medianFilterSize = 15;
  int          step = 15;   /* Pixels to skip in median filtering (speedup) */
  int          startPix, endPix;
  int          idsOrder;
  double       shift;
  double       startLambda, endLambda, stepLambda, lambda;
  double       startIntegral, endIntegral;
  double       norm = 1.0;
  double      *integrals;
  double      *p;
  double      *coeff = NULL;

  cpl_propertylist  *header;

  cpl_table   *detect = NULL;   /* Positions of detected science spectra    */
  cpl_table   *short_strace;    /* Short tracing on science                 */
  cpl_table   *fmodel;          /* Fit of long tracing on flat field        */
  cpl_table   *model;           /* Tracing on flat field aligned to science */
  cpl_table   *smodel;          /* Short tracing aligned to science         */
  cpl_table   *fcoeff;          /* Coefficients of flat field tracings      */
  cpl_table   *short_fcoeff;    /* Coefficients of flat field tracings      */
  cpl_table   *short_scoeff;    /* Coefficients of science tracings         */

  cpl_table   *spectra;         /* Extracted science spectra                */
  cpl_table   *ids;             /* IDS coefficients for each fiber          */
  cpl_table   *trans;           /* Transmission correction for each fiber   */

  cpl_table   *matches;         /* Flat fibers matching science fibers      */
  double       dc0;             /* Mean deviation of matching fibers at row */
  double       dc1;             /* Mean difference of slopes of matching    */

  cpl_image   *science;         /* Science image                            */
  cpl_image   *smo_science;     /* Smoothed science image                   */
  cpl_image   *extracted;       /* Image with extracted science spectra     */
  cpl_image   *fov;             /* Reconstructed image of the FOV           */

/***
  cpl_image   *plane;
  cpl_imagelist   *planes;
***/

#ifdef DEBUG_SHIFTS
  char       tablename[1024];
#endif

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
   * Check if the bad pixels should be corrected.
   */

  cleanBadPixel = pilDfsDbGetBool(parameter, "CleanBadPixel", 1);

  /*
   * Check if the extracted scientific spectra should be also corrected
   * for the fiber transmission variations.
   */

  applyTrans = pilDfsDbGetBool("Parameters", "ApplyTransmission", 1);


  /*
   * Check if the spectro-photometric calibration should be applied.
   */

  calibrateFlux = pilDfsDbGetBool(parameter, "CalibrateFlux", 0);


 /*
  * Check if the wavelength calibration should be refined on skylines.
  */

  useSkylines = pilDfsDbGetBool(parameter, "UseSkylines", 1);


  /*
   * Added by Peter Weilbacher, slightly modified by Carlo Izzo
   * to ensure useSkylines is also set.
   */

  useSkyIndividual = pilDfsDbGetBool(parameter, "UseSkyIndividual", 0);
  if (useSkyIndividual) {
    if (useSkylines) {
      cpl_msg_info(task, "Spectra are aligned to sky individually!");
    }
  }

  /* End of addition by Peter Weilbacher */

  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   * Make sure that the input dataset include the essential datasets.
   */

  scienceCount = (int)pilSofFrameCount(sof, scienceTag);

  if (scienceCount > 1) {
    cpl_msg_error(task, "Too many science exposures in input: "
                "just one is required.");
    return EXIT_FAILURE;
  }

  if (scienceCount == 0) {
    cpl_msg_error(task, "At least one science exposure is required in input.");
    return EXIT_FAILURE;
  }

  traceCount = (int)pilSofFrameCount(sof, traceTag);

  if (traceCount > 1) {
    cpl_msg_error(task, "Too many tracing tables in input: "
                "just one is required.");
    return EXIT_FAILURE;
  }

  if (traceCount == 0) {
    cpl_msg_error(task, "At least one tracing table is required in input.");
    return EXIT_FAILURE;
  }

  currFrame = pilSofLookup(sof, traceTag);
  pilFrmSetType(currFrame, PIL_FRAME_TYPE_CALIB);
  traceName = (char *)pilFrmGetName(currFrame);

  idsCount = (int)pilSofFrameCount(sof, idsTag);

  if (idsCount > 1) {
    cpl_msg_error(task, "Too many inverse dispersion solution tables in input: "
                "just one is required.");
    return EXIT_FAILURE;
  }

  if (idsCount == 0) { 
    cpl_msg_error(task, "At least one inverse dispersion solution table is "
                "required in input.");
    return EXIT_FAILURE;
  }

  currFrame = pilSofLookup(sof, idsTag);
  pilFrmSetType(currFrame, PIL_FRAME_TYPE_CALIB);
  idsName = (char *)pilFrmGetName(currFrame);

  currFrame = pilSofLookup(sof, transTag);
  if (currFrame)
    pilFrmSetType(currFrame, PIL_FRAME_TYPE_CALIB);

  if (applyTrans) {
    transCount = (int)pilSofFrameCount(sof, transTag);

    if (transCount > 1) {
      cpl_msg_error(task, "Too many transmission correction tables in input: "
                  "just one is required.");
      return EXIT_FAILURE;
    }
  
    if (transCount == 0) { 
      cpl_msg_error(task, "At least one transmission correction table is "
                  "required in input.");
      return EXIT_FAILURE;
    }

    transName = (char *)pilFrmGetName(currFrame);
  }


  cpl_msg_info(task, "Loading input frames...");


  /*
   * If the flux calibration was enabled, a spectrophotometric table
   * and an atmospheric extinction table are required in the input set 
   * of frames.
   */

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
   * If bad pixel correction or the quality control check are enabled,
   * a bad pixel table is required in the input set of frames. If no 
   * bad pixel table is present this is an error.
   */

  ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable"));
  if (ccdFrame)
    pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

  if (cleanBadPixel) {
    if (ccdFrame) {
      cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));
      if ((ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
        closeFitsTable(ccdTable, 0);
      }
      else {
        cpl_msg_error(task, "Failure in loading the CCD table");
        return EXIT_FAILURE;
      }
    }
    else {
      cpl_msg_error(task, "No CCD table in input: cannot clean bad pixels!");
      return EXIT_FAILURE;
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
      cpl_msg_error(task, "Failure loading the master bias frame");
  }
  else 
    cpl_msg_error(task, "No master bias in input");

  if (error) {
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


  /*
   * Load the science frame
   */

  currFrame = pilSofLookup(sof, scienceTag);

  if ((scienceImage = openOldFitsFile(pilFrmGetName(currFrame), 1, 0))) {
    pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
    closeFitsImage(scienceImage, 0);
  }
  else {
    cpl_msg_error(task, "Failure loading the science frame");
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


  /*
   * Recreate bias overscans using as reference the science exposure
   */

  if ((mBiasOver = growOverscans(mBiasImage, scienceImage))) {
    if (mBiasImage != mBiasOver) {
      deleteImage(mBiasImage);
      mBiasImage = mBiasOver;
    }
  }
  else {
    cpl_msg_error(task, "Failure in growing overscans in master bias");
    deleteImage(scienceImage);
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

  if (VmSubBias(scienceImage, mBiasImage, biasMethod) == EXIT_FAILURE) {
    cpl_msg_error(task, "Cannot remove bias from science exposure");
    deleteImage(scienceImage);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  deleteImage(mBiasImage);

  grism = getGrism(scienceImage);
  readIntDescriptor(scienceImage->descs, pilTrnGetKeyword("Quadrant"),
                    &quadrant, NULL);
  readStringDescriptor(scienceImage->descs,
                    pilTrnGetKeyword("GrismName", quadrant),
                    grismName, NULL);
  readStringDescriptor(scienceImage->descs,
                    pilTrnGetKeyword("FilterName", quadrant),
                    filterName, NULL);
  readStringDescriptor(scienceImage->descs, pilTrnGetKeyword("IfuMode"),
                    ifuShutter, NULL);

  cpl_msg_info(task, "Processing data from quadrant %d, %s/%s, shutter %s...",
             quadrant, grismName, filterName, ifuShutter);

  error = 0;

  if (grism > 1) {                      /* For HR and MR grisms */
    if (ifuShutter[1] == 'F')
      error = 1;                        /* the IFU Shutter should be ON */
  }

  if (grism == 2) {                     /* For MR grisms */
    if (strncmp(filterName, "GG475", 5))
      error = 1;                        /* the filter must be GG475 */
  }

  if (grism == 0) {                     /* For LR_red grisms */
    if (strncmp(filterName, "OS-red", 6))
      error = 1;                        /* the filter must be OS-red */
  }

  if (grism == 1) {                     /* For LR_blue grisms */
    if (strncmp(filterName, "OS-blue", 7))
      error = 1;                        /* the filter must be OS-blue */
  }

  if (error) {
    cpl_msg_error(task, "Unsupported mode!");
    deleteImage(scienceImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  if (cleanBadPixel) {
    cpl_msg_info(task, "Cleaning bad pixels from science exposure");
    if (EXIT_FAILURE == cleanBadPixels(scienceImage, ccdTable, 0)) {
      cpl_msg_error(task, "Cannot clean bad pixels");
      deleteImage(scienceImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  deleteCcdTable(ccdTable);


  /**************
   *  Here begins the real IFU data processing.
   */

  /*
   *  Allocate the image that will contain the extracted spectra.
   */

  ifuRange(grism, &startLambda, &endLambda, &stepLambda); 

  xlen = (endLambda - startLambda) / stepLambda;

  if (grism > 1) {
    firstSlit = 1;
    lastSlit = 1;
    ylen = 400;
  }
  else {
    firstSlit = 0;
    lastSlit = 3;
    ylen = 1600;
  }

  extracted = cpl_image_new(xlen, ylen, CPL_TYPE_FLOAT);


  /* FIXME:
   *  This one is just to get the reference wavelength. It could be
   *  made a bit more elegant than that... An arbitrary quadrant and
   *  slit numbers are specified here, and the first guess IDS is
   *  immediately rejected because not useful to this recipe.
   */

  coeff = ifuFirstIds(grism, 1, 1, &idsOrder, &lambda);
  cpl_free(coeff);

  /*
   *  The scienceImage is cast to a CPL type, because the IFU DRS 
   *  is based on CPL.
   */

  science = cpl_image_wrap_float(scienceImage->xlen, scienceImage->ylen, 
                                 scienceImage->data);

  for (slit = firstSlit; slit <= lastSlit; slit++) {

    if (grism < 2 && ifuShutter[1] == 'N')
      if (slit != 1)
        continue;

    cpl_msg_info(task, "Processing spectra in pseudo-slit %d:", slit + 1);
    cpl_msg_indent_more();


    /*
     *  Extraction parameters. Note that the FG is not needed by this
     *  recipe.
     */

    ifuExtractionParameters(grism, quadrant, slit, 1,
                            &row, &above, &below, &zero);
    cpl_msg_info(task, "Load flat field tracing table...");

    extension = 2 * (slit - firstSlit) + 1;

    if (grism < 2 && ifuShutter[1] == 'N')
        extension = 1;

    fcoeff = cpl_table_load(traceName, extension, 1);
    short_fcoeff = cpl_table_load(traceName, extension + 1, 1);

    fmodel = ifuComputeTraces(fcoeff, row, above, below);


    /* FIXME:
     * This is the short extraction range, used for the flat field
     * alignment to science tracings. To be entirely consistent,
     * this range should be derived from the "short_fcoeff" table,
     * that was created elsewhere (by the recipe vmifucalib).
     * However, currently the short extraction range is fixed,
     * therefore the following lines of code are (at least for now)
     * appropriate.
     */

    shortRange = above < below ? above : below;
    if (shortRange > 400)
      shortRange = 400;


    /*
     *  Vertical median filter to science.
     */

    smo_science = cpl_image_vertical_median_filter(science, medianFilterSize, 
                                                   row, above, below, step);

    /*
     * Detect traceable spectra on science, and just short-trace them.
     */

    cpl_msg_info(task, "Detect traceable spectra on science...");

    /* Old call, before Peter Weilbacher modification:
    detect = ifuDetect(smo_science, row);
    */

    /* Added by Peter Weilbacher, with some small modification
     * by Carlo Izzo to avoid a memory leak: */

    detectedCount = 0;
    minSignal = 70.0; /* Start from 70.0, as in the original code */
    cpl_msg_indent_more();
    while ((detect == NULL || detectedCount < 50) && minSignal >= minSigLimit) {
      cpl_table_delete(detect); detect = NULL;
      detect = ifuDetect(smo_science, row, minSignal);
      if (detect) {
        detectedCount = cpl_table_get_nrow(detect);
      }
      if (detectedCount > 0) {
        cpl_msg_info(task, "Detection threshold: %f, Detected spectra: %zd",
                     minSignal, detectedCount);
      }
      minSignal -= 3.0;
    }
    cpl_msg_indent_less();


    if (detect == NULL) {
      cpl_msg_indent_more();
      cpl_msg_warning(task, "No traceable spectra on scientific exposure: "
                    "extracting along flat field tracings.");
      cpl_msg_indent_less();
      model = fmodel;
    }
    else {

      /* Removed by Peter Weilbacher:
      detectedCount = cpl_table_get_nrow(detect);
      */

      cpl_msg_indent_more();
      if (detectedCount > 1) {
        cpl_msg_info(task, "%zd spectra detected on scientific exposure.", 
                     detectedCount);
        cpl_msg_indent_less();
        cpl_msg_info(task, "Tracing detected spectra...");
      }
      else {
        cpl_msg_info(task, "One spectrum detected on scientific exposure.");
        cpl_msg_indent_less();
        cpl_msg_info(task, "Tracing detected spectrum...");
      }

      short_strace = ifuTraceDetected(smo_science, row,
                                      shortRange, shortRange, step, detect);

      maxReject = 2 * shortRange - 2 * (shortRange / step)
                + 0.50 * (2 * shortRange) / step;

      short_scoeff = ifuFitDetected(short_strace, 1, tolerance, maxReject);

      cpl_table_delete(short_strace);

      cpl_msg_info(task, "Identify traced spectra on science exposure...");

      matches = ifuMatch(short_fcoeff, short_scoeff, row, &dc0, &dc1);

      /* Added by Peter Weilbacher: */

      cpl_msg_info(task, "Flat tracings alignment: offset: %f, slope: %f",
                   dc0, dc1);

      cpl_table_delete(short_scoeff);
      cpl_table_delete(matches);


      /*
       * Alignment of flat tracings to science frame.
       */

      cpl_msg_info(task, "Align flat field tracings to science exposure...");

      model = ifuAlign(fcoeff, fmodel, dc0, dc1);

    }

    if (extension == 1) {
        vmstrlower(strcpy(outTraceName, outTraceTag));
        strcat(outTraceName, ".fits");
        error = cpl_table_save(fcoeff, NULL, NULL, outTraceName, CPL_IO_DEFAULT);
        if (error)
          cpl_msg_error(task, "Cannot create local product file %s!",
                        outTraceName);

    }
    else {
        error = cpl_table_save(fcoeff, NULL, NULL, outTraceName, CPL_IO_EXTEND);
        if (error)
          cpl_msg_error(task, "Cannot extend local product file %s!",
                        outTraceName);
    }

    smodel = ifuAlign(short_fcoeff, fmodel, dc0, dc1);
    cpl_table_delete(smodel);

    error = cpl_table_save(short_fcoeff, NULL, NULL,
                           outTraceName, CPL_IO_EXTEND);
    if (error)
      cpl_msg_error(task, "Cannot extend local product file %s!", outTraceName);

    cpl_table_delete(fcoeff);
    cpl_table_delete(short_fcoeff);
    cpl_image_delete(smo_science);


    /*
     * Science simple spectral extraction
     */

    cpl_msg_info(task, "Extraction of scientific spectra...");

    spectra = ifuExtraction(science, model);

    if (detect != NULL) {
      cpl_table_delete(model);
      cpl_table_delete(detect); detect = NULL;
    }
    cpl_table_delete(fmodel);


    /*
     *  Load the wavelength calibration
     */

    extension = slit - firstSlit + 1;

    if (grism < 2 && ifuShutter[1] == 'N')
        extension = 1;

    ids = cpl_table_load(idsName, extension, 1);

#if DEBUG_SHIFTS
    sprintf(tablename, "ids_1_%s.fits", useSkyIndividual ? "indi" : "norm");
    cpl_table_save(ids, NULL, NULL, tablename, CPL_IO_CREATE);
#endif

    if (useSkylines) {

      /*
       *  Align to detected sky lines, if requested
       */

      cpl_msg_info(task, 
                   "Align wavelength calibration to detected sky lines...");

      /* Removed by Peter Weilbacher:
      shift = ifuAlignSkylines(spectra, ids, lambda);
      */

      /* Added by Peter Weilbacher: */
      shift = ifuAlignSkylines(spectra, ids, lambda, useSkyIndividual);

      cpl_msg_indent_more();
      cpl_msg_info(task, "Alignment: %.2f CCD pixels", shift);
      cpl_msg_indent_less();

#if DEBUG_SHIFTS
      sprintf(tablename, "ids_1x_%s.fits", useSkyIndividual ? "indi" : "norm");
      cpl_table_save(ids, NULL, NULL, tablename, CPL_IO_CREATE);
#endif

      if (extension == 1) {
          vmstrlower(strcpy(outIdsName, outIdsTag));
          strcat(outIdsName, ".fits");
          error = cpl_table_save(ids, NULL, NULL, outIdsName, CPL_IO_DEFAULT);
          if (error)
            cpl_msg_error(task, "Cannot create local product file %s!",
                          outIdsName);

      }
      else {
          error = cpl_table_save(ids, NULL, NULL, outIdsName, CPL_IO_EXTEND);
          if (error)
            cpl_msg_error(task, "Cannot extend local product file %s!",
                          outIdsName);
      }
      outputFrame = newPilFrame(outIdsName, outIdsTag);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
    }

    cpl_msg_info(task, "Resample scientific spectra at constant wavelength "
               "step (%.2f Angstrom)", stepLambda);
    ifuResampleSpectra(extracted, spectra, ids, slit, lambda, startLambda,
                       stepLambda);

#if DEBUG_SHIFTS
    sprintf(tablename, "ids_2_%s.fits", useSkyIndividual ? "indi" : "norm");
    cpl_table_save(ids, NULL, NULL, tablename, CPL_IO_CREATE);
#endif

    cpl_table_delete(ids);
    cpl_table_delete(spectra);

    cpl_msg_indent_less();

  }

  cpl_image_unwrap(science);


  outputFrame = newPilFrame(outTraceName, outTraceTag);
  
  pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
  pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
  pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
  pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);
  
  pilSofInsert(sof, outputFrame);

  if (applyTrans) {

    /*
     *  Apply transmission correction
     */

    cpl_msg_info(task, "Apply transmission correction to extracted spectra...");

    header = cpl_propertylist_load(transName, 0);
    if (cpl_propertylist_has(header, "ESO QC IFU FLAT FLUX"))
      norm = cpl_propertylist_get_double(header, "ESO QC IFU FLAT FLUX");
    cpl_propertylist_delete(header);
    trans = cpl_table_load(transName, 1, 1);
    ifuApplyTransmission(extracted, trans);
    cpl_table_delete(trans);
  }

  /*
   * Get the integration interval for the reconstructed spatial image.
   */

  ifuRangeTransmission(grism, &startIntegral, &endIntegral);
  startPix = (startIntegral - startLambda) / stepLambda;
  endPix = (endIntegral - startLambda) / stepLambda;


  /*
   *  Spatial image reconstruction.
   */

  data = (float *)cpl_image_get_data(extracted);
  integrals = (double *)cpl_calloc(ylen, sizeof(double));

  for (i = 0; i < ylen; i++, data += xlen) {
    for (j = startPix; j < endPix; j++)
      integrals[i] += data[j];
    integrals[i] /= endPix - startPix;
  }

  fov = cpl_image_new(80, 80, CPL_TYPE_FLOAT);

  p = integrals;
  for (slit = firstSlit; slit <= lastSlit; slit++, p += 400)
    ifuImage(fov, p, quadrant, slit);

  cpl_free(integrals);


  /*
   * Cube reconstruction (ELIMINATED)
   */

/***
  data = (float *)cpl_image_get_data(extracted);
  integrals = (double *)cpl_calloc(ylen, sizeof(double));

  planes = cpl_imagelist_new();

  for (j = 0; j < xlen; j++) { \* Plane by plane = wave by wave *\
    data = (float *)cpl_image_get_data(extracted);   \* Reset *\
    for (i = 0; i < ylen; i++, data += xlen) { \* Spec by spec = row by row *\
      integrals[i] = data[j];
    }
    plane = cpl_image_new(80, 80, CPL_TYPE_FLOAT);
    p = integrals;
    for (slit = firstSlit; slit <= lastSlit; slit++, p += 400)
      ifuImage(plane, p, quadrant, slit);

    cpl_imagelist_set(planes, plane, j);
  }

  cpl_free(integrals);

\***
  cube = newCubeAndAlloc(80, 80, xlen);
  cdata = cube->data;

  for (j = 0; j < xlen; j++) {
    plane = cpl_imagelist_get(planes, j);
    data = (float *)cpl_image_get_data(plane);
    for (i = 0; i < 80*80; i++) {
      cdata[j*80*80 + i] = data[i];
    }
  }
***\

  strlower(strcpy(reducedName, cubeTag));
  strcat(reducedName, ".fits");

  header = cpl_propertylist_new();
  cpl_propertylist_append_double(header, "ESO QC IFU FLAT FLUX", norm);

  if (CPL_ERROR_NONE == cpl_imagelist_save(planes, reducedName, 
                                           CPL_BPP_IEEE_FLOAT, header,
                                           CPL_IO_DEFAULT)) {
    outputFrame = newPilFrame(reducedName, cubeTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", reducedName);
    cpl_propertylist_delete(header);
    cpl_image_delete(fov);
    cpl_image_delete(extracted);
    return EXIT_FAILURE;
  }

  cpl_propertylist_delete(header);
  cpl_imagelist_delete(planes);
***/
  

  /*
   *  To output the CPL images, we convert them into VIMOS format first.
   *  Then we add the descriptor header, and compute the QC1 parameters.
   */

  data = (float *)cpl_image_get_data(fov);
  sciFov = newImage(80, 80, data);
  copyAllDescriptors(scienceImage->descs, &(sciFov->descs));

  //cgarcia: I don't know why ESO ADA keywords where deleted, 
  //but they are important to convey the POSANG information
  //deleteSetOfDescriptors(&(sciFov->descs), "ESO ADA*"); 
  deleteSetOfDescriptors(&(sciFov->descs), "ESO DPR*");
  deleteSetOfDescriptors(&(sciFov->descs), "ESO PRO CCD SKY*");
  deleteSetOfDescriptors(&(sciFov->descs), "ESO PRO SKY CCD*");
  deleteSetOfDescriptors(&(sciFov->descs), "ESO PRO IDS*");
  deleteSetOfDescriptors(&(sciFov->descs), "ESO PRO OPT*");
  deleteSetOfDescriptors(&(sciFov->descs), "ESO PRO CRV*");
  deleteSetOfDescriptors(&(sciFov->descs), "ESO PRO ZERO*");
  deleteSetOfDescriptors(&(sciFov->descs), "CD1_*");  
  deleteSetOfDescriptors(&(sciFov->descs), "CD2_*");  
  //These keywords shouldn't be present in the input raw, but we nevertheless
  //remove them in case they are present (See PIPE-5422)
  deleteSetOfDescriptors(&(sciFov->descs), "PCOUNT");  
  deleteSetOfDescriptors(&(sciFov->descs), "GCOUNT");  

  data = (float *)cpl_image_get_data(extracted);
  sciExtracted = newImage(xlen, ylen, data);
  copyAllDescriptors(scienceImage->descs, &(sciExtracted->descs));

  //cgarcia: I don't know why ESO ADA keywords where deleted, 
  //but they are important to convey the POSANG information
  //deleteSetOfDescriptors(&(sciExtracted->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(sciExtracted->descs), "ESO DPR*");
  deleteSetOfDescriptors(&(sciExtracted->descs), "ESO PRO CCD SKY*");
  deleteSetOfDescriptors(&(sciExtracted->descs), "ESO PRO SKY CCD*");
  deleteSetOfDescriptors(&(sciExtracted->descs), "ESO PRO IDS*");
  deleteSetOfDescriptors(&(sciExtracted->descs), "ESO PRO OPT*");
  deleteSetOfDescriptors(&(sciExtracted->descs), "ESO PRO CRV*");
  deleteSetOfDescriptors(&(sciExtracted->descs), "ESO PRO ZERO*");
  deleteSetOfDescriptors(&(sciExtracted->descs), "CD1_*");  
  deleteSetOfDescriptors(&(sciExtracted->descs), "CD2_*");  
  //These keywords shouldn't be present in the input raw, but we nevertheless
  //remove them in case they are present (See PIPE-5422)
  deleteSetOfDescriptors(&(sciExtracted->descs), "PCOUNT");  
  deleteSetOfDescriptors(&(sciExtracted->descs), "GCOUNT");  

  deleteImage(scienceImage);


  writeIntDescriptor(&(sciExtracted->descs),
                        pilTrnGetKeyword("Naxis", 1), xlen, 
                        pilTrnGetComment("Naxis"));
  writeIntDescriptor(&(sciExtracted->descs),
                        pilTrnGetKeyword("Naxis", 2) , ylen, 
                        pilTrnGetComment("Naxis"));
  writeDoubleDescriptor(&(sciExtracted->descs), 
                        pilTrnGetKeyword("Crval", 1), startLambda, 
                        pilTrnGetComment("Crval"));
  writeDoubleDescriptor(&(sciExtracted->descs), 
                        pilTrnGetKeyword("Crval", 2), 1., 
                        pilTrnGetComment("Crval"));
  writeDoubleDescriptor(&(sciExtracted->descs), 
                        pilTrnGetKeyword("Crpix", 1), 1., 
                        pilTrnGetComment("Crpix"));
  writeDoubleDescriptor(&(sciExtracted->descs), 
                        pilTrnGetKeyword("Crpix", 2), 1., 
                        pilTrnGetComment("Crpix"));
  writeDoubleDescriptor(&(sciExtracted->descs), 
                        pilTrnGetKeyword("Cdelt", 1), stepLambda, 
                        pilTrnGetComment("Cdelt"));
  writeFloatDescriptor(&(sciExtracted->descs), 
                        pilTrnGetKeyword("Cdelt", 2), 1., 
                        pilTrnGetComment("Cdelt"));
  writeStringDescriptor(&(sciExtracted->descs), 
                        pilTrnGetKeyword("Ctype", 1), "LAMBDA", 
                        pilTrnGetComment("Ctype"));
  writeStringDescriptor(&(sciExtracted->descs), 
                        pilTrnGetKeyword("Ctype", 2), "FIBRE", 
                        pilTrnGetComment("Ctype"));


  writeIntDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Naxis", 1), 80,
                        pilTrnGetComment("Naxis"));
  writeIntDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Naxis", 2) , 80,
                        pilTrnGetComment("Naxis"));
  writeDoubleDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Crval", 1), 1.,
                        pilTrnGetComment("Crval"));
  writeDoubleDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Crval", 2), 1.,
                        pilTrnGetComment("Crval"));
  writeDoubleDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Crpix", 1), 1.,
                        pilTrnGetComment("Crpix"));
  writeDoubleDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Crpix", 2), 1.,
                        pilTrnGetComment("Crpix"));
  writeDoubleDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Cdelt", 1), 1.,
                        pilTrnGetComment("Cdelt"));
  writeFloatDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Cdelt", 2), 1.,
                        pilTrnGetComment("Cdelt"));
  writeStringDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Ctype", 1), "FIBRE",
                        pilTrnGetComment("Ctype"));
  writeStringDescriptor(&(sciFov->descs),
                        pilTrnGetKeyword("Ctype", 2), "FIBRE",
                        pilTrnGetComment("Ctype"));


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
      sciFov->data = NULL;
      deleteImage(sciFov);
      cpl_image_delete(fov);
      sciExtracted->data = NULL;
      deleteImage(sciExtracted);
      cpl_image_delete(extracted);
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
      sciFov->data = NULL;
      deleteImage(sciFov);
      cpl_image_delete(fov);
      sciExtracted->data = NULL;
      deleteImage(sciExtracted);
      cpl_image_delete(extracted);
      return EXIT_FAILURE;
    }

    sciFluxExtracted = VmSpApplyPhot(sciExtracted, sphotTable, atmTable);

    if (sciFluxExtracted == NULL) {
      cpl_msg_warning(task, "Spectro-photometric calibration failure: "
                      "the corresponding product, %s, will not be created",
                      fluxReducedTag);
      calibrateFlux = 0;
    }

    deleteTable(sphotTable);
    deleteTable(atmTable);

  }

  if (computeQC) {

    cpl_msg_info(task, "QC1 parameters computation - not yet active.");

  } /* End of QC1 computation. */


  /*
   * Update the extracted spectra header
   */

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciExtracted->descs), 
                                     pilTrnGetKeyword("DataMin"),
                                     imageMinimum(sciExtracted), 
                                     pilTrnGetComment("DataMin"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciExtracted->descs), 
                                     pilTrnGetKeyword("DataMax"),
                                     imageMaximum(sciExtracted), 
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciExtracted->descs), 
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(sciExtracted), 
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciExtracted->descs), 
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(sciExtracted), 
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciExtracted->descs), 
                                     pilTrnGetKeyword("DataMean"), 
                                     imageMean(sciExtracted), 
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciExtracted->descs),
                                     "ESO QC IFU FLAT FLUX",
                                     norm,
                                     "Normalisation of transmission correction",
                                     "ESO PRO*", 1));


  /*
   * Update the extracted flux calibrated spectra header
   */

  if (calibrateFlux) {

    updateOK = (updateOK && 
                insertDoubleDescriptor(&(sciFluxExtracted->descs),
                                       pilTrnGetKeyword("ExposureTime"),
                                       1.0,
                                       pilTrnGetComment("ExposureTime"),
                                       "ESO*", 1));

    updateOK = (updateOK &&
                insertDoubleDescriptor(&(sciFluxExtracted->descs),
                                       pilTrnGetKeyword("DataMin"),
                                       imageMinimum(sciFluxExtracted),
                                       pilTrnGetComment("DataMin"),
                                       "ESO*", 1));

    updateOK = (updateOK &&
                insertDoubleDescriptor(&(sciFluxExtracted->descs),
                                       pilTrnGetKeyword("DataMax"),
                                       imageMaximum(sciFluxExtracted),
                                       pilTrnGetComment("DataMax"),
                                       "ESO*", 1));

    updateOK = (updateOK &&
                insertDoubleDescriptor(&(sciFluxExtracted->descs),
                                       pilTrnGetKeyword("DataMedian"),
                                       imageMedian(sciFluxExtracted),
                                       pilTrnGetComment("DataMedian"),
                                       "ESO PRO*", 1));

    updateOK = (updateOK &&
                insertDoubleDescriptor(&(sciFluxExtracted->descs),
                                       pilTrnGetKeyword("DataStdDeviation"),
                                       imageSigma(sciFluxExtracted),
                                       pilTrnGetComment("DataStdDeviation"),
                                       "ESO PRO*", 1));

    updateOK = (updateOK &&
                insertDoubleDescriptor(&(sciFluxExtracted->descs),
                                       pilTrnGetKeyword("DataMean"),
                                       imageMean(sciFluxExtracted),
                                       pilTrnGetComment("DataMean"),
                                       "ESO PRO*", 1));

    updateOK = (updateOK &&
                insertDoubleDescriptor(&(sciFluxExtracted->descs),
                                       "ESO QC IFU FLAT FLUX", norm,
                                "Normalisation of transmission correction",
                                       "ESO PRO*", 1));

    updateOK = (updateOK &&
                insertDoubleDescriptor(&(sciFluxExtracted->descs),
                                       pilTrnGetKeyword("AirMass"), 0.0,
                                       pilTrnGetComment("AirMass"),
                                       "ESO PRO*", 1));

  }


  /*
   * Update the reconstructed FOV header
   */

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciFov->descs),
                                     pilTrnGetKeyword("DataMin"),
                                     imageMinimum(sciFov),
                                     pilTrnGetComment("DataMin"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciFov->descs),
                                     pilTrnGetKeyword("DataMax"),
                                     imageMaximum(sciFov),
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciFov->descs),
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(sciFov),
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciFov->descs),
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(sciFov),
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciFov->descs),
                                     pilTrnGetKeyword("DataMean"),
                                     imageMean(sciFov),
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(sciFov->descs),
                                     "ESO QC IFU FLAT FLUX",
                                     norm,
                                     "Normalisation of transmission correction",
                                     "ESO PRO*", 1));

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    sciExtracted->data = NULL;
    deleteImage(sciExtracted);
    cpl_image_delete(extracted);
    deleteImage(sciFluxExtracted);
    return EXIT_FAILURE;
  }


  /*
   * Create the product file on disk, set the product attributes and
   * update the set of frames.
   */

  vmstrlower(strcpy(reducedName, reducedTag));
  strcat(reducedName, ".fits");

  if (createFitsImage(reducedName, sciExtracted, reducedTag)) {
    outputFrame = newPilFrame(reducedName, reducedTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", reducedName);
    deleteImage(sciFluxExtracted);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    sciExtracted->data = NULL;
    deleteImage(sciExtracted);
    cpl_image_delete(extracted);
    return EXIT_FAILURE;
  }

  sciExtracted->data = NULL;
  deleteImage(sciExtracted);
  cpl_image_delete(extracted);

  vmstrlower(strcpy(fovName, fovTag));
  strcat(fovName, ".fits");

  if (createFitsImage(fovName, sciFov, fovTag)) {
    outputFrame = newPilFrame(fovName, fovTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", fovName);
    deleteImage(sciFluxExtracted);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    return EXIT_FAILURE;
  }

  sciFov->data = NULL;
  deleteImage(sciFov);
  cpl_image_delete(fov);


  if (calibrateFlux) {

    vmstrlower(strcpy(reducedName, fluxReducedTag));
    strcat(reducedName, ".fits");

    if (createFitsImage(reducedName, sciFluxExtracted, fluxReducedTag)) {
      outputFrame = newPilFrame(reducedName, fluxReducedTag);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
    }
    else {
      cpl_msg_error(task, "Cannot create local product file %s!", reducedName);
      deleteImage(sciFluxExtracted);
      return EXIT_FAILURE;
    }

    deleteImage(sciFluxExtracted);

  }

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
                    "vmifuscience",
    "Extract spectra from an IFU scientific exposure.",
    "This recipe extracts IFU scientific spectra using the input extraction\n"
    "mask, after aligning it to the brightest fiber spectra detected on the\n" 
    "input exposure. The extracted spectra are then resampled at a constant\n"
    "wavelength step, after aligning the input wavelength calibration to\n"
    "the positions of a set of identified sky lines. The extracted spectra\n"
    "are eventually corrected for the relative differences in transmission\n"
    "from fiber to fiber.\n\n"
    "Input files:\n\n"
"  DO category:             Type:       Explanation:             Required:\n"
    "  IFU_SCIENCE              Raw         Scientific spectra          Y\n"
    "  MASTER_BIAS              Calib       Master bias                 Y\n"
    "  IFU_TRACE                Calib       Extraction mask             Y\n"
    "  IFU_IDS                  Calib       Wavelength calibration      Y\n"
    "  IFU_TRANSMISSION         Calib       Transmission correction     Y\n"
    "  EXTINCT_TABLE            Calib       Atmospheric extinction      .\n"
    "  IFU_SPECPHOT_TABLE       Calib       Response curve              .\n"
    "  CCD_TABLE                Calib       Bad pixel table             .\n\n"
    "Output files:\n\n"
    "  DO category:             Data type:  Explanation:\n"
    "  IFU_SCIENCE_REDUCED      FITS image  Reduced scientific spectra\n"
    "  IFU_SCIENCE_FLUX_REDUCED FITS image  Reduced scientific spectra\n"
    "  IFU_FOV                  FITS image  Reconstructed IFU field-of-view\n"
    "  IFU_SKY_TRACE            FITS table  Science aligned extraction mask\n\n"
    "  Only if sky alignment is requested:\n"
    "  IFU_SKY_IDS              FITS table  Sky aligned wavelength solution\n"
    "The extraction mask, the wavelength calibration, and the relative\n"
    "transmission table, are those generated by the recipe vmifucalib. A CCD\n"
    "table must be specified only if a bad pixel cleaning is requested.\n\n"
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

                    vmifuscience_create,
                    vmifuscience_exec,
                    vmifuscience_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
