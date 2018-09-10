/* $Id: vmifustandard.c,v 1.7 2013-08-07 16:45:33 cgarcia Exp $
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
 * $Revision: 1.7 $
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

#include <cpl_memory.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilrecipe.h>
#include <pilfits.h>
#include <pilqc.h>
#include <pilutils.h>

#include "vmimage.h"
#include "vmccdtable.h"
#include "vmmosutils.h"
#include "vmifu.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmextincttable.h"
#include "vmspecphottable.h"
#include "vmstdfluxtable.h"
#include "vmqcutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


static cxint vmifustandard(PilSetOfFrames *);


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
vmifustandard_create(cpl_plugin *plugin)
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
vmifustandard_exec(cpl_plugin *plugin)
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

    if (vmifustandard(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmifustandard");
        
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
vmifustandard_destroy(cpl_plugin *plugin)
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
 *   Extract and calibrate IFU standard star spectrum.
 * 
 * @return EXIT_SUCCESS or EXIT_FAILURE.   
 *
 * @param sof  Set of frames containing one IFU standard star exposure,
 *             a tracing solution, a wavelength calibration, a transmission 
 *             correction, a master bias, and (optionally) a CCD Table. 
 *  
 * @doc
 *   This would subtract the master bias from the standard star exposure, 
 *   and then trace its brightest spectra. Not finding traceable spectra
 *   generates an error. The input tracing solution is aligned to the 
 *   standard star spectral tracings. On the basis of the aligned tracing 
 *   solution, the standard star spectra would be extracted together with 
 *   the sky spectra, resampled at a constant wavelength step, and then 
 *   transmission corrected. The median (sky) spectrum is then subtracted
 *   from all the extracted spectra, and all spectra would be integrated
 *   and used to compose a spatial image of the field-of-view. Finally,
 *   all the spectra corresponding to the position of the standard star
 *   would be summed, and written to a table. A CCD table may also be 
 *   specified, to be used in bad pixel correction.
 *
 * @author C. Izzo, R. Palsa
 */

static cxint 
vmifustandard(PilSetOfFrames *sof)
{

  const char     task[]        = "vmifustandard";

  const char     parameter[]   = "Parameters";

  const char    *idsTag        = pilTrnGetCategory("IfuIds");
  const char    *traceTag      = pilTrnGetCategory("IfuTrace");
  const char    *transTag      = pilTrnGetCategory("IfuTransmission");
  const char    *standardTag   = pilTrnGetCategory("IfuStandard");
  const char    *reducedTag    = pilTrnGetCategory("IfuStandardReduced");
  const char    *fovTag        = pilTrnGetCategory("IfuStdFov");
  const char    *skyTag        = pilTrnGetCategory("IfuScienceSky");
  const char    *starTag       = pilTrnGetCategory("IfuStandardExtracted");
  const char    *sphotTableTag = pilTrnGetCategory("IfuSphotTable");
  const char    *atmTableTag   = pilTrnGetCategory("ExtinctTable");
  const char    *fluxTableTag  = pilTrnGetCategory("StdFluxTable");

  char          *biasMethodTag = NULL;

  char          *idsName;
  char          *traceName;
  char          *transName;
  char           reducedName[PATHNAME_MAX + 1];
  char           fovName[PATHNAME_MAX + 1];
  char           skyName[PATHNAME_MAX + 1];
  char           starName[PATHNAME_MAX + 1];
  char           sphotTableName[PATHNAME_MAX + 1];
  char           parName[30];

  char           grismName[20];
  char           filterName[20];
  char           ifuShutter[10];

  VimosBool      updateOK = VM_TRUE;

  size_t         standardCount, traceCount, idsCount, transCount;
  size_t         detectedCount;

  int            i, j, count, preCount, pos, prePos;
  int            newStart, newEnd;

  unsigned int   cleanBadPixel;
  unsigned int   reduceAnyFrame;
  unsigned int   useSkylines;
  unsigned int   computeQC;
  unsigned int   error;

  int            biasMethodEntry;

  double        *ddata;
  float         *data;
  float         *skydata;
  float         *stardata;

  PilFrame      *ccdFrame, *mBiasFrame;
  PilFrame      *currFrame, *outputFrame;
  PilFrame      *atmFrame, *fluxFrame;

  BiasMethod     biasMethod = BIAS_UNDEF;

  VimosImage    *mBiasImage, *mBiasOver;
  VimosImage    *standardImage = NULL;
  VimosImage    *sciExtracted;
  VimosImage    *skyExtracted;
  VimosImage    *starExtracted;
  VimosImage    *sciFov;
  VimosImage    *tmpImage;
  VimosImage    *atmFile = NULL;
  VimosImage    *fluxFile = NULL;

  VimosImage    *spectrum;
  VimosImage    *extinction;
  VimosImage    *smo_response;
  VimosImage    *response;
  VimosImage    *smo_efficiency;
  VimosImage    *efficiency;
  VimosImage    *stdFlux;
  VimosImage    *convStdFlux;

  VimosTable    *ccdTable = NULL;
  VimosTable    *atmTable = NULL;
  VimosTable    *fluxTable = NULL;
  VimosTable    *sphotTable = NULL;


  int         quadrant;
  int         grism;
  int         xlen, ylen;
  int         extension;
  int         nexp, expno;
  int         row;             /* Reference image row                      */
  int         above;           /* Pixels to trace above reference row      */
  int         below;           /* Pixels to trace below reference row      */
  int         zero;            /* Rough position of zero order contam      */
  float       tolerance = 0.3; /* Max deviation from fit for rejection     */
  int         maxReject;       /* Max rejections to consider a fiber dead  */
  int         shortRange;      /* Radius of short trace                    */
  int         firstSlit, lastSlit, slit;
  int         medianFilterSize = 15;
  int         step = 15;   /* Pixels to skip in median filtering (speedup) */
  int         startPix, endPix;
  int         idsOrder;
  float       gain;
  /* float       ron; */
  double      time, airmass;
  double      shift;
  double      startLambda, endLambda, stepLambda, lambda;
  double      startIntegral, endIntegral;
  double     *integrals;
  double     *p;
  double     *coeff = NULL;
  double      wstart, wstep;
  int         wcount;

  int         responseOrder    = 0;

  cpl_table  *detect = NULL;   /* Positions of detected std star spectra   */
  cpl_table  *short_strace;    /* Short tracing on standard star exposure  */
  cpl_table  *fmodel;          /* Fit of long tracing on flat field        */
  cpl_table  *model;           /* Tracing on flat field aligned to std     */
  cpl_table  *fcoeff;          /* Coefficients of flat field tracings      */
  cpl_table  *short_fcoeff;    /* Coefficients of flat field tracings      */
  cpl_table  *short_scoeff;    /* Coefficients of standard star tracings   */

  cpl_table  *spectra;         /* Extracted standard star spectra          */
  cpl_table  *ids;             /* IDS coefficients for each fiber          */
  cpl_table  *trans;           /* Transmission correction for each fiber   */

  cpl_table  *matches;         /* Flat fibers matching standard exp fibers */
  double      dc0;             /* Mean deviation of matching fibers at row */
  double      dc1;             /* Mean difference of slopes of matching    */

  cpl_image  *standard;        /* Standard star image                      */
  cpl_image  *smo_standard;    /* Smoothed standard star image             */
  cpl_image  *extracted;       /* Image with extracted standard spectra    */
  cpl_image  *cleanExtracted;  /* Image with background spectra to zero    */
  cpl_image  *fov;             /* Reconstructed image of the FOV           */
  cpl_image  *sky;             /* Sky spectrum removed from std spectra    */
  cpl_image  *star;            /* Standard star spectrum                   */


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
  * Check if the wavelength calibration should be refined on skylines.
  */

  useSkylines = pilDfsDbGetBool(parameter, "UseSkylines", 1);


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
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   * Check if only frames with TPL EXPNO = quadrant number should be reduced.
   * This is effective only if the frame is part of a template sequence.
   */

  reduceAnyFrame = pilDfsDbGetBool("Parameters", "ReduceAnyFrame", 0);


  /*
   * Make sure that the input SOF includes the essential datasets.
   */

  if ((atmFrame = pilSofLookup(sof, atmTableTag)))
    pilFrmSetType(atmFrame, PIL_FRAME_TYPE_CALIB);
  else {
    cpl_msg_error(task, "No input atmospheric extinction table found");
    return EXIT_FAILURE;
  }

  if ((fluxFrame = pilSofLookup(sof, fluxTableTag)))
    pilFrmSetType(fluxFrame, PIL_FRAME_TYPE_CALIB);
  else {
    cpl_msg_error(task, "No input standard star flux table found");
    return EXIT_FAILURE;
  }

  standardCount = (int)pilSofFrameCount(sof, standardTag);

  if (standardCount > 1) {
    cpl_msg_error(task, "Too many standard star exposures in input: "
                  "just one is required.");
    return EXIT_FAILURE;
  }

  if (standardCount == 0) {
    cpl_msg_error(task, 
                "At least one standard star exposure is required in input.");
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

  currFrame = pilSofLookup(sof, transTag);
  pilFrmSetType(currFrame, PIL_FRAME_TYPE_CALIB);
  transName = (char *)pilFrmGetName(currFrame);


  cpl_msg_info(task, "Loading input frames...");


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
   * Load the standard star frame
   */

  currFrame = pilSofLookup(sof, standardTag);

  if ((standardImage = openOldFitsFile(pilFrmGetName(currFrame), 1, 0))) {
    pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
    closeFitsImage(standardImage, 0);
  }
  else {
    cpl_msg_error(task, "Failure loading the standard star frame");
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


  if (readDoubleDescriptor(standardImage->descs, 
                           pilTrnGetKeyword("ExposureTime"),
                           &time, NULL) == VM_TRUE) {
    if (time < MIN_DIVISOR) {
      cpl_msg_error(task, "Zero or negative exposure time in input.");
      deleteImage(standardImage);
      deleteImage(mBiasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }
  else {
    cpl_msg_error(task, "Cannot read %s from input exposure.",
                  pilTrnGetKeyword("ExposureTime"));
    deleteImage(standardImage);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }


  if (readDoubleDescriptor(standardImage->descs, pilTrnGetKeyword("Airmass"),
                           &airmass, NULL) == VM_FALSE) {
    if (readDoubleDescriptor(standardImage->descs, "ESO TEL AIRM START",
                           &airmass, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Cannot read %s from input exposure.",
                    pilTrnGetKeyword("Airmass"));
      deleteImage(standardImage);
      deleteImage(mBiasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  if (airmass < 0.0) {
    cpl_msg_error(task, "Negative airmass in input.");
    deleteImage(standardImage);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
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

    if (readIntDescriptor(standardImage->descs, 
                          pilTrnGetKeyword("TplExposures"),
                          &nexp, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Integer descriptor %s not found",
                  pilTrnGetKeyword("TplExposures"));

      deleteImage(standardImage);
      deleteImage(mBiasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }

    if (nexp == 4) {
      if (readIntDescriptor(standardImage->descs,
                            pilTrnGetKeyword("TplExposureNumber"),
                            &expno, NULL) == VM_FALSE) {
        cpl_msg_error(task, "Integer descriptor %s not found",
                    pilTrnGetKeyword("TplExposureNumber"));

        deleteImage(standardImage);
        deleteImage(mBiasImage);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }

      if (readIntDescriptor(standardImage->descs, 
                            pilTrnGetKeyword("Quadrant"),
                            &quadrant, NULL) == VM_FALSE) {
        cpl_msg_error(task, "Integer descriptor %s not found",
                    pilTrnGetKeyword("Quadrant"));

        deleteImage(standardImage);
        deleteImage(mBiasImage);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }

      if (quadrant != expno) {
        cpl_msg_warning(task, "Frame %s is not reduced, because it is not "
                      "expected to contain any standard stars. This recipe "
                      "run will be aborted!", pilFrmGetName(currFrame));
        deleteImage(standardImage);
        deleteImage(mBiasImage);
        deleteCcdTable(ccdTable);
        return EXIT_FAILURE;
      }

    }
  }


  /*
   * Recreate bias overscans using as reference the standard star exposure
   */

  if ((mBiasOver = growOverscans(mBiasImage, standardImage))) {
    if (mBiasImage != mBiasOver) {
      deleteImage(mBiasImage);
      mBiasImage = mBiasOver;
    }
  }
  else {
    cpl_msg_error(task, "Failure in growing overscans in master bias");
    deleteImage(standardImage);
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

  if (VmSubBias(standardImage, mBiasImage, biasMethod) == EXIT_FAILURE) {
    cpl_msg_error(task, "Cannot remove bias from standard star exposure");
    deleteImage(standardImage);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  deleteImage(mBiasImage);

  grism = getGrism(standardImage);
  readIntDescriptor(standardImage->descs, pilTrnGetKeyword("Quadrant"),
                    &quadrant, NULL);
  readStringDescriptor(standardImage->descs,
                    pilTrnGetKeyword("GrismName", quadrant),
                    grismName, NULL);
  readStringDescriptor(standardImage->descs,
                    pilTrnGetKeyword("FilterName", quadrant),
                    filterName, NULL);
  readStringDescriptor(standardImage->descs, pilTrnGetKeyword("IfuMode"),
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
    deleteImage(standardImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  if (cleanBadPixel) {
    cpl_msg_info(task, "Cleaning bad pixels from standard star exposure");
    if (EXIT_FAILURE == cleanBadPixels(standardImage, ccdTable, 0)) {
      cpl_msg_error(task, "Cannot clean bad pixels");
      deleteImage(standardImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  gain = getMeanGainFactor(standardImage);

/*
  cpl_msg_info(task, "Remove cosmic rays hits...");

  ron = computeAverageRon(standardImage);

  VmCosmicClean(standardImage, NULL, 0, 0.0, gain, ron, 4.0, 2.0);
*/

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
   *  The standardImage is cast to a CPL type, because the IFU DRS 
   *  is based on CPL.
   */

  standard = cpl_image_wrap_float(standardImage->xlen, standardImage->ylen, 
                                  standardImage->data);

  cpl_msg_info(task, "Processing data from quadrant %d...", quadrant);

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
     * alignment to standard image tracings. To be entirely consistent,
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
     *  Vertical median filter to standard star exposure.
     */

    smo_standard = cpl_image_vertical_median_filter(standard, medianFilterSize, 
                                                    row, above, below, step);

    /*
     * Detect traceable spectra on standard star exposure, and just 
     * short-trace them.
     */

    cpl_msg_info(task, "Detect traceable spectra on standard star exposure...");

    /* Removed by Peter Weilbacher:
    detect = ifuDetect(smo_standard, row);
    */

    /* 
     * Added by Peter Weilbacher:
     * standard star exposures should contain enough 
     * signal so that a high cutoff works
     */
    detect = ifuDetect(smo_standard, row, 70.0);

    if (detect == NULL) {
      cpl_msg_indent_more();
      cpl_msg_warning(task, "No traceable spectra on standard star exposure: "
                    "extracting along flat field tracings.");
      cpl_msg_indent_less();
      model = fmodel;
    }
    else {
      detectedCount = cpl_table_get_nrow(detect);

      cpl_msg_indent_more();
      if (detectedCount > 1) {
        cpl_msg_info(task, "%zd spectra detected on standard star exposure.", 
                   detectedCount);
        cpl_msg_indent_less();
        cpl_msg_info(task, "Tracing detected spectra...");
      }
      else {
        cpl_msg_info(task, "One spectrum detected on standard star exposure.");
        cpl_msg_indent_less();
        cpl_msg_info(task, "Tracing detected spectrum...");
      }

      short_strace = ifuTraceDetected(smo_standard, row,
                                      shortRange, shortRange, step, detect);

      cpl_table_delete(detect);

      maxReject = 2 * shortRange - 2 * (shortRange / step)
                + 0.50 * (2 * shortRange) / step;

      short_scoeff = ifuFitDetected(short_strace, 1, tolerance, maxReject);

      cpl_table_delete(short_strace);

      cpl_msg_info(task, "Identify traced spectra on standard star exposure...");

      matches = ifuMatch(short_fcoeff, short_scoeff, row, &dc0, &dc1);

      cpl_table_delete(short_scoeff);
      cpl_table_delete(matches);


      /*
       * Alignment of flat tracings to standard star frame.
       */

      cpl_msg_info(task, 
                 "Align flat field tracings to standard star exposure...");

      model = ifuAlign(fcoeff, fmodel, dc0, dc1);

    }

    cpl_table_delete(fcoeff);
    cpl_table_delete(short_fcoeff);
    cpl_image_delete(smo_standard);


    /*
     * Standard star simple spectral extraction
     */

    cpl_msg_info(task, 
                 "Simple extraction of standard star and sky spectra...");

    spectra = ifuSimpleExtraction(standard, model);

    if (detect != NULL)
      cpl_table_delete(model);
    cpl_table_delete(fmodel);


    /*
     *  Load the wavelength calibration
     */

    extension = slit - firstSlit + 1;

    if (grism < 2 && ifuShutter[1] == 'N')
        extension = 1;

    ids = cpl_table_load(idsName, extension, 1);


    if (useSkylines) {

      /*
       *  Align to detected sky lines, if requested
       */

      cpl_msg_info(task, 
                   "Align wavelength calibration to detected sky lines...");

      /* Removed by Peter Weilbacher:
      shift = ifuAlignSkylines(spectra, ids, lambda);
      */

      /* Added by Peter Weilbacher:
       * do not shift individual spectra for standard star exposures 
       */
      shift = ifuAlignSkylines(spectra, ids, lambda, 0);

      cpl_msg_indent_more();
      cpl_msg_info(task, "Alignment: %.2f CCD pixels", shift);
      cpl_msg_indent_less();

    }

    cpl_msg_info(task, "Resample standard star spectra at constant wavelength "
               "step (%.2f Angstrom)", stepLambda);
    ifuResampleSpectra(extracted, spectra, ids, slit, lambda, startLambda,
                       stepLambda);

    cpl_table_delete(ids);
    cpl_table_delete(spectra);

    cpl_msg_indent_less();

  }

  cpl_image_unwrap(standard);


  /*
   *  Apply transmission correction
   */

  cpl_msg_info(task, "Apply transmission correction to extracted spectra...");

  trans = cpl_table_load(transName, 1, 1);
  ifuApplyTransmission(extracted, trans);
  cpl_table_delete(trans);


  /*
   *  Subtract the sky
   */

  sky = ifuSubtractSky(extracted);


  /*
   *  Extract standard star total spectrum.
   */

  cleanExtracted = cpl_image_duplicate(extracted);
  ifuSetZeroLevel(cleanExtracted);
  star = ifuSumSpectrum(cleanExtracted);
  cpl_image_delete(cleanExtracted);


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
   *  To output the CPL images, we convert them into VIMOS format first.
   *  Then we add the descriptor header, and compute the QC1 parameters.
   */

  data = (float *)cpl_image_get_data(fov);
  sciFov = newImage(80, 80, data);
  copyAllDescriptors(standardImage->descs, &(sciFov->descs));

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
  deleteSetOfDescriptors(&(sciFov->descs), "HISTORY*");
  //These keywords shouldn't be present in the input raw, but we nevertheless
  //remove them in case they are present (See PIPE-5422)
  deleteSetOfDescriptors(&(sciFov->descs), "PCOUNT");  
  deleteSetOfDescriptors(&(sciFov->descs), "GCOUNT");  

  data = (float *)cpl_image_get_data(extracted);
  sciExtracted = newImage(xlen, ylen, data);
  copyAllDescriptors(standardImage->descs, &(sciExtracted->descs));

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
  deleteSetOfDescriptors(&(sciExtracted->descs), "HISTORY*");  
  //These keywords shouldn't be present in the input raw, but we nevertheless
  //remove them in case they are present (See PIPE-5422)
  deleteSetOfDescriptors(&(sciExtracted->descs), "PCOUNT");  
  deleteSetOfDescriptors(&(sciExtracted->descs), "GCOUNT");  

  skydata = (float *)cpl_image_get_data(sky);
  skyExtracted = newImage(xlen, 1, skydata);
  copyAllDescriptors(standardImage->descs, &(skyExtracted->descs));

  deleteSetOfDescriptors(&(skyExtracted->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "ESO DPR*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "ESO PRO CCD SKY*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "ESO PRO SKY CCD*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "ESO PRO IDS*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "ESO PRO OPT*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "ESO PRO CRV*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "ESO PRO ZERO*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "CD1_*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "CD2_*");
  deleteSetOfDescriptors(&(skyExtracted->descs), "HISTORY*");
  //These keywords shouldn't be present in the input raw, but we nevertheless
  //remove them in case they are present (See PIPE-5422)
  deleteSetOfDescriptors(&(skyExtracted->descs), "PCOUNT");  
  deleteSetOfDescriptors(&(skyExtracted->descs), "GCOUNT");  

  stardata = (float *)cpl_image_get_data(star);
  starExtracted = newImage(xlen, 1, stardata);
  copyAllDescriptors(standardImage->descs, &(starExtracted->descs));

  deleteSetOfDescriptors(&(starExtracted->descs), "ESO ADA*");
  deleteSetOfDescriptors(&(starExtracted->descs), "ESO DPR*");
  deleteSetOfDescriptors(&(starExtracted->descs), "ESO PRO CCD SKY*");
  deleteSetOfDescriptors(&(starExtracted->descs), "ESO PRO SKY CCD*");
  deleteSetOfDescriptors(&(starExtracted->descs), "ESO PRO IDS*");
  deleteSetOfDescriptors(&(starExtracted->descs), "ESO PRO OPT*");
  deleteSetOfDescriptors(&(starExtracted->descs), "ESO PRO CRV*");
  deleteSetOfDescriptors(&(starExtracted->descs), "ESO PRO ZERO*");
  deleteSetOfDescriptors(&(starExtracted->descs), "CD1_*");
  deleteSetOfDescriptors(&(starExtracted->descs), "CD2_*");
  deleteSetOfDescriptors(&(starExtracted->descs), "HISTORY*");
  //These keywords shouldn't be present in the input raw, but we nevertheless
  //remove them in case they are present (See PIPE-5422)
  deleteSetOfDescriptors(&(starExtracted->descs), "PCOUNT");  
  deleteSetOfDescriptors(&(starExtracted->descs), "GCOUNT");  


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


  writeIntDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Naxis", 1), xlen,
                        pilTrnGetComment("Naxis"));
  writeIntDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Naxis", 2) , 1,
                        pilTrnGetComment("Naxis"));
  writeDoubleDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Crval", 1), startLambda,
                        pilTrnGetComment("Crval"));
  writeDoubleDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Crval", 2), 1.,
                        pilTrnGetComment("Crval"));
  writeDoubleDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Crpix", 1), 1.,
                        pilTrnGetComment("Crpix"));
  writeDoubleDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Crpix", 2), 1.,
                        pilTrnGetComment("Crpix"));
  writeDoubleDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Cdelt", 1), stepLambda,
                        pilTrnGetComment("Cdelt"));
  writeFloatDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Cdelt", 2), 1.,
                        pilTrnGetComment("Cdelt"));
  writeStringDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Ctype", 1), "LAMBDA",
                        pilTrnGetComment("Ctype"));
  writeStringDescriptor(&(skyExtracted->descs),
                        pilTrnGetKeyword("Ctype", 2), " ",
                        pilTrnGetComment("Ctype"));


  writeIntDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Naxis", 1), xlen,
                        pilTrnGetComment("Naxis"));
  writeIntDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Naxis", 2) , 1,
                        pilTrnGetComment("Naxis"));
  writeDoubleDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Crval", 1), startLambda,
                        pilTrnGetComment("Crval"));
  writeDoubleDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Crval", 2), 1.,
                        pilTrnGetComment("Crval"));
  writeDoubleDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Crpix", 1), 1.,
                        pilTrnGetComment("Crpix"));
  writeDoubleDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Crpix", 2), 1.,
                        pilTrnGetComment("Crpix"));
  writeDoubleDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Cdelt", 1), stepLambda,
                        pilTrnGetComment("Cdelt"));
  writeFloatDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Cdelt", 2), 1.,
                        pilTrnGetComment("Cdelt"));
  writeStringDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Ctype", 1), "LAMBDA",
                        pilTrnGetComment("Ctype"));
  writeStringDescriptor(&(starExtracted->descs),
                        pilTrnGetKeyword("Ctype", 2), " ",
                        pilTrnGetComment("Ctype"));


  /*
   * Computation of the efficiency curve and of the response curve
   * begins here.
   */

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
    deleteImage(standardImage);
    sciExtracted->data = NULL;
    deleteImage(sciExtracted);
    cpl_image_delete(extracted);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    skyExtracted->data = NULL;
    deleteImage(skyExtracted);
    cpl_image_delete(sky);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
    return EXIT_FAILURE;
  }


  /*
   * Load the standard star flux table
   */

  error = 1;

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

  if (error) {
    deleteTable(atmTable);
    deleteImage(standardImage);
    sciExtracted->data = NULL;
    deleteImage(sciExtracted);
    cpl_image_delete(extracted);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    skyExtracted->data = NULL;
    deleteImage(skyExtracted);
    cpl_image_delete(sky);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
    return EXIT_FAILURE;
  }


  /*
   * Convert extracted spectrum in electrons per second per Angstrom.
   */

  spectrum = duplicateImage(starExtracted);
  constArithLocal(spectrum, gain / time / stepLambda, VM_OPER_MUL);

  /*
   * Map the atmospheric extinction factors to the same lambda sampling
   * of the extracted spectrum, and convert to actual flux loss.
   */

  extinction = duplicateImage(starExtracted);
  mapTable(extinction, startLambda, stepLambda, atmTable, "WAVE", "EXTINCTION");
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

  stdFlux = duplicateImage(starExtracted);
  mapTable(stdFlux, startLambda, stepLambda, fluxTable, "WAVE", "FLUX");
  deleteTable(fluxTable);

  preCount = count;
  prePos = pos;
  count = 0;
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

  convStdFlux = duplicateImage(starExtracted);

  for (i = 0; i < convStdFlux->xlen; i++) {
    lambda = startLambda + stepLambda * i;
    convStdFlux->data[i] = 0.0026 * lambda * stdFlux->data[i];
  }

  efficiency = duplicateImage(starExtracted);

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
    count = preCount;
    pos = prePos;
/*
    deleteImage(efficiency);
    deleteImage(spectrum);
    deleteImage(stdFlux);
    deleteImage(standardImage);
    sciExtracted->data = NULL;
    deleteImage(sciExtracted);
    cpl_image_delete(extracted);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    skyExtracted->data = NULL;
    deleteImage(skyExtracted);
    cpl_image_delete(sky);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
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

  response = duplicateImage(starExtracted);

  for (i = 0; i < response->xlen; i++) {
    if (efficiency->data[i] > 0.01 && stdFlux->data[i] != 0.0)
      response->data[i] = spectrum->data[i] / stdFlux->data[i];
    else
      response->data[i] = 0.0;
  }

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

  ddata = tblGetDoubleData(sphotTable, "WAVE");
  for (i = 0; i < response->xlen; i++)
    ddata[i] = startLambda + stepLambda * i;

  ddata = tblGetDoubleData(sphotTable, "STD_FLUX");
  for (i = 0; i < response->xlen; i++)
    ddata[i] = stdFlux->data[i];

  deleteImage(stdFlux);

  ddata = tblGetDoubleData(sphotTable, "OBS_FLUX");
  for (i = 0; i < response->xlen; i++)
    ddata[i] = spectrum->data[i];

  deleteImage(spectrum);

  ddata = tblGetDoubleData(sphotTable, "RAW_EFFICIENCY");
  for (i = 0; i < response->xlen; i++)
    ddata[i] = efficiency->data[i];

  deleteImage(efficiency);

  ddata = tblGetDoubleData(sphotTable, "EFFICIENCY");
  for (i = 0; i < response->xlen; i++)
    ddata[i] = smo_efficiency->data[i];

  deleteImage(smo_efficiency);

  ddata = tblGetDoubleData(sphotTable, "RAW_RESPONSE");
  for (i = 0; i < response->xlen; i++)
    ddata[i] = response->data[i];

  deleteImage(response);

  ddata = tblGetDoubleData(sphotTable, "RESPONSE");
  for (i = 0; i < smo_response->xlen; i++)
    ddata[i] = smo_response->data[i];

  deleteImage(smo_response);

  specPhotTableHeader(sphotTable, standardImage->descs);

  vimosDscCopy(&sphotTable->descs, standardImage->descs,
               pilTrnGetKeyword("TargetName"), NULL);

  vimosDscCopy(&sphotTable->descs, standardImage->descs,
               pilTrnGetKeyword("IfuMode"), NULL);

  vimosDscCopy(&sphotTable->descs, standardImage->descs,
               "ESO INS IFUE MAG", NULL);

  if (computeQC) {

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      cpl_msg_info(task, "Computing QC1 parameters...");

      /*
       * Write PRO.CATG, ARCFILE, TPL ID, and OCS CON QUAD, to QC1 group.
       */

      pilQcWriteString("PRO.CATG", sphotTableTag, "Product category");

      qcCopyValue(standardImage->descs, pilTrnGetKeyword("ArchiveFile"),
                  NULL, "Archive File Name");

      qcCopyValue(standardImage->descs, pilTrnGetKeyword("TplId"),
                  NULL, "Template signature ID");

      qcCopyValue(standardImage->descs, pilTrnGetKeyword("Quadrant"),
                  NULL, "Quadrant");

      qcCopyValue(standardImage->descs, 
                  pilTrnGetKeyword("FilterName", quadrant),
                  NULL, "Filter name");

      qcCopyValue(standardImage->descs, pilTrnGetKeyword("GrismName", quadrant),
                  NULL, "Grism name");

      wstart = 3700.;
      wstep = 400.;
      wcount = 15;
      tmpImage = newImageAndAlloc(wcount, 1);
      mapTableDouble(tmpImage, wstart, wstep, sphotTable, "WAVE", "EFFICIENCY");

      for (i = 0; i < wcount; i++) {
        sprintf(parName, "QC.IFU.EFFICIENCY%d.LAMBDA", i + 1);
        qcWriteValueDouble(sphotTable->descs, wstart + wstep * i, parName,
                           "Angstrom", "Wavelength of efficiency evaluation");

        sprintf(parName, "QC.IFU.EFFICIENCY%d", i + 1);
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


  deleteImage(standardImage);


  vmstrlower(strcpy(sphotTableName, sphotTableTag));
  /* strcat(sphotTableName, ".TFITS"); */
  strcat(sphotTableName, ".fits");

  if (writeFitsSpecPhotTable(sphotTableName, sphotTable) != VM_TRUE) {
    cpl_msg_error(task, "Cannot create local product file %s!", sphotTableName);

    deleteTable(sphotTable);
    sciExtracted->data = NULL;
    deleteImage(sciExtracted);
    cpl_image_delete(extracted);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    skyExtracted->data = NULL;
    deleteImage(skyExtracted);
    cpl_image_delete(sky);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
    return EXIT_FAILURE;
  }
  else {
    pilFitsHdrCopy(sphotTableName, 0, NULL, ".*-OBS$", 1);
    pilFitsHdrCopy(sphotTableName, 0, NULL, "^ESO .*", 1);
    for (i = 0; i < 30; i++)
        pilFitsHdrDeleteKeys(sphotTableName, "ESO QC*", 1);

    outputFrame = newPilFrame(sphotTableName, sphotTableTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);
    pilSofInsert(sof, outputFrame);
  }

  deleteTable(sphotTable);

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


  /*
   * Update the sky spectrum header
   */

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(skyExtracted->descs),
                                     pilTrnGetKeyword("DataMin"),
                                     imageMinimum(skyExtracted),
                                     pilTrnGetComment("DataMin"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(skyExtracted->descs),
                                     pilTrnGetKeyword("DataMax"),
                                     imageMaximum(skyExtracted),
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(skyExtracted->descs),
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(skyExtracted),
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(skyExtracted->descs),
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(skyExtracted),
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(skyExtracted->descs),
                                     pilTrnGetKeyword("DataMean"),
                                     imageMean(skyExtracted),
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));


  /*
   * Update the standard star spectrum header
   */

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(starExtracted->descs),
                                     pilTrnGetKeyword("DataMin"),
                                     imageMinimum(starExtracted),
                                     pilTrnGetComment("DataMin"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(starExtracted->descs),
                                     pilTrnGetKeyword("DataMax"),
                                     imageMaximum(starExtracted),
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(starExtracted->descs),
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(starExtracted),
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(starExtracted->descs),
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(starExtracted),
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(starExtracted->descs),
                                     pilTrnGetKeyword("DataMean"),
                                     imageMean(starExtracted),
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    sciExtracted->data = NULL;
    deleteImage(sciExtracted);
    cpl_image_delete(extracted);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    skyExtracted->data = NULL;
    deleteImage(skyExtracted);
    cpl_image_delete(sky);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
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
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", reducedName);
    sciExtracted->data = NULL;
    deleteImage(sciExtracted);
    cpl_image_delete(extracted);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    skyExtracted->data = NULL;
    deleteImage(skyExtracted);
    cpl_image_delete(sky);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
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
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", fovName);
    sciFov->data = NULL;
    deleteImage(sciFov);
    cpl_image_delete(fov);
    skyExtracted->data = NULL;
    deleteImage(skyExtracted);
    cpl_image_delete(sky);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
    return EXIT_FAILURE;
  }

  sciFov->data = NULL;
  deleteImage(sciFov);
  cpl_image_delete(fov);

  vmstrlower(strcpy(skyName, skyTag));
  strcat(skyName, ".fits");

  if (createFitsImage(skyName, skyExtracted, skyTag)) {
    outputFrame = newPilFrame(skyName, skyTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", skyName);
    skyExtracted->data = NULL;
    deleteImage(skyExtracted);
    cpl_image_delete(sky);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
    return EXIT_FAILURE;
  }

  skyExtracted->data = NULL;
  deleteImage(skyExtracted);
  cpl_image_delete(sky);

  vmstrlower(strcpy(starName, starTag));
  strcat(starName, ".fits");

  if (createFitsImage(starName, starExtracted, starTag)) {
    outputFrame = newPilFrame(starName, starTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", starName);
    starExtracted->data = NULL;
    deleteImage(starExtracted);
    cpl_image_delete(star);
    return EXIT_FAILURE;
  }

  starExtracted->data = NULL;
  deleteImage(starExtracted);
  cpl_image_delete(star);

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
                    "vmifustandard",
                    "Extract a standard star spectrum.",
    "This recipe extracts IFU standard star fiber spectra using the input\n"
    "extraction mask, after aligning it to the brightest fiber spectra\n"
    "detected on the input exposure. The extracted spectra are then\n"
    "resampled at a constant wavelength step, after aligning the input\n"
    "wavelength calibration to the positions of a set of identified sky\n"
    "lines. The extracted spectra are finally corrected for the relative\n"
    "differences in transmission from fiber to fiber, and summed together.\n\n"
    "Input files:\n\n"
    "  DO category:            Type:       Explanation:           Required:\n"
    "  IFU_STANDARD            Raw         Standard star spectra     Y\n"
    "  MASTER_BIAS             Calib       Master bias               Y\n"
    "  IFU_TRACE               Calib       Extraction mask           Y\n"
    "  IFU_IDS                 Calib       Wavelength calibration    Y\n"
    "  IFU_TRANSMISSION        Calib       Transmission correction   Y\n"
    "  EXTINCT_TABLE           Calib       Atmospheric extinction    Y\n"
    "  STD_FLUX_TABLE          Calib       Standard star fluxes      Y\n"
    "  CCD_TABLE               Calib       Bad pixel table           .\n\n"
    "Output files:\n\n"
    "  DO category:            Data type:  Explanation:\n"
    "  IFU_STANDARD_REDUCED    FITS image  Reduced standard star spectra\n"
    "  IFU_FOV                 FITS image  Reconstructed field-of-view image\n"
    "  IFU_STANDARD_EXTRACTED  FITS table  Total standard star spectrum\n"
    "  IFU_SCIENCE_SKY         FITS table  Total sky spectrum\n"
    "  IFU_SPECPHOT_TABLE      FITS table  Efficiency and response curves\n\n"
    "The extraction mask, the wavelength calibration, and the relative\n"
    "transmission table, are those generated by the recipe vmifucalib. A CCD\n"
    "table must be specified only if a bad pixel cleaning is requested.\n\n"
    "The spectra extracted from the fibers on the standard star are added\n"
    "and compared to the catalog fluxes of the same star, to obtain the\n"
    "efficiency curve and the response curve to use for flux calibration.\n"
    "The procedure is the same as the one applied in the case of MOS data,\n"
    "by recipe vmmosstandard.\n\n",

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

                    vmifustandard_create,
                    vmifustandard_exec,
                    vmifustandard_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
