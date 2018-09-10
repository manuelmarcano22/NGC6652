/* $Id: vmdet.c,v 1.4 2012-01-26 16:06:47 cgarcia Exp $
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
 * $Date: 2012-01-26 16:06:47 $
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <cxmemory.h>

#include <cpl.h>

#include <pilmemory.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilrecipe.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>
#include <pilfits.h>
#include <pilqc.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmfit.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmcpl.h"
#include "vmqcutils.h"
#include "vimos_dfs.h"


#define  NSIGMA             (4.)
#define  SAMPLE_REGION      (200)
#define  MAXITER            (4)
#define  CELLSIZE           (128)
#define  DEGREE             (2)    /* Degree of polynomial used to fit
                                      a relation between slopes variance
                                      and flux */
#define  MAX_COMMENT_LENGTH (80)

/*
 * Definition of the label strings for all methods the recipe function
 * supports, with their associated method code.
 */

static const char *method_names[] = { "Tolerant", "Intolerant" };
static const int   methods[]      = { 1,           0 };

static unsigned int n_methods = sizeof(methods) / sizeof(int);


static cxint vmdet(PilSetOfFrames *);


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmdet_create(cpl_plugin *plugin)
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

    p = cpl_parameter_new_value("vimos.Parameters.map",
                                CPL_TYPE_BOOL,
                                "Bad pixel map creation",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CreateBadPixelMap");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CreateBadPixelMap");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.error",
                                CPL_TYPE_BOOL,
                                "Response slope uncertainty map creation",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CreateErrorImage");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CreateErrorImage");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.threshold",
                                CPL_TYPE_DOUBLE,
                                "Threshold in determination of bad pixels",
                                "vimos.Parameters",
                                5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "DetectionThreshold");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "DetectionThreshold");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.method",
                               CPL_TYPE_STRING,
                               "Bad pixel detection method",
                               "vimos.Parameters",
                               "Intolerant", 2,
                               "Intolerant", "Tolerant");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "DetectionMode");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "DetectionMode");
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
vmdet_exec(cpl_plugin *plugin)
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

    if (vmdet(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmdet");
        
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
vmdet_destroy(cpl_plugin *plugin)
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
 *   Determine RON, Gain, and create Bad Pixel Map and CCD table.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE.
 *
 * @param sof            Input set of frames. It should contain at
 *                       least five pairs of flat field exposures,
 *                       each pair corresponding to a different 
 *                       exposure time. All frames should belong 
 *                       to the same quadrant, either in imaging 
 *                       or in MOS mode.
 *
 * @doc
 *   A pixel is "bad" when the slope of the linear fit of each
 *   image exposure level versus pixel value deviates from the 
 *   average of all the slopes found for each pixel by more than 
 *   a given threshold. Bad pixels are written to the CCD table 
 *   columns and, if requested, as pixels of value 1.0 in a 0-filled 
 *   image having the same size of the chip (overscans are removed).
 *   For debug purposes, an error image containing the uncertainties
 *   on the fitted slopes can be created.
 *   
 *   The Read Out Noise and the Gain determinations are based on
 *   the paper by L. Mortara and A. Fowler (1981), Solid State 
 *   Imagers for Astronomy, in SPIE Vol. 290, p. 28. In principle, 
 *   RON and Gain can be obtained by least squaring the relation
 *   
 *        variance(S) = (RON/Gain)^2 + average(S)/Gain
 *   
 *   where S is the measured unbiased signal (in ADU), Gain is in 
 *   electron/ADU (that is, ESO DET OUTi CONAD, _not_ ESO DET OUTi GAIN), 
 *   and RON is in electrons. With this method the gain determination 
 *   is accurate, but the RON determination is very poor. In practice, 
 *   the RON is determined here by directly measuring the variance 
 *   within the overscan regions.
 *   
 *   The bias is roughly removed from each image by subtracting the 
 *   overscan regions mean value.
 *   
 *   Since the relation between variance and signal is linear, 
 *   we can use the above relation using the average signal and 
 *   variance determined on just a portion of the chip. The variance 
 *   is evaluated from the difference of frames with same exposure
 *   time. The found bad pixels are excluded from the determination 
 *   of RON and gain.
 *
 * @author C. Izzo
 */

static cxint
vmdet(PilSetOfFrames *sof)
{

  char task[] = "vmdet";

  struct rlimit rlim;

 /*
  *  Data structures
  */

  PilFrame    *frame;
  PilFrame    *outputFrame;
  VimosImage  *firstImage    = NULL;
  VimosImage  *biasImage     = NULL;
  VimosImage  *biasOver      = NULL;
  VimosImage  *badPixelImage = NULL;
  VimosImage  *errorImage    = NULL;
  VimosImage  *meanImage     = NULL;
  VimosImage  *filtered      = NULL;
  VimosImage **flatList;
  VimosImage **flatList1;
  VimosImage **flatList2;
  VimosImage **pairMean;
  VimosImage **pair;
  VimosTable  *ccdTable;
  VimosDpoint *flux_variance[4];        /* Four tables flux versus variance */
  VimosDpoint *time_flux[4];            /* Four tables time versus flux     */
  VimosDpoint *exposure_flux = NULL;    /* Table exposure versus flux       */
  VimosDpoint *flux_sigma    = NULL;    /* Flux versus slopes sigma         */
  VimosDpoint *flux_average  = NULL;    /* Flux versus slopes average       */

 /*
  *  Some auxiliary variables
  */

  char         badPixelImageName[PATHNAME_MAX + 1];
  char         errorImageName[PATHNAME_MAX + 1];
  char         ccdTableName[PATHNAME_MAX + 1];
  const char  *masterBiasCategory = pilTrnGetCategory("MasterBias");
  const char  *flatCategory  = pilTrnGetCategory("DetectorProperties");
  const char  *ccdTableCategory   = pilTrnGetCategory("CcdTable");
  const char  *badPixelCategory   = pilTrnGetCategory("BadPixelMap");
  int          methodEntry;
  int          minFrames = 10;
  int          computeQC;
  int          tolerant = 0;
  int          flatCount, pairCount, newPairCount;
  int          maxFlat = 0;
  int          i, j, k, count, nIter;
  int          x, y, xlen, ylen, npix;
  int          nBoxPix, nSidesX, nSidesY, boxPosX, boxPosY;
  int          winStartX, winStartY, winSizeX, winSizeY;
  int          winPrscStartX, winPrscStartY;
  int          excluded;
  int          quadrant;
  int          error;
  int          keep_flat;
  int         *keep;
  int         *expInt;
  int          maxInt, minInt, maxPos, minPos;
  float        expTime;
  float       *sample1;
  float       *sample2;
  float       *diff;
  float       *cell;
  float        gain[4];                    
  float        delay[4];                    
  float        delay1, delay2;
  float       *noise;                   
  float        value, mean_delay, err_delay, mean_gain, err_gain;
  float        mean_noise, err_noise;
  double       rms;
  double      *fit;
  double       max = 0.0;
  double       ratio = 0.0;
  double       average, newAverage, delta;
  double       tolerance, tmp;
  double      *sigmaCurve   = NULL;
  double      *averageCurve = NULL;
  char        *methodTag = NULL;
  VimosBool    updateOK = VM_TRUE;
  double       badPixelThresh;
  unsigned int createError = 1;
  unsigned int createMap = 1;


  /*
   *  Get task parameters from the recipe database
   */

  methodTag = (char *)pilDfsDbGetString("Parameters", "DetectionMode");

  if ((methodEntry = strselect(methodTag, method_names, n_methods)) < 0) {
    cpl_msg_error(task, "%s: Invalid bad pixel detection mode.", methodTag);
    return EXIT_FAILURE;
  }

  tolerant = methods[methodEntry];


  /* 
   *  Get the threshold used in the determination of bad pixels.
   */
  
  badPixelThresh = pilDfsDbGetDouble("Parameters", "DetectionThreshold", 5.0);

  if (badPixelThresh < 1.0) {
    cpl_msg_error(task,
             "Invalid bad pixel threshold, it should be at least one sigma.");
    return EXIT_FAILURE;
  }


  /* 
   *  Check if a bad pixel image should be created.
   */
  
  createMap = pilDfsDbGetBool("Parameters", "CreateBadPixelMap", 1);


  /*
   *  Check if an image of the uncertainty of the slope fitting pixel 
   *  intensities versus exposures should be created.
   */

  createError = pilDfsDbGetBool("Parameters", "CreateErrorImage", 0);


  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   *  Make sure that there are enough raw flat fields in the input set.
   */

  flatCount = (int) pilSofFrameCount(sof, flatCategory);

  if (flatCount < minFrames) {
    if (flatCount == 0)
      cpl_msg_error(task, "Missing input flat field frames.");
    else
      cpl_msg_error(task, "Input flat fields should be at least %d", minFrames);
    return EXIT_FAILURE;
  } 

  if (flatCount % 2) {
    cpl_msg_error(task, "Incomplete set of input flat field frames.");
    return EXIT_FAILURE;
  }


  /* FIXME:
   * There is a limit to the number of FITS files that can be
   * opened simultaneously. This must be fixed as soon as possible...
   */

  if (flatCount > 24) {
    cpl_msg_warning(task, "Limit of max 25 opened FITS files reached. "
                  "Processing only first 24 flat fields.");
    flatCount = 24;
  }


  /*
   *  Open all images, without mapping the data into memory,
   *  detect pairs with equal exposure times, and store them
   *  in twin lists.
   */

  flatList  = (VimosImage **)cpl_calloc(flatCount, sizeof(VimosImage *));

  frame = pilSofLookupNext(sof, flatCategory);

  for (i = 0; i < flatCount; i++) {
    flatList[i] = openOldFitsFile(pilFrmGetName(frame), 0, 0);
    if ((loadFitsHeader(flatList[i])) == VM_FALSE) {
      cpl_msg_error(task, "Cannot open flat field %s", pilFrmGetName(frame));
      for (j = 0; j < i; j++)
        deleteImage(flatList[j]);
      cpl_free(flatList);
      return EXIT_FAILURE;
    }
    pilFrmSetType(frame, PIL_FRAME_TYPE_RAW);
    frame = pilSofLookupNext(sof, NULL);
  }

  pairCount = flatCount / 2;

  flatList1 = (VimosImage **)cpl_calloc(pairCount, sizeof(VimosImage *));
  flatList2 = (VimosImage **)cpl_calloc(pairCount, sizeof(VimosImage *));
  expInt    = (int *)cpl_calloc(pairCount, sizeof(int));

  for (i = 0; i < flatCount; i++) {
    if (readFloatDescriptor(flatList[i]->descs, 
                            pilTrnGetKeyword("ExposureTime"),
                            &expTime, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Cannot read descriptor %s", 
                  pilTrnGetKeyword("ExposureTime"));
      for (i = 0; i < flatCount; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);
      cpl_free(flatList1);
      cpl_free(flatList2);
      cpl_free(expInt);
      return EXIT_FAILURE;
    }

    for (j = 0; j < pairCount; j++) {
      if (expInt[j] == 0) {                              /* Free slot     */
        expInt[j] = floor(10 * expTime + 0.5);
        flatList1[j] = flatList[i];
        break;
      }
      if (expInt[j] == (int)floor(10 * expTime + 0.5)) { /* Same exposure */
        if (flatList2[j]) {
          cpl_msg_error(task, "Too many frames with same exposure time (%f)",
                      expTime);
          for (i = 0; i < flatCount; i++)
            deleteImage(flatList[i]);
          cpl_free(flatList);
          cpl_free(flatList1);
          cpl_free(flatList2);
          cpl_free(expInt);
          return EXIT_FAILURE;
        }
        flatList2[j] = flatList[i];
        break;
      }
    }

    if (j == pairCount) {            /* A free slot was not found */
      cpl_msg_error(task, "Incomplete pairs of frames with same exposure time. "
                  "Check the exposure times of the input flat fields.");
      for (i = 0; i < flatCount; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);
      cpl_free(flatList1);
      cpl_free(flatList2);
      cpl_free(expInt);
      return EXIT_FAILURE;
    }
    
  }


  /*
   * Get the master bias frame.
   */

  error = 1;

  if ((frame = pilSofLookup(sof, masterBiasCategory))) {
    pilFrmSetType(frame, PIL_FRAME_TYPE_CALIB);
    if ((biasImage = openOldFitsFile(pilFrmGetName(frame), 1, 0))) {
      closeFitsImage(biasImage, 0);
      error = 0;
    }
    else
      cpl_msg_error(task, "Failure opening master bias frame");
  }
  else
    cpl_msg_error(task, "No master bias found in input");

  if (error) {
    for (i = 0; i < flatCount; i++)
      deleteImage(flatList[i]);
    cpl_free(flatList);
    cpl_free(flatList1);
    cpl_free(flatList2);
    cpl_free(expInt);
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
    cpl_free(flatList1);
    cpl_free(flatList2);
    cpl_free(expInt);
    deleteImage(biasImage);
    return EXIT_FAILURE;
  }

  cpl_free(flatList);


  /*
   *  Check available memory
   */

  if (getrlimit(RLIMIT_DATA, &rlim)) {
    cpl_msg_warning(task, "getrlimit() failure: "
                  "set max limit of input flats to infinity.");
  }
  else {
    rlim.rlim_cur = rlim.rlim_max;
    setrlimit(RLIMIT_DATA, &rlim);
    if (rlim.rlim_cur == RLIM_INFINITY)
      keep_flat = 9999;
    else {
      keep_flat = rlim.rlim_cur 
                / (biasImage->xlen * biasImage->ylen * sizeof(float))
                - 6;  /* Extra space for safety, measured in "frames" */
    }
  }

  cpl_msg_debug(task, "Memory limit: %d flat pairs, Input: %d", 
              keep_flat, pairCount);

  /*
   * This part is just to flag the flat fields that will be used
   * for the bad pixel determination. Not more than keep_flat will
   * be used, for reasons of limited memory available.
   */

  keep = (int *)cpl_calloc(pairCount, sizeof(int));

  if (pairCount <= keep_flat) {
    for (i = 0; i < pairCount; i++)
      keep[i] = 1;
  }
  else {
    i = 0;
    while (i < keep_flat) {
      maxInt = 0;
      for (j = 0; j < pairCount; j++) {
        if (keep[j] == 0) {
          if (expInt[j] > maxInt) {
            maxInt = expInt[j];
            maxPos = j;
          }
        }
      }
      keep[maxPos] = 1;
      i++;
      if (i == keep_flat)
        break;
      minInt = maxInt;
      for (j = 0; j < pairCount; j++) {
        if (keep[j] == 0) {
          if (expInt[j] < minInt) {
            minInt = expInt[j];
            minPos = j;
          }
        }
      }
      keep[minPos] = 1;
      i++;
    }
  }

/*
for (i = 0; i < pairCount; i++)
printf("keep[%d] = %d (%d)\n", i, keep[i], expInt[i]);
*/

  /* End of patch for memory limits. To remove it, delete also
   * any reference to keep[]. */

  cpl_free(expInt);


 /*
  *  Compute the mean frame of each pair, and subtract the master bias.
  *  From a central region of each mean frame, determine the median level.
  *  Determine also the variance within this region, using for this the
  *  difference of the images in the pair.
  */

  cpl_msg_info(task, "Creating photon transfer curve...") ;
  cpl_msg_indent_more();

  winSizeX = SAMPLE_REGION;
  winSizeY = SAMPLE_REGION;
  npix = winSizeX * winSizeY / 4;

  for (i = 0; i < 4; i++) {
    flux_variance[i] = newDpoint(pairCount);
    time_flux[i] = newDpoint(pairCount);
  }

  pairMean = (VimosImage **)cpl_calloc(pairCount, sizeof(VimosImage *));
  pair = (VimosImage **)cpl_calloc(2, sizeof(VimosImage *));
  noise = (float *)cpl_calloc(flatCount, sizeof(float));

  newPairCount = 0;

  for (i = 0; i < pairCount; i++) {

    loadFitsData(flatList1[i]);
    loadFitsData(flatList2[i]);

    pair[0] = flatList1[i];
    pair[1] = flatList2[i];

    noise[2 * i] = computeAverageRon(pair[0]);             /* RON in ADU */
    noise[2 * i + 1] = computeAverageRon(pair[1]);

    meanImage = frCombAverage(pair, 2);
    copyAllDescriptors(pair[0]->descs, &(meanImage->descs));

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

    if (VmSubBias(meanImage, biasImage, BIAS_ZMASTER) == EXIT_FAILURE) {
      for (j = i; j < pairCount; j++)
        deleteImage(flatList1[j]);
      cpl_free(flatList1);
      for (j = i; j < pairCount; j++)
        deleteImage(flatList2[j]);
      cpl_free(flatList2);
      cpl_free(pair);
      for (j = 0; j < newPairCount; j++)
        deleteImage(pairMean[j]);
      cpl_free(pairMean);
      deleteImage(meanImage);
      deleteImage(biasImage);
      cpl_free(keep);
      return EXIT_FAILURE;
    }


    /*
     *  Extract four central regions.
     */

    for (j = 0; j < 4; j++) {
      if (j == 0) {
        winStartX = (meanImage->xlen - winSizeX) / 2;
        winStartY = (meanImage->ylen - winSizeY) / 2;
        winPrscStartX = (pair[0]->xlen - winSizeX) / 2;  /* With prescans */
        winPrscStartY = (pair[1]->ylen - winSizeY) / 2;
      }
      if (j == 1) {
        winStartX += winSizeX / 2;
        winPrscStartX += winSizeX / 2;
      }
      if (j == 2) {
        winStartY += winSizeY / 2;
        winPrscStartY += winSizeY / 2;
      }
      if (j == 3) {
        winStartX -= winSizeX / 2;
        winPrscStartX -= winSizeX / 2;
      }
  
      sample1 = extractFloatImage(meanImage->data,
                                  meanImage->xlen, meanImage->ylen,
                                  winStartX, winStartY, 
                                  winSizeX / 2, winSizeY / 2);
  
      flux_variance[j][i].x = medianPixelvalue(sample1, npix);
      time_flux[j][i].y = flux_variance[j][i].x;
  
      cpl_free(sample1);
  
      sample1 = extractFloatImage(pair[0]->data, pair[0]->xlen, pair[0]->ylen,
                                  winPrscStartX, winPrscStartY, 
                                  winSizeX / 2, winSizeY / 2);
  
      sample2 = extractFloatImage(pair[1]->data, pair[1]->xlen, pair[1]->ylen,
                                  winPrscStartX, winPrscStartY, 
                                  winSizeX / 2, winSizeY / 2);
  
      diff = cpl_malloc(npix * sizeof(float));
  
      for (k = 0; k < npix; k++)
        diff[k] = sample1[k] - sample2[k];
  
      flux_variance[j][i].y = computeVarianceFloat2D(diff, 
                                                     winSizeX / 2, 
                                                     winSizeY / 2) / 2;

      readFloatDescriptor(pair[0]->descs, pilTrnGetKeyword("ExposureTime"),
                          &expTime, NULL);

      time_flux[j][i].x = expTime;
  
      if (j == 0)
        cpl_msg_info(task, 
        "Exptime, level, and variance of pair %d:\n"
        "%5.2f s, %7.1f ADU, %7.1f ADU", 
        i + 1, expTime, flux_variance[j][i].x, flux_variance[j][i].y);
  
      cpl_free(diff);
      cpl_free(sample1);
      cpl_free(sample2);
    }

    if (keep[i]) {
      pairMean[newPairCount] = meanImage;
      newPairCount++;
    }
    else {
      deleteImage(meanImage);
    }

    closeFitsImage(flatList1[i], 0);
    deleteImage(flatList1[i]);
    closeFitsImage(flatList2[i], 0);
    deleteImage(flatList2[i]);

  }

  deleteImage(biasImage);
  cpl_free(pair);
  cpl_free(flatList1);
  cpl_free(flatList2);


  /*
   *  It is possible at this point to evaluate the gain by linear 
   *  fitting of the exposure/variance table. The gain is measured 
   *  in e-/ADU.
   *  The fitting of the linearity curve is also performed. In the
   *  absence of non-linear terms, the relation between exposure
   *  time T and measured mean level L is
   *
   *                   L / (T - s) = C
   *  where s is the shutter delay, and C the count rate (constant,
   *  if the lamp is constant...). Therefore, we fit time vs mean level
   *  i.e., the linear relation
   *
   *                   L = CT - sC
   *
   *  Below, fit[1] = C and fit[0] = - sC, therefore s = - fit[0] / fit[1]
   */

  for (j = 0; j < 4; j++) {
    fit = fitSlopeRobust(flux_variance[j], pairCount);
    gain[j] = 1 / fit[1];
    cpl_free(fit);
    fit = fitSlopeRobust(time_flux[j], pairCount);
    cpl_free(fit);
    if (fit[2] > 200) {
      cpl_msg_error(task, "The input data are invalid: either the flat "
                    "lamp was strongly variable during the observation "
                    "(e.g., it switched off before the end of some "
                    "exposures), or the detector response is very far "
                    "from being linear (or some of the input frames are "
                    "saturated).");
      return EXIT_FAILURE;
    }
    fit = fit1DPoly(2, time_flux[j], pairCount, &rms);
    delay1 = -2 * fit[0] 
             / (fit[1] + sqrt(fit[1] * fit[1] - 4 * fit[2] * fit[0]));
    delay2 = -2 * fit[0] 
             / (fit[1] - sqrt(fit[1] * fit[1] - 4 * fit[2] * fit[0]));
    cpl_free(fit);
    if (fabs(delay1) < fabs(delay2))
      delay[j] = delay1;
    else
      delay[j] = delay2;
  }

  mean_gain = computeAverageFloat(gain, 4);
  value = 0.0;
  for (j = 0; j < 4; j++)
    value += (gain[j] - mean_gain) * (gain[j] - mean_gain);

  value /= 3;
  err_gain = sqrt(value);

  mean_delay = computeAverageFloat(delay, 4);
  value = 0.0;
  for (j = 0; j < 4; j++)
    value += (delay[j] - mean_delay) * (delay[j] - mean_delay);

  value /= 3;
  err_delay = sqrt(value);

  mean_noise = computeAverageFloat(noise, flatCount);
  value = 0.0;
  for (j = 0; j < flatCount; j++)
    value += (noise[j] - mean_noise) * (noise[j] - mean_noise);

  value /= flatCount - 1;
  err_noise = sqrt(value);

  mean_noise *= mean_gain;   /* In electrons */
  err_noise *= mean_gain;

  cpl_msg_indent_more();
  cpl_msg_info(task, "Gain          = %f +/- %f e-/ADU", mean_gain, err_gain);
  cpl_msg_info(task, "RON           = %f +/- %f e-", mean_noise, err_noise);
  cpl_msg_info(task, "Shutter delay = %f +/- %f s", mean_delay, err_delay);
  cpl_msg_indent_less();
  cpl_msg_indent_less();


  /*
   *  Here begins the search for bad pixels.
   */

  xlen = pairMean[0]->xlen;
  ylen = pairMean[0]->ylen;

  exposure_flux = newDpoint(newPairCount);


  /*
   *  Initializing the (optional) auxiliary image that will contain the
   *  rms residual for each line fit, and the output bad pixel image:
   */

  if (createError)
    errorImage = newImageAndAlloc(xlen, ylen);

  badPixelImage = newImageAndAlloc(xlen, ylen);


  /*
   * Create initial bad pixel map header from the input flat field
   */

  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("MjdObs"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("DateObs"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("ArchiveFile"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("Instrument"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("OBS.DID"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("OBS.ID"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("PROG.ID"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               "^ESO OCS (DID|CON QUAD)", NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               "^ESO TPL [.DINPSV.]", NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("INS.DID"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("InstrumentMode"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("DET.DID"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("NumberOfWindows", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("SeqWindowStartX", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("SeqWindowStartY", 1), NULL);

  /*
   * The bad pixel map does not contain pre- and overscan areas,
   * i.e. the window size is defined to be the image size.
   */

  writeIntDescriptor(&badPixelImage->descs,
                     pilTrnGetKeyword("SeqWindowSizeX", 1),
                     xlen,
                     pilTrnGetComment("SeqWindowSizeX"));
  writeIntDescriptor(&badPixelImage->descs,
                     pilTrnGetKeyword("SeqWindowSizeY", 1),
                     ylen,
                     pilTrnGetComment("SeqWindowSizeY"));

  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("WINi.BINX", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("WINi.BINY", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("READ.MODE"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("READ.SPEED"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("READ.CLOCK"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("NCHIPS"), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("CHIPi.INDEX", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("CHIPi.ID", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("CHIPi.NAME", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("CHIPi.DATE", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("CHIPi.NX", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("CHIPi.NY", 1), NULL);
  vimosDscCopy(&badPixelImage->descs, pairMean[0]->descs,
               pilTrnGetKeyword("Quadrant"), NULL);


  /* 
   * Store all exposure times (actually, the flatfield median 
   * value is used here as "time").
   */

  j = 0;
  for (i = 0; i < pairCount; i++) {
    if (keep[i]) {
      exposure_flux[j].x = flux_variance[0][i].x;
  
      /*
       * Keep track of the flat field with max exposure:
       */

      if (max < exposure_flux[j].x) {
        max = exposure_flux[j].x;
        maxFlat = j;
      }
      j++;
    }
  }

  cpl_free(keep);

  for (i = 0; i < 4; i++)
    deleteDpoint(flux_variance[i]);

  cpl_msg_info(task, "Analyzing pixels responses...");

  if (keep_flat < pairCount)
    cpl_msg_info(task, 
               "Using %d of %d flat fields pairs for bad pixels detection",
               keep_flat, pairCount);

  for (y = 0; y < ylen; y++) {
    for (x = 0; x < xlen; x++) {

      /* 
       * Getting all values for that position from each flat,
       * determining also the max among these values. 
       */

      for (k = 0; k < newPairCount; k++)
        exposure_flux[k].y = pairMean[k]->data[x + y * xlen];

      if (tolerant) {

        /*
         * Normalize pix values to their value corresponding to the
         * max exposed flat field.
         */

        if (exposure_flux[maxFlat].y > MIN_DIVISOR) {
          ratio = exposure_flux[maxFlat].x / exposure_flux[maxFlat].y;
          for (k = 0; k < newPairCount; k++)
            exposure_flux[k].y *= ratio;
        }

      }

/* debug
for (k = 0; k < newPairCount; k++)
printf("YYYY exposure-flux: %f %f\n", exposure_flux[k].x, exposure_flux[k].y);

printf("--------------------------------------------------------------\n");
*/


      /* 
       *  A line is fit to exposure times versus counts
       */

      fit = fitSlopeRobust(exposure_flux, newPairCount);


      /* 
       * The slope (fit[1]) is temporarily stored in the bad pixel image
       * The rms in the error image.
       */

      badPixelImage->data[x + y * xlen] = fit[1];

      if (errorImage) 
        errorImage->data[x + y * xlen] = fit[2];

      cpl_free(fit);

    }
  }


  /* 
   *  Now do some statistics to catch the bad pixels.
   *  The illumination of the flat fields may not be spatially
   *  uniform. The pixels response may depend on the illumination 
   *  level of the region they belong to. This would modify the
   *  mean response slope, and the standard deviation around this
   *  mean. We have to keep into account this effect, to adjust
   *  accordingly the absolute threshold used to identify the bad 
   *  pixels.
   */

  /*
   *  As an indicator of the flux variations within the flat fields,
   *  we take a smoothed image of any one of the flat fields.
   */

  filtered = VmFrMedFil(pairMean[0], 11, 11, 0);

  /*
   *  Here we define a chessboard of square regions (cells) where 
   *  their central flux will be compared with the local mean slope 
   *  of the pixel response and its variance.
   */

  nSidesX = xlen / CELLSIZE;
  nSidesY = ylen / CELLSIZE;

  flux_sigma = newDpoint(nSidesX * nSidesY);
  flux_average = newDpoint(nSidesX * nSidesY);

  count = 0;        /*  Cells counter */


  /*
   *  Keep far from image edges (usually garbage is there)
   */

  for (x = 1, boxPosX = x * CELLSIZE; 
       x < nSidesX - 1; x++, boxPosX = x * CELLSIZE) {
    for (y = 1, boxPosY = y * CELLSIZE; 
         y < nSidesY - 1; y++, boxPosY = y * CELLSIZE) {

      /*
       *  For each cell we find the mean and the variance of the line 
       *  slopes, then exclude all the values departing from the mean 
       *  more than NSIGMA sigma, repeat until nothing is excluded 
       *  anymore (or MAXITER iterations are done), derive the sigma 
       *  of the filtered distribution, then select bad pixels according 
       *  to the user defined threshold. This can work only if the bad 
       *  pixels are less than the good pixels.
       */

      nBoxPix = CELLSIZE * CELLSIZE;

      cell = extractFloatImage(badPixelImage->data, xlen, ylen,
                               boxPosX, boxPosY, CELLSIZE, CELLSIZE);


      /*
       *  Let the first guess for the average be a median: more robust.
       */

      average = medianPixelvalue(cell, nBoxPix);

      excluded = 1;              /* Just to let it enter the while block */
      nIter = 0;

      while (excluded > 0 && nIter < MAXITER) {

        /* 
         * First calculate the NSIGMA level
         */

        tolerance = 0.0;
        for (j = 0; j < nBoxPix; j++) 
          tolerance += fabs(cell[j] - average);

        tolerance /= nBoxPix;
        tolerance *= MEANDEV_TO_SIGMA * NSIGMA;
    
        cpl_msg_debug(task, "Average q-image: %f", average);
        cpl_msg_debug(task, "Tolerance      : %f", tolerance);


        /* 
         * Then we exclude from sample values departing from average more
         * than NSIGMA sigmas.
         */

        excluded = 0;
        newAverage = 0;
        for (j = 0; j < nBoxPix; j++) {
          delta = fabs(cell[j] - average);
          if (delta < tolerance) {
            cell[j - excluded] = cell[j];
            newAverage = (((double) j) / (j + 1)) * newAverage 
                       + cell[j] / (j + 1);
          }
          else {
            excluded++;
          }
        }
        average = newAverage;
        nBoxPix -= excluded;
        cpl_msg_debug(task,"Excluded       : %d\n",excluded);
        nIter++;
      }

      cpl_free(cell);

      /* 
       *  Setting here the user-specified tolerance.  
       */

      tolerance *= badPixelThresh;
      tolerance /= NSIGMA;

      cpl_msg_debug(task, "Tolerance      : %f", tolerance);


      /*
       *  Take the flux corresponding to current cell from its center,
       *  and put it in relation with the average slope and the
       *  absolute threshold:
       */

      cell = extractFloatImage(pairMean[0]->data, xlen, ylen,
                               boxPosX, boxPosY, CELLSIZE, CELLSIZE);

      flux_sigma[count].x   = medianWirth(cell, nBoxPix);
      flux_average[count].x = flux_sigma[count].x;
      flux_average[count].y = average;
      flux_sigma[count].y   = tolerance;

      cpl_free(cell);

      count++;
    }
  }      /*  End of loop on cells  */


  /*
   *  Here determine the relation between flux and average slope
   *  and slope variance.
   */

  averageCurve = fit1DPoly(DEGREE, flux_average, count, NULL);
  sigmaCurve = fit1DPoly(DEGREE, flux_sigma, count, NULL);

  deleteDpoint(flux_average);
  deleteDpoint(flux_sigma);

  /* 
   * Coding here the Bad Pixel Map
   */

  npix = xlen * ylen;

  for (j = 0; j < npix; j++) {

    /*
     *  For each pixel determine the expected average of the slope, 
     *  and the absolute tolerance of a deviation from this value, 
     *  on the basis of the illumination on that position.
     */

    tolerance = 0.0;
    average = 0.0;
    tmp = 1.0;
    for (x = 0; x <= DEGREE; x++) {
      tolerance += sigmaCurve[x] * tmp;
      average   += averageCurve[x] * tmp;
      tmp *= filtered->data[j];            /* Illumination */
    }

    delta = fabs(badPixelImage->data[j] - average);

    if (delta > tolerance) {
      badPixelImage->data[j] = 1;
    }
    else {
      badPixelImage->data[j] = 0;
    }
  }

  deleteImage(filtered);


  cpl_free(averageCurve);
  cpl_free(sigmaCurve);

  deleteDpoint(exposure_flux);


  /*
   *  Create a CCD bad pixel table from bad pixel image
   */

  if (!(ccdTable = badPixelImage2CcdTable(badPixelImage))) {
    cpl_msg_error(task, "Cannot create CCD table");
    for (i = 0; i < newPairCount; i++)
      deleteImage(pairMean[i]);
    cpl_free(pairMean);
    deleteImage(badPixelImage);
    deleteImage(errorImage);
    return EXIT_FAILURE;
  }

  if (!createMap) {
    deleteImage(badPixelImage);
    badPixelImage = NULL;
  }


  /* 
   * Write the gain + RON values to output CCD table header.
   */

  insertFloatDescriptor(&(ccdTable->descs),
                        pilTrnGetKeyword("ReadNoise", 1),
                        mean_noise,
                        pilTrnGetComment("ReadNoise"),
                        "ESO DET*", 0);

  insertFloatDescriptor(&(ccdTable->descs),
                        pilTrnGetKeyword("Electron2Adu", 1),
                        1 / mean_gain,
                        pilTrnGetComment("Electron2Adu"),
                        "ESO DET*", 0);

  insertFloatDescriptor(&(ccdTable->descs),
                        pilTrnGetKeyword("Adu2Electron", 1),
                        mean_gain,
                        pilTrnGetComment("Adu2Electron"),
                        "ESO DET*", 0);


  /*
   * Update other products header
   */

  if (createMap) {
    updateOK = (updateOK &&
                insertDoubleDescriptor(&(badPixelImage->descs), 
                                       pilTrnGetKeyword("DataMin"),
                                       0.0,
                                       pilTrnGetComment("DataMin"),
                                       "ESO*", 1));
  
    updateOK = (updateOK &&
                insertDoubleDescriptor(&(badPixelImage->descs), 
                                       pilTrnGetKeyword("DataMax"),
                                       1.0,
                                       pilTrnGetComment("DataMax"),
                                       "ESO*", 1));

    updateOK = (updateOK &&
                writeDoubleDescriptor(&(badPixelImage->descs), 
                                      pilTrnGetKeyword("DataMean"),
                                      imageMean(badPixelImage), 
                                      pilTrnGetComment("DataMean")));

    updateOK = (updateOK &&
                writeDoubleDescriptor(&(badPixelImage->descs),
                                      pilTrnGetKeyword("DataStdDeviation"),
                                      imageSigma(badPixelImage), 
                                      pilTrnGetComment("DataStdDeviation")));

    updateOK = (updateOK &&
                writeDoubleDescriptor(&(badPixelImage->descs),
                                      pilTrnGetKeyword("DataMedian"),
                                      0.0,
                                      pilTrnGetComment("DataMedian")));


    if (!updateOK) {
      cpl_msg_error(task, "Failure updating product header");
      for (j = 0; j < i; j++)
        deleteImage(pairMean[j]);
      cpl_free(pairMean);
      deleteImage(badPixelImage);
      deleteImage(errorImage);
      deleteTable(ccdTable);
      return EXIT_FAILURE;
    }

  }


  if (computeQC) {

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      frame = pilSofLookupNext(sof, flatCategory);

      firstImage = openOldFitsFile(pilFrmGetName(frame), 0, 0);
      loadFitsHeader(firstImage);

      /*
       * Write PRO.CATG, ARCFILE, TPL ID, and OCS CON QUAD, to QC1 group.
       */

      readIntDescriptor(firstImage->descs, pilTrnGetKeyword("Quadrant"),
                        &quadrant, NULL);

      pilQcWriteString("PRO.CATG", ccdTableCategory, "Product category");

      qcCopyValue(firstImage->descs, pilTrnGetKeyword("ArchiveFile"),
                  NULL, "Archive File Name");

      qcCopyValue(firstImage->descs, pilTrnGetKeyword("TplId"),
                  NULL, "Template signature ID");

      qcCopyValue(firstImage->descs, pilTrnGetKeyword("Quadrant"),
                  NULL, "Quadrant");
/*
 *    qcCopyValue(firstImage->descs, pilTrnGetKeyword("FilterName", quadrant),
 *                NULL, "Filter name");
 */
      closeFitsImage(firstImage, 0);
      deleteImage(firstImage);

      qcWriteValueDouble(ccdTable->descs, mean_delay, "QC.SHUTTER.DELAY", "s",
                         "Shutter delay");

      qcWriteValueDouble(ccdTable->descs, err_delay, "QC.SHUTTER.DELAY.ERROR", 
                         "s", "Error on shutter delay");

      if (pilQcGroupEnd() == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

    }
    else
      cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

  } /* End of QC1 computation. */


  /*
   * Create the product files on disk, set the products attributes and
   * update the set of frames.
   */

  vmstrlower(strcpy(ccdTableName, ccdTableCategory));
  /* strcat(ccdTableName, ".TFITS"); */
  strcat(ccdTableName, ".fits");

  if (createFitsTable(ccdTableName, ccdTable, ccdTableCategory)) {
    pilFitsHdrCopy(ccdTableName, 0, NULL, ".*-OBS$", 1);
    pilFitsHdrCopy(ccdTableName, 0, NULL, "^ESO .*", 1);

    outputFrame = newPilFrame(ccdTableName, ccdTableCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);

  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!",
                ccdTableName);

    for (i = 0; i < flatCount; i++)
      deleteImage(flatList[i]);
    cpl_free(flatList);
    deleteImage(badPixelImage);
    deleteImage(errorImage);
    deleteTable(ccdTable);

    return EXIT_FAILURE;
  }

  if (createMap) {
    vmstrlower(strcpy(badPixelImageName, badPixelCategory));
    strcat(badPixelImageName, ".fits");

    if (createFitsImage(badPixelImageName, badPixelImage, badPixelCategory)) {
      outputFrame = newPilFrame(badPixelImageName, badPixelCategory);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
    }
    else {
      cpl_msg_error(task, "Cannot create local product file %s!",
                  badPixelImageName);

      for (i = 0; i < flatCount; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);
      deleteImage(badPixelImage);
      deleteImage(errorImage);
      deleteTable(ccdTable);
  
      return EXIT_FAILURE;
    }
  }

  if (createError) {

   /*
    *  This image is created for debug reasons, and it is not inserted
    *  in the set of frames.
    */

    strcpy(errorImageName, "errorImage");
    strcat(errorImageName, ".fits");

    if (!createFitsImage(errorImageName, errorImage, "")) {
      cpl_msg_error(task, "Cannot create local debug product file %s!",
                  errorImageName);
  
      for (i = 0; i < flatCount; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);
      deleteImage(badPixelImage);
      deleteImage(errorImage);
      deleteTable(ccdTable);
  
      return EXIT_FAILURE;
    }
  }

 /*
  * Cleanup.
  */

  if (badPixelImage) 
    deleteImage(badPixelImage);
  deleteImage(errorImage);
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

    cpl_recipe *recipe = cpl_calloc(1, sizeof *recipe);
    cpl_plugin *plugin = &recipe->interface;


    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmdet",

    "Determine RON, gain, and bad pixels from a set of flat field exposures.",

    "This recipe is used to estimate the read-out-noise (RON) and the\n"
    "gain of the CCD, and to determine the positions of the bad pixels.\n"
    "The input SOF should contain at least five pairs of flat field\n"
    "exposures, all belonging to the same quadrant, each pair corresponding\n"
    "to a different exposure time. The flat fields can be produced either\n"
    "in imaging or in MOS mode. In MOS mode a HR grism is used, in order to\n"
    "illuminate the CCD also beyond the central region used in direct\n"
    "imaging mode, but no mask is inserted at the telescope focal plane.\n"
    "This type of exposure cannot really be considered a spectral flat\n"
    "field, because the CCD is exposed to white light (i.e., without\n"
    "a wavelength dependency along the dispersion direction). The flat\n"
    "fields generated for the purpose of determining the detector properties\n"
    "are assigned the DO category DETECTOR_PROPERTIES, to distinguish them\n"
    "from the more common IMG_SCREEN_FLAT or MOS_SCREEN_FLAT that are used\n"
    "to produce master calibrations.\n\n"
    "Input files:\n\n"
    "  DO category:          Type:         Explanation:  Required:\n"
    "  DETECTOR_PROPERTIES   Raw           Flat field       Y\n"
    "  MASTER_BIAS           Calib         Master bias      Y\n\n"
    "Output files:\n\n"
    "  DO category:          Data type:    Explanation:\n"
    "  CCD_TABLE             FITS table    Bad pixel table\n"
    "  (none)                FITS image    Bad pixel image\n"
    "  (none)                FITS image    Error image\n\n"
    "Only the primary product, the bad pixel table, is copied (or moved)\n"
    "to the product directory. Other products are generated only on request\n"
    "(typically for debug purposes) and are not assigned a DO category as\n"
    "they would not be used anywhere in further data processing steps.\n\n"
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

                    vmdet_create,
                    vmdet_exec,
                    vmdet_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
