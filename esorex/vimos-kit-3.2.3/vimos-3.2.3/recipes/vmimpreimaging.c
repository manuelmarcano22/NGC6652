/* $Id: vmimpreimaging.c,v 1.8 2012-11-08 17:56:57 cgarcia Exp $
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
 * $Date: 2012-11-08 17:56:57 $
 * $Revision: 1.8 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <string.h>

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

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmutils.h"
#include "vmimgpreprocessing.h"
#include "vmimgastrometry.h"
#include "vmimgphotcalib.h"
#include "vmimgutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


/*
 * Definition of the label strings for all methods the recipe function
 * supports for removing the bias.
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


static cxint vmimpreimaging(PilSetOfFrames *);


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmimpreimaging_create(cpl_plugin *plugin)
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
                                "Bad pixel correction on reduced science "
                                "image.",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanBadPixel");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanBadPixel");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.clean",
                                CPL_TYPE_BOOL,
                                "Cosmic ray removal from the input "
                                "science image.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CleanCosmic");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CleanCosmic");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.threshold",
                                CPL_TYPE_DOUBLE,
                                "Threshold for cosmic ray candidate "
                                "selection.",
                                "vimos.Parameters",
                                4.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CosmicsThreshold");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CosmicsThreshold");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.cosmics.ratio",
                                CPL_TYPE_DOUBLE,
                                "Ratio for discrimination of objects "
                                "and cosmics.",
                                "vimos.Parameters",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CosmicsRatio");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CosmicsRatio");
    cpl_parameterlist_append(recipe->parameters, p);

#ifdef ONLINE_MODE

    p = cpl_parameter_new_value("vimos.Parameters.comatrix.radius",
                                CPL_TYPE_DOUBLE,
                                "Aperture used for object identification "
                                "during computation of the CO matrix",
                                "vimos.Parameters",
                                6.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "COMatrixSearchRadius");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "COMatrixSearchRadius");
    cpl_parameterlist_append(recipe->parameters, p);

#endif

    p = cpl_parameter_new_value("vimos.Parameters.temperature.check",
                                CPL_TYPE_BOOL,
                                "Check beam temperature when updating "
                                "the world coordinate system.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "TemperatureCheck");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "TemperatureCheck");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.temperature.tolerance",
                                CPL_TYPE_DOUBLE,
                                "Maximum allowed difference between beam "
                                "and ambient temperature.",
                                "vimos.Parameters",
                                3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "TemperatureTolerance");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "TemperatureTolerance");
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
vmimpreimaging_exec(cpl_plugin *plugin)
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

    if (vmimpreimaging(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmimpreimaging");
        
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
vmimpreimaging_destroy(cpl_plugin *plugin)
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
 * @brief
 *   Processing of a preimaging observation from a single quadrant.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise @c EXIT_FAILURE is returned.
 *
 * @param sof   Set of frames containing the references to a preimaging
 *              observation, a master bias and a master flat field. 
 *              Optionally a photometric table and a CCD table may
 *              be added, to apply the photometric calibration and/or
 *              bad pixel cleaning.
 *
 * This function provides the the basic reduction steps for preimaging
 * observations. This includes at least the removal of the bias, and 
 * the flat field correction. Optionally, cleaning of bad pixels and 
 * cosmic ray hits can be performed, and a photometric calibration 
 * applied.
 * 
 * Finally, for science and standard star fields the CD and the CCD to
 * sky transformation matrix are used to create the CO matrix which is
 * written to the image header.
 * 
 * On successful termination the engine creates a reduced image.
 * 
 * Control options and additional parameters are read from the recipe
 * configuration database. The recipe function accepts the following
 * task parameters:
 *
 * @li BiasMethod        Method used for bias removal.
 * @li CleanBadPixel     Bad pixel correction on input image.
 * @li CleanCosmic       Cosmic ray removal from the input image.
 * @li CosmicsRatio      Threshold for cosmic ray candidate selection.
 * @li CosmicsThreshold  Ratio for discrimination of objects and cosmics.
 *
 * If any of these task parameters is not set in the recipe configuration
 * database the recipe function uses the builtin defaults for these
 * task parameters.
 * 
 * @author C.Izzo, R.Palsa
 */   

static cxint
vmimpreimaging(PilSetOfFrames *sof)
{

  const char  task[] = "vmimpreimaging";

  char *biasMethodTag = 0;
  char *rawTag = (char *)pilTrnGetCategory("ImgPreimaging");
  char *productTag = (char *)pilTrnGetCategory("redImgScience");
  char productName[PATHNAME_MAX + 1];

  unsigned int cleanBadPixel, cleanCosmic;
  unsigned int tempCheck;

  int biasMethodEntry;
  int scienceCount = 0;

  float thresholdCosmics, ratioCosmics;

  double coMatrixSearchRadius;
  double tempTolerance;

  PilFrame *ccdFrame, *biasFrame;
  PilFrame *flatFrame, *ipcFrame;
  PilFrame *rawFrame;
  PilFrame *productFrame;

  VimosImage *rawImage;
  VimosImage *biasImage, *flatImage;
  VimosImage *resizedBias = 0;
  VimosImage *rawImageFF = 0;

  VimosTable *ccdTable = 0;
  VimosTable *ipcTable = 0;

  BiasMethod biasMethod = BIAS_UNDEF;


 /*
  * Get task parameters from the recipe database
  */

 /*
  * Retrieve the bias removal method
  */

  biasMethodTag = (char *)pilDfsDbGetString("Parameters", "BiasMethod");

  if ((biasMethodEntry = strselect(biasMethodTag, biasMethodNames,
                                  nBiasMethods)) < 0) {
    cpl_msg_error(task, "%s: Invalid bias removal method!", biasMethodTag);
    return EXIT_FAILURE;
  }

  biasMethod = biasMethods[biasMethodEntry];


 /*
  * Check if bad pixel correction is required.
  */

  cleanBadPixel = pilDfsDbGetBool("Parameters", "CleanBadPixel", 1);


 /*
  * Check if cosmic ray hit correction is required.
  */

  cleanCosmic = pilDfsDbGetBool("Parameters", "CleanCosmic", 1);
  thresholdCosmics = pilDfsDbGetFloat("Parameters", "CosmicsThreshold", 4.0);
  ratioCosmics = pilDfsDbGetFloat("Parameters", "CosmicsRatio", 2.0);


 /*
  * Retrieve the setup for matching the faked stars for distorsions.
  */

  coMatrixSearchRadius = 
          pilDfsDbGetDouble("Parameters", "COMatrixSearchRadius", 6.);

  if (coMatrixSearchRadius < 0.) {
    cpl_msg_error(task, "Search Radius is out of range!");
    return EXIT_FAILURE;
  }


  /*
   * Check whether temperature checks should be used when
   * updating the world coordinate system with a CO matrix
   * and get the temperature tolerance to apply.
   */

  tempCheck = pilDfsDbGetBool("Parameters", "TemperatureCheck", 1);
  tempTolerance = pilDfsDbGetDouble("Parameters", "TemperatureTolerance", 3.);


 /*
  * Make sure that all required input frames are present. Exactly
  * one preimaging observation from a single quadrant is allowed.
  */

  scienceCount = pilSofFrameCount(sof, rawTag);

  if (scienceCount == 0) {
    cpl_msg_error(task, "No image found in input!");
    return EXIT_FAILURE;
  }

  if (scienceCount > 1) {
    cpl_msg_error(task, "More than one image found in input!");
    return EXIT_FAILURE;
  }


 /*
  * Check for the master bias
  */

  if (!(biasFrame = pilSofLookup(sof, pilTrnGetCategory("MasterBias")))) {
    cpl_msg_error(task, "No master bias frame in input!");
    return EXIT_FAILURE;
  }

  pilFrmSetType(biasFrame, PIL_FRAME_TYPE_CALIB);

 /*
  * Check for the master flat field frame
  */

  if (!(flatFrame = pilSofLookup(sof,
                                 pilTrnGetCategory("ImgMasterSkyFlat")))) {
    cpl_msg_error(task, "No master flat field frame in input!");
    return EXIT_FAILURE;
  }

  pilFrmSetType(flatFrame, PIL_FRAME_TYPE_CALIB);

 /*
  * If the bad pixel correction is requested a bad pixel table is
  * required in the input set of frames. If no bad pixel table can
  * be found this is an error. If cosmic ray cleaning is enabled
  * the bad pixel table is used, if it is present in the input set,
  * to avoid bad pixels in the surrounding when correcting for cosmics,
  * but a missing bad pixel table will not cause an error.
  */

  ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable"));
  if (ccdFrame)
    pilFrmSetType(ccdFrame, PIL_FRAME_TYPE_CALIB);

  if (cleanBadPixel || cleanCosmic) {
    if (!ccdFrame) {
      if (cleanBadPixel) {
        cpl_msg_error(task, 
                    "Bad pixel cleaning requires a CCD table in input!");
        return EXIT_FAILURE;
      }
    }
    else {
      cpl_msg_debug(task, "CCD table is %s", pilFrmGetName(ccdFrame));
    }
  }


 /*
  * Check for the photometric calibration data.
  */

  ipcFrame = pilSofLookup(sof, pilTrnGetCategory("PhotometricCoeffTable"));
  if (!ipcFrame)
    //Backwards-compatible. This was the tag used before Nov 2012. Probably it could be removed in the future
    ipcFrame = pilSofLookup(sof, pilTrnGetCategory("PhotometricTable"));
    if (!ipcFrame)
      cpl_msg_warning(task, "No photometric calibration data in input: "
                "no photometric calibration will be applied to the "
                "reduced data.");
  if (ipcFrame)
      pilFrmSetType(ipcFrame, PIL_FRAME_TYPE_CALIB);


 /*
  * Load the observation data.
  */

  if ((rawFrame = pilSofLookup(sof, rawTag))) {
    if ((rawImage = openOldFitsFile(pilFrmGetName(rawFrame), 1, 0)) == NULL) {
      cpl_msg_error(task, "Cannot load image %s!", pilFrmGetName(rawFrame));
      return EXIT_FAILURE;
    }
    closeFitsImage(rawImage, 0);
    pilFrmSetType(rawFrame, PIL_FRAME_TYPE_RAW);
  }
   

 /*
  * Remove the bias from the observation. The master bias does
  * not have overscan areas anymore so they are faked by enlarging
  * the master bias using the observation as reference.
  */

  cpl_msg_info(task, "Bias removal...");

  if (!(biasImage = openOldFitsFile(pilFrmGetName(biasFrame), 1, 0))) {
    cpl_msg_error(task, "Cannot load master bias %s!", pilFrmGetName(biasFrame));
    deleteImage(rawImage);
    return EXIT_FAILURE;
  }

  closeFitsImage(biasImage, 0);

  if (!(resizedBias = growOverscans(biasImage, rawImage))) {
    cpl_msg_error(task, "Restoring overscan regions failed!");
    deleteImage(biasImage);
    deleteImage(rawImage);
    return EXIT_FAILURE;
  }

  if (resizedBias != biasImage) {
    deleteImage(biasImage);
    biasImage = resizedBias;
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

  if (VmSubBias(rawImage, biasImage, biasMethod) == EXIT_FAILURE) {
    cpl_msg_error(task, "Bias correction failed!");
    deleteImage(biasImage);
    deleteImage(rawImage);
    return EXIT_FAILURE;
  }

  deleteImage(biasImage);


 /*
  * Flat field correction.
  */

  cpl_msg_info(task, "Flat field correction ...");

  if (!(flatImage = openOldFitsFile(pilFrmGetName(flatFrame), 1, 0))) {
    cpl_msg_error(task, "Cannot load master flat field %s!",
                pilFrmGetName(flatFrame));
    deleteImage(rawImage);
    return EXIT_FAILURE;
  }

  closeFitsImage(flatImage, 0);

  if (!(rawImageFF = imageArith(rawImage, flatImage, VM_OPER_DIV))) {
    cpl_msg_error(task, "Flat field correction failed!");
    deleteImage(flatImage);
    deleteImage(rawImage);
    return EXIT_FAILURE;
  }

  copyAllDescriptors(rawImage->descs, &(rawImageFF->descs));

  deleteImage(flatImage);
  deleteImage(rawImage);

  if (!findDescriptor(rawImageFF->descs, pilTrnGetKeyword("MaskCcdX0"))) {
    cpl_msg_error(task, "Missing Mask-to-CCD keywords in product!");
    return EXIT_FAILURE;
  }

  if (!findDescriptor(rawImageFF->descs, pilTrnGetKeyword("CcdMaskTemp"))) {
    cpl_msg_error(task, "Missing Mask-to-CCD temperature keyword in product!");
    return EXIT_FAILURE;
  }

  if (!findDescriptor(rawImageFF->descs, pilTrnGetKeyword("MaskCcdY", 0, 0))) {
    cpl_msg_error(task, "Missing Mask-to-CCD keywords in product!");
    return EXIT_FAILURE;
  }

  if (!findDescriptor(rawImageFF->descs, pilTrnGetKeyword("CcdMaskX", 0, 0))) {
    cpl_msg_error(task, "Missing CCD-to-Mask keywords in product!");
    return EXIT_FAILURE;
  }
  

 /*
  * Load bad pixel data in case bad pixel correction and/or cosmic ray
  * cleaning should be done.
  */

  if (cleanBadPixel || (cleanCosmic && ccdFrame)) {
    if (!(ccdTable = openOldFitsTable(pilFrmGetName(ccdFrame), 0))) {
      cpl_msg_error(task, "Cannot load bad pixel data from %s!",
                  pilFrmGetName(ccdFrame));
      return EXIT_FAILURE;
    }
    closeFitsTable(ccdTable, 0);
  }

 
 /*
  *  Bad pixel cleaning
  */
  
  if (cleanBadPixel) {
    cpl_msg_info(task, "Cleaning bad pixels ...");

    if (cleanBadPixels(rawImageFF, ccdTable, 0) == EXIT_FAILURE) {
      cpl_msg_error(task, "Bad pixel cleaning failed!");
      if (ccdTable)
        deleteTable(ccdTable);
      deleteImage(rawImageFF);
      return EXIT_FAILURE;
    }
  }

  
 /*
  * Correct for cosmic ray hits.
  */

  if (cleanCosmic) {

   /* FIXME:
    *  I would prefer the image mode here, assuming there is more 
    *  background than objects, but this requires cleaning up the 
    *  methods for histograms. Let's postpone it a little. (RP)
    */

    float bkgEstimate = imageMedian(rawImageFF);
    float gain = getMeanGainFactor(rawImageFF);
    float ron = computeAverageRon(rawImageFF);

    if (gain < MIN_DIVISOR) {
      cpl_msg_error(task, "Missing or invalid gain factor encountered in %s!",
                  pilFrmGetName(rawFrame));
      if (ccdTable)
        deleteTable(ccdTable);
      deleteImage(rawImageFF);
      return EXIT_FAILURE;
    }

    cpl_msg_info(task, "Removing cosmic ray hits ...");

    if (VmCosmicClean(rawImageFF, ccdTable, 0, bkgEstimate, gain, ron,
                      thresholdCosmics, ratioCosmics) == EXIT_FAILURE) {
      cpl_msg_error(task, "Removal of cosmic ray hits failed!");
      if (ccdTable)
        deleteTable(ccdTable);
      deleteImage(rawImageFF);
      return EXIT_FAILURE;
    }
  }

  if (ccdTable)
    deleteTable(ccdTable);


 /* 
  * If the photometric table is available, photometric calibration is 
  * applied. This just updates header keywords. The image data is not 
  * altered in any way.
  */
  
  if (ipcFrame) {

    cpl_msg_info(task, "Performing photometric calibration ...");

    if (!(ipcTable = openOldFitsTable(pilFrmGetName(ipcFrame), 0))) {
      cpl_msg_error(task, "Cannot load photometric calibration data from %s!",
                  pilFrmGetName(ipcFrame));
      return EXIT_FAILURE;
    }

    closeFitsTable(ipcTable, 0);

    if (!(VmImApplyPhot(rawImageFF, ipcTable))) {
      cpl_msg_error(task, "Photometric calibration failed!");
      if (ipcTable)
        deleteTable(ipcTable);
      deleteImage(rawImageFF);
      return EXIT_FAILURE;
    }

    deleteTable(ipcTable);

  }

/********************  %%% */

  if (findDescriptor(rawImageFF->descs, pilTrnGetKeyword("CD", 1, 1))) {
    if (findDescriptor(rawImageFF->descs, pilTrnGetKeyword("Cdelt", 1)))
      removeDescriptor(&rawImageFF->descs, pilTrnGetKeyword("Cdelt", 1));

    if (findDescriptor(rawImageFF->descs, pilTrnGetKeyword("Cdelt", 2)))
      removeDescriptor(&rawImageFF->descs, pilTrnGetKeyword("Cdelt", 2));
  }
/********************  %%% */

 /*
  * Write the CO matrix to the header of the reduced image. It is
  * generated from the CD matrix and the sky to CCD transformation
  * found in the image header, i.e. no star match table needs to
  * be passed.
  */

 /* FIXME:
  *   The value 10.0 is passed as sigma clipping factor which is
  *   arbitrarily choosen. BG just said it should be a high value
  *   since the tables are created artificially. But in this case it
  *   would be better to remove it completely from the call and
  *   handle this internally choosing a value appropriate for the
  *   create tables. (RP)
  */

  cpl_msg_info(task, "Updating image world coordinate system ...");

  if (VmAstroComputeCO(rawImageFF, coMatrixSearchRadius, 10., 0L,
                       tempCheck, tempTolerance) != VM_TRUE) {
    cpl_msg_error(task, "CO matrix computation failed!");
    deleteImage(rawImageFF);
    return EXIT_FAILURE;
  }


 /*
  * Prepare the products. This includes creation of the temporary product
  * file name from its product category, update of the product header,
  * and the creation of the temporary product file in the current working
  * directory. Finally each product has to be registered with the
  * appropriate attributes in the set of frames.
  */

 /* FIXME:
  * For the moment also keywords which are not task specific
  * are handled here, since this is the last possibility to access
  * the linked list of keywords without reopening the file.
  * This will change in future!
  */


 /*
  * Reduced image
  */

  vmstrlower(strcpy(productName, productTag));
  strcat(productName, ".fits");


  insertDoubleDescriptor(&rawImageFF->descs,
                         pilTrnGetKeyword("DataMin"),
                         imageMinimum(rawImageFF),
                         pilTrnGetComment("DataMin"),
                         "ESO*", 1);

  insertDoubleDescriptor(&rawImageFF->descs,
                         pilTrnGetKeyword("DataMax"),
                         imageMaximum(rawImageFF),
                         pilTrnGetComment("DataMax"),
                         "ESO*", 1);

  insertDoubleDescriptor(&rawImageFF->descs,
                         pilTrnGetKeyword("DataMedian"),
                         imageMedian(rawImageFF), 
                         pilTrnGetComment("DataMedian"),
                         "ESO PRO*", 1);

  insertDoubleDescriptor(&rawImageFF->descs, 
                         pilTrnGetKeyword("DataStdDeviation"),
                         imageSigma(rawImageFF),
                         pilTrnGetComment("DataStdDeviation"),
                         "ESO PRO*", 1);

  insertDoubleDescriptor(&rawImageFF->descs,
                         pilTrnGetKeyword("DataMean"),
                         imageMean(rawImageFF),
                         pilTrnGetComment("DataMean"),
                         "ESO PRO*", 1);

  deleteSetOfDescriptors(&rawImageFF->descs, "ESO DPR*");

  if (createFitsImage(productName, rawImageFF, productTag) != VM_TRUE) {
    cpl_msg_error(task, "Cannot create local product file %s!", productName);
    deleteImage(rawImageFF);
    return EXIT_FAILURE;
  }

  productFrame = newPilFrame(productName, productTag);

  pilFrmSetType(productFrame, PIL_FRAME_TYPE_PRODUCT);
  pilFrmSetFormat(productFrame, PIL_FRAME_FORMAT_IMAGE);
  pilFrmSetProductLevel(productFrame, PIL_PRODUCT_LEVEL_PRIMARY);
  pilFrmSetProductType(productFrame, PIL_PRODUCT_TYPE_REDUCED);

  pilSofInsert(sof, productFrame);

  deleteImage(rawImageFF);

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
                    "vmimpreimaging",
    "Reduce a preimaging science exposure.",
    "This recipe is used to apply basic reduction steps to the imaging\n"
    "observation that is preliminary to a MOS observation of the same\n"
    "field. No source detection is attempted on the image. The image\n"
    "WCS, together with the component describing the instrument optical\n"
    "distortions, is converted into the convention followed by VMMPS\n" 
    "(the VIMOS mask preparation software). This set of coefficients,\n" 
    "the so-called CO-matrix used by the SAO WCSTools package, is written\n"
    "to the header of the reduced image.\n\n"
    "Input files:\n\n"
    "  DO category:              Type:       Explanation:         Required:\n"
    "  IMG_PREIMAGING            Raw         Preimaging exposure     Y\n"
    "  MASTER_BIAS               Calib       Master bias             Y\n"
    "  IMG_MASTER_SKY_FLAT       Calib       Master sky flat         Y\n"
    "  PHOT_COEFF_TABLE          Calib       Photometric table       .\n"
    "  CCD_TABLE                 Calib       Bad pixel table         .\n\n"
    "Output files:\n\n"
    "  DO category:              Data type:  Explanation:\n"
    "  IMG_SCIENCE_REDUCED       FITS image  Reduced preimaging exposure\n\n"
    "If a photometric table is specified, the magnitude zeropoint, the\n"
    "atmospheric extinction coefficient, and the colour term are copied\n"
    "from the photometric table to the header of the reduced image. A CCD\n"
    "table must be specified in input only if a bad pixel cleaning is\n"
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

                    vmimpreimaging_create,
                    vmimpreimaging_exec,
                    vmimpreimaging_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
