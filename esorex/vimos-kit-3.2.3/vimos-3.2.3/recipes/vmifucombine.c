/* $Id: vmifucombine.c,v 1.3 2011-03-14 15:28:58 cizzo Exp $
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


#include <string.h>
#include <math.h>

#include <cxmemory.h>

#include <cpl_recipe.h>
#include <cpl_plugininfo.h>
#include <cpl_parameterlist.h>
#include <cpl_frameset.h>
#include <cpl_propertylist.h>

#include <cpl_memory.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
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
#include "vmcpl.h"
#include "vimos_dfs.h"


static cxint vmifucombine(PilSetOfFrames *);


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmifucombine_create(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;

/*    cpl_parameter *p;   */

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

/*
    p = cpl_parameter_new_value("vimos.Parameters.quality.enable",
                                CPL_TYPE_BOOL,
                                "Compute QC1 parameters",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ComputeQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ComputeQC");
    cpl_parameterlist_append(recipe->parameters, p);
*/


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
vmifucombine_exec(cpl_plugin *plugin)
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

    if (vmifucombine(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmifucombine");
        
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
vmifucombine_destroy(cpl_plugin *plugin)
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
 *   Combine IFU FOV reconstructed from different quadrants.
 * 
 * @return EXIT_SUCCESS or EXIT_FAILURE.   
 *
 * @param sof  Set of frames containing from 2 to 4 IFU_FOV belonging to
 *             different VIMOS quadrants, generated by the vmifuscience
 *             or by the vmifustandard recipes.
 *  
 * @doc
 *   This recipe would simply sum together the images specified in
 *   the input SOF. All images must belong to a different VIMOS
 *   quadrant (and therefore not more than 4 images may be specified
 *   in input). Any input image must have size 80x80. If the inputs
 *   have not the same value for the TPL.START keyword, a warning
 *   is issued, but the images are combined anyway. The PRO.CATG
 *   of the product is still IFU_FOV.
 *
 * @author C. Izzo, R. Palsa
 */

static cxint 
vmifucombine(PilSetOfFrames *sof)
{

  const char     task[] = "vmifucombine";

/*  const char     parameter[] = "Parameters";  */

  const char    *fovTag     = pilTrnGetCategory("IfuFov");
  const char    *fovStdTag  = pilTrnGetCategory("IfuStdFov");
  const char    *fovFullTag = pilTrnGetCategory("IfuFullFov");
  char           fovName[PATHNAME_MAX + 1];
  VimosBool      updateOK = VM_TRUE;
  int            fovCount, fovStdCount;

  int            i, j, q;

  PilFrame      *fovFrame;
  PilFrame      *outputFrame;

  VimosImage   **fovList;

  cpl_propertylist *header;

  char           tplStart[24], atplStart[24];
  char          *p;
  int            quadrant[] = {0, 0, 0, 0};
  double         refnorm = 1;
  double         norm = 1;


  /*
   * Make sure that the input SOF include the right datasets.
   */

  fovCount = (int)pilSofFrameCount(sof, fovTag);
  fovStdCount = (int)pilSofFrameCount(sof, fovStdTag);

  if (fovCount > 0 && fovStdCount > 0) {
    cpl_msg_error(task, "Incompatible IFU FOVs in input.");
    return EXIT_FAILURE;
  }

  if (fovCount == 0) {
    fovCount = fovStdCount;
    fovTag = pilTrnGetCategory("IfuStdFov");
    fovFullTag = pilTrnGetCategory("IfuFullStdFov");
  }

  if (fovCount < 2) {
    cpl_msg_error(task, "Too few IFU FOVs in input: at least 2 are required.");
    return EXIT_FAILURE;
  }

  if (fovCount > 4) {
    cpl_msg_error(task, "Too many IFU FOVs in input: up to 4 may be given.");
    return EXIT_FAILURE;
  }


  cpl_msg_info(task, "Loading input frames...");

  fovList = (VimosImage **)cpl_calloc(fovCount, sizeof(VimosImage *));

  fovFrame = pilSofLookupNext(sof, fovTag);

  for (i = 0; i < fovCount; i++) {
    if ((fovList[i] = openOldFitsFile(pilFrmGetName(fovFrame), 1, 0))) {
      pilFrmSetType(fovFrame, PIL_FRAME_TYPE_CALIB);
      closeFitsImage(fovList[i], 0);
      readIntDescriptor(fovList[i]->descs, pilTrnGetKeyword("Quadrant"),
                        &q, NULL);
      --q;

      if (quadrant[q]) {
        cpl_msg_error(task, "Input frames belonging to the same quadrant.");
        if (i)
          deleteImage(fovList[i]);
        deleteImage(fovList[0]);
        cpl_free(fovList);
        return EXIT_FAILURE;
      }

      quadrant[q] = 1;

      if (fovList[i]->xlen != 80 || fovList[i]->ylen != 80) {
        cpl_msg_error(task, "Wrong input frame %s (not 80x80).", 
                    pilFrmGetName(fovFrame));
        if (i)
          deleteImage(fovList[i]);
        deleteImage(fovList[0]);
        cpl_free(fovList);
        return EXIT_FAILURE;
      }

      if (i)
        p = atplStart;
      else
        p = tplStart;

      if (readStringDescriptor(fovList[i]->descs, 
                               "ESO TPL START", p, NULL) == VM_FALSE) {
        cpl_msg_error(task, "Missing TPL START keyword in file %s",
                    pilFrmGetName(fovFrame));
        if (i)
          deleteImage(fovList[i]);
        deleteImage(fovList[0]);
        cpl_free(fovList);
        return EXIT_FAILURE;
      }

      header = cpl_propertylist_load(pilFrmGetName(fovFrame), 0);
      if (cpl_propertylist_has(header, "ESO QC IFU FLAT FLUX")) {
        if (i) {
          norm = cpl_propertylist_get_double(header, "ESO QC IFU FLAT FLUX");
          norm /= refnorm;
        }
        else {
          refnorm = cpl_propertylist_get_double(header, "ESO QC IFU FLAT FLUX");
        }
      }
      cpl_propertylist_delete(header);

      if (i) {
        if (strcmp(tplStart, atplStart))
          cpl_msg_warning(task, "Input frames are from different exposures!");

        for (j = 0; j < 6400; j++)
          fovList[0]->data[j] += fovList[i]->data[j] / norm;

        deleteImage(fovList[i]);
      }

    }
    else {
      cpl_msg_error(task, "Failure loading IFU_FOV frame %d", i + 1);
      if (i)
        deleteImage(fovList[0]);
      cpl_free(fovList);
      return EXIT_FAILURE;
    }

    fovFrame = pilSofLookupNext(sof, NULL);

  }


  /*
   * Update the reconstructed FOV header
   */

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(fovList[0]->descs),
                                     pilTrnGetKeyword("DataMin"),
                                     imageMinimum(fovList[0]),
                                     pilTrnGetComment("DataMin"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(fovList[0]->descs),
                                     pilTrnGetKeyword("DataMax"),
                                     imageMaximum(fovList[0]),
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(fovList[0]->descs),
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(fovList[0]),
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(fovList[0]->descs),
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(fovList[0]),
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(fovList[0]->descs),
                                     pilTrnGetKeyword("DataMean"),
                                     imageMean(fovList[0]),
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    deleteImage(fovList[0]);
    return EXIT_FAILURE;
  }



  /*
   * Create the product file on disk, set the product attributes and
   * update the set of frames.
   */

  vmstrlower(strcpy(fovName, fovFullTag));
  strcat(fovName, ".fits");

  if (createFitsImage(fovName, fovList[0], fovFullTag)) {
    outputFrame = newPilFrame(fovName, fovFullTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", fovName);
    deleteImage(fovList[0]);
    return EXIT_FAILURE;
  }

  deleteImage(fovList[0]);

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
                    "vmifucombine",
    "Combine partial IFU field-of-views reconstructed by vmifuscience "
    "from different quadrants.",
    "This recipe is used to compose into a single image the reconstructed\n" 
    "images of the IFU field-of-view from different VIMOS quadrants.\n"
    "Such images are created by the recipes vmifuscience and vmifustandard.\n"
    "The input images must belong to different quadrants, so that they\n"
    "cannot be more than 4. The images are simply summed together: in\n"
    "fact, all field-of-view images are always 80x80 pixels, with only\n"
    "the sector corresponding to one quadrant having intensities that\n"
    "differ from zero.\n\n"
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

                    vmifucombine_create,
                    vmifucombine_exec,
                    vmifucombine_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
