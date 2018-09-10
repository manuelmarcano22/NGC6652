/* $Id: vmspphot.c,v 1.3 2011-03-14 15:28:58 cizzo Exp $
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
#include <stdlib.h>
#include <stdio.h>

#include <cpl_recipe.h>
#include <cpl_plugininfo.h>
#include <cpl_parameterlist.h>
#include <cpl_frameset.h>

#include <pilmemory.h>
#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilrecipe.h>
#include <piltranslator.h>
#include <cpl_msg.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmextincttable.h"
#include "vmspecphottable.h"
#include "vmutils.h"
#include "vmmossphotcalib.h"
#include "vmmosutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


static cxint vmspphot(PilSetOfFrames *);


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmspphot_create(cpl_plugin *plugin)
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


    p = cpl_parameter_new_value("vimos.Parameters.flux.response",
                                CPL_TYPE_BOOL,
                                "Apply instrument response.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyResponse");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyResponse");
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
vmspphot_exec(cpl_plugin *plugin)
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

    if (vmspphot(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmspphot");
        
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
vmspphot_destroy(cpl_plugin *plugin)
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
 * @param sof  Set of frames with references to at least one frame
 *             in 1D-stacked format (MOS_SCIENCE_REDUCED) created
 *             by vmmosobsstare, vmmosobsjitter, or vmifuscience, 
 *             and an atmospheric extinction table. A spectro-photometric 
 *             table should also be added if a full flux correction is 
 *             requested.
 *
 * @doc 
 *   Apply the atmospheric absorption and optionally the spectro-photometric
 *   correction to a sequence of at least one frame in 1D-stacked 
 *   format created by vmmosobsstare, vmmosobsjitter, or vmifuscience,
 *   from data obtained with the same instrument mode.
 *
 *   \begin{itemize}
 *
 *     \item ApplyResponse:     Apply also instrument response.
 *
 *   \end{itemize}
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 *
 * @author C.Izzo
 */   

static cxint
vmspphot(PilSetOfFrames *sof)
{
  const char  task[] = "vmspphot";

  const char  parameter[] = "Parameters";

  const char *mosReducedCategory = pilTrnGetCategory("MosScienceReduced");
  const char *ifuReducedCategory = pilTrnGetCategory("IfuScienceReduced");
  const char *reducedCategory;

  const char *mosFluxReducedCategory = 
                                   pilTrnGetCategory("MosScienceFluxReduced");
  const char *ifuFluxReducedCategory = 
                                   pilTrnGetCategory("IfuScienceFluxReduced");
  const char *fluxReducedCategory;

  const char *mosSphotTableCategory = pilTrnGetCategory("MosSphotTable");
  const char *ifuSphotTableCategory = pilTrnGetCategory("IfuSphotTable");
  const char *sphotTableCategory;

  const char *atmTableCategory = pilTrnGetCategory("ExtinctTable");


  char                   output1DName[PATHNAME_MAX + 1];

  VimosBool              updateOK = VM_TRUE;

  unsigned int           calibrateFlux;
  unsigned int           error;

  int                    i, j;

  int                    ifuCount, imaCount, mosCount;
  int                    firstQuadrant, quadrant;
  int                    firstGrism, grism;

  PilFrame              *frame;
  PilFrame              *sphotFrame;
  PilFrame              *atmFrame;
  PilFrame              *outputFrame;

  VimosImage            *sphotFile           = NULL;
  VimosImage            *atmFile             = NULL;
  VimosImage           **imaList             = NULL;
  VimosImage            *tmpImage;
  VimosTable            *sphotTable          = NULL;
  VimosTable            *atmTable            = NULL;



  /*
   * Check if the instrument response should also be applied.
   */

  calibrateFlux = pilDfsDbGetBool(parameter, "ApplyResponse", 1);


  /*
   * Load input frames
   */

  cpl_msg_info(task, "Loading input frames...");


  /*
   * Get 1D-extracted MOS or IFU frames
   */

  ifuCount = pilSofFrameCount(sof, ifuReducedCategory);
  mosCount = pilSofFrameCount(sof, mosReducedCategory);


  if (ifuCount && mosCount) {
    cpl_msg_error(task, "Just MOS or just IFU spectra are allowed in input: "
                  "not both.");
    return EXIT_FAILURE;
  }

  imaCount = ifuCount + mosCount;

  if (imaCount == 0) {
    if (pilSofFrameCount(sof, ifuFluxReducedCategory) 
        || pilSofFrameCount(sof, mosFluxReducedCategory))
      cpl_msg_error(task, "Input frames are already flux calibrated.");
    else
      cpl_msg_error(task, "No 1D-extracted spectral frames found in input.");
    return EXIT_FAILURE;
  }

  if (ifuCount) {
    sphotTableCategory = ifuSphotTableCategory;
    reducedCategory = ifuReducedCategory;
    fluxReducedCategory = ifuFluxReducedCategory;
  }
  else {
    sphotTableCategory = mosSphotTableCategory;
    reducedCategory = mosReducedCategory;
    fluxReducedCategory = mosFluxReducedCategory;
  }

  sphotFrame = pilSofLookup(sof, sphotTableCategory);
  if (sphotFrame)
    pilFrmSetType(sphotFrame, PIL_FRAME_TYPE_CALIB);

  if (calibrateFlux) {

    /*
     * Check that a SPH table is present in input frames.
     */

    if (!sphotFrame) {
      cpl_msg_error(task, "%s not found in input: the instrument "
                    "response correction cannot be applied.",
                    sphotTableCategory);
      return EXIT_FAILURE;
    }
  }


  /*
   * Check that an atmospheric extinction table is present in input frames.
   */

  atmFrame = pilSofLookup(sof, atmTableCategory);
  if (!atmFrame) {
    cpl_msg_error(task, "Missing atmospheric extinction table: "
                  "input frames cannot be corrected.");
    return EXIT_FAILURE;
  }
  else
    pilFrmSetType(atmFrame, PIL_FRAME_TYPE_CALIB);


  /*
   * Load images
   */

  if ((imaList = (VimosImage **)cpl_calloc(imaCount, sizeof(VimosImage *)))) {

    frame = pilSofLookupNext(sof, reducedCategory);

    for (i = 0; i < imaCount; i++) {
      if ((imaList[i] = openOldFitsFile(pilFrmGetName(frame), 1, 0))) {
        closeFitsImage(imaList[i], 0);
        pilFrmSetType(frame, PIL_FRAME_TYPE_CALIB);
      }
      else {
        cpl_msg_error(task, "Failure opening %s image %d", reducedCategory,
                      i + 1);
        for (j = 0; j < i; j++)
          deleteImage(imaList[j]);
        cpl_free(imaList);
        return EXIT_FAILURE;
      }
      frame = pilSofLookupNext(sof, NULL);
    }
  }
  else {
    cpl_msg_error(task, "Failure creating list of input images");
    return EXIT_FAILURE;
  }


  for (i = 0; i < imaCount; i++) {
    if (readIntDescriptor(imaList[i]->descs, pilTrnGetKeyword("Quadrant"), 
                          &quadrant, NULL) == VM_TRUE) {
      if (i) {
        if (firstQuadrant == quadrant) {
          continue;
        }
        cpl_msg_error(task, "Input images do not come from the same quadrant");
      }
      else {
        firstQuadrant = quadrant;
        continue;
      }
    }
    else
      cpl_msg_error(task, "Missing descriptor %s in input images", 
                    pilTrnGetKeyword("Quadrant"));
    for (i = 0; i < imaCount; i++)
      deleteImage(imaList[i]);
    cpl_free(imaList);
    return EXIT_FAILURE;
  }


  for (i = 0; i < imaCount; i++) {
    grism = getGrism(imaList[i]);
    if (grism > -1) {
      if (i) {
        if (firstGrism == grism) {
          continue;
        }
        cpl_msg_error(task, "Input images do not come from the same grism");
      }
      else {
        firstGrism = grism;
        continue;
      }
    }
    else
      cpl_msg_error(task, "Missing grism name in input images.");
    for (i = 0; i < imaCount; i++)
      deleteImage(imaList[i]);
    cpl_free(imaList);
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
    for (j = 0; j < imaCount; j++)
      deleteImage(imaList[j]);
    cpl_free(imaList);
    return EXIT_FAILURE;
  }


  /* 
   * Spectrophotometric calibration to 1D spectra
   */ 

  if (calibrateFlux) {
  
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
      deleteTable(atmTable);
      for (j = 0; j < imaCount; j++)
        deleteImage(imaList[j]);
      cpl_free(imaList);
      return EXIT_FAILURE;
    }

    if (readIntDescriptor(sphotTable->descs, pilTrnGetKeyword("Quadrant"),
                          &quadrant, NULL) == VM_TRUE) {
      if (firstQuadrant != quadrant) {
        cpl_msg_error(task, "The input spectro-photometric table is related "
                      "to quadrant %d instead of %d", quadrant, firstQuadrant);
        deleteTable(atmTable);
        deleteTable(sphotTable);
        for (j = 0; j < imaCount; j++)
          deleteImage(imaList[j]);
        cpl_free(imaList);
        return EXIT_FAILURE;
      }
    }
  
  
    grism = getGrismAgain(sphotTable);
  
    if (grism > -1) {
      if (firstGrism != grism) {
        cpl_msg_error(task, "The input spectro-photometric table is related "
                      "to a different grism ithan input images");
        deleteTable(atmTable);
        deleteTable(sphotTable);
        for (j = 0; j < imaCount; j++)
          deleteImage(imaList[j]);
        cpl_free(imaList);
        return EXIT_FAILURE;
      }
    }
  }


  cpl_msg_info(task, "Applying atmospheric extinction correction...");

  if (calibrateFlux)
    cpl_msg_info(task, "Applying instrument response correction...");

  for (i = 0; i < imaCount; i++) {
    tmpImage = VmSpApplyPhot(imaList[i], sphotTable, atmTable);
    if (!tmpImage) {
      cpl_msg_error(task, "Spectro-photometric calibration failure.");
      deleteTable(atmTable);
      deleteTable(sphotTable);
      for (j = 0; j < imaCount; j++)
        deleteImage(imaList[j]);
      cpl_free(imaList);
      return EXIT_FAILURE;
    }
    deleteImage(imaList[i]);
    imaList[i] = tmpImage;
  }

  deleteTable(atmTable);
  deleteTable(sphotTable);

  for (i = 0; i < imaCount; i++) {

    if (calibrateFlux) {
      updateOK = updateOK && insertDoubleDescriptor(&(imaList[i]->descs),
                             pilTrnGetKeyword("ExposureTime"),
                             1.0,
                             pilTrnGetComment("ExposureTime"),
                             "ESO*", 1);
    }

    updateOK = updateOK && insertDoubleDescriptor(&(imaList[i]->descs),
                           pilTrnGetKeyword("DataMin"),
                           imageMinimum(imaList[i]),
                           pilTrnGetComment("DataMin"),
                           "ESO*", 1);

    updateOK = updateOK && insertDoubleDescriptor(&(imaList[i]->descs),
                           pilTrnGetKeyword("DataMax"),
                           imageMaximum(imaList[i]),
                           pilTrnGetComment("DataMax"),
                           "ESO*", 1);

    updateOK = updateOK && writeDoubleDescriptor(&(imaList[i]->descs),
                           pilTrnGetKeyword("DataMean"),
                           imageMean(imaList[i]),
                           pilTrnGetComment("DataMean"));

    updateOK = updateOK && writeDoubleDescriptor(&(imaList[i]->descs),
                           pilTrnGetKeyword("DataStdDeviation"),
                           imageSigma(imaList[i]),
                           pilTrnGetComment("DataStdDeviation"));

    updateOK = updateOK && writeDoubleDescriptor(&(imaList[i]->descs),
                           pilTrnGetKeyword("DataMedian"),
                           imageMedian(imaList[i]),
                           pilTrnGetComment("DataMedian"));
  
    updateOK = updateOK && writeDoubleDescriptor(&(imaList[i]->descs),
                           pilTrnGetKeyword("AirMass"), 
                           0.0,
                           pilTrnGetComment("AirMass"));

    if (!updateOK) {
      cpl_msg_error(task, "Failure updating product %d header", i+1);
      for (j = 0; j < imaCount; j++)
        deleteImage(imaList[j]);
      cpl_free(imaList);
      return EXIT_FAILURE;
    } 

    if (calibrateFlux)
      reducedCategory = fluxReducedCategory;

    sprintf(output1DName, "%s_%d.fits", reducedCategory, i);
    vmstrlower(output1DName);

    if (createFitsImage(output1DName, imaList[i], reducedCategory)) {
      outputFrame = newPilFrame(output1DName, reducedCategory);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
    }
    else {
      cpl_msg_error(task, "Cannot create local product file %s!", output1DName);
      for (j = 0; j < imaCount; j++)
        deleteImage(imaList[j]);
      cpl_free(imaList);
      return EXIT_FAILURE;
    }

  }

  for (j = 0; j < imaCount; j++)
    deleteImage(imaList[j]);
  cpl_free(imaList);

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
                    "vmspphot",
    "Apply to IFU or MOS 1D-extracted spectra a flux correction.",
    "This recipe can be used to apply a spectro-photometric flux calibration\n"
    "to any number of 1D-extracted spectral frames generated by the recipes\n"
    "vmmosobsstare, vmmosobsjitter, and vmifuscience.\n\n"
"Input files:\n\n"
"  DO category:                Type:       Explanation:            Required:\n"
"  MOS_SCIENCE_REDUCED \n"
"  or IFU_SCIENCE_REDUCED      Product     Extracted MOS or IFU spectra Y\n"
"  EXTINCT_TABLE               Calib       Atmospheric extinction       Y\n"
"  MOS_SPECPHOT_TABLE\n"
"  or IFU_SPECPHOT_TABLE       Calib       Response curve               .\n\n"
"Output files:\n\n"
"  DO category:                Data type:  Explanation:\n"
"  MOS_SCIENCE_FLUX_REDUCED\n"
"  or IFU_SCIENCE_FLUX_REDUCED FITS image  Flux calibrated objects spectra\n\n"
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

                    vmspphot_create,
                    vmspphot_exec,
                    vmspphot_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
