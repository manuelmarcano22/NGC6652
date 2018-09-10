/* $Id: vmifucombinecube.cc,v 1.7 2013-08-23 10:23:38 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2009 European Southern Observatory
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
 * $Date: 2013-08-23 10:23:38 $
 * $Revision: 1.7 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/

#include <string.h>
#include <math.h>
#include <cmath>

#include <cxmemory.h>

#include <cpl_recipe.h>
#include <cpl_plugininfo.h>
#include <cpl_parameterlist.h>
#include <cpl_frameset.h>
#include <cpl_propertylist.h>
#include <cpl_imagelist.h>
#include <cpl_wcs.h>

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
#include "vimos_ifu_wcs.h"
#include "vimos_utils.h"
#include "vimos_dfs.h"

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
cxint cpl_plugin_get_info(cpl_pluginlist *list);

static int vmifucombinecube_create(cpl_plugin *) ;
static int vmifucombinecube_exec(cpl_plugin *) ;
static int vmifucombinecube_destroy(cpl_plugin *) ;
static cxint vmifucombinecube(PilSetOfFrames *);

cpl_propertylist * vimos_ifu_combinecube_add_wcs
(PilSetOfFrames * sof, const char *inputTag, int nx_fiber, int ny_fiber, 
 double spec_crpix1, double spec_crval1, double spec_cdelt1);

/*-----------------------------------------------------------------------------
                            Structure definitions
 -----------------------------------------------------------------------------*/

static char vmifucombinecube_description[] =
"This recipe is used to rearrange into a single cube the images of\n" 
"extracted spectra from different VIMOS quadrants. Such images are\n"
"created by the recipes vmifucalib, vmifuscience, and vmifustandard,\n"
"and the classification tag of the product depends on the classification\n"
"tag of the inputs:\n\n"
"         Input                         Output\n\n"
"   IFU_ARC_SPECTRUM_EXTRACTED    IFU_ARC_SPECTRUM_CUBE\n"
"   IFU_FLAT_SPECTRUM_EXTRACTED   IFU_FLAT_SPECTRUM_CUBE\n"
"   IFU_SCIENCE_REDUCED           IFU_SCIENCE_CUBE\n"
"   IFU_SCIENCE_FLUX_REDUCED      IFU_SCIENCE_FLUX_CUBE\n"
"   IFU_STANDARD_REDUCED          IFU_STANDARD_CUBE\n\n"
"The input images must belong to different quadrants, so that they\n"
"cannot be more than 4. The allocated cubes are the smallest possible,\n"
"depending on the number of quadrants involved and on whether the IFU\n"
"shutter was on or off. The smallest cubes (20x20xN) are produced when\n" 
"just one quadrant is input and the shutter was on. The largest cubes\n"
"(80x80xN) are produced when 3 or more quadrants are given in input\n"
"(and when either quadrants 1 and 3, or 2 and 4, are input), and the\n"
"shutter was off. Pixels belonging to missing quadrants are padded\n" 
"with zeroes if necessary.\n\n";

/*-----------------------------------------------------------------------------
                                Functions code
 -----------------------------------------------------------------------------*/

/*
 * Build table of contents, i.e. the list of available plugins, for
 * this module. This function is exported.
 */

int
cpl_plugin_get_info(cpl_pluginlist *list)
{

    cpl_recipe *recipe = (cpl_recipe*)cx_calloc(1, sizeof *recipe);
    cpl_plugin *plugin = &recipe->interface;


    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmifucombinecube",
                    "Rearrange into cube format images of extracted spectra"
                    "produced by vmifucalib, vmifuscience, and vmifustandard.",
                    vmifucombinecube_description,
                    "ESO VIMOS Pipeline Team and VIMOS Consortium",
                    PACKAGE_BUGREPORT,
                    vimos_get_license(),
                    vmifucombinecube_create,
                    vmifucombinecube_exec,
                    vmifucombinecube_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}

/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmifucombinecube_create(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;

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
     * No parameter list
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
vmifucombinecube_exec(cpl_plugin *plugin)
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

    try
    {
        if (vmifucombinecube(sof) == EXIT_SUCCESS) {
         
            const char * recipe_name = "vmifucombinecube";
       
            /*
             * Post process the product frames registered in the set
             * of frames.
             */

            status = vmCplPostProcessFrames(sof, recipe_name);
          
            if (status == 0) {
  
                /*
                 * Update recipe interface with the product frames found in sof
                 * and destroy it.
                 */

                status = vmCplFramesetImport(recipe->frames, sof);
            }
      }
      else
          status = 1;
    }
    catch (...)
    {
        status = 1;
        cpl_msg_error(cpl_func, "An uncaught error during recipe execution");
    }

    vmCplRecipeTimerStop(NULL);


    /*
     * Release locally acquired resources
     */

    deletePilSetOfFrames(sof);
    
    return status == 0 ? 0 : 1;

}


static cxint
vmifucombinecube_destroy(cpl_plugin *plugin)
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
 *   Rearrange images with extracted IFU spectra into cubes.
 * 
 * @return EXIT_SUCCESS or EXIT_FAILURE.   
 *
 * @param sof  Set of frames containing up to 4 images with extracted
 *             spectra from different quadrants, generated by the 
 *             vmifuscience and vmifucalib recipes.
 *  
 * @doc
 *   This recipe would simply rearrange the extracted spectra
 *   contained in the input images into a cube. The input images
 *   must belong to different VIMOS quadrants (and therefore not 
 *   more than 4 images may be specified in input). Any input image
 *   must have size Nx400 (HR/MR grisms) or Nx1600 (LR grisms). 
 *   If the inputs have not the same value for the TPL.START 
 *   keyword, a warning is issued but the images are combined 
 *   anyway. 
 *
 * @author C. Izzo
 */

static cxint 
vmifucombinecube(PilSetOfFrames *sof)
{

  const char     task[] = "vmifucombinecube";

  const char    *cubeTag;
  const char    *inputTag = NULL;
  const char    *tag;
  const char    *warning;
  char           cubeName[PATHNAME_MAX + 1];
  int            inputCount;

  int            firstSlit, lastSlit, slit;
  int            xlen, ylen;
  int            xref, yref;
  int            xlow, ylow, xhig, yhig;
  int            xl[4], yl[4], xh[4], yh[4]; 
  int            i, j, k, q;
  double         spec_crpix1, spec_crval1, spec_cdelt1;

  cpl_image     *extracted;
  cpl_image    **images;
  cpl_imagelist *planes;
  PilFrame      *frame;
  PilFrame      *outputFrame;

  cpl_propertylist *header;

  char           tplStart[24];
  const char    *astart;
  int            quadrant[] = {0, 0, 0, 0};
  double         refnorm = 1.0;
  double         norm = 1.0;
  double        *integrals;
  double        *p;
  float         *data;


  /*
   * Make sure that the input SOF include the right datasets.
   * All inputs must have the same tag, and only few legal tags
   * are accepted.
   */

  if (pilSofIsEmpty(sof)) {
    cpl_msg_error(task, "No frames in input...");
    return EXIT_FAILURE;
  }

  inputCount = 0;
  frame = pilSofFirst(sof);

  while (frame) {
    if (inputTag == NULL) {
      inputTag = pilFrmGetCategory(frame);
      if (inputTag == NULL) {
        cpl_msg_error(task, "Missing tags in input frames...");
        return EXIT_FAILURE;
      }
    }
    else {
      tag = pilFrmGetCategory(frame);
      if (tag == NULL) {
        cpl_msg_error(task, "Missing tags in input frames...");
        return EXIT_FAILURE;
      }
      if (strcmp(inputTag, tag)) {
        cpl_msg_error(task, "All input frames must be of the same type...");
        return EXIT_FAILURE;
      }
    }
    inputCount++;
    frame = pilSofNext(sof, frame);
  }

  cpl_msg_info(task, "Found %d input %s", inputCount, inputTag);

  if (strcmp(inputTag,        "IFU_ARC_SPECTRUM_EXTRACTED"   ) == 0) {
    cubeTag =                 "IFU_ARC_SPECTRUM_CUBE";
  }
  else if (strcmp(inputTag,   "IFU_FLAT_SPECTRUM_EXTRACTED"  ) == 0) {
    cubeTag =                 "IFU_FLAT_SPECTRUM_CUBE";
  }
  else if (strcmp(inputTag,   "IFU_SCIENCE_REDUCED"          ) == 0) {
    cubeTag =                 "IFU_SCIENCE_CUBE";
  }
  else if (strcmp(inputTag,   "IFU_SCIENCE_FLUX_REDUCED"     ) == 0) {
    cubeTag =                 "IFU_SCIENCE_FLUX_CUBE";
  }
  else if (strcmp(inputTag,   "IFU_STANDARD_REDUCED"         ) == 0) {
    cubeTag =                 "IFU_STANDARD_CUBE";
  }
  else {
    cpl_msg_error(task, "Illegal input: %s. Only extracted IFU spectra "
                  "can be arraged into a cube", inputTag);
    return EXIT_FAILURE;
  }

  if (inputCount > 4) {
    cpl_msg_error(task, "Too many frames in input: up to 4 may be given.");
    return EXIT_FAILURE;
  }

  frame = pilSofLookupNext(sof, inputTag);

  header = cpl_propertylist_load_regexp(pilFrmGetName(frame), 0, 
                                        "C(RVAL|RPIX|DELT|TYPE)[0-9]+", 0);
  spec_crpix1 = cpl_propertylist_get_double(header, "CRPIX1");
  spec_crval1 = cpl_propertylist_get_double(header, "CRVAL1");
  spec_cdelt1 = cpl_propertylist_get_double(header, "CDELT1");
  cpl_propertylist_delete(header);

  for (i = 0; i < inputCount; i++) {

    pilFrmSetType(frame, PIL_FRAME_TYPE_CALIB);

    header = cpl_propertylist_load_regexp(pilFrmGetName(frame), 0, 
                                          pilTrnGetKeyword("Quadrant"), 0);
    q = cpl_propertylist_get_int(header, pilTrnGetKeyword("Quadrant")) - 1;

    cpl_propertylist_delete(header);

    if (quadrant[q]) {
      cpl_msg_error(task, "Input frames belonging to the same quadrant.");
      return EXIT_FAILURE;
    }

    quadrant[q] = 1;

    header = cpl_propertylist_load_regexp(pilFrmGetName(frame), 0, 
                                          "ESO TPL START", 0);
    astart = cpl_propertylist_get_string(header, "ESO TPL START");

    cpl_propertylist_delete(header);

    if (astart == NULL) {
      cpl_msg_error(task, "Missing TPL START keyword in file %s",
                    pilFrmGetName(frame));
      return EXIT_FAILURE;
    }

    if (i) {
      if (strcmp(tplStart, astart)) {
        if (warning) {
          cpl_msg_warning(task, "%s", warning);
          warning = NULL;                       /* Print just once */
        }
      }
    }
    else {
      strncpy(tplStart, astart, 24);
      warning = "Input frames are from different exposures.";
    }

    header = cpl_propertylist_load_regexp(pilFrmGetName(frame), 0, "NAXIS1", 0);
    xlen = cpl_propertylist_get_int(header, "NAXIS1");
    cpl_propertylist_delete(header);
    header = cpl_propertylist_load_regexp(pilFrmGetName(frame), 0, "NAXIS2", 0);
    ylen = cpl_propertylist_get_int(header, "NAXIS2");
    cpl_propertylist_delete(header);

    if (xlen == 0) {
      cpl_msg_error(task, "Missing NAXIS1 keyword in file %s",
                    pilFrmGetName(frame));
      return EXIT_FAILURE;
    }

    if (ylen == 0) {
      cpl_msg_error(task, "Missing NAXIS2 keyword in file %s",
                    pilFrmGetName(frame));
      return EXIT_FAILURE;
    }

    if (i) {
      if (xlen != xref || ylen != yref) {
        cpl_msg_error(task, "Input frames must all have the same size.");
        return EXIT_FAILURE;
      }
    }
    else {
      if (ylen != 400 && ylen != 1600) {
        cpl_msg_error(task, "Input frames have unexpected Y size "
                      "(it should be either 400 or 1600 pixels)");
        return EXIT_FAILURE;
      }
      xref = xlen;
      yref = ylen;
    }

    frame = pilSofLookupNext(sof, NULL);
  }


  /*
   * Determine the true size of a cube spatial slice: images are
   * originally 80x80 (as the full IFU head), but if the shutter was
   * present they may be reduced to 40x40. Also, if in input not
   * all quadrants are present, the size of the cube may be further
   * reduced. Here the span of the used region is defined:
   */ 

  if (ylen == 400) {
    firstSlit = 1;
    lastSlit = 1;
    xl[0] = 40; yl[0] = 40; xh[0] = 60; yh[0] = 60;
    xl[1] = 20; yl[1] = 40; xh[1] = 40; yh[1] = 60;
    xl[2] = 20; yl[2] = 20; xh[2] = 40; yh[2] = 40;
    xl[3] = 40; yl[3] = 20; xh[3] = 60; yh[3] = 40;
  }

  if (ylen == 1600) {
    firstSlit = 0;
    lastSlit = 3;
    xl[0] = 40; yl[0] = 40; xh[0] = 80; yh[0] = 80;
    xl[1] =  0; yl[1] = 40; xh[1] = 40; yh[1] = 80;
    xl[2] =  0; yl[2] =  0; xh[2] = 40; yh[2] = 40;
    xl[3] = 40; yl[3] =  0; xh[3] = 80; yh[3] = 40;
  }

  xlow = ylow = xhig = yhig = 40;

  for (q = 0; q < 4; q++) {
    if (quadrant[q]) {
      if (xlow > xl[q]) 
        xlow = xl[q];
      if (ylow > yl[q]) 
        ylow = yl[q];
      if (xhig < xh[q]) 
        xhig = xh[q];
      if (yhig < yh[q]) 
        yhig = yh[q];
    }
  }

  cpl_msg_info(task, "Reconstructing cube %d x %d x %d "
               "(IFU head from (%d,%d) to (%d,%d))", 
               xhig - xlow, yhig - ylow, xlen, xlow, ylow, xhig, yhig);

  cpl_msg_info(task, "Loading input frames...");


  /*
   * We use an array of images because the imagelist is too rigid
   * (it doesn't allow images of different sizes).
   */

  images = (cpl_image**)cpl_calloc(xlen, sizeof(cpl_image *));
  integrals = (double *)cpl_calloc(ylen, sizeof(double));
  norm = 1.0;

  frame = pilSofLookupNext(sof, inputTag);

  for (i = 0; i < inputCount; i++) {

    extracted = cpl_image_load(pilFrmGetName(frame), CPL_TYPE_FLOAT, 0, 0);

    if (extracted) {

      header = cpl_propertylist_load_regexp(pilFrmGetName(frame), 0,
                                            pilTrnGetKeyword("Quadrant"), 0);
      q = cpl_propertylist_get_int(header, pilTrnGetKeyword("Quadrant"));
      cpl_propertylist_delete(header);

      header = cpl_propertylist_load_regexp(pilFrmGetName(frame), 0,
                                            "ESO QC IFU FLAT FLUX", 0);
      if (header && cpl_propertylist_has(header, "ESO QC IFU FLAT FLUX")) {
        if (i) {
          norm = cpl_propertylist_get_double(header, "ESO QC IFU FLAT FLUX");
          norm /= refnorm;
        }
        else {
          refnorm = cpl_propertylist_get_double(header, "ESO QC IFU FLAT FLUX");
        }
      }
      else {
        norm = 1.0;
      }

      cpl_propertylist_delete(header);

      if (i)
        cpl_image_divide_scalar(extracted, norm);

      for (j = 0; j < xlen; j++) { /* Plane by plane = wave by wave */

        data = (float *)cpl_image_get_data(extracted);

        for (k = 0; k < ylen; k++, data += xlen) /* Spec by spec = row by row */
          integrals[k] = data[j];

        if (i == 0)
          images[j] = cpl_image_new(80, 80, CPL_TYPE_FLOAT);

        p = integrals;
        for (slit = firstSlit; slit <= lastSlit; slit++, p += 400)
          ifuImage(images[j], p, q, slit);
      }

    }
    else {
      cpl_msg_error(task, "Failure loading %s frame %d", inputTag, i + 1);
      cpl_free(integrals);
      cpl_free(images);
      return EXIT_FAILURE;
    }

    cpl_image_delete(extracted);
    frame = pilSofLookupNext(sof, NULL);

  }

  cpl_free(integrals);

  /*
   * Now store result into the smallest possible cube (if necessary)
   */

  planes = cpl_imagelist_new();
  if (xhig - xlow != 80 || yhig - ylow != 80) {
    for (j = 0; j < xlen; j++) { /* Plane by plane = wave by wave */
      cpl_imagelist_set(planes, 
                        cpl_image_extract(images[j], 
                                          xlow + 1, ylow + 1, xhig, yhig), 
                        j);
      cpl_image_delete(images[j]);
    }
    cpl_free(images);
  }
  else {
    for (j = 0; j < xlen; j++) { /* Plane by plane = wave by wave */
      cpl_imagelist_set(planes, images[j], j);
    }
  }

  /* Compute the WCS keywords */
  
  cpl_propertylist * wcsheader = vimos_ifu_combinecube_add_wcs
          (sof, inputTag, xhig - xlow, yhig - ylow, 
           spec_crpix1, spec_crval1, spec_cdelt1);
  

  /*
   * Create the product file on disk, set the product attributes and
   * update the set of frames.
   */

  cpl_msg_info(cpl_func, "Saving the cube...");
  vmstrlower(strcpy(cubeName, cubeTag));
  strcat(cubeName, ".fits");


  if (CPL_ERROR_NONE == cpl_imagelist_save(planes, cubeName,
                                           CPL_BPP_IEEE_FLOAT, wcsheader,
                                           CPL_IO_DEFAULT)) {
    outputFrame = newPilFrame(cubeName, cubeTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", cubeName);
    cpl_propertylist_delete(wcsheader);
    cpl_imagelist_delete(planes);
    return EXIT_FAILURE;
  }

  cpl_propertylist_delete(wcsheader);
  cpl_imagelist_delete(planes);

  return EXIT_SUCCESS;
}

cpl_propertylist * vimos_ifu_combinecube_add_wcs
(PilSetOfFrames * sof, const char *inputTag, int nx_fiber, int ny_fiber, 
 double spec_crpix1, double spec_crval1, double spec_cdelt1)
{
    double alpha_target;
    double delta_target;
    double ra_vimos;
    double dec_vimos;
    double epoch;
    double equinox;
    double plate_scale;
    double pos_angle_ins;
    double pos_angle_wcs;
    cpl_propertylist * wcsheader;
    PilFrame *wcs_ref_frame = pilSofFirst((PilSetOfFrames *)sof);
    
    if((vimos_check_equal_keys(sof, inputTag, "ESO TEL TARG ALPHA", alpha_target) &&
        vimos_check_equal_keys(sof, inputTag, "ESO TEL TARG DELTA", delta_target) &&
        vimos_check_equal_keys(sof, inputTag, "RA", ra_vimos) &&
        vimos_check_equal_keys(sof, inputTag, "DEC", dec_vimos) &&
        vimos_check_equal_keys(sof, inputTag, "ESO TEL TARG EPOCH", epoch) &&
        vimos_check_equal_keys(sof, inputTag, "ESO TEL TARG EQUINOX", equinox) &&
        vimos_check_equal_keys(sof, inputTag, "ESO INS PIXSCALE", plate_scale) &&
        vimos_check_equal_keys(sof, inputTag, "ESO ADA POSANG", pos_angle_ins)))
    {
        cpl_msg_info(cpl_func,"Adding WCS information to the cube");

        cpl_propertylist * header;
        header = cpl_propertylist_load(pilFrmGetName(wcs_ref_frame),0);
        pos_angle_wcs = - (pos_angle_ins + 90);
        
        two_d_linear_wcs spatial_wcs = 
                vimos_ifu_get_2d_wcs_from_pointing(ra_vimos, dec_vimos, 
                pos_angle_wcs, nx_fiber+1, ny_fiber+1, 
                plate_scale, epoch, equinox);

        wcsheader = cpl_propertylist_new();
        cpl_propertylist_update_double(wcsheader, "CRPIX1", spatial_wcs.crpix1());
        cpl_propertylist_update_double(wcsheader, "CRPIX2", spatial_wcs.crpix2());
        cpl_propertylist_update_double(wcsheader, "CRPIX3", spec_crpix1);
        cpl_propertylist_update_double(wcsheader, "CRVAL1", spatial_wcs.crval1());
        cpl_propertylist_update_double(wcsheader, "CRVAL2", spatial_wcs.crval2());
        cpl_propertylist_update_double(wcsheader, "CRVAL3", spec_crval1);
        cpl_propertylist_update_double(wcsheader, "CDELT1", spatial_wcs.cdelt1());
        cpl_propertylist_update_double(wcsheader, "CDELT2", spatial_wcs.cdelt2());
        cpl_propertylist_update_double(wcsheader, "CDELT3", spec_cdelt1);
        cpl_propertylist_update_double(wcsheader, "CD1_1",  spatial_wcs.cd_matrix()[0]);
        cpl_propertylist_update_double(wcsheader, "CD1_2",  spatial_wcs.cd_matrix()[1]);
        cpl_propertylist_update_double(wcsheader, "CD2_1",  spatial_wcs.cd_matrix()[2]);
        cpl_propertylist_update_double(wcsheader, "CD2_2",  spatial_wcs.cd_matrix()[3]);
        cpl_propertylist_update_string(wcsheader, "CTYPE1", spatial_wcs.ctype1().c_str());
        cpl_propertylist_update_string(wcsheader, "CTYPE2", spatial_wcs.ctype2().c_str());
        cpl_propertylist_update_string(wcsheader, "CTYPE3", "WAVE");
        cpl_propertylist_update_string(wcsheader, "CUNIT1", spatial_wcs.cunit1().c_str());
        cpl_propertylist_update_string(wcsheader, "CUNIT2", spatial_wcs.cunit2().c_str());
        cpl_propertylist_update_string(wcsheader, "CUNIT3", "Angstrom");
        
    }
    else
    {
      cpl_msg_warning(cpl_func,"WCS is not available on original images,"
              " or are not compatible. Final image without WCS");
      wcsheader = cpl_propertylist_new();
      cpl_propertylist_update_double(wcsheader, "CRPIX1", 1.0);
      cpl_propertylist_update_double(wcsheader, "CRPIX2", 1.0);
      cpl_propertylist_update_double(wcsheader, "CRPIX3", spec_crpix1);
      cpl_propertylist_update_double(wcsheader, "CRVAL1", nx_fiber+1);
      cpl_propertylist_update_double(wcsheader, "CRVAL2", ny_fiber+1);
      cpl_propertylist_update_double(wcsheader, "CRVAL3", spec_crval1);
      cpl_propertylist_update_double(wcsheader, "CDELT1", 1.0);
      cpl_propertylist_update_double(wcsheader, "CDELT2", 1.0);
      cpl_propertylist_update_double(wcsheader, "CDELT3", spec_cdelt1);
      cpl_propertylist_update_string(wcsheader, "CTYPE1", "FIBRE");
      cpl_propertylist_update_string(wcsheader, "CTYPE2", "FIBRE");
      cpl_propertylist_update_string(wcsheader, "CTYPE3", "WAVE");
    }
 
    
    return wcsheader;
}
