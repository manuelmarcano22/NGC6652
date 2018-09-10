/* $Id: vmimcalphot.c,v 1.5 2012-11-08 17:56:57 cgarcia Exp $
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
 * $Revision: 1.5 $
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
#include <pilstrutils.h>
#include <pilutils.h>
#include <pilfits.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmstarmatchtable.h"
#include "vmphotometrictable.h"
#include "vmfit.h"
#include "vmutils.h"
#include "vmcpl.h"
#include "vimos_dfs.h"


#define COMMENT_LENGTH             (80)
#define MAX_STRING_KEYWORD_LENGTH  (80)
#define MAG_COLUMN_ROOT_NAME       "MAG"


static cxint vmimcalphot(PilSetOfFrames *);


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmimcalphot_create(cpl_plugin *plugin)
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

    p = cpl_parameter_new_value("vimos.Parameters.extinction",
                                CPL_TYPE_BOOL,
                                "Compute also extinction coefficient",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "Extinction");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "Extinction");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.colorterm",
                                CPL_TYPE_BOOL,
                                "Compute also color term",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ColorTerm");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ColorTerm");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.usecolorterm",
                                CPL_TYPE_BOOL,
                                "Use color term in zeropoint computation",
                                "vimos.Parameters",
                                FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "UseColorTerm");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "UseColorTerm");
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
vmimcalphot_exec(cpl_plugin *plugin)
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

    if (vmimcalphot(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmimcalphot");
        
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
vmimcalphot_destroy(cpl_plugin *plugin)
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
 *   Determine a photometric calibration for a given filter/CCD 
 *   combination from a set of star match tables, and a first 
 *   guess photometric calibration.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param sof List of input frames, including one or more star match 
 *            tables, and a standard photometric calibration table.
 *
 * @doc 
 *   A photometric calibration is derived from a set of input
 *   star match tables, and a first guess photometric calibration.
 *   The magnitude night zero point is computed by averaging the
 *   difference between each star catalog magnitude with its 
 *   measured magnitude corrected by the corresponding airmass, 
 *   or, optionally, it is derived together a new extinction
 *   coefficient, or a new color term, or both, modelling the
 *   magnitude differences as a function of airmasses and/or 
 *   stars color indexes.
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
 *
 *   \begin{itemize}
 *
 *     \item Extinction       Determination of the extinction coefficient. 
 *                            If less than 4 input star match tables are 
 *                            found, there are not enough different airmasses 
 *                            to be fitted and therefore this option is 
 *                            disabled. In this case just the magnitude 
 *                            night zero point and, if enabled, the color 
 *                            term, are computed.
 *     \item ColorTerm        Determination of the color term. If less 
 *                            than 4 input stars are found, there are not 
 *                            enough different star colors cannot be fitted 
 *                            and therefore this option is disabled. In this 
 *                            case just the computation of the magnitude 
 *                            night zero point and, if enabled, the extinction 
 *                            coefficient, are computed.
 *   \end{itemize}
 *
 * @author B.Garilli, C.Izzo, adapting from Marco Lolli's (Bologna 
 *   Astronomical Observatory) MIDAS procedures VmImZero.prg, 
 *   VmImExt.prg, and VmImCalPhot.prg
 */   

static cxint
vmimcalphot(PilSetOfFrames *sof)
{
  const char   task[]             = "vmimcalphot";

  const char   parameter[]        = "Parameters";

  char         outPhotTableName[PATHNAME_MAX + 1];
  double       magzero;
  double       airmass            = 1;
  double       colorTerm          = 0;
  double       extinction;
  char         filterName[MAX_STRING_KEYWORD_LENGTH];
  char         filterId[MAX_STRING_KEYWORD_LENGTH];
  char         checkFilterId[MAX_STRING_KEYWORD_LENGTH];
  char         colorSystem[MAX_STRING_KEYWORD_LENGTH];
  const char  *starMatchCat       = pilTrnGetCategory("StarMatchTable");
  const char  *photoCategory      = pilTrnGetCategory("PhotometricTable");
  const char  *photoCoeffCategory = pilTrnGetCategory("PhotometricCoeffTable");
  int          nStarMatchTables   = pilSofFrameCount(sof, starMatchCat);
  int          nValidStarMatchTables;
  int          firstValid;
  int          numStars, nRows;
  float       *fmagValues         = NULL;
  float        sigma, mean, dummy, dummy1;
  double      *extinctionValues   = NULL;
  double      *colorTermValues    = NULL;
  PilFrame    *photometricFrame   = pilSofLookup(sof, photoCategory);
  VimosImage  *photometricFile;
  VimosTable  *photometricTable   = NULL;
  PilFrame    *currentFrame       = pilSofLookup(sof, starMatchCat);
  VimosImage  *currentFile;
  VimosTable **starMatchTable;
  char       **starMatchName;
  PilFrame    *outputFrame;

  int          i, j, k, len, quad, mode, error;
  unsigned int computeExtinction, computeColorTerm, useColorTerm;
  char        *magColumnName[2], col [2];
  char        *magFilColumnName   = NULL;
  double      *magValues          = NULL;
  VimosColumn *nameColumn, *magIdColumn, *zeropColumn;
  VimosColumn *starIdColumn, *magFilColumn;
  VimosColumn *magBestColumn;
  VimosColumn *magColumn[2];
  VimosPixel  *colorValues        = NULL;
  double       rmszerop, rmscolor, rms;
  char         comment[COMMENT_LENGTH], name[80];

  unsigned int doSomethingAnyway = 0;

 /*
  * The variable "doSomethingAnyway" is just a flag to change 
  * the behavior of this task. If set to 1, and not enough input 
  * star match tables are specified in input for fulfilling the 
  * computation of a requested extincion coefficient, or not
  * enough stars are found to permit a computation of a requested 
  * color term, then the requests will be switched off accordingly
  * (up to the point that just a magnitude zero point will be
  * computed), and a photometric calibration produced anyway. 
  * If instead "doSomethingAnyway" is set to 0, then the task
  * stops in case the input is insufficient to fulfill a request,
  * and and error is returned.
  */
  

  if (nStarMatchTables == 0) {
    cpl_msg_error(task, "No input Star Match Tables");
    return EXIT_FAILURE;
  }


 /*
  * Get the standard photometric table
  */

  error = 1;

  if (photometricFrame) {
    pilFrmSetType(photometricFrame, PIL_FRAME_TYPE_CALIB);
    if ((photometricFile = openOldFitsFile(pilFrmGetName(photometricFrame),
                                           0, 0))) {
      if ((photometricTable = newPhotometricTable())) {
        if (readFitsPhotometricTable(photometricTable, photometricFile->fptr) 
                                                              == VM_TRUE) {
          error = 0;

          readIntDescriptor(photometricTable->descs, 
                            pilTrnGetKeyword("Quadrant"), &quad, comment);
          readDoubleDescriptor(photometricTable->descs, 
                               pilTrnGetKeyword("Extinction"), 
                               &extinction, comment);
          readDoubleDescriptor(photometricTable->descs, 
                               pilTrnGetKeyword("ColorTerm"), 
                               &colorTerm, comment);
          readStringDescriptor(photometricTable->descs,
                               pilTrnGetKeyword("Colour"),
                               colorSystem, comment);
          readStringDescriptor(photometricTable->descs, 
                               pilTrnGetKeyword("FilterId", quad), 
                               filterId, comment);
          deleteAllColumns(photometricTable) ;
          vimosDscErase(&photometricTable->descs, "^ESO .*");
        }
        else {
          cpl_msg_error(task, "Failure in reading the standard photometric "
                      "calibration table %s", pilFrmGetName(photometricFrame));
          deleteTable(photometricTable);
        }
      }
      else
        cpl_msg_error(task, "Failure in memory allocation");
    }
    else
      cpl_msg_error(task, 
            "Failure opening the standard photometric calibration table");
  }
  else
    cpl_msg_error(task, "No input standard photometric calibration table");

  if (error)
    return EXIT_FAILURE;

  starMatchTable = 
             (VimosTable **)cpl_malloc(nStarMatchTables * sizeof(VimosTable *));

  if (!starMatchTable) {
    cpl_msg_error(task, "Failure in memory allocation");
    return EXIT_FAILURE;
  }

  starMatchName = (char **)cpl_calloc(nStarMatchTables, sizeof(char *));
  currentFrame = pilSofLookupNext(sof, starMatchCat);

  numStars = 0;
  firstValid = -1;
  nValidStarMatchTables = nStarMatchTables;
  for (i = 0; i < nStarMatchTables; i++) {
    starMatchName[i] = (char *)pilFrmGetName(currentFrame);
    currentFile = openOldFitsFile(pilFrmGetName(currentFrame), 0, 0);
    pilFrmSetType(currentFrame, PIL_FRAME_TYPE_CALIB);
    if ((starMatchTable[i] = newStarMatchTableEmpty())) {
      if (readFitsStarMatchTable(starMatchTable[i], currentFile->fptr) 
                                                           == VM_TRUE) {
        closeFitsImage(currentFile, 0);
        if (starMatchTable[i]->numColumns == 0) {
          deleteTable(starMatchTable[i]);
          starMatchTable[i] = NULL;
          nValidStarMatchTables--;
          currentFrame = pilSofLookupNext(sof, NULL);
          continue;
        }
        else {
          if (firstValid < 0) {
            firstValid = i;
          }
        }
        if (readStringDescriptor(starMatchTable[i]->descs, 
                             pilTrnGetKeyword("FilterId", quad), 
                             checkFilterId, comment) == VM_TRUE) {
          if (!strncmp(checkFilterId, filterId, strlen(filterId))) {
            numStars += starMatchTable[i]->cols->len;
          }
          else {
            cpl_msg_error(task, "Wrong filter ID for star match table %s "
                        "(%s found, %s expected)", 
                        pilFrmGetName(currentFrame), checkFilterId, filterId);
            for (j = 0; j <= i; j++)
              deleteTable(starMatchTable[j]);
            cpl_free(starMatchTable);
            deleteTable(photometricTable);
            return EXIT_FAILURE;
          }
        }
        else {
          cpl_msg_error(task, "Wrong quadrant number for star match table %s "
                      "(%d expected)", pilFrmGetName(currentFrame), quad);
          for (j = 0; j <= i; j++)
            deleteTable(starMatchTable[j]);
          cpl_free(starMatchTable);
          deleteTable(photometricTable);
          return EXIT_FAILURE;
        }
      }
      else {
        cpl_msg_error(task, "Cannot read star match table %s", 
                    pilFrmGetName(currentFrame));
        for (j = 0; j <= i; j++)
          deleteTable(starMatchTable[j]);
        cpl_free(starMatchTable);
        deleteTable(photometricTable);
        return EXIT_FAILURE;
      }
    }
    else {
      cpl_msg_error(task, "Failure in memory allocation");
      for (j = 0; j < i; j++)
        deleteTable(starMatchTable[j]);
      cpl_free(starMatchTable);
      deleteTable(photometricTable);
      return EXIT_FAILURE;
    }
    currentFrame = pilSofLookupNext(sof, NULL);
  }

  if (nValidStarMatchTables == 0) {
    cpl_msg_error(task, "No valid input Star Match Tables");
    return EXIT_FAILURE;
  }

  if (numStars == 0) {
    cpl_msg_error(task, "Empty star match table");
    return EXIT_FAILURE;
  }


 /*
  * Input opened successfully, now get parameters
  */

  cpl_msg_info(task, "Determining the magnitude zero point");

 /*
  * See if also the extinction coefficient should be computed
  */

  error = 0;

  computeExtinction = pilDfsDbGetBool(parameter, "Extinction", 1);

  if (computeExtinction) {
    if (nValidStarMatchTables < 4) {
      if (doSomethingAnyway) {
        cpl_msg_warning(task, "Less than 4 images at different airmasses; "
            "the determination of the extinction coefficient is disabled");
        computeExtinction = 0;
      }
      else {
        cpl_msg_error(task, "Less than 4 images at different airmasses; "
                    "the requested determination of the extinction "
                    "coefficient is impossible!");
        error = 1;
      }
    }
    else {
      cpl_msg_info(task, "together with the extinction coefficient");
    }
  }


 /*
  * See if also the color term should be computed
  */

  computeColorTerm = pilDfsDbGetBool(parameter, "ColorTerm", 1);

  if (computeColorTerm) {
    if (numStars < 4) {
      if (doSomethingAnyway) {
        cpl_msg_warning(task, "Less than 4 stars with different color indexes; "
          "the determination of the color term is disabled");
        computeColorTerm = 0;
      }
      else {
        cpl_msg_error(task, "Less than 4 stars with different color indexes; "
                    "the requested determination of the color term is "
                    "impossible!");
        error = 1;
      }
    }
    else {
      cpl_msg_info(task, "together with the color term");
    }
  }

  if (error) {
    for (j = 0; j < i; j++)
      deleteTable(starMatchTable[j]);
    cpl_free(starMatchTable);
    deleteTable(photometricTable);
    return EXIT_FAILURE;
  }

 /*
  * Compute the recipe running mode: the magnitude zero point is always
  * determined; in addition to that also the extinction coefficient 
  * and/or the color term may be determined. The calibration "mode" 
  * is coded as described in the following table:
  *
  *  Compute      Compute      Calibration
  *  extinction:  color term:  mode:
  *
  *     false       false       0
  *     true        false       1
  *     true        true        2
  *     false       true        3
  */

  mode = computeColorTerm * (3 - 2 * computeExtinction) + computeExtinction;


  /* 
   *  Add some more keywords to the photometric table header
   */

  readStringDescriptor(starMatchTable[firstValid]->descs, 
                       pilTrnGetKeyword("FilterName", quad), 
                       filterName, comment);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("MjdObs"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("DateObs"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("Instrument"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("OBS.DID"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("OBS.ID"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("PROG.ID"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("INS.DID"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("InstrumentMode"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs,
               pilTrnGetKeyword("FilterName", quad), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs,
               pilTrnGetKeyword("FilterId", quad), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("DET.DID"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs,
               pilTrnGetKeyword("Adu2Electron", 1), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs,
               pilTrnGetKeyword("Electron2Adu", 1), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("WINi.BINX", 1), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("WINi.BINY", 1), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("SeqWindowSizeX", 1), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("SeqWindowSizeY", 1), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("READ.MODE"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("READ.SPEED"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("READ.CLOCK"), NULL);

  vimosDscCopy(&photometricTable->descs, starMatchTable[firstValid]->descs, 
               pilTrnGetKeyword("Quadrant"), NULL);

  writeStringDescriptor(&photometricTable->descs, pilTrnGetKeyword("Table"), 
                        VM_IPC, "Photometric Table");

  /* 
   * Get build up column name from filter name.
   */

  len = strlen(MAG_COLUMN_ROOT_NAME);

  magFilColumnName = (char *) cpl_malloc((len+strlen(filterName)+2)*sizeof(char));

  photometricTable->numColumns = 3;
  photometricTable->cols = newStringColumn(numStars, "ImageName");

  /*
   * Now we are back to the original design of the VimosTable where
   * there was no "depth" so newStringColumn does not allocate
   * any more the string arrays. So I have changed the original
   * photometricTable->cols = newStringColumn(numStars,20,"ImageName")
   * (Paola)
   */

  nameColumn = photometricTable->cols;
  photometricTable->cols->next = newStringColumn(numStars, "STAR_ID");
  starIdColumn = photometricTable->cols->next;
  photometricTable->cols->next->next = newFloatColumn(numStars, "zeropoint");
  zeropColumn = photometricTable->cols->next->next ; 


  /*  
   * Determination of color system from 1st starMatchTable, they should 
   * all be the same...
   */

  sprintf(magFilColumnName, "%s_%s", MAG_COLUMN_ROOT_NAME, filterName);

  if (mode != 1) {
    useColorTerm = 1;
    if (mode == 0)
      useColorTerm = pilDfsDbGetBool(parameter, "UseColorTerm", 0);

    if (useColorTerm) {
      col[0] = colorSystem[0];
      col[1] = colorSystem[(strlen(colorSystem)-1)];
      for (i = 0; i < 2; i++) {
        magColumnName[i] = (char *)cpl_malloc((len + 2) * sizeof(char));
        sprintf(magColumnName[i], "%s_%c", MAG_COLUMN_ROOT_NAME, col[i]);
      }
      for (j = 0; j < 2; j++) {
        if (!(magColumn[j] = findColumn(starMatchTable[firstValid]->cols, 
                                        magColumnName[j]))) {
          cpl_msg_error(task, "No Info for %s in starmatch table", colorSystem);
          return EXIT_FAILURE;
        }
      }
    }
  } 
  else {
    strcpy(colorSystem, "none") ;
  }


 /*
  *  Allocate enough space for what is needed, according to "mode":
  */

  switch (mode) {
  case 0 : 
    fmagValues = (float *)cpl_malloc(numStars * sizeof(float)); 
    break;
  case 1 : 
    extinctionValues = (double *)cpl_malloc(numStars * sizeof(double));
    magValues = (double *)cpl_malloc(numStars * sizeof(double)); 
    break;
  case 2 : 
    colorValues = newPixel(numStars);
    break;
  case 3 :
    colorTermValues = (double *)cpl_malloc(numStars * sizeof(double));
    magValues = (double *)cpl_malloc(numStars * sizeof(double)); 
    break;
  default: 
    cpl_msg_debug(task, "Unrecognized mode");
    return EXIT_FAILURE;
  }

  k = 0;
  for (i = 0; i < nStarMatchTables; i++) {

    if (starMatchTable[i] == NULL)
      continue;

   /*  
    *  Airmass of current Star Match Table  
    */

    if (readDoubleDescriptor(starMatchTable[i]->descs, 
                             pilTrnGetKeyword("AirMass"), 
                             &airmass, comment) == VM_FALSE) {
      cpl_msg_error(task, "Cannot read descriptor %s", 
                  pilTrnGetKeyword("AirMass"));
      return EXIT_FAILURE;
    }


   /*  
    * Mag Zero Point of current Star Match Table 
    */

    if (readDoubleDescriptor(starMatchTable[i]->descs,
                             pilTrnGetKeyword("MagZero"), 
                             &magzero, comment) == VM_FALSE) {
      cpl_msg_error(task, "Cannot read descriptor %s", 
                  pilTrnGetKeyword("MagZero"));
      return EXIT_FAILURE;
    }


    /* 
     * File image name to be put in photom table If it doesn't exist, 
     * just put a generic "starMatchTable" string
     */

    if (readStringDescriptor(starMatchTable[i]->descs, 
                             pilTrnGetKeyword("ProductFile"),
                             name, comment) == VM_FALSE) {
      for (j = strlen(starMatchName[i]); j >= 0; j--)
        if (starMatchName[i][j] == '/')
          break;
      j++;
      strncpy(name, starMatchName[i] + j, 80);
    }

    nRows = starMatchTable[i]->cols->len;
    magFilColumn = findColumn(starMatchTable[i]->cols, magFilColumnName);
    magBestColumn = findColumn(starMatchTable[i]->cols, "MAG");
    /* magIdColumn = tblCopyColumn(starMatchTable[i], "ID");  ???? */
    magIdColumn = findColumn(starMatchTable[i]->cols, "ID");


    /* 
     * B.G.: mode 0 & 1 OK, 3 is formally OK, 2 seems also OK but 
     * for output (rms)
     */

    switch (mode) {
    case 0: 
      if (useColorTerm) {
        for (j = 0; j < nRows; j++) {
          if (magFilColumn->colValue->dArray[j] < 50 &&
              magColumn[0]->colValue->dArray[j] < 50 &&
              magColumn[1]->colValue->dArray[j] < 50) {
            fmagValues[k] = magzero
                          + extinction * airmass 
                          + colorTerm * (magColumn[0]->colValue->dArray[j]
                                       - magColumn[1]->colValue->dArray[j])
                          + magFilColumn->colValue->dArray[j]
                          - magBestColumn->colValue->dArray[j];
    
            zeropColumn->colValue->fArray[k] = fmagValues[k];
            starIdColumn->colValue->sArray[k] = 
                   magIdColumn->colValue->sArray[j];
            nameColumn->colValue->sArray[k] = cpl_strdup(name);
            k++;
          }
        }
      }
      else {
        for (j = 0; j < nRows; j++) {
          if (magFilColumn->colValue->dArray[j] < 50) {
            fmagValues[k] = magzero
                          + extinction * airmass
                          + magFilColumn->colValue->dArray[j]
                          - magBestColumn->colValue->dArray[j];

            zeropColumn->colValue->fArray[k] = fmagValues[k];
            starIdColumn->colValue->sArray[k] =
                   magIdColumn->colValue->sArray[j];
            nameColumn->colValue->sArray[k] = cpl_strdup(name);
            k++;
          }
        }
      }
      break;
    case 1: 
      for (j = 0; j < nRows; j++) {

        if (magFilColumn->colValue->dArray[j] < 50) {

          extinctionValues[k] = airmass;

          magValues[k] = magzero
                       + magFilColumn->colValue->dArray[j]
                       - magBestColumn->colValue->dArray[j];
  
          nameColumn->colValue->sArray[k] = cpl_strdup(name);
          zeropColumn->colValue->fArray[k] = (float)magValues[k];
          starIdColumn->colValue->sArray[k] = 
                 magIdColumn->colValue->sArray[j];
          k++;
        }
      }
      break;
    case 2: 
      for (j = 0; j < nRows; j++) {

        if (magFilColumn->colValue->dArray[j] < 50 &&
            magColumn[0]->colValue->dArray[j] < 50 &&
            magColumn[1]->colValue->dArray[j] < 50) {

          colorValues[k].x = airmass;

          colorValues[k].y = magColumn[0]->colValue->dArray[j]
                           - magColumn[1]->colValue->dArray[j];

          colorValues[k].i = magzero
                           + magFilColumn->colValue->dArray[j]
                           - magBestColumn->colValue->dArray[j];

          nameColumn->colValue->sArray[k] =  cpl_strdup(name);
          zeropColumn->colValue->fArray[k] = (float)colorValues[k].i;
          starIdColumn->colValue->sArray[k] =
                 magIdColumn->colValue->sArray[j];
          k++;
        }
      }
      break;
    case 3:
      for (j = 0; j < nRows; j++) {

        if (magFilColumn->colValue->dArray[j] < 50 &&
            magColumn[0]->colValue->dArray[j] < 50 &&
            magColumn[1]->colValue->dArray[j] < 50) {

          colorTermValues[k] = magColumn[0]->colValue->dArray[j]
                             - magColumn[1]->colValue->dArray[j];

          magValues[k] = magzero
                       + magFilColumn->colValue->dArray[j]
                       - magBestColumn->colValue->dArray[j];

          nameColumn->colValue->sArray[k]= cpl_strdup(name);
          zeropColumn->colValue->fArray[k] = (float)magValues[k];
          starIdColumn->colValue->sArray[k] =
                 magIdColumn->colValue->sArray[j];
          k++;
        }
      }
      break;
    }
  }


  /* 
   * Now do the fitting. NOTE: put to 0 rms of whatever has NOT been fitted 
   */

  if (k == 0) {
    cpl_msg_error(task, "No stars selected");
    return EXIT_FAILURE;
  }

  cpl_msg_info(task, " ");

  switch (mode) {
  case 0:
    if (k > 1) {
      xbiwt (fmagValues, k, &mean, &sigma, &dummy, &dummy1 );
      magzero = (double)mean ;
      rmszerop = (double)sigma ;
    } else {
      magzero = (double)fmagValues[0];
      rmszerop = 0. ;
    }
    rms =0;
    rmscolor = 0;
    cpl_msg_info(task, "Zeropoint: %5.2f +/- %5.2f", magzero, rmszerop);
    break;
  case 1:
    stupidLinearFit(extinctionValues, magValues, k, 
                    &magzero, &extinction, &rmszerop, &rms);
    rmscolor = 0;
    cpl_msg_info(task, "Zeropoint : %5.2f +/- %5.2f", magzero, rmszerop);
    cpl_msg_info(task, "Extinction: %5.2f +/- %5.2f", extinction, rms);
    break;
  case 2:   
    fitSurPolErrors(colorValues, k, &magzero,
                    &extinction, &colorTerm, &rmszerop, &rms, &rmscolor);
    cpl_msg_info(task, "Zeropoint : %5.2f +/- %5.2f", magzero, rmszerop);
    cpl_msg_info(task, "Extinction: %5.2f +/- %5.2f", extinction, rms);
    cpl_msg_info(task, "Color term: %5.2f +/- %5.2f", colorTerm, rmscolor);
    break;
  case 3:
    stupidLinearFit(colorTermValues, magValues, k, &magzero, 
                    &colorTerm, &rmszerop, &rmscolor);
    cpl_msg_info(task, "Zeropoint : %5.2f +/- %5.2f", magzero, rmszerop);
    cpl_msg_info(task, "Color term: %5.2f +/- %5.2f", colorTerm, rmscolor);
    rms = 0;
    break;
  }

  cpl_msg_info(task, "Number of stars used: %d", k);
  cpl_msg_info(task, " ");

  nameColumn->colValue->sArray = cpl_realloc(nameColumn->colValue->sArray, 
                                             k * sizeof(char *));
  starIdColumn->colValue->sArray = cpl_realloc(starIdColumn->colValue->sArray,
                                               k * sizeof(char *));
  zeropColumn->colValue->fArray = cpl_realloc(zeropColumn->colValue->fArray,
                                              k * sizeof(float));

  nameColumn->len = k;
  starIdColumn->len = k;
  zeropColumn->len = k;


  /* 
   * Write new photometric calibration coefficients in header keywords 
   */

    writeDoubleDescriptor(&(photometricTable->descs), 
                          pilTrnGetKeyword("MagZero"), magzero,
                          "Zero point magnitude");

    writeDoubleDescriptor(&(photometricTable->descs), 
                          "ESO PRO MAGZERO RMS", rmszerop,
                          "RMS on zero point magnitude");

    writeDoubleDescriptor(&(photometricTable->descs), 
                          pilTrnGetKeyword("Extinction"), extinction,
                          "Atmospheric extinction coefficient");

    writeDoubleDescriptor(&(photometricTable->descs), 
                          pilTrnGetKeyword("ExtinctionRms"), rms,
                          "RMS on extinction");

    writeStringDescriptor(&(photometricTable->descs), 
                          pilTrnGetKeyword("Colour"), colorSystem, 
                          "Color Index");

    writeDoubleDescriptor(&(photometricTable->descs), 
                          pilTrnGetKeyword("ColorTerm"), colorTerm,
                          "Color term for filter");

    writeDoubleDescriptor(&(photometricTable->descs), 
                          pilTrnGetKeyword("ColorTermRms"), rmscolor, 
                          "RMS on color term");

  /*
   * Create the product file on disk, set the product attributes and
   * update the set of frames.
   */

  vmstrlower(strcpy(outPhotTableName, photoCategory));
  /* strcat(outPhotTableName, ".TFITS"); */
  strcat(outPhotTableName, ".fits");

  if (writeFitsPhotometricTable(outPhotTableName, photometricTable)) {
    outputFrame = newPilFrame(outPhotTableName, photoCoeffCategory);

    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("MjdObs"), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("InstrumentMode"), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("READ.CLOCK"), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("WINi.BINX", 1), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("WINi.BINY", 1), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("SeqWindowSizeX", 1), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("SeqWindowSizeY", 1), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("Quadrant"), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("Adu2Electron", 1), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("Electron2Adu", 1), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("FilterId", quad), 1);
    pilFitsHdrCopy(outPhotTableName, 0, NULL, 
                   pilTrnGetKeyword("FilterName", quad), 1);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
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
                    "vmimcalphot",
    "Determine photometric calibration from a set of observed standard stars.",
    "This recipe is used to determine night zeropoints, atmospheric\n"
    "extinction coefficients, and colour terms, from a set of star match\n"
    "tables produced by the recipe vmimstandard. The star match tables may\n"
    "refer to different standard star fields, but they must all be derived\n"
    "from exposures made with the same filter and the same quadrant.\n\n"
    "Input files:\n\n"
    "  DO category:            Type:       Explanation:           Required:\n"
    "  IMG_STAR_MATCH_TABLE    Product     Lust of standard stars    Y\n"
    "  PHOTOMETRIC_TABLE       Calib       Photometric table         Y\n\n"
    "Output files:\n\n"
    "  DO category:            Data type:  Explanation:\n"
    "  PHOTOMETRIC_TABLE       FITS table  Upgraded photometric table\n\n"
    "The only product of this recipe is an upgraded PHOTOMETRIC_TABLE,\n"
    "carrying the newly computed zeropoint, and, if requested, new\n"
    "extinction and colour coefficients.\n\n"
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

                    vmimcalphot_create,
                    vmimcalphot_exec,
                    vmimcalphot_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
