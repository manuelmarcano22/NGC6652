/* $Id: vmifucalib.c,v 1.12 2013-08-07 16:30:28 cgarcia Exp $
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
 * $Date: 2013-08-07 16:30:28 $
 * $Revision: 1.12 $
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


/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

static int vmifucalib_create(cpl_plugin *) ;
static int vmifucalib_exec(cpl_plugin *) ;
static int vmifucalib_destroy(cpl_plugin *) ;

static cxint vmifucalib(PilSetOfFrames *);

/*-----------------------------------------------------------------------------
                            Structure definitions
 -----------------------------------------------------------------------------*/


/*
 * Definition of the label strings for all methods the recipe function
 * supports for combining frames and their associated method code.
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

static const char *biasMethodNames[] = {
  "Master",
  "Zmaster"
};

static const BiasMethod biasMethods[] = {
  BIAS_MASTER,
  BIAS_ZMASTER
};

static char vmifucalib_description[] =
"This recipe produces the extraction mask, the wavelength calibration,\n"
"and the relative transmission correction, from a set of flat field\n"
"and one arc lamp exposures (IFU mode).\n\n"
"Input files:\n\n"
"  DO category:                 Type:       Explanation:       Required:\n"
"  IFU_SCREEN_FLAT              Raw         Flat field            Y\n"
"  IFU_ARC_SPECTRUM             Raw         Arc spectrum          .\n"
"  LINE_CATALOG                 Calib       Arc line catalog      .\n"
"  MASTER_BIAS                  Calib       Master bias           Y\n"
"  IFU_IDENT                    Calib       Fiber identification  .\n"
"  CCD_TABLE                    Calib       Bad pixel table       .\n\n"
"Output files:\n\n"
"  DO category:                 Data type:  Explanation:\n"
"  IFU_MASTER_SCREEN_FLAT       FITS image  Stack of all raw flat fields\n"
"  IFU_ARC_SPECTRUM_EXTRACTED   FITS image  Extracted arc spectra\n"
"  IFU_FLAT_SPECTRUM_EXTRACTED  FITS image  Extracted flat field spectra\n"
"  IFU_TRACE                    FITS table  Extraction mask\n"
"  IFU_IDS                      FITS table  Wavelength calibration\n"
"  IFU_TRANSMISSION             FITS table  Transmission correction\n\n"
"An arc spectrum exposure may not be specified in input, but in that\n"
"case no wavelength calibration nor fiber relative transmission\n"
"correction will be computed, and further reduction of associated\n"
"scientific data will not be possible; the only products in this\n" 
"case will be the stack of all the input raw flat field exposures,\n" 
"and the extraction mask. If the arc spectrum exposure is specified\n" 
"in input, then also a line catalog must be given. The fiber\n" 
"identification table, if given, ensures that all spectra will be\n" 
"safely identified and assigned to the correct spatial position,\n" 
"but in most cases it will be possible to reduce the data even\n" 
"without specifying this table in the input SOF. A CCD table\n"
"must be specified only if a bad pixel cleaning is requested.\n\n"
"For more details, please refer to the VIMOS Pipeline User's Guide.";

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

    cpl_recipe *recipe = cx_calloc(1, sizeof *recipe);
    cpl_plugin *plugin = &recipe->interface;


    cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    VIMOS_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "vmifucalib",
                    "Generate IFU calibrations from a set of flat and arc lamp exposures.",
                    vmifucalib_description,
                    "ESO VIMOS Pipeline Team and VIMOS Consortium",
                    PACKAGE_BUGREPORT,
                    vimos_get_license(),
                    vmifucalib_create,
                    vmifucalib_exec,
                    vmifucalib_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}

/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmifucalib_create(cpl_plugin *plugin)
{
    cpl_recipe    *recipe;
    cpl_parameter *p;

/*    cx_string *path = cx_string_new(); */
    cxint status = 0;

    /* 
     * Check that the plugin is part of a valid recipe 
     */

    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE) 
        recipe = (cpl_recipe *)plugin;
    else 
        return -1;

    /* 
     * Create the parameters list in the cpl_recipe object 
     */

    recipe->parameters = cpl_parameterlist_new(); 
    if (recipe->parameters == NULL)
        return 1;

    /*
     * Fill the parameter list
     */

    p = cpl_parameter_new_value("vimos.Parameters.stacking.singleframes",
                                CPL_TYPE_BOOL,
                                "Flat fields combination method is ignored.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "AllowSingleFrames");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "AllowSingleFrames");
    cpl_parameterlist_append(recipe->parameters, p);

    
    p = cpl_parameter_new_enum("vimos.Parameters.stacking.method",
                               CPL_TYPE_STRING,
                               "Flat fields combination method.",
                               "vimos.Parameters",
                               "MinMax", 5, "Average", "Median", "MinMax",
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
                                0);
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
                                "arc lamp spectra",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ApplyTransmission");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ApplyTransmission");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.ids.maxrms",
                                CPL_TYPE_DOUBLE,
                                "Maximum tolerated RMS of residuals in IDS "
                                "fit",
                                "vimos.Parameters",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MaxIdsRms");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MaxIdsRms");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_enum("vimos.Parameters.line.ident.method",
                               CPL_TYPE_STRING,
                               "Arc line identification method.",
                               "vimos.Parameters",
                               "Blind", 2, "FirstGuess", "Blind");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "LineIdent");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "LineIdent");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.trace.failed",
                                CPL_TYPE_INT,
                                "Maximum percentage of rejected positions "
                                "in fiber spectra tracing",
                                "vimos.Parameters",
                                33);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "MaxTraceRejection");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "MaxTraceRejection");
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
vmifucalib_exec(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;

    cxint status = 0;

    PilSetOfFrames *sof = NULL;


    if (recipe->parameters == NULL || recipe->frames == NULL) {
        return 1;
    }

    /* Issue a banner */
    vimos_print_banner();


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

    if (vmifucalib(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmifucalib");
        
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
vmifucalib_destroy(cpl_plugin *plugin)
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
 *   Create all IFU calibrations needed by the VIMOS IFU DRS.
 * 
 * @return EXIT_SUCCESS or EXIT_FAILURE.   
 *
 * @param sof  Set of frames containing the references to raw flat field
 *             exposures, a raw arc lamp exposure, a master bias, and 
 *             (optionally) a CCD Table. 
 *  
 * @doc
 *   This recipe would combine the input flat field exposures, subtract
 *   the master bias, and then trace the flat field spectra on the
 *   reduced flat field image. On the basis of such tracings, the arc
 *   lamp spectra would be extracted and used to obtain a wavelength 
 *   calibration. Next, also the flat field spectra would be extracted, 
 *   and they would be integrated at a predetermined wavelength range,
 *   providing the transmission correction factors. The products of this
 *   recipe are the tracing solution, the wavelength calibration, and the
 *   transmission correction. If no arc lamp exposure is found in the
 *   input SOF, only the tracing solution would be produced. This recipe
 *   assumes that the flat field and the arc lamp exposures are taken
 *   quasi-simultaneously, so that the spectral traces obtained from
 *   the flat field can be used to extract the arc lamp spectra without
 *   any further alignment operation. A CCD table may also be specified, 
 *   to be used in bad pixel correction.
 *
 * @author C. Izzo, R. Palsa
 */

static cxint 
vmifucalib(PilSetOfFrames *sof)
{

  const char     task[] = "vmifucalib";

  const char     parameter[] = "Parameters";

  const char    *mFlatTag    = pilTrnGetCategory("IfuMasterScreenFlat");
  const char    *flatTag     = pilTrnGetCategory("IfuScreenFlat");
  const char    *arcTag      = pilTrnGetCategory("IfuArcSpectrum");
  const char    *arcExtrTag  = pilTrnGetCategory("IfuArcSpectrumExtracted");
  const char    *flatExtrTag = pilTrnGetCategory("IfuFlatSpectrumExtracted");
  const char    *idsTag      = pilTrnGetCategory("IfuIds");
  const char    *traceTag    = pilTrnGetCategory("IfuTrace");
  const char    *transTag    = pilTrnGetCategory("IfuTransmission");
  const char    *lineCatTag  = pilTrnGetCategory("LineCatalog");
  const char    *identTag    = pilTrnGetCategory("IfuIdent");

  char          *combMethodTag = NULL;
  char          *biasMethodTag = NULL;
  char           masterFlatName[PATHNAME_MAX + 1];
  char           extraName[PATHNAME_MAX + 1];
  char           traceName[PATHNAME_MAX + 1];
  char           transName[PATHNAME_MAX + 1];
  char           idsName[PATHNAME_MAX + 1];
  char           parName[50];
  char           colName[5];

  VimosBool      updateOK = VM_TRUE;

  int            flatCount, arcCount, tracedCount, wavecalCount;
  int            minFrames;

  int            lostFibers[] = {0, 0, 0, 0};
  int            wlostFibers[] = {0, 0, 0, 0};
  int            refRow[] = {0, 0, 0, 0};
  int            centralFiber[] = {0, 0, 0, 0};
  double         traceRms[] = {0.0, 0.0, 0.0, 0.0};
  double         grismAlignment[] = {0.0, 0.0, 0.0, 0.0};
  double        *mcoeff[] = {NULL, NULL, NULL, NULL};
  double        *tcoeff[] = {NULL, NULL, NULL, NULL};
  double         flatFlux = 0.0;
  double         errFlux = 0.0;

  char          *lineIdentificationTag = NULL;
  unsigned int   lineIdentification;

  unsigned int   cleanBadPixel;
  unsigned int   computeQC;
  unsigned int   applyTrans;
  unsigned int   singleFrames;
  unsigned int   error;

  int            biasMethodEntry;
  int            combMethodEntry;

  float         *data;
  float         *refdata;

  PilFrame      *currFrame, *ccdFrame, *mBiasFrame, *lineCatFrame, *identFrame;
  PilFrame      *outputFrame;

  BiasMethod     biasMethod = BIAS_UNDEF;
  CombMethod     combMethod = COMB_UNDEF;
  CombParameters combParameter;

  VimosImage    *mFlat = NULL;
  VimosImage   **flatList, *mBiasImage, *mBiasOver;
  VimosImage    *arcImage = NULL;
  VimosImage    *arcExtracted = NULL;
  VimosImage    *flatExtracted = NULL;

  VimosTable    *ccdTable = NULL;

  char        unit[20];
  char        comment[80];
  char        grismName[20];
  char        filterName[20];
  char        ifuShutter[10];

  int         quadrant;
  int         grism;
  int         aquadrant;
  int         agrism;
  int         xlen, ylen;
  int         reflength;
  int         order;           /* Order of polynomial fit to tracings      */
  int         row;             /* Reference image row                      */
  int         above;           /* Pixels to trace above reference row      */
  int         below;           /* Pixels to trace below reference row      */
  int         zero;            /* Rough position of zero order contam      */
  int         idsOrder;        /* Order of IDS polynomial                  */
  float       tolerance = 0.3; /* Max deviation from fit for rejection     */
  int         maxReject;       /* Max rejections to consider a fiber dead  */
  int         shortRange;      /* Radius of short trace                    */
  int         firstSlit, lastSlit, slit;
  int         medianFilterSize = 15;
  int         step = 15;   /* Pixels to skip in median filtering (speedup) */
  int         startPix, endPix;
  float       lambdaRed, lambdaYel, lambdaBlu;
  float       lambdaHe, lambdaNe, lambdaAr;
  double      meanValue, rmsValue;
  double      meanFwhm = 0;
  double      startLambda, endLambda, stepLambda, lambda;
  double      startTrans, endTrans;
  double      expTime, time;
  double      arcTime[3];
  double      cdelt;
  double     *coeff = NULL;
  double      maxIdsRms;
  double      lambda2pix;
  double      lflux, lflux_err;
  int         maxTraceRejection;
  int         radius = 18;
  int         wradius = 320;              /* 175; */
  int         count;
  int         i, j;

  cpl_table  *linecat = NULL;    /* Line catalog                            */
  cpl_table  *ident = NULL;      /* Fiber identification table              */
  cpl_table  *ftrace = NULL;     /* Long tracing on flat field              */
  cpl_table  *short_ftrace = NULL;
  cpl_table  *fmodel = NULL;     /* Fit of long tracing on flat field       */
  cpl_table  *fcoeff = NULL;     /* Coefficients of flat field tracings     */
  cpl_table  *short_fcoeff = NULL;

  cpl_table  *spectra = NULL;    /* Extracted arc lamp spectra              */
  cpl_table  *fspectra = NULL;   /* Extracted flat spectra                  */
  cpl_table  *ids = NULL;        /* IDS coefficients for each fiber         */
  cpl_table  *trans = NULL;      /* Transmission correction for each fiber  */
  cpl_table  *refident[4];       /* Identif. tables associated to refimage  */
  cpl_table **tables = NULL;     /* Auxiliary pointers to tables            */

  cpl_image  *flat = NULL;       /* Flat field image                        */
  cpl_image  *smo_flat = NULL;   /* Smoothed flat field image               */
  cpl_image  *arc = NULL;        /* Arc lamp spectrum                       */
  cpl_image  *extracted = NULL;  /* Image with extracted arc lamp spectra   */
  cpl_image  *fextracted = NULL; /* Image with extracted flat field spectra */
  cpl_image  *refimage = NULL;   /* Image of rows with identified spectra   */

  cpl_propertylist *palist = NULL;   /* Arc lamp frame header (for CPL)       */
  cpl_propertylist *pflist = NULL;   /* Flat lamp frame header (for CPL)      */
  cpl_propertylist *paheader = NULL; /* Primary header for output IDS table   */
  cpl_propertylist *pfheader = NULL; /* Primary header for output trans table */
  cpl_propertylist *ptheader = NULL; /* Primary header for output trace table */


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
   * Determine the frame stacking method and all method dependent parameters.
   */

  combMethodTag = (char *)pilDfsDbGetString(parameter, "StackMethod");

  if ((combMethodEntry = strselect(combMethodTag, combMethodNames, 
                                   PIL_N_ELEMENTS(combMethods))) < 0) {
    cpl_msg_error(task, "%s: Invalid frame combination method.", combMethodTag);
    return EXIT_FAILURE;
  }

  combMethod = combMethods[combMethodEntry];

  switch (combMethod) {
  case COMB_KSIGMA:
    minFrames = MIN_FRAMES_KSIGMA;
    combParameter.kSigmaLow = pilDfsDbGetDouble(parameter, "KSigmaLow", 5.);
    combParameter.kSigmaHigh = pilDfsDbGetDouble(parameter, "KSigmaHigh", 5.);
    break;

  case COMB_REJECT:
    minFrames = MIN_FRAMES_REJECT;
    combParameter.minRejection = pilDfsDbGetInt(parameter, "MinRejection", 1);
    combParameter.maxRejection = pilDfsDbGetInt(parameter, "MaxRejection", 1);
    break;

  case COMB_MEDIAN:
    minFrames = MIN_FRAMES_MEDIAN;
    break;

  case COMB_AVERAGE:
    minFrames = MIN_FRAMES_AVERAGE;
    break;

  default:
    cpl_msg_warning(task, "Invalid stacking method. Using default: 'Average'");
    combMethod = COMB_AVERAGE;
    minFrames = MIN_FRAMES_AVERAGE;
    break;
  }


  /*
   * Check if a single frame in input should be tolerated. This
   * means that if a single frame is found in input, then the
   * stacking method will be ignored.
   */

  singleFrames = pilDfsDbGetBool("Parameters", "AllowSingleFrames", 1);


  /*
   * Check if the bad pixels should be corrected.
   */

  cleanBadPixel = pilDfsDbGetBool(parameter, "CleanBadPixel", 1);


  /*
   * Check if the extracted arc lamp spectra should be also corrected
   * for the fiber transmission variations.
   */

  applyTrans = pilDfsDbGetBool("Parameters", "ApplyTransmission", 1);


  /*
   * Get maximum rms of residuals in IDS fit.
   */

  maxIdsRms = pilDfsDbGetDouble("Parameters", "MaxIdsRms", 1.0);


  /*
   * Get line identification method
   */

  lineIdentificationTag = (char *)pilDfsDbGetString(parameter, "LineIdent");

  if (strcmp(lineIdentificationTag, "FirstGuess") == 0)
    lineIdentification = 0;
  else
    lineIdentification = 1;


  /*
   * Get max percentage of rejected pixel positions during tracing
   * to consider a fiber "dead".
   */

  maxTraceRejection = pilDfsDbGetInt(parameter, "MaxTraceRejection", 33);

  if (maxTraceRejection < 0 || maxTraceRejection > 99) {
    cpl_msg_error(task, "Invalid percentage for MaxTraceRejection!");
    return EXIT_FAILURE;
  }


  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   * Make sure that the input dataset include the essential datasets.
   */

  arcCount = (int)pilSofFrameCount(sof, arcTag);

  if (arcCount > 1) {
    cpl_msg_error(task, "Too many arc lamp exposures in input: "
                "just one is required.");
    return EXIT_FAILURE;
  }

  flatCount = (int)pilSofFrameCount(sof, flatTag);

  if (flatCount == 0) {
    cpl_msg_error(task, "At least one flat field is required in input.");
    return EXIT_FAILURE;
  }

  if (flatCount < minFrames) {
    if (flatCount == 1 && singleFrames) {
      combMethod = COMB_UNDEF;
    }
    else {
      cpl_msg_error(task, "Not enough flat field frames in input for stacking "
                  "method '%s'!", combMethodTag);
      return EXIT_FAILURE;
    }
  }

  cpl_msg_info(task, "Loading the input frames...");


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
   * Load the arc lamp frame
   */

  if (arcCount) {
    currFrame = pilSofLookup(sof, arcTag);

    if ((arcImage = openOldFitsFile(pilFrmGetName(currFrame), 1, 0))) {
      pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
      closeFitsImage(arcImage, 0);

      /*
       * Load the header of the arc lamp exposure, needed
       * for the CPL based types used in this program.
       */

      palist = cpl_propertylist_load(pilFrmGetName(currFrame), 0);

      /*
       * Get also the exposure times for the individual lamps (for QC1 
       * flux measurement)
       */

      if (getArcLampTimes(arcImage, arcTime)) {
          
        /*
         * Header info is incompatible (e.g., ON lamp with exposure time
         * zero). Assume then that the correct exposure time is the same
         * for all lamps.
         */

        readDoubleDescriptor(arcImage->descs, 
                             pilTrnGetKeyword("ExposureTime"),
                             &time, NULL);

        for (i = 0; i < 3; i++) {
            arcTime[i] = time;
        }
      }
    }
    else {
      cpl_msg_error(task, "Failure loading the arc lamp frame");
      deleteImage(mBiasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }

    agrism = getGrism(arcImage);
    readIntDescriptor(arcImage->descs, pilTrnGetKeyword("Quadrant"),
                      &aquadrant, NULL);
  }
  else
    cpl_msg_warning(task, "No arc lamp exposure in input: no wavelength "
                  "calibration nor transmission correction can be computed.");

  if (arcImage == NULL) {
    if (computeQC) {
      cpl_msg_warning(task, "The QC1 parameters computation is disabled too.");
      computeQC = 0;
    }
  }


  /*
   * Load the flat field frames.
   */

  error = 1;

  if ((flatList = (VimosImage **)cpl_calloc(flatCount, sizeof(VimosImage *)))) {

    error = 0;
    currFrame = pilSofLookupNext(sof, flatTag);

    for (i = 0; i < flatCount; i++) {
      if ((flatList[i] = openOldFitsFile(pilFrmGetName(currFrame), 1, 0))) {
        pilFrmSetType(currFrame, PIL_FRAME_TYPE_RAW);
        closeFitsImage(flatList[i], 0);

        if (i == 0) {

          /*
           * Load the header of the first flat field exposure, needed 
           * for the CPL based types used in this program.
           */

          pflist = cpl_propertylist_load(pilFrmGetName(currFrame), 0);

          /*
           * Get also the exposure time (for QC1 flux measurement)
           */

          readDoubleDescriptor(flatList[i]->descs,
                               pilTrnGetKeyword("ExposureTime"), 
                               &expTime, NULL);

          if (expTime < 0.01) {
            error = 1;
            cpl_msg_error(task, "Invalid exposure time for flat field");
            deleteImage(flatList[0]);
            cpl_free(flatList);
            break;
          }

        }
      }
      else {
        error = 1;
        cpl_msg_error(task, "Failure loading flat field frame %d", i + 1);
        for (j = 0; j < i; j++) 
          deleteImage(flatList[j]);
        cpl_free(flatList);
        break;
      }
      currFrame = pilSofLookupNext(sof, NULL);
    }
  }
  else
    cpl_msg_error(task, "Failure in creating list of input flat fields");

  if (error) {
    cpl_propertylist_delete(pflist);
    cpl_propertylist_delete(palist);
    deleteImage(arcImage);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable); 
    return EXIT_FAILURE;
  }

  grism = getGrism(flatList[0]);
  readIntDescriptor(flatList[0]->descs, pilTrnGetKeyword("Quadrant"),
                    &quadrant, NULL);
  readStringDescriptor(flatList[0]->descs,
                    pilTrnGetKeyword("GrismName", quadrant),
                    grismName, NULL);
  readStringDescriptor(flatList[0]->descs,
                    pilTrnGetKeyword("FilterName", quadrant),
                    filterName, NULL);
  readStringDescriptor(flatList[0]->descs, pilTrnGetKeyword("IfuMode"),
                    ifuShutter, NULL);

  if (arcImage) {
    error = 0;
    if (grism != agrism) {
      error = 1;
      cpl_msg_error(task, "Input data belong to different grisms.");
    }
    if (quadrant != aquadrant) {
      error = 1;
      cpl_msg_error(task, "Input data belong to different quadrants.");
    }
    if (error) {
      for (i = 0; i < flatCount; i++) 
        deleteImage(flatList[i]);
      cpl_free(flatList);
      cpl_propertylist_delete(pflist);
      cpl_propertylist_delete(palist);
      deleteImage(arcImage);
      deleteImage(mBiasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

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
    for (i = 0; i < flatCount; i++) 
      deleteImage(flatList[i]);
    cpl_free(flatList);
    cpl_propertylist_delete(pflist);
    cpl_propertylist_delete(palist);
    deleteImage(arcImage);
    deleteImage(mBiasImage);
    deleteCcdTable(ccdTable);
    return EXIT_FAILURE;
  }

  /*
   * Recreate bias overscans using as reference the first flat field exposure
   */

  if ((mBiasOver = growOverscans(mBiasImage, flatList[0]))) {
    if (mBiasImage != mBiasOver) {
      deleteImage(mBiasImage);
      mBiasImage = mBiasOver;
    }
  }
  else {
    cpl_msg_error(task, "Failure in growing overscans in master bias");
    for (i = 0; i < flatCount; i++) 
      deleteImage(flatList[i]);
    cpl_free(flatList);
    cpl_propertylist_delete(pflist);
    cpl_propertylist_delete(palist);
    deleteImage(arcImage);
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

  if (arcImage) {
    if (VmSubBias(arcImage, mBiasImage, biasMethod) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot remove bias from arc lamp exposure");
      for (i = 0; i < flatCount; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);
      cpl_propertylist_delete(pflist);
      cpl_propertylist_delete(palist);
      deleteImage(arcImage);
      deleteImage(mBiasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  for (i = 0; i < flatCount; i++) {
    if (VmSubBias(flatList[i], mBiasImage, biasMethod) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot remove bias from input flat field %d", i + 1);
      for (i = 0; i < flatCount; i++)
        deleteImage(flatList[i]);
      cpl_free(flatList);
      cpl_propertylist_delete(pflist);
      cpl_propertylist_delete(palist);
      deleteImage(arcImage);
      deleteImage(mBiasImage);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  deleteImage(mBiasImage);


  /*
   * Stack the flat field frames.
   */

  if (flatCount == 1) {
    mFlat = duplicateImage(flatList[0]);
  }
  else {
    cpl_msg_info(task, "Combining %d flat fields with method '%s'", flatCount,
               combMethodNames[combMethodEntry]);
    mFlat = frComb(flatList, flatCount, combMethod, &combParameter, 0);
  }


  /*
   * Create the master flat header
   */

  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               pilTrnGetKeyword("ExposureTime"), NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               pilTrnGetKeyword("MjdObs"), NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               pilTrnGetKeyword("DateObs"), NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               pilTrnGetKeyword("ArchiveFile"), NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^C[A-Z]*[1,2]", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO OBS (DID|PROG ID)", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO OBS ID", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO TPL [.DINPSV.]", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO INS (DID|MODE)", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO INS (GRIS|FILT).*", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO INS IFU.*", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO DET [A-Z].*", NULL);
  vimosDscCopy(&mFlat->descs, flatList[0]->descs,
               "^ESO OCS (DID|CON QUAD)", NULL);

  for (i = 0; i < flatCount; i++)
    deleteImage(flatList[i]);
  cpl_free(flatList);


  if (cleanBadPixel) {
    cpl_msg_info(task, "Cleaning bad pixels from combined flat field");
    if (EXIT_FAILURE == cleanBadPixels(mFlat, ccdTable, 0)) {
      cpl_msg_error(task, "Cannot clean bad pixels");
      cpl_propertylist_delete(pflist);
      cpl_propertylist_delete(palist);
      deleteImage(arcImage);
      deleteImage(mFlat);
      deleteCcdTable(ccdTable);
      return EXIT_FAILURE;
    }
  }

  deleteCcdTable(ccdTable);


  /**************
   *  Here begins the real IFU data processing.
   */

  lineCatFrame = pilSofLookup(sof, lineCatTag);
  if (lineCatFrame)
      pilFrmSetType(lineCatFrame, PIL_FRAME_TYPE_CALIB);

  if (arcImage) {

    /*
     *  Line Catalog
     */

    if (lineCatFrame) {
      linecat = cpl_table_load((char *)pilFrmGetName(lineCatFrame), 1, 0);
    }
    else {
      cpl_msg_error(task, "No input line catalog found");
      cpl_propertylist_delete(pflist);
      cpl_propertylist_delete(palist);
      deleteImage(arcImage);
      deleteImage(mFlat);
      return EXIT_FAILURE;
    }


    /*
     *  Allocate the images that will contain the extracted spectra.
     */

    ifuRange(grism, &startLambda, &endLambda, &stepLambda); 
    ifuRangeTransmission(grism, &startTrans, &endTrans);

    xlen = (endLambda - startLambda) / stepLambda;

    if (grism > 1) {
      idsOrder = order = 4;
      firstSlit = 1;
      lastSlit = 1;
      ylen = 400;
    }
    else {
      idsOrder = order = 3;
      firstSlit = 0;
      lastSlit = 3;
      ylen = 1600;
    }

    if (computeQC)
      for (slit = firstSlit; slit <= lastSlit; slit++)
        mcoeff[slit] = cpl_calloc(idsOrder + 1, sizeof(double));

    if (computeQC) {
      for (slit = firstSlit; slit <= lastSlit; slit++) {
        tcoeff[slit] = cpl_calloc(order + 1, sizeof(double));
      }
    }

    extracted = cpl_image_new(xlen, ylen, CPL_TYPE_FLOAT);
    fextracted = cpl_image_new(xlen, ylen, CPL_TYPE_FLOAT);

  }
  else {
    if (grism > 1) {
      order = 4;
      firstSlit = 1;
      lastSlit = 1;
    }
    else {
      order = 3;
      firstSlit = 0;
      lastSlit = 3;
    }
  }


  /*
   *  arcImage and mFlat are cast to a CPL type, because the IFU DRS 
   *  is based on CPL.
   */

  flat = cpl_image_wrap_float(mFlat->xlen, mFlat->ylen, mFlat->data);
  if (arcImage)
    arc = cpl_image_wrap_float(arcImage->xlen, arcImage->ylen, arcImage->data);

  /*
   *  Build also the FITS (primary) headers for the product tables.
   *  First, the header for the trace table:
   */

  ptheader = cpl_propertylist_new();

  cpl_propertylist_append_string(ptheader, 
                                 pilTrnGetKeyword("DoCategory"), "none");
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("DateObs"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("MjdObs"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("InstrumentMode"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("PROG.ID"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("OBS.ID"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("TplId"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("IfuMode"));
  vm_plist_update(ptheader, pflist, "ESO INS IFUE MAG");
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("Quadrant"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("FilterName", quadrant));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("FilterId", quadrant));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("GrismName", quadrant));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("GrismId", quadrant));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("Adu2Electron", 1));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("Electron2Adu", 1));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("ReadNoise", 1));
  for (i = 1; i <= 4; i++)
    vm_plist_update(ptheader, pflist, pilTrnGetKeyword("BeamTemperature", i));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("READ.CLOCK"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("READ.MODE"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("READ.SPEED"));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("WINi.BINX", 1));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("WINi.BINY", 1));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("SeqWindowSizeX", 1));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("SeqWindowSizeY", 1));
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("MaskId", quadrant));
  vm_plist_update(ptheader, pflist, "^ESO INS LAMP[1-5] .*");
  vm_plist_update(ptheader, pflist, pilTrnGetKeyword("BeamTemperature"));

  cpl_propertylist_delete(pflist);

  if (arcImage) {

    pfheader = cpl_propertylist_duplicate(ptheader);

    if (computeQC) {

      /*
       *  Add placeholders for QC1 parameters (they will be filled later on)
       */

      for (slit = firstSlit; slit <= lastSlit; slit++) {
        sprintf(parName, "ESO QC IFU LOST%d", slit + 1);
        cpl_propertylist_append_int(ptheader, parName, 0);
      }

      for (slit = firstSlit; slit <= lastSlit; slit++) {
        sprintf(parName, "ESO QC IFU TRACE%d RMS", slit + 1);
        cpl_propertylist_append_double(ptheader, parName, 0.0);
      }

      for (slit = firstSlit; slit <= lastSlit; slit++) {
        sprintf(parName, "ESO QC IFU REFROW%d", slit + 1);
        cpl_propertylist_append_int(ptheader, parName, 0);
      }

      for (slit = firstSlit; slit <= lastSlit; slit++) {
        sprintf(parName, "ESO QC IFU TRACE%d CENTRAL", slit + 1);
        cpl_propertylist_append_int(ptheader, parName, 0);
      }

      for (slit = firstSlit; slit <= lastSlit; slit++) {
        sprintf(parName, "ESO QC IFU TRACE%d SLOPE", slit + 1);
        cpl_propertylist_append_double(ptheader, parName, 0.0);
      }

      for (slit = firstSlit; slit <= lastSlit; slit++) {
        for (i = 0; i <= order; i++) {
          sprintf(parName, "ESO QC IFU TRACE%d COEFF%d", slit + 1, i);
          cpl_propertylist_append_double(ptheader, parName, 0);
        }
      }

    }


    /*
     *  Next, the header for the transmission correction table:
     */

    if (computeQC) {

      /*
       *  Add placeholders for QC1 parameters
       */

      cpl_propertylist_append_double(pfheader, "ESO QC IFU FLUX LAMBDA1", 0.0);
      cpl_propertylist_append_double(pfheader, "ESO QC IFU FLUX LAMBDA2", 0.0);
      cpl_propertylist_append_double(pfheader, "ESO QC IFU FLAT FLUX", 0.0);
      cpl_propertylist_append_double(pfheader, "ESO QC IFU FLAT FLUXERR", 0.0);
    }


    /*
     *  Finally, the header for the IDS table
     */

    paheader = cpl_propertylist_new();

    cpl_propertylist_append_string(paheader, 
                                   pilTrnGetKeyword("DoCategory"), "none");
    vm_plist_update(paheader, palist, pilTrnGetKeyword("DateObs"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("MjdObs"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("InstrumentMode"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("PROG.ID"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("OBS.ID"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("TplId"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("IfuMode"));
    vm_plist_update(paheader, palist, "ESO INS IFUE MAG");
    vm_plist_update(paheader, palist, pilTrnGetKeyword("Quadrant"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("FilterName", quadrant));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("FilterId", quadrant));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("GrismName", quadrant));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("GrismId", quadrant));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("Adu2Electron", 1));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("Electron2Adu", 1));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("ReadNoise", 1));
    for (i = 1; i <= 4; i++)
      vm_plist_update(paheader, palist, pilTrnGetKeyword("BeamTemperature", i));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("READ.CLOCK"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("READ.MODE"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("READ.SPEED"));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("WINi.BINX", 1));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("WINi.BINY", 1));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("SeqWindowSizeX", 1));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("SeqWindowSizeY", 1));
    vm_plist_update(paheader, palist, pilTrnGetKeyword("MaskId", quadrant));
    vm_plist_update(paheader, palist, "^ESO INS LAMP[1-5] .*");
    vm_plist_update(paheader, palist, pilTrnGetKeyword("BeamTemperature"));

    cpl_propertylist_delete(palist);

    if (computeQC) {

      /*
       *  Add placeholders for QC1 parameters
       */

      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU IDS RMS", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION1", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION2", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION3", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION1 LAMBDA", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION2 LAMBDA", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION3 LAMBDA", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION1 RMS", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION2 RMS", 0.0);
      cpl_propertylist_append_double(paheader, 
                                     "ESO QC IFU RESOLUTION3 RMS", 0.0);

      for (slit = firstSlit; slit <= lastSlit; slit++) {
        sprintf(parName, "ESO QC IFU WAVELOST%d", slit + 1);
        cpl_propertylist_append_int(paheader, parName, 0);
        for (i = 0; i <= idsOrder; i++) {
          sprintf(parName, "ESO QC IFU WAVECAL%d COEFF%d", slit + 1, i);
          cpl_propertylist_append_double(paheader, parName, 0);
        }
      }
    }
  }


  /*
   *  Load reference for fiber identification, if present.
   */

  for (i = 0; i < 4; i++)
    refident[i] = NULL;

  if ((identFrame = pilSofLookup(sof, identTag))) {
    pilFrmSetType(identFrame, PIL_FRAME_TYPE_CALIB);
    refimage = cpl_image_load(pilFrmGetName(identFrame), CPL_TYPE_FLOAT, 0, 0);
    if (refimage) {
      refdata  = cpl_image_get_data(refimage);
      reflength = cpl_image_get_size_x(refimage);
      refdata += reflength * firstSlit;
      for (slit = firstSlit; slit <= lastSlit; slit++) {

        if (grism < 2 && ifuShutter[1] == 'N')
          if (slit != 1)
            continue;

        refident[slit] = cpl_table_load((char *)pilFrmGetName(identFrame), 
                                        slit - firstSlit + 1, 0);
        if (!refident[slit]) {
          for (i = 0; i < 4; i++)
            if (refident[i])
              cpl_table_delete(refident[i]);
          cpl_image_delete(refimage);
          cpl_propertylist_delete(paheader);
          cpl_propertylist_delete(pfheader);
          cpl_propertylist_delete(ptheader);
          cpl_image_unwrap(flat);
          cpl_image_unwrap(arc);
          deleteImage(arcImage);
          deleteImage(mFlat);
          cpl_image_delete(extracted);
          cpl_image_delete(fextracted);
          cpl_table_delete(linecat);
          for (slit = firstSlit; slit <= lastSlit; slit++)
            cpl_free(mcoeff[slit]);
          for (slit = firstSlit; slit <= lastSlit; slit++)
            cpl_free(tcoeff[slit]);
          cpl_msg_error(task, "Failure in loading the IFU fiber "
                      "identification tables!");
          return EXIT_FAILURE;
        }
      }
    }
    else {
      cpl_msg_error(task, "Failure in loading the IFU fiber "
                  "identification image!");
      cpl_propertylist_delete(paheader);
      cpl_propertylist_delete(pfheader);
      cpl_propertylist_delete(ptheader);
      cpl_image_unwrap(flat);
      cpl_image_unwrap(arc);
      deleteImage(arcImage);
      deleteImage(mFlat);
      cpl_image_delete(extracted);
      cpl_image_delete(fextracted);
      cpl_table_delete(linecat);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(mcoeff[slit]);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(tcoeff[slit]);
      return EXIT_FAILURE;
    }
  }
  else
    cpl_msg_warning(task, "No input fiber identification found - using "
                  "blind identification method.");

  for (slit = firstSlit; slit <= lastSlit; slit++) {

    if (grism < 2 && ifuShutter[1] == 'N')
      if (slit != 1)
        continue;

    cpl_msg_info(task, "Processing spectra in pseudo-slit %d:", slit + 1);
    cpl_msg_indent_more();
    cpl_msg_info(task, "Identify and trace flat field spectra...");

    /*
     *  Extraction parameters. Note that the FG (ident) is currently
     *  ignored (it is not yet fully implemented).
     */

    ifuExtractionParameters(grism, quadrant, slit, 1,
                            &row, &above, &below, &zero);

    refRow[slit] = row;
    cpl_msg_debug(task, "Using reference row %d for slit %d", row, slit + 1);
    cpl_msg_debug(task, "Extracting %d pixels above and %d below", above, below);


    /*
     *  This is the short extraction range, to be used as reference for 
     *  the tracing alignment to the science tracings (recipe vmifuscience).
     */

    shortRange = above < below ? above : below;
    if (shortRange > 400)
      shortRange = 400;


    /*
     *  Vertical median filter to flat.
     */

    smo_flat = cpl_image_vertical_median_filter(flat, medianFilterSize, row,
                                                above, below, step);

    /*
     *  Identify fibers on flat.
     */

    ident = NULL;

    if (identFrame) {
      if (ifuIdentifyUpgrade(smo_flat, row, refdata,
                             refident[slit], radius, wradius)) {
        cpl_msg_warning(task, "Upgrade of fiber identification failed - "
                      "using blind identification method");
        cpl_table_delete(refident[slit]);
      }
      else
        ident = refident[slit];
      refdata += reflength;
    }

    if (ident == NULL) {
      ident = ifuIdentify(smo_flat, row);
    }
// else

    if (ident == NULL) {
      cpl_msg_error(task, "Fiber identification failed!");
      cpl_propertylist_delete(paheader);
      cpl_propertylist_delete(pfheader);
      cpl_propertylist_delete(ptheader);
      cpl_image_delete(extracted);
      cpl_image_delete(fextracted);
      cpl_table_delete(linecat);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(mcoeff[slit]);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(tcoeff[slit]);
      cpl_image_delete(smo_flat);
      cpl_image_unwrap(flat);
      cpl_image_unwrap(arc);
      deleteImage(arcImage);
      deleteImage(mFlat);
      return EXIT_FAILURE;
    }


    /*
     *  Trace all fibers on smoothed flat field.
     */

    tables = ifuTrace(smo_flat, row, above, below, step, ident);

    if (tables == NULL) {
      cpl_msg_error(task, "Fiber tracing failed!");
      cpl_propertylist_delete(paheader);
      cpl_propertylist_delete(pfheader);
      cpl_propertylist_delete(ptheader);
      cpl_image_delete(extracted);
      cpl_image_delete(fextracted);
      cpl_table_delete(linecat);
      cpl_table_delete(ident);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(mcoeff[slit]);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(tcoeff[slit]);
      cpl_image_delete(smo_flat);
      cpl_image_unwrap(flat);
      cpl_image_unwrap(arc);
      deleteImage(arcImage);
      deleteImage(mFlat);
      return EXIT_FAILURE;
    }

    ftrace = tables[0];
    cpl_table_delete(tables[1]);
    cpl_free(tables);


    /*
     * Trace the flat field spectra also in a shorter range.
     */

    tables = ifuTrace(smo_flat, row, shortRange, shortRange, step, ident);

    short_ftrace = tables[0];
    cpl_table_delete(tables[1]);
    cpl_free(tables);
    cpl_table_delete(ident);
    cpl_image_delete(smo_flat);


    /*
     * This is the max number of rejections tolerated during tracing
     * a fiber. The formula below keeps into account the fact that
     * "step" pixels are skipped in the cpl_image_vertical_median_filter(),
     * and therefore less points are used for the tracing. The factor
     * maxTraceRejection is the percentage of positions rejections that
     * is tolerated, before a fiber is flagged as "dead".
     */

    maxReject = above - above / step + below - below / step
              + (maxTraceRejection / 100.) * (above + below) / step;

    tables = ifuFit(ftrace, order, tolerance, maxReject);
    cpl_table_delete(ftrace);

    fcoeff = tables[0];
    fmodel = tables[1];

    cpl_free(tables);

    lostFibers[slit] = cpl_table_count_invalid(fcoeff, "c0");
    tracedCount = 400 - lostFibers[slit];

    error = 1;

    if (tracedCount) {
      cpl_msg_indent_more();
      cpl_msg_info(task, "%d out of 400 spectra were identified and traced.", 
                 tracedCount);
      cpl_msg_indent_less();

      if (computeQC)
        traceRms[slit] = cpl_table_get_column_mean(fcoeff, "rms");

      if ((slit == firstSlit) || (grism < 2 && ifuShutter[1] == 'N')) {
        cpl_propertylist_set_string(ptheader, pilTrnGetKeyword("DoCategory"), 
                                    traceTag);
        vmstrlower(strcpy(traceName, traceTag));
        /* strcat(traceName, ".TFITS"); */
        strcat(traceName, ".fits");
        error = cpl_table_save(fcoeff, ptheader, NULL, 
                               traceName, CPL_IO_DEFAULT);
        cpl_propertylist_delete(ptheader);
        if (error)
          cpl_msg_error(task, "Cannot create local product file %s!", 
                        traceName);
      }
      else {
        error = cpl_table_save(fcoeff, NULL, NULL, traceName, CPL_IO_EXTEND);
        if (error)
          cpl_msg_error(task, "Cannot extend local product file %s!", 
                        traceName);
      }

      if (computeQC) {
        for (i = 0; i <= order; i++) {
          sprintf(colName, "c%d", i);
          tcoeff[slit][i] = cpl_table_get_column_median(fcoeff, colName);
        }
      }

      /*
       *  Fit also short tracings
       */

      maxReject = 2 * shortRange - 2 * (shortRange / step)
                + 0.5 * (2 * shortRange) / step;

      tables = ifuFit(short_ftrace, 1, tolerance, maxReject);

      short_fcoeff = tables[0];
      cpl_table_delete(tables[1]);
      cpl_table_delete(short_ftrace);
      cpl_free(tables);

      if (computeQC) {
        centralFiber[slit] = findCentralFiber(short_fcoeff, row);
        grismAlignment[slit] = cpl_table_get_double(short_fcoeff, "c1", 
                                                    centralFiber[slit], NULL);
      }

      error = cpl_table_save(short_fcoeff, NULL, NULL, 
                             traceName, CPL_IO_EXTEND);
      if (error)
        cpl_msg_error(task, "Cannot extend local product file %s!", traceName);

    }
    else
      cpl_msg_error(task, 
                  "No spectra could be traced on pseudo-slit %d!", slit + 1);

    cpl_table_delete(short_fcoeff);
    cpl_table_delete(fcoeff);

    if (error) {
      cpl_propertylist_delete(paheader);
      cpl_propertylist_delete(pfheader);
      cpl_image_delete(extracted);
      cpl_image_delete(fextracted);
      cpl_table_delete(linecat);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(mcoeff[slit]);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(tcoeff[slit]);
      cpl_table_delete(fmodel);
      cpl_image_unwrap(flat);
      cpl_image_unwrap(arc);
      deleteImage(arcImage);
      deleteImage(mFlat);
      return EXIT_FAILURE;
    }

    if (arcImage) {
      cpl_msg_info(task, "Extract flat field spectra...");
      fspectra = ifuExtraction(flat, fmodel);

      cpl_msg_info(task, "Extract arc lamp spectra...");
      spectra = ifuSimpleExtraction(arc, fmodel);

      cpl_table_delete(fmodel);

      if (lineIdentification == 0) {
        cpl_msg_info(task, "Compute the wavelength calibration...");
        coeff = ifuFirstIds(grism, quadrant, slit, &idsOrder, &lambda);
        ids = ifuComputeIds(spectra, linecat, coeff, idsOrder, lambda, 
                            zero, maxIdsRms);
        cpl_free(coeff);
      }
      else {
        cpl_msg_info(task, "Compute the wavelength calibration (blind)...");
        coeff = ifuFirstIds(grism, quadrant, slit, &idsOrder, &lambda);
        lambda2pix = coeff[1];
        cpl_free(coeff);
        coeff = ifuComputeIdsBlind(spectra, linecat, lambda2pix, idsOrder, 
                                   lambda, maxIdsRms);
        ids = ifuComputeIds(spectra, linecat, coeff, idsOrder, lambda, 
                            zero, maxIdsRms);
        cpl_free(coeff);
      }

      if (!ids) {
        cpl_msg_error(task, "Failure computing the wavelength calibration.");
        cpl_propertylist_delete(paheader);
        cpl_propertylist_delete(pfheader);
        cpl_image_delete(extracted);
        cpl_image_delete(fextracted);
        cpl_table_delete(linecat);
        for (slit = firstSlit; slit <= lastSlit; slit++)
          cpl_free(mcoeff[slit]);
        for (slit = firstSlit; slit <= lastSlit; slit++)
          cpl_free(tcoeff[slit]);
        cpl_table_delete(fspectra);
        cpl_table_delete(spectra);
        cpl_image_unwrap(flat);
        cpl_image_unwrap(arc);
        deleteImage(arcImage);
        deleteImage(mFlat);
        return EXIT_FAILURE;
      }

      wlostFibers[slit] = cpl_table_count_invalid(ids, "c0");
      wavecalCount = 400 - wlostFibers[slit];

      error = 1;

      if (wavecalCount) {
        cpl_msg_indent_more();
        cpl_msg_info(task, "%d out of %d extracted spectra were wavelength "
                   "calibrated.", wavecalCount, tracedCount);
        cpl_msg_indent_less();
  
        if ((slit == firstSlit) || (grism < 2 && ifuShutter[1] == 'N')) {
          cpl_propertylist_set_string(paheader, pilTrnGetKeyword("DoCategory"),
                                      idsTag);
          vmstrlower(strcpy(idsName, idsTag));
          /* strcat(idsName, ".TFITS"); */
          strcat(idsName, ".fits");
          error = cpl_table_save(ids, paheader, NULL, idsName, CPL_IO_DEFAULT);
          cpl_propertylist_delete(paheader);
          if (error)
            cpl_msg_error(task, "Cannot create local product file %s!", 
                          idsName);
        }
        else {
          error = cpl_table_save(ids, NULL, NULL, idsName, CPL_IO_EXTEND);
          if (error) 
            cpl_msg_error(task, "Cannot extend local product file %s!", 
                          idsName);
        }

        if (computeQC) {
          for (i = 0; i <= idsOrder; i++) {
            sprintf(colName, "c%d", i);
            mcoeff[slit][i] = cpl_table_get_column_median(ids, colName);
          }
        }
  
      }
      else
        cpl_msg_error(task, "No spectra could be wavelength calibrated "
                    "on pseudo-slit %d!", slit + 1);

      if (error) {
        cpl_propertylist_delete(pfheader);
        cpl_image_delete(extracted);
        cpl_image_delete(fextracted);
        cpl_table_delete(linecat);
        for (slit = firstSlit; slit <= lastSlit; slit++)
          cpl_free(mcoeff[slit]);
        for (slit = firstSlit; slit <= lastSlit; slit++)
          cpl_free(tcoeff[slit]);
        cpl_table_delete(ids);
        cpl_table_delete(fspectra);
        cpl_table_delete(spectra);
        cpl_image_unwrap(flat);
        cpl_image_unwrap(arc);
        deleteImage(arcImage);
        deleteImage(mFlat);
        return EXIT_FAILURE;
      }


      cpl_msg_info(task, "Resample flat field spectra at constant wavelength "
                 "step (%.2f Angstrom)", stepLambda);
      ifuResampleSpectra(fextracted, fspectra, ids, slit, lambda, startLambda,
                         stepLambda);
      cpl_table_delete(fspectra);

      cpl_msg_info(task, "Resample arc lamp spectra at constant wavelength step "
                 "(%.2f Angstrom)", stepLambda);
      ifuResampleSpectra(extracted, spectra, ids, slit, lambda, startLambda,
                         stepLambda);
      cpl_table_delete(spectra);

      cpl_table_delete(ids);

    }
    else
      cpl_table_delete(fmodel);

    cpl_msg_indent_less();

  }

  cpl_image_unwrap(flat);
  cpl_image_unwrap(arc);


  if (arcImage) {

    /*
     * Determine the transmission correction
     */

    cpl_msg_info(task, "Derive the relative transmission correction "
               "from extracted flat field spectra...");
    startPix = (startTrans - startLambda) / stepLambda;
    endPix = (endTrans - startLambda) / stepLambda;
    trans = ifuTransmission(fextracted, startPix, endPix, &flatFlux, &errFlux);
    flatFlux /= expTime;
    errFlux /= expTime;
    if (applyTrans)
      ifuApplyTransmission(fextracted, trans);
    cpl_propertylist_set_string(pfheader, pilTrnGetKeyword("DoCategory"), 
                                transTag);
    vmstrlower(strcpy(transName, transTag));
    /* strcat(transName, ".TFITS"); */
    strcat(transName, ".fits");
    error = cpl_table_save(trans, pfheader, NULL, transName, CPL_IO_DEFAULT);
    if (error) {
      cpl_msg_error(task, "Cannot create local product file %s!", transName);
      cpl_propertylist_delete(pfheader);
      cpl_image_delete(fextracted);
      cpl_image_delete(extracted);
      cpl_table_delete(linecat);
      cpl_table_delete(trans);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(mcoeff[slit]);
      for (slit = firstSlit; slit <= lastSlit; slit++)
        cpl_free(tcoeff[slit]);
      deleteImage(arcImage);
      deleteImage(mFlat);
      return EXIT_FAILURE;
    }

    if (applyTrans)
      ifuApplyTransmission(extracted, trans);

    cpl_table_delete(trans);
    cpl_propertylist_delete(pfheader);


    /*
     *  To output the CPL images, we convert them into VIMOS format first.
     *  Then we add the descriptor header, and compute the QC1 parameters.
     */

    data = (float *)cpl_image_get_data(extracted);
    arcExtracted = newImage(xlen, ylen, data);

    writeIntDescriptor(&(arcExtracted->descs),
                          pilTrnGetKeyword("Naxis", 1), xlen, 
                          pilTrnGetComment("Naxis"));
    writeIntDescriptor(&(arcExtracted->descs),
                          pilTrnGetKeyword("Naxis", 2) , ylen, 
                          pilTrnGetComment("Naxis"));
    writeDoubleDescriptor(&(arcExtracted->descs), 
                          pilTrnGetKeyword("Crval", 1), startLambda, 
                          pilTrnGetComment("Crval"));
    writeDoubleDescriptor(&(arcExtracted->descs), 
                          pilTrnGetKeyword("Crval", 2), 1., 
                          pilTrnGetComment("Crval"));
    writeDoubleDescriptor(&(arcExtracted->descs), 
                          pilTrnGetKeyword("Crpix", 1), 1., 
                          pilTrnGetComment("Crpix"));
    writeDoubleDescriptor(&(arcExtracted->descs), 
                          pilTrnGetKeyword("Crpix", 2), 1., 
                          pilTrnGetComment("Crpix"));
    writeDoubleDescriptor(&(arcExtracted->descs), 
                          pilTrnGetKeyword("Cdelt", 1), stepLambda, 
                          pilTrnGetComment("Cdelt"));
    writeFloatDescriptor(&(arcExtracted->descs), 
                          pilTrnGetKeyword("Cdelt", 2), 1., 
                          pilTrnGetComment("Cdelt"));
    writeStringDescriptor(&(arcExtracted->descs), 
                          pilTrnGetKeyword("Ctype", 1), "LAMBDA", 
                          pilTrnGetComment("Ctype"));
    writeStringDescriptor(&(arcExtracted->descs), 
                          pilTrnGetKeyword("Ctype", 2), "FIBRE", 
                          pilTrnGetComment("Ctype"));
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("Quadrant"), NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("MjdObs"), NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("DateObs"), NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("ExposureTime"), NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          "^ESO OBS (DID|PROG ID)", NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          "^ESO OBS ID", NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          "^ESO TPL [.DINPSV.]", NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          "^ESO INS (DID|MODE)", NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          "^ESO INS (GRIS|FILT).*", NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          "^ESO INS LAMP.*", NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          "^ESO INS IFU.*", NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          "^ESO DET [A-Z].*", NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("Adu2Electron", 1), NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("Electron2Adu", 1), NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("READ.CLOCK"), NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("READ.MODE"), NULL);
    vimosDscCopy(&arcExtracted->descs, arcImage->descs,
                          pilTrnGetKeyword("READ.SPEED"), NULL);

    data = (float *)cpl_image_get_data(fextracted);
    flatExtracted = newImage(xlen, ylen, data);

    writeIntDescriptor(&(flatExtracted->descs),
                          pilTrnGetKeyword("Naxis", 1), xlen, 
                          pilTrnGetComment("Naxis"));
    writeIntDescriptor(&(flatExtracted->descs),
                          pilTrnGetKeyword("Naxis", 2) , ylen, 
                          pilTrnGetComment("Naxis"));
    writeDoubleDescriptor(&(flatExtracted->descs), 
                          pilTrnGetKeyword("Crval", 1), startLambda, 
                          pilTrnGetComment("Crval"));
    writeDoubleDescriptor(&(flatExtracted->descs), 
                          pilTrnGetKeyword("Crval", 2), 1., 
                          pilTrnGetComment("Crval"));
    writeDoubleDescriptor(&(flatExtracted->descs), 
                          pilTrnGetKeyword("Crpix", 1), 1., 
                          pilTrnGetComment("Crpix"));
    writeDoubleDescriptor(&(flatExtracted->descs), 
                          pilTrnGetKeyword("Crpix", 2), 1., 
                          pilTrnGetComment("Crpix"));
    writeDoubleDescriptor(&(flatExtracted->descs), 
                          pilTrnGetKeyword("Cdelt", 1), stepLambda, 
                          pilTrnGetComment("Cdelt"));
    writeFloatDescriptor(&(flatExtracted->descs), 
                          pilTrnGetKeyword("Cdelt", 2), 1., 
                          pilTrnGetComment("Cdelt"));
    writeStringDescriptor(&(flatExtracted->descs), 
                          pilTrnGetKeyword("Ctype", 1), "LAMBDA", 
                          pilTrnGetComment("Ctype"));
    writeStringDescriptor(&(flatExtracted->descs), 
                          pilTrnGetKeyword("Ctype", 2), "FIBRE", 
                          pilTrnGetComment("Ctype"));
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("Quadrant"), NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("MjdObs"), NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("DateObs"), NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("ExposureTime"), NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          "^ESO OBS (DID|PROG ID)", NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          "^ESO OBS ID", NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          "^ESO TPL [.DINPSV.]", NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          "^ESO INS (DID|MODE)", NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          "^ESO INS (GRIS|FILT).*", NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          "^ESO INS LAMP.*", NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          "^ESO INS IFU.*", NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          "^ESO DET [A-Z].*", NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("Adu2Electron", 1), NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("Electron2Adu", 1), NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("READ.CLOCK"), NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("READ.MODE"), NULL);
    vimosDscCopy(&flatExtracted->descs, mFlat->descs,
                          pilTrnGetKeyword("READ.SPEED"), NULL);

  }
  else
    cpl_propertylist_delete(pfheader);

  if (computeQC) {

    if (arcImage) {
      cpl_msg_info(task, "Computing QC1 parameters...");


      /*
       *  Here are handled the parameters associated to the extracted
       *  arc spectrum (lamps monitoring).
       */
 
      if (pilQcGroupStart() == EXIT_SUCCESS) {
 
        /*
         * Write PRO.CATG, ARCFILE, TPL ID, OCS CON QUAD, etc. to QC1 group.
         */
 
        pilQcWriteString("PRO.CATG", arcExtrTag, "Product category");
 
        qcCopyValue(arcImage->descs, pilTrnGetKeyword("ArchiveFile"),
                    NULL, "Archive File Name");
 
        qcCopyValue(arcImage->descs, pilTrnGetKeyword("TplId"),
                    NULL, "Template signature ID");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("Quadrant"),
                    NULL, "Quadrant");
 
        qcCopyValue(arcImage->descs, pilTrnGetKeyword("GrismName", quadrant),
                    NULL, "Grism name");
 
        qcCopyValue(arcImage->descs, pilTrnGetKeyword("FilterName", quadrant),
                    NULL, "Filter name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("IfuMode"),
                    NULL, "IFU mode");
 
        qcCopyValue(arcImage->descs, "ESO INS IFUE MAG",
                    NULL, "IFU magnification");

        switch (grism) {
        case 0:                         /* LR_red             */
          lambdaHe = 7065.19;
          lambdaNe = 0.0;
          lambdaAr = 7723.80;
          break;
        case 1:                         /* LR_blue            */
          lambdaHe = 5015.68;
          lambdaNe = 6598.96;
          lambdaAr = 0.0;
          break;
        case 2:                         /* MR                 */
          lambdaHe = 7065.19;
          lambdaNe = 7032.41;
          lambdaAr = 7723.80;
          break;
        case 6:                         /* HR_red holographic */
        case 3:                         /* HR_red             */
          lambdaHe = 7065.19;
          lambdaNe = 7032.41;
          lambdaAr = 7723.80;
          break;
        case 4:                         /* HR_orange          */
          lambdaHe = 7065.19;
          lambdaNe = 7032.41;
          lambdaAr = 7383.98;    /* Previously was 7723.80, see DFS04209 */
          break;
        case 7:                         /* HR_blue  holographic */
          lambdaHe = 5015.68;
          lambdaNe = 0.0;
          lambdaAr = 0.0;
          break;
        case 5:                         /* HR_blue            */
          lambdaHe = 5015.68;
          lambdaNe = 5944.83;
          lambdaAr = 0.0;
          break;
        default:
          cpl_msg_error(task, "You should never read this message...");
          return EXIT_FAILURE;
        }

        if (lambdaHe > 1.) {
          extractIfuFlux(extracted, lambdaHe, startLambda, stepLambda,
                         &lflux, &lflux_err);

          if (arcTime[0] > 0.1) {
            lflux     /= arcTime[0];
            lflux_err /= arcTime[0];
          }
          else {
            lflux      = 0.0;
            lflux_err  = 0.0;
          }

          cpl_msg_info(task, "Flux of He %.2f: %.2f +/- %.2f ADU/fiber/s",
                     lambdaHe, lflux, lflux_err);

          qcWriteValueDouble(arcExtracted->descs, lambdaHe,
                             "QC.IFU.HE.LAMBDA", "Angstrom",
                             "He arc lamp line for flux determination");

          qcWriteValueDouble(arcExtracted->descs, lflux,
                             "QC.IFU.HE.FLUX", "ADU/fiber/s",
                             "Flux at chosen He wavelength");

          qcWriteValueDouble(arcExtracted->descs, lflux_err,
                             "QC.IFU.HE.FLUXERR", "ADU/fiber/s",
                             "Error on flux at chosen He wavelength");
        }
        else
          cpl_msg_warning(task, "No He lines in %s spectral range: "
                          "corresponding QC1 parameters are not computed.", 
                          grismName);

        if (lambdaNe > 1.) {
          extractIfuFlux(extracted, lambdaNe, startLambda, stepLambda,
                         &lflux, &lflux_err);

          if (arcTime[1] > 0.1) {
            lflux     /= arcTime[1];
            lflux_err /= arcTime[1];
          }
          else {
            lflux      = 0.0;
            lflux_err  = 0.0;
          }

          cpl_msg_info(task, "Flux of Ne %.2f: %.2f +/- %.2f ADU/fiber/s",
                     lambdaNe, lflux, lflux_err);

          qcWriteValueDouble(arcExtracted->descs, lambdaNe,
                             "QC.IFU.NE.LAMBDA", "Angstrom",
                             "Ne arc lamp line for flux determination");

          qcWriteValueDouble(arcExtracted->descs, lflux,
                             "QC.IFU.NE.FLUX", "ADU/fiber/s",
                             "Flux at chosen Ne wavelength");

          qcWriteValueDouble(arcExtracted->descs, lflux_err,
                             "QC.IFU.NE.FLUXERR", "ADU/fiber/s",
                             "Error on flux at chosen Ne wavelength");
        }
        else
          cpl_msg_warning(task, "No Ne lines in %s spectral range: "
                          "corresponding QC1 parameters are not computed.",
                          grismName);

        if (lambdaAr > 1.) {
          extractIfuFlux(extracted, lambdaAr, startLambda, stepLambda,
                         &lflux, &lflux_err);

          if (arcTime[2] > 0.1) {
            lflux     /= arcTime[2];
            lflux_err /= arcTime[2];
          }
          else {
            lflux      = 0.0;
            lflux_err  = 0.0;
          }

          cpl_msg_info(task, "Flux of Ar %.2f: %.2f +/- %.2f ADU/fiber/s",
                     lambdaAr, lflux, lflux_err);

          qcWriteValueDouble(arcExtracted->descs, lambdaAr,
                             "QC.IFU.AR.LAMBDA", "Angstrom",
                             "Ar arc lamp line for flux determination");

          qcWriteValueDouble(arcExtracted->descs, lflux,
                             "QC.IFU.AR.FLUX", "ADU/fiber/s",
                             "Flux at chosen Ar wavelength");

          qcWriteValueDouble(arcExtracted->descs, lflux_err,
                             "QC.IFU.AR.FLUXERR", "ADU/fiber/s",
                             "Error on flux at chosen Ar wavelength");
        }
        else
          cpl_msg_warning(task, "No Ar lines in %s spectral range: "
                          "corresponding QC1 parameters are not computed.",
                          grismName);

        if (pilQcGroupEnd() == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

      }
      else
        cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");


      /*
       *  Here are handled the parameters associated to the trace table.
       */

      if (pilQcGroupStart() == EXIT_SUCCESS) {

        /*
         * Write PRO.CATG, ARCFILE, TPL ID, OCS CON QUAD, etc. to QC1 group.
         */

        pilQcWriteString("PRO.CATG", traceTag, "Product category");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("ArchiveFile"),
                    NULL, "Archive File Name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("TplId"),
                    NULL, "Template signature ID");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("Quadrant"),
                    NULL, "Quadrant");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("GrismName", quadrant),
                    NULL, "Grism name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("FilterName", quadrant),
                    NULL, "Filter name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("BeamTemperature", 
                    quadrant), NULL, "Beam temperature");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("RotatorAngleStart"), 
                    NULL, "Adaptor initial rotator angle");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("RotatorAngleEnd"), 
                    NULL, "Adaptor final rotator angle");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("IfuMode"), 
                    NULL, "IFU mode");

        qcCopyValue(arcImage->descs, "ESO INS IFUE MAG", 
                    NULL, "IFU magnification");

        pilQcWriteDouble("PRO.WLEN.CEN", lambda, "Angstrom", 
                         "Reference wavelength");


        /* QC IFU REFROWi */

        for (slit = firstSlit; slit <= lastSlit; slit++) {
          sprintf(parName, "QC.IFU.REFROW%d", slit + 1);
          qcWriteValueInt_CPL(traceName, refRow[slit], parName, "pixel",
                        "CCD Y pixel coordinate where the tracing is started");
        }


        /* QC IFU LOSTi */

        for (slit = firstSlit; slit <= lastSlit; slit++) {
          sprintf(parName, "QC.IFU.LOST%d", slit + 1);
          qcWriteValueInt_CPL(traceName, lostFibers[slit], parName, NULL,
                        "Lost fibers on pseudo-slit");
        }


        /* QC IFU TRACEi RMS */

        for (slit = firstSlit; slit <= lastSlit; slit++) {
          sprintf(parName, "QC.IFU.TRACE%d.RMS", slit + 1);
          qcWriteValueDouble_CPL(traceName, traceRms[slit], parName, "pixel", 
                        "Mean RMS of polynomial fits of spectral tracings");
        }


        /* QC IFU TRACEi CENTRAL */

        for (slit = firstSlit; slit <= lastSlit; slit++) {
          sprintf(parName, "QC.IFU.TRACE%d.CENTRAL", slit + 1);
          qcWriteValueInt_CPL(traceName, centralFiber[slit], parName, NULL,
                "Fiber closest to central CCD X pixel at the reference row");
        }


        /* QC IFU TRACEi SLOPE */

        for (slit = firstSlit; slit <= lastSlit; slit++) {
          sprintf(parName, "QC.IFU.TRACE%d.SLOPE", slit + 1);
          qcWriteValueDouble_CPL(traceName, grismAlignment[slit], parName, NULL,
                "Slope of tracing of the central fiber at reference row");
        }

        for (slit = firstSlit; slit <= lastSlit; slit++) {
          for (i = 0; i <= order; i++) {

            sprintf(parName, "QC.IFU.TRACE%d.COEFF%d", slit + 1, i);

            switch (i) {
            case 0:
              sprintf(unit, "pixel"); break;
            case 1:
              sprintf(unit, "pixel/pixel"); break;
            default:
              sprintf(unit, "pixel/pixel^%d", i);
            }

            sprintf(comment, "Median coefficient %d of tracing", i);

            qcWriteValueDouble_CPL(traceName, tcoeff[slit][i],
                                   parName, unit, comment);

          }
          cpl_free(tcoeff[slit]);
        }

        if (pilQcGroupEnd() == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

      }
      else
        cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");


      /*
       *  Here are handled the parameters associated to the trans table.
       */

      if (pilQcGroupStart() == EXIT_SUCCESS) {

        /*
         * Write PRO.CATG, ARCFILE, TPL ID, OCS CON QUAD, etc. to QC1 group.
         */

        pilQcWriteString("PRO.CATG", transTag, "Product category");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("ArchiveFile"),
                    NULL, "Archive File Name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("TplId"),
                    NULL, "Template signature ID");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("Quadrant"),
                    NULL, "Quadrant");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("GrismName", quadrant),
                    NULL, "Grism name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("FilterName", quadrant),
                    NULL, "Filter name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("BeamTemperature", 
                    quadrant), NULL, "Beam temperature");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("RotatorAngleStart"), 
                    NULL, "Adaptor initial rotator angle");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("RotatorAngleEnd"), 
                    NULL, "Adaptor final rotator angle");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("IfuMode"),
                    NULL, "IFU mode");

        qcCopyValue(arcImage->descs, "ESO INS IFUE MAG",         
                    NULL, "IFU magnification");

        pilQcWriteDouble("PRO.WLEN.CEN", lambda, "Angstrom", 
                         "Reference wavelength");


        /* QC IFU FLUX LAMBDAi */

        qcWriteValueDouble_CPL(transName, startTrans, "QC.IFU.FLUX.LAMBDA1", 
                               "Angstrom", "Start wavelength in flat field "
                               "flux determination");

        qcWriteValueDouble_CPL(transName, endTrans, "QC.IFU.FLUX.LAMBDA2", 
                               "Angstrom", "End wavelength in flat field "
                               "flux determination");


        /* QC IFU FLAT FLUX */

        qcWriteValueDouble_CPL(transName, flatFlux, "QC.IFU.FLAT.FLUX", 
                               "ADU/s", "Mean flat field signal within "
                               "given spectral range");

        /* QC IFU FLAT FLUXERR */

        qcWriteValueDouble_CPL(transName, errFlux, "QC.IFU.FLAT.FLUXERR", 
                               "ADU/s", "Error on QC.IFU.FLAT.FLUX");

        if (pilQcGroupEnd() == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

      }
      else
        cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");


      /*
       *  Here are handled the parameters associated to the IDS table.
       */

      if (pilQcGroupStart() == EXIT_SUCCESS) {

        /*
         * Write PRO.CATG, ARCFILE, TPL ID, OCS CON QUAD, etc. to QC1 group.
         */

        pilQcWriteString("PRO.CATG", idsTag, "Product category");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("ArchiveFile"),
                    NULL, "Archive File Name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("TplId"),
                    NULL, "Template signature ID");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("Quadrant"),
                    NULL, "Quadrant");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("GrismName", quadrant),
                    NULL, "Grism name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("FilterName", quadrant),
                    NULL, "Filter name");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("BeamTemperature", 
                    quadrant), NULL, "Beam temperature");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("RotatorAngleStart"), 
                    NULL, "Adaptor initial rotator angle");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("RotatorAngleEnd"), 
                    NULL, "Adaptor final rotator angle");

        qcCopyValue(arcImage->descs, pilTrnGetKeyword("IfuMode"),
                    NULL, "IFU mode");

        qcCopyValue(arcImage->descs, "ESO INS IFUE MAG",         
                    NULL, "IFU magnification");

        pilQcWriteDouble("PRO.WLEN.CEN", lambda, "Angstrom", 
                         "Reference wavelength");


        switch (grism) {
        case 0:                         /* LR_red             */
          lambdaRed = 9122.97;
          lambdaYel = 7635.11;
          lambdaBlu = 5875.62;
          break;
        case 1:                         /* LR_blue            */
          lambdaRed = 6598.95;
          lambdaYel = 5015.68;
          lambdaBlu = 4471.48;
          break;
        case 2:                         /* MR                 */
          lambdaRed = 8264.521;
          lambdaYel = 6678.200;
          lambdaBlu = 5015.675;
          break;
        case 6:                         /* HR_red holographic */
        case 3:                         /* HR_red             */
          lambdaRed = 8264.521;
          lambdaYel = 7438.899;
          lambdaBlu = 6506.528;
          break;
        case 4:                         /* HR_orange          */
          lambdaRed = 7245.167;
          lambdaYel = 6506.528;
          lambdaBlu = 5852.488;
          break;
        case 5:                         /* HR_blue            */
          lambdaRed = 6029.997;
          lambdaYel = 5015.675;
          lambdaBlu = 4471.480;
          break;
        case 7:                         /* HR_blue holographic */
          lambdaRed = 5015.675;
          lambdaYel = 4471.480;
          lambdaBlu = 3888.650;
          break;
        default:
          cpl_msg_error(task, "You should never read this message...");
          return EXIT_FAILURE;
        }


        count = 0;
        if (spectralResolution(arcExtracted, lambdaRed, 
                               &meanValue, &rmsValue, 100000000)
                               == EXIT_SUCCESS) {
          count++;
          meanFwhm += lambdaRed / meanValue;
        }

        cpl_msg_info(task, "Spectral resolution at %.2f: %.2f +/- %.2f",
                   lambdaRed, meanValue, rmsValue);

        qcWriteValueDouble_CPL(idsName, lambdaRed,
                           "QC.IFU.RESOLUTION1.LAMBDA", "Angstrom",
                           "Line used in spectral resolution determination");

        qcWriteValueDouble_CPL(idsName, meanValue,
                           "QC.IFU.RESOLUTION1", NULL,
                           "Mean spectral resolution at red end of spectra");

        qcWriteValueDouble_CPL(idsName, rmsValue,
                           "QC.IFU.RESOLUTION1.RMS", NULL,
                           "RMS of spectral resolution at red end of spectra");

        if (spectralResolution(arcExtracted, lambdaYel,
                               &meanValue, &rmsValue, 100000000)
                               == EXIT_SUCCESS) {
          count++;
          meanFwhm += lambdaYel / meanValue;
        }

        cpl_msg_info(task, "Spectral resolution at %.2f: %.2f +/- %.2f",
                   lambdaYel, meanValue, rmsValue);

        qcWriteValueDouble_CPL(idsName, lambdaYel,
                           "QC.IFU.RESOLUTION2.LAMBDA", "Angstrom",
                           "Line used in spectral resolution determination");

        qcWriteValueDouble_CPL(idsName, meanValue,
                           "QC.IFU.RESOLUTION2", NULL,
                           "Mean spectral resolution at center of spectra");

        qcWriteValueDouble_CPL(idsName, rmsValue,
                           "QC.IFU.RESOLUTION2.RMS", NULL,
                           "RMS of spectral resolution at center of spectra");

        if (spectralResolution(arcExtracted, lambdaBlu,
                               &meanValue, &rmsValue, 100000000)
                               == EXIT_SUCCESS) {
          count++;
          meanFwhm += lambdaBlu / meanValue;
        }

        cpl_msg_info(task, "Spectral resolution at %.2f: %.2f +/- %.2f",
                   lambdaBlu, meanValue, rmsValue);

        qcWriteValueDouble_CPL(idsName, lambdaBlu,
                           "QC.IFU.RESOLUTION3.LAMBDA", "Angstrom",
                           "Line used in spectral resolution determination");

        qcWriteValueDouble_CPL(idsName, meanValue,
                           "QC.IFU.RESOLUTION3", NULL,
                           "Mean spectral resolution at blue end of spectra");

        qcWriteValueDouble_CPL(idsName, rmsValue,
                           "QC.IFU.RESOLUTION3.RMS", NULL,
                           "RMS of spectral resolution at blue end of spectra");

        if (count)
          meanFwhm /= count;          /* In Angstrom */
        else {
          readDoubleDescriptor(arcExtracted->descs, 
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

        rmsValue = distortionsRms_CPL(arcExtracted, linecat, meanFwhm);

        cpl_table_delete(linecat);

        cpl_msg_info(task, "IDS RMS: %.3f pixel", rmsValue);
        qcWriteValueDouble_CPL(idsName, rmsValue, "QC.IFU.IDS.RMS", "pixel", 
                               "Mean IDS rms");

        for (slit = firstSlit; slit <= lastSlit; slit++) {

          sprintf(parName, "QC.IFU.WAVELOST%d", slit + 1);
          qcWriteValueInt_CPL(idsName, wlostFibers[slit], parName, NULL, 
                              "Fibers lost both to trace and wavecalib");

          for (i = 0; i <= idsOrder; i++) {

            sprintf(parName, "QC.IFU.WAVECAL%d.COEFF%d", slit + 1, i);

            switch (i) {
            case 0:
              sprintf(unit, "pixel"); break;
            case 1:
              sprintf(unit, "pixel/Angstrom"); break;
            default:
              sprintf(unit, "pixel/Angstrom^%d", i);
            }

            sprintf(comment, "Median coefficient %d of IDS", i);

            qcWriteValueDouble_CPL(idsName, mcoeff[slit][i],
                                   parName, unit, comment);

          }
          cpl_free(mcoeff[slit]);
        }

        if (pilQcGroupEnd() == EXIT_FAILURE)
          cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

      }
      else
        cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

      deleteImage(arcImage);

    }

  } /* End of QC1 computation. */


  /*
   * Update the master flat field header
   */

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataMin"),
                                     imageMinimum(mFlat), 
                                     pilTrnGetComment("DataMin"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataMax"),
                                     imageMaximum(mFlat), 
                                     pilTrnGetComment("DataMax"),
                                     "ESO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("SummedExposureTime"),
                                     flatCount * expTime, 
                                     pilTrnGetComment("SummedExposureTime"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertIntDescriptor(&(mFlat->descs), 
                                  pilTrnGetKeyword("NFramesCombined"),
                                  flatCount,
                                  pilTrnGetComment("NFramesCombined"),
                                  "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataMedian"),
                                     imageMedian(mFlat), 
                                     pilTrnGetComment("DataMedian"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataStdDeviation"),
                                     imageSigma(mFlat), 
                                     pilTrnGetComment("DataStdDeviation"),
                                     "ESO PRO*", 1));

  updateOK = (updateOK &&
              insertDoubleDescriptor(&(mFlat->descs), 
                                     pilTrnGetKeyword("DataMean"), 
                                     imageMean(mFlat), 
                                     pilTrnGetComment("DataMean"),
                                     "ESO PRO*", 1));

  if (!updateOK) {
    cpl_msg_error(task, "Failure updating product header");
    if (arcImage) {
      arcExtracted->data = NULL;
      deleteImage(arcExtracted);
      cpl_image_delete(extracted);
      flatExtracted->data = NULL;
      deleteImage(flatExtracted);
      cpl_image_delete(fextracted);
    }
    deleteImage(mFlat);
    return EXIT_FAILURE;
  }


  /*
   * Create the product file on disk, set the product attributes and
   * update the set of frames.
   */

  vmstrlower(strcpy(masterFlatName, mFlatTag));
  strcat(masterFlatName, ".fits");

  if (createFitsImage(masterFlatName, mFlat, mFlatTag)) {
    outputFrame = newPilFrame(masterFlatName, mFlatTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", masterFlatName);
    if (arcImage) {
      arcExtracted->data = NULL;
      deleteImage(arcExtracted);
      cpl_image_delete(extracted);
      flatExtracted->data = NULL;
      deleteImage(flatExtracted);
      cpl_image_delete(fextracted);
    }
    deleteImage(mFlat);
    return EXIT_FAILURE;
  }

  deleteImage(mFlat);

  outputFrame = newPilFrame(traceName, traceTag);

  pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
  pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
  pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
  pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

  pilSofInsert(sof, outputFrame);

  if (arcImage) {

    outputFrame = newPilFrame(idsName, idsTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);

    vmstrlower(strcpy(extraName, arcExtrTag));
    strcat(extraName, ".fits");

    if (createFitsImage(extraName, arcExtracted, arcExtrTag)) {
      outputFrame = newPilFrame(extraName, arcExtrTag);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
    }
    else {
      cpl_msg_error(task, "Cannot create local product file %s!", extraName);
      arcExtracted->data = NULL;
      deleteImage(arcExtracted);
      cpl_image_delete(extracted);
      flatExtracted->data = NULL;
      deleteImage(flatExtracted);
      cpl_image_delete(fextracted);
      return EXIT_FAILURE;
    }

    arcExtracted->data = NULL;
    deleteImage(arcExtracted);
    cpl_image_delete(extracted);

    vmstrlower(strcpy(extraName, flatExtrTag));
    strcat(extraName, ".fits");

    if (createFitsImage(extraName, flatExtracted, flatExtrTag)) {
      outputFrame = newPilFrame(extraName, flatExtrTag);

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
    }
    else {
      cpl_msg_error(task, "Cannot create local product file %s!", extraName);
      flatExtracted->data = NULL;
      deleteImage(flatExtracted);
      cpl_image_delete(fextracted);
      return EXIT_FAILURE;
    }

    flatExtracted->data = NULL;
    deleteImage(flatExtracted);
    cpl_image_delete(fextracted);

    outputFrame = newPilFrame(transName, transTag);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_TABLE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }

  return EXIT_SUCCESS;

}


