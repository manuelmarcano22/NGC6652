/* $Id: vmmasktoccd.c,v 1.12 2012-06-12 09:34:40 cgarcia Exp $
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
 * $Date: 2012-06-12 09:34:40 $
 * $Revision: 1.12 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <string.h>
#include <unistd.h>
#include <math.h>

#include <cxmemory.h>
#include <cxstring.h>

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
#include <cpl_array.h>
#include <cpl_table.h>
#include <cpl_image.h>
#include <cpl_fit.h>
#include <piltask.h>
#include <pilutils.h>
#include <pilqc.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmfit.h"
#include "vmutils.h"
#include "vmsextractor.h"
#include "vmimgpreprocessing.h"
#include "vmimgutils.h"
#include "vmqcutils.h"
#include "vmcpl.h"
#include "vmmath.h"


#define SEXTRACTOR_ARGC  20
#define MAX_LINE_LENGTH  520


/*
 * Definition of the label strings for all methods the recipe function
 * supports for combining frames and bias removal, with their associated
 * method code.
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


static SextParameter masktoccd_1[] = {
    {"NUMBER", SEXT_COLUMN_INT},
    {"FLUX_APER(1)", SEXT_COLUMN_FLOAT},
    {"FLUX_MAX", SEXT_COLUMN_FLOAT},
    {"X_IMAGE", SEXT_COLUMN_FLOAT},
    {"Y_IMAGE", SEXT_COLUMN_FLOAT},
    {"ELONGATION", SEXT_COLUMN_FLOAT},
    {"FLUX_RADIUS", SEXT_COLUMN_FLOAT},
    {"FLAGS", SEXT_COLUMN_INT},
    {"VIGNET(15,15)", SEXT_COLUMN_FLOAT},
    {0, SEXT_COLUMN_UNDEF}
};

static SextParameter masktoccd_2[] = {
    {"NUMBER", SEXT_COLUMN_INT},
    {"XPSF_IMAGE", SEXT_COLUMN_FLOAT},
    {"YPSF_IMAGE", SEXT_COLUMN_FLOAT},
    {"FLUX_PSF", SEXT_COLUMN_FLOAT},
    {"X_IMAGE", SEXT_COLUMN_FLOAT},
    {"Y_IMAGE", SEXT_COLUMN_FLOAT},
    {"ELONGATION", SEXT_COLUMN_FLOAT},
    {"MAG_AUTO", SEXT_COLUMN_FLOAT},
    {"FLAGS", SEXT_COLUMN_INT},
    {"ELLIPTICITY", SEXT_COLUMN_FLOAT},
    {"FWHM_WORLD", SEXT_COLUMN_FLOAT},
    {"THETA_IMAGE", SEXT_COLUMN_FLOAT},
    {"A_IMAGE", SEXT_COLUMN_FLOAT},
    {"B_IMAGE", SEXT_COLUMN_FLOAT},
    {0, SEXT_COLUMN_UNDEF}
};


static cxint vmmasktoccd(PilSetOfFrames *);


/*
 * Matches mask spots detected by SExtractor with positions written 
 * in the image header. It does the simplest possible matching, going
 * through the list of sex_pixels to find out those that are within
 * distance r1 from ref_pixel. If more than one match is found, it
 * returs the SExtractor detection with the highest flux associated.
 */

static VimosPixel *
matchSexPosition(VimosPixel *sex_pixel, int nSex, VimosPixel *ref_pixel,
                 double r1)
{

    const char *modName = "matchSexPosition";


    int         matchIndex[10];
    int         match;
    int         matchFound = 0;
    int         nMatch = 0;
    int         i;

    float       matchInt;
    float       minInt = -999.;

    double      inx, iny;
    double      dx,dy,dxy;

    VimosPixel *corrPixel;


    if (sex_pixel == NULL) {
        cpl_msg_error(modName, "No SExtractor detections: cannot match "
                    "positions");
        return NULL;
    }
  
    if (ref_pixel == NULL) {
        cpl_msg_error(modName, "No reference pixel: cannot match positions");
        return NULL;
    } 

    if (r1 < 1.0) {
        cpl_msg_error(modName, "Wrong radius values: %g", r1);
        return NULL;
    }
        
    /* Allocate output pixel */

    corrPixel = newPixel(1);
    inx = ref_pixel->x;
    iny = ref_pixel->y;

    for (i = 0; i < nSex; i++)
    {
        dx = sex_pixel[i].x - inx;
        dy = sex_pixel[i].y - iny;
        dxy = sqrt(dx * dx + dy * dy);
        if (dxy < r1)
        {
            matchIndex[nMatch] = i;
            matchFound=1;
            nMatch++;
        }
    }

    if (matchFound == 0) {
        return NULL;
    }

    for (i = 0; i < nMatch; i++) {
        matchInt = sex_pixel[matchIndex[i]].i;
        if (matchInt > minInt) {
            minInt = matchInt;
            match = matchIndex[i];
        }
    }

    corrPixel->x = sex_pixel[match].x;
    corrPixel->y = sex_pixel[match].y;
    corrPixel->i = sex_pixel[match].i;

    return (corrPixel);

}


/*
 * Create the recipe instance, i.e. setup the parameter list for this
 * recipe and make it available to the application using the interface.
 */

static cxint
vmmasktoccd_create(cpl_plugin *plugin)
{

    cpl_recipe *recipe = (cpl_recipe *)plugin;

    cpl_parameter *p;

    cx_string *path = cx_string_new();

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
                                "Cosmic ray removal from input science image.",
                                "vimos.Parameters",
                                FALSE);
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
                                "Ratio for discrimination of objects and "
                                "cosmics.",
                                "vimos.Parameters",
                                2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "CosmicsRatio");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "CosmicsRatio");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.order.x",
                                CPL_TYPE_INT,
                                "Order of the bivariate polynomial for "
                                "the X coordinate transformation.",
                                "vimos.Parameters",
                                3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PolyOrderX");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PolyOrderX");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.order.y",
                                CPL_TYPE_INT,
                                "Order of the bivariate polynomial for "
                                "the Y coordinate transformation.",
                                "vimos.Parameters",
                                3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PolyOrderY");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PolyOrderY");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.matching.radius",
                                CPL_TYPE_DOUBLE,
                                "Max distance from expected position "
                                "where a pinhole is searched.",
                                "vimos.Parameters",
                                30.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SearchRadius");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SearchRadius");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.iterations",
                                CPL_TYPE_INT,
                                "Max number of times a fitted distorsion "
                                "model is reused in peak search.",
                                "vimos.Parameters",
                                1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "Iterations");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "Iterations");
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("vimos.Parameters.quality.enable",
                                CPL_TYPE_BOOL,
                                "Compute QC1 parameters.",
                                "vimos.Parameters",
                                TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ComputeQC");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "ComputeQC");
    cpl_parameterlist_append(recipe->parameters, p);


    cx_string_set(path, VIMOS_SEXTRACTOR_PATH);
    cx_string_append(path, "/sex");

    p = cpl_parameter_new_value("vimos.SExtractor.SExtractor",
                                CPL_TYPE_STRING,
                                "Path to the sextractor executable.",
                                "vimos.SExtractor",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.SExtractor");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.SExtractor");
    cpl_parameterlist_append(recipe->parameters, p);


    cx_string_set(path, VIMOS_SEXTRACTOR_CONFIG);
    cx_string_append(path, "/vimos.conv");

    p = cpl_parameter_new_value("vimos.SExtractor.FilterName",
                                CPL_TYPE_STRING,
                                "Name of the file containing the filter "
                                "definition.",
                                "vimos.SExtractor",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.FilterName");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.FilterName");
    cpl_parameterlist_append(recipe->parameters, p);


    cx_string_set(path, VIMOS_SEXTRACTOR_CONFIG);
    cx_string_append(path, "/vimos.nnw");

    p = cpl_parameter_new_value("vimos.SExtractor.StarNnwName",
                                CPL_TYPE_STRING,
                                "Name of the neuronal network weights file.",
                                "vimos.SExtractor",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                            "SExtractor.StarNnwName");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, 
                            "SExtractor.StarNnwName");
    cpl_parameterlist_append(recipe->parameters, p);


    cx_string_set(path, VIMOS_SEXTRACTOR_CONFIG);
    cx_string_append(path, "/masktoccd_1.sex");

    p = cpl_parameter_new_value("vimos.SExtractor.Config1",
                                CPL_TYPE_STRING,
                                "SExtractor configuration file.",
                                "vimos.SExtractor",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.Config1");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.Config1");
    cpl_parameterlist_append(recipe->parameters, p);


    cx_string_set(path, VIMOS_SEXTRACTOR_CONFIG);
    cx_string_append(path, "/masktoccd_2.sex");

    p = cpl_parameter_new_value("vimos.SExtractor.Config2",
                                CPL_TYPE_STRING,
                                "SExtractor configuration file.",
                                "vimos.SExtractor",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "SExtractor.Config2");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "SExtractor.Config2");
    cpl_parameterlist_append(recipe->parameters, p);


    cx_string_set(path, VIMOS_SEXTRACTOR_PATH);
    cx_string_append(path, "/psfex");

    p = cpl_parameter_new_value("vimos.PsfEx.PsfEx",
                                CPL_TYPE_STRING,
                                "Path to the PSF modeling executable.",
                                "vimos.PsfEx",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PsfEx.PsfEx");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PsfEx.PsfEx");
    cpl_parameterlist_append(recipe->parameters, p);


    cx_string_set(path, VIMOS_SEXTRACTOR_CONFIG);
    cx_string_append(path, "/masktoccd.psfex");

    p = cpl_parameter_new_value("vimos.PsfEx.Config",
                                CPL_TYPE_STRING,
                                "PSF modeling configuration file.",
                                "vimos.PsfEx",
                                cx_string_get(path));
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "PsfEx.Config");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CFG, "PsfEx.Config");
    cpl_parameterlist_append(recipe->parameters, p);

    cx_string_delete(path);


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
vmmasktoccd_exec(cpl_plugin *plugin)
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

    if (vmmasktoccd(sof) == EXIT_SUCCESS) {
       
        /*
         * Post process the product frames registered in the set
         * of frames.
         */

        status = vmCplPostProcessFrames(sof, "vmmasktoccd");
        
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
vmmasktoccd_destroy(cpl_plugin *plugin)
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

static double get_orientation(cpl_table *table, double *rms)
{
    float *profile;
    float *triprofile;
    float  maxvalue, minvalue, value;
    double stdv = 90.0;
    double mean = 0.0;
    int    nbin, binsize;
    int    nrows = cpl_table_get_nrow(table);
    int    i, j, maxpos;


    if (cpl_table_has_valid(table, "t")) {

        int null;

        /* At least an average of one pinhole per bin */

        binsize = 180 / (nrows - cpl_table_count_invalid(table, "t"));
        if (binsize > 10) 
            binsize = 10;
        if (binsize == 7 || binsize == 8) // Not an exact divisor of 180
            binsize = 9;
        if (binsize == 0) 
            binsize = 1;
        nbin = 180 / binsize;

        profile = cpl_calloc(nbin, sizeof(float));
        triprofile = cpl_calloc(3 * nbin, sizeof(float));
    
        for (i = 0; i < nrows; i++) {
            value = cpl_table_get_double(table, "t", i, &null) + 90.;
            if (null)
                continue;

            j = floor(value / binsize);
            if (j == nbin)
                j = 0;
            if (j >= 0 && j < nbin)
                profile[j]++;
        }

        /*
         *  Find position of max value in histogram
         */

        maxpos = 0;
        minvalue = maxvalue = profile[maxpos];
        for (i = 1; i < nbin; i++) {
            if (profile[i] > maxvalue) {
                maxvalue = profile[i];
                maxpos = i;
            }
            if (profile[i] < minvalue) {
                minvalue = profile[i];
            }
        }

        /*
         * Build shifted profile
         */

        for (j = 0; j < 3; j++)
            for (i = 0; i < nbin; i++)
                triprofile[i + j * nbin] = profile[i];

        for (i = 0, j = nbin / 2 + maxpos; i < nbin; i++, j++)
            profile[i] = triprofile[j];

        if (findPeak1D(profile, nbin, &value, 3) == VM_TRUE) {
            mean = value + nbin / 2 + maxpos;
            if (mean > nbin)
                mean -= nbin;

            stdv = 0.0;
            value = 0.0;
            for (i = 0; i < nbin; i++) {
                stdv += (i - mean) * (i - mean) * (profile[i] - minvalue);
                value += (profile[i] - minvalue);
            }

            stdv /= value;
            stdv = sqrt(stdv / value);

            stdv *= binsize;
            mean *= binsize;
            mean -= 90.0;
        }
        cpl_free(profile);
        cpl_free(triprofile);
    }

    *rms = stdv;
    return mean;

}


/**
 * @memo
 *   Determine mask to CCD mapping
 *
 * @return EXIT_SUCCESS on EXIT_FAILURE
 *
 * @param sof  Set of frames with references to a pinhole mask
 *             image and a master bias. Optionally a CCD table,
 *             a master dark, or a flat field, may be used in the
 *             corresponding reduction steps of the pinhole image.
 *             If bad pixel cleaning is requested, a CCD table
 *             must be given.
 *
 * @doc
 *   The recipe function computes the expected position of each
 *   mask pinhole on the CCD, using the transformation coefficients
 *   stored in the image header, and then searches for the pinhole
 *   illumination peak within a given search radius. The real position
 *   of each pinhole on the CCD is computed using SExtractor, in psf 
 *   fitting mode, which means running SExtractor AND psfex (another 
 *   nice product from E. Bertin) on the image.
 *   Fit a polynomial and find the distorsion coefficients, in a two stage
 *   process, much like IRAF. First only first order fit, to find offset and
 *   rotation, then a full fit to the residuals. Iterate dropping points
 *   that deviate more than 4*rms of the fit.
 *   Produce the MaskToCcd PAF file with the direct and the inverse
 *   transformation coefficients and update the input mask image header.
 *   On successful termination, the output mask image (with the upgraded
 *   header) and the PAF files are added to the set of frames.
 *
 *   Control options and additional parameters are read from the recipe
 *   configuration database. The recipe function accepts the following
 *   task parameters:
 *
 *   \begin{itemize}
 *
 *     \item PolyOrderX:       Order of the bivariate polynomial for the
 *                             X coordinate transformation. Note that
 *                             the MPS require this to be set to 3.
 *
 *     \item PolyOrderY:       Order of the bivariate polynomial for the
 *                             Y coordinate transformation. Note that
 *                             the MPS require this to be set to 3.
 *
 *     \item SearchRadius:     Max distance from expected position where
 *                             a pinhole is searched,
 *
 *     \item Iterations:       Max number of iterations of the optical
 *                             distorsion model fit, repeated rejecting
 *                             points that deviate more than 4 rms.
 *
 *     \item BiasMethod:       Method for bias removal from pinhole mask
 *                             image. Legal settings are:
 *
 *       \subitem Master:        Prescan and overscan regions are trimmed
 *                               away from the pinhole mask image after
 *                               master bias removal.
 *
 *       \subitem Zmaster:       After master bias subtraction the residual
 *                               signal found in the pinhole mask image 
 *                               overscan regions is modeled and subtracted
 *                               from the image. Next, prescan and overscan
 *                               regions are trimmed away.
 *
 *     \item CleanBadPixel:    Bad pixel correction on pinhole mask image.
 *                             If this option is turned on, the recipe 
 *                             expects a CCD\_TABLE in the input set of 
 *                             frames.
 *
 *     \item CleanCosmic:      Cosmic ray events removal from pinhole
 *                             mask image. If a CCD\_TABLE is found in 
 *                             the input set of frames, bad pixels will 
 *                             not be used in computing the values to 
 *                             replace the cosmic rays events.
 *
 *     \item CosmicThreshold:  Threshold for the detection of cosmic rays,
 *                             having the same meaning of the parameter ns
 *                             of the MIDAS command FILTER/COSMIC: it is
 *                             the number of theoretical noise sigmas above
 *                             smoothed level that must be reached to make
 *                             a pixel a cosmic ray event location candidate.
 *                             This parameter is effective when CleanCosmic
 *                             is set to "true".
 *
 *     \item CosmicRatio:      Critical ratio for discrimination of objects
 *                             and cosmic rays, having the same meaning
 *                             of the parameter rc of the MIDAS command
 *                             FILTER/COSMIC: it is the ratio between the
 *                             peak of a cosmic ray event candidate, and the
 *                             average of the 8 nearby pixels that must be
 *                             reached to identify the candidate as a real
 *                             cosmic ray event.
 *
 *   If any of these task parameters is not set in the recipe configuration
 *   database the recipe function uses the corresponding builtin defaults.
 */

static cxint
vmmasktoccd(PilSetOfFrames *sof)
{

  const char task[] = "vmmasktoccd";

  const char *parameter = "Parameters";
  const char *namePAF = "IMG_mask2ccd";
  const char *pipeNamePAF = "IMG_mask2ccd_";

  char svalue[80];
  char path[PATHNAME_MAX + 1];
  char psfex[PATHNAME_MAX + 1];
  char *argv[SEXTRACTOR_ARGC + 1];
  char *biasMethodTag = NULL;
  char *coeffControlString = NULL;
  char *filterName = NULL;
  char *networkName = NULL;
  char *pafFileName = NULL;
  char *imageName = "pinholes.fits";
  char *parameterName = "pinholes.par";
  char *catalogName = "pinholes.cat";
  char *psfName = "pinholes.psf";

  const char *maskCategory = "IMG_MASK_CCD_CONFIG";
  char filename[PATHNAME_MAX + 1];
  int success;

  unsigned int cleanBadPixel, cleanCosmic;
  unsigned int error;

  int quadrant;
  int argc;
  int fitOrderX, fitOrderY;
  int noIterations;
  int computeQC;
  int biasMethodEntry;
  int *exclude;
  int nSex;
  int noRef, noFit, noReject;
  int noFound;
  int noIter;
  int minSlitNo;
  int i, j, k;
  int dummy;          /* Dummy int argument for fitSurfacePolynomial() */
  int nextCycle = 1;  /* flag for the rejection cycle */

  float diffX, diffY;
  float ron = 0.;
  float level = 0.;
  float gain = 0.;

  double searchRadius;
  double xRef, yRef;
  double cosmicThreshold = 0.0;
  double cosmicRatio = 0.0;
  double xrms = 0.;
  double yrms = 0.;
  double deltaXrms, deltaYrms, oldXrms, oldYrms;
  double *coeffX;
  double *coeffY;
  double *cdX;
  double *cdY;

  const double  pscale = 1.0; // 0.205;
  const char   *punit = "pixel"; // "arcsec";

  time_t timeout = sextGetExecutionTimeLimit();

  FILE *fp = NULL;
  FILE *parameterFile = NULL;

  BiasMethod biasMethod = BIAS_UNDEF;

  PilFrame *maskFrame, *ccdFrame, *biasFrame, *darkFrame, *flatFrame;
  PilFrame *outputFrame;

  VimosImage *tmpImage;
  VimosImage *maskImage = NULL;
  VimosImage *biasImage = NULL;
  VimosImage *flatImage = NULL;
  VimosImage *darkImage = NULL;

  VimosTable *ccdTable = NULL;

  VimosPixel *refPos    = NULL;
  VimosPixel *imaPos    = NULL;
  VimosPixel *corrPos   = NULL;
  VimosPixel *surfaceX  = NULL;
  VimosPixel *surfaceY  = NULL;
  VimosPixel *refList   = NULL;
  VimosPixel *itRefList = NULL;
  VimosPixel *sexPos    = NULL;
  VimosPixel *newPos    = NULL;
/*******
  VimosPixel *ccdList   = NULL;
  VimosPixel *maskList  = NULL;
 *******/


  /*
   * Get task parameters from the recipe database
   */

  /*
   * Determine the frame bias removal method first.
   */

  biasMethodTag = (char *)pilDfsDbGetString(parameter, "BiasMethod");

  if ((biasMethodEntry = strselect(biasMethodTag, biasMethodNames,
                                   nBiasMethods)) < 0) {
      cpl_msg_error(task, "%s: Invalid bias removal method.", biasMethodTag);
      return EXIT_FAILURE;
  }

  biasMethod = biasMethods[biasMethodEntry];


  /*
   * Get degrees of the bivariate polynomial used to fit the mask-CCD
   * transformation. They should be set to 3, 3 to match the MPS request.
   */ 

  fitOrderX = pilDfsDbGetInt(parameter, "PolyOrderX", 3);
  fitOrderY = pilDfsDbGetInt(parameter, "PolyOrderY", 3);

  if (fitOrderX < 2 || fitOrderY < 2) {
    cpl_msg_error(task, "Invalid choice of polynomial order for mask to "
                "CCD transformation: order should be at least 2");
    return EXIT_FAILURE;
  }


  /*
   * Get search radius 
   */

  searchRadius = pilDfsDbGetDouble(parameter, "SearchRadius", 30.);


  /*
   * Get max number of model fit iterations.
   */

  noIterations = pilDfsDbGetInt(parameter, "Iterations", 1);


  /*
   * Check if the bad pixels should be corrected.
   */

  cleanBadPixel = pilDfsDbGetBool(parameter, "CleanBadPixel", 0);


  /*
   * Check if the cosmic rays events should be corrected.
   */

  cleanCosmic = pilDfsDbGetBool(parameter, "CleanCosmic", 0);

  if (cleanCosmic) {
      cosmicThreshold = pilDfsDbGetDouble(parameter, "CosmicThreshold", 4.0);
      cosmicRatio = pilDfsDbGetDouble(parameter, "CosmicRatio", 2.0);
      if (cosmicThreshold < 1.0 || cosmicRatio < 1.0) {
          cpl_msg_error(task, "Invalid cosmic ray filtering parameters: both "
                      "CosmicThreshold and CosmicRatio must be greater than "
                      "one");
          return EXIT_FAILURE;
      }
  }


  /*
   * Check if QC1 parameters should be computed
   */

  computeQC = pilDfsDbGetBool("Parameters", "ComputeQC", 1);


  /*
   * Load input frames
   */

  cpl_msg_info(task, "Loading input frames...");


  /*
   * Get the pinhole mask image
   */

  if ((maskFrame = pilSofLookup(sof, pilTrnGetCategory("GridMaskImage")))) {
      if ((maskImage = openOldFitsFile(pilFrmGetName(maskFrame), 1, 0))) {
          pilFrmSetType(maskFrame, PIL_FRAME_TYPE_RAW);
          closeFitsImage(maskImage, 0);
      }
      else {
          cpl_msg_error(task, "Failure opening pinhole mask frame");
          return EXIT_FAILURE;
      }
  }
  else {
      cpl_msg_error(task, "No input pinhole mask frame found");
      return EXIT_FAILURE;
  }


  /*
   * Get the master bias frame
   */

  error = 1;

  if ((biasFrame = pilSofLookup(sof, pilTrnGetCategory("MasterBias")))) {
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
      deleteImage(maskImage);
      return EXIT_FAILURE;
  }


  /*
   * Recreate bias overscans using as a reference the mask frame
   */

  if ((tmpImage = growOverscans(biasImage, maskImage))) {
      if (biasImage != tmpImage) {
          deleteImage(biasImage);
          biasImage = tmpImage;
      }
  }
  else {
      cpl_msg_error(task, "Failure in growing overscans in master bias");
      deleteImage(maskImage);
      deleteImage(biasImage);
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

  if (VmSubBias(maskImage, biasImage, biasMethod) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot remove bias from pinhole mask frame");
      deleteImage(maskImage);
      deleteImage(biasImage);
      return EXIT_FAILURE;
  }

  deleteImage(biasImage);


  /*
   * Get the (optional) master dark frame
   */

  if ((darkFrame = pilSofLookup(sof, pilTrnGetCategory("MasterDark")))) {
      pilFrmSetType(darkFrame, PIL_FRAME_TYPE_CALIB);
      if ((darkImage = openOldFitsFile(pilFrmGetName(darkFrame), 1, 0))) {
          closeFitsImage(darkImage, 0);
      }
      else {
          cpl_msg_error(task, "Failure opening master dark frame");
          deleteImage(maskImage);
          return EXIT_FAILURE;
      }
  }
  else
      cpl_msg_warning(task, "No master dark in input, skipping dark "
                    "subtraction!");


  /*
   * Now go for the (optional) dark subtraction
   */

  if (darkImage) {
      cpl_msg_info(task, "Dark subtraction...");

      if (VmSubDark(maskImage, darkImage) == EXIT_FAILURE) {
          cpl_msg_error(task, "Cannot subtract dark from pinhole mask image");
          deleteImage(maskImage);
          deleteImage(darkImage);
          return EXIT_FAILURE;
      }
      deleteImage(darkImage);
  }


  /*
   * Get the (optional) master flat field frame
   */

  if ((flatFrame = pilSofLookup(sof, pilTrnGetCategory("ImgMasterSkyFlat")))) {
      pilFrmSetType(flatFrame, PIL_FRAME_TYPE_CALIB);
      if ((flatImage = openOldFitsFile(pilFrmGetName(flatFrame), 1, 0))) {
          closeFitsImage(flatImage, 0);
      }
      else {
          cpl_msg_error(task, "Failure opening master flat field frame");
          deleteImage(maskImage);
          return EXIT_FAILURE;
      }
  }
  else {
      cpl_msg_warning(task, "No master flat field in input, skipping flat "
                    "field correction!");
  }


  /*
   * Now go for the (optional) flat fielding
   */

  if (flatImage) {
      cpl_msg_info(task, "Flat fielding...");

      if ((imageArithLocal(maskImage, flatImage, VM_OPER_DIV))) {
          cpl_msg_error(task, "Flat fielding failure");
          deleteImage(maskImage);
          deleteImage(flatImage);
          return EXIT_FAILURE;
      }
      deleteImage(flatImage);
  }


  /*
   * If bad pixel correction is enabled, a bad pixel table is required
   * in the input set of frames. If no bad pixel table is present this
   * is an error. If cleaning of cosmic rays is enabled, the bad pixel
   * table will be used, if present, to avoid bad pixels while smoothing
   * away cosmic ray events. In this case, a missing bad pixel table is
   * not an error.
   */

  error = 0;

  ccdFrame = pilSofLookup(sof, pilTrnGetCategory("CcdTable"));
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
              cpl_msg_error(task, "Failure opening CCD table");
      }
      else {
          if (cleanBadPixel)
              cpl_msg_error(task, "No CCD table in input: cannot clean "
                          "bad pixels.");
          else
              error = 0;
      }
  }

  if (error) {
      deleteImage(maskImage);
      return EXIT_FAILURE;
  }

  if (cleanCosmic) {

      /*
       * Compute median level, gain, and ron of mask frame for cosmic
       * rays cleaning.
       */

      level = imageMedian(maskImage);
      ron = computeAverageRon(maskImage);
      gain =  getMeanGainFactor(maskImage);

      if (gain < MIN_DIVISOR) {
          cpl_msg_error(task, "Wrong or no gain factor for pinhole mask frame");
          deleteImage(maskImage);
          deleteTable(ccdTable);
          return EXIT_FAILURE;
      }

      /*
       * Cosmic rays cleaning
       */

      if (cleanBadPixel) {
          cpl_msg_info(task, "Cleaning bad pixels and cosmic rays in "
                     "pinhole mask frame");
      }
      else {
          cpl_msg_info(task, "Cleaning bad pixels in pinhole mask frame");
      }

      if (VmCosmicClean(maskImage, ccdTable, cleanBadPixel, level, gain, ron, 
                        cosmicThreshold, cosmicRatio) == EXIT_FAILURE)
          error = 1;

  } else if (cleanBadPixel) {

      /*
       * Just bad pixel cleaning
       */

      cpl_msg_info(task, "Cleaning bad pixels in pinhole mask frame");

      if (EXIT_FAILURE == cleanBadPixels(maskImage, ccdTable, 0))
          error = 1;
  }

  deleteTable(ccdTable);

  if (error) {
      cpl_msg_error(task, "Cannot clean pinhole mask frame");
      deleteImage(maskImage);
      return EXIT_FAILURE;
  }


  /*
   * Prepare VimosPixel list for the task Mask to Ccd. Take reference pinhole 
   * positions in mask coordinates (mm) from the image header.
   */

  if (readIntDescriptor(maskImage->descs, pilTrnGetKeyword("NoSlit"), 
                        &noRef, NULL) == VM_FALSE) {
      cpl_msg_error(task, "Descriptor %s not found",
                  pilTrnGetKeyword("NoSlit"));
      deleteImage(maskImage);
      return EXIT_FAILURE;
  }

  if (noRef == 0) {
      cpl_msg_error(task, "Number of reference apertures is zero");
      deleteImage(maskImage);
      return EXIT_FAILURE;
  }

  minSlitNo = MAX((fitOrderX + 1) * (fitOrderX + 1),
                  (fitOrderY + 1) * (fitOrderY + 1));

  if (noRef < minSlitNo) {
      cpl_msg_error(task, "Fit order is too high: %d pinholes, %d X degree "
                  "polynomial %d Y degree polynomial", noRef, fitOrderX,
                  fitOrderY); 
      deleteImage(maskImage);
      return EXIT_FAILURE;
  }
  
  refList = newPixel(noRef);

  for (i = 0; i < noRef; i++) {
      if (readDoubleDescriptor(maskImage->descs,
                               pilTrnGetKeyword("SlitX", i + 1),
                               &xRef, NULL) == VM_FALSE) {
          cpl_msg_error(task, "Descriptor %s not found",
                      pilTrnGetKeyword("SlitX", i + 1));
          return EXIT_FAILURE;
      }

      if (readDoubleDescriptor(maskImage->descs,
                               pilTrnGetKeyword("SlitY", i + 1),
                               &yRef, NULL) == VM_FALSE) {
          cpl_msg_error(task, "descriptor %s not found",
                      pilTrnGetKeyword("SlitY", i + 1));
          return EXIT_FAILURE;
      }

      refList[i].x = xRef;
      refList[i].y = yRef; 
  }


  /*
   * Convert mask coord (mm) into pixel coords using coefficients from
   * the image header 
   */

  corrPos = newPixel(noRef);
  surfaceX = newPixel(noRef);
  surfaceY = newPixel(noRef);
  itRefList = newPixel(noRef);  

  refPos = MaskToCcd(refList, noRef, maskImage->descs);

  if (!refPos) {
      cpl_msg_error(task, "Failure in MaskToCcd");
      deleteImage(maskImage);
      return EXIT_FAILURE;
  }

  for (i = 0; i < noRef; i++) {
      cpl_msg_debug(task, "refPos: %f, %f", (float)refPos[i].x,
                  (float)refPos[i].y);
  }


  /*
   * Find position of pinholes in mask image   
   */ 

  cpl_msg_info(task, "Searching image for pinholes ...");

  if (!(networkName = (char *)sextGetStarNnwName())) {
      cpl_msg_error(task, "Cannot retrieve SExtractor neuronal network "
                  "setup file!");

      deleteImage(maskImage);

      return 0;
  }
  else {
      if (access(networkName, F_OK | R_OK)) {
          cpl_msg_error(task, "Cannot access SExtractor neuronal network "
                      "setup file %s!", networkName);

          deleteImage(maskImage);

          return 0;
      }
  }


  if (!(filterName = (char *)sextGetFilterName())) {
      cpl_msg_error(task, "Cannot retrieve SExtractor filter setup file!");

      deleteImage(maskImage);

      return 0;
  }
  else {
      if (access(filterName, F_OK | R_OK)) {
          cpl_msg_error(task, "Cannot access SExtractor neuronal network "
                      "setup file %s!", filterName);

          deleteImage(maskImage);

          return 0;
      }
  }

  if (!createFitsImage((char *)imageName, maskImage, "UNKNOWN")) {
      cpl_msg_error(task, "Cannot create temporary SExtractor input image "
                  "file!");

      deleteImage(maskImage);

      return 0;
  }


  if (!(parameterFile = fopen(parameterName, "w"))) {
      cpl_msg_error(task, "Cannot create temporary file!");

      deleteImage(maskImage);
      remove(imageName);

      return 0;
  }

  if (sextSaveParameters(parameterFile, masktoccd_1) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot write SExtractor parameter file!");

      fclose(parameterFile);
      remove(parameterName);

      deleteImage(maskImage);
      remove(imageName);

      return 0;
  }

  fclose(parameterFile);


  cpl_msg_info(task, "Initial pinhole detection ...");

  sextGetFileName(path, pilDfsDbGetString("SExtractor", "Config1"),
                  PATHNAME_MAX);

  argc = 0;
  argv[argc++] = (char *)sextGetSextractorPath();
  argv[argc++] = imageName;
  argv[argc++] = "-c";
  argv[argc++] = path;
  argv[argc++] = "-PARAMETERS_NAME";
  argv[argc++] = parameterName;
  argv[argc++] = "-FILTER_NAME";
  argv[argc++] = filterName;
  argv[argc++] = "-STARNNW_NAME";
  argv[argc++] = networkName;
  argv[argc++] = "-CATALOG_NAME";
  argv[argc++] = catalogName;
  argv[argc] = NULL;

  if (pilTaskExecWait(argc, (const char *const *)argv, timeout)) {
      cpl_msg_error(task, "Running SExtractor failed!");

      deleteImage(maskImage);
      remove(imageName);

      remove(parameterName);
      remove(catalogName);

      return EXIT_FAILURE;
  }

  remove(parameterName);


  sextGetFileName(psfex, pilDfsDbGetString("PsfEx", "PsfEx"),
                  PATHNAME_MAX);

  sextGetFileName(path, pilDfsDbGetString("PsfEx", "Config"),
                  PATHNAME_MAX);

  argc = 0;
  argv[argc++] = psfex;
  argv[argc++] = catalogName;
  argv[argc++] = "-c";
  argv[argc++] = path;
  argv[argc++] = "-PSF_NAME";
  argv[argc++] = psfName;
  argv[argc] = NULL;

  if (pilTaskExecWait(argc, (const char *const *)argv, timeout)) {
      cpl_msg_error(task, "Running PSFex failed!");

      deleteImage(maskImage);
      remove(imageName);

      remove(catalogName);
      remove(psfName);

      return EXIT_FAILURE;
  }
  
  remove(catalogName);

  
  cpl_msg_info(task, "Final pinhole detection (refining positions) ...");

  if (!(parameterFile = fopen(parameterName, "w"))) {
      cpl_msg_error(task, "Cannot create temporary file!");

      deleteImage(maskImage);
      remove(imageName);

      remove(psfName);

      return 0;
  }

  if (sextSaveParameters(parameterFile, masktoccd_2) == EXIT_FAILURE) {
      cpl_msg_error(task, "Cannot write SExtractor parameter file!");

      fclose(parameterFile);
      remove(parameterName);

      deleteImage(maskImage);
      remove(imageName);

      remove(psfName);


      return 0;
  }

  fclose(parameterFile);

  sextGetFileName(path, pilDfsDbGetString("SExtractor", "Config2"),
                  PATHNAME_MAX);
  argc = 0;
  argv[argc++] = (char *)sextGetSextractorPath();
  argv[argc++] = imageName;
  argv[argc++] = "-c";
  argv[argc++] = path;
  argv[argc++] = "-PARAMETERS_NAME";
  argv[argc++] = parameterName;
  argv[argc++] = "-FILTER_NAME";
  argv[argc++] = filterName;
  argv[argc++] = "-STARNNW_NAME";
  argv[argc++] = networkName;
  argv[argc++] = "-CATALOG_NAME";
  argv[argc++] = catalogName;
  argv[argc++] = "-PSF_NAME";
  argv[argc++] = psfName;
  argv[argc] = NULL;

  if (pilTaskExecWait(argc, (const char *const *)argv, timeout)) {
      cpl_msg_error(task, "Running SExtractor failed!");

      deleteImage(maskImage);
      remove(imageName);

      remove(parameterName);
      remove(catalogName);
      remove(psfName);

      return EXIT_FAILURE;
  }

  remove(parameterName);
  remove(imageName);
  remove(psfName);


  if (!(fp = fopen(catalogName, "r"))) {
      cpl_msg_warning(task, "Could not open pinholes list `%s'",
                    catalogName);

      deleteImage(maskImage);
      remove(catalogName);

      return EXIT_FAILURE;
  }
  else {
      
      char line[MAX_LINE_LENGTH];

      nSex = 0;
      while (fgets(line, MAX_LINE_LENGTH, fp)) {
          nSex++;
      }  
      rewind(fp);

      sexPos = newPixel(nSex);

      for (i = 0; i < nSex; i++) {

          char c;
          float x, y;

          fscanf(fp, "%*d %f %f %f", &x, &y, &sexPos[i].i); 

          sexPos[i].x = x;
          sexPos[i].y = y;

          while((c = getc(fp)) != '\n' && c != EOF) {
              ;
          }
      }

      fclose(fp);
//      remove(catalogName);
  }

  noFound = 0;
  for (i = 0; i < noRef; i++) {

      imaPos = matchSexPosition(sexPos, nSex, &(refPos[i]), searchRadius);

      if (imaPos != NULL) {
          corrPos[noFound].x = imaPos->x;
          corrPos[noFound].y = imaPos->y;
          itRefList[noFound].x = refList[i].x;
          itRefList[noFound].y = refList[i].y;
          noFound++;
          deletePixel(imaPos);
      }
      else {
          cpl_msg_debug(task, "Pinhole %d not found: refposX = %f, "
                      "refposY = %f", i, refPos[i].x, refPos[i].y); 
      }
  }

  cpl_msg_info(task, "Found %d pinholes", noFound);

  if (noFound < minSlitNo) {
      cpl_msg_info(task, "Too few pinholes found: cannot fit the requested "
                 "polynomial"); 
      return EXIT_FAILURE;
  }

  /*
   * Prepare pixel lists to fit and estimate X and Y coefficients of the 
   * mask to CCD transformation
   */

  for (i = 0; i < noFound; i++) {
      surfaceX[i].x = itRefList[i].x;
      surfaceX[i].y = itRefList[i].y;
      surfaceX[i].i = corrPos[i].x;
      surfaceY[i].x = itRefList[i].x;
      surfaceY[i].y = itRefList[i].y;
      surfaceY[i].i = corrPos[i].y;
  }

  /* 
   * Here again a new part, where we fit with a rejection on the points
   * that deviate from the fit (but we do not touch the matching !!!) -- MS 
   */

  noIter = 0;
  noFit = noFound;
  oldXrms = 99.99;
  oldYrms = 99.99;
  exclude = intVector(0, noFound);

  for (i = 0; i < noFound; i++)
      exclude[i] = 0;

  coeffControlString = createVimosCtrlStr(fitOrderX, fitOrderY); 

  do {

      /* 
       * We try to follow IRAF: first we fit offset and rotation only... 
       */

      cdX = fitSurfacePolynomial(surfaceX, noFit, "(0,0) (1,0) (0,1)", 1, 
                                 &dummy, &xrms);
      cdY = fitSurfacePolynomial(surfaceY, noFit, "(0,0) (1,0) (0,1)", 1, 
                                 &dummy, &yrms);

      /* 
       * ... and then we do a full fit of the residuals. -- MS
       */

      j = 0;

      for (i = 0; i < noFound; i++) {
          if (!exclude[i]) {
              surfaceX[j].i = corrPos[i].x -
                  (cdX[0] + cdX[1] * itRefList[i].x + cdX[2] * itRefList[i].y);
              surfaceY[j].i = corrPos[i].y -
                  (cdY[0] + cdY[1] * itRefList[i].x + cdY[2] * itRefList[i].y);
              j++;
          }
      }

      coeffX = fitSurfacePolynomial(surfaceX, noFit, coeffControlString, 
                                    fitOrderX * 2, &dummy, &xrms);
      coeffY = fitSurfacePolynomial(surfaceY, noFit, coeffControlString, 
                                    fitOrderY * 2, &dummy, &yrms);

      sprintf(svalue, "%.12G", cdX[0]);
      writeStringDescriptor((&(maskImage->descs)),
                            pilTrnGetKeyword("MaskCcdX0"), svalue, "");
      sprintf(svalue, "%.8G", cdX[1]);
      writeStringDescriptor((&(maskImage->descs)),
                            pilTrnGetKeyword("MaskCcdXX"), svalue, "");
      sprintf(svalue, "%.8G", cdX[2]);
      writeStringDescriptor((&(maskImage->descs)),
                            pilTrnGetKeyword("MaskCcdXY"), svalue, "");
      sprintf(svalue, "%.12G", cdY[0]);
      writeStringDescriptor((&(maskImage->descs)),
                            pilTrnGetKeyword("MaskCcdY0"), svalue, "");
      sprintf(svalue, "%.8G", cdY[1]);
      writeStringDescriptor((&(maskImage->descs)),
                            pilTrnGetKeyword("MaskCcdYX"), svalue, "");
      sprintf(svalue, "%.8G", cdY[2]);
      writeStringDescriptor((&(maskImage->descs)),
                            pilTrnGetKeyword("MaskCcdYY"), svalue, "");

      k = 0;
      for (i = 0; i <= fitOrderX; i++) {
          k=i;
          for (j = 0; j <= fitOrderX; j++) {
              sprintf(svalue, "%.8G", coeffX[k]);
              writeStringDescriptor((&(maskImage->descs)),
                                    pilTrnGetKeyword("MaskCcdX", i, j),
                                    svalue, "");
              k = k + fitOrderX + 1;
          }
      }
    
      k = 0;
      for (i = 0; i <= fitOrderY; i++) {
          k=i;
          for (j = 0; j <= fitOrderY; j++) {
              sprintf(svalue, "%.8G", coeffY[k]);
              writeStringDescriptor((&(maskImage->descs)),
                                    pilTrnGetKeyword("MaskCcdY", i, j),
                                    svalue, "");
              k = k + fitOrderY + 1;
          }
      }                        

      /* do we need to iterate some more ?? */

      if (noIter < noIterations) {
          /* 
           * Measure differences between new and old positions, and reject
           * points that are 4-sigma outliers -- MS
           */
          deletePixel(newPos);
          newPos = newPixel(noFit);
          newPos = MaskToCcd(itRefList, noFit, maskImage->descs);
          for (i = 0; i < noFound; i++) {
              diffX = fabs(newPos[i].x - corrPos[i].x);
              diffY = fabs(newPos[i].y - corrPos[i].y);
              if (diffX > 4. * sqrt(xrms) || diffY > 4. * sqrt(yrms))
                  exclude[i] = 1;
          }
          j = 0;
          for (i = 0; i < noFound; i++) {
              if (!exclude[i]) {
                  surfaceX[j].x = itRefList[i].x;
                  surfaceX[j].y = itRefList[i].y;
                  surfaceX[j].i = corrPos[i].x;
                  surfaceY[j].x = itRefList[i].x;
                  surfaceY[j].y = itRefList[i].y;
                  surfaceY[j].i = corrPos[i].y;
                  j++;
              }
          }
      
          deltaXrms = oldXrms - sqrt(xrms);
          deltaYrms = oldYrms - sqrt(yrms);
          oldXrms = sqrt(xrms);
          oldYrms = sqrt(yrms);
          noReject = noFit - j;
          noFit = j;
          noIter++;

          /* if no points are rejected we don't need any more iterations */
          if (noReject == 0) {
              nextCycle = 0;
          }
          /* if the rms doesn't decrease any more, we can stop with a last
             iteration, to make sure we use the latest rejections */
          else
              if (deltaXrms < (0.1 * oldXrms) && deltaYrms < (0.1 * oldYrms)) {
                  noIter = noIterations;
              }
      }
      /* otherwise we just stop the iterations cycle */
      else {
          nextCycle = 0;
      }
  } while (nextCycle);

  cpl_msg_info (task, "Number of fitted points %d", noFit);
  writeDoubleDescriptor(&(maskImage->descs), pilTrnGetKeyword("MaskCcdXrms"),
                        sqrt(xrms), "");
  writeDoubleDescriptor(&(maskImage->descs), pilTrnGetKeyword("MaskCcdYrms"),
                        sqrt(yrms), "");
    
  writeIntDescriptor(&(maskImage->descs), pilTrnGetKeyword("MaskCcdXord"),
                     fitOrderX, "");
  writeIntDescriptor(&(maskImage->descs), pilTrnGetKeyword("MaskCcdYord"),
                     fitOrderY, "");
  
  if ( (xrms >= 1)  || (yrms >= 1) ) {
      cpl_msg_debug(task, "The error in the slit position is > 1 pixel: "
                  "The slit coordinates on the CCD are determined with "
                  "accuracy XRMS = %f, YRMS = %f", sqrt(xrms), sqrt(yrms));
  }


  if (computeQC) {

    cpl_msg_info(task, "Computing QC1 parameters...");

    if (pilQcGroupStart() == EXIT_SUCCESS) {

      /*
       * Write PRO.CATG, ARCFILE, TPL ID, and OCS CON QUAD, to QC1 group.
       */

      readIntDescriptor(maskImage->descs, pilTrnGetKeyword("Quadrant"),
                        &quadrant, NULL);

      pilQcWriteString("PRO.CATG", maskCategory, "Product category");

      qcCopyValue(maskImage->descs, pilTrnGetKeyword("ArchiveFile"),
                  NULL, "Archive File Name");

      qcCopyValue(maskImage->descs, pilTrnGetKeyword("TplId"),
                  NULL, "Template signature ID");

      qcCopyValue(maskImage->descs, pilTrnGetKeyword("Quadrant"),
                  NULL, "Quadrant");

      qcCopyValue(maskImage->descs, pilTrnGetKeyword("FilterName", quadrant),
                  NULL, "Filter name");

      qcCopyValue(maskImage->descs, pilTrnGetKeyword("MaskId", quadrant),
                  NULL, "Mask identification");

      qcCopyValue(maskImage->descs, 
                  pilTrnGetKeyword("BeamTemperature", quadrant), 
                  NULL, "Beam temperature");

      qcWriteValueDouble(maskImage->descs, cdX[0], "QC.MASK.CCD.X0", "pixel",
                         "Mask2ccd coefficient");

      qcWriteValueDouble(maskImage->descs, cdX[1], "QC.MASK.CCD.XX", "pixel/mm",
                         "Mask2ccd coefficient");

      qcWriteValueDouble(maskImage->descs, cdX[2], "QC.MASK.CCD.XY", "pixel/mm",
                         "Mask2ccd coefficient");

      qcWriteValueDouble(maskImage->descs, cdY[0], "QC.MASK.CCD.Y0", "pixel",
                         "Mask2ccd coefficient");

      qcWriteValueDouble(maskImage->descs, cdY[1], "QC.MASK.CCD.YX", "pixel/mm",
                         "Mask2ccd coefficient");

      qcWriteValueDouble(maskImage->descs, cdY[2], "QC.MASK.CCD.YY", "pixel/mm",
                         "Mask2ccd coefficient");


      if (1) {

          /* Deal with new QCs */

          cpl_table *qctable;
          double     rms;
          float      x, y, a, b, e, t;
          char       c;
          int        nsel;


          fp = fopen(catalogName, "r");

          qctable = cpl_table_new(nSex);
          cpl_table_new_column(qctable, "x", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "y", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "a", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "b", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "e", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "t", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "s", CPL_TYPE_INT);

          for (i = 0; i < nSex; i++) {

              fscanf(fp, "%*d %*f %*f %*f %f %f %*f %*f %*d %f %*f %f %f %f", 
                     &x, &y, &e, &t, &a, &b);

              while((c = getc(fp)) != '\n' && c != EOF) { ; }

              cpl_table_set_double(qctable, "x", i, x);
              cpl_table_set_double(qctable, "y", i, y);
              cpl_table_set_double(qctable, "a", i, 2*a);
              cpl_table_set_double(qctable, "b", i, 2*b);
              cpl_table_set_double(qctable, "e", i, e);
              cpl_table_set_double(qctable, "t", i, t);

          }

          fclose(fp);

          for (i = 0; i < noFound; i++) {

              if (exclude[i])
                  continue;

              cpl_table_and_selected_double(qctable, "x", CPL_LESS_THAN,
                                            corrPos[i].x + 2.0);
              cpl_table_and_selected_double(qctable, "x", CPL_GREATER_THAN,
                                            corrPos[i].x - 2.0);
              cpl_table_and_selected_double(qctable, "y", CPL_LESS_THAN,
                                            corrPos[i].y + 2.0);
              nsel = 
              cpl_table_and_selected_double(qctable, "y", CPL_GREATER_THAN,
                                            corrPos[i].y - 2.0);

              if (nsel == 1) {
                  cpl_array *array = cpl_table_where_selected(qctable);
                  cpl_size     pos = cpl_array_get_cplsize(array, 0, NULL);

                  cpl_table_set_int(qctable, "s", (int)pos, 1);
                  cpl_array_delete(array);
              }

              cpl_table_select_all(qctable);

          }

          cpl_table_and_selected_int(qctable, "s", CPL_EQUAL_TO, 1);
          cpl_table_not_selected(qctable);
          cpl_table_erase_selected(qctable);

//          printf("QC PINHOLE COUNT = %d\n", noFit);

          qcWriteValueInt(maskImage->descs, noFit, "QC.PINHOLE.COUNT", NULL,
                          "Number of used pinholes");

          cpl_msg_info(task, "QC IMAGE QUALITY MEAN = %f %s", 
                 pscale * cpl_table_get_column_mean(qctable, "a"), punit);

          qcWriteValueDouble(maskImage->descs, 
                             pscale * cpl_table_get_column_mean(qctable, "a"), 
                             "QC.IMAGE.QUALITY.MEAN", punit,
                             "Mean image quality");

          cpl_msg_info(task, "QC IMAGE QUALITY RMS = %f %s", 
                 pscale * cpl_table_get_column_stdev(qctable, "a"), punit);

          qcWriteValueDouble(maskImage->descs, 
                             pscale * cpl_table_get_column_stdev(qctable, "a"),
                             "QC.IMAGE.QUALITY.RMS", punit,
                             "Image quality RMS");

          cpl_msg_info(task, "QC ELLIPTICITY MEAN = %f", 
                 cpl_table_get_column_mean(qctable, "e"));

          qcWriteValueDouble(maskImage->descs, 
                             cpl_table_get_column_mean(qctable, "e"), 
                             "QC.ELLIPTICITY.MEAN", NULL,
                             "Mean image ellipticity");

          cpl_msg_info(task, "QC ELLIPTICITY RMS = %f", 
                 cpl_table_get_column_stdev(qctable, "e"));

          qcWriteValueDouble(maskImage->descs, 
                             cpl_table_get_column_stdev(qctable, "e"), 
                             "QC.ELLIPTICITY.RMS", NULL,
                             "Image ellipticity RMS");

          cpl_msg_info(task, "QC ORIENTATION MEAN = %f degree", 
                            get_orientation(qctable, &rms));

          qcWriteValueDouble(maskImage->descs, 
                             get_orientation(qctable, &rms),
                             "QC.ORIENTATION.MEAN", "degree",
                             "Mean image orientation");

          cpl_msg_info(task, "QC ORIENTATION RMS = %f degree", rms);

          qcWriteValueDouble(maskImage->descs, rms,
                             "QC.ORIENTATION.RMS", "degree",
                             "Image orientation RMS");

          if (1) {
              cpl_table *sltable;
              char       descName[VM_DESC_LENGTH];
              int        xStart, xEnd, xPos, xStep;
              int        yStart, yEnd, yPos, yStep;
              int        cell;

              yStart = 80;
              yStep = 760;
              cell = 1;

              for (yPos = 0; yPos < 3; yPos++, yStart += yStep) {
                  yEnd = yStart + yStep;
                  xStart = 73;
                  xStep = 634;

                  for (xPos = 0; xPos < 3; xPos++, xStart += xStep, cell++) {
                      xEnd = xStart + xStep;

                      cpl_table_select_all(qctable);

                      cpl_table_and_selected_double(qctable, "x", 
                                                    CPL_LESS_THAN, xEnd);
                      cpl_table_and_selected_double(qctable, "x", 
                                                    CPL_GREATER_THAN, xStart);
                      cpl_table_and_selected_double(qctable, "y", 
                                                    CPL_LESS_THAN, yEnd);
                      nsel = 
                      cpl_table_and_selected_double(qctable, "y", 
                                                    CPL_GREATER_THAN, yStart);

                      if (nsel == 0) {

                          cpl_msg_warning(task, "No pinholes in cell %d\n", 
                                          cell);
                          continue;
                      }

                      cpl_msg_info(task, "-----------------");

                      cpl_msg_info(task, "%d pinholes in cell %d", nsel, cell);

                      sltable = cpl_table_extract_selected(qctable);

                      sprintf(descName, "QC.IMAGE.QUAL%d.MEAN", cell);

                      cpl_msg_info(task, "%s = %f %s", descName,
                           pscale * cpl_table_get_column_mean(sltable, "a"),
                           punit);

                      qcWriteValueDouble(maskImage->descs,
                           pscale * cpl_table_get_column_mean(sltable, "a"),
                           descName, punit, "Mean image quality in cell");

                      sprintf(descName, "QC.IMAGE.QUAL%d.RMS", cell);

                      cpl_msg_info(task, "%s = %f %s", descName,
                           pscale * cpl_table_get_column_stdev(sltable, "a"),
                           punit);

                      qcWriteValueDouble(maskImage->descs,
                           pscale * cpl_table_get_column_stdev(sltable, "a"),
                           descName, punit, "Image quality RMS in cell");

                      sprintf(descName, "QC.ELLIPTIC%d.MEAN", cell);

                      cpl_msg_info(task, "%s = %f", descName,
                             cpl_table_get_column_mean(sltable, "e"));

                      qcWriteValueDouble(maskImage->descs,
                                         cpl_table_get_column_mean(sltable, 
                                                                   "e"),
                                         descName, NULL,
                                         "Mean ellipticity in cell");

                      sprintf(descName, "QC.ELLIPTIC%d.RMS", cell);

                      cpl_msg_info(task, "%s = %f", descName,
                             cpl_table_get_column_stdev(sltable, "e"));

                      qcWriteValueDouble(maskImage->descs,
                                         cpl_table_get_column_stdev(sltable, 
                                                                    "e"),
                                         descName, NULL,
                                         "Ellipticity RMS in cell");

                      sprintf(descName, "QC.ORIENTAT%d.MEAN", cell);

                      cpl_msg_info(task, "%s = %f degree", descName,
                             get_orientation(sltable, &rms));

                      qcWriteValueDouble(maskImage->descs,
                                         get_orientation(sltable, &rms),
                                         descName, "degree",
                                         "Mean orientation in cell");

                      sprintf(descName, "QC.ORIENTAT%d.RMS", cell);

                      cpl_msg_info(task, "%s = %f degree", descName, rms);

                      qcWriteValueDouble(maskImage->descs, rms,
                                         descName, "degree",
                                         "Orientation RMS in cell");

                      cpl_table_delete(sltable);
                  }
              }
          }

          // cpl_table_save(qctable, NULL, NULL, "catalog.fits", CPL_IO_CREATE);
          cpl_table_delete(qctable);

      }

      if (0) {

          /*
           * This part of the code is valid in case you
           * wanted to use the cpl_fit_gaussian(), if it
           * were not slow and returned NaN in all places.
           */

          /* Deal with new QCs */

          cpl_image *qcimage;
          cpl_array *qcarray;
          cpl_table *qctable;
          double     major, minor, angle, rms;


          /*
           *  maskImage is cast to a CPL type, because this part
           * of the code is based on the CPL.
           */

          qcimage = cpl_image_wrap_float(maskImage->xlen, maskImage->ylen, 
                                         maskImage->data);

          qctable = cpl_table_new(noFit);
          cpl_table_new_column(qctable, "x", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "y", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "a", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "b", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "e", CPL_TYPE_DOUBLE);
          cpl_table_new_column(qctable, "t", CPL_TYPE_DOUBLE);

          qcarray = cpl_array_new(7, CPL_TYPE_DOUBLE);

          for (j = 0, i = 0; i < noFound; i++) {

              if (exclude[i])
                  continue;

              cpl_array_set_double(qcarray, 3, corrPos[i].x);
              cpl_array_set_double(qcarray, 4, corrPos[i].y);

              if (cpl_fit_image_gaussian(qcimage, NULL, 
                                         corrPos[i].x, corrPos[i].y,
                                         5, 5, qcarray, 
                                         NULL, NULL, NULL, NULL, NULL,
                                         &major, &minor, &angle, NULL))
                  continue;

              cpl_table_set_double(qctable, "x", j, corrPos[i].x);
              cpl_table_set_double(qctable, "y", j, corrPos[i].y);
              cpl_table_set_double(qctable, "a", j, major);
              cpl_table_set_double(qctable, "b", j, minor);
              cpl_table_set_double(qctable, "e", j, 1 - minor / major);
              cpl_table_set_double(qctable, "t", j, angle * 180 / PI_NUMB);

              j++;
              cpl_array_fill_window_invalid(qcarray, 0, 7);
          }

          cpl_array_delete(qcarray);
          cpl_image_unwrap(qcimage);

          printf("QC PINHOLE COUNT = %d\n", noFit);

          printf("QC IMAGE QUALITY MEAN = %f\n", 
                 cpl_table_get_column_mean(qctable, "a"));

          printf("QC IMAGE QUALITY RMS = %f\n", 
                 cpl_table_get_column_stdev(qctable, "a"));

          printf("QC ELLIPTICITY MEAN = %f\n", 
                 cpl_table_get_column_mean(qctable, "e"));

          printf("QC ELLIPTICITY RMS = %f\n", 
                 cpl_table_get_column_stdev(qctable, "e"));

          printf("QC ORIENTATION MEAN = %f\n", get_orientation(qctable, &rms));

          printf("QC ORIENTATION RMS = %f\n", rms);

          cpl_table_save(qctable, NULL, NULL, "catalog.fits", CPL_IO_CREATE);
          cpl_table_delete(qctable);

      }

      if (pilQcGroupEnd() == EXIT_FAILURE)
        cpl_msg_error(task, "Cannot create QC1 PAF file! (Ignored)");

    }
    else
      cpl_msg_error(task, "Cannot open QC1 logging! (Ignored)");

  } /* End of QC1 computation. */

  
  /* 
   * Compute inverse transformation coefficients SAME MODIFICATION AS FOR
   * THE DIRECT FIT. No rejection needed any more, of course... -- MS
   */

  j = 0;

  for (i = 0; i < noFound; i++) {
      if (!exclude[i]) {
          surfaceX[j].x = corrPos[i].x;
          surfaceX[j].y = corrPos[i].y;
          surfaceX[j].i = itRefList[i].x;
          surfaceY[j].x = corrPos[i].x;
          surfaceY[j].y = corrPos[i].y;
          surfaceY[j].i = itRefList[i].y;
          j++;
      }
  }

  cdX = fitSurfacePolynomial(surfaceX, noFit, "(0,0) (1,0) (0,1)", 
                             1, &dummy, &xrms);
  cdY = fitSurfacePolynomial(surfaceY, noFit, "(0,0) (1,0) (0,1)", 
                             1, &dummy, &yrms);

  j = 0;

  for (i = 0; i < noFound; i++) {
      if (!exclude[i]) {
          surfaceX[j].i = itRefList[i].x - (cdX[0] + cdX[1] * corrPos[i].x 
                                            + cdX[2] * corrPos[i].y);
          surfaceY[j].i = itRefList[i].y - (cdY[0] + cdY[1] * corrPos[i].x
                                            + cdY[2] * corrPos[i].y);
          j++;
      }
  }

  coeffX = fitSurfacePolynomial(surfaceX, noFit, coeffControlString, 
                                fitOrderX * 2, &dummy, &xrms);
  coeffY = fitSurfacePolynomial(surfaceY, noFit, coeffControlString, 
                                fitOrderY * 2, &dummy, &yrms);
  sprintf(svalue, "%.12G", cdX[0]);
  writeStringDescriptor((&(maskImage->descs)),
                        pilTrnGetKeyword("CcdMaskX0"), svalue, "");
  sprintf(svalue, "%.8G", cdX[1]);
  writeStringDescriptor((&(maskImage->descs)),
                        pilTrnGetKeyword("CcdMaskXX"), svalue, "");
  sprintf(svalue, "%.8G", cdX[2]);
  writeStringDescriptor((&(maskImage->descs)),
                        pilTrnGetKeyword("CcdMaskXY"), svalue, "");
  sprintf(svalue, "%.12G", cdY[0]);
  writeStringDescriptor((&(maskImage->descs)),
                        pilTrnGetKeyword("CcdMaskY0"), svalue, "");
  sprintf(svalue, "%.8G", cdY[1]);
  writeStringDescriptor((&(maskImage->descs)),
                        pilTrnGetKeyword("CcdMaskYX"), svalue, "");
  sprintf(svalue, "%.8G", cdY[2]);
  writeStringDescriptor((&(maskImage->descs)),
                        pilTrnGetKeyword("CcdMaskYY"), svalue, "");

  cpl_free(coeffControlString);
    
  k = 0;
  for (i = 0; i <= fitOrderX; i++) {
      k=i;
      for (j = 0; j <= fitOrderX; j++) {
          sprintf(svalue, "%.8G", coeffX[k]);
          writeStringDescriptor((&(maskImage->descs)),
                                pilTrnGetKeyword("CcdMaskX", i, j),
                                svalue, "");
          k = k + fitOrderX + 1;
      }
  }
  
  k = 0;
  for (i = 0; i <= fitOrderY; i++) {
      k=i;
      for (j = 0; j <= fitOrderY; j++) {
          sprintf(svalue, "%.8G", coeffY[k]);
          writeStringDescriptor((&(maskImage->descs)),
                                pilTrnGetKeyword("CcdMaskY", i, j),
                                svalue, "");
          k = k + fitOrderY + 1;
      }
  } 
 
  writeDoubleDescriptor(&maskImage->descs, pilTrnGetKeyword("CcdMaskXrms"),
                        sqrt(xrms), ""); 
  writeDoubleDescriptor(&maskImage->descs, pilTrnGetKeyword("CcdMaskYrms"),
                        sqrt(yrms), "");
  
  writeIntDescriptor(&(maskImage->descs), pilTrnGetKeyword("CcdMaskXord"),
                     fitOrderX,  "");
  writeIntDescriptor(&(maskImage->descs), pilTrnGetKeyword("CcdMaskYord"),
                     fitOrderY,"");

  /*
   * Insert check on ALL reference positions (also the rejected ones)
   * transforming mask coordinates to CCD coordinates, and then back
   * CCD coordinates to mask coordinates, to see what is the real
   * accuracy of the models.
   */

/******* Restore this for check on models self-consistency.
 ******* Do not forget to set the temperature scale factor
 ******* to 1.0 in CcdToMask()    (C.Izzo)

  ccdList = MaskToCcd(refList, noRef, maskImage->descs);
  maskList = CcdToMask(ccdList, noRef, maskImage->descs);
   
  for (i = 0; i < noRef; i++)
    printf("%.5f %.5f %.5f %.5f\n", 
           refList[i].x, refList[i].y, maskList[i].x, maskList[i].y);

 *******
 *******/

  createMaskCcdPAF(maskImage->descs, namePAF, &pafFileName);
  
  if (createMaskCcdPAF(maskImage->descs, pipeNamePAF, &pafFileName) 
      == EXIT_SUCCESS) {
      outputFrame = newPilFrame(pafFileName, pilTrnGetCategory("PAFCategory"));

      pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);  
      pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_PAF);
      pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_PRIMARY);
      pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

      pilSofInsert(sof, outputFrame);
  }
  else {
      cpl_msg_error(task, "Cannot create local product file %s!", pafFileName);
      deleteImage(maskImage);
      return EXIT_FAILURE;
  }


  /*
   * Remove data, just empty primary array is wanted by DFO (24.11.2008 C.Izzo)
   */

  cpl_free(maskImage->data);
  maskImage->data = NULL;
  maskImage->xlen = 0;
  maskImage->ylen = 0;

  vmstrlower(strcpy(filename, maskCategory));
  strcat(filename, ".fits");

  success = 0;
  if (openNewFitsImage(filename, maskImage)) {
    removeDescriptor(&(maskImage->descs), "NAXIS1");
    removeDescriptor(&(maskImage->descs), "NAXIS2");
    writeIntDescriptor(&(maskImage->descs), "NAXIS", 0, "Empty");
    copyFromHeaderToHeader(maskImage->descs, "ESO INS FOCU1 TEMP",
                           &(maskImage->descs), "ESO PRO CCD MASK TEMP");
    if (writeDescsToFitsImage(maskImage->descs, maskImage) == VM_TRUE) {
      if (closeFitsImage(maskImage, 0) == VM_TRUE) {
        cpl_msg_debug(task, "Image %s (%s) created", filename, maskCategory);
        success = 1;
      }
    }
  }

//  if (createFitsImage(filename, maskImage, maskCategory)) {

  if (success) {
    outputFrame = newPilFrame(filename, maskCategory);

    pilFrmSetType(outputFrame, PIL_FRAME_TYPE_PRODUCT);
    pilFrmSetFormat(outputFrame, PIL_FRAME_FORMAT_IMAGE);
    pilFrmSetProductLevel(outputFrame, PIL_PRODUCT_LEVEL_SECONDARY);
    pilFrmSetProductType(outputFrame, PIL_PRODUCT_TYPE_REDUCED);

    pilSofInsert(sof, outputFrame);
  }
  else {
    cpl_msg_error(task, "Cannot create local product file %s!", filename);
    deleteImage(maskImage);
    return EXIT_FAILURE;
  }


  /*
   *  Cleanup
   */
  
  deleteImage(maskImage);
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
                    "vmmasktoccd",
    "Determine CCD to mask transformation and its inverse.",
    "This recipe is used to determine the CCD to Mask transformation and\n"
    "its inverse. SExtractor v2.1.6 is used for determining the positions\n"
    "of the spotlights on a direct imaging exposure of a lamp and a\n"
    "calibration mask containing a regular grid of pinholes. The relation\n"
    "between mask and CCD positions is then determined.\n\n"
    "Input files:\n\n"
    "  DO category:              Type:       Explanation:         Required:\n"
    "  MASK_TO_CCD               Raw         Pinhole mask exposure   Y\n"
    "  MASTER_BIAS               Calib       Master bias             Y\n"
    "  MASTER_DARK               Calib       Master dark             .\n"
    "  IMG_MASTER_SKY_FLAT       Calib       Master sky flat         .\n"
    "  CCD_TABLE                 Calib       Bad pixel table         .\n\n"
    "Output files:\n\n"
    "  DO category:              Data type:  Explanation:\n"
    "  IMG_MASK_CCD_CONFIG       FITS file   Header with new transformation\n"
    "  (none)                    PAF         New transformation for IWS\n\n"
    "In the online processing, the PAF file is created in the working\n"
    "directory and also copied to the product directory. This file is\n"
    "formatted as the IWS configuration file IMG_mask2ccd_Q.cmf (where Q\n"
    "indicates the VIMOS quadrant number). A CCD table needs to be specified\n"
    "in input only if a bad pixel cleaning is requested.\n\n"
    "For more details, please refer to the VIMOS Pipeline User's Guide.",

    "ESO VIMOS Pipeline Team and VIMOS Consortium",

    PACKAGE_BUGREPORT,

    "This file is part of the VIMOS Instrument Pipeline\n"
    "Copyright (C) 2002-2009 European Southern Observatory\n\n"
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

                    vmmasktoccd_create,
                    vmmasktoccd_exec,
                    vmmasktoccd_destroy);

    cpl_pluginlist_append(list, plugin);
    
    return 0;

}
