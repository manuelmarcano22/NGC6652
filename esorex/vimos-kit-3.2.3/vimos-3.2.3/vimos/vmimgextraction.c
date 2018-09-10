/* $Id: vmimgextraction.c,v 1.6 2013-08-07 15:39:39 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2004 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-08-07 15:39:39 $
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#define _DEFAULT_SOURCE // For tempnam() (Obsolete in POSIX)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include <pilmemory.h>
#include <piltask.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilstrutils.h>
#include <pilutils.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmgalaxytable.h"
#include "vmastrometrictable.h"
#include "vmstarmatchtable.h"
#include "vmstartable.h"
#include "vmsextractor.h"
#include "vmimgutils.h"
#include "vmimgextraction.h"
#include "cpl.h"


#define SEXTRACTOR_ARGC  20
#define MIN_DEVIATION    1.e-6

/*
 * Mandatory SExtractor output columns. These columns define the VIMOS
 * galaxy table datatype.
 */

static SextParameter galaxy_table_columns[] = {
  {"NUMBER", SEXT_COLUMN_INT}, 
  {"MAG_ISOCOR", SEXT_COLUMN_FLOAT},
  {"MAGERR_ISOCOR", SEXT_COLUMN_FLOAT},
  {"MAG_APER", SEXT_COLUMN_FLOAT},
  {"MAGERR_APER", SEXT_COLUMN_FLOAT},
  {"MAG_AUTO", SEXT_COLUMN_FLOAT},
  {"MAGERR_AUTO", SEXT_COLUMN_FLOAT},
  {"MAG_BEST", SEXT_COLUMN_FLOAT},
  {"MAGERR_BEST", SEXT_COLUMN_FLOAT},
  {"X_IMAGE", SEXT_COLUMN_FLOAT},
  {"Y_IMAGE", SEXT_COLUMN_FLOAT},
  {"X_WORLD", SEXT_COLUMN_FLOAT},
  {"Y_WORLD", SEXT_COLUMN_FLOAT},
  {"ISOAREA_WORLD", SEXT_COLUMN_FLOAT},
  {"A_IMAGE", SEXT_COLUMN_FLOAT},
  {"B_IMAGE", SEXT_COLUMN_FLOAT},
  {"A_WORLD", SEXT_COLUMN_FLOAT},
  {"B_WORLD", SEXT_COLUMN_FLOAT},
  {"FWHM_IMAGE", SEXT_COLUMN_FLOAT},
  {"FWHM_WORLD", SEXT_COLUMN_FLOAT},
  {"THETA_IMAGE", SEXT_COLUMN_FLOAT},
  {"ERRTHETA_IMAGE", SEXT_COLUMN_FLOAT},
  {"ELLIPTICITY", SEXT_COLUMN_FLOAT},
  {"CLASS_STAR", SEXT_COLUMN_FLOAT},
  {"FLAGS", SEXT_COLUMN_INT},
  {0, SEXT_COLUMN_UNDEF}
};

/*  The following Sex columns are not used in DRS, thus they are NOT  */
/*  part of sex output BG */

/*    {"ISO0", SEXT_COLUMN_INT}, */
/*    {"ISO1", SEXT_COLUMN_INT}, */
/*    {"ISO2", SEXT_COLUMN_INT}, */
/*    {"ISO3", SEXT_COLUMN_INT}, */
/*    {"ISO4", SEXT_COLUMN_INT}, */
/*    {"ISO5", SEXT_COLUMN_INT}, */
/*    {"ISO6", SEXT_COLUMN_INT}, */
/*    {"ISO7", SEXT_COLUMN_INT}, */
/*    {"KRON_RADIUS", SEXT_COLUMN_FLOAT}, */
/*    {"BACKGROUND", SEXT_COLUMN_FLOAT}, */
/*    {"THRESHOLD", SEXT_COLUMN_FLOAT}, */
/*    {"MU_THRESHOLD", SEXT_COLUMN_FLOAT}, */
/*    {"FLUX_MAX", SEXT_COLUMN_FLOAT}, */
/*    {"MU_MAX", SEXT_COLUMN_FLOAT}, */
/*    {"ISOAREA_IMAGE", SEXT_COLUMN_FLOAT}, */
/*    {"A_IMAGE", SEXT_COLUMN_FLOAT}, */
/*    {"B_IMAGE", SEXT_COLUMN_FLOAT}, */
/*    {"THETA_IMAGE", SEXT_COLUMN_FLOAT}, */
/*    {"THETA_WORLD", SEXT_COLUMN_FLOAT}, */
/*    {"ELONGATION", SEXT_COLUMN_FLOAT}, */

/**
 * @name VmImObjectDetection
 *
 * @doc
 *   The modules provides functions for detecting objects in an image,
 *   operations on source lists including correlation with catalog data.
 */

/**@{*/

/*
 * @memo
 *   Append world coordinate system information to a keyword list.
 *
 * @return The function returns #EXIT_SUCCESS# if no error occurred,
 *   otherwise #EXIT_FAILURE# is returned.
 *
 * @param list  Target keyword list.
 * @param wcs   WCS information.
 *
 * @doc
 *   The function appends the world coordinate system information found
 *   in \textbf{wcs} to the list of keywords specified by \textbf{list}.
 *   Using the WCS information from \textbf{wcs} the function appends the
 *   keywords
 *   \begin{itemize}
 *     \item #CRVALi#,
 *     \item #CRPIXi#, and
 *     \item #CTYPEi#
 *   \end{itemize}
 *   for i = 1, 2, the
 *   CD matrix, i.e. the keywords #CDi_j# for i, j = 1, 2, and the keywords
 *   specifying the equinox and the reference frame #EQUINOX# and #RADECSYS#.
 *
 * @author R. Palsa
 */

inline static int
writeWcsInfo(VimosDescriptor **list, struct WorldCoor *wcs)
{

  register int i, j;
  VimosDescriptor *dsc = NULL;


  assert(list != 0 && wcs != 0);

  
  /*
   * Save the table type descriptor and remove it temporarily from the
   * list. It will be restored later.
   */

  if (vimosDscCopy(&dsc, *list, pilTrnGetKeyword("Table"), NULL) ==
      EXIT_FAILURE) {
      return EXIT_FAILURE;
  }
  else {
      vimosDscErase(list, pilTrnGetKeyword("Table"));
  }

  if (writeDoubleDescriptor(list, pilTrnGetKeyword("Crval", 1), wcs->xref,
                            "HH MM SS.SS, RA at ref pixel (deg)") == VM_FALSE) 
    return EXIT_FAILURE;

  if (writeDoubleDescriptor(list, pilTrnGetKeyword("Crval", 2), wcs->yref,
                            "DD MM SS, DEC at ref pixel (deg)") == VM_FALSE) 
    return EXIT_FAILURE;

  if (writeDoubleDescriptor(list, pilTrnGetKeyword("Crpix", 1), wcs->xrefpix,
                            "Ref pixel in X") == VM_FALSE) 
    return EXIT_FAILURE;

  if (writeDoubleDescriptor(list, pilTrnGetKeyword("Crpix", 2), wcs->yrefpix,
                            "Ref pixel in Y") == VM_FALSE) 
    return EXIT_FAILURE;

  if (writeStringDescriptor(list, pilTrnGetKeyword("Ctype", 1), "RA---TAN",
                            "Ref pixel in X") == VM_FALSE) 
    return EXIT_FAILURE;

  if (writeStringDescriptor(list, pilTrnGetKeyword("Ctype", 2), "DEC--TAN",
                            "Ref pixel in Y") == VM_FALSE) 
    return EXIT_FAILURE;


  /*
   * The 4 matrix elements of the CD matrix. Stored as (1, 1), (1, 2),
   * (2, 1) and (2, 2) in wcs->cd.
   */

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      if (writeDoubleDescriptor(list, pilTrnGetKeyword("CD", i + 1, j + 1),
                                wcs->cd[2 * i + j], "Translation matrix "
                                "element") == VM_FALSE)
        return EXIT_FAILURE;
    }
  }

  if (writeDoubleDescriptor(list, pilTrnGetKeyword("Equinox"), wcs->equinox,
                            "Standard FK5") == VM_FALSE)
    return EXIT_FAILURE;

  if (writeStringDescriptor(list, pilTrnGetKeyword("Radecsys"), wcs->radecsys,
                            "FK5") == VM_FALSE)
    return EXIT_FAILURE;


  /*
   * Restore the table type descriptor
   */

  if (vimosDscCopy(list, dsc, pilTrnGetKeyword("Table"), NULL) == 
      EXIT_FAILURE) {
      return EXIT_FAILURE;
  }
  else {
      deleteDescriptor(dsc);
  }

  return EXIT_SUCCESS;

}


/*
 * @memo
 *   Sort a table.
 *
 * @return The function returns a pointer to the sorted table if no error
 *   occurred, otherwise a #NULL# pointer is returned.
 *
 * @param table  Table object.
 * @param name   Name of the reference column.
 *
 * @doc
 *   The function looks for the column named \textbf{name} in the input table
 *   \textbf{table} and uses this column as a reference for sorting the
 *   table. The reference column must be of type \textbf{double} or
 *   \textbf{float} and its values are converted to type float for generating
 *   the index table, i.e you will loose precision (this will be fixed in
 *   future).
 *
 * @author R. Palsa
 */

inline static VimosTable *sortTable(VimosTable *table, const char *name)
{

  const char *fctid = "sortTable";


  int sz;

  struct bucket_t {
    size_t size;
    unsigned char *data;
  } bucket;

  VimosColumn *column;


  assert(table != 0);
  assert(name != 0);


  /*
   * Check for the reference column.
   */

  if (!(column = findColInTab(table, name))) {
    cpl_msg_debug(fctid, "Table column '%s' is missing!", name);
    return 0;
  }

  /*
   * Nothing to do if the table contains less than 2 columns
   */

  if ((sz = (colGetSize(column))) > 1) {
    register int i;
    int *indexTable;
    
    float *buffer;

    void *values;

    VimosColumn *c;


    /*
     * Check if all columns have the same size as the reference column.
     */

    c = table->cols;
    while (c) {
      if (colGetSize(c) != sz) {
        cpl_msg_debug(fctid, "Column sizes do not match!");
        return 0;
      }
      c = c->next;
    }


    /*
     * Generate the index table
     */

    /* FIXME: Is the extra element really needed?
     */

    if ((buffer = cpl_calloc(sz + 1, sizeof(float))) == 0)
      return 0;

    switch (column->colType) {
      case VM_FLOAT:
        values = colGetFloatData(column);

        for (i = 0; i < sz; i++)
          buffer[i] = (float) *((float *)values + i);

        break;

      case VM_DOUBLE:
        cpl_msg_debug(fctid, "Converting double to float. Possible loss of "
                    "precision!");

        values = colGetDoubleData(column);

        for (i = 0; i < sz; i++)
          buffer[i] = (float) *((double *)values + i);

        break;

      default:
        cpl_msg_debug(fctid, "Type of column '%s' is not supported!",
                    name);
        cpl_free(buffer);
        return 0;
    }


    /*
     * Allocate working buffers
     */

    if ((indexTable = cpl_calloc(sz, sizeof(int))) == 0)
      return 0;

    Indexx(sz, buffer, indexTable);
    cpl_free(buffer);


    /*
     * Rearrange the table. There is some memory needed for temporary
     * storage of the column data. This storage place is just a block of
     * bytes which is accessed through the appropriate handle. Its initial
     * size is calculated for a double column assuming that doubles are
     * the largest type, but for being on the safe side the actual size 
     * is checked and the byte buffer is enlarged if necessary.
     */

    bucket.size = sz * sizeof(double);
    bucket.data = (unsigned char *)cpl_calloc(bucket.size, sizeof(unsigned char));
    if (!bucket.data) {
      cpl_free(indexTable);
      return 0;
    }

    column = table->cols;

    while(column) {
      size_t bytes;

      switch (column->colType) {
      case VM_INT: {
        int *handle, *columnData;

        if ((bytes = sz * sizeof(int)) > bucket.size) {
          if (!cpl_realloc(bucket.data, bytes)) {
            cpl_free(bucket.data);
            cpl_free(indexTable);
            return 0;
          }
          else
            bucket.size = bytes;
        }

        columnData = colGetIntData(column);
        
        memcpy(bucket.data, columnData, bytes);
        handle = (int *)bucket.data;

        for (i = 0; i < sz; i++){

          columnData[i] = *(handle + indexTable[i]);
        }
      }
      break;

      case VM_FLOAT: {
        float *handle, *columnData;

        if ((bytes = sz * sizeof(float)) > bucket.size) {
          if (!cpl_realloc(bucket.data, bytes)) {
            cpl_free(bucket.data);
            cpl_free(indexTable);
            return 0;
          }
          else
            bucket.size = bytes;
        }

        columnData = colGetFloatData(column);
        memcpy(bucket.data, columnData, bytes);
        handle = (float *)bucket.data;

        for (i = 0; i < sz; i++)
          columnData[i] = *(handle + indexTable[i]);
          
      }
      break;

      case VM_DOUBLE: {
        double *handle, *columnData;

        if ((bytes = sz * sizeof(double)) > bucket.size) {
          if (!cpl_realloc(bucket.data, bytes)) {
            cpl_free(bucket.data);
            cpl_free(indexTable);
            return 0;
          }
          else
            bucket.size = bytes;
        }

        columnData = colGetDoubleData(column);
        
        memcpy(bucket.data, columnData, bytes);
        handle = (double *)bucket.data;

        for (i = 0; i < sz; i++)
          columnData[i] = *(handle + indexTable[i]);
          
      }
      break;

      case VM_STRING: {
        char **handle, **columnData;

        if ((bytes = sz * sizeof(char *)) > bucket.size) {
          if (!cpl_realloc(bucket.data, bytes)) {
            cpl_free(bucket.data);
            cpl_free(indexTable);
            return 0;
          }
          else
            bucket.size = bytes;
        }

        columnData = colGetStringData(column);
        
        memcpy(bucket.data, columnData, bytes);
        handle = (char **)bucket.data;

        for (i = 0; i < sz; i++) 
          columnData[i] = *(handle + indexTable[i]);
      }
      break;

      default:
        cpl_msg_debug(fctid, "Invalid column type!");
        cpl_free(bucket.data);
        cpl_free(indexTable);
        return 0;
        break;
      }

      column = column->next;

    }

    cpl_free(bucket.data);
    cpl_free(indexTable);

  }

  return table;

}


/*
 * @memo
 *   Prepare an astrometric catalog for processing.
 *
 * @return The function returns a pointer to the updated catalog if no error
 *    occurred, otherwise a #NULL# pointer is returned.
 *
 * @param catalog  Astrometric catalog.
 * @param filter   Filter name.
 *
 * @doc
 *   The function prepares the astrometric catalog \textbf{catalog} for
 *   the object identification task, i.e. the column containing the
 *   magnitudes for the band corresponding to the filter \textbf{filter},
 *   is copied to a new column with the name 'MAG'. Finally the catalog
 *   is sorted with respect to this column.
 *
 * @author R. Palsa
 */

inline static VimosTable *setupAstrometricTable(VimosTable *catalog,
                                                const char *filter)
{

  const char *fctid = "setupAstrometricTable";


  char name[6] = "MAG_";

  VimosColumn *column;


  assert(catalog != 0 && filter != 0);

  
  /*
   * Verify the type of the input table. Just rely on the table tag here.
   */

  if (strncmp(catalog->name, VM_AST, strlen(VM_AST))) {
    cpl_msg_debug(fctid, "Invalid type of input catalog!");
    return 0;
  }


  /*
   * Get the wavelength band designation from the filter name.
   * This seems to be pretty unreliable.
   */

  switch (*filter) {
  case 'U':
    strcat(name, "U");
    break;

  case 'B':
    strcat(name, "B");
    break;

  case 'V':
    strcat(name, "V");
    break;

  case 'R':
    strcat(name, "R");
    break;

  case 'I':
    strcat(name, "I");
    break;

  case 'z':
  case 'Z':
    strcat(name, "z");
    break;

  default:
    cpl_msg_debug(fctid, "Invalid filter name!");
    return 0;
    break;
  }


  /*
   * Copy the appropriate column and change its name and append it to
   * the table.
   */

  if (!(column = tblCopyColumn(catalog, name)))
    return 0;
  else {
    colSetName(column, "MAG");
    tblAppendColumn(catalog, column);
  }


  /*
   * Sort the catalog with respect to magnitude.
   */

  if (!sortTable(catalog, "MAG")) {
    cpl_msg_debug(fctid, "Cannot sort astrometric catalog by magnitude!");
    return 0;
  }

  return catalog;

}


/**
 * @memo
 *   Convert a SExtractor output table into a galaxy table.
 *
 * @return The function returns the updated table object if no error
 *   occurred, otherwise the return value is NULL;
 *
 * @param table  SExtractor output catalog table object.
 * @param image  Source detection image used for the catalog creation.
 *
 * @doc
 *   The function adds the galaxy table tag to the input source list
 *   generated by SExtractor and also all keywords required to make the
 *   input table a valid galaxy table are added to the table header.
 *   Additionally required unit conversions are applied.
 *
 *   The function also applies the sky to CCD transformation matrix to the
 *   coordinates of the detected objects and corrects for temperature
 *   effects.
 *
 * @author R. Palsa
 */

VimosTable *VmImBuildGalaxyTable(VimosTable *table, VimosImage *image)
{

  const char fctid[] = "VmImBuildGalaxyTable";


  const char *labels[] = {    /* Columns which need conversions */
    "A_WORLD",
    "B_WORLD",
    "FWHM_WORLD", 
    "ISOAREA_WORLD",
    "MAG_ISOCOR",
    "MAG_APER",
    "MAG_AUTO",
    "MAG_BEST",
    "X_IMAGE",
    "Y_IMAGE"
  };
  char comment[COMMENT_LENGTH];

  const size_t nc = sizeof labels / sizeof(const char *);

  int i, sz;
  int quadrant;

  double exposureTime, magnitudeCorrection, airmass;

  VimosColumn *columns[nc];
    

  /*
   * Get exposure time from the input observation and compute
   * the magnitude correction to account for the observation
   * duration.
   */

  if (readDoubleDescriptor(image->descs, pilTrnGetKeyword("ExposureTime"),
                           &exposureTime, comment) != VM_TRUE) {
    cpl_msg_error(fctid, "Cannot get exposure time!");
    return 0;
  }
  else
    magnitudeCorrection = 2.5 * log10(exposureTime);


  /*
   * Check if the input table has the properties of a galaxy table,
   * i.e. if all required columns and keywords a present.
   */

  for (i = 0; (size_t)i < nc; i++)
    if (!(columns[i] = findColInTab(table, labels[i]))) {
      cpl_msg_error(fctid, "Table column '%s' is missing!", labels[i]);
      return 0;
    }


  /*
   * Check that the columns all have the same size. Again this is needed
   * since a generic table may be passed.
   */

  sz = colGetSize(columns[0]);

  for (i = 1; (size_t)i < nc; i++) {
    if (colGetSize(columns[i]) != sz) {
      cpl_msg_error(fctid, "Column sizes do not match!");
      return 0;
    }
  }

      
  /*
   * Apply the conversions
   */

  for (i = 0; i < sz; i++) {

    /* FIXME:
     *  It is not obvious what is converted into what for the first 4
     *  columns. Check and add some explanation.
     */

    *(colGetDoubleData(columns[0]) + i) *= 3600.;
    *(colGetDoubleData(columns[1]) + i) *= 3600.;
    *(colGetDoubleData(columns[2]) + i) *= 3600.;
    *(colGetDoubleData(columns[3]) + i) *= 3600. * 3600.;

    *(colGetDoubleData(columns[4]) + i) += magnitudeCorrection;
    *(colGetDoubleData(columns[5]) + i) += magnitudeCorrection;
    *(colGetDoubleData(columns[6]) + i) += magnitudeCorrection;
    *(colGetDoubleData(columns[7]) + i) += magnitudeCorrection;


    /*
     * Shift image pixels by 1. SExtractor pixel coordinates
     * start at (0,0), while the FITS standard requires to
     * start at (1,1). See Greisen & Calabretta, A&A, preprint
     */

    *(colGetDoubleData(columns[8]) + i) += 1.;
    *(colGetDoubleData(columns[9]) + i) += 1.;
  }


  /*
   * Add necessary keywords to the table.
   */

  if (readIntDescriptor(image->descs, pilTrnGetKeyword("Quadrant"),
                        &quadrant, comment) != VM_TRUE)
    return 0;
  
  if (vimosDscCopy(&(table->descs), image->descs,
                   ".*-OBS$", 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("Instrument"), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   "^ESO OBS (DID|ID|PROG ID)", 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("INS.DID"), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("InstrumentMode"), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("FilterId", quadrant), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("FilterName", quadrant), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("DET.DID"), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("Adu2Electron", 1), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("ReadNoise", 1), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("Electron2Adu", 1), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("WINi.BINX", 1), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("WINi.BINY", 1), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("SeqWindowSizeX", 1), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("SeqWindowSizeY", 1), 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   "^ESO DET READ (CLOCK|SPEED|MODE)", 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   "^ESO OCS (CON QUAD|DID)", 0)) {
      return 0;
  }

  if (vimosDscCopy(&(table->descs), image->descs,
                   pilTrnGetKeyword("MagZero"), 0)) {
      return 0;
  }

  if (VmComputeAirmass(image, &airmass) == EXIT_FAILURE)
    return 0;
  else {
    if (writeDoubleDescriptor(&(table->descs), pilTrnGetKeyword("AirMass"),
                              airmass, "Averaged Airmass") != VM_TRUE)
      return 0;
  }


  /*
   * Finally the input table is now dressed up as a galaxy table and
   * deserves its tag.
   */
  
  strcpy(table->name, VM_GAL);
  if (writeStringDescriptor(&(table->descs), pilTrnGetKeyword("Table"),
                            VM_GAL, "Type of table") != VM_TRUE)
                            return 0;
  
  return table;

}

/**
 * @memo
 *   Create a star table from a galaxy table.
 *
 * @return The function returns a pointer to the created star table if no
 *   error occured, otherwise a NULL pointer is returned.
 *
 * @param table      Galaxy table.
 * @param starIndex  Stellarity index threshold.
 * @param magLimit   Magnitude lower limit.
 *
 * @doc
 *   The function creates a star table from a source list \textbf{table}
 *   having the properties of a galaxy table, i.e. the columns
 *   \begin{itemize}
 *     \item #NUMBER#,
 *     \item #X_IMAGE#,
 *     \item #Y_IMAGE#,
 *     \item #X_WORLD# and
 *     \item #Y_WORLD#
 *   \end{itemize}
 *   must be present. All objects in the input source list having a
 *   stellarity index greater than \textbf{starIndex}, which may have
 *   values ranging from 0 to 1, and which are brighter than
 *   \textbf{magLimit} are considered as stars.
 *
 *   A preselection of the detected objects is done on the basis
 *   of the extraction attributes (source list column #FLAGS#). Only
 *   'clean' source detections are considered, i.e. detections have
 *   no flags set (#FLAGS# = 0).
 *
 *   The output star table is created and the selected entries from
 *   \textbf{table} corresponding to stars are copied.
 *
 * @see VmImBuildGalaxyTable, VmImDetectObjects
 *   
 * @author P. Sartoretti, R. Palsa
 */

VimosTable *VmImBuildStarTable(VimosTable *table, float starIndex,
                               float magLimit)
{

  const char  fctid[] = "VmImBuildStarTable";


  const char *labels[] = {     /* Columns needed from the galaxy table */
    "NUMBER",
    "X_IMAGE",
    "Y_IMAGE",
    "X_WORLD",
    "Y_WORLD",
    "FLAGS",
    "CLASS_STAR",
    "MAG_BEST"
  };

  const size_t nc = sizeof labels / sizeof(const char *);
  size_t i, objectCount, starCount;

  int selection;
  int *flags;
  int *selectedEntry;
  int *star;

  double *class;
  double *ximage, *yimage;
  double *xworld, *yworld;
  double *magBest, *magnitude;

  VimosColumn *columns[nc];

  VimosTable  *starTable;


  /*
   * Validate input
   */

  assert(table != 0);

  if((starIndex < 0.) || (starIndex > 1.)) {
    cpl_msg_error(fctid, "Stellarity index is out of range!");
    return 0;
  }


  /*
   * Check if the input table has the properties of a galaxy table,
   * i.e. if all required columns and keywords a present.
   */

  for (i = 0; i < nc; i++)
    if (!(columns[i] = findColInTab(table, labels[i]))) {
      cpl_msg_error(fctid, "Table column '%s' is missing!", labels[i]);
      return 0;
    }

  cpl_msg_debug(fctid, "Stellarity index threshold: %.3f", starIndex);
  cpl_msg_debug(fctid, "Magnitude limit: %.3f", magLimit);

  
  /*
   * Find all bright stars in the input table and record their
   * sequence number
   */

  objectCount = colGetSize(columns[0]);
  selectedEntry = (int *)cpl_calloc(objectCount, sizeof(int));

  /*
   * Check the selection criteria and count the stars.
   *
   * Column 5: FLAGS
   * Column 6: CLASS_STAR
   * Column 7: MAG_BEST
   */

  flags = colGetIntData(columns[5]);
  class = colGetDoubleData(columns[6]);
  magBest = colGetDoubleData(columns[7]);

  for(i = 0, starCount = 0; i < objectCount; i++) {
    if (flags[i] == 0 && class[i] > starIndex && magBest[i] < magLimit) {
      selectedEntry[starCount] = i;
      starCount++;
    }
  }


  /*
   * Beware of empty table
   */

  if (starCount <= 0) {
    cpl_msg_warning(fctid, "No stars found for current settings!");
    starTable = newStarTableEmpty();
    cpl_free(selectedEntry);
    return 0;
  }
  else {
    cpl_msg_info(fctid, "%zd stars have been selected.", starCount);
    starTable = newStarTable(starCount);
  }


  /*
   * Setup the output table and its columns
   */

  if (!starTable) {
    cpl_msg_error(fctid, "Cannot create star table!");
    cpl_free(selectedEntry);
    return 0;
  }
  else {
      char *hint = cpl_strdup(pilTrnGetKeyword("Table"));

      vimosDscCopy(&(starTable->descs), table->descs,
                   ".*-OBS$", hint);
      vimosDscCopy(&(starTable->descs), table->descs,
                   pilTrnGetKeyword("Instrument"), hint);
      vimosDscCopy(&(starTable->descs), table->descs,
                   "^ESO (OBS|INS|DET|OCS)", hint);
      vimosDscCopy(&(starTable->descs), table->descs,
                   "^ESO PRO (MAG ZERO|AIRMASS)", hint);

      cpl_free(hint);
  }
          

  if (starCount > 0) {

    /*
     * Copy NUMBER, X_IMAGE, Y_IMAGE, X_WORLD, Y_WORLD and MAG of the
     * selected stars from the galaxy table to the output star table.
     */

    star = tblGetIntData(starTable, "NUMBER");
    ximage = tblGetDoubleData(starTable, "X_IMAGE");
    yimage = tblGetDoubleData(starTable, "Y_IMAGE");
    xworld = tblGetDoubleData(starTable, "X_WORLD");
    yworld = tblGetDoubleData(starTable, "Y_WORLD");
    magnitude = tblGetDoubleData(starTable, "MAG");

    for (i = 0; i < starCount; i++) {
      selection = selectedEntry[i];
  
      star[i] = *(colGetIntData(columns[0]) + selection);
  
      ximage[i] = *(colGetDoubleData(columns[1]) + selection);
      yimage[i] = *(colGetDoubleData(columns[2]) + selection);
      xworld[i] = *(colGetDoubleData(columns[3]) + selection);
      yworld[i] = *(colGetDoubleData(columns[4]) + selection);
      magnitude[i] = *(colGetDoubleData(columns[7]) + selection);
    }
     
  }

  cpl_free(selectedEntry);

  return starTable;
}


/**
 * @memo
 *   Build a star match table from a list of stars and a catalog.
 *
 * @return The function returns a pointer to the created star match table
 *   if no error occured, otherwise a NULL pointer is returned.
 *
 * @param image          Image with WCS information.
 * @param starTable      List of detected stars.
 * @param astTable       Astrometric catalog data.
 * @param minStars       Minimum number of stars to be found.
 * @param searchRadius   Search radius in arcsec.
 * @param iMagTolerance  Magnitude tolerance for the initial selection
 *                       of matching stars.
 * @param fMagTolerance  Magnitude tolerance for the final selection of
 *                       matching stars.
 * @param sigmaClip
 *
 * @doc
 *   The function identifies the detected stellar objects provided by
 *   the source list \textbf{starTable} with stars listed in a astrometric
 *   catalog \textbf{astTable}. The detected objects are identified using
 *   \textbf{VmImSearchMatches()}. For all identified objects, the 
 *   information from the source list and the astrometric catalog are
 *   merged into the result star match table. For a detailed information
 *   on the remaining parameters see the documentation of
 *   \textbf{VmImSearchMatches()}.
 *
 * @see VmImSearchMatches
 *
 * @author P. Montegriffo, P. Sartoretti, B. Garilli, R. Palsa
 */

VimosTable *VmImBuildStarMatchTable(VimosImage *image, VimosTable *starTable,
                                    VimosTable *astTable, int minStars,
                                    double searchRadius, double iMagTolerance,
                                    double fMagTolerance, float sigmaClip)
{

  const char fctid[] = "VmImBuildStarMatchTable";

  register int i;

  char comment[COMMENT_LENGTH];
  char filter[FLEN_VALUE + 1];
  char **star_name;

  int quadrant, matchCount;
  int *indexTable;
  int *star;

  double *ximage, *yimage;
  double *xworld, *yworld;
  double *ra, *dec;
  double *magnitude;

  struct WorldCoor *wcs;

  VimosTable *stmcTable;


  /*
   * Validate input
   */

  assert(image != 0);
  assert(starTable != 0 && astTable != 0);

/*
  if (starTable->cols == 0 || starTable->cols->len == 0) {
    cpl_msg_error(fctid, "Empty input source list!");
    return 0;
  }
*/

  if (astTable->cols == 0 || astTable->cols->len == 0) {
    cpl_msg_error(fctid, "Empty astrometric catalog!");
    return 0;
  }


  /*
   * Retrieve the filter name used for the observation from the
   * input image.
   */

  if (readIntDescriptor(image->descs, pilTrnGetKeyword("Quadrant"), &quadrant,
                        comment) == VM_FALSE)
    return 0;


  if (readStringDescriptor(image->descs,
                           pilTrnGetKeyword("FilterName", quadrant), filter,
                           comment) == VM_FALSE)
      return 0;


  /*
   * Get world coordinate system from image header
   */

  if (!(wcs = rdimage(image->descs))) {
    cpl_msg_error(fctid, "World coordinate system not found in input image");
    return 0;
  }


  /*
   * Setup the astrometric catalog for the used filter, sort it by
   * magnitude and check that enough astrometric stars are available.
   */

  if (!setupAstrometricTable(astTable, filter)) {
    cpl_msg_error(fctid, "Astrometric table setup failed!");
    return 0;
  }

  if (astTable->cols->len < minStars) {
    cpl_msg_error(fctid, "Too few entries in astrometric catalog for filter %s!",
                filter);
    return 0;
  }


  /*
   * Sort the input source list by magnitude
   */

  if ((starTable->cols) && starTable->cols->len > 0) {
    if (!sortTable(starTable, "MAG")) {
      cpl_msg_error(fctid, "Cannot sort source list by magnitude!");
      return 0;
    }


    /*
     * Project astrometric stars into the image pixel coordinate system
     */

    cpl_msg_debug(fctid, "Projecting astrometric stars into image plane!");

    wcstopix(astTable->cols->len, astTable, wcs);


    /*
     * Convert serch radius from arcsec to pixel.
     */

    searchRadius /= fabs(3600. * wcs->cdelt[0]);


    /*
     * Identify sources from the source list with catalog entries.
     */

    cpl_msg_debug(fctid, "Searching for matching stars.");

    indexTable = VmSearchMatches(starTable, astTable, searchRadius, 
                                 iMagTolerance, fMagTolerance, sigmaClip,
                                 minStars, &matchCount);

    if (!indexTable) {
      char *hint = cpl_strdup(pilTrnGetKeyword("Table"));

      cpl_msg_warning(fctid, "Search for matching stars failed!");

      stmcTable = newStarMatchTableEmpty();

      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   ".*-OBS$", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   pilTrnGetKeyword("Instrument"), hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO (OBS|INS|DET|OCS)", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO PRO (MAG ZERO|AIRMASS)", hint);

      cpl_free(hint);

      return stmcTable;
    }
    else 
      cpl_msg_info(fctid, "%d matching stars found.", matchCount);


    /*
     * Project stars onto sky. World coordinates are needed for the 
     * StarMatch table.
     */

    cpl_msg_debug(fctid, "Projecting stars onto sky!");

    pixtowcs(starTable->cols->len, starTable, wcs);


    /*
     * Build the StarMatch table.
     */

    if (!(stmcTable = newStarMatchTable(matchCount))) {
      cpl_free(indexTable);
      return 0;
    }
  }
  else {
      char *hint = cpl_strdup(pilTrnGetKeyword("Table"));

      stmcTable = newStarMatchTableEmpty();

      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   ".*-OBS$", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   pilTrnGetKeyword("Instrument"), hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO (OBS|INS|DET|OCS)", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO PRO (MAG ZERO|AIRMASS)", hint);

      cpl_free(hint);

      return stmcTable;
  }


  /*
   * Setup the star match table keyword list, i.e. copy the keywords
   * from the star table and insert the world coordinate system keywords.
   */

  if (writeWcsInfo(&stmcTable->descs, wcs) == EXIT_FAILURE) {
    cpl_msg_debug(fctid, "Appending WCS keywords failed!");
    
    cpl_free(indexTable);
    deleteTable(stmcTable);

    return 0;
  }
  else {
      char *hint = cpl_strdup(pilTrnGetKeyword("Table"));

      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   ".*-OBS$", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   pilTrnGetKeyword("Instrument"), hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO (OBS|INS|DET|OCS)", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO PRO (MAG ZERO|AIRMASS)", hint);

      cpl_free(hint);
  }

  /*
   * Copy the data for the matching entries from the star and astrometric
   * table.
   */

  star = tblGetIntData(starTable, "NUMBER");
  star_name = tblGetStringData(astTable, "ID");
  ximage = tblGetDoubleData(starTable, "X_IMAGE");
  yimage = tblGetDoubleData(starTable, "Y_IMAGE");
  xworld = tblGetDoubleData(starTable, "X_WORLD");
  yworld = tblGetDoubleData(starTable, "Y_WORLD");
  magnitude = tblGetDoubleData(starTable, "MAG");
  ra = tblGetDoubleData(astTable, "RA");
  dec = tblGetDoubleData(astTable, "DEC");
  
  for (i = 0; i < matchCount; i++) {
    int j, k;

    j = indexTable[2 * i];
    k = indexTable[2 * i + 1];

    tblSetIntValue(stmcTable, "NUMBER", i, star[j]);
    tblSetStringValue(stmcTable, "ID", i, star_name[k]);
    tblSetDoubleValue(stmcTable, "X_IMAGE", i, ximage[j]);
    tblSetDoubleValue(stmcTable, "Y_IMAGE", i, yimage[j]);
    tblSetDoubleValue(stmcTable, "X_WORLD", i, xworld[j]);
    tblSetDoubleValue(stmcTable, "Y_WORLD", i, yworld[j]);
    tblSetDoubleValue(stmcTable, "MAG", i, magnitude[j]);
    tblSetDoubleValue(stmcTable, "RA", i, ra[k]);
    tblSetDoubleValue(stmcTable, "DEC", i, dec[k]);
  }


  /*
   * Extra magnutude columns from the astrometric table. Not all of
   * them might be present!
   */

  if (!findColInTab(astTable, "MAG_U"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_U"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_U");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_U", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_B"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_B"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_B");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_B", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_V"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_V"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_V");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_V", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_R"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_R"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_R");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_R", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_I"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_I"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_I");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_I", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_z"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_z"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_z");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_z", i, m[indexTable[2 * i + 1]]);
  }


  /*
   * Cleanup
   */

  cpl_free(indexTable);

  return stmcTable;

}

/*PDB - New version of VmImBuildStarMatchTable that has a fix for vmskyccd*/
VimosTable *VmImBuildStarMatchTable_skyccd(VimosImage *image, VimosTable *starTable,
                                    VimosTable *astTable, int minStars,
                                    double searchRadius, double iMagTolerance,
                                    double fMagTolerance, float sigmaClip)
{

  const char fctid[] = "VmImBuildStarMatchTable";

  register int i;

  char comment[COMMENT_LENGTH];
  char filter[FLEN_VALUE + 1];
  char **star_name;

  int quadrant, matchCount;
  int *indexTable;
  int *star;

  double *ximage, *yimage;
  double *xworld, *yworld;
  double *ra, *dec;
  double *magnitude;

  struct WorldCoor *wcs;

  VimosTable *stmcTable;


  /*
   * Validate input
   */

  assert(image != 0);
  assert(starTable != 0 && astTable != 0);

/*
  if (starTable->cols == 0 || starTable->cols->len == 0) {
    cpl_msg_error(fctid, "Empty input source list!");
    return 0;
  }
*/

  if (astTable->cols == 0 || astTable->cols->len == 0) {
    cpl_msg_error(fctid, "Empty astrometric catalog!");
    return 0;
  }


  /*
   * Retrieve the filter name used for the observation from the
   * input image.
   */

  if (readIntDescriptor(image->descs, pilTrnGetKeyword("Quadrant"), &quadrant,
                        comment) == VM_FALSE)
    return 0;


  if (readStringDescriptor(image->descs,
                           pilTrnGetKeyword("FilterName", quadrant), filter,
                           comment) == VM_FALSE)
      return 0;


  /*
   * Get world coordinate system from image header
   */

  if (!(wcs = rdimage(image->descs))) {
    cpl_msg_error(fctid, "World coordinate system not found in input image");
    return 0;
  }


  /*
   * Setup the astrometric catalog for the used filter, sort it by
   * magnitude and check that enough astrometric stars are available.
   */

  if (!setupAstrometricTable(astTable, filter)) {
    cpl_msg_error(fctid, "Astrometric table setup failed!");
    return 0;
  }

  if (astTable->cols->len < minStars) {
    cpl_msg_error(fctid, "Too few entries in astrometric catalog for filter %s!",
                filter);
    return 0;
  }


  /*
   * Sort the input source list by magnitude
   */

  if ((starTable->cols) && starTable->cols->len > 0) {
    if (!sortTable(starTable, "MAG")) {
      cpl_msg_error(fctid, "Cannot sort source list by magnitude!");
      return 0;
    }


    /*
     * Project astrometric stars into the image pixel coordinate system
     */

    cpl_msg_debug(fctid, "Projecting astrometric stars into image plane!");
    wcstopix(astTable->cols->len, astTable, wcs);


    /*
     * Convert serch radius from arcsec to pixel.
     */

    searchRadius /= fabs(3600. * wcs->cdelt[0]);

    /*
     * Identify sources from the source list with catalog entries.
     */

    cpl_msg_debug(fctid, "Searching for matching stars.");

    indexTable = VmSearchMatches(starTable, astTable, searchRadius, 
                                 iMagTolerance, fMagTolerance, sigmaClip,
                                 minStars, &matchCount);

    if (!indexTable) {
      char *hint = cpl_strdup(pilTrnGetKeyword("Table"));

      cpl_msg_warning(fctid, "Search for matching stars failed!");

      stmcTable = newStarMatchTableEmpty();

      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   ".*-OBS$", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   pilTrnGetKeyword("Instrument"), hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO (OBS|INS|DET|OCS)", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO PRO (MAG ZERO|AIRMASS)", hint);

      cpl_free(hint);

      return stmcTable;
    }
    else 
      cpl_msg_info(fctid, "%d matching stars found.", matchCount);


    /*
     * Project stars onto sky. World coordinates are needed for the 
     * StarMatch table.
     */

    cpl_msg_debug(fctid, "Projecting stars onto sky!");

    pixtowcs(starTable->cols->len, starTable, wcs);


    /*
     * Build the StarMatch table.
     */

    if (!(stmcTable = newStarMatchTable(matchCount))) {
      cpl_free(indexTable);
      return 0;
    }
  }
  else {
      char *hint = cpl_strdup(pilTrnGetKeyword("Table"));

      stmcTable = newStarMatchTableEmpty();

      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   ".*-OBS$", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   pilTrnGetKeyword("Instrument"), hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO (OBS|INS|DET|OCS)", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO PRO (MAG ZERO|AIRMASS)", hint);

      cpl_free(hint);

      return stmcTable;
  }


  /*
   * Setup the star match table keyword list, i.e. copy the keywords
   * from the star table and insert the world coordinate system keywords.
   */

  if (writeWcsInfo(&stmcTable->descs, wcs) == EXIT_FAILURE) {
    cpl_msg_debug(fctid, "Appending WCS keywords failed!");
    
    cpl_free(indexTable);
    deleteTable(stmcTable);

    return 0;
  }
  else {
      char *hint = cpl_strdup(pilTrnGetKeyword("Table"));

      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   ".*-OBS$", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   pilTrnGetKeyword("Instrument"), hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO (OBS|INS|DET|OCS)", hint);
      vimosDscCopy(&stmcTable->descs, starTable->descs,
                   "^ESO PRO (MAG ZERO|AIRMASS)", hint);

      cpl_free(hint);
  }

  /*
   * Copy the data for the matching entries from the star and astrometric
   * table.
   */

  star = tblGetIntData(starTable, "NUMBER");
  star_name = tblGetStringData(astTable, "ID");
  /*PDB table is astTable in 2 lines below, this is specific to the
      VmImBuildStarMatchTable_skyccd version of the function and
      is necessary for vmskyccd to work properly*/
  ximage = tblGetDoubleData(astTable, "X_IMAGE");
  yimage = tblGetDoubleData(astTable, "Y_IMAGE");
  xworld = tblGetDoubleData(starTable, "X_WORLD");
  yworld = tblGetDoubleData(starTable, "Y_WORLD");
  magnitude = tblGetDoubleData(starTable, "MAG");
  ra = tblGetDoubleData(astTable, "RA");
  dec = tblGetDoubleData(astTable, "DEC");
  
  for (i = 0; i < matchCount; i++) {
    int j, k;

    j = indexTable[2 * i];
    k = indexTable[2 * i + 1];

    tblSetIntValue(stmcTable, "NUMBER", i, star[j]);
    tblSetStringValue(stmcTable, "ID", i, star_name[k]);
    /*PDB index is [k] in 2 lines below, this is specific to the
      VmImBuildStarMatchTable_skyccd version of the function and
      is necessary for vmskyccd to work properly*/
    tblSetDoubleValue(stmcTable, "X_IMAGE", i, ximage[k]);
    tblSetDoubleValue(stmcTable, "Y_IMAGE", i, yimage[k]);
    tblSetDoubleValue(stmcTable, "X_WORLD", i, xworld[j]);
    tblSetDoubleValue(stmcTable, "Y_WORLD", i, yworld[j]);
    tblSetDoubleValue(stmcTable, "MAG", i, magnitude[j]);
    tblSetDoubleValue(stmcTable, "RA", i, ra[k]);
    tblSetDoubleValue(stmcTable, "DEC", i, dec[k]);
  }


  /*
   * Extra magnutude columns from the astrometric table. Not all of
   * them might be present!
   */

  if (!findColInTab(astTable, "MAG_U"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_U"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_U");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_U", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_B"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_B"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_B");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_B", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_V"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_V"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_V");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_V", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_R"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_R"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_R");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_R", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_I"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_I"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_I");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_I", i, m[indexTable[2 * i + 1]]);
  }

  if (!findColInTab(astTable, "MAG_z"))
    deleteColumn(tblRemoveColumn(stmcTable, "MAG_z"));
  else {
    double *m = tblGetDoubleData(astTable, "MAG_z");

    for (i = 0; i < matchCount; i++)
      tblSetDoubleValue(stmcTable, "MAG_z", i, m[indexTable[2 * i + 1]]);
  }


  /*
   * Cleanup
   */

  cpl_free(indexTable);

  return stmcTable;

}


/**
 * @memo 
 *   Detect objects in an imaging observation.
 *
 * @return The function returns a pointer to the generated source list if
 *   no error occurred, othewise the return value is NULL.
 *
 * @param image      Image which should be searched for objects.
 * @param badPixels  Bad pixel map. (This is not supported. Set it to 0.)
 *
 * @doc 
 *
 * @author R. Palsa
 */

VimosTable *VmImDetectObjects(VimosImage *image, VimosImage *badpixels)
{

  const char fctid[] = "VmImDetectObjects";

  char *configName, *parameterName, *imageName;
  char *networkName, *filterName, *catalogName;
  char *assocName = 0;
  char *checkName = 0;
  char *flagName = 0;
  char *weightName = 0;
  char *argv[SEXTRACTOR_ARGC + 1];
  char cwd[PATHNAME_MAX + 1];

  int argc;

  time_t timeout = sextGetExecutionTimeLimit();

  FILE *configFile, *parameterFile;

  VimosTable *galTable;


  /*
   * Prepare the SExtractor setup files. The input setup is taken
   * from the recipe database. The output is defined here to be
   * exactly what is needed for building a galaxy table. The setup
   * for the galaxy/star classifier (neuronal network setup) as well
   * the convolution mask setup are just read from the files given
   * in the database.
   */

  /*
   * Get the current working directory.
   */

  if (!getcwd(cwd, PATHNAME_MAX)) {
    cpl_msg_error(fctid, "Cannot determine current working directory!");
    return 0;
  }


  /*
   * Create and write the SExtractor setup file. The function tempnam
   * is used to create a name for a temporary file. The name has to
   * be passed to SExtractor so it has to be kept.
   */

  if (!(configName = tempnam(cwd, "sext"))) {
    cpl_msg_error(fctid, "Cannot create unique name for temporary file!");
    return 0;
  }
  else {
    if (!(configFile = fopen(configName, "w"))) {
      cpl_msg_error(fctid, "Cannot create temporary setup file!");
      
      cpl_free(configName);

      return 0;
    }
  }

  if (sextSaveConfiguration(configFile, image) == EXIT_FAILURE) {
    cpl_msg_error(fctid, "Cannot dump SExtractor setup!");

    fclose(configFile);
    remove(configName);

    cpl_free(configName);

    return 0;
  }

  fclose(configFile);


  /*
   * Create and write the SExtractor output configuration file. The 
   * function tempnam is used to create a name for a temporary file.
   * The name has to be passed to SExtractor so it has to be kept.
   */

  if (!(parameterName = tempnam(cwd, "sext"))) {
    cpl_msg_error(fctid, "Cannot create unique name for temporary file!");

    remove(configName);

    cpl_free(configName);

    return 0;
  }
  else {
    if (!(parameterFile = fopen(parameterName, "w"))) {
      cpl_msg_error(fctid, "Cannot create temporary setup file!");
      
      remove(configName);

      cpl_free(configName);
      cpl_free(parameterName);

      return 0;
    }
  }

  if (sextSaveParameters(parameterFile, galaxy_table_columns) == 
      EXIT_FAILURE) {
    cpl_msg_error(fctid, "Cannot write SExtractor parameter file!");

    fclose(parameterFile);

    remove(configName);
    remove(parameterName);

    cpl_free(configName);
    cpl_free(parameterName);

    return 0;
  }

  fclose(parameterFile);


  /*
   * Get the fully expanded path to the neuronal network and convolution
   * mask setup files. Their existance and permissions are checked.
   */

  if (!(networkName = cpl_strdup(sextGetStarNnwName()))) {
    cpl_msg_error(fctid, "Cannot retrieve SExtractor neuronal network "
                "setup file!");

    remove(configName);
    remove(parameterName);

    cpl_free(configName);
    cpl_free(parameterName);

    return 0;
  }
  else 
    if (access(networkName, F_OK | R_OK)) {
      cpl_msg_error(fctid, "Cannot access SExtractor neuronal network "
                  "setup file %s!", networkName);

      remove(configName);
      remove(parameterName);

      cpl_free(configName);
      cpl_free(parameterName);
      cpl_free(networkName);

      return 0;
    }


  if (!(filterName = cpl_strdup(sextGetFilterName()))) {
    cpl_msg_error(fctid, "Cannot retrieve SExtractor filter setup file!");

    remove(configName);
    remove(parameterName);

    cpl_free(configName);
    cpl_free(parameterName);
    cpl_free(networkName);

    return 0;
  }
  else
    if (access(filterName, F_OK | R_OK)) {
      cpl_msg_error(fctid, "Cannot access SExtractor neuronal network "
                  "setup file %s!", filterName);

      remove(configName);
      remove(parameterName);

      cpl_free(configName);
      cpl_free(parameterName);
      cpl_free(networkName);
      cpl_free(filterName);

      return 0;
    }

  
  /* FIXME:
   *  A better solution for generating the output catalog name and the
   *  input image name might be to use one of the file name header keywords
   *  in the image together with the appropriate suffix. Also, it
   *  should be forced that image and catalog file are generated in the
   *  current working directory. (RP)
   */

  /*
   * Generate the output catalog name
   */

  if (!(catalogName = tempnam(cwd, "sext"))) {
    cpl_msg_error(fctid, "Cannot create unique name for temporary file!");

    remove(configName);
    remove(parameterName);

    cpl_free(configName);
    cpl_free(parameterName);
    cpl_free(networkName);
    cpl_free(filterName);

    return 0;
  }

  /*
   * Save input images to local disk files which can be used by
   * SExtractor. Also the CDELT keywords are removed here, which
   * confuse SExtractor if the CD matrix is present too.
   */

  if (findDescriptor(image->descs, pilTrnGetKeyword("CD", 1, 1))) {
    if (findDescriptor(image->descs, pilTrnGetKeyword("Cdelt", 1)))
      removeDescriptor(&image->descs, pilTrnGetKeyword("Cdelt", 1));

    if (findDescriptor(image->descs, pilTrnGetKeyword("Cdelt", 2)))
      removeDescriptor(&image->descs, pilTrnGetKeyword("Cdelt", 2));
  }

  if (!(imageName = tempnam(cwd, "sext"))) {
    cpl_msg_error(fctid, "Cannot create unique name for temporary file!");

    remove(configName);
    remove(parameterName);

    cpl_free(configName);
    cpl_free(parameterName);
    cpl_free(networkName);
    cpl_free(filterName);
    cpl_free(catalogName);

    return 0;
  }
  else {
    if (!createFitsImage(imageName, image, "UNKNOWN")) {
      cpl_msg_error(fctid, "Cannot create temporary SExtractor input image "
                  "file!");

      remove(configName);
      remove(parameterName);

      cpl_free(configName);
      cpl_free(parameterName);
      cpl_free(networkName);
      cpl_free(filterName);
      cpl_free(imageName);
      cpl_free(catalogName);

      return 0;
    }
  }
  

  /*
   * Build the SExtractor command line and execute it.
   */

  argv[0] = cpl_strdup(sextGetSextractorPath());
  argv[1] = imageName;
  argv[2] = "-c";
  argv[3] = configName;
  argv[4] = "-PARAMETERS_NAME";
  argv[5] = parameterName;
  argv[6] = "-CATALOG_NAME";
  argv[7] = catalogName;

  argc = 7;

  if (filterName) {
    argv[++argc] = "-FILTER_NAME";
    argv[++argc] = filterName;
  }

  if (networkName) {
    argv[++argc] = "-STARNNW_NAME";
    argv[++argc] = networkName;
  }

  /*
   * Note that the following assignements refer to an static string
   * of the functions on the right side. Therefore the assignement
   * should be considered valid only within the enclosing 
   * if-statement.
   */

  if ((assocName = (char *)sextGetAssocName())) {
    argv[++argc] = "-ASSOC_NAME";
    argv[++argc] = assocName;
  }
    
  if ((checkName = (char *)sextGetCheckImageName())) {
    argv[++argc] = "-CHECKIMAGE_NAME";
    argv[++argc] = checkName;
  }
    
  if ((flagName = (char *)sextGetFlagImageName())) {
    argv[++argc] = "-FLAG_IMAGE";
    argv[++argc] = flagName;
  }

  if ((weightName = (char *)sextGetWeightImageName())) {
    argv[++argc] = "-WEIGHT_IMAGE";
    argv[++argc] = weightName;
  }

  argv[++argc] = 0;

  if (pilTaskExecWait(argc, (const char *const *)argv, timeout)) {
    cpl_msg_error(fctid, "Running SExtractor failed!");

    remove(configName);
    remove(parameterName);
    remove(imageName);

    cpl_free(configName);
    cpl_free(parameterName);
    cpl_free(networkName);
    cpl_free(filterName);
    cpl_free(imageName);
    cpl_free(catalogName);

    return 0;
  }

  remove(configName);
  remove(parameterName);
  remove(imageName);

  cpl_free(configName);
  cpl_free(parameterName);
  cpl_free(networkName);
  cpl_free(filterName);
  cpl_free(imageName);


  /*
   * Check the validity of SExtractor output source list and build
   * a galaxy table from it.
   */

  if (!(galTable = sextConvertCatalog(catalogName, galaxy_table_columns))) {
    cpl_msg_error(fctid, "SExtractor output catalog conversion failed!");

    remove(catalogName);

    cpl_free(catalogName);

    return 0;
  }
  else {
    remove(catalogName);
    cpl_free(catalogName);
  }

  if (!(galTable = VmImBuildGalaxyTable(galTable, image))) {
    cpl_msg_error(fctid, "Building Galaxy table from SExtractor output catalog "
                "failed!");

    deleteTable(galTable);
    return 0;
  }

  return galTable;
}


/*
 * The function below is basically a copy of searchmatch(). Modifications
 * have been applied to the prototype in the sense that the size of the
 * columns are taken from the input tables instead of being passed, and
 * the order of the arguments has changed. The function body has been 
 * cleaned up with respect to comments and input only!
 */


/**
 * @memo
 *   Search for matching objects.
 *
 * @return The function returns the index table for the matching objects if
 *   no error occurred, otherwise a #NULL# pointer is returned.
 *
 * @param sourceList     Table containing detected objects.
 * @param catalog        Reference catalog.
 * @param radius         Search radius in pixel.
 * @param magTolerance1  Magnitude tolerance for initial selection of
 *                       matching objects.
 * @param magTolerance2  Magnitude tolerance for final selection of
 *                       matching objects.
 * @param sigmaClip      Sigma clipping threshold.
 * @param minStars       Required minimum number of matching pairs.
 * @param matches        Number of matching pairs found.
 *
 * @doc
 *   The function identifies the sources listed in the source list
 *   \textbf{sourceList}. A detected source from \textbf{sourceList} is
 *   identified by first looking for catalog objects listed in 
 *   \textbf{catalog} within \textbf{radius} pixel around the position of
 *   the detected objects and, additionally, discriminating on the difference
 *   in brightness of the detected objects and the catalog sources, using
 *   \textbf{magTolerance1} as an upper limit. During this first selection
 *   at least \textbf{minStars} must be found, or the function returns an
 *   error. In a second step a selection among the surviving objects is
 *   performed applying a sigma clipping method to the positional offsets
 *   from the mean offset of the catalog sources from the detected positions.
 *   The sensitivity of the sigma clipping algorithm is controlled by the
 *   sigma clipping factor \textbf{sigmaClip}. Again the discrimination on
 *   the difference in brightness is applied, this time using
 *   \textbf{magTolerance2} as upper limit. If \textbf{magTolerance1}
 *   or \textbf{magTolerance2} are less than zero, the corresponding
 *   magnitude discrimination is ignored.
 *
 *   The final number of matching pairs is passed back to the caller through
 *   the argument \textbf{matches}.
 *
 * @author P. Montegriffo, P. Sartoretti, R. Palsa
 */

int *VmSearchMatches(VimosTable *sourceList, VimosTable *catalog,
                     double radius, double magTolerance1, double magTolerance2,
                     float sigmaClip, int minStars, int *matches)
{

  const char *fctid = "VmSearchMatches";


  int i, j, jbest;
  int sz, n_ele1, n_ele2;
  int nmatch0;
  int nrej;
  int nmatch = 0;
  int *imatch = 0;     /* returning index table for matching stars */
  int *imatch0 = 0;    /* temporary index table for matching stars */
  int *flag = 0;

  double tol2, dx, dy, dxy, dxys;
  double dmag, dmsum, dmave;
  double dxsum, dysum, dxave, dyave;
  double dxdev, dydev;

  VimosColumn *aXimaCol, *aYimaCol, *aXwldCol;
  VimosColumn *aYwldCol, *aMagCol, *agoffCol;
  VimosColumn *oXimaCol, *oYimaCol, *oMagCol;


  *matches = 0;

  tol2 = pow(radius, 2.);
  cpl_msg_debug(fctid, "Tolerance = %f pixel", radius);


  /*
   * Get the columns from the input source list
   */

  if (!(oXimaCol = findColInTab(sourceList, "X_IMAGE"))) {
    cpl_msg_error(fctid, "Column 'X_IMAGE' not found in source list!");
    return NULL;
  }

  if (!(oYimaCol = findColInTab(sourceList, "Y_IMAGE"))) {
    cpl_msg_error(fctid, "Column 'Y_IMAGE' not found in source list!");
    return NULL;
  }
  
  if (!(oMagCol = findColInTab(sourceList, "MAG"))) {
    cpl_msg_error(fctid, "Column 'MAG' not found in source list!");
    return NULL;
  }


  /*
   * Get columns from the reference catalog
   */

  if (!(aXimaCol = findColInTab(catalog, "X_IMAGE"))) {
    cpl_msg_error(fctid, "Column 'X_IMAGE' not found in reference catalog!");
    return NULL;
  }

  if (!(aYimaCol = findColInTab(catalog, "Y_IMAGE"))) {
    cpl_msg_error(fctid, "Column 'Y_IMAGE' not found in reference catalog!");
    return NULL;
  }

  if (!(aXwldCol = findColInTab(catalog, "RA"))) {
    cpl_msg_error(fctid, "Column 'RA' not found in reference catalog!");
    return NULL;
  }

  if (!(aYwldCol = findColInTab(catalog, "DEC"))) {
    cpl_msg_error(fctid, "Column 'DEC' not found in reference catalog!");
    return NULL;
  }

  if (!(aMagCol = findColInTab(catalog, "MAG"))) {
    cpl_msg_error(fctid, "Column 'MAG' not found in reference catalog!");
    return NULL;
  }

  if (!(agoffCol = findColInTab(catalog, "GOFF"))) {
    cpl_msg_error(fctid, "Column 'GOFF' not found in reference catalog!");
    return NULL;
  }


  /*
   * Get the number of rows in the input tables.
   */

  /* FIXME: This assumes that all columns have the same size. This should
   *        be checked!
   */

  n_ele1 = colGetSize(oXimaCol);
  n_ele2 = colGetSize(aXimaCol);


  /*
   * Allocate working arrays for the index table and the selection flag
   */

  sz = n_ele1 < n_ele2 ? n_ele1 : n_ele2;

  imatch0 = (int *)cpl_calloc((size_t)(2 * sz), sizeof(int));
  if (!imatch0) {
    cpl_msg_error(fctid, "Not enough memory!");
    return NULL;
  }

  
  if ((flag = (int *)cpl_calloc((size_t)n_ele2, sizeof(int))) == 0) {
    cpl_msg_error(fctid, "Not enough memory!");
    cpl_free(imatch0);

    return NULL;
  }


  /*
   * Object idendification starts here
   */

  nmatch0 = 0;
  dxsum = 0.0;
  dysum = 0.0;
  dmsum = 0.0;

  for (i = 0; i < n_ele1; i++) {
    dxys = -1.0;
    jbest = -1;

    for (j = 0; j < n_ele2; j++) {

      /*
       * flags is used to discard problematic object, objects close
       * to the border for instance.
       */

      if ((!(agoffCol->colValue->iArray[j])) && !(*(flag+j))) {
        dx = aXimaCol->colValue->dArray[j] - oXimaCol->colValue->dArray[i];
        dy = aYimaCol->colValue->dArray[j] - oYimaCol->colValue->dArray[i]; 
        dxy = pow(dx, 2.) + pow(dy, 2.);
        dmag = fabs((aMagCol->colValue->dArray[j] - 
                     oMagCol->colValue->dArray[i]));

        if (magTolerance1 > 0.0) {
          if (dxy < tol2 && dmag <= magTolerance1) {
            if (dxys < 0.0) {
              dxys = dxy;
              jbest = j;
            }
            else
              if (dxy < dxys) {
                dxys = dxy;
                jbest = j;
              }
          }
        }
        else {
          if (dxy < tol2) {
            if (dxys < 0.0) {
              dxys = dxy;
              jbest = j;
            }
            else
              if (dxy < dxys) {
                dxys = dxy;
                jbest = j;
              }
          }
        }
      }
    }

    if (jbest > -1) {
      dx = aXimaCol->colValue->dArray[jbest] - oXimaCol->colValue->dArray[i];  
      dy = aYimaCol->colValue->dArray[jbest] - oYimaCol->colValue->dArray[i];
      dxsum += dx;
      dysum += dy;
      dmsum += aMagCol->colValue->dArray[jbest] - oMagCol->colValue->dArray[i];

      *(imatch0+(2 * nmatch0)) = i;
      *(imatch0+(2 * nmatch0 +1)) = jbest;
      *(flag+jbest) = 1;
      nmatch0++;
    }
  }

  cpl_free(flag);

  if (nmatch0 < (minStars < 2 ? minStars : 2)) {
    cpl_msg_error(fctid, "Insufficient number of matches found [%d]", nmatch0);
    cpl_free(imatch0);

    return NULL;
  }
  cpl_msg_debug(fctid, "Found %d matches", nmatch0);

  if (nmatch0 > 1) {
    dxave = dxsum / nmatch0;
    dyave = dysum / nmatch0;
    dmave = dmsum / nmatch0;

    
    /*
     * Get the standard deviations of dx and dy from the mean
     */
    
    dxsum = 0.0;
    dysum = 0.0;
    dxdev = 0.0;
    dydev = 0.0;

    for (i = 0; i < nmatch0; i++) {
      dx = (aXimaCol->colValue->dArray[imatch0[2 * i + 1]]) - 
        (oXimaCol->colValue->dArray[imatch0[2 * i]]) - dxave;
      dy = (aYimaCol->colValue->dArray[imatch0[2 * i + 1]]) - 
        (oYimaCol->colValue->dArray[imatch0[2 * i]]) - dyave;
      dxsum += dx;
      dysum += dy;
      dxdev += dx*dx;
      dydev += dy*dy;
    }
    
    dxdev = sqrt((dxdev - pow(dxsum, 2.) / nmatch0) / (nmatch0 - 1));
    dydev = sqrt((dydev - pow(dysum, 2.) / nmatch0) / (nmatch0 - 1));

    /** FIXME:
     * (Paola) patch here: apply sigma clipping only when dxdev and dydev
     * are significant. If not we play with "numeric noise" and
     * many sources are excluded. This happened in "fitCO"
     * where many artificial sources are excluded due only to noise.
     */

   if(dxdev < MIN_DEVIATION)
     dxdev = MIN_DEVIATION;
   if(dydev < MIN_DEVIATION )
     dydev = MIN_DEVIATION;


    cpl_msg_debug(fctid, "Applying 2-sigma rejection: dxdev=%g; dydev=%g",
                dxdev, dydev);

    
    /*
     * Find best matches applying a Sigma clipping
     */
    
    if (!(imatch = (int *)cpl_calloc((size_t)(2 * nmatch0), sizeof(int)))) {
      cpl_msg_error(fctid, "Not enough memory!");
      cpl_free(imatch0);
      return NULL;
    }

    nmatch = 0;
    nrej = 0;

    for (i = 0; i < nmatch0; i++) {
      dx = (aXimaCol->colValue->dArray[imatch0[2 * i + 1]]) - 
        (oXimaCol->colValue->dArray[imatch0[2 * i]]) - dxave;
      dy = (aYimaCol->colValue->dArray[imatch0[2 * i + 1]]) - 
        (oYimaCol->colValue->dArray[imatch0[2 * i]]) - dyave;
      dmag = (aMagCol->colValue->dArray[imatch0[2 * i + 1]]) - 
        (oMagCol->colValue->dArray[imatch0[2 * i]]) - dmave;
      
      
      if (magTolerance2 > 0.0) {
        if (fabs(dx) <= sigmaClip * dxdev && fabs(dy) <= sigmaClip * dydev && 
            fabs(dmag) <= magTolerance2 ) {
          *(imatch+(2 * nmatch)) = imatch0[2 * i];
          *(imatch+(2 * nmatch +1)) = imatch0[2 * i + 1];
          nmatch++;
        }
        else
          nrej++;
      }
      else {
        if (fabs(dx) <= sigmaClip * dxdev && fabs(dy) <= sigmaClip * dydev) {
          *(imatch+(2 * nmatch)) = imatch0[2 * i];
          *(imatch+(2 * nmatch +1)) = imatch0[2 * i + 1];
          nmatch++;
        }
        else
          nrej++;
      }
    }

    if(nrej > 0)
      cpl_msg_debug(fctid, "Rejected %d pair(s)", nrej);

    cpl_free (imatch0);
    *matches = nmatch;
  }
  else {
    if (!(imatch = (int *)cpl_calloc((size_t)(2 * nmatch0), sizeof(int)))) {
      cpl_msg_error(fctid, "Not enough memory!");
      cpl_free(imatch0);

      return NULL;
    }

    *(imatch) = imatch0[0];
    *(imatch+1) = imatch0[1];
    *matches = 1;
  }  

  return imatch;

}
/**@}*/
