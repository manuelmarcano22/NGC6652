/* $Id: vmsextractor.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Date: 2013-03-25 11:43:04 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <pilmemory.h>
#include <pildfsconfig.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilstrutils.h>
#include <pilfileutils.h>
#include <pilutils.h>

#include "vmtable.h"
#include "vmimgutils.h"
#include "vmsextractor.h"
#include "cpl.h"


#define SEXTRACTOR_GROUP          "SExtractor"
#define SEXTRACTOR_COMMENT_CHARS  "#"

#define SEEING_DEFAULT (1.5)


/**
 * @brief vmsextractor Sextractor Interface
 *
 * The module @b VmSextractor provides a basic interface
 * implementation to allow for the usage of the SExtractor source
 * detection software from within the VIMOS DRS.
 */

/**@{*/

/**
 * @brief
 *   Expand a SExtractor file name and copy it to a buffer.
 *
 * @return The function returns the buffer if no error occurred, otherwise
 *   the return value is @c NULL.
 *
 * @param buffer  Memory area which is filled with the expanded path. It
 *                must be large enough.
 * @param path    Unexpanded path string.
 * @param n       Maximum number of characters copied to the buffer.
 *
 * The function checks if @em path is not @c NULL and has a non-zero
 * length. The string provided by @em path is assumed to be the
 * path to a file. It is expanded to an absolute path and the first
 * @em n characters are copied to the buffer given by @em buffer. The
 * target buffer must be large enough to hold the expanded path. The
 * buffer has to provide space for one additional character, the
 * terminating '\0'. On each call the buffer array is padded with zero's
 * assuming that the buffer's size is @b n + 1!
 */

const char *
sextGetFileName(char *buffer, const char *path, size_t n)
{

  char *s;

  if (!path || strlen(path) == 0)
    return 0;

  memset(buffer, 0, n + 1);

  s = cpl_strdup(pilFileExpandFilePath(path));

  if (strlen(s) > n) {
    cpl_free(s);
    return 0;
  }

  strncpy(buffer, s, n);
  cpl_free(s);

  return buffer;
  
}

/* FIXME:
 *  The following functions should actually print some explanatory error
 *  messages in addition to returning an error status. But note that
 *  in this case they may not be static any longer and should go to their
 *  own module. (RP)
 */

/**
 * @brief
 *   Write the SExtractor runtime configuration to a file.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param cfile  Pointer to the output stream.
 * @param image  Source detection image.
 *
 * The function writes the SExtractor runtime configuration to the
 * file @em cfile. The keyword-value pairs needed for writing
 * the configuration file are taken from the group @b SExtractor,
 * which must be present in the recipe database. The image @em image
 * is used to provide SExtractor with the necessary observation specific
 * information (e.g. seeing). The function does not handle any of the
 * SExtractor filename parameters, i.e. if one of these parameters,
 * @b CatalogName for instance, appears in the recipe database
 * it is simply ignored by this function. Filename parameters have to
 * be treated separately!
 *
 * If the magnitude zero point is not found in the header od the image
 * @em image it is taken from the recipe configuration database. In this
 * case the value is written to the image header.
 */

int sextSaveConfiguration(FILE *cfile, VimosImage *image)
{

  const char fctid[] = "sextSaveConfiguration";

  char *value;

  float zeroPoint, pixelScale;
  float seeingEstimate = 0.;

  double gain;


  /*
   * Check if SExtractor configuration group can be accessed.
   */

  if (!pilDfsDbGroupExists(SEXTRACTOR_GROUP))
    return EXIT_FAILURE;


  /*
   * The current values for pixel scale, detector gain, zeropoint, and
   * the seeing are set from the image header keywords.
   */

  if (readFloatDescriptor(image->descs, pilTrnGetKeyword("PixelScale"),
			  &pixelScale, NULL) != VM_TRUE)
    return EXIT_FAILURE;
  else
    fprintf(cfile, "PIXEL_SCALE %.3f\n", pixelScale);

  /* FIXME:
   *  The following is fully correct only for single port readouts. But
   *  since SExtractor allows only for one gain parameter currently
   *  the average gain is used in case of multi port readouts. Note
   *  also, that the meaning of CONAD and GAIN is exchanged for FORS.
   *  VIMOS does it correctly and ESO DET OUTi GAIN is actually what
   *  should be used here.
   */

  if ((gain = getMeanGainFactor(image)) < 0.)
    return EXIT_FAILURE;
  else 
    fprintf(cfile, "GAIN %.2f\n", gain);


  if (readFloatDescriptor(image->descs, pilTrnGetKeyword("MagZero"), 
                          &zeroPoint, NULL) != VM_TRUE) {
    value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "MagZeropoint");

    if (value) {
      fprintf(cfile, "MAG_ZEROPOINT %s\n", value);

      /* FIXME:
       *   Check if this is the right way to propagate the zeropoint
       *   value, or if it would be better to pass it back to the
       *   caller. (RP)
       */

      if (writeFloatDescriptor(&image->descs,
                               pilTrnGetKeyword("MagZero"), zeroPoint,
                               pilTrnGetComment("MagZero")) != VM_TRUE) {
        return EXIT_FAILURE;
      }
    }
    else
      return EXIT_FAILURE;
  }
  else
    fprintf(cfile, "MAG_ZEROPOINT %.3f\n", zeroPoint);


  if (readFloatDescriptor(image->descs, pilTrnGetKeyword("Seeing"),
			  &seeingEstimate, NULL) != VM_TRUE) {
      float seeingStart = 0.;
      float seeingEnd = 0.;
          
    if (readFloatDescriptor(image->descs, pilTrnGetKeyword("SeeingStart"),
			    &seeingStart, NULL) != VM_TRUE) {
      cpl_msg_warning(fctid, "Missing keyword `%s'",
                    pilTrnGetKeyword("SeeingStart"));
    }
    
    if (readFloatDescriptor(image->descs, pilTrnGetKeyword("SeeingEnd"),
			    &seeingEnd, NULL) != VM_TRUE) {
      cpl_msg_warning(fctid, "Missing keyword `%s'",
                    pilTrnGetKeyword("SeeingEnd"));
    }

    if (seeingStart > 0. && seeingEnd > 0.) {
        seeingEstimate = 0.5 * (seeingStart + seeingEnd);
    }
  }

  if (seeingEstimate <= 0.0) {
    cpl_msg_warning(fctid, "Invalid seeing value found (%.2f); using `%.2f' "
                  "instead", seeingEstimate, SEEING_DEFAULT);
    seeingEstimate = SEEING_DEFAULT;
  }
    
  fprintf(cfile, "SEEING_FWHM %.2f\n", seeingEstimate);


  /*
   * Setup parameters retrieved from the recipe database.
   */

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "AnalysisThresh")))
    fprintf(cfile, "ANALYSIS_THRESH %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "AssocData")))
    fprintf(cfile, "ASSOC_DATA %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "AssocParams")))
    fprintf(cfile, "ASSOC_PARAMS %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "AssocRadius")))
    fprintf(cfile, "ASSOC_RADIUS %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "AssocType")))
    fprintf(cfile, "ASSOC_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "AssocSelecType")))
    fprintf(cfile, "ASSOCSELEC_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "BackFilterSize")))
    fprintf(cfile, "BACK_FILTERSIZE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "BackSize")))
    fprintf(cfile, "BACK_SIZE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "BackType")))
    fprintf(cfile, "BACK_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "BackValue")))
    fprintf(cfile, "BACK_VALUE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "BackPhotoThick")))
    fprintf(cfile, "BACKPHOTO_THICK %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "BackPhotoType")))
    fprintf(cfile, "BACKPHOTO_TYPE %s\n", value);

  /* FIXME:
   *  Currently only the ASCII_HEAD output format is supported! Who wants
   *  to be volunteer for the remaining drivers? (RP)
   */

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "CatalogType"))) {
    if (strncmp(value, "ASCII_HEAD", 10))
      cpl_msg_warning(fctid, "SExtractor output format '%s' is not supported! "
		    "Format reset to ASCII_HEAD.", value);

    fprintf(cfile, "CATALOG_TYPE ASCII_HEAD\n");
  }

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "CheckImageType")))
    fprintf(cfile, "CHECKIMAGE_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "Clean")))
    fprintf(cfile, "CLEAN %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "CleanParam")))
    fprintf(cfile, "CLEAN_PARAM %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "DeblendMinCont")))
    fprintf(cfile, "DEBLEND_MINCONT %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "DeblendNthresh")))
    fprintf(cfile, "DEBLEND_NTHRESH %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "DetectMinArea")))
    fprintf(cfile, "DETECT_MINAREA %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "DetectThresh")))
    fprintf(cfile, "DETECT_THRESH %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "DetectType")))
    fprintf(cfile, "DETECT_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "Filter")))
    fprintf(cfile, "FILTER %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "FilterThresh")))
    fprintf(cfile, "FILTER_THRESH %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "FitsUnsigned")))
    fprintf(cfile, "FITS_UNSIGNED %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "FlagType")))
    fprintf(cfile, "FLAG_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "InterpMaxXlag")))
    fprintf(cfile, "INTERP_MAXXLAG %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "InterpMaxYlag")))
    fprintf(cfile, "INTERP_MAXYLAG %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "InterpType")))
    fprintf(cfile, "INTERP_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "MagGamma")))
    fprintf(cfile, "MAG_GAMMA %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "MaskType")))
    fprintf(cfile, "MASK_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "MemoryBufSize")))
    fprintf(cfile, "MEMORY_BUFSIZE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "MemoryObjStack")))
    fprintf(cfile, "MEMORY_OBJSTACK %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "MemoryPixStack")))
    fprintf(cfile, "MEMORY_PIXSTACK %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "PhotApertures")))
    fprintf(cfile, "PHOT_APERTURES %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "PhotAutoParams")))
    fprintf(cfile, "PHOT_AUTOPARAMS %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "PhotAutoApers")))
    fprintf(cfile, "PHOT_AUTOAPERS %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "PhotFluxFrac")))
    fprintf(cfile, "PHOT_FLUXFRAC %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "SaturLevel")))
    fprintf(cfile, "SATUR_LEVEL %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "ThreshType")))
    fprintf(cfile, "THRESH_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "VerboseType")))
    fprintf(cfile, "VERBOSE_TYPE %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "WeightGain")))
    fprintf(cfile, "WEIGHT_GAIN %s\n", value);

  if ((value = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "WeightType")))
    fprintf(cfile, "WEIGHT_TYPE %s\n", value);

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Write the SExtractor output setup to a file.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param pfile    Pointer the output stream.
 * @param columns  List of SExtractor output column labels.
 *
 * The function writes each SExtractor output column label appearing
 * in the list of column labels @em columns to the output 
 * configuration file @em pfile, using the SExtractor parameter
 * file format.
 */

int sextSaveParameters(FILE *pfile, SextParameter *columns)
{

  register size_t i;

  if (!pfile || !columns)
    return EXIT_FAILURE;

  i = 0;
  while (columns[i].name) {
    fprintf(pfile, "%s\n", columns[i].name);
    i++;
  }

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the name of the SExtractor neuronal network setup file.
 *
 * @return The function returns the filename, or @c NULL if the filename
 *   is not available.
 *
 * The function queries the recipe database for the entry
 * @b StarNnwName in the group @b SExtractor. If
 * the entry is found, a pointer to a static string holding the
 * absolute path is returned.
 */

const char *sextGetStarNnwName(void)
{

  const char *name = pilDfsDbGetString(SEXTRACTOR_GROUP, "StarNnwName");

  static char path[PATHNAME_MAX + 1] = "";


  return sextGetFileName(path, name, PATHNAME_MAX);

}


/**
 * @brief
 *   Get the name of the SExtractor star classifier setup file.
 *
 * @return The function returns the filename, or @c NULL if the filename
 *   is not available.
 *
 * The function queries the recipe database for the entry
 * @em FilterName in the group @b SExtractor. If
 * the entry is found, a pointer to a static string holding the
 * absolute path is returned.
 */

const char *sextGetFilterName(void)
{

  const char *name = pilDfsDbGetString(SEXTRACTOR_GROUP, "FilterName");

  static char path[PATHNAME_MAX + 1] = "";


  return sextGetFileName(path, name, PATHNAME_MAX);

}


/**
 * @brief
 *   Get the name of the SExtractor association setup file.
 *
 * @return The function returns the filename, or @c NULL if the filename
 *   is not available.
 *
 * The function queries the recipe database for the entry
 * @b AssocName in the group @b SExtractor. If the entry is
 * found, a pointer to a static string holding the absolute
 * path is returned.
 */

const char *sextGetAssocName(void)
{

  const char *name = pilDfsDbGetString(SEXTRACTOR_GROUP, "AssocName");

  static char path[PATHNAME_MAX + 1] = "";


  return sextGetFileName(path, name, PATHNAME_MAX);

}


/**
 * @brief
 *   Get the name of the SExtractor check image file.
 *
 * @return The function returns the filename, or @c NULL if the filename
 *   is not available.
 *
 * The function queries the recipe database for the entry
 * @b CheckImageName in the group @b SExtractor. If
 * the entry is found, a pointer to a static string holding the
 * absolute path is returned.
 */

const char *sextGetCheckImageName(void)
{

  const char *name = pilDfsDbGetString(SEXTRACTOR_GROUP, "CheckImageName");

  static char path[PATHNAME_MAX + 1] = "";


  return sextGetFileName(path, name, PATHNAME_MAX);

}


/**
 * @brief
 *   Get the name of the SExtractor flag image file.
 *
 * @return The function returns the filename, or @c NULL if the filename
 *   is not available.
 *
 * The function queries the recipe database for the entry
 * @b FlagImage in the group @b SExtractor. If the entry is
 * found, a pointer to a static string holding the absolute path
 * is returned.
 */

const char *sextGetFlagImageName(void)
{

  const char *name = pilDfsDbGetString(SEXTRACTOR_GROUP, "FlagImage");

  static char path[PATHNAME_MAX + 1] = "";


  return sextGetFileName(path, name, PATHNAME_MAX);

}


/**
 * @brief
 *   Get the name of the SExtractor weight image file.
 *
 * @return The function returns the filename, or @c NULL if the filename
 *   is not available.
 *
 * The function queries the recipe database for the entry
 * @b WeightImage in the group @b SExtractor. If
 * the entry is found, a pointer to a static string holding the
 * absolute path is returned.
 */

const char *sextGetWeightImageName(void)
{

  const char *name = pilDfsDbGetString(SEXTRACTOR_GROUP, "WeightImage");

  static char path[PATHNAME_MAX + 1] = "";


  return sextGetFileName(path, name, PATHNAME_MAX);

}


/**
 * @brief
 *   Get the pathname of the SExtractor executable.
 *
 * @return The absolute pathname to the sextractor executable file, or
 *   @c NULL if the pathname could not be retrieved.
 *
 * The function queries the recipe database for the entry
 * @b SExtractor in the group @b SExtractor. If
 * the entry is found, a pointer to a static string, holding the
 * absolute path to the SExtractor executable, is returned.
 */

const char *sextGetSextractorPath(void)
{

  static char sextractor[PATHNAME_MAX + 1] = "";

  const char *path = pilDfsDbGetString(SEXTRACTOR_GROUP, "SExtractor");


  return sextGetFileName(sextractor, path, PATHNAME_MAX);

}


/**
 * @brief
 *   Get the execution time limit for running SExtractor.
 *
 * @return The function returns the execution time limit.
 *
 * The function queries the recipe database for the entry @b Timeout
 * in the group @b SExtractor. The timeout limit will either be a
 * positive integer to be interpreted as the execution time limit in 
 * seconds or 0.
 */

time_t sextGetExecutionTimeLimit(void)
{

  long t = pilDfsDbGetLong(SEXTRACTOR_GROUP, "Timeout", 0);

  if (t < 0)
    t = 0;

  return t;

}


/**
 * @brief
 *   Convert SExtractor output catalog into a table object.
 *
 * @return The function returns the table object if no error occurred,
 *   otherwise the return value is @c NULL.
 *
 * @param catalog  Name of the SExtractor catalog file.
 * @param columns  Catalog column properties (Not yet available!).
 * @param image    Image used for source detection.
 *
 * The function extracts from the SExtractor output catalog given by
 * @em catalog the columns needed to build a table object with
 * the detected objects. The table object is created and filled with the
 * data from catalog columns listed in @em columns (not yet implemented!),
 * applying the necessary unit conversions. The table header keywords are
 * updated from the source detection image given by @em image.
 */

VimosTable *sextConvertCatalog(const char *catalog, SextParameter *columns)
{

  char catLine[LINE_LENGTH_MAX + 1];

  int n;

  char *window;
  int   x, y, sx, sy, ex, ey;

  size_t i, j, lineCount, winLineCount;

  FILE *catFile;

  VimosTable *table = newTable();
  VimosTable *winTable;

  VimosColumn *tableColumn;
  VimosColumn *winTableColumn;

  VimosColumn *x_image;
  VimosColumn *y_image;



  if (!table)
    return 0;

  /*
   * Get SExtractor window
   */

  window = (char *)pilDfsDbGetString(SEXTRACTOR_GROUP, "Window");
  sscanf(window, "%d,%d,%d,%d", &sx, &sy, &ex, &ey);

  
  /*
   * Count the valid records in the SExtractor output file.
   */

  if (!(catFile = fopen(catalog, "r"))) {
    deleteTable(table);
    return 0;
  }

  lineCount = 0;

  while (fgets(catLine, LINE_LENGTH_MAX, catFile))
    if (!strempty(catLine, SEXTRACTOR_COMMENT_CHARS))
      lineCount++;

  if (ferror(catFile)) {
    deleteTable(table);
    fclose(catFile);
    return 0;
  }


  /*
   * Create the table columns
   */

  i = 0;
  while (columns[i].name) {
    switch (columns[i].type) {
    case SEXT_COLUMN_INT:
      tableColumn = newIntColumn(lineCount, columns[i].name);
      break;

    case SEXT_COLUMN_FLOAT:
      tableColumn = newDoubleColumn(lineCount, columns[i].name);
      break;

    default:
      deleteTable(table);
      fclose(catFile);
      return 0;
      break;
    }


    /*
     * Insert the created column into the table
     */

    if (!tableColumn) {
      deleteTable(table);
      fclose(catFile);
      return 0;
    }
    else 
      tblAppendColumn(table, tableColumn);

    i++;
  }


  /*
   * Read the catalog into the table. Note that here we use the
   * knowledge about the internal layout of the table structure
   * and the order in which the columns were created!
   */

  rewind(catFile);
  clearerr(catFile);

  i = 0;
  lineCount = 0;

  while (fgets(catLine, LINE_LENGTH_MAX, catFile)) {
    if (!strempty(catLine, SEXTRACTOR_COMMENT_CHARS)) {
      char *s;

      /*
       * For the current line read all column entries
       */

      assert(table->numColumns > 0);

      s = catLine;

      for (i = 0; i < (size_t)table->numColumns; i++) {
	int ivalue;
	double dvalue;

	s = strskip(s, isspace);

	switch (columns[i].type) {
	case SEXT_COLUMN_INT:
	  if ((n = sscanf(s, "%d", &ivalue)) == 1)
	    tblSetIntValue(table, columns[i].name, lineCount, ivalue);

	  break;

	case SEXT_COLUMN_FLOAT:
	  if ((n = sscanf(s, "%lf", &dvalue)) == 1)
	    tblSetDoubleValue(table, columns[i].name, lineCount, dvalue);

	  break;

	default:
	  n = 0;
	  break;
	}

	/*
	 * Stop on conversion error
	 */

	if (n != 1) {
	  deleteTable(table);
	  fclose(catFile);
	  return 0;
	}


	/*
	 * Skip entry just read, assuming that entries are searated
	 * by white space characters in the sextractor output file.
	 */

	while(!isspace(*s))
	  s++;

      }

      lineCount++;

    }
  }

  if (ferror(catFile)) {
    deleteTable(table);
    fclose(catFile);
    return 0;
  }
  else
    fclose(catFile);

  x_image = findColumn(table->cols, "X_IMAGE");
  y_image = findColumn(table->cols, "Y_IMAGE");

  winLineCount = 0;

  for (i = 0; i < lineCount; i++) {
    x = x_image->colValue->dArray[i];
    y = y_image->colValue->dArray[i];
    if (x >= sx && x < ex && y >= sy && y < ey)
      winLineCount++;
  }

  if (lineCount == winLineCount)
    return table;

  winTable = newTable();

  /*
   * Create the windowed table columns
   */
   
  i = 0;
  while (columns[i].name) {
    switch (columns[i].type) {
    case SEXT_COLUMN_INT:
      tableColumn = newIntColumn(winLineCount, columns[i].name);
      break;

    case SEXT_COLUMN_FLOAT:
      tableColumn = newDoubleColumn(winLineCount, columns[i].name);
      break;

    default:
      deleteTable(table);
      fclose(catFile);
      return 0;
      break;
    }


    /*
     * Insert the created column into the table
     */

    if (!tableColumn) {
      deleteTable(table);
      deleteTable(winTable);
      fclose(catFile);
      return 0;
    }
    else
      tblAppendColumn(winTable, tableColumn);

    i++;
  }

  for (i = 0; i < (size_t)winTable->numColumns; i++) {

    tableColumn = findColumn(table->cols, columns[i].name);
    winTableColumn = findColumn(winTable->cols, columns[i].name);

    winLineCount = 0;
    for (j = 0; j < lineCount; j++) {
      x = x_image->colValue->dArray[j];
      y = y_image->colValue->dArray[j];
      if (x >= sx && x < ex && y >= sy && y < ey) {
        switch (columns[i].type) {
        case SEXT_COLUMN_INT:
          winTableColumn->colValue->iArray[winLineCount] = 
                                       tableColumn->colValue->iArray[j];
          break;

        case SEXT_COLUMN_FLOAT:
          winTableColumn->colValue->dArray[winLineCount] = 
                                       tableColumn->colValue->dArray[j];
          break;

        default:
          break;
        }
        winLineCount++;

      }
    }
  }

  deleteTable(table);

  return winTable;

}
/**@}*/
