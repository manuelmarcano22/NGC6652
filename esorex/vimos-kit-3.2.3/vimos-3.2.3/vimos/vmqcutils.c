/* $Id: vmqcutils.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <fitsio.h>
#include <fitsio2.h>

#include <pilmemory.h>
#include <pilstrutils.h>
#include <piltranslator.h>
#include <pilerrno.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilqc.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmqcutils.h"
#include "cpl.h"


#define COMMENT_LENGTH    (80)


/**
 * @name vmqcutils Quality Control Utilities
 *
 * The module collects utility functions for quality control operations.
 */

/**@{*/

/*
static int myqfits_is_int(const char * s)
{
    regex_t re_int ;
    int     status ;

    if (s==NULL) return 0 ;
    if (s[0]==0) return 0 ;
    if (regcomp(&re_int, &regex_int[0], REG_EXTENDED|REG_NOSUB)!=0) {
        printf("internal error: compiling int rule");
        abort();
    }
    status = regexec(&re_int, s, 0, NULL, 0) ;
    regfree(&re_int) ;
    return (status) ? 0 : 1 ;
}

static int myqfits_is_float(const char * s)
{
    regex_t re_float;
    int     status ;

    if (s==NULL) return 0 ;
    if (s[0]==0) return 0 ;
    if (regcomp(&re_float, &regex_float[0], REG_EXTENDED|REG_NOSUB)!=0) {
        printf("internal error: compiling float rule");
        abort();
    }
    status = regexec(&re_float, s, 0, NULL, 0) ;
    regfree(&re_float) ;
    return (status) ? 0 : 1 ;
}

static int myqfits_is_boolean(const char * s)
{
    if (s==NULL) return 0 ;
    if (s[0]==0) return 0 ;
    if ((int)strlen(s)>1) return 0 ;
    if (s[0]=='T' || s[0]=='F') return 1 ;
    return 0 ;
}

static int myqfits_is_complex(const char * s)
{
    regex_t re_cmp ;
    int     status ;

    if (s==NULL) return 0 ;
    if (s[0]==0) return 0 ;
    if (regcomp(&re_cmp, &regex_cmp[0], REG_EXTENDED|REG_NOSUB)!=0) {
        printf("internal error: compiling complex rule");
        abort();
    }
    status = regexec(&re_cmp, s, 0, NULL, 0) ;
    regfree(&re_cmp) ;
    return (status) ? 0 : 1 ;
}

#define PRETTY_STRING_STATICBUFS    8
*/
/*----------------------------------------------------------------------------*/
/**
  @brief    Clean out a FITS string value.
  @param    s pointer to allocated FITS value string.
  @return   pointer to statically allocated character string

  From a string FITS value like 'marvin o''hara', remove head and tail
  quotes, replace double '' with simple ', trim blanks on each side,
  and return the result in a statically allocated area.

  Examples:

  - ['o''hara'] becomes [o'hara]
  - ['  H    '] becomes [H]
  - ['1.0    '] becomes [1.0]

 */
/*----------------------------------------------------------------------------*/

/*
static char * myqfits_pretty_string(const char * s)
{
    static char     pretty_buf[PRETTY_STRING_STATICBUFS][81] ;
    static int      flip=0 ;
    char        *   pretty ;
    int             i,j ;

    if (s==NULL) return NULL ;

    pretty = pretty_buf[flip];
    flip++ ;
    if (flip==PRETTY_STRING_STATICBUFS)
        flip=0 ;

    pretty[0] = (char)0 ;
    if (s[0]!='\'') return (char *)s ;

    i=1 ;
    j=0 ;
    while (s[i]==' ') {
        if (i==(int)strlen(s)) break ;
        i++ ;
    }
    if (i>=(int)(strlen(s)-1)) return pretty ;
    while (i<(int)strlen(s)) {
        if (s[i]=='\'') {
            i++ ;
        }
        pretty[j]=s[i];
        i++ ;
        j++ ;
    }
    pretty[j+1]=(char)0;
    j = (int)strlen(pretty)-1;
    while (pretty[j]==' ') j-- ;
    pretty[j+1]=(char)0;
    return pretty;
}
#undef PRETTY_STRING_STATICBUFS
*/

/*
static void keytuple2str(
        char        *   line,
        const char  *   key,
        const char  *   val,
        const char  *   com)
{
    int     len ;
    int     hierarch = 0 ;
    char    cval[81];
    char    cval2[81];
    char    cval_q[81];
    char    ccom[81];
    char    safe_line[512];
    int     i, j ;

    if (line==NULL || key==NULL) return ;

    memset(line, ' ', 80);
    if (key==NULL) return ;

    if (!strcmp(key, "END")) {
        sprintf(line, "END") ;
        return ;
    }
    if (!strcmp(key, "HISTORY") ||
        !strcmp(key, "COMMENT") ||
        !strncmp(key, "        ", 8)) {
        sprintf(line, "%s ", key);
        if (val==NULL) return ;

        len = strlen(val);
        if (len>72) len=72 ;
        strncpy(line+8, val, len);
        return ;
    }

    if (val==NULL) cval[0]=(char)0;
    else if (strlen(val)<1) cval[0]=(char)0;
    else strcpy(cval, val);

    if (com==NULL) strcpy(ccom, "no comment");
    else strcpy(ccom, com);

    if (!strncmp(key, "HIERARCH", 8)) hierarch ++ ;

    if (myqfits_is_int(cval) ||
            myqfits_is_float(cval) ||
            myqfits_is_boolean(cval) ||
            myqfits_is_complex(cval)) {
        if (hierarch) sprintf(safe_line, "%-29s= %s / %s", key, cval, ccom);
        else sprintf(safe_line, "%-8.8s= %20s / %-48s", key, cval, ccom);
        strncpy(line, safe_line, 80);
        line[80]=(char)0;
        return ;
    }

    if (cval[0]==0) {
        if (hierarch) {
            sprintf(safe_line, "%-29s=                    / %s", key, ccom);
        } else {
        sprintf(safe_line, "%-8.8s=                      / %-48s", key, ccom);
        }
        strncpy(line, safe_line, 80);
        line[80]=(char)0;
        return ;
    }

    memset(cval_q, 0, 81);
    strcpy(cval2, myqfits_pretty_string(cval));
    j=0 ;
    i=0 ;
    while (cval2[i] != (char)0) {
        if (cval2[i]=='\'') {
            cval_q[j]='\'';
            j++ ;
            cval_q[j]='\'';
        } else {
            cval_q[j] = cval2[i];
        }
        i++ ;
        j++ ;
    }

    if (hierarch) {
        sprintf(safe_line, "%-29s= '%s' / %s", key, cval_q, ccom);
        if (strlen(key) + strlen(cval_q) + 3 >= 80)
            safe_line[79] = '\'';
    } else {
        sprintf(safe_line, "%-8.8s= '%-8s' / %s", key, cval_q, ccom);
    }
    strncpy(line, safe_line, 80);

    line[80]=(char)0;
    return ;
}
*/


/**
 * @memo
 *   Select images having a consistent signal distribution from a list.
 *
 * @return The number of selected frames if no error occured, otherwise
 *   the return value is 0 and the error flag @c pilErrno is set to
 *   @c P_EGENERIC.
 *
 * @param imageList   List of images.
 * @param imageNoise  Array of image noise values.
 * @param imageCount  Number of images in the list.
 * @param kappa       Kappa-Sigma scale factor.
 *
 * @doc
 *   The function compares the intensity distribution, i.e. the illumination
 *   pattern, of the images provided by \textbf{imageList}. The number
 *   of images in the list has to be provided through \textbf{imageCount}.
 *   Images showing an intensity distribution which is inconsistent with the 
 *   majority of the input images are moved to the end of \textbf{imageList}.
 *   Note that, if not at least 2 images show a consistent intensity
 *   distribution the returned number of consistent images is 0.
 *
 *   The comparison of the signal distribution is done in the following way:
 *   For all possible pairs of images the difference image is computed and
 *   the median value of the differences is determined. The deviation of the
 *   median difference from zero is then compared with \textbf{kappa}
 *   times the expected noise for the difference image. The expected noise
 *   is calculated from the noise found in the individual images which has
 *   to be provided through the array \textbf{imageNoise}. Images having a
 *   median difference smaller than the tolerance are considered consistent.
 *   
 *   The input list is rearranged so that all images which are consistent
 *   are found in the beginning of \textbf{imageList}. Images which are
 *   considered inconsistent are moved to the end of the list.
 *
 *   In case there are no images with a consistent signal distribution in
 *   input list \textbf{imageList} the function returns 0 and sets the error
 *   flag @c pilErrno to @c P_ENOERR. Instead, if an internal error occured
 *   the return value is still 0, but the error flag is set to @c P_EGENERIC.
 *
 *   FIXME: Instead of moving images within the list a new list object
 *   holding the accepted frames should be created. Requires an image
 *   list object!
 *
 * @author P. Sartoretti, R. Palsa
 */

size_t
qcSelectConsistentImages(VimosImage **imageList, float *imageNoise,
                         size_t imageCount, double kappa)
{

  register int i, j, k;

  int acceptedCount = 0, rejectedCount = 0, finalCount = 0;

  size_t imageSize = imageList[0]->xlen * imageList[0]->ylen;

  double sigma;

  VimosMatrix *offset, *threshold;

  VimosImage *referenceImage, *differenceImage;
  VimosImage **acceptedList = NULL;
  VimosImage **rejectedList = NULL;
  VimosImage **finalList = NULL;


  pilErrno = P_ENOERR;


  /*
   * At least 2 images must be present in the input.
   */

  if (imageCount < 2) {
    pilErrno = P_EGENERIC;
    return 0;
  }


  /*
   * Build the symmetric matrix of selection thresholds from the
   * expected noise of the difference image.
   */

  if (!(threshold = newMatrix(imageCount, imageCount))) {
    pilErrno = P_EGENERIC;
    return 0;
  }
  else
    for (i = 0; i < threshold->nr; i++)
      for (j = i + 1; j < threshold->nc; j++) {

        /*
         * Expected uncertainty of the difference of an image pair.
         * Simple error propagation.
         */

        sigma = sqrt(ipow(imageNoise[i], 2) + ipow(imageNoise[j], 2));

        k = i * threshold->nc + j;
        threshold->data[k] = kappa * sigma;
        threshold->data[j * threshold->nc + i] = threshold->data[k];

      }

      
  /*
   * Loop over all possible image pairs (don't care about the order) and
   * compute the median pixel value of the difference image. Build a matrix
   * of the absolute values of the differences.
   */

  if (!(offset = newMatrix(imageCount, imageCount))) {
    deleteMatrix(threshold);

    pilErrno = P_EGENERIC;
    return 0;
  }

  for (i = 0; (size_t)i < imageCount; i++) {
    referenceImage = imageList[i];

    for (j = i + 1; (size_t)j < imageCount; j++) {
      differenceImage = imageArith(imageList[j], referenceImage, VM_OPER_SUB);

      if (!differenceImage) {
        pilErrno = P_EGENERIC;
        return 0;
      }

      /*
       * The absolute value is needed, since we look for the typical
       * offset value and not for the actual median of the differences.
       */

      for (k = 0; (size_t)k < imageSize; k++)
        differenceImage->data[k] = fabs(differenceImage->data[k]);

      k = i * imageCount + j;
      offset->data[k] = imageMean(differenceImage);
      offset->data[j * imageCount + i] = offset->data[k];

      deleteImage(differenceImage);
    }
  }


  /*
   * Analyze the matrix of the median offsets, i.e. compare the offsets
   * from zero with the expected and scaled noise of the difference image.
   * Determine the matrix row with the largest number of images satisfying
   * the selection criteria.
   */

  acceptedList = (VimosImage **)cpl_calloc(imageCount, sizeof(VimosImage *));
  rejectedList = (VimosImage **)cpl_calloc(imageCount, sizeof(VimosImage *));
  finalList = (VimosImage **)cpl_calloc(imageCount, sizeof(VimosImage *));

  if (!acceptedList || !rejectedList || !finalList) {
    deleteMatrix(offset);
    deleteMatrix(threshold);

    if (acceptedList)
      cpl_free(acceptedList);

    if (rejectedList)
      cpl_free(rejectedList);

    if (finalList)
      cpl_free(finalList);
    
    pilErrno = P_EGENERIC;
    return 0;
  }

  for (i = 0; i < offset->nr; i++) {
    acceptedCount = 0;
    rejectedCount = 0;

    for (j = 0; j < offset->nc; j++) {
      k = i * offset->nc + j;
      if (offset->data[k] <= threshold->data[k]) {
        acceptedList[acceptedCount] = imageList[j];
        acceptedCount++;
      }
      else {
        rejectedList[rejectedCount] = imageList[j];
        rejectedCount++;
      }
    }

    if (acceptedCount > finalCount) {
      finalCount = acceptedCount;
      memcpy(finalList, acceptedList, acceptedCount * sizeof(VimosImage *));
      memcpy(finalList + acceptedCount, rejectedList, 
             rejectedCount * sizeof(VimosImage *));
    }
  }

  deleteMatrix(offset);
  deleteMatrix(threshold);


  /*
   * Rearrange the images in the input list.
   */

  for (i = 0; (size_t)i < imageCount; i++)
    imageList[i] = finalList[i];
    

  /*
   * Cleanup
   */

  cpl_free(acceptedList);
  cpl_free(rejectedList);
  cpl_free(finalList);

  return finalCount;

}


/**
 * @brief
 *   Compare the bias level with the nominal value.
 *
 * @return The function returns @c EXIT_SUCCESS if the quality check was
 *   successful, otherwise the return value is @c EXIT_FAILURE.
 *
 * @param bias      Input (master) bias frame.
 * @param mbias     Reference master bias frame.
 * @param maxDev    Maximum allowed deviation from the nominal bias level.
 * @param warnOnly  Only issue a warning if the check fails.
 * @param rCalc     Flag to use/recalculate the bias level.
 *
 * The function determines the median level of the input (master) bias
 * image @em bias. The way the median level is determined can be
 * controlled via the flag @em rCalc. Passing a non-zero flag value
 * forces the internal recalculation of the median image level. If
 * @em rCalc equals zero, the function first checks if the property
 * @b DataMedian exists. If it is found the value is used as median level.
 * If it is not present the function calculates the median bias level
 * internally.
 * 
 * If the median level was calculated by the function the property
 * @b DataMedian is updated.
 * 
 * The median level of the input (master) bias frame is then compared
 * with the nominal bias value, which is taken from the property
 * @b DataMedian of the reference master bias @em mbias. If the
 * difference between the median level of the input image and the
 * nominal bias value is larger than the maximum allowed difference
 * in terms of sigma @em maxDev the function returns with an error and sets
 * the error variable to @c P_ENOERR, indicating that only the check failed.
 * If the flag @b warnOnly is not set to zero the function only issues a
 * warning message but terminates successfully.
 * 
 * If the input image meets the quality control criterium the property
 * @b BiasOffset is set to the difference between the median bias level
 * and the nominal value. The nominal bias value the median level was
 * compared to is stored in the property @b BiasLevel. If the input image
 * does not meet the criterium the header is not updated.
 * 
 * In case of any other kind of failure, the error variable is set to
 * @c P_EGENERIC.
 *
 * @author R. Palsa
 */

int
qcCheckBiasLevel(VimosImage *bias, VimosImage *mbias, double maxDev,
                 unsigned int warnOnly, unsigned int rCalc)
{

  const char fctid[] = "qcCheckBiasLevel";

  char comment[COMMENT_LENGTH];

  double median, nominal, offset, sigma;



  /*
   * Determine median bias level. Use the header entry if computation
   * is not forced or the header entry does not exist. If the median
   * is computed write it back to the header.
   */

  
  if (rCalc || readDoubleDescriptor(bias->descs, pilTrnGetKeyword("BiasLevel"),
                                    &median, comment) == VM_FALSE) {
    cpl_msg_info(fctid, "Calculating image median ...");
    median = imageMedian(bias);
    writeFloatDescriptor(&(bias->descs), pilTrnGetKeyword("BiasLevel"), median,
                         "Median bias level");
  }
  else 
    cpl_msg_info(fctid, "Retrieving image median from header (%s) ...", 
               pilTrnGetKeyword("BiasLevel"));

  cpl_msg_info(fctid, "Median bias level: %.4f", median);


  /*
   * Read nominal bias value from the reference image.
   */

  pilErrno = P_ENOERR;

  if (readDoubleDescriptor(mbias->descs, pilTrnGetKeyword("DataMedian"),
                           &nominal, comment) == VM_FALSE) {
    cpl_msg_error(fctid, "Keword '%s' not found in master bias",
                pilTrnGetKeyword("DataMedian"));

    pilErrno = P_EGENERIC;

    return EXIT_FAILURE;
  }
  else 
    cpl_msg_info(fctid, "Nominal bias level: %.4f", nominal);


  /*
   * Compute offset from the nominal bias constant and check if the image
   * median is within the tolerance.
   */

  offset = median - nominal;

  /*
   * Do not use a call to imageMedSigma() here to speed up
   * computation, otherwise the median pixel value is recomputed.
   */

  sigma = maxDev * imageAverageDeviation(bias, median);

  if (fabs(offset) > sigma) {
    if (!warnOnly) {
      cpl_msg_error(fctid, "Median bias level offset exceeds maximum "
                  "tolerance value of %.2f sigma (%.4f)!", maxDev, sigma);

      pilErrno = P_ENOERR;

      return EXIT_FAILURE;
    }
    else 
      cpl_msg_warning(fctid, "Median bias level offset exceeds maximum "
                    "tolerance value of %.2f sigma (%.4f)!", maxDev, sigma);
  }
  else
    cpl_msg_info(fctid, "Median bias level within tolerance interval "
               "%.4f +/- %.4f (%.2f sigma)", nominal, sigma, maxDev);

  /*
   * Update the header.
   */

  writeFloatDescriptor(&(bias->descs), pilTrnGetKeyword("BiasOffset"),
                       offset, "Offset from nominal bias level");

  writeFloatDescriptor(&(bias->descs), pilTrnGetKeyword("BiasLevel"), nominal,
                       pilTrnGetComment("BiasLevel"));

  writeFloatDescriptor(&(bias->descs), pilTrnGetKeyword("BiasOffset"), offset,
                       pilTrnGetComment("BiasOffset"));

  return EXIT_SUCCESS;

}


/**
 * @memo
 *   Compare the dark level with the nominal value.
 *
 * @return The function returns EXIT_SUCCESS if the quality check was
 *   successful, otherwise the return value is EXIT_FAILURE.
 *
 * @param dark      Input (master) dark frame.
 * @param ccdTable  CCD table.
 * @param maxDev    Maximum allowed deviation from the nominal dark level.
 * @param warnOnly  Only issue a warning if the check fails.
 * @param rCalc     Flag to use/recalculate the dark level.
 *
 * @doc
 *   The function determines the median level of the input (master) dark
 *   image \textbf{dark}. The way the median level is determined can be
 *   controlled via the flag \textbf{rCalc}. Passing a non-zero flag value
 *   forces the internal recalculation of the median image level. If
 *   \textbf{rCalc} equals zero, the function first checks if the
 *   descriptor \textit{DarkLevel} exists. If it is found the descriptor
 *   value is used as median level. If it is not present the function
 *   calculates the median dark level internally.
 *
 *   If the median level was calculated by the function the descriptor
 *   \textit{DarkLevel} is updated.
 *
 *   The median level of the input (master) dark frame is then compared
 *   with the nominal dark constant, which is taken from the descriptor
 *   \textit{DarkLevel} of the CCD table \textbf{ccdTable}. If the
 *   difference between the median level of the input image and the nominal
 *   dark constant is larger than the maximum allowed difference in terms of
 *   sigma \textbf{maxDev} the function fails but sets the error variable to
 *   P_ENOERR, indicating that only the check failed. If the flag
 *   \textbf{warnOnly} is not set to zero the function only issues a
 *   warning message but terminates successfully.
 *
 *   If the input image meets the quality control criterium the difference
 *   between the median dark level and the nominal value is written to the
 *   descriptor \textit{DarkOffset}. If the input image does not meet the
 *   criterium the header is not updated.
 *
 *   In case of any other kind of failure, the error variable is set to
 *   P_EGENERIC.
 *
 * @author R. Palsa, C. Izzo
 */

int
qcCheckDarkLevel(VimosImage *dark, VimosTable *ccdTable, double maxDev,
                 unsigned int warnOnly, unsigned int rCalc)
{

  const char fctid[] = "qcCheckDarkLevel";

  char      *keyname = cpl_strdup(pilTrnGetKeyword("DarkLevel"));
  char       comment[COMMENT_LENGTH];

  double     median, nominal, offset, sigma;


  /*
   * Reset error flag
   */

  pilErrno = P_ENOERR;


  /*
   * Determine median dark level. Use the header entry if computation
   * is not forced or the header entry does not exist. If the median
   * is computed write it back to the header.
   */


  if (rCalc || readDoubleDescriptor(dark->descs, keyname, &median, comment)
      == VM_FALSE) {
    cpl_msg_info(fctid, "Calculating image median...");
    median = imageMedian(dark);
    writeFloatDescriptor(&(dark->descs), keyname, median,
                         "Median dark level");
  }
  else
    cpl_msg_info(fctid, "Retrieving image median from header (%s)...", keyname);

  cpl_msg_info(fctid, "Median dark level: %.4f", median);


  /*
   * Read nominal dark constant from the CCD table.
   */

  if (readDoubleDescriptor(ccdTable->descs, keyname, &nominal, comment)
      == VM_FALSE) {
    cpl_msg_error(fctid, "Invalid CCD table! Descriptor '%s' not found",
                keyname);

    pilErrno = P_EGENERIC;
    cpl_free(keyname);

    return EXIT_FAILURE;
  }

  cpl_msg_info(fctid, "Nominal dark level: %.4f", nominal);


  /*
   * Compute offset from the nominal dark constant and check if the image
   * median is within the tolerance.
   */

  offset = median - nominal;
  sigma = maxDev * imageAverageDeviation(dark, median);

  if (fabs(offset) > sigma) {
    if (!warnOnly) {
      cpl_msg_error(fctid, "Median dark level offset exceeds maximum "
                  "tolerance value of %.2f sigma (%.4f)!", maxDev, sigma);

      pilErrno = P_ENOERR;
      cpl_free(keyname);

      return EXIT_FAILURE;
    }
    else
      cpl_msg_warning(fctid, "Median dark level offset exceeds maximum "
                    "tolerance value of %.2f sigma (%.4f)!", maxDev, sigma);
  }
  else
    cpl_msg_info(fctid, "Median dark level within tolerance interval "
               "%.4f +/- %.4f (%.2f sigma)", nominal, sigma, maxDev);

  /*
   * Update the header.
   */

  writeFloatDescriptor(&(dark->descs), pilTrnGetKeyword("DarkOffset"),
                       offset, "Offset from nominal dark level");

  /*
   * Cleanup
   */

  cpl_free(keyname);

  return EXIT_SUCCESS;

}


/**
 * @memo
 *   Copy a descriptor value to the currently active QC1 PAF object.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param header  Pointer to a descriptor header.
 * @param name    Descriptor name.
 * @param unit    Optional unit to be associated to descriptor value.
 * @param comment Optional comment to be associated to descriptor value.
 *
 * @doc
 *   A descriptor with the specified name is searched in the header. 
 *   Its type is determined, then its value is read with the appropriate 
 *   interface. Only VM_INT, VM_FLOAT, VM_DOUBLE and VM_STRING descriptors 
 *   are supported. From the descriptor name the corresponding PAF keyword 
 *   name is derived by removing any "ESO " at descriptor name beginning, 
 *   and replacing blanks with dots (e.g., "ESO TPL ID" becomes "TPL.ID").
 *   Finally, the new PAF keyword, with the same type as the descriptor, 
 *   is written to the currently active QC1 PAF object. Note that before 
 *   calling this funtion a QC1 PAF object must be created with a call to 
 *   pilQcGroupStart().
 *
 * @author C. Izzo
 */

int
qcCopyValue(VimosDescriptor *header, const char *name, const char *unit,
            const char *comment)
{

  const char fctid[] = "qcCopyValue";

  VimosDescriptor *desc;
  int              status = EXIT_FAILURE;
  int              i;
  char            *pafName;
  char            *keep;
  char            *pos;
  int              ivalue;
  float            fvalue;
  double           dvalue;
  char            *svalue = NULL;


  if (header == NULL) {
    cpl_msg_error(fctid, "Missing header!");
    return status;
  }

  desc = findDescriptor(header, name);

  if (!desc) {
    cpl_msg_error(fctid, "Descriptor %s not found!", name);
    return status;
  }

  switch (desc->descType) {
  case VM_INT :
    ivalue = desc->descValue->i;
    break;
  case VM_FLOAT :
    fvalue = desc->descValue->f;
    break;
  case VM_DOUBLE :
    dvalue = desc->descValue->d;
    break;
  case VM_STRING :
    svalue = cpl_strdup(desc->descValue->s);
    if (!svalue) {
      cpl_msg_error(fctid, "Memory failure!");
      return status;
    }
    break;
  default :
    cpl_msg_error(fctid, "Unsupported descriptor type!");
    return status;
  }

  /*
   *  Construct entry name for PAF
   */

  keep = pafName = cpl_strdup(name);
  if (!pafName) {
    if (svalue)
      cpl_free(svalue);
    cpl_msg_error(fctid, "Memory failure!");
    return status;
  }

  pos = strstr(pafName, "ESO ");

  if (pos == pafName)
    pafName += 4;

  for (i = 0; pafName[i] != '\0'; i++)
    if (pafName[i] == ' ')
      pafName[i] = '.';

  /*
   *  Now write entry to PAF object.
   */

  switch (desc->descType) {
  case VM_INT :
    status = pilQcWriteInt(pafName, ivalue, unit, comment);
    break;
  case VM_FLOAT :
    dvalue = fvalue;
  case VM_DOUBLE :
    status = pilQcWriteDouble(pafName, dvalue, unit, comment);
    break;
  default :    /* VM_STRING */
    status = pilQcWriteString(pafName, svalue, comment);
  }

  if (status)
    cpl_msg_error(fctid, "Could not copy descriptor value to QC1 PAF!");

  if (svalue)
    cpl_free(svalue);
  cpl_free(keep);

  return status;

}


/**
 * @memo
 *   Write an integer value to the active QC1 PAF object and to a header.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param header  Pointer to a descriptor header.
 * @param value   Value to write.
 * @param name    QC1 PAF entry name.
 * @param unit    Optional unit to be associated to value.
 * @param comment Optional comment to be associated to value.
 *
 * @doc
 *   An entry with the specified name is written to the current QC1 PAF 
 *   object. From the entry name, the name of the QC descriptor that
 *   should be written to header is derived prepending the string "ESO "
 *   and replacing all '.' with a blank (e.g., "QC.BIAS.MASTER.MEAN"
 *   becomes "ESO QC BIAS MASTER MEAN"). Finally, the new descriptor
 *   is written to the header. Note that before calling this funtion 
 *   a QC1 PAF object must be created with a call to pilQcGroupStart().
 *
 * @author C. Izzo
 */

int
qcWriteValueInt(VimosDescriptor *header, int value, const char *name,
                const char *unit, const char *comment)
{

  const char fctid[] = "qcWriteValueInt";

  char            *descName;
  int              i, status;


  if (header == NULL) {
    cpl_msg_error(fctid, "Missing header!");
    return EXIT_FAILURE;
  }

  if (pilQcWriteInt(name, value, unit, comment) == EXIT_FAILURE) {
    cpl_msg_error(fctid, "Could not copy value to QC1 PAF!");
    return EXIT_FAILURE;
  }

  descName = cpl_malloc((strlen(name) + 5) * sizeof(char *));

  if (!descName) {
    cpl_msg_error(fctid, "Memory failure!");
    return EXIT_FAILURE;
  }

  strcpy(descName, "ESO ");
  strcat(descName, name);

  for (i = 0; descName[i] != '\0'; i++)
    if (descName[i] == '.')
      descName[i] = ' ';

  status = writeIntDescriptor(&header, descName, value, comment);

  cpl_free(descName);

  if (status == VM_FALSE) {
    cpl_msg_error(fctid, "Could not copy value to descriptor header!");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}


/**
 * @memo
 *   Write a double value to the active QC1 PAF object and to a header.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param header  Pointer to a descriptor header.
 * @param value   Value to write.
 * @param name    QC1 PAF entry name.
 * @param unit    Optional unit to be associated to value.
 * @param comment Optional comment to be associated to value.
 *
 * @doc
 *   An entry with the specified name is written to the current QC1 PAF 
 *   object. From the entry name, the name of the QC descriptor that
 *   should be written to header is derived prepending the string "ESO "
 *   and replacing all '.' with a blank (e.g., "QC.BIAS.MASTER.MEAN"
 *   becomes "ESO QC BIAS MASTER MEAN"). Finally, the new descriptor
 *   is written to the header. Note that before calling this funtion 
 *   a QC1 PAF object must be created with a call to pilQcGroupStart().
 *
 * @author C. Izzo
 */

int
qcWriteValueDouble(VimosDescriptor *header, double value, const char *name,
                   const char *unit, const char *comment)
{

  const char fctid[] = "qcWriteValueDouble";

  char            *descName;
  int              i, status;


  if (header == NULL) {
    cpl_msg_error(fctid, "Missing header!");
    return EXIT_FAILURE;
  }

  if (pilQcWriteDouble(name, value, unit, comment) == EXIT_FAILURE) {
    cpl_msg_error(fctid, "Could not copy value to QC1 PAF!");
    return EXIT_FAILURE;
  }

  descName = cpl_malloc((strlen(name) + 5) * sizeof(char *));

  if (!descName) {
    cpl_msg_error(fctid, "Memory failure!");
    return EXIT_FAILURE;
  }

  strcpy(descName, "ESO ");
  strcat(descName, name);

  for (i = 0; descName[i] != '\0'; i++)
    if (descName[i] == '.')
      descName[i] = ' ';

  status = writeDoubleDescriptor(&header, descName, value, comment);

  cpl_free(descName);

  if (status == VM_FALSE) {
    cpl_msg_error(fctid, "Could not copy value to descriptor header!");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}


/**
 * @memo
 *   Write an integer value to the active QC1 PAF object and to a header.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param filnam  Name of existing FITS file.
 * @param value   Value to write.
 * @param name    QC1 PAF entry name.
 * @param unit    Optional unit to be associated to value.
 * @param comment Optional comment to be associated to value.
 *
 * @doc
 *   This is just identical to the function qcWriteValueInt(), but
 *   it writes the header entries directly to the header of the FITS
 *   file written to disk, using the qfits_replace_card() call.
 *
 *   An entry with the specified name is written to the current QC1 PAF 
 *   object. From the entry name, the name of the QC descriptor that
 *   should be written to header is derived prepending the string "ESO "
 *   and replacing all '.' with a blank (e.g., "QC.BIAS.MASTER.MEAN"
 *   becomes "ESO QC BIAS MASTER MEAN"). Finally, the new descriptor
 *   is written to the header. Note that before calling this funtion 
 *   a QC1 PAF object must be created with a call to pilQcGroupStart().
 *
 * @author C. Izzo
 */

int
qcWriteValueInt_CPL(char *filnam, int value, const char *name,
                const char *unit, const char *comment)
{

  const char fctid[] = "qcWriteValueInt_CPL";

/*
  char             line[80];
  char             val[80];
*/
  char            *descName;
  fitsfile        *fptr;
  int              status = 0;
  int              i;


  if (pilQcWriteInt(name, value, unit, comment) == EXIT_FAILURE) {
    cpl_msg_error(fctid, "Could not copy value to QC1 PAF!");
    return EXIT_FAILURE;
  }

  descName = cpl_malloc((strlen(name) + 15) * sizeof(char *));

  if (!descName) {
    cpl_msg_error(fctid, "Memory failure!");
    return EXIT_FAILURE;
  }

  strcpy(descName, "HIERARCH ESO ");
  strcat(descName, name);

  for (i = 0; descName[i] != '\0'; i++)
    if (descName[i] == '.')
      descName[i] = ' ';

  ffopen(&fptr, filnam, READWRITE, &status);
  fits_update_key(fptr, TINT, descName, &value, (char *)comment, &status);
  fits_close_file(fptr, &status);

/*
  sprintf(val, "%d", value);
  keytuple2str(line, descName, val, (char *)comment);
  qfits_replace_card(filnam, descName, line);
*/

  cpl_free(descName);

  if (status)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @memo
 *   Write a double value to the active QC1 PAF object and to a header.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param filnam  Name of existing FITS file.
 * @param value   Value to write.
 * @param name    QC1 PAF entry name.
 * @param unit    Optional unit to be associated to value.
 * @param comment Optional comment to be associated to value.
 *
 * @doc
 *   This is just identical to the function qcWriteValueDouble(), but
 *   it writes the header entries directly to the header of the FITS
 *   file written to disk, using the qfits_replace_card() call.
 *
 *   An entry with the specified name is written to the current QC1 PAF 
 *   object. From the entry name, the name of the QC descriptor that
 *   should be written to header is derived prepending the string "ESO "
 *   and replacing all '.' with a blank (e.g., "QC.BIAS.MASTER.MEAN"
 *   becomes "ESO QC BIAS MASTER MEAN"). Finally, the new descriptor
 *   is written to the header. Note that before calling this funtion 
 *   a QC1 PAF object must be created with a call to pilQcGroupStart().
 *
 * @author C. Izzo
 */

int
qcWriteValueDouble_CPL(char *filnam, double value, const char *name,
                   const char *unit, const char *comment)
{

  const char fctid[] = "qcWriteValueDouble_CPL";

/*
  char             line[80];
  char             val[80];
*/
  char            *descName;
  fitsfile        *fptr;
  int              status = 0;
  int              i;


  if (pilQcWriteDouble(name, value, unit, comment) == EXIT_FAILURE) {
    cpl_msg_error(fctid, "Could not copy value to QC1 PAF!");
    return EXIT_FAILURE;
  }

  descName = cpl_malloc((strlen(name) + 15) * sizeof(char *));

  if (!descName) {
    cpl_msg_error(fctid, "Memory failure!");
    return EXIT_FAILURE;
  }

  strcpy(descName, "HIERARCH ESO ");
  strcat(descName, name);

  for (i = 0; descName[i] != '\0'; i++)
    if (descName[i] == '.')
      descName[i] = ' ';

  ffopen(&fptr, filnam, READWRITE, &status);
  fits_update_key(fptr, TDOUBLE, descName, &value, (char *)comment, &status);
  fits_close_file(fptr, &status);

/*
  sprintf(val, "%1.6e", value);
  keytuple2str(line, descName, val, (char *)comment);
  qfits_replace_card(filnam, descName, line);
*/

  cpl_free(descName);

  if (status)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}
/**@}*/
