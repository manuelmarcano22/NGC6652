/* $Id: vmutils.c,v 1.17 2013-08-23 10:23:02 cgarcia Exp $
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
 * $Date: 2013-08-23 10:23:02 $
 * $Revision: 1.17 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <fitsio.h>
#include <fitsio2.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <cxstring.h>
#include <cpl_propertylist.h>
#include <cpl_version.h>

#include <pildfsconfig.h>
#include <pilframeset.h>
#include <pilmemory.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <pilrecipe.h>
#include <pilfits.h>
#include <cpl_msg.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmutils.h"
#include "cpl.h"


/* Set parameters for PAF writer */

#define PAF_VALUE_POSITION (30)


/**
 * @name vmutils  Miscellaneous Utilitiy Functions
 *
 * @doc
 *   The module collects high/medium level functions
 *   for common operations in the Imaging DRS.
 */

/**@{*/

static char *_get_base_name(const char *filename)
{
    char *p ;
    p = strrchr (filename, '/');
    return p ? p + 1 : (char *) filename;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the pipeline copyright and license
  @return   The copyright and license string

  The function returns a pointer to the statically allocated license string.
  This string should not be modified using the returned pointer.
 */
/*----------------------------------------------------------------------------*/
const char * vimos_get_license(void)
{
    const char  *   vimos_license = 
        "This file is currently part of the VIMOS Instrument Pipeline\n"
        "Copyright (C) 2002-2011 European Southern Observatory\n\n"
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
        "Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA\n";
    return vimos_license;
}

/**
 * @memo
 *   Write a double precision keyword to a PAF file.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param fp         Pointer to PAF file
 * @param name       Keyword name
 * @param value      Keyword value
 *
 * @doc
 *   The keyword value is written at the position specified
 *   by the parameter PAF_VALUE_POSITION, with the same format
 *   used in FITS file double keywords to avoid losing precision.
 *
 * @author C. Izzo
 */

int
writeDoublePAFEntry(FILE *fp, char *name, double value)
{
  const char  modName[] = "writeDoublePAFEntry";
  int         nBlanks;

  if (name) {
    nBlanks = PAF_VALUE_POSITION - strlen(name);
    if (nBlanks < 1) nBlanks = 1;
  }
  else {
    cpl_msg_debug(modName, "Undefined parameter name");
    return EXIT_FAILURE;
  }
  
  fprintf(fp, "%s%*s\"%.14E\";\n", name, nBlanks, " ", value);

  return EXIT_SUCCESS;
}


/**
 * @memo
 *   Write an integer keyword to a PAF file.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param fp         Pointer to PAF file
 * @param name       Keyword name
 * @param value      Keyword value
 *
 * @doc
 *   The keyword value is written at the position specified
 *   by the parameter PAF_VALUE_POSITION.
 *
 * @author C. Izzo
 */

int
writeIntPAFEntry(FILE *fp, char *name, int value)
{
  const char  modName[] = "writeIntPAFEntry";
  int         nBlanks;

  if (name) {
    nBlanks = PAF_VALUE_POSITION - strlen(name);
    if (nBlanks < 1) nBlanks = 1;
  }
  else {
    cpl_msg_debug(modName, "Undefined parameter name");
    return EXIT_FAILURE;
  }
  
  fprintf(fp, "%s%*s\"%d\";\n", name, nBlanks, " ", value);
  return EXIT_SUCCESS;
}

/**
 * @memo
 *   Write a character string keyword to a PAF file.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE
 *
 * @param fp         Pointer to PAF file
 * @param name       Keyword name
 * @param value      Keyword value
 *
 * @doc
 *   The keyword value is written at the position specified
 *   by the parameter PAF_VALUE_POSITION.
 *
 * @author C. Izzo
 */

int
writeStringPAFEntry(FILE *fp, char *name, char *value)
{
  const char  modName[] = "writeStringPAFEntry";
  int         nBlanks;

  if (name) {
    nBlanks = PAF_VALUE_POSITION - strlen(name);
    if (nBlanks < 1) nBlanks = 1;
  }
  else {
    cpl_msg_debug(modName, "Undefined parameter name");
    return EXIT_FAILURE;
  }

  if (value) fprintf(fp, "%s%*s\"%s\";\n", name, nBlanks, " ", value);
  else       fprintf(fp, "%s;\n", name);

  return EXIT_SUCCESS;
}


/**
 * @memo 
 *   Update descriptor list of a recipe image product with standard keywords
 *
 * @return EXIT_SUCCESS/EXIT_FAILURE
 * 
 * @param image     input image 
 * @param category  category of product (done at this level for the moment)
 *
 * @doc
 *   Write in the descriptor list of a recipe image product the 
 *   standard keywords required by the pipeline 
 *   (ex: data minimum/maximum/mean/median value..) and remove keywords 
 *   like DPR, TEL keywords.   
 *
 * @author P. Sartoretti
 */    

int
UpdateProductDescriptors(VimosImage *image, const char *category)

{
  char modName[] = "UpdateProductDescriptors";
  VimosDescriptor  *minDesc, *maxDesc;


  if ((image == NULL) || (image->descs == NULL)) {
    cpl_msg_error(modName, "No descriptor list to update");
    return EXIT_FAILURE;
  }

  /*
   * DATAMIN and DATAMAX are standard fits keys and are before
   * the hierarch keys
   */

  minDesc = newDoubleDescriptor(pilTrnGetKeyword("DataMin"),
                       imageMinimum(image), "Minimum pixel value");
  if (insertDescriptor(&(image->descs), "ESO*", minDesc, 1)==VM_FALSE) {
    cpl_msg_warning(modName, "Cannot write descriptor %s",
                  pilTrnGetKeyword("DataMin"));
  }
  
  maxDesc = newDoubleDescriptor(pilTrnGetKeyword("DataMax"),
                       imageMaximum(image), "Maximum pixel value");
  if (insertDescriptor(&(image->descs), "ESO*", maxDesc, 1)==VM_FALSE) {
    cpl_msg_warning(modName, "Cannot write descriptor %s",
                  pilTrnGetKeyword("DataMax"));
  }

  if ((writeDoubleDescriptor(&(image->descs), pilTrnGetKeyword("DataMean"),
                     imageMean(image), "Mean pixel value"))==VM_FALSE) {
    cpl_msg_warning(modName, "Cannot write descriptor %s",
                  pilTrnGetKeyword("DataMean"));
  }
  if ((writeDoubleDescriptor(&(image->descs), 
                            pilTrnGetKeyword("DataStdDeviation"),
                            imageSigma(image), 
                             "Standard deviation of pixel"))==VM_FALSE) {
    cpl_msg_warning(modName, "Cannot write descriptor %s",
                  pilTrnGetKeyword("DataStdDeviation"));
  }

  if ((writeDoubleDescriptor(&(image->descs), pilTrnGetKeyword("DataMedian"),
                             imageMedian(image), "Median pixel value"))
      ==VM_FALSE) {
    cpl_msg_warning(modName, "Cannot write descriptor %s",
                  pilTrnGetKeyword("DataMedian"));
  }

  if ((writeStringDescriptor(&(image->descs), pilTrnGetKeyword("DoCategory"),
                      category, "Category of pipeline product")==VM_FALSE)) {
    cpl_msg_warning(modName, "Cannot write descriptor %s",
                  pilTrnGetKeyword("DoCategory"));
  }
  
  if ((deleteSetOfDescriptors(&(image->descs), "ESO DPR*")) == 0) {
    cpl_msg_warning(modName, "Cannot remove descriptors ESO DPR*");
  }
  return EXIT_SUCCESS;
}


/**
 * @memo
 *   Apply a selection criterium to a list of images.
 *
 * @return The number of images in the list matching the given criterium.
 *   In case of an error the return value is -1;
 *
 * @param imageList       List of images.
 * @param imageProperty   Array of characteristic values.
 * @param imageCount      Number of images in the list.
 * @param rangeLow        Lower limit for the tolerance range.
 * @param rangeHigh       Higher limit for the tolerance range.
 * @param accept          Flag indicating the acceptance mode.
 *
 * @doc
 *   The function removes all images from \textbf{imageList} whose
 *   characteristic values are not within the given range. The validity 
 *   range is build from the reference value \textbf{referenceValue}
 *   and the two tolarance values \textbf{daltaLow} and \textbf{deltaHigh}.
 *   To be selected from \textbf{list} an image must fulfill one of the
 *   following conditions:
 *   \[
 *     imageProperty_{i} \in [rangeLow, rangeHigh]
 *   \]
 *   or
 *   \[
 *     imageProperty_{i} \not\in [rangeLow, rangeHigh]
 *   \]
 *   Which relation is applied depends on the setting of the acceptance
 *   flag \textbf{accept}. If \textbf{accept} is set to zero the second
 *   condition is be used, otherwise the first condition is used.
 *
 *   The order of images in the input list is modified so that the list
 *   contains the selected frames first followed by the rejected frames.
 *   The same modification is applied to the array of characteristic values
 *   so that in the end each entry in the array \textbf{imageProperty} still
 *   corresponds to the same image entry in \textbf{imageList}.
 *
 *   FIXME: Instead of moving images within the list a new list object
 *   holding the accepted frames should be created. Requires an image
 *   list object!
 *
 * @author R. Palsa
 */

int
applyListSelection(VimosImage **imageList, float *imageProperty,
                   int imageCount, double rangeLow, double rangeHigh, 
                   unsigned int accept)
{

  const char fctId[] = "applyListSelection";


  register int i, j = 0, k = 0;

  int isValid;

  float *rejectedProperty;

  VimosImage **rejectedList;



  /*
   * Allocate temporary list to save the rejected images and their
   * property value.
   */

  rejectedList = (VimosImage **)cpl_malloc(imageCount * sizeof(VimosImage *));
  if (!rejectedList)
    return -1;

  rejectedProperty = (float *)cpl_malloc(imageCount * sizeof(float));
  if (!rejectedProperty) {
    cpl_free(rejectedList);
    return -1;
  }


  /*
   * Check the frames in the list.
   */

  for (i = 0; i < imageCount; i++) {
    if (!accept)
      isValid = imageProperty[i] < rangeLow || imageProperty[i] > rangeHigh;
    else 
      isValid = rangeLow <= imageProperty[i] && imageProperty[i] <= rangeHigh;

    if (isValid) {
      if (i > j) {
        imageList[j] = imageList[i];
        imageProperty[j] = imageProperty[i];
      }
      j++;
    }
    else {
      cpl_msg_debug(fctId, "Image %d removed from list.", i + 1);
      rejectedList[k] = imageList[i];
      rejectedProperty[k] = imageProperty[i];
      k++;
    }
  }


  /*
   * Insert rejected images back into the original list.
   */

  for (i = j, k = 0; i < imageCount; i++, k++) {
    imageList[i] = rejectedList[k];
    imageProperty[i] = rejectedProperty[k];
  }


  /*
   * Cleanup
   */

  cpl_free(rejectedList);
  cpl_free(rejectedProperty);

  return j;

}


/**
 * @memo
 *   Remap float array according to a resorting of an image list.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE;
 *
 * @param original   Array of VimosImage pointers, as it was before sorting
 * @param sorted     Array of VimosImage pointers, as it is after sorting
 * @param array      Array of floats, to be remapped according to image 
 *                   list remapping
 * @param imageCount Number of input images
 *
 * @doc
 *   Functions like applyListSelection() and qcSelectConsistentImages()
 *   perform a resorting of an input image list, moving the good frames
 *   at the beginning of the list. In case other lists of values were
 *   associated to the original image list, they would no longer match
 *   the corresponding images. To avoid this problem, such lists should 
 *   be remapped in the same way as the rearranged image list. This function
 *   does exactly this: the images array befor and after sorting is passed,
 *   together with the float array to be reordered. Remapping of the float
 *   values is done in place.
 *
 * @see applyListSelection, qcSelectConsistentImages
 *
 * @author C. Izzo
 */

int
remapFloatsLikeImages(VimosImage **original, VimosImage **sorted, 
                      float *array, int imageCount)
{
  const char modName[] = "remapFloatsLikeImages";

  int    i, j;
  float *holder;
  int   *error;

  if (original == NULL) {
    cpl_msg_debug(modName, "NULL input array of images");
    return EXIT_FAILURE;
  }

  if (sorted == NULL) {
    cpl_msg_debug(modName, "NULL input array of images");
    return EXIT_FAILURE;
  }

  if (array == NULL) {
    cpl_msg_debug(modName, "NULL input array of floats");
    return EXIT_FAILURE;
  }

  if (imageCount < 1) {
    cpl_msg_debug(modName, "Wrong number of input images (%d)", imageCount);
    return EXIT_FAILURE;
  }

  if (imageCount == 1)
    return EXIT_SUCCESS;               /* No need for remapping */

  for (i = 0; i < imageCount; i++) {
    if (original[i] == NULL || sorted[i] == NULL) {
      cpl_msg_debug(modName, "NULL images in input");
      return EXIT_FAILURE;
    }
  }

  if (!(holder = (float *)cpl_malloc(imageCount*sizeof(float)))) {
    cpl_msg_debug(modName, "Problems with memory allocation");
    return EXIT_FAILURE;
  }

  if (!(error = (int *)cpl_malloc(imageCount*sizeof(int)))) {
    cpl_free(holder);
    cpl_msg_debug(modName, "Problems with memory allocation");
    return EXIT_FAILURE;
  }

  for (i = 0; i < imageCount; i++)
    error[i] = 1;               /* To catch any mismatch or malfunctioning */

  for (i = 0; i < imageCount;  i++) {  /* Loop on sorted list */
    for (j = 0; j < imageCount; j++) {  /* Loop on unsorted (original) list */
      if (original[j] == sorted[i]) {    /* Found old position j for i */
        holder[i] = array[j];             /* Reassign j-th value of array */
        error[i] = 0;
        break;
      }
    }
  }

  for (i = 0; i < imageCount;  i++) {   /* Error checking */
    if (error[i]) {
      cpl_free(holder);
      cpl_free(error);
      cpl_msg_debug(modName, "Input image arrays are not comparable");
      return EXIT_FAILURE;
    }
  }

  for (i = 0; i < imageCount;  i++)  /* New sorting overwrites old one */
    array[i] = holder[i];

  cpl_free(holder);
  cpl_free(error);
  
  return EXIT_SUCCESS;

}


/**
 * @memo
 *   Remap double array according to a resorting of an image list.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE;
 *
 * @param original   Array of VimosImage pointers, as it was before sorting
 * @param sorted     Array of VimosImage pointers, as it is after sorting
 * @param array      Array of doubles, to be remapped according to image 
 *                   list remapping
 * @param imageCount Number of input images
 *
 * @doc
 *   Functions like applyListSelection() and qcSelectConsistentImages()
 *   perform a resorting of an input image list, moving the good frames
 *   at the beginning of the list. In case other lists of values were
 *   associated to the original image list, they would no longer match
 *   the corresponding images. To avoid this problem, such lists should 
 *   be remapped in the same way as the rearranged image list. This function
 *   does exactly this: the images array befor and after sorting is passed,
 *   together with the double array to be reordered. Remapping of the double
 *   values is done in place.
 *
 * @see applyListSelection, qcSelectConsistentImages
 *
 * @author C. Izzo
 */

int
remapDoublesLikeImages(VimosImage **original, VimosImage **sorted, 
                       double *array, int imageCount)
{
  const char modName[] = "remapDoublesLikeImages";

  int     i, j;
  double *holder;
  int    *error;

  if (original == NULL) {
    cpl_msg_debug(modName, "NULL input array of images");
    return EXIT_FAILURE;
  }

  if (sorted == NULL) {
    cpl_msg_debug(modName, "NULL input array of images");
    return EXIT_FAILURE;
  }

  if (array == NULL) {
    cpl_msg_debug(modName, "NULL input array of doubles");
    return EXIT_FAILURE;
  }

  if (imageCount < 1) {
    cpl_msg_debug(modName, "Wrong number of input images (%d)", imageCount);
    return EXIT_FAILURE;
  }

  if (imageCount == 1)
    return EXIT_SUCCESS;               /* No need for remapping */

  for (i = 0; i < imageCount; i++) {
    if (original[i] == NULL || sorted[i] == NULL) {
      cpl_msg_debug(modName, "NULL images in input");
      return EXIT_FAILURE;
    }
  }

  if (!(holder = (double *)cpl_malloc(imageCount*sizeof(double)))) {
    cpl_msg_debug(modName, "Problems with memory allocation");
    return EXIT_FAILURE;
  }

  if (!(error = (int *)cpl_malloc(imageCount*sizeof(int)))) {
    cpl_free(holder);
    cpl_msg_debug(modName, "Problems with memory allocation");
    return EXIT_FAILURE;
  }

  for (i = 0; i < imageCount; i++)
    error[i] = 1;               /* To catch any mismatch or malfunctioning */

  for (i = 0; i < imageCount;  i++) {  /* Loop on sorted list */
    for (j = 0; j < imageCount; j++) {  /* Loop on unsorted (original) list */
      if (original[j] == sorted[i]) {    /* Found old position j for i */
        holder[i] = array[j];             /* Reassign j-th value of array */
        error[i] = 0;
        break;
      }
    }
  }

  for (i = 0; i < imageCount;  i++) {   /* Error checking */
    if (error[i]) {
      cpl_free(holder);
      cpl_free(error);
      cpl_msg_debug(modName, "Input image arrays are not comparable");
      return EXIT_FAILURE;
    }
  }

  for (i = 0; i < imageCount;  i++)  /* New sorting overwrites old one */
    array[i] = holder[i];

  cpl_free(holder);
  cpl_free(error);
  
  return EXIT_SUCCESS;

}


void sortN (int ncol, float **ra, int sortCol, int fromRow, int forRows)
     /* liberally adapted from NR (BG)*/
     /* sort  matrix ra of ncol columns
        for coloumn sortCol in ascending order, and rearranges all other
        columns sorting is done only over the portion fromRow-forRows */
{
        int i,j,i1,*iwksp;
        float *wksp;

        iwksp=cpl_calloc(forRows, sizeof(int *));
        wksp=cpl_calloc(forRows, sizeof(float *));
        j=0;
        for (i=fromRow;i<fromRow+forRows;i++){
          wksp[j]=ra[i][sortCol];
          j++;
        }
        Indexx(forRows,wksp,iwksp);
        for (j=0; j<ncol; j++) {
          i1=0;
          for (i=fromRow;i<fromRow+forRows;i++) {
            wksp[i1]=ra[i][j];
            i1++;
          }
          i1=0;
          for (i=fromRow;i<fromRow+forRows;i++) {
            ra[i][j]=wksp[iwksp[i1]];
            i1++;
          }
        }
        cpl_free(wksp);
        cpl_free(iwksp);
}


inline static cxint
_vm_plist_append(cpl_propertylist *self, cpl_property *p)
{

    const cxchar *name = cpl_property_get_name(p);
    const cxchar *comment = cpl_property_get_comment(p);


    switch (cpl_property_get_type(p)) {
    case CPL_TYPE_BOOL:
    {
        cxbool value = cpl_property_get_bool(p);

        cpl_propertylist_append_bool(self, name, value);
        break;
    }
                
    case CPL_TYPE_CHAR:
    {
        cxchar value = cpl_property_get_char(p);

        cpl_propertylist_append_char(self, name, value);
        break;
    }

    case CPL_TYPE_INT:
    {
        cxint value = cpl_property_get_int(p);

        cpl_propertylist_append_int(self, name, value);
        break;
    }

    case CPL_TYPE_LONG:
    {
        cxlong value = cpl_property_get_long(p);

        cpl_propertylist_append_long(self, name, value);
        break;
    }

    case CPL_TYPE_FLOAT:
    {
        cxfloat value = cpl_property_get_float(p);

        cpl_propertylist_append_float(self, name, value);
        break;
    }

    case CPL_TYPE_DOUBLE:
    {
        cxdouble value = cpl_property_get_double(p);

        cpl_propertylist_append_double(self, name, value);
        break;
    }

    case CPL_TYPE_STRING:
    {
        const cxchar *value = cpl_property_get_string(p);

        cpl_propertylist_append_string(self, name, value);
        break;
    }

    default:

        /*
         * We should never reach this point! Since the property
         * was a valid property it has a valid type. But this
         * protects against addition of new types.
         */

        return 1;
        break;
    }

    if (comment != NULL) {
        cpl_propertylist_set_comment(self, name, comment);
    }

    return 0;

}

/**
 * @brief
 *   Update a property list.
 *
 * @param self        The property list to update.
 * @param properties  The source property list.
 * @param regexp      A property name pattern.
 *
 * @return The function returns 0 on success and a non-zero value in case an
 *   error occurred.
 *
 * The function updates the target property list with properties from the
 * source list @em properties, which are not present in @em self. If a pattern
 * string is given only properties with names matching the given pattern
 * @em regexp are taken into account when @em self is updated. If a pattern
 * is given, it must be a valid regular expression. If the pattern string is
 * either @c NULL or the empty string, the whole source list is considered
 * during the update operation.
 */

cxint
vm_plist_update(cpl_propertylist *self, cpl_propertylist *properties,
                     const cxchar *regexp)
{

    cxlong i;
    cxlong sz = 0;


    cx_assert(self != NULL);

    if (properties == NULL) {
        return -1;
    }

    sz = cpl_propertylist_get_size(properties);

    if (regexp == NULL || regexp[0] == '\0') {

        for (i = 0; i < sz; ++i) {

            cpl_property *p = cpl_propertylist_get(properties, i);
            const cxchar *name = cpl_property_get_name(p);


            if (!cpl_propertylist_has(self, name)) {

                cxint status = _vm_plist_append(self, p);

                if (status) {
                    return 1;
                }
            }
        }
    }
    else {

        cxint status = 0;

        regex_t re;


        status = regcomp(&re, regexp, REG_EXTENDED | REG_NOSUB);

        if (status) {
            return 1;
        }

        for (i = 0; i < sz; ++i) {

            cpl_property *p = cpl_propertylist_get(properties, i);
            const cxchar *name = cpl_property_get_name(p);

            if (regexec(&re, name, (size_t)0, NULL, 0) == REG_NOMATCH) {
                continue;
            }

            if (!cpl_propertylist_has(self, name)) {

                cxint status = _vm_plist_append(self, p);

                if (status) {
                    return 1;
                }
            }
        }

        regfree(&re);

    }

    return 0;

}

/*
 * Product keywords aliases
 */

#define DATAMD5                     "DATAMD5"
#define PIPEFILE                    "PIPEFILE"
#define PRO_DID                     "ESO PRO DID"
#define DPR_CATG                    "ESO DPR CATG"
#define PRO_CATG                    "ESO PRO CATG"
#define PRO_TYPE                    "ESO PRO TYPE"
#define DPR_TECH                    "ESO DPR TECH"
#define PRO_TECH                    "ESO PRO TECH"
#define PRO_SCIENCE                 "ESO PRO SCIENCE"
#define PRO_DATE                    "ESO PRO DATE"
#define PRO_DATANCOM                "ESO PRO DATANCOM"
#define PRO_REC_ID                  "ESO PRO REC1 ID"
#define PRO_REC_DRS_ID              "ESO PRO REC1 DRS ID"
#define PRO_REC_PIPE_ID             "ESO PRO REC1 PIPE ID"
#define PRO_REC_RAWi_NAME           "ESO PRO REC1 RAW%d NAME"
#define PRO_REC_RAWi_CATG           "ESO PRO REC1 RAW%d CATG"
#define PRO_REC_CALi_NAME           "ESO PRO REC1 CAL%d NAME"
#define PRO_REC_CALi_CATG           "ESO PRO REC1 CAL%d CATG"
#define PRO_REC_CALi_DATAMD5        "ESO PRO REC1 CAL%d DATAMD5"
#define PRO_REC_PARAMi_NAME         "ESO PRO REC1 PARAM%d NAME"
#define PRO_REC_PARAMi_VALUE        "ESO PRO REC1 PARAM%d VALUE"

#define MAX_PLENGTH (64)

/* Little function to insert a new record in a FITS header:
 * If the input keyword is present in header, then
 *     overwrite the corresponding card;
 * Else, if the reference keyword exists, 
 *     insert the input keyword card as indicated;
 * Else
 *     append the input keyword card at the end of header.
 */

/*************

static int _try_insert_key(fitsfile *fptr, const char *keyName, int before, 
                           const char *refKeyName, const char *icard)
{

    int   status = 0;
    int   countKeys = 0;
    int   position = 0;
    int   keyExists = 1;
    int   refKeyExists = 1;
    char  card[FLEN_CARD];

    
    if (!fptr)
        return EXIT_FAILURE;

    fits_read_record(fptr, 0, card, &status);
    fits_find_nextkey(fptr, (char **)&keyName, 1, NULL, 0, card, &status);

    if (status) {
        if (status == KEY_NO_EXIST) {
            keyExists = 0;
            status = 0;
        }
        else
            return status;
    }

    fits_read_record(fptr, 0, card, &status);
    fits_find_nextkey(fptr, (char **)&refKeyName, 1, NULL, 0, card, &status);

    if (status) {
        if (status == KEY_NO_EXIST) {
            refKeyExists = 0;
            status = 0;
        }
        else
            return status;
    }

    if (refKeyExists && !(keyExists)) {
        fits_get_hdrpos(fptr, &countKeys, &position, &status);

        if (before == 0)
            position++;

        fits_insert_record(fptr, position, (char *)icard, &status);
    }
    else
        fits_update_card(fptr, (char *)keyName, (char *)icard, &status);

    return status;

}

*****************/

/* This one counts RAW frames starting from 0, and if seq >= the total
   number of RAW frames in framelist it returns NULL */

static PilFrame *get_raw_frame_time_sequence(PilSetOfFrames *framelist, int seq)
{

    static char *fid = "get_raw_frame_time_sequence";

    PilFrame  *frame = pilSofFirst(framelist);
    int        nframes = 0;
    double    *times;
    PilFrame **sort_frames;

    double       tobs;
    int          i, j, pos;

    char         errText[50];
    fitsfile    *fptr;
    int          status = 0;


    while (frame != NULL) {
        if (pilFrmGetType(frame) == PIL_FRAME_TYPE_RAW
            || pilFrmGetType(frame) == PIL_FRAME_TYPE_UNDEF) {
            if (ffopen(&fptr, pilFrmGetName(frame), READONLY, &status))
            {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Opening file %s: %s",
                              pilFrmGetName(frame), errText);
                return NULL;
            }
            if (fits_read_key(fptr, TDOUBLE, "MJD-OBS", &tobs, NULL, &status)) {
                status = 0;
                fits_close_file(fptr, &status);
                frame = pilSofNext(framelist, frame);
                continue;
            }
            if (fits_close_file(fptr, &status))             {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Closing file %s: %s",
                              pilFrmGetName(frame), errText);
                return NULL;
            }
            nframes++;
        }
        frame = pilSofNext(framelist, frame);
    }

    if (seq > nframes)
        return NULL;

    times = calloc(nframes, sizeof(double));
    sort_frames = calloc(nframes, sizeof(PilFrame *));
    nframes = 0;
    frame = pilSofFirst(framelist);

    while (frame != NULL) {
        if (pilFrmGetType(frame) == PIL_FRAME_TYPE_RAW
            || pilFrmGetType(frame) == PIL_FRAME_TYPE_UNDEF) {
            if (ffopen(&fptr, pilFrmGetName(frame), READONLY, &status))
            {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Opening file %s: %s",
                              pilFrmGetName(frame), errText);
                free(times);
                free(sort_frames);
                return NULL;
            }
            if (fits_read_key(fptr, TDOUBLE, "MJD-OBS", &tobs, NULL, &status)) {
                status = 0;
                fits_close_file(fptr, &status);
                frame = pilSofNext(framelist, frame);
                continue;
            }
            if (fits_close_file(fptr, &status))             {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Closing file %s: %s",
                              pilFrmGetName(frame), errText);
                free(times);
                free(sort_frames);
                return NULL;
            }
            times[nframes] = tobs;
            sort_frames[nframes] = frame;
            nframes++;
        }
        frame = pilSofNext(framelist, frame);
    }

    for (i = 0; i < nframes; i++) {
        pos = 0;
        for (j = 0; j < nframes; j++) {
            if (j == i)
                continue;
            if (times[i] > times[j])
                pos++;
        }
        if (pos == seq) {
            frame = sort_frames[i];
            free(times);
            free(sort_frames);
            return frame;
        }
    }

    free(times);
    free(sort_frames);
    return NULL;
}


int vm_dfs_setup_product_header(PilFrame *proFrame, const char *recipename, 
                                PilSetOfFrames *framelist)
{

    static char *fid = "vm_dfs_setup_product_header";
    char         cval[FLEN_CARD];
    PilFrame    *frame;
    PilFrame    *first_frame = NULL;
    int          nraw;
    int          ncal;
    int          npar;
    int          size;
    int          i;

    char        *productName;
    char         errText[50];
    char         card[FLEN_CARD];
    char         value[FLEN_CARD];
    char         kname[FLEN_CARD];
    char         parName[FLEN_CARD];
    char         parValue[MAX_PLENGTH + 1];
    char       **parString;
    fitsfile    *ifptr;
    fitsfile    *ofptr;
    fitsfile    *fptr;
    int          status = 0;
    char        *md5;
    int          pro_science;


    /*
     * Here is the list of mandatory keywords, first in input header,
     * second in output _image_ header:
     */

    typedef struct {
        char        *key;
        VimosVarType type;
    } klist;

    klist mandatory[] =  {{"ORIGIN", VM_STRING},
                          {"TELESCOP", VM_STRING},
                          {"INSTRUME", VM_STRING},
                          {"OBJECT", VM_STRING},
                          {"RA", VM_DOUBLE},
                          {"DEC", VM_DOUBLE},
                          {"EPOCH", VM_STRING},
                          {"EQUINOX", VM_DOUBLE},
                          {"RADECSYS", VM_STRING},
                          {"DATE-OBS", VM_STRING},
                          {"MJD-OBS", VM_DOUBLE},
                          {"UTC", VM_DOUBLE},
                          {"LST", VM_DOUBLE},
                          {"PI-COI", VM_STRING},
                          {"OBSERVER", VM_STRING}};

    long count_mandatory = sizeof(mandatory) / sizeof(klist);


    if (recipename == 0x0) {
        cpl_msg_error(fid, "Missing recipe name");
        return EXIT_FAILURE;
    }

    if (proFrame == 0x0 || framelist == 0x0) {
        cpl_msg_error(fid, "Null input");
        return EXIT_FAILURE;
    }

    productName = (char *)pilFrmGetName(proFrame);


    md5 = pilFitsMD5Signature(productName);
    if (!md5) {
        cpl_msg_error(fid, "Cannot compute MD5 signature for file %s",
                      productName);
        return EXIT_FAILURE;
    }


    /*
     * Get the first input frame in the input set-of-frames.
     */

    first_frame = get_raw_frame_time_sequence(framelist, 0);

/********* OLD part...

    mintobs = 100000.0;
    frame = pilSofFirst(framelist);

    while (frame != NULL) {
        if (pilFrmGetType(frame) == PIL_FRAME_TYPE_RAW 
            || pilFrmGetType(frame) == PIL_FRAME_TYPE_UNDEF) {
            if (ffopen(&fptr, pilFrmGetName(frame), READONLY, &status))
            {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Opening file %s: %s",
                              pilFrmGetName(frame), errText);
                return EXIT_FAILURE;
            }
            if (fits_read_key(fptr, TDOUBLE, "MJD-OBS", &tobs, NULL, &status)) {
                status = 0;
                fits_close_file(fptr, &status);
                frame = pilSofNext(framelist, frame);
                continue;
            }
            if (fits_close_file(fptr, &status))             {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Closing file %s: %s",
                              pilFrmGetName(frame), errText);
                return EXIT_FAILURE;
            }
            if (tobs < mintobs) {
                mintobs = tobs;
                first_frame = frame;
            }
        }
        frame = pilSofNext(framelist, frame);
    }

****************** End of old part */

    if (first_frame == NULL) {
        frame = pilSofFirst(framelist);
    
        while (frame != NULL) {
            if (pilFrmGetType(frame) == PIL_FRAME_TYPE_CALIB) {
                first_frame = frame;
                break;
            }
            frame = pilSofNext(framelist, frame);
        }
    }

    if (first_frame == NULL) {
        cpl_msg_error(fid, "Data not found in SOF");
        return EXIT_FAILURE;
    }

    /*
     * Now copy all the required entries, if present, from input to 
     * product header.
     */

    if (ffopen(&ifptr, pilFrmGetName(first_frame), READONLY, &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Opening file %s: %s", 
                      pilFrmGetName(first_frame), errText);
        return EXIT_FAILURE;
    }

    if (ffopen(&ofptr, productName, READWRITE, &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Opening file %s: %s", productName, errText);
        status = 0;
        fits_close_file(ifptr, &status);
        return EXIT_FAILURE;
    }

    if (fits_write_date(ofptr, &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing DATE to file %s: %s", productName, errText);
        status = 0;
        fits_close_file(ifptr, &status);
        status = 0;
        fits_close_file(ofptr, &status);
        return EXIT_FAILURE;
    }


    /*
     * Copy only if mandatory keyword is present
     */

    for (i = 0; i < count_mandatory; i++) {
        if (fits_read_card(ifptr, mandatory[i].key, card, &status) == 0) {
            if (fits_update_card(ofptr, mandatory[i].key, card, &status)) {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Writing keyword %s to file %s: %s", 
                              mandatory[i].key, productName, errText);
                status = 0;
                fits_close_file(ifptr, &status);
                status = 0;
                fits_close_file(ofptr, &status);
                return EXIT_FAILURE;
            }
        }
        else
            status = 0;
    }


    /*
     * Here copy all the HIERARCH.ESO._ keywords, excluding the 
     * HIERARCH.ESO.DPR._ keywords, and of the .PRO._ and .DRS._
     * keywords.
     */

    if (fits_get_hdrspace(ifptr, &size, NULL, &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Counting header keywords in file %s: %s",
                      pilFrmGetName(first_frame), errText);
        status = 0;
        fits_close_file(ifptr, &status);
        status = 0;
        fits_close_file(ofptr, &status);
        return EXIT_FAILURE;
    }

    for (i = 1; i <= size; i++) {
        if (fits_read_keyn(ifptr, i, kname, value, NULL, &status))
            status = 0;
        if (strncmp(kname, "ESO", 3) == 0) {
            if (strncmp(kname, "ESO DPR ", 8) != 0 && 
                strncmp(kname, "ESO PRO ", 8) != 0 &&
                strncmp(kname, "ESO DRS ", 8) != 0) {

                fits_read_record(ifptr, i, card, &status);
                fits_update_card(ofptr, kname, card, &status);
                if (status) {
                    fits_get_errstatus(status, errText);
                    cpl_msg_error(fid, "Copying card keyword %s from file "
                                  "%s to file %s: %s", kname,
                                  pilFrmGetName(first_frame),
                                  productName, errText);
                    status = 0;
                    fits_close_file(ifptr, &status);
                    status = 0;
                    fits_close_file(ofptr, &status);
                    return EXIT_FAILURE;
                }
            }
        }
    }


    /*
     * Forbidden keywords are removed from product header.
     */

    fits_delete_key(ofptr, "ESO DPR CATG", &status);
    if (status) {
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ESO.DPR.CATG keyword from file %s: %s",
                          productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }
    fits_delete_key(ofptr, "ESO DPR TYPE", &status);
    if (status) {
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ESO.DPR.TYPE keyword from file %s: %s",
                          productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }
    fits_delete_key(ofptr, "ESO DPR TECH", &status);
    if (status) {
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ESO.DPR.TECH keyword from file %s: %s", 
                          productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }
    fits_delete_key(ofptr, "ARCFILE", &status);
    if (status) {
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ARCFILE keyword from file %s: %s",
                          productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }
    fits_delete_key(ofptr, "ORIGFILE", &status);
    if (status) {
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ORIGFILE keyword from file %s: %s",
                          productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }
    fits_delete_key(ofptr, "CHECKSUM", &status);
    if (status) {         
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing CHECKSUM keyword from file %s: %s",
                          productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;       
    }


    fits_delete_key(ofptr, "ESO DET OUT1 OVSCX", &status);
    if (status) {         
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ESO DET OUT1 OVSCX keyword "
                          "from file %s: %s", productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }


    fits_delete_key(ofptr, "ESO DET OUT1 OVSCY", &status);
    if (status) {
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ESO DET OUT1 OVSCY keyword "
                          "from file %s: %s", productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }


    fits_delete_key(ofptr, "ESO DET OUT1 PRSCX", &status);
    if (status) {
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ESO DET OUT1 PRSCX keyword "
                          "from file %s: %s", productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }


    fits_delete_key(ofptr, "ESO DET OUT1 PRSCY", &status);
    if (status) {
        if (status != KEY_NO_EXIST) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Removing ESO DET OUT1 PRSCY keyword "
                          "from file %s: %s", productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            status = 0;
            fits_close_file(ifptr, &status);
            return EXIT_FAILURE;
        }
        status = 0;
    }


    /*
     * DATAMD5
     */

    if (fits_update_key(ofptr, TSTRING, DATAMD5, md5,
        "MD5 signature of data product", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s", DATAMD5,
                      productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        status = 0;
        fits_close_file(ifptr, &status);
        return EXIT_FAILURE;
    }


    /* 
     * PIPEFILE 
     */

    if (fits_update_key(ofptr, TSTRING, PIPEFILE, productName, 
        "Filename of data product", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s", PIPEFILE,
                      productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        status = 0;
        fits_close_file(ifptr, &status);
        return EXIT_FAILURE;
    }


    /*
     * PRO DID
     */

    if (fits_update_key(ofptr, TSTRING, PRO_DID, "PRO-1.15",
        "Data dictionary for PRO", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s", PRO_DID,
                      productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        status = 0;
        fits_close_file(ifptr, &status);
        return EXIT_FAILURE;
    }


    /* 
     * PRO CATG 
     */

    if (fits_update_key(ofptr, TSTRING, PRO_CATG, 
                        (char *)pilFrmGetCategory(proFrame),
                        "Category of pipeline product frame", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s", PRO_CATG,
                      productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        status = 0;
        fits_close_file(ifptr, &status);
        return EXIT_FAILURE;
    }


    /*
     * PRO TYPE
     */

    if (fits_update_key(ofptr, TSTRING, PRO_TYPE, "REDUCED",
        "Product type", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s", PRO_TYPE,
                      productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        status = 0;
        fits_close_file(ifptr, &status);
        return EXIT_FAILURE;
    }


    /*
     * PRO TECH
     */

    fits_read_key(ifptr, TSTRING, DPR_TECH, value, NULL, &status);
    if (status) {
        status = 0;
        fits_read_key(ifptr, TSTRING, PRO_TECH, value, NULL, &status);
    }

    if (status) {
        fits_get_errstatus(status, errText);
        cpl_msg_debug(fid,
                      "Reading keyword PRO.TECH or DPR.TECH from file %s: %s",
                      pilFrmGetName(first_frame), errText);
        status = 0;
    }
    else if (fits_update_key(ofptr, TSTRING, PRO_TECH, value,
        "Observation Technique", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s",
                      PRO_TECH, productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        status = 0;
        fits_close_file(ifptr, &status);
        return EXIT_FAILURE;
    }

    /*
     * PRO SCIENCE
     */

    fits_read_key(ifptr, TSTRING, DPR_CATG, value, NULL, &status);
    if (status) {
        status = 0;
        fits_read_key(ifptr, TLOGICAL, PRO_SCIENCE, &pro_science, NULL, &status);
    }
    else {
        pro_science = !strncmp(value, "SCIENCE", 7);
    }

    if (status) {
        fits_get_errstatus(status, errText);
        cpl_msg_debug(fid,
               "Reading keyword PRO.SCIENCE or DPR.CATG from file %s: %s",
               pilFrmGetName(first_frame), errText);
        status = 0;
    }
    else if (fits_update_key(ofptr, TLOGICAL, PRO_SCIENCE, &pro_science,
        "Scientific product if T", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s",
                      PRO_SCIENCE, productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        status = 0;
        fits_close_file(ifptr, &status);
        return EXIT_FAILURE;
    }

    if (fits_close_file(ifptr, &status))             {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Closing file %s: %s",
                      pilFrmGetName(first_frame), errText);
        status = 0;
        fits_close_file(ofptr, &status);
        return EXIT_FAILURE;
    }

    /* 
     * PRO REC1 ID 
     */

    if (fits_update_key(ofptr, TSTRING, PRO_REC_ID, recipename,
        "Pipeline recipe (unique) identifier", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s", PRO_REC_ID,
                      productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        return EXIT_FAILURE;
    }


    /* 
     * PRO REC1 DRS ID 
     */

    if (fits_update_key(ofptr, TSTRING, PRO_REC_DRS_ID, 
                        cpl_version_get_version(),
                        "Data Reduction System identifier", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s", PRO_REC_DRS_ID,
                      productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        return EXIT_FAILURE;
    }


    /* 
     * PRO REC1 PIPE ID
     */

    snprintf(cval, FLEN_CARD, "%s/%s", PACKAGE, PACKAGE_VERSION);

    if (fits_update_key(ofptr, TSTRING, PRO_REC_PIPE_ID, cval,
        "Pipeline (unique) identifier", &status)) {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Writing keyword %s to file %s: %s", PRO_REC_PIPE_ID,
                      productName, errText);
        status = 0;
        fits_close_file(ofptr, &status);
        return EXIT_FAILURE;
    }


    /*
     * PRO REC1 RAWi NAME  and  PRO REC1 RAWi CATG
     */

    nraw = 0;
    frame = get_raw_frame_time_sequence(framelist, 0);
    // frame = pilSofFirst(framelist);

    while (frame != NULL) {
/**
        if (pilFrmGetType(frame) == PIL_FRAME_TYPE_RAW
            || pilFrmGetType(frame) == PIL_FRAME_TYPE_UNDEF) {
**/

            ++nraw;

            snprintf(cval, FLEN_CARD, PRO_REC_RAWi_NAME, nraw);

            if (fits_update_key(ofptr, TSTRING, cval, 
                                _get_base_name(pilFrmGetName(frame)), 
                                "File name of raw frame", &status)) {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Writing keyword %s to file %s: %s", cval,
                              productName, errText);
                status = 0;
                fits_close_file(ofptr, &status);
                return EXIT_FAILURE;
            }

            snprintf(cval, FLEN_CARD, PRO_REC_RAWi_CATG, nraw);

            if (fits_update_key(ofptr, TSTRING, cval, 
                                (char *)pilFrmGetCategory(frame), 
                                "Category of raw frame", &status)) {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Writing keyword %s to file %s: %s", cval,
                              productName, errText);
                status = 0;
                fits_close_file(ofptr, &status);
                return EXIT_FAILURE;
            }

/**
        }
**/
        frame = get_raw_frame_time_sequence(framelist, nraw);
        // frame = pilSofNext(framelist, frame);
    }


    /* 
     * PRO DATANCOM (written only if missing)
     */

    if (fits_read_card(ofptr, PRO_DATANCOM, card, &status)) {
        status = 0;
        if (fits_update_key(ofptr, TINT, PRO_DATANCOM, &nraw,
                            "Number of combined frames", &status)) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Writing keyword %s to file %s: %s", 
                          PRO_DATANCOM, productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            return EXIT_FAILURE;
        }
    }


    /*
     * PRO REC1 CALi NAME,  PRO REC1 CALi CATG,  and  PRO REC1 CALi DATAMD5
     */ 
 
    ncal = 0;
    frame = pilSofFirst(framelist);
 
    while (frame != NULL) {
        if (pilFrmGetType(frame) == PIL_FRAME_TYPE_CALIB) {

            ++ncal;
 
            snprintf(cval, FLEN_CARD, PRO_REC_CALi_NAME, ncal);

            if (fits_update_key(ofptr, TSTRING, cval,
                                _get_base_name(pilFrmGetName(frame)),
                                "File name of calibration frame", &status)) {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Writing keyword %s to file %s: %s", cval,
                              productName, errText);
                status = 0;
                fits_close_file(ofptr, &status);
                return EXIT_FAILURE;
            }

            snprintf(cval, FLEN_CARD, PRO_REC_CALi_CATG, ncal);

            if (fits_update_key(ofptr, TSTRING, cval, 
                                (char *)pilFrmGetCategory(frame),
                                "Category of calibration frame", &status)) {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Writing keyword %s to file %s: %s", cval,
                              productName, errText);
                status = 0;
                fits_close_file(ofptr, &status);
                return EXIT_FAILURE;
            }

            snprintf(cval, FLEN_CARD, PRO_REC_CALi_DATAMD5, ncal);

            if (ffopen(&fptr, pilFrmGetName(frame), READONLY, &status))
            {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Opening file %s: %s",
                              pilFrmGetName(frame), errText);
                status = 0;
                fits_close_file(ofptr, &status);
                return EXIT_FAILURE;
            }

            if (fits_read_key(fptr, TSTRING, "DATAMD5", value, NULL, &status))
            {
                fits_get_errstatus(status, errText);
                cpl_msg_debug(fid, "Reading keyword DATAMD5 from file %s: %s",
                              pilFrmGetName(frame), errText);
                status = 0;
            }
            else if (fits_update_key(ofptr, TSTRING, cval, value,
                "MD5 signature of calib frame", &status)) {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Writing keyword DATAMD5 to file %s: %s", 
                              productName, errText);
                status = 0;
                fits_close_file(ofptr, &status);
                status = 0;
                fits_close_file(fptr, &status);
                return EXIT_FAILURE;
            }

            if (fits_close_file(fptr, &status))             {
                fits_get_errstatus(status, errText);
                cpl_msg_error(fid, "Closing file %s: %s",
                              pilFrmGetName(frame), errText);
                status = 0;
                fits_close_file(ofptr, &status);
                return EXIT_FAILURE;
            }

        }
        frame = pilSofNext(framelist, frame);

    }



    parString = pilDfsDumpDBtoString(&npar);

    for (i = 0; i < npar; i++) {

        if (strchr(parString[i], '\'') != NULL)
            sscanf(parString[i], "%[^=]='%[^=']'", parName, parValue);
        else if (strchr(parString[i], '\"') != NULL)
            sscanf(parString[i], "%[^=]=\"%[^=\"]\"", parName, parValue);
        else
            sscanf(parString[i], "%[^=]=%[^=]", parName, parValue);

        snprintf(cval, FLEN_CARD, PRO_REC_PARAMi_NAME, i + 1);

        if (fits_update_key(ofptr, TSTRING, cval, parName,
                            "Input parameter name", &status)) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Writing keyword %s to file %s: %s", cval,
                          productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            for (i = 0; i < npar; i++)
                pil_free(parString[i]);
            pil_free(parString);
            return EXIT_FAILURE;
        }

        snprintf(cval, FLEN_CARD, PRO_REC_PARAMi_VALUE, i + 1);

        if (fits_update_key(ofptr, TSTRING, cval, parValue,
                            "Input parameter value", &status)) {
            fits_get_errstatus(status, errText);
            cpl_msg_error(fid, "Writing keyword %s to file %s: %s", cval,
                          productName, errText);
            status = 0;
            fits_close_file(ofptr, &status);
            for (i = 0; i < npar; i++)
                pil_free(parString[i]);
            pil_free(parString);
            return EXIT_FAILURE;
        }

    }

    for (i = 0; i < npar; i++)
        pil_free(parString[i]);
    pil_free(parString);

    if (fits_close_file(ofptr, &status))             {
        fits_get_errstatus(status, errText);
        cpl_msg_error(fid, "Closing file %s: %s", productName, errText);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

/*
 * Input is assumed to be an arc lamp image - no check is performed
 * to this end.
 * The array "times" should have at least 3 elements: 
 * times[0] is assigned the He lamp exposure time.
 * times[1] is assigned the Ne lamp exposure time.
 * times[2] is assigned the Ar lamp exposure time.
 * This function returns zero on success. Failure occurs in case
 * the related keywords are not found, or a lamp is marked "ON"
 * and yet it has exposure time 0.
 */

int getArcLampTimes(VimosImage *arcLamp, double *times)
{

  char string[80];
  int  time;
  int  i, j;

  for (i = 1; i <= 5; i++) {
    if (VM_FALSE == readStringDescriptor(arcLamp->descs, 
                                         pilTrnGetKeyword("LampName", i), 
                                         string, NULL)) {
        return 1;
    }

    j = 0;
    switch (string[0]) {
    case 'A': j++;         /* Argon  */
    case 'N': j++;         /* Neon   */
    case 'H': j++;         /* Helium */
    default : j--;         /* other  */
    }

    if (j < 0) 
        continue;

    if (VM_FALSE == readStringDescriptor(arcLamp->descs, 
                                         pilTrnGetKeyword("LampState", i), 
                                         string, NULL)) {
        return 1;
    }

    if (string[1] == 'F') {  /* Lamp is OFF */
        times[j] = 0;
        continue;
    }

    if (VM_FALSE == readIntDescriptor(arcLamp->descs, 
                                      pilTrnGetKeyword("LampTime", i), 
                                      &time, NULL)) {
        return 1;
    }

    if (time < 1)
        return 1;

    times[j] = time;

  }

// printf("Helium %f\n", times[0]);
// printf("Neon   %f\n", times[1]);
// printf("Argon  %f\n", times[2]);

  return 0;

}

/*----------------------------------------------------------------------------*/
/**
  @brief    Issue a banner with the pipeline version
 */
/*----------------------------------------------------------------------------*/
void vimos_print_banner(void)
{
    cpl_msg_info(__func__, "*****************************************");
    cpl_msg_info(__func__, "Welcome to VIMOS Pipeline release %s",
                 vimos_get_version());
    cpl_msg_info(__func__, "*****************************************");
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Returns the version of the pipeline
 */
/*----------------------------------------------------------------------------*/
const char * vimos_get_version(void)
{
    static const char version[100] = PACKAGE_VERSION; //Defined in config.h
    return version;
}

/**@}*/
