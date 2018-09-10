/* $Id: vmimage.c,v 1.3 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <vmfit.h>
#include <assert.h>

#include <fitsio.h>
#include <fitsio2.h>

#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmtypes.h"
#include "vmtable.h"
#include "vmimage.h"


#include <pilmessages.h>
#include "cpl.h"


/**
 * @defgroup vmimage Images
 *
 * The module provides the image datatype and the fundamental operations.
 */

/**@{*/

VimosImage *newImage(int xlen, int ylen, float *data) 
{
  const char  modName[] = "newImage";
  VimosImage *tImage;
  
  tImage = (VimosImage *) cpl_malloc(sizeof(VimosImage));

  /* check if space was allocated */
  if (tImage == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
    
  /* should not allocate space for the data!!!!  a VimosImage is only an
     interface to a Midas Image or a Fits Image. The pointer should point to
     the data of the Midas Imiage or of the Fits Image
  */ 
  tImage->data = data;  
  tImage->xlen = xlen;
  tImage->ylen = ylen;
  tImage->descs = NULL;
  tImage->fptr = NULL;
  
  return(tImage);
}

VimosImage *newImageAndAlloc(int xlen, int ylen)
{
  const char  modName[] = "newImageAndAlloc";
  VimosImage *tmpImage;
  
  tmpImage = newImage(xlen, ylen, NULL);
  if (tmpImage == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
    
  tmpImage->data = (float *)cpl_calloc(xlen * ylen, sizeof(float));
  if (tmpImage->data == NULL) {
    /* cleanup */
    deleteImage(tmpImage);
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  
  return(tmpImage);
  
}

void deleteImage(VimosImage *image) 
{
  if (image != NULL) {
    deleteAllDescriptors(image->descs);
    if (image->data != NULL) {
      cpl_free(image->data);
    }
    cpl_free(image);
  }
}

void deleteImageAndAlloc(VimosImage *image)
{
  if (image != NULL) {
    if (image->data != NULL) {
      cpl_free(image->data);
      image->data = NULL;
    }
  }
  deleteImage(image);
}

VimosImage *duplicateImage(VimosImage *imageIn)
{
  VimosImage *imageOut;
  int         i,nPixel;

  imageOut = newImageAndAlloc(imageIn->xlen, imageIn->ylen);
  nPixel = imageIn->xlen * imageIn->ylen;
  for (i=0; i<nPixel; i++) imageOut->data[i] = imageIn->data[i];
  return imageOut;
}

VimosBool readDescsFromFitsImage(VimosDescriptor **desc, VimosImage *image)
{
  const char modName[] = "readDescsFromFitsImage";
  int i,j,len,k;
  int numDescs;
 
  int status, keypos;
  char name[FLEN_KEYWORD], value[FLEN_VALUE], value1[FLEN_VALUE]; 
  char comment[FLEN_COMMENT], comment1[FLEN_COMMENT];
  char type[1];

  
  VimosDescriptor* tDesc=NULL;
  VimosDescriptor* lastDesc;
  

  lastDesc = NULL;
  
  if (!desc) {
    cpl_msg_error(modName, "NULL input descriptor");
    return(VM_FALSE);
  }
  
  status = 0;

  if (image->fptr == NULL) {
    cpl_msg_error(modName, "No pointer to the FITS file");
    return(VM_FALSE);
  }

  /* get no. of keywords */
  if (fits_get_hdrpos(image->fptr, &numDescs, &keypos, &status)) {
    cpl_msg_error(modName, "fits_get_hdrpos returned error code %d", status);
    return(VM_FALSE);
  }
      
  for (i = 1; i <= numDescs; i++) {

    if (fits_read_keyn(image->fptr, i, name, value, comment, &status)) {
      cpl_msg_error(modName, "fits_read_keyn returned error code %d", status); 
      return(VM_FALSE);
    }

/* Check for blank lines */
    if ( strlen(name) == 0 ) {
      status = 0;
      continue;
    }

    if (strncmp("HISTORY",name,7) == 0 || strncmp("COMMENT",name,7) == 0) {
      /* remove quotes from string */
      if (ffc2s(comment, (comment1), &status)) {
	cpl_msg_error(modName, "ffc2s returned error code %d", status); 
        return(VM_FALSE);
      }
      /* remove trailing blanks */
      len = strlen(comment1);
      j = 0;
      while (comment1[j] == ' ') {
	j++;
      }
      for (k=0; k<=len-j; k++) {
	comment1[k] = comment1[k+j];
      }
      if (!(strstr(comment1,"FITS (Flexible Image Transport System) format") 
	  || strstr(comment1,"Astrophysics Supplement Series v4")
	  || strstr(comment1,"Contact the NASA Science")
	  || strstr(comment1,"FITS Definition document"))) {
	tDesc = newStringDescriptor(name, comment1, "");
	if (tDesc == NULL) {
	  cpl_msg_error(modName, "newStringDescriptor returned a NULL");
	  return(VM_FALSE);
	}
      }
    } 
    else 
    {    
      if (fits_get_keytype (value, (type), &status)) {
        cpl_msg_error(modName, "fits_get_keytype returned error code %d", status);
        return(VM_FALSE);
      }

      if (*(type) == 'F') {
	tDesc = newDoubleDescriptor(name, atof(value), comment);
	if (tDesc == NULL) {
	  cpl_msg_error(modName, "newDoubleDescriptor returned a NULL");
	  return(VM_FALSE);
	}
      } else if (*(type) == 'I') {
	tDesc = newIntDescriptor(name, atol(value), comment);
	if (tDesc == NULL) {
	  cpl_msg_error(modName, "newIntDescriptor returned a NULL");
	  return(VM_FALSE);
	}
      } else if (*(type) == 'L') {
	if (*(value) == 'T') {
	  tDesc = newBoolDescriptor(name, VM_TRUE, comment);
	  if (tDesc == NULL) {
	    cpl_msg_error(modName, "newBoolDescriptor returned a NULL");
	    return(VM_FALSE);
	  }
	}
	if (*(value) == 'F') {
	  tDesc = newBoolDescriptor(name,VM_FALSE, comment);
	  if (tDesc == NULL) {
	    cpl_msg_error(modName, "newBoolDescriptor returned a NULL");
	    return(VM_FALSE);
	  }
	}
      } else if (*(type) == 'C') {
        /* remove quotes from string */
	if (ffc2s(value, (value1), &status)) {   
	  cpl_msg_error(modName, "ffc2s returned error code %d", status); 
          return(VM_FALSE);
        }
	tDesc = newStringDescriptor(name, value1, comment);
	if (tDesc == NULL) {
	  cpl_msg_error(modName, "newStringDescriptor returned a NULL");
	  return(VM_FALSE);
	}
      }
      else {
        cpl_msg_error(modName, "Unrecognized key type %s",type);
        return(VM_FALSE);
      }
    }

    if (tDesc != NULL) {
      if (lastDesc == NULL) {
        *desc = tDesc;
      }
      else {
        lastDesc->next = tDesc;
        tDesc->prev = lastDesc;
      }
      lastDesc = tDesc;
    }

  }

  return(VM_TRUE);
}
	


VimosBool writeDescsToFitsImage(VimosDescriptor *desc, VimosImage *image)
{
  const char modName[] = "writeDescsToFitsImage";
  int        tf;
  int        fits_type, status=0;
  
  VimosDescriptor* tDesc;
  
  tDesc = desc;
  
  if (!desc ) {
    cpl_msg_error(modName, "NULL input descriptor");
    return(VM_FALSE);
  }

  if (image->fptr == NULL) {
    cpl_msg_error(modName, "No pointer to fits file");
    return(VM_FALSE);
  }
  
  while (tDesc) {
    
    switch(tDesc->descType) {
      case VM_INT : 
	{
	  fits_type = TINT;
	  if (fits_update_key (image->fptr,fits_type,tDesc->descName,
			       &(tDesc->descValue->i),tDesc->descComment,
                               &status)) {
	    cpl_msg_error(modName, "fits_update_key returned error %d", status);
	    return(VM_FALSE);
	  }   
	  break;
	}
      case VM_BOOL : 
	{
	  fits_type = TLOGICAL;
	  if (tDesc->descValue->b == VM_TRUE) {
	    tf = 1;
	  }
	  else {
	    tf = 0;
	  }
	  if (fits_update_key (image->fptr,fits_type,tDesc->descName,
			       &tf,tDesc->descComment,&status)) {
	    cpl_msg_error(modName, "fits_update_key returned error %d", status);
	    return(VM_FALSE);
	  }      
	  break;
	}
      case VM_FLOAT :
	{
	  fits_type = TFLOAT;
	  if (fits_update_key (image->fptr,fits_type,tDesc->descName,
			       &(tDesc->descValue->f),tDesc->descComment,
			       &status)) {
	    cpl_msg_error(modName, "fits_update_key returned error %d", status);
	    return(VM_FALSE);
	  }   
	  break;
	}
      case VM_DOUBLE :
	{
	  fits_type = TDOUBLE;
	  if (fits_update_key (image->fptr,fits_type,tDesc->descName,
			       &(tDesc->descValue->d),tDesc->descComment,
                               &status)) {
	    cpl_msg_error(modName, "fits_update_key returned error %d", status);
	    return(VM_FALSE);
	  }     
	  break;
	}
      case VM_STRING :
	{
	  fits_type = TSTRING;
	  if (strncmp ((tDesc->descName),"HISTORY",7) == 0){
	    if (fits_write_history (image->fptr, tDesc->descValue->s, 
                                    &status)) {
	      cpl_msg_error(modName, 
                         "fits_write_history returned error %d", status);
	      return(VM_FALSE);
	    }
	    break;
	  } 
	  if (strncmp ((tDesc->descName),"COMMENT",7) == 0){
	    if (fits_write_comment (image->fptr, tDesc->descValue->s, 
				    &status)) {
	      cpl_msg_error(modName, 
                          "fits_write_comment returned error %d", status);
	      return(VM_FALSE);
	    }
	    break;
	  }
          if (fits_update_key (image->fptr,fits_type,tDesc->descName,
                        tDesc->descValue->s,tDesc->descComment,&status)) {
	    cpl_msg_error(modName, "fits_update_key returned error %d", status);
            return(VM_FALSE);
	    }  
	  break;
	}
      default :
	{
	  cpl_msg_error(modName, "Invalid type"); 
          return(VM_FALSE);
	}
    }
    tDesc = tDesc->next;
  
  }  
  
  return(VM_TRUE);
  
}

VimosBool loadFitsHeader(VimosImage *image)
{

  const char modName[] = "loadFitsHeader";
  int status = 0;
  int nfound;
  long naxes[2];


  if (!image)
    return VM_FALSE;

  if (fits_read_keys_lng(image->fptr, "NAXIS", 
                         1, 2, naxes, &nfound, &status)) {
     cpl_msg_debug(modName, "fits_read_keys_lng() returned error %d", status);
     return VM_FALSE;
  }

  image->xlen = naxes[0];
  image->ylen = naxes[1];

  return readDescsFromFitsImage(&(image->descs), image);

}

VimosBool loadFitsData(VimosImage *image)
{

  const char modName[] = "loadFitsData";
  int status = 0;
  int anynull;
  long fpixel = 1;
  long npixels;
  float nullval = 0;


  if (!image)
    return VM_FALSE;

  npixels  = image->xlen * image->ylen;

  cpl_free(image->data);             /* Avoid possible leaks */

  image->data = (float *)cpl_malloc(npixels * sizeof(float));

  if (image->data == NULL) {
    cpl_msg_debug(modName, "Allocation error!");
    return VM_FALSE;
  }

  if (fits_read_img(image->fptr, TFLOAT, fpixel, npixels, &nullval,
                    image->data, &anynull, &status)) {
    cpl_msg_debug(modName, "fits_read_img() returned error %d", status);
    return VM_FALSE;
  }

  return VM_TRUE;

}

/**
 * @memo
 *   Arithmetic operations between 2 images. Result is a new image.
 *
 * @return Pointer to the resulting image
 *
 * @param ima1      pointer to first image
 * @param ima2      pointer to second image
 * @param optype    operation type
 *
 * @doc 
 *   Arithmetic operations between 2 images. The 2 images must have the 
 *   same dimensions. The output is 1 newly allocated image. The aritmetic
 *   operations are VM_OPER_ADD, VM_OPER_SUB, VM_OPER_MUL, VM_OPER_DIV
 *
 * @author P. Sartoretti
 */

VimosImage *imageArith(VimosImage *ima1, VimosImage *ima2, VimosOperator optype)
{
  VimosImage      *ima_out;
  VimosUlong32     i;
  VimosUlong32     nbpix;             /* total number of pixels in images */
  register float  *p1, *p2, *outp;    /* pixel value in ima_out */
  char             modName[] = "imageArith";

/* Error handling */
  
  if ((ima1 == NULL) || (ima2 == NULL)) {
    cpl_msg_debug(modName, "NULL input images");
    return(NULL);
  }

/* Input images must have the same sizes */

  if ((ima1->xlen != ima2->xlen) || (ima1->ylen != ima2->ylen)) {
    cpl_msg_error(modName, "First image is %dx%d, second image is %dx%d:"
               " images of different sizes cannot be combined",
               ima1->xlen, ima1->ylen, ima2->xlen, ima2->ylen);
    return NULL;
  }

/* Allocation */

  ima_out = newImageAndAlloc(ima1->xlen, ima1->ylen);

/* Do the arithmetic operation */

  p1 = ima1->data;
  p2 = ima2->data;
  outp = ima_out->data;
  nbpix = ima1->xlen * ima1->ylen;


  switch(optype){
  
  case VM_OPER_ADD:
    for (i=0; i<nbpix; i++)
      *outp++ = *p1++ + *p2++;
    break;

  case VM_OPER_SUB:
    for (i=0; i<nbpix; i++)
      *outp++ = *p1++ - *p2++;
    break;

  case VM_OPER_MUL:
    for (i=0; i<nbpix; i++)
      *outp++ = *p1++ * *p2++;
    break;

  case VM_OPER_DIV:
    for (i=0; i<nbpix; i++) {
      if (fabs(*p2) < MIN_DIVISOR){
        *outp = MAX_PIXEL_VALUE;
         outp++;
         p1++;
         p2++;
      }
      else { 
        *outp++ = *p1++ / *p2++; 
      }
    }
  break;
          
  default:
    cpl_msg_error(modName, "Unrecognized operator");
    return NULL;
  }
  return(ima_out);
}

/**
 * @memo
 *   Arithmetic operations between an image and a constant. Result is a new
 *   image.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_in    pointer to image
 * @param constant  constant
 * @param optype    operation type
 *
 * @doc 
 *   Arithmetic operations between an image and a constant. 
 *   The result is a new image. The aritmetic
 *   operations are VM_OPER_ADD, VM_OPER_SUB, VM_OPER_MUL, VM_OPER_DIV
 *
 * @author P. Sartoretti
 */


VimosImage *constArith(VimosImage *ima_in, double constant, 
                        VimosOperator optype)
{
  VimosImage   *ima_out;
  VimosUlong32 i;
  VimosUlong32 nbpix;    /* total number of pixels in images */
  double       invconst;
  char         modName[] = "constArith";

/* Error handling   */

  if (ima_in == NULL) {
    cpl_msg_error(modName, "NULL input image");
    return(NULL);
  }
  
  if ((optype == VM_OPER_DIV) && (fabs((double)constant) < MIN_DIVISOR)) {
    cpl_msg_error(modName, "Division by zero");
    return NULL;  
  }

/* Allocations  */
  
  ima_out = newImageAndAlloc(ima_in->xlen, ima_in->ylen);
  nbpix = ima_in->xlen * ima_in->ylen;
    
  switch(optype) {
    
  case VM_OPER_ADD:
    for (i=0; i<nbpix; i++)
      ima_out->data[i] = (float)((double)ima_in->data[i] + constant);
    break;
                 
  case VM_OPER_SUB:
    for (i=0; i<nbpix; i++)
      ima_out->data[i] = (float)((double)ima_in->data[i] - constant);
    break;
                 
  case VM_OPER_MUL:
    for (i=0; i<nbpix; i++)
      ima_out->data[i] = (float)((double)ima_in->data[i] * constant);
    break;
                 
  case VM_OPER_DIV:
    /* Multiplications are faster than divisions !*/
    invconst = (double)1.0 / constant; 
    for (i=0; i<nbpix; i++)
      ima_out->data[i] = (float)((double)ima_in->data[i] * invconst);
    break;

  default:
    cpl_msg_error(modName, "Unrecognized operator");
    return NULL;
  }
 
  return(ima_out);
}

/**
 * @memo
 *   Arithmetic operations between 2 images. Result is stored in the 
 *   first image.
 *
 * @return 0 on success 1 on failure
 *
 * @param ima1      pointer to first image
 * @param ima2      pointer to second image
 * @param optype     operation type
 *
 * @doc 
 *   Arithmetic operations between 2 images. The 2 images must have the 
 *   same dimensions. The output is stored in the first image. The aritmetic
 *   operations are VM_OPER_ADD, VM_OPER_SUB, VM_OPER_MUL, VM_OPER_DIV
 *
 * @author P. Sartoretti
 */
   
int imageArithLocal(VimosImage *ima1, VimosImage *ima2, VimosOperator optype )
{
 
  register float          *p1, *p2;
  register VimosUlong32    i;
  VimosUlong32             nbpix;    
  char                     modName[] = "imageArithLocal";


/* Error handling */
  
  if ((ima1 == NULL) || (ima2 == NULL)) {
    cpl_msg_error(modName, "NULL input images");
    return 1;
  }

  nbpix = ima1->xlen * ima1->ylen; 

/* Input images must have the same sizes */

  if ((ima1->xlen != ima2->xlen) || (ima1->ylen != ima2->ylen)) {
    cpl_msg_error(modName, "First image is %dx%d, second image is %dx%d:"
               " images of different sizes cannot be combined",
               ima1->xlen, ima1->ylen, ima2->xlen, ima2->ylen);
    return 1;
  }

  p1 = ima1->data;
  p2 = ima2->data;

  switch(optype){
  
  case VM_OPER_ADD:
    for (i=0; i<nbpix; i++) *p1++ += *p2++;
    break;

  case VM_OPER_SUB:
    for (i=0; i<nbpix; i++) *p1++ -= *p2++;
    break;

  case VM_OPER_MUL:
    for (i=0; i<nbpix; i++) *p1++ *= *p2++;
    break;

  case VM_OPER_DIV:
    for (i=0; i<nbpix; i++) {
      if (fabs(*p2) < MIN_DIVISOR)
        *p1++ = MAX_PIXEL_VALUE;
      else
        *p1++ /= *p2++; 
    }
    break;
  
  default:
    cpl_msg_error(modName, "Unrecognized operator");
    return 1;
  }
  
  return 0;
}
/**
 * @memo
 *   Arithmetic operations between an image and a constant. Result is stored 
 *   in the input image.
 *
 * @return 0 on success 1 on failure
 *
 * @param ima_in    pointer to image
 * @param constant  constant
 * @param optype          operation type
 *
 * @doc 
 *   Arithmetic operations between an image and a constant. 
 *   The result is stored in the input image. The aritmetic
 *   operations are VM_OPER_ADD, VM_OPER_SUB, VM_OPER_MUL, VM_OPER_DIV
 *
 * @author P. Sartoretti
 */

int constArithLocal(VimosImage *ima_in, double constant, VimosOperator optype)

{
  VimosUlong32    i;
  VimosUlong32    nbpix;    /* total number of pixels in images */
  double          invconst;
  char            modName[] = "constArithLocal";
 
/* Error handling   */

  if (ima_in == NULL) {
    cpl_msg_error(modName, "NULL input image");
    return 1;
  }
 
  if ((fabs((double)constant) < MIN_DIVISOR) && (optype == VM_OPER_DIV)) {
    cpl_msg_error(modName, "division by zero");
    return 1;  
  }

  nbpix = ima_in->xlen * ima_in->ylen;
    
  switch(optype) {
    
  case VM_OPER_ADD:
    for (i=0; i<nbpix; i++)
      ima_in->data[i] = (float)((double)ima_in->data[i] + constant);
    break;
                 
  case VM_OPER_SUB:
    for (i=0; i<nbpix; i++)
      ima_in->data[i] = (float)((double)ima_in->data[i] - constant);
    break;
                 
  case VM_OPER_MUL:
    for (i=0; i<nbpix; i++)
      ima_in->data[i] = (float)((double)ima_in->data[i] * constant);
    break;
                 
  case VM_OPER_DIV:
    /* Multiplications are faster than divisions !*/
    invconst = (double)1.0 / constant; 
    for (i=0; i<nbpix; i++)
      ima_in->data[i] = (float)((double)ima_in->data[i] * invconst);
     break;

  default:
    cpl_msg_error(modName, "Unrecognized operator");
    return 1;
  }
 
 return 0;
}
/**
 * @memo
 *   Apply a average filter to an image.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_in        pointer to input image
 * @param filtsizex     size in x of the filter box
 * @param filtsizey     size in y of the filter box
 * @param excludeCenter flag to exclude the center of the box in the
 *                      computation of the average
 *
 * @doc 
 *   Apply an optimized average filter of sizes filtsizex, filtsixey to
 *   an image. The output is a newly allocated image. Set excludeCenter
 *   flag to 1 to exclude the center of the median kernel from the median
 *   computation. If the median kernel extends over the image edges, the
 *   boundary points are repeated such that the box has always the same
 *   size (this is a modification compared to the original algorithm).
 *
 * @author C.Izzo. The original algorithm is taken from Eclipse (N.Devillard)
 */

VimosImage *VmFrAveFil(VimosImage *ima_in, int filtsizex, int filtsizey, 
            int excludeCenter) 
{
  VimosImage   *filt_img = NULL;
  int           col, row;
  float        *buf = NULL;
  int           avesize, upright_x, loleft_x, upright_y, loleft_y;
  int           uprightuse_x, loleftuse_x;
  int           i, j;
  int           xIsEven = !(filtsizex - (filtsizex/2)*2);
  int           yIsEven = !(filtsizey - (filtsizey/2)*2);
  float        *inpt;
  float        *outpt;
  int           f2x, f2y;
  char          modName[] = "VmFrAveFil";

  if (xIsEven) filtsizex++;
  if (yIsEven) filtsizey++;

  cpl_msg_debug(modName,
  "Filtering image using method AVERAGE, box %dx%d\n", filtsizex, filtsizey);

  if (ima_in->xlen <= filtsizex || ima_in->ylen <= filtsizey) {
    cpl_msg_error(modName, "Average filter size: %dx%d, image size: %d,%d",
           filtsizex, filtsizey, ima_in->xlen, ima_in->ylen);
    return NULL;
  }

  if (excludeCenter) excludeCenter = 1;

  f2x = filtsizex/2;
  f2y = filtsizey/2;
  filt_img = newImageAndAlloc(ima_in->xlen, ima_in->ylen);
  buf = cpl_malloc(filtsizex*filtsizey*sizeof(float));

  for (row=0; row<ima_in->ylen; row++){
    loleft_y = row - f2y;
    upright_y = row + f2y + 1;

    for (col=0; col<ima_in->xlen; col++){
      loleft_x = col - f2x;
      loleftuse_x = MAX(loleft_x,0); /* Lowest x-value on image */
      avesize = filtsizex * filtsizey - excludeCenter;
      upright_x = col + f2x + 1;
      uprightuse_x = MIN(upright_x,ima_in->xlen); /* Highest x-value on image */

      /* Optimized extraction */
      outpt = buf;
      if (excludeCenter) {
        for (j=loleft_y; j<upright_y; j++){

          if (j < 0)
            inpt = ima_in->data+loleftuse_x;
          else if (j > ima_in->ylen-1)
            inpt = ima_in->data+loleftuse_x + (ima_in->ylen-1)*ima_in->xlen;
          else
            inpt = ima_in->data+loleftuse_x + j*ima_in->xlen;

          for (i=loleft_x; i<loleftuse_x; i++)
            *outpt++ = *inpt;

          for (i=loleftuse_x; i<uprightuse_x; i++){
            if (i == col && j == row)
              inpt++;                /*** Skip "central" pixel value ***/
            else
              *outpt++ = *inpt++;
          }

          for (i=uprightuse_x; i<upright_x; i++)
            *outpt++ = *inpt;
        }
      }
      else {
        for (j=loleft_y; j<upright_y; j++){

          if (j < 0)
            inpt = ima_in->data+loleftuse_x;
          else if (j > ima_in->ylen-1)
            inpt = ima_in->data+loleftuse_x + (ima_in->ylen-1)*ima_in->xlen;
          else
            inpt = ima_in->data+loleftuse_x + j*ima_in->xlen;

          for (i=loleft_x; i<loleftuse_x; i++)
            *outpt++ = *inpt;

          for (i=loleftuse_x; i<uprightuse_x; i++)
            *outpt++ = *inpt++;

          for (i=uprightuse_x; i<upright_x; i++)
            *outpt++ = *inpt;
        }
      }

      filt_img->data[col+row*filt_img->xlen] = computeAverageFloat(buf,avesize);
    }
  }
  cpl_free(buf);
  return filt_img;
}

/**
 * @memo
 *   Apply a median filter to an image.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_in        pointer to input image
 * @param filtsizex     size in x of the filter box
 * @param filtsizey     size in y of the filter box
 * @param excludecenter flag to exclude the center of the box in the
 *                      computation of the median
 *
 * @doc 
 *   Apply an optimized median filter of sizes filtsizex, filtsixey 
 *   to an image. The output is a newly allocated image. Set excludeCenter
 *   flag to 1 to exclude the center of the median kernel from the median 
 *   computation. If the median kernel extends over the image edges, the 
 *   boundary points are repeated such that the box has always the same 
 *   size (this is a modification compared to the original algorithm).
 *
 * @author P. Sartoretti, C.Izzo. Algorithm taken from Eclipse (N.Devillard)
 */

VimosImage *VmFrMedFil(VimosImage *ima_in, int filtsizex, int filtsizey, 
            int excludeCenter) 
{
  VimosImage   *filt_img = NULL;
  int           col, row;
  float        *buf = NULL;
  int           medsize, upright_x, loleft_x, upright_y, loleft_y;
  int           uprightuse_x, loleftuse_x;
  int           i, j;
  int           xIsEven = !(filtsizex - (filtsizex/2)*2);
  int           yIsEven = !(filtsizey - (filtsizey/2)*2);
  float        *inpt;
  float        *outpt;
  int           f2x, f2y;
  char          modName[] = "VmFrMedFil";

  if (xIsEven) filtsizex++;
  if (yIsEven) filtsizey++;

  cpl_msg_debug(modName,
  "Filtering image using method MEDIAN, box %dx%d\n", filtsizex, filtsizey);

  if (ima_in->xlen <= filtsizex || ima_in->ylen <= filtsizey) {
    cpl_msg_error(modName, "Median filter size: %dx%d, image size: %d,%d",
           filtsizex, filtsizey, ima_in->xlen, ima_in->ylen);
    return NULL;
  }

  if (excludeCenter) excludeCenter = 1;

  f2x = filtsizex/2;
  f2y = filtsizey/2;
  filt_img = newImageAndAlloc(ima_in->xlen, ima_in->ylen);
  buf = cpl_malloc(filtsizex*filtsizey*sizeof(float));

  for (row=0; row<ima_in->ylen; row++){
    loleft_y = row - f2y;
    upright_y = row + f2y + 1;

    for (col=0; col<ima_in->xlen; col++){
      loleft_x = col - f2x;
      loleftuse_x = MAX(loleft_x,0); /* Lowest x-value on image */
      medsize = filtsizex * filtsizey - excludeCenter;
      upright_x = col + f2x + 1;
      uprightuse_x = MIN(upright_x,ima_in->xlen-1); /* Highest x-value on image */

      /* Optimized extraction */
      outpt = buf;
      if (excludeCenter) {
        for (j=loleft_y; j<upright_y; j++){

          if (j < 0)
            inpt = ima_in->data+loleftuse_x;
          else if (j > ima_in->ylen-1)
            inpt = ima_in->data+loleftuse_x + (ima_in->ylen-1)*ima_in->xlen;
          else
            inpt = ima_in->data+loleftuse_x + j*ima_in->xlen;

          for (i=loleft_x; i<loleftuse_x; i++)
            *outpt++ = *inpt;

          for (i=loleftuse_x; i<uprightuse_x; i++){
            if (i == col && j == row)
              inpt++;                /*** Skip "central" pixel value ***/
            else
              *outpt++ = *inpt++;
          }

          for (i=uprightuse_x; i<upright_x; i++)
            *outpt++ = *inpt;
        }
      }
      else {
        for (j=loleft_y; j<upright_y; j++){

          if (j < 0)
            inpt = ima_in->data+loleftuse_x;
          else if (j > ima_in->ylen-1)
            inpt = ima_in->data+loleftuse_x + (ima_in->ylen-1)*ima_in->xlen;
          else
            inpt = ima_in->data+loleftuse_x + j*ima_in->xlen;

          for (i=loleft_x; i<loleftuse_x; i++)
            *outpt++ = *inpt;

          for (i=loleftuse_x; i<uprightuse_x; i++)
            *outpt++ = *inpt++;

          for (i=uprightuse_x; i<upright_x; i++)
            *outpt++ = *inpt;
        }
      }

      filt_img->data[col+row*filt_img->xlen] = medianPixelvalue(buf,medsize);
    }
  }
  cpl_free(buf);
  return filt_img;
}

/**
 * @memo
 *   Filter an image.
 *
 * @return Pointer to the resulting image
 *
 * @param image         pointer to input image
 * @param xSize         size in x of the filter box
 * @param ySize         size in y of the filter box
 * @param method        filter method
 * @param excludecenter flag to exclude the center of the box in the
 *                      computation of the median
 *
 * @doc
 *   Apply a filter of choice to an image. Available method are
 *   median, average, gaussian filter. An automatic method can
 *   be specified.
 *
 * @see VmFrMedFil, VmFrAveFil
 *
 * @author C.Izzo and P.Sartoretti
 */

VimosImage *VmFrFilter(VimosImage *image, int xSize, int ySize, 
                                 FilterMethod method, int excludeCenter)
{
  const char  modName[] = "VmFrFilter";
  VimosImage *outImage;

  switch (method) {
  case FILTER_AUTO: 
  case FILTER_MEDIAN: 
    outImage = VmFrMedFil(image, xSize, ySize, excludeCenter);
    break;
  case FILTER_AVERAGE:
    outImage = VmFrAveFil(image, xSize, ySize, excludeCenter);
    break;
  default:
    cpl_msg_warning(modName, "Filter method not yet implemented - using MEDIAN");
    outImage = VmFrMedFil(image, xSize, ySize, excludeCenter);
  }
  return outImage;
}

/**
 * @memo
 *   Normalize an image.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_in        pointer to input image
 * @param meth          normalization method
 *
 * @doc 
 *   Normalize an image. Possible normalization methods are: MEAN 
 *   normalize image to its mean value, MEDIAN: normalize image to its 
 *   median value, and MODE: normalize image to its mode. The output 
 *   is a newly allocated image.
 *
 * @author P. Sartoretti
 */

VimosImage *VmImNorm(VimosImage *ima_in, Method meth)
{
  VimosImage  *ima_out;
  float        divisor;
  char         modName[] = "VmImNorm";
  
  /* Error handling */
  
  if (ima_in == NULL) {
    cpl_msg_error(modName, "NULL input image");
    return NULL;
  }

  switch(meth) {
    case MEAN  : divisor = imageMean(ima_in); break;
   
    case MEDIAN: divisor = imageMedian(ima_in); break;
    
    case MODE  : divisor = imageMode(ima_in); break; 
   
    default    : cpl_msg_error(modName, "Unrecognized normalization method");
                 return NULL;
  }

  if (fabs(divisor) < MIN_DIVISOR) {
    cpl_msg_error(modName, "Division by zero");
    return NULL;
  }
  ima_out = constArith(ima_in, divisor, VM_OPER_DIV);

  return ima_out;
}
/**
 * @memo
 *   Combine a list of images.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_list      Pointer to the array of images to combine
 * @param num           Number of images to combine
 * @param combMethod    Combination method
 * @param combParameter Used just when combination the method requires extra
 *                      parameters, as in COMB_REJECT or COMB_KSIGMA.
 * @param flag          Reject or not -32000 values in input
 *
 * @doc
 *   Combine a list of images according to the method. Possible methods are
 *   COMB_KSIGMA, COMB_MEDIAN, COMB_AVERAGE, COMB_SUM, COMB_REJECT, and AUTO.
 *   AUTO automatically selects the combination method based on the number
 *   of frames
 *
 * @author P. Sartoretti, C.Izzo
 */

VimosImage *frComb(VimosImage **ima_list, int num, 
            CombMethod combMethod, CombParameters *combParameter, int flag)
{
  VimosImage *ima_out;
  char        modName[] = "frComb";
 
  if ((combMethod == COMB_KSIGMA)) {
    cpl_msg_debug(modName, "Combination method is k-sigma clipping: "
                        "low, K = %3.1f sigma, high, K = %3.1f sigma", 
                        combParameter->kSigmaLow, combParameter->kSigmaHigh);
    ima_out = frCombKSigma(ima_list, combParameter->kSigmaLow, 
                                     combParameter->kSigmaHigh, num);
  }
  else if ((combMethod == COMB_REJECT)) {
    cpl_msg_debug(modName, "Combination method is MINMAX REJECTION:"
       "rejecting lower %d and upper %d pixel values over %d", 
       (int)(floor(num * combParameter->minRejection/100.)) + 1,
       (int)(floor(num * combParameter->maxRejection/100.)) + 1, num);
    ima_out = frCombMinMaxReject(ima_list, 
              combParameter->minRejection, combParameter->maxRejection, num);
  }
  else if ((combMethod == COMB_MEDIAN)) {
    cpl_msg_debug(modName, "Combination method is MEDIAN");
    ima_out = frCombMedian(ima_list, num, flag);
  }
  else if ((combMethod == COMB_AVERAGE)) {
    cpl_msg_debug(modName, "Combination method is AVERAGE");
    ima_out = frCombAverage(ima_list,num);
  }
  else if ((combMethod == COMB_SUM)) {
    cpl_msg_debug(modName, "Combination method is SUM");
    ima_out = frCombSum(ima_list, num);
  }
  else {
    cpl_msg_error(modName, "Unrecognized combination method");
    ima_out = NULL;
  }
  return ima_out;
}


VimosImage *frComb32000(VimosImage **ima_list, int num,
            CombMethod combMethod, CombParameters *combParameter, int flag)
{
  VimosImage *ima_out;
  char        modName[] = "frComb";
 
  if ((combMethod == COMB_KSIGMA)) {
    cpl_msg_debug(modName, "Combination method is k-sigma clipping: "
                        "low, K = %3.1f sigma, high, K = %3.1f sigma",
                        combParameter->kSigmaLow, combParameter->kSigmaHigh);
    ima_out = frCombKSigma32000(ima_list, combParameter->kSigmaLow,
                                     combParameter->kSigmaHigh, num);
  }
  else if ((combMethod == COMB_REJECT)) {
    cpl_msg_debug(modName, "Combination method is MINMAX REJECTION:"
       "rejecting lower %d and upper %d pixel values over %d",
       (int)(floor(num * combParameter->minRejection/100.)) + 1,
       (int)(floor(num * combParameter->maxRejection/100.)) + 1, num);
    ima_out = frCombMinMaxReject32000(ima_list,
              combParameter->minRejection, combParameter->maxRejection, num);
  }
  else if ((combMethod == COMB_MEDIAN)) {
    cpl_msg_debug(modName, "Combination method is MEDIAN");
    ima_out = frCombMedian(ima_list, num, flag);
  }
  else if ((combMethod == COMB_AVERAGE)) {
    cpl_msg_debug(modName, "Combination method is AVERAGE");
    ima_out = frCombAverage32000(ima_list,num);
  }
  else if ((combMethod == COMB_SUM)) {
    cpl_msg_debug(modName, "Combination method is SUM");
    ima_out = frCombSum(ima_list, num);
  }
  else {
    cpl_msg_error(modName, "Unrecognized combination method");
    ima_out = NULL;
  }
  return ima_out;
}
  
/**
 * @memo
 *   Sum a list of images.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_list        pointer to the array of images to combine
 * @param num             number of images
 *
 * @doc
 *   Sum a list of images. Result is a newly allocated image.
 *
 * @author P. Sartoretti
 */

VimosImage *frCombSum(VimosImage **ima_list, int num)
{ 
  VimosImage *ima_sum; 
  int         xlen, ylen, npix;
  int         i;
  char        modName[] = "frCombSum";
  
  /* Error handling */
  
  if (*ima_list == NULL) {
    cpl_msg_error(modName, "NULL input list");
    return NULL;
  }

/* Check dimensions of images to combine */

  xlen = ima_list[0]->xlen;
  ylen = ima_list[0]->ylen;
  npix = xlen*ylen;
  for (i=1; i<num; i++) {
    if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
      cpl_msg_error(modName, "Images must have the same dimensions");
      return NULL;
    }
  }
  
  ima_sum = newImageAndAlloc(xlen, ylen);

  for (i=0; i<npix; i++) ima_sum->data[i] = 0.;

  for (i=0; i<num; i++) imageArithLocal(ima_sum, ima_list[i], VM_OPER_ADD);
  
  return(ima_sum);
}

/**
 * @memo
 *   Compute average image of a list of images.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_list        pointer to the array of images to combine
 * @param num             number of images
 *
 * @doc
 *   Average a list of images. Result is a newly allocated image.
 *
 * @author P. Sartoretti
 */

VimosImage *frCombAverage(VimosImage **ima_list, int num)
{ 
  VimosImage    *ima_ave;
  float          inv;
  int            xlen,ylen,npix,i;
  char           modName[] = "frCombAverage";
 
 
 /* Error handling */
 
  if (ima_list == NULL) {
    cpl_msg_error(modName,"NULL input list");
    return(NULL);
  }

  xlen = ima_list[0]->xlen;
  ylen = ima_list[0]->ylen;
  npix = xlen*ylen;
  for (i=1; i<num; i++) {
    if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
      cpl_msg_error(modName, "Images must have the same dimensions");
      return NULL;
    }
  }
 
  ima_ave = newImageAndAlloc(xlen, ylen);

  for (i = 0; i < num; i++) 
    imageArithLocal(ima_ave, ima_list[i], VM_OPER_ADD);

  inv = 1.0 / num;
  for (i = 0; i < npix; i++) 
    ima_ave->data[i] *= inv;
  
  return(ima_ave);
}


VimosImage *frCombAverage32000(VimosImage **ima_list, int num)
{
  VimosImage *ima_ave;
  int         xlen, ylen, i, j, n, rej;
  char        modName[] = "frCombAverage32000";
  double     *zvalue;
  float       value;


  /* Error handling */

  if (ima_list == NULL) {
    cpl_msg_error(modName,"NULL input list");
    return(NULL);
  }

  xlen = ima_list[0]->xlen;
  ylen = ima_list[0]->ylen;
  for (i = 1; i < num; i++) {
    if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
      cpl_msg_error(modName, "Images must have the same dimensions");
      return NULL;
    }
  }

  ima_ave = newImageAndAlloc(xlen, ylen);

  zvalue = (double *)cpl_calloc(num, sizeof(double));

  for (j = 0; j < ylen; j++) {
    for (i = 0; i < xlen; i++) {
      rej = 0;
      for (n = 0; n < num; n++) {
        value = ima_list[n]->data[i+j*xlen];
        if (fabs(value + 32000) > 0.001)
          zvalue[n-rej] = value;
        else
          rej++;
      }
      if (rej != num) 
        ima_ave->data[i+j*xlen] = computeAverageDouble(zvalue, num-rej);
      else 
        ima_ave->data[i+j*xlen] = -32000;
    }
  }

  cpl_free(zvalue);

  return(ima_ave);
}


/**
 * @memo
 *   Compute median image of a list of images.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_list        pointer to the array of images to combine
 * @param num             number of images
 *
 * @doc
 *   Compute median image of a list of images. Result is a newly 
 *   allocated image.
 *
 * @author P. Sartoretti
 */

VimosImage *frCombMedian(VimosImage **ima_list, int num, int flag)
{
  VimosImage     *ima_med;
  float          *zvalue, value;
  int             i, j, n, xlen, ylen, rej;
  char            modName[] = "frCombMedian";
  
  if (ima_list == NULL) {
    cpl_msg_error(modName, "NULL input list");
    return(NULL);
  }
  
  /* Check dimensions of images to combine */  
  
  xlen = ima_list[0]->xlen;
  ylen = ima_list[0]->ylen;
  for (i=1; i<num; i++) {
    if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
      cpl_msg_error(modName, "Images must have the same dimensions");
      return NULL;
    }
  }

  if (num < MIN_FRAMES_MEDIAN) {
    cpl_msg_error(modName, 
                "At least %d frames are needed for median computation",
                MIN_FRAMES_MEDIAN);
    return NULL;
  }
  else {
    ima_med = newImageAndAlloc(xlen, ylen); 

    /* Allocate pointer and find pixel intensity of each frame
       corresponding to the position x,y. Compute median value.  
       If flag=1 reject -32000 pixel */   
  
    zvalue =(float*)cpl_calloc(num, sizeof(float));

    if(flag) {
      for (j=0; j<ylen; j++) {
	for (i=0; i<xlen; i++) {
	  rej=0;
	  for (n=0; n<num; n++) {
	    value = ima_list[n]->data[i+j*xlen];
	    if (fabs(value+32000) > 0.001) 
              zvalue[n-rej] = ima_list[n]->data[i+j*xlen];
	    else 
              rej++;
	  }
	  if(rej!=num) ima_med->data[i+j*xlen] = medianPixelvalue(zvalue, num-rej);
	  else ima_med->data[i+j*xlen] = -32000;
	}
      }
    } else {
      for (j=0; j<ylen; j++) {
	for (i=0; i<xlen; i++) {
	  for (n=0; n<num; n++) {
	    zvalue[n] = ima_list[n]->data[i+j*xlen];
	  }
	  ima_med->data[i+j*xlen] = medianPixelvalue(zvalue, num);
	}
      }
    }
    cpl_free(zvalue);
  }
  return ima_med;
}


/**
 * @memo
 *   Combine images with rejection of highest/lowest pixels values.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_list        pointer to the array of images to combine
 * @param minrej          number of pixels to min-reject
 * @param maxrej          number of pixels to max-reject
 * @param num             number of images
 *
 * @doc
 *   A zero value for both minrej and maxrej is illegal, and in that case 
 *   a NULL pointer is returned. If the total number of pixels to be 
 *   rejected is greater then or equal to the number of input images, 
 *   nothing is done and a NULL pointer is returned.
 *
 * @author P.Sartoretti, C.Izzo
 */

VimosImage *frCombMinMaxReject(VimosImage **ima_list, int minrej,
                               int maxrej, int num)
{
  char            modName[] = "frCombMinMaxReject";

  VimosImage     *ima_out;
  float          *zvalue;
  float           average = 0;
  int             i, j, n, xlen, ylen, minn, maxn;

  /* Errors handling */

  if (ima_list == NULL) {
    cpl_msg_error(modName, "NULL input list");
    return(NULL);
  }

  if (num < MIN_FRAMES_REJECT) {
    cpl_msg_error(modName,
                "No rejection with less than %d frames", MIN_FRAMES_REJECT);
    return NULL;
  }
  else {
    if ((minrej + maxrej) >= num) {
      cpl_msg_error(modName, "Max %d values can be rejected", num - 1);
      return NULL;
    }
    if (minrej == 0 && maxrej == 0) {
      cpl_msg_error(modName, "At least one value should be rejected");
      return NULL;
    }

    /* Check dimensions of images to combine */

    xlen = ima_list[0]->xlen;
    ylen = ima_list[0]->ylen;
    for (i = 1; i < num; i++) {
      if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
        cpl_msg_error(modName, "Images must have the same dimensions");
        return NULL;
      }
    }

    /* Allocate output image */

    ima_out = newImageAndAlloc(xlen, ylen);

    /* Find range number of pixel to average */

    minn = minrej;
    maxn = num - maxrej;

    /* Find pixel intensity in each frame.
       Sort pixel intensities, reject min/max and average. */

    zvalue = (float *) cpl_calloc(num, sizeof(float));

    for (j = 0; j < ylen; j++) {
      for (i = 0; i < xlen; i++) {
        for (n = 0; n < num; n++) {
          zvalue[n] = ima_list[n]->data[i + j * xlen];
        }

        sort(num, zvalue);
        average = 0;

        for (n = minn; n < maxn; n++)
          average += zvalue[n];

        ima_out->data[i + j * xlen] = average / (maxn - minn);
      }
    }
    cpl_free(zvalue);
  }
  return ima_out;
}


VimosImage *frCombMinMaxReject32000(VimosImage **ima_list, int minrej,
                               int maxrej, int num)
{
  char            modName[] = "frCombMinMaxReject";

  VimosImage     *ima_out;
  float          *zvalue;
  float           average = 0;
  int             i, j, n, rej, xlen, ylen, minn, maxn;

  /* Errors handling */

  if (ima_list == NULL) {
    cpl_msg_error(modName, "NULL input list");
    return(NULL);
  }

  if (num < MIN_FRAMES_REJECT) {
    cpl_msg_error(modName,
                "No rejection with less than %d frames", MIN_FRAMES_REJECT);
    return NULL;
  }
  else {
    if ((minrej + maxrej) >= num) {
      cpl_msg_error(modName, "Max %d values can be rejected", num - 1);
      return NULL;
    }
    if (minrej == 0 && maxrej == 0) {
      cpl_msg_error(modName, "At least one value should be rejected");
      return NULL;
    }

    /* Check dimensions of images to combine */

    xlen = ima_list[0]->xlen;
    ylen = ima_list[0]->ylen;
    for (i = 1; i < num; i++) {
      if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
        cpl_msg_error(modName, "Images must have the same dimensions");
        return NULL;
      }
    }


    /* Allocate output image */

    ima_out = newImageAndAlloc(xlen, ylen);

    /* Find pixel intensity in each frame.
       Sort pixel intensities, reject min/max and average. */

    zvalue = (float *) cpl_calloc(num, sizeof(float));

    for (j = 0; j < ylen; j++) {
      for (i = 0; i < xlen; i++) {
        rej = 0;
        for (n = 0; n < num; n++) {
          if (fabs(ima_list[n]->data[i+j*xlen]+32000) > 0.001)
            zvalue[n-rej] = ima_list[n]->data[i+j*xlen];
          else
            rej++;
        }

        if (num - rej < MIN_FRAMES_REJECT) {
          if (rej != num)
            ima_out->data[i+j*xlen] = computeAverageFloat(zvalue, num-rej);
          else
            ima_out->data[i+j*xlen] = -32000;
        }
        else {

          sort(num-rej, zvalue);

          /* Find range number of pixel to average */

          minn = minrej;
          maxn = num - rej - maxrej;

          average = 0;
          for (n = minn; n < maxn; n++)
            average += zvalue[n];

          ima_out->data[i + j * xlen] = average / (maxn - minn);

        }
      }
    }
    cpl_free(zvalue);
  }
  return ima_out;
}


/**
 * @memo
 *   Combine images with rejection of high/low pixels values.
 *
 * @return Pointer to the resulting image
 *
 * @param ima_list        pointer to the array of images to combine
 * @param minrej          percent of pixels to min-reject
 * @param maxrej          percent of pixels to max-reject
 * @param num             number of images
 *
 * @doc
 *   Percentages are expressed with numbers from 0 to 100.  
 *   If the percent of pixels to max or min reject is less
 *   then the min percentage step (i.e. 100/num), no image
 *   will be min or max rejected. Not more than 90 percent
 *   of pixel values can be rejected: if such threshold is 
 *   trespassed, a NULL pointer is returned.
 *
 * @author P. Sartoretti
 */

VimosImage *OLDfrCombMinMaxReject(VimosImage **ima_list, double minrej, 
                               double maxrej, int num)
{
  char            modName[] = "frCombMinMaxReject";

  VimosImage     *ima_out; 
  double          maxTotalPercent = 90.0;
  float          *zvalue;
  float           average = 0;
  int             i, j, n, xlen, ylen, minn, maxn;
  
  /* Errors handling */
  
  if (ima_list == NULL) {
    cpl_msg_error(modName, "NULL input list");
    return(NULL);
  }
  
  if (num < MIN_FRAMES_REJECT) {
    cpl_msg_error(modName, 
                "No rejection with less than %d frames", MIN_FRAMES_REJECT);
    return NULL;
  }
  else {
    if ((minrej + maxrej) > maxTotalPercent) {
      cpl_msg_error(modName, 
                    "Rejection should not be over %f2.0%%", maxTotalPercent);
      return NULL; 
    }
    
    /* Check dimensions of images to combine */  
    
    xlen = ima_list[0]->xlen;
    ylen = ima_list[0]->ylen;
    for (i=1; i<num; i++) {
      if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
        cpl_msg_error(modName, "Images must have the same dimensions");
        return NULL;
      }
    }
    
    /* Allocate output image */
    
    ima_out = newImageAndAlloc(xlen, ylen);
    
    /* Find range number of pixel to average */

    minn = (int)(floor(num * minrej/100.));
    maxn = num - (int)(floor(num * maxrej/100.));
    
    /* Find pixel intensity in each frame. 
       Sort pixel intensities, reject min/max and average. */   
    
    zvalue = (float *) cpl_calloc(num, sizeof(float));
    
    for (j = 0; j < ylen; j++) {
      for (i = 0; i < xlen; i++) {
        for (n = 0; n < num; n++) {
          zvalue[n] = ima_list[n]->data[i + j * xlen];
        }
        sort(num, zvalue);
        average = 0;
        for (n = minn; n < maxn; n++) {
          average += zvalue[n];
        }
        ima_out->data[i + j * xlen] = average / (maxn - minn);
      }
    }
    cpl_free(zvalue);
  }
  return ima_out;
}


/**
 * @memo
 *   Combine a list of images with rejection of pixels outside a K*sigma
 *   threshold (sigma-clipping)
 *
 * @return Pointer to the resulting image
 *
 * @param ima_list        pointer to the array of images to combine
 * @param klow            number of sigma for low values sigma-clipping
 * @param khigh           number of sigma for high values sigma-clipping
 * @param num             number of images
 *
 * @doc
 *  Combine frame with rejection of pixels outside
 *  a k*sigma treshold. The threshold is calculated on the median
 *  value, and not on the average to limit the influence of bad pixels.
 *  Input a list of images; klow to set the low treshold (=median-klow*sigma),
 *  and khigh to set the high treshold (=median+khigh*sigma).
 *  Result is a newly allocated image.
 *
 * @author P. Sartoretti
 */

VimosImage *frCombKSigma(VimosImage **ima_list, double dklow, 
                           double dkhigh, int num)
{
  char            modName[] = "frCombKSigma";

  VimosImage     *ima_out;
  float           klow = (float) dklow;
  float           khigh = (float) dkhigh;
  float          *zvalue,sigma,low,high,median,ave;
  int             i, j, n, xlen, ylen, goodpix;
  
  if (ima_list == NULL) {
    cpl_msg_error(modName, "NULL input list");
    return(NULL);
  }
  
  /* Check dimensions of images to combine */  
  
  xlen = ima_list[0]->xlen;
  ylen = ima_list[0]->ylen;
  for (i=1; i<num; i++) {
    if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
      cpl_msg_error(modName, "Images must have the same dimensions");
      return NULL;
    }
  }
  
  if (num < MIN_FRAMES_KSIGMA) {
    cpl_msg_warning(modName,
    "No sigma rejection with less than %d frames", MIN_FRAMES_KSIGMA);
    return NULL;
  }
  else {
    /* 
     * Find pixel intensity of each frame corresponding to the 
     * position x,y. Compute median value and sigma on the median. 
     * Average pixels inside treshold.
     */   
    
    ima_out = newImageAndAlloc(xlen, ylen); 
    
    zvalue = (float*) cpl_calloc(num, sizeof(float));
    
    for (j=0; j<ylen; j++) {
      for (i=0; i<xlen; i++) {
        for (n=0; n<num; n++) {
          zvalue[n] = ima_list[n]->data[i+j*xlen];
        }
        median = medianPixelvalue(zvalue,num);
        sigma=0.;
        for (n=0; n<num; n++) {
          sigma += fabs(zvalue[n]-median);
        }
        sigma /= num;
        sigma *= MEANDEV_TO_SIGMA;

        low  = median - klow*sigma;
        high = median + khigh*sigma;
        ave=0.;
        goodpix=num;
        for (n=0; n<num; n++) {
          if ((zvalue[n] < low) || (zvalue[n] > high)) {
            --goodpix;
          }
          else {
            ave += zvalue[n];
          }
        }
        ima_out->data[i+j*xlen] = ave/goodpix;
      }
    }
    cpl_free(zvalue);
  }
  return ima_out;
}


VimosImage *frCombKSigma32000(VimosImage **ima_list, double dklow,
                           double dkhigh, int num)
{
  char            modName[] = "frCombKSigma32000";

  VimosImage     *ima_out;
  float           klow = (float) dklow;
  float           khigh = (float) dkhigh;
  float          *zvalue, sigma, low, high, median, ave;
  int             i, j, n, rej, xlen, ylen, goodpix;


  if (ima_list == NULL) {
    cpl_msg_error(modName, "NULL input list");
    return(NULL);
  }

  /* Check dimensions of images to combine */

  xlen = ima_list[0]->xlen;
  ylen = ima_list[0]->ylen;
  for (i = 1; i < num; i++) {
    if ((ima_list[i]->xlen != xlen) || (ima_list[i]->ylen != ylen)) {
      cpl_msg_error(modName, "Images must have the same dimensions");
      return NULL;
    }
  }

  if (num < MIN_FRAMES_KSIGMA) {
    cpl_msg_warning(modName,
    "No sigma rejection with less than %d frames", MIN_FRAMES_KSIGMA);
    return NULL;
  }
  else {

    /*
     * Find pixel intensity of each frame corresponding to the
     * position x,y. Compute median value and sigma on the median.
     * Average pixels inside treshold.
     */

    ima_out = newImageAndAlloc(xlen, ylen);

    zvalue = (float *)cpl_calloc(num, sizeof(float));

    for (j = 0; j < ylen; j++) {
      for (i = 0; i < xlen; i++) {
        rej = 0;
        for (n = 0; n < num; n++) {
          if (fabs(ima_list[n]->data[i+j*xlen]+32000) > 0.001)
            zvalue[n-rej] = ima_list[n]->data[i+j*xlen];
          else
            rej++;
        }

        if (num - rej < MIN_FRAMES_KSIGMA) {
          if (rej != num)
            ima_out->data[i+j*xlen] = computeAverageFloat(zvalue, num-rej);
          else
            ima_out->data[i+j*xlen] = -32000;
        }
        else {
          median = medianPixelvalue(zvalue, num);

          sigma = 0.;
          for (n = 0; n < num-rej; n++) {
            sigma += fabs(zvalue[n] - median);
          }
          sigma /= num - rej;
          sigma *= MEANDEV_TO_SIGMA;

          low  = median - klow*sigma;
          high = median + khigh*sigma;
          ave = 0.;
          goodpix = num;
          for (n = 0; n < num - rej; n++) {
            if ((zvalue[n] < low) || (zvalue[n] > high)) {
              --goodpix;
            }
            else {
              ave += zvalue[n];
            }
          }
          ima_out->data[i+j*xlen] = ave/goodpix;
        }
      }
    }
    cpl_free(zvalue);
  }
  return ima_out;
}


/**
 * @memo
 *   Compute the average value of an image
 *
 * @return the average value, or 0 in case of failure
 *
 * @param ima_in   input image
 *
 * @doc 
 *   Compute the average value of an image
 *
 * @author P. Sartoretti
 */

float imageMean(VimosImage *ima_in)
{
  VimosUlong32    i;
  VimosUlong32    nbpix;    /* total number of pixels in images */
  float           mean;
  char            modName[] = "imageMean";
  
  if (ima_in == NULL) {
    cpl_msg_error(modName, "NULL input image");
    return 0.;
  }
  nbpix = ima_in->xlen * ima_in->ylen;
  mean = 0.0;
  for (i=0; i<nbpix; i++) 
    mean += ima_in->data[i];
  
  mean /= nbpix;
  return mean;
} 
/**
 * @memo
 *   Compute the median value of an image
 *
 * @return the median value, or 0 in case of failure
 *
 * @param ima_in   input image
 *
 * @doc 
 *   Compute the median value of an image. Use median-wirth.
 *
 * @author P. Sartoretti
 */

float imageMedian(VimosImage *ima_in)
{
  VimosImage    *copy;/* copy image because median_WIRTH modifies the input */
  VimosUlong32   i,nbpix;  
  float                 median;
  char           modName[] = "imageMedian";
  
  if (ima_in == NULL) {
    cpl_msg_error(modName, "NULL input image");
    return 0.;
  }
  
  nbpix = ima_in->xlen * ima_in->ylen;
  copy = newImageAndAlloc(ima_in->xlen,ima_in->ylen);
  if (copy == NULL) {
    cpl_msg_error(modName, "Cannot copy image: aborting median search");
    return 0.;
  }
  for (i=0; i<nbpix; i++) copy->data[i] = ima_in->data[i]; 
  median = medianWirth(copy->data, nbpix);
  deleteImage(copy);
  return median;
}

/**
 * @memo
 *   Compute average deviation from any given value
 *
 * @return average deviation
 *
 * @param image   input image
 * @param level   reference value
 *
 * @doc 
 *   Compute average deviation from any given value
 *
 * @author C. Izzo
 */

float imageAverageDeviation(VimosImage *image, float level)
{
  VimosUlong32  i, nPix;
  float         averageDeviation = 0.;
  
  if (image == NULL) 
    return -1.;

  nPix = image->xlen * image->ylen;

  for (i = 0; i < nPix; i++)
    averageDeviation += fabs(image->data[i] - level);

  averageDeviation /= nPix;

  return averageDeviation;
}

/**
 * @memo
 *   Compute image sigma
 *
 * @return the sigma value, 
 *
 * @param image   input image
 *
 * @doc 
 *   Compute the sigma value of an image
 *
 * @author P. Sartoretti
 */

float imageSigma(VimosImage *image)
{
  return imageAverageDeviation(image, imageMean(image));
}

/**
 * @memo
 *   Compute image sigma using median value instead of mean estimate
 *
 * @return the sigma value 
 *
 * @param image   input image
 *
 * @doc 
 *   Compute the sigma value of an image using the median value instead of 
 *   the mean
 *
 * @author P. Sartoretti
 */

float imageMedSigma(VimosImage *image)
{
  return imageAverageDeviation(image, imageMedian(image));
}

/**
 * @memo
 *   Compute image minimum value of an image
 *
 * @return the maximum value
 *
 * @param ima_in   input image
 *
 * @doc 
 *   Compute the minimum value of an image 
 *
 * @author P. Sartoretti
 */

float imageMinimum(VimosImage *ima_in)
{
  VimosUlong32  i;
  VimosUlong32  nbpix;
  float         minimum;
  char          modName[] = "imageMinimum";
  
  if (ima_in == NULL) {
    cpl_msg_error(modName, "NULL input image");
    return 0.;
  }
  nbpix = ima_in->xlen * ima_in->ylen;
  minimum = ima_in->data[0];
  for (i=1; i<nbpix; i++) {
    if (ima_in->data[i] < minimum)
      minimum = ima_in->data[i];
  }
  return minimum;
}
/**
 * @memo
 *   Compute image maximum value of an image
 *
 * @return the maximum value
 *
 * @param ima_in   input image
 *
 * @doc 
 *   Compute the maximum value of an image 
 *
 * @author P. Sartoretti
 */

float imageMaximum(VimosImage  *ima_in)
{
  VimosUlong32  i;
  VimosUlong32  nbpix;
  float         maximum;
  char          modName[] = "imageMaximum";
  
  if (ima_in == NULL){
    cpl_msg_error(modName, "NULL input image");
    return 0.;
  }
  nbpix = ima_in->xlen * ima_in->ylen;
  maximum = ima_in->data[0]; 
  for (i=1; i<nbpix; i++) {
    if (ima_in->data[i] > maximum)
      maximum = ima_in->data[i];
  }
  return maximum;
}
/**
 * @memo
 *   Compute image histogram
 *
 * @return the histogram
 *
 * @param ima_in   input image
 * @param nbins    number of bins for the histogram
 *
 * @doc 
 *  Compute an image histogram. Input one image and the number of bins for
 *  the histogram. Output a list of [nbins] VimosDpoints, where x contains 
 *  the central pixelvalue for each bin, and  y contains the number of 
 *  pixels in this bin.
 *
 * @author P. Sartoretti. This algorithm is taken from Eclipse (N.Devillard)
 */

VimosDpoint *imageHistogram(VimosImage *ima_in, unsigned int nbins)
     
{
  float           max,min;
  VimosUlong32    i,nbpix;
  VimosUlong32   *h;      /* frequency of each bin_val */
  int             bin_val;
  double          bin_size;
  VimosDpoint    *histogram;
  
  max = imageMaximum(ima_in);
  min = imageMinimum(ima_in);
  nbpix = ima_in->xlen * ima_in->ylen;
  
  bin_size = (double)(max - min)/(double)nbins;
  h = (VimosUlong32*)cpl_calloc(nbins,sizeof(VimosUlong32));
 
  for (i=0; i<nbpix; i++) { 
    if (ima_in->data[i] >= max) {
      bin_val = nbins - 1;
    }
    else {
      bin_val = (int)((ima_in->data[i] - min) /bin_size);
    }
    h[bin_val]++;
  }
  histogram = newDpoint(nbins);
  for (i=0; i<nbins; i++) {
    histogram[i].x = min + (double)i * bin_size;
    histogram[i].y = (double)h[i];
  }
  cpl_free(h);
  return histogram;
}


/**
 * @memo
 *   Compute image mode value of an image
 *
 * @return the mode, or -1 in case of failure
 *
 * @param ima_in   input image
 *
 * @doc 
 *   Compute the mode value of an image (i.e. the histo peak position) 
 *   It depends on the number of bins (binzise).
 *   Here take the bin size = 1 data unit (i.e. bin number =(int) max - min) 
 *
 * @author P. Sartoretti
 */

float imageMode(VimosImage *ima_in)
{
  const char modName[] = "imageMode";

  VimosDpoint    *hist;

  size_t          bin_num;

  float           mode, imax, imin;
  double          fwhm;
  

  assert(ima_in);

  imax = imageMaximum(ima_in);
  imin = imageMinimum(ima_in);

  if (imin == imax)
    return imin;

  bin_num = (size_t)floor(imax - imin);

  if (bin_num <= 1 ) { 
    cpl_msg_error(modName, "Number of bins is too small");
    return -1.;
  }
  hist = imageHistogram(ima_in, bin_num);
  mode = histogramPeak(hist, &fwhm, bin_num);
  
  return mode;
}


/**
 * @memo
 *   Compute the sum of all pixelvalues in a rectangle.
 *
 * @return the sum, or -1 in case of failure
 *
 * @param ima_in   input image
 * @param x        x lower left corner of the rectangle
 * @param y        y lower left corner of the rectangle
 * @param nx       number of pixels in x
 * @param ny       number of pixels in y    
 *
 * @doc 
 *   Compute the sum of all pixelvalues in a rectangle. The rectangle is given
 *   by the lower left corner coordinates (x,y) and the number of pixels
 *   nx and ny.  
 *
 * @author P. Sartoretti
 */

float sumPixelsInImage(VimosImage *ima_in, int x, int y, int nx, int ny)
{
  float sum;
  int   i,j;
  char  modName[] = "sumPixelsInImage";
  
  if (ima_in == NULL) {
    return 0.;
  }
  
  /* check box coordinates */
  
  if (x<0 || y<0 || x+nx > ima_in->xlen || y+ny > ima_in->ylen 
  || nx < 0 || ny < 0) {
    cpl_msg_error(modName, 
    "Invalid rectangle coordinates: lower left is %d,%d and"
    " upper right is %d,%d", x, y, x+nx-1, y+ny-1);
    return 0.;
  }
  
  sum=0.;
  for (j=y; j<(y+ny); j++) {
    for (i=x; i<(x+nx); i++) {
      sum += ima_in->data[i+j*ima_in->xlen];
    }
  }  
  return(sum);
}    

/**
 * @memo
 *   Collapse a rectangle inside an image into one column/row.
 *
 * @return the collapsed array
 *
 * @param ima_in   input image
 * @param x        x lower left corner of the rectangle
 * @param y        y lower left corner of the rectangle
 * @param nx       number of pixels in x
 * @param ny       number of pixels in y
 * @param collapse direction of collape: ROW or COLUMNS        
 *
 * @doc
 *   collapse a rectangle inside an image into one column. Compute the sum
 *   of the rows or of the columns of a rectangle. Input an Image,
 *   the lower left corner coordinates (x,y) of the rectangle; the
 *   number of pixels (nx, ny) from (x,y), and the collapsing direction:
 *   collapse ROW: sum of rows and collapse COLUMN: sum of columns. 
 *
 * @author P. Sartoretti , C. Izzo
 */

float *collapse2Dto1D(VimosImage *ima_in, int x, int y, int nx, 
                      int ny, Direction collapse)
{
  float       sum, *a1d = NULL;
  int         i;
  char        modName[] = "collapse2Dto1D";
  
  if (ima_in == NULL) {
    return NULL;
  }
  
  /* check box coordinates */
  
  if (x<0 || y<0 || x+nx > ima_in->xlen || y+ny > ima_in->ylen 
  || nx < 0 || ny < 0) {
    cpl_msg_error(modName, 
    "Invalid rectangle coordinates: lower left is %d,%d and"
    " upper right is %d,%d", x, y, x+nx-1, y+ny-1);
    return NULL;
  }
  
  switch(collapse){
  case ROW:
    a1d = cpl_calloc(nx,sizeof(float)); 
    for (i=0; i<nx; i++) {
      sum = sumPixelsInImage(ima_in,x+i,y,1,ny);
      a1d[i]=sum;
    }
    break;
    
  case COLUMN:
    a1d = cpl_calloc(ny,sizeof(float)); 
    for (i=0;i<ny; i++) {
      sum = sumPixelsInImage(ima_in,x,y+i,nx,1);
      a1d[i]=sum;
    }
    break;

  default:
    cpl_msg_error(modName, 
    "Supported directions are COLUMN (sum columns) or ROW (sum rows)");
    break;
  }
  return a1d;
}

/**
 * @memo
 *   Refine the peak centers by computing the baricenter.
 *
 * @return the list of corrected positions
 *
 * @param ima_in   input image
 * @param p_list   list of pixel position to refine
 * @param npix     size of p_list array
 * @param r1       defines a disk in which the barycenter will be computed
 * @param r2       inner radius of the ring where the background is computed
 * @param r3       outer radius of the ring where the background is computed
 *
 * @doc 
 *  Refine the peak centers, by computing a barycenter
 *  The 3 doubles (r1, r2, r3) define 3 circle radiuses:
 *  r2 and r3 define a ring centered on each pixel in the list,
 *  around which the background will be estimated.
 *  r1 defines a disk in which the barycenter will be computed,
 *  i.e. the centroid pixel, weighted by pixel values, from
 *  which the background value has been subtracted:
 *  x_center = (1/sum(pixel[i]-bg))*(sum(x[i] * (pixel[i]-bg)))
 *  y_center = (1/sum(pixel[i]-bg))*(sum(y[i] * (pixel[i]-bg)))
 *
 * @author P. Sartoretti. Algorithm (fine_position_centers) taken from 
 *                        Eclipse (N.Devillard)
 */

VimosPixel *finePosition(VimosImage *ima_in, VimosPixel *p_list, int npix,
                         double r1, double r2, double r3)
{
  VimosPixel *list_out;
  int         pixnum;
  int         i, j;
  double      est_x, est_y;
  double      corr_x, corr_y;
  double      sum_weights;
  double      curpix;
  double      background;
  int              countpix;
  double      sq_r1, sq_r2, sq_r3;
  double      sq_rad;
  char        modName[] = "finePosition";
  
  sq_r1 = (double)(r1 * r1);
  sq_r2 = (double)(r2 * r2);
  sq_r3 = (double)(r3 * r3);
  
  
  /* Check input */
  
  if (ima_in == NULL) {
    cpl_msg_error(modName, "Input NULL image");
    return NULL;
  }
  
  if (p_list == (VimosPixel*)NULL) {
    cpl_msg_error(modName, "No pixel in list: cannot refine centroid positions");
    return NULL;
  }
  
  if ((r1<1.0) || (r2<1.0) || (r3<1.0) || (r2<r1) || (r3<r1) || (r3<r2)) {
    cpl_msg_error(modName, "wrong radius values: %g %g %g", r1, r2, r3);
    return NULL;
  }
  
  
  /* Allocate output list */
  
  list_out = newPixel(npix);
  
  
  for (pixnum=0; pixnum<npix; pixnum++) {
    
    /* First, determine the background value for each position */
    
    countpix = 0;
    background = 0.00;
    est_x = (double)(p_list[pixnum].x - 1);
    est_y = (double)(p_list[pixnum].y - 1);
    
    for (j=0; j<ima_in->ylen; j++) {
       for (i=0; i<ima_in->xlen; i++) {
         sq_rad = ((double)i - est_x) * ((double)i - est_x) +
           ((double)j - est_y) * ((double)j - est_y);
         if ((sq_rad >= sq_r2) && (sq_rad <= sq_r3)) {
           countpix ++;
           background += (double)ima_in->data[i+j*ima_in->xlen];
         }
       }
    }
    
    if (countpix<=0) {
      cpl_msg_error(modName, "No pixels found in background region");
      return NULL;
     } 
    else {
      background /= (double)countpix;
    }
    
    /*  Now determine the barycenter within the circle centered on the
        pixel, of radius r1 */
    
    corr_x = corr_y = 0.00;
    sum_weights = 0.00;
    
    for (j=0; j<ima_in->ylen; j++) {
      for (i=0; i<ima_in->xlen; i++) {
        sq_rad = ((double)i - est_x) * ((double)i - est_x) +
          ((double)j - est_y) * ((double)j - est_y);
         if (sq_rad <= sq_r1) {
           curpix = (double)ima_in->data[i+j*ima_in->xlen] - background;
           sum_weights += curpix; 
           corr_x += (double)i * curpix; 
           corr_y += (double)j * curpix; 
         }
      }
    }
    
    if (fabs(sum_weights)>1e-10) {
       corr_x /= sum_weights;
       corr_y /= sum_weights;
       list_out[pixnum].x = corr_x + 1;
       list_out[pixnum].y = corr_y + 1;
    } 
    else {
      cpl_msg_error(modName, 
      "Cannot compute barycenter: weighting sum is too small");
      return NULL;
    }
  }
  return(list_out);
}
/**
 * @memo
 *   Refine a peak center by computing the baricenter.
 *
 * @return the peak corrected positions or NULL if the peak is not found.
 *
 * @param ima_in   input image
 * @param in_pixel expected position of the peak
 * @param r1       defines a square around expected position where to look
 *                 for the peak 
 *
 * @doc 
 *  Refine a single peak centers, by computing a barycenter. Define a 
 *  squared region of the image of dimension r1 around expected position
 *  and compute baricenter (for baricenter determination see doc 
 *  of "findPeak2D")
 *   
 * @author P. Sartoretti and C.Izzo
 *                      
 */

VimosPixel *finePositionSimple(VimosImage *ima_in, VimosPixel *in_pixel, 
                          double r1)
{
  VimosPixel *corr_pixel;
  int         x_left=0, x_right=0, y_low=0, y_high=0;
  int         xsize, ysize;
  double      est_x, est_y;
  double      corr_y;
  double      inx, iny;
  float       posx, posy;
  float      *sample;
  char        modName[] = "finePositionSimple";

  
  /* Check input */

  if (ima_in == NULL) {
    cpl_msg_error(modName, "Input NULL image");
    return NULL;
  }
  
  if (in_pixel == (VimosPixel*) NULL) {
    cpl_msg_error(modName, "No pixel in list: cannot refine centroid positions");
    return NULL;
  } 
  inx = in_pixel->x;
  iny = in_pixel->y;
  if ((inx < r1) || (iny < r1) || 
      (inx > (ima_in->xlen - r1)) || 
      (iny > (ima_in->ylen - r1))) {
 /* cpl_msg_warning(modName, "Pixel %f, %f is out/too close to image borders\n",
                  inx, iny); */
    return NULL;
  }

  if (r1<1.0) {
    cpl_msg_error(modName, "Wrong radius values: %g", r1);
    return NULL;
  }
        
 /* Allocate output pixel */

  corr_pixel = newPixel(1);

  est_x = inx;
  est_y = iny;

/*  Determine the barycenter within the square with dimension r1 centered 
    on the pixel */

  corr_y = 0.00;
  x_left  = MAX(floor(est_x)-r1, 0.);
  x_right = MIN(ceil(est_x)+r1, ima_in->xlen);
  y_low   = MAX(floor(est_y)-r1, 0.);
  y_high  = MIN(ceil(est_y)+r1, ima_in->ylen);

  xsize   = x_right-x_left;
  ysize   = y_high-y_low;

  sample = extractFloatImage(ima_in->data,
           ima_in->xlen, ima_in->ylen, x_left, y_low, xsize, ysize);

  if (VM_TRUE == findPeak2D(sample, xsize, ysize, &posx, &posy, 3)) {
    corr_pixel->x = posx + x_left;
    corr_pixel->y = posy + y_low;
  }
  else {
    cpl_msg_warning(modName, 
    "Cannot compute baricenter around input pixel %f, %f", inx, iny);
    return NULL;
  }
     
  cpl_free(sample);

  return(corr_pixel);
}

/**
 * @memo
 *   Shift an image in X and/or Y
 *
 * @return Pointer to the resulting image
 *
 * @param imageIn      pointer to input image
 * @param xshift       shift in X coordinate
 * @param yshift       shift in Y coordinate
 * @param outVal       substitution value
 *
 * @doc 
 *   Image shifting. The output is 1 newly allocated image. 
 *
 * @author P. Franzetti 
 */

VimosImage *imageShift(VimosImage *imageIn, float xshift, float yshift, int outXlen, int outYlen, float outVal)
{
  char             modName[] = "imageShift";

  double           xorig,yorig;
  double           cur;    
  double           neighbors[16];
  double           rsc[8];
  double           sumrs;
  double          *kernel; 
   
  int              i,j,k;
  int              inpXlen;     
  int              tabx, taby;   

  VimosImage      *imageOut;

  VimosLong32     *leaps;  

  VimosUlong32     px,py;
  VimosUlong32     pos; 

/* Error handling */
  
  if (imageIn == NULL) {
    cpl_msg_debug(modName, "NULL input image");
    return(NULL);
  }

  /* Setup arrays for the interpolation*/
  inpXlen = imageIn->xlen;
  
  leaps = (VimosLong32 *) cpl_malloc(16*sizeof(VimosLong32));
  
  /* Check if space was allocated */
  if (leaps == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(VM_FALSE);
  }
  
  /* Setup the interpolation kernel and the vector with the leaps */
  if (!setupInterpolation(&kernel, &leaps, inpXlen)) {
    cpl_msg_error(modName, "Function setupInterpolation failure");
    return(VM_FALSE);
  }

  
  /* Output image allocation */

  imageOut = newImageAndAlloc(outXlen,outYlen);


  /* Do the operation */

  for (j=0; j<outXlen; j++){
    for (i=0; i<outYlen; i++){
      xorig=j-xshift;
      yorig=i-yshift;

      /* ------------ Code from specExtract2D -----------------*/


      /* Which is the closest integer  neighbour?    */
      px = (VimosUlong32) xorig;
      py = (VimosUlong32) yorig;
      
      /* if too close to border: set to outVal */
      if ((px < 1) ||
	  (px > (VimosUlong32) (imageIn->xlen-3)) ||
	  (py < 1) ||
	  (py > (VimosUlong32) (imageIn->ylen-3))) {
	imageOut->data[j+i*outXlen] = outVal;
      } else {
	/* Now feed the positions for the closest 16 neighbours  */
	pos = px + py * inpXlen;
	for (k = 0; k < 16; k++) {
	  neighbors[k] =  
	    (double)(imageIn->data[(VimosUlong32)((VimosLong32)pos
						  +leaps[k])]);
	}
	
	/* Which tabulated value index shall we use for the kernel? */
	tabx = (xorig - (double)px) * (double)(TABSPERPIX); 
	taby = (yorig - (double)py) * (double)(TABSPERPIX); 
	
	/* get resampling coefficients from kernel array  */
	/* rsc[0..3] in x, rsc[4..7] in y   */
	
	rsc[0] = kernel[TABSPERPIX + tabx];
	rsc[1] = kernel[tabx];
	rsc[2] = kernel[TABSPERPIX - tabx];
	rsc[3] = kernel[2 * TABSPERPIX - tabx];
	rsc[4] = kernel[TABSPERPIX + taby];
	rsc[5] = kernel[taby];
	rsc[6] = kernel[TABSPERPIX - taby];
	rsc[7] = kernel[2 * TABSPERPIX - taby];
	
	/* compute normalization */
	sumrs = (rsc[0]+rsc[1]+rsc[2]+rsc[3]) *
	  (rsc[4]+rsc[5]+rsc[6]+rsc[7]);
	
	
	/* Compute interpolated pixel */
	cur = rsc[4] * (rsc[0]*neighbors[0]   +
			rsc[1]*neighbors[1]   +
			rsc[2]*neighbors[2]   +
			rsc[3]*neighbors[3])  +
	  rsc[5] * (rsc[0]*neighbors[4]   +
		    rsc[1]*neighbors[5]   +
		    rsc[2]*neighbors[6]   +
		    rsc[3]*neighbors[7])  +
	  rsc[6] * (rsc[0]*neighbors[8]   +
		    rsc[1]*neighbors[9]   +
		    rsc[2]*neighbors[10]  +
		    rsc[3]*neighbors[11]) +
	  rsc[7] * (rsc[0]*neighbors[12]  +
		    rsc[1]*neighbors[13]  +
		    rsc[2]*neighbors[14]  +
		    rsc[3]*neighbors[15]); 
	
	/* write output value */
	imageOut->data[j+i*outXlen] = (float)(cur/sumrs);
      }
    }
  }
  return(imageOut);
}

VimosImage *openOldFitsFile(const char *fileName, int imaFlag, int readFlag)
{
  const char modName[] = "openOldFitsFile";
  int status,  nfound, anynull;
  float nullval;
  long naxes[2]={1,1};
  long fpixel, npixels;
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */

  VimosImage *tImage=NULL;

  status = 0;

  if (readFlag == 0) {
    if (ffopen( &fptr, fileName, READONLY, &status)) {
         cpl_msg_error(modName, "ffopen returned error %d)", status);
	 return(NULL);
    }
  }
  if (readFlag == 1) {
    if (ffopen( &fptr, fileName, READWRITE, &status)) {
         cpl_msg_error(modName, "ffopen returned error %d", status);
	 return(NULL);
    }
  }

  if(imaFlag == 0)
  {
    tImage = newImage(0,0,NULL);
    if (tImage == NULL) {
      cpl_msg_error(modName, "The function newImage returned a NULL");
      return(NULL);
    }
    tImage->fptr = fptr;
  }

  if(imaFlag == 1) 
  {
    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status)) {
       cpl_msg_error(modName, "fits_read_keys_lng returned error %d", status);
       return(NULL);
    }      

    npixels  = naxes[0] * naxes[1];     /* number of pixels in the image */
    tImage = newImageAndAlloc(naxes[0],naxes[1]);
    /* check if space was allocated */
    if (tImage == NULL) {
      cpl_msg_error(modName, "Allocation Error");
      return(NULL);
    }

    fpixel   = 1;
    nullval  = 0;           /* don't check for null values in the image */

    if (fits_read_img(fptr, TFLOAT, fpixel, npixels, &nullval,
		      (tImage->data), &anynull, &status)) {
      cpl_msg_error(modName, "fits_read_img returned error %d", status);      
      deleteImage(tImage);
      return(NULL);
    }
  
    tImage->xlen = naxes[0];
    tImage->ylen = naxes[1];
    tImage->fptr = fptr;

    if (!readDescsFromFitsImage(&(tImage->descs), tImage)) {
      cpl_msg_error(modName, "readDescsFromFitsImage returned an error");
      return(NULL);
    }
  }

  return(tImage);
}

VimosBool  openNewFitsImage(char *imageName, VimosImage *image)
{
  const char modName[] = "openNewFitsImage";
  int        naxis;
  long       npix[2];
  int        status, bitpix;

  if (image->xlen == 0 && image->ylen == 0) {  /* Empty primary array */
    naxis = 0;
  }
  else if (image->xlen == 0 || image->ylen == 0) {
    cpl_msg_error(modName, "Invalid image sizes.");
    return VM_FALSE;
  }
  else {
    naxis = 2;
  }

  npix[0] = image->xlen;
  npix[1] = image->ylen;

  bitpix = -32;
  
 /*
  *  If a file with the same name already exists, overwrite it
  *  (i.e., delete it first).
  */
  status = 0;
  if (!access(imageName, F_OK)) {
    unlink(imageName);
  }
  status = 0;
  if (fits_create_file(&image->fptr, imageName, &status)) {
    cpl_msg_error(modName, "fits_create_file returned error %d", status);
    return VM_FALSE;
  }
  if (fits_create_img(image->fptr,  bitpix, naxis, npix, &status)) {
    cpl_msg_error(modName, "fits_create_img returned error %d", status);
    return VM_FALSE;
  } 
  return VM_TRUE;
}

VimosBool closeFitsImage(VimosImage *image, int flag) 
{
  const char modName[] = "closeFitsImage";
  int status;
  long fpixel, npixels;

  status = 0;

  if (image == NULL) return VM_FALSE;

  npixels  = image->xlen * image->ylen;         
  fpixel   = 1;

  if (flag != 0 ){
    if (fits_write_img(image->fptr, TFLOAT, fpixel, npixels, image->data,
		       &status)) {
      cpl_msg_error(modName, "fits_write_img returned error %d", status);
      return(VM_FALSE); 
    }     
  }

  status = 0;
  if (fits_close_file(image->fptr, &status)) {
    cpl_msg_error(modName, "fits_close_file returned error %d", status);
    return(VM_FALSE);
  }

  return(VM_TRUE);
}

VimosBool appendNewFitsImage(VimosImage *image,fitsfile *fptr,char extname[])
{
  const char modName[] = "appendNewFitsImage";
  int naxis;
  long npix[2],fpixel,npixels;
  int status, bitpix;

  npix[0] = image->xlen;
  npix[1] = image->ylen;

  naxis = 2;
  bitpix = -32;
  npixels  = image->xlen * image->ylen;         
  fpixel   = 1;
  image->fptr = fptr;

  status = 0;

  /* if extension is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, IMAGE_HDU, extname, 0, &status)) {
    if (fits_delete_hdu(fptr, NULL, &status)) {
      cpl_msg_error(modName, "fits_delete_hdu returned error %d", status);
      return (VM_FALSE);
    }
  } else {
    status = 0;
  }
 
  if (fits_create_img(fptr, bitpix, naxis, npix, &status)) {
    cpl_msg_error(modName, "fits_create_img returned error %d", status);
    return(VM_FALSE);
  }
  /*  if (fits_movnam_hdu(fptr, IMAGE_HDU, extname, 0, &status) ) {
      return (VM_FALSE);
      } */

  if (fits_write_img(fptr, TFLOAT, fpixel, npixels, image->data, &status)) {
    cpl_msg_error(modName, "fits_write_img returned error %d", status);
    return(VM_FALSE);      
  }

  if (!writeDescsToFitsImage(image->descs, image)) {
    cpl_msg_error(modName, "writeDescsToFitsImage returned an error");
    return(VM_FALSE);
  }

  if (fits_update_key_str(fptr, "EXTNAME", extname, "", &status)) {
    cpl_msg_error(modName, "fits_update_key_str returned error %d", status);
    return(VM_FALSE);
  } else {
    return(VM_TRUE);
  }
    
}


VimosImage *openFitsImageExtension(fitsfile *fptr, char extname[])
{
  const char modName[] = "openFitsImageExtension";
  int status,  nfound, anynull;
  float nullval;
  long naxes[2], fpixel, npixels;

  VimosImage *tImage;

  status = 0;

  if (fits_movnam_hdu(fptr, IMAGE_HDU, extname, 0, &status)) {
    cpl_msg_error(modName, "fits_movnam_hdu returned error %d", status);
    return(NULL);
  }

  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status)) {
    cpl_msg_error(modName, "fits_read_keys_lng returned error %d", status);
    return(NULL);
  }

  npixels  = naxes[0] * naxes[1];     /* number of pixels in the image */

  tImage = newImageAndAlloc(naxes[0],naxes[1]);
  /* check if space was allocated */
  if (tImage == NULL) {
    cpl_msg_error(modName, "newImageAndAlloc has returned NULL");
    return(NULL);
  }

  fpixel   = 1;
  nullval  = 0;           /* don't check for null values in the image */

  if (fits_read_img(fptr, TFLOAT, fpixel, npixels, &nullval,
		     (tImage->data), &anynull, &status)) {
    cpl_msg_error(modName, "fits_read_img returned error %d", status);
    deleteImage(tImage);
    return(NULL);
  }
  
  tImage->xlen = naxes[0];
  tImage->ylen = naxes[1];
  tImage->fptr = fptr;

  if (!readDescsFromFitsImage(&(tImage->descs), tImage)){
    cpl_msg_error(modName, "readDescsFromFitsImage returned an error");
    return(NULL);
  }

  return(tImage);
}

VimosBool copyFitsExtensions(VimosImage *outImage,VimosImage *inImage)
{
  const char modName[] = "copyFitsExtensions";
  int i;
  int extNum;
  int status = 0;

  if (fits_get_num_hdus(inImage->fptr, &extNum, &status)) {
    cpl_msg_error(modName, "fits_get_num_hdus returned error %d", status);
    return(VM_FALSE);
  }

  if (extNum > 1)
  {
    for (i = 2; i <= extNum; i++)
    {
      if (fits_movabs_hdu(inImage->fptr, i, NULL, &status)) {
        cpl_msg_error(modName, "fits_movabs_hdu returned error %d", status);
        return(VM_FALSE);
      }
      if (fits_copy_hdu(inImage->fptr, outImage->fptr, 0, &status)) {
        cpl_msg_error(modName, "fits_copy_hdu returned error %d", status);
        return(VM_FALSE);
      }
    }
  }

  return (VM_TRUE);
}
/** 
 * @memo 
 *   Get the list of 
 *
 * @return list of extension names. 
 *
 * @param fptr         pointer to fitsfile
 * @param extNum         pointer to extension number 
 *
 * @doc 
 *   Create a list of extension names of a FITS file
 *
 * @author C. Izzo
 */
char **getFitsFileExtensionsNames(fitsfile *fptr, int *extNum)
{
  int    i;
  int    status = 0;
  char **extNames = NULL;

  *extNum = 0;

  fits_get_num_hdus(fptr, extNum, &status);
  extNames = (char **) cpl_malloc(*extNum*sizeof(char *));

  if (*extNum > 1) {
    for (i=2; i <= *extNum; i++) {
      extNames[i-2] = (char *) cpl_malloc(FLEN_VALUE*sizeof(char));
      fits_movabs_hdu(fptr, i, NULL, &status);
      fits_read_key_str(fptr,"EXTNAME",extNames[i-2],NULL,&status);
      if (status) {
        strcpy(extNames[i-2],"Not found");
        status = 0;
      }
    }
    (*extNum)--;
  }
  return extNames;
}
/** 
 * @memo 
 *   Write a VimosImage into a disk FITS file
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param filename         name of the fitsfile
 * @param image            Vimos Image to be written as a FITS file 
 * @param category         Category of image to write
 *
 * @doc 
 *   Create a disk image FITS file from structure VimosImage
 *
 * @author C. Izzo
 */

VimosBool createFitsImage(char *filename, VimosImage *image, 
                          const char *category)
{
  char modName[] = "createFitsImage";

  if (openNewFitsImage(filename, image)) {

   /*
    *  Make sure that the image is float, if the header keywords
    *  were inherited by images of other types.
    */

    if (writeIntDescriptor(&(image->descs), "BITPIX", -32, 
                           "No. of bits per pixel") == VM_TRUE) {
      if (writeIntDescriptor(&(image->descs), "BITPIX", -32, 
                             "No. of bits per pixel") == VM_TRUE) {
        removeDescriptor(&(image->descs), "BSCALE");
        removeDescriptor(&(image->descs), "BZERO");

       /*
        *  Make sure that the image axes have proper dimensions
        */

        if (writeIntDescriptor(&(image->descs), "NAXIS1", 
                               image->xlen, "Pixel in X") == VM_TRUE) {
          if (writeIntDescriptor(&(image->descs), "NAXIS2", 
                                 image->ylen, "Pixel in Y") == VM_TRUE) {

            if (writeDescsToFitsImage(image->descs, image) == VM_TRUE) {
              if (closeFitsImage(image, 1) == VM_TRUE) {
                cpl_msg_debug(modName, "Image %s (%s) created", 
                            filename, category);
                return VM_TRUE;
              }
            }
          }
        }
      }
    }
  }

  return VM_FALSE;
}


int mapTable(VimosImage *image, double start, double step,
             VimosTable *table, char *xname, char *yname)
{
  const char modName[] = "mapTable";

  float  *xdata   = tblGetFloatData(table, xname);
  float  *ydata   = tblGetFloatData(table, yname);
  int     xlength = tblGetSize(table, xname);
  int     length  = image->xlen;
  float   xzero, pos;
  int     i, j, n;


  if (image->ylen != 1) {
    cpl_msg_error(modName, "Input image Y size should be 1");
    return EXIT_FAILURE;
  }

  /*
   * Initialization at 0.0 - this value is left on non-overlapping
   * portions.
   */

  for (i = 0; i < length; i++)
    image->data[i] = 0.0;

  n = 0;
  xzero = xdata[n];

  for (i = 0; i < length; i++) {
    pos = start + step * i;
    if (pos < xzero)
      continue;
    for (j = n; j < xlength; j++) {
      if (xdata[j] > pos) {
        n = j;
        image->data[i] = ydata[j-1]
                       + (ydata[j] - ydata[j-1])
                       * (pos - xdata[j-1]) / (xdata[j] - xdata[j-1]);
        break;
      }
    }
  }

  return EXIT_SUCCESS;

}


int mapTableDouble(VimosImage *image, double start, double step,
                   VimosTable *table, char *xname, char *yname)
{
  const char modName[] = "mapTableDouble";

  double *xdata   = tblGetDoubleData(table, xname);
  double *ydata   = tblGetDoubleData(table, yname);
  int     xlength = tblGetSize(table, xname);
  int     length  = image->xlen;
  float   xzero, pos;
  int     i, j, n;


  if (image->ylen != 1) {
    cpl_msg_error(modName, "Input image Y size should be 1");
    return EXIT_FAILURE;
  }

  /*
   * Initialization at 0.0 - this value is left on non-overlapping
   * portions.
   */

  for (i = 0; i < length; i++)
    image->data[i] = 0.0;

  n = 0;
  xzero = xdata[n];

  for (i = 0; i < length; i++) {
    pos = start + step * i;
    if (pos < xzero)
      continue;
    for (j = n; j < xlength; j++) {
      if (xdata[j] > pos) {
        n = j;
        image->data[i] = ydata[j-1]
                       + (ydata[j] - ydata[j-1])
                       * (pos - xdata[j-1]) / (xdata[j] - xdata[j-1]);
        break;
      }
    }
  }
  
  return EXIT_SUCCESS;
  
}


int polySmooth(VimosImage *image, int ord, int hw)
{

  const char modName[] = "polySmooth";

  VimosDpoint *list     = newDpoint(image->xlen);
  float       *smoo     = malloc(image->xlen * sizeof(float));
  float        mean     = imageMean(image);
  float        dev      = imageAverageDeviation(image, mean);
  double      *c        = NULL;
  double       value;
  double       factor;
  int          box      = 2 * hw + 1;
  float       *row      = malloc(box * sizeof(float));
  int          i, j;


  if (image->ylen > 1)
  {
    free(smoo);
    free(row);
    return EXIT_FAILURE;
  }

  if (box > image->xlen)
  {
    free(smoo);
    free(row);
    return EXIT_FAILURE;
  }

  /*
   * Copy first hw and last hw items
   */

  for (i = 0; i < hw; i++)
    smoo[i] = image->data[i];

  for (i = image->xlen - hw; i < image->xlen; i++)
    smoo[i] = image->data[i];

  /*
   * Median filtering
   */

  for (i = hw; i < image->xlen - hw; i++) {
    for (j = -hw; j <= hw; j++)
      row[j + hw] = image->data[i + j];
    smoo[i] = median(row, box);
  }

  free(row);

  for (i = 0; i < image->xlen; i++)
    image->data[i] = smoo[i];

  free(smoo);

  list = newDpoint(image->xlen);

  for (i = 0; i < image->xlen; i++) {
    list[i].x = i - image->xlen;
    list[i].y = (image->data[i] - mean) / dev;
  }

  c = fit1DPoly(ord, list, image->xlen, NULL);

  if (c) {
    for (i = 0; i < image->xlen; i++) {
      value = 0.0;
      factor = 1.0;
      for (j = 0; j <= ord; j++) {
        value += c[j] * factor;
        factor *= i - image->xlen;
      }
      image->data[i] = value * dev + mean;
    }
  }
  else
    cpl_msg_warning(modName, "No smoothing possible...");


  return EXIT_SUCCESS;

}

/**@}*/

