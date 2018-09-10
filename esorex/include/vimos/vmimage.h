/* $Id: vmimage.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_IMAGE_H
#define VM_IMAGE_H

#include <fitsio.h>

#include <pilmacros.h>

#include <vmtypes.h>
#include <vmtable.h>
#include <vmmath.h>


PIL_BEGIN_DECLS

/* Define the minimum number of frame for different stacking methods */

#define MIN_FRAMES_KSIGMA 2
#define MIN_FRAMES_REJECT 2
#define MIN_FRAMES_MEDIAN 3
#define MIN_FRAMES_AVERAGE 2
#define MIN_FRAMES_SUM 2

#define MAX_PIXEL_VALUE  ((float)LONG_MAX)


typedef struct _VIMOS_IMAGE_
{
  int xlen;
  int ylen;
  float *data;
  VimosDescriptor *descs;
  fitsfile *fptr;
} VimosImage;

typedef enum _COMB_METHOD_
{
  COMB_UNDEF = 0,
  COMB_AUTO, 
  COMB_KSIGMA,
  COMB_REJECT, 
  COMB_MEDIAN, 
  COMB_AVERAGE, 
  COMB_SUM   
} CombMethod;

typedef enum _FILTER_METHOD_
{
  FILTER_UNDEF = 0,
  FILTER_AUTO, 
  FILTER_MEDIAN,
  FILTER_AVERAGE, 
  FILTER_GAUSS
} FilterMethod;

/* enum type for function collapse2Dto1D.c */

typedef enum _IMAGE_DIRECTION_
{
  ROW,
  COLUMN
} Direction;

typedef struct _COMB_PARAMETERS_
{
  double kSigmaLow;
  double kSigmaHigh;
  int minRejection;
  int maxRejection;
} CombParameters;


VimosImage  *newImage(int xlen, int ylen, float *data);
void deleteImage(VimosImage *image);
VimosImage  *newImageAndAlloc(int xlen, int ylen);
void deleteImageAndAlloc(VimosImage *);
VimosImage *duplicateImage(VimosImage *);


/* Arithmetic operations between 2 images. The 2 images must have the 
   same dimensions. The output is 1 newly allocated image */

VimosImage *imageArith(VimosImage *ima1, VimosImage *ima2, 
                        VimosOperator optype);

/* Arithmetic operations between an image and a constant.
   The output is 1 newly allocated image */

VimosImage *constArith(VimosImage *ima_in, double constant, 
                        VimosOperator optype);


/*  Arithmetic operations between 2 images. Stores the result on the 
    first frame. Used in frame combination. */ 
   
int imageArithLocal(VimosImage *ima1, VimosImage *ima2, VimosOperator optype);


/* Arithmetic operations between an Image and a constant. Store
   result in the image */

int constArithLocal(VimosImage *ima_in, double constant, VimosOperator optype);

/* Apply an optimized median filter of sizes filtsizex, filtsixey to an image.
   The median kernel is gradually reduced to (filtsize/2+1)*(filtsize/2+1)
   at borders. The output is a newly allocated image.
   The algorithm is taken from Eclipse (N. Devilla), extended to allow
   the exclusion of the central pixel from median computation */

VimosImage *VmFrMedFil(VimosImage *ima_in, int filtsizex, int filtsizey, int);

/* Apply an optimized average filter of sizes filtsizex, filtsixey to an image.
   The average box is gradually reduced to (filtsize/2+1)*(filtsize/2+1)
   at borders. The output is a newly allocated image. The algorithm is 
   taken from VmFrMedFil(), where the median computation is replaced by
   the average. */

VimosImage *VmFrAveFil(VimosImage *ima_in, int filtsizex, int filtsizey, int);

VimosImage *VmFrFilter(VimosImage *image, 
                                     int xSize, int ySize, FilterMethod, int);

/* Image normalization */

VimosImage *VmImNorm( VimosImage *ima_in, Method meth);
 
/*
 * Combine a list of frame according to the method. Possible methods are
 * COMB_KSIGMA, COMB_MEDIAN, COMB_AVERAGE, COMB_SUM, COMB_REJECT, and AUTO. 
 * AUTO automatically selects the combination method based on the number 
 * of frames
 */

VimosImage *frComb(VimosImage **ima_list, int num, CombMethod combMethod,
            CombParameters *, int flag);
VimosImage *frComb32000(VimosImage **ima_list, int num, CombMethod combMethod,
            CombParameters *, int flag);


/* Combine a list of frames: the output is the sum of all frames */

VimosImage *frCombSum(VimosImage **ima_list, int num);


/* Combine a list of frames: the output is the average of all frames */

VimosImage *frCombAverage(VimosImage **ima_list,int num);
VimosImage *frCombAverage32000(VimosImage **ima_list,int num);


/* Combine a list of frames: the output is the median */

VimosImage *frCombMedian(VimosImage **ima_list, int num, int flag);


/* Combine a list of frame with rejection of pixels outside a
  minimum/maximum treshold.
  Input a list of images, the percentages of pixel to minimum-reject
  (i.e. to set below the minimum treshold), and the percentage
  of pixels to maximum-reject. The percentage range is 0-1  */

VimosImage *OLDfrCombMinMaxReject
       (VimosImage **ima_list, double minrej, double maxrej, int num);

/* Same as above, but with absloute number of images to exclude */

VimosImage *frCombMinMaxReject
       (VimosImage **ima_list, int minrej, int maxrej, int num);
VimosImage *frCombMinMaxReject32000
       (VimosImage **ima_list, int minrej, int maxrej, int num);


/* K-sigma-clipping: Combine frame with rejection of pixels outside
   a k*sigma treshold. The treshold is calculated on the median
   value, and not on the average to limit the influence of bad pixels.
   Input a list of images; klow to set the low treshold (=median-klow*sigma),
   and khigh to set the high treshold
   (=median+khigh*sigma)*/

VimosImage *frCombKSigma
       (VimosImage **ima_list, double klow, double khigh, int num);
VimosImage *frCombKSigma32000
       (VimosImage **ima_list, double klow, double khigh, int num);

/* compute the mean value of an image */

float imageMean(VimosImage *ima_in);


/* compute the median value of an image */

float imageMedian(VimosImage *ima_in);

float imageAverageDeviation(VimosImage *, float);

/* compute image sigma. */

float imageSigma(VimosImage *ima_in);


/* compute image sigma using median_value as mean estimate. Useful for
   kappa/sigma rejection */

float imageMedSigma(VimosImage *ima_in);


/* returns the maximum pixel value of an image */

float imageMaximum(VimosImage  *ima_in);


/* returns the minimum pixel value of an image */

float imageMinimum(VimosImage *ima_in);


/* compute an image histogram. Input one image and the number of bins for
   the histogram. Output a list of [nbins] dpoints. x contains the central
   pixelvalue for each bin, y contains the number of pixels in this bin.
   The algorithm is taken from Eclipse. */

VimosDpoint *imageHistogram(VimosImage *ima_in, unsigned int nbins);


/* compute the mode of an image (=histo peak position).
   It depends on the number of bins (binzise)!
   Here take the bin size = 1 data unit (i.e. bin number =(int) max - min) */

float imageMode (VimosImage *ima_in);


/* compute the sum of all pixelvalues in a rectangle. The rectangle is given
   by the lower left corner coordinates (x,y) and the number of pixels
   nx and ny */

float sumPixelsInImage(VimosImage *ima_in, int x, int y, int nx, int ny);


/* collapse a rectangle inside an image into one column. Compute the sum
   of the rows or of the columns of a rectangle. Input an Image,
  the lower left corner coordinates (x,y) of the rectangle; the
  number of pixels (nx, ny) from (x,y), and the collapsing direction:
  collapse ROW: sum of rows and collapse COLUMN: sum of columns. */

float *collapse2Dto1D(VimosImage *ima_in, int x, int y, int nx, int ny,
                       Direction collapse);

/* Extracted from Eclipse and adapted to Vimos Imaging
   fine_position_centers()
   Input: In an image, a list of pixel positions, a number of pixels
   corresponding to that list, and 3 doubles
   Out: void
   Job: refine the peak centers, by computing a barycenter
   Notice:      the 3 doubles (r1, r2, r3) define 3 circle radiuses:
   r2 and r3 define a ring centered on each pixel in the list,
   around which the background will be estimated.
   r1 defines a disk in which the barycenter will be computed,
   i.e. the centroid pixel, weighted by pixel values, from
   which the background value has been subtracted:
   x_center = (1/sum(pixel[i]-bg))*(sum(x[i] * (pixel[i]-bg)))
   y_center = (1/sum(pixel[i]-bg))*(sum(y[i] * (pixel[i]-bg))) */

VimosPixel *finePosition(VimosImage *ima_in, VimosPixel *p_list, int npix,
        double r1, double r2, double r3);

/*  ... */

VimosPixel *finePositionSimple(VimosImage *ima_in, VimosPixel *in_pixel,
                        double r1);

/*   Shift an image in X and/or Y  */

VimosImage *imageShift(VimosImage *imageIn, float xshift, float yshift,
                       int outXlen, int outYlen,  float outVal);

VimosBool readDescsFromFitsImage(VimosDescriptor **desc, VimosImage *image);
VimosBool writeDescsToFitsImage(VimosDescriptor *desc, VimosImage *image);

/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosImage *openOldFitsFile(char *fileName, int imaFlag, int readFlag)

   Description
   Open an existing FITS file, either in READONLY or in READWRITE mode,
   depending on the value of the readFlag parameter. The image data from the
   primary array are read only if the imaFlag parameter equals 1.

   Input:
   char *fileName
   Name of the FITS file containing the image

   int imaFlag
   Flag to identify wheather the image in the primary array should be read
   (flag = 0 means DON'T READ IMAGE; flag = 1 means READ IMAGE)

   int readFlag
   Flag to identify the opening mode for the image (flag = 0 means READONLY;
   flag = 1 means READWRITE)

   Return Value:
   Pointer to the VimosImage structure where the image stuff is put

   Updates:
   15 Dec 99: Created (BG)
   11 Feb 00: Added the flag parameter to handle different opening modes (MS)
   24 May 00: Modified from openOldFitsImage into openOldFitsFile, to allow
              the opening of a file without reading in the primary array 
	      data (MS)

-----------------------------------------------------------------------------
*/
VimosImage *openOldFitsFile(const char *fileName, int imaFlag, int readFlag);

VimosBool loadFitsHeader(VimosImage *);
VimosBool loadFitsData(VimosImage *);

VimosBool openNewFitsImage(char *imageName, VimosImage *image);

VimosBool closeFitsImage(VimosImage *image, int flag);

VimosBool appendNewFitsImage(VimosImage *image,fitsfile *fptr,char extname[]);

VimosImage *openFitsImageExtension(fitsfile *fptr, char extname[]);

VimosBool copyFitsExtensions(VimosImage *outImage,VimosImage *inImage);

/*
 * Return the list of extensions names
 */

char **getFitsFileExtensionsNames(fitsfile *, int *);

/*
 * Create a disk image FITS file from structure VimosImage
 */

VimosBool createFitsImage(char *filename, VimosImage *, const char *category);

int mapTable(VimosImage *, double, double, VimosTable *, char *, char *);
int mapTableDouble(VimosImage *, double, double, VimosTable *, char *, char *);
int polySmooth(VimosImage *, int, int);

PIL_END_DECLS

#endif /* VM_IMAGE_H */
