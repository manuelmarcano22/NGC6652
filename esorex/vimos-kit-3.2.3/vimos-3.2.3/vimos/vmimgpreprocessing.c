/* $Id: vmimgpreprocessing.c,v 1.3 2013-03-25 11:43:04 cgarcia Exp $
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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>

#include <pilmemory.h>
#include <piltranslator.h>
#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmccdtable.h"
#include "vmdetector.h"
#include "vmimgutils.h"
#include "vmimgpreprocessing.h"
#include "cpl.h"



/* Set the polynomial degree to fit overscan/prescan regions */

#define DEGREE            (2)

/* Set default parameters for cosmic rays cleaning */

#define DEFAULT_THRESHOLD (4.)
#define DEFAULT_RATIO     (2.)
#define SMOOTHING_BOX     (3)


/**
 * @name vmimgpreprocessing Generic Image Preprocessing
 *
 * The module collects high/medium level functions for common operations
 * on images.
 */

/**@{*/

/**
 * @memo 
 *   Master Bias subtraction from image.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE;
 * 
 * @param image  Image to be Master Bias subtracted.
 * @param mbias  Master Bias image.
 * @param method Method used for bias subtraction: BIAS_UNDEF = don't trim 
 *               overscans after master bias subtraction; BIAS_MASTER = trim 
 *               overscans away; BIAS_ZMASTER = apply overscan correction 
 *               before trimming overscans.
 *
 * @doc
 *   Subtract Master Bias from an image, and if method is BIAS_ZMASTER 
 *   then correct the result with the image overscan/prescan region(s)
 *   possible residual trends. Trim away overscan/prescan regions if the
 *   chosen method is not BIAS_UNDEF.
 *
 * @author C. Izzo
 */    

int VmSubBias(VimosImage *image, VimosImage *mbias, BiasMethod method)
{
  char         modName[] = "VmSubBias";

  VimosImage  *outImage = NULL;
  VimosPort   *ports, *currport;
  VimosDpoint *avedpoint;
  int          vertical;
  int          nports, i, j;
  int          averaged;
  int          npoints;
  int          px, py, pnx, pny, ox, oy, onx, ony, wx, wy, wnx, wny;
  float       *presc1d = NULL;
  float       *ovesc1d = NULL;
  int          prescSize = 0, ovescSize = 0;
  double      *coef;

  /* Check input */
  
  cpl_msg_debug(modName, "Subtracting Bias");

  if ((image == NULL) || (mbias == NULL)) {
    cpl_msg_error(modName, "Null input");
    return EXIT_FAILURE;
  }

  /* Subtract Masterbias */
  
  if (imageArithLocal(image, mbias, VM_OPER_SUB) == 0) {
    ports = getPorts(image, &nports);
    if (method != BIAS_UNDEF && ports->prScan->nX + ports->ovScan->nX 
                             && ports->prScan->nY + ports->ovScan->nY) {

      /*
       * Overscans are present
       */

      if (method == BIAS_ZMASTER) {

        /*
         * Apply overscan correction before trimming overscans away.
         */

        for (currport = ports; currport; currport = currport->next) {

          /*
           * Keep regions position and size in short named
           * variables for clarity
           */

          px  = currport->prScan->startX;           /* Prescan  */
          py  = currport->prScan->startY;
          pnx = currport->prScan->nX;
          pny = currport->prScan->nY;
          ox  = currport->ovScan->startX;           /* Overscan */
          oy  = currport->ovScan->startY;
          onx = currport->ovScan->nX;
          ony = currport->ovScan->nY;
          wx  = currport->readOutWindow->startX;    /* Port     */
          wy  = currport->readOutWindow->startY;
          wny = currport->readOutWindow->nY;
          wnx = currport->readOutWindow->nX;
          
          vertical = getReadoutDirection(currport);
          
          /*
           * Collapse prescan and/or overscan in one column or row
           * according to the scan direction.
           */

          averaged = 0;
          if (vertical) {
            if (pny > 0) {
              presc1d = collapse2Dto1D(image, px, py, pnx, pny, ROW);
              prescSize = pnx;
              averaged = pny;
            }
            if (ony > 0) {
              ovesc1d = collapse2Dto1D(image, ox, oy, onx, ony, ROW);
              ovescSize = onx;
              averaged += ony;
            }
          }
          else {
            if (pnx > 0) {
              presc1d = collapse2Dto1D(image, px, py, pnx, pny, COLUMN);
              prescSize = pny;
              averaged = pnx;
            }
            if (onx > 0) {
              ovesc1d = collapse2Dto1D(image, ox, oy, onx, ony, COLUMN);
              ovescSize = ony;
              averaged += onx;
            }
          }

          /*
           * Note that "averaged" is the total number of point
           * BOTH in the collapsed prescan AND the collapsed
           * overscan. This simplifies the direct computation
           * of the average bias in overscan regions.
           */
          
          if (presc1d != NULL && ovesc1d != NULL) {
            if (prescSize != ovescSize) {
              cpl_msg_error(modName, "Lengths of prescan and overscan differ!");
              return EXIT_FAILURE;
            }
          }
          
          /*
           * If the prescan is missing, the overscan exists for sure:
           * get the number of points of collapsed regions from one
           * of them, and allocate the appropriate list of points:
           */

          if (presc1d != NULL) npoints = prescSize;
          else                 npoints = ovescSize;
          
          avedpoint = newDpoint(npoints);
          
          for (i = 0; i < npoints; i++) {
            avedpoint[i].x = i;         /* Local coordinate        */
            avedpoint[i].y = 0;         /* This init is necessary! */
          }
          
          if (presc1d != NULL) {
            for (i = 0; i < npoints; i++) { 
              avedpoint[i].y = presc1d[i] / averaged;
            }
            cpl_free(presc1d);
            presc1d = NULL;
          }
          
          if (ovesc1d != NULL) {
            for (i = 0; i < npoints; i++) { 
              avedpoint[i].y += ovesc1d[i] / averaged;
            }
            cpl_free(ovesc1d);
            ovesc1d = NULL;
          }
          
          /*
           * Fit now the averaged data with a polynomial
           */

          coef = fit1DPoly(DEGREE, avedpoint, npoints, NULL);
          
          /*
           * Now recycle avedpoint to contain model, instead of data points
           */

          for (i = 0; i < npoints; i++) {
            avedpoint[i].y = 0;
            for (j = 0; j <= DEGREE; j++){
              avedpoint[i].y += coef[j] * ipow(i, j);
            }
          } 
          
          /*
           * Subtract the model from image according to the scan
           * direction.
           */

          if (vertical) {
            for (j = 0; j < wny; j++) {
              for (i = 0; i < npoints; i++) {
                image->data[(i + wx) + (j + wy) * image->xlen] 
                                                          -= avedpoint[i].y;
              }
            }
          }
          else {
            for (i = 0; i < npoints; i++) {
              for (j = 0; j < wnx; j++) {
                image->data[(j + wx) + (i + wy) * image->xlen] 
                                                          -= avedpoint[i].y;
              }
            }
          }
          deleteDpoint(avedpoint);
        }
      }
      if (EXIT_FAILURE == trimOverscans(image)) {
        cpl_msg_error(modName, "Cannot trim overscans");
        return EXIT_FAILURE;
      }
    } 
    outImage = image;
    deletePortList(ports);
  }
  else {
    cpl_msg_error(modName, "Cannot subtract Master Bias from original image");
    return EXIT_FAILURE;
  } 
  return EXIT_SUCCESS;
}


/**
 * @memo 
 *   Master Dark subtraction.
 *
 * @return EXIT_SUCCESS / EXIT_FAILURE 
 *
 * @param ima_in Image to be dark subtracted.
 * @param mdark  Master Dark.
 *
 * @doc 
 *   Scale master Dark to image exposure time and subtract. 
 *   Overwrite result on input image. 
 * 
 * @author P. Sartoretti, C. Izzo.
 */    

int VmSubDark(VimosImage *ima_in, VimosImage *mdark)
{

  VimosImage  *scaledDark;
  double       time;
  char         modName[] = "VmSubDark";


  /* check input */
  cpl_msg_debug(modName,"subtracting Dark");

  if (mdark == NULL) {
    cpl_msg_error(modName, "Null master Dark");
    return EXIT_FAILURE;
  }
  if (ima_in == NULL) {
    cpl_msg_error(modName,"Null input image\n");
    return EXIT_FAILURE;
  }
  
  /* scale master dark and subtract */
  
  if (readDoubleDescriptor(ima_in->descs, pilTrnGetKeyword("ExposureTime"), 
                           &time, NULL) == VM_TRUE) {
    scaledDark = constArith(mdark, time, VM_OPER_MUL);
    imageArithLocal(ima_in, scaledDark, VM_OPER_SUB);
    deleteImage(scaledDark);
    return EXIT_SUCCESS;
  }
  else 
    return EXIT_FAILURE;
}


/**
 * @memo
 *   Remove cosmic rays from single CCD exposure.
 * 
 * @return EXIT_SUCCESS / EXIT_FAILURE
 *
 * @param image      Pointer to input image to be cleaned. The image
 *                   pixels values are assumed in ADU.
 * @param ccdTable   Optional CCD table, containing a list of bad
 *                   pixels to be corrected simultaneously with the 
 *                   found cosmic rays events: this is to avoid
 *                   correcting cosmic rays using bad pixel values
 *                   when interpolating hit pixels.
 * @param cleanBad   If set to true a ccdTable must be given, and bad 
 *                   pixels are corrected simultaneously with cosmic
 *                   rays events. If set to false and a ccdTable is 
 *                   given, then cosmic rays events are interpolated
 *                   avoiding bad pixels listed in the ccdTable.
 * @param skyLevel   Sky level (in ADU). If set to a negative number,
 *                   it is internally assigned the median level of
 *                   the input image.
 * @param gain       Inverse gain factor (e-/ADU).
 * @param ron        Read-out-noise (ADU).
 * @param threshold  Threshold for cosmic rays detection, given in
 *                   units of noise sigma. Where the median filtered 
 *                   input image differs from the unfiltered original 
 *                   more than this threshold, a cosmic ray candidate 
 *                   is found. If the threshold is negative, a default 
 *                   value of 4 is assigned.
 * @param ratio      Critical ratio for discrimination between objects
 *                   and cosmic rays. If negative, it is assigned the 
 *                   default value 2.
 *
 * @doc
 *   The algorithm used is the same of the MIDAS command FILTER/COSMIC,
 *   with few extensions, including the possibility of correcting the
 *   cosmic rays simultaneously with existing bad pixels regions, or
 *   avoiding bad pixel regions when interpolating the new values 
 *   replacing the cosmic rays events. 
 *
 * @see cleanBadPixels
 *
 * @author C. Izzo
 */    

int VmCosmicClean(VimosImage *image, VimosTable *ccdTable, int cleanBad,
    float skyLevel, float gain, float ron, float threshold, float ratio)
{
  const char   modName[] = "VmCosmicClean";
  /*** VimosImage  *outputMask; ***/  /* This is just a check image */
  VimosImage  *smoothImage;
  VimosTable  *table;
  VimosColumn *xColumn;
  VimosColumn *yColumn;
  VimosColumn *xBad = NULL;
  VimosColumn *yBad = NULL;
  VimosColumn *colLevel = NULL;
  VimosPort   *port;
  float        sigma, sum, value, smoothValue;
  float        fMax;
  int          iMin, iMax, jMin, jMax, iPosMax, jPosMax;
  int          xLen = image->xlen;
  int          yLen = image->ylen;
  int          nPix = xLen * yLen;
  int          first = 1;  /* position of first cosmic ray candidate
                              encountered while scanning the image */
  int          pos, i, j, k, l, ii, jj, iii = 0, jjj = 0;
  int          numCosmic = 0;
  int          numBad = 0;
  int          numPix = 0;
  int          found, foundContiguousCandidate;
  int         *cosmic;
  int          nports;
  int          shiftX = 0;
  int          shiftY = 0;
  int          cleanJustCosmics = ccdTable && !cleanBad;
  int          cosmicFlag = 1;
  int          badPixelFlag = 0;
  int          exitStatus = EXIT_FAILURE;

  cpl_msg_debug(modName, "Removing Cosmic Rays");

  if (image == NULL) {
    cpl_msg_debug(modName, "Missing input image");
    return(exitStatus);
  }

  if (cleanBad && !ccdTable) {
    cpl_msg_error(modName, "Cannot clean bad pixels if no CCD Table is given");
    return(exitStatus);
  }

  /*** outputMask = duplicateImage(image); ***/

 /*
  *  "cosmic" is a flags holder (initialized to zero):
  *
  *           -1 = candidate for cosmic ray
  *            0 = not a cosmic
  *            1 = a cosmic ray
  *            2 = member of current group of contiguous candidates
  *            3 = examined member of current group
  */
  cosmic = (int *) cpl_malloc(nPix * sizeof(int));
  for (i = 0; i < nPix; i++) {
    cosmic[i] = 0;
    /*** outputMask->data[i] = 0.; ***/
  }

  if (skyLevel < 0.)  skyLevel  = imageMedian(image);
  if (threshold < 0.) threshold = DEFAULT_THRESHOLD;
  if (ratio < 0.)     ratio     = DEFAULT_RATIO;

  cpl_msg_debug(modName, 
    "Sky level: %f, Threshold: %f sigma, Discrimination ratio "
    "cosmic/object: %f, Gain = %f e-/ADU, RON = %f ADU\n", 
    skyLevel, threshold, ratio, gain, ron);

  smoothImage = VmFrMedFil(image, SMOOTHING_BOX, SMOOTHING_BOX, 1);

 /*
  *  Loop on images pixels, searching for cosmic rays candidates.
  *  Border pixels are currently excluded (they cannot contain 
  *  candidates), to avoid that the search for groups of contiguous 
  *  pixels would ever go out of image boundaries. In future we may
  *  overcome this limit, adding an appropriate check when contiguous
  *  pixels are searched.
  */

  for (j = 1; j < yLen - 1; j++) {
    for (i = 1; i < xLen - 1; i++) {
      value = image->data[i + j * xLen];
      smoothValue = smoothImage->data[i + j * xLen];
      if (smoothValue < 1.0) 
        smoothValue = 1.0;
      sigma = sqrt(ron * ron + smoothValue / gain);
      if (value - smoothValue >= threshold * sigma) cosmic[i + j * xLen] = -1;
    }
  }
  deleteImage(smoothImage);

  /*
   *  Search for groups of contiguous cosmic rays candidates.
   */

  do {
    found = 0;
    for (pos = first; pos < nPix; pos++) {
      if (cosmic[pos] == -1) {
        cosmic[pos] = 2;         /*  Candidate found.  */
        i = pos % xLen;            /*  Its coordinates.  */
        j = pos / xLen;
        first = pos;
        first++;      /* ???  really necessary? */
        found = 1;
        break;
      }
    }

    if (found) {

     /*
      *  Determine new group of contiguous cosmic rays candidates.
      *  Initialize the working box boundaries, iMin, iMax, jMin, jMax, and
      *  the value of the max pixel and its position, fMax, iPosMax, jPosMax.
      */

      iMin = iMax = iPosMax = i;
      jMin = jMax = jPosMax = j;
      fMax = image->data[i + j * xLen];

      do {
        foundContiguousCandidate = 0;
        for (l = 0; l <= 1; l++) {
          for (k = 0; k <= 1; k++) {

           /* 
            *  Looping on 4 pixels to North, East, South and West 
            */

            ii = i + k - l;
            jj = j + k + l - 1;
            if (cosmic[ii + jj * xLen] == -1) {
              foundContiguousCandidate = 1;
              cosmic[ii + jj * xLen] = 2; 
                                        /* Candidate belongs to current group */
              iii = ii;                 /* Keep its position */
              jjj = jj;
              if (ii < iMin) iMin = ii; /* Upgrade search box */
              if (ii > iMax) iMax = ii;
              if (jj < jMin) jMin = jj;
              if (jj > jMax) jMax = jj;
              if (image->data[ii + jj * xLen] > fMax) {
                fMax = image->data[ii + jj * xLen];
                iPosMax = ii;
                jPosMax = jj;
              }
            }
          }
        }

       /*
        *  We are done exploring the "cross". Now mark as "examined"
        *  the current candidate (at the center of the cross):
        */

        cosmic[i + j * xLen] = 3; /* It could probably be set to 1 right away */
        if (foundContiguousCandidate) {

         /* 
          * Pass (arbitrarily) the coordinates of the LAST found candidate 
          */

          i = iii;
          j = jjj;
          continue;       /* Skip the rest, continue loop on new candidate */
        }

       /*
        *  Look for leftovers in the (growing!) search box
        */

        for (l = jMin; l <= jMax; l++) {
          for (k = iMin; k <= iMax; k++) {
            if (cosmic[k + l * xLen] == 2) {
              i = k;
              j = l;
              foundContiguousCandidate = 1;
              break;
            }
          }
          if (foundContiguousCandidate) break;
        }
      } while (foundContiguousCandidate);

     /*
      *  No more contiguous candidates are found. Decide now
      *  whether the current group is a cosmic ray or not.
      */

      sum = 0.;                /* Sum of 8 pixels around max position */
      for (l = -1; l <= 1; l++) {
        for (k = -1; k <= 1; k++) {
          if (l != 0 || k != 0) {
            sum += image->data[iPosMax + k + (jPosMax + l) * xLen];
          }
        }
      }
      sum /= 8.;
      sum -= skyLevel;
      if (fMax - skyLevel > ratio * sum) {
        for (l = jMin - 1; l <= jMax + 1; l++) {
          for (k = iMin - 1; k <= iMax + 1; k++) {
            if (cosmic[k + l * xLen] == 3) {
              cosmic[k + l * xLen] = 1;
              /*** outputMask->data[k + l * xLen] = 1.; ***/
              numCosmic++;
            }
          }
        }
      }
      else {
        for (l = jMin - 1; l <= jMax + 1; l++) {
          for (k = iMin - 1; k <= iMax + 1; k++) {
            if (cosmic[k + l * xLen] != -1) {
              if (cosmic[k + l * xLen] == 1) numCosmic--;
              cosmic[k + l * xLen] = 0;
            }
          }
        }
      }
    }
  } while (found);


  /* 
   *  Prepare dummy CCD table (containing bad pixels - i.e., cosmic rays -
   *  coordinates). If a (real) CCD table was given in input, add its
   *  bad pixel coordinates to the new table. Add then also an extra
   *  flag column, to tell cosmic rays (= 1) from bad pixels (= 2).
   */

  table = newTable();
  table->descs = newStringDescriptor(pilTrnGetKeyword("Table"), "COSMIC", " ");

  table->numColumns = 2;

  if (ccdTable) {
    vimosDscCopy(&(table->descs), ccdTable->descs,
                 pilTrnGetKeyword("InstrumentMode"), NULL); 
    if (ccdTable->numColumns > 0) {
      table->numColumns++;
      if ((port = getPorts(image, &nports))) {
        shiftX = port->shiftX - 1;
        shiftY = port->shiftY - 1;
        deletePortList(port);
      }
     /*
      *  Never mind any failure of getPorts, continue anyway (no shift)
      */
      numBad = ccdTable->cols->len;
      xBad = ccdTable->cols;
      yBad = ccdTable->cols->next;
    }
    else {
      numBad = 0;
    }
  }
  else {
    numBad = 0;
  }

  cpl_msg_debug(modName, "Cosmic found: %d", numCosmic);

  numPix = numCosmic + numBad;

  table->cols = newIntColumn(numPix, "X");
  xColumn = table->cols;
  xColumn->next = newIntColumn(numPix, "Y");
  yColumn = xColumn->next;
  if (cleanJustCosmics) {
   /*
    *  If just cosmics must be cleaned, and a list of bad pixels is
    *  given too, we need an extra column to tell cosmics from bads.
    */
    yColumn->next = newIntColumn(numPix, "level");
    colLevel = yColumn->next;
  }

  for (pos = 0, i = 0; pos < nPix; pos++) {
    if (cosmic[pos] == 1) {
      xColumn->colValue->iArray[i] = (pos % xLen) + 1;
      yColumn->colValue->iArray[i] = (pos / xLen) + 1;
      if (cleanJustCosmics) colLevel->colValue->iArray[i] = cosmicFlag;
      i++;
    }
  }

  for (j = 0, i = numCosmic; j < numBad; j++, i++) {
    xColumn->colValue->iArray[i] = xBad->colValue->iArray[j] - shiftX;
    yColumn->colValue->iArray[i] = yBad->colValue->iArray[j] - shiftY;
    if (cleanJustCosmics) colLevel->colValue->iArray[i] = badPixelFlag;
  }

  if (cleanJustCosmics) exitStatus = cleanBadPixels(image, table, cosmicFlag);
  else                  exitStatus = cleanBadPixels(image, table, 0);

  cpl_free(cosmic);
  deleteTable(table);

  return(exitStatus);
}
/**@}*/
