/* $Id: vmmosfringes.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <string.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>

#include "vmimage.h"   
#include "vmwindowtable.h"
#include "cpl.h"


#define MAX_COMMENT_LENGTH (80)


int
VmSpFringCorr(VimosImage **mosImages, VimosWindowTable **winTables,
              int mosCount, int pixels, int interpolate)
{
  char        task[]               = "VmSpFringCorr";

  int                   xStart           = 0;
  int                   xEnd             = 0;
  int                   j, i, imaX, imaY;
  int                   left,right;
  int                   badObjects;
  int                   startFound;

  float                 valEnd,valStart,m,value;



  VimosImage           *meanSkyImage     = NULL;
  VimosImage          **subMosImages     = NULL;

  VimosWindowObject    *object;

  VimosWindowSlit      *wSlit            = NULL;

/*    char        comment[MAX_COMMENT_LENGTH]; */
/*    int                   startPixelSub    = 0; */
/*    float                 startWLenSub=7700.;  */
/*    float                 wLenStart,wLenInc; */
/*    double                dValue; */




  subMosImages = (VimosImage **)cpl_calloc(mosCount, sizeof(VimosImage *));
  if (!subMosImages) { 
    cpl_msg_error(task, "Failure creating list of 2D sky images");
    return EXIT_FAILURE;
  }
    
  for (i = 0; i < mosCount; i++) {
    subMosImages[i] = duplicateImage(mosImages[i]);
    wSlit = winTables[i]->slits;
    
    while(wSlit) {
      object = wSlit->objs;
      while(object) {
	for(imaX=0; imaX<subMosImages[i]->xlen; imaX++){
	  left=right=1;
	  
	  /* starting end ending pixel */
	  xStart = wSlit->specStart + object->objStart - pixels;
	  xEnd = wSlit->specStart + object->objEnd + pixels;
	  
	  /* check pixels */
	  if(xStart-2<0) {
	    left=0;
	    xStart=0;
	  }
	  if(xEnd+2>wSlit->specEnd) {
	    right=0;
	    xEnd=wSlit->specEnd;
	  }
	  
	  if(interpolate) {
	    valStart=valEnd=0.;
	    /* simple linear interpolation */
	    if (left) {
	      for(j=xStart-2; j<=xStart; j++) 
		valStart += subMosImages[i]->data[imaX + j*subMosImages[i]->xlen];
	      valStart/=3;
	      value = valStart;
	      m = 0.;
	    }
	    if (right) {
	      for(j=xEnd+2; j>=xEnd; j--) 
		valEnd += subMosImages[i]->data[imaX + j*subMosImages[i]->xlen];
	      valEnd/=3;
	      value = valEnd;
	      m = 0.;
	    }

	    if(left && right) {
	      m = (valEnd-valStart)/((float)(xEnd)-(float)(xStart));
	      value = valStart;
	    }

	    for(imaY=xStart; imaY<=xEnd; imaY++)
	      subMosImages[i]->data[imaX + imaY*subMosImages[i]->xlen] = m*(imaY-xStart) + value;
	  }
	  else 
	    for(imaY=xStart; imaY<=xEnd; imaY++)
	      subMosImages[i]->data[imaX + imaY*subMosImages[i]->xlen] = -32000.;
	}
	object=object->next;
      }
      if (wSlit->next) wSlit = wSlit->next;
      else {
	while(wSlit->prev) wSlit = wSlit->prev; /* rewind wSlit */
	break;
      }
    }
  }
    
    
  /* compute mean sky image */
    
  meanSkyImage = frCombMedian(subMosImages,mosCount,1);
    

  /* Masking of blue wavelenght. FIXME - should be applicated? */
  
  /*
    readDoubleDescriptor(winTables[0]->descs, pilTrnGetKeyword("WlenStart"),&dValue, comment);
    wLenStart = (float) dValue;
    readDoubleDescriptor(winTables[0]->descs, pilTrnGetKeyword("WlenInc"),&dValue, comment);
    wLenInc = (float) dValue;
    startPixelSub = (startWLenSub - wLenStart) / wLenInc;
    
    for (i=0; i< meanSkyImage->xlen; i++)
    for (j=0; j< meanSkyImage->ylen; j++)
    if(i<startPixelSub-1) meanSkyImage->data[i+j*meanSkyImage->xlen] = 0.0;
  */
    


  if(!interpolate) {
    
    /* eventually interpolate in regions covered by the object in all 
       exposures that contains -32000 values even after the combination 
       FIXME - to be verified */

    badObjects=0;
    while(wSlit) {
      object = wSlit->objs;
      while(object) {
	for(imaX=0; imaX<meanSkyImage->xlen; imaX++){
	  startFound = 0;
	  if(!imaX) {
	    for(j=wSlit->specStart; j<=wSlit->specEnd; j++){
	      /* check if a -32000 region exist */
	      if(meanSkyImage->data[imaX + j*meanSkyImage->xlen] == -32000 && !startFound) {
		xStart=j-1;
		startFound=1;
	      }
	      if(meanSkyImage->data[imaX + j*meanSkyImage->xlen] != -32000 && startFound) { 
		xEnd=j;
		break;
	      }
	    }
	    if(startFound) badObjects++;
	  }  
	  
	  if(startFound) {
	    valStart = meanSkyImage->data[imaX + xStart*meanSkyImage->xlen];
	    valEnd = meanSkyImage->data[imaX + xEnd*meanSkyImage->xlen];
	    
	    /* again simple linear interpolation */
	    m = (valEnd-valStart)/((float)(xEnd)-(float)(xStart)); 
	    for(imaY=xStart+1; imaY<xEnd; imaY++)
	      meanSkyImage->data[imaX + imaY*meanSkyImage->xlen] = m*(imaY-xStart) + valStart;
	  }
	}
	object=object->next;  
      }
      wSlit=wSlit->next;
    }
    if(badObjects) cpl_msg_warning(task,"%d objects has been interpolated in central regions", badObjects);    
  }      
  
  /* finally subtract residual image from input images */
  for (i = 0; i < mosCount; i++)        
    imageArithLocal(mosImages[i], meanSkyImage, VM_OPER_SUB);
  
  /* cleaning */  
  for (j = 0; j < mosCount; j++)
    deleteImage(subMosImages[j]);
  cpl_free(subMosImages);
  
  return EXIT_SUCCESS;
}



