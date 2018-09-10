/* $Id: vmifuutils.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <math.h>

#include <pilmessages.h>
#include <piltranslator.h>

#include "vmimage.h"
#include "vmtable.h"
#include "vmextractiontable.h"


int
ifuExtrShift(VimosImage **stackImages, VimosExtractionTable *extractionTable)  
{
  /* these define where to look for intermodule gap
  to re-shift correctly the ifu extraction table (flexures)
  should be OK for Q2, to be checked for Q1,3,4 */
/*  #define Q1S1pos 320 */
/*  #define Q1S2pos 320 */
/*  #define Q1S3pos 320 */
/*  #define Q1S4pos 320 */
/*  #define Q2S1pos 320 */
/*  #define Q2S2pos 320 */
/*  #define Q2S3pos 320 */
/*  #define Q2S4pos 320 */
/*  #define Q3S1pos 320 */
/*  #define Q3S2pos 320 */
/*  #define Q3S3pos 320 */
/*  #define Q3S4pos 320 */
/*  #define Q4S1pos 320 */
/*  #define Q4S2pos 320 */
/*  #define Q4S3pos 320 */
/*  #define Q4S4pos 320 */

  const char  modName[] = "ifuExtrShift";
  char comment[80];
  int quad, currentSlit, inter, firstSlit, lastSlit;
  int i,j,x,y;
  int startx,nx,starty,ny;
  int    posMin;
  int m, n, leng;

  int xlength, ylength, fib, l, k, startylong, countdead;

  float  maxVal = -99999., minVal;
  float  sigma, offset;

  VimosIntArray *maxPos;
  VimosFloatArray *firstPosX, *firstPosY, *lastPos, *profile, *profilelong;
  VimosExtractionSlit  *slits;
  VimosImage *subimage, *subimagelong;



  if (readIntDescriptor(stackImages[0]->descs, pilTrnGetKeyword("Quadrant"),
                        &quad, comment) == VM_FALSE) {
    pilMsgError(modName, "Cannot read descriptor %s",
                pilTrnGetKeyword("Quadrant"));
    return EXIT_FAILURE;
  }


  /* start loop on pseudoslits: now extract monodim array */
  for (i=0; i<4; i++) {

    /* loop to find first and last pos for profile extraction */
    slits=extractionTable->slits;
    currentSlit = i+1;
    /*    printf("slit: %d \n", currentSlit);*/

    firstPosX = newFloatArray(4);
    firstPosY = newFloatArray(4);
    lastPos = newFloatArray(4);

    /* initialize positon */
    inter = 80;
    /* start loop on the 5 modules of each slit to find start-end positions */
    n = 0;
    while (slits) {
      firstSlit=inter-1;
      lastSlit=inter +2;
      while ((slits->IFUslitNo == currentSlit) && slits) {
	if(slits->IFUfibNo == firstSlit){
	  firstPosX->data[n] = slits->ccdX->data[0];
	  firstPosY->data[n] = slits->ccdY->data[0];
	}
	if(slits->IFUfibNo == lastSlit){
	  lastPos->data[n] = slits->ccdX->data[4];
	  slits = slits->next;

/*  	  printf("%f %f %f \n",
  		 firstPosX->data[n],firstPosY->data[n], lastPos->data[n]);
 */
	  inter += 80;
	  n++;
	  break;
	}
	slits = slits->next;
	if (!slits) break;
      } /* END OF WHILE ON SLITS->IFUSLITNO */
      if (slits) slits = slits->next;
    }  /* END OF WHILE ON SLITS */

    inter = 80;
    /* restart loop on moules for this slit, to do actual computations */
    for (n=0; n<4; n++) {

/*        printf("*********\n"); */

      startx = (int)(firstPosX->data[n]+0.5);
      nx = (int)(lastPos->data[n]+0.5)-startx+1;
      starty= (int)(firstPosY->data[n]+0.5)-150;
      ny = 200;
      subimage = newImageAndAlloc(nx,ny);
      profile = newFloatArray(nx);
      subimage->data=extractFloatImage(stackImages[0]->data,
				       stackImages[0]->xlen,
				       stackImages[0]->ylen,
				       startx,starty, nx, ny );
      for (x = 0; x < nx; x++) {
	for (y = 0; y < ny; y++) {
	  profile->data[x] += subimage->data[x + nx*y];
	}
      }

      /* and find minimum */
      minVal=profile->data[0];
      posMin=0;
      for (x = 1; x < nx; x++) {
	if (profile->data[x] < minVal ) {
	  minVal = profile->data[x];
	  posMin = x;
	}
      }


      xlength = stackImages[0]->xlen -(startx+posMin);
      ylength = 50;
      startylong = (int)(firstPosY->data[n]+0.5)-(ylength/2);

      subimagelong = newImageAndAlloc(xlength,ylength);
      profilelong = newFloatArray(xlength);
      subimagelong->data=extractFloatImage(stackImages[0]->data,
					   stackImages[0]->xlen,
					   stackImages[0]->ylen,
					   (startx+posMin),startylong, 
					   xlength,ylength );

      leng = subimagelong->xlen;
      for (x = 0; x < xlength; x++) {
	for (y = 0; y < ylength; y++) {
	  profilelong->data[x] += subimagelong->data[x + leng*y];
	}
      }
      
      /* from here towards direction, find maximum 
       direction is right if posMin is towards the left of the array
       and viceversa also check to be at least N sigma above minimum
       to avoid gap fluctuations*/

      sigma = 10*sqrt(minVal);
      
      maxVal = -99999.;
      maxPos = newIntArray(80);

      slits = extractionTable->slits;
      while (slits->IFUslitNo < (i+1))
	{
	  slits = slits->next;
	}
      fib = inter + 1;
      while (slits->IFUfibNo != fib) {
	slits = slits->next;
      }
      
      /* count if and howmany fibers are dead at the beginning of the module */
      countdead = 0;
      while(slits->IFUfibTrans < 0.0)
	{
	  countdead++;
	  slits = slits->next;
	}

/*        printf("%f \n",slits->ccdY->data[2]); */

      /* direction is right: find first good peak */
      if (slits->IFUfibTrans >= 0.0)
	{
	  for (j=0; j<profilelong->len; j++) {
	    if (profilelong->data[j] >= maxVal ) {
	      if (profilelong->data[j] > (sigma+minVal)) {
		maxVal = profilelong->data[j];
		maxPos->data[0+countdead] = j;

/*  		printf("%f  %d \n",maxVal, maxPos->data[0+countdead]); */

	      }
	    } else {
	      break;
	    }
	  }
	}
      else{
	
	maxPos->data[0+countdead] = slits->ccdX->data[2];
      }

      /* calc offset for this fiber */
      offset = (float)((int)(slits->ccdX->data[2]) -
		       (startx+posMin+maxPos->data[0+countdead]));
      
      /* apply offset to all positions of this fiber */
      for (l=0; l<slits->ccdX->len; l++) slits->ccdX->data[l] -=offset;
      
      /* go to next fiber in the module */
      slits = slits->next;
      

      /* chech this!!!!! if first fib not good, this can be very offsetted */
      
      for (l=1+countdead; l<80; l++)
	{
	  k = maxPos->data[l-1] + 3;
	  
	  maxVal = -99999.;
	  
	  if (slits->IFUfibTrans >= 0.0)
	    {
	      for (j=k; j<profilelong->len; j++) {
		if (profilelong->data[j] >= maxVal ) {
		  if (profilelong->data[j] > (sigma+minVal)) {
		    maxVal = profilelong->data[j];
		    maxPos->data[l] = j;
		  }
		} else {
		  break;
		}
	      }
/*  	      printf("%d %d \n",(l+1),(startx+posMin)+maxPos->data[l]); */
	      
	    }
	  else {
	    maxPos->data[l] = maxPos->data[l-1] +5;

/*  	    printf("pppp %d %d \n",(l+1),(startx+posMin)+maxPos->data[l]); */

	  }
	  offset = (float)((int)(slits->ccdX->data[2]) - 
			   (startx+posMin+maxPos->data[l]));
	  for (m=0; m<slits->ccdX->len; m++) slits->ccdX->data[m] -=offset;
	  
	  slits = slits->next;
	  
	}

/*   
     for (j=0; j<80; j++) 
     printf("%d %d \n",(j+1),(startx+posMin)+maxPos->data[j]);
*/

      inter += 80;

      deleteIntArray(maxPos);
      deleteImageAndAlloc(subimagelong);
      deleteFloatArray(profilelong);
      deleteImageAndAlloc(subimage);
      deleteFloatArray(profile);
      
    } /* end loop on modules for this pseudoslit */


    /* SPECIAL LOOP FOR FIRST MODULE IN THE SLIT */

    inter = 80;

/* restart loop on moules for this slit, to do actual computations */

/*      printf("*********\n"); */

    startx = (int)(firstPosX->data[0]+0.5);
    nx = (int)(lastPos->data[0]+0.5)-startx+1;
    starty= (int)(firstPosY->data[0]+0.5)-150;
    ny = 200;
    subimage = newImageAndAlloc(nx,ny);
    profile = newFloatArray(nx);
    subimage->data=extractFloatImage(stackImages[0]->data,
				     stackImages[0]->xlen,
				     stackImages[0]->ylen,
				     startx,starty, nx, ny );
    for (x = 0; x < nx; x++) {
      for (y = 0; y < ny; y++) {
	profile->data[x] += subimage->data[x + nx*y];
      }
    }
    
    /* and find minimum */
    minVal=profile->data[0];
    posMin=0;
    for (x = 1; x < nx; x++) {
      if (profile->data[x] < minVal ) {
	minVal = profile->data[x];
	posMin = x;
      }
    }
    
    
    xlength =startx+posMin;
    ylength = 50;
    startylong = (int)(firstPosY->data[0]+0.5)-(ylength/2);
    
    subimagelong = newImageAndAlloc(xlength,ylength);
    profilelong = newFloatArray(xlength);
    subimagelong->data=extractFloatImage(stackImages[0]->data,
					 stackImages[0]->xlen,
					 stackImages[0]->ylen,
					 0,startylong, 
					 xlength, ylength );
    for (x = 0; x < xlength; x++) {
      for (y = 0; y < ylength; y++) {
	profilelong->data[x] += subimagelong->data[x + xlength*y];
      }
    }
    
    /* from here towards direction, find maximum 
       direction is right if posMin is towards the left of the array
       and viceversa also check to be at least N sigma above minimum
       to avoid gap fluctuations*/
    
    sigma = 10*sqrt(minVal);
    
    maxVal = -99999.;
    maxPos = newIntArray(80);
    
    slits = extractionTable->slits;
    while (slits->IFUslitNo < (i+1))
      {
	slits = slits->next;
      }
    fib = inter;
    while (slits->IFUfibNo != fib) {
      slits = slits->next;
    }
    
    /* count if and howmany fibers are dead at the beginning of the module */
    countdead = 0;
    while(slits->IFUfibTrans < 0.0)
      {
	countdead++;
	slits = slits->prev;
      }
    
/*      printf("%f \n",slits->ccdY->data[2]); */

    /* direction is right: find first good peak */
    if (slits->IFUfibTrans >= 0.0)
      {
	for (j=(profilelong->len-1); j>=0; j--) {
	  if (profilelong->data[j] >= maxVal ) {
	    if (profilelong->data[j] > (sigma+minVal)) {
	      maxVal = profilelong->data[j];
	      maxPos->data[79-countdead] = j;
	      /*printf("%f  %d \n",maxVal, maxPos->data[0+countdead]);*/
	    }
	  } else {
	    break;
	  }
	}
      }
    else{
      
      maxPos->data[79-countdead] = slits->ccdX->data[2];
    }

    /* calc offset for this fiber */
    offset = (float)((int)(slits->ccdX->data[2]) -
		     (maxPos->data[79-countdead]));
    
    /* apply offset to all positions of this fiber */
    for (l=0; l<slits->ccdX->len; l++) slits->ccdX->data[l] -=offset;
    
    /* go to next fiber in the module */
    slits = slits->prev;
    
    
    /* chech this!!!!! if first fib not good, this can be very offsetted */
    
    for (l=(78-countdead); l>=0; l--)
      {
	k = maxPos->data[l+1] - 3;
	
	maxVal = -99999.;
	
	if (slits->IFUfibTrans >= 0.0)
	  {
	    for (j=k; j>=0; j--) {
	      if (profilelong->data[j] >= maxVal ) {
		if (profilelong->data[j] > (sigma+minVal)) {
		  maxVal = profilelong->data[j];
		  maxPos->data[l] = j;
		}
	      } else {
		break;
	      }
	    }
/*  	    printf("%d %d \n",(l+1),maxPos->data[l]); */
	    
	  }
	else {
	  maxPos->data[l] = maxPos->data[l+1] - 5;

/*  	  printf("pppp %d %d \n",(l+1),maxPos->data[l]); */

	}
	offset = (float)((int)(slits->ccdX->data[2]) - 
			 (maxPos->data[l]));
	for (m=0; m<slits->ccdX->len; m++) slits->ccdX->data[m] -=offset;
	
	slits = slits->prev;
	
      }
    
    deleteIntArray(maxPos);
    deleteImageAndAlloc(subimagelong);
    deleteFloatArray(profilelong);
    deleteImageAndAlloc(subimage);
    deleteFloatArray(profile);


    deleteFloatArray(firstPosX);
    deleteFloatArray(firstPosY);
    deleteFloatArray(lastPos);

  } /* end loop on pseudoslits */

  return EXIT_SUCCESS;
}


/* this function deletes the first guess for curvature model coefficients
   in image header, setting coefficients to zero.
   This to allow the construction of a good extraction table for IFU spectra */

VimosBool
ifuDeleteCrvMod(VimosImage* image)
{
  const char modName[] = "ifuDeleteCrvMod";
  char      *descName;
  int        i, j, k;
  int orderPol, orderX, orderY;
  VimosBool  wrOK;
  char comment[80];
  
  pilMsgInfo(modName,"Setting to zero curvature model in image header");

  descName = (char *) pilKeyTranslate("CurvatureOrd");

  if ((wrOK = readIntDescriptor(image->descs, descName, &orderPol, comment)) 
      == VM_TRUE) {
    descName = (char *) pilKeyTranslate("CurvatureOrdX");

    if ((wrOK = readIntDescriptor(image->descs, descName, &orderX, comment))
	== VM_TRUE) {
      descName = (char *) pilKeyTranslate("CurvatureOrdY");

      if ((wrOK = readIntDescriptor(image->descs, descName,&orderY, comment))
	  == VM_TRUE) {

	for (i = 0; i <= orderPol; i++) {
	  for (j = 0; j <= orderX; j++) {
	    for (k = 0; k <= orderY; k++) {
	      descName = (char *) pilKeyTranslate("Curvature", i, j, k);
	      if ((wrOK = writeStringDescriptor(&(image->descs), descName,
						"0.0", comment)) == VM_FALSE) {
		pilMsgError(modName, "Cannot set to zero descriptor %s", 
			    descName);
		return VM_FALSE;
	      }
	    }
	  }
	}
      }
      else {
	pilMsgError(modName,"Cannot read orderY");
	return VM_FALSE;
      }
    } 
    else {
      pilMsgError(modName,"Cannot read orderX");
      return VM_FALSE;
    }
  }
		     
  else {
    pilMsgError(modName,"Cannot read orderPol");
    return VM_FALSE;
  }
  
  if (wrOK == VM_FALSE) 
    pilMsgError(modName, "Cannot set to zero descriptor %s", descName);

  return wrOK;
}


 /*
  * We search for two minima in the profile, one on the left side (low
  * pixel number), and one on the right side (high pixel number)
  */

VimosBool
findIfuBorders(VimosFloatArray *profile, double *upper, double *lower)
{
    int    i;
    int    maxPos;
    float  diff, diff1;
    float  maxVal = -99999.;


    /*  Look for the fiber peak. 
     *  Of course, if  the
     *  accuracy of the first guess on the position were to be worse than 3 
     *  pixels) one could still have 2 real maxima in the profile. Then 
     *  this function would be in trouble !!!
     */
    for (i=0; i<profile->len; i++)
    { 
        if (profile->data[i] > maxVal)
        {
            maxVal = profile->data[i];
            maxPos = i;
        }
    }

    /* if peak is at first or last pixel, we are in trouble:
       our first guess is too wrong. This case is NOT handled */

    if (maxPos == 0 || maxPos == (profile->len -1) ) {
        return (VM_FALSE);
    }

    *lower = (double) (maxPos);
    *upper = (double) (maxPos);

    /* now search for the minimum on the left of the peak */
    diff1=-99;
    for (i=maxPos; i>=0; i--) {
        diff = profile->data[maxPos] -  profile->data[i] ;
        if (diff >diff1) {
            *lower = (double) (i);
            diff1=diff;
        }
    } 
  

    /* now search for the minimum on the rigth of the peak */
    diff1=-99;
    for (i=maxPos; i<=profile->len; i++) {
        diff = profile->data[maxPos] -  profile->data[i] ;
        if (diff >diff1) {
            *upper = (double) (i);
            diff1=diff;
        }
    } 

    /* maybe we should check that distance lower-upper is 5 pixels */

    return (VM_TRUE);

}
