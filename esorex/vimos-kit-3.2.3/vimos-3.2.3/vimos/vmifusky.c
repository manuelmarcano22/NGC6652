/* $Id: vmifusky.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <pilmessages.h>

#include "vmimage.h"
#include "vmmatrix.h"
#include "vmifutable.h"
#include "vmobjecttable.h"
#include "vmmath.h"
#include "vmfit.h"
#include "vmifusky.h"


#define MEDIAN  1
#define AVERAGE 2
#define SUB_STEP 10.

#define PSEUDOSLIT 400


VimosImage *
VmIfuSky(VimosImage *imageData, VimosIfuTable *ifuTable,
         VimosObjectTable *objectTable, char *Combine, int quadNo)
{
  int objNumber, groupFibs;
  int j, k, kk, in, li, io, uu, skyObjs;
  int imageXlen, imageYlen;
  int Comb_Meth = MEDIAN;
  int nGroups, gn, xArrayLen;

  int nDeadFibs, nGoodFibs, never, used, nTotGood, nTotDead, i;

  int numBins, pp, index1, last1, tot;

  float ave, sigma, histoStep, histoXMax, histoYMax;
  float minimum, maximum, threshold, wLenInc;
  char comment[80];
  char modName[] = "VmIfuSky";

  VimosIntArray *sortIndex, *checkFibs;
  VimosIntArray *fiberNum, *fiberNumSort;
  VimosFloatArray *skySpec, *spatSpec, *tmpSpec;
  VimosFloatArray *spectralFlux, *spectralFluxSort;
  VimosFloatArray *histoArray, *histoXArray;


  VimosImage *skyData;

  VimosIfuQuad *ifuQuads;
  VimosIfuSlit  *ifuSlits;
  VimosIfuFiber *ifuFibers;
  VimosObjectObject *objects;

  pilMsgInfo (modName, "Computing Sky Spectra");

  if (!strcmp(Combine,"MEDIAN"))
    {
     Comb_Meth = MEDIAN;
    }
  else if (!strcmp(Combine, "AVERAGE"))
    {
     Comb_Meth = AVERAGE;
    }
  else
    {
     pilMsgError(modName, "Unable to set Combine parameter");
     return NULL;
    }

  /* get sizes: imageXlen is the length of spectra in wavelength */  
  imageXlen = imageData->xlen;
  imageYlen = imageData->ylen;

  /* read wavelength increment in spectra */
  readFloatDescriptor(objectTable->descs,"ESO PRO WLEN INC",&wLenInc,comment);

  /* read number of fiber groups (same for each slit), a different sky */
  /* spectrum is to be computed for each group */
  readIntDescriptor(ifuTable->descs,"ESO PRO SKYGROUP",&nGroups,comment);

  /* create output image for sky spectra: X and Y sizes are same as input */
  skyData = newImageAndAlloc(imageXlen, imageYlen);

  /* initialize skyData to zero everywhere, so that for dead fibers the */
  /* sky spectrum will be zero */
  for (i=0; i<(imageXlen*imageYlen); i++) skyData->data[i] = 0.;

  /* explicit initialization */
  skySpec = NULL;
  histoArray = NULL;
  histoXArray = NULL;
  spectralFlux = NULL;
  spectralFluxSort = NULL;
  fiberNum = NULL;
  fiberNumSort = NULL;
  sortIndex = NULL;
  spatSpec = NULL;

  /* define the "mean" sky spectrum */
  skySpec=newFloatArray(imageXlen);

  /* this array will be used to check that a given fiber is not used in more */
  /* than one group...one never knows...Final values will be: 0 for fibers */
  /* use once; <0 for fibers used more than once; 1 for fibers never used */
  /* (dead ones) */
  checkFibs= newIntArray(PSEUDOSLIT);


  /* buffer to be used to compute integrated flux */
  tmpSpec=newFloatArray(imageXlen);

  objNumber = 0;

  /* counters for total n. of good and dead fibers on this image */
  nTotGood = 0;
  nTotDead = 0;

  /* loop on the IFU Table */
  ifuQuads = ifuTable->quads;
  while (ifuQuads)
    {
     /* take the pseudoslits in the quadrant this image refers to */
     if (ifuQuads->quadNo == quadNo)
       {
        ifuSlits = ifuQuads->ifuSlits;
        while(ifuSlits)
          {
	   /* count total number of dead fibs in this slit, for later check */
	   nDeadFibs = 0;
	   ifuFibers = ifuSlits->fibers;
	   while(ifuFibers)
	     {
	      if (ifuFibers->fiberTrans == -1.) nDeadFibs +=1;
	      ifuFibers = ifuFibers->next;
	     }

	   /* set the values in checkFibs to 1, they will be decremented */
	   /* later */
	   for (uu=0; uu<checkFibs->len; uu++) checkFibs->data[uu] = 1;  

	   /* counter for total number of good fibers in this slit */
	   used = 0;

	   /* for each group in the slit */
	   for (gn=1; gn<=nGroups; gn++)
	      {
	       pilMsgInfo(modName,"SLIT: %3d GROUP: %3d \n",
			  ifuSlits->ifuSlitNo,gn);

	       /* cleanup previous allocations, if they exist */
	       deleteIntArray(fiberNum);
	       deleteIntArray(fiberNumSort);
	       deleteIntArray(sortIndex);
	       deleteFloatArray(spatSpec);
	       deleteFloatArray(spectralFlux);
	       deleteFloatArray(spectralFluxSort);
	       deleteFloatArray(histoXArray);
	       deleteFloatArray(histoArray);

	       histoArray = NULL;
	       histoXArray = NULL;
	       spectralFlux = NULL;
	       spectralFluxSort = NULL;
	       fiberNum = NULL;
	       fiberNumSort = NULL;
	       sortIndex = NULL;
	       spatSpec = NULL;


	       /* initialize sky spectrum for this group */
	       for (k=0; k<imageXlen; k++) skySpec->data[k] = 0.;

	       /* counter for number of objects/fibers in the group */
	       groupFibs = 0;

	       /* take GOOD fibers in the slit and count those belonging to */
	       /* this group: this is needed for setting the size of */
	       /* the fiberNum and spectralFlux arrays (see below) */
	       ifuFibers = ifuSlits->fibers;
	       while(ifuFibers)
		 {
		  if ( (ifuFibers->sigmaYGroup == gn) &&
		       (ifuFibers->fiberTrans != -1.) ) groupFibs++;
	          ifuFibers = ifuFibers->next;
		 }
	       pilMsgInfo(modName,"Slit: %d, n. good fibs: %d",
			  ifuSlits->ifuSlitNo, groupFibs);

	       /* initialize arrays for fiber number in slit and integrated */ 
	       /* flux. Set also arrays for SORTED fiber number and */
	       /* integrated flux */

	       fiberNum=newIntArray(groupFibs);
	       fiberNumSort=newIntArray(groupFibs);
	       sortIndex=newIntArray(groupFibs);
	       spectralFlux=newFloatArray(groupFibs);
	       spectralFluxSort=newFloatArray(groupFibs);

	       
	       for (pp=0; pp<groupFibs; pp++)
		 {
		  fiberNum->data[pp] = 0;
		  fiberNumSort->data[pp] = 0;
		  sortIndex->data[pp] = 0;
		  spectralFlux->data[pp] = 0.0;
		  spectralFluxSort->data[pp] = 0.0;
		 }

	       /* begin selection of fibers and objects to compute sky */
	       /* spectrum for the given group */

	       li = 0;
	       ifuFibers = ifuSlits->fibers;
	       while(ifuFibers)
		 {
		  /* do not take dead fibers */
		  if (ifuFibers->fiberTrans != -1.)
		    {

		     /* take fiber if it is in the group */
		     if (ifuFibers->sigmaYGroup == gn)
		       {
			/* look for the object corresponding to the fiber */
			objects = objectTable->objs;

			while(objects)
			  {
			   /* if object is in the group compute integrated */
			   /* flux and save its fiber number */
			   if ( (objects->IFUslitNo == ifuSlits->ifuSlitNo) &&
				(objects->IFUfibNo == ifuFibers->fibNo) )
			     {
			      /* store the fiber number of this object */
			      fiberNum->data[li] = objects->IFUfibNo;

			      /* initialize */
			      for (in = 0; in < imageXlen; in++)
				tmpSpec->data[in] = 0.0;

			      /* take the spectrum of this object */
			      objNumber = objects->rowNum;
			      for (in = 0; in < imageXlen; in++)
				{
				 tmpSpec->data[in] = 
				   imageData->data[in+objNumber*imageXlen];
				}

			      /* compute integrated flux and store value */
			      spectralFlux->data[li] = 
				integrateSpectrum(tmpSpec,wLenInc);

			      /* decrement checkFibs to mark fibers already */
			      /* used for sky computation */
			      checkFibs->data[(objects->IFUfibNo)-1] -= 1;

			      /* increment counter for object */
			      li++;

			     } /* end if object is in the group */

			   objects = objects->next;
			  } /* end loop on objects */

		       } /* end if the fiber is in the group */


		    } /* end if it is a good fiber (not dead one) */

		  ifuFibers = ifuFibers->next;

		 } /* end loop on fibers */

	       /* pilMsgInfo(modName,"LI: %3d",li); */

	       if (li != groupFibs)
		 {
		   pilMsgError(modName,"ERROR!");
		   return NULL;
		 }

	       /* now for this fiber group we have have two arrays, one for */
	       /* integrated flux and one for fiber number, this last to */
	       /* identify objects later */

	       /* start operations on spectralFlux array to select fibers */
	       /* for sky computation : evaluate r.m.s. of integrated fluxes */
	       sigma = computeRMS(spectralFlux->data, spectralFlux->len);

	       /* set step for the flux histogram. SUB_STEP defined  as 10. */
	       histoStep = sigma / SUB_STEP;

	       /* sort arrays and get minimum and maximum flux */
	       Indexx(spectralFlux->len, spectralFlux->data, sortIndex->data);

	       /* sorting */
	       for (io=0; io<spectralFlux->len; io++)
		 {
		  index1 = sortIndex->data[io];
	          spectralFluxSort->data[io] = 
		    spectralFlux->data[index1];
		  fiberNumSort->data[io]= fiberNum->data[index1];
		 }

	       last1 = spectralFlux->len - 1;
	       minimum = spectralFluxSort->data[0];
	       maximum = spectralFluxSort->data[last1];

	       /* compute flux histogram */
	       numBins = (int)((maximum - minimum) / histoStep + 1.);
	       histoArray= newFloatArray(numBins);
	       histoXArray = newFloatArray(numBins);

	       for (io=0; io<numBins; io++)
		 {
		  histoXArray->data[io] = minimum + io*histoStep;
		 }
	       xArrayLen = histoXArray->len;


	       /*
		 if (ifuSlits->ifuSlitNo == 1)
		 {
		   for (pp=300; pp<400; pp++)
		     pilMsgInfo(modName,"%15.5f\n",spectralFlux->data[pp]);
		   pilMsgInfo(modName,"m: %10.5f M: %10.5f S: %10.5f\n",
		   minimum,maximum,histoStep);
		 }
	       */


	       computeHistogram(spectralFlux,xArrayLen,
				histoArray,minimum,maximum,
				histoStep);

	       /*
	       if (ifuSlits->ifuSlitNo == 1)
		 {
		   for (pp=0; pp<histoArray->len; pp++)
		     pilMsgInfo(modName,"%15.5f\n",histoArray->data[pp]);
		 }
	       */

	       tot = 0.;

	       for (pp=0; pp<xArrayLen; pp++)
		 {
		  /* pilMsgInfo(modName,"I: %3d X: %15.5f  Y: %15.5f\n",pp,
		     histoXArray->data[pp],histoArray->data[pp]);
		  */
		  tot += histoArray->data[pp];
		 }
	       /* pilMsgInfo(modName,"TOT: %4d\n", tot); */
	       if (tot != groupFibs)
		 {
		  pilMsgError(modName,"Wrong no. of histogram elements");
		  return NULL;
		 }


	       /* look for maximum in the histogram */
	       histoYMax = histoArray->data[0];
	       histoXMax = histoXArray->data[0];
	       for (io = 1; io < histoArray->len ; io++)
		  {
		   if (histoArray->data[io] > histoYMax)
		     {
		      histoYMax = histoArray->data[io];
		      /* save max flux */
		      histoXMax = histoXArray->data[io];
		     }
		  }

	       /* define threshold as: */
	       /* ("midpoint flux" + (10% of "midpoint flux")) */
	       threshold=(histoXMax+(histoStep/2.)) + 
		           0.1*(histoXMax+(histoStep/2.));

	       /* look into the sorted flux array and find where flux starts */
	       /* to be >= than threshold, to count number of "good" fibers */
	       /* for sky determination */
	       skyObjs = 0;
	       for (io=0; io<spectralFluxSort->len ; io++)
		 {
		  if (spectralFluxSort->data[io] < threshold) skyObjs++;
		 }


	       /* start the actual computation of the sky spectrum for this */
	       /* group of fiber */

	       /* first create the "spatial" sky: each element is the flux */
	       /* at a given pixel from each one of the selected fibers */

	       /* work array for the spatial sky spectrum */
	       spatSpec = newFloatArray(skyObjs);

	       /* for each pixel along wavelength */
	       for (k=0; k<imageXlen; k++)
		  {
		   /* initialize spatial sky spectrum */
		   for (kk=0; kk<skyObjs; kk++) spatSpec->data[kk] = 0.;

		   /* restart loop on objects */
		   objects = objectTable->objs;

		   /* index for spatSpec array, -->MUST<-- go from 0 to */
		   /* (skyObjs-1) */
		   j = 0;
		   while(objects)
		     {
		      /* if the object is in the right slit */
		      if (objects->IFUslitNo == ifuSlits->ifuSlitNo)
			{
			 /* check if it is in the "sky fibers" list */
			 /* as we use sorted arrays we know we must stop */
			 /* at the "skyObjs" fiber in the array */
			 for (io=0; io<skyObjs; io++)
			    {
			     if (objects->IFUfibNo == fiberNumSort->data[io])
			       {
				objNumber = objects->rowNum;
			        spatSpec->data[j] = 
			          imageData->data[k+objNumber*imageXlen];
			        /* increment counter */
			        j++;

			       } /* end if this object is ok for sky */
			    } /* end loop on fibers to be used for sky comp. */
			} /* end if object is in the slit */

		      objects = objects->next;
		    } /* end loop on objects */

		   /* start combining spatSpec to get the sky value at each */
		   /* Y pixel */
		   ave = 0.;
		   switch (Comb_Meth)
		     {
		     case AVERAGE:
		       {
			for (kk=0; kk<skyObjs; kk++) ave += spatSpec->data[kk];
			ave /= ((float)skyObjs);
			break;
		       }
		    default:
		    case MEDIAN:
		      {
		       ave = medianWirth(spatSpec->data, skyObjs);
		       break;
		      }
		     }

		   /* write pixel value in the sky spectrum for this group */
		   skySpec->data[k] = ave;

		  } /* end loop on lambda: sky spectrum completed */

	       pilMsgInfo(modName,"group: %d, n. used for sky: %d",
			  gn,skyObjs);

	       /* restart to write the computed sky spectrum in skyData */
	       /* image at the rows corresponding to each object belonging */
	       /* to this group */
	       objects = objectTable->objs;

	       while(objects)
		 {
		  if (objects->IFUslitNo == ifuSlits->ifuSlitNo)
		    {
		     /* if the object is in the group: this time use the */
		     /* whole range, groupFibs instead of skyObjs, to write */
		     /* this sky spectrum for all the fibers in the group */
		     for (io=0; io<groupFibs; io++)
		        {
			 if (objects->IFUfibNo == fiberNumSort->data[io])
			   {
		            /* write this sky spectrum in skyData at the */
		            /* same row as object */
			    objNumber = objects->rowNum;
			    for (kk=0; kk<imageXlen; kk++)
			        skyData->data[kk + objNumber*imageXlen] = 
				                         skySpec->data[kk];
			   }

			} /* end loop on groupFibs */

		    } /* end if object is in the right slit */

		  objects = objects->next;
		 } /* end loop on objects */
		  

	       for (pp=0; pp<groupFibs; pp++)
		 {
		  spectralFluxSort->data[pp] = 0.0;
		 }

	       /* increment total number of good fibers for this slit */
	       used += groupFibs;

	      } /* close loop on groups for this slit */

	   /* check that every fiber is not used more than once for sky */
	   /* computation in this slit, and count how many fibers are never */
	   /* used (dead ones) */
	   never = 0;
	   nGoodFibs = 0;
	   for (uu=0; uu<checkFibs->len; uu++)
	      {
	       if (checkFibs->data[uu] != 0)
		 {
		  if (checkFibs->data[uu] < 0) 
		    {
		     pilMsgInfo(modName,"Slit %d: Fib %d used more than once",
				ifuSlits->ifuSlitNo, uu+1);
		    }
		  else if (checkFibs->data[uu] > 0)
		    {
		     never += 1;
		     pilMsgInfo(modName,"Slit %3: Fiber no. %d never used",
				ifuSlits->ifuSlitNo, uu+1);
		    }
		 }
	       else if (checkFibs->data[uu] == 0) nGoodFibs += 1;
	      }
	   /* some checks */
	   if (never != nDeadFibs)
	     {
	      pilMsgError(modName,"Wrong n. of dead fibers in pseudoslit %d",
			  ifuSlits->ifuSlitNo);
	      return NULL;
	     }
	   if (used != nGoodFibs)
	     {
	      pilMsgError(modName,"Wrong n. of good fibers in pseudoslit %d",
			  ifuSlits->ifuSlitNo);
	      return NULL;
	     }
	   if ((never+used) != (PSEUDOSLIT))
	     {
	      pilMsgError(modName,"ERROR! n. of good fibers + n. of dead "
                          "fibers in pseudoslit %d != 400",
			  ifuSlits->ifuSlitNo);
	      return NULL;
	     }

	   pilMsgInfo(modName,"slit: %d, Dead: %d, Good: %d",
			  ifuSlits->ifuSlitNo, never, used);
	   pilMsgInfo(modName,"\n");

	   nTotGood += used;
	   nTotDead += never;

	   ifuSlits = ifuSlits->next;
	  }  /* close loop on slits */

       } /* end if is the right quadrant */

     ifuQuads = ifuQuads->next;
    } /* end loop on ifu quadrants */

  /* some more checks */
  if ( (nTotGood+nTotDead) != (PSEUDOSLIT*4) )
    {
     pilMsgError(modName,"Wrong total n. of fibers used in this image");
     return NULL;
    }

  copyAllDescriptors(imageData->descs, &(skyData->descs)); 

  deleteFloatArray(tmpSpec);
  deleteIfuQuad(ifuQuads);
  deleteIfuSlit(ifuSlits);
  deleteIfuFiber(ifuFibers);
  deleteObjectObject(objects);
  deleteFloatArray(skySpec);
  deleteIntArray(checkFibs);

  return skyData;

}
