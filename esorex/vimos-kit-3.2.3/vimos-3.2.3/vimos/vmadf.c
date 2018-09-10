/* $Id: vmadf.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <fitsio.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>
#include <pilerrno.h>

#include "vmtypes.h"
#include "vmtable.h"
#include "vmdistmodels.h"
#include "vmadf.h"
#include "cpl.h"


/* 
 * Allocates a new VimosBezierCurve 
 */

VimosBezierCurve *newBezierCurve()
{
  const char modName[] = "newBezierCurve";
  VimosBezierCurve *tmpCurve;
  
  /* allocate space */
  tmpCurve = (VimosBezierCurve *) cpl_malloc(sizeof(VimosBezierCurve));
  
  /* check if space was allocated */
  if (tmpCurve == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  /* initialize values */
  tmpCurve->x0 = 0.0;
  tmpCurve->a  = 0.0;
  tmpCurve->b  = 0.0;
  tmpCurve->c  = 0.0;

  /* return pointer */
  return(tmpCurve);
  
}

   

/*
  Deletes a VimosBezierCurve
*/
void deleteBezierCurve(VimosBezierCurve *curve)
{
  /* only have to free the space */
  if (curve != NULL) {
    cpl_free(curve);
  }
  
}


/*
  Evaluate a bezier curve
*/
float computeBezierCurve(VimosBezierCurve *curve, float t)
{
  const char modName[] = "computeBezierCurve";
  float result;
  
  pilErrno = 0;

  /* validate input */
  if (curve == NULL) {
    cpl_msg_error(modName, "NULL input pointer");
    pilErrno = 1;
    return(-1.0);
  }
  
  /* bring parameter within limits */
  if (t < 0.0) {
    t = 0.0;
  }
  if (t > 1.0) {
    t = 1.0;
  }
  
  /*compute value */
  result = curve->x0 + curve->c*t + curve->b*t*t + curve->a*t*t*t;
  
  return(result);
    
}

/*
   Allocates a new VimosAdfRectSlit structure
*/
VimosAdfRectSlit *newAdfRectSlit()
{
  const char modName[] = "newAdfRectSlit";
  VimosAdfRectSlit *tmpSlit;
  
  /* allocate space */
  tmpSlit = (VimosAdfRectSlit *) cpl_malloc(sizeof(VimosAdfRectSlit));
  
  /* check if space was allocated */
  if (tmpSlit == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  /* initialize values */
  tmpSlit->slitType = VM_ADF_RECT_SLIT;
  tmpSlit->slitNo   = 0;
  tmpSlit->x        = 0.0;
  tmpSlit->y        = 0.0;
  tmpSlit->dimX     = 0.0;
  tmpSlit->dimY     = 0.0;

  /* return pointer */
  return(tmpSlit);
  
}


/*
  Deletes a VimosAdfRectSlit
*/
void deleteAdfRectSlit(VimosAdfRectSlit *adfSlit)
{
  /* only have to free the space */
  if (adfSlit != NULL) {
    cpl_free(adfSlit);
  }
  
}



/*
   Allocates a new VimosAdfCurvSlit structure
*/
VimosAdfCurvSlit *newAdfCurvSlit()
{
  const char        modName[] = "newAdfCurvSlit";
  VimosAdfCurvSlit *tmpSlit;
  
  /* allocate space */
  tmpSlit = (VimosAdfCurvSlit *) cpl_malloc(sizeof(VimosAdfCurvSlit));
  
  /* check if space was allocated */
  if (tmpSlit == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  /* initialize values */
  tmpSlit->slitType = VM_ADF_CURV_SLIT;
  tmpSlit->slitNo   = 0;
  tmpSlit->deltay = 0.;

  /* Allocate curves  and check allocations */
  /* middle X curve */
  tmpSlit->xMiddle = newBezierCurve();
  if (tmpSlit->xMiddle == NULL) {
    cpl_free(tmpSlit);
    cpl_msg_error(modName, "The function newBezierCurve has returned NULL (x)");
    return(NULL);
  }
  
  /* middle Y curve */
  tmpSlit->yMiddle = newBezierCurve();
  if (tmpSlit->yMiddle == NULL) {
    cpl_free(tmpSlit);
    cpl_msg_error(modName, "The function newBezierCurve has returned NULL (y)");
    return(NULL);
  }

  /* return pointer */
  return(tmpSlit);
  
}


/*
  Deletes a VimosAdfCurvSlit
*/
void deleteAdfCurvSlit(VimosAdfCurvSlit *slit)
{
  if (slit != NULL) {
    /* delete the Bezier Curves */
    deleteBezierCurve(slit->xMiddle);
    deleteBezierCurve(slit->yMiddle);
    
    /* free the object */
    cpl_free(slit);
  }
  
}

/*
   Allocates a new VimosAdfCircSlit structure
*/
VimosAdfCircSlit *newAdfCircSlit()
{
  const char modName[] = "newAdfCircSlit";
  VimosAdfCircSlit *tmpSlit;
  
  /* allocate space */
  tmpSlit = (VimosAdfCircSlit *) cpl_malloc(sizeof(VimosAdfCircSlit));
  
  /* check if space was allocated */
  if (tmpSlit == NULL) {
    cpl_msg_error(modName, "Allocation error");
    return(NULL);
  }

  /* initialize values */
  tmpSlit->slitType = VM_ADF_CIRC_SLIT;
  tmpSlit->slitNo   = 0;
  tmpSlit->x        = 0.0;
  tmpSlit->y        = 0.0;
  tmpSlit->radius   = 0.0;
  tmpSlit->IFUslitNo = 0;
  tmpSlit->IFUfibNo = 0;
  tmpSlit->IFUfibTrans = 0.0;

  /* return pointer */
  return(tmpSlit);
  
}


/*
  Deletes a VimosAdfCircSlit
*/
void deleteAdfCircSlit(VimosAdfCircSlit *adfSlit)
{
  /* only have to free the space */
  if (adfSlit != NULL) {
    cpl_free(adfSlit);
  }
  
}



/*
   Allocates a new VimosAdfRefrSlit structure
*/
VimosAdfRefrSlit *newAdfRefrSlit()
{
  const char modName[] = "newAdfRefrSlit";
  VimosAdfRefrSlit *tmpSlit;
  
  /* allocate space */
  tmpSlit = (VimosAdfRefrSlit *) cpl_malloc(sizeof(VimosAdfRefrSlit));
  
  /* check if space was allocated */
  if (tmpSlit == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  /* initialize values */
  tmpSlit->slitType = VM_ADF_REFR_SLIT;
  tmpSlit->slitNo   = 0;
  tmpSlit->x        = 0.0;
  tmpSlit->y        = 0.0;
  tmpSlit->size     = 0.0;

  /* return pointer */
  return(tmpSlit);
  
}


/*
  Deletes a VimosAdfRefrSlit
*/
void deleteAdfRefrSlit(VimosAdfRefrSlit *adfSlit)
{

  /* only have to free the space */
  if (adfSlit != NULL) {
    cpl_free(adfSlit);
  }
  
}


VimosAdfSlitHolder *newAdfSlitHolder()
{
  const char modName[] = "newAdfSlitHolder";
  VimosAdfSlitHolder *slitHolder;
  
  slitHolder = (VimosAdfSlitHolder *) cpl_malloc(sizeof(VimosAdfSlitHolder));
  
  if (slitHolder == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  slitHolder->slitType = VM_ADF_UNDF_SLIT;
  
  slitHolder->slit = NULL;
  slitHolder->prev = NULL;
  slitHolder->next = NULL;
  
  return(slitHolder);
    
}

void deleteAdfSlitHolder(VimosAdfSlitHolder *slitHolder) 
{

  VimosAdfSlitHolder *tSlit;
  VimosAdfSlitHolder *nextSlit;

  tSlit = slitHolder;
    
  while(tSlit) {
    switch (tSlit->slitType) {
    case VM_ADF_RECT_SLIT : 
      {
        deleteAdfRectSlit((VimosAdfRectSlit *)tSlit->slit);
        break;
      }
    case VM_ADF_CURV_SLIT : 
      {
        deleteAdfCurvSlit((VimosAdfCurvSlit *)tSlit->slit);
        break;
      }
    case VM_ADF_CIRC_SLIT : 
      {
        deleteAdfCircSlit((VimosAdfCircSlit *)tSlit->slit);
        break;
      }
    case VM_ADF_REFR_SLIT : 
      {
        deleteAdfRefrSlit((VimosAdfRefrSlit *)tSlit->slit);
        break;
      }
    default :
      {
        break;
        
      }
    }
    nextSlit = tSlit->next;
    cpl_free(tSlit);
    tSlit = nextSlit;
  }

}





/*
   Allocate a new ADF. This is a generic VimosTable with no name.
   No descriptors or columns are allocated (but this may change in a later
   version).
   This is a new version (MS)
*/
VimosTable *newADF()
{
  const char modName[] = "newADF";
  VimosTable *newTab;
  
  /* allocate new VimosTable */
  newTab = newTable();

  /* if error, return NULL */
  if (newTab == NULL) {
    cpl_msg_error(modName, "The function newTable has returned NULL");
    return(NULL);
  }

  return(newTab);
}

/* 
 *  The following function is doing exactly what readFitsADF() is doing,
 *  but avoiding the dependency from the FITS file. This was in fact
 *  unnecessary, as the descriptors are already contained in the input
 *  VimosImage, and do not need to be read again from disk. This
 *  routine, as well as readFitsADF, need to be improved in the sense
 *  that the ADF table should consist of just a subset of all the
 *  descriptors found in header, that is just the descriptor related
 *  to ADF.
 */

VimosBool readADF(VimosTable *adf, VimosImage *adfImage)
{
  char descVal[81];

  /* validate input */
  if (adf == NULL) {
    return(VM_FALSE);
  }

  copyAllDescriptors(adfImage->descs, &(adf->descs));
  readStringDescriptor(adf->descs, "ESO INS ADF TYPE", descVal, NULL);

  /* make sure ADF types are ok */
  if ( !strncmp(descVal, "MOS", 3) ) {
    strcpy(adf->name, VM_ADF_MOS);
  }
  if ( !strncmp(descVal, "IFU", 3) ) {
    strcpy(adf->name, VM_ADF_IFU);
  }
  if ( !strncmp(descVal, "IMA", 3) ) {
    strcpy(adf->name, VM_ADF_IMA);
  }

  return(VM_TRUE);
}

/*
  Delete an ADF. This is just an esthetic wrapper for deleteTable(adf)
*/
void deleteADF(VimosTable *adf)
{
  if (adf != NULL) {
    /* just a wrapper for deleteTabel() */
    deleteTable(adf);
  }
  
}


/* 
   Extract the slits apertures from the ADF. The function allocates an
   array of *voids of the size required and returns the pointer to this array
   through slitAr. It also allocates an array of VimosAdfSlitType of which the
   pointer is returned through typeAr
   For each slit defined in the ADF, the function creates a VimosAdfRectSlit,
   a VimosAdfCurvSlit or a VimosAdfCircSlit and stores the pointer to this 
   slit in the array **slitAr. It also writes the type of this slit in *typeAr.

One day I should write this more compact...

*/
VimosAdfSlitHolder *extractSlitsFromADF(VimosTable *adf)
{
  const char           modName[] = "extractSlitsFromADF";
  int                  numSlits;
  int                  i;
  char                 strValue[80];
  double               dValue;
  VimosBool            rdOK;
  VimosAdfType         adfType;
  VimosAdfSlitType     slitType=VM_ADF_UNDF_SLIT;
  VimosAdfSlitHolder   *slitHolder;
  VimosAdfSlitHolder   *linkedSlits=NULL;
  VimosAdfSlitHolder   *tSlitHolder;
  VimosAdfRectSlit     *tRectSlit;
  VimosAdfCurvSlit     *tCurvSlit;
  VimosAdfCircSlit     *tCircSlit;
  void                 *tSlit;
  

  /* validate input */
  if (adf == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return(NULL);
  }

  /* initialize to undefined ADF type*/
  adfType = VM_ADF_TYPE_UDF;
  /* only work on ADFs */
  if ( !strcmp(adf->name, VM_ADF_MOS) ) {
    adfType = VM_ADF_TYPE_MOS;
  } 
  if ( !strcmp(adf->name, VM_ADF_IFU) ) {
    adfType = VM_ADF_TYPE_IFU;
  } 
  if ( !strcmp(adf->name, VM_ADF_IMA) ) {
    adfType = VM_ADF_TYPE_IMA;
  } 

  /* read how many slits are defined in the ADF. Every type of ADF should have 
   this descriptor */
  rdOK = readIntDescFromTable(adf, pilTrnGetKeyword("NoSlit"), &numSlits,
                              NULL);
  if (!rdOK) {
    cpl_msg_error(modName, "Cannot read descriptor %s",
                pilTrnGetKeyword("NoSlit"));
    return(NULL);
  }
            
  slitHolder = NULL;
  
  switch (adfType) {
    /* No recognized ADF */
  default : 
    {
      cpl_msg_error(modName, "Unrecognized ADF type");
      return(NULL);
    }
    /* MOS ADF*/
  case VM_ADF_TYPE_MOS : 
    {
      for (i = 0; i < numSlits; i++) {
        rdOK =  readStringDescFromTable(adf, pilTrnGetKeyword("SlitType", i+1),
                                        strValue, NULL);
        if (!rdOK) {
          deleteAdfSlitHolder(slitHolder);
          cpl_msg_error(modName, "Cannot read descriptor %s from table %s", 
                      pilTrnGetKeyword("SlitType", i+1), adf->name);
          return(NULL);
        }
        slitType = VM_ADF_UNDF_SLIT;
        if ( !strcmp(strValue, VM_ADF_RECT_SLIT_NAME) ) {
          slitType = VM_ADF_RECT_SLIT;
        }
        if ( !strcmp(strValue, VM_ADF_CURV_SLIT_NAME) ) {
          slitType = VM_ADF_CURV_SLIT;
        }
        
        switch (slitType) {
        case VM_ADF_RECT_SLIT : 
          {
            /* is a rectangular slit */
            tRectSlit =  newAdfRectSlit();
            if (tRectSlit == NULL) {
              deleteAdfSlitHolder(slitHolder);
              cpl_msg_error(modName, "Function newAdfRectSlit returned NULL");
              return(NULL);
            }
            tRectSlit->slitNo = i+1;
            rdOK = readDoubleDescFromTable(adf,
                                           pilTrnGetKeyword("SlitX", i+1),
                                           &dValue, NULL);
            tRectSlit->x = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitY", i+1),
                                           &dValue, NULL);
            tRectSlit->y = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitDimX", i+1),
                                           &dValue, NULL);
            tRectSlit->dimX = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitDimY", i+1),
                                           &dValue, NULL);
            tRectSlit->dimY = (float) dValue;

            tSlit = (void *) tRectSlit;
            break; 
          } 
        case VM_ADF_CURV_SLIT : 
          {
            /* is a curved slit */
            tCurvSlit =  newAdfCurvSlit();
            if (tCurvSlit == NULL) {
              deleteAdfSlitHolder(slitHolder);
              cpl_msg_error(modName, "Function newAdfCurvSlit returned a NULL");
              return(NULL);
            }        
            tCurvSlit->slitNo = i+1;

            rdOK = readDoubleDescFromTable(adf,
                                           pilTrnGetKeyword("SlitBezierDY",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->deltay = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitBezierXX",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->xMiddle->x0 = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitBezierAX",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->xMiddle->a = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitBezierBX",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->xMiddle->b = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitBezierCX",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->xMiddle->c = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitBezierYY",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->yMiddle->x0 = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitBezierAY",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->yMiddle->a = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitBezierBY",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->yMiddle->b = (float) dValue;

            rdOK = rdOK
                && readDoubleDescFromTable(adf, 
                                           pilTrnGetKeyword("SlitBezierCY",
                                                            tCurvSlit->slitNo),
                                           &dValue, NULL);
            tCurvSlit->yMiddle->c = (float) dValue;

            tSlit = (void *) tCurvSlit;
            break;
          }
        default : 
          {
            cpl_msg_error(modName, "Unrecognized type of slit");
            deleteAdfSlitHolder(slitHolder);
            return(NULL);
          }
        } /* end switch slitType */

        if (i == 0) {
          slitHolder = newAdfSlitHolder();
          if (slitHolder == NULL) {
            cpl_msg_error(modName, "Function newAdfSlitHolder returned a NULL");
            return(NULL);
          }
          linkedSlits = slitHolder;
          linkedSlits->prev = NULL;
          linkedSlits->next = NULL;
        }
        else {
          tSlitHolder = newAdfSlitHolder();
          if (tSlitHolder == NULL) {
            cpl_msg_error(modName, "Function newAdfSlitHolder returned a NULL");
            return(NULL);
          }
          tSlitHolder->prev = linkedSlits;
          linkedSlits->next = tSlitHolder;
          linkedSlits = tSlitHolder;
        }
        linkedSlits->slit = tSlit;
        linkedSlits->slitNo = i+1;
        linkedSlits->slitType = slitType;
        if (!rdOK) {
          deleteAdfSlitHolder(slitHolder);
          cpl_msg_error(modName, "In function readDoubleDescFromTable");
          return(NULL);
        }     

      }
      break;
    }
  case VM_ADF_TYPE_IMA :
  case VM_ADF_TYPE_IFU : 
    {
      for (i = 0; i < numSlits; i++) {
        /* is a circular slit */
        tCircSlit =  newAdfCircSlit();
        if (tCircSlit == NULL) {
          cpl_msg_error(modName, "Function newAdfCircSlit returned a NULL");
          return(NULL);
        }
        tCircSlit->slitNo = i+1;

        rdOK = readDoubleDescFromTable(adf,
                                       pilTrnGetKeyword("SlitX", i+1),
                                       &dValue, NULL);
        tCircSlit->x = (float) dValue;

        rdOK = rdOK
            && readDoubleDescFromTable(adf,
                                       pilTrnGetKeyword("SlitY", i+1),
                                       &dValue, NULL);
        tCircSlit->y = (float) dValue;

        rdOK = rdOK
            && readDoubleDescFromTable(adf, 
                                       pilTrnGetKeyword("SlitRadius", i+1),
                                       &dValue, NULL);
        tCircSlit->radius = (float) dValue;


        if (i == 0) {
          slitHolder = newAdfSlitHolder();
          if (slitHolder == NULL) {
            cpl_msg_error(modName, "Function newAdfSlitHolder returned a NULL");
            return(NULL);
          }
          linkedSlits = slitHolder;
          linkedSlits->prev = NULL;
          linkedSlits->next = NULL;
        }
        else {
          tSlitHolder = newAdfSlitHolder();
          if (tSlitHolder == NULL) {
            cpl_msg_error(modName, "Function newAdfSlitHolder returned a NULL");
            return(NULL);
          }
          tSlitHolder->prev = linkedSlits;
          linkedSlits->next = tSlitHolder;
          linkedSlits = tSlitHolder;
        }
        linkedSlits->slit = tCircSlit;
        linkedSlits->slitNo = i+1;
        linkedSlits->slitType = slitType;
        if (!rdOK) {
          deleteAdfSlitHolder(slitHolder);
          cpl_msg_error(modName, "readDoubleDescFromTable returned an error");
          return(NULL);
        }     
      }
      
      break;
      
    }
  } /* of switch */
  
  return(slitHolder);
}


/*
   Extract the reference apertures from the ADF. The function allocates an
   array of *voids of the size required and returns the pointer to this array
   through slitAr. 
   For each slit defined in the ADF, the function creates a VimosAdfRefrSlit
   and stores the pointer to this slit in the array pointed to by *slitAr.
*/

VimosAdfSlitHolder *extractRefsFromADF(VimosTable *adf)
{
  const char           modName[] = "extractRefsFromADF";
  int                  numSlits;
  int                  i;
  double               dValue;
  VimosBool            rdOK;
  VimosAdfType         adfType;
  VimosAdfSlitHolder   *slitHolder;
  VimosAdfSlitHolder   *linkedSlits=NULL;
  VimosAdfSlitHolder   *tSlitHolder;
  VimosAdfRefrSlit     *tRefrSlit;
  
  /* validate input */
  if (adf == NULL) {
    cpl_msg_error(modName, "NULL input table");
    return(NULL);
  }

  /* initialize to undefined ADF type*/
  adfType = VM_ADF_TYPE_UDF;
  /* only work on MOS ADFs */
  if ( !strcmp(adf->name, VM_ADF_MOS) ) {
    adfType = VM_ADF_TYPE_MOS;
  } 
  
  switch (adfType) {
    /* No recognized ADF */
  default: 
    {
      cpl_msg_error(modName, "Unrecognized ADF type");
      return(NULL);
    }
    /* MOS ADF*/
  case VM_ADF_TYPE_MOS : 
    {
      /* read how many slits are defined in the ADF */
      rdOK = readIntDescFromTable(adf, pilTrnGetKeyword("NoRefSlit"), 
                                  &numSlits, NULL);
      if (!rdOK) {
        cpl_msg_error(modName, "Cannot read descriptor %s from table %s", 
                    pilTrnGetKeyword("NoRefSlit"), adf->name);
        return(NULL);
      }

      slitHolder = NULL;
      
      for (i = 0; i < numSlits; i++) {
        /* create slit object and fill up */
        tRefrSlit = newAdfRefrSlit();
        if (tRefrSlit == NULL) {
          cpl_msg_error(modName, "Function newAdfRefrSlit returned a NULL");
          return(NULL);
        }
        tRefrSlit->slitNo = i+1;

        rdOK = readDoubleDescFromTable(adf, 
                   pilTrnGetKeyword("RefSlitX", i+1), &dValue, NULL);
        tRefrSlit->x = (float) dValue;

        rdOK = rdOK && readDoubleDescFromTable(adf, 
                    pilTrnGetKeyword("RefSlitY", i+1), &dValue, NULL);
        tRefrSlit->y = (float) dValue;

        rdOK = rdOK && readDoubleDescFromTable(adf, 
                    pilTrnGetKeyword("SlitRefDimX", i+1), &dValue, NULL);
        tRefrSlit->size = (float) dValue;

        rdOK = rdOK && readDoubleDescFromTable(adf, 
                    pilTrnGetKeyword("SlitRefObjRA", i+1), &dValue, NULL);
        tRefrSlit->objRA = dValue;

        rdOK = rdOK && readDoubleDescFromTable(adf, 
                    pilTrnGetKeyword("SlitRefObjDec", i+1), &dValue, NULL);
        tRefrSlit->objDec = dValue;

        if (i == 0) {
          slitHolder = newAdfSlitHolder();
          if (slitHolder == NULL) {
            cpl_msg_error(modName, "Function newAdfSlitHolder returned a NULL");
            return(NULL);
          }
          linkedSlits = slitHolder;
          linkedSlits->prev = NULL;
        }
        else {
          tSlitHolder = newAdfSlitHolder();
          if (tSlitHolder == NULL) {
            cpl_msg_error(modName, "Function newAdfSlitHolder returned a NULL");
            return(NULL);
          }
          tSlitHolder->prev = linkedSlits;
          linkedSlits = tSlitHolder;
        }
        linkedSlits->slit = (void *) tRefrSlit;
        linkedSlits->slitNo = i+1;
        linkedSlits->next = NULL;
        linkedSlits->slitType = VM_ADF_REFR_SLIT;
        if (!rdOK) {
          deleteAdfSlitHolder(slitHolder);
          cpl_msg_error(modName, "readDoubleDescFromTable returned an error");
          return(NULL);
        }     
      }
    } 
  } /* of switch on ADF type*/

  return(slitHolder);  
}


/* this routine can surely be written more compact...
One day...

Found a bug: it was computing the inverse dispersion polinomial starting
from the optical distorsion matrix, instead of the inverse dispersion
matrix. 99/12/22 (MS)

  23 Oct 01: added IFUfibPeakX in calcSlitLocationsOnCCD (AZ)
  25 Jul 02: added IFUfibTrans in calcSlitLocationsOnCCD (AZ)
  10 Dec 02: define IFU_NUMPIX to always have 5 pixels per fiber spectra (AZ)
*/

VimosBool calcSlitLocationsOnCCD(void *slit, VimosAdfSlitType slitType,
                          VimosDistModel2D *optModX,VimosDistModel2D *optModY,
                          VimosDistModelFull *crvMod, 
                          VimosDistModelFull *invDispMat, 
                          VimosDistModel2D *contModX, 
                          VimosDistModel2D *contModY,
                          VimosFloatArray **ccdX, VimosFloatArray **ccdY, 
                          VimosFloatArray **maskX, VimosFloatArray **maskY, 
                          VimosDistModel1D ***crvPol, 
			  VimosFloatArray **crvPolRms,
                          VimosDistModel1D ***invDis, 
			  VimosFloatArray **invDisRms, float lambda0, 
                          int *numPix, int *IFUslitNo, int *IFUfibNo,
			  float *IFUfibPeakX, float *IFUfibTrans,
                          VimosFloatArray **zeroX, VimosFloatArray **zeroY,
                          VimosIntArray **invDisQuality)
{

#define IFU_NUMPIX 5

  const char modName[] = "calcSlitLocationsOnCCD";
  float xCen, yCen;
  float xSize;
  float xLow, xHig;
  float yLow, yHig;
  float x1CCD, y1CCD;
  float x2CCD, y2CCD;
  float xSlit, ySlit;
  double scale;
  int   i;
  int   beg, ste;
  
  float tFloat;

  double t, dt;

  pilErrno = 0;

  /* switch on type of slit */
  switch (slitType) {
  default :
    {
      /* unrecognized slit type, return */
      *numPix = 0;
      cpl_msg_error(modName, "Unrecognized slit type");
      return(VM_FALSE);
    }
  case VM_ADF_RECT_SLIT :
    {
      /* rectangular slit */

      /* get mask coordinates of slit centre and size */
      xCen = ((VimosAdfRectSlit *)slit)->x;
      yCen = ((VimosAdfRectSlit *)slit)->y;
      xSize = ((VimosAdfRectSlit *)slit)->dimX;
      /* compute edges of slit on mask */
      xLow = xCen-xSize/2.0;
      xHig = xCen+xSize/2.0;
      /* compute CCD position of edges, we assume dispersion is in Y */
      x1CCD = computeDistModel2D(optModX, xLow, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      y1CCD = computeDistModel2D(optModY, xLow, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      x2CCD = computeDistModel2D(optModX, xHig, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      y2CCD = computeDistModel2D(optModY, xHig, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      /* x1CCD should be smallest */
      if (x2CCD < x1CCD) {
        tFloat = x1CCD;
        x1CCD = x2CCD;
        x2CCD = tFloat;
        tFloat = y1CCD;
        y1CCD = y2CCD;
        y2CCD = tFloat;
      }
      /* number of output rows for this slit */
      *numPix = (x2CCD - x1CCD);
      /* get mm/pixel for this slit from model as average of scale in X */
      scale = fabs((double)(xHig-xLow)/(x2CCD-x1CCD));      
    
      /* allocate output arrays */
      *ccdX = newFloatArray(*numPix);
      if (*ccdX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *ccdY = newFloatArray(*numPix);
      if (*ccdY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *maskX = newFloatArray(*numPix);
      if (*maskX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *maskY = newFloatArray(*numPix);
      if (*maskY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *zeroX = newFloatArray(*numPix);
      if (*zeroX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *zeroY = newFloatArray(*numPix);
      if (*zeroY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *crvPolRms = newFloatArray(*numPix);
      if (*crvPolRms == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *invDisRms = newFloatArray(*numPix);
      if (*invDisRms == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *invDisQuality = newIntArray(*numPix);
      if (*invDisQuality == NULL) {
        cpl_msg_error(modName, "newIntArray has returned NULL");
        return(VM_FALSE);
      }

      *crvPol=(VimosDistModel1D **)cpl_malloc(*numPix*sizeof(VimosDistModel1D *));
      /* check if space was allocated */
      if (*crvPol == NULL) {
        cpl_msg_error(modName, "Allocation Error");
        return(VM_FALSE);
      }

      *invDis=(VimosDistModel1D **)cpl_malloc(*numPix*sizeof(VimosDistModel1D *));
      /* check if space was allocated */
      if (*invDis == NULL) {
        cpl_msg_error(modName, "Allocation Error");
        return(VM_FALSE);
      }      

      /* initialize */
      xSlit = xLow;
      ySlit = yCen;
      
      /* loop to calculate position of slit for central wavelength */
      for (i = 0; i < *numPix; i++) {
        (*maskX)->data[i] = xSlit;
        (*maskY)->data[i] = ySlit;
        (*ccdX)->data[i] = computeDistModel2D(optModX, xSlit, ySlit);
        if (pilErrno) {
          cpl_msg_error(modName, "computeDistModel2D has returned an error");
          return(VM_FALSE);
        }

        (*ccdY)->data[i] = computeDistModel2D(optModY, xSlit, ySlit);
        if (pilErrno) {
          cpl_msg_error(modName, "computeDistModel2D has returned an error");
          return(VM_FALSE);
        }

        if (contModX)
        {
          (*zeroX)->data[i] = computeDistModel2D(contModX, xSlit, ySlit) -
                              (*ccdX)->data[i] ;
          if (pilErrno) {
            cpl_msg_error(modName, "computeDistModel2D has returned an error");
            return(VM_FALSE);
          }
        }
        if (contModY)
        {
          (*zeroY)->data[i] = computeDistModel2D(contModY, xSlit, ySlit) -
                            (*ccdY)->data[i];
          if (pilErrno) {
            cpl_msg_error(modName, "computeDistModel2D has returned an error");
            return(VM_FALSE);
          }
        }

        if (!getDistModel1DFromFull(crvMod, xSlit, ySlit, &((*crvPol)[i]))) {
          cpl_msg_error(modName, "getDistModel1DFromFull has returned an error");
          return(VM_FALSE);
        }

        (*crvPol)[i]->offset = (*ccdY)->data[i];

        if (!getDistModel1DFromFull(invDispMat, xSlit, ySlit, &((*invDis)[i])))
          {
            cpl_msg_error(modName, "getDistModel1DFromFull returned an error");
            return(VM_FALSE);
          }

        (*invDis)[i]->offset = lambda0;
        (*invDisQuality)->data[i] = 1;
	(*crvPolRms)->data[i] = 0.;
	(*invDisRms)->data[i] = 0.;
                
        xSlit += scale;
      }
      *IFUslitNo = 0;
      *IFUfibNo = 0;
      /* added (AZ) */
      *IFUfibPeakX = 0.;
      *IFUfibTrans = 0.0;

      
      break;
    }
  case VM_ADF_CURV_SLIT :
    {
     /* curved slit */

      xLow = ((VimosAdfCurvSlit *)slit)->xMiddle->x0;
      xHig = computeBezierCurve(((VimosAdfCurvSlit *)slit)->xMiddle, 1.0);
      if (pilErrno) {
        cpl_msg_error(modName, "computeBezierCurve has returned an error");
        return(VM_FALSE);
      }

      yLow = ((VimosAdfCurvSlit *)slit)->yMiddle->x0;
      yHig = computeBezierCurve(((VimosAdfCurvSlit *)slit)->yMiddle, 1.0);
      if (pilErrno) {
        cpl_msg_error(modName, "computeBezierCurve has returned an error");
        return(VM_FALSE);
      }
      
      /* compute CCD position of edges, we assume dispersion is in Y */
      x1CCD = computeDistModel2D(optModX, xLow, yLow);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      y1CCD = computeDistModel2D(optModY, xLow, yLow);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      x2CCD = computeDistModel2D(optModX, xHig, yHig);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      y2CCD = computeDistModel2D(optModY, xHig, yHig);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      /* x1CCD show be smallest */
      if (x2CCD < x1CCD) {
        tFloat = x1CCD;
        x1CCD = x2CCD;
        x2CCD = tFloat;
        tFloat = y1CCD;
        y1CCD = y2CCD;
        y2CCD = tFloat;
      }
      /* number of output rows for this slit */
      *numPix = (x2CCD - x1CCD);
      /* get mm/pixel for this slit from model as average of scale in X */
      scale = fabs((xHig-xLow)/(x2CCD-x1CCD));

      /* allocate output arrays */
      *ccdX = newFloatArray(*numPix);
      if (*ccdX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *ccdY = newFloatArray(*numPix);
      if (*ccdY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *maskX = newFloatArray(*numPix);
      if (*maskX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *maskY = newFloatArray(*numPix);
      if (*maskY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *zeroX = newFloatArray(*numPix);
      if (*zeroX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *zeroY = newFloatArray(*numPix);
      if (*zeroY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *crvPolRms = newFloatArray(*numPix);
      if (*crvPolRms == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *invDisRms = newFloatArray(*numPix);
      if (*invDisRms == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *invDisQuality = newIntArray(*numPix);
      if (*invDisQuality == NULL) {
        cpl_msg_error(modName, "newIntArray has returned NULL");
        return(VM_FALSE);
      }

      *crvPol=(VimosDistModel1D **)cpl_malloc(*numPix*sizeof(VimosDistModel1D *));
      /* check if space was allocated */
      if (*crvPol == NULL) {
        cpl_msg_error(modName, "Allocation Error");
        return(VM_FALSE);
      }

      *invDis=(VimosDistModel1D **)cpl_malloc(*numPix*sizeof(VimosDistModel1D *));
      /* check if space was allocated */
      if (*invDis == NULL) {
        cpl_msg_error(modName, "Allocation Error");
        return(VM_FALSE);
      }

      /* initialize */
      xSlit = xLow;
      ySlit = yLow;
      if (xLow > xHig) {
        t = 1.0;
        beg = *numPix - 1;
        ste = -1;
      }
      else {
        t = 0.0;
        beg = 0;
        ste = 1;
      }
      
      /* loop to calculate position of slit for central wavelength */
      for (i = beg; i < *numPix && i >= 0; i += ste) {

        /* we need some way to step through the bezier curves such that we
           step in X in steps of 1 pixel.... */

        (*maskX)->data[i] = xSlit;
        (*maskY)->data[i] = ySlit;
        (*ccdX)->data[i] = computeDistModel2D(optModX, xSlit, ySlit);
        if (pilErrno) {
          cpl_msg_error(modName, "computeDistModel2D has returned an error");
          return(VM_FALSE);
        }
        (*ccdY)->data[i] = computeDistModel2D(optModY, xSlit, ySlit);
	if (pilErrno) {
	  cpl_msg_error(modName, "computeDistModel2D has returned an error");
	  return(VM_FALSE);
	}
	if (contModX)
	{
	  (*zeroX)->data[i] = computeDistModel2D(contModX, xSlit, ySlit) -
	    (*ccdX)->data[i] ;
	  if (pilErrno) {
	    cpl_msg_error(modName, "computeDistModel2D has returned an error");
	    return(VM_FALSE);
	  }
	}
	if (contModY)
	{
	  (*zeroY)->data[i] = computeDistModel2D(contModY, xSlit, ySlit) -
	    (*ccdY)->data[i];
	  if (pilErrno) {
	    cpl_msg_error(modName, "computeDistModel2D has returned an error");
	    return(VM_FALSE);
	  }
	}
        if (!getDistModel1DFromFull(crvMod, xSlit, ySlit, &((*crvPol)[i]))) {
          cpl_msg_error(modName, "getDistModel1DFromFull has returned an error");
          return(VM_FALSE);
        }
        (*crvPol)[i]->offset = (*ccdY)->data[i];
        if (!getDistModel1DFromFull(invDispMat, xSlit, ySlit, &((*invDis)[i])))
        {
            cpl_msg_error(modName, "getDistModel1DFromFull returned an error");
            return(VM_FALSE);
        }
        (*invDis)[i]->offset = lambda0;
        (*invDisQuality)->data[i] = 1;
	(*crvPolRms)->data[i] = 0.;
	(*invDisRms)->data[i] = 0.;

        xSlit = xLow + scale * (i + ste - beg);
        dt = scale / (
                      ((VimosAdfCurvSlit *)slit)->xMiddle->a * 3 * t * t +
                      ((VimosAdfCurvSlit *)slit)->xMiddle->b * 2 * t +
                      ((VimosAdfCurvSlit *)slit)->xMiddle->c
                     );

        t += dt;

        xSlit = computeBezierCurve(((VimosAdfCurvSlit *)slit)->xMiddle, t);
        ySlit = computeBezierCurve(((VimosAdfCurvSlit *)slit)->yMiddle, t);

      }
      *IFUslitNo = 0;
      *IFUfibNo = 0;
      /* added (AZ) */
      *IFUfibPeakX = 0.;
      *IFUfibTrans = 0.0;

      break;
    }
  case VM_ADF_CIRC_SLIT :
    {
     /* circular slit */

      /* get mask coordinates of slit centre and size */
      xCen = ((VimosAdfCircSlit *)slit)->x;
      yCen = ((VimosAdfCircSlit *)slit)->y;
      xSize = 2.0 * ((VimosAdfCircSlit *)slit)->radius;
      /* compute edges of slit on mask */
      xLow = xCen-xSize/2.0;
      xHig = xCen+xSize/2.0;
      /* compute CCD position of edges, we assume dispersion is in Y */
      x1CCD = computeDistModel2D(optModX, xLow, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      y1CCD = computeDistModel2D(optModY, xLow, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      x2CCD = computeDistModel2D(optModX, xHig, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      y2CCD = computeDistModel2D(optModY, xHig, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      /* x1CCD show be smallest */
      if (x2CCD < x1CCD) {
        tFloat = x1CCD;
        x1CCD = x2CCD;
        x2CCD = tFloat;
        tFloat = y1CCD;
        y1CCD = y2CCD;
        y2CCD = tFloat;
      }
      /* number of output rows for this slit */

      *numPix = (x2CCD - x1CCD);


      /* ALEX: if IFU, set numPix to 5 to prevent rounding problems (e.g.
	 4.999 pixels => 4 rows per slit) */

      if ((((VimosAdfCircSlit *)slit)->IFUfibNo != 0) &&
	  (((VimosAdfCircSlit *)slit)->IFUslitNo != 0))
	*numPix = IFU_NUMPIX;


      /* get mm/pixel for this slit from model as average of scale in X */
      scale = fabs((xHig-xLow)/(x2CCD-x1CCD));

      /* allocate output arrays */
      *ccdX = newFloatArray(*numPix);
      if (*ccdX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *ccdY = newFloatArray(*numPix);
      if (*ccdY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *maskX = newFloatArray(*numPix);
      if (*maskX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *maskY = newFloatArray(*numPix);
      if (*maskY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *zeroX = newFloatArray(*numPix);
      if (*zeroX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *zeroY = newFloatArray(*numPix);
      if (*zeroY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *crvPolRms = newFloatArray(*numPix);
      if (*crvPolRms == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *invDisRms = newFloatArray(*numPix);
      if (*invDisRms == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *invDisQuality = newIntArray(*numPix);
      if (*invDisQuality == NULL) {
        cpl_msg_error(modName, "newIntArray has returned NULL");
        return(VM_FALSE);
      }

      *crvPol=(VimosDistModel1D **)cpl_malloc(*numPix*sizeof(VimosDistModel1D *));
      /* check if space was allocated */
      if (*crvPol == NULL) {
        cpl_msg_error(modName, "Allocation Error");
        return(VM_FALSE);
      }
      *invDis=(VimosDistModel1D **)cpl_malloc(*numPix*sizeof(VimosDistModel1D *));
      /* check if space was allocated */
      if (*invDis == NULL) {
        cpl_msg_error(modName, "Allocation Error");
        return(VM_FALSE);
      }

      /* initialize */
      xSlit = xLow;
      ySlit = yCen;
      
      /* loop to calculate position of slit for central wavelength */
      for (i = 0; i < *numPix; i++) {
        (*maskX)->data[i] = xSlit;
        (*maskY)->data[i] = ySlit;
        (*ccdX)->data[i] = computeDistModel2D(optModX, xSlit, ySlit);
        if (pilErrno) {
          cpl_msg_error(modName, "computeDistModel2D has returned an error");
          return(VM_FALSE);
        }
        (*ccdY)->data[i] = computeDistModel2D(optModY, xSlit, ySlit);
	if (pilErrno) {
	  cpl_msg_error(modName, "computeDistModel2D has returned an error");
	  return(VM_FALSE);
	}
	if (contModX)
	{
	  (*zeroX)->data[i] = computeDistModel2D(contModX, xSlit, ySlit);
	  if (pilErrno) {
	    cpl_msg_error(modName, "computeDistModel2D has returned an error");
	    return(VM_FALSE);
	  }
	}
	if (contModY)
	{
	  (*zeroY)->data[i] = computeDistModel2D(contModY, xSlit, ySlit);
	  if (pilErrno) {
	    cpl_msg_error(modName, "computeDistModel2D has returned an error");
	    return(VM_FALSE);
	  }
	}
        if (!getDistModel1DFromFull(crvMod, xSlit, ySlit, &((*crvPol)[i]))) {
          cpl_msg_error(modName, "getDistModel1DFromFull has returned an error");
          return(VM_FALSE);
        }
        (*crvPol)[i]->offset = (*ccdY)->data[i];
        if (!getDistModel1DFromFull(invDispMat, xSlit, ySlit, &((*invDis)[i])))
          {
            cpl_msg_error(modName, "getDistModel1DFromFull returned an error");
            return(VM_FALSE);
          }
        (*invDis)[i]->offset = lambda0;
        (*invDisQuality)->data[i] = 1;
	(*crvPolRms)->data[i] = 0.;
	(*invDisRms)->data[i] = 0.;

        xSlit += scale;
      }
      *IFUslitNo = ((VimosAdfCircSlit *)slit)->IFUslitNo;
      *IFUfibNo = ((VimosAdfCircSlit *)slit)->IFUfibNo;
      *IFUfibTrans = ((VimosAdfCircSlit *)slit)->IFUfibTrans;

      /* 
       * ALEX: TO BE FINISHED: set this value to fiber central pixel
                               (to be modified in the future)
       */
      *IFUfibPeakX = (*ccdX)->data[2];

      
      break;
    }
  case VM_ADF_REFR_SLIT :
    {
     /* square reference slit slit */

      /* get mask coordinates of slit centre and size */
      xCen = ((VimosAdfRefrSlit *)slit)->x;
      yCen = ((VimosAdfRefrSlit *)slit)->y;
      xSize = ((VimosAdfRefrSlit *)slit)->size;
      /* compute edges of slit on mask */
      xLow = xCen-xSize/2.0;
      xHig = xCen+xSize/2.0;
      /* compute CCD position of edges, we assume dispersion is in Y */
      x1CCD = computeDistModel2D(optModX, xLow, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      y1CCD = computeDistModel2D(optModY, xLow, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      x2CCD = computeDistModel2D(optModX, xHig, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      y2CCD = computeDistModel2D(optModY, xHig, yCen);
      if (pilErrno) {
        cpl_msg_error(modName, "computeDistModel2D has returned an error");
        return(VM_FALSE);
      }
      /* x1CCD show be smallest */
      if (x2CCD < x1CCD) {
        tFloat = x1CCD;
        x1CCD = x2CCD;
        x2CCD = tFloat;
        tFloat = y1CCD;
        y1CCD = y2CCD;
        y2CCD = tFloat;
      }
      /* number of output rows for this slit */
      *numPix = (x2CCD - x1CCD);
      /* get mm/pixel for this slit from model as average of scale in X */
      scale = fabs((xHig-xLow)/(x2CCD-x1CCD));

      /* allocate output arrays */
      *ccdX = newFloatArray(*numPix);
      if (*ccdX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *ccdY = newFloatArray(*numPix);
      if (*ccdY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *maskX = newFloatArray(*numPix);
      if (*maskX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *maskY = newFloatArray(*numPix);
      if (*maskY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *zeroX = newFloatArray(*numPix);
      if (*zeroX == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *zeroY = newFloatArray(*numPix);
      if (*zeroY == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *crvPolRms = newFloatArray(*numPix);
      if (*crvPolRms == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *invDisRms = newFloatArray(*numPix);
      if (*invDisRms == NULL) {
        cpl_msg_error(modName, "newFloatArray has returned NULL");
        return(VM_FALSE);
      }
      *invDisQuality = newIntArray(*numPix);
      if (*invDisQuality == NULL) {
        cpl_msg_error(modName, "newIntArray has returned NULL");
        return(VM_FALSE);
      }

      *crvPol=(VimosDistModel1D **)cpl_malloc(*numPix*sizeof(VimosDistModel1D *));
      /* check if space was allocated */
      if (*crvPol == NULL) {
        cpl_msg_error(modName, "Allocation Error");
        return(VM_FALSE);
      }

      *invDis=(VimosDistModel1D **)cpl_malloc(*numPix*sizeof(VimosDistModel1D *));
      /* check if space was allocated */
      if (*invDis == NULL) {
        cpl_msg_error(modName, "Allocation Error");
        return(VM_FALSE);
      }

      /* initialize */
      xSlit = xLow;
      ySlit = yCen;
      
      /* loop to calculate position of slit for central wavelength */
      for (i = 0; i < *numPix; i++) {
        (*maskX)->data[i] = xSlit;
        (*maskY)->data[i] = ySlit;
        (*ccdX)->data[i] = computeDistModel2D(optModX, xSlit, ySlit);
        if (pilErrno) {
          cpl_msg_error(modName, "computeDistModel2D has returned an error");
          return(VM_FALSE);
        }
        (*ccdY)->data[i] = computeDistModel2D(optModY, xSlit, ySlit);
	if (pilErrno) {
	  cpl_msg_error(modName, "computeDistModel2D has returned an error");
	  return(VM_FALSE);
	}
	if (contModX)
	{
	  (*zeroX)->data[i] = computeDistModel2D(contModX, xSlit, ySlit) -
	    (*ccdX)->data[i] ;
	  if (pilErrno) {
	    cpl_msg_error(modName, "computeDistModel2D has returned an error");
	    return(VM_FALSE);
	  }
	}
	if (contModY)
	{
	  (*zeroY)->data[i] = computeDistModel2D(contModY, xSlit, ySlit) -
	    (*ccdY)->data[i];
	  if (pilErrno) {
	    cpl_msg_error(modName, "computeDistModel2D has returned an error");
	    return(VM_FALSE);
	  }
	}
        if (!getDistModel1DFromFull(crvMod, xSlit, ySlit, &((*crvPol)[i]))) {
          cpl_msg_error(modName, "getDistModel1DFromFull has returned an error");
          return(VM_FALSE);
        }
        (*crvPol)[i]->offset = (*ccdY)->data[i];
        if (!getDistModel1DFromFull(invDispMat, xSlit, ySlit, &((*invDis)[i])))
          {
            cpl_msg_error(modName, "getDistModel1DFromFull returned an error");
            return(VM_FALSE);
          }
        (*invDis)[i]->offset = lambda0;
        (*invDisQuality)->data[i] = 1;
	(*crvPolRms)->data[i] = 0.;
	(*invDisRms)->data[i] = 0.;

        xSlit += scale;
      }
      *IFUslitNo = 0;
      *IFUfibNo = 0;
      /* added (AZ) */
      *IFUfibPeakX = 0.;
      *IFUfibTrans = 0.0;

      break;
    }
  } /* of switch (slitType) */

  return(VM_TRUE);  

#undef IFU_NUMPIX

}


VimosAdfType getADFTypeFromDesc(VimosDescriptor *desc)
{
  const char modName[] = "getADFTypeFromDesc";
  VimosDescriptor *tDesc;
  
  tDesc = findDescriptor(desc, "ESO INS ADF TYPE");
  
  if (tDesc == NULL) {
    cpl_msg_error(modName, "Cannot find descriptor %s", "ESO INS ADF TYPE");
    return(VM_ADF_TYPE_UDF);
  } else {
    if (!strncmp("MOS", tDesc->descValue->s, 3)) {
      return(VM_ADF_TYPE_MOS);
    }
    if (!strncmp("IFU", tDesc->descValue->s, 3)) {
      return(VM_ADF_TYPE_IFU);
    }
    if (!strncmp("IMAGE", tDesc->descValue->s, 5)) {
      return(VM_ADF_TYPE_IMA);
    }
  }

  return(VM_ADF_TYPE_UDF);
}


/*
   Read the ADF from the header of the Fits Image. The type of the ADF should 
   match that defined in the image.
 */
VimosBool readFitsADF(VimosTable *adf, VimosImage *adfImage)
{
  char descVal[81];

  /* validate input */
  if (adf == NULL) {
    cpl_msg_error("readFitsADF","NULL input table");
    return(VM_FALSE);
  }
    
  if (!readDescsFromFitsImage(&(adf->descs), adfImage)) {
    cpl_msg_error("readFitsADF","The function readDescsFromFitsImage has returned an error");
    return(VM_FALSE);
  }

  if (!readStringDescriptor(adf->descs, "ESO INS ADF TYPE", descVal, NULL)) 
    {
      cpl_msg_error("readFitsADF", "The function readStringDescriptor has "
                  "returned an error");
      return(VM_FALSE);
    }

  /* make sure ADF types are ok */
  if ( !strncmp(descVal, "MOS", 3) ) {
    strcpy(adf->name, VM_ADF_MOS);
  }
  if ( !strncmp(descVal, "IFU", 3) ) {
    strcpy(adf->name, VM_ADF_IFU);
  }
  if ( !strncmp(descVal, "IMA", 3) ) {
    strcpy(adf->name, VM_ADF_IMA);
  }
        
  return(VM_TRUE);
}
