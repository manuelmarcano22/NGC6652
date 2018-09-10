/* $Id: vmdistmodels.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <math.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <piltranslator.h>
#include <pilerrno.h>

#include "vmtypes.h"
#include "vmfit.h"
#include "vmdistmodels.h"
#include "cpl.h"


/* Constructors and destructors for the models */

VimosDistModel1D *newDistModel1D(const int order)
{
  const char modName[] = "newDistModel1D";
  int i;
  VimosDistModel1D *tmpMod;
  
  /* check order, should be at least zero */
  if (order < 0) {
    cpl_msg_error(modName, "Invalid input order");
    return(NULL);
  }
  
  /* allocated space and check */
  tmpMod = (VimosDistModel1D *) cpl_malloc(sizeof(VimosDistModel1D));
  if (tmpMod == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  
  /* allocate array for coefficients  */
  tmpMod->coefs = (double *) cpl_calloc(order+1, sizeof(double));

  /* check allocation */
  if (tmpMod->coefs == NULL) {
    /* cleanup */
    cpl_free(tmpMod);
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  /* store order */
  tmpMod->order = order;
  tmpMod->offset = 0.0;
  for (i=0; i <= order; i++) {
    tmpMod->coefs[i] = 0.0;
  }
  
  return(tmpMod);
  
}

VimosDistModel2D *newDistModel2D(const int orderX, const int orderY)
{
  const char modName[] = "newDistModel2D";
  VimosDistModel2D *tmpMod;
  int i, j;
  
  /* check order, should be at least zero */
  if ( (orderX < 0) || (orderY < 0) ) {
    cpl_msg_error(modName, "Invalid input order (X or Y)");
    return(NULL);
  }
  
  /* allocated space and check */
  tmpMod = (VimosDistModel2D *) cpl_malloc(sizeof(VimosDistModel2D));
  if (tmpMod == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  /* allocate array for coefficients  */
  tmpMod->coefs = (double **) cpl_calloc(orderX+1, sizeof(double *));
  
  /* check allocation */
  if (tmpMod->coefs == NULL) {
    /* cleanup */
    cpl_free(tmpMod);
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  for (i=0; i <= orderX; i++) {
    tmpMod->coefs[i] = (double *) cpl_calloc(orderY+1, sizeof(double));

    /* check allocation */
    if (tmpMod->coefs[i] == NULL) {
      /* cleanup */
      cpl_free(tmpMod);
      cpl_msg_error(modName, "Allocation Error");
      return(NULL);
    }
  }
  
  /* store orders */
  tmpMod->orderX = orderX;
  tmpMod->orderY = orderY;
  tmpMod->offsetX = 0.0;
  tmpMod->offsetY = 0.0;
  for (i=0; i <= orderX; i++) {
    for (j=0; j<= orderY; j++) {
      tmpMod->coefs[i][j] = 0.0;
    }
  }
  
  return(tmpMod);
  
}

VimosDistModelFull *newDistModelFull(const int orderPol,
                                     const int orderX, const int orderY)
{
  const char modName[] = "newDistModelFull";
  int i;
  VimosDistModelFull *tmpMod;
  
  /* check order, should be at least zero */
  if ( (orderPol < 0) || (orderX < 0) || (orderY < 0) ) {
    cpl_msg_error(modName, "Invalid input order (polynomial, X or Y)");
    return(NULL);
  }
  
  /* allocated space and check */
  tmpMod = (VimosDistModelFull *) cpl_malloc(sizeof(VimosDistModelFull));
  if (tmpMod == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }

  /* allocated space and check */
  tmpMod->coefs = 
    (VimosDistModel2D **) cpl_malloc((orderPol+1)*sizeof(VimosDistModel2D *));
  if (tmpMod->coefs == NULL) {
    /* cleanup */
    cpl_free(tmpMod);
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  
  /* allocate array for coefficients  */
  for (i = 0; i <= orderPol; i++ ) {
    tmpMod->coefs[i] = newDistModel2D(orderX, orderY);
    if (tmpMod->coefs[i] == NULL) {
      /* cleanup */
      cpl_free(tmpMod);
      cpl_msg_error(modName, "The function newDistModel2D has returned NULL");
      return(NULL);
    }
  }

  /* store orders */
  tmpMod->orderPol = orderPol;
  tmpMod->orderX   = orderX;
  tmpMod->orderY   = orderY;
  tmpMod->offsetT  = 0.0;
  tmpMod->offsetX  = 0.0;
  tmpMod->offsetY  = 0.0;

  return(tmpMod);
}

void deleteDistModel1D(VimosDistModel1D *mod)
{
  /* check input */
  if (mod == NULL) {
    return;
  }
  /* delete the array with coefficients */
  cpl_free(mod->coefs);
  /* delete sructure */
  cpl_free(mod);
  
}

void deleteDistModel2D(VimosDistModel2D *mod)
{
  int i;
  
  /* check input */
  if (mod == NULL) {
    return;
  }
  /* delete coefficients */
  for (i=0; i <= mod->orderX; i++) {
    cpl_free(mod->coefs[i]);
  }
  cpl_free(mod->coefs);
  
  /* delete structure */
  cpl_free(mod);
}

void deleteDistModelFull(VimosDistModelFull *mod)
{
  int i;
  
  /* check input */
  if (mod == NULL) {
    return;
  }
  /* delete coefficients */
  for (i=0; i <= mod->orderPol; i++){
    deleteDistModel2D(mod->coefs[i]);
  }
  /* delete structure */
  cpl_free(mod);

}


/* compute model values */
double computeDistModel1D(VimosDistModel1D *mod, float x)
{
  const char modName[] = "computeDistModel1D";
  int i;
  double retval;
  double tmpX;
  
  pilErrno = 0;

  /* check input */
  if (mod == NULL) {
    cpl_msg_error(modName, "NULL input pointer");
    pilErrno = 1;
    return(0.0);
  }
  
  /* initialize */
  retval = 0.0;
  tmpX = 1.0;

  /* loop over coefficients and calculate value */
  for (i=0 ; i <= mod->order; i++) {
    retval += mod->coefs[i]*tmpX;
    /* increase order to evaluate next */
    tmpX *= ((double) x - mod->offset);
  }
  /* return result */
  return(retval);
      
}

double computeDistModel2D(VimosDistModel2D *mod, float x, float y)
{
  const char modName[] = "computeDistModel2D";
  int i,j;
  double retval;
  double tmpX;
  double tmpY;
  
  pilErrno = 0;

  /* check input */
  if (mod == NULL) {
    /* if bad input: we have to return something... */
    cpl_msg_error(modName, "NULL input pointer");
    pilErrno = 1;
    return(0.0);
  }
  
  /* initialize */
  retval = 0.0;
  tmpX = 1.0;

  /* loop over coefficients and calculate value */
  for (i=0; i <= mod->orderX; i++) {
    tmpY = 1.0;
    for (j=0; j <= mod->orderY; j++) {
      retval += mod->coefs[i][j]*tmpX*tmpY;
      /* increase order to evaluate next */
      tmpY *= ((double) y - mod->offsetY);
    }
    /* increase order to evaluate next */
    tmpX *= ((double) x - mod->offsetX);
  }

  /* return result */
  return(retval);
}

double computeDistModelFull(VimosDistModelFull *mod, float t, float x, float y)
{
  const char modName[] = "computeDistModelFull";
  int i;
  double retval;
  double tmpT;
  
  pilErrno = 0;

  /* check input */
  if (mod == NULL) {
    cpl_msg_error(modName, "NULL input pointer");
    pilErrno = 1;
    return(0.0);
  }
  
  /* initialize */
  retval = 0.0;
  tmpT = 1.0;

  /* loop over coefficients and compute value */
  for (i=0; i <= mod->orderPol; i++ ) {
    retval += computeDistModel2D(mod->coefs[i], x, y)*tmpT;
    if (pilErrno) {
      cpl_msg_error(modName, "Function computeDistModel2D returned an error");
      pilErrno = 1;
      return(0.0);
    }
    /* increase order to evaluate next */
    tmpT *= ((double) t - mod->offsetT);
  }
  /* return result */
  return(retval);
}


VimosBool getDistModel1DFromFull(VimosDistModelFull *full, float x, float y,
                                 VimosDistModel1D **mod1D) 
{
  const char modName[] = "getDistModel1DFromFull";
  int i;
  
  pilErrno = 0;

  *mod1D = newDistModel1D(full->orderPol);

  if (*mod1D == NULL) {
    cpl_msg_error(modName, "The function newDistModel1D has returned NULL");
    return(VM_FALSE);
  }

  /* loop over coefficients and compute value */
  for (i=0; i <= full->orderPol; i++ ) {
    (*mod1D)->coefs[i] = computeDistModel2D(full->coefs[i], x, y);
    if (pilErrno) {
      deleteDistModel1D(*mod1D);
      cpl_msg_error(modName, "Function computeDistModel2D returned an error");
      return(VM_FALSE);
    }
  }
  (*mod1D)->offset = full->offsetT;

  return(VM_TRUE);
}


VimosBool readOptDistModel(VimosDescriptor *desc, VimosDistModel2D **optModX,
                           VimosDistModel2D **optModY)
{
  const char modName[] = "readOptDistModel";
  int        iValue;
  double     dValue;
  int        i, j;
  
  
  *optModX = NULL;
  *optModY = NULL;

 /*
  *  Note: the optical distorsion model is made of two 2D distorsion
  *  models, one for the output X coordinate, and the other for the
  *  output Y coordinate. Both models have an equal "order" for the
  *  input x and y. Note also that when an Optical Distorsion Model
  *  is read, the offsets are always left to zero (as in the header
  *  no corresponding keyword is present). This implies that offsets
  *  are _not_ foreseen for the DRS optical distorsion models.  (C.Izzo)
  */

  if (readIntDescriptor(desc, pilTrnGetKeyword("OptDistOrdX"), &iValue, 
                        NULL) == VM_FALSE) {
    cpl_msg_error(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("OptDistOrdX"));
    return VM_TRUE;
  }
  
  *optModX = newDistModel2D(iValue, iValue);
  if (*optModX == NULL) {
    cpl_msg_error(modName, "The function newDistModel2D has returned NULL");
    return(VM_FALSE);
  }
  for (i = 0; i <= iValue; i++) {
    for (j = 0; j <= iValue; j++) {
      if (readDoubleDescriptor(desc, pilTrnGetKeyword("OptDistX", i, j), 
                               &dValue, NULL) == VM_FALSE) {
        deleteDistModel2D(*optModX);
        *optModX = NULL;
        cpl_msg_error(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("OptDistX", i, j)); 
        return VM_FALSE;
      }
      else {
        (*optModX)->coefs[i][j] = dValue;
      }
    }
  }

  if (readIntDescriptor(desc, pilTrnGetKeyword("OptDistOrdY"), 
                        &iValue, NULL) == VM_FALSE) {
    deleteDistModel2D(*optModX);
    *optModX = NULL;
    cpl_msg_error(modName, "Cannot read descriptor %s", 
                pilTrnGetKeyword("OptDistOrdY"));
    return VM_FALSE;
  }
  
  *optModY = newDistModel2D(iValue, iValue);
  if (*optModY == NULL) {
    cpl_msg_error(modName, "The function newDistModel2D has returned NULL");
    return VM_FALSE;
  }  

  for (i = 0; i <= iValue; i++) {
    for (j = 0; j <= iValue; j++) {
      if (readDoubleDescriptor(desc, pilTrnGetKeyword("OptDistY", i, j), 
                               &dValue, NULL) == VM_FALSE) {
        deleteDistModel2D(*optModX);
        deleteDistModel2D(*optModY);
        *optModX = NULL;
        *optModY = NULL;
        cpl_msg_error(modName, "Cannot read descriptor %s", 
                    pilTrnGetKeyword("OptDistY", i, j));
        return VM_FALSE;
      }
      else {
        (*optModY)->coefs[i][j] = dValue;
      }
    }
  }
  
  return VM_TRUE;
}

VimosBool writeOptDistModel(VimosDescriptor **desc, VimosDistModel2D *optModX,
                            VimosDistModel2D *optModY)
{
  const char modName[] = "writeOptDistModel";
  int        i, j;

 /*
  *  Note: the optical distorsion model is made of two 2D distorsion
  *  models, one for the output X coordinate, and the other for the
  *  output Y coordinate. Both models have an equal "order" for the
  *  input x and y, and this is why just optModX/Y->orderX is written
  *  to header. Note also that when an Optical Distorsion Model is 
  *  written to file, any offset is ignored (in the header no corresponding 
  *  keyword is created). This implies that offsets are _not_ foreseen 
  *  for the DRS optical distorsion models.  (C.Izzo)
  */
  
  if (optModX) {
    if (writeIntDescriptor(desc, pilTrnGetKeyword("OptDistOrdX"), 
                           optModX->orderX, "") == VM_FALSE) {
      cpl_msg_error(modName, "Descriptor %s not found", 
                  pilTrnGetKeyword("OptDistOrdX"));
      return VM_FALSE;
    }
  
    for (i = 0; i <= optModX->orderX; i++) {
      for (j = 0; j <= optModX->orderY; j++) {
        if (writeDoubleDescriptor(desc, pilTrnGetKeyword("OptDistX", i, j), 
                                  optModX->coefs[i][j], "") == VM_FALSE) {
          cpl_msg_error(modName, "Descriptor %s not found", 
                      pilTrnGetKeyword("OptDistX", i, j));
          return VM_FALSE;
        }
      }
    }
  }
  
  if (optModY) {
    if (writeIntDescriptor(desc, pilTrnGetKeyword("OptDistOrdY"), 
                           optModY->orderX, "") == VM_FALSE){
      cpl_msg_error(modName, "Descriptor %s not found", 
                  pilTrnGetKeyword("OptDistOrdY"));
      return VM_FALSE;
    }
  
    for (i = 0; i <= optModY->orderX; i++) {
      for (j = 0; j <= optModY->orderY; j++) {
        if (writeDoubleDescriptor(desc, pilTrnGetKeyword("OptDistY", i, j), 
                                  optModY->coefs[i][j], "") == VM_FALSE) {
          cpl_msg_error(modName, "Descriptor %s not found", 
                      pilTrnGetKeyword("OptDistY", i, j));
          return VM_FALSE;
        }
      }
    }
  }

  return VM_TRUE;
}


/*
 * The following ludicrous function is a copy of the writeOptDistModel()
 * (above) that has been made necessary by a bug that is still present
 * in the OSLX: when writing any distorsion model to an image header,
 * its coefficients must be character strings. When one day this problem
 * will be solved, this function should disappear, and any call to it
 * should be replaced by a call to writeOptDistModel().  (C.Izzo)
 */

VimosBool writeOptDistModelString(VimosDescriptor **desc, 
                                  VimosDistModel2D *optModX,
                                  VimosDistModel2D *optModY)
{
  const char modName[] = "writeOptDistModelString";
  int        i, j;
  char       string[80];

 /*
  *  Note: the optical distorsion model is made of two 2D distorsion
  *  models, one for the output X coordinate, and the other for the
  *  output Y coordinate. Both models have an equal "order" for the
  *  input x and y, and this is why just optModX/Y->orderX is written
  *  to header. Note also that when an Optical Distorsion Model is
  *  written to file, any offset is ignored (in the header no corresponding
  *  keyword is created). This implies that offsets are _not_ foreseen
  *  for the DRS optical distorsion models.  (C.Izzo)
  */

  if (writeIntDescriptor(desc, pilTrnGetKeyword("OptDistOrdX"),
                         optModX->orderX, "") == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found",
                pilTrnGetKeyword("OptDistOrdX"));
    return VM_FALSE;
  }

  for (i = 0; i <= optModX->orderX; i++) {
    for (j = 0; j <= optModX->orderY; j++) {
      sprintf(string, "%#.14E", optModX->coefs[i][j]);
      if (writeStringDescriptor(desc, pilTrnGetKeyword("OptDistX", i, j),
                                string, "") == VM_FALSE) {
        cpl_msg_error(modName, "Descriptor %s not found",
                    pilTrnGetKeyword("OptDistX", i, j));
        return VM_FALSE;
      }
    }
  }

  if (writeIntDescriptor(desc, pilTrnGetKeyword("OptDistOrdY"),
                         optModY->orderX, "") == VM_FALSE){
    cpl_msg_error(modName, "Descriptor %s not found",
                pilTrnGetKeyword("OptDistOrdY"));
    return VM_FALSE;
  }

  for (i = 0; i <= optModY->orderX; i++) {
    for (j = 0; j <= optModY->orderY; j++) {
      sprintf(string, "%#.14E", optModY->coefs[i][j]);
      if (writeStringDescriptor(desc, pilTrnGetKeyword("OptDistY", i, j),
                                string, "") == VM_FALSE) {
        cpl_msg_error(modName, "Descriptor %s not found",
                    pilTrnGetKeyword("OptDistY", i, j));
        return VM_FALSE;
      }
    }
  }
  return VM_TRUE;
}


VimosBool readContaminationModel(VimosDescriptor *desc,
                   VimosDistModel2D **contModX, VimosDistModel2D **contModY)
{
  const char modName[] = "readContaminationModel";
  int        iValue;
  double     dValue;
  int        i, j;
  
  *contModX = NULL;
  *contModY = NULL;

  if (readIntDescriptor(desc, pilTrnGetKeyword("ZeroOrdX"), 
                              &iValue, NULL) == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found",  
                pilTrnGetKeyword("ZeroOrdX"));
    return VM_FALSE;
  } 
  
  *contModX = newDistModel2D(iValue, iValue);
  if (*contModX == NULL) {
    cpl_msg_error(modName, "The function newDistModel2D returned NULL");
    return(VM_FALSE);
  }
  for (i = 0; i <= iValue; i++) {
    for (j = 0; j <= iValue; j++) {
      if (readDoubleDescriptor(desc, pilTrnGetKeyword("ZeroX", i, j), 
                               &dValue, NULL) == VM_FALSE) {
        deleteDistModel2D(*contModX);
        *contModX = NULL;
        cpl_msg_error(modName, "Descriptor %s not found", 
                    pilTrnGetKeyword("ZeroX", i, j));
        return VM_FALSE;
      }
      else {
        (*contModX)->coefs[i][j] = dValue;
      }
    }
  }

  if (readIntDescriptor(desc, pilTrnGetKeyword("ZeroOrdY"), 
                        &iValue, NULL) == VM_FALSE) {
    deleteDistModel2D(*contModX);
    *contModX = NULL;
    cpl_msg_error(modName, "Descriptor %s not found",  
                pilTrnGetKeyword("ZeroOrdY"));
    return VM_FALSE;
  }
  
  *contModY = newDistModel2D(iValue, iValue);
  if (*contModY == NULL) {
    cpl_msg_error(modName, "The function newDistModel2D has returned NULL");
    return(VM_FALSE);
  }
  for (i = 0; i <= iValue; i++) {
    for (j = 0; j <= iValue; j++) {
      if (readDoubleDescriptor(desc, pilTrnGetKeyword("ZeroY", i, j),
                               &dValue, NULL) == VM_FALSE) {
        deleteDistModel2D(*contModX);
        deleteDistModel2D(*contModY);
        *contModX = NULL;
        *contModY = NULL;
        cpl_msg_error(modName, "Descriptor %s not found", 
                    pilTrnGetKeyword("ZeroY", i, j));
        return VM_FALSE;
      }
      else {
        (*contModY)->coefs[i][j] = dValue;
      }
    }
  }
  
  return VM_TRUE;
}

VimosBool writeContaminationModel(VimosDescriptor **desc, VimosDistModel2D 
                                  *contModX, VimosDistModel2D *contModY)
{
  const char modName[] = "writeContaminationModel";
  int        i, j;
  
  if (writeIntDescriptor(desc, pilTrnGetKeyword("ZeroOrdX"), 
                         contModX->orderX, "") == VM_FALSE) {
    cpl_msg_error(modName, "Descriptor %s not found", 
                pilTrnGetKeyword("ZeroOrdX"));
    return VM_FALSE;
  }
  
  for (i = 0; i <= contModX->orderX; i++) {
    for (j = 0; j <= contModX->orderX; j++) {
      if (writeDoubleDescriptor(desc, pilTrnGetKeyword("ZeroX", i, j), 
                                contModX->coefs[i][j], "") == VM_FALSE) {
        cpl_msg_error(modName, "Descriptor %s not found", 
                    pilTrnGetKeyword("ZeroX", i, j));
        return VM_FALSE;
      }
    }
  }
  
  if (writeIntDescriptor(desc, pilTrnGetKeyword("ZeroOrdY"), 
                         contModY->orderX, "") == VM_FALSE){
    cpl_msg_error(modName, "Descriptor %s not found", 
                pilTrnGetKeyword("ZeroOrdY"));
    return VM_FALSE;
  }

  for (i = 0; i <= contModY->orderX; i++) {
    for (j = 0; j <= contModY->orderX; j++) {
      if (writeDoubleDescriptor(desc, pilTrnGetKeyword("ZeroY", i, j), 
                                contModY->coefs[i][j], "") == VM_FALSE) {
        cpl_msg_error(modName, "Descriptor %s not found", 
                    pilTrnGetKeyword("ZeroY", i, j));
        return VM_FALSE;
      }
    }
  }

  return VM_TRUE;
}


VimosBool
readCurvatureModel(VimosDescriptor *desc, VimosDistModelFull **crvMod)
{
  const char modName[] = "readCurvatureModel";
  char      *descName;
  int        ordPol, ordX, ordY;
  double     dValue; 
  int        i, j, k;
  VimosBool  rdOK;
  
  *crvMod = NULL;

  descName = (char *) pilKeyTranslate("CurvatureOrd");
  if ((rdOK = readIntDescriptor(desc, descName, &ordPol, NULL)) 
      == VM_TRUE) {
    descName = (char *) pilKeyTranslate("CurvatureOrdX");
    if ((rdOK = readIntDescriptor(desc, descName, &ordX, NULL)) 
        == VM_TRUE) {
      descName = (char *) pilKeyTranslate("CurvatureOrdY");
      if ((rdOK = readIntDescriptor(desc, descName, &ordY, NULL)) 
          == VM_TRUE) {
        *crvMod = newDistModelFull(ordPol, ordX, ordY);
        if (*crvMod == NULL) {
          cpl_msg_error(modName, "Function newDistModelFull failure");
          return VM_FALSE;
        }
        for (i = 0; i <= ordPol; i++) {
          for (j = 0; j <= ordX; j++) {
            for (k = 0; k <= ordY; k++) {
              descName = (char *) pilKeyTranslate("Curvature", i, j, k);
              if ((readDoubleDescriptor(desc, descName, &dValue, NULL)) 
                                                              == VM_TRUE) {
                ((*crvMod)->coefs[i])->coefs[j][k] = dValue;
              }
              else {
                deleteDistModelFull(*crvMod);
                *crvMod = NULL;
                cpl_msg_error(modName, "Cannot read descriptor %s", descName);
                return VM_FALSE;
              }
            }
          }
        }
      }
    }
  }

  if (rdOK == VM_FALSE) 
    cpl_msg_error(modName, "Cannot read descriptor %s", descName);

  return rdOK;
}


VimosBool writeCurvatureModel(VimosDescriptor **desc, 
                              VimosDistModelFull *crvMod)
{
  const char modName[] = "writeCurvatureModel";
  char      *descName;
  int        i, j, k;
  VimosBool  wrOK;
  
  descName = (char *) pilKeyTranslate("CurvatureOrd");
  if ((wrOK = writeIntDescriptor(desc, descName, crvMod->orderPol, ""))
                                                               == VM_TRUE) {
    descName = (char *) pilKeyTranslate("CurvatureOrdX");
    if ((wrOK = writeIntDescriptor(desc, descName, crvMod->orderX, ""))
                                                               == VM_TRUE) {
      descName = (char *) pilKeyTranslate("CurvatureOrdY");
      if ((wrOK = writeIntDescriptor(desc, descName,  crvMod->orderY, ""))
                                                               == VM_TRUE) {
        for (i = 0; i <= crvMod->orderPol; i++) {
          for (j = 0; j <= crvMod->orderX; j++) {
            for (k = 0; k <= crvMod->orderY; k++) {
              descName = (char *) pilKeyTranslate("Curvature", i, j, k);
              if ((wrOK = writeDoubleDescriptor(desc, descName,
                          crvMod->coefs[i]->coefs[j][k],"")) == VM_FALSE) {
                cpl_msg_error(modName, "Cannot write descriptor %s", descName);
                return VM_FALSE;
              }
            }
          }
        }
      }
    }
  }
  
  if (wrOK == VM_FALSE) 
    cpl_msg_error(modName, "Cannot write descriptor %s", descName);

  return wrOK;
}


/*
 * The following ludicrous function is a copy of the writeCurvatureModel()
 * (above) that has been made necessary by a bug that is still present
 * in the OSLX: when writing any distorsion model to an image header,
 * its coefficients must be character strings. When one day this problem
 * will be solved, this function should disappear, and any call to it
 * should be replaced by a call to writeCurvatureModel().  (C.Izzo)
 */

VimosBool writeCurvatureModelString(VimosDescriptor **desc,
                                    VimosDistModelFull *crvMod)
{
  const char modName[] = "writeCurvatureModel";
  char      *descName;
  int        i, j, k;
  char       string[80];
  VimosBool  wrOK;

  descName = (char *) pilKeyTranslate("CurvatureOrd");
  if ((wrOK = writeIntDescriptor(desc, descName, crvMod->orderPol, ""))
                                                               == VM_TRUE) {
    descName = (char *) pilKeyTranslate("CurvatureOrdX");
    if ((wrOK = writeIntDescriptor(desc, descName, crvMod->orderX, ""))
                                                               == VM_TRUE) {
      descName = (char *) pilKeyTranslate("CurvatureOrdY");
      if ((wrOK = writeIntDescriptor(desc, descName,  crvMod->orderY, ""))
                                                               == VM_TRUE) {
        for (i = 0; i <= crvMod->orderPol; i++) {
          for (j = 0; j <= crvMod->orderX; j++) {
            for (k = 0; k <= crvMod->orderY; k++) {
              sprintf(string, "%#.14E", crvMod->coefs[i]->coefs[j][k]);
              descName = (char *) pilKeyTranslate("Curvature", i, j, k);
              if ((wrOK = writeStringDescriptor(desc, descName, string, "")) 
                                                                == VM_FALSE) {
                cpl_msg_error(modName, "Cannot write descriptor %s", descName);
                return VM_FALSE;
              }
            }
          }
        }
      }
    }
  }

  if (wrOK == VM_FALSE)
    cpl_msg_error(modName, "Cannot write descriptor %s", descName);

  return wrOK;
}


VimosBool readInvDispMatrix(VimosDescriptor *desc, 
                            VimosDistModelFull **invDispMat)
{
  const char  modName[] = "readInvDispMatrix";
  char       *descName;
  int         ordPol, ordX, ordY;
  double      dValue;
  int         i, j, k;
  VimosBool   rdOK;
  
  *invDispMat = NULL;

  descName = (char *) pilKeyTranslate("DispersionOrd");
  if ((rdOK = readIntDescriptor(desc, descName, &ordPol, NULL))
      == VM_TRUE) {
    descName = (char *) pilKeyTranslate("DispersionOrdX");
    if ((rdOK = readIntDescriptor(desc, descName, &ordX, NULL))
        == VM_TRUE) {
      descName = (char *) pilKeyTranslate("DispersionOrdY");
      if ((rdOK = readIntDescriptor(desc, descName, &ordY, NULL))
          == VM_TRUE) {
        *invDispMat = newDistModelFull(ordPol, ordX, ordY);
        if (*invDispMat == NULL) {
          cpl_msg_error(modName, "Function newDistModelFull failure");
          return VM_FALSE;
        }
        for (i = 0; i <= ordPol; i++) {
          for (j = 0; j <= ordX; j++) {
            for (k = 0; k <= ordY; k++) {
              descName = (char *) pilKeyTranslate("Dispersion", i, j, k);
              if ((rdOK = readDoubleDescriptor(desc, descName, &dValue, 
                                                    NULL)) == VM_TRUE) {
                ((*invDispMat)->coefs[i])->coefs[j][k] = dValue;
              }
              else {
                deleteDistModelFull(*invDispMat);
                *invDispMat = NULL;
                cpl_msg_error(modName, "Cannot read descriptor %s", descName);
                return VM_FALSE;
              }
            }
          }
        }
      }
    }
  }

  if (rdOK == VM_FALSE) 
    cpl_msg_error(modName, "Cannot read descriptor %s", descName);

  return rdOK;

}

VimosBool writeInvDispMatrix(VimosDescriptor **desc, VimosDistModelFull *idsMat)
{
  const char modName[] = "writeInvDispMatrix";
  char      *descName;
  int        i, j, k;
  VimosBool  wrOK;

  descName = (char *) pilKeyTranslate("DispersionOrd");
  if ((wrOK = writeIntDescriptor(desc, descName, idsMat->orderPol, ""))
      == VM_TRUE) {
    descName = (char *) pilKeyTranslate("DispersionOrdX");
    if ((wrOK = writeIntDescriptor(desc, descName, idsMat->orderX, ""))
        == VM_TRUE) {
      descName = (char *) pilKeyTranslate("DispersionOrdY");
      if ((wrOK = writeIntDescriptor(desc, descName,  idsMat->orderY, ""))
          == VM_TRUE) {
        for (i = 0; i <= idsMat->orderPol; i++) {
          for (j = 0; j <= idsMat->orderX; j++) {
            for (k = 0; k <= idsMat->orderY; k++) {
              descName = (char *) pilKeyTranslate("Dispersion", i, j, k);
              if ((wrOK = writeDoubleDescriptor(desc, descName,
                           idsMat->coefs[i]->coefs[j][k],"")) == VM_FALSE) {
                cpl_msg_error(modName, "Cannot write descriptor %s", descName);
                return VM_FALSE;
              }
            }
          }
        }
      }
    }
  }

  if (wrOK == VM_FALSE) 
    cpl_msg_error(modName, "Cannot write descriptor %s", descName);

  return wrOK;
}

VimosBool writeInvDispMatrixString(VimosDescriptor **desc, 
                                   VimosDistModelFull *idsMat)
{
  const char modName[] = "writeInvDispMatrix";
  char      *descName;
  int        i, j, k;
  char       string[80];
  VimosBool  wrOK;

  descName = (char *) pilKeyTranslate("DispersionOrd");
  if ((wrOK = writeIntDescriptor(desc, descName, idsMat->orderPol, ""))
                                                          == VM_TRUE) {
    descName = (char *) pilKeyTranslate("DispersionOrdX");
    if ((wrOK = writeIntDescriptor(desc, descName, idsMat->orderX, ""))
                                                          == VM_TRUE) {
      descName = (char *) pilKeyTranslate("DispersionOrdY");
      if ((wrOK = writeIntDescriptor(desc, descName,  idsMat->orderY, ""))
                                                          == VM_TRUE) {
        for (i = 0; i <= idsMat->orderPol; i++) {
          for (j = 0; j <= idsMat->orderX; j++) {
            for (k = 0; k <= idsMat->orderY; k++) {
              descName = (char *) pilKeyTranslate("Dispersion", i, j, k);
              sprintf(string, "%#.14E", idsMat->coefs[i]->coefs[j][k]);
              if ((wrOK = writeStringDescriptor(desc, descName,
                           string, "")) == VM_FALSE) {
                cpl_msg_error(modName, "Cannot write descriptor %s", descName);
                return VM_FALSE;
              }
            }
          }
        }
      }
    }
  }

  if (wrOK == VM_FALSE)
    cpl_msg_error(modName, "Cannot write descriptor %s", descName);

  return wrOK;
}

#define rad2deg(x) (57.2957795130823208768*(x))   /* radians -> degrees */
#define deg2rad(x) ( 0.0174532925199432958*(x))   /* degrees -> radians */

VimosGnomonic *newGnomonic(const double alpha0, const double delta0)
{
  const char modName[] = "newGnomonic";
  VimosGnomonic *tmpGnome;
  
  /* allocate space */
  tmpGnome = (VimosGnomonic *) cpl_malloc(sizeof(VimosGnomonic));
  if (tmpGnome == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }
  /* fill up elements. Convert from degrees to radians */
  tmpGnome->alpha0 = deg2rad(alpha0);
  tmpGnome->delta0 = deg2rad(delta0);
  tmpGnome->sina0  = sin(deg2rad(alpha0));
  tmpGnome->cosa0  = cos(deg2rad(alpha0));
  tmpGnome->sind0  = sin(deg2rad(delta0));
  tmpGnome->cosd0  = cos(deg2rad(delta0));
  
  return(tmpGnome);
}



void deleteGnomonic(VimosGnomonic *gnome)
{
  if (gnome != NULL) {
    cpl_free(gnome);
  }
}



VimosBool lm2RaDec(VimosGnomonic *gnome, double l, double m, double *ra, 
                   double *dec)
{
  const char modName[] = "lm2RaDec";
  double s, t, da;
  
  if (gnome == NULL) {
    *ra  = 0.0;
    *dec = 0.0;
    cpl_msg_error(modName, "NULL input pointer");
    return(VM_FALSE);
  }
  
  s = (m * gnome->cosd0 + gnome->sind0);
  t = (gnome->cosd0 - m * gnome->sind0);
  da = atan(l/t);
  *ra = rad2deg(gnome->alpha0 +  da);
  *dec = rad2deg(atan(cos(da) * s/t));

  return(VM_TRUE);  
}

VimosBool raDec2lm(VimosGnomonic *gnome, double ra, double dec, double *l, 
                   double *m)
{
  const char modName[] = "raDec2lm";
  double da, cosa, sina, cosd, sind, t;
  
  if (gnome == NULL) {
    *l = 0.0;
    *m = 0.0;
    cpl_msg_error(modName, "NULL input pointer");
    return(VM_FALSE);
  }


  da = deg2rad(ra) - gnome->alpha0;
  cosa = cos(da);
  sina = sin(da);
  cosd = cos(deg2rad(dec));
  sind = sin(deg2rad(ra));
  t = (sind*gnome->sind0  +  cosd*gnome->cosd0*cosa);
  *l = cosd * sina / t;
  *m = (sind*gnome->cosd0 - cosd*gnome->sind0*cosa) / t;

  return(VM_TRUE);
}

#undef rad2deg
#undef deg2rad

VimosBool fitDistModel2D(VimosPixel *surface, int numPoints, int polyDeg,  
                         double offsetX, double offsetY, 
                         VimosDistModel2D **model, double *rms)
{
  const char  modName[] = "fitDistModel2D";
  char       *controlString;
  VimosPixel *tmpPix;
  double     *tmpCoefs;
  int         i;
  int         degx, degy;
  /* we need a dummy integer argument for fitSurfacePolynomial: */
  int bla;
  
  /* create pixel list that we will use here */
  tmpPix = newPixel(numPoints);
  if (tmpPix == NULL) {
    cpl_msg_error(modName, "Function newPixel failure");
    return(VM_FALSE);
  } 
  
  /* copy data from input pixel list and apply offset in X and Y */
  for (i = 0; i < numPoints; i++) {
    tmpPix[i].x = surface[i].x - offsetX;
    tmpPix[i].y = surface[i].y - offsetY;
    tmpPix[i].i = surface[i].i;
  }
  
  /* fit surface of order polyDeg to surface */
  controlString =
  createVimosCtrlStr(polyDeg, polyDeg);
  tmpCoefs = fitSurfacePolynomial(tmpPix, numPoints, controlString, 
                                  polyDeg+polyDeg, &bla, rms);
  if (tmpCoefs == NULL) {
    cpl_msg_error(modName, "Function fitSurfacePolynomial failure");
    return(VM_FALSE);
  }
  /* create output model */
  *model = newDistModel2D(polyDeg, polyDeg);
  if (*model == NULL) {
    cpl_msg_error(modName, "Function newDistModel2D failure");
    return(VM_FALSE);
  }
  /* copy offsets */
  (*model)->offsetX = offsetX;
  (*model)->offsetY = offsetY;
  
  /* copy coefficients to model */
  i = 0 ;
  for (degx=0; degx <= polyDeg; degx++) {
    for (degy=0; degy <= polyDeg; degy++) {
      (*model)->coefs[degx][degy] = tmpCoefs[i];
      i++ ;
    }
  }
  
  /* free temp array */
  cpl_free(tmpCoefs);

  return(VM_TRUE);

}

