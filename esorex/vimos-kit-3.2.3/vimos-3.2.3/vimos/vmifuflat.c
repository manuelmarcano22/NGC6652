/* $Id: vmifuflat.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <piltranslator.h>

#include "vmtable.h"
#include "vmextractiontable.h"
#include "vmifuflat.h"


#define EXTRA              20.


/**
 * @brief vmifuflat IFU Flat Fields
 *
 * The module provides the IFU flat field related functions
 */

/**@{*/

VimosBool
determineExposedIfuArea(VimosTable *adf,VimosExtractionTable *extTable,
                        int *nSlits, VimosDpoint **minY, VimosDpoint **maxY,
                        VimosDpoint **zeroY)
{

  const char modName[] = "determineExposedIfuArea";
  int        nPixBelow;
  int        nPixAbove;
  float      ySlit;
  float      yZero;
  char       comment[80];

  VimosExtractionSlit *exSlit;

  /* validate input */
  if (extTable == NULL) {
    pilMsgError(modName, "NULL input table");
    return (VM_FALSE);
  }

  exSlit = extTable->slits;

  if (!determineExposedIfuSlit(adf, exSlit, &ySlit, &yZero)) {
    pilMsgError(modName, "Function determineExposedIfuSlit returned an error");
    return(VM_FALSE);
  }
  *nSlits = 1;
  *minY = newDpoint(*nSlits);
  if (*minY == NULL) {
    pilMsgError(modName, "The function newDpoint has returned NULL");
    return(VM_FALSE);
  }
  *maxY = newDpoint(*nSlits);
  if (*maxY == NULL) {
    pilMsgError(modName, "The function newDpoint has returned NULL");
    return(VM_FALSE);
  }
  *zeroY = newDpoint(*nSlits);
  if (*zeroY == NULL) {
    pilMsgError(modName, "The function newDpoint has returned NULL");
    return(VM_FALSE);
  }

  if (!readIntDescriptor(extTable->descs, "ESO PRO SPECT LLEN LO", &nPixBelow,
        comment)) {
    pilMsgError(modName, "Function readIntDescriptor has returned an error");
    return(VM_FALSE);
  }
  if (!readIntDescriptor(extTable->descs, "ESO PRO SPECT LLEN HI", &nPixAbove,
        comment)) {
    pilMsgError(modName, "Function readIntDescriptor has returned an error");
    return(VM_FALSE);
  }

  (*minY)->x = 0.;
  (*minY)->y = (double) (ySlit - nPixBelow - EXTRA);
  (*maxY)->y = (double) (ySlit + nPixAbove + EXTRA);
  (*zeroY)->y = (double) (ySlit + yZero);

  return VM_TRUE;
}


VimosBool
determineExposedIfuSlit(VimosTable *adf, VimosExtractionSlit *exSlit,
                        float *ySlit, float *yZero)
{
  const char  modName[] = "determineExposedIfuSlit";
  int         i, j;
  int         quadNum;
  float       yMin;
  float       yMax;
  int         ifuSlit;
  float       ifuSlitY;
  float       ifuZeroY;
  char        comment[80];

  /* find out the minimum and maximum y coordinate in the illuminated area */

  if (readIntDescriptor(adf->descs, pilTrnGetKeyword("Quadrant"), 
                        &quadNum, comment) == VM_FALSE) {
    pilMsgError(modName, "Keyword %s not found", pilTrnGetKeyword("Quadrant"));
    return VM_FALSE;
  }

  if (readFloatDescriptor(adf->descs, pilTrnGetKeyword("MshuPosH", quadNum), 
                          &yMax, comment) == VM_FALSE) {
    pilMsgError(modName, "Keyword %s not found", 
                pilTrnGetKeyword("MshuPosH", quadNum));
    return VM_FALSE;
  }

  if (readFloatDescriptor(adf->descs, pilTrnGetKeyword("MshuPosL", quadNum), 
                          &yMin, comment) == VM_FALSE) {
    pilMsgError(modName, "Keyword %s not found", 
                pilTrnGetKeyword("MshuPosL", quadNum));
    return VM_FALSE;
  }

  for (i = 0; i <= 3; i++) {
    ifuSlit = exSlit->IFUslitNo;
    ifuSlitY = 0.;
    ifuZeroY = 0.;
    j = 0;

    while (exSlit->IFUslitNo == ifuSlit) {
      ifuSlitY += exSlit->ccdY->data[0];
      ifuZeroY += exSlit->zeroY->data[0];
      j++;
      exSlit = exSlit->next;
    }

    ifuSlitY /= (float) j;
    ifuZeroY /= (float) j;

    if (ifuSlitY > yMin && ifuSlitY < yMax) {
      *ySlit = ifuSlitY;
      *yZero = ifuZeroY;
      return (VM_TRUE);
    }
  }
  return (VM_FALSE);
}
/**@}*/

