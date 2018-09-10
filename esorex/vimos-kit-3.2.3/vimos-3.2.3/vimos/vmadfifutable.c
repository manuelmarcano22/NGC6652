/* $Id: vmadfifutable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
 
#include <pilmessages.h>

#include "vmtable.h"
#include "vmifutable.h"
#include "vmadf.h"
#include "vmadfifutable.h"


VimosAdfSlitHolder *extractSlitsFromIFU(VimosTable *adf, 
			      VimosIfuTable *ifuTab, VimosIfuMode ifuMode)
{
  int                i, quadVal;
  float              radiusVal;
  double             dValue;
  char               comment[80];
  VimosIfuQuad       *ifuQuad;
  VimosIfuSlit       *ifuSlit;
  VimosIfuFiber      *ifuFiber;
  VimosAdfSlitType   slitType;
  VimosAdfCircSlit   *tCircSlit;
  VimosAdfSlitHolder *slitHolder;
  VimosAdfSlitHolder *tSlitHolder;
  VimosAdfSlitHolder *linkedSlits=NULL;


  /* validate input */
  if (ifuTab == NULL) {
    pilMsgError("extractSlitsFromIFU","NULL input table");
    return(NULL);
  }

  if (ifuMode == VM_IFU_MODE_UNDEF) {
    pilMsgError("extractSlitsFromIFU","Undefined input mode");
    return(NULL);
  }

  if (!readIntDescriptor(adf->descs, "ESO OCS CON QUAD", &quadVal, comment)) {
    pilMsgError("extractSlitsFromIFU","The function readIntDescriptor has returned an error");
    return(NULL);
  }

  if (!readDoubleDescriptor(ifuTab->descs, "ESO PRO FIBER RADIUS", &dValue,
			    comment)) {
    pilMsgError("extractSlitsFromIFU","The function readDoubleDescriptor has returned an error");
    return(NULL);
  }

  radiusVal = (float) dValue;

  /* get quadrants in the table */
  ifuQuad = ifuTab->quads;
  i = 0;
  slitHolder = NULL;
  slitType = VM_ADF_CIRC_SLIT;
 
  while (ifuQuad)
    {
      if (quadVal != ifuQuad->quadNo) {
	ifuQuad = ifuQuad->next;
	continue;
      }	

      /* get slits in the quadrant */
     ifuSlit = ifuQuad->ifuSlits;
 
     while (ifuSlit)
       {
	 /* test to see if we are in High Resolution mode. In this
	    case only the slit nr. 2 is illuminated */
	 if ((ifuMode == VM_IFU_MODE_SMALL) && (ifuSlit->ifuSlitNo != 2)) {
	   ifuSlit = ifuSlit->next;
	   continue;
	 }

         /* get fibers in the slit */
        ifuFiber = ifuSlit->fibers;

	/* cycle on all fibers */
        while (ifuFiber)
          {
	    /* here goes the conversion... */
	    tCircSlit =  newAdfCircSlit();
	    if (tCircSlit == NULL) {
	      pilMsgError("extractSlitsFromIFU","The function newAdfCircSlit has returned NULL");
	      return(NULL);
	    }
	    tCircSlit->slitNo = i+1;
	    tCircSlit->x = ifuFiber->fiberx;
	    tCircSlit->y = ifuFiber->fibery;
	    tCircSlit->radius = radiusVal;
	    tCircSlit->IFUslitNo = ifuSlit->ifuSlitNo;
	    tCircSlit->IFUfibNo = ifuFiber->fibNo;
	    tCircSlit->IFUfibTrans = ifuFiber->fiberTrans;

	    if (i == 0) {
	      slitHolder = newAdfSlitHolder();
	      if (slitHolder == NULL) {
		pilMsgError("extractSlitsFromIFU","The function newAdfSlitHolder has returned NULL");
		return(NULL);
	      }

	      linkedSlits = slitHolder;
	      linkedSlits->prev = NULL;
	      linkedSlits->next = NULL;
	    }
	    else {
	      tSlitHolder = newAdfSlitHolder();
	      if (tSlitHolder == NULL) {
		pilMsgError("extractSlitsFromIFU","The function newAdfSlitHolder has returned NULL");
		return(NULL);
	      }

	      tSlitHolder->prev = linkedSlits;
	      linkedSlits->next = tSlitHolder;
	      linkedSlits = tSlitHolder;
	    }
	    linkedSlits->slit = tCircSlit;
	    linkedSlits->slitNo = i+1;
	    linkedSlits->slitType = slitType;

	    i++;
 
	    /* get next fiber */
	    ifuFiber = ifuFiber->next;
 
          }
 
        /* get next slit */
        ifuSlit = ifuSlit->next;
 
       }
 
     /* get next quadrant */
     ifuQuad = ifuQuad->next;
 
    }

  return(slitHolder);

}
