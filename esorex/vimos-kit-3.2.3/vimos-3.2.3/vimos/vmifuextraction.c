/* $Id: vmifuextraction.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <pilmessages.h>

#include "vmifutable.h"
#include "vmextractiontable.h"
#include "vmwindowtable.h"
#include "vmifuextraction.h"


VimosWindowTable *
VmIfuWinTab(VimosIfuTable *ifuTable, VimosExtractionTable *extTable,
            int quadNum)
{
  int i = 0;
  int j = 0;
  int kk;
  int objPos;

  char modName[] = "VmIfuWinTab";

  VimosWindowTable *ifuWinTable;

  VimosIfuQuad *theQuads;
  VimosIfuSlit *theIfuSlits;
  VimosIfuFiber *theFibers;

  VimosWindowSlit *winSlit;
  VimosWindowSlit *lastWinSlit;
  VimosWindowObject *winObject;

  VimosExtractionSlit *theExtSlits;

  pilMsgInfo (modName, "Computing IFU Window Table");

  ifuWinTable = newWindowTable();
  /* copy descriptors of the IFU Extraction Table to the IFU Window Table */
  copyExtTab2WinTab(extTable, ifuWinTable);

  lastWinSlit = NULL;

  theQuads = ifuTable->quads;

  while (theQuads)
    {
     /* take only the quadrant correspondent to the extraction table */
     if (theQuads->quadNo == quadNum)
       {
        /* check against extSlit */
        theExtSlits = extTable->slits;

	/* loop on extSlits for this quadrant */
        while (theExtSlits)
          {
           i = 0;

	   /* for each Ext slit take the IFU slits on this quadrant */
           theIfuSlits = theQuads->ifuSlits;

           while (theIfuSlits)
             {
	      /* for each slit, take fibers */
	      theFibers = theIfuSlits->fibers;

              while (theFibers)
                {
	         i++;

                 if ((theFibers->fibNo == theExtSlits->IFUfibNo) &&
                     (theIfuSlits->ifuSlitNo == theExtSlits->IFUslitNo))
                    {
	             j++;

                     winSlit = newWindowSlit();

		     /* SETTING THE SLIT AS SHORT TYPE BY DEFAULT */
		     winSlit->specLong = VM_FALSE;

                     winSlit->specStart = theExtSlits->y->data[0];
                     winSlit->specEnd = theExtSlits->y->data[0] +
                                        theExtSlits->numRows - 1;

		     /* give here the slit number for IFU slits */
                     winSlit->IFUslitNo = theExtSlits->IFUslitNo;

                     /* give here the fiber number for IFU fibers */
                     winSlit->IFUfibNo =theExtSlits->IFUfibNo;


		     /* set here the fiber transmission */
                     winSlit->IFUfibTrans =theExtSlits->IFUfibTrans;


                     winSlit->slitNo = theExtSlits->slitNo;

		     /* Only 1 obj spectrum for each window slit in IFU case */
                     winSlit->numObj = 1;

                     winObject = newWindowObject();
		     
		     /* Set this to 0, it is the first row of obj spectrum */
		     /* RELATIVE TO SPEC_START: in IFU the obj spectrum spans*/
		     /* the whole fiber (spatially) */
                     winObject->objStart = 0;

		     /* same as for objStart: RELATIVE TO SPEC_START */
                     winObject->objEnd = winSlit->specEnd - 
                                         winSlit->specStart;

		     /* set the object profile */
		     winObject->objProfile=newFloatArray(winObject->objEnd - 
							winObject->objStart+1);

		     /* WARNING!!! for the moment set to 0, to be changed */
		     for (kk=winObject->objStart; kk<=winObject->objEnd;
			  kk++) 
		       {
			winObject->objProfile->data[kk] = 0.0;
			/*
			if ((theFibers->fibNo == 1) &&
			    (theIfuSlits->ifuSlitNo ==1))
			pilMsgInfo(modName," %f ",
				   winObject->objProfile->data[kk]);
			*/
		       }
		     /* Only 1 obj spectrum for each window slit in IFU case */
                     winObject->objNo = 1;

		     /* assign the object position as the average (central) */

		     /*IS IT OK? (in findObiectsInProfile it is said "should */
		     /* make more sofisticated later")*/
	             winObject->objPos = ( winObject->objStart +
                                          winObject->objEnd ) / 2.;

                     objPos = winObject->objPos;
                     winObject->objX = theExtSlits->maskX->data[objPos];
                     winObject->objY =theExtSlits->maskY->data[objPos];
                     winObject->posDef = VM_FALSE;

		     /* WHAT ABOUT THESE PARAMETERS? */

                     winObject->objRA =0.;
                     winObject->objDec =0.;

		     /* these are commented in the windowTable declaration...
                     winObject->objWidth =0.;
                     winObject->objPeak =0.;
                     winObject->skyLevel =0.;
		     */

                     winObject->parDef = winObject->posDef = VM_FALSE;

		     /* we have 1 object per fiber, don't link the obj list */ 
		     /* ok? */
                     winSlit->objs = winObject;

                     /* link Window Slit to list */
                     if (lastWinSlit == NULL)
                       {
			 /* first slit we do */
                        ifuWinTable->slits = winSlit;
                       }
                     else
                       {
			 /* not first, add to tail */
                        lastWinSlit->next = winSlit;
                        winSlit->prev = lastWinSlit;
                       }
                     /* store tail of linked list */
                     lastWinSlit = winSlit;

                    }

                 theFibers = theFibers->next;
                } /* end loop on ifu slit fibers */

              theIfuSlits = theIfuSlits->next;
             } /* end loop in ifu slits */

           theExtSlits = theExtSlits->next;
          } /* end loop on extSlits */

       } /* end "if the quadrant is the right one" */

     theQuads = theQuads->next;
    } /* end loop on quadrants to find the right one */

  return ifuWinTable;

}
