/* $Id: vmifutable.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
#include <string.h>
 
#include <fitsio.h>

#include <pilmemory.h>
#include <pilmessages.h>

#include "vmtable.h"
#include "vmifutable.h"
#include "cpl.h"


/*
 * Returns a pointer to a new IFU Fiber
 */

VimosIfuFiber *newIfuFiber(void)
{
 VimosIfuFiber *newFiber;
  /* allocate memory for IFU Fiber */
  newFiber = (VimosIfuFiber *) cpl_malloc(sizeof(VimosIfuFiber));

  /* check if space was allocated */
  if (newFiber == NULL) {
    pilMsgError("newIfuFiber","Allocation Error");
    return(NULL);
  }

  /* fill up fields with default values */
  newFiber->fibNo       =0;
  newFiber->fiberL      =0;
  newFiber->fiberM      =0;
  newFiber->fiberx      =0.0;
  newFiber->fibery      =0.0;
  newFiber->fiberTrans  =1.0;
  newFiber->fiberPwidth =1.0;
  newFiber->sigmaY      =0.0;
  newFiber->sigmaYGroup =0;

  newFiber->prev = newFiber->next = NULL;
  /* return to caller */
  return(newFiber);

}

/*
Deletes all IFU Fibers contained in the list
*/

void deleteIfuFiber(VimosIfuFiber *fiber)
{
 VimosIfuFiber  *tmpFiber;
 VimosIfuFiber  *nextFiber;

  /* store start of the list */
  tmpFiber = fiber;

  /* traverse  list */
  while (tmpFiber) {
    /* get address of next fiber */
    nextFiber = tmpFiber->next;
    /* free current object */
    cpl_free(tmpFiber);
    /* next one to process */
    tmpFiber = nextFiber;
  }
}

/*
  Returns a pointer to a new IFU Slit
*/

VimosIfuSlit *newIfuSlit()
{
  VimosIfuSlit *newSlit;
  /* allocate memory for IfuSlit */
  newSlit = (VimosIfuSlit *) cpl_malloc(sizeof(VimosIfuSlit));
  
  /* check if space was allocated */
  if (newSlit == NULL) {
    pilMsgError("newIfuSlit","Allocation Error");
    return(NULL);
  }

  /* fill up fields with default values */
  newSlit->ifuSlitNo = 0;
  newSlit->fibers = NULL;
  newSlit->prev = newSlit->next = NULL;

  /* return to caller */
  return(newSlit);
  
}

/*
Deletes all IFU Slits contained in the list
*/

void deleteIfuSlit(VimosIfuSlit *ifuSlit)
{
  VimosIfuSlit *tmpIfuSlit;
  VimosIfuSlit *nextIfuSlit;

  /* store start of the list */
  tmpIfuSlit = ifuSlit;

  /* traverse  list */
  while (tmpIfuSlit) {
    /* delete the fibers in this slit */
    deleteIfuFiber(tmpIfuSlit->fibers);
    /* get address of next IFU slit */
    nextIfuSlit = tmpIfuSlit->next;
    /* free current IFU slit */
    cpl_free(tmpIfuSlit);
    /* next one to process */
    tmpIfuSlit = nextIfuSlit;
  }
}

/*
  Returns a pointer to a new IFU Quadrant
*/

VimosIfuQuad *newIfuQuad()
{
  VimosIfuQuad *ifuQuad;
    /* allocate memory for ifuQuad */
  ifuQuad = (VimosIfuQuad *) cpl_malloc(sizeof(VimosIfuQuad));
  
  /* check if space was allocated */
  if (ifuQuad == NULL) {
    pilMsgError("newIfuQuad","Allocation Error");
    return(NULL);
  }

  /* fill up fields with default values */
  ifuQuad->quadNo = 0;
  ifuQuad->ifuSlits = NULL;

  /* set to 4 the number of slits in each quadrant */
  ifuQuad->numIfuSlits = 4;

  ifuQuad->prev = ifuQuad->next = NULL;

  /* return to caller */
  return(ifuQuad);
  
}


/*
Deletes all IFU quadrants contained in the list
*/

void deleteIfuQuad(VimosIfuQuad *ifuQuad)
{
  VimosIfuQuad *tmpIfuQuad;
  VimosIfuQuad *nextIfuQuad;

  /* store start of the list */
  tmpIfuQuad = ifuQuad;

  /* traverse  list */
  while (tmpIfuQuad) {
    /* delete the slits in this quadrant. This should delete also all */
    /* fibers in each slit */
    deleteIfuSlit(tmpIfuQuad->ifuSlits);
    /* get address of next IFU quadrant */
    nextIfuQuad = tmpIfuQuad->next;
    /* free current IFU quadrant */
    cpl_free(tmpIfuQuad);
    /* next one to process */
    tmpIfuQuad = nextIfuQuad;
  }
}


/*
  Allocate a new IFU Table. 
*/
VimosIfuTable *newIfuTable()
{
  VimosIfuTable *ifuTab;

  /* allocate new VimosIfuTable */
  ifuTab = (VimosIfuTable *) cpl_malloc(sizeof(VimosIfuTable));
 
  /* check if space was allocated */
  if (ifuTab == NULL) {
    pilMsgError("newIfuTable","Allocation Error");
    return(NULL);
  }

  /* copy "IFU" into name of table and initialize fields*/
  strcpy(ifuTab->name, VM_IFU);
  ifuTab->descs = newStringDescriptor("ESO PRO TABLE", "IFU", "");
  if (ifuTab->descs == NULL) {
    /* cleanup */
    cpl_free(ifuTab);
    pilMsgError("newIfuTable","The function newStringDescriptor has returned NULL");
    return(NULL);
  }

  ifuTab->quads   = NULL;

  /* set the number of quadrants to 4 */
  ifuTab->numIfuQuads = 4;

  /* set the total number of fibers in IFU */
  ifuTab->numIfuFibs = 6400;

  ifuTab->fptr = NULL;

  /* return address of new IFU Table */
  return(ifuTab);
  
}


/*
 
*/
void deleteIfuTable(VimosIfuTable *ifuTable)
{
  VimosDescriptor *tmpDesc;
  VimosDescriptor *nxtDesc;
  
  if (ifuTable == NULL) {
    return;
  }

  deleteIfuQuad(ifuTable->quads);

  tmpDesc = ifuTable->descs;  
  while (tmpDesc) {
    nxtDesc = tmpDesc->next;
    deleteDescriptor(tmpDesc);
    tmpDesc = nxtDesc;
  }
 
}


VimosIfuSlit *computeIfuSlit(int l_start, int m_start, int l_step,
                             int m_step, int module_step_m, float x_start,
                             float x_step, float y, float x_module_step)
{
 int k, theFibNum;
 int i,j,l;
 int m_start_ini;
 int l_dir = -1;
 int kk;
 VimosIfuFiber *ifuFiber;
 VimosIfuFiber *lastIfuFiber;
 VimosIfuSlit *theIfuSlit;

/* k says this is the first fiber in the slit */
 k=0;

 m_start_ini = m_start;
 lastIfuFiber = NULL;

 theIfuSlit = newIfuSlit();
 if (theIfuSlit == NULL) {
   pilMsgError("computeIfuSlit","The function newIfuSlit has returned NULL");
   return(NULL);
 }

 /* initialize the fiber number in this slit */
 theFibNum = 1;

 /* build the five modules for each slit */
 for (i=0; i<5; i++)
    {
     /* for each module, 4 lines are read */
     for (j=0; j<4; j++)
        {
        /* kk=0 says this is the first fiber in the module */
         kk=0;

         /* for each line 20 fibers are read */
         for (l=0; l<20; l++)
            {
             ifuFiber = newIfuFiber();
             if (ifuFiber == NULL) {
               pilMsgError("computeIfuSlit","The function newIfuFiber has returned NULL");
             return(NULL);
 }

	     /* set the fiber number in the slit */
             ifuFiber->fibNo = theFibNum;

             /* store first fiber in the module */
             if (kk == 0)
               {
                ifuFiber->fiberL = l_start;
                ifuFiber->fiberM = m_start;
               }
             else
               {
                ifuFiber->fiberL = lastIfuFiber->fiberL + l_step;
                ifuFiber->fiberM = m_start;
               }
             if (k == 0)
               {
                ifuFiber->fiberx = x_start;
               }
             else if (k !=0)
               {
                x_start = x_start + x_step;
               }
             ifuFiber->fiberx = x_start;
             ifuFiber->fibery = y;

             /* printf("FIBER X %6.1f",ifuFiber->fiberx); */
             /* printf("FIBER Y %6.1f\n",ifuFiber->fibery); */

             /* link fiber to list */
             if (lastIfuFiber == NULL) 
               {
                /* first fiber we do */
                theIfuSlit->fibers = ifuFiber;
               }
             else
               {
                /* not first, add to tail */
                lastIfuFiber->next = ifuFiber;
                ifuFiber->prev = lastIfuFiber;
               }
             /* store tail of linked list */
             lastIfuFiber = ifuFiber;

             /* linking done */

             kk = kk+1;
             k = k+1;

             theFibNum++;

            } /* end reading 20 fibers */

         m_start=m_start+m_step;
         l_step=l_step*l_dir;
         l_start = ifuFiber->fiberL;

        } /* end reading 4 lines */
     /* printf("END WRITING MODULE\n"); */
     /* mechanical:  1mm between adjacent 80 fibres modules */
     x_start = x_start + x_module_step;
     m_start=m_start-m_step;
     m_start=m_start_ini+(module_step_m*(i+1));
    } /* end reading 5 modules for each slit */
 /* puts("end computeIfuSlit"); */

 return(theIfuSlit);

}



VimosBool writeTable(VimosIfuTable *anIfuTable)
{

  VimosIfuQuad *theQuads;
  VimosIfuSlit *theSlits;
  VimosIfuFiber *theFib;

  FILE *ofp;
  const char fileName[20] = "ifuTable.dat";
  ofp = fopen(fileName, "w");
  if (ofp == NULL)
    {
      pilMsgError("writeTable","The file cannot be opened");
      return(VM_FALSE);
    }

  theQuads = anIfuTable->quads;
  while (theQuads)
    {
      /* printf("writing quadrant %d: \n",theQuads->quadNo); */
     theSlits = theQuads->ifuSlits;
     while (theSlits)
       {   
	 /* printf("writing slit %d: \n",theSlits->ifuSlitNo); */
        theFib = theSlits->fibers;
        while(theFib)
          {
           fprintf(ofp,
		   "%2d %2d %5d %3d %3d %6.1f %6.1f %8.5f %8.5f %8.5f %2d\n",
		   theQuads->quadNo,theSlits->ifuSlitNo,theFib->fibNo,
		   theFib->fiberL,theFib->fiberM,theFib->fiberx,theFib->fibery,
		   theFib->fiberTrans,theFib->fiberPwidth,theFib->sigmaY,
		   theFib->sigmaYGroup);

           theFib = theFib->next; 
          }
        theSlits = theSlits->next;
       }
     theQuads = theQuads->next; 
    }
  fclose(ofp);

  return(VM_TRUE);
}


VimosBool readFitsIfuTable(VimosIfuTable *ifuTable, fitsfile *fptr)
{
  int    i, ii;
  int    nCols;
  int    nRows;
  int    dummyInt;
  int    null;
  int    j;

  int    numFibCol, lfibCol, mfibCol, quadCol, rowCol;
  int    xfibCol, yfibCol, pwidthCol, transCol;

  int    newQuadNum;
  int    newSlitNum, slitNum;
  int    fibNum;

  int    status;
  char **ttype;
  char   comment[80];

  VimosIfuQuad *ifuQuad;
  VimosIfuQuad *lastIfuQuad;
  VimosIfuSlit *ifuSlit;
  VimosIfuSlit *lastIfuSlit;
  VimosIfuFiber *ifuFiber;
  VimosIfuFiber *lastIfuFiber;

  status = 0;

  /* validate input */
  if (ifuTable == NULL) {
    pilMsgError("readFitsIfuTable","NULL input table");
    return (VM_FALSE);
  }
  
  /* validate input */
  if ( strcmp(ifuTable->name, VM_IFU) ) {
    pilMsgError("readFitsIfuTable","Invalid input table");
    return (VM_FALSE);
  }
  
  /* open Table */
  if (fits_movnam_hdu(fptr, BINARY_TBL, "IFU", 0, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  ifuTable->fptr = fptr;

  /* read Table */
  if (!readDescsFromFitsTable(&(ifuTable->descs), ifuTable->fptr)) {
    pilMsgError("readFitsIfuTable","The function readDescsFromFitsTable has returned an error");
    return (VM_FALSE);
  }
  if (!readIntDescriptor(ifuTable->descs, "TFIELDS", &nCols, comment)) {
    pilMsgError("readFitsIfuTable","The function readIntDescriptor has returned an error");
    return (VM_FALSE);
  }
  if (!readIntDescriptor(ifuTable->descs, "NAXIS2", &nRows, comment)) {
    pilMsgError("readFitsIfuTable","The function readIntDescriptor has returned an error");
    return (VM_FALSE);
  }

 /* allocate space for the column labels */
  ttype = (char **) cpl_malloc(nCols*sizeof(char *));
  for (i = 0; i < nCols; i++) {
    ttype[i] = (char *) cpl_malloc(FLEN_VALUE*sizeof(char));

    /* check if space was allocated */
    if (ttype[i] == NULL) {
      pilMsgError("readFitsIfuTable","Allocation Error");
      return (VM_FALSE);
    } 
  }

  /* read the column names from the TTYPEn keywords */
  if (fits_read_keys_str(ifuTable->fptr, "TTYPE", 1, nCols, ttype, &dummyInt, 
			 &status)) {
    pilMsgError("readFitsIfuTable","The function fits_read_keys_str has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  /*******************************/
  /* ALEX: change "X" to "XPIX" (marcos's format) - 2-Apr-02 */
  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "XPIX", &xfibCol, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  /*******************************/
  /* ALEX: change "Y" to "YPIX" (marcos's format) - 2-Apr-02 */  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "YPIX", &yfibCol, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "L", &lfibCol, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "M", &mfibCol, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "PWIDTH", &pwidthCol, 
		      &status)) {
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "QUAD", &quadCol, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "ROW", &rowCol, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "FIB", &numFibCol, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_get_colnum(ifuTable->fptr, CASEINSEN, "TRANS", &transCol, &status)){
    pilMsgError("readFitsIfuTable","The function fits_get_colnum has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  if (fits_read_col_int(ifuTable->fptr, quadCol, 1, 1, 1, null, &newQuadNum,
			&null, &status)) {
    pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  lastIfuQuad = NULL;

  i=0;

  while (i< nRows)
    {
      /* create new quadrant */
     ifuQuad = newIfuQuad();
     if (ifuQuad == NULL) {
       pilMsgError("readFitsIfuTable","The function newIfuQuad has returned NULL");
       return (VM_FALSE);
     }

     lastIfuSlit = NULL;
     slitNum = 0;
 
     /* set the quadrant number */
     ifuQuad->quadNo = newQuadNum;

     /* get the number of the first slit in the quadrant */
     if (fits_read_col_int(ifuTable->fptr, rowCol, i+1, 1, 1, null, 
			   &newSlitNum, &null, &status)) {
       pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
       return (VM_FALSE);
     }

     /* while we are in the same quadrant */
     while ((i < nRows) && (newQuadNum == ifuQuad->quadNo) )
       {
 
         /* create new slit */
        ifuSlit = newIfuSlit();
        if (ifuSlit == NULL) {
          pilMsgError("readFitsIfuTable","The function newIfuSlit has returned NULL");
          return(VM_FALSE);
        }

        lastIfuFiber = NULL;
 
        /* set the slit number */
        ifuSlit->ifuSlitNo = newSlitNum;

        /* number of fibers in this slit now is initialized to zero */
        fibNum = 0; 
 
        /* count how many fibers are in the slit */
        while ((i+fibNum < nRows) && (newSlitNum == ifuSlit->ifuSlitNo) )
          {
           if (fits_read_col_int(ifuTable->fptr, rowCol, i+fibNum+1, 1, 1, 
				 null, &newSlitNum, &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
             return (VM_FALSE);
            }

           fibNum++;
          }
 
        if (i+fibNum <nRows)
          {
           fibNum--;
          }

        /* read all the fibers in the slit */
        for (j=0; j<fibNum; j++)
          {
           /* create new fiber */
           ifuFiber = newIfuFiber();
	   if (ifuFiber == NULL) {
	     pilMsgError("readFitsIfuTable","The function newIfuFiber has returned NULL");
	     return (VM_FALSE);
	   }
 
           if (fits_read_col_flt(ifuTable->fptr, xfibCol, i+j+1, 1, 1, null,
				 &(ifuFiber->fiberx), &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_flt has returned an error (code %d)", status);
             return (VM_FALSE);
           }
           if (fits_read_col_flt(ifuTable->fptr, yfibCol, i+j+1, 1, 1, null,
				 &(ifuFiber->fibery),& null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_flt has returned an error (code %d)", status);
             return (VM_FALSE);
           } 
           if (fits_read_col_int(ifuTable->fptr, lfibCol, i+j+1, 1, 1, null,
				 &(ifuFiber->fiberL), &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
             return (VM_FALSE);
           }
           if (fits_read_col_int(ifuTable->fptr, mfibCol, i+j+1, 1, 1, null,
				 &(ifuFiber->fiberM), &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
             return (VM_FALSE);
           }
           if (fits_read_col_flt(ifuTable->fptr, pwidthCol, i+j+1, 1, 1, null,
				 &(ifuFiber->fiberPwidth), &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_flt has returned an error (code %d)", status);
             return (VM_FALSE);
           } 
           if (fits_read_col_int(ifuTable->fptr, quadCol, i+j+1, 1, 1, null,
				 &(ifuQuad->quadNo), &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
             return (VM_FALSE);
           } 
           if (fits_read_col_int(ifuTable->fptr, rowCol, i+j+1, 1, 1, null,
				 &(ifuSlit->ifuSlitNo), &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
             return (VM_FALSE);
           } 
           if (fits_read_col_int(ifuTable->fptr, numFibCol, i+j+1, 1, 1, null,
				 &(ifuFiber->fibNo), &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
             return (VM_FALSE);
           } 
           if (fits_read_col_flt(ifuTable->fptr, transCol, i+j+1, 1, 1, null,
				 &(ifuFiber->fiberTrans), &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_flt has returned an error (code %d)", status);
             return (VM_FALSE);
           }

           /* start linking fibers to the linked list */
           /* if first fiber */
           if (lastIfuFiber == NULL)
             {
              lastIfuFiber = ifuFiber;
              ifuSlit->fibers = ifuFiber;
             }
           /* if not */
           else
             {
              lastIfuFiber->next = ifuFiber;
              ifuFiber->prev = lastIfuFiber;
             }
           /* store tail */
           lastIfuFiber = ifuFiber;
 
          } 
 
        /* jump to the next slit */
        i += fibNum;
 
        slitNum++;
 
        /* start linking slits to the linked list */
        /* if first slit */
        if (lastIfuSlit == NULL)
          {
           lastIfuSlit = ifuSlit;
           ifuQuad->ifuSlits = ifuSlit;
          }
        /* if not */
        else
          {
           lastIfuSlit->next = ifuSlit;
           ifuSlit->prev = lastIfuSlit;
          }
        /* store tail */
        lastIfuSlit = ifuSlit;
 
        /* jump to next quadrant ONLY if you are not at the end of the table */
        if (i<nRows)
          {
           if (fits_read_col_int(ifuTable->fptr, quadCol, i+1, 1,1 , null,
				 &newQuadNum, &null, &status)) {
             pilMsgError("readFitsIfuTable","The function fits_read_col_int has returned an error (code %d)", status);
             return (VM_FALSE);
           }
          }
       }
 
     /* start linking quadrants to the linked list */
     /* if first quadrant */
     if (lastIfuQuad == NULL)
       {
        lastIfuQuad = ifuQuad;
        ifuTable->quads = ifuQuad;
       }
     /* if not */
     else
       {
        lastIfuQuad->next = ifuQuad;
        ifuQuad->prev = lastIfuQuad;
       }
     /* store tail */
     lastIfuQuad = ifuQuad;
 
    }

  /* free the memory for the column labels */
  for (ii = 0; ii < 3; ii++)
     {
      cpl_free( ttype[ii] );
     }
 
  return (VM_TRUE);
}


VimosBool writeFitsIfuTable(VimosIfuTable *ifuTable, fitsfile *fptr)
{
  int i;
  int nRows;
  int rowNum;
  int objNum;

/*  int numFibCol,lfibCol,mfibCol,quadCol,rowCol; */
/*  int xfibCol,yfibCol,pwidthCol,transCol; */

  int status;
  char *ttype[84], *tform[84];
 
  VimosIfuQuad *ifuQuad;
  VimosIfuSlit *ifuSlit;
  VimosIfuFiber *ifuFiber;

  /* validate input */
  if (ifuTable == NULL) {
    pilMsgError("writeFitsIfuTable","NULL input table");
    return (VM_FALSE);
  }

  /* validate input */
  if ( strcmp(ifuTable->name, VM_IFU) ) {
    pilMsgError("writeFitsIfuTable","Invalid input table");
    return (VM_FALSE);
  }

  nRows = ifuTable->numIfuFibs;
 
  status = 0;
  ifuTable->fptr = fptr;

  /* if Table is already present, first remove it  */
  if (!fits_movnam_hdu(fptr, BINARY_TBL, "IFU", 0, &status) ) {
    if (fits_delete_hdu(fptr, NULL, &status)) {
      pilMsgError("writeFitsIfuTable","The function fits_delete_hdu has returned an error (code %d)", status);
      return (VM_FALSE);
    }
  } else {
    status = 0;
  }

  /* allocate space for the column labels */
  for (i = 0; i <= 8; i++)
     {
      ttype[i] = (char *) cpl_malloc(FLEN_VALUE);  /* max label length = 69 */
      /* check if space was allocated */
      if (ttype[i] == NULL) {
	pilMsgError("writeFitsIfuTable","Allocation Error");
	return (VM_FALSE);
      }
      tform[i] = (char *) cpl_malloc(FLEN_VALUE);
      /* check if space was allocated */
      if (tform[i] == NULL) {
	pilMsgError("writeFitsIfuTable","Allocation Error");
	return (VM_FALSE);
      }
     }
  ttype[0] = "L";
  tform[0] = "1J";
  ttype[1] = "XPIX";
  tform[1] = "1E";
  ttype[2] = "YPIX";
  tform[2] = "1E";
  ttype[3] = "M";
  tform[3] = "1J";
  ttype[4] = "PWIDTH";
  tform[4] = "1E";
  ttype[5] = "QUAD";
  tform[5] = "1J";
  ttype[6] = "ROW";
  tform[6] = "1J";
  ttype[7] = "FIB";
  tform[7] = "1J";
  ttype[8] = "TRANS";
  tform[8] = "1E";    
/*    ttype[9] = "SIGMAY"; */
/*    tform[9] = "1E";   */
/*      ttype[10] = "SIGMAYGROUP"; */
/*    tform[10] = "1E";     */
/*    ttype[11] = "X_IMA"; */
/*    tform[11] = "1J";   */
/*      ttype[12] = "Y_IMA"; */
/*    tform[12] = "1J";    */ 
  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl(fptr, BINARY_TBL, 0, 9, ttype, tform, NULL, 
		       "IFU", &status)) {
    pilMsgError("writeFitsIfuTable","The function fits_create_tbl has returned an error (code %d)", status);
    return (VM_FALSE);
  }
  if (fits_movnam_hdu(fptr, BINARY_TBL, "IFU", 0, &status)) {
    pilMsgError("writeFitsIfuTable","The function fits_movnam_hdu has returned an error (code %d)", status);
    return (VM_FALSE);
  }

  /* write the table descriptors to the Fits file */
  
  /* BG: if IN the table there are the NAXIS, NAXIS1 and NAXIS2
     descriptors and the various fits standard table descriptors
     (e.g. because it has been created starting from 
     another table read from disk) the write fails with error 241
     NAXIS descriptors must be removed before, if existing 
     Remove also any already existing TTYPE, TFORM, TUNIT. 
     These will be correctly created by cfitsio according to current table
  */

  deleteSetOfDescriptors(&(ifuTable->descs),"NAXIS*");
  deleteSetOfDescriptors(&(ifuTable->descs),"*COUNT");
  deleteSetOfDescriptors(&(ifuTable->descs),"TUNIT*");
  deleteSetOfDescriptors(&(ifuTable->descs),"TFIELDS*");
  deleteSetOfDescriptors(&(ifuTable->descs),"TTYPE*");
  deleteSetOfDescriptors(&(ifuTable->descs),"TFORM*");

/*    if (fits_read_key (ifuTable->fptr,TINT,"NAXIS1",&nbytes,NULL,&status)) { */
/*      pilMsgError("writeFitsIfuTable","The function fits_read_key has returned an error (code %d)", status); */
/*      return (VM_FALSE); */
/*    } */


/*    if (!writeIntDescriptor(&(ifuTable->descs), "NAXIS1", nbytes, "")) { */
/*      pilMsgError("writeFitsIfuTable","The function writeIntDescriptor has returned an error"); */
/*      return (VM_FALSE); */
/*    } */
/*    if (!writeIntDescriptor(&(ifuTable->descs), "NAXIS2", nRows, "")) { */
/*      pilMsgError("writeFitsIfuTable","The function writeIntDescriptor has returned an error"); */
/*      return (VM_FALSE); */
/*    } */
/*    if (!writeIntDescriptor(&(ifuTable->descs),"TFIELDS", 9, "")) { */
/*      pilMsgError("writeFitsIfuTable","The function writeIntDescriptor has returned an error"); */
/*      return (VM_FALSE); */
/*    } */
  if (!writeDescsToFitsTable(ifuTable->descs, ifuTable->fptr)) {
    pilMsgError("writeFitsIfuTable","The function writeDescsToFitsTable has returned an error");
    return (VM_FALSE);
  }

  /* get quadrants in the table */
  ifuQuad = ifuTable->quads;
 
  rowNum = 1;
 
  while (ifuQuad)
    {
      /* get slits in the quadrant */
     ifuSlit = ifuQuad->ifuSlits;
 
     while (ifuSlit)
       {
        /* get fibers in the slit */
        ifuFiber = ifuSlit->fibers;
        objNum = 1;
 
        /* write all fibers in the slit */
        while (ifuFiber)
          {
           if (fits_write_col_int(ifuTable->fptr, 1, rowNum, 1, 1, 
                                  &(ifuFiber->fiberL),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_int has returned an error (code %d)", status);
             return (VM_FALSE);
           }
           if (fits_write_col_flt(ifuTable->fptr, 2, rowNum, 1, 1, 
                                  &(ifuFiber->fiberx),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_flt has returned an error (code %d)", status);
             return (VM_FALSE);   
           }
           if (fits_write_col_flt(ifuTable->fptr, 3, rowNum, 1, 1, 
                                  &(ifuFiber->fibery),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_flt has returned an error (code %d)", status);
             return (VM_FALSE);  
           }
           if (fits_write_col_int(ifuTable->fptr, 4, rowNum, 1, 1, 
                                  &(ifuFiber->fiberM),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_int has returned an error (code %d)", status);
             return (VM_FALSE);   
           }
           if (fits_write_col_flt(ifuTable->fptr, 5, rowNum, 1, 1, 
                                  &(ifuFiber->fiberPwidth),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_flt has returned an error (code %d)", status);
             return (VM_FALSE);   
           }
           if (fits_write_col_int(ifuTable->fptr, 6, rowNum, 1, 1, 
                                  &(ifuQuad->quadNo),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_int has returned an error (code %d)", status);
             return (VM_FALSE);   
           }
           if (fits_write_col_int(ifuTable->fptr, 7, rowNum, 1, 1, 
                                  &(ifuSlit->ifuSlitNo),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_int has returned an error (code %d)", status);
             return (VM_FALSE);   
           }
           if (fits_write_col_int(ifuTable->fptr, 8, rowNum, 1, 1, 
                                  &(ifuFiber->fibNo),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_int has returned an error (code %d)", status);
             return (VM_FALSE);   
           }
           if (fits_write_col_flt(ifuTable->fptr, 9, rowNum, 1, 1, 
                                  &(ifuFiber->fiberTrans),&status)) {
             pilMsgError("writeFitsIfuTable","The function fits_write_col_flt has returned an error (code %d)", status);          
             return (VM_FALSE);    
           }


/*             if (fits_write_col_flt(ifuTable->fptr, 10, rowNum, 1, 1,  */
/*                                    &(ifuFiber->sigmaY),&status)) { */
/*               pilMsgError("writeFitsIfuTable","The function fits_write_col_flt has returned an error (code %d)", status);           */
/*               return (VM_FALSE);     */
/*             } */
/*             if (fits_write_col_flt(ifuTable->fptr, 11, rowNum, 1, 1,  */
/*                                    &(ifuFiber->sigmaYGroup),&status)) { */
/*               pilMsgError("writeFitsIfuTable","The function fits_write_col_flt has returned an error (code %d)", status);           */
/*               return (VM_FALSE);     */
/*             } */

/*             if (fits_write_col_flt(ifuTable->fptr, 12, rowNum, 1, 1,  */
/*                                    &(ifuFiber->x_image),&status)) { */
/*               pilMsgError("writeFitsIfuTable","The function fits_write_col_flt has returned an error (code %d)", status);           */
/*               return (VM_FALSE);     */
/*             } */

/*             if (fits_write_col_flt(ifuTable->fptr, 13, rowNum, 1, 1,  */
/*                                    &(ifuFiber->fiberTrans),&status)) { */
/*               pilMsgError("writeFitsIfuTable","The function fits_write_col_flt has returned an error (code %d)", status);           */
/*               return (VM_FALSE);     */
/*             } */
           rowNum++;
           objNum++;

           /* get next fiber */
           ifuFiber = ifuFiber->next;
 
          }
 
        /* get next slit */
        ifuSlit = ifuSlit->next;
 
       }
 
     /* get next quadrant */
     ifuQuad = ifuQuad->next;
 
    }
 
  return(VM_TRUE);
}
