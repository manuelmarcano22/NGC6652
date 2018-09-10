/* $Id: vmtypes.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>

#include "vmtypes.h"
#include "cpl.h"


/*

  This constructor has, contrary to the other construtors of Vimos Linked
  Lists, an argument that allows to allocate the linked list as an array.  The
  reason for this 'inconsistency' is that VimosDpoint is taken from the
  eclipse library, where dpoints are allocated as arrays. Just to be safe, and
  not run into possible conflicts with eclipse, I allow for VimosDpoints as an
  array.  

*/

VimosDpoint *newDpoint(int number) 
{
  const char   modName[] = "newDpoint";
  int          i;
  VimosDpoint *tPoint;

  if (number <= 0) {
    cpl_msg_error(modName, "Invalid argument");
    return(NULL);
  }
  
  
  tPoint = (VimosDpoint *) cpl_calloc(number, sizeof(VimosDpoint));
  
  /* check if space was allocated */
  if (tPoint == NULL) {
    cpl_msg_error(modName, "Allocation error");
    return(NULL);
  }

  /* initialize links */
  switch (number) {
  case 1: 
    {
      tPoint->prev = tPoint->next = NULL;
      break;
    }
  default: 
    {
      for (i = 1; i < number-1; i++) {
        tPoint[i].next = &(tPoint[i+1]);
        tPoint[i].prev = &(tPoint[i-1]);
      }
      tPoint[0].prev = NULL;
      tPoint[0].next = &(tPoint[1]);
      tPoint[number-1].prev = &(tPoint[number-2]);
      tPoint[number-1].next = NULL;
      break;
    }
  }
    
  return(tPoint);
}

void deleteDpoint(VimosDpoint *dPoint) 
{
  if (dPoint) cpl_free(dPoint);
}

/*

  This constructor has, contrary to the other construtors of Vimos Linked
  Lists, an argument that allows to allocate the linked list as an array.  The
  reason for this 'inconsistency' is that VimosPixel is taken from the eclipse
  library (where it exists as pixel_position), where they are allocated as
  arrays. Just to be safe, and not run into possible conflicts with eclipse, I
  allow for VimosPixels as an array.

*/

VimosPixel *newPixel(int number) 
{
  const char  modName[] = "newPixel";
  int         i;
  VimosPixel *tPixel;
  
  if (number <= 0) {
    cpl_msg_error(modName, "Invalid argument");
    return(NULL);
  }
  
  tPixel = (VimosPixel *) cpl_calloc(number, sizeof(VimosPixel));

  /* check if space was allocated */  
  if (tPixel == NULL) {
    cpl_msg_error(modName, "Allocation Error");
    return(NULL);
  }


  /* initialize links */
  switch (number) {
  case 1: 
    {
      tPixel->prev = tPixel->next = NULL;
      break;
    }
  default: 
    {
      for (i = 1; i < number-1; i++) {
        tPixel[i].next = &(tPixel[i+1]);
        tPixel[i].prev = &(tPixel[i-1]);
      }
      tPixel[0].prev = NULL;
      tPixel[0].next = &(tPixel[1]);
      tPixel[number-1].prev = &(tPixel[number-2]);
      tPixel[number-1].next = NULL;
      break;
    }
  }    
  
  return(tPixel);
}

void deletePixel(VimosPixel *pixel) 
{
  if (pixel) cpl_free(pixel);
}
