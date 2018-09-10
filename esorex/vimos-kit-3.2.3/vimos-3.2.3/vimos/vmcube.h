/* $Id: vmcube.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_CUBE_H
#define VM_CUBE_H

#include <fitsio.h>

#include <pilmacros.h>

#include <vmimage.h>
#include <vmifutable.h>
#include <vmobjecttable.h>


PIL_BEGIN_DECLS

/* 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  stucture VimosCube

  Description:
  Definition of Cube structure. Similar to VimosImage but 3D

  Layout:
     int   xlen
     int   ylen
     int   zlen
     float *data;
     VimosDescriptor *descs
     fitsfile *fptr


   Updates:
   22 Nov 00: Created (AZ)

-------------------------------------------------------------------------------
*/


typedef struct _VIMOS_CUBE_
{
  int xlen;
  int ylen;
  int zlen;
  float *data;
  VimosDescriptor *descs;
  fitsfile *fptr;
} VimosCube;




VimosCube *newCube(int xlen, int ylen, int zlen, float *data);
VimosCube *newCubeAndAlloc(int xlen, int ylen, int zlen);

void deleteCube(VimosCube *Cube);
void deleteCubeAndAlloc(VimosCube *Cube);


/*
  Find (L,M) coordinates of the fiber corresponding to an Object
  in an objectTable
*/

VimosFloatArray *selectFiberForObject(VimosIfuSlit *ifuSlit, 
                                      VimosObjectObject *object,
                                      float *theData, int specLen, 
                                      int objNum, int *theL, int *theM);

VimosBool  openNewFitsCube(char *cubeName, VimosCube *cube);

VimosBool closeFitsCube(VimosCube *cube, int flag);

PIL_END_DECLS

#endif /* VM_CUBE_H */
