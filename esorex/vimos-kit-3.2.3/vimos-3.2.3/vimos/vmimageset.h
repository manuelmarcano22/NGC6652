/* $Id: vmimageset.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_IMAGESET_H
#define VM_IMAGESET_H

#include <pilmacros.h>

#include <vmimage.h>
#include <vmccdtable.h>
#include <vmifutable.h>
#include <vmwindowtable.h>
#include <vmextractiontable.h>
#include <vmobjecttable.h>


PIL_BEGIN_DECLS

struct _VIMOS_SINGLE_IMAGE_ {
  VimosImage *theImage;

  VimosTable *ccdTable;
  VimosTable *sphotStdTable;
  VimosIfuTable *ifuTable;
  VimosExtractionTable *extractionTable;
  VimosWindowTable *windowTable;
  VimosObjectTable *objectTable;

  char *imageName;

  struct _VIMOS_SINGLE_IMAGE_ *next;
  struct _VIMOS_SINGLE_IMAGE_ *prev;
};

typedef struct _VIMOS_SINGLE_IMAGE_ VimosSingleImage;


struct _VIMOS_IMAGE_SET_ {
  VimosSingleImage *images;
  int nImas;
};

typedef struct _VIMOS_IMAGE_SET_ VimosImageSet;


VimosSingleImage *newSingleImage();
void deleteSingleImage(VimosSingleImage *oneImage);

VimosImageSet *newImageSet();
void deleteImageSet(VimosImageSet *imageSet);

void storeImagesAndTables(VimosImageSet *imageSet, char *imageList,
                          char extname[], int ccdT, int sphT, int ifuT,
                          int extT, int winT, int objT);
                                    
void closeFitsImageSet(VimosImageSet *imageSet, int flag);

PIL_END_DECLS

#endif /* VM_IMAGESET_H */
