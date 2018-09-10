/* $Id: pilframe.h,v 1.2 2012-06-14 15:46:12 cgarcia Exp $
 *
 * This file is part of the VIMOS pipeline library
 * Copyright (C) 2000-2004 European Southern Observatory
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2012-06-14 15:46:12 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef _PIL_FRAME_H
#define _PIL_FRAME_H

#include <pilmacros.h>


PIL_BEGIN_DECLS

/* To be replaced by keymap header file */
#define PIL_CATEGORY_MAX  32

typedef struct _PIL_FRAME_ PilFrame;

enum _PIL_FRAME_TYPE_ {
  PIL_FRAME_TYPE_UNDEF = 0,
  PIL_FRAME_TYPE_RAW,
  PIL_FRAME_TYPE_CALIB,
  PIL_FRAME_TYPE_PRODUCT
};

enum _PIL_FRAME_FORMAT_ {
  PIL_FRAME_FORMAT_UNDEF = 0,
  PIL_FRAME_FORMAT_IMAGE,
  PIL_FRAME_FORMAT_TABLE,
  PIL_FRAME_FORMAT_PAF
};

enum _PIL_PRODUCT_LEVEL_ {
  PIL_PRODUCT_LEVEL_UNDEF = 0,
  PIL_PRODUCT_LEVEL_TEMPORARY,
  PIL_PRODUCT_LEVEL_INTERMEDIATE,
  PIL_PRODUCT_LEVEL_SECONDARY,
  PIL_PRODUCT_LEVEL_PRIMARY
};

enum _PIL_PRODUCT_TYPE_ {
  PIL_PRODUCT_TYPE_UNDEF = 0,
  PIL_PRODUCT_TYPE_TEMPORARY,
  PIL_PRODUCT_TYPE_PREPROCESSED,
  PIL_PRODUCT_TYPE_REDUCED,
  PIL_PRODUCT_TYPE_QCPARAM
};

typedef enum _PIL_FRAME_TYPE_ PilFrameType;
typedef enum _PIL_FRAME_FORMAT_ PilFrameFormat;
typedef enum _PIL_PRODUCT_LEVEL_ PilProductLevel;
typedef enum _PIL_PRODUCT_TYPE_ PilProductType;


/*
 * Identifier strings for the different frame types.
 */

#define PIL_FRAME_TYPE_RAW_ID      "RAW"
#define PIL_FRAME_TYPE_CALIB_ID    "CALIB"
#define PIL_FRAME_TYPE_PRODUCT_ID  "PRODUCT"


/*
 * Constructors and Destructors
 */

PilFrame *newPilFrameEmpty(void);
PilFrame *newPilFrame(const char *, const char *);
void deletePilFrame(PilFrame *);


/*
 * Methods
 */

const char *pilFrmGetName(const PilFrame *);
int pilFrmSetName(PilFrame *, const char *);

const char *pilFrmGetCategory(const PilFrame *);
int pilFrmSetCategory(PilFrame *, const char *);

PilFrameType pilFrmGetType(const PilFrame *);
int pilFrmSetType(PilFrame *, PilFrameType);

PilFrameFormat pilFrmGetFormat(const PilFrame *);
int pilFrmSetFormat(PilFrame *, PilFrameFormat);

PilProductLevel pilFrmGetProductLevel(const PilFrame *);
int pilFrmSetProductLevel(PilFrame *, PilProductLevel);

PilProductType pilFrmGetProductType(const PilFrame *);
int pilFrmSetProductType(PilFrame *, PilProductType);

unsigned int pilFrmGetKeepFlag(const PilFrame *);
int pilFrmSetKeepFlag(PilFrame *);
int pilFrmResetKeepFlag(PilFrame *);

unsigned int pilFrmGetIgnoreFlag(const PilFrame *);
int pilFrmSetIgnoreFlag(PilFrame *);
int pilFrmResetIgnoreFlag(PilFrame *);

PilFrame *pilFrmDuplicate(PilFrame *);

PIL_END_DECLS

#endif /* _PIL_FRAME_H */
