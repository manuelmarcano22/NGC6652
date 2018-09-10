/* $Id: pilrecipe.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
 * $Author: cizzo $
 * $Date: 2008-10-21 09:10:13 $
 * $Revision: 1.1.1.1 $
 * $Name: not supported by cvs2svn $
 */

#ifndef _PIL_RECIPE_H
#define _PIL_RECIPE_H

#include <pilmacros.h>
#include <pilframeset.h>
#include <piltimer.h>


PIL_BEGIN_DECLS

/*
 * Structure defining a recipe parameter.
 */

typedef struct _PIL_REC_PARAMETER_ {
  const char *name;
  const char *type;
  const char *unit;
  const char *value;
  const char *comment;
} PilRecParameter;


/*
 * Pipeline runtime environment startup and shutdown
 */

PilSetOfFrames *pilRecStart(const char *, const char *, const char *,
			    PilRecParameter *, int, char **);
void pilRecStop(PilSetOfFrames *);

int pilRecWriteProducts(PilSetOfFrames *);


/*
 * Reset the recipe information
 */

void pilRecInfoClear(void);


/*
 * Set and retrieve recipe information
 */


int pilRecSetName(const char *);
const char *pilRecGetName(void);

int pilRecSetVersion(const char *);
const char *pilRecGetVersion(void);

int pilRecSetInstrument(const char *);
const char *pilRecGetInstrument(void);

int pilRecSetTimeStart(PilTime);
PilTime pilRecGetTimeStart(void);

int pilRecSetTimeStop(PilTime);
PilTime pilRecGetTimeStop(void);

int pilRecSetTimer(PilTimer *);
PilTimer *pilRecGetTimer(void);


/*
 * Additional functions
 */

void pilRecListSet(PilSetOfFrames *);
int pilRecValidateSet(PilSetOfFrames *set);
int pilRecUpdateProductInfo(PilFrame *, const char *, PilSetOfFrames *);

PIL_END_DECLS

#endif /* _PIL_RECIPE_H */
