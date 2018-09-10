/* $Id: vmsextractor.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_SEXTRACTOR_H
#define VM_SEXTRACTOR_H

#include <stdio.h>
#include <time.h>

#include <pilmacros.h>

#include <vmimage.h>
#include <vmtable.h>


PIL_BEGIN_DECLS

/*
 * Column data types known by SExtractor
 */

enum _SEXT_COLUMN_TYPE_ {
  SEXT_COLUMN_UNDEF = 0,
  SEXT_COLUMN_INT,
  SEXT_COLUMN_FLOAT
};

typedef enum _SEXT_COLUMN_TYPE_ SextColumnType;


/*
 * SExtractor output column specification consists just of the
 * column label anf the format string. Both structure members are
 * constants so that this structure can only be initialized!
 */

struct _SEXT_PARAMETER_ {
  const char *name;
  SextColumnType type;
};

typedef struct _SEXT_PARAMETER_ SextParameter;


/*
 * Methods
 */

const char *sextGetSextractorPath(void);

int sextSaveConfiguration(FILE *, VimosImage *);
int sextSaveParameters(FILE *, SextParameter *);

VimosTable *sextConvertCatalog(const char *, SextParameter *);

const char *sextGetStarNnwName(void);
const char *sextGetFilterName(void);
const char *sextGetAssocName(void);
const char *sextGetCheckImageName(void);
const char *sextGetFlagImageName(void);
const char *sextGetWeightImageName(void);

time_t sextGetExecutionTimeLimit(void);

const char *sextGetFileName(char *buffer, const char *path, size_t n);

PIL_END_DECLS

#endif /* VM_SEXTRACTOR_H */
