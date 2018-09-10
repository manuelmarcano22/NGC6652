/* $Id: pilpaf.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
 
#ifndef _PIL_PAF_H
#define _PIL_PAF_H

#include <sys/types.h>

#include <pilmacros.h>


PIL_BEGIN_DECLS

/*
 * Maximum length of a parameter file record, i.e. maximum number of
 * characters per line of a parameter file on disk. This does not include
 * a trailing 0.
 */

#define PAF_RECORD_MAX  (256)


/*
 * PAF value types
 */

enum _PIL_PAF_TYPE_ {
    PAF_TYPE_NONE,
    PAF_TYPE_BOOL,
    PAF_TYPE_INT,
    PAF_TYPE_DOUBLE,
    PAF_TYPE_STRING
};

typedef enum _PIL_PAF_TYPE_ PilPAFType;

/*
 * PAF object
 */

typedef struct _PIL_PAF_ PilPAF;

/*
 * Create, copy and destroy operations
 */

PilPAF *newPilPAF(const char *, const char *, const char *,
                          const char *);
void deletePilPAF(PilPAF *);

/*
 * Nonmodifying operations
 */

int pilPAFIsEmpty(const PilPAF *);
size_t pilPAFGetSize(const PilPAF *);
int pilPAFContains(const PilPAF *, const char *);
size_t pilPAFCount(const PilPAF *, const char *);

/*
 * Header operations
 */

PilPAFType pilPAFType(const PilPAF *, const char *);

const char *pilPAFGetName(const PilPAF *);
const char *pilPAFGetTag(const PilPAF *);
const char *pilPAFGetId(const PilPAF *);
const char *pilPAFGetDescription(const PilPAF *);

int pilPAFSetName(PilPAF *, const char *);
int pilPAFSetTag(PilPAF *, const char *);
int pilPAFSetId(PilPAF *, const char *);
int pilPAFSetDescription(PilPAF *, const char *);

int pilPAFSetHeader(PilPAF *, const char *, const char *, const char *,
                    const char *);

/*
 * Element access
 */

int pilPAFGetValueBool(const PilPAF *, const char *);
int pilPAFGetValueInt(const PilPAF *, const char *);
double pilPAFGetValueDouble(const PilPAF *, const char *);
const char *pilPAFGetValueString(const PilPAF *, const char *);
const char *pilPAFGetComment(const PilPAF *, const char *);

int pilPAFSetValueBool(PilPAF *, const char *, int);
int pilPAFSetValueInt(PilPAF *, const char *, int);
int pilPAFSetValueDouble(PilPAF *, const char *, double);
int pilPAFSetValueString(PilPAF *, const char *, const char *);
int pilPAFSetComment(PilPAF *, const char *, const char *);

/*
 * Inserting and removing elements
 */

int pilPAFInsertBool(PilPAF *, const char *, const char *, int, const char *);
int pilPAFInsertInt(PilPAF *, const char *, const char *, int, const char *);
int pilPAFInsertDouble(PilPAF *, const char *, const char *, double,
                       const char *);
int pilPAFInsertString(PilPAF *, const char *, const char *, const char *,
                       const char *);

int pilPAFInsertAfterBool(PilPAF *, const char *, const char *, int,
                          const char *);
int pilPAFInsertAfterInt(PilPAF *, const char *, const char *, int,
                         const char *);
int pilPAFInsertAfterDouble(PilPAF *, const char *, const char *, double,
                            const char *);
int pilPAFInsertAfterString(PilPAF *, const char *, const char *,
                            const char *, const char *);

int pilPAFPrependBool(PilPAF *, const char *, int, const char *);
int pilPAFPrependInt(PilPAF *, const char *, int, const char *);
int pilPAFPrependDouble(PilPAF *, const char *, double, const char *);
int pilPAFPrependString(PilPAF *, const char *, const char *, const char *);

int pilPAFAppendBool(PilPAF *, const char *, int, const char *);
int pilPAFAppendInt(PilPAF *, const char *, int, const char *);
int pilPAFAppendDouble(PilPAF *, const char *, double, const char *);
int pilPAFAppendString(PilPAF *, const char *, const char *, const char *);

void pilPAFErase(PilPAF *, const char *);
void pilPAFClear(PilPAF *);

/*
 * Read and write operations
 */

int pilPAFWrite(PilPAF *);

/*
 * Miscellaneous utilities
 */

int pilPAFIsValidName(const char *);

PIL_END_DECLS

#endif /* _PIL_PAF_H */
