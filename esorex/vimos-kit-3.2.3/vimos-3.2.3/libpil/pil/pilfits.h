/* $Id: pilfits.h,v 1.2 2011-01-31 08:30:46 cizzo Exp $
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
 * $Date: 2011-01-31 08:30:46 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef _PIL_FITS_H
#define _PIL_FITS_H

#include <fitsio.h>

#include <pilmacros.h>


PIL_BEGIN_DECLS

/*
 * These definitions already take into account the trailing zero
 * terminating the strings.
 */

#define PIL_FITS_CARD_MAX   (FLEN_CARD)
#define PIL_FITS_VALUE_MAX  (FLEN_VALUE)
 

typedef enum _PIL_FITS_IO_MODE_ PilFitsIOMode;

enum _PIL_FITS_IO_MODE_ {
    PIL_FITS_READ,
    PIL_FITS_WRITE,
    PIL_FITS_READWRITE
};

typedef struct _PIL_FITS_FILE_ PilFitsFile;

/*
 * Constructors and Destructors
 */

PilFitsFile *newPilFitsFile(const char *, PilFitsIOMode);
void deletePilFitsFile(PilFitsFile *);

/*
 * Methods
 */

int pilFitsHdrCount(PilFitsFile *);
int pilFitsHdrGoto(PilFitsFile *, int);

int pilFitsHdrReadInt(PilFitsFile *, const char *, int *);
int pilFitsHdrReadLogical(PilFitsFile *, const char *, int *);
int pilFitsHdrReadFloat(PilFitsFile *, const char *, float *);
int pilFitsHdrReadDouble(PilFitsFile *, const char *, double *);
int pilFitsHdrReadString(PilFitsFile *, const char *, char **);

int pilFitsHdrWriteInt(PilFitsFile *, const char *, 
                   int, const char *);
int pilFitsHdrWriteLogical(PilFitsFile *, const char *, 
                   int, const char *);
int pilFitsHdrWriteFloat(PilFitsFile *, const char *, 
                   float, const char *);
int pilFitsHdrWriteDouble(PilFitsFile *, const char *, 
                   double, const char *);
int pilFitsHdrWriteString(PilFitsFile *, const char *, 
                   const char *, const char *);

int pilFitsHdrInsertInt(PilFitsFile *, int, const char *, const char *, 
                    int, const char *);
int pilFitsHdrInsertLogical(PilFitsFile *, int, const char *, const char *, 
                    int, const char *);
int pilFitsHdrInsertFloat(PilFitsFile *, int, const char *, const char *, 
                    float, const char *);
int pilFitsHdrInsertDouble(PilFitsFile *, int, const char *, const char *, 
                    double, const char *);
int pilFitsHdrInsertString(PilFitsFile *, int, const char *, const char *, 
                    const char *, const char *);

int pilFitsHdrWriteDate(PilFitsFile *);
int pilFitsHdrInsertDate(PilFitsFile *, int, const char *);

int pilFitsHdrReadCard(PilFitsFile *, const char *, char **);
int pilFitsHdrWriteCard(PilFitsFile *, const char *, const char *);
int pilFitsHdrInsertCard(PilFitsFile *, int, const char *, const char *);

int pilFitsHdrWriteComment(PilFitsFile *, const char *);
int pilFitsHdrWriteHistory(PilFitsFile *, const char *);

int pilFitsHdrDelete(PilFitsFile *, const char *);
int pilFitsHdrDeleteKeys(const char *, const char *, int);

/*
 * Miscellaneous utilities
 */

int pilFitsHdrCopy(const char *, unsigned int, const char *, const char *,
                   unsigned int);

char *pilFitsMD5Signature(const char *);

PIL_END_DECLS

#endif /* _PIL_FITS_H */
