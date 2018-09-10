/* $Id: pilframeset.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifndef PIL_FRAMESET_H
#define PIL_FRAMESET_H

#include <stdio.h>

#include <pilmacros.h>
#include <pilframe.h>
#include <pildictionary.h>


PIL_BEGIN_DECLS

typedef PilDictionary PilSetOfFrames;
typedef PilDictCapacity PilSetCapacity;

/*
 * Constructors and Destructors
 */

PilSetOfFrames *newPilSetOfFrames(void);
void deletePilSetOfFrames(PilSetOfFrames *);


/*
 * Methods
 */

int pilSofInsert(PilSetOfFrames *, PilFrame *);
void pilSofRemove(PilSetOfFrames *, const char *);

PilFrame *pilSofLookup(PilSetOfFrames *, const char *);
PilFrame *pilSofLookupNext(PilSetOfFrames *, const char *);

PilSetCapacity pilSofFrameCount(PilSetOfFrames *, const char *);
int pilSofIsEmpty(const PilSetOfFrames *);

PilFrame *pilSofFirst(PilSetOfFrames *);
PilFrame *pilSofNext(PilSetOfFrames *, PilFrame *);

PilSetOfFrames *pilSofRead(const char *, PilSetOfFrames *);
int pilSofWrite(PilSetOfFrames *set, const char *);
int pilSofDump(FILE *stream, const char, PilSetOfFrames *);

PIL_END_DECLS

#endif /* PIL_FRAMESET_H */
