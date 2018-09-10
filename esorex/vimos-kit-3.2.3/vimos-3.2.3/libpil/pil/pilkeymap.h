/* $Id: pilkeymap.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifndef _PIL_KEYMAP_H
#define _PIL_KEYMAP_H

#include <pilmacros.h>
#include <pildictionary.h>
#include <pilalias.h>

 
PIL_BEGIN_DECLS

typedef PilDictionary PilKeymap;

/*
 * Constructors and Destructors
 */

PilKeymap *newPilKeymap(void);
void deletePilKeymap(PilKeymap *);


/*
 * Methods
 */

int pilKeymapInsert(PilKeymap *, PilAlias *);
void pilKeymapRemove(PilKeymap *, const char *);
PilAlias *pilKeymapLookup(PilKeymap *, const char *);

const char *pilKeymapGetValue(PilKeymap *, const char *);
const char *pilKeymapGetFormat(PilKeymap *, const char *);
const char *pilKeymapGetComment(PilKeymap *, const char *);

PIL_END_DECLS

#endif /* _PIL_KEYMAP_H */
