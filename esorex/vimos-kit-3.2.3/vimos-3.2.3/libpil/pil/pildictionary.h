/* $Id: pildictionary.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifndef _PIL_DICTIONARY_H
#define _PIL_DICTIONARY_H

#include <dict.h>

#include <pilmacros.h>

 
PIL_BEGIN_DECLS

#define PIL_DICT_CAPACITY_MAX  (DICTCOUNT_T_MAX)

typedef dictcount_t PilDictCapacity;

typedef dict_t PilDictionary;
typedef dnode_t PilDictNode;

typedef dict_comp_t PilDictComparator;
typedef dnode_alloc_t PilDictAllocator;
typedef dnode_free_t PilDictDeallocator;
typedef dnode_process_t PilDictIterator;


/*
 * Node constructors and destructors
 */

extern PilDictNode *newPilDictNode(void *);
extern void deletePilDictNode(PilDictNode *);


/*
 * Node methods
 */

extern const void *pilDictGetKey(PilDictNode *);
extern void *pilDictGetData(PilDictNode *);
extern void pilDictPutData(PilDictNode *, void *);


/*
 * Dictionary constructors and destructors
 */

extern PilDictionary *newPilDictionary(PilDictCapacity, PilDictComparator);
extern void deletePilDictionary(PilDictionary *);


/*
 * Dictionary methods
 */

extern void pilDictAllowDuplicates(PilDictionary *);
extern void pilDictSetAllocator(PilDictionary *, PilDictAllocator,
				PilDictDeallocator, void *);

extern int pilDictVerify(PilDictionary *);
extern int pilDictIsEmpty(PilDictionary *);
extern int pilDictIsFull(PilDictionary *);
extern PilDictCapacity pilDictCapacity(PilDictionary *);

extern PilDictNode *pilDictLookup(PilDictionary *, const void *);
extern int pilDictContains(PilDictionary *, PilDictNode *);

extern int pilDictInsert(PilDictionary *, const void *, void *);
extern void pilDictInsertNode(PilDictionary *, PilDictNode *, const void *);

extern PilDictNode *pilDictRemove(PilDictionary *, PilDictNode *);
extern void pilDictErase(PilDictionary *, PilDictNode *);
extern void pilDictClear(PilDictionary *);

extern PilDictNode *pilDictBegin(PilDictionary *);
extern PilDictNode *pilDictEnd(PilDictionary *);
extern PilDictNode *pilDictNext(PilDictionary *, PilDictNode *);
extern PilDictNode *pilDictPrev(PilDictionary *, PilDictNode *);

PIL_END_DECLS

#endif /* _PIL_DICTIONARY_H */
