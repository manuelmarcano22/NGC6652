/* $Id: pillist.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifndef _PIL_LIST_H
#define _PIL_LIST_H

#include <stddef.h>
#include <list.h>

#include <pilmacros.h>

 
PIL_BEGIN_DECLS

typedef listcount_t PilListSize;
typedef list_t PilList;
typedef lnode_t PilListNode;


/*
 * Node constructors and destructors
 */

PilListNode *newPilListNode(void *);
void deleteListNode(PilListNode *);


/*
 * Node methods
 */

void *pilListNodeGet(const PilListNode *);
void pilListNodePut(PilListNode *, const void *);
int pilListNodeIsListMember(const PilListNode *);


/*
 * List constructors and destructors
 */

PilList *newPilList(void);
void deletePilList(PilList *);
void pilListDestroy(PilList *, void (*)(void *));

/*
 * List methods
 */

/*
 * Nonmodifying list operations
 */

int pilListIsEmpty(const PilList *);
int pilListIsFull(const PilList *);
PilListSize pilListSize(const PilList *);
PilListSize pilListMaxSize(const PilList *);
int pilListContains(const PilList *, const PilListNode *);
int pilListIsSorted(PilList *, int (*)(const void *, const void *));
int pilListVerify(const PilList *);

/*
 * Element access
 */

PilListNode *pilListLookup(PilList *, const void *,
                           int (*)(const void *, const void *));

/*
 * Iterator functions
 */

PilListNode *pilListBegin(const PilList *);
PilListNode *pilListEnd(const PilList *);
PilListNode *pilListNext(const PilList *, const PilListNode *);
PilListNode *pilListPrev(const PilList *, const PilListNode *);

/*
 * Inserting and removing elements
 */

void pilListPushFront(PilList *, PilListNode *);
PilListNode *pilListPopFront(PilList *);
void pilListPushBack(PilList *, PilListNode *);
PilListNode *pilListPopBack(PilList *);

PilListNode *pilListInsert(PilList *, PilListNode *, PilListNode *);
PilListNode *pilListInsertAfter(PilList *, PilListNode *, PilListNode *);
void pilListErase(PilList *, PilListNode *, void (*)(void *));
PilListNode *pilListExtract(PilList *, PilListNode *);
void pilListRemove(PilList *, void *);
void pilListClear(PilList *);

/*
 * Splice functions
 */

void pilListSplice(PilList *, PilListNode *, PilList *, PilListNode *,
                   PilListNode *);
void pilListMerge(PilList *, PilList *, int (*)(const void *, const void *));
void pilListTransfer(PilList *, PilList *, PilListNode *, PilListNode *);

/*
 * Sorting
 */

void pilListSort(PilList *, int(const void *, const void *));

PIL_END_DECLS

#endif /* _PIL_LIST_H */
