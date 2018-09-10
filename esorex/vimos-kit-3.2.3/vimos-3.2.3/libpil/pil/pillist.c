/* $Id: pillist.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>

#include "pillist.h"


 PilListNode *pilListExtract(PilList *, PilListNode *);


/**
 * @defgroup pilList pilList
 *
 * The module @b pilList provides functions to create maintain and
 * destroy a list object.
 *
 * @note
 *   The current implementation is based on the @b kazlib list
 *   implementation. Not all services defined there are available
 *   through this interface.
 */

/**@{*/


/**
 * @brief
 *   Create a new list node.
 *
 * @param data  Pointer to the node data.
 *
 * @return The function returns the pointer to the newly created node if no
 *    error occurs, otherwise the return value is @c NULL.
 *
 * The function allocates the memory for a list node and initializes
 * the node's fields. The data field of the node is initialized to the
 * pointer @em data.
 */

 PilListNode *newPilListNode(void *data)
{

    return lnode_create(data);

}


/**
 * @brief
 *   Destroy a list node.
 *
 * @param node  Pointer to an existing list node.
 *
 * @return Nothing.
 *
 * The function destroys an existing list node object. The user data object
 * pointed to by @em node is nopt deallocated by the destructor!
 */

 void deletePilListNode(PilListNode *node)
{

    if (node)
        lnode_destroy(node);

    return;

}


/**
 * @brief
 *   Get the data stored in a list node.
 *
 * @param node  Pointer to a list node.
 *
 * @return The pointer to the node's user data is returned.
 *
 * The function returns a pointer to the data object referenced by the
 * list node @em node.
 *
 * @see pilListNodePut()
 */

 void *pilListNodeGet(const PilListNode *node)
{

    return lnode_get(node);

}


/**
 * @brief
 *   Assign data to a list node.
 *
 * @param node  Pointer to a list node.
 *
 * @return Nothing.
 *
 * The function sets the list node's data pointer to the provided user
 * data object. 
 *
 * @see pilListNodeGet()
 */

 void pilListNodePut(PilListNode *node, const void *data)
{

    lnode_put(node, (void *)data);
    return;

}


/**
 * @brief
 *   Check if a node is a member of any list.
 *
 * @param node  Pointer to a list node.
 *
 * @return The function returns 1 if @em node is a list member, otherwise the
 *   0 is returned.
 *
 * The function checks if the node @em node is a member of some list
 * by checking whether the node has a predecessor or a successor.
 */

 int pilListNodeIsListMember(const PilListNode *node)
{

    return lnode_is_in_a_list((PilListNode *)node);

}


/**
 * @brief
 *   Get the first node in a list.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The function returns a pointer to the first node in @m list.
 *   If the list is empty a @c NULL pointer is returned.
 *
 * The function retrieves the first node in the list @em list.
 *
 * @see pilListEnd()
 */

 PilListNode *pilListBegin(const PilList *list)
{

    return list_first(list);

}


/**
 * @brief
 *   Get the last node in a list.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The function returns a pointer to the last node in @em list.
 *   If the list is empty a @c NULL pointer is returned.
 *
 * The function retrieves the last node in the list @em list.
 *
 * @see pilListBegin()
 */

 PilListNode *pilListEnd(const PilList *list)
{

    return list_last(list);

}


/**
 * @brief
 *   Get the next node in a list.
 *
 * @param list  Pointer to an existing list.
 * @param node  Pointer to a list node.
 *
 * @return The function returns a pointer to the next node in @em list.
 *   If there is no next node a @c NULL pointer is returned.
 *
 * The function returns a pointer to the successor of the node @em node
 * in the list @em list.
 *
 * @see pilListPrev()
 */

 PilListNode *pilListNext(const PilList *list, const PilListNode *node)
{

    return list_next(list, node);

}


/**
 * @brief
 *   Get the previous node in a list.
 *
 * @param list  Pointer to an existing list.
 * @param node  Pointer to a list node.
 *
 * @return The function returns a pointer to the previous node in a list.
 *   If there is no previous node a @c NULL pointer is returned.
 *
 * The function returns a pointer to the predecessor of the node
 * @em node in the list @em list.
 *
 * @see pilListNext()
 */

 PilListNode *pilListPrev(const PilList *list, const PilListNode *node)
{

    return list_prev(list, node);

}


/**
 * @brief
 *   Remove all nodes from a list.
 *
 * @param list  Pointer to an existing list.
 *
 * @return Nothing.
 *
 * The function deletes all existing list nodes from the list @em list.
 * Each list node is first removed from the list and then destroyed. The
 * data objects pointed to by the list nodes are left untouched, in
 * particular they are not deallocated. It is the resposibility of the
 * caller to ensure that other references to the data object are still
 * present. After calling this function the list is empty.
 *
 * @see pilListErase(), pilListExtract(), pilListPopFront(), pilListPopBack()
 */

 void pilListClear(PilList *list)
{

    if (list)
        list_destroy_nodes(list);

    return;

}


/**
 * @brief
 *   Check whether a list is empty.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The return value is non-zero if the list is empty, otherwise the
 *   returned value is 0.
 *
 * The function tests if the list object @em list contains any data.
 *
 * @see pilListIsFull()
 */

 int pilListIsEmpty(const PilList *list)
{

    return list_isempty((PilList *)list);

}


/**
 * @brief
 *   Create a new list object without any elements.
 *
 * @return Returns the pointer to the newly created list. If the list
 *   could not be created @c NULL is returned.
 *
 * The function allocates memory for a new list object and initializes it
 * to an empty list.
 */

 PilList *newPilList(void)
{

    return list_create(LISTCOUNT_T_MAX);

}


/**
 * @brief
 *   Destroy a list object.
 *
 * @param list  Pointer to an existing list object.
 *
 * @return Nothing.
 *
 * The function destroys the list object @em list. The list must be
 * empty.
 *
 * @see pilListClear(), pilListDestroy()
 */

 void deletePilList(PilList *list)
{

    if (list)
        list_destroy(list);

    return;

}


/**
 * @brief
 *   Destroy a list object and all its nodes.
 *
 * @param list        Pointer to an existing list object.
 * @param deallocate  Data deallocator.
 *
 * @return Nothing.
 *
 * The function destroys the list object @em list. The list must be
 * empty. The data objects associated to the list nodes are deallocated
 * using the deallocation function @em deallocate.
 *
 * @see pilListClear(), deletePilList()
 */

 void pilListDestroy(PilList *list, void (*deallocate)(void *))
{

    PilListNode *node;

    if (!list)
        return;

    assert(deallocate != NULL);

    node = pilListBegin(list);
    while (node) {
        PilListNode *n = node;
        void *data;

        node = pilListNext(list, node);

        pilListExtract(list, n);

        data = pilListNodeGet(n);
        if (data)
            deallocate(data);

        deletePilListNode(n);
    }

    assert(pilListIsEmpty(list));
    list_destroy(list);

    return;

}


/**
 * @brief
 *   Get the actual number of list nodes.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The function returns the number of list nodes.
 *
 * The function retrieves the number of list nodes currently stored in
 * the list @em list.
 */

 PilListSize pilListSize(const PilList *list)
{

    return list_count(list);

}


/**
 * @brief
 *   Get the actual number of list nodes.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The function returns the number of list nodes.
 *
 * The function retrieves the number of list nodes currently stored in
 * the list @em list.
 */

 PilListSize pilListMaxSize(const PilList *list)
{

    assert(list != NULL);
    return LISTCOUNT_T_MAX;

}


/**
 * @brief
 *   Check whether a list is full.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The return value is non-zero if the list is full, otherwise the
 *   returned value is 0.
 *
 * The function tests if the list object @em list is full, i.e.
 * if the number of list nodes equals the list's capacity.
 *
 * @see pilListIsEmpty()
 */

 int pilListIsFull(const PilList *list)
{

    return list_isfull((PilList *)list);

}


/**
 * @brief
 *   Check if a list contains a given node.
 *
 * @param list  Pointer to an existing list.
 * @param node  Pointer to a list node.
 *
 * @return The function returns 1 if the node is found in the list, otherwise
 *   the return value is 0.
 *
 * The function searches the list object @em list for the list node
 * @em node.
 *
 * @see pilListNodeIsMember()
 */

 int pilListContains(const PilList *list, const PilListNode *node)
{

    return list_contains((PilList *)list, (PilListNode *)node);

}


/**
 * @brief
 *   Verify the integrity of a list.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The function returns 1 if the list integrity could be verified,
 *   if the list is corrupted the return value is 0.
 *
 * The internal list structure is verified. The function is provided for
 * debugging purposes only and should be placed in assert statements.
 */

 int pilListVerify(const PilList *list)
{

    return list_verify((PilList *)list);

}


/**
 * @brief
 *   Insert a node into a list before a given node.
 *
 * @param list      Pointer to an existing list.
 * @param position  The list node specifying the location for the insertion.
 * @param node      The list node to be added.
 *
 * @return Pointer to the inserted node.
 *
 * The function inserts the node @em node immediately before the node
 * @em position into the list @em list.
 */

 PilListNode *pilListInsert(PilList *list, PilListNode *position,
                                  PilListNode *node)
{

    list_ins_before(list, node, position);
    return node;

}


/**
 * @brief
 *   Insert a node into a list after a given node.
 *
 * @param list      Pointer to an existing list.
 * @param position  The list node specifying the location for the insertion.
 * @param node      The list node to be added.
 *
 * @return Nothing.
 *
 * The function inserts the node @em node immediately after the node
 * @em position into the list @em list.
 */

 PilListNode *pilListInsertAfter(PilList *list, PilListNode *position,
                                       PilListNode *node)
{

    list_ins_after(list, node, position);
    return node;

}


/**
 * @brief
 *   Add a node at the beginning of a list.
 *
 * @param list  Pointer to an existing list.
 * @param node  The list node to be added.
 *
 * @return Nothing.
 *
 * The function prepends the node @em node to the list @em list. 
 */

 void pilListPushFront(PilList *list, PilListNode *node)
{

    list_prepend(list, node);
    return;

}


/**
 * @brief
 *   Add a node at the end of a list.
 *
 * @param list  Pointer to an existing list.
 * @param node  The list node to be added.
 *
 * @return Nothing.
 *
 * The function appends the node @em node to the list @em list.
 */

 void pilListPushBack(PilList *list, PilListNode *node)
{

    list_append(list, node);
    return;

}


/**
 * @brief
 *   Erase a list element.
 *
 * @param list        Pointer to an existing list.
 * @param node        List node that should be erased.
 * @param deallocate  Data deallocation function.
 *
 * @return Nothing.
 *
 * The function removes the list node @em node from the list @em list and
 * deallocates it. The data stored in the list node is deallocated using
 * the data deallocator @em deallocate.
 *
 * @see pilListClear(), pilListExtract(), pilListPopFront(), pilListPopBack()
 */

 void pilListErase(PilList *list, PilListNode *node,
                         void (*deallocate)(void *))
{

    assert(list != NULL);
    assert(deallocate != NULL);

    if (node) {
        void *data = pilListNodeGet(node);

        list_delete(list, node);
 
        if (data)
            deallocate(data);

        lnode_destroy(node);
    }

    return;

}


/**
 * @brief
 *   Extract a node from a list.
 *
 * @return The function returns a pointer to the removed node.
 *
 * @param list  Pointer to an existing list.
 * @param node  Pointer to the list node which should be removed.
 *
 * The function removes the list node @em node from the list
 * @em list. The deleted node is returned.
 *
 * @see pilListErase(), pilListClear(), pilListPopFront(), pilListPopBack()
 */

 PilListNode *pilListExtract(PilList *list, PilListNode *node)
{

    return list_delete(list, node);

}


/**
 * @brief
 *   Remove the first node from a list.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The function returns a pointer to the removed node.
 *
 * The function removes the first node in the list @em list from this list.
 *
 * @see pilListPopBack(), pilListExtract(), pilListErase(), pilListClear()
 */

 PilListNode *pilListPopFront(PilList *list)
{

    return list_del_first(list);

}


/**
 * @brief
 *   Remove the last node from a list.
 *
 * @param list  Pointer to an existing list.
 *
 * @return The function returns a pointer to the removed node.
 *
 * The function removes the last node in the list @em list from this list.
 *
 * @see pilListPopFront(), pilListExtract(), pilListErase(), pilListClear()
 */

 PilListNode *pilListPopBack(PilList *list)
{

    return list_del_last(list);

}


/**
 * @brief
 *   Transfer a range of a list nodes.
 *
 * @param dst    Pointer to the destination list.
 * @param src    Pointer to the source list.
 * @param begin  First node of the slice.
 * @param end    Last node of the slice.
 *
 * @return Nothing.
 *
 * The function extracts from the source list @em src the sequence
 * of nodes starting with node @em start and ending with node
 * @em end. The extracted part of the source list is then appended to
 * the destination list @em dst.
 */

 void pilListTransfer(PilList *dst, PilList *src, PilListNode *begin,
                            PilListNode *end)
{

    list_extract(dst, src, begin, end);
    return;

}


/**
 * @brief
 *   Merge a list with another list.
 *
 * @param dst      Pointer to the destination list.
 * @param src      Pointer to the source list.
 * @param compare  Pointer to the comparison function.
 *
 * @return Nothing.
 *
 * The function merges the source list @em src into the destination
 * list @em dst. The comparison function @em compare is used to
 * compare data objects referenced by the list nodes. The comparison 
 * function takes two arguments of type @c const @c void @c * pointing to
 * the data objects being compared and returns an integer value less than,
 * equal or greater than zero if its first argument is less than, equal
 * or greater than its second argument. The first argument of the comparator
 * refers to data objects stored in the destination list, while its second
 * argument refers to data stored in the source list.
 *
 * @see pilListSplice()
 */

 void pilListMerge(PilList *dst, PilList *src, 
			 int (*compare)(const void *, const void *))
{

    list_merge(dst, src, compare);
    return;

}


/**
 * @brief
 *  Move a slice of a list to another list.
 *
 * @param dst    Pointer to the destination list.
 * @param pos    Pointer to the destination node.
 * @param src    Pointer to the source list.
 * @param begin  First node in the slice.
 * @param end    Last node of the slice.
 *
 * @return Nothing.
 *
 * The function extract the sequence of list nodes starting with the node
 * @em begin and ending with node @em end from the source list
 * @em src and inserts the extracted sequence into the destination
 * list @em dst before the node @em pos.
 *
 * @see pilListTransfer(), pilListMerge()
 */

 void pilListSplice(PilList *dst, PilListNode *pos, PilList *src,
			  PilListNode *begin, PilListNode *end)
{

    PilList *tmp = newPilList();

    pilListTransfer(tmp, dst, pos, pilListEnd(dst));
    pilListTransfer(dst, src, begin, end);
    pilListTransfer(dst, tmp, pilListBegin(tmp), pilListEnd(tmp));

    assert(pilListIsEmpty(tmp) == 1);
    deletePilList(tmp);

    return;

}


/**
 * @brief
 *   Sort a list.
 *
 * @param list     Pointer to an existing list.
 * @param compare  Comparison function used for sorting the list.
 *
 * @return Nothing.
 *
 * The function sorts the list @em list according to the criterium
 * imposed by the comparison function @em compare. The comparison
 * function @em compare is used to compare data objects referenced by
 * the list nodes. The comparison function takes two arguments of type
 * @c const @c void @c * pointing to the data objects being compared and
 * returns an integer value less than, equal or greater than zero if its
 * first argument is less than, equal or greater than its second argument.
 *
 * @see pilListIsSorted()
 */

 void pilListSort(PilList *list,
                        int (*compare)(const void *, const void *))
{

    list_sort(list, compare);
    return;

}


/**
 * @brief
 *   Test if a list is sorted.
 *
 * @param list     Pointer to an existing list.
 * @param compare  Comparison function used for checking the sorting order.
 *
 * @return The function returns 1 if the list is sorted, otherwise 0 is
 *   returned.
 *
 * The function steps through the list @em list and checks if the
 * order of the list nodes is the one defined by the comparison function
 * @em compare.
 *
 * @see pilListSort()
 */

 int pilListIsSorted(PilList *list,
			   int (*compare)(const void *, const void *))
{

    return list_is_sorted(list, compare);

}


/**
 * @brief
 *   Lookup a list entry.
 *
 * @param list     Pointer to an existing list.
 * @param data     Data to be searched for in the list.
 * @param compare  Data comparison function.
 *
 * @return Pointer to the list node if the entry was found, otherwise
 *   @c NULL is returned.
 *
 * The function locates the data @em data in the list @em list and returns a
 * pointer to the list node pointing to @em data. If such a node does not
 * exist a @c NULL pointer is returned. The data is searched by comparing
 * the data objects stored in the list with @em data using the data
 * comparison function @em compare. The comparison function takes two
 * arguments of type @c const @c void @c *. The first argument refers to the
 * data object referenced by a list node and the data @em data is used as
 * second argument. It is assumed that, if @em data matches a data object
 * in the list, the comparison function returns 0.
 */

 PilListNode *pilListLookup(PilList *list, const void *key,
				  int (*compare)(const void *, const void *))
{

    return list_find(list, key, compare);

}
/**@}*/
