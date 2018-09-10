/* $Id: pildictionary.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#include "pildictionary.h"


/**
 * @defgroup pilDictionary pilDictionary
 *
 * The module @b pilDictionary provides functions to create, destroy and
 * maintain and a dictionary.
 *
 * @note 
 *   The current implementation is based on the @b kazlib dictionary 
 *   implementation. Not all services defined there are available through
 *   this interface.
 */

/**@{*/

/**
 * @brief
 *   Create a new dictionary node.
 *
 * @param data  Pointer to the node data.
 *
 * @return The function returns the pointer to the newly created node if no
 *    error occurs, otherwise the return value is @c NULL.
 *
 * The function allocates the memory for a dictionary node and initializes
 * the node's fields. The data field of the node is initialized to the
 * pointer @em data.
 */

 PilDictNode *newPilDictNode(void *data)
{

  return dnode_create(data);

}


/**
 * @brief
 *   Destroys a dictionary node.
 *
 * @param node  Pointer to an existing dictionary node.
 *
 * @return Nothing.
 *
 * The function destroys an existing dictionary node object.
 */

 void deletePilDictNode(PilDictNode *node)
{

  dnode_destroy(node);
  return;

}


/**
 * @brief
 *   Get the keyword of a dictionary node.
 *
 * @param node  Pointer to a dictionary node.
 *
 * @return The pointer to the node's keyword is returned.
 *
 * The function returns a pointer to the keyword data of the node @em node.
 *
 * @see pilDictGetData(), pilDictPutData()
 */

 const void *pilDictGetKey(PilDictNode *node)
{

  return dnode_getkey(node);

}


/**
 * @brief
 *   Get the user data of a dictionary node.
 *
 * @param node  Pointer to a dictionary node.
 *
 * @return The pointer to the node's user data is returned.
 *
 * The function returns a pointer to the user supplied data of the node
 * @em node.
 *
 * @see pilDictGetData(), pilDictPutData()
 */

 void *pilDictGetData(PilDictNode *node)
{

  return dnode_get(node);

}


/**
 * @brief
 *   Set the user data of a dictionary node.
 *
 * @param node  Pointer to a dictionary node.
 * @param data  Pointer to the user data object.
 *
 * @return Nothing.
 *
 * The function associates the user data object @em data to the dictionary
 * node @em node.
 *
 * @see pilDictGetKey(), pilDictGetData()
 */

 void pilDictPutData(PilDictNode *node, void *data)
{

  dnode_put(node, data);
  return;

}


/**
 * @brief
 *   Create a new dictionary.
 *
 * @param capacity  Maximum number of allowed dictionary entries.
 * @param cmp       Comparator function used for keyword comparision.
 *
 * @return Returns the pointer to the newly created dictionary if no error
 *   occurs, otherwise the return value is @c NULL.
 *
 * The function creates a new dictionary with the given maximum number
 * of possible dictionary nodes and installs the provided comparator
 * function in the dictionary.
 */

 PilDictionary *newPilDictionary(PilDictCapacity capacity,
				       PilDictComparator cmp)
{

  return dict_create(capacity, cmp);

}


/**
 * @brief
 *   Destroys a dictionary object.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return Nothing.
 *
 * The function destroys an existing dictionary object. The dictionary
 * must be empty.
 *
 * @see pilDictClear()
 */

 void deletePilDictionary(PilDictionary *dict)
{

  dict_destroy(dict);
  return;

}


/**
 * @brief
 *   Configure a dictionary for duplicate entries.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return Nothing.
 *
 * The function configures the dictionary @em dict so that multiple entries
 * may have the same keyword.
 */

 void pilDictAllowDuplicates(PilDictionary *dict)
{

  dict_allow_dupes(dict);
  return;

}


/**
 * @brief
 *   Installs the dictionary allocator and deallocator functions.
 *
 * @param dict         Pointer to an existing dictionary object.
 * @param allocator    Allocator function.
 * @param deallocator  Deallocator function.
 * @param context      Dictionary context.
 *
 * @return Nothing.
 *
 * The function replaces the default allocator and deallocator functions
 * of the dictionary @em dict used for the creation and destruction of
 * dictionary nodes by two user supplied functions @em allocator and
 * @em deallocator of the type @c pilDictAllocator and @c pilDictDeallocator
 * respectively. The allocator and deallocator can only be changed for empty
 * dictionaries and both functions have to be replaced simultaneously. A
 * dictionary context @em context may be provided.
 */

 void pilDictSetAllocator(PilDictionary *dict,
				PilDictAllocator allocator,
				PilDictDeallocator deallocator, void *context)
{

  dict_set_allocator(dict, allocator, deallocator, context);
  return;

}


/**
 * @brief
 *   Verify dictionary integrity.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return The function returns 1 if the dictionary could be verified, if the
 *   dictionary is corrupted the return value is 0.
 *
 * The structure of the dictionary @em dict is verified. The function checks
 * only for a limited set of possible corruptions. A return value of 1 does
 * not necessarily imply that the directory structure is not corrupt.
 *
 * @note
 *   The function is provided for debugging purposes only and should be
 *   placed in assert statements.
 */

 int pilDictVerify(PilDictionary *dict)
{

  return dict_verify(dict);

}


/**
 * @brief
 *   Check if a dictionary is empty.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return The return value is different from 0 if the dictionary is empty,
 *   if it is not empty 0 is returned.
 *
 * The function checks if the dictionary @em dict is empty.
 *
 * @see pilDictIsFull()
 */

 int pilDictIsEmpty(PilDictionary *dict)
{

  return dict_isempty(dict);

}


/**
 * @brief
 *   Check if a dictionary is full.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return The return value is different from 0 if the dictionary is full,
 *   if it is not full 0 is returned.
 *
 * The function checks if the dictionary @em dict is full, i.e. if the
 * number of dictionary nodes equals the dictionary's capacity.
 *
 * @see pilDictIsEmpty()
 */

 int pilDictIsFull(PilDictionary *dict)
{

  return dict_isfull(dict);

}


/**
 * @brief
 *   Report the current number of dictionary nodes.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return The function returns the number of dictionary nodes.
 *
 * The function reports the number of dictionary nodes in the dictionary
 * @em dict.
 */

 PilDictCapacity pilDictCapacity(PilDictionary *dict)
{

  return dict_count(dict);

}


/**
 * @brief
 *   Lookup a dictionary entry.
 *
 * @param dict     Pointer to an existing dictionary.
 * @param keyword  Keyword string to be searched in the dictionary.
 *
 * @return Pointer to the dictionary node if the entry was found, otherwise
 *   @c NULL is returned.
 *
 * The function locates a dictionary node in the dictionary @em dict
 * associated to the key @em keyword. If such a node is not found a
 * @c NULL pointer is returned.
 *
 * @see pilDictInsert(), pilDictErase(), pilDictRemove()
 */

 PilDictNode *pilDictLookup(PilDictionary *dict, const void *keyword)
{

  return dict_lookup(dict, keyword);

}


/**
 * @brief
 *   Check if a node is contained in a dictionary.
 *
 * @param dict  Pointer to an existing dictionary.
 * @param node  Dictionary node to be checked.
 *
 * @return The return value is 1, if the node is found in the dictionary,
 *   otherwise 0 is returned.
 *
 * The function checks is the dictionary node @em node is a member of the
 * dictionary @em dict.
 *
 * @see pilDictLookup()
 */

 int pilDictContains(PilDictionary *dict, PilDictNode *node)
{

  return dict_contains(dict, node);

}


/**
 * @brief
 *   Insert a new entry into a dictionary.
 *
 * @param dict     Pointer to an existing dictionary.
 * @param keyword  String to be used as node keyword.
 * @param data     Pointer to user data.
 *
 * @return The function returns 1 if the node insertion was successful,
 *   otherwise the return value is 0.
 *
 * The function allocates a new dictionary node using the dictionary's
 * allocator function and initializes the newly created node with
 * the reference to the user supplied data @em data. The node is then
 * inserted into the dictionary @em dict using the key @em keyword.
 *
 * @note
 *   Only the reference to the keyword data is stored in the dictionary.
 *
 * @see pilDictInsertNode(), pilDictLookup(), pilDictErase(), pilDictRemove()
 */

 int pilDictInsert(PilDictionary *dict, const void *keyword, void *data)
{

  return dict_alloc_insert(dict, keyword, data);

}


/**
 * @brief
 *   Insert an already existing dictionary node into a dictionary.
 *
 * @param dict     Pointer to an existing dictionary.
 * @param node     Pointer to an existing dictionary node.
 * @param keyword  String to be used as the node keyword.
 *
 * @return Nothing
 *
 * The function inserts an already created dictionary node @em node in the
 * dictionary @em dict using the dictionary key @em keyword. The data field
 * of the node is not changed.
 *
 * @note
 *   Only the reference to the keyword data is stored in the dictionary.
 *
 * @see pilDictInsert()
 */

 void pilDictInsertNode(PilDictionary *dict, PilDictNode *node,
			      const void *keyword)
{

  dict_insert(dict, node, keyword);
  return;

}


/**
 * @brief
 *   Remove a node from a dictionary.
 *
 * @param dict  Pointer to an existing dictionary.
 * @param node  Pointer to the dictionary node that should be removed.
 *
 * @return The removed dictionary node.
 *
 * The function removes the dictionary node @em node from the dictionary
 * @em dict. The removed node is returned.
 *
 * @see pilDictLookup(), pilDictInsert(), pilDictErase()
 */

 PilDictNode *pilDictRemove(PilDictionary *dict, PilDictNode *node)
{

  return dict_delete(dict, node);

}


/**
 * @brief
 *   Remove and delete a dictionary node.
 *
 * @param dict  Pointer to an existing dictionary.
 * @param node  Pointer to the dictionary node that should be removed.
 *
 * @return Nothing.
 *
 * The function removes the dictionary node @em node from the dictionary
 * @em dict. The node itself is destroyed using the dictionary's deallocator
 * function.
 *
 * @note
 *   Since the default deallocator does not know anything about the user
 *   supplied data referenced by the node, it cannot deallocate the 
 *   user data. This might result in a memory leak if no appropriate
 *   allocator and deallocator functions were installed!
 *
 * @see pilDictLookup(), pilDictInsert(), pilDictRemove()
 */

 void pilDictErase(PilDictionary *dict, PilDictNode *node)
{

  dict_delete_free(dict, node);
  return;

}


/**
 * @brief
 *   Deletes a dictionary and all its nodes.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return Nothing.
 *
 * The function deletes all existing dictionary nodes in the dictionary
 * @em dict. Each dictionary node is removed from the dictionary and deleted
 * using the dictionary's deallocater function.
 *
 * @note
 *   Since the default deallocator does not know anything about the user
 *   supplied data referenced by a node, it cannot deallocate the 
 *   user data. This might result in a memory leak if no appropriate
 *   allocator and deallocator functions were installed!
 *
 * @see pilDictDelete(), deletePilDictionary(), pilDictSetAllocator()
 */

 void pilDictClear(PilDictionary *dict)
{

  dict_free_nodes(dict);
  return;

}


/**
 * @brief
 *   Get the first node in a dictionary.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return The function returns a pointer to the first node in the
 *   dictionary, if the dictionary is empty a @c NULL pointer is returned.
 *
 * The function looks for the first dictionary node in the dictionary
 * @em dict, i.e. the node having the lowest dictionary keyword.
 *
 * @see pilDictEnd()
 */

 PilDictNode *pilDictBegin(PilDictionary *dict)
{

  return dict_first(dict);

}


/**
 * @brief
 *   Get the last node in a dictionary.
 *
 * @param dict  Pointer to an existing dictionary.
 *
 * @return The function returns a pointer to the last node in the
 *   dictionary, if the dictionary is empty a @c NULL pointer is returned.
 *
 * The function looks for the last dictionary node in the dictionary
 * @em dict, i.e. the node having the highest dictionary keyword.
 *
 * @see pilDictBegin()
 */

 PilDictNode *pilDictEnd(PilDictionary *dict)
{

  return dict_last(dict);

}


/**
 * @brief
 *   Get the next node in a dictionary.
 *
 * @param dict  Pointer to an existing dictionary.
 * @param node  Pointer to a dictionary node.
 *
 * @return The function returns a pointer to the next node in the
 *   dictionary, if there is no next node a @c NULL pointer is returned.
 *
 * The function returns a pointer to the sucessor of the node @em node in
 * the dictionary @em dict.
 *
 * @see pilDictPrev()
 */

 PilDictNode *pilDictNext(PilDictionary *dict, PilDictNode *node)
{

  return dict_next(dict, node);

}


/**
 * @brief
 *   Get the previous node in a dictionary.
 *
 * @param dict  Pointer to an existing dictionary.
 * @param node  Pointer to a dictionary node.
 *
 * @return The function returns a pointer to the previous node in the
 *   dictionary, if there is no previous node a @c NULL pointer is returned.
 *
 * The function returns a pointer to the predecessor of the node @em node
 * in the dictionary @em dict.
 *
 * @see pilDictNext()
 */

 PilDictNode *pilDictPrev(PilDictionary *dict, PilDictNode *node)
{

  return dict_prev(dict, node);

}
/**@}*/
