/* $Id: pilkeymap.c,v 1.2 2013-08-07 18:22:16 cgarcia Exp $
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
 * $Author: cgarcia $
 * $Date: 2013-08-07 18:22:16 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "pilmemory.h"
#include "pilalias.h"
#include "pilstrutils.h"
#include "pilutils.h"
#include "pilkeymap.h"


/**
 * @defgroup pilKeymap pilKeymap
 *
 * The module @b pilKeymap provides functions to create and maintain an 
 * alias table for frame header keywords.
 */

/**@{*/

static int
catgcmp(const void *category1, const void *category2)
{

    return strcmp(category1, category2);

}


static PilDictNode *
_pilKeymapAllocator(void *context)
{

    (void)context;  /* To avoid compiler warnings */

    return pil_malloc(sizeof(PilDictNode));

}


static void
_pilKeymapDeallocator(PilDictNode *node, void *context)
{

    (void)context;  /* To avoid compiler warnings */

    if (node != NULL) {
        char *key = (void *)pilDictGetKey(node);
        PilAlias *data = pilDictGetData(node);

        if (data != NULL) {
            deletePilAlias(data);
        }

        if (key != NULL) {
            pil_free(key);
        }

        pil_free(node);

    }

    return;

}


/**
 * @brief
 *   Initialize the object where the mapping of a set of keywords
 *   will be contained.
 *
 * @return Pointer to the new keyword mapping.
 *
 * The function defines and allocates a dictionary of capacity 
 * @c PIL_DICTCAPACITY_MAX, with a method to compare nodes.
 */

PilKeymap *
newPilKeymap(void)
{

    PilKeymap *self = newPilDictionary(PIL_DICT_CAPACITY_MAX, catgcmp);


    pilDictSetAllocator(self, _pilKeymapAllocator, _pilKeymapDeallocator,
                        NULL);
    return self;

}


/**
 * @brief
 *   Destroys an existing keyword mapping.
 *
 * @param keymap  Pointer to the keyword mapping that shall be destroyed.
 *
 * @return Nothing.
 *
 * The function removes and deallocates the keywords contained
 * in the keyword mapping dictionary, and the dictionary itself.
 */

void
deletePilKeymap(PilKeymap *keymap)
{

    if (!pilDictIsEmpty(keymap))
        pilDictClear(keymap);

    deletePilDictionary(keymap);

    return;

}


/**
 * @brief
 *   Insert a keyword in the keyword mapping definition.
 *
 * @param keymap Pointer to keyword mapping.
 * @param key    Pointer to keyword to be inserted in the keyword mapping
 *               definition.
 *
 * @return @c EXIT_SUCCESS on success, @c EXIT_FAILURE on failure.
 *
 * The function inserts the keyword in the appropriate node
 * of the keyword mapping dictionary.
 */

int
pilKeymapInsert(PilKeymap *keymap, PilAlias *key)
{

    char *name;

    if (key) {
        name = pil_strdup(pilAliasGetName(key));

        if (pilDictInsert(keymap, name, key)) 
            return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;

}


/**
 * @brief
 *   Remove a keyword from the keyword mapping definition.
 *
 * @param keymap  Pointer to keyword mapping.
 * @param name    Name of keyword to be removed.
 *
 * @return Nothing
 *
 * The function removes a keyword from the keyword mapping dictionary
 * and destroys it.
 */

void
pilKeymapRemove(PilKeymap *keymap, const char *name)
{

    char *keyname;


    PilDictNode *node = pilDictLookup(keymap, name);

    if (node) {
        keyname = (char *)pilDictGetKey(node);

        if (keyname)
            pil_free(keyname);

        deletePilAlias((PilAlias *)pilDictGetData(node));
        pilDictErase(keymap, node);
    }

    return;

}


/**
 * @brief
 *   Search a keyword in the keyword mapping definition.
 *
 * @param keymap  Pointer to keyword mapping.
 * @param name    Name of searched keyword.
 *
 * @return Pointer to found keyword, or NULL.
 *
 * The function searches a keyword in the keyword mapping dictionary
 * and returns it if found.
 */

PilAlias *
pilKeymapLookup(PilKeymap *keymap, const char *name)
{

    PilDictNode *node = pilDictLookup(keymap, name);

    if (!node)
        return NULL;

    return (PilAlias *)pilDictGetData(node);

}


/**
 * @brief
 *   Given a keyword name, get the value it was mapped to.
 *
 * @param keymap  Pointer to a keyword mapping.
 * @param name    Name of keyword to be looked for in the keyword mapping.
 *
 * @return Keyword value.
 *
 * The function returns the value a given keyword is mapped to.
 */

const char *
pilKeymapGetValue(PilKeymap *keymap, const char *name)
{

    PilAlias *key = pilKeymapLookup(keymap, name);

    if (key)
        return pilAliasGetValue(key);

    return NULL;

}


/**
 * @brief
 *   Given a keyword name, get its format.
 *
 * @param keymap Pointer to a keyword mapping.
 * @param name Name of keyword to be looked for in the keyword mapping.
 *
 * @return Keyword format.
 *
 * The function returns the format of a given keyword found in the
 * specified keyword mapping.
 */

const char *
pilKeymapGetFormat(PilKeymap *keymap, const char *name)
{

    PilAlias *key = pilKeymapLookup(keymap, name);

    if (key)
        return pil_strdup(pilAliasGetFormat(key));

    return NULL;

}


/**
 * @brief
 *   Given a keyword name, get its comment.
 *
 * @param keymap  Pointer to a keyword mapping.
 * @param name    Name of keyword to be looked for in the keyword mapping.
 *
 * @return Keyword comment.
 *
 * The function returns the comment of a given keyword found in the
 * specified keyword mapping.
 */

const char *
pilKeymapGetComment(PilKeymap *keymap, const char *name)
{

    PilAlias *key = pilKeymapLookup(keymap, name);

    if (key)
        return pilAliasGetComment(key);

    return 0L;

}
/**@}*/
