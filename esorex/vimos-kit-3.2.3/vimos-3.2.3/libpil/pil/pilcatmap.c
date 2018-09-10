/* $Id: pilcatmap.c,v 1.2 2013-08-07 18:22:16 cgarcia Exp $
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
#include "pilstrutils.h"
#include "pilutils.h"
#include "pilcatmap.h"


/**
 * @defgroup pilCatmap pilCatmap
 *
 * The module @b pilCatmap provides functions to create, maintain and
 * destroy a data category and a category mapping definition.
 */

/**@{*/

static int catgcmp(const void *category1, const void *category2)
{

  return strcmp(category1, category2);

}


static PilDictNode *
_pilCatmapAllocator(void *context)
{

    (void)context;  /* To avoid compiler warnings */

    return pil_malloc(sizeof(PilDictNode));

}


static void
_pilCatmapDeallocator(PilDictNode *node, void *context)
{

    (void)context;  /* To avoid compiler warnings */

    if (node != NULL) {
        char *key = (void *)pilDictGetKey(node);
        PilCategory *data = pilDictGetData(node);

        if (data != NULL) {
            deletePilCategory(data);
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
 *   Creates a new empty category.
 *
 * @return Pointer to the newly allocated category.
 *
 * The function allocates and initializes a new, empty category.
 */

PilCategory *
newPilCategoryEmpty(void)
{

  PilCategory *new = (PilCategory *)pil_malloc(sizeof(PilCategory));

  if (new) {
    new->name = NULL;
    new->value = NULL;
  }

  return new;

}


/**
 * @brief
 *   Creates a new category.
 *
 * @param name   Category name.
 * @param value  Category value.
 *
 * @return Pointer to the newly allocated category.
 *
 * The function allocates and initializes a new category with given
 * name and value.
 */

PilCategory *
newPilCategory(const char *name, const char *value)
{

  PilCategory *category = (PilCategory *)pil_malloc(sizeof(PilCategory));

  category->name = NULL;
  category->value = NULL;

  if (pilCatSetName(category, name) == EXIT_FAILURE)
    return NULL;

  if (pilCatSetValue(category, value) == EXIT_FAILURE)
    return NULL;

  return category;

}


/**
 * @brief
 *   Destroys an existing category.
 *
 * @param category  Pointer to the category that shall be destroyed.
 *
 * @return Nothing.
 *
 * The function removes and deallocates the category.
 */

void
deletePilCategory(PilCategory *category)
{

  if (category) {
    if (category->name)
      pil_free(category->name);

    if (category->value)
      pil_free(category->value);

    pil_free(category);
  }

  return;

}


/**
 * @brief
 *   Initialize the object where the mapping of a set of categories
 *   will be contained.
 *
 * @return Pointer to the new category mapping.
 *
 * The function defines and allocates a dictionary of capacity
 * @c PIL_DICT_CAPACITY_MAX, with a method to compare nodes.
 */

PilCatmap *
newPilCatmap(void)
{

    PilCatmap *self = newPilDictionary(PIL_DICT_CAPACITY_MAX, catgcmp);


    pilDictSetAllocator(self, _pilCatmapAllocator, _pilCatmapDeallocator,
                        NULL);
    return self;

}


/**
 * @brief
 *   Destroys an existing category mapping.
 *
 * @param catmap  Pointer to the category mapping that shall be destroyed.
 *
 * @return Nothing.
 *
 * The function removes and deallocates the categories contained
 * in the category mapping dictionary, and the dictionary itself.
 */

void
deletePilCatmap(PilCatmap *catmap)
{

  if (!pilDictIsEmpty(catmap))
    pilDictClear(catmap);

  deletePilDictionary(catmap);

  return;

}


/**
 * @brief
 *   Set a category name.
 *
 * @param category  Pointer to a category.
 * @param name      New category name.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the name of a given category.
 */

int
pilCatSetName(PilCategory *category, const char *name)
{

  if (category) {
    if (name) {
      if (category->name) pil_free(category->name);
      if ((category->name = pil_strdup(name))) {
        return EXIT_SUCCESS;
      }
    }
  }
  return EXIT_FAILURE;

}


/**
 * @brief
 *   Set the value a category is mapped to.
 *
 * @param category  Pointer to a category.
 * @param value     New category value.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the value a given category is mapped to.
 */

int
pilCatSetValue(PilCategory *category, const char *value)
{

  if (category) {
    if (value) {
      if (category->value) pil_free(category->value);
      if ((category->value = pil_strdup(value))) {
        return EXIT_SUCCESS;
      }
    }
  }
  return EXIT_FAILURE;

}


/**
 * @brief
 *   Get a category name.
 *
 * @param category  Pointer to a category.
 *
 * @return Category name.
 *
 * The function returns the name of a given category.
 */

const char *
pilCatGetName(PilCategory *category)
{

  if (category)
    return pil_strdup(category->name);

  return NULL;

}


/**
 * @brief
 *   Get a category value.
 *
 * @param category  Pointer to a category.
 *
 * @return Category value.
 *
 * The function returns the value a given category is mapped to.
 */

const char *
pilCatGetValue(PilCategory *category)
{

  if (category)
    return pil_strdup(category->value);

  return NULL;

}


/**
 * @brief
 *   Insert a category in the category mapping definition.
 *
 * @param catmap    Pointer to category mapping.
 * @param category  Pointer to category to be inserted in the category
 *                  mapping definition.
 *
 * @return 1 on success, 0 on failure.
 *
 * The function inserts the category in the appropriate node
 * of the category mapping dictionary.
 */

int
pilCatmapInsert(PilCatmap *catmap, PilCategory *category)
{

  char *name;

  if (category) {
    name = pil_strdup(category->name);
    if (pilDictInsert(catmap, name, category)) 
      return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;

}


/**
 * @brief
 *   Remove a category from the category mapping definition.
 *
 * @param catmap  Pointer to category mapping.
 * @param name    Name of category to be removed.
 *
 * @return Nothing
 *
 * The function removes a category from the category mapping dictionary
 * and destroys it.
 */

void
pilCatmapRemove(PilCatmap *catmap, const char *name)
{

  char *categoryname;

  PilDictNode *node = pilDictLookup(catmap, name);

  if (node) {
    categoryname = (char *)pilDictGetKey(node);
    if (categoryname)
      pil_free(categoryname);
    deletePilCategory((PilCategory *)pilDictGetData(node));
    pilDictErase(catmap, node);
  }

  return;

}


/**
 * @brief
 *   Search a category in the category mapping definition.
 *
 * @param catmap  Pointer to category mapping.
 * @param name    Name of searched category.
 *
 * @return Pointer to found category, or @c NULL.
 *
 * The function searches a category in the category mapping dictionary
 * and returns it if found.
 */

PilCategory *
pilCatmapLookup(PilCatmap *catmap, const char *name)
{

  PilDictNode *node = pilDictLookup(catmap, name);

  if (!node)
    return NULL;

  return (PilCategory *)pilDictGetData(node);

}


/**
 * @brief
 *   Given a category name, get the value it was mapped to.
 *
 * @param catmap  Pointer to a category mapping.
 * @param name    Name of category to be looked for in the category mapping.
 *
 * @return Category value.
 *
 * The function returns the value a given category is mapped to.
 */

const char *
pilCatmapGetValue(PilCatmap *catmap, const char *name)
{

  PilCategory *category = pilCatmapLookup(catmap, name);

  if (category)
    return category->value;

  return NULL;

}
/**@}*/
