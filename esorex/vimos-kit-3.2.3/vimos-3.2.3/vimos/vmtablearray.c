/* $Id: vmtablearray.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2004 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-03-25 11:43:04 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>

#include <pilmemory.h>

#include "vmtablearray.h"
#include "cpl.h"


/**
 * @defgroup vmtablearray Table Arrays
 *
 *   The module @b VimosTableArray implements a simple, fixed size array
 *   container for tables.
 */

/**@{*/

/*
 * The table array object
 */

struct _VIMOS_TABLE_ARRAY_ {
    int capacity;
    int size;
    VimosTable **data;
};


/**
 * @brief
 *   Create an empty table array of appropriate size.
 *
 * @return The function returns the pointer to the allocated array object
 *   if no error occurred, or @c NULL otherwise.
 *
 * @param size  Size of the table array.
 *
 * The function allocates the memory for an table array of @em size elements.
 * All array elements are initialized to to @c NULL.
 *
 * @author R. Palsa
 */

VimosTableArray *newTableArray(int size)
{

    VimosTableArray *self = 0L;


    assert(size > 0);
    
    if ((self = (VimosTableArray *)cpl_malloc(sizeof *self))) {
      self->data = (VimosTable **)cpl_calloc(size, sizeof(VimosTable *));
      if (!self->data) {
        deleteTableArray(self);
        return 0L;
      }

      self->capacity = size;
      self->size = 0;
    }

    return self;

}


/**
 * @brief
 *   Delete an table array object.
 *
 * @return Nothing.
 *
 * @param array  Table array object to be deallocated.
 *
 * The function deallocates the table array object @em array. The array object
 * must be empty.
 *
 * @note
 *   If @em array is @c NULL nothing is done.
 *
 * @see destroyTableArray
 *
 * @author R. Palsa
 */

void deleteTableArray(VimosTableArray *array)
{

  if (array) {
    assert(tblArrayIsEmpty(array) == VM_TRUE);

    if (array->data)
      cpl_free(array->data);

    cpl_free(array);
  }

  return;

}


/**
 * @brief
 *   Destroy an table array object.
 *
 * @return Nothing.
 *
 * @param array  Table array object to be destroyed.
 *
 * The function works as @b deleteTableArray(), but any table which is stored
 * in the array is deallocated prior to the deallocation of the array object
 * @em array.
 *
 * @note
 *   If @em array is @c NULL nothing is done.
 *
 * @see deleteTableArray
 *
 * @author R. Palsa
 */

void destroyTableArray(VimosTableArray *array)
{

  if (array) {
    if (tblArrayIsEmpty(array) == VM_FALSE) {
      register int i;
      register int sz = tblArraySize(array);

      for (i = 0; i < sz; i++)
        deleteTable(tblArrayRemove(array, i));
    }

    deleteTableArray(array);

  }

  return;

}


/**
 * @brief
 *   Check if an table array object is empty.
 *
 * @return The function returns @c VM_TRUE if the array is empty, otherwise
 *   @c VM_FALSE is returned.
 *
 * @param array  Table array object to be tested.
 *
 * The function checks if any table is stored in the array.
 *
 * @author R. Palsa
 */

int tblArrayIsEmpty(const VimosTableArray *array)
{

  assert(array != 0L);

  return tblArraySize(array) == 0 ? VM_TRUE : VM_FALSE;

}


/**
 * @brief
 *   Get the capacity of an table array.
 *
 * @return The function returns the table array's capacity.
 *
 * @param array  Table array object.
 *
 * The function reports the capacity of the array @em array, i.e. the
 * maximum number of tables which could be stored in the array.
 *
 * @author R. Palsa
 */

int tblArrayCapacity(const VimosTableArray *array)
{

  assert(array != 0L);

  return array->capacity;

}


/**
 * @brief
 *   Get the number of tables stored in an table array.
 *
 * @return The function returns the number of stored tables.
 *
 * @param array  Table array object.
 *
 * The function retrieves the number of tables currently stored in the table
 * array object @em array. The reported number is not the array objects
 * capacity, i.e. the maximum number of tables which could be stuffed into
 * the array.
 *
 * @author R. Palsa
 */

int tblArraySize(const VimosTableArray *array)
{

  assert(array != 0L);

  return array->size;

}


/**
 * @brief
 *   The function retrieves an table from an table array object.
 *
 * @return The function returns a handle to a stored table. If no table
 *   is associated to an array element @c NULL is returned.
 *
 * @param array  Table array object.
 * @param index  Array element offset.
 *
 * The function returns a reference to the table stored in the array
 * @em array. The offset @em index is the usual C array offset, i.e.
 * the first table is accessed by passing 0 as @em index.
 *
 * @see tblArraySet
 *
 * @author R. Palsa
 */

const VimosTable *tblArrayGet(const VimosTableArray *array, int index)
{

  assert(array != 0L);
  assert(index >= 0 && index < tblArrayCapacity(array));

  return array->data[index];

}


/**
 * @brief
 *   Set an element of an table array.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred, or 
 *   @c EXIT_FAILURE otherwise.
 *
 * @param array  Table array object.
 * @param index  Array element offset.
 * @param table  Table object.
 *
 * The functions stores the table @em table into the element with offset
 * @em index of the table array object @em array. To be successfull
 * @em index must be equal or larger than 0, smaller than the array's
 * capacity and the array element must not already reference another
 * table. If a stored table should be replaced by another table first use
 * @b tableArrayRemove() to extract an table from the array object.
 *
 * @note
 *   One cannot use this function to reset an array element, i.e. one cannot
 *   pass @c NULL as table @em table.
 *
 * @see tblArrayGet, tblArrayRemove
 *
 * @author R. Palsa
 */

int tblArraySet(VimosTableArray *array, int index, const VimosTable *table)
{

  assert(array != 0);
  assert(table != 0);
  assert(index >= 0 && index < tblArrayCapacity(array));

  if (!array->data[index]) {
    array->data[index] = (VimosTable *)table;
    array->size++;
  }
  else
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Remove an table from an table array object.
 *
 * @return The function returns a handle to the removed table.
 *
 * @param array  Table array object.
 * @param index  Array element offset.
 *
 * The function removes the table stored at the array offset @em index and
 * returns a reference to the removed table. The array element is reset to
 * @c NULL. If no table was stored in the array element the function succeeds
 * and the returned reference is @c NULL.
 *
 * @author R. Palsa
 */

VimosTable *tblArrayRemove(VimosTableArray *array, int index)
{

  VimosTable *table;


  assert(array != 0L);
  assert(index >= 0 && index < tblArrayCapacity(array));

  if (!array->data[index])
    table = 0L;
  else {
    table = array->data[index];
    array->data[index] = 0L;
    array->size--;
  }

  return table;

}


/**
 * @brief
 *   Get reference to the array objects data.
 *
 * @return The function returns the reference to the array data block. This
 *   function always succeeds.
 *
 * @param array  Table array object.
 *
 * The function returns a handle of the the array's data block.
 *
 * @note
 *   This function should be used with care.
 *
 * @author R. Palsa
 */

VimosTable **tblArrayGetData(const VimosTableArray *array)
{

    return array->data;

}
/**@}*/
