/* $Id: vmimagearray.c,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#include "vmimagearray.h"
#include "cpl.h"


/*
 * The image array object
 */

struct _VIMOS_IMAGE_ARRAY_ {
    int capacity;
    int size;
    VimosImage **data;
};


/**
 * @defgroup VimosImageArray VimosImageArray
 *
 *   The module @b VimosImageArray implements a simple, fixed size array
 *   container for images.
 */

/**@{*/

/**
 * @brief
 *   Create an empty image array of appropriate size.
 *
 * @return The function returns the pointer to the allocated array object
 *   if no error occurred, or @c NULL otherwise.
 *
 * @param size  Size of the image array.
 *
 * The function allocates the memory for an image array of @em size elements.
 * All array elements are initialized to to @c NULL.
 *
 * @author R. Palsa
 */

VimosImageArray *newImageArray(int size)
{

    VimosImageArray *self = 0L;


    assert(size > 0);
    
    if ((self = (VimosImageArray *)cpl_malloc(sizeof *self))) {
      self->data = (VimosImage **)cpl_calloc(size, sizeof(VimosImage *));
      if (!self->data) {
        deleteImageArray(self);
        return 0L;
      }

      self->capacity = size;
      self->size = 0;
    }

    return self;

}


/**
 * @brief
 *   Delete an image array object.
 *
 * @return Nothing.
 *
 * @param array  Image array object to be deallocated.
 *
 * The function deallocates the image array object @em array. The array object
 * must be empty.
 *
 * @note
 *   If @em array is @c NULL nothing is done.
 *
 * @see destroyImageArray
 *
 * @author R. Palsa
 */

void deleteImageArray(VimosImageArray *array)
{

  if (array) {
    assert(imageArrayIsEmpty(array) == VM_TRUE);

    if (array->data)
      cpl_free(array->data);

    cpl_free(array);
  }

  return;

}


/**
 * @brief
 *   Destroy an image array object.
 *
 * @return Nothing.
 *
 * @param array  Image array object to be destroyed.
 *
 * The function works as @b deleteImageArray(), but any image which is stored
 * in the array is deallocated prior to the deallocation of the array object
 * @em array.
 *
 * @note
 *   If @em array is @c NULL nothing is done.
 *
 * @see deleteImageArray
 *
 * @author R. Palsa
 */

void destroyImageArray(VimosImageArray *array)
{

  if (array) {
    if (imageArrayIsEmpty(array) == VM_FALSE) {
      register int i;
      register int sz = imageArraySize(array);

      for (i = 0; i < sz; i++)
        deleteImage(imageArrayRemove(array, i));
    }

    deleteImageArray(array);

  }

  return;

}


/**
 * @brief
 *   Check if an image array object is empty.
 *
 * @return The function returns @c VM_TRUE if the array is empty, otherwise
 *   @c VM_FALSE is returned.
 *
 * @param array  Image array object to be tested.
 *
 * The function checks if any image is stored in the array.
 *
 * @author R. Palsa
 */

int imageArrayIsEmpty(const VimosImageArray *array)
{

  assert(array != 0L);

  return imageArraySize(array) == 0 ? VM_TRUE : VM_FALSE;

}


/**
 * @brief
 *   Get the capacity of an image array.
 *
 * @return The function returns the image array's capacity.
 *
 * @param array  Image array object.
 *
 * The function reports the capacity of the array @em array, i.e. the
 * maximum number of images which could be stored in the array.
 *
 * @author R. Palsa
 */

int imageArrayCapacity(const VimosImageArray *array)
{

  assert(array != 0L);

  return array->capacity;

}


/**
 * @brief
 *   Get the number of images stored in an image array.
 *
 * @return The function returns the number of stored images.
 *
 * @param array  Image array object.
 *
 * The function retrieves the number of images currently stored in the image
 * array object @em array. The reported number is not the array objects
 * capacity, i.e. the maximum number of images which could be stuffed into
 * the array.
 *
 * @author R. Palsa
 */

int imageArraySize(const VimosImageArray *array)
{

  assert(array != 0L);

  return array->size;

}


/**
 * @brief
 *   The function retrieves an image from an image array object.
 *
 * @return The function returns a handle to a stored image. If no image
 *   is associated to an array element @c NULL is returned.
 *
 * @param array  Image array object.
 * @param index  Array element offset.
 *
 * The function returns a reference to the image stored in the array
 * @em array. The offset @em index is the usual C array offset, i.e.
 * the first image is accessed by passing 0 as @em index.
 *
 * @see imageArraySet
 *
 * @author R. Palsa
 */

const VimosImage *imageArrayGet(const VimosImageArray *array, int index)
{

  assert(array != 0L);
  assert(index >= 0 && index < imageArrayCapacity(array));

  return array->data[index];

}


/**
 * @brief
 *   Set an element of an image array.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred, or 
 *   @c EXIT_FAILURE otherwise.
 *
 * @param array  Image array object.
 * @param index  Array element offset.
 * @param image  Image object.
 *
 * The functions stores the image @em image into the element with offset
 * @em index of the image array object @em array. To be successfull
 * @em index must be equal or larger than 0, smaller than the array's
 * capacity and the array element must not already reference another
 * image. If a stored image should be replaced by another image first use
 * @b imageArrayRemove() to extract an image from the array object.
 *
 * @note
 *   One cannot use this function to reset an array element, i.e. one cannot
 *   pass @c NULL as image @em image.
 *
 * @see imageArrayGet, imageArrayRemove
 *
 * @author R. Palsa
 */

int imageArraySet(VimosImageArray *array, int index, const VimosImage *image)
{

  assert(array != 0);
  assert(image != 0);
  assert(index >= 0 && index < imageArrayCapacity(array));

  if (!array->data[index]) {
    array->data[index] = (VimosImage *)image;
    array->size++;
  }
  else
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Remove an image from an image array object.
 *
 * @return The function returns a handle to the removed image.
 *
 * @param array  Image array object.
 * @param index  Array element offset.
 *
 * The function removes the image stored at the array offset @em index and
 * returns a reference to the removed image. The array element is reset to
 * @c NULL. If no image was stored in the array element the function succeeds
 * and the returned reference is @c NULL.
 *
 * @author R. Palsa
 */

VimosImage *imageArrayRemove(VimosImageArray *array, int index)
{

  VimosImage *image;


  assert(array != 0L);
  assert(index >= 0 && index < imageArrayCapacity(array));

  if (!array->data[index])
    image = 0L;
  else {
    image = array->data[index];
    array->data[index] = 0L;
    array->size--;
  }

  return image;

}


/**
 * @brief
 *   Get reference to the array objects data.
 *
 * @return The function returns the reference to the array data block. This
 *   function always succeeds.
 *
 * @param array  Image array object.
 *
 * The function returns a handle of the the array's data block.
 *
 * @note
 *   This function should be used with care.
 *
 * @author R. Palsa
 */

VimosImage **imageArrayGetData(const VimosImageArray *array)
{

    return array->data;

}
/**@}*/
