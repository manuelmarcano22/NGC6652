/* $Id: pilframe.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pilmemory.h"
#include "pilstrutils.h"
#include "pilutils.h"
#include "pilframe.h"


typedef struct _PIL_PRODUCT_ATTRIBUTES_ PilProductAttributes;
typedef struct _PIL_FRAME_FLAGS_ PilFrameFlags;

struct _PIL_PRODUCT_ATTRIBUTES_ {
  PilProductLevel level;
  PilProductType type;
};

struct _PIL_FRAME_FLAGS_ {
  unsigned int keep;
  unsigned int ignore;
};

/**
 * @defgroup pilFrame pilFrame
 *
 * The module @b pilFrame provides functions to create, maintain and
 * destroy a frame.
 */

/**@{*/

/*
 * Definition of the frame data type
 */

struct _PIL_FRAME_ {
  char *filename;
  char *category;
  PilFrameType type;
  PilFrameFormat format;
  PilFrameFlags flags;
  PilProductAttributes attributes;
};


/**
 * @brief
 *   Set the frame keep flag to a particular value.
 *
 * @param frame   Pointer to an existing frame.
 * @param keep    Flag value to be stored.
 *
 * @return Returns @c EXIT_SUCCESS if no error occurred, if the flag could
 *   not be set @c EXIT_FAILURE is returned.
 *
 * The function sets the keep status of an existing frame to the given
 * value.
 */

static int pilFrmChangeKeepFlag(PilFrame *frame, unsigned int keep)
{

  if (!frame)
    return EXIT_FAILURE;

  frame->flags.keep = keep;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Set the frame ignore flag to a particular value.
 *
 * @param frame   Pointer to an existing frame.
 * @param ignore  Flag value to be stored.
 *
 * @return Returns @c EXIT_SUCCESS if no error occurred, if the flag could
 *   not be set @c EXIT_FAILURE is returned.
 *
 * The function sets the ignore status of an existing frame to the given
 * value.
 */

static int pilFrmChangeIgnoreFlag(PilFrame *frame, unsigned int ignore)
{

  if (!frame)
    return EXIT_FAILURE;

  frame->flags.ignore = ignore;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Creates a new frame.
 *
 * The function allocates and initializes a new, empty frame.
 *
 * @return Pointer to the newly allocated frame.
 */

PilFrame *newPilFrameEmpty(void)
{

  PilFrame *new = (PilFrame *)pil_malloc(sizeof(PilFrame));

  if (new) {
    new->filename = NULL;
    new->category = NULL;
    new->type = PIL_FRAME_TYPE_UNDEF;
    new->format = PIL_FRAME_FORMAT_UNDEF;
    new->attributes.level = PIL_PRODUCT_LEVEL_UNDEF;
    new->attributes.type = PIL_PRODUCT_TYPE_UNDEF;
    new->flags.keep = 0;
    new->flags.ignore = 0;
  }

  return new;

}


/**
 * @brief
 *   Creates a new frame.
 *
 * @param name Frame name.
 * @param category Frame category.
 *
 * @return Pointer to the newly allocated frame.
 *
 * The function allocates and initializes a new frame with a given
 * name and category.
 */

PilFrame *newPilFrame(const char *name, const char *category)
{

  PilFrame *new = newPilFrameEmpty();

  if (new) {
    if ((new->filename = pil_strdup(name)) == NULL) {
      deletePilFrame(new);
      return NULL;
    }

    if ((new->category = pil_strdup(category)) == NULL) {
      pil_free(new->filename);
      deletePilFrame(new);
      return NULL;
    }
  }

  return new;

}


/**
 * @brief
 *   Destroys an existing frame.
 *
 * @param frame Pointer to the frame that shall be destroyed.
 *
 * @return Nothing.
 *
 * The function removes and deallocates the frame.
 */

void deletePilFrame(PilFrame *frame)
{
  if (frame) {
    if (frame->filename)
      pil_free(frame->filename);
    if (frame->category)
      pil_free(frame->category);

    pil_free(frame);
  }

  return;

}


/**
 * @brief
 *   Get the frame filename.
 *
 * @param frame Pointer to a frame.
 *
 * @return Pointer to filename string.
 *
 * The function returns a pointer to the filename member of the frame
 * structure. The filename should not be modified directly through this
 * pointer.
 */

const char *pilFrmGetName(const PilFrame *frame)
{

  return frame->filename;

}


/**
 * @brief
 *   Set the frame filename.
 *
 * @param frame Pointer to a frame.
 * @param filename New frame filename.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the filename of a given frame.
 */

int pilFrmSetName(PilFrame *frame, const char *filename)
{

  size_t length;

  if (!frame || !filename)
    return EXIT_FAILURE;

  if ((length = strlen(filename)) > PIL_PATHNAME_MAX)
    return EXIT_FAILURE;
  else {
    if (length++ > strlen(frame->filename))
      frame->filename = (char *)pil_realloc(frame->filename, length);
    
    if (!frame->filename)
      return EXIT_FAILURE;
    else
      memcpy(frame->filename, filename, length);
  }

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the frame category.
 *
 * @param frame Pointer to a frame.
 *
 * @return Frame category.
 *
 * The function returns the category of a given frame.
 */

const char *pilFrmGetCategory(const PilFrame *frame)
{

  return frame->category;

}


/**
 * @brief
 *   Set the frame category.
 *
 * @param frame     Pointer to a frame.
 * @param category  New frame category.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the category of a given frame.
 */

int pilFrmSetCategory(PilFrame *frame, const char *category)
{

  char *s = pil_strdup(category);

  if (!frame || !category || !s)
    return EXIT_FAILURE;

  pil_free(frame->category);
  frame->category = s;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the frame type.
 *
 * @param frame Pointer to a frame.
 *
 * @return Frame type.
 *
 * The function returns the type of a given frame.
 */

PilFrameType pilFrmGetType(const PilFrame *frame)
{

  return frame->type;

}


/**
 * @brief
 *   Set the frame type.
 *
 * @param frame Pointer to a frame.
 * @param type New frame type.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the type of a given frame.
 */

int pilFrmSetType(PilFrame *frame, PilFrameType type)
{

  if (!frame)
    return EXIT_FAILURE;

  frame->type = type;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the data format of a frame.
 *
 * @param frame Pointer to a frame.
 *
 * @return The function returns the data format of @em frame.
 *
 * The function returns the frame's data format.
 */

PilFrameFormat pilFrmGetFormat(const PilFrame *frame)
{

  return frame->format;

}


/**
 * @brief
 *   Set the data format of a frame.
 *
 * @param frame Pointer to a frame.
 * @param format New frame data format.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the data format of a given frame.
 */

int pilFrmSetFormat(PilFrame *frame, PilFrameFormat format)
{

  if (!frame)
    return EXIT_FAILURE;

  frame->format = format;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the product level of a frame.
 *
 * @param frame  Pointer to a frame.
 *
 * @return The function returns the frame product level.
 *
 * The function returns the product level of a given frame. If @em frame
 * is not a product, i.e. the type is different from
 * @c PIL_FRAME_TYPE_PRODUCT, the function returns @c PIL_PRODUCT_LEVEL_UNDEF.
 */

PilProductLevel pilFrmGetProductLevel(const PilFrame *frame)
{

  if (frame->type != PIL_FRAME_TYPE_PRODUCT)
    return PIL_PRODUCT_LEVEL_UNDEF;

  return frame->attributes.level;

}


/**
 * @brief
 *   Set the product level of a frame.
 *
 * @param frame  Pointer to a frame.
 * @param level  Frame product level.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the product level of the frame @em frame. if the
 * input frame @em frame is not of type @c PIL_FRAME_TYPE_PRODUCT the
 * function fails.
 */

int pilFrmSetProductLevel(PilFrame *frame, PilProductLevel level)
{

  if (!frame || frame->type != PIL_FRAME_TYPE_PRODUCT)
    return EXIT_FAILURE;

  frame->attributes.level = level;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the product type of a frame.
 *
 * @param frame  Pointer to a frame.
 *
 * @return The function returns the frame's product type.
 *
 * The function returns the product type of a given frame. If @em frame
 * is not a product, i.e. the type is different from
 * @c PIL_FRAME_TYPE_PRODUCT, the function returns @c PIL_PRODUCT_TYPE_UNDEF.
 */

PilProductType pilFrmGetProductType(const PilFrame *frame)
{

  if (frame->type != PIL_FRAME_TYPE_PRODUCT)
    return PIL_PRODUCT_TYPE_UNDEF;

  return frame->attributes.type;

}


/**
 * @brief
 *   Set the product type of a frame.
 *
 * @param frame  Pointer to a frame.
 * @param typel  Frame product type.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the product type of the frame @em frame. If the
 * input frame @em frame is not of type @c PIL_FRAME_TYPE_PRODUCT the
 * function fails.
 */

int pilFrmSetProductType(PilFrame *frame, PilProductType type)
{

  if (!frame || frame->type != PIL_FRAME_TYPE_PRODUCT)
    return EXIT_FAILURE;

  frame->attributes.type = type;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the value of the keep flag from a frame.
 *
 * @param frame  Pointer to a frame.
 *
 * @return The function returns 1 if the flag is set, otherwise the function
 *   returns 0.
 *
 * The function returns the current setting of the keep flag of @em frame.
 */

unsigned int pilFrmGetKeepFlag(const PilFrame *frame)
{

  return frame->flags.keep == 0 ? 0 : 1;

}


/**
 * @brief
 *   Set the keep flag of a frame.
 *
 * @param frame  Pointer to a frame.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function sets the keep flag of the given frame @em frame.
 */

int pilFrmSetKeepFlag(PilFrame *frame)
{

  return pilFrmChangeKeepFlag(frame, 1);

}


/**
 * @brief
 *   Reset the keep flag of a frame.
 *
 * @param frame  Pointer to a frame.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * The function resets the keep flag of the given frame @em frame to 0. If
 * the flag is 0 the frame may be deleted when discarding temporary data
 * from the working area.
 */

int pilFrmResetKeepFlag(PilFrame *frame)
{

  return pilFrmChangeKeepFlag(frame, 0);

}


/**
 * @brief
 *   Get the the value of the ignore flag from a frame.
 *
 * @param frame  Pointer to an existing frame.
 *
 * @return Returns 0 if the ignore flag is not set, if it is set the return
 *   value is non-zero.
 *
 * The function retrieves the value of the ignore status of the given
 * pipeline frame. If the flag is non-zero for the given frame, the frame
 * is excluded from subsequent processing steps.
 */

unsigned int pilFrmGetIgnoreFlag(const PilFrame *frame)
{

  return frame->flags.ignore == 0 ? 0 : 1;

}


/**
 * @brief
 *   Set the ignore flag of a frame.
 *
 * @param frame  Pointer to an existing frame.
 *
 * @return Returns @c EXIT_SUCCESS if no error occurred, if the flag could
 *   not be set @c EXIT_FAILURE is returned.
 *
 * The function sets the ignore status of the given frame @em frame to
 * exclude it from subsequent processing steps.
 */

int pilFrmSetIgnoreFlag(PilFrame *frame)
{

  return pilFrmChangeIgnoreFlag(frame, 1);

}


/**
 * @brief
 *   Resets the ignore flag of a frame.
 *
 * @param frame  Pointer to an existing frame.
 *
 * @return Returns @c EXIT_SUCCESS if no error occurred, if the flag could
 *   not be reset @c EXIT_FAILURE is returned.
 *
 * The function resets the ignore status of the given frame so that the
 * frame will be treated by the subsequent processing steps.
 */

int pilFrmResetIgnoreFlag(PilFrame *frame)
{

  return pilFrmChangeIgnoreFlag(frame, 0);

}


/**
 * @brief
 *   Duplicate a frame
 *
 * @param frame  Pointer to an existing frame.
 *
 * @return The newly created copy of the input frame @em frame.
 *
 * A new frame is created and initialized with the values of the input frame
 * @em frame.
 */

PilFrame *pilFrmDuplicate(PilFrame *frame)
{

  PilFrame *copy = newPilFrame(frame->filename, frame->category);

  if (copy) {
    copy->type = frame->type;
    copy->format = frame->format;
    copy->attributes.level  = frame->attributes.level;
    copy->attributes.type  = frame->attributes.type;
    copy->flags.keep   = frame->flags.keep;
    copy->flags.ignore = frame->flags.ignore;
  }

  return copy;

}
/**@}*/
