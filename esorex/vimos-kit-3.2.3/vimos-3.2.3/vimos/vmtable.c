/* $Id: vmtable.c,v 1.3 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <regex.h>

#include <fitsio.h>
#include <fitsio2.h>

#include <pilmemory.h>
#include <pilmessages.h>
#include <cpl_msg.h>
#include <pilstrutils.h>
#include <pilutils.h>

#include "vmtypes.h"
#include "vmimage.h"
#include "vmtable.h"
#include "cpl.h"


/**
 * @name vimosTable
 *
 * @doc
 *   The module vimosTable provides functions related to tables
 *   and header keywords (descriptors)
 */

/*@{*/

/* Constructor of VimosDescValue */
VimosDescValue *newDescValue()
{
  VimosDescValue *newValue;
  char            modName[] = "newDescValue";
  
  /* allocate memory */
  newValue = (VimosDescValue *) cpl_malloc(sizeof(VimosDescValue));
  /* check if space was allocated */
  if (newValue == NULL) {
    cpl_msg_debug(modName,"Allocation Error");
    return(NULL);
  }
  
  /* initialize */
  newValue->i = 0;
  return(newValue);  
}


/* Desctructor of VimosDescValue */
void deleteDescValue(VimosDescValue *dValue)
{
  if (dValue) {
    cpl_free(dValue);
  }
  
}

/* Constructor of VimosDescriptor */
VimosDescriptor *newDescriptor()
{
  VimosDescriptor *newDesc;
  char             modName[] = "newDescriptor";
  
  /* allocate space */
  newDesc = (VimosDescriptor *) cpl_malloc(sizeof(VimosDescriptor));
  if (newDesc == NULL) {
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
  /* allocate space for name */
  newDesc->descName = (char *) cpl_malloc(VM_DESC_LENGTH*sizeof(char));
  if (newDesc->descName == NULL) {
    /* cleanup if something went wrong */
    cpl_free(newDesc);
    cpl_msg_debug(modName, "Allocation error");
    return(NULL);
  }
  /* set name to undefined */
  strcpy(newDesc->descName, "Undefined");
  /* set length */
  newDesc->len = 0;

  /* set other defaults */
  newDesc->descType = VM_VARTYPE_UNDEF;
  newDesc->prev = newDesc->next = NULL;
  
  /* allocate value */
  newDesc->descValue = newDescValue();
  if (newDesc->descValue == NULL) {
    /* cleanup */
    deleteDescriptor(newDesc);
    cpl_msg_debug(modName, "The function newDescValue has returned NULL");
    return(NULL);
  }

  /*allocate space for comment */
  newDesc->descComment = (char *) cpl_malloc(VM_DESC_LENGTH*sizeof(char));
  if (newDesc->descComment == NULL) {
    /*cleanup */
    deleteDescriptor(newDesc);
    cpl_msg_debug(modName, "Allocation error");
    return(NULL);
  }

  strcpy(newDesc->descComment, "");

  /* return address of new descriptor */
  return(newDesc);
  
}


/*
 * Desctructor of VimosDescriptor 
 */

void
deleteDescriptor(VimosDescriptor *vDesc)
{
 
  if (vDesc != NULL) {
      cpl_free(vDesc->descName);
      cpl_free(vDesc->descComment);

      switch (vDesc->descType) {
      case VM_STRING:
      case VM_INT_ARRAY:
      case VM_FLOAT_ARRAY:
      case VM_DOUBLE_ARRAY:
          if (vDesc->descValue->s) {
              cpl_free(vDesc->descValue->s);
          }
          break;

      default:
          break;
      }
      deleteDescValue(vDesc->descValue);

      cpl_free(vDesc);
  }

  return;

}

void
deleteAllDescriptors(VimosDescriptor *vDesc)
{

  if (vDesc != NULL) {

      while (vDesc) {
          VimosDescriptor *nextDesc = vDesc->next;

          deleteDescriptor(vDesc);
          vDesc = nextDesc;
      }

  }

  return;

}

int removeDescriptor(VimosDescriptor **vDesc, const char *dscName)
{
  VimosDescriptor *tDesc;
  int              numFound = 0;

 /*
  *  Deleting all descriptors with the same name
  */

  while ((tDesc = findDescriptor(*vDesc,dscName))) {
    numFound++;
    if (tDesc->prev) tDesc->prev->next = tDesc->next;
    if (tDesc->next) tDesc->next->prev = tDesc->prev;
    if (tDesc == *vDesc) {
      *vDesc = tDesc->next;
    }
    deleteDescriptor(tDesc);
  }
  return numFound;
}

VimosDescriptor *
newStringDescriptor(const char *name, const char *value, const char *comment)
{

    size_t sz = 0;

    VimosDescriptor *self = NULL;
  

    self = newDescriptor();

    if (self != NULL) {

        strcpy(self->descName, name);
        strcpy(self->descComment, comment);

        self->descType = VM_STRING;

        if (value != NULL) {
            sz = strlen(value);
        }
  
        self->descValue->s = (char *)cpl_malloc((sz + 1) * sizeof(char));
        if (self->descValue->s == NULL) {
            deleteDescriptor(self);
            return NULL;
        }

        if (sz > 0) {
            strcpy(self->descValue->s, value);
        }
        else {
            self->descValue->s[sz] = '\0';
        }
        self->len = sz + 1;

    }  

    return self;
 
}

VimosDescriptor *newBoolDescriptor(const char *name, VimosBool value,
                                   const char *comment)
{
  VimosDescriptor *tDesc;
  char             modName[] = "newBoolDescriptor";
  
  /* allocate new descriptor */
  tDesc = newDescriptor();
  if (tDesc == NULL) {
    cpl_msg_debug(modName, "The function newDescriptor has returned NULL");
    return(NULL);
  }
  
  /* copy name and comment */
  strcpy(tDesc->descName, name);
  strcpy(tDesc->descComment, comment);
  /* set type and value */
  tDesc->descType = VM_BOOL;
  tDesc->descValue->b = value;
  /* set length */
  tDesc->len = 1;
  
  /* return new descriptor */
  return(tDesc);
  
}

VimosDescriptor *newIntDescriptor(const char *name, int value, 
                                  const char *comment)
{
  VimosDescriptor *tDesc;
  char             modName[] = "newIntDescriptor";
  
  /* allocate new descriptor */
  tDesc = newDescriptor();
  if (tDesc == NULL) {
    cpl_msg_debug(modName, "The function newDescriptor has returned NULL");
    return(NULL);
  }
  
  /* copy name and comment */
  strcpy(tDesc->descName, name);
  strcpy(tDesc->descComment, comment);

  /* set type and value */
  tDesc->descType = VM_INT;
  tDesc->descValue->i = value;
  /* set length */
  tDesc->len = 1;
  
  /* return new descriptor */
  return(tDesc);
  
}

VimosDescriptor *newFloatDescriptor(const char *name, float value,
                                    const char *comment)
{
  VimosDescriptor *tDesc;
  char             modName[] = "newFloatDescriptor";
  
  /* allocate new descriptor */
  tDesc = newDescriptor();
  if (tDesc == NULL) {
    cpl_msg_debug(modName, "The function newDescriptor has returned NULL");
    return(NULL);
  }
  
  /* copy name and comment */
  strcpy(tDesc->descName, name);
  strcpy(tDesc->descComment, comment);

  /* set type and value */
  tDesc->descType = VM_FLOAT;
  tDesc->descValue->f = value;
  /* set length */
  tDesc->len = 1;
    
  /* return new descriptor */
  return(tDesc);
  
}

VimosDescriptor *newDoubleDescriptor(const char *name, double value,
                                     const char *comment)
{
  VimosDescriptor *tDesc;
  char             modName[] = "newDoubleDescriptor";
  
  /* allocate new descriptor */
  tDesc = newDescriptor();
  if (tDesc == NULL) {
    cpl_msg_debug(modName, "The function newDescriptor has returned NULL");
    return(NULL);
  }
  
  /* copy name and comment */
  strcpy(tDesc->descName, name);
  strcpy(tDesc->descComment, comment);

  /* set type and value */
  tDesc->descType = VM_DOUBLE;
  tDesc->descValue->d = value;
  /* set length */
  tDesc->len = 1;
    
  /* return new descriptor */
  return(tDesc);
  
}

VimosDescriptor *newIntArrayDescriptor(const char *name, int *value, 
                                       const char *comment, int len)
{
  int i;
  VimosDescriptor *tDesc;
  char             modName[] = "newIntArrayDescriptor";
  
  /* allocate new descriptor */
  tDesc = newDescriptor();
  if (tDesc == NULL) {
    cpl_msg_debug(modName, "The function newDescriptor has returned NULL");
    return(NULL);
  }
  
  /* copy name and comment */
  strcpy(tDesc->descName, name);
  strcpy(tDesc->descComment, comment);

  /* set type */
  tDesc->descType = VM_INT_ARRAY;
  /* allocate space for value and check */
  tDesc->descValue->iar = (int *) cpl_malloc(len*sizeof(int));
  if (tDesc->descValue->iar == NULL) {
    /* cleanup */
    deleteDescriptor(tDesc);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
  for (i = 0; i < len; i++) {
    tDesc->descValue->iar[i] = value[i];
  }
  
  /* and set length */
  tDesc->len = len;
  
  /* return new descriptor */
  return(tDesc);
  
}

VimosDescriptor *newFloatArrayDescriptor(const char *name, float *value, 
                                         const char *comment, int len)
{
  int i;
  VimosDescriptor *tDesc;
  char             modName[] = "newFloatArrayDescriptor";
  
  /* allocate new descriptor */
  tDesc = newDescriptor();
  if (tDesc == NULL) {
    cpl_msg_debug(modName, "The function newDescriptor has returned NULL");
    return(NULL);
  }
  
  /* copy name and comment */
  strcpy(tDesc->descName, name);
  strcpy(tDesc->descComment, comment);

  /* set type */
  tDesc->descType = VM_FLOAT_ARRAY;
  /* allocate space for value and check */
  tDesc->descValue->far = (float *) cpl_malloc(len*sizeof(float));
  if (tDesc->descValue->far == NULL) {
    /* cleanup */
    deleteDescriptor(tDesc);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
  for (i = 0; i < len; i++) {
    tDesc->descValue->far[i] = value[i];
  }
  
  /* and set length */
  tDesc->len = len;
  
  /* return new descriptor */
  return(tDesc);
  
}
VimosDescriptor *newDoubleArrayDescriptor(const char *name, double *value, 
                                          const char *comment, int len)
{
  int i;
  VimosDescriptor *tDesc;
  char             modName[] = "newDoubleArrayDescriptor";
  
  /* allocate new descriptor */
  tDesc = newDescriptor();
  if (tDesc == NULL) {
    cpl_msg_debug(modName, "The function newDescriptor has returned NULL");
    return(NULL);
  }
  
  /* copy name and comment */
  strcpy(tDesc->descName, name);
  strcpy(tDesc->descComment, comment);

  /* set type */
  tDesc->descType = VM_DOUBLE_ARRAY;
  /* allocate space for value and check */
  tDesc->descValue->dar = (double *) cpl_malloc(len*sizeof(double));
  if (tDesc->descValue->dar == NULL) {
    /* cleanup */
    deleteDescriptor(tDesc);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
  for (i = 0; i < len; i++) {
    tDesc->descValue->dar[i] = value[i];
  }
  
  /* and set length */
  tDesc->len = len;
  
  /* return new descriptor */
  return(tDesc);
  
}

VimosBool writeIntDescriptor(VimosDescriptor **desc, const char *name, 
                             int value, const char *comment)
{
  VimosDescriptor *tDesc;
  VimosDescriptor *tDesc2;
  char             modName[] = "writeIntDescriptor";
  
  tDesc = tDesc2 = findDescriptor(*desc, name);
  /* loop in case descriptor is written more than once, we want to replace
     the last one... */
  while (tDesc2 != NULL)  {
    tDesc = tDesc2;
    tDesc2 = findDescriptor(tDesc2->next, name);
  }

  if (tDesc == NULL) {
    tDesc = newIntDescriptor(name, value, comment);
    if (tDesc == NULL) {
      cpl_msg_debug(modName, "The function newIntDescriptor has returned NULL");
      return(VM_FALSE);
    }
    if (!addDesc2Desc(tDesc, desc)) {
      cpl_msg_debug(modName, "The function addDesc2Desc has returned an error");
      return(VM_FALSE);
    }
  }
  else {
    if (tDesc->len > 1) {
      cpl_free(tDesc->descValue->s);
    }
    tDesc->descType = VM_INT;
    tDesc->len = 1;
    tDesc->descValue->i = value;
    strcpy(tDesc->descComment, comment);

  }
  return(VM_TRUE);
}

VimosBool writeFloatDescriptor(VimosDescriptor **desc, const char *name, 
                               float value, const char *comment)
{
  VimosDescriptor *tDesc;
  VimosDescriptor *tDesc2;
  char             modName[] = "writeFloatDescriptor";
  
  tDesc = tDesc2 = findDescriptor(*desc, name);
  /* loop in case descriptor is written more than once, we want to replace
     the last one... */
  while (tDesc2 != NULL)  {
    tDesc = tDesc2;
    tDesc2 = findDescriptor(tDesc2->next, name);
  }

  if (tDesc == NULL) {
    tDesc = newFloatDescriptor(name, value, comment);
    if (tDesc == NULL) {
      cpl_msg_debug(modName, "The function newFloatDescriptor has returned NULL");
      return(VM_FALSE);
    }
    if (!addDesc2Desc(tDesc, desc)) {
      cpl_msg_debug(modName, "The function addDesc2Desc has returned an error");
      return(VM_FALSE);
    }
  }
  else {
    if (tDesc->len > 1) {
      cpl_free(tDesc->descValue->s);
    }
    tDesc->descType = VM_FLOAT;
    tDesc->len = 1;
    tDesc->descValue->f = value;
    strcpy(tDesc->descComment, comment);
  }
  return(VM_TRUE);
}

VimosBool writeDoubleDescriptor(VimosDescriptor **desc, const char *name, 
                                double value, const char *comment)
{  
  VimosDescriptor *tDesc;
  VimosDescriptor *tDesc2;
  char             modName[] = "writeDoubleDescriptor";
  
  tDesc = tDesc2 = findDescriptor(*desc, name);
  /* loop in case descriptor is written more than once, we want to replace
     the last one... */
  while (tDesc2 != NULL)  {
    tDesc = tDesc2;
    tDesc2 = findDescriptor(tDesc2->next, name);
  }

  if (tDesc == NULL) {
    tDesc = newDoubleDescriptor(name, value, comment);
    if (tDesc == NULL) {
      cpl_msg_debug(modName,"The function newDoubleDescriptor has returned NULL");
      return(VM_FALSE);
    }
    if (!addDesc2Desc(tDesc, desc)) {
      cpl_msg_debug(modName, "The function addDesc2Desc has returned an error");
      return(VM_FALSE);
    }
  }
  else {
    if (tDesc->len > 1) {
      cpl_free(tDesc->descValue->s);
    }
    tDesc->descType = VM_DOUBLE;
    tDesc->len = 1;
    tDesc->descValue->d = value;
    strcpy(tDesc->descComment, comment);
  }
  return(VM_TRUE);
}


VimosBool writeStringDescriptor(VimosDescriptor **desc, const char *name,
                                const char *value, const char *comment)
{
  VimosDescriptor *tDesc;
  VimosDescriptor *tDesc2;
  char             modName[] = "writeStringDescriptor";
  
  tDesc = tDesc2 = findDescriptor(*desc, name);
  /* loop in case descriptor is written more than once, we want to replace
     the last one... */
  while (tDesc2 != NULL)  {
    tDesc = tDesc2;
    tDesc2 = findDescriptor(tDesc2->next, name);
  }

  if (tDesc == NULL) {
    tDesc = newStringDescriptor(name, value, comment);
    if (tDesc == NULL) {
      cpl_msg_debug(modName,"The function newStringDescriptor has returned NULL");
      return(VM_FALSE);
    }
    if (!addDesc2Desc(tDesc, desc)) {
      cpl_msg_debug(modName, "The function addDesc2Desc has returned an error");
      return(VM_FALSE);
    }
  }
  else {
    if (tDesc->len > 1) {
      cpl_free(tDesc->descValue->s);
    }
    tDesc->descType = VM_STRING;
    tDesc->descValue->s = (char *) cpl_malloc((VM_DESC_LENGTH+1)*sizeof(char));
    if (tDesc->descValue->s == NULL) {
      deleteDescriptor(tDesc);
      cpl_msg_debug(modName, "Allocation Error");
      return(VM_FALSE);
    }
    
    strcpy(tDesc->descValue->s, value);
    tDesc->len = strlen(value);
    strcpy(tDesc->descComment, comment);
  }
  return(VM_TRUE);
}



/**
 * @memo
 *   Copy a descriptor from a descriptor list to another.
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param fromDesc Pointer to source descriptor list.
 * @param inName Name of descriptor in the source list.
 * @param toDesc Pointer to target descriptor list.
 * @param outName New name of descriptor in the target list.
 *
 * @doc
 *   Descriptors with the given name(s) are searched in the source 
 *   list and in the target list. If matches are found in both lists, 
 *   issue a warning if types don't match, and overwrite anyway.
 *   If the output name is set to NULL, the input descriptor is copied
 *   with the same name. If a descriptor having the same name of the
 *   output descriptor doesn't exist in the target list, try to insert
 *   it after a descriptor having the same name of the descriptor
 *   preceding the input descriptor in the source list, otherwise
 *   copy the descriptor at the target list end.
 *
 * @author C. Izzo
 */

VimosBool copyFromHeaderToHeader(VimosDescriptor *fromDesc, const char *inName,
                     VimosDescriptor **toDesc, const char *outName)
{
  const char       modName[] = "copyFromHeaderToHeader";
  VimosDescriptor *inDesc;
  VimosDescriptor *outDesc;
  VimosDescriptor *copyDesc;

  if ( toDesc   == NULL) return VM_FALSE;
  if (*toDesc   == NULL) return VM_FALSE;
  if ( inName   == NULL) return VM_FALSE;
  if ( outName  == NULL) outName = inName;

  if ((inDesc = findDescriptor(fromDesc, inName)) == NULL) return VM_FALSE;

 /*
  *  Renaming of input is necessary only when inName and outName
  *  are different - but we do it anyway. The copy is necessary to
  *  avoid modifying the source list.
  */
  copyDesc = copyOfDescriptor(inDesc);
  strcpy(copyDesc->descName, outName);

 
  if ((outDesc = findDescriptor(*toDesc, outName)) != NULL) {
   
   /*
    *  A descriptor having the output name is found in the target list.
    */
    if (inDesc->descType != outDesc->descType) {
      cpl_msg_warning(modName, 
                    "Type mismatch between descriptors %s and %s (ignored)",
                    inName, outName);
    }
    if (outDesc->prev != NULL) {
      insertDescriptor(toDesc, outDesc->prev->descName, copyDesc, 0);
    }
    else if (outDesc->next != NULL) {
      insertDescriptor(toDesc, outDesc->next->descName, copyDesc, 1);
    }
    else {
      removeDescriptor(toDesc, outName);
      *toDesc = copyDesc;
    }
  }
  else if (!strcmp(outName, inName) && (inDesc->prev != NULL)) {
   /*
    *  A descriptor with the name outName was NOT found in the output
    *  header. At this point we may just copy the input descriptor at
    *  the end of the output list: nevertheless, if inName is the same 
    *  as outName (as in this "if" block is), an attempt is made to 
    *  insert the new descriptor at the "same" position that the input 
    *  descriptor had in the source header, i.e., soon after a descriptor 
    *  having the same name of the descriptor preceding the input descriptor.
    */
    
    if (insertDescriptor(toDesc, inDesc->prev->descName, copyDesc, 0) == 
        VM_FALSE) {

        /*
         *  The insertion failed, try to append the descriptor at the end
         */
        return addDesc2Desc(copyDesc, toDesc);

    }
  }
  else {
   /*
    *  A descriptor with the name outName was NOT found in the output
    *  header, and inName is NOT equal to outName. The descriptor is 
    *  appended at the end of the target list.
    */
    return addDesc2Desc(copyDesc, toDesc);
  }
  return VM_TRUE;
}

/**
 * @memo
 *   Insert a descriptor at any position in a descriptor list.
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param desc Pointer to a descriptor list.
 * @param name Name of reference descriptor in the list.
 * @param newDesc Pointer to the descriptor to insert in the list.
 * @param before Flag for inserting the new descriptor just before
 *        or soon after the reference descriptor.
 *
 * @doc
 *   The descriptor is inserted in the list by just changing
 *   the involved ->next and ->prev pointers. If the descriptor 
 *   to insert is part of another list, it is duplicated before
 *   insertion in the new list. If the descriptor to insert is
 *   already part of the list itself, it is moved (not copied)
 *   to the insert position. If a descriptor with the same name
 *   already exists in the list, it is replaced by the new one.
 */

VimosBool
insertDescriptor(VimosDescriptor **desc, const char *name,
                 VimosDescriptor *newDesc, int before)
{

  VimosDescriptor        *tDesc;
  const char modName[]   = "insertDescriptor";

  if ( desc    == NULL) return VM_FALSE;
  if (*desc    == NULL) return VM_FALSE;
  if ( newDesc == NULL) return VM_FALSE;
  if ( name    == NULL) return VM_FALSE;

  if (newDesc->next != NULL || newDesc->prev != NULL) {
   /*
    *  The descriptor to insert in the new list is already part of
    *  a descriptor list.
    */
    tDesc = findDescriptor(*desc, newDesc->descName);
    if (tDesc == newDesc) {
     /*
      *  The descriptor to insert is from the list itself: just move it
      *  to the insert position. Here we just extract it from the list...
      */
      if (newDesc->prev) newDesc->prev->next = newDesc->next;
      if (newDesc->next) newDesc->next->prev = newDesc->prev;
      newDesc->prev = newDesc->next = NULL;
    }
    else {
     /*
      *  The descriptor to insert belongs to another descriptor list,
      *  make a of copy the descriptor, to avoid destroying its 
      *  original list after modification of its ->next and ->prev 
      *  pointers.
      */
      newDesc = copyOfDescriptor(newDesc);
    }
  }

 /*
  *  If a descriptor with the same name is present in the list, destroy it.
  */
  removeDescriptor(desc, newDesc->descName);

 /*
  *  Now insert the new descriptor to the right position.
  */

  if ((tDesc = findDescriptor(*desc, name)) != NULL) {
    if (before) {
      newDesc->next = tDesc;
      newDesc->prev = tDesc->prev;
      if (tDesc->prev) 
	tDesc->prev->next = newDesc;
      tDesc->prev = newDesc;
      if (!newDesc->prev) *desc = newDesc;
    }
    else {
      if (tDesc->next) {
        newDesc->prev = tDesc;
        newDesc->next = tDesc->next;
	if (tDesc->next) 
	  tDesc->next->prev = newDesc;
	tDesc->next = newDesc;
	
      }
      else {
        tDesc->next = newDesc;
        newDesc->prev = tDesc;
      }
    }
    return VM_TRUE;
  }
  else {
    cpl_msg_debug(modName, "Descriptor %s not found (appending).", name);
    return addDesc2Desc(newDesc, desc);
  }

  return VM_FALSE;
}

/**
 * @memo
 *   Create and insert an integer descriptor at any position of a 
 *   descriptor list.
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param desc    Pointer to a descriptor list.
 * @param name    New descriptor name.
 * @param value   Integer descriptor value.
 * @param comment Descriptor comment field.
 * @param refName Name of reference descriptor in the list.
 * @param before  Flag for inserting the new descriptor just before
 *                or soon after the reference descriptor.
 *
 * @doc
 *   The descriptor is created and inserted in the list using
 *   the insertDescriptor function.
 *
 * @see insertDescriptor
 */

VimosBool
insertIntDescriptor(VimosDescriptor **desc, const char *name, int value,
                    const char *comment, const char *refName, int before)
{

    VimosDescriptor *newDesc = newIntDescriptor(name, value, comment);

    if (newDesc) {
        return insertDescriptor(desc, refName, newDesc, before);
    }

    return VM_FALSE;

}

/**
 * @memo
 *   Create and insert a float descriptor at any position of a 
 *   descriptor list.
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param desc    Pointer to a descriptor list.
 * @param name    New descriptor name.
 * @param value   Float descriptor value.
 * @param comment Descriptor comment field.
 * @param refName Name of reference descriptor in the list.
 * @param before  Flag for inserting the new descriptor just before
 *                or soon after the reference descriptor.
 *
 * @doc
 *   The descriptor is created and inserted in the list using
 *   the insertDescriptor function.
 *
 * @see insertDescriptor
 *
 * @author C. Izzo
 */

VimosBool insertFloatDescriptor(VimosDescriptor **desc, const char *name, 
          float value, const char *comment, const char *refName, int before)
{  
  VimosDescriptor *newDesc   = newFloatDescriptor(name, value, comment);

  if (newDesc) return insertDescriptor(desc, refName, newDesc, before);
  return VM_FALSE;
}

/**
 * @memo
 *   Create and insert a double descriptor at any position of a 
 *   descriptor list.
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param desc    Pointer to a descriptor list.
 * @param name    New descriptor name.
 * @param value   Double descriptor value.
 * @param comment Descriptor comment field.
 * @param refName Name of reference descriptor in the list.
 * @param before  Flag for inserting the new descriptor just before
 *                or soon after the reference descriptor.
 *
 * @doc
 *   The descriptor is created and inserted in the list using
 *   the insertDescriptor function.
 *
 * @see insertDescriptor
 *
 * @author C. Izzo
 */

VimosBool insertDoubleDescriptor(VimosDescriptor **desc, const char *name, 
          double value, const char *comment, const char *refName, int before)
{  
  VimosDescriptor *newDesc   = newDoubleDescriptor(name, value, comment);

  if (newDesc) return insertDescriptor(desc, refName, newDesc, before);
  return VM_FALSE;
}

/**
 * @memo
 *   Create and insert a string descriptor at any position of a 
 *   descriptor list.
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param desc    Pointer to a descriptor list.
 * @param name    New descriptor name.
 * @param value   String descriptor value.
 * @param comment Descriptor comment field.
 * @param refName Name of reference descriptor in the list.
 * @param before  Flag for inserting the new descriptor just before
 *                or soon after the reference descriptor.
 *
 * @doc
 *   The descriptor is created and inserted in the list using
 *   the insertDescriptor function.
 *
 * @see insertDescriptor
 *
 * @author C. Izzo
 */

VimosBool insertStringDescriptor(VimosDescriptor **desc, const char *name, 
     const char *value, const char *comment, const char *refName, int before)
{  
  VimosDescriptor *newDesc   = newStringDescriptor(name, value, comment);

  if (newDesc) return insertDescriptor(desc, refName, newDesc, before);
  return VM_FALSE;
}

/**
 * @memo
 *   Delete all descriptors whose name includes a given string.
 *
 * @return Number of found and deleted descriptors.
 *
 * @param  desc      Descriptor list 
 * @param  substring Substring to look for
 *
 * @doc
 *   Descriptors containing the given string are searched along the 
 *   list. When a match is found the descriptor is removed from the 
 *   list. A wildcard "*" (with its usual meaning) may be used at 
 *   the beginning and/or at the end of the substring, if it is not
 *   used the substring is interpreted as a complete descriptor name. 
 *
 * @author P. Sartoretti and C. Izzo
 */

int deleteSetOfDescriptors(VimosDescriptor **desc, char *substring)
{
  char             modName[] = "deleteSetOfDescriptors";
  VimosDescriptor *currDesc  = *desc;
  VimosDescriptor *delDesc   =  NULL;
  int              nchar;
  int              numFound=0;
  int              flag;
  char            *pos, *str, *keep;
  char            *descName;

  keep  = str = cpl_strdup(substring);
  nchar = strlen(str);

  flag = (str[0] == '*') + 2*(str[nchar-1] == '*');

  if (flag == 0) return removeDescriptor(desc, str);
  if (flag  > 1) str[nchar-1] = '\0';
  if (flag != 2) str++;

  while (currDesc) {
    descName = currDesc->descName;
    if ((pos = strstr(descName, str)) != NULL ) {
      switch (flag) {
      case 1 :
        if ((pos+strlen(pos)) == (descName+nchar)) delDesc = currDesc;
        break;
      case 2 :
        if (pos != descName) break;
      default:
        delDesc = currDesc;
      }
    }
    currDesc=currDesc->next;
    if (delDesc) {
      if (delDesc->prev) delDesc->prev->next = delDesc->next;
      if (delDesc->next) delDesc->next->prev = delDesc->prev;
      if (delDesc == *desc) {
        *desc = delDesc->next;
      }
      cpl_msg_debug(modName, "Delete descriptor: %s\n", delDesc->descName);
      deleteDescriptor(delDesc);
      delDesc = NULL;
      numFound++;
    }
  }
  cpl_free(keep);
  return numFound;
}


/*
  Allocates a new VimosDescriptor and copies the value of the input descriptor 
  into the new descriptor. The links prev and next of the new copy are set to
  NULL. 
*/
VimosDescriptor *copyOfDescriptor(VimosDescriptor *desc)
{
  VimosDescriptor *newDesc;
  char             modName[] = "copyOfDescriptor";
  
  if (desc == NULL) {
    cpl_msg_debug(modName, "NULL input descriptor");
    return(NULL);
  }
  
  switch (desc->descType) {
  case VM_INT : {
    newDesc = newIntDescriptor(desc->descName, desc->descValue->i, 
                               desc->descComment);
    if (newDesc == NULL) {
      cpl_msg_debug(modName, "The function newIntDescriptor has returned NULL");
      return(NULL);
    }
    break;
  }
  case VM_INT_ARRAY : {
    newDesc = newIntArrayDescriptor(desc->descName, desc->descValue->iar,
                                    desc->descComment, desc->len);
    if (newDesc == NULL) {
      cpl_msg_debug(modName,
      "The function newIntArrayDescriptor has returned NULL");
      return(NULL);
    }
    break;
  }
  case VM_BOOL : {
    newDesc = newBoolDescriptor(desc->descName, desc->descValue->b,
                                desc->descComment);
    if (newDesc == NULL) {
      cpl_msg_debug(modName, "The function newBoolDescriptor has returned NULL");
      return(NULL);
    }
    break;
  }
  case VM_FLOAT : {
    newDesc = newFloatDescriptor(desc->descName, desc->descValue->f,
                                 desc->descComment);
    if (newDesc == NULL) {
      cpl_msg_debug(modName, "The function newFloatDescriptor has returned NULL");
      return(NULL);
    }
    break;
  }
  case VM_FLOAT_ARRAY : {
    newDesc = newFloatArrayDescriptor(desc->descName, desc->descValue->far,
                                      desc->descComment, desc->len);
    if (newDesc == NULL) {
      cpl_msg_debug(modName,
      "The function newFloatArrayDescriptor has returned NULL");
      return(NULL);
    }
    break;
  }
  case VM_DOUBLE : {
    newDesc = newDoubleDescriptor(desc->descName, desc->descValue->d,
                                  desc->descComment);
    if (newDesc == NULL) {
      cpl_msg_debug(modName, "The function newDoubleDescriptor has returned NULL");
      return(NULL);
    }
    break;
  }
  case VM_DOUBLE_ARRAY : {
    newDesc = newDoubleArrayDescriptor(desc->descName, desc->descValue->dar, 
                                       desc->descComment, desc->len);
    if (newDesc == NULL) {
      cpl_msg_debug(modName, 
      "The function newDoubleArrayDescriptor has returned NULL");
      return(NULL);
    }
    break;
  }
  case VM_STRING : {
    newDesc = newStringDescriptor(desc->descName, desc->descValue->s,
                                  desc->descComment);
    if (newDesc == NULL) {
      cpl_msg_debug(modName, 
      "The function newStringDescriptor has returned NULL");
      return(NULL);
    }
    break;
  }
  default: {
    cpl_msg_debug(modName, "Undefined type of value stored in the descriptor");
    return(NULL);
  }
  
  }  /* of switch */

  /* return address of new copy */
  return(newDesc);
}



/* Constructor of VimosColumnValue */
VimosColumnValue *newColumnValue()
{
  VimosColumnValue *newValue;
  char             modName[] = "newColumnValue";

  /* allocate memory and check */
  newValue = (VimosColumnValue *) cpl_malloc(sizeof(VimosColumnValue));
  if (newValue == NULL) {
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
  /* initialize */
  newValue->p = NULL;
  
  return(newValue);
  
}


/* Desctructor of VimosColumnValue */
void deleteColumnValue(VimosColumnValue *dValue)
{
  if (dValue) {
    cpl_free(dValue->p);
    cpl_free(dValue);
  }
  
}


/* Constructor of VimosColumn */
VimosColumn *newColumn()
{
  VimosColumn *newCol;
  char         modName[] = "newColumn";
  
  newCol = (VimosColumn *) cpl_malloc(sizeof(VimosColumn));
  if (newCol == NULL) {
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
  
  newCol->colName = (char *) cpl_malloc(VM_DESC_LENGTH*sizeof(char));
  if (newCol->colName == NULL) {
    /* cleanup */
    cpl_free(newCol);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
      
  strcpy(newCol->colName, "Undefined");

  newCol->colType = VM_VARTYPE_UNDEF;
  newCol->prev = newCol->next = NULL;
  newCol->len = 0;
 
  
  newCol->colValue = newColumnValue();
  if (newCol->colValue == NULL) {
    /* cleanup */
    deleteColumn(newCol);
    cpl_msg_debug(modName, "The function newColumnValue has returned NULL");
    return(NULL);
  }
  
  return(newCol);
  
}

VimosColumn *newIntColumn(int size, const char *name)
{
  VimosColumn *tCol;
  char         modName[] = "newIntColumn";
  
  tCol = newColumn();
  if (tCol == NULL) {
    cpl_msg_debug(modName, "The function newColumn has returned NULL");
    return(NULL);
  }

  strcpy(tCol->colName, name);
  tCol->len = size;
  tCol->colType = VM_INT;
  tCol->colValue->iArray = (int *) cpl_malloc(size*sizeof(int));
  if (tCol->colValue->iArray == NULL) {
    /* cleanup */
    deleteColumn(tCol);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }

  return(tCol);
}

VimosColumn *newFloatColumn(int size, const char *name)
{
  VimosColumn *tCol;
  char         modName[] = "newFloatColumn";
  
  tCol = newColumn();
  if (tCol == NULL) {
    cpl_msg_debug(modName, "The function newColumn has returned NULL");
    return(NULL);
  }

  strcpy(tCol->colName, name);
  tCol->len = size;
  tCol->colType = VM_FLOAT;
  tCol->colValue->fArray = (float *) cpl_malloc(size*sizeof(float));
  if (tCol->colValue->fArray == NULL) {
    /* cleanup */
    deleteColumn(tCol);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  return(tCol);
}

VimosColumn *newDoubleColumn(int size, const char *name)
{
  VimosColumn *tCol;
  char         modName[] = "newDoubleColumn";
  
  tCol = newColumn();
  if (tCol == NULL) {
    cpl_msg_debug(modName, "The function newColumn has returned NULL");
    return(NULL);
  }

  strcpy(tCol->colName, name);
  tCol->len = size;
  tCol->colType = VM_DOUBLE;
  tCol->colValue->dArray = (double *) cpl_malloc(size*sizeof(double));
  if (tCol->colValue->dArray == NULL) {
    /* cleanup */
    deleteColumn(tCol);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  return(tCol);
}

VimosColumn *newCharacterColumn(int size, const char *name)
{
  VimosColumn *tCol;
  char         modName[] = "newCharacterColumn";
  
  tCol = newColumn();
  if (tCol == NULL) {
    cpl_msg_debug(modName, "The function newColumn has returned NULL");
    return(NULL);
  }

  strcpy(tCol->colName, name);
  tCol->len = size;
  tCol->colType = VM_CHARACTER;
  tCol->colValue->cArray = (char *) cpl_malloc(size*sizeof(char));
  if (tCol->colValue->cArray == NULL) {
    /* cleanup */
    deleteColumn(tCol);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  return(tCol);
}


VimosColumn *newStringColumn(int size, const char *name)
{
  VimosColumn *tCol;
  char         modName[] = "newStringColumn";
  
  tCol = newColumn();
  if (tCol == NULL) {
    cpl_msg_debug(modName, "The function newColumn has returned NULL");
    return(NULL);
  }

  strcpy(tCol->colName, name);
  tCol->len = size;
  tCol->colType = VM_STRING;
  tCol->colValue->sArray = (char **) cpl_calloc(size, sizeof(char *));
  if (tCol->colValue->sArray == NULL) {
    /* cleanup */
    deleteColumn(tCol);
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
  return(tCol);
}


/* Desctructor of VimosColumn */
void deleteColumn(VimosColumn *vCol)
{
  int i;
  
  if (!vCol) {
    return;
  }

  if (vCol->colType == VM_STRING) {
    for (i = 0; i < vCol->len; i++) {
      cpl_free(vCol->colValue->sArray[i]);
    }
  }
  
  cpl_free(vCol->colName);
  deleteColumnValue(vCol->colValue);
  cpl_free(vCol);
  
}



/* 
   Allocates a new VimosTable. The pointers in the table are set to NULL, the
   field numColumns is set to 0, and the name of the table is set to
   CM_EMPTY_TABLE_STRING. 
*/
VimosTable *newTable()
{
  VimosTable *aNewTable;
  char        modName[] = "newTable";
  
  aNewTable = (VimosTable *) cpl_malloc(sizeof(VimosTable));
  if (aNewTable == NULL) {
    cpl_msg_debug(modName, "Allocation Error");
    return(NULL);
  }
  
  strcpy(aNewTable->name, VM_EMPTY_TABLE_STRING);
  
  aNewTable->descs = NULL;
  aNewTable->cols = NULL;
  aNewTable->numColumns = 0;
  aNewTable->fptr = NULL;

  return(aNewTable);
  
}


void deleteTable(VimosTable *table)
{
  VimosDescriptor *tmpDesc;
  VimosDescriptor *nxtDesc;
  VimosColumn     *tmpCol;
  VimosColumn     *nxtCol;
  

  if (!table) {
    return;
  }
  
  tmpDesc = table->descs;
  
  while (tmpDesc) {
    nxtDesc = tmpDesc->next;
    deleteDescriptor(tmpDesc);
    tmpDesc = nxtDesc;
  }

  tmpCol = table->cols;
  
  while (tmpCol) {
    nxtCol = tmpCol->next;
    deleteColumn(tmpCol);
    tmpCol = nxtCol;
  }
  
}


/*
  Removes all columns from a VimosTable
*/
void deleteAllColumns(VimosTable *table)
{
  VimosColumn *tmpCol;
  VimosColumn *nxtCol;
 
  if (table == NULL) {
    return;
  }
  
  tmpCol = table->cols;
  
  while (tmpCol) {
    nxtCol = tmpCol->next;
    deleteColumn(tmpCol);
    tmpCol = nxtCol;
  }
  table->cols = NULL;
  table->numColumns = 0;
  
}



VimosBool copyTableDescriptors(VimosTable *inTable, VimosTable *outTable)
{
  VimosDescriptor *inDesc;
  VimosDescriptor *outDesc;
  VimosDescriptor *tDesc;
  char             modName[] = "copyTableDescriptors";
  
  inDesc = inTable->descs;
  outDesc = outTable->descs;

  tDesc = outDesc;
  while (outDesc) {
    tDesc = outDesc;
    outDesc = outDesc->next;
  }

  outDesc = tDesc;
      
  while (inDesc) {
    tDesc = copyOfDescriptor(inDesc);
    if (tDesc == NULL) {
      cpl_msg_debug(modName, "The function copyOfDescriptor has returned NULL");
      return(VM_FALSE);
    }
    if (outDesc == NULL) {
      outTable->descs = tDesc;
      outDesc = tDesc;
    }
    else {
      outDesc->next =tDesc;
      tDesc->prev = outDesc;
    }
    inDesc = inDesc->next;
    outDesc = outDesc->next;
  }
  return(VM_TRUE);
  
}

VimosBool copyAllDescriptors(VimosDescriptor *inDesc,VimosDescriptor **outDesc)
{
  VimosDescriptor *tDesc1;
  VimosDescriptor *tDesc2;
  char             modName[] = "copyAllDescriptors";
  char descName[80];

  if (outDesc == NULL) {
    cpl_msg_debug(modName, "NULL input descriptor");
    return(VM_FALSE);
  }
  
  tDesc2 = *outDesc;
  
  while (inDesc) {
    strcpy(descName, inDesc->descName);
    if (strncmp(descName, "TTYPE", 5) && strncmp(descName, "TFORM", 5))
    {
      tDesc1 = copyOfDescriptor(inDesc);
      if (tDesc1 == NULL) {
        cpl_msg_debug(modName, "The function copyOfDescriptor has returned NULL");
        return(VM_FALSE);
      }
    
      if (*outDesc == NULL) {
        tDesc2 = tDesc1;
        *outDesc = tDesc1;
      }
      else {
        if (!addDesc2Desc(tDesc1, &tDesc2)) {
          cpl_msg_debug(modName, "The function addDesc2Desc has returned an error");
          return(VM_FALSE);
        }
        tDesc2 = tDesc1;
      }
    }
    inDesc = inDesc->next;
  }
  return(VM_TRUE);
}

/*
  Add a descriptor to the linked list of descriptors of the Table
*/
VimosBool addDesc2Table(VimosDescriptor *desc, VimosTable *table)
{
  char             modName[] = "addDesc2Table";

  if ( (table == NULL) || (desc == NULL) ) {
    cpl_msg_debug(modName, "Invalid input table or descriptor");
    return(VM_FALSE);
  }
  
  if (!addDesc2Desc(desc, &(table->descs))) {
    cpl_msg_debug(modName, "The function addDesc2Desc has returned an error");
    return(VM_FALSE);
  }
  return(VM_TRUE);
}


/*
  Add a linked list of descriptors to a linked list of descriptors
*/
VimosBool addDesc2Desc(VimosDescriptor *inDesc, VimosDescriptor **outDesc)
{
  VimosDescriptor *tDesc1;
  VimosDescriptor *tDesc2;
  char             modName[] = "addDesc2Desc";

  if ( (!outDesc) || (!inDesc) ) {
    cpl_msg_debug(modName, "Invalid input descriptor"); 
    return (VM_FALSE);
  }
  
  tDesc1 = *outDesc;
  tDesc2 = *outDesc;
  
  while (tDesc1) {
    tDesc2 = tDesc1;
    tDesc1 = tDesc1->next;
  }
  if (tDesc2 == NULL) {
    *outDesc = inDesc;
  }
  else {
    tDesc2->next = inDesc;
  }
  inDesc->prev = tDesc2;
  return(VM_TRUE);
}
/**
 * @memo
 *   Find a descriptor in a list.
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param desc      Descriptor list 
 * @param name      Descriptor name - or substring containing name
 *
 * @doc
 *   Descriptor with the given name is searched along the 
 *   list and the pointer to its value is returned. 
 *   A wildcard "*" (with its usual meaning) may be used to find the
 *   first descriptor of a given type (ex. HIERARCH), if it is not
 *   used the given string is interpreted as a complete descriptor name. 
 *
 * @return Pointer to the value of the descriptor found.
 */

VimosDescriptor *findDescriptor(VimosDescriptor *vDesc, const char *name)
{

    int flag;
    int nchar;

    VimosDescriptor *tmpDesc = vDesc;
    VimosDescriptor *result = NULL;


    nchar = strlen(name);
    flag = (name[0] == '*') + 2 * (name[nchar - 1] == '*');

    if (flag == 0) {
        while (!result && tmpDesc) {
            if (!strcmp(tmpDesc->descName, name)) {
                result = tmpDesc;
            }
            tmpDesc = tmpDesc->next;
        }
    }
    else {

        char *str = cpl_strdup(name);

        if (flag > 1) {
            str[nchar - 2] = '\0';
        }

        if (flag != 2) {
            str++;
        }

        while (!result && tmpDesc) {

            char *descName = tmpDesc->descName;
            char *pos = strstr(descName, str);


            if (pos != NULL) {

                switch (flag) {
                case 1:
                    if ((pos + strlen(pos)) == (descName + nchar)) {
                        result = tmpDesc;
                    }
                    break;

                case 2 :
                    if (pos == descName) {
                        result = tmpDesc;
                    }
                    break;

                default:
                    result = tmpDesc;
                    break;
                }
            }

            tmpDesc = tmpDesc->next;
        }

        cpl_free(str);

    }

    return result;

}

VimosDescriptor *findDescInTab(VimosTable *table, const char *name)
{
  char             modName[] = "findDescInTab";

  if (table == NULL) {
    cpl_msg_debug(modName, "Invalid input table");
    return(NULL);
  }
  
  return(findDescriptor(table->descs, name));
}



VimosBool readIntDescriptor(VimosDescriptor *vDesc, const char *name, 
                            int *val, char *comment)
{
  VimosDescriptor *desc;
  char             modName[] = "readIntDescriptor";
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  
  /* not found: */
  if (desc == NULL) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(VM_FALSE);
  }
  
  /* it's there, check type */
  if (desc->descType != VM_INT) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Descriptor %s is not integer", name);
    return(VM_FALSE);
  }
  
  /* descriptor of the right name and type exists, copy value */
  *val = desc->descValue->i;

  if (comment)
    strcpy(comment, desc->descComment);

  return(VM_TRUE);
    
}

VimosBool readIntDescFromTable(VimosTable *tab, const char *name, 
                               int *val, char *comment)
{
  char  modName[] = "readIntDescFromTable";

  if (tab == NULL) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  return(readIntDescriptor(tab->descs, name, val, comment));
  
}



VimosBool readBoolDescriptor(VimosDescriptor *vDesc, const char *name, 
                             VimosBool *val, char *comment)
{
  VimosDescriptor *desc;
  char  modName[] = "readBoolDescriptor";
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  
  /* not found: */
  if (desc == NULL) {
    *val = VM_FALSE;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(VM_FALSE);
  }
  
  /* it's there, check type */
  if (desc->descType != VM_BOOL) {
    *val = VM_FALSE;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Descriptor %s is not boolean", name);
    return(VM_FALSE);
  }
  
  /* descriptor of the right name and type exists, copy value */
  *val = desc->descValue->b;

  if (comment)
    strcpy(comment, desc->descComment);

  return(VM_TRUE);
    
}

VimosBool readBoolDescFromTable(VimosTable *tab, const char *name, 
                                VimosBool *val, char *comment)
{
  char  modName[] = "readBoolDescFromTable";

  if (tab == NULL) {
    *val = VM_FALSE;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  return(readBoolDescriptor(tab->descs, name, val, comment));
  
}



VimosBool readFloatDescriptor(VimosDescriptor *vDesc, const char *name, 
                              float *val, char *comment)
{
  VimosDescriptor *desc;
  char  modName[] = "readFloatDescriptor";
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  
  /* not found: */
  if (desc == NULL) {
    *val = 0.0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(VM_FALSE);
  }
  
  /* it's there, check type 
  In FITS keywords, no way to distinguish between Float and Double, so they
  are all put in Double by default. To avoid mess in the code, we handle
  both cases  */

  if (desc->descType == VM_FLOAT) {
    *val = desc->descValue->f;

    if (comment)
      strcpy(comment, "");

    return(VM_TRUE);
  }
  if (desc->descType == VM_DOUBLE) {
    *val = desc->descValue->d;

    if (comment)
      strcpy(comment, desc->descComment);

    return(VM_TRUE);
  }

  cpl_msg_debug(modName, "Descriptor %s is not float", name);

  return(VM_FALSE);
    
}

VimosBool readFloatDescFromTable(VimosTable *tab, const char *name, 
                                 float *val, char *comment)
{
  char  modName[] = "readFloatDescFromTable";

  if (tab == NULL) {
    *val = 0.0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  return(readFloatDescriptor(tab->descs, name, val, comment));
  
}



VimosBool readDoubleDescriptor(VimosDescriptor *vDesc, const char *name, 
                               double *val, char *comment)
{
  VimosDescriptor *desc;
  char             modName[] = "readDoubleDescriptor";
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  
  /* not found: */
  if (desc == NULL) {
    *val = 0.0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(VM_FALSE);
  }
  
  /* it's there, check type */
  if (desc->descType != VM_DOUBLE) {
    *val = 0.0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Descriptor %s is not double", name);
    return(VM_FALSE);
  }
  
  /* descriptor of the right name and type exists, copy value */
  *val = desc->descValue->d;

  if (comment)
    strcpy(comment, desc->descComment);

  return(VM_TRUE);
    
}

VimosBool readDoubleDescFromTable(VimosTable *tab, const char *name, 
                                  double *val, char *comment)
{
  char             modName[] = "readDoubleDescFromTable";

  if (tab == NULL) {
    *val = 0.0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  return(readDoubleDescriptor(tab->descs, name, val, comment));
  
}




VimosBool readStringDescriptor(VimosDescriptor *vDesc, const char *name, 
                               char *val, char *comment)
{
  VimosDescriptor *desc;
  char             modName[] = "readStringDescriptor";
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  
  /* not found: */
  if (desc == NULL) {
    strcpy(val, "");

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(VM_FALSE);
  }
  
  /* it's there, check type */
  if (desc->descType != VM_STRING) {
    strcpy(val, "");

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Descriptor %s is not a string", name);
    return(VM_FALSE);
  }
  
  /* descriptor of the right name and type exists, copy value */
  strcpy(val, desc->descValue->s);

  if (comment)
    strcpy(comment, desc->descComment);

  return(VM_TRUE);
    
}

VimosBool readStringDescFromTable(VimosTable *tab, const char *name, 
                                  char *val, char *comment)
{
  char             modName[] = "readStringDescFromTable";

  if (tab == NULL) {
    strcpy(val, "");

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  return(readStringDescriptor(tab->descs, name, val, comment));
  
}



int getDescriptorLength(VimosDescriptor *vDesc, const char *name)
{
  VimosDescriptor *desc;
  char             modName[] = "getDescriptorLength";
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  /* not found: */
  if (desc == NULL) {
    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(0);
  } else {
    return(desc->len);
  }
  
}


VimosBool readIntArrayDescriptor(VimosDescriptor *vDesc, const char *name, 
                                 int *val, char *comment, int len)
{
  int i;
  int maxLen;
  
  VimosDescriptor *desc;
  char             modName[] = "readIntArrayDescriptor";
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  
  /* not found: */
  if (desc == NULL) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(VM_FALSE);
  }
  
  /* it's there, check type */
  if (desc->descType != VM_INT_ARRAY) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Descriptor %s is not an array of integers", name);
    return(VM_FALSE);
  }
  
  if (len > desc->len) {
    maxLen = len;
  }
  else {
    maxLen = desc->len;
  }
  
  /* descriptor of the right name and type exists, copy value */
  for (i = 0; i < maxLen; i++) {
    val[i] = desc->descValue->iar[i];
  }

  if (comment)
    strcpy(comment, desc->descComment);
  
  return(VM_TRUE);
    
}

VimosBool readIntArrayDescFromTable(VimosTable *tab, const char *name, 
                                    int *val, char *comment, int len)
{
  char             modName[] = "readIntArrayDescFromTable";

  if (tab == NULL) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  return(readIntArrayDescriptor(tab->descs, name, val, comment, len));
  
}

VimosBool readFloatArrayDescriptor(VimosDescriptor *vDesc, const char *name, 
                                   float *val, char *comment, int len)
{
  int i;
  int maxLen;
  char             modName[] = "readFloatArrayDescriptor";
  
  VimosDescriptor *desc;
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  
  /* not found: */
  if (desc == NULL) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(VM_FALSE);
  }
  
  /* it's there, check type */
  if (desc->descType != VM_FLOAT_ARRAY) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Descriptor %s is not an array of floats", name);
    return(VM_FALSE);
  }
  
  if (len > desc->len) {
    maxLen = len;
  }
  else {
    maxLen = desc->len;
  }
  
  /* descriptor of the right name and type exists, copy value */
  for (i = 0; i < maxLen; i++) {
    val[i] = desc->descValue->far[i];
  }

  if (comment)
    strcpy(comment, desc->descComment);
  
  return(VM_TRUE);
    
}

VimosBool readFloatArrayDescFromTable(VimosTable *tab, const char *name, 
                                      float *val, char *comment, int len)
{
  char             modName[] = "readFloatArrayDescFromTable";

  if (tab == NULL) {
    *val = 0.0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  return(readFloatArrayDescriptor(tab->descs, name, val, comment, len));
  
}


VimosBool readDoubleArrayDescriptor(VimosDescriptor *vDesc, const char *name, 
                                    double *val, char *comment, int len)
{
  int i;
  int maxLen;
  char             modName[] = "readDoubleArrayDescriptor";
  
  VimosDescriptor *desc;
  
  /* do not have to validate input */

  /* look for descriptor of the right name */
  desc = findDescriptor(vDesc, name);
  
  /* not found: */
  if (desc == NULL) {
    *val = 0.0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Cannot find descriptor %s", name);
    return(VM_FALSE);
  }
  
  /* it's there, check type */
  if (desc->descType != VM_DOUBLE_ARRAY) {
    *val = 0.0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "Descriptor %s is not an array of doubles", name);
    return(VM_FALSE);
  }
  
  if (len > desc->len) {
    maxLen = len;
  }
  else {
    maxLen = desc->len;
  }
  
  /* descriptor of the right name and type exists, copy value */
  for (i = 0; i < maxLen; i++) {
    val[i] = desc->descValue->dar[i];
  }

  if (comment)
    strcpy(comment, desc->descComment);
  
  return(VM_TRUE);
    
}

VimosBool readDoubleArrayDescFromTable(VimosTable *tab, const char *name, 
                                       double *val, char *comment, int len)
{
  char             modName[] = "readDoubleArrayDescFromTable";

  if (tab == NULL) {
    *val = 0;

    if (comment)
      strcpy(comment, "");

    cpl_msg_debug(modName, "NULL input table");
    return(VM_FALSE);
  }
  
  return(readDoubleArrayDescriptor(tab->descs, name, val, comment, len));
  
}


/* Get a pointer to the value of the column with name name */
VimosColumn *findColumn(VimosColumn *vCol, const char *name)
{
  VimosColumn  *tmpCol;
  
  tmpCol = vCol;
  
  while (tmpCol) {
    if ( !strcmp(tmpCol->colName, name) ) {
      return(tmpCol);
    }
    tmpCol = tmpCol->next;
  }
  
  return(NULL);
  
}


/* Get a pointer to the column with name name in a VimosTable 
 */
VimosColumn *findColInTab(VimosTable *table, const char *name)
{
  return(findColumn(table->cols, name));
}



/**
 * @memo
 *   Insert a HISTORY/COMMENT descriptor after other History/Comment
 *   descroptors in a descriptor list.
 *
 * @return VM_TRUE or VM_FALSE.
 *
 * @param desc Pointer to a descriptor list.
 * @param name Name of reference descriptor in the list.
 * @param newDesc Pointer to the descriptor to insert in the list.
 *
 * @doc
 *   Similar as InsertStringDescriptor BUT existing descriptors with 
 *   same name are
 *   NOT destroyed
 *
 * @author B. Garilli
 */

VimosBool insertHistoryDescriptor(VimosDescriptor **desc, const char *name, 
     const char *value, const char *comment)
{  
  VimosDescriptor        *tDesc;
  const char modName[]   = "insertHistoryDescriptor";
  int end_of_desc = 0 ;

  VimosDescriptor *newDesc   = newStringDescriptor(name, value, comment);

 if (newDesc) {

  if ( desc    == NULL) return VM_FALSE;
  if (*desc    == NULL) return VM_FALSE;
  if ( newDesc == NULL) return VM_FALSE;
 /*
  *  Now insert the new descriptor to the right position.
  */
  if ((tDesc = findDescriptor(*desc, name)) != NULL) {
    do {
      if (!tDesc->next) {
	end_of_desc = 1;
	break ;
      } else {
	tDesc = tDesc->next;
      }
    } while (strcmp(tDesc->descName,name) == 0);
    if (!end_of_desc) {
      newDesc->prev = tDesc->prev;
      newDesc->next = tDesc;
      tDesc->prev->next = newDesc;
      tDesc->prev = newDesc;
    } else {
      tDesc->next = newDesc;
      newDesc->prev = tDesc;
    }
    return VM_TRUE;
  } else { 
    if (!addDesc2Desc(newDesc, desc)) {
      cpl_msg_debug(modName, "The function addDesc2Desc has returned an error");
      return(VM_FALSE);
    }
    return VM_TRUE;
  }
 }else{
  return VM_FALSE;
 }
}

/**
 * @brief
 *   Locate a descriptor in a list.
 *
 * @param list  List of descriptors to search.
 * @param name  Descriptor name pattern to locate.
 *
 * @return Pointer to the first descriptor in @em list matching the
 *   name pattern @em pattern, or @c NULL if no matching descriptor
 *   is found.
 *
 * The function traverses the given descriptor list @em list looking for
 * the first occurrence of a descriptor with a name matching the name
 * pattern @em name. The name pattern may be either be just an ordinary
 * identifier string, or an extendend POSIX regular expression.
 */

VimosDescriptor *
vimosDscFind(VimosDescriptor *list, const char *name)
{

    int re_status;

    regex_t re;

    VimosDescriptor *dsc = NULL;


    assert(name != NULL);
    assert(list != NULL);

    re_status = regcomp(&re, name, REG_EXTENDED | REG_NOSUB);
    if (re_status) {
        return NULL;
    }

    while (!dsc && list) {
        if (regexec(&re, list->descName, 0, NULL, 0) == 0) {
            dsc = list;
        }

        list = list->next;
    }

    regfree(&re);

    return dsc;

}
    

/**
 * @brief
 *   Copy descriptors.
 *
 * @param tlist  Target descriptor list.
 * @param slist  Source descriptor list.
 * @param name   Descriptor name pattern.
 * @param hint   Position hint for descriptor insertion.
 *
 * @return The function returns @c EXIT_SUCCESS if the descriptors were
 *   copied successfully, or @c EXIT_FAILURE otherwise.
 *
 * The function copies all descriptors of the source list @em slist matching
 * the name pattern @em name to the target descriptor list @tlist. If the
 * position hint @em hint is not @c NULL the descriptors to copy are inserted
 * at this position. If hint is @c NULL the descriptors to copy are simply
 * appended to the target descriptor list.
 */

int
vimosDscCopy(VimosDescriptor **tlist, VimosDescriptor *slist,
             const char *name, const char *hint)
{

    int re_status;

    regex_t re;

    VimosDescriptor *position = NULL;


    assert(name != NULL);
    assert(tlist != NULL);
    assert(slist != NULL);

    re_status = regcomp(&re, name, REG_EXTENDED | REG_NOSUB);
    if (re_status) {
        return EXIT_FAILURE;
    }

    if (hint) {
        position = vimosDscFind(*tlist, hint);
    }

    while (slist) {
        if (regexec(&re, slist->descName, 0, NULL, 0) == 0) {
            VimosDescriptor *clone = copyOfDescriptor(slist);

            if (!clone) {
                return EXIT_FAILURE;
            }

            if (position) {
                if (position->prev) {
                    clone->prev = position->prev;
                    clone->prev->next = clone;
                }
                else {
                    clone->prev = NULL;
                    *tlist = clone;
                }

                clone->next = position;
                position->prev = clone;
            }
            else {
                if (addDesc2Desc(clone, tlist) != VM_TRUE) {
                    return EXIT_FAILURE;
                }
            }
        }

        slist = slist->next;

    }

    regfree(&re);

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Erase descriptors from a list.
 *
 * @param list  Descriptor list.
 * @param name  Descriptor name pattern.
 *
 * @return The function returns the number of descriptors removed from the
 *   list of descriptors or -1 if an error occurred.
 *
 * The function removes all descriptors from the descriptor list @em list
 * matching the name pattern @em name.
 */

int
vimosDscErase(VimosDescriptor **list, const char *name)
{
    int re_status;
    int count = 0;

    regex_t re;

    VimosDescriptor *pos;


    assert(name != NULL);
    assert(list != NULL);


    re_status = regcomp(&re, name, REG_EXTENDED | REG_NOSUB);
    if (re_status) {
        return -1;
    }

    pos = *list;
    while (pos) {
        VimosDescriptor *node = pos;

        pos = pos->next;

        if (regexec(&re, node->descName, 0, NULL, 0) == 0) {
            
            /*
             * Extract the descriptor from the list
             */

            if (node->next) {
                node->next->prev = node->prev;
            }

            if (node->prev) {
                node->prev->next = node->next;
            }

            /*
             * Make node's successor the list head if node was the head
             * of the descriptor list.
             */

            if (node == *list) {
                *list = node->next;
            }

            node->next = NULL;
            node->prev = NULL;

            deleteDescriptor(node);
            ++count;
        }
    }

    regfree(&re);

    return count;
}


/**
 * @brief
 *   Retrieve the size of a column.
 *
 * @return The function returns the size of the column.
 *
 * @param column  Column object.
 *
 * @doc
 *   The function retrieves the number of rows of the column object
 *   \textbf{column}.
 *
 * @author R. Palsa
 */

int colGetSize(VimosColumn *column)
{

  return column->len;

}


/**
 * @brief
 *   Retrieve the size of a column.
 *
 * @return The function returns the size of the column. 
 * 
 * @param table  Table object.
 * @param name   Column name. 
 *
 * @doc
 *   The function retrieves the number of rows of the specified column
 *   (supposedly, the number of rows in the input table).
 *
 * @author C.Izzo
 */
 
int tblGetSize(VimosTable *table, const char *name)
{
 
  return colGetSize(findColInTab(table, name));

}


/**
 * @brief
 *   Set the label of a table column.
 *
 * @return The function returns #EXIT_SUCCESS# if no error occurred, otherwise
 *   the return value is #EXIT_FAILURE#.
 *
 * @param column  Table column object.
 * @param name    Column label.
 *
 * @doc
 *   The function replaces the current label of column \textbf{column}, if
 *   present, with the string \textbf{name}.
 *
 * @author R. Palsa
 */

int colSetName(VimosColumn *column, const char *name)
{

  size_t sz;

  if (!column)
    return EXIT_FAILURE;

  if ((sz = strlen(name)) > VM_DESC_LENGTH - 1)
    return EXIT_FAILURE;

  memcpy(column->colName, name, sz);
  *(column->colName + sz) = '\0';

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the reference to the data of a string column.
 *
 * @return The function returns a pointer to the column data if no error
 *   occurred, otherwise the function returns a #NULL# pointer.
 *
 * @param column Column object.
 *
 * @doc
 *   The function first checks the type of \textbf{column}. If \textbf{column}
 *   is not of type #VM_STRING# the function returns an error. If the column
 *   has the correct type a handle to the data stored in \textbf{column}
 *   is returned, if a data array is present, otherwise an error is returned.
 *
 * @author R. Palsa
 */

char **colGetStringData(VimosColumn *column)
{

  assert(column != 0);

  if (column->colType != VM_STRING || !column->colValue->sArray)
    return 0;

  return column->colValue->sArray;

}


/**
 * @brief
 *   Get the reference to the data of a integer column.
 *
 * @return The function returns a pointer to the column data if no error
 *   occurred, otherwise the function returns a #NULL# pointer.
 *
 * @param column Column object.
 *
 * @doc
 *   The function first checks the type of \textbf{column}. If \textbf{column}
 *   is not of type #VM_INT# the function returns an error. If the column
 *   has the correct type a handle to the data stored in \textbf{column}
 *   is returned, if a data array is present, otherwise an error is returned.
 *
 * @author R. Palsa
 */

int *colGetIntData(VimosColumn *column)
{

  assert(column != 0);
  if (column->colType != VM_INT || !column->colValue->iArray)
    return 0;

  return column->colValue->iArray;

}


/**
 * @brief
 *   Get the reference to the data of a float column.
 *
 * @return The function returns a pointer to the column data if no error
 *   occurred, otherwise the function returns a #NULL# pointer.
 *
 * @param column Column object.
 *
 * @doc
 *   The function first checks the type of \textbf{column}. If \textbf{column}
 *   is not of type #VM_FLOAT# the function returns an error. If the column
 *   has the correct type a handle to the data stored in \textbf{column}
 *   is returned, if a data array is present, otherwise an error is returned.
 *
 * @author R. Palsa
 */

float *colGetFloatData(VimosColumn *column)
{

  assert(column != 0);

  if (column->colType != VM_FLOAT || !column->colValue->fArray)
    return 0;

  return column->colValue->fArray;

}


/**
 * @brief
 *   Get the reference to the data of a double column.
 *
 * @return The function returns a pointer to the column data if no error
 *   occurred, otherwise the function returns a #NULL# pointer.
 *
 * @param column Column object.
 *
 * @doc
 *   The function first checks the type of \textbf{column}. If \textbf{column}
 *   is not of type #VM_DOUBLE# the function returns an error. If the column
 *   has the correct type a handle to the data stored in \textbf{column}
 *   is returned, if a data array is present, otherwise an error is returned.
 *
 * @author R. Palsa
 */

double *colGetDoubleData(VimosColumn *column)
{

  assert(column != 0);

  if (column->colType != VM_DOUBLE || !column->colValue->dArray)
    return 0;

  return column->colValue->dArray;

}


/**
 * @brief
 *   Get the reference to the data of a string column in a table.
 *
 * @return The function returns a pointer to the column data if no error
 *   occurred, otherwise the function returns a #NULL# pointer.
 *
 * @param table  Table object.
 * @param name   Column label.
 *
 * @doc
 *   The function queries the table object \textbf{table} for the column
 *   \textbf{name}. If the column exists and if it is of the type #VM_STRING#,
 *   a handle to the columns string data is returned. In case the column
 *   does not exist, if it is not of the appropriate type, or if there
 *   is no data array assocoated to the column, the function returns
 *   an error.
 *
 *   Note that this function does not provide any size information about
 *   the data array for which the handle is returned.
 *
 * @author R. Palsa
 */

char **tblGetStringData(VimosTable *table, const char *name)
{

  VimosColumn *column;


  assert(table != 0 && name != 0);

  if (!(column = findColInTab(table, name)))
    return 0;

  return colGetStringData(column);

}


/**
 * @brief
 *   Get the reference to the data of an integer column in a table.
 *
 * @return The function returns a pointer to the column data if no error
 *   occurred, otherwise the function returns a #NULL# pointer.
 *
 * @param table  Table object.
 * @param name   Column label.
 *
 * @doc
 *   The function queries the table object \textbf{table} for the column
 *   \textbf{name}. If the column exists and if it is of the type #VM_INT#,
 *   a handle to the columns string data is returned. In case the column
 *   does not exist, if it is not of the appropriate type, or if there
 *   is no data array assocoated to the column, the function returns
 *   an error.
 *
 *   Note that this function does not provide any size information about
 *   the data array for which the handle is returned.
 *
 * @author R. Palsa
 */

int *tblGetIntData(VimosTable *table, const char *name)
{

  VimosColumn *column;


  assert(table != 0 && name != 0);

  if (!(column = findColInTab(table, name))) 
    return 0;
 

  return colGetIntData(column);

}


/**
 * @brief
 *   Get the reference to the data of a float column in a table.
 *
 * @return The function returns a pointer to the column data if no error
 *   occurred, otherwise the function returns a #NULL# pointer.
 *
 * @param table  Table object.
 * @param name   Column label.
 *
 * @doc
 *   The function queries the table object \textbf{table} for the column
 *   \textbf{name}. If the column exists and if it is of the type #VM_FLOAT#,
 *   a handle to the columns string data is returned. In case the column
 *   does not exist, if it is not of the appropriate type, or if there
 *   is no data array assocoated to the column, the function returns
 *   an error.
 *
 *   Note that this function does not provide any size information about
 *   the data array for which the handle is returned.
 *
 * @author R. Palsa
 */

float *tblGetFloatData(VimosTable *table, const char *name)
{

  VimosColumn *column;


  assert(table != 0 && name != 0);

  if (!(column = findColInTab(table, name)))
    return 0;

  return colGetFloatData(column);

}


/**
 * @brief
 *   Get the reference to the data of a double column in a table.
 *
 * @return The function returns a pointer to the column data if no error
 *   occurred, otherwise the function returns a #NULL# pointer.
 *
 * @param table  Table object.
 * @param name   Column label.
 *
 * @doc
 *   The function queries the table object \textbf{table} for the column
 *   \textbf{name}. If the column exists and if it is of the type #VM_DOUBLE#,
 *   a handle to the columns string data is returned. In case the column
 *   does not exist, if it is not of the appropriate type, or if there
 *   is no data array assocoated to the column, the function returns
 *   an error.
 *
 *   Note that this function does not provide any size information about
 *   the data array for which the handle is returned.
 *
 * @author R. Palsa
 */

double *tblGetDoubleData(VimosTable *table, const char *name)
{

  VimosColumn *column;


  assert(table != 0 && name != 0);

  if (!(column = findColInTab(table, name)))
    return 0;

  return colGetDoubleData(column);

}


/**
 * @brief
 *   Assign a string value to a table entry.
 *
 * @return The function returns #EXIT_SUCCESS# if no error occurred,
 *   otherwise the return value is #EXIT_FAILURE#.
 *
 * @param table  Table object.
 * @param name   Column label.
 * @param row    Row number.
 * @param value  Value string to be set.
 *
 * @doc
 *   The function assigns the string \textbf{value} to the table cell 
 *   specified by the column label \textbf{name} and the row number
 *   \textbf{row} of the target table \textbf{table}. The function returns
 *   an error if either the target column is not found in \textbf{table}
 *   or the provided row number is out of range, i.e larger than the
 *   column capacity.
 *
 *   If there is already a string assigned to this table cell, this string
 *   is deallocated before \textbf{value} is assigned. If #NULL# is passed
 *   as \textbf{value}, the table cell will be set to #NULL#. The assignent
 *   is done by copying the value string \textbf{value}.
 *
 * @author R. Palsa
 */

int tblSetStringValue(VimosTable *table, const char *name, unsigned int row,
		      const char *value)
{

  VimosColumn *column;


  assert(table != 0);
  assert(name != 0);


  if (!(column = findColInTab(table, name)) || (int)row > column->len)
    return EXIT_FAILURE;


  /*
   * If the column entry had already a value associated it is deallocated
   * first.
   */

  if (column->colValue->sArray[row])
    cpl_free(column->colValue->sArray[row]);

  if (!value)
    column->colValue->sArray[row] = 0;
  else
    column->colValue->sArray[row] = cpl_strdup(value);

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Assign an integer value to a table entry.
 *
 * @return The function returns #EXIT_SUCCESS# if no error occurred,
 *   otherwise the return value is #EXIT_FAILURE#.
 *
 * @param table  Table object.
 * @param name   Column label.
 * @param row    Row number.
 * @param value  Value to be set.
 *
 * @doc
 *   The function assigns the integer \textbf{value} to the table cell 
 *   specified by the column label \textbf{name} and the row number
 *   \textbf{row} of the target table \textbf{table}. The function returns
 *   an error if either the target column is not found in \textbf{table}
 *   or the provided row number is out of range, i.e larger than the
 *   column capacity.
 *
 * @author R. Palsa
 */

int tblSetIntValue(VimosTable *table, const char *name, unsigned int row,
		   int value)
{

  VimosColumn *column;


  assert(table != 0);
  assert(name != 0);


  if (!(column = findColInTab(table, name)) || (int)row > column->len)
    return EXIT_FAILURE;


  column->colValue->iArray[row] = value;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Assign a float value to a table entry.
 *
 * @return The function returns #EXIT_SUCCESS# if no error occurred,
 *   otherwise the return value is #EXIT_FAILURE#.
 *
 * @param table  Table object.
 * @param name   Column label.
 * @param row    Row number.
 * @param value  Value to be set.
 *
 * @doc
 *   The function assigns the float \textbf{value} to the table cell 
 *   specified by the column label \textbf{name} and the row number
 *   \textbf{row} of the target table \textbf{table}. The function returns
 *   an error if either the target column is not found in \textbf{table}
 *   or the provided row number is out of range, i.e larger than the
 *   column capacity.
 *
 * @author R. Palsa
 */

int tblSetFloatValue(VimosTable *table, const char *name, unsigned int row,
		     float value)
{

  VimosColumn *column;


  assert(table != 0);
  assert(name != 0);


  if (!(column = findColInTab(table, name)) || (int)row > column->len)
    return EXIT_FAILURE;


  column->colValue->fArray[row] = value;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Assign a double value to a table entry.
 *
 * @return The function returns #EXIT_SUCCESS# if no error occurred,
 *   otherwise the return value is #EXIT_FAILURE#.
 *
 * @param table  Table object.
 * @param name   Column label.
 * @param row    Row number.
 * @param value  Value to be set.
 *
 * @doc
 *   The function assigns the double \textbf{value} to the table cell 
 *   specified by the column label \textbf{name} and the row number
 *   \textbf{row} of the target table \textbf{table}. The function returns
 *   an error if either the target column is not found in \textbf{table}
 *   or the provided row number is out of range, i.e larger than the
 *   column capacity.
 *
 * @author R. Palsa
 */

int tblSetDoubleValue(VimosTable *table, const char *name, unsigned int row,
		      double value)
{

  VimosColumn *column;


  assert(table != 0);
  assert(name != 0);


  if (!(column = findColInTab(table, name)) || (int)row > column->len)
    return EXIT_FAILURE;


  column->colValue->dArray[row] = value;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Create a copy of a table column.
 *
 * @return The function returns a pointer to the created copy of the column if
 *   no error occurred, otherwise a #NULL# pointer is returned.
 *
 * @param table  Table object.
 * @param name   Name of the column to be copied.
 *
 * @doc
 *   The function creates a new table column having the same name
 *   \textbf{name} as the source column found in \textbf{table}.
 *   The data of the source columns is copied to the newly created
 *   column object.
 *
 * @author R. Palsa
 */

VimosColumn *tblCopyColumn(VimosTable *table, const char *name)
{

  register int i;

  size_t sz;

  VimosColumn *column, *source;


  if (!table || !name)
    return 0;

  if (!(source = findColInTab(table, name)))
    return 0;


  /*
   * Create new column object with the appropropriate name and type
   * and size. Check the size of the name first. Note that the column
   * constructor creates the name member with a size of VM_DESC_LENGTH
   * characters, i.e. the name string may be made of VM_DESC_LENGTH - 1
   * characters at most!
   *
   * Updating column->next and column->prev with 0 is not needed. The
   * constructor takes care of this.
   */

  if (!(column = newColumn()))
    return 0;

  if ((sz = strlen(name)) > VM_DESC_LENGTH - 1) {
    deleteColumn(column);
    return 0;
  }
  else {
    memcpy(column->colName, source->colName, sz);
    *(column->colName + sz) = '\0';
  }

  column->colType = source->colType;
  column->len = source->len;
 

  /*
   * Allocate memory for the column data
   */

  if (column->len > 0) {
    size_t nbytes = column->len;

    switch (column->colType) {
    case VM_INT:
      nbytes *= sizeof(int);
      break;

    case VM_FLOAT:
      nbytes *= sizeof(float);
      break;

    case VM_DOUBLE:
      nbytes *= sizeof(double);
      break;

    case VM_STRING:
      nbytes *= sizeof(char *);
      break;

    default:
      deleteColumn(column);
      return 0;
      break;
    }

    if ((column->colValue->p = cpl_malloc(nbytes * sizeof(char))) == 0) {
      deleteColumn(column);
      return 0;
    }


    /*
     * Copy the data. Note that a string column is different and needs
     * special treatment.
     */

    memcpy(column->colValue->p, source->colValue->p, nbytes);

    if (source->colType == VM_STRING) {
      char **s = column->colValue->sArray;

      for (i = 0; i < source->len; i++)
	*s++ = cpl_strdup(source->colValue->sArray[i]);
    }
  }

  return column;

}


/**
 * @brief
 *   Remove a column from a table.
 *
 * @return The function returns a pointer to the removed column if no error
 *   occurred, otherwise a #NULL# pointer is returned.
 *
 * @param table  Table object.
 * @param name   Name of the table column to be removed.
 *
 * @doc
 *   The function looks for the column \textbf{name} in the input table
 *   \textbf{table} and extracts it from the table. The extracted column
 *   is not destroyed but retured to the caller. The table object is
 *   updated with the new number of columns.
 *
 * @author R. Palsa
 */

VimosColumn *tblRemoveColumn(VimosTable *table, const char *name)
{

  VimosColumn *column;


  if (!table)
    return 0;

  if (!(column = findColInTab(table, name)))
    return 0;


  /*
   * Extract the column
   */

  if (column->prev)
    column->prev->next = column->next;

  if (column->next) {

    /*
     * Make sure that the column pointer points to the second column
     * if the first one was removed.
     */

    if (!(column->next->prev = column->prev))
      table->cols = column->next;
  }
  column->prev = 0;
  column->next = 0;

  table->numColumns--;

  return column;

}


/**
 * @brief
 *   Append a column to a table.
 *
 * @return The function returns #EXIT_SUCCESS# if no error occurred, otherwise
 *   #EXIT_FAILURE# is returned.
 *
 * @param table   Table object.
 * @param column  Column object.
 *
 * @doc
 *   The function appends the column object \textbf{column} to the table
 *   \textbf{table} and the number of columns in the table is incremented.
 *   
 * @author R. Palsa
 */

int tblAppendColumn(VimosTable *table, VimosColumn *column)
{

  VimosColumn  *lastcol;

  assert(table != 0 && column != 0);
  assert(column->prev == 0 && column->next == 0);

  if (!table->cols)
    table->cols = column;
  else {
    lastcol = table->cols;
    while (lastcol->next) lastcol = lastcol->next;
    lastcol->next = column;
    column->prev = lastcol;
  }


  /*
   * Adjust column counter
   */

  table->numColumns++;

  return EXIT_SUCCESS;

}
      
VimosBool readDescsFromFitsTable(VimosDescriptor **desc, fitsfile *tableID)
{
  int i;
  int numDescs;
  int status, keypos, hdutype;

  char name[FLEN_KEYWORD], value[FLEN_VALUE], value1[FLEN_VALUE]; 
  char comment[FLEN_COMMENT];
  char type[1];
  
  VimosDescriptor* tDesc=NULL;
  VimosDescriptor* lastDesc;
  char             modName[] = "readDescsFromFitsTable";

  lastDesc = NULL;
  
  if (!desc) {
    cpl_msg_debug(modName, "NULL input descriptor");
    return(VM_FALSE);
  }
  
  status = 0;

  if (tableID == NULL) {
    cpl_msg_debug(modName, "No pointer to fits file");
    return(VM_FALSE);
  }

  hdutype = tableID->Fptr->hdutype;
  if (hdutype != BINARY_TBL){
    cpl_msg_debug(modName, "This HDU is not a binary table");
    return(VM_FALSE);
  }

  /* get no. of keywords */
  if (fits_get_hdrpos(tableID, &numDescs, &keypos, &status)) {
    cpl_msg_debug(modName, 
    "The function fits_get_hdrpos has returned an error (code %d)", status);
    return(VM_FALSE);
  }

  for (i = 1; i <= numDescs; i++) {

    if (fits_read_keyn(tableID, i, name, value, comment, &status)) {
      cpl_msg_debug(modName, 
      "The function fits_read_keyn has returned an error (code %d)", status);
      return(VM_FALSE);
    }

/* Check for blank lines */
    if ( strlen(name) == 0 ) {
      status = 0;
      continue;
    }

    if (strncmp("HISTORY",name,7) == 0 || strncmp("COMMENT",name,7) == 0) {
      tDesc = newStringDescriptor(name,comment, "");
      if (tDesc == NULL) {
	cpl_msg_debug(modName,
        "The function newStringDescriptor has returned NULL");
	return(VM_FALSE);
      }
    } else {

      if (fits_get_keytype (value, (type), &status)) {
	cpl_msg_debug(modName, 
        "The function fits_get_keytype returned an error (code %d)", status);
	return(VM_FALSE);
      }

    /* First check on sign, then on integer/real, then on string  */
    /*   to determine type */

      if (*(type) == 'F') {
	tDesc = newDoubleDescriptor(name, atof(value), comment);
	if (tDesc == NULL) {
	  cpl_msg_debug(modName, 
          "The function newDoubleDescriptor has returned NULL");
	  return(VM_FALSE);
	}
      } else if (*(type) == 'I') {
	tDesc = newIntDescriptor(name, atol(value), comment);
	if (tDesc == NULL) {
	  cpl_msg_debug(modName, 
          "The function newIntDescriptor has returned NULL");
	  return(VM_FALSE);
	}
      } else if (*(type) == 'L') {
	if (*(value) == 'T') {
	  tDesc = newBoolDescriptor(name, VM_TRUE, comment);
	  if (tDesc == NULL) {
	    cpl_msg_debug(modName, 
            "The function newBoolDescriptor has returned NULL");
	    return(VM_FALSE);
	  }
	}
	if (*(value) == 'F') {
	  tDesc = newBoolDescriptor(name,VM_FALSE, comment);
	  if (tDesc == NULL) {
	    cpl_msg_debug(modName,
            "The function newBoolDescriptor has returned NULL");
	    return(VM_FALSE);
	  }
	}
      } else if (*(type) == 'C') {
	if (ffc2s(value, (value1), &status)){ /* remove quotes from string */
	  cpl_msg_debug(modName,
          "The function ffc2s has returned an error (code %d)", status);
	  return(VM_FALSE);
	}
	tDesc = newStringDescriptor(name, value1, comment);
	if (tDesc == NULL) {
	  cpl_msg_debug(modName,
          "The function newStringDescriptor has returned NULL");
	  return(VM_FALSE);
	}
      } else {
	cpl_msg_debug(modName, "Unrecognized key type %s",type);
	return(VM_FALSE);
      }
    }
    if (tDesc != NULL) {
      if (lastDesc == NULL) {
	*desc = tDesc;
      }
      else {
	lastDesc->next = tDesc;
	tDesc->prev = lastDesc;
      }
      lastDesc = tDesc;
    }
  }

  return(VM_TRUE);
}
        
VimosBool writeDescsToFitsTable(VimosDescriptor *desc, fitsfile *tableID)
{
  char charBuf[68];
  
  int              fits_type, status;
  int              tf;
  VimosDescriptor* tDesc;
  char             modName[] = "writeDescsToFitsTable";
  
  tDesc = desc;
  
  if (!desc) {
    cpl_msg_debug(modName, "NULL input descriptor");
    return(VM_FALSE);
  }
  
  if (tableID == NULL) {
    cpl_msg_debug(modName, "No pointer to fits file");
    return(VM_FALSE);
  }

  status = 0;

  while (tDesc) {
    status = 0;
    switch(tDesc->descType) {
    case VM_INT : 
      {
        fits_type = TINT;
	if (fits_update_key(tableID,fits_type,tDesc->descName,
			   &(tDesc->descValue->i),tDesc->descComment,&status))
	  {
	    cpl_msg_debug(modName,
            "The function fits_update_key returned an error (code %d)", status);
	    return(VM_FALSE);
	  }
        break;
      }
    case VM_BOOL : 
      {
        fits_type = TLOGICAL;
        if (tDesc->descValue->b == VM_TRUE) {
          tf = 1;
        }
        else {
          tf = 0;
        }
	if (fits_update_key(tableID,fits_type,tDesc->descName,
			    &tf,tDesc->descComment,&status)) {
	  cpl_msg_debug(modName, 
          "The function fits_update_key returned an error (code %d)", status);
          return(VM_FALSE);
        }   
        break;
      }
    case VM_FLOAT :
      {
	fits_type = TFLOAT;
	if (fits_update_key(tableID,fits_type,tDesc->descName,
			    &(tDesc->descValue->f),tDesc->descComment,&status))
	  { 
	    cpl_msg_debug(modName, 
            "The function fits_update_key returned an error (code %d)", status);
	    return(VM_FALSE);
	  }
        break;
      }
    case VM_DOUBLE :
      {
	fits_type = TDOUBLE;
	if (fits_update_key(tableID,fits_type,tDesc->descName,
			   &(tDesc->descValue->d),tDesc->descComment,&status))
	  {
 	    cpl_msg_debug(modName, 
            "The function fits_update_key returned an error (code %d)", status);
	    return(VM_FALSE);
	  }
        break;
      }
    case VM_STRING :
      {
	fits_type = TSTRING;
	if (strncmp ((tDesc->descName),"HISTORY",7) == 0){
	  if (fits_write_history(tableID, tDesc->descValue->s, &status)) {
	    cpl_msg_debug(modName,
            "The function fits_write_history returned error code %d", status);
	    return(VM_FALSE);
	  }
	  break;
	} 
	if (strncmp ((tDesc->descName),"COMMENT",7) == 0){
	  if (fits_write_comment(tableID, tDesc->descValue->s, &status)) {
	    cpl_msg_debug(modName,
            "The function fits_write_comment returned error code %d", status);
	    return(VM_FALSE);
	  }
	  break;
	}
	strcpy(charBuf, (tDesc->descValue->s));
	if (fits_update_key(tableID,fits_type,tDesc->descName,
			    charBuf,tDesc->descComment,&status)) {
	  cpl_msg_debug(modName, 
          "The function fits_update_key returned error code %d", status);
	  return(VM_FALSE);
	}
	break;
      }
    default :
      {
        cpl_msg_debug(modName,
        "Unrecognized type of value stored in the descriptor");
        return(VM_FALSE);
      }
    }
    
    tDesc = tDesc->next;
  }

  return(VM_TRUE);
} 

/** 
 * @memo 
 *   Write a VimosTable into a disk FITS file
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param filename         name of the fitsfile
 * @param table            pointer to VimosTable to be written as a FITS file 
 * @param tabletype        Extension name of table (table type)
 *
 * @doc 
 *   Create a disk table FITS file from structure VimosTable
 *
 * @author C. Izzo
 */

VimosBool createFitsTable(char *filename, VimosTable *table, 
                        const char *tabletype)
{
  VimosImage  *nullImage;
  VimosColumn *currentColumn;

  char       **ttype;
  char       **tform;
  char       **tunit;
  char       **stArray;
  char         modName[] = "createFitsTable";
  int          colNumber=0, nrows=0, i=0;
  int          status = 0, maxLength=0, length=0, num;

  if (table) {

    /*
     * Create the empty primary array first
     */
    
    nullImage = newImage(0, 0, (float *) NULL);
    if (nullImage) {
      if ((openNewFitsImage(filename, nullImage))) {
	
	if (table->numColumns && table->cols) {
	  
	  /*
	   * Get information for constructing the FITS table extension
           */
	  
          nrows = table->cols->len;
          ttype = (char **) cpl_malloc(table->numColumns*sizeof(char *));
          tform = (char **) cpl_malloc(table->numColumns*sizeof(char *));
          tunit = (char **) cpl_malloc(table->numColumns*sizeof(char *));
          currentColumn = table->cols;
          for (colNumber=0; currentColumn; colNumber++) {
            ttype[colNumber] = currentColumn->colName;
            tunit[colNumber] = " ";
            switch (currentColumn->colType) {
            case VM_INT          : tform[colNumber] = "1J"; break;
            case VM_FLOAT        : tform[colNumber] = "1E"; break;
            case VM_DOUBLE       : tform[colNumber] = "1D"; break;
            case VM_STRING       :
              stArray = colGetStringData(currentColumn);
	      for (i = 0; i < nrows; i++) {
		length = strlen(stArray[i])+1;
		if (length > maxLength) maxLength = length;
	      }
	      num = floor(log10(maxLength))+1;
	      tform[colNumber] = (char*) cpl_calloc((num+1+1),sizeof(char));
              sprintf(tform[colNumber],"%dA",maxLength);
              break;

	    default              : 
              cpl_msg_debug(modName, "Unsupported table column type");
              return(VM_FALSE);
            }
            currentColumn = currentColumn->next;
          }
         /*
          * create HDU for BINARY TABLE
          */
          fits_create_tbl(nullImage->fptr, BINARY_TBL,
          (long) nrows, table->numColumns, ttype,
          tform, tunit, (char *) tabletype, &status);
        }
        else {
          fits_create_tbl(nullImage->fptr, BINARY_TBL,
          (long) 0, 0, NULL,
          NULL, NULL, (char *) tabletype, &status);
        }

        if (!status) {

         /*
          * Write descriptor header first (faster)
          */

	  /* BG: if IN the table there are the NAXIS, NAXIS1 and NAXIS2
	     descriptors and the various fits standard table descriptors
	     (e.g. because it has been created starting from 
	     another table read from disk) the write fails with error 241
	     NAXIS descriptors must be removed before, if existing 
	  Remove also EXTNAME, if existing, thus it will be written
	  according to PRO TABLE 
	  Remove also any already existing TTYPE, TFORM, TUNIT. 
	  These will be correctly created by cfitsio according to current table*/

	  deleteSetOfDescriptors(&(table->descs),"NAXIS*");
	  deleteSetOfDescriptors(&(table->descs),"*COUNT");
	  deleteSetOfDescriptors(&(table->descs),"TUNIT*");
	  deleteSetOfDescriptors(&(table->descs),"TFIELDS*");
	  deleteSetOfDescriptors(&(table->descs),"EXTNAME");
	  deleteSetOfDescriptors(&(table->descs),"TTYPE*");
	  deleteSetOfDescriptors(&(table->descs),"TFORM*");
	  deleteSetOfDescriptors(&(table->descs),"TUNIT*");

          if (writeDescsToFitsTable(table->descs, nullImage->fptr)) {

            if (table->numColumns && table->cols) {
             /*
              * Dumping columns into file.
              */
              currentColumn = table->cols;
              for (colNumber=1; currentColumn; colNumber++) {
                switch (currentColumn->colType) {
                case VM_INT    : fits_write_col(nullImage->fptr,TINT,
                                   colNumber, 1, 1, (long) currentColumn->len,
                                   currentColumn->colValue->iArray, &status);
                                   break;
                case VM_FLOAT  : fits_write_col(nullImage->fptr,TFLOAT,
                                   colNumber, 1, 1, (long) currentColumn->len,
                                   currentColumn->colValue->fArray, &status);
                                   break;
                case VM_DOUBLE : fits_write_col(nullImage->fptr,TDOUBLE,
                                   colNumber, 1, 1, (long) currentColumn->len,
                                   currentColumn->colValue->dArray, &status);
                                   break;
                case VM_STRING : fits_write_col(nullImage->fptr,TSTRING,
                                   colNumber, 1, 1, (long) currentColumn->len,
                                   currentColumn->colValue->sArray, &status);
                                   break;
                default        : break;
                }
                currentColumn = currentColumn->next;
              }
            }
            if (closeFitsImage(nullImage, 0)) {
              cpl_msg_info(modName,"Table %s (%s) created.",filename,tabletype);
              return(VM_TRUE);
            }
          }
        }
      }
    }
  }
  return(VM_FALSE);
}

/** 
 * @memo 
 *   Read a fits table file (single extension) into a VimosTable 
 *
 * @return VimosTable structure
 *
 * @param tableName        name of the fitsfile 
 * @param IOflag           flag to only read or also modify the fitsfile 
 *
 * @doc 
 *   Create a VimosTable from a table Fits file on disk
 *
 * @author P. Sartoretti
 * BG: added reading of table Name
 * BG: added numColumns
 */

VimosTable *openOldFitsTable(const char *tableName, int IOflag)
{
  int          status = 0; 
  int          typecode = 0;
  long int     width = 0;
  long int     repeat=0;
  int          null, nCols, nRows, i, j, dummyInt;
  VimosTable  *tTable;
  VimosColumn *currentColumn;
  VimosColumn *holder=NULL;
  char       **ttype;
  char         comment[80];
  char         modName[] = "openOldFitsTable";
 
  tTable = newTable();

  if (IOflag == 0) {
    ffopen(&(tTable->fptr), tableName, READONLY, &status);
  }
  if (IOflag == 1) {
    ffopen(&(tTable->fptr), tableName, READWRITE, &status);
  }
  cpl_msg_debug(modName, "Table opening exit status = %d", status);

  fits_movrel_hdu(tTable->fptr, 1, NULL, &status);
  if (status) {
    cpl_msg_debug(modName, "No table extension found");
    return NULL;
  }
  
  readDescsFromFitsTable(&(tTable->descs), tTable->fptr);
  
  readIntDescriptor(tTable->descs, "TFIELDS", &nCols, comment);
  readIntDescriptor(tTable->descs, "NAXIS2", &nRows, comment);
  readStringDescriptor(tTable->descs, "ESO PRO TABLE", tTable->name, comment);

  if (nCols) {
    tTable->numColumns = nCols;
    ttype = (char **) cpl_malloc(nCols*sizeof(char *));
    for (i = 0; i < nCols; i++) {   
     /* 
      * Allocate space for the column labels 
      */
      ttype[i] = (char *) cpl_malloc(FLEN_VALUE);
    }
   /* 
    * Read the column names from the TTYPEn keywords 
    */
    fits_read_keys_str(tTable->fptr, "TTYPE", 1, nCols, ttype, &dummyInt,
                       &status);
    if(status) {
      cpl_msg_debug(modName, "Problems in reading column names");
      return NULL;
    }
    
   /* 
    * Allocate and fill table columns 
    */ 
    for (i=0; i<nCols; i++) {
     fits_get_coltype(tTable->fptr, i+1, &typecode, &repeat, &width, &status);
     if(status) {
       cpl_msg_debug(modName, "Cannot read column type");
       return NULL;
     }
     
     switch(typecode) {
         case TSHORT: 
         case TINT:
         case TLONG: 
           currentColumn = newIntColumn(nRows, ttype[i]);
           fits_read_col(tTable->fptr, TINT, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->iArray, &null, &status);
           break;
         case TFLOAT:  
           currentColumn = newFloatColumn(nRows, ttype[i]);
           fits_read_col(tTable->fptr, TFLOAT, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->iArray, &null, &status);
           break;
         case TDOUBLE:  
           currentColumn = newDoubleColumn(nRows, ttype[i]);
           fits_read_col(tTable->fptr, TDOUBLE, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->iArray, &null, &status);
           break;
         case TSTRING:
           currentColumn = newStringColumn(nRows, ttype[i]);
           for (j = 0; j < nRows; j++) {
             currentColumn->colValue->sArray[j] =
               (char *) cpl_malloc(((int) repeat + 1)*sizeof(char));
             if (currentColumn->colValue->sArray[j] == NULL) {
               deleteColumn(currentColumn);
               cpl_msg_error(modName, "Allocation Error");
               return(NULL);
             }
           }

           fits_read_col(tTable->fptr, TSTRING, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->sArray, &null, &status);
           break;
         default:
           cpl_msg_debug(modName, "Unsupported table column type");
           return NULL;
           break;
     }  
     if (status) {
       cpl_msg_debug(modName, "Cannot read column data");
       return NULL;
     }
     if (i == 0) {
       tTable->cols = currentColumn;
     }
     else {
       holder->next = currentColumn;
       currentColumn->prev = holder;
     }
     holder = currentColumn;
   }
 }
 return tTable;
}

/** 
 * @memo 
 *   Read a table FITS file (already open) and create a VimosTable
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param tTable         pointer to output VimosTable  
 * @param tableID        pointer to the table fitsfile
 *
 * @doc 
 *   Read a table FITS file (already open) and create a VimosTable
 *     (similar to openOldFitsTable (openOldFitsTable will disapear): 
 *     exept that this function oly 
 *  reads the table fits file (already open) and fills the VimosTabl 
 *
 * @author P. Sartoretti
 * BG: added reading of table Name
 */

VimosBool readFitsTable(VimosTable *tTable, fitsfile *tableID)
{
  int          status = 0; 
  int          typecode = 0;
  long int     width = 0;
  long int     repeat=0;
  int          null, nCols, nRows, i, j, dummyInt, hdutype;
  VimosColumn *currentColumn;
  VimosColumn *holder=NULL;
  char       **ttype;
  char         comment[80];
  char         modName[] = "readFitsTable";
 

  if (tTable == NULL) {
    cpl_msg_error(modName, "Null input Table");
    return VM_FALSE;
  }
  if(!tableID) {
    cpl_msg_error(modName, "Null pointer to fitsfile");
    return VM_FALSE;
  }
   
  hdutype = tableID->Fptr->hdutype;
  if (hdutype != BINARY_TBL){
    cpl_msg_error(modName, "This HDU is not a binary table");
    return(VM_FALSE);
  }

  

  if(!readDescsFromFitsTable(&(tTable->descs), tableID)) {
    cpl_msg_error(modName, "Error in reading descriptors");
    return VM_FALSE;
  }
  
  readIntDescriptor(tTable->descs, "TFIELDS", &nCols, comment);
  readIntDescriptor(tTable->descs, "NAXIS2", &nRows, comment);
  readStringDescriptor(tTable->descs, "ESO PRO TABLE", tTable->name, comment);

  if (nCols) {
    tTable->numColumns = nCols;
    ttype = (char **) cpl_malloc(nCols*sizeof(char *));
    for (i = 0; i < nCols; i++) {   
     /* 
      * Allocate space for the column labels 
      */
      ttype[i] = (char *) cpl_malloc(FLEN_VALUE);
    }
   /* 
    * Read the column names from the TTYPEn keywords 
    */
    fits_read_keys_str(tableID, "TTYPE", 1, nCols, ttype, &dummyInt,
                       &status);
    if(status) {
      cpl_msg_debug(modName, "Problems in reading column names");
      return VM_FALSE;
    }
    
   /* 
    * Allocate and fill table columns 
    */ 
    for (i=0; i<nCols; i++) {
     fits_get_coltype(tableID, i+1, &typecode, &repeat, &width, &status);
     if(status) {
       cpl_msg_debug(modName, "Cannot read column type");
       return VM_FALSE;
     }
     
     switch(typecode) {
         case TSHORT: 
         case TINT:
         case TLONG: 
           currentColumn = newIntColumn(nRows, ttype[i]);
           fits_read_col(tableID, TINT, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->iArray, &null, &status);
           break;
         case TFLOAT:  
           currentColumn = newFloatColumn(nRows, ttype[i]);
           fits_read_col(tableID, TFLOAT, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->iArray, &null, &status);
           break;
         case TDOUBLE:  
           currentColumn = newDoubleColumn(nRows, ttype[i]);
           fits_read_col(tableID, TDOUBLE, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->iArray, &null, &status);
           break;
         case TLOGICAL: 
           currentColumn = newCharacterColumn(nRows, ttype[i]);
           fits_read_col(tableID, TBYTE, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->cArray, &null, &status);
           break;
         case TSTRING:
           currentColumn = newStringColumn(nRows, ttype[i]);
           for (j = 0; j < nRows; j++) {
             currentColumn->colValue->sArray[j] =
               (char *) cpl_malloc(((int)repeat + 1)*sizeof(char));
             if (currentColumn->colValue->sArray[j] == NULL) {
               deleteColumn(currentColumn);
               cpl_msg_error(modName, "Allocation Error");
               return(VM_FALSE);
             }
           }

           fits_read_col(tableID, TSTRING, i+1, 1, 1, (long) nRows, NULL,
                         currentColumn->colValue->sArray, &null, &status);
           break;
         default:
             cpl_msg_warning(cpl_func, "en default");
           cpl_msg_debug(modName, "Unsupported table column type");
           return VM_FALSE;
           break;
     }  
     if (status) {
       cpl_msg_debug(modName, "Cannot read column data");
       return VM_FALSE;
     }
     if (i == 0) {
       tTable->cols = currentColumn;
     }
     else {
       holder->next = currentColumn;
       currentColumn->prev = holder;
     }
     holder = currentColumn;
   }
 }
 return VM_TRUE;
}


/** 
 * @memo 
 *   Close a Fits table 
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param table           pointer to a VimosTable
 * @param flag            flag to write (=1) or not (=0) the
 *                        VimosTable to the fitsFile to close.
 * @doc 
 *   Close a FITS table
 *
 * @author P. Sartoretti
 */

VimosBool closeFitsTable(VimosTable *table, int flag)
{
  int status = 0;
  int ii, nCols;
  VimosColumn *currentColumn;
  
  if (flag ==1) {
    if (writeDescsToFitsTable(table->descs, table->fptr)) {
      
     /*
      * Dump columns into fits file and close it
      */
      nCols = table->numColumns;
      currentColumn = table->cols;
      for (ii=1; ii<=nCols; ii++) {
        switch (currentColumn->colType) {
            case VM_INT:
              fits_write_col(table->fptr, TINT, ii, (long) 1, (long) 1,
                             (long) currentColumn->len,
                             currentColumn->colValue->iArray, 
                             &status);
              break;
            case VM_FLOAT:
              fits_write_col(table->fptr, TFLOAT, ii, (long) 1, (long) 1,
                             (long) currentColumn->len,
                             currentColumn->colValue->fArray, 
                             &status);
              break;
            case VM_DOUBLE:
              fits_write_col(table->fptr, TDOUBLE, ii, (long) 1, (long) 1,
                             (long) currentColumn->len,
                             currentColumn->colValue->dArray, 
                             &status);
              break;
            case VM_STRING:
              fits_write_col_str(table->fptr, ii, (long) 1, (long) 1,
                                 (long) currentColumn->len,
                                 currentColumn->colValue->sArray, 
                                 &status);
              break; 
            default:
              break;
        }
        currentColumn = currentColumn->next;
      }
    }
  }   
  status = 0;
  if(fits_close_file(table->fptr, &status))
    return(VM_FALSE);
  else
    return VM_TRUE;
}

/** 
 * @memo 
 *   Copy relevant informations to primary header
 *
 * @return VM_TRUE or VM_FALSE. 
 *
 * @param fileName         name of the fitsfile
 * @param keyName          name of the keyword to be copied  
 *
 * @doc 
 *   Copy a keyword from secondary to primary header
 *
 * @author C. Izzo
 */

int copyToPrimary(char *fileName, const char *keyName)
{
  const char modName[] = "copyToPrimary";
  int        status = 0;
  char       card[FLEN_CARD];
  fitsfile  *fptr;
  
  if (ffopen(&fptr, fileName, READWRITE, &status)) {
    cpl_msg_error(modName, "Failure in opening file %s", fileName);
    return EXIT_FAILURE;
  }
  if (fits_movabs_hdu(fptr, 2, NULL, &status)) {
    cpl_msg_error(modName, 
    "Failure in moving to first extension of file %s", fileName);
    return EXIT_FAILURE;
  }
  if (fits_read_card(fptr, (char *)keyName, card, &status)) {
    cpl_msg_error(modName, "Keyword %s not found", keyName);
    return EXIT_FAILURE;
  }
  if (fits_movabs_hdu(fptr, 1, NULL, &status)) {
    cpl_msg_error(modName, 
    "Failure in moving to primary array of file %s", fileName);
    return EXIT_FAILURE;
  }
  if (fits_write_record(fptr, card, &status)) {
    cpl_msg_error(modName, 
    "Failure in writing record to primary array of file %s", fileName);
    return EXIT_FAILURE;
  }
  fits_close_file(fptr, &status);
  return EXIT_SUCCESS;
}
/*@}*/
