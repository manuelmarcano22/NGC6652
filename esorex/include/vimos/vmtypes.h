/* $Id: vmtypes.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_TYPES_H
#define VM_TYPES_H

#include <stdlib.h>

#include <pilmacros.h>


PIL_BEGIN_DECLS

#ifdef _DEC_ALPHA
typedef unsigned int VimosUlong32 ;
typedef int VimosLong32 ;
#else
typedef unsigned long VimosUlong32 ;
typedef long VimosLong32 ;
#endif

typedef unsigned short VimosUshort16 ;
typedef short VimosShort16 ;

typedef unsigned char VimosUchar8 ;
typedef char VimosChar8 ;
 
typedef unsigned char VimosByte ;

/* defined in limits.h, redefined here for portability  */

#define LONG32_MIN  (long32)(-2147483647-1) 
#define LONG32_MAX  (long32)(2147483647)
#define ULONG32_MAX (ulong32)(4294967295)

#define SHRT16_MIN  (short16)(-32768)
#define SHRT16_MAX  (short16)(32767)
#define USHRT16_MAX (ushort16)(65535)

/* Boolean type */
typedef enum _VIMOS_BOOL_
{
  VM_FALSE = 0,
  VM_TRUE = 1
} VimosBool ;


typedef enum _VIMOS_QUAD_
{
  VM_QUAD_UNDEF = 0,
  VM_QUAD_1 = 1,
  VM_QUAD_2 = 2,
  VM_QUAD_3 = 3,
  VM_QUAD_4 = 4
} VimosQuad;

typedef enum _VIMOS_OBS_TYPE_ 
{
  VM_OBS_UNDEF = 0,
  VM_OBS_MOS = 1,
  VM_OBS_IFU = 2,
  VM_OBS_IMA = 3
} VimosObsType;


typedef enum _VIMOS_VAR_TYPE_ 
{
  VM_VARTYPE_UNDEF = 0,
  VM_INT,
  VM_BOOL,
  VM_FLOAT,
  VM_DOUBLE,
  VM_CHARACTER,
  VM_STRING,
  VM_POINTER,
  VM_INT_ARRAY,
  VM_FLOAT_ARRAY,
  VM_DOUBLE_ARRAY,
  VM_INT_MATRIX,
  VM_FLOAT_MATRIX,
  VM_DOUBLE_MATRIX,
  VM_INT_CUBE,
  VM_FLOAT_CUBE,
  VM_DOUBLE_CUBE
} VimosVarType;


typedef struct _VIMOS_DPOINT_ {
  double x ;
  double y ;
  struct _VIMOS_DPOINT_ *prev;
  struct _VIMOS_DPOINT_ *next;
} VimosDpoint ;


/*--------------------------------------------------------------------------*/
/* The following structure stores pixel positions in an image               */

typedef struct _VIMOS_PIXEL_
{
  double x;
  double y;
  float  i;
  struct _VIMOS_PIXEL_ *prev;
  struct _VIMOS_PIXEL_ *next;
} VimosPixel;


VimosDpoint *newDpoint(int number);
void deleteDpoint(VimosDpoint *dpoint);

VimosPixel  *newPixel(int number);
void deletePixel(VimosPixel *pixel);

PIL_END_DECLS

#endif /* VM_TYPES_H */
