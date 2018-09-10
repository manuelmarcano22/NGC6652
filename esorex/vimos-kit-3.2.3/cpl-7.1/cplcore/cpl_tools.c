/*
 * This file is part of the ESO Common Pipeline Library
 * Copyright (C) 2001-2017 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cxmemory.h>
#include <cxmessages.h>
#include <cxutils.h>

#include "cpl_propertylist_impl.h"
#include "cpl_memory.h"
#include "cpl_tools.h"

#include "cpl_error_impl.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <signal.h>

/* Needed by strcmp(), strncmp(), strlen() */
#include <string.h>

#include "cpl_tools.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

#ifndef inline
#define inline /* inline */
#endif


#ifndef CPL_PIX_STACK_SIZE
#define CPL_PIX_STACK_SIZE 50
#endif

/* Swap macros */
#define CPL_INT_SWAP(a,b) { register const cpl_size t=(a); (a)=(b); (b)=t; }

#define CONCAT(a,b) a ## _ ## b
#define CONCAT2X(a,b) CONCAT(a,b)


/*----------------------------------------------------------------------------*/
/**
  @enum     cpl_fits_keytype
  @internal
  @brief    Possible FITS key types

  This enum stores all possible types for a FITS keyword. These determine
  the order of appearance in a header, they are a crucial point for
  DICB (ESO) compatibility. This classification is internal to this
  module.
 */
/*----------------------------------------------------------------------------*/
typedef enum _cpl_fits_keytype_ {
    cpl_fits_keytype_undef          =0,

    cpl_fits_keytype_top            =1,

    /* Mandatory keywords */
    /* All FITS files */
    cpl_fits_keytype_bitpix         =2,
    cpl_fits_keytype_naxis          =3,

    cpl_fits_keytype_naxis1         =11,
    cpl_fits_keytype_naxis2         =12,
    cpl_fits_keytype_naxis3         =13,
    cpl_fits_keytype_naxis4         =14,
    cpl_fits_keytype_naxisi         =20,
    /* Random groups only */
    cpl_fits_keytype_group          =30,
    /* Extensions */
    cpl_fits_keytype_pcount         =31,
    cpl_fits_keytype_gcount         =32,
    /* Main header */
    cpl_fits_keytype_extend         =33,
    /* Images */
    cpl_fits_keytype_bscale         =34,
    cpl_fits_keytype_bzero          =35,
    /* Tables */
    cpl_fits_keytype_tfields        =36,
    cpl_fits_keytype_tbcoli         =40,
    cpl_fits_keytype_tformi         =41,

    /* Other primary keywords */
    cpl_fits_keytype_primary        =100,

    /* HIERARCH ESO keywords ordered according to DICB */
    cpl_fits_keytype_hierarch_dpr   =200,
    cpl_fits_keytype_hierarch_obs   =201,
    cpl_fits_keytype_hierarch_tpl   =202,
    cpl_fits_keytype_hierarch_gen   =203,
    cpl_fits_keytype_hierarch_tel   =204,
    cpl_fits_keytype_hierarch_ins   =205,
    cpl_fits_keytype_hierarch_det   =206,
    cpl_fits_keytype_hierarch_log   =207,
    cpl_fits_keytype_hierarch_pro   =208,
    /* Other HIERARCH keywords */
    cpl_fits_keytype_hierarch       =300,

    /* HISTORY and COMMENT */
    cpl_fits_keytype_history        =400,
    cpl_fits_keytype_comment        =500,
    /* END */
    cpl_fits_keytype_end            =1000
} cpl_fits_keytype;

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static int cpl_fits_property_comparison(const cpl_property *, 
                                        const cpl_property *);
static cpl_fits_keytype cpl_fits_property_get_type(const char *) CPL_ATTR_PURE;

/*-----------------------------------------------------------------------------
                           Private variables
 -----------------------------------------------------------------------------*/

static cpl_flops flop_count = 0;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cpl_tools Utility functions
 *
 * This module provides various functions to apply simple operations on
 * data arrays (sorting, median computation).
 *
 * @par Synopsis:
 * @code
 *   #include "cpl_tools.h"
 * @endcode
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                            Function codes
 -----------------------------------------------------------------------------*/

#define CPL_TYPE_IS_TYPE_PIXEL
#define CPL_TYPE_IS_NUM
#define CPL_TYPE int
#define CPL_TYPE_NAME int
#define CPL_TOOLS_SORT_LT(x,y) ((x) < (y))
#include "cpl_tools_body.h"
#undef CPL_TYPE
#undef CPL_TYPE_NAME
#undef CPL_TOOLS_SORT_LT
#undef CPL_TYPE_IS_TYPE_PIXEL

#define CPL_TYPE_IS_NUM
#define CPL_TYPE long
#define CPL_TYPE_NAME long
#define CPL_TOOLS_SORT_LT(x,y) ((x) < (y))
#include "cpl_tools_body.h"
#undef CPL_TYPE
#undef CPL_TYPE_NAME
#undef CPL_TOOLS_SORT_LT

#define CPL_TYPE_IS_NUM
#define CPL_TYPE long long
#define CPL_TYPE_NAME long_long
#define CPL_TOOLS_SORT_LT(x,y) ((x) < (y))
#include "cpl_tools_body.h"
#undef CPL_TYPE
#undef CPL_TYPE_NAME
#undef CPL_TOOLS_SORT_LT

#define CPL_TYPE cpl_size
#define CPL_TYPE_NAME cplsize
#define CPL_TOOLS_SORT_LT(x,y) ((x) < (y))
#include "cpl_tools_body.h"
#undef CPL_TYPE
#undef CPL_TYPE_NAME
#undef CPL_TOOLS_SORT_LT

#define CPL_TYPE_IS_TYPE_PIXEL
#define CPL_TYPE float
#define CPL_TYPE_NAME float
#define CPL_TOOLS_SORT_LT(x,y) ((x) < (y))
#include "cpl_tools_body.h"
#undef CPL_TYPE
#undef CPL_TYPE_NAME
#undef CPL_TOOLS_SORT_LT

#define CPL_TYPE double
#define CPL_TYPE_NAME double
#define CPL_TOOLS_SORT_LT(x,y) ((x) < (y))
#include "cpl_tools_body.h"
#undef CPL_TYPE
#undef CPL_TYPE_NAME
#undef CPL_TOOLS_SORT_LT
#undef CPL_TYPE_IS_NUM
#undef CPL_TYPE_IS_TYPE_PIXEL

#define CPL_TYPE char *
#define CPL_TYPE_NAME string
#define CPL_TOOLS_SORT_LT(x,y) (((x) == NULL && (y) != NULL) ||  \
                                ((x) != NULL && (y) != NULL &&   \
                                 strcmp((x), (y)) < 0))
#include "cpl_tools_body.h"
#undef CPL_TYPE
#undef CPL_TYPE_NAME
#undef CPL_TOOLS_SORT_LT


/*----------------------------------------------------------------------------*/
/**
   @internal
   @brief Initializer with conditional allocation for ifalloc variable
   @param self  The ifalloc variable to be conditionally allocated
   @param size  The number of bytes needed, zero guarantees no heap allocation
   @see cpl_ifalloc_free
   @note Calling this function with a pointer that has previously been passed
         to this function with a non-zero size can lead to a memory leak.

   Example of usage:
    @code

    cpl_ifalloc mybuf;
    cpl_ifalloc_set(&mybuf, 0);

    if (!myerror) {
      double * mynumbers;

      cpl_ifalloc_set(&mybuf, mysize * sizeof(*mynumbers));

      mynumbers = (double*)cpl_ifalloc_get(&mybuf);
    }

    cpl_ifalloc_free(&mybuf);

    @endcode

 */
/*----------------------------------------------------------------------------*/
inline void cpl_ifalloc_set(cpl_ifalloc * self, size_t size)
{
    if (size > (size_t)CPL_IFALLOC_SZ) {
        self->data = self->heap = cpl_malloc(size);
    } else {
        self->data = (void*)(self->stack);
        self->heap = NULL;
    }
}

/*----------------------------------------------------------------------------*/
/**
   @internal
   @brief Get pointer to data area of ifalloc variable
   @param self  The ifalloc variable to accessed
   @see cpl_ifalloc_set
   @note The returned pointer becomes stale after a call to cpl_ifalloc_free()

 */
/*----------------------------------------------------------------------------*/
inline void * cpl_ifalloc_get(cpl_ifalloc * self)
{
    return self->data;
}

/*----------------------------------------------------------------------------*/
/**
   @internal
   @brief Free all memory associated with ifalloc variable
   @param self  The ifalloc variable to be deallocated
   @see cpl_ifalloc_set
 */
/*----------------------------------------------------------------------------*/
inline void cpl_ifalloc_free(cpl_ifalloc * self)
{
    cpl_free(self->heap);
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Find the FITS BITPIX value corresponding to the given CPL type.
  @param    type The CPL type.
  @return   The corresponding FITS BITPIX value, or zero on error.
  @note The input type may have its array/pointer bit set, which is ignored.

  The supported types (and their corresponding BITPIX) are:
  CPL_TYPE_UCHAR:  BYTE_IMG     (8)
  CPL_TYPE_SHORT:  SHORT_IMG   (16)
  CPL_TYPE_USHORT: USHORT_IMG  (20)
  CPL_TYPE_INT:    LONG_IMG    (32)
  CPL_TYPE_FLOAT:  FLOAT_IMG  (-32)
  CPL_TYPE_DOUBLE: DOUBLE_IMG (-64)

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_INVALID_TYPE if no BITPIX value corresponds to the type

 */
/*----------------------------------------------------------------------------*/
int cpl_tools_get_bpp(cpl_type type)
{
    int bpp;

    /* switch on type without POINTER and ARRAY bits */
    switch (type & ~CPL_TYPE_POINTER & ~CPL_TYPE_FLAG_ARRAY) {

    case CPL_TYPE_UCHAR:
        bpp = BYTE_IMG; /* 8 */
        break;

    case CPL_TYPE_SHORT:
        bpp = SHORT_IMG; /* 16 */
        break;

    case CPL_TYPE_USHORT:
        bpp = USHORT_IMG; /* 20 */
        break;

    case CPL_TYPE_INT:
        bpp = LONG_IMG; /* 32 */
        break;

    case CPL_TYPE_FLOAT:
        bpp = FLOAT_IMG; /* -32 */
        break;

    case CPL_TYPE_DOUBLE:
        bpp = DOUBLE_IMG; /* -64 */
        break;

    default:
        bpp = 0;
        (void)cpl_error_set_message_(CPL_ERROR_INVALID_TYPE,
                                     "type=0x%x", (int)type);
    }

    return bpp;
}


/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Determine if an integer is a power of 2.
  @param    p   Integer to check.
  @return   The corresponding power of 2 or -1 if p is not a power of 2
 */
/*----------------------------------------------------------------------------*/
cpl_size cpl_tools_is_power_of_2(cpl_size p)
{
    cpl_size ipow = 0;
    cpl_size done;

    if (p <= 0) return -1;

    /* Count until 1st non-zero bit is found */
    do {
        done = p & 1;
        p >>= 1;
        ipow++;
    } while (!done);

    /* The number is a power iff p is now zero */
    return p == 0 ? ipow-1 : -1;

}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief   Compute x to the power of y
  @param   x  The base
  @param   y  The power
  @return  x to the power of y
  @note    Iff y is negative and x is zero an actual division by zero will occur

  Apart from a possible difference in round-off the result equals pow(x,y).

 */
/*----------------------------------------------------------------------------*/
double cpl_tools_ipow(double x, int y)
{

    double result;
    double pow2 = x;
    int p = abs(y);

    /* Compute the result as a product of powers of 2 of x.
       - this method may produce (slightly) different results than pow(x,y) */

    /* Handle least significant bit in p here in order to avoid an unnecessary
       multiplication of pow2 - which could cause an over- or underflow */
    /* Also, 0^0 is 1, while any other power of 0 is 0 */
    result = p & 1 ? x : 1.0;

    while (p >>= 1) {
        pow2 *= pow2;
        /* Process least significant bit in p */
        if (p & 1) result *= pow2;
    }

    /* It is left to the caller to ensure that no division by zero occurs */
    return y < 0 ? 1.0/result : result;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Generate a valid FITS header for saving functions
  @param    self      The fitsfile to modify
  @param    to_add    The set of keywords to add to the minimal header
  @param    to_rm     The set of keys to remove
  @return   1 newly allocated valid header

  The passed file should contain a minimal header. 
  The propertylist is sorted (DICB) and added after the miniaml header.
  
  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if self is NULL
  - CPL_ERROR_ILLEGAL_INPUT if the propertylist cannot be written
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_fits_add_properties(fitsfile               * self,
                                       const cpl_propertylist * to_add,
                                       const char             * to_rm)
{
    cpl_error_code error = CPL_ERROR_NONE;

    if (to_add != NULL) {
        cpl_propertylist * out;
    
        /* Failure here would indicate a bug in CPL */
        cpl_ensure_code(self, CPL_ERROR_NULL_INPUT);

        if (to_rm != NULL) {
            /* Copy all but the black-listed properties */
            out = cpl_propertylist_new();
            if (cpl_propertylist_copy_property_regexp(out, to_add, to_rm, 1)) {
                error = cpl_error_set_where_();
            }
        } else {
            out = cpl_propertylist_duplicate(to_add);
        }

        if (!error) {
            /* Sort and write the propertylist to the file */
            if (cpl_propertylist_sort(out, cpl_fits_property_comparison) ||
                cpl_propertylist_to_fitsfile(self, out, NULL)) {

                error = cpl_error_set_where_();

            }

        }
        cpl_propertylist_delete(out);
    }

    return error;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Comparison function of two properties
  @param    p1      the first property
  @param    p2      the second property
  @return   -1 if p1<p2, 0 if p1==p2, +1 if p1>p2

  This function implements the sorting of the FITS header keywords to
  insure DICB compliant products
 */
/*----------------------------------------------------------------------------*/
static int cpl_fits_property_comparison(
        const cpl_property  *   p1,
        const cpl_property  *   p2)
{
    const char          *   key1;
    const char          *   key2;
    cpl_fits_keytype        kt1, kt2;
   
    /* Get the keys */
    key1 = cpl_property_get_name(p1);
    key2 = cpl_property_get_name(p2);
    
    /* Get the keywords types */
    kt1 = cpl_fits_property_get_type(key1);
    kt2 = cpl_fits_property_get_type(key2);
    
    /* Compare the types */
    if (kt1 > kt2) return 1;
    else if (kt1 < kt2) return -1;
    else return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Get the FITS keyword type
  @param    key     the keyword
  @return   the keyword type or -1
 */
/*----------------------------------------------------------------------------*/
static cpl_fits_keytype cpl_fits_property_get_type(const char * key)
{
    cpl_fits_keytype        kt;
   
    /* Check entries */
    if (key == NULL) return cpl_fits_keytype_undef;
    
    kt = cpl_fits_keytype_undef;
    if (!strcmp(key, "SIMPLE") || !strcmp(key, "XTENSION")) 
        kt = cpl_fits_keytype_top;
    else if (!strcmp(key, "END"))           kt = cpl_fits_keytype_end;
    else if (!strcmp(key, "BITPIX"))        kt = cpl_fits_keytype_bitpix;
    else if (!strcmp(key, "NAXIS"))         kt = cpl_fits_keytype_naxis;
    else if (!strcmp(key, "NAXIS1"))        kt = cpl_fits_keytype_naxis1;
    else if (!strcmp(key, "NAXIS2"))        kt = cpl_fits_keytype_naxis2;
    else if (!strcmp(key, "NAXIS3"))        kt = cpl_fits_keytype_naxis3;
    else if (!strcmp(key, "NAXIS4"))        kt = cpl_fits_keytype_naxis4;
    else if (!strncmp(key, "NAXIS", 5))     kt = cpl_fits_keytype_naxisi;
    else if (!strcmp(key, "GROUP"))         kt = cpl_fits_keytype_group;
    else if (!strcmp(key, "PCOUNT"))        kt = cpl_fits_keytype_pcount;
    else if (!strcmp(key, "GCOUNT"))        kt = cpl_fits_keytype_gcount;
    else if (!strcmp(key, "EXTEND"))        kt = cpl_fits_keytype_extend;
    else if (!strcmp(key, "BSCALE"))        kt = cpl_fits_keytype_bscale;
    else if (!strcmp(key, "BZERO"))         kt = cpl_fits_keytype_bzero;
    else if (!strcmp(key, "TFIELDS"))       kt = cpl_fits_keytype_tfields;
    else if (!strncmp(key, "TBCOL", 5))     kt = cpl_fits_keytype_tbcoli;
    else if (!strncmp(key, "TFORM", 5))     kt = cpl_fits_keytype_tformi;
    else if (!strncmp(key, "ESO DPR", 7))   kt = cpl_fits_keytype_hierarch_dpr;
    else if (!strncmp(key, "ESO OBS", 7))   kt = cpl_fits_keytype_hierarch_obs;
    else if (!strncmp(key, "ESO TPL", 7))   kt = cpl_fits_keytype_hierarch_tpl;
    else if (!strncmp(key, "ESO GEN", 7))   kt = cpl_fits_keytype_hierarch_gen;
    else if (!strncmp(key, "ESO TEL", 7))   kt = cpl_fits_keytype_hierarch_tel;
    else if (!strncmp(key, "ESO INS", 7))   kt = cpl_fits_keytype_hierarch_ins;
    else if (!strncmp(key, "ESO DET", 7))   kt = cpl_fits_keytype_hierarch_det;
    else if (!strncmp(key, "ESO LOG", 7))   kt = cpl_fits_keytype_hierarch_log;
    else if (!strncmp(key, "ESO PRO", 7))   kt = cpl_fits_keytype_hierarch_pro;
    else if (!strncmp(key, "ESO", 3))       kt = cpl_fits_keytype_hierarch;
    else if (!strcmp(key, "HISTORY"))       kt = cpl_fits_keytype_history;
    else if (!strcmp(key, "COMMENT"))       kt = cpl_fits_keytype_comment;
    else if ((int)strlen(key)<9)            kt = cpl_fits_keytype_primary;
    return kt;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Add to the FLOPS counter.
  @param flops  FLOPS to add
  @return   void
  @note If compiled with -DCPL_ADD_FLOPS this function will incur significant
        thread scheduling and CPU overhead, if not it will not be used.

 */
/*----------------------------------------------------------------------------*/
inline void cpl_tools_add_flops_(cpl_flops flops)
{
#ifdef _OPENMP
#pragma omp atomic
#endif
    flop_count += flops;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Get the FLOPS count.
  @return   The FLOPS count
  @note This function is intended to be used only by the CPL test modules.
  @note If compiled with -DCPL_ADD_FLOPS this function will always return zero.

 */
/*----------------------------------------------------------------------------*/
cpl_flops cpl_tools_get_flops(void)
{
  return flop_count;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Transform: xhat = x- mean(x)
  @param    x   The vector to be transformed
  @param    pm  On return, *pm is the mean of x
  @return   The created transformed vector
  @note No NULL-pointer check in this internal function

 */
/*----------------------------------------------------------------------------*/
cpl_vector * cpl_vector_transform_mean(const cpl_vector * x, double * pm)
{

    cpl_vector * xhat = cpl_vector_duplicate(x);

    *pm = cpl_vector_get_mean(xhat);
    cpl_vector_subtract_scalar(xhat, *pm);

    return xhat;

}


/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief   Count the number of distinct values in the vector
  @param   self       The vector to check, it will be sorted
  @param   pndistinct The number of distinct values
  @return  CPL_ERROR_NONE iff the check is sucessful

 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_vector_count_distinct(cpl_vector * self,
                                         cpl_size * pndistinct)
{

    if (pndistinct == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    } else if (cpl_vector_sort(self, CPL_SORT_ASCENDING)) {
        return cpl_error_set_where_();
    } else {

        const double * dx = cpl_vector_get_data_const(self);
        double         xmin = dx[0];
        const cpl_size n = cpl_vector_get_size(self);
        cpl_size       i, j;

        for (j = i = 1; i < n; i++) {
            if (dx[i] > xmin) {
                xmin = dx[i];
                j++;
            }
        }

        *pndistinct = j;
    }

    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief   Ensure that the vector has the required number of distinct values
  @param   self      The vector to check
  @param   ndistinct The (positive) number of distinct values to require
  @return  CPL_ERROR_NONE iff the check is sucessful
  @see cpl_vector_count_distinct()

 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_vector_ensure_distinct(const cpl_vector * self,
                                          cpl_size ndistinct)
{

    if (ndistinct > 1) {
        const cpl_size n = cpl_vector_get_size(self);

        if (n < 1) {
            return cpl_error_set_where_();
        } else if (ndistinct > n) {
            return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND,
                                          "A %" CPL_SIZE_FORMAT "-element "
                                          "vector cannot have at least %"
                                          CPL_SIZE_FORMAT " distinct values",
                                          n, ndistinct);
        } else {
            cpl_vector * tmp = cpl_vector_duplicate(self);
            cpl_size     idistinct = 0;

            (void)cpl_vector_count_distinct(tmp, &idistinct);

            cpl_vector_delete(tmp);

            if (idistinct < ndistinct) {
                return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%"
                                              CPL_SIZE_FORMAT "-element vector "
                                              "must have at least %"
                                              CPL_SIZE_FORMAT " (not %"
                                              CPL_SIZE_FORMAT ") distinct "
                                              " values", n, ndistinct,
                                              idistinct);
            }
        }
    } else if (ndistinct < 1) {
        return cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
    }

    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Copy a rectangular memory area from one place to another
  @param    self        Preallocated location to copy to
  @param    src         Location to copy from
  @param    size        Size of each element [bytes]
  @param    nx          Number of elements in x direction in source
  @param    ny          Number of elements in y direction in source
  @param    llx         Lower left x position (starting with 1)
  @param    lly         Lower left y position (starting with 1)
  @param    urx         Upper right x position
  @param    ury         Upper right y position
  @return  CPL_ERROR_NONE if OK, otherwise the relevant #_cpl_error_code_.
  @note self must have place for (urx-llx+1) * (ury-lly+1) elements of size size
  @note No NULL-pointer check in this internal function
  @see mempcy()

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_ILLEGAL_INPUT if the window coordinates are not valid
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_tools_copy_window(void *self, const void *src, size_t size,
                                     cpl_size nx,  cpl_size ny,
                                     cpl_size llx, cpl_size lly,
                                     cpl_size urx, cpl_size ury)
{

    /* FIXME: Need to do pointer arithmetic */
    const char * csrc  = (const char*)src;

    cpl_ensure_code(size != 0,  CPL_ERROR_ILLEGAL_INPUT);

    cpl_ensure_code(llx <= urx, CPL_ERROR_ILLEGAL_INPUT);

    cpl_ensure_code(lly > 0,    CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(lly <= ury, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(ury <= ny,  CPL_ERROR_ILLEGAL_INPUT);


    /* assert(sizeof(char) == 1); */


    if (llx == 1 && urx == nx) {
        /* Number of bytes in one row of self */
        const size_t rowsize = size * (size_t)nx;

        (void)memcpy(self, csrc + rowsize * (size_t)(lly-1),
                     rowsize * (size_t)(ury-lly+1));
    } else {
        /* FIXME: Need to do pointer arithmetic */
        char       * cself   = (char*)self;
        /* Number of bytes in one row of self */
        const size_t rowsize = size * (size_t)(urx-llx+1);

        /* Byte just after last byte to write to self */
        const char * cstop   = cself + rowsize * (size_t)(ury-lly+1);

        cpl_ensure_code(llx >  0,   CPL_ERROR_ILLEGAL_INPUT);
        cpl_ensure_code(urx <= nx,  CPL_ERROR_ILLEGAL_INPUT);

        /* Point to first byte to read */
        csrc += size * (size_t)((llx - 1) + (lly-1) * nx);

        for (; cself < cstop; cself += rowsize, csrc += size * (size_t)nx) {
            (void)memcpy(cself, csrc, rowsize);
        }
    }

    return CPL_ERROR_NONE;
} 

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Shift a rectangular memory area
  @param    self    Location of 1st element
  @param    size    Size of each element [bytes]
  @param    nx      Number of elements in x direction in source
  @param    ny      Number of elements in y direction in source
  @param    null    Value to insert into 'empty' zone'
  @param    dx      Shift in X
  @param    dy      Shift in Y
  @return   the #_cpl_error_code_ or CPL_ERROR_NONE
  @note No NULL-pointer check in this internal function
  @see cpl_mask_shift()

  The 'empty zone' in the shifted rectangle is filled with size null's.
  The shift values have to be valid:
  -nx < x_shift < nx and -ny < y_shift < ny
  
  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if in is NULL
  - CPL_ERROR_ILLEGAL_INPUT if the offsets are too big, or size is zero, or
    nx or ny non-positive
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_tools_shift_window(void * self, size_t size,
                                      cpl_size nx, cpl_size ny, int null,
                                      cpl_size dx, cpl_size dy)
{

    /* Need to do pointer arithmetic */
    char * myself = (char*)self;
    /* Number of elements to be set to null */
    const size_t sset  = (size_t)(labs((long int)dy) * nx) * size;
    /* Remainder of buffer will be moved */
    const size_t smove = (size_t)(ny * nx) * size - sset;
    char * pdest = dy < 0 ? myself : myself + sset;
    char * psrc  = dy > 0 ? myself : myself + sset; /* const */
    char * pset  = dy > 0 ? myself : myself + smove;

    cpl_ensure_code(size >  0,    CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(nx   >  0,    CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(ny   >  0,    CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code( dx  <  nx,   CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(-dx  <  nx,   CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code( dy  <  ny,   CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(-dy  <  ny,   CPL_ERROR_ILLEGAL_INPUT);


    if (dx != 0) {
        /* Have to shift one row at a time */

        /* Size of a row */
        const size_t rsize = (size_t)nx * size;
        /* Number of elements to be set to null */
        const size_t srset = (size_t)labs((long int)dx) * size;
        /* Remainder of buffer will be moved */
        const size_t srmove = rsize - srset;

        char * prdest = dx < 0 ? pdest : pdest + srset;
        char * prsrc  = dx > 0 ? psrc  : psrc  + srset;
        char * prset  = dx > 0 ? pdest : pdest + srmove;
        /* source and dest overlap when dy is zero */
        void * (*myshift)(void *, const void *, size_t) = dy ? memcpy : memmove;
        cpl_size j;

        if (dy <= 0) {
            for (j = -dy; j < ny; j++,
                 prdest += rsize, prsrc += rsize, prset += rsize) {
                (void)myshift(prdest, prsrc, srmove);
                (void)memset(prset, null, srset);
            }
        } else {
            /* With a positive dy the rows must be shifted in reverse
               order */
            prdest += smove - rsize;
            prsrc  += smove - rsize;
            prset  += smove - rsize;

            for (j = dy; j < ny; j++,
                 prdest -= rsize, prsrc -= rsize, prset -= rsize) {
                (void)memcpy(prdest, prsrc, srmove);
                (void)memset(prset, null, srset);
            }
        }
    } else if (dy != 0) {
        /* Shift all rows at once */
        (void)memmove(pdest, psrc, smove);
    }

    /* Set all of the 'empty' rows */
    if (dy != 0) (void)memset(pset, null, sset);

    return CPL_ERROR_NONE;
}

/**@}*/
