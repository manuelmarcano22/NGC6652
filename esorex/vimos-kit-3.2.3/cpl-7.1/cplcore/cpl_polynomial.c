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

#include "cpl_polynomial_impl.h"

#include "cpl_tools.h"
#include "cpl_error_impl.h"
#include "cpl_matrix_impl.h"
#include "cpl_memory.h"

#include <cxmessages.h>

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <assert.h>

/*-----------------------------------------------------------------------------
                             Macro definitions
 -----------------------------------------------------------------------------*/

#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
  /* FIXME: Not ready yet :-( */
  #define CPL_POLYNOMIAL_EVAL_FUNCTION cpl_polynomial_eval_horner
#else
  #define CPL_POLYNOMIAL_EVAL_FUNCTION cpl_polynomial_eval_bf
#endif

/* Maximum number of Newton-Raphson iterations */
#ifndef CPL_NR_MAXITE
#define CPL_NR_MAXITE 100
#endif

#define CPL_POLYNOMIAL_CMP                                              \
                /* Verify that it differs within tolerance */           \
                if (fabs(p1->c[i] - p2->c[j]) <= tol) {                 \
                /* Verify that the powers match */                      \
                    for (dim=0; dim < p1->dim; dim++)                   \
                        if (p1->pow[p1->dim * i + dim]                  \
                         != p2->pow[p1->dim * j + dim]) break;          \
                    if (dim == p1->dim) break;   /* - found it */       \
                }

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cpl_polynomial Polynomials
 *
 * This module provides functions to handle uni- and multivariate polynomials.
 *
 * Univariate (1D) polynomials use the Horner rule for evaluation, while
 * multivariate polynomials are evaluated simply as the sum of each term.
 *
 * This means that of the two polynomials
 * @verbatim
 * P1(x) = p0 + p1.x + p4.x^2
 * @endverbatim
   and
 * @verbatim
 * P2(x,y) = p0 + p1.x + p2.y + p3.x.y + p4.x^2 + p5.y^2
 * @endverbatim
 * P1(x) may evaluate to more accurate results than P2(x,0),
 * especially around the roots.
 *
 * Note that a polynomial like P3(z) = p0 + p1.z + p2.z^2 + p3.z^3, z=x^4
 * is preferable to p4(x) = p0 + p1.x^4 + p2.x^8 + p3.x^12.
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                                   Type definition
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
                  Polynomial object.

  If you want to store the following polynomial:
  @verbatim
  p(x,y) = p0 + p1.x + p2.y + p3.x.y + p4.x^2 + p5.y^2
  @endverbatim
  You would internally have:
  @verbatim
  nc  = 6 (from 0 to 5 incl.)
  dim = 2 (x and y)
  c[]  pow[] contains (pow[] is stored row-wise):
  p0   0 0
  p1   1 0
  p2   0 1
  p3   1 1
  p4   2 0
  p5   0 2
  @endverbatim
  i.e. the i'th coefficient is multiplied by the variable in dimension j
  lifted to this power:
  pow[dim * i + j]

  1-dimensional polynomials are stored differently:
  pow[] is not used.  dim == 1.
  Instead all coefficients (including those with zero-value) are stored
  in order of increasing power - in c.
  Thus the constant coefficient is stored in c[0] while the leading coefficient
  is stored in c[nc-1].
  This storage scheme allows for the usage of the Horner rule.

 */
/*----------------------------------------------------------------------------*/

struct _cpl_polynomial_ {

    /* The dimension of the polynomial */
    cpl_size    dim;

    /* if dim == 1: The number of coefficients in c (degree + 1).
       if dim >  1: The number of (non-zero) coefficients in c. */
    cpl_size    nc;

    /* The coefficients. If dim == 1: c0, c_1, ..., c_{nc - 1} */
    double  *   c;

    /* If dim == 1: Not used
       If dim > 1: The nc * dim powers of the variables x1, x2, ..., xdim */
    cpl_size *  pow;

    /* If dim == 1: Not used
       If dim > 1: Temporary storage to speed up evaluation
       max_degree[i] is the largest power used for dimension i.
       eval_pow[i] holds space for 1+max_degree[i] double's. 
       max_degree_alldims is the largest power counting all dims 
       (i.e. 5 for x2y3) */
    cpl_size * max_degree;
    cpl_size   max_degree_alldims;
    double ** eval_pow;
    
#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
    /* Cache variables for Horner evaluation */
    double *   added_horner;
    int    *   added_exist_iterations;
    int        added_exist_index;
    int        added_iter;
#endif
};

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_error_code cpl_polynomial_delete_coeff(cpl_polynomial *,
                                                  const cpl_size *);

static
void cpl_matrix_fill_normal_vandermonde(cpl_matrix *, cpl_matrix *,
                                        const cpl_vector *, cpl_boolean,
                                        cpl_size, const cpl_vector *)
    CPL_ATTR_NONNULL;

static
cpl_error_code cpl_polynomial_fit_2d(cpl_polynomial *, const cpl_bivector *,
                                     const cpl_vector *, cpl_boolean,
                                     const cpl_size *, double *);

static int cpl_vector_is_eqdist(const cpl_vector *);

#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
static double cpl_polynomial_eval_horner(const cpl_polynomial * p,
                                         const cpl_vector * x);

static
double cpl_polynomial_register_horner(cpl_polynomial * p,
                                      int              ireg,
                                      int            * step,
                                      int            * powers_tmp,
                                      const double   * x);

static
void cpl_polynomial_arrange_coeffs_horner(cpl_polynomial  * p,
                                          int               ireg,
                                          int             * step,
                                          int             * powers_tmp);

static
int cpl_polynomial_coeff_index(const cpl_polynomial *,
                               const int *);
#else
static double cpl_polynomial_eval_bf(const cpl_polynomial *, 
                                     const cpl_vector *) CPL_ATTR_NONNULL;
#endif

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a new cpl_polynomial
  @param    dim    The positive polynomial dimension (number of variables)
  @return   1 newly allocated cpl_polynomial, or NULL on error

  The returned object must be deallocated using cpl_polynomial_delete().

  A newly created polynomial has degree 0 and evaluates as 0.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_ILLEGAL_INPUT if dim is negative or zero
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cpl_polynomial_new(cpl_size dim)
{
    cpl_polynomial * p;

    /* Test input */
    cpl_ensure(dim > 0, CPL_ERROR_ILLEGAL_INPUT, NULL);

    /* Allocate struct - and set nc == 0 */
    p = cpl_calloc(1, sizeof(cpl_polynomial));
    p->dim = dim;
    p->nc  = 0;

    return p;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Delete a cpl_polynomial
  @param    p    cpl_polynomial to delete
  @return   void

   The function deallocates the memory used by the polynomial @em p.
   If @em p is @c NULL, nothing is done, and no error is set.

 */
/*----------------------------------------------------------------------------*/
void cpl_polynomial_delete(cpl_polynomial * p)
{

    if (p == NULL) return;

    if (p->nc > 0) {
        cpl_free(p->c);
        if (p->dim > 1) {
            while (p->dim--) {
                cpl_free(p->eval_pow[p->dim]);
            }
            cpl_free(p->pow);
            cpl_free(p->eval_pow);
            cpl_free(p->max_degree);
#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
            cpl_free(p->added_horner);
            cpl_free(p->added_exist_iterations);
#endif
        }
    }

    cpl_free(p);

    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Dump a cpl_polynomial as ASCII to a stream
  @param    p       Input cpl_polynomial to dump
  @param    stream  Output stream, accepts @c stdout or @c stderr
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_

  Each coefficient is preceded by its integer power(s) and
  written on a single line.
  If the polynomial has non-zero coefficients, only those are printed,
  otherwise the (zero-valued) constant term is printed.

  Comment lines start with the hash character.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_FILE_IO if the write operation fails
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_dump(const cpl_polynomial * p,
                                   FILE                 * stream)
{
    const cpl_error_code err = CPL_ERROR_FILE_IO;
    int i, dim;


    cpl_ensure_code(stream != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(p      != NULL, CPL_ERROR_NULL_INPUT);

    cpl_ensure_code( fprintf(stream,
        "#----- %" CPL_SIZE_FORMAT " dimensional polynomial -----\n", p->dim)
                     > 0, err);

    for (i=0; i < p->dim; i++)
        cpl_ensure_code(fprintf(stream, "%d.dim.power  ", i+1) > 0, err);

    cpl_ensure_code(fprintf(stream, "coefficient\n") > 0,  err);

    if (p->nc == 0) {
        for (dim=0; dim < p->dim; dim++)
            cpl_ensure_code(fprintf(stream, "  %5d      ", 0) > 0, err);

        cpl_ensure_code(fprintf(stream, "0\n") > 0,  err);
    }

    for (i=0; i < p->nc; i++) {
        if (p->dim != 1) {
            for (dim=0; dim < p->dim; dim++)
                cpl_ensure_code(fprintf(stream, "  %5" CPL_SIZE_FORMAT "      ",
                                        p->pow[p->dim * i + dim]) > 0, err);
        } else if (p->c[i] != 0.0) {
            /* Only non-zero coefficients */
            cpl_ensure_code(fprintf(stream, "  %5d      ", i) > 0,  err);
        } else {
            continue;
        }

        cpl_ensure_code(fprintf(stream, "%g\n", p->c[i]) > 0, err);

    }

    cpl_ensure_code(
        fprintf(stream, "#------------------------------------\n") > 0, err);

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    This function duplicates an existing polynomial
  @param    p   The input cpl_polynomial
  @return   A newly allocated cpl_polynomial or NULL on error

  Notice that the returned object is a newly allocated cpl_polynomial that
  must be deallocated by the caller using cpl_polynomial_delete().

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
*/
/*----------------------------------------------------------------------------*/
cpl_polynomial * cpl_polynomial_duplicate(const cpl_polynomial * p)
{
    cpl_polynomial * out;


    /* Test inputs */
    cpl_ensure(p != NULL,  CPL_ERROR_NULL_INPUT, NULL);

    if (p->nc == 0) return cpl_polynomial_new(p->dim);

    /* Create a new cpl_polynomial and copy the non-pointer members */
    out = cpl_malloc(sizeof(cpl_polynomial));
    memcpy(out, p, sizeof(cpl_polynomial));

    /* Copy the coefficients */
    out->c = cpl_malloc((size_t)out->nc * sizeof(double));
    memcpy(out->c, p->c, (size_t)out->nc * sizeof(double));

    if (out->dim > 1) {
        /* Copy the powers - and eval-workspace */
        int dim;

        out->pow = cpl_malloc((size_t)out->nc * (size_t)out->dim
                              * sizeof(cpl_size));
        memcpy(out->pow, p->pow, (size_t)out->nc * (size_t)out->dim
               * sizeof(cpl_size));

        out->max_degree = cpl_malloc((size_t)out->dim * sizeof(cpl_size));
        memcpy(out->max_degree, p->max_degree,
               (size_t)out->dim * sizeof(cpl_size));

        /* out->max_degree_alldims = p->max_degree_alldims;  Unnecesary.. */

        out->eval_pow = cpl_malloc((size_t)out->dim * sizeof(double*));

        for (dim=0; dim < out->dim; dim++) {

            out->eval_pow[dim] = cpl_malloc((size_t)(1+out->max_degree[dim])
                    * sizeof(double));
            out->eval_pow[dim][0] = 1.0; /* Always 1 */

        }
#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
        out->added_horner = 
            cpl_malloc(out->nc * sizeof(*out->added_horner));
        out->added_exist_iterations = 
            cpl_malloc(out->nc * sizeof(*out->added_exist_iterations));
#endif
    }


    return out;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    This function copies contents of a polynomial into another one
  @param    out  Pre-allocated output cpl_polynomial
  @param    in   Input cpl_polynomial
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_

  in and out must point to different polynomials.

  If out already contains coefficients, they are overwritten.

  This is the only function that can modify the dimension of a polynomial.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INCOMPATIBLE_INPUT if in and out point to the same polynomial
*/
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_copy(cpl_polynomial * out,
                                   const cpl_polynomial * in)
{


    /* Check input */
    cpl_ensure_code(out != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(in  != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(in  != out,  CPL_ERROR_INCOMPATIBLE_INPUT);


    /* Verify: Should this be done only when in->pow == out->pow ?? */ 
    if (out->dim == in->dim && in->nc > 0 && out->nc > 0) {
        /* The polynomials have equal dimension
          - realloc() if needed and just copy all buffers */

        if (in->dim > 1) {
            /* Copy the powers */
            int dim;

            if (out->nc != in->nc) out->pow =
                cpl_realloc(out->pow, (size_t)in->dim * (size_t)in->nc
                            * sizeof(cpl_size));
            memcpy(out->pow, in->pow, (size_t)in->dim * (size_t)in->nc
                   * sizeof(cpl_size));

            memcpy(out->max_degree, in->max_degree,
                   (size_t)in->dim * sizeof(cpl_size));

            out->max_degree_alldims = in->max_degree_alldims; 

            for (dim=0; dim < in->dim; dim++) {
                out->eval_pow[dim] =
                    cpl_realloc(out->eval_pow[dim],
                            (size_t)(1+out->max_degree[dim]) * sizeof(double));
            }
#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
            cpl_free(out->added_horner);
            cpl_free(out->added_exist_iterations);
            out->added_horner           = 
                cpl_malloc((size_t)out->nc * sizeof(*out->added_horner));
            out->added_exist_iterations = 
                cpl_malloc((size_t)out->nc
                           * sizeof(*out->added_exist_iterations));
#endif
        }

        /* Resize the coefficient buffer  */
        if (out->nc != in->nc) {
            out->nc = in->nc;
            out->c = cpl_realloc(out->c, (size_t)in->nc * sizeof(double));
        }

        memcpy(out->c, in->c, (size_t)in->nc * sizeof(double));


    } else {
        /* When the dimensions are different,
           or if one of the polynomials have no coefficients,
           avoiding the malloc()/free()'s is too complicated
           - instead just duplicate the source and swap it
           into place */
        cpl_polynomial * tmp  = cpl_polynomial_duplicate(in);
        void           * swap = cpl_malloc(sizeof(cpl_polynomial));


        memcpy(swap, out,  sizeof(cpl_polynomial));
        memcpy(out,  tmp,  sizeof(cpl_polynomial));
        memcpy(tmp,  swap, sizeof(cpl_polynomial));

        /* Delete the buffers originally used by out - and the swap space */
        cpl_polynomial_delete(tmp);
        cpl_free(swap);

    }

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compare the coefficients of two polynomials
  @param    p1   the 1st polynomial
  @param    p2   the 2nd polynomial
  @param    tol  The absolute (non-negative) tolerance
  @return   0 when equal, positive when different, negative on error.

  The two polynomials are considered equal iff they have identical
  dimensions and the absolute difference between their coefficients
  does not exceed the tolerance.

  This means that the following pair of polynomials per definition are
  considered different:
  P1(x1,x2) = 3*x1 different from P2(x1) = 3*x1.

  If all parameters are valid and p1 and p2 point to the same polynomial the
  functions returns 0.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_ILLEGAL_INPUT if tol is negative
 */
/*----------------------------------------------------------------------------*/
int cpl_polynomial_compare(const cpl_polynomial * p1, const cpl_polynomial * p2,
                           double tol)
{

    cpl_size i;

    /* Test parameters */
    cpl_ensure(p1 != NULL, CPL_ERROR_NULL_INPUT, -1);
    cpl_ensure(p2 != NULL, CPL_ERROR_NULL_INPUT, -2);

    cpl_ensure(tol >= 0,   CPL_ERROR_ILLEGAL_INPUT, -5);


    if (p1 == p2) return 0;

    if (p1->dim != p2->dim) return 1;

    if (p1->dim == 1) {
        /* If the degrees differ, the leading coefficient(s) must be
           smaller than tol */
        for (i = p1->nc; i > p2->nc; i--) if (fabs(p1->c[i]) > tol) return 2;
        for (i = p2->nc; i > p1->nc; i--) if (fabs(p2->c[i]) > tol) return 3;

        /* At this point i is equal to min(p1->nc, p2->nc) */
        while (i--) if (fabs(p1->c[i]-p2->c[i]) > tol) return 4;

    } else {
        /* For each 'non-zero' coefficient in p1, the same one must be in p2
           and vice-versa */
        cpl_size * found2 = cpl_calloc((size_t)p2->nc, sizeof(cpl_size));
        cpl_size dim;
        cpl_size j;

        /* Verify that each coefficient in p1 is either below tol or
           present in p2. Remember which of the coefficients in p2 were
           found */
        for (i = 0; i < p1->nc; i++) {
            if (fabs(p1->c[i]) <= tol) continue;
            for (j = 0; j < p2->nc; j++) {
                if (found2[j]) continue; /* Use a coefficient just once */
                CPL_POLYNOMIAL_CMP;
            }
            /* Verify that one was found */
            if (j == p2->nc) break;
            found2[j] = TRUE;
        }
        if (i < p1->nc) {
            cpl_free(found2);
            return 5;
        }
        for (j = 0; j < p2->nc; j++) {
            if (found2[j] || fabs(p2->c[j]) <= tol) continue;
            for (i=0; i < p1->nc; i++) {
                CPL_POLYNOMIAL_CMP;
            }
            /* Verify that one was found */
            if (i == p1->nc) break;
        }
        cpl_free(found2);

        if (j < p2->nc) return 6;
    }

    return 0;

}

/*----------------------------------------------------------------------------*/
/**
  @brief    The degree of the polynomial
  @param    p   the polynomial
  @return   The degree or negative on error

  The degree is the highest sum of exponents (with a non-zero coefficient).

  If there are no non-zero coefficients the degree is zero.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
 */
/*----------------------------------------------------------------------------*/
cpl_size cpl_polynomial_get_degree(const cpl_polynomial * p)
{
    cpl_size degree;


    /* Test parameters */
    cpl_ensure(p != NULL, CPL_ERROR_NULL_INPUT, -1);

    if (p->nc == 0) return 0;

    if (p->dim == 1) {
        return p->nc-1;
    } else {
        cpl_size dim;
        cpl_size sum;
        cpl_size i;

        degree = -1;
        for (i = 0; i < p->nc; i++) {
            if (p->c[i] == 0) continue;
            sum = 0;
            for (dim = 0; dim < p->dim; dim++) sum += p->pow[p->dim * i + dim];
            if (sum > degree) degree = sum;
        }
    }

    /* User may have reset leading coefficient(s) to zero */
    if (degree < 0) return 0;

    return degree;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    The dimension of the polynomial
  @param    p   the polynomial
  @return   The dimension or negative on error

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
 */
/*----------------------------------------------------------------------------*/
cpl_size cpl_polynomial_get_dimension(const cpl_polynomial * p)
{


    /* Test entries */
    cpl_ensure(p != NULL, CPL_ERROR_NULL_INPUT, -1);

    return p->dim;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get a coefficient of the polynomial
  @param    in      the input polynomial
  @param    pows    the powers of the different variables
  @return   The coefficient or undefined on error

  The array pows must have the size of the polynomial dimension and have
  non-negative elements.

  It is allowed to specify a (set of) power(s) for which no coefficient
  has previously been set. In this case the function returns zero.

  Example of usage:
    @code

    const cpl_size power       = 3;
    double         coefficient = cpl_polynomial_get_coeff(poly1d, &power);

    @endcode

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_ILLEGAL_INPUT if pows contains negative values
 */
/*----------------------------------------------------------------------------*/
double cpl_polynomial_get_coeff(const cpl_polynomial * in,
                                const cpl_size       * pows)
{
    double coef = 0.0; /* Default value */
    cpl_size dim;


    /* Test entries */
    cpl_ensure(in   != NULL, CPL_ERROR_NULL_INPUT, coef);
    cpl_ensure(pows != NULL, CPL_ERROR_NULL_INPUT, coef);

    for (dim = 0;  dim < in->dim; dim++) {
        if (pows[dim] < 0) break;
    }

    if (dim < in->dim) {
        (void)cpl_error_set_message_(CPL_ERROR_ILLEGAL_INPUT, "Dimension %"
                                     CPL_SIZE_FORMAT " of %" CPL_SIZE_FORMAT
                                     " has negative power %" CPL_SIZE_FORMAT,
                                     dim+1, in->dim, pows[dim]);
    } else if (in->dim == 1) {

        if (pows[0] < in->nc) coef = in->c[pows[0]];

    } else {
        cpl_size i;
        for (i=0; i < in->nc; i++) {
            if (!memcmp(in->pow + (size_t)(in->dim * i), pows,
                        (size_t)in->dim * sizeof(*pows)))
                break; /* Found the right combination of powers */
        }
        if (i < in->nc) coef = in->c[i];
    }

    return coef;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Set a coefficient of the polynomial
  @param    in      the input polynomial
  @param    pows    the powers of the different variables
  @param    c       the coefficient
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_

  The array pows must have the size of the polynomial dimension and have
  non-negative elements.

  If the coefficient is already there, it is overwritten, if not, a new
  coefficient is added to the polynomial. This may cause the degree of the
  polynomial to be increased.

  Setting the coefficient of x1^4 * x3^2 in the 4-dimensional polynomial
  poly4d to 12.3 would be performed by:

    @code

    const cpl_size pows[] = {4, 0, 2, 0};
    cpl_error_code error  = cpl_polynomial_set_coeff(poly4d, pows, 12.3);

    @endcode

  Setting the coefficient of x^3 in the 1-dimensional polynomial poly1d to
  12.3 would be performed by:
    @code

    const cpl_size power = 3;
    cpl_error_code error = cpl_polynomial_set_coeff(poly1d, &power, 12.3);

    @endcode


  For efficiency reasons the coefficients of a 1D-polynomial are best
  inserted with that of the highest power first.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_ILLEGAL_INPUT if pows contains negative values
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_set_coeff(cpl_polynomial * in,
                                        const cpl_size * pows,
                                        double           c)
{
    cpl_size     newnc;
    cpl_size     dim;
    cpl_size     ind  = -1; /* Assume coefficient is absent */


    /* Test entries */
    cpl_ensure_code(in,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(pows, CPL_ERROR_NULL_INPUT);
    for (dim=0; dim < in->dim; dim++)
        cpl_ensure_code(pows[dim] >= 0, CPL_ERROR_ILLEGAL_INPUT);

    /* Handle zero-valued coefficients as a special case */
    /* FIXME: 0 is a valid value */
    if (c == 0.0) return cpl_polynomial_delete_coeff(in, pows); 

    /* Find the coeff */
    if (in->dim == 1) {
        if (pows[0] < in->nc) ind = pows[0];
    } else {
        cpl_size i;
        for (i=0; i < in->nc; i++) {
            if (!memcmp(in->pow + (size_t)in->dim * (size_t)i, pows,
                        (size_t)in->dim * sizeof(cpl_size)))
                break; /* Found the right combination of powers */
        }
        if (i < in->nc) ind = i;
    }

    /* The coefficient already exists : replace it and return */
    if (ind >= 0) {

        in->c[ind] = c;

        return CPL_ERROR_NONE;
    }


    /* The coefficient is a new one */
    /* Extend the arrays and update counter */

    /* New size. 1D-polynomials can grow with more than one coefficient */
    newnc = in->dim == 1 ? 1 + pows[0] : 1 + in->nc;

    in->c = in->nc ? cpl_realloc(in->c, (size_t)newnc * sizeof(double))
                   : cpl_malloc((size_t)newnc * sizeof(double));

    if (in->dim == 1) {
        /* Initialize coefficients between existing and new one to zero */
        if (newnc - in->nc > 1)
            memset(in->c + (size_t)in->nc, 0, (size_t)(newnc - in->nc - 1)
                   * sizeof(double));
    } else {
        /* Extend and copy the power arrays */

        in->pow = newnc == 1
            ? cpl_malloc((size_t)in->dim * (size_t)newnc * sizeof(cpl_size))
            : cpl_realloc(in->pow, (size_t)in->dim * (size_t)newnc
                          * sizeof(cpl_size));

        memcpy(in->pow + (size_t)in->dim * (size_t)in->nc, pows,
               (size_t)in->dim * sizeof(cpl_size));

        if (in->nc == 0) {

            in->max_degree = cpl_malloc((size_t)in->dim * sizeof(cpl_size));
            in->eval_pow   = cpl_malloc((size_t)in->dim * sizeof(double*));
            in->max_degree_alldims = 0; 

            for (dim=0; dim < in->dim; dim++) {
                in->max_degree[dim] = pows[dim];
                in->max_degree_alldims += in->max_degree[dim];
                in->eval_pow[dim] = cpl_malloc((size_t)(1+in->max_degree[dim])
                                               * sizeof(double));
                in->eval_pow[dim][0] = 1.0; /* Always 1 */
            }
#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
            /* FIXME: This is going to be deleted a few lines below, anyway */
            in->added_horner           = 
                cpl_malloc(in->nc * sizeof(*in->added_horner));
            in->added_exist_iterations = 
                cpl_malloc(in->nc * sizeof(*in->added_exist_iterations));
#endif            

        } else {
            for (dim = 0; dim < in->dim; dim++) {
                /* Verify */
                if (pows[dim] <= in->max_degree[dim]) continue;
                in->max_degree_alldims += pows[dim] - in->max_degree[dim];
                in->max_degree[dim] = pows[dim];
                in->eval_pow[dim] =
                    cpl_realloc(in->eval_pow[dim],
                              (size_t)(1+in->max_degree[dim]) * sizeof(double));

            }
        }
    }

    /* Set the new number of coefficients */
    in->nc = newnc;

    /* Set the new coefficient */
    in->c[in->nc-1] = c;

#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
    /* Arrange the added terms for Horner scheme */
    if(in->dim != 1)
    {
        /* Variables for Horner speed-up */
        int     * horner_step = cpl_malloc(sizeof(int) * dim);
        int     * powers_tmp  = cpl_malloc(sizeof(int) * dim);
        const int ireg        = 0;

        for(int istep = 0; istep < dim; ++istep) {
            horner_step[istep] = in->max_degree_alldims;
        }
        cpl_free(in->added_horner);
        cpl_free(in->added_exist_iterations);
        in->added_horner           = 
            cpl_malloc(in->nc * sizeof(*in->added_horner));
        in->added_exist_iterations = 
            cpl_malloc(in->nc * sizeof(*in->added_exist_iterations));
        
        in->added_exist_index = 0;
        in->added_iter = 0;
        cpl_polynomial_arrange_coeffs_horner
            (in, ireg, horner_step, powers_tmp);
        cpl_free(horner_step);
        cpl_free(powers_tmp);
    }
#endif
        
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Evaluate the polynomial at the given point 
  @param    p  The polynomial
  @param    x  Point of evaluation
  @return   The computed value or undefined on error.

  The length of x must be the polynomial dimension.

  A polynomial with no non-zero coefficients evaluates as 0.

  If the classical evaluation method is used, the computational cost is:  

  For 1-dimensional polynomials the call requires 2n FLOPs
  where n+1 is the number of coefficients in p, see also
  cpl_polynomial_eval_1d().

  For multivariate polynomials the call requires
  n*(1+dim) + d_1 + d_2 + ... + d_dim FLOPs,
  where dim is the dimension, n is the number of coefficients in p and d_i is
  the highest power used in dimension i, i = 1, 2, ..., dim.

  If the Horner evaluation method is used the complexity has not been studied
  yet.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INCOMPATIBLE_INPUT if the length of x differs from the dimension
    of the polynomial
 */
/*----------------------------------------------------------------------------*/
double cpl_polynomial_eval(const cpl_polynomial * p, const cpl_vector * x)
{
    /* Test entries */
    cpl_ensure(p, CPL_ERROR_NULL_INPUT, -1);
    cpl_ensure(x, CPL_ERROR_NULL_INPUT, -2);
    cpl_ensure(p->dim == cpl_vector_get_size(x), CPL_ERROR_INCOMPATIBLE_INPUT,
               -3);

    if (p->dim == 1) {
        return cpl_polynomial_eval_1d(p, cpl_vector_get(x, 0), NULL);
    } else if (p->nc == 0) {
        return 0;
    } else {
        return CPL_POLYNOMIAL_EVAL_FUNCTION(p,x);
    }
}

    
#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
/*----------------------------------------------------------------------------*/
/**
  @brief    Evaluate the polynomial at the given point
  @param    p  The polynomial
  @param    x  Point of evaluation
  @return   The computed value or undefined on error.

  The length of x must be the polynomial dimension.

  See http://www.uni-giessen.de/tomas.sauer/Publ/runhorner.pdf.gz
  for a study of the Horner multivariate evaluation algorithm. 
  

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INCOMPATIBLE_INPUT if the length of x differs from the dimension
    of the polynomial
 */
/*----------------------------------------------------------------------------*/
double cpl_polynomial_eval_horner(const cpl_polynomial * p,
                                  const cpl_vector     * x)
{
    /* This variable is used by the Horner recursive algorithm to know in
     * which step of the multivariate evaluation we are. */
    int            * horner_step;
    /* The ireg is another variable which specifies in which Horner "register"
     * we are interested. */
    int              ireg;
    
    int              istep;
    int              ndim;
    double           val;
    const double   * x_ptr;
    int            * powers_tmp;
    cpl_polynomial * p_nonconst;


    ndim = p->dim;

    /* Have a non-const copy */
    CPL_DIAG_PRAGMA_PUSH_IGN(-Wcast-qual);
    p_nonconst = (cpl_polynomial*)p;
    CPL_DIAG_PRAGMA_POP;
    
    /* Fill the step array */
    /* Consider making this a member (cast to non-const before) */
    horner_step = cpl_malloc(sizeof(int) * ndim);
    /* Consider a memcpy from a reference copy */
    for(istep = 0; istep < ndim; ++istep)
        horner_step[istep] = p->max_degree_alldims;
    
    /* The final value is the first "register" in the Horner tree */
    ireg = 0;
    x_ptr = cpl_vector_get_data_const(x); //For performance improvement
    powers_tmp = cpl_malloc(sizeof(int) * ndim);
    p_nonconst->added_exist_index = 0;
    p_nonconst->added_iter = 0;
    val = cpl_polynomial_register_horner(p_nonconst, ireg, horner_step, 
            powers_tmp, x_ptr);

    /* Housekeeping */
    cpl_free(horner_step);
    cpl_free(powers_tmp);
    
    return val;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Arrange the coefficients so that they are suitable to be used
            efficiently by the Horner algorithm
  @param    p  The polynomial

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
 */
/*----------------------------------------------------------------------------*/
void cpl_polynomial_arrange_coeffs_horner(cpl_polynomial       * p,
        int                    ireg,
        int                  * step,
        int                  * powers_tmp)
{
    const int       ndim = p->dim;

    /* If we are in ireg == ndim - 1 the "added" term is one of the coefficients */
    if(ireg == p->dim - 1) {
        int    i;
        int    coeff_index;

        /* Fill the powers for this iteration */
        powers_tmp[0] = p->max_degree_alldims - step[0];
        for(i = 1; i < ndim; ++i)
            powers_tmp[i] = step[i-1] - step[i];
        coeff_index = cpl_polynomial_coeff_index(p, powers_tmp);
        if(coeff_index != -1)
        {
            const double added = p->c[coeff_index];
            p->added_horner[p->added_exist_index] = added;
            p->added_exist_iterations[p->added_exist_index] = p->added_iter;
            p->added_exist_index++;
        }
        p->added_iter++;
    }
    /* If not, call recursively to get the next register */
    else
        cpl_polynomial_arrange_coeffs_horner(p, ireg + 1, step, 
                powers_tmp);

    /* If we are in the last step for this register, the multiply term is 0 */
    if(step[ireg] != 0) {
        /* We are not, recursively call with different step */

        /* Allocate and copy the new step */
        int * new_step = cpl_malloc(sizeof(int) * ndim);
        memcpy(new_step, step, sizeof(int) * ndim);
        for(int i = ireg; i<ndim; ++i)
            (new_step[i])--;
        cpl_polynomial_arrange_coeffs_horner(p, ireg, new_step, 
                powers_tmp);
        /* Housekeeping */
        cpl_free(new_step);
    }
}
    
/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief   Returns a given register in the Horner multivariate algorithm 
  @param   p          The polynomial
  @param   ireg       The index of the register of interest
  @param   step       The step within the Horner algorithm
  @param   powers_tmp Already allocated array of dim ndim for performance 
                      improvements
  @param   x          The evaluation point

  The Horner method described in 
  http://www.uni-giessen.de/tomas.sauer/Publ/runhorner.pdf.gz
  defines several steps that cover all the possible powers of the polynomial.
  Note that the variable step does not contain the powers (but it is
  closely related).
  For each step there are a set of "registers" that are used to compute other
  registers.
  A given register is defined as:
   register = multiply * x[ireg] + added
  where "multiply" and "added" are defined in terms of other "registers" (or
  the polynomial coefficients in the case of ireg == ndim - 1)
  Note that this is a recursive function and therefore it is not easy to follow. 

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
 */
/*----------------------------------------------------------------------------*/
double cpl_polynomial_register_horner(cpl_polynomial * p,
        int                    ireg,
        int                  * step,
        int                  * powers_tmp,
        const double         * x)
{
    double          multiply;
    double          added;
    int             ndim;

    /* Allocate and copy the new step */
    ndim = p->dim;

    /* If we are in ireg == ndim - 1 the "added" term is one of the coefficients */
    if(ireg == p->dim - 1)
    {
        /* Verify that this iteration has an existing coefficient */
        if(p->added_exist_iterations[p->added_exist_index] == p->added_iter)
        {
            added = p->added_horner[p->added_exist_index];
            p->added_exist_index++;
        }
        else
            added = 0;
        p->added_iter++;
    }
    /* If not, call recursively to get the next register */
    else
        added = cpl_polynomial_register_horner(p, ireg + 1, step, 
                powers_tmp, x);

    /* If we are in the last step for this register, the multiply term is 0 */
    if(step[ireg] == 0)
        multiply = 0;
    /* If not, recursively call with different step */
    else
    {
        int * new_step = cpl_malloc(sizeof(int) * ndim);
        memcpy(new_step, step, sizeof(int) * ndim);
        for(int i = ireg; i<ndim; ++i)
            (new_step[i])--;
        multiply = cpl_polynomial_register_horner(p, ireg, new_step, 
                powers_tmp, x);
        /* Housekeeping */
        cpl_free(new_step);
    }


    /* Return the value of this register */
    return (multiply * x[ireg] + added);
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Gives the index in the coefficient array  
  @param    in      the input polynomial
  @param    pows    the powers of the different variables
  @return   The index in the coefficient array if it exists, -1 elsewhere
  
  The implementation is the same as cpl_polynomial_get_coeff,
  except that it only returns the index in the array 

 */
/*----------------------------------------------------------------------------*/
static
int cpl_polynomial_coeff_index(const cpl_polynomial    *   in,
                               const int               *   pows)
{
    int dim;
    int idx = -1;


    /* Test entries */
    cpl_ensure(in,  CPL_ERROR_NULL_INPUT, idx);
    cpl_ensure(pows, CPL_ERROR_NULL_INPUT, idx);

    for (dim = 0;  dim < in->dim; dim++) {
        cpl_error_ensure(pows[dim] >= 0, CPL_ERROR_ILLEGAL_INPUT, return(idx),
                         "Dimension %d of %lld has negative power %d", dim+1,
                         in->dim, pows[dim]);
    }

    /* Find the coeff */
    if (in->dim == 1) {
        if (pows[0] < in->nc) idx = pows[0];
    } else {
        int i;
        for (i=0; i < in->nc; i++) {
            if (!memcmp(in->pow + in->dim * i, pows, in->dim * sizeof(int)))
                break; /* Found the right combination of powers */
        }
        if (i < in->nc) idx = i;
    }
    return idx;
}

#else

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Evaluate the polynomial at the given point (brute force algorithm)
  @param    p  The polynomial
  @param    x  Point of evaluation
  @return   The computed value or undefined on error.

  The length of x must be the polynomial dimension.

  A polynomial with no non-zero coefficients evaluates as 0.

  For 1-dimensional polynomials the call requires 2n FLOPs
  where n+1 is the number of coefficients in p, see also
  cpl_polynomial_eval_1d().

  For multivariate polynomials the call requires
  n*(1+dim) + d_1 + d_2 + ... + d_dim FLOPs,
  where dim is the dimension, n is the number of coefficients in p and d_i is
  the highest power used in dimension i, i = 1, 2, ..., dim.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INCOMPATIBLE_INPUT if the length of x differs from the dimension
    of the polynomial
 */
/*----------------------------------------------------------------------------*/
static
double cpl_polynomial_eval_bf(const cpl_polynomial * p, const cpl_vector * x)
{

    const double * x_val = cpl_vector_get_data_const(x);
    double         z     = 0.0;
    cpl_size       i, j;

    /* Build the table of powers of the evaluation point */
    for (j=0; j < p->dim; j++) {
        /* Allow changes to (the eval-buffer of) p!
               This is acceptable, because the data that is modified
               here is not accessed outside this else-block.
               FIXME: Revise this approach for shared memory parallelism.
         */

        double * eval_powj  = (double *) p->eval_pow[j];

        /* Already done: eval_powj[0] = 1.0; */
        for (i=1; i <= p->max_degree[j]; i++)
            eval_powj[i] = eval_powj[i-1] * x_val[j];

        cpl_tools_add_flops( p->max_degree[j]);

    }

    /* Evaluate the polynomial using the table.
           If the p->eval_pow does not fit in the cache,
           its memory layout should be changed */
    for (i=0; i < p->nc; i++) {
        /* dim is at least 2 */
        double val = p->eval_pow[0][p->pow[p->dim * i + 0]]
                   * p->eval_pow[1][p->pow[p->dim * i + 1]];
        for (j = 2; j < p->dim; j++) {
            val *= p->eval_pow[j][p->pow[p->dim * i + j]];
        }
        z += p->c[i] * val;
    }
    cpl_tools_add_flops( p->nc * ( 1 + p->dim ) );
    return z;

}

#endif

/*----------------------------------------------------------------------------*/
/**
  @brief   Collapse one dimension of a multi-variate polynomial by composition
  @param   self  The multi-variate polynomial
  @param   dim   The dimension to collapse (zero for first dimension)
  @param   other The polynomial to replace dimension dim of self
  @return  The collapsed polynomial or NULL on error

  The dimension of the polynomial self must be one greater than that of the
  other polynomial. Given these two polynomials the dimension dim of self is
  collapsed by creating a new polynomial from
    self(x0, x1, ..., x{dim-1},
         other(x0, x1, ..., x{dim-1}, x{dim+1}, x{dim+2}, ..., x{n-1}),
         x{dim+1}, x{dim+2}, ..., x{n-1}).

  The created polynomial thus has a dimension which is one less than the
  polynomial self and which is equal to that of the other polynomial.
  Collapsing one dimension of a 1D-polynomial is equivalent to evaluating it,
  which can be done with cpl_polynomial_eval_1d().

  FIXME: The other polynomial must currently have a degree of zero, i.e. it must
  be a constant.

  Currently,
  the call requires dn + p FLOPs, where d the dimension of the polynomial self,
  p is the largest power of dimension dim and n the number of (non-zero)
  coefficients of the polynomial self.

  The returned object is a newly allocated cpl_polynomial that
  must be deallocated by the caller using cpl_polynomial_delete().

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INVALID_TYPE if the polynomial is uni-variate.
  - CPL_ERROR_ILLEGAL_INPUT if dim is negative.
  - CPL_ERROR_ACCESS_OUT_OF_RANGE if dim exceeds the dimension of self.
  - CPL_ERROR_INCOMPATIBLE_INPUT if other has the wrong dimension.
  - CPL_ERROR_UNSUPPORTED_MODE if other is not of degree 0.

 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cpl_polynomial_extract(const cpl_polynomial * self,
                                        cpl_size dim,
                                        const cpl_polynomial * other)
{

    cpl_polynomial * collapsed;
    cpl_size         newdim;


    cpl_ensure(self  != NULL,    CPL_ERROR_NULL_INPUT,          NULL);
    cpl_ensure(other != NULL,    CPL_ERROR_NULL_INPUT,          NULL);
    cpl_ensure(self->dim > 1,    CPL_ERROR_INVALID_TYPE,        NULL);
    cpl_ensure(dim >= 0,         CPL_ERROR_ILLEGAL_INPUT,       NULL);
    cpl_ensure(dim < self->dim,  CPL_ERROR_ACCESS_OUT_OF_RANGE, NULL);

    newdim = self->dim - 1;

    cpl_ensure(other->dim == newdim, CPL_ERROR_INCOMPATIBLE_INPUT, NULL);

    /* FIXME: Generalize this */
    cpl_ensure(cpl_polynomial_get_degree(other) == 0,
               CPL_ERROR_UNSUPPORTED_MODE, NULL);

    collapsed = cpl_polynomial_new(newdim);

    if (self->nc > 0) {
        double     * eval_powj;
        cpl_size   * pows = (cpl_size*)cpl_malloc((size_t)newdim
                                                  * sizeof(*pows));
        const double x = other->nc == 0 ? 0.0 : other->c[0];
        cpl_size     i;


        /* Build the table of powers of the evaluation point */
        /* Allow changes to (the eval-buffer of) self!
           This is acceptable, because the data that is modified
           here is not accessed outside this else-block.
           FIXME: Revise this approach for shared memory parallelism.
        */

        eval_powj  = (double *) self->eval_pow[dim];

        /* Already done: eval_powj[0] = 1.0; */
        for (i=1; i <= self->max_degree[dim]; i++)
            eval_powj[i] = eval_powj[i-1] * x;

        for (i=0; i < self->nc; i++) {
            double coeff = self->c[i];
            cpl_size j, k;

            for (k=0, j=0; j < self->dim; j++) {
                if (j == dim) {
                    coeff *= eval_powj[self->pow[self->dim * i + j]];
                } else {
                    pows[k++] = self->pow[self->dim * i + j];
                }
            }

            coeff += cpl_polynomial_get_coeff(collapsed, pows);

            (void)cpl_polynomial_set_coeff(collapsed, pows, coeff);

        }

        cpl_tools_add_flops( self->max_degree[dim] + self->nc * self->dim );

        cpl_free(pows);
    }

    return collapsed;

}


/*----------------------------------------------------------------------------*/
/**
  @brief    Compute a first order partial derivative
  @param    self  The polynomial to be modified in place
  @param    dim   The dimension to differentiate (zero for first dimension)
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_

  The dimension of the polynomial is preserved, even if the operation may cause
  the polynomial to become independent of the dimension dim of the variable.

  The call requires n FLOPs, where n is the number of (non-zero) polynomial
  coefficients whose power in dimension dim is at least 1.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_ILLEGAL_INPUT if dim is negative.
  - CPL_ERROR_ACCESS_OUT_OF_RANGE if dim exceeds the dimension of self.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_derivative(cpl_polynomial * self, cpl_size dim)
{

    cpl_size i;


    cpl_ensure_code(self != NULL,    CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(dim >= 0,        CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(dim < self->dim, CPL_ERROR_ACCESS_OUT_OF_RANGE);

    if (self->nc == 0) return CPL_ERROR_NONE;

    if (self->dim == 1) {
        for (i=0; i < self->nc-1; i++) {
            self->c[i] = self->c[i+1] * (double) (i+1);
        }
        cpl_polynomial_delete_coeff(self, &i);
    } else {

        /* Remove the terms that are constant with respect to dim */
        for (i=0; i < self->nc; ) {
            if (self->pow[self->dim * i + dim] == 0) {
                /* FIXME: Memory aliasing occurs here - the contents of the
                   const pointer is actually modified */
                cpl_polynomial_delete_coeff(self, self->pow + self->dim * i);
            } else {
                i++;
            }
        }

        /* At this point self contains the correct number of terms
           which now need to be differentiated */
        if (self->nc > 0) {
            for (i=0; i < self->nc; i++) {
                assert(self->pow[self->dim * i + dim] > 0);
                self->c[i] *= (double) self->pow[self->dim * i + dim]--;
            }

            if (self->max_degree[dim] > 0) {
                /* Decrement the size of the work-array */

                self->eval_pow[dim] =
                    cpl_realloc(self->eval_pow[dim],
                                (size_t)self->max_degree[dim] * sizeof(double));

                self->max_degree[dim]--;
                self->max_degree_alldims--;

            }
        }
    }

    cpl_tools_add_flops( self->nc );

    return CPL_ERROR_NONE;

}



/* A nested Horner algorithm can be used to evaluate both differences and
   derivatives.
   Caveat: All four arguments are evaluated more than once */
#define CPL_HORNER_NESTED(A, B, PA, PB)         \
    do {                                        \
        PB = PA;                                \
        while (n> 1) {                          \
            PA = PA * A + self->c[--n];         \
            PB = PB * B + PA;                   \
        }                                       \
        PA = PA * A + self->c[0];               \
                                                \
    } while (0)

/*----------------------------------------------------------------------------*/
/**
  @brief    Evaluate a univariate (1D) polynomial using Horners rule.
  @param    self  The 1D-polynomial
  @param    x     The point of evaluation
  @param    pd    Iff pd is non-NULL, the derivative evaluated at x
  @return   The result or undefined on error.

  A polynomial with no non-zero coefficents evaluates to 0 with a
  derivative that does likewise.

  The result is computed as p_0 + x * ( p_1 + x * ( p_2 + ... x * p_n ))
  and requires 2n FLOPs where n+1 is the number of coefficients.

  If the derivative is requested it is computed using a nested Horner rule.
  This requires 4n FLOPs in total.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INVALID_TYPE if the polynomial is not (1D) univariate
 */
/*----------------------------------------------------------------------------*/
double cpl_polynomial_eval_1d(const cpl_polynomial * self, double x,
                              double * pd)
{

    cpl_long_double result;
    cpl_size        n;

    cpl_ensure(self      != NULL, CPL_ERROR_NULL_INPUT,   -1);
    cpl_ensure(self->dim == 1,    CPL_ERROR_INVALID_TYPE, -3);

    if (self->nc == 0) {
        if (pd != NULL) *pd = 0;
        return 0;
    }

    n = self->nc;

    result = self->c[--n];

    if (pd == NULL) {
        cpl_tools_add_flops( 2 * n );
        while (n) result = x * result + self->c[--n];
    } else {
        cpl_long_double d = 0;

        cpl_tools_add_flops( 4 * n );
        if (n) CPL_HORNER_NESTED(x, x, result, d);
        *pd = d;
    }

    return result;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Evaluate p(a) - p(b) using Horners rule.
  @param    self  The 1D-polynomial
  @param    a     The evaluation point of the minuend
  @param    b     The evaluation point of the subtrahend
  @param    ppa   Iff ppa is not NULL, p(a)
  @return   The difference or undefined on error

  The call requires about 4n FLOPs where n is the number of coefficients in
  self, which is the same as that required for two separate polynomial
  evaluations. cpl_polynomial_eval_1d_diff() is however more accurate.

  ppa may be NULL. If it is not, *ppa is set to self(a), which is calculated at
  no extra cost.

  The underlying algorithm is the same as that used in cpl_polynomial_eval_1d()
  when the derivative is also requested.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INVALID_TYPE if the polynomial has the wrong dimension
 */
/*----------------------------------------------------------------------------*/
double cpl_polynomial_eval_1d_diff(const cpl_polynomial * self, double a,
                                   double b, double * ppa)
{

    cpl_long_double diff, pa;
    cpl_size        n;


    cpl_ensure(self      != NULL, CPL_ERROR_NULL_INPUT,   -1);
    cpl_ensure(self->dim == 1,    CPL_ERROR_INVALID_TYPE, -3);

    if (self->nc == 0) {
        if (ppa != NULL) *ppa = 0;
        return 0;
    }

    n = self->nc;

    pa = self->c[--n];

    cpl_tools_add_flops( 4 * n );
    CPL_HORNER_NESTED(a, b, pa, diff);

    if (ppa != NULL) *ppa = pa;

    return diff * (a - b);

}
#undef CPL_HORNER_NESTED

/*----------------------------------------------------------------------------*/
/**
  @brief    Evaluate a 1D-polynomial on equidistant points using Horners rule
  @param    v  Preallocated vector to contain the result
  @param    p  The 1D-polynomial
  @param    x0 The first point of evaluation
  @param    d  The increment between points of evaluation
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_
  @see cpl_vector_fill

  The evaluation points are x_i = x0 + i * d, i=0, 1, ..., n-1,
  where n is the length of the vector.

  If d is zero it is preferable to simply use
  cpl_vector_fill(v, cpl_polynomial_eval_1d(p, x0, NULL)).

  The call requires about 2nm FLOPs, where m+1 is the number of coefficients in
  p.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INVALID_TYPE if the polynomial has the wrong dimension
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_vector_fill_polynomial(cpl_vector * v,
                                          const cpl_polynomial * p,
                                          double x0, double d)
{
    cpl_size i  = cpl_vector_get_size(v);
    double * dv = cpl_vector_get_data(v);

    cpl_ensure_code(v != NULL,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(p != NULL,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(p->dim == 1, CPL_ERROR_INVALID_TYPE);

    do {
        i--;
        dv[i] = cpl_polynomial_eval_1d(p, x0 + (double)i * d, NULL);
    } while (i > 0);

    cpl_tools_add_flops( 2 * cpl_vector_get_size(v) );

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    A real solution to p(x) = 0 using Newton-Raphsons method
  @param    p    The 1D-polynomial
  @param    x0   First guess of the solution
  @param    px   The solution, on error see below
  @param    mul  The root multiplicity (or 1 if unknown)
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_

  Even if a real solution exists, it may not be found if the first guess is
  too far from the solution. But a solution is guaranteed to be found if all
  roots of p are real. If the constant term is zero, the solution 0 will be
  returned regardless of the first guess.

  No solution is found when the iterative process stops because:
  1) It can not proceed because p`(x) = 0 (CPL_ERROR_DIVISION_BY_ZERO).
  2) Only a finite number of iterations are allowed (CPL_ERROR_CONTINUE).
  Both cases may be due to lack of a real solution or a bad first guess.
  In these two cases *px is set to the value where the error occurred.
  In case of other errors *px is unmodified.

  The accuracy and robustness deteriorates with increasing multiplicity
  of the solution. This is also the case with numerical multiplicity,
  i.e. when multiple solutions are located close together.

  mul is assumed to be the multiplicity of the solution. Knowledge of the
  root multiplicity often improves the robustness and accuracy. If there
  is no knowledge of the root multiplicity mul should be 1.
  Setting mul to a too high value should be avoided.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INVALID_TYPE if the polynomial has the wrong dimension
  - CPL_ERROR_ILLEGAL_INPUT if the multiplicity is non-positive
  - CPL_ERROR_DIVISION_BY_ZERO if a division by zero occurs
  - CPL_ERROR_CONTINUE if the algorithm does not converge
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_solve_1d(const cpl_polynomial * p,
                                       double x0, double * px, cpl_size mul)
{

    return cpl_polynomial_solve_1d_(p, x0, px, mul, CPL_FALSE)
        ? cpl_error_set_where_() : CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Fit a polynomial to a set of samples in a least squares sense
  @param  self    Polynomial of dimension d to hold the coefficients
  @param  samppos Matrix of p sample positions, with d rows and p columns
  @param  sampsym NULL, or d booleans, true iff the sampling is symmetric
  @param  fitvals Vector of the p values to fit
  @param  fitsigm Uncertainties of the sampled values, or NULL for all ones
  @param  dimdeg  True iff there is a fitting degree per dimension
  @param  mindeg  Pointer to 1 or d minimum fitting degree(s), or NULL
  @param  maxdeg  Pointer to 1 or d maximum fitting degree(s), at least mindeg
  @return CPL_ERROR_NONE on success, else the relevant #_cpl_error_code_
  @note Currently only uni- and bi-variate polynomials are supported,
        fitsigm must be NULL. For all but uni-variate polynomials mindeg must
        be zero.
  @see cpl_vector_fill_polynomial_fit_residual()

  Any pre-set non-zero coefficients in self are overwritten or reset by the fit.

  For 1D-polynomials N = 1 + maxdeg - mindeg coefficients are fitted.

  For multi-variate polynomials the fit depends on dimdeg:

  If dimdeg is false, an n-degree coefficient is fitted iff
  mindeg <= n <= maxdeg. For a 2D-polynomial this means that
  N * (N + 1) / 2 coefficients are fitted.

  If dimdeg is true, nci = 1 + maxdeg[i] + mindeg[i] coefficients are fitted
  for dimension i, i.e. for a 2D-polynomial N = nc1 * nc2 coefficients are
  fitted.

  The number of distinct samples should exceed the number of coefficients to
  fit. The number of distinct samples may also equal the number of coefficients
  to fit, but in this case the fit has another meaning (any non-zero residual
  is due to rounding errors, not a fitting error).
  It is an error to try to fit more coefficients than there are distinct
  samples.

  If the relative uncertainties of the sampled values are known, they may be
  passed via fitsigm. NULL means that all uncertainties equals one.

  sampsym is ignored if mindeg is nonzero, otherwise
  the caller may use sampsym to indicate an a priori knowledge that the sampling
  positions are symmetric. NULL indicates no knowledge of such symmetry.
  sampsym[i] may be set to true iff the sampling is symmetric around u_i, where
  u_i is the mean of the sampling positions in dimension i.

  In 1D this implies that the sampling points as pairs average u_0 (with an odd
  number of samples one sample must equal u_0). E.g. both x = (1, 2, 4, 6, 7)
  and x = (1, 6, 4, 2, 7) have sampling symmetry, while x = (1, 2, 4, 6) does
  not.

  In 2D this implies that the sampling points are symmetric in the 2D-plane.
  For the first dimension sampling symmetry means that the sampling is line-
  symmetric around y = u_1, while for the second dimension, sampling symmetry
  implies line-symmetry around x = u_2. Point symmetry around
  (x,y) = (u_1, u_2) means that both sampsym[0] and sampsym[1] may be set to
  true.

  Knowledge of symmetric sampling allows the fit to be both faster and
  eliminates certain round-off errors.

  For higher order fitting the fitting problem known as "Runge's phenomenon"
  is minimized using the socalled "Chebyshev nodes" as sampling points.
  For Chebyshev nodes sampsym can be set to CPL_TRUE.

  Warning: An increase in the polynomial degree will normally reduce the
  fitting error. However, due to rounding errors and the limited accuracy
  of the solver of the normal equations, an increase in the polynomial degree
  may at some point cause the fitting error to _increase_. In some cases this
  happens with an increase of the polynomial degree from 8 to 9.

  The fit is done in the following steps:
  1) If mindeg is zero, the sampling positions are first transformed into
     Xhat_i = X_i - mean(X_i), i=1, .., dimension.
  2) The Vandermonde matrix is formed from Xhat.
  3) The normal equations of the Vandermonde matrix is solved.
  4) If mindeg is zero, the resulting polynomial in Xhat is transformed
      back to X.

  For a univariate (1D) fit the call requires 6MN + N^3/3 + 7/2N^2 + O(M) FLOPs
  where M is the number of data points and where N is the number of polynomial
  coefficients to fit, N = 1 + maxdeg - mindeg.

  For a bivariate fit the call requires MN^2 + N^3/3 + O(MN) FLOPs where M
  is the number of data points and where N is the number of polynomial
  coefficients to fit.

  Examples of usage:
    @code

    cpl_polynomial  * fit1d     = cpl_polynomial_new(1);
    cpl_matrix      * samppos1d = my_sampling_points_1d(); // 1-row matrix
    cpl_vector      * fitvals   = my_sampling_values();
    const cpl_boolean sampsym   = CPL_TRUE;
    const int         maxdeg1d  = 4; // Fit 5 coefficients
    cpl_error_code    error1d
        = cpl_polynomial_fit(fit1d, samppos1d, &sampsym, fitvals, NULL,
                             CPL_FALSE, NULL, &maxdeg1d);
    @endcode


    @code

    cpl_polynomial  * fit2d      = cpl_polynomial_new(2);
    cpl_matrix      * samppos2d  = my_sampling_points_2d(); // 2-row matrix
    cpl_vector      * fitvals    = my_sampling_values();
    const int         maxdeg2d[] = {2, 1}; // Fit 6 coefficients
    cpl_error_code    error2d
        = cpl_polynomial_fit(fit2d, samppos2d, NULL, fitvals, NULL, CPL_FALSE,
                             NULL, maxdeg2d);

    @endcode

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_ILLEGAL_INPUT if a mindeg value is negative, or if a maxdeg value
    is less than the corresponding mindeg value.
  - CPL_ERROR_DATA_NOT_FOUND if the number of columns in samppos is less than
    the number of coefficients to be determined.
  - CPL_ERROR_INCOMPATIBLE_INPUT if samppos, fitvals or fitsigm have
    incompatible sizes, or if samppos, self or sampsym have incompatible sizes.
  - CPL_ERROR_SINGULAR_MATRIX if samppos contains too few distinct values
  - CPL_ERROR_DIVISION_BY_ZERO if an element in fitsigm is zero
  - CPL_ERROR_UNSUPPORTED_MODE if the polynomial dimension exceeds two, or if
    there is a non-zero value in the mindeg array.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_fit(cpl_polynomial    * self,
                                  const cpl_matrix  * samppos,
                                  const cpl_boolean * sampsym,
                                  const cpl_vector  * fitvals,
                                  const cpl_vector  * fitsigm,
                                  cpl_boolean         dimdeg,
                                  const cpl_size    * mindeg,
                                  const cpl_size    * maxdeg)
{

    const cpl_size mdim = cpl_polynomial_get_dimension(self);
    const cpl_size np   = cpl_vector_get_size(fitvals);

    cpl_ensure_code(self    != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(samppos != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(fitvals != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(maxdeg  != NULL, CPL_ERROR_NULL_INPUT);

    if (cpl_matrix_get_ncol(samppos) != np)
        return cpl_error_set_message_(CPL_ERROR_INCOMPATIBLE_INPUT, "Number of "
                                      "fitting values = %" CPL_SIZE_FORMAT
                                      " <=> % " CPL_SIZE_FORMAT " = samppos col"
                                      "umns", np, cpl_matrix_get_ncol(samppos));

    if (cpl_matrix_get_nrow(samppos) != mdim)
        return cpl_error_set_message_(CPL_ERROR_INCOMPATIBLE_INPUT, "Fitting "
                                      "dimension = %" CPL_SIZE_FORMAT " <=> %"
                                      CPL_SIZE_FORMAT " = samppos rows", mdim,
                                      cpl_matrix_get_nrow(samppos));

    if (fitsigm != NULL && cpl_vector_get_size(fitsigm) != np)
        return cpl_error_set_message_(CPL_ERROR_INCOMPATIBLE_INPUT, "Number of "
                                      "fitting values = %" CPL_SIZE_FORMAT " <="
                                      "> % " CPL_SIZE_FORMAT " = number of er"
                                      "rors", np, cpl_vector_get_size(fitsigm));

    if (fitsigm != NULL)
        return cpl_error_set_message_(CPL_ERROR_UNSUPPORTED_MODE, "fitsigm = "
                                      "%p != NULL is not supported (yet)",
                                      (const void*)fitsigm);

CPL_DIAG_PRAGMA_PUSH_IGN(-Wcast-qual);
    if (mdim == 1) {
        cpl_vector * x_pos
            = cpl_vector_wrap(np, cpl_matrix_get_data((cpl_matrix*)samppos));
        const cpl_boolean isampsym = sampsym ? sampsym[0] : CPL_FALSE;
        const int         mindeg0  = mindeg ? mindeg[0] : 0;
        const cpl_error_code error = cpl_polynomial_fit_1d(self, x_pos, fitvals,
                                                           mindeg0, maxdeg[0],
                                                           isampsym, NULL);
CPL_DIAG_PRAGMA_POP;

        (void)cpl_vector_unwrap(x_pos);
        if (error) return cpl_error_set_where_();

    } else if (mdim > 2) {
        return cpl_error_set_message_(CPL_ERROR_UNSUPPORTED_MODE, "The fitting "
                                      "dimension %d > 2 is not supported",
                                      (int)mdim);
CPL_DIAG_PRAGMA_PUSH_IGN(-Wcast-qual);
    } else if ((!mindeg || (mindeg[0] == 0 && (!dimdeg || mindeg[1] == 0)))) {
        cpl_vector * x_pos
            = cpl_vector_wrap(np, cpl_matrix_get_data((cpl_matrix *)samppos));
        cpl_vector * y_pos = cpl_vector_wrap(np, np
                              + cpl_matrix_get_data((cpl_matrix *)samppos));
        cpl_bivector * xy_pos = cpl_bivector_wrap_vectors(x_pos, y_pos);
        const cpl_error_code error = cpl_polynomial_fit_2d(self, xy_pos,
                                                           fitvals, dimdeg,
                                                           maxdeg, NULL);
CPL_DIAG_PRAGMA_POP;
        assert( mdim == 2 );

        (void)cpl_vector_unwrap(x_pos);
        (void)cpl_vector_unwrap(y_pos);
        cpl_bivector_unwrap_vectors(xy_pos);

        if (error) return cpl_error_set_where_();

    } else {
        assert( mdim == 2 );
        return cpl_error_set_message_(CPL_ERROR_UNSUPPORTED_MODE, "In a 2D-fit "
                                      "mindeg must be NULL or contain zero(s)");
    }

    return CPL_ERROR_NONE;
}



/*----------------------------------------------------------------------------*/
/**
  @brief  Compute the residual of a polynomial fit
  @param  self    Vector to hold the fitting residuals, fitvals may be used
  @param  fitvals Vector of the p fitted values
  @param  fitsigm Uncertainties of the sampled values or NULL for a uniform
                  uncertainty
  @param  fit     The fitted polynomial
  @param  samppos Matrix of p sample positions, with d rows and p columns
  @param  rechisq If non-NULL, the reduced chi square of the fit
  @return CPL_ERROR_NONE on success, else the relevant #_cpl_error_code_
  @note If necessary, self is resized to the length of fitvals.
  @see cpl_polynomial_fit()

  It is allowed to pass the same vector as both fitvals and as self,
  in which case fitvals is overwritten with the residuals.

  If the relative uncertainties of the sampled values are known, they may be
  passed via fitsigm. NULL means that all uncertainties equal one. The
  uncertainties are taken into account when computing the reduced
  chi square value.

  If rechisq is non-NULL, the reduced chi square of the fit is computed as
  well.

  The mean square error, which was computed directly by the former CPL functions
  cpl_polynomial_fit_1d_create() and cpl_polynomial_fit_2d_create() can be
  computed from the fitting residual like this:

  @code
    const double mse = cpl_vector_product(fitresidual, fitresidual)
                     / cpl_vector_get_size(fitresidual);
  @endcode

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer (other than fitsigm) is NULL
  - CPL_ERROR_INCOMPATIBLE_INPUT if samppos, fitvals, fitsigm or fit
    have incompatible sizes
  - CPL_ERROR_DIVISION_BY_ZERO if an element in fitsigm is zero
  - CPL_ERROR_DATA_NOT_FOUND if the number of columns in samppos is less than
    the number of coefficients in the fitted polynomial.
 */
/*----------------------------------------------------------------------------*/
cpl_error_code
cpl_vector_fill_polynomial_fit_residual(cpl_vector           * self,
                                        const cpl_vector     * fitvals,
                                        const cpl_vector     * fitsigm,
                                        const cpl_polynomial * fit,
                                        const cpl_matrix     * samppos,
                                        double               * rechisq)
{
    const cpl_size mdim = cpl_polynomial_get_dimension(fit);
    const cpl_size np = cpl_vector_get_size(fitvals);
    cpl_vector   * sampoint;
    double       * dsampoint;
    double       * dself;
    cpl_size       i, j;

    cpl_ensure_code(self    != NULL,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(fitvals != NULL,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(fit     != NULL,   CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(samppos != NULL,   CPL_ERROR_NULL_INPUT);

    cpl_ensure_code(cpl_matrix_get_ncol(samppos) == np,
                    CPL_ERROR_INCOMPATIBLE_INPUT);

    cpl_ensure_code(cpl_matrix_get_nrow(samppos) == mdim,
                    CPL_ERROR_INCOMPATIBLE_INPUT);

    /* May not be under-determined */
    cpl_ensure_code(np >= fit->nc, CPL_ERROR_DATA_NOT_FOUND);

    cpl_vector_set_size(self, np);
    dself = cpl_vector_get_data(self);

    sampoint = cpl_vector_new(mdim);
    dsampoint = cpl_vector_get_data(sampoint);

    /* FIXME: Use direct pointer access for matrix ? */
    for (i = 0; i < np; i++) {
        for (j = 0; j < mdim; j++) {
            dsampoint[j] = cpl_matrix_get(samppos, j, i);
        }
        dself[i] = cpl_vector_get(fitvals, i)
            - cpl_polynomial_eval(fit, sampoint);
    }

    cpl_vector_delete(sampoint);

    cpl_tools_add_flops( np );

    if (rechisq != NULL) {
        /* Assuming that the np sampling points are distinct! */
        const int nfree = np - fit->nc;

        cpl_ensure_code(nfree > 0, CPL_ERROR_DATA_NOT_FOUND);

        if (fitsigm == NULL) {
            *rechisq = cpl_vector_product(self, self) / nfree;
        } else {
            const double * dsigm = cpl_vector_get_data_const(fitsigm);
            double dot = 0.0;

            cpl_ensure_code(cpl_vector_get_size(fitsigm) == np,
                            CPL_ERROR_INCOMPATIBLE_INPUT);

            for (i = 0; i < np; i++) {
                double delta;
                /* Sigmas may be negative, the sign is ignored */
                cpl_ensure_code(dsigm[i] != 0.0, CPL_ERROR_DIVISION_BY_ZERO);
                delta = dself[i] / dsigm[i];
                dot += delta * delta;
            }
            *rechisq = dot / nfree;
            cpl_tools_add_flops( 3 * np );
        }
    }

    return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
  @brief  Modify p, p(x0, x1, ..., xi, ...) := (x0, x1, ..., xi+u, ...)
  @param  p    The polynomial to be modified in place
  @param  i    The dimension to shift (0 for first)
  @param  u    The shift
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_
  @note Currently, only dimensions 1 and 2 are supported

  Shifting the polynomial p(x) = x^n with u = 1 will generate the binomial
  coefficients for n.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_ILLEGAL_INPUT if i is negative
  - CPL_ERROR_ACCESS_OUT_OF_RANGE if i exceeds the dimension of p
  - CPL_ERROR_UNSUPPORTED_MODE if the polynomial has a dimension larger than 2
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_shift_1d(cpl_polynomial * p, cpl_size i, double u)
{

    const cpl_size ndim = cpl_polynomial_get_dimension(p);

    cpl_ensure_code(p != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(i >= 0,    CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(i < ndim,  CPL_ERROR_ACCESS_OUT_OF_RANGE);

    if (ndim == 1) {
        /* Should not be able to fail now */
        if (p->nc) cpl_polynomial_shift_double( p->c, p->nc, u);
    } else if (ndim == 2) {
        const cpl_size   ndeg = cpl_polynomial_get_degree(p);
        cpl_polynomial * p1d  = cpl_polynomial_new(1);
        cpl_size         powers[2];
        cpl_size       * pi2  = powers + 1 - i; /* Other power */

        cpl_ensure_code(p1d != NULL, CPL_ERROR_NULL_INPUT);
        /* Should not be able to fail now */

        for (*pi2 = ndeg; *pi2 >=0; (*pi2)--) {
            cpl_size * pshift = powers + i;

            /* Copy to a 1D all coefficients with this power*/
            /* - *pi2 decreases so all coeffs in p1d are overwritten */
            /* - *pshift decreases to avoid realloc() in _set_coeff() */
            for (*pshift = ndeg - *pi2; *pshift >= 0; (*pshift)--) {
                assert(powers[i] == *pshift);
                cpl_polynomial_set_coeff(p1d, pshift,
                                         cpl_polynomial_get_coeff(p, powers));
            }

            if (p1d->nc) cpl_polynomial_shift_double( p1d->c, p1d->nc, u);

            /* Copy back, overwriting any old coeffs */
            for (*pshift = ndeg - *pi2; *pshift >= 0; (*pshift)--) {
                assert(powers[i] == *pshift);
                cpl_polynomial_set_coeff(p, powers,
                                         cpl_polynomial_get_coeff(p1d, pshift));
            }
        }

        cpl_polynomial_delete(p1d);
    } else {
        return cpl_error_set_(CPL_ERROR_UNSUPPORTED_MODE);
    }

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Fit a 1D-polynomial to a 1D-signal in a least squares sense
  @param    x_pos   Vector of positions of the signal to fit.
  @param    values  Vector of values of the signal to fit.
  @param    degree  Non-negative polynomial degree.
  @param    mse     Iff mse is not null, the mean squared error on success
  @return   The fitted polynomial or NULL on error
  @see cpl_polynomial_fit()
  @deprecated Replace this call with cpl_polynomial_fit() and
     optionally cpl_vector_fill_polynomial_fit_residual().

 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cpl_polynomial_fit_1d_create(const cpl_vector    *   x_pos,
                                              const cpl_vector    *   values,
                                              cpl_size                degree,
                                              double              *   mse)
{

    cpl_polynomial * self = cpl_polynomial_new(1);
    const cpl_boolean is_eqdist = cpl_vector_is_eqdist(x_pos) == 1;
    const cpl_error_code error = cpl_polynomial_fit_1d(self, x_pos, values, 0,
                                                       degree, is_eqdist, mse);

    if (error != CPL_ERROR_NONE) {
        cpl_polynomial_delete(self);
        self = NULL;
        cpl_error_set_(error);
    }

    return self;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Fit a 2D-polynomial to a 2D-surface in a least squares sense
  @param    xy_pos  Bivector  positions of the surface to fit.
  @param    values  Vector of values of the surface to fit.
  @param    degree  Non-negative polynomial degree.
  @param    mse     Iff mse is not null, the mean squared error on success
  @return   The fitted polynomial or NULL on error
  @see cpl_polynomial_fit()
  @deprecated Replace this call with cpl_polynomial_fit() and
     optionally cpl_vector_fill_polynomial_fit_residual().

 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cpl_polynomial_fit_2d_create(cpl_bivector     *  xy_pos,
                                              cpl_vector       *  values,
                                              cpl_size            degree,
                                              double           *  mse)
{
    cpl_polynomial * self = cpl_polynomial_new(2);
    const cpl_error_code error = cpl_polynomial_fit_2d(self, xy_pos, values,
                                                       CPL_FALSE, &degree, mse);

    if (error != CPL_ERROR_NONE) {
        cpl_polynomial_delete(self);
        self = NULL;
        cpl_error_set_(error);
    }

    return self;

}

/*----------------------------------------------------------------------------*/
/**
  @brief    Add two polynomials of the same dimension
  @param    self      The polynomial to hold the result
  @param    first     The 1st polynomial to add
  @param    second    The 2nd polynomial to add
  @return   CPL_ERROR_NONE or the relevant CPL error code
  @note self may be passed also as first and/or second

  Possible CPL error code set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INCOMPATIBLE_INPUT if the polynomials do not have identical
    dimensions
  - CPL_ERROR_UNSUPPORTED_MODE if the dimension is not 1 (FIXME)
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_add(cpl_polynomial * self,
                                  const cpl_polynomial * first,
                                  const cpl_polynomial * second)
{
    cpl_size       degree0 = cpl_polynomial_get_degree(self);
    const cpl_size degree1 = cpl_polynomial_get_degree(first);
    const cpl_size degree2 = cpl_polynomial_get_degree(second);
    const cpl_size maxdeg  = degree1 > degree2 ? degree1 : degree2;


    cpl_ensure_code(self   != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(first  != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(second != NULL, CPL_ERROR_NULL_INPUT);

    cpl_ensure_code(cpl_polynomial_get_dimension(self) ==
                    cpl_polynomial_get_dimension(first),
                    CPL_ERROR_INCOMPATIBLE_INPUT);
    cpl_ensure_code(cpl_polynomial_get_dimension(self) ==
                    cpl_polynomial_get_dimension(second),
                    CPL_ERROR_INCOMPATIBLE_INPUT);

    /* FIXME: */
    cpl_ensure_code(cpl_polynomial_get_dimension(self) == 1,
                    CPL_ERROR_UNSUPPORTED_MODE);

    if (degree0 < maxdeg) {
        degree0 = maxdeg;
    } else {
        /* Reset coefficients in self as needed */
        for (; degree0 > maxdeg; degree0--) {
            cpl_polynomial_delete_coeff(self, &degree0);
        }
    }

    assert( degree0 == maxdeg );

    for (; degree0 >= 0; degree0--) {
        const double val1 = cpl_polynomial_get_coeff(first, &degree0);
        const double val2 = cpl_polynomial_get_coeff(second, &degree0);
        cpl_polynomial_set_coeff(self, &degree0, val1 + val2);
    }

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Subtract two polynomials of the same dimension
  @param    self      The polynomial to hold the result
  @param    first     The polynomial to subtract from
  @param    second    The polynomial to subtract
  @return   CPL_ERROR_NONE or the relevant CPL error code
  @note self may be passed also as first and/or second

  Possible CPL error code set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_INCOMPATIBLE_INPUT if the polynomials do not have identical
    dimensions
  - CPL_ERROR_UNSUPPORTED_MODE if the dimension is not 1 (FIXME)
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_subtract(cpl_polynomial * self,
                                       const cpl_polynomial * first,
                                       const cpl_polynomial * second)
{
    cpl_size       degree0 = cpl_polynomial_get_degree(self);
    const cpl_size degree1 = cpl_polynomial_get_degree(first);
    const cpl_size degree2 = cpl_polynomial_get_degree(second);
    const cpl_size maxdeg  = degree1 > degree2 ? degree1 : degree2;


    cpl_ensure_code(self   != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(first  != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(second != NULL, CPL_ERROR_NULL_INPUT);

    cpl_ensure_code(cpl_polynomial_get_dimension(self) ==
                    cpl_polynomial_get_dimension(first),
                    CPL_ERROR_INCOMPATIBLE_INPUT);
    cpl_ensure_code(cpl_polynomial_get_dimension(self) ==
                    cpl_polynomial_get_dimension(second),
                    CPL_ERROR_INCOMPATIBLE_INPUT);

    /* FIXME: */
    cpl_ensure_code(cpl_polynomial_get_dimension(self) == 1,
                    CPL_ERROR_UNSUPPORTED_MODE);

    if (degree0 < maxdeg) {
        degree0 = maxdeg;
    } else {
        /* Reset coefficients in self as needed */
        for (; degree0 > maxdeg; degree0--) {
            cpl_polynomial_delete_coeff(self, &degree0);
        }
    }

    assert( degree0 == maxdeg );

    for (; degree0 >= 0; degree0--) {
        const double val1 = cpl_polynomial_get_coeff(first, &degree0);
        const double val2 = cpl_polynomial_get_coeff(second, &degree0);
        cpl_polynomial_set_coeff(self, &degree0, val1 - val2);
    }

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Multiply a polynomial with a scalar
  @param    self      The polynomial to hold the result
  @param    other     The polynomial to scale, may equal self
  @param    factor    The factor to multiply with
  @return   CPL_ERROR_NONE or the relevant CPL error code

  Possible CPL error code set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_UNSUPPORTED_MODE if the dimension is not 1 (FIXME)
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_multiply_scalar(cpl_polynomial * self,
                                              const cpl_polynomial * other,
                                              double factor)
{

    const cpl_size maxdeg  = cpl_polynomial_get_degree(other);
    const cpl_size zerodeg = cpl_polynomial_get_degree(self);
    cpl_size degree;

    cpl_ensure_code(self  != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(other != NULL, CPL_ERROR_NULL_INPUT);

    cpl_ensure_code(cpl_polynomial_get_dimension(self) == 1,
                    CPL_ERROR_UNSUPPORTED_MODE);
    cpl_ensure_code(cpl_polynomial_get_dimension(other) == 1,
                    CPL_ERROR_UNSUPPORTED_MODE);

    for (degree = 0; degree <= maxdeg; degree++) {
        const double val = factor * cpl_polynomial_get_coeff(other, &degree);
        cpl_polynomial_set_coeff(self, &degree, val);
    }

    /* Reset coefficients in self as needed */
    for (; degree <= zerodeg; degree++) {
        cpl_polynomial_delete_coeff(self, &zerodeg);
    }

    return CPL_ERROR_NONE;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief  Fit a 1D-polynomial to a 1D-signal in a least squares sense
  @param  self    1D-polynomial to hold the fit
  @param  x_pos   Vector of positions of the signal to fit.
  @param  values  Vector of values of the signal to fit.
  @param  mindeg  The non-negative minimum fitting degree
  @param  degree  The polynomial fitting degree, at least mindeg
  @param  symsamp True iff the x_pos values are symmetric around their mean
  @param  mse     Iff mse is not null, the mean squared error on success
  @return The fitted polynomial or NULL on error
  @see cpl_polynomial_fit_1d_create()

  symsamp is ignored if mindeg is nonzero, otherwise
  symsamp may to be set to CPL_TRUE if and only if the values in x_pos are
  known a-priori to be symmetric around their mean, e.g. (1, 2, 4, 6, 10,
  14, 16, 18, 19), but not (1, 2, 4, 6, 10, 14, 16). Setting symsamp to
  CPL_TRUE while mindeg is zero eliminates certain round-off errors.
  For higher order fitting the fitting problem known as "Runge's phenomenon"
  is minimized using the socalled "Chebyshev nodes" as sampling points.
  For Chebyshev nodes symsamp can be set to CPL_TRUE.

 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_fit_1d(cpl_polynomial * self,
                                     const cpl_vector * x_pos,
                                     const cpl_vector * values,
                                     cpl_size mindeg, cpl_size degree,
                                     cpl_boolean symsamp, double * mse)
{

    /* Number of unknowns to determine */
    const cpl_size     pdeg = cpl_polynomial_get_degree(self);
    const cpl_size     nc = 1 + degree - mindeg;
    const cpl_size     np = cpl_vector_get_size(x_pos);
    double             mean;
    double             delta, xval;
    cpl_vector       * xhat;
    const cpl_vector * xuse;
    cpl_matrix       * mh;   /* Hankel matrix */
    cpl_matrix       * mx;
    cpl_size           i, j;
    cpl_error_code     error;


    cpl_ensure_code(np > 0,      cpl_error_get_code());
    cpl_ensure_code(values,      CPL_ERROR_NULL_INPUT);

    cpl_ensure_code(cpl_vector_get_size(values) == np,
                    CPL_ERROR_INCOMPATIBLE_INPUT);

    cpl_ensure_code(mindeg >= 0,      CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(degree >= mindeg, CPL_ERROR_ILLEGAL_INPUT);

    symsamp = symsamp && (mindeg == 0); /* symsamp not usable with mindeg > 0 */

    /* Reset coefficients below mindeg */
    for (i = 0; i < mindeg; i++) {
        cpl_polynomial_delete_coeff(self, &i);
    }
    /* Reset coefficients above degree */
    for (i = degree + 1; i <= pdeg; i++) {
        cpl_polynomial_delete_coeff(self, &i);
    }

    if (degree == mindeg) {
        /* Handle this one-coefficient polynomial as a special case */

        if (degree == 0) {
            /* Handle this zero-degree polynomial as a special case */

            /* If requested, compute mean squared error */
            if (mse != NULL) *mse =
                cpl_tools_get_variance_double(cpl_vector_get_data_const(values),
                                              np, &mean);
            else
                mean = cpl_vector_get_mean(values);

        } else {
            /* A polynomial with just one coefficient and a positive degree
               is requested. The coefficient must therefore be non-zero. */

            /* Raise values to the power of mindeg */
            const double * xpos = cpl_vector_get_data_const(x_pos);
            const double * dval = cpl_vector_get_data_const(values);
            double h = 0.0; /* Hankel = Transpose(Vandermonde) * Vandermonde */
            double vtv = 0.0; /* Transpose(Vandermonde) * values */

            for (i=0; i < np; i++) {
                const double xn = cpl_tools_ipow(xpos[i], (int)mindeg);

                vtv += xn * dval[i];
                h   += xn * xn;
            }

            if (h > 0.0) {
                mean = vtv / h;
            } else {
                return cpl_error_set_message_(CPL_ERROR_DIVISION_BY_ZERO,
                                             "mindeg=%" CPL_SIZE_FORMAT ". "
                                              "degree=%" CPL_SIZE_FORMAT ". "
                                              "nc=%" CPL_SIZE_FORMAT ". "
                                              "np=%" CPL_SIZE_FORMAT ". "
                                              "Coeff = %g / %g",
                                              mindeg, degree,
                                              nc, np, vtv, h);
            }
        }
        /* Should not be able to fail now, nevertheless propagate */
        return cpl_error_set_(cpl_polynomial_set_coeff(self, &degree, mean));
    }
    cpl_ensure_code(np >= nc,   CPL_ERROR_DATA_NOT_FOUND);

    /* The Hankel matrix may be singular in such a fashion, that the pivot
       points in its Cholesky decomposition are positive due to rounding errors.
       To ensure that such singular systems are robustly detected, the number of
       distinct sampling points is counted.
    */

    cpl_ensure_code(!cpl_vector_ensure_distinct(x_pos, nc),
               CPL_ERROR_SINGULAR_MATRIX);

    /* Try to detect if the x-points are equidistant
       - in which every other skew diagonal of the Hankel matrix is zero */
    xval = cpl_vector_get(x_pos, 1);
    delta = xval - cpl_vector_get(x_pos, 0);
    for (i=1; i < np-1; i++) {
        const double dprev = delta;
        const double xprev = xval;

        xval = cpl_vector_get(x_pos, i+1);
        delta = xval - xprev;
        if (delta != dprev) break;
    }

    if (mindeg == 0) {
        /* Transform: xhat = x - mean(x) */
        xhat = cpl_vector_transform_mean(x_pos, &mean);
        xuse = xhat;
    } else {
        mean = 0.0;
        xhat = NULL;
        xuse = x_pos;
    }

    assert( xuse != NULL );

    /* Generate Hankel matrix, H = V' * V, where V is the Vandermonde matrix */
    /* FIXME: It is faster and likely more accurate to compute the QR
       factorization of the Vandermonde matrix, QR = V, see
       C.J. Demeure: Fast QR Factorization of Vandermonde Matrices */

    /* mh is initialized only if xuse is not equidistant */
    mh = symsamp ? cpl_matrix_new(nc, nc) :
        cpl_matrix_wrap(nc, nc, cpl_malloc((size_t)(nc * nc) * sizeof(double)));
    mx = cpl_matrix_wrap(nc, 1, cpl_malloc((size_t)(nc *  1) * sizeof(double)));

    cpl_matrix_fill_normal_vandermonde(mh, mx, xuse, symsamp, mindeg, values);
#ifdef CPL_POLYNOMIAL_FIT_DEBUG
    cpl_msg_warning(cpl_func, "MINDEG=%" CPL_SIZE_FORMAT ". degree=%"
                    CPL_SIZE_FORMAT ". nc=%" CPL_SIZE_FORMAT ". np=%"
                    CPL_SIZE_FORMAT ". mean=%g", mindeg, degree, nc, np, mean);
    cpl_matrix_dump(mh, stdout);
    cpl_matrix_dump(mx, stdout);
#endif

    cpl_vector_delete(xhat);

    error = cpl_matrix_solve_spd(mh, mx);

    cpl_matrix_delete(mh);

    if (error) {
        cpl_matrix_delete(mx);
        return cpl_error_set_message_(error, "mindeg=%" CPL_SIZE_FORMAT
                                      ". degree=%" CPL_SIZE_FORMAT ". nc=%"
                                      CPL_SIZE_FORMAT ". np=%" CPL_SIZE_FORMAT
                                      ". mean=%g", mindeg, degree, nc, np,
                                      mean);
    }


    /* Scale back - and store coefficients - with leading coefficient first,
       see doxygen of cpl_polynomial_set_coeff() */
    for (j = nc-1; j >= 0; j--) {
        const double coeff = cpl_matrix_get(mx, j, 0);
        const cpl_size k = j + mindeg;

        cpl_polynomial_set_coeff(self, &k, coeff);

    }

    cpl_matrix_delete(mx);

    if (mindeg == 0) {
        /* Shift back */
        cpl_polynomial_shift_1d(self, 0, -mean);
    }

    /* If requested, compute mean squared error */
    if (mse != NULL) {
        *mse = 0;
        for (i=0; i < np; i++) {
            /* Subtract from the true value, square, accumulate */
            const double residue = cpl_vector_get(values, i)
                - cpl_polynomial_eval_1d(self, cpl_vector_get(x_pos, i), NULL);
            *mse += residue * residue;
        }
        /* Average the error term */
        *mse /= (double)np;

        cpl_tools_add_flops( 3 * np + 1 );
    }

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief  Fit a 2D-polynomial to a 2D-surface in a least squares sense
  @param  xy_pos  Bivector  positions of the surface to fit.
  @param  values  Vector of values of the surface to fit.
  @param  dimdeg  True iff there is a fitting degree per dimension
  @param  maxdeg  Pointer to 1 or d maximum fitting degree(s), at least mindeg
  @param  mse     Iff mse is not null, the mean squared error on success
  @return The fitted polynomial or NULL on error
  @see cpl_polynomial_fit_2d_create
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code cpl_polynomial_fit_2d(cpl_polynomial * self,
                                            const cpl_bivector * xy_pos,
                                            const cpl_vector * values,
                                            cpl_boolean dimdeg,
                                            const cpl_size * maxdeg,
                                            double * mse)
{

    const cpl_size     np = cpl_bivector_get_size(xy_pos);
    cpl_size           degree; /* The degree of the fitted polynomial */
    /* Number of unknowns to determine */
    cpl_size           nc;
    cpl_matrix       * mv;   /* The transpose of the Vandermonde matrix */
    cpl_matrix       * mh;   /* Block-Hankel matrix, V'*V */
    cpl_matrix       * mb;
    cpl_matrix       * mx;
    const double     * coeffs1d;
    double           * dmv;
    cpl_vector       * xhat;
    cpl_vector       * yhat;
    double             xmean;
    double             ymean;
    cpl_size           powers[2];
    cpl_size           degx, degy;
    cpl_size           i, j;
    cpl_error_code     error;


    cpl_ensure_code(self   != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(cpl_polynomial_get_dimension(self) == 2,
                    CPL_ERROR_INVALID_TYPE);
    cpl_ensure_code(np > 0,         cpl_error_get_code());
    cpl_ensure_code(values != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(maxdeg != NULL, CPL_ERROR_NULL_INPUT);

    cpl_ensure_code(cpl_vector_get_size(values) == np,
                    CPL_ERROR_INCOMPATIBLE_INPUT);

    cpl_ensure_code(maxdeg[0] >= 0, CPL_ERROR_ILLEGAL_INPUT);

    if (dimdeg) {

        cpl_ensure_code(maxdeg[1] >= 0, CPL_ERROR_ILLEGAL_INPUT);

        degree = maxdeg[0] + maxdeg[1];
        nc = (maxdeg[0] + 1) * (maxdeg[1] + 1);

    } else {

        degree = maxdeg[0];
        nc = (maxdeg[0] + 1) * (maxdeg[0] + 2) / 2;

    }

    cpl_ensure_code(np >= nc, CPL_ERROR_DATA_NOT_FOUND);

    /* Reset coefficients above degree, using easy direct access. */
    for (i = 0; i < self->nc;) {
        if (dimdeg ?
            self->pow[2 * i + 0] > maxdeg[0] ||
            self->pow[2 * i + 1] > maxdeg[1]
            : self->pow[2 * i + 0] + self->pow[2 * i + 1] > degree) {
            cpl_polynomial_delete_coeff(self, self->pow + 2 * i);
        } else {
            i++;
        }
    }

    if (degree == 0) {
        /* Handle this as a special case */
        powers[0] = powers[1] = 0;
        /* Copy-paste from 1D 0-degree */

        /* If requested, compute mean squared error */
        if (mse != NULL)
            *mse = cpl_tools_get_variance_double(cpl_vector_get_data_const
                                                 (values), np, &xmean);
        else
            xmean = cpl_vector_get_mean(values);

        return cpl_polynomial_set_coeff(self, powers, xmean);
    }

    /* Transform: xhat = x - mean(x) */
    xhat = cpl_vector_transform_mean(cpl_bivector_get_x_const(xy_pos), &xmean);
    assert( xhat != NULL );

    /* Transform: yhat = y - mean(y) */
    yhat = cpl_vector_transform_mean(cpl_bivector_get_y_const(xy_pos), &ymean);
    assert( yhat != NULL );

    /* Initialize matrices */
    /* mv contains the polynomial terms in the order described */
    /* above in each row, for each input point. */
    dmv = (double*)cpl_malloc((size_t)(nc * np) * sizeof(double));
    mv = cpl_matrix_wrap(nc, np, dmv);

    /* Has redundant FLOPs, appears to improve accuracy */
    for (i=0; i < np; i++) {
        const double x = cpl_vector_get(xhat, i);
        const double y = cpl_vector_get(yhat, i);
        double yvalue = 1.0;
        j = 0;
        for (degy = 0; degy <= (dimdeg ? maxdeg[1] : degree); degy++) {
            double xvalue = 1.0;
            for (degx = 0; degx <= (dimdeg ? maxdeg[0] : degree - degy);
                 degx++, j++) {
                dmv[np * j + i] = xvalue * yvalue;
                xvalue *= x;
            }
            yvalue *= y;
        }
        assert( j == nc );
    }
    cpl_tools_add_flops( np * (nc * 2 + 1 + (dimdeg ? maxdeg[1] : degree)));

    cpl_vector_delete(xhat);
    cpl_vector_delete(yhat);

    /* mb contains the values, it is _not_ modified */
CPL_DIAG_PRAGMA_PUSH_IGN(-Wcast-qual);
    mb = cpl_matrix_wrap(np, 1, cpl_vector_get_data((cpl_vector*)values));
CPL_DIAG_PRAGMA_POP;

    /* Form the right hand side of the normal equations */
    mx = cpl_matrix_product_create(mv, mb);

    (void)cpl_matrix_unwrap(mb);

    /* Form the matrix of the normal equations */
    mh = cpl_matrix_product_normal_create(mv);
    cpl_matrix_delete(mv);

    /* Solve XA=B by a least-square solution (aka pseudo-inverse). */
    error = cpl_matrix_solve_spd(mh, mx);

    cpl_matrix_delete(mh);

    if (error) {
        cpl_matrix_delete(mx);
        return cpl_error_set_(error);
    }

    /* Store coefficients for output */

    coeffs1d = cpl_matrix_get_data(mx);

    for (j = 0; j < nc; j++) {
        if (isnan(coeffs1d[j])) {
            cpl_matrix_delete(mx);
            return cpl_error_set_(CPL_ERROR_DIVISION_BY_ZERO);
        }
    }

    j = 0;
    for (degy = 0; degy <= (dimdeg ? maxdeg[1] : degree); degy++) {
        powers[1] = degy;
        for (degx = 0; degx <= (dimdeg ? maxdeg[0] : degree - degy);
             degx++, j++) {
            powers[0] = degx;
            if (coeffs1d[j] != cpl_matrix_get(mx, j, 0)) break;
            cpl_polynomial_set_coeff(self, powers, cpl_matrix_get(mx, j, 0));
        }
        if (degx <= (dimdeg ? maxdeg[0] : degree - degy)) break;
    }
    cpl_matrix_delete(mx);

    cpl_ensure_code(j == nc, CPL_ERROR_UNSPECIFIED);


    /* Transform the polynomial back */
    cpl_polynomial_shift_1d(self, 0, -xmean);
    cpl_polynomial_shift_1d(self, 1, -ymean);

    /* If requested, compute mean squared error */
    if (mse != NULL) {
        const cpl_vector * x_pos = cpl_bivector_get_x_const(xy_pos);
        const cpl_vector * y_pos = cpl_bivector_get_y_const(xy_pos);
        cpl_vector * x_val = cpl_vector_new(2);

        *mse = 0;
        for (i=0; i<np; i++) {
            double residue;
            cpl_vector_set(x_val, 0, cpl_vector_get(x_pos, i));
            cpl_vector_set(x_val, 1, cpl_vector_get(y_pos, i));
            /* Subtract from the true value, square, accumulate */
            residue = cpl_vector_get(values, i)
                - cpl_polynomial_eval(self, x_val);
            *mse += residue * residue;
        }
        cpl_vector_delete(x_val);
        /* Average the error term */
        *mse /= (double)np;
        cpl_tools_add_flops( 3 * np + 1 );
    }

    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Delete (set to zero) a coefficient of a polynomial
  @param    in      The polynomial to modify
  @param    pows    The power(s) of the different variables
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_
  @see cpl_polynomial_set_coeff()
  @note Passing a row from self->pow as pows causes memory aliasing. This is
        allowed, but pows[] will be modified (or even unallocated) by the call.

  Possible #_cpl_error_code_ set in this function:
  - CPL_ERROR_NULL_INPUT if an input pointer is NULL
  - CPL_ERROR_ILLEGAL_INPUT if pows contains negative values
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code cpl_polynomial_delete_coeff(cpl_polynomial * self,
                                                  const cpl_size * pows)
{

    cpl_size i, dim;
    cpl_size ind;


    cpl_ensure_code(self != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(pows != NULL, CPL_ERROR_NULL_INPUT);
    for (dim=0; dim < self->dim; dim++)
        cpl_ensure_code(pows[dim] >= 0, CPL_ERROR_ILLEGAL_INPUT);

    if (self->dim == 1) {
        /* Handle 1D as a special case */
        if (pows[0] == self->nc-1) {
            /* Leading 1D-coefficient is zero - find the new leading */
            do {
                self->nc--;
            } while (self->nc > 0 && self->c[self->nc - 1] == 0.0);

            if (self->nc > 0)
                self->c = cpl_realloc(self->c,
                                      (size_t)self->nc * sizeof(double));
            else
                cpl_free(self->c);

        } else if (pows[0] < self->nc-1) {
            self->c[pows[0]] = 0.0;
        }
        return CPL_ERROR_NONE;
    }

    /* Find the coeff of the multi-variate polynomial */
    for (i=0; i < self->nc; i++) {
        if (!memcmp(self->pow + (size_t)(self->dim * i), pows,
                    (size_t)self->dim * sizeof(cpl_size)))
            break; /* Found the right combination of powers */
    }
    if (i == self->nc) return CPL_ERROR_NONE;

    /* The coefficient exists: Set it to zero and return */
    ind = i;

    /* Shrink polynomium */
    self->nc--;

    /* Shrink array of powers */
    if (self->nc > 0) {

        /* self->nc now points to the last coefficient in the polynomial */

        /* Reduce (if possible) the length of the power-table */
        for (dim = 0; dim < self->dim; dim++) {
            cpl_size new_max;

            if (self->pow[self->dim * ind + dim] < self->max_degree[dim])
                continue;

            /* The max-power of this dim may decrease */

            new_max = 0;
            for (i=0; i < self->nc && new_max < self->max_degree[dim]; i++) {
                if (self->pow[self->dim * i + dim] > new_max)
                    new_max = self->pow[self->dim * i + dim];
            }

            if (new_max == self->max_degree[dim]) continue;

            /* The max-power of this dim decreases */
            self->max_degree_alldims -= self->max_degree[dim] - new_max;
            self->max_degree[dim] = new_max;
            self->eval_pow[dim] =
                cpl_realloc(self->eval_pow[dim],
                            (size_t)(1+self->max_degree[dim]) * sizeof(double));
        }

        if (ind < self->nc) {
            /* Move last coefficient to place of zero-valued one */
            self->c[ind] = self->c[self->nc];

            /* Copy last powers to place of those of the zero-valued one */
            memcpy(self->pow + (size_t)(self->dim * ind), self->pow
                   + (size_t)(self->dim * self->nc),
                   (size_t)self->dim * sizeof(cpl_size));

        }

        self->c = cpl_realloc(self->c, (size_t)self->nc * sizeof(double));
        self->pow = cpl_realloc(self->pow, (size_t)self->dim
                                * (size_t)self->nc * sizeof(cpl_size));

#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
        cpl_free(self->added_horner);
        cpl_free(self->added_exist_iterations);
        self->added_horner           = 
            cpl_malloc(self->nc * sizeof(*self->added_horner));
        self->added_exist_iterations = 
            cpl_malloc(self->nc * sizeof(*self->added_exist_iterations));
#endif
    } else {
        cpl_free(self->c);
        cpl_free(self->pow);
        cpl_free(self->max_degree);
#ifdef CPL_POLYNOMIAL_USE_MULTI_HORNER
        cpl_free(self->added_horner);
        cpl_free(self->added_exist_iterations);
#endif
        self->max_degree_alldims = 0;
        for (dim=0; dim < self->dim; dim++) {
            cpl_free(self->eval_pow[dim]);
        }
        cpl_free(self->eval_pow);
    }

    return CPL_ERROR_NONE;

}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Fill the Hankel Matrix H=V'*V, where V is a 1D-Vandermonde matrix
  @param    self      The matrix H
  @param    mx        A right multiplication with V', mx = V' * values
  @param    xhat      The mean-transformed x-values
  @param    is_eqdist True iff xhat contains equidistant points
  @param    mindeg    The non-negative minimum fitting degree
  @param    values    The values to be interpolated
  @return   void
  @note self must have its elements initialized to zero iff is_eqdist is true.

 */
/*----------------------------------------------------------------------------*/
static void cpl_matrix_fill_normal_vandermonde(cpl_matrix * self,
                                               cpl_matrix * mx,
                                               const cpl_vector * xhat,
                                               cpl_boolean is_eqdist,
                                               cpl_size mindeg,
                                               const cpl_vector * values)
{


    const double * dval = cpl_vector_get_data_const(values);
    const double * xval = cpl_vector_get_data_const(xhat);
    cpl_vector   * phat = cpl_vector_duplicate(xhat); /* Powers of xhat */
    cpl_vector   * qhat = NULL;                       /* mindeg Power of xhat */
    double       * dhat = cpl_vector_get_data(phat);
    double       * ehat = NULL;
    const cpl_size nc   = cpl_matrix_get_ncol(self);
    const cpl_size np   = cpl_vector_get_size(xhat);
    cpl_size       i,j;


    assert( nc == cpl_matrix_get_nrow(self) );
    assert( nc == cpl_matrix_get_nrow(mx) );
    assert( 1  == cpl_matrix_get_ncol(mx) );
    assert( np == cpl_vector_get_size(values) );


    /* Fill Hankel matrix from top-left to main skew diagonal
       - on and above (non-skew) main diagonal */
    /* Also compute transpose(V) * b */
    /* Peel off 1st iteration */
    if (mindeg > 0) {
        double hsum = 0.0;
        cpl_size k;

        qhat = mindeg == 1 ? cpl_vector_duplicate(xhat) : cpl_vector_new(np);
        ehat = cpl_vector_get_data(qhat);

        /* Raise xhat to the power of mindeg */
        for (k=0; k < np; k++) {
            const double x = xval[k];

            if (mindeg > 1) ehat[k] = cpl_tools_ipow(x, (int)mindeg);
            dhat[k] *= ehat[k];

            hsum += ehat[k] * ehat[k];
        }
        cpl_matrix_set(self, 0, 0, hsum);
    } else {
        cpl_matrix_set(self, 0, 0, (double)np);
    }
    /* qhat is xhat to the power of mindeg, iff mindeg > 0 */
    /* dhat is xhat to the power of 1+mindeg, iff mindeg > 0 */
    for (j=1; j < 2; j++) {
        double vsum0 = 0.0;
        double hsum = 0.0;
        double vsum = 0.0;
        cpl_size k;

        for (k=0; k < np; k++) {
            const double y = dval[k];

            hsum += mindeg > 0 ? ehat[k] * dhat[k] : dhat[k];
            vsum += y * dhat[k];
            vsum0 += mindeg > 0 ? ehat[k] * y : y;
        }
        cpl_matrix_set(mx, 0, 0, vsum0);
        cpl_matrix_set(mx, j, 0, vsum);
        if (is_eqdist) continue;
        k = j;
        for (i=0; i <= k; i++, k--) {
            cpl_matrix_set(self, i, k, hsum);
        }
    }
    for (; j < nc; j++) {
        double   hsum = 0.0;
        double   vsum = 0.0;
        cpl_size k;

        for (k=0; k < np; k++) {
            const double x = xval[k];
            const double y = dval[k];

            dhat[k] *= x;
            hsum += mindeg > 0 ? ehat[k] * dhat[k] : dhat[k];
            vsum += y * dhat[k];
        }
        cpl_matrix_set(mx, j, 0, vsum);
        if (is_eqdist && (j&1)) continue;
        k = j;
        for (i=0; i <= k; i++, k--) {
            cpl_matrix_set(self, i, k, hsum);
        }
    }
    /* Fill remaining Hankel matrix - on and above (non-skew) main diagonal */

    if (mindeg > 0) cpl_vector_multiply(phat, qhat);
    cpl_vector_delete(qhat);
    for (i = 1; i < nc; i++) {
        cpl_size k;
        double   hsum = 0.0;

        if (is_eqdist && ((i+nc)&1)==0) {
            cpl_vector_multiply(phat, xhat);
            continue;
        }

        for (k=0; k < np; k++) {
            const double x = xval[k];

            dhat[k] *= x;
            hsum += dhat[k];
        }
        k = i;
        for (j = nc-1; k <= j; k++, j--) {
            cpl_matrix_set(self, k, j, hsum);
        }
    }

    cpl_tools_add_flops( 6 * np * ( nc - 1) );

    cpl_vector_delete(phat);

}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Try to detect if the x-points are equidistant
  @param    self  The points to check
  @return   1 if equidistant, 0 if not, -1 on error.
  @note If yes, then every other skew diagonal of the Hankel matrix is zero,
        H = V' * V, where V is the 1D-Vandermonde matrix from the x-points.

 */
/*----------------------------------------------------------------------------*/
static int cpl_vector_is_eqdist(const cpl_vector * self)
{

    double xval, delta;
    const cpl_size np = cpl_vector_get_size(self);
    cpl_size i;

    cpl_ensure(self, CPL_ERROR_NULL_INPUT, -1);

    if (cpl_vector_get_size(self) == 1) return 1;

    xval = cpl_vector_get(self, 1);
    delta = xval - cpl_vector_get(self, 0);
    for (i=1; i < np-1; i++) {
        const double dprev = delta;
        const double xprev = xval;

        xval = cpl_vector_get(self, i+1);
        delta = xval - xprev;
        if (delta != dprev) break;
    }

    return i == np-1 ? 1 : 0;

}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Given p and u, modify the polynomial to p(x) := p(x+u)
  @param    p  The polynomial coefficients to be modified in place
  @param    n  The number of coefficients
  @param    u  The shift
  @return   void
  @see      cpl_polynomial_shift_1d
  @note     The function will seg-fault on NULL input.

 */
/*----------------------------------------------------------------------------*/
void cpl_polynomial_shift_double(double * coeffs, cpl_size n, double u)
{
    for (size_t j = 0; j < (size_t)n-1; j++)
        for (size_t i = 1; i < (size_t)n - j; i++ )
            coeffs[n-1-i] += coeffs[n-i] * u;

    cpl_tools_add_flops( n * ( n - 1) );
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    A real solution to p(x) = 0 using Newton-Raphsons method
  @param    p    The 1D-polynomial
  @param    x0   First guess of the solution
  @param    px   The solution, on error see above
  @param    mul  The root multiplicity (or 1 if unknown)
  @param    bpos Iff CPL_TRUE, then fail if a derivative is non-positive
  @return   CPL_ERROR_NONE or the relevant #_cpl_error_code_
  @see cpl_polynomial_solve_1d()

 */
/*----------------------------------------------------------------------------*/
cpl_error_code cpl_polynomial_solve_1d_(const cpl_polynomial * p,
                                        double x0, double * px, cpl_size mul,
                                        cpl_boolean bpos)
{

    /* Initialize to ensure at least one iteration */
    double       r = 1;
    double       d = 0;
    double       xprev = 2 * x0 + 1;
    const double mm = (double)mul; /* Root multiplicity */
    cpl_size     mite;
    cpl_size     i;


    cpl_ensure_code(px != NULL,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(p  != NULL,  CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(p->dim == 1, CPL_ERROR_INVALID_TYPE);
    cpl_ensure_code(mul > 0,     CPL_ERROR_ILLEGAL_INPUT);

    /* Iterating towards zero is not as simple as it sounds, so don't */
    if (p->nc == 0 || p->c[0] == 0.0) {
       *px = 0.0;
       return CPL_ERROR_NONE;
    }

    *px = x0;

    mite = p->nc * CPL_NR_MAXITE;

    for (i = 0; i < mite; i++, xprev = *px) {

        const double rprev = r;
        const double dprev = d;

        /* Compute residual, r = p(x) and derivative, d = p`(x) */
        r = cpl_polynomial_eval_1d(p, *px, &d);

        /* Stop if:
           0) If bpos is true and d <= 0.0. This indicates a failure to solve
           1) Correction did not decrease - unless p`(x) changed sign
           2) p`(x) == 0. It is insufficient to implement this as d == 0,
              because some non-zero divisors can still result in inf. */
        if (bpos && d <= 0.0) break;
        if (d * dprev >= 0.0 && fabs(r * dprev) >= fabs(rprev * d)) break;

        /* Compute and apply the accelerated Newton-Raphson correction */
        *px -= mm * r / d;

        /* *px can become NaN - in which case the iteration goes on
           until the maximum number of iterations is reached.
           In one case this happens because a p'(x) == 0 is pertubed by
           round-off. */

        /* Stop also if:
           3) The correction did not change the solution
              - will typically save at most one Horner-evaluation  */
        if (fabs(*px-xprev) < fabs(*px) * DBL_EPSILON) break;

        /* if x_i == x_j for i > j > 0 the iteration cannot converge.
           This can only happen with at least two sign-changes for p`(x_k)
           i >= k >= j. This is not checked for. */

    }

    cpl_tools_add_flops( i * 12 );

    if (i == mite) return
        cpl_error_set_message_(CPL_ERROR_CONTINUE, "x0=%g, mul=%"
                               CPL_SIZE_FORMAT ", degree=%" CPL_SIZE_FORMAT
                               ", i=%d, p(x=%g)=%g", x0, mul, p->nc-1, (int)i,
                               *px, cpl_polynomial_eval_1d(p, *px, NULL));

    /* At this point:
       In absence of rounding r or d is zero.
       Due to rounding r or d is zero or close to zero.
       If there is no solution only d is (close to) zero.
       If there is a single solution only r is (close to) zero.
       If there is a multiple solution both r and d are (close to) zero
       - in this case |r| cannot be bigger than |d| because p is one
       degree higher than p'. */

    if (bpos && d <= 0.0) return
        cpl_error_set_message_(CPL_ERROR_ILLEGAL_INPUT, "x0=%g, mul=%"
                               CPL_SIZE_FORMAT ", degree=%" CPL_SIZE_FORMAT
                               ", i=%d, p(x=%g)=%g, p'(x)=%g <= 0.0",
                               x0, mul, p->nc-1, (int)i, *px, r, d);

    if (fabs(r) > fabs(d)) {
        /* When d is computed in long double precision at a multiple root,
           |r| _can_ be bigger than |d|.

           Quick fix: Assume that solution is still OK, if |r| is small
           compared to the largest coefficient */

        /* Since r is non-zero, at least one coefficient must be non-zero */
        cpl_size n = p->nc;
        double max = 0.0;
        while (n--) if (fabs(p->c[n]) > max) max = fabs(p->c[n]);

        if (fabs(r) > max * DBL_EPSILON) return
            cpl_error_set_message_(CPL_ERROR_DIVISION_BY_ZERO, "x0=%g, mul=%"
                                   CPL_SIZE_FORMAT ", degree=%" CPL_SIZE_FORMAT
                                   ", i=%d, p(x=%g)=%g, p'(x)=%g",
                                   x0, mul, p->nc-1, (int)i, *px, r, d);
    }

    return CPL_ERROR_NONE;

}
