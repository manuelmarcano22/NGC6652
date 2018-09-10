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

#ifndef CPL_VECTOR_FIT_IMPL_H
#define CPL_VECTOR_FIT_IMPL_H


/*
 * FIXME: The code in this file is a copy of the cpl_fit module and has to
 *        stay here until cpl_vector_fit_gaussian() is moved to the cpl_fit
 *        module, which can be done only before the release of the next major
 *        version, in order not to break the library hierarchy!
 *
 *        When the code in this file is finally moved to the cpl_fit module,
 *        this file has to be removed!
 */

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include "cpl_fit.h"

#include "cpl_vector.h"
#include "cpl_matrix.h"
#include "cpl_memory.h"
#include "cpl_error_impl.h"
#include "cpl_errorstate.h"


#include <assert.h>
#include <math.h>
/*----------------------------------------------------------------------------*/
/*
   @brief   Get new position in parameter space (L-M algorithm)
   @param   a       Current fit parameters.
   @param   ia      Non-NULL array defining with non-zero values which
                    parameters participate in the fit.
   @param   M       Number of fit parameters
   @param   N       Number of positions
   @param   D       Dimension of x-positions
   @param   lambda  Lambda in L-M algorithm.
   @param   f       Function that evaluates the fit function.
   @param   dfda    Function that evaluates the partial derivaties
                    of the fit function w.r.t. fit parameters.
   @param   x       The input positions (pointer to MxD matrix buffer).
   @param   y       The N values to fit.
   @param   sigma   A vector of size N containing the uncertainties of the
                    y-values. If NULL, a constant uncertainty equal to 1 is
            assumed.
   @param   partials The partial derivatives (work space).
   @param   alpha   Alpha in L-M algorithm (work space).
   @param   beta    Beta in L-M algorithm (work space).
   @param   a_da    (output) Candidate position in parameter space.

   @return  0 iff okay.

   This function computes a potentially better set of parameters @em a + @em da,
   where @em da solves the equation @em alpha(@em lambda) * @em da = @em beta .

   Possible #_cpl_error_code_ set in this function:
   - CPL_ERROR_ILLEGAL_INPUT if the fit function or its derivative could
   not be evaluated.
   - CPL_ERROR_SINGULAR_MATRIX if @em alpha is singular.

*/
/*----------------------------------------------------------------------------*/
inline static int
get_candidate(const double *a, const int ia[],
          cpl_size M, cpl_size N, cpl_size D,
          double lambda,
          int    (*f)(const double x[], const double a[], double *result),
          int (*dfda)(const double x[], const double a[], double result[]),
          const double *x,
          const double *y,
          const double *sigma,
          double *partials,
          cpl_matrix *alpha,
          cpl_matrix *beta,
          double *a_da)
{
    cpl_size Mfit = 0;    /* Number of non-constant fit parameters */
    cpl_matrix *da;       /* Solution of   alpha * da = beta */
    double *alpha_data;
    double *beta_data;
    double *da_data;
    cpl_size i, j;
    int imfit = 0;
    int jmfit = 0;
    int k = 0;

    /* For efficiency, don't check input in this static function */

    Mfit = cpl_matrix_get_nrow(alpha);

    alpha_data    = cpl_matrix_get_data(alpha);
    beta_data     = cpl_matrix_get_data(beta);

    /* Build alpha, beta:
     *
     *  alpha[i,j] = sum_{k=1,N} (sigma_k)^-2 * df/da_i * df/da_j  *
     *                           (1 + delta_ij lambda) ,
     *
     *   beta[i]   = sum_{k=1,N} (sigma_k)^-2 * ( y_k - f(x_k) ) * df/da_i
     *
     * where (i,j) loop over the non-constant parameters (0 to Mfit-1),
     * delta is Kronecker's delta, and all df/da are evaluated in x_k
     */

    cpl_matrix_fill(alpha, 0.0);
    cpl_matrix_fill(beta , 0.0);

    for (k = 0; k < N; k++)
    {
        double sm2 = 0.0;                /* (sigma_k)^-2 */
        double fx_k = 0.0;               /* f(x_k)       */
        const double *x_k = &(x[0+k*D]); /* x_k          */

        if (sigma == NULL)
        {
            sm2 = 1.0;
        }
        else
        {
            sm2 = 1.0 / (sigma[k] * sigma[k]);
        }

        /* Evaluate f(x_k) */
        cpl_ensure( f(x_k, a, &fx_k) == 0, CPL_ERROR_ILLEGAL_INPUT, -1);

        /* Evaluate (all) df/da (x_k) */
        cpl_ensure( dfda(x_k, a, partials) == 0,
            CPL_ERROR_ILLEGAL_INPUT, -1);

        for (i = 0, imfit = 0; i < M; i++)
        {
            if (ia[i] != 0)
            {
                /* Beta */
                beta_data[imfit] +=
                sm2 * (y[k] - fx_k) * partials[i];

                /* Alpha is symmetrical, so compute
                   only lower-left part */
                for (j = 0, jmfit = 0; j < i; j++)
                {
                    if (ia[j] != 0)
                    {
                        alpha_data[jmfit + imfit*Mfit] +=
                        sm2 * partials[i] *
                        partials[j];

                        jmfit += 1;
                    }
                }

                /* Alpha, diagonal terms */
                j = i;
                jmfit = imfit;

                alpha_data[jmfit + imfit*Mfit] +=
                sm2 * partials[i] *
                partials[j] * (1 + lambda);

                imfit += 1;
            }
        }

        assert( imfit == Mfit );
    }

    /* Create upper-right part of alpha */
    for (i = 0, imfit = 0; i < M; i++)
    {
        if (ia[i] != 0)
        {
            for (j = i+1, jmfit = imfit+1; j < M; j++)
            {
                if (ia[j] != 0)
                {
                    alpha_data[jmfit + imfit*Mfit] =
                    alpha_data[imfit + jmfit*Mfit];

                    jmfit += 1;
                }
            }
            assert( jmfit == Mfit );

            imfit += 1;
        }
    }
    assert( imfit == Mfit );

    da = cpl_matrix_solve(alpha, beta);

    cpl_ensure(da != NULL, cpl_error_get_code(), -1);

    /* Create a+da vector by adding a and da */
    da_data   = cpl_matrix_get_data(da);

    for (i = 0, imfit = 0; i < M; i++)
    {
        if (ia[i] != 0)
        {
            a_da[i] = a[i] + da_data[0 + imfit*1];

            imfit += 1;
        }
        else
        {
            a_da[i] = a[i];
        }
    }

    assert( imfit == Mfit );

    cpl_matrix_delete(da);

    return 0;
}

/*----------------------------------------------------------------------------*/
/*
   @brief   Compute chi square
   @param   N       Number of positions
   @param   D       Dimension of x-positions
   @param   f       Function that evaluates the fit function.
   @param   a       The fit parameters.
   @param   x       Where to evaluate the fit function (N x D matrix).
   @param   y       The N values to fit.
   @param   sigma   A vector of size N containing the uncertainties of the
                    y-values. If NULL, a constant uncertainty equal to 1 is
            assumed.

   @return  chi square, or a negative number on error.

   This function calculates chi square defined as
   sum_i (y_i - f(x_i, a))^2/sigma_i^2

   Possible #_cpl_error_code_ set in this function:
   - CPL_ERROR_ILLEGAL_INPUT if the fit function could not be evaluated
*/
/*----------------------------------------------------------------------------*/

inline static double
get_chisq(cpl_size N, cpl_size D,
      int (*f)(const double x[], const double a[], double *result),
      const double *a,
      const double *x,
      const double *y,
      const double *sigma)
{
    double chi_sq = 0.0;     /* Result */
    cpl_size i;

    /* For efficiency, don't check input in this static function */

    for (i = 0; i < N; i++)
    {
        double fx_i;
        double residual;             /* Residual in units of uncertainty */
        const double *x_i = &(x[0+i*D]);

        /* Evaluate */
        cpl_ensure( f(x_i, a, &fx_i) == 0, CPL_ERROR_ILLEGAL_INPUT, -1.0);

        /* Accumulate */
        if (sigma == NULL)
        {
            residual = (fx_i - y[i]);
        }
        else
        {
            residual = (fx_i - y[i]) / sigma[i];
        }

        chi_sq += residual*residual;

    }

    return chi_sq;
}


/*----------------------------------------------------------------------------*/
/*
   @brief   Fit a function to a set of data
   @param   x        N x D matrix of the positions to fit.
                     Each matrix row is a D-dimensional position.
   @param   sigma_x  Uncertainty (one sigma, gaussian errors assumed)
                     assosiated with @em x. Taking into account the
             uncertainty of the independent variable is currently
             unsupported, and this parameter must therefore be set
             to NULL.
   @param   y        The N values to fit.
   @param   sigma_y  Vector of size N containing the uncertainties of
                     the y-values. If this parameter is NULL, constant
             uncertainties are assumed.
   @param   a        Vector containing M fit parameters. Must contain
                     a guess solution on input and contains the best
             fit parameters on output.
   @param   ia       Array of size M defining which fit parameters participate
                     in the fit (non-zero) and which fit parameters are held
             constant (zero). At least one element must be non-zero.
             Alternatively, pass NULL to fit all parameters.
   @param   f        Function that evaluates the fit function
                     at the position specified by the first argument (an array
             of size D) using the fit parameters specified by the second
             argument (an array of size M). The result must be output
             using the third parameter, and the function must return zero
             iff the evaluation succeded.
   @param   dfda     Function that evaluates the first order partial
                     derivatives of the fit function with respect to the fit
             parameters at the position specified by the first argument
             (an array of size D) using the parameters specified by the
             second argument (an array of size M). The result must
             be output using the third parameter (array of size M), and
             the function must return zero iff the evaluation succeded.
   @param relative_tolerance
                     The algorithm converges by definition if the relative
                     decrease in chi squared is less than @em tolerance
                     @em tolerance_count times in a row. Recommended default:
                     CPL_FIT_LVMQ_TOLERANCE
   @param tolerance_count
                     The algorithm converges by definition if the relative
                     decrease in chi squared is less than @em tolerance
                     @em tolerance_count times in a row. Recommended default:
                     CPL_FIT_LVMQ_COUNT
   @param max_iterations
                     If this number of iterations is reached without
                     convergence, the algorithm diverges, by definition.
                     Recommended default: CPL_FIT_LVMQ_MAXITER
   @param mse        If non-NULL, the mean squared error of the best fit is
                     computed.
   @param red_chisq  If non-NULL, the reduced chi square of the best fit is
                     computed. This requires @em sigma_y to be specified.
   @param covariance If non-NULL, the formal covariance matrix of the best
                     fit parameters is computed (or NULL on error). On success
             the diagonal terms of the covariance matrix are guaranteed
             to be positive. However, terms that involve a constant
             parameter (as defined by the input array @em ia) are
             always set to zero. Computation of the covariacne matrix
             requires @em sigma_y to be specified.


   @return  CPL_ERROR_NONE iff OK.

   This function makes a minimum chi squared fit of the specified function
   to the specified data set using a Levenberg-Marquardt algorithm.

   Possible #_cpl_error_code_ set in this function:
   - CPL_ERROR_NULL_INPUT if an input pointer other than @em sigma_x, @em
     sigma_y, @em mse, @em red_chisq or @em covariance is NULL.
   - CPL_ERROR_ILLEGAL_INPUT if an input matrix/vector is empty, if @em ia
     contains only zero values, if any of @em relative_tolerance,
     @em tolerance_count or max_iterations @em is non-positive, if N <= M
     and @em red_chisq is non-NULL, if any element of @em sigma_x or @em sigma_y
     is non-positive, or if evaluation of the fit function or its derivative
     failed.
   - CPL_ERROR_INCOMPATIBLE_INPUT if the dimensions of the input
     vectors/matrices do not match, or if chi square or covariance computation
     is requested and @em sigma_y is NULL.
   - CPL_ERROR_ILLEGAL_OUTPUT if memory allocation failed.
   - CPL_ERROR_CONTINUE if the Levenberg-Marquardt algorithm failed to converge.
   - CPL_ERROR_SINGULAR_MATRIX if the covariance matrix could not be computed.

*/
/*----------------------------------------------------------------------------*/
inline static cpl_error_code
cpl_fit_lvmq_(const cpl_matrix *x, const cpl_matrix *sigma_x,
         const cpl_vector *y, const cpl_vector *sigma_y,
         cpl_vector *a, const int ia[],
         int    (*f)(const double x[], const double a[], double *result),
         int (*dfda)(const double x[], const double a[], double result[]),
         double relative_tolerance,
         int tolerance_count,
         int max_iterations,
         double *mse,
         double *red_chisq,
         cpl_matrix **covariance)
{
    const double *x_data     = NULL; /* Pointer to input data                 */
    const double *y_data     = NULL; /* Pointer to input data                 */
    const double *sigma_data = NULL; /* Pointer to input data                 */
    cpl_size N    = 0;               /* Number of data points                 */
    cpl_size D    = 0;               /* Dimension of x-points                 */
    cpl_size M = 0;                  /* Number of fit parameters              */
    cpl_size Mfit = 0;               /* Number of non-constant fit
                                        parameters                            */

    double lambda    = 0.0;          /* Lambda in L-M algorithm               */
    double MAXLAMBDA = 10e40;        /* Parameter to control the graceful exit
                    if steepest descent unexpectedly fails */
    double chi_sq    = 0.0;          /* Current  chi^2                        */
    int count        = 0;            /* Number of successive small improvements
                    in chi^2 */
    int iterations   = 0;

    cpl_matrix *alpha  = NULL;       /* The MxM ~curvature matrix used in L-M */
    cpl_matrix *beta   = NULL;       /* Mx1 matrix = -.5 grad(chi^2)          */
    double *a_data     = NULL;       /* Parameters, a                         */
    double *a_da       = NULL;       /* Candidate position a+da               */
    double *part       = NULL;       /* The partial derivatives df/da         */
    int *ia_local      = NULL;       /* non-NULL version of ia                */

    /* If covariance computation is requested, then either
     * return the covariance matrix or return NULL.
     */
    if (covariance != NULL) *covariance = NULL;

    /* Validate input */
    cpl_ensure_code(x       != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(sigma_x == NULL, CPL_ERROR_UNSUPPORTED_MODE);
    cpl_ensure_code(y       != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(a       != NULL, CPL_ERROR_NULL_INPUT);
    /* ia may be NULL */
    cpl_ensure_code(f       != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(dfda    != NULL, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(relative_tolerance > 0, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(tolerance_count    > 0, CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(max_iterations     > 0, CPL_ERROR_ILLEGAL_INPUT);

    /* Chi^2 and covariance computations require sigmas to be known */
    cpl_ensure_code( sigma_y != NULL ||
                     (red_chisq == NULL && covariance == NULL),
                     CPL_ERROR_INCOMPATIBLE_INPUT);

    D = cpl_matrix_get_ncol(x);
    N = cpl_matrix_get_nrow(x);
    M = cpl_vector_get_size(a);
    cpl_ensure_code(N > 0 && D > 0 && M > 0, CPL_ERROR_ILLEGAL_INPUT);

    cpl_ensure_code( cpl_vector_get_size(y) == N,
             CPL_ERROR_INCOMPATIBLE_INPUT);

    x_data = cpl_matrix_get_data_const(x);
    y_data = cpl_vector_get_data_const(y);
    a_data = cpl_vector_get_data(a);

    if (sigma_y != NULL)
    {
        cpl_ensure_code( cpl_vector_get_size(sigma_y) == N,
                 CPL_ERROR_INCOMPATIBLE_INPUT);
        /* Sigmas must be positive */
        cpl_ensure_code( cpl_vector_get_min (sigma_y) > 0,
                 CPL_ERROR_ILLEGAL_INPUT);
        sigma_data = cpl_vector_get_data_const(sigma_y);
    }

    ia_local = cpl_malloc((size_t)M * sizeof(int));
    cpl_ensure_code(ia_local != NULL, CPL_ERROR_ILLEGAL_OUTPUT);

    /* Count non-constant fit parameters, copy ia */
    if (ia != NULL)
    {
        cpl_size i;

        Mfit = 0;
        for (i = 0; i < M; i++)
        {
            ia_local[i] = ia[i];

            if (ia[i] != 0)
            {
                Mfit += 1;
            }
        }

        if (! (Mfit > 0))
        {
            cpl_free(ia_local);
            return cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
        }
    }
    else
    {
        /* All parameters participate */
        cpl_size i;

        Mfit = M;

        for (i = 0; i < M; i++)
        {
            ia_local[i] = 1;
        }
    }

    /* To compute reduced chi^2, we need N > Mfit */
    if (! ( red_chisq == NULL || N > Mfit ) )
    {
        cpl_free(ia_local);
        return cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
    }

    /* Create alpha, beta, a_da, part  work space */
    alpha = cpl_matrix_new(Mfit, Mfit);
    if (alpha == NULL)
    {
        cpl_free(ia_local);
        return cpl_error_set_(CPL_ERROR_ILLEGAL_OUTPUT);
    }

    beta = cpl_matrix_new(Mfit, 1);
    if (beta == NULL)
    {
        cpl_free(ia_local);
        cpl_matrix_delete(alpha);
        return cpl_error_set_(CPL_ERROR_ILLEGAL_OUTPUT);
    }

    a_da = cpl_malloc((size_t)M * sizeof(double));
    if (a_da == NULL)
    {
        cpl_free(ia_local);
        cpl_matrix_delete(alpha);
        cpl_matrix_delete(beta);
        return cpl_error_set_(CPL_ERROR_ILLEGAL_OUTPUT);
    }

    part = cpl_malloc((size_t)M * sizeof(double));
    if (part == NULL)
    {
        cpl_free(ia_local);
        cpl_matrix_delete(alpha);
        cpl_matrix_delete(beta);
        cpl_free(a_da);
        return cpl_error_set_(CPL_ERROR_ILLEGAL_OUTPUT);
    }

    /* Initialize loop variables */
    lambda = 0.001;
    count = 0;
    iterations = 0;
    if( (chi_sq = get_chisq(N, D, f, a_data, x_data, y_data, sigma_data)) < 0)
    {
        cpl_free(ia_local);
        cpl_matrix_delete(alpha);
        cpl_matrix_delete(beta);
        cpl_free(a_da);
        cpl_free(part);
        return cpl_error_set_where_();
    }

    /* Iterate until chi^2 didn't improve significantly many
       times in a row (where 'many' is defined by tolerance_count) */
    while (count < tolerance_count &&
       lambda < MAXLAMBDA &&
       iterations < max_iterations)
    {
        /* In each iteration lambda increases, or chi^2 decreases or
           count increases. Because chi^2 is bounded from below
           (and lambda and count from above), the loop will terminate */

        double chi_sq_candidate = 0.0;
        double chi_sq_diff;
        int canderror;
        int cont = 3; /* Try to get a candidate this many times */



        /* Get candidate position in parameter space = a+da,
         * where  alpha * da = beta .
         * Increase lambda until alpha is non-singular
         */

        do {
            cpl_errorstate prevstate = cpl_errorstate_get();
            canderror = get_candidate(a_data, ia_local,
                                      M, N, D,
                                      lambda, f, dfda,
                                      x_data, y_data, sigma_data,
                                      part, alpha, beta, a_da);
            cont--;

            if (canderror || !cpl_errorstate_is_equal(prevstate)) {
                if (cont > 0 && !cpl_errorstate_is_equal(prevstate) &&
                    cpl_error_get_code() == CPL_ERROR_SINGULAR_MATRIX &&
                    lambda < MAXLAMBDA) {
                    /* Recover since lambda did not diverge */
                    lambda *= 9.0;
                    cpl_errorstate_set(prevstate);
                    canderror = 1; /* Ensure another try */
                } else {
                    cpl_free(ia_local);
                    cpl_matrix_delete(alpha);
                    cpl_matrix_delete(beta);
                    cpl_free(a_da);
                    cpl_free(part);

                    /* Set error either if lambda diverged or if 
                       get_candidate() failed in some other way */
                    return cpl_error_set_message_(cpl_error_get_code() ==
                                                  CPL_ERROR_SINGULAR_MATRIX
                                                  ? CPL_ERROR_CONTINUE
                                                  : cpl_error_get_code(),
                                                  "get_candidate()=%d. cont=%d. "
                                                  "lambda = %g", canderror, cont,
                                                  lambda);
                }
            }
        } while (canderror);

        /* Get chi^2(a+da) */
        chi_sq_candidate = get_chisq(N, D, f, a_da, x_data, y_data, sigma_data);
        /* Check for invalid candidate, including NaN */
        if (chi_sq_candidate < 0.0 || chi_sq_candidate != chi_sq_candidate)
        {
            cpl_free(ia_local);
            cpl_matrix_delete(alpha);
            cpl_matrix_delete(beta);
            cpl_free(a_da);
            cpl_free(part);
            return cpl_error_get_code()
                ? cpl_error_set_where_()
                : cpl_error_set_(CPL_ERROR_CONTINUE);
        }

	chi_sq_diff = chi_sq_candidate - chi_sq;
        if (chi_sq_diff > sqrt(DBL_EPSILON) )
        {
            /* Move towards steepest descent */
            lambda *= 9.0;
        }
        else
        {
            /* Move towards Newton's algorithm */
            lambda /= 10.0;

            /* Count the number of successive improvements in chi^2 of
               less than 0.01 (default) relative */
            if ( chi_sq < sqrt(DBL_EPSILON) ||
             (chi_sq - chi_sq_candidate)/chi_sq < relative_tolerance)
            {
                count += 1;
            }
            else
            {
                /* Chi^2 improved by a significant amount,
                   reset counter */
                count = 0;
            }

            /* chi^2 improved, update a and chi^2 */
            {
            cpl_size i;
            for (i = 0; i < M; i++) a_data[i] = a_da[i];
            }
            chi_sq = chi_sq_candidate;
        }
        iterations++;

    }

    /* Set error if we didn't converge */
    if ( !( lambda < MAXLAMBDA && iterations < max_iterations ) )
    {
        cpl_free(ia_local);
        cpl_matrix_delete(alpha);
        cpl_matrix_delete(beta);
        cpl_free(a_da);
        cpl_free(part);
        return cpl_error_set_(CPL_ERROR_CONTINUE);
    }

    /* Compute mse if requested */
    if (mse != NULL)
    {
        cpl_size i;

        *mse = 0.0;

        for(i = 0; i < N; i++)
        {
            double fx_i = 0.0;
            double residual = 0.0;

            /* Evaluate f(x_i) at the best fit parameters */
            if( f(&(x_data[i*D]),
              a_data,
              &fx_i) != 0)
            {
                cpl_free(ia_local);
                cpl_matrix_delete(alpha);
                cpl_matrix_delete(beta);
                cpl_free(a_da);
                cpl_free(part);
                return cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
            }

            residual = y_data[i] - fx_i;
            *mse += residual * residual;
        }
        *mse /= (double)N;
    }

    /* Compute reduced chi^2 if requested */
    if (red_chisq != NULL)
    {
        /* We already know the optimal chi^2 (and that N > Mfit)*/
        *red_chisq = chi_sq / (double)(N-Mfit);
    }

    /* Compute covariance matrix if requested
     * cov = alpha(lambda=0)^-1
     */
    if (covariance != NULL)
    {
        cpl_matrix *cov;

        if( get_candidate(a_data, ia_local,
                  M, N, D, 0.0, f, dfda,
                  x_data, y_data, sigma_data,
                  part, alpha, beta, a_da)
        != 0)
        {
            cpl_free(ia_local);
            cpl_matrix_delete(alpha);
            cpl_matrix_delete(beta);
            cpl_free(a_da);
            cpl_free(part);
            return cpl_error_set_where_();
        }

        cov = cpl_matrix_invert_create(alpha);
        if (cov == NULL)
        {
            cpl_free(ia_local);
            cpl_matrix_delete(alpha);
            cpl_matrix_delete(beta);
            cpl_free(a_da);
            cpl_free(part);
            return cpl_error_set_where_();
        }

        /* Make sure that variances are positive */
        {
        cpl_size i;
        for (i = 0; i < Mfit; i++)
            {
            if ( !(cpl_matrix_get(cov, i, i) > 0) )
                {
                cpl_free(ia_local);
                cpl_matrix_delete(alpha);
                cpl_matrix_delete(beta);
                cpl_free(a_da);
                cpl_free(part);
                cpl_matrix_delete(cov);
                *covariance = NULL;
                return cpl_error_set_(CPL_ERROR_SINGULAR_MATRIX);
                }
            }
        }

        /* Expand covariance matrix from Mfit x Mfit
           to M x M. Set rows/columns corresponding to fixed
           parameters to zero */

        *covariance = cpl_matrix_new(M, M);
        if (*covariance == NULL)
        {
            cpl_free(ia_local);
            cpl_matrix_delete(alpha);
            cpl_matrix_delete(beta);
            cpl_free(a_da);
            cpl_free(part);
            cpl_matrix_delete(cov);
            return cpl_error_set_(CPL_ERROR_ILLEGAL_OUTPUT);
        }

        {
            cpl_size j, jmfit;

        for (j = 0, jmfit = 0; j < M; j++)
            if (ia_local[j] != 0)
            {
                cpl_size i, imfit;

                for (i = 0, imfit = 0; i < M; i++)
                if (ia_local[i] != 0)
                    {
                    cpl_matrix_set(*covariance, i, j,
                               cpl_matrix_get(
                               cov, imfit, jmfit));
                    imfit += 1;
                    }

                assert( imfit == Mfit );

                jmfit += 1;
            }

        assert( jmfit == Mfit );
        }

        cpl_matrix_delete(cov);
    }


    cpl_free(ia_local);
    cpl_matrix_delete(alpha);
    cpl_matrix_delete(beta);
    cpl_free(a_da);
    cpl_free(part);

    return CPL_ERROR_NONE;
}
#endif /* CPL_VECTOR_FIT_IMPL_H */
