/*
 * This file is part of the HDRL
 * Copyright (C) 2013,2014 European Southern Observatory
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

#if HDRL_OPERATION == HDRL_IMLIST_BASIC_IMLIST

    cpl_size i;

    /* Check input */
    cpl_ensure_code(himlist1, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(himlist2, CPL_ERROR_NULL_INPUT);

    /* Check image sets compatibility     */
    cpl_ensure_code( himlist1->ni == himlist2->ni, CPL_ERROR_ILLEGAL_INPUT);

    /* Loop on the planes and apply the operation */
    for (i=0; i<himlist1->ni; i++) {
        const cpl_error_code error_code = 
            HDRL_OPERATOR(himlist1->images[i], himlist2->images[i]);
        cpl_ensure_code(!error_code, error_code);
    }

    return CPL_ERROR_NONE;

#elif HDRL_OPERATION == HDRL_IMLIST_BASIC_IMAGE

    cpl_size i;

    /* Check input */
    cpl_ensure_code(himlist, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(himg, CPL_ERROR_NULL_INPUT);

    /* Loop on the planes and apply the operation */
    for (i=0; i<himlist->ni; i++) {
        const cpl_error_code error_code = HDRL_OPERATOR(himlist->images[i], 
                himg);
        cpl_ensure_code(!error_code, error_code);
    }

    return CPL_ERROR_NONE;

#elif HDRL_OPERATION == HDRL_IMLIST_BASIC_SCALAR

    cpl_size i;

    /* Check input */
    cpl_ensure_code(himlist, CPL_ERROR_NULL_INPUT);

    /* Loop on the planes and apply the operation */
    for (i=0; i<himlist->ni; i++) {
        const cpl_error_code error_code = HDRL_OPERATOR(himlist->images[i],
                                                        value);
        cpl_ensure_code(!error_code, error_code);
    }

    return CPL_ERROR_NONE;


#endif

