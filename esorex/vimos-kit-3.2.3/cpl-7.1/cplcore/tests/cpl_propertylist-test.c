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

#include <complex.h>

#include "cpl_init.h"
#include "cpl_image_io.h"
#include "cpl_error.h"
#include "cpl_propertylist_impl.h"
#include "cpl_table.h"
#include "cpl_msg.h"
#include "cpl_test.h"
#include "cpl_memory.h"

#include "cpl_fits.h"
#include "cpl_io_fits.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <math.h>

#include <fitsio.h>

/*----------------------------------------------------------------------------
                 Defines
 ----------------------------------------------------------------------------*/

#define BASE "cpl_propertylist-test"

#define NOKEY "Non-existing key"

#define BADNUM "l.O"
/* This complex format is not supported */
#define COMPLEXVAL "1D1,2E2"

#define LONGNAME80 "0123456789012345678901234567890123456789"   \
    "0123456789012345678901234567890123456789"

/*----------------------------------------------------------------------------
                 Functions prototypes
 ----------------------------------------------------------------------------*/

static int cpl_test_property_compare_name(const cpl_property  *,
                                          const cpl_property  *);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{

    const char *keys[] = {
        "a", "b", "c", "d", "e", "f", "g", "h", "i",
        "A", "B", "C", "D", "E", "F", "G", "H", "I"
    };

    const char *comments[] = {
        "A character value",
        "A boolean value",
        "A integer value",
        "A long integer value",
        "A floating point number",
        "A double precision number",
        "A string value",
        "A floating point complex number",
        "A double precision complex number"
    };

    cpl_type types[] = {
        CPL_TYPE_CHAR,
        CPL_TYPE_BOOL,
        CPL_TYPE_INT,
        CPL_TYPE_LONG,
        CPL_TYPE_FLOAT,
        CPL_TYPE_DOUBLE,
        CPL_TYPE_STRING,
        CPL_TYPE_FLOAT_COMPLEX,
        CPL_TYPE_DOUBLE_COMPLEX
    };

    int status, j;

    long i;
    long naxes[2] = {256,256};

    const float fval0 = -1.23456789;
    float fval1, fval2;
    const double dval0 = -1.23456789;
    double dval1, dval2;
    const float complex zf0 = fval0 + fval0 * fval0 * _Complex_I;
    float  complex zf1, zf2;
    const double complex zd0 = zf0;
    double complex zd1, zd2;
    float f0,f1,f2;
    double d0,d1,d2;
    const int nprops = sizeof(types)/sizeof(types[0]);

    const char *filename3 = BASE ".fits";
    /*const char regex[] = "^(DATE-OBS|ESO|COMMENT|HISTORY)$";*/
    const char to_rm[] = "^HIERARCH ESO |^NBXIS1$|^NBXIS2$|^HISTORY$";
    cpl_propertylist *plist, *_plist, *plist2;
    cpl_property *pro, *pro2;
    const cpl_property *pro1;
    cpl_table *t;
    cpl_type t1;
    cpl_error_code code;
    int value;
    cpl_boolean   do_bench;
    int bench_size;
    const char * strval;

    fitsfile *fptr = NULL;
    fitsfile *_fptr;

    FILE *stream;

#ifdef CPL_PROPERTYLIST_TEST_FCARD

    struct fcard {
        const char *key;
        const char *val;
        const char *com;
        cpl_type type;
    };


    struct fcard hdr[] = {
        {"SIMPLE", "T",
         "Standard FITS format (NOST-100.0)",
         CPL_TYPE_BOOL},
        {"BITPIX", "16",
         "# of bits storing pix values",
         CPL_TYPE_INT},
        {"NAXIS", "2",
         "# of axes in frame",
         CPL_TYPE_INT},
        {"NAXIS1", "2148",
         "# pixels/axis",
         CPL_TYPE_INT},
        {"NAXIS2", "2340",
         "# pixels/axis",
         CPL_TYPE_INT},
        {"ORIGIN", "ESO",
         "European Southern Observatory",
         CPL_TYPE_STRING},
        {"DATE", "2002-03-08T04:27:21.420",
         "Date this file was written (dd/mm/yyyy)",
         CPL_TYPE_STRING},
        {"MJD-OBS", "52341.17813019",
         "Obs start 2002-03-08T04:16:30.448",
         CPL_TYPE_DOUBLE},
        {"DATE-OBS", "2002-03-08T04:16:30.448",
         "Date of observation",
         CPL_TYPE_STRING},
        {"EXPTIME", "600.000",
         "Total integration time. 00:10:00.000",
         CPL_TYPE_DOUBLE},
        {"TELESCOP", "VLT",
         "ESO <TEL>",
         CPL_TYPE_STRING},
        {"RA", "181.41734",
         "12:05:40.1 RA (J2000) pointing",
         CPL_TYPE_DOUBLE},
        {"DEC", "-7.65555",
         "-07:39:19.9  DEC (J2000) pointing",
         CPL_TYPE_DOUBLE},
        {"EQUINOX", "2000.",
         "Standard FK5 (years)",
         CPL_TYPE_DOUBLE},
        {"RADECSYS", "FK5",
         "Coordinate reference frame",
         CPL_TYPE_STRING},
        {"LST", "38309.370",
         "10:38:29.370 LST at start",
         CPL_TYPE_DOUBLE},
        {"UTC", "15438.000",
         "04:17:18.000 UT at start",
         CPL_TYPE_DOUBLE},
        {"OBSERVER", "UNKNOWN",
         "Name of observer",
         CPL_TYPE_STRING},
        {"INSTRUME", "UNKNOWN",
         "Instrument used",
         CPL_TYPE_STRING},
        {"PI-COI", "'555555555'",
         "Name of PI and COI",
         CPL_TYPE_STRING},
        {"OBJECT", "None",
         "Original target",
         CPL_TYPE_STRING},
        {"PCOUNT", "0",
         "Number of parameters per group",
         CPL_TYPE_INT},
        {"GCOUNT", "1",
         "Number of groups",
         CPL_TYPE_INT},
        {"CRVAL1", "181.41734",
         "12:05:40.1, RA at ref pixel",
         CPL_TYPE_DOUBLE},
        {"CRPIX1", "2341.8585366",
         "Reference pixel in X",
         CPL_TYPE_DOUBLE},
        {"CDELT1", "0.20500000",
         "SS arcsec per pixel in RA",
         CPL_TYPE_DOUBLE},
        {"CTYPE1", "RA---TAN",
         "pixel coordinate system",
         CPL_TYPE_STRING},
        {"CRVAL2", "-7.65555",
         "-07:39:19.9, DEC at ref pixel",
         CPL_TYPE_DOUBLE},
        {"CRPIX2", "2487.8585366",
         "Reference pixel in Y",
         CPL_TYPE_DOUBLE},
        {"CDELT2", "0.20500000",
         "SS arcsec per pixel in DEC",
         CPL_TYPE_DOUBLE},
        {"CTYPE2", "DEC--TAN",
         "pixel coordinate system",
         CPL_TYPE_STRING},
        {"BSCALE", "1.0",
         "pixel=FITS*BSCALE+BZERO",
         CPL_TYPE_DOUBLE},
        {"BZERO", "32768.0",
         "pixel=FITS*BSCALE+BZERO",
         CPL_TYPE_DOUBLE},
        {"CD1_1", "0.000057",
         "Translation matrix element",
         CPL_TYPE_DOUBLE},
        {"CD1_2", "0.000000",
         "Translation matrix element",
         CPL_TYPE_DOUBLE},
        {"CD2_1", "0.000000",
         "Translation matrix element",
         CPL_TYPE_DOUBLE},
        {"CD2_2", "0.000057",
         "Translation matrix element",
         CPL_TYPE_DOUBLE},
        {"HIERARCH ESO OBS DID", "ESO-VLT-DIC.OBS-1.7",
         "OBS Dictionary",
         CPL_TYPE_STRING},
        {"HIERARCH ESO OBS OBSERVER", "UNKNOWN",
         "Observer Name",
         CPL_TYPE_STRING},
        {"HIERARCH ESO OBS PI-COI NAME", "UNKNOWN",
         "PI-COI name",
         CPL_TYPE_STRING},
        {"HIERARCH ESO INS GRAT NAME", "HR",
         "Grating name",
         CPL_TYPE_STRING},
        {"HIERARCH ESO PRO CATG", "X",
         "Product category",
         CPL_TYPE_STRING},
        {"HIERARCH ESO TPL NEXP", "5",
         "Number of exposures",
         CPL_TYPE_INT},
        {"HISTORY", "1st history record", NULL, CPL_TYPE_STRING},
        {"COMMENT", "1st comment record", NULL, CPL_TYPE_STRING},
        {"HISTORY", "2st history record", NULL, CPL_TYPE_STRING},
        {"COMMENT", "2st comment record", NULL, CPL_TYPE_STRING},
        {"COMMENT", "3st comment record", NULL, CPL_TYPE_STRING},
        {"HISTORY", "3st history record", NULL, CPL_TYPE_STRING},
        {"END", NULL, NULL, CPL_TYPE_STRING}
    };
#endif

    const char *longname = LONGNAME80;

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    do_bench = cpl_msg_get_level() <= CPL_MSG_INFO ? CPL_TRUE : CPL_FALSE;

    /* Always test the _dump functions, but produce no output
       when the message level is above info */
    stream = cpl_msg_get_level() > CPL_MSG_INFO
        ? fopen("/dev/null", "a") : stdout;

    cpl_test_nonnull(stream);

    cpl_test_eq(sizeof(comments)/sizeof(comments[0]), nprops);
    cpl_test_eq(sizeof(keys)/sizeof(keys[0]), 2 * nprops);

    /*
     * Test 1: Create a property list and check its validity.
     */

    plist = cpl_propertylist_new();

    cpl_test_nonnull(plist);
    cpl_test(cpl_propertylist_is_empty(plist));
    cpl_test_zero(cpl_propertylist_get_size(plist));

    pro1 = cpl_propertylist_get_const(plist, 100);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_null(pro1);

    /*
     * Test 2: Append properties to the list created in the previous test
     *         and verify the data.
     */

    cpl_propertylist_append_char(plist, keys[0], 'a');
    cpl_propertylist_set_comment(plist, keys[0], comments[0]);

    cpl_propertylist_append_bool(plist, keys[1], 1);
    cpl_propertylist_set_comment(plist, keys[1], comments[1]);

    cpl_propertylist_append_int(plist, keys[2], -1);
    cpl_propertylist_set_comment(plist, keys[2], comments[2]);

    cpl_propertylist_append_long(plist, keys[3], 32768);
    cpl_propertylist_set_comment(plist, keys[3], comments[3]);

    cpl_propertylist_append_float(plist, keys[4], fval0);
    cpl_propertylist_set_comment(plist, keys[4], comments[4]);

    cpl_propertylist_append_double(plist, keys[5], dval0);
    cpl_propertylist_set_comment(plist, keys[5], comments[5]);

    cpl_propertylist_append_string(plist, keys[6], comments[6]);
    cpl_propertylist_set_comment(plist, keys[6], comments[6]);

    cpl_propertylist_append_float_complex(plist, keys[7], zf0);
    cpl_propertylist_set_comment(plist, keys[7], comments[7]);

    cpl_propertylist_append_double_complex(plist, keys[8], zd0);
    cpl_propertylist_set_comment(plist, keys[8], comments[8]);

    cpl_test_zero(cpl_propertylist_is_empty(plist));

    cpl_test_eq(cpl_propertylist_get_size(plist), nprops);

    pro1 = cpl_propertylist_get_const(plist, nprops-1);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(pro1);

    pro1 = cpl_propertylist_get_const(plist, nprops);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_null(pro1);

    code = cpl_propertylist_save(plist, BASE "_2.fits", CPL_IO_CREATE);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_2.fits");

    _plist = cpl_propertylist_load(BASE "_2.fits", 0);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(_plist);
    /* Some properties are added */
    cpl_test_leq(nprops, cpl_propertylist_get_size(_plist));

    for (i = 0; i < nprops; i++) {
        const cpl_property * p1 =
            cpl_propertylist_get_property_const(plist, keys[i]);
        const cpl_property * p2 =
            cpl_propertylist_get_property_const(_plist, keys[i + nprops]);

        cpl_test_eq(i, i);

        cpl_test_error(CPL_ERROR_NONE);
        cpl_test_nonnull(p1);
        cpl_test_nonnull(p2);

        if (types[i] == CPL_TYPE_CHAR) {
            /* FITS I/O promotes this type to a single character string */
            cpl_test_eq(cpl_property_get_type(p2), CPL_TYPE_STRING);
            cpl_test_eq(cpl_property_get_string(p2)[0], cpl_property_get_char(p1));
        } else if (types[i] == CPL_TYPE_LONG) {
            /* FITS I/O casts this type to int */
            cpl_test_eq(cpl_property_get_type(p2), CPL_TYPE_INT);
            cpl_test_eq(cpl_property_get_int(p2), cpl_property_get_long(p1));
        } else if (types[i] == CPL_TYPE_FLOAT) {
            /* FITS I/O promotes this type to double */
            cpl_test_eq(cpl_property_get_type(p2), CPL_TYPE_DOUBLE);
            cpl_test_abs(cpl_property_get_double(p2),
                         cpl_property_get_float(p1), 2.0 * FLT_EPSILON);
        } else if (types[i] == CPL_TYPE_FLOAT_COMPLEX) {
            /* FITS I/O promotes this type to double complex */
            cpl_test_eq(cpl_property_get_type(p2), CPL_TYPE_DOUBLE_COMPLEX);
            zd1 = cpl_property_get_float_complex(p1);
            zd2 = cpl_property_get_double_complex(p2);
            cpl_test_abs_complex(zd1, zd2, 2.0 * FLT_EPSILON);
        } else {
            cpl_test_eq(cpl_property_get_type(p2), types[i]);
            if (types[i] == CPL_TYPE_BOOL) {
                cpl_test_eq(cpl_property_get_bool(p2),
                            cpl_property_get_bool(p1));
            } else if (types[i] == CPL_TYPE_INT) {
                cpl_test_eq(cpl_property_get_int(p2), cpl_property_get_int(p1));
            } else if (types[i] == CPL_TYPE_DOUBLE) {
                cpl_test_abs(cpl_property_get_double(p2),
                             cpl_property_get_double(p1), DBL_EPSILON);
            } else if (types[i] == CPL_TYPE_STRING) {
                cpl_test_eq_string(cpl_property_get_string(p2),
                                   cpl_property_get_string(p1));
            } else if (types[i] == CPL_TYPE_DOUBLE_COMPLEX) {
                zd1 = cpl_property_get_double_complex(p1);
                zd2 = cpl_property_get_double_complex(p2);
                cpl_test_abs_complex(zd1, zd2, 32.0 * DBL_EPSILON);
            }
        }

        cpl_test_zero(strncmp(cpl_property_get_comment(p2), comments[i],
                              strlen(cpl_property_get_comment(p2))));

    }

    cpl_propertylist_delete(_plist);

    for (i = 0; i < nprops; i++) {
        const cpl_property *p = cpl_propertylist_get_const(plist, i);

        cpl_test_eq_string(cpl_property_get_name(p), keys[i]);
        cpl_test_eq_string(cpl_property_get_comment(p), comments[i]);
        cpl_test_eq(cpl_property_get_type(p), types[i]);

        cpl_test(cpl_propertylist_has(plist, keys[i]));
        cpl_test_eq_string(cpl_propertylist_get_comment(plist, keys[i]),
                           comments[i]);
        cpl_test_eq(cpl_propertylist_get_type(plist, keys[i]), types[i]);

    }

    cpl_test_eq(cpl_propertylist_get_char(plist, keys[0]), 'a');
    cpl_test_eq(cpl_propertylist_get_bool(plist, keys[1]), 1);
    cpl_test_eq(cpl_propertylist_get_int(plist, keys[2]), -1);
    cpl_test_eq(cpl_propertylist_get_long(plist, keys[3]), 32768);

    fval1 = cpl_propertylist_get_float(plist, keys[4]);
    cpl_test_abs(fval0, fval1, 0.0);

    dval1 = cpl_propertylist_get_double(plist, keys[5]);
    cpl_test_abs(dval0, dval1, 0.0);

    cpl_test_eq_string(cpl_propertylist_get_string(plist, keys[6]),
                       comments[6]);

    /*
     * Test 3: Modify the values of the property list entries
     *         and verify the data.
     */

    cpl_test_zero(cpl_propertylist_set_char(plist, keys[0], 'b'));
    cpl_test_eq(cpl_propertylist_get_char(plist, keys[0]), 'b');

    cpl_test_zero(cpl_propertylist_set_bool(plist, keys[1], 0));
    cpl_test_zero(cpl_propertylist_get_bool(plist, keys[1]));

    cpl_test_zero(cpl_propertylist_set_int(plist, keys[2], -1));
    cpl_test_eq(cpl_propertylist_get_int(plist, keys[2]), -1);

    cpl_test_zero(cpl_propertylist_set_long(plist, keys[3], 1));
    cpl_test_eq(cpl_propertylist_get_long(plist, keys[3]), 1);

    fval1 = 9.87654321;
    cpl_test_noneq(fval0, fval1);
    cpl_test_zero(cpl_propertylist_set_float(plist, keys[4], fval1));
    fval2 = cpl_propertylist_get_float(plist, keys[4]);
    cpl_test_abs(fval1, fval2, 0.0);

    dval1 = -9.87654321;
    cpl_test_noneq(dval0, dval1);
    cpl_test_zero(cpl_propertylist_set_double(plist, keys[5], dval1));
    dval2 = cpl_propertylist_get_double(plist, keys[5]);
    cpl_test_abs(dval1, dval2, 0.0);

    cpl_test_zero(cpl_propertylist_set_string(plist, keys[6], comments[0]));
    cpl_test_eq_string(cpl_propertylist_get_string(plist, keys[6]),
                       comments[0]);

    zf1 = 9.87654321 * zf0;
    cpl_test_noneq(zf0, zf1);
    cpl_test_zero(cpl_propertylist_set_float_complex(plist, keys[7], zf1));
    zf2 = cpl_propertylist_get_float_complex(plist, keys[7]);
    cpl_test_abs_complex(zf1, zf2, 0.0);

    zd1 = -9.87654321 * zd0;
    cpl_test_noneq(zd0, zd1);
    cpl_test_zero(cpl_propertylist_set_double_complex(plist, keys[8], zd1));
    zd2 = cpl_propertylist_get_double_complex(plist, keys[8]);
    cpl_test_abs_complex(zd1, zd2, 0.0);


    /*
     * Test 4: Check that trying to modify an entry with a different
     *         type is properly failing.
     */

    code = cpl_propertylist_set_char(plist, keys[1], 'a');
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);

    code = cpl_propertylist_set_bool(plist, keys[2], 1);
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);

    code = cpl_propertylist_set_int(plist, keys[3], 1);
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);

    code = cpl_propertylist_set_long(plist, keys[4], 1);
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);

    code = cpl_propertylist_set_float(plist, keys[5], 1.);
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);

    code = cpl_propertylist_set_double(plist, keys[6], 1.);
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);

    code = cpl_propertylist_set_string(plist, keys[0], comments[0]);
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);

    code = cpl_propertylist_set_float(plist, keys[7], 1.);
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);

    code = cpl_propertylist_set_double(plist, keys[8], 1.);
    cpl_test_eq_error(code, CPL_ERROR_TYPE_MISMATCH);


    /*
     * Test 5: Verify that values are inserted correctly into the property
     *         list.
     */

    cpl_test_eq(cpl_propertylist_insert_char(plist, keys[0],
                                           keys[0 + nprops], 'a'), 0);

    cpl_test_eq(cpl_propertylist_insert_after_char(plist, keys[0],
                                                 keys[0 + nprops], 'c'), 0);

    cpl_test_eq(cpl_propertylist_insert_bool(plist, keys[1],
                                           keys[1 + nprops], 0), 0);

    cpl_test_eq(cpl_propertylist_insert_after_bool(plist, keys[1],
                                                 keys[1 + nprops], 1), 0);


    cpl_test_eq(cpl_propertylist_insert_int(plist, keys[2],
                                          keys[2 + nprops], 0), 0);

    cpl_test_eq(cpl_propertylist_insert_after_int(plist, keys[2],
                                                keys[2 + nprops], 1), 0);


    cpl_test_eq(cpl_propertylist_insert_long(plist, keys[3], keys[3 + nprops],
                                           123456789), 0);

    cpl_test_eq(cpl_propertylist_insert_after_long(plist, keys[3],
                                                   keys[3 + nprops],
                                                   123456789), 0);


    cpl_test_eq(cpl_propertylist_insert_float(plist, keys[4], keys[4 + nprops],
                                            fval0), 0);

    cpl_test_eq(cpl_propertylist_insert_after_float(plist, keys[4],
                                                    keys[4 + nprops],
                                                  -fval0), 0);


    cpl_test_eq(cpl_propertylist_insert_double(plist, keys[5], keys[5 + nprops],
                                             dval0), 0);

    cpl_test_eq(cpl_propertylist_insert_after_double(plist, keys[5],
                                                     keys[5 + nprops],
                                                   -dval0), 0);


    cpl_test_eq(cpl_propertylist_insert_string(plist, keys[6],
                                             keys[6 + nprops], ""), 0);

    cpl_test_eq(cpl_propertylist_insert_after_string(plist, keys[6],
                                                   keys[6 + nprops], ""), 0);

    cpl_test_eq(cpl_propertylist_insert_float(plist, keys[7], keys[7 + nprops],
                                            fval0), 0);

    cpl_test_eq(cpl_propertylist_insert_after_float(plist, keys[7],
                                                    keys[7 + nprops],
                                                  -fval0), 0);


    cpl_test_eq(cpl_propertylist_insert_double(plist, keys[8], keys[8 + nprops],
                                             dval0), 0);

    cpl_test_eq(cpl_propertylist_insert_after_double(plist, keys[8],
                                                     keys[8 + nprops],
                                                   -dval0), 0);


    for (i = 0; i < nprops; i++) {
        cpl_property *p0 = cpl_propertylist_get(plist, 3 * i);
        cpl_property *p1 = cpl_propertylist_get(plist, 3 * i + 1);
        cpl_property *p2 = cpl_propertylist_get(plist, 3 * i + 2);


        cpl_test_eq_string(cpl_property_get_name(p0), keys[i + nprops]);
        cpl_test_eq_string(cpl_property_get_name(p1), keys[i]);
        cpl_test_eq_string(cpl_property_get_name(p2), keys[i + nprops]);

        switch (cpl_property_get_type(p0)) {
        case CPL_TYPE_CHAR:
            cpl_test_eq(cpl_property_get_char(p0), 'a');
            cpl_test_eq(cpl_property_get_char(p2), 'c');
            break;

        case CPL_TYPE_BOOL:
            cpl_test_zero(cpl_property_get_bool(p0));
            cpl_test_eq(cpl_property_get_bool(p2), 1);
            break;

        case CPL_TYPE_INT:
            cpl_test_zero(cpl_property_get_int(p0));
            cpl_test_zero(cpl_property_get_int(p0));
            break;

        case CPL_TYPE_LONG:
            cpl_test_eq(cpl_property_get_long(p0), 123456789);
            cpl_test_eq(cpl_property_get_long(p0), 123456789);
            break;

        case CPL_TYPE_FLOAT:
            fval1 = cpl_property_get_float(p0);
            cpl_test_abs(fval0, fval1, 0.0);

            fval1 = -cpl_property_get_float(p2);
            cpl_test_abs(fval0, fval1, 0.0);
            break;

        case CPL_TYPE_DOUBLE:
            dval1 = cpl_property_get_double(p0);
            cpl_test_abs(dval0, dval1, 0.0);

            dval1 = -cpl_property_get_double(p2);
            cpl_test_abs(dval0, dval1, 0.0);
            break;

        case CPL_TYPE_STRING:
            cpl_test_eq_string(cpl_property_get_string(p0), "");
            cpl_test_eq_string(cpl_property_get_string(p2), "");
            break;

        default:
            /* This point should never be reached */
            cpl_test(0);
            break;
        }
    }


    /*
     * Test 6: Verify that modification of or insertion at/after a non
     *         existing elements is reported correctly.
     */

    cpl_test_zero(cpl_propertylist_has(plist, NOKEY));

    code = cpl_propertylist_set_char(plist, NOKEY, 'a');
    cpl_test_eq_error(code, CPL_ERROR_DATA_NOT_FOUND);

    code = cpl_propertylist_set_bool(plist, NOKEY, 1);
    cpl_test_eq_error(code, CPL_ERROR_DATA_NOT_FOUND);

    code = cpl_propertylist_set_int(plist, NOKEY, 1);
    cpl_test_eq_error(code, CPL_ERROR_DATA_NOT_FOUND);

    code = cpl_propertylist_set_long(plist, NOKEY, 1);
    cpl_test_eq_error(code, CPL_ERROR_DATA_NOT_FOUND);

    code = cpl_propertylist_set_float(plist, NOKEY, 1);
    cpl_test_eq_error(code, CPL_ERROR_DATA_NOT_FOUND);

    code = cpl_propertylist_set_double(plist, NOKEY, 1);
    cpl_test_eq_error(code, CPL_ERROR_DATA_NOT_FOUND);

    code = cpl_propertylist_set_string(plist, NOKEY, "");
    cpl_test_eq_error(code, CPL_ERROR_DATA_NOT_FOUND);

    code = cpl_propertylist_insert_char(plist, NOKEY, "h", 'a');
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_bool(plist, NOKEY, "h", 1);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_int(plist, NOKEY, "h", 1);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_long(plist, NOKEY, "h", 1);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_float(plist, NOKEY, "h", 1.0);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_double(plist, NOKEY, "h", 1.0);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_string(plist, NOKEY, "h", "");
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_after_char(plist, NOKEY, "h", 'a');
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_after_bool(plist, NOKEY, "h", 1);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_after_int(plist, NOKEY, "h", 1);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_after_long(plist, NOKEY, "h", 1);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_after_float(plist, NOKEY, "h", 1.0);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_after_double(plist, NOKEY, "h", 1.0);
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    code = cpl_propertylist_insert_after_string(plist, NOKEY, "h", "");
    cpl_test_eq_error(code, CPL_ERROR_UNSPECIFIED);

    /*
     * Test 7: Create a copy of the property list and verify that original
     *         and copy are identical but do not share any resources.
     */

    _plist = cpl_propertylist_duplicate(plist);
    plist2 = cpl_propertylist_duplicate(plist);
    cpl_test_nonnull(_plist);
    cpl_test_assert(_plist != plist);
    cpl_test_nonnull(plist2);
    cpl_test_assert(plist2 != plist);

    cpl_test_error(CPL_ERROR_NONE);

    for (i = 0; i < cpl_propertylist_get_size(plist); i++) {
        cpl_property *p = cpl_propertylist_get(plist, i);
        cpl_property *_p = cpl_propertylist_get(_plist, i);

        cpl_test_assert(cpl_property_get_name(p) != cpl_property_get_name(_p));
        cpl_test_eq_string(cpl_property_get_name(p),
                           cpl_property_get_name(_p));
        cpl_test_assert(cpl_property_get_comment(p) == NULL ||
                  (cpl_property_get_comment(p) != cpl_property_get_comment(_p)));
        cpl_test_assert(cpl_property_get_comment(p) == NULL ||
                  !strcmp(cpl_property_get_comment(p), cpl_property_get_comment(_p)));

        switch (cpl_property_get_type(p)) {
        case CPL_TYPE_CHAR:
            cpl_test_eq(cpl_property_get_char(p),
                      cpl_property_get_char(_p));
            break;

        case CPL_TYPE_BOOL:
            cpl_test_eq(cpl_property_get_bool(p),
                      cpl_property_get_bool(_p));
            break;

        case CPL_TYPE_INT:
            cpl_test_eq(cpl_property_get_int(p),
                      cpl_property_get_int(_p));
            break;

        case CPL_TYPE_LONG:
            cpl_test_eq(cpl_property_get_long(p),
                      cpl_property_get_long(_p));
            break;

        case CPL_TYPE_FLOAT:
            fval1 = cpl_property_get_float(p);
            fval2 = cpl_property_get_float(_p);
            cpl_test_abs(fval1, fval2, 0.0);
            break;

        case CPL_TYPE_DOUBLE:
            dval1 = cpl_property_get_double(p);
            dval2 = cpl_property_get_double(_p);
            cpl_test_abs(dval1, dval2, 0.0);
            break;

        case CPL_TYPE_STRING:
            cpl_test_eq_string(cpl_property_get_string(p),
                               cpl_property_get_string(_p));
            break;

        case CPL_TYPE_FLOAT_COMPLEX:
            zf1 = cpl_property_get_float_complex(p);
            zf2 = cpl_property_get_float_complex(_p);
            cpl_test_abs_complex(zf1, zf2, 0.0);
            break;

        case CPL_TYPE_DOUBLE_COMPLEX:
            zd1 = cpl_property_get_double_complex(p);
            zd2 = cpl_property_get_double_complex(_p);
            cpl_test_abs_complex(zd1, zd2, 0.0);
            break;

        default:
            /* This point should never be reached */
            cpl_test(0);
            break;
        }
    }

    cpl_propertylist_delete(_plist);

    cpl_test_error(CPL_ERROR_NONE);

    /*
     * Test 8: Erase elements from the property list and verify the list
     *         structure and the data. Each key exists twice in capital
     *         letters and once in lowercase letters.
     */

    for (i = 0; i < nprops; i++) {
        cpl_propertylist_erase(plist, keys[i + nprops]);
        cpl_test_eq(cpl_propertylist_has(plist, keys[i + nprops]), 1);

        cpl_propertylist_erase(plist, keys[i + nprops]);
        cpl_test_zero(cpl_propertylist_has(plist, keys[i + nprops]));
    }
    cpl_test_eq(cpl_propertylist_get_size(plist), nprops);

    for (i = 0; i < nprops; i++) {
        cpl_property *p = cpl_propertylist_get(plist, i);
        cpl_test_eq_string(cpl_property_get_name(p), keys[i]);
    }

    cpl_test_eq(cpl_propertylist_get_char(plist, keys[0]), 'b');
    cpl_test_zero(cpl_propertylist_get_bool(plist, keys[1]));
    cpl_test_eq(cpl_propertylist_get_int(plist, keys[2]), -1);
    cpl_test_eq(cpl_propertylist_get_long(plist, keys[3]), 1);

    fval1 = 9.87654321;
    cpl_test_noneq(fval0, fval1);
    fval2 = cpl_propertylist_get_float(plist, keys[4]);
    cpl_test_abs(fval1, fval2, 0.0);

    dval1 = -9.87654321;
    cpl_test_noneq(dval0, dval1);
    dval2 = cpl_propertylist_get_double(plist, keys[5]);
    cpl_test_abs(dval1, dval2, 0.0);

    cpl_test_eq_string(cpl_propertylist_get_string(plist, keys[6]),
                       comments[0]);

    /*
     * Test 9: Erase all elements from the property list and verify that
     *         the list is empty.
     */

    cpl_propertylist_empty(plist);

//    cpl_test_assert(cpl_propertylist_is_empty(plist));
    cpl_test(cpl_propertylist_is_empty(plist));
   cpl_test_zero(cpl_propertylist_get_size(plist));

    cpl_propertylist_delete(plist);

     /*
     * Test 10: Write a propertylist to a FITS file
     *          and save it in disk
     */

    code = cpl_propertylist_update_string(plist2,"ESO OBS DID", "BLABLA");
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_update_string(plist2,"ESO OBS OBSERVER", "NAME");
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_update_string(plist2,"ESO OBS PI-COI NAME", "NAME");
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_update_string(plist2,"ESO INS GRAT NAME", "DODO");
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_update_string(plist2,"ESO PRO CATG", "PRODUCT");
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_update_int(plist2,"ESO TPL NEXP", 4);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_update_int(plist2, "TOFITS", 6);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_set_comment(plist2, "TOFITS", "Test TOFITS function");
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    (void)remove(filename3);
    status = 0;

    fits_create_diskfile(&fptr, filename3, &status);
    cpl_test_zero(status);

    /*
    * Create simple HDU
    */
    fits_create_img(fptr, 8, 0, naxes, &status);
    cpl_test_zero(status);

    /*Save a propertylist to a FITS file*/
    code = cpl_propertylist_to_fitsfile(fptr, plist2, NULL);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    fits_close_file(fptr, &status);
    cpl_test_zero(status);
    cpl_test_fits(filename3);

    cpl_propertylist_delete(plist2);

    /* Test 11: Load back the propertylist from a FITS file using CFITSIO
    *          and compare it with a plist loaded using CPL.
    *         Retrieve a property usign its name
    */

    plist2 = cpl_propertylist_load(filename3, 0);
    cpl_test_nonnull(plist2);

    /* Try to load a non-existent extension and
       check that the proper error code is set */
    plist = cpl_propertylist_load(filename3, 10);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(plist);

    pro2 = cpl_propertylist_get_property(plist2, "TOFITS");
    cpl_test_nonnull(pro2);
    cpl_test_eq_string("TOFITS", cpl_property_get_name(pro2));
    cpl_test_eq(cpl_property_get_int(pro2), 6);

    status = 0;
    fits_open_diskfile(&_fptr, filename3, READONLY, &status);
    cpl_test_zero(status);

    _plist = cpl_propertylist_from_fitsfile(_fptr);
    cpl_test_nonnull(_plist);

    /* Check that the property of both lists are the same */
    pro = cpl_propertylist_get_property(_plist, "TOFITS");
    cpl_test_nonnull(pro);

    cpl_test_eq(cpl_property_get_int(pro), cpl_property_get_int(pro2));
    //    cpl_propertylist_delete(plist2);

    fits_close_file(_fptr, &status);
    cpl_test_zero(status);
    cpl_test_fits(filename3);

    /*
    * Test 12: Compare the properties of both lists
    */

    cpl_test_eq(cpl_propertylist_get_size(plist2),
              cpl_propertylist_get_size(_plist));

    for (i = 0; i < cpl_propertylist_get_size(plist2); i++) {
        const cpl_property *p = cpl_propertylist_get_const(plist2, i);
        cpl_property *_p = cpl_propertylist_get(_plist, i);

        cpl_test_eq_string(cpl_property_get_name(p),
                           cpl_property_get_name(_p));
        cpl_test_eq_string(cpl_property_get_comment(p),
                           cpl_property_get_comment(_p));
        cpl_test_eq(cpl_property_get_type(p), cpl_property_get_type(_p));

        switch (cpl_property_get_type(p)) {
        case CPL_TYPE_BOOL:
            cpl_test_eq(cpl_property_get_bool(p), cpl_property_get_bool(_p));
            break;

        case CPL_TYPE_INT:
            cpl_test_eq(cpl_property_get_int(p), cpl_property_get_int(_p));
            break;

        case CPL_TYPE_DOUBLE:
            cpl_test_eq(cpl_property_get_double(p),
                      cpl_property_get_double(_p));
            break;

        case CPL_TYPE_STRING:
            cpl_test_eq_string(cpl_property_get_string(p),
                             cpl_property_get_string(_p));
            break;

        default:
            /* This point should never be reached */
            cpl_test(0);
            break;
        }
    }

    cpl_propertylist_delete(_plist);
    plist = cpl_propertylist_duplicate(plist2);

    cpl_propertylist_delete(plist2);

    /*
     * Test 13: Copy all properties matching a given pattern from one
     *          property list to another.
     */

    _plist = cpl_propertylist_new();

    cpl_propertylist_copy_property_regexp(_plist, plist, "^ESO .*", 0);

    cpl_test(cpl_propertylist_has(_plist, "ESO OBS DID"));
    cpl_test(cpl_propertylist_has(_plist, "ESO OBS OBSERVER"));
    cpl_test(cpl_propertylist_has(_plist, "ESO OBS PI-COI NAME"));
    cpl_test(cpl_propertylist_has(_plist, "ESO INS GRAT NAME"));
    cpl_test(cpl_propertylist_has(_plist, "ESO PRO CATG"));
    cpl_test(cpl_propertylist_has(_plist, "ESO TPL NEXP"));

    cpl_propertylist_empty(_plist);
    cpl_test(cpl_propertylist_is_empty(_plist));

    cpl_propertylist_copy_property_regexp(_plist, plist, "^ESO .*", 1);
    cpl_test_zero(cpl_propertylist_has(_plist, "ESO OBS DID"));


    /*
     * Test 14a: Erase all properties matching the given pattern from the
     *          property list.
     */

    cpl_propertylist_empty(_plist);
    cpl_test(cpl_propertylist_is_empty(_plist));

    code = cpl_propertylist_copy_property_regexp(NULL, plist, ".", 0);
    cpl_test_eq_error(code, CPL_ERROR_NULL_INPUT);

    code = cpl_propertylist_copy_property_regexp(_plist, NULL, ".", 0);
    cpl_test_eq_error(code, CPL_ERROR_NULL_INPUT);

    code = cpl_propertylist_copy_property_regexp(_plist, plist, NULL, 0);
    cpl_test_eq_error(code, CPL_ERROR_NULL_INPUT);

    code = cpl_propertylist_copy_property_regexp(_plist, plist, ")|(", 1);
    cpl_test_eq_error(code, CPL_ERROR_ILLEGAL_INPUT);

    code = cpl_propertylist_copy_property_regexp(_plist, plist, "^ESO .*", 0);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_eq(cpl_propertylist_get_size(_plist), 6);

    value  = cpl_propertylist_erase_regexp(_plist, "^ESO OBS .*", 0);
    cpl_test_eq(value, 3);
    cpl_test_eq(cpl_propertylist_get_size(_plist), 3);

    value = cpl_propertylist_erase_regexp(_plist, "ESO TPL NEXP", 0);
    cpl_test_eq(value, 1);
    cpl_test_eq(cpl_propertylist_get_size(_plist), 2);

    /*
    *  Test 14b: Erase a non-existing property and check return.
                 When regexp is null, the return should be -1
    */
    value = cpl_propertylist_erase_regexp(_plist, "ESO TPL NEXP", 0);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_zero(value);

    value = cpl_propertylist_erase_regexp(_plist, ")|(", 1);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_eq(value, -1);

    value = cpl_propertylist_erase_regexp(NULL, ".", 1);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(value, -1);

    value = cpl_propertylist_erase_regexp(_plist, NULL, 0);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(value, -1);

    cpl_propertylist_delete(_plist);
    cpl_propertylist_delete(plist);


    /*
     * Test 15: Create a property list from a file. Only properties matching
     *          the given pattern are loaded.
     */

    plist = cpl_propertylist_load_regexp(filename3, 0, "^ESO .*", 0);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(plist);
    cpl_test_zero(cpl_propertylist_is_empty(plist));
    cpl_test_eq(cpl_propertylist_get_size(plist), 6);
    cpl_test(cpl_propertylist_has(plist, "ESO OBS DID"));
    cpl_test(cpl_propertylist_has(plist, "ESO OBS OBSERVER"));
    cpl_test(cpl_propertylist_has(plist, "ESO OBS PI-COI NAME"));
    cpl_test(cpl_propertylist_has(plist, "ESO INS GRAT NAME"));
    cpl_test(cpl_propertylist_has(plist, "ESO PRO CATG"));
    cpl_test(cpl_propertylist_has(plist, "ESO TPL NEXP"));

    cpl_propertylist_delete(plist);

    /*
     * Test 15b: Failure tests of cpl_propertylist_load{,_regexp}().
     */

    remove("no.fits");

    value = cpl_fits_count_extensions(filename3);
    cpl_test_leq(0, value);

    plist = cpl_propertylist_load(NULL, 0);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load_regexp(NULL, 0, ".", 0);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load_regexp(filename3, 0, NULL, 0);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load("no.fits", 0);
    cpl_test_error(CPL_ERROR_FILE_IO);
    cpl_test_null(plist);

    plist = cpl_propertylist_load_regexp("no.fits", 0, ".", 0);
    cpl_test_error(CPL_ERROR_FILE_IO);
    cpl_test_null(plist);

    plist = cpl_propertylist_load("/dev/null", 0);
    cpl_test_error(CPL_ERROR_BAD_FILE_FORMAT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load_regexp("/dev/null", 0, ".", 0);
    cpl_test_error(CPL_ERROR_BAD_FILE_FORMAT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load(filename3, -1);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load_regexp(filename3, -1, ".", 0);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load(filename3, -1);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load_regexp(filename3, -1, ".", 0);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(plist);

    plist = cpl_propertylist_load(filename3, value + 1);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(plist);

    plist = cpl_propertylist_load_regexp(filename3, value + 1, ".", 0);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(plist);

    plist = cpl_propertylist_load_regexp(filename3, 0, "\\", 0);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(plist);


    /*
     * Test 16: Append a property list to another.
     */

    plist = cpl_propertylist_new();
    _plist = cpl_propertylist_new();

    cpl_propertylist_append_char(plist, keys[0], 'a');
    cpl_propertylist_set_comment(plist, keys[0], comments[0]);

    cpl_propertylist_append_bool(plist, keys[1], 1);
    cpl_propertylist_set_comment(plist, keys[1], comments[1]);

    cpl_propertylist_append_int(plist, keys[2], -1);
    cpl_propertylist_set_comment(plist, keys[2], comments[2]);

    cpl_propertylist_append_long(plist, keys[3], 32768);
    cpl_propertylist_set_comment(plist, keys[3], comments[3]);

    cpl_propertylist_append_float(_plist, keys[4], fval0);
    cpl_propertylist_set_comment(_plist, keys[4], comments[4]);

    cpl_propertylist_append_double(_plist, keys[5], dval0);
    cpl_propertylist_set_comment(_plist, keys[5], comments[5]);

    cpl_propertylist_append_string(_plist, keys[6], comments[6]);
    cpl_propertylist_set_comment(_plist, keys[6], comments[6]);

    cpl_propertylist_append_float_complex(_plist, keys[7], zf0);
    cpl_propertylist_set_comment(_plist, keys[7], comments[7]);

    cpl_propertylist_append_double_complex(_plist, keys[8], zd0);
    cpl_propertylist_set_comment(_plist, keys[8], comments[8]);

    cpl_test_zero(cpl_propertylist_is_empty(plist));
    cpl_test_eq(cpl_propertylist_get_size(plist), 4);

    cpl_test_zero(cpl_propertylist_is_empty(_plist));
    cpl_test_eq(cpl_propertylist_get_size(_plist), 5);

    code = cpl_propertylist_append(plist, _plist);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    cpl_test_zero(cpl_propertylist_is_empty(plist));
    cpl_test_eq(cpl_propertylist_get_size(plist), nprops);

    cpl_test_zero(cpl_propertylist_is_empty(_plist));
    cpl_test_eq(cpl_propertylist_get_size(_plist), 5);

    for (i = 0; i < cpl_propertylist_get_size(plist); i++) {
        cpl_property *p = cpl_propertylist_get(plist, i);

        cpl_test_eq_string(cpl_property_get_name(p), keys[i]);
        cpl_test_eq_string(cpl_property_get_comment(p), comments[i]);
        cpl_test_eq(cpl_property_get_type(p), types[i]);

        cpl_test(cpl_propertylist_has(plist, keys[i]));
        cpl_test_eq_string(cpl_propertylist_get_comment(plist, keys[i]),
                           comments[i]);
        cpl_test_eq(cpl_propertylist_get_type(plist, keys[i]), types[i]);
    }

    /*
    * Test 17: Sequentially access the elements of a property list
    */

    /* First element */

    /*    pro = cpl_propertylist_get_first(plist);
    cpl_test_nonnull(pro);

    cpl_test_eq_string(cpl_property_get_name(pro), keys[0]);
    cpl_test_eq_string(cpl_property_get_comment(pro), comments[0]);
    cpl_test_eq(cpl_property_get_type(pro), types[0]);*/

    /*
    * Test 18: Loop through all elements except the first,
    *           and check the results.
    */
    /*    for (i = 1; i < cpl_propertylist_get_size(plist); i++) {

        cpl_property *p = cpl_propertylist_get_next(plist);

        cpl_test_eq_string(cpl_property_get_name(p), keys[i]);
        cpl_test_eq_string(cpl_property_get_comment(p), comments[i]);
        cpl_test_eq(cpl_property_get_type(p), types[i]);

        cpl_test_assert(cpl_propertylist_has(plist, keys[i]));
        cpl_test_zero(strcmp(cpl_propertylist_get_comment(plist, keys[i]),
                          comments[i]));
        cpl_test_eq(cpl_propertylist_get_type(plist, keys[i]), types[i]);
     }*/

     /*
    * Test 19: Try to access the next element and see
    *           if the end of list is handled properly.
    */

    /*    pro = cpl_propertylist_get_next(plist);
    cpl_test_eq(pro, NULL);*/


/*
     * Test 20: Create a FITS header using a list containing a property with
     *          a name of length 80 characters (the length of a FITS card)
     */

    cpl_propertylist_empty(plist);

    code = cpl_propertylist_append_string(plist, longname, comments[6]);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    cpl_propertylist_delete(_plist);
    _plist = NULL;

    cpl_propertylist_delete(plist);
    plist = NULL;

    /*
    * Test 21:  Sort a property list
    */

    plist = cpl_propertylist_new();

    cpl_test_nonnull(plist);
    cpl_test(cpl_propertylist_is_empty(plist));
    cpl_test_zero(cpl_propertylist_get_size(plist));


    for(i = 0; i < 1; i++) {

        cpl_propertylist_append_int(plist, "FDEBACGH", -1);
        cpl_propertylist_set_comment(plist, "FDEBACGH", comments[2]);

        cpl_propertylist_append_char(plist, "CDEFGB1", 'a');
        cpl_propertylist_set_comment(plist, "CDEFGB1", comments[0]);

        cpl_propertylist_append_bool(plist, "ABCDEFGH", 1);
        cpl_propertylist_set_comment(plist, "ABCDEFGH", comments[1]);

        cpl_propertylist_append_float(plist, "ZZZGDBCA", fval0);
        cpl_propertylist_set_comment(plist, "ZZZGDBCA", comments[4]);

        cpl_propertylist_append_string(plist, "BBBBB2", comments[6]);
        cpl_propertylist_set_comment(plist, "BBBBB2", comments[6]);

        cpl_propertylist_append_long(plist, "HISTORY", 32768);
        cpl_propertylist_set_comment(plist, "HISTORY", comments[3]);

        cpl_propertylist_append_int(plist, "HISTORY", 0);
        cpl_propertylist_set_comment(plist, "HISTORY", comments[3]);

        cpl_propertylist_append_double(plist, "HISTORY", dval0);
        cpl_propertylist_set_comment(plist, "HISTORY", comments[5]);

    }

    cpl_test_zero(cpl_propertylist_is_empty(plist));
    cpl_test_eq(cpl_propertylist_get_size(plist), 8);

    code = cpl_propertylist_sort(plist, cpl_test_property_compare_name);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    cpl_propertylist_delete(plist);

    /*
    * Test 21a:  Save a propertylist which contains a HISTORY
    *            field without a value
    */

    plist = cpl_propertylist_new();
    cpl_propertylist_append_string(plist, "HISTORY", "HELLO");
    cpl_propertylist_append_string(plist, "HISTORY", "WORLD");
    cpl_propertylist_append_string(plist, "HISTORY", "");
    cpl_propertylist_append_string(plist, "HISTORY", "HELLO");
    cpl_propertylist_append_string(plist, "HISTORY", "AGAIN");

    cpl_propertylist_append_string(plist, "COMMENT", "HELLO");
    cpl_propertylist_append_string(plist, "COMMENT", "WORLD");
    cpl_propertylist_append_string(plist, "COMMENT", "");
    cpl_propertylist_append_string(plist, "COMMENT", "HELLO");
    cpl_propertylist_append_string(plist, "COMMENT", "AGAIN");
    cpl_test_error(CPL_ERROR_NONE);

    t = cpl_table_new(1);
    code = cpl_table_save(t, plist, NULL, BASE "_21a.fits", CPL_IO_CREATE);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_21a.fits");


    cpl_propertylist_delete(plist);
    cpl_table_delete(t);

    /*
     * Test 21b: load plist from an image and sort it.
     * This test was compared with result from cpl_propertylist_sort().
     * Results are expected to be the same.
     */

    plist = cpl_propertylist_load(BASE "_image.fits", 0);

    if (plist == NULL) {
        cpl_test_error(CPL_ERROR_FILE_IO);
    } else {

        /* Inactive per default */

        code = cpl_propertylist_sort(plist, cpl_test_property_compare_name);
        cpl_test_eq_error(code, CPL_ERROR_NONE);


        code = cpl_image_save(NULL, BASE "_sorted.fits", CPL_TYPE_UCHAR,
                              plist, CPL_IO_CREATE);
        cpl_test_eq_error(code, CPL_ERROR_NONE);
        cpl_test_fits(BASE "_sorted.fits");
        cpl_test_zero(remove(BASE "_sorted.fits"));

        cpl_propertylist_delete(plist);
    }


    /*
     * Test 21c:  Save a very long property list into a FITS file
     */

    plist = cpl_propertylist_new();

    cpl_test_nonnull(plist);
    cpl_test(cpl_propertylist_is_empty(plist));
    cpl_test_zero(cpl_propertylist_get_size(plist));

    bench_size = do_bench ? 100000 : 10000;

    for(i = 0; i < bench_size; i++) {
        cpl_propertylist_append_float(plist, keys[4], fval0);
        cpl_propertylist_set_comment(plist, keys[4], comments[4]);

        cpl_propertylist_append_int(plist, keys[2], -1);
        cpl_propertylist_set_comment(plist, keys[2], comments[2]);

        cpl_propertylist_append_char(plist, keys[0], 'a');
        cpl_propertylist_set_comment(plist, keys[0], comments[0]);

        cpl_propertylist_append_bool(plist, keys[1], 1);
        cpl_propertylist_set_comment(plist, keys[1], comments[1]);

        cpl_propertylist_append_string(plist, keys[6], comments[6]);
        cpl_propertylist_set_comment(plist, keys[6], comments[6]);

        cpl_propertylist_append_long(plist, keys[3], 32768);
        cpl_propertylist_set_comment(plist, keys[3], comments[3]);

        cpl_propertylist_append_double(plist, keys[5], dval0);
        cpl_propertylist_set_comment(plist, keys[5], comments[5]);

        cpl_propertylist_append_float_complex(plist, keys[7], fval0);
        cpl_propertylist_set_comment(plist, keys[7], comments[7]);

        cpl_propertylist_append_double_complex(plist, keys[8], dval0);
        cpl_propertylist_set_comment(plist, keys[8], comments[8]);

    }

    cpl_test_zero(cpl_propertylist_is_empty(plist));
    cpl_test_eq(cpl_propertylist_get_size(plist), nprops * bench_size);

    remove(filename3);
    status = 0;

    fits_create_diskfile(&fptr, filename3, &status);
    cpl_test_zero(status);

    fits_create_img(fptr, 8, 0, naxes, &status);
    cpl_test_zero(status);

    code = cpl_propertylist_to_fitsfile(fptr, plist, NULL);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    fits_close_file(fptr, &status);
    cpl_test_zero(status);
    cpl_test_fits(filename3);

    cpl_propertylist_delete(plist);

    /*
     *   Test 22: Save a property list into a FITS file using
     *            cpl_propertylist_save()
     */

    plist = cpl_propertylist_new();
    cpl_propertylist_append_char(plist, keys[0], 'a');
    cpl_propertylist_set_comment(plist, keys[0], comments[0]);

    cpl_propertylist_append_bool(plist, keys[1], 1);
    cpl_propertylist_set_comment(plist, keys[1], comments[1]);

    cpl_propertylist_append_int(plist, keys[2], -1);
    cpl_propertylist_set_comment(plist, keys[2], comments[2]);

    cpl_propertylist_append_long(plist, keys[3], 32768);
    cpl_propertylist_set_comment(plist, keys[3], comments[3]);

    cpl_propertylist_append_float(plist, keys[4], fval0);
    cpl_propertylist_set_comment(plist, keys[4], comments[4]);

    cpl_propertylist_append_double(plist, keys[5], dval0);
    cpl_propertylist_set_comment(plist, keys[5], comments[5]);

    cpl_propertylist_append_string(plist, keys[6], comments[6]);
    cpl_propertylist_set_comment(plist, keys[6], comments[6]);

    cpl_propertylist_append_float_complex(plist, keys[7], zf0);
    cpl_propertylist_set_comment(plist, keys[7], comments[7]);

    cpl_propertylist_append_double_complex(plist, keys[8], zd0);
    cpl_propertylist_set_comment(plist, keys[8], comments[8]);

    cpl_test_zero(cpl_propertylist_is_empty(plist));
    cpl_test_eq(cpl_propertylist_get_size(plist), nprops);


    code = cpl_propertylist_save(plist, BASE "_22.fits", CPL_IO_CREATE);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_22.fits");

    code = cpl_propertylist_save(plist, BASE "_22.fits", CPL_IO_APPEND);
    cpl_test_eq_error(code, CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_fits(BASE "_22.fits");

    cpl_propertylist_update_string(plist, keys[6], "updated string");
    code = cpl_propertylist_save(plist, BASE "_22.fits", CPL_IO_EXTEND);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_22.fits");

    cpl_propertylist_delete(plist);

    /*
     * Save a NULL property list to an extension
     */

    remove(BASE "_null.fits");
    plist = cpl_propertylist_new();
    code = cpl_propertylist_save(plist, BASE "_null.fits", CPL_IO_CREATE);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_null.fits");

    code = cpl_propertylist_save(NULL, BASE "_null.fits", CPL_IO_EXTEND);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_null.fits");
    cpl_propertylist_delete(plist);

    /*
     * Check that XTENSION is set to IMAGE
     */

    plist = cpl_propertylist_load(BASE "_null.fits",1);
    cpl_test_eq_string("IMAGE", cpl_propertylist_get_string(plist, "XTENSION"));
    cpl_propertylist_delete(plist);

    /*
    *   Test 23: Save a property list to a FITS file using
    *            a regexp to remove keys
    */
    remove(filename3);
    status = 0;

    fits_create_diskfile(&fptr, filename3, &status);
    cpl_test_zero(status);

    /*
    * Create simple HDU
    */
    fits_create_img(fptr, 8, 0, naxes, &status);
    cpl_test_zero(status);

    plist = cpl_propertylist_new();
    cpl_propertylist_append_char(plist, keys[0], 'a');
    cpl_propertylist_set_comment(plist, keys[0], comments[0]);

    cpl_propertylist_append_bool(plist, keys[1], 1);
    cpl_propertylist_set_comment(plist, keys[1], comments[1]);

    cpl_propertylist_append_int(plist, keys[2], -1);
    cpl_propertylist_set_comment(plist, keys[2], comments[2]);

    cpl_propertylist_append_long(plist, "NBXIS1", 32768);
    cpl_propertylist_set_comment(plist, "NBXIS1", comments[3]);

    cpl_propertylist_append_float(plist, keys[4], fval0);
    cpl_propertylist_set_comment(plist, keys[4], comments[4]);

    cpl_propertylist_append_double(plist, keys[5], dval0);
    cpl_propertylist_set_comment(plist, keys[5], comments[5]);

    cpl_propertylist_append_string(plist, "HIERARCH ESO TEST1", "One string without any comment");
    cpl_propertylist_append_string(plist, "HIERARCH ESO TEST2", "Two string without any comment");


    cpl_propertylist_append_float_complex(plist, keys[7], zf0);
    cpl_propertylist_set_comment(plist, keys[7], comments[7]);

    cpl_propertylist_append_double_complex(plist, keys[8], zd0);
    cpl_propertylist_set_comment(plist, keys[8], comments[8]);

    cpl_propertylist_append_long(plist, "NBXIS2", 32769);
    cpl_propertylist_set_comment(plist, "NBXIS2", comments[3]);

    cpl_propertylist_append_string(plist, "HISTORY", "One history string without any comment");
    cpl_propertylist_append_string(plist, "HISTORY", "Two history string without any comment");
    cpl_propertylist_append_string(plist, "HISTORY", "Three history string without any comment");

#ifdef CPL_FITS_TEST_NON_STRING_HISTORY_COMMENT
    /* This will produce a pretty strange FITS header */

    /* Check handling of HISTORY of non-string type */
    cpl_propertylist_append_float(plist, "HISTORY", 4.0);

    /* Check handling of COMMENT of non-string type */
    cpl_propertylist_append_float(plist, "COMMENT", 5.0);

    cpl_test_eq(cpl_propertylist_get_size(plist), 16);
#else
    cpl_test_eq(cpl_propertylist_get_size(plist), 14);
#endif
    cpl_test_zero(cpl_propertylist_is_empty(plist));

    /* Remove the following keys from the saving
       HIERARCH ESO |NBXIS1|NBXIS2|HISTORY
    */

    code = cpl_propertylist_to_fitsfile(fptr, plist, ")|(");
    cpl_test_eq_error(code, CPL_ERROR_ILLEGAL_INPUT);

    code = cpl_propertylist_to_fitsfile(fptr, plist, to_rm);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    /* Pass a null regexp to the function*/
    code = cpl_propertylist_to_fitsfile(fptr, plist, NULL);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    /* Pass an empty property list to the function */
    cpl_propertylist_empty(plist);
    code = cpl_propertylist_to_fitsfile(fptr, plist, to_rm);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_to_fitsfile(fptr, NULL, to_rm);
    cpl_test_eq_error(code, CPL_ERROR_NULL_INPUT);

    code = cpl_propertylist_to_fitsfile(NULL, plist, to_rm);
    cpl_test_eq_error(code, CPL_ERROR_NULL_INPUT);

    fits_close_file(fptr, &status);
    cpl_test_zero(status);
    cpl_test_fits(filename3);

    cpl_propertylist_delete(plist);

    plist = cpl_propertylist_load(filename3, 0);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(plist);

    /*
    *   Test 24: Create a property list with compressed
    *            keywords and save it in disk and in a table.
    *            Check if compressed
    *            keywords are removed after saving it.
    */

    remove(BASE "_24.fits");
    cpl_propertylist_empty(plist);
    cpl_propertylist_append_string(plist,"ZTENSION", "COMPRESSED_EMPTY");
    cpl_propertylist_append_long(plist, "ZNAXIS", 2);
    cpl_propertylist_set_comment(plist, "ZNAXIS", "compressed NAXIS");
    cpl_propertylist_append_long(plist, "ZNAXIS1", 2048);
    cpl_propertylist_set_comment(plist, "ZNAXIS1", "compressed NAXIS1");
    cpl_propertylist_append_long(plist, "ZNAXIS2", 1024);
    cpl_propertylist_set_comment(plist, "ZNAXIS2", "compressed NAXIS2");
    cpl_test_zero(cpl_propertylist_is_empty(plist));
    cpl_test_eq(cpl_propertylist_get_size(plist), 4);


    /* Save it in disk and check if compression keywords are removed */
    code = cpl_propertylist_save(plist, BASE "_24.fits", CPL_IO_CREATE);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_24.fits");
    code = cpl_propertylist_save(plist, BASE "_24.fits", CPL_IO_EXTEND);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_24.fits");
//    cpl_propertylist_delete(plist);

    /* primary header */
    plist2 = cpl_propertylist_load(BASE "_24.fits", 0);
    cpl_test_zero(cpl_propertylist_has(plist2, "ZTENSION"));
    cpl_propertylist_delete(plist2);

    /* extension header */
    plist2 = cpl_propertylist_load(BASE "_24.fits", 1);
    cpl_test_zero(cpl_propertylist_has(plist2, "ZNAXIS"));
    cpl_propertylist_delete(plist2);


    (void)remove(BASE "_24.fits");
    t = cpl_table_new(1);

    /* Compressed keywords should be removed when saving table*/
    code = cpl_table_save(t, plist, NULL, BASE "_24.fits", CPL_IO_CREATE);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_24.fits");
    code = cpl_table_save(t, NULL, plist, BASE "_24.fits", CPL_IO_EXTEND);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_fits(BASE "_24.fits");
    cpl_propertylist_delete(plist);

    /* Check if compressed keywords were removed from headers*/
    plist = cpl_propertylist_load(BASE "_24.fits",1);
    cpl_test_zero(cpl_propertylist_has(plist, "ZTENSION"));
    cpl_test_zero(cpl_propertylist_has(plist, "ZNAXIS"));
    cpl_test_zero(cpl_propertylist_has(plist, "ZNAXIS1"));
    cpl_test_zero(cpl_propertylist_has(plist, "ZNAXIS2"));
    cpl_propertylist_delete(plist);

    plist = cpl_propertylist_load(BASE "_24.fits",2);
    cpl_test_zero(cpl_propertylist_has(plist, "ZTENSION"));
    cpl_test_zero(cpl_propertylist_has(plist, "ZNAXIS"));
    cpl_test_zero(cpl_propertylist_has(plist, "ZNAXIS1"));
    cpl_test_zero(cpl_propertylist_has(plist, "ZNAXIS2"));

    /* Test propertylist dump */
    cpl_propertylist_dump(plist, stream);

    cpl_table_delete(t);
    cpl_propertylist_delete(plist);

    /*
     *  Test 25: Test the new error message function
     */

    plist = cpl_propertylist_load(BASE "_null.fits", 2);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);

    /*
     *  Test 26: Verify that NANs cannot be saved
     */

    cpl_propertylist_delete(plist); /* Could be non-NULL */
    plist = cpl_propertylist_new();

    cpl_test_zero(cpl_propertylist_append_double(plist, "NANKEY", 0.0/0.0));

    /* Allow verification that the value is really NAN */
    cpl_propertylist_dump(plist, stream);

    code = cpl_propertylist_save(plist, BASE "_26.fits", CPL_IO_CREATE);
    cpl_test_eq_error(code, CPL_ERROR_ILLEGAL_INPUT);

    cpl_propertylist_empty(plist);

    /*
     *  Test 26A: Verify that a (non-FITS) NaN cannot be loaded
     */
    cpl_test_zero(cpl_propertylist_append_string(plist, "BADKEY", BADNUM));

    cpl_test_zero(cpl_propertylist_save(plist, BASE "_26.fits",
                                        CPL_IO_CREATE));
    cpl_test_fits(BASE "_26.fits");
    cpl_propertylist_delete(plist);

#ifdef CPL_NO_PERL
    plist = cpl_propertylist_load(BASE "_26.fits", 0);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(cpl_propertylist_get_type(plist, "BADKEY"), CPL_TYPE_STRING);
    cpl_propertylist_delete(plist);
#else

    /* Need to close file due to non CPL access */
    cpl_test_zero(cpl_io_fits_close(BASE "_26.fits", &status));
    cpl_test_zero(status);

    cpl_test_eq(strlen(BADNUM), 3); /* For the below substitution to work */
    /* Remove quotes from value in FITS card */
    if (system("perl -pi -e 's/." BADNUM "     ./ " BADNUM "      /' "
               BASE "_26.fits") == 0) {

        /* The file is no longer valid FITS, and should not be loadable */
        plist = cpl_propertylist_load(BASE "_26.fits", 0);
        cpl_test_error(CPL_ERROR_BAD_FILE_FORMAT);
        cpl_test_null(plist);
    }

#endif

    /*
     *  Test 26B: Verify that a complex number cannot be loaded
     */

    plist = cpl_propertylist_new();
    code = cpl_propertylist_append_string(plist, "BADKEY", COMPLEXVAL);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_zero(cpl_propertylist_save(plist, BASE "_26.fits",
                                        CPL_IO_CREATE));
    cpl_test_fits(BASE "_26.fits");
    cpl_propertylist_delete(plist);

#ifdef CPL_NO_PERL
    plist = cpl_propertylist_load(BASE "_26.fits", 0);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(cpl_propertylist_get_type(plist, "BADKEY"), CPL_TYPE_STRING);
#else
    /* Need to close file due to non CPL access */
    cpl_test_zero(cpl_io_fits_close(BASE "_26.fits", &status));
    cpl_test_zero(status);

    cpl_test_eq(strlen(COMPLEXVAL), 7); /* For the below substitution to work */
    /* Replace quotes with parenthesis to form a complex number */
    cpl_test_zero(system("perl -pi -e 's/\\047" COMPLEXVAL " \\047/("
                         COMPLEXVAL ") /' " BASE "_26.fits"));

    /* The file should still be valid FITS, but no longer loadable by CPL */
    cpl_test_fits(BASE "_26.fits");

    plist = cpl_propertylist_load(BASE "_26.fits", 0);
    cpl_test_error(CPL_ERROR_BAD_FILE_FORMAT);
    cpl_test_null(plist);

    /*  - Except if the complex card is ignored */
    plist = cpl_propertylist_load_regexp(BASE "_26.fits", 0, "BADKEY", 1);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_zero(cpl_propertylist_has(plist, "BADKEY"));
#endif
    cpl_propertylist_delete(plist);


    /*
     * Test 27: Append a property to a property list
     */

    plist = cpl_propertylist_new();

    /* Catch error code */
    code = cpl_propertylist_append_property(plist, NULL);
    cpl_test_eq_error(code, CPL_ERROR_NULL_INPUT);

    for(j = 0; j < 6; j++) {
        char *name = cpl_sprintf("HIERARCH ESO TEST%d APPEND", j);
        cpl_property *p = cpl_property_new(name, CPL_TYPE_INT);
        cpl_free(name);

        cpl_test_nonnull(p);
        cpl_test_zero(cpl_property_set_int(p, j+2));
        cpl_test_zero(cpl_propertylist_append_property(plist, p));
        cpl_property_delete(p);
    }

    cpl_test_eq(cpl_propertylist_get_size(plist),6);
    cpl_test_eq(cpl_propertylist_get_int(plist,"HIERARCH ESO TEST1 APPEND"),3);

     /*
      * Test 28: Prepend a property to a property list
      */


     /* Catch error code */
    code = cpl_propertylist_prepend_property(plist, NULL);
    cpl_test_eq_error(code, CPL_ERROR_NULL_INPUT);

    for(j = 0; j < 3; j++) {
        char *name = cpl_sprintf("HIERARCH ESO TEST%d PREPEND",j);
        cpl_property *p = cpl_property_new(name,CPL_TYPE_STRING);
           cpl_free(name);
        cpl_test_nonnull(p);
        cpl_test_zero(cpl_property_set_string(p,"test prepend"));
        cpl_test_zero(cpl_propertylist_prepend_property(plist,p));
        cpl_property_delete(p);
    }

    cpl_test_eq(cpl_propertylist_get_size(plist),9);


     /* Test 29: Insert a property into a property list */

    code = cpl_propertylist_insert_property(plist, "HIERARCH ESO TEST0 APPEND",
                                            NULL);
    cpl_test_eq_error(code, CPL_ERROR_NULL_INPUT);

    for(j = 0; j < 1; j++){
        const char *name = "HIERARCH ESO TEST0 APPEND";
        cpl_property *p = cpl_property_new("INSERT1", CPL_TYPE_FLOAT);
        cpl_property *p1 = cpl_property_new("INSERT2", CPL_TYPE_FLOAT);
        cpl_test_nonnull(p);
        cpl_test_nonnull(p1);
        cpl_test_zero(cpl_property_set_float(p, 1.0));
        cpl_test_zero(cpl_propertylist_insert_property(plist,name,p));
        cpl_test_zero(cpl_property_set_float(p1, 2.0));
        cpl_test_zero(cpl_propertylist_insert_after_property(plist,name,p1));
        cpl_property_delete(p);
        cpl_property_delete(p1);
    }

    cpl_test_eq(cpl_propertylist_get_size(plist),11);

//    cpl_propertylist_dump(plist,stdout);
    cpl_propertylist_delete(plist);


    /*
     * Test 30: Test the casting on the float and double accessors
     */

    plist = cpl_propertylist_new();
    cpl_test_zero(cpl_propertylist_update_string(plist,"C","String"));
    cpl_test_zero(cpl_propertylist_update_int(plist,"I",8));
    f0 = 235.89;
    cpl_test_zero(cpl_propertylist_update_float(plist,"F",f0));
    d0 = 1234.9994048;
    cpl_test_zero(cpl_propertylist_update_double(plist,"D",d0));
    cpl_test_zero(cpl_propertylist_save(plist,BASE "_30.fits",CPL_IO_CREATE));
    cpl_test_fits(BASE "_30.fits");

    cpl_propertylist_delete(plist);

    plist = cpl_propertylist_load(BASE "_30.fits",0);
    cpl_test_nonnull(plist);

    /* The float keyword is casted to double when loaded from disk */
    t1 = cpl_propertylist_get_type(plist,"F");
    cpl_test_eq(t1,CPL_TYPE_DOUBLE);

    /* Test that the casting in cpl works */
    f1 = cpl_propertylist_get_float(plist, "F");
    cpl_test_error(CPL_ERROR_NONE);

    d1 = cpl_propertylist_get_double(plist, "F");
    cpl_test_error(CPL_ERROR_NONE);

    f2 = cpl_propertylist_get_float(plist,"D");
    cpl_test_error(CPL_ERROR_NONE);

    d2 = cpl_propertylist_get_double(plist, "D");
    cpl_test_error(CPL_ERROR_NONE);

    /* Test the values */
    cpl_test_abs(f0, f1, 0.0);
    cpl_test_abs(d0, d2, 0.0);
    cpl_test_rel(d1, f1, FLT_EPSILON);
    cpl_test_rel(d2, f2, FLT_EPSILON);


    /*
     * Test 31: Verify the truncation of a long string in I/O
     */

    cpl_propertylist_empty(plist);

    code = cpl_propertylist_append_string(plist, "ESO PRO REC1 PIPE ID",
                                          LONGNAME80);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_save(plist, BASE "_31.fits", CPL_IO_CREATE);
    cpl_test_fits(BASE "_31.fits");

    plist2 = cpl_propertylist_load(BASE "_31.fits", 0);
    cpl_test_eq_error(code, CPL_ERROR_NONE);
    cpl_test_nonnull(plist2);

    strval = cpl_propertylist_get_string(plist2, "ESO PRO REC1 PIPE ID");
    cpl_test_nonnull(strval);

    if (strval != NULL) {
        value = strlen(strval);
        cpl_test_leq(value, 80);

        cpl_test_zero(strncmp(strval, LONGNAME80, value));
    }

    cpl_propertylist_delete(plist2);
    cpl_propertylist_empty(plist);

    /*
     * Test 32: Verify the handling of setting a string property to NULL
     */

    code = cpl_propertylist_update_string(plist, keys[6], NULL);
    cpl_test_eq_error(CPL_ERROR_NULL_INPUT, code);

    code = cpl_propertylist_append_string(plist, keys[6], comments[6]);
    cpl_test_eq_error(code, CPL_ERROR_NONE);

    code = cpl_propertylist_set_string(plist, keys[6], NULL);
    cpl_test_eq_error(CPL_ERROR_NULL_INPUT, code);

    code = cpl_propertylist_insert_string(plist, keys[6], keys[6 + nprops],
                                          NULL);
    cpl_test_eq_error(CPL_ERROR_NULL_INPUT, code);

    code = cpl_propertylist_insert_after_string(plist, keys[6],
                                          keys[6 + nprops], NULL);
    cpl_test_eq_error(CPL_ERROR_NULL_INPUT, code);

    code = cpl_propertylist_prepend_string(plist, keys[6], NULL);
    cpl_test_eq_error(CPL_ERROR_NULL_INPUT, code);



    cpl_propertylist_delete(plist);

    if (stream != stdout) cpl_test_zero( fclose(stream) );

    /*
     * All tests done
     */


    return cpl_test_end(0);

}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Compare two properties
  @param    p1  First property to compare
  @param    p2  Second property to compare
  @return  An integer less than,  equal to, or greater than zero if
           the name of p1 is found, respectively, to be less than,
           to match, or be greater than that of p2.
  @see cpl_property_get_name(), strcmp()

 */
/*----------------------------------------------------------------------------*/
static int cpl_test_property_compare_name(
        const cpl_property  *   p1,
        const cpl_property  *   p2)
{

    /*
    * Compare the two properties
    */
    return strcmp(cpl_property_get_name(p1),
                  cpl_property_get_name(p2));
}

