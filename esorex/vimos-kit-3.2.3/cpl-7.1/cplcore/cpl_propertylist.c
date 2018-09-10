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
#  include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <regex.h>
#include <complex.h>

#include <fitsio.h>

#include <cxmacros.h>
#include <cxmemory.h>
#include <cxdeque.h>
#include <cxmessages.h>
#include <cxstrutils.h>
#include <cxutils.h>

#include "cpl_io_fits.h"

#include "cpl_errorstate.h"
#include "cpl_error_impl.h"
#include "cpl_tools.h"
#include "cpl_memory.h"
#include "cpl_io.h"
#include "cpl_propertylist_impl.h"


/**
 * @defgroup cpl_propertylist Property Lists
 *
 * This module implements a container for @em properties (see
 * @ref cpl_property) which can be used to store auxiliary values related to
 * another data object, an image or a table for instance. The property values
 * can be set and retrieved by their associated name and properties can be
 * added and removed from the list. The property list container is an ordered
 * sequence of properties.
 *
 * @par Synopsis:
 * @code
 *   #include <cpl_propertylist.h>
 * @endcode
 */

/**@{*/

enum {
    FITS_STDKEY_MAX = 8,
    FITS_SVALUE_MAX = 68
};


/*
 * The property list type.
 */

struct _cpl_propertylist_ {
    cx_deque *properties;
};


/*
 * Regular expresion filter type
 */

struct _cpl_regexp_ {
    regex_t re;
    cxbool invert;
};

typedef struct _cpl_regexp_ cpl_regexp;


/*
 * FITS keyword type mapping for floating point WCS keywords.
 */

struct _cpl_wcskeys_ {
    const cxchar *key;
    cxint sz;
    cxchar type;
};

typedef struct _cpl_wcskeys_ cpl_wcskeys;


static cpl_wcskeys FITS_KEYMAP_WCS[] = {
    {"^CRPIX[0-9]+",     -1, 'F'},
    {"^CRVAL[0-9]+",     -1, 'F'},
    {"^CDELT[0-9]+",     -1, 'F'},
    {"^CRDER[0-9]+",     -1, 'F'},
    {"^CSYER[0-9]+",     -1, 'F'},
    {"^PC[0-9]+_[0-9]+", -1, 'F'},
    {"^PV[0-9]+_[0-9]+", -1, 'F'},
    {"^PS[0-9]+_[0-9]+", -1, 'F'},
    {"EQUINOX",           7, 'F'},
    {"EPOCH",             5, 'F'},
    {"MJD-OBS",           7, 'F'},
    {"LONGPOLE",          8, 'F'},
    {"LATPOLE",           7, 'F'}
};


/*
 * Private methods
 */

inline static cxint
_cpl_wcskeys_find(const cxchar *key, const cpl_wcskeys *keymap,
                  cxint mapsize)
{

    cxint i   = 0;
    cxint pos = mapsize;

    for (i = 0; i < mapsize; ++i) {

        if (keymap[i].sz < 0) {

            cpl_regexp pattern;
            cxint status = regcomp(&pattern.re, keymap[i].key,
                                   REG_EXTENDED | REG_NOSUB);
            if (status) {
                return -1;
            }

            if (regexec(&pattern.re, key, (size_t)0, NULL, 0) == 0) {
                pos = i;
                regfree(&pattern.re);
                break;
            }
            regfree(&pattern.re);

        }
        else {

            if (strncmp(key, keymap[i].key, keymap[i].sz) == 0) {
                pos = i;
                break;
            }

        }

    }

    return pos;

}


inline static cxint
_cpl_propertylist_filter_regexp(cxcptr key, cxcptr filter)
{

    const cpl_regexp *_filter = (const cpl_regexp *)filter;

    if (regexec(&_filter->re, key, (size_t)0, NULL, 0) == REG_NOMATCH)
        return _filter->invert == TRUE ? TRUE : FALSE;

    return _filter->invert == TRUE ? FALSE : TRUE;

}


inline static cxbool
_cpl_propertylist_compare(const cpl_property *property, const char *name)
{
    const cxchar *key = cpl_property_get_name(property);

    return strcmp(key, name) == 0 ? TRUE : FALSE;

}


inline static cxbool
_cpl_propertylist_compare_start(const cpl_property *property,
                                const char *part_name)
{

    const cxchar *key = cpl_property_get_name(property);

    if (strstr(key, part_name) == key)
        return TRUE;

    return FALSE;

}


inline static cxbool
_cpl_propertylist_compare_regexp(const cpl_property *property, cpl_regexp *re)
{

    const cxchar *key = cpl_property_get_name(property);

    return _cpl_propertylist_filter_regexp(key, re);

}


inline static cx_deque_iterator
_cpl_propertylist_find(const cpl_propertylist *self, const char *name)
{

    cx_deque_iterator first, last;
    cpl_property *p;

    first = cx_deque_begin(self->properties);
    last = cx_deque_end(self->properties);

    while (first != last) {
        p = cx_deque_get(self->properties, first);

        if (_cpl_propertylist_compare(p, name))
            break;

        first = cx_deque_next(self->properties, first);
    }

    return first;

}

inline static cpl_property *
_cpl_propertylist_get(const cpl_propertylist *self, const char *name)
{

    cx_deque_iterator pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties))
        return NULL;

    return cx_deque_get(self->properties, pos);

}


inline static int
_cpl_propertylist_insert(cpl_propertylist *self, const cxchar *where,
                         cxbool after, const cxchar *name, cpl_type type,
                         cxcptr value)
{

    cx_deque_iterator pos;
    cpl_property *property;


    /*
     * Find the position where value should be inserted.
     */

    pos = _cpl_propertylist_find(self, where);

    if (pos == cx_deque_end(self->properties)) {
        return 1;
    }

    if (after) {
        pos = cx_deque_next(self->properties, pos);
    }


    /*
     * Create the property for value and fill it.
     */

    property = cpl_property_new(name, type);
    if (!property) {
        return 1;
    }


    /*
     * Map property type to the driver function's argument type.
     */

    switch (type) {
        case CPL_TYPE_CHAR:
            cpl_property_set_char(property, *((const cxchar *)value));
            break;

        case CPL_TYPE_BOOL:
            cpl_property_set_bool(property, *((const cxint *)value));
            break;

        case CPL_TYPE_INT:
            cpl_property_set_int(property, *((const cxint *)value));
            break;

        case CPL_TYPE_LONG:
            cpl_property_set_long(property, *((const cxlong *)value));
            break;

        case CPL_TYPE_LONG_LONG:
            cpl_property_set_long_long(property, *((const cxllong *)value));
            break;

        case CPL_TYPE_FLOAT:
            cpl_property_set_float(property, *((const cxfloat *)value));
            break;

        case CPL_TYPE_DOUBLE:
            cpl_property_set_double(property, *((const cxdouble *)value));
            break;

        case CPL_TYPE_STRING:
            cpl_property_set_string(property, ((const cxchar *)value));
            break;

        case CPL_TYPE_FLOAT_COMPLEX:
            cpl_property_set_float_complex(property,
                                           *((const float complex *)value));
            break;

        case CPL_TYPE_DOUBLE_COMPLEX:
            cpl_property_set_double_complex(property,
                                            *((const double complex *)value));
            break;

        default:
            return 1;
            break;
    }


    /*
     * Insert it into the deque
     */

    cx_deque_insert(self->properties, pos, property);

    return 0;

}


/*
 * @brief Insert cards from a FITS CHU into a propertylist
 * @param self   The propertylist to insert
 * @param file   The CFITSIO file object
 * @param hdumov Absolute extension number to move to first (0 for primary)
 * @param filter An optional compare function for filtering the properties
 * @param data   An optional regexp w. invert flag for filtering the properties
 * @return CPL_ERROR_NONE, or the relevant CPL error on failure
 * @see cpl_propertylist_to_fitsfile
 *
 * The function converts the current FITS header referenced by the cfitsio
 * file pointer to a property list. If the header cannot be accessed, i.e.
 * the number of keywords cannot be determined, or a keyword cannot be read,
 * -1 is returned. If the header contains FITS keywords whose type cannot be
 * determined the function returns 2. If a keyword type is not supported the
 * return value is 3.
 */

static cpl_error_code
_cpl_propertylist_fill_from_fits(cpl_propertylist *self, fitsfile *file,
                                int hdumov,
                                cx_compare_func filter, cxptr data)
{

    int status = 0;
    int ncards = 0;
    int i;


    if (hdumov >= 0 && fits_movabs_hdu(file, 1+hdumov, NULL, &status)) {
        return cpl_error_set_fits(CPL_ERROR_DATA_NOT_FOUND, status,
                                 fits_movabs_hdu, "hdumov=%d", hdumov);
    }

    if (fits_get_hdrspace(file, &ncards, NULL, &status)) {
        return cpl_error_set_fits(CPL_ERROR_FILE_IO, status,
                                  fits_get_hdrspace, "hdumov=%d", hdumov);
    }


    for (i = 1; i <= ncards; i++) {

        /*
         * In CFITSIO the length of a FITS comment card's value is larger
         * than what fits into a FITS value string. The FITS comment
         * card's 'value' is returned as comment by CFITSIO, but needs
         * to be copied to the value to build a correct property.
         * Therefore the value character array must have at least
         * the size FLEN_COMMENT.
         */

        char type                  = '\0';
        char key[FLEN_KEYWORD]     = "";
        char value[FLEN_COMMENT]   = "";
        char comment[FLEN_COMMENT] = "";

        cpl_property *property;


        if (fits_read_keyn(file, i, key, value, comment, &status)) {
            return cpl_error_set_fits(CPL_ERROR_FILE_IO, status,
                                      fits_read_keyn, "Card %d of %d",
                                      i, ncards);
        }

        if (filter != NULL && filter(key, data) == FALSE) {
            /* Card is filtered out */
            continue;
        }


        if (fits_get_keytype(value, &type, &status)) {

            if (status != VALUE_UNDEFINED) {
                return cpl_error_set_fits(CPL_ERROR_FILE_IO, status,
                                          fits_get_keytype, "card %d of %d "
                                          "has key='%s', value='%s', "
                                          "comment='%s'", i, ncards,
                                          key, value, comment);
            }


            /*
             *  Card has no value
             */

            status = 0;

            // FIXME: Check the following statement

            /*
             * Skip totally empty records.
             */

            if (key[0] == '\0' && comment[0] == '\0') {
                continue;
            }


            /*
             *  Make the comment of a FITS comment card the value of
             *  a string property.
             */

            strncpy(value, comment, FLEN_COMMENT);
            value[FLEN_COMMENT - 1] = '\0';

            comment[0] = '\0';
            type = 'C';

        }


        /*
         * Promote WCS keywords to double properties if the determined
         * type is integer.
         */

        if (type == 'I') {

            cxuint nwcs = CX_N_ELEMENTS(FITS_KEYMAP_WCS);
            cxuint j = _cpl_wcskeys_find(key, FITS_KEYMAP_WCS, nwcs);

            if ((j < nwcs) && (type != FITS_KEYMAP_WCS[j].type)) {
                type = FITS_KEYMAP_WCS[j].type;
            }

        }


        /*
         * Create the property from the parsed FITS card.
         */

        switch (type) {
        case 'L':
            {
                property = cpl_property_new(key, CPL_TYPE_BOOL);
                cpl_property_set_bool(property, *value == 'T' ? 1 : 0);

                break;
            }

        case 'I':
            {

                int itype = 0;
                int isgn  = 0;


                fits_get_inttype(value, &itype, &isgn, &status);

                if ((itype == TUINT) || (itype == TLONGLONG)) {

                    long long lval = 0LL;

                    const int nscan = sscanf(value, "%lld", &lval);


                    if (nscan != 1) {

                        return cpl_error_set_message_(CPL_ERROR_BAD_FILE_FORMAT,
                                                      "card %d of %d has "
                                                      "key='%s', "
                                                      "integer value='%s', "
                                                      "comment='%s'. "
                                                      "sscanf()=%d",
                                                      i, ncards,
                                                      key, value, comment,
                                                      nscan);

                    }

                    property = cpl_property_new(key, CPL_TYPE_LONG_LONG);
                    cpl_property_set_long_long(property, lval);


                }
                else {

                    int ival = 0;
                    const int nscan = sscanf(value, "%d", &ival);

                    if (nscan != 1) {

                        return cpl_error_set_message_(CPL_ERROR_BAD_FILE_FORMAT,
                                                      "card %d of %d has "
                                                      "key='%s', "
                                                      "integer value='%s', "
                                                      "comment='%s'. "
                                                      "sscanf()=%d",
                                                      i, ncards,
                                                      key, value, comment,
                                                      nscan);

                    }

                    property = cpl_property_new(key, CPL_TYPE_INT);
                    cpl_property_set_int(property, ival);

                }

                break;

            }

        case 'F':
            {

                double dval;

                const int nscan = sscanf(value, "%lf", &dval);


                if (nscan != 1) {

                    return cpl_error_set_message_(CPL_ERROR_BAD_FILE_FORMAT,
                                                  "card %d of %d has key='%s', "
                                                  "floating-point value='%s', "
                                                  "comment='%s'. sscanf()=%d",
                                                  i, ncards, key, value,
                                                  comment, nscan);
                }

                property = cpl_property_new(key, CPL_TYPE_DOUBLE);
                cpl_property_set_double(property, dval);

                break;

            }

        case 'C':

            /*
             * FITS standard: blank keywords may be followed by
             * any ASCII text as it is for COMMENT and HISTORY.
             *
             * In order to preserve this header record it is
             * changed into COMMENT record, so that it can be
             * stored in the property list.
             */

            property = cpl_property_new(key[0] != '\0' ? key : "COMMENT",
                                        CPL_TYPE_STRING);

            /*
             * Remove leading and trailing single quotes (string value
             * indicators) only it the keyword is not a comment or a history
             * keyword or a blank keyword.
             */

            if ((strncmp(key, "COMMENT", 7) != 0) &&
                    (strncmp(key, "HISTORY", 7) != 0)) {

                size_t sz = strlen(value);

                if (sz > 0) {

                    if (value[0] == '\'') {
                        value[0] = ' ';
                    }

                    if (value[sz - 1] == '\'') {
                        value[sz - 1] = ' ';
                    }

                    cx_strstrip(value);

                }

            }

            cpl_property_set_string(property, value);

            break;

        case 'X':
            {

                double dreal = 0.0;
                double dimag = 0.0;

                // FIXME: Find a CFITSIO routine to parse the string value
                const int nscan = sscanf(value, "(%lf,%lf", &dreal, &dimag);
                const double complex zval = dreal + dimag * _Complex_I;

                if (nscan != 2) {

                    return cpl_error_set_message_(CPL_ERROR_BAD_FILE_FORMAT,
                                                  "card %d of %d has key='%s', "
                                                  "complex floating-point "
                                                  "value='%s', comment='%s'. "
                                                  "sscanf()=%d", i, ncards, key,
                                                  value, comment, nscan);

                }

                property = cpl_property_new(key, CPL_TYPE_DOUBLE_COMPLEX);
                cpl_property_set_double_complex(property, zval);

                break;

            }

        default:
            return cpl_error_set_message_(CPL_ERROR_UNSUPPORTED_MODE,
                                          "card %d of %d has key='%s', "
                                          "value='%s', type='%c', comment='%s'",
                                          i, ncards, key, value, type,
                                          comment);
        }

        cpl_property_set_comment(property, comment);

        cx_deque_push_back(self->properties, property);

    }

    return 0;

}


/*
 * @brief
 *   Save the current property list to a FITS file using cfitsio.
 *
 * @param file   The CFITSIO file object
 * @param self   The propertylist to be saved
 * @param filter An optional compare function for filtering the properties
 * @param data   An optional regexp w. invert flag for filtering the properties
 *
 * @return
 *  CPL_ERROR_NONE, or the relevant CPL error on failure
 *
 * @note
 *  It will not save keys that match a regular expression given in filter.
 *
 * @see cpl_propertylist_to_fitsfile
 */

static cpl_error_code
_cpl_propertylist_to_fitsfile(fitsfile *file, const cpl_propertylist *self,
                              cx_compare_func filter, cxptr data)
{

    cx_deque_iterator first, last;

    const cpl_boolean do_filter = filter != NULL && data != NULL;

#ifdef CPL_USE_LONGSTRN
    /*
      Used for the OGIP Long String Keyword Convention:
      http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/ofwg_recomm/r13.html
    */
    cpl_boolean has_longstr = CPL_FALSE;
#endif

    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    first = cx_deque_begin(self->properties);
    last  = cx_deque_end(self->properties);

    while (first != last) {

        const cpl_property *p = cpl_propertylist_get_const(self, first);

        const cxchar *name    = cpl_property_get_name(p);
        const cxchar *comment = cpl_property_get_comment(p);

        cxint error = 0;

        const cpl_type type = cpl_property_get_type(p);


        if (do_filter && filter(name, data) == TRUE) {
            first = cx_deque_next(self->properties, first);
            continue;
        }

        switch (type) {
        case CPL_TYPE_CHAR:
            {
                const cxchar c = cpl_property_get_char(p);
                cxchar value[2];

                /*
                 * Character properties should be represented as a single
                 * character string, not as its numerical equivalent.
                 */

                value[0] = c;
                value[1] = '\0';

                (void)fits_update_key(file, TSTRING, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        case CPL_TYPE_BOOL:
            {
                const cxint b = cpl_property_get_bool(p);
                const cxint value = b == TRUE ? 1 : 0;

                (void)fits_update_key(file, TLOGICAL, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        case CPL_TYPE_INT:
            {
                const cxint value = cpl_property_get_int(p);

                (void)fits_update_key(file, TINT, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        case CPL_TYPE_LONG:
            {
                const cxlong value = cpl_property_get_long(p);

                (void)fits_update_key(file, TLONG, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        case CPL_TYPE_LONG_LONG:
            {
                const long long value = cpl_property_get_long_long(p);

                (void)fits_update_key(file, TLONGLONG, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        case CPL_TYPE_FLOAT:
            {
                const cxfloat value = cpl_property_get_float(p);

                (void)fits_update_key(file, TFLOAT, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        case CPL_TYPE_DOUBLE:
            {
                const cxdouble value = cpl_property_get_double(p);

                (void)fits_update_key(file, TDOUBLE, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        case CPL_TYPE_STRING:
            {
                const cxchar *value = cpl_property_get_string(p);

                if (strcmp(name, "COMMENT") == 0) {

                    if (fits_write_comment(file, (strlen(value) < 1
                                                  ? " " : value),
                                           &error)) {
                        return cpl_error_set_fits(CPL_ERROR_ILLEGAL_INPUT,
                                                  error, fits_write_comment,
                                                  "name='%s', type=%d ('%s'), "
                                                  "comment='%s'", name,
                                                  type, cpl_type_get_name(type),
                                                  comment);
                    }

                }
                else if (strcmp(name, "HISTORY") == 0) {

                    if (fits_write_history(file, (strlen(value) < 1
                                                  ? " " : value),
                                           &error)) {
                        return cpl_error_set_fits(CPL_ERROR_ILLEGAL_INPUT,
                                                  error, fits_write_history,
                                                  "name='%s', type=%d ('%s'), "
                                                  "comment='%s'", name,
                                                  type, cpl_type_get_name(type),
                                                  comment);
                    }

                }
                else {

#ifdef CPL_USE_LONGSTRN
                    /* In v. 3.24 fits_update_key_longstr() has a buffer
                       overflow triggered by a long comment. :-(((((((((((( */
                    /* Try to work around this by using a truncated comment. */
                    cxchar shortcomment[FLEN_COMMENT];

                    if (comment != NULL) {
                        (void)strncpy(shortcomment, comment, FLEN_COMMENT-1);
                        shortcomment[FLEN_COMMENT-1] = '\0';
                    }

                    if (!has_longstr) {

                        has_longstr = CPL_TRUE;

                        if (fits_write_key_longwarn(file, &error)) {
                            return cpl_error_set_fits(CPL_ERROR_FILE_IO, error,
                                                      fits_write_key_longwarn,
                                                      " ");
                        }

                    }

                    if (fits_update_key_longstr(file, (cxchar*)name,
                                                (cxptr)value,
                                                comment != NULL ?
                                                (cxchar*)shortcomment : NULL,
                                                &error)) {
                        return cpl_error_set_fits(CPL_ERROR_ILLEGAL_INPUT,
                                                  error,fits_update_key_longstr,
                                                  "name='%s', value='%s', "
                                                  "comment='%s'",
                                                  name, value, comment);
                    }

#else

                    if (fits_update_key_str(file, name, (cxptr)value,
                                            (cxchar*)comment, &error)) {
                        return cpl_error_set_fits(CPL_ERROR_ILLEGAL_INPUT,
                                                  error, fits_update_key_str,
                                                  "name='%s', value='%s', "
                                                  "comment='%s'",
                                                  name, value, comment);
                    }

#endif

                }

            }
            break;

        case CPL_TYPE_FLOAT_COMPLEX:
            {
                const float complex value = cpl_property_get_float_complex(p);

                (void)fits_update_key(file, TCOMPLEX, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        case CPL_TYPE_DOUBLE_COMPLEX:
            {
                const double complex value = cpl_property_get_double_complex(p);

                (void)fits_update_key(file, TDBLCOMPLEX, name,
                                      (cxptr)&value, (cxchar*)comment, &error);
            }
            break;

        default:
            return cpl_error_set_message_(CPL_ERROR_UNSUPPORTED_MODE,
                                         "name='%s', type=%d ('%s'), "
                                         "comment='%s'", name, type,
                                         cpl_type_get_name(type), comment);

        }

        if (error) {
            return cpl_error_set_fits(CPL_ERROR_ILLEGAL_INPUT, error,
                                      fits_update_key, "name='%s', "
                                      "type=%d ('%s'), comment='%s'", name,
                                      type, cpl_type_get_name(type), comment);
        }

        first = cx_deque_next(self->properties, first);

    }

    return CPL_ERROR_NONE;

}


/*
 * Public methods
 */

/**
 * @brief
 *    Create an empty property list.
 *
 * @return
 *   The newly created property list.
 *
 * The function creates a new property list and returns a handle for it.
 * To destroy the returned property list object use the property list
 * destructor @b cpl_propertylist_delete().
 *
 * @see cpl_propertylist_delete()
 */

cpl_propertylist *
cpl_propertylist_new(void)
{

    cpl_propertylist *self = cx_malloc(sizeof *self);

    self->properties = cx_deque_new();

    return self;

}


/**
 * @brief
 *   Create a copy of the given property list.
 *
 * @param self  The property list to be copied.
 *
 * @return
 *   The created copy or @c NULL in case an error occurred.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a deep copy of the given property list @em self,
 * i.e the created copy and the original property list do not share any
 * resources.
 */

cpl_propertylist *
cpl_propertylist_duplicate(const cpl_propertylist *self)
{


    cx_deque_iterator first, last;

    cpl_propertylist *copy = NULL;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    cx_assert(self->properties != NULL);


    copy = cpl_propertylist_new();

    first = cx_deque_begin(self->properties);
    last = cx_deque_end(self->properties);

    while (first != last) {
        cpl_property *tmp = cx_deque_get(self->properties, first);

        cx_deque_push_back(copy->properties, cpl_property_duplicate(tmp));
        first = cx_deque_next(self->properties, first);
    }

    return copy;

}


/**
 * @brief
 *    Destroy a property list.
 *
 * @param self  The property list to .
 *
 * @return
 *   Nothing.
 *
 * The function destroys the property list @em self and its whole
 * contents. If @em self is @c NULL, nothing is done and no error is set.
 */

void
cpl_propertylist_delete(cpl_propertylist *self)
{

    if (self) {
        cx_deque_destroy(self->properties, (cx_free_func)cpl_property_delete);
        cx_free(self);
    }

    return;

}


/**
 * @brief
 *    Get the current size of a property list.
 *
 * @param self  A property list.
 *
 * @return
 *   The property list's current size, or 0 if the list is empty. If an
 *   error occurs the function returns 0 and sets an appropriate error
 *   code.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function reports the current number of elements stored in the property
 * list @em self.
 */

cpl_size
cpl_propertylist_get_size(const cpl_propertylist *self)
{



    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0L;
    }

    return (cpl_size) cx_deque_size(self->properties);

}


/**
 * @brief
 *   Check whether a property list is empty.
 *
 * @param self  A property list.
 *
 * @return
 *   The function returns 1 if the list is empty, and 0 otherwise.
 *   In case an error occurs the function returns -1 and sets an
 *   appropriate error code.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function checks if @em self contains any properties.
 */

int
cpl_propertylist_is_empty(const cpl_propertylist *self)
{



    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return -1;
    }

    return cx_deque_empty(self->properties);

}


/**
 * @brief
 *   Get the the type of a property list entry.
 *
 * @param self  A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The type of the stored value. If an error occurs the function returns
 *   @c CPL_TYPE_INVALID and sets an appropriate error code.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns the type of the value stored in @em self with the
 * name @em name. If there is more than one property with the same @em name,
 * it takes the first one from the list.
 */

cpl_type
cpl_propertylist_get_type(const cpl_propertylist *self, const char *name)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return CPL_TYPE_INVALID;
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return CPL_TYPE_INVALID;
    }

    return cpl_property_get_type(property);

}


/**
 * @brief
 *   Check whether a property is present in a property list.
 *
 * @param self  A property list.
 * @param name  The property name to look up.
 *
 * @return
 *   The function returns 1 if the property is present, or 0 otherwise.
 *   If an error occurs the function returns 0 and sets an appropriate
 *   error code.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property with
 * the name @em name and reports whether it was found or not.
 */

int
cpl_propertylist_has(const cpl_propertylist *self, const char *name)
{



    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    return _cpl_propertylist_get(self, name) != NULL ? 1 : 0;

}


/**
 * @brief
 *   Modify the comment field of the given property list entry.
 *
 * @param self    A property list.
 * @param name     The property name to look up.
 * @param comment  New comment string.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its comment is replaced by
 * the string @em comment. The provided comment string may be @c NULL.
 * In this case an already existing comment is deleted. If there is more
 * than one property with the same @em name, it takes the first one from the
 * list.
 */

cpl_error_code
cpl_propertylist_set_comment(cpl_propertylist *self, const char *name,
                             const char *comment)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    cpl_property_set_comment(property, comment);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of the given character property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New character value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its character value is replaced
 * with the character @em value. If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_char(cpl_propertylist *self, const char *name,
                          char value)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_char(property, value);

}


/**
 * @brief
 *   Set the value of the given boolean property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New boolean value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its boolean value is replaced
 * with the boolean @em value.  If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_bool(cpl_propertylist *self, const char *name, int value)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_bool(property, value);

}


/**
 * @brief
 *   Set the value of the given integer property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New integer value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its integer value is replaced
 * with the integer @em value. If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_int(cpl_propertylist *self, const char *name, int value)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_int(property, value);

}


/**
 * @brief
 *   Set the value of the given long property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New long value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its long value is replaced with
 * the long @em value. If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_long(cpl_propertylist *self, const char *name,
                          long value)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_long(property, value);

}


/**
 * @brief
 *   Set the value of the given long long property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New long long value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its long long value is replaced
 * with @em value. If there is more than one property with the same @em name,
 * it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_long_long(cpl_propertylist *self, const char *name,
                               long long value)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_long_long(property, value);

}


/**
 * @brief
 *   Set the value of the given float property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New float value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its float value is replaced with
 * the float @em value. If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_float(cpl_propertylist *self, const char *name,
                           float value)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_float(property, value);

}


/**
 * @brief
 *   Set the value of the given double property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New double value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its double value is replaced with
 * the double @em value. If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_double(cpl_propertylist *self, const char *name,
                            double value)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_double(property, value);

}


/**
 * @brief
 *   Set the value of the given string property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New string value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>name</i> or <i>value</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its string value is replaced with
 * the string @em value. If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_string(cpl_propertylist *self, const char *name,
                            const char *value)
{


    cpl_property *property;


    if (self == NULL || name == NULL || value == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_string(property, value);

}


/**
 * @brief
 *   Set the value of the given float complex property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New float complex value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its float complex value is replaced
 * with the float complex @em value. If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_float_complex(cpl_propertylist *self, const char *name,
                                   float complex value)
{
    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_float_complex(property, value)
        ? cpl_error_set_where_() : CPL_ERROR_NONE;

}


/**
 * @brief
 *   Set the value of the given double complex property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 * @param value  New double complex value.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its double complex value is replaced
 * with the double complex @em value. If there is more than one property with
 * the same @em name, it takes the first one from the list.
 */

cpl_error_code
cpl_propertylist_set_double_complex(cpl_propertylist *self, const char *name,
                                    double complex value)
{

    cpl_property *property;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = _cpl_propertylist_get(self, name);

    if (property == NULL) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    return cpl_property_set_double_complex(property, value)
        ? cpl_error_set_where_() : CPL_ERROR_NONE; 

}


/**
 * @brief
 *   Access property list elements by index.
 *
 * @param self      The property list to query.
 * @param position  Index of the element to retrieve.
 *
 * @return
 *   The function returns the property with index @em position, or @c NULL
 *   if @em position is out of range.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns a handle for the property list element, the property,
 * with the index @em position. Numbering of property list elements extends from
 * 0 to @b cpl_propertylist_get_size() - 1. If @em position is less than 0 or
 * greater equal than @b cpl_propertylist_get_size() the function returns
 * @c NULL.
 */

const cpl_property *
cpl_propertylist_get_const(const cpl_propertylist *self, long position)
{


#ifdef CPL_PROPERTYLIST_ENABLE_LOOP_CHECK
    cxsize i = 0;
#endif

    cx_deque_iterator first, last;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (position < 0) {
        return NULL;
    }

    first = cx_deque_begin(self->properties);
    last = cx_deque_end(self->properties);

#ifdef CPL_PROPERTYLIST_ENABLE_LOOP_CHECK
    /* FIXME: Except from never stopping in a unit test, what does this do? */
    while (i < (cxsize)position && first != last) {
        first = cx_deque_next(self->properties, first);
        i++;
    }
#endif

    if (first == last) {
        return NULL;
    }

    if ((cx_deque_const_iterator)position >= last) {
        return NULL;
    }

    return cx_deque_get(self->properties, (cx_deque_const_iterator)position);

}

/**
 * @brief
 *   Access property list elements by index.
 *
 * @param self      The property list to query.
 * @param position  Index of the element to retrieve.
 *
 * @return
 *   The function returns the property with index @em position, or @c NULL
 *   if @em position is out of range.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns a handle for the property list element, the property,
 * with the index @em position. Numbering of property list elements extends from
 * 0 to @b cpl_propertylist_get_size() - 1. If @em position is less than 0 or
 * greater equal than @b cpl_propertylist_get_size() the function returns
 * @c NULL.
 */

cpl_property *
cpl_propertylist_get(cpl_propertylist *self, long position)
{

    cpl_errorstate prestate = cpl_errorstate_get();
    cpl_property *property = (cpl_property *)cpl_propertylist_get_const(self, position);

    if (!cpl_errorstate_is_equal(prestate))
        (void)cpl_error_set_where_();

    return property;

}
/**
 * @brief
 *   Get the comment of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The comment of the property list entry, or @c NULL.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its comment string is returned.
 * If an entry with the name @em name is not found, or if the entry has
 * no comment the function returns @c NULL. If there is more than one property
 * with the same @em name, it takes the first one from the list.
 */

const char *
cpl_propertylist_get_comment(const cpl_propertylist *self, const char *name)
{


    cpl_property *property;


    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return NULL;
    }

    return cpl_property_get_comment(property);

}


/**
 * @brief
 *   Get the character value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The character value stored in the list entry. The function returns '\\0'
 *   if an error occurs and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         @c CPL_TYPE_CHAR.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its character value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list.
 */

char
cpl_propertylist_get_char(const cpl_propertylist *self, const char *name)
{


    cxchar result;

    cpl_property *property;

    cpl_errorstate prevstate;

    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return '\0';
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return '\0';
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_char(property);

    /*
     * If an error occurred change any possible set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        cpl_error_set_where_();
        return '\0';
    }


    return result;

}


/**
 * @brief
 *   Get the boolean value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The integer representation of the boolean value stored in
 *   the list entry. @c TRUE is represented as non-zero value while
 *   0 indicates @c FALSE. The function returns 0 if an error occurs
 *   and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         @c CPL_TYPE_BOOL.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its boolean value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list.
 */

int
cpl_propertylist_get_bool(const cpl_propertylist *self, const char *name)
{


    cxbool result;

    cpl_property *property;

    cpl_errorstate prevstate;

    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return 0;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_bool(property);

    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        cpl_error_set_where_();
        return 0;
    }


    return result == TRUE ? 1 : 0;

}


/**
 * @brief
 *   Get the integer value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The integer value stored in the list entry. The function returns 0 if
 *   an error occurs and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         @c CPL_TYPE_INT.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its integer value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list.
 */

int
cpl_propertylist_get_int(const cpl_propertylist *self, const char *name)
{


    cxint result;

    cpl_property *property;

    cpl_errorstate prevstate;

    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return 0;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_int(property);

    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        cpl_error_set_where_();
        return 0;
    }

    return result;

}


/**
 * @brief
 *   Get the long value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The long value stored in the list entry. The function returns 0 if
 *   an error occurs and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         <tt>CPL_TYPE_LONG</tt> or <tt>CPL_TYPE_INT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its long value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list.
 *
 * The function may be used to access the value of all integer type properties
 * whose integer value has a rank less or equal to the functions return type.
 * If the value of a compatible property is retrieved, it is promoted to
 * the return type of the function.
 */

long
cpl_propertylist_get_long(const cpl_propertylist *self, const char *name)
{


    cxlong result;

    cpl_property *property;

    cpl_errorstate prevstate;


    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return 0;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_long(property);

    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        cpl_error_set_where_();
        return 0;
    }


    return result;

}


/**
 * @brief
 *   Get the long long value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The long value stored in the list entry. The function returns 0 if
 *   an error occurs and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         <tt>CPL_TYPE_LONG_LONG</tt>, <tt>CPL_TYPE_LONG</tt>, or
 *         <tt>CPL_TYPE_INT</tt>
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its long long value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list.
 *
 * The function may be used to access the value of all integer type properties
 * whose integer value has a rank less or equal to the functions return type.
 * If the value of a compatible property is retrieved, it is promoted to
 * the return type of the function.
 */

long long
cpl_propertylist_get_long_long(const cpl_propertylist *self, const char *name)
{


    cxlong result;

    cpl_property *property;

    cpl_errorstate prevstate;


    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return 0;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_long_long(property);


    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        cpl_error_set_where_();
        return 0;
    }


    return result;

}


/**
 * @brief
 *   Get the float value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The float value stored in the list entry. The function returns 0 if
 *   an error occurs and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         @c CPL_TYPE_FLOAT or @c CPL_TYPE_DOUBLE.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its float value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list. If the value is of type double, the function
 * casts it to float before returning it.
 */

float
cpl_propertylist_get_float(const cpl_propertylist *self, const char *name)
{


    cxfloat result;

    cpl_property *property;

    cpl_errorstate prevstate;

    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return 0;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_float(property);

    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        cpl_error_set_where_();
        return 0;
    }

    return result;

}


/**
 * @brief
 *   Get the double value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The double value stored in the list entry. The function returns 0 if
 *   an error occurs and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         <tt>CPL_TYPE_DOUBLE</tt> or <tt>CPL_TYPE_FLOAT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its double value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list. If the value is of type float, the function
 * casts it to double before returning it.
 *
 * The function may be used to access the value of all floating-point type
 * properties whose floating-point value has a rank less or equal to the
 * functions return type. If the value of a compatible property is retrieved,
 * it is promoted to the return type of the function.
 */

double
cpl_propertylist_get_double(const cpl_propertylist *self, const char *name)
{


    cxdouble result;

    cpl_property *property;

    cpl_errorstate prevstate;


    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return 0;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_double(property);

    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        cpl_error_set_where_();
        return 0;
    }


    return result;

}


/**
 * @brief
 *   Get the string value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   A handle to the string value stored in the list entry. The
 *   function returns @c NULL if an error occurs and an appropriate
 *   error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         @c CPL_TYPE_STRING.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, a handle to its string value
 * is returned. If there is more than one property with the same @em name,
 * it takes the first one from the list.
 */

const char *
cpl_propertylist_get_string(const cpl_propertylist *self, const char *name)
{


    const cxchar *result;

    cpl_property *property;

    cpl_errorstate prevstate;

    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return NULL;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_string(property);

    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        cpl_error_set_where_();
        return NULL;
    }

    return result;

}


/**
 * @brief
 *   Get the float complex value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The float complex value stored in the list entry. The function returns 0 if
 *   an error occurs and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         @c CPL_TYPE_FLOAT_COMPLEX or @c CPL_TYPE_DOUBLE_COMPLEX.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its float complex value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list. If the value is of type double, the function
 * casts it to float complex before returning it.
 */

float complex
cpl_propertylist_get_float_complex(const cpl_propertylist *self,
                                   const char *name)
{

    float complex result;

    cpl_property *property;

    cpl_errorstate prevstate;

    if (self == NULL || name == NULL) {
        (void)cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0.0;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        (void)cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return 0.0;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_float_complex(property);

    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        (void)cpl_error_set_where_();
        return 0.0;
    }

    return result;

}


/**
 * @brief
 *   Get the double complex value of the given property list entry.
 *
 * @param self   A property list.
 * @param name   The property name to look up.
 *
 * @return
 *   The double complex value stored in the list entry. The function returns 0 if
 *   an error occurs and an appropriate error code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>self</i> does not contain a property with
 *         the name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The sought-after property <i>name</i> is not of type
 *         @c CPL_TYPE_DOUBLE_COMPLEX or @c CPL_TYPE_FLOAT_COMPLEX.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property list @em self for a property named
 * @em name. If it is present in the list, its double complex value is returned.
 * If there is more than one property with the same @em name, it takes the
 * first one from the list. If the value is of type float, the function
 * casts it to double complex before returning it.
 *
 * The function may be used to access the value of all complex type
 * properties whose complex value has a rank less or equal to the
 * functions return type. If the value of a compatible property is retrieved,
 * it is promoted to the return type of the function.
 */

double complex
cpl_propertylist_get_double_complex(const cpl_propertylist *self,
                                    const char *name)
{

    double complex result;

    cpl_property *property;

    cpl_errorstate prevstate;


    if (self == NULL || name == NULL) {
        (void)cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0.0;
    }

    property = _cpl_propertylist_get(self, name);

    if (!property) {
        (void)cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
        return 0.0;
    }

    prevstate = cpl_errorstate_get();

    result = cpl_property_get_double_complex(property);

    /*
     * If an error occurred change any possibly set location to this
     * function.
     */

    if (!cpl_errorstate_is_equal(prevstate)) {
        (void)cpl_error_set_where_();
        return 0.0;
    }


    return result;

}

/**
 * @brief
 *   Insert a character value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The character value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new character property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_char(cpl_propertylist *self, const char *here,
                             const char *name, char value)
{


    cxint status = 0;


    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, here, FALSE, name,
                                      CPL_TYPE_CHAR, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a boolean value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The boolean value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new boolean property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_bool(cpl_propertylist *self, const char *here,
                             const char *name, int value)
{


    cxint status = 0;


    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status =  _cpl_propertylist_insert(self, here, FALSE, name,
                                       CPL_TYPE_BOOL, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a integer value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The integer value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new integer property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_int(cpl_propertylist *self, const char *here,
                            const char *name, int value)
{


    cxint status = 0;


    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, here, FALSE, name,
                                      CPL_TYPE_INT, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a long value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new long property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_long(cpl_propertylist *self, const char *here,
                             const char *name, long value)
{


    cxint status = 0;


    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, here, FALSE, name,
                                      CPL_TYPE_LONG, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a long long value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The long long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new long long property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_long_long(cpl_propertylist *self, const char *here,
                                  const char *name, long long value)
{


    cxint status = 0;


    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, here, FALSE, name,
                                      CPL_TYPE_LONG_LONG, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a float value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The float value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new float property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_float(cpl_propertylist *self, const char *here,
                              const char *name, float value)
{


    cxint status = 0;


    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, here, FALSE, name,
                                      CPL_TYPE_FLOAT, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a double value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The double value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new double property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_double(cpl_propertylist *self, const char *here,
                               const char *name, double value)
{


    cxint status = 0;


    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, here, FALSE, name,
                                      CPL_TYPE_DOUBLE, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a string value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The string value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i>, <i>name</i> or <i>value</i>
 *         or <i>name</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new string property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_string(cpl_propertylist *self, const char *here,
                               const char *name, const char *value)
{


    cxint status = 0;


    if (self == NULL || here == NULL || name == NULL || value == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, here, FALSE, name,
                                      CPL_TYPE_STRING, (cxcptr)value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a float complex value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The float complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new float complex property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_float_complex(cpl_propertylist *self, const char *here,
                                      const char *name, float complex value)
{

    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    return _cpl_propertylist_insert(self, here, FALSE, name,
                                    CPL_TYPE_FLOAT_COMPLEX,
                                    (cxcptr)&value)
        ? cpl_error_set_(CPL_ERROR_UNSPECIFIED) : CPL_ERROR_NONE;
}


/**
 * @brief
 *   Insert a double complex value into a property list at the given position.
 *
 * @param self   A property list.
 * @param here   Name indicating the position at which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The double complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>here</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new double complex property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_double_complex(cpl_propertylist *self, const char *here,
                                       const char *name, double complex value)
{

    if (self == NULL || here == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    return _cpl_propertylist_insert(self, here, FALSE, name,
                                    CPL_TYPE_DOUBLE_COMPLEX,
                                    (cxcptr)&value)
        ? cpl_error_set_(CPL_ERROR_UNSPECIFIED) : CPL_ERROR_NONE;
}

/**
 * @brief
 *   Insert a character value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The character value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new character property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_char(cpl_propertylist *self, const char *after,
                                   const char *name, char value)
{


    cxint status = 0;


    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, after, TRUE, name,
                                      CPL_TYPE_CHAR, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a boolean value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The boolean value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new boolean property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_bool(cpl_propertylist *self, const char *after,
                                   const char *name, int value)
{


    cxint status = 0;


    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, after, TRUE, name,
                                      CPL_TYPE_BOOL, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a integer value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The integer value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new integer property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_int(cpl_propertylist *self, const char *after,
                                  const char *name, int value)
{


    cxint status = 0;


    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, after, TRUE, name,
                                      CPL_TYPE_INT, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a long value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new long property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_long(cpl_propertylist *self, const char *after,
                                   const char *name, long value)
{


    cxint status = 0;


    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, after, TRUE, name,
                                      CPL_TYPE_LONG, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a long long value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The long long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new long long property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_long_long(cpl_propertylist *self,
                                        const char *after, const char *name,
                                        long long value)
{


    cxint status = 0;


    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, after, TRUE, name,
                                      CPL_TYPE_LONG_LONG, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a float value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The float value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new float property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_float(cpl_propertylist *self, const char *after,
                                    const char *name, float value)
{


    cxint status = 0;


    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, after, TRUE, name,
                                      CPL_TYPE_FLOAT, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a double value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The double value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new double property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_double(cpl_propertylist *self,
                                     const char *after, const char *name,
                                     double value)
{


    cxint status = 0;


    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status = _cpl_propertylist_insert(self, after, TRUE, name,
                                      CPL_TYPE_DOUBLE, (cxcptr)&value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Insert a string value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The string value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i>, <i>name</i> or <i>value</i>
 *         or <i>name</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new string property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_string(cpl_propertylist *self,
                                     const char *after, const char *name,
                                     const char *value)
{


    cxint status = 0;


    if (self == NULL || after == NULL || name == NULL || value == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    status =  _cpl_propertylist_insert(self, after, TRUE, name,
                                       CPL_TYPE_STRING, (cxcptr)value);

    if (status) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    return CPL_ERROR_NONE;

}

/**
 * @brief
 *   Insert a float complex value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The float complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new float complex property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_float_complex(cpl_propertylist *self,
                                            const char *after, const char *name,
                                            float complex value)
{

    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    return _cpl_propertylist_insert(self, after, TRUE, name,
                                    CPL_TYPE_FLOAT_COMPLEX, (cxcptr)&value)
        ? cpl_error_set_(CPL_ERROR_UNSPECIFIED) : CPL_ERROR_NONE;
}


/**
 * @brief
 *   Insert a double complex value into a property list after the given position.
 *
 * @param self   A property list.
 * @param after  Name of the property after which the value is inserted.
 * @param name   The property name to be assigned to the value.
 * @param value  The double complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>after</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_UNSPECIFIED</td>
 *       <td class="ecr">
 *         A property with the name <i>name</i> could not be inserted
 *         into <i>self</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new double complex property with name @em name and
 * value @em value. The property is inserted into the property list
 * @em self after the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_double_complex(cpl_propertylist *self,
                                             const char *after,
                                             const char *name,
                                             double complex value)
{

    if (self == NULL || after == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    return _cpl_propertylist_insert(self, after, TRUE, name,
                                    CPL_TYPE_DOUBLE_COMPLEX,
                                    (cxcptr)&value)
        ? cpl_error_set_(CPL_ERROR_UNSPECIFIED) : CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a character value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The character value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new character property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_char(cpl_propertylist *self, const char *name,
                              char value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_CHAR);
    cx_assert(property != NULL);

    cpl_property_set_char(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a boolean value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The boolean value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new boolean property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_bool(cpl_propertylist *self, const char *name,
                              int value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_BOOL);
    cx_assert(property != NULL);

    cpl_property_set_bool(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a integer value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The integer value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new integer property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_int(cpl_propertylist *self, const char *name,
                             int value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_INT);
    cx_assert(property != NULL);

    cpl_property_set_int(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a long value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new long property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_long(cpl_propertylist *self, const char *name,
                              long value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_LONG);
    cx_assert(property != NULL);

    cpl_property_set_long(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a long long value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The long long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new long long property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_long_long(cpl_propertylist *self, const char *name,
                                   long long value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_LONG_LONG);
    cx_assert(property != NULL);

    cpl_property_set_long_long(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a float value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The float value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new float property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_float(cpl_propertylist *self, const char *name,
                               float value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_FLOAT);
    cx_assert(property != NULL);

    cpl_property_set_float(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a double value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The double value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new double property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_double(cpl_propertylist *self, const char *name,
                                double value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_DOUBLE);
    cx_assert(property != NULL);

    cpl_property_set_double(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a string value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The string value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>name</i> or <i>value</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new string property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_string(cpl_propertylist *self, const char *name,
                                const char *value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL || value == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_STRING);
    cx_assert(property != NULL);

    cpl_property_set_string(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a float complex value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The float complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new float complex property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_float_complex(cpl_propertylist *self, const char *name,
                                       float complex value)
{

    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_FLOAT_COMPLEX);
    cx_assert(property != NULL);

    cpl_property_set_float_complex(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Prepend a double complex value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The double complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new double complex property with name @em name and
 * value @em value. The property is prepended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_prepend_double_complex(cpl_propertylist *self, const char *name,
                                        double complex value)
{

    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_DOUBLE_COMPLEX);
    cx_assert(property != NULL);

    cpl_property_set_double_complex(property, value);
    cx_deque_push_front(self->properties, property);

    return CPL_ERROR_NONE;

}

/**
 * @brief
 *   Append a character value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The character value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new character property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_char(cpl_propertylist *self, const char *name,
                             char value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_CHAR);
    cx_assert(property != NULL);

    cpl_property_set_char(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append a boolean value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The boolean value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new boolean property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_bool(cpl_propertylist *self, const char *name,
                             int value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_BOOL);
    cx_assert(property != NULL);

    cpl_property_set_bool(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append an integer value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The integer value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new integer property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_int(cpl_propertylist *self, const char *name,
                            int value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_INT);
    cx_assert(property != NULL);

    cpl_property_set_int(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append a long value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new long property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_long(cpl_propertylist *self, const char *name,
                             long value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_LONG);
    cx_assert(property != NULL);

    cpl_property_set_long(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append a long long value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The long long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new long long property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_long_long(cpl_propertylist *self, const char *name,
                                  long long value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_LONG_LONG);
    cx_assert(property != NULL);

    cpl_property_set_long_long(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append a float value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The float value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new float property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_float(cpl_propertylist *self, const char *name,
                              float value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_FLOAT);
    cx_assert(property != NULL);

    cpl_property_set_float(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append a double value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The double value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new double property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_double(cpl_propertylist *self, const char *name,
                               double value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_DOUBLE);
    cx_assert(property != NULL);

    cpl_property_set_double(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append a string value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The string value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>name</i> or <i>value</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new string property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_string(cpl_propertylist *self, const char *name,
                               const char *value)
{


    cpl_property *property = NULL;


    if (self == NULL || name == NULL || value == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_STRING);
    cx_assert(property != NULL);

    cpl_property_set_string(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append a float complex value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The float complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new float complex property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_float_complex(cpl_propertylist *self, const char *name,
                                      float complex value)
{

    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_FLOAT_COMPLEX);
    cx_assert(property != NULL);

    cpl_property_set_float_complex(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Append a double complex value to a property list.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The double complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new double complex property with name @em name and
 * value @em value. The property is appended to the property list @em self.
 */

cpl_error_code
cpl_propertylist_append_double_complex(cpl_propertylist *self, const char *name,
                                       double complex value)
{

    cpl_property *property = NULL;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    property = cpl_property_new(name, CPL_TYPE_DOUBLE_COMPLEX);
    cx_assert(property != NULL);

    cpl_property_set_double_complex(property, value);
    cx_deque_push_back(self->properties, property);

    return CPL_ERROR_NONE;

}

/**
 * @brief
 *   Append a property list..
 *
 * @param self   A property list.
 * @param other  The property to append.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function appends the property list @em other to the property list
 * @em self.
 */

cpl_error_code
cpl_propertylist_append(cpl_propertylist *self,
                        const cpl_propertylist *other)
{


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    if (other != NULL) {

        cx_deque_const_iterator pos = cx_deque_begin(other->properties);

        while (pos != cx_deque_end(other->properties)) {

            const cpl_property *p = cx_deque_get(other->properties, pos);

            cx_deque_push_back(self->properties, cpl_property_duplicate(p));
            pos = cx_deque_next(other->properties, pos);

        }

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Erase the given property from a property list.
 *
 * @param self   A property list.
 * @param name   Name of the property to erase.
 *
 * @return
 *   On success the function returns the number of erased entries. If
 *   an error occurs the function returns 0 and an appropriate error
 *   code is set.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches the property with the name @em name in the property
 * list @em self and removes it. The property is destroyed. If @em self
 * contains multiple duplicates of a property named @em name, only the
 * first one is erased.
 */

int
cpl_propertylist_erase(cpl_propertylist *self, const char *name)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return 0;
    }

    pos = _cpl_propertylist_find(self, name);
    if (pos == cx_deque_end(self->properties)) {
        return 0;
    }

    cx_deque_erase(self->properties, pos, (cx_free_func)cpl_property_delete);

    return 1;

}

/**
 * @brief
 *   Erase all properties with name matching a given regular expression.
 *
 * @param self    A property list.
 * @param regexp  Regular expression.
 * @param invert  Flag inverting the sense of matching.
 *
 * @return
 *   On success the function returns the number of erased entries or 0 if
 *   no entries are erased. If an error occurs the function returns -1 and an
 *   appropriate error code is set. In CPL versions earlier than 5.0, the
 *   return value in case of error is 0.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>regexp</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_ILLEGAL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>regexp</i> is the empty string.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function searches for all the properties matching in the list
 * @em self and removes them. Whether a property matches or not depends on
 * the given regular expression @em regexp, and the flag @em invert. If
 * @em invert is @c 0, all properties matching @em regexp are removed from
 * the list. If @em invert is set to a non-zero value, all properties which
 * do not match @em regexp are erased. The removed properties are destroyed.
 *
 * The function expects POSIX 1003.2 compliant extended regular expressions.
 */

int
cpl_propertylist_erase_regexp(cpl_propertylist *self, const char *regexp,
                              int invert)
{


    cxint status = 0;
    cxint count = 0;

    cx_deque_iterator first, pos;

    cpl_property    *p;

    cpl_regexp filter;


    if (self == NULL || regexp == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return -1;
    }

    if (regexp[0] == '\0') {
        cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
        return -1;
    }

    status = regcomp(&filter.re, regexp, REG_EXTENDED | REG_NOSUB);
    if (status) {
        (void)cpl_error_set_regex(CPL_ERROR_ILLEGAL_INPUT, status,
                                  &filter.re, "regexp='%s', invert=%d",
                                  regexp, invert);
        return -1;
    }

    filter.invert = invert == 0 ? FALSE : TRUE;

    first = cx_deque_begin(self->properties);

    while (first < cx_deque_end(self->properties)) {
        pos = first;

        p = cx_deque_get(self->properties, pos);
        if (_cpl_propertylist_compare_regexp(p, &filter) == TRUE) {
            cx_deque_erase(self->properties, pos,
                          (cx_free_func)cpl_property_delete);
            count++;
        }
        else {
            first = cx_deque_next(self->properties, first);
        }
    }

    regfree(&filter.re);

    return count;

}


/**
 * @brief
 *   Remove all properties from a property list.
 *
 * @param self  A property list.
 *
 * @return Nothing.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function removes all properties from @em self. Each property
 * is properly deallocated. After calling this function @em self is
 * empty.
 */

void
cpl_propertylist_empty(cpl_propertylist *self)
{


    cx_deque_iterator first;


    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return;
    }

    first = cx_deque_begin(self->properties);

    /* cx_deque_end changes its value everytime cx_deque_erase()
    * is called. The elements are shifted from end to begin every
    * time an element is erased, therefore we always erase element
    * first with a value of 0.
    */
    while (first < cx_deque_end(self->properties)) {
        cx_deque_iterator pos = first;

        cx_deque_erase(self->properties, pos,
                      (cx_free_func)cpl_property_delete);

    }

    return;

}


/**
 * @brief
 *   Update a property list with a character value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The character value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type <tt>CPL_TYPE_CHAR</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the character value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type @c CPL_TYPE_CHAR.
 */

cpl_error_code
cpl_propertylist_update_char(cpl_propertylist *self, const char *name,
                             char value)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_CHAR);


        cx_assert(property != NULL);

        cpl_property_set_char(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_get_type(property) != CPL_TYPE_CHAR) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        cpl_property_set_char(property, value);

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a boolean value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The boolean value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type <tt>CPL_TYPE_BOOL</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the boolean value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type @c CPL_TYPE_BOOL.
 */

cpl_error_code
cpl_propertylist_update_bool(cpl_propertylist *self, const char *name,
                             int value)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_BOOL);


        cx_assert(property != NULL);

        cpl_property_set_bool(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_get_type(property) != CPL_TYPE_BOOL) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        cpl_property_set_bool(property, value);

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a integer value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The integer value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type <tt>CPL_TYPE_INT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the integer value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type @c CPL_TYPE_INT.
 */

cpl_error_code
cpl_propertylist_update_int(cpl_propertylist *self, const char *name,
                            int value)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_INT);


        cx_assert(property != NULL);

        cpl_property_set_int(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_get_type(property) != CPL_TYPE_INT) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        cpl_property_set_int(property, value);

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a long value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type <tt>CPL_TYPE_LONG</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the long value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type @c CPL_TYPE_LONG.
 */

cpl_error_code
cpl_propertylist_update_long(cpl_propertylist *self, const char *name,
                             long value)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_LONG);


        cx_assert(property != NULL);

        cpl_property_set_long(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_get_type(property) != CPL_TYPE_LONG) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        cpl_property_set_long(property, value);

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a long long value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The long long value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type <tt>CPL_TYPE_LONG_LONG</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the long long value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type @c CPL_TYPE_LONG_LONG.
 */

cpl_error_code
cpl_propertylist_update_long_long(cpl_propertylist *self, const char *name,
                                  long long value)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_LONG_LONG);


        cx_assert(property != NULL);

        cpl_property_set_long_long(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_get_type(property) != CPL_TYPE_LONG_LONG) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        cpl_property_set_long_long(property, value);

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a float value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The float value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type <tt>CPL_TYPE_FLOAT</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the float value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type @c CPL_TYPE_FLOAT.
 */

cpl_error_code
cpl_propertylist_update_float(cpl_propertylist *self, const char *name,
                              float value)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_FLOAT);


        cx_assert(property != NULL);

        cpl_property_set_float(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_get_type(property) != CPL_TYPE_FLOAT) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        cpl_property_set_float(property, value);

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a double value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The double value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type <tt>CPL_TYPE_DOUBLE</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the double value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type @c CPL_TYPE_DOUBLE.
 */

cpl_error_code
cpl_propertylist_update_double(cpl_propertylist *self, const char *name,
                               double value)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_DOUBLE);


        cx_assert(property != NULL);

        cpl_property_set_double(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_get_type(property) != CPL_TYPE_DOUBLE) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        cpl_property_set_double(property, value);

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a string value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The string value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>name</i> or <i>value</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type <tt>CPL_TYPE_STRING</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the string value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type @c CPL_TYPE_STRING.
 */

cpl_error_code
cpl_propertylist_update_string(cpl_propertylist *self, const char *name,
                               const char *value)
{


    cx_deque_iterator pos;


    if (self == NULL || name == NULL || value == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_STRING);


        cx_assert(property != NULL);

        cpl_property_set_string(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_get_type(property) != CPL_TYPE_STRING) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        cpl_property_set_string(property, value);

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a float complex value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The float complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type
           <tt>CPL_TYPE_FLOAT_COMPLEX</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the float complex value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type
 *  @c CPL_TYPE_FLOAT_COMPLEX.
 */

cpl_error_code
cpl_propertylist_update_float_complex(cpl_propertylist *self, const char *name,
                                      float complex value)
{

    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name, CPL_TYPE_FLOAT_COMPLEX);


        cx_assert(property != NULL);

        cpl_property_set_float_complex(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_set_float_complex(property, value)) {
            return cpl_error_set_where_();
        }

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Update a property list with a double complex value.
 *
 * @param self   A property list.
 * @param name   The property name to be assigned to the value.
 * @param value  The double complex value to store.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>name</i> is a <tt>NULL</tt>
 *         pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of type
 *         <tt>CPL_TYPE_DOUBLE_COMPLEX</tt>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function updates the property list @em self with the double complex value
 * @em value. This means, if a property with the name @em name exists already
 * its value is updated, otherwise a property with the name @em name is
 * created and added to @em self. The update will fail if a property with
 * the name @em name exists already which is not of type
 * @c CPL_TYPE_DOUBLE_COMPLEX.
 */

cpl_error_code
cpl_propertylist_update_double_complex(cpl_propertylist *self, const char *name,
                                       double complex value)
{

    cx_deque_iterator pos;


    if (self == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties)) {

        cpl_property *property = cpl_property_new(name,
                                                  CPL_TYPE_DOUBLE_COMPLEX);


        cx_assert(property != NULL);

        cpl_property_set_double_complex(property, value);
        cx_deque_push_back(self->properties, property);
    }
    else {

        cpl_property *property = cx_deque_get(self->properties, pos);


        cx_assert(property != NULL);

        if (cpl_property_set_double_complex(property, value)) {
            cpl_error_set_where_();
        }

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Copy a property from another property list.
 *
 * @param self   A property list.
 * @param other  The property list from which a property is copied.
 * @param name   The name of the property to copy.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>other</i> or <i>name</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The property list <i>other</i> does not contain a property with the
 *         name <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of the same type as the property
 *         which should be copied from <i>other</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function copies the property @em name from the property list
 * @em other to the property list @em self. If the property list @em self
 * does not already contain a property @em name the property is appended
 * to @em self. If a property @em name exists already in @em self the
 * function overwrites the contents of this property if and only if this
 * property is of the same type as the property to be copied from @em other.
 */

cpl_error_code
cpl_propertylist_copy_property(cpl_propertylist *self,
                               const cpl_propertylist *other,
                               const char *name)
{


    cx_deque_iterator spos;
    cx_deque_iterator tpos;


    if (self == NULL || other == NULL || name == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    spos = _cpl_propertylist_find(other, name);

    if (spos == cx_deque_end(other->properties)) {
        return cpl_error_set_message_(CPL_ERROR_DATA_NOT_FOUND, "%s", name);
    }

    tpos = _cpl_propertylist_find(self, name);

    if (tpos == cx_deque_end(self->properties)) {

        cpl_property *p = cpl_property_duplicate(cx_deque_get(other->properties,
                                                             spos));
        cx_deque_push_back(self->properties, p);

    }
    else {

        cpl_property *p = cx_deque_get(self->properties, tpos);
        cpl_property *_p = cx_deque_get(other->properties, spos);


        if (cpl_property_get_type(p) != cpl_property_get_type(_p)) {
            return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);
        }

        switch (cpl_property_get_type(_p)) {

        case CPL_TYPE_CHAR:
            cpl_property_set_char(p, cpl_property_get_char(_p));
            break;

        case CPL_TYPE_BOOL:
            cpl_property_set_bool(p, cpl_property_get_bool(_p));
            break;

        case CPL_TYPE_INT:
            cpl_property_set_int(p, cpl_property_get_int(_p));
            break;

        case CPL_TYPE_LONG:
            cpl_property_set_long(p, cpl_property_get_long(_p));
            break;

        case CPL_TYPE_FLOAT:
            cpl_property_set_float(p, cpl_property_get_float(_p));
            break;

        case CPL_TYPE_DOUBLE:
            cpl_property_set_double(p, cpl_property_get_double(_p));
            break;

        case CPL_TYPE_STRING:
            cpl_property_set_string(p, cpl_property_get_string(_p));
            break;

        case CPL_TYPE_FLOAT_COMPLEX:
            cpl_property_set_float_complex(p,
                                           cpl_property_get_float_complex(_p));
            break;

        case CPL_TYPE_DOUBLE_COMPLEX:
            cpl_property_set_double_complex(p,
                                            cpl_property_get_double_complex(_p));
            break;

        default:
            /* This point should never be reached */
            return cpl_error_set_(CPL_ERROR_UNSUPPORTED_MODE);
            break;
        }

        cpl_property_set_comment(p, cpl_property_get_comment(_p));

    }

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Copy matching properties from another property list.
 *
 * @param self    A property list.
 * @param other   The property list from which a property is copied.
 * @param regexp  The regular expression used to select properties.
 * @param invert  Flag inverting the sense of matching.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success or a CPL error
 *   code otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>other</i> or <i>regexp</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_ILLEGAL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>regexp</i> is an invalid regular expression,
 *         including the empty stringr.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_TYPE_MISMATCH</td>
 *       <td class="ecr">
 *         The property list <i>self</i> contains a property with the
 *         name <i>name</i> which is not of the same type as the property
 *         which should be copied from <i>other</i>.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function copies all properties with matching names from the property
 * list @em other to the property list @em self. If the flag @em invert is
 * zero, all properties whose names match the regular expression @em regexp
 * are copied. If @em invert is set to a non-zero value, all properties with
 * names not matching @em regexp are copied rather. The function expects
 * POSIX 1003.2 compliant extended regular expressions.
 *
 * If the property list @em self does not already contain one of the
 * properties to be copied this property is appended to @em self. If a
 * property to be copied exists already in @em self the function overwrites
 * the contents of this property.
 *
 * Before properties are copied from the property list @em other to @em self
 * the types of the properties are checked and if any type mismatch is
 * detected the function stops processing immediately. The property list
 * @em self is not at all modified in this case.
 *
 * @see cpl_propertylist_copy_property(), cpl_propertylist_append()
 */

cpl_error_code
cpl_propertylist_copy_property_regexp(cpl_propertylist *self,
                                      const cpl_propertylist *other,
                                      const char *regexp,
                                      int invert)
{


    cxint status;

    cxsize i;
    cxsize count = 0;

    cx_deque_const_iterator first, last;

    typedef struct _property_pair_ {
        cpl_property *s;
        cpl_property *t;
    } property_pair;

    property_pair *pairs = NULL;

    cpl_regexp filter;


    if (self == NULL || other == NULL || regexp == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    if (regexp[0] == '\0') {
        return cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
    }

    status = regcomp(&filter.re, regexp, REG_EXTENDED | REG_NOSUB);
    if (status) {
        (void)cpl_error_set_regex(CPL_ERROR_ILLEGAL_INPUT, status,
                                  &filter.re, "regexp='%s', invert=%d",
                                  regexp, invert);
        return CPL_ERROR_ILLEGAL_INPUT;
    }

    filter.invert = invert == 0 ? FALSE : TRUE;


    count = cx_deque_size(other->properties);

    if (count == 0) {
        regfree(&filter.re);
        return CPL_ERROR_NONE;
    }

    pairs = cx_malloc(count * sizeof(property_pair));
    cx_assert(pairs != NULL);

    count = 0;


    first = cx_deque_begin(other->properties);
    last  = cx_deque_end(other->properties);

    while (first != last) {

        cpl_property *p = cx_deque_get(other->properties, first);


        if (_cpl_propertylist_compare_regexp(p, &filter) == TRUE) {

            const cxchar *name = cpl_property_get_name(p);

            cx_deque_const_iterator pos = _cpl_propertylist_find(self, name);

            cpl_property *_p = NULL;


            if (pos != cx_deque_end(self->properties)) {

                _p = cx_deque_get(self->properties, pos);

                if (cpl_property_get_type(p) != cpl_property_get_type(_p)) {

                    regfree(&filter.re);

                    cx_free(pairs);
                    pairs = NULL;

                    return cpl_error_set_(CPL_ERROR_TYPE_MISMATCH);

                }

            }

            pairs[count].s = p;
            pairs[count].t = _p;
            ++count;

        }

        first = cx_deque_next(other->properties, first);

    }

    regfree(&filter.re);


    for (i = 0; i < count; i++) {

        if (pairs[i].t == NULL) {

            cpl_property *p = cpl_property_duplicate(pairs[i].s);
            cx_deque_push_back(self->properties, p);

        }
        else {

            switch (cpl_property_get_type(pairs[i].s)) {

            case CPL_TYPE_CHAR:
                cpl_property_set_char(pairs[i].t,
                                      cpl_property_get_char(pairs[i].s));
                break;

            case CPL_TYPE_BOOL:
                cpl_property_set_bool(pairs[i].t,
                                      cpl_property_get_bool(pairs[i].s));
                break;

            case CPL_TYPE_INT:
                cpl_property_set_int(pairs[i].t,
                                     cpl_property_get_int(pairs[i].s));
                break;

            case CPL_TYPE_LONG:
                cpl_property_set_long(pairs[i].t,
                                      cpl_property_get_long(pairs[i].s));
                break;

            case CPL_TYPE_FLOAT:
                cpl_property_set_float(pairs[i].t,
                                       cpl_property_get_float(pairs[i].s));
                break;

            case CPL_TYPE_DOUBLE:
                cpl_property_set_double(pairs[i].t,
                                        cpl_property_get_double(pairs[i].s));
                break;

            case CPL_TYPE_STRING:
                cpl_property_set_string(pairs[i].t,
                                        cpl_property_get_string(pairs[i].s));
                break;

            case CPL_TYPE_FLOAT_COMPLEX:
                cpl_property_set_float_complex(pairs[i].t,
                                    cpl_property_get_float_complex(pairs[i].s));
                break;

            case CPL_TYPE_DOUBLE_COMPLEX:
                cpl_property_set_double_complex(pairs[i].t,
                                   cpl_property_get_double_complex(pairs[i].s));
                break;

            default:
                /* This point should never be reached */
                cx_free(pairs);
                return cpl_error_set_(CPL_ERROR_UNSUPPORTED_MODE);
                break;
            }

        }

    }

    cx_free(pairs);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Sort a property list.
 *
 * @param self     The property list to sort.
 * @param compare  The function used to compare two properties.
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success, or a CPL error code
 *   otherwise.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function sorts the property list @em self in place, using the function
 * @em compare to determine whether a property is less, equal or greater than
 * another one.
 *
 * The function @em compare must be of the type cpl_propertylist_compare_func.
 *
 * @see
 *   cpl_propertylist_compare_func
 */

/*cpl_error_code
cpl_propertylist_sort(cpl_propertylist* self,
                     int (*compare)(const void*, const void*))*/
cpl_error_code
cpl_propertylist_sort(cpl_propertylist* self,
                     cpl_propertylist_compare_func compare)
{



    cx_compare_func _compare = (cx_compare_func)compare;


    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_deque_sort(self->properties, _compare);

    return CPL_ERROR_NONE;

}


/**
 * @brief
 *   Create a property list from a file.
 *
 * @param name      Name of the input file.
 * @param position  Index of the data set to read.
 *
 * @return
 *   The function returns the newly created property list or @c NULL if an
 *   error occurred.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>name</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_ILLEGAL_INPUT</td>
 *       <td class="ecr">
 *         The <i>position</i> is less than 0 or the properties cannot be
 *         read from the file <i>name</i>.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_FILE_IO</td>
 *       <td class="ecr">
 *         The file <i>name</i> does not exist.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_BAD_FILE_FORMAT</td>
 *       <td class="ecr">
 *         The file <i>name</i> is not a valid FITS file.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The requested data set at index <i>position</i> does not exist.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function reads the properties of the data set with index @em position
 * from the file @em name.
 *
 * Currently only the FITS file format is supported. The property list is
 * created by reading the FITS keywords from extension @em position. The
 * numbering of the data sections starts from 0.
 * When creating the property list from a FITS header, any keyword without
 * a value such as undefined keywords, are not transformed into
 * a property. In the case of float or double (complex) keywords, there is no
 * way to identify the type returned by CFITSIO, therefore this function will
 * always load them as double (complex).
 *
 * @see cpl_propertylist_load_regexp()
 */

cpl_propertylist *
cpl_propertylist_load(const char *name, cpl_size position)
{

    cxint status = 0;

    cpl_propertylist *self;

    cpl_error_code code;

    fitsfile *file = NULL;



    if (name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if ((position < 0) || (position > CX_MAXINT)) {
        cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }

    if (cpl_io_fits_open_diskfile(&file, name, READONLY, &status)) {

        (void)cpl_error_set_fits(status == FILE_NOT_OPENED ?
                                 CPL_ERROR_FILE_IO : CPL_ERROR_BAD_FILE_FORMAT,
                                 status, fits_open_diskfile, "filename='%s', "
                                 "position=%d", name, (cxint)position);
        return NULL;

    }

    self = cpl_propertylist_new();

    code = _cpl_propertylist_fill_from_fits(self, file, (cxint)position,
                                            NULL, NULL);

    if (cpl_io_fits_close_file(file, &status)) {

        code = cpl_error_set_fits(CPL_ERROR_BAD_FILE_FORMAT, status,
                                  fits_close_file, "filename='%s', "
                                  "position=%d", name, (cxint)position);

    } else if (code) {

        /*
         *  Propagate error
         */

        cpl_error_set_where_();

    }

    if (code) {
        cpl_propertylist_delete(self);
        self = NULL;
    }

    return self;

}


/**
 * @brief
 *   Create a filtered property list from a file.
 *
 * @param name      Name of the input file.
 * @param position  Index of the data set to read.
 * @param regexp    Regular expression used to filter properties.
 * @param invert    Flag inverting the sense of matching property names.
 *
 * @return
 *   The function returns the newly created property list or @c NULL if an
 *   error occurred.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>name</i> or the parameter <i>regexp</i> is a
 *         <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_ILLEGAL_INPUT</td>
 *       <td class="ecr">
 *         The <i>position</i> is less than 0, the properties cannot be
 *         read from the file <i>name</i>, or <i>regexp</i> is not a valid
*          extended regular expression, including the empty string..
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_FILE_IO</td>
 *       <td class="ecr">
 *         The file <i>name</i> does not exist.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_BAD_FILE_FORMAT</td>
 *       <td class="ecr">
 *         The file <i>name</i> is not a valid FITS file.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_DATA_NOT_FOUND</td>
 *       <td class="ecr">
 *         The requested data set at index <i>position</i> does not exist.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function reads all properties of the data set with index @em position
 * with matching names from the file @em name. If the flag @em invert is zero,
 * all properties whose names match the regular expression @em regexp are
 * read. If @em invert is set to a non-zero value, all properties with
 * names not matching @em regexp are read rather. The function expects
 * POSIX 1003.2 compliant extended regular expressions.
 *
 * Currently only the FITS file format is supported. The property list is
 * created by reading the FITS keywords from extension @em position.
 * The numbering of the data sections starts from 0.
 *
  * When creating the property list from a FITS header, any keyword without
 * a value such as undefined keywords, are not transformed into
 * a property. In the case of float or double (complex) keywords, there is no
 * way to identify the type returned by CFITSIO, therefore this function will
 * always load them as double (complex).
 *
 * FITS format specific keyword prefixes (e.g. @c HIERARCH) must
 * not be part of the given pattern string @em regexp, but only the actual
 * FITS keyword name may be given.
 *
 * @see cpl_propertylist_load()
 */

cpl_propertylist *
cpl_propertylist_load_regexp(const char *name, cpl_size position,
                             const char *regexp, int invert)
{

    cxint status = 0;

    cpl_propertylist *self;

    cpl_regexp filter;

    cpl_error_code code;

    fitsfile *file = NULL;


    if (name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if ((position < 0) || (position > CX_MAXINT)) {
        cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }


    if (regexp == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (regexp[0] == '\0') {
        cpl_error_set_(CPL_ERROR_ILLEGAL_INPUT);
        return NULL;
    }


    /*
     * Set up the regular expression filter
     */

    status = regcomp(&filter.re, regexp, REG_EXTENDED | REG_NOSUB);

    if (status) {

        (void)cpl_error_set_regex(CPL_ERROR_ILLEGAL_INPUT, status,
                                  &filter.re, "regexp='%s', invert=%d",
                                  regexp, invert);
        return NULL;

    }

    filter.invert = invert == 0 ? FALSE : TRUE;

    if (cpl_io_fits_open_diskfile(&file, name, READONLY, &status)) {

        (void)cpl_error_set_fits(status == FILE_NOT_OPENED ?
                                 CPL_ERROR_FILE_IO : CPL_ERROR_BAD_FILE_FORMAT,
                                 status, fits_open_diskfile, "filename='%s', "
                                 "position=%d, regexp='%s'",
                                 name, (cxint)position, regexp);

        status = 0;
        cpl_io_fits_close_file(file, &status);
        regfree(&filter.re);

        return NULL;

    }

    self = cpl_propertylist_new();

    code = _cpl_propertylist_fill_from_fits(self, file, (cxint)position,
                                           _cpl_propertylist_filter_regexp,
                                           &filter);
    regfree(&filter.re);

    if (cpl_io_fits_close_file(file, &status)) {

        code = cpl_error_set_fits(CPL_ERROR_BAD_FILE_FORMAT, status,
                                  fits_close_file, "filename='%s', position="
                                  "%d, regexp='%s'", name, (cxint)position,
                                  regexp);

    } else if (code) {

        /*
         *  Propagate error
         */

        cpl_error_set_(code);

    }

    if (code) {
        cpl_propertylist_delete(self);
        self = NULL;
    }

    return self;

}


/**
 * @internal
 * @brief
 *   Write a property list to a FITS file.
 *
 * @param file         The FITS file to write to.
 * @param properties  The property list to write to the file.
 * @param to_rm       The regular expression used to filter properties
 *
 * @return
 *   The function returns @c CPL_ERROR_NONE on success, or an appropriate
 *   error code otherwise.
 *
 * This function takes a sorted property list and appends it to the
 * provided FITS file. All properties with names matching the regular expression
 * @c to_rm will not be copied to the output. The function expects
 * POSIX 1003.2 compliant extended regular expressions.
 */

cpl_error_code
cpl_propertylist_to_fitsfile(fitsfile *file,
                             const cpl_propertylist *self,
                             const char *to_rm)
{

    cpl_regexp filter;

    cpl_error_code error;

    if (file == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    if (self == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    if (to_rm != NULL) {

        /* A regular expression must be applied */

        const cxint rstatus = regcomp(&filter.re, to_rm,
                                      REG_EXTENDED | REG_NOSUB);

        if (rstatus) {
            return cpl_error_set_regex(CPL_ERROR_ILLEGAL_INPUT, rstatus,
                                       &filter.re, "to_rm='%s'", to_rm);
        }

        filter.invert = 0;

    }

    error = _cpl_propertylist_to_fitsfile(file, self,
                                          _cpl_propertylist_filter_regexp,
                                          to_rm ? &filter : NULL);
    if (to_rm != NULL) {
        regfree(&filter.re);
    }

    return error ? cpl_error_set_where_() : CPL_ERROR_NONE;

}


/**
 * @internal
 * @brief
 *   Create a property list from a cfitsio FITS file.
 *
 * @param file  FITS file pointing to the header to convert.
 *
 * @return The function returns the created property list, or @c NULL
 *   in case of an error.
 *
 * The function  converts the FITS keywords from the current HDU of the FITS
 * file @em file into properties and puts them into a property list.
 *
 * The special FITS keyword END indicating the end of a FITS header is not
 * transformed into a property, but simply ignored.
 *
 * In case of an error, an appropriate error code is set. If a FITS header
 * card cannot be parsed the error code is set to @c CPL_ERROR_ILLEGAL_INPUT
 * or to @c CPL_ERROR_INVALID_TYPE if a FITS keyword type is not supported.
 * If @em file is @c NULL the error code is set to @c CPL_ERROR_NULL_INPUT.
 */

cpl_propertylist *
cpl_propertylist_from_fitsfile(fitsfile *file)
{

    cpl_propertylist *self;

    cpl_ensure(file != NULL, CPL_ERROR_NULL_INPUT, NULL);

    self = cpl_propertylist_new();

    if (_cpl_propertylist_fill_from_fits(self, file, 0, NULL, NULL)) {
        cpl_propertylist_delete(self);
        self = NULL;
        cpl_error_set_where_();
    }

    return self;

}


/**
 * @brief
 *   Access property list elements by property name.
 *
 * @param self      The property list to query.
 * @param name  The name of the property to retrieve.
 *
 * @return
 *   The function returns the property with name @em name, or @c NULL
 *   if it does not exist.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or the <i>name</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns a handle to the property list element, the property,
 * with the name @em name. If more than one property exist with the same
 * @em name, then the first one found will be returned.
 */

const cpl_property *
cpl_propertylist_get_property_const(const cpl_propertylist *self, const char *name)
{
    cx_deque_iterator pos;

    if (self == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    if (name == NULL) {
        cpl_error_set_(CPL_ERROR_NULL_INPUT);
        return NULL;
    }

    pos = _cpl_propertylist_find(self, name);

    if (pos == cx_deque_end(self->properties))
        return NULL;

    return cx_deque_get(self->properties, pos);

}


/**
 * @brief
 *   Access property list elements by property name.
 *
 * @param self      The property list to query.
 * @param name  The name of the property to retrieve.
 *
 * @return
 *   The function returns the property with name @em name, or @c NULL
 *   if it does not exist.
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or the <i>name</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function returns a handle to the property list element, the property,
 * with the name @em name. If more than one property exist with the same
 * @em name, then the first one found will be returned.
 */

cpl_property *
cpl_propertylist_get_property(cpl_propertylist *self, const char *name)
{

    cpl_errorstate prestate = cpl_errorstate_get();
    cpl_property *property = (cpl_property *)cpl_propertylist_get_property_const(self, name);

    if (!cpl_errorstate_is_equal(prestate))
        (void)cpl_error_set_where_();

    return property;

}


/**
 * @brief
 *   Save a property list to a FITS file
 *
 * @param self      The property list to save or NULL if empty
 * @param filename  Name of the file to write
 * @param mode      The desired output options (combined with bitwise or)
 *
 * @return
 *   CPL_ERROR_NONE or the relevant #_cpl_error_code_ on error
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The <i>filename</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_ILLEGAL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>mode</i> is invalid.
 *       </td>
 *     </tr>
 *     <tr>
 *       <td class="ecl">CPL_ERROR_FILE_IO</td>
 *       <td class="ecr">
 *         The file cannot be written or accessed.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * This function saves a property list to a FITS file, using cfitsio.
 * The data unit is empty.
 *
 * Supported output modes are CPL_IO_CREATE (create a new file) and
 * CPL_IO_EXTEND  (append to an existing file)
 */

cpl_error_code
cpl_propertylist_save(const cpl_propertylist *self, const char *filename,
                      unsigned mode)
{

    const cxchar *badkeys = mode & CPL_IO_EXTEND ?
                            CPL_FITS_BADKEYS_EXT  "|" CPL_FITS_COMPRKEYS :
                            CPL_FITS_BADKEYS_PRIM "|" CPL_FITS_COMPRKEYS;

    cxint error = 0;

    cpl_error_code code = CPL_ERROR_NONE;

    fitsfile *fptr;



    /*
     *  Check entries
     */

    cpl_ensure_code(filename, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code(mode == CPL_IO_CREATE || mode == CPL_IO_EXTEND,
                    CPL_ERROR_ILLEGAL_INPUT);


    /* Switch on the mode */
    if (mode == CPL_IO_EXTEND) {
        /* Append */

        /* Open the file */
        if (cpl_io_fits_open_diskfile(&fptr, filename, READWRITE, &error)) {
            return cpl_error_set_fits(CPL_ERROR_FILE_IO, error,
                                      fits_open_diskfile, "filename='%s', "
                                      "mode=%u", filename, mode);
        }

    }
    else {
        /* Main HDU */

        /* Create the file */
        cxchar *sval = cpl_sprintf("!%s", filename);

        cpl_io_fits_create_file(&fptr, sval, &error);
        cpl_free(sval);

        if (error != 0) {
            return cpl_error_set_fits(CPL_ERROR_FILE_IO, error,
                                      fits_create_file, "filename='%s', "
                                      "mode=%u", filename, mode);
        }

    }


    /* Create empty header */
    if (fits_create_img(fptr, BYTE_IMG, 0, NULL, &error)) {

        cxint error2 = 0;

        cpl_io_fits_close_file(fptr, &error2);

        return cpl_error_set_fits(CPL_ERROR_FILE_IO, error, fits_create_img,
                                  "filename='%s', mode=%u", filename, mode);

    }


    /* Add DATE */
    if (mode != CPL_IO_EXTEND && fits_write_date(fptr, &error)) {

        cxint error2 = 0;

        cpl_io_fits_close_file(fptr, &error2);

        return cpl_error_set_fits(CPL_ERROR_FILE_IO, error, fits_write_date,
                                  "filename='%s', mode=%u", filename, mode);

    }


    /* Add the property list */
    if (cpl_fits_add_properties(fptr, self, badkeys)) {
        code = cpl_error_set_where_();
    }


    /* Close (and write to disk) */
    if (cpl_io_fits_close_file(fptr, &error)) {
        return cpl_error_set_fits(CPL_ERROR_FILE_IO, error, fits_close_file,
                                  "filename='%s', mode=%u", filename, mode);
    }

    return code;

}


/**
 * @brief
 *   Print a property list.
 *
 * @param self    Pointer to property list
 * @param stream   The output stream
 *
 * @return Nothing.
 *
 * This function is mainly intended for debug purposes.
 * If the specified stream is @c NULL, it is set to @em stdout.
 * The function used for printing is the standard C @c fprintf().
 */

void cpl_propertylist_dump(const cpl_propertylist *self, FILE *stream)
{

    cxchar c = '\0';

    cpl_size i = 0;
    cpl_size sz = cpl_propertylist_get_size(self);


    if (stream == NULL)
        stream = stdout;

    if (self == NULL) {
        fprintf(stream, "NULL property list\n\n");
        return;
    }


    fprintf(stream, "Property list at address %p:\n", (cxcptr)self);

    for (i = 0; i < sz; i++) {
        const cpl_property *p = cpl_propertylist_get_const(self, i);
        const cxchar *name = cpl_property_get_name(p);
        const cxchar *comment = cpl_property_get_comment(p);
        cpl_size size = cpl_property_get_size(p);
        cpl_type type = cpl_property_get_type(p);
        const cxchar *typestr = cpl_type_get_name(type);

        fprintf(stream, "Property at address %p\n", (cxcptr)p);
        fprintf(stream, "\tname   : %p '%s'\n", name, name);
        fprintf(stream, "\tcomment: %p '%s'\n", comment, comment);
        fprintf(stream, "\ttype   : %#09x '%s'\n", type, typestr);
        fprintf(stream, "\tsize   : %" CPL_SIZE_FORMAT "\n", size);
        fprintf(stream, "\tvalue  : ");


        switch (type) {
        case CPL_TYPE_CHAR:
            c = cpl_property_get_char(p);
            if (!c)
                fprintf(stream, "''");
            else
                fprintf(stream, "'%c'", c);
            break;

        case CPL_TYPE_BOOL:
            fprintf(stream, "%d", cpl_property_get_bool(p));
            break;

        case CPL_TYPE_INT:
            fprintf(stream, "%d", cpl_property_get_int(p));
            break;

        case CPL_TYPE_LONG:
            fprintf(stream, "%ld", cpl_property_get_long(p));
            break;

        case CPL_TYPE_LONG_LONG:
            fprintf(stream, "%lld", cpl_property_get_long_long(p));
            break;

        case CPL_TYPE_FLOAT:
            fprintf(stream, "%.7g", cpl_property_get_float(p));
            break;

        case CPL_TYPE_DOUBLE:
            fprintf(stream, "%.15g", cpl_property_get_double(p));
            break;

        case CPL_TYPE_STRING:
            fprintf(stream, "'%s'", cpl_property_get_string(p));
            break;

        case CPL_TYPE_FLOAT_COMPLEX: {
            /* Print with FITS header format */
            const float complex z = cpl_property_get_float_complex(p);
            fprintf(stream, "(%.7g,%.7g)", crealf(z), cimagf(z));
            break;
        }
        case CPL_TYPE_DOUBLE_COMPLEX: {
            /* Print with FITS header format */
            const double complex z = cpl_property_get_double_complex(p);
            fprintf(stream, "(%.15g,%.15g)", creal(z), cimag(z));
            break;
        }
        default:
            fprintf(stream, "unknown.");
            break;

        }

        fprintf(stream, "\n");

    }

    return;

}


/**
 * @brief
 *    Append a property to a property list
 *
 * @param   self      Property list to append to
 * @param   property  The property to append
 *
 * @return
 *    The function returns @c CPL_ERROR_NONE on success or a CPL error
 *    code otherwise.
 *
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>property</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * This function creates a new property and appends it to the end of a property list.
 * It will not check if the property already exists.
 */

cpl_error_code
cpl_propertylist_append_property(cpl_propertylist *self,
                                 const cpl_property *property)
{

    if((self == NULL) || (property == NULL)){
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_deque_push_back(self->properties, cpl_property_duplicate(property));

    return CPL_ERROR_NONE;
}


/**
 * @brief
 *    Prepend a property to a property list
 *
 * @param   self      Property list to prepend to
 * @param   property  The property to prepend
 *
 * @return
 *    The function returns @c CPL_ERROR_NONE on success or a CPL error
 *    code otherwise.
 *
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i> or <i>property</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * This function creates a new property and prepends it to the beginning of a property list.
 * It will not check if the property already exists.
 */

cpl_error_code
cpl_propertylist_prepend_property(cpl_propertylist *self,
                                  const cpl_property *property)
{

    if((self == NULL) || (property == NULL)){
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    cx_deque_push_front(self->properties, cpl_property_duplicate(property));


    return CPL_ERROR_NONE;
}


/**
 * @brief
 *    Insert a property into a property list at the given position
 *
 * @param   self      Property list
 * @param   here      Name indicating the position where to insert the property
 * @param   property  The property to insert
 *
 * @return
 *    The function returns @c CPL_ERROR_NONE on success or a CPL error
 *    code otherwise.
 *
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>property</i> or <i>here</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new property and inserts it into the property list
 * @em self at the position of the property named @em here.
 */

cpl_error_code
cpl_propertylist_insert_property(cpl_propertylist *self,
                                 const char *here,
                                 const cpl_property *property)
{

    cx_deque_iterator pos;

    if (self == NULL || here == NULL || property == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    /*
     * Find the position where property should be inserted.
     */

    pos = _cpl_propertylist_find(self, here);

    if (pos == cx_deque_end(self->properties)) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }


    /*
     * Insert it into the deque
     */

    cx_deque_insert(self->properties, pos, cpl_property_duplicate(property));


    return CPL_ERROR_NONE;

}


/**
 * @brief
 *    Insert a property into a property list after the given position
 *
 * @param   self      Property list
 * @param   after     Name of the property after which to insert the property
 * @param   property  The property to insert
 *
 * @return
 *    The function returns @c CPL_ERROR_NONE on success or a CPL error
 *    code otherwise.
 *
 *
 * @error
 *   <table class="ec" align="center">
 *     <tr>
 *       <td class="ecl">CPL_ERROR_NULL_INPUT</td>
 *       <td class="ecr">
 *         The parameter <i>self</i>, <i>property</i> or <i>after</i> is a <tt>NULL</tt> pointer.
 *       </td>
 *     </tr>
 *   </table>
 * @enderror
 *
 * The function creates a new property and inserts it into the property list
 * @em self after the position of the property named @em after.
 */

cpl_error_code
cpl_propertylist_insert_after_property(cpl_propertylist *self,
                                 const char *after,
                                 const cpl_property *property)
{

    cx_deque_iterator pos;

    if (self == NULL || after == NULL || property == NULL) {
        return cpl_error_set_(CPL_ERROR_NULL_INPUT);
    }

    /*
     * Find the position where property should be inserted.
     */
    pos = _cpl_propertylist_find(self, after);

    if (pos == cx_deque_end(self->properties)) {
        return cpl_error_set_(CPL_ERROR_UNSPECIFIED);
    }

    /*
     * Get the position after it
     */
    pos = cx_deque_next(self->properties, pos);

    /*
     * Insert it into the deque
     */
    cx_deque_insert(self->properties, pos, cpl_property_duplicate(property));


    return CPL_ERROR_NONE;

}
/**@}*/
