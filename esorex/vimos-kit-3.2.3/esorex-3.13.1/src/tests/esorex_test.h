/*
 * This file is part of the ESO Recipe Execution Tool
 * Copyright (C) 2017 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

#include <cpl_test.h>

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* Tests if the two string pointers are both NULL or otherwise the strings
   themselves are both equal. */
#define test_equal_string(first, second) \
    if (! (first == NULL && second == NULL)) { \
        cpl_test_eq_string(first, second); \
    }

/* Tests if two CPL parameter objects are equal or not, i.e. are all their
   attributes the same or not. */
#define test_equal_param(first, second) \
    cpl_test_eq(cpl_parameter_get_class(first), \
                cpl_parameter_get_class(second)); \
    cpl_test_eq(cpl_parameter_get_type(first), \
                cpl_parameter_get_type(second)); \
    cpl_test_eq_string(cpl_parameter_get_name(first), \
                       cpl_parameter_get_name(second)); \
    cpl_test_eq_string(cpl_parameter_get_help(first), \
                       cpl_parameter_get_help(second)); \
    cpl_test_eq_string(cpl_parameter_get_context(first), \
                       cpl_parameter_get_context(second)); \
    test_equal_string(cpl_parameter_get_tag(first), \
                      cpl_parameter_get_tag(second)); \
    cpl_test_eq(cpl_parameter_get_default_flag(first), \
                cpl_parameter_get_default_flag(second)); \
    cpl_test_eq(cpl_parameter_is_enabled(first, CPL_PARAMETER_MODE_CLI), \
                cpl_parameter_is_enabled(second, CPL_PARAMETER_MODE_CLI)); \
    cpl_test_eq(cpl_parameter_is_enabled(first, CPL_PARAMETER_MODE_ENV), \
                cpl_parameter_is_enabled(second, CPL_PARAMETER_MODE_ENV)); \
    cpl_test_eq(cpl_parameter_is_enabled(first, CPL_PARAMETER_MODE_CFG), \
                cpl_parameter_is_enabled(second, CPL_PARAMETER_MODE_CFG)); \
    test_equal_string(cpl_parameter_get_alias(first, CPL_PARAMETER_MODE_CLI), \
                      cpl_parameter_get_alias(second, CPL_PARAMETER_MODE_CLI));\
    test_equal_string(cpl_parameter_get_alias(first, CPL_PARAMETER_MODE_ENV), \
                      cpl_parameter_get_alias(second, CPL_PARAMETER_MODE_ENV));\
    test_equal_string(cpl_parameter_get_alias(first, CPL_PARAMETER_MODE_CFG), \
                      cpl_parameter_get_alias(second, CPL_PARAMETER_MODE_CFG));\
    switch (cpl_parameter_get_type(first)) { \
    case CPL_TYPE_BOOL: \
        cpl_test_eq(cpl_parameter_get_bool(first), \
                    cpl_parameter_get_bool(second)); \
        cpl_test_eq(cpl_parameter_get_default_bool(first), \
                    cpl_parameter_get_default_bool(second)); \
        break; \
    case CPL_TYPE_INT: \
        cpl_test_eq(cpl_parameter_get_int(first), \
                    cpl_parameter_get_int(second)); \
        cpl_test_eq(cpl_parameter_get_default_int(first), \
                    cpl_parameter_get_default_int(second)); \
        break; \
    case CPL_TYPE_DOUBLE: \
        cpl_test_abs(cpl_parameter_get_double(first), \
                     cpl_parameter_get_double(second), DBL_EPSILON); \
        cpl_test_abs(cpl_parameter_get_default_double(first), \
                     cpl_parameter_get_default_double(second), DBL_EPSILON); \
        break; \
    case CPL_TYPE_STRING: \
        cpl_test_eq_string(cpl_parameter_get_string(first), \
                           cpl_parameter_get_string(second)); \
        cpl_test_eq_string(cpl_parameter_get_default_string(first), \
                           cpl_parameter_get_default_string(second)); \
        break; \
    default: \
        cpl_test_assert( \
            cpl_parameter_get_type(first) == CPL_TYPE_BOOL || \
            cpl_parameter_get_type(first) == CPL_TYPE_INT || \
            cpl_parameter_get_type(first) == CPL_TYPE_DOUBLE || \
            cpl_parameter_get_type(first) == CPL_TYPE_STRING \
        ); \
    } \
    switch (cpl_parameter_get_class(first)) { \
    case CPL_PARAMETER_CLASS_RANGE: \
        switch (cpl_parameter_get_type(first)) { \
        case CPL_TYPE_INT: \
            cpl_test_eq(cpl_parameter_get_range_min_int(first), \
                        cpl_parameter_get_range_min_int(second)); \
            cpl_test_eq(cpl_parameter_get_range_max_int(first), \
                        cpl_parameter_get_range_max_int(second)); \
            break; \
        case CPL_TYPE_DOUBLE: \
            cpl_test_abs(cpl_parameter_get_range_min_double(first), \
                         cpl_parameter_get_range_min_double(second), \
                         DBL_EPSILON); \
            cpl_test_abs(cpl_parameter_get_range_max_double(first), \
                         cpl_parameter_get_range_max_double(second), \
                         DBL_EPSILON); \
            break; \
        default: \
            cpl_test_assert( \
                cpl_parameter_get_type(first) == CPL_TYPE_INT || \
                cpl_parameter_get_type(first) == CPL_TYPE_DOUBLE \
            ); \
        } \
        break; \
    case CPL_PARAMETER_CLASS_ENUM: \
        cpl_test_eq(cpl_parameter_get_enum_size(first), \
                    cpl_parameter_get_enum_size(second)); \
        switch (cpl_parameter_get_type(first)) { \
        case CPL_TYPE_INT: \
            { \
                int n = 0; \
                for (n = 0; n < cpl_parameter_get_enum_size(first); ++n) { \
                    cpl_test_eq(cpl_parameter_get_enum_int(first, n), \
                                cpl_parameter_get_enum_int(second, n)); \
                } \
            } \
            break; \
        case CPL_TYPE_DOUBLE: \
            { \
                int n = 0; \
                for (n = 0; n < cpl_parameter_get_enum_size(first); ++n) { \
                    cpl_test_abs(cpl_parameter_get_enum_double(first, n), \
                                 cpl_parameter_get_enum_double(second, n), \
                                 DBL_EPSILON); \
                } \
            } \
            break; \
        case CPL_TYPE_STRING: \
            { \
                int n = 0; \
                for (n = 0; n < cpl_parameter_get_enum_size(first); ++n) { \
                    cpl_test_eq_string( \
                            cpl_parameter_get_enum_string(first, n), \
                            cpl_parameter_get_enum_string(second, n) \
                        ); \
                } \
            } \
            break; \
        default: \
            cpl_test_assert( \
                cpl_parameter_get_type(first) == CPL_TYPE_INT || \
                cpl_parameter_get_type(first) == CPL_TYPE_DOUBLE || \
                cpl_parameter_get_type(first) == CPL_TYPE_STRING \
            ); \
        } \
        break; \
    default: \
        cpl_test_assert( \
            cpl_parameter_get_class(first) == CPL_PARAMETER_CLASS_VALUE || \
            cpl_parameter_get_class(first) == CPL_PARAMETER_CLASS_RANGE || \
            cpl_parameter_get_class(first) == CPL_PARAMETER_CLASS_ENUM \
        ); \
    } \
    cpl_test_eq_error(cpl_error_get_code(), CPL_ERROR_NONE)

/* Tests if two CPL frame objects are equal or not, i.e. are all their
   attributes the same or not. */
#define test_equal_frame(first, second) \
    cpl_test_eq_string(cpl_frame_get_filename(first), \
                       cpl_frame_get_filename(second)); \
    cpl_test_eq_string(cpl_frame_get_tag(first), \
                       cpl_frame_get_tag(second)); \
    cpl_test_eq(cpl_frame_get_type(first), \
                cpl_frame_get_type(second)); \
    cpl_test_eq(cpl_frame_get_group(first), \
                cpl_frame_get_group(second)); \
    cpl_test_eq(cpl_frame_get_level(first), \
                cpl_frame_get_level(second))
