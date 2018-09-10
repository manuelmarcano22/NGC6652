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

/* Author: Artur Szostak <aszostak@partner.eso.org> */

/*
 * These are additional tests relying on overloading certain functions to mock
 * various failure modes that are difficult to trigger in any other way.
 * The overloading is not compatible with esorex_json_parse_errors*-test.c,
 * therefore we needed to compile them into their own test binary.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef ENABLE_PYTHON_RECIPES

/* Undefine NDEBUG to make sure assert() macros are compiled in for these
   tests. */
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cpl_test.h>
#include "er_json.c"

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void test_add_recipeconfig_tag(void);

/*-----------------------------------------------------------------------------
                            Overloaded functions
 -----------------------------------------------------------------------------*/

/* Overload the following functions to produce an error. */

int cpl_recipeconfig_set_tag(cpl_recipeconfig* self, const char* tag,
                             cpl_size min_count, cpl_size max_count)
{
    (void) self;
    (void) tag;
    (void) min_count;
    (void) max_count;
    return cpl_error_set_message("cpl_recipeconfig_set_tag",
                                 CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
}

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    test_add_recipeconfig_tag();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

static void test_add_recipeconfig_tag(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;
    cpl_recipeconfig * config = NULL;

    /* Test error handling in add_recipeconfig_tag for an edge case where
       the cpl_recipeconfig_set_tag function fails. */
    json =
        "{\n"
        "    \"tag\": \"RAW\",\n"
        "    \"min\": 1,\n"
        "    \"max\": 1,\n"
        "    \"inputs\": [\n"
        "         {\n"
        "             \"tag\": \"CALIB\",\n"
        "             \"min\": 2,\n"
        "             \"max\": 3\n"
        "         }\n"
        "      ],\n"
        "    \"outputs\": [\"PROD\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
