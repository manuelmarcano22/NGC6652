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

static void test_json_to_plugin(void);
static void test_add_recipeconfig_inputs(void);
static void test_add_recipeconfig_outputs(void);

/*-----------------------------------------------------------------------------
                            Overloaded functions
 -----------------------------------------------------------------------------*/

/* Overload the following functions to produce an error. */

cpl_error_code cpl_plugin_init(cpl_plugin *self, unsigned int api,
                               unsigned long version, unsigned long type,
                               const char *name, const char *synopsis,
                               const char *description, const char *author,
                               const char *email, const char *copyright,
                               cpl_plugin_func initialize,
                               cpl_plugin_func execute,
                               cpl_plugin_func deinitialize)
{
    (void) self;
    (void) api;
    (void) version;
    (void) type;
    (void) name;
    (void) synopsis;
    (void) description;
    (void) author;
    (void) email;
    (void) copyright;
    (void) initialize;
    (void) execute;
    (void) deinitialize;
    return cpl_error_set_message("cpl_plugin_init", CPL_ERROR_ILLEGAL_OUTPUT,
                                 "dummy error");
}

int
cpl_recipeconfig_set_input(cpl_recipeconfig* self,
                           const char* tag, const char* input,
                           cpl_size min_count, cpl_size max_count)
{
    (void) self;
    (void) tag;
    (void) input;
    (void) min_count;
    (void) max_count;
    return cpl_error_set_message("cpl_recipeconfig_set_input",
                                 CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
}

int
cpl_recipeconfig_set_output(cpl_recipeconfig* self,
                            const char* tag, const char* output)
{
    (void) self;
    (void) tag;
    (void) output;
    return cpl_error_set_message("cpl_recipeconfig_set_output",
                                 CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
}

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    test_json_to_plugin();
    test_add_recipeconfig_inputs();
    test_add_recipeconfig_outputs();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

static void test_json_to_plugin(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;

    /* Test edge case where failure in cpl_plugin_init is handled correctly by
       er_json_to_plugin. */
    json =
        "{\n"
        "    \"class\": \"testclass\",\n"
        "    \"name\": \"testrecipe\",\n"
        "    \"version\": 123,\n"
        "    \"synopsis\": \"test plugin\",\n"
        "    \"description\": \"description ...\",\n"
        "    \"author\": \"Some Author\",\n"
        "    \"email\": \"someone@eso.org\",\n"
        "    \"copyright\": \"copyright ...\",\n"
        "    \"parameters\": [],\n"
        "    \"frames\": []\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);
    er_json_node_delete(parsetree);
}

static void test_add_recipeconfig_inputs(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;
    cpl_recipeconfig * config = NULL;

    /* Test error handling in add_recipeconfig_inputs for an edge case where
       the cpl_recipeconfig_set_input function fails. */
    json =
        "[\n"
        "  {\n"
        "    \"tag\": \"CALIB\",\n"
        "    \"max\": 1,\n"
        "    \"min\": 1\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_inputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);
}

static void test_add_recipeconfig_outputs(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;
    cpl_recipeconfig * config = NULL;

    /* Test error handling in add_recipeconfig_outputs for an edge case where
       the cpl_recipeconfig_set_output function fails. */
    json = "[\"PROD\"]";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_outputs(parsetree, json, config, "RAW"),
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
