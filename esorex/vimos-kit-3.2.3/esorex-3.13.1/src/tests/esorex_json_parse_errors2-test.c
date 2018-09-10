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
static void test_json_to_parameter(void);
static void test_json_to_parameterlist(void);
static void test_json_to_frameset(void);

/*-----------------------------------------------------------------------------
                            Overloaded functions
 -----------------------------------------------------------------------------*/

/* Overload the following function to return a specific value that can be
   controlled by these tests, or to produce an error. */

static cpl_parameter * cpl_parameter_new_value_return_value = NULL;

cpl_parameter *
cpl_parameter_new_value(const char *name, cpl_type type,
                        const char *description, const char *context, ...)
{
    (void) name;
    (void) type;
    (void) description;
    (void) context;
    if (cpl_parameter_new_value_return_value == NULL)
    {
        cpl_error_set_message("cpl_parameter_new_value",
                              CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
    }
    return cpl_parameter_new_value_return_value;
}

/* Overload the following functions to produce an error. */

cpl_error_code cpl_parameterlist_append(cpl_parameterlist *self,
                                        cpl_parameter *parameter)
{
    (void) self;
    (void) parameter;
    return cpl_error_set_message("cpl_parameterlist_append",
                                 CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
}

cpl_error_code cpl_frameset_insert(cpl_frameset *self, cpl_frame *frame)
{
    (void) self;
    (void) frame;
    return cpl_error_set_message("cpl_frameset_insert",
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
    test_json_to_parameter();
    test_json_to_parameterlist();
    test_json_to_frameset();

    /* Setting this to NULL helps in finding memory leaks with Valgrind,
       because the memory will not be marked as still reachable. */
    cpl_parameter_new_value_return_value = NULL;

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

    /* Test edge case where failure in cpl_parameter_new_value is handled
       appropriately by er_json_to_plugin. */
    cpl_parameter_new_value_return_value = NULL;
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

    /* Test edge case where failure in cpl_parameterlist_append is handled
       appropriately by er_json_to_plugin. */
    cpl_parameter_new_value_return_value = cpl_parameter_new_enum("__class__",
            CPL_TYPE_STRING, "Python class name", "__python__", "testclass",
            1, "testclass");
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

static void test_json_to_parameter(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;

    /* Test edge case where failure in cpl_parameter_new_value is handled
       appropriately by json_to_parameter. */
    cpl_parameter_new_value_return_value = NULL;
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);
    er_json_node_delete(parsetree);
}

static void test_json_to_parameterlist(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;

    /* Test edge case where failure in cpl_parameterlist_append is handled
       appropriately by json_to_parameterlist. */
    json =
        "[\n"
        "  {\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"choices\": [3]\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameterlist(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);
    er_json_node_delete(parsetree);
}

static void test_json_to_frameset(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;

    /* Test edge case where failure in cpl_frameset_insert is handled
       appropriately by json_to_frameset. */
    json =
        "[\n"
        "  {\n"
        "    \"filename\": \"test1.fits\",\n"
        "    \"tag\": \"RAW\",\n"
        "    \"type\": 1,\n"
        "    \"group\": 0,\n"
        "    \"level\": 0\n"
        "  },\n"
        "  {\n"
        "    \"filename\": \"test2.fits\",\n"
        "    \"tag\": \"CALIB\",\n"
        "    \"type\": 1,\n"
        "    \"group\": 0,\n"
        "    \"level\": 0\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_frameset(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);
    er_json_node_delete(parsetree);
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
