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
 * These tests rely on overloading certain standard library and 3rd party
 * functions to mock various failure modes that are difficult to trigger in any
 * other way.
 */

/* Undefine NDEBUG to make sure assert() macros are compiled in for these
   tests. */
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <errno.h>
#include <cpl_test.h>
#include "er_json.c"

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void test_parse_number(void);
static void test_parameter_to_json(void);
static void test_plugin_to_json(void);

#ifdef ENABLE_PYTHON_RECIPES
static void test_json_to_pluginlist(void);
#endif

/*-----------------------------------------------------------------------------
                            Overloaded functions
 -----------------------------------------------------------------------------*/

/* Overload the strtod function so that we can control if it should generate
   and error or behave as normal. */

static cpl_boolean strtod_produce_error = CPL_FALSE;

double strtod(const char * text, char ** endptr)
{
    if (strtod_produce_error)
    {
        if (endptr != NULL) *endptr = (char *) text;
        errno = EINVAL;
        return 0.0;
    }
    else
    {
        /* Use strtof as a surrogate to reproduce the original behaviour of
           this function. Since we do not use very large or small numbers in
           these tests this approximation is good enough. */
        double result = strtof(text, endptr);
        return result;
    }
}

/* Overload the following functions to return a specific value that can be
   controlled by these tests. */

static cpl_type cpl_parameter_get_type_return_value = CPL_TYPE_INVALID;

cpl_type cpl_parameter_get_type(const cpl_parameter *self)
{
    (void) self;
    return cpl_parameter_get_type_return_value;
}

static cpl_parameter_class
cpl_parameter_get_class_return_value = CPL_PARAMETER_CLASS_INVALID;

cpl_parameter_class cpl_parameter_get_class(const cpl_parameter *self)
{
    (void) self;
    return cpl_parameter_get_class_return_value;
}

static int cpl_parameter_get_enum_size_return_value = 0;

int cpl_parameter_get_enum_size(const cpl_parameter *self)
{
    (void) self;
    return cpl_parameter_get_enum_size_return_value;
}

static char** cpl_recipeconfig_get_tags_return_value = NULL;

char** cpl_recipeconfig_get_tags(const cpl_recipeconfig* self)
{
    (void) self;
    /* If the value to return is NULL then return NULL and set and error.
       Otherwise depp copy the values to reproduce the behaviour of the original
       function. */
    if (cpl_recipeconfig_get_tags_return_value == NULL)
    {
        cpl_error_set_message("cpl_recipeconfig_get_tags",
                              CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
        return NULL;
    }
    cpl_size count = 0;
    char ** str = NULL;
    for (str = cpl_recipeconfig_get_tags_return_value; *str != NULL; ++str)
    {
        ++count;
    }
    char ** result = cpl_calloc(count + 1, sizeof(char*));
    cpl_size n = 0;
    for (n = 0; n < count; ++n)
    {
        result[n] = cpl_strdup(cpl_recipeconfig_get_tags_return_value[n]);
    }
    return result;
}

/* Overload the following functions to return NULL. */

char** cpl_recipeconfig_get_inputs(const cpl_recipeconfig* self,
                                   const char* tag)
{
    (void) self;
    (void) tag;
    return NULL;
}

char** cpl_recipeconfig_get_outputs(const cpl_recipeconfig* self,
                                    const char* tag)
{
    (void) self;
    (void) tag;
    return NULL;
}

/* Overload the following function to produce an error. */
cpl_error_code
cpl_pluginlist_append(cpl_pluginlist *self, const cpl_plugin *plugin)
{
    (void) self;
    (void) plugin;
    return cpl_error_set_message("cpl_pluginlist_append",
                                 CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
}

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    test_parse_number();
    test_parameter_to_json();
    test_plugin_to_json();

#ifdef ENABLE_PYTHON_RECIPES
    test_json_to_pluginlist();
#endif

    /* Setting this to NULL helps in finding memory leaks with Valgrind,
       because the memory will not be marked as still reachable. */
    cpl_recipeconfig_get_tags_return_value = NULL;

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

static void test_parse_number(void)
{
    /* Test that the parse_number fuction produces an error code when parsing
       invalid input. We specifically check the two branches in parse_number for
       short and long input strings. */
    strtod_produce_error = CPL_TRUE;
    cpl_test_null(er_json_parse("123x4"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[123x4, \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
}

static void test_parameter_to_json(void)
{
    cpl_parameter * param1 = cpl_parameter_new_value("par1", CPL_TYPE_INT,
                                                    "desc","ctx", 10);
    cpl_parameter * param2 = cpl_parameter_new_range("par2", CPL_TYPE_INT,
                                                    "desc","ctx", 5, 1, 10);
    cpl_parameter * param3 = cpl_parameter_new_enum("par3", CPL_TYPE_INT,
                                                    "desc","ctx", 2,
                                                    3, 1, 2, 3);
    strtod_produce_error = CPL_FALSE;

    /* Test that the parameter_to_json fuction produces an error code when
       trying to handle a parameter that has the wrong class or type. */
    cpl_parameter_get_type_return_value = CPL_TYPE_INVALID;
    cpl_parameter_get_class_return_value = CPL_PARAMETER_CLASS_INVALID;
    cpl_parameter_get_enum_size_return_value = 3;
    cpl_test_null(parameter_to_json(param1));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    cpl_test_null(parameter_to_json(param2));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    cpl_test_null(parameter_to_json(param3));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);

    cpl_parameter_get_type_return_value = CPL_TYPE_INT;
    cpl_parameter_get_class_return_value = CPL_PARAMETER_CLASS_INVALID;
    cpl_parameter_get_enum_size_return_value = 3;
    cpl_test_null(parameter_to_json(param1));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    cpl_test_null(parameter_to_json(param2));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    cpl_test_null(parameter_to_json(param3));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);

    cpl_parameter_get_type_return_value = CPL_TYPE_CHAR;
    cpl_parameter_get_class_return_value = CPL_PARAMETER_CLASS_VALUE;
    cpl_test_null(parameter_to_json(param1));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    cpl_parameter_get_type_return_value = CPL_TYPE_BOOL;
    cpl_parameter_get_class_return_value = CPL_PARAMETER_CLASS_RANGE;
    cpl_test_null(parameter_to_json(param2));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    cpl_parameter_get_type_return_value = CPL_TYPE_BOOL;
    cpl_parameter_get_class_return_value = CPL_PARAMETER_CLASS_ENUM;
    cpl_parameter_get_enum_size_return_value = 3;
    cpl_test_null(parameter_to_json(param3));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);

    /* Test an edge case where an enumeration parameter indicates zero
       choices. */
    cpl_parameter_get_type_return_value = CPL_TYPE_INT;
    cpl_parameter_get_class_return_value = CPL_PARAMETER_CLASS_ENUM;
    cpl_parameter_get_enum_size_return_value = 0;
    cpl_test_null(parameter_to_json(param3));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_parameter_delete(param1);
    cpl_parameter_delete(param2);
    cpl_parameter_delete(param3);
}

static void test_plugin_to_json(void)
{
    cpl_plugin * plugin = NULL;
    cpl_recipe2 * recipe = NULL;
    char * json = NULL;
    er_json_node * parsetree = NULL;
    const er_json_node * node = NULL;
    const er_json_node * node2 = NULL;
    char* reference_result[2] = {"RAW", NULL};

    /* Must restore normal behaviour of strtod function which is used in the
       JSON parsing. */
    strtod_produce_error = CPL_FALSE;

    /* Prepare the plugin object to parse. */
    recipe = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));

    recipe->config = cpl_recipeconfig_new();
    cpl_recipeconfig_set_tag(recipe->config, "RAW", 2, 3);
    cpl_test_assert(recipe->config != NULL);

    plugin = &recipe->base.interface;
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE_V2, "test",
                                    "simple test", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);

    /* Test that NULL is returned if cpl_recipeconfig_get_tags returns a
       NULL value, which indicates an error. */
    cpl_recipeconfig_get_tags_return_value = NULL;
    cpl_test_null(er_plugin_to_json(plugin));
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);

    /* Test an edge case where the functions cpl_recipeconfig_get_inputs and
       cpl_recipeconfig_get_outputs return NULL. */
    cpl_recipeconfig_get_tags_return_value = reference_result;
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(json);

    /* Parse the produced JSON and make sure we get a valid result from the
       edge case. */
    parsetree = er_json_parse(json);
    cpl_test_nonnull(parsetree);
    cpl_test_error(CPL_ERROR_NONE);

    node = er_json_node_object_get(parsetree, "recipeconfig");
    cpl_test_nonnull(node);
    cpl_test_eq(er_json_node_array_size(node), 1);
    er_json_array_iterator iter = er_json_node_array_begin(node);
    node = er_json_node_array_get(node, iter);
    node2 = er_json_node_object_get(node, "tag");
    cpl_test_eq_string(er_json_node_get_string(node2), "RAW");
    node2 = er_json_node_object_get(node, "min");
    cpl_test_eq(er_json_node_get_number(node2), 2);
    node2 = er_json_node_object_get(node, "max");
    cpl_test_eq(er_json_node_get_number(node2), 3);

    cpl_free(json);
    cpl_recipeconfig_delete(recipe->config);
    cpl_plugin_delete(plugin);
    er_json_node_delete(parsetree);
}

#ifdef ENABLE_PYTHON_RECIPES

static void test_json_to_pluginlist(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;

    /* Test edge case where failure in cpl_pluginlist_append is handled
       appropriately by er_json_to_pluginlist. */
    json =
        "[\n"
        "  {\n"
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
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_pluginlist(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);
    er_json_node_delete(parsetree);
}

#endif /* ENABLE_PYTHON_RECIPES */
