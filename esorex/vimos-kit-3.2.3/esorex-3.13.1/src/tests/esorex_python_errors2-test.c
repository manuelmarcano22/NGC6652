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
 * various failure modes that are difficult to trigger in any other way. These
 * tests are not compatible with other tests from esorex_python_errors*-test.c
 * because commonly use functions are being overloaded. Thus, they must be
 * compiled into their own binary.
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

#include <cpl.h>
#include <cpl_test.h>

/*-----------------------------------------------------------------------------
                             Overload functions
 -----------------------------------------------------------------------------*/

/* Overload the following functions to allow producing errors on demand. */

cpl_boolean cpl_parameterlist_append_must_produce_error = CPL_FALSE;

cpl_error_code _cpl_parameterlist_append(cpl_parameterlist *self,
                                         cpl_parameter *parameter)
{
    if (cpl_parameterlist_append_must_produce_error)
    {
        return cpl_error_set_message("cpl_parameterlist_append",
                                     CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
    }
    else
    {
        return cpl_parameterlist_append(self, parameter);
    }
}

cpl_boolean cpl_recipeconfig_get_tags_must_produce_error = CPL_FALSE;

char** _cpl_recipeconfig_get_tags(const cpl_recipeconfig* self)
{
    if (cpl_recipeconfig_get_tags_must_produce_error)
    {
        cpl_error_set_message("cpl_recipeconfig_get_tags",
                              CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
        return NULL;
    }
    else
    {
        return cpl_recipeconfig_get_tags(self);
    }
}

cpl_boolean cpl_recipeconfig_get_inputs_must_produce_error = CPL_FALSE;

char** _cpl_recipeconfig_get_inputs(const cpl_recipeconfig* self,
                                    const char* tag)
{
    if (cpl_recipeconfig_get_inputs_must_produce_error)
    {
        cpl_error_set_message("cpl_recipeconfig_get_inputs",
                              CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
        return NULL;
    }
    else
    {
        return cpl_recipeconfig_get_inputs(self, tag);
    }
}

cpl_boolean cpl_recipeconfig_get_outputs_must_produce_error = CPL_FALSE;

char** _cpl_recipeconfig_get_outputs(const cpl_recipeconfig* self,
                                     const char* tag)
{
    if (cpl_recipeconfig_get_outputs_must_produce_error)
    {
        cpl_error_set_message("cpl_recipeconfig_get_outputs",
                              CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
        return NULL;
    }
    else
    {
        return cpl_recipeconfig_get_outputs(self, tag);
    }
}

cpl_boolean cpl_recipeconfig_set_tag_must_produce_error = CPL_FALSE;

int _cpl_recipeconfig_set_tag(cpl_recipeconfig* self, const char* tag,
                             cpl_size min_count, cpl_size max_count)
{
    if (cpl_recipeconfig_set_tag_must_produce_error)
    {
        return cpl_error_set_message("cpl_recipeconfig_set_tag",
                                     CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
    }
    else
    {
        return cpl_recipeconfig_set_tag(self, tag, min_count, max_count);
    }
}

cpl_boolean cpl_recipeconfig_set_input_must_produce_error = CPL_FALSE;

int _cpl_recipeconfig_set_input(cpl_recipeconfig* self, const char* tag,
                                const char* input, cpl_size min_count,
                                cpl_size max_count)
{
    if (cpl_recipeconfig_set_input_must_produce_error)
    {
        return cpl_error_set_message("cpl_recipeconfig_set_input",
                                     CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
    }
    else
    {
        return cpl_recipeconfig_set_input(self, tag, input, min_count,
                                          max_count);
    }
}

cpl_boolean cpl_recipeconfig_set_output_must_produce_error = CPL_FALSE;

int _cpl_recipeconfig_set_output(cpl_recipeconfig* self, const char* tag,
                                 const char* output)
{
    if (cpl_recipeconfig_set_output_must_produce_error)
    {
        return cpl_error_set_message("cpl_recipeconfig_set_output",
                                     CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
    }
    else
    {
        return cpl_recipeconfig_set_output(self, tag, output);
    }
}

/* Replace calls in the code we want to test with our overloaded functions. */
#define cpl_parameterlist_append _cpl_parameterlist_append
#define cpl_recipeconfig_get_tags _cpl_recipeconfig_get_tags
#define cpl_recipeconfig_get_inputs _cpl_recipeconfig_get_inputs
#define cpl_recipeconfig_get_outputs _cpl_recipeconfig_get_outputs
#define cpl_recipeconfig_set_tag _cpl_recipeconfig_set_tag
#define cpl_recipeconfig_set_input _cpl_recipeconfig_set_input
#define cpl_recipeconfig_set_output _cpl_recipeconfig_set_output
#include "er_python.c"
#undef cpl_parameterlist_append
#undef cpl_recipeconfig_get_tags
#undef cpl_recipeconfig_get_inputs
#undef cpl_recipeconfig_get_outputs
#undef cpl_recipeconfig_set_tag
#undef cpl_recipeconfig_set_input
#undef cpl_recipeconfig_set_output

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void test_plugin_initialize(void);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    test_plugin_initialize();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

static void test_plugin_initialize(void)
{
    char * module_name = NULL;
    cpl_recipe2 * recipe = NULL;
    cpl_recipe2 * recipe2 = NULL;
    cpl_pluginlist * pluginlist = NULL;
    cpl_parameter * param = NULL;

    /* Loading a simulated plugin with one parameter and a simple recipe config
       structure. */
    recipe = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_init(&recipe->base.interface, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE_V2, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);

    recipe->base.parameters = cpl_parameterlist_new();
    cpl_test_assert(recipe->base.parameters != NULL);
    param = cpl_parameter_new_value("test.par1", CPL_TYPE_INT, "int param",
                                    "test", 5);
    cpl_test_assert(param != NULL);
    cpl_test_assert(cpl_parameterlist_append(recipe->base.parameters, param)
                    == CPL_ERROR_NONE);
    param = cpl_parameter_new_value("__class__", CPL_TYPE_STRING,
                                    "Python class name", "__python__",
                                    "recipeclass");
    cpl_test_assert(param != NULL);
    cpl_test_assert(cpl_parameterlist_append(recipe->base.parameters, param)
                    == CPL_ERROR_NONE);

    recipe->config = cpl_recipeconfig_new();
    cpl_test_assert(recipe->config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(recipe->config, "RAW", 1, 1) == 0);
    cpl_test_assert(cpl_recipeconfig_set_input(recipe->config, "RAW",
                                               "CALIB", 1, 1)
                    == 0);
    cpl_test_assert(cpl_recipeconfig_set_output(recipe->config, "RAW", "PROD")
                    == 0);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, &recipe->base.interface)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, recipe, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    /* Test error handling of plugin_initialize if cpl_parameterlist_append
       produces and error. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_copy(&recipe2->base.interface,
                                    &recipe->base.interface)
                    == CPL_ERROR_NONE);
    cpl_parameterlist_append_must_produce_error = CPL_TRUE;
    cpl_recipeconfig_get_tags_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_inputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_outputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_tag_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_input_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_output_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(plugin_initialize(&recipe2->base.interface),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(&recipe2->base.interface);

    /* Test error handling of plugin_initialize if cpl_recipeconfig_get_tags
       produces and error. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_copy(&recipe2->base.interface,
                                    &recipe->base.interface)
                    == CPL_ERROR_NONE);
    cpl_parameterlist_append_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_tags_must_produce_error = CPL_TRUE;
    cpl_recipeconfig_get_inputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_outputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_tag_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_input_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_output_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(plugin_initialize(&recipe2->base.interface),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(&recipe2->base.interface);

    /* Test error handling of plugin_initialize if cpl_recipeconfig_get_inputs
       produces and error. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_copy(&recipe2->base.interface,
                                    &recipe->base.interface)
                    == CPL_ERROR_NONE);
    cpl_parameterlist_append_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_tags_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_inputs_must_produce_error = CPL_TRUE;
    cpl_recipeconfig_get_outputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_tag_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_input_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_output_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(plugin_initialize(&recipe2->base.interface),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(&recipe2->base.interface);

    /* Test error handling of plugin_initialize if cpl_recipeconfig_get_outputs
       produces and error. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_copy(&recipe2->base.interface,
                                    &recipe->base.interface)
                    == CPL_ERROR_NONE);
    cpl_parameterlist_append_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_tags_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_inputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_outputs_must_produce_error = CPL_TRUE;
    cpl_recipeconfig_set_tag_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_input_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_output_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(plugin_initialize(&recipe2->base.interface),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(&recipe2->base.interface);

    /* Test error handling of plugin_initialize if cpl_recipeconfig_set_tag
       produces and error. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_copy(&recipe2->base.interface,
                                    &recipe->base.interface)
                    == CPL_ERROR_NONE);
    cpl_parameterlist_append_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_tags_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_inputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_outputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_tag_must_produce_error = CPL_TRUE;
    cpl_recipeconfig_set_input_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_output_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(plugin_initialize(&recipe2->base.interface),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(&recipe2->base.interface);

    /* Test error handling of plugin_initialize if cpl_recipeconfig_set_input
       produces and error. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_copy(&recipe2->base.interface,
                                    &recipe->base.interface)
                    == CPL_ERROR_NONE);
    cpl_parameterlist_append_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_tags_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_inputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_outputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_tag_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_input_must_produce_error = CPL_TRUE;
    cpl_recipeconfig_set_output_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(plugin_initialize(&recipe2->base.interface),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(&recipe2->base.interface);

    /* Test error handling of plugin_initialize if cpl_recipeconfig_set_output
       produces and error. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_copy(&recipe2->base.interface,
                                    &recipe->base.interface)
                    == CPL_ERROR_NONE);
    cpl_parameterlist_append_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_tags_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_inputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_get_outputs_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_tag_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_input_must_produce_error = CPL_FALSE;
    cpl_recipeconfig_set_output_must_produce_error = CPL_TRUE;
    cpl_test_eq_error(plugin_initialize(&recipe2->base.interface),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(&recipe2->base.interface);

    /* Cleanup simulated plugin object. */
    cpl_parameterlist_delete(recipe->base.parameters);
    cpl_frameset_delete(recipe->base.frames);
    cpl_recipeconfig_delete(recipe->config);
    er_python_cleanup();
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
