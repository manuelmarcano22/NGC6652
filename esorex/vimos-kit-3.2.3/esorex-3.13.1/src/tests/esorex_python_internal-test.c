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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef ENABLE_PYTHON_RECIPES

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <cpl_test.h>
#include "er_python.c"

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void mock_python(const char * script);
static void test_plugin_initialize(void);
static void test_plugin_execute(void);
static void test_update_plugin(void);
static void test_plugin_compare(void);
static void test_run_python_command(void);
static void test_stop_python_interpreter(void);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Update the PATH environment variable to make sure we pick up our mock
       python script. */
    const char * path = "./mock_python_internal_test:/usr/bin:/bin";
    cpl_test_assert(setenv("PATH", path, 0x1) == 0);

    /* Insert tests below */
    test_plugin_initialize();
    test_plugin_execute();
    test_update_plugin();
    test_plugin_compare();
    test_run_python_command();
    test_stop_python_interpreter();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

/*
 * This function is used to create a shell script that will mock up a Python
 * interpreter process. The input is the shell script to write to the python
 * script file.
 */
static void mock_python(const char * script)
{
    /* Do not overwrite the python script if a previous test failed. This will
       make it easier to debug by having the failing script available. */
    cpl_test_assert(cpl_test_get_failed() == 0);

    /* Create the directory for the python script. We put this in its own
       directory so that this unit test can be made concurrency safe. */
    if (mkdir("./mock_python_internal_test", 0777) != 0)
    {
        cpl_test_assert(errno == EEXIST);
    }

    FILE* file = fopen("./mock_python_internal_test/python", "w");
    cpl_test_assert(file != NULL);
    size_t bytes_to_write = strlen(script);
    size_t bytes_written = fwrite(script, sizeof(char), bytes_to_write, file);
    cpl_test_assert(bytes_written == bytes_to_write);
    int fclose_result = fclose(file);
    cpl_test_assert(fclose_result == 0);
    mode_t mode = S_IRUSR | S_IWUSR | S_IXUSR |
                  S_IRGRP | S_IXGRP |
                  S_IROTH | S_IXOTH;
    cpl_test_assert(chmod("./mock_python_internal_test/python", mode) == 0);
}

static void test_plugin_initialize(void)
{
    char * module_name = NULL;
    cpl_recipe * recipe = NULL;
    cpl_recipe * recipe2 = NULL;
    cpl_recipe2 * recipe_v2 = NULL;
    cpl_recipe2 * recipe2_v2 = NULL;
    cpl_pluginlist * pluginlist = NULL;

    /* Test an edge case where plugin_initialize is given a plugin type that it
       cannot handle.
       We first have to prepare the module lookup tables as would be existing in
       the normal case just before a call to the function under test. Only then
       can we make the invocation to perform the actual test.
      */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_init(&recipe->interface, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_NONE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    recipe2 = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_copy(&recipe2->interface, &recipe->interface)
                    == CPL_ERROR_NONE);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, &recipe->interface)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, recipe, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    /* Invoke the function we are trying to test. */
    cpl_test_eq_error(plugin_initialize(&recipe2->interface),
                      CPL_ERROR_UNSUPPORTED_MODE);

    cpl_parameterlist_delete(recipe2->parameters);
    cpl_frameset_delete(recipe2->frames);
    cpl_plugin_delete(&recipe2->interface);
    er_python_cleanup();

    /* Test an edge case where plugin_initialize is given a plugin that has no
       parameters structure set.
      */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_init(&recipe->interface, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    recipe2 = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_copy(&recipe2->interface, &recipe->interface)
                    == CPL_ERROR_NONE);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, &recipe->interface)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, recipe, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    cpl_test_eq_error(plugin_initialize(&recipe2->interface), CPL_ERROR_NONE);

    cpl_parameterlist_delete(recipe2->parameters);
    cpl_frameset_delete(recipe2->frames);
    cpl_plugin_delete(&recipe2->interface);
    er_python_cleanup();

    /* Check an edge case when a recipe v2 plugin is given a plugin with an
       empty recipe config.
      */
    recipe_v2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_init(&recipe_v2->base.interface, CPL_PLUGIN_API,
                                    1, CPL_PLUGIN_TYPE_RECIPE_V2, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    recipe2_v2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_assert(cpl_plugin_copy(&recipe2_v2->base.interface,
                                    &recipe_v2->base.interface)
                    == CPL_ERROR_NONE);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist,
                                          &recipe_v2->base.interface)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, recipe_v2, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    cpl_test_eq_error(plugin_initialize(&recipe2_v2->base.interface),
                      CPL_ERROR_NONE);

    cpl_recipeconfig_delete(recipe2_v2->config);
    cpl_parameterlist_delete(recipe2_v2->base.parameters);
    cpl_frameset_delete(recipe2_v2->base.frames);
    cpl_plugin_delete(&recipe2_v2->base.interface);
    er_python_cleanup();
}

static void test_plugin_execute(void)
{
    char * module_name = NULL;
    cpl_recipe * recipe = NULL;
    cpl_recipe * recipe2 = NULL;
    cpl_pluginlist * pluginlist = NULL;
    cpl_parameter * param = NULL;

    /* Test an edge case where the plugin_execute function is given a plugin
       type that they cannot handle.
       We first have to prepare the module lookup tables as would be existing in
       the normal case just before a call to the functions under test. Only then
       can we make the invocations to perform the actual test.
      */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_init(&recipe->interface, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_NONE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    recipe2 = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_copy(&recipe2->interface, &recipe->interface)
                    == CPL_ERROR_NONE);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, &recipe->interface)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, recipe, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    /* Invoke the function we are trying to test. */
    cpl_test_eq_error(plugin_execute(&recipe2->interface),
                      CPL_ERROR_UNSUPPORTED_MODE);

    cpl_parameterlist_delete(recipe2->parameters);
    cpl_frameset_delete(recipe2->frames);
    cpl_plugin_delete(&recipe2->interface);
    er_python_cleanup();

    /* Test an edge case where the plugin_execute function is given a plugin
       that has no parameters structure set.
      */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_init(&recipe->interface, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    recipe2 = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_copy(&recipe2->interface, &recipe->interface)
                    == CPL_ERROR_NONE);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, &recipe->interface)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, recipe, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    cpl_test_eq_error(plugin_execute(&recipe2->interface),
                      CPL_ERROR_DATA_NOT_FOUND);

    cpl_parameterlist_delete(recipe2->parameters);
    cpl_frameset_delete(recipe2->frames);
    cpl_plugin_delete(&recipe2->interface);
    er_python_cleanup();

    /* Test error handling of plugin_execute if the __class__ parameter in the
       simulate plugin is missing. */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_init(&recipe->interface, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    recipe->parameters = cpl_parameterlist_new();

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, &recipe->interface)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, recipe, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    recipe2 = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_copy(&recipe2->interface,
                                    &recipe->interface)
                    == CPL_ERROR_NONE);

    cpl_test_eq_error(plugin_execute(&recipe2->interface),
                      CPL_ERROR_DATA_NOT_FOUND);

    cpl_parameterlist_delete(recipe2->parameters);
    cpl_frameset_delete(recipe2->frames);
    cpl_plugin_delete(&recipe2->interface);
    cpl_parameterlist_delete(recipe->parameters);
    cpl_frameset_delete(recipe->frames);
    er_python_cleanup();

    /* Test error handling of plugin_execute if the __class__ parameter in the
       simulate plugin is wrong. */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_init(&recipe->interface, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);

    recipe->parameters = cpl_parameterlist_new();
    param = cpl_parameter_new_value("__class__", CPL_TYPE_INT,
                                    "Python class name", "__python__", 123);
    cpl_test_assert(param != NULL);
    cpl_test_assert(cpl_parameterlist_append(recipe->parameters, param)
                    == CPL_ERROR_NONE);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, &recipe->interface)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, recipe, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    recipe2 = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_copy(&recipe2->interface,
                                    &recipe->interface)
                    == CPL_ERROR_NONE);

    cpl_test_eq_error(plugin_initialize(&recipe2->interface), CPL_ERROR_NONE);
    cpl_test_eq_error(plugin_execute(&recipe2->interface),
                      CPL_ERROR_TYPE_MISMATCH);

    cpl_parameterlist_delete(recipe2->parameters);
    cpl_frameset_delete(recipe2->frames);
    cpl_plugin_delete(&recipe2->interface);
    cpl_parameterlist_delete(recipe->parameters);
    cpl_frameset_delete(recipe->frames);
    er_python_cleanup();
}

static void test_update_plugin(void)
{
    cpl_recipe * recipe = NULL;
    cpl_recipe * recipe2 = NULL;
    cpl_frame * frame = NULL;

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_init(&recipe->interface, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    recipe2 = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_assert(cpl_plugin_copy(&recipe2->interface, &recipe->interface)
                    == CPL_ERROR_NONE);

    /* Test error handling if update_plugin is given a plugin with the wrong API
       version. */
    recipe2->interface.api = 0;
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_TYPE_MISMATCH);
    recipe2->interface.api = recipe->interface.api;   /* restore value. */

    /* Test error handling if update_plugin is given a plugin with a version
       number that does not match. */
    recipe2->interface.version = 0;
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_TYPE_MISMATCH);
    recipe2->interface.version = recipe->interface.version; /* restore value. */

    /* Test error handling if update_plugin is given a plugin with a type that
       does not match. */
    recipe2->interface.type = CPL_PLUGIN_TYPE_RECIPE_V2;
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_TYPE_MISMATCH);
    recipe2->interface.type = recipe->interface.type;   /* restore value. */

    /* Test error handling if update_plugin is given a plugin with a name that
       does not match. */
    const char * oldname = recipe2->interface.name;
    recipe2->interface.name = "something";
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_TYPE_MISMATCH);
    recipe2->interface.name = oldname;   /* restore value. */

    /* Test error handling if update_plugin is given a plugins with a type that
       it cannot handle. */
    recipe->interface.type = CPL_PLUGIN_TYPE_NONE;
    recipe2->interface.type = CPL_PLUGIN_TYPE_NONE;
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_UNSUPPORTED_MODE);
    recipe->interface.type = CPL_PLUGIN_TYPE_RECIPE;
    recipe2->interface.type = CPL_PLUGIN_TYPE_RECIPE;

    /* Test error handling if update_plugin is given plugin objects that do not
       have an allocated frameset. */
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_ILLEGAL_INPUT);

    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_frameset_delete(recipe->frames);
    recipe->frames = NULL;

    recipe2->frames = cpl_frameset_new();
    cpl_test_assert(recipe2->frames != NULL);
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_frameset_delete(recipe2->frames);
    recipe2->frames = NULL;

    /* Test edge case where frame values have been modified. */
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    frame = cpl_frame_new();
    cpl_test_assert(frame != NULL);
    cpl_test_assert(cpl_frame_set_filename(frame, "input.fits")
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_tag(frame, "RAW") == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_type(frame, CPL_FRAME_TYPE_NONE)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_group(frame, CPL_FRAME_GROUP_NONE)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_level(frame, CPL_FRAME_LEVEL_NONE)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frameset_insert(recipe->frames, frame)
                    == CPL_ERROR_NONE);

    recipe2->frames = cpl_frameset_new();
    cpl_test_assert(recipe2->frames != NULL);
    frame = cpl_frame_new();
    cpl_test_assert(frame != NULL);
    cpl_test_assert(cpl_frame_set_filename(frame, "input_moved.fits")
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_tag(frame, "RAW2") == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_type(frame, CPL_FRAME_TYPE_IMAGE)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_group(frame, CPL_FRAME_GROUP_RAW)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_level(frame, CPL_FRAME_LEVEL_FINAL)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frameset_insert(recipe2->frames, frame)
                    == CPL_ERROR_NONE);

    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_NONE);

    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);
    cpl_frameset_delete(recipe2->frames);
    cpl_plugin_delete(&recipe2->interface);
}

static void test_plugin_compare(void)
{
    /* Test that the comparison function returns the correct results. */
    cpl_plugin a, b;
    a.api = 1;
    a.version = 1;
    a.type = 1;
    a.name = "a";
    a.author = "a";
    a.email = "a";
    a.synopsis = "a";
    a.description = "a";
    a.copyright = "a";
    b.api = 2;
    b.version = 2;
    b.type = 2;
    b.name = "b";
    b.author = "b";
    b.email = "b";
    b.synopsis = "b";
    b.description = "b";
    b.copyright = "b";
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.api = 3;
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.api = b.api;
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.version = 3;
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.version = b.version;
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.type = 3;
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.type = b.type;
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.name = "c";
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.name = b.name;
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.author = "c";
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.author = b.author;
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.email = "c";
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.email = b.email;
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.synopsis = "c";
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.synopsis = b.synopsis;
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.description = "c";
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.description = b.description;
    cpl_test_eq(plugin_compare(&a, &b), -1);

    a.copyright = "c";
    cpl_test_eq(plugin_compare(&a, &b), 1);

    a.copyright = b.copyright;
    cpl_test_eq(plugin_compare(&a, &b), 0);
}

static void test_run_python_command(void)
{
    char * output = NULL;

    /* Mockup the python command that will just echo the input to output. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)
            " >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
    );

    /* Test for success of run_python_command when given an empty input and the
       simple echoing python mockup script. */
    output = run_python_command("x=1", "");
    cpl_test_eq_string(output, "");
    cpl_test_error(CPL_ERROR_NONE);
    free(output);

    /* Test that the input propagates to the output for the given mockup. */
    output = run_python_command("x=1", "some data");
    cpl_test_eq_string(output, "some data");
    cpl_test_error(CPL_ERROR_NONE);
    free(output);

    /* Test error handling if the python mockup exits with an error code. */
    mock_python(
        "#!/bin/sh\n"
        "exit 99\n"
    );
    cpl_test_null(run_python_command("x=1", "some data"));
    cpl_test_error(CPL_ERROR_FILE_IO);
}

static void test_stop_python_interpreter(void)
{
    /* Test error handling for edge case where we try stop the python
       interpreter when the child process was never started. */
    cpl_test_eq_error(stop_python_interpreter(1, -1, -1), CPL_ERROR_FILE_IO);

    /* Test that if a simulated python process hangs for a long time then the
       stop_python_interpreter function will send a terminate signal to bring
       the process down. i.e. avoid waiting for the process forever.
       An error code will be set since the process did not terminate normally.
     */
    mock_python(
        "#!/bin/sh\n"
        "sleep 601\n"
    );

    int start_time = time(NULL);

    pid_t pid;
    int input, output;
    cpl_test_assert(start_python_interpreter("x=1", &pid, &input, &output)
                    == CPL_ERROR_NONE);
    cpl_test_eq_error(stop_python_interpreter(pid, input, output),
                      CPL_ERROR_FILE_IO);

    int end_time = time(NULL);
    cpl_test(end_time - start_time < 600);
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
