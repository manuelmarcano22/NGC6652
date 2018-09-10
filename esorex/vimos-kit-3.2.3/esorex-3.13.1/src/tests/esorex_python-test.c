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
#include <sys/stat.h>
#include <cpl_test.h>
#include <cpl.h>
#include "er_python.h"

/*-----------------------------------------------------------------------------
                                 Defines
 -----------------------------------------------------------------------------*/

#define PYTHON_INPUT_FD   3
#define PYTHON_OUTPUT_FD  4

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void mock_python(const char * script);
static void test_error_from_interpreter(void);
static void test_load_modules_function(void);
static void test_module_set_function(void);
static void test_get_plugin_list(void);
static void test_plugin_execution(void);
static void test_plugin_execution_with_simple_name(void);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Update the PATH environment variable to make sure we pick up our mock
       python script. */
    const char * path = "./mock_python_test:/usr/bin:/bin";
    cpl_test_assert(setenv("PATH", path, 0x1) == 0);

    /* Insert tests below */
    test_error_from_interpreter();
    test_load_modules_function();
    test_module_set_function();
    test_get_plugin_list();
    test_plugin_execution();
    test_plugin_execution_with_simple_name();

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
    if (mkdir("./mock_python_test", 0777) != 0)
    {
        cpl_test_assert(errno == EEXIST);
    }

    FILE* file = fopen("./mock_python_test/python", "w");
    cpl_test_assert(file != NULL);
    size_t bytes_to_write = strlen(script);
    size_t bytes_written = fwrite(script, sizeof(char), bytes_to_write, file);
    cpl_test_assert(bytes_written == bytes_to_write);
    int fclose_result = fclose(file);
    cpl_test_assert(fclose_result == 0);
    mode_t mode = S_IRUSR | S_IWUSR | S_IXUSR |
                  S_IRGRP | S_IXGRP |
                  S_IROTH | S_IXOTH;
    cpl_test_assert(chmod("./mock_python_test/python", mode) == 0);
}

static void test_error_from_interpreter(void)
{
    /* Test handling of errors from the external mocked up Python process. */
    er_stringarray_t * modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "./test.py");

    mock_python("#!/bin/sh\n"
                "exit 1\n");

    cpl_test_eq_error(er_python_load_modules(modulelist), CPL_ERROR_FILE_IO);

    er_python_cleanup();
    er_stringarray_delete(modulelist);
}

static void test_load_modules_function(void)
{
    er_stringarray_t * modulelist = NULL;

    /* Check for success if an empty module list is given. */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{}\n"
        "EOF\n"
    );

    cpl_test_eq_error(er_python_load_modules(modulelist), CPL_ERROR_NONE);

    er_python_cleanup();
    er_stringarray_delete(modulelist);

    /* Check error handling if NULL is given. */
    cpl_test_eq_error(er_python_load_modules(NULL), CPL_ERROR_NULL_INPUT);

    /* The simplest case of python returning an empty JSON dictionary structure
       should succeed with no errors. */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "./test.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{}\n"
        "EOF\n"
    );

    cpl_test_eq_error(er_python_load_modules(modulelist), CPL_ERROR_NONE);
    cpl_test(! er_stringarray_present(modulelist, 0));

    er_python_cleanup();
    er_stringarray_delete(modulelist);

    /* Mock up a recipe by returning expected output from the mocked python
       interpreter and see that this succeeds. */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "./test.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"./test.py\": [\n"
        "    {\n"
        "      \"class\": \"testclass\",\n"
        "      \"name\": \"testrecipe\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin\",\n"
        "      \"description\": \"description ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"test.par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"id\": 0,\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"value\": 12,\n"
        "            \"default\": 5,\n"
        "            \"present\": false,\n"
        "            \"tag\": \"myparam\",\n"
        "            \"cli_enabled\": true,\n"
        "            \"cli_alias\": \"par1_cli\",\n"
        "            \"env_enabled\": false,\n"
        "            \"env_alias\": \"par1_env\",\n"
        "            \"cfg_enabled\": true,\n"
        "            \"cfg_alias\": \"par1_cfg\"\n"
        "         }\n"
        "      ]\n"
        "    }\n"
        "  ]\n"
        "}\n"
        "EOF\n"
    );

    cpl_test_eq_error(er_python_load_modules(modulelist), CPL_ERROR_NONE);
    cpl_test(er_stringarray_present(modulelist, 0));

    er_python_cleanup();
    er_stringarray_delete(modulelist);

    /* Check error handling if we got corrupted JSON output.  */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "./test.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{ ... , xyz\n"
        "EOF\n"
    );

    cpl_test_eq_error(er_python_load_modules(modulelist),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_test(er_stringarray_present(modulelist, 0));

    er_python_cleanup();
    er_stringarray_delete(modulelist);

    /* Check error handling if we did not get an array of plugins.  */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "./test.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"./test.py\": \"xyz\"\n"
        "}\n"
        "EOF\n"
    );

    cpl_test_eq_error(er_python_load_modules(modulelist),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_test(! er_stringarray_present(modulelist, 0));

    er_python_cleanup();
    er_stringarray_delete(modulelist);

    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "./test.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "\"./test.py\"\n"
        "EOF\n"
    );

    cpl_test_eq_error(er_python_load_modules(modulelist),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_test(er_stringarray_present(modulelist, 0));

    er_python_cleanup();
    er_stringarray_delete(modulelist);

    /* Check error handling if the plugin structure is invalid.  */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "./test.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"./test.py\": [\n"
        "    {\"dummy\": 345}\n"
        "  ]\n"
        "}\n"
        "EOF\n"
    );

    cpl_test_eq_error(er_python_load_modules(modulelist),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_test(! er_stringarray_present(modulelist, 0));

    er_python_cleanup();
    er_stringarray_delete(modulelist);

    /* Check error handling if a second invalid plugin structure is present. The
       er_python_load_modules function should still be able to identify the
       plugin that is OK.  */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "./test1.py");
    er_stringarray_append(modulelist, "./test2.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"./test1.py\": [\n"
        "    {\n"
        "      \"class\": \"testclass\",\n"
        "      \"name\": \"testrecipe\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin\",\n"
        "      \"description\": \"description ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"test.par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"id\": 0,\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"value\": 12,\n"
        "            \"default\": 5,\n"
        "            \"present\": false,\n"
        "            \"tag\": \"myparam\",\n"
        "            \"cli_enabled\": true,\n"
        "            \"cli_alias\": \"par1_cli\",\n"
        "            \"env_enabled\": false,\n"
        "            \"env_alias\": \"par1_env\",\n"
        "            \"cfg_enabled\": true,\n"
        "            \"cfg_alias\": \"par1_cfg\"\n"
        "         }\n"
        "      ]\n"
        "    }\n"
        "  ],\n"
        "  \"./test2.py\": [\n"
        "    {\"dummy\": 345}\n"
        "  ]\n"
        "}\n"
        "EOF\n"
    );

    cpl_test_eq_error(er_python_load_modules(modulelist),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_test(er_stringarray_present(modulelist, 0));
    cpl_test(! er_stringarray_present(modulelist, 1));

    er_python_cleanup();
    er_stringarray_delete(modulelist);
}

static void test_module_set_function(void)
{
    er_stringarray_t * modulelist = NULL;

    /* Perform basic tests of the Python module selection function. */
    cpl_test_eq_error(er_python_select_module(NULL), CPL_ERROR_NONE);
    cpl_test_eq_error(er_python_select_module("test.py"),
                      CPL_ERROR_DATA_NOT_FOUND);

    /* Test that the er_python_select_module function finds a module that has
       been loaded correctly.
      */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "test.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"test.py\": [\n"
        "    {\n"
        "      \"class\": \"testclass\",\n"
        "      \"name\": \"testrecipe\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin\",\n"
        "      \"description\": \"description ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": []\n"
        "    }\n"
        "  ]\n"
        "}\n"
        "EOF\n"
    );

    cpl_test_assert(er_python_load_modules(modulelist) == CPL_ERROR_NONE);
    er_stringarray_delete(modulelist);

    cpl_test_eq_error(er_python_select_module("test.py"), CPL_ERROR_NONE);

    /* Check that an unloaded module is still not found. */
    cpl_test_eq_error(er_python_select_module("other.py"),
                      CPL_ERROR_DATA_NOT_FOUND);

    er_python_cleanup();
}

static void test_get_plugin_list(void)
{
    er_stringarray_t * modulelist = NULL;
    cpl_pluginlist * list = NULL;
    int result = -1;

    /* Check error handling if NULL is given. */
    cpl_test(er_python_get_plugin_list(NULL) != 0);
    cpl_test_error(CPL_ERROR_NULL_INPUT);

    /* Check error handling if trying to get a plugin list if nothing was
       selected. */
    list = cpl_pluginlist_new();
    cpl_test_assert(list != NULL);
    cpl_test(er_python_get_plugin_list(list) != 0);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_pluginlist_delete(list);

    /* Test that a plugin list is returned successfully if a plugin was
       instantiated and selected correctly. */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    list = cpl_pluginlist_new();
    cpl_test_assert(list != NULL);
    er_stringarray_append(modulelist, "test.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"test.py\": [\n"
        "    {\n"
        "      \"class\": \"testclass\",\n"
        "      \"name\": \"testrecipe\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin\",\n"
        "      \"description\": \"description ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": []\n"
        "    }\n"
        "  ]\n"
        "}\n"
        "EOF\n"
    );

    cpl_test_assert(er_python_load_modules(modulelist) == CPL_ERROR_NONE);
    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    result = er_python_get_plugin_list(list);
    cpl_test_eq(result, 0);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_pluginlist_delete(list);
    er_stringarray_delete(modulelist);
    er_python_cleanup();
}

static void test_plugin_execution(void)
{
    er_stringarray_t * modulelist = NULL;
    cpl_pluginlist * list = NULL;
    int result = -1;
    cpl_recipe * recipe = NULL;
    cpl_recipe2 * recipe2 = NULL;
    cpl_plugin_func initfunc = NULL;
    cpl_plugin_func execfunc = NULL;
    cpl_plugin_func deinitfunc = NULL;
    cpl_frame * frame = NULL;

    /* Prepare a reference input frame. */
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

    /* First prepare a python mockup that will load a version 1 and 2 recipe. */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "recipe_v1.py");
    er_stringarray_append(modulelist, "recipe_v2.py");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"recipe_v1.py\": [\n"
        "    {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ]\n"
        "    }\n"
        "  ],\n"
        "  \"recipe_v2.py\": [\n"
        "    {\n"
        "      \"class\": \"testclass2\",\n"
        "      \"name\": \"testrecipe2\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v2\",\n"
        "      \"description\": \"description v2 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"range\",\n"
        "            \"description\": \"range param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 6,\n"
        "            \"min\": 1,\n"
        "            \"max\": 9\n"
        "         }\n"
        "      ],\n"
        "      \"recipeconfig\": [\n"
        "         {\n"
        "            \"tag\": \"RAW\",\n"
        "            \"min\": 1,\n"
        "            \"max\": 1,\n"
        "            \"inputs\": [\n"
        "               {\n"
        "                  \"tag\": \"CALIB\",\n"
        "                  \"min\": 1,\n"
        "                  \"max\": 1\n"
        "               }\n"
        "            ],\n"
        "            \"outputs\": [\"PROD\"]\n"
        "         }\n"
        "      ]\n"
        "    }\n"
        "  ]\n"
        "}\n"
        "EOF\n"
    );

    /* Load the mocked recipe plugins and make sure this succeeds. */
    cpl_test_assert(er_python_load_modules(modulelist) == CPL_ERROR_NONE);
    cpl_test(er_stringarray_present(modulelist, 0));
    cpl_test(er_stringarray_present(modulelist, 1));
    cpl_test_assert(er_python_select_module("recipe_v1.py") == CPL_ERROR_NONE);
    cpl_test_assert(er_python_select_module("recipe_v2.py") == CPL_ERROR_NONE);

    /* Now select recipe v1 for the following tests. */
    cpl_test_assert(er_python_select_module("recipe_v1.py") == CPL_ERROR_NONE);
    list = cpl_pluginlist_new();
    cpl_test_assert(list != NULL);
    result = er_python_get_plugin_list(list);
    cpl_test_eq(result, 0);
    cpl_test_error(CPL_ERROR_NONE);

    /* Prepare a python mocking script to simulate the response from a valid
       recipe execution and test the recipe invokation mechanics. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"input.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        },\n"
        "        {\n"
        "            \"filename\": \"product.fits\",\n"
        "            \"tag\": \"PROD\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 3,\n"
        "            \"level\": 3\n"
        "        }\n"
        "      ]\n"
        "  },\n"
        "  \"result\": 0,\n"
        "  \"error\": null\n"
        "}\n"
        "EOF\n"
    );

    /* Test that invoking the recipe works successfully. We have to also mock up
       how EsoRex will create a new recipe object. */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq(execfunc(&recipe->interface), 0);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* The following tests check error handling when invoking the plugin's init,
       exec and deinit functions. */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);

    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);

    /* Check for correct error handling when a NULL is given */
    cpl_test_eq_error(initfunc(NULL), CPL_ERROR_NULL_INPUT);
    cpl_test_eq_error(execfunc(NULL), CPL_ERROR_NULL_INPUT);
    cpl_test_eq_error(deinitfunc(NULL), CPL_ERROR_NULL_INPUT);

    /* Check for correct error handling if the wrong API version is given. */
    recipe->interface.api = 0;
    cpl_test_eq_error(initfunc(&recipe->interface), CPL_ERROR_UNSUPPORTED_MODE);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_UNSUPPORTED_MODE);
    cpl_test_eq_error(deinitfunc(&recipe->interface),
                      CPL_ERROR_UNSUPPORTED_MODE);
    recipe->interface.api = CPL_PLUGIN_API;

    /* Check for correct error handling if the plugin type is wrong. */
    recipe->interface.type = CPL_PLUGIN_TYPE_NONE;
    cpl_test_eq_error(initfunc(&recipe->interface), CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_eq_error(deinitfunc(&recipe->interface),
                      CPL_ERROR_UNSUPPORTED_MODE);
    recipe->interface.type = CPL_PLUGIN_TYPE_RECIPE;

    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Test the error handling if the recipe produces a truncated frameset. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "      ]\n"
        "  },\n"
        "  \"result\": 0,\n"
        "  \"error\": null\n"
        "}\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_ILLEGAL_INPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Test the error handling if the recipe does not produce a plugin
       structure in the output JSON. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"result\": 0,\n"
        "  \"error\": null\n"
        "}\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_ILLEGAL_INPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Test the error handling if the recipe does not produce a result field
       in the output JSON. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"input.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        },\n"
        "        {\n"
        "            \"filename\": \"product.fits\",\n"
        "            \"tag\": \"PROD\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 3,\n"
        "            \"level\": 3\n"
        "        }\n"
        "      ]\n"
        "  },\n"
        "  \"error\": null\n"
        "}\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_ILLEGAL_INPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Check for correct handling when the recipe indicates an error. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"input.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        }\n"
        "      ]\n"
        "  },\n"
        "  \"result\": 99,\n"
        "  \"error\": \"FAIL FROM RECIPE\"\n"
        "}\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq(execfunc(&recipe->interface), 99);
    cpl_test(strstr(cpl_error_get_message(), "FAIL FROM RECIPE") != NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Another test for error handling, but check the case when the error
       message is set to null. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"input.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        }\n"
        "      ]\n"
        "  },\n"
        "  \"result\": 999,\n"
        "  \"error\": null\n"
        "}\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq(execfunc(&recipe->interface), 999);
    cpl_test(strstr(cpl_error_get_message(), "FAIL FROM RECIPE") == NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Make sure error handling is as expected even if the error message string
       is not set by the recipe. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"input.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        }\n"
        "      ]\n"
        "  },\n"
        "  \"result\": 88\n"
        "}\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq(execfunc(&recipe->interface), 88);
    cpl_test_error(CPL_ERROR_ILLEGAL_OUTPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Check error handling if the 'result' field is not a number. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"input.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        }\n"
        "      ]\n"
        "  },\n"
        "  \"result\": null,\n"
        "  \"error\": \"FAIL FROM RECIPE\"\n"
        "}\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_ILLEGAL_INPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Test error handling when the recipe produces an invalid 'plugin'
       structure. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\"dummy\": 123},\n"
        "  \"result\": 0,\n"
        "  \"error\": null\n"
        "}\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_ILLEGAL_INPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Test error handling when the recipe produces the wrong format. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "[]\n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_ILLEGAL_INPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Test error handling when the recipe produces garbage. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{ ... % \" , \n"
        "EOF\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_ILLEGAL_INPUT);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Test the error handling when we simulate a complete failure of the
       python process. */
    mock_python(
        "#!/bin/sh\n"
        "exit 7\n"
    );

    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq_error(execfunc(&recipe->interface), CPL_ERROR_FILE_IO);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Now select recipe v2 and test the invocation mechanics. */
    cpl_test_assert(er_python_select_module("recipe_v2.py") == CPL_ERROR_NONE);
    cpl_pluginlist_delete(list);
    list = cpl_pluginlist_new();
    cpl_test_assert(list != NULL);
    result = er_python_get_plugin_list(list);
    cpl_test_eq(result, 0);
    cpl_test_error(CPL_ERROR_NONE);

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass2\",\n"
        "      \"name\": \"testrecipe2\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v2\",\n"
        "      \"description\": \"description v2 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"input.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        },\n"
        "        {\n"
        "            \"filename\": \"product.fits\",\n"
        "            \"tag\": \"PROD\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 3,\n"
        "            \"level\": 3\n"
        "        }\n"
        "      ],\n"
        "      \"recipeconfig\": [\n"
        "         {\n"
        "            \"tag\": \"RAW\",\n"
        "            \"min\": 1,\n"
        "            \"max\": 1,\n"
        "            \"inputs\": [\n"
        "               {\n"
        "                  \"tag\": \"CALIB\",\n"
        "                  \"min\": 1,\n"
        "                  \"max\": 1\n"
        "               }\n"
        "            ],\n"
        "            \"outputs\": [\"PROD\"]\n"
        "         }\n"
        "      ]\n"
        "  },\n"
        "  \"result\": 0,\n"
        "  \"error\": null\n"
        "}\n"
        "EOF\n"
    );

    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    cpl_test_eq_error(cpl_plugin_copy(&recipe2->base.interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe2->base.frames = cpl_frameset_new();
    cpl_test_assert(recipe2->base.frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe2->base.frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    cpl_test_error(CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe2->base.interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe2->base.interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe2->base.interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq(execfunc(&recipe2->base.interface), 0);
    deinitfunc = cpl_plugin_get_deinit(&recipe2->base.interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe2->base.interface), 0);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_plugin_delete(&recipe2->base.interface);

    /* Cleanup common objects for the recipe execution tests. */
    cpl_pluginlist_delete(list);
    er_stringarray_delete(modulelist);
    er_python_cleanup();
    cpl_frame_delete(frame);
}

static void test_plugin_execution_with_simple_name(void)
{
    /* This is a test for an edge case where the python module name used does
       not end in a .py extension. */
    er_stringarray_t * modulelist = NULL;
    cpl_pluginlist * list = NULL;
    int result = -1;
    cpl_recipe * recipe = NULL;
    cpl_plugin_func initfunc = NULL;
    cpl_plugin_func execfunc = NULL;
    cpl_plugin_func deinitfunc = NULL;
    cpl_frame * frame = NULL;

    /* Prepare a reference input frame. */
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

    /* First prepare a python mockup that will simulate loading of a recipe
       plugin. */
    modulelist = er_stringarray_new();
    cpl_test_assert(modulelist != NULL);
    er_stringarray_append(modulelist, "testrecipe");

    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"testrecipe\": [\n"
        "    {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ]\n"
        "    }\n"
        "  ]\n"
        "}\n"
        "EOF\n"
    );

    /* Load the mocked recipe plugin and make sure this succeeds. */
    cpl_test_assert(er_python_load_modules(modulelist) == CPL_ERROR_NONE);
    cpl_test(er_stringarray_present(modulelist, 0));
    cpl_test_assert(er_python_select_module("testrecipe") == CPL_ERROR_NONE);
    list = cpl_pluginlist_new();
    cpl_test_assert(list != NULL);
    result = er_python_get_plugin_list(list);
    cpl_test_eq(result, 0);
    cpl_test_error(CPL_ERROR_NONE);

    /* Prepare a python mocking script to simulate the recipe response. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)" > /dev/null\n"
        "cat <<EOF >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
        "{\n"
        "  \"plugin\": {\n"
        "      \"class\": \"testclass1\",\n"
        "      \"name\": \"testrecipe1\",\n"
        "      \"version\": 123,\n"
        "      \"synopsis\": \"test plugin v1\",\n"
        "      \"description\": \"description v1 ...\",\n"
        "      \"author\": \"Some Author\",\n"
        "      \"email\": \"someone@eso.org\",\n"
        "      \"copyright\": \"copyright ...\",\n"
        "      \"parameters\": [\n"
        "         {\n"
        "            \"name\": \"par1\",\n"
        "            \"class\": \"value\",\n"
        "            \"description\": \"int param\",\n"
        "            \"context\": \"test\",\n"
        "            \"default\": 3\n"
        "         }\n"
        "      ],\n"
        "      \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"input.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        },\n"
        "        {\n"
        "            \"filename\": \"product.fits\",\n"
        "            \"tag\": \"PROD\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 3,\n"
        "            \"level\": 3\n"
        "        }\n"
        "      ]\n"
        "  },\n"
        "  \"result\": 0,\n"
        "  \"error\": null\n"
        "}\n"
        "EOF\n"
    );

    /* Prepare the recipe plugin as EsoRex would and invoke the recipe. */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    cpl_test_eq_error(cpl_plugin_copy(&recipe->interface,
                                      cpl_pluginlist_get_first(list)),
                      CPL_ERROR_NONE);
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    cpl_test_assert(cpl_frameset_insert(recipe->frames,
                                        cpl_frame_duplicate(frame))
                    == CPL_ERROR_NONE);
    initfunc = cpl_plugin_get_init(&recipe->interface);
    cpl_test_assert(initfunc != NULL);
    cpl_test_eq(initfunc(&recipe->interface), 0);
    execfunc = cpl_plugin_get_exec(&recipe->interface);
    cpl_test_assert(execfunc != NULL);
    cpl_test_eq(execfunc(&recipe->interface), 0);
    deinitfunc = cpl_plugin_get_deinit(&recipe->interface);
    cpl_test_assert(deinitfunc != NULL);
    cpl_test_eq(deinitfunc(&recipe->interface), 0);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);

    /* Cleanup common objects. */
    cpl_pluginlist_delete(list);
    er_stringarray_delete(modulelist);
    er_python_cleanup();
    cpl_frame_delete(frame);
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
