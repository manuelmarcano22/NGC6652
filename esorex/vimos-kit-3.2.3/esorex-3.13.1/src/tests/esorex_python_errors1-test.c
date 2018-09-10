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
 * These are tests relying on overloading certain functions to mock various
 * failure modes that are difficult to trigger in any other way.
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

#include <errno.h>
#include <cpl.h>
#include <cpl_test.h>

/*-----------------------------------------------------------------------------
                            Overloaded functions
 -----------------------------------------------------------------------------*/

/* Overload the following function to produce a desired error. */

int errno_value_to_use = EIO;

int close(int fd)
{
    (void) fd;
    errno = errno_value_to_use;
    return -1;
}

/* Overload cpl_plugin_set_init to produce an error on demand. */

cpl_boolean cpl_plugin_set_init_must_produce_error = CPL_FALSE;

cpl_error_code _cpl_plugin_set_init(cpl_plugin * self, cpl_plugin_func func)
{
    if (cpl_plugin_set_init_must_produce_error)
    {
        return cpl_error_set_message("cpl_plugin_set_init",
                                     CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
    }
    else
    {
        return cpl_plugin_set_init(self, func);
    }
}

/* Overload cpl_pluginlist_append to either produce an error or call the
   original CPL function. */

cpl_boolean cpl_pluginlist_append_must_produce_error = CPL_FALSE;

cpl_error_code _cpl_pluginlist_append(cpl_pluginlist * self,
                                     const cpl_plugin * plugin)
{
    if (cpl_pluginlist_append_must_produce_error)
    {
        return cpl_error_set_message("cpl_pluginlist_append",
                                     CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
    }
    else
    {
        return cpl_pluginlist_append(self, plugin);
    }
}

cpl_boolean cpl_frameset_insert_must_produce_error = CPL_FALSE;

cpl_error_code _cpl_frameset_insert(cpl_frameset *self, cpl_frame *frame)
{
    if (cpl_pluginlist_append_must_produce_error)
    {
        return cpl_error_set_message("cpl_frameset_insert",
                                     CPL_ERROR_ILLEGAL_OUTPUT, "dummy error");
    }
    else
    {
        return cpl_frameset_insert(self, frame);
    }
}

/* Replace calls in the code we want to test with our overloaded functions. */
#define cpl_plugin_set_init _cpl_plugin_set_init
#define cpl_pluginlist_append _cpl_pluginlist_append
#define cpl_frameset_insert _cpl_frameset_insert
#include "er_python.c"
#undef cpl_plugin_set_init
#undef cpl_pluginlist_append
#undef cpl_frameset_insert

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void test_close_file(void);
static void test_get_plugin_list(void);
static void test_update_plugin(void);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    test_close_file();
    test_get_plugin_list();
    test_update_plugin();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

static void test_close_file(void)
{
    /* Test edge case where close() fails with EIO. This should be ignored. */
    errno_value_to_use = EIO;
    close_file(123);
    cpl_test_error(CPL_ERROR_NONE);

    /* Test edge case where close() fails with EINTR. This should be ignored. */
    errno_value_to_use = EINTR;
    close_file(123);
    cpl_test_error(CPL_ERROR_NONE);
}

static void test_get_plugin_list(void)
{
    char * module_name = NULL;
    cpl_plugin * plugin = NULL;
    cpl_pluginlist * pluginlist = NULL;
    cpl_pluginlist * list = NULL;

    /* Test and edge case where the cpl_plugin_set_init function fails when
       invoking er_python_get_plugin_list. We first have to prepare the module
       lookup tables as would be existing in the normal case just before a call
       to the er_python_get_plugin_list function. Only then can we set the
       cpl_plugin_set_init function to produce an error and make the call to
       er_python_get_plugin_list which is being tested.
      */
    plugin = cpl_plugin_new();
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, plugin)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, plugin, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    list = cpl_pluginlist_new();
    cpl_test_assert(list != NULL);

    cpl_pluginlist_append_must_produce_error = CPL_FALSE;
    cpl_plugin_set_init_must_produce_error = CPL_TRUE;
    cpl_test_eq_error(er_python_get_plugin_list(list),
                      CPL_ERROR_ILLEGAL_OUTPUT);

    cpl_pluginlist_delete(list);
    er_python_cleanup();

    /* This second edge case is similar to the above one, except we now test
       the error behaviour of er_python_get_plugin_list when the function
       cpl_pluginlist_append fails instead. */
    plugin = cpl_plugin_new();
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);

    pluginlist = cpl_pluginlist_new();
    cpl_test_assert(pluginlist != NULL);
    cpl_test_assert(cpl_pluginlist_append(pluginlist, plugin)
                    == CPL_ERROR_NONE);

    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    module_name = cpl_strdup("test.py");
    cx_map_insert(module_list, module_name, pluginlist);
    cx_map_insert(module_lut, plugin, module_name);

    cpl_test_assert(er_python_select_module("test.py") == CPL_ERROR_NONE);

    list = cpl_pluginlist_new();
    cpl_test_assert(list != NULL);

    cpl_pluginlist_append_must_produce_error = CPL_TRUE;
    cpl_plugin_set_init_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(er_python_get_plugin_list(list),
                      CPL_ERROR_ILLEGAL_OUTPUT);

    cpl_pluginlist_delete(list);
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
    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);

    recipe2->frames = cpl_frameset_new();
    cpl_test_assert(recipe2->frames != NULL);
    frame = cpl_frame_new();
    cpl_test_assert(frame != NULL);
    cpl_test_assert(cpl_frame_set_filename(frame, "output.fits")
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_tag(frame, "PROD") == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_type(frame, CPL_FRAME_TYPE_IMAGE)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_group(frame, CPL_FRAME_GROUP_PRODUCT)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frame_set_level(frame, CPL_FRAME_LEVEL_FINAL)
                    == CPL_ERROR_NONE);
    cpl_test_assert(cpl_frameset_insert(recipe2->frames, frame)
                    == CPL_ERROR_NONE);

    /* Test error handling if the cpl_frameset_insert function returns an
       error while calling update_plugin. */
    cpl_frameset_insert_must_produce_error = CPL_TRUE;
    cpl_test_eq_error(update_plugin(&recipe->interface, &recipe2->interface),
                      CPL_ERROR_ILLEGAL_OUTPUT);
    cpl_frameset_insert_must_produce_error = CPL_FALSE;

    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(&recipe->interface);
    cpl_frameset_delete(recipe2->frames);
    cpl_plugin_delete(&recipe2->interface);
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
