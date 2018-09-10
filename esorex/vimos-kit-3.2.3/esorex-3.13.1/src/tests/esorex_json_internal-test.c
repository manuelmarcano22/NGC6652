/*
 * This file is part of the ESO Recipe Execution Tool
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/* Unit tests for testing various edge cases of internal static functions. */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef ENABLE_PYTHON_RECIPES

#include <string.h>
#include <float.h>

#include "er_json.c"
#include "esorex_test.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/* Tests if two CPL recipe configuration objects are equal or not, i.e. are all
   their attributes the same or not.
 */
#define test_equal_config(first, second) \
    { \
        char ** tag = NULL; \
        char ** ptr = NULL; \
        char ** first_tags = cpl_recipeconfig_get_tags(first); \
        cpl_test_assert(first_tags != NULL); \
        char ** second_tags = cpl_recipeconfig_get_tags(second); \
        cpl_test_assert(second_tags != NULL); \
        char ** first_tag = first_tags; \
        char ** second_tag = second_tags; \
        while (*first_tag != NULL && *second_tag != NULL) { \
            cpl_test_eq_string(*first_tag, *second_tag); \
            ++first_tag; \
            ++second_tag; \
        } \
        cpl_test_assert(*first_tag == NULL && *second_tag == NULL); \
        for (tag = first_tags; *tag != NULL; ++tag) { \
            cpl_test_eq( \
                    cpl_recipeconfig_get_min_count(first, *tag, *tag), \
                    cpl_recipeconfig_get_min_count(second, *tag, *tag) \
                ); \
            cpl_test_eq( \
                    cpl_recipeconfig_get_max_count(first, *tag, *tag), \
                    cpl_recipeconfig_get_max_count(second, *tag, *tag) \
                ); \
            char ** first_inputs = cpl_recipeconfig_get_inputs(first, *tag); \
            cpl_test_assert(first_inputs != NULL); \
            char ** second_inputs = cpl_recipeconfig_get_inputs(second, *tag); \
            cpl_test_assert(second_inputs != NULL); \
            char ** first_outputs = cpl_recipeconfig_get_outputs(first, *tag); \
            cpl_test_assert(first_outputs != NULL); \
            char ** second_outputs = cpl_recipeconfig_get_outputs(second,*tag);\
            cpl_test_assert(second_outputs != NULL); \
            char ** first_input = first_inputs; \
            char ** second_input = second_inputs; \
            while (*first_input != NULL && *second_input != NULL) { \
                cpl_test_eq_string(*first_input, *second_input); \
                ++first_input; \
                ++second_input; \
            } \
            cpl_test_assert(*first_input == NULL && *second_input == NULL); \
            char ** first_output = first_outputs; \
            char ** second_output = second_outputs; \
            while (*first_output != NULL && *second_output != NULL) { \
                cpl_test_eq_string(*first_output, *second_output); \
                ++first_output; \
                ++second_output; \
            } \
            cpl_test_assert(*first_output == NULL && *second_output == NULL); \
            char ** input = NULL; \
            for (input = first_inputs; *input != NULL; ++input) { \
                cpl_test_eq( \
                        cpl_recipeconfig_get_min_count(first, *tag, *input), \
                        cpl_recipeconfig_get_min_count(second, *tag, *input) \
                    ); \
                cpl_test_eq( \
                        cpl_recipeconfig_get_max_count(first, *tag, *input), \
                        cpl_recipeconfig_get_max_count(second, *tag, *input) \
                    ); \
            } \
            for (ptr = first_inputs; *ptr != NULL; ++ptr) cpl_free(*ptr); \
            cpl_free(first_inputs); \
            for (ptr = second_inputs; *ptr != NULL; ++ptr) cpl_free(*ptr); \
            cpl_free(second_inputs); \
            for (ptr = first_outputs; *ptr != NULL; ++ptr) cpl_free(*ptr); \
            cpl_free(first_outputs); \
            for (ptr = second_outputs; *ptr != NULL; ++ptr) cpl_free(*ptr); \
            cpl_free(second_outputs); \
        } \
        for (ptr = first_tags; *ptr != NULL; ++ptr) cpl_free(*ptr); \
        cpl_free(first_tags); \
        for (ptr = second_tags; *ptr != NULL; ++ptr) cpl_free(*ptr); \
        cpl_free(second_tags); \
    }

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void test_get_value_from_object(void);
static void test_get_string_from_object(void);
static void test_get_ulong_from_object(void);
static void test_get_choices_array(void);
static void test_node_is_integer(void);
static void test_json_to_frame(void);
static void test_add_recipeconfig_inputs(void);
static void test_add_recipeconfig_outputs(void);
static void test_add_recipeconfig_tag(void);
static void test_json_to_parameter(void);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    test_get_value_from_object();
    test_get_string_from_object();
    test_get_ulong_from_object();
    test_get_choices_array();
    test_node_is_integer();
    test_json_to_frame();
    test_add_recipeconfig_inputs();
    test_add_recipeconfig_outputs();
    test_add_recipeconfig_tag();
    test_json_to_parameter();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

static void test_get_value_from_object(void)
{
    const char * json = "{}";
    er_json_node * node = er_json_parse(json);

    /* Check handling of error conditions. */
    cpl_test_null(get_value_from_object(node, "key", json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    er_json_node_delete(node);
}

static void test_get_string_from_object(void)
{
    const char * json = NULL;
    er_json_node * node = NULL;

    /* Check handling of error conditions. */
    json = "{}";
    node = er_json_parse(json);
    cpl_test_null(get_string_from_object(node, "key", json, CPL_TRUE));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(get_string_from_object(node, "key", json, CPL_FALSE));
    cpl_test_error(CPL_ERROR_NONE);
    er_json_node_delete(node);

    /* Check handling of error if the child key node is not a string. */
    json = "{\"key\": 123}";
    node = er_json_parse(json);
    cpl_test_assert(node != NULL);
    cpl_test_null(get_string_from_object(node, "key", json, CPL_TRUE));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(node);
}

static void test_get_ulong_from_object(void)
{
    const char * json = NULL;
    er_json_node * node = NULL;

    /* Check handling of error conditions. */
    json = "{}";
    node = er_json_parse(json);
    cpl_test_eq(get_ulong_from_object(node, "key", json, CPL_TRUE), 0);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_eq(get_ulong_from_object(node, "key", json, CPL_FALSE), 0);
    cpl_test_error(CPL_ERROR_NONE);
    er_json_node_delete(node);

    /* Check handling of error if the child key node is not a number. */
    json = "{\"key\": null}";
    node = er_json_parse(json);
    cpl_test_assert(node != NULL);
    cpl_test_eq(get_ulong_from_object(node, "key", json, CPL_TRUE), 0);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(node);
}

static void test_get_choices_array(void)
{
    const char * json = NULL;
    er_json_node * node = NULL;

    /* Check error handling when the 'choices' key is missing. */
    json = "{}";
    node = er_json_parse(json);
    cpl_test_null(get_choices_array(node, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(node);

    /* Check error handling when 'choices' is not an array. */
    json = "{\"choices\": \"value\"}";
    node = er_json_parse(json);
    cpl_test_null(get_choices_array(node, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(node);
}


static void test_node_is_integer(void)
{
    const char * json = NULL;
    er_json_node * node = NULL;

    /* Test that node_is_integer returns the correct value for various number
       strings. */
    json = "12";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_TRUE);
    er_json_node_delete(node);

    json = "+12";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_TRUE);
    er_json_node_delete(node);

    json = "-12";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_TRUE);
    er_json_node_delete(node);

    json = "12 ]";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_TRUE);
    er_json_node_delete(node);

    json = "12}";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_TRUE);
    er_json_node_delete(node);

    json = "12,  ";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_TRUE);
    er_json_node_delete(node);

    json = "1.0";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_FALSE);
    er_json_node_delete(node);

    json = "123.";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_FALSE);
    er_json_node_delete(node);

    json = "1.23e4";
    node = er_json_node_new_object(json);
    cpl_test_eq(node_is_integer(node), CPL_FALSE);
    er_json_node_delete(node);
}

static void test_json_to_frame(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;

    /* Test that json_to_frame produces an error if keywords are missing in
       the JSON input text. */
    json =
        "{\n"
        "    \"tag\": \"RAW\",\n"
        "    \"type\": 1,\n"
        "    \"group\": 0,\n"
        "    \"level\": 0\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_frame(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"filename\": \"test.fits\",\n"
        "    \"type\": 1,\n"
        "    \"group\": 0,\n"
        "    \"level\": 0\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_frame(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"filename\": \"test.fits\",\n"
        "    \"tag\": \"RAW\",\n"
        "    \"group\": 0,\n"
        "    \"level\": 0\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_frame(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"filename\": \"test.fits\",\n"
        "    \"tag\": \"RAW\",\n"
        "    \"type\": 1,\n"
        "    \"level\": 0\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_frame(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"filename\": \"test.fits\",\n"
        "    \"tag\": \"RAW\",\n"
        "    \"type\": 1,\n"
        "    \"group\": 0\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_frame(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test for correct error handling if the supposed frame structure is not
       a JSON object (i.e. dictionary). */
    json = "[\"dummy\", 123]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_frame(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);
}

static void test_add_recipeconfig_inputs(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;
    cpl_recipeconfig * config = NULL;
    cpl_recipeconfig * ref_config = NULL;

    ref_config = cpl_recipeconfig_new();
    cpl_test_assert(ref_config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(ref_config, "RAW", 1, 1) == 0);

    /* Test that add_recipeconfig_inputs produces an error it does not get an
       array structure in the JSON text. */
    json = "{}";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_inputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Test error handling if the wrong format is used in an array item. */
    json = "[\"dummy\"]";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_inputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Test that we get an error if trying to add an input that has the same
       tag value as the primary top level tag. */
    json =
        "[\n"
        "  {\n"
        "    \"tag\": \"RAW\",\n"
        "    \"min\": 1,\n"
        "    \"max\": 1\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_inputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Check error handling if the 'tag' field is an empty string. */
    json =
        "[\n"
        "  {\n"
        "    \"tag\": \"\",\n"
        "    \"min\": 1,\n"
        "    \"max\": 1\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_inputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Check error handling if keywords are missing. */
    json =
        "[\n"
        "  {\n"
        "    \"min\": 1,\n"
        "    \"max\": 1\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_inputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    json =
        "[\n"
        "  {\n"
        "    \"tag\": \"CALIB\",\n"
        "    \"max\": 1\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_inputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    json =
        "[\n"
        "  {\n"
        "    \"tag\": \"CALIB\",\n"
        "    \"min\": 1\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_inputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    cpl_recipeconfig_delete(ref_config);
}

static void test_add_recipeconfig_outputs(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;
    cpl_recipeconfig * config = NULL;
    cpl_recipeconfig * ref_config = NULL;

    ref_config = cpl_recipeconfig_new();
    cpl_test_assert(ref_config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(ref_config, "RAW", 1, 1) == 0);

    /* Test that add_recipeconfig_outputs produces an error it does not get an
       array structure in the JSON text. */
    json = "{}";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_outputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Test error handling if the wrong format is used in an array item. */
    json = "[123]";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_outputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Test error handling if an output entry is an empty string. */
    json = "[\"\"]";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_assert(cpl_recipeconfig_set_tag(config, "RAW", 1, 1) == 0);
    cpl_test_eq_error(add_recipeconfig_outputs(parsetree, json, config, "RAW"),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    cpl_recipeconfig_delete(ref_config);
}

static void test_add_recipeconfig_tag(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;
    cpl_recipeconfig * config = NULL;
    cpl_recipeconfig * ref_config = NULL;

    ref_config = cpl_recipeconfig_new();
    cpl_test_assert(ref_config != NULL);

    /* Test that add_recipeconfig_tag produces an error it does not get a
       JSON object (dictionary) structure. */
    json = "[]";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Test error handling if the parsed JSON has the wrong format. */
    json = "{\"dummy\": 123}";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Check error handling if required fields are missing. */
    json =
        "{\n"
        "    \"min\": 1,\n"
        "    \"max\": 1\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"tag\": \"RAW\",\n"
        "    \"max\": 1\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"tag\": \"RAW\",\n"
        "    \"min\": 1\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Test error handling if the tag has an empty string. */
    json =
        "{\n"
        "    \"tag\": \"\",\n"
        "    \"min\": 1,\n"
        "    \"max\": 1\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_INPUT);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Check for correct error handling if the inputs or outputs are invalid.
       Note that we do not compare to the reference structure with the
       test_equal_config macro, because we know that the 'config' object will
       be only partially constructed. */
    json =
        "{\n"
        "    \"tag\": \"RAW\",\n"
        "    \"min\": 1,\n"
        "    \"max\": 1,\n"
        "    \"inputs\": [\"dummy\"],\n"
        "    \"outputs\": [\"PROD\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

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
        "    \"outputs\": [{\"dummy\": 123}]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_ILLEGAL_INPUT);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Check for correct parsing of a valid recipe config structure. */
    cpl_test_assert(cpl_recipeconfig_set_tag(ref_config, "RAW", 1, 1) == 0);
    json =
        "{\n"
        "    \"tag\": \"RAW\",\n"
        "    \"min\": 1,\n"
        "    \"max\": 1\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    config = cpl_recipeconfig_new();
    cpl_test_assert(config != NULL);
    cpl_test_eq_error(add_recipeconfig_tag(parsetree, json, config),
                      CPL_ERROR_NONE);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    /* Check for correct parsing of a recipe config structure with inputs and
       outputs. */
    cpl_test_assert(cpl_recipeconfig_set_input(ref_config, "RAW", "CALIB", 2, 3)
                    == 0);
    cpl_test_assert(cpl_recipeconfig_set_output(ref_config, "RAW", "PROD")
                    == 0);
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
                      CPL_ERROR_NONE);
    test_equal_config(config, ref_config);
    cpl_recipeconfig_delete(config);
    er_json_node_delete(parsetree);

    cpl_recipeconfig_delete(ref_config);
}

static void test_json_to_parameter(void)
{
    const char * json = NULL;
    er_json_node * parsetree = NULL;
    cpl_parameter * param = NULL;
    cpl_parameter * ref_param = NULL;

    /* Test that an integer value parameter is parsed correctly. */
    ref_param = cpl_parameter_new_value("testpar", CPL_TYPE_INT,
                                        "test param", "test", 3);
    cpl_test_assert(ref_param != NULL);
    cpl_parameter_set_id(ref_param, 2);
    cpl_parameter_set_int(ref_param, 5);
    cpl_parameter_set_default_flag(ref_param, 1);
    cpl_parameter_set_tag(ref_param, "mypar");
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CLI, "par_cli");
    cpl_parameter_enable(ref_param, CPL_PARAMETER_MODE_CLI);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_ENV, "par_env");
    cpl_parameter_enable(ref_param, CPL_PARAMETER_MODE_ENV);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CFG, "par_cfg");
    cpl_parameter_disable(ref_param, CPL_PARAMETER_MODE_CFG);
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"id\": 2,\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5,\n"
        "    \"default\": 3,\n"
        "    \"present\": true,\n"
        "    \"tag\": \"mypar\",\n"
        "    \"cli_enabled\": true,\n"
        "    \"cli_alias\": \"par_cli\",\n"
        "    \"env_enabled\": true,\n"
        "    \"env_alias\": \"par_env\",\n"
        "    \"cfg_enabled\": false,\n"
        "    \"cfg_alias\": \"par_cfg\"\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Test that a floating point value parameter is parsed correctly. */
    ref_param = cpl_parameter_new_value("testpar", CPL_TYPE_DOUBLE,
                                        "test param", "test", 3.5);
    cpl_test_assert(ref_param != NULL);
    cpl_parameter_set_id(ref_param, 3);
    cpl_parameter_set_double(ref_param, 5.2);
    cpl_parameter_set_default_flag(ref_param, 0);
    cpl_parameter_set_tag(ref_param, "mypar");
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CLI, "par_cli");
    cpl_parameter_disable(ref_param, CPL_PARAMETER_MODE_CLI);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_ENV, "par_env");
    cpl_parameter_disable(ref_param, CPL_PARAMETER_MODE_ENV);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CFG, "par_cfg");
    cpl_parameter_enable(ref_param, CPL_PARAMETER_MODE_CFG);
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"id\": 3,\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5.2,\n"
        "    \"default\": 3.5,\n"
        "    \"present\": false,\n"
        "    \"tag\": \"mypar\",\n"
        "    \"cli_enabled\": false,\n"
        "    \"cli_alias\": \"par_cli\",\n"
        "    \"env_enabled\": false,\n"
        "    \"env_alias\": \"par_env\",\n"
        "    \"cfg_enabled\": true,\n"
        "    \"cfg_alias\": \"par_cfg\"\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Check that a missing description in the JSON does not cause errors. */
    ref_param = cpl_parameter_new_value("testpar", CPL_TYPE_INT, "", "test", 3);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Check that a missing context in the JSON does not cause errors. */
    ref_param = cpl_parameter_new_value("testpar", CPL_TYPE_INT,
                                        "test param", "", 3);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Check error handling if an invalid JSON structure is given to the
       json_to_parameter function. */
    json = "\"dummy\"";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if class, name or default keywords are missing. */
    json = "{}";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\"\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check for correct error handling if 'default' is invalid. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": {}\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check for correct error handling if there is a type mismatch between
       'value' and 'default'. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": true,\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 3,\n"
        "    \"default\": true\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": \"string\",\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 3,\n"
        "    \"default\": \"string\"\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 3.12,\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": \"string\",\n"
        "    \"default\": 3.12\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": {},\n"
        "    \"default\": 123\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    /* Check for correct error handling if 'class' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"unknown\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if 'present' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"present\": 123\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if 'id' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"id\": true,\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if 'tag' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"tag\": 1.234\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check that a null value for 'tag' does not produce an error. */
    ref_param = cpl_parameter_new_value("testpar", CPL_TYPE_INT,
                                        "test param", "test", 3);
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"tag\": null\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Test error handling if 'cli_enabled' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"cli_enabled\": 123\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if 'cli_alias' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"cli_alias\": 123\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if 'env_enabled' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"env_enabled\": 123\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if 'env_alias' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"env_alias\": 123\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if 'cfg_enabled' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"cfg_enabled\": 123\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if 'cfg_alias' has an invalid value. */
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 3,\n"
        "    \"cfg_alias\": 123\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check that no error occurs if 'cli_alias' is null. */
    ref_param = cpl_parameter_new_value("testpar", CPL_TYPE_INT,
                                        "test param", "test", 4);
    cpl_test_assert(ref_param != NULL);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CLI, NULL);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_ENV, "abc");
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CFG, "def");
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 4,\n"
        "    \"cli_alias\": null,\n"
        "    \"env_alias\": \"abc\",\n"
        "    \"cfg_alias\": \"def\"\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Check that no error occurs if 'env_alias' is null. */
    ref_param = cpl_parameter_new_value("testpar", CPL_TYPE_INT,
                                        "test param", "test", 4);
    cpl_test_assert(ref_param != NULL);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CLI, "abc");
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_ENV, NULL);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CFG, "def");
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 4,\n"
        "    \"cli_alias\": \"abc\",\n"
        "    \"env_alias\": null,\n"
        "    \"cfg_alias\": \"def\"\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Check that no error occurs if 'cfg_alias' is null. */
    ref_param = cpl_parameter_new_value("testpar", CPL_TYPE_INT,
                                        "test param", "test", 4);
    cpl_test_assert(ref_param != NULL);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CLI, "abc");
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_ENV, "def");
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CFG, NULL);
    json =
        "{\n"
        "    \"class\": \"value\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"test param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 4,\n"
        "    \"cli_alias\": \"abc\",\n"
        "    \"env_alias\": \"def\",\n"
        "    \"cfg_alias\": null\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Test that an integer range parameter is parsed correctly. */
    ref_param = cpl_parameter_new_range("testpar", CPL_TYPE_INT,
                                        "range param", "test", 3, 2, 7);
    cpl_test_assert(ref_param != NULL);
    cpl_parameter_set_id(ref_param, 6);
    cpl_parameter_set_int(ref_param, 5);
    cpl_parameter_set_default_flag(ref_param, 1);
    cpl_parameter_set_tag(ref_param, "myrange");
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CLI, "par_cli");
    cpl_parameter_enable(ref_param, CPL_PARAMETER_MODE_CLI);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_ENV, "par_env");
    cpl_parameter_enable(ref_param, CPL_PARAMETER_MODE_ENV);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CFG, "par_cfg");
    cpl_parameter_disable(ref_param, CPL_PARAMETER_MODE_CFG);
    json =
        "{\n"
        "    \"class\": \"range\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"id\": 6,\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5,\n"
        "    \"default\": 3,\n"
        "    \"min\": 2,\n"
        "    \"max\": 7,\n"
        "    \"present\": true,\n"
        "    \"tag\": \"myrange\",\n"
        "    \"cli_enabled\": true,\n"
        "    \"cli_alias\": \"par_cli\",\n"
        "    \"env_enabled\": true,\n"
        "    \"env_alias\": \"par_env\",\n"
        "    \"cfg_enabled\": false,\n"
        "    \"cfg_alias\": \"par_cfg\"\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Test error handling if the 'min' or 'max' keywords are missing. */
    json =
        "{\n"
        "    \"class\": \"range\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5,\n"
        "    \"default\": 3,\n"
        "    \"max\": 7\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"range\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5,\n"
        "    \"default\": 3,\n"
        "    \"min\": 2\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if the 'min' or 'max' keywords are invalid. */
    json =
        "{\n"
        "    \"class\": \"range\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5,\n"
        "    \"default\": 3,\n"
        "    \"min\": {},\n"
        "    \"max\": 7\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"range\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5,\n"
        "    \"default\": 3,\n"
        "    \"min\": 2,\n"
        "    \"max\": {}\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    /* Test error handling if the 'min' or 'max' have types that do not match.
     */
    json =
        "{\n"
        "    \"class\": \"range\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5,\n"
        "    \"default\": 3,\n"
        "    \"min\": 2.3,\n"
        "    \"max\": 7\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"range\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 5,\n"
        "    \"default\": 3,\n"
        "    \"min\": 2,\n"
        "    \"max\": 7.1\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    /* Test for error if a value type is used that is not compatible with a
       range parameter. */
    json =
        "{\n"
        "    \"class\": \"range\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": true,\n"
        "    \"default\": true,\n"
        "    \"min\": true,\n"
        "    \"max\": false\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test that an integer enum parameter is parsed correctly. */
    ref_param = cpl_parameter_new_enum("testpar", CPL_TYPE_INT,
                                       "enum param", "test", 3, 3, 2, 3, 7);
    cpl_test_assert(ref_param != NULL);
    cpl_parameter_set_id(ref_param, 3);
    cpl_parameter_set_int(ref_param, 7);
    cpl_parameter_set_default_flag(ref_param, 1);
    cpl_parameter_set_tag(ref_param, "myenum");
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CLI, "par_cli");
    cpl_parameter_enable(ref_param, CPL_PARAMETER_MODE_CLI);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_ENV, "par_env");
    cpl_parameter_enable(ref_param, CPL_PARAMETER_MODE_ENV);
    cpl_parameter_set_alias(ref_param, CPL_PARAMETER_MODE_CFG, "par_cfg");
    cpl_parameter_disable(ref_param, CPL_PARAMETER_MODE_CFG);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"id\": 3,\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 7,\n"
        "    \"default\": 3,\n"
        "    \"choices\": [2, 3, 7],\n"
        "    \"present\": true,\n"
        "    \"tag\": \"myenum\",\n"
        "    \"cli_enabled\": true,\n"
        "    \"cli_alias\": \"par_cli\",\n"
        "    \"env_enabled\": true,\n"
        "    \"env_alias\": \"par_env\",\n"
        "    \"cfg_enabled\": false,\n"
        "    \"cfg_alias\": \"par_cfg\"\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    /* Test for error if a value type is used that is not compatible with an
       enumeration parameter. */
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": true,\n"
        "    \"default\": true,\n"
        "    \"choices\": [true, false]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if 'choices' does not contain the default value. */
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 4,\n"
        "    \"default\": 3,\n"
        "    \"choices\": [2, 4]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 4.5,\n"
        "    \"default\": 3.5,\n"
        "    \"choices\": [4.5, 5.5]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": \"x\",\n"
        "    \"default\": \"y\",\n"
        "    \"choices\": [\"x\", \"z\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if 'choices' is invalid for different types. */
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 4,\n"
        "    \"default\": 3,\n"
        "    \"choices\": [2.0, 3.0, 4.0]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 4,\n"
        "    \"default\": 3,\n"
        "    \"choices\": [\"x\", \"y\", \"z\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 4.0,\n"
        "    \"default\": 3.0,\n"
        "    \"choices\": [\"x\", \"y\", \"z\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": \"x\",\n"
        "    \"default\": \"y\",\n"
        "    \"choices\": [2, 3, 4]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_TYPE_MISMATCH);
    er_json_node_delete(parsetree);

    /* Check error handling if 'choices' is invalid for different types. */
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 4,\n"
        "    \"default\": 3,\n"
        "    \"choices\": {}\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": 4.1,\n"
        "    \"default\": 3.2,\n"
        "    \"choices\": {}\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"testpar\",\n"
        "    \"description\": \"range param\",\n"
        "    \"context\": \"test\",\n"
        "    \"value\": \"x\",\n"
        "    \"default\": \"y\",\n"
        "    \"choices\": {}\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(json_to_parameter(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check correct parsing of enumeration parameters with different number of
       choices. We check up to the minimum number of choices supported. */
    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 1, 1);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 2, 1, 2);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 3, 1, 2, 3);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2, 3]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 4, 1, 2, 3, 4);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2, 3, 4]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 5, 1, 2, 3, 4, 5);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2, 3, 4, 5]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 6, 1, 2, 3, 4, 5, 6);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2, 3, 4, 5, 6]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 7, 1, 2, 3, 4, 5, 6, 7);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2, 3, 4, 5, 6, 7]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 8, 1, 2, 3, 4, 5, 6, 7, 8);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2, 3, 4, 5, 6, 7, 8]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2, 3, 4, 5, 6, 7, 8, 9]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_INT, "enum param",
                                       "test", 1,
                                       10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1,\n"
        "    \"choices\": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 1, 1.1);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 2, 1.1, 2.2);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 3, 1.1, 2.2, 3.3);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2, 3.3]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 4, 1.1, 2.2, 3.3, 4.4);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2, 3.3, 4.4]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 5, 1.1, 2.2, 3.3, 4.4, 5.5);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2, 3.3, 4.4, 5.5]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 6, 1.1, 2.2, 3.3, 4.4, 5.5,
                                       6.6);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2, 3.3, 4.4, 5.5, 6.6]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 7, 1.1, 2.2, 3.3, 4.4, 5.5,
                                       6.6, 7.7);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 8, 1.1, 2.2, 3.3, 4.4, 5.5,
                                       6.6, 7.7, 8.8);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 9, 1.1, 2.2, 3.3, 4.4, 5.5,
                                       6.6, 7.7, 8.8, 9.9);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_DOUBLE, "enum param",
                                       "test", 1.1, 10, 1.1, 2.2, 3.3, 4.4, 5.5,
                                       6.6, 7.7, 8.8, 9.9, 10.0);
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": 1.1,\n"
        "    \"choices\": [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.0]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 1, "a");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 2, "a", "b");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 3, "a", "b", "c");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\", \"c\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 4, "a", "b", "c", "d");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\", \"c\", \"d\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 5, "a", "b", "c", "d", "e");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\", \"c\", \"d\", \"e\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 6, "a", "b", "c", "d", "e",
                                       "f");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\", \"c\", \"d\", \"e\", \"f\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 7, "a", "b", "c", "d", "e",
                                       "f", "h");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"h\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 8, "a", "b", "c", "d", "e",
                                       "f", "h", "i");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"h\",\n"
        "                  \"i\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 9, "a", "b", "c", "d", "e",
                                       "f", "h", "i", "j");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"h\",\n"
        "                  \"i\", \"j\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);

    ref_param = cpl_parameter_new_enum("par", CPL_TYPE_STRING, "enum param",
                                       "test", "a", 10, "a", "b", "c", "d", "e",
                                       "f", "h", "i", "j", "k");
    cpl_test_assert(ref_param != NULL);
    json =
        "{\n"
        "    \"class\": \"enum\",\n"
        "    \"name\": \"par\",\n"
        "    \"description\": \"enum param\",\n"
        "    \"context\": \"test\",\n"
        "    \"default\": \"a\",\n"
        "    \"choices\": [\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"h\",\n"
        "                  \"i\", \"j\", \"k\"]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    param = json_to_parameter(parsetree, json);
    cpl_test_error(CPL_ERROR_NONE);
    test_equal_param(param, ref_param);
    cpl_parameter_delete(param);
    cpl_parameter_delete(ref_param);
    er_json_node_delete(parsetree);
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
