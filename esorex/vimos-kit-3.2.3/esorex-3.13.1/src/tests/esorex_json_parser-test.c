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


/* Undefine NDEBUG to make sure assert() macros are compiled in for these
   tests. */
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <math.h>
#include <cpl_test.h>
#include "er_json.c"

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void test_json_find_line_column(void);
static void test_json_node_null(void);
static void test_json_node_bool(void);
static void test_json_node_number(void);
static void test_json_node_string(void);
static void test_json_node_array(void);
static void test_json_node_object(void);
static void test_json_node_invalid_access(void);
static void test_simple_create_and_delete_node(void);
static void test_parsing_whitespace(void);
static void test_parsing_null(void);
static void test_parsing_boolean(void);
static void test_parsing_number(void);
static void test_parsing_string(void);
static void test_parsing_array(void);
static void test_parsing_object(void);
static void test_parse_errors(void);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    test_json_find_line_column();
    test_json_node_null();
    test_json_node_bool();
    test_json_node_number();
    test_json_node_string();
    test_json_node_array();
    test_json_node_object();
    test_json_node_invalid_access();
    test_simple_create_and_delete_node();
    test_parsing_whitespace();
    test_parsing_null();
    test_parsing_boolean();
    test_parsing_number();
    test_parsing_string();
    test_parsing_array();
    test_parsing_object();
    test_parse_errors();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

static void test_json_find_line_column(void)
{
    int line, column;
    const char * json = "a\nbc\nd";

    /* Test the calculation of the line and column numbers for various positions
       within the JSON text. */
    cpl_test_eq_error(er_json_find_line_column(json, json, &line, &column),
                      CPL_ERROR_NONE);
    cpl_test_eq(line, 1);
    cpl_test_eq(column, 1);

    cpl_test_eq_error(er_json_find_line_column(json, json+1, &line, &column),
                      CPL_ERROR_NONE);
    cpl_test_eq(line, 1);
    cpl_test_eq(column, 2);

    cpl_test_eq_error(er_json_find_line_column(json, json+2, &line, &column),
                      CPL_ERROR_NONE);
    cpl_test_eq(line, 2);
    cpl_test_eq(column, 1);

    cpl_test_eq_error(er_json_find_line_column(json, json+3, &line, &column),
                      CPL_ERROR_NONE);
    cpl_test_eq(line, 2);
    cpl_test_eq(column, 2);

    cpl_test_eq_error(er_json_find_line_column(json, json+4, &line, &column),
                      CPL_ERROR_NONE);
    cpl_test_eq(line, 2);
    cpl_test_eq(column, 3);

    cpl_test_eq_error(er_json_find_line_column(json, json+5, &line, &column),
                      CPL_ERROR_NONE);
    cpl_test_eq(line, 3);
    cpl_test_eq(column, 1);

    cpl_test_eq_error(er_json_find_line_column(json, json+6, &line, &column),
                      CPL_ERROR_NONE);
    cpl_test_eq(line, 3);
    cpl_test_eq(column, 2);

    /* Check error handling when invalid pointers are passed to the function. */
    cpl_test_eq_error(er_json_find_line_column(NULL, json, &line, &column),
                      CPL_ERROR_NULL_INPUT);
    cpl_test_eq_error(er_json_find_line_column(json, NULL, &line, &column),
                      CPL_ERROR_NULL_INPUT);
    cpl_test_eq_error(er_json_find_line_column(json, json, NULL, &column),
                      CPL_ERROR_NULL_INPUT);
    cpl_test_eq_error(er_json_find_line_column(json, json, &line, NULL),
                      CPL_ERROR_NULL_INPUT);
}

static void test_json_node_null(void)
{
    /* Test construction and destruction of null type JSON value objects. */
    const char * json = "null";
    er_json_node * obj = er_json_node_new_null(json);
    cpl_test_nonnull(obj);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(er_json_node_type(obj), JSON_NULL);
    cpl_test_eq_ptr(er_json_node_location(obj), json);
    cpl_test_eq_ptr(obj->object, NULL);
    cpl_test_error(CPL_ERROR_NONE);

    er_json_node_delete(obj);
    cpl_test_error(CPL_ERROR_NONE);
}

static void test_json_node_bool(void)
{
    /* Test construction and destruction of boolean type JSON value objects. */
    const char * json = "true";
    er_json_node * obj = er_json_node_new_bool(TRUE, json);
    cpl_test_nonnull(obj);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(er_json_node_type(obj), JSON_BOOL);
    cpl_test_eq_ptr(er_json_node_location(obj), json);
    cpl_test_eq(obj->bool, TRUE);
    cpl_test_error(CPL_ERROR_NONE);

    er_json_node_delete(obj);
    cpl_test_error(CPL_ERROR_NONE);

    json = "false";
    obj = er_json_node_new_bool(FALSE, json);
    cpl_test_nonnull(obj);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(er_json_node_type(obj), JSON_BOOL);
    cpl_test_eq_ptr(er_json_node_location(obj), json);
    cpl_test_eq(obj->bool, FALSE);
    cpl_test_error(CPL_ERROR_NONE);

    er_json_node_delete(obj);
    cpl_test_error(CPL_ERROR_NONE);
}

static void test_json_node_number(void)
{
    /* Test construction and destruction of number type JSON value objects. */
    const char * json = "123.4";
    er_json_node * obj = er_json_node_new_number(123.4, json);
    cpl_test_nonnull(obj);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(er_json_node_type(obj), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(obj), json);
    cpl_test_abs(obj->number, 123.4, DBL_EPSILON);
    cpl_test_error(CPL_ERROR_NONE);

    er_json_node_delete(obj);
    cpl_test_error(CPL_ERROR_NONE);
}

static void test_json_node_string(void)
{
    /* Test construction and destruction of string type JSON value objects. */
    const char * json = "\"test\"";
    er_json_node * obj = er_json_node_new_string("test", json);
    cpl_test_nonnull(obj);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(er_json_node_type(obj), JSON_STRING);
    cpl_test_eq_ptr(er_json_node_location(obj), json);
    cpl_test_nonnull(obj->string);
    cpl_test_eq_string(obj->string, "test");
    cpl_test_error(CPL_ERROR_NONE);

    er_json_node_delete(obj);
    cpl_test_error(CPL_ERROR_NONE);
}

static void test_json_node_array(void)
{
    /* Test construction and destruction of array type JSON value objects. */
    const char * json = "[1.2, true, \"test\"]";
    er_json_node * obj = er_json_node_new_array(json);
    cpl_test_nonnull(obj);
    cpl_test_error(CPL_ERROR_NONE);

    /* Created the empty array object. */
    cpl_test_eq(er_json_node_type(obj), JSON_ARRAY);
    cpl_test_eq_ptr(er_json_node_location(obj), json);
    cpl_test_nonnull(obj->array);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_array_empty(obj), TRUE);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_array_size(obj), 0);
    cpl_test_error(CPL_ERROR_NONE);

    /* Add the items to the array. */
    er_json_node_array_append(obj, er_json_node_new_number(1.2, json+1));
    cpl_test_error(CPL_ERROR_NONE);
    er_json_node_array_append(obj, er_json_node_new_bool(TRUE, json+7));
    cpl_test_error(CPL_ERROR_NONE);
    er_json_node_array_append(obj, er_json_node_new_string("test", json+12));
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_array_empty(obj), FALSE);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_array_size(obj), 3);
    cpl_test_error(CPL_ERROR_NONE);

    /* Now iterate through the array, checking the values as we go.
       Check the first item is the correct number. */
    er_json_array_iterator iter = er_json_node_array_begin(obj);
    cpl_test_nonnull(iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_noneq_ptr(iter, er_json_node_array_end(obj));
    cpl_test_error(CPL_ERROR_NONE);
    const er_json_node * item = er_json_node_array_get(obj, iter);
    cpl_test_nonnull(item);
    cpl_test_eq(er_json_node_type(item), JSON_NUMBER);
    cpl_test_abs(er_json_node_get_number(item), 1.2, DBL_EPSILON);

    /* Check the second item is a boolean. */
    iter = er_json_node_array_next(obj, iter);
    cpl_test_nonnull(iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_noneq_ptr(iter, er_json_node_array_end(obj));
    cpl_test_error(CPL_ERROR_NONE);
    item = er_json_node_array_get(obj, iter);
    cpl_test_nonnull(item);
    cpl_test_eq(er_json_node_type(item), JSON_BOOL);
    cpl_test_eq(er_json_node_get_bool(item), TRUE);

    /* Check the third item is the correct string. */
    iter = er_json_node_array_next(obj, iter);
    cpl_test_nonnull(iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_noneq_ptr(iter, er_json_node_array_end(obj));
    cpl_test_error(CPL_ERROR_NONE);
    item = er_json_node_array_get(obj, iter);
    cpl_test_nonnull(item);
    cpl_test_eq(er_json_node_type(item), JSON_STRING);
    cpl_test_eq_string(er_json_node_get_string(item), "test");

    /* Check we reached the end of the array. */
    iter = er_json_node_array_next(obj, iter);
    cpl_test_nonnull(iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq_ptr(iter, er_json_node_array_end(obj));

    /* Check that we can go back one item and get the string again. */
    iter = er_json_node_array_previous(obj, iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_noneq_ptr(iter, er_json_node_array_end(obj));
    cpl_test_eq_string(er_json_node_get_string(item), "test");

    er_json_node_delete(obj);
    cpl_test_error(CPL_ERROR_NONE);
}

static void test_json_node_object(void)
{
    /* Test construction and destruction of object type JSON value objects. */
    const char * json = "{\"num\": 1.2, \"str\": \"test\"}";
    er_json_node * obj = er_json_node_new_object(json);
    cpl_test_nonnull(obj);
    cpl_test_error(CPL_ERROR_NONE);

    /* Created the empty JSON object. */
    cpl_test_eq(er_json_node_type(obj), JSON_OBJECT);
    cpl_test_eq_ptr(er_json_node_location(obj), json);
    cpl_test_nonnull(obj->object);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_object_empty(obj), TRUE);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_object_size(obj), 0);
    cpl_test_error(CPL_ERROR_NONE);

    /* Add the key/value pairs to the object. */
    er_json_node_object_insert(obj, "num", er_json_node_new_number(1.2, json+1));
    cpl_test_error(CPL_ERROR_NONE);
    er_json_node_object_insert(obj, "str", er_json_node_new_string("test",json+13));
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_object_empty(obj), FALSE);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_object_size(obj), 2);
    cpl_test_error(CPL_ERROR_NONE);

    /* Check that the keys can be found and their values are correct.
       We first use the er_json_node_object_get query method. */
    const er_json_node * item = er_json_node_object_get(obj, "num");
    cpl_test_nonnull(item);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_type(item), JSON_NUMBER);
    cpl_test_abs(er_json_node_get_number(item), 1.2, DBL_EPSILON);

    item = er_json_node_object_get(obj, "str");
    cpl_test_nonnull(item);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(er_json_node_type(item), JSON_STRING);
    cpl_test_eq_string(er_json_node_get_string(item), "test");

    item = er_json_node_object_get(obj, "unknown");
    cpl_test_null(item);
    cpl_test_error(CPL_ERROR_NONE);

    /* Check the iterator methods. Start with the first item. */
    er_json_object_iterator iter = er_json_node_object_begin(obj);
    cpl_test_nonnull(iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_noneq_ptr(iter, er_json_node_object_end(obj));
    cpl_test_error(CPL_ERROR_NONE);

    /* Get the second item. */
    iter = er_json_node_object_next(obj, iter);
    cpl_test_nonnull(iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_noneq_ptr(iter, er_json_node_object_end(obj));
    cpl_test_error(CPL_ERROR_NONE);

    /* Check that we hit the end of the items. */
    iter = er_json_node_object_next(obj, iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq_ptr(iter, er_json_node_object_end(obj));
    cpl_test_error(CPL_ERROR_NONE);

    /* Check that we can move back to get to the beginning. */
    iter = er_json_node_object_previous(obj, iter);
    cpl_test_error(CPL_ERROR_NONE);
    iter = er_json_node_object_previous(obj, iter);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq_ptr(iter, er_json_node_object_begin(obj));
    cpl_test_error(CPL_ERROR_NONE);

    /* Now we explicitly iterate through the object's entries to check them.
       Note: there is no guarantees on the order of the items, we only know
       that there must be 2, so we will have to check the types of the items
       dynamically as we go. */
    for (iter = er_json_node_object_begin(obj);
         iter != er_json_node_object_end(obj);
         iter = er_json_node_object_next(obj, iter))
    {
        cpl_test_nonnull(er_json_node_object_get_key(obj, iter));
        cpl_test_error(CPL_ERROR_NONE);
        item = er_json_node_object_get_value(obj, iter);
        cpl_test_nonnull(item);
        cpl_test_error(CPL_ERROR_NONE);
        if (strcmp(er_json_node_object_get_key(obj, iter), "num") == 0)
        {
            cpl_test_eq(er_json_node_type(item), JSON_NUMBER);
            cpl_test_abs(er_json_node_get_number(item), 1.2, DBL_EPSILON);
        }
        else
        {
            cpl_test_eq_string(er_json_node_object_get_key(obj, iter), "str");
            cpl_test_eq(er_json_node_type(item), JSON_STRING);
            cpl_test_eq_string(er_json_node_get_string(item), "test");
        }
    }

    er_json_node_delete(obj);
    cpl_test_error(CPL_ERROR_NONE);
}

static void test_json_node_invalid_access(void)
{
    /* Create some valid dummy objects for testing. */
    const char * nulljson = "null";
    er_json_node * nullnode = er_json_node_new_null(nulljson);
    cpl_test_assert(nullnode != NULL);
    const char * arrayjson = "[]";
    er_json_node * arraynode = er_json_node_new_array(arrayjson);
    er_json_array_iterator arrayiter = er_json_node_array_begin(arraynode);
    cpl_test_assert(arrayiter != NULL);
    const char * objectjson = "{}";
    er_json_node * objectnode = er_json_node_new_object(objectjson);
    er_json_object_iterator objectiter = er_json_node_object_begin(objectnode);
    cpl_test_assert(objectiter != NULL);

    /* Check error handling when NULL pointers or invalid objects are passed to
       the JSON node API methods. */

    cpl_test_eq(er_json_node_type(NULL), JSON_NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);

    cpl_test_eq_ptr(er_json_node_location(NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);

    cpl_test_eq(er_json_node_get_bool(NULL), CPL_FALSE);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(er_json_node_get_bool(nullnode), CPL_FALSE);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test(isnan(er_json_node_get_number(NULL)));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test(isnan(er_json_node_get_number(nullnode)));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_get_string(NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_get_string(nullnode), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq(er_json_node_array_size(NULL), -1);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(er_json_node_array_size(nullnode), -1);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq(er_json_node_array_empty(NULL), CPL_TRUE);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(er_json_node_array_empty(nullnode), CPL_TRUE);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_array_begin(NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_array_begin(nullnode), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_array_end(NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_array_end(nullnode), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_array_next(NULL, arrayiter), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_array_next(arraynode, NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_array_next(nullnode, arrayiter), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_array_previous(NULL, arrayiter), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_array_previous(arraynode, NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_array_previous(nullnode, arrayiter), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_array_get(NULL, arrayiter), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_array_get(arraynode, NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_array_get(nullnode, arrayiter), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq(er_json_node_object_size(NULL), -1);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(er_json_node_object_size(nullnode), -1);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq(er_json_node_object_empty(NULL), CPL_TRUE);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(er_json_node_object_empty(nullnode), CPL_TRUE);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_object_begin(NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_begin(nullnode), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_object_end(NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_end(nullnode), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_object_next(NULL, objectiter), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_next(objectnode, NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_next(nullnode, objectiter), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_object_previous(NULL, objectiter), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_previous(objectnode, NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_previous(nullnode, objectiter), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_object_get_key(NULL, objectiter), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_get_key(objectnode, NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_get_key(nullnode, objectiter), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_object_get_value(NULL, objectiter), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_get_value(objectnode, NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_get_value(nullnode, objectiter), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    cpl_test_eq_ptr(er_json_node_object_get(NULL, "testkey"), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_get(objectnode, NULL), NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq_ptr(er_json_node_object_get(nullnode, "testkey"), NULL);
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);

    er_json_node_delete(objectnode);
    er_json_node_delete(arraynode);
    er_json_node_delete(nullnode);
}

static void test_simple_create_and_delete_node(void)
{
    /* Test parsing an empty string returns a null JSON object. */
    const char * json = "";
    er_json_node * tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NULL);
    er_json_node_delete(tree);
    cpl_test_error(CPL_ERROR_NONE);

    /* Test for error codes if the parser gets a NULL. */
    tree = er_json_parse(NULL);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(tree);

    /* Test that the delete method can handle a NULL input. */
    er_json_node_delete(NULL);
    cpl_test_error(CPL_ERROR_NONE);
}

static void test_parsing_whitespace(void)
{
    /* Test parsing of various white space characters. */
    const char * json = "  ";
    er_json_node * tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NULL);
    er_json_node_delete(tree);

    json = " \n\t\r \v ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NULL);
    er_json_node_delete(tree);
}

static void test_parsing_null(void)
{
    /* Test parsing a null value. */
    const char * json = "null";
    er_json_node * tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NULL);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    er_json_node_delete(tree);

    json = " null";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NULL);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    er_json_node_delete(tree);

    json = "null  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NULL);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    er_json_node_delete(tree);

    json = "  null  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NULL);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    er_json_node_delete(tree);
}

static void test_parsing_boolean(void)
{
    /* Test parsing a boolean value. */
    const char * json = "true";
    er_json_node * tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_BOOL);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test_eq(er_json_node_get_bool(tree), TRUE);
    er_json_node_delete(tree);

    json = "  true";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_BOOL);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    cpl_test_eq(er_json_node_get_bool(tree), TRUE);
    er_json_node_delete(tree);

    json = "true ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_BOOL);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test_eq(er_json_node_get_bool(tree), TRUE);
    er_json_node_delete(tree);

    json = " true ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_BOOL);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test_eq(er_json_node_get_bool(tree), TRUE);
    er_json_node_delete(tree);

    json = " false ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_BOOL);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test_eq(er_json_node_get_bool(tree), FALSE);
    er_json_node_delete(tree);
}

static void test_parsing_number(void)
{
    /* Test parsing a number value. */
    const char * json = "123";
    er_json_node * tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test_abs(er_json_node_get_number(tree), 123.0, DBL_EPSILON);
    er_json_node_delete(tree);

    json = " 12.3";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test_abs(er_json_node_get_number(tree), 12.3, DBL_EPSILON);
    er_json_node_delete(tree);

    json = "12.3  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test_abs(er_json_node_get_number(tree), 12.3, DBL_EPSILON);
    er_json_node_delete(tree);

    json = "  12.3  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    cpl_test_abs(er_json_node_get_number(tree), 12.3, DBL_EPSILON);
    er_json_node_delete(tree);

    json = "  12.3e4";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    cpl_test_abs(er_json_node_get_number(tree), 12.3e4, DBL_EPSILON);
    er_json_node_delete(tree);

    json = "  -12.3e4 ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    cpl_test_abs(er_json_node_get_number(tree), -12.3e4, DBL_EPSILON);
    er_json_node_delete(tree);

    json = " nan ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test(isnan(er_json_node_get_number(tree)));
    er_json_node_delete(tree);

    json = " NaN ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test(isnan(er_json_node_get_number(tree)));
    er_json_node_delete(tree);

    json = "inf ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test(isinf(er_json_node_get_number(tree)));
    cpl_test(er_json_node_get_number(tree) > 0.0);
    er_json_node_delete(tree);

    json = " inf ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test(isinf(er_json_node_get_number(tree)));
    cpl_test(er_json_node_get_number(tree) > 0.0);
    er_json_node_delete(tree);

    json = " Inf ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test(isinf(er_json_node_get_number(tree)));
    cpl_test(er_json_node_get_number(tree) > 0.0);
    er_json_node_delete(tree);

    json = " Infinity ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test(isinf(er_json_node_get_number(tree)));
    cpl_test(er_json_node_get_number(tree) > 0.0);
    er_json_node_delete(tree);

    json = " -inf ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test(isinf(er_json_node_get_number(tree)));
    cpl_test(er_json_node_get_number(tree) < 0.0);
    er_json_node_delete(tree);

    json = " -Inf ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test(isinf(er_json_node_get_number(tree)));
    cpl_test(er_json_node_get_number(tree) < 0.0);
    er_json_node_delete(tree);

    json = " -Infinity ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test(isinf(er_json_node_get_number(tree)));
    cpl_test(er_json_node_get_number(tree) < 0.0);
    er_json_node_delete(tree);

    /* Make sure number overflows are handled in some reasonable manner, but
       do not cause an error. */
    json = "  1234567890123456789012345678901234567890123456789.0  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    cpl_test(er_json_node_get_number(tree) > 1.234567890123e48);
    cpl_test(er_json_node_get_number(tree) < 1.234567890124e48);
    er_json_node_delete(tree);

    json = "  1234567890123456789012345678901234567890123456789e999999999  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_NUMBER);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    cpl_test(isinf(er_json_node_get_number(tree)));
    cpl_test(er_json_node_get_number(tree) > 0.0);
    er_json_node_delete(tree);
}

static void test_parsing_string(void)
{
    /* Test parsing strings. */
    const char * json = "\"test\"";
    er_json_node * tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_STRING);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test_eq_string(er_json_node_get_string(tree), "test");
    er_json_node_delete(tree);

    json = "  \"test\"";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_STRING);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    cpl_test_eq_string(er_json_node_get_string(tree), "test");
    er_json_node_delete(tree);

    json = " \"test\"  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_STRING);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test_eq_string(er_json_node_get_string(tree), "test");
    er_json_node_delete(tree);

    json = " \"test\\n\\t\\\\\\\"\\b\\f\\r\\u0023\\/end\"  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_STRING);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test_eq_string(er_json_node_get_string(tree),
                       "test\n\t\\\"\b\f\r#/end");
    er_json_node_delete(tree);
}

static void test_parsing_array(void)
{
    /* Test parsing an array. */
    const char * json = "[]";
    er_json_node * tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_ARRAY);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test_eq(er_json_node_array_size(tree), 0);
    er_json_node_delete(tree);

    json = "[1]";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_ARRAY);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test_eq(er_json_node_array_size(tree), 1);
    er_json_array_iterator iter = er_json_node_array_begin(tree);
    /* Check the item. */
    const er_json_node * item = er_json_node_array_get(tree, iter);
    cpl_test_eq(er_json_node_type(item), JSON_NUMBER);
    cpl_test_abs(er_json_node_get_number(item), 1.0, DBL_EPSILON);
    cpl_test_eq_ptr(er_json_node_location(item), json+1);
    er_json_node_delete(tree);

    json = "  [ 1 , \" test \" ]";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_ARRAY);
    cpl_test_eq_ptr(er_json_node_location(tree), json+2);
    cpl_test_eq(er_json_node_array_size(tree), 2);
    /* Check first item. */
    iter = er_json_node_array_begin(tree);
    item = er_json_node_array_get(tree, iter);
    cpl_test_eq(er_json_node_type(item), JSON_NUMBER);
    cpl_test_abs(er_json_node_get_number(item), 1.0, DBL_EPSILON);
    cpl_test_eq_ptr(er_json_node_location(item), json+4);
    /* Check second item. */
    iter = er_json_node_array_next(tree, iter);
    item = er_json_node_array_get(tree, iter);
    cpl_test_eq(er_json_node_type(item), JSON_STRING);
    cpl_test_eq_string(er_json_node_get_string(item), " test ");
    cpl_test_eq_ptr(er_json_node_location(item), json+8);
    er_json_node_delete(tree);

    json = " [1,\" test \"]  ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_ARRAY);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test_eq(er_json_node_array_size(tree), 2);
    /* Check first item. */
    iter = er_json_node_array_begin(tree);
    item = er_json_node_array_get(tree, iter);
    cpl_test_eq(er_json_node_type(item), JSON_NUMBER);
    cpl_test_abs(er_json_node_get_number(item), 1.0, DBL_EPSILON);
    cpl_test_eq_ptr(er_json_node_location(item), json+2);
    /* Check second item. */
    iter = er_json_node_array_next(tree, iter);
    item = er_json_node_array_get(tree, iter);
    cpl_test_eq(er_json_node_type(item), JSON_STRING);
    cpl_test_eq_string(er_json_node_get_string(item), " test ");
    cpl_test_eq_ptr(er_json_node_location(item), json+4);
    er_json_node_delete(tree);
}

static void test_parsing_object(void)
{
    /* Test parsing an object. */
    const char * json = "{}";
    er_json_node * tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_OBJECT);
    cpl_test_eq_ptr(er_json_node_location(tree), json);
    cpl_test_eq(er_json_node_object_size(tree), 0);
    er_json_node_delete(tree);

    json = " {\"num\":1.2} ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_OBJECT);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test_eq(er_json_node_object_size(tree), 1);
    /* Check value. */
    const er_json_node * item = er_json_node_object_get(tree, "num");
    cpl_test_nonnull(item);
    cpl_test_eq(er_json_node_type(item), JSON_NUMBER);
    cpl_test_abs(er_json_node_get_number(item), 1.2, DBL_EPSILON);
    cpl_test_eq_ptr(er_json_node_location(item), json+8);
    er_json_node_delete(tree);

    /* Test nested structure. */
    json = " {\"array\": [1.2 , {\"str\":\"hello\\nworld\"}] } ";
    tree = er_json_parse(json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(tree);
    cpl_test_eq(er_json_node_type(tree), JSON_OBJECT);
    cpl_test_eq_ptr(er_json_node_location(tree), json+1);
    cpl_test_eq(er_json_node_object_size(tree), 1);
    /* Check value of top level object. */
    const er_json_node * array = er_json_node_object_get(tree, "array");
    cpl_test_nonnull(item);
    cpl_test_eq(er_json_node_type(array), JSON_ARRAY);
    cpl_test_eq_ptr(er_json_node_location(array), json+11);
    cpl_test_eq(er_json_node_array_size(array), 2);
    /* Check the first array item. */
    er_json_array_iterator iter = er_json_node_array_begin(array);
    cpl_test_nonnull(iter);
    item = er_json_node_array_get(array, iter);
    cpl_test_nonnull(item);
    cpl_test_eq(er_json_node_type(item), JSON_NUMBER);
    cpl_test_abs(er_json_node_get_number(item), 1.2, DBL_EPSILON);
    cpl_test_eq_ptr(er_json_node_location(item), json+12);
    /* Check the second array item. */
    iter = er_json_node_array_next(array, iter);
    cpl_test_nonnull(iter);
    item = er_json_node_array_get(array, iter);
    cpl_test_nonnull(item);
    cpl_test_eq(er_json_node_type(item), JSON_OBJECT);
    cpl_test_eq_ptr(er_json_node_location(item), json+18);
    cpl_test_eq(er_json_node_object_size(item), 1);
    /* Check the second array item's entry. */
    item = er_json_node_object_get(item, "str");
    cpl_test_nonnull(item);
    cpl_test_eq(er_json_node_type(item), JSON_STRING);
    cpl_test_eq_ptr(er_json_node_location(item), json+25);
    cpl_test_eq_string(er_json_node_get_string(item), "hello\nworld");
    er_json_node_delete(tree);
}

static void test_parse_errors(void)
{
    /* Test handling of invalid JSON input. */
    cpl_test_null(er_json_parse("random"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("\""));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("["));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[1"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[1,"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[,"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{1:2}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"1\""));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"1\":"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"1\":}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{:123}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"1\":123,}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"1\":123"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"1\":123,"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("\"\\"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("\"\\x\""));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("\"\\u\""));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("\"\\u00\""));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("\"\\u30      \""));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("\"\\u7FFF\""));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("\"\\u0000\""));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("nul"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("na"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("nana"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("-"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("Infin"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("-infi"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("tru"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("truee"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("fal"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("falses"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("123x4"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[nul, \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[tru, \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[fal, \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[123x4, \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[\"value, \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[\"\\u30    \", \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{ : \"check with long input string . . . . ."
                                " . . . . . . . . . . . . . . . . .\"}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{12345: \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"key\\u13\": \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"key\" \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"}"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("{\"key\": 123 check with long input string"
                                " . . . . . . . . . . . . . . . . ."));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_null(er_json_parse("[random, \"check with long input string"
                                " . . . . . . . . . . . . . . . . .\"]"));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
}
