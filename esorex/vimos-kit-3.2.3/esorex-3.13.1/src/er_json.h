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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifndef ER_JSON_H
#define ER_JSON_H

#include <cpl.h>
#include "er_stringarray.h"

CPL_BEGIN_DECLS

/*
 * Functions for converting CPL objects to JSON.
 */

cpl_error_code er_frameset_to_json(const cpl_frameset *frameset,
                                   const char *json_file);
cpl_error_code er_frameset_to_text(const cpl_frameset *frameset,
                                   const char *text_file);
cpl_error_code er_recipe_parameterlist_to_json(const cpl_parameterlist *plist,
                                               const char *recipe_name,
                                               const char *json_file);
char *er_json_escape_string(const char *str);
char *er_plugin_to_json(const cpl_plugin *plugin);

/*
 * JSON parsing function and parse tree navigation functions.
 */

typedef enum _er_json_type_
{
    JSON_NULL,
    JSON_OBJECT,
    JSON_ARRAY,
    JSON_STRING,
    JSON_NUMBER,
    JSON_BOOL
} er_json_type;

typedef struct _er_json_node_ er_json_node;
typedef const struct _cx_lnode_ * er_json_array_iterator;
typedef const struct _cx_tnode_ * er_json_object_iterator;

er_json_node * er_json_parse(const char * json);

void er_json_node_delete(er_json_node * obj);
er_json_type er_json_node_type(const er_json_node * obj);
const char * er_json_node_location(const er_json_node * obj);
cpl_boolean er_json_node_get_bool(const er_json_node * obj);
double er_json_node_get_number(const er_json_node * obj);
const char * er_json_node_get_string(const er_json_node * obj);
cpl_size er_json_node_array_size(const er_json_node * obj);
cpl_boolean er_json_node_array_empty(const er_json_node * obj);
er_json_array_iterator er_json_node_array_begin(const er_json_node * obj);
er_json_array_iterator er_json_node_array_end(const er_json_node * obj);
er_json_array_iterator er_json_node_array_next(const er_json_node * obj,
                                               er_json_array_iterator current);
er_json_array_iterator er_json_node_array_previous(
                                                const er_json_node * obj,
                                                er_json_array_iterator current);
const er_json_node * er_json_node_array_get(const er_json_node * obj,
                                            er_json_array_iterator iterator);
cpl_size er_json_node_object_size(const er_json_node * obj);
cpl_boolean er_json_node_object_empty(const er_json_node * obj);
er_json_object_iterator er_json_node_object_begin(const er_json_node * obj);
er_json_object_iterator er_json_node_object_end(const er_json_node * obj);
er_json_object_iterator er_json_node_object_next(
                                            const er_json_node * obj,
                                            er_json_object_iterator current);
er_json_object_iterator er_json_node_object_previous(
                                            const er_json_node * obj,
                                            er_json_object_iterator current);
const char * er_json_node_object_get_key(const er_json_node * obj,
                                         er_json_object_iterator iterator);
const er_json_node * er_json_node_object_get_value(
                                            const er_json_node * obj,
                                            er_json_object_iterator iterator);
const er_json_node * er_json_node_object_get(const er_json_node * obj,
                                             const char * key);

cpl_error_code er_json_find_line_column(const char * json,
                                        const char * location,
                                        int * line, int * column);

/*
 * Utility functions for converting a parsed JSON node into actual objects.
 */

er_stringarray_t * er_json_to_string_array(const er_json_node * parsetree,
                                           const char * json);

#ifdef ENABLE_PYTHON_RECIPES

cpl_plugin * er_json_to_plugin(const er_json_node * parsetree,
                               const char * json);
cpl_pluginlist * er_json_to_pluginlist(const er_json_node * parsetree,
                                       const char * json);
void er_json_pluginlist_delete(cpl_pluginlist * list);

#endif

CPL_END_DECLS

#endif /* ER_JSON_H */
