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

#include <string.h>
#include <float.h>

#include "er_json.h"
#include "esorex_test.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

#ifndef STR_LENGTH    
#define STR_LENGTH 80
#endif

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void esorex_json_parameterlist(void);
static void esorex_json_frameset(void);
static void esorex_json_to_string_array(void);
static void esorex_json_escape_string(void);

#ifdef ENABLE_PYTHON_RECIPES

static void esorex_json_plugin_conversion(void);
static void esorex_test_plugin_parsing(void);
static void esorex_test_pluginlist_parsing(void);

#endif

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/
int main(void)
{

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    esorex_json_parameterlist();
    esorex_json_frameset();
    esorex_json_to_string_array();
    esorex_json_escape_string();

#ifdef ENABLE_PYTHON_RECIPES

    esorex_json_plugin_conversion();
    esorex_test_plugin_parsing();
    esorex_test_pluginlist_parsing();

#endif

    /* End of tests */
    return cpl_test_end(0);
}

static void esorex_json_parameterlist(void)
{
    cpl_parameter     * par1;
    cpl_parameter     * par2;
    cpl_parameter     * par3;
    cpl_parameter     * par4;
    cpl_parameter     * par5;
    cpl_parameter     * par6;
    cpl_parameterlist * parlist;
    
    /* Export a parameterlist to a JSON file */

    par1 = cpl_parameter_new_value("recipe.par1", CPL_TYPE_DOUBLE,
                                   "Description 1","rec1", 10.0);
    
    cpl_test_nonnull( par1 );

    par2 = cpl_parameter_new_value("recipe.par2", CPL_TYPE_STRING,
                                   "Description 2","recipe", "String");
    
    cpl_test_nonnull( par2 );
    
    par3 = cpl_parameter_new_value("recipe.par3", CPL_TYPE_BOOL,
                                   "Description 3","recipe", CPL_TRUE);
    
    cpl_test_nonnull( par3 );

    par4 = cpl_parameter_new_value("recipe.par4", CPL_TYPE_INT,
                                   "Escaped chars \" / \f \n \r \t \\","recipe", -10);
    
    par5 = cpl_parameter_new_enum("recipe.par5", CPL_TYPE_STRING,
                                  "Description 5","recipe", "Str1", 
                                  3,"Str1", "Str2", "Str3");
    
    par6 = cpl_parameter_new_range("recipe.par6", CPL_TYPE_INT,
                                  "Description 6","recipe", 0, -1, 1);

    cpl_test_nonnull( par4 );

    parlist = cpl_parameterlist_new();
    
    cpl_test_nonnull( parlist );

    cpl_parameterlist_append(parlist, par1);
    
    /* Dumping to json */
    er_recipe_parameterlist_to_json(parlist, "recipe_name", "parlist.json");
    
    cpl_test_error(CPL_ERROR_NONE);

    /* Read the file and compare it with a static version */
    {
        char * whole_file = 0;
        long length;
        FILE * fjson = fopen("parlist.json", "rb");
        const char * expected_json="[\n"
                "  {\n"
                "    \"name\": \"recipe.par1\",\n"
                "    \"value\": 10,\n"
                "    \"valtype\": \"double\",\n"
                "    \"partype\": \"value\",\n"
                "    \"display_name\": \"recipe.par1\",\n"
                "    \"description\": \"Description 1\",\n"
                "    \"recipe\": \"recipe_name\"\n"
                "  }\n"
                "]\n\0";
        
        cpl_test_nonnull(fjson);
        fseek(fjson, 0, SEEK_END);
        length = ftell(fjson);
        rewind(fjson);
        whole_file = cpl_malloc(length+1);
        
        cpl_test_nonnull(whole_file);
        fread(whole_file, 1, length, fjson);
        whole_file[length]='\0';
        fclose(fjson);
        
        cpl_test_zero( strcmp(whole_file, expected_json) );
        
        cpl_test_zero( remove("parlist.json") );
        cpl_free(whole_file);
    }

    /* Now with no recipe name */
    er_recipe_parameterlist_to_json(parlist, "", "parlist11.json");
    
    cpl_test_error(CPL_ERROR_NONE);

    /* Read the file and compare it with a static version */
    {
        char * whole_file = 0;
        long length;
        FILE * fjson = fopen("parlist11.json", "rb");
        const char * expected_json="[\n"
                "  {\n"
                "    \"name\": \"recipe.par1\",\n"
                "    \"value\": 10,\n"
                "    \"valtype\": \"double\",\n"
                "    \"partype\": \"value\",\n"
                "    \"display_name\": \"recipe.par1\",\n"
                "    \"description\": \"Description 1\"\n"
                "  }\n"
                "]\n\0";
        
        cpl_test_nonnull(fjson);
        fseek(fjson, 0, SEEK_END);
        length = ftell(fjson);
        rewind(fjson);
        whole_file = cpl_malloc(length+1);
        
        cpl_test_nonnull(whole_file);
        fread(whole_file, 1, length, fjson);
        whole_file[length]='\0';
        fclose(fjson);
        
        cpl_test_zero( strcmp(whole_file, expected_json) );
        
        cpl_test_zero( remove("parlist11.json") );
        cpl_free(whole_file);
    }
    
    /* Do the same with a list that contains more than 1 parameter, including
     * a parameter with control characters */
    cpl_parameterlist_append(parlist, par2);
    cpl_parameterlist_append(parlist, par3);
    cpl_parameterlist_append(parlist, par4);
    cpl_parameterlist_append(parlist, par5);
    cpl_parameterlist_append(parlist, par6);
    
    cpl_test_error(CPL_ERROR_NONE);
            
    /* Dumping to json */
    er_recipe_parameterlist_to_json(parlist, "recipe_name", "parlist2.json");

    cpl_test_error(CPL_ERROR_NONE);

    /* Read the file and compare it with a static version */
    {
        char * whole_file = 0;
        long length;
        FILE * fjson = fopen("parlist2.json", "rb");
        const char * expected_json="[\n"
                "  {\n"
                "    \"name\": \"recipe.par1\",\n"
                "    \"value\": 10,\n"
                "    \"valtype\": \"double\",\n"
                "    \"partype\": \"value\",\n"
                "    \"display_name\": \"recipe.par1\",\n"
                "    \"description\": \"Description 1\",\n"
                "    \"recipe\": \"recipe_name\"\n"
                "  },\n"
                "  {\n"
                "    \"name\": \"recipe.par2\",\n"
                "    \"value\": \"String\",\n"
                "    \"valtype\": \"string\",\n"
                "    \"partype\": \"value\",\n"
                "    \"display_name\": \"recipe.par2\",\n"
                "    \"description\": \"Description 2\",\n"
                "    \"recipe\": \"recipe_name\"\n"
                "  },\n"
                "  {\n"
                "    \"name\": \"recipe.par3\",\n"
                "    \"value\": true,\n"
                "    \"valtype\": \"bool\",\n"
                "    \"partype\": \"value\",\n"
                "    \"display_name\": \"recipe.par3\",\n"
                "    \"description\": \"Description 3\",\n"
                "    \"recipe\": \"recipe_name\"\n"
                "  },\n"
                "  {\n"
                "    \"name\": \"recipe.par4\",\n"
                "    \"value\": -10,\n"
                "    \"valtype\": \"int\",\n"
                "    \"partype\": \"value\",\n"
                "    \"display_name\": \"recipe.par4\",\n"
                "    \"description\": \"Escaped chars \\\" \\/ \\f \\n \\r \\t \\\\\",\n"
                "    \"recipe\": \"recipe_name\"\n"
                "  },\n"
                "  {\n"
                "    \"name\": \"recipe.par5\",\n"
                "    \"value\": \"Str1\",\n"
                "    \"valtype\": \"string\",\n"
                "    \"partype\": \"enum\",\n"
                "    \"valenum\": [ \"Str1\" , \"Str2\" , \"Str3\" ],\n"
                "    \"display_name\": \"recipe.par5\",\n"
                "    \"description\": \"Description 5\",\n"
                "    \"recipe\": \"recipe_name\"\n"
                "  },\n"
                "  {\n"
                "    \"name\": \"recipe.par6\",\n"
                "    \"value\": 0,\n"
                "    \"valtype\": \"int\",\n"
                "    \"partype\": \"range\",\n"
                "    \"valmin\": -1,\n"
                "    \"valmax\": 1,\n"
                "    \"display_name\": \"recipe.par6\",\n"
                "    \"description\": \"Description 6\",\n"
                "    \"recipe\": \"recipe_name\"\n"
                "  }\n"
                "]\n\0";

        cpl_test_nonnull(fjson);
        fseek(fjson, 0, SEEK_END);
        length = ftell(fjson);
        rewind(fjson);
        whole_file = cpl_malloc(length+1);

        cpl_test_nonnull(whole_file);
        fread(whole_file, 1, length, fjson);
        whole_file[length]='\0';
        fclose(fjson);

        cpl_test_zero( strcmp(whole_file, expected_json) );

        cpl_test_zero( remove("parlist2.json") );
        cpl_free(whole_file);
    }

    /* Try to save in a non-existent path */
    {
        cpl_error_code errcode = 
                er_recipe_parameterlist_to_json(parlist, "recipe_name",
                        "/non-existent/path/parlist.json");
        cpl_test_error(CPL_ERROR_FILE_NOT_CREATED);
        cpl_test_eq(errcode , CPL_ERROR_FILE_NOT_CREATED);
    }
    
    /* Trying with null inputs */
    {
        cpl_error_code errcode = 
                er_recipe_parameterlist_to_json(NULL, "recipe_name",
                        "parlist.json");
        cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
        cpl_test_eq(errcode , CPL_ERROR_ILLEGAL_INPUT);
        
        errcode =  er_recipe_parameterlist_to_json(parlist, NULL,
                         "parlist.json");
        cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
        cpl_test_eq(errcode , CPL_ERROR_ILLEGAL_INPUT);

        errcode =  er_recipe_parameterlist_to_json(parlist, "recipe_name",
                         NULL);
        cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
        cpl_test_eq(errcode , CPL_ERROR_ILLEGAL_INPUT);
    }
    
    cpl_parameterlist_delete(parlist);
}

static void esorex_json_frameset(void)
{
    cpl_frame    * fset1;
    cpl_frame    * fset2;
    cpl_frame    * fset3;
    cpl_frameset * frameset;
    
    /* Export a frameset to a JSON file */

    fset1 = cpl_frame_new();
    cpl_frame_set_filename(fset1, "filename1.fits");
    cpl_frame_set_tag(fset1, "FLAT");
    
    cpl_test_nonnull( fset1 );

    fset2 = cpl_frame_new();
    cpl_frame_set_filename(fset2, "filename2.fits");
    cpl_frame_set_tag(fset2, "BIAS");
    
    cpl_test_nonnull( fset2 );
    
    fset3 = cpl_frame_new();
    cpl_frame_set_filename(fset3, "filename3.fits");
    cpl_frame_set_tag(fset3, "SCIENCE");
    
    cpl_test_nonnull( fset3 );

    frameset = cpl_frameset_new();
    
    cpl_test_nonnull( frameset );

    cpl_frameset_insert(frameset, fset1);
    
    /* Dumping to json */
    er_frameset_to_json(frameset, "frameset.json");
    
    cpl_test_error(CPL_ERROR_NONE);

    /* Read the file and compare it with a static version */
    {
        char * whole_file = 0;
        long length;
        FILE * fjson = fopen("frameset.json", "rb");
        const char * expected_json="[\n"
                "  {\n"
                "    \"name\": \"filename1.fits\",\n"
                "    \"category\": \"FLAT\"\n"
                "  }\n"
                "]\n\0";
        
        cpl_test_nonnull(fjson);
        fseek(fjson, 0, SEEK_END);
        length = ftell(fjson);
        rewind(fjson);
        whole_file = cpl_malloc(length+1);
        
        cpl_test_nonnull(whole_file);
        fread(whole_file, 1, length, fjson);
        whole_file[length]='\0';
        fclose(fjson);
        
        cpl_test_zero( strcmp(whole_file, expected_json) );
        
        cpl_test_zero( remove("frameset.json") );
        cpl_free(whole_file);
    }

    /* Do the same with a list that contains more than 1 frame */
    cpl_frameset_insert(frameset, fset2);
    cpl_frameset_insert(frameset, fset3);
    
    cpl_test_error(CPL_ERROR_NONE);
            
    /* Dumping to json */
    er_frameset_to_json(frameset, "frameset2.json");

    cpl_test_error(CPL_ERROR_NONE);
    /* Read the file and compare it with a static version */
    {
        char * whole_file = 0;
        long length;
        FILE * fjson = fopen("frameset2.json", "rb");
        const char * expected_json="[\n"
                "  {\n"
                "    \"name\": \"filename1.fits\",\n"
                "    \"category\": \"FLAT\"\n"
                "  },\n"
                "  {\n"
                "    \"name\": \"filename2.fits\",\n"
                "    \"category\": \"BIAS\"\n"
                "  },\n"
                "  {\n"
                "    \"name\": \"filename3.fits\",\n"
                "    \"category\": \"SCIENCE\"\n"
                "  }\n"
                "]\n\0";

        cpl_test_nonnull(fjson);
        fseek(fjson, 0, SEEK_END);
        length = ftell(fjson);
        rewind(fjson);
        whole_file = cpl_malloc(length+1);

        cpl_test_nonnull(whole_file);
        fread(whole_file, 1, length, fjson);
        whole_file[length]='\0';
        fclose(fjson);

        cpl_test_zero( strcmp(whole_file, expected_json) );

        cpl_test_zero( remove("frameset2.json") );
        cpl_free(whole_file);
    }

    /* Try to save in a non-existent path */
    {
        cpl_error_code errcode =
                er_frameset_to_json(frameset, 
                        "/non-existent/path/frameset.json");
        cpl_test_error(CPL_ERROR_FILE_NOT_CREATED);
        cpl_test_eq(errcode , CPL_ERROR_FILE_NOT_CREATED);
    }
    
    /* Trying with null inputs */
    {
        cpl_error_code errcode = 
                er_frameset_to_json(NULL, "frameset.json");
        cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
        cpl_test_eq(errcode , CPL_ERROR_ILLEGAL_INPUT);
        
        errcode = er_frameset_to_json(frameset, NULL);
        cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
        cpl_test_eq(errcode , CPL_ERROR_ILLEGAL_INPUT);
    }
    
    cpl_frameset_delete(frameset);
}

static void esorex_json_to_string_array(void)
{
    /* Test conversion of a parsed JSON array of strings to er_stringarray_t. */
    const char * json = "[\"foo\",\"bar\"]";
    er_json_node * tree = er_json_parse(json);
    cpl_test_assert(tree != NULL);
    er_stringarray_t * array = er_json_to_string_array(tree, json);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(array);
    cpl_test_eq(er_stringarray_size(array), 2);
    cpl_test_eq_string(er_stringarray_get(array, 0), "foo");
    cpl_test_eq_string(er_stringarray_get(array, 1), "bar");

    /* Check handling NULL inputs. */
    cpl_test_null(er_json_to_string_array(NULL, json));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(er_json_to_string_array(tree, NULL));
    cpl_test_error(CPL_ERROR_NULL_INPUT);

    er_stringarray_delete(array);
    er_json_node_delete(tree);

    /* Check handling of parsing error conditions. */
    json = "{\"foo\":\"bar\"}";
    tree = er_json_parse(json);
    cpl_test_assert(tree != NULL);
    cpl_test_null(er_json_to_string_array(tree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(tree);

    json = "[\"foo\",123]";
    tree = er_json_parse(json);
    cpl_test_assert(tree != NULL);
    cpl_test_null(er_json_to_string_array(tree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(tree);
}

static void esorex_json_escape_string(void)
{
    /* Test the er_json_escape_string function. */
    char * result = er_json_escape_string(" / \t \n \r \f \b \" \\ ");
    cpl_test_nonnull(result);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq_string(result, " \\/ \\t \\n \\r \\f \\b \\\" \\\\ ");
    cpl_free(result);

    /* Test error handling. */
    cpl_test_null(er_json_escape_string(NULL));
    cpl_test_error(CPL_ERROR_NONE);
}

#ifdef ENABLE_PYTHON_RECIPES

static void esorex_json_plugin_conversion(void)
{
    cpl_plugin * plugin = NULL;
    cpl_plugin * parsed_plugin = NULL;
    cpl_recipe * recipe = NULL;
    cpl_recipe * parsed_recipe = NULL;
    cpl_recipe2 * recipe2 = NULL;
    cpl_recipe2 * parsed_recipe2 = NULL;
    cpl_parameter * param1 = NULL;
    cpl_parameter * param2 = NULL;
    cpl_parameter * param3 = NULL;
    cpl_parameter * param4 = NULL;
    cpl_parameter * param5 = NULL;
    cpl_parameter * param6 = NULL;
    cpl_parameter * param7 = NULL;
    cpl_parameter * param8 = NULL;
    cpl_parameter * param9 = NULL;
    cpl_parameter * param10 = NULL;
    cpl_parameter * parsed_param0 = NULL;
    cpl_parameter * parsed_param1 = NULL;
    cpl_parameter * parsed_param2 = NULL;
    cpl_parameter * parsed_param3 = NULL;
    cpl_parameter * parsed_param4 = NULL;
    cpl_parameter * parsed_param5 = NULL;
    cpl_parameter * parsed_param6 = NULL;
    cpl_parameter * parsed_param7 = NULL;
    cpl_parameter * parsed_param8 = NULL;
    cpl_parameter * parsed_param9 = NULL;
    cpl_parameter * parsed_param10 = NULL;
    cpl_frame * frame1 = NULL;
    cpl_frame * frame2 = NULL;
    cpl_frame * parsed_frame1 = NULL;
    cpl_frame * parsed_frame2 = NULL;
    char * json = NULL;
    er_json_node * parsetree = NULL;
    char** tags = NULL;
    char** inputs = NULL;
    char** outputs = NULL;
    char** charptr = NULL;

    /* Test some basic API error handling for the er_plugin_to_json and
       er_json_to_plugin functions. */
    cpl_test_null(er_plugin_to_json(NULL));
    cpl_test_error(CPL_ERROR_NULL_INPUT);

    parsetree = er_json_parse("null");
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(NULL, NULL));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(er_json_to_plugin(parsetree, NULL));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(er_json_to_plugin(NULL, "null"));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    er_json_node_delete(parsetree);

    /* Test for error if an invalid plugin type is given. */
    plugin = cpl_plugin_new();
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_NONE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_UNSUPPORTED_MODE);
    cpl_test_null(json);
    cpl_plugin_delete(plugin);

    /* Construct a simple recipe object, convert it to JSON text and check that
       the set strings are present in the output. */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
    plugin = &recipe->interface;
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(json);
    cpl_test(strstr(json, "\"test_recipe\"") != NULL);
    cpl_test(strstr(json, "\"simple test recipe\"") != NULL);
    cpl_test(strstr(json, "\"description ...\"") != NULL);
    cpl_test(strstr(json, "\"some author\"") != NULL);
    cpl_test(strstr(json, "\"author@eso.org\"") != NULL);
    cpl_test(strstr(json, "\"copyright ...\"") != NULL);

    /* Parse the produced JSON back to a CPL plugin object and check that
       the structure is correct. */
    parsetree = er_json_parse(json);
    cpl_test_nonnull(parsetree);
    cpl_test_error(CPL_ERROR_NONE);
    parsed_plugin = er_json_to_plugin(parsetree, json);
    cpl_test_nonnull(parsed_plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(cpl_plugin_get_version(parsed_plugin), 1);
    cpl_test_eq(cpl_plugin_get_type(parsed_plugin), CPL_PLUGIN_TYPE_RECIPE);
    cpl_test_eq_string(cpl_plugin_get_name(parsed_plugin), "test_recipe");
    cpl_test_eq_string(cpl_plugin_get_synopsis(parsed_plugin),
                       "simple test recipe");
    cpl_test_eq_string(cpl_plugin_get_description(parsed_plugin),
                       "description ...");
    cpl_test_eq_string(cpl_plugin_get_author(parsed_plugin), "some author");
    cpl_test_eq_string(cpl_plugin_get_email(parsed_plugin), "author@eso.org");
    cpl_test_eq_string(cpl_plugin_get_copyright(parsed_plugin),
                       "copyright ...");
    cpl_free(json);
    cpl_plugin_delete(plugin);
    parsed_plugin->deinitialize(parsed_plugin);
    cpl_plugin_delete(parsed_plugin);
    er_json_node_delete(parsetree);

    /* Construct a recipe v2 object, convert it to JSON text and check that
       the set strings are present in the output. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    plugin = &recipe2->base.interface;
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE_V2, "test_recipe2",
                                    "simple test recipe v2", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(json);
    cpl_test(strstr(json, "\"test_recipe2\"") != NULL);
    cpl_test(strstr(json, "\"simple test recipe v2\"") != NULL);
    cpl_test(strstr(json, "\"description ...\"") != NULL);
    cpl_test(strstr(json, "\"some author\"") != NULL);
    cpl_test(strstr(json, "\"author@eso.org\"") != NULL);
    cpl_test(strstr(json, "\"copyright ...\"") != NULL);
    cpl_free(json);
    cpl_plugin_delete(plugin);

    /* Test a more complex example of a recipe plugin with parameters and an
       input frameset. */
    recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));

    recipe->parameters = cpl_parameterlist_new();
    cpl_test_assert(recipe->parameters != NULL);
    param1 = cpl_parameter_new_value("test.par1", CPL_TYPE_BOOL,
                                     "bool value", "test_1", CPL_TRUE);
    cpl_test_assert(param1 != NULL);
    cpl_parameterlist_append(recipe->parameters, param1);
    param2 = cpl_parameter_new_value("test.par2", CPL_TYPE_INT,
                                     "int value", "test_2", 2);
    cpl_test_assert(param2 != NULL);
    cpl_parameterlist_append(recipe->parameters, param2);
    param3 = cpl_parameter_new_value("test.par3", CPL_TYPE_DOUBLE,
                                     "float value", "test_3", 4.0);
    cpl_test_assert(param3 != NULL);
    cpl_parameterlist_append(recipe->parameters, param3);
    param4 = cpl_parameter_new_value("test.par4", CPL_TYPE_STRING,
                                     "string value", "test_4", "something");
    cpl_test_assert(param4 != NULL);
    cpl_parameterlist_append(recipe->parameters, param4);
    param5 = cpl_parameter_new_range("test.par5", CPL_TYPE_INT,
                                     "range value1", "test_5", 3, 1, 5);
    cpl_test_assert(param5 != NULL);
    cpl_parameterlist_append(recipe->parameters, param5);
    param6 = cpl_parameter_new_range("test.par6", CPL_TYPE_DOUBLE,
                                     "range value2", "test_6", 3.0, 1.0, 5.0);
    cpl_test_assert(param6 != NULL);
    cpl_parameterlist_append(recipe->parameters, param6);
    param7 = cpl_parameter_new_enum("test.par7", CPL_TYPE_INT, "enum value1",
                                    "test_7", 2, 3, 1, 2, 4);
    cpl_test_assert(param7 != NULL);
    cpl_parameterlist_append(recipe->parameters, param7);
    param8 = cpl_parameter_new_enum("test.par8", CPL_TYPE_DOUBLE, "enum value2",
                                    "test_8", 4.4, 3, 1.1, 2.2, 4.4);
    cpl_test_assert(param8 != NULL);
    cpl_parameterlist_append(recipe->parameters, param8);
    param9 = cpl_parameter_new_enum("test.par9", CPL_TYPE_STRING, "enum value3",
                                    "test_9", "y", 3, "x", "y", "z");
    cpl_test_assert(param9 != NULL);
    cpl_parameter_set_default_flag(param9, 1);
    cpl_parameter_set_alias(param9, CPL_PARAMETER_MODE_CLI, "par9_cli");
    cpl_parameter_set_alias(param9, CPL_PARAMETER_MODE_ENV, "par9_env");
    cpl_parameter_set_alias(param9, CPL_PARAMETER_MODE_CFG, "par9_cfg");
    cpl_parameterlist_append(recipe->parameters, param9);
    param10 = cpl_parameter_new_value("test.par10", CPL_TYPE_BOOL, "disabled",
                                     "test_10", CPL_FALSE);
    cpl_test_assert(param10 != NULL);
    cpl_parameter_set_default_flag(param10, 0);
    cpl_parameter_disable(param10, CPL_PARAMETER_MODE_CLI);
    cpl_parameter_disable(param10, CPL_PARAMETER_MODE_ENV);
    cpl_parameter_disable(param10, CPL_PARAMETER_MODE_CFG);
    cpl_parameterlist_append(recipe->parameters, param10);

    recipe->frames = cpl_frameset_new();
    cpl_test_assert(recipe->frames != NULL);
    frame1 = cpl_frame_new();
    cpl_test_assert(frame1 != NULL);
    cpl_frame_set_filename(frame1, "test1.fits");
    cpl_frame_set_tag(frame1, "RAW");
    cpl_frame_set_type(frame1, CPL_FRAME_TYPE_IMAGE);
    cpl_frame_set_group(frame1, CPL_FRAME_GROUP_RAW);
    cpl_frame_set_level(frame1, CPL_FRAME_LEVEL_FINAL);
    cpl_frameset_insert(recipe->frames, frame1);
    frame2 = cpl_frame_new();
    cpl_test_assert(frame2 != NULL);
    cpl_frame_set_filename(frame2, "test2.fits");
    cpl_frame_set_tag(frame2, "CALIB");
    cpl_frameset_insert(recipe->frames, frame2);

    plugin = &recipe->interface;
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE, "test_recipe",
                                    "simple test recipe", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(json);
    cpl_test(strstr(json, "\"test_recipe\"") != NULL);
    cpl_test(strstr(json, "\"simple test recipe\"") != NULL);
    cpl_test(strstr(json, "\"description ...\"") != NULL);
    cpl_test(strstr(json, "\"some author\"") != NULL);
    cpl_test(strstr(json, "\"author@eso.org\"") != NULL);
    cpl_test(strstr(json, "\"copyright ...\"") != NULL);

    /* Parse the produced JSON back into a CPL plugin object and check that
       the recipe structure is correct, including the parameters and frameset.
     */
    parsetree = er_json_parse(json);
    cpl_test_nonnull(parsetree);
    cpl_test_error(CPL_ERROR_NONE);
    parsed_plugin = er_json_to_plugin(parsetree, json);
    parsed_recipe = (cpl_recipe *) parsed_plugin;
    cpl_test_nonnull(parsed_plugin);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(cpl_plugin_get_version(parsed_plugin), 1);
    cpl_test_eq(cpl_plugin_get_type(parsed_plugin), CPL_PLUGIN_TYPE_RECIPE);
    cpl_test_eq_string(cpl_plugin_get_name(parsed_plugin), "test_recipe");
    cpl_test_eq_string(cpl_plugin_get_synopsis(parsed_plugin),
                       "simple test recipe");
    cpl_test_eq_string(cpl_plugin_get_description(parsed_plugin),
                       "description ...");
    cpl_test_eq_string(cpl_plugin_get_author(parsed_plugin), "some author");
    cpl_test_eq_string(cpl_plugin_get_email(parsed_plugin), "author@eso.org");
    cpl_test_eq_string(cpl_plugin_get_copyright(parsed_plugin),
                       "copyright ...");

    /* Check that there are 11 parameters in the parameters list. One of them
       must be a specially added __class__ parameters. The rest should be the
       same parameters as declared in the original plugin object that was
       serialised. */
    cpl_test_eq(cpl_parameterlist_get_size(parsed_recipe->parameters), 11);
    parsed_param0 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "__class__");
    cpl_test_nonnull(parsed_param0);
    cpl_test_eq(cpl_parameter_get_type(parsed_param0), CPL_TYPE_STRING);
    cpl_test_eq_string(cpl_parameter_get_string(parsed_param0), "unknown");
    cpl_test_eq_string(cpl_parameter_get_default_string(parsed_param0),
                       "unknown");
    parsed_param1 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par1");
    cpl_test_nonnull(parsed_param1);
    test_equal_param(param1, parsed_param1);
    parsed_param2 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par2");
    cpl_test_nonnull(parsed_param2);
    test_equal_param(param2, parsed_param2);
    parsed_param3 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par3");
    cpl_test_nonnull(parsed_param3);
    test_equal_param(param3, parsed_param3);
    parsed_param4 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par4");
    cpl_test_nonnull(parsed_param4);
    test_equal_param(param4, parsed_param4);
    parsed_param5 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par5");
    cpl_test_nonnull(parsed_param5);
    test_equal_param(param5, parsed_param5);
    parsed_param6 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par6");
    cpl_test_nonnull(parsed_param6);
    test_equal_param(param6, parsed_param6);
    parsed_param7 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par7");
    cpl_test_nonnull(parsed_param7);
    test_equal_param(param7, parsed_param7);
    parsed_param8 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par8");
    cpl_test_nonnull(parsed_param8);
    test_equal_param(param8, parsed_param8);
    parsed_param9 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par9");
    cpl_test_nonnull(parsed_param9);
    test_equal_param(param9, parsed_param9);
    parsed_param10 = cpl_parameterlist_find(parsed_recipe->parameters,
                                           "test.par10");
    cpl_test_nonnull(parsed_param10);
    test_equal_param(param10, parsed_param10);

    /* Check the two parsed frames are the same as the original. */
    cpl_test_eq(cpl_frameset_get_size(parsed_recipe->frames), 2);
    parsed_frame1 = cpl_frameset_get_position(parsed_recipe->frames, 0);
    cpl_test_nonnull(parsed_frame1);
    test_equal_frame(frame1, parsed_frame1);
    parsed_frame2 = cpl_frameset_get_position(parsed_recipe->frames, 1);
    cpl_test_nonnull(parsed_frame2);
    test_equal_frame(frame2, parsed_frame2);

    cpl_free(json);
    cpl_parameterlist_delete(recipe->parameters);
    cpl_frameset_delete(recipe->frames);
    cpl_plugin_delete(plugin);
    parsed_plugin->deinitialize(parsed_plugin);
    cpl_plugin_delete(parsed_plugin);
    er_json_node_delete(parsetree);

    /* Construct a more complex example of a recipe v2 object, including its
       recipeconfig object. We then try produce JSON out of it and parse it
       back, to see that we get the same structure. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));

    recipe2->base.parameters = cpl_parameterlist_new();
    cpl_test_assert(recipe2->base.parameters != NULL);
    param1 = cpl_parameter_new_value("test.par1", CPL_TYPE_INT,
                                     "int value", "test", 123);
    cpl_test_assert(param1 != NULL);
    cpl_parameterlist_append(recipe2->base.parameters, param1);

    recipe2->base.frames = cpl_frameset_new();
    cpl_test_assert(recipe2->base.frames != NULL);
    frame1 = cpl_frame_new();
    cpl_frame_set_filename(frame1, "test.fits");
    cpl_frame_set_tag(frame1, "RAW");
    cpl_test_assert(frame1 != NULL);
    cpl_frameset_insert(recipe2->base.frames, frame1);

    recipe2->config = cpl_recipeconfig_new();
    cpl_test_assert(recipe2->config != NULL);
    cpl_recipeconfig_set_tag(recipe2->config, "RAW", 1, 1);
    cpl_recipeconfig_set_input(recipe2->config, "RAW", "CALIB", 2, 3);
    cpl_recipeconfig_set_output(recipe2->config, "RAW", "PROD");
    cpl_recipeconfig_set_tag(recipe2->config, "RAW2", 1, 1);
    cpl_recipeconfig_set_input(recipe2->config, "RAW2", "CALIB", 2, 2);
    cpl_recipeconfig_set_input(recipe2->config, "RAW2", "DARK", 1, 1);
    cpl_recipeconfig_set_output(recipe2->config, "RAW2", "PROD2");
    cpl_recipeconfig_set_output(recipe2->config, "RAW2", "PROD3");
    cpl_test_error(CPL_ERROR_NONE);

    plugin = &recipe2->base.interface;
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE_V2, "test_recipe2",
                                    "simple test recipe v2", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(json);

    parsetree = er_json_parse(json);
    cpl_test_nonnull(parsetree);
    cpl_test_error(CPL_ERROR_NONE);
    parsed_plugin = er_json_to_plugin(parsetree, json);
    parsed_recipe2 = (cpl_recipe2 *) parsed_plugin;
    cpl_test_nonnull(parsed_plugin);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(cpl_plugin_get_version(parsed_plugin), 1);
    cpl_test_eq(cpl_plugin_get_type(parsed_plugin), CPL_PLUGIN_TYPE_RECIPE_V2);
    cpl_test_eq_string(cpl_plugin_get_name(parsed_plugin), "test_recipe2");
    cpl_test_eq_string(cpl_plugin_get_synopsis(parsed_plugin),
                       "simple test recipe v2");
    cpl_test_eq_string(cpl_plugin_get_description(parsed_plugin),
                       "description ...");
    cpl_test_eq_string(cpl_plugin_get_author(parsed_plugin), "some author");
    cpl_test_eq_string(cpl_plugin_get_email(parsed_plugin), "author@eso.org");
    cpl_test_eq_string(cpl_plugin_get_copyright(parsed_plugin),
                       "copyright ...");

    cpl_test_eq(cpl_parameterlist_get_size(parsed_recipe2->base.parameters), 2);
    parsed_param0 = cpl_parameterlist_find(parsed_recipe2->base.parameters,
                                           "__class__");
    cpl_test_nonnull(parsed_param0);
    cpl_test_eq(cpl_parameter_get_type(parsed_param0), CPL_TYPE_STRING);
    cpl_test_eq_string(cpl_parameter_get_string(parsed_param0), "unknown");
    cpl_test_eq_string(cpl_parameter_get_default_string(parsed_param0),
                       "unknown");
    parsed_param1 = cpl_parameterlist_find(parsed_recipe2->base.parameters,
                                           "test.par1");
    cpl_test_nonnull(parsed_param1);
    test_equal_param(param1, parsed_param1);

    cpl_test_eq(cpl_frameset_get_size(parsed_recipe2->base.frames), 1);
    parsed_frame1 = cpl_frameset_get_position(parsed_recipe2->base.frames, 0);
    cpl_test_nonnull(parsed_frame1);
    test_equal_frame(frame1, parsed_frame1);

    tags = cpl_recipeconfig_get_tags(parsed_recipe2->config);
    cpl_test_nonnull(tags);
    cpl_test_eq_string(tags[0], "RAW");
    cpl_test_eq_string(tags[1], "RAW2");
    cpl_test_null(tags[2]);
    for (charptr = tags; *charptr != NULL; ++charptr) cpl_free(*charptr);
    cpl_free(tags);

    inputs = cpl_recipeconfig_get_inputs(parsed_recipe2->config, "RAW");
    cpl_test_nonnull(inputs);
    cpl_test_eq_string(inputs[0], "CALIB");
    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW", "CALIB"),
                2);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW", "CALIB"),
                3);
    cpl_test_null(inputs[1]);
    for (charptr = inputs; *charptr != NULL; ++charptr) cpl_free(*charptr);
    cpl_free(inputs);

    outputs = cpl_recipeconfig_get_outputs(parsed_recipe2->config, "RAW");
    cpl_test_nonnull(outputs);
    cpl_test_eq_string(inputs[0], "PROD");
    cpl_test_null(outputs[1]);
    for (charptr = outputs; *charptr != NULL; ++charptr) cpl_free(*charptr);
    cpl_free(outputs);

    inputs = cpl_recipeconfig_get_inputs(parsed_recipe2->config, "RAW2");
    cpl_test_nonnull(inputs);
    cpl_test_eq_string(inputs[0], "CALIB");
    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW2", "CALIB"),
                2);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW2", "CALIB"),
                2);
    cpl_test_eq_string(inputs[1], "DARK");
    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW2", "DARK"),
                1);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW2", "DARK"),
                1);
    cpl_test_null(inputs[2]);
    for (charptr = inputs; *charptr != NULL; ++charptr) cpl_free(*charptr);
    cpl_free(inputs);

    outputs = cpl_recipeconfig_get_outputs(parsed_recipe2->config, "RAW2");
    cpl_test_nonnull(outputs);
    cpl_test_eq_string(inputs[0], "PROD2");
    cpl_test_eq_string(inputs[1], "PROD3");
    cpl_test_null(outputs[2]);
    for (charptr = outputs; *charptr != NULL; ++charptr) cpl_free(*charptr);
    cpl_free(outputs);

    cpl_free(json);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(plugin);
    parsed_plugin->deinitialize(parsed_plugin);
    cpl_plugin_delete(parsed_plugin);
    er_json_node_delete(parsetree);

    /* Test parsing of the recipeconfig object if negative values are used for
       the minimum and maximum.
     */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));
    recipe2->base.parameters = cpl_parameterlist_new();
    cpl_test_assert(recipe2->base.parameters != NULL);
    recipe2->base.frames = cpl_frameset_new();
    cpl_test_assert(recipe2->base.frames != NULL);

    recipe2->config = cpl_recipeconfig_new();
    cpl_test_assert(recipe2->config != NULL);
    cpl_recipeconfig_set_tag(recipe2->config, "RAW1", -1, 5);
    cpl_recipeconfig_set_input(recipe2->config, "RAW1", "CALIB1", -1, 3);
    cpl_recipeconfig_set_input(recipe2->config, "RAW1", "CALIB2", 2, -1);
    cpl_recipeconfig_set_input(recipe2->config, "RAW1", "CALIB3", -1, -1);
    cpl_recipeconfig_set_tag(recipe2->config, "RAW2", 4, -1);
    cpl_recipeconfig_set_tag(recipe2->config, "RAW3", -1, -1);
    cpl_test_error(CPL_ERROR_NONE);

    plugin = &recipe2->base.interface;
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE_V2, "test_recipe2",
                                    "simple test recipe v2", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(json);

    parsetree = er_json_parse(json);
    cpl_test_nonnull(parsetree);
    cpl_test_error(CPL_ERROR_NONE);
    parsed_plugin = er_json_to_plugin(parsetree, json);
    parsed_recipe2 = (cpl_recipe2 *) parsed_plugin;
    cpl_test_nonnull(parsed_plugin);
    cpl_test_error(CPL_ERROR_NONE);

    tags = cpl_recipeconfig_get_tags(parsed_recipe2->config);
    cpl_test_nonnull(tags);
    cpl_test_eq_string(tags[0], "RAW1");
    cpl_test_eq_string(tags[1], "RAW2");
    cpl_test_eq_string(tags[2], "RAW3");
    cpl_test_null(tags[3]);
    for (charptr = tags; *charptr != NULL; ++charptr) cpl_free(*charptr);
    cpl_free(tags);

    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW1", "RAW1"),
                -1);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW1", "RAW1"),
                5);
    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW2", "RAW2"),
                4);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW2", "RAW2"),
                -1);
    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW3", "RAW3"),
                -1);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW3", "RAW3"),
                -1);

    inputs = cpl_recipeconfig_get_inputs(parsed_recipe2->config, "RAW1");
    cpl_test_nonnull(inputs);
    cpl_test_eq_string(inputs[0], "CALIB1");
    cpl_test_eq_string(inputs[1], "CALIB2");
    cpl_test_eq_string(inputs[2], "CALIB3");
    cpl_test_null(inputs[3]);
    for (charptr = inputs; *charptr != NULL; ++charptr) cpl_free(*charptr);
    cpl_free(inputs);

    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW1", "CALIB1"),
                -1);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW1", "CALIB1"),
                3);
    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW1", "CALIB2"),
                2);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW1", "CALIB2"),
                -1);
    cpl_test_eq(cpl_recipeconfig_get_min_count(parsed_recipe2->config,
                                               "RAW1", "CALIB3"),
                -1);
    cpl_test_eq(cpl_recipeconfig_get_max_count(parsed_recipe2->config,
                                               "RAW1", "CALIB3"),
                -1);

    cpl_free(json);
    cpl_parameterlist_delete(recipe2->base.parameters);
    cpl_frameset_delete(recipe2->base.frames);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(plugin);
    parsed_plugin->deinitialize(parsed_plugin);
    cpl_plugin_delete(parsed_plugin);
    er_json_node_delete(parsetree);

    /* Test parsing of the recipe configuration structure when there is only
       one tag and it has no inputs or outputs. */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));

    recipe2->config = cpl_recipeconfig_new();
    cpl_test_assert(recipe2->config != NULL);
    cpl_recipeconfig_set_tag(recipe2->config, "RAW", 1, 1);
    cpl_test_error(CPL_ERROR_NONE);

    plugin = &recipe2->base.interface;
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE_V2, "test_recipe2",
                                    "simple test recipe v2", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(json);

    parsetree = er_json_parse(json);
    cpl_test_nonnull(parsetree);
    cpl_test_error(CPL_ERROR_NONE);
    parsed_plugin = er_json_to_plugin(parsetree, json);
    parsed_recipe2 = (cpl_recipe2 *) parsed_plugin;
    cpl_test_nonnull(parsed_plugin);
    cpl_test_error(CPL_ERROR_NONE);

    tags = cpl_recipeconfig_get_tags(parsed_recipe2->config);
    cpl_test_nonnull(tags);
    cpl_test_eq_string(tags[0], "RAW");
    cpl_test_null(tags[1]);
    for (charptr = tags; *charptr != NULL; ++charptr) cpl_free(*charptr);
    cpl_free(tags);

    inputs = cpl_recipeconfig_get_inputs(parsed_recipe2->config, "RAW");
    cpl_test_nonnull(inputs);
    cpl_test_null(inputs[0]);
    cpl_free(inputs);

    outputs = cpl_recipeconfig_get_outputs(parsed_recipe2->config, "RAW");
    cpl_test_nonnull(outputs);
    cpl_test_null(outputs[0]);
    cpl_free(outputs);

    cpl_free(json);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(plugin);
    parsed_plugin->deinitialize(parsed_plugin);
    cpl_plugin_delete(parsed_plugin);
    er_json_node_delete(parsetree);

    /* Test parsing of the recipe configuration structure when it is empty.
       */
    recipe2 = (cpl_recipe2 *) cpl_calloc(1, sizeof(cpl_recipe2));

    recipe2->config = cpl_recipeconfig_new();
    cpl_test_assert(recipe2->config != NULL);

    plugin = &recipe2->base.interface;
    cpl_test_assert(cpl_plugin_init(plugin, CPL_PLUGIN_API, 1,
                                    CPL_PLUGIN_TYPE_RECIPE_V2, "test_recipe2",
                                    "simple test recipe v2", "description ...",
                                    "some author", "author@eso.org",
                                    "copyright ...", NULL, NULL, NULL)
                    == CPL_ERROR_NONE);
    json = er_plugin_to_json(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_nonnull(json);

    parsetree = er_json_parse(json);
    cpl_test_nonnull(parsetree);
    cpl_test_error(CPL_ERROR_NONE);
    parsed_plugin = er_json_to_plugin(parsetree, json);
    parsed_recipe2 = (cpl_recipe2 *) parsed_plugin;
    cpl_test_nonnull(parsed_plugin);
    cpl_test_error(CPL_ERROR_NONE);

    tags = cpl_recipeconfig_get_tags(parsed_recipe2->config);
    cpl_test_nonnull(tags);
    cpl_test_null(tags[0]);
    cpl_free(tags);

    cpl_free(json);
    cpl_recipeconfig_delete(recipe2->config);
    cpl_plugin_delete(plugin);
    parsed_plugin->deinitialize(parsed_plugin);
    cpl_plugin_delete(parsed_plugin);
    er_json_node_delete(parsetree);
}

static void esorex_test_plugin_parsing(void)
{
    const char * json = NULL;
    cpl_plugin * plugin = NULL;
    cpl_recipe * recipe = NULL;
    cpl_parameter * ref_param1 = NULL;
    cpl_parameter * parsed_param0 = NULL;
    cpl_parameter * parsed_param1 = NULL;
    cpl_frame * ref_frame = NULL;
    cpl_frame * parsed_frame = NULL;
    er_json_node * parsetree = NULL;

    /* Test er_json_to_plugin produces a correct object for the corresponding
       JSON text. */
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
        "    \"parameters\": [\n"
        "        {\n"
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
        "        }\n"
        "    ],\n"
        "    \"frames\": [\n"
        "        {\n"
        "            \"filename\": \"test.fits\",\n"
        "            \"tag\": \"RAW\",\n"
        "            \"type\": 2,\n"
        "            \"group\": 1,\n"
        "            \"level\": 3\n"
        "        }\n"
        "    ]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    plugin = er_json_to_plugin(parsetree, json);
    cpl_test_nonnull(plugin);
    cpl_test_error(CPL_ERROR_NONE);

    /* Prepare reference objects to compare to. */
    ref_param1 = cpl_parameter_new_value("test.par1", CPL_TYPE_INT,
                                         "int param", "test", 5);
    cpl_test_assert(ref_param1 != NULL);
    cpl_parameter_set_int(ref_param1, 12);
    cpl_parameter_set_tag(ref_param1, "myparam");
    cpl_parameter_enable(ref_param1, CPL_PARAMETER_MODE_CLI);
    cpl_parameter_set_alias(ref_param1, CPL_PARAMETER_MODE_CLI, "par1_cli");
    cpl_parameter_disable(ref_param1, CPL_PARAMETER_MODE_ENV);
    cpl_parameter_set_alias(ref_param1, CPL_PARAMETER_MODE_ENV, "par1_env");
    cpl_parameter_enable(ref_param1, CPL_PARAMETER_MODE_CFG);
    cpl_parameter_set_alias(ref_param1, CPL_PARAMETER_MODE_CFG, "par1_cfg");
    cpl_test_error(CPL_ERROR_NONE);

    ref_frame = cpl_frame_new();
    cpl_test_assert(ref_frame != NULL);
    cpl_frame_set_filename(ref_frame, "test.fits");
    cpl_frame_set_tag(ref_frame, "RAW");
    cpl_frame_set_type(ref_frame, CPL_FRAME_TYPE_IMAGE);
    cpl_frame_set_group(ref_frame, CPL_FRAME_GROUP_RAW);
    cpl_frame_set_level(ref_frame, CPL_FRAME_LEVEL_FINAL);

    /* Check the produced plugin structure. */
    cpl_test_eq(cpl_plugin_get_version(plugin), 123);
    cpl_test_eq(cpl_plugin_get_type(plugin), CPL_PLUGIN_TYPE_RECIPE);
    cpl_test_eq_string(cpl_plugin_get_name(plugin), "testrecipe");
    cpl_test_eq_string(cpl_plugin_get_synopsis(plugin), "test plugin");
    cpl_test_eq_string(cpl_plugin_get_description(plugin), "description ...");
    cpl_test_eq_string(cpl_plugin_get_author(plugin), "Some Author");
    cpl_test_eq_string(cpl_plugin_get_email(plugin), "someone@eso.org");
    cpl_test_eq_string(cpl_plugin_get_copyright(plugin), "copyright ...");

    recipe = (cpl_recipe *) plugin;
    cpl_test_eq(cpl_parameterlist_get_size(recipe->parameters), 2);
    parsed_param0 = cpl_parameterlist_find(recipe->parameters, "__class__");
    cpl_test_nonnull(parsed_param0);
    cpl_test_eq(cpl_parameter_get_type(parsed_param0), CPL_TYPE_STRING);
    cpl_test_eq_string(cpl_parameter_get_string(parsed_param0), "testclass");
    cpl_test_eq_string(cpl_parameter_get_default_string(parsed_param0),
                       "testclass");
    parsed_param1 = cpl_parameterlist_find(recipe->parameters, "test.par1");
    cpl_test_nonnull(parsed_param1);
    test_equal_param(ref_param1, parsed_param1);

    cpl_test_eq(cpl_frameset_get_size(recipe->frames), 1);
    parsed_frame = cpl_frameset_get_position(recipe->frames, 0);
    cpl_test_nonnull(parsed_frame);
    test_equal_frame(ref_frame, parsed_frame);

    cpl_parameter_delete(ref_param1);
    cpl_frame_delete(ref_frame);

    plugin->deinitialize(plugin);
    cpl_plugin_delete(plugin);

    /* Test for correct error handling if NULL pointers are passed to the
       er_json_to_plugin function. */
    cpl_test_null(er_json_to_plugin(NULL, json));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(er_json_to_plugin(parsetree, NULL));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    er_json_node_delete(parsetree);

    /* Check that defaults are used if certain keywords are missing. */
    json =
        "{\n"
        "    \"class\": \"testclass\",\n"
        "    \"parameters\": [],\n"
        "    \"frames\": []\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    plugin = er_json_to_plugin(parsetree, json);
    cpl_test_nonnull(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq(cpl_plugin_get_version(plugin), 0);
    cpl_test_eq(cpl_plugin_get_type(plugin), CPL_PLUGIN_TYPE_RECIPE);
    cpl_test_eq_string(cpl_plugin_get_name(plugin), "testclass");
    cpl_test_eq_string(cpl_plugin_get_synopsis(plugin), "unknown");
    cpl_test_eq_string(cpl_plugin_get_description(plugin), "unknown");
    cpl_test_eq_string(cpl_plugin_get_author(plugin), "unknown");
    cpl_test_eq_string(cpl_plugin_get_email(plugin), "unknown");
    cpl_test_eq_string(cpl_plugin_get_copyright(plugin), "unknown");
    plugin->deinitialize(plugin);
    cpl_plugin_delete(plugin);
    er_json_node_delete(parsetree);

    /* Test error handling if the JSON does not contain an object/dictionary. */
    json = "[\"dummy\", 123]";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if the JSON object does not contain valid plugin
       structure. */
    json = "{\"dummy\": 123}";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if the 'class' keyword is missing. */
    json =
        "{\n"
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
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if the 'parameters' keyword is missing. */
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
        "    \"frames\": []\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if 'parameters' is not an array. */
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
        "    \"parameters\": {\"dummy\": 123},\n"
        "    \"frames\": []\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if the 'parameters' array contains invalid
       structures. */
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
        "    \"parameters\": [{\"dummy\": 123}],\n"
        "    \"frames\": []\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Make sure the parsing works if the 'frames' keyword is missing. */
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
        "    \"parameters\": []\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    plugin = er_json_to_plugin(parsetree, json);
    cpl_test_nonnull(plugin);
    cpl_test_error(CPL_ERROR_NONE);
    plugin->deinitialize(plugin);
    cpl_plugin_delete(plugin);
    er_json_node_delete(parsetree);

    /* Check error handling if 'frames' is not an array. */
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
        "    \"frames\": {\"dummy\": 123}\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if the 'frames' array contains invalid
       structures. */
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
        "    \"frames\": [{\"dummy\": 123}]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if a 'recipeconfig' is not an array. */
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
        "    \"frames\": [],\n"
        "    \"recipeconfig\": {\"dummy\": 123}\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Check error handling if a 'recipeconfig' array contains invalid
       structures. */
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
        "    \"frames\": [],\n"
        "    \"recipeconfig\": [{\"dummy\": 123}]\n"
        "}\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_plugin(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);
}

static void esorex_test_pluginlist_parsing(void)
{
    const char * json = NULL;
    cpl_pluginlist * list = NULL;
    cpl_plugin * plugin = NULL;
    cpl_recipe * recipe = NULL;
    cpl_parameter * param = NULL;
    er_json_node * parsetree = NULL;

    /* Test er_json_to_pluginlist produces a correct plugin list for the
       corresponding JSON text. */
    json =
        "[\n"
        "  {\n"
        "    \"class\": \"testclass1\",\n"
        "    \"name\": \"testrecipe1\",\n"
        "    \"version\": 123,\n"
        "    \"synopsis\": \"test plugin\",\n"
        "    \"description\": \"description ...\",\n"
        "    \"author\": \"Some Author\",\n"
        "    \"email\": \"someone@eso.org\",\n"
        "    \"copyright\": \"copyright ...\",\n"
        "    \"parameters\": [],\n"
        "    \"frames\": []\n"
        "  },\n"
        "  {\n"
        "    \"class\": \"testclass2\",\n"
        "    \"name\": \"testrecipe2\",\n"
        "    \"version\": 345,\n"
        "    \"synopsis\": \"test plugin 2\",\n"
        "    \"description\": \"description 2 ...\",\n"
        "    \"author\": \"Some Author 2\",\n"
        "    \"email\": \"someone2@eso.org\",\n"
        "    \"copyright\": \"copyright 2 ...\",\n"
        "    \"parameters\": [],\n"
        "    \"frames\": []\n"
        "  }\n"
        "]\n";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    list = er_json_to_pluginlist(parsetree, json);
    cpl_test_nonnull(list);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_test_eq(cpl_pluginlist_get_size(list), 2);

    plugin = cpl_pluginlist_get_first(list);
    cpl_test_nonnull(plugin);
    cpl_test_eq(cpl_plugin_get_version(plugin), 123);
    cpl_test_eq(cpl_plugin_get_type(plugin), CPL_PLUGIN_TYPE_RECIPE);
    cpl_test_eq_string(cpl_plugin_get_name(plugin), "testrecipe1");
    cpl_test_eq_string(cpl_plugin_get_synopsis(plugin), "test plugin");
    cpl_test_eq_string(cpl_plugin_get_description(plugin), "description ...");
    cpl_test_eq_string(cpl_plugin_get_author(plugin), "Some Author");
    cpl_test_eq_string(cpl_plugin_get_email(plugin), "someone@eso.org");
    cpl_test_eq_string(cpl_plugin_get_copyright(plugin), "copyright ...");

    recipe = (cpl_recipe *) plugin;
    cpl_test_eq(cpl_parameterlist_get_size(recipe->parameters), 1);
    param = cpl_parameterlist_find(recipe->parameters, "__class__");
    cpl_test_nonnull(param);
    cpl_test_eq(cpl_parameter_get_type(param), CPL_TYPE_STRING);
    cpl_test_eq_string(cpl_parameter_get_string(param), "testclass1");
    cpl_test_eq_string(cpl_parameter_get_default_string(param), "testclass1");

    plugin = cpl_pluginlist_get_last(list);
    cpl_test_nonnull(plugin);
    cpl_test_eq(cpl_plugin_get_version(plugin), 345);
    cpl_test_eq(cpl_plugin_get_type(plugin), CPL_PLUGIN_TYPE_RECIPE);
    cpl_test_eq_string(cpl_plugin_get_name(plugin), "testrecipe2");
    cpl_test_eq_string(cpl_plugin_get_synopsis(plugin), "test plugin 2");
    cpl_test_eq_string(cpl_plugin_get_description(plugin), "description 2 ...");
    cpl_test_eq_string(cpl_plugin_get_author(plugin), "Some Author 2");
    cpl_test_eq_string(cpl_plugin_get_email(plugin), "someone2@eso.org");
    cpl_test_eq_string(cpl_plugin_get_copyright(plugin), "copyright 2 ...");

    recipe = (cpl_recipe *) plugin;
    cpl_test_eq(cpl_parameterlist_get_size(recipe->parameters), 1);
    param = cpl_parameterlist_find(recipe->parameters, "__class__");
    cpl_test_nonnull(param);
    cpl_test_eq(cpl_parameter_get_type(param), CPL_TYPE_STRING);
    cpl_test_eq_string(cpl_parameter_get_string(param), "testclass2");
    cpl_test_eq_string(cpl_parameter_get_default_string(param), "testclass2");

    /* Make sure the delete method does not produce an error. */
    er_json_pluginlist_delete(list);
    cpl_test_error(CPL_ERROR_NONE);

    /* Make sure er_json_pluginlist_delete does not produce an error, even if
       it is given a NULL pointer. */
    er_json_pluginlist_delete(NULL);
    cpl_test_error(CPL_ERROR_NONE);

    /* Test for correct error handling if NULL pointers are passed to the
       er_json_to_pluginlist function. */
    cpl_test_null(er_json_to_pluginlist(NULL, json));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_null(er_json_to_pluginlist(parsetree, NULL));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if the JSON does not contain an array. */
    json = "{\"dummy\": 123}";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_pluginlist(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);

    /* Test error handling if the JSON array does not contain valid plugin
       structures. */
    json = "[\"dummy\", 123]";
    parsetree = er_json_parse(json);
    cpl_test_assert(parsetree != NULL);
    cpl_test_null(er_json_to_pluginlist(parsetree, json));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    er_json_node_delete(parsetree);
}

#endif /* ENABLE_PYTHON_RECIPES */
