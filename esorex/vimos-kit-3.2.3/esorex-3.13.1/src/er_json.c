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
#include <config.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <er_json.h>
#include <cxmap.h>
#include <cxlist.h>

#ifdef HAVE_LIBAVCALL
#include <avcall.h>
#else
  #ifdef HAVE_LIBFFI
  #include <ffi.h>
  #endif
#endif


/**
 * @defgroup esorex_json_utils JSON Utilities
 *
 * This module allows to write some CPL objects to JSON format.
 * There is a recursive descent parser provided that also allows to parse JSON
 * text into a parse tree.
 * The parser should be compatible with UTF-8 and has an expanded number range,
 * specifically NaN and infinity are accepted.
 * Additional utility functions provide conversion of the parse tree back into
 * CPL objects.
 *
 * @par Synopsis:
 * @code
 *   #include <er_json.h>
 * @endcode
 */

/**@{*/

/**
 * @brief
 *   Write a JSON file with the contents of a frameset.
 *
 * @param frameset  The frameset to export
 * @param json_file Name of the JSON file to write 
 *
 * @return @c CPL_ERROR_NONE on success.
 *
 * This function will export the content of a frameset to a JSON format
 * which can be machine readable. The format looks like:
 * <pre>
 * [
 *  {
 *    "name": "ifu_trace.fits",
 *    "category": "IFU_TRACE"
 *  },
 *  {
 *    "name": "ifu_transmission.fits",
 *    "category": "IFU_TRANSMISSION"
 *  }
 * ]
 * </pre>
 */
cpl_error_code er_frameset_to_json(const cpl_frameset * frameset, 
                                   const char* json_file)
{
    FILE             * fjson;
    cpl_size           nframes;
    cpl_size           iframe;

    /* Sanity checks  */ 
    if(!frameset || !json_file)
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Null input");

    /* Open the output file */
    fjson = fopen(json_file, "w");
    if(!fjson)
    {
        errno = 0;
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                     "Cannot open JSON file for writing");
    }

    /* Open a JSON list */
    fprintf(fjson, "[\n");

    /* Loop on frames */
    nframes = cpl_frameset_get_size(frameset);
    iframe = 0;

    for(iframe = 0 ; iframe < nframes; iframe++)
    {
        const cpl_frame  * frame;
        char             * filename;
        char             * category;

        /* Get this frame */
        frame = cpl_frameset_get_position_const(frameset, iframe);
        filename = er_json_escape_string(cpl_frame_get_filename(frame));
        category = er_json_escape_string(cpl_frame_get_tag(frame));

        /* Open a JSON item */
        fprintf(fjson, "  {\n");

        /* Write this frame fields */
        if(strlen(filename) > 0)
            fprintf(fjson, "    \"name\": \"%s\",\n",filename);
        if(strlen(category) > 0)
            fprintf(fjson, "    \"category\": \"%s\"\n",category);

        /* Close JSON item */
        if(iframe == cpl_frameset_get_size(frameset) - 1)
            fprintf(fjson, "  }\n");
        else
            fprintf(fjson, "  },\n");
        cpl_free(filename);
        cpl_free(category);
    }
    /* Close JSON list */
    fprintf(fjson, "]\n");
    fclose(fjson);

    return CPL_ERROR_NONE;
}

/**
 * @brief
 *   Write a text file with the contents of a frameset.
 *
 * @param frameset  The frameset to export
 * @param text_file Name of the JSON file to write 
 *
 * @return @c CPL_ERROR_NONE on success.
 *
 * This function will export the contents of the frameset to a text human
 * readable file, like the one used for the input sof. 
 * The output will look like this:
 * <pre>
 * ifu_trace.fits     IFU_TRACE
 * ifu_transmission.fits     IFU_TRANSMISSION
 * </pre>
 */
cpl_error_code er_frameset_to_text(const cpl_frameset * frameset, 
                                   const char* text_file)
{
    FILE             * ftext;
    cpl_size           nframes;
    cpl_size           iframe;

    /* Sanity checks  */ 
    if(!frameset || !text_file)
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Null input");

    /* Open the output file */
    ftext = fopen(text_file, "w");
    if(!ftext)
    {
        errno = 0;
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                     "Cannot open file %s for writing", text_file);
    }

    /* Loop on frames */
    nframes = cpl_frameset_get_size(frameset);
    iframe = 0;

    for(iframe = 0 ; iframe < nframes; iframe++)
    {
        const cpl_frame  * frame;
        const char       * filename;
        const char       * category;

        /* Get this frame */
        frame = cpl_frameset_get_position_const(frameset, iframe);
        filename = cpl_frame_get_filename(frame);
        category = cpl_frame_get_tag(frame);

        /* Writing a frame*/
        fprintf(ftext, "%s     %s", filename, category);
        if(iframe != cpl_frameset_get_size(frameset) - 1)
            fprintf(ftext, "\n");

        iframe++;
    }
    fclose(ftext);

    return CPL_ERROR_NONE;
}

/**
 * @brief
 *   Write a JSON file with the contents of a parameterlist.
 *
 * @param plist       The parameter list  to export
 * @param recipe_name The recipe this parameter list applies to
 * @param json_file   Name of the JSON file to write 
 *
 * @return @c CPL_ERROR_NONE on success.
 *
 * This function exports the content of a parameter list to a machine readable
 * JSON format. The output looks like:
 * <pre>
 * [
 *  {
 *    "name": "vimos.Parameters.stacking.singleframes",
 *    "value": true,
 *    "display_name": "AllowSingleFrames",
 *    "description": "Frame combination method is ignored.",
 *    "recipe": "vmbias"
 *  },
 *  {
 *    "name": "vimos.Parameters.stacking.method",
 *    "value": "Median",
 *    "display_name": "StackMethod",
 *    "description": "Frame combination method",
 *    "recipe": "vmbias"
 *  }
 * ]
 * </pre>
 */
cpl_error_code er_recipe_parameterlist_to_json(const cpl_parameterlist * plist,
                                               const char * recipe_name,
                                               const char * json_file)
{
    FILE                * fjson;
    const cpl_parameter * param;

    /* Sanity checks  */ 
    if(!plist || !recipe_name || !json_file)
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Null input");

    /* Open the output file */
    fjson = fopen(json_file, "w");
    if(!fjson)
    {
        errno = 0;
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_NOT_CREATED,
                                     "Cannot open JSON file for writing");
    }

    /* Open a JSON list */
    fprintf(fjson, "[\n");

    /* Loop on parameters */
    param = cpl_parameterlist_get_first_const(plist);

    while (param != NULL)
    {
        if (cpl_parameter_is_enabled(param, CPL_PARAMETER_MODE_CFG))
        {
            char              * description;
            char              * name;
            char              * display_name;
            cpl_type            valtype;
            cpl_parameter_class partype;
            int                 ivalue;
            int                 bvalue;
            double              dvalue;
            char              * svalue;
            int                 enum_size;
            int                 ienum;

            description = er_json_escape_string(cpl_parameter_get_help(param));
            name = er_json_escape_string(cpl_parameter_get_name(param));
            display_name = er_json_escape_string(
                    cpl_parameter_get_alias(param, CPL_PARAMETER_MODE_CLI));
            if(display_name == NULL)
            {
                cpl_error_reset();
                display_name = er_json_escape_string(cpl_parameter_get_name(param));
            }
            valtype = cpl_parameter_get_type(param);
            partype = cpl_parameter_get_class(param);

            /* Open a JSON item */
            fprintf(fjson, "  {\n");

            /* Write this parameter fields */
            if(strlen(name) > 0)
                fprintf(fjson, "    \"name\": \"%s\",\n",name);
            /* Write value and value type */
            switch (valtype)
            {
                case CPL_TYPE_BOOL:
                    bvalue = cpl_parameter_get_bool(param);
                    if (bvalue == 0)
                        fprintf(fjson, "    \"value\": false,\n");
                    else
                        fprintf(fjson, "    \"value\": true,\n");
                    fprintf(fjson, "    \"valtype\": \"bool\",\n");
                    break;

                case CPL_TYPE_INT:
                    ivalue = cpl_parameter_get_int(param);
                    fprintf(fjson, "    \"value\": %d,\n",ivalue);
                    fprintf(fjson, "    \"valtype\": \"int\",\n");
                    break;

                case CPL_TYPE_DOUBLE:
                    dvalue = cpl_parameter_get_double(param);
                    fprintf(fjson, "    \"value\": %.17g,\n",dvalue); //IEEE 754 double
                    fprintf(fjson, "    \"valtype\": \"double\",\n");
                    break;

                case CPL_TYPE_STRING:
                    svalue = er_json_escape_string(cpl_parameter_get_string(param));
                    fprintf(fjson, "    \"value\": \"%s\",\n",svalue);
                    fprintf(fjson, "    \"valtype\": \"string\",\n");
                    cpl_free(svalue);
                    break;

                default:
                    return cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                                 "Type not supported for cpl_propertylist");
                    break;
            }

            /* Write parameter type (called class inside CPL) */
            switch (partype)
            {
                case CPL_PARAMETER_CLASS_VALUE:
                    fprintf(fjson, "    \"partype\": \"value\",\n");
                    break;

                    /* If a range, write the minimum and maximum.
               Only integers and doubles are supported */
                case CPL_PARAMETER_CLASS_RANGE:
                    fprintf(fjson, "    \"partype\": \"range\",\n");
                    switch (valtype)
                    {
                        case CPL_TYPE_INT:
                            ivalue = cpl_parameter_get_range_min_int(param);
                            fprintf(fjson, "    \"valmin\": %d,\n",ivalue);
                            ivalue = cpl_parameter_get_range_max_int(param);
                            fprintf(fjson, "    \"valmax\": %d,\n",ivalue);
                            break;

                        case CPL_TYPE_DOUBLE:
                            dvalue = cpl_parameter_get_range_min_double(param);
                            fprintf(fjson, "    \"valmin\": %.17g,\n",dvalue);
                            dvalue = cpl_parameter_get_range_max_double(param);
                            fprintf(fjson, "    \"valmax\": %.17g,\n",dvalue);
                            break;
                        default:
                            return cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                                         "Range parameters not supported for this type");
                            break;
                    }
                    break;
                    /* If an enum, write the list of alternatives as a JSON list.
               Only integers, doubles and strings are supported */
                        case CPL_PARAMETER_CLASS_ENUM:
                            fprintf(fjson, "    \"partype\": \"enum\",\n");
                            enum_size = cpl_parameter_get_enum_size(param);
                            fprintf(fjson, "    \"valenum\": [");
                            for(ienum = 0; ienum < enum_size; ienum++)
                            {
                                switch (valtype)
                                {
                                    case CPL_TYPE_INT:
                                        ivalue = cpl_parameter_get_enum_int(param, ienum);
                                        fprintf(fjson, " %d ",ivalue);
                                        break;

                                    case CPL_TYPE_DOUBLE:
                                        dvalue = cpl_parameter_get_enum_double(param, ienum);
                                        fprintf(fjson, " %.17g ",dvalue);
                                        break;

                                    case CPL_TYPE_STRING:
                                        svalue = er_json_escape_string(
                                                cpl_parameter_get_enum_string(param, ienum));
                                        fprintf(fjson, " \"%s\" ",svalue);
                                        cpl_free(svalue);
                                        break;
                                    default:
                                        return cpl_error_set_message(cpl_func,
                                                                     CPL_ERROR_TYPE_MISMATCH,
                                                                     "Enum parameters not supported for this type");
                                        break;
                                }
                                if(ienum != enum_size - 1)
                                    fprintf(fjson, ",");
                            }
                            fprintf(fjson, "],\n");
                            break;
                        default:
                            return cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                                         "Parameter class not supported for cpl_propertylist");
                            break;
            }
            if(strlen(display_name) > 0)
                fprintf(fjson, "    \"display_name\": \"%s\",\n",display_name);
            if(strlen(description) > 0)
            {
                if(strlen(recipe_name) > 0) //Description is not the last field
                    fprintf(fjson, "    \"description\": \"%s\",\n",description);
                else                        //Description is the last field
                    fprintf(fjson, "    \"description\": \"%s\"\n",description);
            }
            if(strlen(recipe_name) > 0)
                fprintf(fjson, "    \"recipe\": \"%s\"\n",recipe_name);

            /* Close JSON item */
            if(param == cpl_parameterlist_get_last_const(plist))
                fprintf(fjson, "  }\n");
            else
                fprintf(fjson, "  },\n");

            cpl_free(description);
            cpl_free(name);
            cpl_free(display_name);
        }

        param = cpl_parameterlist_get_next_const(plist);

    }
    /* Close JSON list */
    fprintf(fjson, "]\n");
    fclose(fjson);

    return CPL_ERROR_NONE;
}

char * er_json_escape_string(const char * str)
{
    char * escaped_str;
    size_t i, j;
    if(str == NULL)
        return NULL;

    escaped_str = cpl_malloc(2*strlen(str) * sizeof(char)+1);
    for(i = 0, j = 0; i < strlen(str); i++)
    {
        if(str[i] == 0x22 || str[i] == 0x5C || str[i] == 0x2F )
        {
            escaped_str[j++] = '\\';
            escaped_str[j++] = str[i];
        }
        else if(str[i] == 0x08)
        {
            escaped_str[j++] = '\\';
            escaped_str[j++] = 'b';
        }
        else if(str[i] == 0x09)
        {
            escaped_str[j++] = '\\';
            escaped_str[j++] = 't';
        }
        else if(str[i] == 0x0A)
        {
            escaped_str[j++] = '\\';
            escaped_str[j++] = 'n';
        }
        else if(str[i] == 0x0C)
        {
            escaped_str[j++] = '\\';
            escaped_str[j++] = 'f';
        }
        else if(str[i] == 0x0D)
        {
            escaped_str[j++] = '\\';
            escaped_str[j++] = 'r';
        }
        else
            escaped_str[j++] = str[i];
    }
    escaped_str[j++]='\0';    
    escaped_str = cpl_realloc(escaped_str, j);

    return escaped_str;
}

/*
 * Creates and returns a properly escaped JSON string surrounded by quotes.
 * The result must be deleted with cpl_free().
 */
static char * prepare_escaped_string(const char * string)
{
    if (string != NULL)
    {
        char * tmp = er_json_escape_string(string);
        char * result = cpl_sprintf("\"%s\"", tmp);
        cpl_free(tmp);
        return result;
    }
    else
    {
        return cpl_strdup("null");
    }
}

/*
 * Serialises a CPL parameter object to properly formatted JSON text and
 * returned as a null terminated string.
 * If an error occurred then NULL is returned. Otherwise the result must be
 * cleaned up by the caller with cpl_free().
 */
static char * parameter_to_json(const cpl_parameter * param)
{
    assert(param != NULL);

    const char * present = cpl_parameter_get_default_flag(param) ? "true"
                                                                 : "false";
    int param_id = cpl_parameter_get_id(param);
    const char * tag = cpl_parameter_get_tag(param);
    const char * cli_enabled
        = cpl_parameter_is_enabled(param, CPL_PARAMETER_MODE_CLI) ? "true"
                                                                  : "false";
    const char * env_enabled
        = cpl_parameter_is_enabled(param, CPL_PARAMETER_MODE_ENV) ? "true"
                                                                  : "false";
    const char * cfg_enabled
        = cpl_parameter_is_enabled(param, CPL_PARAMETER_MODE_CFG) ? "true"
                                                                  : "false";
    const char * cli_alias = cpl_parameter_get_alias(param,
                                                     CPL_PARAMETER_MODE_CLI);
    const char * env_alias = cpl_parameter_get_alias(param,
                                                     CPL_PARAMETER_MODE_ENV);
    const char * cfg_alias = cpl_parameter_get_alias(param,
                                                     CPL_PARAMETER_MODE_CFG);
    if (cpl_error_get_code() != CPL_ERROR_NONE) return NULL;

    char * tmp = NULL;
    char * value = NULL;
    char * default_value = NULL;
    cpl_type value_type = cpl_parameter_get_type(param);
    switch (value_type)
    {
        case CPL_TYPE_BOOL:
            value = (cpl_parameter_get_bool(param) == TRUE)
                    ? cpl_strdup("true") : cpl_strdup("false");
            default_value = (cpl_parameter_get_default_bool(param) == TRUE)
                            ? cpl_strdup("true") : cpl_strdup("false");
            break;

        case CPL_TYPE_INT:
            value = cpl_sprintf("%d", cpl_parameter_get_int(param));
            default_value =
                cpl_sprintf("%d", cpl_parameter_get_default_int(param));
            break;

        case CPL_TYPE_DOUBLE:
            value = cpl_sprintf("%.17e", cpl_parameter_get_double(param));
            default_value =
                cpl_sprintf("%.17e", cpl_parameter_get_default_double(param));
            break;

        case CPL_TYPE_STRING:
            tmp = er_json_escape_string(cpl_parameter_get_string(param));
            value = cpl_sprintf("\"%s\"", tmp);
            cpl_free(tmp);
            tmp = er_json_escape_string(
                                    cpl_parameter_get_default_string(param));
            default_value = cpl_sprintf("\"%s\"", tmp);
            cpl_free(tmp);
            break;

        default:
            cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                  "Type not supported for cpl_parameter.");
            return NULL;
    }

    const char * class_name = NULL;
    char * range_min = NULL;
    char * range_max = NULL;
    char * choices = NULL;
    int n = 0;
    switch (cpl_parameter_get_class(param))
    {
        case CPL_PARAMETER_CLASS_VALUE:
            class_name = "value";
            break;

        case CPL_PARAMETER_CLASS_RANGE:
            /* If parameter is a range then prepare the minimum and maximum
               values. Only integers and doubles are supported */
            class_name = "range";
            switch (value_type)
            {
                case CPL_TYPE_INT:
                    range_min = cpl_sprintf("%d",
                                    cpl_parameter_get_range_min_int(param));
                    range_max = cpl_sprintf("%d",
                                    cpl_parameter_get_range_max_int(param));
                    break;

                case CPL_TYPE_DOUBLE:
                    range_min = cpl_sprintf("%.17e",
                                    cpl_parameter_get_range_min_double(param));
                    range_max = cpl_sprintf("%.17e",
                                    cpl_parameter_get_range_max_double(param));
                    break;
                default:
                    cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                          "Type not supported for a range"
                                          " cpl_parameter.");
                    cpl_free(value);
                    cpl_free(default_value);
                    return NULL;
            }
            break;

        case CPL_PARAMETER_CLASS_ENUM:
            /* If parameter is an enum then prepare the list of choices.
               Only integers, doubles and strings are supported. */
            class_name = "enum";
            for (n = 0; n < cpl_parameter_get_enum_size(param); ++n)
            {
                /* For each enum we print the value and add it to the 'choices'
                   string. These must be comma separated. */
                switch (value_type)
                {
                    case CPL_TYPE_INT:
                        if (n == 0)
                        {
                            int val = cpl_parameter_get_enum_int(param, n);
                            choices = cpl_sprintf("[\n            %d", val);
                        }
                        else
                        {
                            int val = cpl_parameter_get_enum_int(param, n);
                            tmp = cpl_sprintf("%s,\n            %d", choices,
                                              val);
                            cpl_free(choices);
                            choices = tmp;
                        }
                        break;

                    case CPL_TYPE_DOUBLE:
                        if (n == 0)
                        {
                            double val = cpl_parameter_get_enum_double(param,
                                                                       n);
                            choices = cpl_sprintf("[\n            %.17e", val);
                        }
                        else
                        {
                            double val = cpl_parameter_get_enum_double(param,
                                                                       n);
                            tmp = cpl_sprintf("%s,\n            %.17e", choices,
                                              val);
                            cpl_free(choices);
                            choices = tmp;
                        }
                        break;

                    case CPL_TYPE_STRING:
                        if (n == 0)
                        {
                            char * val = er_json_escape_string(
                                    cpl_parameter_get_enum_string(param, n));
                            choices = cpl_sprintf("[\n            \"%s\"", val);
                            cpl_free(val);
                        }
                        else
                        {
                            char * val = er_json_escape_string(
                                    cpl_parameter_get_enum_string(param, n));
                            tmp = cpl_sprintf("%s,\n            \"%s\"",
                                                      choices, val);
                            cpl_free(choices);
                            cpl_free(val);
                            choices = tmp;
                        }
                        break;
                    default:
                        cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                              "Type not supported for an"
                                              " enumeration cpl_parameter.");
                        cpl_free(choices);
                        cpl_free(value);
                        cpl_free(default_value);
                        return NULL;
                }
            }
            if (choices == NULL)
            {
                cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                      "At least one enumeration choice must be"
                                      " given for cpl_parameter.");
                cpl_free(value);
                cpl_free(default_value);
                return NULL;
            }
            else
            {
                /* End the JSON list of choices with a ']'. */
                tmp = cpl_sprintf("%s\n          ]", choices);
                cpl_free(choices);
                choices = tmp;
            }
            break;

        default:
            cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                  "Parameter class not supported for"
                                  " cpl_parameter.");
            cpl_free(value);
            cpl_free(default_value);
            return NULL;
    }

    assert(class_name != NULL);

    /* Prepare the tag and alias strings. We have to take into account that
       these might be NULL. */
    char * tagstr = prepare_escaped_string(tag);
    char * cli_alias_str = prepare_escaped_string(cli_alias);
    char * env_alias_str = prepare_escaped_string(env_alias);
    char * cfg_alias_str = prepare_escaped_string(cfg_alias);

    /* Create the full JSON output text for the parameter object. */
    char * json = NULL;
    if (choices != NULL)
    {
        json = cpl_sprintf(
                "{\n"
                "        \"name\": \"%s\",\n"
                "        \"class\": \"%s\",\n"
                "        \"id\": %d,\n"
                "        \"description\": \"%s\",\n"
                "        \"context\": \"%s\",\n"
                "        \"value\": %s,\n"
                "        \"default\": %s,\n"
                "        \"choices\": %s,\n"
                "        \"present\": %s,\n"
                "        \"tag\": %s,\n"
                "        \"cli_enabled\": %s,\n"
                "        \"cli_alias\": %s,\n"
                "        \"env_enabled\": %s,\n"
                "        \"env_alias\": %s,\n"
                "        \"cfg_enabled\": %s,\n"
                "        \"cfg_alias\": %s\n"
                "      }",
                cpl_parameter_get_name(param), class_name, param_id,
                cpl_parameter_get_help(param), cpl_parameter_get_context(param),
                value, default_value, choices, present, tagstr,
                cli_enabled, cli_alias_str, env_enabled, env_alias_str,
                cfg_enabled, cfg_alias_str
            );
    }
    else if (range_min != NULL)
    {
        assert(range_max != NULL);
        json = cpl_sprintf(
                "{\n"
                "        \"name\": \"%s\",\n"
                "        \"class\": \"%s\",\n"
                "        \"id\": %d,\n"
                "        \"description\": \"%s\",\n"
                "        \"context\": \"%s\",\n"
                "        \"value\": %s,\n"
                "        \"default\": %s,\n"
                "        \"min\": %s,\n"
                "        \"max\": %s,\n"
                "        \"present\": %s,\n"
                "        \"tag\": %s,\n"
                "        \"cli_enabled\": %s,\n"
                "        \"cli_alias\": %s,\n"
                "        \"env_enabled\": %s,\n"
                "        \"env_alias\": %s,\n"
                "        \"cfg_enabled\": %s,\n"
                "        \"cfg_alias\": %s\n"
                "      }",
                cpl_parameter_get_name(param), class_name, param_id,
                cpl_parameter_get_help(param), cpl_parameter_get_context(param),
                value, default_value, range_min, range_max, present,
                tagstr, cli_enabled, cli_alias_str, env_enabled, env_alias_str,
                cfg_enabled, cfg_alias_str
            );
    }
    else
    {
        json = cpl_sprintf(
                "{\n"
                "        \"name\": \"%s\",\n"
                "        \"class\": \"%s\",\n"
                "        \"id\": %d,\n"
                "        \"description\": \"%s\",\n"
                "        \"context\": \"%s\",\n"
                "        \"value\": %s,\n"
                "        \"default\": %s,\n"
                "        \"present\": %s,\n"
                "        \"tag\": %s,\n"
                "        \"cli_enabled\": %s,\n"
                "        \"cli_alias\": %s,\n"
                "        \"env_enabled\": %s,\n"
                "        \"env_alias\": %s,\n"
                "        \"cfg_enabled\": %s,\n"
                "        \"cfg_alias\": %s\n"
                "      }",
                cpl_parameter_get_name(param), class_name, param_id,
                cpl_parameter_get_help(param), cpl_parameter_get_context(param),
                value, default_value, present, tagstr,
                cli_enabled, cli_alias_str, env_enabled, env_alias_str,
                cfg_enabled, cfg_alias_str
            );
    }

    /* Cleaning up allocated memory. NOTE: it is safe to pass NULL to cpl_free
       since it will check for these. */
    cpl_free(tagstr);
    cpl_free(cli_alias_str);
    cpl_free(env_alias_str);
    cpl_free(cfg_alias_str);
    cpl_free(range_min);
    cpl_free(range_max);
    cpl_free(choices);
    cpl_free(value);
    cpl_free(default_value);
    return json;
}

/*
 * Serialises a CPL frame object to properly formatted JSON text and returned
 * as a null terminated string.
 * If an error occurred then NULL is returned. Otherwise the result must be
 * cleaned up by the caller with cpl_free().
 */
static char * frame_to_json(const cpl_frame * frame)
{
    assert(frame != NULL);

    char * filename = er_json_escape_string(cpl_frame_get_filename(frame));
    char * tag = er_json_escape_string(cpl_frame_get_tag(frame));
    unsigned int type = (unsigned int) cpl_frame_get_type(frame);
    unsigned int group = (unsigned int) cpl_frame_get_group(frame);
    unsigned int level = (unsigned int) cpl_frame_get_level(frame);

    char * json = cpl_sprintf(
            "{\n"
            "        \"filename\": \"%s\",\n"
            "        \"tag\": \"%s\",\n"
            "        \"type\": %u,\n"
            "        \"group\": %u,\n"
            "        \"level\": %u\n"
            "      }",
            filename, tag, type, group, level
        );

    cpl_free(tag);
    cpl_free(filename);
    return json;
}

/**
 * @brief Serialises a CPL plugin object to JSON text.
 *
 * @param plugin The CPL plugin object to convert.
 *
 * @return The JSON test corresponding to the CPL plugin, otherwise @c NULL if
 *      an error occurred.
 *
 * This function will produce JSON text that corresponds to a CPL plugin object.
 * The following is an example of the produced JSON format:
 * <pre>
 * {
 *   "class": "Recipe",
 *   "name": "somerecipe",
 *   "version": 10205,
 *   "synopsis": "recipe",
 *   "description": "This is a recipe",
 *   "author": "Some Author",
 *   "email": "someone@eso.org",
 *   "copyright": "copyright message",
 *   "parameters": [
 *     {
 *       "name": "somerecipe.par1",
 *       "class": "value",
 *       "id": 0,
 *       "description": "integer parameter",
 *       "context": "somerecipe",
 *       "value": 12,
 *       "default": 5,
 *       "present": false,
 *       "tag": null,
 *       "cli_enabled": true,
 *       "cli_alias": "somerecipe.par1",
 *       "env_enabled": true,
 *       "env_alias": "somerecipe.par1",
 *       "cfg_enabled": true,
 *       "cfg_alias": "somerecipe.par1"
 *     }
 *   ],
 *   "frames": [
 *     {
 *       "filename": "test.fits",
 *       "tag": "RAW",
 *       "type": 2,
 *       "group": 1,
 *       "level": 3
 *     }
 *   ]
 * }
 * </pre>
 *
 * @note The caller must free the result with @c cpl_free() once it is no longer
 *      needed.
 */
char * er_plugin_to_json(const cpl_plugin * plugin)
{
    cpl_error_ensure(plugin != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Plugin is NULL.");

    /* Identify the needed parts of the recipe data that have to be serialised
       into JSON and sent to the Python interpreter. */
    cpl_parameterlist * params = NULL;
    cpl_frameset * frames = NULL;
    cpl_recipeconfig * config = NULL;
    unsigned long plugin_type = cpl_plugin_get_type(plugin);
    switch (plugin_type)
    {
        case CPL_PLUGIN_TYPE_RECIPE:
            {
                cpl_recipe * recipe = (cpl_recipe *) plugin;
                params = recipe->parameters;
                frames = recipe->frames;
            }
            break;

        case CPL_PLUGIN_TYPE_RECIPE_V2:
            {
                cpl_recipe2 * recipe = (cpl_recipe2 *) plugin;
                params = recipe->base.parameters;
                frames = recipe->base.frames;
                config = recipe->config;
            }
            break;

        default:
            cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                                  "Plugin type %lu is not supported.",
                                  plugin_type);
            return NULL;
    }

    /* Construct a properly formatted JSON array of parameters from the plugin's
       CPL parameter list. If no parameters are given then '[]' is used for the
       JSON text. */
    char * itemsjson = NULL;
    if (params != NULL)
    {
        const cpl_parameter * firstpar
                                    = cpl_parameterlist_get_first_const(params);
        const cpl_parameter * p = NULL;
        for (p = firstpar;
             p != NULL;
             p = cpl_parameterlist_get_next_const(params))
        {
            char * paramjson = parameter_to_json(p);
            if (p != firstpar)
            {
                char * tmp = cpl_sprintf("%s,\n      %s", itemsjson, paramjson);
                cpl_free(paramjson);
                cpl_free(itemsjson);
                itemsjson = tmp;
            }
            else
            {
                itemsjson = paramjson;
            }
        }
    }
    char * parlistjson;
    if (itemsjson != NULL)
    {
        parlistjson = cpl_sprintf("[\n      %s\n    ]", itemsjson);
    }
    else
    {
        parlistjson = cpl_strdup("[]");
    }
    cpl_free(itemsjson);

    /* Construct a properly formatted JSON array of frames from the plugin's CPL
       frameset object. If there are not frames then '[]' is used for the JSON
       text. */
    itemsjson = NULL;
    if (frames != NULL)
    {
        cpl_frameset_iterator * frameiter = cpl_frameset_iterator_new(frames);
        const cpl_frame * firstframe
                                = cpl_frameset_iterator_get_const(frameiter);
        const cpl_frame * frame = NULL;
        for (frame = firstframe;
             (frame = cpl_frameset_iterator_get_const(frameiter)) != NULL;
             cpl_frameset_iterator_advance(frameiter, 1))
        {
            char * framejson = frame_to_json(frame);
            if (frame != firstframe)
            {
                char * tmp = cpl_sprintf("%s,\n      %s", itemsjson, framejson);
                cpl_free(framejson);
                cpl_free(itemsjson);
                itemsjson = tmp;
            }
            else
            {
                itemsjson = framejson;
            }
        }
        cpl_frameset_iterator_delete(frameiter);
    }
    char * framelistjson;
    if (itemsjson != NULL)
    {
        framelistjson = cpl_sprintf("[\n      %s\n    ]", itemsjson);
    }
    else
    {
        framelistjson = cpl_strdup("[]");
    }
    cpl_free(itemsjson);

    /* Convert the recipe configuration object into JSON text if we are dealing
       with a recipe v2 structure. */
    char * configjson = NULL;
    if (config != NULL)
    {
        itemsjson = NULL;

        char ** tags = cpl_recipeconfig_get_tags(config);
        if (tags == NULL)
        {
            cpl_free(parlistjson);
            cpl_free(framelistjson);
            return NULL;
        }
        char ** tag = NULL;
        for (tag = tags; *tag != NULL; ++tag)
        {
            char * tmp;
            cpl_size min, max;

            /* Encode the primary tag's input list as a JSON array of
               objects. */
            char * inputitems = NULL;
            char ** inputs = cpl_recipeconfig_get_inputs(config, *tag);
            if (inputs != NULL)
            {
                inputitems = cpl_strdup("[");
                char ** input = NULL;
                for (input = inputs; *input != NULL; ++input)
                {
                    min = cpl_recipeconfig_get_min_count(config, *tag, *input);
                    max = cpl_recipeconfig_get_max_count(config, *tag, *input);
                    char * inputstr = er_json_escape_string(*input);
                    cpl_boolean first_item = (input == inputs);
                    tmp = cpl_sprintf(
                                "%s%s\n"
                                "          {\n"
                                "            \"tag\": \"%s\",\n"
#if defined CPL_SIZE_BITS && CPL_SIZE_BITS == 32
                                "            \"min\": %d,\n"
                                "            \"max\": %d\n"
#else
                                "            \"min\": %lld,\n"
                                "            \"max\": %lld\n"
#endif
                                "          }",
                                inputitems, (first_item ? "" : ","), inputstr,
                                min, max
                            );
                    cpl_free(inputstr);
                    cpl_free(inputitems);
                    inputitems = tmp;
                    cpl_free(*input);
                }
                cpl_free(inputs);
                tmp = cpl_sprintf("%s\n        ]", inputitems);
                cpl_free(inputitems);
                inputitems = tmp;
            }

            /* Encode the primary tag's output list as a JSON array of
               strings. */
            char * outputitems = NULL;
            char ** outputs = cpl_recipeconfig_get_outputs(config, *tag);
            if (outputs != NULL)
            {
                outputitems = cpl_strdup("[");
                char ** output = NULL;
                for (output = outputs; *output != NULL; ++output)
                {
                    char * outputstr = er_json_escape_string(*output);
                    cpl_boolean first_item = (output == outputs);
                    tmp = cpl_sprintf("%s%s\n          \"%s\"", outputitems,
                                      (first_item ? "" : ","), outputstr);
                    cpl_free(outputstr);
                    cpl_free(outputitems);
                    outputitems = tmp;
                    cpl_free(*output);
                }
                cpl_free(outputs);
                tmp = cpl_sprintf("%s\n        ]", outputitems);
                cpl_free(outputitems);
                outputitems = tmp;
            }

            /* If not the first tag entry in the list then add a comma and
               space prefix. Otherwise we create a new empty string for
               itemsjson. */
            if (tag != tags)
            {
                tmp = cpl_sprintf("%s,\n      ", itemsjson);
                cpl_free(itemsjson);
                itemsjson = tmp;
            }
            else
            {
                itemsjson = cpl_strdup("");
            }

            min = cpl_recipeconfig_get_min_count(config, *tag, *tag);
            max = cpl_recipeconfig_get_max_count(config, *tag, *tag);
            char * tagstr = er_json_escape_string(*tag);

            /* Start a new JSON object for the current tag. */
            tmp = cpl_sprintf(
                        "%s{\n"
                        "        \"tag\": \"%s\",\n"
#if defined CPL_SIZE_BITS && CPL_SIZE_BITS == 32
                        "        \"min\": %d,\n"
                        "        \"max\": %d",
#else
                        "        \"min\": %lld,\n"
                        "        \"max\": %lld",
#endif
                        itemsjson, tagstr, min, max
                    );
            cpl_free(tagstr);
            cpl_free(itemsjson);
            itemsjson = tmp;

            /* Add the inputs if any are available. */
            if (inputitems != NULL)
            {
                tmp = cpl_sprintf("%s,\n        \"inputs\": %s",
                                  itemsjson, inputitems);
                cpl_free(itemsjson);
                itemsjson = tmp;
                cpl_free(inputitems);
            }

            /* Add the outputs if any are available. */
            if (outputitems != NULL)
            {
                tmp = cpl_sprintf("%s,\n        \"outputs\": %s",
                                  itemsjson, outputitems);
                cpl_free(itemsjson);
                itemsjson = tmp;
                cpl_free(outputitems);
            }

            /* Complete the JSON object. */
            tmp = cpl_sprintf("%s\n      }", itemsjson);
            cpl_free(itemsjson);
            itemsjson = tmp;

            cpl_free(*tag);
        }
        cpl_free(tags);

        /* Create the final top level array if any tag items were available.
           Otherwise just use and empty JSON array. */
        if (itemsjson != NULL)
        {
            configjson = cpl_sprintf("[\n      %s\n    ]", itemsjson);
        }
        else
        {
            configjson = cpl_strdup("[]");
        }
        cpl_free(itemsjson);
    }

    char * name = er_json_escape_string(cpl_plugin_get_name(plugin));
    unsigned long version = cpl_plugin_get_version(plugin);
    char * synopsis = er_json_escape_string(cpl_plugin_get_synopsis(plugin));
    char * desc = er_json_escape_string(cpl_plugin_get_description(plugin));
    char * author = er_json_escape_string(cpl_plugin_get_author(plugin));
    char * email = er_json_escape_string(cpl_plugin_get_email(plugin));
    char * copyright = er_json_escape_string(cpl_plugin_get_copyright(plugin));

    char * json;
    if (configjson == NULL)
    {
        json = cpl_sprintf(
                "{\n"
                "    \"class\": \"unknown\",\n"
                "    \"name\": \"%s\",\n"
                "    \"version\": %lu,\n"
                "    \"synopsis\": \"%s\",\n"
                "    \"description\": \"%s\",\n"
                "    \"author\": \"%s\",\n"
                "    \"email\": \"%s\",\n"
                "    \"copyright\": \"%s\",\n"
                "    \"parameters\": %s,\n"
                "    \"frames\": %s\n"
                "  }",
                name, version, synopsis, desc, author, email, copyright,
                parlistjson, framelistjson
            );
    }
    else
    {
        json = cpl_sprintf(
                "{\n"
                "    \"class\": \"unknown\",\n"
                "    \"name\": \"%s\",\n"
                "    \"version\": %lu,\n"
                "    \"synopsis\": \"%s\",\n"
                "    \"description\": \"%s\",\n"
                "    \"author\": \"%s\",\n"
                "    \"email\": \"%s\",\n"
                "    \"copyright\": \"%s\",\n"
                "    \"parameters\": %s,\n"
                "    \"frames\": %s,\n"
                "    \"recipeconfig\": %s\n"
                "  }",
                name, version, synopsis, desc, author, email, copyright,
                parlistjson, framelistjson, configjson
            );
    }

    cpl_free(copyright);
    cpl_free(email);
    cpl_free(author);
    cpl_free(desc);
    cpl_free(synopsis);
    cpl_free(name);
    cpl_free(parlistjson);
    cpl_free(framelistjson);
    cpl_free(configjson);

    return json;
}


/**
 * @struct er_json_node
 * A JSON node object that forms part of the parse tree produced by
 * er_json_parse().
 */
struct _er_json_node_
{
    er_json_type type;    /* Identifies which union member is active. */
    union
    {
        cx_map * object;
        cx_list * array;
        char * string;
        double number;
        cpl_boolean bool;
    };
    const char * location;   /* Corresponding position in the JSON text. */
};

/**
 * @struct er_json_array_iterator
 * Iterator pointer for navigating an er_json_node structure that corresponds to
 * a JSON array, i.e. er_json_node_type() returns @c JSON_ARRAY for the node.
 */

/**
 * @struct er_json_object_iterator
 * Iterator pointer for navigating an er_json_node structure that corresponds to
 * a JSON object, i.e. er_json_node_type() returns @c JSON_OBJECT for the node.
 */

/*
 * NOTE: in the following functions we only use cpl_error_ensure() in the ones
 * that are exported, i.e. that are accessible through the public API.
 * In the static functions we instead use assert() to check preconditions. This
 * will only test the code in debug builds. This is fine, since the static
 * functions are not reachable outside this file.
 */

/*
 * Create new null type JSON node.
 * @param location  The corresponding position in the JSON text this object was
 *                  parsed from.
 * @return New node object. cpl_malloc() will abort if out of memory.
 */
static er_json_node * er_json_node_new_null(const char * location)
{
    assert(location != NULL);
    er_json_node * obj = cpl_malloc(sizeof(er_json_node));
    obj->type = JSON_NULL;
    obj->object = NULL;
    obj->location = location;
    return obj;
}


static cxint _string_key_compare(cxcptr a, cxcptr b)
{
    return (strcmp(a, b) < 0) ? TRUE : FALSE;
}

/*
 * Create new object type JSON node.
 * @param location  The corresponding position in the JSON text this object was
 *                  parsed from.
 * @return New node corresponding to a mapping/dictionary. cpl_malloc() will
 *         abort if out of memory.
 */
static er_json_node * er_json_node_new_object(const char * location)
{
    assert(location != NULL);
    er_json_node * obj = cpl_malloc(sizeof(er_json_node));
    obj->type = JSON_OBJECT;
    obj->object = cx_map_new(_string_key_compare, cpl_free,
                             (cx_free_func)er_json_node_delete);
    obj->location = location;
    return obj;
}

/*
 * Create new array type JSON node.
 * @param location  The corresponding position in the JSON text this object was
 *                  parsed from.
 * @return New node corresponding to a JSON array. cpl_malloc() will abort if
 *         out of memory.
 */
static er_json_node * er_json_node_new_array(const char * location)
{
    assert(location != NULL);
    er_json_node * obj = cpl_malloc(sizeof(er_json_node));
    obj->type = JSON_ARRAY;
    obj->array = cx_list_new();
    obj->location = location;
    return obj;
}

/*
 * Create new string type JSON node.
 * @param string    The null terminated string value for this node.
 * @param location  The corresponding position in the JSON text this object was
 *                  parsed from.
 * @return New node object. cpl_malloc() will abort if out of memory.
 */
static er_json_node * er_json_node_new_string(const char * string,
                                              const char * location)
{
    assert(location != NULL);
    er_json_node * obj = cpl_malloc(sizeof(er_json_node));
    obj->type = JSON_STRING;
    obj->string = cpl_strdup(string);
    obj->location = location;
    return obj;
}

/*
 * Create new number type JSON node.
 * @param number    Floating point number for this node.
 * @param location  The corresponding position in the JSON text this object was
 *                  parsed from.
 * @return New node storing a double. cpl_malloc() will abort if out of memory.
 */
static er_json_node * er_json_node_new_number(double number,
                                              const char * location)
{
    assert(location != NULL);
    er_json_node * obj = cpl_malloc(sizeof(er_json_node));
    obj->type = JSON_NUMBER;
    obj->number = number;
    obj->location = location;
    return obj;
}

/*
 * Create new boolean type JSON node.
 * @param bool      Boolean value, either CPL_TRUE or CPL_FALSE.
 * @param location  The corresponding position in the JSON text this object was
 *                  parsed from.
 * @return New boolean node object. cpl_malloc() will abort if out of memory.
 */
static er_json_node * er_json_node_new_bool(cpl_boolean bool,
                                            const char * location)
{
    assert(location != NULL);
    er_json_node * obj = cpl_malloc(sizeof(er_json_node));
    obj->type = JSON_BOOL;
    obj->bool = bool;
    obj->location = location;
    return obj;
}

/**
 * @brief Delete a parse tree JSON node returned by er_json_parse().
 *
 * @param obj  The parse tree object to be deleted.
 *
 * This function must be called on the resultant object from a call to the
 * er_json_parse() function when it is no longer needed.
 * This will free the allocated buffers.
 */
void er_json_node_delete(er_json_node * obj)
{
    if (obj == NULL) return;
    switch (obj->type)
    {
        case JSON_NULL:
        case JSON_NUMBER:
        case JSON_BOOL:
            /* Nothing to be done. */
            break;
        case JSON_OBJECT:
            cx_map_delete(obj->object);
            break;
        case JSON_ARRAY:
            cx_list_destroy(obj->array, (cx_free_func)er_json_node_delete);
            break;
        case JSON_STRING:
            cpl_free(obj->string);
            break;
        default:
            assert(
                obj->type == JSON_NULL ||
                obj->type == JSON_OBJECT ||
                obj->type == JSON_ARRAY ||
                obj->type == JSON_STRING ||
                obj->type == JSON_NUMBER ||
                obj->type == JSON_BOOL
            );
    }
    cpl_free(obj);
}

/**
 * @brief Return the type of the JSON node.
 *
 * @param obj  The JSON node to query.
 *
 * @return The type code for this JSON node or @c JSON_NULL if @p obj is
 *         invalid. Must use @c cpl_error_get_code() to confirm if it was
 *         invalid.
 *
 * This returns the type code that can be used to figure out what kind of JSON
 * node one is dealing with. This is important to decide which functions can be
 * used with this node. e.g. the array functions or object functions etc.
 * Valid type codes returned by this function are one of the following:
 * <ul>
 *   <li>JSON_NULL</li>
 *   <li>JSON_BOOL</li>
 *   <li>JSON_NUMBER</li>
 *   <li>JSON_STRING</li>
 *   <li>JSON_ARRAY</li>
 *   <li>JSON_OBJECT</li>
 * </ul>
 */
er_json_type er_json_node_type(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return JSON_NULL,
                     "JSON node is NULL.");
    return obj->type;
}

/**
 * @brief Returns the node's location in the original JSON text.
 *
 * @param obj  The JSON node to query.
 *
 * @return The string position within the original JSON text or @c NULL if this
 *         node is not valid. A CPL error is set if an error occurs.
 *
 * Returns the string location within the original JSON text that this JSON node
 * was parsed from. Allows to calculate the line and column numbers for error
 * messages using er_json_find_line_column().
 */
const char * er_json_node_location(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    return obj->location;
}

/**
 * @brief Return the boolean value for a JSON node.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_BOOL.
 *
 * @return The @c CPL_TRUE or @c CPL_FALSE value. If this was not a valid node
 *         then @c CPL_FALSE is returned by default, @c cpl_error_get_code()
 *         can be used to check for such errors.
 *
 * For boolean JSON nodes this function will return its value as a
 * @c cpl_boolean.
 */
cpl_boolean er_json_node_get_bool(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return CPL_FALSE,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_BOOL, CPL_ERROR_ILLEGAL_INPUT,
                     return CPL_FALSE, "JSON node is not a boolean.");
    return obj->bool;
}

/**
 * @brief Returns the number value for a JSON node.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_NUMBER.
 *
 * @return The floating point value or NaN (Not a Number) if this was not a
 *         valid node.  @c cpl_error_get_code() can be used to check if the node
 *         really contains a NaN value or if it was in fact invalid.
 *
 * For a number type JSON node this function will return its value as a double
 * precision floating point number.
 */
double er_json_node_get_number(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NAN,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_NUMBER, CPL_ERROR_ILLEGAL_INPUT,
                     return NAN, "JSON node is not a number.");
    return obj->number;
}

/**
 * @brief Returns the string value for a JSON node.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_STRING.
 *
 * @return The string value or @c NULL if the node is not valid.
 *         A CPL error is set if an error occurs.
 *
 * For a string type JSON node this function will return its value as a null
 * terminated character string.
 */
const char * er_json_node_get_string(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_STRING, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not a string.");
    return obj->string;
}

/**
 * @brief Returns the size of an array node.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_ARRAY.
 *
 * @return The number of elements in the array or @c -1 if the node is not
 *         valid. A CPL error is set if an error occurs.
 *
 * For array type JSON nodes this will return the number of elements in the
 * array.
 */
cpl_size er_json_node_array_size(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return -1,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_ARRAY, CPL_ERROR_ILLEGAL_INPUT,
                     return -1, "JSON node is not an array.");
    return cx_list_size(obj->array);
}

/**
 * @brief Checks if the array node is empty or not.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_ARRAY.
 *
 * @return @c CPL_TRUE if the array is empty or @c CPL_FALSE otherwise. If the
 *         node is not valid @c CPL_FALSE is returned and one should check this
 *         error condition with @c cpl_error_get_code().
 *
 * For array type JSON nodes this function will test if the array is empty.
 */
cpl_boolean er_json_node_array_empty(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return CPL_TRUE,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_ARRAY, CPL_ERROR_ILLEGAL_INPUT,
                     return CPL_TRUE, "JSON node is not an array.");
    return cx_list_empty(obj->array);
}

/**
 * @brief Gets an iterator to the first element of the array.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_ARRAY.
 *
 * @return An iterator pointer or @c NULL if this node is not valid.
 *         A CPL error is set if an error occurs.
 *
 * For an array type node this function returns an iterator to the first element
 * of the array.
 */
er_json_array_iterator er_json_node_array_begin(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_ARRAY, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an array.");
    return cx_list_begin(obj->array);
}

/**
 * @brief Gets an iterator pointing to one past the last element of the array.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_ARRAY.
 *
 * @return An iterator pointer or @c NULL if this node is not valid.
 *         A CPL error is set if an error occurs.
 *
 * For an array type node this function returns an iterator marking a position
 * that is one element past the end of the array. This allows to construct the
 * following kind of iteration loop:
 *
 * @code
 *   for (er_json_array_iterator iter = er_json_node_array_begin(obj);
 *        iter != er_json_node_array_end(obj);
 *        iter = er_json_node_array_next(obj, iter))
 *   {
 *      ...
 *   }
 * @endcode
 */
er_json_array_iterator er_json_node_array_end(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_ARRAY, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an array.");
    return cx_list_end(obj->array);
}

/**
 * @brief Moves to the next position in an array.
 *
 * @param obj  The JSON node to from which @p current was created.
 *             Must be a node with er_json_node_type() returning @c JSON_ARRAY.
 * @param current  The iterator object being incremented.
 *
 * @return The incremented iterator or @c NULL if the current iterator or array
 *         node is invalid. A CPL error is set if an error occurs.
 *
 * For array type nodes this will increment an iterator object that was created
 * from the node. i.e. move the iterator to the next position in the array.
 */
er_json_array_iterator er_json_node_array_next(const er_json_node * obj,
                                               er_json_array_iterator current)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(current != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Iterator is NULL.");
    cpl_error_ensure(obj->type == JSON_ARRAY, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an array.");
    return cx_list_next(obj->array, current);
}

/**
 * @brief Moves to the previous position in an array.
 *
 * @param obj  The JSON node to from which @p current was created.
 *             Must be a node with er_json_node_type() returning @c JSON_ARRAY.
 * @param current  The iterator object being decremented.
 *
 * @return The decremented iterator or @c NULL if the current iterator or array
 *         node is invalid. A CPL error is set if an error occurs.
 *
 * For array type nodes this will decrement an iterator object that was created
 * from the node. i.e. move the iterator to the previous position in the array.
 */
er_json_array_iterator er_json_node_array_previous(
                                                const er_json_node * obj,
                                                er_json_array_iterator current)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(current != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Iterator is NULL.");
    cpl_error_ensure(obj->type == JSON_ARRAY, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an array.");
    return cx_list_previous(obj->array, current);
}

/**
 * @brief Returns the corresponding array item pointed to by an iterator.
 *
 * @param obj  The JSON node to from which @p iterator was created.
 *             Must be a node with er_json_node_type() returning @c JSON_ARRAY.
 * @param iterator  The iterator object pointing to the item of interest.
 *
 * @return The child JSON node item or @c NULL if the iterator or parent array
 *         node is invalid. A CPL error is set if an error occurs.
 *
 * For array type nodes this function returns the child item pointed to by the
 * given iterator. The iterator must have been created from the given array
 * node.
 */
const er_json_node * er_json_node_array_get(const er_json_node * obj,
                                            er_json_array_iterator iterator)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(iterator != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Iterator is NULL.");
    cpl_error_ensure(obj->type == JSON_ARRAY, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an array.");
    return (er_json_node *) cx_list_get(obj->array, iterator);
}

/*
 * Adds a JSON node to an existing array node object.
 * @param obj   The array node to modify. Must be of type JSON_ARRAY.
 * @param value The node being added.
 *
 * NOTE: Takes ownership of value. i.e. the caller should no longer call
 * er_json_node_delete() for value after invoking this function.
 */
static void er_json_node_array_append(er_json_node * obj, er_json_node * value)
{
    assert(obj != NULL);
    assert(obj->type == JSON_ARRAY);
    cx_list_push_back(obj->array, value);
}

/**
 * @brief Returns the size of an object node.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 *
 * @return The number of attributes in the object or @c -1 if the node is not
 *         valid. A CPL error is set if an error occurs.
 *
 * For object type JSON nodes this will return the number of attributes for the
 * object.
 */
cpl_size er_json_node_object_size(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return -1,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return -1, "JSON node is not an object.");
    return cx_map_size(obj->object);
}

/**
 * @brief Checks if the object node is empty or not.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 *
 * @return @c CPL_TRUE if the object has no attributes or @c CPL_FALSE
 *         otherwise. If the node is not valid @c CPL_FALSE is returned and one
 *         should check this error condition with @c cpl_error_get_code().
 *
 * For object type JSON nodes this function will test if the object contains any
 * attributes.
 */
cpl_boolean er_json_node_object_empty(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return CPL_TRUE,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return CPL_TRUE, "JSON node is not an object.");
    return cx_map_empty(obj->object);
}

/**
 * @brief Gets an iterator to the first attribute of the object.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 *
 * @return An iterator pointer or @c NULL if this node is not valid.
 *         A CPL error is set if an error occurs.
 *
 * For an object type node this function returns an iterator to the first
 * attribute of the object.
 */
er_json_object_iterator er_json_node_object_begin(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an object.");
    return cx_map_begin(obj->object);
}

/**
 * @brief Gets an iterator pointing to one past the last attribute of the
 *        object.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 *
 * @return An iterator pointer or @c NULL if this node is not valid.
 *         A CPL error is set if an error occurs.
 *
 * For an object type node this function returns an iterator marking a position
 * that is one attribute past the end of the object's attribute list. This
 * allows to construct the following kind of iteration loop:
 *
 * @code
 *   for (er_json_object_iterator iter = er_json_node_object_begin(obj);
 *        iter != er_json_node_object_end(obj);
 *        iter = er_json_node_object_next(obj, iter))
 *   {
 *      ...
 *   }
 * @endcode
 */
er_json_object_iterator er_json_node_object_end(const er_json_node * obj)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an object.");
    return cx_map_end(obj->object);
}

/**
 * @brief Moves to the next attribute in an object.
 *
 * @param obj  The JSON node to from which @p current was created.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 * @param current  The iterator object being incremented.
 *
 * @return The incremented iterator or @c NULL if the current iterator or object
 *         node is invalid. A CPL error is set if an error occurs.
 *
 * For object type nodes this will increment an iterator object that was created
 * from the node. i.e. move the iterator to the next position in the object's
 * attribute list.
 */
er_json_object_iterator er_json_node_object_next(
                                            const er_json_node * obj,
                                            er_json_object_iterator current)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(current != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Iterator is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an object.");
    return cx_map_next(obj->object, current);
}

/**
 * @brief Moves to the previous attribute in an object.
 *
 * @param obj  The JSON node to from which @p current was created.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 * @param current  The iterator object being decremented.
 *
 * @return The decremented iterator or @c NULL if the current iterator or object
 *         node is invalid. A CPL error is set if an error occurs.
 *
 * For object type nodes this will decrement an iterator object that was created
 * from the node. i.e. move the iterator to the previous position in the
 * object's attribute list.
 */
er_json_object_iterator er_json_node_object_previous(
                                            const er_json_node * obj,
                                            er_json_object_iterator current)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(current != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Iterator is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an object.");
    return cx_map_previous(obj->object, current);
}

/**
 * @brief Returns the corresponding attribute name pointed to by an iterator.
 *
 * @param obj  The JSON node to from which @p iterator was created.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 * @param iterator  The iterator object pointing to the attribute of interest.
 *
 * @return The attribute name as a string, or @c NULL if the iterator or parent
 *         object node is invalid. A CPL error is set if an error occurs.
 *
 * For object type nodes this function returns the key or attribute name pointed
 * to by the given iterator. The iterator must have been created from the given
 * object node.
 */
const char * er_json_node_object_get_key(const er_json_node * obj,
                                         er_json_object_iterator iterator)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(iterator != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Iterator is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an object.");
    return (const char *) cx_map_get_key(obj->object, iterator);
}

/**
 * @brief Returns the corresponding attribute value pointed to by an iterator.
 *
 * @param obj  The JSON node to from which @p iterator was created.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 * @param iterator  The iterator object pointing to the attribute of interest.
 *
 * @return The attribute's value node object, or @c NULL if the iterator or
 *         parent node is invalid. A CPL error is set if an error occurs.
 *
 * For object type nodes this function returns the attribute's value, pointed to
 * by the given iterator. The iterator must have been created from the given
 * object node.
 */
const er_json_node * er_json_node_object_get_value(
                                            const er_json_node * obj,
                                            er_json_object_iterator iterator)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(iterator != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Iterator is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an object.");
    return (er_json_node *) cx_map_get_value(obj->object, iterator);
}

/**
 * @brief Finds the value for a particular attribute in an object node.
 *
 * @param obj  The JSON node to query.
 *             Must be a node with er_json_node_type() returning @c JSON_OBJECT.
 * @param key  The name of the attribute to search for in the object node.
 *
 * @return The attribute's value node or @c NULL if no such attribute exists.
 *         If the parent node being queries is not valid then @c NULL is also
 *         returned and an error code is set, which can be queried with
 *         @c cpl_error_get_code().
 *
 * This function will try to find an attribute in an object type JSON node. The
 * attribute is named by the given key. If such an attribute is found it's
 * corresponding value is returned or @c NULL if no such attribute can be found.
 */
const er_json_node * er_json_node_object_get(const er_json_node * obj,
                                             const char * key)
{
    cpl_error_ensure(obj != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON node is NULL.");
    cpl_error_ensure(key != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Key is NULL.");
    cpl_error_ensure(obj->type == JSON_OBJECT, CPL_ERROR_ILLEGAL_INPUT,
                     return NULL, "JSON node is not an object.");
    er_json_object_iterator item = cx_map_find(obj->object, key);
    if (item != cx_map_end(obj->object))
    {
        return (er_json_node *) cx_map_get_value(obj->object, item);
    }
    else
    {
        return NULL;
    }
}

/*
 * Adds a new key/value pair to the JSON object node.
 * @param obj    The node to modify. Must be of type JSON_OBJECT.
 * @param key    The key string for the new entry.
 * @param value  The value node for the new entry.
 *
 * NOTE: Takes ownership of value. i.e. the caller should no longer call
 * er_json_node_delete() for value after invoking this function.
 */
static void er_json_node_object_insert(er_json_node * obj,
                                       const char * key, er_json_node * value)
{
    assert(obj != NULL);
    assert(obj->type == JSON_OBJECT);
    cx_map_insert(obj->object, cpl_strdup(key), value);
}


/* Structure to store the context for the parsing routines. It is easier to pass
   one pointer as argument to each function rather than all the values below.
 */
typedef struct _json_parser_context_
{
    const char * json;     /* Original JSON text to parse. */
    char * buffer;         /* Working buffer for parsing functions.
                              Will be big enough to fit the remainder of the
                              JSON string pointed to by 'current'. */
    const char * current;  /* The current parse position. */
    char * error;          /* The error message to use if JSON is invalid. */
    er_json_node * tree;     /* The parse tree built by the parser. */
} json_parser_context;

/*
 * Skips any white space from the current parsing position.
 * @param ctx  Pointer to the parsing context data for the current parse.
 */
static void parse_whitespace(json_parser_context * ctx)
{
    assert(ctx != NULL);
    while (*ctx->current != '\0' && isspace(*ctx->current)) ++ctx->current;
}

/*
 * Parses a 'null' JSON terminal from the current parsing position.
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return A new node object corresponding to the null terminal. NULL is
 *         returned on parse errors.
 */
static er_json_node * parse_null(json_parser_context * ctx)
{
    assert(ctx != NULL);
    er_json_node * value = NULL;
    const char * nullstr = "null";
    if (strncmp(ctx->current, nullstr, strlen(nullstr)) == 0)
    {
        value = er_json_node_new_null(ctx->current);
        ctx->current += strlen(nullstr);
    }
    else
    {
        if (strlen(ctx->current) > 20)
        {
            ctx->error = cpl_sprintf("Expected null but found '%.20s...'",
                                     ctx->current);
        }
        else
        {
            ctx->error = cpl_sprintf("Expected null but found '%.20s'",
                                     ctx->current);
        }
    }
    return value;
}

/*
 * Parses a boolean JSON terminal from the current parse position.
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return A new node object corresponding to the boolean terminal. NULL is
 *         returned on parse errors.
 */
static er_json_node * parse_boolean(json_parser_context * ctx)
{
    assert(ctx != NULL);
    er_json_node * value = NULL;
    const char * truestr = "true";
    const char * falsestr = "false";
    if (strncmp(ctx->current, truestr, strlen(truestr)) == 0)
    {
        value = er_json_node_new_bool(CPL_TRUE, ctx->current);
        ctx->current += strlen(truestr);
    }
    else if (strncmp(ctx->current, falsestr, strlen(falsestr)) == 0)
    {
        value = er_json_node_new_bool(CPL_FALSE, ctx->current);
        ctx->current += strlen(falsestr);
    }
    else
    {
        if (strlen(ctx->current) > 20)
        {
            ctx->error = cpl_sprintf("Expected true or false but found"
                                     " '%.20s...'", ctx->current);
        }
        else
        {
            ctx->error = cpl_sprintf("Expected true or false but found '%.20s'",
                                     ctx->current);
        }
    }
    return value;
}

/*
 * Parses a JSON number from the current parse position.
 * This actually accepts a larger set of numbers than strict JSON. NaN and
 * infinity are also accepted.
 *
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return A new number node corresponding to the null terminal. NULL is
 *         returned on parse errors.
 */
static er_json_node * parse_number(json_parser_context * ctx)
{
    assert(ctx != NULL);
    er_json_node * value = NULL;
    const char * special_value_str[8] = {
        /* NOTE: The order of these strings matters. Must have substrings after
           the longer strings. i.e. Infinity after Inf etc. */
        "nan", "NaN", "inf", "Infinity", "Inf", "-inf", "-Infinity", "-Inf"
    };
    double special_value[8] = {
        NAN, NAN, INFINITY, INFINITY, INFINITY, -INFINITY, -INFINITY, -INFINITY
    };
    int i = 0;
    for (i = 0; i < 8; ++i)
    {
        const char * str = special_value_str[i];
        if (strncmp(ctx->current, str, strlen(str)) == 0)
        {
            value = er_json_node_new_number(special_value[i], ctx->current);
            ctx->current += strlen(str);
            return value;
        }
    }

    errno = 0;  /* Clear errno before checking strtod for errors. */
    char * endptr = NULL;
    double number = strtod(ctx->current, &endptr);
    /* NOTE: we ignore floating point overflow and underflows (ERANGE). */
    if (errno == 0 || errno == ERANGE)
    {
        value = er_json_node_new_number(number, ctx->current);
        ctx->current = endptr;
    }
    else
    {
        if (strlen(ctx->current) > 20)
        {
            ctx->error = cpl_sprintf("Expected a number but found '%.20s...'",
                                     ctx->current);
        }
        else
        {
            ctx->error = cpl_sprintf("Expected number but found '%.20s'",
                                     ctx->current);
        }
    }
    return value;
}

/*
 * Parses the hexadecimal number for a Unicode escape sequence from the current
 * parse position.
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return The corresponding character code or -1 on a parse error.
 */
static long parse_unicode_character(json_parser_context * ctx, char * buffer)
{
    assert(ctx != NULL);
    assert(buffer != NULL);

    /* Write a zero to the buffer to ensure there will be a null byte
       terminating the string. Then copy the 4 characters from the current JSON
       parse location into the buffer. */
    buffer[4] = '\0';
    strncpy(buffer, ctx->current, 4);

    if (strlen(buffer) != 4)
    {
        ctx->error = cpl_strdup("A unicode escape sequence requires four"
                                " hexadecimal digits");
        return -1;
    }

    /* Decode the hex digits. */
    errno = 0;  /* Clear errno before checking strtol for errors. */
    char * endptr = NULL;
    long charcode = strtol(buffer, &endptr, 16);
    if (errno != 0 || endptr != buffer + 4)
    {
        if (strlen(ctx->current) > 20)
        {
            ctx->error = cpl_sprintf("Expected a four character hexadecimal"
                                     " number but found '%.20s...'",
                                     ctx->current);
        }
        else
        {
            ctx->error = cpl_sprintf("Expected a four character hexadecimal"
                                     " number but found '%.20s'", ctx->current);
        }
        return -1;
    }

    if (charcode == 0 || charcode > 255)
    {
        ctx->error = cpl_strdup("Unicode value is outside the supported range");
        return -1;
    }

    /* NOTE: Because this will be called in the for loop inside the
       parse_string_characters function we only increment by 3 characters,
       since it will be again incremented in the for loop. */
    ctx->current += 3;
    return charcode;
}

/*
 * Parses and converts escape sequences in a JSON string from the current parse
 * position. Will stop at the end of string character ("), i.e. the one that is
 * not part of an escape sequence.
 *
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return The converted null terminated string or NULL on parse errors.
 */
static const char * parse_string_characters(json_parser_context * ctx)
{
    assert(ctx != NULL);
    assert(ctx->buffer != NULL);

    char * result = ctx->buffer;

    /* The following variable keeps track of the current character position to
       write to in the output buffer. */
    char * output = result;

    for (; *ctx->current != '\0' && *ctx->current != '"'; ++ctx->current)
    {
        char char_to_use;

        /* Check if we found the start of an escape sequence. */
        if (*ctx->current != '\\')
        {
            char_to_use = *ctx->current;
        }
        else
        {
            long charcode;
            ++ctx->current;  /* Move past the '\' character. */
            if (*ctx->current == '\0')
            {
                ctx->error = cpl_strdup("Missing escape character after '\\'");
                return NULL;
            }
            switch (*ctx->current)
            {
                case '"': char_to_use = '"'; break;
                case '\\': char_to_use = '\\'; break;
                case '/': char_to_use = '/'; break;
                case 'b': char_to_use = '\b'; break;
                case 'f': char_to_use = '\f'; break;
                case 'n': char_to_use = '\n'; break;
                case 'r': char_to_use = '\r'; break;
                case 't': char_to_use = '\t'; break;

                case 'u':
                    ++ctx->current;  /* Move past the 'u' character. */
                    /* NOTE: we let parse_unicode_character use the buffer space
                       just after out current output location, which will not
                       jet contain any usable data. */
                    charcode = parse_unicode_character(ctx, output);
                    if (charcode == -1) return NULL;
                    char_to_use = (char)(0xFF & charcode);
                    break;

                default:
                    ctx->error = cpl_sprintf("Invalid escape sequence '\\%c'",
                                             *ctx->current);
                    return NULL;
            }
        }

        *(output++) = char_to_use;
    }

    *(output++) = '\0'; /* Make sure the resultant string is null terminated. */

    /* Move the working buffer position so that multiple calls to
       parse_string_characters do not clobber each other's results. */
    ctx->buffer = output;

    return result;
}

/*
 * Parses a JSON string from the current parse position.
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return The null terminated string value or NULL if there were parse errors.
 */
static const char * parse_string(json_parser_context * ctx)
{
    assert(ctx != NULL);

    /* Move past the first '"' character. */
    assert(*ctx->current == '"');
    ++ctx->current;

    const char * string = parse_string_characters(ctx);
    if (string == NULL) return NULL;

    /* Make sure we reached the end of the string. */
    if (*ctx->current == '"')
    {
        ++ctx->current;
    }
    else
    {
        /* NOTE: the parse_string_characters function will always get to the
           end of the input if there is no '"' character to terminate the
           string. Therefore there is no point in trying to print an error
           message that uses ctx->current. */
        ctx->error = cpl_sprintf("Expected '\"' character to end a string");
        return NULL;
    }

    return string;
}

/*
 * Parses a JSON string from the current parse position and returns a JSON
 * string node object.
 *
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return A new string node or NULL if there was a parsing error.
 */
static er_json_node * parse_string_as_object(json_parser_context * ctx)
{
    assert(ctx != NULL);
    const char * location = ctx->current;  /* Start location of string. */

    const char * string = parse_string(ctx);
    if (string == NULL) return NULL;

    return er_json_node_new_string(string, location);
}

/* Need to forward declare this function which is called recursively in the
   parse_array and parse_object functions.
 */
static er_json_node * parse_value(json_parser_context * ctx);

/*
 * Parses a JSON array from the current parse position.
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return A new array node object or NULL if there was a parsing error.
 */
static er_json_node * parse_array(json_parser_context * ctx)
{
    assert(ctx != NULL);

    er_json_node * array = er_json_node_new_array(ctx->current);

    /* Move past the first '[' character. */
    assert(*ctx->current == '[');
    ++ctx->current;

    /* Check for white space in case this is an empty array. */
    parse_whitespace(ctx);

    if (*ctx->current != ']')  /* Check if the array was empty. */
    {
        do
        {
            if (*ctx->current == ',') ++ctx->current;

            parse_whitespace(ctx);  /* White space before a JSON value. */
            if (*ctx->current == '\0')
            {
                ctx->error = cpl_sprintf("Expected the start of a new value");
                er_json_node_delete(array);
                return NULL;
            }
            er_json_node * item = parse_value(ctx);
            if (item == NULL)
            {
                er_json_node_delete(array);
                return NULL;
            }
            er_json_node_array_append(array, item);
            parse_whitespace(ctx);  /* White space after a JSON value. */
        }
        while (*ctx->current == ',');
    }

    /* Make sure we reached the end of the array. */
    if (*ctx->current == ']')
    {
        ++ctx->current;
    }
    else
    {
        if (*ctx->current == '\0')
        {
            ctx->error = cpl_sprintf("Expected ']' character to end the array");
        }
        else
        {
            ctx->error = cpl_sprintf("Expected ']' character to end the array,"
                                     " but found '%c' instead,", *ctx->current);
        }
        er_json_node_delete(array);
        return NULL;
    }

    return array;
}

/*
 * Parses a JSON object from the current parse position.
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return A new object node or NULL if there was a parsing error.
 */
static er_json_node * parse_object(json_parser_context * ctx)
{
    assert(ctx != NULL);

    er_json_node * object = er_json_node_new_object(ctx->current);

    /* Move past the first '{' character. */
    assert(*ctx->current == '{');
    ++ctx->current;

    /* Check for white space in case this is an empty object. */
    parse_whitespace(ctx);

    if (*ctx->current != '}')  /* Check if the JSON object was empty. */
    {
        do
        {
            if (*ctx->current == ',') ++ctx->current;

            /* First parse the key string. */
            parse_whitespace(ctx);  /* White space before the key string. */
            if (*ctx->current != '"')
            {
                if (*ctx->current == '\0')
                {
                    ctx->error = cpl_sprintf("Expected a key string");
                }
                else if (strlen(ctx->current) > 20)
                {
                    ctx->error = cpl_sprintf("Expected a key string, but found"
                                             " '%.20s...' instead,",
                                             ctx->current);
                }
                else
                {
                    ctx->error = cpl_sprintf("Expected a key string, but found"
                                             " '%.20s' instead,", ctx->current);
                }
                er_json_node_delete(object);
                return NULL;
            }
            const char * key = parse_string(ctx);
            if (key == NULL)
            {
                er_json_node_delete(object);
                return NULL;
            }
            parse_whitespace(ctx);  /* White space after the key string. */

            /* The next character must be the separator ':'. */
            if (*ctx->current != ':')
            {
                if (strlen(ctx->current) > 20)
                {
                    ctx->error = cpl_sprintf("Expected ':' but found '%.20s...'"
                                             " instead,", ctx->current);
                }
                else
                {
                    ctx->error = cpl_sprintf("Expected ':' but found '%.20s'"
                                             " instead,", ctx->current);
                }
                er_json_node_delete(object);
                return NULL;
            }
            ++ctx->current;  /* Move past the ':' character. */

            /* Now the value part should follow. */
            parse_whitespace(ctx);  /* White space before the JSON value. */
            if (*ctx->current == '\0')
            {
                ctx->error = cpl_sprintf("Expected the start of a new value");
                er_json_node_delete(object);
                return NULL;
            }
            er_json_node * item = parse_value(ctx);
            if (item == NULL)
            {
                er_json_node_delete(object);
                return NULL;
            }
            er_json_node_object_insert(object, key, item);
            parse_whitespace(ctx);  /* White space after the JSON value. */
        }
        while (*ctx->current == ',');

        /* Check for left over white space before the last '}'. */
        parse_whitespace(ctx);
    }

    /* Make sure we reached the end of the JSON object. */
    if (*ctx->current == '}')
    {
        ++ctx->current;
    }
    else
    {
        if (*ctx->current == '\0')
        {
            ctx->error = cpl_sprintf("Expected '}' character to end a object");
        }
        else
        {
            ctx->error = cpl_sprintf("Expected '}' character to end a object,"
                                     " but found '%c' instead,", *ctx->current);
        }
        er_json_node_delete(object);
        return NULL;
    }

    return object;
}

/*
 * Parses a JSON value from the current parse position. This can be either a
 * null, boolean, number, string, array or object.
 *
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return A new JSON node or NULL if there was a parsing error.
 */
static er_json_node * parse_value(json_parser_context * ctx)
{
    assert(ctx != NULL);
    assert(*ctx->current != '\0');

    /* If the current character is a digit then this should be a number. */
    if (isdigit(*ctx->current))
    {
        return parse_number(ctx);
    }

    /* Identify the type of value we are dealing with by the current character.
       We then invoke the appropriate parsing function. */
    switch (*ctx->current)
    {
        case 'n':
            /* If the first character is a 'n' we might be dealing with either
               a 'null' terminal or 'nan' number terminal. We need to check what
               the character after the 'n' is to decide which one we are dealing
               with. */
            if (*(ctx->current+1) == 'u')
            {
                return parse_null(ctx);
            }
            else
            {
                return parse_number(ctx);
            }
        case 't': case 'f':
            return parse_boolean(ctx);
        case '-': case 'i': case 'I': case 'N':
            return parse_number(ctx);
        case '"':
            return parse_string_as_object(ctx);
        case '[':
            return parse_array(ctx);
        case '{':
            return parse_object(ctx);
        default:
            if (strlen(ctx->current) > 20)
            {
                ctx->error = cpl_sprintf("Expected the start of a JSON value,"
                                         " but found '%.20s...' instead,",
                                         ctx->current);
            }
            else
            {
                ctx->error = cpl_sprintf("Expected the start of a JSON value,"
                                         " but found '%.20s' instead,",
                                         ctx->current);
            }
            return NULL;
    }
}

/*
 * Implements a recursive descent parser to parse the JSON text pointed to by
 * the current parse position. Any errors will be marked by setting the error
 * message in the ctx->error pointer.
 *
 * @param ctx  Pointer to the parsing context data for the current parse.
 * @return CPL_TRUE if the parse was successful or CPL_FALSE otherwise.
 */
static cpl_boolean parse_json(json_parser_context * ctx)
{
    parse_whitespace(ctx);
    if (*ctx->current == '\0')
    {
        /* Special case, where the text only contains white space. */
        ctx->tree = er_json_node_new_null(ctx->json);
        return CPL_TRUE;
    }
    ctx->tree = parse_value(ctx);
    if (ctx->tree == NULL) return CPL_FALSE;
    parse_whitespace(ctx);
    if (*ctx->current != '\0')
    {
        ctx->error = cpl_sprintf("Unexpected character '%c' after a JSON value"
                                 " node", *ctx->current);
        return CPL_FALSE;
    }
    return CPL_TRUE;
}

/**
 * @brief Finds line and column numbers for a character position within a text.
 *
 * @param[in]  json      Original JSON text.
 * @param[in]  location  A character position in the text, i.e. within @p json.
 * @param[out] line      Pointer to an integer that will contain the computed
 *                       line number.
 * @param[out] column    Pointer to an integer that will contain the computed
 *                       column number.
 *
 * @return @c CPL_ERROR_NONE on success, or an appropriate error code if the
 *         input pointers were not valid.
 *
 * This function will compute the line and column number corresponding to the
 * given character location within the JSON text.
 */
cpl_error_code er_json_find_line_column(const char * json,
                                        const char * location,
                                        int * line, int * column)
{
    cpl_error_ensure(json != NULL && location != NULL &&
                     line != NULL && column != NULL,
                     CPL_ERROR_NULL_INPUT, return CPL_ERROR_NULL_INPUT,
                     "Used NULL pointer.");

    *line = 1;
    *column = 1;

    /* The following sanity checks makes sure that the location pointer is
       within the JSON text. */
    assert(location >= json);
    assert(location <= json + strlen(json) + 1);

    /* Find the line, column position of the error. */
    const char * c = NULL;
    for (c = json; c != location; ++c)
    {
        if (*c == '\n')
        {
            ++(*line);
            *column = 1;
        }
        else
        {
            ++(*column);
        }
    }

    return CPL_ERROR_NONE;
}

/**
 * @brief Parses a JSON text and produces a corresponding parse tree.
 *
 * @param json  The JSON null terminated text string to parse.
 *
 * @return A JSON node object to the root of the parse tree, or @c NULL if there
 *         was a parse error. In that case a CPL error code with an appropriate
 *         message is set.
 *
 * This function will parse the input JSON text string and produce a parse tree.
 * The parse tree is a tree of er_json_node objects representing the parsed
 * information. The tree can be navigated to extract the converted string,
 * boolean and floating point number data.
 */
er_json_node * er_json_parse(const char * json)
{
    cpl_error_ensure(json != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON text is NULL.");

    /* Initialise the parser context, allocating a working buffer.
       NOTE: The extra 8 bytes guarantee there is space for a NULL terminating
       character or room for the parse_unicode_character to copy too in case we
       are near the end of the JSON string.  */
    char * buffer = cpl_malloc(strlen(json) + 8);
    json_parser_context ctx;
    ctx.json = json;
    ctx.buffer = buffer;
    ctx.current = json;
    ctx.error = NULL;
    ctx.tree = NULL;

    /* Perform the parse. */
    if (! parse_json(&ctx))
    {
        /* If an error occurred then delete the parse tree, generate the
           appropriate error code/message and return NULL. */
        er_json_node_delete(ctx.tree);
        ctx.tree = NULL;
        int line, col;
        er_json_find_line_column(ctx.json, ctx.current, &line, &col);
        assert(ctx.error != NULL);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "%s at line %d column %d.", ctx.error, line, col);
    }

    /* Cleanup the working buffer and any error message. NOTE: cpl_free already
       checks if the pointer is NULL so its safe to have ctx.error == NULL. */
    cpl_free(ctx.error);
    cpl_free(buffer);

    return ctx.tree;
}

/**
 * @brief Converts a parsed JSON array node into a er_stringarray_t object.
 *
 * @param parsetree  The parse tree for the JSON text, as produced by the
 *                   er_json_parse() function. This can be a subnode of the
 *                   result produced by er_json_parse().
 * @param json       The original JSON text from which the @p parsetree node was
 *                   produced.
 *
 * @return @c CPL_ERROR_NONE on success or an appropriate error code otherwise.
 *
 * Converts a JSON array node that is returned by the er_json_parse() function,
 * into a er_stringarray_t object, which contains a list of strings. One of
 * the subnodes of the resultant parse tree from er_json_parse() can be used
 * instead, if we only want to convert a child element from the JSON text.
 */
er_stringarray_t * er_json_to_string_array(const er_json_node * parsetree,
                                           const char * json)
{
    cpl_error_ensure(parsetree != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Input JSON parse tree is NULL.");
    cpl_error_ensure(json != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON text is NULL.");

    /* Make sure the type of the top level node is an array. */
    if (er_json_node_type(parsetree) != JSON_ARRAY)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an array of strings at line %d column"
                              " %d.", line, col);
        return NULL;
    }

    er_stringarray_t * array = er_stringarray_new();

    /* Add all the items in the parsed array to the er_stringarray_t object,
       making sure that what we parsed for each item is actually a string. */
    er_json_array_iterator iter = NULL;
    for (iter =  er_json_node_array_begin(parsetree);
         iter != er_json_node_array_end(parsetree);
         iter = er_json_node_array_next(parsetree, iter))
    {
        const er_json_node * item = er_json_node_array_get(parsetree, iter);
        if (er_json_node_type(item) != JSON_STRING)
        {
            int line, col;
            er_json_find_line_column(json, er_json_node_location(item),
                                     &line, &col);
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Expected a string item in the array at"
                                  " line %d column %d.", line, col);
            er_stringarray_delete(array);
            return NULL;
        }
        er_stringarray_append(array, er_json_node_get_string(item));
    }

    return array;
}

#ifdef ENABLE_PYTHON_RECIPES

/*
 * Attempts to fetch a child value node from a JSON object node.
 *
 * @param object  The JSON object node to query.
 * @param key     The name of the child node to find.
 * @param json    The original JSON text string from which the object node was
 *                parsed.
 *
 * @return  The JSON value node corresponding to the given key. @c NULL is
 *      returned if the child node could not be found.
 */
static const er_json_node * get_value_from_object(const er_json_node * object,
                                                  const char * key,
                                                  const char * json)
{
    assert(object != NULL);
    assert(key != NULL);
    assert(json != NULL);

    const er_json_node * value = er_json_node_object_get(object, key);
    if (value == NULL)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(object),
                                    &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Failed to find key '%s' in object at line %d"
                              " column %d.", key, line, col);
    }
    return value;
}

/*
 * Attempts to extract a string value corresponding to a given key from a JSON
 * object node.
 *
 * @param object  The JSON object node to query.
 * @param key     The name of the child node to find.
 * @param json    The original JSON text string from which the object node was
 *                parsed.
 * @param must_exist  Flag indicating if an error should be set if the key is
 *                    missing.
 *
 * @return  The string value corresponding to the given key. @c NULL is returned
 *      if the given key does not exist in the JSON object or the value is not a
 *      string.
 */
static const char * get_string_from_object(const er_json_node * object,
                                           const char * key, const char * json,
                                           cpl_boolean must_exist)
{
    assert(object != NULL);
    assert(key != NULL);
    assert(json != NULL);

    const er_json_node * value = er_json_node_object_get(object, key);
    if (value == NULL)
    {
        if (must_exist)
        {
            int line, col;
            er_json_find_line_column(json, er_json_node_location(object),
                                     &line, &col);
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Failed to find key '%s' in object at line %d"
                                  " column %d.", key, line, col);
        }
        return NULL;
    }
    if (er_json_node_type(value) != JSON_STRING)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(value),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected a string for key '%s' at line %d"
                              " column %d.", key, line, col);
        return NULL;
    }
    return er_json_node_get_string(value);
}

/*
 * Attempts to extract a number value as a floating-point double, corresponding
 * to a given key from a JSON object node.
 *
 * @param object  The JSON object node to query.
 * @param key     The name of the child node to find.
 * @param json    The original JSON text string from which the object node was
 *                parsed.
 * @param must_exist  Flag indicating if an error should be set if the key is
 *                    missing.
 *
 * @return  The double value corresponding to the given key. @c NULL is returned
 *      if the given key does not exist in the JSON object or the value is not a
 *      number.
 */
static double get_number_from_object(const er_json_node * object,
                                     const char * key, const char * json,
                                     cpl_boolean must_exist)
{
    assert(object != NULL);
    assert(key != NULL);
    assert(json != NULL);

    const er_json_node * value = er_json_node_object_get(object, key);
    if (value == NULL)
    {
        if (must_exist)
        {
            int line, col;
            er_json_find_line_column(json, er_json_node_location(object),
                                     &line, &col);
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Failed to find key '%s' in object at line %d"
                                  " column %d.", key, line, col);
        }
        return 0;
    }
    if (er_json_node_type(value) != JSON_NUMBER)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(value),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected a number for key '%s' at line %d"
                              " column %d.", key, line, col);
        return 0;
    }
    return er_json_node_get_number(value);
}

/*
 * Attempts to extract an unsigned long integer value corresponding to a given
 * key from a JSON object node.
 *
 * @param object  The JSON object node to query.
 * @param key     The name of the child node to find.
 * @param json    The original JSON text string from which the object node was
 *                parsed.
 * @param must_exist  Flag indicating if an error should be set if the key is
 *                    missing.
 *
 * @return  The integer value corresponding to the given key. @c NULL is
 *      returned if the given key does not exist in the JSON object or the value
 *      is not a number.
 */
static unsigned long get_ulong_from_object(const er_json_node * object,
                                           const char * key, const char * json,
                                           cpl_boolean must_exist)
{
    return (unsigned long) get_number_from_object(object, key, json,
                                                  must_exist);
}

/*
 * Attempts to extract a cpl_size type integer value corresponding to a given
 * key from a JSON object node.
 *
 * @param object  The JSON object node to query.
 * @param key     The name of the child node to find.
 * @param json    The original JSON text string from which the object node was
 *                parsed.
 *
 * @return  The integer value corresponding to the given key as cpl_size.
 *      @c NULL is returned if the given key does not exist in the JSON object
 *      or the value is not a number.
 */
static cpl_size get_cpl_size_from_object(const er_json_node * object,
                                         const char * key, const char * json)
{
    return (cpl_size) get_number_from_object(object, key, json, CPL_TRUE);
}

/*
 * Finds and returns the JSON array node called 'choices' from a JSON object.
 *
 * @param object  The JSON object node to query.
 * @param json    The original JSON text string from which the object node was
 *                parsed.
 *
 * @return  The JSON array node called 'choices'. @c NULL is returned if the
 *      node could not be found or does not correspond to an array type.
 */
static const er_json_node * get_choices_array(const er_json_node * object,
                                              const char * json)
{
    assert(object != NULL);
    assert(json != NULL);

    /* Fetch the "choices" JSON node from the given object parse tree and make
       sure it is an array type. */
    const char * key = "choices";
    const er_json_node * array_node = get_value_from_object(object, key, json);
    if (array_node == NULL) return NULL;
    if (er_json_node_type(array_node) != JSON_ARRAY)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(array_node),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an array for key '%s' at line %d"
                              " column %d.", key, line, col);
        return NULL;
    }
    return array_node;
}

/*
 * Returns true if the representation of the JSON object is compatible with a
 * pure integer. i.e. this is to check if a number was represented as an integer
 * or floating point.
 */
static cpl_boolean node_is_integer(const er_json_node * object)
{
    const char * text = er_json_node_location(object);
    cpl_boolean result = CPL_TRUE;
    const char * c = NULL;
    for (c = text; *c != '\0' && *c != ']' && *c != '}' && *c != ','; ++c)
    {
        if (isspace(*c) || isdigit(*c) || *c == '-' || *c == '+') continue;
        result = CPL_FALSE;
        break;
    }
    return result;
}

/*
 * Extracts the values from the 'choices' array found in a JSON object as
 * integers. These are written to the output pointers passed to this function.
 *
 * @param[in] object  The JSON object node to extract the choices from.
 * @param[in] json    The original JSON text string from which the object node
 *                    was parsed.
 * @param[in] default_value  A default value that must exist in the choices
 *                           array. An error is returned if this value could not
 *                           be found.
 * @param[out] choices_count  A pointer that will store the number of choices
 *                            that was written to the @c choices array.
 * @param[out] choices  A pointer that will store the new array of choices.
 * @param[out] mismatch  A pointer where the node that has a type mismatch is
 *                       recorded if any. This is only filled if this function
 *                       returns @c CPL_ERROR_TYPE_MISMATCH.
 *
 * @return  CPL_ERROR_NONE is returned on success and an appropriate error code
 *      otherwise.
 *
 * @note The caller must call @c cpl_free on the @p choices array if this
 *      function returned a @c CPL_ERROR_NONE error code.
 */
static cpl_error_code get_ints_from_choices(const er_json_node * object,
                                            const char * json,
                                            int default_value,
                                            int * choices_count,
                                            int ** choices,
                                            const er_json_node ** mismatch)
{
    assert(object != NULL);
    assert(json != NULL);
    assert(choices_count != NULL);
    assert(choices != NULL);
    assert(mismatch != NULL);

    const er_json_node * array_node = get_choices_array(object, json);
    if (array_node == NULL) return cpl_error_get_code();

    /* Check every element in the array to make sure that it is a integer.
       We also check that at least one of the choices corresponds to the default
       value. */
    cpl_boolean default_not_found = CPL_TRUE;
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(array_node);
         iter != er_json_node_array_end(array_node);
         iter = er_json_node_array_next(array_node, iter))
    {
        const er_json_node * item = er_json_node_array_get(array_node, iter);
        if (er_json_node_type(item) != JSON_NUMBER)
        {
            /* Mark the JSON node that has a type mismatch, this should be later
               handled by the called. */
            *mismatch = item;
            return CPL_ERROR_TYPE_MISMATCH;
        }
        double num = er_json_node_get_number(item);
        if (! node_is_integer(item))  /* Check if num is an integer. */
        {
            *mismatch = item;
            return CPL_ERROR_TYPE_MISMATCH;
        }
        int inum = (int)num;
        if (inum == default_value)
        {
            default_not_found = CPL_FALSE;
        }
    }
    if (default_not_found)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(array_node),
                                 &line, &col);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "The default value %d was not found in"
                                     " the array of choices at line %d column"
                                     " %d.", default_value, line, col);
    }

    /* We now allocate memory and construct the output array. */
    cpl_size size = er_json_node_array_size(array_node);
    *choices = cpl_malloc(size * sizeof(int));
    *choices_count = (int)size;
    int n;
    for (iter = er_json_node_array_begin(array_node), n = 0;
         iter != er_json_node_array_end(array_node);
         iter = er_json_node_array_next(array_node, iter), ++n)
    {
        const er_json_node * item = er_json_node_array_get(array_node, iter);
        (*choices)[n] = (int)er_json_node_get_number(item);
    }
    return CPL_ERROR_NONE;
}

/*
 * This function behaves the same way as get_ints_from_choices, but handles
 * choices that are of floating point type.
 * See get_ints_from_choices for details.
 */
static cpl_error_code get_doubles_from_choices(const er_json_node * object,
                                               const char * json,
                                               double default_value,
                                               int * choices_count,
                                               double ** choices,
                                               const er_json_node ** mismatch)
{
    assert(object != NULL);
    assert(json != NULL);
    assert(choices_count != NULL);
    assert(choices != NULL);
    assert(mismatch != NULL);

    const er_json_node * array_node = get_choices_array(object, json);
    if (array_node == NULL) return cpl_error_get_code();

    /* Check every element in the array to make sure that it is a double
       floating point. We also check that at least one of the choices
       corresponds to the default value. */
    cpl_boolean default_not_found = CPL_TRUE;
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(array_node);
         iter != er_json_node_array_end(array_node);
         iter = er_json_node_array_next(array_node, iter))
    {
        const er_json_node * item = er_json_node_array_get(array_node, iter);
        if (er_json_node_type(item) != JSON_NUMBER)
        {
            /* Mark the JSON node that has a type mismatch, this should be later
               handled by the called. */
            *mismatch = item;
            return CPL_ERROR_TYPE_MISMATCH;
        }
        if (er_json_node_get_number(item) == default_value)
        {
            default_not_found = CPL_FALSE;
        }
    }
    if (default_not_found)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(array_node),
                                 &line, &col);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "The default value %.17g was not found in"
                                     " the array of choices at line %d column"
                                     " %d.", default_value, line, col);
    }

    /* We now allocate memory and construct the output array. */
    cpl_size size = er_json_node_array_size(array_node);
    *choices = cpl_malloc(size * sizeof(double));
    *choices_count = (int)size;
    int n;
    for (iter = er_json_node_array_begin(array_node), n = 0;
         iter != er_json_node_array_end(array_node);
         iter = er_json_node_array_next(array_node, iter), ++n)
    {
        const er_json_node * item = er_json_node_array_get(array_node, iter);
        (*choices)[n] = er_json_node_get_number(item);
    }
    return CPL_ERROR_NONE;
}

/*
 * This function behaves the same way as get_ints_from_choices, but handles
 * choices that are of string type.
 * See get_ints_from_choices for details.
 */
static cpl_error_code get_strings_from_choices(const er_json_node * object,
                                               const char * json,
                                               const char * default_value,
                                               int * choices_count,
                                               const char *** choices,
                                               const er_json_node ** mismatch)
{
    assert(object != NULL);
    assert(json != NULL);
    assert(default_value != NULL);
    assert(choices_count != NULL);
    assert(choices != NULL);
    assert(mismatch != NULL);

    const er_json_node * array_node = get_choices_array(object, json);
    if (array_node == NULL) return cpl_error_get_code();

    /* Check every element in the array to make sure that it is a string.
       We also check that at least one of the choices corresponds to the default
       value. */
    cpl_boolean default_not_found = CPL_TRUE;
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(array_node);
         iter != er_json_node_array_end(array_node);
         iter = er_json_node_array_next(array_node, iter))
    {
        const er_json_node * item = er_json_node_array_get(array_node, iter);
        if (er_json_node_type(item) != JSON_STRING)
        {
            /* Mark the JSON node that has a type mismatch, this should be later
               handled by the called. */
            *mismatch = item;
            return CPL_ERROR_TYPE_MISMATCH;
        }
        if (strcmp(er_json_node_get_string(item), default_value) == 0)
        {
            default_not_found = CPL_FALSE;
        }
    }
    if (default_not_found)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(array_node),
                                 &line, &col);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "The default value '%s' was not found in"
                                     " the array of choices at line %d column"
                                     " %d.", default_value, line, col);
    }

    /* We now allocate memory and construct the output array. */
    cpl_size size = er_json_node_array_size(array_node);
    *choices = cpl_malloc(size * sizeof(char *));
    *choices_count = (int)size;
    int n;
    for (iter = er_json_node_array_begin(array_node), n = 0;
         iter != er_json_node_array_end(array_node);
         iter = er_json_node_array_next(array_node, iter), ++n)
    {
        const er_json_node * item = er_json_node_array_get(array_node, iter);
        (*choices)[n] = er_json_node_get_string(item);
    }
    return CPL_ERROR_NONE;
}

/*
 * Enables or disables a CPL parameter based on the boolean value of a given
 * JSON sub node in the parse tree. This is done for a particular CPL parameter
 * mode i.e. CLI, ENV or CFG.
 *
 * @param[out] param  The CPL parameter object to update.
 * @param[in]  parsetree  The JSON parse tree to query. This must point to a
 *                        JSON object type node.
 * @param[in]  json   The original JSON text string from which the parse tree
 *                    was generated.
 * @param[in]  key    The name of the boolean node to use to decide if the CPL
 *                    parameter should be enabled or not.
 * @param[in]  mode   The CPL parameter mode to update.
 *
 * @return  @c CPL_ERROR_NONE on success and an appropriate error code
 *      otherwise.
 */
static cpl_error_code set_param_enable_flag(cpl_parameter * param,
                                            const er_json_node * parsetree,
                                            const char * json,
                                            const char * key,
                                            cpl_parameter_mode mode)
{
    assert(param != NULL);
    assert(parsetree != NULL);
    assert(json != NULL);
    assert(key != NULL);

    const er_json_node * node = er_json_node_object_get(parsetree, key);
    if (node != NULL)
    {
        if (er_json_node_type(node) != JSON_BOOL)
        {
            int line, col;
            er_json_find_line_column(json, er_json_node_location(node),
                                     &line, &col);
            return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                         "Expected a boolean at line %d column"
                                         " %d.", line, col);
        }
        if (er_json_node_get_bool(node))
        {
            cpl_parameter_enable(param, mode);
        }
        else
        {
            cpl_parameter_disable(param, mode);
        }
    }
    return CPL_ERROR_NONE;
}

/*
 * Sets the CPL parameter's alias string for a particular mode, based on the
 * value of a child node in a JSON parse tree.
 *
 * @param[out] param  The CPL parameter object to update.
 * @param[in]  parsetree  The JSON parse tree to query. This must point to a
 *                        JSON object type node.
 * @param[in]  json   The original JSON text string from which the parse tree
 *                    was generated.
 * @param[in]  key    The name of a string node who's value will be used for
 *                    the alias.
 * @param[in]  mode   The CPL parameter mode to update.
 *
 * @return  @c CPL_ERROR_NONE on success and an appropriate error code
 *      otherwise.
 */
static cpl_error_code set_param_alias(cpl_parameter * param,
                                      const er_json_node * parsetree,
                                      const char * json,
                                      const char * key,
                                      cpl_parameter_mode mode)
{
    assert(param != NULL);
    assert(parsetree != NULL);
    assert(json != NULL);
    assert(key != NULL);

    const er_json_node * node = er_json_node_object_get(parsetree, key);
    if (node != NULL)
    {
        if (er_json_node_type(node) == JSON_STRING)
        {
            cpl_parameter_set_alias(param, mode, er_json_node_get_string(node));
        }
        else if (er_json_node_type(node) == JSON_NULL)
        {
            cpl_parameter_set_alias(param, mode, NULL);
        }
        else
        {
            int line, col;
            er_json_find_line_column(json, er_json_node_location(node),
                                     &line, &col);
            return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                         "Expected a string or null at line %d"
                                         " column %d.", line, col);
        }
    }
    return CPL_ERROR_NONE;
}

/*
 * Converts a JSON parse tree to a CPL parameter object.
 *
 * @param  parsetree  The JSON parse tree or sub tree returned by
 *                    er_json_parse().
 * @param  json       The original JSON text string from which the parse tree
 *                    was generated.
 *
 * @return A new CPL parameter object corresponding to the parsed JSON or
 *      @c NULL if an error occurred. Errors indicate that the parsed JSON does
 *      not correspond to a proper CPL parameter object.
 */
static cpl_parameter * json_to_parameter(const er_json_node * parsetree,
                                         const char * json)
{
    int line, col;

    if (er_json_node_type(parsetree) != JSON_OBJECT)
    {
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an object at line %d column %d.",
                              line, col);
        return NULL;
    }

    /* Find the common attributes for the parameter object in the JSON parse
       tree. For the description and context we use empty strings as default
       values if they were not found in the parse tree. */
    const char * name = get_string_from_object(parsetree, "name", json,
                                               CPL_TRUE);
    const char * class_name = get_string_from_object(parsetree, "class", json,
                                                     CPL_TRUE);
    const char * desc = get_string_from_object(parsetree, "description", json,
                                               CPL_FALSE);
    if (desc == NULL) desc = "";
    const char * context = get_string_from_object(parsetree, "context", json,
                                                  CPL_FALSE);
    if (context == NULL) context = "";

    const er_json_node * node = get_value_from_object(parsetree, "default",
                                                      json);
    if (node == NULL) return NULL;
    if (cpl_error_get_code() != CPL_ERROR_NONE) return NULL;

    /* Construct the CPL parameter object. We have to perform a custom
       invocation of one of cpl_parameter_new_value, cpl_parameter_new_range or
       cpl_parameter_new_enum, depending on the class indicated in the JSON
       parse tree (i.e. the value of class_name). */
    const er_json_node * mismatch_node = NULL;
    cpl_parameter * param = NULL;
    if (strcmp(class_name, "value") == 0)
    {
        double num = NAN;
        /* We have to select the appropriate CPL type for the parameter object
           based on the JSON type of the default value's JSON node. */
        switch (er_json_node_type(node))
        {
            case JSON_BOOL:
                param = cpl_parameter_new_value(name, CPL_TYPE_BOOL,
                                                desc, context,
                                                er_json_node_get_bool(node));
                break;

            case JSON_NUMBER:
                num = er_json_node_get_number(node);
                if (node_is_integer(node))  /* Check if num is an integer. */
                {
                    int inum = (int)num;
                    param = cpl_parameter_new_value(name, CPL_TYPE_INT, desc,
                                                    context, inum);
                }
                else
                {
                    param = cpl_parameter_new_value(name, CPL_TYPE_DOUBLE, desc,
                                                    context, num);
                }
                break;

            case JSON_STRING:
                param = cpl_parameter_new_value(name, CPL_TYPE_STRING,
                                                desc, context,
                                                er_json_node_get_string(node));
                break;

            default:
                er_json_find_line_column(json, er_json_node_location(node),
                                         &line, &col);
                cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                      "Type not supported for value at line %d"
                                      " column %d.", line, col);
                return NULL;
        }
    }
    else if (strcmp(class_name, "range") == 0)
    {
        double num = NAN;
        double min = NAN;
        double max = NAN;
        /* Fetch the JSON nodes corresponding to the minimum and maximum values
           of the range. */
        const er_json_node * min_node = get_value_from_object(parsetree, "min",
                                                              json);
        const er_json_node * max_node = get_value_from_object(parsetree, "max",
                                                              json);
        if (min_node == NULL || max_node == NULL) return NULL;

        switch (er_json_node_type(node))
        {
            case JSON_NUMBER:
                if (er_json_node_type(min_node) != JSON_NUMBER)
                {
                    mismatch_node = min_node;
                    break;
                }
                if (er_json_node_type(max_node) != JSON_NUMBER)
                {
                    mismatch_node = max_node;
                    break;
                }
                num = er_json_node_get_number(node);
                min = er_json_node_get_number(min_node);
                max = er_json_node_get_number(max_node);
                if (node_is_integer(node))  /* Check if num is an integer. */
                {
                    /* Check that min and max are also integers. */
                    if (! node_is_integer(min_node))
                    {
                        mismatch_node = min_node;
                        break;
                    }
                    if (! node_is_integer(max_node))
                    {
                        mismatch_node = max_node;
                        break;
                    }

                    /* Convert to an actual int type and create the CPL
                       parameter object */
                    int inum = (int)num;
                    int imin = (int)min;
                    int imax = (int)max;
                    param = cpl_parameter_new_range(name, CPL_TYPE_INT, desc,
                                                    context, inum, imin, imax);
                }
                else
                {
                    /* The number is already a double in this case so no
                       typecasting is needed. */
                    param = cpl_parameter_new_range(name, CPL_TYPE_DOUBLE, desc,
                                                    context, num, min, max);
                }
                break;

            default:
                er_json_find_line_column(json, er_json_node_location(node),
                                         &line, &col);
                cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                      "Type not supported for range at line %d"
                                      " column %d.", line, col);
                return NULL;
        }
    }
    else if (strcmp(class_name, "enum") == 0)
    {
        int type = CPL_TYPE_INVALID;
        int choices_count = 0;
        double num = NAN;
        int default_int = 0;
        int * int_choices = NULL;
        double default_num = NAN;
        double * num_choices = NULL;
        const char * default_str = NULL;
        const char ** str_choices = NULL;

        /* Have to parse the choices list appropriately depending on the type
           of the default value (i.e. depending on whether is was an integer,
           double or string. */
        switch (er_json_node_type(node))
        {
            case JSON_NUMBER:
                num = er_json_node_get_number(node);
                if (node_is_integer(node))
                {
                    type = CPL_TYPE_INT;
                    default_int = (int)num;
                    if (get_ints_from_choices(parsetree, json, default_int,
                                              &choices_count, &int_choices,
                                              &mismatch_node)
                        != CPL_ERROR_NONE)
                    {
                        if (mismatch_node == NULL) return NULL;
                    }
                }
                else
                {
                    type = CPL_TYPE_DOUBLE;
                    default_num = num;
                    if (get_doubles_from_choices(parsetree, json, default_num,
                                                 &choices_count, &num_choices,
                                                 &mismatch_node)
                        != CPL_ERROR_NONE)
                    {
                        if (mismatch_node == NULL) return NULL;
                    }
                }
                break;

            case JSON_STRING:
                type = CPL_TYPE_STRING;
                default_str = er_json_node_get_string(node);
                if (get_strings_from_choices(parsetree, json, default_str,
                                             &choices_count, &str_choices,
                                             &mismatch_node)
                    != CPL_ERROR_NONE)
                {
                    if (mismatch_node == NULL) return NULL;
                }
                break;

            default:
                er_json_find_line_column(json, er_json_node_location(node),
                                         &line, &col);
                cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                      "Type not supported for value at line %d"
                                      " column %d.", line, col);
                return NULL;
        }

        if (mismatch_node == NULL)
        {
#ifdef HAVE_LIBAVCALL
            /* Prepare the dynamic argument list and make a call to the
               cpl_parameter_new_enum function to construct the CPL object. */
            av_alist alist;
            av_start_ptr(alist, &cpl_parameter_new_enum,
                            cpl_parameter *, &param);
            av_ptr(alist, const char *, name);
            av_int(alist, type);
            av_ptr(alist, const char *, desc);
            av_ptr(alist, const char *, context);
            switch (type)
            {
                case CPL_TYPE_INT:
                    av_int(alist, default_int);
                    break;
                case CPL_TYPE_DOUBLE:
                    av_double(alist, default_num);
                    break;
                case CPL_TYPE_STRING:
                    av_ptr(alist, const char *, default_str);
                    break;
                default:
                    assert(type == CPL_TYPE_INT || type == CPL_TYPE_DOUBLE
                           || type == CPL_TYPE_STRING);
                    break;
            }
            av_int(alist, choices_count);
            int n = 0;
            switch (type)
            {
                case CPL_TYPE_INT:
                    for (n = 0; n < choices_count; ++n)
                    {
                        av_int(alist, int_choices[n]);
                    }
                    break;
                case CPL_TYPE_DOUBLE:
                    for (n = 0; n < choices_count; ++n)
                    {
                        av_double(alist, num_choices[n]);
                    }
                    break;
                case CPL_TYPE_STRING:
                    for (n = 0; n < choices_count; ++n)
                    {
                        av_ptr(alist, const char *, str_choices[n]);
                    }
                    break;
                default:
                    assert(type == CPL_TYPE_INT || type == CPL_TYPE_DOUBLE
                           || type == CPL_TYPE_STRING);
                    break;
            }
            av_call(alist);
#else /* HAVE_LIBAVCALL */
#ifdef HAVE_LIBFFI
            /* Prepare the dynamic argument list and make a call to the
               cpl_parameter_new_enum function to construct the CPL object. */
            ffi_cif cif;
            unsigned int nargs = 6 + choices_count;
            ffi_type ** argtypes = cpl_malloc(nargs*sizeof(ffi_type *));
            const void ** args = cpl_malloc(nargs*sizeof(const void *));
            argtypes[0] = &ffi_type_pointer;
            args[0] = &name;
            argtypes[1] = &ffi_type_sint;
            args[1] = &type;
            argtypes[2] = &ffi_type_pointer;
            args[2] = &desc;
            argtypes[3] = &ffi_type_pointer;
            args[3] = &context;
            switch (type)
            {
                case CPL_TYPE_INT:
                    argtypes[4] = &ffi_type_sint;
                    args[4] = &default_int;
                    break;
                case CPL_TYPE_DOUBLE:
                    argtypes[4] = &ffi_type_double;
                    args[4] = &default_num;
                    break;
                case CPL_TYPE_STRING:
                    argtypes[4] = &ffi_type_pointer;
                    args[4] = &default_str;
                    break;
                default:
                    assert(type == CPL_TYPE_INT || type == CPL_TYPE_DOUBLE
                           || type == CPL_TYPE_STRING);
                    break;
            }
            argtypes[5] = &ffi_type_sint;
            args[5] = &choices_count;
            int n = 0;
            switch (type)
            {
                case CPL_TYPE_INT:
                    for (n = 0; n < choices_count; ++n)
                    {
                        argtypes[6+n] = &ffi_type_sint;
                        args[6+n] = &int_choices[n];
                    }
                    break;
                case CPL_TYPE_DOUBLE:
                    for (n = 0; n < choices_count; ++n)
                    {
                        argtypes[6+n] = &ffi_type_double;
                        args[6+n] = &num_choices[n];
                    }
                    break;
                case CPL_TYPE_STRING:
                    for (n = 0; n < choices_count; ++n)
                    {
                        argtypes[6+n] = &ffi_type_pointer;
                        args[6+n] = &str_choices[n];
                    }
                    break;
                default:
                    assert(type == CPL_TYPE_INT || type == CPL_TYPE_DOUBLE
                           || type == CPL_TYPE_STRING);
                    break;
            }
            if (ffi_prep_cif(&cif, FFI_DEFAULT_ABI, nargs, &ffi_type_pointer,
                             argtypes)
                == FFI_OK)
            {
                ffi_call(&cif, (void*)&cpl_parameter_new_enum, &param,
                         (void**)args);
            }
            else
            {
                cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                      "Failed to make dynamic call to"
                                      " cpl_parameter_new_enum.");
            }
            cpl_free(args);
            cpl_free(argtypes);
#else /* HAVE_LIBFFI */
    /* Since we do not have either libavcall or libffi available generate a
       compiler error. */
#error Must have libavcall or libffi available.
#endif /* HAVE_LIBFFI */
#endif /* HAVE_LIBAVCALL */

            /* NOTE: cpl_free already checks if the pointer is NULL, so the
            following is safe if the array was not allocated, as long as the
            pointer was initialised to NULL when declared. */
            cpl_free(int_choices);
            cpl_free(num_choices);
            cpl_free(str_choices);
        }
    }
    else
    {
        node = er_json_node_object_get(parsetree, "class");
        er_json_find_line_column(json, er_json_node_location(node),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Invalid value for 'class' at line %d column %d."
                              " Must be one of 'value', 'range' or 'enum'.",
                              line, col);
        return NULL;
    }

    if (mismatch_node != NULL)
    {
        er_json_find_line_column(json, er_json_node_location(mismatch_node),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                "Type mismatch between value at line %d"
                                " column %d and the type used for the"
                                " parameter's default value.", line, col);
        return NULL;
    }
    if (param == NULL) return NULL;

    /* Parse the value node and make sure that its type is compatible with the
       default value's type. The parameter is updated with the parsed value if
       one is found. Otherwise it is left as the default value. */
    node = er_json_node_object_get(parsetree, "value");
    if (node != NULL)
    {
        cpl_boolean type_mismatch = CPL_TRUE;
        switch (cpl_parameter_get_type(param))
        {
            case CPL_TYPE_BOOL:
                if (er_json_node_type(node) == JSON_BOOL)
                {
                    cpl_parameter_set_bool(param, er_json_node_get_bool(node));
                    type_mismatch = CPL_FALSE;
                }
                break;

            case CPL_TYPE_INT:
                if (er_json_node_type(node) == JSON_NUMBER)
                {
                    double num = er_json_node_get_number(node);
                    if (node_is_integer(node))
                    {
                        int inum = (int)num;
                        cpl_parameter_set_int(param, inum);
                        type_mismatch = CPL_FALSE;
                    }
                }
                break;

            case CPL_TYPE_DOUBLE:
                if (er_json_node_type(node) == JSON_NUMBER)
                {
                    cpl_parameter_set_double(param,
                                             er_json_node_get_number(node));
                    type_mismatch = CPL_FALSE;
                }
                break;

            case CPL_TYPE_STRING:
                if (er_json_node_type(node) == JSON_STRING)
                {
                    cpl_parameter_set_string(param,
                                             er_json_node_get_string(node));
                    type_mismatch = CPL_FALSE;
                }
                break;

            default:
                type_mismatch = CPL_TRUE;
        }
        if (type_mismatch)
        {
            er_json_find_line_column(json, er_json_node_location(node),
                                     &line, &col);
            cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                  "Type mismatch between value at line %d"
                                  " column %d and the type used for the"
                                  " parameter's default value.", line, col);
            cpl_parameter_delete(param);
            return NULL;
        }
    }

    /* Update the flag indicating presence if it is available in the parse
       tree. */
    node = er_json_node_object_get(parsetree, "present");
    if (node != NULL)
    {
        if (er_json_node_type(node) != JSON_BOOL)
        {
            er_json_find_line_column(json, er_json_node_location(node),
                                     &line, &col);
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Expected a boolean at line %d column %d.",
                                  line, col);
            cpl_parameter_delete(param);
            return NULL;
        }
        cpl_parameter_set_default_flag(param, er_json_node_get_bool(node));
    }

    /* Update the ID if it is available in the parse tree. */
    node = er_json_node_object_get(parsetree, "id");
    if (node != NULL)
    {
        if (er_json_node_type(node) != JSON_NUMBER)
        {
            er_json_find_line_column(json, er_json_node_location(node),
                                     &line, &col);
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Expected an integer at line %d column %d.",
                                  line, col);
            cpl_parameter_delete(param);
            return NULL;
        }
        int id = (int)er_json_node_get_number(node);
        cpl_parameter_set_id(param, id);
    }

    /* Update the tag if it is available in the parse tree and not NULL. */
    node = er_json_node_object_get(parsetree, "tag");
    if (node != NULL)
    {
        if (er_json_node_type(node) == JSON_STRING)
        {
            cpl_parameter_set_tag(param, er_json_node_get_string(node));
        }
        else if (er_json_node_type(node) == JSON_NULL)
        {
            /* Nothing to do in this case. The tag is set by default and
               calling cpl_parameter_set_tag() with a NULL value for the new
               tag will generate a CPL error. */
        }
        else
        {
            er_json_find_line_column(json, er_json_node_location(node),
                                     &line, &col);
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "Expected a string or null at line %d column"
                                  " %d.", line, col);
            cpl_parameter_delete(param);
            return NULL;
        }
    }

    /* Update the enable flags if they are available in the parse tree. */
    if (set_param_enable_flag(param, parsetree, json, "cli_enabled",
                              CPL_PARAMETER_MODE_CLI)
        != CPL_ERROR_NONE)
    {
        cpl_parameter_delete(param);
        return NULL;
    }
    if (set_param_enable_flag(param, parsetree, json, "env_enabled",
                              CPL_PARAMETER_MODE_ENV)
        != CPL_ERROR_NONE)
    {
        cpl_parameter_delete(param);
        return NULL;
    }
    if (set_param_enable_flag(param, parsetree, json, "cfg_enabled",
                              CPL_PARAMETER_MODE_CFG)
        != CPL_ERROR_NONE)
    {
        cpl_parameter_delete(param);
        return NULL;
    }

    /* Update the aliases if they are available in the parse tree. */
    if (set_param_alias(param, parsetree, json, "cli_alias",
                        CPL_PARAMETER_MODE_CLI)
        != CPL_ERROR_NONE)
    {
        cpl_parameter_delete(param);
        return NULL;
    }
    if (set_param_alias(param, parsetree, json, "env_alias",
                        CPL_PARAMETER_MODE_ENV)
        != CPL_ERROR_NONE)
    {
        cpl_parameter_delete(param);
        return NULL;
    }
    if (set_param_alias(param, parsetree, json, "cfg_alias",
                        CPL_PARAMETER_MODE_CFG)
        != CPL_ERROR_NONE)
    {
        cpl_parameter_delete(param);
        return NULL;
    }

    return param;
}

/*
 * Converts a JSON parse tree to a CPL parameter list.
 *
 * @param  parsetree  The JSON parse tree or sub tree returned by
 *                    er_json_parse().
 * @param  json       The original JSON text string from which the parse tree
 *                    was generated.
 *
 * @return A new CPL parameter list corresponding to the parsed JSON or @c NULL
 *      if an error occurred. Errors indicate that the parsed JSON does not
 *      correspond to a proper CPL parameter list.
 */
static cpl_parameterlist * json_to_parameterlist(const er_json_node * parsetree,
                                                 const char * json)
{
    if (er_json_node_type(parsetree) != JSON_ARRAY)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an array at line %d column %d.",
                              line, col);
        return NULL;
    }

    cpl_parameterlist * list = cpl_parameterlist_new();

    /* Traverse every element in the JSON array and convert it into a CPL
       parameter object. Add the new parameter object to the output list and
       return the result once done. */
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(parsetree);
         iter != er_json_node_array_end(parsetree);
         iter = er_json_node_array_next(parsetree, iter))
    {
        const er_json_node * item = er_json_node_array_get(parsetree, iter);
        cpl_parameter * param = json_to_parameter(item, json);
        if (param == NULL ||
            cpl_parameterlist_append(list, param) != CPL_ERROR_NONE)
        {
            cpl_parameter_delete(param);
            cpl_parameterlist_delete(list);
            return NULL;
        }
    }

    return list;
}

/*
 * Converts a JSON parse tree to a CPL frame object.
 *
 * @param  parsetree  The JSON parse tree or sub tree returned by
 *                    er_json_parse().
 * @param  json       The original JSON text string from which the parse tree
 *                    was generated.
 *
 * @return A new CPL frame object corresponding to the parsed JSON or @c NULL
 *      if an error occurred. Errors indicate that the parsed JSON does not
 *      correspond to a proper CPL frame object.
 */
static cpl_frame * json_to_frame(const er_json_node * parsetree,
                                 const char * json)
{
    int line, col;

    if (er_json_node_type(parsetree) != JSON_OBJECT)
    {
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an object at line %d column %d.",
                              line, col);
        return NULL;
    }

    const char * filename = get_string_from_object(parsetree, "filename", json,
                                                   CPL_TRUE);
    if (filename == NULL) return NULL;
    const char * tag = get_string_from_object(parsetree, "tag", json, CPL_TRUE);
    if (tag == NULL) return NULL;
    unsigned long type = get_ulong_from_object(parsetree, "type", json,
                                               CPL_TRUE);
    unsigned long group = get_ulong_from_object(parsetree, "group", json,
                                                CPL_TRUE);
    unsigned long level = get_ulong_from_object(parsetree, "level", json,
                                                CPL_TRUE);
    if (cpl_error_get_code() != CPL_ERROR_NONE)
    {
        /* Stop here if we had any parse errors. */
        return NULL;
    }

    cpl_frame * frame = cpl_frame_new();
    cpl_frame_set_filename(frame, filename);
    cpl_frame_set_tag(frame, tag);
    cpl_frame_set_type(frame, (cpl_frame_type) type);
    cpl_frame_set_group(frame, (cpl_frame_group) group);
    cpl_frame_set_level(frame, (cpl_frame_level) level);
    return frame;
}

/*
 * Converts a JSON parse tree to a CPL frame set.
 *
 * @param  parsetree  The JSON parse tree or sub tree returned by
 *                    er_json_parse().
 * @param  json       The original JSON text string from which the parse tree
 *                    was generated.
 *
 * @return A new CPL frame set corresponding to the parsed JSON or @c NULL
 *      if an error occurred. Errors indicate that the parsed JSON does not
 *      correspond to a proper CPL frame set.
 */
static cpl_frameset * json_to_frameset(const er_json_node * parsetree,
                                       const char * json)
{
    if (er_json_node_type(parsetree) != JSON_ARRAY)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an array at line %d column %d.",
                              line, col);
        return NULL;
    }

    cpl_frameset * frameset = cpl_frameset_new();

    /* Traverse every element in the JSON array and convert it into a CPL
       frame object. Add the new frame to the output frameset and return the
       result once done. */
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(parsetree);
         iter != er_json_node_array_end(parsetree);
         iter = er_json_node_array_next(parsetree, iter))
    {
        const er_json_node * item = er_json_node_array_get(parsetree, iter);
        cpl_frame * frame = json_to_frame(item, json);
        if (frame == NULL ||
            cpl_frameset_insert(frameset, frame) != CPL_ERROR_NONE)
        {
            cpl_frame_delete(frame);
            cpl_frameset_delete(frameset);
            return NULL;
        }
    }

    return frameset;
}

/*
 * Adds all input entries found in a JSON parse tree to a CPL recipe
 * configuration object.
 *
 * @param[in]  parsetree  The JSON parse tree or sub tree returned by
 *                        er_json_parse().
 * @param[in]  json       The original JSON text string from which the parse
 *                        tree was generated.
 * @param[out] config     The recipe configuration object being updated.
 * @param[in]  tag        The primary tag for which to add the input entries in
 *                        @p config.
 *
 * @return @c CPL_ERROR_NONE on success and an appropriate error code otherwise.
 *      Errors indicate that the parsed JSON does not correspond to a proper
 *      list of input entries.
 */
static cpl_error_code add_recipeconfig_inputs(const er_json_node * parsetree,
                                              const char * json,
                                              cpl_recipeconfig * config,
                                              const char * tag)
{
    assert(parsetree != NULL);
    assert(json != NULL);
    assert(config != NULL);
    assert(tag != NULL);

    int line, col;
    if (er_json_node_type(parsetree) != JSON_ARRAY)
    {
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Expected an array at line %d column %d.",
                                     line, col);
    }

    /* Traverse the JSON array and each input entry to the configuration. */
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(parsetree);
         iter != er_json_node_array_end(parsetree);
         iter = er_json_node_array_next(parsetree, iter))
    {
        const er_json_node * item = er_json_node_array_get(parsetree, iter);

        if (er_json_node_type(item) != JSON_OBJECT)
        {
            er_json_find_line_column(json, er_json_node_location(item),
                                     &line, &col);
            return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                         "Expected an object at line %d column"
                                         " %d.", line, col);
        }

        /* Extract the input's tag string, check that it is not an empty string
           and that it does not equal the primary tag. We do these extra checks,
           since the recipeconfig object does not allow empty strings for the
           inputs or tag == input. */
        const char * input = get_string_from_object(item, "tag", json, CPL_TRUE);
        if (input == NULL) return cpl_error_get_code();
        if (strcmp(input, "") == 0)
        {
            const er_json_node * node = er_json_node_object_get(item, "tag");
            er_json_find_line_column(json, er_json_node_location(node),
                                    &line, &col);
            return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                         "Input tag is an empty string at line"
                                         " %d column %d.", line, col);
        }
        if (strcmp(input, tag) == 0)
        {
            const er_json_node * node = er_json_node_object_get(item, "tag");
            er_json_find_line_column(json, er_json_node_location(node),
                                     &line, &col);
            return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                         "Input tag at line %d column %d cannot"
                                         " be the same as the primary tag.",
                                         line, col);
        }

        /* Extract the minimum and maximum frame counts for the input. */
        cpl_size min = get_cpl_size_from_object(item, "min", json);
        cpl_size max = get_cpl_size_from_object(item, "max", json);
        if (cpl_error_get_code() != CPL_ERROR_NONE)
        {
            /* Stop here if we had any parse errors. */
            return cpl_error_get_code();
        }

        if (cpl_recipeconfig_set_input(config, tag, input,
                                       (cpl_size)min, (cpl_size)max) != 0)
        {
            return cpl_error_get_code();
        }
    }

    return CPL_ERROR_NONE;
}

/*
 * Adds all output entries found in a JSON parse tree to a CPL recipe
 * configuration object.
 *
 * @param[in]  parsetree  The JSON parse tree or sub tree returned by
 *                        er_json_parse().
 * @param[in]  json       The original JSON text string from which the parse
 *                        tree was generated.
 * @param[out] config     The recipe configuration object being updated.
 * @param[in]  tag        The primary tag for which to add the output entries in
 *                        @p config.
 *
 * @return @c CPL_ERROR_NONE on success and an appropriate error code otherwise.
 *      Errors indicate that the parsed JSON does not correspond to a proper
 *      list of output entries.
 */
static cpl_error_code add_recipeconfig_outputs(const er_json_node * parsetree,
                                               const char * json,
                                               cpl_recipeconfig * config,
                                               const char * tag)
{
    assert(parsetree != NULL);
    assert(json != NULL);
    assert(config != NULL);
    assert(tag != NULL);

    int line, col;
    if (er_json_node_type(parsetree) != JSON_ARRAY)
    {
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Expected an array at line %d column %d.",
                                     line, col);
    }

    /* Traverse the JSON array of strings . */
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(parsetree);
         iter != er_json_node_array_end(parsetree);
         iter = er_json_node_array_next(parsetree, iter))
    {
        const er_json_node * item = er_json_node_array_get(parsetree, iter);

        if (er_json_node_type(item) != JSON_STRING)
        {
            er_json_find_line_column(json, er_json_node_location(item),
                                     &line, &col);
            return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                         "Expected a string at line %d column"
                                         " %d.", line, col);
        }
        const char * output = er_json_node_get_string(item);
        /* Make sure the output tag is not an empty string since the
           cpl_recipeconfig object does not accept empty strings. */
        if (strcmp(output, "") == 0)
        {
            er_json_find_line_column(json, er_json_node_location(item),
                                     &line, &col);
            return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                         "Output tag is an empty string at line"
                                         " %d column %d.", line, col);
        }

        if (cpl_recipeconfig_set_output(config, tag, output) != 0)
        {
            return cpl_error_get_code();
        }
    }

    return CPL_ERROR_NONE;
}

/*
 * Adds a new tag to a CPL recipe configuration object based on the information
 * found in a given JSON parse tree.
 *
 * @param[in]  parsetree  The JSON parse tree or sub tree returned by
 *                        er_json_parse().
 * @param[in]  json       The original JSON text string from which the parse
 *                        tree was generated.
 * @param[out] config     The recipe configuration object being updated.
 *
 * @return @c CPL_ERROR_NONE on success and an appropriate error code otherwise.
 *      Errors indicate that the parsed JSON does not correspond to a proper
 *      tag entry.
 */
static cpl_error_code add_recipeconfig_tag(const er_json_node * parsetree,
                                           const char * json,
                                           cpl_recipeconfig * config)
{
    assert(parsetree != NULL);
    assert(json != NULL);
    assert(config != NULL);

    const er_json_node * node;
    int line, col;

    if (er_json_node_type(parsetree) != JSON_OBJECT)
    {
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Expected an object at line %d column %d.",
                                     line, col);
    }

    /* Extract the tag string and check that it is not an empty string. We do
       this check since the recipeconfig object does not allow empty strings
       for the tags. */
    const char * tag = get_string_from_object(parsetree, "tag", json, CPL_TRUE);
    if (tag == NULL) return cpl_error_get_code();
    if (strcmp(tag, "") == 0)
    {
        node = er_json_node_object_get(parsetree, "tag");
        er_json_find_line_column(json, er_json_node_location(node),
                                 &line, &col);
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Tag is an empty string at line %d column"
                                     " %d.", line, col);
    }

    /* Extract the minimum and maximum frame counts corresponding to the tag,
       then actually add the new tag entry to the configuration. */
    cpl_size min = get_cpl_size_from_object(parsetree, "min", json);
    cpl_size max = get_cpl_size_from_object(parsetree, "max", json);
    if (cpl_error_get_code() != CPL_ERROR_NONE)
    {
        /* Stop here if we had any parse errors. */
        return cpl_error_get_code();
    }

    if (cpl_recipeconfig_set_tag(config, tag, (cpl_size)min, (cpl_size)max)
        != 0)
    {
        return cpl_error_get_code();
    }

    /* Convert and add the input and output arrays if they were given. */
    node = er_json_node_object_get(parsetree, "inputs");
    if (node != NULL)
    {
        if (add_recipeconfig_inputs(node, json, config, tag) != CPL_ERROR_NONE)
        {
            return cpl_error_get_code();
        }
    }
    node = er_json_node_object_get(parsetree, "outputs");
    if (node != NULL)
    {
        if (add_recipeconfig_outputs(node, json, config, tag) != CPL_ERROR_NONE)
        {
            return cpl_error_get_code();
        }
    }

    return CPL_ERROR_NONE;
}

/*
 * Converts a JSON parse tree to a CPL recipe configuration object.
 *
 * @param  parsetree  The JSON parse tree or sub tree returned by
 *                    er_json_parse().
 * @param  json       The original JSON text string from which the parse tree
 *                    was generated.
 *
 * @return A new CPL recipe configuration object corresponding to the parsed
 *      JSON or @c NULL if an error occurred. Errors indicate that the parsed
 *      JSON does not correspond to a proper recipe configuration object.
 */
static cpl_recipeconfig * json_to_recipeconfig(const er_json_node * parsetree,
                                               const char * json)
{
    int line, col;
    if (er_json_node_type(parsetree) != JSON_ARRAY)
    {
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an array at line %d column %d.",
                              line, col);
        return NULL;
    }

    cpl_recipeconfig * config = cpl_recipeconfig_new();

    /* Traverse the JSON array and convert each tag, adding the information to
       the configuration object as we go. */
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(parsetree);
         iter != er_json_node_array_end(parsetree);
         iter = er_json_node_array_next(parsetree, iter))
    {
        const er_json_node * item = er_json_node_array_get(parsetree, iter);
        if (add_recipeconfig_tag(item, json, config) != CPL_ERROR_NONE)
        {
            cpl_recipeconfig_delete(config);
            return NULL;
        }
    }

    return config;
}

/*
 * This is a destructor method for a recipe v1 type plugin object.
 * The er_json_to_plugin function will assign the plugin's destructor method to
 * this function when it is creating a new v1 type recipe.
 */
static int recipe_deinitialize(cpl_plugin * plugin)
{
    assert(plugin != NULL);
    cpl_recipe * recipe = (cpl_recipe *) plugin;
    cpl_parameterlist_delete(recipe->parameters);
    cpl_frameset_delete(recipe->frames);
    return 0;
}

/*
 * This is a destructor method for a recipe v2 type plugin object.
 * The er_json_to_plugin function will assign the plugin's destructor method to
 * this function when it is creating a new v2 type recipe.
 */
static int recipe2_deinitialize(cpl_plugin * plugin)
{
    assert(plugin != NULL);
    cpl_recipe2 * recipe = (cpl_recipe2 *) plugin;
    cpl_parameterlist_delete(recipe->base.parameters);
    cpl_frameset_delete(recipe->base.frames);
    cpl_recipeconfig_delete(recipe->config);
    return 0;
}

/**
 * @brief Converts a JSON parse tree to a CPL plugin object.
 *
 * This function will convert a parse tree or sub tree returned by the
 * er_json_parse() function to an appropriate CPL plugin object. Depending on
 * the structure of the parsed JSON a recipe v1 or v2 object is created.
 *
 * A special parameter will also be added to the recipe's parameter list called
 * @a \__class__. This will correspond to the value of the @a class entry found
 * in the JSON parse tree. It should hold the name of the corresponding Python
 * class that implements the plugin.
 *
 * @param  parsetree  The JSON parse tree or sub tree returned by
 *                    er_json_parse().
 * @param  json       The original JSON text string from which the parse tree
 *                    was generated.
 *
 * @return  A new CPL plugin object for the parsed recipe. @c NULL is returned
 *      if an error is detected. Errors typically indicate that the parsed JSON
 *      does not actually correspond to a valid CPL plugin object.
 *
 * @note The returned object must be released by the caller with the following
 *      code when it is no longer needed:
 *      @code
 *          cpl_plugin * plugin = er_json_to_plugin(parsetree, json);
 *          plugin->deinitialize(plugin);
 *          cpl_plugin_delete(plugin);
 *      @endcode
 */
cpl_plugin * er_json_to_plugin(const er_json_node * parsetree,
                               const char * json)
{
    cpl_error_ensure(parsetree != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Input JSON parse tree is NULL.");
    cpl_error_ensure(json != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON text is NULL.");

    /* Make sure the top level node in the parse tree is an object. */
    if (er_json_node_type(parsetree) != JSON_OBJECT)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an object at line %d column %d.",
                              line, col);
        return NULL;
    }

    const char * classname = get_string_from_object(parsetree, "class", json,
                                                    CPL_TRUE);
    const char * name = get_string_from_object(parsetree, "name", json,
                                               CPL_FALSE);
    unsigned long version = get_ulong_from_object(parsetree, "version", json,
                                                  CPL_FALSE);
    const char * synopsis = get_string_from_object(parsetree, "synopsis", json,
                                                   CPL_FALSE);
    const char * desc = get_string_from_object(parsetree, "description", json,
                                               CPL_FALSE);
    const char * author = get_string_from_object(parsetree, "author", json,
                                                 CPL_FALSE);
    const char * email = get_string_from_object(parsetree, "email", json,
                                                CPL_FALSE);
    const char * cpright = get_string_from_object(parsetree, "copyright", json,
                                                  CPL_FALSE);
    if (cpl_error_get_code() != CPL_ERROR_NONE)
    {
        /* Stop here if we had any parse errors. */
        return NULL;
    }

    /* Should be ensured by error code checks above. */
    assert(classname != NULL);

    /* Set default values if they were not found in the object node. */
    if (name == NULL) name = classname;
    if (synopsis == NULL) synopsis = "unknown";
    if (desc == NULL) desc = "unknown";
    if (author == NULL) author = "unknown";
    if (email == NULL) email = "unknown";
    if (cpright == NULL) cpright = "unknown";

    const er_json_node * params = er_json_node_object_get(parsetree,
                                                          "parameters");
    if (params == NULL)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Object starting at line %d column %d does not"
                              " have a key called 'parameters'.", line, col);
        return NULL;
    }

    /* Convert the parameter list to a cpl_parameterlist object. */
    cpl_parameterlist * paramlist = json_to_parameterlist(params, json);
    if (paramlist == NULL) return NULL;

    /* Convert the frames list to a cpl_frameset object. */
    const er_json_node * frames = er_json_node_object_get(parsetree, "frames");
    cpl_frameset * frameset = NULL;
    if (frames != NULL)
    {
        frameset = json_to_frameset(frames, json);
        if (frameset == NULL)
        {
            cpl_parameterlist_delete(paramlist);
            return NULL;
        }
    }

    /* Create and add the class name parameter to the parameter list, so that we
       know which Python class to use when invoking this plugin.
       NOTE: this parameter will have to be stripped later on before giving this
       plugin object to the upper EsoRex code.
       */
    cpl_parameter * param = cpl_parameter_new_value("__class__",
                                                    CPL_TYPE_STRING,
                                                    "Python class name",
                                                    "__python__",
                                                    classname);
    if (param == NULL)
    {
        cpl_parameterlist_delete(paramlist);
        cpl_frameset_delete(frameset);
        return NULL;
    }
    if (cpl_parameterlist_append(paramlist, param) != CPL_ERROR_NONE)
    {
        cpl_parameter_delete(param);
        cpl_parameterlist_delete(paramlist);
        cpl_frameset_delete(frameset);
        return NULL;
    }

    cpl_plugin * plugin = NULL;
    unsigned long plugin_type = CPL_PLUGIN_TYPE_NONE;
    cpl_plugin_func deinitfunc = NULL;

    /* Decide what kind of recipe plugin we need to create. If the parse tree
       contains 'recipeconfig' then its a version 2 recipe, otherwise it must
       be a version 1 recipe. */
    const er_json_node * recipeconfig = er_json_node_object_get(parsetree,
                                                                "recipeconfig");
    if (recipeconfig == NULL)
    {
        cpl_recipe * recipe = cpl_calloc(1, sizeof(cpl_recipe));
        plugin = &recipe->interface;
        plugin_type = CPL_PLUGIN_TYPE_RECIPE;
        deinitfunc = recipe_deinitialize;
        recipe->parameters = paramlist;
        recipe->frames = frameset;
    }
    else
    {
        /* Parse the recipeconfig object first before allocating memory for the
           plugin object. In case there are parsing errors, will have less to
           cleanup. */
        cpl_recipeconfig * config = json_to_recipeconfig(recipeconfig, json);
        if (config == NULL)
        {
            cpl_parameterlist_delete(paramlist);
            /* The following is safe because NULL is implicitly checked for. */
            cpl_frameset_delete(frameset);
            return NULL;
        }

        cpl_recipe2 * recipe = cpl_calloc(1, sizeof(cpl_recipe2));
        plugin = &recipe->base.interface;
        plugin_type = CPL_PLUGIN_TYPE_RECIPE_V2;
        deinitfunc = recipe2_deinitialize;
        recipe->base.parameters = paramlist;
        recipe->base.frames = frameset;
        recipe->config = config;
    }

    cpl_error_code result = cpl_plugin_init(plugin, CPL_PLUGIN_API, version,
                                            plugin_type, name, synopsis, desc,
                                            author, email, cpright,
                                            NULL, NULL, deinitfunc);
    if (result != CPL_ERROR_NONE)
    {
        cpl_parameterlist_delete(paramlist);
        cpl_frameset_delete(frameset);
        /* Do not call deinitialize since we never got to complete
           initialisation successfully. Just delete the allocated memory. */
        cpl_free(plugin);
        return NULL;
    }

    return plugin;
}


/**
 * @brief Converts a JSON parse tree to a CPL plugin list.
 *
 * This function will convert a parse tree or sub tree returned by the
 * er_json_parse() function to an appropriate CPL plugin list.
 *
 * @param  parsetree  The JSON parse tree or sub tree returned by
 *                    er_json_parse().
 * @param  json       The original JSON text string from which the parse tree
 *                    was generated.
 *
 * @return  A new CPL plugin list for the parsed JSON. @c NULL is returned if
 *      an error is detected. Errors typically indicate that the parsed JSON
 *      does not actually correspond to a valid CPL plugin list.
 *
 * @note The returned object must be freed with er_json_pluginlist_delete()
 *      by the caller when it is no longer needed.
 */
cpl_pluginlist * er_json_to_pluginlist(const er_json_node * parsetree,
                                       const char * json)
{
    cpl_error_ensure(parsetree != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "Input JSON parse tree is NULL.");
    cpl_error_ensure(json != NULL, CPL_ERROR_NULL_INPUT, return NULL,
                     "JSON text is NULL.");

    /* Make sure the top level node in the parse tree is an array. */
    if (er_json_node_type(parsetree) != JSON_ARRAY)
    {
        int line, col;
        er_json_find_line_column(json, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an array at line %d column %d.",
                              line, col);
        return NULL;
    }

    cpl_pluginlist * plugins = cpl_pluginlist_new();

    /* Convert and add every plugin found in the array to the plugins list. */
    er_json_array_iterator iter = NULL;
    for (iter = er_json_node_array_begin(parsetree);
         iter != er_json_node_array_end(parsetree);
         iter = er_json_node_array_next(parsetree, iter))
    {
        const er_json_node * item = er_json_node_array_get(parsetree, iter);
        cpl_plugin * plugin = er_json_to_plugin(item, json);
        if (plugin == NULL)
        {
            er_json_pluginlist_delete(plugins);
            return NULL;
        }
        if (cpl_pluginlist_append(plugins, plugin) != CPL_ERROR_NONE)
        {
            /* Ignore result from deinitialize because we are already handling
               another error. */
            int result = plugin->deinitialize(plugin);
            (void) result;  /* Avoid warnings if the NDEBUG macro is set. */
            cpl_plugin_delete(plugin);
            er_json_pluginlist_delete(plugins);
            return NULL;
        }
    }

    return plugins;
}

/**
 * @brief Deletes a CPL plugin list and all plugins that it contains.
 *
 * @param list  The list of CPL plugins to delete.
 */
void er_json_pluginlist_delete(cpl_pluginlist * list)
{
    cpl_plugin * plugin = NULL;
    if (list == NULL) return;
    for (plugin = cpl_pluginlist_get_first(list);
         plugin != NULL;
         plugin = cpl_pluginlist_get_next(list))
    {
        /* Call the plugin's destructor if it has one. NOTE: we ignore the
           output because in our context it should always succeed. */
        if (plugin->deinitialize != NULL)
        {
            int result = plugin->deinitialize(plugin);
            (void) result;  /* Avoid warnings if the NDEBUG macro is set. */
            assert(result == 0);
        }
    }
    cpl_pluginlist_delete(list);
}

#endif /* ENABLE_PYTHON_RECIPES */

/**@}*/
