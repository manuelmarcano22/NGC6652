/* $Id: fors_dfs.c,v 1.6 2013-07-11 11:22:27 cgarcia Exp $
 *
 * This file is part of the FORS library
 * Copyright (C) 2002-2006 European Southern Observatory
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

/*
 * $Author: cgarcia $
 * $Date: 2013-07-11 11:22:27 $
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vimos_dfs.h>

#include <cpl.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

/*----------------------------------------------------------------------------*/
/**
 * @defgroup vimosdfs DFS Utilities
 *
 *   The module vimos_dfs collects medium level functions related to 
 *   DFS data IO.
 */
/*----------------------------------------------------------------------------*/

/**@{*/

#define WCS_KEYS "^((CRVAL|CRPIX|CTYPE|CDELT)[0-9]|RADECSYS|CD[0-9]_[0-9])$"

/*------------------------------------------------------------------------------
    Prototypes
 -----------------------------------------------------------------------------*/
static void
errorstate_dump_one(unsigned self, unsigned first, unsigned last);

/*------------------------------------------------------------------------------
    Implementation
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
char *vmstrlower(char *s)
{

  char *t = s;

  while (*t) {
    *t = tolower(*t);
    t++;
  }

  return s;

}

char *dfs_generate_filename_tfits(const char *category)
{
    char *filename = cpl_calloc(strlen(category) + 6, sizeof(char));

    vmstrlower(strcpy(filename, category));

    strcat(filename, ".fits");

    return filename;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Dump a single CPL error
  @param    self      The number of the current error to be dumped
  @param    first     The number of the first error to be dumped
  @param    last      The number of the last error to be dumped

  This is a callback function to cpl_errorstate_dump()

  @see cpl_errorstate_dump
 */
/*----------------------------------------------------------------------------*/
static void
errorstate_dump_one(unsigned self, unsigned first, unsigned last)
{

    const cpl_boolean is_reverse = first > last ? CPL_TRUE : CPL_FALSE;
    const unsigned    newest     = is_reverse ? first : last;
    const unsigned    oldest     = is_reverse ? last : first;
    const char      * revmsg     = is_reverse ? " in reverse order" : "";

    /* Cannot use internal CPL functions
       cx_assert( oldest <= self );
       cx_assert( newest >= self ); */

    if (newest == 0) {
        cpl_msg_info(cpl_func, "No error(s) to dump");
        /* cx_assert( oldest == 0); */
    } else {
        /*  cx_assert( oldest > 0);
            cx_assert( newest >= oldest); */
        if (self == first) {
            if (oldest == 1) {
                cpl_msg_debug(cpl_func, "Dumping all %u error(s)%s:", newest,
                              revmsg);
            } else {
                cpl_msg_error(cpl_func, "Dumping the %u most recent error(s) "
                              "out of a total of %u errors%s:",
                              newest - oldest + 1, newest, revmsg);
            }
        }

        const char *message_from_cpl = cpl_error_get_message();

        if (message_from_cpl == NULL) {
            /* This should never happen */
            cpl_msg_error(cpl_func, "Unspecified error");
        }

        /* Skip the standard (non-informative) CPL message string, 
           which usually terminates with ': '

           If no user-defined error message is given, 
           print the CPL standard message
        */
        while (*message_from_cpl != '\0' && *message_from_cpl != ':') {
            message_from_cpl += 1;
        }
        
        if (*message_from_cpl != '\0') {
            message_from_cpl += 1;
            
            if (*message_from_cpl == ' ') message_from_cpl++;
            
            if (*message_from_cpl != '\0') {
                /* Still something left of the string */
                cpl_msg_error(cpl_func, "%s [%s]", message_from_cpl,
                              cpl_error_get_where());
            }
            else {
                cpl_msg_error(cpl_func, "%s [%s]",
                              cpl_error_get_message(), cpl_error_get_where());
        }
        }
        else {
            /* Found no ':' is CPL message string */
            cpl_msg_error(cpl_func, "%s [%s]",
                          cpl_error_get_message(), cpl_error_get_where());
        }
    }

    return;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Reading a recipe integer parameter value
 *
 * @param parlist   The input parameter list
 * @param name      The parameter name
 * @param defaults  The defaults table
 *
 * @return The parameter value
 *
 * This function is just a wrapper to the basic CPL function
 * @c cpl_parameter_get_int(), but if a @em defaults table is 
 * passed then the parameter value is searched in that table.
 * If it is found, the parameter value will be modified also
 * on the input parameter list before being returned, so that
 * it would appear on the recipe products headers. If the
 * parameter is not found, then the parameter value is simply
 * read from the input parameter list. If a @em defaults table
 * is not specified, then this function works exactly as the
 * function @c cpl_parameter_get_int().
 */
/*----------------------------------------------------------------------------*/
int dfs_get_parameter_int(cpl_parameterlist *parlist, const char *name,
                          const cpl_table *defaults)
{
    const char *func = "dfs_get_parameter_int";

    const char    *alias;
    cpl_parameter *param;


    if (parlist == NULL) {
        cpl_msg_error(func, "Missing input parameter list");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    if (name == NULL) {
        cpl_msg_error(func, "Missing input parameter name");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    param = cpl_parameterlist_find(parlist, name);

    if (param == NULL) {
        cpl_msg_error(func, "Wrong parameter name: %s", name);
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return 0;
    }

    if (cpl_parameter_get_type(param) != CPL_TYPE_INT) {
        cpl_msg_error(func, "Unexpected type for parameter "
                      "\"%s\": it should be integer", name);
        cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
        return 0;
    }

    alias = cpl_parameter_get_alias(param, CPL_PARAMETER_MODE_CLI);

//    if (defaults && cpl_parameter_get_default_flag(param) == 0) {
    if (defaults && 
        cpl_parameter_get_default_int(param) == cpl_parameter_get_int(param)) {
    
        if (cpl_table_has_column(defaults, alias)) {
            if (cpl_table_get_column_type(defaults, alias) != CPL_TYPE_INT) {
                cpl_msg_error(func, "Unexpected type for CONFIG_TABLE "
                              "column \"%s\": it should be integer", alias);
                cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
                return 0;
            }
            if (cpl_table_is_valid(defaults, alias, 0)) {
                cpl_parameter_set_int(param, cpl_table_get_int(defaults, 
                                                               alias, 0, NULL));
            }
            else {
                cpl_msg_error(func, "Invalid parameter value in table "
                              "column \"%s\"", alias);
                cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
                return 0;
            }
        }
        else {
            cpl_msg_warning(func, "Parameter \"%s\" not found in CONFIG_TABLE "
                            "- using recipe default", alias);
        }
    }

    cpl_msg_info(func, "%s: %d", alias, cpl_parameter_get_int(param));

    return cpl_parameter_get_int(param);

}

/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Reading a recipe double parameter value
 *
 * @param parlist   The input parameter list
 * @param name      The parameter name
 * @param defaults  The defaults table
 *
 * @return The parameter value
 *
 * This function is just a wrapper to the basic CPL function
 * @c cpl_parameter_get_double(), but if a @em defaults table is 
 * passed then the parameter value is searched in that table.
 * If it is found, the parameter value will be modified also
 * on the input parameter list before being returned, so that
 * it would appear on the recipe products headers. If the
 * parameter is not found, then the parameter value is simply
 * read from the input parameter list. If a @em defaults table
 * is not specified, then this function works exactly as the
 * function @c cpl_parameter_get_double().
 */
/*----------------------------------------------------------------------------*/
double dfs_get_parameter_double(cpl_parameterlist *parlist, 
                                const char *name, const cpl_table *defaults)
{
    const char *func = "dfs_get_parameter_double";

    const char    *alias;
    cpl_parameter *param;


    if (parlist == NULL) {
        cpl_msg_error(func, "Missing input parameter list");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    if (name == NULL) {
        cpl_msg_error(func, "Missing input parameter name");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    param = cpl_parameterlist_find(parlist, name);

    if (param == NULL) {
        cpl_msg_error(func, "Wrong parameter name: %s", name);
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return 0;
    }

    if (cpl_parameter_get_type(param) != CPL_TYPE_DOUBLE) {
        cpl_msg_error(func, "Unexpected type for parameter "
                      "\"%s\": it should be double", name);
        cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
        return 0;
    }

    alias = cpl_parameter_get_alias(param, CPL_PARAMETER_MODE_CLI);

//    if (defaults && cpl_parameter_get_default_flag(param) == 0) {
    if (defaults && 
        cpl_parameter_get_default_double(param) == 
        cpl_parameter_get_double(param)) {
    
        if (cpl_table_has_column(defaults, alias)) {
            if (cpl_table_get_column_type(defaults, alias) != CPL_TYPE_DOUBLE) {
                cpl_msg_error(func, "Unexpected type for GRISM_TABL "
                              "column \"%s\": it should be double", alias);
                cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
                return 0;
            }
            if (cpl_table_is_valid(defaults, alias, 0)) {
                cpl_parameter_set_double(param, cpl_table_get_double(defaults, 
                                                               alias, 0, NULL));
            }
            else {
                cpl_msg_error(func, "Invalid parameter value in table "
                              "column \"%s\"", alias);
                cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
                return 0;
            }
        }
        else {
            cpl_msg_warning(func, "Parameter \"%s\" not found in CONFIG_TABLE "
                            "- using recipe default", alias);
        }
    }

    cpl_msg_info(func, "%s: %f", alias, cpl_parameter_get_double(param));

    return cpl_parameter_get_double(param);

}

/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Reading a recipe string parameter value
 *
 * @param parlist   The input parameter list
 * @param name      The parameter name
 * @param defaults  The defaults table
 *
 * @return The parameter value
 *
 * This function is just a wrapper to the basic CPL function
 * @c cpl_parameter_get_string(), but if a @em defaults table is 
 * passed then the parameter value is searched in that table.
 * If it is found, the parameter value will be modified also
 * on the input parameter list before being returned, so that
 * it would appear on the recipe products headers. If the
 * parameter is not found, then the parameter value is simply
 * read from the input parameter list. If a @em defaults table
 * is not specified, then this function works exactly as the
 * function @c cpl_parameter_get_string().
 */
/*----------------------------------------------------------------------------*/
const char *dfs_get_parameter_string(cpl_parameterlist *parlist, 
                                     const char *name, 
                                     const cpl_table *defaults)
{
    const char *func = "dfs_get_parameter_string";

    const char    *alias;
    cpl_parameter *param;


    if (parlist == NULL) {
        cpl_msg_error(func, "Missing input parameter list");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    if (name == NULL) {
        cpl_msg_error(func, "Missing input parameter name");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    param = cpl_parameterlist_find(parlist, name);

    if (param == NULL) {
        cpl_msg_error(func, "Wrong parameter name: %s", name);
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return 0;
    }

    if (cpl_parameter_get_type(param) != CPL_TYPE_STRING) {
        cpl_msg_error(func, "Unexpected type for parameter "
                      "\"%s\": it should be string", name);
        cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
        return 0;
    }

    alias = cpl_parameter_get_alias(param, CPL_PARAMETER_MODE_CLI);

//    if (defaults && cpl_parameter_get_default_flag(param) == 0) {
    if (defaults && 
        strcmp(cpl_parameter_get_default_string(param), 
               cpl_parameter_get_string(param)) == 0) {
    
        if (cpl_table_has_column(defaults, alias)) {
            if (cpl_table_get_column_type(defaults, alias) != CPL_TYPE_STRING) {
                cpl_msg_error(func, "Unexpected type for CONFIG_TABLE "
                              "column \"%s\": it should be string", alias);
                cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
                return 0;
            }
            if (cpl_table_is_valid(defaults, alias, 0)) {
                cpl_parameter_set_string(param, cpl_table_get_string(defaults, 
                                                             alias, 0));
            }
            else {
                cpl_msg_error(func, "Invalid parameter value in table "
                              "column \"%s\"", alias);
                cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
                return 0;
            }
        }
        else {
            cpl_msg_warning(func, "Parameter \"%s\" not found in CONFIG_TABLE "
                            "- using recipe default", alias);
        }
    }

    cpl_msg_info(func, "%s: %s", alias, cpl_parameter_get_string(param));

    return cpl_parameter_get_string(param);

}

/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Reading a recipe boolean parameter value
 *
 * @param parlist   The input parameter list
 * @param name      The parameter name
 * @param defaults  The defaults table
 *
 * @return The parameter value
 *
 * This function is just a wrapper to the basic CPL function
 * @c cpl_parameter_get_bool(), but if a @em defaults table is 
 * passed then the parameter value is searched in that table.
 * If it is found, the parameter value will be modified also
 * on the input parameter list before being returned, so that
 * it would appear on the recipe products headers. If the
 * parameter is not found, then the parameter value is simply
 * read from the input parameter list. If a @em defaults table
 * is not specified, then this function works exactly as the
 * function @c cpl_parameter_get_bool().
 */
/*----------------------------------------------------------------------------*/
int dfs_get_parameter_bool(cpl_parameterlist *parlist, const char *name,
                           const cpl_table *defaults)
{
    const char *func = "dfs_get_parameter_bool";

    const char    *alias;
    cpl_parameter *param;
    int            value;


    if (parlist == NULL) {
        cpl_msg_error(func, "Missing input parameter list");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    if (name == NULL) {
        cpl_msg_error(func, "Missing input parameter name");
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    param = cpl_parameterlist_find(parlist, name);

    if (param == NULL) {
        cpl_msg_error(func, "Wrong parameter name: %s", name);
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return 0;
    }

    if (cpl_parameter_get_type(param) != CPL_TYPE_BOOL) {
        cpl_msg_error(func, "Unexpected type for parameter "
                      "\"%s\": it should be boolean", name);
        cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
        return 0;
    }

    alias = cpl_parameter_get_alias(param, CPL_PARAMETER_MODE_CLI);

//    if (defaults && cpl_parameter_get_default_flag(param) == 0) {
    if (defaults && 
        cpl_parameter_get_default_bool(param) == 
        cpl_parameter_get_bool(param)) {
    
        if (cpl_table_has_column(defaults, alias)) {
            if (cpl_table_get_column_type(defaults, alias) != CPL_TYPE_INT) {
                cpl_msg_error(func, "Unexpected type for CONFIG_TABLE "
                              "column \"%s\": it should be integer", alias);
                cpl_error_set(func, CPL_ERROR_INVALID_TYPE);
                return 0;
            }
            if (cpl_table_is_valid(defaults, alias, 0)) {
                value = cpl_table_get_int(defaults, alias, 0, NULL);
                if (value < 0 || value > 1) {
                    cpl_msg_error(func, "Illegal parameter value in table "
                                  "column \"%s\": it should be either 0 or 1", 
                                  alias);
                    cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
                    return 0;
                }
                cpl_parameter_set_bool(param, value);
            }
            else {
                cpl_msg_error(func, "Invalid parameter value in table "
                              "column \"%s\"", alias);
                cpl_error_set(func, CPL_ERROR_ILLEGAL_INPUT);
                return 0;
            }
        }
        else {
            cpl_msg_warning(func, "Parameter \"%s\" not found in CONFIG_TABLE "
                            "- using recipe default", alias);
        }
    }

    value = cpl_parameter_get_bool(param);

    if (value) {
        cpl_msg_info(func, "%s: TRUE", alias);
    }
    else {
        cpl_msg_info(func, "%s: FALSE", alias);
    }

    return value;

}

/*----------------------------------------------------------------------------*/
/**
 * @see dfs_get_parameter_bool
 */
/*----------------------------------------------------------------------------*/
int dfs_get_parameter_bool_const(const cpl_parameterlist *parlist, const char *name)
{
    return dfs_get_parameter_bool((cpl_parameterlist *)parlist, name, NULL);
}

/*----------------------------------------------------------------------------*/
/**
 * @see dfs_get_parameter_int
 */
/*----------------------------------------------------------------------------*/
int dfs_get_parameter_int_const(const cpl_parameterlist *parlist, const char *name)
{
    return dfs_get_parameter_int((cpl_parameterlist *)parlist, name, NULL);
}

/*----------------------------------------------------------------------------*/
/**
 * @see dfs_get_parameter_double
 */
/*----------------------------------------------------------------------------*/
double dfs_get_parameter_double_const(const cpl_parameterlist *parlist, const char *name)
{
    return dfs_get_parameter_double((cpl_parameterlist *)parlist, name, NULL);
}

/*----------------------------------------------------------------------------*/
/**
 * @see dfs_get_parameter_string
 */
/*----------------------------------------------------------------------------*/
const char *dfs_get_parameter_string_const(const cpl_parameterlist *parlist, const char *name)
{
    return dfs_get_parameter_string((cpl_parameterlist *)parlist, name, NULL);
}

/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Loading image data of given category.
 *
 * @param frameset The input set-of-frames
 * @param category The category of the image to load
 * @param type     The data type of the loaded image
 * @param ext      The FITS file extension to access (first = 0)
 * @param calib    1 = calibration file, 0 = raw file
 *
 * @return The loaded image
 *
 * This function is just a wrapper to the basic CPL functions
 * @c cpl_frameset_find() and @c cpl_image_load(), as they 
 * typically are called every time an image should be loaded 
 * by a recipe. Error checking and proper messaging are also
 * included here, to give a more readable look to the main 
 * recipe code.
 *
 * In case of any error, a @c NULL pointer is returned. The 
 * error codes that are set in this case are the same set by
 * the above mentioned CPL functions. The "where" string
 * (accessible via a call to @c cpl_error_get_where() ) is
 * not modified by this function, and therefore the function
 * where the failure occurred is also reported.
 *
 */
/*----------------------------------------------------------------------------*/
cpl_image *dfs_load_image(cpl_frameset *frameset, const char *category, 
                          cpl_type type, int ext, int calib)
{
    const char *func = "dfs_load_image";

    cpl_frame  *frame = NULL;
    cpl_image  *image = NULL;


    frame = cpl_frameset_find(frameset, category);

    if (frame) {
        image = cpl_image_load(cpl_frame_get_filename(frame), type, 0, ext);
        if (image == NULL) {
            cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
            cpl_msg_error(func, "Cannot load image %s",
                          cpl_frame_get_filename(frame));
        }
        else {
            if (calib) 
                cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
            else
                cpl_frame_set_group(frame, CPL_FRAME_GROUP_RAW);
        }
    }

    return image;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Loading table data of given category.
 *   
 * @param frameset The input set-of-frames
 * @param category The category of the table to load
 * @param ext      The FITS file extension to access (first = 0)
 * 
 * @return The loaded table
 * 
 * This function is just a wrapper to the basic CPL functions
 * @c cpl_frameset_find() and @c cpl_table_load(), as they 
 * typically are called every time a table should be loaded
 * by a recipe. Error checking and proper messaging are also
 * included here, to give a more readable look to the main
 * recipe code.
 *
 * In case of any error, a @c NULL pointer is returned. The
 * error codes that are set in this case are the same set by
 * the above mentioned CPL functions. The "where" string
 * (accessible via a call to @c cpl_error_get_where() ) is
 * not modified by this function, and therefore the function
 * where the failure occurred is also reported.
 *
 */
/*----------------------------------------------------------------------------*/
cpl_table *dfs_load_table(cpl_frameset *frameset, const char *category, int ext)
{
    const char *func = "dfs_load_table";

    cpl_frame  *frame = NULL;
    cpl_table  *table = NULL;


    frame = cpl_frameset_find(frameset, category);

    if (frame) {
        table = cpl_table_load(cpl_frame_get_filename(frame), ext, 1);
        if (table == NULL) {
            cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
            cpl_msg_error(func, "Cannot load table %s",
                          cpl_frame_get_filename(frame));
        }
    }

    return table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Loading header associated to data of given category.
 *
 * @param frameset The input set-of-frames
 * @param category The category of the frame containing the header
 * @param ext      The FITS file extension to access (first = 0)
 *
 * @return The loaded property list
 *
 * This function is just a wrapper to the basic CPL functions
 * @c cpl_frameset_find() and @c cpl_propertylist_load(), as they
 * typically are called every time a header should be loaded
 * by a recipe. Error checking and proper messaging are also
 * included here, to give a more readable look to the main
 * recipe code.
 *
 * In case of any error, a @c NULL pointer is returned. The
 * error codes that are set in this case are the same set by
 * the above mentioned CPL functions. The "where" string
 * (accessible via a call to @c cpl_error_get_where() ) is
 * not modified by this function, and therefore the function
 * where the failure occurred is also reported.
 */
/*----------------------------------------------------------------------------*/
cpl_propertylist *dfs_load_header(cpl_frameset *frameset, 
                                  const char *category, int ext)
{
    const char *func = "dfs_load_header";

    cpl_frame         *frame = NULL;
    cpl_propertylist  *plist = NULL;


    frame = cpl_frameset_find(frameset, category);

    if (frame) {
        plist = cpl_propertylist_load(cpl_frame_get_filename(frame), ext);
        if (plist == NULL) {
            cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
            cpl_msg_error(func, "Cannot load header from %s",
                          cpl_frame_get_filename(frame));
        }
    }

    return plist;
}

    
/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Saving image data of given category.
 *
 * @param frameset   The input set-of-frames (to be upgraded)
 * @param image      The image to save
 * @param category   The category of the image to save
 * @param header     Header to input to cpl_dfs_setup_product_header()
 * @param parlist    The recipe parameter list
 * @param recipename The name of the recipe
 * @param version    The version of the pipeline
 *
 * @return 0 in case of success.
 *
 * This function is just a wrapper to the basic CPL functions
 * that are routinely called every time an image product must 
 * be saved to disk by a recipe. Error checking and proper 
 * messaging are also included here, to give a more readable 
 * look to the main recipe code.
 *
 * The output file name will be derived from the specified
 * category by lowercasing it and by appending the suffix ".fits".
 * The new image is properly logged in the input set-of-frames
 * in case of success.
 *
 * The error codes that are set in this case are the same set 
 * by the above mentioned CPL functions. The "where" string
 * (accessible via a call to @c cpl_error_get_where() ) is
 * not modified by this function, and therefore the function
 * where the failure occurred is also reported.
 */
/*----------------------------------------------------------------------------*/
int dfs_save_image(cpl_frameset *frameset, const cpl_image *image, 
                   const char *category, cpl_propertylist *header,
                   const cpl_parameterlist *parlist, const char *recipename, 
                   const char *version)
{
    const char       *func = "dfs_save_image";

    char             *filename;
    cpl_frame        *frame;
    cpl_propertylist *plist;


    if (category == NULL || frameset == NULL || image == NULL) {
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return -1;
    }

    cpl_msg_info(func, "Saving %s image to disk...", category);

    filename = cpl_calloc(strlen(category) + 6, sizeof(char));

    vmstrlower(strcpy(filename, category));
    strcat(filename, ".fits");

    frame = cpl_frame_new();

    cpl_frame_set_filename(frame, filename);
    cpl_frame_set_tag(frame, category);
    cpl_frame_set_type(frame, CPL_FRAME_TYPE_IMAGE);
    cpl_frame_set_group(frame, CPL_FRAME_GROUP_PRODUCT);
    cpl_frame_set_level(frame, CPL_FRAME_LEVEL_FINAL);
    if (cpl_error_get_code()) {
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        cpl_msg_error(func, "Cannot initialise the product frame");
        cpl_frame_delete(frame);
        cpl_free(filename);
        return -1;
    }


    /*
     * Produce DFS compliant FITS header for image
     */

    if (header == NULL)
        plist = cpl_propertylist_new();
    else
        plist = cpl_propertylist_duplicate(header);

    if (cpl_dfs_setup_product_header(plist, frame, frameset, parlist,
                                     recipename, version, "PRO-1.15", NULL)) {
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        cpl_msg_error(func, "Problem with product %s FITS header definition",
                      category);
        cpl_propertylist_delete(plist);
        cpl_frame_delete(frame);
        cpl_free(filename);
        return -1;
    }

/* CPL3.0 
    cpl_propertylist_erase_regexp(plist, "^ESO DET OUT1 OVSC*", 0);
    cpl_propertylist_erase_regexp(plist, "^ESO DET OUT1 PRSC*", 0);
*/
/* CPL2.0 */
    cpl_propertylist_erase(plist, "ESO DET OUT1 OVSCX");
    cpl_propertylist_erase(plist, "ESO DET OUT1 PRSCX");
    cpl_propertylist_erase(plist, "ESO DET OUT1 OVSCY");
    cpl_propertylist_erase(plist, "ESO DET OUT1 PRSCY");
    cpl_propertylist_erase_regexp(plist, 
"^ESO PRO CRV |^ESO PRO IDS |^ESO PRO ZERO |^ESO PRO OPT |^ESO PRO CCD |^ESO PRO SKY ", 0);

    /*
     * Write image to disk
     */

    if (cpl_image_save(image, filename, CPL_BPP_IEEE_FLOAT, plist,
                       CPL_IO_DEFAULT)) {
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        cpl_msg_error(func, "Cannot save product %s to disk", filename);
        cpl_propertylist_delete(plist);
        cpl_frame_delete(frame);
        cpl_free(filename);
        return -1;
    }

    cpl_propertylist_delete(plist);

    cpl_free(filename);

    cpl_frameset_insert(frameset, frame);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Saving table data of given category.
 *
 * @param frameset   The input set-of-frames (to be upgraded)
 * @param table      The table to save
 * @param category   The category of the table to save
 * @param header     Header to input to cpl_dfs_setup_product_header()
 * @param parlist    The recipe parameter list
 * @param recipename The name of the recipe
 * @param version    The version of the pipeline
 *
 * @return 0 in case of success.
 *
 * This function is just a wrapper to the basic CPL functions
 * that are routinely called every time a table product must 
 * be saved to disk by a recipe. Error checking and proper 
 * messaging are also included here, to give a more readable 
 * look to the main recipe code.
 *
 * The output file name will be derived from the specified
 * category by lowercasing it and by appending the suffix ".fits".
 * The new table is properly logged in the input set-of-frames
 * in case of success.
 *
 * The error codes that are set in this case are the same set 
 * by the above mentioned CPL functions. The "where" string
 * (accessible via a call to @c cpl_error_get_where() ) is
 * not modified by this function, and therefore the function
 * where the failure occurred is also reported.
 */
/*----------------------------------------------------------------------------*/
int dfs_save_table(cpl_frameset *frameset, const cpl_table *table, 
                   const char *category, cpl_propertylist *header,
                   const cpl_parameterlist *parlist, const char *recipename, 
                   const char *version)
{
    const char       *func = "dfs_save_table";

    char             *filename;
    cpl_frame        *frame;
    cpl_propertylist *plist;


    if (category == NULL || frameset == NULL || table == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        return -1;
    }

    cpl_msg_info(func, "Saving %s table to disk...", category);

    /*
    filename = cpl_calloc(strlen(category) + 7, sizeof(char));
    */
    filename = cpl_calloc(strlen(category) + 6, sizeof(char));

    vmstrlower(strcpy(filename, category));

    /*
    strcat(filename, ".tfits");
    */
    strcat(filename, ".fits");

    frame = cpl_frame_new();

    cpl_frame_set_filename(frame, filename);
    cpl_frame_set_tag(frame, category);
    cpl_frame_set_type(frame, CPL_FRAME_TYPE_TABLE);
    cpl_frame_set_group(frame, CPL_FRAME_GROUP_PRODUCT);
    cpl_frame_set_level(frame, CPL_FRAME_LEVEL_FINAL);
    if (cpl_error_get_code()) {
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        cpl_msg_error(func, "Cannot initialise the product frame");
        cpl_frame_delete(frame);
        cpl_free(filename);
        return -1;
    }


    /*
     * Produce DFS compliant FITS header for table
     */

    if (header == NULL)
        plist = cpl_propertylist_new();
    else
        plist = header;

    if (cpl_dfs_setup_product_header(plist, frame, frameset, parlist,
                                     recipename, version, "PRO-1.15", NULL)) {
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        cpl_msg_error(func, "Problem with product %s FITS header definition",
                      category);
        if (header == NULL)
            cpl_propertylist_delete(plist);
        cpl_frame_delete(frame);
        cpl_free(filename);
        return -1;
    }

/* For CPL3.0
    cpl_propertylist_erase_regexp(plist, "^ESO DET OUT1 OVSC*", 0);
    cpl_propertylist_erase_regexp(plist, "^ESO DET OUT1 PRSC*", 0);
*/
/* For CPL2.0 */
    cpl_propertylist_erase(plist, "ESO DET OUT1 OVSCX");
    cpl_propertylist_erase(plist, "ESO DET OUT1 PRSCX");
    cpl_propertylist_erase(plist, "ESO DET OUT1 OVSCY");
    cpl_propertylist_erase(plist, "ESO DET OUT1 PRSCY");

    /*
     * Write table to disk
     */
    
    
    if (cpl_table_save(table, plist, NULL, filename, CPL_IO_DEFAULT)) {
        cpl_msg_error(cpl_error_get_where(), "%s", cpl_error_get_message());
        cpl_msg_error(func, "Cannot save product %s to disk", filename);
        if (header == NULL)
            cpl_propertylist_delete(plist);
        cpl_frame_delete(frame);
        cpl_free(filename);
        return -1;
    }

    if (header == NULL)
        cpl_propertylist_delete(plist);
    cpl_free(filename);

    cpl_frameset_insert(frameset, frame);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief
 *   Saving table data of given category.
 *
 * @param frameset The input set-of-frames
 * @param keyword  The keyword that should be identical in all frames
 *
 * @return 1 if consistent, 0 if inconsistent
 *
 * The specified FITS header @em keyword should be identical
 * in all the input frames. Only the FITS primary header is
 * examined. If @em keyword is missing this is considered a 
 * case of identical keywords. Only integer and string keywords 
 * are compared: in case of other types 0 is always returned.
 * If a file is not FITS, it is ignored.
 */
/*----------------------------------------------------------------------------*/
int dfs_equal_keyword(cpl_frameset *frameset, const char *keyword)
{
    const char       *func = "dfs_equal_keyword";

    cpl_frame        *frame;
    cpl_propertylist *reference;
    cpl_type          rtype;
    cpl_type          type;
    const char       *rstring;
    const char       *string;
    int               rintero;
    int               intero;
    int               found;


    if (frameset == NULL || keyword == NULL) {
        cpl_error_set(func, CPL_ERROR_NULL_INPUT);
        return 0;
    }

    if (cpl_frameset_is_empty(frameset)) {
        cpl_error_set(func, CPL_ERROR_DATA_NOT_FOUND);
        return 0;
    }
        
    found = 0;

    for(int i = 0; i < cpl_frameset_get_size(frameset); ++i) {

        frame = cpl_frameset_get_position(frameset, i);

        reference = cpl_propertylist_load(cpl_frame_get_filename(frame), 0);
        if (cpl_error_get_code() == CPL_ERROR_BAD_FILE_FORMAT) {
            cpl_error_reset();
            continue;
        }

        if (cpl_propertylist_has(reference, keyword)) {
            rtype = cpl_propertylist_get_type(reference, keyword);

            if (rtype == CPL_TYPE_STRING) {
                found = 1;
                rstring = cpl_strdup(cpl_propertylist_get_string(reference, 
                                                                 keyword));
                cpl_propertylist_delete(reference);
                break;
            }

            if (rtype == CPL_TYPE_INT) {
                found = 1;
                rintero = cpl_propertylist_get_int(reference, keyword);
                cpl_propertylist_delete(reference);
                break;
            }

            cpl_propertylist_delete(reference);
            return 0;
        }

        cpl_propertylist_delete(reference);

    }


    if (!found)
        return 1;

    for(int i = 0; i < cpl_frameset_get_size(frameset); ++i) {

        frame = cpl_frameset_get_position(frameset, i);

        reference = cpl_propertylist_load(cpl_frame_get_filename(frame), 0);
        if (cpl_error_get_code() == CPL_ERROR_BAD_FILE_FORMAT) {
            cpl_error_reset();
            continue;
        }

        if (cpl_propertylist_has(reference, keyword)) {

            type = cpl_propertylist_get_type(reference, keyword);

            if (rtype != type) {
                cpl_propertylist_delete(reference);
                return 0;
            }

            if (rtype == CPL_TYPE_STRING) {
                string = cpl_propertylist_get_string(reference, 
                                                     keyword);
                if (strncmp(rstring, string, 15)) {
                    cpl_propertylist_delete(reference);
                    return 0;
                }
            }

            if (rtype == CPL_TYPE_INT) {
                intero = cpl_propertylist_get_int(reference, keyword);
                if (rintero - intero) {
                    cpl_propertylist_delete(reference);
                    return 0;
                }
            }
        }

        cpl_propertylist_delete(reference);
    }

    if (rtype == CPL_TYPE_STRING)
        cpl_free((void *)rstring);

    return 1;

}

/**
 * @brief Save a table in a extension (different from the first one)
 *
 * @param table      Table to save
 * @param tag        PRO.CATG of the table
 * @param extheader  Header for the extension or NULL
 *
 * @return CPL_ERROR_NONE of corresponding cpl_error_code on error.
 */
cpl_error_code dfs_save_table_ext(cpl_table        * table,
				  const char             * tag,
				  cpl_propertylist * extheader)
{
    char * filename = cpl_calloc(strlen(tag) + 6, sizeof(char));
    cpl_propertylist * header;

    if (extheader) { 
	header = cpl_propertylist_duplicate(extheader);

	cpl_propertylist_erase_regexp(header, 
                                      "^ESO DPR |^ARCFILE$|^ORIGFILE$", 0);
    } else {
	header = NULL;
    }

    vmstrlower(strcpy(filename, tag));
    strcat(filename, ".fits");
    
    if (cpl_table_save(table, NULL, header, filename, CPL_IO_EXTEND)) {
	cpl_free(filename);
	cpl_ensure_code(0, CPL_ERROR_FILE_IO);
    }
    
    cpl_propertylist_delete(header);
    cpl_free(filename);

    return CPL_ERROR_NONE;
}

/**
 * @brief Save an image in a extension 
 *
 * @param image      Image to save
 * @param tag        PRO.CATG of the image
 * @param extheader  Header for the extension or NULL
 *
 * @return CPL_ERROR_NONE of corresponding cpl_error_code on error.
 */
cpl_error_code dfs_save_image_ext(cpl_image        * image,
				  const char       * tag,
				  cpl_propertylist * extheader)
{
    char * filename = cpl_calloc(strlen(tag) + 6, sizeof(char));
    
    cpl_propertylist * header = NULL;

    if (extheader) {
	header = cpl_propertylist_duplicate(extheader);

	cpl_propertylist_erase_regexp(header, 
        "^ESO DPR |^ARCFILE$|^ORIGFILE$|^ESO PRO CRV |^ESO PRO IDS |^ESO PRO ZERO |^ESO PRO OPT |^ESO PRO CCD |^ESO PRO SKY ", 0);
    } else {
	header = NULL;
    }

    vmstrlower(strcpy(filename, tag));
    strcat(filename, ".fits");
    
    if (cpl_image_save(image, filename, CPL_BPP_IEEE_FLOAT,
		       header, CPL_IO_EXTEND)) {
	cpl_free(filename);
	cpl_ensure_code(0, CPL_ERROR_FILE_IO);
    }
    
    if (header != NULL)
        cpl_propertylist_delete(header);
    cpl_free(filename);

    return CPL_ERROR_NONE;
}

/**
 * @brief Save a product with an empty primary extension
 *
 * @param frameset   Frameset
 * @param parlist    Parlist
 * @param tag        PRO.CATG of the product
 * @param recipename Name of the recipe
 * @param version    Version
 *
 * @return CPL_ERROR_NONE of corresponding cpl_error_code on error.
 */
cpl_error_code dfs_save_image_null(cpl_frameset * frameset,
                                   cpl_propertylist * extra_header, 
                                   cpl_parameterlist * parlist,
				   const char * tag,
				   const char * recipename,
				   const char * version)
{
    const char * regexp = "ESO DET OUT1 OVSCX|"
	                  "ESO DET OUT1 PRSCX|"
                          "ESO DET OUT1 OVSCY|"
                          "ESO DET OUT1 PRSCY";

    char * filename = cpl_calloc(strlen(tag) + 6, sizeof(char));

    cpl_error_code error;

    cpl_propertylist * pro = cpl_propertylist_new();
		
    cpl_propertylist_append_string(pro, "ESO PRO CATG", tag);
    if(extra_header != NULL)
        cpl_propertylist_append(pro, extra_header);
    
    vmstrlower(strcpy(filename, tag));
    strcat(filename, ".fits");

    error = cpl_dfs_save_image(frameset, NULL, parlist, frameset, NULL, NULL,
			       CPL_BPP_IEEE_FLOAT, recipename, pro,
			       regexp, version, filename);
    
    cpl_free(filename);
    cpl_propertylist_delete(pro);

    return error;
}

void
vimos_dfs_set_groups(cpl_frameset * set)
{
    cpl_frame *f;
    int i, n;
    
    if(set == NULL)
        return;
    
    n = cpl_frameset_get_size(set);
    for (i = 0 ; i<n; i++)
    {
        f = cpl_frameset_get_position(set, i);

        const char *tag = cpl_frame_get_tag(f);

        if (tag != NULL) {
            if (strcmp(tag, "BIAS") == 0 ) {
                cpl_frame_set_group(f, CPL_FRAME_GROUP_RAW);
            }
            else if (strcmp(tag, "MASTER_BIAS") == 0) {
                cpl_frame_set_group(f, CPL_FRAME_GROUP_CALIB);
            }
            else {
                cpl_msg_warning(cpl_func, "Unrecognized frame tag: '%s'",
                                tag);
            }
        }
    }
    
    return;
}

/*----------------------------------------------------------------------------*/
/** 
 * @brief  Extract frames with given tag from frameset
 * @param  frames      frame set
 * @param  tag         to search for
 * @return newly allocated, possibly empty, frameset, or NULL on error
 */
/*----------------------------------------------------------------------------*/
cpl_frameset *
vimos_frameset_extract(const cpl_frameset *frames,
                      const char *tag)
{
    cpl_frameset *subset = NULL;
    const cpl_frame *f;

    if( frames == NULL)
        return NULL;
    if( tag    == NULL)
        return NULL;
    
    subset = cpl_frameset_new();

    for (f = cpl_frameset_find_const(frames, tag);
         f != NULL;
         f = cpl_frameset_find_const(frames, NULL)) {

        cpl_frameset_insert(subset, cpl_frame_duplicate(f));
    }

    return subset;
}

/**@}*/
