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
#   include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef HAVE_GETOPT_LONG_ONLY
#   include <getopt.h>
#else
#   include "getopt.h"
#endif

#include "ltdl.h"


#include "er_macros.h"
#include "er_help.h"
#include "er_stringarray.h"
#include "er_plugin.h"
#include "er_pluginlist.h"
#include "er_params.h"
#include "er_paramutils.h"
#include "er_fileutils.h"
#include "er_stringutils.h"


/**
 * @defgroup er_params_process EsoRex Parameter Processing Functions
 *
 * EsoRex Parameter Processing Functions
 *
 */

/**@{*/

/**********************************************************************/
/**
 * @brief
 *   Parse esorex caller configuration file
 *
 * @param param_list  List of Parameters
 * @param file_name   Name of File to parse
 *
 * @return ?????
 *
 * Function searches in file @em file_name for keyword, value pairs
 * If a keyword is found that is present in the @em param_list the value
 * found is stored as current value of the parameter in the parameter list
 * @em param_list.
 */
/**********************************************************************/

int params_parse_config_file (cpl_parameterlist *param_list,
                              const char *file_name)

{
    const char *f_file_kw;
    char input[MAXSTRLENCONF];
    char temp[MAXSTRLENCONF];                 /* Temporary copy of "input" */
    char **kv_pair;

    cpl_parameter *p;

    size_t value_length, i;

    int f_file_active;
    int num_values;             /* The number of values read with sscanf() */
    int maj_cfg_version;        /* Maj.version No. from config file */
    int maj_pkg_version;        /* Maj.version No. of this executable */

    FILE *file_descriptor;


    if (file_name == NULL) return -1;

    file_descriptor = fopen (file_name, "r");
    if (file_descriptor == NULL) return -1;


    /* Initialize buffer */

    temp[0] = '\0';
    memset (input, 0, MAXSTRLENCONF);


    /* Start reading lines from file */

    while (fgets (input, MAXSTRLENCONF, file_descriptor) != NULL)
    {
        /* Strip leading and trailing whitespace from entire line */
        (void) cx_strstrip (input);

        /* Check if we have an old configuration file */

        if (strncmp (input, ER_HELP_CFG_VER, strlen (ER_HELP_CFG_VER)) == 0)
        {
            /*
             * If we get to this part, then we have the comment
             * about the version
             */
            strcpy (temp, input);

            /* Convert all the non-numbers to white space */
            for (i = 0; i < strlen (temp); i++)
            {
                if (!isdigit ((int) temp[i]))
                    temp[i] = ' ';
            }

            num_values = sscanf (temp, "%d", &maj_cfg_version);
            if (num_values == 1)
            {
                num_values = sscanf (VERSION, "%d", &maj_pkg_version);
                if (num_values != 1)
                    cpl_msg_error (er_func, "Unidentifed " PACKAGE " version");
                if (maj_pkg_version > maj_cfg_version)
                {
                    cpl_msg_warning (er_func,
                                     "The configuration file '%s' appears "
                                     "to have been created with an older "
                                     "version of EsoRex.", file_name);
                }
                if (maj_pkg_version < maj_cfg_version)
                {
                    cpl_msg_warning (er_func,
                                     "The configuration file '%s' appears "
                                     "to have been created with a more "
                                     "recent version of EsoRex.", file_name);
                }
            }
        }                         /* End of configuration version checking */


        /*
         *      Pos 1 = '#', line is remark,
         *      Pos 1 = '\0', line is empty/blank,
         *      otherwise try to find keyword=value
         */

        if ((input[0] != '#') && (input[0] != '\0'))
        {
            /* Try to split this line into a key-value pair at an "=" sign */
            if ((kv_pair = cx_strsplit (input, "=", 2)))
            {
                /*
                 * Strip leading and trailing whitespace from keyword and
                 * value
                 */
                if (kv_pair[0]) (void) cx_strstrip (kv_pair[0]);
                if (kv_pair[1]) (void) cx_strstrip (kv_pair[1]);
            }

            /* Ensure that key-value splitting was successful */
            /* Ensure that both the key and the value are set to something */
            if ((kv_pair == NULL)
                    || (kv_pair[0] == NULL) || (kv_pair[1] == NULL)
                    || (kv_pair[0][0] == '\0') || (kv_pair[1][0] == '\0'))
            {
                cpl_msg_error (er_func, "'%s' does not fit 'key=value' "
                               "format in configuration file '%s'", input,
                               file_name);

                cpl_msg_error (er_func, "Unable to continue parsing "
                               "configuration file '%s'",
                               file_name);

                fclose (file_descriptor);

                if (kv_pair) cx_strfreev (kv_pair);     /* free array */
                return -2;
            }

            /*
             * Found key-value pair, try to match it with accepted file
             * keywords
             */

            int f_unknown_param = -1;

            p = cpl_parameterlist_get_first(param_list);

            while (p != NULL)
            {
                /* Over-ride them with values */
                f_file_active = cpl_parameter_is_enabled(p, CPL_PARAMETER_MODE_CFG);
                f_file_kw = cpl_parameter_get_name(p);

                if (strcmp(kv_pair[0], f_file_kw) == 0)
                {
                    if (f_file_active != 0)
                    {
                        f_unknown_param = 0;
                        value_length = strlen (kv_pair[1]);

                        /* Handle quoted values */
                        if (kv_pair[1][0] == '"')
                        {
                            if ((kv_pair[1][value_length - 1] == '"')
                                    && (kv_pair[1][value_length - 2] != '\\'))
                            {
                                /* Matching quotes at start and end - remove them */
                                memmove (&kv_pair[1][0], &kv_pair[1][1],
                                         value_length - 2);
                                kv_pair[1][value_length - 2] = '\0';
                                value_length -= 2;
                            }
                            else
                            {
                                /* Quotes are mismatched */
                                cpl_msg_error (er_func,
                                               "Quoting mismatch within value "
                                               "'%s' for parameter '%s' in "
                                               "configuration file '%s'",
                                               kv_pair[1], kv_pair[0], file_name);

                                cpl_msg_error (er_func, "Unable to continue "
                                               "parsing configuration file '%s'",
                                               file_name);

                                fclose (file_descriptor);

                                cx_strfreev (kv_pair);
                                return -3;
                            }
                        }

                        /* Strip backslash-escaping out of the value string */
                        for (i=0; i<value_length; ++i)
                        {
                            if (kv_pair[1][i] == '\\')
                            {               /* Backslash found - remove it */
                                memmove (&kv_pair[1][i], &kv_pair[1][i + 1],
                                         value_length - 1 - i);
                                kv_pair[1][value_length - 1] = '\0';
                                value_length--;
                            }
                            else if (kv_pair[1][i] == '"')
                            {               /* Unescaped quote found */
                                cpl_msg_error (er_func,
                                               "Unescaped quote within value '%s' for "
                                               "parameter '%s' in configuration file '%s'",
                                               kv_pair[1], kv_pair[0], file_name);

                                cpl_msg_error (er_func,
                                               "Unable to continue parsing "
                                               "configuration file '%s'",
                                               file_name);

                                fclose (file_descriptor);

                                cx_strfreev (kv_pair);
                                return -4;
                            }
                        }

                        /* Assign the value to the parameter */

                        if (paramutils_set_from_string(p, kv_pair[1], file_name) != 0)
                        {
                            cpl_msg_error(er_func,
                                          "Error parsing/interpreting entry '%s'"
                                          " in configuration file '%s'",
                                           cpl_parameter_get_name(p), file_name);
                            fclose(file_descriptor);

                            cx_strfreev(kv_pair);
                            return -5;
                        }
                    }
                    else
                    {
                        f_unknown_param = 1;
                    }
                }

                p = cpl_parameterlist_get_next(param_list);
            }

            if (f_unknown_param < 0)
            {
                cpl_msg_error(er_func, "Parameter '%s' not recognised "
                              "in configuration file '%s'", kv_pair[0],
                              file_name);

                cpl_msg_error(er_func, "Unable to continue parsing "
                              "configuration file '%s'", file_name);

                fclose(file_descriptor);

                cx_strfreev(kv_pair);
                return -6;
            }
            else if (f_unknown_param > 0)
            {
                cpl_msg_warning(er_func, "Parameter '%s' should not appear "
                                "in a configuration file, but is present in "
                                "'%s'", kv_pair[0], file_name);
            }

            cx_strfreev(kv_pair);
        }

        /* Reinitialize buffer */
        memset(input, '\0', MAXSTRLENCONF);
    }

    fclose(file_descriptor);

    return 0;

}                               /* End of params_parse_config_file() */


/**********************************************************************/
/**
 * @brief
 *   Parse caller configuration environment
 *
 * @param param_list  List of Parameters
 *
 * @return None
 *
 * Function searches shell environment for keyword, value pairs
 * If a keyword is found that is present in the @em param_list the value
 * found is stored as current value of the parameter in the parameter list
 * @em param_list.
 */
/**********************************************************************/

int params_parse_config_environment (cpl_parameterlist * param_list)

{
    char *value;

    cpl_parameter *p;

    int f_env_active;


    p = cpl_parameterlist_get_first (param_list);

    while (p != NULL)
    {
        f_env_active = 0;
        f_env_active = cpl_parameter_is_enabled (p, CPL_PARAMETER_MODE_ENV);

        if (f_env_active != 0)
        {                   /* we have an env-var for that parameter */
            const char  *envvar;

            envvar = cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_ENV);
            value = getenv (envvar);        /* get value of that variable */

            if (value != NULL)
            {
                int  err;
                char   str[256];

                err = (int) strlen(envvar);
                if (err > 220)
                {
                    cpl_msg_warning (er_func, "Name of env. variable > "
                                     "220 chars - skipped");
                    return err;
                }

                (void) sprintf (str, "environment variable '%s'", envvar);
                err = paramutils_set_from_string (p, value, str);
                if (err != 0) return err;
            }
        }

        p = cpl_parameterlist_get_next (param_list);
    }

    return 0;

} /* End of params_parse_config_environment() */


/**********************************************************************/
/**
 * @brief
 *   Convert cpl parameter list to genopt options list (PRIVATE FUNCTION)
 *
 * @param param_list  List of Parameters
 * @return List of getopt options. This pointer must be deallocated by 
 *         the caller.
 */
/**********************************************************************/

static struct option *
param_list_to_cmdl_opt_list (cpl_parameterlist * param_list)

{
        int countparam = 0, i = 0;
        int  alias_flag = 0;

        char  *cptr;

        size_t strlength = 0;

        struct option *newoptions = NULL;
        cpl_parameter *p;


        /*
         * if env. var ESOREX_USE_PARAM_NAME is defined, we use not the
         *  alias but the full name of parameters
         */

        cptr = getenv ("ESOREX_USE_PARAM_NAME");
        if (cptr == NULL)
            alias_flag = 1;
        else
            alias_flag = 0;

        countparam = cpl_parameterlist_get_size (param_list);

        newoptions = (struct option *)
                     cpl_malloc ((size_t) sizeof (struct option) * (countparam + 1));

        if (newoptions != NULL)
        {
            p = cpl_parameterlist_get_first (param_list);

            while (p != NULL)
            {
                if (alias_flag == 1)
                    strlength =
                            strlen (cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_CLI));
                else
                    strlength = strlen (cpl_parameter_get_name (p));

                newoptions[i].name =
                        (char *) cpl_malloc ((size_t) (sizeof (char) * (strlength + 1)));

                if (alias_flag == 1)
                    (void) strcpy ((char *) newoptions[i].name,
                                   cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_CLI));
                else
                    (void) strcpy ((char *) newoptions[i].name,
                                   cpl_parameter_get_name (p));

                newoptions[i].has_arg = optional_argument;
                newoptions[i].flag = NULL;
                newoptions[i].val = cpl_parameter_get_id (p);

                p = cpl_parameterlist_get_next (param_list);

                i++;
            }

            /* add closing entry to option list */
            newoptions[i].name = 0;
            newoptions[i].has_arg = 0;
            newoptions[i].flag = 0;
            newoptions[i].val = 0;
        }

        return newoptions;

}                               /* End of param_list_to_cmdl_opt_list() */


/**********************************************************************/
/**
 * @brief
 *   Parse esorex caller commandline
 *
 * @param param_list  List of Parameters
 * @param argc        Count of arguments given on command line
 * @param argv        Array of arguments given on command line
 * @param f_caller    Boolean indicating if parsing should begin from the
 *                    recipe name argument or not.
 *
 * @retval recipe_name Name of the recipe given on the command line
 * @retval set_of_frames_filenames Set of strings containing filenames of
 * set of frames
 * @return None
 *
 * Function parses command file arguments given in @em argv for keyword-value
 * pairs. If a keyword is found that is present in the @em param_list the
 * value found is stored as current value of the parameter in the parameter
 * list @em param_list.
 * The recipe name found is returned in @em recipe_name.
 * @em set_of_frames_filenames contains a set of filenames, filenames
 * point to a set-of-frames file.
 */
/**********************************************************************/

int
params_parse_config_commandline(cpl_parameterlist *param_list,
                                char *recipe_name,
                                er_stringarray_t *set_of_frames_filenames,
                                int argc, char *argv[],
                                int f_caller)

{
    struct option *getopt_input_list = NULL;

    cpl_parameter *p;
    cpl_type ptype;

    int f_unknown_param = 0;
    int err = 0, option_id = 0;
    int option_index = 0, i, n, iparam = 0, countparam = 0;


    /*  Number all parameters uniquely so we can find them using getopt */
    p = cpl_parameterlist_get_first (param_list);

    i = 0;
    while (p != NULL)
    {
        cpl_parameter_set_id (p, i++);
        p = cpl_parameterlist_get_next (param_list);
    }

    /*  Generate temporary opt array so we can use getopt */

    getopt_input_list = param_list_to_cmdl_opt_list (param_list);

    opterr = 0;      /* suppress getopt_warnings */
    optind = 0;      /* make sure we can enter this routine multiple times */

    /*
     * If called for plugin options, advance argv and reduce argc so that
     * getopt() will start processing from the plugin name.
     */

    if (f_caller == 0)
    {
        const char* tmp_name_ptr = recipe_name;

        for (i=0; i<argc; ++i)
        {
            n = (int) strlen(argv[i]);
            if (n > FILEMAX)
            {
                cpl_msg_error (er_func, "length of recipe name "
                               "(= %d chars)\n(name: %44.44s...)\nis > max. "
                               "size (= %d)", n, recipe_name,FILEMAX);

                return -1;
            }

            if (!strcmp (argv[i], tmp_name_ptr))
            {
                argv += i;
                argc -= i;
                break;
            }
        }
    }

    while (1)
    {
        option_id = getopt_long_only(argc, argv, "+",
                                     getopt_input_list, &option_index);

        if (option_id == -1) break;         /* end of --option list */

        f_unknown_param = 1;

        /* = 63 ('?') and 0 => mistyped arguments (Error) */
        if (option_id != option_index)
            goto check_getopt;

        p = cpl_parameterlist_get_first (param_list);

        while (p != NULL)
        {                       /* Found keyword in list */
            if ((int) cpl_parameter_get_id (p) == option_id)
            {
                f_unknown_param = 0;

                /*
                 * Cmdl            Optarg
                 * --------------- ------------------------------------------
                 * --blah         optarg = NULL
                 * --blah=        optarg = "" (string length 0)
                 * --blah=xxxxx   optarg = "xxxxx" (string length > 0)
                 */

                if (optarg == NULL)
                {
                    ptype = cpl_parameter_get_type (p);
                    if (ptype == CPL_TYPE_BOOL)
                    {
                        (void) paramutils_set_from_string (p, "TRUE",
                                                           "command line");
                    }
                    else if (ptype == CPL_TYPE_STRING)
                    {
                        (void) paramutils_set_from_string (p, "",
                                                           "command line");
                    }
                    else
                    {
                        cpl_msg_error (er_func,
                                       "'--%s' needs an argument but none was specified...",
                                       cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_CLI));
                    }
                }
                else if (strlen (optarg) == 0)
                {
                    cpl_msg_error (er_func,
                                   "'--%s' needs an argument but none was "
                                   "specified after the equals sign (\"=\").",
                                   cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_CLI));

                    countparam = cpl_parameterlist_get_size (param_list);
                    for (iparam = 0; iparam < countparam; iparam++)
                    {
                        cpl_free ((char *) getopt_input_list[iparam].name);
                    }
                    cpl_free (getopt_input_list);

                    cpl_msg_error (er_func, "'%s' during processing of '%s'",
                                   cpl_error_get_message (),
                                   cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_CLI));

                    return -1;

                }
                else
                {
                    err = paramutils_set_from_string (p, optarg, "command line");
                    if (err != 0)
                    {
                        countparam = cpl_parameterlist_get_size (param_list);
                        for (iparam = 0; iparam < countparam; iparam++)
                        {
                            cpl_free ((char *) getopt_input_list[iparam].name);
                        }
                        cpl_free (getopt_input_list);

                        cpl_msg_error (er_func,
                                       "CPL error '%s' occurred during the processing of the "
                                       "option '--%s'",
                                       cpl_error_get_message (),
                                       cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_CLI));

                        cpl_msg_warning (er_func,
                                         "Check that the supplied value is of a "
                                         "suitable type for the parameter (e.g. enum, string, "
                                         "boolean, etc.) and is within any required range.");

                        return err;
                    }

                }
            }

            p = cpl_parameterlist_get_next (param_list);
        }


        check_getopt:
        if (f_unknown_param)
        {
            countparam = cpl_parameterlist_get_size (param_list);
            for (iparam = 0; iparam < countparam; iparam++)
            {
                cpl_free ((char *) getopt_input_list[iparam].name);
            }

            cpl_free (getopt_input_list);
            cpl_msg_error(er_func,"Command line parameter '%s' not "
                          "recognized", argv[optind - 1]);
            cpl_msg_error (er_func, "Unable to continue parsing command "
                           "line options");

            return -1;
        }
    }

    /* Deallocate the parameter list */

    countparam = cpl_parameterlist_get_size (param_list);
    for (iparam = 0; iparam < countparam; iparam++)
    {
        cpl_free ((char *) getopt_input_list[iparam].name);
    }
    cpl_free (getopt_input_list);


    /* Start checking for the plugin name */

    if (optind < argc)
    {
        if (f_caller != 0)
        {           /* Assume first non-cmdl parameter is plugin name */
            n = (int) strlen(argv[optind]);
            if (n > FILEMAX)
            {
                cpl_msg_error (er_func,
                               "length of recipe name (= %d chars)\n"
                               "(name: %44.44s...)\nis > max. size (= %d)",
                               n,argv[optind],FILEMAX);
                return -1;
            }

            (void) strcpy(recipe_name, argv[optind]);
        }
        else
        { /* Assume all remaining cmdl parameter are set of frame filenames */
            while (optind < argc)
            {
                er_stringarray_append (set_of_frames_filenames, argv[optind]);
                optind++;
            }
        }

    }

    return 0;

} /* End of params_parse_config_commandline() */


/**********************************************************************/
/**
 * @brief
 *   Postprocess parameter list after reading configuration data from file,
 *   environment and commandline
 *
 * @param param_list  List of Parameters
 *
 * @return None
 *
 * Function searches @em param_list for parameter for which no argument was
 * given either in a file, in a environment variable or on the command line.
 * If such a parameter is found the default value is copied to the current
 * value.
 *
 */
/**********************************************************************/

void
params_parse_config_postprocess(cpl_parameterlist *param_list)

{
    const char *s_default;

    cpl_parameter *p = NULL;
    cpl_type ptype;

    int b_default;
    int i_default;
    int f_is_present = 0;

    double d_default;


    p = cpl_parameterlist_get_first (param_list);

    while (p != NULL)
    {
        f_is_present = cpl_parameter_get_default_flag (p);

        if (!f_is_present)
        {
            ptype = cpl_parameter_get_type (p);

            switch (ptype)
            {
                case CPL_TYPE_BOOL:
                    b_default = cpl_parameter_get_default_bool (p);
                    cpl_parameter_set_bool (p, b_default);
                    break;

                case CPL_TYPE_INT:
                    i_default = cpl_parameter_get_default_int (p);
                    cpl_parameter_set_int (p, i_default);
                    break;

                case CPL_TYPE_DOUBLE:
                    d_default = cpl_parameter_get_default_double (p);
                    cpl_parameter_set_double (p, d_default);
                    break;

                case CPL_TYPE_STRING:
                    s_default = cpl_parameter_get_default_string (p);
                    if (s_default == NULL)
                    {
                        cpl_parameter_set_string (p, "");
                    }
                    else
                    {
                        cpl_parameter_set_string (p, s_default);
                    }
                    break;

                default:
                    break;

            }                      /* End of switch(ptype) */

        }

        p = cpl_parameterlist_get_next (param_list);
    }

} /* End of params_parse_config_postprocess() */


/**********************************************************************/
/**
 * @brief
 *   Validates and handles parameters
 *
 * @param param_list    List of Parameters
 * @param plugin_name   Name of recipe specified on command line
 *
 * @retval 0 if successfull, !=0 otherwise
 *
 * Function performs the necessary operations based on the values
 * stored in the caller parameter list. All parameters are validated and
 * those which can immediately be processed will be. All others are postponed
 * until the Plugin is executed.
 */
/**********************************************************************/

int
params_handle_parameters(char *plugin_name, cpl_parameterlist *param_list)

{
    cpl_parameter *p = NULL;

    cpl_msg_severity msg_level = CPL_MSG_ERROR;

    er_stringarray_t *plugin_list = NULL, *pllib_list = NULL;

    const char *val_string_tmp = NULL;
    const char *val_create = NULL;
    const char *val_config = NULL;
    const char *val_plugin_config = NULL;
    const char *val_output_dir = NULL;
    const char *val_output_mask = NULL;
    const char *val_link_dir = NULL;
    const char *val_log_file = NULL;
    const char *val_log_dir = NULL;
    const char *val_paf_config = NULL;
    const char *val_products_sof = NULL;

    char *val_plugin_dir = NULL;
    char **val_recipe_dirs = NULL;

    int err = 0;
    int f_val_config = 0,
            f_val_recipe_dirs = 0,
            f_val_plugin_config = 0,
            f_val_version = 0,
            f_val_plugins = 0,
            f_val_help = 0,
            f_val_create = 0,
            f_val_params = 0,
            f_val_nolink = 0,
            f_val_noprefix = 0,
            f_val_output_dir = 0,
            f_val_link_dir = 0,
            f_val_output_mask = 0,
            f_val_log_file = 0,
            f_val_log_dir = 0,
            f_val_paf_config = 0,
            f_val_products_sof = 0,
            f_continue = 1;


    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".create-config");
    if (cpl_parameter_get_default_flag (p))
    {
        val_create = cpl_parameter_get_string (p);
        ForceNull(val_create)
        /* --create-config=bla.bla  found */
        if (val_create != NULL)
        {
            /* handle FALSE and TRUE from boolean history of this param */
            if (strcmp(val_create,"TRUE") == 0)
            {
                f_val_create = 1;
            }
            else if (strcmp(val_create,"FALSE") != 0)
            {
                f_val_create = 11;      /* indicates filename is given */
            }
        }
        else
        {
            /* --create-config  is interpreted as ...=TRUE */
            f_val_create = 1;
        }
    }


    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".config");
    if (cpl_parameter_get_default_flag (p))
    {
        val_config = cpl_parameter_get_string (p);
        ForceNull(val_config)
        if (val_config != NULL)
        {
            f_val_config = fileutils_file_exists (val_config);
        }
    }
    else
        f_val_config = 1;

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".recipe-dir");
    val_string_tmp = cpl_parameter_get_string (p);
    if ((val_string_tmp != NULL) && (*val_string_tmp != '\0'))
    {
        char **dir_str;
        int  nodirs, single;

        /* extract the individual directories from command-line argument */
        val_recipe_dirs = cx_strsplit (val_string_tmp, ":", -1);

        nodirs = single = 0;                /* get total no. of recipe dirs */
        for (dir_str = val_recipe_dirs; *dir_str != NULL; dir_str++)
        {
            single ++;          /* total count of dirs */
            f_val_recipe_dirs = fileutils_directory_exists (*dir_str);
            if (f_val_recipe_dirs != 0)
            {               /* one more existing directory */
                nodirs ++;
            }
        }
        if (nodirs < 1)
        {
            char  blastr[16];

            if (single < 2)
                (void) strcpy(blastr,"directory");
            else
                (void) strcpy(blastr,"directories");
            cpl_msg_warning (er_func, "Only non-existent %s, '%s', specified "
                             "as argument to '--recipe-dir'. "
                             "(Hint: use the '--params' option to check on how "
                             "EsoRex or the shell interpreted any environment "
                             "variables, etc.)", blastr,val_string_tmp);
        }
    }


    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".recipe-config");
    if (cpl_parameter_get_default_flag (p))
    {
        val_plugin_config = cpl_parameter_get_string (p);
        ForceNull(val_plugin_config)
        if (val_plugin_config != NULL)
        {
            f_val_plugin_config = fileutils_file_exists (val_plugin_config);
        }
    }
    else
        f_val_plugin_config = 1;


    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".version");
    f_val_version = cpl_parameter_get_bool (p);

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".recipes");
    f_val_plugins = cpl_parameter_get_bool (p);

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".help");
    f_val_help = cpl_parameter_get_bool (p);

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".params");
    f_val_params = cpl_parameter_get_bool (p);

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".output-dir");
    val_output_dir = cpl_parameter_get_string (p);
    ForceNull(val_output_dir)
    if (val_output_dir != NULL)
    {
        f_val_output_dir = fileutils_directory_exists ( val_output_dir);
    }

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".output-prefix");
    val_output_mask = cpl_parameter_get_string (p);
    ForceNull(val_output_mask)
    if (val_output_mask != NULL)
    {
        f_val_output_mask = 1;
    }


    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".suppress-prefix");
    f_val_noprefix = cpl_parameter_get_bool (p);

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".link-dir");
    val_link_dir = cpl_parameter_get_string (p);
    ForceNull(val_link_dir)
    if (val_link_dir != NULL)
    {
        f_val_link_dir = fileutils_directory_exists (val_link_dir);
    }

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".suppress-link");
    f_val_nolink = cpl_parameter_get_bool (p);

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".log-file");
    val_log_file = cpl_parameter_get_string (p);
    ForceNull(val_log_file)
    if (val_log_file != NULL)
    {
        f_val_log_file = 1;
    }

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".log-dir");
    val_log_dir = cpl_parameter_get_string (p);
    ForceNull(val_log_dir)
    if (val_log_dir != NULL)
    {
        f_val_log_dir = fileutils_directory_exists (val_log_dir);
    }

    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".paf-config");
    if (cpl_parameter_get_default_flag (p))
    {
        val_paf_config = cpl_parameter_get_string (p);
        ForceNull(val_paf_config)
        if (val_paf_config != NULL)
        {
            f_val_paf_config = fileutils_file_exists (val_paf_config);
        }
    }
    else
        f_val_paf_config = 1;


    p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".products-sof");
    if (cpl_parameter_get_default_flag (p))
    {
        val_products_sof = cpl_parameter_get_string (p);
        ForceNull(val_products_sof )
        if (val_products_sof  != NULL)
        {
            f_val_products_sof = fileutils_file_exists (val_products_sof );
        }
    }
    else
        f_val_products_sof = 1;

    /* Now deal with the parameters individually */

    /* esorex.caller.msg-level */

    msg_level = message_severity (param_list,2);
    cpl_msg_set_level (msg_level);


    /* esorex.caller.version */

    if ((f_continue != 0) && (f_val_version != 0))
    {
        const char *cdescr;

        f_continue = 0;
        err = -99999;

        cdescr = cpl_get_description(CPL_DESCRIPTION_DEFAULT);
        (void) printf("Libraries used: %s\n\n", cdescr);

        (void) printf("Copyright 2003-2018, European Southern Observatory.\n");
        (void) printf("Report bugs to <cpl-help@eso.org>.\n\n");
    }

    /* esorex.caller.help */

    if ((f_continue != 0) && (f_val_help != 0))
    {
        if (*plugin_name == '\0')
        {
            er_help_display ("recipe", param_list, 0);
        }
    }


    /*
     *  esorex.caller.params
     */

    if ((f_continue != 0) && (f_val_params != 0))
    {
        (void) er_paramutils_print_list (param_list, "Caller Parameters");
    }

    /*
     *  esorex.caller.config
     */

    if ((f_continue != 0) && (f_val_config == 0))
    {
        if (val_config == NULL)
        {
            cpl_msg_error (er_func, "'--config' was specified, but no argument "
                           "was given.");
        }
        else
        {
            cpl_msg_error (er_func, "'--config=%s' was specified, but the file "
                           "'%s' could not be found.",
                           val_config, val_config);
        }
        f_continue = 0;
        err = -1;
        goto clean_up;
    }

    /*
     *  esorex.caller.recipe-config
     */

    if ((f_continue != 0) && (f_val_plugin_config == 0))
    {
        if (val_plugin_config == NULL)
        {
            cpl_msg_error (er_func,
                           "'--recipe-config' was specified, but no "
                           "argument given.");
        }
        else
        {
            cpl_msg_error (er_func,
                           "'--recipe-config=%s' was specified, but the file "
                           "'%s' could not be found.", val_plugin_config,
                           val_plugin_config);
        }
        f_continue = 0;
        err = -1;
        goto clean_up;
    }

    /*
     *  esorex.caller.log-dir
     */

    if ((f_continue != 0) && (f_val_log_dir == 0))
    {
        if (val_log_dir == NULL)
        {
            cpl_msg_error (er_func, "'--log-dir' is a mandatory option, but "
                           "no argument was given. Check that there is only "
                           "an equals sign (\"=\") between the '--log-dir' "
                           "and the directory path, and that all environment "
                           "variables are resolved.");
        }
        else
        {
            cpl_msg_error (er_func, "'--log-dir=%s' specified, but '%s' is "
                           "not a valid directory.", val_log_dir,
                           val_log_dir);
        }
        f_continue = 0;
        err = -1;
        goto clean_up;
    }

    /*
     *  esorex.caller.log-file
     */

    if ((f_continue != 0) && (f_val_log_file == 0))
    {
        err = -99999;
        if (val_log_file == NULL)
        {
            cpl_msg_error (er_func, "'--log-file' is a mandatory option, but "
                           "no argument was given.  Check that there is only "
                           "an equals sign (\"=\") between the '--log-file' "
                           "and the file path, and that all environment "
                           "variables are resolved.");
            f_continue = 0;
            err = -1;
            goto clean_up;
        }
    }

    /*
     *  esorex.caller.time
     */

    /*
     *  esorex.caller.overwrite
     */

    /*
     *  esorex.caller.recipes
     */

    if ((f_continue != 0) && (f_val_plugins != 0))
    {
        plugin_list = er_stringarray_new ();
        pllib_list = NULL;
        pllib_list = er_pluginlist_create_list (val_recipe_dirs);

        er_pluginlist_create_cache (pllib_list, plugin_list);
        er_pluginlist_print_list (pllib_list);

        er_stringarray_delete (plugin_list);
        er_stringarray_delete (pllib_list);

        err = -99999;
        f_continue = 0;

    }

    /*
     *  esorex.caller.output-dir
     *  esorex.caller.output-prefix
     *  esorex.caller.suppress-prefix
     *  esorex.caller.output-rbs
     */

    if ((f_continue != 0) && (f_val_help == 0))
    {
        if (f_val_output_dir == 0)
        {
            if ((val_output_dir != NULL) && (strlen (val_output_dir) > 0))
            {
                cpl_msg_error (er_func, "The argument to '--output-dir' [%s] "
                               "is not a valid directory", val_output_dir);
            }
            else
            {
                cpl_msg_error (er_func, "No destination given for program "
                               "output (use the option '--output-dir').");
            }
            f_continue = 0;
            err = -1;
            goto clean_up;
        }

        if ((f_val_output_mask == 0) && (f_val_noprefix == 0))
        {
            cpl_msg_warning (er_func,
                             "Output file prefixes have been disabled using "
                             "the '--suppress-prefix' option.");
            cpl_msg_error (er_func,
                           "However, no valid prefix has been given using the "
                           "'--output-prefix' option.");
            cpl_msg_error (er_func,
                           "Use '--help' (without specifying any recipe) "
                           "for help on these command options");
            f_continue = 0;
            err = -1;
            goto clean_up;
        }

    }

    /*
     * esorex.caller.link-dir
     * esorex.caller.suppress-link
     */

    if ((f_continue != 0) && (f_val_help == 0) && (f_val_nolink == 0))
    {
        if (f_val_link_dir == 0)
        {
            if (val_link_dir != NULL)
            {
                cpl_msg_error (er_func,
                               "The argument to '--link-dir' [%s] is not a "
                               "valid directory", val_link_dir);
            }
            cpl_msg_error (er_func,
                           "Non-valid directory given for symbolic links.");
            f_continue = 0;
            err = -1;
            goto clean_up;
        }
    }

    /*
     *  esorex.caller.create-config
     */

    if ((f_continue != 0) && (f_val_create != 0))
    {
        if (*plugin_name == '\0')
        {                   /* no plugin (recipe) given */
            cpl_parameter *param;
            size_t  k;
            char *saved_name;

            /* Deactivate the "create-config" parameter while we do this */
            param = cpl_parameterlist_find (param_list,
                                            PACKAGE_RESOURCE ".create-config");

            if (f_val_create == 11)     /* Generate the configuration file */
            {
                k = strlen(val_create)+1;
                saved_name = malloc(k);
                (void) strcpy(saved_name,val_create);
                cpl_parameter_set_string (param, "FALSE");
                er_help_create_config (f_val_create,plugin_name, saved_name,
                                       param_list, param_list);
            }
            else
            {
                saved_name = malloc(8);
                (void) strcpy(saved_name,"TRUE");
                cpl_parameter_set_string (param, "FALSE");
                er_help_create_config (f_val_create,PACKAGE, "",
                                       param_list, param_list);
            }

            cpl_parameter_set_string (param, saved_name);
            free(saved_name);
        }
    }

    /*
     *  esorex.caller.paf-config
     */

    if ((f_continue != 0) && (f_val_paf_config == 0))
    {
        if (val_paf_config == NULL)
        {
            cpl_msg_error (er_func, "'--paf-config' was specified, but no "
                           "argument given. Check that there is only an "
                           "equals sign (\"=\") between the '--paf-config' "
                           "and the directory path, and that all environment "
                           "variables are resolved.");
        }
        else
        {
            cpl_msg_error (er_func, "'--paf-config=%s' specified, but '%s' "
                           "could not be opened for reading.", val_paf_config,
                           val_paf_config);
        }
        f_continue = 0;
        err = -1;
        goto clean_up;
    }

    /*
     *  esorex.caller.products-sof
     */

    if ((f_continue != 0) && (f_val_products_sof == 0))
    {
        if (val_products_sof == NULL)
        {
            cpl_msg_error (er_func, "'--products-sof' was specified, but no "
                           "argument given. Check that there is only an "
                           "equals sign (\"=\") between the '--products-sof' "
                           "and the directory path, and that all environment "
                           "variables are resolved.");
            f_continue = 0;
            err = -1;
            goto clean_up;
        }
    }

    /*
     * Clean-up
     */

    clean_up:
    if (val_recipe_dirs != NULL) cx_strfreev (val_recipe_dirs);
    if (val_plugin_dir != NULL) cpl_free (val_plugin_dir);

    return err;

} /* End of params_handle_parameters() */


/**********************************************************************/
/**
 * @brief
 *   Processes configuration data for caller
 *
 * @param caller_parameter_list  List of Parameters necessary for caller
 * @param global_conf_file       FQFN for global configuration file
 * @param local_conf_file        FQFN for local i.e. user specific
 *                               configuration file
 * @param argc                   Count of commandline arguments
 * @param argv                   Handle to commandline aguments
 *
 * @retval plugin_name             Name of recipe specified on command line
 * @retval set_of_frames_filenames Array of SOF filenames specified on
 *                                 commandline
 *
 * @returns 0 if successfull, !=0 otherwise
 */
/**********************************************************************/

int
params_process_configuration(cpl_parameterlist *caller_parameter_list,
                             char *global_conf_file,
                             char *local_conf_file,
                             int argc,
                             char *argv[],
                             char *plugin_name,
                             er_stringarray_t *set_of_frames_filenames)

{
    int err = 0, err2 = 0;


    /*  Global, local and environment config are not necessary as such,
    so we discard the error code! */

    /* Parse the global configuration file */
    err2 = params_parse_config_file (caller_parameter_list, global_conf_file);
    if (err2 == -5) return err;

    er_paramutils_tilde_convert (caller_parameter_list);

    /* Parse the local configuration file */
    err2 = params_parse_config_file (caller_parameter_list, local_conf_file);
    if (err2 == -5) return err;

    er_paramutils_tilde_convert (caller_parameter_list);

    /* Parse any enviroment variables */
    err2 = params_parse_config_environment (caller_parameter_list);
    if (err2 != 0) return err;

    er_paramutils_tilde_convert (caller_parameter_list);

    /* Now parse the command line variables */
    err = params_parse_config_commandline (caller_parameter_list,
                                           plugin_name,
                                           set_of_frames_filenames,
                                           argc, argv, 1);

    /* Check the error code from parsing the command line */
    if (err != 0) return err;

    er_paramutils_tilde_convert (caller_parameter_list);

    params_parse_config_postprocess (caller_parameter_list);
    er_paramutils_tilde_convert (caller_parameter_list);


    /*
  Now actually handle all these parameters 
  and return a flag indicating whether or not we should run the plugin 
     */

    err = params_handle_parameters (plugin_name, caller_parameter_list);
    return err;

}                               /* End of params_process_configuration() */

/**@}*/

/* End of file */
