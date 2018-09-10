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

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>


#include "cpl.h"

#include "er_macros.h"
#include "er_help.h"
#include "er_plugin.h"
#include "er_paramutils.h"
#include "er_stringarray.h"
#include "er_stringutils.h"
#include "er_fileutils.h"
#include "er_json.h"

#define SIZE_X    1024 

static  char  *wdsc = NULL;         /* space for help description */
static  char  *entrx = NULL;        /* space for help entry */

static int wdsclen = 0;
static int entrlen = 0;




/**
 * @defgroup er_help  EsoRex help functions.
 *
 * Generates help information for EsoRex functionality.
 */

/**@{*/

/**********************************************************************/
/**
 * @brief               Fills a string with a "usage" description.
 * @param param         Parameter for which a description is to be created.
 * @param print_default  TRUE = display default values
 *
 * Generates a string containing the help information about a parameter.
 * The specification of the default value is option and is controlled by
 * the @a print_default boolean argument.
 */
/**********************************************************************/

static void er_help_create_description (cpl_parameter  *param,
                                        int   print_default) 
{
    const char  *cpp, *cstr = NULL;
    char  *cptr;

    int   indx;
    int   countr, n, str_len;

    cpl_type      type;            /* Parameter type */

    int        f_env_active;    /* Flag for if an ENV variable is available */


    if (param == NULL)
    {
        cpl_msg_warning(er_func, "A NULL-parameter was found");
        return;
    }

    if (wdsclen == 0)
    {
        wdsc = (char *) cpl_malloc((size_t) SIZE_X);
        if (wdsc == NULL)
        {
            cpl_msg_warning (er_func,
                             "Could not allocate initial %d bytes for help description",
                             SIZE_X);
            return;
        }
        wdsclen = SIZE_X;
    }

    cpp = cpl_parameter_get_help(param);
    if (cpp == NULL)
    {
        cpl_msg_warning(er_func, "cpl_parameter_get_help failed ");
        return;
    }

    str_len = (int) strlen(cpp);
    n = str_len + 1;
    if (n > wdsclen)
    {                                      /* allocate already more space */
        wdsclen = n + 160;
        er_enlarge("help_create",&wdsc,wdsclen);
    }


    (void) strcpy(wdsc, cpp);
    cptr = wdsc + str_len - 1;
    for (n=str_len; n>0; n--)
    {
        if (isspace((int)*cptr))
            cptr --;
        else
            break;
    }

    if (*cptr != '.') *++cptr = '.';    /* append a full-stop, if none */
    *++cptr = '\0';
    countr = (int) strlen(wdsc);        /* length counter for wdsc */


    /* Add environment variable, if any */

    f_env_active = FALSE;
    f_env_active = cpl_parameter_is_enabled(param, CPL_PARAMETER_MODE_ENV);

    if (f_env_active == TRUE)
    {
        const char  *env_par;
        char        *pt1, *pt2;

        /* get name for eventual env. variable (a.b.c if nothing defined) */
        env_par = cpl_parameter_get_alias(param, CPL_PARAMETER_MODE_ENV);

        pt1 = strchr(env_par,'.');      /* get first '.' */
        pt2 = strrchr(env_par,'.');     /* get last '.' */
        if ((pt1 != NULL) && (pt2 != NULL) && (pt1 < pt2)) {
            ;       /* it's the a.b.c string - no variable defined */
        }
        else
        {
            char  tmp[512];

            (void)strcpy(tmp,
                         " This option may also be set using the environment variable ");
            (void)strcat(tmp, env_par);
            (void)strcat(tmp, ".");
            n = (int) strlen(tmp) + 1;
            countr += n;
            if (countr > wdsclen)
            {                                      /* allocate already more space */
                wdsclen = countr + 160;
                er_enlarge("help_create",&wdsc,wdsclen);
            }
            (void) strcat(wdsc,tmp);
        }
    }

    /* Add options, depending on the type */

    type = cpl_parameter_get_type(param);
    switch(type)
    {
        case CPL_TYPE_BOOL:
        case CPL_TYPE_INT:  /* ????? - Need to add enumeration/range handling */
        case CPL_TYPE_DOUBLE:   /* ????? - Need to add enumeration/range handling */
            break;

        case CPL_TYPE_STRING:
            if (cpl_parameter_get_class(param) == CPL_PARAMETER_CLASS_ENUM)
            {
                int  klen;
                char  tmp[SIZE_X];

                (void)strcpy(tmp, " <");
                klen = 4;               /* the < and > chars */
                for(indx = 0; indx < cpl_parameter_get_enum_size(param); indx++)
                {
                    if (indx != 0)
                    {
                        klen += 3;
                        (void)strcat(tmp, " | ");
                    }
                    cstr = cpl_parameter_get_enum_string(param, indx);
                    if (cstr == NULL)
                    {
                        cpl_msg_warning(er_func, "cpl_parameter_get_enum_string failed ");
                        return;
                    }
                    klen += 4 + (int) strlen(cstr);     /* for the ` | ' above */
                    if (klen > SIZE_X)
                    {
                        cpl_msg_warning(er_func, "internal options overflow...");
                        return;
                    }
                    else
                        (void)strcat(tmp, cstr);
                }

                (void)strcat(tmp, ">");
                n = (int) strlen(tmp) + 1;
                countr += n;
                if (countr > wdsclen)
                {                                      /* allocate more space */
                    wdsclen = countr + 160;
                    er_enlarge("help_create",&wdsc,wdsclen);
                }
                (void) strcat(wdsc,tmp);
            }
            break;

        default:
            cpl_msg_warning(er_func, "Unrecognized parameter type");
            break;
    }


    /* Add defaults, if required */

    if (print_default == TRUE)
    {
        char  tmp[MAXSTRLENCONF];
        const char *cppp;

        switch(type)
        {
            case CPL_TYPE_BOOL:
                switch(cpl_parameter_get_bool(param))
                {
                    case TRUE:
                        (void)strcpy(tmp, " [TRUE]");
                        countr += 7;
                        break;

                    case FALSE:
                        (void)strcpy(tmp, " [FALSE]");
                        countr += 8;
                        break;

                    default:
                        (void)strcpy(tmp, " [???]");
                        countr += 6;
                        break;
                }
                break;

                    case CPL_TYPE_INT:
                        n = snprintf(tmp,(size_t)80, " [%d]", cpl_parameter_get_int(param));
                        countr += n;
                        break;

                    case CPL_TYPE_DOUBLE:
                        (void)strcpy(tmp, " [");
                        cppp = er_strutils_dblstr(cpl_parameter_get_double(param));
                        if ((int)strlen(cppp) > (MAXSTRLENCONF - 4))
                        {
                            cpl_msg_warning(er_func, "cpl_parameter_get_double(param) > 1000 chars...");
                            break;
                        }
                        (void)strcat(tmp, cppp);
                        (void)strcat(tmp, "]");
                        countr += (int) strlen(tmp);
                        break;

                    case CPL_TYPE_STRING:
                        (void)strcpy(tmp, " [");
                        cppp = cpl_parameter_get_string(param);
                        if ((int)strlen(cppp) > (MAXSTRLENCONF - 4))
                        {
                            cpl_msg_warning(er_func, "cpl_parameter_get_string(param) > 1000 chars...");
                            break;
                        }
                        (void)strcat(tmp, cppp);
                        (void)strcat(tmp, "]");
                        countr += (int) strlen(tmp);
                        break;

                    default:
                        cpl_msg_warning(er_func, "Unrecognized parameter type");
                        break;
        }

        if (countr > wdsclen)
        {                                      /* allocate more space */
            wdsclen = countr + 160;
            er_enlarge(er_func,&wdsc,wdsclen);
        }
        (void) strcat(wdsc,tmp);
    }
    return;

}


/**********************************************************************/
/**
 * @brief               Fills a string with a "usage" description.
 * @param entry         An existing string, into which the keyword and
 *                      default value shall be placed.
 * @param param         Parameter for which a description is to be created.
 *
 * Generates a string containing the full parameter name and default value.
 * The specification of the default value is option and is controlled by
 * the @a print_default boolean argument.
 */
/**********************************************************************/

static void er_help_create_entry (cpl_parameter  *param)

{
    char  value[SIZE_X-80];     /* Default value of an option */
    const char *cpp;

    int  n;

    cpl_type  type;         /* Parameter type */


    value[0] = '\0';

    if (param == NULL)
    {
        cpl_msg_warning(er_func, "A NULL-parameter was found");
        return;
    }

    cpp = cpl_parameter_get_name(param);
    if (cpp == NULL)
    {
        cpl_msg_warning(er_func, "cpl_parameter_get_name failed ");
        return;
    }
    n = (int) strlen(cpp) + 24;     /* for the stuff below, except TYPE_STRING */
    if (n > entrlen)
    {
        entrlen = n;
        er_enlarge("help_create_entry",&entrx,entrlen);
    }

    (void) strcpy(entrx, cpp);
    (void) strcat(entrx, "=");


    /* Add the default value, depending on type */

    type = cpl_parameter_get_type(param);
    switch(type)
    {
        case CPL_TYPE_BOOL:
            switch(cpl_parameter_get_bool(param))
            {
                case TRUE:
                    (void)strcat(entrx, "TRUE");
                    break;

                case FALSE:
                    (void)strcat(entrx, "FALSE");
                    break;

                default:
                    cpl_msg_warning(er_func, "Unable to set one of the parameter values");
                    (void)strcat(entrx, "???");
                    break;
            }
            break;

                case CPL_TYPE_INT:
                    (void)snprintf(value,(size_t)80, "%d", cpl_parameter_get_int(param));
                    (void)strcat(entrx,value);
                    break;

                case CPL_TYPE_DOUBLE:
                    (void)strcat(entrx,
                                 er_strutils_dblstr(cpl_parameter_get_double(param)) );
                    break;

                case CPL_TYPE_STRING:
                    cpp = cpl_parameter_get_string(param);
                    ForceNull(cpp)
                    if (cpp == NULL) break;

                    n = (int) strlen(cpp) + 1;
                    if (n > entrlen)
                    {
                        entrlen = n;
                        er_enlarge("help_create_entry",&entrx,entrlen);
                    }
                    (void)strcat(entrx, cpp);
                    break;

                default:
                    cpl_msg_warning(er_func, "Unrecognized parameter type");
                    break;
    }

    return;
}


/**********************************************************************/
/**
 * @brief             Provides help information on EsoRex
 * @param param_list  List of Parameters
 * @param show_all    Flag to display disabled ('hidden') parameters
 *
 * Prints the help information for EsoRex (but not the plugin, if any).
 */
/**********************************************************************/

static void er_help_print_paramlist(cpl_parameterlist *param_list,
                                    int show_all)
{

    /*
     * Loop through all the parameters in the list and create the
     * description text.
     */

    cpl_parameter *param = cpl_parameterlist_get_first(param_list);

    while (param != NULL)
    {
        const char *alias = cpl_parameter_get_alias(param,
                                                    CPL_PARAMETER_MODE_CLI);

        int enabled = cpl_parameter_is_enabled(param, CPL_PARAMETER_MODE_CLI);

        if ((enabled || show_all) && (alias != NULL)) {
            er_help_create_description(param, TRUE);
            er_paramutils_print_key_desc(COMMAND_LINE_PREFIX, alias, wdsc);
        }

        /*
         * Get the next parameter
         */

        param = cpl_parameterlist_get_next(param_list);
    }

    /* Pad with an additional blank line */
    printf( "\n" );

} 


/**********************************************************************/
/**
 * @brief             Counts the number of parameters in a list.
 * @param param_list  List of Parameters
 *
 * This function is used to count the number of parameters that is
 * within a specified parameter list.
 */
/**********************************************************************/

static int er_help_count_paramlist ( cpl_parameterlist *param_list)

{
    int    num_params = 0;   /* Counter */

    cpl_parameter   *param  = NULL;   /* Parameter within the list */


    /* Loop through all the parameters in the list */

    param = cpl_parameterlist_get_first(param_list);
    while (param != NULL)
    {
        num_params++;
        param = cpl_parameterlist_get_next(param_list);
    }

    /* Return the number of parameters counted */
    return (num_params);

}


/**********************************************************************/
/**
 * @brief             Provides help information on EsoRex
 * @param plugin_name     The name of the plugin
 * @param param_list      List of the plugin parameters
 * @param flag_show_all   Flag to force displaying disabled parameters.
 *
 * Prints the help information for EsoRex (but not the plugin, if any).
 */
/**********************************************************************/

void er_help_display(const char *plugin_name, cpl_parameterlist *param_list,
                     int flag_show_all)

{
    /* Print some header text */
    (void) printf("Usage: %s [%s-options] %s [%s-options] sof\n\n",
                  PACKAGE, PACKAGE, plugin_name, plugin_name);

    (void) printf("Options:\n\n");

    er_help_print_paramlist(param_list, flag_show_all);
}

/**********************************************************************/
/**
 * @brief             Prints a simulated man page for a plugin.
 * @param recipe          The recipe plugin for which the page should be generated.
 * @param flag_show_all   Flag to force displaying disabled parameters.
 *
 * Prints a "man page" for a plugin.
 */
/**********************************************************************/

void er_help_manpage(cpl_recipe  *recipe, int flag_show_all)

{
    int left_margin;     /* Column of left margin */
    int right_margin;    /* Column of right margin */
    int margin = 2;      /* Size of margins from either side */
    int num_params;         /* Number of options provided by the recipe */
    int  n, strsize;

    char *version_str;      /* Plugin version (as a string) */
    char *str = NULL;
    const char  *cpp;


    /* Check that we can actually do something */

    if (recipe == NULL)
    {
        cpl_msg_warning(er_func, "Unable to print documentation for the recipe");
        return;
    }


    /* Determine the right hand margin */

    right_margin = er_strutils_termwidth();
    if (right_margin < 1) right_margin = 80;
    right_margin -= margin;

    str = (char *) cpl_malloc((size_t) SIZE_X);
    if ( str == NULL)
    {
        cpl_msg_error (er_func, "Could not allocate %d bytes for str",SIZE_X);
        return ;
    }
    strsize = SIZE_X;           /* dynamic size */


    /* Plugin name */

    (void) printf("\nNAME\n");
    (void) strcpy(str, "  ");

    cpp = cpl_plugin_get_name((cpl_plugin *)recipe);
    if (cpp == NULL)
    {
        cpl_msg_error (er_func, "cpl_plugin_get_name() failed");
        goto free_str;
    }

    n = 12 + (int) strlen(cpp);     /* add the constant char stuff before + after */
    if (n > strsize)
    {
        strsize = n;
        er_enlarge("help_manpage",&str,strsize);
    }
    (void) strcat(str, cpp);

    (void) strcat(str, " -- ");
    left_margin = (int) strlen(str);    /* new length of string */

    cpp = cpl_plugin_get_synopsis((cpl_plugin *)recipe);
    if (cpp == NULL)
    {
        cpl_msg_error (er_func, "cpl_plugin_get_synopsis() failed");
        goto free_str;
    }

    n = (int) strlen(cpp) + left_margin;
    if (n > strsize)
    {
        strsize = n;
        er_enlarge("help_manpage",&str,strsize);
    }
    (void) strcat(str, cpp);

    (void) printf("%s\n\n", er_strutils_split(str , left_margin, right_margin) );
    left_margin = margin;

    /* Print some header text */

    (void) printf("SYNOPSIS\n");
    (void) printf("  %s [%s-options] %s [%s-options] sof\n\n",
                  PACKAGE, PACKAGE, cpl_plugin_get_name((cpl_plugin *)recipe),
                  cpl_plugin_get_name((cpl_plugin *)recipe));

    /* Description of the plugin */

    (void) printf("DESCRIPTION\n");
    (void) strcpy(str, "  ");

    cpp = cpl_plugin_get_description((cpl_plugin *)recipe);
    if (cpp == NULL)
    {
        cpl_msg_error (er_func, "cpl_plugin_get_description() failed");
        goto free_str;
    }

    n = 4 + (int) strlen(cpp);
    if (n > strsize)
    {
        strsize = n;
        er_enlarge("help_manpage",&str,strsize);
    }
    (void) strcat(str, cpp);

    (void) printf("%s\n\n", er_strutils_indent( str, left_margin) );


    /* Plugin options */

    printf("OPTIONS\n");
    num_params = er_help_count_paramlist(recipe->parameters);

    if (num_params < 1)
        printf("  No options are provided by this recipe.\n\n");
    else
    {
        char  tmp[256];         /* actually 194 are needed... */

        if (num_params > 1)
            (void) printf("  The following options are provided by this recipe.\n\n");
        else
            (void) printf("  A single option is provided by this recipe.\n\n");

        er_help_print_paramlist(recipe->parameters, flag_show_all);

        (void) strcpy(tmp, "  Note that it is also possible to create a "
                      "configuration file containing these options, along with suitable "
                      "default values. Please refer to the "
                      "details provided by the '" PACKAGE " --help' command.");
        printf("%s\n\n", er_strutils_split(tmp , left_margin, right_margin) );
    }


    /* Author's name */

    printf("AUTHOR\n");
    strcpy(str, "  ");

    cpp = cpl_plugin_get_author((cpl_plugin *)recipe);
    if (cpp == NULL)
    {
        cpl_msg_error (er_func, "cpl_plugin_get_author() failed");
        goto free_str;
    }

    n = 4 + (int) strlen(cpp);
    if (n > strsize)
    {
        strsize = n;
        er_enlarge("help_manpage",&str,strsize);
    }
    (void) strcat(str, cpp);
    printf("%s\n\n", er_strutils_split(str , left_margin, right_margin) );


    /* Bugs and E-mail address */

    (void) printf("BUGS\n");
    (void) strcpy(str, "  Please report any problems with this recipe to ");
    n = 2 + (int) strlen(str);              /* don't forget the \0 */

    cpp = cpl_plugin_get_email((cpl_plugin *)recipe);
    if (cpp == NULL)
    {
        cpl_msg_error (er_func, "cpl_plugin_get_email() failed");
        goto free_str;
    }

    n += (int) strlen(cpp);
    if (n > strsize)
    {
        strsize = n;
        er_enlarge("help_manpage",&str,strsize);
    }
    (void) strcat(str, cpp);
    (void) printf("%s\n\n", er_strutils_split(str , left_margin, right_margin) );


    /* Version */

    (void) printf("RELEASE\n");
    (void) strcpy(str, "  ");

    cpp = cpl_plugin_get_name((cpl_plugin *)recipe);    /* already checked above */

    n = 24 + (int) strlen(cpp);     /* const chars before + after ... */
    if (n > strsize)
    {
        strsize = n;
        er_enlarge("help_manpage",&str,strsize);
    }
    (void) strcat(str, cpp);
    (void) strcat(str, " -- version ");

    version_str = cpl_plugin_get_version_string((cpl_plugin *)recipe);
    if (version_str == NULL)
    {
        cpl_msg_error (er_func, "cpl_plugin_get_version_string() failed");
        goto free_str;
    }

    n += (int) strlen(version_str);
    if (n > strsize)
    {
        strsize = n;
        er_enlarge("help_manpage",&str,strsize);
    }
    (void) strcat(str, version_str);
    cpl_free(version_str);

    (void) printf("%s\n\n", er_strutils_split( str , left_margin, right_margin) );


    /* License */

    (void) printf("LICENSE\n");
    cpp = cpl_plugin_get_copyright((cpl_plugin *)recipe);
    if (cpp == NULL)
    {
        cpl_msg_error (er_func, "cpl_plugin_get_copyright() failed");
        goto free_str;
    }
    (void) printf("  %s\n", er_strutils_indent(cpp , left_margin));


    free_str:
    cpl_free(str);
    return;

}


/**********************************************************************/
/**
 * @brief            Creates header text for a configuration file.
 * @param fp         File pointer, into which the text shall be written.
 * @param file_name  Name of the file
 *
 * Writes to the file pointer @a fp a block of header comments to mark
 * the start of a EsoRex created configure file.
 */
/**********************************************************************/

static void er_help_create_config_header (FILE *fp, const char *file_name)

{
    time_t epoch;
    struct tm *timestruct;        /* Time structure */

    char timestr[MAXSTRLENCONF];  /* String with date/time */


    fprintf(fp, "# File: %s\n", file_name);
    fprintf(fp, "#\n");
    fprintf(fp, "# Note: This configuration file has been automatically\n");
    fprintf(fp, "%s%s) program.\n", ER_HELP_CFG_VER, VERSION);
    fprintf(fp, "#\n");

    epoch = time(NULL);
    timestruct = localtime(&epoch);

    (void)strftime(timestr, MAXSTRLENCONF, "%d-%b-%Y %H:%M:%S", timestruct);
    fprintf(fp, "# Date: %s\n", timestr);
    fprintf(fp, "#\n");

    fprintf(fp, "#\n\n");

}


/**********************************************************************/
/**
 * @brief             Creates a list of config-file parameter entries.
 * @param param_list  List of Parameters
 *
 * Generates a sequence of parameter entries suitable for a 
 * configuration file.
 */
/**********************************************************************/

static void er_help_create_paramlist (FILE *fp, cpl_parameterlist *param_list)

{
    cpl_parameter *param = NULL;        /* Parameter */

    char  *split_str = NULL;        /* Split string */
    char  **value_str = NULL;       /* String of just the parameter value */

    int  c;                 /* Loop counter */


    /* Loop through all the parameters in the list */

    if (entrlen == 0)
    {
        entrx = (char *) cpl_malloc((size_t) SIZE_X);
        if (entrx == NULL)
        {
            cpl_msg_warning (er_func, "Could not allocate initial %d bytes for help entry",
                             SIZE_X);
            return;
        }
        entrlen = SIZE_X;
    }


    param = cpl_parameterlist_get_first(param_list);
    while (param != NULL)
    {
        /*
         * Only parameters which are not disabled are written to
         * the configuration file.
         */

        if (cpl_parameter_is_enabled(param, CPL_PARAMETER_MODE_CFG))
        {

            /*
             * Create the comment and entry text
             */

            er_help_create_description(param, FALSE);
            er_help_create_entry(param);

            /*
             * Display it on the terminal
             */

            fprintf(fp, "# %s%s\n", COMMAND_LINE_PREFIX,
                    cpl_parameter_get_alias(param, CPL_PARAMETER_MODE_CLI));

            split_str = er_strutils_split((const char *) wdsc, 2,
                                          DEFAULT_TERM_WIDTH-2);

            /*
             * Parse the string and replace any space following a
             * new-line with a #
             */

            for (c = 0; c < (int)(strlen(split_str)-1); c++)
            {
                if (*(split_str+c) == '\n')
                {
                    if (*(split_str+c+1) == ' ') *(split_str+c+1) = '#';
                }
            }

            /*
             * Write the actual comment (including the initial #)
             * into the file
             */

            fprintf(fp, "# %s\n", split_str);

            /*
             * Write the actual entry - but comment out blank ones!
             */

            value_str = cx_strsplit(entrx, "=", 0);
            if ( ( value_str == NULL ) ||
                    ( *(value_str+1) == NULL ) ||
                    ( **(value_str+1) == '\0'))
            {
                const char *alias =
                        cpl_parameter_get_alias(param, CPL_PARAMETER_MODE_CLI);

                cpl_msg_warning(er_func,
                                "No sensible default available for '%s%s'. "
                                "The entry will be inserted into the "
                                "configuration file, but will be commented "
                                "out.", COMMAND_LINE_PREFIX, alias);

                fprintf(fp, "#");
            }
            fprintf(fp, "%s\n\n", entrx);

            /* Get the next parameter */

            cx_strfreev(value_str); value_str = NULL;
        }

        param = cpl_parameterlist_get_next(param_list);
    }

    /* Pad with an additional blank line and end of file message */
    fprintf(fp, "#\n");
    fprintf(fp, "# End of file\n");
}


/**********************************************************************/
/**
 * @brief  Generates a default configuration file.
 * @param def_flag    = 1 for default , = 11 for specific config file
 * @param context_name    Name of the application context
 * @param file_name       Name of the file to be created (without 
 *                        extension)
 * @param caller_parlist  List of Parameters from the caller, which
 *                        will govern the behaviour of this function.
 * @param config_parlist  List of the parameters that shall be 
 *                        written into the configuration file.
 *
 * Generates a configuration file containing all the parameters from
 * the specified list, along with relevant comments and the default
 * values.
 */
/**********************************************************************/

void er_help_create_config ( int def_flag, 
                             const char  *context_name, 
                             const char  *file_name, 
                             cpl_parameterlist *caller_parlist, 
                             cpl_parameterlist *config_parlist)
{
    char final_file_name[FILEMAX+1];
    char command_str[2*FILEMAX];
    char  *cpp;

    int  file_exists = TRUE;

    FILE *fp = NULL;

    int status;

    CPL_UNUSED(caller_parlist);


    /*  check if default (=> $HOME/.esorex/esorex.rc) or user specified config file   */

    if (def_flag == 11)
    {
        (void) strcpy(final_file_name, file_name);
        goto filename_complete;     /* skip the default stuff */
    }

    /* First, generate the global directory name, and establish the directory */

    cpp = getenv("HOME");
    if (cpp == NULL)
        (void)strcpy( final_file_name, "/" );
    else
    {
        (void) strcpy(final_file_name,cpp);
        (void) strcat( final_file_name, "/" );
    }
    (void)strcat( final_file_name, GLOBAL_RC_DIR );

    if (fileutils_directory_exists(final_file_name) == FALSE)
    {
        cpl_msg_warning(er_func, "No global configuration directory exists. "
                        "Creating global configuration directory '%s'",
                        final_file_name);

        (void) strcpy(command_str,"mkdir ");
        (void) strcat(command_str,final_file_name);
        status = system(command_str);
        if(status != 0) cpl_msg_error(er_func, "Error creating directory '%s'",
                                      command_str);
        if (fileutils_directory_exists(final_file_name) == FALSE)
        {
            cpl_msg_error(er_func, "Unable to create global configuration "
                          "directory '%s'\n", final_file_name );
            return;
        }
    }

    /* Now complete the filename */

    (void)strcat( final_file_name, "/" );
    (void)strcat( final_file_name, context_name);
    (void)strcat( final_file_name, GLOBAL_RC_EXTENSION);


    filename_complete:

    /* Check for file existence */
    file_exists = fileutils_file_exists(final_file_name);

    /* Make a copy of the old configuration file */
    if (file_exists == TRUE)
    {
        (void) strcpy(command_str,"cp ");
        (void)strcat(command_str,final_file_name);
        (void)strcat(command_str," ");
        (void)strcat(command_str,final_file_name);
        (void)strcat(command_str, GLOBAL_RC_BACKUP);
        if (system(command_str) != 0)
        {
            cpl_msg_error(er_func, "Unable to make a backup copy of the "
                          "configuration file.");
            cpl_msg_warning(er_func, "No attempt made to create a "
                            "new configuration file.");
            return;
        }
        else
        {
            cpl_msg_info(er_func, "Old configuration file copied to '%s%s'",
                         final_file_name, GLOBAL_RC_BACKUP);
        }
    }

    /* If we're clear to go, then report it, and being the generation process */
    cpl_msg_info(er_func, "Creating configuration file '%s'", final_file_name);

    if(strlen(final_file_name)>5 &&
            strncmp(final_file_name + strlen(final_file_name) - 5, ".json", 5) == 0 )
    {
        /* We write a JSON format configuration file */
        cpl_msg_info(cpl_func, "Writing configuration file in JSON format");
        er_recipe_parameterlist_to_json(config_parlist, context_name,
                                        final_file_name);
    }
    else
    {
        /* We write a plain text configuration file */

        /* Open the file */
        fp = fopen(final_file_name, "w");
        if (fp == NULL)
        {
            cpl_msg_error(er_func, "Unable to create configuration file '%s'",
                          final_file_name);
            return;
        }

        /* Generate the header text */
        er_help_create_config_header(fp, final_file_name);

        /* Generate the actually configuration contents */
        er_help_create_paramlist(fp, config_parlist);

        /* Close the file and remove the filename string */
        fclose(fp);
    }

}

void er_help_free(void)

{
    if (entrlen > 0) cpl_free(entrx);
    if (wdsclen > 0) cpl_free(wdsc);
}

/**@}*/

/* End of file */
