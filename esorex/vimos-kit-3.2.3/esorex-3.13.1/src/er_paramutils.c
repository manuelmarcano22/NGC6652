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
#  include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <wordexp.h>


#include <cpl.h>

#include "er_macros.h"
#include "er_plugin.h"
#include "er_paramutils.h"
#include "er_stringutils.h"

#define  SIZE_X  200
#define  MAXSTR  MAXSTRLENCONF - 4

/**
 * @defgroup er_params_utils EsoRex Parameter Utility Functions
 *
 * EsoRex Parameter Utility Functions
 *
 */

/**@{*/

/**********************************************************************/
/**
 * @brief  Replaces any tilde in parameters within a parameter list.
 *
 * @param param_list  List of Parameters
 *
 * @return None
 *
 * Function searches @em param_list for parameters whose values begin 
 * with a tilde (~) and if found, replace the tilde with the user's
 * home directory (using the $HOME environment variable).
 * 
 */
/**********************************************************************/

void er_paramutils_tilde_convert (cpl_parameterlist * param_list)

{
    cpl_parameter *p = NULL;
    cpl_type ptype;

    int f_is_present = 0;       /* FALSE */

    const char *old_str = NULL;


    p = cpl_parameterlist_get_first (param_list);

    while (p != NULL)
    {
        f_is_present = cpl_parameter_get_default_flag (p);

        if (f_is_present)
        {
            ptype = cpl_parameter_get_type (p);
            if (ptype == CPL_TYPE_STRING)
            {
                old_str = cpl_parameter_get_string (p);
                if (*old_str == '~')
                {
                    char  *cpp;
                    wordexp_t  result;

                    /*
                     * use C library (glibc) routine wordexp() to expand
                     *  ~/ or ~name/
                     */

                    if (wordexp(old_str,&result,0) != 0)
                    {
                        cpl_msg_error (er_func,
                                       "Expansion of ~ (tilde) failed...");
                        return;
                    }

                    /* use 1st and only result word */
                    cpp = *result.we_wordv;
                    if (cpp == NULL)
                    {
                        cpl_msg_error (er_func, "Env. variable HOME not set");
                        return;
                    }

                    cpl_parameter_set_string (p, cpp);
                    wordfree(&result);
                }
            }
        }

        /* Get the next parameter in the list */
        p = cpl_parameterlist_get_next (param_list);
    }           /* End of loop through all parameters in the list */

}




/**********************************************************************/
/**
 * @brief Neatly print a parameter and its description.
 *
 * @param flag_prefix This is the prefix used to flag certain types of
 *                    switches/parameters. For example a double hypen
 *                    "--" may be used before the long form of command
 *                    line options.
 * @param keyword     This is the keyword itself. It is appended to the
 *                    flag-prefix, if any.
 * @param description This is the descriptive text that follows the 
 *                    keyword. It is put at a tab-stop, so that all the
 *                    descriptions are neatly aligned.
 *
 * This function neatly prints a keyword and its associated description
 * to the normal output channel. The formatting is handled automatically
 * to get neat alignment of the output.
 */
/**********************************************************************/

void
er_paramutils_print_key_desc(const char *flag_prefix, const char *keyword,
                             const char *description)
{
    int i, ii, tab = COMMENT_TAB_POSITION;

    char *mystr, *ptr, *mykeyword;


    mystr = (char *) cpl_malloc((size_t) SIZE_X);
    if (mystr == NULL)
    {
        cpl_msg_warning (er_func, "Could not allocate initial %d bytes "
                         "for help description", SIZE_X);
        return;
    }

    /* store the indent, optional flag ("--") and keyword */

    i = (int) strlen(flag_prefix);
    ii = (int) strlen(keyword);
    if ((i+ii) > (SIZE_X-4))
    {
        cpl_msg_warning (er_func, "Size of keyword (= %d) too large - will "
                         "be truncated...", ii);
        mykeyword = (char *) cpl_malloc((size_t) SIZE_X);
        (void) strncpy(mykeyword,keyword,(SIZE_X-6));
        mykeyword[SIZE_X-6] = '\0';
    }
    else
        mykeyword = (char *) keyword;

    i = sprintf(mystr,"  %s%s",flag_prefix,mykeyword);
    if (tab < i) tab = i + 2;     /* if keyword quite long, move the tab */


    /* pad out to the tab-point, then append ": " */

    ptr = mystr + i;
    for (ii = i; ii < tab; ii++) *ptr++ = ' ';
    *ptr++ = ':';
    *ptr++ = ' ';
    *ptr = '\0';
    tab += 2;

    /* append description  */

    i = tab + (int) strlen(description) + 2;
    if (i > SIZE_X)
    {
        er_enlarge("paramutils_print_key_desc",&mystr,i);
    }
    (void) strcpy(mystr+tab, description);      /* append full description */


    /* finally, print the nicely wrapped string */

    (void) printf("%s\n",
                  er_strutils_split(mystr, tab, er_strutils_termwidth()));

    cpl_free(mystr);
} 


/**********************************************************************/
/**
 * @brief Manage buffer of sources of all used parameters (Esorex + recipes)
 *
 * @param flag       Control flag:
 *                   = 1, for given param. save source in internal buffer
 *                   = 2, for given param. pull out source from internal buffer
 *                   = 3, free allocated memory
 * @param name       Name of parameter
 * @param source     desription of source of a param.
 *           maybe "default", config-file name, "command line", etc
 *
 * This function manages a buffer of sources of all used parameters 
 * in Esorex + all called recipes.
 * Max MAX_SRC can be handled, i.e. for at most MAX_SRC different params.
 * Esorex has already 25 parameters...
 * Max length of parameter name is PAR_LEN chars.
 */
/**********************************************************************/

int
er_manage_sources(const int flag, const char *name, char **source)

#define MAX_SRC 60
#define PAR_LEN 40

{
    int offset, n, m, indx, lp;
    static int no_sources = -1;

    static char  *ptr, *p_sources[MAX_SRC];
    static char  p_names[MAX_SRC*PAR_LEN];

    size_t  slen;




    /* write parameter + source into save buffer */

    if (flag == 1)
    {
        if (no_sources < 0)             /* first time */
        {
            no_sources = 1;
            indx = 0;
            offset = 0;
        }
        else
        {                       /* loop through existing list */
            indx = -1;
            lp = PAR_LEN-1;
            for (n=0; n<no_sources; n++)
            {
                offset = n * PAR_LEN;
                m = strncmp(name,p_names+offset,lp);
                if (m == 0)
                {
                    indx = n;
                    cpl_free(p_sources[indx]); /* free space for new source */
                    break;
                }
            }
            if (indx == -1)  /* we need to add new entry */
            {
                if (no_sources < MAX_SRC)
                {
                    no_sources ++;
                    indx = no_sources - 1;
                    offset = indx * PAR_LEN;
                }
                else
                    return 1;               /* buffer full... */
            }
        }

        (void) strncpy(p_names+offset,name,PAR_LEN-1); /* save parameter name */
        *(p_names+offset+PAR_LEN-1) = '\0';

        slen = strlen(*source) + 1;      /* get memory for source */
        ptr = (char *) cpl_malloc(slen);
        if (ptr == NULL) return (-2);    /* no space? just leave... */

        (void) strcpy(ptr,*source);
        p_sources[indx] = ptr;
        return 0;
    }


    /* for given parameter read source out of save buffer */

    else if (flag == 2)
    {
        lp = PAR_LEN-1;
        for (n=0; n<no_sources; n++)
        {
            offset = n*PAR_LEN;
            m = strncmp(name,p_names+offset,lp);
            if (m == 0)
            {
                *source = p_sources[n];
                return 0;
            }
        }
        return (-1);                    /* par. name not found */
    }


    /* free all the allocated memory */

    else
    {
        for (n=0; n<no_sources; n++)
        {
            cpl_free(p_sources[n]);
        }
        no_sources = -1;
        return 0;
    }
}


/**********************************************************************/
/**
 * @brief Neatly print auxiliary information for a parameter
 *
 * @param flag_prefix This is the prefix used to flag certain types of
 *                    switches/parameters. For example a double hypen
 *                    "--" may be used before the long form of command
 *                    line options. It is only used for spacing.
 * @param keyword     This is the keyword itself. It is only used for 
 *                    spacing.
 * @param description This is the descriptive text that is to be 
 *                    printed. It is aligned to a tab stop to ensure 
 *                    that the output is neat.
 *
 * This function neatly prints an auxiliary line of information 
 * regarding a parameter. The infomation is printed to the normal output 
 * stream. The formatting is handled automatically to get neat alignment 
 * of the output.
 */
/**********************************************************************/

void
er_paramutils_print_aux_info(const char *flag_prefix, const char *keyword,
                             const char *description)

{
    int i, ii, tab = COMMENT_TAB_POSITION;

    char  spc[88];              /* to determine string padding */
    char  str[MAXSTRLENCONF];


    i = (int) strlen(flag_prefix);
    ii = (int) strlen(keyword);
    if ((i+ii) > 80)
    {
        cpl_msg_error (er_func, "Size of prefix (= %d), keyword (= %d) too "
                       "large", i, ii);
        cpl_error_set (er_func, CPL_ERROR_ILLEGAL_INPUT);
        return;
    }

    i = sprintf (spc, "  %s%s", flag_prefix, keyword);
    /* if the keyword is quite long, move the tab */
    if (tab < i) tab = i + 2;


    /* for the real string, pad out to the tab-point */

    memset ((void *) str, 32,(size_t) tab);     /* int 32 = ' ' */
    str[tab] = '\0';

    /* build up string with the description */
    (void) strcat (str, "  (source = ");
    i = (int) strlen(description);
    if (i > 900)
    {
        ii = (int) strlen(str);
        /* truncate description */
        (void) sprintf (&str[ii], "%900.900s", description);
    }
    else
        (void) strcat (str, description);
    (void) strcat (str, ")");

    /* finally, print the nicely wrapped string */

    (void) printf("%s\n",
                  er_strutils_split(str, (tab + 2), er_strutils_termwidth()));

}                      



/**********************************************************************/
/**
 * @brief  Pretty-print a parameter list with a given header text.
 *
 * @param  param_list   List of Parameters
 * @param  header       Header Text
 *
 * @returns Zero on sucess, or non-zero in the case of error.
 *
 * This function takes a parameter list, and a title and prints them.
 * The function makes use of the @c COMMENT_TAB_POSITION constant that
 * is defined in @c er_macros.h .
 */
/**********************************************************************/

int
er_paramutils_print_list (cpl_parameterlist *param_list, const char *header)

{
    cpl_parameter *p;
    cpl_msg_severity msg_level = CPL_MSG_ERROR;
    cpl_type ptype;

    const char *tmp_cmdl_keywd;
    const char *tmp_string;
    char tmp_cmdl_descr[MAXSTRLENCONF];

    int i, tmp_int;

    double tmp_double;


    tmp_cmdl_descr[0] = '\0';

    /* Check that the inputs are valid */

    if (param_list == NULL) return -1;
    if (header == NULL) return -1;

    /* Determine the message level (for printing auxiliary information) */
    msg_level = cpl_msg_get_level ();

    /* Print some header text */
    (void) printf ("%s :\n\n", header);

    /* Get the first parameter from the list */
    p = cpl_parameterlist_get_first (param_list);

    /* Loop through all the parameters in the list */
    while (p != NULL)
    {                   /* Get and print the keyword */
        tmp_cmdl_keywd = cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_CLI);

        /* Determine the type of the parameter */
        ptype = cpl_parameter_get_type (p);

        /* Depending on the parameter type, print the parameter appropriately */
        switch (ptype)
        {
            case CPL_TYPE_BOOL:
                (void) sprintf (tmp_cmdl_descr, "%s",
                                cpl_parameter_get_bool (p) ? "TRUE" : "FALSE");
                break;

            case CPL_TYPE_INT:
                tmp_int = cpl_parameter_get_int (p);
                (void) sprintf (tmp_cmdl_descr, "%d", tmp_int);
                break;

            case CPL_TYPE_DOUBLE:
                tmp_double = cpl_parameter_get_double (p);
                (void) sprintf (tmp_cmdl_descr, "%s",
                                er_strutils_dblstr(tmp_double));
                break;

            case CPL_TYPE_STRING:
                tmp_string = cpl_parameter_get_string (p);
                i = (int) strlen(tmp_string);
                if (i > MAXSTR)
                    (void) snprintf (tmp_cmdl_descr, MAXSTR, "%s",
                                     tmp_string);
                else
                    (void) sprintf (tmp_cmdl_descr, "%s",
                                    tmp_string ? tmp_string : "-");
                break;

            default:
                break;
        }

        /* Actually print the parameter pair */
        er_paramutils_print_key_desc ("", tmp_cmdl_keywd, tmp_cmdl_descr);

        /*
         * If we are set to DEBUG level, then print the source
         * of the value too
         */
        if (msg_level == CPL_MSG_DEBUG)
        {
            char  *myptr;

            i = er_manage_sources(2,cpl_parameter_get_name (p),&myptr);
            if (i == 0)
                er_paramutils_print_aux_info ("", tmp_cmdl_keywd, myptr);
            else
                er_paramutils_print_aux_info ("", tmp_cmdl_keywd,
                                              "ambiguous...");
        }

        /* Determine the next parameter in the list */
        p = cpl_parameterlist_get_next (param_list);
    }

    printf ("\n");

    return 0;

}                               /* End of er_paramutils_print_list() */



/**********************************************************************/
/**
 * @brief  Set the value of a parameter from that of a given string.
 *
 * @param p          Parameter that will have its value set.
 * @param value      A string containing the value to be assigned to the 
 *                   parameter.
 * @param source     The source of @c value
 *
 * @returns 0 if successfull, !=0 otherwise
 *
 * Function sets the current value of the given parameter. Converts
 * the string value to the necessary type if needed.
 */
/**********************************************************************/

int
paramutils_set_from_string(cpl_parameter *p, const char *value,
                           const char *source)

{
    cpl_parameter_class pclass;
    cpl_type ptype;

    int f_found = 0;    /* FALSE */
    int scanlen = 0, tmp_int = 0, i = 0;

    double tmp_double = 0.0;

    char **pptr;


    if (p == NULL)
    {
        cpl_error_set (er_func, CPL_ERROR_NULL_INPUT);
        return 1;
    }


    /* save for later where the value came from */

    pptr = (char **) &source;
    er_manage_sources(1,cpl_parameter_get_name(p),pptr);

    cpl_parameter_set_default_flag (p, TRUE);

    pclass = cpl_parameter_get_class (p);
    ptype = cpl_parameter_get_type (p);

    /*
     * ????? - This is a mess. The following switch blocks should be a
     * single switch block, calling a function fn_class(p, value); depending
     * on the class. Each of those functions would determine the type and then
     * do the actual processing.
     */

    switch (pclass)
    {
        case CPL_PARAMETER_CLASS_VALUE:
            switch (ptype)
            {
                case CPL_TYPE_BOOL:
                    if ((strcmp (value, "TRUE") == 0) ||
                            (strcmp (value, "true") == 0))
                    {
                        cpl_parameter_set_bool (p, TRUE);
                    }
                    else if ((strcmp (value, "FALSE") == 0) ||
                            (strcmp (value, "false") == 0))
                    {
                        cpl_parameter_set_bool (p, FALSE);
                    }
                    else
                    {
                        cpl_error_set (er_func, CPL_ERROR_TYPE_MISMATCH);
                        return -1;
                    }
                    break;

                case CPL_TYPE_INT:
                    scanlen = sscanf (value, " %d ", &tmp_int);
                    if (scanlen == 1)
                    {
                        cpl_parameter_set_int (p, tmp_int);
                    }
                    else
                    {
                        cpl_error_set (er_func, CPL_ERROR_TYPE_MISMATCH);
                        return -1;
                    }
                    break;

                case CPL_TYPE_DOUBLE:
                    scanlen = sscanf (value, " %le ", &tmp_double);
                    if (scanlen == 1)
                    {
                        cpl_parameter_set_double (p, tmp_double);
                    }
                    else
                    {
                        cpl_error_set (er_func, CPL_ERROR_TYPE_MISMATCH);
                        return -1;
                    }
                    break;

                case CPL_TYPE_STRING:
                    cpl_parameter_set_string (p, value);
                    break;

                default:
                    break;
            }

            break;


                case CPL_PARAMETER_CLASS_RANGE:

                    switch (ptype)
                    {
                        case CPL_TYPE_BOOL:
                            break;

                        case CPL_TYPE_INT:
                            scanlen = sscanf (value, " %d ", &tmp_int);
                            if (scanlen == 1)
                            {
                                if ((tmp_int <= cpl_parameter_get_range_max_int (p)) &&
                                        (tmp_int >= cpl_parameter_get_range_min_int (p)))
                                {
                                    cpl_parameter_set_int (p, tmp_int);
                                }
                                else
                                {
                                    cpl_error_set(er_func, CPL_ERROR_ILLEGAL_INPUT);
                                    cpl_msg_error(er_func,
                                                  "Specified 'int' value (%d) for the parameter "
                                                  "'%s' was outside the permitted range %d to %d\n",
                                                  tmp_int,
                                                  cpl_parameter_get_alias(p, CPL_PARAMETER_MODE_CLI),
                                                  cpl_parameter_get_range_min_int(p),
                                                  cpl_parameter_get_range_max_int(p));
                                    return -1;
                                }
                            }
                            else
                            {
                                cpl_error_set (er_func, CPL_ERROR_TYPE_MISMATCH);
                                return -1;
                            }
                            break;

                        case CPL_TYPE_DOUBLE:
                            scanlen = sscanf (value, " %le ", &tmp_double);
                            if (scanlen == 1)
                            {
                                if ((tmp_double <= cpl_parameter_get_range_max_double (p)) &&
                                        (tmp_double >= cpl_parameter_get_range_min_double (p)))
                                {
                                    cpl_parameter_set_double (p, tmp_double);
                                }
                                else
                                {
                                    char  *tmp;


                                    cpl_error_set (er_func, CPL_ERROR_ILLEGAL_INPUT);


                                    tmp = (char *) cpl_malloc((size_t) 1024);
                                    if ( tmp == NULL)
                                    {
                                        cpl_msg_error (er_func, "Could not allocate 1024 bytes for tmp");
                                        return -1;
                                    }

                                    /*
                                     * formulate the error message in steps...
                                     * er_strutils_dblstr() uses an internal static buffer!
                                     */

                                    (void) strcpy (tmp, "Specified 'double' value (");
                                    (void) strcat (tmp, er_strutils_dblstr (tmp_double));
                                    (void) strcat (tmp, ") for the parameter '");
                                    (void) strcat (tmp, cpl_parameter_get_alias (p, CPL_PARAMETER_MODE_CLI));
                                    (void) strcat (tmp, "' was outside the permitted range ");
                                    (void) strcat (tmp, er_strutils_dblstr (cpl_parameter_get_range_min_double (p)));
                                    (void) strcat (tmp, " to ");
                                    (void) strcat (tmp, er_strutils_dblstr (cpl_parameter_get_range_max_double (p)));
                                    (void) strcat (tmp, ".\n");
                                    cpl_msg_error (er_func, "%s", tmp);

                                    cpl_free(tmp);
                                    return -1;
                                }
                            }
                            else
                            {
                                cpl_error_set (er_func, CPL_ERROR_TYPE_MISMATCH);
                                return -1;
                            }
                            break;

                        case CPL_TYPE_STRING:
                            break;

                        default:
                            break;
                    }

                    break;


                        case CPL_PARAMETER_CLASS_ENUM:

                            switch (ptype)
                            {

                                case CPL_TYPE_BOOL:
                                    break;

                                case CPL_TYPE_INT:
                                    scanlen = sscanf (value, " %d ", &tmp_int);
                                    if (scanlen == 1)
                                    {
                                        for (i=0; i<cpl_parameter_get_enum_size (p); i++)
                                        {
                                            if (tmp_int == cpl_parameter_get_enum_int (p, i))
                                            {
                                                cpl_parameter_set_int (p, tmp_int);
                                                f_found = 1;
                                            }
                                        }
                                        if (!f_found)
                                        {
                                            cpl_error_set (er_func, CPL_ERROR_INCOMPATIBLE_INPUT);
                                            return -1;
                                        }
                                    }
                                    else
                                    {
                                        cpl_error_set (er_func, CPL_ERROR_TYPE_MISMATCH);
                                        return -1;
                                    }
                                    break;

                                case CPL_TYPE_DOUBLE:
                                    scanlen = sscanf (value, " %le ", &tmp_double);
                                    if (scanlen == 1)
                                    {

                                        /* WARNING --- The following code may need fabs() */

                                        for (i=0; i<cpl_parameter_get_enum_size (p); i++)
                                        {
                                            if ((cpl_parameter_get_enum_double (p, i) - tmp_double) < 0.000001)
                                            {
                                                cpl_parameter_set_double (p, tmp_double);
                                                f_found = 1;
                                            }
                                        }
                                        if (!f_found )
                                        {
                                            cpl_error_set (er_func, CPL_ERROR_INCOMPATIBLE_INPUT);
                                            return -1;
                                        }

                                    }
                                    else
                                    {
                                        cpl_error_set (er_func, CPL_ERROR_TYPE_MISMATCH);
                                        return -1;
                                    }
                                    break;

                                case CPL_TYPE_STRING:
                                    for (i=0; i<cpl_parameter_get_enum_size (p); i++)
                                    {
                                        if (strcmp (cpl_parameter_get_enum_string (p, i), value) == 0)
                                        {
                                            cpl_parameter_set_string (p, value);
                                            f_found = 1;
                                        }
                                    }
                                    if (! f_found )
                                    {
                                        cpl_error_set (er_func, CPL_ERROR_INCOMPATIBLE_INPUT);
                                        return -1;
                                    }

                                    break;

                                default:
                                    break;
                            }

                            break;


                                default:
                                    break;

    }                            /* End of switch(pclass) */

    return 0;

}
/**@}*/

/* End of file */
