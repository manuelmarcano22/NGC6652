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

#include <dirent.h>
#include <fcntl.h>

#include <unistd.h>

#include <ctype.h>


#include <termios.h>
#include <sys/ioctl.h>
#include <stdlib.h>
#include <string.h>


#include <stdio.h>
#include <float.h>
#include <math.h>

#include <cpl.h>

#include "er_stringutils.h"

#define MAX_DBL_STR  16
#define COMMENT_TAB_POSITION   24



/**
 * @defgroup er_stringutils EsoRex String Utility Functions
 *
 * Utility functions for handling strings.
 *
 */

/**@{*/


/**********************************************************************/
/*
 * @brief         Indent a string.
 *
 * @param str     String to be processed.
 * @param indent  Number of characters to indent
 *
 * @return Pointer to the modified character string, or @c NULL in case of
 *   failure. The modified character string is statically allocated within
 *   this function and it cannot exceed 1024 characters.
 *
 */
/**********************************************************************/

char *er_strutils_indent (const char * str, int indent)

{
    static char output[MAX_MSG_LENGTH];
    const char *warning = "\n    ... [ TEXT TRUNCATED - MAXIMUM STRING LENGTH EXCEEDED ]\n";

    int c, i, s;

    memset (output, 0, MAX_MSG_LENGTH);

    c = 0;
    for (s=0; s<(int)strlen(str); s++)
    {
        if (c >= (MAX_MSG_LENGTH - 2 - indent))
        {
            (void) strcpy (&output[MAX_MSG_LENGTH - 2 - strlen (warning)], warning);
            /* Put a "safety-nul" on the end and return */
            output[MAX_MSG_LENGTH - 1] = '\0';
            return (output);
        }

        output[c++] = str[s];
        if (str[s] == '\n')
        {
            for (i = 0; i < indent; i++) output[c++] = ' ';
        }
    }

    /* Put a "safety-nul" on the end and return */
    output[MAX_MSG_LENGTH - 1] = '\0';
    return (output);

}                  

/**********************************************************************/
/*
 * @brief
 *   Split a string according to the max allowed page width.
 *
 * @param s       String to be processed.
 * @param blanks  Number of blanks to be inserted at every split point.
 * @param width   Max number of characters between split points.
 *
 * @return Pointer to the modified character string, or @c NULL in case of
 *   failure. The modified character string is statically allocated within
 *   this function and it cannot exceed 1024 characters.
 *
 * This private function is used for splitting a string avoiding to exceed
 * a maximum width (as for instance the width of the terminal where the 
 * string is going to be printed). The splitting is performed without 
 * breaking words, i.e. by replacing with a newline character ('\n') 
 * the last blank character before the maximum allowed width. Newline
 * characters already present in the input string are preserved.
 * Single words that exceed the max allowed width would not be split, 
 * just in this case long lines are tolerated. A number of blanks to 
 * be inserted at every split point must be specified, setting the 
 * left indentation level for the printed string. This number must 
 * not exceed the maximum allowed width.
 *
 * @note This function was originally a copy of the function @c strsplit from
 *       @c cpl_msg.c in the CPL (Derek).
 */
/**********************************************************************/

char *er_strutils_split (const char *ss, int blanks, int width)

{
    static char split[MAX_MSG_LENGTH];
    char *s, *savptr;

    int i, j, k;
    int cuti = 0;
    int cutj = 0;
    int limit;
    int nn, m;


    /* if the keyword name (alias) itself is large, we split it up as well and start
   description on the next line (the original split function failed on that) */

    if (blanks < (width/2))
    {
        s = (char *) ss;            /* take full string as it is */
    }

    else
    {
        char  csav;

        savptr = (char *) malloc((size_t)(blanks+2));
        (void) strncpy(savptr,ss,blanks);   /* isolate keyword name in savptr */
        *(savptr+blanks) = '\0';

        limit = width;
        m = 0, nn = blanks;         /* nn = length of keyword string */
        splitting:
        s = savptr+m;
        if (nn <= limit)
        {                   /* end of string reached */
            if (limit == width)
                (void) printf("%s\n",s);
            else
                (void) printf("  --%s\n",s);
        }
        else
        {
            for (i=limit; i>0; i--)
            {
                if (s[i] == ' ' || s[i] == '\t' || s[i] == '\n')
                {
                    s[i] = '\0';
                    if (m == 0)         /* 1st time */
                    {
                        (void) printf("%s\n",s);
                        limit = width - 4;  /* because of "  --" in the beginning */
                    }
                    else
                    {
                        (void) printf("  --%s\n",s);
                    }
                    m += (i+1); nn -= i;
                    goto splitting;
                }
            }
            csav = s[limit];            /* full line without spaces ... */
            s[limit] = '\0';
            if (m == 0)         /* 1st time */
            {
                (void) printf("%s\n",s);
                limit = width - 4;      /* because of "  --" in the beginning */
            }
            else
            {
                (void) printf("  --%s\n",s);
            }
            s[limit] = csav;
            m += limit; nn -= limit;
            goto splitting;
        }

        m = blanks;         /* offset description to default tab size */
        blanks = COMMENT_TAB_POSITION + 2;

        nn = strlen(ss)-m;      /* replace original name by new no. of blanks */
        s = (char *) malloc((size_t)(nn+blanks+2));
        for (i=0; i<blanks; i++) s[i] = ' ';
        (void) strcpy(s+blanks,ss+m);   /* and append original description */
    }

    /* now continue as original code */

    limit = width;
    for (i=0, j=0; i<MAX_MSG_LENGTH && j<MAX_MSG_LENGTH; i++, j++)
    {
        split[j] = s[i];

        if (s[i] == ' ' || s[i] == '\t' || s[i] == '\0' || s[i] == '\n')
        {
            if (i > limit)
            {
                /* Go back to the previous cuttable position, if possible */
                if (limit - cuti < width - blanks)
                {
                    j = cutj;
                    i = cuti;
                }
                else
                {
                    if (s[i] == '\0')
                        break;
                }

                /* Split here, and insert blanks */
                split[j] = '\n';

                for (k=0, j++; k<blanks && j<MAX_MSG_LENGTH; k++, j++)
                {
                    split[j] = ' ';
                }
                j--;

                limit = width - blanks + i;
            }
            else
            {
                if (s[i] == '\0')
                    break;

                if (s[i] == '\n')
                {
                    /* Split point already present in input string */
                    /* just add the secified number of blanks */

                    i++;
                    if (s[i] == '\0')
                    {
                        split[j] = '\0';
                        break;
                    }

                    for (k = 0, j++; k < blanks && j < MAX_MSG_LENGTH; k++, j++)
                    {
                        split[j] = ' ';
                    }
                    j--;

                    limit = width - blanks + i;
                }

                /* Keep track of the last cuttable position */

                cutj = j;
                cuti = i;

            }
        }
    }

    split[MAX_MSG_LENGTH - 1] = '\0';       /* safety belt */
    return split;

}                               /* End of er_strutils_split() */


/**********************************************************************/
/**
 * @brief   Determines the number of columns in the terminal
 *
 * @returns The number of columns in the terminal, or 80 if this cannot
 *          be determined.
 *
 * This function determines the number of columns that are present in
 * the terminal window at the time that it is called. If it is not 
 * possible to determine this, then a deafult of 80 characters is 
 * assumed.
 */
/**********************************************************************/

int er_strutils_termwidth (void)

{
    struct winsize win;

    int width;
    int fd = STDOUT_FILENO;                          /* File descriptor */


    /* Determine the width of the terminal (assume 80 if it cannot be found) */
    if (ioctl (fd, TIOCGWINSZ, &win) < 0 || win.ws_col < 1)
    {
        width = DEFAULT_TERM_WIDTH;
    }
    else
    {
        width = (int) win.ws_col;
    }

    return (width);

}                               /* End of er_strutils_termwidth() */


/**********************************************************************/
/**
 * @brief   Determines the number of rows in the terminal
 *
 * @returns The number of rows in the terminal, or 24 if this cannot
 *          be determined.
 *
 * This function determines the number of rows that are present in
 * the terminal window at the time that it is called. If it is not 
 * possible to determine this, then a deafult of 24 lines is assumed.
 */
/**********************************************************************/

int er_strutils_termheight (void)

{
    struct winsize win;

    int fd = STDOUT_FILENO;                          /* File descriptor */
    int height;


    /* Determine the width of the terminal (assume 23 if it cannot be found) */
    if (ioctl (fd, TIOCGWINSZ, &win) < 0 || win.ws_col < 1)
    {
        height = DEFAULT_TERM_WIDTH;
    }
    else
    {
        height = (int) win.ws_row;
    }

    return (height);

} /* End of er_strutils_termwidth() */


/**********************************************************************/
/**
 * @brief   Fills a string with a double precision value.
 *
 * @param   value A double precision number that is to be converted
 *          into a string.
 * @returns A pointer to a static buffer, containing the string
 *          representation of the given double.
 *
 * @warning This function uses a static buffer, and so should only be
 *          called once per operational element.
 *
 * This function is used to generate a string-representation of a
 * supplied double precision number, with the minimum amount of
 * precision in order to be able to accurately reconstruct the original
 * number. The string is stored locally as a static variable, and a pointer
 * to it is returned. Thus, the function must not be called multiple
 * times in a single operational element.
 */
/**********************************************************************/

const char *
er_strutils_dblstr (double value)

{
    const char *ptr;                /* Pointer to the output string */
    char  er_fmt[24];                       /* Printing format */
    static char  er_buf[256];                 /* Buffer to hold the created string */

    int prec;                                        /* Precision */
    int numread;                                     /* Number of values read using sscanf */

    double result;                                   /* The number written to the string */




    for (prec = 3; prec < MAX_DBL_STR; prec++)
    {
        snprintf (er_fmt,(size_t)20, "%c.%dg", '%', prec);
        snprintf (er_buf,(size_t)255, er_fmt, value);
        numread = sscanf (er_buf, "%lf", &result);
        if (numread != 1)
        {
            (void) strcpy (er_buf, "***");
            break;
        }

        /*
         * Admittedly, the next line is a direct comparison of two
         * double-precision numbers. However, this is actually what we want
         * to do, given that we are trying to match precision. And as the
         * number was originally read using the same system as that here, it
         * should generate the same internal representation.
         *
         * And "no", one can not use DBL_EPSILON. We are not testing
         * variance from 1.0.
         */

        if (result == value)
            break;

    }

    if (prec >= MAX_DBL_STR)
    {
        cpl_msg_warning ("strutils_dblstr", "Some precision may have been lost in the "
                         "conversion of the value '%s...'", er_buf);
    }

    /* Now check for the case that there is no decimal point */
    ptr = er_buf;
    if (strstr (ptr, ".") == NULL)
    {               /* Now check if it is scientific notation */
        if (strstr (ptr, "e") == NULL)
        {               /* The representation is decimal - simply append */
            (void) strcat (er_buf, ".0");
        }
    }

    return (er_buf);
}

/**@}*/

/* End of file */
