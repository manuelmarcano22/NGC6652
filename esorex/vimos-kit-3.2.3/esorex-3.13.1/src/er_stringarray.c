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
#include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <cpl.h>

#include "er_macros.h"
#include "er_stringarray.h"

/**
 * @defgroup  er_stringarray  Array of Strings
 *
 * This module provides a collection of functions that create,
 * manipulate and destroy arrays of character strings. 
 *
 */

/**@{*/

/**
 * @struct er_stringarray_t
 * A simple array of strings
 */
struct er_stringarray_s
{
    char  **strings;
    int     nstrings;
};


/**********************************************************************/
/**
 * @brief   Creates a new er_stringarray
 *
 * @returns A newly allocated er_stringarray
 */
/**********************************************************************/

er_stringarray_t * er_stringarray_new(void)

{
    er_stringarray_t *csa;


    csa = (er_stringarray_t *) cpl_malloc((size_t)sizeof(er_stringarray_t));

    if (csa != NULL)
    {
        csa->strings  = NULL;
        csa->nstrings = 0;
    }

    return csa;

}  /* End of er_stringarray_new() */


/**********************************************************************/
/**
 * @brief
 *   Resizes a er_stringarray
 *
 * @param csa  Stringarray to resize
 * @param size New size of string array
 *
 * Resizes the stringarray so it can can contain @em size strings.
 * If the new size is equal to the old size nothing happens.
 * If the new size is less, all strings at the far side of the new size are
 * properly deleted and the er_stringarray is shrunk to the new size.
 * If the new size is larger than old size the array is grown to
 * the new size, all string present are kept.
 */
/**********************************************************************/

void er_stringarray_resize ( er_stringarray_t *csa, int size)

{
    int     i;

    char **tmp;


    if (csa==NULL) return;

    if (size==csa->nstrings) return;


    tmp = (char **) cpl_calloc((size_t)size, (size_t)sizeof(char *));

    if (size < csa->nstrings)
    {
        for (i=0; i<size; i++)
        {
            tmp[i] = csa->strings[i];
        }

        for (i=size; i<csa->nstrings; i++)
        {
            cpl_free(csa->strings[i]);
        }

    }
    else
    {
        for (i=0; i<csa->nstrings; i++)
        {
            tmp[i] = csa->strings[i];
        }

        for (i=csa->nstrings; i<size; i++)
        {
            tmp[i] = NULL;
        }
    }

    cpl_free(csa->strings);         /* get rid of old pointer array */

    csa->nstrings = size;           /* fill structure with new values */
    csa->strings  = tmp;

} 

/**********************************************************************/
/**
 * @brief
 *   Deletes a er_stringarray
 *
 * @param csa Stringarray to delete
 *
 * Deletes all information stored in a er_stringarray. All strings
 * contained in it are deleted. All memory associated with csa is properly
 * deallocated.
 */
/**********************************************************************/

void er_stringarray_delete ( er_stringarray_t *csa)

{
    int i;


    if (csa==NULL) return;

    if (csa->strings!=NULL)
    {
        for (i=0; i<csa->nstrings; i++)
        {
            if (csa->strings[i] != NULL ) cpl_free(csa->strings[i]);
        }
        cpl_free(csa->strings);
    }


    cpl_free(csa);

} 

/**********************************************************************/
/**
 * @brief
 *   Get a reference to the string at position index
 *
 * @param csa Stringarray to operate on
 * @param indx Position in Stringarray
 *
 * Returns A reference to the string stored at position @em indx in @em csa.
 */
/**********************************************************************/

char* er_stringarray_get ( er_stringarray_t *csa, int indx)

{

    if (csa == NULL) return NULL;
    if (csa->strings == NULL) return NULL;
    if ( (indx < 0) || (indx > csa->nstrings) ) return NULL;

    return csa->strings[indx];

}


/**********************************************************************/
/**
 * @brief
 *   Sets position index in stringarray to a copy of the string given
 *
 * @param csa Stringarray to operate on
 * @param instring String to store
 * @param indx Position in Stringarray
 *
 * Stores a copy of @em instring at position @em indx in @em csa.
 */
/**********************************************************************/

void er_stringarray_set (er_stringarray_t *csa, const char *instring, int indx)

{
    char  *cpt;


    if (csa == NULL) return;
    if (csa->strings == NULL) return;
    if ( (indx < 0) || (indx > csa->nstrings) ) return;

    if (csa->strings[indx] != NULL )
    {
        cpl_free(csa->strings[indx]);
    }

    cpt = cpl_malloc((size_t)((int)strlen(instring)+1));
    if (cpt != NULL) csa->strings[indx] = cpt;
    (void) strcpy(csa->strings[indx],instring);

} 



/**********************************************************************/
/**
 * @brief
 *   Adds a copy to given string at end of stringarray
 *
 * @param csa      Stringarray to operate on
 * @param instring String to store
 *
 * Stores a copy to @em instring at the last position in @em csa.
 */
/**********************************************************************/

void er_stringarray_append ( er_stringarray_t *csa, const char *instring)

{
    int new_size = 0;


    if (csa==NULL) return;

    new_size = csa->nstrings + 1;
    er_stringarray_resize(csa, new_size);

    er_stringarray_set(csa, instring, new_size-1);

}  


/**********************************************************************/
/**
 * @brief
 *   Removes a string at position index
 *
 * @param csa Stringarray to operate on
 * @param indx Position in stringarray
 *
 * Removes a @em string at position @em indx.
 * Effectively creates a 'hole' in the stringarray.
 */
/**********************************************************************/

void er_stringarray_remove ( er_stringarray_t *csa, int indx)

{

    if (csa == NULL) return;
    if (csa->strings == NULL) return;
    if ( (indx < 0) || (indx > csa->nstrings) ) return;

    if (csa->strings[indx] != NULL)
    {
        cpl_free(csa->strings[indx]);
        csa->strings[indx] = NULL;
    }

    return;

}



/**********************************************************************/
/**
 * @brief
 *   Checks whether there is a string present at position index
 *
 * @param csa Stringarray to operate on
 * @param indx Position in stringarray
 *
 * @returns 1 if there is, 0 if there is not.
 *
 */
/**********************************************************************/

int er_stringarray_present ( er_stringarray_t *csa, int indx)

{

    if (csa == NULL) return 0;
    if (csa->strings == NULL) return 0;
    if ((indx < 0) || (indx > csa->nstrings) ) return 0;

    if (csa->strings[indx] != NULL)
        return 1;

    return 0;

} 



/**********************************************************************/
/**
 * @brief    Returns current size aka capacity of stringarray.
 *
 * @param csa Stringarray to operate on
 *
 * @returns  Maximum number of strings the array can store
 *
 */
/**********************************************************************/

int er_stringarray_size ( er_stringarray_t *csa)

{

    if (csa == NULL) return 0;
    if (csa->strings == NULL) return 0;

    return(csa->nstrings);

}



/**********************************************************************/
/**
 * @brief     Prints a stringarray to standardformat
 *
 * @param csa Stringarray to operate on
 *
 * The format is "index : stringcontents".
 *
 */
/**********************************************************************/

void er_stringarray_print ( er_stringarray_t *csa)

{
    int i;


    if (csa == NULL) return;

    if (csa->strings != NULL)
    {
        for (i = 0; i < csa->nstrings; i++)
        {
            if (csa->strings[i] != NULL )
            {
                printf("%d : %s\n",i,csa->strings[i]);
            }
        }
    }

} 

/**@}*/

/* End of file */
