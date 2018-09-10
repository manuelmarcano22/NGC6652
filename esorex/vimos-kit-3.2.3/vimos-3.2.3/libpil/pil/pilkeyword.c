/* $Id: pilkeyword.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
 *
 * This file is part of the VIMOS pipeline library
 * Copyright (C) 2000-2004 European Southern Observatory
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: cizzo $
 * $Date: 2008-10-21 09:10:13 $
 * $Revision: 1.1.1.1 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pilmemory.h"
#include "pilstrutils.h"
#include "pilkeyword.h"


/**
 * @defgroup pilKeyword pilKeyword
 *
 * The module @b pilKeyword implements a keyword class.
 *
 * A keyword has keyword name, an associated string value and additionally
 * may have a comment.
 */

/**@{*/

/*
 * The keyword object type
 */

struct _PIL_KEYWORD_ {
    char *name;
    char *value;
    char *comment;
};

 void deletePilKeyword(PilKeyword *);


/*
 * @brief
 *   Copy a string value to a keyword member.
 *
 * @param member  Keyword member of type pointer to @c char.
 * @param value   String value to be copied.
 *
 * @return The function returns @c EXIT_SUCCESS on success and
 *   @c EXIT_FAILURE otherwise.
 *
 * The function copies the string @em value into the keyword member
 * @em member. The keyword member must be a pointer to @c char. If
 * @em member in not @c NULL, i.e. a string was already assigned to
 * @em member, this string is deallocated before @em value is copied
 * and assigned to @em member. If @em value is @c NULL, member will
 * be set to @c NULL.
 */

 static int pilKeyCopyString(char **member, const char *value)
{

    if (*member)
	pil_free(*member);

    if (!value)
	*member = 0L;
    else
	if (!(*member = pil_strdup(value)))
	    return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Create an empty keyword.
 *
 * @return The function returns a pointer to the allocated keyword if no error
 *   occurred, otherwise a @c NULL pointer is returned.
 *
 * The function allocates memory for a keyword object. The keyword fields
 * are initialized with @c NULL.
 */

 PilKeyword *newPilKeywordEmpty(void)
{

    PilKeyword *key = (PilKeyword *)pil_malloc(sizeof *key);

    if (key) {
	key->name = 0L;
	key->value = 0L;
	key->comment = 0L;
    }

    return key;

}


/**
 * @brief
 *   Create a new keyword.
 *
 * @param name     Keyword identifier.
 * @param value    Keyword value string.
 * @param comment  Description of the keyword.
 *
 * @return The function returns a pointer to the allocated keyword if no error
 *   occurred, otherwise a @c NULL pointer is returned.
 *
 * The function allocates the memory for a new keyword object. The keyword
 * identifier @em name and its associated value string @em value
 * are used to initialize the structure and have to be given. The keyword
 * description @em comment is optional. If @em comment is @c NULL
 * the keyword description is left empty.
 * 
 * The strings passed for @em name, @em value and @em comment
 * have to be 0 terminated.
 * 
 * The keyword fields are filled by copying the passed input strings.
 */

 PilKeyword *newPilKeyword(const char *name, const char *value,
				 const char *comment)
{

    PilKeyword *self;


    assert(name != 0 && value !=0);

    if ((self = newPilKeywordEmpty())) {

	if (pilKeyCopyString(&self->name, name) == EXIT_FAILURE) {
	    deletePilKeyword(self);
	    return 0L;
	}

	if (pilKeyCopyString(&self->value, value) == EXIT_FAILURE) {
	    deletePilKeyword(self);
	    return 0L;
	}

	if (comment) {
	    if (pilKeyCopyString(&self->comment, comment) == EXIT_FAILURE) {
		deletePilKeyword(self);
		return 0L;
	    }
	}

    }

    return self;

}


/**
 * @brief
 *   Destroy a keyword object.
 *
 * @param keyword  Keyword object.
 *
 * @return Nothing.
 *
 * The function first deallocates the memory used for storing the
 * keyword identifier, its value string and the keyword description.
 * Finally the keyword object @em keyword itself is destroyed.
 * 
 * If @em keyword is @c NULL, no operation is performed.
 */

 void deletePilKeyword(PilKeyword *keyword)
{

    if (!keyword)
	return;

    if (keyword->name)
	pil_free(keyword->name);

    if (keyword->value)
	pil_free(keyword->value);

    if (keyword->comment)
	pil_free(keyword->comment);

    pil_free(keyword);

    return;

}


/**
 * @brief
 *   Retrieve the name of a keyword.
 *
 * @param keyword  Keyword structure to be examined.
 *
 * @return The function returns the name string of a keyword.
 *
 * The function retrieves the keyword identifier, i.e. its name, from
 * the keyword structure @em keyword. The keyword structure must exist.
 */

 const char *pilKeyGetName(PilKeyword *keyword)
{

    assert(keyword != 0L);

    return keyword->name;

}


/**
 * @brief
 *   Set the name of a keyword.
 *
 * @param keyword  Keyword structure to be modified.
 * @param name     Name string.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise @c EXIT_FAILURE is returned.
 *
 * The function stores the string @em name as the keyword identifier into
 * the keyword structure @em keyword. The string @em name is copied.
 */

 int pilKeySetName(PilKeyword *keyword, const char *name)
{

    assert(keyword != 0L && name != 0L);

    return pilKeyCopyString(&keyword->name, name);

}


/**
 * @brief
 *   Retrieve a keyword value.
 *
 * @param keyword  Keyword structure to be examined.
 *
 * @return The function returns the value string of a keyword.
 *
 * The function retrieves the keyword value string from the keyword
 * structure @em keyword. The keyword structure must exist.
 */

 const char *pilKeyGetValue(PilKeyword *keyword)
{

    assert(keyword != 0L);

    return keyword->value;

}


/**
 * @brief
 *   Set the value of a keyword.
 *
 * @param keyword  Keyword structure to be modified.
 * @param name     Value string.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise @c EXIT_FAILURE is returned.
 *
 * The function stores the string @em value as the keyword value into
 * the keyword structure @em keyword. The value string @em value is copied.
 */

 int pilKeySetValue(PilKeyword *keyword, const char *value)
{

    assert(keyword != 0L && value != 0L);

    return pilKeyCopyString(&keyword->value, value);

}


/**
 * @brief
 *   Retrieve the comment of a keyword.
 *
 * @param keyword  Keyword structure to be examined.
 *
 * @return The function returns the comment string of a keyword.
 *
 * The function retrieves the keyword comment from the keyword 
 * structure @em keyword. The keyword structure must exist.
 */

 const char *pilKeyGetComment(PilKeyword *keyword)
{

    assert(keyword != 0L);

    return keyword->comment;

}


/**
 * @brief
 *   Set the comment of a keyword.
 *
 * @param keyword  Keyword structure to be modified.
 * @param comment  Comment string.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise @c EXIT_FAILURE is returned.
 *
 * The function stores the string @em comment as the keyword comment into
 * the keyword structure @em keyword. The string @em comment is copied.
 */

 int pilKeySetComment(PilKeyword *keyword, const char *comment)
{

    assert(keyword != 0L && comment != 0L);

    return pilKeyCopyString(&keyword->comment, comment);

}


/**
 * @brief
 *   Set the name, value and optionally the comment of a keyword.
 *
 * @param keyword  Keyword structure to be modified.
 * @param name     Name string.
 * @param value    Value string.
 * @param comment  Comment string or @c NULL.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise @c EXIT_FAILURE is returned.
 *
 * The function copies the strings @em name, @em value and, if it is
 * provided, the @em comment into the appropriate members of the keyword
 * structure @em keyword.
 */

 int pilKeySet(PilKeyword *keyword, const char *name, const char *value,
		     const char *comment)
{

    assert(keyword != 0L && name != 0L && value != 0L);

    if (pilKeyCopyString(&keyword->name, name) == EXIT_FAILURE)
	return EXIT_FAILURE;

    if (pilKeyCopyString(&keyword->value, value) == EXIT_FAILURE)
	return EXIT_FAILURE;

    if (pilKeyCopyString(&keyword->comment, comment) == EXIT_FAILURE)
	return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Clear the contents of a keyword.
 *
 * @param keyword  Keyword to be reset.
 *
 * @return Nothing
 *
 * The function resets the value of all members of the keyword 
 * @em keyword to @c NULL.
 */

 void pilKeyClear(PilKeyword *keyword)
{

    assert(keyword != 0L);

    if (keyword->name)
	pil_free(keyword->name);

    if (keyword->value)
	pil_free(keyword->value);

    if (keyword->comment)
	pil_free(keyword->comment);

    keyword->name = 0L;
    keyword->value = 0L;
    keyword->comment = 0L;

    return;

}
/**@}*/
