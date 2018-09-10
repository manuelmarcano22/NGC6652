/* $Id: pilalias.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
#include "pilalias.h"


/*
 * The alias object type
 */

struct _PIL_ALIAS_ {
    PilKeyword *keyword;
    char *format;
};

 void deletePilAlias(PilAlias *);


/**
 * @defgroup pilAlias pilAlias
 *
 * The module @b pilAlias implements an alias class.
 *
 * An alias consists of a name alias string, i.e. its name and the
 * translation of the alias string, i.e. its value. In addition to this
 * an alias must have format string used for printing and an optional
 * comment.
 */

/**@{*/

/**
 * @brief
 *   Create a new alias.
 *
 * @return The function returns a pointer to the allocated alias object if
 *   no error occurred, otherwise a @c NULL pointer is returned.
 *
 * @param name     Name alias.
 * @param value    Translated alias.
 * @param format   Format string.
 * @param comment  Comment string.
 *
 * The function allocates the memory for a new alias object. The name alias
 * string @em name, its associated mapping @em value and the format string 
 * @em format are used to initialize the structure and have to be given.
 * The alias description @em comment is optional. If @em comment is @c NULL
 * the alias description is left empty, i.e. set to @c NULL.
 * 
 * The strings passed for @em name, @em value, @em format and @em comment
 * have to be 0 terminated.
 * 
 * The keyword fields are filled by copying the passed input strings.
 */

 PilAlias *newPilAlias(const char *name, const char *value,
			     const char *format, const char *comment)
{

    PilAlias *self;


    assert(name != 0 && value !=0 && format != 0);

    if ((self = (PilAlias *)pil_malloc(sizeof *self))) {

	if (!(self->keyword = newPilKeyword(name, value, comment))) {
	    deletePilAlias(self);
	    return 0L;
	}

	if (!(self->format = pil_strdup(format))) {
	    deletePilAlias(self);
	    return 0L;
	}

    }

    return self;

}


/**
 * @brief
 *   Destroy an alias object.
 *
 * @return Nothing.
 *
 * @param alias  Alias object.
 *
 * The function first deallocates all previously allocated alias object
 * members and finally the alias object @em alias itself.
 * 
 * If @em alias is @c NULL, no operation is performed.
 */

 void deletePilAlias(PilAlias *alias)
{

    if (!alias)
	return;

    if (alias->keyword)
	deletePilKeyword(alias->keyword);

    if (alias->format)
	pil_free(alias->format);

    pil_free(alias);

    return;

}


/**
 * @brief
 *   Retrieve the name alias from an alias object.
 *
 * @return The function returns the name alias string on success and @c NULL
 *   if an error occurred.
 *
 * @param alias  Alias object.
 *
 * The function queries the alias object @em alias for the name alias
 * string. The alias object must exist.
 */

 const char *pilAliasGetName(PilAlias *alias)
{

    assert(alias != 0L);

    return pilKeyGetName(alias->keyword);

}


/**
 * @brief
 *   Set the name alias string of an alias object.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param alias  Alias object.
 * @param name   Name string.
 *
 * The function sets the name alias string @em name of the alias object
 * @em alias. The alias object must exist.
 */

 int pilAliasSetName(PilAlias *alias, const char *name)
{

    assert(alias != 0L && name != 0L);

    return pilKeySetName(alias->keyword, name);

}


/**
 * @brief
 *   Retrieve the name value from an alias object.
 *
 * @return The function returns the name value string on success and @c NULL
 *   if an error occurred.
 *
 * @param alias  Alias object.
 *
 * The function queries the alias object @em alias for the alias translation,
 * i.e. the name value string. The alias object must exist.
 */

 const char *pilAliasGetValue(PilAlias *alias)
{

    assert(alias != 0L);

    return pilKeyGetValue(alias->keyword);

}


/**
 * @brief
 *   Set the value string of an alias object.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param alias  Alias object.
 * @param value  Value string.
 *
 * The function sets the value string @em value of the alias object
 * @em alias. The alias object must exist.
 */

 int pilAliasSetValue(PilAlias *alias, const char *value)
{

    assert(alias != 0L && value != 0L);

    return pilKeySetValue(alias->keyword, value);

}


/**
 * @brief
 *   Retrieve the format string from an alias object.
 *
 * @return The function returns the format string on success and @c NULL
 *   if an error occurred.
 *
 * @param alias  Alias object.
 *
 * The function queries the alias object @em alias for the format string.
 * The alias object must exist.
 */

 const char *pilAliasGetFormat(PilAlias *alias)
{

    assert(alias != 0L);

    return alias->format;

}


/**
 * @brief
 *   Set the format string of an alias object.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param alias   Alias object.
 * @param format  Format string.
 *
 * The function sets the format string @em format of the alias object
 * @em alias. The alias object must exist.
 *
 * For details on the format string refer to the @b fprintf() documentation.
 */

 int pilAliasSetFormat(PilAlias *alias, const char *format)
{

    assert(alias != 0L && format != 0L);

    if (alias->format)
	pil_free(alias->format);

    if (!(alias->format = pil_strdup(format)))
	return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Retrieve the comment from an alias object.
 *
 * @return The function returns the comment string on success and @c NULL
 *   if an error occurred.
 *
 * @param alias  Alias object.
 *
 * The function queries the alias object @em alias for the comment string.
 * The alias object must exist.
 */

 const char *pilAliasGetComment(PilAlias *alias)
{

    assert(alias != 0L);

    return pilKeyGetComment(alias->keyword);

}


/**
 * @brief
 *   Set the comment string of an alias object.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param alias    Alias object.
 * @param comment  Comment string.
 *
 * The function sets the comment string @em comment of the alias object
 * @em alias. The alias object must exist.
 */

 int pilAliasSetComment(PilAlias *alias, const char *comment)
{

    assert(alias != 0L && comment != 0L);

    return pilKeySetComment(alias->keyword, comment);

}


/**
 * @brief
 *   Set the name, value, format and optionally the comment of an alias
 *   object.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param alias    Alias object.
 * @param name     Name string.
 * @param value    Value string.
 * @param format   Format string.
 * @param comment  Comment string or @c NULL.
 *
 * The function sets the @em name, @em value, @em format and, if it is
 * provided, the @em comment of the alias object @em alias.
 */

 int pilAliasSet(PilAlias *alias, const char *name, const char *value,
		       const char *format, const char *comment)
{

    assert(alias != 0L && name != 0L && value != 0L && format != 0L);

    if (pilKeySet(alias->keyword, name, value, comment) == EXIT_FAILURE)
	return EXIT_FAILURE;

    if (pilAliasSetFormat(alias, format) == EXIT_FAILURE)
	return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Clear the contents of an alias object.
 *
 * @return Nothing.
 *
 * @param alias  Alias object to be reset.
 *
 * The function resets the value of all members of the alias object 
 * @em alias to @c NULL.
 */

 void pilAliasClear(PilAlias *alias)
{

    assert(alias != 0L);

    pilKeyClear(alias->keyword);

    if (alias->format)
	pil_free(alias->format);

    alias->format = 0L;

    return;

}
/**@}*/
