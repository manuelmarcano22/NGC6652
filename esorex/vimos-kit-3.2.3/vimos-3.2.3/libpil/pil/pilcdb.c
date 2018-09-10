/* $Id: pilcdb.c,v 1.2 2013-08-07 16:25:37 cgarcia Exp $
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
 * $Author: cgarcia $
 * $Date: 2013-08-07 16:25:37 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <float.h>

#include "pilmemory.h"
#include "pildictionary.h"
#include "pilstrutils.h"
#include "pilutils.h"
#include "pilcdb.h"


#define PILCFG_GRP_DEFAULT    "<top>"
#define PILCFG_IFS_DEFAULT    '.'
#define PILCFG_COMMENT_CHARS  "#;"


/** 
 * @defgroup pilCdb pilCdb
 *
 * The module @b pilCdb provides functions to create a configuration
 * database from a file and to maintain it. The module also provides
 * methods to test the existence of database entries, to query the
 * database, retrieve and modify the stored information. There are no
 * other methods provided to create new database entries rather than
 * reading a configuration file.
 * 
 * The configuration database object keeps keyword/value pairs.
 * Regardless of the actual datatype of a value, the database stores
 * always its string representation. The value string is only interpreted
 * during retrieval through the dedicated methods.
 */

/**@{*/

/*
 * The database object type
 */

struct _PIL_CDB_ {
  char grp_ifs;
  PilCdbKeyCase key_case;
  PilDictionary *groups;
};

typedef struct _PIL_CDB_VALUE_ {
  char *value;
  PilCdbKeyMode read_only;
} PilCdbValue;
  
typedef PilDictNode PilCdbGroup;
typedef PilDictNode PilCdbEntry;

 static void deletePilCdbGroup(PilCdbGroup *, void *);
 static void deletePilCdbEntry(PilCdbEntry *, void *);
 static int pilCdbAddGroup(PilCdb *, const char *);
 PilCdbKeyCase pilCdbGetKeyCase(const PilCdb *);
 char pilCdbGetGroupIFS(const PilCdb *db);


 static int keycmp(const void *keyword1, const void *keyword2)
{
  return strcmp(keyword1, keyword2);
}


/*
 * @brief
 *   Build the configuration database internal keyword string.
 *
 * @return Pointer to the configuration database keyword string, or @c NULL
 *   in the case an error occurred.
 *
 * @param name   Configuration database keyword string.
 *
 * The function allocates memory for the fully qualified keyword name
 * used to uniquely identify a configuration database entry. The allocated
 * memory can be deallocated using pil_free(). It is a private function.
 */

 static char *db_create_name(PilCdb *db, const char *name)
{

  char *db_key;


  if (!name || strlen(name) == 0)
    return NULL;

  if (!(db_key = pil_strdup(name)))
    return NULL;
  else 
    if (pilCdbGetKeyCase(db) == IGNORE_CASE)
      strlower(db_key);

  return db_key;

}


/*
 * @brief
 *   Create a new configuration database value object.
 *
 * @return Pointer to the newly allocated database value object.
 *
 * @param value      The value string to be stored in the value object.
 * @param read_only  Flag indicating if the the configuration database value
 *                   is read only.
 *
 * Allocate the memory needed for a configuration database value object.
 * The new value object is initialized with the given value string and
 * read only flag. It is a private function.
 */

 static PilCdbValue *newPilCdbValue(const char *value, 
					  PilCdbKeyMode read_only)
{

  PilCdbValue *db_data = (PilCdbValue *)pil_malloc(sizeof *db_data);

  if (db_data) {
    if (value)
      db_data->value = pil_strdup(value);
    else
      db_data->value = NULL;

    db_data->read_only = read_only;
  }

  return db_data;

}


/*
 * @brief
 *   Destroy a configuration database value object.
 *
 * @return Nothing.
 *
 * @param db_data  Pointer to an existing configuration database value object.
 *
 * The function deallocates a configuration database value, making sure
 * that the value objects contents is properly deallocated before the
 * value object itself is destroyed. If @c NULL is passed to the function
 * nothing is done. It is a private function.
 */

 static void deletePilCdbValue(PilCdbValue *db_data)
{

  if (db_data) {
    if (db_data->value)
      pil_free(db_data->value);
    pil_free(db_data);
  }

  return;

}


/*
 * @brief
 *   Create a new configuration database entry.
 *
 * @return Pointer to the newly allocated database entry.
 *
 * Allocate the memory needed for a configuration database entry. The
 * memory needed to actually store the keyword string and the configuration
 * database value object is not allocated. It is a private function.
 */

 static PilCdbEntry *newPilCdbEntry(void *context)
{
  (void)context;
  return newPilDictNode(NULL);

}


/*
 * @brief
 *   Destroys a configuration database entry.
 *
 * @return Nothing.
 *
 * Destroys a configuration database entry including the entry's keyword
 * and value field, if the fields are used.
 */

 static void deletePilCdbEntry(PilCdbEntry *db_entry, void *context)
{

  (void)context;
  if (db_entry) {

    /*
     * Deallocate the database entry's keyword and value. The value's
     * destructor takes care that the value string is properly deallocated
     * before the entry itself is destroyed.
     */

    pil_free((void *)pilDictGetKey(db_entry));
    deletePilCdbValue(pilDictGetData(db_entry));


    /*
     * Destroy the database node
     */

    pil_free(db_entry);

  }

  return;

}


/*
 * @brief
 *   Create a new configuration database group.
 *
 * @return Pointer to the newly allocated database group.
 *
 * Allocate the memory needed for a configuration database group. Also
 * memory needed for the object that maintains the actual group entries
 * is created. It is a private function.
 */

 static PilCdbGroup *newPilCdbGroup(void *context)
{

  (void)context;
  PilCdbGroup *new = newPilDictNode(NULL);
  PilDictionary *entries;


  if (new) {
    if ((entries = newPilDictionary(PIL_DICT_CAPACITY_MAX, keycmp)) == NULL) {
      deletePilCdbGroup(new, NULL);
      return NULL;
    }
    else {
      pilDictSetAllocator(entries, newPilCdbEntry, deletePilCdbEntry, NULL);
      pilDictPutData(new, entries);
    }
  }

  return new;

}


/*
 * @brief
 *   Destroys a configuration database group.
 *
 * @return Nothing.
 *
 * Destroys a configuration database group including the group's keyword
 * and all group entries, if they are present.
 */

 static void deletePilCdbGroup(PilCdbGroup *db_group, void *context)
{

  (void)context;
  PilDictionary *entries;


  if (db_group) {

    /*
     * Deallocate the database group's keyword and value. The value
     * is another dictionary, i.e. if it is not already empty it has
     * to be cleaned up before it can be destroyed.
     */

    pil_free((void *)pilDictGetKey(db_group));

    entries = pilDictGetData(db_group);
    if (!pilDictIsEmpty(entries))
      pilDictClear(entries);
    deletePilDictionary(entries);


    /*
     * Destroy the database node
     */

    pil_free(db_group);

  }

  return;

}


/*
 * @brief
 *   Check if a database value is read only.
 *
 * @return The function returns @c READ_ONLY if the database value cannot
 *   be modified, otherwise @c READ_WRITE.
 *
 * @param db_value  Pointer to an existing configuration database value.
 *
 * The function checks if the read only flag is set for the given
 * database value. It is a private function.
 */

 static PilCdbKeyMode pilCdbGetMode(PilCdbValue *db_value)
{

  return db_value->read_only;

}


/*
 * @brief
 *   Set the database value access mode for the given keyword.
 *
 * @return Nothing.
 *
 * @param db_value  Pointer to an existing configuration database value.
 * @param mode      Database value access mode.
 *
 * The function sets the database values access mode according to 
 * @em mode. Allowed values are @c READ_WRITE or @c READ_ONLY.
 */

 static void pilCdbSetMode(PilCdbValue *db_value, PilCdbKeyMode mode)
{

  db_value->read_only = mode;
  return;

}


/*
 * @brief
 *   Set the value field of a database value object.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param db_value  Pointer to an existing database value object.
 * @param value     The string that is set as the value object's contents.
 *
 * The function replaces the contents of the database value object
 * by the given string, resizing the memory portion accordingly.
 * It is a private function.
 */

 static int pilCdbSetValue(PilCdbValue *db_value, const char *value,
				 PilCdbKeyMode mode)
{

  size_t sz = strlen(value);

  
  if (sz++ > strlen(db_value->value)) {
    db_value->value = (char *)pil_realloc(db_value->value, sz * sizeof(char));
    if (!db_value->value) {
      return EXIT_FAILURE;
    }
  }

  if (!memcpy(db_value->value, value, sz)) {
    return EXIT_FAILURE;
  }


  /*
   * Set keyword access mode
   */

  db_value->read_only = mode;

  return EXIT_SUCCESS;

}


/*
 * @brief
 *   Lookup the database group in the configuration database for a given group
 *   name.
 *
 * @return Pointer to the configuration database group or @c NULL if the group
 *   is not found.
 *
 * @param db     Pointer to an existing configuration database.
 * @param group  Name of the group where which should be searched for
 *               the configuration database keyword.
 *
 * The function searches for a configuration database group having the given
 * group name. If no such entry exists  or the database is empty a @c NULL
 * pointer is returned. It is a private function.
 */

 static PilCdbGroup *pilCdbFindGroup(PilCdb *db, const char *group)
{

  char *db_group_name;
  PilCdbGroup *db_group;

  if (pilDictIsEmpty(db->groups))
    return NULL;

  if (!(db_group_name = db_create_name(db, group)))
    return NULL;
  else 
    db_group = pilDictLookup(db->groups, db_group_name);

  pil_free(db_group_name);
  return db_group;

}


/*
 * @brief
 *   Create a database group in the configuration database for the given
 *   group name.
 *
 * @return The function returns @c EXIT_FAILURE if an error occurred,
 *   otherwise @c EXIT_SUCCESS is returned.
 *
 * @param db     Pointer to an existing configuration database.
 * @param group  Name of the group where which should be created
 *
 * The function creates a new group entry in the configuration database
 * provided by the caller. To be successful the database must not be
 * full and the group must either be a top level group, i.e. the group
 * name must not contain an internal field separator character, or if
 * the group is not a top level group the parent group must already exist
 * in the database. Adding an already existing group is an error. It is a
 * private function.
 */

 static int pilCdbAddGroup(PilCdb *db, const char *group)
{

  PilCdbGroup *db_group;
  char *db_group_name, *db_group_ifs, *db_parent_name;


  if (pilDictIsFull(db->groups))
    return EXIT_FAILURE;

  if (!(db_group_name = db_create_name(db, group)))
    return EXIT_FAILURE;

  if (pilDictLookup(db->groups, db_group_name)) {
    pil_free(db_group_name);
    return EXIT_FAILURE;
  }


  /*
   * Check if the group is a top level group. If not, i.e. the group name
   * contains at least one internal field separator character, the existence
   * of the parent group is checked and only if it is there the given group
   * is inserted. If the group simply is a top level group it is just
   * inserted.
   */

  if ((db_group_ifs = strrchr(db_group_name, pilCdbGetGroupIFS(db)))) {

    db_parent_name = pil_strdup(db_group_name);
    *(db_parent_name + (db_group_ifs - db_group_name)) = '\0';

    if (!pilDictLookup(db->groups, db_parent_name)) {
      pil_free(db_parent_name);
      pil_free(db_group_name);
      return EXIT_FAILURE;
    }

    pil_free(db_parent_name);

  }
      

  /*
   * Create a new group and insert it into the dictionary.
   */

  if (!(db_group = newPilCdbGroup(NULL))) {
    pil_free(db_group_name);
    return EXIT_FAILURE;
  }
  else 
    pilDictInsertNode(db->groups, db_group, db_group_name);

  return EXIT_SUCCESS;

}


/*
 * @brief
 *   Lookup the given database entry in the given configuration database
 *   group.
 *
 * @return Pointer to the configuration database entry or @c NULL if the
 *   entry is not found.
 *
 * @param db_group  Pointer to the database group to be searched for
 *                  the entry.
 * @param name      Keyword string of the entry to be looked up.
 *
 * The function looks for a group entry with the given keyword string. As
 * database group a pointer to a database group entry is required. Such
 * a pointer is, for instance, returned by pilCdbFindGroup(). It is
 * a private function.
 */

 static PilCdbEntry *pilCdbFindEntry(PilCdb *db, PilCdbGroup *db_group,
					   const char *name)
{

  char *db_entry_name;
  PilCdbEntry *db_entry;
  PilDictionary *entries = pilDictGetData(db_group);


  if (pilDictIsEmpty(entries))
    return NULL;

  if (!(db_entry_name = db_create_name(db, name)))
    return NULL;
  else 
    db_entry = pilDictLookup(entries, db_entry_name);

  pil_free(db_entry_name);
  return db_entry;

}


/*
 * @brief
 *   Lookup the database value for a given group and keyword in the
 *   configuration database.
 *
 * @return Pointer to the configuration database value or @c NULL if the entry
 *   is not found.
 *
 * @param db     Pointer to an existing configuration database.
 * @param group  Name of the group where which should be searched for
 *               the configuration database keyword.
 * @param name   Name of the configuration database keyword.
 *
 * The function searches for a configuration database entry belonging to the
 * given database group and having the given keyword name. It is a low level
 * function, that expects correct group and keyword name. If no such entry
 * exists  or the database is empty a @c NULL pointer is returned. It is a
 * private function.
 */

static PilCdbValue *pilCdbFindValue(PilCdb *db, const char *group,
				    const char *name)
{

  char *db_key;
  PilCdbEntry *db_entry;
  PilCdbGroup *db_group = pilCdbFindGroup(db, group);


  if (!db_group)
    return NULL;
  else {
    db_key = db_create_name(db, name);
    if (!(db_entry = pilDictLookup(pilDictGetData(db_group), db_key))) {
      pil_free(db_key);
      return NULL;
    }

    pil_free(db_key);
  }

  return (PilCdbValue *)pilDictGetData(db_entry);

}


/*
 * @brief
 *   Add a new paramter entry to the configuration database or modify
 *   an existing entry.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise @c EXIT_FAILURE.
 *
 * @param db     Pointer to an existing configuration database.
 * @param group  Configuration group to which the entry is added.
 * @param name   Keyword for the database entry.
 * @param value  Value string of the database entry.
 * @param mode   Indicator of keyword access mode.
 *
 * The function looks for a database entry in the given group having
 * the keyword name provided by the caller. If no such entry is found the
 * function creates a new database entry with the specified keyword name
 * within the given group, i.e. the group must already exist. The value of
 * the database entry is changed to the value given by the caller. An access
 * mode for the database entry must be specified, which is either
 * @c READ_WRITE or @c READ_ONLY. It is a private function.
 */

static int pilCdbAddEntry(PilCdb *db, const char *group, const char *name,
			  const char *value, PilCdbKeyMode mode)
{

  PilCdbValue *db_data;
  PilCdbEntry *db_entry;
  PilCdbGroup *db_group = pilCdbFindGroup(db, group);

  char *db_name;


  /*
   * Check if the database group exists.
   */

  if (!db_group)
    return EXIT_FAILURE;


  /*
   * Check if the database keyword is already in use.
   */

  if ((db_entry = pilCdbFindEntry(db, db_group, name))) {
    db_data = pilDictGetData(db_entry);
    if (pilCdbGetMode(db_data) == READ_ONLY)
      return EXIT_FAILURE;

    /*
     * The keyword already exists in the database and the value is
     * writable. Only the value string has to be modified.
     */

    pilCdbSetValue(db_data, value, mode);

  }
  else {

    /*
     * The keyword is not yet used. Create the new value object with the
     * given keyword mode and add it to the database.
     */

    if (!(db_name = db_create_name(db, name)))
      return EXIT_FAILURE;

    if ((db_data = newPilCdbValue(value, mode)) == NULL)
      return EXIT_FAILURE;
    else 
      if (!pilDictInsert(pilDictGetData(db_group), db_name, db_data)) {
	deletePilCdbValue(db_data);
	return EXIT_FAILURE;
      }
  }

  return EXIT_SUCCESS;

}


/*
 * @brief
 *   Dump the contents of a configuration group to a file.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param db_group  The group which should be dumped.
 * @param file      File pointer to the output stream.
 *
 * Walk through the list of entries in @em db_group and write the
 * entry name and its value to the output stream @em file.
 * 
 * If either @em db_group or @em file are equal to @c NULL the
 * function returns an error.
 */

 static int pilCdbDumpGroup(PilCdbGroup *db_group, FILE *file)
{
  PilCdbEntry *db_entry;
  PilCdbValue *db_data;
  char *db_key;
  PilDictionary *entries;


  if (!db_group || !file)
    return EXIT_FAILURE;


  /*
   * If there are no entries in this group nothing has to be written to
   * the output stream.
   */

  if (!(entries = pilDictGetData(db_group)) || pilDictIsEmpty(entries))
    return EXIT_FAILURE;

  db_entry = pilDictBegin(entries);
  while (db_entry) {
    db_key = (char *)pilDictGetKey(db_entry);
    db_data = pilDictGetData(db_entry);

    /*
     * If the keyword is read only also the prefix must be written.
     */
      
    if (pilCdbGetMode(db_data) == READ_ONLY)
      fprintf(file, "const ");

    /*
     * Use the correct representation for empty values and values
     * containing whitespace characters.
     */

    if (strempty(db_data->value, NULL))
      fprintf(file, "%s=\"\"\n", db_key);
    else {
      if (strchr(db_data->value, ' ') || strchr(db_data->value, '\t') ||
	  strchr(db_data->value, '\v') || strchr(db_data->value, '\n') ||
	  strchr(db_data->value, '\r') || strchr(db_data->value, '\f'))
	fprintf(file, "%s=\"%s\"\n", db_key, db_data->value);
      else 
	fprintf(file, "%s=%s\n", db_key, db_data->value);
    }
    db_entry = pilDictNext(entries, db_entry);
  }

  return EXIT_SUCCESS;

}


 static char **pilCdbDumpGrouptoString(PilCdbGroup *db_group, int *n)
{
  PilCdbEntry *db_entry;
  PilCdbValue *db_data;
  char *db_key;
  PilDictionary *entries;
  char **parString;
  int count = 0;


  *n = 0;

  if (!db_group)
    return NULL;


  /*
   * If there are no entries in this group nothing has to be written to
   * the output stream.
   */

  if (!(entries = pilDictGetData(db_group)) || pilDictIsEmpty(entries))
    return NULL;

  count = 0;
  db_entry = pilDictBegin(entries);
  while (db_entry) {
    db_entry = pilDictNext(entries, db_entry);
    count++;
  }

  parString = (char **)pil_calloc(count, sizeof(char *));

  count = 0;
  db_entry = pilDictBegin(entries);
  while (db_entry) {
    db_key = (char *)pilDictGetKey(db_entry);
    db_data = pilDictGetData(db_entry);

    /*
     * Use the correct representation for empty values and values
     * containing whitespace characters.
     */

    if (strempty(db_data->value, NULL)) {
      parString[count] = (char *)pil_calloc(strlen(db_key) + 5, sizeof(char));
      sprintf(parString[count], "%s=\"\"", db_key);
    }
    else {
      if (strchr(db_data->value, ' ') || strchr(db_data->value, '\t') ||
          strchr(db_data->value, '\v') || strchr(db_data->value, '\n') ||
          strchr(db_data->value, '\r') || strchr(db_data->value, '\f')) {
        parString[count] = (char *)pil_calloc(strlen(db_key) 
                                              + strlen(db_data->value) + 5, 
                                              sizeof(char));
        sprintf(parString[count], "%s=\"%s\"", db_key, db_data->value);
      }
      else {
        parString[count] = (char *)pil_calloc(strlen(db_key) 
                                              + strlen(db_data->value) + 2, 
                                              sizeof(char));
        sprintf(parString[count], "%s=%s", db_key, db_data->value);
      }
    }
    db_entry = pilDictNext(entries, db_entry);
    count++;
  }

  *n = count;

  return parString;

}


/**
 * @brief
 *   Creates a new configuration database.
 *
 * @return Pointer to the newly allocated configuration database.
 *
 * The function allocates and initializes a new, empty configuration
 * database. The new database is configured for a case sensitive
 * keyword comparison and uses the '.' character as internal field
 * separator for configuration group names.
 *
 * @see pilCdbSetKeyCase(), pilCdbSetGroupIFS()
 */

PilCdb *newPilCdb(void)
{

  PilCdb *new = (PilCdb *)pil_malloc(sizeof(PilCdb));


  if (new) {
    if (!(new->groups = newPilDictionary(PIL_DICT_CAPACITY_MAX, keycmp))) {
      pil_free(new);
      return NULL;
    }
    else {
      pilDictSetAllocator(new->groups, newPilCdbGroup, deletePilCdbGroup,
			  NULL);
      new->key_case = USE_CASE;
      new->grp_ifs = PILCFG_IFS_DEFAULT;
    }
  }

  return new;

}


/**
 * @brief
 *   Delete an existing configuration data base.
 *
 * @return Nothing.
 *
 * @param db  Pointer to an existing configuration database.
 *
 * The function deallocates the configuration database object and all
 * remaining database entries.
 */

void deletePilCdb(PilCdb *db)
{

  /*
   * If the configuration database is not empty, all database entries
   * are destroyed first. Then the dictionary is deallocated.
   */

  if (!pilDictIsEmpty(db->groups))
    pilDictClear(db->groups);
  deletePilDictionary(db->groups);

  pil_free(db);
  return;

}


/**
 * @brief
 *   Retrieve the group name internal field separator character from a
 *   database object.
 *
 * @return The function returns the internal field separator character.
 *
 * @param db    Pointer to an existing configuration database.
 *
 * The function retrieves the character from an object of type
 * PilCdb pointed to by @em db which is used as an internal
 * field separator for the group name strings and separates one subgroup
 * name from the other.
 *
 * @see pilCdbSetGroupIFS()
 */

 char pilCdbGetGroupIFS(const PilCdb *db)
{

  return db->grp_ifs;

}


/**
 * @brief
 *   Configure the group name internal field separator character.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise @c EXIT_FAILURE is returned.
 *
 * @param db    Pointer to an existing configuration database.
 * @param ifs   Internal field separator character to be used.
 *
 * The function stores the given character in the database object pointed
 * to by @em db. The character is used as an internal field separator
 * for the group names, i.e. to separate subgroups from each other.
 * To configure the database the database must be empty.
 *
 * @see pilCdbGetGroupIFS()
 */

 int pilCdbSetGroupIFS(PilCdb *db, char ifs)
{

  if (!pilDictIsEmpty(db->groups))
    return EXIT_FAILURE;
  else
    db->grp_ifs = ifs;

  return EXIT_SUCCESS;
}


/**
 * @brief
 *   Retrieve the keyword comparison mode the database is configured for.
 *
 * @return The function returns USE_CASE if the database is configured for
 *   case sensitive keyword comparison. If keywords should be compared 
 *   case insensitive the function returns IGNORE_CASE.
 *
 * @param db    Pointer to an existing configuration database.
 *
 * The function retrieves the the current mode for keyword comparisons,
 * which is either case sensitive or case insensitive. The setting is
 * stored in the database object pointed to by @em db.
 *
 * @see pilCdbSetKeyCase()
 */

 PilCdbKeyCase pilCdbGetKeyCase(const PilCdb *db)
{

  return db->key_case;

}


/**
 * @brief
 *   Configure the database to use the given mode when comparing keywords.
 *
 * @return The function returns @c EXIT_SUCCESS if the given mode could be
 *   set, otherwise it returns @c EXIT_FAILURE.
 *
 * @param db    Pointer to an existing configuration database.
 * @param mode  Keyword comparison mode.
 *
 * The function configures the configuration database @em db to use
 * the given mode when comparing database keywords. To be able to setup the
 * database correctly for the given mode the configuration database
 * must be empty. If this is not the case the function returns an error.
 *
 * @see newPilCdb()
 */

 int pilCdbSetKeyCase(PilCdb *db, PilCdbKeyCase mode)
{

  if (!pilDictIsEmpty(db->groups))
    return EXIT_FAILURE;
  else 
    db->key_case = mode;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Parse a configuration file and add its contents to the configuration
 *   database.
 *
 * @return The function returns @c EXIT_SUCCESS if the file could be parsed
 *   and its contents could be added to the database. If either the file
 *   could not be read or has a wrong format, or if the file entries could
 *   not be installed in the database the function returns @c EXIT_FAILURE.
 *
 * @param db       Pointer to an existing configuration database.
 * @param cfgfile  File pointer.
 *
 * The function checks if the file specified by @em filename exists
 * and is readable. If the file can be read the keyword-value pairs
 * it contains are added to the appropriate group of the configuration
 * database @em db.
 */

int pilCdbParseFile(PilCdb *db, FILE *cfgfile)
{

  char cfgline[PIL_LINE_LENGTH_MAX + 1];
  char grpname[PIL_LINE_LENGTH_MAX + 1];
  char keyname[PIL_LINE_LENGTH_MAX + 1], keyvalue[PIL_LINE_LENGTH_MAX + 1];
  char *lstart;

  PilCdbKeyMode keymode;



  /*
   * Check if the configuration database exists and the 
   * configuration file is given.
   */

  if (!db || !cfgfile)
    return EXIT_FAILURE;


  /*
   * Initially the built in top level group is used
   */

  strcpy(grpname, PILCFG_GRP_DEFAULT);


   /*
    * Read configuration file ignoring empty lines and comment lines.
    */

  clearerr(cfgfile);

  while (fgets(cfgline, PIL_LINE_LENGTH_MAX, cfgfile)) {
    if (!strempty(cfgline, PILCFG_COMMENT_CHARS)) {
      lstart = strskip(cfgline, isspace);

      if (sscanf(lstart, "[%[^]]", grpname)) {

	/*
	 * A group name was found. If it is already in use there is
	 * nothing to do. Otherwise the group is created.
	 */

	if (!pilCdbFindGroup(db, grpname))
	  if (pilCdbAddGroup(db, grpname) == EXIT_FAILURE)
	    return EXIT_FAILURE;
      }
      else {

	/*
	 * If the current keyword-value pair does not belong to a
	 * particular group it is added to the default group. If
	 * needed the default group has to be created. This is done
	 * here.
	 */

	if (!strcmp(grpname, PILCFG_GRP_DEFAULT))
	  if (!pilCdbFindGroup(db, grpname))
	    if (pilCdbAddGroup(db, grpname) == EXIT_FAILURE)
	      return EXIT_FAILURE;

	/*
	 * Here we handle the keyword-value pairs. If the keyword string
	 * starts with 'const' the keyword is added as READ_ONLY, otherwise
	 * it is added as READ_WRITE. The value might be enclosed in either
	 * double or single quotes. If the value is not quoted parsing
	 * stops at any of the comment characters if present.
	 */

	if (sscanf(lstart, "const %[^=] = \"%[^\"]\"", keyname, keyvalue) == 2
	    || sscanf(lstart, "const %[^=] = '%[^']'", keyname, keyvalue) == 2
	    || sscanf(lstart, "const %[^=] = %[^" PILCFG_COMMENT_CHARS "]'",
		      keyname, keyvalue) == 2)
	  keymode = READ_ONLY;
	else 
	  if (sscanf(lstart, "%[^=] = \"%[^\"]\"", keyname, keyvalue) == 2
	      || sscanf(lstart, "%[^=] = '%[^']'", keyname, keyvalue) == 2 
	      || sscanf(lstart, "%[^=] = %[^" PILCFG_COMMENT_CHARS "]'",
			keyname, keyvalue) == 2)
	    keymode = READ_WRITE;
	  else 
	    return EXIT_FAILURE;

	/*
	 * Leading and trailing whitespace characters are removed from
	 * the keyword name and the keyword value. If the value is two
	 * successive double or single quotes the value is handled as
	 * an empty string value.
	 */

	strtrim(keyname, 2);
	strtrim(keyvalue, 2);

	if (!strcmp(keyvalue, "\"\"") || !strcmp(keyvalue, "''"))
	  keyvalue[0] = '\0';
	  

	/*
	 * Add the entry to the database
	 */

	if (pilCdbAddEntry(db, grpname, keyname, keyvalue,
			   keymode) == EXIT_FAILURE)
	  return EXIT_FAILURE;
      }
    }
  }


  /*
   * Check if reading the configuration file did not stop due to an error
   */

  if (!feof(cfgfile) || ferror(cfgfile))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Dump the contents of a configuration data base to a file.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param db       Pointer to an existing configuration database.
 * @param cfgfile  File pointer.
 *
 * The function write the contents of the configuration database @em db
 * to the stream pointed to by @em cfgfile. The output format is the
 * same as for input configuration files.
 * 
 * The function loops the group dictionary. For each entry the function
 * writes the group name to the output file followed by all parameters
 * defined in that group.
 */

int pilCdbDumpDB(PilCdb *db, FILE *cfgfile)
{

  PilCdbGroup *db_group;


  /*
   * If there is no database, no output file or if the group dictionary
   * is empty nothing is done.
   */

  if (!db || !cfgfile || pilDictIsEmpty(db->groups))
    return EXIT_FAILURE;


  /*
   * The default group is special. In the output a configuration file
   * all entries appearing before the first group header are added to
   * the default group. Therefore the contents of the default group
   * must be dumped first.
   */

  if ((db_group = pilCdbFindGroup(db, PILCFG_GRP_DEFAULT))) {
    if (pilCdbDumpGroup(db_group, cfgfile) == EXIT_FAILURE)
      return EXIT_FAILURE;
    else 
      fprintf(cfgfile, "\n");
  }


  /*
   * The remaining database groups must be dumped. Walk through the
   * group dictionary and write all groups which are not the default
   * group to the output stream. For all groups the group header must
   * be printed.
   */

  db_group = pilDictBegin(db->groups);
  while (db_group) {
    if (strcmp(pilDictGetKey(db_group), PILCFG_GRP_DEFAULT)) {
      fprintf(cfgfile, "[%s]\n", (char *)pilDictGetKey(db_group));

      if (pilCdbDumpGroup(db_group, cfgfile) == EXIT_FAILURE)
	return EXIT_FAILURE;
      else 
	if (db_group != pilDictEnd(db->groups))
	  fprintf(cfgfile, "\n");
    }
    db_group = pilDictNext(db->groups, db_group);
  }

  return EXIT_SUCCESS;

}


char **pilCdbDumpDBtoString(PilCdb *db, int *n)
{

  PilCdbGroup *db_group;

  *n = 0;

  if ((db_group = pilCdbFindGroup(db, "Parameters")))
    return pilCdbDumpGrouptoString(db_group, n);
  else
    return NULL;
}


/**
 * @brief
 *   Check if a configuration group exists.
 *
 * @return The function returns 1 if the group exists, otherwise 0 is
 *   returned.
 *
 * @param db     Pointer to the database to be queried.
 * @param group  Name of the configuration group.
 *
 * The function searches the configuration database pointed to by
 * @em db for the given configuration group name. If it is found
 * the group exists.
 */

int pilCdbGroupExists(const PilCdb *db, const char *group)
{

  if (pilCdbFindGroup((PilCdb *)db, group))
    return 1;

  return 0;

}


/**
 * @brief
 *   Check if a configuration parameter exists.
 *
 * @return The function returns 1 if the parameter exists, otherwise 0 is
 *   returned.
 *
 * @param db     Pointer to the database to be queried.
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 *
 * The function searches the configuration database pointed to by
 * @em db for the given configuration group name. If the group
 * exists, this group is searched for the given configuration parameter
 * name.
 */

int pilCdbEntryExists(const PilCdb *db, const char *group,
			 const char *name)
{

  if (pilCdbFindValue((PilCdb *)db, group, name))
    return 1;

  return 0;

}


/**
 * @brief
 *   Retrieve the access mode of the value associated to a database entry.
 *
 * @return The function returns the database values access mode.
 *
 * @param db     Pointer to the database to be queried.
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 *
 * The function locates the database entry and retrieves the current
 * setting of the values access mode flag. The function does not check
 * if the requested database entry exists or not.
 */

PilCdbKeyMode pilCdbGetKeyMode(const PilCdb *db, const char *group,
			       const char *name)
{

  return pilCdbGetMode(pilCdbFindValue((PilCdb *)db, group, name));

}


/**
 * @brief
 *   Set the access mode of the value associated to a database entry.
 *
 * @return The function returns @c EXIT_SUCCESS if setting the access mode
 *   was successful, otherwise @c EXIT_FAILURE is returned.
 *
 * @param db     Pointer to the database to be queried.
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 * @param mode   Access mode setting to be used for the value.
 *
 * The function locates the database entry and sets the access mode flag
 * of the associated value to @em mode. If the database entry does
 * not exist the function fails.
 */

int pilCdbSetKeyMode(PilCdb *db, const char *group, const char *name,
		     PilCdbKeyMode mode)
{

  PilCdbValue *db_value = pilCdbFindValue((PilCdb *)db, group, name);

  if (!db_value)
    return EXIT_FAILURE;

  pilCdbSetMode(db_value, mode);

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Add a new group  to the database
 *
 * @return The function returns @c EXIT_SUCCESS if the group was created
 *   successfully, otherwise @c EXIT_FAILURE is returned.
 *
 * @param db     Pointer to the database to be updated.
 * @param group  Name of the configuration group.
 *
 * The function creates the configuration group @em group in the
 * database @em db. The group is created if it is a top level
 * group. If it is a subgroup its parent group must exist. Otherwise
 * the function fails.
 */

int pilCdbCreateGroup(PilCdb *db, const char *group)
{

  if (!db)
    return EXIT_FAILURE;


  /*
   * Create the database group.
   */

  return pilCdbAddGroup(db, group);
  
}


/**
 * @brief
 *   Add a new database entry to a database.
 *
 * @return The function returns @c EXIT_SUCCESS if the entry was created
 *   successfully, otherwise @c EXIT_FAILURE is returned.
 *
 * @param db     Pointer to the database to be updated.
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 * @param value  Value string of the parameter.
 *
 * The function creates a new configuration parameter @em name within
 * the configuration group @em group of the database @em db. If
 * the group does not exist it is created as long as it is a top level group
 * or its parent group already exists. The parameter's value field will
 * contain a copy of the value string @em value and the access mode
 * is set to @c READ_WRITE.
 * 
 * If the database contains already an entry @em name within the group
 * @em group the function will fail. To modify an existing entry
 * the function @b pilCdbModifyValue() must be used.
 * 
 * Note that the group and entry identifiers which are actually inserted
 * into the database are duplicates of the strings @em group and
 * @em name. This applies also to the entry's value string
 * @em value.
 *
 * @see pilCdbModifyValue()
 */

int pilCdbCreateEntry(PilCdb *db, const char *group, const char *name,
		     const char *value)
{

  if (!db || pilCdbFindValue(db, group, name))
    return EXIT_FAILURE;


  /*
   * Check if the database group exists. If it is not yet there it is
   * created.
   */

  if (!pilCdbFindGroup(db, group))
    if (pilCdbAddGroup(db, group) == EXIT_FAILURE)
      return EXIT_FAILURE;


  /*
   * Add the parameter entry
   */

  if (pilCdbAddEntry(db, group, name, value, READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get a configuration database value as string.
 *
 * @return The value string.
 *
 * @param db     Pointer to the database to be queried.
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 *
 * The database pointed to by @em db is queried for the parameter
 * @em name stored in the configuration group @em group. If
 * the parameter is found its value is returned as a string.
 * 
 * Actually the function returns the pointer to the value string which
 * is stored in the database. Therefore the value string referenced by
 * the returned pointer must not be modified through this pointer directly.
 * To modify the value string in the database use the interface function
 * @b pilCdbModifyValue().
 *
 * @see pilCdbModifyValue()
 */

const char *pilCdbGetString(const PilCdb *db, const char *group,
			    const char *name)
{

  PilCdbValue *db_value = pilCdbFindValue((PilCdb *)db, group, name);

  if (!db_value)
    return NULL;

  return db_value->value;

}


/**
 * @brief
 *   Get a configuration database value as boolean.
 *
 * @return The function returns 1 if the value is true and 0 if the value
 *   is false. If an error occurs the value of the argument errval is
 *   returned.
 *
 * @param db      Pointer to the database to be queried.
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * The database pointed to by @em db is queried for the parameter
 * @em name stored in the configuration group @em group. If
 * the parameter is found its value is returned as a boolean. The string
 * stored in the database is interpreted as a boolean value, if it matches
 * any of the strings given in the following lists.
 * 
 * The strings "true", "T" and "1" are interpreted as true. The
 * strings "false", "F" and "0" are interpreted as false.
 * 
 * If the configuration group or the parameter can not be found, or if
 * the string stored in the database cannot be converted, or any other
 * error occurs, the function returns a user provided value. This value
 * has to be passed as argument @em errval when calling the function.
 *
 * @see pilCdbModifyValue()
 */

int pilCdbGetBool(const PilCdb *db, const char *group, const char *name,
		  int errval)
{

  PilCdbValue *db_value = pilCdbFindValue((PilCdb *)db, group, name);

  if (db_value) {
    if (!strncmp(db_value->value, "true", 4) || *db_value->value == 'T' ||
	*db_value->value == '1')
      return 1;
    else 
      if (!strncmp(db_value->value, "false", 5) || *db_value->value == 'F' ||
	  *db_value->value == '0')
	return 0;
  }

  return errval;

}


/**
 * @brief
 *   Get a configuration database value as integer.
 *
 * @return The numerical value of the database value string as int. If an
 *   error occurs the value of the argument errval is returned.
 *
 * @param db      Pointer to the database to be queried.
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * The database pointed to by @em db is queried for the parameter
 * @em name stored in the configuration group @em group. If
 * the parameter is found its value is returned as an integer value.
 * To be recognized as an integer the value must fulfil the relation
 * @f$ @c INT_MIN <= value <= @c INT_MAX @f$.
 * 
 * If the configuration group or the parameter can not be found, or if
 * the string stored in the database cannot be converted, or any other
 * error occurs, the function returns a user provided value. This value
 * has to be passed as argument @em errval when calling the function.
 *
 * @see pilCdbModifyValue()
 */

int pilCdbGetInt(const PilCdb *db, const char *group, const char *name,
		 int errval)
{

  PilCdbValue *db_value = pilCdbFindValue((PilCdb *)db, group, name);
  char *last;
  long val;


  if (db_value) {
    val = strtol(db_value->value, &last, 10);

    /*
     * The full string must have been converted and the value must be
     * a valid integer number. No range violation may have happened.
     */

    if (*last == '\0')
      if (INT_MIN <= val && val <= INT_MAX)
	return (int)val;
  }

  return errval;

}


/**
 * @brief
 *   Get a configuration database value as long.
 *
 * @return The numerical value of the database value string as a long. If
 *   an error occurs the value of the argument errval is returned.
 *
 * @param db      Pointer to the database to be queried.
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * The database pointed to by @em db is queried for the parameter
 * @em name stored in the configuration group @em group. If
 * the parameter is found its value is returned as a long value. To be
 * recognized as a valid long the value must fulfil the relation
 * @f$ @c LONG_MIN <= value <= @c LONG_MAX @f$.
 * 
 * If the configuration group or the parameter can not be found, or if
 * the string stored in the database cannot be converted, or any other
 * error occurs, the function returns a user provided value. This value
 * has to be passed as argument @em errval when calling the function.
 *
 * @see pilCdbModifyValue()
 */

long pilCdbGetLong(const PilCdb *db, const char *group, const char *name,
		   long errval)
{

  PilCdbValue *db_value = pilCdbFindValue((PilCdb *)db, group, name);
  char *last;
  long val;


  if (db_value) {
    errno = 0;
    val = strtol(db_value->value, &last, 10);

    /*
     * The full string must have been converted and the value must be
     * a valid integer number. No range violation may have happened.
     */

    if (*last == '\0' && errno == 0)
	return val;
  }

  return errval;

}


/**
 * @brief
 *   Get a configuration database value as float.
 *
 * @return The numerical value of the database value string as a float. If
 *   an error occurs the value of the argument errval is returned.
 *
 * @param db      Pointer to the database to be queried.
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * The database pointed to by @em db is queried for the parameter
 * @em name stored in the configuration group @em group. If
 * the parameter is found its value is returned as a float value. To be
 * recognized as a valid float the value must fulfil the relation
 * @f$ @c FLT_MIN <= value <= @c FLT_MAX @f$.
 * 
 * If the configuration group or the parameter can not be found, or if
 * the string stored in the database cannot be converted, or any other
 * error occurs, the function returns a user provided value. This value
 * has to be passed as argument @em errval when calling the function.
 *
 * @see pilCdbModifyValue()
 */

float pilCdbGetFloat(const PilCdb *db, const char *group, const char *name,
		     float errval)
{

  PilCdbValue *db_value = pilCdbFindValue((PilCdb *)db, group, name);
  char *last;
  double val;


  if (db_value) {
    val = strtod(db_value->value, &last);

    /*
     * The full string must have been converted and the value must be
     * a valid float number. No range violation may have happened.
     */

    if (*last == '\0')
      if (FLT_MIN <= val && val <= FLT_MAX)
	return (float)val;
  }

  return errval;

}


/**
 * @brief
 *   Get a configuration database value as double.
 *
 * @return The numerical value of the database value string as a double. If
 *   an error occurs the value of the argument errval is returned.
 *
 * @param db      Pointer to the database to be queried.
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * The database pointed to by @em db is queried for the parameter
 * @em name stored in the configuration group @em group. If
 * the parameter is found its value is returned as a double value. To be
 * recognized as a valid double the value must fulfil the relation
 * @f$ @c DBL_MIN <= value <= @c DBL_MAX @f$.
 * 
 * If the configuration group or the parameter can not be found, or if
 * the string stored in the database cannot be converted, or any other
 * error occurs, the function returns a user provided value. This value
 * has to be passed as argument @em errval when calling the function.
 *
 * @see pilCdbModifyValue()
 */

double pilCdbGetDouble(const PilCdb *db, const char *group, const char *name,
		       double errval)
{

  PilCdbValue *db_value = pilCdbFindValue((PilCdb *)db,	group, name);
  char *last;
  double val;


  if (db_value) {
    errno = 0;
    val = strtod(db_value->value, &last);

    /*
     * The full string must have been converted and the value must be
     * a valid float number. No range violation may have happened.
     */

    if (*last == '\0' && errno == 0)
	return val;
  }

  return errval;

}


/**
 * @brief
 *   Modify an existing value string of a database parameter entry.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param db     Pointer to an existing database.
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 * @param value  Value string of the parameter.
 *
 * The function updates the value string of an already existing parameter
 * entry, given by @em group and @em name, in the configuration
 * database pointed to by @em db. The function cannot be used to
 * create new database entries. If the parameter entry does not
 * exist, i.e. either the group or the parameter does not exist, the
 * function will fail.
 * 
 * The function can be used for any kind of the retrievable datatypes.
 * For the numerical types the string representation of the numerical
 * value has to be passed. 
 *
 * @see pilCdbGetString(), pilCdbGetBool(), pilCdbGetInt(),
 *   pilCdbGetLong(), pilCdbGetFloat(), pilCdbGetDouble()
 */

int pilCdbModifyValue(PilCdb *db, const char *group, const char *name,
			 const char *value)
{


  if (!db || !pilCdbFindValue(db, group, name))
    return EXIT_FAILURE;

  if (pilCdbAddEntry(db, group, name, value, READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}
/**@}*/
