/* $Id: pildfsconfig.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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


#include <pwd.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include "pilmemory.h"
#include "pilcdb.h"
#include "pilstrutils.h"
#include "pilfileutils.h"
#include "pilutils.h"
#include "pildfsconfig.h"


#define PIL_ROOTDIR         "PIPE_HOME"
#define PIL_EXPORT_FLAG     "DFS_PIPE_ALLPRODUCTS"
#define PIL_OVERWRITE_FLAG  "DFS_PIPE_OVERWRITE"
#define PIL_DISPLAY         "DFS_PIPE_DISPLAY"
#define DFS_LOGDIR          "DFS_LOG"
#define DFS_PRODIR          "DFS_PRODUCT"
#define DFS_RBDIR           "DFS_REDBLOCK"
#define DFS_EXPORTDIR       "DFS_DATA_REDUCED_OLAS"
#define PIL_CFGFILESUFFIX   ".ini"
#define PIL_CFGFILEPATH     "config"
#define PIL_USERPREFIX      ".pipeline"


/**
 * @defgroup pilDfsConfig pilDfsConfig
 *
 * The module @b pilDfsConfig provides functions to create the recipe
 * configuration database. This is an internal database structure used
 * for storing dataflow system and pipeline specific information retrieved
 * from the runtime environment and configuration files when the internal
 * structure is created. 
 * 
 * The module also implements methods to retrieve data from the recipe
 * configuration database, to check the existance of database entries,
 * to modify existing entries, to load the contents of a configuration file
 * into the internal structure and to query the runtime environment.
 */

/**@{*/

/*
 * This is the recipe configuration data base. It is hidden to the outside
 * world and can only be accessed by the interface defined in this module.
 */

static PilCdb *cfgDB = NULL;


/*
 * @brief
 *   Get a normalized pathname from an environment variable.
 *
 * @param env       Source environment variable.
 *
 * @return The function returns the string containing the path if no error
 *   occurred, otherwise it returns @c NULL.
 *
 * The function creates a copy of the string value of @em env and
 * replaces a, possibly present, trailing slash by a '\0'.
 */

 static char *getpath(const char *env)
{

  register char *path = getenv(env);

  if (!path)
    return NULL;

  return pilFileTrimPath(pil_strdup(path));

}


/*
 * @brief
 *   Build the full configuration filename from the format string
 *   and the argument list.
 *
 * @param fmt  Format string used for printing to the output string
 * 
 * @return The function returns a string, where the filename is stored, if
 *   no error occurred, otherwise it returns @c NULL.
 *
 * The function uses the format string and the variable-length argument
 * facilities to create the full filename of a pipeline configuration
 * file. The full filename is printed into the internally created output
 * string.
 *
 * The memory allocated for the returned string can be deallocated using
 * @b pil_free().
 */

 static char *rcfilepath(const char *fmt, ...)
{

  size_t sz;
  char path[PIL_PATHNAME_MAX + 1];
  va_list ap;

  va_start(ap, fmt);

#ifdef HAVE_VSNPRINTF
  sz = vsnprintf(path, PIL_PATHNAME_MAX, fmt, ap);
#else
  sz = vsprintf(path, fmt, ap);
#endif

  va_end(ap);

  if (sz != strlen(path))
    return NULL;

  return pil_strdup(path);

}


/*
 * @brief
 *   Retrieve the root directory of the pipeline installation from the
 *   runtime environment.
 *
 * @param instrument  Instrument name.
 *
 * @return The function returns a string, where the pathname is stored, if
 *   no error occurred, otherwise it returns @c NULL.
 *
 * The function checks if the environment variable PIPE_HOME exists. If it
 * exists it is checked whether the stored pathname contains the instrument
 * name given by @em instrument. If the instrument name is the last
 * part of the path, the path is copied to the output string which is
 * created internally. If @em instrument is not already part of the
 * pathname the root directory name is build by concatenating the pathname
 * and @em instrument and copying the result into the, internally
 * created, output string.
 *
 * The memory allocated for the returned string can be deallocated using
 * @b pil_free().
 */

 static char *pilDfsGetRootDir(const char *instrument)
{

  register char *s;
  register char *troot;

  char *root;
  size_t sz;
  size_t offset;


  if (!instrument || !(troot = getenv(PIL_ROOTDIR)))
    return NULL;
  
  root = pilFileTrimPath(pil_strdup(troot));

  s = strstr(root, instrument);
  if (!s || *(s + strlen(instrument)) != '\0') {
    sz = strlen(root);
    offset = sz;

    /* Add an extra byte for the slash and memory for the instrument name */
    sz += strlen(instrument) + 2;
    
    if ((root = (char *)pil_realloc(root, sz * sizeof(char))) == NULL) {
      pil_free(root);
      return NULL;
    }
    else {
      s = root + offset;
      *s++ = '/';
      memcpy(s, instrument, strlen(instrument) + 1);
    }
  }

  return root;

}


/*
 * @brief
 *   Retrieve the user's home directory from the runtime environment.
 *
 * @return The function returns a string containing the user's home directory
 *   if no error occurred, otherwise it returns @c NULL.
 *
 * Retrieves the user's home directory from the system and copies the path
 * name to the, internally created, output string.
 *
 * The memory allocated for the returned string can be deallocated using
 * @b pil_free().
 */

 static char *pilDfsGetHomeDir(void)
{

  register char *path;


#if defined HAVE_GETUID && defined HAVE_GETPWUID
  struct passwd *pw;

  if (!(pw = getpwuid(getuid())))
    return NULL;
  else
    path = pw->pw_dir;
#else
  if ((path = getenv("HOME")) == NULL)
    return NULL;
#endif

  return pilFileTrimPath(pil_strdup(path));

}


/*
 * @brief
 *   Retrieve the path to the dataflow system logfile directory from the
 *   runtime environment.
 *
 * @return The function returns a string containing the pathname if no error
 *   occurred, otherwise it returns @c NULL.
 *
 * The function retrieves the path of a dataflow system specific 
 * directory from the environment. This directory specifies the location
 * where log files must be placed. The path is copied to the, internally
 * created, output string.
 *
 * The memory allocated for the returned string can be deallocated using
 * @b pil_free().
 */

 static char *pilDfsGetLogDir(void)
{

  return getpath(DFS_LOGDIR);

}
  

/*
 * @brief
 *   Retrieve the path to the pipeline product directory from the
 *   runtime environment.
 *
 * @return The function returns a string containing the pathname if no error
 *   occurred, otherwise it returns @c NULL.
 *
 * The function retrieves the path of a dataflow system specific 
 * directory from the environment. This directory specifies the location
 * where the pipeline products must be placed. The path is copied to the,
 * internally, created output string.
 *
 * The memory allocated for the returned string can be deallocated using
 * @b pil_free().
 */

 static char *pilDfsGetProductDir(void)
{

  return getpath(DFS_PRODIR);

}


/*
 * @brief
 *   Retrieve the path to the online archive system link directory from the
 *   runtime environment.
 *
 * @return The function returns a string containing the pathname if no error
 *   occurred, otherwise it returns @c NULL.
 *
 * The function retrieves the path of a dataflow system specific 
 * directory from the environment. This directory specifies the location
 * where the symbolic links to the pipeline products will be placed for all
 * products being exported to the archive. The path is copied to the,
 * internally created, output string.
 *
 * The memory allocated for the returned string can be deallocated using
 * @b pil_free().
 */

 static char *pilDfsGetExportDir(void)
{

  return getpath(DFS_EXPORTDIR);

}


/*
 * @brief
 *   Retrieve the pipeline product export configuration from the runtime
 *   environment.
 *
 * @return The function returns a string indicating the current configuration.
 *
 * The function queries the runtime environment for the variable
 * DFS_PIPE_ALLPRODUCTS. If the environment variable is set to 'NO' a
 * symbolic link is created only for the main data product. If it is set to
 * 'YES' or any other value a link is created for all data products. The
 * returned, internally created, string indicates the current configuration
 * using the strings "MainOnly" and "AllProducts" respectively.
 *
 * The memory allocated for the returned string can be deallocated using
 * @b pil_free().
 */

 static char *pilDfsGetExportFlag(void)
{

  char *value, *flag;
  register char *tvalue = getenv(PIL_EXPORT_FLAG);


  if (!tvalue)
    flag = NULL;
  else {
    value = strlower(pil_strdup(tvalue));

    if (strncmp(value, "no", 2) != 0)
      flag = pil_strdup("AllProducts");
    else 
      flag = pil_strdup("MainOnly");

    pil_free(value);
  }

  return flag;

}

  
/*
 * @brief
 *   Retrieve the configuration setting for overwriting existing pipeline
 *   products.
 *
 * @return The function returns the string "true" if products may be
 *   overwritten, otherwise the string "false" is returned.
 *
 * The function queries the runtime environment for the variable
 * DFS_PIPE_OVERWRITE. If its value is set to 'YES' the pipeline will
 * overwrite already existing pipeline products in the pipeline product
 * directory if necessary. If the variable is set to 'NO' or does
 * not exist at all, the pipeline will not overwrite any file in the
 * output directory. The returned, internally created, string indicates
 * the current configuration. The returned string can be used as a
 * boolean value in the configuration database.
 *
 * The memory allocated for the returned string can be deallocated using
 * @b pil_free().
 */

 static char *pilDfsGetOverwriteFlag(void)
{

  char *value, *flag;
  register char *tvalue = getenv(PIL_OVERWRITE_FLAG);


  if (!tvalue)
    flag = NULL;
  else {
    value = strlower(pil_strdup(tvalue));

    if (strncmp(value, "yes", 3) != 0)
      flag = pil_strdup("true");
    else
      flag = pil_strdup("false");

    pil_free(value);
  }

  return flag;

}



/*
 * @brief
 *   Fill the recipe configuration database with initial defaults.
 *
 * @return The function returns @c EXIT_SUCCESS if the default values could be
 *   inserted, otherwise the return value is @c EXIT_FAILURE.
 *
 * The function inserts the minimum set of recipe configuration defaults
 * into the recipe configuration database. This set is needed by
 * the higher level functions. Only dataflow system specific defaults
 * are treated here.
 *
 * The set of database entries inserted into the group @b DfsConfig, are:
 *   @li PipelineMode=Online
 *   @li AllowUserConfiguration=true
 *   @li LogDir=cwd()
 *   @li Verbosity=Info
 *   @li LogLevel=Off
 *   @li ProductDir=cwd()
 *   @li ProductPrefix=recipe()
 *   @li OverwriteProducts=false
 *   @li CopyProducts=false
 *   @li ExportDir=cwd()
 *   @li ExportProducts=NoExport
 *
 * The set of database entries inserted into the group @b Visualization,
 * are:
 *   @li EnableDisplays=false
 *   @li EnableGraphics=false
 */

 static int pilDfsInitDB(void)
{


  /*
   * Defaults for the database group 'DfsConfig'
   */

  if (pilDfsDbCreateEntry("DfsConfig", "PipelineMode", "Online", 
			  READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "AllowUserConfiguration", "true",
			  READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "LogDir", ".", READ_WRITE) ==
      EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "Verbosity", "Info", READ_WRITE) ==
      EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "LogLevel", "Off", READ_WRITE) ==
      EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "ProductDir", ".", READ_WRITE) ==
      EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "ProductPrefix", "recipe()",
			  READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "OverwriteProducts", "false",
			  READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "CopyProducts", "false", 
			  READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "ExportDir", ".", READ_WRITE) ==
      EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("DfsConfig", "ExportProducts", "NoExport",
			  READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;


  /*
   * Defaults for the database group 'Visualization'
   */

  if (pilDfsDbCreateEntry("Visualization", "EnableDisplays", "false",
			  READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (pilDfsDbCreateEntry("Visualization", "EnableGraphics", "false",
			  READ_WRITE) == EXIT_FAILURE)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/*
 * @brief
 *   Retrieve a dataflow specific environment variable and update the
 *   corresponding database entry.
 *
 * @param agent  Function retrieving the environment variables value.
 * @param group  Database group name.
 * @param name   Database entry name.
 *
 * @return The function returns 1 if the database could be updated,
 *   otherwise the return value is 0.
 *
 * The function calls @em agent to retrieve the value of the
 * environment variable from the system. The function @em agent
 * must return the variables value string, or, if the variable is
 * not defined, @c NULL. The database entry specified through @em group
 * and @em name is updated with the retrieved value string. If
 * @em agent returns @c NULL the database entry is left untouched.
 */

 static int pilDfsUpdateFromEnv(char *(*agent)(void), const char *group,
				      const char *name)
{

  register char *value;

  if ((value = agent())) {
    if (pilCdbModifyValue(cfgDB, group, name, value) == EXIT_FAILURE) {
      pil_free(value);
      return 0;
    }
    else 
      pil_free(value);
  }

  return 1;

}


/**
 * @brief
 *   Setup the recipe configuration database.
 *
 * @param grp_ifs   Group name internal field separator.
 * @param key_case  Keyword case sensitivity flag.
 *
 * @return The function returns @c EXIT_SUCCESS if the database could be
 *   created, if an error occures @c EXIT_FAILURE is returned.
 *
 * A new recipe configuration database is created with the specified
 * attributes. The attributes allow for a user specified field separator
 * character which is used to separate subgroup names from each other and to
 * configure the case sensitivity of the database.
 *
 * Allowed  values for the separator character @em grp_ifs are all
 * printable characters except characters belonging to [a-z,A-Z,0-9,_].
 * To use the built in default, i.e. the '.' character, 0 has to be passed
 * when calling the function.
 *
 * To make the database queries case sensitive @em key_case must be
 * @c USE_CASE. If @em key_case is @c IGNORE_CASE the database queries
 * will be case insensitive.
 *
 * Finally the recipe configuration database is filled with the default
 * entries.
 */

int pilDfsCreateDB(int grp_ifs, PilCdbKeyCase key_case)
{


  /*
   * Do not create a new database if there is already one.
   */

  if (cfgDB)
    return EXIT_FAILURE;
  else { 

    /*
     * Create a new database
     */

    if ((cfgDB = newPilCdb()) == NULL)
      return EXIT_FAILURE;

    pilCdbSetKeyCase(cfgDB, key_case);


    /*
     * The internal field separator is configured only if grp_ifs is not 0.
     * In this case only, characters which do not belong to [a-z,A-Z,0-9,_]
     * are allowed as internal field separators.
     */

    if (grp_ifs != 0) {
      
      /*
       * Note: The characer class functions have to be enclosed in
       *   parantheses to use the actual function on HP-UX where these
       *   functions are also provided as macros, which are taken by
       *   default and may lead to compiler warnings.
       */

      if ((isspace)(grp_ifs) || !(ispunct)(grp_ifs)) {
	deletePilCdb(cfgDB);
	return EXIT_FAILURE;
      }

      if (pilCdbSetGroupIFS(cfgDB, (char)grp_ifs) == EXIT_FAILURE) {
	deletePilCdb(cfgDB);
	return EXIT_FAILURE;
      }
    }


    /*
     * Initialize the database with internal defaults
     */

    if (pilDfsInitDB() == EXIT_FAILURE) {
      deletePilCdb(cfgDB);
      return EXIT_FAILURE;
    }
  }

  
  return EXIT_SUCCESS;
    
}


/**
 * @brief
 *   Destroy the recipe configuration database.
 *
 * @return Nothing.
 *
 * Public destructor of the recipe configuration database.
 */

void pilDfsFreeDB(void)
{

  if (cfgDB) {
    deletePilCdb(cfgDB);
    cfgDB = NULL;
  }
  
  return;

}


/**
 * @brief
 *   Dump the contents of the recipe configuration database to a file.
 *
 * @param filename  Path and name of the file to which the database is
 *                  written.
 *
 * @return The function returns @c EXIT_SUCCESS if writing the database was
 *   successful, otherwise the return value is @c EXIT_FAILURE.
 *
 * The given file is created and and the current contents of the recipe
 * configuration database is dumped to the file the same format of the
 * recipe configuration files. The dumped file can then be used as
 * a recipe configuration file. If @em filename is set to @c NULL, the
 * contents of the configuration database is dumped to @b stdout.
 */

int pilDfsDumpDB(const char *filename)
{

  int status = EXIT_SUCCESS;

  FILE *cfgfile = stdout;


  if (filename && strlen(filename) > 0) {

    if (!(cfgfile = fopen(filename, "w")))
      return EXIT_FAILURE;
  }
   
  if (pilCdbDumpDB(cfgDB, cfgfile) == EXIT_FAILURE || ferror(cfgfile))
    status = EXIT_FAILURE;

  if (cfgfile != stdout)
    fclose(cfgfile);

  return status;

}

char **pilDfsDumpDBtoString(int *n)
{

  return pilCdbDumpDBtoString(cfgDB, n);

}


/**
 * @brief
 *   Read the pipeline configuration files and initialize the recipe
 *   configuration database.
 *
 * @param instrument  Instrument to be configured
 * @param recipe      Recipe to be configured
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred when
 *   parsing the configuration files or updating the database. If an error
 *   occurred the return value is @c EXIT_FAILURE.
 */

int pilDfsReadSetupFiles(const char *instrument, const char * recipe)
{

  char *prefix, *filename;
  FILE *cfgfile;


  /*
   * Do not overwrite an existing database
   */

  if (!instrument || !recipe || !cfgDB)
    return EXIT_FAILURE;


  /*
   * Check if the root directory of the pipeline installation is defined.
   */

  if (!(prefix = pilDfsGetRootDir(instrument)))
    return EXIT_FAILURE;


  /*
   * Parse the instrument specific configuration file.
   */

  filename = rcfilepath("%s/%s/%s%s", prefix, PIL_CFGFILEPATH, instrument,
		     PIL_CFGFILESUFFIX);
  if (!filename) {
    pil_free(prefix);
    return EXIT_FAILURE;
  }
  else {
    if (!(cfgfile = fopen(filename, "r"))) {      
      pil_free(filename);
      pil_free(prefix);
      return EXIT_FAILURE;
    }
    else {
      if (pilCdbParseFile(cfgDB, cfgfile) == EXIT_FAILURE) {
	fclose(cfgfile);
	pil_free(filename);
	pil_free(prefix);
	return EXIT_FAILURE;
      }
      else {
	fclose(cfgfile);
        pil_free(filename);
      }
    }
  }


  /*
   * Parse the system configuration file for the recipe
   */

  filename = rcfilepath("%s/%s/%s%s", prefix, PIL_CFGFILEPATH, recipe,
		     PIL_CFGFILESUFFIX);
  if (!filename) {
    pil_free(prefix);
    return EXIT_FAILURE;
  }
  else {
    if (!(cfgfile = fopen(filename, "r"))) {
      pil_free(filename);
      pil_free(prefix);
      return EXIT_FAILURE;
    }
    else {
      if (pilCdbParseFile(cfgDB, cfgfile) == EXIT_FAILURE) {
	fclose(cfgfile);
	pil_free(filename);
	pil_free(prefix);
	return EXIT_FAILURE;
      }
      else {
	pil_free(filename);
	pil_free(prefix);
      }
    }
  }


  /*
   * If user specific configuration files exist, try to parse them if
   * the user is allowed to.
   */

  if (pilCdbGetBool(cfgDB, "DfsConfig", "AllowUserConfiguration", 0))
    if ((prefix = pilDfsGetHomeDir())) {
      filename = rcfilepath("%s/%s/%s/%s%s", prefix, PIL_USERPREFIX,
			    instrument, instrument, PIL_CFGFILESUFFIX);
      if (filename && (cfgfile = fopen(filename, "r"))) {
	pilCdbParseFile(cfgDB, cfgfile);
	fclose(cfgfile);
	pil_free(filename);
      }

      filename = rcfilepath("%s/%s/%s/%s%s", prefix, PIL_USERPREFIX, 
			    instrument, recipe, PIL_CFGFILESUFFIX);
      if (filename && (cfgfile = fopen(filename, "r"))) {
	pilCdbParseFile(cfgDB, cfgfile);
	fclose(cfgfile);
	pil_free(filename);
      }

      pil_free(prefix);
    }

  
  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Retrieve the pipeline specific configuration of the dataflow system
 *   from the runtime environment and updates the recipe configuration
 *   database.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise it returns @c EXIT_FAILURE;
 *
 * The function queries the applications runtime environment for the 
 * following environment variables:
 *   @li  DFS_PIPE_ALLPRODUCTS,
 *   @li  DFS_PIPE_OVERWRITE,
 *   @li  DFS_PIPE_DISPLAY (not yet implemented),
 *   @li  DFS_LOG,
 *   @li  DFS_PRODUCT, and
 *   @li  DFS_DATA_REDUCED_OLAS
 * and assigns the values of the defined variables to their corresponding
 * entries in the configuration database. If a variable is not set the
 * database entry is not modified.
 *
 * The function requires that the database contains the configuration 
 * group @b DfsConfig and that the database entries
 *   @li  ExportProducts,
 *   @li  OverwriteProducts,
 *   @li  Display (not yet implemented),
 *   @li  LogDir,
 *   @li  ProductDir, and
 *   @li  ExportDir,
 * corresponding to the environment variables, are already existent in the
 * database. If the configuration group or any of the database entries is
 * missing the function will return an error. 
 */

int pilDfsGetEnv(void)
{

  if (!pilDfsUpdateFromEnv(pilDfsGetLogDir, "DfsConfig", "LogDir"))
    return EXIT_FAILURE;

  if (!pilDfsUpdateFromEnv(pilDfsGetProductDir, "DfsConfig", "ProductDir"))
    return EXIT_FAILURE;

  if (!pilDfsUpdateFromEnv(pilDfsGetExportDir, "DfsConfig", "ExportDir"))
    return EXIT_FAILURE;

  if (!pilDfsUpdateFromEnv(pilDfsGetExportFlag, "DfsConfig",
			   "ExportProducts"))
    return EXIT_FAILURE;

  if (!pilDfsUpdateFromEnv(pilDfsGetOverwriteFlag, "DfsConfig",
			   "OverwriteProducts"))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Check if a configuration group exists in the recipe configuration
 *   database.
 *
 * @param group  Name of the configuration group.
 *
 * @return The function returns 1 if the group exists, otherwise 0 is
 *   returned.
 *
 * The function searches the recipe configuration database for the given
 * configuration group name. If it is found the group exists.
 */

int pilDfsDbGroupExists(const char *group)
{

  return pilCdbGroupExists(cfgDB, group);

}


/**
 * @brief
 *   Check if a parameter exists in the recipe configuration database.
 *
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 *
 * @return The function returns 1 if the parameter exists, otherwise 0 is
 *   returned.
 *
 * The function searches the recipe configuration database for the
 * given configuration group name. If the group exists, this group is
 * searched for the given parameter name.
 */

int pilDfsDbEntryExists(const char *group, const char *name)
{

  return pilCdbEntryExists(cfgDB, group, name);

}


/**
 * @brief
 *   Retrieve the access mode for an entry in the recipe configuration
 *   database.
 *
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 *
 * @return The function returns the entry's access mode setting.
 *
 * The function searches the recipe configuration database for the
 * given configuration database entry and returns READ_WRITE if the entry
 * may be modified or READ_ONLY if the entry cannot be changed.
 */

PilCdbKeyMode pilDfsDbGetKeyMode(const char *group, const char *name)
{

  return pilCdbGetKeyMode(cfgDB, group, name);

}

/**
 * @brief
 *   Create a new group in the recipe configuration database.
 *
 * @param group   Name of the configuration group.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function creates a recipe configuration group in the recipe
 * configuration database with the groupname @em group. If
 * a group with this name is already present the function fails.
 * The function will also fail if the group to be created is a
 * subgroup entry and its parent group does not yet exist.
 *
 * @internal: Maybe at this interface level it should not be an error
 *            if the parent group does not exist, but instead the full
 *            group hierarchy should be created?
 */

int pilDfsDbCreateGroup(const char *group)
{

  return pilCdbCreateGroup(cfgDB, group);

}


/**
 * @brief
 *   Create a new entry in the recipe configuration database.
 *
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param value   Value string of the parameter.
 * @param access  Access mode of the database entry.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function creates a new entry in the recipe configuration database.
 * The entry is created within the group @em group with the database
 * keyword name @em name. The string @em value is copied to the value
 * field of the new entry and its access mode is set to @em access. If
 * the database group @em group does not exist it is created as long as
 * the parent group already exists or the group itself is a top level group.
 *
 * In case that a database parameter with the name @em name already
 * exists within the group @em group the function will fail. To modify
 * an existing database value @b pilDfsDbModifyValue() has to be used.
 *
 * The function can be used for any kind of the retrievable datatypes.
 * For the numerical types the string representation of the numerical
 * value has to be passed. 
 *
 * @see pilDfsDbModifyValue(), pilDfsDbGetString(), pilDfsDbGetBool(), 
 *   pilDfsDbGetInt(), pilDfsDbGetLong(), pilDfsDbGetFloat(),
 *   pilDfsDbGetDouble()
 */

int pilDfsDbCreateEntry(const char *group, const char *name,
			const char *value, PilCdbKeyMode access)
{

  if (pilCdbCreateEntry(cfgDB, group, name, value) == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (access == READ_ONLY)
    if (pilCdbSetKeyMode(cfgDB, group, name, access) == EXIT_FAILURE)
      return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Retrieve a string value from the recipe configuration database.
 *
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 *
 * @return The value string.
 *
 * The recipe configuration database is queried for the parameter
 * @em name stored in the configuration group @em group. If
 * the parameter is found its value is returned as a string.
 *
 * Actually the function returns the pointer to the value string which
 * is stored in the database. Therefore the value string referenced by
 * the returned pointer must not be modified through this pointer directly.
 * To modify the value string in the database use the interface function
 * @b pilDfsDbModifyValue().
 *
 * @see pilDfsDbModifyValue()
 */

const char *pilDfsDbGetString(const char *group, const char *name)
{

  return pilCdbGetString(cfgDB, group, name);

}


/*
 * @brief
 *   Get a recipe configuration database value as boolean.
 *
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * @return The function returns 1 if the value is true and 0 if the value
 *   is false. If an error occurs the value of the argument errval is
 *   returned.
 *
 * The recipe configuration database is queried for the parameter
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
 * @see pilDfsDbModifyValue()
 */

int pilDfsDbGetBool(const char *group, const char *name, int errval)
{

  return pilCdbGetBool(cfgDB, group, name, errval);

}


/**
 * @brief
 *   Get a recipe configuration database value as integer.
 *
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * @return The numerical value of the database value string as int. If an
 *   error occurs the value of the argument errval is returned.
 *
 * The recipe configuration database is queried for the parameter
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
 * @see pilDfsDbModifyValue()
 */

int pilDfsDbGetInt(const char *group, const char *name, int errval)
{

  return pilCdbGetInt(cfgDB, group, name, errval);

}


/**
 * @brief
 *   Get a recipe configuration database value as long.
 *
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * @return The numerical value of the database value string as a long. If
 *   an error occurs the value of the argument errval is returned.
 *
 * The recipe configuration database is queried for the parameter
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
 * @see pilDfsDbModifyValue()
 */

long pilDfsDbGetLong(const char *group, const char *name, long errval)
{

  return pilCdbGetLong(cfgDB, group, name, errval);

}


/**
 * @brief
 *   Get a recipe configuration database value as float.
 *
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * @return The numerical value of the database value string as a float. If
 *   an error occurs the value of the argument errval is returned.
 *
 * The recipe configuration database is queried for the parameter
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
 * @see pilDfsDbModifyValue()
 */

float pilDfsDbGetFloat(const char *group, const char *name, float errval)
{

  return pilCdbGetFloat(cfgDB, group, name, errval);

}


/**
 * @brief
 *   Get a recipe configuration database value as double.
 *
 * @param group   Name of the configuration group.
 * @param name    Name of the configuration parameter.
 * @param errval  Value to be returned in case of failure.
 *
 * @return The numerical value of the database value string as a double. If
 *   an error occurs the value of the argument errval is returned.
 *
 * The recipe configuration database is queried for the parameter
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
 * @see pilDfsDbModifyValue()
 */

double pilDfsDbGetDouble(const char *group, const char *name,
			     double errval)
{

  return pilCdbGetDouble(cfgDB, group, name, errval);

}


/**
 * @brief
 *   Modify an existing value string of a recipe configuration database
 *   parameter.
 *
 * @param group  Name of the configuration group.
 * @param name   Name of the configuration parameter.
 * @param value  Value string of the parameter.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function updates the value string of an already existing parameter
 * entry, given by @em group and @em name, in the recipe
 * configuration database. The function cannot be used to create new
 * database entries. If the parameter entry does not exist, i.e. either the
 * group or the parameter does not exist, the function returns an error.
 *
 * The function can be used for any kind of the retrievable datatypes.
 * For the numerical types the string representation of the numerical
 * value has to be passed. 
 *
 * @see pilDfsDbGetString(), pilDfsDbGetBool(), pilDfsDbGetInt(),
 *   pilDfsDbGetLong(), pilDfsDbGetFloat(), pilDfsDbGetDouble()
 */

int pilDfsDbModifyValue(const char *group, const char *name,
			const char *value)
{

  return pilCdbModifyValue(cfgDB, group, name, value);

}
/**@}*/
