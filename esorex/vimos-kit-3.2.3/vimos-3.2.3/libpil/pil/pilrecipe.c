/* $Id: pilrecipe.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#ifdef HAVE_GETOPT_LONG
#  include <getopt.h>
#else
#  include "getopt.h"
#endif

#include "pilrecipe.h"
#include "pildfsconfig.h"
#include "pilframeset.h"
#include "pilmemory.h"
#include "pilmessages.h"
#include "piltranslator.h"
#include "pilstrutils.h"
#include "pilfileutils.h"
#include "pilfits.h"
#include "piltimer.h"
#include "pildate.h"
#include "pilutils.h"


#ifndef PRODUCT_DID
#  define PRODUCT_DID  "Unknown"
#endif


/**
 * @defgroup pilRecipe pilRecipe
 *
 * The module @b pilRecipe provides functions to start and initialize the
 * recipe runtime environment, including the internally created recipe
 * configuration database, messaging system, etc. This also includes the
 * parsing of the command line arguments.
 * 
 * Also provided are services to retrieve recipe related information from
 * the system and write the generated products to their output directories.
 * 
 * Finally a function to shutdown the recipe runtime environment is
 * provided.
 */

/**@{*/

/*
 * Parameter group name separator
 */

#define GROUP_SEPARATOR  '.'


/*
 * Default parameter group
 */

static char *group_default = "Parameters";


/*
 * Define the options common to all recipes.
 */

typedef struct _PIL_REC_OPTION_ {
  struct option opt;
  const char *comment;
} PilRecOption;

static PilRecOption common_options[] = {
  {{"copy", no_argument, NULL, 'c'}, 
   "Copy product file instead of moving them."},
  {{"display", optional_argument, NULL, 'd'}, 
   "Turn on displays. The display name may be given."},
  {{"dump", optional_argument, NULL, 'D'},
  "Dump the recipe setup to the given file, or stdout."},
  {{"export", required_argument, NULL, 'e'},
  "Export the product files."},
  {{"export-dir", required_argument, NULL, 'E'},
  "Destination directory for exported product files."},
  {{"graphics", optional_argument, NULL, 'g'}, 
   "Enable graphics output. The device name may be given."},
  {{"help", no_argument, NULL, 'h'}, 
  "Print this message, then exit."},
  {{"log-dir", required_argument, NULL, 'L'},
  "Directory where log files will be stored."},
  {{"log-level", required_argument, NULL, 'l'},
  "Sets the log level to the given value."},
  {{"mode", required_argument, NULL, 'M'},
  "The recipe mode (\"Engineering\", \"Online\" or \"Offline\")."},
  {{"overwrite", no_argument, NULL, 'o'},
  "Overwrite existing product files."},
  {{"product-root", required_argument, NULL, 'R'},
  "Sets the output directory and basename for products."},
  {{"time", no_argument, NULL, 't'},
   "Measure recipe execution time."},
  {{"verbose", optional_argument, NULL, 'v'},
  "Controls the amount of terminal output."},
  {{"version", no_argument, NULL, 'V'},
  "Print version information, then exit."},
  {{0, 0, 0, 0},
   ""}
};


/*
 * String of common recipe options. The leading colon is needed to make
 * getopt_long() return a ':' instead of a '?' if a required argument
 * is not present.
 */

static char common_optstr[] = ":cd::D::e:E:g::hL:l:M:oR:tv::V";


/*
 * Recipe info structure storing the recipe name, version and the
 * start of execution.
 */

typedef struct _PIL_REC_INFO_ {
  char *name;
  char *version;
  char *instrument;
  PilTime start;
  PilTime stop;
  PilTimer *timer;
    
} PilRecInfo;

static PilRecInfo recipe_info = {NULL, NULL, NULL, 0, 0, NULL};


/*
 * @brief
 *   Test if a string is a command line option element.
 *
 * @param string  A string to test.
 *
 * @return The function returns 1 if the string is an option, otherwise
 *   0 is returned.
 *
 * The string @em string is recognized as an option element if it's first
 * character is a '-' either followed by a lower or upper case letter, or
 * another dash which must be followed by a lower or upper case letter.
 */

 static int
is_option(const char *string)
{

  int opt = 0;

  /*
   * Note: The characer class functions have to be enclosed in
   *   parantheses to use the actual function on HP-UX where these
   *   functions are also provided as macros, which are taken by
   *   default and may lead to compiler warnings.
   */

  if (string[0] == '-' && string[1] != '\0') {
    if ((isalpha)(string[1])) {
      opt = 1;
    }
    else {
      if (string[1] == '-' && string[2] != '\0') {
        if ((isalpha)(string[2])) {
        opt = 1;
        }
      }
    }
  }

  return opt;

}


/*
 * @brief
 *   Copy a file.
 *
 * @param srcpath  Source file name.
 * @param dstpath  Destination file name.
 *
 * @return The function returns 0 if no error occurred, otherwise -1 is
 *    returned and errno is set appropriately.
 *
 * The function provides a very basic way to copy the source file 
 * @em srcpath to a destination file @em dstpath.
 *
 * The implementation just uses @b read() and @b write() to do the job.
 * It is by far not comparable with the @b cp shell command. Actually it
 * just writes the source file contents into a new file with the
 * appropriate output name.
 * 
 * If an error occurs the destination file is removed before the function
 * returns.
 */

 static int
copy(const char *srcpath, const char *dstpath)
{

  int src, dst;
  struct stat sb, db;
  char *buf;
  ssize_t rbytes = 0, wbytes = 0, blksize = 4096;


  if ((src = open(srcpath, O_RDONLY)) == -1)
    return -1;
  else {
    if (fstat(src, &sb) == -1 || !S_ISREG(sb.st_mode)) {
      close(src);
      return -1;
    }

    if ((dst = open(dstpath, O_CREAT | O_WRONLY | O_TRUNC,
		    sb.st_mode)) == -1) {
      close(src);
      return -1;
    }
    else 
      if (fstat(dst, &db) == -1 || !S_ISREG(db.st_mode)) {
	close(src);
	close(dst);
	unlink(dstpath);
	return -1;
      }
  }

#ifdef HAVE_ST_BLOCKS
  blksize = db.st_blksize;
#else
#ifdef DEV_BSIZE
  blksize = DEV_BSIZE;
#endif
#endif

  if ((buf = (char *)pil_malloc(blksize)) == NULL) {
    close(src);
    close(dst);
    unlink(dstpath);
    return -1;
  }

  while ((rbytes = read(src, buf, blksize)) > 0)
    if ((wbytes = write(dst, buf, rbytes)) != rbytes) {
      wbytes = -1;
      break;
    }

  close(src);
  close(dst);
  pil_free(buf);
  
  if (rbytes == -1 || wbytes == -1) {
    unlink(dstpath);
    return -1;
  }
  
  return 0;

}


/*
 * @brief
 *   Move a file.
 *
 * @param srcpath  Source file name.
 * @param dstpath  Destination file name.
 *
 * @return The function returns 0 if no error occurred, otherwise -1 is
 *    returned and errno is set appropriately.
 *
 * The function moves the source file @em srcpath to the given
 * destination path @em dstpath.
 * 
 * If an error occurs the destination file is removed before the function
 * returns and the source file is left untouched.
 */

 static int
move(const char *srcpath, const char *dstpath)
{

  struct stat sb;


  if (copy(srcpath, dstpath) == -1)
    return -1;


  /*
   * Remove the source, but if the source file cannot be checked or is not
   * writable revert to the original state, i.e. remove the file copy.
   */ 

  if (stat(srcpath, &sb) == -1 || !(sb.st_mode & S_IWUSR)) {
    unlink(dstpath);
    return -1;
  }

  unlink(srcpath);
  return 0;

}


/*
 * @brief
 *   Destructor for the recipe information structure.
 *
 * @return Nothing.
 *
 * The function just frees the memory blocks used for the recipe name
 * and version string. It also assigns @c NULL to the pointers.
 */

 static void
deleteRecInfo(void)
{

  if (recipe_info.name) {
    pil_free(recipe_info.name);
    recipe_info.name = NULL;
  }

  if (recipe_info.version) {
    pil_free(recipe_info.version);
    recipe_info.version = NULL;
  }

  if (recipe_info.instrument) {
    pil_free(recipe_info.instrument);
    recipe_info.instrument = NULL;
  }

  if (recipe_info.timer) {
      deletePilTimer(recipe_info.timer);
      recipe_info.timer = NULL;
  }

  return;

}


/*
 * @brief
 *   Parse out the group from a parameter name.
 *
 * @param name  Parameter name string to parse.
 *
 * @return The function returns the part of the parameter name specifying the
 *   parameter group as a new string. If the parameter name does not have
 *   a group part the function returns @c NULL.
 *
 * The function scans the string @em name for the last occurrence of the
 * group separator character @c GROUP_SEPARATOR and returns a string
 * containing the full group specification, i.e. all the characters up to,
 * but not including, the last group separator. The string containing the
 * group name can be deallocated using @b pil_free().
 */

 static char *
pilRecGetGroupName(const char *name) {

  char *group = NULL;
  char *s = strrchr(name, GROUP_SEPARATOR);

  if (s) {
    size_t sz = s - name;

    group = (char *)pil_malloc((sz + 1) * sizeof(char));
    if (!group) {
      errno = ENOMEM;
      return NULL;
    }
    else {
      strncpy(group, name, sz);
      group[sz] = '\0';
    }
  }

  return group;

}


/*
 * @brief
 *   Merge common recipe options with recipe parameters.
 *
 * @param ropt  Common recipe options
 * @param rpar  Recipe parameters
 *
 * @return The function returns a pointer to the merged list of common options
 *   and recipe parameters if no error occurred, otherwise the return value
 *   is @c NULL.
 *
 * The function creates the full list of command line options and puts
 * them into a structure of type @c struct option appropriate for
 * @b getopt_long().
 * 
 * The recipe parameters are inserted into the created output structure
 * assuming that all recipe parameters have to be followed by an argument.
 */

 static struct option *
pilRecMergeOptions(PilRecOption *ropt, PilRecParameter *rpar)
{

  register int i, j;

  int nc = sizeof(common_options) / sizeof(PilRecOption) - 1;
  int np = 0, nopt = 0;

  struct option *options;



  /*
   * Compute the total number of allowed command line options.
   */

  if (rpar)
    while (rpar[np].name)
      np++;

  nopt = nc + np;


  /*
   * Create the list of merged options. Add 1 for the terminating entry.
   */

  if ((options = (struct option *)pil_malloc((nopt + 1) * 
                                             sizeof(struct option))) == NULL)
    return NULL;
  else {
    for (i = 0; i < nc; i++) {
      options[i].name = pil_strdup(ropt[i].opt.name);
      options[i].has_arg = ropt[i].opt.has_arg;
      options[i].flag = ropt[i].opt.flag;
      options[i].val = ropt[i].opt.val;
    }

    if (np > 0)
      for (i = 0, j = nc; i < np; i++, j++) {
	options[j].name = pil_strdup(rpar[i].name);
	options[j].has_arg = required_argument;
	options[j].flag = NULL;
	options[j].val = 0;
      }

    /*
     * Write terminating entry
     */

    options[nopt].name = NULL;
    options[nopt].has_arg = no_argument;
    options[nopt].flag = NULL;
    options[nopt].val = 0;
  }
      
  return options;

}


/*
 * @brief
 *   Destroy a merged option structure.
 *
 * @param options  Merged option structure.
 *
 * @return Nothing.
 *
 * The function deallocates the name field fore each option found in the
 * structure @em options and finally deallocates the structure itself.
 */

 static void
pilRecFreeMergedOptions(struct option *options)
{

    register int i = 0;

    while (options[i].name) {
        pil_free((void *)options[i].name);
        i++;
    }

    pil_free(options);

    return;
    
}


/*
 * @brief
 *   Print recipe usage to stderr.
 *
 * @param recipe  The recipe name.
 * @param parg    Structure defining recipe parameters.
 *
 * @return Nothing.
 *
 * The function prints a standard usage message to @b stderr. It does
 * not make use of the messaging subsystem. If @em parg is not @c NULL
 * the list of allowed recipe parameters is also printed, taking the 
 * parameter descriptions from the parameter objects comment field.
 */

 static void
pilRecUsage(PilRecParameter *parg)
{

  register int i;

  int terminalWidth = 80;               /* Just a patch */

  static unsigned short flag = 0;
  static size_t sz = 0;
  size_t indent = 11;                   /* Minimum indentation */
  size_t length;

  char optfmt[64], optname[4096];
  char message[4096];
  const char *recipe = pilRecGetName();


  /*
   * IF the recipe name is not set use a default.
   */

  if (!recipe)
    recipe = "Unknown";


  /*
   * On the first call get the size of the longest option name.
   */

  if (!flag) {
    flag = 1;

    i = 0;
    while (common_options[i].opt.name) {
      length = strlen(common_options[i++].opt.name);
      sz = sz < length ? length : sz;
    }

    if (parg) {
      i = 0;
      while (parg[i].name) {
	length = strlen(parg[i++].name);
	sz = sz < length ? length : sz;
      }
    }

    indent += sz;
  }


  /*
   * Print the message
   */

  fprintf(stderr, "Usage: %s [option...] [recipe parameter...] "
	  "set-of-frames...\n\n", recipe);
  fprintf(stderr, "Recipe control options:\n");

  i = 0;
  while (common_options[i].opt.name) {

    /*
     * Note: The characer class functions have to be enclosed in
     *   parantheses to use the actual function on HP-UX where these
     *   functions are also provided as macros, which are taken by
     *   default and may lead to compiler warnings.
     */

    if ((isalnum)(common_options[i].opt.val) && 
	strchr(common_optstr, common_options[i].opt.val)) {

      sprintf(optname, "%s, -%c", common_options[i].opt.name,
	      common_options[i].opt.val);
      sprintf(optfmt, "  --%%-%lus   %%s\n", (unsigned long)sz + 4);
      sprintf(message, optfmt, optname, common_options[i].comment);
      fprintf(stderr, "%s\n", strsplit(message, indent, terminalWidth));
    }
    else {
      sprintf(optfmt, "  --%%-%lus       %%s\n", (unsigned long)sz);
      sprintf(message, optfmt, common_options[i].opt.name,
              common_options[i].comment);
      fprintf(stderr, "%s\n", strsplit(message, indent, terminalWidth));

    }
    i++;
  }


  /*
   * This is printed only if there is at least one parameter
   */

  if (parg && parg[0].name) {

    fprintf(stderr, "\nRecipe parameters:\n");

    i = 0;
    sprintf(optfmt, "  --%%-%lus       %%s [%%s]\n", (unsigned long)sz);
    while (parg[i].name) {
      sprintf(message, optfmt, parg[i].name, parg[i].comment, parg[i].value);
      fprintf(stderr, "%s\n", strsplit(message, indent, terminalWidth));
      i++;
    }

    fprintf(stderr, "\nParameter built-in defaults are given in square "
	    "brackets.\n");
  }

  return;

}


/*
 * @brief
 *   Print the recipe/pipeline version string to stderr.
 *
 * @return Nothing.
 *
 * Prints the standard recipe/pipeline version information to 
 * @b stderr. The recipe information is taken from the
 * global recipe info structure @c recipe_info.
 */

 static void
pilRecVersion(void)
{

  const char *recipe = pilRecGetName();
  const char *version = pilRecGetVersion();
  const char *instrument = pilRecGetInstrument();


  /*
   * Use defaults if the name, version or instrument are not set
   */

  if (!recipe)
    recipe = "Unknown";

  if (!version)
    version = "Unknown";

  if (!instrument)
    instrument = "Unknown";

  fprintf(stderr, "%s %s\n", recipe, version);
  fprintf(stderr, "The program is part of the %s Instrument pipeline.\n\n",
	  instrument);
  fprintf(stderr, "Copyright (C) 2002-2004 European Southern Observatory\n");
  fprintf(stderr, "This program is free software; you can redistribute it "
	  "and/or\n");
  fprintf(stderr, "modify it under the terms of the GNU General Public "
	  "License as\n");
  fprintf(stderr, "published by the Free Software Foundation; either "
	  "version 2 of\n");
  fprintf(stderr, "the License, or (at your option) any later version.\n\n");
  fprintf(stderr, "This program is distributed in the hope that it will be "
	  "useful,\n");
  fprintf(stderr, "but WITHOUT ANY WARRANTY; without even the implied "
	  "warranty of\n");
  fprintf(stderr, "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
  fprintf(stderr, "See the GNU General Public License for more details.\n");

  return;

}


/*
 * @brief
 *   Write recipe start time to the standard output.
 *
 * @return Nothing.
 *
 * The function gets the recipe name and start time from the recipe
 * information structure and writes a message indicating the recipe's name 
 * and starting time to the standard output.
 */

 static void
pilRecTimerStart()
{

  fprintf(stdout, "Recipe: %s started at: %s\n", recipe_info.name, 
          pilTimerGetCalendarTime(recipe_info.start));

  return;

}


/*
 * @brief
 *   Computes execution time and print it to the standard output.
 *
 * @return Nothing.
 *
 * The function computes recipe execution time using the timer of the
 * recipe. The recipes termination time is computed from the elapsed time
 * and the recipes start time. The stopping time and the number of seconds
 * elapsed since recipe start are written to the standard output.
 */

 static void
pilRecTimerStop(void)
{

  PilTime elapsed;


  if (recipe_info.timer) {
    if (pilTimerIsActive(recipe_info.timer))
      pilTimerStop(recipe_info.timer, NULL);

    elapsed = pilTimerElapsed(recipe_info.timer, NULL);

    fprintf(stdout, "Recipe %s terminated at: %s\n", recipe_info.name,
            pilTimerGetCalendarTime(recipe_info.start + elapsed));
    fprintf(stdout, "Elapsed time: %.3f seconds\n", elapsed);
  }

  return;
}


/*
 * @brief
 *   Set the log level of the recipe.
 *
 * @param log_level  Log level identifier.
 *
 * @return Nothing.
 *
 * The function converts the log level identifier string @em log_level
 * into the appropriate code number for the messaging subsystem and
 * sets the log level through @b pilMsgEnableLog().
 * 
 * If the identifier @em log_level is not recognized as a valid
 * level for the messaging system the current setup of the
 * messaging subsystem is not changed.
 *
 * @see pilMsgEnableLog()
 */

 static void
pilRecSetLogLevel(const char *log_level)
{

  if (!strcmp(log_level, "Off"))
    pilMsgEnableLog(PIL_MSG_OFF);
  else if (!strcmp(log_level, "Debug"))
    pilMsgEnableLog(PIL_MSG_DEBUG);
  else if (!strcmp(log_level, "Info"))
    pilMsgEnableLog(PIL_MSG_INFO);
  else if (!strcmp(log_level, "Warning"))
    pilMsgEnableLog(PIL_MSG_WARNING);
  else if (!strcmp(log_level, "Error"))
    pilMsgEnableLog(PIL_MSG_ERROR);

  return;

}


/*
 * @brief
 *   Set the verbosity of the recipe.
 *
 * @param verbosity  Verbosity identifier.
 *
 * @return Nothing.
 *
 * The function converts the verbosity identifier string @em verbosity
 * into the appropriate code number for the messaging subsystem and
 * sets the verbosity through @b pilMsgEnableTerminal().
 * 
 * If the identifier @em verbosity is not recognized as a valid
 * level for the messaging system the current setup of the messaging
 * subsystem is not changed.
 *
 * @see pilMsgEnableTerminal()
 */

 static void
pilRecSetVerbosity(const char *verbosity)
{

  if (!strcmp(verbosity, "Off"))
    pilMsgEnableTerminal(PIL_MSG_OFF);
  else if (!strcmp(verbosity, "Debug"))
    pilMsgEnableTerminal(PIL_MSG_DEBUG);
  else if (!strcmp(verbosity, "Info"))
    pilMsgEnableTerminal(PIL_MSG_INFO);
  else if (!strcmp(verbosity, "Warning"))
    pilMsgEnableTerminal(PIL_MSG_WARNING);
  else if (!strcmp(verbosity, "Error"))
    pilMsgEnableTerminal(PIL_MSG_ERROR);

  return;

}


/*
 * @brief
 *   Update the recipe configuration database with the builtin recipe
 *   parameter defaults.
 *
 * @param rpar  Recipe task parameters
 *
 * @return The function returns @c EXIT_SUCCESS if the database was updated
 *   successfully, otherwise the return value is @c EXIT_FAILURE.
 *
 * The function adds the name of the recipe task parameters and its
 * default value provided by the set of task parameters @em rpar
 * to the recipe configuration database. The task parameters are
 * added to the configuration group @b Parameters. If the 
 * configuration group does not exist it is created.
 */

 static int
pilRecParameterDefaults(PilRecParameter *rpar)
{

  size_t np = 0;
    
  if (!rpar)
    return EXIT_FAILURE;

  while (rpar[np].name) {
    int _errno;

    char *group;
    char *name = (char *)rpar[np].name;


    _errno = errno;
    errno = 0;

    group = pilRecGetGroupName(rpar[np].name);
    if (!group) {
      if (!errno) {
        errno = _errno;
        group = group_default;
      }
      else {
        return EXIT_FAILURE;
      }
    }

    if (group != group_default)
      name = name + strlen(group) + 1;

    if (!pilDfsDbGroupExists(group)) {
      if (pilDfsDbCreateGroup(group) == EXIT_FAILURE)
        return EXIT_FAILURE;
    }

    if (pilDfsDbEntryExists(group, name)) {
      if (pilDfsDbModifyValue(group, name,
                              rpar[np].value) == EXIT_FAILURE) {
        return EXIT_FAILURE;
      }
    }
    else {
      if (pilDfsDbCreateEntry(group, name, rpar[np].value,
                              READ_WRITE) == EXIT_FAILURE) {
        return EXIT_FAILURE;
      }
    }

    if (group != group_default)
      pil_free(group);

    np++;
  }

  return EXIT_SUCCESS;

}


/*
 * @brief
 *   build up the full product file name.
 *
 * @param name    Destination string for the fully resolved product name.
 * @param path    Path to the product directory.
 * @param prefix  Unresolved prefix of the product filename.
 * @param count   Output product file running number.
 * @param frame   Source frame for which the product name is generated.
 *
 * @return The function returns the pointer to the fully resolved product
 *   filename string. If the product name cannot be resolved the function
 *   returns @c NULL.
 *
 * The function generates a file name for the pipeline product frame
 * @em frame from the path to the output directory @em path, the file
 * name root @em prefix and the output file running number @em count.
 * The filename extension is determined from the frame's datatype.
 * 
 * If there is no path the argument @em path has to be set to @c NULL.
 * Otherwise the provided path name should not have a slash as the last
 * character, since it is inserted by @b pilRecBuildProductName().
 * 
 * If no running number should be appended to the output file name
 * the argument @em count must be smaller than 0.
 * 
 * The file name root @em prefix may be a simple string which is used
 * literally as basename, but it may also be one of the special strings:
 *
 *   @li @b recipe(): The recipe's name is used as the product file
 *                    name root.
 *   @li @b category(): The frame category identifier is retrieved from the
 *                      source frame, converted to lower case and used as
 *                      file name root. This method suppresses the product
 *                      file running number in any case.
 */

 static char *
pilRecBuildProductName(char *name, const char *path, const char *prefix,
                       int count, const PilFrame *frame)
{

  size_t sz, ndigit = 1, failed = 0;

  char  category[PIL_CATEGORY_MAX + 1];
  const char dummy_path[] = "./";


  if (!path)
    path = dummy_path;


  /*
   * If the running number should be used get the number of digits needed.
   */

  if (count > 0)
    ndigit = log10(count) + 1;
  ndigit = ndigit > 4 ? ndigit : 4;


  /*
   * Resolve the file name root. If the prefix string does not contain
   * the token '()' it is used literally.
   */

  if (!strstr(prefix, "()")) {
    sz = strlen(path) + strlen(prefix);

    if (count < 0) {
      if (sz + 1 <= PIL_PATHNAME_MAX)
	sprintf(name, "%s/%s", path, prefix);
      else
	return NULL;
    }
    else {
      if (sz + ndigit + 2 <= PIL_PATHNAME_MAX)
	sprintf(name, "%s/%s_%04d", path, prefix, count);
      else
	return NULL;
    }
  }
  else {

    /*
     * Check if the special string corresponds to an existing conversion
     * function. If it is not a supported conversion abort.
     */

    if (strcmp(prefix, "recipe()") && strcmp(prefix, "category()"))
      return NULL;

    /*
     * Do the conversion.
     */

    if (!strcmp(prefix, "recipe()")) {
      sz = strlen(path) + strlen(pilRecGetName());

      if (count < 0) {
	if (sz + 1 <= PIL_PATHNAME_MAX)
	  sprintf(name, "%s/%s", path, pilRecGetName());
	else
	  return NULL;
      }
      else {
	if (sz + ndigit + 2 <= PIL_PATHNAME_MAX)
	  sprintf(name, "%s/%s_%04d", path, pilRecGetName(), count);
	else
	  return NULL;
      }
    }

    if (!strcmp(prefix, "category()")) {
      sz = strlen(path) + strlen(pilFrmGetCategory(frame));

      if (sz + 1 <= PIL_PATHNAME_MAX) {
	strncpy(category, pilFrmGetCategory(frame), PIL_CATEGORY_MAX);
	strlower(category);
	sprintf(name, "%s/%s", path, category);
      }
      else
	return NULL;
    }
  }


  /*
   * Determine and append the file name extension from the frame's
   * datatype.
   */

  sz = strlen(name);

  switch (pilFrmGetFormat(frame)) {
  case PIL_FRAME_FORMAT_IMAGE:
  case PIL_FRAME_FORMAT_TABLE:
    if (sz + strlen(".fits") <= PIL_PATHNAME_MAX)
      strcat(name, ".fits");
    else 
      failed = 1;
    break;

/******
  case PIL_FRAME_FORMAT_TABLE:
    if (sz + strlen(".tfits") <= PIL_PATHNAME_MAX)
      strcat(name, ".tfits");
    else 
      failed = 1;
    break;
 ******/

  case PIL_FRAME_FORMAT_PAF:
    if (sz + strlen(".paf") <= PIL_PATHNAME_MAX)
      strcat(name, ".paf");
    else 
      failed = 1;
    break;

  default:
    failed = 1;
    break;
  }

  if (failed) {
    *name = '\0';
    return NULL;
  }

  return name;

}


/*
 * @brief
 *   Copy a pipeline product file to the product directory.
 *
 * @param frame  Frame that should be saved.
 * @param path   Absolute product path name.

 * @param path   Path to the product directory.
 * @param name   Basename of the output product file.
 * @param count  Product sequence number.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function writes the pipeline product given by @em frame to its final
 * destination.
 *
 * If the product destination file given by @em path already exists in the
 * pipeline product directory the function's behaviour depends on the value
 * of the recipe database entry @b OverwriteProducts in the group
 * @b DfsConfig. If products may be overwritten the file permissions are
 * checked and modified if necessary, otherwise the function returns an error.
 * 
 * If the database entry @b CopyProducts in the group @b DfsConfig is
 * set to @em true the function will copy the local product file to the
 * target directory. Otherwise the local products are moved.
 */

 static int
pilRecCopyProduct(PilFrame *frame, const char *path)
{

  const char fctid[] = "pilRecCopyProduct";

  char srcfile[PIL_PATHNAME_MAX + 1];


  if (!frame)
    return EXIT_FAILURE;


  /*
   * Check if the input file is readable
   */

  strcpy(srcfile, pilFrmGetName(frame));

  if (access(srcfile, R_OK)) {
    pilMsgError(fctid, "%s: Permission denied", srcfile);
    return EXIT_FAILURE;
  }


  /*
   * Check if the output file exists already
   */

  if (!access(path, F_OK)) {
    if (pilDfsDbGetBool("DfsConfig", "OverwriteProducts", 0)) {
      if (access(path, R_OK | W_OK)) {
	if (chmod(path, 0644) == -1) {
	  pilMsgError(fctid, "Setting permissions failed for %s", path);
	  return EXIT_FAILURE;
	}
      }
    }
    else {
      pilMsgError(fctid, "Product file '%s' already exists!", path);
      return EXIT_FAILURE;
    }
  }

  if (!pilDfsDbGetBool("DfsConfig", "CopyProducts", 0)) {
    if (move(srcfile, path) == -1) {
      pilMsgError(fctid, "Cannot move local file: %s", srcfile);
      return EXIT_FAILURE;
    }
  }
  else {
    if (copy(srcfile, path) == -1) {
      pilMsgError(fctid, "Cannot copy local file: %s", srcfile);
      return EXIT_FAILURE;
    }
  }


  /*
   * Make file read only for everyone.
   */

  if (chmod(path, 0444) == -1)
    pilMsgWarning(fctid, "Cannot change file permissions: %s", path);

  pilMsgInfo(fctid, "Product file %s (%s) written to %s", srcfile,
	     pilFrmGetCategory(frame), path);

  return EXIT_SUCCESS;
	   
}


/*
 * @brief
 *   Export a pipeline product file to the archive.
 *
 * @param e_path  Path to the export directory.
 * @param p_path  Absolute path to the product file to export.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The input file given by @em p_path is located in the product
 * directory and it is verified that the file exists and is readable.
 * The product file is exported by creating a link (see the note below)
 * to the product file in the pipeline export directory given by @em e_path.
 * The name of the link will be the same as the name of the product file.
 * 
 * The export directory name @em e_path should not contain a trailing slash.
 * 
 * @note
 *   In case symbolic links are not supported by the system the product
 *   file is simply copied. To copy the file is safer, since hard
 *   links cannot be used across file systems.
 */

 static int
pilRecExportProduct(const char *e_path, const char *p_path)
{

  const char fctid[] = "pilRecExportProduct";

  char linkname[PIL_PATHNAME_MAX + 1];
  char *s;

  struct stat lb;


  s = strrchr(p_path, '/');
  s++;

  sprintf(linkname, "%s/%s", e_path, s);


  /*
   * Check if the input file is readable
   */

  if (access(p_path, R_OK)) {
    pilMsgError(fctid, "Product file is unreadable: %s", p_path);
    return EXIT_FAILURE;
  }
   

  /*
   * To export the file to the archive usually a symbolic link to the
   * product file is created. But symbolic links are not conforming to
   * the ANSI C standard. If symbolic links are not available the file
   * is just copied. Hard links might be an alternative but they cannot
   * be used across file systems, so a real copy is created here.
   */

#if defined HAVE_LSTAT && defined HAVE_SYMLINK
  if (lstat(linkname, &lb)) {
    if (errno != ENOENT) {
      pilMsgError(fctid, "Cannot get file status: %s", linkname);
      return EXIT_FAILURE;
    }
  }
  else {
    if (S_ISLNK(lb.st_mode))
      unlink(linkname);
    else {
      pilMsgError(fctid, "Not a symbolic link: %s", linkname);
      return EXIT_FAILURE;
    }
  }

  if (symlink(p_path, linkname) != 0) {
    pilMsgError(fctid, "Cannot create symbolic link for %s", p_path);
    return EXIT_FAILURE;
  }
#else
  if (stat(linkname, &lb)) {
    if (errno != ENOENT) {
      pilMsgError(fctid, "Cannot get file status: %s", linkname);
      return EXIT_FAILURE;
    }
  }
  else 
    unlink(linkname);

  if (copy(p_path, linkname) != 0) {
    pilMsgError(fctid, "Cannot copy %s", p_path);
    return EXIT_FAILURE;
  }
#endif

  return EXIT_SUCCESS;
	   
}


/*
 * @brief
 *   Save the recipe logfile.
 *
 * @param filename  Name of the local logfile.
 * @param path      Path to the logfile's target location.
 * @param prefix    Basename of the output logfile.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The output file name for the logfile is generated from the target
 * location given by @em path and the basename of the output logfile
 * @em prefix. The output logfile name is completed by appending the
 * suffix ".log" internally. The target directory name @em path should
 * not contain a trailing slash.
 * 
 * If this file already exists in the target directory the function's
 * behaviour depends on the value of the recipe database entry 
 * @b OverwriteProducts in the group @b DfsConfig. If products may be
 * overwritten the file permissions are checked and modified if necessary,
 * otherwise the function returns an error.
 * 
 * If the database entry @b CopyProducts in the group @b DfsConfig is set
 * to @em true the function will copy the locally created logfile to its
 * target directory. Otherwise the local logfile is moved to that directory.
 */

 static int
pilRecSaveLog(const char *filename, const char *path, const char *prefix)
{

  const char fctid[] = "pilRecSaveLog";

  char cwd[PIL_PATHNAME_MAX + 1];
  char srcfile[PIL_PATHNAME_MAX + 1], dstfile[PIL_PATHNAME_MAX + 1];


  /*
   * Build the fully qualified source and target logfile names.
   */

  getcwd(cwd, PIL_PATHNAME_MAX);
  sprintf(srcfile, "%s/%s", cwd, filename);

  sprintf(dstfile, "%s/%s.log", path, prefix);


  /*
   * Check if the local logfile is readable.
   */

  if (access(filename, R_OK)) {
    pilMsgError(fctid, "Local logfile is unreadable: %s", filename);
    return EXIT_FAILURE;
  }


  /*
   * Check if the target logfile exists already
   */

  if (!access(dstfile, F_OK)) {
    if (pilDfsDbGetBool("DfsConfig", "OverwriteProducts", 0)) {
      if (access(dstfile, R_OK | W_OK)) {
	if (chmod(dstfile, 0644) == -1) {
	  pilMsgError(fctid, "Setting permissions failed for %s", dstfile);
	  return EXIT_FAILURE;
	}
      }
    }
    else {
      pilMsgError(fctid, "Product file '%s' already exists!", dstfile);
      return EXIT_FAILURE;
    }
  }

  if (!pilDfsDbGetBool("DfsConfig", "CopyProducts", 0)) {
    if (move(srcfile, dstfile) == -1) {
      pilMsgError(fctid, "Cannot move local logfile: %s", filename);
      return EXIT_FAILURE;
    }
  }
  else {
    if (copy(srcfile, dstfile) == -1) {
      pilMsgError(fctid, "Cannot copy local logfile: %s", filename);
      return EXIT_FAILURE;
    }
  }


  /*
   * Make logfile read only for everyone.
   */

  if (chmod(dstfile, 0444) == -1)
    pilMsgWarning(fctid, "Cannot change file permissions: %s", dstfile);

  pilMsgDebug(fctid, "Recipe logfile %s saved as %s", filename, dstfile);

  return EXIT_SUCCESS;
	   
}


/**
 * @brief
 *   Initialize pipeline runtime environment at recipe start.
 *
 * @param instrument  Instrument name for which the pipeline environment
 *                    should be configured.
 * @param recipe      Recipe identifier string
 * @param version     Recipe version string
 * @param parg        Recipe task parameters
 * @param argv        Command line arguments
 * @param argc        Number of command line parameters
 *
 * @return The function returns a pointer to a newly created set of frames if
 *   no error occurred, otherwise the function returns @c NULL.
 *
 * The function initializes the recipe log file and retrieves then 
 * the pipeline configuration settings from the runtime environment,
 * if such information is present.
 * It then parses the site configuration files and, if present, the user
 * specific configuration files. Finally the command line is parsed for any 
 * pipeline configuration option. From this information the pipeline
 * runtime environment is setup using the following scheme: if the user
 * is allowed to override the site configuration the user specific
 * configuration takes precedence over the site configuration. In the same
 * way, environment variables will override settings from any configuration
 * file. Command line options have the highest precedence.
 * 
 * If the setup of the pipeline runtime environment was successful, the
 * function creates a new set of frames from the input files given on the
 * command line. The command line is also parsed for recipe parameters,
 * which are set within the pipeline runtime environment.
 * 
 * All command line options processed by the function are removed from
 * the command line. Unprocessed commandline options are left in the
 * argument list and the argument count is adjusted accordingly, i.e.
 * additional options may be passed to the recipe and can be processed
 * directly by the application.
 */

PilSetOfFrames *
pilRecStart(const char *instrument, const char *recipe, const char *version,
            PilRecParameter *parg, int argc, char **argv)
{

  const char fctid[] = "pilRecStart";

  register int i;

  int c;
  int skip = 0;
  int option_index = 0;
  int last_optind = optind;
  int dump = 0, timer = 0, status = 0;

  char *s, *path, *dump_file = NULL;

  PilSetOfFrames *set;

  struct option *options = pilRecMergeOptions(common_options, parg);



  if (!options)
    return NULL;

  /*
   * Disable error message from getopt itself
   */

  opterr = 0;


  /*
   * Activate swapping mechanism
   */

  /*
  pil_enable_swap();
  */

  /*
   * Put the recipe name and version into the recipe information structure
   */

  atexit(deleteRecInfo);

  if (pilRecSetName(recipe) == EXIT_FAILURE)
    return NULL;

  if (pilRecSetVersion(version) == EXIT_FAILURE)
    return NULL;

  if (pilRecSetInstrument(instrument) == EXIT_FAILURE)
    return NULL;


  /*
   * Create the configuration database.
   */

  if (pilDfsCreateDB('.', USE_CASE) == EXIT_FAILURE) {
    fprintf(stderr, "%s: %s: Configuration database creation failed!\n",
	    recipe, fctid);
    return NULL;
  }
  else
    atexit(pilDfsFreeDB);


  /*
   * Initialize the logging subsystem with the built in defaults.
   */

  if (pilMsgStart() == EXIT_FAILURE) {
    fprintf(stderr, "%s: %s: Cannot startup messaging system!\n",
	    recipe, fctid);
    return NULL;
  }
  else {
    pilMsgSetRecipeName(recipe);
    pilRecSetLogLevel((char *)pilDfsDbGetString("DfsConfig", "LogLevel"));
    pilRecSetVerbosity((char *)pilDfsDbGetString("DfsConfig", "Verbosity"));
  }


  /*
   * Update the database with the built in defaults of the recipe
   * parameters.
   */

  assert(pilRecParameterDefaults(parg) != EXIT_FAILURE);


  /*
   * Update the database with the contents of the configuration files.
   */

  if (pilDfsReadSetupFiles(instrument, recipe) == EXIT_FAILURE)
    fprintf(stderr, "%s: %s: Cannot read configuration files, starting "
	    "from built-in defaults!\n", recipe, fctid);

  /*
   * Update the database with environment variables
   */

  assert(pilDfsGetEnv() != EXIT_FAILURE);


  /*
   * Parse the command line
   */

  while ((c = getopt_long(argc, argv, common_optstr, options, 
			  &option_index)) != -1) {

    switch (c) {
    case 'h':
      pilRecUsage(parg);
      return NULL;
      break;

    case 'V':
      pilRecVersion();
      return NULL;
      break;

    case 'D':
      dump = 1;
      if (optarg) {
	if (!(s =  pilFileExpandFilePath(optarg))) {
	  fprintf(stderr, "%s: %s: Invalid file name!\n", recipe, optarg);
	  return NULL;
	}
	else 
	  dump_file = pil_strdup(s);
      }
      break;

    case 'e':
      if (!strncmp(optarg, "NoExport", 8) || !strncmp(optarg, "MainOnly", 8) ||
	  !strncmp(optarg, "AllProducts", 11))
	pilDfsDbModifyValue("DfsConfig", "ExportProducts", optarg);
      else {
	fprintf(stderr, "%s: Invalid product export mode!\n", recipe);
	return NULL;
      }
      break;

    case 'E':
      if (is_option(optarg)) {
	fprintf(stderr, "%s: Option `%s' requires an argument!\n",
		recipe, options[option_index].name);
	return NULL;
      }
        
      pilDfsDbModifyValue("DfsConfig", "ExportDir", pilFileTrimPath(optarg));
      break;

    case 'v':
      if (optarg) {
	if (!strncmp(optarg, "Off", 3) || !strncmp(optarg, "Debug", 5) ||
	    !strncmp(optarg, "Info", 4) || !strncmp(optarg, "Warning", 7) ||
	    !strncmp(optarg, "Error", 5))
	  pilDfsDbModifyValue("DfsConfig", "Verbosity", optarg);
	else {
	  fprintf(stderr, "%s: Invalid verbosity level!\n", recipe);
	  return NULL;
	}
      }
      else
	pilDfsDbModifyValue("DfsConfig", "Verbosity", "Debug");
      break;

    case 'L':
      if (is_option(optarg)) {
	fprintf(stderr, "%s: Option `%s' requires an argument!\n",
		recipe, options[option_index].name);
	return NULL;
      }
        
      pilDfsDbModifyValue("DfsConfig", "LogDir", pilFileTrimPath(optarg));
      break;

    case 'l':
	if (!strncmp(optarg, "Off", 3) || !strncmp(optarg, "Debug", 5) ||
	    !strncmp(optarg, "Info", 4) || !strncmp(optarg, "Warning", 7) ||
	    !strncmp(optarg, "Error", 5))
	  pilDfsDbModifyValue("DfsConfig", "LogLevel", optarg);
	else {
	  fprintf(stderr, "%s: Invalid logging level!\n", recipe);
	  return NULL;
	}
	break;

    case 'c':
      pilDfsDbModifyValue("DfsConfig", "CopyProducts", "true");
      break;

    case 'o':
      pilDfsDbModifyValue("DfsConfig", "OverwriteProducts", "true");
      break;

    case 'R':
      if (is_option(optarg)) {
	fprintf(stderr, "%s: Option `%s' requires an argument!\n",
		recipe, options[option_index].name);
	return NULL;
      }
        
      s = strrchr(optarg, '/');

      /*
       * Note: The characer class functions have to be enclosed in
       *   parantheses to use the actual function on HP-UX where these
       *   functions are also provided as macros, which are taken by
       *   default and may lead to compiler warnings.
       */

      if (s && (isalnum)(*(s + 1))) {
	pilDfsDbModifyValue("DfsConfig", "ProductPrefix", s + 1);
	*s = '\0';
      }

      pilDfsDbModifyValue("DfsConfig", "ProductDir", pilFileTrimPath(optarg));
      break;

    case 'M':
      if (!strcmp(optarg, "Engineering") || !strcmp(optarg, "Online") ||
	  !strcmp(optarg, "Offline"))
	pilDfsDbModifyValue("DfsConfig", "PipelineMode", optarg);
      else {
	fprintf(stderr, "%s: Invalid pipeline mode!\n", recipe);
	return NULL;
      }
      break;

    case 'd':
      if (optarg)
	if (!strncmp(optarg, "disabled", 8))
	  pilDfsDbModifyValue("Visualization", "EnableDisplays", "false");
	else {
	  pilDfsDbModifyValue("Visualization", "EnableDisplays", "true");
	  if (pilDfsDbEntryExists("Visualization", "Display"))
	    pilDfsDbModifyValue("Visualization", "Display", optarg);
	  else
	    pilDfsDbCreateEntry("Visualization", "Display", optarg,
				READ_WRITE);
	}
      else 
	pilDfsDbModifyValue("Visualization", "EnableDisplays", "true");
      break;

    case 'g':
      if (optarg)
	if (!strncmp(optarg, "disabled", 8))
	  pilDfsDbModifyValue("Visualization", "EnableGraphics", "false");
	else {
	  pilDfsDbModifyValue("Visualization", "EnableGraphics", "true");
	  if (pilDfsDbEntryExists("Visualization", "GraphicsDevice"))
	    pilDfsDbModifyValue("Visualization", "GraphicsDevice", optarg);
	  else
	    pilDfsDbCreateEntry("Visualization", "GraphicsDevice", optarg,
				READ_WRITE);
	}
      else 
	pilDfsDbModifyValue("Visualization", "EnableGraphics", "true");
      break;

    case 't':
      timer = 1;
      break;

    case 0:
      if (is_option(optarg)) {
	fprintf(stderr, "%s: Option `%s' requires an argument!\n",
		recipe, options[option_index].name);
	return NULL;
      }
      else {
        int _errno;
        char *group, *name;

        
        /*
         * At this point the default parameter settings are already in
         * the database, i.e. no need for existance tests or creation
         * of new entries. The existing defaults are just modified.
         */

        _errno = errno;
        errno = 0;

        name = (char *)options[option_index].name;
        group = pilRecGetGroupName(name);

        if (!group) {
          if (!errno) {
            errno = _errno;
            group = group_default;
          }
          else {
            return NULL;
          }
        }

        if (group != group_default)
          name = name + strlen(group) + 1;


        if (pilDfsDbModifyValue(group, name, optarg) == EXIT_FAILURE) {
          fprintf(stderr, "%s: Cannot change recipe parameter entry %s!\n",
                  recipe, options[option_index].name);

          if (group != group_default)
            pil_free(group);

          return NULL;
        }

        if (group != group_default)
          pil_free(group);

      }
      break;

    case ':':
      fprintf(stderr, "%s: Option `%s' requires an argument!\n", recipe,
	      argv[optind - 1]);
      return NULL;
      break;

    default:

      /*
       * Take into account that optind is does not change if getopt
       * processes a list of short options (multiple short options
       * dash starting with a single dash).
       */

      if (last_optind == optind)
	fprintf(stderr, "%s: Invalid option -- %c!\n", recipe, optopt);
      else {
	s = argv[optind - 1];

	/*
	 * Skip leading dash(es)
	 */

	while (*s == '-' && skip < 2) {
	  skip++;
	  s++;
	}

	fprintf(stderr, "%s: Invalid option -- %s!\n", recipe, s);
      }

      fprintf(stderr, "Try `%s --help' for more information.\n", recipe);

      return NULL;
      break;
    }

    if (last_optind != optind)
      last_optind = optind;

  }


  /*
   * Destroy merged options structure
   */

  pilRecFreeMergedOptions(options);


  /*
   * If option dump was used, the contents of the final configuration
   * database is dumped to either stdout or the given file. After the
   * dump the execution stops.
   */

  if (dump) {
    if (pilDfsDumpDB(dump_file) == EXIT_FAILURE) {
      fprintf(stderr, "%s: Error saving recipe configuration!\n", recipe);
      status = 1;
    }
    else
      status = 0;

    if (dump_file)
      pil_free(dump_file);

    return NULL;
  }


  /*
   * Initialize the system
   */

  recipe_info.timer = newPilTimer();
  recipe_info.start = pilTimerStart(recipe_info.timer, NULL);

  if (timer) {
    pilRecTimerStart();
    atexit(pilRecTimerStop);
  }


  /*
   * Setup the logging subsystem for the requested log level and verbosity.
   */

  pilRecSetLogLevel((char *)pilDfsDbGetString("DfsConfig", "LogLevel"));
  pilRecSetVerbosity((char *)pilDfsDbGetString("DfsConfig", "Verbosity"));


  /*
   * Show additional information in debug or engineering mode.
   */

  if (!strcmp(pilDfsDbGetString("DfsConfig", "LogLevel"), "Debug") ||
      !strcmp(pilDfsDbGetString("DfsConfig", "Verbosity"), "Debug")) {
    pilMsgEnableTimeTag();
    pilMsgEnableComponentTag();
  }
  else 
    if (!strcmp(pilDfsDbGetString("DfsConfig", "PipelineMode"),
		"Engineering")) {
      pilMsgEnableComponentTag();
    }


  /*
   * Initialize the build-in header keyword and frame category translation
   * tables. Load external tables if requiested.
   */

  pilTrnInitKeywordMap();
  if (pilDfsDbEntryExists("Dictionaries", "FitsKeywordMap")) {
    path = pilFileExpandDirPath(pilDfsDbGetString("Dictionaries", 
                                "FitsKeywordMap"));

    if (path)
      pilTrnLoadKeywordMap(path);
    else {
      fprintf(stderr, "%s: Could not expand path: %s\n", recipe,
	      pilDfsDbGetString("Dictionaries", "FitsKeywordMap"));
      return NULL;
    }
  }

  pilTrnInitCategoryMap();
  if (pilDfsDbEntryExists("Dictionaries", "FrameCategoryMap")) {
    path = pilFileExpandDirPath(pilDfsDbGetString("Dictionaries", 
				"FrameCategoryMap"));

    if (path)
      pilTrnLoadCategoryMap(path);
    else {
      fprintf(stderr, "%s: Could not expand path: %s\n", recipe,
	      pilDfsDbGetString("Dictionaries", "FrameCategoryMap"));
      return NULL;
    }
  }


  /*
   * There must be at least one set of frames on the command line.
   */

  if (argc <= optind) {
    fprintf(stderr, "%s: Missing set of frames!\n", recipe);
    pilRecUsage(parg);
    return NULL;
  }


  /*
   * Create the set of frames.
   */

  if (!(set = newPilSetOfFrames()))
    return NULL;

  for (i = optind; i < argc; i++)
    if (!pilSofRead(argv[i], set)) {
      fprintf(stderr, "%s: %s: Invalid set of frames!\n", recipe, argv[i]);
      deletePilSetOfFrames(set);
      return NULL;
    }


  /*
   * Check that all frames in the set can be accessed for reading.
   */

  if (!pilRecValidateSet(set)) {
    deletePilSetOfFrames(set);
    return NULL;
  }

  return set;

}


void
pilRecStop(PilSetOfFrames *set)
{

  const char fctid[] = "pilRecStop";


  const char *logfile_prefix = pilDfsDbGetString("DfsConfig", "ProductPrefix");

  PilFrame *frame;

  char *tdir;
  char logfile_dir[PIL_PATHNAME_MAX + 1];



  if (!set || !(frame = pilSofFirst(set)))
    return;


  /*
   * Get the absolute path to the logfile destination path
   * from the recipe configuration database entry 'LogDir'.
   */

  tdir = pilFileExpandDirPath(pilDfsDbGetString("DfsConfig", "LogDir"));

  if (tdir)
    strcpy(logfile_dir, tdir);
  else
    return;


  /*
   * If the product prefix is one of the special functions
   * the recipe name is used as log file name prefix.
   */

  if (strstr(logfile_prefix, "()"))
    logfile_prefix = pilRecGetName();


  /*
   * Remove intermediate and temporary product files from the current
   * working directory.
   */

  pilMsgDebug(fctid, "Removing temporary product files from disk ...");

  frame = pilSofFirst(set);
  while (frame) {
    if (pilFrmGetProductLevel(frame) == PIL_PRODUCT_LEVEL_TEMPORARY ||
	pilFrmGetProductLevel(frame) == PIL_PRODUCT_LEVEL_INTERMEDIATE) {
      if (!pilFrmGetKeepFlag(frame)) {
	if (unlink(pilFrmGetName(frame)) == -1) {
	  pilMsgWarning(fctid, "Cannot remove local file %s", 
			pilFrmGetName(frame));
	}
      }
    }

    frame = pilSofNext(set, frame);

  }


  /*
   * Shutdown the logging system and save the logfile.
   */

  
  if (pilMsgLogLevel() != PIL_MSG_OFF) {
    pilMsgCloseLog();
    pilRecSaveLog(pilMsgGetLogFile(), logfile_dir, logfile_prefix);
  }

  pilMsgStop();


  /*
   * Destroy the set of frames.
   */

  deletePilSetOfFrames(set);

  return;

}


/**
 * @brief
 *   List the contents of the set of frames.
 *
 * @param set  Pointer to an existing set of frames.
 *
 * @return Nothing.
 *
 * The function uses the messaging subsystem to print the contents of
 * @em set to the terminal and/or the log file (depending to the
 * current setup of the messaging system). The conents of the set is
 * listed in a standard format.
 */

void
pilRecListSet(PilSetOfFrames *set)
{

  PilFrame *frame = pilSofFirst(set);

  while (frame) {
    pilMsgInfo(pilRecGetName(), "%s (%s)", pilFrmGetName(frame),
	       pilFrmGetCategory(frame));
    frame = pilSofNext(set, frame);
  }

  return;

}

  
/**
 * @brief
 *   Write pipeline products to disk.
 *
 * @param set  Set of frames being searched for pipeline products
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function scans the given set of frames for pipeline products and
 * writes them to the pipeline product directory given by the pipeline
 * environment setup. If product export to the archive is enabled also
 * the symbolic links to the product file(s) are created in the archive
 * directory.
 */

int
pilRecWriteProducts(PilSetOfFrames *set)
{

  const char fctid[] = "pilRecWriteProducts";


  const char *product_prefix = pilDfsDbGetString("DfsConfig", "ProductPrefix");
  const char *export_mode = pilDfsDbGetString("DfsConfig", "ExportProducts");

  PilFrame *frame, *main_product = NULL;

  size_t frame_count = 0, disable_export = 0;

  struct stat pstat, estat;

  char *tdir;
  char product_dir[PIL_PATHNAME_MAX + 1], export_dir[PIL_PATHNAME_MAX + 1];



  if (!set || !(frame = pilSofFirst(set)))
    return EXIT_FAILURE;


  /*
   * Get the absolute path to the product directory
   */

  tdir = pilFileExpandDirPath(pilDfsDbGetString("DfsConfig", "ProductDir"));

  if (tdir)
    strcpy(product_dir, tdir);
  else 
    return EXIT_FAILURE;


  /*
   * Check if the product directory exists and has the appropriate file
   * permissions.
   */

  if (access(product_dir, F_OK)) {
    pilMsgError(fctid, "Pipeline product directory %s does not exist!",
		product_dir);
    return EXIT_FAILURE;
  }

  if (access(product_dir, R_OK | W_OK)) {
    pilMsgError(fctid, "Pipeline product directory %s is not writable!",
		product_dir);
    return EXIT_FAILURE;
  }


  /*
   * If product export is enabled get the absolute path to the export
   * directory check that it has the correct file permissions.
   */

  if (!strcmp(export_mode, "MainOnly") || !strcmp(export_mode,
						  "AllProducts")) {
    tdir = pilFileExpandDirPath(pilDfsDbGetString("DfsConfig", "ExportDir"));

    if (tdir)
      strcpy(export_dir, tdir);
    else 
      return EXIT_FAILURE;

    if (access(export_dir, F_OK)) {
      pilMsgError(fctid, "Product export directory %s does not exist!",
		  export_dir);
      return EXIT_FAILURE;
    }

    if (access(export_dir, R_OK | W_OK)) {
      pilMsgError(fctid, "Pipeline export directory %s is not writable!",
		  export_dir);
      return EXIT_FAILURE;
    }

    /*
     * Check if the product directory and the export directory
     * are identical. If they are exporting the product is
     * disabled to avoid creating a link on itself. Also since
     * the product directory and the export directory are the
     * same creating a link is not necessary at all.
     */

    stat(product_dir, &pstat);
    stat(export_dir, &estat);

    if ((pstat.st_dev == estat.st_dev) && (pstat.st_ino == estat.st_ino)) {
      pilMsgDebug(fctid, "Product and export directory are identical. "
		  "Export disabled internally!");
      disable_export = 1;
    }

  }


  /*
   * Lookup the main product, since it must be written first. Make sure
   * that the main product is defined uniquely. Write the main product
   * to disk and create the archive links if required.
   */

  pilMsgInfo(fctid, "Writing product files ...");

  while (frame) {
    char product_name[PIL_PATHNAME_MAX + 1];

    if (pilFrmGetType(frame) == PIL_FRAME_TYPE_PRODUCT && 
	pilFrmGetProductLevel(frame) == PIL_PRODUCT_LEVEL_PRIMARY) {
      if (!main_product) {
	main_product = frame;

        if (!pilRecBuildProductName(product_name, product_dir,
                                    product_prefix, frame_count, frame)) {
          pilMsgError(fctid, "Cannot build product name for %s",
                      pilFrmGetName(frame));
          return EXIT_FAILURE;
        }

        if (pilFrmGetFormat(frame) != PIL_FRAME_FORMAT_PAF) {
          if (pilRecUpdateProductInfo(frame, product_name, set))
              return EXIT_FAILURE;
        }

	if (pilRecCopyProduct(frame, product_name) ==  EXIT_FAILURE) {
	  return EXIT_FAILURE;
	}

	/*
	 * Create the archive link if product export is enabled and
	 * possible.
	 */

	if (!disable_export && (!strcmp(export_mode, "MainOnly") || 
	    !strcmp(export_mode, "AllProducts"))) {
	  if (pilRecExportProduct(export_dir, product_name) == EXIT_FAILURE) {
	    return EXIT_FAILURE;
	  }
	}

	frame_count++;

      }
      else {
	pilMsgError(fctid, "Multiply defined primary recipe product!");
	return EXIT_FAILURE;
      }
    }

    frame = pilSofNext(set, frame);

  }


  /*
   * The main product must have been defined.
   */

  if (!main_product) {
    pilMsgError(fctid, "Missing primary pipeline product!");
    return EXIT_FAILURE;
  }


  /*
   * Write the remaining products to disk and export them to the archive
   * if necessary.
   */

  frame = pilSofFirst(set);
  while (frame) {
    char product_name[PIL_PATHNAME_MAX + 1];

    if (pilFrmGetType(frame) == PIL_FRAME_TYPE_PRODUCT && 
	pilFrmGetProductLevel(frame) == PIL_PRODUCT_LEVEL_SECONDARY) {

      if (!pilRecBuildProductName(product_name, product_dir,
                                  product_prefix, frame_count, frame)) {
        pilMsgError(fctid, "Cannot build product name for %s",
                    pilFrmGetName(frame));
        return EXIT_FAILURE;
      }

      if (pilFrmGetFormat(frame) != PIL_FRAME_FORMAT_PAF) {
        if (pilRecUpdateProductInfo(frame, product_name, set))
            return EXIT_FAILURE;
      }

      if (pilRecCopyProduct(frame, product_name) == EXIT_FAILURE) {
	return EXIT_FAILURE;
      }

      /*
       * Create the archive link if product export is enabled for auxiliary
       * products too.
       */

      if (!disable_export && !strcmp(export_mode, "AllProducts")) {
	if (pilRecExportProduct(export_dir, product_name) == EXIT_FAILURE) {
	  return EXIT_FAILURE;
	}
      }

      frame_count++;

    }

    frame = pilSofNext(set, frame);

  }

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Reset the recipe information.
 *
 * @return Nothing.
 *
 * The recipe information data currently stored in the recipe information
 * structure is cleared. After this call the contents of the recipe
 * information structure is undefined. Before it is used again it must be
 * re-initialized using pilRecSetName(), pilRecSetVersion() and
 * pilRecSetInstrument()
 */

void
pilRecInfoClear(void)
{

    deleteRecInfo();
    return;

}


/**
 * @brief
 *   Set the recipe name.
 *
 * @param recipe  Recipe name.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function stores the recipe name string in the recipe information
 * structure.
 */

int
pilRecSetName(const char *recipe)
{

  if (recipe_info.name)
    pil_free(recipe_info.name);

  if (!(recipe_info.name = pil_strdup(recipe)))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Retrieve the recipe name.
 *
 * @return Pointer to the recipe identifier string.
 *
 * The function gives back the location of the recipe identification
 * information. The string must not be modified through the returned
 * pointer.
 */

const char *
pilRecGetName(void)
{

  return recipe_info.name;

}


/**
 * @brief
 *   Set the recipe version string.
 *
 * @param recipe  Recipe version.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function stores the recipe version string in the recipe information
 * structure.
 */

int
pilRecSetVersion(const char *version)
{

  if (recipe_info.version)
    pil_free(recipe_info.version);

  if (!(recipe_info.version = pil_strdup(version)))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Retrieve the recipe version information.
 *
 * @return Pointer to the recipe version string.
 *
 * The function gives back the location of the recipe version information
 * string. The string must not be modified through the returned pointer.
 */

const char *
pilRecGetVersion(void)
{

  return recipe_info.version;

}


/**
 * @brief
 *   Set the instrument name.
 *
 * @param instrument  Instrument name.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function stores the instrument name string in the recipe information
 * structure.
 */

int
pilRecSetInstrument(const char *instrument)
{

  if (recipe_info.instrument)
    pil_free(recipe_info.instrument);

  if (!(recipe_info.instrument = pil_strdup(instrument)))
    return EXIT_FAILURE;


  /*
   * Instrument names should always be uppercase.
   */

  strupper(recipe_info.instrument);

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Retrieve the instrument name.
 *
 * @return Pointer to the instrument name string.
 *
 * The function gives back the location of the instrument name string.
 * The string must not be modified through the returned pointer.
 */

const char *
pilRecGetInstrument(void)
{

  return recipe_info.instrument;

}


/**
 * @brief
 *   Set the start time of the recipe.
 * 
 * @param time  The time to set as recipe start time.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function stores the time @em time, as the time the recipe execution
 * started, to the recipe information structure.
 */

int
pilRecSetTimeStart(PilTime time)
{

    recipe_info.start = time;
    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the start time of the recipe.
 * 
 * @return The function returns the recipe's start time.
 *
 * The function retrieves the time when the recipe execution started from the
 * recipe information structure.
 */

PilTime
pilRecGetTimeStart(void)
{

    return recipe_info.start;

}


/**
 * @brief
 *   Set the stop time of the recipe.
 * 
 * @param time  The time to set as recipe stop time.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function stores the time @em time, as the time the recipe execution
 * stopped, to the recipe information structure.
 */

int
pilRecSetTimeStop(PilTime time)
{

    recipe_info.stop = time;
    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the stop time of the recipe.
 * 
 * @return The function returns the recipe's stop time.
 *
 * The function retrieves the time when the recipe execution stopped from the
 * recipe information structure.
 */

PilTime
pilRecGetTimeStop(void)
{

    return recipe_info.stop;

}


/**
 * @brief
 *   Set the timer of a recipe.
 *
 * @param timer  The timer to set as the recipe timer.
 *
 * @return The function returns @c EXIT_SUCCESS on success or @c EXIT_FAILURE
 *   in case an error occurred.
 *
 * The function installs the given timer @em timer as the recipe's execution
 * timer. If there is already a timer installed it is stopped first and then
 * replaced by @em timer.
 */

int
pilRecSetTimer(PilTimer *timer)
{

    if (recipe_info.timer != NULL) {
        if (pilTimerIsActive(recipe_info.timer)) {
            pilTimerStop(recipe_info.timer, NULL);
        }

        deletePilTimer(recipe_info.timer);
        recipe_info.timer = NULL;
    }

    recipe_info.timer = timer;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Get the timer of a recipe.
 *
 * @return The function returns the currently installed recipe timer, or
 *   @c NULL if no timer is installed.
 *
 * The function retrieves the currently installed recipe timer from the
 * recipe's information structure.
 */

PilTimer *
pilRecGetTimer(void)
{

    return recipe_info.timer;

}


/**
 * @brief
 *   Validate the set of input frames.
 *
 * @param set  set of frames to be validated.
 *
 * @return The function returns 1 if the input set of frames is
 *   valid, otherwise it returns 0.
 *
 * The validity of the given set of frames is checked. To be a valid set
 * the frames listed in the set must exist and they have to be readable.
 * Furthermore the must be valid FITS files.
 */

int
pilRecValidateSet(PilSetOfFrames *set)
{

  const char fctid[] = "pilRecValidateSet";

  PilFrame *frame;
  char *name;


  if (!set)
    return 0;

  frame = pilSofFirst(set);

  while (frame) {
    name = (char *)pilFrmGetName(frame);
    if (access(name, F_OK)) {
      pilMsgError(fctid, "%s: No such file or directory!", name);
      return 0;
    }

    if (access(name, R_OK)) {
      pilMsgError(fctid, "%s: Permission denied!", name);
      return 0;
    }

    if (!pilFileIsFits(name)) {
      pilMsgError(fctid, "File '%s' is not a standard FITS file!", name);
      return 0;
    }

    frame = pilSofNext(set, frame);
  }

  return 1;

}


/**
 * @brief
 *   Update the product header with generic product information.
 *
 * @param frame  The frame to update.
 * @param name   The fully expanded name of the product file.
 *
 * @return The function returns 0 if the product information was updated
 *   successfully, or 1 otherwise.
 *
 * The function writes the generic product information to the products. The
 * information written to the product frame @em frame are the:
 *   @li MD5 signature for the products data sections.
 *   @li file name of @em frame as it was created by the recipe.
 *   @li the recipe name and start time.
 *   @li the pipeline version.
 *
 * If @em name is not a NULL pointer, the name given is expected to be the
 * absolute path and file name of the final product file, and it is used
 * as the product name to be written to the header instead of the frame's
 * local name. 
 */

int
pilRecUpdateProductInfo(PilFrame *frame, const char *name, PilSetOfFrames *set)
{

  char *key, *value;
  char *md5, *type;
  char *date;

  register size_t sz = PIL_FITS_VALUE_MAX - 1;

  PilFitsFile *fitsfile;


  assert(frame != NULL);


  md5 = pilFitsMD5Signature(pilFrmGetName(frame));

  if (!md5)
    return EXIT_FAILURE;


  key = (char *)pil_calloc(PIL_FITS_CARD_MAX, sizeof(char));
  value  = (char *)pil_calloc(PIL_FITS_CARD_MAX, sizeof(char));

  if (!key || !value)
      return EXIT_FAILURE;


  fitsfile = newPilFitsFile(pilFrmGetName(frame), PIL_FITS_READWRITE);

  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("Instrument"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("Origin"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("Date"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("DataMD5"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("DprCategory"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("DprType"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("DprTechnique"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("OriginalFile"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("ArchiveFile"));
  pilFitsHdrDelete(fitsfile, pilTrnGetKeyword("Checksum"));


  pilFitsHdrInsertString(fitsfile, 1, "ESO*",
                         pilTrnGetKeyword("Origin"),
                         "ESO",
                         pilTrnGetComment("Origin"));

  date = pilDateGetISO8601();
  if (!date) {
      date = "";
  }

  pilFitsHdrInsertString(fitsfile, 1, "ESO*",
                         pilTrnGetKeyword("Date"),
                         date,
                         pilTrnGetComment("Date"));

  pilFitsHdrInsertString(fitsfile, 1, "ESO*",
                         pilTrnGetKeyword("Instrument"),
                         recipe_info.instrument,
                         pilTrnGetComment("Instrument"));

  pilFitsHdrInsertString(fitsfile, 1, "ESO*",
                         pilTrnGetKeyword("DataMD5"),
                         md5,
                         pilTrnGetComment("DataMD5"));

  snprintf(key, sz, "%s", pilTrnGetKeyword("DataMD5"));

  if (name) {
    char *s = pilFileBaseName(name);

    if (!s) {
      pil_free(key);
      pil_free(value);

      return EXIT_FAILURE;
    }

    pilFitsHdrInsertString(fitsfile, 0, key,
                           pilTrnGetKeyword("ProductFile"),
                           s,
                           pilTrnGetComment("ProductFile"));
  }
  else {
    pilFitsHdrInsertString(fitsfile, 0, key,
                           pilTrnGetKeyword("ProductFile"),
                           pilFrmGetName(frame),
                           pilTrnGetComment("ProductFile"));
  }

  if (pilFitsHdrInsertString(fitsfile, 1, "ESO PRO*", 
                             pilTrnGetKeyword("ProductDID"),
                             PRODUCT_DID,
                             pilTrnGetComment("ProductDID"))
      == EXIT_FAILURE) {
    pilFitsHdrWriteString(fitsfile,
                          pilTrnGetKeyword("ProductDID"),
                          PRODUCT_DID,
                          pilTrnGetComment("ProductDID"));

  }

  snprintf(key, sz, "%s", pilTrnGetKeyword("ProductDID"));
  pilFitsHdrInsertString(fitsfile, 0, key,
                         pilTrnGetKeyword("DoCategory"),
                         pilFrmGetCategory(frame),
                         pilTrnGetComment("DoCategory"));

  switch (pilFrmGetProductType(frame)) {
      case PIL_PRODUCT_TYPE_TEMPORARY:
        type = "TEMPORARY";
        break;

      case PIL_PRODUCT_TYPE_PREPROCESSED:
        type = "PREPROCESSED";
        break;

      case PIL_PRODUCT_TYPE_REDUCED:
        type = "REDUCED";
        break;

      case PIL_PRODUCT_TYPE_QCPARAM:
        type = "QCPARAM";
        break;

      default:
        type = "UNKNOWN";
        break;
  }

  snprintf(key, sz, "%s", pilTrnGetKeyword("DoCategory"));
  pilFitsHdrInsertString(fitsfile, 0, key,
                         pilTrnGetKeyword("ProductType"),
                         type,
                         pilTrnGetComment("ProductType"));

  snprintf(key, sz, "%s", pilTrnGetKeyword("ProductType"));
  pilFitsHdrInsertString(fitsfile, 0, key,
                         pilTrnGetKeyword("RecipeId", 1),
                         recipe_info.name,
                         pilTrnGetComment("RecipeId"));

  snprintf(key, sz, "%s", pilTrnGetKeyword("RecipeId", 1));
  snprintf(value, sz, "%s/%s", recipe_info.instrument, recipe_info.version);
  pilFitsHdrInsertString(fitsfile, 0, key,
                         pilTrnGetKeyword("PipelineId", 1),
                         value,
                         pilTrnGetComment("PipelineId"));
  
  snprintf(key, sz, "%s", pilTrnGetKeyword("PipelineId", 1));
  pilFitsHdrInsertString(fitsfile, 0, key,
                         pilTrnGetKeyword("RecipeStart", 1),
                         pilTimerGetTimeISO8601(recipe_info.start),
                         pilTrnGetComment("RecipeStart"));


  /*
   * Add the names of the raw and calibration frames used to the product
   * header if an input set was given.
   */

  snprintf(key, sz, "%s", pilTrnGetKeyword("RecipeStart", 1));

  if (set) {
    register size_t nr = 0;
    register size_t nc = 0;

    PilFrame *source = pilSofFirst(set);


    while (source) {
      char *s = (char *)pilFrmGetName(source);

      PilFitsFile *cal_file;

      PilFrameType frame_type = pilFrmGetType(source);


      switch (frame_type) {
          case PIL_FRAME_TYPE_RAW:
            ++nr;

            pilFitsHdrInsertString(fitsfile, 0, key,
                                   pilTrnGetKeyword("RawFrameId", 1, nr),
                                   pilFileBaseName(s),
                                   pilTrnGetComment("RawFrameId"));
              
            snprintf(key, sz, "%s", pilTrnGetKeyword("RawFrameId", 1, nr));
            pilFitsHdrInsertString(fitsfile, 0, key,
                                   pilTrnGetKeyword("RawFrameCategory", 1, nr),
                                   pilFrmGetCategory(source),
                                   pilTrnGetComment("RawFrameCategory"));

            snprintf(key, sz, "%s",
                     pilTrnGetKeyword("RawFrameCategory", 1, nr));
            break;

          case PIL_FRAME_TYPE_CALIB:
            ++nc;

            pilFitsHdrInsertString(fitsfile, 0, key,
                                   pilTrnGetKeyword("CalFrameId", 1, nc),
                                   pilFileBaseName(s),
                                   pilTrnGetComment("CalFrameId"));
              
            snprintf(key, sz, "%s", pilTrnGetKeyword("CalFrameId", 1, nc));
            pilFitsHdrInsertString(fitsfile, 0, key,
                                   pilTrnGetKeyword("CalFrameCategory", 1, nc),
                                   pilFrmGetCategory(source),
                                   pilTrnGetComment("CalFrameCategory"));

            snprintf(key, sz, "%s",
                     pilTrnGetKeyword("CalFrameCategory", 1, nc));

            cal_file = newPilFitsFile(s, PIL_FITS_READ);

            if (cal_file) {
              char *cal_md5;

              if (pilFitsHdrReadString(cal_file, pilTrnGetKeyword("DataMD5"),
                                       &cal_md5) == EXIT_SUCCESS) {
                pilFitsHdrInsertString(fitsfile, 0, key,
                                       pilTrnGetKeyword("CalFrameMD5", 1, nc),
                                       cal_md5,
                                       pilTrnGetComment("CalFrameMD5"));
                snprintf(key, sz, "%s",
                         pilTrnGetKeyword("CalFrameMD5", 1, nc));
                pil_free(cal_md5);
              }

              deletePilFitsFile(cal_file);
            }
            break;

          default:
            break;
      }

      source = pilSofNext(set, source);

    }
  }

  deletePilFitsFile(fitsfile);

  pil_free(key);
  pil_free(value);

  return EXIT_SUCCESS;

}
/**@}*/
