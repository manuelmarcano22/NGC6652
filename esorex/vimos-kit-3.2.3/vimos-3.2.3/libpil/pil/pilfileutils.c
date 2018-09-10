/* $Id: pilfileutils.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <wordexp.h>

#include "pilfileutils.h"
#include "pilutils.h"


#ifndef FITS_RECORD_SIZE
#define FITS_RECORD_SIZE  80
#endif

/*
 * According to the ESO Data Interface Control Document
 * (GEN-SPE-ESO-19400-794/1.1/0) the maximum allowed record length for PAF
 * files is 256 characters, but this might change.
 */

#ifndef PAF_RECORD_SIZE
#define PAF_RECORD_SIZE  256
#endif

#ifndef PAF_COMMENT_CHARS
#define PAF_COMMENT_CHARS  "#"
#endif


/**
 * @defgroup pilFile pilFile
 *
 * The module @b pilFile provides utility functions for processing file
 * and directory names respectively. Currently it also provides two
 * utility functions for checking for the type of a file.
 * 
 * @note
 *   The two filetype utility functions provided here may be renamed
 *   and moved to dedicated modules in the future.
 *
 */

/**@{*/

 static int directory_up(char *);


/*
 * @brief
 *   Remove trailing slashes from directory names.
 *
 * @param path  Path name.
 *
 * @return The function returns a pointer to the original but truncated
 *   path name.
 *
 * If @em path has one or more '/' as last character(s), they are
 * replaced by a '@\0' if it is not the filesystem root.
 */

 static char *path_trim(char *path)
{

  char *p = path;

  if (path) {
    p += strlen(path) - 1;
    while (p > path && *p == '/')
      *p-- = '\0';
  }

  return path;

}


/*
 * @brief
 *   Resolve the special tokens '.' and '..' in a path string.
 *
 * @param dest  Buffer receiving the resolved path string.
 * @param src   Source path string.
 *
 * @return The function returns @c EXIT_SUCCESS if the path could be
 *   resolved, otherwise the function returns @c EXIT_FAILURE.
 *
 * The function resolves the special tokens '.' and '..' if they are
 * present in the input path string @em src. Any other combination
 * of dots is considered as a directory name. The input string may not
 * consist of more than @c PIL_PATHNAME_MAX characters. The output buffer
 * must be large enough to receive @c PIL_PATHNAME_MAX characters.
 * 
 * If @em src starts with '.' or '..' the tokens are replaced by the
 * current working directory or its parent directory respectively. If the
 * token '.' appears in place of a directory name within the input path it
 * is simply ignored.
 * 
 * If '..' appears within the input string the last directory component
 * in the output path is removed, i.e. move up one directory level relative
 * to the directory after which the '..' would appear in the output path
 * string.
 * 
 * If there are so many '..' tokens that the resolved path string would go
 * beyond the filesystem root the resolved path is the filesystem root,
 * i.e. '/'.
 * 
 * If the input path does not contain any '/' the input path gets
 * prefixed by the current working directory.
 * 
 * The fully resolved path does not have a trailing slash.
 */

 static int path_resolve(char *dest, const char *src)
{

  size_t sz;

  char *s;
  char spath[PIL_PATHNAME_MAX + 1];
  char tpath[2 * (PIL_PATHNAME_MAX + 1)];  /* Avoid buffer overflow */



  /*
   * Save the input path string to a local buffer since strtok()
   * destroys its input argument.
   */

  if (strlen(src) > PIL_PATHNAME_MAX)
    return EXIT_FAILURE;
  else
    strncpy(spath, src, PIL_PATHNAME_MAX);


  /*
   * Relative and absolute paths are treated differently.
   */

  if (*spath == '/')
    strcpy(tpath, "/");
  else {
    if (!getcwd(tpath, PIL_PATHNAME_MAX))
      return EXIT_FAILURE;

    if (strcmp(tpath, "/"))
      strcat(tpath, "/");
  }


  /*
   * Parse the remaining input string into tokens, resolve the special tokens
   * '.' and '..'. Normal directory names are just appended to the output
   * path string.
   */

  s = strtok(spath, "/");
  sz = strlen(tpath);

  while (s) {
    if (!strcmp(s, "..")) {
      directory_up(tpath);
      sz = strlen(tpath);
    }
    else {
      if (strcmp(s, ".")) {
	sz += strlen(s) + 1;

	if (sz > 2 * PIL_PATHNAME_MAX + 1)
	  return EXIT_FAILURE;
        else {
	  strcat(tpath, s);
	  strcat(tpath, "/");
	}
      }
    }
    s = strtok(NULL, "/");
  }


  /*
   * Remove the trailing slash
   */

  path_trim(tpath);


  /*
   * Check the length of the fully expanded path string and write
   * it to the output buffer.
   */

  if (strlen(tpath) > PIL_PATHNAME_MAX)
    return EXIT_FAILURE;
  else 
    strncpy(dest, tpath, PIL_PATHNAME_MAX);

  return EXIT_SUCCESS;

}


/*
 * @brief
 *   Move one directory up in a path string.
 *
 * @param path  Path string.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function moves up one directory in the input path string by 
 * removing the path component, assuming that it specifies a 
 * subdirectory. The input path may have a trailing slash. The function
 * does not do any kind of expansion. The input path must be terminated.
 * The resulting path string will always have a trailing slash.
 */

 static int directory_up(char *path)
{

  char *s;

  path_trim(path);
  if ((s = strrchr(path, '/')))
    *(s + 1) = '\0';
  else
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}


/**
 * @brief 
 *   Check if the given file complies to the  standard FITS format.
 *
 * @param filename  The name of the file to be checked.
 *
 * @return The function returns 1 if the given file is a standard FITS file,
 *   otherwise the return value is 0.
 *
 * The file format is determined from the first bytes of the file. To
 * be a standard FITS file the first 80 bytes must contain the mandatory
 * FITS keyword SIMPLE and its value must be 'T'.
 * 
 * The function does not check if the file specified by @em filename
 * exists or if the process has sufficient privileges to read it.
 */

extern int pilFileIsFits(const char *filename)
{

  FILE *file;
  register int k, isfits = 0;
  char fits_record[FITS_RECORD_SIZE];


  if ((file = fopen(filename, "r"))) {
    if (fread(fits_record, 1, FITS_RECORD_SIZE, file) == FITS_RECORD_SIZE) {
      if (memcmp(fits_record, "SIMPLE  =", 9) == 0) {
	for (k = 10; fits_record[k] == ' ' && k < FITS_RECORD_SIZE; k++)
	  ;
	if (k < FITS_RECORD_SIZE && fits_record[k] == 'T')
	  isfits = 1;
      }
    }

    fclose(file);

  }

  return isfits;

}


/**
 * @brief 
 *   Check if the given file complies to the VLT parameter file format (PAF).
 *
 * @param filename  The name of the file to be checked.
 *
 * @return The function returns 1 if the given file is a valid PAF file,
 *   otherwise the return value is 0.
 *
 * The file format is determined from the first bytes of the file. To
 * be a VLT parameter file the first line, which is not a comment (first
 * character is a '#'), must start with the string "PAF.HDR.START;".
 * 
 * The function does not check if the file specified by @em filename
 * exists or if the process has sufficient privileges to read it.
 */

extern int pilFileIsPaf(const char *filename)
{

  FILE *file;
  register int ispaf = 0;
  char paf_record[PAF_RECORD_SIZE + 1];


  if ((file = fopen(filename, "r"))) {
    while (fgets(paf_record, PAF_RECORD_SIZE, file))
      if (!strempty(paf_record, PAF_COMMENT_CHARS))
	if (!strncmp(paf_record, "PAF.HDR.START", 13)) {
	  ispaf = 1;
	  break;
	}

    fclose(file);
  }
    
  return ispaf;

}


/**
 * @brief
 *   Remove trailing slash from directory names.
 *
 * @param path  Path name.
 *
 * @return The function returns a pointer to the original but truncated
 *   path name.
 *
 * If @em path has one or more '/' as last character(s) , they are
 * replaced by a '@\0'. The slash is not replaced if the input path is a
 * single '/', i.e. the filesystem root. The function assumes the 
 * @em path points to a valid path string that is properly
 * terminated. Note that the function modifies its input.
 */

char *pilFileTrimPath(char *path)
{

  return path_trim(path);

}


/*
 * @brief
 *   Strip directory and suffix from a filename.
 *
 * @param filename  Path to the file.
 *
 * @return Pointer to the stripped filename, or @c NULL in case an error
 *    occurred.
 *
 * The function removes any leading directory components from
 * @em filename. The input string must not contain more than
 * @c PIL_PATHNAME_MAX characters.
 * 
 * The result is stored in a static, internal buffer, i.e. it is
 * overwritten by the next call to @b pilFileBaseName().
 */

 char *pilFileBaseName(const char *filename)
{

  static char basename[PIL_PATHNAME_MAX + 1];

  char *s;
  char tfile[PIL_PATHNAME_MAX + 1];


  /*
   * Create a working copy
   */

  if (strlen(filename) > PIL_PATHNAME_MAX)
    return NULL;
  else 
    strncpy(tfile, filename, PIL_PATHNAME_MAX);


  /*
   * Remove trailing slashes.
   */

  path_trim(tfile);


  /*
   * Strip leading directories
   */

  if ((s = strrchr(tfile, '/')))
    strncpy(basename, s + 1, PIL_PATHNAME_MAX);
  else
    strncpy(basename, tfile, PIL_PATHNAME_MAX);

  return basename;
}


/*
 * @brief
 *   Strip non-directory suffix from a filename.
 *
 * @param filename  Path to the file.
 *
 * @return Pointer to the stripped filename, or @c NULL in case an error
 *    occurred.
 *
 * The function removes the trailing directory component from
 * @em filename. If filename does not contain any slashes the
 * result is the current working directory. The input string must
 *  not contain more than @c PIL_PATHNAME_MAX characters.
 * 
 * The result is stored in a static, internal buffer, i.e. it is
 * overwritten by the next call to @b pilFileBaseName().
 */

 char *pilFileDirName(const char *filename)
{

  static char dirname[PIL_PATHNAME_MAX + 1];

  size_t sz;

  char *s;


  /*
   * Is there any slash?
   */

  if (!(s = strchr(filename, '/'))) {
    if (!getcwd(dirname, PIL_PATHNAME_MAX))
      return NULL;
  }
  else {
    if ((sz = s - filename) > PIL_PATHNAME_MAX)
      return NULL;
    else {
      strncpy(dirname, filename, sz);
      dirname[sz + 1] = '\0';
    }
  }

  return dirname;

}


/*
 * @brief
 *   Do a complete expansion of a directory path string.
 *
 * @param path  Directory path to be expanded.
 *
 * @return The function returns a pointer to the expanded path string if
 *   expanding the input string was successful. Otherwise the function
 *   returns @c NULL.
 *
 * The function does an expansion of the directory path string
 * @em path. The result is @em path converted into an absolute
 * file path.
 * 
 * If @em path is a relative path, i.e. it starts with '.' or '..'
 * these tokens are replaced by the absolute path to the current
 * working directory or its parent directory respectively.
 * 
 * If @em path contains the tilde character (~) as the first
 * character it is replaced by the users home directory if the second
 * character is a '/'. If the tilde is followed by a valid user name
 * it is replaced by the path to this users's home directory.
 * 
 * If @em path contains the reference to an environment variable,
 * the environment variable is replaced by its contents. If the environment
 * variable does not exist the function returns an error.
 */

 char *pilFileExpandDirPath(const char *path)
{

  static char expanded_path[PIL_PATHNAME_MAX + 1];

  char tpath[PIL_PATHNAME_MAX + 1];

  wordexp_t pwordexp = {0, 0, 0};


  /*
   * Initialize the output path
   */

  *expanded_path = '\0';

  /*
   * Perform tilde, variable expansion. The expanded path must
   * be unique.
   */

  if (wordexp(path, &pwordexp, WRDE_NOCMD | WRDE_UNDEF) || 
      pwordexp.we_wordc > 1) {
    if (pwordexp.we_wordc > 0)
      wordfree(&pwordexp);
    return NULL;
  }
  else {
    if (!pwordexp.we_wordv[0] ||
        strlen(pwordexp.we_wordv[0]) > PIL_PATHNAME_MAX) {
      wordfree(&pwordexp);
      return NULL;
    }
    strncpy(tpath, pwordexp.we_wordv[0], PIL_PATHNAME_MAX);
    wordfree(&pwordexp);
  }
    
  /*
   * Convert the path string into an absolute path resolving the
   * special tokens '.' and '..'.
   */

  if (path_resolve(expanded_path, tpath) == EXIT_FAILURE)
    return NULL;

  return expanded_path;
  
}


/*
 * @brief
 *   Do a complete expansion of a file path string.
 *
 * @param filepath  File path to be expanded.
 *
 * @return The function returns a pointer to the expanded path string if
 *   expanding the input string was successful. Otherwise the function
 *   returns @c NULL.
 *
 * The function does an expansion of the directory part of the path string
 * @em filepath. The result is @em filepath converted into an absolute file
 * path. The last component of @em filepath is considered as the name of
 * the file.
 * 
 * To do the expansion the function uses @b pilFileExpandDirPath()
 * internally. The input file path must not have a trailing slash. A slash
 * as the last character of the input string causes an error.
 *
 * @see pilFileExpandDirPath()
 */

 char *pilFileExpandFilePath(const char *filepath)
{

  if (*(filepath + strlen(filepath) - 1) == '/')
    return NULL;

  return pilFileExpandDirPath(filepath);

}
/**@}*/
