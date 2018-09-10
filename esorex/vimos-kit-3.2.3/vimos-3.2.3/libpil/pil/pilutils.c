/* $Id: pilutils.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2004 European Southern Observatory
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

#include <ctype.h>
#include <string.h>

/**
 * @defgroup pilUtils pilUtils
 *
 * The module @b pilUtils provides functions related to string manipulation.
 */

/**@{*/

/**
 * @brief 
 *   Convert all uppercase characters in a string into lowercase
 *   characters.
 *
 * @return Returns a pointer to the converted string.
 *
 * @param s  Pointer to the string to be converted.
 *
 * Walks through the given string and turns uppercase characters into
 * lowercase characters using @b tolower().
 *
 * @see strupper
 *
 * @author R. Palsa
 */

 char *strlower(char *s)
{

  char *t = s;

  while (*t) {
    *t = tolower(*t);
    t++;
  }

  return s;

}


/**
 * @brief 
 *   Convert all lowercase characters in a string into uppercase
 *   characters.
 *
 * @return Returns a pointer to the converted string.
 *
 * @param s  Pointer to the string to be converted.
 *
 * Walks through the given string and turns lowercase characters into
 * uppercase characters using @b toupper().
 *
 * @see strlower
 *
 * @author R. Palsa
 */

 char *strupper(char *s)
{

  char *t = s;

  while (*t) {
    *t = toupper(*t);
    t++;
  }

  return s;

}


/**
 * @brief
 *   Remove leading and/or trailing whitespace characters from a string.
 *
 * @return The function returns a pointer to the modified string if no
 *   error occurred, otherwise @c NULL.
 *
 * @param s     String to be processed.
 * @param mode  Processing option.
 *
 * The function removes leading or trailing whitespace characters, or 
 * both from the string @em s, according to the processing option
 * @em mode. Whitespace characters are recognized by @b isspace().
 * 
 * The following processing options are valid:
 *   @li 0: Remove leading whitespace characters.
 *   @li 1: Remove trailing whitespace characters.
 *   @li 2: Remove leading and trailing whitespace characters.
 *
 * @author R. Palsa
 */

 char *strtrim(char *s, unsigned int mode)
{

  register char *t, *r;

  if (mode > 2)
    return NULL;

  /*
   * Note: The characer class functions have to be enclosed in
   *   parantheses to use the actual function on HP-UX where these
   *   functions are also provided as macros, which are taken by
   *   default and may lead to compiler warnings.
   */

  if (mode == 1 || mode == 2) {
    t = s + strlen(s) - 1;
    while ((isspace)((unsigned char)*t))
      t--;
    *++t = '\0';
  }
    
  if (mode == 0 || mode == 2) {
    t = r = s;
    while ((isspace)((unsigned char)*t))
      t++;
    while ((*r++ = *t++))
      ;
  }

  return s;

}


/**
 * @brief
 *   Locate first character in a string that does not belong to the
 *   specified character class.
 *
 * @return Pointer to the first character that is not a member of
 *   the given character class.
 *
 * @param s      String to be processed.
 * @param ctype  Function used to test if a character belongs to a particular
 *               character class.
 *
 * Searches the string @em s for the first occurence of a character
 * which does not belong to a certain character class. The character
 * class is represented through a function that returns a non zero vlaue
 * if a character belongs to that class and 0 otherwise. Such functions
 * are the character classification routines like @b isspace() for
 * instance. It is expected that the input string is properly terminated.
 * In case the whole string consists of characters of the specified class
 * the function will return the location of the terminating '@\0'.
 *
 * @author R. Palsa
 */

 char *strskip(const char *s, int (*ctype)(int))
{

  register char *t = (char *)s;

  while (ctype((unsigned char)*t))
    t++;

  return t;

}


/**
 * @brief
 *   Split a string according to a max allowed width.
 *
 * @return Pointer to the modified character string, or @c NULL in case of
 *         failure. The modified character string is statically allocated
 *         within this function and it cannot exceed 1024 characters.
 *
 * @param s       String to be processed.
 * @param blanks  Number of blanks to be inserted at every split point.
 * @param width   Max number of characters between split points.
 *
 * This function is typically used in splitting a string
 * avoiding to exceed a given maximum width (as for instance 
 * the width of a terminal where the string is going to be
 * printed). The splitting is performed preserving words, 
 * i.e. by replacing with a newline character ('@\n') the last 
 * blank character before the maximum allowed width. Newline 
 * characters already present in the input string are preserved. 
 * Single words that are long enough to exceed the max allowed 
 * width would not be split, just in this case long lines are 
 * tolerated. A number of blanks to be inserted at every split 
 * point must be specified, setting the left indentation level 
 * for the printed string. This number must not exceed the
 * maximum allowed width.
 *
 * @author C. Izzo
 */

#define MAX_LENGTH (1024)

char *strsplit(char *s, unsigned int blanks, unsigned int width)
{

  unsigned int  i, j, k;
  unsigned int  cuti  = 0;
  unsigned int  cutj  = 0;
  unsigned int  limit = width;

  static char split[MAX_LENGTH];

  if (blanks >= width) {
   /*
    * Give up indentation
    */
    blanks = 0;
  }

  for (i = 0, j = 0; i < MAX_LENGTH && j < MAX_LENGTH; i++, j++) {

    split[j] = s[i];

    if (s[i] == ' ' || s[i] == '\0' || s[i] == '\n') {

      if (i > limit) {

       /*
        *  Go back to the previous cuttable position, if possible
        */

        if (limit - cuti < width - blanks) {
          j = cutj;
          i = cuti;
        }
        else {
          if (s[i] == '\0')
            break;
        }

       /*
        *  Split here, and insert blanks
        */

        split[j] = '\n';

        for (k = 0, j++; k < blanks && j < MAX_LENGTH; k++, j++) {
          split[j] = ' ';
        }
        j--;

        limit += width - blanks - (limit - i);
      }
      else {
        if (s[i] == '\0') break;
        if (s[i] == '\n') {

         /*
          *  Split point already present in input string: just add
          *  the specified number of blanks
          */

          i++;
          if (s[i] == '\0') {
            split[j] = '\0';
            break;
          }

          for (k = 0, j++; k < blanks && j < MAX_LENGTH; k++, j++) {
            split[j] = ' ';
          }
          j--;

          limit += width - blanks - (limit - i);
        }

       /*
        *  Keep track of the last cuttable position
        */

        cutj = j;
        cuti = i;

      }
    }
  }


 /*
  *  Safety belt!
  */

  split[MAX_LENGTH - 1] = '\0';

  return split;

}

#undef MAX_LENGTH


/**
 * @brief
 *   Test if a string represents an empty or comment line.
 *
 * @return The function returns 1 if the string is found to be empty or if
 *   the first non--whitespace character is one out of the set of provided
 *   comment characters. Otherwise the function returns 0.
 *
 * @param s  String to be tested.
 * @param p  String containing all allowed comment characters.
 *
 * The function skips all leading whitespace characters in the string
 * @em s. Whitespace characters are recognized by @b isspace().
 * If the first character which is not a whitespace character is either
 * '@\0' or one out of the pattern string @em p, the string is
 * considered as empty and the function returns 1.
 * 
 * If @em p is set to @c NULL there is no checking for special
 * characters that should be considered as whitespaces.
 *
 * @author R. Palsa
 */

 int strempty(const char *s, const char *p)
{

  register char *t = strskip(s, isspace);

  if (*t) {
    if (!p || !strchr(p, *t))
      return 0;
  }

  return 1;

}


/**
 * @brief
 *   Select a string pattern from a list.
 *
 * @return The function returns the index of the pattern string if the
 *   string is found in the list, otherwise the return value is -1.
 *
 * @param pattern  Pattern string.
 * @param list     List of strings to be searched.
 * @param sz       Number of elements the list contains.
 *
 * The function sequentially walks through the list until the pattern
 * @em pattern is found in @em list or the size @em sz of the list is
 * reached. The function internally uses @b strcmp() to find the pattern
 * in the list.
 *
 * @author R. Palsa
 */

size_t strselect(const char *pattern, const char **list, size_t sz)
{

  size_t i = 0;

  while (i < sz) {
    if (!strcmp(pattern, list[i]))
      return i;
    i++;
  }

  return -1;

}
/**@}*/
