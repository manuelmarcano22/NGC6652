/*
 * This file is part of the ESO Recipe Execution Tool
 * Copyright (C) 2017 European Southern Observatory
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

/* Author: Artur Szostak <aszostak@partner.eso.org> */

/*
 * This unit test will simulate execution of a very simple Python based recipe.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef ENABLE_PYTHON_RECIPES

#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include <string.h>
#include <cpl.h>
#include <cpl_test.h>

#ifndef TEST_RECIPE_DIR
#define TEST_RECIPE_DIR  "."
#endif

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void write_file(const char * filename, const char * text);
static char * read_file(const char * filename);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    if (system("python --version > /dev/null") != 0)
    {
        /* Indicate this unit test is skipped if Python is not available. */
        return 77;
    }

    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    cpl_test_assert(argc >= 1);
    char * argv0 = cpl_strdup(argv[0]);
    const char * dir = dirname(argv0);

    /* Create a dummy input file that contains the integer 12. */
    char * filename = cpl_sprintf("%s/python_recipe_test_input.fits", dir);
    write_file(filename, "12\n");
    cpl_free(filename);

    /* Create the Set-of-Frames file. */
    filename = cpl_sprintf("%s/python_recipe_test.sof", dir);
    char * value = cpl_sprintf("%s/python_recipe_test_input.fits RAW\n", dir);
    write_file(filename, value);
    cpl_free(value);
    cpl_free(filename);

    /* Prepare the command to execute esorex  */
    char * esorex_command = cpl_sprintf(
            "%s/../esorex --suppress-prefix=TRUE --recipe-dir=%s"
            " testrecipe --par1=7 %s/python_recipe_test.sof"
            " > %s/python_recipe_test_recipe_output.log",
            dir, TEST_RECIPE_DIR, dir, dir
        );
    cpl_test_eq(system(esorex_command), 0);

    /* Read the console output from EsoRex and the output file produced by the
       recipe. */
    filename = cpl_sprintf("%s/python_recipe_test_recipe_output.log", dir);
    char * console_output = read_file(filename);
    cpl_test_assert(console_output != NULL);
    cpl_free(filename);
    filename = cpl_sprintf("%s/python_recipe_test_output.fits", dir);
    char * output = read_file(filename);
    cpl_test_assert(output != NULL);
    cpl_free(filename);

    /* Check that the output is correct, i.e. input + par1 = 12 + 7 = 19 */
    cpl_test_eq_string(output, "19\n");

    /* Check that EsoRex registered the output file. */
    cpl_test(strstr(console_output, "1 product created") != NULL);

    cpl_free(output);
    cpl_free(console_output);
    cpl_free(esorex_command);
    cpl_free(argv0);
    return cpl_test_end(0);
}


/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

/*
 * This function will write the given text to a file.
 */
static void write_file(const char * filename, const char * text)
{
    FILE* file = fopen(filename, "w");
    cpl_test_assert(file != NULL);
    size_t bytes_to_write = strlen(text);
    size_t bytes_written = fwrite(text, sizeof(char), bytes_to_write, file);
    cpl_test_assert(bytes_written == bytes_to_write);
    int fclose_result = fclose(file);
    cpl_test_assert(fclose_result == 0);
}

/*
 * Read a text file and return the contents as a string.
 * The caller must free the returned string with cpl_free().
 */
static char * read_file(const char * filename)
{
    size_t buffersize = 1024*4;
    char * buffer = cpl_malloc(buffersize);
    FILE* file = fopen(filename, "r");
    cpl_test_assert(file != NULL);
    size_t bytes_read = fread(buffer, sizeof(char), buffersize-1, file);
    int fclose_result = fclose(file);
    cpl_test_assert(fclose_result == 0);
    buffer[bytes_read] = '\0';   /* Make sure to null terminate the string. */
    return buffer;
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
