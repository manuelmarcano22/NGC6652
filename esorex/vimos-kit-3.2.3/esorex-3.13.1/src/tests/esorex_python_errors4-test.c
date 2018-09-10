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
 * These are tests relying on overloading certain functions to mock various
 * failure modes that are difficult to trigger in any other way.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef ENABLE_PYTHON_RECIPES

/* Undefine NDEBUG to make sure assert() macros are compiled in for these
   tests. */
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <errno.h>
#include <cpl.h>
#include <cpl_test.h>

/*-----------------------------------------------------------------------------
                            Overloaded functions
 -----------------------------------------------------------------------------*/

/* Overload the following function to produce a desired error. */

int errno_value_to_use = EIO;

int _overloaded_execlp(const char *file, const char *arg, ...)
{
    (void) file;
    (void) arg;
    errno = errno_value_to_use;
    return -1;
}

/* Replace calls in the code we want to test with our overloaded functions. */
#define execlp _overloaded_execlp
#include "er_python.c"
#undef execlp

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void test_start_python_interpreter(void);
static void test_run_python_command(void);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Insert tests below */
    test_start_python_interpreter();
    test_run_python_command();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

static void test_start_python_interpreter(void)
{
    /* Test edge case where execlp() fails. In this case the function
       start_python_interpreter should not be able to detect the error since
       it happens in the child process. */
    pid_t pid = -1;
    int input = -1;
    int output = -1;
    cpl_test_eq_error(start_python_interpreter("x=1", &pid, &input, &output),
                      CPL_ERROR_NONE);
}

static void test_run_python_command(void)
{
    /* Test edge case where execlp() fails. The run_python_command function
       should detect the error. */
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
