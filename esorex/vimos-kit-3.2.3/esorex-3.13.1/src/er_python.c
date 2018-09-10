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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <time.h>
#include <sys/wait.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <libgen.h>
#include <cpl.h>
#include <cxmap.h>
#include "er_json.h"
#include "er_python.h"

/**
 * @defgroup esorex_python Python Interfacing Functions
 *
 * This module provides functions needed to interface with a Python interpreter
 * running as an external process that will execute a Python based plugin.
 * The CPL plugin API is preserved as much as possible, both on the EsoRex side
 * and within the Python plugin code. Thus, the EsoRex core code will invoke the
 * Python based plugin through the usual calls to the @c cpl_plugin class (more
 * specifically, derived classes such as @c cpl_recipe).
 * This module provides the mechanism to load Python modules, identify the
 * Python based plugins, select the desired module, and then actually perform
 * the necessary communication with the Python interpreter to execute the
 * plugin.
 *
 * The exported functions of this interface deal with the differences in loading
 * Python module plugins versus C/C++ based shared library plugins. Only the
 * loading and unloading needs special treatment from EsoRex, since the actual
 * invocation of the plugin is through the usual @c cpl_plugin API. Thus,
 * starting/stopping the Python interpreter and the communication protocol is
 * all hidden behind the interface from EsoRex's point of view inside the
 * private functions of this module.
 *
 * Loading of a Python module involves starting the Python interpreter, loading
 * the named Python module file and finding all Python classes that derive from
 * a class called @c CplPlugin. This is accomplished with the function
 * er_python_load_modules(). All such classes are registered and can be
 * subsequently selected with er_python_select_module(). Selecting the desired
 * module is necessary, so that er_python_get_plugin_list() knows for which
 * Python module to return the list of available plugins. The function
 * er_python_get_plugin_list() can be considered as the equivalent of the
 * @c cpl_plugin_get_info() entry point function that must be implemented by a
 * C/C++ based shared library plugin.
 * Once er_python_get_plugin_list() has successfully returned a list of CPL
 * plugin objects, there is no difference in the invocation procedure to execute
 * the plugin between a Python module or shared library based plugin.
 * Once the invocation of the plugin is complete, one must cleanup by calling
 * the er_python_cleanup() function.
 * A simple code example for the call sequence to invoke a Python based plugin
 * will look as follows (error handling is ignored for simplicity):
 * @code
 * er_stringarray_t * modules = er_stringarray_new();
 * cpl_pluginlist * list = cpl_pluginlist_new();
 * cpl_plugin_func plugin_init, plugin_exec, plugin_deinit;
 *
 * // Python specific module initialisation.
 * er_stringarray_append(modules, "./recipe.py");
 * er_python_load_modules(modules);
 * er_python_select_module("./recipe.py");
 * er_python_get_plugin_list(list);
 *
 * // Common plugin invocation sequence.
 * cpl_recipe * recipe = (cpl_recipe *) cpl_calloc(1, sizeof(cpl_recipe));
 * recipe->frames = cpl_frameset_new();
 * plugin_init = cpl_plugin_get_init(&recipe->interface);
 * if (plugin_init != NULL) {
 *   plugin_init(&recipe->interface);
 * }
 * plugin_exec = cpl_plugin_get_exec(&recipe->interface);
 * if (plugin_exec != NULL) {
 *   plugin_exec(&recipe->interface);
 * }
 * plugin_deinit = cpl_plugin_get_deinit(&recipe->interface);
 * if (plugin_deinit != NULL) {
 *   plugin_deinit(&recipe->interface);
 * }
 * cpl_frameset_delete(recipe->frames);
 * cpl_plugin_delete(&recipe->interface);
 *
 * // Python specific cleanup.
 * er_python_cleanup();
 * @endcode
 *
 * Under the hood, communication with the Python interpreter is handled by the
 * @c execute function, which is registered with the @c cpl_plugin objects
 * returned by er_python_get_plugin_list(). The communication protocol uses JSON
 * as the exchange format that is sent over Unix pipes. The @c execute function
 * is a stub that will implement the following steps:
 * <ol>
 *   <li>
 *      Start an external Python interpreter process with the @c python shell
 *      command and give it the initial Python stub code to execute.
 *   </li>
 *   <li>
 *      Encode the input @c cpl_plugin structure as JSON format.
 *   </li>
 *   <li>
 *      Send the JSON text to the Python interpreter over a Unix pipe.
 *   </li>
 *   <li>
 *      Close the input pipe to indicate end of input.
 *   </li>
 *   <li>
 *      Read the output from a second pipe, to which the Python process will
 *      write the results in JSON format.
 *   </li>
 *   <li>
 *      Wait for the output pipe to be closed, which indicates end of output.
 *   </li>
 *   <li>
 *      Decode the output JSON text back into a @c cpl_plugin structure and
 *      update the @c cpl_plugin argument of the @c execute function.
 *   </li>
 *   <li>
 *      Wait for the Python interpreter to terminate (join with the process).
 *   </li>
 *   <li>
 *      Return the plugin's return code received from the Python plugin.
 *   </li>
 * </ol>
 * The other end of the protocol is implemented in the initial Python stub code
 * passed to the Python interpreter to execute. The steps it will perform are:
 * <ol>
 *   <li>
 *      Import the desired Python module containing the plugin.
 *   </li>
 *   <li>
 *      Instantiate the selected plugin object (a class that derives from
 *      @c CplPlugin).
 *   </li>
 *   <li>
 *      Read the JSON input from the Unix pipe.
 *   </li>
 *   <li>
 *      Wait for the input pipe to be closed, which indicates the end of input.
 *   </li>
 *   <li>
 *      Decode the JSON into a Python object.
 *   </li>
 *   <li>
 *      Call the plugin object's @c execute method, passing the decoded Python
 *      object as the method's argument.
 *   </li>
 *   <li>
 *      Once the Python @c execute method completes, the updated argument and
 *      return value are encoded as JSON again.
 *   </li>
 *   <li>
 *      The encoded JSON result is send over the output pipe.
 *   </li>
 *   <li>
 *      The output pipe is closed to indicate the end of output.
 *   </li>
 *   <li>
 *      The Python interpreter ends execution and terminates.
 *   </li>
 * </ol>
 * Note that the JSON is not send over the stdin and stdout pipes, but rather
 * two new dedicated pipes are created to handle the interprocess communication.
 * This allows the stdin and stdout pipes to behave as one would normally
 * expect.
 * For example, they will be attached to the same terminal or file as the parent
 * process (i.e. EsoRex) was attached to. This means that print statements in
 * the Python code will write the output to the same location as one would
 * expect if the Python code was executed manually in the Python interpreter.
 * Separation of the EsoRex to Python communication protocol from stdin/stdout
 * avoids unexpected interference or having stdin/stdout redirected in an
 * unexpected manner.
 *
 * Developers of Python based plugins must implement a plugin class that has the
 * following requirements:
 * <ul>
 *   <li>
 *      Must derive from a class called @c CplPlugin or it must itself be called
 *      @c CplPlugin.
 *   </li>
 *   <li>
 *      The @c \__init__ constructor must take no arguments or all arguments
 *      must have defaults.
 *   </li>
 *   <li>
 *      The following class attributes must be defined (unless indicated as
 *      optional) or set in the @c \__init__ function:
 *     <dl>
 *       <dt>name</dt>
 *       <dd>Name of the plugin/recipe. If not given then the Python class's
 *          name is used instead. (type @c basestring, optional)</dd>
 *       <dt>version</dt>
 *       <dd>The version number of the plugin. (type @c int)</dd>
 *       <dt>synopsis</dt>
 *       <dd>A summary text giving a brief description of the plugin. If not
 *          given then the first line of the doc string for the class is used
 *          instead. (type @c basestring, optional)</dd>
 *       </li>
 *       <dt>description</dt>
 *       <dd>A description text for the plugin. If not given then the doc string
 *          for the class is used instead. (type @c basestring, optional)</dd>
 *       </li>
 *       <dt>author</dt>
 *       <dd>The author of the plugin. An empty string is used if not given.
 *          (type @c basestring, optional)</dd>
 *       <dt>email</dt>
 *       <dd>The email address of plugin's author. An empty string is used if
 *          not given. (type @c basestring, optional)</dd>
 *       <dt>copyright</dt>
 *       <dd>The copyright notice for the plugin. An empty string is used if
 *          not given. (type @c basestring, optional)</dd>
 *       <dt>parameters</dt>
 *       <dd>A list of input parameters accepted by the plugin (type @c list,
 *          optional). This must be a list of dictionaries with the following
 *          keys:
 *          </dd>
 *       <dl>
 *         <dt>class</dt>
 *         <dd>The type of parameter. Must be one of @c 'value', @c 'range' or
 *          @c 'enum'. (type @c basestring)</dd>
 *         <dt>name</dt>
 *         <dd>The name of the parameter. (type @c basestring)</dd>
 *         <dt>default</dt>
 *         <dd>The default value of the parameter if no value is assigned.
 *          (type can be one of @c bool, @c int, @c float or @c basestring)</dd>
 *         <dt>min</dt>
 *         <dd>The minimum allowed value of the range. This only needs to be set
 *          if @a class is @c 'range'. (type @c int or @c float, but @a default
 *          must have the same matching type)</dd>
 *         <dt>max</dt>
 *         <dd>The maximum allowed value of the range. This only needs to be set
 *          if @a class is @c 'range'. (type @c int or @c float, but @a default
 *          must have the same matching type)</dd>
 *         <dt>choices</dt>
 *         <dd>The list of allowed values. This only needs to be set if @a class
 *          is @c 'enum'. The individual elements in the list can be one of
 *          @c int, @c float or @c basestring; but must all have the same type
 *          and must match the type used for @a default. (type @c list)</dd>
 *         <dt>description</dt>
 *         <dd>Text describing the parameter. (type @c basestring,
 *          optional)</dd>
 *         <dt>context</dt>
 *         <dd>The application context for the parameter, typically the plugin
 *          name. (type @c basestring, optional)</dd>
 *         <dt>tag</dt>
 *         <dd>A free form tag associated with a parameter. (type @c basestring,
 *          optional)</dd>
 *         <dt>cli_enabled</dt>
 *         <dd>Flag indicating if the parameter is available on the EsoRex
 *          command line. (type @c bool, optional)</dd>
 *         <dt>cli_alias</dt>
 *         <dd>The name of the parameter when given on the EsoRex command line.
 *          If not set then the full parameter name is used.
 *          (type @c basestring, optional)</dd>
 *         <dt>env_enabled</dt>
 *         <dd>Flag indicating if the parameter can be set as an environment
 *          variable. (type @c bool, optional)</dd>
 *         <dt>env_alias</dt>
 *         <dd>The name of the environment variable to use as this parameter.
 *          If not set then the full parameter name is used as the environment
 *          variable. (type @c basestring, optional)</dd>
 *         <dt>cfg_enabled</dt>
 *         <dd>Flag indicating if the parameter can be set from a configuration
 *          file. (type @c bool, optional)</dd>
 *         <dt>cfg_alias</dt>
 *         <dd>The name of the parameter to use in the configuration file.
 *          If not set then the full parameter name is used.
 *          (type @c basestring, optional)</dd>
 *       </dl>
 *       <dt>recipeconfig</dt>
 *       <dd>A list of frame tags and their associations to input and output
 *          frames. (type @c list, optional). This must be a list of
 *          dictionaries with the following keys:
 *          </dd>
 *       <dl>
 *         <dt>tag</dt>
 *         <dd>The frame's tag value, typically the tag used in the
 *          Set-of-Frames. (type @c basestring)</dd>
 *         <dt>min</dt>
 *         <dd>The minimum number of frames of type @a tag. (type @c int)</dd>
 *         <dt>max</dt>
 *         <dd>The maximum number of frames of type @a tag. (type @c int)</dd>
 *         <dt>inputs</dt>
 *         <dd>The input frames associated with the frame of type @a tag
 *          (type @c list). This must be a list of dictionaries with the
 *          following keys:</dd>
 *         <dl>
 *           <dt>tag</dt>
 *           <dd>The input frame's tag value. (type @c basestring)</dd>
 *           <dt>min</dt>
 *           <dd>The minimum number of associated input frames.
 *               (type @c int)</dd>
 *           <dt>max</dt>
 *           <dd>The maximum number of associated input frames.
 *               (type @c int)</dd>
 *         </dl>
 *         <dt>outputs</dt>
 *         <dd>The output frames associated with the frame of type @a tag.
 *          (type @c list of @c basestring elements)</dd>
 *       </dl>
 *     </dl>
 *   </li>
 *   <li>
 *      Must have a method called @c execute defined as follows:
 *
 *          def execute(self, plugin):
 *
 *      It will take a single dictionary argument @a plugin and return an
 *      integer return code.
 *      @c 0 will indicate success and any other value will indicate an error.
 *      If an error does occur, the @c execute method should assign an
 *      appropriate string message to the attribute @c self.error_message.
 *
 *      Input parameters can be accessed as @c plugin['parameters']. Input
 *      frames are accessed as @c plugin['frames']. Any new output frames must
 *      be appended to the @c plugin['frames'] list.
 *   </li>
 * </ul>
 *
 * The following is a simple example of a Python class that implements the
 * required interface:
 * @code{.py}
 *  class CplPlugin(object):
 *      name = "test"
 *      version = 102030
 *      synopsis = "short description"
 *      description = "long description"
 *      author = "Some Author"
 *      email = "someone@example.org"
 *      copyright = "copyright 2017"
 *      parameters = [
 *          {
 *              'class': 'value',
 *              'name': 'test.par1',
 *              'description': 'parameter 1',
 *              'context': 'test',
 *              'default': 3
 *          }
 *      ]
 *
 *      def execute(self, plugin):
 *          # Fetch the parameter value.
 *          parameters = plugin['parameters']
 *          par1 = None
 *          for par in parameters:
 *              if par['name'] == 'par1':
 *                  par1 = par['value']
 *
 *          # Fetch the first frame's file name.
 *          frames = plugin['frames']
 *          filename = frames[0]['filename']
 *
 *          try:
 *              # ... process frames ...
 *          except:
 *              # Set the error message and indicate an error.
 *              self.error_message = "Recipe failed"
 *              return 1
 *
 *          # Append new frames as output.
 *          frames.append(
 *              {
 *                  'filename': 'output.fits',
 *                  'tag': 'PROD',
 *                  'type': 2,
 *                  'group': 3,
 *                  'level': 3,
 *              }
 *          )
 *
 *          # Indicate success.
 *          return 0
 * @endcode
 */

/**@{*/

/* The following are file descriptor numbers for the pipes that the Python
   process should use to communicate with EsoRex. These will be prepared just
   before executing the Python interpreter binary. */
#define PYTHON_INPUT_FD   3
#define PYTHON_OUTPUT_FD  4


/* The following is a global map of Python module file names to cpl_pluginlist
   objects created by the er_python_load_modules function. */
static cx_map * module_list = NULL;

/* Identifies the current module to be handled by er_python_get_plugin_list.
   This is set with er_python_select_module. */
static cx_map_iterator current_module = NULL;

/* The following is a map between cpl_plugin objects found within the plugin
   lists in module_list, and the corresponding python module names. In effect,
   this is the inverse lookup of module_list. */
static cx_map * module_lut = NULL;

/* The following command is passed to the -c option to the Python interpreter.
   It bootstraps the more complex script in the JSON input from standard input.
   This keeps the command line length reasonably short. */
static const char * python_start_command =
    "import sys, os, json;"
    " _data = json.load(os.fdopen("CPL_STRINGIFY(PYTHON_INPUT_FD)",'r'));"
    " exec(_data['script'])";


/*
 * Attempts to close a file descriptor and prints an error if an error occurs
 * that is not a signal interrupt. The errno is reset to make sure we ignore
 * any errors from close, since there is nothing more we can do about them.
 */
static void close_file(int file)
{
    if (close(file) == 0) return;
    if (errno != EINTR)
    {
        cpl_msg_error(cpl_func,
                      "Failed to close file descriptor for pipe to Python"
                      " interpreter: %s", strerror(errno));
    }
    errno = 0;
}

/*
 * Starts a new Python interpreter process.
 *
 * @param[in]  python_command  The Python command string to pass to the
 *                             interpreter with the '-c' option for evaluation.
 * @param[out] python_pid  A pointer where the new process ID for the Python
 *                         interpreter process will be written.
 * @param[out] python_input  A pointer where the input file descriptor is
 *                           written. This descriptor is for a pipe to which
 *                           input data should be written to the Python process.
 * @param[out] python_output A pointer where the output file descriptor is
 *                           written. This descriptor is for a pipe on which
 *                           output data should be read from the Python process.
 *
 * @return  @c CPL_ERROR_NONE is returned when the new Python interpreter
 *      process has been started. Otherwise an appropriate error code is
 *      returned.
 */
static cpl_error_code start_python_interpreter(const char * python_command,
                                               pid_t * python_pid,
                                               int * python_input,
                                               int * python_output)
{
    assert(python_command != NULL);
    assert(python_pid != NULL);
    assert(python_input != NULL);
    assert(python_output != NULL);

    int rpipes[2];
    int wpipes[2];

    /* Prepare two streams that will be connected to the Python process for
       input and output. */
    if (pipe(rpipes) != 0)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                     "Failed to create read pipe for Python"
                                     " interpreter: %s", strerror(errno));
    }
    if (pipe(wpipes) != 0)
    {
        close_file(rpipes[0]);
        close_file(rpipes[1]);
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                     "Failed to create write pipe for Python"
                                     " interpreter: %s", strerror(errno));
    }

    *python_pid = fork();
    if (*python_pid == -1)
    {
        close_file(rpipes[0]);
        close_file(rpipes[1]);
        close_file(wpipes[0]);
        close_file(wpipes[1]);
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                     "Failed to start Python process: %s",
                                     strerror(errno));
    }
    if (*python_pid == 0)
    {
        /* In child process after the fork. Here we reconnect the read/write
           pipes we created to use file descriptors that the Python process will
           expect to communicate with EsoRex on. Once the pipes are ready we
           then call exec to start executing the Python interpreter in the child
           process. NOTE: before actually using dup2 to set the required file
           descriptor numbers, we need to duplicate the two descriptors that
           will be redirected. This is important to avoid really closing the
           underlying pipe by mistake in the dup2 call if the required file
           descriptor number overlaps with one of the existing pipe descriptors.
           Having a duplicate will avoid really closing the pipe until it is
           redirected. We also have to do the duplication before closing any
           ends of the pipes to make sure that we do not end up with the same
           file descriptor number as before, since the POSIX file descriptor
           allocation scheme (always use the lowest available number) mandates
           that this will happen. */
        int rpipe = dup(rpipes[1]);
        if (rpipe == -1)
        {
            cpl_msg_debug(cpl_func, "Failed to duplicate the output pipe for"
                          " the Python process: %s", strerror(errno));
            exit(EXIT_FAILURE);
        }
        int wpipe = dup(wpipes[0]);
        if (wpipe == -1)
        {
            cpl_msg_debug(cpl_func, "Failed to duplicate the input pipe for"
                          " the Python process: %s", strerror(errno));
            exit(EXIT_FAILURE);
        }

        /* The following assertions are here to make sure that PYTHON_INPUT_FD
           and PYTHON_OUTPUT_FD were chosen in such a way that they will not
           overlap with the values of rpipe or wpipe when these get duplicated.
           For example, with PYTHON_INPUT_FD == 3 and PYTHON_OUTPUT_FD == 4 this
           will be correct, since the POSIX file descriptor allocation scheme
           will have rpipe == 7 and wpipe == 8. */
        assert(rpipe != PYTHON_INPUT_FD);
        assert(rpipe != PYTHON_OUTPUT_FD);
        assert(wpipe != PYTHON_INPUT_FD);
        assert(wpipe != PYTHON_OUTPUT_FD);

        /* Now we can close all the ends of the pipes we do not need. */
        close_file(rpipes[0]);
        close_file(rpipes[1]);
        close_file(wpipes[0]);
        close_file(wpipes[1]);

        /* Redirect the saved (duplicated) ends of the pipes to the required
           file descriptor numbers. */
        if (dup2(rpipe, PYTHON_OUTPUT_FD) == -1)
        {
            cpl_msg_debug(cpl_func, "Failed to redirect output pipe for the"
                          " Python process: %s", strerror(errno));
            exit(EXIT_FAILURE);
        }
        if (dup2(wpipe, PYTHON_INPUT_FD) == -1)
        {
            cpl_msg_debug(cpl_func, "Failed to redirect input pipe for the"
                          " Python process: %s", strerror(errno));
            exit(EXIT_FAILURE);
        }
        close_file(rpipe);
        close_file(wpipe);

        /* execlp() will either not return, or we know we have an error so stop
           the process. */
        (void) execlp("python", "python", "-c", python_command, (char *)0);
        cpl_msg_debug(cpl_func, "Failed to execute python command: %s",
                      strerror(errno));
        exit(EXIT_FAILURE);
    }
    else
    {
        /* In parent process after the fork. Close the ends of the pipes we do
           not need, so that we avoid any possible deadlocks or other problems.
           */
        close_file(rpipes[1]);
        close_file(wpipes[0]);
        *python_input = wpipes[1];
        *python_output = rpipes[0];
    }
    return CPL_ERROR_NONE;
}

/*
 * Stops an existing Python interpreter process and joins the process.
 *
 * @param[out] python_pid  The process ID for the Python interpreter.
 * @param[out] python_input  The input file descriptor to close. If this value
 *                           is already set to -1 then it is ignored.
 * @param[out] python_output The output file descriptor to close.
 *
 * @return  @c CPL_ERROR_NONE is returned when the interpreter is stopped and
 *      it returned a zero exit status. Otherwise an appropriate error code is
 *      returned if the interpreter returned a non-zero value.
 */
static cpl_error_code stop_python_interpreter(pid_t python_pid,
                                              int python_input,
                                              int python_output)
{
    int status = -1;
    int result = -1;

    /* Close the input/output pipes to the Python command. Start with
       python_input to let the interpreter know it needs to terminate.
       python_input might have been closed before, so only close it if it
       appears to be valid. */
    if (python_input != -1)
    {
        close_file(python_input);
    }
    close_file(python_output);

    time_t start_time = time(NULL);
    do
    {
        result = waitpid(python_pid, &status, WNOHANG);
        if (result == -1)
        {
            if (errno != EINTR)
            {
                /* waitpid() was not interrupted, thus there was an error when
                   joining the child process. */
                return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                             "Failed to join with the Python"
                                             " process: %s", strerror(errno));
            }
        }
        if (result == 0)
        {
            /* The child process did not complete yet so sleep a little.
               If we have been waiting too long, then also send a termination
               signal.
              */
            time_t now = time(NULL);
            if (now - start_time > 10)  /* Waiting more than 10 seconds? */
            {
                (void) kill(python_pid, SIGTERM);
            }
            struct timespec ts = {0, 10000000};  /* 10 millisecond wait. */
            nanosleep(&ts, NULL);
        }
    } while (result != python_pid);
    if (WIFEXITED(status) && WEXITSTATUS(status) == EXIT_SUCCESS)
    {
        return CPL_ERROR_NONE;
    }
    else
    {
        /* Python terminated with an error. */
        return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                     "Python process existed with error code"
                                     " %d", WEXITSTATUS(status));
    }
}

/*
 * Writes data to the pipe connected to the Python process.
 *
 * @param python_input  The file descriptor for the input pipe to the Python
 *                      process.
 * @param data          Pointer to the data buffer to write to the pipe.
 * @param data_size     The number of bytes in the data buffer.
 *
 * @return  @c CPL_ERROR_NONE is returned on success and an appropriate error
 *      code otherwise.
 */
static cpl_error_code write_to_pipe(int python_input,
                                    const char * data,
                                    size_t data_size)
{
    assert(data != NULL);

    size_t total_written = 0;
    do
    {
        ssize_t bytes_written = write(python_input,
                                      data + total_written,
                                      data_size - total_written);
        if (bytes_written == -1)
        {
            if (errno == EINTR) continue;  /* Got interrupted, try again. */
            return cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                         "Failed to write commands to the"
                                         " Python interpreter: %s",
                                         strerror(errno));
        }
        total_written += bytes_written;
    }
    while (total_written < data_size);

    return CPL_ERROR_NONE;
}

/*
 * Performs blocking reads until all available data has been read from the
 * given pipe.
 *
 * @param python_output  The file descriptor for the output pipe to the Python
 *                       process from which we read.
 *
 * @return  A null terminated string containing the result read from the pipe.
 *
 * @note The caller must release the returned buffer with @c free() from the
 *       C standard library.
 */
static char * read_from_pipe(int python_output)
{
    char * newptr = NULL;
    ssize_t bytes_read = 0;
    size_t total_read = 0;

    /* NOTE: We do not use cpl_malloc or cpl_realloc because those functions
       terminate immediately if a memory allocation failed. But we do not want
       that, since we will still need to join with the Python interpreter child
       process and cleanup properly. */
    size_t data_size = 1024*4;
    char * data = (char *) malloc(data_size);
    if (data == NULL)
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                              "Buffer allocation failed: %s", strerror(errno));
        return NULL;
    }

    do
    {
        /* If the buffer is 50% full then double its size before reading. */
        if (total_read > data_size / 2)
        {
            data_size *= 2;
            newptr = (char *) realloc(data, data_size);
            if (newptr != NULL)
            {
                data = newptr;
            }
            else
            {
                cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                      "Buffer reallocation failed: %s",
                                      strerror(errno));
                free(data);
                return NULL;
            }
        }

        /* Read from the Python command's output stream and check for errors. */
        bytes_read = read(python_output, data + total_read,
                          data_size - 1 - total_read);
        if (bytes_read < 0)
        {
            if (errno == EINTR) continue;  /* Got interrupted, try again. */
            cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                  "Failed to read output from the Python"
                                  " interpreter: %s", strerror(errno));
            free(data);
            return NULL;
        }
        total_read += bytes_read;
    }
    /* Keep reading until no more bytes can be read out. */
    while (bytes_read != 0);

    /* Reduce the memory requirement to what we actually read out of the pipe
       and also make sure the result is a null terminated character string
       before returning it. */
    newptr = (char *) realloc(data, data_size);
    if (newptr != NULL)
    {
        data = newptr;
    }
    else
    {
        cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                              "Buffer reallocation failed: %s",
                              strerror(errno));
        free(data);
        return NULL;
    }
    data[total_read] = '\0';

    return data;
}

/*
 * This helper function produces a debug message for each line for the given
 * text. Each line is also numbered. This is necessary since the CPL messaging
 * facility truncated lines to 1k characters.
 *
 * @param funcname  The name of the function to print in the debug messages.
 *                  This should typically be the name of the function in which
 *                  print_multiline_debug_message() has been used.
 * @param header  An initial message to print before the line numbered text.
 * @param text    The text which is printed with one debug message per line.
 */
static void print_multiline_debug_message(const char * funcname,
                                          const char * header,
                                          const char * text)
{
    /* Print text over multiple debug message lines to produce line numbers
       and thereby also try avoid the max line length limitation of the logging
       subsystem. */
    cpl_msg_debug(funcname, "%s", header);
    char * strbuf = cpl_strdup(text);
    char * line = strbuf;
    int lineno = 1;
    cpl_msg_indent_more();
    char * pchar = line;
    while (*pchar != '\0')
    {
        if (*pchar == '\n')
        {
            *(pchar++) = '\0';
            cpl_msg_debug(funcname, "%.4d %s", lineno, line);
            line = pchar;
            ++lineno;
        }
        else
        {
            ++pchar;
        }
    }
    if (*line != '\0')
    {
        cpl_msg_debug(funcname, "%.4d %s", lineno, line);
    }
    cpl_msg_indent_less();
    cpl_free(strbuf);
}

/*
 * This function starts a Python interpreter invoking it with the given command
 * with the '-c' option. Input text is then written to an input pipe connected
 * to the Python process. This process will then block while trying to read the
 * results written by the interpreter to the output pipe. When the Python
 * process completes an EOF will be read by this process at which point all the
 * results would have been read out. The interpreter is then waits and joins
 * with it child process and getting its exit status.
 * If all goes well the results read out of the output pipe are returned.
 *
 * @param command  The Python command use in the '-c' command line argument.
 * @param input    The input text to write to the input pipe connected to the
 *                 Python process.
 *
 * @return  The output generated by the python process, or @c NULL if an error
 *          occurred.
 *
 * @note The caller must release the returned buffer with @c free() from the
 *       C standard library.
 */
static char * run_python_command(const char * command,
                                 const char * input)
{
    assert(command != NULL);
    assert(input != NULL);

    cpl_msg_debug(cpl_func, "Command sent to Python interpreter:\n%s", command);
    print_multiline_debug_message(cpl_func, "Input to Python interpreter:",
                                  input);

    pid_t python_pid = -1;
    int python_input = -1, python_output = -1;
    if (start_python_interpreter(command, &python_pid,
                                 &python_input, &python_output)
        != CPL_ERROR_NONE)
    {
        /* Error already set in start_python_interpreter, so just return. */
        return NULL;
    }

    if (write_to_pipe(python_input, input, strlen(input)) != CPL_ERROR_NONE)
    {
        /* Cleanup, but preserve the CPL error state from write_to_pipe. */
        cpl_errorstate state = cpl_errorstate_get();
        stop_python_interpreter(python_pid, python_input, python_output);
        cpl_errorstate_set(state);
        return NULL;
    }

    /* First close the input stream to the Python interpreter to let it know
       that no more input is coming and give it a change to run the script. */
    if (close(python_input) != 0)
    {
        if (errno != EINTR)
        {
            cpl_error_set_message(cpl_func, CPL_ERROR_FILE_IO,
                                  "Failed to close input stream to the Python"
                                  " interpreter: %s", strerror(errno));
        }
        /* Cleanup, but preserve the CPL error state.
           NOTE: since we had a problem closing the pipe this might hang the
           python process. Thus, send a kill signal just in case. */
        cpl_errorstate state = cpl_errorstate_get();
        (void) kill(python_pid, SIGTERM);
        stop_python_interpreter(python_pid, -1, python_output);
        cpl_errorstate_set(state);
        return NULL;
    }

    /* Let stop_python_interpreter know that it does not have to try close the
       input pipe, because it has now already been closed. */
    python_input = -1;

    char * output = read_from_pipe(python_output);
    if (output == NULL)
    {
        /* Cleanup, but preserve the CPL error state from read_from_pipe. */
        cpl_errorstate state = cpl_errorstate_get();
        stop_python_interpreter(python_pid, python_input, python_output);
        cpl_errorstate_set(state);
        return NULL;
    }

    if (stop_python_interpreter(python_pid, python_input, python_output)
            != CPL_ERROR_NONE)
    {
        free(output);
        /* Error already set in stop_python_interpreter, so just return. */
        return NULL;
    }

    if (strlen(output) == 0)
    {
        cpl_msg_debug(cpl_func, "Output from Python interpreter: (none)");
    }
    else
    {
        print_multiline_debug_message(cpl_func,
                                      "Output from Python interpreter:",
                                      output);
    }

    return output;
}

/*
 * A comparison function to compare string keys for the cx_map class.
 */
static cxint string_key_compare(cxcptr a, cxcptr b)
{
    return (strcmp(a, b) < 0) ? TRUE : FALSE;
}

/*
 * This comparison function for two plugin objects will return an integer as
 * follows:
 *   -1 if pa < pb
 *    0 if pa == pb
 *    1 if pa > pb
 */
static int plugin_compare(const cpl_plugin * pa, const cpl_plugin * pb)
{
    int result = 0;
    /* NOTE: we order the comparison from what should typically take the least
       amount of time to the most. This should allow us to navigate to the
       matching node within cx_map more quickly. */
    if (pa->api < pb->api) return -1;
    if (pa->api > pb->api) return 1;
    if (pa->version < pb->version) return -1;
    if (pa->version > pb->version) return 1;
    if (pa->type < pb->type) return -1;
    if (pa->type > pb->type) return 1;
    result = strcmp(pa->name, pb->name);
    if (result < 0) return -1;
    if (result > 0) return 1;
    result = strcmp(pa->author, pb->author);
    if (result < 0) return -1;
    if (result > 0) return 1;
    result = strcmp(pa->email, pb->email);
    if (result < 0) return -1;
    if (result > 0) return 1;
    result = strcmp(pa->synopsis, pb->synopsis);
    if (result < 0) return -1;
    if (result > 0) return 1;
    result = strcmp(pa->description, pb->description);
    if (result < 0) return -1;
    if (result > 0) return 1;
    result = strcmp(pa->copyright, pb->copyright);
    if (result < 0) return -1;
    if (result > 0) return 1;
    return 0;
}

/*
 * A comparison function to compare cpl_plugin like keys for the cx_map class.
 */
static cxint plugin_key_compare(cxcptr a, cxcptr b)
{
    const cpl_plugin * pa = (const cpl_plugin *) a;
    const cpl_plugin * pb = (const cpl_plugin *) b;
    return (plugin_compare(pa, pb) < 0) ? TRUE : FALSE;
}

/**
 * @brief Loads relevant information for Python modules that contain plugins.
 *
 * This function will go through a list of module paths that point to Python
 * code files and attempt to load these. If the module can be loaded
 * successfully and it contains any classes that derive from a base class called
 * @c CplPlugin (including @c CplPlugin itself) then these are recorded, such
 * that er_python_select_module() can then be used to select a particular Python
 * module to run. Only classes directly within the Python module's namespace are
 * considered, i.e. secondary imported modules will not have their classes added
 * unless they are imported into the namespace with a
 * \"<tt>from ... import ...</tt>\" Python statement.
 *
 * @param[in,out] paths  The list of module paths that need to be checked if
 *                       they can be loaded. Any paths that do not contain
 *                       loadable plugins will be removed from the list by this
 *                       function.
 *
 * @return @c CPL_ERROR_NONE on success and an appropriate error code if a
 *      severe error occurred.
 *
 * @note This function will still succeed if a Python module could not be
 *      loaded. All such modules are simply ignored. However, the user may get
 *      error messages printed by the Python process on stderr.
 */
cpl_error_code er_python_load_modules(er_stringarray_t * paths)
{
    cpl_error_ensure(paths != NULL, CPL_ERROR_NULL_INPUT,
                     return CPL_ERROR_NULL_INPUT, "paths is NULL.");

    er_python_cleanup();
    module_list = cx_map_new(string_key_compare, cpl_free,
                             (cx_free_func)er_json_pluginlist_delete);
    /* NOTE: we do not set a destructor for the keys or values of the module_lut
       mapping since these are pointing to objects inside module_list. */
    module_lut = cx_map_new(plugin_key_compare, NULL, NULL);

    /* Prepare a JSON compatible string containing a list of module names. */
    char * paths_string = NULL;
    int i = 0;
    for (i = 0; i < er_stringarray_size(paths); ++i)
    {
        char * path = er_json_escape_string(er_stringarray_get(paths, i));
        if (i == 0)
        {
            paths_string = cpl_sprintf("\"%s\"", path);
        }
        else
        {
            char * tmpstr = cpl_sprintf("%s,\n    \"%s\"", paths_string, path);
            cpl_free(paths_string);
            paths_string = tmpstr;
        }
        cpl_free(path);
    }
    if (paths_string == NULL)
    {
        /* No paths to load. */
        return CPL_ERROR_NONE;
    }

    /* The following script attempts to find Python based plugins for EsoRex
       in the module paths identified in _data['paths']. */
    char * script = er_json_escape_string(
        "import inspect\n"
        "valid_modules = {}\n"
        "for path in _data['paths']:\n"
        /* First try import the python module file. */
        "  module = os.path.splitext(os.path.basename(path))[0]\n"
        "  sys.path.append(os.path.dirname(path))\n"
        "  try:\n"
        "    __import__(module)\n"
        /*   Look for any classes that inherit from a class called CplPlugin.
             If such a class is found then add the name under the corresponding
             module path key in the valid modules dictionary. */
        "    for name, cls in inspect.getmembers(sys.modules[module],"
                                               " predicate=inspect.isclass):\n"
        "      if len([base for base in cls.__mro__"
                     " if base.__name__ == 'CplPlugin']) > 0:\n"
        "        if path not in valid_modules:\n"
        "          valid_modules[path] = set()\n"
        "        valid_modules[path].add(cls)\n"
        "  except:\n"
        "    import traceback\n"
        "    traceback.print_exc()\n"
        "  del sys.path[-1]\n"
        "results = {}\n"
        "for path, clslist in valid_modules.items():\n"
        "  results[path] = []\n"
        "  for cls in clslist:\n"
        "    plugin = {'class': cls.__name__}\n"
        /*   Create an instance of the class and find all the relevant
             information to fill into the plugin dictionary. */
        "    try:\n"
        "      inst = cls()\n"
        "    except:\n"
        "      continue\n"
        "    try:\n"
        "      name = inst.name\n"
        "    except:\n"
        "      name = cls.__name__\n"
        "    try:\n"
        "      version = inst.version\n"
        "    except:\n"
        "      version = None\n"
        "    try:\n"
        "      synopsis = inst.synopsis\n"
        "    except:\n"
        "      try:\n"
        "        synopsis = inspect.getdoc(inst).splitlines()[0]\n"
        "      except:\n"
        "        synopsis = None\n"
        "    try:\n"
        "      description = inst.description\n"
        "    except:\n"
        "      try:\n"
        "        description = inspect.getdoc(inst)\n"
        "      except:\n"
        "        description = None\n"
        "    try:\n"
        "      author = inst.author\n"
        "    except:\n"
        "      author = None\n"
        "    try:\n"
        "      email = inst.email\n"
        "    except:\n"
        "      email = None\n"
        "    try:\n"
        "      copyright = inst.copyright\n"
        "    except:\n"
        "      copyright = None\n"
        "    plugin['name'] = name\n"
        "    plugin['version'] = version\n"
        "    plugin['synopsis'] = synopsis\n"
        "    plugin['description'] = description\n"
        "    plugin['author'] = author\n"
        "    plugin['email'] = email\n"
        "    plugin['copyright'] = copyright\n"
        "    try:\n"
        "      plugin['parameters'] = inst.parameters\n"
        "    except:\n"
        "      pass\n"
        "    try:\n"
        "      plugin['recipeconfig'] = inst.recipeconfig\n"
        "    except:\n"
        "      pass\n"
        "    results[path].append(plugin)\n"
        /* Generate and write out the results as JSON. */
        "json.dump(results,"
                 " os.fdopen("CPL_STRINGIFY(PYTHON_OUTPUT_FD)",'w'),"
                 " indent=2)\n"
    );

    /* Create the input for the above script.
       NOTE: cpl_sprintf will create a valid buffer or terminate the process. */
    char * input = cpl_sprintf(
            "{\n"
            "  \"script\": \"%s\",\n"
            "  \"paths\": [\n"
            "    %s\n"
            "  ]\n"
            "}",
            script, paths_string
        );
    cpl_free(paths_string);

    char * output = run_python_command(python_start_command, input);
    cpl_free(script);
    cpl_free(input);
    if (output == NULL) return cpl_error_get_code();

    /* Parse the output JSON into a parse tree of JSON nodes. These can then
       be converted into appropriate CPL structures. */
    er_json_node * parsetree = er_json_parse(output);
    if (parsetree == NULL)
    {
        free(output);
        return cpl_error_get_code();
    }

    /* Make sure the top level node in the parse tree is an object. */
    if (er_json_node_type(parsetree) != JSON_OBJECT)
    {
        int line, col;
        er_json_find_line_column(output, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an object at line %d column %d.",
                              line, col);
        er_json_node_delete(parsetree);
        free(output);
        return cpl_error_get_code();
    }

    er_stringarray_t * valid_paths = er_stringarray_new();

    er_json_object_iterator iter  = NULL;
    for (iter = er_json_node_object_begin(parsetree);
         iter != er_json_node_object_end(parsetree);
         iter = er_json_node_object_next(parsetree, iter))
    {
        const char * module = er_json_node_object_get_key(parsetree, iter);
        const er_json_node * item = er_json_node_object_get_value(parsetree,
                                                                  iter);
        cpl_pluginlist * pluginlist = er_json_to_pluginlist(item, output);
        /* NOTE: If no plugin list could be created then skip the module,
           but continue processing the other possible modules. */
        if (pluginlist != NULL)
        {
            er_stringarray_append(valid_paths, module);

            /* Need to add a copy to the map object because it will be later
               deleted by it in er_python_cleanup. */
            char * module_name = cpl_strdup(module);

            cx_map_insert(module_list, module_name, pluginlist);

            /* Add all the plugin entries to the reverse lookup table. */
            cpl_plugin * plugin = NULL;
            for (plugin = cpl_pluginlist_get_first(pluginlist);
                 plugin != NULL;
                 plugin = cpl_pluginlist_get_next(pluginlist))
            {
                cx_map_insert(module_lut, plugin, module_name);
            }
        }
    }

    er_json_node_delete(parsetree);
    free(output);

    /* For every path in the input paths list, we check if it was found in
       valid_paths. Any path not found in valid_paths is removed. */
    for (i = 0; i < er_stringarray_size(paths); ++i)
    {
        const char * path = er_stringarray_get(paths, i);
        int found = 0;
        int j = 0;
        for (j = 0; j < er_stringarray_size(valid_paths); ++j)
        {
            if (strcmp(path, er_stringarray_get(valid_paths, j)) == 0)
            {
                found = 1;
                break;
            }
        }
        if (! found)
        {
            er_stringarray_remove(paths, i);
        }
    }
    er_stringarray_delete(valid_paths);

    return cpl_error_get_code();
}

/**
 * @brief Selects a specific Python module to use for execution.
 *
 * This function is used to select the Python module that will be used when
 * the er_python_get_plugin_list() function is invoked and attempts to load
 * all available plugin classes.
 *
 * @param module_name  The full path name of the python module to use.
 *
 * @return @c CPL_ERROR_NONE on success or an appropriate error code if a
 *      severe error occurs.
 */
cpl_error_code er_python_select_module(const char * module_name)
{
    if (module_name == NULL)
    {
        current_module = NULL;
        return CPL_ERROR_NONE;
    }

    if (module_list == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "Python modules need to be loaded first"
                                     " with er_python_load_modules.");
    }
    cx_map_iterator module = cx_map_find(module_list, module_name);
    if (module != cx_map_end(module_list))
    {
        current_module = module;
        return CPL_ERROR_NONE;
    }
    else
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "Failed to find preloaded Python module:"
                                     " %s", module_name);
    }
}

#ifdef HAVE_BROKEN_CPL_PARAMETER_DUPLICATE

/*
 * The following is a workaround for a bug in cpl_parameter_duplicate, which
 * incorrectly copies ranges.
 * This function will make a duplicate copy of a CPL parameter and return it.
 */
static cpl_parameter * fixed_parameter_duplicate(const cpl_parameter * par)
{
    if (cpl_parameter_get_class(par) == CPL_PARAMETER_CLASS_RANGE)
    {
        /* Create a new parameter range object if we were dealing with a range
           to begin with. */
        cpl_parameter * newpar;
        switch (cpl_parameter_get_type(par))
        {
            case CPL_TYPE_INT:
                newpar = cpl_parameter_new_range(
                                        cpl_parameter_get_name(par),
                                        CPL_TYPE_INT,
                                        cpl_parameter_get_help(par),
                                        cpl_parameter_get_context(par),
                                        cpl_parameter_get_default_int(par),
                                        cpl_parameter_get_range_min_int(par),
                                        cpl_parameter_get_range_max_int(par)
                                    );
                cpl_parameter_set_int(newpar, cpl_parameter_get_int(par));
                break;

            case CPL_TYPE_DOUBLE:
                newpar = cpl_parameter_new_range(
                                        cpl_parameter_get_name(par),
                                        CPL_TYPE_DOUBLE,
                                        cpl_parameter_get_help(par),
                                        cpl_parameter_get_context(par),
                                        cpl_parameter_get_default_double(par),
                                        cpl_parameter_get_range_min_double(par),
                                        cpl_parameter_get_range_max_double(par)
                                    );
                cpl_parameter_set_double(newpar, cpl_parameter_get_double(par));
                break;

            default:
                return cpl_parameter_duplicate(par);
        }

        /* We need to copy over additional information that is not settable by
           the cpl_parameter_new_range function. */
        cpl_parameter_set_default_flag(
                        newpar, cpl_parameter_get_default_flag(par));
        cpl_parameter_set_id(newpar, cpl_parameter_get_id(par));
        if (cpl_parameter_get_tag(par) != NULL)
        {
            cpl_parameter_set_tag(newpar, cpl_parameter_get_tag(par));
        }

        cpl_parameter_mode mode[3] = {
                CPL_PARAMETER_MODE_CLI,
                CPL_PARAMETER_MODE_ENV,
                CPL_PARAMETER_MODE_CFG
            };
        int n = 0;
        for (n = 0; n < 3; ++n)
        {
            cpl_parameter_set_alias(newpar, mode[n],
                                    cpl_parameter_get_alias(par, mode[n]));
            if (cpl_parameter_is_enabled(par, mode[n]))
            {
                cpl_parameter_enable(newpar, mode[n]);
            }
            else
            {
                cpl_parameter_disable(newpar, mode[n]);
            }
        }

        return newpar;
    }
    else
    {
        /* If dealing with any other parameter class then just use the normal
           duplicate method which should work fine. */
        return cpl_parameter_duplicate(par);
    }
}

#endif /* HAVE_BROKEN_CPL_PARAMETER_DUPLICATE */

/*
 * Releases a NULL terminated array of strings and all the strings that it
 * contains.
 */
static void free_string_array(char ** array)
{
    if (array == NULL) return;
    char ** str = NULL;
    for (str = array; *str != NULL; ++str)
    {
        cpl_free(*str);
    }
    cpl_free(array);
}

/*
 * This is the constructor method that is set by er_python_get_plugin_list()
 * for the new plugin objects that it creates.
 * It will attempt to find the preloaded data corresponding to the given plugin.
 * The preloaded data should have been filled in the global lookup tables by
 * er_python_load_modules() before hand.
 * Once the corresponding information is found it is copied over to the given
 * plugin to initialise it. Specifically things like the recipe parameters and
 * recipe configuration information is copied over.
 */
static int plugin_initialize(cpl_plugin * plugin)
{
    if (plugin == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
                                     "NULL received for the plugin.");
    }

    unsigned int plugin_api = cpl_plugin_get_api(plugin);
    if (plugin_api != CPL_PLUGIN_API)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                                     "Plugin API version %u is not supported.",
                                     plugin_api);
    }

    /* Find the corresponding plugin object we should have already loaded
       beforehand by er_python_load_modules. We will need this originally loaded
       object to be able to correctly fill in the parameters and configuration
       for the plugin structure received by this function. */
    cx_map_iterator plgiter = cx_map_find(module_lut, plugin);
    if (plgiter == cx_map_end(module_lut))
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "Failed to find preloaded Python plugin"
                                     " object for '%s'.",
                                     cpl_plugin_get_name(plugin));
    }
    cpl_plugin * preloaded_plugin = (cpl_plugin *) cx_map_get_key(module_lut,
                                                                  plgiter);
    assert(cpl_plugin_get_type(plugin)
           == cpl_plugin_get_type(preloaded_plugin));

    /* Allocate the necessary data structures for the recipe plugin, based on
       the type of plugin we are dealing with. Also find the corresponding
       data strutures in the preloaded plugin. */
    cpl_parameterlist * params = NULL;
    cpl_recipeconfig * config = NULL;
    cpl_parameterlist * preloaded_params = NULL;
    cpl_recipeconfig * preloaded_config = NULL;
    unsigned long plugin_type = cpl_plugin_get_type(plugin);
    switch (plugin_type)
    {
        case CPL_PLUGIN_TYPE_RECIPE:
            {
                cpl_recipe * recipe = (cpl_recipe *) plugin;
                recipe->parameters = cpl_parameterlist_new();
                params = recipe->parameters;

                recipe = (cpl_recipe *) preloaded_plugin;
                preloaded_params = recipe->parameters;
            }
            break;

        case CPL_PLUGIN_TYPE_RECIPE_V2:
            {
                cpl_recipe2 * recipe = (cpl_recipe2 *) plugin;
                recipe->base.parameters = cpl_parameterlist_new();
                params = recipe->base.parameters;
                recipe->config = cpl_recipeconfig_new();
                config = recipe->config;

                recipe = (cpl_recipe2 *) preloaded_plugin;
                preloaded_params = recipe->base.parameters;
                preloaded_config = recipe->config;
            }
            break;

        default:
            return cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                                         "Plugin type %lu is not supported.",
                                         plugin_type);
    }

    /* Copy default parameters into the received plugin object. NOTE: we have
       to skip the __class__ parameter from the preloaded parameter list since
       this was added internally when parsing the JSON output from the Python
       interpreter. */
    if (preloaded_params != NULL)
    {
        const cpl_parameter * last =
                            cpl_parameterlist_get_last_const(preloaded_params);
        assert(strcmp(cpl_parameter_get_name(last), "__class__") == 0);
        const cpl_parameter * par = NULL;
        for (par = cpl_parameterlist_get_first_const(preloaded_params);
             par != last;
             par = cpl_parameterlist_get_next_const(preloaded_params))
        {
#ifdef HAVE_BROKEN_CPL_PARAMETER_DUPLICATE
            cpl_parameter * newpar = fixed_parameter_duplicate(par);
#else
            cpl_parameter * newpar = cpl_parameter_duplicate(par);
#endif
            if (cpl_parameterlist_append(params, newpar) != CPL_ERROR_NONE)
            {
                cpl_parameter_delete(newpar);
                return cpl_error_get_code();
            }
        }
    }

    /* Copy the recipe configuration structure if we are dealing with a version
       2 type recipe structure. */
    if (preloaded_config != NULL)
    {
        char ** inputs = NULL;
        char ** outputs = NULL;
        cpl_size min, max;

        char ** tags = cpl_recipeconfig_get_tags(preloaded_config);
        if (tags == NULL) return cpl_error_get_code();

        char ** tag = NULL;
        for (tag = tags; *tag != NULL; ++tag)
        {
            inputs = cpl_recipeconfig_get_inputs(preloaded_config, *tag);
            if (inputs == NULL) break;
            outputs = cpl_recipeconfig_get_outputs(preloaded_config, *tag);
            if (outputs == NULL) break;

            min = cpl_recipeconfig_get_min_count(preloaded_config, *tag, *tag);
            max = cpl_recipeconfig_get_max_count(preloaded_config, *tag, *tag);

            if (cpl_recipeconfig_set_tag(config, *tag, min, max) != 0) break;

            char ** input = NULL;
            for (input = inputs; *input != NULL; ++input)
            {
                min = cpl_recipeconfig_get_min_count(preloaded_config,
                                                     *tag, *input);
                max = cpl_recipeconfig_get_max_count(preloaded_config,
                                                     *tag, *input);

                if (cpl_recipeconfig_set_input(config, *tag, *input, min, max)
                    != 0)
                {
                    goto cleanup;
                }
            }

            char ** output = NULL;
            for (output = outputs; *output != NULL; ++output)
            {
                if (cpl_recipeconfig_set_output(config, *tag, *output) != 0)
                {
                    goto cleanup;
                }
            }

            free_string_array(outputs);
            free_string_array(inputs);
            inputs = outputs = NULL;
        }

    cleanup:
        free_string_array(outputs);
        free_string_array(inputs);
        free_string_array(tags);

        /* Check if we had an error during the copy. */
        if (cpl_error_get_code() != CPL_ERROR_NONE) return cpl_error_get_code();
    }

    return CPL_ERROR_NONE;
}

/*
 * This function will update the frameset data in 'plugin' from the frameset
 * found in 'result'. Values for old entries in 'plugin' are updated and any
 * new frames in 'result' are copied over (a deep copy is performed).
 *
 * @return @c CPL_ERROR_NONE on success and an appropriate error code otherwise.
 */
static cpl_error_code update_plugin(cpl_plugin * plugin, cpl_plugin * result)
{
    assert(plugin != NULL);
    assert(result != NULL);

    /* Perform some sanity checks to make sure that the two plugin structures,
       are compatible. i.e. they must have the same API number, name, version
       and type code. */
    if (cpl_plugin_get_api(plugin) != cpl_plugin_get_api(result))
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                     "Plugin API number %u for the returned"
                                     " plugin does not correspond to the"
                                     " original value %u.",
                                     cpl_plugin_get_api(result),
                                     cpl_plugin_get_api(plugin));
    }

    if (cpl_plugin_get_version(plugin) != cpl_plugin_get_version(result))
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                     "Plugin version %lu for the returned"
                                     " plugin does not correspond to the"
                                     " original value %lu.",
                                     cpl_plugin_get_version(result),
                                     cpl_plugin_get_version(plugin));
    }

    unsigned long plugin_type = cpl_plugin_get_type(plugin);
    if (plugin_type != cpl_plugin_get_type(result))
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                     "Plugin type code %lu for the returned"
                                     " plugin does not correspond to the"
                                     " original value %lu.",
                                     cpl_plugin_get_type(result), plugin_type);
    }

    if (strcmp(cpl_plugin_get_name(plugin), cpl_plugin_get_name(result)) != 0)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_TYPE_MISMATCH,
                                     "Plugin name '%s' for the returned"
                                     " plugin does not correspond to the"
                                     " original value '%s'.",
                                     cpl_plugin_get_name(result),
                                     cpl_plugin_get_name(plugin));
    }

    /* Prepare the pointers to the frameset objects in each of the plugins.
       This will then make it easier to update the data structures. */
    cpl_frameset * plugin_frames = NULL;
    cpl_frameset * result_frames = NULL;
    switch (plugin_type)
    {
        case CPL_PLUGIN_TYPE_RECIPE:
            {
                cpl_recipe * recipe = (cpl_recipe *) plugin;
                plugin_frames = recipe->frames;
                recipe = (cpl_recipe *) result;
                result_frames = recipe->frames;
            }
            break;

        case CPL_PLUGIN_TYPE_RECIPE_V2:
            {
                cpl_recipe2 * recipe = (cpl_recipe2 *) plugin;
                plugin_frames = recipe->base.frames;
                recipe = (cpl_recipe2 *) result;
                result_frames = recipe->base.frames;
            }
            break;

        default:
            return cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                                         "Plugin type %lu is not supported.",
                                         plugin_type);
    }

    if (result_frames == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Received and empty frameset from the"
                                     " recipe plugin as output.");
    }
    if (plugin_frames == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                     "Plugin does not have its frameset"
                                     " initialised for updating.");
    }

    /* Copy the first N frame values of from the resultant plugin to the
       existing plugin. N is the number of frames passed to the plugin recipe
       as input. This is necessary since the recipe plugin is expected to update
       the frame information. */
    cpl_frameset_iterator * ip = cpl_frameset_iterator_new(plugin_frames);
    cpl_frameset_iterator * ir = cpl_frameset_iterator_new(result_frames);
    cpl_frame * fp = cpl_frameset_iterator_get(ip);
    cpl_frame * fr = cpl_frameset_iterator_get(ir);
    while ((fp = cpl_frameset_iterator_get(ip)) != NULL)
    {
        fr = cpl_frameset_iterator_get(ir);
        if (fr == NULL)
        {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                                  "The frameset returned from the plugin recipe"
                                  " is shorter than what was given to it.");
            break;
        }

        /* We only copy the filename and tag if these were really changed.
           This avoids mallocs and frees. */
        if (strcmp(cpl_frame_get_filename(fp), cpl_frame_get_filename(fr)) != 0)
        {
            cpl_frame_set_filename(fp, cpl_frame_get_filename(fr));
        }
        if (strcmp(cpl_frame_get_tag(fp), cpl_frame_get_tag(fr)) != 0)
        {
            cpl_frame_set_tag(fp, cpl_frame_get_tag(fr));
        }

        cpl_frame_set_type(fp, cpl_frame_get_type(fr));
        cpl_frame_set_group(fp, cpl_frame_get_group(fr));
        cpl_frame_set_level(fp, cpl_frame_get_level(fr));

        cpl_frameset_iterator_advance(ip, 1);
        cpl_frameset_iterator_advance(ir, 1);
    }
    cpl_frameset_iterator_delete(ip);

    /* Insert copies of all the new frames to the output plugin structure. */
    while ((fr = cpl_frameset_iterator_get(ir)) != NULL)
    {
        cpl_frame * newframe = cpl_frame_duplicate(fr);
        if (cpl_frameset_insert(plugin_frames, newframe)
            != CPL_ERROR_NONE)
        {
            cpl_frame_delete(newframe);
            break;
        }
        cpl_frameset_iterator_advance(ir, 1);
    }
    cpl_frameset_iterator_delete(ir);

    return cpl_error_get_code();
}

/*
 * This is the recipe execution method that is set by
 * er_python_get_plugin_list() for the new plugin objects that it creates.
 * This will execute a Python based recipe plugin in a Python interpreter
 * process.
 *
 * Requirements for the Python plugin class:
 *
 * 1) It must derive from a class called CplPlugin.
 *
 * 2) __init__ must take no parameters or all parameters must have defaults.
 *
 * 3) Must have one method called execute with the following format:
 *
 *      execute(self, plugin) -> int
 *
 * 4) It should have the following class attributes defined (alternatively these
 *    can also be set by the __init__ function):
 *
 *      name (string) - else taken from the Python class's name.
 *      version (int)
 *      synopsis (string) - else taken from first line of doc string.
 *      description (string) - else taken from doc string.
 *      author (string)
 *      email (string)
 *      copyright (string)
 *      parameters (list) - a list of dictionaries with the following format:
 *          {
 *              'class': (string),
 *              'name': (string),
 *              'id': (integer),
 *              'description': (string),
 *              'context': (string),
 *              'value': (string) | (integer) | (double) | (boolean),
 *              'default': (string) | (integer) | (double) | (boolean),
 *              'min': (integer) | (double),   # Only defined for ranges.
 *              'max': (integer) | (double),   # Only defined for ranges.
 *              'choices': [                   # Only defined for enumerations.
 *                      (string) | (integer) | (double),
 *                      (string) | (integer) | (double),
 *                      ...
 *                  ],
 *              'present': True | False,
 *              'tag': (string),
 *              'cli_enabled': True | False,
 *              'cli_alias': (string) | None,
 *              'env_enabled': True | False,
 *              'env_alias': (string) | None,
 *              'cfg_enabled': True | False,
 *              'cfg_alias': (string) | None
 *          }
 *      recipeconfig (list) - this is optional, only needed for v2 type recipes.
 *                            If given then this is a list of dictionaries with
 *                            the following format:
 *          {
 *              'tag': (string),
 *              'min': (integer),
 *              'max': (integer),
 *              'inputs': [
 *                      {'tag': (string), 'min': (integer), 'max': (integer)},
 *                      {'tag': (string), 'min': (integer), 'max': (integer)},
 *                      ...
 *                  ],
 *              'outputs': [(string), (string), ...]
 *          }
 */
static int plugin_execute(cpl_plugin * plugin)
{
    if (plugin == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
                                     "NULL received for the plugin.");
    }

    unsigned int plugin_api = cpl_plugin_get_api(plugin);
    if (plugin_api != CPL_PLUGIN_API)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                                     "Plugin API version %u is not supported.",
                                     plugin_api);
    }

    /* Find the corresponding plugin object we should have already loaded
       beforehand by er_python_load_modules. We will need this to be able to
       correctly identify the Python module and class name to use when invoking
       the Python recipe code. */
    cx_map_iterator plgiter = cx_map_find(module_lut, plugin);
    if (plgiter == cx_map_end(module_lut))
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "Failed to find preloaded Python plugin"
                                     " object for '%s'.",
                                     cpl_plugin_get_name(plugin));
    }
    cpl_plugin * preloaded_plugin = (cpl_plugin *) cx_map_get_key(module_lut,
                                                                  plgiter);
    const char * module_path = (const char *) cx_map_get_value(module_lut,
                                                               plgiter);
    assert(cpl_plugin_get_type(plugin)
           == cpl_plugin_get_type(preloaded_plugin));

    /* Identify the needed parts of the recipe data that for the preloaded
       plugin object. */
    cpl_parameterlist * preloaded_params = NULL;
    unsigned long plugin_type = cpl_plugin_get_type(plugin);
    switch (plugin_type)
    {
        case CPL_PLUGIN_TYPE_RECIPE:
            {
                cpl_recipe * recipe = (cpl_recipe *) preloaded_plugin;
                preloaded_params = recipe->parameters;
            }
            break;

        case CPL_PLUGIN_TYPE_RECIPE_V2:
            {
                cpl_recipe2 * recipe = (cpl_recipe2 *) preloaded_plugin;
                preloaded_params = recipe->base.parameters;
            }
            break;

        default:
            return cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                                         "Plugin type %lu is not supported.",
                                         plugin_type);
    }

    /* Identify the Python class that we will have to instantiate. */
    if (preloaded_params == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "Missing parameter list in preloaded"
                                     " Python plugin object for '%s'.",
                                     cpl_plugin_get_name(plugin));
    }
    const cpl_parameter * param =
        cpl_parameterlist_get_last_const(preloaded_params);
    if (param == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "Parameter list in preloaded Python plugin"
                                     " object for '%s' is empty.",
                                     cpl_plugin_get_name(plugin));
    }
    assert(strcmp(cpl_parameter_get_name(param), "__class__") == 0);
    const char * class_name = cpl_parameter_get_string(param);
    if (class_name == NULL) return cpl_error_get_code();

    /* Need to convert the plugin object to JSON text to send to the Python
       interpreter. */
    char * plugin_data = er_plugin_to_json(plugin);

    /* Extract the module name to be used in the Python import statement.
       We also extract the directory to the module while we are at it.
       NOTE: basename takes a "char *" string, not a "const char *" string.
       Thus, we create a modifiable string, just in case. */
    char * modpath = cpl_strdup(module_path);
    char * module_name = cpl_strdup(basename(modpath));
    cpl_free(modpath);
    modpath = cpl_strdup(module_path);
    char * module_dir = cpl_strdup(dirname(modpath));
    cpl_free(modpath);
    char * pos = strchr(module_name, '.');
    if (pos != NULL) *pos = '\0';

    /* The following script attempts to load the Python module, instantiate the
       appropriate Python class and invoke it's execute method. The updated data
       is finally converted to JSON and written to the pipe connected to the
       EsoRex process. */
    char * script = cpl_sprintf(
        "sys.path.append(_data['modulepath'])\n"
        "import %s\n"
        "plugin = %s.%s()\n"
        "result = plugin.execute(_data['plugin'])\n"
        "if hasattr(plugin, 'error_message'):\n"
        "  error_message = getattr(plugin, 'error_message')\n"
        "else:"
        "  error_message = None\n"
        "json.dump({'result': result, 'error': error_message,"
                  " 'plugin': _data['plugin']},"
                 " os.fdopen("CPL_STRINGIFY(PYTHON_OUTPUT_FD)",'w'),"
                 " indent=2)\n",
        module_name, module_name, class_name
    );
    cpl_free(module_name);
    char * excaped_script = er_json_escape_string(script);
    cpl_free(script);

    /* Create the input for the above script.
       NOTE: cpl_sprintf will create a valid buffer or terminate the process. */
    char * input = cpl_sprintf(
            "{\n"
            "  \"script\": \"%s\",\n"
            "  \"modulepath\": \"%s\",\n"
            "  \"plugin\": %s\n"
            "}",
            excaped_script, module_dir, plugin_data);
    cpl_free(excaped_script);
    cpl_free(module_dir);
    cpl_free(plugin_data);

    char * output = run_python_command(python_start_command, input);
    cpl_free(input);
    if (output == NULL) return cpl_error_get_code();

    /* Parse the output JSON and convert to appropriate CPL objects. */
    er_json_node * parsetree = er_json_parse(output);
    if (parsetree == NULL)
    {
        free(output);
        return cpl_error_get_code();
    }

    /* Make sure the top level node in the parse tree is an object. */
    if (er_json_node_type(parsetree) != JSON_OBJECT)
    {
        int line, col;
        er_json_find_line_column(output, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected an object at line %d column %d.",
                              line, col);
        er_json_node_delete(parsetree);
        free(output);
        return cpl_error_get_code();
    }

    /* Get the result field from the parsed output. */
    const er_json_node * node = er_json_node_object_get(parsetree, "result");
    if (node == NULL)
    {
        int line, col;
        er_json_find_line_column(output, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Missing a 'result' key in output object at line"
                              " %d column %d.", line, col);
        er_json_node_delete(parsetree);
        free(output);
        return cpl_error_get_code();
    }
    if (er_json_node_type(node) != JSON_NUMBER)
    {
        int line, col;
        er_json_find_line_column(output, er_json_node_location(node),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Expected a number for key 'result' in output"
                              " object at line %d column %d.", line, col);
        er_json_node_delete(parsetree);
        free(output);
        return cpl_error_get_code();
    }
    int result = (int) er_json_node_get_number(node);
    if (result != 0)
    {
        /* An error was marked by returning non-zero from the Python plugin's
           execute method call. Now attempt to get the error message. If there
           was not error message as a atring object then just produce a defualt
           error message instead. */
        const char * error_msg = NULL;
        node = er_json_node_object_get(parsetree, "error");
        if (node != NULL && er_json_node_type(node) == JSON_STRING)
        {
            error_msg = er_json_node_get_string(node);
        }
        if (error_msg != NULL)
        {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT, "%s",
                                  error_msg);
        }
        else
        {
            cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
                                  "Plugin execution failed with error code %d",
                                  result);
        }
        er_json_node_delete(parsetree);
        free(output);
        return result;
    }

    /* Get the plugin structure from the parsed output and update the plugin
       object that was passed to this function, by copying all the new frames
       back to the plugin object. */
    node = er_json_node_object_get(parsetree, "plugin");
    if (node == NULL)
    {
        int line, col;
        er_json_find_line_column(output, er_json_node_location(parsetree),
                                 &line, &col);
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "Missing a 'plugin' key in output object at line"
                              " %d column %d.", line, col);
        er_json_node_delete(parsetree);
        free(output);
        return cpl_error_get_code();
    }
    cpl_plugin * plugin_result = er_json_to_plugin(node, output);
    if (plugin_result == NULL)
    {
        er_json_node_delete(parsetree);
        free(output);
        return cpl_error_get_code();
    }
    cpl_error_code errorcode = update_plugin(plugin, plugin_result);

    /* Cleanup the returned plugin structure, parsetree and output test. */
    plugin_result->deinitialize(plugin_result);
    cpl_plugin_delete(plugin_result);
    er_json_node_delete(parsetree);
    free(output);

    return errorcode;
}

/*
 * This is the destructor method that is set by er_python_get_plugin_list()
 * for the new plugin objects that it creates. It will free memory by deleting
 * any objects created and added to the plugin by plugin_initialize().
 */
static int plugin_deinitialize(cpl_plugin * plugin)
{
    if (plugin == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
                                     "NULL received for the plugin.");
    }

    unsigned int plugin_api = cpl_plugin_get_api(plugin);
    if (plugin_api != CPL_PLUGIN_API)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                                     "Plugin API version %u is not supported.",
                                     plugin_api);
    }

    unsigned long plugin_type = cpl_plugin_get_type(plugin);
    switch (plugin_type)
    {
        case CPL_PLUGIN_TYPE_RECIPE:
            {
                cpl_recipe * recipe = (cpl_recipe *) plugin;
                cpl_parameterlist_delete(recipe->parameters);
            }
            break;

        case CPL_PLUGIN_TYPE_RECIPE_V2:
            {
                cpl_recipe2 * recipe = (cpl_recipe2 *) plugin;
                cpl_recipeconfig_delete(recipe->config);
                cpl_parameterlist_delete(recipe->base.parameters);
            }
            break;

        default:
            return cpl_error_set_message(cpl_func, CPL_ERROR_UNSUPPORTED_MODE,
                                         "Plugin type %lu is not supported.",
                                         plugin_type);
    }

    return CPL_ERROR_NONE;
}

/**
 * @brief Fills a list with all available plugins for a selected Python module.
 *
 * This function is the equivalent of the @c cpl_plugin_get_info() function
 * found in a compiled plugin library. It will fill a list of CPL plugin objects
 * for all plugins found in the currently selected Python module.
 * er_python_select_module() must be called before this function is invoked.
 *
 * @param list  The plugin list that will be filled with CPL plugin objects.
 *
 * @return @c 0 is returned on success and an appropriate CPL error code
 *         otherwise.
 */
int er_python_get_plugin_list(cpl_pluginlist * list)
{
    cpl_error_ensure(list != NULL, CPL_ERROR_NULL_INPUT,
                     return CPL_ERROR_NULL_INPUT,
                     "NULL pointer given for list.");

    if (current_module == NULL)
    {
        return cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                     "Python module has not been selected yet"
                                     " with er_python_select_module.");
    }

    cpl_pluginlist * pluginlist =
        (cpl_pluginlist *) cx_map_get_value(module_list, current_module);
    assert(pluginlist != NULL);

    /* Copy all the plugins we found in er_python_load_modules into the output
       list. But we only copy the base plugin information, not the parameters.
       The parameters need to be initialised in the plugin->initialize function.
       */
    cpl_plugin * rawplugin = NULL;
    for (rawplugin = cpl_pluginlist_get_first(pluginlist);
         rawplugin != NULL;
         rawplugin = cpl_pluginlist_get_next(pluginlist))
    {
        cpl_plugin * plugin = cpl_plugin_new();
        cpl_error_code result = CPL_ERROR_NONE;
        result |= cpl_plugin_copy(plugin, rawplugin);

        /* We need to override the init, exec and deinit functions to ones that
           will appropiately interact with the EsoRex API interface routines.
           */
        result |= cpl_plugin_set_init(plugin, plugin_initialize);
        result |= cpl_plugin_set_exec(plugin, plugin_execute);
        result |= cpl_plugin_set_deinit(plugin, plugin_deinitialize);

        if (result != CPL_ERROR_NONE)
        {
            cpl_plugin_delete(plugin);
            return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                         "Failed to create a copy of a plugin"
                                         " structure and append it to the"
                                         " output list.");
        }

        if (cpl_pluginlist_append(list, plugin) != CPL_ERROR_NONE)
        {
            cpl_plugin_delete(plugin);
            return cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                         "Failed to create a copy of a plugin"
                                         " structure and append it to the"
                                         " output list.");
        }
    }

    return CPL_ERROR_NONE;
}

/**
 * @brief Resets loaded Python module information.
 *
 * This function should be called when one no longer needs to deal with any
 * Python based plugins.
 */
void er_python_cleanup(void)
{
    current_module = NULL;

    /* Do not need to check for NULL since we set it to that value initially
       and cx_map_delete anyway checks for NULL. */
    cx_map_delete(module_lut);
    module_lut = NULL;
    cx_map_delete(module_list);
    module_list = NULL;
}

/**@}*/
