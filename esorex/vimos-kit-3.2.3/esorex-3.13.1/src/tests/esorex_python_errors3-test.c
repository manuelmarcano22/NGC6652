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
 * These are additional tests relying on overloading certain functions to mock
 * various failure modes that are difficult to trigger in any other way. These
 * tests are not compatible with other tests from esorex_python_errors*-test.c
 * because commonly use functions are being overloaded. Thus, they must be
 * compiled into their own binary.
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

#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <errno.h>
#include <cpl.h>
#include <cpl_test.h>

/*-----------------------------------------------------------------------------
                             Overload functions
 -----------------------------------------------------------------------------*/

/* Overload the following functions to allow producing errors on demand. */

cpl_boolean pipe_must_produce_error = CPL_FALSE;
int pipe_error_countdown = 0;

int _overloaded_pipe(int fildes[2])
{
    --pipe_error_countdown;
    if (pipe_must_produce_error && pipe_error_countdown < 0)
    {
        errno = EIO;
        return -1;
    }
    else
    {
        return pipe(fildes);
    }
}

cpl_boolean fork_must_produce_error = CPL_FALSE;

pid_t _overloaded_fork(void)
{
    if (fork_must_produce_error)
    {
        errno = EIO;
        return -1;
    }
    else
    {
        return fork();
    }
}

cpl_boolean dup_must_produce_error = CPL_FALSE;
int dup_error_countdown = 0;

int _overloaded_dup(int oldfd)
{
    --dup_error_countdown;
    if (dup_must_produce_error && dup_error_countdown < 0)
    {
        errno = EIO;
        return -1;
    }
    else
    {
        return dup(oldfd);
    }
}

cpl_boolean dup2_must_produce_error = CPL_FALSE;
int dup2_error_countdown = 0;

int _overloaded_dup2(int oldfd, int newfd)
{
    --dup2_error_countdown;
    if (dup2_must_produce_error && dup2_error_countdown < 0)
    {
        errno = EIO;
        return -1;
    }
    else
    {
        return dup2(oldfd, newfd);
    }
}

cpl_boolean close_must_produce_error = CPL_FALSE;
int close_failure_code = EIO;
int close_error_countdown = 0;
int close_on_which_pid = -1;

int _overloaded_close(int fd)
{
    if (getpid() != close_on_which_pid)
    {
        return close(fd);
    }
    --close_error_countdown;
    if (close_must_produce_error && close_error_countdown < 0)
    {
        errno = close_failure_code;
        return -1;
    }
    else
    {
        return close(fd);
    }
}

cpl_boolean read_must_produce_error = CPL_FALSE;
int read_failure_code = EIO;
int read_call_countdown = 0;
int read_nbyte_limit = -1;

ssize_t _overloaded_read(int fildes, void *buf, size_t nbyte)
{
    --read_call_countdown;
    if (read_must_produce_error || read_call_countdown >= 0)
    {
        errno = read_failure_code;
        return -1;
    }
    else
    {
        if (read_nbyte_limit > 0)
        {
            return read(fildes, buf, read_nbyte_limit);
        }
        else
        {
            return read(fildes, buf, nbyte);
        }
    }
}

cpl_boolean write_must_produce_error = CPL_FALSE;
int write_failure_code = EIO;
int write_call_countdown = 0;
int write_nbyte_limit = -1;

ssize_t _overloaded_write(int fildes, const void * buf, size_t nbyte)
{
    --write_call_countdown;
    if (write_must_produce_error || write_call_countdown >= 0)
    {
        errno = write_failure_code;
        return -1;
    }
    else
    {
        if (write_nbyte_limit > 0)
        {
            return write(fildes, buf, write_nbyte_limit);
        }
        else
        {
            return write(fildes, buf, nbyte);
        }
    }
}

cpl_boolean malloc_must_produce_error = CPL_FALSE;

void * _overloaded_malloc(size_t size)
{
    if (malloc_must_produce_error)
    {
        return NULL;
    }
    else
    {
        return malloc(size);
    }
}

cpl_boolean realloc_must_produce_error = CPL_FALSE;

void * _overloaded_realloc(void *ptr, size_t size)
{
    if (realloc_must_produce_error)
    {
        return NULL;
    }
    else
    {
        return realloc(ptr, size);
    }
}

/* Overload waitpid to simulate a signal interruption, and delay when the real
   waitpid is called by a certain number of invocations. */
int waitpid_call_countdown = 0;

pid_t _overloaded_waitpid(pid_t pid, int * wstatus, int options)
{
    --waitpid_call_countdown;
    if (waitpid_call_countdown >= 0)
    {
        errno = EINTR;
        return -1;
    }
    else
    {
        return waitpid(pid, wstatus, options);
    }
}

/* Replace calls in the code we want to test with our overloaded functions. */
#define pipe _overloaded_pipe
#define fork _overloaded_fork
#define dup _overloaded_dup
#define dup2 _overloaded_dup2
#define close _overloaded_close
#define read _overloaded_read
#define write _overloaded_write
#define malloc _overloaded_malloc
#define realloc _overloaded_realloc
#define waitpid _overloaded_waitpid
#include "er_python.c"
#undef pipe
#undef fork
#undef dup
#undef dup2
#undef close
#undef read
#undef write
#undef malloc
#undef realloc
#undef waitpid

/*-----------------------------------------------------------------------------
                        Private function prototypes
 -----------------------------------------------------------------------------*/

static void mock_python(const char * script);
static void test_start_python_interpreter(void);
static void test_run_python_command(void);
static void test_stop_python_interpreter(void);
static void test_pipe_read_write(void);

/*-----------------------------------------------------------------------------
                                  Main
 -----------------------------------------------------------------------------*/

int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    /* Update the PATH environment variable to make sure we pick up our mock
       python script. */
    const char * path = "./mock_python_error3_test:/usr/bin:/bin";
    cpl_test_assert(setenv("PATH", path, 0x1) == 0);

    /* Need to restrict overloading close() to this process. */
    close_on_which_pid = getpid();

    /* Insert tests below */
    test_start_python_interpreter();
    test_run_python_command();
    test_stop_python_interpreter();
    test_pipe_read_write();

    /* End of tests */
    return cpl_test_end(0);
}

/*-----------------------------------------------------------------------------
                             Private functions
 -----------------------------------------------------------------------------*/

/*
 * This function is used to create a shell script that will mock up a Python
 * interpreter process. The input is the shell script to write to the python
 * script file.
 */
static void mock_python(const char * script)
{
    /* Do not overwrite the python script if a previous test failed. This will
       make it easier to debug by having the failing script available. */
    cpl_test_assert(cpl_test_get_failed() == 0);

    /* Create the directory for the python script. We put this in its own
       directory so that this unit test can be made concurrency safe. */
    if (mkdir("./mock_python_error3_test", 0777) != 0)
    {
        cpl_test_assert(errno == EEXIST);
    }

    FILE* file = fopen("./mock_python_error3_test/python", "w");
    cpl_test_assert(file != NULL);
    size_t bytes_to_write = strlen(script);
    size_t bytes_written = fwrite(script, sizeof(char), bytes_to_write, file);
    cpl_test_assert(bytes_written == bytes_to_write);
    int fclose_result = fclose(file);
    cpl_test_assert(fclose_result == 0);
    mode_t mode = S_IRUSR | S_IWUSR | S_IXUSR |
                  S_IRGRP | S_IXGRP |
                  S_IROTH | S_IXOTH;
    cpl_test_assert(chmod("./mock_python_error3_test/python", mode) == 0);
}

static void test_start_python_interpreter(void)
{
    pid_t pid = -1;
    int input = -1;
    int output = -1;

    close_must_produce_error = CPL_FALSE;
    read_must_produce_error = CPL_FALSE;
    write_must_produce_error = CPL_FALSE;
    malloc_must_produce_error = CPL_FALSE;
    realloc_must_produce_error = CPL_FALSE;

    /* Test error handling for edge case where pipe() fails. Note that pipe() is
       called twice, so check both cases. */
    pipe_must_produce_error = CPL_TRUE;
    pipe_error_countdown = 0;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(start_python_interpreter("x=1", &pid, &input, &output),
                      CPL_ERROR_FILE_IO);

    pipe_must_produce_error = CPL_TRUE;
    pipe_error_countdown = 1;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(start_python_interpreter("x=1", &pid, &input, &output),
                      CPL_ERROR_FILE_IO);

    /* Test error handling for edge case where fork() fails. */
    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_TRUE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(start_python_interpreter("x=1", &pid, &input, &output),
                      CPL_ERROR_FILE_IO);

    /* Test edge case where dup() fails. In this case start_python_interpreter
       will not be able to detect this since it happens in the child process.
       Note that dup() is called twice, so check both cases. */
    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_TRUE;
    dup_error_countdown = 0;
    dup2_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(start_python_interpreter("x=1", &pid, &input, &output),
                      CPL_ERROR_NONE);

    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_TRUE;
    dup_error_countdown = 1;
    dup2_must_produce_error = CPL_FALSE;
    cpl_test_eq_error(start_python_interpreter("x=1", &pid, &input, &output),
                      CPL_ERROR_NONE);

    /* Test edge case where dup2() fails. As for dup(), start_python_interpreter
       will not be able to detect this since it happens in the child process. */
    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_TRUE;
    dup2_error_countdown = 0;
    cpl_test_eq_error(start_python_interpreter("x=1", &pid, &input, &output),
                      CPL_ERROR_NONE);

    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_TRUE;
    dup2_error_countdown = 1;
    cpl_test_eq_error(start_python_interpreter("x=1", &pid, &input, &output),
                      CPL_ERROR_NONE);
}

static void test_run_python_command(void)
{
    close_must_produce_error = CPL_FALSE;
    read_must_produce_error = CPL_FALSE;
    write_must_produce_error = CPL_FALSE;
    malloc_must_produce_error = CPL_FALSE;
    realloc_must_produce_error = CPL_FALSE;

    /* Test error handling for edge case where pipe() fails. */
    pipe_must_produce_error = CPL_TRUE;
    pipe_error_countdown = 0;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_FALSE;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);

    /* Test edge case where dup() fails. In this case run_python_command should
       be able to detect the error. Note that dup() is called twice, so check
       both cases. */
    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_TRUE;
    dup_error_countdown = 0;
    dup2_must_produce_error = CPL_FALSE;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);

    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_TRUE;
    dup_error_countdown = 1;
    dup2_must_produce_error = CPL_FALSE;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);

    /* Test also edge cases for dup2() failures. */
    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_TRUE;
    dup2_error_countdown = 0;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);

    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_TRUE;
    dup2_error_countdown = 1;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);

    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_FALSE;

    /* Mockup the python command that will just echo the input to output. */
    mock_python(
        "#!/bin/sh\n"
        "cat <&"CPL_STRINGIFY(PYTHON_INPUT_FD)
            " >&"CPL_STRINGIFY(PYTHON_OUTPUT_FD)"\n"
    );

    /* Simulate a failure in the write() function and check error handling. */
    write_must_produce_error = CPL_TRUE;
    write_failure_code = EIO;
    write_call_countdown = 0;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);
    write_must_produce_error = CPL_FALSE;

    /* Simulate a failure in the read() function and check error handling. */
    read_must_produce_error = CPL_TRUE;
    read_failure_code = EIO;
    read_call_countdown = 0;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);
    read_must_produce_error = CPL_FALSE;

    /* Simulate a failure in the close() function called in run_python_command
       and check error handling. */
    close_must_produce_error = CPL_TRUE;
    close_failure_code = EIO;
    close_error_countdown = 2;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_FILE_IO);
    close_must_produce_error = CPL_FALSE;

    /* Simulate an interrupt in the close() function and check that NULL is
       returned. No error code should be set in this case however. */
    close_must_produce_error = CPL_TRUE;
    close_failure_code = EINTR;
    close_error_countdown = 2;
    cpl_test_null(run_python_command("x=1", "input data"));
    cpl_test_error(CPL_ERROR_NONE);
    close_must_produce_error = CPL_FALSE;
}

static void test_stop_python_interpreter(void)
{
    close_must_produce_error = CPL_FALSE;
    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_FALSE;
    read_must_produce_error = CPL_FALSE;
    write_must_produce_error = CPL_FALSE;
    malloc_must_produce_error = CPL_FALSE;
    realloc_must_produce_error = CPL_FALSE;

    /* Test edge case where waitpid does not return immediately with the status
       of the exited child, but is interrupted. */
    waitpid_call_countdown = 3;
    int input = dup(0);
    int output = dup(1);
    pid_t pid = fork();
    if (pid == 0)
    {
        exit(EXIT_SUCCESS);
    }
    else
    {
        cpl_test_eq_error(stop_python_interpreter(pid, input, output),
                          CPL_ERROR_NONE);
    }

    /* Test edge case where the child is terminated forcefully with a signal,
       i.e. abnormal termination and we should have WIFEXITED(status) != true.
     */
    waitpid_call_countdown = 0;
    input = dup(0);
    output = dup(1);
    pid = fork();
    if (pid == 0)
    {
        sleep(10);
        exit(EXIT_FAILURE);
    }
    else
    {
        cpl_test_assert(kill(pid, SIGTERM) == 0);
        cpl_test_eq_error(stop_python_interpreter(pid, input, output),
                          CPL_ERROR_FILE_IO);
    }
}

static void test_pipe_read_write(void)
{
    int fildes[2];
    const char * input = "some data";
    char * output = NULL;
    char * hugeinput = NULL;
    int hugesize = 0;
    int n = 0;
    pid_t pid = -1;
    int status = -1;

    close_must_produce_error = CPL_FALSE;
    pipe_must_produce_error = CPL_FALSE;
    fork_must_produce_error = CPL_FALSE;
    dup_must_produce_error = CPL_FALSE;
    dup2_must_produce_error = CPL_FALSE;
    read_must_produce_error = CPL_FALSE;
    write_must_produce_error = CPL_FALSE;
    malloc_must_produce_error = CPL_FALSE;
    realloc_must_produce_error = CPL_FALSE;
    read_nbyte_limit = -1;
    write_nbyte_limit = -1;

    /* Prepare a 32 MB input data buffer for some of the following tests. */
    hugesize = 32*1024*1024;
    hugeinput = cpl_malloc(hugesize);
    const char * pattern = "0123456789abcdef";
    for (n = 0; n < hugesize; ++n)
    {
        hugeinput[n] = pattern[n % 16];
    }
    hugeinput[hugesize-1] = '\0';

    /* Test that writing to a Unix pipe succeeds. */
    cpl_test_assert(pipe(fildes) == 0);
    cpl_test_eq_error(write_to_pipe(fildes[1], input, strlen(input)),
                      CPL_ERROR_NONE);
    (void) close(fildes[1]);

    /* Test that reading from the other end succeeds and that we get the same
       result as what was written to the pipe. */
    output = read_from_pipe(fildes[0]);
    cpl_test_nonnull(output);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq_string(input, output);
    free(output);
    (void) close(fildes[0]);

    /* Test a read write cycle if we only allow one byte to be written or read
       at a time. */
    read_nbyte_limit = 1;
    write_nbyte_limit = 1;
    cpl_test_assert(pipe(fildes) == 0);
    cpl_test_eq_error(write_to_pipe(fildes[1], input, strlen(input)),
                      CPL_ERROR_NONE);
    (void) close(fildes[1]);
    output = read_from_pipe(fildes[0]);
    cpl_test_nonnull(output);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq_string(input, output);
    free(output);
    (void) close(fildes[0]);
    read_nbyte_limit = -1;
    write_nbyte_limit = -1;

    /* Do the read write test with 2 bytes at a time. */
    read_nbyte_limit = 2;
    write_nbyte_limit = 2;
    cpl_test_assert(pipe(fildes) == 0);
    cpl_test_eq_error(write_to_pipe(fildes[1], input, strlen(input)),
                      CPL_ERROR_NONE);
    (void) close(fildes[1]);
    output = read_from_pipe(fildes[0]);
    cpl_test_nonnull(output);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq_string(input, output);
    free(output);
    (void) close(fildes[0]);
    read_nbyte_limit = -1;
    write_nbyte_limit = -1;

    /* Test error handling if write() fails. */
    cpl_test_assert(pipe(fildes) == 0);
    write_must_produce_error = CPL_TRUE;
    write_failure_code = EIO;
    write_call_countdown = 0;
    cpl_test_eq_error(write_to_pipe(fildes[1], input, strlen(input)),
                      CPL_ERROR_FILE_IO);
    (void) close(fildes[1]);
    write_must_produce_error = CPL_FALSE;

    /* Test error handling if read() fails. */
    read_must_produce_error = CPL_TRUE;
    read_failure_code = EIO;
    read_call_countdown = 0;
    cpl_test_null(read_from_pipe(fildes[0]));
    cpl_test_error(CPL_ERROR_FILE_IO);
    (void) close(fildes[0]);
    read_must_produce_error = CPL_FALSE;

    /* Simulate signal interruptions in write() calls and see that writing to
       the pipe still completes correctly. */
    cpl_test_assert(pipe(fildes) == 0);
    write_must_produce_error = CPL_FALSE;
    write_failure_code = EINTR;
    write_call_countdown = 2;
    cpl_test_eq_error(write_to_pipe(fildes[1], input, strlen(input)),
                      CPL_ERROR_NONE);
    (void) close(fildes[1]);
    write_must_produce_error = CPL_FALSE;

    /* Simulate interruptions for read() and see that read_from_pipe still
       completed successfully. */
    read_must_produce_error = CPL_FALSE;
    read_failure_code = EINTR;
    read_call_countdown = 2;
    output = read_from_pipe(fildes[0]);
    cpl_test_nonnull(output);
    cpl_test_error(CPL_ERROR_NONE);
    cpl_test_eq_string(input, output);
    free(output);
    (void) close(fildes[0]);
    read_must_produce_error = CPL_FALSE;

    /* Test error handling in read_from_pipe if malloc() fails. */
    cpl_test_assert(pipe(fildes) == 0);
    cpl_test_assert(write_to_pipe(fildes[1], input, strlen(input))
                    == CPL_ERROR_NONE);
    (void) close(fildes[1]);
    malloc_must_produce_error = CPL_TRUE;
    cpl_test_null(read_from_pipe(fildes[0]));
    cpl_test_error(CPL_ERROR_FILE_IO);
    (void) close(fildes[0]);
    malloc_must_produce_error = CPL_FALSE;

    /* Test error handling in read_from_pipe if realloc() fails. Note that
       we need to check the edge case when realloc is called at the end of the
       read_from_pipe function to shrink wrap the buffer space allocated to
       exactly the number of bytes read, and also for huge input when it is
       called to increase the internal buffer. */
    cpl_test_assert(pipe(fildes) == 0);
    cpl_test_assert(write_to_pipe(fildes[1], input, strlen(input))
                    == CPL_ERROR_NONE);
    (void) close(fildes[1]);
    realloc_must_produce_error = CPL_TRUE;
    cpl_test_null(read_from_pipe(fildes[0]));
    cpl_test_error(CPL_ERROR_FILE_IO);
    (void) close(fildes[0]);
    realloc_must_produce_error = CPL_FALSE;

    cpl_test_assert(pipe(fildes) == 0);
    /* NOTE: We have to run in a separate process to avoid a deadlock, since
       something has to read from the pipe so that it does not go full. */
    pid = fork();
    if (pid == 0)
    {
        (void) close(fildes[0]);
        write_to_pipe(fildes[1], hugeinput, strlen(hugeinput));
        (void) close(fildes[1]);
        exit(EXIT_SUCCESS);
    }
    else
    {
        (void) close(fildes[1]);
        realloc_must_produce_error = CPL_TRUE;
        cpl_test_null(read_from_pipe(fildes[0]));
        cpl_test_error(CPL_ERROR_FILE_IO);
        (void) close(fildes[0]);
        realloc_must_produce_error = CPL_FALSE;
        wait(&status);
    }

    /* Test that processing a large volume of data through the pipe produces
       correct output. */
    cpl_test_assert(pipe(fildes) == 0);
    pid = fork();
    if (pid == 0)
    {
        (void) close(fildes[0]);
        write_to_pipe(fildes[1], hugeinput, strlen(hugeinput));
        (void) close(fildes[1]);
        exit(EXIT_SUCCESS);
    }
    else
    {
        (void) close(fildes[1]);
        output = read_from_pipe(fildes[0]);
        cpl_test_nonnull(output);
        cpl_test_error(CPL_ERROR_NONE);
        cpl_test_eq_string(hugeinput, output);
        free(output);
        (void) close(fildes[0]);
        wait(&status);
    }

    cpl_free(hugeinput);
}


#else /* ENABLE_PYTHON_RECIPES */

int main(void)
{
    /* Indicate to the automake test runner that the tests are skipped. */
    return 77;
}

#endif  /* ENABLE_PYTHON_RECIPES */
