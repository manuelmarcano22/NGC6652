/* $Id: piltask.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#define _XOPEN_SOURCE
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>
#ifdef HAVE_VFORK_H
# include <vfork.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/time.h>

#include "pilmemory.h"
#include "pilmessages.h"
#include "pilstrutils.h"
#include "pilfileutils.h"
#include "piltask.h"


/*
 * Status of the child process returned to the caller if
 * execve() failed. The same value as system() is used.
 */

#define PIL_PROC_CHILD_FAILED  127


/**
 * @defgroup pilTask pilTask
 *
 * The module @b pilTask provides functions for child process management.
 * Currently the provided functions are rather basic.
 */
 
/**@{*/

/*
 * This holds the process id pilTaskExecWait() is waiting for.
 * Since suspending caller execution until the child has finished
 * implies that there can only be one child forked at a time, a
 * global variable for the child pid is an appropritate implementation
 * here.
 */

static pid_t exec_wait_pid = -1;


/*
 * @brief
 *   Kill the child process waited for by @b pilTaskExecWait().
 *
 * @param signum  Signal number (not used).
 *
 * @return Nothing.
 *
 * Kill the child process with the pid, stored in the global variable
 * @b exec_wait_pid. This is the child process created by calling
 * @b pilTaskExecWait().
 */

static void exec_wait_kill_child(int signum)
{

    signum = SIGKILL;

    if (exec_wait_pid >= 0)
	kill(exec_wait_pid, signum);

    return;

}


/*
 * @brief
 *   Kill the child process waited for and then commit suicide.
 *
 * @param signum  Signal number (not used).
 *
 * @return Nothing.
 *
 * Kill the child process with the pid, stored in the global variable
 * @b exec_wait_pid. Then get the pid of the current process and kill it.
 */

static void exec_wait_kill_process(int signum)
{

    pid_t pid = getpid();


    signum = SIGKILL;

    if (exec_wait_pid >= 0)
	kill(exec_wait_pid, signum);

    kill(pid, signum);
    return;

}


/**
 * @brief
 *   Execute a command as a child process and wait for it.
 *
 * @param argc     Number of arguments in the list.
 * @param argv     Argument list specifying the command being executed
 *                 and its options.
 * @param timeout  Execution time limit for the child process in seconds.
 *
 * @return If the command launched as a child process exited normaly, the
 *   function returns the commands exit code. If launching the command
 *   failed, the function returns 127, like @b system(), does. If the
 *   command terminated because of receiving a signal, or if any other
 *   error occurs the return value is -1.
 *
 * The function launches the command provided by the argument list
 * @em argv as a child process and waits until the command is
 * finished. A execution time limit for the launched command can be
 * provided by @em timeout. In case the child exceeds this execution
 * time limit it is killed and the function returns an error. The timer
 * is not used if @em timeout is set to 0.
 * 
 * The layout of @em argv has to be the same as the argument list
 * which is received by a C @b main(), i.e @em argv[0] contains the
 * full path to the executable and the remaining elements contain
 * the command line arguments. Note that @em argv has to have one
 * more element than given by @em argc, since the last element must
 * be @c NULL.
 * 
 * An alternative to passing ordinary command line arguments (not the
 * executable path or the program name) as individual elements of the
 * argument list is to pass them together in a single element string.
 * 
 * The function removes the directory prefix from the executable path
 * and used only the filename as name for the child process.
 */

int pilTaskExecWait(int argc, const char *const argv[], time_t timeout)
{

    const char fctid[] = "pilTaskExecWait";


    char *cmd, **args;

    int i;
    int status;

    pid_t p;

    struct sigaction alarm, act, oact[9];

    struct stat sb;

    struct itimerval exec_timer, previous_timer;


    /*
     * Basic check of the argument vector. At least the command itself
     * should be present and it has to fit into the path buffer used by
     * some functions.
     */

    if (argc < 1 || !argv || !argv[0] || strlen(argv[0]) > PIL_PATHNAME_MAX)
	return -1;

    assert(argv[argc] == 0);


    /*
     * Check file permissions. Executables having the set UID bit and/or
     * set GID bit set are not executed.
     */

    if (access(argv[0], F_OK | X_OK))
	return -1;
    else {
	if (stat(argv[0], &sb) || !S_ISREG(sb.st_mode))
	    return -1;

	if (sb.st_mode & S_ISUID)
	    return -1;

	/*
	 * Take into account special uses for the set GID bit.
	 */

	if ((sb.st_mode & S_ISGID) && (sb.st_mode & S_IXGRP))
	    return -1;
    }

      
    /*
     * Setup the command and the argument list passed to execve
     */

    cmd = pil_strdup(argv[0]);

    if (!(args = (char **)pil_calloc(argc + 1, sizeof(char *)))) {
	pil_free(cmd);
	return -1;
    }

    args[0] = pil_strdup(pilFileBaseName(cmd));
    args[argc] = 0;

    for (i = 1; i < argc; i++)
	args[i] = pil_strdup(argv[i]);


    /*
     * Setup signal handlers used for killing the child process after a
     * timeout and any terminating signal received by the parent, during
     * execution of the child process.
     */

    alarm.sa_handler = exec_wait_kill_child;
    sigemptyset(&alarm.sa_mask);
    alarm.sa_flags = 0;

    act.sa_handler = exec_wait_kill_process;
    sigemptyset(&act.sa_mask);
    act.sa_flags = 0;


    /*
     * Setup the timer and make sure the process dies after exceeding
     * the time limit.
     */

    if (timeout > 0) {
	exec_timer.it_interval.tv_sec = timeout;
	exec_timer.it_interval.tv_usec = 0;
	exec_timer.it_value.tv_sec = timeout;
	exec_timer.it_value.tv_usec = 0;

	setitimer(ITIMER_REAL, &exec_timer, &previous_timer);

	/*
	 * Install the handler for the SIGALRM signal.
	 */
 
	sigaction(SIGALRM, &alarm, &oact[0]);
    }


    /*
     * Catch signals with a handler that kills the child process, if
     * there is one and then kills the current process.
     */

    sigaction(SIGHUP, &act, &oact[1]);
    sigaction(SIGINT, &act, &oact[2]);
    sigaction(SIGBUS, &act, &oact[3]);
    sigaction(SIGFPE, &act, &oact[4]);
    sigaction(SIGQUIT, &act, &oact[5]);
    sigaction(SIGABRT, &act, &oact[6]);
    sigaction(SIGTERM, &act, &oact[7]);
    sigaction(SIGSEGV, &act, &oact[8]);


    /*
     * Create the child process
     */

    exec_wait_pid = fork();

    switch (exec_wait_pid) {
	case -1: 

	    /* fork failed */

	    status = -1;
	    break;

	case 0:
 
	    /* Child side */

	    /*
	     * Restore signals
	     */

	    sigaction(SIGHUP, &oact[1], 0L);
	    sigaction(SIGINT, &oact[2], 0L);
	    sigaction(SIGBUS, &oact[3], 0L);
	    sigaction(SIGFPE, &oact[4], 0L);
	    sigaction(SIGQUIT, &oact[5], 0L);
	    sigaction(SIGABRT, &oact[6], 0L);
	    sigaction(SIGTERM, &oact[7], 0L);
	    sigaction(SIGSEGV, &oact[8], 0L);

	    execve(cmd, (char *const *)args, 0);

	    return PIL_PROC_CHILD_FAILED;
	    break;

	default:

	    /* Parent side */

	    /*
	     * Wait for the child to finish.
	     */

#ifndef HAVE_WAITPID
	    do {
		p = wait(&status);
	    } while (p != exec_wait_pid);
#else
	    p = waitpid(exec_wait_pid, &status, 0);
#endif

	    if (p != exec_wait_pid) {
		kill(exec_wait_pid, SIGKILL);
		status = -1;
	    }
	    else {
		if (WIFSIGNALED(status)) {
		    if (timeout > 0 && WTERMSIG(status) == SIGALRM)
			pilMsgDebug(fctid, "Execution time limit exceeded! "
				    "Process %d killed!", exec_wait_pid);
		    else 
			pilMsgDebug(fctid, "Process %d received signal %d. "
				    "Terminated!", exec_wait_pid,
				    WTERMSIG(status));
		    status = -1;
		}

		if (WIFEXITED(status))
		    status = WEXITSTATUS(status);
	    }
	    break;
    }


    /*
     * Reset the child pid.
     */

    exec_wait_pid = -1;


    /*
     * Restore the signals and the timer.
     */

    if (timeout > 0)
	sigaction(SIGALRM, &oact[0], 0L);

    sigaction(SIGHUP, &oact[1], 0L);
    sigaction(SIGINT, &oact[2], 0L);
    sigaction(SIGBUS, &oact[3], 0L);
    sigaction(SIGFPE, &oact[4], 0L);
    sigaction(SIGQUIT, &oact[5], 0L);
    sigaction(SIGABRT, &oact[6], 0L);
    sigaction(SIGTERM, &oact[7], 0L);
    sigaction(SIGSEGV, &oact[8], 0L);

    setitimer(ITIMER_REAL, &previous_timer, 0);


    /*
     * Cleanup
     */

    pil_free(cmd);

    for (i = 0; i < argc; i++)
	pil_free(args[i]);
    pil_free(args);

    return status;

}
/**@}*/
