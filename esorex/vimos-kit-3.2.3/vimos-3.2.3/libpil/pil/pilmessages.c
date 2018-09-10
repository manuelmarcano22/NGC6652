/* $Id: pilmessages.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
#include <stdio.h>
#include <string.h>
#include <unistd.h>
/*
#include <termio.h>
*/
#include <signal.h>
/*
#include <stropts.h>
*/
#include <time.h>

#include "pildate.h"
#include "pilmessages.h"

/*
 *  This is the max length for text lines that are written to the logfile.
 *  It is also the max length for text lines sent to the terminal, in case 
 *  the window size cannot be determined by the appropriate call to ioctl().
 */

#define PIL_TERM_WIDTH (80)

static struct sigaction act, oact;


/**
 * @defgroup pilMessages pilMessages
 *
 * The module @b pilMessages provides functions to display and
 * logging messages from the DRS.
 */

/**@{*/

/* FIXME:
 * The following private members might be organized in one (or more)
 * structures, in case a complete description of the messaging system
 * configuration is requested by the client.
 */

static PilMsgSeverity logMinLevel                     = PIL_MSG_OFF;
static PilMsgSeverity terminalMinLevel                = PIL_MSG_ERROR;
static int            timeTag                         = 0;
static int            recipeTag                       = 1;
static int            componentTag                    = 0;

static char           recipeName[MAX_RECIPE_NAME]     = "Undefined";
static char           logfileName[MAX_LOGFILE_NAME]   = ".logfile";
static FILE          *logfilePointer                  = NULL;

static int            pageWidth                       = PIL_TERM_WIDTH;

/*
 * Enable this if you want split lines in log file:
 *
 * static int            logWidth                        = PIL_TERM_WIDTH;
 */

static int            indentStep                      = 2;
static int            indentValue                     = 0;

static FILE          *msgStdout;
static FILE          *msgStderr;
static int            outStream;
static int            errStream;


/*
 * Installable print handlers
 */

static PilPrintFunc pil_message_printer = NULL;
static PilPrintFunc pil_error_printer = NULL;

/*
 * To save the default handlers at messaging startup, so that they
 * can be restored when the system is shut down.
 */

static PilPrintFunc _pil_message_printer = NULL;
static PilPrintFunc _pil_error_printer = NULL;


 static void
_pil_message_print(const char *message)
{

    fputs(message, msgStdout);
    fflush(msgStdout);

    return;

}


 static void
_pil_error_print(const char *message)
{

    fputs(message, msgStderr);
    fflush(msgStderr);

    return;

}


static void
_pil_print(const char *format, ...)
{

    va_list args;
    char string[MAX_MESSAGE_LENGTH + 1];
    PilPrintFunc printer;


    if (format == NULL) {
        return;
    }

    va_start(args, format);

#ifdef HAVE_VSNPRINTF
    vsnprintf(string, (size_t) MAX_MESSAGE_LENGTH, format, args);
#else
    vsprintf(string, format, args);
#endif
    
    string[MAX_MESSAGE_LENGTH] = '\0';

    printer = pil_message_printer;

    if (printer != NULL) {
        printer(string);
    }
    else {
        fputs(string, stdout);
        fflush(stdout);
    }

    return;

}


static void
_pil_printerr(const char *format, ...)
{

    va_list args;
    char string[MAX_MESSAGE_LENGTH + 1];
    PilPrintFunc printer;


    if (format == NULL) {
        return;
    }

    va_start(args, format);

#ifdef HAVE_VSNPRINTF
    vsnprintf(string, (size_t) MAX_MESSAGE_LENGTH, format, args);
#else
    vsprintf(string, format, args);
#endif
    
    string[MAX_MESSAGE_LENGTH] = '\0';

    printer = pil_error_printer;

    if (printer != NULL) {
        printer(string);
    }
    else {
        fputs(string, stderr);
        fflush(stderr);
    }

    return;

}


/*
 * @brief
 *  Signal handler for signal @c SIGWINCH
 *
 * @param i  Dummy argument (not used!)
 *
 * @return Nothing.
 *
 * The function accomodates the output line width of the messaging
 * subsystem to the new window size on arrival of the signal @c SIGWINCH.
 */

/*
 static void pilMsgChangeWidth(int i)
{
  struct winsize win;

  i = 0;

  if (ioctl(outStream, TIOCGWINSZ, &win) < 0 || win.ws_col < 1)
    pageWidth = PIL_TERM_WIDTH;
  else
    pageWidth = win.ws_col;
}
*/


/*
 * @brief
 *   Format and output message string.
 *
 * @param severity      Severity level of the incoming message. 
 * @param functionName  Name of the function generating the message.
 * @param format        Format string in the usual C convention.
 * @param al            Variable argument list associated to the format string.
 *
 * @return Nothing.
 *
 * This function is used internally to actually display/add the
 * message to terminal and/or logfile. Messages with severity 
 * level equal to "warning" or greater would be sent to stderr, 
 * the other messages would go to stdout.
 *
 * If the severity level is lower than the levels set by 
 * @b pilMsgEnableTerminal() and @b pilMsgEnableLog(), then the message
 * is not displayed.
 *
 * @see pilMsgEnableTerminal(), pilMsgEnableLog()
 */

 static void pilMsgOut(PilMsgSeverity severity,
			     const char *functionName, char *format,
			     va_list al)
{
  time_t seconds;
  char   messageText[MAX_MESSAGE_LENGTH];
  char   messageLog[MAX_MESSAGE_LENGTH];
  char   messageTerminal[MAX_MESSAGE_LENGTH];
/*
 * Enable this if you want split lines in log file:
 * int    startLogLine;
 */
  int    startTerminalLine;
  int    i;


  if (severity < terminalMinLevel && severity < logMinLevel) return;
  if (severity == PIL_MSG_OFF) return;

  seconds = time((time_t *)0);
   
#ifdef HAVE_VSNPRINTF
  vsnprintf(messageText, (size_t) MAX_MESSAGE_LENGTH, format, al);
#else
  vsprintf(messageText, format, al);
#endif


  /* 
   *  For safety, truncate the message to MAX_MESSAGE_LENGTH 
   */

  messageText[MAX_MESSAGE_LENGTH-1] = '\0';


  /* 
   *  Date and time. Note that time tag and severity field are not 
   *  affected by indentation. Date and time are always present in 
   *  logfile, optional in terminal output.
   */

  strftime(messageLog, MAX_MESSAGE_LENGTH, "%H:%M:%S ", localtime(&seconds));
  if (timeTag) {
    strftime(messageTerminal,MAX_MESSAGE_LENGTH,"%H:%M:%S ",
	     localtime(&seconds));
  }
  else {
    messageTerminal[0] = '\0';
  }


  /*
   *  Severity label
   */

  if (severity == PIL_MSG_ERROR) {
    strcat(messageLog, ERROR_STRING);
    strcat(messageTerminal, ERROR_STRING);
  }
  else if (severity == PIL_MSG_WARNING) {
    strcat(messageLog, WARNING_STRING);
    strcat(messageTerminal, WARNING_STRING);
  }
  else if (severity == PIL_MSG_INFO) {
    strcat(messageLog, INFO_STRING);
    strcat(messageTerminal, INFO_STRING);
  }
  else if (severity == PIL_MSG_DEBUG) {
    strcat(messageLog, DEBUG_STRING);
    strcat(messageTerminal, DEBUG_STRING);
  }


  /*
   *  Recipe, function name, and message appended:
   */

  if (recipeTag) {
    strcat(messageTerminal, recipeName);
    strcat(messageTerminal, ": ");
  }

/*
 * Enable this if you want split lines in log file:
 *
 * startLogLine = strlen(messageLog);
 */

  startTerminalLine = strlen(messageTerminal);

  if (componentTag) {
    strcat(messageTerminal, functionName);
    strcat(messageTerminal, "()  ");
  }

  strcat(messageLog, functionName);
  strcat(messageLog, "()  ");

  /*
   *  Message indentation
   */

  for (i=0; i<indentValue; i++) {
    strcat(messageLog, " ");
    strcat(messageTerminal, " ");
  }

  strcat(messageLog, messageText);
  strcat(messageTerminal, messageText);

  if (severity >= logMinLevel) {
/*
 * Enable this (and disable the next line) if you want split lines in log file:
 *
 *  fprintf(logfilePointer, "%s\n", strsplit(messageLog, startLogLine,
 *                                           logWidth));
 */
    fprintf(logfilePointer, "%s\n", messageLog);
  }

  if (severity >= terminalMinLevel) {
    if (severity > PIL_MSG_WARNING) {
      _pil_printerr("%s\n", strsplit(messageTerminal, startTerminalLine,
                                     pageWidth));
    }
    else {
      _pil_print("%s\n", strsplit(messageTerminal, startTerminalLine,
                                  pageWidth));
    }
  }
}


/**
 * @brief
 *   Initialize the messaging system
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE
 *
 * Currently just the terminal width is determined (if possible), 
 * and the resized window signal handler is deployed.
 */

int pilMsgStart(void)
{
/*
  struct winsize win;
*/

  /*
   *  First duplicate stdout and stderr streams
   */

  if (!(outStream = dup(fileno(stdout))))
    return EXIT_FAILURE;

  if (!(errStream = dup(fileno(stderr))))
    return EXIT_FAILURE;

  if (!(msgStdout = fdopen(outStream, "a")))
    return EXIT_FAILURE;

  if (!(msgStderr = fdopen(errStream, "a")))
    return EXIT_FAILURE;

  _pil_message_printer = pilMsgSetPrintHandler(_pil_message_print);
  _pil_error_printer = pilMsgSetErrorHandler(_pil_error_print);


  /*
   *  Get the terminal window size, and if successful deploy the handler 
   *  for any image resizing at runtime.
   */

/*

  if (ioctl(outStream, TIOCGWINSZ, &win) < 0 || win.ws_col < 1)
    return EXIT_SUCCESS;

  pageWidth = win.ws_col;

  act.sa_handler = pilMsgChangeWidth;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;             * Probably more appropriate flags         *
                                * initialization should be inserted here. * 

  act.sa_flags &= ~SA_SIGINFO;  * Eliminates SA_SIGINFO from any setting  *
                                * above.                                  * 

  sigaction(SIGWINCH, &act, &oact);

*/

  return EXIT_SUCCESS;
}


/**
 * @brief
 *   Shutdown the messaging system
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE
 *
 * The original SIGWINCH signal handler is restored.
 */

void pilMsgStop(void)
{
/*
  if (act.sa_handler == pilMsgChangeWidth) 
    sigaction(SIGWINCH, &oact, NULL);
*/

  fclose(msgStdout);
  fclose(msgStderr);

  pilMsgSetPrintHandler(_pil_message_printer);
  pilMsgSetErrorHandler(_pil_error_printer);

  return;

}


/**
 * @brief
 *   Open and initialize the log file.
 *
 * @param severity  Specification of the lowest severity level a message
 *                  should have to be written to log file.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE
 *
 * The logfile name is the same as the recipe name, with
 * a time tag and the extension .log appended. Tipically
 * this function would be called at pipeline initialization,
 * when a given recipe is run. Calling it again would
 * close the previous logfile, and start a new one with
 * the same name and a different time tag.
 */

int pilMsgEnableLog(PilMsgSeverity severity)
{
  char *timeLabel;

  if (logfilePointer) {

    /*
     *  If a log file was already open, close it before initializing
     *  a new one.
     */

    if (pilMsgCloseLog() == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
  }
  if (severity != PIL_MSG_OFF) {
    logMinLevel = severity;

    if ((logfilePointer = fopen(logfileName, "w")) == NULL) {
      return EXIT_FAILURE;
    }

    /*
     *  Write log file header
     */

    timeLabel = pilDateGetISO8601();
    fprintf(logfilePointer, "\n");
    fprintf(logfilePointer, "Start time     : %s\n", timeLabel);
    fprintf(logfilePointer, "Recipe name    : %s\n", recipeName);
    fprintf(logfilePointer, "Severity level : ");
    switch(severity) {
      case PIL_MSG_DEBUG   : fprintf(logfilePointer, DEBUG_STRING); break;
      case PIL_MSG_INFO    : fprintf(logfilePointer, INFO_STRING); break;
      case PIL_MSG_WARNING : fprintf(logfilePointer, WARNING_STRING); break;
      case PIL_MSG_ERROR   : fprintf(logfilePointer, ERROR_STRING); break;
      default              : ;
    }
    fprintf(logfilePointer, "\n\n");
  }
  return EXIT_SUCCESS;
}


/**
 * @brief
 *   Close the current log file.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * Typically this routine should be called at pipeline end.
 */

int pilMsgCloseLog(void)
{
  if (logMinLevel != PIL_MSG_OFF) {
    logMinLevel = PIL_MSG_OFF;
    if (fclose(logfilePointer)) return EXIT_FAILURE;
    logfilePointer = NULL;
  }
  return EXIT_SUCCESS;
}


/**
 * @brief
 *   Get the logfile name.
 *
 * @return Logfile name
 *
 * Read access to internal static string.
 */

const char *pilMsgGetLogFile(void)
{
  return logfileName;
}


/**
 * @brief
 *   Enable output to terminal.
 *
 * @param severity  Specification of the lowest severity level a message
 *                  should have to be displayed to terminal.
 *
 * @return Nothing.
 *
 * Typically this function would be called at pipeline
 * initialization, when a given recipe is run. Calling it
 * again would have no effect but changing the severity 
 * level.
 */

void pilMsgEnableTerminal(PilMsgSeverity severity)
{
  terminalMinLevel = severity;

  if (severity == PIL_MSG_DEBUG) pilMsgEnableComponentTag();
}


/**
 * @brief
 *   Get current log messaging level.
 *
 * @return Current severity level.
 *
 * Get current log messaging level.
 */

PilMsgSeverity pilMsgLogLevel(void)
{
  return logMinLevel;
}


/**
 * @brief
 *   Get current terminal messaging level.
 *
 * @return Current severity level.
 *
 * Get current terminal messaging level.
 */

PilMsgSeverity pilMsgTerminalLevel(void) 
{
  return terminalMinLevel;
}


/**
 * @brief
 *   Enable a Time Tag to output messages.
 *
 * @return Nothing.
 *
 * There is no analogous routine for the log file, as a
 * time tag would be always added to a message written
 * to the log file. If this routine is not called, the
 * Time Tag is not used.
 */

void pilMsgEnableTimeTag(void)
{
  timeTag = 1;
}


/**
 * @brief
 *   Disable Time Tag to output messages.
 *
 * @return Nothing.
 *
 * There is no analogous routine for the log file, as a
 * time tag would be always added to a message written
 * to the log file.
 */

void pilMsgDisableTimeTag(void)
{
  timeTag = 0;
}


/**
 * @brief
 *   Display Recipe Name before output messages.
 *
 * @return Nothing.
 *
 * There is no analogous routine for the log file, as a
 * recipe name would be always added to a message written
 * to the log file.
 */

void pilMsgEnableRecipeTag(void)
{
  recipeTag = 1;
}


/**
 * @brief
 *   Disable display of Recipe Name before output messages.
 *
 * @return Nothing.
 *
 * There is no analogous routine for the log file, as a
 * recipe name would be always added to a message written
 * to the log file. If this routine is not called, the
 * recipe name is displayed before messages.
 */

void pilMsgDisableRecipeTag(void)
{
  recipeTag = 0;
}


/**
 * @brief
 *   Display Component Name (i.e., function) before output messages.
 *
 * @return Nothing.
 *
 * There is no analogous routine for the log file, as a
 * component name would be always added to a message written
 * to the log file.
 */

void pilMsgEnableComponentTag(void)
{
  componentTag = 1;
}


/**
 * @brief
 *   Disable display of Component Name before output messages.
 *
 * @return Nothing.
 *
 * There is no analogous routine for the log file, as a
 * recipe name would be always added to a message written
 * to the log file. If this routine is not called, the
 * recipe name is displayed before messages.
 */

void pilMsgDisableComponentTag(void)
{
  componentTag = 0;
}


/**
 * @brief
 *   Communicate the recipe name to the messaging system.
 *
 * @param name Any task identification string, typically the
 *   DRS recipe name..
 *
 * @return Nothing.
 *
 * This routine should be called at the recipe start,
 * before a possible call to the function pilMsgEnableLog. 
 * The recipe name will be displayed at the beginning 
 * of any message to terminal, and in the header of 
 * the logfile.
 */

void pilMsgSetRecipeName(const char *name)
{
  if (strlen(name) >= MAX_RECIPE_NAME) {
    strncpy(recipeName, name, MAX_RECIPE_NAME);
    recipeName[MAX_RECIPE_NAME-1] = '\0';
  }
  else {
    strcpy(recipeName, name);
  }
}


/**
 * @brief
 *   Set the max width of the displayed text.
 *
 * @param width Max width of the displayed text.
 *
 * @return Nothing.
 *
 * If a message is longer than this width, it would be broken into 
 * shorter lines before display to terminal. This function is called 
 * by the messaging system every time the terminal window is resized,
 * and the width is set equal to the new width of the terminal window.
 */

void pilMsgSetWidth(int width)
{
  pageWidth = width;
}


/**
 * @brief
 *   Set the indentation step.
 *
 * @param step Indentation step.
 *
 * @return Nothing.
 *
 * To maintain a constant message display style, this routine
 * should be called at most once, and just at pipeline start.
 * A message line might be moved leftward or rightward by a
 * number of characters that is a multiple of the indentation 
 * step. Setting the indentation step to zero would eliminate 
 * messages indentation. If this function is not called, the 
 * indentation step is set to 2.
 */

void pilMsgSetIndentStep(int step)
{
  indentStep = step;
}


/**
 * @brief
 *   Set the indentation level.
 *
 * @param level Indentation level.
 *
 * @return Nothing.
 *
 * Any new message line will be indented by a number of character
 * equal to the indentation level mutiplied by the indentation
 * step.
 */

void pilMsgIndent(int level)
{
  indentValue = level*indentStep;
}


/**
 * @brief
 *   Increase messages indentation by one indentation step.
 *
 * @return Nothing.
 *
 * Calling this function is equivalent to increase the indentation
 * level by 1.
 */

void pilMsgIndentMore(void)
{
  indentValue += indentStep;
}


/**
 * @brief
 *   Decrease messages indentation by one indentation step.
 *
 * @return Nothing.
 *
 * Calling this function is equivalent to decrease the indentation
 * level by 1.
 */

void pilMsgIndentLess(void)
{
  indentValue -= indentStep;
}


/**
 * @brief
 *   Display an error message.
 *
 * @param functionName  Name of the function generating the message.
 * @param format        Format string.
 * @param ...           Variable argument list associated to the format string.
 *
 * @return Nothing.
 *
 * Any new message line will be indented by a number of character
 * equal to the indentation level mutiplied by the indentation
 * step.
 *
 * The format string @em format should follow the usual C convention. New
 * line characters should not be used, as the message would be split
 * automatically according to the width specified with @ b pilMsgSetWidth().
 * Inserting a new line character would enforce breaking a line of 
 * text even before the current row is filled, and should be indented 
 * only as the begin of a new paragraph of text.
 *
 * @see pilMsgSetWidth()
 */

void pilMsgError(const char *functionName, char *format, ...)
{
  va_list al;

  va_start(al, format);
  pilMsgOut(PIL_MSG_ERROR, functionName, format, al);
  va_end(al);
}


/**
 * @brief
 *   Display a warning message.
 *
 * @param functionName  Name of the function generating the message.
 * @param format        Format string.
 * @param ...           Variable argument list associated to the format string.
 *
 * @return Nothing.
 *
 * Any new message line will be indented by a number of character
 * equal to the indentation level mutiplied by the indentation
 * step.
 *
 * The format string @em format should follow the usual C convention. New
 * line characters should not be used, as the message would be split
 * automatically according to the width specified with @ b pilMsgSetWidth().
 * Inserting a new line character would enforce breaking a line of 
 * text even before the current row is filled, and should be indented 
 * only as the begin of a new paragraph of text.
 *
 * @see pilMsgSetWidth()
 */

void pilMsgWarning(const char *functionName, char *format, ...)
{
  va_list al;

  va_start(al, format);
  pilMsgOut(PIL_MSG_WARNING, functionName, format, al);
  va_end(al);
}


/**
 * @brief
 *   Display an informative message.
 *
 * @param functionName  Name of the function generating the message.
 * @param format        Format string.
 * @param ...           Variable argument list associated to the format string.
 *
 * @return Nothing.
 *
 * Any new message line will be indented by a number of character
 * equal to the indentation level mutiplied by the indentation
 * step.
 *
 * The format string @em format should follow the usual C convention. New
 * line characters should not be used, as the message would be split
 * automatically according to the width specified with @ b pilMsgSetWidth().
 * Inserting a new line character would enforce breaking a line of 
 * text even before the current row is filled, and should be indented 
 * only as the begin of a new paragraph of text.
 *
 * @see pilMsgSetWidth()
 */

void pilMsgInfo(const char *functionName, char *format, ...)
{
  va_list al;

  va_start(al, format);
  pilMsgOut(PIL_MSG_INFO, functionName, format, al);
  va_end(al);
}


/**
 * @brief
 *   Display a debug message.
 *
 * @param functionName  Name of the function generating the message.
 * @param format        Format string.
 * @param ...           Variable argument list associated to the format string.
 *
 * @return Nothing.
 *
 * Any new message line will be indented by a number of character
 * equal to the indentation level mutiplied by the indentation
 * step.
 *
 * The format string @em format should follow the usual C convention. New
 * line characters should not be used, as the message would be split
 * automatically according to the width specified with @ b pilMsgSetWidth().
 * Inserting a new line character would enforce breaking a line of 
 * text even before the current row is filled, and should be indented 
 * only as the begin of a new paragraph of text.
 *
 * @see pilMsgSetWidth()
 */

void pilMsgDebug(const char *functionName, char *format, ...)
{
  va_list al;

  va_start(al, format);
  pilMsgOut(PIL_MSG_DEBUG, functionName, format, al);
  va_end(al);
}


/**
 * @brief
 *   Set handler for message output.
 *
 * @param func New handler function.
 *
 * @return The previously set print handler.
 *
 * The function @em func is installed as the new message printing function.
 * Any ordinary message is printed using this handler. The default print
 * handler just outputs the message text to @c stdout.
 */

PilPrintFunc
pilMsgSetPrintHandler(PilPrintFunc func)
{

    PilPrintFunc previous;

    previous = pil_message_printer;
    pil_message_printer = func;

    return previous;

}


/**
 * @brief
 *   Set handler for error message output.
 *
 * @param func New handler function.
 *
 * @return The previously set error message handler.
 *
 * The function @em func is installed as the new error message printing
 * function. Any error message is printed using this handler. The default
 * print handler just outputs the error message text to @c stderr.
 */

PilPrintFunc
pilMsgSetErrorHandler(PilPrintFunc func)
{

    PilPrintFunc previous;

    previous = pil_error_printer;
    pil_error_printer = func;

    return previous;

}
/**@}*/
