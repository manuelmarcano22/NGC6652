/* $Id: pilmessages.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifndef _PIL_MESSAGES_H
#define _PIL_MESSAGES_H

#include <stdarg.h>

#include <pilmacros.h>
#include <pilutils.h>

 
PIL_BEGIN_DECLS

/*
 *  Constants
 */

#define MAX_MESSAGE_LENGTH  (1024)
#define MAX_FUNCTION_NAME   (30)
#define MAX_RECIPE_NAME     (40)
#define MAX_LOGFILE_NAME    (52)


/* 
 *  Strings used for the severity field in the message:
 */

#define ERROR_STRING   "[ERR] "
#define WARNING_STRING "[WAR] "
#define INFO_STRING    "[INF] "
#define DEBUG_STRING   "[DBG] "

/*
 *  Message severity
 */

typedef enum _PIL_MSG_SEVERITY_ {
  PIL_MSG_DEBUG   = 0,
  PIL_MSG_INFO,
  PIL_MSG_WARNING,
  PIL_MSG_ERROR,
  PIL_MSG_OFF
} PilMsgSeverity;

/*
 * Message handler
 */

typedef void (*PilPrintFunc)(const char *);


/*
 *  Methods
 */

int pilMsgStart(void);
void pilMsgStop(void);

int pilMsgEnableLog(PilMsgSeverity);
int pilMsgCloseLog(void);
const char *pilMsgGetLogFile(void);
void pilMsgEnableTerminal(PilMsgSeverity);
PilMsgSeverity pilMsgLogLevel(void);
PilMsgSeverity pilMsgTerminalLevel(void);

void pilMsgEnableTimeTag(void);
void pilMsgDisableTimeTag(void);
void pilMsgEnableComponentTag(void);
void pilMsgDisableComponentTag(void);
void pilMsgEnableRecipeTag(void);
void pilMsgDisableRecipeTag(void);
void pilMsgSetRecipeName(const char *);

void pilMsgSetWidth(int);
void pilMsgSetIndentStep(int);
void pilMsgIndentMore(void);
void pilMsgIndentLess(void);
void pilMsgIndent(int);

void pilMsgError(const char *functionName, char *format, ...);
void pilMsgWarning(const char *functionName, char *format, ...);
void pilMsgInfo(const char *functionName, char *format, ...);
void pilMsgDebug(const char *functionName, char *format, ...);

/*
 * Installable print handlers
 */

PilPrintFunc pilMsgSetPrintHandler(PilPrintFunc);
PilPrintFunc pilMsgSetErrorHandler(PilPrintFunc);

PIL_END_DECLS

#endif /* _PIL_MESSAGES_H */
