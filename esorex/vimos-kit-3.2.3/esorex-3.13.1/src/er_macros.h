/*
 * This file is part of the ESO Recipe Execution Tool
 * Copyright (C) 2001-2017 European Southern Observatory
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

#ifndef ER_MACROS_H
#define ER_MACROS_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#define MAXSTRLENCONF     1024

#define ForceNull(x)  if (*x == '\0') x = NULL;

#if defined HAVE_DECL___FUNC__ && HAVE_DECL___FUNC__
#define er_func __func__
#else
#define er_func ""
#endif

/*
 * Global resource files. Uncomment to allow the automatic detection and 
 * use of global configuration files.
 */

#define GLOBAL_RC_EXTENSION     ".rc"
#define GLOBAL_RC_BACKUP        ".bak"
#define GLOBAL_RC_DIR           "." PACKAGE
#define GLOBAL_RC_NAME          PACKAGE GLOBAL_RC_EXTENSION

/*
 * Name of resource-prefix to be used in EsoRex configuration files.
 */

#define PACKAGE_RESOURCE   PACKAGE ".caller"

/* 
 * Prefix to use for environment variables containing configuration settings 
 */ 

#define PACKAGE_ENV     "ESOREX"

/*
 * Presentation definitions
 */

#define COMMENT_TAB_POSITION   24
#define COMMAND_LINE_PREFIX    "--"

/* needed because of the er_fileutils replace functions  */

#define PATHSET_MAX     1024

#endif /* ER_MACROS_H */
