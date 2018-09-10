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

#ifndef ER_STRINGUTILS_H
#define ER_STRINGUTILS_H

/* This constant is from cpl_msg.c */
#define MAX_MSG_LENGTH      (16384)
#define DEFAULT_TERM_WIDTH     (80)

CPL_BEGIN_DECLS

extern char * er_strutils_indent( const char *str, int );
extern char * er_strutils_split ( const char *, int, int );
extern int er_strutils_termwidth ( void );
extern int er_strutils_termheight ( void );
extern const char * er_strutils_dblstr( double );  
extern void er_strutils_scrtest(void);

CPL_END_DECLS

#endif /* ER_STRINGUTILS_H */
