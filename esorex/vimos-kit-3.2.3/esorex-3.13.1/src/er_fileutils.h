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

#ifndef EP_FILEUTILS_H
#define EP_FILEUTILS_H


CPL_BEGIN_DECLS

const char *er_fileutils_tilde_replace(const char *);
const char *er_fileutils_dot_replace(const char *name);
int fileutils_directory_exists(const char *);
int fileutils_file_exists(const char *);
char *fileutils_create_fqfname(char *, char *);
char *fileutils_fqfname_filename(const char *);
char *fileutils_fqfname_dirname(const char *);
int fileutils_copy(const char *, const char *);
int fileutils_move(const char *, const char *);
int er_fileutils_link(const char *, const char *);
int er_fileutils_file_is_fits(const char * self);

#define  FILEMAX 4096

CPL_END_DECLS

#endif /* EP_FILEUTILS_H */
