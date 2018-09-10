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

#ifndef ER_PLUGIN_H
#define ER_PLUGIN_H

#include <cxstrutils.h>

#include <cpl.h>

#include <er_stringarray.h>


CPL_BEGIN_DECLS

extern int plugin_process_plugin(cpl_parameterlist *, char *,
                                 er_stringarray_t *, int, char **);

extern cpl_msg_severity message_severity(cpl_parameterlist *param_list,
                                         int flag);

extern void er_enlarge(const char *fn, char **pptr, int msize);

extern int mysscanf(char *myline, char *path, char *tag, char *group);

extern cpl_frameset *er_frameset_load(const char *name, cpl_frameset *set,
                                      int check_input_files_flag);

CPL_END_DECLS

#endif  /* ER_PLUGIN_H */
