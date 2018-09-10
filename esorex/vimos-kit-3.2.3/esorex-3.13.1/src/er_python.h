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

#ifndef ER_PYTHON_H
#define ER_PYTHON_H

#include <cpl_pluginlist.h>
#include "er_stringarray.h"

CPL_BEGIN_DECLS

extern cpl_error_code er_python_load_modules(er_stringarray_t *);
extern cpl_error_code er_python_select_module(const char *);
extern int er_python_get_plugin_list(cpl_pluginlist *);
extern void er_python_cleanup(void);

CPL_END_DECLS

#endif /* ER_PYTHON_H */
