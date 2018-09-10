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

#ifndef EP_PARAM_UTILS_H
#define EP_PARAM_UTILS_H

#include <cpl.h>


CPL_BEGIN_DECLS

void er_paramutils_print_key_desc(const char *, const char *, const char *);
int er_paramutils_print_list(cpl_parameterlist*, const char *);
int paramutils_set_from_string(cpl_parameter *, const char *, const char *);
void er_paramutils_tilde_convert(cpl_parameterlist *);
int er_manage_sources(const int , const char *, char **);

CPL_END_DECLS

#endif /* EP_PARAM_UTILS_H */
