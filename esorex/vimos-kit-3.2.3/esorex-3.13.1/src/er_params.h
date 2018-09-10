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

#ifndef ER_PARAMS_H
#define ER_PARAMS_H

#include <cpl.h>

CPL_BEGIN_DECLS


extern int params_parse_config_file(cpl_parameterlist *, const char *);

extern int params_parse_config_environment(cpl_parameterlist *);

extern int params_parse_config_commandline(cpl_parameterlist *, char *,
                                           er_stringarray_t *, int,
                                           char **argv, int);

extern void params_parse_config_postprocess(cpl_parameterlist *);

extern int params_handle_parameters(char *, cpl_parameterlist *);

extern int params_process_configuration(cpl_parameterlist *, char *,
                                        char *, int, char **, char *,
                                        er_stringarray_t *);

CPL_END_DECLS

#endif
