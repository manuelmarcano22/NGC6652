/* $Id: fors_dfs.h,v 1.4 2013-04-24 14:05:15 cgarcia Exp $
 *
 * This file is part of the FORS Pipeline
 * Copyright (C) 2002-2006 European Southern Observatory
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

/*
 * $Author: cgarcia $
 * $Date: 2013-04-24 14:05:15 $
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VIMOS_DFS_H
#define VIMOS_DFS_H

#include <cpl.h>

CPL_BEGIN_DECLS

char *vmstrlower(char *s);
int dfs_save_image(cpl_frameset *, const cpl_image *, const char *, 
                   cpl_propertylist *, const cpl_parameterlist *, 
                   const char *, const char *);
int dfs_save_table(cpl_frameset *, const cpl_table *, const char *, 
                   cpl_propertylist *, const cpl_parameterlist *, 
                   const char *, const char *);
cpl_image *dfs_load_image(cpl_frameset *, const char *, cpl_type, int, int);
cpl_table *dfs_load_table(cpl_frameset *, const char *, int);
cpl_propertylist *dfs_load_header(cpl_frameset *, const char *, int);
int dfs_equal_keyword(cpl_frameset *frameset, const char *keyword);
int dfs_get_parameter_bool(cpl_parameterlist *, const char *, 
                           const cpl_table *);
int dfs_get_parameter_int(cpl_parameterlist *, const char *, 
                          const cpl_table *);
double dfs_get_parameter_double(cpl_parameterlist *, const char *, 
                                const cpl_table *);
const char *dfs_get_parameter_string(cpl_parameterlist *, const char *, 
                                     const cpl_table *);
int dfs_get_parameter_bool_const(const cpl_parameterlist *, const char *);
int dfs_get_parameter_int_const(const cpl_parameterlist *, const char *);
double dfs_get_parameter_double_const(const cpl_parameterlist *, const char *);
const char *dfs_get_parameter_string_const(const cpl_parameterlist *, const char *);
char *dfs_generate_filename_tfits(const char *);

cpl_error_code dfs_save_table_ext(cpl_table *, const char *,
				  cpl_propertylist *);
cpl_error_code dfs_save_image_ext(cpl_image *, const char *,
				  cpl_propertylist *);
cpl_error_code dfs_save_image_null(cpl_frameset *, 
                                   cpl_propertylist *, cpl_parameterlist *,
				   const char *, const char *,
				   const char *);

void vimos_dfs_set_groups(cpl_frameset * set);

cpl_frameset *
vimos_frameset_extract(const cpl_frameset *frames,
                       const char *tag);
CPL_END_DECLS

#endif   /* VIMOS_DFS_H */
