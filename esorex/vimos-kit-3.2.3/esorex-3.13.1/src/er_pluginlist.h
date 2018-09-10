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

#ifndef ER_PLUGINLIST_H
#define ER_PLUGINLIST_H


CPL_BEGIN_DECLS

extern er_stringarray_t *er_pluginlist_create_list(char **listname);

extern int er_pluginlist_get_libpath(er_stringarray_t *list_of_pllibs,
                                     const char *plugin_name, char *libpath);

extern cpl_plugin *er_pluginlist_get_plugin(const char *library_name,
                                            const char *plugin_name,
                                            lt_dlhandle *module);

extern void er_pluginlist_create_cache(er_stringarray_t *nama,
                                       er_stringarray_t *namb);

extern void er_pluginlist_print_list(er_stringarray_t *list_of_pllibs);

CPL_END_DECLS

#endif /* ER_PLUGIN_LIST_H */
