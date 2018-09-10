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

#ifndef ER_STRINGARRAY_H
#define ER_STRINGARRAY_H

#include <cpl_macros.h>

CPL_BEGIN_DECLS

typedef struct er_stringarray_s er_stringarray_t;


extern er_stringarray_t *er_stringarray_new(void);
extern void er_stringarray_resize(er_stringarray_t *, int);
extern void er_stringarray_delete(er_stringarray_t *);
extern char *er_stringarray_get(er_stringarray_t *, int);
extern void er_stringarray_set(er_stringarray_t *, const char *, int);
extern void er_stringarray_remove(er_stringarray_t *, int);
extern int er_stringarray_present(er_stringarray_t *, int);
extern int er_stringarray_size(er_stringarray_t *);
extern void er_stringarray_print(er_stringarray_t *);
extern void er_stringarray_append(er_stringarray_t *, const char *);

CPL_END_DECLS

#endif /* ER_STRINGARRAY_H */
