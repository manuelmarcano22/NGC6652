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

#ifndef ER_PAF_H
#define ER_PAF_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cpl.h>


struct er_paf_config_item
{
    char             *pro_catg_regexp;
    er_stringarray_t *requested_keywords;
};

struct er_paf_config
{
    struct er_paf_config_item *items;
    unsigned int               nitems;
};

extern int er_create_recipe_pafs(const cpl_frameset *frames,
                                 const char *recipe_name,
                                 const char *paf_config_filename);

extern int er_paf_single_hdu_fill_write(const cpl_frame *frame,
                                        struct er_paf_config_item *config_item,
                                        int idx_paf,
                                        const char *recipe_name,
                                        char *arcfile);

extern int er_paf_multi_hdu_fill_write(const cpl_frame *frame,
                                       struct er_paf_config_item *config_item,
                                       int idx_paf,
                                       const char *recipe_name,
                                       char *arcfile);

extern struct er_paf_config *er_paf_read_config_file
(
    const char  *val_paf_config
);

extern cpl_propertylist *
er_paf_fill_paf_proplist_single_hdu(const cpl_frame *frame,
                                    struct er_paf_config_item *config_item,
                                    char *arcfile);

extern cpl_propertylist *
er_paf_fill_paf_proplist_multi_hdu(const cpl_frame *frame,
                                   struct er_paf_config_item *config_item,
                                   int iext,
                                   char *arcfile);

extern int er_paf_parse_config(er_stringarray_t *file_stripped,
                               struct er_paf_config *paf_configuration);

extern struct er_paf_config_item *
er_paf_get_matching_item(struct er_paf_config *paf_configuration,
                         const char *this_pro_catg);

extern er_stringarray_t *
er_paf_read_whole_file(const char *val_paf_config);

extern int er_paf_config_delete(struct er_paf_config *paf_configuration);

extern char *er_paf_get_arcfile(const cpl_frameset *frames);

#endif /* ER_PAF_H*/
