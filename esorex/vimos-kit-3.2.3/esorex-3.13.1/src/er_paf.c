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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <regex.h>

#include "cxstrutils.h"
#include "cpl_dfs.h"

#include "er_macros.h"
#include "er_stringarray.h"
#include "er_fileutils.h"
#include "er_paf.h"

#define ER_PAF_MANDATORY_KEYS_MAIN_EXT  "PRO CATG|INSTRUME|MJD-OBS|ESO TPL ID|DATE-OBS|ESO OBS ID" 
#define ER_PAF_MANDATORY_KEYS_EXT  "EXTNAME" 
#define ER_PAF_QC_KEYS  "ESO QC" 

/**
 * @defgroup  er_paf  PAF files handling
 *
 * This module provides a collection of functions that deal
 * with the PAF files creation  
 *
 */

/**@{*/

/**********************************************************************/
/**
 * @brief Creates the PAF files
 *
 * @param frames        The frames to create the PAF from
 * @param recipe_name   The name of the corresponding recipe
 * @param paf_config_filename    The configuration file for PAF creation
 *
 * @returns 0 if successful, !=0 otherwise
 *
 * Create a PAF file using the important keywords from the header
 * of the products. The important keywords are specified in the configuration
 * file.
 *
 */
/**********************************************************************/

int er_create_recipe_pafs
(const cpl_frameset * frames,
 const char *         recipe_name,
 const char *         paf_config_filename)
{
    const cpl_frame *      frame = NULL;
    struct er_paf_config * paf_configuration = NULL;
    int                    paf_err = CPL_ERROR_NONE;
    char *                 arcfile = NULL;
    int                    idx_paf = 0;

    cpl_msg_info(er_func, "Creating PAF files");

    /* Read the configuration file */
    paf_configuration = er_paf_read_config_file(paf_config_filename);
    if(paf_configuration == NULL)
    {
        cpl_msg_error(er_func, "Could not parse the PAF configuration file");
        return CPL_ERROR_BAD_FILE_FORMAT;
    }

    /* Loop into all the products */
    cpl_msg_debug(__func__,"Looping into all the products");

    cpl_frameset_iterator *it = cpl_frameset_iterator_new(frames);

    while ((frame = cpl_frameset_iterator_get_const(it)) != NULL)
    {
        const char * filename;
        filename = cpl_frame_get_filename(frame);
        if (cpl_frame_get_group(frame) == CPL_FRAME_GROUP_PRODUCT &&
                er_fileutils_file_is_fits(filename))
        {
            struct er_paf_config_item * config_item = NULL;
            const char *                this_pro_catg = NULL;
            cpl_propertylist *          procatg_proplist = NULL;

            /* Get PRO CATG of this product */
            procatg_proplist = cpl_propertylist_load(filename, 0);
            this_pro_catg = cpl_propertylist_get_string(procatg_proplist, CPL_DFS_PRO_CATG);

            if(this_pro_catg == NULL)
            {
                cpl_msg_error(er_func,"Could not get PRO CATG in file %s",
                              filename);
                cpl_frameset_iterator_delete(it);
                er_paf_config_delete(paf_configuration);
                cpl_free(arcfile);
                paf_err = cpl_error_get_code();
                return paf_err;
            }

            ++idx_paf;
            cpl_msg_debug(er_func, "Product with PRO CATG=%s", this_pro_catg);
            cpl_msg_indent_more();

            /* Get the corresponding matching in the config file (if there is) */
            config_item = er_paf_get_matching_item
                    (paf_configuration, this_pro_catg);
            if(config_item != NULL)
            {
                int next;

                cpl_msg_debug(er_func, "Matched with regexp =%s",
                              config_item->pro_catg_regexp);

                /*
                 *  Try to get ARCFILE from the product itself. If it is not present
                 *  fall back to retrieving it from the first raw file in the frame
                 *  set.
                 */

                arcfile = cpl_strdup(cpl_propertylist_get_string(procatg_proplist,
                                                                 "ESO PRO ANCESTOR"));

                if(arcfile == NULL)
                {
                    arcfile = er_paf_get_arcfile(frames);
                    if(arcfile == NULL)
                    {
                        cpl_msg_error(er_func,
                                      "Could not read ARCFILE information from "
                                      "raw frames");
                        cpl_frameset_iterator_delete(it);
                        er_paf_config_delete(paf_configuration);
                        return CPL_ERROR_DATA_NOT_FOUND;
                    }
                }

                /* Decide what to do depending on number of extensions */
                next = cpl_frame_get_nextensions(frame);
                if(next == 0)
                {
                    paf_err = er_paf_single_hdu_fill_write
                            (frame, config_item, idx_paf, recipe_name, arcfile);
                }
                else if(next != 0)
                {
                    paf_err = er_paf_multi_hdu_fill_write
                            (frame, config_item, idx_paf, recipe_name, arcfile);
                }
                if(paf_err != CPL_ERROR_NONE)
                {
                    cpl_propertylist_delete(procatg_proplist);
                    cpl_frameset_iterator_delete(it);
                    er_paf_config_delete(paf_configuration);
                    cpl_free(arcfile);
                    return paf_err;
                }
            }
            cpl_msg_indent_less();
            cpl_propertylist_delete(procatg_proplist);
        }
        else {
            cpl_msg_debug(__func__,"File is not a FITS product %s",
                          filename);
        }

        cpl_frameset_iterator_advance(it, 1);

    }

    cpl_frameset_iterator_delete(it);

    /* Cleanup */
    er_paf_config_delete(paf_configuration);
    if(arcfile != NULL)
        cpl_free(arcfile);

    return paf_err;
}

/**********************************************************************/
/**
 * @brief Creates a list of keywords and saves the corresponding PAF file
 *        for single HDU products.
 *
 * @param frame     The frame product where to extract the keywords from
 * @param config_item   The configuration contains the keywords requested
 * @param idx_paf       The index number for this PAF file (it will be added to
 *          the filename
 * @param recipe_name   The recipe name (this is saved in the PAF file)    
 * @param arcfile   This is added as the ARCFILE keyword
 * @returns CPL_ERROR_NONE if success or the error code otherwise
 *
 */
/**********************************************************************/

int er_paf_single_hdu_fill_write
(const cpl_frame *           frame,
 struct er_paf_config_item * config_item,
 int                         idx_paf,
 const char *                recipe_name,
 char *                      arcfile)
{
    cpl_propertylist *  paf_proplist = NULL;
    char *              paf_filename;

    /* Get the keywords to write */
    paf_proplist = 
            er_paf_fill_paf_proplist_single_hdu(frame, config_item, arcfile);
    if(paf_proplist == NULL)
    {
        cpl_msg_debug(__func__,"No keywords found for file %s. "
                      "Skipping writing of paf", cpl_frame_get_filename(frame));
        return CPL_ERROR_NONE;
    }

    /* Get the PAF filename */
    paf_filename = cpl_sprintf("qc_%04d.paf",idx_paf);

    /* Writing PAF file */
    cpl_dfs_save_paf(cpl_propertylist_get_string
                     (paf_proplist,"INSTRUME"),
                     recipe_name,
                     paf_proplist,
                     paf_filename);

    /* Cleaning */
    cpl_free(paf_filename);
    cpl_propertylist_delete(paf_proplist);

    return CPL_ERROR_NONE;
}

/**********************************************************************/
/**
 * @brief Creates a list of keywords and saves the corresponding PAF file
 *        for multiple HDU products.
 *
 * @param frame     The frame product where to extract the keywords from
 * @param config_item   The configuration contains the keywords requested
 * @param idx_paf       The index number for these PAF files (it will be added 
 *          to the filenames)
 * @param recipe_name   The recipe name (this is saved in the PAF file)    
 * @param arcfile   This is added as the ARCFILE keyword
 * @returns CPL_ERROR_NONE if success or the error code otherwise
 *
 */
/**********************************************************************/

int er_paf_multi_hdu_fill_write
(const cpl_frame *           frame,
 struct er_paf_config_item * config_item,
 int                         idx_paf,
 const char *                recipe_name,
 char *                      arcfile)
{
    int                 iext;
    int                 next;

    next = cpl_frame_get_nextensions(frame);
    for(iext = 0 ; iext <= next ; ++iext)
    {
        cpl_propertylist *  paf_proplist;
        char *              paf_filename;

        /* Get the keywords to write */
        paf_proplist =
                er_paf_fill_paf_proplist_multi_hdu(frame,config_item, iext, arcfile);
        if(paf_proplist == NULL)
        {
            cpl_msg_warning(__func__,"Problem reading headers in file %s, extension %d",
                            cpl_frame_get_filename(frame), iext);
            return CPL_ERROR_ILLEGAL_INPUT;
        }
        if(cpl_propertylist_get_size(paf_proplist) == 0)
        {
            cpl_msg_debug(__func__,"Skipping writing of paf for file %s, "
                          "extension %d (starting in 0).",
                          cpl_frame_get_filename(frame), iext);
            cpl_propertylist_delete(paf_proplist);
            continue;
        }

        /* Get the PAF filename */
        paf_filename = cpl_sprintf("qc_%04d_%04d.paf",idx_paf, iext);

        /* Writing PAF file */
        cpl_dfs_save_paf(cpl_propertylist_get_string
                         (paf_proplist,"INSTRUME"),
                         recipe_name,
                         paf_proplist,
                         paf_filename);
        /* Cleaning */
        cpl_free(paf_filename);
        cpl_propertylist_delete(paf_proplist);
    }

    return CPL_ERROR_NONE;
}


/**********************************************************************/
/**
 * @brief Get the matching configuration section 
 *
 * @param paf_configuration    The whole configuration
 * @param this_pro_catg        The procatg to match
 * @returns a pointer to a er_paf_config_item structure if there is
 *  a match, NULL otherwise
 *
 * This function loops through the config pro catg regular expressions and
 * return the first match 
 */
/**********************************************************************/
struct er_paf_config_item * er_paf_get_matching_item
(struct er_paf_config * paf_configuration, const char * this_pro_catg)
{
        struct er_paf_config_item * config_item = NULL;
        unsigned int                item;

        cpl_msg_indent_more();
        /* Loop through the config pro catg regular expressions */
        for(item = 0; item < paf_configuration->nitems; ++item)
        {
            regex_t pattern;
            int     regstatus;

            cpl_msg_debug(er_func, "Comparing PRO CATG with %s",
                          paf_configuration->items[item].pro_catg_regexp);

            regstatus = regcomp(&pattern,
                                paf_configuration->items[item].pro_catg_regexp,
                                REG_EXTENDED|REG_NOSUB);
            if(regstatus)
            {
                cpl_msg_error(er_func, "Bad formatted PRO CATG regular expression");
                cpl_error_set(er_func, CPL_ERROR_ILLEGAL_INPUT);
                cpl_msg_indent_less();
                regfree(&pattern);
                return NULL;
            }
            regstatus =
                    regexec(&pattern, this_pro_catg, (size_t)0, NULL, 0);
            if (regstatus != REG_NOMATCH)
            {
                config_item = paf_configuration->items+item;
                cpl_msg_indent_less();
                regfree(&pattern);
                return config_item;
            }
            regfree(&pattern);
        }
        cpl_msg_indent_less();

        /* Return */
        return config_item;
}


/**********************************************************************/
/**
 * @brief Reads the configuration to write PAF files
 *
 * @param paf_config_filename    Path to the configuration file
 * @returns a pointer to a er_paf_config structure if successfull,
 *  NULL otherwise
 *
 * This function reads the configuration file for PAF writing functionality.
 * It returns an allocated er_paf_config object, which must be deallocated
 * with er_paf_config_delete.
 */
/**********************************************************************/

struct er_paf_config * er_paf_read_config_file
(const char * paf_config_filename)
{
        struct er_paf_config * paf_configuration;
        er_stringarray_t * file_buffer;
        er_stringarray_t * file_stripped;
        int iline;

        cpl_msg_debug(__func__,"Reading configuration file %s",paf_config_filename);

        /* Allocate and init the structure */
        paf_configuration = cpl_malloc(sizeof(struct er_paf_config));
        paf_configuration->items = NULL;
        paf_configuration->nitems = 0;

        /* Read the whole file */
        file_buffer = er_paf_read_whole_file(paf_config_filename);
        if(file_buffer == NULL)
            return NULL;

        /* Strip the comments and blank lines */
        cpl_msg_debug(__func__,"Stripping comments and blank lines");
        file_stripped = er_stringarray_new();
        for (iline = 0; iline < er_stringarray_size(file_buffer); ++iline)
        {
            char * line;
            line = er_stringarray_get(file_buffer, iline);
            cx_strstrip(line);
            if (strcmp(line, "") != 0 && line[0] != '#' )
                er_stringarray_append(file_stripped, line);
        }

        /* Parse the configuration lines and store it in a convenient structure */
        er_paf_parse_config(file_stripped, paf_configuration);

        /* clean and return */
        er_stringarray_delete(file_buffer);
        er_stringarray_delete(file_stripped);
        return paf_configuration;
}

/**********************************************************************/
/**
 * @brief Creates a propertylist with all keywords to be saved for single
 *        extension products
 *
 * @param frame     The frame product where to extract the keywords from
 * @param config_item   The configuration contains the keywords requested
 * @param arcfile   This is added as the ARCFILE keyword
 * @returns a pointer to the propertylist if sucess,
 *  NULL otherwise
 *
 * This function will read all the keywords that are going to be
 * saved in the PAF file. It includes the mandatory keywords, the QC 
 * keywords and the specified keywords in the config file and the ARCFILE
 * keyword which comes from a external frame (the first raw frame from the 
 * frameset).
 * This function should be used only for single extension products  
 */
/**********************************************************************/

cpl_propertylist* er_paf_fill_paf_proplist_single_hdu
(const cpl_frame *           frame,
 struct er_paf_config_item * config_item,
 char *                      arcfile)
{
    cpl_propertylist * paf_proplist;
    cpl_propertylist * all_keywords;
    const char *       filename;
    int                ikey;
    int                nkeys;

    /* Get All keywords from FITS header */
    filename = cpl_frame_get_filename(frame);
    all_keywords = cpl_propertylist_load(filename, 0);
    paf_proplist = cpl_propertylist_new();
    if(all_keywords == NULL)
    {
        cpl_msg_error(er_func,"Could not read the main header of %s", filename);
        return NULL;
    }

    /* Retrieve mandatory keywords */
    cpl_propertylist_copy_property_regexp
    (paf_proplist, all_keywords, ER_PAF_MANDATORY_KEYS_MAIN_EXT, 0);
    cpl_propertylist_append_string(paf_proplist, "ARCFILE", arcfile);

    /* Retrieve requested keywords in configuration */
    nkeys = er_stringarray_size(config_item->requested_keywords);
    for(ikey = 0 ; ikey < nkeys; ++ikey)
    {
        char * key = er_stringarray_get(config_item->requested_keywords, ikey);
        cpl_propertylist *  matched_keys;

        matched_keys = cpl_propertylist_new();
        cpl_propertylist_copy_property_regexp(matched_keys, all_keywords,key,0);
        if(cpl_propertylist_get_size(matched_keys) == 0)
            cpl_msg_debug(er_func, "Cannot find keywords matching %s", key);
        cpl_propertylist_append(paf_proplist, matched_keys);
        cpl_propertylist_delete(matched_keys);
    }

    /* Retrieve QC keywords */
    cpl_propertylist_copy_property_regexp
    (paf_proplist, all_keywords, ER_PAF_QC_KEYS, 0);

    /* Bad return if there are no keywords at all */
    if(cpl_propertylist_get_size(paf_proplist) == 0)
    {
        cpl_msg_error(er_func, "Could not get any keyword from file %s",
                      filename);
        cpl_propertylist_delete(paf_proplist);
        cpl_propertylist_delete(all_keywords);
        return NULL;
    }

    /* Clean and return */
    cpl_propertylist_delete(all_keywords);
    return paf_proplist;
}

/**********************************************************************/
/**
 * @brief Creates a propertylist with all keywords to be saved for multi
 *        extension products
 *
 * @param frame     The frame product where to extract the keywords from
 * @param config_item   The configuration contains the keywords requested
 * @param iext      The number of the extension to get the keywords from
 * @param arcfile   This is added as the ARCFILE keyword
 * @returns a pointer to the propertylist if success,
 *  NULL otherwise
 *
 * This function will read all the keywords that are going to be
 * saved in the PAF file. It includes the mandatory keywords, the QC 
 * keywords and the specified keywords in the config file and the ARCFILE
 * keyword which comes from a external frame (the first raw frame from the 
 * frameset). All the keywords are retrieved from both the main extension
 * and the requested extension. If there are duplicated keywords in the 
 * main header, the extension ones will take precedence 
 * 
 * This function should be used only for multi extension products  
 */
/**********************************************************************/

cpl_propertylist* er_paf_fill_paf_proplist_multi_hdu
(const cpl_frame *           frame,
 struct er_paf_config_item * config_item,
 int                         iext,
 char *                      arcfile)
{
    cpl_propertylist * paf_proplist;
    cpl_propertylist * all_keywords_hdu_main;
    cpl_propertylist * all_keywords_hdu_ext;
    cpl_propertylist * qc_keywords;
    const char *       filename;
    int                ikey;
    int                nkeys;

    /* Read All keywords from FITS headers */
    paf_proplist = cpl_propertylist_new();
    filename = cpl_frame_get_filename(frame);
    all_keywords_hdu_main = cpl_propertylist_load(filename, 0);
    if(all_keywords_hdu_main == NULL)
    {
        cpl_msg_error(er_func,"Could not read the main header of %s", filename);
        return NULL;
    }
    all_keywords_hdu_ext = cpl_propertylist_load(filename, iext);
    if(all_keywords_hdu_ext == NULL)
    {
        cpl_msg_error(er_func,"Could not read header of extension %d in %s",
                      iext, filename);
        return NULL;
    }

    /* Retrieve mandatory keywords from main header */
    cpl_propertylist_copy_property_regexp
    (paf_proplist, all_keywords_hdu_main,
     ER_PAF_MANDATORY_KEYS_MAIN_EXT, 0);
    cpl_propertylist_append_string(paf_proplist, "ARCFILE", arcfile);

    /* Retrieve mandatory keywords from extension */
    cpl_propertylist_copy_property_regexp
    (paf_proplist, all_keywords_hdu_ext, ER_PAF_MANDATORY_KEYS_EXT, 0);

    /* Retrieve requested keywords in configuration */
    nkeys = er_stringarray_size(config_item->requested_keywords);
    for(ikey = 0 ; ikey < nkeys; ++ikey)
    {
        char * key = er_stringarray_get(config_item->requested_keywords, ikey);
        cpl_propertylist *  matched_keys;

        matched_keys = cpl_propertylist_new();
        cpl_propertylist_copy_property_regexp(matched_keys,
                                              all_keywords_hdu_main,key,0);
        cpl_propertylist_copy_property_regexp(matched_keys,
                                              all_keywords_hdu_ext,key,0);
        if(cpl_propertylist_get_size(matched_keys) == 0)
            cpl_msg_debug(er_func, "Cannot find keywords matching %s", key);
        cpl_propertylist_append(paf_proplist, matched_keys);
        cpl_propertylist_delete(matched_keys);
    }

    /* Retrieve QC keywords */
    qc_keywords = cpl_propertylist_new();
    cpl_propertylist_copy_property_regexp
    (qc_keywords, all_keywords_hdu_ext, ER_PAF_QC_KEYS, 0);
    if(cpl_propertylist_get_size(qc_keywords) == 0)
    {
        cpl_msg_debug(er_func, "No QC keywords found in file %s, extension %d",
                      filename, iext);
        cpl_propertylist_delete(paf_proplist);
        cpl_propertylist_delete(all_keywords_hdu_main);
        cpl_propertylist_delete(all_keywords_hdu_ext);
        return qc_keywords;
    }
    cpl_propertylist_append(paf_proplist, qc_keywords);
    cpl_propertylist_delete(qc_keywords);

    /* Bad return if there are no keywords at all */
    if(cpl_propertylist_get_size(paf_proplist) == 0)
    {
        cpl_msg_warning(er_func, "Could not get any keyword from file %s",
                        filename);
        cpl_propertylist_delete(all_keywords_hdu_main);
        cpl_propertylist_delete(all_keywords_hdu_ext);
        return paf_proplist;
    }

    /* Clean and return */
    cpl_propertylist_delete(all_keywords_hdu_main);
    cpl_propertylist_delete(all_keywords_hdu_ext);
    return paf_proplist;
}

char * er_paf_get_arcfile(const cpl_frameset * frames)
{
    const cpl_frame *  frame;
    char *             arcfile = NULL;
    cpl_propertylist * all_keywords;
    const char *       filename;

    /* Loop into all the products */

    cpl_frameset_iterator *it = cpl_frameset_iterator_new(frames);

    while ((frame = cpl_frameset_iterator_get_const(it)) != NULL)
    {
        if (cpl_frame_get_group(frame) == CPL_FRAME_GROUP_RAW)
        {
            filename = cpl_frame_get_filename(frame);
            all_keywords = cpl_propertylist_load(filename, 0);
            if(all_keywords == NULL)
            {
                cpl_frameset_iterator_delete(it);
                cpl_msg_error(er_func,"Could not read the main header of %s",
                              filename);
                return NULL;
            }

            const char *_arcfile = cpl_propertylist_get_string(all_keywords,
                                                               "ESO PRO ANCESTOR");

            if (!_arcfile && !cpl_propertylist_has(all_keywords, "ESO PRO CATG")) {
                _arcfile = cpl_propertylist_get_string(all_keywords, "ARCFILE");
            }

            if (!_arcfile) {
                cpl_msg_error(er_func,"Could not get ARCFILE information from "
                              "the main header of %s", filename);
                cpl_frameset_iterator_delete(it);
                cpl_propertylist_delete(all_keywords);
                return NULL;
            }

            arcfile = cpl_strdup(_arcfile);

            cpl_propertylist_delete(all_keywords);
            cpl_frameset_iterator_delete(it);

            return arcfile;
        }

        cpl_frameset_iterator_advance(it, 1);

    }

    cpl_frameset_iterator_delete(it);

    cpl_msg_error(er_func,"There are no RAW frames");
    return arcfile;

}


/**********************************************************************/
/**
 * @brief Reads all the file in a string buffer
 *
 * @param paf_config_filename  Path to the configuration file
 * @returns a pointer to the lines buffer if successfull, NULL otherwise
 *
 */
/**********************************************************************/
er_stringarray_t * er_paf_read_whole_file(const char * paf_config_filename)
{
    FILE * file_descriptor;
    char line[MAXSTRLENCONF];
    er_stringarray_t * file_buffer;

    /* Allocate string array */
    file_buffer = er_stringarray_new();

    /*  Go through all the file */
    file_descriptor = fopen(paf_config_filename, "r");
    if (file_descriptor == NULL)
        return NULL;
    memset(line, 0, MAXSTRLENCONF);

    while (fgets(line, MAXSTRLENCONF -1, file_descriptor) != NULL)
    {
        /* Fill the array with this new line */
        er_stringarray_append(file_buffer, line);
        /* Reinitialize buffer */
        memset(line, '\0', MAXSTRLENCONF);
    }

    fclose(file_descriptor);
    return file_buffer;
}

/**********************************************************************/
/**
 * @brief Actually parses the configuration into items, one per PRO CATG regexp 
 *
 * @param config_lines     The lines of text to convert
 * @param paf_configuration     The configuration
 * @returns 0 if successfull, !=0 otherwise
 *
 */
/**********************************************************************/
int er_paf_parse_config
(er_stringarray_t *     config_lines,
 struct er_paf_config * paf_configuration)
{
    unsigned int                iline = 0;
    unsigned int                nlines;
    unsigned int                nitem = 0;
    struct er_paf_config_item * items;

    cpl_msg_debug(__func__,"Parsing the paf configuration file");

    /* Alias for the items structure */
    items = paf_configuration->items;

    /* Parse all the lines */
    nlines = er_stringarray_size(config_lines);
    cpl_msg_debug(__func__,"Number of lines in pkd file: %d", nlines);
    while(iline < nlines)
    {
        char * current_line = er_stringarray_get(config_lines, (int)iline);
        if(strncmp(current_line, "PRO CATG",8) == 0)
        {
            struct er_paf_config_item* new_item;
            char *                     pro_catg_regexp = NULL;

            /* Allocate space for the new config item */
            ++nitem;
            items = cpl_realloc(items,
                                nitem*sizeof(struct er_paf_config_item));
            new_item = items + nitem -1;
            new_item->pro_catg_regexp = cpl_malloc(MAXSTRLENCONF*sizeof(char));
            new_item->requested_keywords = er_stringarray_new();
            paf_configuration->nitems = nitem;

            /* Get the regexp for PRO CATG */
            pro_catg_regexp = strchr(current_line, '=');
            if(pro_catg_regexp == NULL || strlen(pro_catg_regexp) < 2)
            {
                cpl_msg_error(er_func,"Parse error in PAF config. "
                              "No '=' delimiter");
                return -1;
            }
            strncpy(new_item->pro_catg_regexp,
                    pro_catg_regexp + 1, MAXSTRLENCONF);
            /* Remove leading and trailing characters */
            cx_strstrip(new_item->pro_catg_regexp);

            ++iline;
            while(iline < nlines)
            {
                char * keyword = er_stringarray_get(config_lines, (int)iline);
                if(strncmp(keyword, "PRO CATG",8) == 0)
                    break;
                /* Remove leading and trailing characters */
                cx_strstrip(keyword);

                /* If it is a keyword containing only spaces, do not count it */
                if(strcmp(keyword,"") != 0)
                {
                    /* Add this keyword to the requested keywords */
                    er_stringarray_append(new_item->requested_keywords, keyword);
                }
                ++iline;
            }
        }
        else
        {
            cpl_msg_error(er_func, "Error parsing PAF configuration file");
            return CPL_ERROR_BAD_FILE_FORMAT;
        }
    }
    cpl_msg_debug(__func__,"Number of pro catg definitions found in file: %d",
                  paf_configuration->nitems);

    /* Get the alias back */
    paf_configuration->items = items;

    return 0;
}

/**********************************************************************/
/**
 * @brief Deallocates the configuration to write PAF files
 *
 * @param paf_configuration     The configuration
 * @returns 0 if successfull, !=0 otherwise
 *
 * This function deallocates the configuration of a PAF writing process,
 * that is, the er_paf_config structure. 
 */
/**********************************************************************/

int er_paf_config_delete(struct er_paf_config* paf_configuration)
{
    unsigned int item;

    /* Deallocate the items */
    for (item = 0; item < paf_configuration->nitems; ++item)
    {
        er_stringarray_delete
        (paf_configuration->items[item].requested_keywords);
        cpl_free(paf_configuration->items[item].pro_catg_regexp);
    }
    cpl_free(paf_configuration->items);

    /* Deallocate the configuration itself */
    cpl_free(paf_configuration);

    return CPL_ERROR_NONE;
}


/**@}*/

/* End of file */
