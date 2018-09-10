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
#  include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <wordexp.h>

#include <cpl_dfs.h>
#include <cpl.h>

#include "ltdl.h"

#include "er_macros.h"
#include "er_help.h"
#include "er_stringarray.h"
#include "er_fileutils.h"
#include "er_paramutils.h"
#include "er_params.h"
#include "er_plugin.h"
#include "er_pluginlist.h"
#include "er_paf.h"
#include "er_json.h"
#include "er_stringutils.h"

static cpl_frameset *holdme;

#define SIZE_A  2048
#define SIZE_AA 2 * SIZE_A
#define SIZE_B  128
#define SIZE_C  1024
#define SIZE_D  1024

#define Xspace(x)  (x == ' ') || (x == '\t') || (x == '\n') 
#define Xnospace(x)  (x != ' ') && (x != '\t')  && (x != '\n')


/**
 * @defgroup esorex_plugin_process EsoRex Plugin Processing Functions
 *
 * EsoRex Plugin Processing Functions
 *
 */

/**@{*/



/**********************************************************************/
/**
 * @brief  Sets the message severity level for the terminal
 * @param  param_list  A list of the command line parameters
 * @param  flag  = 1, get level for logging
 *               = 2, get level for terminal messages
 * 
 * This function takes the list of all the command line parameters, and
 * checks for the existence of one to set the terminal message reporting
 * level. If it exists, the log-level is set to the requested value.
 *
 */
/**********************************************************************/

cpl_msg_severity message_severity (cpl_parameterlist *param_list, int flag)

{
    const char *msg_level_value = NULL;

    cpl_parameter *p;


    if (param_list == NULL) return (CPL_MSG_OFF);


    if (flag == 1)
        p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".log-level");
    else
        p = cpl_parameterlist_find (param_list, PACKAGE_RESOURCE ".msg-level");

    msg_level_value = cpl_parameter_get_string (p);
    ForceNull(msg_level_value)

    if (msg_level_value != NULL)
    {
        if (strcmp (msg_level_value, "debug") == 0)
            return CPL_MSG_DEBUG;

        else if (strcmp (msg_level_value, "info") == 0)
            return CPL_MSG_INFO;

        else if (strcmp (msg_level_value, "warning") == 0)
            return CPL_MSG_WARNING;

        else if (strcmp (msg_level_value, "error") == 0)
            return CPL_MSG_ERROR;

        else if (strcmp (msg_level_value, "off") == 0)
            return CPL_MSG_OFF;

        else
        {
            char  tmp[16];

            if (flag == 1)
                (void) strcpy(tmp,"logfile");
            else
                (void) strcpy(tmp,"terminal");
            cpl_msg_error
            (er_func, "%s messaging level '%s' is not recognized", tmp,msg_level_value);
        }
    }

    return (CPL_MSG_OFF);

}                               /* End of message_severity() */



int  add_size(const char *path, double *fsz)
{
    struct stat buf;


    if (stat(path,&buf) == -1)      /* get file structure info */
    {
        return(-1);
    }

    *fsz += (double) buf.st_size;       /* add current files size in bytes */

    return(0);
}

int  mysscanf (char *myline, char *path, char * tag,char * group)
{
    int   iout, keep;
    register int i;

    char  cc;
    char  *mpt[3];



    mpt[0] = path;
    mpt[1] = tag;
    mpt[2] = group;

    i = -1;
    keep = 0;

    for (iout=1; iout <4; iout++)
    {
        while(1)                /* find blank space */
        {
            cc = myline[++i];
            if (cc == '\0')
            {
                (void) strcpy(mpt[iout-1],&myline[keep]);
                return iout;
            }
            if (Xspace(cc)) break;
        }

        myline[i] = '\0';
        (void) strcpy(mpt[iout-1],&myline[keep]);

        if (iout == 3) return iout;     /* not more than 3 items... */

        while(1)                /* find non-space char */
        {
            cc = myline[++i];
            if (cc == '\0') return iout;
            if (Xnospace(cc)) break;
        }

        keep = i;
    }
    return iout;
}



/**********************************************************************/
/**  
 * @brief
 *   Create a new frame set from a @em set @em of @em frames file.
 * 
 * @param name  Input file path.
 * @param set   Frame set to be updated with the contents of @em name.
 * @param flag_check_sof_exist  Boolean indicating if an error will be produced
 *                              if the the Set-of-Frames is missing.
 * 
 * @return Pointer to the newly created frame set if @em set was @c NULL,
 *   or the updated set @em set. In case an error occurred the return value
 *   is @c NULL.
 *
 * The function reads the given input file @em filename and either, 
 * if @c NULL is passed for @em set, creates a new frame set from its
 * contents, or updates @em set, if @em set points to an already existing   
 * frame set.
 *
 * @note
 *   The current format of the @em set @em of @em frames file is as follows:
 *     - One frame description per line.
 *     - Each frame description consists of an absolute file path followed
 *       by the frame tag and, optionally, the group the frame belongs to.
 *     - Currently the only defined group tags are RAW, CALIB and PRODUCT.
 *     - The individual fields are separated by whitespace.
 *     - Blank lines are ignored.
 *     - Lines that begin with a hash (#) are treated as comment lines.
 */
/**********************************************************************/

cpl_frameset *er_frameset_load
(const char *name, cpl_frameset * set, int flag_check_sof_exist)
{
    FILE *fp;

    char  *line, *path, *group, *tag, *xpath;
    char  estr[240];
    register char  cc;

    int line_number = 0, created = 0;
    int  ii, mlen, n, no_comnt;

    int load_sucesfull = 1;

    cpl_frame *frame = NULL;

    cpl_msg_severity msg_level;


    if (name == NULL)
    {
        cpl_msg_error (er_func, "No SOF name provided, when one was expected.");
        return NULL;
    }

    if (!(fp = fopen (name, "r")))
    {
        cpl_msg_error (er_func, "Unable to open SOF file '%s'", name);
        return NULL;
    }

    /* If set was NULL, create a new frame set that we can fill. */

    if (!set)
    {
        set = cpl_frameset_new ();
        created = 1;
    }


    line = (char *) cpl_calloc((size_t) SIZE_C, (size_t) 1);
    path = (char *) cpl_calloc((size_t) SIZE_C, (size_t) 1);
    xpath = (char *) cpl_calloc((size_t) SIZE_C, (size_t) 1);
    group = (char *) cpl_calloc((size_t) SIZE_C, (size_t) 1);
    tag = (char *) cpl_calloc((size_t) SIZE_C, (size_t) 1);

    msg_level = cpl_msg_get_level();          /* get the message level */


    /* Loop over all the lines in the set-of-frames file */

    while( fgets (line, (SIZE_C-1), fp))
    {
        line_number++;
        no_comnt = 0;       /* check for comments + empty records */
        line[SIZE_C-1] = '\0';
        mlen = (int)strlen(line);

        /* printf("line no. %d has length %d\n",line_number,mlen);  */

        for (ii=0; ii<mlen; ii++)
        {
            cc = line[ii];
            if (Xnospace(cc))           /* look for 1st non-blank char */
            {
                if (cc != '#') no_comnt = 1;    /* # indicates comment line */
                break;
            }
        }

        if (no_comnt == 1)
        {
            n = mysscanf (&line[ii], xpath, tag, group);      /* split up line */

            if (n >= 1)
            {
                if (xpath[0] == '$')        /* check for env. variables */
                {
                    int   wer;
                    char  *wp;
                    wordexp_t  wresult;

                    wer = wordexp(xpath,&wresult,0);    /* expand the $VAR/file */
                    if (wer == 0)
                    {
                        wp = *wresult.we_wordv;             /* use 1st and only result word */
                        if (wp != NULL)
                            (void) strcpy(path,wp);
                        else
                            wer = 1;
                        wordfree(&wresult);
                    }

                    if (wer != 0)
                    {
                        (void) snprintf(estr,(size_t)238,"Expansion of %s failed...",xpath);
                        cpl_msg_error (er_func, "%s", estr);
                        goto dealloc;
                    }
                }
                else
                {
                    (void) strcpy(path,xpath);
                }

                if (!(fileutils_file_exists (path)))    /* Ensure that the path exists */
                {                   /* Invalid file name? */

                    /* allow a single FITS file instead of SOF */

                    if (strncmp(path,"SIMPLE",6) == 0)
                    {               /* it's a single FITS file */
                        n = (int) strlen(name);
                        if (n > (SIZE_C-1))
                        {
                            n = SIZE_C - 1;
                            (void) strncpy(path,name,(size_t)n);
                            path[n] = '\0';
                            cpl_msg_warning (er_func, "FITS file name truncated to %d chars", n);
                        }
                        else
                        {
                            (void) strcpy(path,name);
                        }
                        (void) strcpy(tag,"COMMAND_LINE");
                        n = -1;             /* to remember that later on */
                    }
                    else                /* Yes, bad filename */
                    {
                        char * tagmsg = n >= 2
                                ? cpl_sprintf("tag '%s'", tag)
                                        : cpl_sprintf("no tag");
                        char * groupmsg = n >= 3
                                ? cpl_sprintf("group '%s'", group)
                                        : cpl_sprintf("no group");

                        if (flag_check_sof_exist)
                        {
                            cpl_msg_error (er_func,
                                           "Could not open the input file '%s' with %s and %s in "
                                           "line %d of the SOF '%s'", path, tagmsg, groupmsg,
                                           line_number, name);
                            load_sucesfull = 0;
                        }
                        else
                        {
                            cpl_msg_debug (er_func,
                                           "Could not open the input file '%s' with %s and %s in "
                                           "line %d of the SOF '%s'", path, tagmsg, groupmsg,
                                           line_number, name);
                        }
                        cpl_free(tagmsg);
                        cpl_free(groupmsg);

                    }
                }

                frame = cpl_frame_new ();           /* allocate a new frame */
                cpl_frame_set_filename (frame, path);       /* and provide the filename */

                if (n == 1)             /* no tag in SOF line */
                    cpl_frame_set_tag (frame, "");
                else
                    cpl_frame_set_tag (frame, tag);

                /* Set the group component of the frame (or set a default) */
                if (n > 2)              /* so, omitted for single FITS file */
                {
                    cpl_frame_group grp;

                    if (!strcmp (group, CPL_FRAME_GROUP_RAW_ID))
                        grp = CPL_FRAME_GROUP_RAW;
                    else
                    {
                        if (!strcmp (group, CPL_FRAME_GROUP_CALIB_ID))
                            grp = CPL_FRAME_GROUP_CALIB;
                        else
                        {
                            if (!strcmp (group, CPL_FRAME_GROUP_PRODUCT_ID))
                                grp = CPL_FRAME_GROUP_PRODUCT;
                            else
                                grp = CPL_FRAME_GROUP_NONE;
                        }
                    }

                    cpl_frame_set_group (frame, grp);
                }

                cpl_frameset_insert (set, frame);
                if (n == -1) goto dealloc;          /* was single FITS file */
            }

            else
            {           /* Invalid frame description. */
                if (msg_level == CPL_MSG_DEBUG)
                {
                    cpl_msg_debug (er_func,
                                   "Invalid frame description '%s', specified on line %d of the SOF "
                                   "file, '%s'.", path, line_number, name);
                }
            }
        }

    }                            /* End of while() */


    dealloc:                /* Deallocate resources */
    cpl_free (line);
    cpl_free (path);
    cpl_free (xpath);
    cpl_free (group);
    cpl_free (tag);
    fclose (fp);

    /* Check if the loading was successfully */
    if(!load_sucesfull)
    {
        if(created)
            cpl_frameset_delete(set);
        return NULL;
    }

    /* Return the pointer to the frameset */
    return set;

}                               /* End of er_frameset_load() */



/**********************************************************************/
/**
 * @brief Process the Plugin
 *
 * @param frameset         Set of frames to be moved
 * @param output_prefix    Character string with the prefix to be used
 *
 * @returns 0 
 *
 * This function updates the DATAMD5 keyword for any product already created
 * in the case a recipe returned an error status (failed)
 *
 */
/**********************************************************************/

inline static int
_upda_products (cpl_frameset *frameset, unsigned int flags)
{
    cpl_error_code reto;

    cpl_size pcount = 0;
    cpl_frameset_iterator *it = cpl_frameset_iterator_new(frameset);

    const cpl_frame *frame = NULL;

    while ((frame = cpl_frameset_iterator_get_const(it)) != NULL) {

        if (cpl_frame_get_group(frame) == CPL_FRAME_GROUP_PRODUCT) {
            ++pcount;
        }

        cpl_frameset_iterator_advance(it, 1);

    }

    cpl_frameset_iterator_delete(it);
    it = NULL;

    if ((pcount != 0) && (flags != CPL_DFS_SIGNATURE_NONE)) {
        cpl_msg_info(er_func, "Calculating product checksums");
    }

    reto = cpl_dfs_sign_products(frameset, flags);

    if (reto != CPL_ERROR_NONE)
    {
        cpl_msg_error (er_func, "could not update the product header...");
    }

    return 0;
}




/**********************************************************************/
/**
 * @brief Process the Plugin
 *
 * @param plugin_name       Name of Plugin to process
 * @param frameset         Set of frames to be moved
 * @param output_directory Character string with the output directory
 * @param output_prefix    Character string with the prefix to be used
 * @param flag_noprefix    Boolean indicating if prefix/number system 
 *                         should be ignored, in which case the recipe
 *                         output filenames are used.
 * @param flag_readonly    Boolean indicating whether or not output files
 *                         should be set to oct(444) for Paranal.
 *
 * @returns 0 if successfull, !=0 otherwise
 *
 * This function moves the products from the temporary working
 * directory, as used by the recipe, to the output directory requested 
 * by the user.
 *
 * This function also sets the PIPEFILE keyword in the FITS header if
 * appropriate. If this cannot be done, or the FITS card has 
 * insufficient size to then a warning is printed.
 */
/**********************************************************************/

inline static int
_move_products (cpl_frameset * frameset,
                char * output_directory,
                char * output_prefix,
                char * link_directory,
                unsigned int signature_mask,
                int flag_noprefix, int flag_nolink, int flag_readonly)

{
    cpl_frame *frame = NULL;       /* The current frame being processed */

    int filenum = 0;               /* Number of the output product */
    int ii, move_err;
    int target_path_is_cwd = 0;    /* assume that output path and CWD differ */

    const char *fileutils_pntr;
    const char *cpp;
    const char *input_name = NULL;       /* Name of the input file */

    char *output_name = NULL;            /* Name of the output file */
    char *suffix = NULL;                 /* The file suffix */
    char current_dir[FILEMAX];           /* Name of current (work) directory */
    char  command_str[SIZE_AA];
    char link_dir[PATHSET_MAX];          /* Name of the symbolic link dir. */
    char output_dir[PATHSET_MAX];        /* Name of the output dir. */

    struct stat sb_cdir;
    struct stat sb_tdir;

    cpl_msg_severity msg_level;


    /* allocate space */

    output_name = (char *) cpl_malloc((size_t) SIZE_A);

    current_dir[0] = '\0';


    /* Start by getting the "cleaned" output directory name */

    fileutils_pntr =
            er_fileutils_dot_replace (er_fileutils_tilde_replace (output_directory));
    if (fileutils_pntr == NULL)
    {
        cpl_msg_warning (er_func, "missing output directory ...");
    }
    (void) strcpy (output_dir, fileutils_pntr);

    fileutils_pntr = er_fileutils_tilde_replace (link_directory);
    if (fileutils_pntr == NULL)
    {
        cpl_msg_warning (er_func, "missing directory for symbolic link ignored...");
    }
    (void) strcpy (link_dir, fileutils_pntr);


    /* Determine the name of the current directory */

    if (getcwd (current_dir, (size_t) FILEMAX) == NULL)
    {
        current_dir[0] = '.';
        current_dir[1] = '\0';
    }

    // FIXME: This is a hot fix to avoid a failing move of the product
    //        files in case the output directory is identical to the
    //        current directory but referenced through a link so that
    //        the two path names are different.
    //        Finally this section should use the file utilities function,
    //        which need to be reviewed first.

    if ((stat(current_dir, &sb_cdir) != 0) || (stat(output_dir, &sb_tdir) != 0))
    {
        cpl_msg_error(er_func, "cannot get directory status information");
        move_err = 65;
        goto report_products;
    }
    else
    {

        /*
         * Check whether the current working directory and the final
         * location of the product files are identical, i.e. they are
         * located on the same device and they have the the same inode.
         */

        if ((sb_cdir.st_dev == sb_tdir.st_dev) &&
                (sb_cdir.st_ino == sb_tdir.st_ino))
        {
            target_path_is_cwd = 1;
        }

    }

    /* MD5 errors ignored at the moment... */
    ii = _upda_products(frameset, signature_mask);
    holdme = frameset;


    /* Loop through the frameset */

    move_err = 0;               /* init possible error from moving files */

    cpl_frameset_iterator *it = cpl_frameset_iterator_new(frameset);

    while ((frame = cpl_frameset_iterator_get(it)) != NULL)
    {
        if (cpl_frame_get_group (frame) == CPL_FRAME_GROUP_PRODUCT)
        {

            /*
             *  Get the filename as specified by the frame
             */

            input_name = cpl_frame_get_filename (frame);

            /*
             * Temporarily remove the suffix (find the first '.' and
             * replace it)
             */

            suffix = strrchr (input_name, '.');
            if (suffix) *suffix = '\0';

            /*
             *  Generate the output filename (keep original if this
             * is suppressed)
             */

            if (flag_noprefix != 0)
            {
                (void) strcpy (output_name, input_name);
            }
            else
            {
                (void) sprintf (output_name, "%s_%04d", output_prefix,
                                filenum);
            }

            /*
             * Now restore the suffix (note difference between
             *  suffix and *suffix)
             */

            if (suffix) *suffix = '.';
            if (suffix) (void) strcat (output_name, suffix);

            /*
             *  Check if any move-command is required
             */

            cpp = cpl_frame_get_filename (frame);
            if (cpp == NULL)
            {
                cpl_frameset_iterator_delete(it);
                cpl_msg_error (er_func, "cpl_frame_get_filename failed");
                move_err = 63;
                goto report_products;
            }


            /*
             * Move product files only if the current directory and product
             * target path are different, or if the local product name and
             * the final product name are different.
             */

            if (target_path_is_cwd && (strcmp(cpp, output_name) == 0))
            {

                /*
                 *  No move. Report that the product has been created
                 */

                cpl_msg_info (er_func, "Created product %s (in place)",
                              output_name);

            }
            else
            {

                /*
                 *  Create the full path
                 */

                char  tmp[SIZE_A];

                if (output_name[0] == '\0')
                    (void) strcpy(tmp,"dummy_output");
                else
                    (void) strcpy(tmp,output_name);

                (void) strcpy (output_name, output_dir);
                ii = (int) strlen(output_name) - 1;
                if (output_name[ii] != '/') output_name[++ii] = '/';
                (void) strcpy (&output_name[ii+1],tmp);

                strcpy(command_str,"mv ");
                strcat(command_str,cpp);
                strcat(command_str," ");
                strcat(command_str,output_name);
                msg_level = cpl_msg_get_level();       /* if in in debug mode */
                if (msg_level == CPL_MSG_DEBUG)
                {                           /* print the command string */
                    printf("move command:\n%s\n",command_str);
                }
                ii = system(command_str);
                if (ii != 0)
                {                   /* now, product file should be there... */
                    cpl_frameset_iterator_delete(it);
                    cpl_msg_error (er_func,
                                   "Unable to move product file to final output path (%s)",
                                   output_name);
                    move_err = 61;
                    goto report_products;
                }

                /* Update the name in the sof */
                cpl_frame_set_filename(frame, output_name);
                cpl_msg_info (er_func, "Created product %s", output_name);
            }

            /* If required, generate the symbolic links for the archive system */
            if (flag_nolink == 0)
            {
                if (er_fileutils_link (link_dir, output_name) != EXIT_SUCCESS)
                {
                    cpl_frameset_iterator_delete(it);

                    cpl_msg_error (er_func,
                                   "Unable to create symbolic link for %s in "
                                   "directory %s", output_name, link_dir);
                    move_err = 62;
                    goto report_products;
                }
            }

            /* If required, lock the file permissions to read-only */
            if (flag_readonly != 0) chmod (output_name, 0444);

            filenum++;      /* Increment the product number */
        }

        /* Get the next frame in the list */

        cpl_frameset_iterator_advance(it, 1);

    }               /* End of (while) loop through all frames in frameset */

    cpl_frameset_iterator_delete(it);

    /* Report how many products were created (getting the plural right!) */

    report_products:
    cpl_msg_info (er_func, "%d product%s created", filenum,
                  (filenum == 1) ? "" : "s");
    cpl_free (output_name);             /* free the allocated memory */

    return move_err;
}                            


/**********************************************************************/
/**
 * @brief Process the Plugin
 *
 * @param caller_parameters Caller parameters
 * @param plugin_name       Name of Plugin to process
 * @param sof_filename_list List of strings with SOF filenames
 * @param argc              Count of remaining arguments from cmdl
 * @param argv              Handle to remaining arguments from cmdl
 *
 * @returns 0 if successfull, !=0 otherwise
 *
 * Processes a Plugin.
 * Write more here...
 *
 */
/**********************************************************************/

int plugin_process_plugin (cpl_parameterlist *caller_parameters,
                           char *plugin_name,
                           er_stringarray_t * sof_filename_list,
                           int argc, char * argv[])

{
    const char *val_output_mask,
    *val_link_dir,
    *val_log_dir,
    *val_log_file,
    *val_output_dir, *cpp;
    const char *val_string_tmp = NULL;
    const char *cdescr;

    const char *val_paf_config;

    const char *val_products_sof;

    const char *val_create = NULL;

    char **val_recipe_dirs = NULL,
            *log_file_name_full = NULL;

    char *cptr,
    *plugin_conf_file_global=NULL,
    *plugin_conf_file_local=NULL,
    library_path[MAXSTRLENCONF],
    cmdline[MAXSTRLENCONF];

    char *help_str = NULL;

    int e_code = 0, e_code2 = 0, i = 0, ii = 0;
    int countr, n, sz;

    int flag_time_plugin = 0;
    int flag_plugin_found = 0;
    int flag_help = 0;
    int flag_man_page = 0;
    int flag_show_all = 0;
    int flag_create_config = 0;
    int flag_params = 0;
    int flag_readonly = 0;
    int flag_noprefix = 0;
    int flag_nolink = 0;
    int f_val_unload_plugin = 1;        /* default is TRUE */


    unsigned int signature_mask = CPL_DFS_SIGNATURE_DATAMD5 | CPL_DFS_SIGNATURE_CHECKSUM;

    unsigned long  uil;

    double plugin_time_start=0.0, plugin_time_end=0.0;

    cpl_msg_severity msg_level = CPL_MSG_ERROR,
            msg_sev_logfile, msg_sev_terminal;
    cpl_parameter *p = NULL;
    cpl_plugin *tplugin = NULL;
    lt_dlhandle module = NULL;
    cpl_plugin_func plugin_func_init = NULL,
            plugin_func_exec = NULL,
            plugin_func_deinit = NULL;
    cpl_recipe *trecipe = NULL;
    cpl_recipe2 *t2recipe = NULL;

    cpl_error_code log_err;

    er_stringarray_t *list_of_pllib_names = NULL;


    library_path[0] = '\0';

    /* check, if we want to take times */
    p = cpl_parameterlist_find (caller_parameters, PACKAGE_RESOURCE ".time");
    flag_time_plugin = cpl_parameter_get_bool (p);
    if (flag_time_plugin != 0)          /* --time option */
    {
        plugin_time_start = cpl_test_get_walltime();
    }

    /*
     * check whether hidden parameters should be visible in the help and
     * man-page output
     */

    p = cpl_parameterlist_find(caller_parameters, PACKAGE_RESOURCE ".show-hidden");
    flag_show_all = cpl_parameter_get_bool(p);


    msg_sev_logfile = message_severity (caller_parameters,1);
    msg_sev_terminal = message_severity (caller_parameters,2);

    p = cpl_parameterlist_find (caller_parameters, PACKAGE_RESOURCE ".log-dir");
    val_log_dir = cpl_parameter_get_string (p);
    ForceNull(val_log_dir)

    p = cpl_parameterlist_find (caller_parameters, PACKAGE_RESOURCE ".log-file");
    val_log_file = cpl_parameter_get_string (p);
    ForceNull(val_log_file)
    log_file_name_full =
            fileutils_create_fqfname ((char *) val_log_dir, (char *) val_log_file);


    log_err = cpl_msg_set_log_level (msg_sev_logfile);
    if (log_err != CPL_ERROR_NONE)
    {
        (void) printf ("WARNING : EsoRex is unable to establish the message log "
                "(Error = %d)\n", log_err);
        (void) printf ("          (Check write permission for temporary files.)\n");
    }

    /*
     * Always save the library dependencies and the invocation command line to
     * the logfile (only)
     */

    cpl_msg_set_level (CPL_MSG_OFF);
    cdescr = cpl_get_description(CPL_DESCRIPTION_DEFAULT);

    (void) cpl_msg_info(er_func,"This is EsoRex, version %s", PACKAGE_VERSION);
    (void) cpl_msg_info(er_func,"using the libraries: %s",cdescr);

    sz = 0;
    for (i=0; i < argc - 1; ++i) {
        int nc = snprintf(&cmdline[sz], MAXSTRLENCONF, "%s ", argv[i]);
        sz += nc;
    }
    sz = snprintf(&cmdline[sz], MAXSTRLENCONF, "%s", argv[argc - 1]);
    (void) cpl_msg_info(er_func,"Invocation command line was: %s", cmdline);
    cpl_msg_set_level (msg_sev_terminal);


    /* look for recipe directory */

    p = cpl_parameterlist_find (caller_parameters, PACKAGE_RESOURCE ".recipe-dir");
    val_string_tmp = cpl_parameter_get_string (p);
    ForceNull(val_string_tmp)
    if (val_string_tmp != NULL)
        val_recipe_dirs = cx_strsplit (val_string_tmp, ":", -1);

    /* get the plugin  */

    if ((plugin_name[0] == ' ') || (plugin_name[0] == '\0'))
        sz = 0;
    else
        sz = 1;

    if ((val_recipe_dirs != NULL) && (sz > 0))
    {
        list_of_pllib_names = er_pluginlist_create_list (val_recipe_dirs);

        sz = er_pluginlist_get_libpath (list_of_pllib_names, plugin_name, library_path);
        if (sz > 0)
        {
            tplugin = er_pluginlist_get_plugin (library_path, plugin_name, &module);
            if (tplugin != NULL)
            {
                flag_plugin_found = 1;
                goto next_step;
            }
        }
        cpl_msg_error (er_func, "Unable to find recipe '%s'."
                       " Check that the recipe is in the path specified by the"
                       " '--recipe-dir' option.", plugin_name);
        e_code = CPL_ERROR_INCOMPATIBLE_INPUT;
        goto cleanup;               /* we cannot do anything */
    }


    /* ------------------------------------- */
    /* here we process the plugin (if found) */
    /* ------------------------------------- */

    next_step:
    uil = cpl_plugin_get_type(tplugin);

    /*  initialize the plugin/recipe structure: trecipe->parameters is handled by plugin  */

    if (uil < 2)
    {
        n = sizeof (cpl_recipe);
        trecipe = (cpl_recipe *) cpl_calloc (1,(size_t) n);
    }
    else
    {
        n = sizeof (cpl_recipe2);
        t2recipe = (cpl_recipe2 *) cpl_calloc (1,(size_t) n);
        trecipe = &t2recipe->base;
    }

    trecipe->frames = cpl_frameset_new ();


    if (flag_plugin_found != 0)         /*  = TRUE */
    {
        cpl_plugin_copy ((cpl_plugin *) & trecipe->interface, tplugin);


        /*  Run Plugin Initialization...  */

        plugin_func_init = cpl_plugin_get_init ((cpl_plugin *) trecipe);
        if (plugin_func_init != NULL)
        {
            cpl_msg_set_domain (cpl_plugin_get_name ((cpl_plugin *) trecipe));
            e_code = plugin_func_init ((cpl_plugin *) trecipe);
            cpl_msg_set_domain (PACKAGE);
            if (e_code != 0)
            {
                cpl_msg_error (er_func, "Init of recipe failed...");
                goto plugin_deinit;
            }
        }

        /* loop through all parameters in the list */

        p = cpl_parameterlist_get_first (trecipe->parameters);
        while (p != NULL)
        {                           /* Set the tag */
            char recip_def[] = "recipe default";
            char  *myptr;

            myptr = recip_def;
            er_manage_sources(1,cpl_parameter_get_name(p),&myptr);

            /* Get the next parameter in the list */
            p = cpl_parameterlist_get_next (trecipe->parameters);
        }


        cpp = cpl_plugin_get_name ((cpl_plugin *) trecipe);


        /* Parse any global plugin configuration file */

        cptr = getenv ("HOME");
        if (cptr == NULL)
            n = 0;
        else
            n = (int) strlen(cptr);
        if (cpp != NULL) n += (int)strlen(cpp);

        n = n + ((int)strlen(cpp)) +
                ((int) strlen(GLOBAL_RC_DIR)) + ((int) strlen(GLOBAL_RC_EXTENSION));
        n += 4;                 /* for the additionial const chars */

        plugin_conf_file_global = (char *) cpl_malloc((size_t) n);
        if ( plugin_conf_file_global == NULL)
        {
            cpl_msg_error (er_func,
                           "Could not allocate %d bytes for plugin_conf_file_global",n);
            e_code = CPL_ERROR_ILLEGAL_OUTPUT;
            goto plugin_deinit;
        }


        if (cptr == NULL)
            (void) strcpy(plugin_conf_file_global, "/");
        else
        {
            (void) strcpy(plugin_conf_file_global,cptr);
            (void) strcat(plugin_conf_file_global, "/");
        }
        (void) strcat(plugin_conf_file_global, GLOBAL_RC_DIR);
        (void) strcat(plugin_conf_file_global, "/");
        if (cpp != NULL) (void) strcat(plugin_conf_file_global, cpp);
        (void) strcat(plugin_conf_file_global, GLOBAL_RC_EXTENSION);

        if (e_code == 0)
        {
            if (fileutils_file_exists(plugin_conf_file_global) != 0)
            {
                e_code = params_parse_config_file (trecipe->parameters,
                                                   plugin_conf_file_global);
                if (e_code != 0) goto plugin_deinit;
            }
        }

        /* Parse any specifically specified plugin configuration file */
        if (e_code == 0)
        {
            p = cpl_parameterlist_find (caller_parameters,
                                        PACKAGE_RESOURCE ".recipe-config");
            val_string_tmp = cpl_parameter_get_string (p);

            if ((val_string_tmp != (char *) 0) && (strlen (val_string_tmp) > 0))
                e_code = params_parse_config_file (trecipe->parameters, val_string_tmp);

        }

        if (e_code == 0)
        {
            char *str_tmp = NULL;

            e_code = params_parse_config_commandline
                    (trecipe->parameters, plugin_name,
                     sof_filename_list, argc, argv, 0);
            if (e_code != 0) goto plugin_deinit;

            params_parse_config_postprocess (trecipe->parameters);
            ii = er_stringarray_size(sof_filename_list);

            for (i=0; i<ii; i++)
            {
                int flag_check_sof_exist;
                /* Get the check input files parameter */
                p = cpl_parameterlist_find (caller_parameters,
                                            PACKAGE_RESOURCE ".check-sof-exist");
                flag_check_sof_exist = cpl_parameter_get_bool (p);

                str_tmp = er_stringarray_get (sof_filename_list, i);
                if (er_frameset_load (str_tmp, trecipe->frames, flag_check_sof_exist) == NULL)
                {
                    cpl_msg_error (er_func, "Problem occurred loading frameset "
                                   "from the SOF file '%s'.\n"
                                   "If you want to ignore errors of missing files in "
                                   "the sof, set '--check-sof-exist=false'.", str_tmp);

                    e_code = CPL_ERROR_FILE_NOT_FOUND;
                }
            }
        }

        /* Check if we want to create the plugin (recipe) configuration file */

        p = cpl_parameterlist_find (caller_parameters, PACKAGE_RESOURCE ".create-config");
        if (cpl_parameter_get_default_flag (p))
        {
            val_create = cpl_parameter_get_string (p);
            ForceNull(val_create)
            if (val_create != NULL)                      /* --create-config=bla.bla  found */
            {                 /* handle FALSE and TRUE from boolean history of this param */
                if (strcmp(val_create,"TRUE") == 0)
                {
                    flag_create_config = 1;
                }
                else if (strcmp(val_create,"FALSE") != 0)
                {
                    flag_create_config = 11;             /* indicates filename is given */
                }
            }
            else                                 /* --create-config  is interpreted as ...=TRUE */
            {
                flag_create_config = 1;
            }
        }
        if (flag_create_config == 11)
        {
            er_help_create_config (flag_create_config, plugin_name, val_create,
                                   caller_parameters, trecipe->parameters);
            e_code = -99999;
        }
        else if (flag_create_config == 1)
        {
            er_help_create_config (flag_create_config, plugin_name, NULL,
                                   caller_parameters, trecipe->parameters);
            e_code = -99999;
        }

        /* Check if we want to display the plugin help */

        p = cpl_parameterlist_find (caller_parameters, PACKAGE_RESOURCE ".man-page");
        flag_man_page = cpl_parameter_get_bool (p);
        if (flag_man_page != 0)
        {
            er_help_manpage(trecipe, flag_show_all);
            e_code = -99999;
        }

        /* Check if we want to display the plugin help */

        p = cpl_parameterlist_find (caller_parameters, PACKAGE_RESOURCE ".help");
        flag_help = cpl_parameter_get_bool (p);
        if (flag_help != 0)
        {
            int  h_size = 512;
            char tmp_buf[512];

            help_str = (char *) cpl_malloc((size_t) h_size);
            if (help_str == NULL)
            {
                cpl_msg_error (er_func, "Could not allocate %d bytes for help_str",h_size);
                e_code = CPL_ERROR_ILLEGAL_OUTPUT;
                goto plugin_deinit;
            }


            (void) strcpy (help_str, "Recipe: ");
            countr = 9;
            cpp = cpl_plugin_get_name ((cpl_plugin *) trecipe);
            n = (int) strlen(cpp);
            countr += (5 + n);
            if (countr > h_size)
            {                   /* allocate more space */
                h_size = countr + 120;
                er_enlarge(er_func,&help_str,h_size);
            }
            (void) strcat (help_str, cpp);
            (void) strcat (help_str, " -- ");

            cpp = cpl_plugin_get_synopsis ((cpl_plugin *) trecipe);
            n = (int) strlen(cpp) + countr;
            if (n > h_size)
            {                   /* allocate more space */
                h_size = n + 120;
                er_enlarge(er_func,&help_str,h_size);
            }
            (void) strcat (help_str, cpp);
            printf ("%s\n\n", er_strutils_split (help_str, 2, er_strutils_termwidth ()));

            msg_level = cpl_msg_get_level();        /* test, if we are in debug mode */
            if (msg_level == CPL_MSG_DEBUG)
            {
                (void) strcpy (help_str, "Library: ");
                countr = 12;                /* length of above */
                n = countr + (int) strlen(library_path);
                if (n > h_size)
                {                       /* allocate more space */
                    h_size = n + 120;
                    er_enlarge(er_func,&help_str,h_size);
                }

                (void) strcat (help_str, library_path);
                printf ("%s\n\n", er_strutils_split (help_str, 2, er_strutils_termwidth ()));
            }

            /* Display the actual help for the plugin */

            er_help_display(plugin_name, trecipe->parameters, flag_show_all);

            /* Explain why the esorex help doesn't appear */

            (void) strcpy (tmp_buf,"For help on the options of " PACKAGE
                           " itself, please use the command '" PACKAGE " --help' "
                           "(that is, without specifying any recipe name). "
                           "For more information about the recipe, one can also use "
                           "the command '" PACKAGE " --man-page ");
            countr = (int) strlen(tmp_buf);

            cpp = cpl_plugin_get_name ((cpl_plugin *) trecipe);
            n = (int) strlen(cpp) + countr + 4;

            if (n > h_size)
            {
                h_size = n;
                er_enlarge("plugin_process_plugin",&help_str,h_size);
            }

            (void) strcpy (help_str, tmp_buf);
            (void) strcat (help_str, cpp);
            (void) strcat (help_str, "'.");
            (void) printf
                    ("%s\n", er_strutils_split (help_str, 0, er_strutils_termwidth ()));

            goto plugin_deinit;         /* avoid plugin execution */
        }

        /* Check if we want to display the plugin parameters */

        p = cpl_parameterlist_find (caller_parameters,PACKAGE_RESOURCE ".params");
        flag_params = cpl_parameter_get_bool (p);
        if (flag_params != 0)
        {
            if (er_paramutils_print_list(trecipe->parameters,"Recipe Parameters") != 0)
            {
                cpl_msg_error (er_func,"Unable to print the recipe parameter list\n");
                e_code = CPL_ERROR_INCOMPATIBLE_INPUT;
            }
            else
                goto plugin_deinit;         /* avoid plugin execution */
        }

        /* Check if we don't want to cleanup after processing (for debugging) */

        p = cpl_parameterlist_find(caller_parameters, PACKAGE_RESOURCE ".unload-plugin");
        f_val_unload_plugin = cpl_parameter_get_bool(p);


        /*
         *  Set product signature bit mask
         */

        p = cpl_parameterlist_find(caller_parameters, PACKAGE_RESOURCE ".no-checksum");

        if (cpl_parameter_get_bool(p) == TRUE) {
            signature_mask &= ~CPL_DFS_SIGNATURE_CHECKSUM;
        }

        p = cpl_parameterlist_find(caller_parameters, PACKAGE_RESOURCE ".no-datamd5");

        if (cpl_parameter_get_bool(p) == TRUE) {
            signature_mask &= ~CPL_DFS_SIGNATURE_DATAMD5;
        }


        /*
         *  Run Plugin Execute - if all went well until here...
         */

        if (e_code == 0)
        {
            plugin_func_exec = cpl_plugin_get_exec ((cpl_plugin *) trecipe);

            if (plugin_func_exec != NULL)
            {           /* We have a pointer, so run the plugin */
                const char *recipe_name;

                recipe_name = cpl_plugin_get_name ((cpl_plugin *) trecipe);
                cpl_msg_set_domain (recipe_name);
                e_code = plugin_func_exec ((cpl_plugin *) trecipe);
                cpl_msg_set_domain (PACKAGE);
                if (e_code != 0)
                {
                    cpl_msg_error (er_func,"Execution of recipe '%s' failed, status = %d",
                                   recipe_name,e_code);
                }
            }
            else
            {
                /* NULL-pointer, so we simply set an error */
                e_code = CPL_ERROR_INCOMPATIBLE_INPUT;
            }
        }

        if (e_code == 0)            /* successful execution of recipe */
        {
            p = cpl_parameterlist_find (caller_parameters,
                                        PACKAGE_RESOURCE ".output-dir");
            val_output_dir = cpl_parameter_get_string (p);
            ForceNull(val_output_dir)

            p = cpl_parameterlist_find (caller_parameters,
                                        PACKAGE_RESOURCE ".link-dir");
            val_link_dir = cpl_parameter_get_string (p);
            ForceNull(val_link_dir)

            p = cpl_parameterlist_find (caller_parameters,
                                        PACKAGE_RESOURCE ".output-prefix");
            val_output_mask = cpl_parameter_get_string (p);
            ForceNull(val_output_mask)

            p = cpl_parameterlist_find (caller_parameters,
                                        PACKAGE_RESOURCE ".output-readonly");
            flag_readonly = cpl_parameter_get_bool (p);

            p = cpl_parameterlist_find (caller_parameters,
                                        PACKAGE_RESOURCE ".suppress-prefix");
            flag_noprefix = cpl_parameter_get_bool (p);

            p = cpl_parameterlist_find (caller_parameters,
                                        PACKAGE_RESOURCE ".suppress-link");
            flag_nolink = cpl_parameter_get_bool (p);

            /* Move all product files to output directory */

            e_code = _move_products(trecipe->frames, (char *)val_output_dir,
                                    (char *)val_output_mask, (char *)val_link_dir,
                                    signature_mask,
                                    flag_noprefix, flag_nolink, flag_readonly);

            if (e_code != 0) cpl_msg_error (er_func,
                                            "An error occurred while trying to move the output products");
        }

        else            /* recipe failed, only get correct MD5 sum */
        {           /* for all product files created anyway */
            if (e_code != -99999)
            {
                ii = _upda_products (trecipe->frames, signature_mask);
            }
        }

        /* Run PAF creation */
        p = cpl_parameterlist_find (caller_parameters,
                                    PACKAGE_RESOURCE ".paf-config");
        val_paf_config = cpl_parameter_get_string (p);
        ForceNull(val_paf_config)
        if(val_paf_config != NULL)
        {
            if (e_code == 0)    /* successful execution of recipe and product move*/
            {
                e_code = er_create_recipe_pafs(trecipe->frames,
                                               trecipe->interface.name,
                                               val_paf_config);
                if (e_code != 0)
                {
                    cpl_msg_error (er_func,"Cannot create paf files, status = %d",
                                   e_code);
                }
            }
            else
            {
                cpl_msg_warning (er_func,"Writing of paf files omitted"
                                 " due to previous errors");
            }
        }

        /* Write output sof */
        p = cpl_parameterlist_find (caller_parameters,
                                    PACKAGE_RESOURCE ".products-sof");
        val_products_sof = cpl_parameter_get_string (p);
        ForceNull(val_products_sof)
        if(val_products_sof != NULL)
        {
            if (e_code == 0)    /* successful execution of recipe and product move*/
            {
                cpl_frameset * output_frames;
                cpl_frame    * frame;
                cpl_frameset_iterator *it =
                        cpl_frameset_iterator_new(trecipe->frames);

                output_frames = cpl_frameset_new();

                while ((frame = cpl_frameset_iterator_get(it)) != NULL)
                {
                    if (cpl_frame_get_group (frame) == CPL_FRAME_GROUP_PRODUCT)
                    {
                        cpl_frameset_insert(output_frames,
                                            cpl_frame_duplicate(frame));
                    }
                    cpl_frameset_iterator_advance(it, 1);
                }

                cpl_frameset_iterator_delete(it);

                if(strlen(val_products_sof)>5 &&
                        strncmp(val_products_sof + strlen(val_products_sof) - 5,
                                ".json", 5) == 0 )
                    e_code = er_frameset_to_json(output_frames, val_products_sof);
                else
                    e_code = er_frameset_to_text(output_frames, val_products_sof);

                cpl_frameset_delete(output_frames);
                if (e_code != 0)
                {
                    cpl_msg_error (er_func,"Cannot create output sof file, "
                                   "status = %d",
                                   e_code);
                }
            }
            else
            {
                cpl_msg_warning (er_func,"Writing of output sof omitted"
                                 " due to previous errors");
            }
        }

        /* Run Plugin Deinitialisation...  */

        plugin_deinit:

        plugin_func_deinit = cpl_plugin_get_deinit ((cpl_plugin *) trecipe);

        if (plugin_func_deinit != NULL)
        {
            cpl_msg_set_domain (cpl_plugin_get_name ((cpl_plugin *) trecipe));
            e_code2 = plugin_func_deinit ((cpl_plugin *) trecipe);
            cpl_msg_set_domain (PACKAGE);
        }

        if (e_code2 == 0)
        {
            /* For now deinitialization is empty */
        }

        /* If out main err.code is "OKAY", then use the deinit one, */
        /* (otherwise, preserve the original error). */

        if (e_code == 0) e_code = e_code2;
    }


    if ((flag_time_plugin != 0) && (flag_plugin_found != 0))
    {
        double fsz = 0.0;
        int infile_count = 0;
        cpl_frame *frame = NULL;
        const char *input_name = NULL;

        plugin_time_end = cpl_test_get_walltime();
        if (plugin_time_end > plugin_time_start)
        {
            plugin_time_end -= plugin_time_start;   /* reuse for time difference */
        }
        else
        {
            plugin_time_end = -1.0;
        }

        /* loop again over frameset - now, all tags should be set by the recipe */

        if (holdme != NULL)
        {
            cpl_frameset_iterator *it = cpl_frameset_iterator_new(holdme);

            while ((frame = cpl_frameset_iterator_get(it)) != NULL)
            {
                if (cpl_frame_get_group (frame) == CPL_FRAME_GROUP_RAW)
                {            /* Get the filename as specified by the frame */
                    input_name = cpl_frame_get_filename (frame);
                    if (add_size (input_name,&fsz) != 0)        /* accumulate size (in bytes) of raw input frames */
                    {
                        cpl_msg_warning(er_func,"could not get size of %s\n",input_name);
                    }
                    else
                    {
                        infile_count ++;
                    }
                }
                cpl_frameset_iterator_advance(it, 1);
            }

            cpl_frameset_iterator_delete(it);
        }

        if ((plugin_time_end > 0.0) && (fsz > 0.0))
        {
            cpl_msg_info (er_func, "Recipe operation(s) took %14.3g seconds to complete.",
                          plugin_time_end);
        }

        if (fsz > 0.0)
        {
            fsz *= 0.000001;
            if (infile_count > 1)
                cpl_msg_info (er_func, "Total size of %d raw input frames  = %8.2f MB\n",infile_count,fsz);
            else
                cpl_msg_info (er_func, "Size of single raw input frame  = %8.2f MB\n",fsz);

            if (plugin_time_end > 0.0)
            {
                fsz /= plugin_time_end;
                cpl_msg_info (er_func, "=> processing rate of %8.2f MB/sec \n",fsz);
            }
        }

    }

    /* Terminate CPL messaging */

    cleanup:

    /* currently cpl_msg_stop_log () always returns CPL_ERROR_NONE ... 080613 */
    if (cpl_msg_stop_log () != CPL_ERROR_NONE)
        cpl_msg_error (er_func, "An error was encountered while closing the logfile");

    else
    {                       /*  Move log file...  */
        char  tmpbuf[FILEMAX+12];

        if (getcwd (tmpbuf, (size_t) FILEMAX) == NULL)
        {
            tmpbuf[0] = '.';
            tmpbuf[1] = '\0';
        }
        (void) strcat (tmpbuf, "/.logfile");

        if (fileutils_file_exists(tmpbuf) != 0)
        {                       /* .logfile exists ... */
            n = fileutils_copy(tmpbuf,log_file_name_full);  /* copy .logfile to "real" log file */
            if (n < 0)
            {
                (void) printf("we could not copy .logfile to %s (err-code = %d)\n",
                              log_file_name_full,n);
            }
            else
            {                       /* get rid of .logfile  */
                if (n == 0) (void) unlink(tmpbuf);
            }
        }
    }


    /* free memory again */

    if (help_str != NULL) cpl_free (help_str);

    if (log_file_name_full != NULL) cpl_free (log_file_name_full);
    if (plugin_conf_file_global != NULL) cpl_free (plugin_conf_file_global);
    if (plugin_conf_file_local != NULL) cpl_free (plugin_conf_file_local);

    if (val_recipe_dirs != NULL) cx_strfreev (val_recipe_dirs);
    if (list_of_pllib_names != NULL) er_stringarray_delete (list_of_pllib_names);

    if (tplugin != NULL) cpl_plugin_delete (tplugin);

    if ((module != NULL) && (f_val_unload_plugin != 0)) {
        lt_dlclose (module);
    }

    if (t2recipe != NULL)
    {                   /* we have a version 2 recipe! */
        cpl_frameset_delete (trecipe->frames);
        cpl_plugin_delete ((cpl_plugin *) t2recipe);
    }
    else if (trecipe != NULL)
    {                   /* trecipe->parameters handled by plugin */
        cpl_frameset_delete (trecipe->frames);
        cpl_plugin_delete ((cpl_plugin *) trecipe);
    }

    if (e_code == -99999) e_code = 0;
    return e_code;
}                               /* End of plugin_process_plugin() */




/**********************************************************************/
/**
 * @brief  Enlarge memory buffer
 * @param  fn          name of caller
 * @param  pptr        addr. of pointer to allocated memory buffer
 * @param  msize       new (increased) size for memory buffer
 *
 * @returns 0 if successfull, !=0 otherwise
 * 
 * This function frees the currently allocated space and uses the pointer
 * 'mem_pntr' to point to newly allocated memory
 */
/**********************************************************************/

void er_enlarge (const char *fn, char **pptr, int msize)

{
    void  *newptr, *ptr;


    ptr = (void *) *pptr;           /* get address to free */

    newptr = cpl_realloc(ptr,(size_t) msize);
    if (newptr == NULL)             /* couldn't get the memory... */
    {
        cpl_msg_error (fn, "Could not allocate %d bytes - fatal error", msize);
        exit (EXIT_FAILURE);
    }

    *pptr = (char *) newptr;
}

/**@}*/

/* End of file */


