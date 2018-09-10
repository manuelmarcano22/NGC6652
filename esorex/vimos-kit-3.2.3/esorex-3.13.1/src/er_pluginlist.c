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
#   include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include <dirent.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include "ltdl.h"

#include <cpl.h>

#include "er_macros.h"
#include "er_fileutils.h"
#include "er_stringarray.h"
#include "er_plugin.h"
#include "er_pluginlist.h"
#include "er_paramutils.h"
#ifdef ENABLE_PYTHON_RECIPES
#include "er_python.h"
#endif

extern int  MyOS;       /* 1 = Linux, 0 = Mac OSX */

/*
 * Symbolic constants used in the search_directory() function
 */

#define INDENT_INCR  2          /* Amount to increment per indentation */
#define INDENT_STR   "| "       /* String to use for each indentation */
#define ENTRY_STR    "+-"       /* String to use at the start of each entry */
#define SEP_STR      "/"        /* Directory separator ("/" on UNIX) */

#ifdef ENABLE_PYTHON_RECIPES

/* File extention string for Python modules. */
#define PYTHON_MODULE_FILE_EXT ".py"

#endif


/**
 * @defgroup esorex_pluginlist EsoRex Plugin Listing Functions
 *
 * EsoRex Plugin Listing Functions.
 *
 */

/**@{*/

/**********************************************************************/
/**
 * @brief
 *   Test whether a file name indicates a plugin and what kind it is.
 *
 * @param full_filename
 *   This is the full filename. That is, it includes the full path
 *   specification and any extension too.
 *
 * @returns
 *   1 if the specified file name could contain a dynamically loadable
 *   shared library plugin, 2 if it could be a Python module, or 0 if
 *   it does not appear to be a valid file name at all.
 *
 * This function tests a fully specified filename to see whether or not
 * it may contain a dynamic library or Python module.
 */
/**********************************************************************/

static int plugin_valid_filename (const char *full_filename)

{
    char *ext;                            /* File extension */

    /* If there was no filename, then return */
    if (full_filename == NULL) return 0;

    /* Get the file extension */
    ext = strrchr (full_filename, '.');

    /* Reject any file which doesn't have the correct file extension */
    if (ext == NULL) return 0;

    /* Check whether the filename extension indicates a valid module */
    if (strcmp(ext, LT_MODULE_EXT) == 0)
    {
        return 1;
    }

#ifdef ENABLE_PYTHON_RECIPES

    else if (strcmp(ext, PYTHON_MODULE_FILE_EXT) == 0)
    {
        return 2;
    }

#endif

    else
    {
        /* Invalid file extension. */
        return 0;
    }
}


/**********************************************************************/
/**
 * @brief
 *   Test whether a file is a loadable dynamic library plugin file
 *
 * @param full_filename
 *   This is the full filename. That is, it includes the full path
 *   specification and any extension too.
 *
 * @returns
 *   TRUE if the specified module is valid, or FALSE if not.
 *
 * This function tests a fully specified filename to see whether or not
 * it may be loaded as a dynamic library module.
 */
/**********************************************************************/

static int plugin_valid_sharedlib (const char *full_filename)

{
    lt_dlhandle module = NULL;            /* Handle to the dynamic library */

    int (*get_plugin_list) (cpl_pluginlist *) = NULL;

    /*
     * Attempt to open the file as a dynamic library module
     */
    module = lt_dlopen (full_filename);

    /* If this didn't load properly, then return */
    if (!module)
    {
        const char  *dl_errstring;

        dl_errstring = lt_dlerror();
        cpl_msg_warning (er_func,
                         "lt_dlopen (%s) returned NULL,\nlt_dlerror() = %s\n",
                         full_filename,dl_errstring);
        return 0;
    }

    /* Check that it has the required plugin function */
    get_plugin_list = (int (*)()) lt_dlsym (module, "cpl_plugin_get_info");
    lt_dlclose (module);            /* Close the module */

    if (get_plugin_list == NULL) return 0;

    /* If we get this far, then it was obviously a valid plugin file */
    return 1;
}

#ifdef ENABLE_PYTHON_RECIPES

/**********************************************************************/
/**
 * @brief Checks if a file name corresponds to a Python module.
 *
 * @param[in] filename   The file name string to check.
 *
 * Returns 1 if the file name has a file extension corresponding to a Python
 * module and 0 otherwise.
 */
/**********************************************************************/

static int is_python_module (const char * filename)
{
    assert(filename != NULL);

    size_t len = strlen(filename);
    size_t extlen = strlen(PYTHON_MODULE_FILE_EXT);
    if (len < extlen) return 0;
    if (strcmp(filename + (len - extlen), PYTHON_MODULE_FILE_EXT) == 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

#endif

/**********************************************************************/
/**
 * @brief Recursively searches a directory tree for plugin file names
 *
 * @param dirname          The name of the directory to search.
 * @param list_of_pllibs   A list of filenames of possible plugin shared
 *                         libraries. If this function finds any more,
 *                         then they are added to this list.
 * @param list_of_pypllibs A list of filenames of possible Python plugin
 *                         modules. If this function finds any more,
 *                         then they are added to this list.
 *
 * Scans all the files in the specified directory. If any of these are
 * regular files, then another function is called to check whether or
 * not they are a valid plugin file name. If so, then this filename is
 * added to one of the provided @a list_of_pllibs or @a list_of_pypllibs,
 * depending on the type of plugin found. In the case that one of the
 * scanned files is actually a directory, then this function is
 * (recursively) called on that function. That way, eventually, all
 * files in or below the original directory are searched.
 *
 *
 *
 * @todo Implement some check to prevent endless looping, for example,
 *       in the event of an unforeseen error during the recursive calling
 *       of this function.
 *
 */
/**********************************************************************/

static void search_directory (const char * dirname, 
                              er_stringarray_t * list_of_pllibs,
                              er_stringarray_t * list_of_pypllibs )

{
    struct stat finfo;      /* File information */
    struct dirent *entry;       /* Entry in a given directory listing */

    DIR *dp;            /* Pointer to a directory */
    char *fqfn;                     /* Fully qualified filename (path+file+ext) */

    static int indent = 0;      /* Level of indentation */
    int valid;          /* Is the specified file a valid plugin file */

    int filename_len = 4096;    /* the "usual" FILENAME_MAX */
    int status;
    int ldir, lentry;


    assert(dirname != NULL);
    assert(list_of_pllibs != NULL);
    assert(list_of_pypllibs != NULL);


    dp = opendir (dirname);         /* Open the directory */
    if (dp == NULL)
    {
        cpl_msg_warning (er_func, "Unable to open directory %s", dirname);
        return;
    }


    ldir = (int) strlen(dirname);
    fqfn = cpl_malloc ((size_t) filename_len);
    if (fqfn == NULL)
    {
        closedir(dp);
        return;
    }

    while ((entry = readdir (dp)) != NULL)
    {
        /* Ignore hidden files, current directory and parent directory */
        if (entry->d_name[0] == '.') continue;

        /* We have a directory, so recursively search it with this function */
        lentry = (int) strlen(entry->d_name);
        if ((lentry+ldir+4) > filename_len) /* ensure enough memory is allocated */
        {
            filename_len = (lentry+ldir+4) * 2;
            er_enlarge("search_directory",&fqfn,filename_len);
        }
        (void) sprintf (fqfn, "%s%s%s", dirname, SEP_STR, entry->d_name);

        /* Get information about the current entry (file) */
        status = stat (fqfn, &finfo);
        if (status == 0)
        {
            if (S_ISREG (finfo.st_mode))
            {
                valid = plugin_valid_filename (fqfn);
                if (valid == 1)
                {
                    /* If we get this far, then we've possibly found a shared
                       library plugin */
                    /* Add it to the list of plugins that we have so far */

                    er_stringarray_append (list_of_pllibs, fqfn);
                }

#ifdef ENABLE_PYTHON_RECIPES

                else if (valid == 2)
                {
                    /* Found a Python module file instead. */
                    er_stringarray_append (list_of_pypllibs, fqfn);
                }

#endif

            }
            else if (S_ISDIR (finfo.st_mode))
            {
                /* Recursively search any subdirectory that we find */
                search_directory (fqfn, list_of_pllibs, list_of_pypllibs);
            }

        }                         /* End if valid status */
    }                            /* End of while */

    /* Release the memory allocated to the filename */

    cpl_free (fqfn);

    indent -= INDENT_INCR;      /* Decrement the indentation again */
    closedir (dp);

}                               /* End of search_directory() */



/**********************************************************************/
/**
 * @brief Adds items to a list such that they are not duplicated.
 *
 * @param src  The list of items to add.
 * @param dest New items will be appended to this list.
 *
 * @returns Pointer to stringarray
 *
 * This function will add all items in the @a src string list to the
 * @a dest list. But only those items are added that are not already in
 * the @a dest list. i.e. it keeps the items in @a dest unique.
 */
/**********************************************************************/

static void append_list_unique(er_stringarray_t * src, er_stringarray_t * dest)
{
    int i1 = 0, i2 = 0;
    for (i1=0; i1<er_stringarray_size (src); i1++)
    {
        char * src_item = er_stringarray_get (src, i1);
        if (src_item != NULL)   /* Skip any empty slots in the list. */
        {
            int duplicate = 0;
            for (i2=0; i2<er_stringarray_size (dest); i2++)
            {
                if (strcmp (src_item, er_stringarray_get (dest, i2)) == 0)
                {
                    duplicate = 1;
                }
            }
            if (duplicate == 0)
            {
                er_stringarray_append (dest, src_item);
            }
        }
    }
}


/**********************************************************************/
/**
 * @brief Creates a string array of plugin libraries in a given directory
 *
 * @param directories_to_search  An array of directory names.
 *
 * @returns Pointer to stringarray
 *
 * This function establishes new stringarrays, @c list_of_pllibs and
 * @c list_of_pypllibs. It then calls the @c search_directory function to
 * recursively search the specified @a directory_to_search .
 * @c search_directory will add the filenames of any found plugin files to
 * @c list_of_pllibs or @c list_of_pypllibs depending on the type of plugin.
 * When this operation is complete, the function will concolidate the two
 * lists of plugins into a single list and return this as a stringarray.
 *
 */
/**********************************************************************/

er_stringarray_t *er_pluginlist_create_list (char ** directories_to_search)

{
    er_stringarray_t *list_of_pllibs;   /* List of plugin shared libraries */
    er_stringarray_t *list_of_pypllibs; /* List of Python plugin modules */
    er_stringarray_t *reduced_list;  /* List of all plugin library file names */

    char **dir_str;                         /* Pointer to the list of directories */

    int i = 0;


    /* Obviously, we need a directory to search, so bail out if there isn't one */

    if (directories_to_search == NULL) return NULL;


    /* Create new arrays of strings to contain any plugin filenames */

    if ((list_of_pllibs = er_stringarray_new ()) == NULL) return NULL;
    if ((list_of_pypllibs = er_stringarray_new ()) == NULL) return NULL;

    for (dir_str = directories_to_search; *dir_str != NULL; dir_str++)
    {
        /*
         * Recursively search the specified directory, adding any found
         * plugins to the list
         */
        search_directory (er_fileutils_tilde_replace (*dir_str),
                          list_of_pllibs, list_of_pypllibs);
    }

    /* Double check the shared libraries, making sure they really are valid
       plugins for EsoRex. Remove any file names that are not valid plugins. */
    for (i=0; i<er_stringarray_size(list_of_pllibs); ++i)
    {
        if (! plugin_valid_sharedlib(er_stringarray_get(list_of_pllibs, i)))
        {
            er_stringarray_remove(list_of_pllibs, i);
        }
    }

#ifdef ENABLE_PYTHON_RECIPES

    /* We also have to check the Python modules. But this is done implicitly
       by er_python_load_modules(). It will also remove file names that are
       not valid Python modules for EsoRex. */
    if (er_python_load_modules(list_of_pypllibs) != CPL_ERROR_NONE)
    {
        /* At this point we actually ignore errors and just print the
           information as debug messages instead. This is so that we do not
           block handling of compiled recipies if there are problems with the
           Python ones. */
        cpl_msg_debug(er_func, "%s", cpl_error_get_message());
        cpl_error_reset();
    }

#endif

    /* Create a new array of strings for the consolidated list */
    reduced_list = er_stringarray_new ();
    if (reduced_list != NULL)
    {
        /* Consolidate the list - there may be duplicate entries */
        append_list_unique(list_of_pllibs, reduced_list);
        append_list_unique(list_of_pypllibs, reduced_list);
    }

    er_stringarray_delete (list_of_pypllibs);
    er_stringarray_delete (list_of_pllibs);

    /* Return the (possibly empty) list of plugin filenames */

    return reduced_list;

}                               /* End of er_pluginlist_create_list() */


/**********************************************************************/
/**
 * @brief  Gets the full pathname of the library containing the named plugin
 *
 * @param list_of_pllibs String array of plugin libraries
 * @param plugin_name    Name to look for in pllibs given in string array
 * @param[out] libpath   Filled with the library path string
 *
 * @returns the length of the library path found,
 *          or 0, in the case that no suitable library was found
 *
 * The function takes a list of plugin library path names, and opens
 * all of them searching for the named plugin: @c plugin_name . While 
 * doing this, the funciton also checks for the existence of multiple 
 * plugins with the same name. Should multiples be found, then the 
 * function will print a warning message to this effect, indicating 
 * which one it chose to use.
 *
 * The choice is made based on the version number. If the same version 
 * number exists multiple times, then the first one it finds in the 
 * library path is used.
 *
 */
/**********************************************************************/

int er_pluginlist_get_libpath (er_stringarray_t * list_of_pllibs,
                               const char * plugin_name, char *libpath)

{
    char *cptr, fully_qualified_library_name[MAXSTRLENCONF];

    lt_dlhandle module = NULL;

    cpl_pluginlist *pl_list = NULL;
    cpl_plugin *pl = NULL;

    int num_copies = 0, i = 0;
    int num_best = 0;
    int  m, e_code;

#ifdef ENABLE_PYTHON_RECIPES
    int is_shared_lib = 0;
#endif

    int (*get_plugin_list) (cpl_pluginlist *) = NULL;


    unsigned long this_version = 0;
    unsigned long best_version = 0;


    fully_qualified_library_name[0] = '\0';

    if (list_of_pllibs == NULL)
    {
        cpl_msg_warning (er_func, "Empty plugin list");
        return 0;
    }


    for (i=0; i< er_stringarray_size (list_of_pllibs); i++)
    {
        cptr = er_stringarray_get (list_of_pllibs, i);
        m = (int) strlen(cptr);
        if (m > (MAXSTRLENCONF-1))
        {
            cpl_msg_error (er_func, "size of plugin lib > %d ...", MAXSTRLENCONF);
            return 0;
        }

        (void) strcpy(fully_qualified_library_name, cptr);

#ifdef ENABLE_PYTHON_RECIPES

        is_shared_lib = !is_python_module(fully_qualified_library_name);

        if (is_shared_lib)
        {

#endif

            /* try to open library */
            module = lt_dlopen (fully_qualified_library_name);
            if (!module)
            {
                cpl_msg_error (er_func, "Could not open %s: %s",
                               fully_qualified_library_name, lt_dlerror ());
                return 0;
            }

            /* get plugin list  */

            get_plugin_list = (int (*)()) lt_dlsym (module, "cpl_plugin_get_info");

            /* Cannot use lt_dlerror due to bug in 1.4.2 version of libtool */

            if (get_plugin_list == NULL)
            {
                cpl_msg_error (er_func, "Could not find cpl_plugin_get_info() "
                               "index function for plugin library (%s): %s",
                               fully_qualified_library_name, lt_dlerror ());
                lt_dlclose(module);
                return 0;
            }

#ifdef ENABLE_PYTHON_RECIPES

        }
        else
        {
            /* In this case we are dealing with a Python module. Need to select
               the Python module and assign the get_plugin_list function pointer
               to the special interfacing function that will communicate with
               the Python plugins. */
            if (er_python_select_module(fully_qualified_library_name)
                    != CPL_ERROR_NONE)
            {
                cpl_msg_error(er_func, "%s", cpl_error_get_message());
                return 0;
            }
            get_plugin_list = er_python_get_plugin_list;
        }

#endif

        /* Make a pluginlist containing all the plugins in the library */

        pl_list = cpl_pluginlist_new ();
        e_code = get_plugin_list (pl_list);
        if (e_code != 0)
        {
            cpl_msg_warning (er_func,
                             "Unexpected error (%d) in recovering a pluginlist from "
                             "library '%s'", e_code,
                             fully_qualified_library_name);
        }

        /* Check the "one-plugin-per-library" assumption */

        m = cpl_pluginlist_get_size (pl_list);
        if (m == 0)
        {
            cpl_msg_warning (er_func, "No plugins contained within "
                             "library '%s'", fully_qualified_library_name);
        }
        else if (m > 1)
        {
            cpl_msg_warning (er_func, "Multiple plugins contained within "
                             "library '%s'", fully_qualified_library_name);
        }

        /* Once we have a valid library, fetch the plugin out of it */

        pl = cpl_pluginlist_find (pl_list, plugin_name);

        /* Double check again, to make sure we got it */

        if (pl != NULL)
        {       /* Increment a counter of all instances that we've found */
            num_copies++;

            /* Get the version number of this plugin */

            this_version = cpl_plugin_get_version (pl);

            /* Make sure the version number is value (i.e. non-zero) */

            if (this_version == 0)
            {
                cpl_msg_warning (er_func,
                                 "Recipe '%s' (in library '%s') has an invalid version number",
                                 plugin_name, fully_qualified_library_name);
            }

            /* Now see if this version is later than the best we have so far */

            if (this_version > best_version)
            {                   /* size check done above ... */
                (void) strcpy(libpath, fully_qualified_library_name);
                best_version = this_version;
                num_best = 1;
            }
            else if (this_version == best_version)
            {
                num_best++;
            }
        }
        cpl_pluginlist_delete (pl_list);

#ifdef ENABLE_PYTHON_RECIPES

        if (is_shared_lib)
        {
            lt_dlclose(module);
        }

#else

        lt_dlclose(module);

#endif

    }

    if (num_best > 1)
    {
        cpl_msg_error (er_func,
                       "Multiple copies (%d) of the latest version of recipe '%s' were found. "
                       "Using path '%s' ('%s' version %lu)", num_best, plugin_name, libpath,
                       plugin_name, best_version);
    }
    else if (num_copies > 1)
    {
        cpl_msg_warning (er_func,
                         "Older copies (%d) of recipe '%s' were also found. "
                         "Using latest version, path '%s' ('%s' version %lu)",
                         (num_copies - 1), plugin_name, libpath, plugin_name, best_version);
    }

    return ((int)strlen(libpath));

}                               /* End of er_pluginlist_get_libpath() */


/**********************************************************************/
/**
 * @brief  Gets a plugin, of a given name, from a given library.
 *
 * @param library_name   Path/filename of the shared library
 * @param plugin_name    Name of plugin to recover.
 * @param module         (output) A handle to the plugin library. This
 *                       handle should be closed (with lt_dlclose) when
 *                       the returned plugin is no longer used. If the
 *                       requested plugin could not be found or we are
 *                       dealing with a Python plugin, then this is set
 *                       to NULL.
 *
 * @returns Pointer to Plugin, or NULL if the requested plugin could not 
 *          be found
 *
 */
/**********************************************************************/

cpl_plugin *er_pluginlist_get_plugin (const char * library_name,
                                      const char * plugin_name,
                                      lt_dlhandle *module)

{
    const char *fn = "er_pluginlist_get_plugin";
    cpl_pluginlist *pl_list = NULL;
    cpl_plugin *pl = NULL, *pl_copy = NULL;

    int mm, e_code = 0;
    int (*get_plugin_list) (cpl_pluginlist *) = NULL;


    if (library_name == NULL)
    {
        cpl_msg_error (fn, "NULL pointer received instead of library_name");
        return NULL;
    }
    if (plugin_name == NULL)
    {
        cpl_msg_error (fn, "NULL pointer received instead of plugin_name");
        return NULL;
    }

    pl_list = cpl_pluginlist_new ();

#ifdef ENABLE_PYTHON_RECIPES

    int is_shared_lib = !is_python_module(library_name);

    if (is_shared_lib)
    {

#endif

        /* Try to open the module */

        *module = lt_dlopen (library_name);
        if (!(*module))
        {
            cpl_msg_error (fn, "Could not open %s: %s", library_name, lt_dlerror ());
            return NULL;
        }


        /* get plugin list */

        get_plugin_list = (int (*)()) lt_dlsym (*module, "cpl_plugin_get_info");


        /* lt_dlclose (module) will be done in caller routine plugin_process_plugin() */

        if (get_plugin_list == NULL)
        {
            cpl_msg_error (fn, "Could not find cpl_plugin_get_info() "
                           "index function for plugin library (%s): %s",
                           library_name, lt_dlerror ());
            return NULL;
        }

#ifdef ENABLE_PYTHON_RECIPES

    }
    else
    {
        /* In this case we are dealing with a Python module. Need to select the
           module and assign the get_plugin_list function pointer to the special
           interfacing function that will communicate with the Python plugins.
           */
        module = NULL;
        if (er_python_select_module(library_name) != CPL_ERROR_NONE)
        {
            cpl_msg_error(er_func, "%s", cpl_error_get_message());
            return NULL;
        }
        get_plugin_list = er_python_get_plugin_list;
    }

#endif

    pl_copy = NULL;
    e_code = get_plugin_list (pl_list);

    if (e_code != 0)
    {
        cpl_msg_warning (fn,
                         "Unexpected error (%d) in recovering a pluginlist from "
                         "library '%s'", e_code, library_name);
        goto end_of_it;
    }

    /* Check the "one-plugin-per-library" assumption */

    mm = cpl_pluginlist_get_size (pl_list);
    if (mm > 1)
    {
        cpl_msg_warning (fn, "Multiple plugins contained within "
                         "library '%s'", library_name);
    }
    else if (mm < 1)
    {
        cpl_msg_warning (fn, "No plugins contained within "
                         "library '%s'", library_name);
        goto end_of_it;
    }


    /* Now search for the plugin in the specified library */

    pl = cpl_pluginlist_find (pl_list, plugin_name);


    /* If found, make a copy of it (remains NULL otherwise) */

    if (pl != NULL)
    {
        pl_copy = cpl_plugin_new ();
        e_code = cpl_plugin_copy (pl_copy, pl);
        if (e_code != CPL_ERROR_NONE)
        {
            cpl_msg_error (fn, "Failed to make internal copy of plugin");
        }
    }

    /* Delete the plugin list, and return the pointer to the copy */

    end_of_it:
    cpl_pluginlist_delete (pl_list);

    return pl_copy;
}                               /* End of er_pluginlist_get_plugin() */


/**********************************************************************/
/**
 * @brief  Creates a string array of the names of available plugins.
 *
 * @param  list_of_pllibs        String array of pllibs
 * @retval list_of_plugin_names  String array of plugin names
 *
 * This function creates a string array of all plugins contained in list
 * of pllibs specified in a string array. The function does not display
 * these in any way.
 */
/**********************************************************************/

void er_pluginlist_create_cache (er_stringarray_t * list_of_pllibs,
                                 er_stringarray_t * list_of_plugin_names)

{
    const char *fn = "er_pluginlist_create_cache";

    char *cptr, fully_qualified_library_name[MAXSTRLENCONF];

    cpl_plugin *tplugin = NULL;

    lt_dlhandle module = NULL;

    int m, e_code, i = 0;
    int (*get_plugin_list) (cpl_pluginlist *) = NULL;

#ifdef ENABLE_PYTHON_RECIPES
    int is_shared_lib = 0;
#endif


    if (list_of_pllibs == NULL) return;

    fully_qualified_library_name[0] = '\0';

    for (i = 0; i < er_stringarray_size (list_of_pllibs); i++)
    {
        cpl_pluginlist *pl_list = cpl_pluginlist_new ();
        cptr = er_stringarray_get (list_of_pllibs, i);

        m = (int) strlen(cptr);
        if (m > (MAXSTRLENCONF-1))
        {
            cpl_msg_error (er_func, "size of plugin lib > %d ...", MAXSTRLENCONF);
            return;
        }
        (void) strcpy(fully_qualified_library_name,cptr);

#ifdef ENABLE_PYTHON_RECIPES

        is_shared_lib = !is_python_module(fully_qualified_library_name);

        if (is_shared_lib)
        {

#endif

            /* Try to open the module */

            module = lt_dlopen (fully_qualified_library_name);
            if (!module)
            {
                cpl_msg_error (fn, "Could not open %s: %s",
                               fully_qualified_library_name, lt_dlerror ());
                cpl_pluginlist_delete (pl_list);
                return;
            }

            /* get plugin list */

            get_plugin_list = (int (*)()) lt_dlsym (module, "cpl_plugin_get_info");

            if (get_plugin_list == NULL)
            {
                cpl_msg_error (fn, "Could not find cpl_plugin_get_info()"
                               " of plugin library %s",
                               fully_qualified_library_name);
                cpl_pluginlist_delete (pl_list);
                lt_dlclose(module);
                return;
            }

#ifdef ENABLE_PYTHON_RECIPES

        }
        else
        {
            /* In this case we are dealing with a Python module. Need to select
               the Python module and assign the get_plugin_list function pointer
               to the special interfacing function that will communicate with
               the Python plugins. */
            if (er_python_select_module(fully_qualified_library_name)
                    != CPL_ERROR_NONE)
            {
                cpl_msg_error(er_func, "%s", cpl_error_get_message());
                cpl_pluginlist_delete (pl_list);
                return;
            }
            get_plugin_list = er_python_get_plugin_list;
        }

#endif

        e_code = get_plugin_list (pl_list);
        (void) e_code;  /* suppress compiler warning. */

        tplugin = cpl_pluginlist_get_first (pl_list);
        while (tplugin != NULL)
        {
            if (cpl_plugin_get_api(tplugin) == CPL_PLUGIN_API)
            {
                if (cpl_plugin_get_type(tplugin) & CPL_PLUGIN_TYPE_RECIPE)
                {
                    er_stringarray_append(list_of_plugin_names,
                                          (char *)cpl_plugin_get_name(tplugin));
                }
            }

            tplugin = cpl_pluginlist_get_next(pl_list);

        }                         /* End of while() */

        cpl_pluginlist_delete(pl_list);

#ifdef ENABLE_PYTHON_RECIPES

        if (is_shared_lib)
        {
            lt_dlclose(module);
        }

#else

        lt_dlclose(module);

#endif

    }                            /* End of for(each name in list_of_pllibs) */

    return;

}                               /* End of er_pluginlist_create_cache() */


/**********************************************************************/
/**
 * @brief  Neatly print a list of all the plugins.
 *
 * @param  list_of_pllibs   String array of plugin library names.
 *
 * This function will take a stringarray, @a list_of_pllibs , which
 * contains the names of all the plugin files (i.e. the dynamic libraries)
 * which are visible within the current path. It then extracts from these
 * files the plugins themselves, and makes use of CPL calls to access
 * the name and synopsis of each individual file, before printing them.
 */
/**********************************************************************/

void er_pluginlist_print_list (er_stringarray_t * list_of_pllibs)

{
    char *ccptr, fully_qualified_library_name[MAXSTRLENCONF];

    cpl_plugin *tplugin = NULL;

    lt_dlhandle module = NULL;

    int (*get_plugin_list) (cpl_pluginlist *) = NULL;

    int m, e_code, i = 0;
    int num_plugins = 0;

#ifdef ENABLE_PYTHON_RECIPES
    int is_shared_lib = 0;
#endif


    fully_qualified_library_name[0] = '\0';

    if (list_of_pllibs == NULL) return;


    (void) printf("List of Available Recipes :\n\n");

    /* Loop through each dynamic library file that we have */

    for (i=0; i < er_stringarray_size(list_of_pllibs); i++)
    {
        cpl_pluginlist *pl_list = cpl_pluginlist_new ();

        ccptr = er_stringarray_get (list_of_pllibs, i);
        m = (int) strlen(ccptr);
        if (m > (MAXSTRLENCONF-1))
        {
            cpl_msg_error (er_func, "size of plugin lib > %d ...", MAXSTRLENCONF);
            return;
        }

        (void) strcpy(fully_qualified_library_name,ccptr);

#ifdef ENABLE_PYTHON_RECIPES

        is_shared_lib = !is_python_module(fully_qualified_library_name);

        if (is_shared_lib)
        {

#endif

            /* try to open library */
            module = lt_dlopen (fully_qualified_library_name);
            if (!module)
            {
                cpl_msg_error (er_func, "Could not open %s: %s",
                               fully_qualified_library_name, lt_dlerror ());
                cpl_pluginlist_delete (pl_list);
                return;
            }

            get_plugin_list = (int (*)()) lt_dlsym (module,
                                                    "cpl_plugin_get_info");

            if (get_plugin_list == NULL)
            {
                cpl_msg_error (er_func,
                               "Could not find the cpl_plugin_get_info()"
                               " function of plugin library %s",
                               fully_qualified_library_name);
                cpl_pluginlist_delete (pl_list);

                lt_dlclose(module);
                return;
            }

#ifdef ENABLE_PYTHON_RECIPES

        }
        else
        {
            /* In this case we are dealing with a Python module. Need to select
               the Python module and assign the get_plugin_list function pointer
               to the special interfacing function that will communicate with
               the Python plugins. */
            if (er_python_select_module(fully_qualified_library_name)
                    != CPL_ERROR_NONE)
            {
                cpl_msg_error(er_func, "%s", cpl_error_get_message());
                cpl_pluginlist_delete (pl_list);
                return;
            }
            get_plugin_list = er_python_get_plugin_list;
        }

#endif

        e_code = get_plugin_list (pl_list);
        if (e_code != CPL_ERROR_NONE)
        {
            cpl_msg_error (er_func, "Unable to read plugin list");
            cpl_pluginlist_delete (pl_list);

#ifdef ENABLE_PYTHON_RECIPES

            if (is_shared_lib)
            {
                lt_dlclose(module);
            }

#else

            lt_dlclose(module);

#endif

            return;
        }

        /* Get the first plugin in the list */

        tplugin = cpl_pluginlist_get_first (pl_list);

        /* Loop through all the plugins that we can find */

        while (tplugin != NULL)
        {       /* Check that we are printing appropriate plugins only */
            if (cpl_plugin_get_api (tplugin) != CPL_PLUGIN_API) continue;
            if (!(cpl_plugin_get_type (tplugin) & CPL_PLUGIN_TYPE_RECIPE)) continue;

            er_paramutils_print_key_desc ("", cpl_plugin_get_name (tplugin),
                                          cpl_plugin_get_synopsis (tplugin));

            num_plugins++;      /* Increment the plugin counter */

            /* Get the next plugin in the list */
            tplugin = cpl_pluginlist_get_next (pl_list);
        }

        cpl_pluginlist_delete (pl_list);

#ifdef ENABLE_PYTHON_RECIPES

        if (is_shared_lib)
        {
            lt_dlclose(module);
        }

#else

        lt_dlclose(module);

#endif

    }                            /* End of for(each name in list_of_pllibs) */

    /* Print a warning message if no plugins were listed */

    if (num_plugins == 0)
        (void) printf("  No recipes were found in the specified recipe directory.\n");

    (void) printf("\n");

    return;

}                               /* End of er_pluginlist_print_list() */

/**@}*/


/* End of file */
