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
#include <config.h>
#endif

#ifdef HAVE_WORDEXP
#include <wordexp.h>
#endif


#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include <ctype.h>

#include <cpl.h>

#include "er_macros.h"
#include "er_plugin.h"
#include "er_fileutils.h"

#define DEV_BLOCKSIZE   4096

#define Xspace(x)  (x == ' ') || (x == '\t') || (x == '\n')


/**
 * @defgroup er_fileutils EsoRex File Utility Functions
 *
 * Utility functions for hadling files.
 *
 */

/**@{*/

#ifdef HAVE_WORDEXP

/**********************************************************************/
/**
 * @brief  Replaces any environment variables with literal strings
 *
 * @param  path      Input path
 *
 * @return A string with the replaced variables.
 *
 * This function performs the replacement of any environment
 * variables. It only works if the WORDEXP package is available
 * (detected at "configure" time and indicated as the HAVE_WORDEXP
 * symbolic constant). If this packages doesn't exist then this
 * function does nothing.
 */
/**********************************************************************/

char *er_fileutils_env_replace (const char *path)

{

    static char expanded_path[PATHSET_MAX + 1];
    char tpath[PATHSET_MAX + 1];

    wordexp_t pwordexp;


    (void) strncpy (expanded_path, path, PATHSET_MAX);


    /* Perform tilde and environment variable expansion */

    if (wordexp (path, &pwordexp, WRDE_NOCMD | WRDE_UNDEF) != 0)
    {
        cpl_msg_warning (er_func, "Unable to expand path '%s'\n", path);

        if (pwordexp.we_wordc > 0) wordfree (&pwordexp);
        return expanded_path;
    }


    if (!pwordexp.we_wordv[0] || strlen (pwordexp.we_wordv[0]) > PATHSET_MAX)
    {
        cpl_msg_warning (er_func, "Buffer overflow while expanding "
                         "path '%s'\n", path);
        wordfree (&pwordexp);
        return expanded_path;
    }

    (void) strcpy (expanded_path, pwordexp.we_wordv[0]);
    wordfree (&pwordexp);


    /* Return a pointer to the amended (if at all) string */
    return expanded_path;
}

#endif


/**********************************************************************/
/**
 * @brief  Replaces the leading tilde (~) in a name.
 *
 * @param name  The input file name/path string
 *
 * @return A string with the replaced tilde.
 *
 * This function takes an input string. If the string begins with a tilde
 * (~), then the function returns this string, but with the tilde replaced
 * by the home directory of the current user (using the $HOME environment
 * variable). If not, then the string is returned unmodified. This function
 * also cleans up any multiple or trailing forward-slashes.
 *
 * This function also performs the replacement of any environment
 * variables.
 *
 * OJO: this routine returns a pointer to a static string,
 *      if you need the returned filename later on, you must copy it
 *      to a string of PATHSET_MAX length
 *      otherwise the returned pointer may point later on to something
 *      very different   :-(
 */
/**********************************************************************/

const char *er_fileutils_tilde_replace (const char * name)

{
    static char str[PATHSET_MAX];
    char *marker, *cpp;

    int  len, lname;


    if (name == NULL) return NULL;

    lname = (int) strlen(name);

    /* Only do this, if we start with a tilde (~) */

    if (*name != '~')
    {
        if (lname > (PATHSET_MAX-1))
        {
            cpl_msg_error(er_func,
                          "Buffer overflow in filename '%s' - fatal error",
                          name);
            exit (EXIT_FAILURE);
        }
        return (name);
    }

    else
    {       /* get HOME env. var., if we start with a tilde */
        cpp = getenv ("HOME");
        if (cpp == NULL)
        {
            cpl_msg_error(er_func,
                          "Env. variable HOME not set, could not replace `~'");
            exit (EXIT_FAILURE);
        }

        len = (int) strlen(cpp);
        if (len > (PATHSET_MAX-1))
        {
            cpl_msg_error (er_func,
                           "Buffer overflow in filename '%s' - fatal error",
                           name);
            exit (EXIT_FAILURE);
        }
        (void) strcpy (str,cpp);

        len += lname;               /* +1 for \0 */
        if (len > PATHSET_MAX)          /* -1 for (name+1) below */
        {
            cpl_msg_error(er_func,
                          "Buffer overflow in filename '%s' - fatal error",
                          name);
            exit (EXIT_FAILURE);
        }

        (void) strcat (str, name + 1);
    }


    /* Remove any multiple slashes */


    marker = strstr (str, "//");
    while (marker != NULL)
    {
        memmove (marker, marker + 1, strlen (marker));
        marker = strstr (str, "//");
    }


    /* Remove any trailing slashes */

    marker = str - 1 + (int) strlen(str);       /* last char of `str' */
    if ((*marker) == '/') *marker = '\0';

#ifdef HAVE_WORDEXP
    return (er_fileutils_env_replace (str));
#else

    return (str);
#endif

}                            


/**********************************************************************/
/**
 * @brief  Replaces the leading dot (.) in a name.
 *
 * @param name  The input file name/path string
 *
 * @return A string with the replaced ".".
 *
 * This function takes an input string. If the string begins with a
 * single dot (.), then the function returns this string, but with the
 * dot replaced by the current working directory (as an absolute path)
 * using the $PWD environment variable. If not, then the string is
 * returned unmodified.
 *
 * OJO: this routine returns a pointer to a static string,
 *      if you need the returned filename later on, you must copy it
 *      to a string of PATHSET_MAX length
 *      otherwise the returned pointer may point later on to something
 *      very different   :-(
 */
/**********************************************************************/

const char *er_fileutils_dot_replace (const char * name)

{
    static char str[PATHSET_MAX];
    char  *cpp;

    int  len, lname;


    if (name == NULL) return NULL;

    lname = (int) strlen(name);

    /* if we don't start with a dot (.) just pass the name on */

    if (*name != '.')
    {
        if (lname > (PATHSET_MAX-1))
        {
            cpl_msg_error (er_func,
                           "Buffer overflow in filename '%s' - fatal error",
                           name);
            exit (EXIT_FAILURE);
        }

        return (name);
    }


    /* Get the environment variable PWD */

    cpp = getcwd(str, PATHSET_MAX);

    if (cpp == NULL)
    {
        cpl_msg_error(er_func, "Buffer overflow in filename '%s' - fatal error",
                      name);
        exit(EXIT_FAILURE);
    }

    len = (int) strlen(cpp);


    /* Handle the ".." case (add a "/." as the "." replacement) */

    if (*(name + 1) == '.')
    {
        if ((len+2) > (PATHSET_MAX-1))
        {
            cpl_msg_error (er_func,
                           "Buffer overflow in filename '%s' - fatal error",
                           name);
            exit (EXIT_FAILURE);
        }
        (void) strcat (str, "/.");
        len += 2;
    }

    len += lname;                   /* +1 for \0 */
    if (len > PATHSET_MAX)              /* -1 for (name+1) below */
    {
        cpl_msg_error (er_func, "Buffer overflow in filename '%s'", name);
        cpl_msg_error (er_func, "Fatal error replacing current working "
                       "directory symbol due to buffer overflow");
        exit (EXIT_FAILURE);
    }

    (void) strcat(str, name + 1);       /* Add the rest of the file name */
    return (str);
}                            


/**********************************************************************/
/**
 * @brief
 *   Checks whether the given directory actually exists
 *
 * @param directory_name Directory
 *
 * @returns TRUE if directory exists, FALSE otherwise
 */
/**********************************************************************/

int fileutils_directory_exists (const char * directory_name)

{
    DIR *dp;


    if (directory_name != NULL)
    {               /* no need to save the pointer, here! */
        dp = opendir (er_fileutils_tilde_replace (directory_name));
        if (dp != NULL)
        {
            (void) closedir (dp);
            return 1;           /* that's TRUE */
        }
    }

    return 0;           /* that's FALSE */
}       


/**********************************************************************/
/**
 * @brief
 *   Checks whether the given file actually exists
 *
 * @param file_name Filename
 *
 * @returns 1 if file exists, 0 otherwise
 */
/**********************************************************************/

int fileutils_file_exists (const char * file_name)

{
    int file_descriptor;

    const char *resname;


    if (file_name != NULL)
    {               /* no need to save the pointer, here! */
        resname = er_fileutils_dot_replace (er_fileutils_tilde_replace (file_name));

        file_descriptor = open (resname, O_RDONLY);
        if (file_descriptor > -1)
        {
            (void) close (file_descriptor);
            return 1;
        }
    }

    return 0;           /* that's FALSE */
} 


/**********************************************************************/
/**
 * @brief
 *   Creates a Fully Qualified File Name based on a directory name
 *   and a filename
 *
 * @param dir_name  Directory
 * @param file_name File
 *
 * @returns A newly allocated string containing FQFN
 */
/**********************************************************************/

char *fileutils_create_fqfname (char * dir_name, char * file_name)

{
    char *fqfn = NULL;

    int n;
    int ldn = 0, lfn = 0, lextra = 2;


    if (file_name == NULL) return NULL;

    if (dir_name != NULL)
    {
        ldn = (int) strlen (dir_name);
        if (dir_name[ldn] == '/') lextra = 1;
    }

    lfn = (int) strlen (file_name);
    n = sizeof (char) * (ldn + lfn + lextra);
    fqfn = (char *) cpl_malloc ((size_t) n);
    if (fqfn == NULL) return NULL;

    (void) strcpy (fqfn, dir_name);
    if (lextra == 2) (void) strcat (fqfn, "/");
    (void) strcat (fqfn, file_name);

    return fqfn;

}                               /* End of fileutils_create_fqfname() */


/**********************************************************************/
/**
 * @brief
 *   Extracts filename from a Fully Qualified File Name
 *
 * @param path FQFN
 *
 * @returns A newly allocated string containing Filename
 */
/**********************************************************************/

char *fileutils_fqfname_filename (const char * path)

{
    char  *fname;
    const char *pc;

    int len, j;


    if (path == NULL) return NULL;

    len = (int) strlen (path);
    j = len;

    pc = path + len;

    while ((j >= 0) && (*pc != '/'))
    {
        pc--;
        j--;
    }

    fname = (char *) cpl_calloc ((size_t) (len - j), (size_t) sizeof (char));
    if (fname == NULL) return NULL;

    (void) strncpy (fname, &path[j + 1], (size_t) (len-j-1));

    return fname;

}                               /* End of fileutils_fqfname_filename() */


/**********************************************************************/
/**
 * @brief
 *   Extracts directoryname from a Fully Qualified File Name
 *
 * @param path FQFN
 *
 * @returns A newly allocated string containing directoryname
 */
/**********************************************************************/

char *fileutils_fqfname_dirname (const char * path)

{
    char  *dname;
    const char *pc;

    int len, j;


    if (path == NULL) return NULL;

    len = (int) strlen (path);
    j = len;

    pc = path + len;

    while ((j >= 0) && (*pc != '/'))
    {
        pc--;
        j--;
    }

    dname = (char *) cpl_calloc ((size_t) (j+1), (size_t) sizeof (char));
    if (dname == NULL) return NULL;

    (void) strncpy (dname, path, (size_t) j);

    return dname;

} 


/*
 * @brief
 *   Copy a file.
 *
 * @param srcpath  Source file name.
 * @param dstpath  Destination file name.
 *
 * @return The function returns 0 if no error occurred, otherwise -1 is
 *    returned and errno is set appropriately.
 *
 * The function provides a very basic way to copy the source file
 * @em srcpath to a destination file @em dstpath.
 *
 * The implementation just uses @b read() and @b write() to do the job.
 * It is by far not comparable with the @b cp shell command. Actually it
 * just writes the source file contents into a new file with the
 * appropriate output name.
 *
 * If an error occurs the destination file is removed before the function
 * returns.
 */

int fileutils_copy (const char * srcpath, const char * dstpath)

{
    char *buf;

    int src, dst;
    int rbytes = 0;
    int wbytes = 0;
    int blksize = DEV_BLOCKSIZE;

    struct stat sb, db;


    if ((stat(srcpath,&sb) == 0) && (stat(dstpath,&db) == 0))
    {
        if (sb.st_ino == db.st_ino)
            return 99;  /* if inodes are the same we are done already... */
    }

    if ((src = open (srcpath, O_RDONLY)) == -1) return (-1);

    if ((fstat (src, &sb) == -1) || (!S_ISREG (sb.st_mode)))
    {
        (void) close (src);
        return -2;
    }

    if ((dst = open (dstpath, O_CREAT | O_WRONLY | O_TRUNC, sb.st_mode)) == -1)
    {
        (void) close (src);
        return -3;
    }

    if ((fstat (dst, &db) == -1) || (!S_ISREG (db.st_mode)))
    {
        (void) close (src);
        (void) close (dst);
        (void) unlink (dstpath);
        return -4;
    }

#ifdef HAVE_ST_BLKSIZE
    blksize = db.st_blksize;
#else
#   ifdef DEV_BSIZE
    blksize = DEV_BSIZE;
#   endif
#endif

    if ((buf = (char *) cpl_malloc ((size_t)blksize)) == NULL)
    {
        (void) close (src);
        (void) close (dst);
        (void) unlink (dstpath);
        return -5;
    }

    while ((rbytes = (int) read (src, buf, (size_t)blksize)) > 0)
    {
        if ((wbytes = (int) write (dst, buf, (size_t)rbytes)) != rbytes)
        {
            wbytes = -1;
            break;
        }
    }

    (void) close (src);
    (void) close (dst);
    cpl_free (buf);


    if ((rbytes == -1) || (wbytes == -1))
    {
        (void) unlink (dstpath);
        return -6;
    }

    return 0;

}


/*
 * @brief
 *   Move a file.
 *
 * @param srcpath  Source file name.
 * @param dstpath  Destination file name.
 *
 * @return The function returns 0 if no error occurred, otherwise -1 is
 *    returned and errno is set appropriately.
 *
 * The function moves the source file @em srcpath to the given
 * destination path @em dstpath.
 *
 * If an error occurs the destination file is removed before the function
 * returns and the source file is left untouched.
 */



int fileutils_move (const char *srcpath, const char *dstpath)

{
    int  ii;

    struct stat sb;


    if ((ii = fileutils_copy (srcpath, dstpath)) != 0)
    {
        if (ii == 99)       /* different path name, but same inode */
            return 99;      /* => it's the same file - do nothing */
        else
            return (-2);
    }


    /*
     * Remove the source, but if the source file cannot be checked or is not
     * writable revert to the original state, i.e. remove the file copy.
     */

    if (stat (srcpath, &sb) == -1 || !(sb.st_mode & S_IWUSR))
    {
        (void) unlink (dstpath);
        return -1;
    }

    (void) unlink (srcpath);
    return 0;

} 


/*
 * @brief
 *   Export a pipeline product file to the archive.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * @param link_path  Path to the export directory.
 * @param file_path  Absolute path to the product file to export.
 *   
 * The input file given by @em file_path is located in the product
 * directory and it is verified that the file exists and is readable.
 * The product file is exported by creating a link (see the note below)
 * to the product file in the pipeline export directory given by @em link_path.
 * The name of the link will be the same as the name of the product file.
 * 
 * The export directory name @em link_path should not contain a trailing slash.
 * 
 * @warning
 *   In case symbolic links are not supported by the system the product
 *   file is simply copied. To copy the file is safer, since hard
 *   links cannot be used across file systems.
 *
 * @note
 *   This function is a copy of the function pilRecExportProduct() from the
 *   VIMOS system, written by R. Palsa.
 *
 */

int er_fileutils_link (const char *link_path, const char *file_path)

{
    char linkname[FILEMAX + 1];
    const char *s;

    int  len;

    struct stat lb;


    /* Check if the input file is readable */

    if (access (file_path, R_OK))
    {
        cpl_msg_error (er_func, "Product file is unreadable: %s", file_path);
        return EXIT_FAILURE;
    }


    /* pull out the file name */

    s = strrchr (file_path, '/');           /* points to last '/' */
    if (s == NULL)
    {
        s = file_path;
    }
    else
    {
        ++s;
    }
    len = (int) strlen(s);
    len = (int) strlen(link_path) + len + 1;
    if (len > FILEMAX)
    {
        cpl_msg_error (er_func, "Buffer overflow - fatal error");
        return EXIT_FAILURE;
    }
    (void) snprintf (linkname, (size_t)FILEMAX, "%s/%s", link_path, s);


    /* check if link already there... */

    /* if (fileutils_file_is_there (linkname) == 1)  */
    if (fileutils_file_exists (linkname) == 1)
    {
        cpl_msg_warning (er_func, "the link %s already exists - will be overwritten!",linkname);
    }


    /*
     * To export the file to the archive usually a symbolic link to the
     * product file is created. But symbolic links are not conforming to
     * the ANSI C standard. If symbolic links are not available the file
     * is just copied. Hard links might be an alternative but they cannot
     * be used across file systems, so a real copy is created here.
     */

#if defined HAVE_LSTAT && defined HAVE_SYMLINK

    if (lstat (linkname, &lb))
    {
        if (errno != ENOENT)
        {
            cpl_msg_error (er_func, "Cannot get file status: %s", linkname);
            return EXIT_FAILURE;
        }
    }
    else
    {
        if (S_ISLNK (lb.st_mode))
            (void) unlink (linkname);
        else
        {
            cpl_msg_error (er_func, "Not a symbolic link: %s", linkname);
            return EXIT_FAILURE;
        }
    }

    if (symlink (file_path, linkname) != 0)
    {
        cpl_msg_error (er_func, "Cannot create symbolic link for %s", file_path);
        return EXIT_FAILURE;
    }

#else

    if (stat (linkname, &lb))
    {
        if (errno != ENOENT)
        {
            cpl_msg_error (er_func, "Cannot get file status: %s", linkname);
            return EXIT_FAILURE;
        }
    }
    else
        (void) unlink (linkname);

    if (fileutils_copy (file_path, linkname) < 0)
    {
        cpl_msg_error (er_func, "Cannot copy %s", file_path);
        return EXIT_FAILURE;
    }
#endif

    return EXIT_SUCCESS;

}                               /* End of er_fileutils_link() */


/*----------------------------------------------------------------------------*/
/**
  @brief  Determine whether a given file is a FITS file
  @param  name     The name of the fits
  @return 1 if the file is FITS, 0 if not, negative on error

 */
/*----------------------------------------------------------------------------*/

int er_fileutils_file_is_fits(const char *name)

{
    FILE *fp;

    char  *line, *keyw, *eqchar, *kval;

    int is_fits = 0;
    int  i, j, n;

    int mysscanf();




    if (!(fp = fopen (name, "r")))
    {
        cpl_msg_error (er_func, "Unable to open FITS file '%s'", name);
        return (-2);
    }

    line = (char *) cpl_calloc((size_t) 82, (size_t) 1);
    keyw = (char *) cpl_calloc((size_t) 80, (size_t) 1);
    eqchar = (char *) cpl_calloc((size_t) 80, (size_t) 1);
    kval = (char *) cpl_calloc((size_t) 80, (size_t) 1);


    /* Loop over first few lines in the file */

    for (i=0; i<10; i++)
    {
        line[0] = '\0';
        if (fgets (line, 80, fp) == NULL)
        {
            cpl_msg_error (er_func, "Failed to read file '%s'", name);
            goto close_file;
        }

        j = 0;
        while (Xspace(line[j]))     /* omit leading space */
        {
            j ++;
            if (line[j] == '\0')
            {
                cpl_msg_error (er_func, "%s not a FITS file!", name);
                goto close_file;
            }
        }

        /*
         * look for pattern 'SIMPLE  =                    T' at the
         * beginning of the file which indicates that the file is a
         * FITS file.
         */

        n = mysscanf (&line[j], keyw, eqchar, kval);        /* split up line */
        /* printf("we get n = %d, keyw = %s, %s, %s\n",n,keyw,kval, eqchar); */

        if ((n == 3) && (strcmp(keyw,"SIMPLE") == 0)
                && (*eqchar == '=') && (*kval == 'T'))
        {
            is_fits = 1;
            goto close_file;
        }
    }

    /* first 10 lines indicate no FITS header ... */
    cpl_msg_error (er_func, "%s not a FITS file!", name);


    close_file:
    cpl_free(line);
    cpl_free(keyw);
    cpl_free(eqchar);
    cpl_free(kval);

    fclose (fp);
    return is_fits;
}

/**@}*/

/* End of file */
