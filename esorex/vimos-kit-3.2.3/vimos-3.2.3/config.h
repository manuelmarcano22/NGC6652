/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define if you have the `fits_get_hduaddrll' function */
/* #undef HAVE_FITS_GET_HDUADDRLL */

/* Define if you have the `fits_get_hduoff' function */
#define HAVE_FITS_GET_HDUOFF 1

/* Define to 1 if you have the `floor' function. */
#define HAVE_FLOOR 1

/* Define if you have the `getpwuid' function */
#define HAVE_GETPWUID 1

/* Define if you have the `getuid' function */
#define HAVE_GETUID 1

/* Define to 1 iff you have GSL */
#define HAVE_GSL 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 iff you have GSL */
#define HAVE_LIBGSL 1

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define to 1 if you have <pthread.h>. */
#define HAVE_PTHREAD_H 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define if you have the `strdup' function */
/* #undef HAVE_STRDUP */

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1


/*
 * The maximum length, in bytes of an input line. It is set to the
 * POSIX2 value.
 */

#ifndef LINE_LENGTH_MAX
#  define LINE_LENGTH_MAX  2048
#endif
            

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define if online support tools should be built */
/* #undef ONLINE_MODE */

/* Name of package */
#define PACKAGE "vimos"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "usd-help@eso.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "VIMOS Instrument Pipeline"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "VIMOS Instrument Pipeline 3.2.3"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "vimos"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "3.2.3"


/*
 * Define the maximum number of characters a filename including the
 * path may have, excluding the string terminator.
 */

#ifndef PATHNAME_MAX
#  define PATHNAME_MAX  4096
#endif
            

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "3.2.3"

/* VIMOS binary age */
#define VIMOS_BINARY_AGE 203

/* VIMOS binary version number */
#define VIMOS_BINARY_VERSION 30203

/* VIMOS interface age */
#define VIMOS_INTERFACE_AGE 0

/* VIMOS major version number */
#define VIMOS_MAJOR_VERSION 3

/* VIMOS micro version number */
#define VIMOS_MICRO_VERSION 3

/* VIMOS minor version number */
#define VIMOS_MINOR_VERSION 2

/* Plugin directory tree prefix */
#define VIMOS_PLUGIN_DIR "esopipes-plugins"

/* Absolute path to the plugin directory tree */
#define VIMOS_PLUGIN_PATH "NONE/lib/esopipes-plugins"

/* Absolute path to the sextractor configuration files */
#define VIMOS_SEXTRACTOR_CONFIG "/home/mmarcano/Documents/VIMOS/NGC6652/esorex/share/esopipes/vimos-3.2.3/config"

/* Absolute path to the sextractor executable */
#define VIMOS_SEXTRACTOR_PATH "/home/mmarcano/Documents/VIMOS/NGC6652/esorex/lib/vimos-3.2.3/bin"

/* Define if using the dmalloc debugging malloc package */
/* #undef WITH_DMALLOC */

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */


#ifndef HAVE_STRDUP
#  define strdup  cx_strdup
#endif
              
