/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Width of type 'char' */
/* #undef CHAR_BIT */

/* Log domain internally used by the library */
#define CX_LOG_DOMAIN "CxLib"

/* Define if thread support is enabled. */
#define CX_THREADS_ENABLED 1

/* A `va_copy()' style function */
#define CX_VA_COPY va_copy

/* Define if you have the `asprintf' function */
/* #undef HAVE_ASPRINTF */

/* Define if CHAR_BIT is defined in limits.h and equals 8 */
#define HAVE_CHAR_BIT 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <float.h> header file. */
#define HAVE_FLOAT_H 1

/* Define if you have the `fpathconf' function */
#define HAVE_FPATHCONF 1

/* Define to 1 if the system has the type `intmax_t'. */
#define HAVE_INTMAX_T 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the `localeconv' function. */
#define HAVE_LOCALECONV 1

/* Define to 1 if you have the <locale.h> header file. */
#define HAVE_LOCALE_H 1

/* Define to 1 if the system has the type `long double'. */
#define HAVE_LONG_DOUBLE 1

/* Define to 1 if the system has the type `long long int'. */
#define HAVE_LONG_LONG_INT 1

/* Define to 1 if `lstat' has the bug that it succeeds when given the
   zero-length file name argument. */
/* #undef HAVE_LSTAT_EMPTY_STRING_BUG */

/* Define to 1 if you have the `memcpy' function. */
#define HAVE_MEMCPY 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define if you have the `pathconf' function */
#define HAVE_PATHCONF 1

/* Define if printf default precision for format `g' is 1 (ISO C standard) or
   6 */
/* #undef HAVE_PRINTF_FLT_FMT_G_STD */

/* Define if printf format `%p' produces the same output as `%#x' or `%#lx' */
/* #undef HAVE_PRINTF_PTR_FMT_ALTERNATE */

/* Define if printf outputs `(nil)' when printing NULL using `%p' */
#define HAVE_PRINTF_PTR_FMT_NIL 1

/* Define if printf treats pointers as signed when using a sign flag */
#define HAVE_PRINTF_PTR_FMT_SIGNED 1

/* Define if printf outputs `(null)' when printing NULL using `%s' */
#define HAVE_PRINTF_STR_FMT_NULL 1

/* Define to 1 if you have <pthread.h>. */
#define HAVE_PTHREAD_H 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#define HAVE_PTRDIFF_T 1

/* Define if you have the `snprintf' function */
#define HAVE_SNPRINTF 1

/* Define to 1 if you have the `sprintf' function. */
#define HAVE_SPRINTF 1

/* Define to 1 if `stat' has the bug that it succeeds when given the
   zero-length file name argument. */
/* #undef HAVE_STAT_EMPTY_STRING_BUG */

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define if you have the `stpcpy' function */
/* #undef HAVE_STPCPY */

/* Define to 1 if you have the `strchr' function. */
#define HAVE_STRCHR 1

/* Define if you have the `strdup' function */
#define HAVE_STRDUP 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if `decimal_point' is a member of `struct lconv'. */
#define HAVE_STRUCT_LCONV_DECIMAL_POINT 1

/* Define to 1 if `thousands_sep' is a member of `struct lconv'. */
#define HAVE_STRUCT_LCONV_THOUSANDS_SEP 1

/* Define if you have the `sysconf' function */
#define HAVE_SYSCONF 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if the system has the type `uintmax_t'. */
#define HAVE_UINTMAX_T 1

/* Define to 1 if the system has the type `uintptr_t'. */
#define HAVE_UINTPTR_T 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if the system has the type `unsigned long long int'. */
#define HAVE_UNSIGNED_LONG_LONG_INT 1

/* Define to 1 if you have the <values.h> header file. */
#define HAVE_VALUES_H 1

/* Define if you have the `vasprintf' function */
/* #undef HAVE_VASPRINTF */

/* Define if you have an implementation of `va_copy()'. */
#define HAVE_VA_COPY 1

/* Define if you have an implementation of a `va_copy()' style function. */
#define HAVE_VA_COPY_STYLE_FUNCTION 1

/* Define if `va_lists' can be copied by value */
/* #undef HAVE_VA_LIST_COPY_BY_VALUE */

/* Define if you have the `vsnprintf' function */
#define HAVE_VSNPRINTF 1

/* Define if you have the C99 `vsnprintf' function. */
#define HAVE_VSNPRINTF_C99 1

/* Define if realloc(NULL,) works */
#define HAVE_WORKING_REALLOC 1

/* Define if you have an implementation of `__va_copy()'. */
/* #undef HAVE___VA_COPY */

/* Define to 1 if `lstat' dereferences a symlink specified with a trailing
   slash. */
#define LSTAT_FOLLOWS_SLASHED_SYMLINK 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Name of package */
#define PACKAGE "cext"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "cpl-help@eso.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "C Extension Library"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "C Extension Library 1.2.2"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "cext"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.2.2"

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* Define to 1 if the `S_IS*' macros in <sys/stat.h> do not work properly. */
/* #undef STAT_MACROS_BROKEN */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.2.2"

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

/* Define to the widest signed integer type if <stdint.h> and <inttypes.h> do
   not define. */
/* #undef intmax_t */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to the widest unsigned integer type if <stdint.h> and <inttypes.h>
   do not define. */
/* #undef uintmax_t */

/* Define to the type of an unsigned integer type wide enough to hold a
   pointer, if such a type exists, and if the system does not define it. */
/* #undef uintptr_t */


#ifndef HAVE_STRDUP
#  define strdup  cx_strdup
#endif
              
