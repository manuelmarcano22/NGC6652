# Process this file with autoconf to produce a configure script.

AC_INIT([VIMOS Pipeline Library], [0.5.0], [usg-help@eso.org], [libpil])
AC_PREREQ([2.59])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR([admin])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AM_MAINTAINER_MODE

LIBPIL_SET_PREFIX(${PIPE_HOME:-/usr/local})

# Immediately before every release do:
#-------------------------------------
#   if (the interface is totally unchanged from previous release)
#       REVISION++;
#   else {
#       /* interfaces have been added, removed or changed */
#       REVISION = 0;
#       CURRENT++;
#       if (any interfaces have been _added_ since last release)
#           AGE++;
#       if (any interfaces have been _removed_ or incompatibly changed)
#           AGE = 0;
#   }

LIBPIL_SET_VERSION_INFO([$VERSION])

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_LN_S

ESO_ENABLE_DEBUG(no)
ESO_ENABLE_STRICT(no)
ESO_PROG_CC_FLAG([fno-builtin], [CFLAGS="$CFLAGS -fno-builtin"])

ESO_CHECK_DOCTOOLS

AC_ENABLE_STATIC(no)
AC_ENABLE_SHARED(yes)

AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# Checks for libraries.
CPL_CHECK_CFITSIO
KAZLIB_CONVENIENCE_LIBRARY
AC_SUBST(INCKAZ)
AC_SUBST(LIBKAZ)
AC_SUBST(CFITSIO_INCLUDES)
AC_SUBST(CFITSIO_LDFLAGS)

AC_CHECK_LIB(m, pow, [LIBS="$LIBS -lm"])
AC_CHECK_LIB(socket, socket, [LIBS="$LIBS -lsocket"])
AC_CHECK_LIB(nsl, inet_ntoa, [LIBS="$LIBS -lnsl"])

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_HEADER_STAT
AC_CHECK_HEADERS([alloca.h errno.h limits.h fcntl.h math.h sys/time.h \
                  termio.h wordexp.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_UID_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SIGNAL
AC_STRUCT_ST_BLKSIZE

# Checks for library functions.
AC_CHECK_FUNCS([atexit getcwd memset stat strchr strrchr strstr strtod strtol],
               [], [AC_MSG_ERROR([Required function is not available!])])
AC_FUNC_ALLOCA
AC_FUNC_VFORK
AC_FUNC_MEMCMP
AC_FUNC_STRFTIME

ESO_FUNC_GETOPT
ESO_FUNC_GETPWUID
ESO_FUNC_GETUID
ESO_FUNC_LSTAT
ESO_FUNC_STRDUP
ESO_FUNC_STPCPY
ESO_FUNC_VSNPRINTF
ESO_FUNC_VASPRINTF
ESO_FUNC_SYSCONF
ESO_FUNC_FPATHCONF
ESO_FUNC_SYMLINK
ESO_FUNC_WORDEXP
ESO_FUNC_GETTIMEOFDAY
ESO_FUNC_REALLOC_SANITY


LIBPIL_SET_DID(PRODUCT_DID, [ESO-VLT-DIC.PRO-1.15])
LIBPIL_SET_PATHS

AM_WITH_DMALLOC

AC_CONFIG_SUBDIRS(kazlib)

AC_CONFIG_FILES(Makefile
                doxygen/Doxyfile
                pil/Makefile
                tests/Makefile)
AC_OUTPUT
