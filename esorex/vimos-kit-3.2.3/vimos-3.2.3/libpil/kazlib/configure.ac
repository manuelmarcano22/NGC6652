# Process this file with autoconf to produce a configure script.

AC_INIT([kazlib], [1.20])
AC_PREREQ([2.59])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR([admin])

AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AM_MAINTAINER_MODE

KAZLIB_ENABLE_CONVENIENCE

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O

KAZLIB_PROG_CC_FLAG([Wall], [CFLAGS="$CFLAGS -Wall -fPIC"])
KAZLIB_PROG_CC_FLAG([Wmissing-prototypes], 
                    [CFLAGS="$CFLAGS -Wmissing-prototypes"])
KAZLIB_PROG_CC_FLAG([ansi], [CFLAGS="$CFLAGS -ansi"])
KAZLIB_PROG_CC_FLAG([pedantic], [CFLAGS="$CFLAGS -pedantic"])

AC_ENABLE_STATIC(yes)
AC_ENABLE_SHARED(no)

AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

KAZLIB_CHECK_DOCTOOLS

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h stddef.h])
AC_CHECK_HEADER([pthread.h],[kaz_have_pthread_h=yes],[kaz_have_pthread_h=no])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([strchr vprintf], [],
               [AC_MSG_ERROR([Required function is not available!])])

KAZLIB_SET_EXTRA_DEFS
KAZLIB_ENABLE_EXCEPTIONS(no)

AC_CONFIG_FILES(Makefile)
AC_OUTPUT
