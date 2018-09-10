# Process this file with autoconf to produce a configure script.

AC_INIT([VIMOS Instrument Pipeline], [3.2.3], [usd-help@eso.org], [vimos])
AC_PREREQ([2.59])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR([admin])
AC_CONFIG_MACRO_DIR([m4macros])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AM_MAINTAINER_MODE([enable])

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

VIMOS_SET_VERSION_INFO([$VERSION])

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_PROG_LN_S

ESO_ENABLE_DEBUG(no)
ESO_ENABLE_STRICT(no)
ESO_PROG_CC_FLAG([std=c99], [CFLAGS="$CFLAGS -std=c99"])

ESO_CHECK_DOCTOOLS

AC_ENABLE_STATIC(no)
AC_ENABLE_SHARED(yes)

AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# Checks for libraries.
AC_CHECK_LIB(m, pow, [LIBS="$LIBS -lm"])
AC_CHECK_LIB(socket, socket, [LIBS="$LIBS -lsocket"])
AC_CHECK_LIB(nsl, inet_ntoa, [LIBS="$LIBS -lnsl"])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])
ESO_FUNC_STRDUP

# Check for CPL presence and usability
CPL_CHECK_CFITSIO
CPL_CHECK_CEXT
CPL_CHECK_LIBS

# Check for GSL 
CPL_CHECK_GSL

# mark hdrl folder as containing a configurable package
AC_CONFIG_SUBDIRS([mosca/hdrl])
# check and define all required variables to build and
# link hdrl external located in the hdrl folder
HDRL_CHECK([mosca/hdrl])

VIMOS_SET_PATHS
VIMOS_CREATE_SYMBOLS

VIMOS_ENABLE_ONLINE

AM_WITH_DMALLOC

#Check for functions used by the PAF creation code:
ESO_FUNC_GETPWUID
ESO_FUNC_GETUID

# Add supplementary preprocessor symbols to config.h.in
AH_VERBATIM([LINE_LENGTH_MAX],
            [
/*
 * The maximum length, in bytes of an input line. It is set to the
 * POSIX2 value.
 */

#ifndef LINE_LENGTH_MAX
#  define LINE_LENGTH_MAX  2048
#endif
            ])

AH_VERBATIM([PATHNAME_MAX],
            [
/*
 * Define the maximum number of characters a filename including the
 * path may have, excluding the string terminator.
 */

#ifndef PATHNAME_MAX
#  define PATHNAME_MAX  4096
#endif
            ])

AC_CONFIG_SUBDIRS(libpil external/libwcs)

VIMOS_SET_SEX_PREFIX

VIMOS_SET_PSFEX_PREFIX

AC_CONFIG_SUBDIRS(external/sextractor external/psfex)

AC_CONFIG_FILES(Makefile
                doxygen/Doxyfile
                mosca/Makefile
                mosca/libmosca/Makefile
                mosca/libmosca/tests/Makefile
                irplib/Makefile
                irplib/tests/Makefile
                config/Makefile
                external/Makefile
                vimos/Makefile
                recipes/Makefile
                tests/Makefile
                regtests/Makefile
                regtests/tests/Makefile
                reflex/Makefile
                reflex/vimos_ifu.xml
                reflex/vimos_mos.xml)
AC_OUTPUT
