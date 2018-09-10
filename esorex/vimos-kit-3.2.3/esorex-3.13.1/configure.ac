# Process this file with autoconf to produce a configure script.

AC_INIT([ESO Recipe Execution Tool], [3.13.1], [cpl-help@eso.org], [esorex])
AC_PREREQ([2.59])
LT_PREREQ([2.2.6])

AC_CONFIG_SRCDIR([src/er_main.cpp])
AC_CONFIG_AUX_DIR([admin])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
LT_CONFIG_LTDL_DIR([libltdl])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])

AM_MAINTAINER_MODE

LT_INIT([disable-static])
AC_SUBST(LIBTOOL_DEPS)

LTDL_INIT([subproject], [convenience])

ESOREX_SET_PREFIX(${CPLDIR:-/usr/local})

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP

AC_PROG_CXX
AC_PROG_CXXCPP

ESO_ENABLE_DEBUG(no)
ESO_ENABLE_STRICT(no)

AC_PROG_CC_C_O
AC_PROG_CXX_C_O

ESO_PROG_PKGCONFIG
ESO_CHECK_DOCTOOLS

# Checks for libraries.
AC_CHECK_LIB(socket, socket, [LIBS="$LIBS -lsocket"])
AC_CHECK_LIB(nsl, inet_ntoa, [LIBS="$LIBS -lnsl"])
ESOREX_LIBAVCALL_OR_LIBFFI

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_DIRENT
AC_HEADER_STAT
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h unistd.h sys/stat.h sys/types.h pwd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_STRUCT_ST_BLKSIZE
AC_SYS_LARGEFILE

# Checks for library functions.
AC_FUNC_CLOSEDIR_VOID
AC_FUNC_STAT
AC_CHECK_FUNCS([getcwd])

# Checks for symbolic link functions
ESO_FUNC_SYMLINK
ESO_FUNC_LSTAT

ESOREX_FUNC_GETOPT

# Check for CPL presence and usability
CPL_CHECK_CEXT
CPL_CHECK_LIBS
ESO_CHECK_EXTRA_LIBS
ESOREX_CHECK_CPL_PARAMETER_DUPLICATE

ESOREX_SET_PATHS

# Determine the shared library suffix on this platform
LT_SYS_MODULE_EXT


AC_CONFIG_FILES(Makefile
                Doxyfile
                src/Makefile
                src/tests/Makefile
                etc/Makefile
                etc/esorex.rc)
AC_OUTPUT
