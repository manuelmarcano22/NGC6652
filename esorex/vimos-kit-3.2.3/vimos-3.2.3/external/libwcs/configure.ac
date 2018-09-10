# Process this file with autoconf to produce a configure script.

AC_INIT([libvimoswcs], [2.9.3])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR(admin)

AM_INIT_AUTOMAKE

AM_MAINTAINER_MODE

# Checks for programs.
AC_PROG_CPP
AC_PROG_CC
AM_PROG_CC_C_O

ESO_PROG_CC_FLAG(fno-builtin, [CFLAGS="$CFLAGS -fno-builtin"])

AC_ENABLE_STATIC(yes)
AC_ENABLE_SHARED(yes)
AC_PROG_LIBTOOL

AC_CONFIG_FILES(Makefile)
AC_OUTPUT
