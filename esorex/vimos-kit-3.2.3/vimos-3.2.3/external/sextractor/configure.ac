# Process this file with autoconf to produce a configure script.

AC_INIT([sextractor],[2.1.6])
AC_PREREQ([2.59])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR([admin])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE

AM_MAINTAINER_MODE

# Checks for programs.
SEXTRACTOR_PROG_CC

ESO_PROG_CC_FLAG(fno-builtin, [CFLAGS="$CFLAGS -fno-builtin"])

ESO_PROG_AR
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB(m, pow, [LIBM="-lm"])
AC_SUBST(LIBM)

#Set the path to where the sextractor binary will be installed:
AC_ARG_VAR([SEXBINPATH], [path to install sextractor binary])
eval sext_bindir="${prefix}/bin"
if test -z "$SEXBINPATH"; then
    eval sext_bindir="${SEXBINPATH}"
fi
eval sext_bindir="${SEXBINPATH}"
AC_SUBST(sext_bindir)

AC_CONFIG_FILES(Makefile
                source/Makefile
                source/wcs/Makefile)
AC_OUTPUT
