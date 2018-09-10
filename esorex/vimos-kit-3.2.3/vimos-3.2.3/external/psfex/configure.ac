# Process this file with autoconf to produce a configure script.

AC_INIT([psfex],[1.8.1])
AC_PREREQ([2.59])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR([admin])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE

AM_MAINTAINER_MODE

# Checks for programs.
PSFEX_PROG_CC

ESO_PROG_CC_FLAG(fno-builtin, [CFLAGS="$CFLAGS -fno-builtin"])

ESO_PROG_AR
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB(m, pow, [LIBM="-lm"])
AC_SUBST(LIBM)

#Set the path to where the psfex binary will be installed:
AC_ARG_VAR([SEXBINPATH], [path to install psfex binary])
eval psfex_bindir="${prefix}/bin"
if test -z "$PSFEXBINPATH"; then
    eval psfex_bindir="${PSFEXBINPATH}"
fi
eval psfex_bindir="${PSFEXBINPATH}"
AC_SUBST(psfex_bindir)


AC_CONFIG_FILES(Makefile
                source/Makefile)
AC_OUTPUT
