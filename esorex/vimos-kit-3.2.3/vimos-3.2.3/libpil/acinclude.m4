# LIBPIL_SET_PREFIX(PREFIX)
#---------------------------
AC_DEFUN([LIBPIL_SET_PREFIX],
[
    unset CDPATH
    # make $PIPE_HOME the default for the installation
    AC_PREFIX_DEFAULT($1)

    if test "x$prefix" = "xNONE"; then
        prefix=$ac_default_prefix
        ac_configure_args="$ac_configure_args --prefix $prefix"
    fi

    if test "x$exec_prefix" = "xNONE"; then
        exec_prefix=$prefix
    fi

])


# LIBPIL_SET_VERSION_INFO(VERSION, [CURRENT], [REVISION], [AGE])
#---------------------------------------------------------------
# Setup various version information, especially the libtool versioning
AC_DEFUN([LIBPIL_SET_VERSION_INFO],
[
    libpil_version=`echo "$1" | sed -e 's/[[a-z,A-Z]].*$//'`

    libpil_major_version=`echo "$libpil_version" | \
        sed 's/\([[0-9]]*\).\(.*\)/\1/'`
    libpil_minor_version=`echo "$libpil_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\)\(.*\)/\2/'`
    libpil_micro_version=`echo "$libpil_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

    if test -z "$libpil_major_version"; then
        libpil_major_version=0
    fi

    if test -z "$libpil_minor_version"; then
        libpil_minor_version=0
    fi

    if test -z "$libpil_micro_version"; then
        libpil_micro_version=0
    fi

    LIBPIL_VERSION="$libpil_version"
    LIBPIL_MAJOR_VERSION=$libpil_major_version
    LIBPIL_MINOR_VERSION=$libpil_minor_version
    LIBPIL_MICRO_VERSION=$libpil_micro_version

    if test -z "$4"; then
        LIBPIL_INTERFACE_AGE=0
    else
        LIBPIL_INTERFACE_AGE="$4"
    fi

    LIBPIL_BINARY_AGE=`expr 100 '*' $LIBPIL_MINOR_VERSION + $LIBPIL_MICRO_VERSION`
    LIBPIL_BINARY_VERSION=`expr 10000 '*' $LIBPIL_MAJOR_VERSION + \
                          $LIBPIL_BINARY_AGE`

    AC_SUBST(LIBPIL_VERSION)
    AC_SUBST(LIBPIL_MAJOR_VERSION)
    AC_SUBST(LIBPIL_MINOR_VERSION)
    AC_SUBST(LIBPIL_MICRO_VERSION)
    AC_SUBST(LIBPIL_INTERFACE_AGE)
    AC_SUBST(LIBPIL_BINARY_VERSION)
    AC_SUBST(LIBPIL_BINARY_AGE)

    AC_DEFINE_UNQUOTED(LIBPIL_MAJOR_VERSION, $LIBPIL_MAJOR_VERSION,
                       [LIBPIL major version number])
    AC_DEFINE_UNQUOTED(LIBPIL_MINOR_VERSION, $LIBPIL_MINOR_VERSION,
                       [LIBPIL minor version number])
    AC_DEFINE_UNQUOTED(LIBPIL_MICRO_VERSION, $LIBPIL_MICRO_VERSION,
                       [LIBPIL micro version number])
    AC_DEFINE_UNQUOTED(LIBPIL_INTERFACE_AGE, $LIBPIL_INTERFACE_AGE,
                       [LIBPIL interface age])
    AC_DEFINE_UNQUOTED(LIBPIL_BINARY_VERSION, $LIBPIL_BINARY_VERSION,
                       [LIBPIL binary version number])
    AC_DEFINE_UNQUOTED(LIBPIL_BINARY_AGE, $LIBPIL_BINARY_AGE,
                       [LIBPIL binary age])

    ESO_SET_LIBRARY_VERSION([$2], [$3], [$4])
])


# LIBPIL_SET_DID(SYMBOL, NAME)
#-----------------------------
# Adds a DID identifier string to config.h
AC_DEFUN([LIBPIL_SET_DID],
[

    if test -n "$1" || test -n "$2"; then
        AC_DEFINE_UNQUOTED($1, "$2",
                  [Product DID identifier the library complies to])
        AC_SUBST(PRODUCT_DID)
    fi

])


# LIBPIL_SET_PATHS
#---------------------
# Define auxiliary directories of the installed directory tree.
AC_DEFUN([LIBPIL_SET_PATHS],
[

    if test -z "$pilprivatelibdir"; then
        pilprivatelibdir='${libdir}/${PACKAGE}-${VERSION}'
    fi

    if test x"$includedir" = x'${prefix}/include'; then
        includedir='${prefix}/include/pil'
    fi

    if test -z "$pipedocsdir"; then
        pipedocsdir='${datadir}/doc/esopipes/${PACKAGE}-${VERSION}/'
    fi

    if test -z "$apidocdir"; then
        apidocdir='${pipedocsdir}/html'
    fi

    if test -z "$configdir"; then
       configdir='${datadir}/${PACKAGE}/config'
    fi

    AC_SUBST(configdir)
    AC_SUBST(pilprivatelibdir)
    AC_SUBST(pipedocsdir)
    AC_SUBST(apidocdir)

])
