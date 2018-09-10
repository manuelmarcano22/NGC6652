AC_DEFUN([VIMOS_SET_SEX_PREFIX],
[
    if test "x$exec_prefix" != "xNONE"; then
        sex_exec_prefix=$exec_prefix/lib/${PACKAGE}-${VERSION}/bin
        ac_configure_args="$ac_configure_args SEXBINPATH=$sex_exec_prefix"
    elif test "x$prefix" != "xNONE"; then
        sex_exec_prefix=$prefix/lib/${PACKAGE}-${VERSION}/bin
        ac_configure_args="$ac_configure_args SEXBINPATH=$sex_exec_prefix"
    else
        sex_exec_prefix=$ac_default_prefix/lib/${PACKAGE}-${VERSION}/bin
        ac_configure_args="$ac_configure_args SEXBINPATH=$sex_exec_prefix"
    fi
])

AC_DEFUN([VIMOS_SET_PSFEX_PREFIX],
[
    if test "x$exec_prefix" != "xNONE"; then
        psfex_exec_prefix=$exec_prefix/lib/${PACKAGE}-${VERSION}/bin
        ac_configure_args="$ac_configure_args PSFEXBINPATH=$psfex_exec_prefix"
    elif test "x$prefix" != "xNONE"; then
        psfex_exec_prefix=$prefix/lib/${PACKAGE}-${VERSION}/bin
        ac_configure_args="$ac_configure_args PSFEXBINPATH=$psfex_exec_prefix"
    else
        psfex_exec_prefix=$ac_default_prefix/lib/${PACKAGE}-${VERSION}/bin
        ac_configure_args="$ac_configure_args PSFEXBINPATH=$psfex_exec_prefix"
    fi
])


# VIMOS_SET_VERSION_INFO(VERSION, [CURRENT], [REVISION], [AGE])
#--------------------------------------------------------------
# Setup various version information, especially the libtool versioning
AC_DEFUN([VIMOS_SET_VERSION_INFO],
[
    vimos_version=`echo "$1" | sed -e 's/[[a-z,A-Z]].*$//'`

    vimos_major_version=`echo "$vimos_version" | \
        sed 's/\([[0-9]]*\).\(.*\)/\1/'`
    vimos_minor_version=`echo "$vimos_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\)\(.*\)/\2/'`
    vimos_micro_version=`echo "$vimos_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

    if test -z "$vimos_major_version"; then
        vimos_major_version=0
    fi

    if test -z "$vimos_minor_version"; then
        vimos_minor_version=0
    fi

    if test -z "$vimos_micro_version"; then
        vimos_micro_version=0
    fi

    VIMOS_VERSION="$vimos_version"
    VIMOS_MAJOR_VERSION=$vimos_major_version
    VIMOS_MINOR_VERSION=$vimos_minor_version
    VIMOS_MICRO_VERSION=$vimos_micro_version

    if test -z "$4"; then
        VIMOS_INTERFACE_AGE=0
    else
        VIMOS_INTERFACE_AGE="$4"
    fi

    VIMOS_BINARY_AGE=`expr 100 '*' $VIMOS_MINOR_VERSION + $VIMOS_MICRO_VERSION`
    VIMOS_BINARY_VERSION=`expr 10000 '*' $VIMOS_MAJOR_VERSION + \
                          $VIMOS_BINARY_AGE`

    AC_SUBST(VIMOS_VERSION)
    AC_SUBST(VIMOS_MAJOR_VERSION)
    AC_SUBST(VIMOS_MINOR_VERSION)
    AC_SUBST(VIMOS_MICRO_VERSION)
    AC_SUBST(VIMOS_INTERFACE_AGE)
    AC_SUBST(VIMOS_BINARY_VERSION)
    AC_SUBST(VIMOS_BINARY_AGE)

    AC_DEFINE_UNQUOTED(VIMOS_MAJOR_VERSION, $VIMOS_MAJOR_VERSION,
                       [VIMOS major version number])
    AC_DEFINE_UNQUOTED(VIMOS_MINOR_VERSION, $VIMOS_MINOR_VERSION,
                       [VIMOS minor version number])
    AC_DEFINE_UNQUOTED(VIMOS_MICRO_VERSION, $VIMOS_MICRO_VERSION,
                       [VIMOS micro version number])
    AC_DEFINE_UNQUOTED(VIMOS_INTERFACE_AGE, $VIMOS_INTERFACE_AGE,
                       [VIMOS interface age])
    AC_DEFINE_UNQUOTED(VIMOS_BINARY_VERSION, $VIMOS_BINARY_VERSION,
                       [VIMOS binary version number])
    AC_DEFINE_UNQUOTED(VIMOS_BINARY_AGE, $VIMOS_BINARY_AGE,
                       [VIMOS binary age])

    ESO_SET_LIBRARY_VERSION([$2], [$3], [$4])
])


# VIMOS_SET_PATHS
#----------------
# Define auxiliary directories of the installed directory tree.
AC_DEFUN([VIMOS_SET_PATHS],
[

    if test -z "$plugindir"; then
        plugindir='${libdir}/esopipes-plugins/${PACKAGE}-${VERSION}'
    fi

    if test -z "$privatelibdir"; then
        privatelibdir='${libdir}/${PACKAGE}-${VERSION}'
    fi

    if test -z "$pipedocsdir"; then
        pipedocsdir='${datadir}/doc/esopipes/${PACKAGE}-${VERSION}/'
    fi

    htmldir='${pipedocsdir}/html'

    if test -z "$apidocdir"; then
        apidocdir='${pipedocsdir}/html'
    fi

    if test -z "$configdir"; then
       configdir='${prefix}/share/esopipes/${PACKAGE}-${VERSION}/config'
    fi

    if test -z "$wkfextradir"; then
        wkfextradir='${datadir}/esopipes/${PACKAGE}-${VERSION}/reflex'
    fi

    if test -z "$wkfcopydir"; then
        wkfcopydir='${datadir}/reflex/workflows/${PACKAGE}-${VERSION}'
    fi

    AC_SUBST(plugindir)
    AC_SUBST(privatelibdir)
    AC_SUBST(pipedocsdir)
    AC_SUBST(apidocdir)
    AC_SUBST(configdir)
    AC_SUBST(wkfextradir)
    AC_SUBST(wkfcopydir)


    # Define a preprocesor symbol for the plugin search paths

    AC_DEFINE_UNQUOTED(VIMOS_PLUGIN_DIR, "esopipes-plugins",
                       [Plugin directory tree prefix])

    eval plugin_dir="$plugindir"
    plugin_path=`eval echo $plugin_dir | \
                sed -e "s/\/${PACKAGE}-${VERSION}.*$//"`

    AC_DEFINE_UNQUOTED(VIMOS_PLUGIN_PATH, "$plugin_path",
                       [Absolute path to the plugin directory tree])

    # Define the preprocessor symbols for the sextractor executable
    # and the configuration files.

    eval sext_bindir="${prefix}/lib/${PACKAGE}-${VERSION}/bin"

    AC_DEFINE_UNQUOTED(VIMOS_SEXTRACTOR_PATH, "$sext_bindir",
                       [Absolute path to the sextractor executable])

    AC_SUBST(sext_bindir)

    eval d="$configdir"
    eval sext_configdir="$d"

    AC_DEFINE_UNQUOTED(VIMOS_SEXTRACTOR_CONFIG, "$sext_configdir",
                       [Absolute path to the sextractor configuration files])

])


# VIMOS_CREATE_SYMBOLS
#---------------------
# Define include and library related makefile symbols
AC_DEFUN([VIMOS_CREATE_SYMBOLS],
[

    # Symbols for package include file and library search paths

    VIMOS_INCLUDES='-I$(top_srcdir)/vimos'
    VIMOS_LDFLAGS='-L$(top_builddir)/vimos'
    IRPLIB_INCLUDES='-I$(top_srcdir)/irplib'
    MOSCA_INCLUDES='-I$(top_srcdir)/mosca/libmosca'
    # No -L for IRPLIB which is statically linked

    LIBVIMOSWCS_INCLUDES='-I$(top_srcdir)/external/libwcs'
    LIBVIMOSWCS_LDFLAGS='-L$(top_builddir)/external/libwcs'

    LIBPIL_INCLUDES='-I$(top_srcdir)/libpil/pil -I$(top_srcdir)/libpil/kazlib'
    LIBPIL_LDFLAGS='-L$(top_builddir)/libpil/pil -L$(top_builddir)/libpil/kazlib'

#    CFITSIO_INCLUDES='-I$(top_srcdir)/libpil/cfitsio'
#    CFITSIO_LDFLAGS='-L$(top_srcdir)/libpil/cfitsio'

#    CFITSIO_INCLUDES='-I$(CPLDIR)/include'
#    CFITSIO_LDFLAGS='-L$(CPLDIR)/lib'

    # Library aliases

    LIBVIMOS='$(top_builddir)/vimos/libvimos.la'
    LIBPIL='$(top_builddir)/libpil/pil/libpil.la'
#    LIBCFITSIO='$(top_srcdir)/libpil/cfitsio/libcfitsio.la'
    LIBVIMOSWCS='$(top_builddir)/external/libwcs/libvimoswcs.la'
    LIBIRPLIB='$(top_builddir)/irplib/libirplib.la'
    LIBMOSCA='$(top_builddir)/mosca/libmosca/libmosca.la'


    # Substitute the defined symbols

    AC_SUBST(CFITSIO_INCLUDES)
    AC_SUBST(CFITSIO_LDFLAGS)

    AC_SUBST(LIBVIMOSWCS_INCLUDES)
    AC_SUBST(LIBVIMOSWCS_LDFLAGS)

    AC_SUBST(LIBPIL_INCLUDES)
    AC_SUBST(LIBPIL_LDFLAGS)

    AC_SUBST(VIMOS_INCLUDES)
    AC_SUBST(VIMOS_LDFLAGS)

    AC_SUBST(LIBVIMOS)
    AC_SUBST(LIBPIL)
    AC_SUBST(LIBCFITSIO)
    AC_SUBST(LIBVIMOSWCS)

    AC_SUBST(IRPLIB_INCLUDES)
    AC_SUBST(LIBIRPLIB)

    AC_SUBST(MOSCA_INCLUDES)
    AC_SUBST(LIBMOSCA)

    # Check for CPL and user defined libraries
    AC_REQUIRE([CPL_CHECK_LIBS])
    AC_REQUIRE([ESO_CHECK_EXTRA_LIBS])

    all_includes='$(VIMOS_INCLUDES) $(MOSCA_INCLUDES) $(LIBPIL_INCLUDES) $(LIBVIMOSWCS_INCLUDES) $(CFITSIO_INCLUDES) $(IRPLIB_INCLUDES) $(CPL_INCLUDES) $(CX_INCLUDES) $(EXTRA_INCLUDES)'
    all_ldflags='$(VIMOS_LDFLAGS) $(LIBPIL_LDFLAGS) $(LIBVIMOSWCS_LDFLAGS) $(CFITSIO_LDFLAGS) $(CPL_LDFLAGS) $(CX_LDFLAGS) $(EXTRA_LDFLAGS)'

    AC_SUBST(all_includes)
    AC_SUBST(all_ldflags)

])


# VIMOS_ENABLE_ONLINE
#--------------------
# Enable the building of extra tools for PSO.
AC_DEFUN([VIMOS_ENABLE_ONLINE],
[

    AH_TEMPLATE([ONLINE_MODE],
                [Define if online support tools should be built])

    AC_ARG_ENABLE(online,
                  AC_HELP_STRING([--enable-online],
                                 [enable online support for PSO [[default=yes]]]),
                  vimos_enable_online=$enableval, vimos_enable_online=no)

    AC_CACHE_CHECK([whether an online support should be enabled],
                   vimos_cv_enable_online,
                   vimos_cv_enable_online=$vimos_enable_online)

    if test x"$vimos_cv_enable_online" = xyes; then

        PSFEXDIR=psfex
        PSFEX_CONFIG="masktoccd.psfex masktoccd_1.sex masktoccd_2.sex"

        AC_DEFINE(ONLINE_MODE)

    else

        PSFEXDIR=""
        PSFEX_CONFIG=""

    fi

    AM_CONDITIONAL([ONLINE_MODE], [test x$vimos_cv_enable_online = xyes])

    AC_SUBST(PSFEXDIR)
    AC_SUBST(PSFEX_CONFIG)

])
