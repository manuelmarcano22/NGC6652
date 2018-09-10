# KAZLIB_CONVENIENCE_LIBRARY(dir=kazlib)
#---------------------------------------
AC_DEFUN([KAZLIB_CONVENIENCE_LIBRARY],
[

    case $kazlib_enable_convenience in
        no) AC_MSG_ERROR([this package needs a convenience libkaz])
            ;;
        "") kazlib_enable_convenience=yes
            ac_configure_args="$ac_configure_args --enable-kazlib-convenience"
            ;;
    esac

    LIBKAZ='$(top_builddir)/'ifelse($#,1,[$1],['kazlib'])/libkazc.la
    INCKAZ='-I$(top_srcdir)/'ifelse($#,1,[$1],['kazlib'])

])

# KAZLIB_INSTALLABLE_LIBRARY(dir=kazlib)
#---------------------------------------
AC_DEFUN([KAZLIB_INSTALLABLE_LIBRARY],
[

    AC_CHECK_LIB(kaz, dict_create,
                 [
                     test x"${kazlib_enable_install}" != xyes && \
                          kazlib_enable_install=no
                 ],
                 [
                     if test x"${kazlib_enable_install}" = xno; then
                         AC_MSG_WARN([libkaz not installed but installation disabled])
                     else
                         kazlib_enable_install=yes
                     fi
                 ])

    if test x"${kazlib_enable_install}" = xyes; then
        ac_configure_args="$ac_configure_args --enable-kazlib-install"

        LIBKAZ='$(top_builddir)/'ifelse($#, 1, [$1], ['kazlib'])/libkaz.la
        INCKAZ='-I$(top_srcdir)/'ifelse($#, 1, [$1], ['kazlib'])
    else
        ac_configure_args="$ac_configure_args --enable-kazlib-install=no"
        LIBKAZ="-lkaz"
        INCKAZ=
    fi
])
