# ESOREX_SET_PREFIX(PREFIX)
#---------------------------
AC_DEFUN([ESOREX_SET_PREFIX],
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


# ESOREX_SET_PATHS
#------------------
# Define auxiliary directories of the installed directory tree.
AC_DEFUN([ESOREX_SET_PATHS],
[
    AC_REQUIRE(ESO_PROG_PKGCONFIG)

    if test -z "$apidocdir"; then
        apidocdir='${datadir}/doc/${PACKAGE}/html'
    fi
    AC_SUBST(apidocdir)

    if test -z "$configdir"; then
       configdir='${datadir}/${PACKAGE}/config'
    fi
    AC_SUBST(configdir)

    if test -z "$recipedir"; then
        recipedir="`eval echo $libdir`"
    fi
    AC_SUBST(recipedir)

    if test -z "$bashcompdir"; then
        if test -n "$PKGCONFIG" && \
          pkg-config --exists bash-completion && \
          pkg-config --define-variable=prefix=${prefix} \
            --variable=completionsdir bash-completion 2>&1 > /dev/null; then
        bashcompdir="$(pkg-config --define-variable=prefix=${prefix} --variable=completionsdir bash-completion)"
        else
            bashcompdir="${datarootdir}/bash-completion/completions"
        fi
    fi
    AC_SUBST(bashcompdir)

])


# ESOREX_FUNC_GETOPT
#--------------------
# Checks for GNU getopt_long declaration and function.
AC_DEFUN([ESOREX_FUNC_GETOPT],
[

    AH_TEMPLATE([HAVE_GETOPT_LONG_ONLY],
                [Define if you have the `getopt_long_only' function])

    ESO_CHECK_FUNC(getopt_long_only, [#include <getopt.h>],
                   HAVE_GETOPT_LONG_ONLY)

    if test x"$ac_cv_func_getopt_long_only" = xno ||
       test x"$eso_cv_have_decl_getopt_long_only" = xno; then
        if test -n "$LIBTOOL"; then
            GETOPT="getopt.lo getopt1.lo"
        else
            GETOPT="getopt.$ac_objext getopt1.$ac_objext"
        fi
    fi

    AC_SUBST(GETOPT)

])


# ESOREX_AVCALL
#--------------------
# Checks for the libavcall.a static library.
AC_DEFUN([ESOREX_AVCALL],
[

    AC_ARG_WITH(avcall-include,
                AC_HELP_STRING([--with-avcall-include],
                               [location of the avcall.h header file]),
                esorex_with_avcall_include=$withval)

    AC_ARG_WITH(avcall-lib,
                AC_HELP_STRING([--with-avcall-lib],
                               [full path to the libavcall.a library]),
                esorex_with_avcall_lib=$withval)

    AH_TEMPLATE([HAVE_LIBAVCALL],
                [Define if libavcall.a was found])

    AC_MSG_CHECKING([for libavcall])

    # Directories given as arguments replace a standard system installation
    # setup if they are given.

    esorex_avcall_cflags=""
    if test -n "$esorex_with_avcall_include"; then
        esorex_avcall_cflags="-I$esorex_with_avcall_include"
    else
        for N in /usr/include /opt/local/include ; do
            if test -f "$N/avcall.h" ; then
                esorex_avcall_cflags="-I$N"
            fi
        done
    fi

    esorex_avcall_lib=""
    if test -n "$esorex_with_avcall_lib"; then
        esorex_avcall_lib="$esorex_with_avcall_lib"
    else
        for N in /lib /lib64 /usr/lib /usr/lib64 /opt/local/lib ; do
            if test -f "$N/libavcall.a" ; then
                esorex_avcall_lib="$N/libavcall.a"
            fi
        done
    fi

    # Check whether the header file and the library are present and whether
    # they can be used.

    esorex_have_avcall_header="no"
    esorex_have_avcall_library="no"

    AC_LANG_PUSH(C)

    esorex_avcall_cflags_save="$CFLAGS"
    esorex_avcall_ldflags_save="$LDFLAGS"
    esorex_avcall_libs_save="$LIBS"

    CFLAGS="$esorex_avcall_cflags"
    LDFLAGS=""
    LIBS="$esorex_avcall_lib"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                       [[
                       #include <avcall.h>
                       ]],
                       [
                       av_alist alist;
                       ])],
                      [
                      esorex_have_avcall_header="yes"
                      ],
                      [
                      esorex_have_avcall_header="no"
                      ])

    AC_LINK_IFELSE([AC_LANG_PROGRAM(
                    [[
                    #include <avcall.h>
                    int func(int x) { return x+3; }
                    ]],
                    [[
                    int result;
                    av_alist alist;
                    av_start_int(alist, &func, &result);
                    av_int(alist, (int)2);
                    av_call(alist);
                    ]])],
                   [
                   esorex_have_avcall_library="yes"
                   ],
                   [
                   esorex_have_avcall_library="no"
                   ])

    AC_LANG_POP(C)

    CFLAGS="$esorex_avcall_cflags_save"
    LDFLAGS="$esorex_avcall_ldflags_save"
    LIBS="$esorex_avcall_libs_save"

    if test x"$esorex_have_avcall_library" = xno || \
       test x"$esorex_have_avcall_header" = xno; then

        esorex_found_libavcall=no
        AC_MSG_RESULT([no])

        AVCALL_INCLUDES=""
        LIBAVCALL=""

    else

        esorex_found_libavcall=yes
        AC_MSG_RESULT([yes])

        AC_DEFINE(HAVE_LIBAVCALL)
        AVCALL_INCLUDES="$esorex_avcall_cflags"
        LIBAVCALL="$esorex_avcall_lib"

    fi

    AC_SUBST(AVCALL_INCLUDES)
    AC_SUBST(LIBAVCALL)
    AM_CONDITIONAL([HAVE_LIBAVCALL], [test x$esorex_found_libavcall = xyes])

])


# ESOREX_LIBFFI
#--------------------
# Checks for the libffi.so library.
AC_DEFUN([ESOREX_LIBFFI],
[

    AC_ARG_WITH(libffi,
                AC_HELP_STRING([--with-libffi],
                               [location of where libffi is installed]),
                esorex_with_libffi=$withval)

    AC_ARG_WITH(libffi-includes,
                AC_HELP_STRING([--with-libffi-includes],
                               [location of the libffi header files]),
                esorex_with_libffi_includes=$withval)

    AC_ARG_WITH(libffi-libs,
                AC_HELP_STRING([--with-libffi-libs],
                               [location of the libffi libraries]),
                esorex_with_libffi_libs=$withval)

    AH_TEMPLATE([HAVE_LIBFFI],
                [Define if libffi.so was found])

    AC_MSG_CHECKING([for libffi])

    # Directories given as arguments replace a standard system installation
    # setup if they are given.

    esorex_libffi_cflags=""
    esorex_libffi_ldflags=""
    esorex_libffi_libs="-lffi"

    if test -n "$esorex_with_libffi"; then
        esorex_libffi_cflags="-I$esorex_with_libffi/include"
        esorex_libffi_ldflags="-L$esorex_with_libffi/lib"
        if test -e "$esorex_with_libffi/lib64"; then
            esorex_libffi_ldflags="$esorex_libffi_ldflags -L$esorex_with_libffi/lib64"
        fi
    fi

    if test -n "$esorex_with_libffi_includes"; then
        esorex_libffi_cflags="-I$esorex_with_libffi_includes"
    fi

    if test -n "$esorex_with_libffi_libs"; then
        esorex_libffi_ldflags="-L$esorex_with_libffi_libs"
    fi

    # Check whether the header file and the library are present and whether
    # they can be used.

    esorex_have_libffi_header="no"
    esorex_have_libffi_library="no"

    AC_LANG_PUSH(C)

    esorex_libffi_cflags_save="$CFLAGS"
    esorex_libffi_ldflags_save="$LDFLAGS"
    esorex_libffi_libs_save="$LIBS"

    CFLAGS="$esorex_libffi_cflags"
    LDFLAGS="$esorex_libffi_ldflags"
    LIBS="$esorex_libffi_libs"

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                       [[
                       #include <ffi.h>
                       ]],
                       [
                       ffi_cif cif;
                       ])],
                      [
                      esorex_have_libffi_header="yes"
                      ],
                      [
                      esorex_have_libffi_header="no"
                      ])

    AC_LINK_IFELSE([AC_LANG_PROGRAM(
                    [[
                    #include <ffi.h>
                    int func(int x) { return x+3; }
                    ]],
                    [[
                    int result;
                    int value = 2;
                    ffi_cif cif;
                    ffi_type * argtypes[1];
                    void * values[1];
                    argtypes[0] = &ffi_type_sint;
                    values[0] = &value;
                    ffi_prep_cif(&cif, FFI_DEFAULT_ABI, 1,
                                 &ffi_type_sint, argtypes);
                    ffi_call(&cif, func, &result, values);
                    ]])],
                   [
                   esorex_have_libffi_library="yes"
                   ],
                   [
                   esorex_have_libffi_library="no"
                   ])

    AC_LANG_POP(C)

    CFLAGS="$esorex_libffi_cflags_save"
    LDFLAGS="$esorex_libffi_ldflags_save"
    LIBS="$esorex_libffi_libs_save"

    if test x"$esorex_have_libffi_library" = xno || \
       test x"$esorex_have_libffi_header" = xno; then

        esorex_found_libffi=no
        AC_MSG_RESULT([no])

        FFI_INCLUDES=""
        FFI_LDFLAGS=""
        LIBFFI=""

    else

        esorex_found_libffi=yes
        AC_MSG_RESULT([yes])

        AC_DEFINE(HAVE_LIBFFI)
        FFI_INCLUDES="$esorex_libffi_cflags"
        FFI_LDFLAGS="$esorex_libffi_ldflags"
        LIBFFI="$esorex_libffi_libs"

    fi

    AC_SUBST(FFI_INCLUDES)
    AC_SUBST(FFI_LDFLAGS)
    AC_SUBST(LIBFFI)
    AM_CONDITIONAL([HAVE_LIBFFI], [test x$esorex_found_libffi = xyes])

])


# ESOREX_LIBAVCALL_OR_LIBFFI
#-----------------------------
# Checks for the libavcall.a or libffi.so libraries.
AC_DEFUN([ESOREX_LIBAVCALL_OR_LIBFFI],
[
    AH_TEMPLATE([ENABLE_PYTHON_RECIPES],
                [Define if Python based recipes should be supported.])

    AC_ARG_ENABLE(python-recipes,
                  AC_HELP_STRING([--enable-python-recipes],
                                 [enables Python based recipes [default=yes]]),
                  esorex_enable_pyrecipes=$enableval,
                  esorex_enable_pyrecipes=default)

    if test x"$esorex_enable_pyrecipes" = xdefault; then
        esorex_enable_pyrecipes=yes
        esorex_warn_disabling_pyrecipes=no
    elif test x"$esorex_enable_pyrecipes" = xyes; then
        esorex_warn_disabling_pyrecipes=yes
    else
        esorex_warn_disabling_pyrecipes=no
    fi

    AC_CACHE_CHECK([whether Python based recipes can be supported],
                   esorex_cv_enable_pyrecipes,
                   [if test x"$esorex_enable_pyrecipes" = xyes; then
                        AC_REQUIRE([ESOREX_AVCALL])
                        AC_REQUIRE([ESOREX_LIBFFI])
                        if test x"$esorex_found_libavcall" = xno && \
                           test x"$esorex_found_libffi" = xno ; then

                            esorex_cv_enable_pyrecipes=no
                        else
                            esorex_cv_enable_pyrecipes=yes
                        fi
                    else
                        esorex_cv_enable_pyrecipes=no
                    fi
                   ])

    if test x"$esorex_warn_disabling_pyrecipes" = xyes && \
       test x"$esorex_cv_enable_pyrecipes" = xno; then

        AC_MSG_WARN([Did not find libavcall or libffi. Disabling support for \
Python based recipe plugins.])

    fi

    if test x"$esorex_cv_enable_pyrecipes" = xyes; then
        AC_DEFINE(ENABLE_PYTHON_RECIPES)
    fi
    AM_CONDITIONAL([ENABLE_PYTHON_RECIPES],
                   [test x"$esorex_cv_enable_pyrecipes" = xyes])
])


# ESOREX_CHECK_CPL_PARAMETER_DUPLICATE
#--------------------------------------
# Checks to see if cpl_parameter_duplicate is broken or not.
# If it is broken then HAVE_BROKEN_CPL_PARAMETER_DUPLICATE will be defined.
#
AC_DEFUN([ESOREX_CHECK_CPL_PARAMETER_DUPLICATE],
[
    AC_REQUIRE([CPL_CHECK_LIBS])

    AC_MSG_CHECKING([if cpl_parameter_duplicate needs workaround])

    AH_TEMPLATE([HAVE_BROKEN_CPL_PARAMETER_DUPLICATE],
                [Define if cpl_parameter_duplicate is broken for ranges])

    AC_LANG_PUSH(C)

    esorex_check_cpl_parameter_duplicate_cflags_save="$CFLAGS"
    esorex_check_cpl_parameter_duplicate_ldflags_save="$LDFLAGS"
    esorex_check_cpl_parameter_duplicate_libs_save="$LIBS"

    CFLAGS="$CPL_CFLAGS"
    LDFLAGS="$CPL_LDFLAGS"
    LIBS="$LIBCPLCORE $LIBCPLUI"

    AC_RUN_IFELSE([AC_LANG_PROGRAM(
                   [[
                   #include <cpl.h>
                   ]],
                   [[
                    cpl_parameter * a = cpl_parameter_new_range(
                            "a", CPL_TYPE_INT, "test", "test", 3, 1, 5);
                    cpl_parameter * b = cpl_parameter_duplicate(a);
                    if (cpl_parameter_get_range_min_int(a)
                        != cpl_parameter_get_range_min_int(b)) return 1;
                    if (cpl_parameter_get_range_max_int(a)
                        != cpl_parameter_get_range_max_int(b)) return 1;
                    return 0;
                   ]])
                  ],
                  [esorex_have_broken_cpl_parameter_duplicate=no],
                  [esorex_have_broken_cpl_parameter_duplicate=yes],
                  [esorex_have_broken_cpl_parameter_duplicate=unknown])

    AC_LANG_POP(C)

    CFLAGS="$esorex_check_cpl_parameter_duplicate_cflags_save"
    LDFLAGS="$esorex_check_cpl_parameter_duplicate_ldflags_save"
    LIBS="$esorex_check_cpl_parameter_duplicate_libs_save"

    if test x$esorex_have_broken_cpl_parameter_duplicate = xyes; then
        AC_DEFINE(HAVE_BROKEN_CPL_PARAMETER_DUPLICATE)
        AC_MSG_RESULT([yes])
    else
        if test x$esorex_have_broken_cpl_parameter_duplicate = xunknown; then
            AC_MSG_RESULT([unknown])
            AC_MSG_WARN([Cannot check if cpl_parameter_duplicate needs a \
workaround. Will assume it does not.])
        else
            AC_MSG_RESULT([no])
        fi
    fi
])

# vim et:sw=4:ts=4
