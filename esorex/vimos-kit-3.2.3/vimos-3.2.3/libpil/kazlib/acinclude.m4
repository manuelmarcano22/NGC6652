# KAZLIB_ENABLE_CONVENIENCE
#--------------------------
AC_DEFUN([KAZLIB_ENABLE_CONVENIENCE],
[

    AC_ARG_ENABLE(kazlib-convenience,
                  AC_HELP_STRING([--enable-kazlib-convenience],
                                 [build convenience library (default is no)]),
                  kazlib_enable_convenience=$enableval,
                  kazlib_enable_convenience=no)

    if test x"${kazlib_enable_convenience}" = xno; then
        kazlib_enable_install=yes
    fi

    AM_CONDITIONAL(KAZLIB_CONVENIENCE,
                   test x"${kazlib_enable_convenience}" = xyes)

])

        
# KAZLIB_ENABLE_EXCEPTIONS(except=no)
#------------------------------------
AC_DEFUN([KAZLIB_ENABLE_EXCEPTIONS],
[

    AC_ARG_ENABLE(kazlib-exceptions,
                  AC_HELP_STRING([--enable-kazlib-exceptions],
                                 [build exceptions module (default is no)]),
                  kazlib_enable_exceptions=$enableval,
                  kazlib_enable_exceptions=no)

    AC_MSG_CHECKING([whether exceptions should be compiled])

    if test x"${kazlib_enable_exceptions}" = xno; then
        EXCEPTIONS=""
    else
        EXCEPTIONS="except.lo"
    fi

    AC_SUBST(EXCEPTIONS)

    AC_MSG_RESULT([$kazlib_enable_exceptions])

])

        
# KAZLIB_SET_EXTRA_DEFS
#----------------------
AC_DEFUN([KAZLIB_SET_EXTRA_DEFS],
[
    AC_ARG_ENABLE(obsolete,
                  AC_HELP_STRING([--enable-obsolete],
                                 [Enable obsolete functions (default is no)]),
                  kazlib_use_obsolete=$enableval, kazlib_use_obsolete=no)

    AC_ARG_ENABLE(threads,
                  AC_HELP_STRING([--enable-threads],
                                 [Enable thread support (default is no)]),
                  kazlib_use_threads=$enableval, kazlib_use_threads=no)

    AC_ARG_ENABLE(rcsid,
                  AC_HELP_STRING([--enable-rcsid],
                                 [Compile in the RCS id (default is yes)]),
                  kazlib_use_rcsid=$enableval, kazlib_use_rcsid=yes)

    AC_ARG_ENABLE(ndebug,
                  AC_HELP_STRING([--enable-ndebug],
                                 [Define symbol NDEBUG (default is yes)]),
                  kazlib_use_ndebug=$enableval, kazlib_use_ndebug=yes)


    if test x$kazlib_use_obsolete = xyes; then
        AC_DEFINE(KAZLIB_OBSOLECENT_DEBUG, 1, [Disables obsolete functions])
    fi

    if test x$kazlib_use_threads = xyes; then
        if test x$kazlib_have_pthread_h = xyes; then
            AC_DEFINE(KAZLIB_POSIX_THREADS, 1, 
                      [Compile with threads support. Requires POSIX threads!])
        else
            AC_MSG_WARN([Missing header file pthread.h.])
            AC_MSG_WARN([Threads support was turned off!])
        fi
    fi

    if test x$kazlib_use_rcsid = xyes; then
        AC_DEFINE(KAZLIB_RCSID, 1, [Compile the RCS id into the library])
    fi

    if test x$kazlib_use_ndebug = xyes; then
        AC_DEFINE(NDEBUG, 1, [Suppress assertions in dict.c])
    fi
])

   
# KAZLIB_PROG_CC_FLAG(FLAG, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#--------------------------------------------------------------------
AC_DEFUN([KAZLIB_PROG_CC_FLAG],
[
    AC_REQUIRE([AC_PROG_CC])

    flag=`echo $1 | sed 'y%.=/+-%___p_%'`
    AC_CACHE_CHECK([whether $CC supports -$1],
                   [kazlib_cv_prog_cc_$flag],
                   [
                       eval "kazlib_cv_prog_cc_$flag=no"
                       AC_LANG_PUSH(C)

                       echo 'int main() { return 0; }' >conftest.$ac_ext

                       try_compile="`$CC -$1 -c conftest.$ac_ext 2>&1`"
                       if test -z "$try_compile"; then
                           try_link="`$CC -$1 -o conftest$ac_exeext \
                                    conftest.$ac_ext 2>&1`"
                           if test -z "$try_link"; then
                               eval "kazlib_cv_prog_cc_$flag=yes"
                           fi
                       fi
                       rm -f conftest*

                       AC_LANG_POP(C)
                   ])

    if eval "test \"`echo '$kazlib_cv_prog_cc_'$flag`\" = yes"; then
        :
        $2
    else
        :
        $3
    fi
])


# KAZLIB_CHECK_DOCTOOLS
#----------------------
AC_DEFUN([KAZLIB_CHECK_DOCTOOLS],
[
    AC_ARG_VAR([LATEX], [latex command])
    AC_PATH_PROG([LATEX], [latex])

    AC_ARG_VAR([MAKEINDEX], [makeindex command])
    AC_PATH_PROG([MAKEINDEX], [makeindex])


    if test -z "${LATEX}"; then
        LATEX=":"
    fi

    if test -z "${MAKEINDEX}"; then
        MAKEINDEX=":"
    fi

])
