# ESO_PROG_CC_FLAG(FLAG, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#-----------------------------------------------------------------
AC_DEFUN([ESO_PROG_CC_FLAG],
[
    AC_REQUIRE([AC_PROG_CC])

    flag=`echo $1 | sed 'y%.=/+-%___p_%'`
    AC_MSG_CHECKING([whether $CC supports -$1])
    AC_CACHE_VAL([eso_cv_prog_cc_$flag],
                 [
                     eval "eso_cv_prog_cc_$flag=no"
                     AC_LANG_SAVE
                     AC_LANG_C
		 
                     echo 'int main() { return 0; }' >conftest.$ac_ext
		 
                     try_compile="`$CC -$1 -c conftest.$ac_ext 2>&1`"
                     if test -z "$try_compile"; then
                         try_link="`$CC -$1 -o conftest$ac_exeext \
                                  conftest.$ac_ext 2>&1`"
                         if test -z "$try_link"; then
                             eval "eso_cv_prog_cc_$flag=yes"
                         fi
                     fi
                     rm -f conftest*
		 
                     AC_LANG_RESTORE
                 ])

    if eval "test \"`echo '$eso_cv_prog_cc_'$flag`\" = yes"; then
        AC_MSG_RESULT(yes)
        $2
    else
        AC_MSG_RESULT(no)
        $3
    fi
])


# ESO_PROG_AR
#------------
# Checks if ar is in the path
AC_DEFUN([ESO_PROG_AR],
[
    AC_CHECK_PROG(AR, ar, ar, NONE)

    if test x"$AR" = xNONE; then
        AC_MSG_ERROR([Cannot find \'ar\'])
    fi

])
