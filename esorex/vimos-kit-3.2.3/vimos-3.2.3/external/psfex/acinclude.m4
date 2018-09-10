# PSFEX_PROG_CC
#--------------
# Check for a C compiler the psfex way.
AC_DEFUN([PSFEX_PROG_CC],
[
    AC_REQUIRE([AC_CANONICAL_TARGET])
    AC_REQUIRE([AC_PROG_CC])

    AC_ARG_WITH(machine,
                [  --with-machine=MACHINE  configure psfex for building on MACHINE],
                psfex_host=$withval, psfex_host=auto)

    AC_MSG_CHECKING([for psfex target host configuration])

    AC_CACHE_VAL(psfex_cv_host,
    [
        if test x"$psfex_host" = xauto; then
            case "${target}" in
                *ibm-aix*)    psfex_host_tag="ibm";;
                *dec-osf*)    psfex_host_tag="alpha";;
                *dec-ultrix*) psfex_host_tag="ultrix";;
                *hp-hpux*)    psfex_host_tag="hpux";;
                *solaris*)    psfex_host_tag="solaris";;
                *sunos*)      psfex_host_tag="sunos";;
                *sgi-irix*)   psfex_host_tag="sgi";;
                *linux-gnu*)  psfex_host_tag="linuxpc";;
                *darwin*)     psfex_host_tag="linuxpc";;
                *)
                    AC_MSG_ERROR([Unsupported host type $host])
                    ;;
            esac
        else
            psfex_host_tag=$psfex_host
        fi

        psfex_cv_host="psfex_host_tag=$psfex_host_tag"
    ])

    eval "$psfex_cv_host"
    AC_MSG_RESULT($psfex_host_tag)

    #Check extra compiler flags and add them to CFLAGS
    case "${psfex_host_tag}" in
        ibm)
            CPPFLAGS="$CPPFLAGS -DIBM_AIX"
            if test x"$GCC" != xyes; then
                ESO_PROG_CC_FLAG(O, [CFLAGS="$CFLAGS -O"])
            fi
            ;;
        alpha)
            CPPFLAGS="$CPPFLAGS -DDEC_ALPHA"
            if test x"$GCC" != xyes; then
                ESO_PROG_CC_FLAG(O, [CFLAGS="$CFLAGS -O"])
            fi
            ;;
        ultrix)
            CPPFLAGS="$CPPFLAGS -DBSWAP"
            if test x"$GCC" != xyes; then
                ESO_PROG_CC_FLAG(O, [CFLAGS="$CFLAGS -O"])
            fi
            ;;
        hpux)
            if test x"$GCC" != xyes; then
                CPPFLAGS="$CPPFLAGS -DHP_UX"
            fi
            ;;
        solaris)
            CPPFLAGS="$CPPFLAGS -DSUN_SOLARIS"
            if test x"$GCC" != xyes; then
                ESO_PROG_CC_FLAG(fast, [CFLAGS="$CFLAGS -fast"])
                ESO_PROG_CC_FLAG(xO5, [CFLAGS="$CFLAGS -xO5"])
            fi
            ;;
        sunos)
            CPPFLAGS="$CPPFLAGS -DSUN_OS"
            ;;
        sgi)
            ;;
        linuxpc)
            CPPFLAGS="$CPPFLAGS -DPC_LINUX"
            if test x"$GCC" = xyes; then
                ESO_PROG_CC_FLAG(malign-double,
                                 [CFLAGS="$CFLAGS -malign-double"])
                ESO_PROG_CC_FLAG(finline-functions,
                                 [CFLAGS="$CFLAGS -finline-functions"])
                ESO_PROG_CC_FLAG(funroll-loops,
                                 [CFLAGS="$CFLAGS -funroll-loops"])
            fi
            ;;
        *)
            if test x"$psfex_host" != xauto; then
                AC_MSG_ERROR([Unsupported host type $psfex_host])
            fi
            ;;
    esac

    SEXMACHINE="$psfex_host_tag"
    AC_SUBST(SEXMACHINE)
])
