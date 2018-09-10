# SEXTRACTOR_PROG_CC
#-------------------
# Check for a C compiler the sextractor way.
AC_DEFUN([SEXTRACTOR_PROG_CC],
[
    AC_REQUIRE([AC_CANONICAL_TARGET])
    AC_REQUIRE([AC_PROG_CC])

    AC_ARG_WITH(machine,
                [  --with-machine=MACHINE  configure sextractor for building on MACHINE],
                sextractor_host=$withval, sextractor_host=auto)

    AC_MSG_CHECKING([for sextractor target host configuration])

    AC_CACHE_VAL(sextractor_cv_host,
    [
        if test x"$sextractor_host" = xauto; then
            case "${target}" in
                *ibm-aix*)    sextractor_host_tag="ibm";;
                *dec-osf*)    sextractor_host_tag="alpha";;
                *dec-ultrix*) sextractor_host_tag="ultrix";;
                *hp-hpux*)    sextractor_host_tag="hpux";;
                *solaris*)    sextractor_host_tag="solaris";;
                *sunos*)      sextractor_host_tag="sunos";;
                *sgi-irix*)   sextractor_host_tag="sgi";;
                *linux-gnu*)  sextractor_host_tag="linuxpc";;
                *darwin*)     sextractor_host_tag="linuxpc";;
                *)
                    AC_MSG_ERROR([Unsupported host type $host])
                    ;;
            esac
        else
            sextractor_host_tag=$sextractor_host
        fi

        sextractor_cv_host="sextractor_host_tag=$sextractor_host_tag"
    ])

    eval "$sextractor_cv_host"
    AC_MSG_RESULT($sextractor_host_tag)

    #Check extra compiler flags and add them to CFLAGS
    case "${sextractor_host_tag}" in
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
                ESO_PROG_CC_FLAG(funroll-loops,
                                 [CFLAGS="$CFLAGS -O1 -funroll-loops"])
            fi
            ;;
        *)
            if test x"$sextractor_host" != xauto; then
                AC_MSG_ERROR([Unsupported host type $sextractor_host])
            fi
            ;;
    esac

    SEXMACHINE="$sextractor_host_tag"
    AC_SUBST(SEXMACHINE)
])
