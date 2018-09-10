## Process this file with automake to produce Makefile.in

##   This file is part of the ESO Common Pipeline Library
##   Copyright (C) 2001-2008 European Southern Observatory

##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.

##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = $(srcdir)/Makefile.in

endif

AM_CPPFLAGS = $(all_includes) -DCX_LOG_DOMAIN=\"CasuLib\" -I$(top_srcdir) -DHDRL_USE_PRIVATE
AM_LDFLAGS = $(CPL_LDFLAGS)
LDADD = ../libcatalogue.la  $(top_builddir)/libhdrl.la $(LIBCPLCORE) $(LIBCPLDRS) $(LIBCEXT) -lm

check_PROGRAMS = imcore_conf-test imcore_background-test \
	imcore_apinit-test imcore_classify-test imcore_radii-test \
	imcore_moments-test imcore_areals-test imcore_seeing-test \
	imcore_polynm-test casu_imcore-test hdrl_casubkg_badpatch-test \
	hdrl_casuaddmul-test hdrl_simulerror_montecarlo-test

imcore_conf_test_SOURCES = imcore_conf-test.c

imcore_background_test_SOURCES = imcore_background-test.c

imcore_apinit_test_SOURCES = imcore_apinit-test.c

imcore_classify_test_SOURCES = imcore_classify-test.c

imcore_radii_test_SOURCES = imcore_radii-test.c

imcore_moments_test_SOURCES = imcore_moments-test.c

imcore_areals_test_SOURCES = imcore_areals-test.c

imcore_seeing_test_SOURCES = imcore_seeing-test.c

imcore_polynm_test_SOURCES = imcore_polynm-test.c

casu_imcore_test_SOURCES = casu_imcore-test.c

hdrl_casubkg_badpatch_test_SOURCES = hdrl_casubkg_badpatch-test.c

hdrl_casuaddmul_test_SOURCES = hdrl_casuaddmul-test.c

hdrl_simulerror_montecarlo_test_SOURCES = hdrl_simulerror_montecarlo-test.c

TESTS = imcore_conf-test imcore_background-test \
	imcore_apinit-test imcore_classify-test imcore_radii-test \
	imcore_moments-test imcore_areals-test imcore_seeing-test \
	imcore_polynm-test casu_imcore-test hdrl_casubkg_badpatch-test \
	hdrl_casuaddmul-test hdrl_simulerror_montecarlo-test

# Be sure to reexport important environment variables.
TESTS_ENVIRONMENT = MAKE="$(MAKE)" CC="$(CC)" CFLAGS="$(CFLAGS)" \
	CPPFLAGS="$(CPPFLAGS)" LD="$(LD)" LDFLAGS="$(LDFLAGS)" \
	LIBS="$(LIBS)" LN_S="$(LN_S)" NM="$(NM)" RANLIB="$(RANLIB)" \
	OBJEXT="$(OBJEXT)" EXEEXT="$(EXEEXT)" \
	MALLOC_PERTURB_=31 MALLOC_CHECK_=2 \
	GOMP_SPINCOUNT=10000

XFAIL_TESTS =

# We need to remove any files that the above tests created.
clean-local:
	$(RM) *.fits *.log
