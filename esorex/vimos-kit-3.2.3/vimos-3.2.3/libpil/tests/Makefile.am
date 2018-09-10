# Process this file with automake to produce Makefile.in

##   This file is part of the VIMOS Pipeline Library
##   Copyright (C) 2001-2004 European Southern Observatory

##   This library is free software; you can redistribute it and/or modify
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

MAINTAINERCLEANFILES = Makefile.in

endif


AM_CPPFLAGS = -I$(top_srcdir)/pil -I$(top_srcdir)/kazlib
LIBADD = $(top_builddir)/pil/libpil.la

check_PROGRAMS = paf-test

paf_test_SOURCES = paf-test.c
paf_test_LDFLAGS = $(CFITSIO_LDFLAGS)
paf_test_LDADD = $(top_builddir)/pil/libpil.la $(LIBCFITSIO)

# Be sure to reexport important environment variables.
TESTS_ENVIRONMENT = MAKE="$(MAKE)" CC="$(CC)" CFLAGS="$(CFLAGS)" \
        CPPFLAGS="$(CPPFLAGS)" LD="$(LD)" LDFLAGS="$(LDFLAGS)" \
        LIBS="$(LIBS)" LN_S="$(LN_S)" NM="$(NM)" RANLIB="$(RANLIB)" \
        OBJEXT="$(OBJEXT)" EXEEXT="$(EXEEXT)" MALLOC_PERTURB_=231 \
        MALLOC_CHECK_=2

TESTS = paf-test

# We need to remove any files that the above tests created.
clean-local:
	$(RM) piltest.paf
