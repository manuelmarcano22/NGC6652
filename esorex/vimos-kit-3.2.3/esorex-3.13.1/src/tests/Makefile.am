## Process this file with automake to produce Makefile.in

##  This file is part of the ESO Recipe Execution Tool
##  Copyright (C) 2002-2017 European Southern Observatory

##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.

##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.

##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = $(srcdir)/Makefile.in

endif


if HAVE_LIBAVCALL
AM_CPPFLAGS = -I$(top_srcdir)/src $(CPL_INCLUDES) $(CX_INCLUDES)
AM_LDFLAGS = $(CPL_LDFLAGS)
LDADD = $(LIBESOREX) $(LIBAVCALL) $(LIBCPLUI) $(LIBCPLCORE) $(LIBCEXT)
else
if HAVE_LIBFFI
AM_CPPFLAGS = -I$(top_srcdir)/src $(CPL_INCLUDES) $(CX_INCLUDES)
AM_LDFLAGS = $(CPL_LDFLAGS) $(FFI_LDFLAGS)
LDADD = $(LIBESOREX) $(LIBFFI) $(LIBCPLUI) $(LIBCPLCORE) $(LIBCEXT)
else
AM_CPPFLAGS = -I$(top_srcdir)/src $(CPL_INCLUDES) $(CX_INCLUDES)
AM_LDFLAGS = $(CPL_LDFLAGS)
LDADD = $(LIBESOREX) $(LIBCPLUI) $(LIBCPLCORE) $(LIBCEXT)
endif
endif

LIBESOREX = $(top_builddir)/src/libesorex.la

check_PROGRAMS = esorex_json-test \
    esorex_json_internal-test \
    esorex_json_parser-test \
    esorex_json_parse_errors1-test \
    esorex_json_parse_errors2-test \
    esorex_json_parse_errors3-test \
    esorex_json_parse_errors4-test \
    esorex_python-test \
    esorex_python_internal-test \
    esorex_python_errors1-test \
    esorex_python_errors2-test \
    esorex_python_errors3-test \
    esorex_python_errors4-test \
    esorex_python_recipe-test

EXTRA_DIST = testrecipe.py
noinst_HEADERS = esorex_test.h

esorex_json_test_SOURCES = esorex_json-test.c
esorex_json_test_LDFLAGS = $(AM_LDFLAGS)
esorex_json_test_LDADD = $(LDADD)

esorex_json_internal_test_SOURCES = esorex_json_internal-test.c
esorex_json_internal_test_LDFLAGS = $(AM_LDFLAGS)
esorex_json_internal_test_LDADD = $(LDADD)

esorex_json_parser_test_SOURCES = esorex_json_parser-test.c
esorex_json_parser_test_LDFLAGS = $(AM_LDFLAGS)
esorex_json_parser_test_LDADD = $(LDADD)

esorex_json_parse_errors1_test_SOURCES = esorex_json_parse_errors1-test.c
esorex_json_parse_errors1_test_LDFLAGS = $(AM_LDFLAGS)
esorex_json_parse_errors1_test_LDADD = $(LDADD)

esorex_json_parse_errors2_test_SOURCES = esorex_json_parse_errors2-test.c
esorex_json_parse_errors2_test_LDFLAGS = $(AM_LDFLAGS)
esorex_json_parse_errors2_test_LDADD = $(LDADD)

esorex_json_parse_errors3_test_SOURCES = esorex_json_parse_errors3-test.c
esorex_json_parse_errors3_test_LDFLAGS = $(AM_LDFLAGS)
esorex_json_parse_errors3_test_LDADD = $(LDADD)

esorex_json_parse_errors4_test_SOURCES = esorex_json_parse_errors4-test.c
esorex_json_parse_errors4_test_LDFLAGS = $(AM_LDFLAGS)
esorex_json_parse_errors4_test_LDADD = $(LDADD)

esorex_python_test_SOURCES = esorex_python-test.c
esorex_python_test_LDFLAGS = $(AM_LDFLAGS)
esorex_python_test_LDADD = $(LDADD)

esorex_python_internal_test_SOURCES = esorex_python_internal-test.c
esorex_python_internal_test_LDFLAGS = $(AM_LDFLAGS)
esorex_python_internal_test_LDADD = $(LDADD)

esorex_python_errors1_test_SOURCES = esorex_python_errors1-test.c
esorex_python_errors1_test_LDFLAGS = $(AM_LDFLAGS)
esorex_python_errors1_test_LDADD = $(LDADD)

esorex_python_errors2_test_SOURCES = esorex_python_errors2-test.c
esorex_python_errors2_test_LDFLAGS = $(AM_LDFLAGS)
esorex_python_errors2_test_LDADD = $(LDADD)

esorex_python_errors3_test_SOURCES = esorex_python_errors3-test.c
esorex_python_errors3_test_LDFLAGS = $(AM_LDFLAGS)
esorex_python_errors3_test_LDADD = $(LDADD)

esorex_python_errors4_test_SOURCES = esorex_python_errors4-test.c
esorex_python_errors4_test_LDFLAGS = $(AM_LDFLAGS)
esorex_python_errors4_test_LDADD = $(LDADD)

esorex_python_recipe_test_SOURCES = esorex_python_recipe-test.c
esorex_python_recipe_test_CPPFLAGS = $(AM_CPPFLAGS) \
    -DTEST_RECIPE_DIR=\"$(abs_srcdir)\"
esorex_python_recipe_test_LDFLAGS = $(AM_LDFLAGS)
esorex_python_recipe_test_LDADD = $(LDADD)


# Be sure to reexport important environment variables.
TESTS_ENVIRONMENT = MAKE="$(MAKE)" CC="$(CC)" CFLAGS="$(CFLAGS)" \
	CPPFLAGS="$(CPPFLAGS)" LD="$(LD)" LDFLAGS="$(LDFLAGS)" \
	LIBS="$(LIBS)" LN_S="$(LN_S)" NM="$(NM)" RANLIB="$(RANLIB)" \
	OBJEXT="$(OBJEXT)" EXEEXT="$(EXEEXT)"


TESTS = $(check_PROGRAMS)

XFAIL_TESTS =

# We need to remove any files that the above tests created.
clean-local:
	$(RM) *.fits *.json *.log *.pyc *.sof
	$(RM) -r mock_python_*/
