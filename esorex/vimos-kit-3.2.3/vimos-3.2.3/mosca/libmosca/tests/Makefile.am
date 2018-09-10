## This file is part of the Refran Library
## Copyright (C) 2009 César Enrique García
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
##

#AUTOMAKE_OPTIONS = 1.8 foreign color-tests
AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~ *.testlog *.stdout test-reports/*

EXTRA_DIST = test_driver data

AM_CPPFLAGS = -DDATADIR=\"$(srcdir)/\"
AM_CXXFLAGS = -Wall -Wextra -pedantic \
	$(MOSCA_INCLUDES) $(CPL_INCLUDES) $(HDRL_INCLUDES)

if MAINTAINER_MODE

MAINTAINERCLEANFILES = $(srcdir)/Makefile.in

endif

AM_LDFLAGS = $(HDRL_LDFLAGS) $(CPL_LDFLAGS)
LDADD = $(HDRL_LIBS) $(LIBMOSCA) $(LIBCPLCORE) $(LIBCPLUI) \
	$(BOOST_UNIT_TEST_FRAMEWORK_LIB)

check_PROGRAMS = image_smooth-test flat_combine-test fiera_config-test \
	image_normalisation-test


image_smooth_test_SOURCES = image_smooth-test.cpp

flat_combine_test_SOURCES = flat_combine-test.cpp

fiera_config_test_SOURCES = fiera_config-test.cpp

image_normalisation_test_SOURCES = image_normalisation-test.cpp

LOG_COMPILER = $(srcdir)/test_driver

TESTS = $(check_PROGRAMS)
