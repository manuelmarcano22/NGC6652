## Process this file with automake to produce Makefile.in

##   This file is part of the HDRL Pipeline
##   Copyright (C) 2012,2013 European Southern Observatory

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
##   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

SUBDIRS = catalogue . catalogue/tests tests

EXTRA_DIST =  m4macros/eso.m4 m4macros/cpl-pipeline.m4

if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif

ACLOCAL_AMFLAGS = -I m4macros

AM_CPPFLAGS = -DCX_LOG_DOMAIN=\"HdrlLib\" -DHDRL_USE_PRIVATE $(CPL_INCLUDES)

hdrlHEADERS = hdrl_overscan.h \
		 hdrl_sigclip.h \
		 hdrl_utils.h \
		 hdrl_combine.h \
		 hdrl_bpm_2d.h \
		 hdrl_bpm_3d.h \
		 hdrl_bpm_fit.h \
		 hdrl_bpm_utils.h \
		 hdrl_iter.h \
		 hdrl_frameiter.h \
		 hdrl_multiiter.h \
		 hdrl_collapse.h \
		 hdrl_elemop.h \
		 hdrl_flat.h \
		 hdrl_prototyping.h \
		 hdrl_image.h \
		 hdrl_image_math.h \
		 hdrl_imagelist.h \
		 hdrl_imagelist_io.h \
		 hdrl_imagelist_basic.h \
		 hdrl_imagelist_view.h \
		 hdrl_parameter.h \
		 hdrl_parameter_defs.h \
		 hdrl_strehl.h \
		 hdrl_types.h \
		 hdrl_test.h \
		 hdrl_buffer.h \
		 hdrl_overscan_defs.h \
		 hdrl_image_defs.h \
		 hdrl_imagelist_defs.h \
		 hdrl_imagelist_basic_body.h \
		 hdrl_lacosmics.h \
		 hdrl_fit.h \
		 hdrl_fringe.h \
		 hdrl_catalogue.h \
		 hdrl_random.h \
		 hdrl.h

libhdrl_la_SOURCES = hdrl_overscan.c \
		    hdrl_sigclip.c \
		    hdrl_utils.c \
		    hdrl_combine.c \
		    hdrl_bpm_2d.c \
		    hdrl_bpm_3d.c \
		    hdrl_bpm_fit.c \
		    hdrl_bpm_utils.c \
		    hdrl_iter.c \
		    hdrl_frameiter.c \
		    hdrl_multiiter.c \
		    hdrl_collapse.c \
		    hdrl_elemop.c \
		    hdrl_flat.c \
		    hdrl_fit.c \
		    hdrl_strehl.c \
		    hdrl_prototyping.c \
		    hdrl_image.c \
		    hdrl_image_math.c \
		    hdrl_imagelist_io.c \
		    hdrl_imagelist_basic.c \
		    hdrl_imagelist_view.c \
		    hdrl_parameter.c \
		    hdrl_lacosmics.c \
		    hdrl_buffer.c \
		    hdrl_fringe.c \
		    hdrl_catalogue.c \
		    hdrl_random.c \
		    igam.c

libhdrl_la_CFLAGS = $(AM_CFLAGS)
libhdrl_la_LDFLAGS = $(CPL_LDFLAGS) -static
libhdrl_la_LIBADD = $(LIBCPLDFS) $(LIBCPLUI) $(LIBCPLDRS) $(LIBCPLCORE) catalogue/libcatalogue.la

if VLTSW
lib_LTLIBRARIES = libhdrl.la
pkginclude_HEADERS = $(hdrlHEADERS)
else
noinst_LTLIBRARIES = libhdrl.la
noinst_HEADERS = $(hdrlHEADERS)
endif
