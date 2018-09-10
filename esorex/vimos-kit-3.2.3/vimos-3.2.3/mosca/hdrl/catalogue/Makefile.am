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

# NOTE: Only build this directory. The tests subdir is built from the top level
# Makefile to maintain the correct build order and break the circular dependency
# between libcatalogue, libhdrl and the unit test binaries.
SUBDIRS = .

if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif

AM_CPPFLAGS = -DCX_LOG_DOMAIN=\"HdrlLib\"  -DHDRL_USE_PRIVATE $(CPL_INCLUDES) -I$(top_srcdir)

noinst_LTLIBRARIES = libcatalogue.la

noinst_HEADERS = \
		 casu_utilfunctions.h  \
		 ap.h  \
		 classify.h  \
		 floatmath.h  \
		 imcore.h  \
		 imcore_radii.h  \
		 imcore_version.h  \
		 util.h

libcatalogue_la_SOURCES = \
		    apclust.c \
		    apinit.c \
		    apline.c \
		    areals.c \
		    classify.c \
		    create_table_6.c \
		    create_table.c \
		    imcore_background.c \
		    imcore_conf.c \
		    imcore_extend.c \
		    imcore_filter.c \
		    imcore_overlp.c \
		    imcore_phopt.c \
		    imcore_radii.c \
		    moments.c \
		    casu_utilfunctions.c \
		    polynm.c \
		    seeing.c \
		    solve.c \
		    terminate.c
