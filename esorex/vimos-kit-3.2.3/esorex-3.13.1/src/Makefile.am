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
##  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA


AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILE = *~

SUBDIRS = . tests


if MAINTAINER_MODE

MAINTAINERCLEANFILES = $(srcdir)/Makefile.in

endif


LIBESOREX = $(builddir)/libesorex.la

if HAVE_LIBAVCALL
AM_CPPFLAGS = -DCX_LOG_DOMAIN=\"EsoRex\" -DSYSCONFDIR=\"@sysconfdir@\" \
              $(CPL_INCLUDES) $(CX_INCLUDES) $(INCLTDL) $(AVCALL_INCLUDES)
EXTRA_FFI_LIBS = $(LIBAVCALL)
else
if HAVE_LIBFFI
AM_CPPFLAGS = -DCX_LOG_DOMAIN=\"EsoRex\" -DSYSCONFDIR=\"@sysconfdir@\" \
              $(CPL_INCLUDES) $(CX_INCLUDES) $(INCLTDL) $(FFI_INCLUDES)
EXTRA_FFI_LIBS = $(LIBFFI)
else
AM_CPPFLAGS = -DCX_LOG_DOMAIN=\"EsoRex\" -DSYSCONFDIR=\"@sysconfdir@\" \
              $(CPL_INCLUDES) $(CX_INCLUDES) $(INCLTDL)
EXTRA_FFI_LIBS =
endif
endif


bin_PROGRAMS = esorex

noinst_LTLIBRARIES = libesorex.la


noinst_HEADERS = \
    er_fileutils.h \
    er_help.h \
    er_macros.h \
    er_params.h \
    er_paramutils.h \
    er_plugin.h \
    er_pluginlist.h \
    er_stringarray.h \
    er_stringutils.h \
    er_paf.h \
    er_json.h \
    er_python.h

libesorex_la_SOURCES = \
    er_fileutils.c \
    er_help.c er_params.c \
    er_paramutils.c \
    er_plugin.c \
    er_pluginlist.c \
    er_stringarray.c \
    er_stringutils.c \
    er_paf.c \
    er_json.c

if ENABLE_PYTHON_RECIPES
libesorex_la_SOURCES += er_python.c
else
EXTRA_FFI_LIBS :=
endif


esorex_SOURCES = er_main.cpp 

EXTRA_esorex_SOURCES = getopt.c getopt1.c getopt.h

config_DATA =

esorex_LDFLAGS = $(CPL_LDFLAGS) $(LIBLTDL)
esorex_LDADD = @GETOPT@ $(LIBESOREX) $(EXTRA_FFI_LIBS) \
               $(LIBCPLDFS) $(LIBCPLUI) $(LIBCPLCORE) $(LIBCEXT) $(LIBLTDL)
esorex_DEPENDENCIES = @GETOPT@ $(LIBESOREX)
