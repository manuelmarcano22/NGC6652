## Process this file with automake to produce Makefile.in

##   This file is part of the VIMOS Pipeline
##   Copyright (C) 2001-2004 European Southern Observatory

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

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 config.h.in config.log \
        configure stamp-h.in

endif


EXTRA_DIST = Files Readme $(LIBVIMOSWCS_EXTRAS)


noinst_LTLIBRARIES = libvimoswcs.la

noinst_HEADERS = fitsfile.h fitshead.h imio.h lvimoswcs.h vimoswcs.h vimoswcscat.h vimoswcslib.h

libvimoswcs_la_SOURCES = actread.c binread.c catutil.c cel.c ctgread.c daoread.c \
	dateutil.c dsspos.c fileutil.c findstar.c fitsfile.c fitsvimoswcs.c \
	gscread.c hget.c hput.c iget.c imgetvimoswcs.c imhfile.c imio.c \
	imrotate.c imsetvimoswcs.c lin.c matchstar.c platefit.c platepos.c \
	proj.c slasubs.c sortstar.c sph.c tabread.c tnxpos.c \
	ty2read.c uacread.c ujcread.c vimoswcs.c vimoswcscon.c vimoswcsinit.c vimoswcslib.c \
	vimoswcstrig.c webread.c worldpos.c

LIBVIMOSWCS_EXTRAS = agascread.c catread.c fortcat.c fortvimoswcs.c imvimoswcs.c median.c \
	polyfit.c tnxfit.c vimoswcsfort.c vimoswcsinit0.c vimoswcsinit1.c vimoswcsinitm.c \
	vimoswcspoly.c
