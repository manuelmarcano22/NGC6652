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

AUTOMAKE_OPTIONS = 1.8 foreign no-dependencies

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


EXTRA_DIST = $(SEXTRACTOR_EXTRAS)


noinst_LIBRARIES = libwcs_c.a

noinst_HEADERS = cel.h lin.h proj.h wcs.h wcstrig.h

libwcs_c_a_SOURCES = cel.c lin.c proj.c sph.c wcs.c wcstrig.c

DEFS = -I. -I$(srcdir) -DCOPYSIGN -DSIGNBIT

SEXTRACTOR_EXTRAS = wcs_ori README.SEX Makefile.orig tcel.c tlin.c tpgc.c
	 tproj1.c tproj2.c tsph.c twcs1.c twcs2.c tca.map
