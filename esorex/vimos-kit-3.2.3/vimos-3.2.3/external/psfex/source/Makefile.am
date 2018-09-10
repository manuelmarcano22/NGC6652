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


EXTRA_DIST = Macrodefs. Macrodefs.aix Macrodefs.alpha Macrodefs.debug \
	Macrodefs.gcc Macrodefs.hpux Macrodefs.insure Macrodefs.linux \
	Macrodefs.linuxp2 Macrodefs.linuxpc Macrodefs.profile \
	Macrodefs.purify Macrodefs.sgi Macrodefs.solaris Macrodefs.sunos \
	Macrodefs.ultrix Makefile.orig


psfex_bin_PROGRAMS = psfex

noinst_HEADERS = check.h define.h fft.h fitscat.h fitscat_defs.h fitshead.h \
	fitsutil.h globals.h key.h poly.h prefs.h psf.h random.h sample.h \
	types.h vignet.h

psfex_SOURCES = basis.c check.c fft.c fitscat.c fitsconv.c fitsfield.c \
	fitshead.c fitskey.c fitsmisc.c fitsread.c fitstab.c fitsutil.c \
	fitswrite.c main.c makeit.c misc.c pc.c poly.c prefs.c psf.c \
	random.c sample.c vignet.c

psfex_LDADD = $(LIBM)


COMPILE = $(CC) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
