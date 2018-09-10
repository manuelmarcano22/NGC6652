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


SUBDIRS = wcs

EXTRA_DIST = Macrodefs. Macrodefs.aix Macrodefs.alpha Macrodefs.debug \
	Macrodefs.gcc Macrodefs.hpux Macrodefs.insure Macrodefs.linux \
	Macrodefs.linuxp2 Macrodefs.linuxpc Macrodefs.profile \
	Macrodefs.purify Macrodefs.sgi Macrodefs.solaris Macrodefs.sunos \
	Macrodefs.ultrix Makefile.orig fitsfield.c sexchanges.log


sext_bin_PROGRAMS = sex

noinst_HEADERS = assoc.h astrom.h back.h bpro.h check.h clean.h define.h \
	extract.h field.h filter.h fitscat.h fitscat_defs.h fitshead.h \
	fitsutil.h flag.h globals.h growth.h image.h interpolate.h key.h \
	neurro.h param.h photom.h plist.h poly.h prefs.h psf.h retina.h \
	sexhead.h sexhead1.h sexheadsc.h som.h types.h weight.h

sex_SOURCES = analyse.c assoc.c astrom.c back.c bpro.c catout.c check.c \
	clean.c extract.c field.c filter.c fitscat.c fitsconv.c fitshead.c \
	fitskey.c fitsmisc.c fitsread.c fitstab.c fitsutil.c fitswrite.c \
	flag.c graph.c growth.c image.c interpolate.c main.c makeit.c \
	manobjlist.c misc.c neurro.c pc.c photom.c plist.c poly.c prefs.c \
	psf.c readimage.c refine.c retina.c scan.c som.c weight.c

sex_LDFLAGS = -L./wcs 
sex_LDADD = -lwcs_c $(LIBM)


STDH  = define.h globals.h types.h
STDO  =

COMPILE = $(CC) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)


analyse.o: $(STDH) $(STDO) analyse.c assoc.h back.h check.h growth.h image.h \
	photom.h plist.h growth.h psf.h

assoc.o: $(STDH) $(STDO) assoc.c assoc.h

astrom.o: $(STDH) $(STDO) astrom.c astrom.h

back.o: $(STDH) $(STDO) back.c back.h field.h

bpro.o: $(STDH) $(STDO) bpro.c bpro.h fitscat.h

catout.o: $(STDH) $(STDO) catout.c param.h sexhead.h sexhead1.h sexheadsc.h

check.o: $(STDH) $(STDO) check.c check.h

clean.o: $(STDH) $(STDO) clean.c clean.h check.h

extract.o: $(STDH) $(STDO) extract.c extract.h plist.h

field.o: $(STDH) $(STDO) field.c field.h interpolate.h

filter.o: $(STDH) $(STDO) filter.c bpro.h filter.h image.h

fitscat.o: $(STDH) $(STDO) fitscat.c fitscat_defs.h fitscat.h

fitsconv.o: $(STDH) $(STDO) fitsconv.c fitscat_defs.h fitscat.h

fitshead.o: $(STDH) $(STDO) fitshead.c fitscat_defs.h fitscat.h

fitskey.o: $(STDH) $(STDO) fitskey.c fitscat_defs.h fitscat.h

fitsmisc.o: $(STDH) $(STDO) fitsmisc.c fitscat_defs.h fitscat.h

fitsread.o: $(STDH) $(STDO) fitsread.c fitscat_defs.h fitscat.h

fitstab.o: $(STDH) $(STDO) fitstab.c fitscat_defs.h fitscat.h

fitsutil.o: $(STDH) $(STDO) fitsutil.c fitscat_defs.h fitscat.h

fitswrite.o: $(STDH) $(STDO) fitswrite.c fitscat_defs.h fitscat.h

flag.o: $(STDH) $(STDO) flag.c flag.h plist.h

graph.o: $(STDH) $(STDO) graph.c

growth.o: $(STDH) $(STDO) growth.c growth.h

image.o: $(STDH) $(STDO) image.c

interpolate.o: $(STDH) $(STDO) interpolate.c interpolate.h

main.o: $(STDH) $(STDO) main.c

makeit.o: $(STDH) $(STDO) makeit.c check.h growth.h interpolate.h psf.h \
	weight.h

manobjlist.o: $(STDH) $(STDO) manobjlist.c plist.h

misc.o: $(STDH) $(STDO) misc.c

neurro.o: $(STDH) $(STDO) neurro.c neurro.h

pc.o: $(STDH) $(STDO) pc.c image.h psf.h

photom.o: $(STDH) $(STDO) photom.c photom.h

plist.o: $(STDH) $(STDO) plist.c plist.h

poly.o: $(STDH) $(STDO) poly.c poly.h

prefs.o: $(STDH) $(STDO) prefs.c prefs.h

psf.o: $(STDH) $(STDO) psf.c image.h psf.h

readimage.o: $(STDH) $(STDO) readimage.c check.h interpolate.h weight.h

refine.o: $(STDH) $(STDO) refine.c extract.h plist.h

retina.o: $(STDH) $(STDO) retina.c image.h retina.h fitscat.h

scan.o: $(STDH) $(STDO) scan.c check.h clean.h extract.h plist.h

som.o: $(STDH) $(STDO) som.c som.h

weight.o: $(STDH) $(STDO) weight.c weight.h field.h
