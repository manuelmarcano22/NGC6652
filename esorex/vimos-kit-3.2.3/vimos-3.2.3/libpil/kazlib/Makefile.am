## Process this file with automake to produce Makefile.in

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

EXTRA_DIST = m4macros/kazlib.m4 CHANGES MUST_READ Makefile.gcc Makefile.vc \
	blast.pl docs.ist docs.ltx


## Build library either for installation or as convenience library for some
## software project that wants to use the library just internally.
if KAZLIB_CONVENIENCE

noinst_HEADERS = dict.h except.h hash.h list.h sfx.h
noinst_LTLIBRARIES = libkazc.la

libkazc_la_SOURCES = dict.c hash.c list.c

EXTRA_libkazc_la_SOURCES = except.c sfx.c

libkazc_la_LIBADD = @EXCEPTIONS@
libkazc_la_DEPENDENCIES = @EXCEPTIONS@

else

include_HEADERS = dict.h except.h hash.h list.h sfx.h
lib_LTLIBRARIES = libkaz.la

libkaz_la_SOURCES = dict.c hash.c list.c

EXTRA_libkaz_la_SOURCES = except.c sfx.c

libkaz_la_LIBADD = @EXCEPTIONS@
libkaz_la_DEPENDENCIES = @EXCEPTIONS@

endif


check_PROGRAMS = tlist tdict thash teh tsfx

tlist_SOURCES = list.c
tlist_CFLAGS = -DKAZLIB_TEST_MAIN

tdict_SOURCES = dict.c
tdict_CFLAGS = -DKAZLIB_TEST_MAIN

thash_SOURCES = hash.c
thash_CFLAGS = -DKAZLIB_TEST_MAIN

teh_SOURCES = except.c
teh_CFLAGS = -DKAZLIB_TEST_MAIN

tsfx_SOURCES = sfx.h sfx.c
tsfx_CFLAGS = -DKAZLIB_TEST_MAIN
tsfx_LDADD = hash.o except.o


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 config.h.in config.log \
	configure stamp-h.in docs.aux docs.log docs.dvi docs.toc docs.idx \
	docs.ind docs.ilg

dvi-local: docs.dvi

docs.dvi: docs.ltx docs.toc docs.ind
	@echo $(LATEX) docs.ltx; \
	$(LATEX) docs.ltx; $(LATEX) docs.ltx

docs.ind: docs.idx docs.ist
	$(MAKEINDEX) -s docs.ist docs.idx

docs.toc: docs.ltx
	$(LATEX) docs.ltx
	rm -f docs.dvi

docs.idx: docs.ltx
	$(LATEX) docs.ltx
	rm -f docs.dvi

dist-hook: dvi
	@if test -f $(srcdir)/docs.dvi; then \
	  echo "cp -p $(srcdir)/docs.dvi $(distdir)"; \
	  cp -p $(srcdir)/docs.dvi $(distdir); \
	fi
	find $(distdir) -type d ! -perm -222 -exec chmod u+w {} \; -o \
	-type f ! -perm -222 -exec chmod u+w {} \; || chmod -R u+w $(distdir)

endif
