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

ACLOCAL_AMFLAGS = -I m4macros

DISTCLEANFILES = *~


COMPILE_FIRST = irplib config external mosca

TOPSUBDIRS =  libpil vimos regtests reflex

COMPILE_LAST = recipes tests


DOXYGEN_SUBDIRS = 

SUBDIRS = $(COMPILE_FIRST) $(TOPSUBDIRS) $(COMPILE_LAST)

EXTRA_DIST = m4macros/eso.m4 m4macros/cpl.m4 admin/doxygen.am \
	doxygen/Doxyfile.in


pipedocs_DATA = ChangeLog AUTHORS NEWS README TODO COPYING


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 config.h.in configure \
	stamp-h.in external/Makefile.in

clean-local: clean-doxygen

dist-hook: doxygen
	@if test -d $(top_srcdir)/html; then \
	  echo "cp -pr $(top_srcdir)/html $(distdir)"; \
	  cp -pr $(top_srcdir)/html $(distdir); \
	fi
	find $(distdir) -type d ! -perm -222 -exec chmod u+w {} \; -o \
	-type f ! -perm -222 -exec chmod u+w {} \; || chmod -R u+w $(distdir)

else

clean-local:

endif

uninstall-local: uninstall-doxygen

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

include $(top_srcdir)/admin/doxygen.am
