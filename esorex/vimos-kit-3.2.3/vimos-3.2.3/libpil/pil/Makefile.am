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


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


# CFITSIO_INCLUDES = -I$(CPLDIR)/include
LIBCFITSIO = -lcfitsio

AM_CPPFLAGS = $(INCKAZ) $(CFITSIO_INCLUDES) 

noinst_HEADERS = pilalias.h pilastroutils.h pilcatmap.h pilcdb.h \
	pildate.h pildfsconfig.h pildictionary.h pilerrno.h pilfileutils.h \
	pilfits.h pilframe.h pilframeset.h pilkeymap.h pilkeyword.h \
	pillist.h pilmacros.h pilmemory.h pilmessages.h pilpaf.h pilqc.h \
	pilrecipe.h pilstrutils.h piltask.h piltimer.h piltranslator.h \
	pilutils.h md5.h

noinst_LTLIBRARIES = libpil.la

libpil_la_SOURCES = pilalias.c pilastroutils.c pilcatmap.c pilcdb.c \
	pildate.c pildfsconfig.c pildictionary.c pilerrno.c pilfileutils.c \
	pilfits.c pilframe.c pilframeset.c pilkeymap.c pilkeyword.c \
	pillist.c pilmemory.c pilmessages.c pilpaf.c pilqc.c pilrecipe.c \
	pilstrutils.c piltask.c piltimer.c piltranslator.c pilutils.c \
	md5.c #xmemory.c

EXTRA_libpil_la_SOURCES = getopt.c getopt1.c getopt.h

libpil_la_LDFLAGS = $(CFITSIO_LDFLAGS) -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)
libpil_la_LIBADD = @GETOPT@ $(LIBCFITSIO) $(LIBKAZ)
libpil_la_DEPENDENCIES = @GETOPT@ $(LIBKAZ) 
