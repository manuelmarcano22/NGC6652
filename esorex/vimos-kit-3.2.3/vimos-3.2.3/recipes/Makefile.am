## Process this file with automake to produce Makefile.in

##   This file is part of the VIMOS Pipeline
##   Copyright (C) 2001-2004 European Southern Observatory

##   This library is free software; you can redistribute it and/or
##   modify it under the terms of the GNU Library General Public
##   License as published by the Free Software Foundation; either
##   version 2 of the License, or (at your option) any later version.

##   This library is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##   Library General Public License for more details.

##   You should have received a copy of the GNU Library General Public License
##   along with this library; see the file COPYING.LIB.  If not, write to
##   the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
##   Boston, MA 02111-1307, USA.

AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


AM_CPPFLAGS = $(HDRL_INCLUDES) $(MOSCA_INCLUDES) $(GSL_CFLAGS) $(all_includes)


noinst_HEADERS = 

if ONLINE_MODE
plugin_LTLIBRARIES = vmbias.la vmdark.la vmdet.la vmimcalphot.la \
                     vmimflatscreen.la vmimflatsky.la vmimobsjitter.la \
                     vmimobsstare.la vmimpreimaging.la vmimstandard.la \
                     vmmasktoccd.la vmmosobsjitter.la vmmosobsstare.la \
                     vmskyccd.la vmspcaldisp.la vmspflat.la vmifucalib.la \
                     vmifuscience.la vmifustandard.la vmifucombine.la \
                     vmmosstandard.la vmmoscombine.la vmspphot.la \
                     vmifucombinecube.la vmmoscalib.la \
                     vmmosscience.la
else
plugin_LTLIBRARIES = vmbias.la vmdark.la vmdet.la vmimcalphot.la \
                     vmimflatscreen.la vmimflatsky.la vmimobsjitter.la \
                     vmimobsstare.la vmimpreimaging.la vmimstandard.la \
                     vmmosobsjitter.la vmmosobsstare.la  vmspcaldisp.la \
                     vmspflat.la vmifucalib.la vmifuscience.la \
                     vmifustandard.la vmifucombine.la vmmosstandard.la \
                     vmmoscombine.la vmspphot.la \
                     vmifucombinecube.la vmmoscalib.la vmmosscience.la
endif

vmbias_la_SOURCES = vmbias.cc
vmbias_la_LIBADD = $(HDRL_LIBS) $(LIBMOSCA) $(LIBVIMOS) $(GSL_LIBS) $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmbias_la_LDFLAGS = -module -avoid-version -no-undefined $(HDRL_LDFLAGS) $(CPL_LDFLAGS)
vmbias_la_DEPENDENCIES = $(LIBVIMOS) 

vmdark_la_SOURCES = vmdark.c
vmdark_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmdark_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmdark_la_DEPENDENCIES = $(LIBVIMOS) 

vmifucalib_la_SOURCES = vmifucalib.c
vmifucalib_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmifucalib_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmifucalib_la_DEPENDENCIES = $(LIBVIMOS) 

vmifucombine_la_SOURCES = vmifucombine.c
vmifucombine_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmifucombine_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmifucombine_la_DEPENDENCIES = $(LIBVIMOS) 

vmifucombinecube_la_SOURCES = vmifucombinecube.cc
vmifucombinecube_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmifucombinecube_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmifucombinecube_la_DEPENDENCIES = $(LIBVIMOS) 

vmifuscience_la_SOURCES = vmifuscience.c
vmifuscience_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmifuscience_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmifuscience_la_DEPENDENCIES = $(LIBVIMOS) 

vmifustandard_la_SOURCES = vmifustandard.c
vmifustandard_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmifustandard_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmifustandard_la_DEPENDENCIES = $(LIBVIMOS) 

vmmosstandard_la_SOURCES = vmmosstandard.c
vmmosstandard_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmmosstandard_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmmosstandard_la_DEPENDENCIES = $(LIBVIMOS) 

vmmoscombine_la_SOURCES = vmmoscombine.c
vmmoscombine_la_LIBADD = $(LIBVIMOS)  $(LIBCFITSIO) $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmmoscombine_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmmoscombine_la_DEPENDENCIES = $(LIBVIMOS) 

vmspphot_la_SOURCES = vmspphot.c
vmspphot_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmspphot_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmspphot_la_DEPENDENCIES = $(LIBVIMOS) 

vmdet_la_SOURCES = vmdet.c
vmdet_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmdet_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmdet_la_DEPENDENCIES = $(LIBVIMOS) 

vmimcalphot_la_SOURCES = vmimcalphot.c
vmimcalphot_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmimcalphot_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmimcalphot_la_DEPENDENCIES = $(LIBVIMOS) 

vmimflatscreen_la_SOURCES = vmimflatscreen.c
vmimflatscreen_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmimflatscreen_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmimflatscreen_la_DEPENDENCIES = $(LIBVIMOS) 

vmimflatsky_la_SOURCES = vmimflatsky.c
vmimflatsky_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmimflatsky_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmimflatsky_la_DEPENDENCIES = $(LIBVIMOS) 

vmimobsjitter_la_SOURCES = vmimobsjitter.c
vmimobsjitter_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmimobsjitter_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmimobsjitter_la_DEPENDENCIES = $(LIBVIMOS) 

vmimobsstare_la_SOURCES = vmimobsstare.c
vmimobsstare_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmimobsstare_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmimobsstare_la_DEPENDENCIES = $(LIBVIMOS) 

vmimpreimaging_la_SOURCES = vmimpreimaging.c
vmimpreimaging_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmimpreimaging_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmimpreimaging_la_DEPENDENCIES = $(LIBVIMOS) 

vmimstandard_la_SOURCES = vmimstandard.c
vmimstandard_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmimstandard_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmimstandard_la_DEPENDENCIES = $(LIBVIMOS) 

if ONLINE_MODE
vmmasktoccd_la_SOURCES = vmmasktoccd.c
else
vmmasktoccd_la_SOURCES =
endif
vmmasktoccd_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmmasktoccd_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmmasktoccd_la_DEPENDENCIES = $(LIBVIMOS) 

vmmosobsjitter_la_SOURCES = vmmosobsjitter.c
vmmosobsjitter_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmmosobsjitter_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmmosobsjitter_la_DEPENDENCIES = $(LIBVIMOS) 

vmmosobsstare_la_SOURCES = vmmosobsstare.c
vmmosobsstare_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmmosobsstare_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmmosobsstare_la_DEPENDENCIES = $(LIBVIMOS) 

if ONLINE_MODE
vmskyccd_la_SOURCES = vmskyccd.c
else
vmskyccd_la_SOURCES =
endif
vmskyccd_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmskyccd_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmskyccd_la_DEPENDENCIES = $(LIBVIMOS) 

vmspcaldisp_la_SOURCES = vmspcaldisp.c
vmspcaldisp_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmspcaldisp_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmspcaldisp_la_DEPENDENCIES = $(LIBVIMOS) 

vmspflat_la_SOURCES = vmspflat.c
vmspflat_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmspflat_la_LDFLAGS = -module -avoid-version -no-undefined $(CPL_LDFLAGS)
vmspflat_la_DEPENDENCIES = $(LIBVIMOS) 

vmmoscalib_la_SOURCES = vimos_calib.cc
vmmoscalib_la_LIBADD = $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmmoscalib_la_LDFLAGS = -module -avoid-version  -no-undefined $(CPL_LDFLAGS)
vmmoscalib_la_DEPENDENCIES = $(LIBVIMOS) 

vmmosscience_la_SOURCES = vimos_science.cc
vmmosscience_la_LIBADD = $(HDRL_LIBS) $(LIBMOSCA) $(LIBVIMOS)  $(LIBCEXT) $(LIBCPLCORE) $(LIBCPLDFS) $(LIBCPLUI)
vmmosscience_la_LDFLAGS = -module -avoid-version -no-undefined $(HDRL_LDFLAGS) $(CPL_LDFLAGS)
vmmosscience_la_DEPENDENCIES = $(LIBVIMOS) 
