## Process this file with automake to produce Makefile.in

##   This file is part of the VIMOS Pipeline
##   Copyright (C) 2000-2004 European Southern Observatory

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


AM_CPPFLAGS = $(HDRL_INCLUDES) $(MOSCA_INCLUDES) $(GSL_CFLAGS) $(all_includes)

pkginclude_HEADERS = vmadf.h vmadfifutable.h vmastrometrictable.h \
	vmccdtable.h vmcpl.h vmcube.h vmdetector.h vmdistmodels.h \
	vmextractiontable.h vmfit.h vmgalaxytable.h vmgrismtable.h \
	vmidstable.h vmifuextraction.h vmifufibers.h vmifuflat.h \
	vmifuimage.h vmifusky.h vmifusphotcalib.h vmifutable.h \
	vmifutypes.h vmifuutils.h vmimage.h vmimagearray.h vmimageset.h \
	vmimgastrometry.h vmimgphotcalib.h vmimgpreprocessing.h \
	vmimgextraction.h vmimgresampling.h vmimgutils.h vmlinecatalog.h \
	vmmath.h vmmatrix.h vmmosextraction.h vmmosflat.h vmmosfringes.h \
	vmmosmodels.h vmmossky.h vmmossphotcalib.h vmmosutils.h \
	vmmoswavecalib.h vmobjecttable.h vmphotometrictable.h \
	vmqcutils.h vmsextractor.h vmsphottable.h vmstarmatchtable.h \
	vmstartable.h vmtable.h vmtablearray.h vmtypes.h vmutils.h \
	vmwcsutils.h vmwindowtable.h vmifu.h vmspecphottable.h \
	vmstdfluxtable.h vmextincttable.h moses.h fors_tools.h \
    vimos_dfs.h vimos_calib_impl.h vimos_science_impl.h \
    vimos_calib_mult.h vimos_ifu_wcs.h two_d_linear_wcs.h vimos_utils.h \
	vimos_overscan.h vimos_detmodel.h vimos_detected_slits.h \
    vimos_calibrated_slits.h vimos_flat_normalise.h vimos_response.h \
	vimos_lines.h vimos_grism.h


privatelib_LTLIBRARIES = libvimos.la

libvimos_la_SOURCES = vmadf.c vmadfifutable.c vmastrometrictable.c \
	vmccdtable.c vmcpl.c vmcube.c vmdetector.c vmdistmodels.c \
	vmextractiontable.c vmfit.c vmgalaxytable.c vmgrismtable.c \
	vmidstable.c vmifuextraction.c vmifufibers.c vmifuflat.c \
	vmifuimage.c vmifusky.c vmifusphotcalib.c vmifutable.c \
	vmifutypes.c vmifuutils.c vmimage.c vmimagearray.c vmimageset.c \
	vmimgastrometry.c vmimgphotcalib.c vmimgpreprocessing.c \
	vmimgextraction.c vmimgresampling.c vmimgutils.c vmlinecatalog.c \
	vmmath.c vmmatrix.c vmmosextraction.c vmmosflat.c vmmosfringes.c \
	vmmosmodels.c vmmossky.c vmmossphotcalib.c vmmosutils.c \
	vmmoswavecalib.c vmobjecttable.c vmphotometrictable.c \
	vmqcutils.c vmsextractor.c vmsphottable.c vmstarmatchtable.c \
	vmstartable.c vmtable.c vmtablearray.c vmtypes.c vmutils.c \
	vmwcsutils.c vmwindowtable.c vmifu.c vmspecphottable.c \
	vmstdfluxtable.c vmextincttable.c moses.c fors_tools.c \
	vimos_dfs.c vimos_ifu_wcs.cc vimos_utils.cc vimos_utils.tcc \
	two_d_linear_wcs.cc vimos_overscan.cc vimos_detmodel.cc \
	vimos_calibrated_slits.cc vimos_detected_slits.cc \
	vimos_flat_normalise.cc vimos_response.cc vimos_calib_impl.cc \
	vimos_science_impl.cc vimos_calib_mult.cc vimos_lines.cc vimos_grism.cc

libvimos_la_LDFLAGS = $(HDRL_LDFLAGS) $(CPL_LDFLAGS) -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)
libvimos_la_LIBADD = $(HDRL_LIBS) $(LIBMOSCA) $(LIBIRPLIB) $(LIBPIL) $(LIBCFITSIO) $(LIBCPLDFS) $(LIBCPLDRS) $(LIBVIMOSWCS) $(LIBCPLUI) $(LIBCPLCORE) $(LIBCEXT)
libvimos_la_DEPENDENCIES = $(LIBMOSCA) $(LIBIRPLIB) $(LIBPIL) $(LIBVIMOSWCS)
