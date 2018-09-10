## Process this file with automake to produce Makefile.in

##   This file is part of the Refran Library
##   Copyright (C) 2009 César Enrique García

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
##   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA


AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~

SUBDIRS = tests

if MAINTAINER_MODE

MAINTAINERCLEANFILES = $(srcdir)/Makefile.in

endif

#AM_CXXFLAGS = -Wall -Wextra -pedantic -std=c++0x
AM_CXXFLAGS = -Wall -Wextra 
AM_CFLAGS = -Wall -Wextra -std=c99 


AM_CPPFLAGS = $(HDRL_INCLUDES) $(CPL_INCLUDES) $(GSL_CFLAGS)

pkginclude_HEADERS = ccd_config.h calibrated_slit.h \
                  detected_slits.h detected_slit.h \
                  extinction.h fiera_config.h \
                  flat_combine.h flat_combine.tcc grism_config.h \
                  global_distortion.h \
                  image_utils.h image_utils.tcc \
                  image_normalisation.h image_normalisation.tcc \
                  image_smooth.h image_smooth.cpp \
                  image_spline_fit.h image_spline_fit.tcc \
                  mosca_image.h mosca_image.tcc \
                  rect_region.h reduce_method.h response.h \
                  spectrum.h spec_std_star.h spatial_distortion.h \
                  slit_trace_distortion.h \
                  statistics.h statistics.cpp \
                  type_traits.h vector_utils.h vector_utils.tcc \
                  wavelength_calibration.h 
 
noinst_HEADERS =  

privatelib_LTLIBRARIES = libmosca.la

libmosca_la_SOURCES = ccd_config.cpp calibrated_slit.cpp detected_slits.cpp \
                      detected_slit.cpp extinction.cpp \
                      fiera_config.cpp global_distortion.cpp grism_config.cpp \
                      mosca_image.cpp rect_region.cpp \
                      reduce_method.cpp response.cpp slit_trace_distortion.cpp \
                      spectrum.cpp spatial_distortion.cpp \
                      spec_std_star.cpp vector_utils.cpp \
                      wavelength_calibration.cpp 

#The -static is needed because in VLT machines, with LD_LIBRARY_PATH=/usrlib:/vlt..., the /usr/lib/libstdc++.so is chosen at runtime, which is the wrong version (it was compiled with /vlt/VLT2010/....). This will compile everything static for this library. Once DFS machines are VLT-free this can be removed.
libmosca_la_LDFLAGS = -version-info 0:1:0 $(HDRL_LDFLAGS) $(CPL_LDFLAGS)
libmosca_la_LIBADD =  -lm $(HDRL_LIBS) $(LIBCPLCORE) $(GSL_LIBS)
libmosca_la_DEPENDENCY =  $(LIBHDRL)
