/* $Id: moses.h,v 1.10 2011-03-28 12:36:03 cizzo Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2006 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
 * $Author: cizzo $
 * $Date: 2011-03-28 12:36:03 $
 * $Revision: 1.10 $
 * $Name: not supported by cvs2svn $
 */

#ifndef MOSES_H
#define MOSES_H

#include <cpl.h>

CPL_BEGIN_DECLS

cpl_table *mos_hough_table(cpl_table *, const char *, const char *);
cpl_image *mos_remove_bias(cpl_image *, cpl_image *, cpl_table *);
cpl_image *mos_normalise_flat(cpl_image *, cpl_image *, cpl_table *, 
                              cpl_table *, double, double, double, double, 
                              int, int);
cpl_image *mos_normalise_longflat(cpl_image *, int, int, int);
cpl_error_code mos_arc_background_1D(float *, float *, int, int, int);
cpl_image *mos_arc_background(cpl_image *, int, int);
int mos_lines_width(const float *, int);
cpl_vector *mos_peak_candidates(const float *, int, float, float);
cpl_vector *mos_refine_peaks(const float *, int, cpl_vector *, int);
void mos_set_multiplex(int);
cpl_bivector *mos_identify_peaks(cpl_vector *, cpl_vector *,
                                 double, double, double);
cpl_bivector *mos_identify_peaks_new(cpl_vector *, cpl_vector *,
                                     double, double, double);
cpl_bivector *mos_find_peaks(const float *, int, cpl_vector *,
                             cpl_polynomial *, double, int);
double mos_eval_dds(cpl_polynomial *, double, double, double, double);
cpl_polynomial *mos_poly_wav2pix(cpl_bivector *, int, double, int,
                                 int *, double *,
                                 cpl_bivector **pixwav_used);
cpl_polynomial *mos_poly_pix2wav(cpl_bivector *, int, double, int,
                                 int *, double *);
cpl_image *mos_wavelength_calibration_raw(const cpl_image *, cpl_vector *,
                                          double, float, int, int, double, 
                                          double, double *, double *, int *, 
                                          double *, cpl_table *, cpl_image *, 
                                          cpl_image *, cpl_table *, cpl_mask *,
                                          cpl_table *detected_lines, 
                                          double disp_tolerance,
                                          double ratio_tolerance);
cpl_error_code mos_interpolate_wavecalib_slit(cpl_table *, cpl_table *, 
                                              int, int);
cpl_error_code mos_interpolate_wavecalib(cpl_table *, cpl_image *, int, int);
int mos_clean_outliers(cpl_table *, const char *);
cpl_table *mos_locate_spectra(cpl_mask *);
cpl_error_code mos_validate_slits(cpl_table *);
cpl_error_code mos_rotate_slits(cpl_table *, int, int, int);
cpl_table *mos_identify_slits(cpl_table *, cpl_table *, cpl_table *);
cpl_table *mos_identify_slits_fast(cpl_table *, cpl_table *, cpl_table *);
cpl_table *mos_identify_slits_linear(cpl_table *, cpl_table *);
cpl_table *mos_trace_flat(cpl_image *, cpl_table *, double,
                          double, double, double);
cpl_table *mos_poly_trace(cpl_table *, cpl_table *, int);
cpl_error_code mos_global_trace(cpl_table *, cpl_table *, int);
cpl_image *mos_spatial_calibration(cpl_image *, cpl_table *, cpl_table *, 
                                   double, double, double, double, int,
                                   cpl_image *);
cpl_image *mos_wavelength_calibration_final(cpl_image *, cpl_table *,
                                            cpl_vector *, double, float, int, 
                                            int, double, double, double *, 
                                            double *, int *, double *, 
                                            cpl_table *, cpl_image *, 
                                            cpl_image *, cpl_table *,
                                            cpl_table *detected_lines,  
                                            double disp_tolerance,
                                            double ratio_tolerance);
cpl_table *mos_global_distortion(cpl_table *, cpl_table *,
                                 cpl_table *, cpl_table *, double);
cpl_table *mos_average_global_distortion(cpl_table **, int, double, double);
cpl_table *mos_build_slit_location(cpl_table *, cpl_table *, int);
cpl_table *mos_build_curv_coeff(cpl_table *, cpl_table *, cpl_table *);
cpl_table *mos_build_disp_coeff(cpl_table *, cpl_table *);
cpl_table *mos_wavelength_align(cpl_image *, cpl_table *, double, double,
                                double, cpl_table *, cpl_vector *, int, 
                                int, cpl_image *, int);
cpl_table *mos_wavelength_align_lss(cpl_image *, double, double, double,
                              cpl_table *, cpl_vector *, int, int, 
                              cpl_image *, int);
cpl_image *mos_wavelength_calibration(cpl_image *, double, double, double, 
                                      double, cpl_table *, int);
cpl_image *mos_map_pixel(cpl_table *, double, double, double, double, int);
cpl_image *mos_map_idscoeff(cpl_table *, int, double, double, double);
cpl_image *mos_map_wavelengths(cpl_image *, cpl_image *, cpl_table *, 
                               cpl_table *, double, double, double, double);
cpl_image *mos_map_spectrum(cpl_image *, cpl_image *, cpl_image *, cpl_table *,
                            cpl_table *, double, double, double, double, int);
cpl_image *mos_subtract_sky(cpl_image *, cpl_table *, cpl_table *, double,
                            double, double, double);
cpl_table *mos_sky_map(cpl_image *, cpl_image *, double, cpl_image *);
cpl_table *mos_sky_map_super(cpl_image *, cpl_image *, double, double, int,
                             cpl_image *);
cpl_image *mos_sky_local_old(cpl_image *spectra, cpl_table *slits);
cpl_image *mos_sky_local(cpl_image *, cpl_table *, int);

cpl_error_code mos_clean_cosmics(cpl_image *, float, float, float);
cpl_error_code mos_clean_bad_pixels(cpl_image *, cpl_table *, int);

double mos_distortions_rms(cpl_image *, cpl_vector *, double, double, int, int);

cpl_image *mos_spatial_map(cpl_image *, cpl_table *, cpl_table *, double,
                           double, double, double);
cpl_image *mos_detect_objects(cpl_image *, cpl_table *, int, int, int, double);
cpl_image **mos_extract_objects(cpl_image *science, cpl_image *science_var, 
                                cpl_image *sky,
                                cpl_table *objects, int extraction, double ron,
                                double gain, int ncombined);
int mos_spectral_resolution(cpl_image *, double, double, double, int,
                            double *, double *, double *, double *, int *);
cpl_table *mos_resolution_table(cpl_image *image, double startwave,
                                double dispersion, int saturation,
                                cpl_vector *lines);
double mos_integrate_signal(cpl_image *, cpl_image *, int, int, double, double);

cpl_polynomial *mos_montecarlo_polyfit(cpl_table *points, cpl_table *evaluate,
                                       int samples, int order);

cpl_error_code mos_randomise_image(cpl_image *, double ron,
                                   double gain, double bias);

/* Instrument dependent!  */

double mos_get_gain_vimos(cpl_propertylist *);
cpl_table *mos_load_overscans_vimos(const cpl_propertylist *, int);
cpl_table *mos_load_overscans_fors(const cpl_propertylist *);
cpl_table *mos_load_slits_vimos(cpl_propertylist *, int include_ref);
int mos_check_multiplex(cpl_table *);
int mos_check_multiplex_old(cpl_table *);
int mos_assign_multiplex_group(cpl_table *);
cpl_table *mos_load_slits_fors_mxu(cpl_propertylist *);
cpl_table *mos_load_slits_fors_mos(cpl_propertylist *);
cpl_table *mos_load_slits_fors_lss(cpl_propertylist *);
cpl_table *mos_load_slits_fors_pmos(cpl_propertylist *);

cpl_error_code mos_refmask_find_gaps(cpl_mask  *, cpl_image *, double);
cpl_error_code mos_saturation_process(cpl_image *);
cpl_error_code mos_subtract_background(cpl_image *);

cpl_error_code mos_object_intersect(cpl_table **, cpl_table *, int, float);

int mos_get_maxobjs_per_slit(cpl_table *);
int mos_get_nobjects(cpl_table *);

int mos_check_slits(cpl_table *, float);
int *fors_get_nobjs_perslit(cpl_table *);
int mos_rebin_signal(cpl_image **, int);
int mos_rebin_error(cpl_image **, int);
cpl_table *mos_photometric_calibration(cpl_image *, double, double, double,
                                 double, cpl_table *, double, cpl_table *,
                                 int);
int map_table(cpl_image *, double, double, cpl_table *, 
              const char *, const char *);
cpl_image *mos_ksigma_stack(cpl_imagelist *imlist,
                            double klow, double khigh, int kiter,
                            cpl_image **);
cpl_image *mos_apply_photometry(cpl_image *, cpl_table *response,
                                cpl_table *ext_table, double startwave,
                                double dispersion, double gain,
                                double exptime, double airmass);
cpl_image *mos_propagate_photometry_error(cpl_image *, cpl_image *,
                                cpl_table *response,
                                cpl_table *ext_table, double startwave,
                                double dispersion, double gain,
                                double exptime, double airmass);
int mos_check_polarisation(cpl_image *q_image, cpl_image *q_error,
                           cpl_image *u_image, cpl_image *u_error,
                           double startwave, double dispersion,
                           double band, cpl_table *pol_sta,
                           double ra, double dec, char *filter,
                           int *polarisation,
                           double *p_offset, double *p_error,
                           double *a_offset, double *a_error);
int mos_compute_offset(cpl_table *, cpl_table *, double *);
cpl_error_code mos_image_shift(cpl_image *, double dx, double dy);
int mos_slit_closest_to_center(cpl_table *slits, int nx, int ny);
cpl_error_code mos_extract_flux(cpl_image *, cpl_table *, double, double,
                                int, double, double *, double *);
cpl_error_code mos_extract_flux_mapped(cpl_image *, cpl_table *,
                                       double, double,
                                       double lambda, double startwave,
                                       double dispersion, int dx, double gain,
                                       double *o_flux, double *o_err);
int mos_median_in_slit(cpl_table *, cpl_table *, int slit,
                       char *label, double *mvalue);
cpl_image *mos_image_filter_median(cpl_image *image, int nx, int ny);

cpl_error_code mos_slit_wavemap(cpl_image *wavemap, int slit, cpl_table *slits,
                                cpl_table *polytraces, double reference,
                                double blue, double red, double dispersion);

CPL_END_DECLS

#endif   /* MOSES_H */
