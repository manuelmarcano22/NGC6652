/*
 * This file is part of the VIMOS Data Reduction Pipeline
 * Copyright (C) 2002-2010 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 * vimos_response.cpp
 *
 *  Created on: 2014 4 2
 *      Author: cgarcia
 */

#include <stdexcept>
#include <vector>
#include <sstream>
#include <cmath>
#include <vimos_response.h>
#include "extinction.h"
#include "spec_std_star.h"
#include "response.h"
#include "vimos_dfs.h"
#include "vimos_flat_normalise.h"

#define MAX_COLNAME      (80)

cpl_table * vimos_compute_response
(cpl_image *spectra, cpl_image * mapped_flat_sed,
 cpl_propertylist *flat_sed_header, cpl_table *objects,
 double startwave, double dispersion, double gain,
 double exptime, cpl_table *ext_table, double airmass, cpl_table *flux_table,
 const std::vector<double>& ignored_waves,
 const std::vector<std::pair<double, double> >& ignored_wave_ranges,
 int nknots, int degree, cpl_table *& response_interp,
 double& flat_sed_norm_factor,
 const vimos::detected_slits& det_slits)
{

    cpl_image *spectrum       = NULL; // Extracted standard star spectrum
    cpl_table *response_table;
    int        nx, ny;


    if (spectra == NULL || ext_table == NULL || flux_table == NULL) 
        throw std::invalid_argument("Empty spectra, ext_table or flux_table");

    if (!cpl_table_has_column(ext_table, "WAVE")) 
        throw std::invalid_argument("Column WAVE in atmospheric extinction table");

    if (!cpl_table_has_column(ext_table, "EXTINCTION")) 
        throw std::invalid_argument("Column EXTINCTION in atmospheric extinction table");

    if (!cpl_table_has_column(flux_table, "WAVE")) 
        throw std::invalid_argument("Column WAVE in standard star flux table");

    if (!cpl_table_has_column(flux_table, "FLUX")) 
        throw std::invalid_argument("Column FLUX in standard star flux table");

    if (gain < 0.1) 
        throw std::invalid_argument("Invalid gain factor (<0.1)");

    if (exptime < 0.001) 
        throw std::invalid_argument("Invalid exposure time (<0.001)");

    if (dispersion < 0.001) 
        throw std::invalid_argument("Invalid dispersion (<0.001)");

    if (nknots < 2 && degree < 0 ) 
        throw std::invalid_argument("Number of knots in spline fitting the " 
                                    "instrument response must be at least 2");

    if (degree < 1 && nknots < 0)
        throw std::invalid_argument("Order of the polynomial fitting the " 
                                    "instrument response must be at least 1");

    nx = cpl_image_get_size_x(spectra);
    ny = cpl_image_get_size_y(spectra);

    /*
     * Find brightest spectrum and duplicate it.
     */
    cpl_size obj_slit = -1;
    {
        cpl_size        maxpos_x, maxpos_y;
        cpl_image *brights = cpl_image_collapse_create(spectra, 1);

        cpl_image_get_maxpos(brights, &maxpos_x, &maxpos_y);

        if(mapped_flat_sed != NULL)
        {
            /* Identify the slit of this object which will be stored in obj_slit*/
            //TODO: Move it away from here. It is also repeated in the science photometric correction
            size_t nslits = cpl_table_get_nrow(objects);
            int   maxobjects = 1;
            char  name[MAX_COLNAME];
            snprintf(name, MAX_COLNAME, "object_%d", maxobjects);
            while (cpl_table_has_column(objects, name)) {
                maxobjects++;
                snprintf(name, MAX_COLNAME, "object_%d", maxobjects);
            }
            maxobjects--;
            obj_slit = -1;
            for (cpl_size i_slit = 0; i_slit < (cpl_size)nslits; i_slit++) {
                for (int i_obj = 1; i_obj <= maxobjects; i_obj++) {
                    snprintf(name, MAX_COLNAME, "row_%d", i_obj);
                    if (cpl_table_is_valid(objects, name, i_slit))
                    {
                        int null;
                        int idx_obj = 
                            cpl_table_get_int(objects, name, i_slit, &null);
                        if (maxpos_y == idx_obj+1)
                            obj_slit = i_slit;
                    }
                }
            }
        }
        
        
        cpl_image_delete(brights);
        spectrum = cpl_image_extract(spectra, 1, maxpos_y, nx, maxpos_y);
    }
    /* Divide the target spectrum by the corresponding slit profile */
    //TODO: This modifies the observed standard spectra. The correction should be applied only at response computation time.
    //It has an impact, since the resp.observed_flux() is saved in the response product.
    cpl_image * spectrum_sedcorr = NULL;
    if(mapped_flat_sed != NULL)
    {
        spectrum_sedcorr = cpl_image_duplicate(spectrum);
        for (int i_pix = 0; i_pix < nx; i_pix++)
        {
            int   null;
            double profile_val  = cpl_image_get(mapped_flat_sed, i_pix+1, obj_slit+1, &null);
            if(profile_val != 0)
                cpl_image_set(spectrum_sedcorr, i_pix+1, 1, 
                        cpl_image_get(spectrum_sedcorr, i_pix+1, 1, &null) / profile_val);
        }
        std::ostringstream norm_key;
        norm_key<< "ESO QC FLAT SED_"<<det_slits[obj_slit].slit_id()<<" NORM";
        flat_sed_norm_factor = 
          cpl_propertylist_get_double(flat_sed_header, norm_key.str().c_str());
        if(cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND)
        {
            cpl_image_delete(spectrum_sedcorr);
            std::string error_msg("Could not find keyword in flat sed: ");
            throw std::runtime_error(error_msg+norm_key.str());
        }
    }

    /* Prepare response computation */
    mosca::spec_std_star std_star(flux_table);
    mosca::response resp;
    mosca::response resp_sedcorr;

    /*
     * Convert standard star spectrum in electrons per second per Angstrom.
     */
    cpl_image_multiply_scalar(spectrum, gain / exptime / dispersion);

    /* Cast to a mosca::spectrum */
    mosca::spectrum std_obs_spectrum(spectrum, startwave, dispersion);

    /* Correct from atmospheric extinction */
    mosca::extinction atm_extinction(ext_table);
    mosca::spectrum std_extcorrect = 
            atm_extinction.correct_spectrum(std_obs_spectrum, airmass);

    /* Compute the normal response */
    resp.compute_response(std_extcorrect, std_star);

    /* Do the same for the SED correction */
    if(mapped_flat_sed != NULL)
    {
        /* Convert standard star spectrum in electrons per second per Angstrom. */
        cpl_image_multiply_scalar(spectrum_sedcorr, gain / exptime / dispersion);

        /* Cast to a mosca::spectrum */
        mosca::spectrum std_obs_spectrum_sedcorr(spectrum_sedcorr, startwave, dispersion);

        /* Correct from atmospheric extinction */
        mosca::spectrum std_extcorrect_sedcorr = 
                atm_extinction.correct_spectrum(std_obs_spectrum_sedcorr, airmass);

        //Compute the response with the SED correction
        resp_sedcorr.compute_response(std_extcorrect_sedcorr, std_star);
    }

    /* Fit the response */
    try
    {
        if(nknots > 0 )
            resp.fit_response_spline(nknots, ignored_waves, ignored_wave_ranges);
        else if(degree > 0 )
            resp.fit_response_pol(degree, ignored_waves, ignored_wave_ranges);
        if(mapped_flat_sed != NULL)
        {
            if(nknots > 0 )
                resp_sedcorr.fit_response_spline(nknots, ignored_waves, ignored_wave_ranges);
            else if(degree > 0 )
                resp_sedcorr.fit_response_pol(degree, ignored_waves, ignored_wave_ranges);
        }
    }
    catch (std::length_error& ex)
    {
        cpl_msg_error(cpl_func, "Too few points in response fitting");
        return NULL;
    }
    
    if(nknots > 0 && resp.nknots_used_response() != (size_t)nknots)
        cpl_msg_warning(cpl_func, "Number of nknots in response fitting too high. "
                "Changed to maximum: %zd", resp.nknots_used_response());
    if(nknots > 0 && resp.nknots_used_efficiency() != (size_t)nknots)
        cpl_msg_warning(cpl_func, "Number of nknots in efficiency fitting too high. "
                "Changed to maximum: %zd", resp.nknots_used_efficiency());
    if(degree > 0 && resp.degree_used_response() != (size_t)degree)
        cpl_msg_warning(cpl_func, "Degree in response fitting too high. "
                "Changed to maximum: %zd", resp.degree_used_response());
    if(degree > 0 && resp.degree_used_efficiency() != (size_t)degree)
        cpl_msg_warning(cpl_func, "Degree in efficiency fitting too high. "
                "Changed to maximum: %zd", resp.degree_used_efficiency());

    cpl_image_delete(spectrum); spectrum = NULL;

    /*
     * Assemble the product spectrophotometric response_table.
     */
    response_table = cpl_table_new(resp.wave_tab().size());

    cpl_table_new_column(response_table, "WAVE", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(response_table, "WAVE", "Angstrom");
    cpl_table_copy_data_double(response_table, "WAVE", &(resp.wave_tab()[0]));

    cpl_table_new_column(response_table, "STD_FLUX", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(response_table, "STD_FLUX", 
                              "10^(-16) erg/(cm^2 s Angstrom)");
    cpl_table_copy_data_double(response_table, "STD_FLUX", &(resp.flux_tab()[0]));

    cpl_table_new_column(response_table, "OBS_FLUX", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(response_table, "OBS_FLUX", "electron/(s Angstrom)");
    cpl_table_copy_data_double(response_table, "OBS_FLUX", &(resp.flux_obs()[0]));

    cpl_table_new_column(response_table, "RAW_EFFICIENCY", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(response_table, "RAW_EFFICIENCY", "electron/photon");
    cpl_table_copy_data_double(response_table, "RAW_EFFICIENCY", &(resp.efficiency_raw()[0]));

    cpl_table_new_column(response_table, "EFFICIENCY", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(response_table, "EFFICIENCY", "electron/photon");
    cpl_table_copy_data_double(response_table, "EFFICIENCY", &(resp.efficiency_fit()[0]));

    cpl_table_new_column(response_table, "RAW_RESPONSE", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(response_table, "RAW_RESPONSE", 
                              "10^(-16) erg/(cm^2 electron)");
    cpl_table_copy_data_double(response_table, "RAW_RESPONSE", &(resp.response_raw()[0]));

    cpl_table_new_column(response_table, "USED_FIT", CPL_TYPE_INT);
    int * used_fit = new int[resp.wave_tab().size()];
    for(size_t i_used = 0; i_used < resp.wave_tab().size(); ++i_used)
    {
        used_fit[i_used] = 1;
        for(size_t i_ignore = 0 ; i_ignore < resp.ignored_waves().size(); ++i_ignore)
        {
            if(resp.wave_tab()[i_used] == resp.ignored_waves()[i_ignore])
                used_fit[i_used] = 0;
        }
    }
    cpl_table_copy_data_int(response_table, "USED_FIT", used_fit);
    delete[] used_fit;

    bool inc_extrapolation = false;
    response_interp = cpl_table_new(resp.wave_obs(inc_extrapolation).size());
    cpl_table_new_column(response_interp, "WAVE", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(response_interp, "WAVE", "Angstrom");
    cpl_table_copy_data_double(response_interp, "WAVE", &(resp.wave_obs(inc_extrapolation)[0]));

    cpl_table_new_column(response_interp, "EFFICIENCY", CPL_TYPE_DOUBLE);
    cpl_table_set_column_unit(response_interp, "EFFICIENCY", "electron/photon");
    cpl_table_copy_data_double(response_interp, "EFFICIENCY", &(resp.efficiency_fit_obs(inc_extrapolation)[0]));

    if(mapped_flat_sed != NULL)
    {
        cpl_table_new_column(response_table, "OBS_FLUX_FFSED", CPL_TYPE_DOUBLE);
        cpl_table_set_column_unit(response_table, "OBS_FLUX_FFSED", "electron/(s Angstrom)");
        cpl_table_copy_data_double(response_table, "OBS_FLUX_FFSED", &(resp_sedcorr.flux_obs()[0]));

        cpl_table_new_column(response_table, "RAW_RESPONSE_FFSED", CPL_TYPE_DOUBLE);
        cpl_table_set_column_unit(response_table, "RAW_RESPONSE_FFSED", 
                                  "10^(-16) erg/(cm^2 electron)");
        cpl_table_copy_data_double(response_table, "RAW_RESPONSE_FFSED", &(resp_sedcorr.response_raw()[0]));

        cpl_table_new_column(response_table, "RESPONSE_FFSED", CPL_TYPE_DOUBLE);
        cpl_table_set_column_unit(response_table, 
                                  "RESPONSE_FFSED", "10^(-16) erg/(cm^2 electron)");
        cpl_table_copy_data_double(response_table, "RESPONSE_FFSED", &(resp_sedcorr.response_fit()[0]));

        cpl_table_new_column(response_interp, "RESPONSE_FFSED", CPL_TYPE_DOUBLE);
        cpl_table_set_column_unit(response_interp, 
                                  "RESPONSE_FFSED", "10^(-16) erg/(cm^2 electron)");
        cpl_table_copy_data_double(response_interp, "RESPONSE_FFSED", &(resp_sedcorr.response_fit_obs(inc_extrapolation)[0]));
    }
    else
    {
        cpl_table_new_column(response_table, "RESPONSE", CPL_TYPE_DOUBLE);
        cpl_table_set_column_unit(response_table, 
                                  "RESPONSE", "10^(-16) erg/(cm^2 electron)");
        cpl_table_copy_data_double(response_table, "RESPONSE", &(resp.response_fit()[0]));

        cpl_table_new_column(response_interp, "RESPONSE", CPL_TYPE_DOUBLE);
        cpl_table_set_column_unit(response_interp, 
                                  "RESPONSE", "10^(-16) erg/(cm^2 electron)");
        cpl_table_copy_data_double(response_interp, "RESPONSE", &(resp.response_fit_obs(inc_extrapolation)[0]));
    }

    return response_table;
}

int vimos_science_correct_flat_sed
(cpl_image *spectra,  cpl_table *objects,
 cpl_image * mapped_flat_sed,
 cpl_propertylist * flat_sed_header,
 cpl_propertylist * specphot_header,
 const vimos::detected_slits& det_slits)
{

    cpl_size nx = cpl_image_get_size_x(spectra);
    cpl_size nslits = cpl_table_get_nrow(objects);
    int   maxobjects = 1;
    char  name[MAX_COLNAME];
    snprintf(name, MAX_COLNAME, "object_%d", maxobjects);
    while (cpl_table_has_column(objects, name)) {
        maxobjects++;
        snprintf(name, MAX_COLNAME, "object_%d", maxobjects);
    }
    maxobjects--;
    //TODO: Move it away from here. It is also repeated in the science photometric correction
    for (cpl_size i_slit = 0; i_slit < nslits; i_slit++) {
        /* Get the proper normalisation factor. The flat sed has to be
         * normalises by the same normalisation factor applied in the
         * creation of specphot.
         * So the original factor applied to the flat sed has to be removed 
         * and then aaply the one in the specphot */
        std::ostringstream norm_key;
        norm_key<< "ESO QC FLAT SED_"<<det_slits[i_slit].slit_id()<<" NORM";
        double flat_sed_norm_orig =
            cpl_propertylist_get_double(flat_sed_header, norm_key.str().c_str());
        double specphot_flat_sed_norm = 
            cpl_propertylist_get_double(specphot_header, "ESO QC RESP FLAT SED_NORM");
        double flat_sed_final_norm_factor = specphot_flat_sed_norm / flat_sed_norm_orig ;

        for (int i_obj = 1; i_obj <= maxobjects; i_obj++) {
            snprintf(name, MAX_COLNAME, "row_%d", i_obj);
            if (cpl_table_is_valid(objects, name, i_slit))
            {
                int null;
                int idx_obj = cpl_table_get_int(objects, name, i_slit, &null);
                /* Divide the target spectrum by the corresponding slit profile */
                for (cpl_size i_pix = 0; i_pix < nx; i_pix++)
                {
                    double profile_val  = cpl_image_get(mapped_flat_sed, i_pix+1, i_slit+1, &null);
                    if(profile_val != 0)
                        cpl_image_set(spectra, i_pix+1, idx_obj+1, 
                                cpl_image_get(spectra, i_pix+1, idx_obj+1, &null) / profile_val *
                                flat_sed_final_norm_factor);
                    else    
                        cpl_image_set(spectra, i_pix+1, idx_obj+1, 0.);
                }
            }
        }
    }
    
    if(cpl_propertylist_get_bool(specphot_header, "ESO QC RESP FLAT SED CORR_SLITWID") &&
       !cpl_propertylist_get_bool(flat_sed_header, "ESO QC FLAT SED CORR_SLITWID"))
        cpl_msg_warning(cpl_func, "The flat SED used to compute the response "
                "includes in its normalisation factors the slit width. "
                "However, the flat SED used to correct the science doesn't. "
                "The flux calibration in this case cannot be performed, "
                "therefore stopping.");

    return 0;
}

int vimos_science_correct_flat_sed_mapped
(cpl_image *mapped_image,  cpl_table *objects,
 cpl_image * mapped_flat_sed,
 cpl_propertylist * flat_sed_header,
 cpl_propertylist * specphot_header,
 const vimos::detected_slits& det_slits)
{

    cpl_size nx = cpl_image_get_size_x(mapped_image);
    cpl_size nslits = cpl_table_get_nrow(objects);
    //TODO: Move it away from here. It is also repeated in the science photometric correction
    for (cpl_size i_slit = 0; i_slit < nslits; i_slit++) {
        /* Get the proper normalisation factor. The flat sed has to be
         * normalises by the same normalisation factor applied in the
         * creation of specphot.
         * So the original factor applied to the flat sed has to be removed 
         * and then aaply the one in the specphot */
        int null;
        std::ostringstream norm_key;
        norm_key<< "ESO QC FLAT SED_"<<det_slits[i_slit].slit_id()<<" NORM";
        double flat_sed_norm_orig =
            cpl_propertylist_get_double(flat_sed_header, norm_key.str().c_str());
        double specphot_flat_sed_norm = 
            cpl_propertylist_get_double(specphot_header, "ESO QC RESP FLAT SED_NORM");
        double flat_sed_final_norm_factor = specphot_flat_sed_norm / flat_sed_norm_orig ;

        int slit_start = cpl_table_get_int(objects, "position", i_slit, &null);
        int slit_length = cpl_table_get_int(objects, "length", i_slit, &null);
        for (int j_pix = slit_start; j_pix < slit_start + slit_length; j_pix++) 
        {
            /* Divide the target spectrum by the corresponding slit profile */
            for (cpl_size i_pix = 0; i_pix < nx; i_pix++)
            {
                double profile_val  = cpl_image_get(mapped_flat_sed, i_pix+1, i_slit+1, &null);
                if(profile_val != 0)
                    cpl_image_set(mapped_image, i_pix+1, j_pix+1, 
                            cpl_image_get(mapped_image, i_pix+1, j_pix+1, &null) / profile_val *
                            flat_sed_final_norm_factor);
                else    
                    cpl_image_set(mapped_image, i_pix+1, j_pix+1, 0.);
            }
        }
    }
    
    if(cpl_propertylist_get_bool(specphot_header, "ESO QC RESP FLAT SED CORR_SLITWID") &&
       !cpl_propertylist_get_bool(flat_sed_header, "ESO QC FLAT SED CORR_SLITWID"))
        cpl_msg_warning(cpl_func, "The flat SED used to compute the response "
                "includes in its normalisation factors the slit width. "
                "However, the flat SED used to correct the science doesn't. "
                "The flux calibration in this case cannot be performed, "
                "therefore stopping.");

    return 0;
}

