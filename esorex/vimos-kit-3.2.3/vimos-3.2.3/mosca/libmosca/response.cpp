/*
 * This file is part of the FORS Data Reduction Pipeline
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
 * response.cpp
 *
 *  Created on: 2014 3 28
 *      Author: cgarcia
 */

#include <functional>
#include <response.h>
#include <vector_utils.h>

mosca::response::response() :
m_nknots_response(0), m_nknots_eff(0), m_degree_response(0), m_degree_eff(0),
m_start_valid_wave_tab(0.0), m_stop_valid_wave_tab(0.0)
{
}

mosca::response::~response()
{
}

void mosca::response::compute_response
(const mosca::spectrum& observed_spectrum, mosca::spec_std_star& std_star)
{
    std::vector<double> std_flux = std_star.flux();
    std::vector<double> std_wave = std_star.wave();
    std::vector<double> std_binsize = std_star.binsize();

    /* Reset all the outputs */
    m_ignored_waves.clear();
    m_wave_tab.clear();
    m_efficiency_raw.clear();
    m_efficiency_fit.clear();
    m_response_fit.clear();
    m_response_raw.clear();
    m_wave_tab.clear();
    m_flux_obs.clear();

    //The response and efficiency are computed in the bins of the
    //tabulated star spectrum
    for(size_t ibin = 0; ibin < std_flux.size() ; ++ibin)
    {
        //TODO: the 0.0026 [cm^2/(erg Angstrom)] is hardcoded, and depends on telescope area 
        //see moses.c  
        //Get the wavelength range
        double wave_bin_start = std_wave[ibin] - std_binsize[ibin] / 2.;
        double wave_bin_end   = std_wave[ibin] + std_binsize[ibin] / 2.;
        
        //Compute the integrated flux of the observed star in this wave range
        double obs_std_flux_bin =
            observed_spectrum.integrate(wave_bin_start, wave_bin_end, true, 0.5) /
               std_binsize[ibin];
        double wavelength = std_wave[ibin];
        double std_flux_bin_phys = std_flux[ibin] * 0.0026 * wavelength; //[photons/s]
        double efficiency = obs_std_flux_bin / std_flux_bin_phys;
        double respon;
        if(obs_std_flux_bin != 0)
            respon =  std_flux[ibin] / obs_std_flux_bin;
        else
            respon = 0;
        
        m_efficiency_raw.push_back(efficiency);
        m_response_raw.push_back(respon);
        m_wave_tab.push_back(wavelength);
        m_wave_tab_bin.push_back(std_binsize[ibin]);
        m_flux_obs.push_back(obs_std_flux_bin);
        m_flux_tab.push_back(std_flux[ibin]);
    }
    
    /* Get the wavelength values of the observed flux */ 
    m_wave_obs = observed_spectrum.wave();
}

void mosca::response::fit_response_spline
(size_t spline_knots,
 const std::vector<double>& waves_to_ignore,
 const std::vector<std::pair<double, double> >& wave_ranges_to_ignore)
{
    double threshold = 0.001;

    m_prepare_fit(waves_to_ignore, wave_ranges_to_ignore);
    
    /* Fitting a cubic spline */
    m_nknots_response = spline_knots;
    m_nknots_eff = spline_knots;
    mosca::vector_cubicspline splfit;
    
    /* Get the knots range */
    double min_x_knot = *std::min_element(m_wave_obs.begin(), m_wave_obs.end());
    double max_x_knot = *std::max_element(m_wave_obs.begin(), m_wave_obs.end());
    
    /* Fit response */
    
    /* Get the threshold in terms of the maximum */
    double max_value = *std::max_element(m_response_fit.begin(), m_response_fit.end());
    double thres_value = threshold * (double)max_value;
    
    /* Create a "mask" of pixels to use */ //TODO: Consolidate the creation of masks
    std::vector<bool> mask_resp;
    std::transform(m_response_fit.begin(), m_response_fit.end(), 
                   std::back_inserter(mask_resp), 
                   std::bind1st(std::less_equal<double>(),thres_value));

    splfit.fit(m_wave_tab, m_response_fit, mask_resp, m_nknots_response,
               min_x_knot, max_x_knot);
    
    /* Get the response fit at the observed binning */
    for(size_t i = 0; i < m_wave_obs.size(); i++)
        m_response_fit_obs.push_back(splfit.eval(m_wave_obs[i]));

    /* Get the response fit at the observed binning, 
     * removing the bins outside the response fit*/
    for(size_t i = 0; i < m_wave_obs_response_valid.size(); i++)
        m_response_fit_obs_response_valid.push_back(splfit.eval(m_wave_obs_response_valid[i]));
    
    /* Fit efficiency */

    /* Get the threshold in terms of the maximum */
    max_value = *std::max_element(m_efficiency_fit.begin(), m_efficiency_fit.end());
    thres_value = threshold * (double)max_value;
    
    /* Create a "mask" of pixels to use */
    std::vector<bool> mask_eff;
    std::transform(m_efficiency_fit.begin(), m_efficiency_fit.end(), 
                   std::back_inserter(mask_eff), 
                   std::bind1st(std::less_equal<double>(),thres_value));
    splfit.fit(m_wave_tab, m_efficiency_fit, mask_eff, m_nknots_eff ,
               min_x_knot, max_x_knot);
    /* Get the efficiency fit at the observed binning */
    for(size_t i = 0; i < m_wave_obs.size(); i++)
        m_efficiency_fit_obs.push_back(splfit.eval(m_wave_obs[i]));
    
    /* Get the efficiency fit at the observed binning, 
     * removing the bins outside the response fit*/
    for(size_t i = 0; i < m_wave_obs_response_valid.size(); i++)
        m_efficiency_fit_obs_response_valid.push_back(splfit.eval(m_wave_obs_response_valid[i]));
}

void mosca::response::fit_response_pol
(size_t pol_degree,
 const std::vector<double>& waves_to_ignore,
 const std::vector<std::pair<double, double> >& wave_ranges_to_ignore)
{
    double threshold = 0.001;

    m_prepare_fit(waves_to_ignore, wave_ranges_to_ignore);

    /* Fitting a polynomial */
    m_degree_response = pol_degree;
    m_degree_eff = pol_degree;
    mosca::vector_polynomial polfit;
    
    /* Fit response */
    /* Get the threshold in terms of the maximum */
    double max_value = *std::max_element(m_response_fit.begin(), m_response_fit.end());
    double thres_value = threshold * (double)max_value;

    /* Create a "mask" of pixels to use */
    std::vector<bool> mask_resp;
    std::transform(m_response_fit.begin(), m_response_fit.end(), 
                   std::back_inserter(mask_resp), 
                   std::bind1st(std::less_equal<double>(),thres_value));
    polfit.fit(m_wave_tab, m_response_fit, mask_resp, m_degree_response);

    /* Get the response fit at the observed binning */
    for(size_t i = 0; i < m_wave_obs.size(); i++)
        m_response_fit_obs.push_back(polfit.eval(m_wave_obs[i]));

    /* Get the response fit at the observed binning, 
     * removing the bins outside the response fit*/
    for(size_t i = 0; i < m_wave_obs_response_valid.size(); i++)
        m_response_fit_obs_response_valid.push_back(polfit.eval(m_wave_obs_response_valid[i]));

    /* Fit efficiency */

    /* Get the threshold in terms of the maximum */
    max_value = *std::max_element(m_efficiency_fit.begin(), m_efficiency_fit.end());
    thres_value = threshold * (double)max_value;

    /* Create a "mask" of pixels to use */
    std::vector<bool> mask_eff;
    std::transform(m_response_fit.begin(), m_response_fit.end(), 
                   std::back_inserter(mask_eff), 
                   std::bind1st(std::less_equal<double>(),thres_value));
    polfit.fit(m_wave_tab, m_efficiency_fit, mask_eff, m_degree_eff);

    /* Get the efficiency fit at the observed binning */
    for(size_t i = 0; i < m_wave_obs.size(); i++)
        m_efficiency_fit_obs.push_back(polfit.eval(m_wave_obs[i]));

    /* Get the efficiency fit at the observed binning, 
     * removing the bins outside the response fit*/
    for(size_t i = 0; i < m_wave_obs_response_valid.size(); i++)
        m_efficiency_fit_obs_response_valid.push_back(polfit.eval(m_wave_obs_response_valid[i]));
}

void mosca::response::m_prepare_fit
(const std::vector<double>& waves_to_ignore,
 const std::vector<std::pair<double, double> >& wave_ranges_to_ignore)
{
    /* Filling the data to fit with the proper data */
    m_response_fit = m_response_raw;
    m_efficiency_fit = m_efficiency_raw;
    std::vector<bool> valid_bins(m_response_raw.size(), true);
    for(size_t ibin = 0; ibin < m_response_raw.size() ; ++ibin)
    {
        //Get the wavelength range
        double wave_bin_start = m_wave_tab[ibin] - m_wave_tab_bin[ibin] / 2.;
        double wave_bin_end   = m_wave_tab[ibin] + m_wave_tab_bin[ibin] / 2.;

        //Check whether this bin should be ignored or not
        bool ignore = false;
        //Is this bin affected by the ignored lines?
        for(size_t iline=0; iline < waves_to_ignore.size(); iline++)
            if(waves_to_ignore[iline]>= wave_bin_start && 
               waves_to_ignore[iline]<= wave_bin_end)
                ignore = true;
        //Is this bin affected by the ignored ranges?
        for(size_t irange=0; irange < wave_ranges_to_ignore.size(); irange++)
        {
            double range_min = std::min(wave_ranges_to_ignore[irange].first,
                                        wave_ranges_to_ignore[irange].second);
            double range_max = std::max(wave_ranges_to_ignore[irange].first,
                                        wave_ranges_to_ignore[irange].second);
            
            if(wave_bin_start <= range_max && wave_bin_end >= range_min)
            {
                ignore = true;
                break;
            }
        }
        //If there is no flux from the standard star, ignore also this bin
        if(m_flux_obs[ibin] == 0)
            ignore = true;
        
        if(ignore)
        {
            m_response_fit[ibin] = 0;
            m_efficiency_fit[ibin] = 0;
            m_ignored_waves.push_back(m_wave_tab[ibin]);
            valid_bins[ibin] = false;
        }
    }

    //Copy the range without extrapolation into m_wave_obs_response_valid 
    std::vector<bool>::iterator first_valid = 
            std::find(valid_bins.begin(), valid_bins.end(), true);
    std::vector<bool>::reverse_iterator last_valid = 
            std::find(valid_bins.rbegin(), valid_bins.rend(), true);
    m_start_valid_wave_tab =
            *(m_wave_tab.begin()+std::distance(valid_bins.begin(),first_valid));
    m_stop_valid_wave_tab =
            *(m_wave_tab.begin()+std::distance(valid_bins.begin(),--last_valid.base()));
    
    for(size_t i = 0; i < m_wave_obs.size(); i++)
        if(m_wave_obs[i] >= m_start_valid_wave_tab &&
           m_wave_obs[i] <= m_stop_valid_wave_tab)
            m_wave_obs_response_valid.push_back(m_wave_obs[i]);
    
}

std::vector<double>& mosca::response::wave_tab()
{
    return m_wave_tab;
}

const std::vector<double>& mosca::response::wave_tab() const
{
    return m_wave_tab;
}


std::vector<double>& mosca::response::efficiency_raw()
{
    return m_efficiency_raw;
}

std::vector<double>& mosca::response::efficiency_fit()
{
    return m_efficiency_fit;
}

std::vector<double>& mosca::response::response_fit()
{
    return m_response_fit;
}

std::vector<double>& mosca::response::response_raw()
{
    return m_response_raw;
}

std::vector<double>& mosca::response::flux_obs()
{
    return m_flux_obs;
}

std::vector<double>& mosca::response::flux_tab()
{
    return m_flux_tab;
}

const std::vector<double>& mosca::response::flux_obs() const
{
	return m_flux_obs;
}

const std::vector<double>& mosca::response::flux_tab() const
{
	return m_flux_tab;
}


std::vector<double>& mosca::response::wave_obs(bool include_extrapolation)
{
    if(include_extrapolation)
        return m_wave_obs;
    else 
        return m_wave_obs_response_valid;
}

std::vector<double>& mosca::response::response_fit_obs(bool include_extrapolation)
{
    if(include_extrapolation)
        return m_response_fit_obs;
    else 
        return m_response_fit_obs_response_valid;
}

std::vector<double>& mosca::response::efficiency_fit_obs(bool include_extrapolation)
{
    if(include_extrapolation)
        return m_efficiency_fit_obs;
    else 
        return m_efficiency_fit_obs_response_valid;
}

std::vector<double>& mosca::response::ignored_waves()
{
    return m_ignored_waves;
}

size_t mosca::response::nknots_used_response() const
{
    return m_nknots_response;
}

size_t mosca::response::nknots_used_efficiency() const
{
    return m_nknots_eff;
}

size_t mosca::response::degree_used_response() const
{
    return m_degree_response;
}

size_t mosca::response::degree_used_efficiency() const
{
    return m_degree_eff;
}

double mosca::response::start_valid_wave_tab() const
{
    return m_start_valid_wave_tab;
}

double mosca::response::stop_valid_wave_tab() const
{
    return m_stop_valid_wave_tab;
}
