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
 * response.h
 *
 *  Created on: 2014 3 28
 *      Author: cgarcia
 */

#ifndef RESPONSE_H_
#define RESPONSE_H_

#include <vector>
#include <utility>
#include "spectrum.h"
#include "spec_std_star.h"

namespace mosca {

class response
{
public:

    typedef std::pair<double, double> wave_range;

    response();

    virtual ~response();
    
    void compute_response
       (const spectrum& observed_spectrum, spec_std_star& std_star);

    void fit_response_spline
        (size_t spline_knots,
         const std::vector<double>& waves_to_ignore  = std::vector<double>(),
         const std::vector<wave_range>& wave_ranges_to_ignore = 
                         std::vector<wave_range>());    

    void fit_response_pol
        (size_t pol_degree,
         const std::vector<double>& waves_to_ignore  = std::vector<double>(),
         const std::vector<wave_range>& wave_ranges_to_ignore = 
                         std::vector<wave_range>());    
    
    std::vector<double>& wave_tab();

    const std::vector<double>& wave_tab() const;

    std::vector<double>& efficiency_raw();

    std::vector<double>& efficiency_fit();

    std::vector<double>& response_fit();

    std::vector<double>& response_raw();

    std::vector<double>& flux_obs();

    std::vector<double>& flux_tab();
    
    const std::vector<double>& flux_obs() const;

    const std::vector<double>& flux_tab() const;

    std::vector<double>& ignored_waves();
    
    //@param include_extrapolation: include points outside of the fitted region
    std::vector<double>& wave_obs(bool include_extrapolation);
    
    std::vector<double>& response_fit_obs(bool include_extrapolation);

    std::vector<double>& efficiency_fit_obs(bool include_extrapolation);

    size_t nknots_used_response() const;

    size_t nknots_used_efficiency() const;

    size_t degree_used_response() const;

    size_t degree_used_efficiency() const;

    //First non-extrapolated wavelength point in wave_tab
    double start_valid_wave_tab() const;

    //Last non-extrapolated wavelength point in wave_tab
    double stop_valid_wave_tab() const;

private:
    
    void m_prepare_fit
        (const std::vector<double>& waves_to_ignore,
         const std::vector<wave_range>& wave_ranges_to_ignore);    

    //The list of wavelengths as tabulated in the standard star 
    std::vector<double> m_wave_tab;

    //The wavelengths width for each bin as tabulated in the standard star 
    std::vector<double> m_wave_tab_bin;
    
    //The computed efficiency at each of the bins in m_wave
    std::vector<double> m_efficiency_raw;

    //The fitted efficiency at each of the bins in m_wave
    std::vector<double> m_efficiency_fit;

    //The computed response at each of the bins in m_wave
    std::vector<double> m_response_raw;

    //The fitted response at each of the bins in m_wave
    std::vector<double> m_response_fit;

    //The integrated observed flux within the bin in m_wave - m_wave_bin, m_wave + m_wave_bin,  
    std::vector<double> m_flux_obs;

    //The tabulated standard star flux at each of the bins in m_wave
    std::vector<double> m_flux_tab;

    //The list of wavelengths that have been ignored for the fit 
    std::vector<double> m_ignored_waves;
    
    //The list of wavelengths of the observed flux 
    std::vector<double> m_wave_obs;
    
    //The fitted response at each of the bins in m_wave_obs
    std::vector<double> m_response_fit_obs;

    //The fitted efficiency at each of the bins in m_wave
    std::vector<double> m_efficiency_fit_obs;

    //The list of wavelengths of the observed flux, removing those wavelengths
    //for which the response fit is being extrapolated
    std::vector<double> m_wave_obs_response_valid;
    
    //The fitted response at each of the bins in m_wave_obs_response_valid
    std::vector<double> m_response_fit_obs_response_valid;

    //The fitted efficiency at each of the bins in m_wave_obs_response_valid
    std::vector<double> m_efficiency_fit_obs_response_valid;

    size_t m_nknots_response;
    
    size_t m_nknots_eff;

    size_t m_degree_response;
    
    size_t m_degree_eff;

    double m_start_valid_wave_tab;

    double m_stop_valid_wave_tab;
};

} /* namespace mosca */
#endif /* RESPONSE_H_ */
