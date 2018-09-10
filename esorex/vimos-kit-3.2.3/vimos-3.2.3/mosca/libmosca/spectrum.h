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
 * spectrum.h
 *
 *  Created on: 2014 3 28
 *      Author: cgarcia
 */

#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include <vector>
#include <string>
#include <cpl.h>
#include <gsl/gsl_interp.h>

namespace mosca {

class spectrum
{
public:

    //TODO: Add the error of the spectrum?
    spectrum(const cpl_image * spec, double start_wave, double step_wave);

    spectrum(std::vector<double>& flux, std::vector<double>& wave);

    spectrum(const spectrum& other);
    
    spectrum();
    
    virtual ~spectrum();
    
    std::vector<double> flux() const;

    std::vector<double> wave() const;

    //Returns the flux at a given wave, performing linear interpolation if needed
    double flux_at_wave(double wave) const;
    
    spectrum rebin
      (double start_wave, double end_wave, double step_wave);

    //Returns the integrated flux in a wave range. It uses linear interpolation
    //If ignore_neg_flux is true, all the bins in which the standard star
    //observation has no flux, will be ignored 
    //The ignore_threshold means that if a bin has less than that fraction with
    //positive values, then it is ignored
    double integrate(double start_wave, double end_wave, 
                     bool ignore_neg_flux = false,
                     float ignore_threshold = 0) const;
    
    bool regular_sampling();
    
    //0 if not regular sampling
    double dispersion();
    double average_dispersion();
    
    std::string units();
    
    
private:
    
    std::vector<double> m_flux;

    std::vector<double> m_wave;
    
    //Useful for caching
    mutable std::vector<double> m_flux_nonzero;

    //Useful for caching
    mutable std::vector<double> m_wave_nonzero;
    
    //Useful for caching
    mutable gsl_interp_accel * m_gsl_acc;

    //Useful for caching
    mutable gsl_interp * m_gsl_interp;

    void m_create_filtered_flux() const;

};

} /* namespace mosca */
#endif /* SPECTRUM_H_ */
