/*
 * This file is part of the MOSCA library
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
 * extinction.h
 *
 *  Created on: 2014 4 2
 *      Author: cgarcia
 */

#ifndef EXTINCTION_H_
#define EXTINCTION_H_

#include <cpl.h>
#include <gsl/gsl_interp.h>
#include "spectrum.h"

namespace mosca {

class extinction
{
public:
    extinction(const cpl_table * ext);
    
    extinction(const extinction& other);

    extinction();
 
    spectrum correct_spectrum(const spectrum& spec, double airmass);
    
    virtual ~extinction();

private:

    double eval_at_wave(double wavelength);

    cpl_table * m_extinction;
    
    gsl_interp_accel * m_gsl_acc;

    gsl_interp * m_gsl_interp;

};

} /* namespace mosca */
#endif /* EXTINCTION_H_ */
