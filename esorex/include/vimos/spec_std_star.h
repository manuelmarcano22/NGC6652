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
 * spec_std_star.h
 *
 *  Created on: 2014 3 28
 *      Author: cgarcia
 */

#ifndef SPEC_STD_STAR_H_
#define SPEC_STD_STAR_H_

#include "cpl.h"

namespace mosca {

class spec_std_star
{
public:

    spec_std_star(cpl_table * tabulated_data);

    spec_std_star(const spec_std_star& other);

    spec_std_star(spec_std_star& other);

    spec_std_star();

    virtual ~spec_std_star();
    
    std::vector<double> flux() const; 

    std::vector<double> wave() const; 

    std::vector<double> binsize() const; 

private:

    cpl_table * m_std_star_data;
};

} /* namespace mosca */
#endif /* SPEC_STD_STAR_H_ */
