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
 * spec_std_star.cpp
 *
 *  Created on: 2014 3 28
 *      Author: cgarcia
 */

#include <vector>
#include <spec_std_star.h>

mosca::spec_std_star::spec_std_star(cpl_table * tabulated_data)
{
    m_std_star_data = cpl_table_duplicate(tabulated_data);
    cpl_table_cast_column(m_std_star_data, "WAVE", "WAVE_D", CPL_TYPE_DOUBLE);
    cpl_table_cast_column(m_std_star_data, "FLUX", "FLUX_D", CPL_TYPE_DOUBLE);
    cpl_table_cast_column(m_std_star_data, "BIN", "BIN_D", CPL_TYPE_DOUBLE);
}

mosca::spec_std_star::spec_std_star(const mosca::spec_std_star& other) :
 m_std_star_data(NULL)
{
    m_std_star_data = cpl_table_duplicate(other.m_std_star_data);
}

mosca::spec_std_star::spec_std_star() :
 m_std_star_data(NULL)
{
}

mosca::spec_std_star::~spec_std_star()
{
    if(m_std_star_data != NULL)
        cpl_table_delete(m_std_star_data);
}

mosca::spec_std_star::spec_std_star(spec_std_star& other) :
         m_std_star_data(NULL)
{
    if(other.m_std_star_data != NULL)
        m_std_star_data = cpl_table_duplicate(other.m_std_star_data);
}

std::vector<double> mosca::spec_std_star::flux() const 
{
    std::vector<double> flux_vec;
    double * flux_p =
        cpl_table_get_data_double(m_std_star_data, "FLUX_D");
    cpl_size size = cpl_table_get_nrow(m_std_star_data);
    flux_vec.insert(flux_vec.end(), flux_p, flux_p + size);
    return flux_vec;
}

std::vector<double> mosca::spec_std_star::wave() const 
{
    std::vector<double> wave_vec;
    double * wave_p =
        cpl_table_get_data_double(m_std_star_data, "WAVE_D");
    cpl_size size = cpl_table_get_nrow(m_std_star_data);
    wave_vec.insert(wave_vec.end(), wave_p, wave_p + size);
    return wave_vec;
}

std::vector<double> mosca::spec_std_star::binsize() const
{
    std::vector<double> binsize_vec;
    double * binsize_p =
        cpl_table_get_data_double(m_std_star_data, "BIN_D");
    cpl_size size = cpl_table_get_nrow(m_std_star_data);
    binsize_vec.insert(binsize_vec.end(), binsize_p, binsize_p + size);
    return binsize_vec;
}


