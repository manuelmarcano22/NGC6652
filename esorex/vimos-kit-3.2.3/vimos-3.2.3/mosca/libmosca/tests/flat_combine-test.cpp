/* $Id: flat_combine-test.cpp,v 1.5 2013-08-14 13:36:49 cgarcia Exp $
 *
 * This file is part of the MOSCA library
 * Copyright (C) 2013 European Southern Observatory
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013-08-14 13:36:49 $
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#include "config.h"

#if defined HAVE_BOOST_UNIT_TEST_FRAMEWORK && HAVE_CXX11

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE flat_combine-test
#define _GLIBCXX_USE_NANOSLEEP
#include <memory>
#include <random>
#include <functional>
#include <algorithm>
#include <functional>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <cpl.h>
#include "flat_combine.h"
#include "reduce_method.h"
#include "type_traits.h"

BOOST_AUTO_TEST_SUITE(flat_combine_types)

typedef boost::mpl::list<float, double> test_types;
typedef std::default_random_engine rnd_generator;

BOOST_AUTO_TEST_CASE_TEMPLATE(flat_combine_types, T, test_types)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    //cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    cpl_size   nx = 1000;
    cpl_size   ny = 500;
    
    
    std::vector<mosca::image> flats;
    flats.push_back(mosca::image(nx, ny));
    flats.push_back(mosca::image(nx, ny));
    std::vector<mosca::calibrated_slit> slits;
    mosca::wavelength_calibration wave_cal;
    double smooth_size =5;
    mosca::reduce_median reduce_method;
    mosca::grism_config grism_cfg(1., 4000., 7000., 5500.);
    std::auto_ptr<mosca::image> master_flat = 
        mosca::flat_combine<T, mosca::reduce_median>
            (flats, slits, wave_cal, grism_cfg, smooth_size, reduce_method);
    
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

typedef boost::mpl::list<float> test_default_type;


BOOST_AUTO_TEST_SUITE_END()

#else

int main(void)
{
    return 0;
}

#endif
