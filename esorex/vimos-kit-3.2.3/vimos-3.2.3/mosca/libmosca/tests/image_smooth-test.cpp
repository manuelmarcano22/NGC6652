/* $Id: image_smooth-test.cpp,v 1.15 2013-08-08 19:37:40 cgarcia Exp $
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
 * $Date: 2013-08-08 19:37:40 $
 * $Revision: 1.15 $
 * $Name: not supported by cvs2svn $
 */

#include "config.h"

#if defined HAVE_BOOST_UNIT_TEST_FRAMEWORK && HAVE_CXX11

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE image_smooth-test
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
#include "image_smooth.h"
#include "type_traits.h"

BOOST_AUTO_TEST_SUITE(image_smooth_types)

typedef boost::mpl::list<float, double> test_types;
typedef std::minstd_rand0 rnd_generator; //std::default_random_engine with a modern C++11 (not available in psdlin1) 

BOOST_AUTO_TEST_CASE_TEMPLATE(image_smooth_types, T, test_types)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    cpl_size   nx = 100;
    cpl_size   ny = 100;
    
    rnd_generator rnd; 
    std::normal_distribution<T> gauss(5.0,2.0);

    for(cpl_size half_width = 1; half_width < nx / 2 - 1; half_width += 5)
    {

        mosca::image image(nx, ny, type);
        
        //Fill the image with random data
        std::generate(image.get_data<T>(), image.get_data<T>() + nx*ny,
                      std::bind(gauss, rnd));
        mosca::image_smooth_1d_median<T>(image, half_width, 
                                      mosca::DISPERSION_AXIS);
    }
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

typedef boost::mpl::list<float> test_default_type;

BOOST_AUTO_TEST_CASE_TEMPLATE(image_smooth_default_type, T, test_default_type)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    cpl_size   nx = 100;
    cpl_size   ny = 100;
    
    rnd_generator rnd;
    std::normal_distribution<T> gauss(5.0,2.0);

    for(cpl_size half_width = 1; half_width < nx / 2 - 1; half_width += 5)
    {

        mosca::image image(nx, ny, type);
        
        //Fill the image with random data
        std::generate(image.get_data<T>(), image.get_data<T>() + nx*ny,
                      std::bind(gauss,rnd));
        mosca::image_smooth_1d_median<T>(image, half_width, 
                                      mosca::DISPERSION_AXIS);
    }
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(image_smooth_errors)

BOOST_AUTO_TEST_CASE(wrong_half_width_size_disp)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_size   nx = 100;
    cpl_size   ny = 1000;
    
    mosca::image image(nx, ny, CPL_TYPE_FLOAT);
    BOOST_REQUIRE_THROW(mosca::image_smooth_1d_median<float>(
            image, 60, mosca::DISPERSION_AXIS), std::out_of_range);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(wrong_half_width_size_spa)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_size   nx = 1000;
    cpl_size   ny = 100;
    
    mosca::image image(nx, ny, CPL_TYPE_FLOAT);
    BOOST_REQUIRE_THROW(mosca::image_smooth_1d_median<float>(
            image, 60, mosca::SPATIAL_AXIS), std::out_of_range);
    cpl_end();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(image_smooth_axes)

typedef std::minstd_rand0 rnd_generator; //std::default_random_engine with a modern C++11 (not available in psdlin1) 

BOOST_AUTO_TEST_CASE(smooth_disp_spat_x_y)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    cpl_size   nx = 100;
    cpl_size   ny = 100;
    mosca::axis axes[] = {mosca::DISPERSION_AXIS, mosca::SPATIAL_AXIS, 
            mosca::X_AXIS, mosca::Y_AXIS};
    
    rnd_generator rnd;
    std::normal_distribution<double> gauss(50.0,2.0);

    for(mosca::axis axis : axes)
    {
        for(cpl_size kernel_size = 1; kernel_size*2  < nx; kernel_size+=5)
        {
            mosca::image image(nx, ny, CPL_TYPE_FLOAT);
            //Fill the image with random data
            std::generate(image.get_data<float>(), image.get_data<float>() + nx*ny,
                          std::bind(gauss,rnd));
            mosca::image_smooth_1d_median<float>(image, kernel_size, axis);
        }
    }
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_AUTO_TEST_CASE(smooth_xy)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    cpl_size   nx = 100;
    cpl_size   ny = 100;
    

    for(cpl_size kernel_size = 1; kernel_size*2 < nx; kernel_size+=5)
    {
        mosca::image image(nx, ny, CPL_TYPE_FLOAT);
        //Fill the image with random data
        mosca::image_smooth_1d_median<float>(image, kernel_size, 
                                      mosca::X_AXIS);
        mosca::image_smooth_1d_median<float>(image, kernel_size, 
                                      mosca::Y_AXIS);
    }
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_AUTO_TEST_SUITE_END()

struct random_image_fixture
{
    typedef std::minstd_rand0 rnd_generator; //std::default_random_engine with a modern C++11 (not available in psdlin1) 

    random_image_fixture()
    {
        BOOST_TEST_MESSAGE( "setup fixture" ); 
        cpl_init(CPL_INIT_DEFAULT);
        half_width = 5;
        nx = 1000;
        ny = 1000;
        image.reset(new mosca::image(nx, ny, CPL_TYPE_FLOAT));
        rnd_generator rnd;
        std::normal_distribution<double> gauss(5.0,2.0);
        std::generate(image->get_data<float>(), image->get_data<float>() + nx*ny,
                      std::bind(gauss,rnd));
        cpl_image * image_tmp =
                cpl_image_wrap_float(nx, ny, image->get_data<float>());
        image_cpl = cpl_image_duplicate(image_tmp);
        cpl_image_unwrap(image_tmp);
    }
    
    ~random_image_fixture()
    {
        BOOST_TEST_MESSAGE( "teardown fixture" );
        cpl_image_delete(image_cpl);
        cpl_end();
    }
    
    cpl_image *                   image_cpl;
    std::unique_ptr<mosca::image> image;
    cpl_size                      half_width;
    cpl_size                      nx;
    cpl_size                      ny;
};

BOOST_FIXTURE_TEST_SUITE(image_smooth_benchmark, random_image_fixture)


BOOST_AUTO_TEST_CASE(benchmark_mosca_smooth_x)
{    
    mosca::image_smooth_1d_median<float>(*image, half_width,mosca::X_AXIS);
}

BOOST_AUTO_TEST_CASE(benchmark_mosca_smooth_y)
{
    
    mosca::image_smooth_1d_median<float>(*image, half_width,mosca::Y_AXIS);
}

BOOST_AUTO_TEST_SUITE_END()

#else

int main(void)
{
    return 0;
}

#endif
