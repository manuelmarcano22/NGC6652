/* 
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


#include "config.h"

#if defined HAVE_BOOST_UNIT_TEST_FRAMEWORK && HAVE_CXX11

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE image_normalisation-test
#define _GLIBCXX_USE_NANOSLEEP
#include <random>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include "cpl.h"
#include "image_normalisation.h"

BOOST_AUTO_TEST_SUITE(image_normalisation_exceptions)

BOOST_AUTO_TEST_CASE(empty_image)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    
    size_t nx = 100;
    size_t ny = 10;
    
    mosca::image empty(nx, ny);
    mosca::image weight(nx, ny);
    std::vector<double> spa_profile;
    std::vector<double> spec_profile;
    
    BOOST_CHECK_NO_THROW(image_normalise(empty, weight, 0, 0, 0, 0, 0.,
                                      spa_profile, spec_profile));
    
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_AUTO_TEST_CASE(size_mismatch)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    
    size_t nx = 100;
    size_t ny = 10;
    
    mosca::image empty(nx, ny);
    mosca::image weight(ny, nx);
    std::vector<double> spa_profile;
    std::vector<double> spec_profile;
    
    BOOST_CHECK_THROW(image_normalise(empty, weight, 0, 0, 0, 0, 0., 
                                      spa_profile, spec_profile), 
                      std::invalid_argument);
    
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(image_normalisation_templates)

typedef boost::mpl::list<float, double, int> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(smooth, T, test_types)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    
    size_t nx = 100;
    size_t ny = 10;
    T flux = 1000;
    cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    
    mosca::image flat_image(nx, ny, type);
    mosca::image weight(nx, ny,type);
    std::vector<T> spa_profile;
    std::vector<T> spec_profile;
    
    std::fill(flat_image.get_data<T>(), flat_image.get_data<T>() + nx*ny, flux);
    std::fill(weight.get_data<T>(), weight.get_data<T>() + nx*ny, 1);
    
    //Loop on various smoothing values in both directions
    //TODO: If spa_smooth_radius is larger than half the size in spatial dir,
    //then it fails the test. Check it
    for (int spa_smooth_radius = 0; spa_smooth_radius <= 3; spa_smooth_radius++)
        for (int disp_smooth_radius = 0; disp_smooth_radius <= 12; disp_smooth_radius++)
        {
            mosca::image smoothed = 
                image_normalise(flat_image, weight, spa_smooth_radius, 
                                disp_smooth_radius, 0, 0, 0., 
                                spa_profile, spec_profile);
            //When doing smoothing and the input is already 
            //flat we should be getting the same thing
            //Actually even if not doing normalisation, 
            //should get something flat everywhere
            BOOST_CHECK_EQUAL_COLLECTIONS(flat_image.get_data<T>(), 
                                          flat_image.get_data<T>() + nx*ny, 
                                          smoothed.get_data<T>(), 
                                          smoothed.get_data<T>() + nx*ny);
        }
    
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

typedef boost::mpl::list<float, double> test_types_fp;
BOOST_TEST_DECORATOR( * boost::unit_test::tolerance(1e-6))
BOOST_TEST_DECORATOR( * boost::unit_test::tolerance((float)1e-6))

BOOST_AUTO_TEST_CASE_TEMPLATE(fit, T, test_types_fp)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    
    size_t nx = 100;
    size_t ny = 10;
    T flux = 10;
    cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    
    mosca::image flat_image(nx, ny, type);
    mosca::image weight(nx, ny,type);
    std::vector<T> spa_profile;
    std::vector<T> spec_profile;
    
    std::fill(flat_image.get_data<T>(), flat_image.get_data<T>() + nx*ny, flux);
    std::fill(weight.get_data<T>(), weight.get_data<T>() + nx * ny, 1.0);
    
    //Loop on various fitting values in both directions
    for (int spa_fit_degree = 0; spa_fit_degree <= 3; spa_fit_degree++)
        for (int disp_fit_nknots = 2; disp_fit_nknots <= 10; disp_fit_nknots++)
        {
            mosca::image smoothed = 
                image_normalise(flat_image, weight, 0, 0, 
                                spa_fit_degree, disp_fit_nknots, 0., 
                                spa_profile, spec_profile);
            //Check for close enough fits
            std::vector<T> differences;
            std::transform (flat_image.get_data<T>(),
                            flat_image.get_data<T>() + nx*ny,
                            smoothed.get_data<T>(),
                            std::back_inserter(differences), std::minus<T>());
            std::vector<T> zero(differences.size());
            BOOST_TEST(differences == zero, boost::test_tools::per_element());
        }
    
    //Check for case disp_fit_nknots==1 which is not support by gsl:
    BOOST_CHECK_THROW(image_normalise(flat_image, weight, 0, 0, 
                                      0, 1, 0., 
                                      spa_profile, spec_profile), std::exception);
    
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_TEST_DECORATOR( * boost::unit_test::tolerance(1e-4))
BOOST_TEST_DECORATOR( * boost::unit_test::tolerance((float)1e-4))

BOOST_AUTO_TEST_CASE_TEMPLATE(disp_y, T, test_types_fp)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    
    size_t ny = 100;
    size_t nx = 10;
    T av_flux = 100;
    T x_span_factor = 0.1;
    T y_span_factor = 0.1;
    cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    
    mosca::image flat_image(nx, ny, type, mosca::Y_AXIS);
    mosca::image weight(nx, ny, type, mosca::Y_AXIS);
    std::vector<T> spa_profile;
    std::vector<T> spec_profile;
    
    std::fill(weight.get_data<T>(), weight.get_data<T>() + nx*ny, 1);
    
    size_t ipix = 0;
    std::generate(flat_image.get_data<T>(), flat_image.get_data<T>() + nx*ny, 
                  [nx, ny, av_flux, x_span_factor, y_span_factor, &ipix]()
                  {//Linear gradient
                      size_t i = ipix % nx;
                      size_t j = ipix / nx;
                      T fx = std::sqrt(av_flux) * (1 + x_span_factor * ((i-nx/2.)/nx));
                      T fy = std::sqrt(av_flux) * (1 + y_span_factor * ((j-ny/2.)/ny));
                      ipix++;
                      return fx * fy ;
                  });
    //Loop on various normalisation parameters
    for (int spa_smooth_radius = 0; spa_smooth_radius <= 3; spa_smooth_radius++)
        for (int disp_smooth_radius = 0; disp_smooth_radius <= 3; disp_smooth_radius++)
            for (int spa_fit_degree = 1; spa_fit_degree <= 3; spa_fit_degree++)
                for (int disp_fit_nknots = 2; disp_fit_nknots <= 4; disp_fit_nknots++)
                {
                    mosca::image smoothed = 
                            image_normalise(flat_image, weight, 
                                    spa_smooth_radius, disp_smooth_radius, 
                                    spa_fit_degree, disp_fit_nknots, 0., 
                                    spa_profile, spec_profile);
                    //Check for close enough fits
                    std::vector<T> differences;
                    std::transform (flat_image.get_data<T>(),
                                    flat_image.get_data<T>() + nx*ny,
                                    smoothed.get_data<T>(),
                                    std::back_inserter(differences), std::minus<T>());
                    std::vector<T> zero(differences.size());
                    BOOST_TEST(differences == zero, boost::test_tools::per_element());
                }
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_TEST_DECORATOR( * boost::unit_test::tolerance(1e-4))
BOOST_TEST_DECORATOR( * boost::unit_test::tolerance((float)1e-4))

BOOST_AUTO_TEST_CASE_TEMPLATE(gradient, T, test_types_fp)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    
    size_t nx = 100;
    size_t ny = 10;
    T av_flux = 100;
    cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    
    mosca::image flat_image(nx, ny, type);
    mosca::image weight(nx, ny,type);
    std::vector<T> spa_profile;
    std::vector<T> spec_profile;
    
    std::fill(weight.get_data<T>(), weight.get_data<T>() + nx*ny, 1);
    
    //Loop on various gradient values in both directions
    for (T x_span_factor = 0.1; x_span_factor <= 0.3; x_span_factor+=0.1)
        for (T y_span_factor = 0.1; y_span_factor <= 0.3; y_span_factor+=0.1)
        {
            size_t ipix = 0;
            std::generate(flat_image.get_data<T>(), flat_image.get_data<T>() + nx*ny, 
                          [nx, ny, av_flux, x_span_factor, y_span_factor, &ipix]()
                          {//Linear gradient
                            size_t i = ipix % nx;
                            size_t j = ipix / nx;
                            T fx = std::sqrt(av_flux) * (1 + x_span_factor * ((i-nx/2.)/nx));
                            T fy = std::sqrt(av_flux) * (1 + y_span_factor * ((j-ny/2.)/ny));
                            ipix++;
                            return fx * fy ;
                          });
            //Loop on various normalisation parameters
            for (int spa_smooth_radius = 1; spa_smooth_radius <= 3; spa_smooth_radius++)
                for (int disp_smooth_radius = 1; disp_smooth_radius <= 5; disp_smooth_radius++)
                    for (int spa_fit_degree = 1; spa_fit_degree <= 3; spa_fit_degree++)
                        for (int disp_fit_nknots = 2; disp_fit_nknots <= 6; disp_fit_nknots++)
                        {
                            mosca::image smoothed = 
                                    image_normalise(flat_image, weight, 
                                            spa_smooth_radius, disp_smooth_radius, 
                                            spa_fit_degree, disp_fit_nknots, 0., 
                                            spa_profile, spec_profile);
                            //Check for close enough fits
                            std::vector<T> differences;
                            std::transform (flat_image.get_data<T>(),
                                            flat_image.get_data<T>() + nx*ny,
                                            smoothed.get_data<T>(),
                                            std::back_inserter(differences), std::minus<T>());
                            std::vector<T> zero(differences.size());
                            BOOST_TEST(differences == zero, boost::test_tools::per_element());
                        }
        }
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(image_normalisation_weights)

typedef boost::mpl::list<float, double> test_types_fp;

BOOST_AUTO_TEST_CASE_TEMPLATE(weight_flat_smooth, T, test_types_fp)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    
    size_t nx = 100;
    size_t ny = 10;
    T flux = 10;
    cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    
    mosca::image flat_image(nx, ny, type);
    mosca::image weight(nx, ny,type);
    std::vector<T> spa_profile;
    std::vector<T> spec_profile;
    
    std::fill(flat_image.get_data<T>(), flat_image.get_data<T>() + nx*ny, flux);
    std::fill(weight.get_data<T>(), weight.get_data<T>() + nx*ny, 1);
    //Set the outer borders of the weights to 0 and the image to something bad 
    std::fill(weight.get_data<T>(), weight.get_data<T>() + nx, 0);
    std::fill(weight.get_data<T>() + nx*(ny-1), weight.get_data<T>() + nx*ny, 0);
    std::fill(flat_image.get_data<T>(), flat_image.get_data<T>() + nx, 100);
    std::fill(flat_image.get_data<T>() + nx*(ny-1), flat_image.get_data<T>() + nx*ny, 50);
    for(int iy=1;iy<=ny;iy++)
    {
        cpl_image_set(weight.get_cpl_image(), 1, iy, 0.);
        cpl_image_set(weight.get_cpl_image(), nx, iy, 0.);
        cpl_image_set(flat_image.get_cpl_image(), 1, iy, 20.);
        cpl_image_set(flat_image.get_cpl_image(), nx, iy, 40.);
    }
    //One pixel in the corner is actually valid...
    cpl_image_set(weight.get_cpl_image(), nx, ny, 1.);
    cpl_image_set(flat_image.get_cpl_image(), nx, ny, flux);
    
    //Loop on various smoothing values in both directions
    //TODO: If spa_smooth_radius is larger than half the size in spatial dir,
    //then it fails the test. Check it
    for (int spa_smooth_radius = 0; spa_smooth_radius <= 3; spa_smooth_radius++)
        for (int disp_smooth_radius = 0; disp_smooth_radius <= 12; disp_smooth_radius++)
        {
            mosca::image smoothed = 
                image_normalise(flat_image, weight, spa_smooth_radius, 
                                disp_smooth_radius, 0, 0, 0., 
                                spa_profile, spec_profile);
            //When doing smoothing and the input is already 
            //flat we should be getting the same thing
            //Actually even if not doing normalisation, 
            //should get something flat everywhere
            //Here we check only the regions with weights != 0
            for(int iy=1;iy<ny-1;iy++)
            {
                BOOST_CHECK_EQUAL_COLLECTIONS(flat_image.get_data<T>() + 1 + iy * nx, 
                                              flat_image.get_data<T>() + iy * nx + nx -1, 
                                              smoothed.get_data<T>() + 1 + iy * nx, 
                                              smoothed.get_data<T>() + iy * nx + nx -1);
            }
        }
    
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_TEST_DECORATOR( * boost::unit_test::tolerance(0.7))
BOOST_TEST_DECORATOR( * boost::unit_test::tolerance((float)0.7))

BOOST_AUTO_TEST_CASE_TEMPLATE(random_weight_gradient_smooth, T, test_types_fp)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_errorstate prestate =  cpl_errorstate_get();
    
    size_t nx = 100;
    size_t ny = 10;
    T av_flux = 100;
    T random_stddev  = 1000;
    float bad_frac = 0.20;
    cpl_type type = mosca::type_trait<T>::cpl_eq_type;
    
    mosca::image flat_image(nx, ny, type);
    mosca::image weight(nx, ny,type);
    std::vector<T> spa_profile;
    std::vector<T> spec_profile;
    
    std::fill(weight.get_data<T>(), weight.get_data<T>() + nx*ny, 1);
    
    //Loop on various gradient values in both directions
    for (T x_span_factor = 0.0; x_span_factor <= 0.3; x_span_factor+=0.1)
        for (T y_span_factor = 0.0; y_span_factor <= 0.3; y_span_factor+=0.1)
        {
            size_t ipix = 0;
            std::generate(flat_image.get_data<T>(), flat_image.get_data<T>() + nx*ny, 
                          [nx, ny, av_flux, x_span_factor, y_span_factor, &ipix]()
                          {//Linear gradient
                            size_t i = ipix % nx;
                            size_t j = ipix / nx;
                            T fx = std::sqrt(av_flux) * (1 + x_span_factor * ((i-nx/2.)/nx));
                            T fy = std::sqrt(av_flux) * (1 + y_span_factor * ((j-ny/2.)/ny));
                            ipix++;
                            return fx * fy;
                          });
            //Distribute randomly invalid pixels
            std::fill(weight.get_data<T>(), weight.get_data<T>() + nx*ny, 1);
            std::default_random_engine rnd;
            std::uniform_int_distribution<> xpix_dist(1, nx);        
            std::uniform_int_distribution<> ypix_dist(1, ny);        
            std::normal_distribution<> val_dist(av_flux, random_stddev);
            for(size_t ibad = 0; (float)ibad < nx*ny*bad_frac; ++ibad)
            {
                size_t xpix = xpix_dist(rnd); 
                size_t ypix = ypix_dist(rnd);
                T im_val = val_dist(rnd);
                cpl_image_set(weight.get_cpl_image(), xpix, ypix, 0.);
                cpl_image_set(flat_image.get_cpl_image(), xpix, ypix, im_val);
            }
            //Loop on various normalisation parameters
            //Since this is not a flat image, it doesn't make much sense
            //to try with the smoothing parameters, since they are not 
            //really able to cope with the gradient.
            for (int spa_fit_degree = 1; spa_fit_degree <= 2; spa_fit_degree++)
                for (int disp_fit_nknots = 3; disp_fit_nknots <= 6; disp_fit_nknots++)
                {
                    mosca::image smoothed = 
                            image_normalise(flat_image, weight, 
                                    0, 0, 
                                    spa_fit_degree, disp_fit_nknots, 0., 
                                    spa_profile, spec_profile);
                    //Test in the pixels where the mask is valid
                    mosca::image flat_w = flat_image;
                    std::transform (flat_image.get_data<T>(), flat_image.get_data<T>() + nx*ny,
                                    weight.get_data<T>(),
                                    flat_w.get_data<T>(), std::multiplies<T>());
                    mosca::image smooth_w = smoothed;
                    std::transform (smooth_w.get_data<T>(), smooth_w.get_data<T>() + nx*ny,
                                    weight.get_data<T>(),
                                    smooth_w.get_data<T>(), std::multiplies<T>());
                    //Check for close enough fits. The tolerance is
                    //actually quite high 0.7, but that's the best that can
                    //be done with 20% of the points masked
                    std::vector<T> differences;
                    std::transform (flat_w.get_data<T>(),
                                    flat_w.get_data<T>() + nx*ny,
                                    smooth_w.get_data<T>(),
                                    std::back_inserter(differences), std::minus<T>());
                    std::vector<T> zero(differences.size());
                    BOOST_TEST(differences == zero, boost::test_tools::per_element());
                }
        }
    BOOST_REQUIRE_EQUAL(prestate, cpl_errorstate_get() );
    cpl_end();
}

BOOST_AUTO_TEST_SUITE_END()

#else

int main(void)
{
    return 0;
}

#endif
