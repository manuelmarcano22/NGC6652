/*
 * This file is part of the CASU Pipeline utilities
 * Copyright (C) 2015,2016 European Southern Observatory
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include "../casu_utilfunctions.h"
#include "../../hdrl_catalogue.h"
#include "../imcore.h"
#include <cpl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NTEST 10
#define NCOLS 63

static void fill_plist(cpl_propertylist * pl)
{
    cpl_propertylist_update_string(pl,"CTYPE1","RA---TAN");
    cpl_propertylist_update_string(pl,"CTYPE2","DEC--TAN");
    cpl_propertylist_update_double(pl,"CRVAL1",30.0);
    cpl_propertylist_update_double(pl,"CRVAL2",12.0);
    cpl_propertylist_update_double(pl,"CRPIX1",512.0);
    cpl_propertylist_update_double(pl,"CRPIX2",512.0);
    cpl_propertylist_update_double(pl,"CD1_1",-1.0/3600);
    cpl_propertylist_update_double(pl,"CD1_2",0.0);
    cpl_propertylist_update_double(pl,"CD2_1",0.0);
    cpl_propertylist_update_double(pl,"CD2_2",1.0/3600);
}

static
cpl_image * create_gauss(double dx, double dy, double fwhm)
{
    /* sample gaussian on large grid and sum it down to requested fwhm */
    const size_t factor = 16;
    double sigma = fwhm * factor / (2 * sqrt(2 * log(2)));
    size_t nx = fwhm * 20;
    size_t ny = fwhm * 20;
    size_t nnx = nx * factor;
    size_t nny = ny * factor;
    dx *= factor;
    dy *= factor;
    cpl_image * g = cpl_image_new(nnx, nny, CPL_TYPE_FLOAT);
    cpl_image_fill_gaussian(g, nnx / 2 + dx, nny / 2 + dx,
                            2.0 * CPL_MATH_PI * sigma * sigma,
                            sigma, sigma);
    cpl_image * r = cpl_image_new(nx, ny, CPL_TYPE_FLOAT);
    float * dg = cpl_image_get_data_float(g);
    float * dr = cpl_image_get_data_float(r);
    for (size_t y = 0; y < ny; y++) {
        for (size_t x = 0; x < nx; x++) {
            for (size_t i = y * factor; i < y * factor + factor; i++) {
                for (size_t j = x * factor; j < x * factor + factor; j++) {
                    dr[y * nx + x] += dg[i * nnx + j];
                }
            }
        }
    }
    for (size_t y = 0; y < nny; y++) {
        for (size_t x = 0; x < nnx; x++) {
            dr[y / factor * nx + x / factor] += dg[y * nx + x];
        }
    }
    cpl_image_divide_scalar(r, factor * factor);

    cpl_image_delete(g);
    return r;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Check imcore_conf() in various conditions
  @return cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrl_catalogue_imcore_conf_basic(void)
{
    int retval,i,nrow,ncol,nl;
    cpl_image *im,*bkg,*cnf;
    cpl_table *tab;
    double sigma=2.0,norm2,tot[NTEST],sky=500.0,diff;
    double xpos[] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,
                     1000.0};
    double ypos[] = {100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,
                     1000.0};
    double norm[] = {1000.0,100.0,200.0,500.0,550.0,600.0,650.0,700.0,
                     750.0,800.0};
    casu_fits *inf,*inconf;
    cpl_propertylist *pl;

    /* Generate a field with some stars and a confidence map */

    bkg = cpl_image_new(1024,1024,CPL_TYPE_FLOAT);
    im = cpl_image_new(1024,1024,CPL_TYPE_FLOAT);
    cnf = cpl_image_new(1024,1024,CPL_TYPE_FLOAT);
    norm2 = 2.0*CPL_MATH_PI*sigma*sigma;
    cpl_image_fill_noise_uniform(bkg,-10.0,10.0);
    cpl_image_add_scalar(bkg,sky);
    cpl_image_fill_noise_uniform(cnf,99.9,100.1);
    for (i = 0; i < NTEST; i++) {
        cpl_image_fill_gaussian(im,xpos[i],ypos[i],norm[i]*norm2,sigma,sigma);
        tot[i] = cpl_image_get_flux(im);
        cpl_image_add(bkg,im);
    }
    pl = cpl_propertylist_new();
    inf = casu_fits_wrap(bkg,NULL,NULL,pl);
    inconf = casu_fits_wrap(cnf,NULL,NULL,NULL);
    cpl_image_delete(im);
    cpl_propertylist_delete(pl);

    /* Give it a WCS */

    fill_plist(casu_fits_get_ehu(inf));

    /* Run imcore */

    hdrl_imcore_result * res = cpl_calloc(sizeof(*res), 1);

    retval = imcore_conf(inf,inconf,5,1.5,0,5.0,1,64,6,3.0,1.0,
                    HDRL_SATURATION_INIT,res);
    cpl_test_eq(retval,CASU_OK);
    cpl_test_nonnull(res->catalogue);

    cpl_image_delete(res->segmentation_map); // TODO check results
    cpl_image_delete(res->background);

    /* Check the results. Start by checking the number of rows and columns. 
       Sort the table by X */

    tab = casu_tfits_get_table(res->catalogue);
    cpl_test_nonnull(tab);
    ncol = cpl_table_get_ncol(tab);
    cpl_test_eq(ncol, NCOLS);
    nrow = cpl_table_get_nrow(tab);
    cpl_test_eq(nrow,NTEST);
    pl = cpl_propertylist_new();
    cpl_propertylist_append_bool(pl,"X_coordinate",0);
    cpl_table_sort(tab,pl);
    cpl_propertylist_delete(pl);

    /* Test the column content of the table */

    for (i = 0; i < NTEST; i++) {
        cpl_test_abs(xpos[i],cpl_table_get_float(tab,"X_coordinate",(cpl_size)i,
                                                 &nl),0.2);
        cpl_test_abs(ypos[i],cpl_table_get_float(tab,"Y_coordinate",(cpl_size)i,
                                                 &nl),0.2);
        diff = fabs(cpl_table_get_float(tab,"Aper_flux_5",(cpl_size)i,&nl) - 
                    tot[i]);
        diff /= cpl_table_get_float(tab,"Aper_flux_5_err",(cpl_size)i,&nl);
        cpl_test_lt(diff,1.6);
    }

    /* Get out of here */

    casu_fits_delete(inf);
    casu_fits_delete(inconf);
    casu_tfits_delete(res->catalogue);
    cpl_free(res);
    return cpl_error_get_code();

}
/*----------------------------------------------------------------------------*/
/**
  @brief Check background subtraction of imcore_conf() in various conditions
  @return cpl_error_code
 */
/*----------------------------------------------------------------------------*/
static cpl_error_code hdrl_catalogue_imcore_conf_backsub(void)
{

    int retval, nrow, nl;
    cpl_image *im, *cnf;
    cpl_table *tab;
    casu_fits *inf, *inconf;
    cpl_propertylist *pl;

    im = cpl_image_new(1024,1024,CPL_TYPE_FLOAT);
    cnf = cpl_image_new(1024,1024,CPL_TYPE_FLOAT);
    pl = cpl_propertylist_new();

    cpl_image_fill_noise_uniform(im,-10.0,10.0);
    //cpl_image_fill_noise_uniform(cnf,99.9,100.1);
    cpl_image_add_scalar(cnf,100.);

    /*Set a region arround the object to zero*/
    for (int x = 400; x < 600; x++) {
        for (int y = 400; y < 600; y++) {
            cpl_image_set(im, x, y, 0.);
        }
    }

    /* Insert an tophat object */
    for (int x = 500; x < 505; x++) {
        for (int y = 500; y < 505; y++) {
            cpl_image_set(im, x, y, 4990.);
        }
    }
    cpl_image_add_scalar(im, 10.);
    /*Now the object has a value of 5000 == 4990 +10 */

    //cpl_image_save(im, "armin.fits",CPL_TYPE_FLOAT,NULL,CPL_IO_CREATE);
    inf = casu_fits_wrap(im,NULL,NULL,pl);
    inconf = casu_fits_wrap(cnf,NULL,NULL,NULL);

    /* Run imcore */

    hdrl_imcore_result * res = cpl_calloc(sizeof(*res), 1);
    retval = imcore_conf(inf,inconf,5,3.0,0,1.0,0,16,7,3.0,1.0,
                    HDRL_SATURATION_INIT,res);
    cpl_test_eq(retval,CASU_OK);
    tab = casu_tfits_get_table(res->catalogue);
    cpl_test_nonnull(tab);

/*
    cpl_table_save(tab, NULL, NULL, "catalogue_nobacksub.fits", CPL_IO_CREATE);
    cpl_image_save(res->background, "background.fits", CPL_TYPE_FLOAT, NULL,
                   CPL_IO_CREATE);
*/

    /* Check the results. */

    nrow = cpl_table_get_nrow(tab);
    cpl_test_eq(nrow,1);

    /* Test the column content of the table */
    cpl_test_abs(5000., cpl_table_get_float(tab, "Aper_flux_1", 0, &nl),0.2);



    casu_fits_delete(inf);
    casu_fits_delete(inconf);
    casu_tfits_delete(res->catalogue);
    cpl_image_delete(res->segmentation_map);
    cpl_image_delete(res->background);

    cpl_free(res);
    cpl_propertylist_delete(pl);
    return cpl_error_get_code();
}

static cpl_error_code test_gaussians(double dx, double dy, double fwhm)
{
    cpl_propertylist *pl;

    /* Generate a field with some stars and a confidence map */
    cpl_image * im = create_gauss(dx, dy, fwhm);
    pl = cpl_propertylist_new();
    casu_fits * inf = casu_fits_wrap(im, NULL, NULL, pl);
    cpl_propertylist_delete(pl);

    /* Give it a WCS */

    fill_plist(casu_fits_get_ehu(inf));

    /* Run imcore */

    hdrl_imcore_result * res = cpl_calloc(sizeof(*res), 1);

    int retval = imcore_conf(inf,NULL,5,2.5,0,fwhm,1,fwhm*3,6,3.0,1.0,
                             HDRL_SATURATION_INIT,res);
    cpl_test_eq(retval,CASU_OK);
    cpl_test_nonnull(res->catalogue);

    cpl_image_delete(res->segmentation_map); // TODO check results
    cpl_image_delete(res->background);

    /* Check the results. Start by checking the number of rows and columns.
       Sort the table by X */

    cpl_table * tab = casu_tfits_get_table(res->catalogue);
    cpl_test_nonnull(tab);
    cpl_size ncol = cpl_table_get_ncol(tab);
    cpl_test_eq(ncol, NCOLS);
    cpl_size nrow = cpl_table_get_nrow(tab);
    cpl_test_eq(nrow, 1);

    int nl;
    cpl_test_abs(fwhm, cpl_table_get_float(tab, "FWHM", 0, &nl), 0.006);

    /* Get out of here */

    casu_fits_delete(inf);
    casu_tfits_delete(res->catalogue);
    cpl_free(res);
    return cpl_error_get_code();

}


/*----------------------------------------------------------------------------*/
/**
  @brief   Unit tests of imcore_conf in the catalogue module
 **/
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);
    hdrl_catalogue_imcore_conf_basic();
    hdrl_catalogue_imcore_conf_backsub();
    test_gaussians(.5, .5, 3.);
    test_gaussians(.0, .5, 4.);
    test_gaussians(.5, .5, 5.);
    test_gaussians(.8, .5, 6.);
    test_gaussians(.1, .2, 7.);
    cpl_test_error(CPL_ERROR_NONE);

    return cpl_test_end(0);
}
