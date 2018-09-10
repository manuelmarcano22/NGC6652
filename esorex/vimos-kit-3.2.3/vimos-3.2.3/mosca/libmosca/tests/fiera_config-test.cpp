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
#define BOOST_TEST_MODULE fiera_config-test
#define _GLIBCXX_USE_NANOSLEEP
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <cpl.h>
#include "fiera_config.h"

BOOST_AUTO_TEST_SUITE(fiera_config_fors_headers)

BOOST_AUTO_TEST_CASE(fors_marlene_out1_bin2)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/marlene_out1_bin2.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 4);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 1029);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 1033);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 5);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 1028);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(fors_tk2048eb4_1_160_out1_bin1)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/tk2048eb4_1_160_out1_bin1.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 15);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 2047);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 2064);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2079);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 2047);

    //Check valid pixels region port 1
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).llx(), 16);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).urx(), 2063);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).ury(), 2047);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(fors_tk2048eb4_1_160_out4_bin1)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/tk2048eb4_1_160_out4_bin1.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 15);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 1023);

    //Check overscan region port 1
    BOOST_CHECK(ccd_config.overscan_region(0).is_empty());

    //Check valid pixels region port 1
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).llx(), 16);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).urx(), 1039);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).ury(), 1023);

    //Check prescan region port 2
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(1).llx(), 2064);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(1).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(1).urx(), 2079);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(1).ury(), 1023);

    //Check overscan region port 2
    BOOST_CHECK(ccd_config.overscan_region(1).is_empty());

    //Check valid pixels region port 2
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(1).llx(), 1040);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(1).lly(), 0);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(1).urx(), 2063);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(1).ury(), 1023);

    //Check prescan region port 3
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(2).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(2).lly(), 1024);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(2).urx(), 15);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(2).ury(), 2047);

    //Check overscan region port 3
    BOOST_CHECK(ccd_config.overscan_region(2).is_empty());

    //Check valid pixels region port 3
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(2).llx(), 16);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(2).lly(), 1024);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(2).urx(), 1039);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(2).ury(), 2047);

    //Check prescan region port 4
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(3).llx(), 2064);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(3).lly(), 1024);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(3).urx(), 2079);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(3).ury(), 2047);

    //Check overscan region port 4
    BOOST_CHECK(ccd_config.overscan_region(3).is_empty());

    //Check valid pixels region port 4
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(3).llx(), 1040);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(3).lly(), 1024);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(3).urx(), 2063);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(3).ury(), 2047);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(fors_tk2048eb4_1_160_out4_bin2)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/tk2048eb4_1_160_out4_bin2.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 7);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 511);

    //Check overscan region port 1
    BOOST_CHECK(ccd_config.overscan_region(0).is_empty());

    //Check valid pixels region port 1
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).llx(), 8);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).urx(), 511);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(0).ury(), 511);

    //Check prescan region port 2
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(1).llx(), 1032);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(1).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(1).urx(), 1039);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(1).ury(), 511);

    //Check overscan region port 2
    BOOST_CHECK(ccd_config.overscan_region(1).is_empty());

    //Check valid pixels region port 2
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(1).llx(), 528);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(1).lly(), 0);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(1).urx(), 1031);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(1).ury(), 511);

    //Check prescan region port 3
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(2).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(2).lly(), 512);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(2).urx(), 7);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(2).ury(), 1023);

    //Check overscan region port 3
    BOOST_CHECK(ccd_config.overscan_region(2).is_empty());

    //Check valid pixels region port 3
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(2).llx(), 8);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(2).lly(), 512);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(2).urx(), 511);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(2).ury(), 1023);

    //Check prescan region port 4
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(3).llx(), 1032);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(3).lly(), 512);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(3).urx(), 1039);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(3).ury(), 1023);

    //Check overscan region port 4
    BOOST_CHECK(ccd_config.overscan_region(3).is_empty());

    //Check valid pixels region port 4
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(3).llx(), 528);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(3).lly(), 512);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(3).urx(), 1031);
    BOOST_REQUIRE_EQUAL(ccd_config.validpix_region(3).ury(), 1023);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(fors_ccid20_14_5_6_out1_bin2)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/ccid20_14_5_6_out1_bin2.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 4);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 1029);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 1033);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 5);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 1028);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(fors_ccid20_14_5_3_out1_bin2)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/ccid20_14_5_3_out1_bin2.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 1029);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 1033);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 4);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 5);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2047);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 1028);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(fors_ccid20_14_5_3_out1_bin1)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/ccid20_14_5_3_out1_bin1.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 2058);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 4095);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 2067);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 4095);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 9);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 10);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 4095);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 2057);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(vimos_eev_ccd_44__ccd_59a)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/eev_ccd_44__ccd_59a.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 49);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 4095);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 2098);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2147);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 4095);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 50);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2097);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 4095);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(vimos_eev_ccd_44__ccd_59b)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/eev_ccd_44__ccd_59b.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 49);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 4095);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 2098);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2147);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 4095);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 50);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2097);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 4095);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(vimos_eev_ccd_44__ccd_60a)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/eev_ccd_44__ccd_60a.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 49);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 4095);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 2098);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2147);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 4095);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 50);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2097);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 4095);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(vimos_eev_ccd_44__ccd_60b)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/eev_ccd_44__ccd_60b.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 49);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 4095);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 2098);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2147);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 4095);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 50);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2097);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 4095);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(vimos_eev_ccd_44_82__david)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
           cpl_propertylist_load(DATADIR "data/eev_ccd_44_82__david.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 49);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 2439);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 2098);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2147);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 2439);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 50);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2097);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 2439);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(vimos_eev_ccd_44_82__tom)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
            cpl_propertylist_load(DATADIR "data/eev_ccd_44_82__tom.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 49);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 2439);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 2098);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2147);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 2439);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 50);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2097);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 2439);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_CASE(vimos_eev_ccd_44_82__tom_ny_4096)
{
    cpl_init(CPL_INIT_DEFAULT);
    cpl_propertylist * header = 
     cpl_propertylist_load(DATADIR "data/eev_ccd_44_82__tom_ny_4096.fits.gz",0);
    mosca::fiera_config ccd_config(header);

    //Check prescan region port 1
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).llx(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).urx(), 49);
    BOOST_CHECK_EQUAL(ccd_config.prescan_region(0).ury(), 4095);

    //Check overscan region port 1
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).llx(), 2098);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).urx(), 2147);
    BOOST_CHECK_EQUAL(ccd_config.overscan_region(0).ury(), 4095);

    //Check valid pixels region port 1
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).llx(), 50);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).lly(), 0);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).urx(), 2097);
    BOOST_CHECK_EQUAL(ccd_config.validpix_region(0).ury(), 4095);

    cpl_propertylist_delete(header);
    cpl_end();
}

BOOST_AUTO_TEST_SUITE_END()

#else

int main(void)
{
    return 0;
}

#endif
