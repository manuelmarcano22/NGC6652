/* $Id: moses.h,v 1.41 2013/09/09 12:19:20 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
 * $Author: cgarcia $
 * $Date: 2013/09/09 12:19:20 $
 * $Revision: 1.41 $
 * $Name:  $
 */

#ifndef VIMOS_OVERSCAN_H
#define VIMOS_OVERSCAN_H

#include <cpl.h>

#ifdef __cplusplus

#include "fiera_config.h"
#include "mosca_image.h"

class vimos_preoverscan
{
public:

    mosca::image subtract_prescan(mosca::image& image, 
                                  const mosca::ccd_config& ccd_config);

    std::vector<mosca::image> subtract_prescan(std::vector<mosca::image>& ima_list, 
            const mosca::ccd_config& ccd_config);

    mosca::image subtract_overscan(mosca::image& image, 
                                  const mosca::ccd_config& ccd_config);

    std::vector<mosca::image> subtract_overscan(std::vector<mosca::image>& ima_list, 
            const mosca::ccd_config& ccd_config);

    mosca::image trimm_preoverscan(mosca::image& image, 
                                   const mosca::ccd_config& ccd_config);

    std::vector<mosca::image> trimm_preoverscan
         (std::vector<mosca::image>& ima_list, 
          const mosca::ccd_config& ccd_config);

    void fix_wcs_trimm(cpl_propertylist * header);

    double get_median_correction() const;
private:

    double m_median_correction;
};

#endif

//TODO: Remove this if all the recipes can call C++ code.
CPL_BEGIN_DECLS

cpl_image * vimos_subtract_prescan(cpl_image * image, 
                                   cpl_image * image_var, 
                                   cpl_propertylist * header);

cpl_image * vimos_subtract_overscan(cpl_image * image, 
                                   cpl_image * image_var, 
                                   cpl_propertylist * header);

cpl_image * vimos_trimm_preoverscan(cpl_image * image, 
                                    cpl_propertylist * header);

CPL_END_DECLS

//void vimos_trimm_fill_info(cpl_propertylist * header,
                          //const mosca::ccd_config& ccd_config);

#endif   /* VIMOS_OVERSCAN_H */
