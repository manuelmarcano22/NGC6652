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
 * fiera_ccd.h
 *
 *  Created on: 2013 11 25
 *      Author: cgarcia
 */

#ifndef MOSCA_FIERA_CONFIG_H
#define MOSCA_FIERA_CONFIG_H

#include <string>
#include "ccd_config.h"
#include "cpl_propertylist.h" 

namespace mosca 
{

class fiera_config : public ccd_config
{
public:
    
    fiera_config(const cpl_propertylist * header);
    
    fiera_config();

    virtual ~fiera_config();
    
protected:
    
    std::string m_chip_id; 
};

} /* namespace mosca */
#endif /* MOSCA_FIERA_CONFIG_H */
