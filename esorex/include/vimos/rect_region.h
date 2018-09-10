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
 * rect_region.h
 *
 *  Created on: 2013 11 25
 *      Author: cgarcia
 */

#ifndef RECT_REGION_H
#define RECT_REGION_H

#include <vector>
#include "hdrl.h"

namespace mosca 
{

//A  rectangular region.
class rect_region
{
public:
    rect_region(int llx_pix, int lly_pix, int urx_pix, int ury_pix);
    
    rect_region();
    
    int llx() const;
    
    int lly() const;
    
    int urx() const;
    
    int ury() const;
    
    void set_llx(int llx_pix);
    
    void set_lly(int lly_pix);
    
    void set_urx(int urx_pix);
    
    void set_ury(int ury_pix);
    
    int length_x() const;
    
    int length_y() const;

    int pix_inc_x() const;
    
    int pix_inc_y() const;

    int area_geom() const;

    int area_pix_inc() const;
    
    bool is_empty() const;

    rect_region coord_0to1(void) const;

    rect_region coord_1to0(void) const;
    
    hdrl_parameter * hdrl_param();
    
    virtual ~rect_region();
    
private:
    
    int m_llx;
    
    int m_lly;
    
    int m_urx;
    
    int m_ury;
    
    hdrl_parameter * m_hdrl_param;

    bool m_is_empty;
    
    friend bool operator ==
      (const rect_region& lhs, const rect_region& rhs);
    
};

rect_region rect_region_minenclose
(const rect_region & reg1, const rect_region & reg2);

rect_region rect_region_minenclose
(const rect_region & reg1, const rect_region & reg2, const rect_region & reg3);

rect_region rect_region_minenclose
(const std::vector<rect_region> & regions);

bool operator == (const mosca::rect_region& lhs, const mosca::rect_region& rhs);

bool operator != (const mosca::rect_region& lhs, const mosca::rect_region& rhs);

} /* namespace mosca */

#endif /* RECT_REGION_H_ */
