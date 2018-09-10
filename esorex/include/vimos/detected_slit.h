/* $Id: detected_slit.h,v 1.4 2013-10-24 08:52:49 cgarcia Exp $
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
 * $Date: 2013-10-24 08:52:49 $
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifndef DETECTED_SLIT_H
#define DETECTED_SLIT_H

#include <vector>
#include "cpl_polynomial.h"

namespace mosca
{

/**
 * This class represents a detected slit in a spectroscopic image.
 * TODO: Separate this in two: the pure detected slit (in terms of
 * physical extent) and a traced slit.
 */
class detected_slit
{
public:
    /* 
     * Constructor from the extents
     * @param slit_id          The identification of the slit
     * @param disp_bottom      The coordinate of the bottom of the slit in the 
     *                         spectral axis
     * @param spa_bottom       The coordinate of the bottom of the slit in the 
     *                         spatial axis
     * @param disp_top         The coordinate of the top of the slit in the 
     *                         spectral axis
     * @param spa_top          The coordinate of the top of the slit in the 
     *                         spatial axis
     * @param trace_bottom_pol_coeffs The coefficients of the polynomail fitted to the
     *                                bottom of the slit trace
     * @param trace_top_pol_coeffs    The coefficients of the polynomail fitted to the
     *                                top of the slit trace
     * 
     */
    detected_slit(int slit_id, double disp_bottom, double spa_bottom, 
                  double disp_top, double spa_top,
                  int position_spatial_corrected,
                  int length_spatial_corrected,
                  const std::vector<double>& trace_bottom_pol_coeffs,                  
                  const std::vector<double>& trace_top_pol_coeffs);
    
    detected_slit();

    detected_slit(const detected_slit& rhs);

    virtual ~detected_slit();
    
    void get_extent(double& disp_bottom, double& spa_bottom, 
                    double& disp_top,    double& spa_top);
    
    void get_extent_pix(int& disp_bottom, int& spa_bottom, 
                        int& disp_top,    int& spa_top) const;
    
    int get_position_spatial_corrected() const;

    int get_length_spatial_corrected() const;
    
    bool within_trace(double dispersion_pos, double spatial_pos) const;
    
    double spatial_correct(double dispersion_pos,
                           double spatial_pos) const;
    
    int slit_id() const;

private:
    
    //The identification of the slit
    int m_slit_id;
    
    //The coordinate of the bottom of the slit in the spectral axis
    double m_disp_bottom;
    
    //The coordinate of the bottom of the slit in the spatial axis
    double m_spa_bottom;
    
    //The coordinate of the top of the slit in the spectral axis
    double m_disp_top;
    
    //The coordinate of the top of the slit in the spatial axis
    double m_spa_top;
  
    //Position of the bottom of the slit in an image already corrected from spatial distortion
    int m_position_spatial_corrected;
    
    //Length of the slit in an image already corrected from spatial distortion
    int m_length_spatial_corrected;
    
    //The coefficients of the polynomial fit to the bottom of the trace 
    std::vector<double> m_trace_bottom_pol_coeffs;
    
    //The coefficients of the polynomial fit to the top of the trace 
    std::vector<double> m_trace_top_pol_coeffs;

    //The polynomial with the fit to the bottom of the trace
    cpl_polynomial * m_trace_bottom_pol;

    //The polynomial with the fit to the top of the trace
    cpl_polynomial * m_trace_top_pol;
};

}

#endif
