/* $Id: reduce_method.h,v 1.2 2013-10-16 11:30:33 cgarcia Exp $
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
 * $Date: 2013-10-16 11:30:33 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef REDUCE_METHODS_H
#define REDUCE_METHODS_H

#include "hdrl.h"

namespace mosca
{

    /**
     * This reduce method will simply compute the mean over a list of values  
     */
    class reduce_mean
    {
    public:
        reduce_mean();
        
        /**
         * Gets the HDRL reduce function
         * @return
         */
        hdrl_parameter  * hdrl_reduce();
    };

    /**
     * This reduce method will simply compute the median over a list of values  
     */
    class reduce_median
    {
    public:
        reduce_median();
        
        /**
         * Gets the HDRL reduce function
         * @return
         */
        hdrl_parameter  * hdrl_reduce();
    };

    /**
     * This reduce method will compute a kappa-sigma clipping algorithm
     * over a list of values.  
     */
    class reduce_sigma_clipping
    {
    public:
        reduce_sigma_clipping(double kappa_low, double kappa_high, int iter);
        
        /**
         * Gets the HDRL reduce function
         * @return
         */
        hdrl_parameter  * hdrl_reduce();
        
    private:
        double m_kappa_high;
        double m_kappa_low;
        int m_iter;
    };

    /**
     * This reduce method will compute a weighted mean algorithm
     * over a list of values.  
     */
    class reduce_weighted_mean
    {
    public:
        reduce_weighted_mean();
        
        /**
         * Gets the HDRL reduce function
         * @return
         */
        hdrl_parameter  * hdrl_reduce();
    };
}

#endif
