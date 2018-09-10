/* $Id: reduce_method.cpp,v 1.2 2013-10-16 11:30:33 cgarcia Exp $
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


#include "reduce_method.h"

mosca::reduce_mean::reduce_mean()
{
}

hdrl_parameter * mosca::reduce_mean::hdrl_reduce()
{
    return hdrl_collapse_mean_parameter_create();
}

mosca::reduce_median::reduce_median()
{
}

hdrl_parameter  * mosca::reduce_median::hdrl_reduce()
{
    return hdrl_collapse_median_parameter_create();
}

mosca::reduce_weighted_mean::reduce_weighted_mean()
{
}

hdrl_parameter  * mosca::reduce_weighted_mean::hdrl_reduce()
{
    return hdrl_collapse_weighted_mean_parameter_create();
}

mosca::reduce_sigma_clipping::reduce_sigma_clipping
(double kappa_low, double kappa_high, int iter) :
m_kappa_high(kappa_high), m_kappa_low(kappa_low), m_iter(iter)
{
}

hdrl_parameter  * mosca::reduce_sigma_clipping::hdrl_reduce()
{
    return hdrl_collapse_sigclip_parameter_create
               (m_kappa_low, m_kappa_high, m_iter);
}
