/*
 * This file is part of the MOSCA library
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

#ifndef STATISTICS_H
#define STATISTICS_H

namespace mosca 
{

template <typename ITER>
double mean(const ITER& begin, const ITER& end);

template <typename ITER>
double variance(const ITER& begin, const ITER& end);

template <typename ITER>
double robust_variance(const ITER& begin, const ITER& end);

template <typename ITER>
void quartile(const ITER& begin, const ITER& end,
              double& first_quartile, double& median, 
              double& third_quartile);

} /* namespace mosca */

#endif /* STATISTICS_H */

#include "statistics.cpp"
