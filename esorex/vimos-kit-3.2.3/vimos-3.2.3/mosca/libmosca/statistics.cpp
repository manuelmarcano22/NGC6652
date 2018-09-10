/* $Id: mosca_image.tcc,v 1.6 2013-07-24 07:40:54 cgarcia Exp $
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
 * $Date: 2013-07-24 07:40:54 $
 * $Revision: 1.6 $
 * $Name: not supported by cvs2svn $
 */

#ifndef STATISTICS_CPP
#define STATISTICS_CPP

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include "statistics.h"


template <typename ITER>
double mosca::mean(const ITER& begin, const ITER& end)
{
    long double mean = 0;
    ITER iter;
    size_t count = 0;

    for(iter = begin; iter != end; ++iter, ++count)
        mean += (*iter -mean) / (count+1);

    return mean;
}

template <typename ITER>
double mosca::variance(const ITER& begin, const ITER& end)
{
    long double variance = 0 ;
    ITER iter;
    size_t count = 0;

    double mean = mosca::mean(begin, end);
    
    /* find the sum of the squares */
    for(iter = begin; iter != end; ++iter, ++count)
    {
        const long double delta = (*iter - mean);
        variance += (delta * delta - variance) / (count + 1);
    }

    return variance;
}

template <typename ITER>
double mosca::robust_variance(const ITER& begin, const ITER& end)
{
    double first_quartile, median, third_quartile;
    mosca::quartile(begin, end, first_quartile, median, third_quartile);
    //https://en.wikipedia.org/wiki/Interquartile_range#Interquartile_range_of_distributions
    double sigma = (third_quartile - first_quartile) / 1.349;
    double variance = sigma * sigma;
    return variance;
}

template <typename ITER>
void mosca::quartile(const ITER& begin, const ITER& end,
                     double& first_quartile, double& median, 
                     double& third_quartile)
{
    ITER iter;
    int k;
    int size = 0;
    
    for (iter = begin, k = 0; iter != end ;++iter, ++k)
      size++;

    double * data = new double[size];

    for (iter = begin, k = 0; iter != end ;++iter, ++k)
        data[k] = double(*iter);

    gsl_sort(data, 1, size);
    median = gsl_stats_median_from_sorted_data(data, 1, size);
    // sigma = (third - first) / 1.35
    first_quartile = gsl_stats_quantile_from_sorted_data(data, 1, size, 0.25);
    third_quartile = gsl_stats_quantile_from_sorted_data(data, 1, size, 0.75);
    delete[] data;
}




#endif
