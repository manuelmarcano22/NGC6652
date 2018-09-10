/* $Id: vimos_calib_impl.c,v 1.22 2013-10-22 16:57:09 cgarcia Exp $
 *
 * This file is part of the VIMOS Data Reduction Pipeline
 * Copyright (C) 2006 European Southern Observatory
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

#include <string>
#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>
#include "vimos_lines.h"
#include <vimos_dfs.h>

//Get the list of lines after proper filtering
cpl_vector * vimos_lines_filter(cpl_frameset * frameset, 
                                const char * ignore_lines,
                                const char * used_linesets)
{
    cpl_table        *wavelengths  = NULL;
    cpl_size          nlines_all;
    cpl_size          n_selected = 0;
    cpl_size          i;
    cpl_vector       *lines;
    double            lambda;
    int               null;
    const char *      wcolumn = "WLEN";
    const char *      linesetcolumn = "LINE_SET";

    /*
     * Read the wavelengths table 
     */
    wavelengths = dfs_load_table(frameset, "LINE_CATALOG", 1);

    if (wavelengths == NULL)
    {
        cpl_msg_error(cpl_func, "Cannot load line catalog");
        return NULL;
    }


    nlines_all = cpl_table_get_nrow(wavelengths);

    if (nlines_all == 0)
    {
        cpl_msg_error(cpl_func, "Empty input line catalog");
        cpl_table_delete(wavelengths);
        return NULL;
    }

    if (cpl_table_has_column(wavelengths, wcolumn) != 1 || 
        cpl_table_has_column(wavelengths, linesetcolumn) != 1)
    {
        cpl_msg_error(cpl_func, "One or several of the following columns "
                      "are not present in input line catalog: %s %s",
                      wcolumn, linesetcolumn);
        cpl_table_delete(wavelengths);
        return NULL;
    }

    /*
     * Deselect lines which are not in the lineset 
     */
    //Parse command line option used_linesets
    std::stringstream used_linesets_ss(used_linesets);
    std::vector<std::string> linesets;
    std::string lineset;
    //Parsing used_linesets (values are separated by comma)
    while(std::getline(used_linesets_ss, lineset, ','))
        linesets.push_back(lineset);

    //Deselect lines which are not in the lineset
    for(cpl_size iline = 0; iline < nlines_all; iline++)
    {
        if(cpl_table_get_string(wavelengths, linesetcolumn, iline) == NULL)
        {
            cpl_table_unselect_row(wavelengths, iline);
            break;
        }
        std::string table_lineset = 
                cpl_table_get_string(wavelengths, linesetcolumn, iline);
        if(std::find(linesets.begin(), linesets.end(), table_lineset) == linesets.end())
            cpl_table_unselect_row(wavelengths, iline);
    }

    /*
     * Deselect lines which are present in ignore_lines 
     */
    std::string ignore_lines_str(ignore_lines);
    while(ignore_lines_str.length() > 0)
    {
        //Parsing ignore_lines (values are separated by comma)
        std::string::size_type found = ignore_lines_str.find(',');
        std::string lambda_str;
        if(found != std::string::npos)
        {
            lambda_str = ignore_lines_str.substr(0, found);
            ignore_lines_str = ignore_lines_str.substr(found+1);
        }
        else
        {
            lambda_str = ignore_lines_str;
            ignore_lines_str = "";
        }
        std::istringstream iss(lambda_str);
        if ( !(iss >> lambda) || !(iss.eof() || (iss >> std::ws && iss.eof())) )
        {
            cpl_msg_error(cpl_func, "Cannot interpret number in ignored_lines");
            cpl_table_delete(wavelengths);
            return NULL;
        }

        //Search for closest line in catalog. The line is unselected but
        //it will be checked again against the next ignored line. In this way,
        //if a value appears many times in the ignored_lines, only one line
        //will be removed
        cpl_size i_ignore = 0;
        double min_lambda_dif = 
             std::fabs(lambda - cpl_table_get(wavelengths, wcolumn, 0, &null));
        for (i = 1; i < nlines_all; i++)
        {
            double lambda_dif = 
              std::fabs(lambda - cpl_table_get(wavelengths, wcolumn, i, &null));
            if(lambda_dif < min_lambda_dif)
            {
                min_lambda_dif = lambda_dif;
                i_ignore = i;
            }
         }
        cpl_table_unselect_row(wavelengths, i_ignore);
    } 
    
    /* Create the final list of reference lines */
    n_selected = cpl_table_count_selected(wavelengths);
    lines = cpl_vector_new(n_selected);
    cpl_size i_line = 0;
    for (i = 0; i < nlines_all; i++)
    {
        lambda = cpl_table_get(wavelengths, wcolumn, i, &null);
        if(cpl_table_is_selected(wavelengths, i))
        {
            cpl_vector_set(lines, i_line, lambda);
            i_line++;
        }
    }

    cpl_table_delete(wavelengths);

    return lines;
}

/**@}*/
 
