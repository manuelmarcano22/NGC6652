/* $Id: VimosUtils.tcc,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2012 European Southern Observatory
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
 * $Date: 2013-03-25 11:43:04 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VIMOSUTILS_TCC
#define VIMOSUTILS_TCC


#include <cmath>
#include "vimos_utils.h"
#include "cpl_propertylist.h"


//Forward declaration.
//This is used to get the proper function to read from the header a 
//keyword of type T (using memeber method get_key_func()).  
template <typename T>
struct header_trait;

template<typename T> 
bool vimos_check_equal_keys
(PilSetOfFrames * frames, 
 const std::string& tag, 
 const std::string&  keyname,
 T& keyvalue)
{
    bool first = true;
    T     value = T();
    
    if(pilSofFrameCount(frames, tag.c_str()) < 2)
       return true;

    PilFrame *frame = pilSofFirst((PilSetOfFrames *)frames);
    while (frame != NULL) {
        cpl_propertylist * header;
        header = cpl_propertylist_load(pilFrmGetName(frame),0);
        if(first == true)
        {
            value = header_trait<T>::get_key_func(header, keyname);
            first = false;
        }
        else
            if(value != header_trait<T>::get_key_func(header, keyname))
            {
                cpl_propertylist_delete(header);
                return false;
            }
        cpl_propertylist_delete(header);
        frame = pilSofNext((PilSetOfFrames *)frames, frame);
    }
    keyvalue = value;
    return true;
}

template <typename T>
struct FitsTypes_trait
{
};

template <>
struct header_trait<double>
{
    static double get_key_func
    (cpl_propertylist * header, const std::string& keyname)
    {
        return cpl_propertylist_get_double(header, keyname.c_str());
    }
};

template <>
struct header_trait<int>
{
    static int get_key_func
    (cpl_propertylist * header, const std::string& keyname)
    {
        return cpl_propertylist_get_int(header, keyname.c_str());
    }
};


#endif
