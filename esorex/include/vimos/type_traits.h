/* $Id: type_traits.h,v 1.4 2013-03-25 14:23:56 cgarcia Exp $
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
 * $Date: 2013-03-25 14:23:56 $
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifndef TYPE_TRAITS_H
#define TYPE_TRAITS_H

#include "cpl_type.h"

namespace mosca
{
template<typename T>
struct type_trait
{    
};

template<>
struct type_trait<double>
{    
    static const cpl_type cpl_eq_type = CPL_TYPE_DOUBLE;
};

template<>
struct type_trait<float>
{    
    static const cpl_type cpl_eq_type = CPL_TYPE_FLOAT;
};

template<>
struct type_trait<int>
{    
    static const cpl_type cpl_eq_type = CPL_TYPE_INT;
};

template<>
struct type_trait<long>
{    
    static const cpl_type cpl_eq_type = CPL_TYPE_LONG;
};

template<>
struct type_trait<bool>
{    
    static const cpl_type cpl_eq_type = CPL_TYPE_BOOL;
};

}

#endif
