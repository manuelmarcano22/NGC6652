/* $Id: vimos_calib_impl.h,v 1.1 2010-06-16 11:07:49 cizzo Exp $
 *
 * This file is part of the FORS Library
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
 * $Author: cizzo $
 * $Date: 2010-06-16 11:07:49 $
 * $Revision: 1.1 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VIMOS_CALIB_IMPL_H
#define VIMOS_CALIB_IMPL_H

#include <cpl.h>

CPL_BEGIN_DECLS

int vimos_calib_impl(cpl_frameset *frames, cpl_parameterlist *parameters);

CPL_END_DECLS

#endif
