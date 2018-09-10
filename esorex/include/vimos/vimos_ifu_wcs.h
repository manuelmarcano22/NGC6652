/* $Id: VimosIfuWCS.h,v 1.3 2013-03-25 11:43:04 cgarcia Exp $
 *
 * This file is part of the VIMOS Pipeline
 * Copyright (C) 2002-2004 European Southern Observatory
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
 * $Revision: 1.3 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VIMOS_IFU_WCS_H
#define VIMOS_IFU_WCS_H

#include "two_d_linear_wcs.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup VimosIfuWCS       IFU WCS routines
 */
/*----------------------------------------------------------------------------*/

/** 
 *  
 * Gets the WCS IFU keywords based on the 
 * pointing of the instrument. The correction is needed because of 
 * IFU offset wrt the instrument optical axis.   
 *  Algorithm provided by P. Hammersley
 *  Originally programmed in Fortran by M. Rejkuba
 *  Testing and OB data: M. Rejkuba
 * @param raVimos     RA coordinate where VIMOS is pointing
 * @param decVimos    DEC coordinate where VIMOS is pointing
 * @param posAngle    Position angle of the rotator
 * @param dimX        The number of fibers in X in the IFU FOV 
 * @param dimY        The number of fibers in Y in the IFU FOV
 * @param fiberSize   The space between two fibers
 * @param epoch       Epoch of the WCS solution
 * @param equinox     Equinox to which the coordinates are refered.
 * @return  The WCS solution
 */
two_d_linear_wcs vimos_ifu_get_2d_wcs_from_pointing
(double raVimos,  double decVimos, 
 double posAngle,
 cpl_size nxFibers, cpl_size nyFibers, double fiberSize, 
 double epoch, double equinox);


#endif /* VIMOS_IFU_WCS_H */
