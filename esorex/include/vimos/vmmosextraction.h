/* $Id: vmmosextraction.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VM_MOSEXTRACTION_H
#define VM_MOSEXTRACTION_H

#include <pilmacros.h>

#include <vmextractiontable.h>
#include <vmwindowtable.h>
#include <vmobjecttable.h>


PIL_BEGIN_DECLS

enum _VIMOS_SPEC_EDGE_METHOD_ 
{
  VM_SP_EDGE_UNDEF = 0,
  VM_SP_EDGE_GRAD  = 1
};

typedef enum _VIMOS_SPEC_EDGE_METHOD_ VimosSpecEdgeMethod;


enum _VIMOS_SPEC_SAMPLE_OPTION_
{
  VM_SP_LIN_LAMBDA = 0,
  VM_SP_LOG_LAMBDA = 1
};

typedef enum _VIMOS_SPEC_SAMPLE_OPTION_ VimosSpecSampleOption;


VimosWindowTable *VmSpDetObj(VimosImage *, VimosExtractionTable *, int, int,
                             float, float, float, int, int, float, float);
VimosImage **VmSpEx2D(VimosImage **, VimosExtractionTable *,
		      VimosSpecSampleOption);
VimosImage **VmSpEx1D(VimosImage **, VimosWindowTable *,
		      VimosObjectTable *, int, int);
VimosImage *VmSpStack2D(VimosImage **, VimosWindowTable **,
                        VimosExtractionTable *, VimosIntArray *, int,
                        CombMethod, CombParameters *, VimosFloatArray *,
                        int, int, int);
VimosExtractionTable *VmSpExTab(VimosImage *, VimosTable *, VimosIfuTable *,
				VimosExtractionTable *);
int VmSpMatch(VimosImage **, int, VimosExtractionTable *, VimosTable **, int);

PIL_END_DECLS

#endif /* VM_MOSEXTRACTION_H */
