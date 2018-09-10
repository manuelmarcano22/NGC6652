/* $Id: vmdetector.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_DETECTOR_H
#define VM_DETECTOR_H

#include <pilmacros.h>

#include <vmmath.h>
#include <vmfit.h>
#include <vmimage.h>
#include <vmtable.h>


PIL_BEGIN_DECLS

#define  MAX_RDPORTS    (9)

typedef struct _VIMOS_WINDOW_
{
  int            startX;
  int            startY;
  int            nX;
  int            nY;
  struct _VIMOS_WINDOW_ *prev;
  struct _VIMOS_WINDOW_ *next;
} VimosWindow;

typedef struct _VIMOS_PORT_
{
  VimosWindow         *window;
  VimosWindow         *prScan;
  VimosWindow         *ovScan;
  VimosWindow         *readOutWindow;
  int                  shiftX;    /* Coordinate shift with respect lower */
  int                  shiftY;    /* left corner of chip                 */
  struct _VIMOS_PORT_ *prev;
  struct _VIMOS_PORT_ *next;
} VimosPort;

/*
 * Constructor of VimosWindow
 */

VimosWindow *newWindow();

/*
 * Destructor of VimosWindow                                               */

void deleteWindow(VimosWindow *window);

/*
 * Truncating a VimosWindow list at any point (including beginning)
 */

void deleteWindowList(VimosWindow *window);

/*
 * Constructor of VimosPort
 */

VimosPort *newPort();

/*
 * Destructor of VimosPort
 */

void deletePort(VimosPort *port);

/*
 * Truncating a VimosPort list at any point (including beginning)
 */

void deletePortList(VimosPort *port);

/*
 * This one to get list of ports, readout windows, and pre/overscan
 * regions, from any image header
 */

VimosPort  *getPorts(VimosImage *image, int *nports);

int getTotalReadoutWindow(VimosPort *, int *, int *, int *, int *);

/*
 * Estimate average bias for each port, from overscan regions.
 */

VimosFloatArray *estimateImageBias(VimosImage *, VimosPort *port);

/*
 * Subtract bias from the overscan regions of a given image.
 * (also from the overscan regions themselves).
 */

VimosBool subtractOverscan(float a[], int nx, int ny, VimosPort *port);

/*
 * Get scan direction (1 = vertical, 0 = horizontal)
 */
int getReadoutDirection(VimosPort *port);

/*
 * Get chip size
 */
int getChipSize(VimosImage *image, int *chipSizeX, int *chipSizeY);

/*
 * Compute RON from an image overscans.
 */

VimosFloatArray *estimateImageRon(VimosImage *, VimosPort *);

/*
 * Read RON from an image keywords header.
 */

VimosFloatArray *getImageRon(VimosImage *);

PIL_END_DECLS

#endif /* VM_DETECTOR_H */
