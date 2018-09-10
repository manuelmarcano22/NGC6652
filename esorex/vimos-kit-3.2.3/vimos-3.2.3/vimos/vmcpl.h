/* $Id: vmcpl.h,v 1.4 2013-08-23 10:22:37 cgarcia Exp $
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
 * $Date: 2013-08-23 10:22:37 $
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VM_CPL_H
#define VM_CPL_H

#include <pilmacros.h>
#include <pilframeset.h>
#include <pildfsconfig.h>
#include <pilmessages.h>
#include <piltimer.h>

#include <cpl_frameset.h>
#include <cpl_parameterlist.h>
#include <cpl_msg.h>


PIL_BEGIN_DECLS

int vmCplFramesetImport(cpl_frameset *, const PilSetOfFrames *);
int vmCplFramesetExport(const cpl_frameset *, PilSetOfFrames *);

int vmCplParlistExport(const cpl_parameterlist *);

PilMsgSeverity vmCplMsgSeverityExport(cpl_msg_severity);

/*
 * Recipe timer functions
 */

int vmCplRecipeTimerStart(PilTime *);
int vmCplRecipeTimerStop(PilTime *);

/*
 * High-level functions
 */

int vmCplRecipeStart(const char *, const char *);
int vmCplRecipeStop(void);

/*
 * Miscellaneous utilities
 */

int vmCplPostProcessFrames(PilSetOfFrames *, const char *);

PIL_END_DECLS

#endif /* VM_CPL_H */
