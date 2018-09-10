/* $Id: pildfsconfig.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
 *
 * This file is part of the VIMOS pipeline library
 * Copyright (C) 2000-2004 European Southern Observatory
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: cizzo $
 * $Date: 2008-10-21 09:10:13 $
 * $Revision: 1.1.1.1 $
 * $Name: not supported by cvs2svn $
 */

#ifndef _PIL_DFSCONFIG_H
#define _PIL_DFSCONFIG_H

#include <pilmacros.h>
#include <pilcdb.h>


PIL_BEGIN_DECLS

/*
 * Constructor and destructor for the pipeline configuration database
 */

extern int pilDfsCreateDB(int, PilCdbKeyCase);
extern void pilDfsFreeDB(void);


/*
 * Methods
 */

extern int pilDfsReadSetupFiles(const char *, const char *);
extern int pilDfsDumpDB(const char *);
extern char **pilDfsDumpDBtoString(int *);
extern int pilDfsGetEnv(void);

extern int pilDfsDbGroupExists(const char *);
extern int pilDfsDbEntryExists(const char *, const char *);

extern PilCdbKeyMode pilDfsDbGetKeyMode(const char *, const char *);

extern int pilDfsDbCreateGroup(const char *);
extern int pilDfsDbCreateEntry(const char *, const char *, const char *,
			       PilCdbKeyMode);

extern int pilDfsDbModifyValue(const char *, const char *, const char *);

extern const char *pilDfsDbGetString(const char *, const char *);
extern int pilDfsDbGetBool(const char *, const char *, int);
extern int pilDfsDbGetInt(const char *, const char *, int);
extern long pilDfsDbGetLong(const char *, const char *, long);
extern float pilDfsDbGetFloat(const char *, const char *, float);
extern double pilDfsDbGetDouble(const char *, const char *, double);

PIL_END_DECLS

#endif /* _PIL_DFSCONFIG_H */
