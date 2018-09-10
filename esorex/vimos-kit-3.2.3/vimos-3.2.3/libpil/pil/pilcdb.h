/* $Id: pilcdb.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifndef _PIL_CDB_H
#define _PIL_CDB_H

#include <stdio.h>

#include <pilmacros.h>
#include <pildictionary.h>


PIL_BEGIN_DECLS

/*
 * Keyword access modes
 */

typedef enum _PIL_CDB_KEYMODE_ {
  READ_WRITE = 0,
  READ_ONLY
} PilCdbKeyMode;


/*
 * Case sensitivity flag for keyword comparison
 */

typedef enum _PIL_CONFIG_KEYCASE_ {
  IGNORE_CASE = 0,
  USE_CASE
} PilCdbKeyCase;


/*
 * The opaque database object type
 */

typedef struct _PIL_CDB_ PilCdb;


/*
 * Database constructor and destructor
 */

extern PilCdb *newPilCdb(void);
extern void deletePilCdb(PilCdb *db);


/*
 * Configuration file I/O
 */

extern int pilCdbParseFile(PilCdb *, FILE *);
extern int pilCdbDumpDB(PilCdb *, FILE *);
extern char **pilCdbDumpDBtoString(PilCdb *, int *);


/*
 * Methods
 */

extern PilCdbKeyCase pilCdbGetKeyCase(const PilCdb *);
extern int pilCdbSetKeyCase(PilCdb *, PilCdbKeyCase);

extern char pilCdbGetGroupIFS(const PilCdb *);
extern int pilCdbSetGroupIFS(PilCdb *, char );

extern int pilCdbGroupExists(const PilCdb *, const char *);
extern int pilCdbEntryExists(const PilCdb *, const char *, const char *);

extern PilCdbKeyMode pilCdbGetKeyMode(const PilCdb *, const char *,
				      const char *);
extern int pilCdbSetKeyMode(PilCdb *, const char *, const char *, 
			    PilCdbKeyMode);

extern int pilCdbCreateGroup(PilCdb *, const char *);
extern int pilCdbCreateEntry(PilCdb *, const char *, const char *, 
			     const char *);

extern int pilCdbModifyValue(PilCdb *, const char *, const char *, 
			     const char *);

extern const char *pilCdbGetString(const PilCdb *, const char *,
				   const char *);
extern int pilCdbGetBool(const PilCdb *, const char *,const char *, int);
extern int pilCdbGetInt(const PilCdb *, const char *, const char *, int);
extern long pilCdbGetLong(const PilCdb *, const char *, const char *, long);
extern float pilCdbGetFloat(const PilCdb *, const char *, const char *,
			    float);
extern double pilCdbGetDouble(const PilCdb *, const char *, const char *,
			      double);

PIL_END_DECLS

#endif /* _PIL_CDB_H */
