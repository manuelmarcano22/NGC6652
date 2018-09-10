/* $Id: pilalias.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
 
#ifndef _PIL_ALIAS_H
#define _PIL_ALIAS_H

#include <pilmacros.h>

 
PIL_BEGIN_DECLS

typedef struct _PIL_ALIAS_ PilAlias;


/*
 * Alias constructor and destructor
 */

extern PilAlias *newPilAlias(const char *, const char *, const char *,
			     const char *);

extern void deletePilAlias(PilAlias *);


/*
 * Methods
 */

extern const char *pilAliasGetName(PilAlias *);
extern int pilAliasSetName(PilAlias *, const char *);

extern const char *pilAliasGetValue(PilAlias *);
extern int pilAliasSetValue(PilAlias *, const char *);

extern const char *pilAliasGetFormat(PilAlias *);
extern int pilAliasSetFormat(PilAlias *, const char *);

extern const char *pilAliasGetComment(PilAlias *);
extern int pilAliasSetComment(PilAlias *, const char *);

extern int pilAliasSet(PilAlias *, const char *, const char *, const char *,
		       const char *);
extern void pilAliasClear(PilAlias *);

PIL_END_DECLS

#endif /* _PIL_ALIAS_H */
