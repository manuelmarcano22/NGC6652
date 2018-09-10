/* $Id: piltranslator.h,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifndef _PIL_TRANSLATOR_H_
#define _PIL_TRANSLATOR_H_

#include <stdarg.h>

#include <pilmacros.h>
#include <pilkeymap.h>
#include <pilcatmap.h>


PIL_BEGIN_DECLS

/*
 * Keyword mapping
 */

int pilTrnInitKeywordMap(void);
int pilTrnClearKeywordMap(void);
int pilTrnLoadKeywordMap(const char *);

int pilTrnAddKey(const char *name, const char *value,
		 const char *format, const char *comment);

const char *pilKeyTranslate(const char *, ...);
const char *pilTrnGetKeyword(const char *alias, ...);
const char *pilTrnGetComment(const char *alias);


/*
 * Frame category mapping
 */

int pilTrnInitCategoryMap(void);
int pilTrnClearCategoryMap(void);
int pilTrnLoadCategoryMap(const char *);

const char *pilTrnGetCategory(const char *);
int pilTrnAddCategory(const char *name, const char *value);

const char *pilCatTranslate(const char *);

PIL_END_DECLS

#endif /* _PIL_TRANSLATOR_H_ */