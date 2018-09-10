/* $Id: pilmemory.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>

#include "pilmemory.h"


/*
 * Define our own version of realloc() in case NULL cannot be passed
 * to realloc() as first argument.
 */

#ifdef HAVE_WORKING_REALLOC
#  define pil_memory_realloc realloc
#else /* !HAVE_WORKING_REALLOC */

static void *pil_memory_realloc(void *memory, size_t nbytes)
{
    if (!memory)
        return malloc(nbytes);
    else
        return realloc(memory, nbytes);
}

#endif /* !HAVE_WORKING_REALLOC */


/*
 * Internal error message utility
 */

 static void
pil_memory_print_error(const char *pos, size_t nbytes)
{

    fprintf(stderr, "%s: failed to allocate %d bytes", pos, nbytes);
    fprintf(stderr, "\naborting...\n");
    abort();

}


/*
 * Allocate `nbytes' bytes.
 */

void *
pil_malloc(size_t nbytes)
{

    if (nbytes) {
        void *mblk = malloc(nbytes);

        if (mblk)
            return mblk;

        pil_memory_print_error(PIL_CODE_POS, nbytes);
    }

    return NULL;

}


/*
 * Allocate `nbytes' bytes and clear them
 */

void *
pil_malloc_clear(size_t nbytes)
{

    if (nbytes) {
        void *mblk = calloc(1, nbytes);

        if (mblk)
            return mblk;

        pil_memory_print_error(PIL_CODE_POS, nbytes);
    }

    return NULL;

}


/*
 * Allocate memory for `natoms' elements of size `size'.
 */

void *
pil_calloc(size_t natoms, size_t nbytes)
{

   if (natoms && nbytes) {
        void *mblk = calloc(natoms, nbytes);

        if (mblk)
            return mblk;

        pil_memory_print_error(PIL_CODE_POS, natoms * nbytes);
    }

    return NULL;

}


/*
 * Change the size of a memory block
 */

void *
pil_realloc(void *memory, size_t nbytes)
{

    if (nbytes) {
        void *mblk = pil_memory_realloc(memory, nbytes);

        if (mblk)
            return mblk;

        pil_memory_print_error(PIL_CODE_POS, nbytes);
    }

    return NULL;

}


/*
 * Memory block deallocation
 */

void
pil_free(void *memory)
{

    if (memory) {
        free(memory);
    }

    return;

}
