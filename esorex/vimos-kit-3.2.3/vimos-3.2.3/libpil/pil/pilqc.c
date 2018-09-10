/* $Id: pilqc.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
#include <string.h>
#include <assert.h>

#include "pilmemory.h"
#include "pillist.h"
#include "pildate.h"
#include "pilerrno.h"
#include "pilpaf.h"
#include "pilqc.h"

#define INSTRUMENT          "[VIMOS]"
#define DICT_LINE_LENGTH    (80)
#define MAX_PAF_NAME_LENGTH (80)
#define PAF_ROOT_NAME       "qc"

/** 
 * @defgroup pilQC QC1 PAF File Utilities
 *
 * TBD
 */

/**@{*/

static PilPAF *pafFile = NULL;
static char    pafName[MAX_PAF_NAME_LENGTH];
static int     pafIndex = 0;

/*
 *  Returns version of QC1 dictionary, as indicated in its "Revision:" field,
 *  or NULL in case of failure (i.e. dictionary or "Revision:" field are not 
 *  found).
 */

static char *
_pilQcDictVersion(const char *name)
{

    static char      version[DICT_LINE_LENGTH];
    char             line[DICT_LINE_LENGTH];
    FILE            *fp;


    if (version[0] != '\0')
        return version;

    if (!(fp = fopen(name, "r")))
        return NULL;

    while(fgets(line, DICT_LINE_LENGTH, fp)) {
        if (line == strstr(line, "Revision:")) {
            sscanf(line, "Revision: %s", version);
            break;
        }
    }

    fclose(fp);

    if (line)
      return version;

    return NULL;

}


/**
 * @brief
 *   Initiate a new QC1 group.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE.
 *
 * A new PAF object is initiated, with the name depending on the root
 * PAF_ROOT_NAME and the current pafIndex. If the previous QC1 PAF file
 * it's found open, this is an error: pilQcGroupEnd() should be called 
 * first.
 */

int
pilQcGroupStart()
{

    if (pafFile) 
        return EXIT_FAILURE;

    sprintf(pafName, "%s%.4d.paf", PAF_ROOT_NAME, pafIndex);

    if (!(pafFile = newPilPAF(pafName, "QC1 parameters", NULL, NULL)))
        return EXIT_FAILURE;

   /* FIXME:
    * Here a call to _pilQcDictVersion() should be used, to get the
    * dictionary version.
    */

    pilQcWriteString("QC.DID", "VIMOS_QC-1.0", "QC1 dictionary");

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Close current QC1 PAF file.
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE.
 *
 * The current QC1 PAF object is written to disk file. If no PAF object
 * is present, this is an error: pilQcGroupStart() should be called
 * first. If the PAF file is empty, the PAF object is destroyed, but
 * no PAF file is created.
 */

int
pilQcGroupEnd()
{

    if (!pafFile)
        return EXIT_FAILURE;

    if (!pilPAFIsEmpty(pafFile)) {
        pilPAFWrite(pafFile);
        pafIndex++;
    }

    deletePilPAF(pafFile);
    pafFile = NULL;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Add string parameter to current QC1 group.
 *
 * @param name     Parameter name
 * @param value    Parameter value
 * @param comment  Parameter comment;
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE.
 *
 * To the current QC1 PAF object is appended a string parameter. If no 
 * PAF object is present, this is an error: pilQcGroupStart() should be 
 * called first.
 */

int
pilQcWriteString(const char *name, const char *value, const char *comment)
{

    int   status = EXIT_FAILURE;
    int   length = strlen(INSTRUMENT) + 1;
    char *allComment;

    assert(comment != 0x0);

    length += strlen(comment) + 1;

    if (!(allComment = pil_malloc(length * sizeof(char))))
      return status;

    sprintf(allComment, "%s %s", comment, INSTRUMENT);

    status = pilPAFAppendString(pafFile, name, value, allComment);

    pil_free(allComment);

    return status;

}


/**
 * @brief
 *   Add double parameter to current QC1 group.
 *
 * @param name     Parameter name
 * @param value    Parameter value
 * @param unit     Parameter unit
 * @param comment  Parameter comment;
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE.
 *
 * To the current QC1 PAF object is appended a double parameter. 
 * The comment string is mandatory. The parameter unit must be 
 * specified, unless the specified parameter is adimensional, otherwise
 * a null pointer should be passed. To the comment string the unit 
 * string (if present) will be appended, enclosed in round brackets, 
 * and then the string "[VIMOS]". If no PAF object is present, this 
 * is an error: pilQcGroupStart() should be called first.
 */

int
pilQcWriteDouble(const char *name, double value, const char *unit, 
              const char *comment)
{

    int   status = EXIT_FAILURE;
    int   length = strlen(INSTRUMENT) + 1;
    char *allComment;

    assert(comment != 0x0);

    length += strlen(comment) + 1;

    if (unit)
      length += strlen(unit) + 3;

    if (!(allComment = pil_malloc(length * sizeof(char))))
      return status;

    if (unit)
      sprintf(allComment, "%s (%s) %s", comment, unit, INSTRUMENT);
    else
      sprintf(allComment, "%s %s", comment, INSTRUMENT);

    status = pilPAFAppendDouble(pafFile, name, value, allComment);

    pil_free(allComment);

    return status;

}


/**
 * @brief
 *   Add integer parameter to current QC1 group.
 *
 * @param name     Parameter name
 * @param value    Parameter value
 * @param unit     Parameter unit
 * @param comment  Parameter comment;
 *
 * @return EXIT_SUCCESS or EXIT_FAILURE.
 *
 * To the current QC1 PAF object is appended a double parameter. If no
 * PAF object is present, this is an error: pilQcGroupStart() should be
 * called first.
 */

int
pilQcWriteInt(const char *name, int value, const char *unit, 
              const char *comment)
{

    int   status = EXIT_FAILURE;
    int   length = strlen(INSTRUMENT) + 1;
    char *allComment;

    assert(comment != 0x0);

    length += strlen(comment) + 1;

    if (unit)
      length += strlen(unit) + 3;

    if (!(allComment = pil_malloc(length * sizeof(char))))
      return status;

    if (unit)
      sprintf(allComment, "%s (%s) %s", comment, unit, INSTRUMENT);
    else
      sprintf(allComment, "%s %s", comment, INSTRUMENT);

    status = pilPAFAppendInt(pafFile, name, value, allComment);

    pil_free(allComment);

    return status;

}
/**@}*/
