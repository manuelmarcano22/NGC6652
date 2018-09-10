/* $Id: pilpaf.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
 * Foundation,  Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
#include <ctype.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/types.h>
#include <assert.h>

#include "pilmemory.h"
#include "pilstrutils.h"
#include "pillist.h"
#include "pildate.h"
#include "pilerrno.h"
#include "pilpaf.h"


/*
 * Reserved PAF header keywords
 */

#define PAF_HDR_START      "PAF.HDR.START"
#define PAF_TYPE           "PAF.TYPE"
#define PAF_ID             "PAF.ID"
#define PAF_NAME           "PAF.NAME"
#define PAF_DESC           "PAF.DESC"
#define PAF_CRTE_NAME      "PAF.CRTE.NAME"
#define PAF_CRTE_TIME      "PAF.CRTE.DAYTIM"
#define PAF_LCHG_NAME      "PAF.LCHG.NAME"
#define PAF_LCHG_TIME      "PAF.LCHG.DAYTIM"
#define PAF_CHCK_NAME      "PAF.CHCK.NAME"
#define PAF_CHCK_TIME      "PAF.CHCK.DAYTIM"
#define PAF_CHCK_CHECKSUM  "PAF.CHCK.CHECKSUM"
#define PAF_HDR_END        "PAF.HDR.END"


/*
 * Value and comment field start position
 */

#define PAF_FIELD_OFFSET_VALUE    20
#define PAF_FIELD_OFFSET_COMMENT  45


/** 
 * @defgroup pilPaf PAF File Utilities
 *
 * TBD
 */

/**@{*/

/*
 * PAF record definition. This corresponds to one line of a parameter file.
 */

struct _PIL_PAF_RECORD_ {
    char *name;
    char *comment;
    PilPAFType type;

    union {
        int *bval;
        int *ival;
        double *dval;
        char *sval;
    } data;
};

typedef struct _PIL_PAF_RECORD_ PilPAFRecord;


/*
 * The PAF object
 */

struct _PIL_PAF_ {
    PilList *header;
    PilList *records;
};


/*
 * Compute record type size in bytes.
 */

 static size_t
_pilPAFValueSize(PilPAFType type, const void *value)
{
    size_t sz;

    switch (type) {
        case PAF_TYPE_BOOL:
            sz = sizeof(int);
            break;
            
        case PAF_TYPE_INT:
            sz = sizeof(int);
            break;

        case PAF_TYPE_DOUBLE:
            sz = sizeof(double);
            break;

        case PAF_TYPE_STRING:
            sz = (strlen((char *)value) + 1) * sizeof(char);
            break;

        default:
                sz = 0;
                break;
        }

    return sz;

}


/*
 * Comparison function for looking up records by their name.
 */

 static int
_pilPAFRecordCompare(const void *record, const void *name)
{

    PilPAFRecord *r = (PilPAFRecord *)record;

    assert(name != NULL);
    assert(record != NULL);

    return strcmp(r->name, name);

}


/*
 * Destroy a PAF record
 */

 static void
_pilPAFRecordDestroy(PilPAFRecord *record)
{

    if (record) {
        if (record->name)
            pil_free(record->name);

        if ((void *)record->data.sval)
            pil_free((void *)record->data.sval);

        if (record->comment)
            pil_free(record->comment);

        pil_free(record);
    }

    return;

}


/*
 * Create a new PAF record
 */

 static PilPAFRecord *
_pilPAFRecordCreate(const char *name, PilPAFType type, const void *value,
                    const char *comment)
{

    size_t sz;

    PilPAFRecord *record = (PilPAFRecord *)pil_malloc(sizeof *record);


    if (record) {
        record->name = pil_strdup(name);
        record->comment = comment ? pil_strdup(comment) : NULL;
        record->type = type;

        sz = _pilPAFValueSize(type, value);

        if (sz == 0) {
            record->data.sval = NULL;
        }
        else {
            record->data.sval = (char *)pil_malloc(sz);
            if (!record->data.sval) {
                _pilPAFRecordDestroy(record);
                return NULL;
            }
        }

        memcpy(record->data.sval, value, sz);
    }

    return record;

}

/*
 * Set name, value and comment of a PAF record
 */

 static void
_pilPAFRecordSet(PilPAFRecord *record, const char *name, PilPAFType type,
                 const void *value, const char *comment)
{

    if (name) {
        if (record->name)
            pil_free(record->name);
        record->name = pil_strdup(name);
    }

    if (comment) {
        if (record->comment)
            pil_free(record->comment);
        record->comment = pil_strdup(comment);
    }

    if (value) {
        size_t sz = _pilPAFValueSize(type, value);

        if (record->data.sval) {
            size_t size = _pilPAFValueSize(record->type, record->data.sval);

            if (sz != size)
                record->data.sval = (char *)pil_realloc(record->data.sval, sz);
        }
        else
            record->data.sval = (char *)pil_malloc(sz);

        memcpy(record->data.sval, value, sz);
        record->type = type;
    }

    return;

}


/*
 * Create a record from name, value and comment and insert the record
 * into the PAF header or record list at the position given by key
 */

 static int
_pilPAFInsert(PilList *list, const char *key, const char *name,
              PilPAFType type, const void *value, const char *comment)
{

    PilPAFRecord *record;
    PilListNode *node, *position;


    position = pilListLookup(list, key, _pilPAFRecordCompare);
    if (!position)
        return 1;

    record = _pilPAFRecordCreate(name, type, value, comment);
    if (!record)
        return 1;

    node = newPilListNode(record);
    if (!node)
        return 1;

    pilListInsert(list, position, node);

    return 0;

}


/*
 * Create a record from name, value and comment and insert the record
 * into the PAF header or record list after the position given by key
 */

 static int
_pilPAFInsertAfter(PilList *list, const char *key, const char *name,
                   PilPAFType type, const void *value, const char *comment)
{

    PilPAFRecord *record;
    PilListNode *node, *position;


    position = pilListLookup(list, key, _pilPAFRecordCompare);
    if (!position)
        return 1;

    record = _pilPAFRecordCreate(name, type, value, comment);
    if (!record)
        return 1;

    node = newPilListNode(record);
    if (!node)
        return 1;

    pilListInsertAfter(list, position, node);

    return 0;

}


/*
 * Create a record from name, value and comment and prepend the record to
 * the PAF header or record list.
 */

 static int
_pilPAFPrepend(PilList *list, const char *name, PilPAFType type,
               const void *value, const char *comment)
{

    PilListNode *node;
    PilPAFRecord *record;


    record = _pilPAFRecordCreate(name, type, value, comment);
    if (!record)
        return 1;

    node = newPilListNode(record);
    if (!node)
        return 1;

    pilListPushFront(list, node);

    return 0;

}


/*
 * Create a record from name, value and comment and append the record to
 * the PAF header or record list.
 */

 static int
_pilPAFAppend(PilList *list, const char *name, PilPAFType type,
              const void *value, const char *comment)
{

    PilListNode *node;
    PilPAFRecord *record;


    record = _pilPAFRecordCreate(name, type, value, comment);
    if (!record)
        return 1;

    node = newPilListNode(record);
    if (!node)
        return 1;

    pilListPushBack(list, node);

    return 0;

}


/*
 * Create a new PAF header.
 */

 static PilList *
_pilPAFHeaderCreate(const char *name, const char *type, const char *id,
                    const char *desc)
{

    PilList *hdr = newPilList();

    if (hdr) {
        _pilPAFAppend(hdr, PAF_HDR_START, PAF_TYPE_NONE, NULL, NULL);
        _pilPAFAppend(hdr, PAF_TYPE, PAF_TYPE_STRING, type,
                      "Type of parameter file");

        if (id)
            _pilPAFAppend(hdr, PAF_ID, PAF_TYPE_STRING, id, NULL);
        else
            _pilPAFAppend(hdr, PAF_ID, PAF_TYPE_STRING, "", NULL);

        _pilPAFAppend(hdr, PAF_NAME, PAF_TYPE_STRING, name, "Name of PAF");

        if (desc)
            _pilPAFAppend(hdr, PAF_DESC, PAF_TYPE_STRING, desc,
                          "Short description of PAF");
        else
            _pilPAFAppend(hdr, PAF_DESC, PAF_TYPE_STRING, "",
                          "Short description of PAF");

        _pilPAFAppend(hdr, PAF_CRTE_NAME, PAF_TYPE_NONE, NULL,
                      "Name of creator");
        _pilPAFAppend(hdr, PAF_CRTE_TIME, PAF_TYPE_NONE, NULL,
                      "Civil time for creation");
        _pilPAFAppend(hdr, PAF_LCHG_NAME, PAF_TYPE_NONE, NULL,
                      "Author of par. file");
        _pilPAFAppend(hdr, PAF_LCHG_TIME, PAF_TYPE_NONE, NULL,
                      "Timestamp for last change");
        _pilPAFAppend(hdr, PAF_CHCK_NAME, PAF_TYPE_STRING, "",
                      "Name of appl. checking");
        _pilPAFAppend(hdr, PAF_CHCK_TIME, PAF_TYPE_STRING, "",
                      "Time for checking");
        _pilPAFAppend(hdr, PAF_CHCK_CHECKSUM, PAF_TYPE_STRING, "",
                      "Checksum for the PAF");
        _pilPAFAppend(hdr, PAF_HDR_END, PAF_TYPE_NONE, NULL, NULL);
    }

    return hdr;

}


/*
 * Get the value for a given PAF record name.
 */

 static int
_pilPAFGetValueBool(PilList *list, const char *name)
{

    PilPAFRecord *record;
    PilListNode *node = pilListLookup(list, name, _pilPAFRecordCompare);

    if (!node) {
        pilErrno = PIL_ENODATA;
        return 0;
    }

    record = pilListNodeGet(node);
    if (record->type != PAF_TYPE_BOOL) {
        pilErrno = PIL_EBADTYPE;
        return 0;
    }

    return *record->data.bval;

}


 static int
_pilPAFGetValueInt(PilList *list, const char *name)
{

    PilPAFRecord *record;
    PilListNode *node = pilListLookup(list, name, _pilPAFRecordCompare);

    if (!node) {
        pilErrno = PIL_ENODATA;
        return 0;
    }

    record = pilListNodeGet(node);
    if (record->type != PAF_TYPE_INT) {
        pilErrno = PIL_EBADTYPE;
        return 0;
    }

    return *record->data.ival;

}


 static double
_pilPAFGetValueDouble(PilList *list, const char *name)
{

    PilPAFRecord *record;
    PilListNode *node = pilListLookup(list, name, _pilPAFRecordCompare);

    if (!node) {
        pilErrno = PIL_ENODATA;
        return 0.;
    }

    record = pilListNodeGet(node);
    if (record->type != PAF_TYPE_DOUBLE) {
        pilErrno = PIL_EBADTYPE;
        return 0.;
    }

    return *record->data.dval;

}


 static const char *
_pilPAFGetValueString(PilList *list, const char *name)
{

    PilPAFRecord *record;
    PilListNode *node = pilListLookup(list, name, _pilPAFRecordCompare);

    if (!node) {
        pilErrno = PIL_ENODATA;
        return NULL;
    }

    record = pilListNodeGet(node);
    if (record->type != PAF_TYPE_STRING) {
        pilErrno = PIL_EBADTYPE;
        return NULL;
    }

    return record->data.sval;

}


/*
 * Modify the value for a given PAF record name.
 */

 static int
_pilPAFSetValue(PilList *list, const char *name, PilPAFType type,
                const void *value)
{

    size_t sz = _pilPAFValueSize(type, value);

    PilPAFRecord *record;
    PilListNode *node;


    /*
     * Invalid type
     */

    if (sz == 0)
        return 1;

    node = pilListLookup(list, name, _pilPAFRecordCompare);
    if (!node) {
        pilErrno = PIL_ENODATA;
        return 1;
    }

    record = pilListNodeGet(node);

    if (record->data.sval) {
        if (record->type == type) {
            size_t size = _pilPAFValueSize(record->type, record->data.sval);

            if (sz != size) {
                record->data.sval = (char *)pil_realloc(record->data.sval, sz);
                if (!record->data.sval)
                    return 1;
            }
        }
        else {
            pilErrno = PIL_EBADTYPE;
            return 1;
        }
    }
    else {
        record->data.sval = (char *)pil_malloc(sz);
        if (!record->data.sval)
            return 1;
    }
    memcpy(record->data.sval, value, sz);
       
    record->type = type;

    return 0;

}


/*
 * Format a record so that it can be written to a parameter file on disk.
 * The formatted record is written to the given output buffer.
 */

 static const char *
_pilPAFFormatRecord(PilPAFRecord *record)
{

    static char buffer[PAF_RECORD_MAX + 1];
    char value[PAF_RECORD_MAX + 1];

    int pos, sz;


    memset(buffer, ' ', PAF_RECORD_MAX);


    /*
     * Verify that the record name fits into the buffer. The extra
     * character is for the semicolon which has to be present.
     */

    if (strlen(record->name) + 1 > PAF_RECORD_MAX)
        return NULL;


    /*
     * Build the formatted string from the record structure
     */

    sz = strlen(record->name);
    strncpy(buffer, record->name, sz);

    pos = sz;
    if (record->data.sval) {
        if (pos < PAF_FIELD_OFFSET_VALUE)
            pos = PAF_FIELD_OFFSET_VALUE;
        else
            pos++;

        switch (record->type) {
            case PAF_TYPE_BOOL:
                snprintf(value, PAF_RECORD_MAX, "%c",
                         *record->data.bval ? 'T' : 'F');
                break;

            case PAF_TYPE_INT:
                snprintf(value, PAF_RECORD_MAX, "%d", *record->data.ival);
                break;

            case PAF_TYPE_DOUBLE:
                snprintf(value, PAF_RECORD_MAX, "%.15G", *record->data.dval);
                if (!strchr(value, '.')) {
                    if (strchr(value, 'E'))
                        snprintf(value, PAF_RECORD_MAX, "%.1E",
                                 *record->data.dval);
                    else
                        strcat(value, ".");
                }
                break;

            case PAF_TYPE_STRING:
                snprintf(value, PAF_RECORD_MAX, "\"%s\"", record->data.sval);
                break;

            case PAF_TYPE_NONE:

                /* 
                 * Should not reach this point. If type is PAF_TYPE_NONE
                 * the data pointer should always be NULL.
                 */
                
                break;
        }

        sz = strlen(value);

        /* 
         * Verify that writing the value string does not overflow the buffer.
         */
        
        if (sz > PAF_RECORD_MAX - pos + 1)
            return NULL;

        strncpy(&buffer[pos], value, sz);
        pos += sz;
    }

    buffer[pos++] = ';';


    /*
     * Comments are not printed if there is room in the buffer for at least 3
     * characters, so that not only the hash and/or the following blank
     * could be stored because of the finite record size.
     */

    if (record->comment && (PAF_RECORD_MAX - pos) >= 2) {
        if (pos < PAF_FIELD_OFFSET_COMMENT)
            pos = PAF_FIELD_OFFSET_COMMENT;
        else
            pos++;

        strncpy(&buffer[pos], "# ", 2);
        pos += 2;
        sz = strlen(record->comment);
        strncpy(&buffer[pos], record->comment, sz);
        pos += sz;
    }

    buffer[pos] = '\0';
    
    return buffer;
}


/**
 * @brief
 *   Destroy a PAF object.
 *
 * @param paf  PAF object to destroy.
 *
 * @return Nothing.
 *
 * The function deallocates all memory used for the PAF object @em paf.
 */

 void
deletePilPAF(PilPAF *paf)
{

    if (paf) {
        pilListDestroy(paf->header, (void (*)(void *))_pilPAFRecordDestroy);
        pilListDestroy(paf->records, (void (*)(void *))_pilPAFRecordDestroy);

        pil_free(paf);
    }

    return;

}


/**
 * @brief
 *   Create a new PAF object.
 *
 * @param name  Parameter file name.
 * @param type  Parameter file type.
 * @param id    Parameter file identifier string.
 * @param desc  Short description for the parameter file.
 *
 * @return The handle for the newly created PAF object.
 *
 * The function allocates the memory for a PAF object and initializes the
 * PAF object with the strings @em name, @em type, @em id and @em desc passed
 * to the function, where @em id and @em desc may be omitted, i.e. @c NULL
 * maybe passed for @em id and/or @em desc. The argument @em name will be
 * used as the file name in case the PAF object is written to a disk file.
 */

PilPAF *
newPilPAF(const char *name, const char *type, const char *id,
          const char *desc)
{

    PilPAF *paf;


    if (!name || !type)
        return NULL;

    paf = (PilPAF *)pil_malloc(sizeof(PilPAF));
    if (paf) {
        paf->header = _pilPAFHeaderCreate(name, type, id, desc);
        paf->records = newPilList();

        if (!paf->header || !paf->records) {
            deletePilPAF(paf);
            return NULL;
        }
    }

    return paf;

}


/**
 * @brief
 *   Check whether a PAF object is empty.
 *
 * @param paf  PAF object to query.
 *
 * @return The function returns 1 if @em paf is empty, otherwise the return
 *   value is 0.
 *
 * The function checks whether the record list of the PAF object @em paf
 * contains any records.
 */

int
pilPAFIsEmpty(const PilPAF *paf)
{

    assert(paf != NULL);
    assert(paf->records != NULL);

    return pilListIsEmpty(paf->records);

}


/**
 * @brief
 *   Get the actual size of the given PAF object.
 *
 * @param paf  PAF object to query.
 *
 * @return The number of PAF records currently stored in the record list, or
 *   0 if the PAF object is empty.
 *
 * The function computes the actual size, i.e. the number of parmeter records,
 * stored in the record list of the PAF object @em paf. Header entries do not
 * contribute the size of a PAF object.
 */

size_t
pilPAFGetSize(const PilPAF *paf)
{
    assert(paf != NULL);
    assert(paf->records != NULL);

    return (size_t)pilListSize(paf->records);

}


/**
 * @brief
 *   Check whether a PAF object contains a parameter with the given name.
 *
 * @param paf   PAF object to query.
 * @param name  Parameter name to look for.
 *
 * @return The function returns 1 if the name is present in the record list,
 *   or 0 if it was not found.
 *
 * The function searches the record list of the PAF object @em paf for an
 * entry with the parameter name @name and returns 1 in case such a record
 * was found. The PAF header is not searched.
 */

int
pilPAFContains(const PilPAF *paf, const char *name)
{

    PilListNode *node;


    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    node = pilListBegin(paf->records);
    while (node) {
        if (_pilPAFRecordCompare(pilListNodeGet(node), name) == 0)
            return 1;

        node = pilListNext(paf->records, node);
    }

    return 0;

}


/**
 * @brief
 *   Get the number of records with the given name.
 *
 * @param paf   PAF object to query.
 * @param name  Parameter name to look for.
 *
 * @return The number of PAF records stored in the record list with
 *   the parameter name @em name.
 *
 * The function searches the record list of the PAF object @em paf for all
 * entries with the parameter name @name and returns the number of entries
 * found.
 */

size_t
pilPAFCount(const PilPAF *paf, const char *name)
{

    size_t count = 0;

    PilListNode *node;


    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    node = pilListBegin(paf->records);
    while (node) {
        if (_pilPAFRecordCompare(pilListNodeGet(node), name) == 0)
            ++count;

        node = pilListNext(paf->records, node);
    }

    return count;

}


/**
 * @brief
 *   Retrieve the name of a PAF object.
 *
 * @param paf  PAF object to be queried.
 *
 * @return The function returns the name of the PAF object @em paf.
 *
 * The function queries the PAF object @em paf for the name set in its
 * header and returns the string.
 */

const char *
pilPAFGetName(const PilPAF *paf)
{

    assert(paf != NULL);
    assert(paf->header != NULL);

    return _pilPAFGetValueString(paf->header, PAF_NAME);

}


/**
 * @brief
 *   Retrieve a PAF object's type tag.
 *
 * @param paf  PAF object to be queried.
 *
 * @return The function returns the type string of the PAF object @em paf.
 *
 * The function queries the PAF object @em paf for the type set in its
 * header and returns the string.
 */

const char *
pilPAFGetTag(const PilPAF *paf)
{

    assert(paf != NULL);
    assert(paf->header != NULL);

    return _pilPAFGetValueString(paf->header, PAF_TYPE);

}


/**
 * @brief
 *   Retrieve a PAF object's id tag.
 *
 * @param paf  PAF object to be queried.
 *
 * @return The function returns the id string of the PAF object @em paf.
 *
 * The function queries the PAF object @em paf for the id tag set in its
 * header and returns the string.
 */

const char *
pilPAFGetId(const PilPAF *paf)
{

    assert(paf != NULL);
    assert(paf->header != NULL);

    return _pilPAFGetValueString(paf->header, PAF_ID);

}


/**
 * @brief
 *   Retrieve a PAF object's short description.
 *
 * @param paf  PAF object to be queried.
 *
 * @return The function returns the short description of the PAF object
 *   @em paf.
 *
 * The function queries the PAF object @em paf for the short description field
 * set in its header and returns the string.
 */

const char *
pilPAFGetDescription(const PilPAF *paf)
{

    assert(paf != NULL);
    assert(paf->header != NULL);

    return _pilPAFGetValueString(paf->header, PAF_DESC);

}


/**
 * @brief
 *   Modify a PAF object's name field.
 *
 * @param paf   A PAF object.
 * @param name  New name string
 *
 * @return The function returns EXIT_SUCCESS if the name field was
 *   successfully modified, otherwise EXIT_FAILURE is returned.
 *
 * The function replaces the contents of the name field in the PAF header
 * with the string @em name.
 *
 * @note
 *   Changing the name of a PAF object implies a change of the file name
 *   used when the PAF object is written to disk.
 *
 * @see pilPAFWrite()
 */

int
pilPAFSetName(PilPAF *paf, const char *name)
{

    assert(paf != NULL);
    assert(paf->header != NULL);

    return _pilPAFSetValue(paf->header, PAF_NAME, PAF_TYPE_STRING, name);

}


/**
 * @brief
 *   Modify a PAF object's type tag.
 *
 * @param paf   A PAF object.
 * @param type  New type string
 *
 * @return The function returns EXIT_SUCCESS if the type field was
 *   successfully modified, otherwise EXIT_FAILURE is returned.
 *
 * The function replaces the contents of the type field in the PAF header
 * with the string @em type.
 */

int
pilPAFSetTag(PilPAF *paf, const char *type)
{

    assert(paf != NULL);
    assert(paf->header != NULL);

    return _pilPAFSetValue(paf->header, PAF_TYPE, PAF_TYPE_STRING, type);

}


/**
 * @brief
 *   Modify a PAF object's id string field.
 *
 * @param paf  A PAF object.
 * @param id   New id string
 *
 * @return The function returns EXIT_SUCCESS if the id string was
 *   successfully modified, otherwise EXIT_FAILURE is returned.
 *
 * The function replaces the contents of the id field in the PAF header
 * with the string @em id.
 */

int
pilPAFSetId(PilPAF *paf, const char *id)
{

    assert(paf != NULL);
    assert(paf->header != NULL);

    return _pilPAFSetValue(paf->header, PAF_ID, PAF_TYPE_STRING, id);

}


/**
 * @brief
 *   Modify a PAF object's short description.
 *
 * @param paf   A PAF object.
 * @param desc  New description.
 *
 * @return The function returns EXIT_SUCCESS if the short description was
 *   successfully modified, otherwise EXIT_FAILURE is returned.
 *
 * The function replaces the contents of the short description in the PAF
 * header with the string @em desc.
 */

int
pilPAFSetDescription(PilPAF *paf, const char *desc)
{

    assert(paf != NULL);
    assert(paf->header != NULL);

    return _pilPAFSetValue(paf->header, PAF_DESC, PAF_TYPE_STRING,
                           desc);

}


/**
 * @brief
 *   Modify the PAF header fields.
 *
 * @param paf   A PAF object.
 * @param name  New name.
 * @param type  New type.
 * @param id    New id string.
 * @param desc  New short description.
 *
 * @return The function returns EXIT_SUCCESS if the header was successfully
 *   modified, or EXIT_FAILURE otherwise.
 *
 * The function sets the header fields name, type, id and desc. If NULL
 * is passed for any of these arguments the corresponding header field is
 * left untouched.
 *
 * @note
 *   Changing the name of a PAF object implies a change of the file name
 *   used when the PAF object is written to disk.
 *
 * @see pilPAFWrite()
 */


int
pilPAFSetHeader(PilPAF *paf, const char *name, const char *type,
                const char *id, const char *desc)
{

    if (name)
        if (_pilPAFSetValue(paf->header, PAF_NAME, PAF_TYPE_STRING, name))
            return EXIT_FAILURE;

    if (type)
        if (_pilPAFSetValue(paf->header, PAF_TYPE, PAF_TYPE_STRING, type))
            return EXIT_FAILURE;

    if (id)
        if (_pilPAFSetValue(paf->header, PAF_ID, PAF_TYPE_STRING, id))
            return EXIT_FAILURE;

    if (desc)
        if (_pilPAFSetValue(paf->header, PAF_DESC, PAF_TYPE_STRING, desc))
            return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Verify that the given string is a valid PAF keyword.
 *
 * @param name  Name string to verify.
 *
 * @return The function returns 1 if the given name is valid, if not 0 is
 *   returned.
 *
 * The function does a simple check on the given string. To be a valid PAF
 * keyword the name should start with a letter, it should only contain
 * upper case letters and it may not contain blanks.
 */

 int
pilPAFIsValidName(const char *name)
{

    register size_t i, sz;


    assert(name != NULL);

    if (strchr(name, ' '))
        return 0;

    sz = strlen(name);
    for (i = 0; i <sz; i++) {
        char c = name[i];

        /*
         * Names may be composed from uppercase letters, digits, the dot
         * and the underscore only.
         */

        /*
         * Note: The characer class functions have to be enclosed in
         *   parantheses to use the actual function on HP-UX where these
         *   functions are also provided as macros, which are taken by
         *   default and may lead to compiler warnings.
         */

        if (!(isupper)(c) && !(isdigit)(c) && c != '.' && c != '_' && c != '-')
            return 0;
    }

    return 1;

}


/**
 * @brief
 *   Determine the type of a parameter value.
 *
 * @param paf   PAF object to query.
 * @param name  Parameter name.
 *
 * @return The function returns the type of the value associated with
 *   the parameter @em name.  If no such entry exists the returned type
 *   is @c PAF_TYPE_NONE and the error variable @b pilErrno is set to
 *   @c PIL_ENODATA.
 *
 * The function locates the first parameter with the name @em name in record
 * list of the PAF object @em paf and returns the type of the value.
 */

PilPAFType
pilPAFType(const PilPAF *paf, const char *name)
{

    PilPAFRecord *record;
    PilListNode *node;


    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    node = pilListLookup(paf->records, name, _pilPAFRecordCompare);
    if (!node) {
        pilErrno = PIL_ENODATA;
        return PAF_TYPE_NONE;
    }

    record = pilListNodeGet(node);

    return record->type;

}
    

/**
 * @brief
 *    Retrieve a boolean value from a PAF object.
 *
 * @param paf   A PAF object.
 * @param name  Parameter name to look up.
 *
 * @return The function returns the boolean value of the parameter with the
 *   name @em name. If a parameter with the name @em name does not exist,
 *   or if the parameter has not the expected type 0 is returned and the
 *   error variable @b pilErrno is set to @c PIL_ENODATA or @c PIL_EBADTYPE,
 *   respectively.
 *
 * The function looks for the first occurrence of a parameter with the name
 * @em name in the record list of the PAF object @em paf and returns its
 * value.
 */

int
pilPAFGetValueBool(const PilPAF *paf, const char *name)
{

    assert(paf != NULL);
    assert(paf->records != NULL);

    return _pilPAFGetValueBool(paf->records, name);

}


/**
 * @brief
 *    Retrieve a integer value from a PAF object.
 *
 * @param paf   A PAF object.
 * @param name  Parameter name to look up.
 *
 * @return The function returns the integer value of the parameter with the
 *   name @em name. If a parameter with the name @em name does not exist,
 *   or if the parameter has not the expected type 0 is returned and the
 *   error variable @b pilErrno is set to @c PIL_ENODATA or @c PIL_EBADTYPE,
 *   respectively.
 *
 * The function looks for the first occurrence of a parameter with the name
 * @em name in the record list of the PAF object @em paf and returns its
 * value.
 */

int
pilPAFGetValueInt(const PilPAF *paf, const char *name)
{

    assert(paf != NULL);
    assert(paf->records != NULL);

    return _pilPAFGetValueInt(paf->records, name);

}


/**
 * @brief
 *    Retrieve a double value from a PAF object.
 *
 * @param paf   A PAF object.
 * @param name  Parameter name to look up.
 *
 * @return The function returns the double value of the parameter with the
 *   name @em name. If a parameter with the name @em name does not exist,
 *   or if the parameter has not the expected type 0 is returned and the
 *   error variable @b pilErrno is set to @c PIL_ENODATA or @c PIL_EBADTYPE,
 *   respectively.
 *
 * The function looks for the first occurrence of a parameter with the name
 * @em name in the record list of the PAF object @em paf and returns its
 * value.
 */

double
pilPAFGetValueDouble(const PilPAF *paf, const char *name)
{

    assert(paf != NULL);
    assert(paf->records != NULL);

    return _pilPAFGetValueDouble(paf->records, name);

}


/**
 * @brief
 *    Retrieve a string value from a PAF object.
 *
 * @param paf   A PAF object.
 * @param name  Parameter name to look up.
 *
 * @return The function returns the string value of the parameter with the
 *   name @em name. If a parameter with the name @em name does not exist,
 *   or if the parameter has not the expected type 0 is returned and the
 *   error variable @b pilErrno is set to @c PIL_ENODATA or @c PIL_EBADTYPE,
 *   respectively.
 *
 * The function looks for the first occurrence of a parameter with the name
 * @em name in the record list of the PAF object @em paf and returns its
 * value.
 */

const char *
pilPAFGetValueString(const PilPAF *paf, const char *name)
{

    assert(paf != NULL);
    assert(paf->records != NULL);

    return _pilPAFGetValueString(paf->records, name);

}


/**
 * @brief
 *    Retrieve the comment of a parameter.
 *
 * @param paf   A PAF object.
 * @param name  Parameter name to look up.
 *
 * @return The function returns the parameter's comment string, or @c NULL
 *   if no comment is present. If a parameter @em name does not exist
 *   @c NULL is returned and the error variable @b pilErrno is set to
 *   @c PIL_ENODATA.
 *
 * The function looks for the first occurrence of a parameter with the name
 * @em name in the record list of the PAF object @em paf and returns the
 * comment string of this parameter, or @c NULL if no comment string is
 * present.
 */

const char *
pilPAFGetComment(const PilPAF *paf, const char *name)
{

    PilPAFRecord *record;
    PilListNode *node;

    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    node = pilListLookup(paf->records, name, _pilPAFRecordCompare);
    if (!node) {
        pilErrno = PIL_ENODATA;
        return NULL;
    }

    record = pilListNodeGet(node);

    return record->comment;

}


/**
 * @brief
 *   Assign a boolean value to a parameter.
 *
 * @param paf    A PAF object.
 * @param name   Parameter name.
 * @param value  Parameter value.
 *
 * @return The function returns @c EXIT_SUCCESS if the parameter's value
 *   was successfully changed, or @c EXIT_FAILURE otherwise.
 *
 * The function searches the record list of the PAF object for the first
 * occurrence of a parameter with the name @em name. If such a parameter
 * exists its value is set to @em value.
 */

int
pilPAFSetValueBool(PilPAF *paf, const char *name, int value)
{

    register int status;


    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    status = _pilPAFSetValue(paf->records, name, PAF_TYPE_BOOL, &value);

    return status == 0 ? EXIT_SUCCESS : EXIT_FAILURE;

}


/**
 * @brief
 *   Assign an integer value to a parameter.
 *
 * @param paf    A PAF object.
 * @param name   Parameter name.
 * @param value  Parameter value.
 *
 * @return The function returns @c EXIT_SUCCESS if the parameter's value
 *   was successfully changed, or @c EXIT_FAILURE otherwise.
 *
 * The function searches the record list of the PAF object for the first
 * occurrence of a parameter with the name @em name. If such a parameter
 * exists its value is set to @em value.
 */

int
pilPAFSetValueInt(PilPAF *paf, const char *name, int value)
{

    register int status;


    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    status = _pilPAFSetValue(paf->records, name, PAF_TYPE_INT, &value);

    return status == 0 ? EXIT_SUCCESS : EXIT_FAILURE;

}


/**
 * @brief
 *   Assign a double value to a parameter.
 *
 * @param paf    A PAF object.
 * @param name   Parameter name.
 * @param value  Parameter value.
 *
 * @return The function returns @c EXIT_SUCCESS if the parameter's value
 *   was successfully changed, or @c EXIT_FAILURE otherwise.
 *
 * The function searches the record list of the PAF object for the first
 * occurrence of a parameter with the name @em name. If such a parameter
 * exists its value is set to @em value.
 */

int
pilPAFSetValueDouble(PilPAF *paf, const char *name, double value)
{

    register int status;


    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    status = _pilPAFSetValue(paf->records, name, PAF_TYPE_DOUBLE, &value);

    return status == 0 ? EXIT_SUCCESS : EXIT_FAILURE;

}


/**
 * @brief
 *   Assign a string value to a parameter.
 *
 * @param paf    A PAF object.
 * @param name   Parameter name.
 * @param value  Parameter value.
 *
 * @return The function returns @c EXIT_SUCCESS if the parameter's value
 *   was successfully changed, or @c EXIT_FAILURE otherwise.
 *
 * The function searches the record list of the PAF object for the first
 * occurrence of a parameter with the name @em name. If such a parameter
 * exists its value is set to @em value.
 */

int
pilPAFSetValueString(PilPAF *paf, const char *name, const char *value)
{

    register int status;


    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    status = _pilPAFSetValue(paf->records, name, PAF_TYPE_STRING, value);

    return status == 0 ? EXIT_SUCCESS : EXIT_FAILURE;

}


/**
 * @brief
 *   Set the comment string of a parameter.
 *
 * @param paf      A PAF object.
 * @param name     Parameter name.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the parameter's comment
 *   was successfully set, or @c EXIT_FAILURE otherwise.
 *
 * The function searches the record list of the PAF object for the first
 * occurrence of a parameter with the name @em name. If such a parameter
 * exists its comment string is replaced by @em comment.
 */

int
pilPAFSetComment(PilPAF *paf, const char *name, const char *comment)
{

    PilPAFRecord *record;
    PilListNode *node;

    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    node = pilListLookup(paf->records, name, _pilPAFRecordCompare);
    if (!node)
        return EXIT_FAILURE;

    record = pilListNodeGet(node);
    if (record->comment) {
        size_t sz = strlen(comment);
        
        if (sz != strlen(record->comment)) {
            record->comment = (char *)pil_realloc(record->comment, 
                                                  (sz + 1) * sizeof(char));
            if (!record->comment)
                return EXIT_FAILURE;
        }

        memcpy(record->comment, comment, sz + 1);
    }
    else
        record->comment = pil_strdup(comment);

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a boolean value into a PAF object at a given position.
 *
 * @param paf      PAF object into which the record is inserted.
 * @param key      Position where the value is inserted.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully inserted, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it, at the position given by the parameter with the name
 * @em key, into the record list of the PAF object @em paf. If @em paf
 * contains multiple duplicates of @em key, the new value is inserted at
 * the position of the first occurence of @em key. If @em key does not
 * exist, the function fails and the value is not inserted.
 */

 int
pilPAFInsertBool(PilPAF *paf, const char *key, const char *name, int value,
                 const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFInsert(paf->records, key, name, PAF_TYPE_BOOL, &value,
                      comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a integer value into a PAF object at a given position.
 *
 * @param paf      PAF object into which the record is inserted.
 * @param key      Position where the value is inserted.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully inserted, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it, at the position given by the parameter with the name
 * @em key, into the record list of the PAF object @em paf. If @em paf
 * contains multiple duplicates of @em key, the new value is inserted at
 * the position of the first occurence of @em key. If @em key does not
 * exist, the function fails and the value is not inserted.
 */

 int
pilPAFInsertInt(PilPAF *paf, const char *key, const char *name, int value,
                const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFInsert(paf->records, key, name, PAF_TYPE_INT, &value,
                      comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a double value into a PAF object at a given position.
 *
 * @param paf      PAF object into which the record is inserted.
 * @param key      Position where the value is inserted.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully inserted, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it, at the position given by the parameter with the name
 * @em key, into the record list of the PAF object @em paf. If @em paf
 * contains multiple duplicates of @em key, the new value is inserted at
 * the position of the first occurence of @em key. If @em key does not
 * exist, the function fails and the value is not inserted.
 */

 int
pilPAFInsertDouble(PilPAF *paf, const char *key, const char *name,
                   double value, const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFInsert(paf->records, key, name, PAF_TYPE_DOUBLE, &value,
                      comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a string value into a PAF object at a given position.
 *
 * @param paf      PAF object into which the record is inserted.
 * @param key      Position where the value is inserted.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully inserted, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it, at the position given by the parameter with the name
 * @em key, into the record list of the PAF object @em paf. If @em paf
 * contains multiple duplicates of @em key, the new value is inserted at
 * the position of the first occurence of @em key. If @em key does not
 * exist, the function fails and the value is not inserted.
 */

 int
pilPAFInsertString(PilPAF *paf, const char *key, const char *name,
                   const char *value, const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFInsert(paf->records, key, name, PAF_TYPE_STRING, value,
                      comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a boolean value into a PAF object after a given position.
 *
 * @param paf      PAF object into which the record is inserted.
 * @param key      Position where the value is inserted.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully inserted, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it, after the position given by the parameter with the name
 * @em key, into the record list of the PAF object @em paf. If @em paf
 * contains multiple duplicates of @em key, the new value is inserted after
 * the position of the first occurence of @em key. If @em key does not
 * exist, the function fails and the value is not inserted.
 */

 int
pilPAFInsertAfterBool(PilPAF *paf, const char *key, const char *name,
                      int value, const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFInsertAfter(paf->records, key, name, PAF_TYPE_BOOL, &value,
                           comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a integer value into a PAF object after a given position.
 *
 * @param paf      PAF object into which the record is inserted.
 * @param key      Position where the value is inserted.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully inserted, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it, after the position given by the parameter with the name
 * @em key, into the record list of the PAF object @em paf. If @em paf
 * contains multiple duplicates of @em key, the new value is inserted after
 * the position of the first occurence of @em key. If @em key does not
 * exist, the function fails and the value is not inserted.
 */

 int
pilPAFInsertAfterInt(PilPAF *paf, const char *key, const char *name,
                     int value, const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFInsertAfter(paf->records, key, name, PAF_TYPE_INT, &value,
                      comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a double value into a PAF object after a given position.
 *
 * @param paf      PAF object into which the record is inserted.
 * @param key      Position where the value is inserted.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully inserted, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it, after the position given by the parameter with the name
 * @em key, into the record list of the PAF object @em paf. If @em paf
 * contains multiple duplicates of @em key, the new value is inserted after
 * the position of the first occurence of @em key. If @em key does not
 * exist, the function fails and the value is not inserted.
 */

 int
pilPAFInsertAfterDouble(PilPAF *paf, const char *key, const char *name,
                        double value, const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFInsertAfter(paf->records, key, name, PAF_TYPE_DOUBLE, &value,
                           comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a string value into a PAF object after a given position.
 *
 * @param paf      PAF object into which the record is inserted.
 * @param key      Position where the value is inserted.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully inserted, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it, after the position given by the parameter with the name
 * @em key, into the record list of the PAF object @em paf. If @em paf
 * contains multiple duplicates of @em key, the new value is inserted after
 * the position of the first occurence of @em key. If @em key does not
 * exist, the function fails and the value is not inserted.
 */

 int
pilPAFInsertAfterString(PilPAF *paf, const char *key, const char *name,
                        const char *value, const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFInsertAfter(paf->records, key, name, PAF_TYPE_STRING, value,
                           comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Prepend a boolean value to a PAF object.
 *
 * @param paf      PAF object to which the record is prepended.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully prepended, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and inserts it as first element to the record list of the PAF object
 * @em paf.
 */

 int
pilPAFPrependBool(PilPAF *paf, const char *name, int value,
                  const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFPrepend(paf->records, name, PAF_TYPE_BOOL, &value, comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Prepend a integer value to a PAF object.
 *
 * @param paf      PAF object to which the record is prepended.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully prepended, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and prepends it to the record list of the PAF object @em paf.
 */

 int
pilPAFPrependInt(PilPAF *paf, const char *name, int value, const char *comment)
{

    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFPrepend(paf->records, name, PAF_TYPE_INT, &value, comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Prepend a double value to a PAF object.
 *
 * @param paf      PAF object to which the record is prepended.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully prepended, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and prepends it to the record list of the PAF object @em paf.
 */

 int
pilPAFPrependDouble(PilPAF *paf, const char *name, double value,
                   const char *comment)
{

    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFPrepend(paf->records, name, PAF_TYPE_DOUBLE, &value, comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Prepend a string value to a PAF object.
 *
 * @param paf      PAF object to which the record is prepended.
 * @param name     Parameter name.
 * @param value    Parameter value string.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully prepended, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and prepends it to the record list of the PAF object @em paf.
 */

 int
pilPAFPrependString(PilPAF *paf, const char *name, const char *value,
                   const char *comment)
{

    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFPrepend(paf->records, name, PAF_TYPE_STRING, value, comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Append a boolean value to a PAF object.
 *
 * @param paf      PAF object to which the record is appended.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully appended, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and appends it to the record list of the PAF object @em paf.
 */

 int
pilPAFAppendBool(PilPAF *paf, const char *name, int value, const char *comment)
{
    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFAppend(paf->records, name, PAF_TYPE_BOOL, &value, comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Append a integer value to a PAF object.
 *
 * @param paf      PAF object to which the record is appended.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully appended, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and appends it to the record list of the PAF object @em paf.
 */

 int
pilPAFAppendInt(PilPAF *paf, const char *name, int value, const char *comment)
{

    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFAppend(paf->records, name, PAF_TYPE_INT, &value, comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Append a double value to a PAF object.
 *
 * @param paf      PAF object to which the record is appended.
 * @param name     Parameter name.
 * @param value    Parameter value.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully appended, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and appends it to the record list of the PAF object @em paf.
 */

 int
pilPAFAppendDouble(PilPAF *paf, const char *name, double value,
                   const char *comment)
{

    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFAppend(paf->records, name, PAF_TYPE_DOUBLE, &value, comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Append a string value to a PAF object.
 *
 * @param paf      PAF object to which the record is appended.
 * @param name     Parameter name.
 * @param value    Parameter value string.
 * @param comment  Parameter description.
 *
 * @return The function returns @c EXIT_SUCCESS if the record was
 *   successfully appended, or @c EXIT_FAILURE otherwise.
 *
 * The function creates a new PAF record for the given name, value and comment
 * and appends it to the record list of the PAF object @em paf.
 */

 int
pilPAFAppendString(PilPAF *paf, const char *name, const char *value,
                   const char *comment)
{

    assert(paf != NULL);
    assert(name != NULL);
    
    if (!pilPAFIsValidName(name) && name[0] != '#' && name [0] != '\0')
        return EXIT_FAILURE;

    assert(paf->records != NULL);

    if (_pilPAFAppend(paf->records, name, PAF_TYPE_STRING, value, comment))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Erase the PAF record with the given name.
 *
 * @param paf   PAF object.
 * @param name  Parameter name of the parameter to remove.
 *
 * @return Nothing.
 *
 * The function finds the first occurrence of the parameter with the name
 * @em name in the list of PAF records and erases it from the list, destroying
 * the record data. If no such parameter is found nothing is done to the 
 * record list. If there are multiple duplicates of @em name present in the
 * record list only the first one is erased.
 */

 void
pilPAFErase(PilPAF *paf, const char *name)
{

    PilListNode *node;

    assert(paf != NULL);
    assert(paf->records != NULL);
    assert(name != NULL);

    node = pilListLookup(paf->records, name, _pilPAFRecordCompare);
    pilListErase(paf->records, node, (void (*)(void *))_pilPAFRecordDestroy);

    return;

}


/**
 * @brief
 *   Erase all PAF records.
 *
 * @param paf   PAF object.
 *
 * @return Nothing.
 *
 * The function removes all records of the PAF object @em paf from the record
 * list and destroys them, i.e. deallocates each removed record. After calling
 * this function the record list is empty. The PAF header is not modified in
 * any way by this function.
 */

 void
pilPAFClear(PilPAF *paf)
{

    PilListNode *node;

    if (!paf)
        return;

    assert(paf->records != NULL);

    node = pilListBegin(paf->records);
    while (node) {
        PilListNode *n = node;
        
        node = pilListNext(paf->records, node);
        pilListErase(paf->records, n, (void (*)(void *))_pilPAFRecordDestroy);
    }

    assert(pilListIsEmpty(paf->records));

    return;

}


/**
 * @brief
 *   Write a PAF object to a disk file.
 *
 * @param paf  PAF object to save.
 *
 * @return The function returns @c EXIT_SUCCESS if the object was
 *   successfully written to the file, otherwise @c EXIT_FAILURE is
 *   returned.
 *
 * The function converts the PAF object into an ASCII file on disk, using
 * the name used to create the PAF object as a file name for the disk file.
 */

int
pilPAFWrite(PilPAF *paf)
{

    const char *record;
    const char *user, *timestamp;

    PilListNode *node;

    FILE *stream;

#if defined HAVE_GETUID && defined HAVE_GETPWUID
    struct passwd *pw;
#endif


    if (!paf)
        return EXIT_FAILURE;

    assert(paf->header != NULL);


    /*
     * Finalize the header. If necessary the creator and the creation time
     * are set and the timestamp for the last modification.
     */

    /* Get user id */

#if defined HAVE_GETUID && defined HAVE_GETPWUID
    pw = getpwuid(getuid());

    if (!pw)
        return EXIT_FAILURE;

    user = pw->pw_name;
#else
    user = getenv("USER");
    user = user == NULL ? getenv("LOGNAME") : user;

    if (!user)
        return EXIT_FAILURE;
#endif

    /* Get timestamp */

    timestamp = pilDateGetISO8601();
    if (!timestamp)
        return EXIT_FAILURE;


    if (!_pilPAFGetValueString(paf->header, PAF_CRTE_NAME))
        _pilPAFSetValue(paf->header, PAF_CRTE_NAME, PAF_TYPE_STRING, user);

    if (!_pilPAFGetValueString(paf->header, PAF_CRTE_TIME))
        _pilPAFSetValue(paf->header, PAF_CRTE_TIME, PAF_TYPE_STRING,
                        timestamp);

    _pilPAFSetValue(paf->header, PAF_LCHG_NAME, PAF_TYPE_STRING, user);
    _pilPAFSetValue(paf->header, PAF_LCHG_TIME, PAF_TYPE_STRING, timestamp);
    

    /*
     * Create output file
     */

    stream = fopen(_pilPAFGetValueString(paf->header, PAF_NAME), "wb");
    if (!stream)
        return EXIT_FAILURE;


    /*
     * Write the PAF header
     */

    node = pilListBegin(paf->header);


    /*
     * Fail on empty header. This should never happen!
     */

    if (!node) {
        fclose(stream);
        return EXIT_FAILURE;
    }
        
    while (node) {
        record = _pilPAFFormatRecord(pilListNodeGet(node));
        if (!record) {
            fclose(stream);
            return EXIT_FAILURE;
        }
            
        fprintf(stream, "%s\n", record);
        node = pilListNext(paf->header, node);
    }


    /*
     * Write the PAF records
     */

    node = pilListBegin(paf->records);

    if (node) {
        char buffer[PAF_RECORD_MAX];

        buffer[0] = '#';
        memset(&buffer[1], '-', 78);
        buffer[79] = '\0';
        fprintf(stream, "%s\n", buffer);
    }

    while (node) {
        record = _pilPAFFormatRecord(pilListNodeGet(node));
        if (!record) {
            fclose(stream);
            return EXIT_FAILURE;
        }

        fprintf(stream, "%s\n", record);
        node = pilListNext(paf->records, node);
    }

    fclose(stream);

    return EXIT_SUCCESS;
    
}
/**@}*/
