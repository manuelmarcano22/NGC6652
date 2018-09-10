/* $Id: pilfits.c,v 1.5 2013-08-07 16:26:39 cgarcia Exp $
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
 * $Author: cgarcia $
 * $Date: 2013-08-07 16:26:39 $
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <regex.h>

#include "fitsio.h"
/* #include "fitsio2.h" */

#include "md5.h"

#include "pilmemory.h"
#include "pilstrutils.h"
#include "pilmessages.h"
#include "pilfits.h"


/*
 * Size of an MD5 hash in bytes (32 bytes are 128 bits) 
 */

#define MD5HASHSZ  32


/*
 * Here is a number of definitions strictly related to the FITS format. 
 */

#define PIL_FITS_BLOCK_SIZE (2880)
#define PIL_FITS_NCARDS     ((PIL_FITS_BLOCK_SIZE) / ((PIL_FITS_CARD_MAX) - 1))


/*
 * The PilFitsFile object type
 */

struct _PIL_FITS_FILE_ {
    fitsfile *fptr;
};


/**
 * @defgroup pilFits pilFits
 *
 * The module @b pilFits currently provides just functions to create,
 * remove, and maintain FITS header keywords.
 */

/**@{*/

/*
 * @brief
 *   Get type of an existing keyword.
 *
 * @param fitsFile     FITS file descriptor.
 *
 * @return Keyword type: 1 = character string
 *                       2 = logical
 *                       3 = integer
 *                       4 = float
 *                       5 = complex
 *                    or 0 = keyword not found
 *
 * The function returns the type of an existing keyword, otherwise
 * returns zero.
 */

 static int
pilFitsKeyType(PilFitsFile *fitsFile, const char *keyName)
{

    int  status = 0;

    char card[PIL_FITS_CARD_MAX];
    char comment[PIL_FITS_CARD_MAX];
    char value[PIL_FITS_VALUE_MAX];
    char type[2];


    if (ffgcrd(fitsFile->fptr, (char *)keyName, card, &status))
        return 0;

    ffpsvc(card, value, comment, &status);
    if (status || value[0] == '\0')
        return 0;

    ffdtyp(value, type, &status);

    switch (type[0]) {
        case 'C': return 1;
        case 'L': return 2;
        case 'I': return 3;
        case 'F': return 4;
        case 'X': return 5;
        default : return 0;
    }

}


/**
 * @brief
 *   Destroy a FITS file object.
 *
 * @param fitsFile     FITS file descriptor.
 *
 * @return Nothing.
 *
 * The function destroys the given FITS file object, i.e. the opened
 * stream the object refers to is first closed and then the object is
 * deallocated.
 */

void
deletePilFitsFile(PilFitsFile *fitsFile)
{

    if (fitsFile) {
        if (fitsFile->fptr) {
            int status = 0;

            ffclos(fitsFile->fptr, &status);
        }

        pil_free(fitsFile);
    }

    return;

}


/**
 * @brief
 *   Create a new FITS file object.
 *
 * @param filename  File name.
 * @param mode      IO mode.
 *
 * @return Newly created FITS file object, or NULL in case of failure.
 *
 * The function creates a FITS file object from the named file @em filename
 * for the given I/O mode @em mode. If the function returns successfully
 * the file @em filename has been opened for the given mode.
 *
 * @see deletePilFitsFile()
 */

PilFitsFile *
newPilFitsFile(const char *filename, PilFitsIOMode mode)
{

  PilFitsFile *fitsFile = (PilFitsFile *)pil_malloc(sizeof *fitsFile);


  if (fitsFile) {
      int io_mode;
      int status = 0;

      fitsFile->fptr = NULL;

      switch (mode) {
          case PIL_FITS_READ:
              io_mode = READONLY;
              break;

          case PIL_FITS_WRITE:
              io_mode = READWRITE;
              break;

          case PIL_FITS_READWRITE:
              io_mode = READWRITE;
              break;

          default:
              deletePilFitsFile(fitsFile);
              return NULL;
              break;
      }

      if (ffopen(&fitsFile->fptr, filename, io_mode, &status)) {
          deletePilFitsFile(fitsFile);
          return NULL;
      }

  }
   
  return fitsFile;

}


/**
 * @brief
 *   Get number of header+data sections in FITS file.
 *
 * @param fitsFile     FITS file descriptor.
 *
 * @return Number of header+data units in FITS file, or 0 in case of failure.
 *
 * The function returns the number of header+data sections contained
 * in the FITS file, including the primary array.
 */

int
pilFitsHdrCount(PilFitsFile *fitsFile)
{

    int status = 0;
    int nHdu   = 0;

    if (fitsFile) 
        if (!ffthdu(fitsFile->fptr, &nHdu, &status))
            return nHdu;

    return 0;

}


/**  
 * @brief
 *   Move to a given header+data section in FITS file.
 *
 * @param fitsFile     FITS file descriptor.
 * @param section      Sequence number of header+data section (first = 0).
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * The function is used to access a given header+data section in a FITS file.
 * Sections are counted from 0, that corresponds to the so-called primary 
 * array. 
 */

int pilFitsHdrGoto(PilFitsFile *fitsFile, int section)
{

    int status = 0;


    section++;

    if (fitsFile) 
        if (!ffmahd(fitsFile->fptr, section, NULL, &status))
            return EXIT_SUCCESS;

    return EXIT_FAILURE;

}


/**
 * @brief
 *   Read a keyword value from header as an integer.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Returned value.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function is used to read the value of a given integer keyword.
 * Data type conversion is performed for numerical values if the
 * keyword is not of type integer.
 */

int pilFitsHdrReadInt(PilFitsFile *fitsFile, const char *keyName, int *value)
{

    int status = 0;


    if (fitsFile)
        if (!ffgky(fitsFile->fptr, TINT, (char *)keyName, value,
                   NULL, &status))
            return EXIT_SUCCESS;

    return EXIT_FAILURE;

}


/**
 * @brief
 *   Read a keyword value from header as a float.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Returned value.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function is used to read the value of a given float keyword.
 * Data type conversion is performed for numerical values if the
 * keyword is not of type float.
 */

int pilFitsHdrReadFloat(PilFitsFile *fitsFile, 
                        const char *keyName, float *value)
{

    int status = 0;


    if (fitsFile)
        if (!ffgky(fitsFile->fptr, TFLOAT, (char *)keyName, value,
                   NULL, &status))
            return EXIT_SUCCESS;

    return EXIT_FAILURE;

}


/**
 * @brief
 *   Read a keyword value from header as a double.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Returned value.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function is used to read the value of a given double keyword.
 * Data type conversion is performed for numerical values if the
 * keyword is not of type double.
 */

int pilFitsHdrReadDouble(PilFitsFile *fitsFile,
                         const char *keyName, double *value)
{

    int status = 0;


    if (fitsFile)
        if (!ffgky(fitsFile->fptr, TDOUBLE, (char *)keyName, value,
                   NULL, &status))
            return EXIT_SUCCESS;

    return EXIT_FAILURE;

}


/**
 * @brief
 *   Read a keyword value from header as a string.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Returned pointer to string.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function is used to read the value of a given string keyword.
 * Memory is allocated, therefore the returned string should be 
 * deallocated when no longer needed.
 */
 
int pilFitsHdrReadString(PilFitsFile *fitsFile,
                         const char *keyName, char **value)
{

    int status = 0;


    if (fitsFile)
        *value = (char *)pil_malloc(FLEN_VALUE * sizeof(char));

    if (*value)
        if (!ffgky(fitsFile->fptr, TSTRING, (char *)keyName, *value,
                   NULL, &status))
            return EXIT_SUCCESS;

    pil_free(*value);
      
    return EXIT_FAILURE;
  
}


/**
 * @brief
 *   Write or overwrite an integer keyword to a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Keyword value.
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function overwrites an existing keyword value and comment, or 
 * if a keyword with the same name doesn't exist it appends the new 
 * keyword at the end of the header. A NULL pointer can be given instead
 * of the comment string, and in that case the comment will be left
 * blank or untouched.
 *
 * @see pilFitsHdrInsertInt()
 */

int pilFitsHdrWriteInt(PilFitsFile *fitsFile, 
                       const char *keyName, int value, const char *comment)
{

    int status = 0;
    int type;


    if (!fitsFile)
        return EXIT_FAILURE;

    type = pilFitsKeyType(fitsFile, (char *)keyName);

    if (type && type != 3)
        return EXIT_FAILURE;

    if (ffuky(fitsFile->fptr, 
              TINT, (char *)keyName, &value, (char *)comment, &status))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Write or overwrite a logical keyword to a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Keyword value (0 = false, else true)
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function overwrites an existing keyword value and comment, or
 * if a keyword with the same name doesn't exist it appends the new
 * keyword at the end of the header. A NULL pointer can be given instead
 * of the comment string, and in that case the comment will be left
 * blank or untouched.
 * 
 * @see pilFitsHdrInsertLogical()
 */

int pilFitsHdrWriteLogical(PilFitsFile *fitsFile, 
                           const char *keyName, int value, const char *comment)
{
  
    int   status = 0;
    int type;
    char *f[] = {"F", "T"};

    value = value != 0;

    if (!fitsFile)
        return EXIT_FAILURE;

    type = pilFitsKeyType(fitsFile, (char *)keyName);

    if (type && type != 2)
        return EXIT_FAILURE;

    if (ffuky(fitsFile->fptr, 
              TLOGICAL, (char *)keyName, f[value], (char *)comment, &status))
        return EXIT_FAILURE;
  
    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Write or overwrite a float keyword to a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Keyword value.
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function overwrites an existing keyword value and comment, or
 * if a keyword with the same name doesn't exist it appends the new
 * keyword at the end of the header. A NULL pointer can be given instead
 * of the comment string, and in that case the comment will be left
 * blank or untouched.
 *
 * @see pilFitsHdrInsertFloat()
 */

int pilFitsHdrWriteFloat(PilFitsFile *fitsFile,
                         const char *keyName, float value, const char *comment)
{

    int status = 0;
    int type;


    if (!fitsFile)
        return EXIT_FAILURE;

    type = pilFitsKeyType(fitsFile, (char *)keyName);

    if (type && type != 4)
        return EXIT_FAILURE;

    if (ffuky(fitsFile->fptr, 
              TFLOAT, (char *)keyName, &value, (char *)comment, &status))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Write or overwrite a double keyword to a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Keyword value.
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function overwrites an existing keyword value and comment, or
 * if a keyword with the same name doesn't exist it appends the new
 * keyword at the end of the header. A NULL pointer can be given instead
 * of the comment string, and in that case the comment will be left
 * blank or untouched.
 *
 * @see pilFitsHdrInsertDouble()
 */

int pilFitsHdrWriteDouble(PilFitsFile *fitsFile, const char *keyName, 
                          double value, const char *comment)
{

    int   status = 0;
    int type;


    if (!fitsFile)
        return EXIT_FAILURE;

    type = pilFitsKeyType(fitsFile, (char *)keyName);

    if (type && type != 4)
        return EXIT_FAILURE;

    if (ffuky(fitsFile->fptr,
              TDOUBLE, (char *)keyName, &value, (char *)comment, &status))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Write or overwrite a string keyword to a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Keyword value.
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function overwrites an existing keyword value and comment, or
 * if a keyword with the same name doesn't exist it appends the new
 * keyword at the end of the header. A NULL pointer can be given instead
 * of the comment string, and in that case the comment will be left
 * blank or untouched.
 *
 * @see pilFitsHdrInsertString()
 */

int pilFitsHdrWriteString(PilFitsFile *fitsFile, const char *keyName,
                          const char *value, const char *comment)
{

    int   status = 0;
    int type;


    if (!fitsFile)
        return EXIT_FAILURE;

    type = pilFitsKeyType(fitsFile, (char *)keyName);

    if (type && type != 1)
        return EXIT_FAILURE;

    if (ffuky(fitsFile->fptr, TSTRING, (char *)keyName, 
              (char *)value, (char *)comment, &status))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert an integer keyword in a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param before       0 = insert after reference keyword, else insert before.
 * @param refKeyName   Name of reference keyword.
 * @param keyName      Name of keyword to insert.
 * @param value        Keyword value.
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function inserts the given integer keyword before or after
 * the specified reference keyword. Wildcards can be inserted in the
 * reference keyword name: ? will match any single character, # will
 * match any sequence of consecutive digits, and * will match any 
 * string of characters. If the reference keyword is not found, a
 * failure is returned. A NULL pointer can be given instead of the 
 * comment string, and in that case the keyword comment will be left 
 * blank.
 *
 * @see pilFitsHdrWriteInt()
 */

int pilFitsHdrInsertInt(PilFitsFile *fitsFile, int before, 
                        const char *refKeyName, const char *keyName,
                        int value, const char *comment)
{

    int   status = 0;
    int   countKeys = 0;
    int   position = 0;
    char  card[FLEN_CARD];


    if (!fitsFile)
        return EXIT_FAILURE;

    ffgrec(fitsFile->fptr, 0, card, &status);
    ffgnxk(fitsFile->fptr, (char **)&refKeyName, 1, NULL, 0, card, &status);

    if (before) {

        /*
         * Move the current header position to the correct insertion position.
         */

        ffghps(fitsFile->fptr, &countKeys, &position, &status);
        position -= 2;
        ffgrec(fitsFile->fptr, position, card, &status);

    }

    ffikyj(fitsFile->fptr, (char *)keyName, value, (char *)comment, &status);

    if (status)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a float keyword in a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param before       0 = insert after reference keyword, else insert before.
 * @param refKeyName   Name of reference keyword.
 * @param keyName      Name of keyword to insert.
 * @param value        Keyword value.
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function inserts the given float keyword before or after
 * the specified reference keyword. Wildcards can be inserted in the
 * reference keyword name: ? will match any single character, # will
 * match any sequence of consecutive digits, and * will match any
 * string of characters. If the reference keyword is not found, a
 * failure is returned. A NULL pointer can be given instead of the
 * comment string, and in that case the keyword comment will be left
 * blank.
 *
 * @see pilFitsHdrWriteFloat()
 */

int pilFitsHdrInsertFloat(PilFitsFile *fitsFile, int before,
                          const char *refKeyName, const char *keyName,
                          float value, const char *comment)
{

    int   status = 0;
    int   countKeys = 0;
    int   position = 0;
    char  card[FLEN_CARD];

    static int decimals = 6;


    if (!fitsFile)
        return EXIT_FAILURE;

    ffgrec(fitsFile->fptr, 0, card, &status);
    ffgnxk(fitsFile->fptr, (char **)&refKeyName, 1, NULL, 0, card, &status);

    if (before) {

        /*
         * Move the current header position to the correct insertion position.
         */

        ffghps(fitsFile->fptr, &countKeys, &position, &status);
        position -= 2;
        ffgrec(fitsFile->fptr, position, card, &status);

    }

    ffikyd(fitsFile->fptr, (char *)keyName,
           value, decimals, (char *)comment, &status);

    if (status)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a double keyword in a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param before       0 = insert after reference keyword, else insert before.
 * @param refKeyName   Name of reference keyword.
 * @param keyName      Name of keyword to insert.
 * @param value        Keyword value.
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function inserts the given double keyword before or after
 * the specified reference keyword. Wildcards can be inserted in the
 * reference keyword name: ? will match any single character, # will
 * match any sequence of consecutive digits, and * will match any
 * string of characters. If the reference keyword is not found, a
 * failure is returned. A NULL pointer can be given instead of the
 * comment string, and in that case the keyword comment will be left
 * blank.
 *
 * @see pilFitsHdrWriteDouble()
 */

int pilFitsHdrInsertDouble(PilFitsFile *fitsFile, int before,
                           const char *refKeyName, const char *keyName,
                           double value, const char *comment)
{

    int   status = 0;
    int   countKeys = 0;
    int   position = 0;
    char  card[FLEN_CARD];

    static int decimals = 11;


    if (!fitsFile)
        return EXIT_FAILURE;

    ffgrec(fitsFile->fptr, 0, card, &status);
    ffgnxk(fitsFile->fptr, (char **)&refKeyName, 1, NULL, 0, card, &status);

    if (before) {

        /*
         * Move the current header position to the correct insertion position.
         */

        ffghps(fitsFile->fptr, &countKeys, &position, &status);
        position -= 2;
        ffgrec(fitsFile->fptr, position, card, &status);

    }

    ffikyd(fitsFile->fptr, (char *)keyName,
           value, decimals, (char *)comment, &status);

    if (status)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a character keyword in a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param before       0 = insert after reference keyword, else insert before.
 * @param refKeyName   Name of reference keyword.
 * @param keyName      Name of keyword to insert.
 * @param value        Keyword value.
 * @param comment      Keyword comment.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function inserts the given character keyword before or after
 * the specified reference keyword. Wildcards can be inserted in the
 * reference keyword name: ? will match any single character, # will
 * match any sequence of consecutive digits, and * will match any
 * string of characters. If the reference keyword is not found, a
 * failure is returned. A NULL pointer can be given instead of the
 * comment string, and in that case the keyword comment will be left
 * blank.
 *
 * @see pilFitsHdrWriteString()
 */

int pilFitsHdrInsertString(PilFitsFile *fitsFile, int before,
                           const char *refKeyName, const char *keyName,
                           const char *string, const char *comment)
{

    int   status = 0;
    int   countKeys = 0;
    int   position = 0;
    char  card[FLEN_CARD];


    if (!fitsFile)
        return EXIT_FAILURE;

    ffgrec(fitsFile->fptr, 0, card, &status);
    ffgnxk(fitsFile->fptr, (char **)&refKeyName, 1, NULL, 0, card, &status);

    if (before) {

        /*
         * Move the current header position to the correct insertion position.
         */

        ffghps(fitsFile->fptr, &countKeys, &position, &status);
        position -= 2;
        ffgrec(fitsFile->fptr, position, card, &status);

    }

    ffikys(fitsFile->fptr, (char *)keyName, 
           (char *)string, (char *)comment, &status);

    if (status)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Write or overwrite the DATE keyword.
 *
 * @param fitsFile     FITS file descriptor.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function overwrites the existing DATE keyword with the current
 * system date, or it appends it to the header if a DATE keyword doesn't 
 * already exist. The format "yyyy-mm-ddThh:mm:ss" is used.
 *
 * @see pilFitsHdrInsertDate()
 */

int pilFitsHdrWriteDate(PilFitsFile *fitsFile)
{

    int status = 0;

    if (fitsFile)
        if (!ffpdat(fitsFile->fptr, &status))
            return EXIT_SUCCESS;

    return EXIT_FAILURE;

}


/**
 * @brief
 *   Inserts the DATE keyword at a header position.
 *
 * @param fitsFile     FITS file descriptor.
 * @param before       0 = insert after reference keyword, else insert before.
 * @param refKeyName   Name of reference keyword.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * Inserts a DATE keyword with the current system date at the specified 
 * header position. If a DATE keyword already exists, it is removed.
 *
 * @see pilFitsHdrWriteDate()
 */

int pilFitsHdrInsertDate(PilFitsFile *fitsFile, 
                         int before, const char *refKeyName)
{

    int   status = 0;
    int   countKeys = 0;
    int   position = 0;
    char  datecard[FLEN_CARD];
    char  card[FLEN_CARD];


    if (!fitsFile)
        return EXIT_FAILURE;

    ffpdat(fitsFile->fptr, &status);
    ffgcrd(fitsFile->fptr, "DATE", datecard, &status);
    ffdkey(fitsFile->fptr, "DATE", &status);
  
    ffgrec(fitsFile->fptr, 0, card, &status);
    ffgnxk(fitsFile->fptr, (char **)&refKeyName, 1, NULL, 0, card, &status);
    ffghps(fitsFile->fptr, &countKeys, &position, &status);

    if (before) {

        /*
         * Move the current header position to the correct insertion position.
         */

        position--;
        ffgrec(fitsFile->fptr, position, card, &status);

    }

    ffirec(fitsFile->fptr, position, datecard, &status);

    if (status)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Read a record from FITS header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param value        Returned pointer to string.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function is used to read a given FITS header record, identifyed 
 * by its keyword name. Memory is allocated, therefore the returned string 
 * should be deallocated when no longer needed.
 */

int pilFitsHdrReadCard(PilFitsFile *fitsFile, const char *keyName, char **card)
{

    int status = 0;


    if (!fitsFile)
        return EXIT_FAILURE;

    *card = (char *)pil_malloc(FLEN_CARD * sizeof(char));

    if (*card)
        if (!(ffgcrd(fitsFile->fptr, (char *)keyName, *card, &status)))
            return EXIT_SUCCESS;

    return EXIT_FAILURE;

}


/**
 * @brief
 *   Write or overwrite a record to a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Keyword name.
 * @param card         Record to write or overwrite.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function overwrites an existing record in the FITS header, or
 * if a record with a keyword with the same name doesn't exist it appends 
 * the new record at the end of the header. 
 *
 * @see pilFitsHdrReadCard(), pilFitsHdrInsertCard()
 */

int pilFitsHdrWriteCard(PilFitsFile *fitsFile, 
                        const char *keyName, const char *card)
{

    int status = 0;

    if (!fitsFile)
        return EXIT_FAILURE;

    if (ffucrd(fitsFile->fptr, (char *)keyName, (char *)card, &status))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Insert a record in a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param before       0 = insert after reference keyword, else insert before.
 * @param refKeyName   Name of reference keyword.
 * @param card         Record to insert.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function inserts the given record before or after the specified 
 * reference keyword. Wildcards can be inserted in the reference keyword 
 * name: ? will match any single character, # will match any sequence of 
 * consecutive digits, and * will match any string of characters. If the 
 * reference keyword is not found, a failure is returned.
 *
 * @see pilFitsHdrReadRecord()
 */

int pilFitsHdrInsertCard(PilFitsFile *fitsFile, int before, 
                         const char *refKeyName, const char *card)
{

    int  status = 0;
    int  countKeys = 0;
    int  position = 0;
    char dummy[FLEN_CARD];


    if (!fitsFile)
        return EXIT_FAILURE;

    ffgrec(fitsFile->fptr, 0, dummy, &status);
    ffgnxk(fitsFile->fptr, (char **)&refKeyName, 1, NULL, 0, dummy, &status);
    ffghps(fitsFile->fptr, &countKeys, &position, &status);

    if (before) {

        /*
         * Move the current header position to the correct insertion position.
         */

        position--;
        ffgrec(fitsFile->fptr, position, dummy, &status);

    }

    ffirec(fitsFile->fptr, position, (char *)card, &status);

    if (status)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Write a comment to a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param comment      Comment to write.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function appends comments to the FITS header. The comment
 * string will be split over multiple comment lines if it is longer
 * than 70 characters.
 */

int pilFitsHdrWriteComment(PilFitsFile *fitsFile, const char *comment)
{

    int status = 0;


    if (fitsFile)
        if (comment)
            if (!(ffpcom(fitsFile->fptr, (char *)comment, &status)))
                return EXIT_SUCCESS;

    return EXIT_FAILURE;
  
}


/**
 * @brief
 *   Write history to a FITS file header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param history      History to write.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function appends history comments to the FITS header. The history
 * string will be split over multiple lines if it is longer than 70 
 * characters.
 */

int pilFitsHdrWriteHistory(PilFitsFile *fitsFile, const char *history)
{

    int status = 0;


    if (fitsFile)
        if (history)
            if (!(ffphis(fitsFile->fptr, (char *)history, &status)))
                return EXIT_SUCCESS;

    return EXIT_FAILURE;

}


/**
 * @brief
 *   Delete a keyword from header.
 *
 * @param fitsFile     FITS file descriptor.
 * @param keyName      Name of keyword to delete.
 *
 * @return EXIT_SUCCESS of EXIT_FAILURE;
 *
 * This function deletes an existing keyword from the FITS header.
 * If a keyword with the specified name is not found, a failure is
 * returned.
 */


int pilFitsHdrDelete(PilFitsFile *fitsFile, const char *keyName)
{

    int status = 0;


    if (fitsFile)
        if (!ffdkey(fitsFile->fptr, (char *)keyName, &status))
            return EXIT_FAILURE;

    return EXIT_SUCCESS;

}

int pilFitsHdrDeleteKeys(const char *filename, const char *keyname, int ext) 
{

    fitsfile *file;
    int status = 0;
    const char *fid = "pilFitsHdrDeleteKeys";

    assert(filename != NULL);
    assert(keyname != NULL);

    if (ffopen(&file, filename, READWRITE, &status)) {
        pilMsgError(fid, "Cannot open file `%s'", filename);
        return EXIT_FAILURE;
    }

    if (fits_movabs_hdu(file, ext + 1, NULL, &status)) {
        status = 0;
        pilMsgError(fid, "Cannot access extension header %d", ext);
        fits_close_file(file, &status);
        return EXIT_FAILURE;
    }

    fits_delete_key(file, keyname, &status);
    fits_close_file(file, &status);

    return EXIT_SUCCESS;
}


/**
 * @brief
 *   Copy a set of keywords from one extension header to another.
 *
 * @param filename  Name of the FITS file.
 * @param target    Target extension number.
 * @param hint      Position hint keyword insertion. (Not yet implemented!)
 * @param name      Source keyword name pattern
 * @param source    Source extension number.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occured, or
 *   @c EXIT_FAILURE otherwise.
 *
 * The function copies all keyword matching the name pattern @em name
 * of the extension header @em source to the extension header @em target.
 * The pattern @em name may be a regular expression. If either the source
 * or target header are not present the function returns an error. If no
 * keyword matching the name pattern is found in the source header the
 * function exits successfully. The numbering of the extensions is
 * assumend to start from 0, where 0 indicates the file's primary header.
 */

int pilFitsHdrCopy(const char *filename, unsigned int target, const char *hint,
                   const char *name, unsigned int source)
{

    const char *fctid = "pilFitsHdrCopy";

    char **key_cache;
    char key[FLEN_CARD];

    register int i, n;

    int sz;
    int re_status;
    int status = 0;
    int key_count = 0;
    (void)hint;

    regex_t re;

    fitsfile *file;


    assert(filename != NULL);
    assert(name != NULL);

    re_status = regcomp(&re, name, REG_EXTENDED | REG_NOSUB);
    if (re_status) {
        return EXIT_FAILURE;
    }

    if (ffopen(&file, filename, READWRITE, &status)) {
        pilMsgError(fctid, "Cannot open file `%s'", filename);
        return EXIT_FAILURE;
    }


    /*
     * Check whether we can access the source and target extension headers
     */

    n = MIN(source, target) + 1;
    for (i = 0; i < 2; i++) {
        if (fits_movabs_hdu(file, n, NULL, &status)) {
            status = 0;
            pilMsgError(fctid, "Cannot access extension header %d", n);
            fits_close_file(file, &status);
            return EXIT_FAILURE;
        }

        n = MAX(source, target) + 1;
    }


    /*
     * Create a cache for the source keywords, traverse the source
     * header and add an entry to the cache for each matching keyword
     * name.
     */

    fits_movabs_hdu(file, source + 1, NULL, &status);
    fits_get_hdrspace(file, &key_count, NULL, &status);

    key_cache = (char **)pil_calloc(key_count + 1, sizeof(char *));

    n = 0;
    for (i = 0; i < key_count; i++) {
        char card[FLEN_CARD];

        fits_read_record(file, i + 1, card, &status);
        fits_get_keyname(card, key, &sz, &status);

        if (regexec(&re, key, 0, NULL, 0) == 0) {
            key_cache[n] = pil_strdup(card);
            ++n;
        }
    }
    regfree(&re);
    

    /*
     * Write the contents of the cache to the target header. We
     * stop at the first empty slot in the cache.
     */

    fits_movabs_hdu(file, target + 1, NULL, &status);

    i = 0;
    while (key_cache[i]) {

        fits_get_keyname(key_cache[i], key, &sz, &status);

        /*
         * Try to modify the card with the keyword key. If this
         * fails it does not exist and we just append it.
         */

        fits_modify_card(file, key, key_cache[i], &status);
        if (status != 0) {
            status = 0;
            fits_write_record(file, key_cache[i], &status);
            if (status != 0) {
                pilMsgWarning(fctid, "Cannot write keyword `%s' to header %d",
                              key, target);
            }
        }
        ++i;
    }

    status = 0;
    fits_close_file(file, &status);

    /*
     * Destroy the cache
     */

    i = 0;
    while (key_cache[i]) {
        pil_free(key_cache[i]);
        ++i;
    }

    pil_free(key_cache);

    return EXIT_SUCCESS;

}

/**
 * @brief  
 *   Compute the MD5 hash of data zones in a FITS file.
 *
 * @param   filename  Name of the FITS file to examine.
 *
 * @return  Statically allocated character string, or NULL.
 *
 * This function expects the name of a FITS file.
 * It will compute the MD5 hash on all data blocks in the main data section
 * and possibly extensions (including zero-padding blocks if necessary) and
 * return it as a string suitable for inclusion into a FITS keyword.
 *
 * The returned string is statically allocated inside this function,
 * so do not free it or modify it. This function returns NULL in case
 * of error.
 */

char *
pilFitsMD5Signature(const char *filename)
{

    const char        modName[] = "pilFitsMD5Signature";

    static char       datamd5[MD5HASHSZ+1];
    struct MD5Context ctx;
    unsigned char     digest[16];
    FILE             *in;
    char              buf[PIL_FITS_BLOCK_SIZE];
    char             *buf_c;
    int               i;
    int               in_header;
    int               check_fits;


    /* 
     * Check entries 
     */

    if (filename == NULL) 
        return NULL;

    if ((in = fopen(filename, "r")) == NULL) {
        pilMsgDebug(modName, "Cannot open file %s", filename);
        return NULL;
    }

    /* 
     * Initialize all variables 
     */

    MD5Init(&ctx);

    in_header = 1;
    check_fits = 0;

    /* 
     * Loop over input file 
     */

    while (fread(buf, 1, PIL_FITS_BLOCK_SIZE, in) == PIL_FITS_BLOCK_SIZE) {
        if (check_fits == 0) {
            check_fits = 1;

            /* 
             * First time in the loop: check the file is FITS.
             * Examine first characters in block.
             */

            if (buf[0] != 'S' ||
                buf[1] != 'I' ||
                buf[2] != 'M' ||
                buf[3] != 'P' ||
                buf[4] != 'L' ||
                buf[5] != 'E' ||
                buf[6] != ' ' ||
                buf[7] != ' ' ||
                buf[8] != '='   ) {
                pilMsgDebug(modName, "File [%s] is not FITS", filename);
                fclose(in);
                return NULL;
            }
        }
        if (in_header) {
            buf_c = buf;
            for (i = 0; i < PIL_FITS_NCARDS; i++) {
                if (buf_c[0] == 'E' &&
                    buf_c[1] == 'N' &&
                    buf_c[2] == 'D' &&
                    buf_c[3] == ' ') {
                    in_header=0;
                    break;
                }
                buf_c += PIL_FITS_CARD_MAX - 1;
            }
        } else {
            /* 
             * If current block is a data block
             * try to locate an extension header 
             */
            if (buf[0] == 'X' &&
                buf[1] == 'T' &&
                buf[2] == 'E' &&
                buf[3] == 'N' &&
                buf[4] == 'S' &&
                buf[5] == 'I' &&
                buf[6] == 'O' &&
                buf[7] == 'N' &&
                buf[8] == '='    ) {
                in_header = 1;
                buf_c = buf;
                for (i = 0; i < PIL_FITS_NCARDS; i++) {
                    /* 
                     * Try to find an END marker in this block 
                     */
                    if (buf_c[0] == 'E' &&
                        buf_c[1] == 'N' &&
                        buf_c[2] == 'D' &&
                        buf_c[3] == ' ') {

                        /* 
                         * Found END marker in same block as XTENSION 
                         */

                        in_header = 0;
                        break;
                    }
                    buf_c += PIL_FITS_CARD_MAX - 1;
                }
            } else {
                MD5Update(&ctx, (unsigned char *)buf, PIL_FITS_BLOCK_SIZE);
            }
        }
    }
    fclose(in);
    if (check_fits == 0) {

        /* 
         * Never went through the read loop: file is not FITS 
         */

        pilMsgDebug(modName, "file [%s] is not FITS", filename);
        return NULL;
    }
    /* 
     * Got to the end of file: summarize 
     */

    MD5Final(digest, &ctx);

    /* 
     * Write digest into a string 
     */

    sprintf(datamd5,
            "%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",
            digest[ 0],
            digest[ 1],
            digest[ 2],
            digest[ 3],
            digest[ 4],
            digest[ 5],
            digest[ 6],
            digest[ 7],
            digest[ 8],
            digest[ 9],
            digest[10],
            digest[11],
            digest[12],
            digest[13],
            digest[14],
            digest[15]);

    return datamd5;

}
/**@}*/
