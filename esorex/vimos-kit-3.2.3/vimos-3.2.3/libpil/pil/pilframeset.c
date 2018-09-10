/* $Id: pilframeset.c,v 1.1.1.1 2008-10-21 09:10:13 cizzo Exp $
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
#  include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

#include "pilmemory.h"
#include "pilstrutils.h"
#include "pilutils.h"
#include "pilframeset.h"


/** 
 * @defgroup pilSetOfFrames pilSetOfFrames
 *
 * The module @b pilSetOfFrames provides functions to create, maintain
 * and destroy a collection of frames.
 */

/**@{*/

static int catgcmp(const void *category1, const void *category2)
{

    return strcmp(category1, category2);

}


/** 
 * @brief
 *   Creates a new set of frames.
 *
 * @return Pointer to the newly allocated set of frames.
 *
 * The function allocates and initializes a new, empty frame collection.
 */

PilSetOfFrames *newPilSetOfFrames(void)
{

    PilSetOfFrames *new = newPilDictionary(PIL_DICT_CAPACITY_MAX, catgcmp);

    if (new)
        pilDictAllowDuplicates(new);

    return new;

}


/**
 * @brief
 *   Destroys an existing set of frames.
 *
 * @param set Pointer to the set of frames that shall be destroyed.
 *
 * @return Nothing.
 *
 * The function removes and deallocates all remaining frames from the frame
 * collection and then deallocates the set itself.
 */

void deletePilSetOfFrames(PilSetOfFrames *set)
{

    PilDictNode *node;

    if (!pilDictIsEmpty(set)) {

        /*
         * For all remaining nodes deallocate the keyword and the data first
         */

        node = pilDictBegin(set);
        while (node) {
            deletePilFrame((PilFrame *)pilDictGetData(node));
            pil_free((void *)pilDictGetKey(node));
            node = pilDictNext(set, node);
        }

        /*
         * Remove and deallocate all the nodes in the set
         */

        pilDictClear(set);
    }

    /*
     * Destroy the set itself
     */

    deletePilDictionary(set);

    return;

}


/**
 * @brief
 *   Add a frame to an existing set of frames.
 *
 * @param set Pointer to an existing set of frames.
 * @param frame Pointer to the frame to be added to the set.
 *
 * @return Returns 1 on success or 0 otherwise.
 *
 * An existing frame is inserted into the set using its frame category as
 * criteria where the frame is inserted.
 */

int pilSofInsert(PilSetOfFrames *set, PilFrame *frame)
{

    return pilDictInsert(set, pil_strdup(pilFrmGetCategory(frame)), frame);

}


/**
 * @brief 
 *   Destructive removal of a single frame from a set of frames.
 *
 * @param set       Pointer to an existing set of frames.
 * @param category  Frame category name used for searching the set of frames.
 *
 * @return Nothing.
 *
 * The function searches the set of frames for the first entry with the
 * given frame category and removes it from the set. After its removal
 * from the frame collection, the frame itself is also deallocated.
 */

void pilSofRemove(PilSetOfFrames *set, const char *category)
{

    PilDictNode *frame_node = pilDictLookup(set, category);

    if (frame_node) {
        deletePilFrame((PilFrame *)pilDictGetData(frame_node));
        pil_free((void *)pilDictGetKey(frame_node));
        pilDictErase(set, frame_node);
    }

    return;

}


/**
 * @brief 
 *   Lookup a frame by its frame category name.
 *
 * @param set       Pointer to an existing set of frames.
 * @param category  Frame category name used for searching the set of frames.
 *
 * @return Pointer to the frame in the set or @c NULL.
 *
 * The function searches the set of frames for an entry with the
 * given frame category and returns a pointer to this frame. In the case
 * of having multiple entries with the same category name always the first
 * entry is returned.
 *
 * @see pilSofLookupNext()
 */

PilFrame *pilSofLookup(PilSetOfFrames *set, const char *category)
{

    PilDictNode *frame_node;

    assert(category != 0x0);

    frame_node = pilDictLookup(set, category);

    if (!frame_node)
        return NULL;

    return (PilFrame *)pilDictGetData(frame_node);
  
}


/**
 * @brief 
 *   Lookup the next frame in a chain by its frame category name.
 *
 * @param set       Pointer to an existing set of frames.
 * @param category  Frame category name used for searching the set of frames.
 *
 * @return Pointer to the frame in the set or @c NULL.
 *
 * The function searches the set of frames for an entry with the
 * given frame category and returns a pointer to this frame. In the case
 * of having multiple entries with the same category name the function
 * behaves like pilSofLookup and returns the pointer to the first frame
 * in the chain, but it also stores the current location within the chain
 * internally. Using @c NULL as category name in subsequent calls to the
 * function returns the pointer to the next frame in the chain. 
 * 
 * If no more frames can be found in the chain or if an error occurs 
 * the function returns @c NULL.
 *
 * @see pilSofLookup()
 */

PilFrame *pilSofLookupNext(PilSetOfFrames *set, const char *category)
{

    static PilDictionary *frame_set = NULL;
    static PilDictNode *frame_node = NULL;
    PilDictNode *previous;

    if (!set)
        return NULL;

    if (category) {
        frame_set = set;
        frame_node = pilDictLookup(set, category);
        if (!frame_node)
            return NULL;
    }
    else {

        /*
         * Check if the current set is the same as the set of the previous
         * call to the function.
         */

        if (frame_set != set) {
            frame_set = NULL;
            return NULL;
        }


        /*
         * Check if the frame node from the previous call is an existing node
         * and it is still a member of the set under consideration.
         */

        if (!frame_node || !pilDictContains(set, frame_node))
            return NULL;

        previous = frame_node;
        if (!(frame_node = pilDictNext(set, frame_node)))
            return NULL;

        if (strcmp(pilDictGetKey(frame_node), pilDictGetKey(previous)))
            return NULL;
    }

    return (PilFrame *)pilDictGetData(frame_node);

}


/**
 * @brief
 *   Counts the number of frames in the set with a given frame category.
 *
 * @param set       Pointer to an existing set of frames.
 * @param category  Frame category name to be searched in the frame 
 *                  collection.
 *
 * @return The number of frames.
 *
 * The function loops through the whole set of frames, comparing the
 * given frame category with the category entries of the set and counting
 * the matches.
 */

PilSetCapacity pilSofFrameCount(PilSetOfFrames *set, const char *category)
{

    PilDictNode *frame_node;
    PilSetCapacity count = 0;


    if ((frame_node = pilDictLookup(set, category))) {
        count++;

        while((frame_node = pilDictNext(set, frame_node))) {
            if (strcmp(category, pilDictGetKey(frame_node)))
                break;

            count++;
        }
    }

    return count;

}

/**
 * @brief
 *   Check whether a set of frames is empty.
 *
 * @param set  The set of frames to check. 
 *
 * @return The function returns 1 if the set is empty, and 0 otherwise.
 *
 * The function checks if @em set contains any frames.
 */

int pilSofIsEmpty(const PilSetOfFrames *set)
{
    
    return pilDictIsEmpty((PilDictionary *)set);

}


/**
 * @brief
 *   Get the first frame in the given set.
 *
 * @param set  Set of frames to query.
 *
 * @return The pointer to the first frame in the set. If the set is empty or
 *   the function returns @c NULL.
 *
 * The function finds and returns the first frame in the input frame
 * collection @em set.
 *
 * @see pilSofNext()
 */

PilFrame *pilSofFirst(PilSetOfFrames *set)
{

    PilDictNode *node = NULL;


    if (set)
        node = pilDictBegin(set);

    return node ? pilDictGetData(node) : NULL;

}


/**
 * @brief
 *   Get the next frame a set.
 *
 * @param set    Set of frames to query for the next frame.
 * @param frame  Predceding frame.
 *
 * @return The pointer to the next frame in a set. If there are no more
 *   frames in the set the function returns @c NULL.
 *
 * The function finds and returns the successor of the frame @em frame in
 * the set of frames @em set. If @em frame is the last frame in @em set, i.e.
 * there is no successor frame, the function returns @c NULL.
 *
 * @see pilSofFirst()
 */

PilFrame *pilSofNext(PilSetOfFrames *set, PilFrame *frame)
{

    PilFrame *next = NULL;


    assert(frame != NULL);

    if (set) {
        PilDictNode *node = pilDictBegin(set);

        while (node) {
            if (pilDictGetData(node) == frame) {
                node = pilDictNext(set, node);

                if (node)
                    next = pilDictGetData(node);

                break;
            }

            node = pilDictNext(set, node);
        }
    }

    return next;

}


/**
 * @brief
 *   Create a set of frames from an ASCII file.
 *
 * @param filename  Name of the input ASCII file.
 * @param set       Pointer to an existing set of frames.
 *
 * @return Pointer to the newly created or updated frame collection. If the
 *   conversion of the ASCII file into the set of frames fails the return
 *   value is @c NULL.
 *
 * The function reads the provided ASCII file and adds the contents of
 * the file to an existing set of frames. The function expects that the
 * input file has the following file format:
 *   @li One file specification per line.
 *   @li Each file specification consists of a file name followed by
 *       the frame category string optionally followed by a frame type
 *       string. The fields are seperated by whitespaces.
 */

PilSetOfFrames *pilSofRead(const char *filename, PilSetOfFrames *set)
{


    if (set) {
        char sofline[PIL_LINE_LENGTH_MAX + 1];
        FILE *fp;


        if (!(fp = fopen(filename, "r")))
            return NULL;

        while (fgets(sofline, PIL_LINE_LENGTH_MAX, fp)) {
            if (!strempty(sofline, NULL)) {
                char path[PIL_LINE_LENGTH_MAX + 1];
                char group[PIL_LINE_LENGTH_MAX + 1];
                char tag[PIL_LINE_LENGTH_MAX + 1];

                int n;


                n = sscanf(sofline, "%s %s %s", path, tag, group);

                if (n >= 2) {
                    PilFrame *frame = newPilFrame(path, tag);


                    pilFrmSetType(frame, PIL_FRAME_TYPE_UNDEF);

                    if (n > 2) {
                        if (!strcmp(group, PIL_FRAME_TYPE_RAW_ID)) {
                            pilFrmSetType(frame, PIL_FRAME_TYPE_RAW);
                        }
                        else {
                            if (!strcmp(group, PIL_FRAME_TYPE_CALIB_ID)) {
                                pilFrmSetType(frame, PIL_FRAME_TYPE_CALIB);
                            }
                            else {
                                if (!strcmp(group,
                                            PIL_FRAME_TYPE_PRODUCT_ID)) {
                                    pilFrmSetType(frame,
                                                  PIL_FRAME_TYPE_PRODUCT);
                                }
                            }
                        }
                    }

                    /*
                     * Insert the created frame in the set.
                     */

                    pilSofInsert(set, frame);

                }
                else {

                    /*
                     * Invalid frame description. Cleanup and return.
                     */

                    fclose(fp);
                    return NULL;

                }
            }
        }

        fclose(fp);

    }

    return set;

}


/**
 * @brief
 *   Write a set of frames to an ASCII file.
 *
 * @param set       Pointer to an existing set of frames.
 * @param filename  Name of the output ASCII file.
 *
 * @return The function returns @c EXIT_SUCCESS if no error occurred,
 *   otherwise the return value is @c EXIT_FAILURE.
 *
 * The function writes the input set of frames to the specified ASCII
 * file. The output file has the following file format:
 *   @li One file specification per line.
 *   @li Each file specification consits of a file name followed by
 *       the frame category string. The fields are seperated by
 *       whitespaces.
 */

int pilSofWrite(PilSetOfFrames *set, const char *filename)
{

    PilDictNode *node;
    PilFrame *frame;

    FILE *fp = fopen(filename, "w");


    if (!fp)
        return EXIT_FAILURE;

    for (node = pilDictBegin(set); node; node = pilDictNext(set, node)) {
        frame = (PilFrame *)pilDictGetData(node);

        if (!frame) {
            fclose(fp);
            return EXIT_FAILURE;
        }
        else {
            if (pilFrmGetName(frame)) {
                fprintf(fp, "%s", pilFrmGetName(frame));

                if (pilFrmGetCategory(frame))
                    fprintf(fp, " %s", pilFrmGetCategory(frame));
                else
                    fprintf(fp, " UNKNOWN");

                switch (pilFrmGetType(frame)) {
                    case PIL_FRAME_TYPE_RAW:
                        fprintf(fp, " %s", PIL_FRAME_TYPE_RAW_ID);
                        break;

                    case PIL_FRAME_TYPE_CALIB:
                        fprintf(fp, " %s", PIL_FRAME_TYPE_CALIB_ID);
                        break;

                    case PIL_FRAME_TYPE_PRODUCT:
                        fprintf(fp, " %s", PIL_FRAME_TYPE_PRODUCT_ID);
                        break;

                    default:
                        break;
                }

                fprintf(fp, "\n");

            }
        }
    }

    fclose(fp);
    return EXIT_SUCCESS;

}


/**
 * @brief
 *   Print the contents of the frame collection to the specified stream.
 *
 * @param stream  Output stream.
 * @param format  Single character indicating the output format.
 * @param set     Pointer to an existing set of frames.
 *
 * @return The function returns the number of processed set entries.
 *
 * The function prints the contents of the frame collection to the
 * output stream. The output stream must be of type #FILE *#.
 * The output format is given by a single character. The available
 * format characters are:
 *   @li @b B: Prints the basic information, namely the file name and
 *             the frame category.
 *   @li @b X: Prints the basic information and all additional frame
 *             attributes in a short format.
 *   @li @b I: Prints the frame number within the set, the set keyword and
 *             all frame information in a pretty format.
 * 
 * @note
 *   This function is only provided for debugging purposes. Its public
 *   usage is not recommended, since it is likely to be changed or even
 *   removed.
 */

int pilSofDump(FILE *stream, const char format, PilSetOfFrames *set)
{

    PilDictNode *node;
    PilDictCapacity count = 0, max = pilDictCapacity(set);
    PilFrame *frame;
    char *category;


    for (node = pilDictBegin(set); node; node = pilDictNext(set, node)) {
        ++count;
        category = (char *)pilDictGetKey(node);
        frame = pilDictGetData(node);

        switch (format) {
            case 'B':
                fprintf(stream, "%s\t%s\n", pilFrmGetName(frame),
                        pilFrmGetCategory(frame));
                break;

            case 'X':
                fprintf(stream, "%s\t%s\n", pilFrmGetName(frame),
                        pilFrmGetCategory(frame));
                fprintf(stream, "type = %d, level = %d, keep = %d, ignore = "
                        "%d\n", pilFrmGetType(frame),
                        pilFrmGetProductLevel(frame), pilFrmGetKeepFlag(frame),
                        pilFrmGetIgnoreFlag(frame));
                break;

            case 'I':
                fprintf(stream, "Frame %ld of %ld:\n", count, max);
                fprintf(stream, "  Keyword:\t%s\n", category);
                fprintf(stream, "  File:\t%s\n", pilFrmGetName(frame));
                fprintf(stream, "  Category:\t%s\n", pilFrmGetCategory(frame));
                fprintf(stream, "  Type:\t\t%d\n", pilFrmGetType(frame));
                fprintf(stream, "  Level:\t%d\n",
                        pilFrmGetProductLevel(frame));
                fprintf(stream, "  Keep:\t\t%d\n", pilFrmGetKeepFlag(frame));
                fprintf(stream, "  Ignore:\t%d\n", pilFrmGetIgnoreFlag(frame));
                break;

            default:
                return count;
                break;
        }
    }

    return count;

}
/**@}*/
