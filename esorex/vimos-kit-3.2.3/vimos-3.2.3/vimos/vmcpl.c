/* $Id: vmcpl.c,v 1.5 2013-08-23 10:22:37 cgarcia Exp $
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
 * $Revision: 1.5 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdlib.h>

#include <cxstring.h>

#include <vmutils.h>
#include <piltranslator.h>
#include <pilrecipe.h>

#include "vmcpl.h"


/* FIXME: This is used to compensate for a deficiency in CPL 1.0. Change
 *        request has been submitted, so this symbol must be removed when
 *        we switch to the next version of CPL (1.0.x or 1.x)!
 */

#ifndef CPL_FRAME_TYPE_PAF
#  define CPL_FRAME_TYPE_PAF  (1 << 4)
#endif


/**
 * @defgroup vmcpl CPL interface
 * 
 * This module provides the necessary adaptor functionality in order
 * to execute data reduction tasks using the CPL recipe interface.
 */

/**@{*/

/*
 * Transform a CPL frame into a PilFrame. This function does not handle
 * product specific frame attributes because this function should only
 * be called for the transformation of raw and calibration frames.
 *
 * This relies on the fact that VIMOS does not have recipe chains, so that
 * product frames from previous data reduction tasks should not appear
 * in a set of frames passed to a recipe. Set of frames are 'always'
 * created from a set of frames file, which provides only the filename,
 * the tag (frame category) and the frame group information.
 */

inline static PilFrame *
_vmCplFrameExport(const cpl_frame *self)
{

    PilFrame *frame = NULL;

    if (self != NULL) {

        const char *name = cpl_frame_get_filename(self);
        const char *tag = cpl_frame_get_tag(self);

        cpl_frame_group group = cpl_frame_get_group(self);


        if (name != NULL && tag != NULL) {

            frame = newPilFrame(name, tag);

            switch (group) {
                case CPL_FRAME_GROUP_NONE:
                    pilFrmSetType(frame, PIL_FRAME_TYPE_UNDEF);
                    break;

                case CPL_FRAME_GROUP_RAW:
                    pilFrmSetType(frame, PIL_FRAME_TYPE_RAW);
                    break;

                case CPL_FRAME_GROUP_CALIB:
                    pilFrmSetType(frame, PIL_FRAME_TYPE_CALIB);
                    break;

                case CPL_FRAME_GROUP_PRODUCT:
                    pilFrmSetType(frame, PIL_FRAME_TYPE_PRODUCT);
                    break;

                default:
                    deletePilFrame(frame);
                    return NULL;
                    break;
            }

        }

    }

    return frame;

}


inline static cpl_frame *
_vmCplFrameImport(const PilFrame *frame)
{

    cpl_frame *self = NULL;


    if (frame != NULL) {

        const char *name = pilFrmGetName(frame);
        const char *tag  = pilFrmGetCategory(frame);

        PilFrameType group = pilFrmGetType(frame);
        PilFrameFormat type = pilFrmGetFormat(frame);
        PilProductLevel level = pilFrmGetProductLevel(frame);


        self = cpl_frame_new();

        if (cpl_frame_set_filename(self, name) != 0) {
            cpl_frame_delete(self);
            return NULL;
        }

        if (cpl_frame_set_tag(self, tag) != 0) {
            cpl_frame_delete(self);
            return NULL;
        }

        switch (group) {
            case PIL_FRAME_TYPE_UNDEF:
                cpl_frame_set_group(self, CPL_FRAME_GROUP_NONE);
                break;

            case PIL_FRAME_TYPE_RAW:
                cpl_frame_set_group(self, CPL_FRAME_GROUP_RAW);
                break;

            case PIL_FRAME_TYPE_CALIB:
                cpl_frame_set_group(self, CPL_FRAME_GROUP_CALIB);
                break;

            case PIL_FRAME_TYPE_PRODUCT:
                cpl_frame_set_group(self, CPL_FRAME_GROUP_PRODUCT);
                break;

            default:
                cpl_frame_delete(self);
                return NULL;
                break;
        }

        switch (type) {
            case PIL_FRAME_FORMAT_UNDEF:
                cpl_frame_set_type(self, CPL_FRAME_TYPE_NONE);
                break;

            case PIL_FRAME_FORMAT_IMAGE:
                cpl_frame_set_type(self, CPL_FRAME_TYPE_IMAGE);
                break;

            case PIL_FRAME_FORMAT_TABLE:
                cpl_frame_set_type(self, CPL_FRAME_TYPE_TABLE);
                break;

            case PIL_FRAME_FORMAT_PAF:
                cpl_frame_set_type(self, CPL_FRAME_TYPE_PAF);
                break;

            default:
                cpl_frame_delete(self);
                return NULL;
                break;
        }

        switch (level) {
            case PIL_PRODUCT_LEVEL_UNDEF:
                cpl_frame_set_level(self, CPL_FRAME_LEVEL_NONE);
                break;

            case PIL_PRODUCT_LEVEL_TEMPORARY:
                cpl_frame_set_level(self, CPL_FRAME_LEVEL_TEMPORARY);
                break;

            case PIL_PRODUCT_LEVEL_INTERMEDIATE:
                cpl_frame_set_level(self, CPL_FRAME_LEVEL_INTERMEDIATE);
                break;

            case PIL_PRODUCT_LEVEL_SECONDARY:
            case PIL_PRODUCT_LEVEL_PRIMARY:
                cpl_frame_set_level(self, CPL_FRAME_LEVEL_FINAL);
                break;

            default:
                cpl_frame_delete(self);
                return NULL;
                break;
        }

    }

    return self;

}


/**
 * @brief
 *   Export a CPL frameset.
 *
 * @param set  The CPL frameset to export.
 * @param sof  The set of frames to fill.
 * 
 * @return The function returns 0 on success or a non-zero value otherwise.
 *
 * The function fills an existing set of frames object @em sof with the
 * converted contents of the source CPL frame set @em set. If set is empty,
 * or if @em set is @c NULL, the target set of frames is left unchanged.
 *
 * The intended use of this function is to convert an input CPL frameset
 * into a set of frames structure. The function therefore does not deal
 * with other frame and product attributes than the frame group. Calling
 * this function for a frameset containing products with extended
 * attributes will propagate only the frame group to the created set of
 * frames.
 *
 * In most case this is exactly what is needed, since the product attributes
 * are only used to help the calling application to write the frames to
 * the disk. On input to a subsequent data reduction task the file name,
 * the tag and the frame group are sufficient.
 *
 * @see vmCplImportSof()
 */

int
vmCplFramesetExport(const cpl_frameset *set, PilSetOfFrames *sof)
{

    if (sof == NULL) {
        return 1;
    }

    if (set != NULL || !cpl_frameset_is_empty(set)) {

        for (int i =0; i< cpl_frameset_get_size(set); i ++) 
        {
            const cpl_frame * frame = cpl_frameset_get_position_const(set, i);

            if (cpl_frame_get_group(frame) != CPL_FRAME_GROUP_PRODUCT) {

                PilFrame *_frame = _vmCplFrameExport(frame);


                if (_frame == NULL) {
                    return 2;
                }

                if (pilSofInsert(sof, _frame) == 0) {
                    return 3;
                }

            }

        }

    }

    return 0;

}


/**
 * @brief
 *   Imports products to a CPL frameset.
 *
 * @param set  The CPL frameset to be updated.
 * @param sof  The set of frames from which products are imported.
 * 
 * @return The function returns 0 on success or a non-zero value otherwise.
 *
 * The function updates an existing CPL frame set @em set with the product
 * frames found in the set of frames @em sof. The target CPL frame set must
 * exist. If @em set is @c NULL the function returns an error. If the source
 * set of frames @em sof is empty, or if it is @c NULL the target frame set
 * is left unchanged, otherwise all product frames found in @em sof are
 * converted to CPL frames and inserted in @em set.
 *
 * @see vmCplExportFrameset()
 */

int
vmCplFramesetImport(cpl_frameset *set, const PilSetOfFrames *sof)
{

    if (set == NULL) {
        return 1;
    }

    if (sof != NULL || !pilSofIsEmpty(sof)) {

        PilFrame *frame = pilSofFirst((PilSetOfFrames *)sof);

        while (frame != NULL) {
            PilFrameType group = pilFrmGetType(frame);
            const char *name = pilFrmGetName(frame);
            const char *_name;
            cpl_frame *_frame;
            cpl_frame_group _group;
            int nframes = cpl_frameset_get_size(set);
            int i;

            switch (group) {
                case PIL_FRAME_TYPE_UNDEF:
                    _group = CPL_FRAME_GROUP_NONE;
                    break;

                case PIL_FRAME_TYPE_RAW:
                    _group = CPL_FRAME_GROUP_RAW;
                    break;

                case PIL_FRAME_TYPE_CALIB:
                    _group = CPL_FRAME_GROUP_CALIB;
                    break;

                case PIL_FRAME_TYPE_PRODUCT:
                    _group = CPL_FRAME_GROUP_PRODUCT;
                    break;

                default:
                    break;
            }

            for (i = 0; i < nframes; i++) {
                _frame = cpl_frameset_get_position(set, i);
                _name = cpl_frame_get_filename(_frame);
                if (strcmp(name, _name) == 0) {
                    cpl_frame_set_group(_frame, _group);
                    break;
                }
            }

            frame = pilSofNext((PilSetOfFrames *)sof, frame);
        }

        frame = pilSofFirst((PilSetOfFrames *)sof);

        while (frame != NULL) {

            if (pilFrmGetType(frame) == PIL_FRAME_TYPE_PRODUCT) {

                cpl_frame *_frame = _vmCplFrameImport(frame);


                if (_frame == NULL) {
                    return 2;
                }

                cpl_frameset_insert(set, _frame);

            }

            frame = pilSofNext((PilSetOfFrames *)sof, frame);

        }

    }

    return 0;

}


/**
 * @brief
 *   Export a CPL parameter list.
 *
 * @param list  The CPL parameter list to export.
 * 
 * @return The function returns 0 on success or a non-zero value otherwise.
 *
 * The function fills the recipe configuration database with the contents
 * of the CPL parameter list @em list. The recipe configuration database
 * must be initialized before this function is called. If the parameter
 * list is empty, or if @em list is @c NULL, the recipe configuration
 * database is left unchanged.
 */

int
vmCplParlistExport(const cpl_parameterlist *list)
{

    if (list != NULL 
        || cpl_parameterlist_get_size((cpl_parameterlist *)list) == 0) {

        cpl_parameter *parameter = 
                       cpl_parameterlist_get_first((cpl_parameterlist *)list);


        while (parameter) {

            const char *name;
            const char *group;

            char *context = (char *)cpl_parameter_get_context(parameter);

            int status;

            cx_string *value;


            /*
             * Verify that the parameter context starts with 'vimos.'
             */

            if (strstr(context, "vimos.") != context) {
                return -1;
            }
            else {
                group = context + strlen("vimos.");
            }


            /*
             * The parameter's command line alias is used as the name
             * of the corresponding recipe configuration database entry.
             */

            name = cpl_parameter_get_alias(parameter, CPL_PARAMETER_MODE_CLI);

            if (name == NULL) {
                return -2;
            }
            else {
                char *s = strrchr(name, (int)'.');

                if (s != NULL) {
                    name = s + 1;
                }
            }


            /*
             * Create a string representation of the parameter value to
             * be inserted into the recipe configuration data base.
             */

            value = cx_string_new();

            switch (cpl_parameter_get_type(parameter)) {
                case CPL_TYPE_BOOL:
                    if (cpl_parameter_get_bool(parameter) != 0) {
                        cx_string_set(value, "true");
                    }
                    else {
                        cx_string_set(value, "false");
                    }
                    break;

                case CPL_TYPE_INT:
                {
                    int i = cpl_parameter_get_int(parameter);

                    cx_string_sprintf(value, "%i", i);
                    break;
                }

                case CPL_TYPE_DOUBLE:
                {
                    double d = cpl_parameter_get_double(parameter);

                    cx_string_sprintf(value, "%g", d);
                    break;
                }

                case CPL_TYPE_STRING:
                    cx_string_set(value, cpl_parameter_get_string(parameter));
                    break;

                default:
                    return -3;
                    break;
            }


            /*
             * Create the new recipe configuration database entry.
             */

            status = pilDfsDbCreateEntry(group, name, cx_string_get(value),
                                         READ_WRITE);

            if (status != EXIT_SUCCESS) {
                cx_string_delete(value);
                return 1;
            }

            cx_string_delete(value);
            parameter = cpl_parameterlist_get_next((cpl_parameterlist *)list);

        }

    }

    return 0;

}


/**
 * @brief
 *   Convert the messaging severity level.
 *
 * @param level  Message severity level to convert.
 *
 * @return The converted message severity level.
 *
 * The function translates a CPL message severity level to a VIMOS message
 * severity level.
 */

PilMsgSeverity
vmCplMsgSeverityExport(cpl_msg_severity level)
{

    PilMsgSeverity _level;


    switch (level) {
    case CPL_MSG_DEBUG:
        _level = PIL_MSG_DEBUG;
        break;

    case CPL_MSG_INFO:
        _level = PIL_MSG_INFO;
        break;

    case CPL_MSG_WARNING:
        _level = PIL_MSG_WARNING;
        break;

    case CPL_MSG_ERROR:
        _level = PIL_MSG_ERROR;
        break;

    case CPL_MSG_OFF:
        _level = PIL_MSG_OFF;
        break;

    default:
        break;
    }

    return _level;

}


/**
 * @brief
 *   Start the recipe execution timer
 *
 * @param time  Address of the variable where the start time will be stored.
 *
 * @return The function returns 0 on success and a non-zero values in case of
 *   an error.
 *
 * The function starts the recipe's execution timer and saves the start time
 * to the recipe information structure. In addition, if @em time is not
 * @c NULL, the recipe's start time is stored in the variable pointed to by
 * @em time.
 */

int
vmCplRecipeTimerStart(PilTime *time)
{

    PilTime start;
    PilTimer *timer = pilRecGetTimer();


    if (timer == NULL) {
        timer = newPilTimer();

        if (timer == NULL) {
            return 1;
        }

        pilRecSetTimer(timer);
    }
    else {
        if (pilTimerIsActive(timer)) {
            pilTimerStop(timer, NULL);
        }
    }

    start = pilTimerStart(timer, NULL);
    pilRecSetTimeStart(start);

    if (time != NULL) {
        *time = start;
    }

    return 0;

}


/**
 * @brief
 *   Stop the recipe execution timer
 *
 * @param time  Address of the variable where the stop time will be stored.
 *
 * @return The function returns 0 on success and a non-zero values in case of
 *   an error.
 *
 * The function stops the recipe's execution timer and saves the stop time
 * to the recipe information structure. In addition, if @em time is not
 * @c NULL, the recipe's stop time is stored in the variable pointed to by
 * @em time.
 */

int
vmCplRecipeTimerStop(PilTime *time)
{

    PilTimer *timer = pilRecGetTimer();


    if (timer == NULL) {
        return 1;
    }

    if (pilTimerIsActive(timer)) {

        PilTime elapsed;
        PilTime start = pilRecGetTimeStart();

        pilTimerStop(timer, NULL);
        elapsed = pilTimerElapsed(timer, NULL);

        pilRecSetTimeStop(start + elapsed);

        if (time != NULL) {
            *time = start + elapsed;
        }

        return 0;
    }
        
    return 2;

}


/**
 * @brief
 *   Initialize the VIMOS recipe subsystems.
 *
 * @param recipe  The name of the recipe to initialize.
 *
 * @return The function returns 0 on success or  a non-zero value otherwise.
 *
 * The function uses the current CPL subsystem settings to initialize
 * the corresponding VIMOS recipe subsystems in the same way.
 */

int
vmCplRecipeStart(const char *recipe, const char *version)
{

    PilMsgSeverity level;

    cx_print_func printer;


    /*
     * Initialize the recipe information structure
     */

    if (pilRecSetName(recipe) == EXIT_FAILURE) {
        return 1;
    }

    if (pilRecSetVersion(version) == EXIT_FAILURE) {
        return 1;
    }

    if (pilRecSetInstrument("vimos") == EXIT_FAILURE) {
        return 1;
    }


    /*
     * Initialize the messaging subsystem
     */

    pilMsgStart();
    pilMsgSetRecipeName(recipe);

    printer = cx_print_set_handler(NULL);
    cx_print_set_handler(printer);
    pilMsgSetPrintHandler(printer);

    printer = cx_printerr_set_handler(NULL);
    cx_printerr_set_handler(printer);
    pilMsgSetErrorHandler(printer);

    level = PIL_MSG_OFF; /* vmCplMsgSeverityExport(cpl_msg_get_log_level()); */
    pilMsgEnableLog(level);

    if (level == PIL_MSG_DEBUG) {
        pilMsgEnableTimeTag();
        pilMsgEnableComponentTag();
    }

 /* level = vmCplMsgSeverityExport(cpl_msg_get_level()); */
    pilMsgEnableTerminal(level);

    if (level == PIL_MSG_DEBUG) {
        pilMsgEnableTimeTag();
        pilMsgEnableComponentTag();
    }


    /*
     * Initialize the recipe configuration database
     */

    if (pilDfsCreateDB('.', USE_CASE) == EXIT_FAILURE) {
        return 2;
    }


    /*
     * Initialize keyword and category translation tables
     */

    /* FIXME: Add suport for loading translation tables from files
     *        since the old system supports it too. (RP)
     */

    if (pilTrnInitKeywordMap() == EXIT_FAILURE) {
        return 3;
    }

    if (pilTrnInitCategoryMap() == EXIT_FAILURE) {
        return 4;
    }

    return 0;

}


/**
 * @brief
 *   Shutdown the VIMOS recipe subsystems.
 *
 * @return The function returns 0 on success or  a non-zero value otherwise.
 *
 * The function stops the VIMOS recipe subsystems initialized by a call to
 * @b vmCplRecipeStart().
 * the corresponding VIMOS recipe subsystems in the same way.
 */

int
vmCplRecipeStop(void)
{

    /*
     * Services should be disabled in the reverse order of
     * initialization. See vmCplRecipeStart() for details.
     */

    pilTrnClearCategoryMap();
    pilTrnClearKeywordMap();

    pilDfsFreeDB();

    if (pilMsgCloseLog() != EXIT_SUCCESS) {
        return 1;
    }
    pilMsgStop();

    pilRecInfoClear();

    return 0;

}


/**
 * @brief
 *   Apply post processing steps to all frames in a set of frames.
 *
 * @param sof  Set of frames to process.
 *
 * @return The function returns 0 on success, or a non-zero value otherwise.
 *
 * The function finalizes the set of frames after the data reduction task
 * has completed. For registered product frames, for instance, the FITS
 * header entries are completed with generic product information.
 */

int
vmCplPostProcessFrames(PilSetOfFrames *sof, const char *recipename)
{

    // const char *fctid = "vmCplPostProcessFrames";

    PilFrame *frame = NULL;


    if (sof == NULL) {
        return 1;
    }


    /*
     * Process all frames in the set. Apply the type specific
     * post processing task (if any).
     */

    frame = pilSofFirst(sof);

    if (frame == NULL) {
        return 2;
    }


    while (frame) {

        if (pilFrmGetType(frame) == PIL_FRAME_TYPE_PRODUCT) {
            if (pilFrmGetFormat(frame) != PIL_FRAME_FORMAT_PAF) {
/*                if (pilRecUpdateProductInfo(frame, NULL, sof) != 0) {   */
                if (vm_dfs_setup_product_header(frame, recipename, sof) != 0) {
                    return 3;
                }
            }
        }

        frame = pilSofNext(sof, frame);

    }

    return 0;

}
/**@}*/
