/*
 * This file is part of the CASU Pipeline utilities
 * Copyright (C) 2015,2016 European Southern Observatory
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

#include "imcore.h"
#include "hdrl.h"
#include <stdio.h>
#include <string.h>


/**@{*/

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Initialise catalogues

    \par Name:
        imcore_tabinit
    \par Purpose:
        Initialise the output table.
    \par Description:
        Wrapper routine to call the relevant initialisation routine for
        each of the allowed types of catalogues.
    \par Language:
        C
    \param ap
        The current ap structure
    \param xcol
        TODO
    \param ycol
        TODO
    \param cattype
        The type of catalogue to be produced
    \param tab
        TODO
    \param res
        TODO
    \returns
        Nothing
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern void imcore_tabinit(ap_t *ap, intptr_t *xcol, intptr_t *ycol,
                           hdrl_catalogue_options cattype,
                           cpl_table **tab, hdrl_imcore_result * res)
{
    imcore_tabinit_6(xcol,ycol,tab);
    if (cattype & HDRL_CATALOGUE_SEGMAP) {
        res->segmentation_map = cpl_image_new(ap->lsiz,ap->csiz,CPL_TYPE_INT);
    }
    else {
        res->segmentation_map = NULL;
    }
    if (cattype & HDRL_CATALOGUE_BKG) {
        res->background = cpl_image_new(ap->lsiz,ap->csiz,CPL_TYPE_FLOAT);
    }
    else {
        res->background = NULL;
    }
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do seeing estimate

    \par Name:
        imcore_do_seeing
    \par Purpose:
        Do the seeing estimate
    \par Description:
        Wrapper routine to call the relevant routine to work out the seeing
        for each of the allowed types of catalogues
    \par Language:
        C
    \param ap
        The current ap structure
    \param cattype
        The type of catalogue to be produced
    \param nobjects
        TODO
    \param tab
        TODO
    \retval CASU_OK
        If all went well
    \retval CASU_FATAL
        If catalogue type is unrecognised
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_do_seeing(ap_t *ap, int cattype, intptr_t nobjects,
                            cpl_table *tab) {
    int status;

    status = imcore_do_seeing_6(ap,nobjects,tab);
    return(status);
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Process results

    \par Name:
        imcore_process_results
    \par Purpose:
        Process the results for each object and store them in the table
    \par Description:
        Wrapper routine to call the relevant routine to work out the results
        for each of the allowed types of catalogues
    \par Language:
        C
    \param ap
        The current ap structure
    \param cattype
        The type of catalogue to be produced
    \param gain
        The header keyword with the gain in e-/ADU
    \param nobjects
        TODO
    \param tab
        TODO
    \param res
        TODO
    \retval CASU_OK
        If all went well
    \retval CASU_FATAL
        If catalogue type is unrecognised
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_process_results(ap_t *ap, int cattype, float gain,
                                  intptr_t *nobjects, cpl_table *tab,
                                  hdrl_imcore_result * res) {
    int status;

    status = imcore_process_results_6(ap,gain,nobjects,tab, res);

    return(status);
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Close the table structure

    \par Name:
        imcore_tabclose
    \par Purpose:
        Close the table structure
    \par Description:
        Wrapper routine to call the relevant routine to close the table
        for each of the allowed types of catalogues
    \par Language:
        C
    \param ap
        The current ap structure
    \param cattype
        The type of catalogue to be produced
    \retval CASU_OK
        If all went well
    \retval CASU_FATAL
        If catalogue type is unrecognised
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_tabclose(ap_t *ap, int cattype) {
    return CASU_OK;
}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Initialise tables (generic)

    \par Name:
        imcore_tabinit_gen
    \par Purpose:
        Initialise tables (generic)
    \par Description:
        Generic routine to create FITS tables for the output catalogues
    \par Language:
        C
    \param ncols
        The number of columns in the table
    \param ttype
        Array of column names for FITS table
    \param tunit
        Array of units for each of the columns
    \param tform
        Array of formats for each of the columns as defined in the FITS
        standard
    \param tab
        TODO
    \returns
        Nothing
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern void imcore_tabinit_gen(intptr_t ncols, const char *ttype[],
                               const char *tunit[], cpl_type tform[],
                               cpl_table **tab) {
    intptr_t i;
    const char *fctid = "imcore_tabinit_gen";

    /* First, create the table with a default number of rows. */

    if ((*tab = cpl_table_new(0)) == NULL) {
        cpl_msg_error(fctid,"Unable to open cpl table!");
        return;
    }

    /* Now define all of the columns */

    for (i = 0; i < ncols; i++) {
        cpl_table_new_column(*tab,ttype[i],tform[i]);
        cpl_table_set_column_unit(*tab,ttype[i],tunit[i]);
    }

}

/*---------------------------------------------------------------------------*/
/**
    \ingroup cataloguemodules
    \brief Do seeing estimate (generic)

    \par Name:
        imcore_do_seeing_gen
    \par Purpose:
        Do seeing estimate (generic)
    \par Description:
        Wrapper routine for doing the seeing estimate
    \par Language:
        C
    \param ap
        The current ap structure
    \param col_ellipt
        The name of the column for ellipticity
    \param col_pkht
        The name of the column for the peak height
    \param col_areals
        The array of names of the areal profile columns
    \param nobjects
        TODO
    \param tab
        TODO
    \retval CASU_OK
        If all is ok. This is currently the only value.
    \par QC headers:
        None
    \par DRS headers:
        None
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int imcore_do_seeing_gen(ap_t *ap, const char *col_ellipt,
                         const char *col_pkht, char *col_areals[NAREAL],
                         intptr_t nobjects, cpl_table *tab) {
    intptr_t i;
    float fwhm,*areal[NAREAL],*ellipt,*pkht,*work;

    /* Get some space and read the relevant columns */
    /* CONSTANTS: 3 */
    work = NULL;
    if (nobjects >= 3) {
        ellipt = cpl_table_get_data_float(tab,col_ellipt);
        pkht = cpl_table_get_data_float(tab,col_pkht);
        work = cpl_malloc(nobjects*sizeof(*work));
        for (i = 0; i < NAREAL; i++)
            areal[i] = cpl_table_get_data_float(tab,col_areals[i]);

        /* Do the seeing calculation */

        imcore_seeing(ap,nobjects,ellipt,pkht,areal,work,&fwhm);
    } else {
        fwhm = 0.0;
    }
    ap->fwhm = fwhm;

    /* Get out of here */

    freespace(work);
    return(CASU_OK);
}

/**@}*/
