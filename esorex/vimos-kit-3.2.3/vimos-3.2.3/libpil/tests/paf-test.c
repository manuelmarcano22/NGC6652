/* $Id: paf-test.c,v 1.1.1.1 2008-10-21 09:10:12 cizzo Exp $
 * ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * COPYRIGHT (c) 2000 European Southern Observatory
 * LICENSE: GNU General Public License version 2 or later
 *
 * PROJECT:   VLT Data Flow System
 * AUTHOR:    Ralf Palsa -- ESO/DMD/DPG
 * SUBSYSTEM: Instrument pipelines
 *
 * PURPOSE:
 *   PAF module tests.
 *
 * DESCRIPTION:
 *
 * $Name: not supported by cvs2svn $
 * $Revision: 1.1.1.1 $
 * ----------------------------------------------------------------------------
 */

#undef NDEBUG

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pilpaf.h"


int
main()
{

    PilPAF *paf;


    /*
     * Test 1: Create and destroy a PAF header.
     */

    paf = newPilPAF("piltest.paf", "Pipeline Result", NULL, NULL);

    assert(pilPAFIsEmpty(paf));
    assert(pilPAFGetSize(paf) == 0);

    pilPAFAppendString(paf, "QC1.DID", "VIMOS-1.0",
                       "VIMOS QC1 dictionary version");
    pilPAFAppendDouble(paf, "DATAMIN", 1.23456789123456789e-10, NULL);
    pilPAFAppendInt(paf, "DET.OUTPUTS", 2, NULL);
    pilPAFAppendBool(paf, "SIMPLE", 1, NULL);

    pilPAFPrependString(paf, "QC1.DID", "VIMOS-1.0",
                       "VIMOS QC1 dictionary version");
    pilPAFPrependDouble(paf, "DATAMIN", 300., NULL);
    pilPAFPrependInt(paf, "DET.OUTPUTS", 2, NULL);
    pilPAFPrependBool(paf, "SIMPLE", 1, NULL);

    pilPAFInsertBool(paf, "QC1.DID", "SIMPLE", 0, NULL);
    pilPAFInsertInt(paf, "QC1.DID", "DET.CHIPS", 1, NULL);
    pilPAFInsertDouble(paf, "QC1.DID", "MJD-OBS", 52341.98213595,
                       "Obs start 2002-03-08T23:34:16.546");
    pilPAFInsertString(paf, "QC1.DID", "DATE-OBS", "2002-03-08T23:34:16.546",
                       NULL);

    pilPAFInsertAfterBool(paf, "QC1.DID", "SIMPLE", 0, NULL);
    pilPAFInsertAfterInt(paf, "QC1.DID", "DET.CHIPS", 1, NULL);
    pilPAFInsertAfterDouble(paf, "QC1.DID", "MJD-OBS", 52341.98213595,
                       "Obs start 2002-03-08T23:34:16.546");
    pilPAFInsertAfterString(paf, "QC1.DID", "DATE-OBS",
                            "2002-03-08T23:34:16.546", NULL);

    assert(pilPAFSetValueBool(paf, "SIMPLE", 0) == EXIT_SUCCESS);
    assert(pilPAFGetValueBool(paf, "SIMPLE") == 0);
    assert(pilPAFSetValueInt(paf, "DET.OUTPUTS", 4) == EXIT_SUCCESS);
    assert(pilPAFGetValueInt(paf, "DET.OUTPUTS") == 4);
    assert(pilPAFSetValueDouble(paf, "DATAMIN", 195.) == EXIT_SUCCESS);
    assert(pilPAFGetValueDouble(paf, "DATAMIN") == 195.);
    assert(pilPAFSetValueString(paf, "QC1.DID", "VIMOS-2.0") == EXIT_SUCCESS);
    assert(!strcmp(pilPAFGetValueString(paf, "QC1.DID"), "VIMOS-2.0"));

    pilPAFWrite(paf);
    deletePilPAF(paf);

    /*
     * All tests succeeded
     */

    return 0;

}
