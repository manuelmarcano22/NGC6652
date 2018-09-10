#include "casu_utilfunctions.h"

#include "imcore.h"

#include "hdrl_utils.h"

#include <cpl.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


//raw copy from 'casu_wcsutils.c':
   /* WCS keywords that should be removed from FITS tables*/
#define SZKEY 9  //seems not used
#define NNOTABKEYS 6
static const char *notabkeys[NNOTABKEYS] = {"^CRVAL[1-2]*$","^CRPIX[1-2]*",
                                            "^CD[1-2]*_[1-2]*","^CDELT[1-2]*",
                                            "^CTYPE[1-2]*","^PV[1-9]*_[1-9]*"};

static void hdrl_update_radec(casu_tfits* tcat, const cpl_wcs * wcs);

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_get_ehu
    \par Purpose:
        Get the propertylist for the extension header for a given casu_fits 
        image.
    \par Description:
        Get the propertylist for the extension header for a given casu_fits 
        image. This is the extension that is relevant of the image.
        This should only need to be read once and then can be used to add
        things to the primary header.
    \par Language:
        C
    \param p
        The input casu_fits object
    \return 
        The propertylist representing the extension header of the input image
        (NULL if there is an error).
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern cpl_propertylist *casu_fits_get_ehu(casu_fits *p) {

    /* Check for nonsense input */

    if (p == NULL)
        return(NULL);

    /* If the propertylist hasn't already been loaded, then do it now */

    if (p->ehu == NULL) 
        p->ehu = cpl_propertylist_load(p->fname,(cpl_size)(p->nexten));
    
    /* Return it */

    return(p->ehu);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_tfits_get_ehu
    \par Purpose:
        Get the propertylist for the extension header for a given casu_tfits 
        image.
    \par Description:
        Get the propertylist for the extension header for a given casu_tfits 
        image. This is the extension that is relevant of the image.
        This should only need to be read once and then can be used to add
        things to the primary header.
    \par Language:
        C
    \param p
        The input casu_tfits object
    \return 
        The propertylist represeting the extension header of the input table
        (NULL if there is an error).
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern cpl_propertylist *casu_tfits_get_ehu(casu_tfits *p) {

    /* Check for nonsense input */

    if (p == NULL)
        return(NULL);

    /* If the propertylist hasn't already been loaded, then do it now */

    if (p->ehu == NULL) 
        p->ehu = cpl_propertylist_load(p->fname,(cpl_size)(p->nexten));
    
    /* Return it */

    return(p->ehu);
}

/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_get_image
    \par Purpose:
        Get the CPL image from the casu_fits object
    \par Description:
        Return the CPL image from the input casu_fits object. This image is
        suitable for use in all cpl_image routines.
    \par Language:
        C
    \param p
        The input casu_fits object
    \return 
        The cpl_image object. NULL if there was an error.
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern cpl_image *casu_fits_get_image(casu_fits *p) {
    const char *fctid = "casu_fits_get_image";
    cpl_image *im2;
    
    /* Check for nonsense input */

    if (p == NULL)
        return(NULL);
    
    /* If the image is already loaded, then return it now */

    if (p->image != NULL)
        return(p->image);

    /* If the image hasn't been loaded, then do that now */
    
    if (p->image == NULL) 
        p->image = cpl_image_load(p->fname,p->type,0,(cpl_size)(p->nexten));
    if (p->image == NULL) {
        cpl_msg_error(fctid,"Unable to load %s[%" CPL_SIZE_FORMAT "] -- %s\n",
                      p->fname,(cpl_size)(p->nexten),
                      cpl_error_get_message());
        cpl_error_reset();
        return(NULL);
    } 

    /* If this is an unspecified type, the force it to be float */

    if (p->type == CPL_TYPE_UNSPECIFIED && 
        cpl_image_get_type(p->image) != CPL_TYPE_FLOAT) {
        im2 = cpl_image_cast(p->image,CPL_TYPE_FLOAT);
        cpl_image_delete(p->image);
        p->image = im2;
    }

    /* Return it */

    return(p->image);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_get_nexten
    \par Purpose:
        Get the FITS extension number for the current image in a casu_fits 
        object
    \par Description:
        Get the FITS extension number for the current image in a casu_fits 
        object
    \par Language:
        C
    \param p
        The input casu_fits object
    \return 
        The extension number (-1 in case of error)
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern intptr_t casu_fits_get_nexten(casu_fits *p) {
    
    /* Check for nonsense input */

    if (p == NULL)
        return(-1);

    /* Return it */

    return(p->nexten);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_get_phu
    \par Purpose:
        Get the propertylist for the primary header for a given casu_fits image.
    \par Description:
        Get the propertylist for the primary header for a given casu_fits image.
        This should only need to be read once and then can be used to add
        things to the primary header.
    \par Language:
        C
    \param p
        The input casu_fits object
    \return 
        The propertylist representing the primary header of the input image
        (NULL if there is an error).
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern cpl_propertylist *casu_fits_get_phu(casu_fits *p) {

    /* Check for nonsense input */

    if (p == NULL)
        return(NULL);

    /* If the propertylist hasn't already been loaded, then do it now */

    if (p->phu == NULL) {
        if (p->casufitstype == CASU_FITS_MEF) 
            p->phu = cpl_propertylist_load(p->fname,0);
        else
            p->phu = cpl_propertylist_load(p->fname,p->nexten);
    }
    
    /* Return it */

    return(p->phu);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_tfits_get_phu
    \par Purpose:
        Get the propertylist for the primary header for a given casu_tfits 
        image.
    \par Description:
        Get the propertylist for the primary header for a given casu_tfits 
        image. This should only need to be read once and then can be used 
        to add things to the primary header.
    \par Language:
        C
    \param p
        The input casu_tfits object
    \return 
        The propertylist represeting the primary header of the input table
        (NULL if there is an error).
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern cpl_propertylist *casu_tfits_get_phu(casu_tfits *p) {

    /* Check for nonsense input */

    if (p == NULL)
        return(NULL);

    /* If the propertylist hasn't already been loaded, then do it now */

    if (p->phu == NULL) 
        p->phu = cpl_propertylist_load(p->fname,0);
    
    /* Return it */

    return(p->phu);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_tabwcs
    \par Purpose:
        Remove FITS image WCS keywords from a propertylist and replace with
        tabular keywords
    \par Description:
        Replace the FITS image WCS keywords in a propertylist with the
        relevant table FITS keyword. This is not very general
    \par Language:
        C
    \param p
        The input propertylist
    \param xcol
        The column number for the X position
    \param ycol
        The column number for the Y position
    \param status
        Standard input and output casu status variable
    \return 
        Standard casu status variable
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int casu_tabwcs(cpl_propertylist *p, intptr_t xcol, intptr_t ycol,
			int *status) {
    intptr_t i;
    char key[9],key2[9];
    const char *fctid="casu_tabwcs";

    /* Inherited status */

    if (*status != CASU_OK)
        return(*status);
    if (p == NULL) {
        cpl_msg_error(fctid,"Propertylist passed is NULL\nProgramming error");
        FATAL_ERROR
    }

    /* If either of the columns are -1, then just get rid of the image WCS 
       and get out of here */

    if (xcol == -1 || ycol == -1) {
        casu_removewcs(p,status);
        GOOD_STATUS
    }

    /* Go through the standard WCS header keywords one by one
       and translate them.  Start with CTYPE */

    (void)snprintf(key,8,"TCTYP%zd",xcol);
    casu_rename_property(p,"CTYPE1",key);
    (void)snprintf(key,8,"TCTYP%zd",ycol);
    casu_rename_property(p,"CTYPE2",key);


    /* Now CRVAL */

    (void)snprintf(key,8,"TCRVL%zd",xcol);
    casu_rename_property(p,"CRVAL1",key);
    (void)snprintf(key,8,"TCRVL%zd",ycol);
    casu_rename_property(p,"CRVAL2",key);

    /* Now CRPIX */

    (void)snprintf(key,8,"TCRPX%zd",xcol);
    casu_rename_property(p,"CRPIX1",key);
    (void)snprintf(key,8,"TCRPX%zd",ycol);
    casu_rename_property(p,"CRPIX2",key);

    /* Now PV matrix */

    for (i = 1; i <= 5; i++) {
        (void)snprintf(key2,8,"PV2_%zd",i);
        (void)snprintf(key,8,"TV%zd_%zd",ycol,i);
        if (cpl_propertylist_has(p,key2))
            casu_rename_property(p,key2,key);
    }

    /* Now the CD matrix */

    (void)snprintf(key,8,"TC%zd_%zd",xcol,xcol);
    casu_rename_property(p,"CD1_1",key);
    (void)snprintf(key,8,"TC%zd_%zd",xcol,ycol);
    casu_rename_property(p,"CD1_2",key);
    (void)snprintf(key,8,"TC%zd_%zd",ycol,xcol);
    casu_rename_property(p,"CD2_1",key);
    (void)snprintf(key,8,"TC%zd_%zd",ycol,ycol);
    casu_rename_property(p,"CD2_2",key);

    /* Now get out of here */

    GOOD_STATUS
    
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_removewcs
    \par Purpose:
        Remove FITS image WCS keywords from a propertylist
    \par Description:
        Remove FITS WCS keywords from a propertylist. This is sometimes
        necessary if a FITS table header has been based on an image header.
    \par Language:
        C
    \param p
        The input propertylist
    \param status
        Standard input and output casu status variable
    \return 
        Standard casu status variable
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern int casu_removewcs(cpl_propertylist *p, int *status) {
    intptr_t i;
    const char *fctid = "casu_removewcs";

    /* Inherited status */

    if (*status != CASU_OK)
        return(*status);
    if (p == NULL) {
        cpl_msg_error(fctid,"Propertylist passed is NULL\nProgramming error");
        FATAL_ERROR
    }

    /* Loop through all the template keywords and remove them */

    for (i = 0; i < NNOTABKEYS; i++) 
        cpl_propertylist_erase_regexp(p,notabkeys[i],0);

    GOOD_STATUS
}

/*---------------------------------------------------------------------------*/
/**
    \brief Rename a property in a given propertylist

    \par Name:
        casu_rename_property
    \par Purpose:
        Rename a property in a propertylist
    \par Description:
        Rename a property in a propertylist
    \par Language:
        C
    \param p
        The input propertylist for the table
    \param oldname
        The old property name
    \param newname
        The new property name
    \returns
        Nothing
    \author
        Jim Lewis, CASU
*/
/*---------------------------------------------------------------------------*/

extern void casu_rename_property(cpl_propertylist *p, const char *oldname,
                                 char *newname) {
    cpl_propertylist *temp;
    cpl_property *property;

    /* First get the old property. Note that we have to do this in a 
       particularly silly way since you cannot reference an individual property
       in a propertylist by its name. You can only do it by its index 
       number. Remeber to change this when CPL comes to it's senses... */

    if (! cpl_propertylist_has(p,oldname))
        return;
    temp = cpl_propertylist_new();
    cpl_propertylist_copy_property(temp,p,oldname);
    property = cpl_propertylist_get(temp,0);

    /* Now change its name */

    cpl_property_set_name(property,newname);

    /* Now insert this into the propertylist and delete the old one */

    cpl_propertylist_append(p,temp);
    cpl_propertylist_erase(p,oldname);
    cpl_propertylist_delete(temp);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_tfits_get_table
    \par Purpose:
        Get the CPL table from the casu_tfits object
    \par Description:
        Return the CPL table from the input casu_tfits object. This table is
        suitable for use in all cpl_table routines.
    \par Language:
        C
    \param p
        The input casu_tfits object
    \return 
        The cpl_image object. NULL if there was an error.
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern cpl_table *casu_tfits_get_table(casu_tfits *p) {
    
    /* Check for nonsense input */

    if (p == NULL)
        return(NULL);

    /* Return it */

    return(p->table);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_tfits_wrap
    \par Purpose:
        Wrap an table in a casu_tfits wrapper
    \par Description:
        The input table is inserted into a casu_tfits wrapper. A model
        casu_tfits object may be provided to give the new object 
        headers. If the phu and ehu parameters are not null then they will
        be used as the propertylists for the new object. If not, then
        an attempt will be made to copy the propertylists from the model.
    \par Language:
        C
    \param tab
        The input cpl_table
    \param model
        The input casu_tfits model object
    \param phu
        The input propertylist for the primary header for the new object.
    \param ehu
        The input propertylist for the extension header for the new object.
    \return
        The new casu_tfits structure.
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern casu_tfits *casu_tfits_wrap(cpl_table *tab, casu_tfits *model,
                                   cpl_propertylist *phu,
                                   cpl_propertylist *ehu) {
    casu_tfits *p;

    /* Check for nonsense input */

    if (tab == NULL)
        return(NULL);

    /* Get the casu_tfits structure */

    p = cpl_malloc(sizeof(casu_tfits));

    /* Load stuff in */

    p->table = tab;
    p->nexten = -1;
    if (phu != NULL)
        p->phu = phu;
    else if (model != NULL)
        p->phu = cpl_propertylist_duplicate(casu_tfits_get_phu(model));
    else
        p->phu = cpl_propertylist_new();
    if (ehu != NULL)
        p->ehu = ehu;
    else if (model != NULL)
        p->ehu = cpl_propertylist_duplicate(casu_tfits_get_ehu(model));
    else
        p->ehu = cpl_propertylist_new();
    p->fname = NULL;
    p->status = CASU_OK;
    p->extname = NULL;
    p->fullname = NULL;

    /* Get out of here */

    return(p);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_wrap
    \par Purpose:
        Wrap an image in a casu_fits wrapper
    \par Description:
        The input image is inserted into a casu_fits wrapper. A model
        casu_fits object may be provided to give the new object
        headers. If the phu and ehu parameters are not null then they will
        be used as the propertylists for the new object. If not, then
        an attempt will be made to copy the propertylists from the model.
        If the model and the propertylists are both NULL, then empty
        propertylists are given
    \par Language:
        C
    \param im
        The input cpl_image
    \param model
        The input casu_fits model object
    \param phu
        The input propertylist for the extension header for the new object.
    \param ehu
        The input propertylist for the extension header for the new object.
    \return 
        The new casu_fits structure.
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern casu_fits *casu_fits_wrap(cpl_image *im, casu_fits *model, 
                                 cpl_propertylist *phu, cpl_propertylist *ehu) {
    casu_fits *p;

    /* Check for nonsense input */

    if (im == NULL)
        return(NULL);

    /* Get the casu_fits structure */

    p = cpl_malloc(sizeof(casu_fits));

    /* Load stuff in */

    p->image = im;
    p->nexten = -1;
    if (phu != NULL) 
        p->phu = cpl_propertylist_duplicate(phu);
    else if (model != NULL) 
        p->phu = cpl_propertylist_duplicate(casu_fits_get_phu(model));
    else
        p->phu = cpl_propertylist_new();
    if (ehu != NULL)
        p->ehu = cpl_propertylist_duplicate(ehu);
    else if (model != NULL) 
        p->ehu = cpl_propertylist_duplicate(casu_fits_get_ehu(model));
    else 
        p->ehu = cpl_propertylist_new();
    p->fname = NULL;
    p->status = CASU_OK;
    p->extname = NULL;
    p->fullname = NULL;
    if (model != NULL) 
        p->casufitstype = model->casufitstype;
    else 
        p->casufitstype = CASU_FITS_MEF;
    p->type = cpl_image_get_type(im);
   
    /* Get out of here */

    return(p);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_delete
    \par Purpose:
        Free all the workspace associated with a casu_fits object
    \par Description:
        Free all the workspace associated with a casu_fits object
    \par Language:
        C
    \param p
        The input casu_fits object
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern void casu_fits_delete(casu_fits *p) {

    /* Check for nonsense input */

    if (p == NULL)
        return;

    /* Free up workspace if it's been used */

    freeimage(p->image);
    freepropertylist(p->phu);
    freepropertylist(p->ehu);
    freespace(p->fname);
    freespace(p->extname);
    freespace(p->fullname);
    cpl_free(p);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_tfits_delete
    \par Purpose:
        Free all the workspace associated with a casu_tfits object
    \par Description:
        Free all the workspace associated with a casu_tfits object
    \par Language:
        C
    \param p
        The input casu_tfits object
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern void casu_tfits_delete(casu_tfits *p) {

    /* Check for nonsense input */

    if (p == NULL)
        return;

    /* Free up workspace if it's been used */

    freetable(p->table);
    freepropertylist(p->phu);
    freepropertylist(p->ehu);
    freespace(p->fname);
    freespace(p->extname);
    freespace(p->fullname);
    cpl_free(p);
}
/*---------------------------------------------------------------------------*/
/**
    \ingroup casu_modules
    \brief Generate object catalogues from input images

    \par Name:
        casu_imcore
    \par Purpose:
        Generate object catalogues from input images
    \par Description:
        A frame and its confidence map are given. Detection thresholds and
        various other parameters are also given. Output is a table with all
        the extracted objects with object classifications included.
    \par Language:
        C
    \param infile
        The input frame with the image to be analysed
    \param conf
        The input frame with the confidence map
    \param wcs
        TODO
    \param ipix
        The minimum allowable size of an object
    \param threshold
        The detection threshold in sigma above sky
    \param icrowd
        If set then the deblending software will be used
    \param rcore
        The core radius in pixels
    \param nbsize
        The smoothing box size for background map estimation
    \param cattype
        The type of catalogue to be produced
    \param filtfwhm
        The FWHM of the smoothing kernel in the detection algorithm
    \param outtab
        The output table of object
    \param gainloc
        The detector gain in e-/ADU
    \param res
        TODO
    \param status
        The input/output status that has the same return value a below
    \retval CASU_OK 
        if everything is ok
    \retval CASU_WARN,CASU_FATAL
        errors in the called routines
    \par QC headers:
        The following values will go into the extension propertylist
        - \b SATURATION
            Saturation level in ADU
        - \b MEAN_SKY
            Mean sky brightness in ADU
        - \b SKY_NOISE
            Pixel noise at sky level in ADU
        - \b IMAGE_SIZE
            The average FWHM of stellar objects in pixels
        - \b ELLIPTICITY
            The average stellar ellipticity (1 - b/a)
        - \b APERTURE_CORR
            The stellar aperture correction for 1x core flux
        - \b NOISE_OBJ
            The number of noise objects
    \par Other headers:
        The following values will go into the extension propertylist
        - \b APCORxx
            A series of aperture correction values for each of the core
            radius apertures.
        - \b SYMBOLx
            A series of keywords to be used by GAIA for plotting ellipses
    \author
        Jim Lewis, CASU

 */
/*---------------------------------------------------------------------------*/

extern int hdrl_casu_imcore(casu_fits *infile, casu_fits *conf, const cpl_wcs * wcs,
                       intptr_t ipix,
                       float threshold, intptr_t icrowd, float rcore,
                       int bkg_subtr, intptr_t nbsize,
                       hdrl_catalogue_options cattype, float filtfwhm, float gainloc,
                       float saturation, hdrl_imcore_result * res,
                       int *status)
{
    int retval;
    const char *fctid = "casu_imcore";
    cpl_propertylist *plist;
    casu_fits *in;

    /* Inherited status */

    res->catalogue = NULL;
    if (*status != CASU_OK)
        return(*status);

    /* Copy the input, the background is subtracted in-place */

    in = casu_fits_duplicate(infile);

    /* Call the main processing routine and get the catalogue */


    retval = imcore_conf(in, conf, ipix, threshold, icrowd, rcore, bkg_subtr,
			 nbsize, cattype,
                         filtfwhm, gainloc, saturation, res);

    if (retval != CASU_OK)
        FATAL_ERROR;
    if ((intptr_t)cpl_table_get_nrow(casu_tfits_get_table(res->catalogue)) == 0) {
        cpl_msg_warning(fctid,"No objects found in image");
        casu_fits_delete(in);
        WARN_RETURN
    }

   /* Get the property list from the input frame */

    plist = casu_fits_get_phu(infile);
    if (plist == NULL) {
        cpl_msg_error(fctid,"Unable to open propertylist %s",
                      casu_fits_get_filename(infile));
        casu_fits_delete(in);
        FATAL_ERROR
    }

    /* Do the classification */
    if (cattype & HDRL_CATALOGUE_CAT_COMPLETE) {
        retval = imcore_classify(res->catalogue,plist,16.0,cattype);
        if (retval != CASU_OK) {
            casu_fits_delete(in);
            WARN_RETURN
        }

        /* Update the RA and DEC of the objects in the object catalogue */
        if (wcs) {
            hdrl_update_radec(res->catalogue, wcs);
        }

        cpl_propertylist_set_comment(casu_tfits_get_ehu(res->catalogue),
                                     "ESO QC IMAGE_SIZE",
                                     "[pixel] Average FWHM of stellar objects");
    }
    else {
        cpl_table_select_all(casu_tfits_get_table(res->catalogue));
        cpl_table_erase_selected(casu_tfits_get_table(res->catalogue));
    }
    casu_fits_delete(in);
    GOOD_STATUS
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_duplicate
    \par Purpose:
        Copy a casu_fits structure into another one.
    \par Description:
        An input casu_fits structure is duplcated and returned
    \par Language:
        C
    \param in
        The input casu_fits object
    \return
        The output casu_fits object.
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern casu_fits *casu_fits_duplicate(casu_fits *in) {
    casu_fits *p;

    /* Check for nonsense input */

    if (in == NULL)
        return(NULL);

    /* Get the casu_fits structure */

    p = cpl_malloc(sizeof(casu_fits));

    /* Now copy everything over */

    if (in->image != NULL)
        p->image = cpl_image_duplicate(in->image);
    else
        p->image = NULL;
    p->phu = cpl_propertylist_duplicate(casu_fits_get_phu(in));
    p->ehu = cpl_propertylist_duplicate(casu_fits_get_ehu(in));
    p->fname = cpl_strdup(in->fname);
    p->extname = cpl_strdup(in->extname);
    p->fullname = cpl_strdup(in->fullname);
    p->nexten = in->nexten;
    p->status = in->status;
    p->casufitstype = in->casufitstype;
    p->type = in->type;
   
    /* Get out of here */

    return(p);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_get_fullname
    \par Purpose:
        Get the fullname of the FITS extension from which the current 
        casu_fits object originated
    \par Description:
        Get the fullname of the FITS extension  from which the current 
        casu_fits object originated. If this is null, then the image didn't 
        originate in an FITS file.
    \par Language:
        C
    \param p
        The input casu_fits object
    \return 
        The fullname name of the file from which this image originated
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern char *casu_fits_get_fullname(casu_fits *p) {

    /* Check for nonsense input */

    if (p == NULL)
        return(NULL);

    /* Return it */

    return(p->fullname);
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_fits_get_filename
    \par Purpose:
        Get the filename from which the current casu_fits object originated
    \par Description:
        Get the filename from which the current casu_fits object originated. If
        this is null, then the image didn't originate in an FITS file.
    \par Language:
        C
    \param p
        The input casu_fits object
    \return 
        The name of the file from which this image originated. This is
        not a copy of the name, so be careful not to modify the string.
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern char *casu_fits_get_filename(casu_fits *p) {

    /* Check for nonsense input */

    if (p == NULL)
        return(NULL);

    /* Return it */

    return(p->fname);
}

/* ---------------------------------------------------------------------------*/
/**
 * @brief update RA and DEC by using the passed wcs information
 * @param tcat The input casu_tfits object
 * @param wcs The wcs paramters for x,y to ra,dec transformation
 * @return
 */
/* ---------------------------------------------------------------------------*/
static void hdrl_update_radec(casu_tfits* tcat, const cpl_wcs * wcs)
{

    /* Update the RA and DEC of the objects in the object catalogue */
    double r, d;
    cpl_table* cat = casu_tfits_get_table(tcat);
    //cpl_table_dump_structure(cat, stdout);
    intptr_t n = (intptr_t) cpl_table_get_nrow(cat);

    float* x = cpl_table_get_data_float(cat, "X_coordinate");
    float* y = cpl_table_get_data_float(cat, "Y_coordinate");
    double* ra = cpl_table_get_data_double(cat, "RA");
    double* dec = cpl_table_get_data_double(cat, "DEC");
    for (intptr_t i = 0; i < n; i++) {
        casu_xytoradec(wcs, (double) x[i], (double) y[i], &r, &d);
        ra[i] = (double) r;
        dec[i] = (double) d;
    }
}
/*---------------------------------------------------------------------------*/
/**
    \par Name:
        casu_xytoradec
    \par Purpose:
        Convert x,y -> ra,dec
    \par Description:
        A WCS structure is used to convert input x,y coordinates
        to equatorial coordinates.
    \par Language:
        C
    \param wcs
        Input WCS structure
    \param x
        Input X
    \param y
        Input Y
    \param ra
        Output RA
    \param dec
        Output Dec
    \return
        Nothing
    \author
        Jim Lewis, CASU
 */
/*---------------------------------------------------------------------------*/

extern void casu_xytoradec(const cpl_wcs *wcs, double x, double y, double *ra,
                           double *dec) {
    double *xy, *radec;
    cpl_matrix *from, *to;
    cpl_array *status;

    /* Load up the information */
    
    from = cpl_matrix_new(1, 2);
    xy = cpl_matrix_get_data(from);
    xy[0] = x;
    xy[1] = y;

    /* Call the conversion routine */

    cpl_wcs_convert(wcs, from, &to, &status, CPL_WCS_PHYS2WORLD);

    /* Pass it back now */

    radec = cpl_matrix_get_data(to);
    *ra = radec[0];
    *dec = radec[1];

    /* Tidy and exit */

    cpl_matrix_delete(from);
    cpl_matrix_delete(to);
    cpl_array_delete(status);
    return;
}



