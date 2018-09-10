
#ifndef CASU_UTILFUNCTIONS_H
#define CASU_UTILFUNCTIONS_H

#include "hdrl.h"
#include <cpl.h>

/* Define error return values */
#define CASU_OK 0
#define CASU_WARN 1
#define CASU_FATAL 2

/* Useful error handling macros for routines with inherited status */
#define FATAL_ERROR {*status = CASU_FATAL; return(*status);}
#define WARN_RETURN {*status = CASU_WARN; return(*status);}
#define WARN_CONTINUE {*status = CASU_WARN;}
#define GOOD_STATUS {*status = CASU_OK; return(*status);}

/* Define size for recipe description character variables */

#define SZ_ALLDESC 4096

/* Macros to define the type of FITS file we might have here */
#define CASU_FITS_MEF        0
#define CASU_FITS_SIMPLE     1
#define CASU_FITS_SIMPLE_CMP 2
#define CASU_FITS_MEF_NOPHU  3

/* !!! Macros !!! */
#define casu_nint(_x) ((intptr_t)(_x + ((_x) < 0.0 ? -0.5 : 0.5)))
#define min(_x,_y) (_x < _y ? _x : _y)
#define max(_x,_y) (_x > _y ? _x : _y) 

#define freespace(_p) if (_p != NULL) {cpl_free(_p); _p = NULL;}
#define freespace2(_p,_n) if (_p != NULL) {int _nfs2; for (_nfs2 = 0; _nfs2 < _n; _nfs2++) {freespace(_p[_nfs2])}}; freespace(_p);
#define freeframe(_p) if (_p != NULL) {cpl_frame_delete(_p); _p = NULL;}
#define freeimage(_p) if (_p != NULL) {cpl_image_delete(_p); _p = NULL;}
#define freeframeset(_p) if (_p != NULL) {cpl_frameset_delete(_p); _p = NULL;}
#define freeplist(_p) if (_p != NULL) {cpl_propertylist_delete(_p); _p = NULL;}
#define freetable(_p) if (_p != NULL) {cpl_table_delete(_p); _p = NULL;}
#define freepropertylist(_p) if (_p != NULL) {cpl_propertylist_delete(_p); _p = NULL;}
#define freeimagelist(_p) if (_p != NULL) {cpl_imagelist_delete(_p); _p = NULL;}
#define freearray(_p) if (_p != NULL) {cpl_array_delete(_p); _p = NULL;}
#define freewcs(_p) if (_p != NULL) {cpl_wcs_delete(_p); _p = NULL;}
#define closefile(_p) if (_p != NULL) {fclose(_p); _p = NULL;}
#define freefits(_p) if (_p != NULL) {casu_fits_delete(_p); _p = NULL;}
#define freefitslist(_p,_n) if (_p != NULL) {casu_fits_delete_list(_p,_n); _p = NULL;}
#define freetfitslist(_p,_n) if (_p != NULL) {casu_tfits_delete_list(_p,_n); _p = NULL;}
#define freetfits(_p) if (_p != NULL) {casu_tfits_delete(_p); _p = NULL;}
#define freemask(_p) if (_p != NULL) {casu_mask_delete(_p); _p = NULL;}

#define FATAL_ERR(_a) {freetable(tab); cpl_msg_error(fctid,_a); tidy(); return(CASU_FATAL);}
//check above for safety


typedef struct {
    cpl_image        *image;
    cpl_propertylist *phu;
    cpl_propertylist *ehu;
    char             *fname;
    char             *extname;
    char             *fullname;
    int              nexten;
    int              status;
    int              casufitstype;
    cpl_type         type;
} casu_fits;

typedef struct {
    cpl_table        *table;
    cpl_propertylist *phu;
    cpl_propertylist *ehu;
    char             *fname;
    char             *extname;
    char             *fullname;
    int              nexten;
    int              status;
} casu_tfits;

typedef struct {
    casu_tfits * catalogue; /* cpl_table and property list */
    cpl_image * segmentation_map;
    cpl_image * background;
} hdrl_imcore_result;


extern cpl_propertylist  *casu_fits_get_ehu(casu_fits *p);
extern cpl_propertylist *casu_tfits_get_ehu(casu_tfits *p);
extern cpl_image *casu_fits_get_image(casu_fits *p);
extern intptr_t casu_fits_get_nexten(casu_fits *p);
extern cpl_propertylist  *casu_fits_get_phu(casu_fits *p);
extern cpl_propertylist *casu_tfits_get_phu(casu_tfits *p);
extern int casu_tabwcs(cpl_propertylist *p, intptr_t xcol, intptr_t ycol,
			int *status);
    extern void casu_rename_property(cpl_propertylist *p, const char *oldname,char *newname);
    extern int casu_removewcs(cpl_propertylist *p, int *status);
extern cpl_table *casu_tfits_get_table(casu_tfits *p);
extern casu_fits *casu_fits_wrap(cpl_image *im, casu_fits *model, 
                                 cpl_propertylist *phu, cpl_propertylist *ehu);
extern casu_tfits *casu_tfits_wrap(cpl_table *tab, casu_tfits *model,
                                   cpl_propertylist *phu,
                                   cpl_propertylist *ehu);
extern void casu_fits_delete(casu_fits *p);
extern void casu_tfits_delete(casu_tfits *p);
extern int hdrl_casu_imcore(casu_fits *infile, casu_fits *conf, const cpl_wcs *,
                            intptr_t ipix,
                            float threshold, intptr_t icrowd, float rcore,
                            int bkg_subtr, intptr_t nbsize,
                            hdrl_catalogue_options cattype, float filtfwhm, float gain,
                            float saturation,
                            hdrl_imcore_result * res, int *status);
    extern casu_fits *casu_fits_duplicate(casu_fits *in);
    //extern int imcore_conf(casu_fits *infile, casu_fits *conf, int ipix, 
                       //float threshold, int icrowd, float rcore, int nbsize, 
                       //int cattype, float filtfwhm, float gain,
                       //casu_tfits **outcat);
    extern char *casu_fits_get_fullname(casu_fits *p);
    extern char *casu_fits_get_filename(casu_fits *p);
extern void casu_xytoradec(const cpl_wcs *wcs, double x, double y, double *ra,
                           double *dec);
#endif
