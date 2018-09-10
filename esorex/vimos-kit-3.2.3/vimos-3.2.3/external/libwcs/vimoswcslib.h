#ifndef vimoswcslib_h_
#define vimoswcslib_h_

/*=============================================================================
*
*   WCSLIB - an implementation of the FITS WCS proposal.
*   Copyright (C) 1995-1999, Mark Calabretta
*
*   This library is free software; you can redistribute it and/or modify it
*   under the terms of the GNU Library General Public License as published
*   by the Free Software Foundation; either version 2 of the License, or (at
*   your option) any later version.
*
*   This library is distributed in the hope that it will be useful, but
*   WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library
*   General Public License for more details.
*
*   You should have received a copy of the GNU Library General Public License
*   along with this library; if not, write to the Free Software Foundation,
*   Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*   Correspondence concerning WCSLIB may be directed to:
*      Internet email: mcalabre@atnf.csiro.au
*      Postal address: Dr. Mark Calabretta,
*                      Australia Telescope National Facility,
*                      P.O. Box 76,
*                      Epping, NSW, 2121,
*                      AUSTRALIA
*
*   Author: Mark Calabretta, Australia Telescope National Facility
*   $Id: vimoswcslib.h,v 1.1.1.1 2008-10-21 09:10:12 cizzo Exp $
*===========================================================================*/

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__STDC__) && !defined(__cplusplus)
#ifndef const
#define const
#endif
#endif

struct prjprm {
   int flag;
   int n;
   double r0;
   double p[10];
   double w[10];
};

#if __STDC__ || defined(__cplusplus)
   int vimosazpset(struct prjprm *);
   int azpfwd(const double, const double, struct prjprm *, double *, double *);
   int azprev(const double, const double, struct prjprm *, double *, double *);
   int vimostanset(struct prjprm *);
   int tanfwd(const double, const double, struct prjprm *, double *, double *);
   int tanrev(const double, const double, struct prjprm *, double *, double *);
   int vimossinset(struct prjprm *);
   int sinfwd(const double, const double, struct prjprm *, double *, double *);
   int sinrev(const double, const double, struct prjprm *, double *, double *);
   int vimosstgset(struct prjprm *);
   int stgfwd(const double, const double, struct prjprm *, double *, double *);
   int stgrev(const double, const double, struct prjprm *, double *, double *);
   int vimosarcset(struct prjprm *);
   int arcfwd(const double, const double, struct prjprm *, double *, double *);
   int arcrev(const double, const double, struct prjprm *, double *, double *);
   int vimoszpnset(struct prjprm *);
   int zpnfwd(const double, const double, struct prjprm *, double *, double *);
   int zpnrev(const double, const double, struct prjprm *, double *, double *);
   int vimoszeaset(struct prjprm *);
   int zeafwd(const double, const double, struct prjprm *, double *, double *);
   int zearev(const double, const double, struct prjprm *, double *, double *);
   int vimosairset(struct prjprm *);
   int airfwd(const double, const double, struct prjprm *, double *, double *);
   int airrev(const double, const double, struct prjprm *, double *, double *);
   int vimoscypset(struct prjprm *);
   int cypfwd(const double, const double, struct prjprm *, double *, double *);
   int cyprev(const double, const double, struct prjprm *, double *, double *);
   int vimoscarset(struct prjprm *);
   int carfwd(const double, const double, struct prjprm *, double *, double *);
   int carrev(const double, const double, struct prjprm *, double *, double *);
   int vimosmerset(struct prjprm *);
   int merfwd(const double, const double, struct prjprm *, double *, double *);
   int merrev(const double, const double, struct prjprm *, double *, double *);
   int vimosceaset(struct prjprm *);
   int ceafwd(const double, const double, struct prjprm *, double *, double *);
   int cearev(const double, const double, struct prjprm *, double *, double *);
   int vimoscopset(struct prjprm *);
   int copfwd(const double, const double, struct prjprm *, double *, double *);
   int coprev(const double, const double, struct prjprm *, double *, double *);
   int vimoscodset(struct prjprm *);
   int codfwd(const double, const double, struct prjprm *, double *, double *);
   int codrev(const double, const double, struct prjprm *, double *, double *);
   int vimoscoeset(struct prjprm *);
   int coefwd(const double, const double, struct prjprm *, double *, double *);
   int coerev(const double, const double, struct prjprm *, double *, double *);
   int vimoscooset(struct prjprm *);
   int coofwd(const double, const double, struct prjprm *, double *, double *);
   int coorev(const double, const double, struct prjprm *, double *, double *);
   int vimosbonset(struct prjprm *);
   int bonfwd(const double, const double, struct prjprm *, double *, double *);
   int bonrev(const double, const double, struct prjprm *, double *, double *);
   int vimospcoset(struct prjprm *);
   int pcofwd(const double, const double, struct prjprm *, double *, double *);
   int pcorev(const double, const double, struct prjprm *, double *, double *);
   int glsset(struct prjprm *);
   int glsfwd(const double, const double, struct prjprm *, double *, double *);
   int glsrev(const double, const double, struct prjprm *, double *, double *);
   int vimosparset(struct prjprm *);
   int parfwd(const double, const double, struct prjprm *, double *, double *);
   int parrev(const double, const double, struct prjprm *, double *, double *);
   int vimosaitset(struct prjprm *);
   int aitfwd(const double, const double, struct prjprm *, double *, double *);
   int aitrev(const double, const double, struct prjprm *, double *, double *);
   int vimosmolset(struct prjprm *);
   int molfwd(const double, const double, struct prjprm *, double *, double *);
   int molrev(const double, const double, struct prjprm *, double *, double *);
   int vimoscscset(struct prjprm *);
   int cscfwd(const double, const double, struct prjprm *, double *, double *);
   int cscrev(const double, const double, struct prjprm *, double *, double *);
   int vimosqscset(struct prjprm *);
   int qscfwd(const double, const double, struct prjprm *, double *, double *);
   int qscrev(const double, const double, struct prjprm *, double *, double *);
   int vimostscset(struct prjprm *);
   int tscfwd(const double, const double, struct prjprm *, double *, double *);
   int tscrev(const double, const double, struct prjprm *, double *, double *);
#else
   int vimosazpset(), azpfwd(), azprev();
   int vimostanset(), tanfwd(), tanrev();
   int vimossinset(), sinfwd(), sinrev();
   int vimosstgset(), stgfwd(), stgrev();
   int vimosarcset(), arcfwd(), arcrev();
   int vimoszpnset(), zpnfwd(), zpnrev();
   int vimoszeaset(), zeafwd(), zearev();
   int vimosairset(), airfwd(), airrev();
   int vimoscypset(), cypfwd(), cyprev();
   int vimoscarset(), carfwd(), carrev();
   int vimosmerset(), merfwd(), merrev();
   int vimosceaset(), ceafwd(), cearev();
   int vimoscopset(), copfwd(), coprev();
   int vimoscodset(), codfwd(), codrev();
   int vimoscoeset(), coefwd(), coerev();
   int vimoscooset(), coofwd(), coorev();
   int vimosbonset(), bonfwd(), bonrev();
   int vimospcoset(), pcofwd(), pcorev();
   int glsset(), glsfwd(), glsrev();
   int vimosparset(), parfwd(), parrev();
   int vimosaitset(), aitfwd(), aitrev();
   int vimosmolset(), molfwd(), molrev();
   int vimoscscset(), cscfwd(), cscrev();
   int vimosqscset(), qscfwd(), qscrev();
   int vimostscset(), tscfwd(), tscrev();
#endif

extern const char *prjset_errmsg[];
extern const char *prjfwd_errmsg[];
extern const char *prjrev_errmsg[];

#define PRJSET 137


extern int npcode;
extern char pcodes[25][4];

struct celprm {
   int flag;
   double ref[4];
   double euler[5];

#if __STDC__  || defined(__cplusplus)
   int (*prjfwd)(const double, const double,
                 struct prjprm *,
                 double *, double *);
   int (*prjrev)(const double, const double,
                 struct prjprm *,
                 double *, double *);
#else
   int (*prjfwd)();
   int (*prjrev)();
#endif
};

#if __STDC__  || defined(__cplusplus)
   int vimoscelset(const char *, struct celprm *, struct prjprm *);
   int celfwd(const char *,
              const double, const double,
              struct celprm *,
              double *, double *,
              struct prjprm *,
              double *, double *);
   int celrev(const char *,
              const double, const double,
              struct prjprm *,
              double *, double *,
              struct celprm *,
              double *, double *);
#else
   int vimoscelset(), celfwd(), celrev();
#endif

extern const char *vimoscelset_errmsg[];
extern const char *celfwd_errmsg[];
extern const char *celrev_errmsg[];

#define CELSET 137

struct linprm {
   int flag;
   int naxis;
   double *crpix;
   double *pc;
   double *cdelt;

   /* Intermediates. */
   double *piximg;
   double *imgpix;
};

#if __STDC__  || defined(__cplusplus)
   int vimoslinset(struct linprm *);
   int linfwd(const double[], struct linprm *, double[]);
   int linrev(const double[], struct linprm *, double[]);
   int vimosmatinv(const int, const double [], double []);
#else
   int vimoslinset(), linfwd(), linrev(), vimosmatinv();
#endif

extern const char *vimoslinset_errmsg[];
extern const char *linfwd_errmsg[];
extern const char *linrev_errmsg[];

#define LINSET 137


struct vimoswcsprm {
   int flag;
   char pcode[4];
   char lngtyp[5], lattyp[5];
   int lng, lat;
   int cubeface;
};

#if __STDC__ || defined(__cplusplus)
   int vimoswcsset(const int,
              const char[][9],
              struct vimoswcsprm *);

   int vimoswcsfwd(const char[][9],
              struct vimoswcsprm *,
              const double[],
              const double[],
              struct celprm *, 
              double *,
              double *, 
              struct prjprm *, 
              double[], 
              struct linprm *,
              double[]);

   int vimoswcsrev(const char[][9],
              struct vimoswcsprm *,
              const double[], 
              struct linprm *,
              double[], 
              struct prjprm *, 
              double *,
              double *, 
              const double[], 
              struct celprm *, 
              double[]);

   int vimoswcsmix(const char[][9],
              struct vimoswcsprm *,
              const int,
              const int,
              const double[],
              const double,
              int,
              double[],
              const double[],
              struct celprm *,
              double *,
              double *,
              struct prjprm *,
              double[], 
              struct linprm *,
              double[]);

#else
   int vimoswcsset(), vimoswcsfwd(), vimoswcsrev(), vimoswcsmix();
#endif

extern const char *vimoswcsset_errmsg[];
extern const char *vimoswcsfwd_errmsg[];
extern const char *vimoswcsrev_errmsg[];
extern const char *vimoswcsmix_errmsg[];

#define VIMOSWCSSET 137


#if __STDC__  || defined(__cplusplus)
   int sphfwd(const double, const double,
              const double [],
              double *, double *);
   int sphrev(const double, const double,
              const double [],
              double *, double *);
#else
   int sphfwd(), sphrev();
#endif

#ifdef PI
#undef PI
#endif

#ifdef D2R
#undef D2R
#endif

#ifdef R2D
#undef R2D
#endif

#ifdef SQRT2
#undef SQRT2
#endif

#ifdef SQRT2INV
#undef SQRT2INV
#endif

#define PI 3.141592653589793238462643
#define D2R PI/180.0
#define R2D 180.0/PI
#define SQRT2 1.4142135623730950488
#define SQRT2INV 1.0/SQRT2

#if !defined(__STDC__) && !defined(__cplusplus)
#ifndef const
#define const
#endif
#endif

#if __STDC__ || defined(__cplusplus)
   double cosdeg(const double);
   double sindeg(const double);
   double tandeg(const double);
   double acosdeg(const double);
   double asindeg(const double);
   double atandeg(const double);
   double atan2deg(const double, const double);
#else
   double cosdeg();
   double sindeg();
   double tandeg();
   double acosdeg();
   double asindeg();
   double atandeg();
   double atan2deg();
#endif

/* Domain tolerance for asin and acos functions. */
#define VIMOSWCSTRIG_TOL 1e-10

#ifdef __cplusplus
};
#endif

#endif /* vimoswcslib_h_ */

/* Feb  3 2000	Doug Mink - Make cplusplus ifdefs for braces all-inclusive
 *
 * Feb 15 2001	Doug Mink - Undefine math constants if already defined
 */
