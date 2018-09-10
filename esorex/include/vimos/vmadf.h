/* $Id: vmadf.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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
 * $Date: 2013-03-25 11:43:04 $
 * $Revision: 1.2 $
 * $Name: not supported by cvs2svn $
 */

#ifndef VM_ADF_H
#define VM_ADF_H

#include <pilmacros.h>

#include <vmimage.h>
#include <vmmatrix.h>
#include <vmdistmodels.h>


PIL_BEGIN_DECLS

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   defines of ADF table names

   Description:
   Allowed names for VimosTables that contain an ADF

   Values:
     VM_ADF_MOS     MOS ADF
     VM_ADF_IFU     IFU ADF
     VM_ADF_IMA     Imaging ADF

   Updates:
   11 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
#define    VM_ADF_MOS     "ADF MOS"
#define    VM_ADF_IFU     "ADF IFU"
#define    VM_ADF_IMA     "ADF IMA"

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   defines of ADF type names

   Description:
   Allowed string for the image keyword ESO INS ADF TYPE

   Values:
     ADF_TYPE_MOS     MOS
     ADF_TYPE_IFU     IFU 
     ADF_TYPE_IMA     IMAGE

   Updates:
   21 Dec 99: Created (MS)

--------------------------------------------------------------------------------
*/
#define    ADF_TYPE_MOS     "MOS"
#define    ADF_TYPE_IFU     "IFU"
#define    ADF_TYPE_IMA     "IMAGE"


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   enum VimosAdfType

   Description:
   Enumeration of type of possible ADFs

   Values:
    VM_ADF_TYPE_UNDF     ADF type undefined
    VM_ADF_TYPE_MOS      MOS ADF
    VM_ADF_TYPE_IFU      IFU ADF
    VM_ADF_TYPE_IMA      Imaging ADF

   Updates:
   11 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
typedef enum _VIMOS_ADF_TYPE_ 
{
  VM_ADF_TYPE_UDF = 0,
  VM_ADF_TYPE_MOS = 1,
  VM_ADF_TYPE_IFU = 2,
  VM_ADF_TYPE_IMA = 3
} VimosAdfType;


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   define of IFU Mode names

   Description: 
   Allowed strings in the image header keyword that describes the IFU Mode

   Values:
     IFU_MODE_SMALL   SMALL
     IFU_MODE_LARGE   LARGE

   Updates:
   13 Dec 99: Created (MS)

--------------------------------------------------------------------------------
*/
#define    IFU_MODE_SMALL     "SMALL"
#define    IFU_MODE_LARGE     "LARGE"

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   enum VimosIfuMode

   Description:
   Enumeration of type of possible IFU Modes

   Values:
    VM_IFU_MODE_UNDEF     IFU Mode undefined
    VM_IFU_MODE_SMALL     IFU small field (only IFU slit 2 fibers)
    VM_IFU_MODE_LARGE     IFU large field (all IFU slits fibers)

   Updates:
   13 Dec 99: Created (MS)

--------------------------------------------------------------------------------
*/
typedef enum _VIMOS_IFU_MODE_ 
{
  VM_IFU_MODE_UNDEF = 0,
  VM_IFU_MODE_SMALL = 1,
  VM_IFU_MODE_LARGE = 2
} VimosIfuMode;


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  enum VimosIfuMode

  Description:
  Enumeration of type of possible IFU Modes

  Values:
  VM_MSHU_MODE_UNDEF     Mask Shutter Mode undefined
  VM_MSHU_MODE_ON        Mask Shutter Mode ON (mask shutters are being used)
  VM_MSHU_MODE_OFF       Mask Shutter Mode OFF (mask shutters are NOT being used)
  
  Updates:
  18 Apr 02: Created (MS)

--------------------------------------------------------------------------------
*/
typedef enum _VIMOS_MSHU_MODE_ 
{
  VM_MSHU_MODE_UNDEF = 0,
  VM_MSHU_MODE_ON = 1,
  VM_MSHU_MODE_OFF = 2
} VimosMshuMode;


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   defines of MOS ADF name for slit types

   Description:
   Allowed names for slit types in a MOS ADF

   Values:
     VM_ADF_RECT_SLIT_NAME     For rectangular slit
     VM_ADF_CURV_SLIT_NAME     For curved slit

   Updates:
   11 Nov 98: Created (TAO)
   30 Nov 99: Changed the curved slit name into CURVE (MS)
--------------------------------------------------------------------------------
*/
#define    VM_ADF_RECT_SLIT_NAME     "RECTANGLE"
#define    VM_ADF_CURV_SLIT_NAME     "CURVE"


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   enum VimosAdfRectSlit

   Description:
   Enumeration of slit types that can occur in a MOS ADF

   Values:
     VM_ADF_UNDF_SLIT     Slit type undefined
     VM_ADF_RECT_SLIT     Rectangular slit
     VM_ADF_CURV_SLIT     Curved slit
     VM_ADF_CIRC_SLIT     Circular aperture
     VM_ADF_REFR_SLIT     Reference aperture (assumed to be square)

   Updates:
   17 Nov 98: Added circular and reference apertures (TAO)
   11 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
typedef enum _VIMOS_ADF_SLIT_TYPE_ 
{
  VM_ADF_UNDF_SLIT = 0,
  VM_ADF_RECT_SLIT = 1,
  VM_ADF_CURV_SLIT = 2,
  VM_ADF_CIRC_SLIT = 3,
  VM_ADF_REFR_SLIT = 4
} VimosAdfSlitType;




/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  structure VimosBezierCurve

   Description:
   Structure containing the parameters of Bezier curve.
   A Bezier curve is:
     x(t) = a t^3 + b t^2 + c t + x_0
   with t running from 0 to 1.

   Layout:
     float x0          Bezier coefficient
     float a           Bezier coefficient
     float b           Bezier coefficient
     float c           Bezier coefficient

   Updates:
   19 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_BEZIER_CURVE_ 
{
  float x0;
  float a;
  float b;
  float c;
} VimosBezierCurve;



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   structure VimosAdfRectSlit

   Description:
   Structure containing the parameters of a slit as defined for a rectangular
   ADF slit

   Layout:
     VimosAdfSlitType  slitType  Type of slit: VM_ADF_RECT_SLIT 
     int               slitNo    Sequential number of slit in ADF
     float             x         X mask coordinate (mm) of slit centre 
     float             y         Y mask coordinate (mm) of slit centre 
     float             dimX      Size of slit in spatial direction (mm) 
     float             dimY      Size of slit in dispersion direction (mm) 

   Updates:
   11 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_ADF_RECT_SLIT_
{
  VimosAdfSlitType    slitType; /* Type of slit: VM_ADF_RECT_SLIT            */
  int                 slitNo;   /* Sequential number of slit in ADF          */
  float               x;        /* X mask coordinate (mm) of slit centre     */
  float               y;        /* Y mask coordinate (mm) of slit centre     */
  float               dimX;     /* Size of slit in spatial direction (mm)    */
  float               dimY;     /* Size of slit in dispersion direction (mm) */
} VimosAdfRectSlit;



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   structure VimosAdfCurvSlit

   Description:
   Structure containing the parameters of a slit as defined for a curved
   ADF slit. A curved slit is parametrized by a Bezier curve:
     x(t) = a_x t^3 + b_x t^2 + c_x t + x_0
     y(t) = a_y t^3 + b_y t^2 + c_y t + y_0
   where the parameter t runs from 0 to 1

   Layout:
   VimosAdfSlitType    slitType   Type of slit: VM_ADF_TYPE_CURVE    
   int                 slitNo     Sequential number of slit in ADF       
   float               deltay     Separation between the 2 edges along Y axis
   VimosBezierCurve    *xMiddle;
   VimosBezierCurve    *yMiddle;
   
   Updates:
   11 Nov 98: Created (TAO)
   30 Nov 99: Updated the description of the slit, to follow the new
              header keyword structure. Only 1 curve is given (the one
	      describing the middle line of the slit), plus the separation
	      between the two edges along the Y axis. (MS)
--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_ADF_CURV_SLIT_
{
  VimosAdfSlitType    slitType; /* Type of slit: VM_ADF_TYPE_CURVE        */
  int                 slitNo;   /* Sequential number of slit in ADF       */
  float               deltay;
  VimosBezierCurve    *xMiddle;
  VimosBezierCurve    *yMiddle;
} VimosAdfCurvSlit;


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   structure VimosAdfCircSlit

   Description:
   Structure containing the parameters of a slit as defined for a circular
   ADF slit (e.g. for IFU)

   Layout:
     VimosAdfSlitType  slitType  Type of slit: VM_ADF_CIRC_SLIT 
     int               slitNo    Sequential number of slit in ADF
     float             x         X mask coordinate (mm) of slit centre 
     float             y         Y mask coordinate (mm) of slit centre 
     float             radius    Radius of aperture  (mm) 
     int               IFUslitNo IFU Slit Number
     int               IFUfibNo  Ifu Fiber number within Slit
     float             IFUfibTrans  Ifu Fiber Relative Transmission

   Updates:
   17 Nov 98: Created (TAO)
   14 Dec 99: Included IFU slit and fiber info (MS)
   25 Jul 02: Added IFUfibTrans (AZ)

--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_ADF_CIRC_SLIT_
{
  VimosAdfSlitType    slitType;  /* Type of slit: VM_ADF_CIRC_SLIT           */
  int                 slitNo;    /* Sequential number of slit in ADF         */
  float               x;         /* X mask coordinate (mm) of slit centre    */
  float               y;         /* Y mask coordinate (mm) of slit centre    */
  float               radius;    /* Radius of aperture (mm)                 */
  int                 IFUslitNo; /* IFU Slit number */
  int                 IFUfibNo;  /* IFU Fiber number within IFU Slit */
  float               IFUfibTrans;  /* IFU Fiber Transmission */
} VimosAdfCircSlit;


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   structure VimosAdfRefrSlit

   Description:
   Structure containing the parameters of a slit as defined for reference
   aperture in a mask (e.g. for mask alignment). A reference aperture is
   assumed to be a square aperture.

   Layout:
     VimosAdfSlitType  slitType  Type of slit: VM_ADF_REFR_SLIT 
     int               slitNo    Sequential number of slit in ADF
     float             x         X mask coordinate (mm) of slit centre 
     float             y         Y mask coordinate (mm) of slit centre 
     float             size      Size of aperture  (mm) 
     double            objRA     RA of reference object (J2000)       
     double            objDec    Dec of reference object (J2000)  

   Updates:
   17 Nov 98: Created (TAO)
   10 Jan 00: Modified RA and Dec to the the reference object coords (MS)

--------------------------------------------------------------------------------
*/
typedef struct _VIMOS_ADF_REFR_SLIT_
{
  VimosAdfSlitType    slitType; /* Type of slit: VM_ADF_REFR_SLIT         */
  int                 slitNo;   /* Sequential number of slit in ADF       */
  float               x;        /* X mask coordinate (mm) of slit centre  */
  float               y;        /* Y mask coordinate (mm) of slit centre  */
  float               size;     /* Size of aperture (mm)                  */
  double              objRA;    /* RA of reference object (J2000)       */
  double              objDec;   /* Dec of reference object (J2000)      */
} VimosAdfRefrSlit;


typedef struct _VIMOS_ADF_SLIT_HOLDER_ 
{
  VimosAdfSlitType               slitType;
  int                            slitNo;
  void                           *slit; /* pointer to one of the slit objects, 
                                           should be cast */
  struct _VIMOS_ADF_SLIT_HOLDER_ *prev;
  struct _VIMOS_ADF_SLIT_HOLDER_ *next;
} VimosAdfSlitHolder;



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBezierCurve  *newBezierCurve()

   Description: 
   Allocates a new VimosBezierCurve structure
   
   Input:
   void   
   
   Return Value (succes):
   Pointer to a newly allocated Bezier structure. All fields are
   initialized to 0.

   Return Value (error):
   NULL
   
   Updates:
   14 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   19 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosBezierCurve *newBezierCurve();



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void  deleteBezierCurve(VimosBezierCurve *curve)

   Description: 
   Deletes a VimosBezierCurve
   
   Input:
   VimosBezierCurve *curve
   pointer to VimosBezierCurve to delete
   
   Return Value:
   void
   
   Updates:
   19 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteBezierCurve(VimosBezierCurve *curve);





/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   float   computeBezierCurve(VimosBezierCurve *curve, float t)

   Description: 
   Computes the value of a Bezier Curve. If t is less than 0, 0 is use in
   computing the value, if t larger than 1, 1 is used. 
   
   Input:
   VimosBezierCurve *curve
   pointer to VimosBezierCurve to use in computation
   
   Return Value (succes):
   Value of Bezier curve
   
   Return Value (error):
   -1.0

   Updates:
   19 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
float   computeBezierCurve(VimosBezierCurve *curve, float t);




/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosAdfRectSlit  *newAdfRectSlit()

   Description: 
   Allocates a new VimosAdfRectSlit structure
   
   Input:
   void   
   
   Return Value (succes):
   Pointer to a newly allocated VimosAdfRectSlit structure. All fields are
   initialized to 0, except the slitType field is set to VM_ADF_RECT_SLIT.

   Return Value (error):
   NULL
     
   Updates:
   14 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   11 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosAdfRectSlit  *newAdfRectSlit();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteAdfRectSlit(VimosAdfRectSlit *adfSlit)
   
   Description: 

   Deletes a VimosAdfRectSlit
   
   Input:
   VimosAdfRectSlit *adfSlit
   Pointer to VimosAdfRectSlit to be deleted
   
   Return Value:
   void
   
   Updates:
   11 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteAdfRectSlit(VimosAdfRectSlit *adfSlit);




/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosAdfCurvSlit  *newAdfCurvSlit()

   Description: 
   Allocates  a new VimosAdfCurvSlit structure.
   
   Input:
   void   
   
   Return Value (succes):
   Pointer to a newly allocated VimosAdfCurvSlit structure. All fields are
   initialized to 0, except the slitType field is set to VM_ADF_CURV_SLIT.

   Return Value (error):
   NULL
      
   Updates:
   14 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   20 Nov 98: Adopted for new layout with VimosBezierCurves
   11 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosAdfCurvSlit  *newAdfCurvSlit();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteAdfCurvSlit(VimosAdfCurvSlit *adfSlit)
   
   Description: 

   Deletes a VimosAdfCurvSlit
   
   Input:
   VimosAdfCurvSlit *adfSlit
   Pointer to VimosAdfCurvSlit to be deleted
   
   Return Value:
   void
   
   Updates:
   20 Nov 98: Adopted for new layout with VimosBezierCurves
   11 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteAdfCurvSlit(VimosAdfCurvSlit *adfSlit);




/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosAdfCircSlit  *newAdfCircSlit()

   Description: 
   Allocates a new VimosAdfCircSlit structure
   
   Input:
   void   
   
   Return Value (succes):
   Pointer to a newly allocated VimosAdfCircSlit structure. All fields are
   initialized to 0, except the slitType field is set to VM_ADF_CIRC_SLIT.

   Return Value (error):
   NULL

   
   Updates:
   14 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   17 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosAdfCircSlit  *newAdfCircSlit();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteAdfCircSlit(VimosAdfCircSlit *adfSlit)
   
   Description: 

   Deletes a VimosAdfCircSlit
   
   Input:
   VimosAdfCircSlit *adfSlit
   Pointer to VimosAdfCircSlit to be deleted
   
   Return Value:
   void
   
   Updates:
   17 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteAdfCircSlit(VimosAdfCircSlit *adfSlit);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosAdfRefrSlit  *newAdfRefrSlit()

   Description: 
   Allocates a new VimosAdfRefrSlit structure
   
   Input:
   void   
   
   Return Value (succes):
   Pointer to a newly allocated VimosAdfRefrSlit structure. All fields are
   initialized to 0, except the slitType field is set to VM_ADF_REFR_SLIT.

   Return Value (error):
   NULL

   
   Updates:
   14 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   17 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosAdfRefrSlit  *newAdfRefrSlit();


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void deleteAdfRefrSlit(VimosAdfRefrSlit *adfSlit)
   
   Description: 

   Deletes a VimosAdfRefrSlit
   
   Input:
   VimosAdfRefrSlit *adfSlit
   Pointer to VimosAdfRefrSlit to be deleted
   
   Return Value:
   void
   
   Updates:
   17 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteAdfRefrSlit(VimosAdfRefrSlit *adfSlit);


VimosAdfSlitHolder *newAdfSlitHolder();
void deleteAdfSlitHolder(VimosAdfSlitHolder *holder);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosTable *newADF()

   Description:
   Allocate a new ADF. This is a generic VimosTable with no name.
   No descriptors or columns are allocated (but this may change in a later
   version).
   
   Input:
   none
   
   Return Value (success):
   Pointer to the new VimosTable

   Return Value (error):
   NULL

   
   Updates:
   14 Jun 00: Return NULL on error instead of exiting (Maura)
   21 Dec 99: modified to have no input type (MS)
   24 Nov 98: call exit(-1) on error (TAO)
   11 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosTable  *newADF();

/*
 *  The following function is doing exactly what readFitsADF() is doing,
 *  but avoiding the dependency from the FITS file. This was in fact
 *  unnecessary, as the descriptors are already contained in the input
 *  VimosImage, and do not need to be read again from disk. This
 *  routine, as well as readFitsADF, need to be improved in the sense
 *  that the ADF table should consist of just a subset of all the
 *  descriptors found in header, that is just the descriptor related
 *  to ADF.
 */

VimosBool readADF(VimosTable *adf, VimosImage *adfImage);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   void deleteADF(VimosTable *adf)
   
   Description:
   Delete an ADF. This is just an esthetic wrapper for deleteTable(eTable)
   
   Input: 
   VimosTable *adf
   Pointer of ADF to be deleted
   
   Return Value:
   void
   
   Updates:
   11 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- 
*/
void deleteADF(VimosTable *adf);


    
/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   vimosAdfType getADFTypeFromDesc(VimosDescriptor *desc)
   
   Description:
   Determine which type of ADF is defined in the Descriptors
   
   Input: 
   VimosDescriptor *desc
   Descriptor to search
   
   Return Value:
   VimosAdfType
   Type code of the type of ADF defined. If no correct ADF type is defined,
   VM_ADF_TYPE_UDF is returned
   
   Updates:
   22 Mar 99: Created (TAO)

-------------------------------------------------------------------------------- 
*/
VimosAdfType getADFTypeFromDesc(VimosDescriptor *desc);




/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosAdfSlitHolder *extractSlitsFromADF(VimosTable *adf)

   Description:
   
   Input:
   VimosTable *adf
   Pointer to the ADF Table to read from. 
   
   Return Value (success):

   Return Value (error):
   NULL
   
   Updates:
   14 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   20 Nov 98: Changed signature: output is now a VimosAdfSlitHolder and is
              is the return value. Adapted implementation to use
              VimosAdfSlitHolder (TAO)
   19 Nov 98: Added typeAr array, changed signature (TAO)
              Accommodated new layout of curved slits (TAO)
   17 Nov 98: Accommodated circular apertures for IFU and IMA
   11 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- 
*/
VimosAdfSlitHolder *extractSlitsFromADF(VimosTable *adf);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosAdfSlitHolder *extractRefsFromADF(VimosTable *adf)

   Description:

   Input:
   VimosTable *adf
   Pointer to the ADF Table to read from. Should be a MOS ADF
   
   Return Value (success):

   Return Value (error):
   NULL
      
   Updates:
   14 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: call exit(-1) on error (TAO)
   20 Nov 98: Changed signature: output is now a VimosAdfSlitHolder and is
              is the return value. Adapted implementation to use
              VimosAdfSlitHolder (TAO)
   19 Nov 98: Changed signature to stay conform extractSlitsFromAdf (TAO)
   17 Nov 98: Created (TAO)

-------------------------------------------------------------------------------- 
*/
VimosAdfSlitHolder *extractRefsFromADF(VimosTable *adf);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  20 Dec 99: added return of IFU slit and fiber number for circular slit (MS)
  22 Dec 99: added the handling of the contamination model (MS)
  23 Oct 01: added IFUfibPeakX parameter in calcSlitLocationsOnCCD (AZ)
  23 May 02: added rms values for the curvature polynomial and inverse
             dispersion relation fit (MS)
  25 Jul 02: added IFUfibTrans parameter in calcSlitLocationsOnCCD (AZ)
-------------------------------------------------------------------------------- 
*/

VimosBool calcSlitLocationsOnCCD(void *slit, VimosAdfSlitType slitType,
                          VimosDistModel2D *optModX,VimosDistModel2D *optModY,
                          VimosDistModelFull *crvMod, 
                          VimosDistModelFull *invDispMat,
			  VimosDistModel2D *contModX, 
			  VimosDistModel2D *contModY,
                          VimosFloatArray **ccdX, VimosFloatArray **ccdY, 
                          VimosFloatArray **maskX, VimosFloatArray **maskY, 
 			  VimosDistModel1D ***crvPol, 
			  VimosFloatArray **crvPolRms, 
			  VimosDistModel1D ***invDis, 
			  VimosFloatArray **invDisRms, float lamda0,
                          int *numPix, int *IFUslitNo, int *IFUfibNo,
			  float *IFUfibPeakX, float *IFUfibTrans,
			  VimosFloatArray **zeroX, VimosFloatArray **zeroY,
			  VimosIntArray **invDisQuality);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   VimosBool readFitsADF(VimosTable *adf, VimosImage *adfImage)

   Description:
   Read the ADF from the header of the FITS Image. The type of the ADF should 
   match that defined in the image.
   
   Input:
   VimosTable *adf  
   Pointer to the ADF to copy the FITS header into.

   VimosImage *adfImage
   Structure holding the image 

   Return Value:
   VimosBool (true if read is OK, otherwise false).
   
   Updates:
   15 Dec 99: Created (BG)
   03 Apr 00: modified input, so that the image is not opened twice (MS)
--------------------------------------------------------------------------------
*/
VimosBool readFitsADF(VimosTable *adf, VimosImage *adfImage);

PIL_END_DECLS

#endif /* VM_ADF_H */
