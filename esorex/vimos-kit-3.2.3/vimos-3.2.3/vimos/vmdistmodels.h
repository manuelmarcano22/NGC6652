/* $Id: vmdistmodels.h,v 1.2 2013-03-25 11:43:04 cgarcia Exp $
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

#ifndef VM_DISTMODELS_H
#define VM_DISTMODELS_H

#include <pilmacros.h>

#include <vmtable.h>


PIL_BEGIN_DECLS

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   structure VimosDistModel1D

   Description:
   Structure containing the coefficients of a 1D distortion model:
     f(x) = sum a_i (x-offset)^i
   where i runs from 0 to order of the model

   Layout:
     int     order      Order of the model
     double   *coefs     Array of length order with the coefficients of the 
                        polynomial
     double   offset     Offset term in polynomial
   Updates:
   13 Nov 98: Created (TAO)
   18 Nov 98: Added offset term (TAO)

--------------------------------------------------------------------------------
*/
typedef struct _VM_DIST_MODEL_1D_ 
{
  int       order;    /* Order of the polynomial of the model          */
  double     *coefs;   /* Array of length order+1 with the coefficients 
                         of the polynomial                             */
  double     offset;   /* Offset term in polynomial                     */
} VimosDistModel1D;


/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   structure VimosDistModel2D

   Description:
   Structure containing the coefficients of a 2D distortion model:
     f(x,y) = sum a_ij (x-offsetX)^i (y-offsetY)^j
   where i runs from 0 to orderX and j from 0 to orderY of the model

   Layout:
     int     orderX     Order of the model in X
     int     orderY     Order of the model in Y
     double   **coefs    Coefficients of the polynomial
     double   offsetX    X offset term in polynomial
     double   offsetY    Y offset term in polynomial

   Updates:
   13 Nov 98: Created (TAO)
   18 Nov 98: Offset term added (TAO)

--------------------------------------------------------------------------------
*/
typedef struct _VM_DIST_MODEL_2D_ 
{
  int       orderX;   /* Order of the polynomial of the model in X      */
  int       orderY;   /* Order of the polynomial of the model in Y      */
  double     **coefs;  /* Matrix with the coefficients of the polynomial */
  double     offsetX;  /* X offset term in polynomial                    */
  double     offsetY;  /* Y offset term in polynomial                    */

} VimosDistModel2D;




/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   structure VimosDistModelFull

   Description:
   Structure containing the coefficients of a full distortion model:
     f(t; x,y) = sum A_i (t-offsetT)^i
   with
     A_i = sum a_i,jk (x-offsetX)^i (y-offsetY)^j.
   i runs from 0 to orderPol,  j from 0 to orderX and k from 0 to
   order Y of the model. The models a_i,jk are stored as an array of
   VimosDistModel2D

   Layout:
     int                orderPol  Order of the model polynomial (range of i)
     int                orderX    Order of the model in X (range of j)
     int                orderY    Order of the model in Y (range of k)
     VimosDistModel2D   *coefs    Models for the coefficients a_i,jk
     double              offsetT   t offset term in polynomial
     double              offsetX   X offset term in polynomial
     double              offsetY   Y offset term in polynomial

   Updates:
   13 Nov 98: Created (TAO)
   18 Nov 98: Offset terms added (TAO)

--------------------------------------------------------------------------------
*/
typedef struct _VM_DIST_MODEL_FULL_ 
{
  int               orderPol; 
  int               orderX;   
  int               orderY;   
  VimosDistModel2D  **coefs;   
  double             offsetT;
  double             offsetX;
  double             offsetY;
  
} VimosDistModelFull;



/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   structure VimosGnomonic

   Description:
   Structure containing the coefficients of a Gnomonic projection

   Layout:
    double         alpha0    RA of projection centre in radians!!
    double         delta0    Dec of projection centre in radians!!
    double         sina0     sin(alpha0)
    double         cosa0     cos(alpha0)
    double         sind0     sin(delta0)
    double         cosd0     cos(delta0)

   Updates:
   25 Feb 99: Created (TAO)


--------------------------------------------------------------------------------
*/
typedef struct _VM_GNOMONIC_ 
{
  double         alpha0;
  double         delta0;
  double         sina0;
  double         cosa0;
  double         sind0;
  double         cosd0;
} VimosGnomonic;



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosDistModel1D    *newDistModel1D(const int order)

   Description: 
   Allocates a new VimosDistModel1D of order order
   
   Input:
   const int   order     order of polynomial of the model

   Return Value (succes):
   Pointer to a newly allocated VimosDistModel1D structure. All coefficients 
   are initialized to 0.

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: Call exit(-1) on error
   18 Nov 98: Offset term added (TAO)
   13 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosDistModel1D   *newDistModel1D(const int order);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosDistModel2D    *newDistModel2D(const int orderX, const int orderY)

   Description: 
   Allocates a new VimosDistModel2D of order orderX in X and orderY in Y
   
   Input:
   const int   orderX    maximum order of X terms of the model
   const int   orderY    maximum order of Y terms of the model
   
   Return Value (succes):
   Pointer to a newly allocated VimosDistModel2D structure. All coefficients 
   are initialized to 0.

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: Call exit(-1) on error
   18 Nov 98: Offset terms added 
   13 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosDistModel2D   *newDistModel2D(const int orderX, const int orderY);


 
/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosDistModelFull  *newDistModelFull(const int orderPol, 
                         const int orderX, const int orderY)

   Description: 
   Allocates a new VimosDistModelFull of order orderPol in t, 
   orderX in X and orderY in Y
   
   Input:
   const int   orderPol  maximum order of t terms of the model
   const int   orderX    maximum order of X terms of the model
   const int   orderY    maximum order of Y terms of the model
   
   Return Value (succes):
   Pointer to a newly allocated VimosDistModelFull structure. All coefficients 
   are initialized to 0.

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead of exiting (Maura)
   24 Nov 98: Call exit(-1) on error
   18 Nov 98: Offset term sadded (TAO)
   13 Nov 98: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosDistModelFull *newDistModelFull(const int orderPol, 
                                     const int orderX, const int orderY);


 
/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void  deleteDistModel1D(VimosDistModel1D *model)
   void  deleteDistModel2D(VimosDistModel2D *model)
   void  deleteDistModelFull(VimosDistModelFull *model)

   Description: 
   Desctructors of distortion model structures
   
   Input:
   VimosDistModel1D   *model        
   Pointer of structure to be deleted
   VimosDistModel2D   *model
   Pointer of structure to be deleted
   VimosDistModelFull *model
   Pointer of structure to be deleted
   
   Return Value:
   void
   
   Updates:
   13 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteDistModel1D(VimosDistModel1D *mod);
void deleteDistModel2D(VimosDistModel2D *mod);
void deleteDistModelFull(VimosDistModelFull *mod);

/* compute model values */
double computeDistModel1D(VimosDistModel1D *mod, float x);
double computeDistModel2D(VimosDistModel2D *mod, float x, float y);
double computeDistModelFull(VimosDistModelFull *mod,float t,float x,float y);




/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBool getDistModel1DFromFull(VimosDistModelFull *full, float x, 
                                    float y, VimosDistModel1D **mod1D)

   Description: 
   Extract the a 1 Distortion Model that is modelled using a Full Distortion
   Model.
   A Full Distortion Model gives
     f(t; x,y) = sum A_i (t-offsetT)^i
   with
     A_i = sum a_i,jk (x-offsetX)^i (y-offsetY)^j,
   where the a_i,jk are the coefficients of the Full Model. This function
   computes the coefficients A_i so that f(t) can be computed using a 1D
   Distorion Model

   Input:
   VimosDistModelFull
   Full Distortion Model
   float x, y
   position to compute 1D Distortion Model for

   Output:
   VimosDistModel1D **mod1D
   the computed 1D Distortion Model

   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   
   Updates:
   24 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosBool getDistModel1DFromFull(VimosDistModelFull *full, float x, float y, 
                                 VimosDistModel1D **mod1D);
  
 
/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBool readOptDistModel(VimosDescriptor *desc, 
                              VimosDistModel2D **optModX,
                              VimosDistModel2D **optModY)

   Description: 
   Read the Optocal Distortion Models from the descriptors
   
   Input:
   VimosDescriptor *desc
   Linked list of descriptors that contani the coefficients of the Optical
   Distortion Models
   
   Output:
   VimosDistModel2D **optModX
   Contains the Optocal Distortion Model for X
   VimosDistModel2D **optModY
   Contains the Optocal Distortion Model for Y

   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   24 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosBool readOptDistModel(VimosDescriptor *desc, VimosDistModel2D **optModX,
                           VimosDistModel2D **optModY);
/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBool writeOptDistModel(VimosDescriptor **desc, 
                               VimosDistModel2D *optModX,
                               VimosDistModel2D *optModY)

   Description: 
   Write the Optical Distortion Models to the descriptors
   
   Input:
   VimosDistModel2D *optModX
   Contains the Optical Distortion Model for X
   VimosDistModel2D *optModY
   Contains the Optical Distortion Model for Y

   Output:   
   VimosDescriptor **desc
   Address of pointer to Linked list of descriptors to write to
   
   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   30 Mar 99: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosBool writeOptDistModel(VimosDescriptor **desc, VimosDistModel2D *optModX,
                            VimosDistModel2D *optModY);
VimosBool writeOptDistModelString(VimosDescriptor **, VimosDistModel2D *,
                                  VimosDistModel2D *);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBool readContaminationModel(VimosDescriptor *desc, 
                                    VimosDistModel2D **contModX, 
                                    VimosDistModel2D **contModY)

   Description: 
   Read the Contamination Models from the descriptors
   
   Input:
   VimosDescriptor *desc
   Linked list of descriptors that contain the coefficients of the 
   Contamination Models
   
   Output:
   VimosDistModel2D **contModX
   Contains the Contamination Model for X
   VimosDistModel2D **contModY
   Contains the Contamination Model for Y

   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   22 Dec 99: Created (MS)
   
--------------------------------------------------------------------------------
*/
VimosBool readContaminationModel(VimosDescriptor *desc, 
		 VimosDistModel2D **contModX, VimosDistModel2D **contModY);

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBool writeContaminationModel(VimosDescriptor **desc, 
                               VimosDistModel2D *contModX,
                               VimosDistModel2D *contModY)

   Description: 
   Write the Optical Distortion Models to the descriptors
   
   Input:
   VimosDistModel2D *contModX
   Contains the Contamination Model for X
   VimosDistModel2D *contModY
   Contains the Contamination Model for Y

   Output:   
   VimosDescriptor **desc
   Address of pointer to Linked list of descriptors to write to
   
   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   30 Nov 00: Created (MS)
   
--------------------------------------------------------------------------------
*/
VimosBool writeContaminationModel(VimosDescriptor **desc, VimosDistModel2D 
				  *contModX, VimosDistModel2D *contModY);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   
   VimosBool readCurvatureModel(VimosDescriptor *desc, 
                                VimosDistModelFull **crvMod)

   Description: 
   Read the Curvature Model from the descriptors
   
   Input:
   VimosDescriptor *desc
   Pointer to linked list of descriptors that contain the coefficients of the 
   Curvature Model
   
   Output:
   VimosDistModelFull **crvMod
   Contains the Curvature Model

   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   24 Nov 98: Created (TAO)
      
--------------------------------------------------------------------------------
*/
VimosBool readCurvatureModel(VimosDescriptor *desc, 
                             VimosDistModelFull **crvMod);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   
   VimosBool writeCurvatureModel(VimosDescriptor **desc, 
                                 VimosDistModelFull *crvMod)

   Description: 
   Write a  Curvature Model to the descriptors
   
   Input:
   VimosDescriptor **desc
   Address of pointer to Linked list of descriptors to write to
   
   Output:
   VimosDistModelFull *crvMod
   Contains the Curvature Model

   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   24 Nov 98: Created (TAO)
      
--------------------------------------------------------------------------------
*/
VimosBool writeCurvatureModel(VimosDescriptor **desc, 
                              VimosDistModelFull *crvMod);
VimosBool writeCurvatureModelString(VimosDescriptor **desc, 
                                    VimosDistModelFull *crvMod);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBool readInvDispMatrix(VimosDescriptor *desc, 
                               VimosDistModelFull **invDispMat)

   Description: 
   Read the Inverse Dispersion Matrix from the descriptors
   
   Input:
   VimosDescriptor *desc
   Linked list of descriptors that contain the coefficients of the 
   Inverse Dispersion Matrix
   
   Output:
   VimosDistModelFull **invDispMat
   Contains the Inverse Dispersion Matrix

   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   24 Nov 98: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosBool readInvDispMatrix(VimosDescriptor *desc, 
                            VimosDistModelFull **invDispMat);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   
   VimosBool writeInvDispMatrix(VimosDescriptor **desc, 
                                VimosDistModelFull *idsMat)

   Description: 
   Write an Inverse Dispersion Relation Matrix to the descriptors
   
   Input:
   VimosDescriptor **desc
   Address of pointer to Linked list of descriptors to write to
   
   Output:
   VimosDistModelFull *idsMat
   Contains the Inverse Dispersion Relation Matrix

   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   24 Nov 98: Created (TAO)
      
--------------------------------------------------------------------------------
*/
VimosBool writeInvDispMatrix(VimosDescriptor **desc, 
                             VimosDistModelFull *idsMat);
VimosBool writeInvDispMatrixString(VimosDescriptor **desc, 
                             VimosDistModelFull *idsMat);



/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosGnomonic    *newGnomonic(const double alpha0, const double delta0)

   Description: 
   Allocates a new VimosGnomonic structure for projection centre (alpha0, delat0)
   
   Input:
   const double alpha0    RA of projection centre (in degrees!!)
   const double delta0    Dec of projection centre (in degrees!!)

   Return Value (succes):
   Pointer to a newly allocated VimosGnomonic structure with members set to
   correct values for projection centre given as input.

   Return Value (error):
   NULL
   
   Updates:
   13 Jun 00: Return NULL on error instead exiting (Maura)
   25 feb 99: Created (TAO)

--------------------------------------------------------------------------------
*/
VimosGnomonic   *newGnomonic(const double alpha0, const double delta0);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void  *deleteGnomonic(VimosGnomonic *gnome)

   Description: 
   Desctructor of VimosGnomonic
   
   Input:
   VimosGnomonic *gnome
   Pointer of structure to be deleted
   
   Return Value:
   void
   
   Updates:
   25 feb 99: Created (TAO)
   
--------------------------------------------------------------------------------
*/
void deleteGnomonic(VimosGnomonic *gnome);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBool  lm2RaDec(VimosGnomonic *gnome, double l, double m, 
                       double *ra, double *dec)

   Description: 
   Convert gnomonic coordinates in Ra and Dec. The coordinate system of Virmos 
   is modeled as having a gnomonic coordinate system in the focal plane (the
   l,m coordinates). This gnomonic coordinate system serevs as an intermediate 
   system. The relation between the gnomonic coordinates and pixel coordiantes 
   on the CCD or the mask coordinates is modeled as a VimosDistModel2D.
   
   Input:
   VimosGnomonic *gnome
   Pointer of VimosGnomonic that contains the RA and Dec of the projection
   centre (which should be de position on the sky of the telescope axis, not
   the centre of the CCD) 

   double l
   double m
   Input l,m coordinates

   Output
   double *ra
   double *dec
   Output Ra and Dec
   
   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   25 feb 99: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosBool lm2RaDec(VimosGnomonic *gnome, double l, double m, 
                   double *ra, double *dec);


/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   void  deleteGnomonic(VimosGnomonic *gnome)


   Description: 
   Convert Ra and Dec in gnomonic coordinates. The coordinate system of Virmos 
   is modeled as having a gnomonic coordinate system in the focal plane (the
   l,m coordinates). This gnomonic coordinate system serevs as an intermediate 
   system. The relation between the gnomonic coordinates and pixel coordiantes 
   on the CCD or the mask coordinates is modeled as a VimosDistModel2D.
   
   Input:
   VimosGnomonic *gnome
   Pointer of VimosGnomonic that contains the RA and Dec of the projection
   centre (which should be de position on the sky of the telescope axis, not
   the centre of the CCD) 

   double ra
   double dec
   Input Ra and Dec

   Output
   double *l
   double *m
   Output l,m coordinates
   
   Return Value:
   void
   
   Updates:
   25 feb 99: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosBool raDec2lm(VimosGnomonic *gnome, double ra, double dec, 
                   double *l, double *m);
     

/* 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   VimosBool  fitDistModel2D(VimosPixel *surface, int numPoints, int polyDeg,  
                             double offsetX, double offsetY, 
                             VimosDistModel2D **model, double *rms)

   Description: 
   Determine the coefficients of the VimosDistModel2D that connects e.g. one 
   mask coordinate (e.g. in X)  to the gnomonic coordinates of Virmos
   
   Input:
   VimosPixel *surface
   Pointer to VimosPixel that contains the X or Y position of calibration
   positions 

   int numPoints
   Number of input points (i.e. length of *surface)
   
   int polyDeg
   Degree to fit to surface. Only those terms of x^i*y^j with i+j <= polyDeg
   are fitted
   
   double offsetX
   double offsetY
   Offset to be applied to the x and y fields of the surface before fitting. 
   This typically is te position of the telescope axis.

   Output
   VimosDistModel **model
   Model containing the coefficients of the fit

   double *rms
   RMS deviation from fit

   Return Value (succes):
   VM_TRUE

   Return Value (error):
   VM_FALSE
   
   Updates:
   13 Jun 00: Return VimosBool instead of void (Maura)
   25 feb 99: Created (TAO)
   
--------------------------------------------------------------------------------
*/
VimosBool fitDistModel2D(VimosPixel *surface, int numPoints, int polyDeg,  
                         double offsetX, double offsetY, 
                         VimosDistModel2D **model, double *rms);

PIL_END_DECLS

#endif /* VM_DISTMODELS_H */

