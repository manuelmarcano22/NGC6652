/* $Id: piltranslator.c,v 1.4 2012-11-08 17:57:37 cgarcia Exp $
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
 * $Date: 2012-11-08 17:57:37 $
 * $Revision: 1.4 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "piltranslator.h"
#include "pilmemory.h"
#include "pilmessages.h"
#include "pilutils.h"


#define COMMENT_CHARS   "#"
#define KEYMAP_ALIAS    "Alias:"
#define KEYMAP_NAME     "Parameter Name:"
#define KEYMAP_FORM     "Value Format:"
#define KEYMAP_COMM     "Comment Field:"
#define CATMAP_ALIAS    "Alias:"
#define CATMAP_NAME     "Category Name:"
#define MAX_TRN_LENGTH  (1024)


static PilKeymap *keymap;
static PilCatmap *catmap;


/**
 * @defgroup pilTranslator pilTranslator
 *
 * The module @b pilTranslator is used to define a dictionary of
 * data categories and keywords, and to return the translation
 * from names to their mapped values. Names can be indexed, and
 * in that case their mapped value would contain a C-like integer 
 * format (%d) where the index(es) would be inserted.
 */

/**@{*/

/**
 * @brief
 *   Add to the keyword mapping a keyword with given name, value,
 *   format and mapping. If a keyword with the same name is already
 *   in the dictionary, replace its value, format and mapping with 
 *   the new ones.
 *
 * @param name Keyword name.
 * @param value Keyword value.
 * @param format Keyword format.
 * @param comment Keyword comment.
 *
 * @return @c EXIT_SUCCESS on success, @c EXIT_FAILURE on failure.
 *
 * The function creates and adds the new keyword to the mapping dictionary.
 */

int
pilTrnAddKey(const char *name, const char *value, const char *format,
             const char *comment)
{

  PilAlias  *alias;
  
  if ((alias = pilKeymapLookup(keymap, name))) {
    pilAliasSetValue(alias, value);
    pilAliasSetFormat(alias, format);
    pilAliasSetComment(alias, comment);

    return EXIT_SUCCESS;
  }
  else {
    alias = newPilAlias(name, value, format, comment);

    return pilKeymapInsert(keymap, alias);
  }

}


/**
 * @brief
 *   Create a default keyword mapping dictionary.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * This function would be called if no configuration file defining
 * the keyword mapping were found. Just one keyword mapping definition
 * is allowed in one pipeline run, an error condition will occur
 * otherwise.
 */

int
pilTrnInitKeywordMap(void)
{
  char modName[] = "pilTrnInitKeywordMap";

  if (keymap == NULL) {
    keymap = newPilKeymap();
  }
  else {
    pilMsgWarning(modName, 
    "Double definition of keywords mapping - the first one is taken");
    return EXIT_FAILURE;
  }


  /*
   * Basic FITS keywords
   */

  pilTrnAddKey("Alpha", "RA", "%.6f", "Telescope pointing RA (J2000)");
  pilTrnAddKey("Delta", "DEC", "%.5f", "Telescope pointing DEC (J2000)");
  pilTrnAddKey("Instrument", "INSTRUME", "%s", "Instrument used");
  pilTrnAddKey("ExposureTime", "EXPTIME", "%.3f", "Integration time (s)"); 
  pilTrnAddKey("SiderialTime", "LST", "%.3f", "Local siderial time (s)"); 
  pilTrnAddKey("Origin", "ORIGIN", "%30s", "European Southern Observatory");
  pilTrnAddKey("Date", "DATE", "%s", "Date this file was written");
  pilTrnAddKey("MjdObs", "MJD-OBS", "%f", "Start of observation");
  pilTrnAddKey("DateObs", "DATE-OBS", "%30s", "Date of observation");
  pilTrnAddKey("Equinox", "EQUINOX", "%f", "equinox");
  pilTrnAddKey("Radecsys", "RADECSYS", "%s", "reference frame");
  pilTrnAddKey("Longpole", "LONGPOLE", "%f", "Longitude of Npole");
  pilTrnAddKey("CO1", "CO1_%d", "%f", "X Coeffs of the wcs distorsion model"); 
  pilTrnAddKey("CO2", "CO2_%d", "%f", "Y Coeffs of the wcs distorsion model");
  pilTrnAddKey("Naxis", "NAXIS%d", "%f", "number of pixels");
  pilTrnAddKey("Crpix", "CRPIX%d", "%f", "Pix coords at ref point");
  pilTrnAddKey("Crval", "CRVAL%d", "%f", "RA/Dec at ref point");
  pilTrnAddKey("Ctype", "CTYPE%d", "%s", "Projection Type");
  pilTrnAddKey("Cdelt", "CDELT%d", "%f", "Coord increment at ref point");
  pilTrnAddKey("CD", "CD%d_%d", "%s", "Rotation matrix");
  pilTrnAddKey("DataMin", "DATAMIN", "%.3f", "Minimum pixel value");
  pilTrnAddKey("DataMax", "DATAMAX", "%.3f", "Maximum pixel value");

  /*
   * Observation
   */

  pilTrnAddKey("OBS.DID", "ESO OBS DID", "%s", "OBS Dictionary");
  pilTrnAddKey("OBS.ID", "ESO OBS ID", "%d", "Observation block ID");
  pilTrnAddKey("PROG.ID", "ESO OBS PROG ID", "%20s",
               "ESO program identification");
  pilTrnAddKey("TargetName", "ESO OBS TARG NAME", "%20s", "Target name");

  /*
   * Template
   */

  pilTrnAddKey("TplId", "ESO TPL ID", "%d", "Template ID");
  pilTrnAddKey("TplExposures", "ESO TPL NEXP", "%d",
               "Number of exposures in template");
  pilTrnAddKey("TplExposureNumber", "ESO TPL EXPNO", "%d",
               "Exposure number within template");

  /*
   * Data Product keywords
   */

  pilTrnAddKey("DprCategory", "ESO DPR CATG", "%32s",
               "Data Product category");
  pilTrnAddKey("DprType", "ESO DPR TYPE", "%32s",
               "Data Product type");
  pilTrnAddKey("DprTechnique", "ESO DPR TECH", "%32s",
               "Data Product technique");

  /*
   * Archive
   */

  pilTrnAddKey("OriginalFile", "ORIGFILE", "%40s",
               "Original file name");
  pilTrnAddKey("ArchiveFile", "ARCFILE", "%40s",
               "Archive file name");
  pilTrnAddKey("Checksum", "CHECKSUM", "%16s",
               "ASCII 1's complement checksum");

  /*
   * Process
   */

  pilTrnAddKey("DataMD5", "DATAMD5", "%32s",
               "MD5 data signature");
  pilTrnAddKey("ProductFile", "PIPEFILE", "%64s",
               "Filename of data product");
  pilTrnAddKey("ProductDID", "ESO PRO DID", "%30s",
               "Data dictionary for PRO");
  pilTrnAddKey("DoCategory", "ESO PRO CATG", "%30s",
               "Pipeline product category");
  pilTrnAddKey("DrsId", "ESO PRO REC%d DRS ID", "%20s",
               "Data reduction system identifier");
  pilTrnAddKey("PipelineId", "ESO PRO REC%d PIPE ID", "%32s",
               "Pipeline (unique) identifier");
  pilTrnAddKey("RecipeId", "ESO PRO REC%d ID", "%32s",
               "Pipeline recipe (unique) identifier");
  pilTrnAddKey("RecipeStart", "ESO PRO REC%d START", "%30s",
               "Date when recipe execution started.");
  pilTrnAddKey("RawFrameId", "ESO PRO REC%d RAW%d NAME", "%45s",
               "File name of raw frame.");
  pilTrnAddKey("RawFrameCategory", "ESO PRO REC%d RAW%d CATG", "%32s",
               "Frame category of raw frame.");
  pilTrnAddKey("CalFrameId", "ESO PRO REC%d CAL%d NAME", "%45s",
               "File name of calibration frame.");
  pilTrnAddKey("CalFrameCategory", "ESO PRO REC%d CAL%d CATG", "%32s",
               "Frame category of calibration frame.");
  pilTrnAddKey("CalFrameMD5", "ESO PRO REC%d CAL%d DATAMD5", "%32s",
               "MD5 data signature of calibration frame.");
  pilTrnAddKey("ProductType","ESO PRO TYPE", "%20s",
               "Product Type");
  pilTrnAddKey("ProductSource","ESO PRO ARCFILE", "%48s",
               "Archive file name of raw data frame.");
  pilTrnAddKey("DataMean", "ESO PRO DATAAVG", "%.3f",
               "Mean pixel value");
  pilTrnAddKey("DataStdDeviation", "ESO PRO DATARMS", "%.3f",
               "Standard deviation of pixel values");
  pilTrnAddKey("DataMedian", "ESO PRO DATAMED", "%.3f",
               "Median pixel value");
  pilTrnAddKey("NFramesCombined","ESO PRO DATANCOM", "%d",
               "Number of frames combined");
  pilTrnAddKey("SummedExposureTime","ESO PRO EXPTTOT", "%.2f",
               "Total exposure time of all frames combined");

  /*
   * Telescope
   */

  pilTrnAddKey("Longitude", "ESO TEL GEOLON", "%.4f", "Telescope geographic "
               " longitude (+=East) in degrees");
  pilTrnAddKey("Latitude", "ESO TEL GEOLAT", "%.4f", "Telescope geographic "
               " latitude (+=North) in degrees");
  pilTrnAddKey("AmbientTemperature", "ESO TEL AMBI TEMP", "%.2f",
               "Observatory ambient temperature queried from ASM");
  pilTrnAddKey("SeeingStart", "ESO TEL AMBI FWHM START", "%.2f",
               "Site seeing at observation start queried from ASM");
  pilTrnAddKey("SeeingEnd", "ESO TEL AMBI FWHM END", "%.2f",
               "Site seeing at observation end queried from ASM");
  pilTrnAddKey("FocalScale", "ESO TEL FOCU SCALE", "%.4f",
               "arcsecs per mm on the focal plane (mask)");

  /*
   * Adapter
   */

  /*
   * Instrument
   */

  pilTrnAddKey("INS.DID", "ESO INS DID", "%30s",
               "Data dictionary for INS");
  pilTrnAddKey("BeamTemperature", "ESO INS FOCU%d TEMP", "%f",
               "Temperature(C) for the beam");
  pilTrnAddKey("MaskId", "ESO INS MASK%d ID", "%d",
               "Mask Identifier");
  pilTrnAddKey("NoRefSlit", "ESO INS REF NO", "%d",
               "No of ref apertures"); 
  pilTrnAddKey("RefSlitX", "ESO INS REF%d X", "%f",
               "X coord in mm");
  pilTrnAddKey("RefSlitY", "ESO INS REF%d Y", "%f",
               "Y coord in mm");
  pilTrnAddKey("SlitRefDimX", "ESO INS REF%d DIMX", "%f",
               "X size in mm");
  pilTrnAddKey("SlitRefDimY", "ESO INS REF%d DIMY", "%f",
               "Y size in mm");
  pilTrnAddKey("SlitRefObjRA", "ESO INS REF%d OBJ RA", "%f",
               "?");
  pilTrnAddKey("SlitRefObjDec", "ESO INS REF%d OBJ DEC", "%f",
               "?");
  pilTrnAddKey("NoSlit", "ESO INS SLIT NO", "%d",
               "No of slits"); 
  pilTrnAddKey("SlitDimX", "ESO INS SLIT%d DIMX", "%f",
               "X size in mm");
  pilTrnAddKey("SlitDimY", "ESO INS SLIT%d DIMY", "%f",
               "Y size in mm");
  pilTrnAddKey("SlitRadius", "ESO INS SLIT%d RADIUS", "%f",
               "Slit radius in mm");
  pilTrnAddKey("SlitX", "ESO INS SLIT%d X", "%f",
               "X position in mm");
  pilTrnAddKey("SlitY", "ESO INS SLIT%d Y", "%f",
               "Y position in mm");
  pilTrnAddKey("SlitType", "ESO INS SLIT%d TYPE", "%s",
               "Type of slit");
  pilTrnAddKey("SlitId", "ESO INS SLIT%d ID", "%s",
               "ID of slit");
  pilTrnAddKey("SlitObjRA", "ESO INS SLIT%d OBJ RA", "%f",
               "RA (J2000) of object");
  pilTrnAddKey("SlitObjDec", "ESO INS SLIT%d OBJ DEC", "%f",
               "DEC (J2000) of object");
  pilTrnAddKey("SlitBezierDY", "ESO INS SLIT%d BEZIER DY", "%f",
               "Bezier slit width");
  pilTrnAddKey("SlitBezierXX", "ESO INS SLIT%d BEZIER XX", "%f",
               "Bezier slit coefficient");
  pilTrnAddKey("SlitBezierAX", "ESO INS SLIT%d BEZIER AX", "%f",
               "Bezier slit coefficient");
  pilTrnAddKey("SlitBezierBX", "ESO INS SLIT%d BEZIER BX", "%f",
               "Bezier slit coefficient");
  pilTrnAddKey("SlitBezierCX", "ESO INS SLIT%d BEZIER CX", "%f",
               "Bezier slit coefficient");
  pilTrnAddKey("SlitBezierYY", "ESO INS SLIT%d BEZIER YY", "%f",
               "Bezier slit coefficient");
  pilTrnAddKey("SlitBezierAY", "ESO INS SLIT%d BEZIER AY", "%f",
               "Bezier slit coefficient");
  pilTrnAddKey("SlitBezierBY", "ESO INS SLIT%d BEZIER BY", "%f",
               "Bezier slit coefficient");
  pilTrnAddKey("SlitBezierCY", "ESO INS SLIT%d BEZIER CY", "%f",
               "Bezier slit coefficient");

  pilTrnAddKey("MshuRefL", "ESO OCS CON SHUL POS%d","%f",
               "Shutter reference pos (L)");
  pilTrnAddKey("MshuRefH", "ESO OCS CON SHUH POS%d","%f",
               "Shutter reference pos (H)");
  pilTrnAddKey("MshuPosL", "ESO INS MSHU%d POSL","%f",
               "Shutter pos (L)");
  pilTrnAddKey("MshuPosH", "ESO INS MSHU%d POSH","%f",
               "Shutter pos (H)");
  pilTrnAddKey("MshuMode", "ESO INS MSHU%d MODE", "%s",
               "Shutter mode: ON or OFF"); 
  pilTrnAddKey("PixelScale", "ESO INS PIXSCALE", "%f",
               "Pixelscale arcsec/pix");
  pilTrnAddKey("InstrumentMode", "ESO INS MODE", "%s",
               "Instrument mode");
  pilTrnAddKey("IfuMode", "ESO INS IFUS MODE", "%s",
               "IFU mode");
  pilTrnAddKey("AdfId", "ESO INS ADF ID", "%d",
               "ADF Identifier");
  pilTrnAddKey("Quadrant", "ESO OCS CON QUAD", "%d",
               "Quadrant used");

  /*
   * Detector
   */

  pilTrnAddKey("DET.DID", "ESO DET DID", "%30s",
               "Dictionary Name and Revision");
  pilTrnAddKey("NCHIPS", "ESO DET CHIPS", "%d",
               "# of chips in detector array");
  pilTrnAddKey("CHIPi.INDEX", "ESO DET CHIP%d INDEX", "%d",
               "Chip index");
  pilTrnAddKey("CHIPi.ID", "ESO DET CHIP%d ID", "%d",
               "Chip index");
  pilTrnAddKey("CHIPi.NAME", "ESO DET CHIP%d NAME", "%16s",
               "Detector chip name");
  pilTrnAddKey("CHIPi.DATE", "ESO DET CHIP%d DATE", "%10s",
               "Date of installation [YYYY-MM-DD]");
  pilTrnAddKey("CHIPi.NX", "ESO DET CHIP%d NX", "%d",
               "# of pixels along X");
  pilTrnAddKey("CHIPi.NY", "ESO DET CHIP%d NY", "%d",
               "# of pixels along Y");

  /* The 2 following entries are deprecated don't use them any more */
  pilTrnAddKey("CHIP1.NX", "ESO DET CHIP NX","%d",
               "# of pixels along X");
  pilTrnAddKey("CHIP1.NY", "ESO DET CHIP NY","%d",
               "# of pixels along Y");

  pilTrnAddKey("CHIPi.PSZX", "ESO DET CHIP%d PSZX", "%d",
               "Size of pixel in X");
  pilTrnAddKey("CHIPi.PSZY", "ESO DET CHIP%d PSZY", "%d",
               "Size of pixel in Y");
  pilTrnAddKey("NumberOfPorts", "ESO DET OUTPUTS", "%d",
               "Total number of ports");
  pilTrnAddKey("PortStartX", "ESO DET OUT1 X", "%d",
               "X start of port");
  pilTrnAddKey("PortStartY", "ESO DET OUT1 Y", "%d",
               "Y start of port");
  pilTrnAddKey("PortSizeX", "ESO DET OUT1 NX", "%d",
               "X port size");
  pilTrnAddKey("PortSizeY", "ESO DET OUT1 NY", "%d",
               "Y port size");
  pilTrnAddKey("PortOvscX", "ESO DET OUT1 OVSCX", "%d",
               "Overscan in X");
  pilTrnAddKey("PortOvscY", "ESO DET OUT1 OVSCY", "%d",
               "Overscan in Y");
  pilTrnAddKey("PortPrscX", "ESO DET OUT1 PRSCX", "%d",
               "Prescan in X");
  pilTrnAddKey("PortPrscY", "ESO DET OUT1 PRSCY", "%d",
               "Prescan in Y");
  pilTrnAddKey("SeqPortStartX", "ESO DET OUT%d X", "%d",
               "X start of port");
  pilTrnAddKey("SeqPortStartY", "ESO DET OUT%d Y", "%d",
               "Y start of port");
  pilTrnAddKey("SeqPortSizeX", "ESO DET OUT%d NX", "%d",
               "X port size");
  pilTrnAddKey("SeqPortSizeY", "ESO DET OUT%d NY", "%d",
               "Y port size");
  pilTrnAddKey("SeqPortOvscX", "ESO DET OUT%d OVSCX", "%d",
               "Overscan in X");
  pilTrnAddKey("SeqPortOvscY", "ESO DET OUT%d OVSCY", "%d",
               "Overscan in Y");
  pilTrnAddKey("SeqPortPrscX", "ESO DET OUT%d PRSCX", "%d",
               "Prescan in X");
  pilTrnAddKey("SeqPortPrscY", "ESO DET OUT%d PRSCY", "%d",
               "Prescan in Y");
  pilTrnAddKey("NumberOfWindows", "ESO DET WINDOWS", "%d",
               "Number of readout windows");
  pilTrnAddKey("SeqWindowStartX", "ESO DET WIN%d STRX", "%d",
               "X start of window");
  pilTrnAddKey("SeqWindowStartY", "ESO DET WIN%d STRY", "%d",
               "Y start of window");
  pilTrnAddKey("SeqWindowSizeX", "ESO DET WIN%d NX", "%d",
               "X window size");
  pilTrnAddKey("SeqWindowSizeY", "ESO DET WIN%d NY", "%d",
               "Y window size");
  pilTrnAddKey("WINi.BINX", "ESO DET WIN%d BINX", "%d",
               "Binning factor along X");
  pilTrnAddKey("WINi.BINY", "ESO DET WIN%d BINY", "%d",
               "Binning factor along Y");
  pilTrnAddKey("READ.CLOCK", "ESO DET READ CLOCK", "%32s",
               "Readout clock pattern used");
  pilTrnAddKey("READ.MODE", "ESO DET READ MODE", "%10s",
               "Readout method");
  pilTrnAddKey("READ.SPEED", "ESO DET READ SPEED", "%10s",
               "Readout speed");
  pilTrnAddKey("SeqBiasLevel", "ESO DET OUT%d BIAS", "%.2f",
               "Mean bias level");
  pilTrnAddKey("SeqDarkLevel", "ESO DET OUT%d DARK", "%.2f",
               "Mean dark level");
  pilTrnAddKey("ReadNoise", "ESO DET OUT%d RON", "%.2f",
               "Read out noise per output (e-)");
  pilTrnAddKey("Adu2Electron", "ESO DET OUT%d CONAD", "%.2f",
               "Conversion from ADUs to electrons");
  pilTrnAddKey("Electron2Adu", "ESO DET OUT%d GAIN", "%.2f",
               "Conversion from electrons to ADU");

  /*
   * VIMOS DRS specific keywords
   */

  pilTrnAddKey("Name","ESO PRO NAME", "%s", 
               "Product name");

  pilTrnAddKey("IdsFlag","ESO PRO IDS FLAG", "%d", 
               "IDS creation: 0 = based on first guess distorsion models");

  pilTrnAddKey("IdsDate", "ESO PRO IDS MAT DAYTIM", "%d", 
               "Inverse Dispersion Solution Time");
  pilTrnAddKey("IdsTemp", "ESO PRO IDS MAT TEMP", "%d", 
               "Inverse Dispersion Solution Temperature");
  pilTrnAddKey("DispersionOrd",  "ESO PRO IDS REL ORD", "%d",
               "IDS relation polynomial order");
  pilTrnAddKey("DispersionOrdX", "ESO PRO IDS MAT XORD", "%d",
               "X order of the model of the IDS relation polynomial "
               "coefficients");
  pilTrnAddKey("DispersionOrdY", "ESO PRO IDS MAT YORD", "%d",
               "Y order of the model of the IDS relation polynomial "
               "coefficients");
  pilTrnAddKey("Dispersion", "ESO PRO IDS MAT_%d_%d_%d", "%d",
               "Coefficient of the model of the IDS relation coefficients");
  pilTrnAddKey("IdsXrms", "ESO PRO IDS MAT X_RMS", "%d", 
               "Inverse Dispersion Solution X RMS");
  pilTrnAddKey("IdsYrms", "ESO PRO IDS MAT Y_RMS", "%d", 
               "Inverse Dispersion Solution Y RMS");

  pilTrnAddKey("OptDate", "ESO PRO OPT DIS DAYTIM", "%d", 
               "Optical Distortion Model Time");
  pilTrnAddKey("OptTemp","ESO PRO OPT DIS TEMP", "%d", 
               "Optical Distorsion Temperature");
  pilTrnAddKey("OptDistOrdX", "ESO PRO OPT DIS XORD", "%d",
               "Optical distorsion X polynomial order");
  pilTrnAddKey("OptDistOrdY", "ESO PRO OPT DIS YORD", "%d",
               "Optical distorsion Y polynomial order");
  pilTrnAddKey("OptDistX", "ESO PRO OPT DIS X_%d_%d", "%d",
               "Optical distorsion coefficient of X polynomial");
  pilTrnAddKey("OptDistY", "ESO PRO OPT DIS Y_%d_%d", "%d",
               "Optical distorsion coefficient of Y polynomial");
  pilTrnAddKey("OptXrms", "ESO PRO OPT DIS X_RMS", "%d", 
               "Optical distortion model X RMS");
  pilTrnAddKey("OptYrms", "ESO PRO OPT DIS Y_RMS", "%d", 
               "Optical distortion model Y RMS");

  pilTrnAddKey("CurvatureOrd", "ESO PRO CRV POL ORD", "%d",
               "Curvature polynomial order");
  pilTrnAddKey("CurvatureOrdX", "ESO PRO CRV MOD XORD", "%d",
               "X order of the model of the curvature polynomial "
               "coefficients");
  pilTrnAddKey("CurvatureOrdY", "ESO PRO CRV MOD YORD", "%d",
               "Y order of the model of the curvature polynomial "
               "coefficients");
  pilTrnAddKey("Curvature", "ESO PRO CRV MOD_%d_%d_%d", "%d",
               "Coefficient of the model of the curvature polynomial "
               "coefficients");
  pilTrnAddKey("CurvatureXrms", "ESO PRO CRV MOD X_RMS", "%d", 
               "Curvature model X RMS");
  pilTrnAddKey("CurvatureYrms", "ESO PRO CRV MOD Y_RMS", "%d", 
               "Curvature model Y RMS");

  pilTrnAddKey("ZeroOrderFlag", "ESO PRO ZERO", "%d",
               "Contamination flag");
  pilTrnAddKey("ZeroOrdX", "ESO PRO ZERO XORD", "%d", 
               "X order of contamination model");
  pilTrnAddKey("ZeroOrdY", "ESO PRO ZERO YORD", "%d", 
               "Y order of contamination model");
  pilTrnAddKey("ZeroX", "ESO PRO ZERO X_%d_%d", "%d", 
               "X branch coefficients of the contamin. model");
  pilTrnAddKey("ZeroY", "ESO PRO ZERO Y_%d_%d", "%d", 
               "Y branch coefficients of the contamin. model");
  pilTrnAddKey("ZeroXrms", "ESO PRO ZERO X_RMS", "%d", 
               "Contamination model X RMS");
  pilTrnAddKey("ZeroYrms", "ESO PRO ZERO Y_RMS", "%d", 
               "Contamination model Y RMS");
  pilTrnAddKey("ZeroOrderWidth", "ESO PRO ZERO PWIDTH", "%f", 
               "Width of the Zero Order contamination");

  pilTrnAddKey("ObsType", "ESO PRO OBS TYPE", "%s",
               "Observation type");

  pilTrnAddKey("NumPixBelow", "ESO PRO SPECT LLEN LO", "%d",
               "Number of pixels below position of central wavelength");
  pilTrnAddKey("NumPixAbove", "ESO PRO SPECT LLEN HI", "%d",
               "Number of pixels above position of central wavelength");
  pilTrnAddKey("Table", "ESO PRO TABLE", "%s",
               "Type of table");
  pilTrnAddKey("BiasLevel", "ESO PRO BIAS", "%.2f",
               "Bias level");
  pilTrnAddKey("DarkLevel", "ESO PRO DARK", "%.2f",
               "Mean dark level");
  pilTrnAddKey("GrismId", "ESO INS GRIS%d ID", "%s",
               "Grism Identifier");
  pilTrnAddKey("GrismName", "ESO INS GRIS%d NAME", "%s",
               "Grism Name");
  pilTrnAddKey("FilterId", "ESO INS FILT%d ID", "%s",
               "Filter Identifier");
  pilTrnAddKey("FilterName", "ESO INS FILT%d NAME", "%s",
               "Filter Name");
  pilTrnAddKey("Colour", "ESO PRO COLOUR", "%s",
               "Color system");

  pilTrnAddKey("RotatorAngleStart", "ESO ADA ABSROT START", "%s",
               "Adaptor initial absolute rotator angle");
  pilTrnAddKey("RotatorAngleEnd", "ESO ADA ABSROT END", "%s",
               "Adaptor final absolute rotator angle");

  pilTrnAddKey("CcdMaskTime","ESO PRO CCD MASK DAYTIM", "%d", 
               "CCD To Mask Matrix time stamp");
  pilTrnAddKey("CcdMaskTemp","ESO PRO CCD MASK TEMP", "%d", 
               "CCD To Mask Matrix Temperature");
  pilTrnAddKey("MaskCcdX0", "ESO PRO MASK CCD X0", "%s",
               "X offset coefficient");
  pilTrnAddKey("MaskCcdXX", "ESO PRO MASK CCD XX", "%s",
               "X scale coefficient");
  pilTrnAddKey("MaskCcdXY", "ESO PRO MASK CCD XY", "%s",
               "X rotation coefficient");
  pilTrnAddKey("MaskCcdY0", "ESO PRO MASK CCD Y0", "%s",
               "Y offset coefficient");
  pilTrnAddKey("MaskCcdYY", "ESO PRO MASK CCD YY", "%s",
               "Y scale coefficient");
  pilTrnAddKey("MaskCcdYX", "ESO PRO MASK CCD YX", "%s",
               "Y rotation coefficient");
  pilTrnAddKey("MaskCcdX", "ESO PRO MASK CCD X_%d_%d", "%s",
               "X mask to CCD");
  pilTrnAddKey("MaskCcdY", "ESO PRO MASK CCD Y_%d_%d", "%s",
               "Y mask to CCD");
  pilTrnAddKey("MaskCcdXord", "ESO PRO MASK CCD XORD", "%d",
               "Order X mask to CCD");
  pilTrnAddKey("MaskCcdYord", "ESO PRO MASK CCD YORD", "%d",
               "Order Y mask to CCD");
  pilTrnAddKey("MaskCcdXrms", "ESO PRO MASK CCD XRMS", "%.2f",
               "rms error in X");
  pilTrnAddKey("MaskCcdYrms", "ESO PRO MASK CCD YRMS", "%.2f",
               "rms error in Y");
  pilTrnAddKey("CcdMaskX0", "ESO PRO CCD MASK X0", "%s",
               "X offset coefficient");
  pilTrnAddKey("CcdMaskXX", "ESO PRO CCD MASK XX", "%s",
               "X scale coefficient");
  pilTrnAddKey("CcdMaskXY", "ESO PRO CCD MASK XY", "%s",
               "X rotation coefficient");
  pilTrnAddKey("CcdMaskY0", "ESO PRO CCD MASK Y0", "%s",
               "Y offset coefficient");
  pilTrnAddKey("CcdMaskYY", "ESO PRO CCD MASK YY", "%s",
               "Y scale coefficient");
  pilTrnAddKey("CcdMaskYX", "ESO PRO CCD MASK YX", "%s",
               "Y rotation coefficient");
  pilTrnAddKey("CcdMaskX", "ESO PRO CCD MASK X_%d_%d", "%s",
               "X CCD to mask");
  pilTrnAddKey("CcdMaskY", "ESO PRO CCD MASK Y_%d_%d", "%s",
               "Y CCD to mask");
  pilTrnAddKey("CcdMaskXrms", "ESO PRO CCD MASK XRMS", "%.2f",
               "rms error in X");
  pilTrnAddKey("CcdMaskYrms", "ESO PRO CCD MASK YRMS", "%.2f",
               "rms error in Y");
  pilTrnAddKey("CcdMaskXord", "ESO PRO CCD MASK XORD", "%d",
               "Order X CCD to mask");
  pilTrnAddKey("CcdMaskYord", "ESO PRO CCD MASK YORD", "%d",
               "Order Y CCD to mask");

  pilTrnAddKey("InvCO1rms","ESO PRO INVCO1 RMS", "%d",       
               "Inverse CO1 RMS");                           
  pilTrnAddKey("InvCO2rms","ESO PRO INVCO2 RMS", "%d",       
               "Inverse CO2 RMS");                           
  pilTrnAddKey("CO1rms","ESO PRO CO1 RMS", "%d",             
               "CO1 RMS");                                   
  pilTrnAddKey("CO2rms","ESO PRO CO2 RMS", "%d",             
               "CO2 RMS");                                   

  pilTrnAddKey("CcdSkyTime","ESO PRO CCD SKY DAYTIM", "%d", 
               "Sky to CCD Matrix Time");
  pilTrnAddKey("CcdSkyTemp","ESO PRO CCD SKY TEMP", "%d", 
               "Sky to CCD Matrix Temperature");
  pilTrnAddKey("CcdSkyXord", "ESO PRO CCD SKY XORD", "%d", 
               "No of coef X CCD to Sky");   
  pilTrnAddKey("CcdSkyYord", "ESO PRO CCD SKY YORD", "%d", 
               "No of coef Y CCD to Sky");
  pilTrnAddKey("CcdSkyX", "ESO PRO CCD SKY X_%d_%d", "%.2f",
               "X CCD to Sky");
  pilTrnAddKey("CcdSkyY", "ESO PRO CCD SKY Y_%d_%d", "%.2f",
               "Y CCD to Sky");
  pilTrnAddKey("CcdSkyXrms", "ESO PRO CCD SKY XRMS", "%.2f", 
               "RMS in X sky fit (arcsec)");
  pilTrnAddKey("CcdSkyYrms", "ESO PRO CCD SKY YRMS", "%.2f", 
               "RMS in Y sky fit (arcsec)");
  pilTrnAddKey("SkyCcdXord", "ESO PRO SKY CCD XORD", "%d", 
               "No of coef X Sky to CCD");   
  pilTrnAddKey("SkyCcdYord", "ESO PRO SKY CCD YORD", "%d", 
               "No of coef Y Sky to CCD");
  pilTrnAddKey("SkyCcdX", "ESO PRO SKY CCD X_%d_%d", "%.2f",
               "X Sky to CCD");
  pilTrnAddKey("SkyCcdY", "ESO PRO SKY CCD Y_%d_%d", "%.2f",
               "Y Sky to CCD");
  pilTrnAddKey("SkyCcdXrms", "ESO PRO SKY CCD XRMS", "%.2f", 
               "RMS in X ccd fit (pix)");
  pilTrnAddKey("SkyCcdYrms", "ESO PRO SKY CCD YRMS", "%.2f", 
               "RMS in Y ccd fit (pix)"); 

  pilTrnAddKey("SkyMaskXord", "ESO PRO TEL DIS XORD", "%d", 
               "Order X Sky to Mask");
  pilTrnAddKey("SkyMaskYord", "ESO PRO TEL DIS YORD", "%d", 
               "Order Y Sky to Mask");
  pilTrnAddKey("SkyMaskYrms", "ESO PRO TEL DIS YRMS", "%.2f",
               "rms error in Y");
  pilTrnAddKey("SkyMaskXrms", "ESO PRO TEL DIS XRMS", "%.2f",
               "rms error in X");
  pilTrnAddKey("SkyMaskX", "ESO PRO TEL DIS X_%d_%d", "%.2f",
               "X Sky to mask");
  pilTrnAddKey("SkyMaskY", "ESO PRO TEL DIS Y_%d_%d", "%.2f",
               "Y Sky to mask");
  pilTrnAddKey("MaskPlateScale", "ESO PRO PLATSCAL", "%f", 
               "plate scale of telescope in mask plane (arcsec/mm)");
  pilTrnAddKey("MaskOptAxisX", "ESO PRO TEL X", "%f",
               "X position of Tel optical axis in the quadrant "
               "coordinate system (mm)");
  pilTrnAddKey("MaskOptAxisY", "ESO PRO TEL Y", "%f",
               "Y position of Tel optical axis in the quadrant "
               "coordinate system (mm)");
  pilTrnAddKey("TabPixelScale", "ESO PRO CDELT", "%f",
               "Pixel Scale in FilTab");
  pilTrnAddKey("TabCrpix", "ESO PRO CRPIX%d", "%f", 
               "CRPIXi in Filter/Stmc Tables");
  /**************************************************************
    Add keyword "ShiftCrpix" in Filter Table. This keys is 
    not listed in the FDR doc, but is necessary to reset center of images
  *************************************************************/
  pilTrnAddKey("ShiftCrpix", "ESO PRO SHIFT CRPIX%d", "%f", 
               "Shift of CRPIXi needed to reset image center");
  pilTrnAddKey("TabCrval", "ESO PRO CRVAL%d", "%f", 
               "CRVALi in Filter/Stmc Tables");
  pilTrnAddKey("TabCDMatrix", "ESO PRO CD%d_%d", "%.2f",
               "CD Matrix in Filter/Stmc Tables");
  pilTrnAddKey("PiezoVolt", "ESO PRO PIEZO%d VOLT", "%f",
               "Piezo Voltage");
  pilTrnAddKey("NoPiezo", "ESO PRO PIEZO NO", "%d",
               "No of piezos");
  pilTrnAddKey("RotAngle", "ESO PRO ROT", "%f",
               "Rotator Angle");
  pilTrnAddKey("OffX", "ESO PRO OFF X", "%f",
               "Average offset in X (pixel) of Mask Image");
  pilTrnAddKey("OffY", "ESO PRO OFF Y", "%f",
               "Average offset in Y (pixel) of Mask Image");
  pilTrnAddKey("HistoBinNo", "ESO PRO BIN NO", "%d",
               "Number of bins of histogram - recipe VmMskChk");
  pilTrnAddKey("HistoBinWidth", "ESO PRO BIN WIDTH", "%f", 
               "width of bins of histogram - recipe VmMskChk");
  pilTrnAddKey("HistoBinStart", "ESO PRO BIN START", "%f", 
               "value of first bin - rec VmMskChk");
  pilTrnAddKey("HistoBinX", "ESO PRO BIN X_%d", "%d",
               "Number in bin i for offset in X - recipe VmMskChk");
  pilTrnAddKey("HistoBinY", "ESO PRO BIN Y_%d", "%d",
               "Number in bin i for offset in Y - recipe VmMskChk");
  pilTrnAddKey("FlatImgRMS", "ESO PRO FLAT RMS", "%f",
               "RMS of normalized master flat image");
  pilTrnAddKey("DarkOffset", "ESO PRO DARK DIFF", "%f",
               "Deviation from nominal dark level");
  pilTrnAddKey("BiasOffset", "ESO PRO BIAS DIFF", "%f",
               "Offset from nominal bias level");
  pilTrnAddKey("Seeing", "ESO PRO FWHM", "%f",
               "Estimated seeing in arcsec");
  pilTrnAddKey("SkyLevel", "ESO PRO SKY", "%f", 
               "Estimated sky brigthness in mag/arcsec2");
  pilTrnAddKey("LimitPoint", "ESO PRO LIMIT POINT", "%f", 
               "Detection limit for point sources in mag");
  pilTrnAddKey("LimitAperture", "ESO PRO LIMIT APERT", "%f", 
               "Detection limit for 3arcsec aperture in mag");
  pilTrnAddKey("MagZero", "ESO PRO MAG ZERO", "%f",
               "Zero point magnitude");
  pilTrnAddKey("MagZeroRms", "ESO PRO MAGZERO RMS", "%f",
               "Zero point magnitude standard deviation");
  pilTrnAddKey("ColorTerm", "ESO PRO COLTERM", "%f",
               "Color term for filter");
  pilTrnAddKey("ColorTermRms", "ESO PRO COLTERM RMS", "%f",
               "Color term standard deviation");
  pilTrnAddKey("Extinction", "ESO PRO EXTINCT", "%f",
               "Atmospheric extinction");
  pilTrnAddKey("ExtinctionRms", "ESO PRO EXTINC RMS", "%f",
               "Atmospheric extinction standard deviation");
  pilTrnAddKey("AirMass", "ESO PRO AIRMASS", "%f", "Software airmass");
  pilTrnAddKey("Airmass", "AIRMASS", "%f", "Airmass");
  pilTrnAddKey("WlenStart", "ESO PRO WLEN START", "%f",
               "Start wave for filter");
  pilTrnAddKey("WlenEnd", "ESO PRO WLEN END", "%f",
               "End wave for filter");
  pilTrnAddKey("WlenCen", "ESO PRO WLEN CEN", "%f",
               "Central wave of spectrum: Lambda zero");
  pilTrnAddKey("WlenInc", "ESO PRO WLEN INC", "%f",
               "Increment in spectrum");
  pilTrnAddKey("SphotOrder", "ESO PRO PHO MOS ORD", "%d",
               "Order of the spectrophotometric sensitivity function "
               "polynomial fit");
  pilTrnAddKey("SphotModel", "ESO PRO PHO MOS_%d", "%f",
               "Spectrophotometric sensitivity function polynomial fit "
               "coefficient");
  pilTrnAddKey("LampName", "ESO INS LAMP%d NAME", "%s", "Lamp Name");
  pilTrnAddKey("LampState", "ESO INS LAMP%d STATE", "%s", "Lamp State");
  pilTrnAddKey("LampTime", "ESO INS LAMP%d TIME", "%d", "Lamp Time");

  pilTrnAddKey("MatchNstars", "ESO PRO MATCH NSTARS", "%d", 
               "Number of matched stars");
  pilTrnAddKey("MatrixFittedPoints", "ESO PRO MATRIX NFIT", "%d", 
               "Number of points used for the matrix fit");
  
 /*
  *  PAF entries
  */
  pilTrnAddKey("PafHeaderStart", "PAF.HDR.START", "%s",
               "PAF header start");
  pilTrnAddKey("PafType", "PAF.TYPE", "%s",
               "PAF type");
  pilTrnAddKey("PafId", "PAF.ID", "%s",
               "PAF ID");
  pilTrnAddKey("PafName", "PAF.NAME", "%s",
               "PAF NAME");
  pilTrnAddKey("PafDesc", "PAF.DESC", "%s",
               "PAF DESC");
  pilTrnAddKey("PafCrteName", "PAF.CRTE.NAME", "%s",
               "PAF CRTE name");
  pilTrnAddKey("PafCrteDaytim", "PAF.CRTE.DAYTIM", "%s",
               "PAF CRTE day time");
  pilTrnAddKey("PafLchgName", "PAF.LCHG.NAME", "%s",
               "PAF LCHG name");
  pilTrnAddKey("PafLchgDaytim", "PAF.LCHG.DAYTIM", "%s",
               "PAF LCHG day time");
  pilTrnAddKey("PafChckName", "PAF.CHCK.NAME", "%s",
               "PAF CHCK name");
  pilTrnAddKey("PafChckDaytim", "PAF.CHCK.DAYTIM", "%s",
               "PAF CHCK day time");
  pilTrnAddKey("PafChecksum", "PAF.CHCK.CHECKSUM", "%s",
               "PAF Checksum");
  pilTrnAddKey("PafHeaderEnd", "PAF.HDR.END", "%s",
               "PAF header end");

  pilTrnAddKey("PAFMaskCcdX", "PRO.MASK.CCD.X_%d_%d", "%s", 
               "Mask to CCD coordinate model coefficient");
  pilTrnAddKey("PAFMaskCcdY", "PRO.MASK.CCD.Y_%d_%d", "%s", 
               "Mask to CCD coordinate model coefficient");
  pilTrnAddKey("PAFMaskCcdXord", "PRO.MASK.CCD.XORD", "%s", 
               "Mask to CCD coordinate model coefficient");  
  pilTrnAddKey("PAFMaskCcdYord", "PRO.MASK.CCD.YORD", "%s", 
               "Mask to CCD coordinate model coefficient");  
  pilTrnAddKey("PAFMaskCcdXrms", "PRO.MASK.CCD.XRMS", "%s", 
               "Mask to CCD coordinate model coefficient");  
  pilTrnAddKey("PAFMaskCcdYrms", "PRO.MASK.CCD.YRMS", "%s", 
               "Mask to CCD coordinate model coefficient");  
  pilTrnAddKey("PAFMaskCcdX0", "PRO.MASK.CCD.X0", "%s", 
               "Mask to CCD coordinate model coefficient");
  pilTrnAddKey("PAFMaskCcdY0", "PRO.MASK.CCD.Y0", "%s", 
               "Mask to CCD coordinate model coefficient");
  pilTrnAddKey("PAFMaskCcdXX", "PRO.MASK.CCD.XX", "%s", 
               "Mask to CCD coordinate model coefficient");
  pilTrnAddKey("PAFMaskCcdYY", "PRO.MASK.CCD.YY", "%s", 
               "Mask to CCD coordinate model coefficient");
  pilTrnAddKey("PAFMaskCcdXY", "PRO.MASK.CCD.XY", "%s", 
               "Mask to CCD coordinate model coefficient");
  pilTrnAddKey("PAFMaskCcdYX", "PRO.MASK.CCD.YX", "%s", 
               "Mask to CCD coordinate model coefficient");

  pilTrnAddKey("PAFCcdMaskDate","PRO.CCD.MASK.DAYTIM", "%s", 
               "Mask to CCD Date");
  pilTrnAddKey("PAFCcdMaskTemp","PRO.CCD.MASK.TEMP", "%s", 
               "Mask to CCD Temperature");

  pilTrnAddKey("PAFCcdMaskX0", "PRO.CCD.MASK.X0", "%s", 
               "CCD to Mask coordinate model coefficient");
  pilTrnAddKey("PAFCcdMaskY0", "PRO.CCD.MASK.Y0", "%s", 
               "CCD to Mask coordinate model coefficient");
  pilTrnAddKey("PAFCcdMaskXX", "PRO.CCD.MASK.XX", "%s", 
               "CCD to Mask coordinate model coefficient");
  pilTrnAddKey("PAFCcdMaskYY", "PRO.CCD.MASK.YY", "%s", 
               "CCD to Mask coordinate model coefficient");
  pilTrnAddKey("PAFCcdMaskXY", "PRO.CCD.MASK.XY", "%s", 
               "CCD to Mask coordinate model coefficient");
  pilTrnAddKey("PAFCcdMaskYX", "PRO.CCD.MASK.YX", "%s", 
               "CCD to Mask coordinate model coefficient");
  pilTrnAddKey("PAFCcdMaskX", "PRO.CCD.MASK.X_%d_%d", "%s", 
               "CCD to Mask coordinate model coefficient");
  pilTrnAddKey("PAFCcdMaskY", "PRO.CCD.MASK.Y_%d_%d", "%s", 
               "CCD to Mask coordinate model coefficient");
  pilTrnAddKey("PAFCcdMaskXord", "PRO.CCD.MASK.XORD", "%s", 
               "CCD to Mask coordinate model X order");  
  pilTrnAddKey("PAFCcdMaskYord", "PRO.CCD.MASK.YORD", "%s", 
               "CCD to Mask coordinate model Y order");  
  pilTrnAddKey("PAFCcdMaskXrms", "PRO.CCD.MASK.XRMS", "%s", 
               "CCD to Mask coordinate X model rms");  
  pilTrnAddKey("PAFCcdMaskYrms", "PRO.CCD.MASK.YRMS", "%s", 
               "CCD to Mask coordinate Y model rms");  

  pilTrnAddKey("PAFCcdSkyDate", "PRO.CCD.SKY.DAYTIM", "%s",
               "Time stamp for CCD to Sky matrix");
  pilTrnAddKey("PAFCcdSkyTemp", "PRO.CCD.SKY.TEMP", "%s",
               "Temperature CCD to Sky matrix");
  pilTrnAddKey("PAFCcdSkyXord", "PRO.CCD.SKY.XORD", "%s", 
               "Order in X of the CCD to Sky polynomial");  
  pilTrnAddKey("PAFCcdSkyYord", "PRO.CCD.SKY.YORD", "%s", 
               "Order in Y of the CCD to Sky polynomial");  
  pilTrnAddKey("PAFCcdSkyX", "PRO.CCD.SKY.X_%d_%d", "%s", 
               "X CCD to Sky matrix coefficient");
  pilTrnAddKey("PAFCcdSkyY", "PRO.CCD.SKY.Y_%d_%d", "%s", 
               "Y CCD to Sky matrix coefficient");
  pilTrnAddKey("PAFCcdSkyXrms", "PRO.CCD.SKY.XRMS", "%s", 
               "RMS of CCD to Sky model fit in X");  
  pilTrnAddKey("PAFCcdSkyYrms", "PRO.CCD.SKY.YRMS", "%s", 
               "RMS of CCD to Sky model fit in Y");  
  pilTrnAddKey("PAFSkyCcdXord", "PRO.SKY.CCD.XORD", "%s", 
               "Order in X of the Sky to CCD polynomial");  
  pilTrnAddKey("PAFSkyCcdYord", "PRO.SKY.CCD.YORD", "%s", 
               "Order in Y of the Sky to CCD polynomial");  
  pilTrnAddKey("PAFSkyCcdX", "PRO.SKY.CCD.X_%d_%d", "%s", 
               "X Sky to CCD matrix coefficient");
  pilTrnAddKey("PAFSkyCcdY", "PRO.SKY.CCD.Y_%d_%d", "%s", 
               "Y Sky to CCD matrix coefficient");
  pilTrnAddKey("PAFSkyCcdXrms", "PRO.SKY.CCD.XRMS", "%s", 
               "RMS of Sky to CCD model fit in X");  
  pilTrnAddKey("PAFSkyCcdYrms", "PRO.SKY.CCD.YRMS", "%s", 
               "RMS of Sky to CCD model fit in Y");  

  pilTrnAddKey("PAFCrvOptDate","PRO.CRV.OPT.DAYTIM", "%s", 
               "Spectral curvature and optical distorsion models date");
  pilTrnAddKey("PAFCrvOptTemp","PRO.CRV.OPT.TEMP", "%s", 
               "Instrument temperature");

  pilTrnAddKey("PAFIdsDate","PRO.IDS.MAT.DAYTIM", "%s", 
               "Inverse Dispersion Relation date");
  pilTrnAddKey("PAFIdsTemp","PRO.IDS.MAT.TEMP", "%s", 
               "Inverse Dispersion Relation temperature");
  pilTrnAddKey("PAFIdsModOrd", "PRO.IDS.REL.ORD", "%s", 
               "Inverse Dispersion model order");  
  pilTrnAddKey("PAFIdsModXord", "PRO.IDS.MAT.XORD", "%s", 
               "Inverse Dispersion X model order");  
  pilTrnAddKey("PAFIdsModYord", "PRO.IDS.MAT.YORD", "%s", 
               "Inverse Dispersion Y model order");  
  pilTrnAddKey("PAFIdsMod", "PRO.IDS.MAT_%d_%d_%d", "%s", 
               "Inverse Dispersion model coefficient");
  pilTrnAddKey("PAFIdsXrms", "PRO.IDS.MAT.X_RMS", "%s", 
               "RMS of inverse dispersion model fit in X");  
  pilTrnAddKey("PAFIdsYrms", "PRO.IDS.MAT.Y_RMS", "%s", 
               "RMS of inverse dispersion model fit in Y");  

  pilTrnAddKey("PAFOptDate","PRO.OPT.DIS.DAYTIM", "%s", 
               "Optical distorsion date");
  pilTrnAddKey("PAFOptTemp","PRO.OPT.DIS.TEMP", "%s", 
               "Instrument temperature");
  pilTrnAddKey("PAFOptDisXord", "PRO.OPT.DIS.XORD", "%s", 
               "Optical X distorsion model order");  
  pilTrnAddKey("PAFOptDisYord", "PRO.OPT.DIS.YORD", "%s", 
               "Optical Y distorsion model order");  
  pilTrnAddKey("PAFOptDisX", "PRO.OPT.DIS.X_%d_%d", "%s", 
               "Optical X distorsion model coefficient");
  pilTrnAddKey("PAFOptDisY", "PRO.OPT.DIS.Y_%d_%d", "%s", 
               "Optical Y distorsion model coefficient");
  pilTrnAddKey("PAFOptDisXrms", "PRO.OPT.DIS.X_RMS", "%s", 
               "RMS of optical distortion model fit in X");  
  pilTrnAddKey("PAFOptDisYrms", "PRO.OPT.DIS.Y_RMS", "%s", 
               "RMS of optical distortion model fit in Y");  

  pilTrnAddKey("PAFCrvModOrd", "PRO.CRV.POL.ORD", "%s", 
               "Curvature model order");  
  pilTrnAddKey("PAFCrvModXord", "PRO.CRV.MOD.XORD", "%s", 
               "Curvature X model order");  
  pilTrnAddKey("PAFCrvModYord", "PRO.CRV.MOD.YORD", "%s", 
               "Curvature Y model order");  
  pilTrnAddKey("PAFCrvMod", "PRO.CRV.MOD_%d_%d_%d", "%s", 
               "Curvature model coefficient");
  pilTrnAddKey("PAFCrvXrms", "PRO.CRV.MOD.X_RMS", "%s", 
               "RMS of curvature model fit in X");  
  pilTrnAddKey("PAFCrvYrms", "PRO.CRV.MOD.Y_RMS", "%s", 
               "RMS of curvature model fit in Y");  

  pilTrnAddKey("PAFZeroXord", "PRO.ZERO.XORD", "%s", 
               "Contamination model X order");  
  pilTrnAddKey("PAFZeroYord", "PRO.ZERO.YORD", "%s", 
               "Contamination model Y order");  
  pilTrnAddKey("PAFZeroX", "PRO.ZERO.X_%d_%d", "%s", 
               "X contamination model coefficient");
  pilTrnAddKey("PAFZeroY", "PRO.ZERO.Y_%d_%d", "%s", 
               "Y contamination model coefficient");
  pilTrnAddKey("PAFZeroXrms", "PRO.ZERO.X_RMS", "%s", 
               "RMS of contamination model fit in X");  
  pilTrnAddKey("PAFZeroYrms", "PRO.ZERO.Y_RMS", "%s", 
               "RMS of contamination model fit in Y");  
  pilTrnAddKey("PAFZeroWidth", "PRO.ZERO.PWIDTH", "%s", 
               "Width of the Zero Order contamination");

  return EXIT_SUCCESS;
}


/**
 * @brief
 *   Clear the keyword mapping.
 *
 * @return The function returns 0 on success and a non-zero value
 *    otherwise.
 *
 * The function removes a entries from the recipe keyword map. After
 * calling this function the recipe's keyword map does not exist any
 * longer. Before keywords are accessed again the recipe keyword map
 * must be (re-)initialized using the function pilTrnInitKeywordMap().
 */

int
pilTrnClearKeywordMap(void)
{

    if (keymap != NULL) {
        deletePilKeymap(keymap);
        keymap = NULL;
    }

    return 0;

}


/**
 * @brief
 *   Load keyword mapping from file.
 *
 * @param filename Name of file containing keywords aliases definition.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * This function would be called if a configuration file defining
 * the keyword mapping were found. Just one keyword mapping definition
 * is allowed in one pipeline run, an error condition will occur
 * otherwise.
 */

int pilTrnLoadKeywordMap(const char *filename)
{
  char  modName[] = "pilTrnLoadKeywordMap";
  char  line[PIL_LINE_LENGTH_MAX];
  char  alias[PIL_LINE_LENGTH_MAX];
  char  name[PIL_LINE_LENGTH_MAX];
  char  form[PIL_LINE_LENGTH_MAX];
  char  comment[PIL_LINE_LENGTH_MAX];
  int   nameFlag, formFlag, commentFlag, aliasFlag;
  FILE *fp;

  if ((fp = fopen (filename, "r")) == NULL) {
    pilMsgWarning(modName, 
                  "Problems in opening keyword map file %s", filename);
    if (keymap == NULL) {
      pilMsgWarning(modName, "No default keyword map was loaded");
      return EXIT_FAILURE;
    } 
    else {
      pilMsgWarning(modName, "Using default keyword mapping only");
      return EXIT_SUCCESS;
    }
  }

  if (keymap == NULL) {
    pilMsgWarning(modName, 
      "No default keynames mapping loaded: "
                  "relying just on mapping from file %s", filename);
    keymap = newPilKeymap();
  }

  nameFlag = 0;
  formFlag = 0;
  commentFlag = 0;
  aliasFlag = 0;
  
  while (fgets(line, PIL_LINE_LENGTH_MAX, fp)) {

   /* 
    * Empty line: assume keyword definition is complete
    */

    if (strempty(line, COMMENT_CHARS)) {
      pilMsgDebug(modName,"Empty line");
      
     /* 
      * Add key to table if all needed items are provided
      */

      if (nameFlag && formFlag && commentFlag && aliasFlag) {
        if (pilTrnAddKey(alias, name, form, comment) == EXIT_FAILURE) {
          fclose(fp);
          return EXIT_FAILURE;
        }
        pilMsgDebug(modName, "Alias '%s' added to keyword map\n", alias);
      }
      else if (nameFlag || formFlag || commentFlag || aliasFlag)
        pilMsgWarning(modName, "A keyword definition in keyword map "
                      "file %s is incomplete", filename);
 
      nameFlag = 0;
      formFlag = 0;
      commentFlag = 0;
      aliasFlag = 0;
      continue;
    }

    if (sscanf(line, KEYMAP_NAME "%[^\n]\n", name)) {
      strtrim(name,2);
      nameFlag = 1;
      pilMsgDebug(modName, "Name: %s\n", name);
      continue;
    }

    if (sscanf(line, KEYMAP_FORM "%[^\n]\n", form)) {
      strtrim(form,2);
      formFlag = 1;
      pilMsgDebug(modName, "Form: %s\n", form);
      continue;
    }

    if (sscanf(line, KEYMAP_COMM "%[^\n]\n", comment)) {
      strtrim(comment,2);
      commentFlag = 1;
      pilMsgDebug(modName, "Comment: %s\n", comment);
      continue;
    }

    if (sscanf(line, KEYMAP_ALIAS "%[^\n]\n", alias)) {
      strtrim(alias,2);
      aliasFlag = 1;
      pilMsgDebug(modName, "Alias: %s\n", alias);
      continue;
    }

  }
   
  fclose(fp);

 /* 
  * Add last key to table if all needed items are provided
  */

  if (nameFlag && formFlag && commentFlag && aliasFlag) {
    if ((pilTrnAddKey(alias, name, form, comment) == EXIT_FAILURE)) {
      return EXIT_FAILURE;
    }
    pilMsgDebug(modName, "Alias '%s' added to keyword map\n", alias);
  }
  else if (nameFlag || formFlag || commentFlag || aliasFlag)
    pilMsgWarning(modName, 
      "A keyword definition in keyword map file %s is incomplete", filename);

  return EXIT_SUCCESS;
}

/**
 * @brief
 *   Translate a keyword name into its value, inserting indexes
 *   where requested.
 *
 * @param name Keyword name.
 * @param ... Variable argument list of integer indexes.
 *
 * @return Keyword translation.
 *
 * The number of specified integer indexes must agree with the
 * number of %d-like formats in the keyword value string.
 * The keyword translation is a character string that is
 * newly allocated (not a pointer to the dictionary component)
 * and it is responsibility of the caller to deallocate it
 * if and when it is found appropriate. Ideally, the translation
 * of the same keyword should be requested just once, and kept
 * in the calling routine.
 */

const char *pilKeyTranslate(const char *name, ...)
{
  const char  modName[] = "pilKeyTranslate";
  const char *value = pilKeymapGetValue(keymap, name);
  const char *p;
  char       *translation;
  va_list     indexes;
  int         ind;
  int         n = 0;
  int         i = 0;
  int         j = 0;
  int         l = 0;

  if (value == NULL) {
    pilMsgError(modName, "Translation of alias %s not found", name);
    return NULL;
  }

  va_start(indexes, name);

  p = value;

 /*
  * l is the original length of the string: it will be modified to
  * accomodate the indexes according to the specified "%d" formats.
  */
  l = strlen(p);
  while ((p = strstr(p,"%"))) {
   /*
    * n = number found within %d format (set to ZERO if no number found).
    * i = number of characters occupied by %d format specification.
    * p = pointer moving along the character string.
    */
    for (n=0, i=2, p++; *p != 'd'; i++, p++) {
      if (*p >= '0' && *p <= '9') {
        n = atoi(p);
        i += strstr(p,"d") - p;
        break;
      }
    }
   /*
    *  Get integer value to be displayed with current %d format, and find
    *  j (its number of digits).
    */
    ind = va_arg(indexes,int);

    if (ind > 0) {
      for (j=0; ind; ind /= 10, j++);
    }
    else if (ind == 0) {
      j = 1;
    }
    else {
      return NULL;
    }

   /*
    *  The space required to display the current index is at least its 
    *  own number of digits, or more if the %d format specifies more.
    */
    n = MAX(n,j);
   /*
    *  The length of the string must be modified. i spaces were already
    *  taken by the %d format specification, n spaces are actually required
    *  by the current index.
    */
    l += n-i;
  }
  l++;
  translation = (char *)pil_malloc(l * sizeof(char));
  va_end(indexes);
 /*
  *  Reset list of indexes, and apply translation.
  */
  va_start(indexes, name);
  n = vsprintf(translation,value,indexes);

  va_end(indexes);

  return translation;

}


/**
 * @brief
 *   Translate a keyword alias into its value, inserting numerical indexes
 *   where requested.
 *
 * @param name Keyword alias.
 * @param ...  Variable argument list of integer indexes.
 *
 * @return Keyword translation if no error occurred, otherwise the 
 *   function returns @c NULL.
 *
 * The number of specified integer indexes must agree with the
 * number of %d-like formats in the keyword value string.
 * The keyword translation is a character string that is
 * statically allocated (not a pointer to the dictionary 
 * component), which is overwritten every time a new translation
 * is requested.
 */

const char *pilTrnGetKeyword(const char *alias, ...)
{
  const char  modName[] = "pilTrnGetKeyword";

  static char translation[MAX_TRN_LENGTH];
  const char *value;
  const char *p;
  va_list     indexes;
  int         ind;
  int         n = 0;
  int         i = 0;
  int         j = 0;
  int         l = 0;

  if (!(value = pilKeymapGetValue(keymap, alias))) {
    pilMsgError(modName, "Translation of alias %s not found", alias);
    return NULL;
  }

  va_start(indexes, alias);

  p = value;

 /*
  * l is the original length of the string: it will be modified
  * to count also the space taken by the indexes according to the 
  * specified "%d" formats.
  */

  l = strlen(p);

  while ((p = strstr(p,"%"))) {

   /*
    * n = number found within %d format (set to ZERO if no number found).
    * i = number of characters occupied by %d format specification.
    * p = pointer moving along the character string.
    */

    for (n = 0, i = 2, p++; *p != 'd'; i++, p++) {
      if (*p >= '0' && *p <= '9') {
        n = atoi(p);
        i += strstr(p, "d") - p;
        break;
      }
    }

   /*
    *  Get integer value to be displayed with current %d format, and find
    *  j (its number of digits).
    */

    ind = va_arg(indexes,int);

    if (ind > 0) {
      for (j = 0; ind; ind /= 10, j++);
    }
    else if (ind == 0) {
      j = 1;
    }
    else {
      return NULL;
    }

   /*
    *  The space required to display the current index is at least its 
    *  own number of digits, or more if the %d format specifies more.
    */

    n = MAX(n, j);

   /*
    *  The length of the string must be modified. i spaces were already
    *  taken by the %d format specification, n spaces are actually required
    *  by the current index.
    */

    l += n - i;
  }

  l++;                         /* Space for the terminating null char      */

  if (l > MAX_TRN_LENGTH)      /* Not enough space for holding translation */
    return NULL;

  va_end(indexes);

 /*
  *  Reset list of indexes, and apply translation.
  */

  va_start(indexes, alias);

  n = vsprintf(translation, value, indexes);

  va_end(indexes);

  return translation;

}


/**
 * @brief
 *   Get comment string of an alias translation.
 *
 * @param alias  Alias string.
 *
 * @return The function returns the reference to the comment string if the
 *   alias was found, if it was not foun, or in case of an error the return
 *   value is @c NULL.
 *
 * The function searches the internal keyword map for the alias string
 * @em alias. If the alias string is found and if there is a comment
 * associated to its translation the reference to the comment is returned.
 * If the alias string is not found or there is no comment available
 * the function returns @c NULL.
 *
 * @note
 *   The comment string must not be modified through the returned reference.
 */

const char *pilTrnGetComment(const char *alias)
{

    return pilKeymapGetComment(keymap, alias);
    
}


/**
 * @brief
 *   Add to the category mapping a category with given name and value.
 *
 * @param name Category name.
 * @param value Category value.
 *
 * @return @c EXIT_SUCCESS on success, @c EXIT_FAILURE on failure.
 *
 * The function creates and adds the new category to the mapping dictionary.
 */

int pilTrnAddCategory(const char *name, const char *value)
{

  PilCategory *category;

  if ((category = pilCatmapLookup(catmap, name))) {
    return pilCatSetValue(category, value);
  }
  else {
    category = newPilCategory(name, value);
    return pilCatmapInsert(catmap, category);
  }

}

/**
 * @brief
 *   Create a default category mapping dictionary.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * This function would be called if no configuration file defining
 * the category mapping were found. Just one category mapping definition
 * is allowed in one pipeline run, an error condition will occur
 * otherwise.
 */

int
pilTrnInitCategoryMap(void)
{
  char modName[] = "pilTrnInitCategoryMap";

  if (catmap == NULL) {
    catmap = newPilCatmap();
  }
  else {
    pilMsgWarning(modName, 
    "Double definition of categories - the first one is taken");
    return EXIT_FAILURE;
  }

  /*
   * General supported frame categories
   */

  pilTrnAddCategory("BadPixelMap", "BAD_PIXEL_MAP");
  pilTrnAddCategory("DetectorProperties", "DETECTOR_PROPERTIES");

  pilTrnAddCategory("Bias", "BIAS");
  pilTrnAddCategory("MasterBias", "MASTER_BIAS");
  pilTrnAddCategory("Dark", "DARK");
  pilTrnAddCategory("MasterDark", "MASTER_DARK");

  pilTrnAddCategory("ImgScreenFlat", "IMG_SCREEN_FLAT");
  pilTrnAddCategory("ImgMasterScreenFlat", "IMG_MASTER_SCREEN_FLAT");
  pilTrnAddCategory("ImgCombScreenFlat", "IMG_COMBINED_SCREEN_FLAT");
  pilTrnAddCategory("ImgSkyFlat", "IMG_SKY_FLAT");
  pilTrnAddCategory("ImgMasterSkyFlat", "IMG_MASTER_SKY_FLAT");

  pilTrnAddCategory("ImgPreimaging", "IMG_PREIMAGING");
  pilTrnAddCategory("ImgScience", "IMG_SCIENCE");
  pilTrnAddCategory("ImgStandard", "IMG_STANDARD");
  pilTrnAddCategory("ImgAstrometry", "IMG_ASTROMETRY");
  pilTrnAddCategory("ImgFringes", "IMG_FRINGES");

  pilTrnAddCategory("redImgScience", "IMG_SCIENCE_REDUCED");
  pilTrnAddCategory("redImgStandard", "IMG_STANDARD_REDUCED");
  pilTrnAddCategory("redImgAstrometry","IMG_ASTROMETRY_REDUCED");
  
  pilTrnAddCategory("ImgStarMatchAstrometry","IMG_ASTROMETRY_STAR_MATCH");

  pilTrnAddCategory("MosScreenFlat", "MOS_SCREEN_FLAT");
  pilTrnAddCategory("MosMasterScreenFlat", "MOS_MASTER_SCREEN_FLAT");
  pilTrnAddCategory("MosCombScreenFlat", "MOS_COMBINED_SCREEN_FLAT");
  pilTrnAddCategory("MosArcSpectrum", "MOS_ARC_SPECTRUM");
  pilTrnAddCategory("MosArcSpectrumExtracted", "MOS_ARC_SPECTRUM_EXTRACTED");
  pilTrnAddCategory("MosScience", "MOS_SCIENCE");
  pilTrnAddCategory("MosScienceFluxReduced", "MOS_SCIENCE_FLUX_REDUCED");
  pilTrnAddCategory("MosScienceReduced", "MOS_SCIENCE_REDUCED");
  pilTrnAddCategory("MosScienceExtracted", "MOS_SCIENCE_EXTRACTED");
  pilTrnAddCategory("MosScienceSky", "MOS_SCIENCE_SKY");
  pilTrnAddCategory("MosSkyReduced", "MOS_SKY_REDUCED");
  pilTrnAddCategory("MosFringesSky", "MOS_FRINGES_SKY");
  pilTrnAddCategory("MosFringes", "MOS_FRINGES");
  pilTrnAddCategory("MosStandard", "MOS_STANDARD");
  pilTrnAddCategory("MosStandardReduced", "MOS_STANDARD_REDUCED");
  pilTrnAddCategory("MosStandardExtracted", "MOS_STANDARD_EXTRACTED");
  pilTrnAddCategory("MosStandardSkyReduced", "MOS_STANDARD_SKY_EXTRACTED");
  pilTrnAddCategory("MosStandardSky", "MOS_STANDARD_SKY");

  pilTrnAddCategory("IfuScreenFlat", "IFU_SCREEN_FLAT");
  pilTrnAddCategory("IfuMasterScreenFlat", "IFU_MASTER_SCREEN_FLAT");
  pilTrnAddCategory("IfuArcSpectrum", "IFU_ARC_SPECTRUM");
  pilTrnAddCategory("IfuArcSpectrumExtracted", "IFU_ARC_SPECTRUM_EXTRACTED");
  pilTrnAddCategory("IfuFlatSpectrumExtracted", "IFU_FLAT_SPECTRUM_EXTRACTED");
  pilTrnAddCategory("IfuFov", "IFU_FOV");
  pilTrnAddCategory("IfuStdFov", "IFU_STD_FOV");
  pilTrnAddCategory("IfuFullFov", "IFU_FULL_FOV");
  pilTrnAddCategory("IfuFullStdFov", "IFU_FULL_STD_FOV");
  pilTrnAddCategory("IfuIds", "IFU_IDS");
  pilTrnAddCategory("IfuTrace", "IFU_TRACE");
  pilTrnAddCategory("IfuTransmission", "IFU_TRANSMISSION");
  pilTrnAddCategory("IfuIdent", "IFU_IDENT");
  pilTrnAddCategory("IfuScience", "IFU_SCIENCE");
  pilTrnAddCategory("IfuScienceReduced", "IFU_SCIENCE_REDUCED");
  pilTrnAddCategory("IfuScienceFluxReduced", "IFU_SCIENCE_FLUX_REDUCED");
  pilTrnAddCategory("IfuStandard", "IFU_STANDARD");
  pilTrnAddCategory("IfuStandardReduced", "IFU_STANDARD_REDUCED");
  pilTrnAddCategory("IfuStandardExtracted", "IFU_STANDARD_EXTRACTED");
  pilTrnAddCategory("IfuScienceSky", "IFU_SCIENCE_SKY");

  pilTrnAddCategory("LineCatalog", "LINE_CATALOG");
  pilTrnAddCategory("AtmosphericExtinction", "ATMOSPHERIC_EXTINCTION");

  /* FIXME:
   *   VIMOS specific frame categories should go to a instrument
   *   specific module. Only instrument independent categories should
   *   go here.
   */

  pilTrnAddCategory("WindowTable", "WINDOW_TABLE");
  pilTrnAddCategory("ObjectTable", "OBJECT_TABLE");
  pilTrnAddCategory("CcdTable", "CCD_TABLE");
  pilTrnAddCategory("GrismTable", "GRISM_TABLE");
  pilTrnAddCategory("ExtractTable", "EXTRACT_TABLE");
  pilTrnAddCategory("SphotTable", "SPECPHOT_TABLE"); 
  pilTrnAddCategory("MosSphotTable", "MOS_SPECPHOT_TABLE"); 
  pilTrnAddCategory("IfuSphotTable", "IFU_SPECPHOT_TABLE"); 
  pilTrnAddCategory("StdFluxTable", "STD_FLUX_TABLE"); 
  pilTrnAddCategory("ExtinctTable", "EXTINCT_TABLE"); 
  pilTrnAddCategory("IfuTable", "IFU_TABLE");
  pilTrnAddCategory("FilterTable", "FILTER_TABLE");
  pilTrnAddCategory("TelescopeTable", "TELESCOPE_TABLE");
  pilTrnAddCategory("PhotometricTable", "PHOTOMETRIC_TABLE");
  pilTrnAddCategory("PhotometricCoeffTable", "PHOT_COEFF_TABLE");
  pilTrnAddCategory("PhotometricCatalog", "PHOTOMETRIC_CATALOG");
  pilTrnAddCategory("GalaxyTable", "IMG_GALAXY_TABLE");
  pilTrnAddCategory("StarTable", "IMG_STAR_TABLE");
  pilTrnAddCategory("AstrometricTable", "ASTROMETRIC_TABLE");
  pilTrnAddCategory("StarMatchTable", "IMG_STAR_MATCH_TABLE");
  pilTrnAddCategory("FlexureTable", "FLEXURES_TABLE");
  pilTrnAddCategory("PiezoTable", "PIEZO_TABLE");
  pilTrnAddCategory("GridMaskImage", "MASK_TO_CCD");
  pilTrnAddCategory("AstroMaskImage", "MASK_COORDINATES");
  pilTrnAddCategory("FlexureCompensation", "FLEXURE_COMPENSATION");
  pilTrnAddCategory("InstrumentFlexure", "INSTRUMENT_FLEXURE");
  pilTrnAddCategory("ImgScienceReducedSequence", 
                    "IMG_SCIENCE_REDUCED_SEQUENCE");
  pilTrnAddCategory("MosScienceReducedSequence", 
                    "MOS_SCIENCE_REDUCED_SEQUENCE");
  pilTrnAddCategory("IfuScienceReducedSequence",
                    "IFU_SCIENCE_REDUCED_SEQUENCE");
  pilTrnAddCategory("Stack2dSpectra", "STACK_2D_SPECTRA");
  pilTrnAddCategory("Stack1dSpectra", "STACK_1D_SPECTRA");
  pilTrnAddCategory("MosZeroOrder", "MOS_ZERO_ORDER");

 /*
  * PAF files generic category
  */

  pilTrnAddCategory("PAFCategory", "PAF");

  return EXIT_SUCCESS;

}


/**
 * @brief
 *   Clear the frame category mapping.
 *
 * @return The function returns 0 on success and a non-zero value
 *    otherwise.
 *
 * The function removes a entries from the recipe frame category map.
 * After calling this function the recipe's frame map does not exist any
 * longer. Before category aliases are accessed again the recipe frame
 * category map must be (re-)initialized using the function 
 * pilTrnInitCategoryMap().
 */

int
pilTrnClearCategoryMap(void)
{

    if (catmap != NULL) {
        deletePilCatmap(catmap);
        catmap = NULL;
    }

    return 0;

}


/**
 * @brief
 *   Load category names mapping from file.
 *
 * @param filename Name of file containing categories aliases definition.
 *
 * @return @c EXIT_SUCCESS or @c EXIT_FAILURE.
 *
 * This function would be called if a configuration file defining
 * the category names mapping were found. Just one category mapping 
 * definition is allowed in one pipeline run, an error condition 
 * will occur otherwise.
 */

int pilTrnLoadCategoryMap(const char *filename)
{
  char  modName[] = "pilTrnLoadCategoryMap";
  char  line[PIL_LINE_LENGTH_MAX];
  char  alias[PIL_LINE_LENGTH_MAX];
  char  name[PIL_LINE_LENGTH_MAX];
  int   nameFlag, aliasFlag;
  FILE *fp;

  if ((fp = fopen (filename, "r")) == NULL) {
    pilMsgWarning(modName, "Problems opening category map file %s", filename);
    if (catmap == NULL) {
      pilMsgWarning(modName, "No default category map was loaded");
      return EXIT_FAILURE;
    } 
    else {
      pilMsgWarning(modName, "Using default category mapping only");
      return EXIT_SUCCESS;
    }
  }

  if (catmap == NULL) {
    pilMsgWarning(modName, "No default category names mapping loaded: "
                  "relying just on mapping from file %s", filename);
    catmap = newPilCatmap();
  }

  nameFlag = 0;
  aliasFlag = 0;
  
  while (fgets(line, PIL_LINE_LENGTH_MAX, fp)) {

   /* 
    * Empty line: assume category definition is complete
    */

    if (strempty(line, COMMENT_CHARS)) {
      pilMsgDebug(modName,"Empty line");
      
     /* 
      * Add category to table if all needed items are provided
      */

      if (nameFlag && aliasFlag) {
        if (pilTrnAddCategory(alias, name) == EXIT_FAILURE) {
          fclose(fp);
          return EXIT_FAILURE;
        }
        pilMsgDebug(modName, "Alias '%s' added to category map\n", alias);
      }
      else if (nameFlag || aliasFlag)
        pilMsgWarning(modName, "A category definition in category map "
                      "file %s is incomplete", filename);

      nameFlag = 0;
      aliasFlag = 0;
      continue;
    }

    if (sscanf(line, CATMAP_NAME "%[^\n]\n", name)) {
      strtrim(name,2);
      nameFlag = 1;
      pilMsgDebug(modName, "Name: %s\n", name);
      continue;
    }

    if (sscanf(line, CATMAP_ALIAS "%[^\n]\n", alias)) {
      strtrim(alias,2);
      aliasFlag = 1;
      pilMsgDebug(modName, "Alias: %s\n", alias);
      continue;
    }

  }
   
  fclose(fp);

 /* 
  * Add last category to table if all needed items are provided
  */

  if (nameFlag && aliasFlag) {
    if ((pilTrnAddCategory(alias, name) == EXIT_FAILURE)) {
       return EXIT_FAILURE;
    }
    pilMsgDebug(modName, "Alias '%s' added to category map\n", alias);
  }
  else if (nameFlag || aliasFlag)
    pilMsgWarning(modName, "A category definition in category map "
                  "file %s is incomplete", filename);

  return EXIT_SUCCESS;
}


/**
 * @brief
 *   Resolve a category alias.
 *
 * @param alias  Category alias
 *
 * @return The function returns a reference to the fully resolved category
 *   alias if no error occurred, otherwise the function returns @c NULL.
 *
 * The function searches the category map for the alias string
 * @em alias and returns a reference to the actual category name
 * corresponding to the given alias. The returned reference points into
 * the category map. No new memory is allocated for the category name
 * string.
 */

const char *pilTrnGetCategory(const char *alias)
{

  const char modName[] = "pilTrnGetCategory";
  const char *name = pilCatmapGetValue(catmap, alias);

  if (!name)
    pilMsgError(modName, "Translation of alias %s not found", alias);

  return name;

}
/**@}*/
