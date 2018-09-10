try:
  import numpy
  from pipeline_display import *
  from pipeline_product import *
  from reflex_plot_widgets import *
  from numpy.polynomial import Polynomial
  numpy.seterr(invalid='ignore')
except ImportError:
  donothing=1


class PlotableReducedArc :
  def __init__(self, fits):
    arc = PipelineProduct(fits)
    self.multiext = False
    self.ext_sel = 0
    if len(arc.all_hdu) > 1:
      self.multiext = True
    self.arcs = []
    if(self.multiext) :
      for iext in range(1, len(arc.all_hdu)) :
        self.arcs.append(PipelineProduct(fits))
    else :
        self.arcs.append(arc)
    self.arcdisp = ImageDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    if(self.multiext) :
      read_ext = 1
    else:
      read_ext = 0
    for arc in self.arcs :
      arc.readImage(read_ext)
      arc.read2DLinearWCS()
      read_ext = read_ext + 1

  def plot(self, subplot, title, tooltip):
    self.arcdisp.setLabels('Lambda [Angstrom]', 'Y [pix]')
    image_clean = self.arcs[self.ext_sel].image[numpy.isfinite(self.arcs[self.ext_sel].image)]
    avg = numpy.median(image_clean)
    quartiles = numpy.percentile(image_clean, [25, 75])
    self.z_lim = avg -  30 * (avg - quartiles[0]), avg + 30 * (quartiles[1] - avg)
    self.arcdisp.setZLimits(self.z_lim)
    self.arcdisp.setXLinearWCSAxis(self.arcs[self.ext_sel].crval1,
                                   self.arcs[self.ext_sel].cdelt1, 
                                   self.arcs[self.ext_sel].crpix1)
    self.arcdisp.display(subplot, title, tooltip, self.arcs[self.ext_sel].image)
    
  def selectMultiplexing(self, multiplexing_idx):
    self.ext_sel = multiplexing_idx - 1
    
class PlotableWavelengthMap :
  def __init__(self, fits):
    map = PipelineProduct(fits)
    self.multiext = False
    self.ext_sel = 0
    if len(map.all_hdu) > 1:
      self.multiext = True
    self.maps = []
    if(self.multiext) :
      for iext in range(1, len(map.all_hdu)) :
        self.maps.append(PipelineProduct(fits))
    else :
        self.maps.append(map)
    self.mapdisp = ImageDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    if(self.multiext) :
      read_ext = 1
    else:
      read_ext = 0
    for map in self.maps :
      map.readImage(read_ext)
      read_ext = read_ext + 1

  def plot(self, subplot, title, tooltip):
    self.mapdisp.setZLimits((3000, 10000.))
    self.mapdisp.display(subplot, title, tooltip, self.maps[self.ext_sel].image)
    
  def selectMultiplexing(self, multiplexing_idx):
    self.ext_sel = multiplexing_idx - 1
    
class PlotableFlat(object) :
  def __init__(self, fits_flat, fits_slittrace):
    flat      = PipelineProduct(fits_flat)
    self.multiext = False
    self.ext_sel = 0
    if len(flat.all_hdu) > 2:
      self.multiext = True
    self.flats = []
    if(self.multiext) :
      for iext in range(1, len(flat.all_hdu), 2) :
        self.flats.append(PipelineProduct(fits_flat))
    else :
        self.flats.append(flat)
    self.slittraces = []
    if fits_slittrace is not None:
      slittrace = PipelineProduct(fits_slittrace)
      if(self.multiext) :
        for iext in range(1, len(slittrace.all_hdu)) :
          self.slittraces.append(PipelineProduct(fits_slittrace))
      else :
          self.slittraces.append(slittrace)
    self.flatdisp  = ImageDisplay()
    self.loadFromFits()
    

  def loadFromFits(self) :
    #Reading the flat image
    if(self.multiext) :
      read_ext = 1
    else:
      read_ext = 0
    for flat in self.flats :
      flat.readImage(read_ext)
      read_ext = read_ext + 2
    
    #Reading the polinomial traces
    if len(self.slittraces) >= 1:
      self.nslits = []
      self.ypos_traces = []
      self.xpos_top_traces = []
      self.xpos_bottom_traces = []
      for i in range(0, len(self.flats)) :
        slittrace = self.slittraces[i]
        ndegree = slittrace.getTableNcols(i+1) - 1
        nslits = slittrace.getTableNrows(i+1) / 2
        self.nslits.append(nslits)
        degreecols = []
        for deg in range(ndegree):
          colname = 'c%d'%deg
          slittrace.readTableColumn(i+1, colname)
          degreecols.append(slittrace.column)
    
        top_trace_polynomials = []
        bottom_trace_polynomials = []
        for slit in range(nslits) :
          top_trace_coeff = []
          bottom_trace_coeff = []
          for deg in range(ndegree) :
            top_trace_coeff.append(degreecols[deg][2*slit])
            bottom_trace_coeff.append(degreecols[deg][2*slit + 1])
        
          top_trace_pol = Polynomial(top_trace_coeff)  
          bottom_trace_pol = Polynomial(bottom_trace_coeff)
          top_trace_polynomials.append(top_trace_pol) 
          bottom_trace_polynomials.append(bottom_trace_pol) 

        #Creating the points to plot based on the polynomial traces
        #We interchange X and Y since the traces are rotated
        ypos_traces = []
        xpos_top_traces = []
        xpos_bottom_traces = []
        for slit in range(self.nslits[i]) :
          xpos_top = []
          xpos_bottom = []
          ypos = []
          for ypix in range(self.flats[i].image.shape[0]) :
            ypos.append(ypix+1) 
            xpos_top.append(self.flats[i].image.shape[1] - top_trace_polynomials[slit](ypix)+1) 
            xpos_bottom.append(self.flats[i].image.shape[1] - bottom_trace_polynomials[slit](ypix)+1)
          ypos_traces.append(ypos)
          xpos_top_traces.append(xpos_top) 
          xpos_bottom_traces.append(xpos_bottom)
        self.ypos_traces.append(ypos_traces)
        self.xpos_top_traces.append(xpos_top_traces)
        self.xpos_bottom_traces.append(xpos_bottom_traces)

  def plot(self, subplot, title, tooltip):
    self.flatdisp.setLabels('X [pix]', 'Y [pix]')
    self.flatdisp.display(subplot, title, tooltip, self.flats[self.ext_sel].image)
    
    if len(self.slittraces) >= 1:
      subplot.autoscale(enable=False)
      for slit in range(self.nslits[self.ext_sel]) :
        subplot.plot(self.xpos_top_traces[self.ext_sel][slit], self.ypos_traces[self.ext_sel][slit],
                     linestyle='solid',color='red')
        subplot.plot(self.xpos_bottom_traces[self.ext_sel][slit], self.ypos_traces[self.ext_sel][slit],
                     linestyle='solid',color='darkred')

  def selectMultiplexing(self, multiplexing_idx):
    self.ext_sel = multiplexing_idx - 1
    
class PlotableNormFlat (PlotableFlat) :
  def __init__(self, fits_flat, fits_slittrace):
    super(PlotableNormFlat, self).__init__(fits_flat, fits_slittrace)

  def plot(self, subplot, title, tooltip):
    self.flatdisp.setZLimits((0.9, 1.1))
    self.flats[self.ext_sel].image[self.flats[self.ext_sel].image > 5.] = 0
    super(PlotableNormFlat, self).plot(subplot, title, tooltip)

class PlotableRawFlat (PlotableFlat) :
  def __init__(self, fits_flat_raw, fits_master_flat, fits_slittrace):
    super(PlotableRawFlat, self).__init__(fits_flat_raw, fits_slittrace)
    if len(self.slittraces) >= 1 and fits_master_flat is not None :
      master_flat = PipelineProduct(fits_master_flat)
      self.trimm_lly = master_flat.all_hdu[0].header.get('HIERARCH ESO QC TRIMM LLY')

      #Change the traces by the amount of overscan in Y that has been removed
      #TODO: Fix this
      for xpos_top in self.xpos_top_traces[self.ext_sel]:
        for j, xpos in enumerate(xpos_top):
          xpos_top[j] = xpos + self.trimm_lly - 1
      for xpos_bottom in self.xpos_bottom_traces[self.ext_sel]:
        for j, xpos in enumerate(xpos_bottom):
          xpos_bottom[j] = xpos + self.trimm_lly -1
    

  def plot(self, subplot, title, tooltip):
    super(PlotableRawFlat, self).plot(subplot, title, tooltip)

class PlotableSpatialMap :
  def __init__(self, fits_spatialmap):
    self.spatialmap      = PipelineProduct(fits_spatialmap)
    self.spatialmapdisp  = ImageDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the flat image
    self.spatialmap.readImage()

  def plot(self, subplot, title, tooltip):
    self.spatialmapdisp.setLabels('X', 'Y')
    self.spatialmapdisp.setZLimits((0., 100))
    self.spatialmapdisp.display(subplot, title, tooltip, self.spatialmap.image)

class PlotableMappedScience :
  def __init__(self, fits_mappedscience, fits_objecttable):
    self.ext_sel = 0
    mappedscience      = PipelineProduct(fits_mappedscience)
    self.multiext = False
    self.ext_sel = 0
    if len(mappedscience.all_hdu) > 2:
      self.multiext = True
    self.mappedsciences = []
    if(self.multiext) :
      for iext in range(1, len(mappedscience.all_hdu), 2) :
        self.mappedsciences.append(PipelineProduct(fits_mappedscience))
    else :
        self.mappedsciences.append(mappedscience)

    self.mappedsciencedisp  = ImageDisplay()
    if fits_objecttable is not None:
      self.objecttables = []
      for iext in range(0, len(self.mappedsciences)) :
        self.objecttables.append(PipelineProduct(fits_objecttable))
    else :
      self.objecttables  = None
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the object images
    if(self.multiext) :
      read_ext = 1
    else:
      read_ext = 0
    for mappedscience in self.mappedsciences :
      mappedscience.readImage(read_ext)
      read_ext = read_ext + 2
    
    #Reading the object table
    if self.objecttables is not None:
      self.nobjects = []
      self.ybottom_obj_extracts = []
      self.ytop_obj_extracts = []
      for iext in range(0, len(self.objecttables)) :
        objecttable = self.objecttables[iext]
        nslit = objecttable.getTableNrows(iext+1)
        n_skip_columns = 10
        for col in objecttable.all_hdu[iext+1].data.columns:
          if(col.name == 'multiplex'):
            n_skip_columns = 12
        maxobjectperslit = (objecttable.getTableNcols(iext+1) - n_skip_columns ) / 4
        start_extracted_cols = []
        end_extracted_cols = []
        for obj in range(maxobjectperslit):
          colname = 'start_%d'%(obj+1)
          objecttable.readTableColumn(iext+1, colname)
          start_extracted_cols.append(objecttable.column)
          colname = 'end_%d'%(obj+1)
          objecttable.readTableColumn(iext+1, colname)
          end_extracted_cols.append(objecttable.column)

        ybottom_obj_extract = []
        ytop_obj_extract = []
        for slit in range(nslit) :
          for obj in range(maxobjectperslit) :
            ybottom = start_extracted_cols[obj][slit]
            ytop    = end_extracted_cols[obj][slit]
            if ybottom != -1 :
              ybottom_obj_extract.append(ybottom)
              ytop_obj_extract.append(ytop)
        self.ybottom_obj_extracts.append(ybottom_obj_extract)
        self.ytop_obj_extracts.append(ytop_obj_extract)
        self.nobjects.append(len(ybottom_obj_extract))

  def plot(self, subplot, title, tooltip):
    self.mappedsciencedisp.setLabels('X [pix]', 'Y [pix]')
    self.mappedsciencedisp.display(subplot, title, tooltip, 
                                   self.mappedsciences[self.ext_sel].image)
    if self.objecttables is not None:
      subplot.autoscale(enable=False)
      for obj in range(self.nobjects[self.ext_sel]) :
        subplot.axhline(self.ytop_obj_extracts[self.ext_sel][obj], linestyle='solid',color='red')
        subplot.axhline(self.ybottom_obj_extracts[self.ext_sel][obj], linestyle='solid',color='yellow')
        
  def getObjectInPosition(self, ypos) :
    for obj in range(self.nobjects[self.ext_sel]) :
      if ypos > self.ybottom_obj_extracts[self.ext_sel][obj] and \
         ypos < self.ytop_obj_extracts[self.ext_sel][obj] :
        return self.nobjects[self.ext_sel] - obj
    return -1

  def selectMultiplexing(self, multiplexing_idx):
    self.ext_sel = multiplexing_idx - 1

class PlotableDispResiduals :
  def __init__(self, fits_dispresiduals):
    self.multiext = False
    self.ext_sel = 0
    disp = PipelineProduct(fits_dispresiduals)
    self.next = len(disp.all_hdu) - 1
    if len(disp.all_hdu) > 2:
      self.multiext = True
    self.dispresiduals = []
    for iext in range(1, len(disp.all_hdu)) :
      self.dispresiduals.append(PipelineProduct(fits_dispresiduals))
    self.resdisplay  = ScatterDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    self.residuals = []
    self.allwave = []
    self.allypos = []
    self.allresiduals = []
    for iext in range(0, self.next) :
      #Reading the residuals table
      self.dispresiduals[iext].readTableColumn(iext + 1, 'wavelength')
      wave     = self.dispresiduals[iext].column
      nwave    = self.dispresiduals[iext].getTableNrows(iext+1)
      ncolumns = self.dispresiduals[iext].getTableNcols(iext+1)
      nselectedrows = (ncolumns - 1) // 3
      this_residuals = []
      this_allwave = []
      this_allypos = []
      this_allresiduals = []
      for i in range(nselectedrows) :
        #TODO: Currently the residuals are computed every 10 rows. 
        #This is hard-coded in the pipeline. It would be better just to detect the
        #columns whose name start with 'r' 
        colname = 'r%d'%(i*10) 
        self.dispresiduals[iext].readTableColumn(iext + 1, colname)
        row_residuals = self.dispresiduals[iext].column
        this_residuals.append(row_residuals)
        this_allwave.extend(wave)
        this_allresiduals.extend(row_residuals)
        ypos = i*10.
        this_allypos.extend([ypos] * nwave)
      self.residuals.append(this_residuals)
      self.allwave.append(this_allwave)
      self.allypos.append(this_allypos)
      self.allresiduals.append(this_allresiduals)

  def plotResVsWave(self, subplot, title, tooltip, excluded_lines = None):
    self.resdisplay.setLabels('Wavelength [Ang]','Residual [pix]')
    self.resdisplay.display(subplot, title, tooltip, self.allwave[self.ext_sel],
                            self.allresiduals[self.ext_sel])
    if excluded_lines is not None :
      for line in excluded_lines :
        subplot.axvline(line, linestyle='solid',color='red')


  def plotResVsY(self, subplot, title, tooltip):
    self.resdisplay.setLabels('Ypos [pix]','Residual [pix]')
    self.resdisplay.display(subplot, title, tooltip, self.allypos[self.ext_sel],
                            self.allresiduals[self.ext_sel])

  def selectMultiplexing(self, multiplexing_idx):
    self.ext_sel = multiplexing_idx - 1

  def getClosestLine(self, wave_selected) :
    distance = numpy.fabs(self.allwave[self.ext_sel] - wave_selected)
    idx = numpy.nanargmin(distance)
    return self.allwave[self.ext_sel][idx]

class PlotableDetectedLines :
  def __init__(self, fits_detectedlines):
    self.multiext = False
    self.ext_sel = 0
    detect = PipelineProduct(fits_detectedlines)
    self.next = len(detect.all_hdu) - 1
    if len(detect.all_hdu) > 2:
      self.multiext = True
    self.detectedlines = []
    for iext in range(1, len(detect.all_hdu)) :
      self.detectedlines.append(PipelineProduct(fits_detectedlines))

    self.xydisplay     = ScatterDisplay()
    self.resdisplay    = ScatterDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the residuals table
    self.x_pix = []
    self.y_pix = []
    self.x_pix_iter = []
    self.y_pix_iter = []
    self.wave  = []
    self.wave_iter  = []
    self.res_xpos  = []
    self.fit_used = []

    for iext in range(0, self.next) :
    
      try :
        self.detectedlines[iext].readTableColumn(iext + 1, 'xpos_rectified')
        x_pix = self.detectedlines[iext].column
        self.detectedlines[iext].readTableColumn(iext + 1, 'ypos_rectified')
        y_pix = self.detectedlines[iext].column
        self.detectedlines[iext].readTableColumn(iext + 1, 'xpos_rectified_iter')
        x_pix_iter = self.detectedlines[iext].column
        self.detectedlines[iext].readTableColumn(iext + 1, 'ypos_rectified_iter')
        y_pix_iter = self.detectedlines[iext].column
      except KeyError:
        self.detectedlines[iext].readTableColumn(iext + 1, 'xpos')
        x_pix = self.detectedlines[iext].column
        self.detectedlines[iext].readTableColumn(iext + 1, 'ypos')
        y_pix = self.detectedlines[iext].column
        self.detectedlines[iext].readTableColumn(iext + 1, 'xpos_iter')
        x_pix_iter = self.detectedlines[iext].column
        self.detectedlines[iext].readTableColumn(iext + 1, 'ypos_iter')
        y_pix_iter = self.detectedlines[iext].column

      
      self.detectedlines[iext].readTableColumn(iext + 1, 'wave_ident')
      wave  = self.detectedlines[iext].column
      self.detectedlines[iext].readTableColumn(iext + 1, 'wave_ident_iter')
      wave_iter  = self.detectedlines[iext].column
      self.detectedlines[iext].readTableColumn(iext + 1, 'res_xpos')
      res_xpos  = self.detectedlines[iext].column
      self.detectedlines[iext].readTableColumn(iext + 1, 'fit_used')
      fit_used = self.detectedlines[iext].column

      self.x_pix.append(x_pix)
      self.y_pix.append(y_pix)
      self.x_pix_iter.append(x_pix_iter)
      self.y_pix_iter.append(y_pix_iter)
      self.wave.append(wave)
      self.wave_iter.append(wave_iter)
      self.res_xpos.append(res_xpos)
      self.fit_used.append(fit_used)

  def plotXVsY(self, subplot, title, tooltip):
    #We first plot all the detected lines
    self.xydisplay.setLabels('Xpos [pix]','Ypos [pix]')
    self.xydisplay.setColor('black')
    self.xydisplay.display(subplot, title, tooltip, self.x_pix[self.ext_sel],
                           self.y_pix[self.ext_sel])
    #We then overplot the identified lines in the second iteration
    self.xydisplay.setColor('lightgreen')
    self.xydisplay.display(subplot, title, tooltip,
                           self.x_pix_iter[self.ext_sel][numpy.isfinite(self.wave_iter[self.ext_sel])],
                           self.y_pix_iter[self.ext_sel][numpy.isfinite(self.wave_iter[self.ext_sel])])
    #And then we overplot the identified lines in the first iteration
    self.xydisplay.setColor('green')
    self.xydisplay.display(subplot, title, tooltip, 
                           self.x_pix[self.ext_sel][numpy.isfinite(self.wave[self.ext_sel])],
                           self.y_pix[self.ext_sel][numpy.isfinite(self.wave[self.ext_sel])])
    #We then overplot the identified lines which have been rejected in the fit
    self.xydisplay.setColor('red')
    identified = numpy.logical_or(numpy.isfinite(self.wave_iter[self.ext_sel]), numpy.isfinite(self.wave[self.ext_sel]))
    not_used_fit = numpy.logical_and(identified, self.fit_used[self.ext_sel] == 0)
    self.xydisplay.display(subplot, title, tooltip,
                           self.x_pix_iter[self.ext_sel][not_used_fit],
                           self.y_pix_iter[self.ext_sel][not_used_fit])

  def plotResVsWave(self, subplot, title, tooltip, excluded_lines = None):
    self.resdisplay.setLabels('Wavelength [Ang]','Residual [pix]')
    self.resdisplay.setColor('black')
    self.resdisplay.display(subplot, title, tooltip, 
                            self.wave[numpy.isfinite(self.res_xpos)],
                            self.res_xpos[numpy.isfinite(self.res_xpos)])
    if excluded_lines is not None :
      for line in excluded_lines : 
        subplot.axvline(line, linestyle='solid',color='red')

  def selectMultiplexing(self, multiplexing_idx):
    self.ext_sel = multiplexing_idx - 1


class PlotableSkylinesOffsets :
  def __init__(self, fits_skylines_off):
    self.skylines_off  = PipelineProduct(fits_skylines_off)
    self.resdisplay    = ScatterDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the slylines offset table
    nslits = self.skylines_off.getTableNcols(1) - 1

    skylines_wave = self.skylines_off.readTableColumn(1, 'wave')
    self.allskylines_wave = list()
    self.allwave_res = list()
    
    for col in range(nslits) :
      self.allskylines_wave.extend(skylines_wave)
      wave_res = self.skylines_off.readTableColumn(1, col + 1)
      self.allwave_res.extend(wave_res)
      
  def plot(self, subplot, title, tooltip):
    self.resdisplay.setLabels('Wavelength [Ang]','Residual [Ang]')
    self.resdisplay.setColor('black')
    self.resdisplay.setPointSize(7)
    self.resdisplay.display(subplot, title, tooltip, 
                            self.allskylines_wave, self.allwave_res)

class PlotableExtractedScience :
  def __init__(self, fits_extractedscience):
    self.obj_id = -1
    self.ext_sel = 0
    extractedscience      = PipelineProduct(fits_extractedscience)
    self.multiext = False
    self.ext_sel = 0
    if len(extractedscience.all_hdu) > 1:
      self.multiext = True
    self.extractedsciences = []
    if(self.multiext) :
      for iext in range(1, len(extractedscience.all_hdu)) :
        self.extractedsciences.append(PipelineProduct(fits_extractedscience))
    else :
        self.extractedsciences.append(extractedscience)
    self.spectrumdisplay   = SpectrumDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #Reading the object images
    if(self.multiext) :
      read_ext = 1
    else:
      read_ext = 0
      

    self.nobj = []
    self.bunit = []
    self.wave = []
    for i in range(0, len(self.extractedsciences)) :
      extractedscience = self.extractedsciences[i]
      extractedscience.readImage(read_ext)
      read_ext = read_ext + 1

      self.nobj.append(extractedscience.image.shape[0])
      crpix1  = extractedscience.readKeyword('CRPIX1', 0)
      crval1  = extractedscience.readKeyword('CRVAL1', 0)
      cdelt1  = extractedscience.readKeyword('CD1_1', 0)
      self.bunit.append(extractedscience.readKeyword('BUNIT', 0))
      nwave  =   extractedscience.image.shape[1]
      self.wave.append(numpy.arange(1, nwave+1, 1))
      self.wave[i] = (self.wave[i] - crpix1) * cdelt1 + crval1
      
    if(self.obj_id == -1) : # Select brightest
      self.selectBrightest()
    self.setFluxSelected()

  def selectBrightest(self):
    if self.nobj[self.ext_sel] == 1:
      self.obj_id = 1
    median = 0
    for obj in range(self.nobj[self.ext_sel]) :
      new_median = numpy.median(self.extractedsciences[self.ext_sel].image[obj,:]) 
      if new_median > median :
        median = new_median
        self.obj_id = obj + 1
  
  def setFluxSelected(self) :
    self.flux = self.extractedsciences[self.ext_sel].image[self.obj_id-1,:]

  def selectObject(self, obj_id):
    self.obj_id = obj_id
    self.setFluxSelected()

  def selectMultiplexing(self, multiplexing_idx):
    self.ext_sel = multiplexing_idx - 1
    self.selectBrightest()
    self.setFluxSelected()

  def plot(self, subplot, title, tooltip):

    self.spectrumdisplay.setLabels('Lambda', 'Total Flux ['+self.bunit[self.ext_sel]+']')
    self.spectrumdisplay.display(subplot, title, tooltip, self.wave[self.ext_sel], self.flux,
                                autolimits = True)

class PlotableSpecPhot :
  def __init__(self, fits):
    self.resp     = PipelineProduct(fits)
    self.respdisp = SpectrumDisplay()
    self.tabdisp  = ScatterDisplay()
    self.flat_sed = False
    self.loadFromFits()

  def loadFromFits(self) :
    self.wave         = self.resp.readTableColumn(1, 'WAVE')
    self.wave_obs     = self.resp.readTableColumn(2, 'WAVE')
    self.std_ref_flux = self.resp.readTableColumn(1, 'STD_FLUX')
    self.std_obs_flux = self.resp.readTableColumn(1, 'OBS_FLUX')
    if 'RESPONSE' in self.resp.all_hdu[2].columns.names :
      self.fit_response = self.resp.readTableColumn(2, 'RESPONSE')
      self.raw_response = self.resp.readTableColumn(1, 'RAW_RESPONSE')
    else :
      self.fit_response = self.resp.readTableColumn(2, 'RESPONSE_FFSED')
      self.raw_response = self.resp.readTableColumn(1, 'RAW_RESPONSE_FFSED')
      self.flat_sed = True
    self.used_fit     = self.resp.readTableColumn(1, 'USED_FIT')
    self.raw_response_nonnull = self.raw_response[self.raw_response > 0]
    self.wave_nonnull = self.wave[self.raw_response > 0]
    self.wave_used = self.wave[self.used_fit > 0]
    self.raw_response_used = self.raw_response[self.used_fit > 0]    

  def plotResponse(self, subplot, title, tooltip):
    self.respdisp.setLabels('$\lambda\, [\AA]$','$10^{-16} erg\, cm^{-2} e-^{-1}$')
    self.respdisp.flux_lim = 0., numpy.max(self.raw_response_nonnull) * 1.1
    self.respdisp.display(subplot, title, tooltip, self.wave_obs, self.fit_response, autolimits = False)
    subplot.scatter(self.wave_nonnull, self.raw_response_nonnull, color='darkblue')
    subplot.scatter(self.wave_used, self.raw_response_used, color='lightgreen')

  def plotStdExtracted(self, subplot, title, tooltip):
    self.respdisp.setLabels('$\lambda\, [\AA]$','$e-\, s^{-1} \AA^{-1}$')
    std_obs_flux_nonnull = self.std_obs_flux[self.std_obs_flux > 0]
    wave_nonnull = self.wave[self.std_obs_flux > 0]
    self.respdisp.display(subplot, title, tooltip,
                          wave_nonnull, std_obs_flux_nonnull, autolimits = True)

  def plotStdTabulated(self, subplot, title, tooltip):
    self.tabdisp.setLabels('$\lambda\, [\AA]$','$10^{-16} erg\,  cm^{-2} s^{-1} \AA-^{-1}$')
    self.tabdisp.display(subplot, title, tooltip, self.wave, self.std_ref_flux)

class PlotableStdTabRedFlux :
  def __init__(self, reducedfluxstd_fits, reducedstd_fits, specphot_fits):
    self.reducedfluxstd = PlotableExtractedScience(reducedfluxstd_fits)
    self.reducedstd     = PlotableExtractedScience(reducedstd_fits)
    self.specphot       = PlotableSpecPhot(specphot_fits)
    self.tabstddisp     = ScatterDisplay()
    self.stdreddisp     = ScatterDisplay()
    self.loadFromFits()

  def loadFromFits(self) :
    #This will select the brightest spectrum, which is the criteria
    #used to extract the standard star
    self.reducedfluxstd.loadFromFits()
    self.reducedstd.loadFromFits()
    self.specphot.loadFromFits()
    self.std_ref_flux_nonnull = self.specphot.std_ref_flux[self.specphot.raw_response > 0]
    self.std_ref_flux_used = self.specphot.std_ref_flux[self.specphot.used_fit > 0]
    

  def plotStdTabRedFlux(self, subplot, title, tooltip) :
    self.tabstddisp.setLabels('$\lambda\, [\AA]$','$10^{-16} erg\, cm^{-2} s^{-1} \AA^{-1}$')
    self.tabstddisp.setLimits(self.reducedfluxstd.wave[self.reducedfluxstd.ext_sel][0],
                              self.reducedfluxstd.wave[self.reducedfluxstd.ext_sel][len(self.reducedfluxstd.wave[self.reducedfluxstd.ext_sel])-1],
                              0.,
                              numpy.max(self.reducedfluxstd.flux) * 1.1)
    self.tabstddisp.setColor('red') 
    self.tabstddisp.display(subplot, title, tooltip, 
                            self.reducedfluxstd.wave[self.reducedfluxstd.ext_sel], self.reducedfluxstd.flux)
    self.tabstddisp.setColor('darkblue') 
    self.tabstddisp.setPointSize(20) 
    subplot.scatter(self.specphot.wave, self.specphot.std_ref_flux)
    subplot.scatter(self.specphot.wave_used, self.std_ref_flux_used, color='lightgreen')


  def plotStdRed(self, subplot, title, tooltip) :
    self.stdreddisp.setLabels('$\lambda\, [\AA]$','$ADU/s$')
    self.stdreddisp.setLimits(self.reducedstd.wave[self.reducedstd.ext_sel][0],
                              self.reducedstd.wave[self.reducedstd.ext_sel][len(self.reducedstd.wave[self.reducedstd.ext_sel])-1],
                              0.,
                              numpy.max(self.reducedstd.flux) * 1.1) 
    self.stdreddisp.setColor('red') 
    self.stdreddisp.display(subplot, title, tooltip, 
                            self.reducedstd.wave[self.reducedstd.ext_sel], self.reducedstd.flux)



