# import the needed modules
try:
  import matplotlib.gridspec as gridspec
  import reflex
  import_sucess = True

#NOTE for developers: 
# -If you want to modify the current script to cope
#  with different parameters, this is the function to modify:
#  setInteractiveParameters()
# -If you want to modify the current script to read different data from
#  the input FITS, this is the function to modify:
#  readFitsData()                  (from class DataPlotterManager) 
# -If you want to modify the current script to modify the plots (using the same
#  data),  this is the function to modify:
#  plotProductsGraphics()          (from class DataPlotterManager)
# -If you want to modify the text that appears in the "Help" button,
#  this is the function to modify:
#  setWindowHelp()
# -If you want to modify the title of the window, modify this function:
#  setWindowTitle()


  #This class deals with the specific details of data reading and final plotting.
  class DataPlotterManager:
    # This function will read all the columns, images and whatever is needed
    # from the products. The variables , self.plot_x, self.plot_y, etc...
    # are used later in function plotProductsGraphics().
    # Add/delete these variables as you need (only that plotProductsGraphics()
    # has to use the same names).
    # You can also create some additional variables (like statistics) after
    # reading the files.
    # If you use a control variable (self.xxx_found), you can modify 
    # later on the layout of the plotting window based on the presence of 
    # given input files. 
    # sof contains all the set of frames
    def readFitsData(self, fitsFiles):
      #Initialise the objects to read/display
      self.lamp_reduced    = None
      self.slit_map        = None
      self.flat_norm       = None
      self.flat_mscreen    = None
      self.disp_residuals  = None
      self.detected_lines  = None
      self.wave_map        = None
      self.excluded_lines  = list()

      #Read all the products
      frames = dict()
      for frame in fitsFiles:
        if frame == '' :
          continue
        category = frame.category
        frames[category] = frame


      if 'MOS_ARC_SPECTRUM_EXTRACTED' in frames:
        self.lamp_reduced = PlotableReducedArc(frames["MOS_ARC_SPECTRUM_EXTRACTED"])

      if 'MOS_SPATIAL_MAP_' in frames:
        self.slit_map = PlotableSpatialMap(frames["MOS_SPATIAL_MAP"])

      if 'MOS_DISP_RESIDUALS_TABLE' in frames:
        self.disp_residuals = PlotableDispResiduals(frames["MOS_DISP_RESIDUALS_TABLE"])

      if 'MOS_WAVELENGTH_MAP' in frames:
        self.wave_map = PlotableWavelengthMap(frames["MOS_WAVELENGTH_MAP"])

      if 'MOS_DETECTED_LINES' in frames:
        self.detected_lines = PlotableDetectedLines(frames["MOS_DETECTED_LINES"])

      curv_frame = None
      if 'MOS_CURV_COEFF' in frames :
          curv_frame = frames['MOS_CURV_COEFF']
                  
      flat_norm = None
      if 'MOS_MASTER_SCREEN_FLAT' in frames :
        flat_norm = frames['MOS_MASTER_SCREEN_FLAT']
        self.flat_norm  = PlotableNormFlat(flat_norm, curv_frame)
                  
      if 'MOS_SCREEN_FLAT' in frames :
        self.flat_mscreen  = PlotableRawFlat(frames['MOS_SCREEN_FLAT'],
#                                             flat_norm, curv_frame)
                                             None, None)
                  
      
    # This function creates all the subplots. It is responsible for the plotting 
    # layouts. 
    # There can different layouts, depending on the availability of data
    # Note that subplot(I,J,K) means the Kth plot in a IxJ grid 
    # Note also that the last one is actually a box with text, no graphs.
    def addSubplots(self, figure):
      if self.flat_norm is not None and self.lamp_reduced is not None \
         and self.disp_residuals is not None and self.flat_mscreen is not None \
         and self.detected_lines is not None and self.wave_map is not None:
        #self.subplot_slit_map       = figure.add_subplot(3,2,2)
        self.subplot_wave_map       = figure.add_subplot(2,3,1)
        self.subplot_flat_raw       = figure.add_subplot(2,3,2)
        self.subplot_flat_norm      = figure.add_subplot(2,3,3)
        gs = gridspec.GridSpec(80, 11)
        self.subplot_line_res_wave  = figure.add_subplot(gs[40:53,0:5])
        self.subplot_line_x_y       = figure.add_subplot(gs[40:53,6:11])
        self.subplot_lamp_reduced   = figure.add_subplot(gs[60:75,0:10])
        self.subplot_mplex_select   = figure.add_subplot(gs[60:75,10:11])
        self.subplot_txtinfo        = figure.add_subplot(gs[79:80,0:11])
      else : 
        self.subtext_nodata      = figure.add_subplot(1,1,1)
          
    # This is the function that makes the plots.
    # Add new plots or delete them using the given scheme.
    # The data has been already stored in self.plot_x, self.plot_xdif, etc ...
    # It is mandatory to add a tooltip variable to each subplot.
    # One might be tempted to merge addSubplots() and plotProductsGraphics().
    # There is a reason not to do it: addSubplots() is called only once at
    # startup, while plotProductsGraphics() is called always there is a resize.
    def plotProductsGraphics(self):
      if self.flat_norm is not None and self.lamp_reduced is not None \
         and self.disp_residuals is not None and self.flat_mscreen is not None \
         and self.wave_map is not None :

        #Reduced lamp
        self.plotReducedLamp()

        #Spatial map
        self.plotSpatialMap()

        #Master screen flat 
        self.plotScreenFlat()

        #Master flat normalise
        self.plotNormalisedFlat()
        
        #Dispersion residuals vs wave
        self.plotResiduals()

        #Wavelength maps
        self.plotWaveMap()

        #Line positions X vs Y
        self.plotLinePositions()
        
        #Additional text info
        self.showTextInfo()

      else :
        #Data not found info
        self.showNoData()
 
    def plotReducedLamp(self) :
      title_lamp_reduced   = 'Wavelength-calibrated arc lamp frame'
      tooltip_lamp_reduced ="""Wavelength-calibrated arc lamp frame.
Arc lines should be vertical without scatter or tilt. 
There should be no regions with no arc lines at all. 
Arc line coverage may vary with slit position."""
      self.lamp_reduced.plot(self.subplot_lamp_reduced, title_lamp_reduced,
                             tooltip_lamp_reduced)

    def plotSpatialMap(self) :
      if self.slit_map is not None:
        title_slit_map   = 'Slit spatial map'
        tooltip_slit_map ="""Reduced arc lamp frame."""
        self.slit_map.plot(self.subplot_slit_map, title_slit_map,
                           tooltip_slit_map)

    def plotScreenFlat(self) :
      title_flat_mscreen   = 'First raw flat field frame'
      tooltip_flat_mscreen ="""First raw flat field frame.
This frame is the first raw flat frame. It serves to verify that for MOS/MXU data all slits have been detected (compare with spatial map displayed above)."""
      self.flat_mscreen.plot(self.subplot_flat_raw, title_flat_mscreen,
                         tooltip_flat_mscreen)

    def plotNormalisedFlat(self) :
      title_flat_norm   = 'Normalised master flat frame'
      tooltip_flat_norm ="""Normalised master flat frame.
This is the result of the flat field normalisation. 
For LSS data it may make sense to keep the relative
flux variation along the slit to correct for slit illumination. 
For MOS/MXU data this generally does not work well."""
      self.flat_norm.plot(self.subplot_flat_norm, title_flat_norm,
                         tooltip_flat_norm)

    def plotResiduals(self) :
      title_line_res_wave   = 'Residuals of wavelength calibration'
      tooltip_line_res_wave ="""Residuals of wavelength calibration.
The residuals should not show any trends.  
Outliers may be tolerated, but at the edges 
they may indicate a questionable result. 
Clicking with middle button will add that line 
to the --ignore_lines recipe parameter."""
      self.subplot_line_res_wave.clear()
      self.disp_residuals.plotResVsWave(self.subplot_line_res_wave,
                                        title_line_res_wave,
                                        tooltip_line_res_wave,
                                        self.excluded_lines)
    def plotWaveMap(self) :
      title_wave_map   = 'Wavelength map'
      tooltip_wave_map ="""Wavelength map.
The value in each pixel represents the wavelength derived from the 
wavelength calibration for each slit."""
      self.wave_map.plot(self.subplot_wave_map, title_wave_map,
                         tooltip_wave_map)


    def plotLinePositions(self) :
      title_line_x_y   = 'Detected / Identified arc lines'
      tooltip_line_x_y ="""Detected/identified arc lines.
Black points are detected lines. 
Green points are identified lines.
Light green points are recovered identified lines after pattern-matching iteration
(only if wradius >0).
Red points are rejected lines in the wave vs pixel fit (only if wradius >0)
Most detected lines should also be identified. 
If lines at the edges are not identified this means that
the dispersion relation there will be extrapolated."""
      self.detected_lines.plotXVsY(self.subplot_line_x_y,
                                   title_line_x_y, tooltip_line_x_y)

    def showTextInfo(self) :
      self.subplot_txtinfo.set_axis_off()
      grism = 'unkown'
      for i in range(1,5):
        keyname = 'HIERARCH ESO INS GRIS%1d NAME' % (i)
        try:
          grism = self.lamp_reduced.arcs[0].readKeyword(keyname)
        except:
          pass
      filter_name = 'free'
      for i in range(1,5):
        keyname = 'HIERARCH ESO INS FILT%1d NAME' % (i)
        try:
          filter_name = self.lamp_reduced.arcs[0].readKeyword(keyname)
        except:
          pass

      self.subplot_txtinfo.text(0.0, -1.9, 'Grism/Filter: '+grism+"/"+
                                filter_name, 
                                ha='left', va='center', weight='bold')
      mask_id = 'free'
      for i in range(1,5):
        keyname = 'HIERARCH ESO INS MASK%1d ID' % (i)
        try:
          mask_id = self.lamp_reduced.arcs[0].readKeyword(keyname)
        except:
          pass
      self.subplot_txtinfo.text(0.6, -1.9, 'Mask id: '+str(mask_id), 
                               ha='left', va='center', weight='bold')

    def showNoData(self) :
      self.subtext_nodata.set_axis_off()
      self.text_nodata = 'Calibrations not found in the products'
      self.subtext_nodata.text(0.1, 0.6, self.text_nodata, color='#11557c', fontsize=18,
                               ha='left', va='center', alpha=1.0)
      self.subtext_nodata.tooltip='Calibrations  not found in the products'

  
    def plotWidgets(self) :
      widgets = list()
      if self.flat_norm is not None and self.lamp_reduced is not None \
         and self.disp_residuals is not None and self.flat_mscreen is not None:
        # Multiplexing Selector radiobutton
        keys = []
        for i in range(1, len(self.flat_norm.flats)+1) :
          keys.append(str(i))
        self.radiomplex =  InteractiveRadioButtons(self.subplot_mplex_select, 
            self.setMultiplexing, keys, 0, title='Mplex Select')
        widgets.append(self.radiomplex)
        # Clickable subplot
        self.clickableresidual = InteractiveClickableSubplot(
          self.subplot_line_res_wave, self.setExcludedLine)
        widgets.append(self.clickableresidual)

      return widgets

    def setMultiplexing(self, multiplexing_idx) :
      multiplexing_number = int(multiplexing_idx)
      self.flat_norm.selectMultiplexing(multiplexing_number)
      self.lamp_reduced.selectMultiplexing(multiplexing_number)
      self.disp_residuals.selectMultiplexing(multiplexing_number)
      self.wave_map.selectMultiplexing(multiplexing_number)
      self.detected_lines.selectMultiplexing(multiplexing_number)
      self.subplot_flat_norm.cla()
      self.subplot_lamp_reduced.cla()
      self.subplot_line_res_wave.cla()
      self.subplot_wave_map.cla()
      self.subplot_line_x_y.cla()
      self.plotNormalisedFlat()
      self.plotReducedLamp()
      self.plotResiduals()
      self.plotWaveMap()
      self.plotLinePositions()
      
    def setExcludedLine(self, point) :
      excluded_line = self.disp_residuals.getClosestLine(point.xdata)
      if excluded_line in self.excluded_lines :
        self.excluded_lines.remove(excluded_line)
      else :
        self.excluded_lines.append(excluded_line)
      self.plotResiduals()

      param_value = str()
      for line in self.excluded_lines :
        param_value = param_value + '%.4f, ' % line
      if len(param_value) > 1 :
        param_value = param_value[:-2]
      new_params = list()
      new_params.append(reflex.RecipeParameter('vmmoscalib','ignore_lines',
                                               value=param_value))
      return new_params
                
    # This function specifies which are the parameters that should be presented
    # in the window to be edited.
    # Note that the parameter has to be also in the in_sop port (otherwise it 
    # won't appear in the window) 
    # The descriptions are used to show a tooltip. They should match one to one
    # with the parameter list 
    # Note also that parameters have to be prefixed by the 'recipe name:'
    def setInteractiveParameters(self):
      paramList = list()
      paramList.append(reflex.RecipeParameter('vmmoscalib','sradius',group='flat field norm',description='Smooth box radius for flat field along spatial direction'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','dradius',group='flat field norm',description='Smooth box radius for flat field along dispersion direction (not used for long-slit-like data)'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','s_degree',group='flat field norm',description='Degree of flat field fitting polynomial along spatial direction (used for LSS-like data only)'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','d_nknots',group='flat field norm',description='Number of knots in flat field fitting splines along dispersion direction'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','slit_ident',group='master flat / slit traces',description='Attempt slit identification. For multiplexing data slit identification is always performed. This is must be on for a proper flux normalization of FLAT_SED.'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','cdegree',group='master flat / slit traces',description='Degree of spectral curvature polynomial'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','reference',group='master flat / slit traces',description='Reference wavelength for slit map'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','wreject',group='wave calib / distorsions',description='Rejection threshold in dispersion relation fit (pixel)'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','peakdetection',group='wave calib / distorsions',description='Initial peak detection threshold (ADU)'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','used_linesets',group='wave calib / distorsions',description='Linesets to use. Valid are standard and extended (see column LINE_SET in the line catalogue)'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','ignore_lines',group='wave calib / distorsions',description='Catalog lines nearest to wavelengths in this list will be ignored for wavelength calibration'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','wdegree',group='wave calib / distorsions',description='Degree of wavelength calibration polynomial'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','dispersion',group='wave calib / distorsions',description='Expected spectral dispersion (Angstrom/pixel)'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','wradius',group='wave calib / distorsions',description='Search radius if iterating pattern-matching with first-guess method'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','startwavelength',group='wave calib / distorsions',description='Start wavelength in spectral extraction'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','endwavelength',group='wave calib / distorsions',description='End wavelength in spectral extraction'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','wmodelss',group='wave calib / distorsions',description='Interpolation mode of wavelength solution (0 = no interpolation, 1 = fill gaps, 2 = global model)'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','wmodemos',group='wave calib / distorsions',description='Interpolation mode of wavelength solution (0 = no interpolation, 1 = local (slit) solution, 2 = global model)'))
      paramList.append(reflex.RecipeParameter('vmmoscalib','line_ident_tol',group='wave calib / distorsions',description='Tolerance for the ratio of detected lines vs reference lines. This is used during for arc line identification'))
      return paramList

    def setWindowHelp(self):
      help_text = """
In this window, the user will interact with the Vimos calibration (flat, slit traces and wavelength calibration)"""
      return help_text

    def setWindowTitle(self):
      title = 'Vimos Interactive Calibration'
      return title

except ImportError:
  import_sucess = 'false'
  print "Error importing modules pyfits, wx, matplotlib, numpy"

#This is the 'main' function
if __name__ == '__main__':

  # import reflex modules
  import reflex_interactive_app
  import sys

  # import UVES reflex modules
  from vimos_plot_common import *

  # Create interactive application
  interactive_app = reflex_interactive_app.PipelineInteractiveApp(enable_init_sop=True)

  #Check if import failed or not
  if import_sucess == 'false' :
    interactive_app.setEnableGUI('false')

  #Open the interactive window if enabled
  if interactive_app.isGUIEnabled() :
    #Get the specific functions for this window
    dataPlotManager = DataPlotterManager()
    interactive_app.setPlotManager(dataPlotManager)
    interactive_app.showGUI()
  else :
    interactive_app.passProductsThrough()

  # print outputs
  interactive_app.print_outputs()

  sys.exit()
