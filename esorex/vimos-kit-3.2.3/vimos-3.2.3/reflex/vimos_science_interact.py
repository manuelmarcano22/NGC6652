# import the needed modules
try:
  import reflex
  import matplotlib.gridspec
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
      self.sci_mapped    = None
      self.sci_extracted = None
      self.skylines_off  = None
      self.object_table  = None
      self.objid = 1

      #Read all the products
      frames = dict()
      for frame in fitsFiles:
        if frame == '' :
          continue
        category = frame.category
        frames[category] = frame

      if 'OBJECT_SCI_TABLE' in frames :
        self.object_table = frames['OBJECT_SCI_TABLE']

      if 'MOS_SCIENCE_EXTRACTED' in frames:
        self.sci_mapped   = PlotableMappedScience(frames["MOS_SCIENCE_EXTRACTED"],
                                                  self.object_table)
          
      if 'MOS_SCIENCE_REDUCED' in frames:
        self.sci_extracted   = PlotableExtractedScience(frames["MOS_SCIENCE_REDUCED"])
          
      if 'SKY_SHIFTS_SLIT_SCI_'in frames:
        self.skylines_off   = PlotableSkylinesOffsets(frames["SKY_SHIFTS_SLIT_SCI_"])
        if numpy.isnan(self.skylines_off.allwave_res).all() :
          self.skylines_off = None #All the offsets are NULL, not plotting
          
    # This function creates all the subplots. It is responsible for the plotting 
    # layouts. 
    # There can different layouts, depending on the availability of data
    # Note that subplot(I,J,K) means the Kth plot in a IxJ grid 
    # Note also that the last one is actually a box with text, no graphs.
    def addSubplots(self, figure):
      if self.sci_mapped is not None :
        if self.skylines_off is not None :
          gs = matplotlib.gridspec.GridSpec(7,7)
          self.subplot_sci_mapped    = figure.add_subplot(gs[0:2,:])
          self.subplot_sci_extracted = figure.add_subplot(gs[3:5,:])
          self.subplot_skylines_off  = figure.add_subplot(gs[6,:])
        else :
          gs = matplotlib.gridspec.GridSpec(4,6)
          self.subplot_sci_mapped    = figure.add_subplot(gs[0:2,0:6])
          self.subplot_sci_extracted = figure.add_subplot(gs[2:4,0:5])
          self.subplot_mplex_select   = figure.add_subplot(4,6,24)
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
      if self.sci_mapped is not None :

        #Reduced science
        self.plotSciMapped()

        #Extracted science
        if self.sci_extracted is not None :
            self.plotSciExtracted()
        else:
            self.showNoExtracted()

        #Sky lines offset
        self.plotSkyLinesOffsets()

      else :
        #Data not found info
        self.showNoData()
  
    def plotSciMapped(self) :
        title_sci_mapped   = 'Mapped science (not flux-calibrated)'
        tooltip_sci_mapped = """Mapped science (not flux-calibrated).
The regions of extracted spectra are marked by red 
(upper boundary of extraction region) 
and yellow (lower boundary) lines. 
By middle-clicking within one of these regions,
you select an extracted spectrum,
which will be plotted below."""
        self.sci_mapped.plot(self.subplot_sci_mapped, title_sci_mapped,
                             tooltip_sci_mapped)

    def plotSciExtracted(self) :
        title_sci_extracted = 'Extracted science spectrum (not flux-calibrated)'
        tooltip_sci_extracted ="""Extracted science spectrum (not flux-calibrated).
This is the sky-subtracted spectrum as extracted by the pipeline."""
        self.subplot_sci_extracted.clear()
        self.sci_extracted.plot(self.subplot_sci_extracted, title_sci_extracted,
                                tooltip_sci_extracted)

    def showNoExtracted(self) :
      #Data not found info
      self.subplot_sci_extracted.set_axis_off()
      self.text_nodata = 'No objects found. There are no extracted spectra.'
      self.subplot_sci_extracted.text(0.1, 0.6, self.text_nodata, color='#11557c', fontsize=18,
                               ha='left', va='center', alpha=1.0)
      self.subplot_sci_extracted.tooltip='No objects found. There are no extracted spectra'


    def plotSkyLinesOffsets(self) :
      if self.skylines_off is not None :
        title_skylines_off   = 'Sky lines offsets'
        tooltip_skylines_off ="""Sky lines offsets.
This plots shows the offsets of certain sky lines
(defined in the MASTER_SKYLINECAT table or 
pipeline values as given in the User Manual) 
relative to their reference position. 
The offsets are adjusted by default by 
a constant shift (--skyalign parameter)."""
        if self.sci_extracted is not None:
          wave_minlim = self.sci_extracted.spectrumdisplay.wave_lim[0]
          wave_maxlim = self.sci_extracted.spectrumdisplay.wave_lim[1]
          res_minlim = numpy.nanmin(self.skylines_off.allwave_res)
          res_maxlim = numpy.nanmax(self.skylines_off.allwave_res)
          if (res_minlim == res_maxlim) :
            res_minlim = 0
            res_maxlim = 2 * res_maxlim
          self.skylines_off.resdisplay.setLimits(wave_minlim, wave_maxlim,
                                                 res_minlim, res_maxlim)
        self.skylines_off.plot(self.subplot_skylines_off, title_skylines_off,
                               tooltip_skylines_off)

    def showNoData(self) :
      #Data not found info
      self.subtext_nodata.set_axis_off()
      self.text_nodata = 'Science data not found in the products (PRO.CATG=)'
      self.subtext_nodata.text(0.1, 0.6, self.text_nodata, color='#11557c', fontsize=18,
                               ha='left', va='center', alpha=1.0)
      self.subtext_nodata.tooltip='Science data not found in the products'

    def plotWidgets(self) :
      widgets = list()
      
      if self.sci_mapped is not None :
        #Object selector
        if self.sci_mapped is not None :
          # Clickable subplot
          self.clickablesmapped = InteractiveClickableSubplot(
            self.subplot_sci_mapped, self.setExtractedObject)
          widgets.append(self.clickablesmapped)
        
        # Multiplexing Selector radiobutton
        keys = []
        for i in range(1, len(self.sci_mapped.mappedsciences)+1) :
          keys.append(str(i))
        self.radiomplex =  InteractiveRadioButtons(self.subplot_mplex_select, 
            self.setMultiplexing, keys, 0, title='Mplex Selection')
        widgets.append(self.radiomplex)
      
      return widgets
    
    def setExtractedObject(self, point) :
      obj = self.sci_mapped.getObjectInPosition(point.ydata)
      if obj != -1 :
        self.sci_extracted.selectObject(obj)
      self.plotSciExtracted()

    def setMultiplexing(self, multiplexing_idx) :
      multiplexing_number = int(multiplexing_idx)
      self.sci_mapped.selectMultiplexing(multiplexing_number)
      self.sci_extracted.selectMultiplexing(multiplexing_number)
      self.subplot_sci_mapped.cla()
      self.subplot_sci_extracted.cla()
      self.plotSciMapped()
      self.plotSciExtracted()
        
                
    # This function specifies which are the parameters that should be presented
    # in the window to be edited.
    # Note that the parameter has to be also in the in_sop port (otherwise it 
    # won't appear in the window) 
    # The descriptions are used to show a tooltip. They should match one to one
    # with the parameter list 
    # Note also that parameters have to be prefixed by the 'recipe name:'
    def setInteractiveParameters(self):
      paramList = list()
      paramList.append(reflex.RecipeParameter('vmmosscience','skyglobal',group='Sky subtraction',description='Subtract global sky spectrum from CCD'))
      paramList.append(reflex.RecipeParameter('vmmosscience','skymedian',group='Sky subtraction',description='Sky subtraction from extracted slit spectra'))
      paramList.append(reflex.RecipeParameter('vmmosscience','skylocal',group='Sky subtraction',description='Sky subtraction from CCD slit spectra'))
      paramList.append(reflex.RecipeParameter('vmmosscience','cosmics',group='Sky subtraction',description='Eliminate cosmic rays hits (only if global sky subtraction is also requested)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','slit_margin',group='Spectra extraction',description='Number of pixels to exclude at each slit in object detection and extraction.'))
      paramList.append(reflex.RecipeParameter('vmmosscience','ext_radius',group='Spectra extraction',description='Maximum extraction radius for detected objects (pixel)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','cont_radius',group='Spectra extraction',description='Minimum distance at which two objects of equal luminosity do not contaminate each other (pixel)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','ext_mode',group='Spectra extraction',description='Object extraction method: 0 = aperture, 1 = Horne optimal extraction'))
      paramList.append(reflex.RecipeParameter('vmmosscience','skyalign',group='Wave calib',description='Polynomial order for sky lines alignment, or -1 to avoid alignment'))
      paramList.append(reflex.RecipeParameter('vmmosscience','detection',group='Spectra extraction',description='Object detection threshold (ADU)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','alignment',group='Dithering',description='Type of alignment of dithered frames, either to the nearest neighbour pixel or to fractions of pixel.  <integer | float>'))
      paramList.append(reflex.RecipeParameter('vmmosscience','dither',group='Dithering',description='Align dithered frames before stacking(for multiple input frames)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','compute',group='Dithering',description='Compute offsets of dithered images from detected objects (true), or read offsets from header (false)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','stack_method',group='Stacking',description='Frames combination method. <average | median | minmax | ksigma>'))
      paramList.append(reflex.RecipeParameter('vmmosscience','minrejection',group='Stacking',description='Number of lowest values to be rejected'))
      paramList.append(reflex.RecipeParameter('vmmosscience','maxrejection',group='Stacking',description='Number of highest values to be rejected'))
      paramList.append(reflex.RecipeParameter('vmmosscience','klow',group='Stacking',description='Low threshold in ksigma method'))
      paramList.append(reflex.RecipeParameter('vmmosscience','khigh',group='Stacking',description='High threshold in ksigma method'))
      paramList.append(reflex.RecipeParameter('vmmosscience','kiter',group='Stacking',description='Max number of iterations in ksigma method'))
      return paramList

    def setWindowHelp(self):
      help_text = """
In this window, the user will interact with the Vimos science reduction"""
      return help_text

    def setWindowTitle(self):
      title = 'Vimos Interactive Science Reduction'
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
