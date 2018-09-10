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
      self.specphot   = None
      self.stdflux    = None

      #Read all the products
      frames = dict()
      for frame in fitsFiles:
        if frame == '' :
          continue
        category = frame.category
        frames[category] = frame

      if 'MOS_SPECPHOT_TABLE' in frames:
        self.specphot   = PlotableSpecPhot(frames["MOS_SPECPHOT_TABLE"])

      reducedfluxstd = None
      reducedstd = None
      if 'MOS_STANDARD_FLUX_REDUCED' in frames :
        reducedfluxstd = frames['MOS_STANDARD_FLUX_REDUCED']
      
      if 'MOS_STANDARD_REDUCED' in frames :
        reducedstd = frames['MOS_STANDARD_REDUCED']
      
      if 'STD_FLUX_TABLE' in frames :
        if reducedfluxstd is not None :
          self.stdflux    = PlotableStdTabRedFlux(reducedfluxstd, reducedstd,
                                                  frames["MOS_SPECPHOT_TABLE"])
          
    # This function creates all the subplots. It is responsible for the plotting 
    # layouts. 
    # There can different layouts, depending on the availability of data
    # Note that subplot(I,J,K) means the Kth plot in a IxJ grid 
    # Note also that the last one is actually a box with text, no graphs.
    def addSubplots(self, figure):
      if self.specphot is not None and self.stdflux is not None :
        self.subplot_std_extracted = figure.add_subplot(3,1,1)
        self.subplot_resp          = figure.add_subplot(3,1,2)
        self.subplot_std_tabulated = figure.add_subplot(3,1,3)
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
      if self.specphot is not None and self.stdflux is not None :

        #Extracted standard star
        self.plotStdExtracted()

        #Response
        self.plotResponse()

        #Tabulated data standard star
        self.plotStdTabulated()

      else :
        #Data not found info
        self.showNoData()
  
    def plotStdExtracted(self) :
        title_std_extracted = 'Extracted standard star spectrum'
        tooltip_std_extracted ="""Extracted standard star spectrum."""
        self.stdflux.plotStdRed(self.subplot_std_extracted, 
                                title_std_extracted, tooltip_std_extracted)

    def plotResponse(self) :
        if self.specphot.flat_sed :
          title_resp   = 'Raw and fitted Response (Flat SED correction applied)'
        else :
          title_resp   = 'Raw and fitted Response'
        tooltip_resp = """Raw and fitted Response. Points in green have been used to compute the fit"""
        self.specphot.plotResponse(self.subplot_resp, title_resp,
                               tooltip_resp)

    def plotStdTabulated(self) :
        title_std_tabulated = 'Flux calibrated and tabulated standard star spectrum'
        tooltip_std_tabulated ="""Tabulated standard star spectrum. Points in green have been used to compute the fit"""
        self.stdflux.plotStdTabRedFlux(self.subplot_std_tabulated, 
                                       title_std_tabulated, tooltip_std_tabulated)

    def showNoData(self) :
      #Data not found info
      self.subtext_nodata.set_axis_off()
      self.text_nodata = 'Response data and calibrated std not found in the products'
      self.subtext_nodata.text(0.1, 0.6, self.text_nodata, color='#11557c', fontsize=18,
                               ha='left', va='center', alpha=1.0)
      self.subtext_nodata.tooltip='Response data not found in the products'

    # This function specifies which are the parameters that should be presented
    # in the window to be edited.
    # Note that the parameter has to be also in the in_sop port (otherwise it 
    # won't appear in the window) 
    # The descriptions are used to show a tooltip. They should match one to one
    # with the parameter list 
    # Note also that parameters have to be prefixed by the 'recipe name:'
    def setInteractiveParameters(self):
      paramList = list()
      paramList.append(reflex.RecipeParameter('vmmosscience','resp_fit_nknots',group='Response',description='Number of knots in spline fitting of the instrument response'))
      paramList.append(reflex.RecipeParameter('vmmosscience','resp_fit_degree',group='Response',description='Degree for polynomial fitting of the instrument response'))
      paramList.append(reflex.RecipeParameter('vmmosscience','resp_ignore_mode',group='Response',description='Types of lines/regions to ignore in response. Valid ones are stellar_absorption, telluric and command_line (from parameter resp_ignore_points'))
      paramList.append(reflex.RecipeParameter('vmmosscience','resp_ignore_points',group='Response',description='Extra lines/regions to ignore in response. Use a comma separated list of values. A range can be specified like 4500.0-4600.0'))
      paramList.append(reflex.RecipeParameter('vmmosscience','resp_use_flat_sed',group='Response',description='Use the flat SED to normalise the observed spectra. Value are true, false, grism_table'))
      paramList.append(reflex.RecipeParameter('vmmosscience','resp_shift',group='Response',description='The extracted standard star will be shifted these many angstroms before using it to compute the response.  This is useful for observed std stars not centered in the slits. Positive values will shift the spectrum to the red. Shift is given in Angstroms but no fraction of pixels will be shifted'))
      paramList.append(reflex.RecipeParameter('vmmosscience','skyglobal',group='Sky subtraction',description='Subtract global sky spectrum from CCD'))
      paramList.append(reflex.RecipeParameter('vmmosscience','skymedian',group='Sky subtraction',description='Sky subtraction from extracted slit spectra'))
      paramList.append(reflex.RecipeParameter('vmmosscience','skylocal',group='Sky subtraction',description='Sky subtraction from CCD slit spectra'))
      paramList.append(reflex.RecipeParameter('vmmosscience','cosmics',group='Sky subtraction',description='Eliminate cosmic rays hits (only if global sky subtraction is also requested)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','slit_margin',group='Spectra extraction',description='Number of pixels to exclude at each slit in object detection and extraction.'))
      paramList.append(reflex.RecipeParameter('vmmosscience','ext_radius',group='Spectra extraction',description='Maximum extraction radius for detected objects (pixel)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','cont_radius',group='Spectra extraction',description='Minimum distance at which two objects of equal luminosity do not contaminate each other (pixel)'))
      paramList.append(reflex.RecipeParameter('vmmosscience','ext_mode',group='Spectra extraction',description='Object extraction method: 0 = aperture, 1 = Horne optimal extraction'))
      paramList.append(reflex.RecipeParameter('vmmosscience','skyalign',group='Wave calib',description='Polynomial order for sky lines alignment, or -1 to avoid alignment'))


      return paramList

    def setWindowHelp(self):
      help_text = """
In this window, the user will interact with the VIMOS response computation"""
      return help_text

    def setWindowTitle(self):
      title = 'Vimos Interactive Response Reduction'
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
