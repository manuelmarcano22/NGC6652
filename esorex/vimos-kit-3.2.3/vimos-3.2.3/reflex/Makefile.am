
AUTOMAKE_OPTIONS = foreign

#Put here all the workflows. They must end with extension .in. They should
#be like a normal workflow but using @prefix@ for paths 
#to the OCA rules (grep OCA in the xml file) or Python scripts
WORKFLOWS = vimos_ifu.xml vimos_mos.xml
#Put here the OCA rules used by the workflow
OCAWKF =  vimos_ifu_wkf.oca vimos_ifu_wkf.dvd.oca vimos_mos_wkf.oca
#Put Python scripts used by the workflows
PYTHONWKF = vimos_calib_interact.py vimos_science_interact.py vimos_response_interact.py vimos_plot_common.py

#This installs the workflow and the OCA rules in reflexaux directory
wkfextra_DATA = $(WORKFLOWS) $(OCAWKF) $(PYTHONWKF)

#Pack in the distribution the original workflow
EXTRA_DIST = $(WORKFLOWS).in $(OCAWKF) $(PYTHONWKF)

#This also installs the workflow in the common reflex directory
wkfcopy_DATA = $(WORKFLOWS)

