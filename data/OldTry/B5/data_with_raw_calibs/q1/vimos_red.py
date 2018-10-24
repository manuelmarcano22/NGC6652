from astropy.io import fits
import glob, os, shutil

files = glob.glob('VIMOS*.fits')

quadrants = ['q1', 'q2', 'q3', 'q4']

for i in quadrants:                                                                                                                                                                                                 
      os.mkdir(str(i)) 


for f in files:
	hd = fits.getheader(f)
	if hd['ESO OCS CON QUAD'] == 1:
		shutil.move(str(f),"q1/")
	if hd['ESO OCS CON QUAD'] == 2:
		shutil.move(str(f),"q2/")
	if hd['ESO OCS CON QUAD'] == 3:
		shutil.move(str(f),"q3/")
	if hd['ESO OCS CON QUAD'] == 4:
		shutil.move(str(f),"q4/")


