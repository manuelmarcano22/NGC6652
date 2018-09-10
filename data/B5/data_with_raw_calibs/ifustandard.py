from astropy.io import fits
import glob, os



files = glob.glob('*fits')

std = []
for f in files:
    hd = fits.getheader(f)
    try:
        if hd['ESO DPR TYPE']=='STD':
            std.append(f)
    except:
    	pass

sof = open('ifustandard.sof', 'w')

for i in std:
    sof.write('./'+str(i)+'\tIFU_STANDARD\n')


sof.write('master_bias.fits\tMASTER_BIAS\n')
sof.write('ifu_ids.fits\tIFU_IDS\n')
sof.write('ifu_trace.fits\tIFU_TRACE\n')
sof.write('ifu_transmission.fits\tIFU_TRANSMISSION\n')
sof.write('/home/mmarcano/Documents/VIMOS/NGC6652/esorex/calibrations/calib/vimos-3.2.3/cal/extinct_table.tfits\tEXTINCT_TABLE\n')
sof.write('/home/mmarcano/Documents/VIMOS/NGC6652/esorex/calibrations/calib/vimos-3.2.3/cal/ltt4816.tfits\tSTD_FLUX_TABLE\n')

sof.close()
