from astropy.io import fits
import glob, os



files = glob.glob('*fits')

obj = []
for f in files:
    hd = fits.getheader(f)
    try:
        if hd['ESO DPR TYPE']=='OBJECT':
            obj.append(f)
    except:
        pass


sof = open('ifuscience.sof', 'w')

#only one science is required
for i in obj[0:1]:
    sof.write('./'+str(i)+'\tIFU_SCIENCE\n')


sof.write('master_bias.fits\tMASTER_BIAS\n')
sof.write('ifu_ids.fits\tIFU_IDS\n')
sof.write('ifu_trace.fits\tIFU_TRACE\n')
sof.write('ifu_transmission.fits\tIFU_TRANSMISSION\n')
sof.write('/home/mmarcano/Documents/VIMOS/NGC6652/esorex/calibrations/calib/vimos-3.2.3/cal/extinct_table.tfits\tEXTINCT_TABLE\n')
sof.write('ifu_specphot_table.fits\tIFU_SPECPHOT_TABLE\n')

sof.close()
