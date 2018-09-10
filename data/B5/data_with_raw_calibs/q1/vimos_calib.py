from astropy.io import fits
import glob, os



files = glob.glob('*fits')


flat = []
for f in files:
    hd = fits.getheader(f)
    try:
        if hd['ESO DPR TYPE']=='FLAT,LAMP':
            flat.append(f)
    except:
        pass

arc = []
for f in files:
    hd = fits.getheader(f)
    try:
        if hd['ESO DPR TYPE']=='WAVE,LAMP':
            arc.append(f)
    except:
        pass

arc = []
for f in files:
    hd = fits.getheader(f)
    try:
        if hd['ESO DPR TYPE']=='WAVE,LAMP':
            arc.append(f)
    except:
        pass








sof = open('calib.sof', 'w')        
for i in flat:
    sof.write('./'+str(i)+'\tIFU_SCREEN_FLAT\n')

        
#Only one flat needed
for i in arc[-1:]:
    sof.write('./'+str(i)+'\tIFU_ARC_SPECTRUM\n')




sof.write('master_bias.fits\tMASTER_BIAS\n')
sof.write('/home/mmarcano/Documents/VIMOS/NGC6652/esorex/calibrations/calib/vimos-3.2.3/cal/lcat_HR_blue.tfits\tLINE_CATALOG\n')
sof.write('/home/mmarcano/Documents/VIMOS/NGC6652/esorex/calibrations/calib/vimos-3.2.3/cal/ifu_ident_HR_blue.1.fits\tIFU_IDENT\n')

sof.close()




