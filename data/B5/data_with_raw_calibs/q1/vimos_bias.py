from astropy.io import fits
import glob, os


files = glob.glob('*.fits')

bias = []
for f in files:
    hd = fits.getheader(f)
    if hd['ESO DPR TYPE']=='BIAS':
            bias.append(f)



sof = open('bias.sof', 'w')        
for i in bias:
    sof.write('./'+str(i)+'\tBIAS\n')


sof.close()