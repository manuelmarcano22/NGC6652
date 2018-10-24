import glob
import os
import re
import itertools
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import numpy as np
import subprocess

def flatten(listOfLists):
    "Flatten one level of nesting"
    return itertools.chain.from_iterable(listOfLists)

def reduceifu(datadirectory, rawdatadictionary):
    """Something"""
    datadir = datadirectory
    rawdata = rawdatadictionary



    # Put to each the data dir and print hte quadran
    for k,v in rawdata.items():
        try:
            rawdata[k] = datadir+v
        except:
            temp = [datadir+i for i in v]
            rawdata[k] = temp


    #Get the quadrant of SCIENCE
    t = rawdata['IFU_SCIENCE']
    if isinstance(t,list):
        t2 = fits.open(t[0])
        quadrant = t2[0].header['HIERARCH ESO OCS CON QUAD']

    print('Quadrant {}'.format(quadrant))

    
    
    ##Create bias
    
    outdir = datadir+'quadrant'+str(quadrant)

    #Create outdir if it doesnt exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    filenamebias = outdir+'/bias.sof'
    with open(filenamebias,'w') as file:
        for bias in rawdata['BIAS'][0:]:
            cat = 'BIAS'
            string = '{} {}\n'.format(bias,cat)
            file.write(string)
            
            
                
    stackmethod = '--StackMethod=Average'
    #esorexpath = '/bin/esorex'
    esorexpath = '/home/mmarcano/Documents/VIMOS/NGC6652/esorex/bin/esorex'
    routine = 'vmbias'
    outputdir = '--output-dir={}'.format(outdir)
    logdir = '--log-dir={}'.format(outdir)
    
    commandbias = '{} {} {} {}  {}'.format(esorexpath,outputdir,logdir,routine,  filenamebias)
    print(commandbias)
    #os.system(commandbias)
    subprocess.run(commandbias,shell=True,check=True)
    
    os.rename(outdir+'/esorex.log',outdir+'/'+routine+'.log')
    
    rawdata['MASTER_BIAS'] = outdir+'/master_bias.fits'
    
    
    ##Calib
    
    inputfilescalib = ['IFU_SCREEN_FLAT','IFU_ARC_SPECTRUM','LINE_CATALOG','IFU_IDENT','MASTER_BIAS']

    filenamecalib = outdir+'/calib.sof'
    with open(filenamecalib,'w') as file:
        for i in inputfilescalib:
            tmpfile = rawdata[i]
            if isinstance(tmpfile, list):
                if i == 'IFU_ARC_SPECTRUM':
                    string = '{} {}\n'.format(tmpfile[0],i)
                    file.write(string)
                else:
                    for j in tmpfile:
                        string = '{} {}\n'.format(j,i)
                        file.write(string)
            else:
                string = '{} {}\n'.format(tmpfile,i)
                file.write(string)


    # Define the routine and call it
    routine = 'vmifucalib'

    commandcalib = '{} {} {} {}  {}'.format(esorexpath,outputdir,logdir,routine,  filenamecalib)
    print(commandcalib)
    subprocess.run(commandcalib,shell=True)
    
    os.rename(outdir+'/esorex.log',outdir+'/'+routine+'.log')

    rawdata['IFU_IDS'] = outdir+'/ifu_ids.fits'
    rawdata['IFU_TRACE'] = outdir+'/ifu_trace.fits'
    rawdata['IFU_TRANSMISSION'] = outdir+'/ifu_transmission.fits'
    
    
    ##Standard
    inputfilesstandard = ['IFU_STANDARD','MASTER_BIAS','IFU_IDS','IFU_TRACE','IFU_TRANSMISSION','EXTINCT_TABLE','STD_FLUX_TABLE']

    filenamestandard = outdir+'/ifustandard.sof'
    with open(filenamestandard,'w') as file:
        for i in inputfilesstandard:
            tmpfile = rawdata[i]
            if isinstance(tmpfile, list):    
                for j in tmpfile:
                    string = '{} {}\n'.format(j,i)
                    file.write(string)
            else:
                string = '{} {}\n'.format(tmpfile,i)
                file.write(string)

                routine = 'vmifustandard'

    commandstandard = '{} {} {} {}  {}'.format(esorexpath,outputdir,logdir,routine,  filenamestandard)
    print(commandstandard)
    subprocess.run(commandstandard,shell=True)
    os.rename(outdir+'/esorex.log',outdir+'/'+routine+'.log')

    #Science Soif
    
    rawdata['IFU_SPECPHOT_TABLE'] = outdir+'/ifu_specphot_table.fits'

    inputfilesscience = ['IFU_SCIENCE','MASTER_BIAS','IFU_IDS','IFU_TRACE','IFU_TRANSMISSION','EXTINCT_TABLE','IFU_SPECPHOT_TABLE']

    filenamescience = outdir+'/ifuscience.sof'
    with open(filenamescience,'w') as file:
        for i in inputfilesscience:
            tmpfile = rawdata[i]
            if isinstance(tmpfile, list):    
                for j in tmpfile:
                    string = '{} {}\n'.format(j,i)
                    file.write(string)
            else:
                string = '{} {}\n'.format(tmpfile,i)
                file.write(string)
                
                
                routine = 'vmifuscience'

    calibrateflux = '--CalibrateFlux=true'

    commandscience = '{} {} {} {} {} {}'.format(esorexpath,outputdir,logdir,routine, calibrateflux ,filenamescience)
    print(commandscience)
    subprocess.run(commandscience,shell=True)
    os.rename(outdir+'/esorex.log',outdir+'/'+routine+'.log')

    
    
    
def reducecombinecube(datadirectory,typecombine='ifu_science_flux_reduced',quadrants=[1,2,3,4]):
    """Or it can be 'IFU_SCIENCE_REDUCED' """
    #cubes = ['ifu_science_reduced','ifu_science_flux_reduced']
    datadir = datadirectory
    
    
    combinedir = datadir+'Combine'

    #Create outdir if it doesnt exist
    if not os.path.exists(combinedir):
        os.makedirs(combinedir)

    filenamescombinecube = combinedir+'/ifucombinefov{}.sof'.format(typecombine.lower())

    with open(filenamescombinecube,'w') as file:
        for i in quadrants:
            fitstemp='{}quadrant{}/{}.fits'.format(datadir,i,typecombine.lower())
            string = '{} {}\n'.format(fitstemp,typecombine.upper())
            file.write(string)

    #esorexpath = '/bin/esorex'
    esorexpath = '/home/mmarcano/Documents/VIMOS/NGC6652/esorex/bin/esorex'
    outputdir = '--output-dir={}'.format(combinedir)
    logdir = '--log-dir={}'.format(combinedir)
    routine = 'vmifucombinecube'
    command = '{} {} {} {} {} '.format(esorexpath,outputdir,logdir,routine, filenamescombinecube)
    print(command)
    subprocess.run(command,shell=True)
    
def changewcs(original,destination):
    originalcube = original
    newcube = destination
    
    copyfile = 'cp {} {}'.format(originalcube,newcube)
    os.system(copyfile)
    fits.info(newcube)
    oricube = fits.open(originalcube)
    oricubehearder = oricube[0].header
    ifura = oricubehearder['HIERARCH ESO INS IFU RA']
    ifudec = oricubehearder['HIERARCH ESO INS IFU DEC']
    
    fits.setval(newcube, 'CRVAL1', value=ifura)
    fits.setval(newcube, 'CRVAL2', value=ifudec)