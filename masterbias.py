# astropy
from astropy.io import fits
from astropy import units as u
from astropy.nddata import CCDData
from astropy.modeling import models
# others
import os
import glob
import ccdproc

#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#

bias_list_default = glob.glob('allframe/biases/frame*.fits', recursive = True)

def create_bias(path, bias_list=bias_list_default):
    '''create masterbias file for 90prime
         - bias_list_default: bias frames from 90prime box folder'''
    amp = 16 # number of amps on 90prime
    extnames = [int(fits.open(bias_list[0])[i].header['EXTNAME'][2:]) for i in range(1,17)]
    gains = [float(fits.open(bias_list[0])[0].header['GAIN'+str(extnames[i])]) for i in range(16)]
    readnoises = [float(fits.open(bias_list[0])[0].header['RDNOIS'+str(extnames[i])]) for i in range(16)]

    for i,bias in enumerate(bias_list):
        with fits.open(bias) as hdr:
            header = hdr[0].header
        raw = [CCDData.read(bias, hdu=x+1) for x in range(amp)]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], 
                                   oscan_model=models.Chebyshev1D(3), 
                                   trim=x.header['DATASEC'], 
                                   gain=gains[j]*u.electron/u.adu, 
                                   readnoise=readnoises[j]*u.electron) for j,x in enumerate(raw)]
        bias_hdu = fits.HDUList([fits.PrimaryHDU(header=header)])
        for x in red: bias_hdu.append(fits.ImageHDU(x.data,header=x.header))
        bias_hdu.writeto(bias.replace('.fits','_red.fits'),overwrite=True)
        bias_list[i] = bias.replace('.fits','_red.fits')
    mbias = [ccdproc.combine(bias_list,hdu=x+1,unit=u.electron) for x in range(amp)]
    mbias_hdu = fits.HDUList([fits.PrimaryHDU()])
    #for x in mbias: mbias_hdu.append(fits.ImageHDU(x.data.T,header=x.header)) 
    #mbias_hdu.writeto(path+'mbias_90prime.T.fits',overwrite=True)
    for x in mbias: mbias_hdu.append(fits.ImageHDU(x.data,header=x.header))
    mbias_hdu.writeto(path+'mbias_90prime.fits',overwrite=True)
    print('master bias created')
    #possible: rearrange amps so [2] and [3] are reversed......
    for bias in bias_list: os.remove(bias)
    return