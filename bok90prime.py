#!/usr/bin/env python

'''
Telescope setting file for 90prime telescope.
Based on css.py file, to push to SAGUARO-MMA/Pipeline
'''

__version__ = "1.0" #last updated 29/01/2022

import numpy as np
import datetime
from astropy.io import fits
import os
from astropy.utils.exceptions import AstropyWarning
from slack import WebClient
import warnings
import glob
import subprocess
import gc
import saguaro_pipe
import ccdproc
from astropy import units as u
from astropy.nddata import CCDData
from astropy.modeling import models


warnings.simplefilter('ignore', category=AstropyWarning)
gc.enable()

def read_path(date):
    '''
    Returns the absolute path where raw files are stored.
    '''    
    read_date = datetime.datetime.strptime(date,'%Y/%m/%d').strftime('%Y/%y%b%d')
    return '/home/saguaro/data/90prime/raw/'+read_date

def write_path():
    '''
    Returns the absolute path where reduced files are written to.
    '''    
    return '/home/saguaro/data/90prime/'

def work_path(date):
    '''
    Returns the working directory for the pipeline.
    '''
    return '/home/saguaro/data/90prime/tmp/'+date+'/'

def log_path():
    '''
    Returns the relative path where log files are written to.
    '''
    return 'log/'

def red_path(date):
    '''
    Returns the relative path where reduced science images are written to.
    '''
    return 'red/'+date+'/'

def file_name():
    '''
    Returns name pattern of input files.
    '''
    return 'G96_*_med.fits.fz'

def create_mask():
    '''
    Returns if a mask should be created for each image.
    '''
    return True

def bad_pixel_mask():
    '''
    Returns the full path of the bad pixel mask.
    '''
    return '/home/saguaro/software/Pipeline/90prime_bpm.fits'

def fieldID(header):
    '''
    Returns field ID of image.
    '''
    fieldID = header['OBJECT']
    return fieldID

def saturation():
    '''
    Returns the saturation level in electrons.
    " will change these later"
    '''
    return 57000

def binning():
    '''
    Returns the image binning used during the determination of the satellite trail mask.
    " will change these later"
    '''
    return 2

def slack_client():
    '''
    Returns the slackclient of the saguaro pipeline channel.
    '''
    with open('/home/saguaro/software/saguaro_slack.txt','r') as f:
        slack_token = f.readline().rstrip()
    return WebClient(token=slack_token)

def mask_edge_pixels(data,mask):
    '''
    Returns the image with edge pixels masked.
    '''
    return data

def mask_bp():
    '''
    Returns if bad pixels were masked.
    '''
    return 'F'

def mask_cr():
    '''
    Returns if cosmic rays were masked.
    '''
    return 'T'

def mask_sp():
    '''
    Returns if saturated pixels were masked.
    '''
    return 'T'

def mask_scp():
    '''
    Returns if pixels connected to saturated pixles were masked.
    '''
    return 'T'

def mask_sat():
    '''
    Returns if satellite trails were masked.
    '''
    return 'F'

def mask_ep():
    '''
    Returns if edge pixels were masked.
    '''
    return 'T'

def cosmic():
    '''
    Returns if cosmic ray correction should be done.
    '''
    return 'F'

def sat():
    '''
    Returns if satellite trail correction should be done.'''
    return 'T'

def sat_readnoise(header):
    '''
    Returns the readnoise value from the header.
    '''
    return header['RDNOISE']

def sat_gain(header):
    '''
    Returns the gain value from the header.
    '''
    return header['GAIN']

def sat_objlim():
    '''
    Returns the value used in the mask creation for determining saturated stars.
    '''
    return 20

def sat_sigclip():
    '''
    Returns the vaule used for sigma clipping in the satellite trail fitting.
    '''
    return 6

def sat_repeat():
    '''
    Return how many times the satellite fitting should be done.
    '''
    return 0

def time_zone():
    '''
    Returns the time zone of the pipeline.
    '''
    return -7

def tel_zone():
    '''
    Returns the adjustment to the time zone for the telescope.
    '''
    return 0

def tel_delta():
    '''
    Returns the adjustment to the time zone for the telescope.
    '''
    return 0

# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

masterbias_string = 'allframe/biases/mbias_90prime.fits'

def remove_bias(f):
    ''' reduce using ccdproc (bias, trim, readnoise, gain, overscan),
        then write to 4 files (1 for each ccd) with names '...nobias.fits'
    '''

    # open file and save information from header
    hdul = fits.open(f)
    extnames = [int(hdul[i].header['EXTNAME'][2:]) for i in range(1,17)]
    gains = [float(hdul[0].header['GAIN'+str(extnames[i])]) for i in range(16)]
    readnoises = [float(hdul[0].header['RDNOIS'+str(extnames[i])]) for i in range(16)]
    mbias = [CCDData.read(masterbias_string, hdu=x+1, unit=u.electron) for x in range(16)]
    
    # open and reduce all 16 extensions using ccdproc, add reduced arrays to red_data list
    raw = [CCDData.read(f, hdu=x+1) for x in range(16)] #removed unit='adu' bc of error
    red_data = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'],oscan_model=models.Chebyshev1D(3), 
                               trim=x.header['DATASEC'], gain=gains[k]*u.electron/u.adu, 
                               readnoise=readnoises[k]*u.electron, master_bias=mbias[k],
                                    gain_corrected=True).data for k,x in enumerate(raw)]

    # lists of each extension's x and y ranges (in context of full image)
    yranges = [hdul[i].header['DETSEC'][1:-1].split(',')[0] for i in range(1,len(hdul))]
    xranges = [hdul[i].header['DETSEC'][1:-1].split(',')[1] for i in range(1,len(hdul))]    
    
    # flip each extension according to xranges and yranges
    for i in range(len(yranges)):
        y1, y2 = yranges[i].split(':')
        x1, x2 = xranges[i].split(':')
        if y1 > y2 and x1 > x2:
            red_data[i] = np.flip(red_data[i].T, 1)[::-1]
        if y1 < y2 and x1 > x2:
            red_data[i] = np.flip(red_data[i].T, 1)
        if y1 > y2 and x1 < x2:
            red_data[i] = red_data[i].T[::-1]
        if y1 < y2 and x1 < x2:
            red_data[i] = red_data[i].T

    # generate a large empty array to fit all 16 extensions
    emptysize = hdul[0].header['DETSIZE'][1:-1].split(',')
    empty = np.zeros(shape=(int(emptysize[0][2:]),int(emptysize[1][2:])))
    # create lists with each extension's lower bound on x/y ranges
    lowys = [int(min(yranges[i].split(':'))) for i in range(len(yranges))]
    lowxs = [int(min(xranges[i].split(':'))) for i in range(len(xranges))]
    
    # assign each extension's data to the correct section of empty array
    for i in range(len(red_data)):
        highy = lowys[i] + int(hdul[1].header['DATASEC'][1:-1].split(',')[0][2:])-1
        highx = lowxs[i] + int(hdul[1].header['DATASEC'][1:-1].split(',')[1][2:])-1
        empty[lowys[i]-1:highy, lowxs[i]-1:highx] = red_data[i]
    
    # split the large empty array into 4 arrays (one for each CCD, "quad data")
    ycut, xcut = int(int(emptysize[0][2:])/2), int(int(emptysize[1][2:])/2)
    qdata = [empty[:ycut,:xcut].T, empty[:ycut,xcut:].T, 
             empty[ycut:,:xcut].T, empty[ycut:,xcut:].T]
    
    # create 4 headers (one for each CCD, "quad headers")
    qheads = [hdul[0].header.copy(), hdul[0].header.copy(), hdul[0].header.copy(), hdul[0].header.copy()]
    extheads = [hdul[4].header.copy(), hdul[9].header.copy(), hdul[8].header.copy(), hdul[13].header.copy()]
    kw_list = ['EQUINOX','WCSDIM','CTYPE1','CTYPE2','CRVAL2',
        'CRVAL1','CRPIX1','CRPIX2', 'CD1_1','CD1_2','CD2_1','CD2_2']
    for h in kw_list:
        for i in range(len(qheads)):      
            qheads[i][h] = extheads[i][h]
            qheads[i]['DATASEC']='[1:4032,1:4096]'
    
    # write to 4 fits files, nobias
    for i in range(4):
        fits.writeto(f.replace('.fits','_q'+str(i+1)+'nobias.fits'), qdata[i], qheads[i], overwrite=True)

    return

def flatfield(f, masterflat_string):
    ''' takes in a 'nobias' file and writes a 'red.fits' file, one quadrant (ccd)
    '''
    master_flat = CCDData(fits.open(masterflat_string)[0].data, unit=u.electron)
    opened = fits.open(f)[0]
    nobias = CCDData(opened.data, header=opened.header, unit=u.electron)
    red = ccdproc.flat_correct(nobias, master_flat)
    fits.writeto(f.replace('nobias.fits', 'red.fits'), red,
                            nobias[0].header, overwrite=True) # 'red' as in reduced
    return

def science_process(science_file):
    '''Function to process science images. In two steps: remove bias, flatfield'''

    # open file, determine masterflat string (filter-dependent)
    hdul = fits.open(science_file)
    filt = hdul[0].header['FILTER']
    masterflat_quads = ['allframe/flat'+filt+'/mflat'+str(i+1)+'_90prime.fits' for i in range(4)]

    # reduce!
    print('removing bias...')
    remove_bias(science_file)
    nobiasfiles = [science_file.replace('.fits', '_q'+str(i+1)+'nobias.fits') for i in range(4)]

    for i in range(4):
        print('flatfielding quad no.', (i+1))
        flatfield(nobiasfiles[i], masterflat_quads[i])

    reduced_quads = [science_file.replace('.fits', '_q'+str(i+1)+'red.fits') for i in range(4)]
    comment = 'data is reduced!'

    return reduced_quads, comment

# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

def output(f,date):
    '''
    Return the full path of the output used to determine if a file has been processed.
    '''
    return os.path.exists(write_path()+red_path(date)+os.path.basename(f.replace('.fits','_trans.fits')))

def find_ref(reduced):
    '''
    Function to find a reference for image subtraction.
    '''
    with fits.open(reduced) as hdr:
        header = hdr[0].header
    field = fieldID(header)
    ref_files = glob.glob(write_path()+'ref/'+field+'/*')
    for f in ref_files:
        subprocess.call(['cp',f,'.'])
        if '.fz' in f:
            f = saguaro_pipe.funpack_file(os.path.basename(f))
    ref_file = field+'_wcs.fits'
    if os.path.exists(ref_file):
        return ref_file
    else:
        return None
