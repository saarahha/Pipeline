#!/usr/bin/env python

'''
Telescope setting file for CSS 1.5m Mt Lemmon telescope.
'''

__version__ = "1.0" #last updated 28/09/2021

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

warnings.simplefilter('ignore', category=AstropyWarning)
gc.enable()

def read_path(date):
    '''
    Returns the absolute path where raw files are stored.
    '''    
    read_date = datetime.datetime.strptime(date,'%Y/%m/%d').strftime('%Y/%y%b%d')
    return '/home/saguaro/data/css/raw/'+read_date

def write_path():
    '''
    Returns the absolute path where reduced files are written to.
    '''    
    return '/home/saguaro/data/css/'

def work_path(date):
    '''
    Returns the working directory for the pipeline.
    '''
    return '/home/saguaro/data/css/tmp/'+date+'/'

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
    return '/home/saguaro/software/Pipeline/css_bpm.fits'

def fieldID(header):
    '''
    Returns field ID of image.
    '''
    fieldID = header['OBJECT']
    return fieldID

def saturation():
    '''
    Returns the saturation level in electrons.
    '''
    return 57000

def binning():
    '''
    Returns the image binning used during the determination of the satellite trail mask.
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

def science_process(science_file,unique_dir,log_file_name):
    '''
    Function to process science images. As CSS data has already been processed, so only a mask is created.
    '''
    subprocess.call(['cp',science_file,science_file.replace('.fits','_mask.fits'),'.'])
    science_file = saguaro_pipe.funpack_file(os.path.basename(science_file))
    mask_file = saguaro_pipe.funpack_file(os.path.basename(science_file).replace('.fits','_mask.fits'))
    with fits.open(science_file) as hdr:
        Red = hdr[0].data
        header = hdr[0].header
    comment = 'No reduction needed. Creating mask.\n '
    with fits.open(mask_file) as hdr:
        mask_bp = hdr[0].data
    Red, header, comment = saguaro_pipe.mask_create(science_file,'css',unique_dir,Red,mask_bp,header,comment,saturation(),log_file_name)
    fits.writeto(science_file,Red,header,overwrite=True)
    return science_file, comment

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
