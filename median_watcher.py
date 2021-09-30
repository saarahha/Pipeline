#!/usr/bin/env python

'''
Script to create median images from the 4 CSS images per field.
'''

__version__ = "1.0" #last updated 30/09/2021

import sys
import numpy as np
import os
import datetime
import time
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
import warnings
import glob
import subprocess
import gc
import uuid
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import fnmatch as fn
import css
import shutil
import logging
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import requests
import saguaro_pipe

warnings.simplefilter('ignore', category=AstropyWarning)
gc.enable()

def check_field(event):
    '''
    Check if this field has already been processed.
    '''
    try:
        try:
            file = str(event.dest_path)
        except AttributeError:
            file = str(event.src_path) #get name of new file
        if fn.fnmatch(os.path.basename(file),'G96_*.calb.fz'): #only continue if event is a fits file
            saguaro_pipe.copying(file) #check to see if write is finished writing
    except AttributeError: #if event is a file
        file = event
    field = file.split('_2B_')[1].split('_00')[0]
    if field not in field_list:
        field_list.append(field)
        if 'N' in field or 'S' in field:
            return field

def action(event,date,read_path,write_path,field):
    '''
    Waits for all 4 images (max time 5 mins) to create a median image.
    '''
    t1 = datetime.datetime.utcnow()
    while True:
        files = sorted(glob.glob(read_path+'/G96_*'+field+'*.calb.fz'))
        if len(files) == 4:
            break
        else:
            t2 = datetime.datetime.utcnow()
            if (t2-t1).total_seconds() > 300:
                break
            else:
                time.sleep(1)
    t1 = datetime.datetime.utcnow()
    while True:
        headers = sorted(glob.glob(read_path+'/G96_*'+field+'*.arch_h'))
        if len(headers) == len(files):
            break
        else:
            t2 = datetime.datetime.utcnow()
            if (t2-t1).total_seconds() > 300:
                break
            else:
                time.sleep(1)
    logger.info('Number of files for field = '+str(len(files)))
    combine = []
    mjd = []
    back = []
    zp = []
    out_file = os.path.basename(files[0]).split('_00')[0]+'_med.fits'
    unique_dir = '/home/saguaro/data/css/tmp/'+date+'/'+uuid.uuid1().hex+'/'
    os.makedirs(unique_dir)
    for i,f in enumerate(files):
        subprocess.call(['cp',f,unique_dir])
        try:
            header = fits.open(f.replace('calb.fz','arch_h'))[0].header
            c = unique_dir+os.path.basename(f)
            con = True
        except:
            logging.error('No header available for '+f)
            con = False
        if con:
            with fits.open(c) as hdr:
                hdr.verify('fix+ignore')
                head0 = hdr[0].header
                head1 = hdr[1].header
                header['CTYPE1'] = 'RA---TPV'
                header['CTYPE2'] = 'DEC--TPV'
                data = hdr[1].data
                hdul = fits.HDUList([fits.PrimaryHDU(data,header)])
                hdul.writeto(c.replace('.calb.fz','.fits'),output_verify='fix+ignore')
                try:
                    stars = header['WCSMATCH'] #check image type
                    t = header['MJD']
                    if stars > 100: #only continue if science image
                         combine.append(c.replace('.calb.fz','.fits'))
                         mjd.append(t)
                         back.append(np.median(data))
                         zp.append(header['MAGZP'])
                    else:
                        bad_images += 1
                except:
                    logging.critical('Error with file '+f)
                with fits.open('/home/saguaro/software/css_bpm.fits') as hdr:
                    mask_header = hdr[0].header
                    data = hdr[0].data
                    fits.writeto(c.replace('.calb.fz','_mask.fits'),data,mask_header+header,output_verify='fix+ignore')
                with open(unique_dir+out_file.replace('.fits','.head'),'w') as swarp_head:
                    for card in header.cards:
                        swarp_head.write(str(card)+'\n')
                shutil.copy(unique_dir+out_file.replace('.fits','.head'),unique_dir+out_file.replace('.fits','_mask.head'))
    if len(combine) > 1:    
        masks = [x.replace('.fits','_mask.fits') for x in combine]
        subprocess.call(['swarp']+combine+['-c']+[os.environ['ZOGYHOME']+'/Config/swarp_css.config']+['-IMAGE_SIZE']+['5280,5280']+['-IMAGEOUT_NAME']+[unique_dir+out_file]+['-SUBTRACT_BACK']+['YES']+['-GAIN_KEYWORD']+['GAIN']+['-BACK_SIZE']+['256']+['-BACK_FILTERSIZE']+['3']+['-FSCALASTRO_TYPE']+['VARIABLE']+['-FSCALE_KEYWORD']+['FLXSCALE'])
        subprocess.call(['swarp']+masks+['-c']+[os.environ['ZOGYHOME']+'Config/swarp_css.config']+['-IMAGE_SIZE']+['5280,5280']+['-IMAGEOUT_NAME']+[unique_dir+out_file.replace('.fits','_mask.fits')]+['-SUBTRACT_BACK']+['NO']+['-GAIN_DEFAULT']+['1']+['-COMBINE_TYPE']+['SUM'])
        with fits.open(unique_dir+out_file,mode='update') as hdr:
            header_swarp = hdr[0].header
            data = hdr[0].data
            data /= header_swarp['FLXSCALE']
            data += np.median(back)
            header_swarp['EXPTIME'] = 30
            header_swarp['RDNOISE'] = 11.6/np.sqrt(len(combine))
            header_swarp['GAIN'] = 3.1
            try:
                del header_swarp['CROTA1']
                del header_swarp['CROTA2']
            except KeyError:
                pass
            header_swarp['MJD'] = mjd[0]+((mjd[-1]-mjd[0])/2)
            header_swarp['NCOMBINE'] = (len(combine), 'Number of files used to create median.')
            for i,c in enumerate(combine):
                header_swarp['FILE'+str(i+1)] = (c)
        with fits.open(unique_dir+out_file.replace('.fits','_mask.fits')) as mask_hdr:
            mask_header = mask_hdr[0].header
            mask_data = mask_hdr[0].data
        mask_data = (mask_data+0.5).astype(np.uint8)
        mask_data[mask_data!=0] = 1
        fits.writeto(unique_dir+out_file.replace('.fits','_mask.fits'),mask_data,mask_header,overwrite=True)
        subprocess.call(['fpack','-D','-Y','-g','-q0','16',unique_dir+out_file])
        subprocess.call(['fpack','-D','-Y','-g',unique_dir+out_file.replace('.fits','_mask.fits')])
        subprocess.call(['mv',unique_dir+out_file+'.fz',write_path])
        subprocess.call(['mv',unique_dir+out_file.replace('.fits','_mask.fits')+'.fz',write_path])
        os.chdir('/home/saguaro/')
        subprocess.call(['rm','-r',unique_dir])
        print('Created '+out_file)
    ncombine.append(len(combine))
    logger.info('Number of files used in median = '+str(len(combine)))
        

class FileWatcher(FileSystemEventHandler,object):
    '''
    Monitors directory for new files.
    '''
    
    def __init__(self, date, read_path, write_path): #parameters needed for action
        self._date = date
        self._read_path = read_path
        self._write_path = write_path
        self._field_list = field_list

    def on_created(self, event):
        '''
        Action to take for new files.
        '''
        field = check_field(event)
        if field:
            logger.info('Found field '+field)
            action(event,self._date,self._read_path,self._write_path,field)

    def on_moved(self, event):
        '''
        Action to take for renamed files.
        '''
        field = check_field(event)
        if field:
            action(event,self._date,self._read_path,self._write_path,field)

global field_list, ncombine, logger, bad_images, missing_head

try:
    mode = sys.argv[1]
except:
    mode = False

try:
    date = sys.argv[2]
except:
    date = (datetime.datetime.utcnow()).strftime('%Y/%m/%d')

log_stream = StringIO() #create log stream for upload to slack
log_file_name = '/home/saguaro/data/log/median_watcher_'+(datetime.datetime.utcnow()).strftime('%Y%m%d_%H%M%S') #create log file name
log = logging.getLogger(log_file_name) #create logger
log.setLevel(logging.INFO) #set level of logger
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s") #set format of logger
logging.Formatter.converter = time.gmtime #convert time in logger to UCT
filehandler = logging.FileHandler(log_file_name+'.log', 'w+') #create log file
filehandler.setFormatter(formatter) #add format to log file
log.addHandler(filehandler) #link log file to logger
streamhandler_slack = logging.StreamHandler(log_stream) #add log stream to logger
streamhandler_slack.setFormatter(formatter) #add format to log stream
log.addHandler(streamhandler_slack) #link logger to log stream
logger = saguaro_pipe.MyLogger(log,log_stream,'css') #load logger handler

logger.critical('Median watcher started.')

read_path = '/home/data/css/G96/'+datetime.datetime.strptime(date,'%Y/%m/%d').strftime('%Y/%y%b%d')
read_dir = False
while read_dir is False:
    if not os.path.exists(read_path):
        print('waiting for directory to be created...')
        done = saguaro_pipe.scheduled_exit(datetime.datetime.utcnow(),'css')
        if done:
            logger.critical('Scheduled time reached.')
            logger.critical('No data ingested.')
            logging.shutdown()
            sys.exit()
        else:
            time.sleep(1)
    else:
        read_dir = True
write_path = '/home/saguaro/data/css/raw/'+datetime.datetime.strptime(date,'%Y/%m/%d').strftime('%Y/%y%b%d')
if not os.path.exists(write_path):
    os.makedirs(write_path)

field_list = []
ncombine = []
bad_images = 0
missing_head = 0

if mode: #to rerun all data
    files = glob.glob(read_path+'/G96*.calb.fz')
    for f in files:
        field = check_field(f)
        if field:
            action(f,date,read_path,write_path,field)
else: #normal real-time reduction
    observer = Observer()
    observer.schedule(FileWatcher(date,read_path,write_path), read_path, recursive=False)
    observer.start()
    while True:
        done = saguaro_pipe.scheduled_exit(datetime.datetime.utcnow(),'css')
        if done:
            logger.critical('Scheduled time reached.')
            observer.stop()
            observer.join()
            break
        else:
            time.sleep(1)

ncombine = np.asarray(ncombine)
logger.critical('Median watcher summary: %i fields observed, %i medians made with 4 images, %i medians made with 3 images, %i medians made with 2 images, %i medians made with 1 image, %i medians not made. %i images not used due to bad weather.' %(len(ncombine),len(ncombine[ncombine==4]),len(ncombine[ncombine==3]),len(ncombine[ncombine==2]),len(ncombine[ncombine==1]),len(ncombine[ncombine==0]),bad_images))
files_raw = len(glob.glob(read_path+'/G96*_N*.calb.fz'))+len(glob.glob(read_path+'/G96*_S*.calb.fz'))
files_sex = len(glob.glob(read_path+'/G96*_N*.sext.gz'))+len(glob.glob(read_path+'/G96*_S*.sext.gz'))
files_head = len(glob.glob(read_path+'/G96*_N*.arch_h'))+len(glob.glob(read_path+'/G96*_S*.arch_h')) 
logger.critical('There seems to be %i, %i images missing based on the SExtractor files, header files respectively.'%(files_sex-files_raw,files_head-files_raw))
files_trans = glob.glob(css.write_path()+css.red_path(date)+'/G96*Scorr.fits.fz')
candidates = []
for f in files_trans:
    with fits.open(f) as hdr:
        candidates.append(hdr[1].header['T-NTRANS'])
plt.hist(candidates,bins=np.arange(0, 9000, 100))
plt.title('Candidate summary for '+date)
plt.xlabel('Number of candidates per field')
plt.savefig(log_file_name+'.pdf')
my_file = {'file' : (log_file_name+'.pdf', open(log_file_name+'.pdf', 'rb'), 'pdf')}
with open('/home/saguaro/software/saguaro_slack.txt','r') as f:
    slack_token = f.readline().rstrip()
payload={"filename":log_file_name+'.pdf', "token":slack_token, "channels":['#pipeline']}
r = requests.post("https://slack.com/api/files.upload", params=payload, files=my_file)
logging.shutdown()
sys.exit()
