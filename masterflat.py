# astropy
from astropy.io import fits
from astropy import units as u
from astropy.nddata import CCDData
from astropy.modeling import models
# others
import glob
import ccdproc
import numpy as np

#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#

flat_list_g = glob.glob('allframe/flatg/frame*.fits', recursive = True)
flat_list_u = glob.glob('allframe/flatu/frame*.fits', recursive = True)
flat_list_r = glob.glob('allframe/flatr/frame*.fits', recursive = True)

# array-flipping functions, conveniently named
def flipy(img):
    return img[::-1]
def flipx(img):
    return np.flip(img,1)
def flipxy(img):
    return np.flip(img,1)[::-1]

def arrange_flat_quad(red, header, quadnum):
    '''flip and and fill in the 4 reduced arrays, based on quadnum
        - wastes time by filling in the full size empty array
        - this function can definitely be improved!!!'''

    # quadnum condition for all 16 amps
    if quadnum == 1 or quadnum == 2:
        red[0].data = flipxy(red[0].data.T)
        red[1].data = flipx(red[1].data.T)
        red[2].data = flipy(red[2].data.T)
        red[3].data = red[3].data.T
    if quadnum == 3 or quadnum ==4:
        red[0].data = red[0].data.T
        red[1].data = flipy(red[1].data.T)
        red[2].data = flipx(red[2].data.T)
        red[3].data = flipxy(red[3].data.T)

    # build empty array to fit all 4 quads        
    emptysize = header['DETSIZE'][1:-1].split(',')
    empty = np.zeros(shape=(int(emptysize[0][2:]),int(emptysize[1][2:])))

    # indeces to map reduced arrays to correct place in empty array
    lowys = [int(min(red[i].header['DETSEC'][1:-1].split(',')[0].split(':')))
             for i in range(len(red))]
    lowxs = [int(min(red[i].header['DETSEC'][1:-1].split(',')[1].split(':')))
             for i in range(len(red))]
    
    # assign each extension's data to the correct section of empty array
    for i in range(len(red)):
        highy = lowys[i] + int(red[0].header['DATASEC'][1:-1].split(',')[0][2:])-1
        highx = lowxs[i] + int(red[0].header['DATASEC'][1:-1].split(',')[1][2:])-1
        empty[lowys[i]-1:highy, lowxs[i]-1:highx] = red[i].data
    
    # split the large empty array into 4 arrays (one for each CCD, "quad data")
    ycut, xcut = int(int(emptysize[0][2:])/2), int(int(emptysize[1][2:])/2)
    qdata = [empty[:ycut,:xcut].T, empty[:ycut,xcut:].T, 
             empty[ycut:,:xcut].T, empty[ycut:,xcut:].T] #this fills only one
    
    # select the non-empty array from qdata, set equal to flat_full
    nonzero_array = qdata[np.argmax([np.sum(qdata[i]) for i in range(len(qdata))])]
    flat_full = CCDData((nonzero_array), header=header, unit=u.electron/u.second)

    return flat_full

def create_flat_quad(flat_list, quadnum, red_path, mbias=None):
    '''use ccdproc to reduce each flat, combine, and writeto masterflat file
        ** takes a little while to run, only works on one quadrant at a time **
       '''
    scale = []
    flats = []
    extnames = [int(fits.open(flat_list[0])[i].header['EXTNAME'][2:]) for i in range(1,17)]
    gains = [float(fits.open(flat_list[0])[0].header['GAIN'+str(extnames[i])]) for i in range(16)]
    readnoises = [float(fits.open(flat_list[0])[0].header['RDNOIS'+str(extnames[i])]) for i in range(16)]
    mbias = [CCDData.read('allframe/biases/mbias_16.fits', hdu=x+1, unit=u.electron) for x in range(16)]
    #mbias = [CCDData.read('allframe/biases/mbias_16.T.fits', hdu=x+1, unit=u.electron) for x in range(16)]
    pick_amps ={1:(0,4), 2:(4,8), 3:(8,12), 4:(12,16)}
    low,high = pick_amps[quadnum][0],pick_amps[quadnum][1]
        
    for i, flat in enumerate(flat_list):
        #print('processing flat no. ', i, '...')
        with fits.open(flat) as hdr:
            header = hdr[0].header

        raw = [CCDData.read(flat, hdu=x+1) for x in range(low,high)]
        red = [ccdproc.ccd_process(x, oscan=x.header['BIASSEC'], 
                                   oscan_model=models.Chebyshev1D(3), 
                                   trim=x.header['DATASEC'], 
                                   gain=gains[low+k]*u.electron/u.adu, 
                                   readnoise=readnoises[low+k]*u.electron, 
                                   master_bias=mbias[low+k],
                                   gain_corrected=True) for k,x in enumerate(raw)]
        flat_full = arrange_flat_quad(red, header, quadnum)

    
        norm = 1/np.median(flat_full[200:800,1200:1800])
        scale.append(norm)
        flats.append(flat_full)

    print('creating master flat for quadrant ' + str(quadnum) + '...')
    mflat = ccdproc.combine(flats,method='average',scale=scale,sigma_clip=True)
    # something weird happened with naming, this temporary condition fixes the issue
    if quadnum == 2:
        mflat.write(red_path+'mflat3_90prime.fits',overwrite=True)
    if quadnum == 3:
        mflat.write(red_path+'mflat2_90prime.fits',overwrite=True)
    else:
        mflat.write(red_path+'mflat'+str(quadnum)+'_90prime.fits',overwrite=True)
    print('master flat created for quadrant ' + str(quadnum))
    return

def create_mflats(filt, red_path):
    '''make 4 master flats (one for each ccd) in a given filter'''
    if filt == 'g':
        flat_list = flat_list_g
    if filt == 'u':
        flat_list = flat_list_u
    if filt == 'r':
        flat_list = flat_list_r

    for i in range(1,5):
        create_flat_quad(flat_list, i, red_path)

