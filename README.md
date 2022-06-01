# SAGUARO Pipeline for 90prime Data

SAGUARO’s main pipeline handles data reduction and image subtraction for incoming telescope data. The scripts above were created to work with the newly incorporated 2.3m Bok Telescope. 

## Setup

To run these scripts, you first need the “master bias” and “master flat” files on your machine. You have two options:

1. Download the master files which were created using **masterbias.py** and **masterflat.py,** or,
2. Download individual flat and bias frames taken by 90prime, then run these through **masterbias.py** and **masterflat.py** yourself (see Usage steps 1 and 2).

In either case, the master files should be saved to a folder named ‘allframe’ in the same location as these scripts. Contact saarah@sas.upenn.edu for a link to the Box drive containing all relevant files (as they are too large to upload to github).

## Usage

1. make a master bias (one satisfies all filters) using **masterbias.py**
    
    ```python
    import masterbias as mb
    mb.create_bias('allframe/biases/')
    ```
    
2. make master flats (filter dependent) using **masterflat.py**
    
    ```python
    import masterflat as mf
    mf.create_mflats('g', 'allframe/flatg/') # do this for g, r, and u
    ```
    
3. reduce your data from 90prime using **bok90prime.py**
    
    ```python
    import bok90prime as bok
    science_file = 'data/d0424.0064.fits' # raw image taken by 90prime
    bok.science_process(science_file)
    # makes 4 reduced files (one for each quadrant):
    # 'data/d0424.0064_q#red.fits' for # 1-4
    ```