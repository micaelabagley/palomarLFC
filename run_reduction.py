#! /usr/bin/env python
########################################################################
## Micaela Bagley
## November 2014
## 
## usage: reduction.py [-h] [--biasdir BIASDIR] [--darkdir DARKDIR]
##                     [--flatdir FLATDIR] [--SaveSteps] 
##                     [--LogSteps LOGSTEPS] [--DarkSubtract] 
##                     [--UseBias USEBIAS] [--UseDark USEDARK]
##                     [--UseFlat [USEFLAT [USEFLAT ...]]]
##                     directory binning
## 
## Reduce a night of Palomar 200" data
##
## positional arguments:
##   directory             Directory to be reduced. Ex: 23mar12
##   binning               Binning of data to be reduced. ('binned' or
##                        'unbinned')
##
## optional arguments:
##   -h, --help            show this help message and exit
##   --biasdir BIASDIR     Directory containing the biases. Default is biases/.
##   --darkdir DARKDIR     Directory containing the darks. Default is darks/.
##   --flatdir FLATDIR     Directory containing the flats. Default is flats/.
##   --SaveSteps           Save each reduction step as a new image, with 
##                         strings ('bs', 'ds', 'ff') appended to 
##                         filenames to indicate reduction steps taken.
##   --LogSteps LOGSTEPS   Create a log file to record the creation of
##                         calibration files and the science images 
##                         included in reduction.
##   --DarkSubtract        Dark subtract the data.
##   --UseBias USEBIAS     Use a previously-created Master Bias.
##   --UseDark USEDARK     Use a previously-created Master Dark.
##   --UseFlat [USEFLAT [USEFLAT ...]]
##                         Use a previously-created Master Flat. List a 
##                         Flat for as many filters as necessary.
########################################################################
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table
from glob import glob
import time
import os
from reduction import biasproc
from reduction import darkproc
from reduction import flatproc
from reduction import create_log


def get_filter(image):
    hdr = fits.getheader(image)
    return hdr['FILTER'].rstrip("'")


def check_header(image, keyword):
    '''Check header for presence of keyword'''
    if keyword in fits.getheader(image):
        return True
    else:
        return False


def expt_norm(image):
    '''Normalize an image by its exposure time'''
    im,hdr = fits.getdata(image, header=True)
    exptime = hdr['EXPTIME']
    new = im / exptime
    # add a header keyword
    hdr['EXPTNORM'] = (exptime, 'Image normalized by EXPTIME')
    fits.writeto(image, new, header=hdr, clobber=True)


def main():
    parser = argparse.ArgumentParser(description=
        'Reduce a night of Palomar 200" data')
    parser.add_argument('directory', type=str, nargs=1, 
        help='Directory to be reduced. Ex: 23mar12 ')
    parser.add_argument('binning', type=str, nargs=1, 
        help="Binning of data to be reduced. (binned or unbinned)")
    parser.add_argument('--biasdir', type=str, default='biases',
        help='Directory containing the biases. Default is biases/.')
    parser.add_argument('--darkdir', type=str, default='darks',
        help='Directory containing the darks. Default is darks/.')
    parser.add_argument('--flatdir', type=str, default='flats',
        help='Directory containing the flats. Default is flats/.')
    parser.add_argument('--Overwrite', action='store_true',
        help="Overwrite files after each reduction step. If not set, " +\
             "images are saved as separate files, with strings " +\
             "('bs', 'ds', 'ff') appended to filenames to indicate " +\
             "reduction steps taken.")
    parser.add_argument('--LogSteps', type=str,
        help='Create a log file to record the creation of calibration ' +\
             'files and the science images included in reduction.')
    parser.add_argument('--DarkSubtract', action='store_true',
        help='Dark subtract the data.')
    parser.add_argument('--UseBias', type=str,
        help='Use a previously-created Master Bias.')
    parser.add_argument('--UseDark', type=str,
        help='Use a previously-created Master Dark.')
    parser.add_argument('--UseFlat', type=str, nargs='*',
        help='Use a previously-created Master Flat. List a Flat for ' +\
             'as many filters as necessary.')
    args = parser.parse_args()
    directory = args.directory[0]
    binning = args.binning[0]
    biasdir = args.biasdir
    darkdir = args.darkdir
    flatdir = args.flatdir
    Overwrite = args.Overwrite
    SaveSteps = False if Overwrite else True
    LogSteps = args.LogSteps
    DarkSubtract = args.DarkSubtract
    UseBias = args.UseBias
    UseDark = args.UseDark
    UseFlat = args.UseFlat

    # does user wish to create a log file for the calibration files?
    if LogSteps :
        # initiate log file
        logfile = LogSteps
        create_log.logfile_init(logfile, directory, binning)

    # set up binning
    if binning == 'binned':
        bin1,bin2 = (2,2)
    if binning == 'unbinned':
        bin1,bin2 = (1,1)

    # filenames - extensions for SaveSteps option
    bstr = 'bs'
    dstr = 'ds'
    fstr = 'ff'

    ############################
    '''GET LISTS OF ALL FILES'''
    ############################
#    if not UseBias:
    # list of all biases in biasdir
    biaslist_all = glob(os.path.join(biasdir,'ccd*0.fits'))
    # find all the biases with the proper binning
    biaslist = []
    for bias in biaslist_all:
        hdr = fits.getheader(bias)
        if (hdr['CCDBIN1'],hdr['CCDBIN2']) == (bin1,bin2):
            biaslist.append(bias)

#    if not UseFlat:
    # list of all flats in flatdir
    # base files
    flatlist_all = glob(os.path.join(flatdir,'ccd*0.fits'))
    # find all the flats with proper binning
    flatlist = []
    for flat in flatlist_all:
        hdr = fits.getheader(flat)
        if (hdr['CCDBIN1'],hdr['CCDBIN2']) == (bin1,bin2):
            flatlist.append(flat)
   
#    if DarkSubtract:
    if not UseDark:
        # list of all darks in darkdir
        darklist_all = glob(os.path.join(darkdir,'ccd*0.fits'))
        # find all the darks with proper binning
        darklist = []
        for dark in darklist_all:
            hdr = fits.getheader(dark)
            if (hdr['CCDBIN1'],hdr['CCDBIN2']) == (bin1,bin2):
                darklist.append(dark)
   
    # list of all science images in scidir
    scilist_all = glob(os.path.join(directory,'ccd*0.fits'))
    # find all the science images with the proper binning
    scilist = []
    for sci in scilist_all:
        hdr = fits.getheader(sci)
        if (hdr['CCDBIN1'],hdr['CCDBIN2']) == (bin1,bin2):
            scilist.append(sci)

    ##########################
    ''' COMBINE THE BIASES '''
    ##########################
    # is a previously-created Master Bias to be used?
    if UseBias:
        MasterBias = UseBias
        # add to log file?
        if LogSteps:
            create_log.logfile_cals(logfile, 'bias', [], binning, 
                UseCal=MasterBias)
    
    else:
        print '\nCombining the biases'
        MasterBias = biasproc.bias_combine(biaslist, binning)
        # add to log file?
        if LogSteps:
            create_log.logfile_cals(logfile, 'bias', biaslist, binning)
    
    ########################
    ''' BIAS SUBTRACTION '''
    ########################
    # make list of all images that need to be bias subtracted
    biasSub = []
    bs_extension = bstr + '.fits'
    if DarkSubtract:
        if not UseDark:
            # bias subtract the darks
            for dark in darklist:
                bs_dark = dark.replace('fits',bs_extension)
                check_biassub = check_header(dark, 'BIASSUB')
                if (os.path.exists(bs_dark) == 0) & (check_biassub == 0):
                    biasSub.append(dark)
    if not UseFlat:
        # bias subtract the flats
        for flat in flatlist:
            bs_flat = flat.replace('fits',bs_extension)
            check_biassub = check_header(flat, 'BIASSUB')
            if (os.path.exists(bs_flat) == 0) & (check_biassub == 0):
                biasSub.append(flat)
    # bias subtract the science images
    for sci in scilist:
        bs_sci = sci.replace('fits',bs_extension)
        check_biassub = check_header(sci, 'BIASSUB')
        if (os.path.exists(bs_sci) == 0) & (check_biassub == 0):
            biasSub.append(sci)

    if not biasSub:
        print '\nAll images already bias-subtracted. Skipping.'
    else:
        print '\nBias subtracting images'
        biasproc.bias_subtract(biasSub, MasterBias, SaveSteps=SaveSteps)

    if SaveSteps:
        flatlist = [x.replace('fits', bs_extension) for x in flatlist]
        scilist = [x.replace('fits', bs_extension) for x in scilist]
        if DarkSubtract:
            darklist = [x.replace('fits', bs_extension) for x in darklist]

    if DarkSubtract:
        ##########################
        ''' COMBINE THE DARKS '''
        ##########################
        # is a previously-created Master Dark to be used?
        if UseDark:
            MasterDark = UseDark
            # add to log file?
            if LogSteps:
                create_log.logfile_cals(logfile, 'dark', [], binning, 
                    UseCal=MasterDark)
    
        else:
            print '\nCombining the darks'
            MasterDark = darkproc.dark_combine(darklist, binning)
            # add to log file?
            if LogSteps:
                create_log.logfile_cals(logfile, 'dark', darklist, binning)
    
        ########################
        ''' DARK SUBTRACTION '''
        ########################
        # make list of all images that need to be dark subtracted
        darkSub = []
        ds_extension = dstr + '.fits'
        if not UseFlat:
            # dark subtract the flats
            for flat in flatlist:
                ds_flat = flat.replace('fits',ds_extension)
                check_darksub = check_header(flat, 'DARKSUB')
                if (os.path.exists(ds_flat) == 0) & (check_darksub == 0):
                    darkSub.append(flat)

        # dark subtract the science images
        for sci in scilist:
            ds_sci = sci.replace('fits',ds_extension)
            check_darksub = check_header(sci, 'DARKSUB')
            if (os.path.exists(ds_sci) == 0) & (check_darksub == 0):
                darkSub.append(sci)

        if not darkSub:
            print '\nAll images already dark-subtracted. Skipping.'
        else:
            print '\nDark subtracting images'
            darkproc.dark_subtract(darkSub, MasterDark, SaveSteps=SaveSteps)
    
        if SaveSteps:
            flatlist = [x.replace('fits', ds_extension) for x in flatlist]
            scilist = [x.replace('fits', ds_extension) for x in scilist]
    
    ###########################################################
    ''' COMBINE THE FLATS AND FLAT FIELD THE SCIENCE IMAGES '''
    ###########################################################
    filters = [get_filter(x) for x in scilist]
    filts = np.unique(filters)
    for f in filts:
        if UseFlat:
            MasterFlat = [x for x in UseFlat if get_filter(x) == f]
            if MasterFlat: 
                MasterFlat = MasterFlat[0]
            # add to log file?
            if LogSteps:
                create_log.logfile_cals(logfile, 'flat', [], binning, filt=f, 
                    UseCal=MasterFlat)
            
        else:
            print '\nCombining the %s band flats' % f
            # get list of flats for this filter
            flatlist_filt = [x for x in flatlist if get_filter(x) == f]
            MasterFlat = flatproc.flat_combine(flatlist_filt, f, binning)
            # add to log file?
            if LogSteps:
                create_log.logfile_cals(logfile, 'flat', flatlist_filt, 
                    binning, filt=f)
    
        #####################
        ''' FLAT FIELDING '''
        #####################
        ff_extension = fstr + '.fits'
        flatField = []
        filtlist = [x for x in scilist if get_filter(x) == f]
        for sci in filtlist:
            ff_sci = sci.replace('fits',ff_extension)
            check_flatfield = check_header(sci, 'FLATFLD')
            if (os.path.exists(ff_sci) == 0) & (check_flatfield == 0):
                flatField.append(sci)
        
        if not flatField:
            print '\nScience images already flat-fielded. Skipping.'
        else:
            print '\nFlat fielding the science images'
            flatproc.flat_field(flatField, MasterFlat, SaveSteps=SaveSteps)

    if SaveSteps:
        scilist = [x.replace('fits', ff_extension) for x in scilist]
    
    # add science images to log file?
    if LogSteps:
        create_log.logfile_sci(logfile, scilist, SaveSteps=SaveSteps)

    # normalize images by their exptimes
    # check that science images have not already been exptime normalized
    exptNorm = []
    for sci in scilist:
        if 'EXPTNORM' not in fits.getheader(sci):
            exptNorm.append(sci)
    
    if not exptNorm:
        print '\nScience images already exptime-normalized. Skipping.'
    else:
        print '\nNormalizing images by their exposure times.'
        for im in exptNorm:
            expt_norm(im)
    
    
if __name__ == '__main__':
    main()
