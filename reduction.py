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
import pyfits
from astropy.table import Table
from glob import glob
import time
import os
import biasproc
import darkproc
import flatproc
import create_log


def get_filter(image):
    hdr = pyfits.getheader(image)
    return hdr['FILTER'].rstrip("'")


def expt_norm(image):
    '''Normalize an image by its exposure time'''
    im,hdr = pyfits.getdata(image, header=True)
    exptime = hdr['EXPTIME']
    new = im / exptime
    # add a header keyword
    hdr['EXPTNORM'] = (exptime, 'Image normalized by EXPTIME')
    pyfits.writeto(image, new, header=hdr, clobber=True)


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
    parser.add_argument('--SaveSteps', action='store_true',
        help="Save each reduction step as a new image, with strings " +\
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
    SaveSteps = args.SaveSteps
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

    ############################
    '''GET LISTS OF ALL FILES'''
    ############################
    if not UseBias:
        # list of all biases in biasdir
        biaslist_all = glob(os.path.join(biasdir,'*.fits'))
        # find all the biases with the proper binning
        biaslist = []
        for bias in biaslist_all:
            hdr = pyfits.getheader(bias)
            if (hdr['CCDBIN1'],hdr['CCDBIN2']) == (bin1,bin2):
                biaslist.append(bias)

    if not UseFlat:
        # list of all flats in flatdir
        flatlist_all = glob(os.path.join(flatdir,'*.fits'))
        # find all the flats with proper binning
        flatlist = []
        for flat in flatlist_all:
            hdr = pyfits.getheader(flat)
            if (hdr['CCDBIN1'],hdr['CCDBIN2']) == (bin1,bin2):
                flatlist.append(flat)

    if DarkSubtract:
        if not UseDarks:
            # list of all darks in darkdir
            darklist_all = glob(os.path.join(darkdir,'*.fits'))
            # find all the darks with proper binning
            darklist = []
            for dark in darklist_all:
                hdr = pyfits.getheader(dark)
                if (hdr['CCDBIN1'],hdr['CCDBIN2']) == (bin1,bin2):
                    darklist.append(dark)

    # list of all science images in scidir
    scilist_all = glob(os.path.join(directory,'*.fits'))
    # find all the science images with the proper binning
    scilist = []
    for sci in scilist_all:
        hdr = pyfits.getheader(sci)
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
        print 'Combining the biases'
        MasterBias = biasproc.bias_combine(biaslist, binning)
        #MasterBias = 'Bias_%s.fits'%binning
        # add to log file?
        if LogSteps:
            create_log.logfile_cals(logfile, 'bias', biaslist, binning)
    
    ########################
    ''' BIAS SUBTRACTION '''
    ########################
    # make list of all images that need to be bias subtracted
    biasSub = []
    # bias subtract the darks
    if DarkSubtract:
        if not UseDark:
            for dark in darklist:
                if 'BIASSUB' not in pyfits.getheader(dark):
                    biasSub.append(dark)
    if not UseFlat:
        # bias subtract the flats
        for flat in flatlist:
            if 'BIASSUB' not in pyfits.getheader(flat):
                biasSub.append(flat)
    # bias subtract the science images
    for sci in scilist:
        if 'BIASSUB' not in pyfits.getheader(sci):
            biasSub.append(sci)
    
    if not biasSub:
        print 'All images already bias-subtracted. Skipping.'
    else:
        print 'Bias subtracting images'
        biasproc.bias_subtract(biasSub, MasterBias, SaveSteps=SaveSteps)
    
    
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
            print 'Combining the darks'
            if SaveSteps:
                darklist_save = [''.join([os.path.splitext(x)[0], '.bs.fits']) \
                    for x in darklist]
                MasterDark = darkproc.dark_combine(darklist_save, binning)
            else:
                MasterDark = darkproc.dark_combine(darklist, binning)
            #MasterDark = 'Dark_%s.fits'%binning
            # add to log file?
            if LogSteps:
                create_log.logfile_cals(logfile, 'dark', darklist, binning)
    
        ########################
        ''' DARK SUBTRACTION '''
        ########################
        # make list of all images that need to be dark subtracted
        darkSub = []
        if not UseFlat:
            # dark subtract the flats
            for flat in flatlist:
                if 'DARKSUB' not in pyfits.getheader(flat):
                    darkSub.append(flat)
        # dark subtract the science images
        for sci in scilist:
            if 'DARKSUB' not in pyfits.getheader(sci):
                darkSub.append(sci)
    
        if not darkSub:
            print 'All images already dark-subtracted. Skipping.'
        else:
            print 'Dark subtracting images'
            darkproc.dark_subtract(darkSub, MasterDark, SaveSteps=SaveSteps)
    
    
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
            print 'Combining the %s band flats' % f
            # get list of flats for this filter
            flatlist_filt = [x for x in flatlist if get_filter(x) == f]
####### NEEDS WORK:
#######     creating lists of files to either combine or process
#######     when SaveingSteps. Don't want to flat-field all *.fits,
#######     just the ones that are *bs.fits or *ds.fits
#######
            if SaveSteps:
                flatlist_save = [x for x in \
                                 glob(os.path.join(flatdir,'*.bs.fits')) \
                                 if x in flatlist_filt]
                MasterFlat = flatproc.flat_combine(flatlist_save, f, binning)
            else:
                MasterFlat = flatproc.flat_combine(flatlist_filt, f, binning)
            #MasterFlat = 'Flat%s_%s.fits' % (f, binning)
            # add to log file?
            if LogSteps:
                create_log.logfile_cals(logfile, 'flat', flatlist_filt, 
                    binning, filt=f)
    
        # check that science images have not already been flat fielded
        flatField = []
        filtlist = [x for x in scilist if get_filter(x) == f]
        for sci in filtlist:
            if 'FLATFLD' not in pyfits.getheader(sci):
                flatField.append(sci)
        if not flatField:
            print 'Science images already flat-fielded. Skipping.'
        else:
            print 'Flat fielding the science images'
            flatproc.flat_field(flatField, MasterFlat, SaveSteps=SaveSteps)
        
    
    # add science images to log file?
    if LogSteps:
        create_log.logfile_sci(logfile, scilist, SaveSteps=SaveSteps)

    # normalize images by their exptimes
    # check that science images have not already been exptime normalized
    exptNorm = []
    for sci in scilist:
        if 'EXPTNORM' not in pyfits.getheader(sci):
            exptNorm.append(sci)
    if not exptNorm:
        print 'Science images already exptime-normalized. Skipping.'
    else:
        print 'Normalizing images by their exposure times.'
        for im in exptNorm:
            expt_norm(im)
    
    
if __name__ == '__main__':
    main()
