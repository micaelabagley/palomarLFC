#! /usr/bin/env python
import time
from astropy.io import fits
import os
from glob import glob

def get_filter(image):
    hdr = fits.getheader(image)
    return hdr['FILTER'].rstrip("'")


def logfile_init(logfile, night, binning):
    '''Start a reduction log

       logfile  - Log file name
       night    - String identifying the night of data to be reduced
       binning  - String identifying binning of images: 'binned' or 'unbinned
    '''
    log = open(logfile, 'w')
    logstr = '###    Reduction log for %s, %s data   ### \n' % (night, binning)
    log.write('#' * len(logstr) + '\n')
    log.write(logstr)
    log.write('#' * len(logstr) + '\n')

    # get the date and time
    now = time.strftime('%c')
    log.write('# ' + now + '\n')
    log.write('\n')
    log.close()


def logfile_cals(logfile, caltype, imlist, binning, filt=None, UseCal=None):
    '''Add to reduction log info on calibration files.
       
       logfile  - Log file name
       caltype  - String identifying type of calibration file: 
                    'bias', 'dark', 'flat'. 
       imlist   - List of files combined to create master frame
       binning  - String identifying binning of images: 'binned' or 'unbinned'
       filt     - If caltype='flat', provide the filter as well (string)
       UseCal   - If using a previously-created calibration file, provide
                  the name of the file
    '''
    nfiles = len(imlist)
    log = open(logfile, 'a')


    # basic info about the calibration file
    if caltype == 'flat':
        # get filename of master calibration file created
        calfile = '%s%s%s_%s.fits' % (caltype[0].upper(), caltype[1:], 
            filt, binning)
        
        # add filter information to logfile
        log.write('# Master %s band %s%s \n' % 
            (filt, caltype[0].upper(), caltype[1:]))
    else:
        calfile = '%s%s_%s.fits' % (caltype[0].upper(), caltype[1:], binning)
        log.write('# Master %s%s \n' % (caltype[0].upper(), caltype[1:]))

    # using a previously-made calibration file?    
    if UseCal is not None:
        log.write('Using previously made file %s\n' % (UseCal))
        log.write('\n')
        log.close()

    else:
        log.write('%s, NCOMBINE=%i \n' % (calfile, nfiles))

        # list the files that were combined to form master frame
        for im in sorted(imlist):
            log.write('    %s \n' % os.path.basename(im))
        log.write('\n')

        # plural form of caltype
        if caltype[-1] == 's':
            plural = '%ses'%caltype
        else:
            plural = '%ss'%caltype

        # add unused files to logfile
        log.write('    Unused %s: \n'%plural)
        caldir = os.path.dirname(imlist[0])
        if caltype == 'flat':
            # add only unused flats of the same filter
            unused = [x for x in glob(os.path.join(caldir,'unused','*.fits')) \
                if get_filter(x) == filt]
        else:
            unused = glob(os.path.join(caldir, 'unused', '*.fits'))
        # desired binning of images
        if binning == 'binned':
            bin1,bin2 = (2,2)
        elif binning == 'unbinned':
            bin1,bin2 = (1,1)
        for im in sorted(unused):
            # check that unused image has the correct binning
            hdr = fits.getheader(im)
            if (hdr['CCDBIN1'],hdr['CCDBIN2']) == (bin1,bin2):
                log.write('    %s  - \n' % (os.path.basename(im)))
        log.write('\n')
        log.write('\n')
        log.close()


def logfile_sci(logfile, imlist, SaveSteps=False):
    '''Add science images and info to logfile
        
       logfile  - Log file name
       imlist   - List of science files reduced 
     SaveSteps  - Set to True if each step in the reduction is saved
                  as a new file
    '''
    imlist.sort()
    nfiles = len(imlist)
    log = open(logfile, 'a')
    log.write('# Science Images\n')
    if SaveSteps:   
        # if user is saving interim reduction steps to new files, 
        # include the extensions for each caliration type to header
        log.write('#           *.fits - raw image \n')
        log.write('#        *.bs.fits - bias-subtracted image \n')
        log.write('#     *.bs.ff.fits - bias-subtracted and flat-fielded \n')
    log.write('# Image  Object  Exptime  Filter\n')
    # for each file, print the image name, object type, exptime, and filter
    for i in range(nfiles):
        imname = os.path.splitext(os.path.basename(imlist[i]))[0]
        obj = fits.getheader(imlist[i])['OBJECT']
        exptime = fits.getheader(imlist[i])['EXPTIME']
        filt = fits.getheader(imlist[i])['FILTER']
        log.write('%s    %s   %.3f   %s \n' % (imname, obj, exptime, filt))
    log.write('\n')
    log.write('\n')
    log.write('# Comments: \n')
    log.write('\n')
    log.close()



def main():
    pass


if __name__ == '__main__':
    main()
