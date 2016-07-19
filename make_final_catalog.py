#! /usr/bin/env python
import argparse
from glob import glob
import os
import shutil
import re
import numpy as np
from astropy.io import fits
import pywcs
from astropy.table import Table
from astropy.io import ascii
from match_cats import match_cats
from run_SE import run_SE
from maglim import find_maglim
from trim_images import trim_image

def get_filter(image):
    hdr = fits.getheader(image)
    return hdr['FILTER'].rstrip("'")


def read_cat(catfile):
    '''Read in fits tables or ascii text files'''
    extension = os.path.splitext(catfile)[1]
    if extension == '.fits':
        f = fits.open(catfile)
        cat = f[1].data
        f.close()
    if extension == '.cat':
        cat = Table.read(catfile, format='ascii')
    return cat


def make_reg(RA, Dec, reg, radius, color, width=1, nsource=np.array([])):
    """Add circular regions in fk5 coords to a ds9 region file"""
    if nsource.shape[0] != 0:
        for r,d,nsrc in zip(RA, Dec, nsource):
            reg.write('circle(%.6f,%.6f,%f") # color=%s width=%i text={%i}\n'%\
                (r,d,radius,color,width,nsrc))
    else:
        for r,d in zip(RA, Dec):
            reg.write('circle(%.6f,%.6f,%f") # color=%s width=%i\n' % \
                (r,d,radius,color,width))



def fix_header(image, wispfield):
    '''Fix the DATASEC keyword in the image's header to extend
       out to the full size of the image. After IRAF's imcombine,
       the DATASEC keyword should have been updated to include 
       the new size of the image
       
       Add the total exposure time of the combined images to 
       the header
    '''
    im,hdr = fits.getdata(image, header=True)
    nx1 = hdr['NAXIS1']
    nx2 = hdr['NAXIS2']
    hdr['DATASEC'] = '[1:%i,1:%i]' % (nx1-1, nx2-1)

    if 'TEXPTIME' not in hdr:
        ims = hdr['IMCMB*']
        Nim = len(ims)
        texptime = 0.
        # add up all the exposure times
        for i in range(Nim):
            if os.path.splitext(ims[i])[1] != '.fits':
                ims[i] = ims[i] + '.fits'
            thdr = fits.getheader(os.path.join(wispfield,ims[i]))
            texptime += thdr['exptime']
        hdr.insert(19, ('TEXPTIME', texptime,
                               'Total EXPTIME of combined images'))
    fits.writeto(image, im, hdr, clobber=True)
    

def image_limits(image, wispfield):
    '''Find the '''
    # read in combined image
    hdr = fits.getheader(image)
    hdr_wcs = pywcs.WCS(hdr)
    # list of images that were combined
    ims = hdr['IMCMB*']
    Nim = len(ims)
    xx = np.zeros((Nim,2), dtype=int)
    yy = np.zeros((Nim,2), dtype=int)
    for i in range(Nim):
        if os.path.splitext(ims[i])[1] != '.fits':
            ims[i] = ims[i] + '.fits'
        hd = fits.getheader(os.path.join(wispfield,ims[i]))
        nx1 = hd['NAXIS1']
        nx2 = hd['NAXIS2']
        hd_wcs = pywcs.WCS(hd)
        # find WCS coordinates of input image corners
        coo = hd_wcs.wcs_pix2sky([[1,1],[1,nx1-1],[nx2-1,1],[nx2-1,nx1-1]],1)
        # find corresponding pixel coords in combined image
        pix = hdr_wcs.wcs_sky2pix(coo,1)
        xx[i,0] = pix[0][1]
        xx[i,1] = pix[1][1]
        yy[i,0] = pix[0][0]
        yy[i,1] = pix[2][0]

    xmin,xmax = np.max(xx[:,0]), np.min(xx[:,1])
    ymin,ymax = np.max(yy[:,0]), np.min(yy[:,1])
    return xmin,xmax,ymin,ymax


def setup(wispfield):
    '''Perform miscellaneous tasks to set up for the construction 
       of the final catalog. These include:
          - Copy any WISP IR catalogs 
          - Find how many Palomar filters are present
          - Determine which SExtractor mode to use
          - Run SExtractor
          - Read in calibration info
          - Collect all necessary info about each image in a dict
    '''
    # copy any WISP catalogs that exist to this directory
    ID = re.search('\d+', wispfield)
    parID = ID.group(0)
    WISPdir = '/data/highzgal/PUBLICACCESS/WISPS/data/Par'+parID+ \
              '/DATA/DIRECT_GRISM/'
    WISPfiles = [os.path.join(WISPdir,'fin_F110.cat'), \
                 os.path.join(WISPdir,'fin_F160.cat'), \
                 os.path.join(WISPdir,'fin_F140.cat') ]
    for wispf in WISPfiles:
        if os.path.exists(wispf):
            shutil.copy(wispf, wispfield)
    WISPimages = glob(os.path.join(WISPdir, 'F*W_drz.fits'))

    # get list of images for this field
    images = [x for x in [os.path.join(wispfield,'%s_g.fits'%wispfield), \
               os.path.join(wispfield,'%s_i.fits'%wispfield)] \
               if x in glob(os.path.join(wispfield,'%s_*.fits'%wispfield))]
    images.sort()

    # Run SE on all Palomar images, dual mode if both filters are present
    if len(images) == 1:
        SEmode = 'single'
    elif len(images) == 2:
        SEmode = 'dual'
    else:
        print "Don't know what to do with %i images" % len(images)

    print '\n--------------------------------------------'
    print 'Running %s in %s image mode'%([x for x in images],SEmode)
    print '--------------------------------------------\n'
    run_SE(images, 'Catalog', mode=SEmode)

    # read in calibration info for this field: 
    #   filter, color term, zero point, 1 sigma maglimit, sigma
    cal = np.genfromtxt(os.path.join(wispfield, 'sdss_calibration.dat'),
        dtype=[('filts','S10'), ('cterm',float), ('zp',float),
               ('maglim',float), ('sigma',float)])
    filts = cal['filts']
    cterm = cal['cterm']
    zp = cal['zp']
    maglim = cal['maglim']
    sigma = cal['sigma']

    # dict of info for each filter
    info = {}
    for image in images:
        filt = get_filter(image)
        # fix the DATASEC keyword in the headers, 
        # and add TEXPTIME if necessary
        fix_header(image, wispfield)
        # find limits of each image in pixel coordinates
        xmin,xmax,ymin,ymax = image_limits(image, wispfield)    
        # fill in info for image
        cat = os.path.splitext(image)[0] + '_final_cat.fits'
        info[filt] = {'cat':cat, 'cterm':cterm[filts == filt], \
                      'zp':zp[filts == filt], 'maglim':maglim[filts == filt], \
                      'sigma':sigma[filts == filt], \
                      'xmin':xmin, 'xmax':xmax, 'ymin':ymin, 'ymax':ymax}

    if WISPimages:
        info['WISPim'] = WISPimages[0]
    else:
        info['WISPim'] = []

    return info


def calibrate(mag, emag, zp, maglimit, cterm=None, color=None):
    '''Calibrate photometry and set magnitude limit'''
    # calibrate photometry
    if cterm:
        cal = mag + cterm*color + zp
    else:
        cal = mag + zp

    # set maglimit
    cal[cal >= maglimit] = maglimit
    emag[cal >= maglimit] = 0.0
    return cal,emag


def make_cat(info, WISPfilters, threshold, wispfield):
    '''Construct catalog by:
        1) Palomar photometry, AUTO magnitudes, run in dual image
           mode if both filters exist with i band for detection
        2) Calibrate photometry using parameters from sdss_calibration.dat
        3) Determining 1sigma limiting magnitude for each image
        4) Set all photometry fainter than 1sigma to 1sigma
    '''
    Palcats = [info[f]['cat'] for f in ['g','i'] if info.has_key(f)]
    Palcats.sort()
    filts = info.keys()
    filts = [x for x in filts if x != 'WISPim']

    # region file of Palomar objects
    reg = open(os.path.join(wispfield,'Palomar.reg'), 'w')
    reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
            'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 '
            'move=1 delete=1 include=1 source=1 \n')
    reg.write('fk5 \n')

    # read in Palomar catalogs
    pal1 = read_cat(Palcats[0])
    pal1x = pal1['X_IMAGE']
    pal1y = pal1['Y_IMAGE']
    # take only sources within full overlap
    f = 'i' if 'i' in filts else 'g'
    w = np.where((pal1x > info[f]['xmin']) & (pal1x < info[f]['xmax']) & 
                 (pal1y > info[f]['ymin']) & (pal1y < info[f]['ymax']))
    pal1x = pal1x[w]
    pal1y = pal1y[w]
    pal1RA = pal1['X_WORLD'][w]
    pal1Dec = pal1['Y_WORLD'][w]
    ID = pal1['NUMBER'][w]
    if len(Palcats) == 2:
        # 2 filters used
        # g in dual image mode
        AUTO_g = pal1['MAG_AUTO'][w]
        eAUTO_g = pal1['MAGERR_AUTO'][w]
        # i 
        pal2 = read_cat(Palcats[1])
        AUTO_i = pal2['MAG_AUTO'][w]
        eAUTO_i = pal2['MAGERR_AUTO'][w]

        # calibrate photometry
        color = AUTO_g - AUTO_i
        mag_g,emag_g = calibrate(AUTO_g, eAUTO_g, info['g']['zp'],
            info['g']['maglim'], cterm=info['g']['cterm'], color=color)
        mag_i,emag_i = calibrate(AUTO_i, eAUTO_i, info['i']['zp'],
            info['i']['maglim'], cterm=info['i']['cterm'], color=color)
        
    elif len(Palcats) == 1:
        if filts[0] == 'g':
            AUTO_g = pal1['MAG_AUTO'][w]
            eAUTO_g = pal1['MAGERR_AUTO'][w]
            mag_i = np.array([99.99]*pal1RA.shape[0])
            emag_i = np.array([-9.99]*pal1RA.shape[0])
            # calibrate photometry
            mag_g,emag_g = calibrate(AUTO_g, eAUTO_g, info['g']['zp'],
                info['g']['maglim'])

        if filts[0] == 'i':
            mag_g = np.array([99.99]*pal1RA.shape[0])
            emag_g = np.array([-9.99]*pal1RA.shape[0])
            AUTO_i = pal1['MAG_AUTO'][w]
            eAUTO_i = pal1['MAGERR_AUTO'][w]
            # calibrate photometry
            mag_i,emag_i = calibrate(AUTO_i, eAUTO_i, info['i']['zp'],
                info['i']['maglim'])

    # add all sources to region file
    make_reg(pal1RA, pal1Dec, reg, 2, 'blue', nsource=ID)
    reg.close()

    # construct Palomar catalog
    cat = Table([ID, pal1RA, pal1Dec, pal1x, pal1y, mag_g,emag_g,mag_i,emag_i],
                names=('PalID', 'PalRA', 'PalDec', 'Palxx', 'Palyy', 
                       'g', 'err_g', 'i', 'err_i'))
    cat['PalRA'].format = '{:.6f}'
    cat['PalDec'].format = '{:.6f}'
    cat['Palxx'].format = '{:.3f}'
    cat['Palyy'].format = '{:.3f}'
    cat['g'].format = '{:.3f}'
    cat['err_g'].format = '{:.3f}'
    cat['i'].format = '{:.3f}'
    cat['err_i'].format = '{:.3f}'
    cat.sort('PalID')
    ascii.write(cat, output=os.path.join(wispfield,'Palomar.cat'), 
                format='fixed_width_two_line', position_char='=')


    # if WISP filters are present:
    #   1) match Palomar and WISP catalogs
    #   2) create a combined catalog
    #   3) trim palomar images?
    if WISPfilters:
        # read in first WISP catalog for RA and Dec
        t1 = read_cat(WISPfilters[0])
        wispRA = t1['X_WORLD']
        wispDec = t1['Y_WORLD']
        ID = t1['NUMBER']
        # get photometry
        if os.path.join(wispfield,'fin_F110.cat') in WISPfilters:
            t110 = read_cat(os.path.join(wispfield,'fin_F110.cat'))
            mag_110 = t110['MAG_F1153W']
            emag_110 = t110['MAGERR_AUTO']
        else:
            mag_110 = np.array([99.99]*wispRA.shape[0])
            emag_110 = np.array([-9.99]*wispRA.shape[0])
        if os.path.join(wispfield,'fin_F140.cat') in WISPfilters:
            t140 = read_cat(os.path.join(wispfield,'fin_F140.cat'))
            mag_140 = t140['MAG_F1392W']
            emag_140 = t140['MAGERR_AUTO']
        else:
            mag_140 = np.array([99.99]*wispRA.shape[0])
            emag_140 = np.array([-9.99]*wispRA.shape[0])
        if os.path.join(wispfield,'fin_F160.cat') in WISPfilters:
            t160 = read_cat(os.path.join(wispfield,'fin_F160.cat'))
            mag_160 = t160['MAG_F1573W']
            emag_160 = t160['MAGERR_AUTO']
        else:
            mag_160 = np.array([99.99]*wispRA.shape[0])
            emag_160 = np.array([-9.99]*wispRA.shape[0])

        # match Palomar to WISP
        '''
        idx,separc = match_cats(wispRA, wispDec, pal1RA, pal1Dec)
        match = (separc.value*3600. <= threshold)
        nomatch = (separc.value*3600. > threshold)
        print '\n%i WISP objs matched to Palomar sources' % idx[match].shape[0]
        print '%i WISP objs NOT matched to Palomar sources\n' % \
            idx[nomatch].shape[0]
        '''
        idx,separc = match_cats(pal1RA, pal1Dec, wispRA, wispDec)
        match = (separc.value*3600. <= threshold)
        nomatch = (separc.value*3600. > threshold)
        print '\n%i Palomar objs matched to WISP sources' % idx[match].shape[0]

        # add successfully matched sources to region file with width=4
        reg = open(os.path.join(wispfield,'WISP.reg'), 'w')
        reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
                'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 '
                'move=1 delete=1 include=1 source=1 \n')
        reg.write('fk5 \n')
        make_reg(wispRA[idx[match]], wispDec[idx[match]], reg, 1.5, 
                 'green', width=4)
        make_reg(wispRA, wispDec, reg, 1.5, 'green')
        reg.close()

        # create arrays for info from Palomar catalog
        palRA = np.array([99.99]*wispRA.shape[0])
        palDec = np.array([99.99]*wispRA.shape[0])
        pal_g = np.array([99.99]*wispRA.shape[0])
        epal_g = np.array([-9.99]*wispRA.shape[0])
        pal_i = np.array([99.99]*wispRA.shape[0])
        epal_i = np.array([-9.99]*wispRA.shape[0])
        wisparrays = [palRA, palDec, pal_g, epal_g, pal_i, epal_i] 
        palarrays = [pal1RA, pal1Dec, mag_g, emag_g, mag_i, emag_i]   

        for warr,parr in zip(wisparrays,palarrays):
            warr[idx[match]] = parr[match]

        wispCat = Table([ID, wispRA, wispDec, palRA, palDec, pal_g, epal_g,
                         pal_i, epal_i, mag_110, emag_110, mag_140, emag_140,
                         mag_160, emag_160],
                        names=('ID','RA_WISP', 'Dec_WISP', 'RA_Pal', 'Dec_Pal',
                               'g', 'eg', 'i', 'ei', 'F110W', 'eF110W',
                               'F140W', 'eF140W', 'F160W', 'eF160W'))
        wispCat['RA_Pal'].format = '{:.6f}'
        wispCat['Dec_Pal'].format = '{:.6f}'
        wispCat['g'].format = '{:.3f}'
        wispCat['eg'].format = '{:.3f}'
        wispCat['i'].format = '{:.3f}'
        wispCat['ei'].format = '{:.3f}'
        wispCat.sort('ID')
        ascii.write(wispCat,output=os.path.join(wispfield,'Palomar_WISPS.cat'),
                    format='fixed_width_two_line', position_char='=')

        # trim images
        WISPimage = info['WISPim']
        WISPim,WISPhdr = fits.getdata(WISPimage, header=True)
        cRA = WISPhdr['CRVAL1'] 
        cDec = WISPhdr['CRVAL2']
        xsize,ysize = 3.,3.
        for filt in filts:
            Palim = os.path.join(wispfield, '%s_%s.fits'%(wispfield,filt))
            outname = os.path.splitext(Palim)[0] + '_WISPFOV.fits'
            trim_image(Palim, cRA, cDec, xsize, ysize, outname)


def main():
    parser = argparse.ArgumentParser(description=
        'Run SE on Palomar images, calibrate photometry, and create catalog.')
    parser.add_argument('wispfield', type=str, nargs=1,
        help='WISP field for which to construct a catalog. ' +\
             'Must match name of directory containing all relevant files')
    parser.add_argument('--threshold', type=float, default=0.75,
        help='Radius in arcsec to be used for matching Palomar and WISP '+\
             'catalogs. Default is 0.75"')
    args = parser.parse_args()
    wispfield = args.wispfield[0].strip('/')
    threshold = args.threshold

    # setup for creation of catalog
    # info contains colorterms, zp's, maglimits, the background sigma
    # and the name of the SE catalog for each filter
    info = setup(wispfield)
    
    # construct catalog
    WISPfilters = glob(os.path.join(wispfield,'fin_*.cat'))
    WISPfilters.sort()
    make_cat(info, WISPfilters, threshold, wispfield)



if __name__ == '__main__':
    main()


