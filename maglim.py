#! /usr/bin/env python
import argparse
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit
from PyPhot import aper
import os

class CircularAperture():
    def __init__(self, im, seg, xc, yc, radius):
        self.im = im
        self.seg = seg
        self.xc = xc
        self.yc = yc
        self.radius = radius
        self.clean = self.circular_aperture()
#        self.aperture = self.circular_aperture()
#        self.clean = self.check_seg()

    def circular_aperture(self):
        '''Create a circular mask'''
        clean = np.zeros(self.xc.shape, dtype=bool)
        for i,v in enumerate(self.xc):
            # get x and y vectors of image
            y,x = np.ogrid[:self.seg.shape[0], :self.seg.shape[1]]
            # distance of every pixel from the center
            r2 = (x - self.xc[i])**2 + (y - self.yc[i])**2
            # find all pixels inside a circular aperture
            mask = r2 <= (self.radius)**2
            clean[i] = self.check_seg(mask)
            # maks = r2 <= (radius - 0.5)**2
        return clean

    def check_seg(self, mask):
        '''Check that all pixels in aperture are 0. 
           Non-zero pixels indicate contamination from source flux'''
        segcheck = self.seg[mask]
        imgcheck = self.im[mask]
        if (np.any(segcheck) != 0) | (np.any(imgcheck) == 0):
            # aperture is contaminated 
            return False
        else:
            return True


def make_reg(xc, yc, reg, radius):
    '''Add circular regions in image coords to a ds9 region file'''
    for x,y in zip(xc,yc):
        reg.write('circle(%.3f,%3.f,%f)\n' % (x,y,radius))


def photometry(image, segmap, naper, rad, zp):
    '''Perform photometry'''
    reg = open('test.reg', 'w')
    reg.write('global color=green dashlist=8 3 width=1 font="helvetica 10 '
              'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 '
              'move=1 delete=1 include=1 source=1 \n')
    reg.write('image \n')

    im,hdr = pyfits.getdata(image, header=True)
    seg,shdr = pyfits.getdata(segmap, header=True)
    xx = np.random.randint(hdr['NAXIS1']-20, size=naper) + 10
    yy = np.random.randint(hdr['NAXIS2']-20, size=naper) + 10
 
    Aperture = CircularAperture(im, seg, xx, yy, rad)
    # perform photometry on uncontaminated aperture
    w = np.where(Aperture.clean)
    # aperture photometry
    phot = aper.aper(im, xx[w], yy[w], phpadu=1, apr=rad, zeropoint=zp, 
        exact=True, badpix=[0,0], setskyval=0.)
    mag,magerr,flux,fluxerr = phot[0],phot[1],phot[2],phot[3]
    make_reg(xx[w], yy[w], reg, rad)
#    fluxes = flux[fluxes != 0.0]
    reg.close()
    return flux


def fit_gaussian(flux, ax, left=False):
    '''Fit a Gaussian to the histogram of binned fluxes'

       The right side of the histogram may be contaminated by flux
       from sources (especially the wings). Setting left=True fits a 
       Gaussian to only the left side of the histogram. 
    '''
    # bin the fluxes
    nbins = 50 
    hist,bins = np.histogram(flux, bins=nbins)
    bincenters = 0.5 * (bins[1:]+bins[:-1])
    width = 1.01 * (bins[1] - bins[0])
    ax.bar(bincenters, hist, align='center', width=width, alpha=0.4,
           linewidth=0)

    # mirror left side to avoid contamination from sources
    if left:
        # find the peak of the histogram
        peak = np.argmax(hist)
        # left side of histogram
        left = hist[:peak+1]
        # duplicate and flip the left side
        right = left[:-1][::-1]
        # full histogram
        tot = np.append([left], [right])
        ntot = tot.shape[0]
        # if necessary, shorten bincenters to match new length of tot
        nbins = hist.shape[0]
        if ntot < nbins:
            bincenters = bincenters[0:ntot]
        # or lengthen bincenters
        if ntot > nbins:
            binsize = bins[1] - bins[0]
            bincenters = np.append(bincenters,np.arange(bincenters[-1]+binsize,
                                   bincenters[-1]+binsize*(ntot-nbins),binsize))
        # plot the new bins 
        width = 1.01 * (bincenters[1] - bincenters[0])
        ax.bar(bincenters, tot, align='center', width=width, alpha=0.5, 
               color='k', linewidth=0)
        hist = tot

    # fit a Gaussian
    gaussfunc = lambda x,a,mu,sigma: a * np.exp(-(x-mu)**2 / (2.*sigma**2))
    # initial estimates
    p0 = [np.max(hist), bincenters[len(bincenters)/2.], \
          (np.max(bincenters)-np.min(bincenters))/4.]
    popt,pcov = curve_fit(gaussfunc, bincenters, hist, p0=p0)
    # add to plot
    ax.plot(bincenters, gaussfunc(bincenters, *popt), 'k', alpha=0.5,
            linewidth=2)
    ax.set_xticks([np.min(bincenters),np.median(bincenters),np.max(bincenters)])
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    return popt


def find_maglim(image, segmap, zp, wispfield, rad=2.5, naper=5000):
    ''' '''
    print 'Finding limiting magnitude of %s'%image
    filt = pyfits.getheader(image)['FILTER'].rstrip("'")

    # set up plot
    fig = plt.figure(figsize=(8,6.5))
    gs = gridspec.GridSpec(2,2)
    gs.update(left=0.1, right=0.95, top=0.9, bottom=0.1,
              hspace=0.3, wspace=0.3)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
#    ax5 = fig.add_subplot(gs[1,1])
#    ax6 = fig.add_subplot(gs[1,2])
#    ax7 = fig.add_subplot(gs[2,0])
#    ax8 = fig.add_subplot(gs[2,1])
#    ax9 = fig.add_subplot(gs[2,2])
    axs = [ax1, ax2, ax3, ax4]#, ax5, ax6, ax7, ax8, ax9]

    fig.suptitle(r'%s %s band,  Histogram of fluxes'%(wispfield,filt))
    ax3.set_xlabel(r'Sky Flux [counts/s]', fontsize=15)
    ax3.set_ylabel(r'Number', fontsize=15)

    sigmas = np.zeros(4, dtype=float)
    for run in range(4):
        print '  Run %i'%(run+1)
        ax = axs[run]

        # photometry
        flux = photometry(image, segmap, naper, rad, zp)
        # remove outliers
        lowsig,upsig = np.percentile(flux, 15.9),np.percentile(flux,84.1)
        dist = upsig - lowsig
        flux = flux[(np.abs(flux-np.median(flux)) <= dist)]

        # bin the fluxes and fit a Gaussian
        popt = fit_gaussian(flux, ax, left=True)

        # popt[2] is the sigma of the distribution
        sigmas[run] = popt[2]

    fig.savefig(os.path.join(wispfield,'%s_%s_maglim.pdf'%(wispfield,filt)))
    sig = np.abs(np.median(sigmas))
    maglim = -2.5*np.log10(sig) + zp

    return maglim,sig

    
def main():
    parser = argparse.ArgumentParser(description=
        'Calculate the limiting 1sigma magnitude for an image')
    parser.add_argument('image', type=str,
        help='Image for which to calculate the limiting magnitude')
    parser.add_argument('zp', type=float,
        help='Zero point shift for the image')
    parser.add_argument('segmap', type=str,
        help='Segmentation map for image')
    parser.add_argument('wispfield', type=str,
        help='')
    parser.add_argument('--rad', type=float, default=2.5,
        help='Aperture radius')
    parser.add_argument('--naper', type=int, default=5000,
        help='Number of apertures to place randomly on image')
    args = parser.parse_args()
    image = args.image
    zp = args.zp
    segmap = args.segmap
    wispfield = args.wispfield
    rad = args.rad
    naper = args.naper

    find_maglim(image, segmap, zp, wispfield, rad, naper)


if __name__ == '__main__':
    main()
