#!/usr/bin/env python
# -*- coding:utf-8 -*-

""" Program that shifts a FITS spectrum based on an line-profile Gaussian fit 
(optimized to work within a few tenths of km/s).
"""

from glob import glob
import numpy as np
import sys
import pyhdust.spectools as spt
import pyhdust.phc as phc
from argparse import ArgumentParser
import pyfits as pf
from scipy.optimize import curve_fit
from warnings import warn

__version__ = "0.92"
__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


class MyParser(ArgumentParser): 
    def error(self, message):
        sys.stderr.write('# ERROR! %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser(description=__doc__)
parser.add_argument('--version', action='version', 
    version='%(prog)s {0}'.format(__version__))
parser.add_argument("INPUT", help=("String to retrive the spectrum(a) "
    "filename(s) (wildcards accepted wrapped with \"\")"), type=str)
parser.add_argument("-l", "--line", action="store", dest="line", 
    help=("Wavelength of the reference line [default: %(default)s]"), 
    type=float, default=6562.79)

group2 = parser.add_argument_group('other arguments (overwrite power)')
group2.add_argument("-r", "--remove", action="store_true", dest="rm", 
    help=("If this flag is enabled, it restores the spectrum(a) to the "
        "original value (i.e., VELSHIFT = 0)"), default=False)
group2.add_argument("-f", "--force-dv", action="store", dest="vs", 
    help=("Force VELSHIFT to the non-zero VS value (in km/s) "
        "[default: None]"), type=float, default=0.)

args = parser.parse_args()


def gauss(x, a, x0, sigma, y0=1.):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+y0


def find_max_peak(vl, nflx, idx=-1):
    """ Based on max flux. 
    """
    if idx == -1:
        idx = len(vl)/2
    idx_ = phc.find_nearest(nflx[:idx], np.max(nflx[:idx]), idx=True)
    nfV, vlV = nflx[idx_], vl[idx_]
    idx_ = phc.find_nearest(nflx[idx:], np.max(nflx[idx:]), idx=True)
    nfR, vlR = nflx[idx+idx_], vl[idx+idx_]
    return nfV, vlV, nfR, vlR


def has_central_absorp(nf0, vlV, vlR):
    """ Return ``True`` or ``False`` if the line has a central absorption 
    (i.e., it is an absorption line or it is an double peak emission line).
    """
    if nf0 < 1:
        # print('# F(lb=0) < 1: A double-peak emission line was found!')
        return True
    if vlR-vlV > 100:
        return True
    else:
        return False


def gauss_fit(x, y, a0=None, x0=None, sig0=None, y0=None, emission=True):
    """ Return ``curve_fit``, i.e., ``popt, pcov``.
    """
    q = 5
    func = np.max
    if not emission:
        func = np.min
        q = 95
    if a0 is None:
        a0 = np.percentile(y, q) - np.median(y)
    if x0 is None:
        x0 = x[np.where(y == func(y))]
    if sig0 is None:
        sig0 = (np.max(x[-1])-np.min(x[0]))/30.
    if y0 is None:
        y0 = np.percentile(y, q)
    return curve_fit(gauss, x, y, p0=[a0, x0, sig0, y0])


def apply_shift(imfits, vshift, wl=None):
    if wl is None:
        # print('watch!', imfits[0].header['CDELT1'], len(imfits[0].data))
        wl = np.arange(len(imfits[0].data)) * imfits[0].header['CDELT1'] + \
            imfits[0].header['CRVAL1']

    new_wl = wl*(1 + vshift/phc.c.cgs*1e5)

    imfits[0].header['CDELT1'] = (new_wl[-1]-new_wl[0])/(len(new_wl)-1)
    imfits[0].header['CRVAL1'] = new_wl[0]
    # imfits[0].data = np.interp(new_wl, wl, flux)

    old_sh = 0
    if 'VELSHIFT' in imfits[0].header:
        old_sh = imfits[0].header['VELSHIFT'] 
    imfits[0].header['VELSHIFT'] = (vshift + old_sh, 
        'Lambda shift in km/s')

    print('# {0} updated with {1:.1f} km/s shift (accum. {2:.0f} km/s)'.
        format(fitsfile, vshift, vshift+old_sh))
    return

if __name__ == '__main__':

    lbc = args.line
    linput = args.INPUT

    for fitsfile in glob(linput):
        imfits = pf.open(fitsfile, mode='update')

        if args.rm is True:
            if 'VELSHIFT' in imfits[0].header:
                vshift = -imfits[0].header['VELSHIFT']
                apply_shift(imfits, vshift)
            else: 
                warn.warn('No VELSHIFT for {0}'.format(fitsfile))

        elif args.vs != 0:
            v0 = 0
            if 'VELSHIFT' in imfits[0].header:
                v0 = imfits[0].header['VELSHIFT']
            vshift = args.vs-v0
            apply_shift(imfits, vshift)

        else:
            try:
                flux = imfits[0].data
                wl = np.arange(len(flux)) * imfits[0].header['CDELT1'] + \
                    imfits[0].header['CRVAL1']

                vl, nflx = spt.lineProf(wl, flux, lbc=lbc)
                v0 = np.median( vl[nflx > np.percentile(nflx, 80)] )
                idx0 = phc.find_nearest(vl, v0, idx=True)
                nf0 = nflx[idx0]
                nfV, vlV, nfR, vlR = find_max_peak(vl, nflx, idx=idx0)

                if has_central_absorp(nf0, vlV, vlR):
                    idxV = phc.find_nearest(vl, vlV, idx=True)
                    idxR = phc.find_nearest(vl, vlR, idx=True)
                    popt, pcov = gauss_fit(vl[idxV:idxR], nflx[idxV:idxR], 
                        emission=False)
                else:
                    popt, pcov = gauss_fit(vl, nflx)

                vshift = -popt[1]
                apply_shift(imfits, vshift, wl)
            except RuntimeError:
                warn.warn('Gaussian fit did not work for {0}'.format(fitsfile))
                vshift = phc.user_input('Enter a vshift value (km/s): ')
                try:
                    vshift = float(vshift)
                    apply_shift(imfits, vshift, wl)
                except ValueError:
                    warn.warn("You entered a invalid value! File unchanged.")

        imfits.close()

    if len(glob(linput)) == 0:
        warn.warn('{0} files were not found!'.format(linput))
