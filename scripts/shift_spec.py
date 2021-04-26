#!/usr/bin/env python3
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
# from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from warnings import warn
from lmfit import Model

__version__ = "1.0"
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
parser.add_argument("-w", "--half-width", action="store", dest="hw", 
    help=("Half-width in km/s of the reference line [default: %(default)s]"), 
    type=float, default=6562.79)

group2 = parser.add_argument_group('other arguments (overwrite power)')
group2.add_argument("-r", "--remove", action="store_true", dest="rm", 
    help=("If this flag is enabled, it restores the spectrum(a) to the "
        "original value (i.e., VELSHIFT = 0)"), default=False)
group2.add_argument("-f", "--force-dv", action="store", dest="vs", 
    help=("Force VELSHIFT to the non-zero VS value (in km/s) "
        "[default: None]"), type=float, default=0.)

args = parser.parse_args()


def gauss(x, a, x0, sigma):
    y0=1.
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


def has_central_absorp(nf0, vlV, vlR, hw):
    """ Return ``True`` or ``False`` if the line has a central absorption 
    (i.e., it is an absorption line or it is an double peak emission line).
    """
    if nf0 < 1:
        # print('# F(lb=0) < 1: A double-peak emission line was found!')
        return True
    if vlR-vlV > hw/10:
        return True
    else:
        return False


def gauss_fit(x, y, a0=None, x0=None, sig0=None, emission=True):
    """ Return ``curve_fit``, i.e., ``popt, pcov``.

    def gauss_fit(x, y, a0=None, x0=None, sig0=None, emission=True, ssize=0.05):
    
    # def gauss(x, a, x0, sigma):
    #     y0=1.
    #     return a*_np.exp(-(x-x0)**2/(2*sigma**2))+y0
    if ssize < 0 or ssize > .5:
        _warn.warn('Invalid ssize value...', stacklevel=2)
        ssize = 0
    ssize = int(ssize * len(y))
    if ssize == 0:
        ssize = 1

    q = 95
    func = _np.max
    if not emission:
        func = _np.min
        q = 5
    if a0 is None:
        a0 = _np.abs(_np.percentile(y, q)) - _np.median(y)
    if x0 is None:
        x0 = x[_np.where(y == func(y))]
    if sig0 is None:
        sig0 = (_np.max(x)-_np.min(x))/10.
    # if y0 is None:
    #     y0 = np.median(y)
    # gmodel = _Model(gauss)
    # gmodel.set_param_hint('a', min=0.2, max=20)
    # if not emission:
    #     gmodel.set_param_hint('a', min=-0.2, max=-4)
    # gmodel.set_param_hint('sigma', min=50, max=1000)
    # result = gmodel.fit(y, x=x, a=a0, x0=x0, sigma=sig0)
    # print(result.params['a'], result.params['sigma'], result.params['x0'])
    # return result.params['a']*_np.sqrt(_np.pi*2)*result.params['sigma']
    medx0, medx1 = _np.average(x[:ssize]), _np.average(x[-ssize:])
    if ssize > 9:
        medy0, medy1 = _np.median(y[:ssize]), _np.median(y[-ssize:])
    else:
        medy0, medy1 = _np.average(y[:ssize]), _np.average(y[-ssize:])
    new_y = medy0 + (medy1 - medy0) * (x - medx0) / (medx1 - medx0)
    g_init = _models.Gaussian1D(amplitude=a0, mean=x0, stddev=sig0)
    fit_g = _fitting.LevMarLSQFitter()
    g = fit_g(g_init, x, y-new_y)
    """
    q = 95
    func = np.max
    if not emission:
        func = np.min
        q = 5
    if a0 is None:
        a0 = np.abs(np.percentile(y, q)) - np.median(y)
    if x0 is None:
        x0 = x[np.where(y == func(y))]
    if sig0 is None:
        sig0 = (np.max(x)-np.min(x))/10.
    # if y0 is None:
    #     y0 = np.median(y)
    gmodel = Model(gauss)
    gmodel.set_param_hint('a', min=0.2, max=4)
    if not emission:
        gmodel.set_param_hint('a', min=-0.2, max=-4)
    gmodel.set_param_hint('sigma', min=50, max=1000)
    result = gmodel.fit(y, x=x, a=a0, x0=x0, sigma=sig0)

    fig, (ax0, ax1) = plt.subplots(2, 1)
    ax0.plot(x, y, 'bo')
    ax0.plot(x, result.init_fit, 'k--')
    ax0.plot(x, result.best_fit, 'r-')
    ax0.set_title(fitsfile)
    print(lbc, np.min(x))
    idx = np.where((wl > lbc*0.98) & (wl < lbc*1.03))
    ax1.plot(wl[idx], flux[idx], 'o')
    plt.show(block=False)
    cmd = phc.user_input('# Problem (y/other)? ')
    if cmd.lower().startswith('y'):
        raise ValueError
    # phc.savefig(fig, figname=fitsfile)
    return result.params['x0']


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
                warn('No VELSHIFT for {0}'.format(fitsfile))

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

                vl, nflx = spt.lineProf(wl, flux, lbc=lbc, hwidth=args.hw)
                nfxsig = np.std(nflx)
                emission = True
                if np.percentile(nflx, 5) + nfxsig < 1:
                    emission = False
                    if np.percentile(nflx, 95) - 1.5*nfxsig > 1:
                        emission = True

                vshift = -gauss_fit(vl, nflx, emission=emission)
                if np.abs(vshift) > args.hw:
                    raise ValueError('VelShift out of bounds!')
                apply_shift(imfits, vshift, wl)
            except:
                warn('Gaussian fit did not work for {0}'.format(fitsfile))
                vshift = phc.user_input('Enter a vshift value (km/s): ')
                try:
                    vshift = float(vshift)
                    apply_shift(imfits, vshift, wl)
                except ValueError:
                    warn("You entered a invalid value! File unchanged.")

        imfits.close()

    if len(glob(linput)) == 0:
        warn('{0} files were not found!'.format(linput))
