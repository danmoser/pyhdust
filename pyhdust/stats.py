# -*- coding:utf-8 -*-

"""PyHdust *stats* module: statistical tools

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import numpy as _np
import pyhdust.phc as _phc
import warnings as _warn

try:
    import matplotlib.pyplot as _plt
    from scipy.stats import mode as _mode
except ImportError:
    _warn.warn('Matplotlib and/or Scipy module not installed!!!')


__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def mad(data, axis=None):
    """ Return 1.48xMAD (median absolute deviation) 

    The MAD is a robust statistic, being more resilient to outliers in a data 
    set than the standard deviation."""
    return 1.4826*_np.median(_np.abs(data - _np.median(data, axis)), axis)


def summary(x, verbose=False):
    """ Returns the summary of the variable: "median", "minus sigma" and 
    "plus sigma" ROBUST values (i.e., median and [15.9, 84.1] percentiles). 

    Example:

    .. code::

        import pyhdust.stats as stt

        for i in range(8):
            a = _np.random.randn(10**i)+2
            print(np.average(a), np.std(a), stt.summary(a))
    """
    data = _np.hstack((_np.median(x), _np.percentile(x, (15.87, 84.13))))
    if verbose:
        print('# median and [15.9, 84.1] percentiles: ')
        print('# {0} {1} {2}'.format(*data))
    return data


def cdf(x, xlim=None, savefig=False):
    """ Display the CDF (Cumulative Density Distribution) of a sample `x`.

    A comparison with a gaussian and a linear one are made.
    """
    n = len(x)
    # probability
    p = _np.arange(n)/(n-1.)
    sortedx = _np.sort(x)
    if xlim is None:
        xlim = [sortedx[0], sortedx[-1]]
    linx = _np.linspace(xlim[0], xlim[1], n)
    madx = mad(x)
    gausx = _np.random.randn(n)*madx/2. + (xlim[0] + xlim[1])/2.
    gausx = _np.sort(gausx)

    fig, ax = _plt.subplots()
    ax.plot(sortedx, p, label='Data')
    ax.plot(linx, p, ls=':', label='linearized', color='gray')
    ax.plot(gausx, p, ls="--", label='Gauss equiv.', color='gray')
    ax.set_xlim(xlim)
    ax.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=8, 
            labelspacing=0.05)
    ax.set_ylabel('c.d.f.')
    if savefig:
        _phc.savefig(fig)  # figname='outname')
    return


def means(inarr, wtharr=None, quiet=False):
    """ Calculate many "means" for a given input array `inarr`.

    `wtharr` is the weights array (e.g., inverse of the uncertainty).

    Return simple, geom, harm, rms, median, mode
    """
    inarr = _np.array(inarr).astype(float)
    if wtharr is None:
        wtharr = _np.ones(len(inarr))
    simp = _np.sum(wtharr*inarr)/_np.sum(wtharr)
    geom = _np.exp(_np.sum(wtharr*_np.log(inarr))/_np.sum(wtharr))
    harm = _np.sum(wtharr)/_np.sum(wtharr/inarr)
    rms = _np.sqrt(_np.sum(wtharr*inarr**2)/_np.sum(wtharr))
    medi = _np.median(inarr)
    mode = _mode(inarr)[0][0]
    if not quiet:
        print('# Simple: {}'.format(simp))
        print('# Geometric: {}'.format(geom))
        print('# Harmonic: {}'.format(harm))
        print('# RMS: {}'.format(harm))
        print('# Median: {}'.format(medi))
        print('# Mode: {}'.format(mode))
    return simp, geom, harm, rms, medi, mode


def snr(count_rate, texp=1., nexp=1, npix=10., bg=10., dk=0., ron=2., var=0.):
    """Calcute the Signal-to-Noise ratio based on Poisson statistics.

    :param count_rate: = rate of counts (e-/time)
    :param npix: = number os pixels for the given count
    :param bg: = background rate per pixel (e-/time)
    :param dk: = dark rate per pixel (e-/time)
    :param ron: = readout noise (single pixel, in e-)
    :param var: = variance on the source erroes (e-)
    """
    return count_rate*texp*nexp / \
        _np.sqrt(texp*nexp*(count_rate+bg*npix+dk*npix) +ron**2*npix*nexp+ var)
