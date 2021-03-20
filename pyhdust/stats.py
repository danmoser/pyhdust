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
    ax.plot(gausx, p, ls="-.", label='Gauss equiv.', color='gray')
    ax.set_xlim(xlim)
    ax.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=8, 
            labelspacing=0.05)
    ax.set_ylabel('c.d.f.')
    ax.set_ylim([0, 1])
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


def corr_coef(x, y, clear_nan=True):
    """ Pearson correlation coefficient for two ``x`` and ``y`` arrays (same 
    length).

    See also ``scipy.stats.pearsonr()``
    """
    if len(x) != len(y):
        _warn('# length of x and y arrays do not match!')
        return None
    x = _np.array(x)
    y = _np.array(y)
    if clear_nan:
        idx = _np.isnan(x)
        idy = _np.isnan(y)
        idxy = idx+idy
        x = x[~idxy]
        y = y[~idxy]
    avgx = _np.average(x)
    avgy = _np.average(y)
    return (_np.sum((x-avgx)*(y-avgy)))/(
        _np.sqrt(_np.sum((x-avgx)**2)*_np.sum((y-avgy)**2)) )


def corr_coef_spearman(x, y, clear_nan=True):
    """ Spearman's correlation coefficient for two ``x`` and ``y`` arrays (same 
    length).

    See also ``scipy.stats.spearmanr()``
    """
    if len(x) != len(y):
        _warn('# length of x and y arrays do not match!')
        return None
    x = _np.array(x)
    y = _np.array(y)
    if clear_nan:
        idx = _np.isnan(x)
        idy = _np.isnan(y)
        idxy = idx+idy
        x = x[~idxy]
        y = y[~idxy]
    idx = _np.argsort(x)
    x = x[idx]
    y = y[idx]
    n = len(x)
    diff = (_np.arange(n)-_np.argsort(y))
    diff2 = diff**2
    return 1 - 6*_np.sum(diff2)/(n*(n**2-1))


def corr_coef_cov(x, y, clear_nan=True):
    r""" Correlation coefficient based on the Covariance of two ``x`` and ``y`` 
    arrays (same length).

    :math:`\rho(x,y)=Cov(x,y)/sqrt(Var(x)*Var(y))`

    If :math:`\rho(x,y)= 0` we say that X and Y are "uncorrelated." If two 
    variables are independent, then their correlation will be 0. However, 
    like with covariance. it doesn't go the other way. A correlation of 0
    does not imply independence.
    """
    if len(x) != len(y):
        _warn('# length of x and y arrays do not match!')
        return None
    x = _np.array(x)
    y = _np.array(y)
    if clear_nan:
        idx = _np.isnan(x)
        idy = _np.isnan(y)
        idxy = idx+idy
        x = x[~idxy]
        y = y[~idxy]
    idx = _np.argsort(x)
    c, cov = _np.polyfit(x, y, 1, cov=True)
    return cov[0, 1]/_np.sqrt(cov[0, 0]*cov[1, 1])


def corr_coef_cov_with_err(x, y, yerr, xerr=None, clear_nan=True, 
    nsample=1000):
    r""" TO BE DONE
    Correlation coefficient based on the Covariance of two ``x`` and ``y`` 
    arrays (same length).

    :math:`\rho(x,y)=Cov(x,y)/sqrt(Var(x)*Var(y))`

    If :math:`\rho(x,y)= 0` we say that X and Y are "uncorrelated." If two 
    variables are independent, then their correlation will be 0. However, 
    like with covariance. it doesn't go the other way. A correlation of 0
    does not imply independence.
    """
    if xerr is None:
        xerr = _np.zeros(len(x))
    if len(x) != len(y) or len(x) != len(xerr) or len(xerr) != len(yerr):
        _warn('# length of input arrays do not match!')
        return None
    x = _np.array(x)
    xerr = _np.array(xerr)
    y = _np.array(y)
    yerr = _np.array(yerr)
    if clear_nan:
        idx = _np.isnan(x)
        idxe = _np.isnan(xerr)
        idy = _np.isnan(y)
        idye = _np.isnan(yerr)
        idxy = idx+idy+idxe+idye
        x = x[~idxy]
        xerr = xerr[~idxy]
        y = y[~idxy]
        yerr = yerr[~idxy]
    mx = _np.random.multivariate_normal(x, _np.diag(xerr**2), nsample)
    my = _np.random.multivariate_normal(y, _np.diag(yerr**2), nsample)
    return corr_coef_cov(mx.reshape(-1), my.reshape(-1))
