# -*- coding:utf-8 -*-

"""PyHdust *stats* module: statistical tools

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import numpy as _np
import pyhdust.phc as _phc
import matplotlib.pyplot as _plt

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
