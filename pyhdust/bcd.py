# -*- coding:utf-8 -*-

"""PyHdust *bcd* auxiliary module: PyHdust BCD module.

.. code: python

    result = bcd.analyze_all()
    print '-' * 30
    for k in result.keys():
        print k, result[k]
        # _np.savetxt('result',(k))

:co-author: Bruno Mota; Daniel Moser; Ahmed Elshaer 
:license: Copyright 2015 Antoine Merand
"""
from __future__ import print_function
import os as os
import numpy as _np
import sys as _sys


def eprint(*args, **kwargs):
    print(*args, file=_sys.stderr, **kwargs)
    return

try:
    import matplotlib.pyplot as _plt
except ImportError:
    eprint('# Warning! Matplotlib module not installed!!!')

__author__ = "Antoine Merand"
__email__ = "amerand@eso.org"


def balmer_jump(filename):
    r""" Calculate the Balmer_jump of a given spectral. 

    INPUT: `filename` is a file with 4 columns, as [wav (:math:`\AA`), 
    log10(wav), norm_flux, log(flux)] 

    OUTPUT: offset, intersect """
    data = _np.loadtxt(filename)
    # -- first range:
    w1 = (data[:, 1] >= 3.5797) * (data[:, 1] <= 3.6989) * \
        (data[:, 2] >= .98) * (data[:, 2] <= 1.0)
    # -- second range
    w2 = (data[:, 1] >= 3.53) * (data[:, 1] <= 3.5658) * \
        (data[:, 2] >= .30) * (data[:, 2] <= .64)
    # -- linear fit for each range
    c1 = _np.polyfit(data[:, 1][w1], data[:, 3][w1], 1)
    c2 = _np.polyfit(data[:, 1][w2], data[:, 3][w2], 1)
    # -- computing offset at 3700A:
    x0 = _np.log10(3700.)
    offset = _np.polyval(c1, x0) - _np.polyval(c2, x0)
    # print 'offset at %4.0fA = %.3f' % (10**x0, offset)
    # -- balmer series, using the Rydberg formula
    B = 3645.0682
    m = 2.0  # for Balmer
    n = _np.arange(10, 18)
    wl = B * (n**2 / (n**2 - m**2))
    # -- maxima as middle points between minima
    wlMax = 0.5 * (wl[1:] + wl[:-1])
    # -- distance between minima:
    wlStep = _np.abs(_np.diff(wl))
    # -- fit all maxima
    maxima = []
    for i in range(len(wlMax)):
        wl0 = wlMax[i]
        dwl = wlStep[i] / 6.
        x = data[:, 1][_np.abs(data[:, 0] - wl0) < dwl]
        y = data[:, 3][_np.abs(data[:, 0] - wl0) < dwl]
        c = _np.polyfit(x, y, 2)
        # -- actual position of maximum
        wl1 = -0.5 * c[1] / c[0]
        # -- (wl at max, max, poly coef, width of the fit)
        maxima.append((wl1, _np.polyval(c, wl1), c, dwl))
    # -- plot data
    _plt.figure(0)
    _plt.clf()
    _plt.plot(data[:, 1][data[:, 1] > 3.5051], data[:, 3][data[:, 1] > 3.5051],
             alpha=0.3, color='k', label=filename)
    _plt.plot(data[:, 1][w1], data[:, 3][w1], '.r', label='range 1', alpha=0.1)
    _plt.plot(data[:, 1][w2], data[:, 3][w2], '.b', label='range 2', alpha=0.1)
    _plt.xlabel("log (wavelength A)")
    _plt.ylabel("log (flux)")
    _plt.title('HD91373')
    # -- plot linear fit:
    x = _np.linspace(x0 - 0.06, x0 + 0.2, 1000)
    # x = _np.linspace(x0-500, x0+1500, 1000)

    _plt.plot(x, _np.polyval(c1, x), '-r')
    _plt.plot(x, _np.polyval(c2, x), '-b')
    # _plt.legend()

    # -- Balmer Jump
    _plt.plot([x0, x0], [_np.polyval(c1, x0), _np.polyval(c2, x0)], '-g')
    # (red line+blue line)/2
    intersect = 0.5 * (_np.polyval(c1, x) + _np.polyval(c2, x))

    maximaCur = _np.interp(
        x, [m[0] for m in maxima][::-1], [m[1] for m in maxima][::-1])
    _plt.plot([m[0] for m in maxima][::-1], [m[1] for m in maxima][::-1], 'y')
    # _plt.plot(x, maximaCur, '-g')
    wlIntersect = x[_np.argmin(_np.abs(intersect - maximaCur))]
    # print 'intersection = %8.3fA' % (10**wlIntersect)

    # -- position of maxima
    _plt.plot(x, intersect, color='orange', linestyle='dashed')

    # -- top of the line polynomial fit:
    for m in maxima:
        x = _np.linspace(m[0] - m[3], m[0] + m[3], 1.)
        _plt.plot(x, _np.polyval(m[2], x), '-', color='orange')
    _plt.xlim(3.5, 3.8)
    # _plt.ylim(-12.9, -11.4)

    return offset, 10**wlIntersect


def analyze_all(fmt=['png']):
    """ Run `balmer_jump` function to all files starting with 'HD'. """
    filenames = os.listdir('./')
    filenames = filter(
        lambda x: x.startswith('HD') and not x.endswith('.pdf'), filenames)
    # print(filenames)
    res = {}
    for f in filenames:
        print('*' * 5, f, '*' * 5)
        res[f] = balmer_jump(f)
        _plt.figure(0)
        # _plt.savefig(f+'.png')
        for ext in fmt:
            _plt.savefig('{0}.{1}'.format(f, ext), dpi=600, transparent=True)
    return res


def bcd(obj, wav, nflx, logflx, elogflx=None, label='Spec', folder=None, 
    doplot=False):
    r""" Calculate the Balmer_jump of a given spectral.

    INPUT: wav (:math:`\AA`), nflux (Normalized flux), log10(flux) and
    log10(err_flux).

    Note that the continuum of nflux after the Balmer Jump is ~1.0 and before
    is 0.30 < flux < 0.64!

    TODO: errors

    OUTPUT: offset, intersect """
    # data = _np.loadtxt(filename)
    logwav = _np.log10(wav)

    # print('logwav')
    # print(logwav)

    # print('nflx')
    nflx = _np.array(nflx)
    # print(nflx)

    # -- first range:
    w1 = (logwav >= 3.5797) & (logwav <= 3.6989) 
    logwav_w1 = logwav[w1]
    logflx_w1 = logflx[w1]
    nflx_w1 = nflx[w1]
    # print('log_w1')
    # print(len(logwav_w1))
    # print(len(logflx_w1))

    w1 = (nflx_w1 >= .98) & (nflx_w1 <= 1.02) & _np.logical_not(_np.isnan(
        logflx_w1))
    logwav_w1 = logwav_w1[w1]
    logflx_w1 = logflx_w1[w1]
    # print('log_w1')
    # print(len(logwav_w1))
    # print(len(logflx_w1))

    w2 = (logwav >= 3.53) & (logwav <= 3.5658)
    logwav_w2 = logwav[w2]
    logflx_w2 = logflx[w2]
    nflx_w2 = nflx[w2]

    # -- second range
    w2 = (nflx_w2 >= .98) & (nflx_w2 <= 1.02) & _np.logical_not(_np.isnan(
        logflx_w2))
    logwav_w2 = logwav_w2[w2]
    logflx_w2 = logflx_w2[w2]

    #  -- linear fit for each range
    # print(len(logwav[w1]))
    # print(len(logflx[w1]))
    # c1 = _np.polyfit(logwav[w1], logflx[w1], 1)
    # print(logwav_w1)
    # print(logflx_w1)

    c1 = _np.polyfit(logwav_w1, logflx_w1, 1)
    c2 = _np.polyfit(logwav_w2, logflx_w2, 1)
    # -- computing offset at 3700A:
    x0 = _np.log10(3700.)
    offset = _np.polyval(c1, x0) - _np.polyval(c2, x0)
    # print 'D (flux diff) at %4.0f A = %.3f' % (10**x0, offset)
    # -- balmer series, using the Rydberg formula
    B = 3645.0682
    m = 2.0  # for Balmer
    n = _np.arange(10, 18)
    wl = B * (n**2 / (n**2 - m**2))
    # -- maxima as middle points between minima
    wlMax = 0.5 * (wl[1:] + wl[:-1])
    # -- distance between minima:
    wlStep = _np.abs(_np.diff(wl))
    # -- fit all maxima
    maxima = []
    for i in range(len(wlMax)):
        wl0 = wlMax[i]
        dwl = wlStep[i] / 6.
        x = logwav[_np.abs(wav - wl0) < dwl]
        y = logflx[_np.abs(wav - wl0) < dwl]
        c = _np.polyfit(x, y, 2)
        # -- actual position of maximum
        wl1 = -0.5 * c[1] / c[0]
        # -- (wl at max, max, poly coef, width of the fit)
        maxima.append((wl1, _np.polyval(c, wl1), c, dwl))

    x = _np.linspace(x0 - 0.06, x0 + 0.2, 1000)
    # (red line+blue line)/2
    intersect = 0.5 * (_np.polyval(c1, x) + _np.polyval(c2, x))
    maximaCur = _np.interp(
        x, [m[0] for m in maxima][::-1], [m[1] for m in maxima][::-1])
    wlIntersect = x[_np.argmin(_np.abs(intersect - maximaCur))]
    # print 'intersection = %8.3f A' % (10**wlIntersect)

    if doplot:
    # -- plot data
        _plt.figure(0)
        _plt.clf()
        _plt.plot(logwav[logwav > 3.5051], logflx[logwav > 3.5051],
             alpha=0.3, color='k', label=label)
        _plt.plot(logwav[w1], logflx[w1], '.r', label='range 1', alpha=0.1)
        _plt.plot(logwav[w2], logflx[w2], '.b', label='range 2', alpha=0.1)
        _plt.xlabel(r'$\log \lambda [\AA]$')
        _plt.ylabel(r'$\log F_\lambda$')
        # _plt.title('HD91373')
        # -- plot linear fit:
        # x = _np.linspace(x0-500, x0+1500, 1000)
        _plt.plot(x, _np.polyval(c1, x), '-r')
        _plt.plot(x, _np.polyval(c2, x), '-b')
        # _plt.legend()
        # -- Balmer Jump
        _plt.plot([x0, x0], [_np.polyval(c1, x0), _np.polyval(c2, x0)], '-g')
        _plt.plot([m[0] for m in maxima][::-1], [m[1] for m in maxima][::-1],
            'y')
        # _plt.plot(x, maximaCur, '-g')
        # -- position of maxima
        _plt.plot(x, intersect, color='orange', linestyle='dashed')
        # -- top of the line polynomial fit:
        for m in maxima:
            x = _np.linspace(m[0] - m[3], m[0] + m[3], 1.)
            _plt.plot(x, _np.polyval(m[2], x), '-', color='orange')
        _plt.xlim(3.5, 3.8)
        # _plt.ylim(-12.9, -10.4)
        _plt.ylim(min(logflx[logwav > 3.5051]), max(logflx[logwav > 3.5051]))
        _plt.minorticks_on()
        _plt.tight_layout()
        # _plt.autoscale(axis='both')
        # _plt.autoscale(axis='y')
        # _plt.show()
        _plt.savefig(str(folder)+obj+'_adjust_bcd.png')
        _plt.close()

    return offset, 10**wlIntersect

# MAIN
if __name__ == '__main__':
    pass
