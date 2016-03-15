#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
PyHdust auxiliary module: PyHdust BCD module.

:license: ?
"""
import numpy as _np
import pyhdust.phc as _phc

try:
    import matplotlib.pyplot as _plt
except:
    print('# Warning! Matplotlib module not installed!!!')

__author__ = "Daniel Moser; Ahmed Elshaer; Antoine Merand"
__email__ = "dmfaes@gmail.com"


def bcd(wav, nflx, logflx, elogflx=None, addsuf=None, doplot=False, 
    fmt=['png']):
    """ Calculate the Balmer_jump of a given spectral. 

    INPUT: wav (math:`\AA`), nflux (Normalized flux), log10(flux) and 
    log10(err_flux).

    Note that the continuum of nflux after the Balmer Jump is ~1.0 and before 
    is 0.30 < flux < 0.64!

    TODO: errors

    OUTPUT: offset, intersect """
    # data = _np.loadtxt(filename)
    logwav = _np.log10(wav)
    # -- first range:
    w1 = (logwav >= 3.5797) & (logwav <= 3.6989) & \
        (nflx >= .98) & (nflx <= 1.0)
    # -- second range
    w2 = (logwav >= 3.53) & (logwav <= 3.5658) & \
        (nflx >= .30) & (nflx <= .64)
    # -- linear fit for each range
    c1 = _np.polyfit(logwav[w1], logflx[w1], 1)
    c2 = _np.polyfit(logwav[w2], logflx[w2], 1)
    # -- computing offset at 3700A:
    x0 = _np.log10(3700.)
    offset = _np.polyval(c1, x0) - _np.polyval(c2, x0)
    print 'D (flux diff) at %4.0fA = %.3f' % (10**x0, offset)
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
    print 'intersection = %8.3fA' % (10**wlIntersect)

    if doplot:
    # -- plot data
        _plt.figure(0)
        _plt.clf()
        _plt.plot(logwav[logwav > 3.5051], logflx[logwav > 3.5051],
             alpha=0.3, color='k', label=addsuf)
        _plt.plot(logwav[w1], logflx[w1], '.r', label='range 1', alpha=0.1)
        _plt.plot(logwav[w2], logflx[w2], '.b', label='range 2', alpha=0.1)
        _plt.xlabel("log (wavelength A)")
        _plt.ylabel("log (flux)")
        if addsuf is not None:
            _plt.title(addsuf)
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
        _plt.ylim(-12.9, -11.4)

        if addsuf is None:
            addsuf = _phc.dtflag()
        for f in fmt:
            print('# Saved bcd{0}.{1}'.format(addsuf, f))
            _plt.savefig('bcd{0}.{1}'.format(addsuf, f), transparent=True)
        _plt.close()

    return offset, 10**wlIntersect


# MAIN
if __name__ == '__main__':
    pass
