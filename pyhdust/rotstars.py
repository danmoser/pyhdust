# -*- coding:utf-8 -*-

"""
PyHdust module: Rotating stars tools.

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
import re as _re
import numpy as _np
import pyhdust.phc as _phc

# try:
#     import matplotlib.pyplot as _plt
#     from scipy import interpolate as _interpolate
# except:
#     print('# Warning! matplotlib and/or scipy module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def readscr(scrfile):
    '''
    Read source generated with `ref_estrela.txt`.

    OUTPUT: M, Req and TP (2*solar units and K).
    '''
    f0 = open(scrfile)
    lines = f0.readlines()
    f0.close()

    n = int(_phc.fltTxtOccur('STAR =', lines, n=1))
    M = _phc.fltTxtOccur('M =', lines, n=n)
    Rp = _phc.fltTxtOccur('R_pole =', lines, n=n)
    if n == 2:
        ob = _phc.fltTxtOccur('R_eq/R_pole =', lines, n=1)
        Tp = _phc.fltTxtOccur('Teff_pole =', lines, n=1)
    else:
        W = _phc.fltTxtOccur('W =', lines, n=1)
        bet = _phc.fltTxtOccur('Beta_GD =', lines, n=1)
        L = _phc.fltTxtOccur('L =', lines, n=n)
        wfrac = _np.sqrt(27. / 8 * (1 + 0.5 * W**2)**-3 * W**2)
        ob, Tp, A = rotStar(Tp=L, M=M, rp=Rp, beta=bet, wfrac=wfrac, 
            quiet=True, LnotTp=True)
    # print M,Rp*ob,Tp
    return M, Rp * ob, Tp


def vrot_scr(scrfile):
    """ Returns the `vrot` value of a given source star.

    OUTPUT: vrot in km/s. """
    M, Req, Tp = readscr(scrfile)
    # Be_M04.80_ob1.40_H0.30_Z0.014_bE_Ell
    rule = '(?<=_ob)(.+)(?=_H)'
    ob = float(_re.findall(rule, scrfile)[0])
    vrot = wrot(ob, is_ob=True) * \
        _np.sqrt(_phc.G.cgs * _phc.Msun.cgs * M / Req / _phc.Rsun.cgs)
    return vrot*1e-5


def wrot(par, is_ob=False):
    """ Converts math:`w_{\rm frac} = \Omega/\Omega_c` into 
    math:`W = vrot/vorb`.

    If `is_ob == True`, it consider the param as ob (instead of 
    math:`w_{\rm frac}`). """
    if is_ob:
        wfrac = (1.5 ** 1.5) * _np.sqrt(2. * (par - 1.) / par ** 3)
    else: 
        wfrac = par

    gam = 2. * _np.cos((_np.pi + _np.arccos(wfrac)) / 3.)
    W = _np.sqrt(gam ** 3 / wfrac)

    return W


def beta(par, is_ob=False):
    """ Calculate the math:`\beta` value from Espinosa-Lara for a given 
    rotation rate math:`w_{\rm frac} = \Omega/\Omega_c`

    If `is_ob == True`, it consider the param as ob (instead of
    math:`w_{\rm frac}`). """

    # Ekstrom et al. 2008, Eq. 9
    if is_ob:
        wfrac = (1.5 ** 1.5) * _np.sqrt(2. * (par - 1.) / par ** 3)
    else: 
        wfrac = par

    # avoid exceptions
    if wfrac == 0:
        return .25
    elif wfrac == 1:
        return 0.13535
    elif wfrac < 0 or wfrac > 1:
        print('# Warning! Invalid value of wfrac.')
        return 0.

    # Espinosa-Lara VLTI-School 2013 lecture, slide 18...
    delt = 1.
    omega1 = 0.
    omega = wfrac
    while delt >= 1e-5:
        f = (3. / (2. + omega**2))**3 * omega**2 - wfrac**2
        df = -108. * omega * (omega**2 - 1.) / (omega**2 + 2.)**4
        omega1 = omega - f / df
        delt = _np.abs(omega1 - omega) / omega
        omega = omega1

    nthe = 100
    theta = _np.linspace(0, _np.pi / 2, nthe + 1)[1:]
    grav = _np.zeros(nthe)
    teff = _np.zeros(nthe)
    corr = _np.zeros(nthe)
    beta = 0.

    for ithe in range(nthe):

        delt = 1.
        r1 = 0.
        r = 1.
        while delt >= 1e-5:
            f = omega**2 * r**3 * \
                _np.sin(theta[ithe])**2 - (2. + omega**2) * r + 2.
            df = 3. * omega**2 * r**2 * \
                _np.sin(theta[ithe])**2 - (2. + omega**2)
            r1 = r - f / df
            delt = _np.abs(r1 - r) / r
            r = r1

        delt = 1.
        n1 = 0.
        ftheta = 1. / 3. * omega**2 * r**3 * _np.cos(theta[ithe])**3 + \
            _np.cos(theta[ithe]) + _np.log(_np.tan(theta[ithe] / 2.))
        n = theta[ithe]
        while delt >= 1e-5:
            f = _np.cos(n) + _np.log(_np.tan(n / 2.)) - ftheta
            df = -_np.sin(n) + 1. / _np.sin(n)
            n1 = n - f / df
            delt = abs(n1 - n) / n
            n = n1

        grav[ithe] = _np.sqrt(1. / r**4 + omega**4 * r**2 * _np.sin(
            theta[ithe])**2 - 2. * omega**2 * _np.sin(theta[ithe])**2 / r)

        corr[ithe] = _np.sqrt(_np.tan(n) / _np.tan(theta[ithe]))
        teff[ithe] = corr[ithe] * grav[ithe]**0.25

    u = ~_np.isnan(teff)
    coef = _np.polyfit(_np.log(grav[u]), _np.log(teff[u]), 1)
    beta = coef[0]

    return beta


def rotStar(Tp=20000., M=10.3065, rp=5.38462, star='B', beta=0.25, wfrac=0.8,
            th_res=5001, quiet=False, LnotTp=False):
    """ Return the photospheric parameters of a rotating star.

    `LnotTp`: the value of "Tp" is the Luminosity (in solar units).

    Calculation of Von Zeipel's Beta parameter as function of W: see math...

    INPUT: th_res (theta resolution, integer)...

    OUTPUT: printed status + (ob, Tp values, Area[cm2])
    """
    Rsun = _phc.Rsun.cgs
    Msun = _phc.Msun.cgs
    Lsun = _phc.Lsun.cgs
    G = _phc.G.cgs
    # AU = _phc.au.cgs
    # pc = _phc.pc.cgs
    sigma = _phc.sigma.cgs
    M = M * Msun
    rp = rp * Rsun
    if wfrac == 0.:
        wfrac = 1e-9
    if LnotTp:
        Tp = (Tp * Lsun / 4. / _np.pi / rp**2 / sigma)**.25

    # DEFS ###
    def rt(th, wfrac):
        if th == 0:
            r = 1.
        else:
            r = (-3. * _np.cos((_np.arccos(wfrac * _np.sin(th)) + 4 *
                _np.pi) / 3)) / (wfrac * _np.sin(th))
        return r

    def area(wfrac):
        ths = _np.linspace(_np.pi / 2, 0, th_res)
        a = 0.
        for i in range(len(ths)):
            a = a + 2 * _np.pi * rt(ths[i], wfrac) ** 2 * _np.sin(ths[i])
        return 2 * a * ths[-2]

    def g(wfrac, M, rp, th):
        wcrit = _np.sqrt(8 * G * M / (27 * rp ** 3))
        g = (wcrit * wfrac) ** 2 * rp * rt(th, wfrac) * \
            _np.sin(th) ** 2 - G * M / (rp * rt(th, wfrac)) ** 2
        return g

    def lum(wfrac, Tp, rp, M, C, beta):
        ths = _np.linspace(_np.pi / 2, 0, th_res)
        l = 0.
        for i in range(len(ths)):
            l = l + rt(ths[i], wfrac) ** 2 * _np.sin(ths[i]) * \
                (abs(g(wfrac, M, rp, ths[i]))) ** (4 * beta)
        return 2 * 2 * _np.pi * ths[-2] * sigma * rp ** 2 * C ** (4 * beta) * l

    def lumf(wfrac, Tp, rp, M, beta):
        ths = _np.linspace(_np.pi / 2, 0, th_res)
        l = 0.
        for i in range(len(ths)):
            l = l + rt(ths[i], wfrac) ** 2 * _np.sin(ths[i]) * \
                abs(g(wfrac, M, rp, ths[i])) ** (4 * beta)
        return l * ths[-2] * rp ** 2

    Bstars = _np.array(_phc.bestars, dtype=str)
    if star in Bstars:
        i = _np.where(Bstars[:, 0] == star)
        i = i[0][0]
        print Bstars[i][0]
        Tp = float(Bstars[i][1])
        M = float(Bstars[i][2]) * Msun
        rp = float(Bstars[i][3]) * Rsun
        # comentar linha abaixo se 1a. rodada:
        # Tp = 27438.63 #K

    wcrit = _np.sqrt(8 * G * M / (27 * rp ** 3))
    C = Tp ** (1. / beta) / abs(G * M / rp ** 2)

    vrot = wcrit * wfrac * rp * rt(_np.pi / 2, wfrac)
    lum0 = 4 * _np.pi * rp ** 2 * sigma * Tp ** 4 / Lsun

    # a = rp**2*Tp**4*abs(g(wfrac,M,rp,0.))**(4*beta)
    # print('Teff_pol* = %.2f' % ( (a/b)**beta ) )
    b = lumf(wfrac, Tp, rp, M, beta)
    c = lumf(0.0001, Tp, rp, M, beta)
    Cw = (c / b) ** (1. / (4. * beta)) * C
    ob = rt(_np.pi / 2, wfrac)  # /(rp / Rsun)

    # OUTPUT ###
    if not quiet:
        print('# Parameters:')
        print('wfrac     = %.4f' % (wfrac))
        print('W         = %.4f' % (_np.sqrt(2 * (ob - 1))))
        print('Star Mass = %.2f Msun' % (M / Msun))
        print('Rpole     = %.2f Rsun' % (rp / Rsun))
        print('Req       = %.2f Rpole' % (rt(_np.pi / 2, wfrac)))
        print('Teff_pol  = %.1f' % (Tp))

        print('Star Area = %.2f' % (area(wfrac)))
        print('Star Lum. = %.1f' % (lum(wfrac, Tp, rp, C, M, beta) / Lsun))
        print('Star Lum.*= %.1f' % (lum0))

        print('vrot(km/s)= %.1f' % (vrot / 1e5))
        print('vorb(km/s)= %.1f' %
              (_np.sqrt(G * M / rp / rt(_np.pi / 2, wfrac)) / 1e5) )
        print('vcrt(km/s)= %.1f' % (wcrit * rp * rt(_np.pi / 2, 1.) / 1e5))

        print('log(g)pole= %.2f' % (_np.log10(abs(g(wfrac, M, rp, 0.))) ))
        print('log(g)eq  = %.2f' %
              (_np.log10(abs(g(wfrac, M, rp, _np.pi / 2))) ))
        print('Teff_eq   = %.1f' %
              ( (C * abs(g(wfrac, M, rp, _np.pi / 2))) ** beta) )
        print('Teff_eq*  = %.1f' %
              ( (Cw * abs(g(wfrac, M, rp, _np.pi / 2))) ** beta) )

        print('Teff_pol* = %.2f' % ( (Cw * abs(g(wfrac, M, rp, 0.))) ** beta) )
        print('T_pol/eq* = %.4f' % ((Cw * abs(g(wfrac, M, rp, 0.))) ** beta / 
            (Cw * abs(g(wfrac, M, rp, _np.pi / 2))) ** beta) )

        print('# \"*\" == case where L is constant!')
    return ob, (Cw * abs(g(wfrac, M, rp, 0.))) ** beta, area(wfrac) * (rp**2)


# MAIN ###
if __name__ == "__main__":
    pass
