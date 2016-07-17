# -*- coding:utf-8 -*-

"""PyHdust *rotstars* module: Rotating stars tools.

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import re as _re
import numpy as _np
import pyhdust.phc as _phc
import warnings as _warn

# try:
#     import matplotlib.pyplot as _plt
#     from scipy import interpolate as _interpolate
# except:
#     print('# Warning! matplotlib and/or scipy module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def readscr(scrfile):
    ''' Read source generated with *ref_estrela.txt*.

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
    """ Returns the ``vrot`` value of a given source star.

    OUTPUT: vrot in km/s. """
    M, Req, Tp = readscr(scrfile)
    # Be_M04.80_ob1.40_H0.30_Z0.014_bE_Ell
    rule = '(?<=_ob)(.+)(?=_H)'
    ob = float(_re.findall(rule, scrfile)[0])
    vrot = wrot(ob, is_ob=True) * \
        _np.sqrt(_phc.G.cgs * _phc.Msun.cgs * M / Req / _phc.Rsun.cgs)
    return vrot*1e-5


def wrot(par, is_ob=False):
    r""" Converts :math:`w_{\rm frac} = \Omega/\Omega_c` into 
    :math:`W = vrot/vorb`.

    If ``is_ob == True``, it considers the param as the oblateness (instead of 
    :math:`w_{\rm frac}`). """
    if is_ob:
        wfrac = (1.5 ** 1.5) * _np.sqrt(2. * (par - 1.) / par ** 3)
    else: 
        wfrac = par

    gam = 2. * _np.cos((_np.pi + _np.arccos(wfrac)) / 3.)
    W = _np.sqrt(gam ** 3 / wfrac)

    return W


def wfrac_rot(W):
    """ Returns wfrac (Omega/Omega_crit) value from a W value.

    Equation 1.23 de Faes (2015).
    """
    return _np.sqrt(27/8.*W**2/(1+.5*W**2)**3)


def beta(par, is_ob=False):
    r""" Calculate the :math:`\beta` value from Espinosa-Lara for a given 
    rotation rate :math:`w_{\rm frac} = \Omega/\Omega_c`

    If ``is_ob == True``, it consider the param as ob (instead of
    :math:`w_{\rm frac}`). """

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
        _warn.warn('Invalid value of wfrac.')
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

    ``LnotTp``: the value of "Tp" is the Luminosity (in solar units).

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

    if star.startswith('B'):
        Bstars = _np.array(bestarsHarm1988, dtype=str)
        if star in Bstars:
            i = _np.where(Bstars[:, 0] == star)
            i = i[0][0]
            print(Bstars[i][0])
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


def rochearea(wfrac, isW=False):
    """ Calculate the Roche area of a rigid rotator.

    Equation 4.23 from Cranmer 1996 (thesis).

    Area in (squared) radial unit (it must be multiplied to Rpole**2 to a 
    physical size).
    """
    if isW:
        w = wfrac_rot(wfrac)
    else:
        w = wfrac
    return 4*_np.pi*(1+.19444*w**2+0.28053*w**2-1.9014*w**6+6.8298*w**8-
        9.502*w**10+4.6631*w**12)


bestarsHarm1988 = [
    # The numbers below are based on Harmanec 1988
    # B1.5 and B2.5 interpolated by Faes.
    # Teff fixed: Rp2 from Lum1; Lum2 from Rp1.
    # SpType  Teff    Mass     Rp    Lum   Rp2     Lum2
    ['B0.0', 29854., 14.57, 05.80, 27290., 6.19, 23948.8487173],
    ['B0.5', 28510., 13.19, 05.46, 19953., 5.80, 17651.9502267],
    ['B1.0', 26182., 11.03, 04.91, 11588., 5.24, 10152.9628687],
    ['B1.5', 24599., 09.72, 04.58, 07768., 4.87, 6883.65832266],
    ['B2.0', 23121., 08.62, 04.28, 05297., 4.55, 4691.72482578],
    ['B2.5', 20980., 07.18, 03.90, 02931., 4.11, 2641.00783143],
    ['B3.0', 19055., 06.07, 03.56, 01690., 3.78, 1497.45695726],
    ['B4.0', 17179., 05.12, 03.26, 00946., 3.48, 829.555139678],
    ['B5.0', 15488., 04.36, 03.01, 00530., 3.21, 467.232334920],
    ['B6.0', 14093., 03.80, 02.81, 00316., 2.99, 279.154727515],
    ['B7.0', 12942., 03.38, 02.65, 00200., 2.82, 176.569574061],
    ['B8.0', 11561., 02.91, 02.44, 00109., 2.61, 95.3190701227],
    ['B9.0', 10351., 02.52, 02.25, 0059.1, 2.39, 52.0850169839]]
    # ['B9.5', 09886., 02.38, 02.17, 00046., 2.32, 40.3107085348]]

bestarsSK1982 = [
    # Schmidt-Kaller1982. Used (and interpolated) by Porter1996, Townsedn2004, 
    # SpType Teff   Mass    Rp    Lum
    ['B0.0', 30105., 17.5, 7.70, 43651.],
    ['B0.5', 27859., 14.6, 6.90, 25703.],
    ['B1.0', 25985., 12.5, 6.30, 16218.],
    ['B1.5', 24347., 10.8, 5.70, 10232.],
    ['B2.0', 22813., 09.6, 5.40, 07079.],
    ['B2.5', 21498., 08.6, 5.00, 04786.],
    ['B3.0', 20222., 07.7, 4.70, 03311.],
    ['B4.0', 18206., 06.4, 4.20, 01737.],
    ['B5.0', 16673., 05.5, 3.80, 01000.],
    ['B6.0', 15302., 04.8, 3.50, 00602.],
    ['B7.0', 14103., 04.2, 3.20, 00363.],
    ['B8.0', 13202., 03.8, 3.00, 00245.],
    ['B9.0', 12246., 03.4, 2.80, 00158.]]

bestarsdJN1987 = [
    # Derived by de Jager & Niewuwenhuijzen 1987 to the main sequence (b=5.) 
    #  lum class IV (b=4.); Used by Cranmer2005
    # Conclusion: 5 and 4 apper do be ZAMS and mid-MS; 3 late MS
    # Conclusion: SpTypes appear to be shifhed by -1.0 here (cooler stars)
    # SpType b-val Teff_V Mass_V Rp_5 Lum_V   Teff_4 Mass_4 Rp_4 Lum_4
    ['B0.0', 1.200, 26841, 13.8, 6.58, 20134., 26911, 15.11, 7.84, 28919.],
    ['B0.5', 1.350, 24944, 11.4, 5.82, 11742., 24809, 12.30, 6.90, 16183.],
    ['B1.0', 1.500, 23213, 9.63, 5.16, 06917., 22915, 10.17, 6.11, 09222.],
    ['B1.5', 1.650, 21629, 8.17, 4.58, 04118., 21204, 08.54, 5.44, 05355.],
    ['B2.0', 1.800, 20178, 7.01, 4.08, 02478., 19655, 07.27, 4.87, 03171.],
    ['B2.5', 1.875, 19498, 6.51, 3.86, 01930., 18935, 06.74, 4.62, 02458.],
    ['B3.0', 1.950, 18846, 6.07, 3.65, 01508., 18250, 06.27, 4.39, 01915.],
    ['B4.0', 2.100, 17621, 5.31, 3.28, 00928., 16972, 05.48, 3.99, 01181.],
    ['B5.0', 2.250, 16493, 4.69, 2.95, 00578., 15810, 04.84, 3.64, 00743.],
    ['B6.0', 2.400, 15452, 4.18, 2.67, 00364., 14749, 04.33, 3.36, 00478.],
    ['B7.0', 2.550, 14491, 3.75, 2.42, 00232., 13780, 03.91, 3.12, 00314.],
    ['B8.0', 2.700, 13601, 3.40, 2.21, 00150., 12893, 03.57, 2.92, 00211.],
    ['B9.0', 2.850, 12778, 3.10, 2.03, 00098., 12080, 03.29, 2.76, 00145.]]    

bestarsdJN1987_3 = [
    # Derived by de Jager & Niewuwenhuijzen 1987 to the main sequence (b=5.) 
    #  lum class IV (b=4.); Used by Cranmer2005
    # Conclusions with Geneva models: class III is still in the main sequence!
    #  (but leaving, ~Achernar)
    # Conclusion: SpTypes appear to be shifhed by -1 step here (cooler stars)
    # SpType b-val Teff_3 Mass_3 Rp_3 Lum_3
    ['B0.0', 1.200, 25030, 14.8, 9.93, 34661.],
    ['B0.5', 1.350, 23009, 12.2, 8.92, 19969.],
    ['B1.0', 1.500, 21198, 10.2, 8.05, 11731.],
    ['B1.5', 1.650, 19570, 8.65, 7.31, 07032.],
    ['B2.0', 1.800, 18105, 7.43, 6.69, 04305.],
    ['B2.5', 1.875, 17427, 6.93, 6.41, 03396.],
    ['B3.0', 1.950, 16782, 6.48, 6.16, 02693.],
    ['B4.0', 2.100, 15586, 5.71, 5.71, 01723.],
    ['B5.0', 2.250, 14502, 5.10, 5.33, 01128.],
    ['B6.0', 2.400, 13519, 4.60, 5.03, 00756.],
    ['B7.0', 2.550, 12624, 4.20, 4.78, 00520.],
    ['B8.0', 2.700, 11809, 3.86, 4.58, 00366.],
    ['B9.0', 2.850, 11065, 3.59, 4.43, 00264.]]   

bestarsBeAtlas = [
    # For ob=1.10
    # SpType Tpole    Teff    Mass   Rp    Lum
    ['B0.0', _np.NaN, _np.NaN, _np.NaN, _np.NaN, _np.NaN],
    ['B0.5', 28905.8, 26765.7, 14.6, 7.50, 31183.26],
    ['B1.0', 26945.8, 24950.9, 12.5, 6.82, 19471.38],
    ['B1.5', 25085.2, 23228.2, 10.8, 6.23, 12204.70],
    ['B2.0', 23629.3, 21879.9, 09.6, 5.80, 08327.67],
    ['B2.5', 22296.1, 20645.4, 08.6, 5.43, 05785.96],
    ['B3.0', 20919.7, 19370.9, 07.7, 5.11, 03971.25],
    ['B4.0', 18739.3, 17351.9, 06.4, 4.62, 02090.08],
    ['B5.0', 17063.8, 15800.5, 05.5, 4.26, 01221.76],
    ['B6.0', 15587.7, 14433.6, 04.8, 4.02, 00757.60],
    ['B7.0', 14300.3, 13241.6, 04.2, 3.72, 00459.55],
    ['B8.0', 13329.9, 12343.0, 03.8, 3.55, 00315.96],
    ['B9.0', 12307.1, 11395.9, 03.4, 3.37, 00206.89]]

bestarsBeAtlas_N = [
    # For ob=1.10
    # SpType Tpole    Teff    Mass   Rp    Lum
    ['B0.0', 28905.8, 26765.7, 14.6, 7.50, 31183.26],
    ['B0.5', 26945.8, 24950.9, 12.5, 6.82, 19471.38],
    ['B1.0', 25085.2, 23228.2, 10.8, 6.23, 12204.70],
    ['B1.5', 23629.3, 21879.9, 09.6, 5.80, 08327.67],
    ['B2.0', 22296.1, 20645.4, 08.6, 5.43, 05785.96],
    ['B2.5', 20919.7, 19370.9, 07.7, 5.11, 03971.25],
    ['B3.0', 18739.3, 17351.9, 06.4, 4.62, 02090.08],
    ['B4.0', 17063.8, 15800.5, 05.5, 4.26, 01221.76],
    ['B5.0', 15587.7, 14433.6, 04.8, 4.02, 00757.60],
    ['B6.0', 14300.3, 13241.6, 04.2, 3.72, 00459.55],
    ['B7.0', 13329.9, 12343.0, 03.8, 3.55, 00315.96],
    ['B8.0', 12307.1, 11395.9, 03.4, 3.37, 00206.89],
    ['B9.0', _np.NaN, _np.NaN, _np.NaN, _np.NaN, _np.NaN]]


# MAIN ###
if __name__ == "__main__":
    pass
