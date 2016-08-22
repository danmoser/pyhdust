# -*- coding:utf-8 -*-

"""PyHdust *spectools* module: spectroscopy tools

Algumas definicoes: para todas as rotinas funcionarem, todos os espectros devem
estar agrupados num mesmo caminho (`path`), em estrutura de
noite/estrelas/espec.

Neste momento, as rotinas somente leem arquivos `*.cal.fits`. Para receber este
sufixo `.cal`, algumas informacoes no header sao necessarias:

    * 'MJD-OBS' ou 'MJD' ou 'JD' ou 'DATE-OBS'
    * 'CRVAL1' + 'CDELT1'

IMPORTANT NOTE: after the version 0.981, the "analline" function returns 
FWHM instead of `depthcent`.

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import os as _os
import numpy as _np
import datetime as _dt
# import time as _time
from glob import glob as _glob
# from itertools import product as _iproduct
import pyhdust.phc as _phc
import pyhdust.jdcal as _jdcal
import pyhdust.input as _inp
import pyhdust as _hdt
from six import string_types as _strtypes
import warnings as _warn


try:    
    import pyfits as _pyfits
    import matplotlib as _mpl
    import matplotlib.pyplot as _plt
    import matplotlib.patches as _mpatches
    from matplotlib.ticker import MaxNLocator as _MaxNLocator
    import matplotlib.gridspec as _gridspec
    import scipy.interpolate as _interpolate 
    from scipy.optimize import curve_fit as _curve_fit
except ImportError:
    _warn.warn('matplotlib, scipy, and/or pyfits module not installed!!!')

try:
    import pyqt_fit.nonparam_regression as _smooth
    from pyqt_fit import npr_methods as _npr_methods
except ImportError:
    _warn.warn('pyqt_fit module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"

_outfold = ''


class Spec(object):

    """Definicao de classe espectro para conter toda a informacao util
    para plots e analises.

    EW in km/s

    Para configurar uma ou mais linhas:

    >>> spdtb = Spec()
    >>> spdtb.lbc == 0
    >>> #significa que vetor wl eh vetor velocidades, e nao comprimento de 
    >>> # onda.
    >>> spdtb.lbc = 6564.
    >>> spdtb2 = Spec()
    >>> spdtb2.lbc = 4863.

    Como usar (hard way):

    >>> spdtb = Spec()
    >>> #read spec `flux` and `wl` for a given `lbc`
    >>> (spdtb.EW, spdtb.EC, spdtb.VR, spdtb.peaksep, spdtb.depthcent,\\ 
    >>> spdtb.F0) = analline(wl, flux, lbc)
    >>> spdtb.MJD = 1
    >>> spdtb.file = file

    And then:

    >>> #to record it to the database:
    >>> spdtb.addspec()

    Para carregar uma tabela anterior, faca:

    >>> spdtb = Spec()
    >>> #(...) read new specs and then join with previous ones
    >>> spdtb.data = _np.vstack((spdtb.data, _np.loadtxt('hdt/datafile.txt')))
    >>> spdtb.metadata = _np.vstack(( spdtb.metadata, \\
    >>> _np.loadtxt('hdt/metafile.txt') ))
    >>> spdtb.updatecount() #to update the counter

    Ou simplesmente (nome de arquivos default):

    >>> spdtb.loaddata()
    """

    def __init__(self, wl=None, flux=None, lbc=None, hwidth=1000., EW=_np.NaN,
        EC=_np.NaN, VR=_np.NaN, peaksep=_np.NaN, depthcent=_np.NaN, F0=_np.NaN,
        dateobs='', MJD=0., datereduc='', file='', gaussfit=False):
        self.wl = wl
        self.flux = flux
        self.lbc = lbc
        self.hwidth = hwidth 
        self.EW = EW
        self.EC = EC
        self.VR = VR
        self.peaksep = peaksep
        self.depthcent = depthcent
        self.F0 = F0
        self.file = file
        self.datereduc = datereduc
        self.dateobs = dateobs
        self.MJD = MJD
        self.count = 0
        self.data = _np.empty(0)
        self.metadata = _np.empty(0)
        self.gaussfit = gaussfit

    def reset(self):
        """Reset the class parameters
        """
        self.wl = None
        self.flux = None
        self.EW = _np.NaN
        self.EC = _np.NaN
        self.VR = _np.NaN
        self.peaksep = _np.NaN
        self.depthcent = _np.NaN
        self.F0 = _np.NaN
        self.file = ''
        self.datereduc = ''
        self.dateobs = ''
        self.MJD = 0.        

    def clear(self):
        """Clear the class parameters
        """
        self.__init__()

    def addspec(self):
        """Record the class parameters into the database
        """
        self.count += 1
        if self.count == 1:
            self.data = _np.array( self.lastinfo() )
            self.metadata = _np.array( self.lastmeta() )
        else:
            self.data = _np.vstack(( self.data, self.lastinfo() ))
            self.metadata = _np.vstack(( self.metadata, self.lastmeta() ))
        # if self.flux != None and self.wl != None and self.lbc != None:
        #    self.savespec()

    def lastinfo(self):
        """Print the current class parameters (last spec)
        """
        return self.MJD, self.EW, self.EC, self.VR, self.peaksep, \
            self.depthcent, self.F0 

    def lastmeta(self):
        """Print the current class parameters (last spec)
        """
        return self.MJD, self.dateobs, self.datereduc, self.file

    def savedata(self, datafile=_outfold + '/datafile.txt',
        metafile=_outfold + '/metafile.txt'):
        """Save current table
        """
        header = ['MJD', 'EW', 'EC', 'VR', 'peaksep', 'depthcent', 'F0']
        _np.savetxt(datafile, self.data, fmt='%12.6f',
        header=(len(header) * '{:>12s}').format(*header))
        _np.savetxt(metafile, self.metadata, fmt='%s', delimiter=',')
        return

    def loaddata(self, datafile=_outfold + '/datafile.txt',
        metafile=_outfold + '/metafile.txt'):
        """Function to load a previous table
        Usage:

        >>> spdtb = Spec()
        >>> spdtb.loaddata()
        """
        self.data = _np.loadtxt(datafile)
        if _os.path.exists(metafile):
            self.metadata = _np.loadtxt(metafile, dtype='str', delimiter=',')
        self.updatecount()
        return

    def updatecount(self, num=0):
        if num > 0:
            self.count = num
        else:
            self.count = len(self.data)
        return

    def loadspec(self, file):
        """Load a fits file (parameters `wl`, `flux`, `MJD`, `dateobs`,
        `datareduc` and `file`).

        Currently, only compatible for standard fits.
        """
        if file.find('.fit') == -1:
            _warn.warn("# ERROR! `loadspec` unrecognized format!")
            return 
        (self.wl, self.flux, self.MJD, self.dateobs, self.datereduc,
            self.file) = loadfits(file)
        (self.EW, self.EC, self.VR, self.peaksep, self.depthcent, self.F0) = \
            analline(self.wl, self.flux, self.lbc, hwidth=self.hwidth, 
            verb=False, gaussfit=self.gaussfit)
        return

    def plotspec(self, outname=''):
        """Export current spec into a PNG file.
        """
        if self.wl is None or self.flux is None:
            _warn.warn('wrong Spec() parameters! {0}'.format(self.file))
            return
        if outname == '':
            path, file = _phc.trimpathname(self.file)
            outname = _phc.rmext(file)
        # Normalization:
        flux = linfit(self.wl, self.flux)
        wl = self.wl
        fig = _plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(wl, flux)
        ax.set_ylabel('norm. flux')
        ax.set_xlabel('wavelength (arb. units)')
        ax.set_title(outname)
        _plt.savefig('{0}/{1:.2f}_{2}.png'.format(_outfold, self.MJD, outname))
        if self.lbc > 0:
            vels = (self.wl - self.lbc) / self.lbc * _phc.c.cgs * 1e-5
            idx = _np.where(_np.abs(vels) <= self.hwidth)
            flux = linfit(vels[idx], flux[idx])
            vels = vels[idx]
            _plt.clf()
            ax = fig.add_subplot(111)
            ax.plot(vels, flux)
            ax.set_ylabel('norm. flux')
            ax.set_xlabel('vel. (km/s)')
            ax.set_title('{0:.2f} {1} {2:.2f}'.format(self.MJD, outname, 
                self.lbc))
            _plt.savefig('{0}/{1:.2f}_{2}_{3:.2f}.png'.format(_outfold, 
                self.MJD, outname, self.lbc))
        _plt.close()
        return


def shiftfits(fitsfile, newsh=None):
    """ Update FITS spec header for a given shift value. """
    imfits = _pyfits.open(fitsfile, mode='update')
    if 'WLSHIFT' in imfits[0].header:
        print('# WLSHIFT = {0} for {1}'.format(imfits[0].header['WLSHIFT'],
        _phc.trimpathname(fitsfile)[1]))
    else:
        print('# No WLSHIFT available for {0}'.format(
            _phc.trimpathname(fitsfile)[1]))
    if newsh is None:
        newsh = _phc.user_input('Type the new shift: ')
    if newsh != '':
        imfits[0].header['WLSHIFT'] = float(newsh)
        imfits.close()
    return 


def checkshiftfits(fitslist, lbc=6562.8):
    """ Do *shiftfits* sistematically 

    INPUT: list of files

    OUTPUT: fits files header updated with WLSHIFT.
    """
    fig, ax = _plt.subplots()
    for f in fitslist:
        print(f)
        data = loadfits(f)
        vel, flx = lineProf(data[0], data[1], lbc=lbc)
        good = False
        imfits = _pyfits.open(f)
        if 'WLSHIFT' in imfits[0].header:
            shift0 = float(imfits[0].header['WLSHIFT'])
        else:
            shift0 = 0.
        shift = 0
        while not good:
            ax.plot([0, 0], [0.7, 1.2], ls='--', color='gray')
            veli = vel + shift*3e5/lbc
            ax.plot(veli, flx)
            _plt.show()
            _plt.draw()
            ri = _phc.user_input('\n# Is it good?(y/other): ')
            if ri != 'y':
                try:
                    shift = float(_phc.user_input('Type shift: '))
                except:
                    shift = 0.
            else:
                good = True
            ax.cla()
        if shift != 0:
            shiftfits(f, newsh=shift+shift0)
    _plt.close(fig)
    return


def loadfits(fitsfile):
    """load FITS spec

    Out: wl, flux, MJD, dateobs, datereduc, fitsfile
    """
    imfits = _pyfits.open(fitsfile)
    flux = imfits[0].data
    wl = _np.arange(len(flux)) * imfits[0].header['CDELT1'] +\
        imfits[0].header['CRVAL1']
    (MJD, dateobs, datereduc) = (0., '', '')
    if 'MJD-OBS' in imfits[0].header:
        MJD = float(imfits[0].header['MJD-OBS'])
    elif 'MJD' in imfits[0].header:
        MJD = float(imfits[0].header['MJD'])
    elif 'JD' in imfits[0].header:
        MJD = float(imfits[0].header['JD']) - 2400000.5
    elif 'DATE-OBS' in imfits[0].header:
        dtobs = imfits[0].header['DATE-OBS']
        dtobs, tobs = check_dtobs(dtobs)
        MJD = _jdcal.gcal2jd(*dtobs)[1] + tobs
    elif 'FRAME' in imfits[0].header:
        dtobs = imfits[0].header['FRAME']
        dtobs, tobs = check_dtobs(dtobs)
        MJD = _jdcal.gcal2jd(*dtobs)[1] + tobs
    else:
        MJD = _jdcal.MJD_JD2000
        _warn.warn('No DATE-OBS information is available! {0}\nAssuming '
            'MJD_JD2000'.format(fitsfile))
    if 'DATE-OBS' in imfits[0].header:
        dateobs = imfits[0].header['DATE-OBS']
    elif 'FRAME' in imfits[0].header:
        dateobs = imfits[0].header['FRAME']
    if 'IRAF-TLM' in imfits[0].header:
        datereduc = imfits[0].header['IRAF-TLM']
    elif 'DATE' in imfits[0].header:
        datereduc = imfits[0].header['DATE']
    if 'WLSHIFT' in imfits[0].header:
        shift = float(imfits[0].header['WLSHIFT'])
        wl += shift
    imfits.close()
    return wl, flux, MJD, dateobs, datereduc, fitsfile


def vac2air(wl):
    """The IAU standard for conversion from air to vacuum wavelengths is given
    in Morton (1991, ApJS, 77, 119). For vacuum wavelengths (VAC) in Angstroms,
    convert to air wavelength (AIR) via:
    AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4 )
    """
    return wl / (1.0 + 2.735182E-4 + 131.4182 / wl**2 + 2.76249E8 / wl**4 )


def air2vac(wl):
    """The IAU standard for conversion from air to vacuum wavelengths is given
    in Morton (1991, ApJS, 77, 119). For vacuum wavelengths (VAC) in Angstroms,
    convert to air wavelength (AIR) via:
    AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4 )

    Fitting the inverse curve:
    VAC = AIR / (1.0 - 2.73443407E-4 - 1.31275255E2 / AIR^2 - 2.75708212E8 / 
    AIR^4 )
    """
    return wl / (1.0 - 2.73443407e-04 - 1.31275255e+02 / wl**2 - 
        2.75708212e+08 / wl**4)


def vel2wl(vel, lbc):
    """ Vel. to wavelength. Vel must be in km/s and output is in `lbc` units. 
    """
    wl = (vel / _phc.c.cgs * 1e5 + 1) * lbc
    return wl


def wl2vel(wl, lbc):
    """ Wavelength to vel., in km/s. `wl` and `lbc` units must be the same. """
    vels = (wl - lbc) / lbc * _phc.c.cgs * 1e-5
    return vels


def hydrogenlinewl(ni, nf):
    """Generate H line transitions wavelengths in meters for VACUUM

    Rydberg constant `R` was manually adjusted to fit Halpha and Hbeta lines.
    """
    return (10967850. * (1. / nf**2 - 1. / ni**2))**-1.


def calcres_R(hwidth=1350, nbins=108):
    """
    (h)Width in km/s.
    *WARNING*: `width` in HDUST input is only half.
    To HDUST effective R, multiple the input width by 2 he_re.

    # R = lbd/Dlbd = _phc.c/Dv = _phc.c*nbins/width 
    # nbins = R*width/_phc.c
    """
    return round(_phc.c.cgs * nbins / hwidth / 1e5)


def calcres_nbins(R=12000, hwidth=1350):
    """
    (h)Width in km/s.
    *WARNING*: `width` in HDUST input is only half.
    To HDUST effective R, multiple the input width by 2 he_re.

    # R = lbd/Dlbd = _phc.c/Dv = _phc.c*nbins/width 
    # nbins = R*width/_phc.c
    """
    return round(R * hwidth * 1e5 / _phc.c.cgs)


def lineProf(x, flx, lbc, flxerr=_np.empty(0), hwidth=1000., ssize=0.05):
    '''
    lineProf() - retorna um array (flx) normalizado e um array x em 
    VELOCIDADES. `lbc` deve fornecido em mesma unidade de x para conversão 
    lambda -> vel. Se vetor x jah esta em vel., usar funcao linfit().

    x eh importante, pois y pode ser nao igualmente amostrado.
    x e y devem estar em ordem crescente.

    ssize = % do tamanho de y; numero de pontos usados nas extremidades
    para a media do contínuo. 'ssize' de .5 à 0 (exclusive).

    OUTPUT: vel (array), flx (array)
    '''    
    x = (x - lbc) / lbc * _phc.c.cgs * 1e-5  # km/s
    idx = _np.where(_np.abs(x) <= 1.001 * hwidth)
    if len(flxerr) == 0:
        flux = linfit(x[idx], flx[idx], ssize=ssize)  # yerr=flxerr,
        if len(x[idx]) == 0:
            _warn.warn('Wrong `lbc` in the lineProf function')
        return x[idx], flux
    else:
        flux, flxerr = linfit(x[idx], flx[idx], yerr=flxerr[idx], ssize=ssize)
        if len(x[idx]) == 0:
            _warn.warn('Wrong `lbc` in the lineProf function')
        return x[idx], flux, flxerr


def linfit(x, y, ssize=0.05, yerr=_np.empty(0)):
    '''
    linfit() - retorna um array (y) normalizado, em posicoes de x

    x eh importante, pois y pode ser nao igualmente amostrado.
    x e y devem estar em ordem crescente.

    ssize = % do tamanho de y; numero de pontos usados nas extremidades
    para a media do contínuo. 'ssize' de .5 à 0 (exclusive).

    OUTPUT: y, yerr (if given)

    .. code:: python

        #Example:
        import numpy as np
        import matplotlib.pyplot as plt
        import pyhdust.phc as phc
        import pyhdust.spectools as spt

        wv = _np.linspace(6500, 6600, 101)
        flx = (np.arange(101)[::-1])/100.+1+phc.normgauss(4, x=wv, 
        xc=6562.79)*5

        plt.plot(wv, flx)
        normflx = linfit(wv, flx)
        plt.plot(wv, normflx, ls='--')

        plt.xlabel(r'$\lambda$ ($\AA$)')
        plt.ylabel('Flux (arb. unit)')

    .. image:: _static/spt_linfit.png
        :align: center
        :width: 500
    '''
    if ssize < 0 or ssize > .5:
        _warn.warn('Invalid ssize value...', stacklevel=2)
        ssize = 0
    ssize = int(ssize * len(y))
    if ssize == 0:
        ssize = 1
    medx0, medx1 = _np.average(x[:ssize]), _np.average(x[-ssize:])
    if ssize > 9:
        medy0, medy1 = _np.median(y[:ssize]), _np.median(y[-ssize:])
    else:
        medy0, medy1 = _np.average(y[:ssize]), _np.average(y[-ssize:])
    new_y = medy0 + (medy1 - medy0) * (x - medx0) / (medx1 - medx0)
    idx = _np.where(new_y != 0)
    y[idx] = y[idx] / new_y[idx]
    if len(yerr) == 0.:
        return y
    else:
        yerr = yerr / _np.average(new_y)
        return y, yerr


def EWcalc(vels, flux, vw=1000):
    """
    Supoe que o fluxo jah estah normalizado, e vetores ordenad_os.

    Devolve o valor EW, alem dos vetores cortados em vw.
    """
    idx = _np.where(_np.abs(vels) <= vw)
    outvels = vels[idx]
    normflux = flux[idx]
    ew = 0.
    if len(outvels) < 3:
        # normflux = _np.ones(len(outvels))
        return ew
    for i in range(len(outvels) - 1):
        dl = outvels[i + 1] - outvels[i]
        ew += (1. - (normflux[i + 1] + normflux[i]) / 2.) * dl
    return ew


def ECcalc(vels, flux, ssize=.05, gaussfit=False, doublegf=True):
    """
    Supoe que o fluxo jah estah normalizado, e vetores ordenados.

    Calcula o topo da emissao da linha, e retorna em que velocidade ela
    ocorre.
    """
    vels = _np.array(vels)
    flux = _np.array(flux)
    # if lncore > 0:
    #     idx = _np.where(_np.abs(vels) < lncore)
    #     vels = vels[idx]
    #     flux = flux[idx]
    if not gaussfit:
        idx = _np.where(_np.max(flux) == flux)
        if flux[idx][0] < 1:
            return _np.NaN, 0.
        if len(idx[0]) > 1:
            idx = idx[0][0]
        return flux[idx][0], vels[idx][0]
    else:
        # check if there is a peak
        ssize = int(ssize * len(vels))
        if ssize == 0:
            ssize = 1
        contmax = _np.max(_np.append(flux[:ssize], flux[-ssize:]))
        fluxmax = _np.max(flux)
        if fluxmax < 1.01 * contmax:
            return _np.NaN, vels[-1]

        # Define model function to be used to fit to the data above
        def gauss(x, *p):
            A, mu, sigma = p
            return A * _np.exp(-(x - mu)**2 / (2. * sigma**2)) + 1
        # 
        ivc = _np.abs(vels - 0).argmin()
        if doublegf:
            i0 = _np.abs(flux[:ivc] - _np.max(flux[:ivc])).argmin()
            i1 = _np.abs(flux[ivc:] - _np.max(flux[ivc:])).argmin() + ivc
            try: 
                p0 = [1., vels[i0], 40.]
                coeff0, tmp = _curve_fit(gauss, vels[:ivc], flux[:ivc], p0=p0)
                p1 = [1., vels[i1], 40.]
                coeff1, tmp = _curve_fit(gauss, vels[ivc:], flux[ivc:], p0=p1)
                EC = _np.max([coeff0[0] + 1., coeff1[0] + 1.])
                vel = _np.abs(coeff0[1] / 2) + _np.abs(coeff1[1] / 2)
                return EC, vel
            except:
                return 1., vels[-1]
        else:
            try:
                p0 = [1., 0, 40.]
                coeff0, tmp = _curve_fit(gauss, vels[:ivc], flux[:ivc], p0=p0)
                EC = coeff0[0] + 1.
                return EC, coeff0[1]
            except:
                return 1., vels[-1]


def VRcalc(vels, flux, vw=1000, gaussfit=False, ssize=0.05):
    """
    Calcula o PICO para os dois lados (azul/vermelho) da linha, ajustando
    a velocidade de repouso (TBD). 
    """
    # calcula e aplica correcao de vel. repousp
    vc = 0.
    vels += vc
    # corta em vw, e faz o teste de tamanho
    if len(vels) < 5:
        vw = 0
    if vw > 0:
        idx = _np.where(_np.abs(vels) <= vw)
        outvels = vels[idx]
        normflux = flux[idx]
    else:
        ew0, ew1 = 0.
        return ew0, ew1, vc
    #
    ivc = _np.abs(outvels - 0).argmin()
    if not gaussfit:
        V = _np.max(normflux[:ivc])
        R = _np.max(normflux[ivc:])
    else:
        # check if there is a peak
        ssize = int(ssize * len(vels))
        if ssize == 0:
            ssize = 1
        contmax = _np.max(_np.append(flux[:ssize], flux[-ssize:]))
        fluxmax = _np.max(flux)
        if fluxmax < 1.01 * contmax:
            # print('# Bad profile!')
            return 0, 0, vc

        # Define model function to be used to fit to the data above
        def gauss(x, *p):
            A, mu, sigma = p
            return A * _np.exp(-(x - mu)**2 / (2. * sigma**2)) + 1.
        # 
        ivc = _np.abs(vels - 0).argmin()
        i0 = _np.abs(flux[:ivc] - _np.max(flux[:ivc])).argmin()
        i1 = _np.abs(flux[ivc:] - _np.max(flux[ivc:])).argmin() + ivc
        try:
            p0 = [1., vels[i0], 40.]
            coeff0, tmp = _curve_fit(gauss, vels[:ivc], flux[:ivc], p0=p0)
            p1 = [1., vels[i1], 40.]
            coeff1, tmp = _curve_fit(gauss, vels[ivc:], flux[ivc:], p0=p1)
            V = coeff0[0] + 1.
            R = coeff1[0] + 1.
        except:
            return 1., 1., vc
    return V, R, vc


def PScalc(vels, flux, vc=0., ssize=.05, gaussfit=False):
    """
    Calcula peak_separation

    `doublegaussfit` = True, do it before and after zero velocity. False, use
    maximum (default).
    """
    # check if there is a peak
    ssize = int(ssize * len(vels))
    if ssize == 0:
        ssize = 1
    contmax = _np.max(_np.append(flux[:ssize], flux[-ssize:]))
    fluxmax = _np.max(flux)
    if fluxmax < 1.01 * contmax:
        return _np.NaN, _np.NaN
    vels += vc
    ivc = _np.abs(vels - 0).argmin()
    i0 = _np.abs(flux[:ivc] - _np.max(flux[:ivc])).argmin()
    i1 = _np.abs(flux[ivc:] - _np.max(flux[ivc:])).argmin() + ivc
    if not gaussfit:
        return vels[i0], vels[i1]
    else:
        # Define model function to be used to fit to the data above
        def gauss(x, *p):
            A, mu, sigma = p
            return A * _np.exp(-(x - mu)**2 / (2. * sigma**2)) + 1.
        # 
        try:
            p0 = [1., vels[i0], 20.]
            coeff0, tmp = _curve_fit(gauss, vels[:ivc], flux[:ivc], p0=p0)
            p1 = [1., vels[i1], 20.]
            coeff1, tmp = _curve_fit(gauss, vels[ivc:], flux[ivc:], p0=p1)
            return coeff0[1], coeff1[1]
        except:
            # print vels[i0], flux[i0], vels[i1], flux[i1]
            return 0, 0


def FWHM(vels, flux, halfmax, vmax=350., flxincr=.01):
    """ Calc. FWHM (Full-Width at Half Maximum) based on the value of the 
    Half Maximum

    TODO: Gaussfit"""
    vels = _np.array(vels)
    flux = _np.array(flux)
    # remove vels bigger than maxvel
    idx = _np.where(_np.abs(vels) < vmax)
    vels = vels[idx]
    flux = flux[idx]
    difflx = _np.abs(flux - halfmax)
    # remove diff bigger than hmf*halfmax
    i = 0
    idx = _np.where(difflx < halfmax * flxincr*i)
    while len(vels[idx]) < 2:
        i += 1
        idx = _np.where(difflx < halfmax * flxincr*i)
    vels = vels[idx]
    difflx = difflx[idx]
    #
    # difvels: ordered vels based on the flux difference
    # idx = _np.argsort(difflx)
    # difvels = vels[idx][:4]
    #
    # difvels: ordered vels closest to the 0 vel.
    idx = _np.argsort(_np.abs(vels))
    difvels = vels[idx][:2]
    return _np.sum(_np.abs(difvels))


def DCcalc(vels, flux, vmax=None, vc=0., ssize=0.05):
    """
    Calculo, na presenca de emissao, da profundidade do reverso central.

    TODO: gauss fit
    """
    vels += vc
    ivc = _np.abs(vels - 0).argmin()
    # check if there is a peak
    ssize = int(ssize * len(vels))
    if ssize == 0:
        ssize = 1
    contmax = _np.max(_np.append(flux[:ssize], flux[-ssize:]))
    fluxmax = _np.max(flux)
    if fluxmax < 1.01 * contmax:
        return flux[ivc], flux[ivc]
    # if a vmax is not given...
    if not isinstance(vmax, (int, long, float)):
        vmax = _np.abs(flux - _np.max(flux)).argmin()
        vmax = vels[vmax]
    ivmax = _np.abs(vels - vmax).argmin()
    return flux[ivmax], flux[ivc]


def analline(lbd, flux, lbdc, hwidth=1000, verb=True, gaussfit=False,
    doublegf=True):
    """
    Return the analysis of a line.

    Both lbd and flux need to be ordered (a normalization IS FORCED).
    lbd,lbdc must have the same unit, and width in km/s is required.
    The line will be cutted so that the total DeltaLambda will be 2*width

    if `lbdc` <= 0, lbd array is assumed to be a velocity array (in km/s)!

    | EXAMPLE: Using sed2data. lbc = 0.6565 (halpha), obs = 1 (width==1000)
    |     analline(lbd=sed2data[obs,:,2], flux=sed2data[obs,:,3], lbc=lbc)

    The EW is the equivalent width in km/s, 
    EC is the Emission over Continuum ratio, 
    VR ratio, 
    peaksep in km/s, 
    FWHM is the Full-Width at Half Maximum (emission as maximum)
    F0 is the depth of rest wavelength normalized to the continuum 

    OUTPUT: EW, EC, VR, peaksep, FWHM, F0 
    """
    if lbdc > 0:
        vels = (lbd - lbdc) / lbdc * _phc.c.cgs * 1e-5
    else:
        vels = lbd
    # check if the file have the desired info.
    if vels[0] > -hwidth * .95 or vels[-1] < hwidth * .95:
        if verb:
            _warn.warn('spec out of range (wavelength)! Check hwidth!')
        return _np.NaN, _np.NaN, _np.NaN, _np.NaN, _np.NaN, _np.NaN

    idx = _np.where(_np.abs(vels) <= hwidth)
    vels = vels[idx]
    flux = flux[idx]
    # Normalization:
    flux = linfit(vels, flux)
    # Output:
    EW = EWcalc(vels, flux, vw=hwidth)
    EC, velEC = ECcalc(vels, flux, gaussfit=gaussfit, doublegf=doublegf)
    ew0, ew1, vc = VRcalc(vels, flux, vw=hwidth, gaussfit=gaussfit)
    if ew1 == 0 or EC is _np.NaN:
        VR = 1
    else:
        VR = ew0 / ew1
    if EC is _np.NaN:
        peaksep = _np.NaN
    else:
        vel0, vel1 = PScalc(vels, flux, gaussfit=gaussfit)
        peaksep = vel1 - vel0
    if peaksep is _np.NaN:
        EC = peaksep
        VR = peaksep
    EC2, F0 = DCcalc(vels, flux, vmax=velEC)
    # depthcent = EC2 - F0
    if EC2 < 1:
        EC2 = 1.
    fwhm = FWHM(vels, flux, (EC2 + F0) / 2., vmax=_np.abs(velEC))
    return EW, EC, VR, peaksep, fwhm, F0


def kurlog(file=None, output=None):
    """ Generate a list of teff and logg present in a Kurucz file.

    If output is not specified, it is saved as `file`+.log """
    if file is None:
        file = _os.path.join(_hdt.hdtpath(), 'refs', 'fp00k0.pck')
    teffs = []
    loggs = []
    fp = open(file)
    for i, line in enumerate(fp):
        if line.find('TEFF') > -1:
            teffs += [float(line.split()[1])]
            loggs += [float(line.split()[3])]
    fp.close()
    return teffs, loggs


def kuruczflux(teff, logg, wavrange=None):
    """ Return fluxes from a Kurucz model.

    Fluxes are in ergs/cm**2/s/hz/ster and wavelength in nm (wavrange must be 
    in nm).

    OUTPUT: wv, flux, info"""
    kurfile = _os.path.join(_hdt.hdtpath(), 'refs', 'fp00k0.pck')
    kurwvlines = (174 - 22)
    kurflxcol = 10
    # wave
    read = _phc.readrange(kurfile, 22, 22 + kurwvlines)
    wave = _np.array([val for line in read for val in line.split()], 
        dtype=float)
    # choose best
    bestT = _np.inf
    bestg = _np.inf
    fp = open(kurfile)
    for i, line in enumerate(fp):
        if line.find('TEFF') > -1:
            readT = float(line.split()[1])
            if _np.abs(readT - teff) <= _np.abs(bestT - teff):
                bestT = readT
                readg = float(line.split()[3])
                if _np.abs(readg - logg) <= _np.abs(bestg - logg):
                    bestg = readg
                    i0 = i + 1
    fp.close()
    best = [bestT, bestg]
    # read best flux
    read = _phc.readrange(kurfile, i0, i0 + kurwvlines)
    flux = _np.array([val for line in read for val in 
        (line[i:i + kurflxcol] for i in range(0, len(line) - 1, kurflxcol))], 
        dtype=float)
    # cut range
    if wavrange is None:
        return wave, flux, best
    else:
        idx = _np.where((wave > wavrange[0]) & (wave < wavrange[-1]))
        return wave[idx], flux[idx], best    


def plot_all(fs2list, obsl=None, fmt=['png'], out=None, lbc=.6564, 
    hwidth=1000., solidfiles=True, xax=0, philist=[0], figname=None, 
    nolabels=False, obsidx=False):
    r""" plot_all-like routine

    ``obsl`` list, in degrees. It will find the closest values. It the find 
    :math:`\Delta\theta > 3^\circ`, a warning message is displayed. The ``obs``
    index can be used if ``obsidx = True``.

    ``solinefiles`` keep solid lines for files (changes only colors), and 
    change line shapes between observers. If ``False``, do the opposite.
    """
    if isinstance(fs2list, _strtypes):
        fs2list = [fs2list]
    if not isinstance(obsl, list) and obsl is not None:
        _warn.warn('Wrong `obsl` format (None or list)', stacklevel=2)
        return

    fig = _plt.figure(figsize=(9, 9))
    lins, cols = (3, 2)
    gs = _gridspec.GridSpec(lins, cols)
    gs.update(hspace=0.25)

    axt = _plt.subplot(gs[0, 1])
    ax0 = _plt.subplot(gs[1, 0])
    ax1 = _plt.subplot(gs[1, 1])
    ax2 = _plt.subplot(gs[2, 0])
    ax3 = _plt.subplot(gs[2, 1])

    xtitle = 'radial scale'
    for f in fs2list:
        m = _inp.HdustMod(f)
        tfile = _os.path.join(m.proj, m.modn, m.modn+m.suf+'*avg.temp')
        tfile = _glob(tfile)
        if len(tfile) > 0:
            npt, rplus, lev = (0, 0, 0)
            tfile.sort()
            tfile = tfile[-1]
            ncr, ncmu, ncphi, nLTE, nNLTE, Rstar, Ra, beta, data, pcr, pcmu, \
                pcphi = _hdt.readtemp(tfile)
            for phiidx in range(0, len(philist)):
                icphi = philist[phiidx]
                x = data[0, :, 0, icphi]
                if (xax == 0):
                    x = _np.log10(x / Rstar - 1.)
                    xtitle = r'$\log_{10}(r/R_*-1)$'
                elif (xax == 1):
                    x = x / Rstar
                    xtitle = r'$r/R_*$'
                elif (xax == 2):
                    x = 1. - Rstar / x
                    xtitle = r'$1-R_*/r$'
                y = data[3 + lev, :, ncmu / 2 + npt + rplus, icphi]
                y = y / 1000.
                axt.plot(x, y, 'o-')

        fs2d = _hdt.readfullsed2(f)

        iobs = range(len(fs2d))
        if obsl is not None:
            if not obsidx:
                iobs = [_phc.find_nearest(_np.arccos(fs2d[:, 0, 0])*180/_np.pi, 
                    ob, idx=True) for ob in obsl]
            else:
                iobs = obsl

        for ob in iobs:
            obfmt = r'{:.1f}$^\circ$, {:.1f}$^\circ$'.format(_np.arccos(
                fs2d[ob, 0, 0])*180/_np.pi, _np.arccos(fs2d[ob, 0, 1]))
            if solidfiles:
                pdict = {'color': _phc.cycles(fs2list.index(f)), 
                    'dashes': _phc.dashes(iobs.index(ob))}
            else:
                pdict = {'dashes': _phc.dashes(fs2list.index(f)), 
                    'color': _phc.cycles(iobs.index(ob))}

            ax0.plot(fs2d[ob, :, 2], fs2d[ob, :, 3], 
                label=_os.path.basename(f), **pdict)
            ax1.plot(fs2d[ob, :, 2], fs2d[ob, :, 3], 
                label=obfmt, **pdict)
            ax2.plot(fs2d[ob, :, 2], fs2d[ob, :, 7]*100, **pdict)
            ax3.plot(*lineProf(fs2d[ob, :, 2], fs2d[ob, :, 3], lbc=lbc, 
                hwidth=hwidth), **pdict)

    axt.set_xlabel(xtitle, labelpad=1)
    axt.set_ylabel(r'Temperature (10$^3$ K)')
    ax0.set_xlim([.37, 1.])
    ax0.autoscale(axis='y', tight=True)
    ax0.set_yscale('log')
    ax0.set_xlabel(r'$\mu$m')
    ax0.set_ylabel(r'$\lambda F_\lambda/F$')
    ax1.set_xlim([1., 100.])
    ax1.autoscale(axis='y', tight=True)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$\mu$m', labelpad=1)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.set_ticks_position('both')
    ax1.set_ylabel(r'$\lambda F_\lambda/F$')
    ax2.set_xlim([.37, .9])
    ax2.autoscale(axis='y', tight=True)
    ax2.set_xlabel(r'$\mu$m')
    ax2.set_ylabel('P (%)')
    ax3.set_xlim([-hwidth, hwidth])
    ax3.set_xlabel(r'km/s')
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.set_ticks_position('both')
    ax3.set_ylabel('Normalized Flux')

    if not nolabels:
        ax1.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=9,
        labelspacing=0.05)
    if len(fs2list) > 1 and not nolabels:
        ax0.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=8,
            labelspacing=0.05)

    _phc.savefig(fig, fmt=fmt, figname=figname)  # figname='outname')
    return


def splitKurucz(filen, path=None):
    """
    Split atmospheric Kurucz file (e.g., 'ap00k0.dat') into individual models.

    INPUT: file, path (strings)

    OUTPUT: *files written
    """
    if path is None:
        path = _os.getcwd()
    allk = _np.loadtxt(filen, dtype=str, delimiter='\n')
    src = _os.path.splitext(_os.path.split(filen)[1])[0]
    if not _os.path.exists(src):
        _os.mkdir(src)
    src = _os.path.join(src, src)

    for i in range(0, len(allk) - 1):
        if 'EFF' in allk[i]:
            iref = i
            teff = int(allk[i].split()[1][:-1])
            logg = float(allk[i].split()[3][:-3])
        elif 'DECK6 72' in allk[i]:
            allk[i] = allk[i].replace('DECK6 72', 'DECK6 71')
        elif 'EFF' in allk[i + 1]:
            _np.savetxt(src+'tef%05dg%.1f.dat' % (teff, logg), 
                allk[iref:i + 1], fmt='%s')

    _np.savetxt(src+'tef%05dg%.1f.dat' % (teff, logg), allk[iref:], fmt='%s')
    return


def writeFits(flx, lbd, extrahead=None, savename=None, quiet=False, path=None,
    lbdc=None, externhd=None):
    """ Write a 1D spectra FITS.

    | INPUT: flux array, lbd array, extrahead flag+info, save name.
    | - lbd array: if len(lbd)==2: lbd = [CRVAL1, CDELT1]
    |              else: CDELT1 = (lbd[-1]-lbd[0])/(len(lbd)-1)
    |                    CRVAL1 = lbd[0]
    |   WARNING: lbd must be in ANGSTROMS (FITS default). It can also be 
    |   velocities. In this case, it must be in km/s and lbdc is given in 
    | ANGSTROM.
    | - extrahead: matrix (n,2). Example: [['OBJECT','Achernar'], ['COMMENT',
    | 'value']]

    `externhd` = copy the header from an external file.

    OUTPUT: write FITS file.
    """
    if path is None or path == '':
        path = _os.getcwd()
    if path[-1] != ['/']:
        path += '/'
    if lbdc is not None:
        lbd = (lbd / _phc.c.cgs * 1e5 + 1) * lbdc
    hdu = _pyfits.PrimaryHDU(flx)
    hdulist = _pyfits.HDUList([hdu])
    if externhd is not None:
        extf = _pyfits.open(externhd)
        hdulist[0].header = extf[0].header
        hdulist[0].header['BZERO'] = 0.
    hdulist[0].header['CRVAL1'] = lbd[0]
    if len(lbd) == 2:
        hdulist[0].header['CDELT1'] = lbd[1]
    else:
        hdulist[0].header['CDELT1'] = (lbd[-1] - lbd[0]) / (len(lbd) - 1)
    if extrahead is not None:
        for e in extrahead:
            hdulist[0].header[e[0]] = e[1]
    if savename is None:
        savename = 'spec_{0}'.format(_phc.dtflag())
    if savename.find('.fit') == -1:
        savename += '.fits'
    hdu.writeto(path + savename, clobber=True)
    if not quiet:
        print('# FITS file {0}{1} saved!'.format(path, savename))
    return


def averagespecs(speclist, n=999, path='', objname='OBJECT'):
    """ Average specs taken in the same MJD, in groups of approx. `n` 
    elements.

    OUTPUT: Files written. """
    if len(path) > 0 and path[-1] != '/':
        path += '/'
    speclist = _np.array(speclist)
    obsdates = []
    for sp in speclist:
        data = loadfits(sp)
        obsdates.append(data[2])
    obsdates = _np.array(obsdates)
    # Sorting things
    idx = _np.argsort(obsdates)
    speclist = speclist[idx]
    obsdates = obsdates[idx]
    # Same day
    iMJD = []
    for m in obsdates:
        iMJD.append(divmod(m, 1)[0])
    idxMJD = _np.unique(iMJD)
    # Do the avgs based on the MJD
    for i in idxMJD:
        idx = _np.where(iMJD == i)
        N = len(speclist[idx])
        for j in _phc.splitequal(N/n, N):
            fidx = speclist[idx][j[0]:j[1]]
            data = loadfits(fidx[0])
            wl = data[0]
            newdate = _np.average( obsdates[idx][j[0]:j[1]] )
            MJD = int(divmod(newdate, 1)[0])
            MJDfrac = int(round( divmod(newdate, 1)[1]*10000 ))
            fluxes = _np.zeros(len(wl))
            for f in fidx:
                data = loadfits(f)
                fluxes += _np.interp(wl, data[0], data[1])
            flx = fluxes/len(fidx)
            outname = 'alpEri_PUCHEROS_VIS_{0}_{1:04d}_avg.fits'.format(MJD, 
                MJDfrac)
            writeFits( flx, wl, savename=outname, path=path, extrahead=[ 
                ['OBJECT', objname], ['Comment', 'Averaged from {0} spectra'.
                format(len(fidx))], ['MJD-OBS', newdate] ] )
    return


def cardelli(lbd, flux, ebv=0., Rv=3.1):
    """
    Milky Way Extinction law from Cardelli et al. 1989

    `lbd` must be in microns.

    OUTPUT: Corrected flux.
    """
    x = 1. / _np.array(lbd)  # CCM x is 1/microns
    a, b = _np.ndarray(x.shape, x.dtype), _np.ndarray(x.shape, x.dtype)

    if any((x < 0.3) | (10 < x)):
        raise ValueError('Some wavelengths outside CCM 89 extinction curve ' + 
            'range')

    irs = (0.3 <= x) & (x <= 1.1)
    opts = (1.1 <= x) & (x <= 3.3)
    nuv1s = (3.3 <= x) & (x <= 5.9)
    nuv2s = (5.9 <= x) & (x <= 8)
    fuvs = (8 <= x) & (x <= 10)

    # CCM Infrared
    a[irs] = .574 * x[irs]**1.61
    b[irs] = -0.527 * x[irs]**1.61

    # CCM NIR/optical
    a[opts] = _np.polyval((.32999, -.7753, .01979, .72085, -.02427, -.50447, 
        .17699, 1), x[opts] - 1.82)
    b[opts] = _np.polyval((-2.09002, 5.3026, -.62251, -5.38434, 1.07233, 
        2.28305, 1.41338, 0), x[opts] - 1.82)

    # CCM NUV
    a[nuv1s] = 1.752 - .316 * x[nuv1s] - 0.104 / ((x[nuv1s] - 4.67)**2 + .341)
    b[nuv1s] = -3.09 + 1.825 * x[nuv1s] + 1.206 / ((x[nuv1s] - 4.62)**2 + .263)

    y = x[nuv2s] - 5.9
    Fa = -.04473 * y**2 - .009779 * y**3
    Fb = -.2130 * y**2 - .1207 * y**3
    a[nuv2s] = 1.752 - .316 * x[nuv2s] - 0.104 / \
        ((x[nuv2s] - 4.67)**2 + .341) + Fa
    b[nuv2s] = -3.09 + 1.825 * x[nuv2s] + \
        1.206 / ((x[nuv2s] - 4.62)**2 + .263) + Fb

    # CCM FUV
    a[fuvs] = _np.polyval((-.070, .137, -.628, -1.073), x[fuvs] - 8)
    b[fuvs] = _np.polyval((.374, -.42, 4.257, 13.67), x[fuvs] - 8)

    AlbAv = a + b / Rv
    return flux * 10**(-AlbAv * Rv * ebv / 2.5)


def fitzpatrick(wave, flux, ebv, Rv=3.1, LMC2=False, AVGLMC=False):
    """
    Deredden a flux vector using the Fitzpatrick (1999) parameterization

    Parameters
    ----------
    wave :   array
             Wavelength in Angstrom
    flux :   array
             Calibrated flux vector, same number of elements as wave.
    ebv  :   float, optional
             Color excess E(B-V). If a positive ebv is supplied,
             then fluxes will be dereddened rather than reddened.
             The default is 3.1.
    AVGLMC : boolean
             If True, then the default fit parameters c1,c2,c3,c4,gamma,x0 
             are set to the average values determined for reddening in the 
             general Large Magellanic Cloud (LMC) field by
             Misselt et al. (1999, ApJ, 515, 128). The default is
             False.
    LMC2 :   boolean
             If True, the fit parameters are set to the values determined
             for the LMC2 field (including 30 Dor) by Misselt et al.
             Note that neither `AVGLMC` nor `LMC2` will alter the default value 
             of Rv, which is poorly known for the LMC.

    Returns
    -------             
    new_flux : array 
               Dereddened flux vector, same units and number of elements
               as input flux.

    Notes
    -----

    .. note::

        This function was ported from the IDL Astronomy User's Library.

        The following five parameters allow the user to customize
        the adopted extinction curve.    For example, see Clayton et al. (2003,
        ApJ, 588, 871) for examples of these parameters in different 
        interstellar environments.

        x0 - Centroid of 2200 A bump in microns (default = 4.596)
        gamma - Width of 2200 A bump in microns (default  =0.99)
        c3 - Strength of the 2200 A bump (default = 3.23)
        c4 - FUV curvature (default = 0.41)
        c2 - Slope of the linear UV extinction component 
            (default = -0.824 + 4.717/R)
        c1 - Intercept of the linear UV extinction component 
            (default = 2.030 - 3.007*c2
    """
    # x = 10000./ wave # Convert to inverse microns
    x = 1. / wave  # microns
    curve = x * 0.

    # Set some standard values:
    x0 = 4.596
    gamma = 0.99
    c3 = 3.23      
    c4 = 0.41    
    c2 = -0.824 + 4.717 / Rv
    c1 = 2.030 - 3.007 * c2

    if LMC2:
        x0 = 4.626
        gamma = 1.05   
        c4 = 0.42   
        c3 = 1.92      
        c2 = 1.31
        c1 = -2.16
    elif AVGLMC:   
        x0 = 4.596  
        gamma = 0.91
        c4 = 0.64  
        c3 = 2.73      
        c2 = 1.11
        c1 = -1.28

    # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function  
    # and R-dependent coefficients
    xcutuv = _np.array([10000.0 / 2700.0])
    xspluv = 10000.0 / _np.array([2700.0, 2600.0])

    iuv = _np.where(x >= xcutuv)[0]
    N_UV = len(iuv)
    iopir = _np.where(x < xcutuv)[0]
    Nopir = len(iopir)
    if (N_UV > 0): 
        xuv = _np.concatenate((xspluv, x[iuv]))
    else: 
        xuv = xspluv

    yuv = c1 + c2 * xuv
    yuv = yuv + c3 * xuv**2 / ((xuv**2 - x0**2)**2 + (xuv * gamma)**2)
    yuv = yuv + c4 * (0.5392 * (_np.maximum(xuv, 5.9) - 5.9)**2 + 0.05644 * (
        _np.maximum(xuv, 5.9) - 5.9)**3)
    yuv = yuv + Rv
    yspluv = yuv[0:2]  # save spline points

    if (N_UV > 0): 
        curve[iuv] = yuv[2::]  # remove spline points

    # Compute optical portion of A(lambda)/E(B-V) curve
    # using cubic spline anchored in UV, optical, and IR
    xsplopir = _np.concatenate(([0], 10000.0 / _np.array([26500.0, 12200.0,
        6000.0, 5470.0, 4670.0, 4110.0])))
    ysplir = _np.array([0.0, 0.26469, 0.82925]) * Rv / 3.1 
    ysplop = _np.array((_np.polyval([-4.22809e-01, 1.00270, 2.13572e-04][::-1],
        Rv ), _np.polyval([-5.13540e-02, 1.00216, -7.35778e-05][::-1], Rv ), 
        _np.polyval([ 7.00127e-01, 1.00184, -3.32598e-05][::-1], Rv ), 
        _np.polyval([ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, 
            -4.45636e-05][::-1], Rv ) ))
    ysplopir = _np.concatenate((ysplir, ysplop))

    if (Nopir > 0): 
        tck = _interpolate.splrep(_np.concatenate((xsplopir, xspluv)),
            _np.concatenate((ysplopir, yspluv)), s=0)
        curve[iopir] = _interpolate.splev(x[iopir], tck)

    # Now apply extinction correction to input flux vector
    curve *= -ebv

    return flux * 10.**(0.4 * curve)


def sort_specs(specs, path=None):
    """ Specs in an (N,2) array, where specs[:,0] are the files paths and 
    specs[:,1] the instrument name. 

    Return ordered_specs"""
    if path is not None:
        if path[-1] != '/': 
            path += '/'
    else:
        path = ''
    nsp = _np.shape(specs)[0]
    MJDs = _np.zeros(nsp)
    specs = _np.array(specs)
    lims = [_np.inf, -_np.inf]
    for i in range(nsp):
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(path + 
            specs[i][0])
        MJDs[i] = MJD
        if MJDs[i] < lims[0]:
            lims[0] = MJDs[i]
        if MJDs[i] > lims[1]:
            lims[1] = MJDs[i]
    return specs[MJDs.argsort()], lims


def convgaussFunc(wl, flx, lbc, hwidth=1000., convgauss=0., frac=0., ssize=.05, 
    wlout=False):
    """ Do a Gaussian convolution of a given Line Profile with a Gaussian. 

    `wl`, `flx`, `lbc`: wavelenght and flux of the spectrum containing the 
    line, and its value.

    `hwidth`, `ssize`: width to be keeped around the line (km/s), and the 
    region (in percentage) where the continuum level will be evaluted around 
    the selected region.

    `convgauss`: if bigger then 0., do the convolution. Its values is the sigma 
    of the gaussian conv. profile (in km/s).

    `frac`: controls the intensity of the convolution. `frac`=0 means pure 
    profile output and `frac`=1 a pure gaussian output with the same EW value.

    `wlout`: returns a wavelength array instead of a velocity array (standar)

    OUTPUT: vel/wl, flux (arrays)
    """
    (x, yo) = lineProf(wl, flx, lbc=lbc, hwidth=hwidth + 3 * convgauss, 
        ssize=ssize)
    y1 = yo
    y2 = 0.
    if convgauss > 0 and frac > 0:
        step = _np.min([x[j + 1] - x[j] for j in range(len(x) - 1)])
        xn = _np.arange(-hwidth - 3 * convgauss,
                        hwidth + 3 * convgauss + step, step)
        cf = _phc.normgauss(convgauss, x=xn)
        yo = _np.interp(xn, x, yo)
        x = xn
        y1 = yo * (1 - frac)
        y2 = _np.convolve(yo * frac, cf / _np.trapz(cf), 'same')
    if wlout:
        x = (x / _phc.c.cgs * 1e5 + 1) * lbc
    return x, y1 + y2


def cutpastrefspec(ivl, iflx, irefvl, ireflx, hwidth, ssize=.05):
    """ Cut and paste a given line profile into a reference line profile. 

    Both profiles (with any resolution) must be normalized and given in vel. 

    It was designed to solve the problem of Achernar's Halpha line wings 
    problem and it works like this: given a reference profile (`refvl`, 
    `reflx`), the selected profile will be cutted at the `hwidth` position 
    and them pasted in the corresponding position (and intensity level) of 
    the reference spectrum.

    OUTPUT: refvl, reflx
    """
    flx = _np.interp(irefvl, ivl, iflx)
    i0 = _np.abs(irefvl + hwidth).argmin()
    i1 = _np.abs(irefvl - hwidth).argmin()
    ssize = int(ssize * len(flx))
    if ssize == 0:
        ssize = 1
    refav = _np.average( ireflx[i0 - ssize / 2:i0 + ssize / 2 + 1] ) / 2. + \
        _np.average( ireflx[i1 - ssize / 2:i1 + ssize / 2 + 1] ) / 2.
    av = _np.average( flx[i0 - ssize / 2:i0 + ssize / 2 + 1] ) / 2. + \
        _np.average( flx[i1 - ssize / 2:i1 + ssize / 2 + 1] ) / 2.
    flx += refav - av
    reflx = _np.array(ireflx).copy()
    reflx[i0:i1 + 1] = flx[i0:i1 + 1]
    return irefvl, reflx


def load_specs_fits(speclist, ref, lbc, lncore=None, hwidth=None, 
    gaussfit=False, plotcut=0):
    """ Load a list of specs and do the *line core cut & paste*

    `lncore`: cut and paste hwidth of the line center. It can be None, and 
    must be < hwidth. If hwidth is None, it is assumed to be 1000 km/s. 

    `speclist` : [ ['path+file.fits', 'INSTRUMENT'], ... ]

    `ref`: reference spectra to do the cut & paste 

    `plotcut`: if plotcut > 0, save the cutted spectra in steps of this 
    variable.

    OUTPUT: dtb_obs
    """
    if hwidth is None:
        hwidth = 1000.
    # do core cut?
    docore = lncore < hwidth
    if lncore is None:
        docore = False
    # load ref
    refwl, reflx = loadfits(ref[0])[0:2]
    refvl, reflx = lineProf(refwl, reflx, lbc=lbc)
    # load specs
    dtb_obs = Spec(lbc=lbc, hwidth=hwidth, gaussfit=gaussfit)
    for i in range(_np.shape(speclist)[0]):
        dtb_obs.loadspec(speclist[i][0])
        vl, flx = lineProf(dtb_obs.wl, dtb_obs.flux, lbc=lbc)
        if docore:
            cuted = cutpastrefspec(vl, flx, refvl, reflx, lncore)
            dtb_obs.flux = cuted[1]
            dtb_obs.wl = (cuted[0]/_phc.c.cgs*1e5+1)*lbc
            (dtb_obs.EW, dtb_obs.EC, dtb_obs.VR, dtb_obs.peaksep, 
            dtb_obs.depthcent, dtb_obs.F0) = analline(dtb_obs.wl, 
            dtb_obs.flux, dtb_obs.lbc, hwidth=lncore, verb=False, 
            gaussfit=dtb_obs.gaussfit)
        else:
            (dtb_obs.EW, dtb_obs.EC, dtb_obs.VR, dtb_obs.peaksep, 
            dtb_obs.depthcent, dtb_obs.F0) = analline(dtb_obs.wl, 
            dtb_obs.flux, dtb_obs.lbc, hwidth=hwidth, verb=False, 
            gaussfit=dtb_obs.gaussfit)
        dtb_obs.addspec()
    # complementary plot
    if plotcut > 0 and docore:
        fig0, ax = _plt.subplots()
        for i in range(_np.shape(speclist)[0]):
            dtb_obs.loadspec(speclist[i][0])
            vl, flx = lineProf(dtb_obs.wl, dtb_obs.flux, lbc=lbc)
            cuted = cutpastrefspec(vl, flx, refvl, reflx, lncore)
            if i % plotcut == 0:
                ax.plot(cuted[0], cuted[1])
        _phc.savefig(fig0)
    return dtb_obs


def plot_spec_info(speclist, dtb_obs, mAEW=False, mgray=None):
    """ Standard plot of the Spec class (EW, E/C, V/R, peak-sep., FWHM, F0) 

    OUTPUT: figure (fig pyplot)
    """
    if mAEW:
        dtb_obs.data[:, 1] *= 1000*dtb_obs.lbc/_phc.c.cgs*1e5
    # Legend, Markers and colors idx...
    instm = list(_np.unique(speclist[:, 1]))
    # coridx = [ phc.cycles(instm.index(i)) for i in speclist[:, 1]]
    cores = _phc.gradColor(range(len(instm)), cmapn='inferno')
    coridx = [ cores[instm.index(i)] for i in speclist[:, 1] ]
    coridx = _np.array(coridx)
    mkidx = [ _phc.cycles(instm.index(i), 'mk') for i in speclist[:, 1]]
    mkidx = _np.array(mkidx)
    # Plots
    fig = _plt.figure()
    lins, cols = (7, 1)
    gssteps = [slice(0, 2), 2, 3, 4, 5, 6]
    gs = _gridspec.GridSpec(lins, cols)
    axs = [_plt.subplot(gs[g, :]) for g in gssteps]
    # EW
    axs[0].invert_yaxis()
    axs[-1].set_xlabel('Julian date - 2400000.5')

    ylabels = [u'EW (m\u00c5)', 'E/C', 'V/R', ('pk. sep.'+'\n'+'(km/s)'), 
        'FWHM'+'\n'+'(km/s)', r'F${\lambda 0}$']
    for i, ax in enumerate(axs):
        # binned
        x, y = _phc.bindata(dtb_obs.data[:, 0], dtb_obs.data[:, i+1])
        # yi = _savgol(y, 3, 1)
        ax.plot(x, y, color='gray', zorder=0)
        # points
        for uniquem in set(mkidx):
            idx = _np.where(mkidx == uniquem)
            ax.plot(dtb_obs.data[:, 0][idx], dtb_obs.data[:, i+1][idx], 
                color=coridx[idx][0], marker=uniquem, ls='')
        ax.set_ylabel(ylabels[i])
    #
    xlim = axs[0].get_xlim()
    axs[2].plot(xlim, [1, 1], ls=":", color='k', zorder=1)
    for i in range(1, len(axs)):
        # ax.locator_params(axis='y', nbins=4)
        axs[i].yaxis.set_major_locator(_MaxNLocator(nbins=4, prune='upper')) 
        if i in [1, 2, 3]:
            axs[i].get_yticklabels()[-1].set_visible(False)
    for ax in axs[:-1]:
        ax.set_xticklabels([])
    # Legend
    for i in range(len(instm)):
        # axs[0].plot([np.NaN], [np.NaN], label=instm[i], color=phc.cycles(i), 
            # marker=phc.cycles(i, 'mk'), ls='')
        axs[0].plot([_np.NaN], [_np.NaN], label=instm[i], color=cores[i], 
            marker=_phc.cycles(i, 'mk'), ls='')
    axs[0].legend(loc='best', fancybox=True, framealpha=0.5, fontsize=8, 
        labelspacing=0.05, ncol=2)
        # bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.
    fig.subplots_adjust(hspace=0.01)
    # Gray
    for ax in axs:
        ax.set_xlim(xlim)
        if mgray is not None:
            ylim = ax.get_ylim()
            rect = _mpatches.Rectangle([mgray[0], ylim[0]], 
                mgray[1]-mgray[0], ylim[1]-ylim[0], ec="gray", fc='gray', 
                alpha=0.5, zorder=1)
            ax.add_patch(rect)
        if len(mgray) == 4:
            if mgray is not None:
                ylim = ax.get_ylim()
                rect = _mpatches.Rectangle([mgray[2], ylim[0]], 
                    mgray[3]-mgray[2], ylim[1]-ylim[0], ec="gray", fc='gray', 
                    alpha=0.5, zorder=1, hatch='//')
                ax.add_patch(rect)
    return fig


# TODO: Check if obsolete
def normalize_range(lb, spec, a, b):
    """This function is obsolete and must be removed.

    Still here for compatibility issues.
    """
    a2 = (spec[b] - spec[a]) / (lb[b] - lb[a])
    a1 = spec[a] - a2 * lb[a]
    return spec / (a1 + a2 * lb)


def normalize_spec(lb, flx, q=2, diff=0.03, perc=0, nlbp=50):
    """ Normalize a spectrum using the non-parametric regression algorithm of
    Local Polynomial Kernel (order=``q``). 

    If perc > 0, a "percentile filter" is applyed to the spectrum (divided in
    nlbp bins).

    For details, see http://pythonhosted.org/PyQt-Fit/NonParam_tut.html .

    INPUT: lb, flx

    OUTPUT: norm_flx
    """
    if perc <= 0:
        k1 = _smooth.NonParamRegression(lb, flx, 
            method=_npr_methods.LocalPolynomialKernel(q=1))
        k1.fit()

        idx0 = _np.where(flx != 0)
        ilb = lb[idx0]
        iflx = flx[idx0]
        idxi = _np.where(_np.abs(k1(ilb)/iflx-1) < diff)
        xsi = ilb[idxi]
        ysi = iflx[idxi]
    else:
        xsi, ysi = _phc.bindata(lb, flx, nbins=nlbp, perc=perc)

    k2 = _smooth.NonParamRegression(xsi, ysi, 
        method=_npr_methods.LocalPolynomialKernel(q=q))
    k2.fit()
    return flx/k2(lb)


def renorm(vl, y):
    """ Renormalize ``y`` so that the equivalent width is preserved when the 
    continuum is shifted to 1. 
    """
    ext = _np.mean([y[0], y[-1]])
    a0 = _np.trapz(y, vl)
    A = ((a0-_np.trapz(_np.tile(1, len(vl)), vl))/
        (a0-_np.trapz(_np.tile(ext, len(vl)), vl)))
    B = 1-A*ext
    return A*y+B


def normEW(vl, y, area=None):
    """ Normalize ``y`` curve to have a specific area. If ``area is None``, 
    then the normalized equivalent width is preserved.
    """
    if area is None:
        area = _np.trapz(linfit(vl, y), vl)
    y0 = linfit(vl, y)-1
    a1 = _np.trapz(y0, vl)
    a0 = _np.trapz(_np.tile([1], len(vl)), vl)
    f = (area-a0)/a1
    return f*y0+1


def checksubdirs(path, star, lbc, hwidth=1000, showleg=True, plots=False):
    """
    Faz o que tem que fazer.
    """
    if not _os.path.exists('{0}/{1}'.format(path, star)):
        _os.system('mkdir {0}/{1}'.format(path, star))

    nights = [o for o in _os.listdir(path) if _os.path.isdir('{0}/{1}'.
        format(path, o))]

    fig = _plt.figure()
    ax = fig.add_subplot(111)
    spdtb = Spec()
    spdtb.lbc = lbc
    spdtb.hwidth = 1000.
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path, night)) if
            _os.path.isdir('%s/%s/%s' % (path, night, o))]
        for target in targets:
            if target.find(star) > -1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path, night, target))
                if len(scal) > 0:
                    for cal in scal:
                        spdtb.loadspec(cal)
                        spdtb.addspec()
                        if not _np.isnan(spdtb.EW):
                            if plots:
                                spdtb.plotspec()
                            vels = (spdtb.wl - lbc) / lbc * _phc.c.cgs * 1e-5
                            idx = _np.where(_np.abs(vels) <= hwidth)
                            flux = linfit(vels[idx], spdtb.flux[idx])
                            vels = vels[idx]
                            leg = spdtb.MJD
                            ax.plot(vels, flux, label=leg, alpha=0.7, 
                                color=_phc.colors[_np.mod(spdtb.count, 
                                len(_phc.colors))])
                else:   
                    print('# Data not reduced for %s at %s!' % (star, night))
    ax.set_xlim([-hwidth, hwidth])
    ax.set_ylim([-1, 5])
    if showleg:
        legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
        _plt.setp(legend.get_texts(), fontsize='small')    
    _plt.savefig('{0}/{1}_at_{2}.png'.format(_outfold, star, lbc))
    _plt.close()
    spdtb.savedata(datafile='{0}/{1}.txt'.format(_outfold, star),
        metafile='{0}/meta_{1}.txt'.format(_outfold, star))
    return


def VREWcalc(vels, flux, vw=1000):
    """
    Supoe que o fluxo jah estah normalizado, e vetores ordenad_os.

    Calcula o ew para os dois lados (azul/vermelho) da linha, ajustando
    a velocidade de repouso (TBD).
    """
    # calcula e aplica correcao de vel. repousp
    vc = 0.
    vels += vc
    # corta em vw, e faz o teste de tamanho
    if len(vels) < 5:
        vw = 0
    if vw > 0:
        idx = _np.where(_np.abs(vels) <= vw)
        outvels = vels[idx]
        normflux = flux[idx]
    else:
        ew0 = 0.
        ew1 = 0.
        return ew0, ew1, vc
    #
    ivc = _np.abs(outvels - 0).argmin()
    ew0 = 0.
    for i in range(0, ivc):
        dl = outvels[i + 1] - outvels[i]
        ew0 += (1. - (normflux[i + 1] + normflux[i]) / 2.) * dl
    ew1 = 0.
    for i in range(ivc, len(outvels) - 1):
        dl = outvels[i + 1] - outvels[i]
        ew1 += (1. - (normflux[i + 1] + normflux[i]) / 2.) * dl
    return ew0, ew1, vc


def plotSpecData(dtb, limits=None, civcfg=[1, 'm', 2013, 1, 1],
    fmt=['png'], ident=None, lims=None, setylim=False, addsuf=''):
    """ Plot spec class database `vs` MJD e civil date

    Plot originally done to London, Canada, 2014.

    INPUT: civcfg = [step, 'd'/'m'/'y', starting year, month, day]

    `lims` sequence: 'EW', 'E/C', 'V/R', 'Pk. sep. (km/s)', 'E-F0', 'F0'

    `lims` = [[-2,4+2,2],[1.,1.4+.1,0.1],[.6,1.4+.2,.2],[0,400+100,100],
    [.30,.45+.05,.05],[0.6,1.20+.2,.2]]

    If `lims` is defined, `setylim` can be set to True.

    OUTPUT: Written image."""
    if isinstance(dtb, _strtypes):
        print('# Loading dtb {0}'.format(dtb))
        dtb = _np.loadtxt(dtb)
    if ident is not None:
        idref = _np.unique(ident)

    ylabels = ['EW', 'E/C', 'V/R', 'Pk. sep. (km/s)', 'E-F0', 'F0']
    fig, ax = _plt.subplots(6, 1, sharex=True, figsize=(9.6, 8))

    icolor = 'blue'
    for i in range(1, len(ylabels) + 1):
        ax[i - 1].plot(*_phc.bindata(dtb[:, 0], dtb[:, i], 20))
        for j in range(len(dtb[:, 0])):
            if ident is not None:
                idx = _np.where(ident[j] == idref)[0]
                icolor = _phc.colors[idx]
            ax[i - 1].plot(dtb[j, 0], dtb[j, i], 'o', color=icolor)
        ax[i - 1].set_ylabel(ylabels[i - 1])
        if lims is not None:
            if lims[i - 1][-1] != 0:
                ax[i - 1].set_yticks(_np.arange(*lims[i - 1]))
            if setylim:
                ax[i - 1].set_ylim([ lims[i - 1][0], lims[i - 1][1] ])

    if ident is not None:
        for id in idref:
            idx = _np.where(id == idref)[0]
            icolor = _phc.colors[idx]
            ax[0].plot([], [], 'o', color=icolor, label=id)
        ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., 
            prop={'size': 6})

    if limits is None:
        # limits = ax[0].get_xlim()
        limits = [dtb[0, 0], dtb[-1, 0]]
    else:
        ax[0].set_xlim(limits)
    mjd0, mjd1 = limits
    ax[5].set_xlabel('MJD')
    ticks = _phc.gentkdates(mjd0, mjd1, civcfg[0], civcfg[1], 
        dtstart=_dt.datetime(civcfg[2], civcfg[3], civcfg[4]).date())
    mjdticks = [_jdcal.gcal2jd(date.year, date.month, date.day)[1] for date in 
        ticks]
    # ticks = [dt.datetime(*jdcal.jd2gcal(jdcal.MJD_0, date)[:3]).date() for \
    # date in ax[0].get_xticks()]
    # mjdticks = ax[0].get_xticks()
    for i in range(1, 6 + 1):
        ax2 = ax[i - 1].twiny()
        ax2.set_xlim(limits)
        ax2.set_xticks(mjdticks)
        ax2.set_xticklabels(['' for date in ticks])
        if i == 1:
            ax2.set_xlabel('Civil date')
            ax2.set_xticklabels([date.strftime("%d %b %y") for date in ticks])
            _plt.setp( ax2.xaxis.get_majorticklabels(), rotation=45 )
    _plt.subplots_adjust(left=0.13, right=0.8, top=0.88, bottom=0.06, 
        hspace=.15)
    for f in fmt:
        print ('SpecQ{1}.{0}'.format(f, addsuf))
        _plt.savefig('SpecQ{1}.{0}'.format(f, addsuf), transparent=True)
    _plt.close()
    return


def din_spec(metadata, lbc=6562.86, hwidth=1500., res=50, interv=None,
    fmt=['png'], outname='din_spec', pxsize=8, vmin=None, vmax=None, avg=True,
    cmapn='inferno', refspec=None, figsize=None):
    """ Plot dynamical specs. from metadata table of the Spec class.

    `interv` controls the interval between specs (in days).

    `res` is the resolution in km/s.

    By default (`avg`=True), the average of spectra in that bin is show. If 
    `avg`=False, the nearest bin-centered (in time) spectra will be shown.

    if `refspec` is not None, them it will be a difference spectra.
    """
    # Define MJD and bins
    dates = _np.array(metadata[:, 0], dtype=float)
    t0 = _np.min(dates)
    tf = _np.max(dates)
    if interv is None:
        interv = _np.linspace(t0, tf, 21)
    else:
        interv = _np.arange(t0, tf + interv, interv)
    dt = interv[1] - interv[0]
    # Select specs 
    wl0 = _np.arange(-hwidth, hwidth + res, res)
    # Load refspec, if required
    baselevel = 1.
    if refspec is not None:
        wl, flux, tmp, tmp, tmp, tmp = loadfits(refspec)
        wl, flux = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        refflx = _np.interp(wl0, wl, flux)
        baselevel = 0
    fluxes = _np.zeros(( len(wl0), len(interv) )) + baselevel
    for i in range(len(interv)):
        # method 1
        if not avg:
            date = _phc.find_nearest(dates, interv[i])
            if date < interv[i] + dt / 2 and date > interv[i] - dt / 2:
                j = list(dates).index(date)
                wl, flux, tmp, tmp, tmp, tmp = loadfits(metadata[j, 3])
                wl, flux = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
                if refspec is None:
                    fluxes[:, i] = _np.interp(wl0, wl, flux)
                else:
                    flux = _np.interp(wl0, wl, flux)
                    fluxes[:, i] = flux - refflx
        # method 2
        else:
            k = 0
            for j in range(len(dates)):
                if dates[j] < interv[i] + dt / 2 and dates[j] > interv[i] - \
                    dt / 2:
                    wl, flux, tmp, tmp, tmp, tmp = loadfits(metadata[j, 3])
                    wl, flux = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
                    fluxes[:, i] += _np.interp(wl0, wl, flux)
                    k += 1
            if k > 0:
                # fluxes[:,i]/= k
                wl = vel2wl(wl0, lbc)
                tmp, fluxes[:, i] = lineProf(wl, fluxes[:, i], lbc=lbc, 
                    hwidth=hwidth)
                if refspec is not None:
                    fluxes[:, i] = fluxes[:, i] - refflx
        if all(fluxes[:, i] == baselevel):
            fluxes[:, i] = _np.NaN
    # Create image
    img = _np.empty((pxsize * len(interv), len(wl0)))
    for i in range(len(interv)):
        img[i * pxsize:(i + 1) * pxsize] = _np.tile(fluxes[:, i], pxsize).\
            reshape(pxsize, len(wl0))
    # Save image
    if figsize is None:
        fig, ax = _plt.subplots(figsize=(len(wl0) / 16, pxsize * 
            len(interv) / 16), dpi=80)
    else:
        fig, ax = _plt.subplots(figsize=figsize)
    # _plt.figure(figsize=(len(wl0) / 16, pxsize * len(interv) / 16), dpi=80)
    # print _np.min(img), _np.max(img)
    cmapn = _plt.get_cmap(cmapn)
    cmapn.set_bad('k', 1.)
    ax.imshow(img, vmin=vmin, vmax=vmax, cmap=cmapn, origin='lower')
    ax.set_xlabel(r'Velocity (km s$^{-1}$)')
    ax.set_ylabel(r'Julian Day - 2400000.5')
    # ax.set_xlim([-hwidth, hwidth])
    ax.set_yticks(_np.linspace(pxsize*len(interv)*.1, pxsize*len(interv)*.9, 
        8))
    ax.set_yticklabels([int(round((tf-t0)*t/(pxsize*len(interv))+t0)) 
        for t in ax.get_yticks()], rotation='vertical')
    ax.set_xticklabels([int(round(t*2.*hwidth/(len(wl0)-1)-hwidth)) for 
        t in ax.get_xticks()])  # , rotation='vertical')
    # fig.tight_layout()
    ax.xaxis.set_tick_params(color='gray', width=1.1)
    ax.yaxis.set_tick_params(color='gray', width=1.1)
    _phc.savefig(fig, fmt=fmt, figname=outname)
    return


def plot_line_str(fig, ax, lbc='', ylabel='', fs=14, xlim=None, dlim=None, 
    cmapn='gnuplot', lfs=10, ylim=None):
    """ Line plotting structure """
    if lbc is not '':
        ax.set_title(r'$\lambda_c$ = {0:.1f} $\AA$'.format(lbc), size=fs)
    if ylabel is not '':
        ax.set_ylabel(ylabel, size=fs)

    if xlim is not None:
        ax.xlims = ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel(r'Velocity (km s$^{-1}$)', size=fs)
    # reverse to keep order consistent
    ax.legend()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper right', labelspacing=0.1, 
        fancybox=True, framealpha=0.5, fontsize=lfs)  # loc=(1.05, .01)

    rect = _mpatches.Rectangle([0.835, 0.01], 0.15, 0.44, ec="black", 
        fc='white', transform=ax.transAxes, zorder=10, alpha=0.5)
    ax.add_patch(rect)

    ax3 = fig.add_axes([0.82, 0.12, 0.025, 0.35])
    # ax3.set_axis_bgcolor('white')
    cmap = _plt.get_cmap(cmapn)
    norm = _mpl.colors.Normalize(vmin=dlim[0], vmax=dlim[1])
    cb = _mpl.colorbar.ColorbarBase(ax3, cmap=cmap, norm=norm, 
        orientation='vertical')
    cb.set_label('MJD', size=fs) 
    fig.subplots_adjust(left=0.1, right=0.95, top=0.94, bottom=0.1)  
    # , hspace=0.3, wspace=.3)  
    return fig, ax


def spec_time(speclist, lbc=6562.8, fmt=['png', 'pdf'], outname=None, 
    cmapn='inferno', hwidth=1000., outpath='', figsize=(5, 15), ysh=0.01):
    """ Plot specs over time as suggested by Rivi """
    if outname is None or outname is "":
        outname = _phc.dtflag()
    MJDs = [_np.inf, 0]
    for sp in speclist:
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(sp)
        if MJD < MJDs[0]:
            MJDs[0] = MJD
        if MJD > MJDs[1]:
            MJDs[1] = MJD
    MJDref = 56245
    if MJDs[0] > MJDref:
        MJDs[0] = MJDref
    # Plot
    extrem = [_np.inf, 0]
    fig, ax = _plt.subplots()
    for sp in speclist:
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(sp)
        vel, flux = lineProf(wl, flux, lbc, hwidth=hwidth)
        if len(flux) == 0:
            raise NameError('Wrong lbc in spt.spe')
        if cmapn is not None:
            cor = _phc.gradColor([MJD], min=MJDs[0], max=(MJDs[1]+
                0.1*(MJDs[1]-MJDs[0])), cmapn=cmapn)[0]
        else:
            cor = 'k'
        print(MJD, MJDs, extrem, ysh, (MJD-MJDs[0])*ysh, flux, sp)
        ax.plot(vel, flux+(MJD-MJDs[0])*ysh, color=cor)
        if _np.max(flux+(MJD-MJDs[0])*ysh) > extrem[1]:
            extrem[1] = _np.max(flux+(MJD-MJDs[0])*ysh)
        if _np.min(flux+(MJD-MJDs[0])*ysh) < extrem[0]:
            extrem[0] = _np.min(flux+(MJD-MJDs[0])*ysh)
    wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits('/data/Dropbox/work'
        '/sci_16-15aeri/alpEri_FEROS_2000AVE.mt')
    vel, flux = lineProf(wl, flux, 6561.8, hwidth=hwidth)
    ax.text(650., 0.8, 'photospheric ref.', horizontalalignment='center', 
        verticalalignment='center')  # , transform=ax.transAxes)
    ax.plot(vel, flux+(MJDref-MJDs[0])*ysh, color='k', ls=':')
    if _np.min(flux+(MJDref-MJDs[0])*ysh) < extrem[0]:
        extrem[0] = _np.min(flux+(MJDref-MJDs[0])*ysh)
    s2d = _hdt.readfullsed2('/data/Dropbox/work/sci_16-15aeri/'
        'fullsed_mod03_VDDn0_1p4e12_Be_aeri2014.sed2')
    vel, flux = lineProf(s2d[4, :, 2], s2d[4, :, 3], .656461, hwidth=hwidth)
    ax.plot(vel, flux+(56910-MJDs[0])*ysh, color='k', ls='--')
    ax.text(800, 1.06+(56910-MJDs[0])*ysh, 'model', 
        horizontalalignment='center', verticalalignment='center') 
    ax.set_xlabel(r'Velocity (km s$^{-1}$)')
    ax.set_ylabel(r'Julian Day - 2400000.5')
    ax.set_ylim(extrem)
    ax.set_xlim([-hwidth, hwidth])
    # ax.set_yticks(_np.arange(56300, 57000+100, 100))
    yref = [1., 1+_np.diff(MJDs)*ysh]
    yMJDs = _np.arange(56300, 57100, 100)
    ax.set_yticks(list(_phc.renormvals(yMJDs, MJDs, yref)))
    ax.set_yticklabels(yMJDs, rotation='vertical')
    fig.set_size_inches(figsize)
    fig.subplots_adjust(left=0.1, right=0.94, top=0.99, bottom=0.04)
    ax.minorticks_on()
    ax3 = ax.twinx()
    ax3.set_yticks(list(_phc.renormvals(yMJDs, MJDs, yref)))
    ax3.set_yticklabels([])
    ax3.minorticks_on()
    ax2 = ax.twinx()
    ax2.spines['right'].set_position(('axes', 1.05))
    ax2.set_ylabel('Civil date')
    dtminticks = _phc.gentkdates(56201., 57023., 1, 'm')
    i = 1
    dtticks = _phc.gentkdates(56201., 57023., 3, 'm')
    mjdticks = [_jdcal.gcal2jd(date.year, date.month, date.day)[1] for date in 
        dtticks]
    while dtticks[0] not in dtminticks:
        dtminticks = _phc.gentkdates(yMJDs[0]+i, yMJDs[-1], 1, 'm')
        i += 1
    minjdticks = [_jdcal.gcal2jd(date.year, date.month, date.day)[1] for date 
        in dtminticks]
    ax2.set_yticks(list(_phc.renormvals(mjdticks, MJDs, yref)))
    ax2.set_yticks(list(_phc.renormvals(minjdticks, MJDs, yref)), minor=True)
    xlabs = [date.strftime('%Y-%m-%d') for date in dtticks]
    # xlabs[1::2] = ['']*len(xlabs[1::2])
    ax2.set_yticklabels(xlabs, rotation='vertical')
    ax2.set_ylim(extrem)
    ax3.set_ylim(extrem)
    ax.xaxis.set_tick_params(length=8, width=1.5)
    ax.xaxis.set_tick_params(length=6, which='minor')
    ax.yaxis.set_tick_params(length=4, which='minor')
    ax.yaxis.set_tick_params(length=8, width=1.5)
    ax2.yaxis.set_tick_params(length=4, which='minor')
    ax2.yaxis.set_tick_params(length=8, width=1.5)
    ax3.yaxis.set_tick_params(length=4, which='minor')
    ax3.yaxis.set_tick_params(length=8, width=1.5)
        # , fontsize=10)
    _phc.savefig(fig, figname=outpath+outname, fmt=fmt)
    return


def extractfromsplot(file, splot):
    """Ce = center; Co = core
    #LcCe, LcCo, lcGW, lcEW, lvCe, lcCo, lvEW, lrCe, LrCo, lrEW
    """
    out = _np.array(10 * [_np.NaN])
    readflag = False
    for line in splot:
        if line.find(']:') > 0 and readflag:
            readflag = False
        if line.find(file) > 0:
            readflag = True
        if readflag:
            info = line.split()
            # if _re.match("^\d+?\.\d+?$", info[0]) is not None:
            try:
                float(info[0])
                info = _np.array(info, dtype=float)
                if info[0] > 6556 and info[0] < 6556 + 4.33:
                    if len(info) == 4:
                        out[6] = float(info[3])
                    elif len(info) == 7:
                        out[4] = float(info[0])
                        out[5] = float(info[4])
                elif info[0] > 6556 + 4.33 and info[0] < 6556 + 2 * 4.33:
                    if len(info) == 4:
                        out[3] = float(info[3])
                    elif len(info) == 7:
                        out[0] = float(info[0])
                        out[1] = float(info[4])
                        out[2] = float(info[5])
                elif info[0] > 6556 + 2 * 4.33 and info[0] < 6556 + 3 * 4.33:
                    if len(info) == 4:
                        out[9] = float(info[3])
                    elif len(info) == 7:
                        out[7] = float(info[0])
                        out[8] = float(info[4])
            except:
                pass
    return out


def check_dtobs(dtobs):
    """ Check if the dtobs fits the float format. Required for MJD calc. """
    if 'T' in dtobs:
        dtobs = dtobs.replace('.', '')
        tobs, dtobs = dtobs.split('T')
        if len(tobs) == 10:
            dtobs, tobs = tobs, dtobs
        tobs = tobs.split(':')
        tobs = float(tobs[0]) * 3600 + float(tobs[1]) * 60 + float(tobs[2])
        tobs /= (24 * 3600)
    else:
        tobs = 0.
    if dtobs[4] == '-':
        dtobs = dtobs.split('-')
    elif dtobs[2] == '/':
        dtobs = dtobs.split('/')[::-1]
    else:
        _warn.warn('Wrong "DATE-OBS" in header! {0}'.format(dtobs))
        raise SystemExit(1)
    dtobs = _np.array(dtobs, dtype='int32')
    return dtobs, tobs


# TODO: Check if obsolete
def overplotsubdirs(path, star, limits=(6540, 6600), showleg=True):
    """
    Realiza o plot de espectros da estrela `star` dentre do diretorio `path`.
    Atualmente, faz o plot entre os valores `limits` (Angstroms).

    Gera os arquivos `path/star/star.log` e `path/star/star_specs.png`.
    """
    # path = _os.getcwd()
    # star = _phc.user_input('Type the star name: ')
    # ref0 = 6540
    # ref1 = 6600

    ref0, ref1 = limits

    if not _os.path.exists('{0}/{1}'.format(path, star)):
        _os.system('mkdir {0}/{1}'.format(path, star))
    f0 = open('{0}/{1}/{1}.log'.format(path, star), 'w')

    nights = [o for o in _os.listdir(path) if _os.path.isdir('{0}/{1}'.
        format(path, o))]

    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path, night)) if
            _os.path.isdir('%s/%s/%s' % (path, night, o))]
        for target in targets:
            if target.find(star) > -1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path, night, target))
                if len(scal) > 0:            
                    srv = _glob('%s/%s/%s/*.rv.fits' % (path, night, target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    # legendl += (night,)
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec)) * \
                            imfits[0].header['CDELT1'] + \
                            imfits[0].header['CRVAL1']
                        # a = _phc.user_input('type to continue: ')
                        if lbda[-1] > 6560:  # and flag == '1':
                            min_dif = min(abs(lbda - ref0))
                            a0 = _np.where(abs(lbda - ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda - ref1))
                            a1 = _np.where(abs(lbda - ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda, spec, a0, a1)
                            msg = '{0}, {1}, {2}'.format((0.1 * i), night, cal)
                            print(msg)
                            f0.writelines(msg + '\n')
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            _plt.plot(lbda, spec, label=leg, alpha=0.7,
                                color=_phc.colors[_np.mod(i, 
                                    len(_phc.colors))])
                            i += 1
                else:
                    print('# Data not reduced for %s at %s!' % (star, night))
                    msg = '{0}, {1}, {2}'.format('NC', night, 'None')
                    f0.writelines(msg + '\n')

    if showleg:
        legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
        _plt.setp(legend.get_texts(), fontsize='small')
    _plt.xlim([ref0, ref1])
    _plt.ylim([-1, 5])
    # _plt.xlabel('vel. (km/s)')
    _plt.savefig('{0}/{1}/{1}_specs.png'.format(path, star))
    _plt.close()
    f0.close()
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = _pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = _np.arange(len(specs[i]))*imfits[0].header['CDELT1']+
    #       imfits[0].header['CRVAL1']
    #     if Ha:
    #         if i == 0:
    #             lbds[i] = (lbds[i]-6561.5)/6561.5*3e5
    #         else:
    #             lbds[i] = (lbds[i]-6562.8)/6562.8*3e5
    #     else:
    #         if i == 0:
    #             lbds[i] = (lbds[i]-6676.8)/6676.8*3e5
    #         else:        
    #             lbds[i] = (lbds[i]-6678.)/6678.*3e5
    #     
    #     a = _np.where( abs(lbds[i]+1000) ==  min(abs(lbds[i]+1000)) )
    #     b = _np.where( abs(lbds[i]-1000) ==  min(abs(lbds[i]-1000)) )
    #     
    #     specs[i] = normalize_range(lbds[i],specs[i],a,b)
    #     
    #     legendl += [imfits[0].header['DATE-OBS']]
    # 
    # figure(2)
    # for i in range(len(specs)):
    #     plot(lbds[i], specs[i], label=legendl[i])
    # 
    # legend(legendl, 'lower right')
    # legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #'lower right'
    # xlim([-1000,1000])
    # if Ha:
    #     title('Halpha profile from LNA-Janot for Achernar')
    #     ylim([.65,1.1])
    # else:
    #     title('HeI 6678 profile from LNA-Janot for Achernar')
    #     ylim([.9,1.05])
    # 
    # legend = _plt.legend(legendl, loc=(0.75, .05), labelspacing=0.1)
    # _plt.setp(legend.get_texts(),  fontsize='small')
    # 
    # xlabel('vel. (km/s)')

    print('# Plot done!')
    return


def diffplotsubdirs(path, star, limits=(6540, 6600)):
    """
    Realiza o plot de espectros da estrela `star` dentre do diretorio `path`.
    Atualmente, faz o plot entre os valores `limits` (Angstroms).

    Gera os arquivos `path/star/star.log` e `path/star/star_specs_dif.png`.
    """
    ref0, ref1 = limits

    if not _os.path.exists('{0}/{1}'.format(path, star)):
        _os.system('mkdir {0}/{1}'.format(path, star))
    # f0 = open('{0}/{1}/{1}.log'.format(path, star), 'w')

    nights = [o for o in _os.listdir(path) if _os.path.isdir('{0}/{1}'.
        format(path, o))]    

    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path, night)) if 
            _os.path.isdir('%s/%s/%s' % (path, night, o))]
        for target in targets:
            if target.find(star) > -1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path, night, target))
                if len(scal) > 0:            
                    srv = _glob('%s/%s/%s/*.rv.fits' % (path, night, target))
                    if len(srv) != len(scal):
                        print('# Specs with dopcor at %s!' % night)
                        srv = scal
                    # legendl += (night,)
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec)) * imfits[0].\
                            header['CDELT1'] + imfits[0].header['CRVAL1']
                        # a = _phc.user_input('type to continue: ')
                        if lbda[0] > 5500:  # and flag == '1':
                            min_dif = min(abs(lbda - ref0))
                            a0 = _np.where(abs(lbda - ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda - ref1))
                            a1 = _np.where(abs(lbda - ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda, spec, a0, a1) + \
                                (0.1 * i)
                            print (0.1 * i)
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            _plt.plot([ref0, ref1], [1 + 0.1 * i, 1 + 0.1 * i], 
                                'k--', alpha=0.5)
                            _plt.plot(lbda, spec, label=leg, 
                                color=_phc.colors[i])
                            i += 1
                else:
                    print('# Data not reduced for %s at %s!' % (star, night))

    legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
    _plt.setp(legend.get_texts(), fontsize='small')
    _plt.xlim([ref0, ref1])
    # _plt.xlabel('vel. (km/s)')
    _plt.savefig('{0}/{1}/{1}_specs_dif.png'.format(path, star))

    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = _pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = _np.arange(len(specs[i]))*imfits[0].header['CDELT1']+
    #       imfits[0].header['CRVAL1']
    #     if Ha:
    #         if i == 0:
    #             lbds[i] = (lbds[i]-6561.5)/6561.5*3e5
    #         else:
    #             lbds[i] = (lbds[i]-6562.8)/6562.8*3e5
    #     else:
    #         if i == 0:
    #             lbds[i] = (lbds[i]-6676.8)/6676.8*3e5
    #         else:        
    #             lbds[i] = (lbds[i]-6678.)/6678.*3e5
    #     
    #     a = _np.where( abs(lbds[i]+1000) ==  min(abs(lbds[i]+1000)) )
    #     b = _np.where( abs(lbds[i]-1000) ==  min(abs(lbds[i]-1000)) )
    #     
    #     specs[i] = normalize_range(lbds[i],specs[i],a,b)
    #     
    #     legendl += [imfits[0].header['DATE-OBS']]
    # 
    # figure(2)
    # for i in range(len(specs)):
    #     plot(lbds[i], specs[i], label=legendl[i])
    # 
    # legend(legendl, 'lower right')
    # legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #'lower right'
    # xlim([-1000,1000])
    # if Ha:
    #     title('Halpha profile from LNA-Janot for Achernar')
    #     ylim([.65,1.1])
    # else:
    #     title('HeI 6678 profile from LNA-Janot for Achernar')
    #     ylim([.9,1.05])
    # 
    # legend = _plt.legend(legendl, loc=(0.75, .05), labelspacing=0.1)
    # _plt.setp(legend.get_texts(),  fontsize='small')
    # 
    # xlabel('vel. (km/s)')
    print('# Plot done!')
    return


def refplotsubdirs(path, star, limits=(6540, 6600)):
    """
    Realiza o plot de espectros da estrela `star` dentre do diretorio `path`.
    Atualmente, faz o plot entre os valores `limits` (Angstroms).

    Gera os arquivos `path/star/star.log` e 
    `path/star/star_specs_REFERENCIA.png`.
    """
    ref0, ref1 = limits

    if not _os.path.exists('{0}/{1}'.format(path, star)):
        _os.system('mkdir {0}/{1}'.format(path, star))
    f0 = open('{0}/{1}/{1}.log'.format(path, star), 'w')

    nights = [o for o in _os.listdir(path) if 
        _os.path.isdir('{0}/{1}'.format(path, o))]

    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path, night)) if 
            _os.path.isdir('%s/%s/%s' % (path, night, o))]
        for target in targets:
            if target.find(star) > -1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path, night, target))
                if len(scal) > 0:            
                    srv = _glob('%s/%s/%s/*.rv.fits' % (path, night, target))
                    if len(srv) != len(scal):
                        srv = scal
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec)) * imfits[0].\
                            header['CDELT1'] + imfits[0].header['CRVAL1']
                        # a = _phc.user_input('type to continue: ')
                        if lbda[0] > 5500:  # and flag == '1':
                            min_dif = min(abs(lbda - ref0))
                            a0 = _np.where(abs(lbda - ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda - ref1))
                            a1 = _np.where(abs(lbda - ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda, spec, a0, a1)
                            print (0.1 * i)
                            leg = imfits[0].header['DATE-OBS']
                            refleg = '2012-11-20T23:51:37.392'
                            refleg = '2008-06-13'
                            if leg == refleg:
                                f0 = open('{0}/{1}/ref.txt'.format(path, star),
                                    'w')
                                f0.writelines([str(x) + '\t' for x in lbda])
                                f0.writelines('\n')
                                f0.writelines([str(x) + '\t' for x in spec])
                                f0.writelines('\n')
                                f0.close()
                            i += 1
                else:
                    print('# Data not reduced for %s at %s!' % (star, night))

    f0 = open('{0}/{1}/ref.txt'.format(path, star))
    lines = f0.readlines()
    f0.close()
    specref = _np.array(lines[1].split(), dtype=float)
    lbdaref = _np.array(lines[0].split(), dtype=float)
    func = _interpolate.interp1d(lbdaref, specref)  # , kind='cubic')
    lbdaref = _np.linspace(ref0, ref1, 5000)
    specref = func(lbdaref)

    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path, night)) if
            _os.path.isdir('%s/%s/%s' % (path, night, o))]
        for target in targets:
            if target.find(star) > -1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path, night, target))
                if len(scal) > 0:            
                    srv = _glob('%s/%s/%s/*.rv.fits' % (path, night, target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    # legendl += (night,)
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec)) * imfits[0].\
                            header['CDELT1'] + imfits[0].header['CRVAL1']
                        # a = _phc.user_input('type to continue: ')
                        if lbda[0] > 5500:  # and flag == '1':
                            min_dif = min(abs(lbda - ref0))
                            a0 = _np.where(abs(lbda - ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda - ref1))
                            a1 = _np.where(abs(lbda - ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda, spec, a0, a1)
                            func = _interpolate.interp1d(lbda, spec)  
                            # , kind='cubic')
                            # Tive problemas de 'out-of-bounds'... um espectro 
                            # estava desordenado:
                            # print imfits[0].header['CDELT1'], 
                            # imfits[0].header['CRVAL1'], cal
                            spec = func(lbdaref)
                            print (0.1 * i)
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            if i < 130:
                                _plt.plot(lbdaref, spec - specref, label=leg, 
                                    alpha=0.8, color=_phc.colors[i])
                            i += 1
                else:
                    print('# Data not reduced for %s at %s!' % (star, night))

    legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
    _plt.setp(legend.get_texts(), fontsize='small')
    _plt.xlim([ref0, ref1])
    _plt.title('Ref.= %s' % refleg)
    # _plt.xlabel('vel. (km/s)')
    _plt.savefig('{0}/{1}/{1}_specs_{2}.png'.format(path, star, refleg[:10]))

    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = _pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = _np.arange(len(specs[i]))*imfits[0].header['CDELT1']+\
    #          imfits[0].header['CRVAL1']
    #     if Ha:
    #         if i == 0:
    #             lbds[i] = (lbds[i]-6561.5)/6561.5*3e5
    #         else:
    #             lbds[i] = (lbds[i]-6562.8)/6562.8*3e5
    #     else:
    #         if i == 0:
    #             lbds[i] = (lbds[i]-6676.8)/6676.8*3e5
    #         else:        
    #             lbds[i] = (lbds[i]-6678.)/6678.*3e5
    #     
    #     a = _np.where( abs(lbds[i]+1000) ==  min(abs(lbds[i]+1000)) )
    #     b = _np.where( abs(lbds[i]-1000) ==  min(abs(lbds[i]-1000)) )
    #     
    #     specs[i] = normalize_range(lbds[i],specs[i],a,b)
    #     
    #     legendl += [imfits[0].header['DATE-OBS']]
    # 
    # figure(2)
    # for i in range(len(specs)):
    #     plot(lbds[i], specs[i], label=legendl[i])
    # 
    # legend(legendl, 'lower right')
    # legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #'lower right'
    # xlim([-1000,1000])
    # if Ha:
    #     title('Halpha profile from LNA-Janot for Achernar')
    #     ylim([.65,1.1])
    # else:
    #     title('HeI 6678 profile from LNA-Janot for Achernar')
    #     ylim([.9,1.05])
    # 
    # legend = _plt.legend(legendl, loc=(0.75, .05), labelspacing=0.1)
    # _plt.setp(legend.get_texts(),  fontsize='small')
    # 
    # xlabel('vel. (km/s)')
    print('# Plot done!')
    return


def overplotsubdirs2(path, star, limits=(6540, 6600)):
    """
    Realiza o plot de espectros da estrela `star` dentre do diretorio `path`.
    Atualmente, faz o plot entre os valores `limits` (Angstroms).

    Ha' um criterio de escolha de espectros aqui (rudimentar).

    Gera os arquivos `path/star/star.log` e `path/star/star_specs2.png`.
    """
    ref0, ref1 = limits

    if not _os.path.exists('{0}/{1}'.format(path, star)):
        _os.system('mkdir {0}/{1}'.format(path, star))
    f0 = open('{0}/{1}/{1}.log'.format(path, star), 'w')

    nights = [o for o in _os.listdir(path) if _os.path.isdir('{0}/{1}'.
        format(path, o))]

    ax = _plt.figure()
    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path, night)) if 
            _os.path.isdir('%s/%s/%s' % (path, night, o))]
        for target in targets:
            if target.find(star) > -1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path, night, target))
                if len(scal) > 0:            
                    srv = _glob('%s/%s/%s/*.rv.fits' % (path, night, target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    # legendl += (night,)
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec)) * imfits[0].\
                            header['CDELT1'] + imfits[0].header['CRVAL1']
                        # a = _phc.user_input('type to continue: ')
                        if lbda[0] > 5500:  # and flag == '1':
                            min_dif = min(abs(lbda - ref0))
                            a0 = _np.where(abs(lbda - ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda - ref1))
                            a1 = _np.where(abs(lbda - ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda, spec, a0, a1)
                            print (0.1 * i), night
                            prtcolor = _phc.colors[i]
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            check = False
                            if leg.find('2012-11-20T23:51:37.392') != -1:
                                leg = '2012-11-20'
                                prtcolor = _phc.colors[0]
                                check = True
                            elif leg.find('22/01/2013') != -1:
                                leg = '2013-01-22'
                                check = True
                            # elif leg.find('03/07/2013') != -1:
                            #     leg = '2013-07-03'
                            #     check = True
                            elif leg.find('28/07/2013') != -1:
                                leg = '2013-07-28'
                                check = True
                            elif leg.find('2013-11-12T01:30:38.938') != -1:
                                leg = '2013-11-12'
                                check = True
                            else:
                                print(leg)
                            if check:
                                print(cal)
                                _plt.plot(lbda, spec, label=leg, alpha=0.7, 
                                    color=prtcolor)
                            i += 1
                else:
                    msg = '# Data not reduced for %s at %s!' % (star, night)
                    print(msg)
                    f0.writelines(msg)

    font = { 'size': 16, }

    legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
    # _plt.setp(legend.get_texts(),  fontsize='small')
    _plt.xlim([ref0, ref1])
    _plt.ylim([.58, 1.2])
    _plt.xlabel(r'wavelength ($\AA$)', fontdict=font)
    _plt.ylabel('Normalized flux', fontdict=font)
    # _plt.xlabel('vel. (km/s)')
    _plt.savefig('{0}/{1}/{1}_specs2.png'.format(path, star))
    _plt.close()
    f0.close()
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = _pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = _np.arange(len(specs[i]))*imfits[0].header['CDELT1']+
    #           imfits[0].header['CRVAL1']
    #     if Ha:
    #         if i == 0:
    #             lbds[i] = (lbds[i]-6561.5)/6561.5*3e5
    #         else:
    #             lbds[i] = (lbds[i]-6562.8)/6562.8*3e5
    #     else:
    #         if i == 0:
    #             lbds[i] = (lbds[i]-6676.8)/6676.8*3e5
    #         else:        
    #             lbds[i] = (lbds[i]-6678.)/6678.*3e5
    #     
    #     a = _np.where( abs(lbds[i]+1000) ==  min(abs(lbds[i]+1000)) )
    #     b = _np.where( abs(lbds[i]-1000) ==  min(abs(lbds[i]-1000)) )
    #     
    #     specs[i] = normalize_range(lbds[i],specs[i],a,b)
    #     
    #     legendl += [imfits[0].header['DATE-OBS']]
    # 
    # figure(2)
    # for i in range(len(specs)):
    #     plot(lbds[i], specs[i], label=legendl[i])
    # 
    # legend(legendl, 'lower right')
    # legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #'lower right'
    # xlim([-1000,1000])
    # if Ha:
    #     title('Halpha profile from LNA-Janot for Achernar')
    #     ylim([.65,1.1])
    # else:
    #     title('HeI 6678 profile from LNA-Janot for Achernar')
    #     ylim([.9,1.05])
    # 
    # legend = _plt.legend(legendl, loc=(0.75, .05), labelspacing=0.1)
    # _plt.setp(legend.get_texts(),  fontsize='small')
    # 
    # xlabel('vel. (km/s)')
    print('# Plot done!')
    return


def overPlotLineSeries(fullseds, obsers=[0], lbc=.6564606, fmt=['png'],
    convgauss=0., frac=0., addsuf='', labels=None, hwidth=1000., ssize=.05,
    outpath='', ylim=[.7, 2.2], cmapn='gnuplot'):
    """Generate overplot spec. line from a HDUST mod list, separated by
    observers.

    Observers config. must be the same between models in `fullseds` list.

    If `convgauss` > 0, do a gaussian convolution.
    """
    if labels is None:
        labels = [''] * len(fullseds)
    for obs in obsers:
        fig, ax = _plt.subplots()
        fig2, ax2 = _plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = _hdt.readfullsed2(file)
            obsdegs = (_np.arccos(sed2data[:, 0, 0]) * 180 / _np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x, yo) = lineProf(sed2data[obs, :, 2], sed2data[obs, :, 3], 
                lbc=lbc, hwidth=hwidth + 3 * convgauss, ssize=ssize)
            y1 = yo
            y2 = 0.
            if convgauss > 0:
                step = _np.min([x[j + 1] - x[j] for j in range(len(x) - 1)])
                xn = _np.arange(-hwidth-3*convgauss, hwidth+3*convgauss+step,
                    step)
                cf = _phc.normgauss(convgauss, x=xn)
                yo = _np.interp(xn, x, yo)
                x = xn
                y1 = yo * (1 - frac)
                y2 = _np.convolve(yo * frac, cf / _np.trapz(cf), 'same')
                ax2.plot(x, y1, color=_phc.colors[_np.mod(i, 
                    len(_phc.colors))])
                ax2.plot(x, y2, color=_phc.colors[_np.mod(i, 
                    len(_phc.colors))])
            y = y1 + y2
            # y = linfit(x, y1+y2)
            if file == fullseds[0]:
                ax.plot(x, y, label='{0:02.1f} deg. {1}'.format(obsdegs[k], 
                labels[i]), color=_phc.colors[_np.mod(i, len(_phc.colors))])
                # ew0 = EWcalc(x, y, vw=hwidth)
            else:
                ax.plot(x, y, color=_phc.colors[_np.mod(i, len(_phc.colors))],
                    label=labels[i])
            # ewf = EWcalc(x, y, vw=hwidth)

        plot_line_str(fig, ax, lbc=lbc, ylim=ylim, cmapn=cmapn, xlim=[-hwidth, 
            hwidth])
        figname = outpath + 'modsover_lbc{1:.4f}_obs{0:02.1f}{2}'.\
            format(obsdegs[k], lbc, addsuf)
        _phc.savefig(fig, figname, fmt)

        plot_line_str(fig2, ax2, lbc=lbc, ylim=ylim, cmapn=cmapn, 
            xlim=[-hwidth, hwidth])
        figname = outpath + 'modsover_lbc{1:.4f}_obs{0:02.1f}{2}Extra'.\
            format(obsdegs[k], lbc, addsuf)
        _phc.savefig(fig, figname, fmt)

    return


def overPlotLineFits(specs, lbc=.6564606, fmt=['png'], hwidth=1500., 
    ylim=None, yzero=False, addsuf='', dlim=None, cmapn='jet', xlim=None,
    outpath=''):
    """Generate overplot spec. line from a FITS file list.
    """
    fig, ax = _plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(_phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x, y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        if dlim is None:
            cor = _phc.colors[_np.mod(i, len(_phc.colors))]
        else:
            cor = _phc.gradColor([MJD], min=dlim[0], max=dlim[1], 
                cmapn=cmapn)[0]
        ax.plot(x, y, label='{0}'.format(dateobs), color=cor)

    ylabel = 'Overplotted spectra'
    fig, ax = plot_line_str(fig, ax, lbc=lbc, ylabel=ylabel, xlim=xlim, 
        dlim=dlim, cmapn=cmapn, ylim=ylim)

    figname = outpath + 'fitsover_lbc{1:.4f}{0}'.format(addsuf, lbc)
    _phc.savefig(fig, figname, fmt)
    return


def incrPlotLineSeries(fullseds, obsers=[0], lbc=.6564606, fmt=['png'], 
    addsuf='', outpath=''):
    """Generate incremented spec. line from a HDUST mod list, separated by
    observers. The increment is 0.1 for each file in fullseds sequence.

    Observers config. must be the same between models in `fullseds` list.
    """
    for obs in obsers:
        fig, ax = _plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = _hdt.readfullsed2(file)
            obsdegs = (_np.arccos(sed2data[:, 0, 0]) * 180 / _np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x, y) = lineProf(sed2data[obs, :, 2], sed2data[obs, :, 3], 
                lbc=lbc)
            if file == fullseds[0]:
                ax.plot(x, y + 0.1 * i, label='{0:02.1f} deg.'.format(
                    obsdegs[k]), color=_phc.colors[_np.mod(i, 
                    len(_phc.colors))])
            else:
                ax.plot(x, y + 0.1 * i, color=_phc.colors[_np.mod(i, 
                    len(_phc.colors))])
        ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
        ax.legend(loc='best', fancybox=True, framealpha=0.5)
        figname = outpath + 'modsincr_lbc{1:.4f}_obs{0:02.1f}{2}'.\
            format(obsdegs[k], lbc, addsuf)
        for f in fmt:
            print('# Saved {1}.{0}'.format(f, figname))
            fig.savefig(figname + '.{0}'.format(f), transparent=True)
        _plt.close()
    return


def incrPlotLineFits(specs, lbc=.6564606, fmt=['png'], hwidth=1500., 
    yzero=False, addsuf='', dlim=None, cmapn='jet', xlim=None, outpath='',
    ylim=None):
    """Generate incremented spec. line from FITS files list.
    The increment is 0.1 for each file in sequence.
    """
    fig, ax = _plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(_phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x, y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        if dlim is None:
            cor = _phc.colors[_np.mod(i, len(_phc.colors))]
        else:
            cor = _phc.gradColor([MJD], min=dlim[0], max=dlim[1], 
                cmapn=cmapn)[0]
        ax.plot(x, y + 0.1 * i, label='{0}'.format(dateobs), color=cor)
    if yzero:
        ylim = ax.get_ylim()
        ax.plot([0, 0], ylim, ls='-', color='Gray')

    ylabel = 'Spaced spectra'
    fig, ax = plot_line_str(fig, ax, lbc=lbc, ylabel=ylabel, xlim=xlim, 
        dlim=dlim, cmapn=cmapn, ylim=ylim)

    figname = outpath + 'fitsincr_lbc{1:.4f}{0}'.format(addsuf, lbc)
    _phc.savefig(fig, figname, fmt)
    return


def diffPlotLineSeries(fullseds, obsers=[0], lbc=.6564606, fmt=['png'], 
    rvel=None, rflx=None, hwidth=1000., outpath='', addsuf=''):
    """Generate overplot of DIFFERENCE spec. line from a HDUST mod list.
    The model will be linearly interpolated
    with the reference spec. If none is given as reference, 
    then it assumes the first of the list.

    It is recommend to run first (rvel, rflx) =  lineProf(rvel, rflx,
    lbc=lbc, hwidth=hwidth).

    Observers config. must be the same between models in
    `fullseds` list.
    """
    for obs in obsers:
        fig, ax = _plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = _hdt.readfullsed2(file)
            obsdegs = (_np.arccos(sed2data[:, 0, 0]) * 180 / _np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x, y) = lineProf(sed2data[obs, :, 2], sed2data[obs, :, 3], 
                lbc=lbc, width=hwidth)
            if rvel is None or rflx is None:
                refspec = _hdt.readfullsed2(fullseds[0])
                (vel, flx) = lineProf(refspec[obs, :, 2], refspec[obs, :, 3], 
                    lbc=lbc, hwidth=hwidth)
            else:
                flx = _np.interp(x, rvel, rflx)
            if file == fullseds[0]:
                ax.plot(x, y - flx, label='{0:02.1f} deg.'.format(obsdegs[k]), 
                    color=_phc.colors[_np.mod(i, len(_phc.colors))])
            else:
                ax.plot(x, y - flx, color=_phc.colors[_np.mod(i, 
                    len(_phc.colors))])
        ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
        ax.set_ylabel('Difference spectra (spec - ref.)')
        ax.legend(fontsize=8, loc='best', fancybox=True, framealpha=0.5)
        figname = outpath + 'modsdiff_lbc{1:.4f}_obs{0:02.1f}{2}'.\
            format(obsdegs[k], lbc, addsuf)
        for f in fmt:
            print('# Saved {1}.{0}'.format(f, figname))
        _plt.close()
    return


def diffPlotLineFits(specs, lbc=.6564606, fmt=['png'], xlim=None,
    rvel=None, rflx=None, hwidth=1500., addsuf='', cmapn='jet', dlim=None,
    outpath='', ylim=None):
    """Generate overplot of DIFFERENCE spec. line from a FITS files list.
    The observations will be linearly interpolated
    with the reference spec. If none is given as reference, 
    then it assumes the first of the list.

    It is recommend to run first (rvel, rflx) =  lineProf(rvel, rflx,
    lbc=lbc, hwidth=hwidth).

    If `cmap` is None or empty, the phc.colors vector is read.
    """
    fig, ax = _plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(_phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x, y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if rvel is None or rflx is None:
            # wl0, flux0, MJD, dateobs0, datereduc, fitsfile = \
            #   loadfits(specs[0])
            # (rvel,flx) = lineProf(wl0, flux0, lbc=lbc, hwidth=hwidth)
            # flx = _np.interp(x, rvel, rflx)
            rvel = x
            rflx = y
            flx = y[:]
        else:
            flx = _np.interp(x, rvel, rflx)
        # if spec == specs[0]:
            # ax.plot(x, y-flx, label='{0}'.format(dateobs), \
            # color= _phc.colors[_np.mod(i, len(_phc.colors))])
        # else:
            # ax.plot(x, y-flx, color= _phc.colors[_np.mod(i, 
                # len(_phc.colors))])
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        if dlim is None:
            cor = _phc.colors[_np.mod(i, len(_phc.colors))]
        else:
            cor = _phc.gradColor([MJD], min=dlim[0], max=dlim[1], 
                cmapn=cmapn)[0]
        ax.plot(x, y - flx, label='{0}'.format(dateobs), color=cor)

    ylabel = 'Difference spectra'
    fig, ax = plot_line_str(fig, ax, lbc=lbc, ylabel=ylabel, xlim=xlim, 
        dlim=dlim, cmapn=cmapn, ylim=ylim)

    figname = outpath + 'fitsdiff_lbc{1:.4f}{0}'.format(addsuf, lbc)
    _phc.savefig(fig, figname, fmt)
    return


def diffPlotLineObs(fullseds, obsers=[0], lbc=.6564606, fmt=['png'], 
    rvel=None, rflx=None, hwidth=1000., addsuf='', outpath=''):
    """Generate overplot of DIFFERENCE spec. line from a HDUST OBSERVERS list.
    The model will be linearly interpolated
    with the reference spec. If none is given as reference, 
    then it assumes the first observer of the list.

    It is recommend to run first (rvel, rflx) =  lineProf(rvel, rflx,
    lbc=lbc, hwidth=hwidth).

    Observers config. must be the same between models in
    `fullseds` list.
    """
    for file in fullseds:
        fig, ax = _plt.subplots()
        sed2data = _hdt.readfullsed2(file)
        obsdegs = (_np.arccos(sed2data[:, 0, 0]) * 180 / _np.pi)[obsers]
        obsdegs = list(obsdegs)
        for obs in obsers:
            i = obsers.index(obs)
            (x, y) = lineProf(sed2data[obs, :, 2], sed2data[obs, :, 3], 
                lbc=lbc, hwidth=hwidth)
            if rvel is None or rflx is None:
                (vel, flx) = lineProf(sed2data[obsers[0], :, 2], 
                sed2data[obsers[0], :, 3], lbc=lbc, hwidth=hwidth)
            else:
                flx = _np.interp(x, rvel, rflx)
            ax.plot(x, y - flx, label='{0:02.1f} deg.'.format(obsdegs[i]), 
                color=_phc.colors[_np.mod(i, len(_phc.colors))])
        ax.set_title(u'lbc={0:.5f}$\mu$m, {1}'.format(lbc, 
            _phc.trimpathname(file)[1]))
        ax.set_ylabel('Difference spectra (spec - ref.)')
        ax.legend(fontsize=8, loc='best', fancybox=True, framealpha=0.5)
        figname = outpath + 'modsdiff_lbc{1:.4f}{0}'.format(addsuf, lbc)
        for f in fmt:
            print('# Saved {1}.{0}'.format(f, figname))
            fig.savefig(figname + '.{0}'.format(f), transparent=True)
        _plt.close()
    return


def max_func_pts(x, y, ws=0.01, avgbox=3):
    """ `ws` window size where the maximum will be evaluated. Example: `ws=0.02`
    corresponds to 2% of the length of the input. """
    x, y = (_np.array(x), _np.array(y))
    N = len(x)
    parts = _phc.splitequal(N*ws, N)
    n = len(parts)
    xout, yout = (_np.zeros(n), _np.zeros(n))
    for i in range(n):
        p = parts[i]
        Y = y[p[0]:p[1]]
        X = x[p[0]:p[1]]
        idx = _np.argsort(Y)
        xout[i] = _np.average(X[idx][-avgbox:])
        yout[i] = _np.average(Y[idx][-avgbox:])
    return xout, yout


def sum_ec(fwl, fflx):
    dmin = _np.inf
    wlmin = _np.inf
    wlmax = 0
    for f in fwl:
        if _np.min(_np.diff(f)) < dmin:
            dmin = _np.min(_np.diff(f))
        if _np.min(f) < wlmin:
            wlmin = _np.min(f)
        if _np.max(f) > wlmax:
            wlmax = _np.max(f)
    swl = _np.arange(wlmin, wlmax, dmin)
    sflx = _np.zeros(len(swl))
    for i in range(len(fwl)):
        idx = _np.where( (swl > _np.min(fwl[i])) & (swl < _np.max(fwl[i])) )
        sflx[idx] += _np.interp(swl[idx], fwl[i], fflx[i])
    return swl, sflx


def lbdc2range(lbdc):
    """ Function doc

    """
    dl = lbdc[1] - lbdc[0]
    return _np.linspace(lbdc[0] - dl / 2, lbdc[-1] + dl / 2, len(lbdc) + 1)

# MAIN ###
if __name__ == "__main__":
    pass
