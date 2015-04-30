#-*- coding:utf-8 -*-

"""
PyHdust *spectools* module: spectroscopy tools

Algumas definicoes: para todas as rotinas funcionarem, todos os espectros devem
estar agrupados num mesmo caminho (`path`), em estrutura de
noite/estrelas/espec.

Neste momento, as rotinas somente leem arquivos `*.cal.fits`. Para receber este
sufixo `.cal`, algumas informacoes no header sao necessarias:

    * 'MJD-OBS' ou 'MJD' ou 'JD' ou 'DATE-OBS'
    * 'CRVAL1' + 'CDELT1'

:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""

import os as _os
import numpy as _np
import datetime as _dt
from glob import glob as _glob
import pyhdust.phc as _phc
import pyhdust.jdcal as _jdcal
import pyhdust as _hdt

try:
    import pyfits as _pyfits
    import matplotlib.pyplot as _plt
    from scipy.interpolate import interp1d as _interp1d
except:
    print('# Warning! matplotlib, scipy and/or pyfits module not installed!!!')

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
    >>> #significa que vetor wl eh vetor velocidades, e nao comprimento de onda.
    >>> spdtb.lbc = 6564.
    >>> spdtb2 = Spec()
    >>> spdtb2.lbc = 4863.

    Como usar (hard way):

    >>> spdtb = spt.Spec()
    >>> #read spec `flux` and `wl` for a given `lbc`
    >>> (spdtb.EW, spdtb.EC, spdtb.VR, spdtb.peaksep, spdtb.depthcent,\\ 
    >>> spdtb.F0) = spt.analline(wl, flux, lbc)
    >>> spdtb.MJD = 1
    >>> spdtb.file = file

    And then:

    >>> #to record it to the database:
    >>> spdtb.addspec()

    Para carregar uma tabela anterior, faca:

    >>> spdtb = spt.Spec()
    >>> #(...) read new specs and then join with previous ones
    >>> spdtb.data = _np.vstack(( spdtb.data, _np.loadtxt('hdt/datafile.txt') ))
    >>> spdtb.metadata = _np.vstack(( spdtb.metadata, \\
    >>> _np.loadtxt('hdt/metafile.txt') ))
    >>> spdtb.updatecount() #to update the counter

    Ou simplesmente (nome de arquivos default):

    >>> spdtb.loaddata()
    """
    def __init__(self, wl=None, flux=None, lbc=None, hwidth=None, EW=_np.NaN,
        EC=_np.NaN, VR=_np.NaN,
        peaksep=_np.NaN, depthcent=_np.NaN, F0=_np.NaN, dateobs='', MJD=0.,
        datereduc='', file=''):
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
        self.MJD= MJD
        self.count = 0
        self.data = _np.empty(0)
        self.metadata = _np.empty(0)

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
        self.MJD= 0.        

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
        #if self.flux != None and self.wl != None and self.lbc != None:
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

    def savedata(self, datafile=_outfold+'/datafile.txt',\
        metafile=_outfold+'/metafile.txt'):
        """Save current table
        """
        header = ['MJD', 'EW', 'EC', 'VR', 'peaksep', 'depthcent', 'F0']
        _np.savetxt(datafile, self.data, fmt='%12.6f',
        header=(len(header)*'{:>12s}').format(*header))
        _np.savetxt(metafile, self.metadata, fmt='%s', delimiter=',')
        return

    def loaddata(self,  datafile=_outfold+'/datafile.txt',\
        metafile=_outfold+'/metafile.txt'):
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
        if file.find('.fit') ==-1:
            print("# ERROR! `loadspec` unrecognized format!")
            return
        (self.wl, self.flux, self.MJD, self.dateobs, self.datereduc, self.file)\
        = loadfits(file)
        (self.EW, self.EC, self.VR, self.peaksep, self.depthcent, self.F0) = \
        analline(self.wl, self.flux, self.lbc, verb=False)
        return

    def plotspec(self, outname=''):
        """Export current spec into a PNG file.
        """
        if self.wl == None or self.flux == None:
            print('# ERROR: wrong Spec() parameters! {}'.format(self.file))
            return
        if outname == '':
            path, file = _phc.trimpathname(self.file)
            outname = _phc.rmext(file)
        #Normalization:
        flux = linfit(self.wl, self.flux)
        wl = self.wl
        fig = _plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(wl, flux)
        ax.set_ylabel('norm. flux')
        ax.set_xlabel('wavelength (arb. units)')
        ax.set_title(outname)
        _plt.savefig('{}/{:.2f}_{}.png'.format(_outfold,self.MJD,outname))
        if self.lbc > 0:
            vels = (self.wl-self.lbc)/self.lbc*_phc.c.cgs*1e-5
            idx = _np.where(_np.abs(vels) <= self.hwidth)
            flux = linfit(vels[idx], flux[idx])
            vels = vels[idx]
            _plt.clf()
            ax = fig.add_subplot(111)
            ax.plot(vels, flux)
            ax.set_ylabel('norm. flux')
            ax.set_xlabel('vel. (km/s)')
            ax.set_title('{:.2f} {} {:.2f}'.format(self.MJD,outname,self.lbc))
            _plt.savefig('{}/{:.2f}_{}_{:.2f}.png'.format(_outfold,self.MJD,\
            outname,self.lbc))
        _plt.close()
        return


def extractfromsplot(file, splot):
    """Ce = center; Co = core
    #LcCe, LcCo, lcGW, lcEW, lvCe, lcCo, lvEW, lrCe, LrCo, lrEW
    """
    out = _np.array(10*[_np.NaN])
    readflag = False
    for line in splot:
        if line.find(']:') > 0 and readflag:
            readflag = False
        if line.find(file) > 0:
            readflag = True
        if readflag:
            info = line.split()
            #if _re.match("^\d+?\.\d+?$", info[0]) is not None:
            try:
                float(info[0])
                info = _np.array(info, dtype=float)
                if info[0] > 6556 and info[0] < 6556+4.33:
                    if len(info) == 4:
                        out[6] = float(info[3])
                    elif len(info) == 7:
                        out[4] = float(info[0])
                        out[5] = float(info[4])
                elif info[0] > 6556+4.33 and info[0] < 6556+2*4.33:
                    if len(info) == 4:
                        out[3] = float(info[3])
                    elif len(info) == 7:
                        out[0] = float(info[0])
                        out[1] = float(info[4])
                        out[2] = float(info[5])
                elif info[0] > 6556+2*4.33 and info[0] < 6556+3*4.33:
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
        dtobs = dtobs.replace('.','')
        tobs, dtobs = dtobs.split('T')
        if len(tobs) == 10:
            dtobs, tobs = tobs, dtobs
        tobs = tobs.split(':')
        tobs = float(tobs[0])*3600+float(tobs[1])*60+float(tobs[2])
        tobs /= (24*3600)
    else:
        tobs = 0.
    if dtobs[4] == '-':
        dtobs = dtobs.split('-')
    elif dtobs[2] == '/':
        dtobs = dtobs.split('/')[::-1]
    else:
        print('# ERROR! Wrong "DATE-OBS" in header! {}'.format(fitsfile))
        raise systemExit(1)
    dtobs = _np.array(dtobs, dtype='int32')
    return dtobs, tobs


def shiftfits(fitsfile, newsh=None):
    """ Update FITS spec header for a given shift value. """
    imfits = _pyfits.open(fitsfile, mode='update')
    if 'WLSHIFT' in imfits[0].header:
        print('# WLSHIFT = {0} for {1}'.format(imfits[0].header['WLSHIFT'], \
        _phc.trimpathname(fitsfile)[1]))
    else:
        print('# No WLSHIFT available for {0}'.format( \
        _phc.trimpathname(fitsfile)[1]))
    if newsh is None:
        newsh = raw_input('Type the new shift: ')
    if newsh != '':
        imfits[0].header['WLSHIFT'] = float(newsh)
        imfits.close()
    return


def loadfits(fitsfile):
    """load FITS spec

    Out: wl, flux, MJD, dateobs, datereduc, fitsfile
    """
    imfits = _pyfits.open(fitsfile)
    flux = imfits[0].data
    wl = _np.arange(len(flux))*imfits[0].header['CDELT1']+\
        imfits[0].header['CRVAL1']
    (MJD, dateobs, datereduc) = (0., '', '')
    if 'MJD-OBS' in imfits[0].header:
        MJD = float(imfits[0].header['MJD-OBS'])
    elif 'MJD' in imfits[0].header:
        MJD = float(imfits[0].header['MJD'])
    elif 'JD' in imfits[0].header:
        MJD = float(imfits[0].header['JD'])-2400000.5
    elif 'DATE-OBS' in imfits[0].header:
        dtobs = imfits[0].header['DATE-OBS']
        dtobs,tobs = check_dtobs(dtobs)
        MJD = _jdcal.gcal2jd(*dtobs)[1]+tobs
    elif 'FRAME' in imfits[0].header:
        dtobs = imfits[0].header['FRAME']
        dtobs,tobs = check_dtobs(dtobs)
        MJD = _jdcal.gcal2jd(*dtobs)[1]+tobs
    else:
        MJD = _jdcal.MJD_JD2000
        print('# Warning! No DATE-OBS information is available! {}'.\
        format(fitsfile))
        print('# Assuming MJD_JD2000')
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
    VAC = AIR / (1.0 - 2.73443407E-4 - 1.31275255E2 / AIR^2 - 2.75708212E8 / AIR^4 )
    """
    return wl / (1.0 -2.73443407e-04 -1.31275255e+02 / wl**2\
    -2.75708212e+08 / wl**4)


def checksubdirs(path, star, lbc, hwidth=1000, showleg=True, plots=False):
    """
    Faz o que tem que fazer.
    """
    if _os.path.exists('{}/{}'.format(path,star)) == False:
        _os.system('mkdir {}/{}'.format(path,star))

    nights = [o for o in _os.listdir(path) if _os.path.isdir('{}/{}'.format(path,\
    o))]

    fig = _plt.figure()
    ax = fig.add_subplot(111)
    spdtb = Spec()
    spdtb.lbc = lbc
    spdtb.hwidth = 1000.
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path,night)) if
            _os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:
                    for cal in scal:
                        spdtb.loadspec(cal)
                        spdtb.addspec()
                        if _np.isnan(spdtb.EW) == False:
                            if plots:
                                spdtb.plotspec()
                            vels = (spdtb.wl-lbc)/lbc*_phc.c.cgs*1e-5
                            idx = _np.where(_np.abs(vels) <= hwidth)
                            flux = linfit(vels[idx], spdtb.flux[idx])
                            vels = vels[idx]
                            leg = spdtb.MJD
                            ax.plot(vels, flux, label=leg, alpha=0.7, color=\
                            _phc.colors[_np.mod(spdtb.count, len(_phc.colors))])
                else:   
                    print('# Data not reduced for %s at %s!' % (star,night))
    ax.set_xlim([-hwidth,hwidth])
    ax.set_ylim([-1,5])
    if showleg:
        legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
        _plt.setp(legend.get_texts(),  fontsize='small')    
    _plt.savefig('{}/{}_at_{}.png'.format(_outfold,star,lbc))
    _plt.close()
    spdtb.savedata(datafile='{}/{}.txt'.format(_outfold,star),\
    metafile='{}/meta_{}.txt'.format(_outfold,star))
    return


def hydrogenlinewl(ni, nf):
    """Generate H line transitions wavelengths in meters for VACUUM

    Rydberg constant `R` was manually adjusted to fit Halpha and Hbeta lines.
    """
    return (10967850.*(1./nf**2-1./ni**2))**-1.


def calcres_R(hwidth=1350, nbins=108):
    """
    (h)Width in km/s.
    *WARNING*: `width` in HDUST input is only half.
    To HDUST effective R, multiple the input width by 2 he_re.

    # R = lbd/Dlbd = _phc.c/Dv = _phc.c*nbins/width 
    # nbins = R*width/_phc.c
    """
    return round(_phc.c.cgs*nbins/hwidth/1e5)


def calcres_nbins(R=12000, hwidth=1350):
    """
    (h)Width in km/s.
    *WARNING*: `width` in HDUST input is only half.
    To HDUST effective R, multiple the input width by 2 he_re.

    # R = lbd/Dlbd = _phc.c/Dv = _phc.c*nbins/width 
    # nbins = R*width/_phc.c
    """
    return round(R*hwidth*1e5/_phc.c.cgs)


def lineProf(x, flx, lbc, flxerr=_np.empty(0), hwidth=1000., ssize=0.05):
    '''
    lineProf() - retorna um array (flx) normalizado e um array x em VELOCIDADES.
    `lbc` deve fornecido em mesma unidade de x para conversão lambda -> vel.
    Se vetor x jah esta em vel., usar funcao linfit().

    x eh importante, pois y pode ser nao igualmente amostrado.
    x e y devem estar em ordem crescente.

    ssize = % do tamanho de y; numero de pontos usados nas extremidades
    para a media do contínuo. 'ssize' de .5 à 0 (exclusive).

    OUTPUT: vel (array), flx (array)
    '''    
    x = (x-lbc)/lbc*_phc.c.cgs*1e-5 #km/s
    idx = _np.where(_np.abs(x) <= hwidth)
    flux = linfit(x[idx], flx[idx], yerr=flxerr, ssize=ssize)
    return x[idx], flux


def linfit(x, y, ssize=0.05, yerr=_np.empty(0)):
    '''
    linfit() - retorna um array (y) normalizado, em posicoes de x

    x eh importante, pois y pode ser nao igualmente amostrado.
    x e y devem estar em ordem crescente.

    ssize = % do tamanho de y; numero de pontos usados nas extremidades
    para a media do contínuo. 'ssize' de .5 à 0 (exclusive).

    OUTPUT: y, yerr (if given)
    '''
    if ssize < 0 or ssize > .5:
        print('# Invalid ssize value...')
        raise SystemExit(1)
    ssize = int(ssize*len(y))
    if ssize == 0:
        ssize = 1
    medx0, medx1 = _np.average(x[:ssize]),_np.average(x[-ssize:])
    if ssize > 20:
        medy0, medy1 = _np.median(y[:ssize]),_np.median(y[-ssize:])
    else:
        medy0, medy1 = _np.average(y[:ssize]),_np.average(y[-ssize:])
    new_y = medy0 + (medy1 - medy0) * (x - medx0) / (medx1 - medx0)
    idx = _np.where(new_y != 0)
    y[idx] = y[idx]/new_y[idx]
    if len(yerr) == 0.:
        return y
    else:
        yerr = yerr/_np.average(new_y)
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
        #normflux = _np.ones(len(outvels))
        return ew
    for i in range(len(outvels)-2):
        dl = outvels[i+1]-outvels[i]
        ew += (1.-(normflux[i+1]+normflux[i])/2.)*dl
    return ew


def ECcalc(vels, flux):
    """
    Supoe que o fluxo jah estah normalizado, e vetores ordenad_os.

    Calcula o topo da emissao da linha, e retorna em que velocidade ela
    ocorre.
    """
    idx = _np.where(_np.max(flux) == flux)
    return flux[idx], vels[idx]


def VRcalc(vels, flux, vw=1000):
    """
    Supoe que o fluxo jah estah normalizado, e vetores ordenad_os.

    Calcula o ew para os dois lados (azul/vermelho) da linha, ajustando
    a velocidade de repouso. 
    """
    #calcula e aplica correcao de vel. repousp
    vc = 0.
    vels += vc
    #corta em vw, e faz o teste de tamanho
    if vw > 0:
        idx = _np.where(_np.abs(vels) <= vw)
        outvels = vels[idx]
        normflux = flux[idx]
    if len(vels) < 3:
        ew0, ew1 = 0.
        return ew0, ew1, vc
    #
    ivc = _np.abs(outvels-0).argmin()
    ew0 = 0.
    for i in range(0,ivc+1-2):
        dl = outvels[i+1]-outvels[i]
        ew0 += (1.-(normflux[i+1]+normflux[i])/2.)*dl
    ew1 = 0.
    for i in range(ivc,len(outvels)-2):
        dl = outvels[i+1]-outvels[i]
        ew1 += (1.-(normflux[i+1]+normflux[i])/2.)*dl
    return ew0, ew1, vc


def PScalc(vels, flux, vc=0., ssize=.05):
    """
    Calcula peak_separation
    """
    #check if there is a peak
    ssize = int(ssize*len(vels))
    if ssize == 0:
        ssize = 1
    contmax = _np.max(_np.append(flux[:ssize],flux[-ssize:]))
    fluxmax = _np.max(flux)
    if fluxmax < 1.01*contmax:
        return 0, 0
    vels += vc
    ivc = _np.abs(vels-0).argmin()
    i0 = _np.abs(flux[:ivc]-_np.max(flux[:ivc])).argmin()
    i1 = _np.abs(flux[ivc+1:]-_np.max(flux[ivc+1:])).argmin()+ivc+1
    return vels[i0], vels[i1]


def DCcalc(vels, flux, vmax=None, vc=0., ssize=0.05):
    """
    Calculo, na presenca de emissao, da profundidade do reverso central.
    
    """
    vels += vc
    ivc = _np.abs(vels-0).argmin()
    #check if there is a peak
    ssize = int(ssize*len(vels))
    if ssize == 0:
        ssize = 1
    contmax = _np.max(_np.append(flux[:ssize],flux[-ssize:]))
    fluxmax = _np.max(flux)
    if fluxmax < 1.01*contmax:
        return flux[ivc], flux[ivc]
    #if a vmax is not given...
    if isinstance(vmax, (int, long, float)) == False:
        vmax = _np.abs(flux-_np.max(flux)).argmin()
        vmax = vels[vmax]
    ivmax = _np.abs(vels-vmax).argmin()
    return flux[ivmax], flux[ivc]


def analline(lbd, flux, lbdc, hwidth=1000, verb=True):
    """
    Return the analysis of a line.

    Both lbd and flux need to be ordered (a normalization IS FORCED).
    lbd,lbdc must have the same unit, and width in km/s is required.
    The line will be cutted so that the total DeltaLambda will be 2*width

    if lbdc <= 0, lbd array is assumed to be a velocity array (in km/s)!

    | EXAMPLE: Using sed2data. lbc = 0.6565 (halpha), obs = 1 (width==1000)
    |     analline(lbd=sed2data[obs,:,2], flux=sed2data[obs,:,3], lbc=lbc)

    OUTPUT: EW, EC, VR, peaksep, depthcent, F0 
    """
    if lbdc > 0:
        vels = (lbd-lbdc)/lbdc*_phc.c.cgs*1e-5
    else:
        vels = lbd
    #check if the file have the desired info.
    if vels[0] > -hwidth or vels[-1] < hwidth:
        if verb:
            print('# ERROR: spec out of range (wavelength)! Check hwidth!')
        return _np.NaN, _np.NaN, _np.NaN, _np.NaN, _np.NaN, _np.NaN

    idx = _np.where(_np.abs(vels) <= hwidth)
    vels = vels[idx]
    flux = flux[idx]
    #Normalization:
    flux = linfit(vels, flux)
    #Output:
    EW = EWcalc(vels, flux, vw=hwidth)
    EC, velEC = ECcalc(vels, flux)
    ew0, ew1, vc = VRcalc(vels, flux, vw=hwidth)
    if ew1 == 0:
        VR = 1
    else:
        VR = ew0/ew1
    vel0, vel1 = PScalc(vels, flux)
    peaksep = vel1-vel0
    EC2, F0 = DCcalc(vels, flux, vmax=velEC)
    depthcent = EC2-F0
    
    return EW, EC, VR, peaksep, depthcent, F0


def normalize_range(lb,spec,a,b):
    """This function is obsolete and must be removed.

    Still here for compatibility issues.
    """
    a2 = (spec[b]-spec[a])/(lb[b]-lb[a])
    a1 = spec[a]-a2*lb[a]
    return spec/(a1+a2*lb)


def overplotsubdirs(path, star, limits=(6540,6600), showleg=True):
    """
    Realiza o plot de espectros da estrela `star` dentre do diretorio `path`.
    Atualmente, faz o plot entre os valores `limits` (Angstroms).

    Gera os arquivos `path/star/star.log` e `path/star/star_specs.png`.
    """
    #path = _os.getcwd()
    #star = raw_input('Type the star name: ')
    #ref0 = 6540
    #ref1 = 6600

    ref0,ref1 = limits
    
    if _os.path.exists('{}/{}'.format(path,star)) == False:
        _os.system('mkdir {}/{}'.format(path,star))
    f0 = open('{0}/{1}/{1}.log'.format(path,star), 'w')
    
    nights = [o for o in _os.listdir(path) if _os.path.isdir('{}/{}'.format(path,\
        o))]
    
    legendl = ()
    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path,night)) if
            _os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = _glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    #legendl += (night,)
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec))*imfits[0].header['CDELT1']+\
                            imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[-1] > 6560: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = _np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = _np.where(abs(lbda-ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda,spec,a0,a1)
                            msg = '{}, {}, {}'.format((0.1*i), night, cal)
                            print(msg)
                            f0.writelines(msg+'\n')
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            _plt.plot(lbda, spec, label=leg, alpha=0.7,
                                color=_phc.colors[_np.mod(i, len(_phc.colors))])
                            i+=1
                else:
                    print('# Data not reduced for %s at %s!' % (star,night))
                    msg =  '{}, {}, {}'.format('NC', night, 'None')
                    f0.writelines(msg+'\n')

    if showleg:
        legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
        _plt.setp(legend.get_texts(),  fontsize='small')
    _plt.xlim([ref0,ref1])
    _plt.ylim([-1, 5])
    #_plt.xlabel('vel. (km/s)')
    _plt.savefig('{0}/{1}/{1}_specs.png'.format(path,star))
    _plt.close()
    f0.close()
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = _pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = _np.arange(len(specs[i]))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
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
    # #legend(legendl, 'lower right')
    # #legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #'lower right'
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


def diffplotsubdirs(path, star, limits=(6540,6600)):
    """
    Realiza o plot de espectros da estrela `star` dentre do diretorio `path`.
    Atualmente, faz o plot entre os valores `limits` (Angstroms).

    Gera os arquivos `path/star/star.log` e `path/star/star_specs_dif.png`.
    """
    ref0,ref1 = limits
    
    if _os.path.exists('{}/{}'.format(path,star)) == False:
        _os.system('mkdir {}/{}'.format(path,star))
    f0 = open('{0}/{1}/{1}.log'.format(path,star), 'w')
    
    nights = [o for o in _os.listdir(path) if _os.path.isdir('{}/{}'.format(path,\
        o))]    

    legendl = ()
    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path,night)) if _os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = _glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        print('# Specs with dopcor at %s!' % night)
                        srv = scal
                    #legendl += (night,)
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[0] > 5500: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = _np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = _np.where(abs(lbda-ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda,spec,a0,a1)+(0.1*i)
                            print (0.1*i)
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            _plt.plot([ref0,ref1], [1+0.1*i,1+0.1*i], 'k--', alpha=0.5)
                            _plt.plot(lbda, spec, label=leg, color=_phc.colors[i])
                            i+=1
                else:
                    print('# Data not reduced for %s at %s!' % (star,night))
    
    legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
    _plt.setp(legend.get_texts(),  fontsize='small')
    _plt.xlim([ref0,ref1])
    #_plt.xlabel('vel. (km/s)')
    _plt.savefig('{0}/{1}/{1}_specs_dif.png'.format(path,star))
    
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = _pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = _np.arange(len(specs[i]))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
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
    # #legend(legendl, 'lower right')
    # #legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #'lower right'
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


def refplotsubdirs(path, star, limits=(6540,6600)):
    """
    Realiza o plot de espectros da estrela `star` dentre do diretorio `path`.
    Atualmente, faz o plot entre os valores `limits` (Angstroms).

    Gera os arquivos `path/star/star.log` e `path/star/star_specs_REFERENCIA.png`.
    """
    ref0,ref1 = limits
    
    if _os.path.exists('{}/{}'.format(path,star)) == False:
        _os.system('mkdir {}/{}'.format(path,star))
    f0 = open('{0}/{1}/{1}.log'.format(path,star), 'w')
    
    nights = [o for o in _os.listdir(path) if _os.path.isdir('{}/{}'.format(path,\
        o))]

    legendl = ()
    i = 0
    
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path,night)) if _os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = _glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        srv = scal
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[0] > 5500: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = _np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = _np.where(abs(lbda-ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda,spec,a0,a1)
                            print (0.1*i)
                            leg = imfits[0].header['DATE-OBS']
                            refleg = '2012-11-20T23:51:37.392'
                            refleg = '2008-06-13'
                            if leg == refleg:
                                f0=open('{0}/{1}/ref.txt'.format(path,star),'w')
                                f0.writelines([str(x)+'\t' for x in lbda])
                                f0.writelines('\n')
                                f0.writelines([str(x)+'\t' for x in spec])
                                f0.writelines('\n')
                                f0.close()
                            i+=1
                else:
                    print('# Data not reduced for %s at %s!' % (star,night))
    
    
    f0=open('{0}/{1}/ref.txt'.format(path,star))
    lines = f0.readlines()
    f0.close()
    specref = _np.array(lines[1].split(), dtype=float)
    lbdaref = _np.array(lines[0].split(), dtype=float)
    func = interp1d(lbdaref, specref)#, kind='cubic')
    lbdaref = _np.linspace(ref0,ref1,5000)
    specref = func(lbdaref)
    
    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path,night)) if _os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = _glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    #legendl += (night,)
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[0] > 5500: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = _np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = _np.where(abs(lbda-ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda,spec,a0,a1)
                            func = interp1d(lbda, spec)#, kind='cubic')
                            #Tive problemas de 'out-of-bounds'... um espectro estava desordenado:
                            #print imfits[0].header['CDELT1'], imfits[0].header['CRVAL1'], cal
                            spec = func(lbdaref)
                            print (0.1*i)
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            if i <130:
                                _plt.plot(lbdaref, spec-specref, label=leg, alpha=0.8, color=_phc.colors[i])
                            i+=1
                else:
                    print('# Data not reduced for %s at %s!' % (star,night))
    
    legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
    _plt.setp(legend.get_texts(),  fontsize='small')
    _plt.xlim([ref0,ref1])
    _plt.title('Ref.= %s' % refleg)
    #_plt.xlabel('vel. (km/s)')
    _plt.savefig('{0}/{1}/{1}_specs_{2}.png'.format(path,star,refleg[:10]))
    
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = _pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = _np.arange(len(specs[i]))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
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
    # #legend(legendl, 'lower right')
    # #legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #'lower right'
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


def overplotsubdirs2(path, star, limits=(6540,6600)):
    """
    Realiza o plot de espectros da estrela `star` dentre do diretorio `path`.
    Atualmente, faz o plot entre os valores `limits` (Angstroms).

    Ha' um criterio de escolha de espectros aqui (rudimentar).

    Gera os arquivos `path/star/star.log` e `path/star/star_specs2.png`.
    """
    ref0,ref1 = limits
    
    if _os.path.exists('{}/{}'.format(path,star)) == False:
        _os.system('mkdir {}/{}'.format(path,star))
    f0 = open('{0}/{1}/{1}.log'.format(path,star), 'w')
    
    nights = [o for o in _os.listdir(path) if _os.path.isdir('{}/{}'.format(path,\
        o))]

    legendl = ()
    ax = _plt.figure()
    i = 0
    for night in nights:
        targets = [o for o in _os.listdir('%s/%s' % (path,night)) if _os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = _glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = _glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    #legendl += (night,)
                    for cal in scal:              
                        imfits = _pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = _np.arange(len(spec))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[0] > 5500: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = _np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = _np.where(abs(lbda-ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda,spec,a0,a1)
                            print (0.1*i), night
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
                            #elif leg.find('03/07/2013') != -1:
                            #    leg = '2013-07-03'
                            #    check = True
                            elif leg.find('28/07/2013') != -1:
                                leg = '2013-07-28'
                                check = True
                            elif leg.find('2013-11-12T01:30:38.938') != -1:
                                leg = '2013-11-12'
                                check = True
                            else:
                                print leg       
                            if check:
                                print cal
                                _plt.plot(lbda, spec, label=leg, alpha=0.7, color=prtcolor)
                            i+=1
                else:
                    msg = '# Data not reduced for %s at %s!' % (star,night)
                    print(msg)
                    f0.writelines(msg)
                    
    
    font = { 'size' : 16, }
    
    legend = _plt.legend(loc=(0.75, .05), labelspacing=0.1)
    #_plt.setp(legend.get_texts(),  fontsize='small')
    _plt.xlim([ref0,ref1])
    _plt.ylim([.58,1.2])
    _plt.xlabel(r'wavelength ($\AA$)', fontdict=font)
    _plt.ylabel('Normalized flux', fontdict=font)
    #_plt.xlabel('vel. (km/s)')
    _plt.savefig('{0}/{1}/{1}_specs2.png'.format(path,star))
    _plt.close()
    f0.close()
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = _pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = _np.arange(len(specs[i]))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
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
    # #legend(legendl, 'lower right')
    # #legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.) #'lower right'
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


def overPlotLineSeries(fullseds, obsers=[0], lbc=.6564606, formats=['png']):
    """Generate overplot spec. line from a HDUST mod list, separated by
    observers.

    Observers config. must be the same between models in `fullseds` list.
    """
    for obs in obsers:
        fig, ax = _plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = hdt.readfullsed2(file)
            obsdegs = (_np.arccos(sed2data[:,0,0])*180/_np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x,y) = lineProf(sed2data[obs,:,2], sed2data[obs,:,3], lbc=lbc)
            if file == fullseds[0]:
                ax.plot(x, y, label='{0:02.1f} deg.'.format(obsdegs[k]), \
                color=_phc.colors[_np.mod(i, len(_phc.colors))])
            else:
                ax.plot(x, y, color=_phc.colors[_np.mod(i, len(_phc.colors))])
        ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
        ax.legend()
        for fmt in formats:
            fig.savefig('modsover_lbc{2:.4f}_obs{0:02.1f}.{1}'.format(obsdegs[k], fmt,\
            lbc), transparent=True)
        _plt.close()
    return

def overPlotLineFits(specs, lbc=.6564606, formats=['png'], hwidth=1500., \
    ylim=None, yzero=False, addsuf=''):
    """Generate overplot spec. line from a FITS file list.
    """
    fig, ax = _plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(_phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x,y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        ax.plot(x, y, label='{0}'.format(dateobs), \
        color=_phc.colors[_np.mod(i, len(_phc.colors))])
    if ylim != None:
        ax.set_ylim(ylim)
    if yzero:
        ylim = ax.get_ylim()
        ax.plot([0,0], ylim, ls='-', color='Gray')
    ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
    ax.set_ylabel('Spec - Ref.')
    legend = ax.legend(loc=(1.05, .01), labelspacing=0.1)
    _plt.setp(legend.get_texts(),  fontsize='small')
    _plt.subplots_adjust(left=0.1, right=0.78, top=0.9, bottom=0.1)#, hspace=0.3, wspace=.3)   
    for fmt in formats:
        fig.savefig('fitsover_lbc{1:.4f}{2}.{0}'.format(\
        fmt, lbc, addsuf), transparent=True)
    _plt.close()
    return


def incrPlotLineSeries(fullseds, obsers=[0], lbc=.6564606, formats=['png']):
    """Generate incremented spec. line from a HDUST mod list, separated by
    observers. The increment is 0.1 for each file in fullseds sequence.

    Observers config. must be the same between models in `fullseds` list.
    """
    for obs in obsers:
        fig, ax = _plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = hdt.readfullsed2(file)
            obsdegs = (_np.arccos(sed2data[:,0,0])*180/_np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x,y) = lineProf(sed2data[obs,:,2], sed2data[obs,:,3], lbc=lbc)
            if file == fullseds[0]:
                ax.plot(x, y+0.1*i, label='{0:02.1f} deg.'.format(obsdegs[k]), \
                color=_phc.colors[_np.mod(i, len(_phc.colors))])
            else:
                ax.plot(x, y+0.1*i, color=_phc.colors[_np.mod(i, len(_phc.colors))])
        ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
        ax.legend()
        for fmt in formats:
            fig.savefig('modsincr_lbc{2:.4f}_obs{0:02.1f}.{1}'.format(obsdegs[k], fmt,\
            lbc), transparent=True)
        _plt.close()
    return


def incrPlotLineFits(specs, lbc=.6564606, formats=['png'], hwidth=1500., \
    yzero=False, addsuf=''):
    """Generate incremented spec. line from FITS files list.
    The increment is 0.1 for each file in sequence.
    """
    fig, ax = _plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(_phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x,y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        ax.plot(x, y+0.1*i, label='{0}'.format(dateobs), \
        color=_phc.colors[_np.mod(i, len(_phc.colors))])
    if yzero:
        ylim = ax.get_ylim()
        ax.plot([0,0], ylim, ls='-', color='Gray')
    ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
    ax.set_ylabel('Spec - Ref.')
    legend = ax.legend(loc=(1.05, .01), labelspacing=0.1)
    _plt.setp(legend.get_texts(),  fontsize='small')
    _plt.subplots_adjust(left=0.1, right=0.78, top=0.9, bottom=0.1)#, hspace=0.3, wspace=.3)   
    for fmt in formats:
        fig.savefig('fitsincr_lbc{1:.4f}{2}.{0}'.format(\
        fmt, lbc, addsuf), transparent=True)
    _plt.close()
    return


def diffPlotLineSeries(fullseds, obsers=[0], lbc=.6564606, formats=['png'], \
    rvel=None, rflx=None, hwidth=1000.):
    """Generate overplot of DIFFERENCE spec. line from a HDUST mod list.
    The model will be linearly interpolated
    with the reference spec. If none is given as reference, 
    then it assumes the first of the list.

    It is recommend to run first (rvel, rflx) =  spt.lineProf(rvel, rflx,
    lbc=lbc, hwidth=hwidth).

    Observers config. must be the same between models in
    `fullseds` list.
    """
    for obs in obsers:
        fig, ax = _plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = hdt.readfullsed2(file)
            obsdegs = (_np.arccos(sed2data[:,0,0])*180/_np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x,y) = lineProf(sed2data[obs,:,2], sed2data[obs,:,3], lbc=lbc, \
            hwidth=hwidth)
            if rvel == None or rflx == None:
                refspec = hdt.readfullsed2(fullseds[0])
                (vel,flx) = lineProf(refspec[obs,:,2], refspec[obs,:,3], lbc=lbc, \
                hwidth=hwidth)
            else:
                flx = _np.interp(x, rvel, rflx)
            if file == fullseds[0]:
                ax.plot(x, y-flx, label='{0:02.1f} deg.'.format(obsdegs[k]), \
                color=_phc.colors[_np.mod(i, len(_phc.colors))])
            else:
                ax.plot(x, y-flx, color=_phc.colors[_np.mod(i, len(_phc.colors))])
        ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
        ax.set_ylabel('Spec - Ref.')
        ax.legend()
        for fmt in formats:
            fig.savefig('modsdiff_lbc{2:.4f}_obs{0:02.1f}.{1}'.format(obsdegs[k], fmt,\
            lbc), transparent=True)
        _plt.close()
    return


def diffPlotLineFits(specs, lbc=.6564606, formats=['png'], \
    rvel=None, rflx=None, hwidth=1500., addsuf='', cmapn='jet'):
    """Generate overplot of DIFFERENCE spec. line from a FITS files list.
    The observations will be linearly interpolated
    with the reference spec. If none is given as reference, 
    then it assumes the first of the list.

    It is recommend to run first (rvel, rflx) =  spt.lineProf(rvel, rflx,
    lbc=lbc, hwidth=hwidth).

    If `cmap` is None or empty, the phc.colors vector is read.
    """
    fig, ax = _plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(_phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x,y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if rvel == None or rflx == None:
            #~ wl0, flux0, MJD, dateobs0, datereduc, fitsfile = loadfits(specs[0])
            #~ (rvel,flx) = lineProf(wl0, flux0, lbc=lbc, hwidth=hwidth)
            #~ flx = _np.interp(x, rvel, rflx)
            rvel = x
            rflx = y
            flx = y[:]
        else:
            flx = _np.interp(x, rvel, rflx)
        #~ if spec == specs[0]:
            #~ ax.plot(x, y-flx, label='{0}'.format(dateobs), \
            #~ color=_phc.colors[_np.mod(i, len(_phc.colors))])
        #~ else:
            #~ ax.plot(x, y-flx, color=_phc.colors[_np.mod(i, len(_phc.colors))])
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        if cmapn == '' or cmapn == None:
            color=_phc.colors[_np.mod(i, len(_phc.colors))]
        else:
            color=_phc.gradColor(range(len(specs)), cmapn='jet')[i]
        ax.plot(x, y-flx, label='{0}'.format(dateobs), color=color)
    ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
    ax.set_ylabel('Spec - Ref.')
    legend = ax.legend(loc=(1.05, .01), labelspacing=0.1)
    _plt.setp(legend.get_texts(),  fontsize='small')
    _plt.subplots_adjust(left=0.1, right=0.78, top=0.9, bottom=0.1)#, hspace=0.3, wspace=.3)   
    for fmt in formats:
        fig.savefig('fitsdiff_lbc{1:.4f}{2}.{0}'.format(\
        fmt, lbc, addsuf), transparent=True)
    _plt.close()
    return


def diffPlotLineObs(fullseds, obsers=[0], lbc=.6564606, formats=['png'], \
    rvel=None, rflx=None, hwidth=1000.):
    """Generate overplot of DIFFERENCE spec. line from a HDUST OBSERVERS list.
    The model will be linearly interpolated
    with the reference spec. If none is given as reference, 
    then it assumes the first observer of the list.

    It is recommend to run first (rvel, rflx) =  spt.lineProf(rvel, rflx,
    lbc=lbc, hwidth=hwidth).

    Observers config. must be the same between models in
    `fullseds` list.
    """
    for file in fullseds:
        fig, ax = _plt.subplots()
        sed2data = hdt.readfullsed2(file)
        obsdegs = (_np.arccos(sed2data[:,0,0])*180/_np.pi)[obsers]
        obsdegs = list(obsdegs)
        for obs in obsers:
            i = obsers.index(obs)
            (x,y) = lineProf(sed2data[obs,:,2], sed2data[obs,:,3], lbc=lbc, \
            hwidth=hwidth)
            if rvel == None or rflx == None:
                (vel,flx) = lineProf(sed2data[obsers[0],:,2], \
                sed2data[obsers[0],:,3], lbc=lbc, hwidth=hwidth)
            else:
                flx = _np.interp(x, rvel, rflx)
            ax.plot(x, y-flx, label='{0:02.1f} deg.'.format(obsdegs[i]), \
            color=_phc.colors[_np.mod(i, len(_phc.colors))])
        ax.set_title(u'lbc={0:.5f}$\mu$m, {1}'.format(lbc, \
        _phc.trimpathname(file)[1]))
        ax.set_ylabel('Spec - Ref.')
        ax.legend()
        for fmt in formats:
            fig.savefig('modsdiff_lbc{2:.4f}_{0}.{1}'.format(_phc.trimpathname(file)[1], \
            fmt, lbc), transparent=True)
        _plt.close()
    return


def kurlog(file=None, output=None):
    """ Generate a list of teff and logg present in a Kurucz file.

    If output is not specified, it is saved as `file`+.log """
    if file == None:
        file = '{0}/refs/fp00k0.pck'.format(_hdt.hdtpath())
    teffs = []
    loggs = []
    fp = open(file)
    for i, line in enumerate(fp):
        if line.find('TEFF') > -1:
            teffs+= [float(line.split()[1])]
            loggs+= [float(line.split()[3])]
    fp.close()
    return teffs, loggs


def kuruczflux(teff, logg, range=None):
    """ Return fluxes from a Kurucz model.

    Fluxes are in ergs/cm**2/s/hz/ster and wavelength in nm (range must be in
    nm).

    OUTPUT: wv, flux, info"""
    kurfile = '{0}/refs/fp00k0.pck'.format(_hdt.hdtpath())
    kurwvlines = (174-22)
    kurflxcol = 10
    #wave
    read = _phc.readrange(kurfile, 22, 22+kurwvlines)
    wave = _np.array([val for line in read for val in line.split()], dtype=float)
    #choose best
    bestT = _np.inf
    bestg = _np.inf
    fp = open(kurfile)
    for i, line in enumerate(fp):
        if line.find('TEFF') > -1:
            readT = float(line.split()[1])
            if _np.abs(readT-teff) <= _np.abs(bestT-teff):
                bestT = readT
                readg = float(line.split()[3])
                if _np.abs(readg-logg) <= _np.abs(bestg-logg):
                    bestg = readg
                    i0 = i+1
    fp.close()
    best = [bestT, bestg]
    #read best flux
    read = _phc.readrange(kurfile, i0, i0+kurwvlines)
    flux = _np.array([val for line in read for val in \
    (line[i:i+kurflxcol] for i in xrange(0, len(line)-1, kurflxcol))], \
    dtype=float)
    #cut range
    if range == None:
        return wave, flux, best
    else:
        idx = _np.where((wave > range[0]) & (wave < range[1]))
        return wave[idx], flux[idx], best    


def plotAll(files, obs=None, boxes=['s'], range=[[0,-1]], ncol=None, fmt=['png'],
    out=None):
    """ PlotAll routine for fullsed files.

    | `files` and `obs` are lists that define files and observers. 
    |
    | The number of boxes is defined by the `boxes` list.
    | Valid boxes are: 's' for SED, 'p' for polarimetry, 'l' for line.
    | 
    | `Range` is a list of lists and defines the interval for 's' and 'p'
    | (in microns). For automatic selection, leave [0,-1] for the respective box.
    | 'l' option must be [hwidht, lbd0] (mandatory).
    |
    | `ncol` fix the number of columns of the output file.

    The standard output name will be `date_time.png`, but it can be change
    with the variables `out` and `fmt`.

    TDB: normalize line, yrange and (x,y) log_scale !
    """
    for interv in range:
        if len(interv) != 2:
            print('# ERROR! Check range variable configuration!!!')
            return
    if len(files) != len(boxes) or len(files) != len(range):
        print('# ERROR! Check input configuration!!!')
        return
    #
    nbox = len(files)
    if ncol == None:
        nlin = int(round(_np.sqrt(nbox)))
        if nlin**2 <= nbox:
            ncol = nlin
        else:
            ncol = nlin+1
    else:
        nlin = nbox/ncol
        if nlin*ncol < nbox:
            ncol += 1
    #
    boxn = {'s':3, 'p':7, 'l':3}
    if obs == None:
        obs = [0]
    fig, axs = _plt.subplots(nlin, ncol)#, sharex=True, sharey=True)
    for file in files:
        sed2data = hdt.readfullsed2(file)
        for ob in obs:
            for I in _product(nlin, ncol):
                i,j = I
                type = boxn[boxes[i+j]]
                y = hdt.sed2data[ob,:,type]
                if boxes[i+j] == 'l':
                    print('linha')
                axs[i,j].plot(sed2data[ob,:,2], y)
    for I in _product(nlin, ncol):
        i,j = I
        if boxes[i+j] == 'l':
            axs[i,j].set_ylabel('Norm. flux')
        if boxes[i+j] == 'p':
            axs[i,j].set_ylabel('Q (%)')
        if boxes[i+j] == 's':
            axs[i,j].set_ylabel('Flux (arb. unit)')
        axs[i,j].set_xlabel(r'$\lambda$ ($\mu$m)')
        if list(range(i+j)) != [0,-1] and boxes[i+j] != 'l':
            axs[i,j].set_xlim(range(i+j))
    #
    if out == None:
        out = 'plotall_{0}{1}{2}_{3}{4}'.format(*time.localtime[:5])
    for format in fmt:
        _plt.savefig('{0}.{1}'.format(out, format))
        print('# File {0}.{1} saved!'.format(out, format))
    _plt.close()
    return


def splitKurucz(file, path=None):
    """
    Split atmospheric Kurucz file (e.g., 'ap00k0.dat') into individual models.

    INPUT: file, path (strings)

    OUTPUT: *files written
    """
    if path == None:
        path = _os.getcwd()
    allk = _np.loadtxt(file, dtype=str, delimiter='\n')
    
    for i in range(0,len(allk)-1):
        if 'EFF' in allk[i]:
            iref = i
            teff = int(allk[i].split()[1][:-1])
            logg = float(allk[i].split()[3][:-3])
        elif 'DECK6 72' in allk[i]:
            allk[i] = allk[i].replace('DECK6 72','DECK6 71')
        elif 'EFF' in allk[i+1]:
            _np.savetxt('ap00k0tef%05dg%.1f.dat' % (teff,logg), allk[iref:i+1], fmt='%s')
    
    _np.savetxt('ap00k0tef%05dg%.1f.dat' % (teff,logg), allk[iref:], fmt='%s')
    return


def din_spec(refspec, metadata, lbc=6562.86, hwidth=1500., res=50, interv=None,
    fmt=['png'], outname='din_spec', pxsize=8):
    """ Plot dynamical specs. from metadata table of the Spec class.

    `interv` controls the interval between specs.

    `res` is the resolution in km/s.

    TBD: average, binned specs.
    """
    #Define MJD and bins
    dates = _np.array(metadata[:,0], dtype=float)
    t0 = _np.min(dates)
    tf = _np.max(dates)
    if interv == None:
        interv = _np.linspace(t0,tf,21)
    else:
        interv = _np.arange(t0,tf+interv,interv)
    dt = interv[1]-interv[0]
    #Select specs 
    interv_names = _np.zeros(len(interv), dtype="|S512")
    wl0 = _np.arange(-hwidth,hwidth+res,res)
    fluxes = _np.ones(( len(wl0),len(interv) ))
    for i in range(len(interv)):
        date = _phc.find_nearest(dates,interv[i])
        #check if it is inside bin
        if date > interv[i]-dt/2 and date < interv[i]+dt/2:
            j = list(dates).index(date)
            interv_names = metadata[j,3]
            wl, flux, tmp, tmp, tmp, tmp = loadfits(refspec)
            wl, flux = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
            fluxes[:,i] = _np.interp(wl0, wl, flux)
    #Create image
    img = _np.empty((pxsize*len(interv),len(wl0)))
    for i in range(len(interv)):
        img[i*pxsize:(i+1)*pxsize] = _np.tile(fluxes[:,i],pxsize).reshape(pxsize,len(wl0))
    #Save image
    _plt.imshow(img)
    for ext in fmt:
        print(_os.getcwd(), 'hdt/{}.{}'.format(outname, ext))
        _plt.savefig('hdt/{}.{}'.format(outname, ext))
    _plt.close()
    return


def writeFits(flx, lbd, extrahead=None, savename=None, quiet=False, path=None,
    lbdc=None):
    """ Write a 1D spectra FITS.

    | INPUT: flux array, lbd array, extrahead flag+info, save name.
    | - lbd array: if len(lbd)==2: lbd = [CRVAL1, CDELT1]
    |              else: CDELT1 = (lbd[-1]-lbd[0])/(len(lbd)-1)
    |                    CRVAL1 = lbd[0]
    |   WARNING: lbd must be in ANGSTROMS (FITS default). It can also be 
    |   velocities. In this case, it must be in km/s and lbdc is given in ANGSTROM.
    | - extrahead: matrix (n,2). Example: [['OBJECT','Achernar'], ['COMMENT','value']]

    OUTPUT: write FITS file.
    """
    if path is None:
        path = _os.getcwd()
    if path[-1] != ['/']:
        path+= '/'
    if lbdc is not None:
        lbd = (lbd/_phc.c.cgs*1e5+1)*lbdc
    hdu = _pyfits.PrimaryHDU(flx)
    hdulist = _pyfits.HDUList([hdu])
    hdulist[0].header['CRVAL1'] = lbd[0]
    if len(lbd) == 2:
        hdulist[0].header['CDELT1'] = lbd[1]
    else:
        hdulist[0].header['CDELT1'] = (lbd[-1]-lbd[0])/(len(lbd)-1)
    if savename is None:
        savename = 'spec_{0}'.format(_phc.dtflag())
    if savename.find('.fit') == -1:
        savename+= '.fits'
    hdu.writeto(path+savename, clobber=True)
    if not quiet:
        print('# FITS file {0}{1} saved!'.format(path,savename))
    return


def plotSpecData(dtb, limits=[56470.,56720.], civcfg=[1,'m',2013,1,1],
    fmt=['png'], ident=None, lims=None, setylim=False, addsuf=''):
    """ Plot spec class database `vs` MJD e civil date

    Plot originally done to London, Canada, 2014.

    INPUT: civcfg = [step, 'd'/'m'/'y', starting year, month, day]

    `lims` sequence: 'EW', 'E/C', 'V/R', 'Pk. sep. (km/s)', 'E-F0', 'F0'

    `lims` = [[-2,4+2,2],[1.,1.4+.1,0.1],[.6,1.4+.2,.2],[0,400+100,100],
    [.30,.45+.05,.05],[0.6,1.20+.2,.2]]

    If `lims` is defined, `setylim` can be set to True.

    OUTPUT: Written image."""
    if isinstance(dtb, basestring):
        dtb = _np.loadtxt(dtb)
    if ident is not None:
        idref = _np.unique(ident)

    ylabels = ['EW', 'E/C', 'V/R', 'Pk. sep. (km/s)', 'E-F0', 'F0']
    fig, ax = _plt.subplots(6,1, sharex=True, figsize=(9.6,8))

    icolor = 'blue'
    for i in range(1,len(ylabels)+1):
        ax[i-1].plot(*_phc.bindata(dtb[:,0], dtb[:,i], 20))
        for j in range(len(dtb[:,0])):
            if ident is not None:
                idx = _np.where(ident[j] == idref)[0]
                icolor = _phc.colors[idx]
            ax[i-1].plot(dtb[j,0], dtb[j,i], 'o', color=icolor)
        ax[i-1].set_ylabel(ylabels[i-1])
        if lims is not None:
            if lims[i-1][-1] != 0:
                ax[i-1].set_yticks(_np.arange(*lims[i-1]))
            if setylim:
                ax[i-1].set_ylim([ lims[i-1][0],lims[i-1][1] ])

    if ident is not None:
        patch = []
        for id in idref:
            idx = _np.where(id == idref)[0]
            icolor = _phc.colors[idx]
            print icolor, id
            ax[0].plot([], [], 'o', color=icolor, label=id)
        ax[0].legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0., prop={'size':6})

    if limits is None:
        limits = ax[0].get_xlim()
    else:
        ax[0].set_xlim(limits)
    mjd0, mjd1 = limits
    ax[5].set_xlabel('MJD')
    ticks = _phc.gentkdates(mjd0, mjd1, civcfg[0], civcfg[1], dtstart=\
    _dt.datetime(civcfg[2],civcfg[3],civcfg[4]).date())
    mjdticks = [_jdcal.gcal2jd(date.year,date.month,date.day)[1] for date in \
    ticks]
    #ticks = [dt.datetime(*jdcal.jd2gcal(jdcal.MJD_0, date)[:3]).date() for \
    #date in ax[0].get_xticks()]
    #mjdticks = ax[0].get_xticks()
    for i in range(1,6+1):
        ax2 = ax[i-1].twiny()
        ax2.set_xlim(limits)
        ax2.set_xticks(mjdticks)
        ax2.set_xticklabels(['' for date in ticks])
        if i == 1:
            ax2.set_xlabel('Civil date')
            ax2.set_xticklabels([date.strftime("%d %b %y") for date in ticks])
            _plt.setp( ax2.xaxis.get_majorticklabels(), rotation=45 )
    _plt.subplots_adjust(left=0.13, right=0.8, top=0.88, bottom=0.06, hspace=.15)
    for f in fmt:
        _plt.savefig('SpecQ{1}.{0}'.format(f, addsuf), transparent=True)
    _plt.close()
    return


def cardelli(lbd, flux, ebv=0., Rv=3.1):
    """
    Milky Way Extinction law from Cardelli et al. 1989

    `lbd` must be in microns.

    OUTPUT: Corrected flux.
    """
    x=1./_np.array(lbd) #CCM x is 1/microns
    a, b=_np.ndarray(x.shape,x.dtype),_np.ndarray(x.shape,x.dtype)

    if any((x<0.3)|(10<x)):
        raise ValueError('Some wavelengths outside CCM 89 extinction curve range')

    irs = (0.3 <= x) & (x <= 1.1)
    opts = (1.1 <= x) & (x <= 3.3)
    nuv1s = (3.3 <= x) & (x <= 5.9)
    nuv2s = (5.9 <= x) & (x <= 8)
    fuvs = (8 <= x) & (x <= 10)

    #CCM Infrared
    a[irs]=.574*x[irs]**1.61
    b[irs]=-0.527*x[irs]**1.61

    #CCM NIR/optical
    a[opts]=_np.polyval((.32999,-.7753,.01979,.72085,-.02427,-.50447,.17699,1),x[opts]-1.82)
    b[opts]=_np.polyval((-2.09002,5.3026,-.62251,-5.38434,1.07233,2.28305,1.41338,0),x[opts]-1.82)

    #CCM NUV
    a[nuv1s]=1.752-.316*x[nuv1s]-0.104/((x[nuv1s]-4.67)**2+.341)
    b[nuv1s]=-3.09+1.825*x[nuv1s]+1.206/((x[nuv1s]-4.62)**2+.263)

    y=x[nuv2s]-5.9
    Fa=-.04473*y**2-.009779*y**3
    Fb=-.2130*y**2-.1207*y**3
    a[nuv2s]=1.752-.316*x[nuv2s]-0.104/((x[nuv2s]-4.67)**2+.341)+Fa
    b[nuv2s]=-3.09+1.825*x[nuv2s]+1.206/((x[nuv2s]-4.62)**2+.263)+Fb

    #CCM FUV
    a[fuvs]=_np.polyval((-.070,.137,-.628,-1.073),x[fuvs]-8)
    b[fuvs]=_np.polyval((.374,-.42,4.257,13.67),x[fuvs]-8)

    AlbAv = a+b/Rv
    return flux*10**(-AlbAv*Rv*ebv/2.5)


### MAIN ###
if __name__ == "__main__":
    pass
