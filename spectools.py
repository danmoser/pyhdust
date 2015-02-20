##!/usr/bin/env python
#-*- coding:utf-8 -*-
#Modified by D. Moser in 2014-10-17

"""
SPECTROSCOPY tools

includes None
"""

import matplotlib.pyplot as plt
import numpy as np
import pyhdust.phc as phc
import os
import glob
import pyfits
from scipy.interpolate import interp1d
import pyhdust.jdcal as jdcal
import pyhdust as hdt
from pyhdust import hdtpath

outfold = 'hdt/'

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
    >>> spdtb.data = np.vstack(( spdtb.data, np.loadtxt('hdt/datafile.txt') ))
    >>> spdtb.metadata = np.vstack(( spdtb.metadata, \\
    >>> np.loadtxt('hdt/metafile.txt') ))
    >>> spdtb.updatecount() #to update the counter

    Ou simplesmente (nome de arquivos default):

    >>> spdtb.loaddata()
    """
    def __init__(self, wl=None, flux=None, lbc=None, hwidth=None, EW=np.NaN,
        EC=np.NaN, VR=np.NaN,
        peaksep=np.NaN, depthcent=np.NaN, F0=np.NaN, dateobs='', MJD=0.,
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
        self.data = np.empty(0)
        self.metadata = np.empty(0)

    def reset(self):
        """Reset the class parameters
        """
        self.wl = None
        self.flux = None
        self.EW = np.NaN
        self.EC = np.NaN
        self.VR = np.NaN
        self.peaksep = np.NaN
        self.depthcent = np.NaN
        self.F0 = np.NaN
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
            self.data = np.array( self.lastinfo() )
            self.metadata = np.array( self.lastmeta() )
        else:
            self.data = np.vstack(( self.data, self.lastinfo() ))
            self.metadata = np.vstack(( self.metadata, self.lastmeta() ))
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

    def savedata(self, datafile=outfold+'/datafile.txt',\
        metafile=outfold+'/metafile.txt'):
        """Save current table
        """
        np.savetxt(datafile, self.data, fmt='%12.6f')
        np.savetxt(metafile, self.metadata, fmt='%s', delimiter=',')
        return

    def loaddata(self,  datafile=outfold+'/datafile.txt',\
        metafile=outfold+'/metafile.txt'):
        """Function to load a previous table
        Usage:

        >>> spdtb = Spec()
        >>> spdtb.loaddata()
        """
        self.data = np.loadtxt(datafile)
        self.metadata = np.loadtxt(metafile, fmt='%s', delimiter=',')
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
        if file[-4:] != 'fits':
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
            path, file = phc.trimpathname(self.file)
            outname = phc.rmext(file)
        #Normalization:
        flux = linfit(self.wl, self.flux)
        wl = self.wl
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(wl, flux)
        ax.set_ylabel('norm. flux')
        ax.set_xlabel('wavelength (arb. units)')
        ax.set_title(outname)
        plt.savefig('{}/{:.2f}_{}.png'.format(outfold,self.MJD,outname))
        if self.lbc > 0:
            vels = (self.wl-self.lbc)/self.lbc*phc.c.cgs*1e-5
            idx = np.where(np.abs(vels) <= self.hwidth)
            flux = linfit(vels[idx], flux[idx])
            vels = vels[idx]
            plt.clf()
            ax = fig.add_subplot(111)
            ax.plot(vels, flux)
            ax.set_ylabel('norm. flux')
            ax.set_xlabel('vel. (km/s)')
            ax.set_title('{:.2f} {} {:.2f}'.format(self.MJD,outname,self.lbc))
            plt.savefig('{}/{:.2f}_{}_{:.2f}.png'.format(outfold,self.MJD,\
            outname,self.lbc))
        plt.close()
        return

###
def extractfromsplot(file, splot):
    """Ce = center; Co = core
    #LcCe, LcCo, lcGW, lcEW, lvCe, lcCo, lvEW, lrCe, LrCo, lrEW
    """
    out = np.array(10*[np.NaN])
    readflag = False
    for line in splot:
        if line.find(']:') > 0 and readflag:
            readflag = False
        if line.find(file) > 0:
            readflag = True
        if readflag:
            info = line.split()
            #if re.match("^\d+?\.\d+?$", info[0]) is not None:
            try:
                float(info[0])
                info = np.array(info, dtype=float)
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
    dtobs = np.array(dtobs, dtype='int32')
    return dtobs, tobs

def shiftfits(fitsfile, newsh=''):
    """ Update FITS spec header for a given shift value. """
    imfits = pyfits.open(fitsfile, mode='update')
    if 'WLSHIFT' in imfits[0].header:
        print('# WLSHIFT = {0} for {1}'.format(imfits[0].header['WLSHIFT'], \
        phc.trimpathname(fitsfile)[1]))
    else:
        print('# No WLSHIFT available for {0}'.format( \
        phc.trimpathname(fitsfile)[1]))
    if newsh == '':
        newsh = raw_input('Type the new shift: ')
    if newsh != '':
        imfits[0].header['WLSHIFT'] = newsh
        imfits.close()
    return

def loadfits(fitsfile):
    """load FITS spec

    Out: wl, flux, MJD, dateobs, datereduc, fitsfile
    """
    imfits = pyfits.open(fitsfile)
    flux = imfits[0].data
    wl = np.arange(len(flux))*imfits[0].header['CDELT1']+\
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
        MJD = jdcal.gcal2jd(*dtobs)[1]+tobs
    else:
        MJD = 0.
        print('# ERROR! No DATE-OBS information is available! {}'.\
        format(fitsfile))
    if 'DATE-OBS' in imfits[0].header:
        dateobs = imfits[0].header['DATE-OBS']
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
    if os.path.exists('{}/{}'.format(path,star)) == False:
        os.system('mkdir {}/{}'.format(path,star))

    nights = [o for o in os.listdir(path) if os.path.isdir('{}/{}'.format(path,\
    o))]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    spdtb = Spec()
    spdtb.lbc = lbc
    spdtb.hwidth = 1000.
    for night in nights:
        targets = [o for o in os.listdir('%s/%s' % (path,night)) if
            os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = glob.glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:
                    for cal in scal:
                        spdtb.loadspec(cal)
                        spdtb.addspec()
                        if np.isnan(spdtb.EW) == False:
                            if plots:
                                spdtb.plotspec()
                            vels = (spdtb.wl-lbc)/lbc*phc.c.cgs*1e-5
                            idx = np.where(np.abs(vels) <= hwidth)
                            flux = linfit(vels[idx], spdtb.flux[idx])
                            vels = vels[idx]
                            leg = spdtb.MJD
                            ax.plot(vels, flux, label=leg, alpha=0.7, color=\
                            phc.colors[np.mod(spdtb.count, len(phc.colors))])
                else:   
                    print('# Data not reduced for %s at %s!' % (star,night))
    ax.set_xlim([-hwidth,hwidth])
    ax.set_ylim([-1,5])
    if showleg:
        legend = plt.legend(loc=(0.75, .05), labelspacing=0.1)
        plt.setp(legend.get_texts(),  fontsize='small')    
    plt.savefig('{}/{}_at_{}.png'.format(outfold,star,lbc))
    plt.close()
    spdtb.savedata(datafile='{}/{}.txt'.format(outfold,star),\
    metafile='{}/meta_{}.txt'.format(outfold,star))
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
    To HDUST effective R, multiple the input width by 2 here.

    # R = lbd/Dlbd = phc.c/Dv = phc.c*nbins/width 
    # nbins = R*width/phc.c
    """
    return round(phc.c.cgs*nbins/hwidth/1e5)

def calcres_nbins(R=12000, hwidth=1350):
    """
    (h)Width in km/s.
    *WARNING*: `width` in HDUST input is only half.
    To HDUST effective R, multiple the input width by 2 here.

    # R = lbd/Dlbd = phc.c/Dv = phc.c*nbins/width 
    # nbins = R*width/phc.c
    """
    return round(R*hwidth*1e5/phc.c.cgs)

def lineProf(x, flx, lbc, flxerr=np.empty(0), hwidth=1000., ssize=0.05):
    '''
    lineProf() - retorna um array (flx) normalizado e um array x em VELOCIDADES.
    `lbc` deve fornecido em mesma unidade de x para conversão lambda -> vel.
    Se vetor x jah esta em vel., usar funcao linfit().

    x eh importante, pois y pode ser nao igualmente amostrado.
    x e y devem estar em ordem crescente.

    ssize = % do tamanho de y; numero de pontos usados nas extremidades
    para a media do contínuo. 'ssize' de .5 à 0 (exclusive).
    '''    
    x = (x-lbc)/lbc*phc.c.cgs*1e-5 #km/s
    idx = np.where(np.abs(x) <= hwidth)
    flux = linfit(x[idx], flx[idx], yerr=flxerr, ssize=ssize)
    return x[idx], flux

def linfit(x, y, ssize=0.05, yerr=np.empty(0)):
    '''
    linfit() - retorna um array (y) normalizado, em posicoes de x

    x eh importante, pois y pode ser nao igualmente amostrado.
    x e y devem estar em ordem crescente.

    ssize = % do tamanho de y; numero de pontos usados nas extremidades
    para a media do contínuo. 'ssize' de .5 à 0 (exclusive).
    '''
    if ssize < 0 or ssize > .5:
        print('# Invalid ssize value...')
        raise SystemExit(1)
    ssize = int(ssize*len(y))
    if ssize == 0:
        ssize = 1
    medx0, medx1 = np.average(x[:ssize]),np.average(x[-ssize:])
    if ssize > 20:
        medy0, medy1 = np.median(y[:ssize]),np.median(y[-ssize:])
    else:
        medy0, medy1 = np.average(y[:ssize]),np.average(y[-ssize:])
    new_y = medy0 + (medy1 - medy0) * (x - medx0) / (medx1 - medx0)
    idx = np.where(new_y != 0)
    y[idx] = y[idx]/new_y[idx]
    if len(yerr) == 0.:
        return y
    else:
        yerr = yerr/np.average(new_y)
        return y, yerr

def EWcalc(vels, flux, vw=1000):
    """
    Supoe que o fluxo jah estah normalizado, e vetores ordenados.

    Devolve o valor EW, alem dos vetores cortados em vw.
    """
    idx = np.where(np.abs(vels) <= vw)
    outvels = vels[idx]
    normflux = flux[idx]
    ew = 0.
    if len(outvels) < 3:
        #normflux = np.ones(len(outvels))
        return ew
    for i in range(len(outvels)-2):
        dl = outvels[i+1]-outvels[i]
        ew += (1.-(normflux[i+1]+normflux[i])/2.)*dl
    return ew

def ECcalc(vels, flux):
    """
    Supoe que o fluxo jah estah normalizado, e vetores ordenados.

    Calcula o topo da emissao da linha, e retorna em que velocidade ela
    ocorre.
    """
    idx = np.where(np.max(flux) == flux)
    return flux[idx], vels[idx]

def VRcalc(vels, flux, vw=1000):
    """
    Supoe que o fluxo jah estah normalizado, e vetores ordenados.

    Calcula o ew para os dois lados (azul/vermelho) da linha, ajustando
    a velocidade de repouso. 
    """
    #calcula e aplica correcao de vel. repousp
    vc = 0.
    vels += vc
    #corta em vw, e faz o teste de tamanho
    if vw > 0:
        idx = np.where(np.abs(vels) <= vw)
        outvels = vels[idx]
        normflux = flux[idx]
    if len(vels) < 3:
        ew0, ew1 = 0.
        return ew0, ew1, vc
    #
    ivc = np.abs(outvels-0).argmin()
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
    contmax = np.max(np.append(flux[:ssize],flux[-ssize:]))
    fluxmax = np.max(flux)
    if fluxmax < 1.01*contmax:
        return 0, 0
    vels += vc
    ivc = np.abs(vels-0).argmin()
    i0 = np.abs(flux[:ivc]-np.max(flux[:ivc])).argmin()
    i1 = np.abs(flux[ivc+1:]-np.max(flux[ivc+1:])).argmin()+ivc+1
    return vels[i0], vels[i1]

def DCcalc(vels, flux, vmax=None, vc=0., ssize=0.05):
    """
    Calculo, na presenca de emissao, da profundidade do reverso central.
    
    """
    vels += vc
    ivc = np.abs(vels-0).argmin()
    #check if there is a peak
    ssize = int(ssize*len(vels))
    if ssize == 0:
        ssize = 1
    contmax = np.max(np.append(flux[:ssize],flux[-ssize:]))
    fluxmax = np.max(flux)
    if fluxmax < 1.01*contmax:
        return flux[ivc], flux[ivc]
    #if a vmax is not given...
    if isinstance(vmax, (int, long, float)) == False:
        vmax = np.abs(flux-np.max(flux)).argmin()
        vmax = vels[vmax]
    ivmax = np.abs(vels-vmax).argmin()
    return flux[ivmax], flux[ivc]

def analline(lbd, flux, lbdc, hwidth=1000, verb=True):
    """
    Return the analysis of a line.

    Both lbd and flux need to be ordered (a normalization IS FORCED).
    lbd,lbdc must have the same unit, and width in km/s is required.
    The line will be cutted so that the total DeltaLambda will be 2*width

    if lbdc <= 0, lbd array is assumed to be a velocity array (in km/s)!


    EXAMPLE: Using sed2data. lbc = 0.6565 (halpha), obs = 1 (width==1000)

        analline(lbd=sed2data[obs,:,2], flux=sed2data[obs,:,3], lbc=lbc)
    """
    if lbdc > 0:
        vels = (lbd-lbdc)/lbdc*phc.c.cgs*1e-5
    else:
        vels = lbd
    #check if the file have the desired info.
    if vels[0] > -hwidth or vels[-1] < hwidth:
        if verb:
            print('# ERROR: spec out of range (wavelength)!')
        return np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN

    idx = np.where(np.abs(vels) <= hwidth)
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
    #path = os.getcwd()
    #star = raw_input('Type the star name: ')
    #ref0 = 6540
    #ref1 = 6600

    ref0,ref1 = limits
    
    if os.path.exists('{}/{}'.format(path,star)) == False:
        os.system('mkdir {}/{}'.format(path,star))
    f0 = open('{0}/{1}/{1}.log'.format(path,star), 'w')
    
    nights = [o for o in os.listdir(path) if os.path.isdir('{}/{}'.format(path,\
        o))]
    
    legendl = ()
    i = 0
    for night in nights:
        targets = [o for o in os.listdir('%s/%s' % (path,night)) if
            os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = glob.glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = glob.glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    #legendl += (night,)
                    for cal in scal:              
                        imfits = pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = np.arange(len(spec))*imfits[0].header['CDELT1']+\
                            imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[-1] > 6560: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = np.where(abs(lbda-ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda,spec,a0,a1)
                            msg = '{}, {}, {}'.format((0.1*i), night, cal)
                            print(msg)
                            f0.writelines(msg+'\n')
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            plt.plot(lbda, spec, label=leg, alpha=0.7,
                                color=phc.colors[np.mod(i, len(phc.colors))])
                            i+=1
                else:
                    print('# Data not reduced for %s at %s!' % (star,night))
                    msg =  '{}, {}, {}'.format('NC', night, 'None')
                    f0.writelines(msg+'\n')

    if showleg:
        legend = plt.legend(loc=(0.75, .05), labelspacing=0.1)
        plt.setp(legend.get_texts(),  fontsize='small')
    plt.xlim([ref0,ref1])
    plt.ylim([-1, 5])
    #plt.xlabel('vel. (km/s)')
    plt.savefig('{0}/{1}/{1}_specs.png'.format(path,star))
    plt.close()
    f0.close()
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = np.arange(len(specs[i]))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
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
    #     a = np.where( abs(lbds[i]+1000) ==  min(abs(lbds[i]+1000)) )
    #     b = np.where( abs(lbds[i]-1000) ==  min(abs(lbds[i]-1000)) )
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
    # legend = plt.legend(legendl, loc=(0.75, .05), labelspacing=0.1)
    # plt.setp(legend.get_texts(),  fontsize='small')
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
    
    if os.path.exists('{}/{}'.format(path,star)) == False:
        os.system('mkdir {}/{}'.format(path,star))
    f0 = open('{0}/{1}/{1}.log'.format(path,star), 'w')
    
    nights = [o for o in os.listdir(path) if os.path.isdir('{}/{}'.format(path,\
        o))]    

    legendl = ()
    i = 0
    for night in nights:
        targets = [o for o in os.listdir('%s/%s' % (path,night)) if os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = glob.glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = glob.glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        print('# Specs with dopcor at %s!' % night)
                        srv = scal
                    #legendl += (night,)
                    for cal in scal:              
                        imfits = pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = np.arange(len(spec))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[0] > 5500: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = np.where(abs(lbda-ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda,spec,a0,a1)+(0.1*i)
                            print (0.1*i)
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            plt.plot([ref0,ref1], [1+0.1*i,1+0.1*i], 'k--', alpha=0.5)
                            plt.plot(lbda, spec, label=leg, color=phc.colors[i])
                            i+=1
                else:
                    print('# Data not reduced for %s at %s!' % (star,night))
    
    legend = plt.legend(loc=(0.75, .05), labelspacing=0.1)
    plt.setp(legend.get_texts(),  fontsize='small')
    plt.xlim([ref0,ref1])
    #plt.xlabel('vel. (km/s)')
    plt.savefig('{0}/{1}/{1}_specs_dif.png'.format(path,star))
    
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = np.arange(len(specs[i]))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
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
    #     a = np.where( abs(lbds[i]+1000) ==  min(abs(lbds[i]+1000)) )
    #     b = np.where( abs(lbds[i]-1000) ==  min(abs(lbds[i]-1000)) )
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
    # legend = plt.legend(legendl, loc=(0.75, .05), labelspacing=0.1)
    # plt.setp(legend.get_texts(),  fontsize='small')
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
    
    if os.path.exists('{}/{}'.format(path,star)) == False:
        os.system('mkdir {}/{}'.format(path,star))
    f0 = open('{0}/{1}/{1}.log'.format(path,star), 'w')
    
    nights = [o for o in os.listdir(path) if os.path.isdir('{}/{}'.format(path,\
        o))]

    legendl = ()
    i = 0
    
    for night in nights:
        targets = [o for o in os.listdir('%s/%s' % (path,night)) if os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = glob.glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = glob.glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        srv = scal
                    for cal in scal:              
                        imfits = pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = np.arange(len(spec))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[0] > 5500: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = np.where(abs(lbda-ref1) == min_dif)[0][0]
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
    specref = np.array(lines[1].split(), dtype=float)
    lbdaref = np.array(lines[0].split(), dtype=float)
    func = interp1d(lbdaref, specref)#, kind='cubic')
    lbdaref = np.linspace(ref0,ref1,5000)
    specref = func(lbdaref)
    
    i = 0
    for night in nights:
        targets = [o for o in os.listdir('%s/%s' % (path,night)) if os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = glob.glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = glob.glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    #legendl += (night,)
                    for cal in scal:              
                        imfits = pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = np.arange(len(spec))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[0] > 5500: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = np.where(abs(lbda-ref1) == min_dif)[0][0]
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
                                plt.plot(lbdaref, spec-specref, label=leg, alpha=0.8, color=phc.colors[i])
                            i+=1
                else:
                    print('# Data not reduced for %s at %s!' % (star,night))
    
    legend = plt.legend(loc=(0.75, .05), labelspacing=0.1)
    plt.setp(legend.get_texts(),  fontsize='small')
    plt.xlim([ref0,ref1])
    plt.title('Ref.= %s' % refleg)
    #plt.xlabel('vel. (km/s)')
    plt.savefig('{0}/{1}/{1}_specs_{2}.png'.format(path,star,refleg[:10]))
    
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = np.arange(len(specs[i]))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
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
    #     a = np.where( abs(lbds[i]+1000) ==  min(abs(lbds[i]+1000)) )
    #     b = np.where( abs(lbds[i]-1000) ==  min(abs(lbds[i]-1000)) )
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
    # legend = plt.legend(legendl, loc=(0.75, .05), labelspacing=0.1)
    # plt.setp(legend.get_texts(),  fontsize='small')
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
    
    if os.path.exists('{}/{}'.format(path,star)) == False:
        os.system('mkdir {}/{}'.format(path,star))
    f0 = open('{0}/{1}/{1}.log'.format(path,star), 'w')
    
    nights = [o for o in os.listdir(path) if os.path.isdir('{}/{}'.format(path,\
        o))]

    legendl = ()
    ax = plt.figure()
    i = 0
    for night in nights:
        targets = [o for o in os.listdir('%s/%s' % (path,night)) if os.path.isdir('%s/%s/%s' % (path,night,o))]
        for target in targets:
            if target.find(star)>-1:
                scal = glob.glob('%s/%s/%s/*.cal.fits' % (path,night,target))
                if len(scal) > 0:            
                    srv  = glob.glob('%s/%s/%s/*.rv.fits'  % (path,night,target))
                    if len(srv) != len(scal):
                        print('# Specs without dopcor at %s!' % night)
                        srv = scal
                    #legendl += (night,)
                    for cal in scal:              
                        imfits = pyfits.open(cal)
                        spec = imfits[0].data
                        lbda = np.arange(len(spec))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
                        #a = raw_input('type to continue: ')
                        if lbda[0] > 5500: #and flag == '1':
                            min_dif = min(abs(lbda-ref0))
                            a0 = np.where(abs(lbda-ref0) == min_dif)[0][0]
                            min_dif = min(abs(lbda-ref1))
                            a1 = np.where(abs(lbda-ref1) == min_dif)[0][0]
                            spec = normalize_range(lbda,spec,a0,a1)
                            print (0.1*i), night
                            prtcolor = phc.colors[i]
                            try:
                                leg = imfits[0].header['DATE-OBS']
                            except:
                                leg = imfits[0].header['FRAME']
                            check = False
                            if leg.find('2012-11-20T23:51:37.392') != -1:
                                leg = '2012-11-20'
                                prtcolor = phc.colors[0]
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
                                plt.plot(lbda, spec, label=leg, alpha=0.7, color=prtcolor)
                            i+=1
                else:
                    msg = '# Data not reduced for %s at %s!' % (star,night)
                    print(msg)
                    f0.writelines(msg)
                    
    
    font = { 'size' : 16, }
    
    legend = plt.legend(loc=(0.75, .05), labelspacing=0.1)
    #plt.setp(legend.get_texts(),  fontsize='small')
    plt.xlim([ref0,ref1])
    plt.ylim([.58,1.2])
    plt.xlabel(r'wavelength ($\AA$)', fontdict=font)
    plt.ylabel('Normalized flux', fontdict=font)
    #plt.xlabel('vel. (km/s)')
    plt.savefig('{0}/{1}/{1}_specs2.png'.format(path,star))
    plt.close()
    f0.close()
    # 
    # Ha = False # False do HeI 6678
    # 
    # for i in range(len(ifits)):
    #     imfits = pyfits.open(ifits[i])
    #     print imfits[0].header[3]
    #     specs[i][:len(imfits[0].data)] = imfits[0].data
    #     lbds[i] = np.arange(len(specs[i]))*imfits[0].header['CDELT1']+imfits[0].header['CRVAL1']
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
    #     a = np.where( abs(lbds[i]+1000) ==  min(abs(lbds[i]+1000)) )
    #     b = np.where( abs(lbds[i]-1000) ==  min(abs(lbds[i]-1000)) )
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
    # legend = plt.legend(legendl, loc=(0.75, .05), labelspacing=0.1)
    # plt.setp(legend.get_texts(),  fontsize='small')
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
        fig, ax = plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = hdt.readfullsed2(file)
            obsdegs = (np.arccos(sed2data[:,0,0])*180/np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x,y) = lineProf(sed2data[obs,:,2], sed2data[obs,:,3], lbc=lbc)
            if file == fullseds[0]:
                ax.plot(x, y, label='{0:02.1f} deg.'.format(obsdegs[k]), \
                color=phc.colors[np.mod(i, len(phc.colors))])
            else:
                ax.plot(x, y, color=phc.colors[np.mod(i, len(phc.colors))])
        ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
        ax.legend()
        for fmt in formats:
            fig.savefig('modsover_lbc{2:.4f}_obs{0:02.1f}.{1}'.format(obsdegs[k], fmt,\
            lbc), transparent=True)
        plt.close()
    return

def overPlotLineFits(specs, lbc=.6564606, formats=['png'], hwidth=1500., \
    ylim=None, yzero=False, addsuf=''):
    """Generate overplot spec. line from a FITS file list.
    """
    fig, ax = plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x,y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        ax.plot(x, y, label='{0}'.format(dateobs), \
        color=phc.colors[np.mod(i, len(phc.colors))])
    if ylim != None:
        ax.set_ylim(ylim)
    if yzero:
        ylim = ax.get_ylim()
        ax.plot([0,0], ylim, ls='-', color='Gray')
    ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
    ax.set_ylabel('Spec - Ref.')
    legend = ax.legend(loc=(1.05, .01), labelspacing=0.1)
    plt.setp(legend.get_texts(),  fontsize='small')
    plt.subplots_adjust(left=0.1, right=0.78, top=0.9, bottom=0.1)#, hspace=0.3, wspace=.3)   
    for fmt in formats:
        fig.savefig('fitsover_lbc{1:.4f}{2}.{0}'.format(\
        fmt, lbc, addsuf), transparent=True)
    plt.close()
    return

def incrPlotLineSeries(fullseds, obsers=[0], lbc=.6564606, formats=['png']):
    """Generate incremented spec. line from a HDUST mod list, separated by
    observers. The increment is 0.1 for each file in fullseds sequence.

    Observers config. must be the same between models in `fullseds` list.
    """
    for obs in obsers:
        fig, ax = plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = hdt.readfullsed2(file)
            obsdegs = (np.arccos(sed2data[:,0,0])*180/np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x,y) = lineProf(sed2data[obs,:,2], sed2data[obs,:,3], lbc=lbc)
            if file == fullseds[0]:
                ax.plot(x, y+0.1*i, label='{0:02.1f} deg.'.format(obsdegs[k]), \
                color=phc.colors[np.mod(i, len(phc.colors))])
            else:
                ax.plot(x, y+0.1*i, color=phc.colors[np.mod(i, len(phc.colors))])
        ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
        ax.legend()
        for fmt in formats:
            fig.savefig('modsincr_lbc{2:.4f}_obs{0:02.1f}.{1}'.format(obsdegs[k], fmt,\
            lbc), transparent=True)
        plt.close()
    return

def incrPlotLineFits(specs, lbc=.6564606, formats=['png'], hwidth=1500., \
    yzero=False, addsuf=''):
    """Generate incremented spec. line from FITS files list.
    The increment is 0.1 for each file in sequence.
    """
    fig, ax = plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x,y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        ax.plot(x, y+0.1*i, label='{0}'.format(dateobs), \
        color=phc.colors[np.mod(i, len(phc.colors))])
    if yzero:
        ylim = ax.get_ylim()
        ax.plot([0,0], ylim, ls='-', color='Gray')
    ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
    ax.set_ylabel('Spec - Ref.')
    legend = ax.legend(loc=(1.05, .01), labelspacing=0.1)
    plt.setp(legend.get_texts(),  fontsize='small')
    plt.subplots_adjust(left=0.1, right=0.78, top=0.9, bottom=0.1)#, hspace=0.3, wspace=.3)   
    for fmt in formats:
        fig.savefig('fitsincr_lbc{1:.4f}{2}.{0}'.format(\
        fmt, lbc, addsuf), transparent=True)
    plt.close()
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
        fig, ax = plt.subplots()
        k = obsers.index(obs)
        for file in fullseds:
            i = fullseds.index(file)
            sed2data = hdt.readfullsed2(file)
            obsdegs = (np.arccos(sed2data[:,0,0])*180/np.pi)[obsers]
            obsdegs = list(obsdegs)
            (x,y) = lineProf(sed2data[obs,:,2], sed2data[obs,:,3], lbc=lbc, \
            hwidth=hwidth)
            if rvel == None or rflx == None:
                refspec = hdt.readfullsed2(fullseds[0])
                (vel,flx) = lineProf(refspec[obs,:,2], refspec[obs,:,3], lbc=lbc, \
                hwidth=hwidth)
            else:
                flx = np.interp(x, rvel, rflx)
            if file == fullseds[0]:
                ax.plot(x, y-flx, label='{0:02.1f} deg.'.format(obsdegs[k]), \
                color=phc.colors[np.mod(i, len(phc.colors))])
            else:
                ax.plot(x, y-flx, color=phc.colors[np.mod(i, len(phc.colors))])
        ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
        ax.set_ylabel('Spec - Ref.')
        ax.legend()
        for fmt in formats:
            fig.savefig('modsdiff_lbc{2:.4f}_obs{0:02.1f}.{1}'.format(obsdegs[k], fmt,\
            lbc), transparent=True)
        plt.close()
    return

def diffPlotLineFits(specs, lbc=.6564606, formats=['png'], \
    rvel=None, rflx=None, hwidth=1500., addsuf=''):
    """Generate overplot of DIFFERENCE spec. line from a FITS files list.
    The observations will be linearly interpolated
    with the reference spec. If none is given as reference, 
    then it assumes the first of the list.

    It is recommend to run first (rvel, rflx) =  spt.lineProf(rvel, rflx,
    lbc=lbc, hwidth=hwidth).
    """
    fig, ax = plt.subplots()
    for spec in specs:
        i = specs.index(spec)
        print("# Reading {0}...".format(phc.trimpathname(spec)[1]))
        wl, flux, MJD, dateobs, datereduc, fitsfile = loadfits(spec)
        (x,y) = lineProf(wl, flux, lbc=lbc, hwidth=hwidth)
        if rvel == None or rflx == None:
            wl0, flux0, MJD, dateobs0, datereduc, fitsfile = loadfits(specs[0])
            (rvel,flx) = lineProf(wl0, flux0, lbc=lbc, hwidth=hwidth)
            flx = np.interp(x, rvel, rflx)
        else:
            flx = np.interp(x, rvel, rflx)
        #~ if spec == specs[0]:
            #~ ax.plot(x, y-flx, label='{0}'.format(dateobs), \
            #~ color=phc.colors[np.mod(i, len(phc.colors))])
        #~ else:
            #~ ax.plot(x, y-flx, color=phc.colors[np.mod(i, len(phc.colors))])
        if dateobs.find('-') > 0:
            dateobs = dateobs[:10]
        elif dateobs.find('/') > 0:
            dtobs = dateobs.split('/')[::-1]
            dateobs = "-".join(dtobs)
        ax.plot(x, y-flx, label='{0}'.format(dateobs), \
        color=phc.colors[np.mod(i, len(phc.colors))])
    ax.set_title(u'lbc = {0:.5f} $\mu$m'.format(lbc))
    ax.set_ylabel('Spec - Ref.')
    legend = ax.legend(loc=(1.05, .01), labelspacing=0.1)
    plt.setp(legend.get_texts(),  fontsize='small')
    plt.subplots_adjust(left=0.1, right=0.78, top=0.9, bottom=0.1)#, hspace=0.3, wspace=.3)   
    for fmt in formats:
        fig.savefig('fitsdiff_lbc{1:.4f}{2}.{0}'.format(\
        fmt, lbc, addsuf), transparent=True)
    plt.close()
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
        fig, ax = plt.subplots()
        sed2data = hdt.readfullsed2(file)
        obsdegs = (np.arccos(sed2data[:,0,0])*180/np.pi)[obsers]
        obsdegs = list(obsdegs)
        for obs in obsers:
            i = obsers.index(obs)
            (x,y) = lineProf(sed2data[obs,:,2], sed2data[obs,:,3], lbc=lbc, \
            hwidth=hwidth)
            if rvel == None or rflx == None:
                (vel,flx) = lineProf(sed2data[obsers[0],:,2], \
                sed2data[obsers[0],:,3], lbc=lbc, hwidth=hwidth)
            else:
                flx = np.interp(x, rvel, rflx)
            ax.plot(x, y-flx, label='{0:02.1f} deg.'.format(obsdegs[i]), \
            color=phc.colors[np.mod(i, len(phc.colors))])
        ax.set_title(u'lbc={0:.5f}$\mu$m, {1}'.format(lbc, \
        phc.trimpathname(file)[1]))
        ax.set_ylabel('Spec - Ref.')
        ax.legend()
        for fmt in formats:
            fig.savefig('modsdiff_lbc{2:.4f}_{0}.{1}'.format(phc.trimpathname(file)[1], \
            fmt, lbc), transparent=True)
        plt.close()
    return


def kurlog(file=None, output=None):
    """ Generate a list of teff and logg present in a Kurucz file.

    If output is not specified, it is saved as `file`+.log """
    if file == None:
        file = '{0}/refs/fp00k0.pck'.format(hdtpath())
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
    nm)."""
    kurfile = '{0}/refs/fp00k0.pck'.format(hdtpath())
    kurwvlines = (174-22)
    kurflxcol = 10
    #wave
    read = phc.readrange(kurfile, 22, 22+kurwvlines)
    wave = np.array([val for line in read for val in line.split()], dtype=float)
    #choose best
    bestT = np.inf
    bestg = np.inf
    fp = open(kurfile)
    for i, line in enumerate(fp):
        if line.find('TEFF') > -1:
            readT = float(line.split()[1])
            if np.abs(readT-teff) <= np.abs(bestT-teff):
                bestT = readT
                readg = float(line.split()[3])
                if np.abs(readg-logg) <= np.abs(bestg-logg):
                    bestg = readg
                    i0 = i+1
    fp.close()
    best = [bestT, bestg]
    #read best flux
    read = phc.readrange(kurfile, i0, i0+kurwvlines)
    flux = np.array([val for line in read for val in \
    (line[i:i+kurflxcol] for i in xrange(0, len(line)-1, kurflxcol))], \
    dtype=float)
    #cut range
    if range == None:
        return wave, flux, best
    else:
        idx = np.where((wave > range[0]) & (wave < range[1]))
        return wave[idx], flux[idx], best    


def plotSed2(files, obslist, tags=None, fmt='png', obsls=False, xlog=False,\
    ylog=False, xlim=None, ylim=None, norm=False, pol=False):
    """ Plot flux from Hdust models.

    `files` works as following: it must be (outs,,) matrix, where `n`
    controls the linestyle type  (i.e., how many files on the graph) and `m`
    the models which will make use of it.

    `obsls` [==True] sets observers as colors. Set it to False to observers as
    linestyles.

    `tags` is the output name for each list entry in `files`.
    """
    #
    if isinstance(files[0], list) == False:
        files = [files]
    if isinstance(fmt, list) == False:
        fmt = [fmt]
    if isinstance(obslist, list) == False:
        obslist = [obslist]
    obslist = list(obslist)
    counter = 0
    lflc = 0
    #~ color=phc.colors[np.mod(spdtb.count, len(phc.colors))
    for lfiles in files:
        flc = 0
        if tags is not None:
            outname = tags[lflc]
            lflc+= 1
        else:
            outname = phc.trimpathname(lfiles[0])[1].replace('.sed2','')
        if pol:
            outname='pol_'+outname
        fig, ax = plt.subplots()
        ax.set_title(outname)
        for file in lfiles:
            obsc = 0
            flc+= 1
            sed2data = readfullsed2(file)
            legname = phc.trimpathname(file)[1].replace('.sed2','')
            obsdeg = list(np.arccos(sed2data[:,0,0][obslist])*180/np.pi)
            for obs in obslist:
                obsc+= 1
                if obsls:
                    i = phc.colors[np.mod(flc-1,  len(phc.colors))]
                    j = phc.ls[np.mod(obsc-1, len(phc.ls))]
                else:
                    i = phc.colors[np.mod(obsc-1, len(phc.colors))]
                    j = phc.ls[np.mod(flc-1,  len(phc.ls))]
                #
                if pol:
                    k = 7
                    yk = 100.
                    ax.set_ylabel(r'$Q$(%)')
                else:
                    k = 3
                    yk = 1.
                    ax.set_ylabel(r'Flux')
                if norm and (xlim is not None):
                    idx = np.where( (sed2data[obs,:,2] >= xlim[0]) & \
                    (sed2data[obs,:,2] <= xlim[1]) )
                    x = sed2data[obs,:,2][idx]
                    y = linfit(x, sed2data[obs,:,k][idx]*yk)
                else:
                    x = sed2data[obs,:,2]
                    y = sed2data[obs,:,k]*yk
                ax.plot( x, y, \
                label='{0:.1f}o. {1}'.format(obsdeg[obslist.index(obs)],legname), \
                color=i, ls=j )
            #ax.set_title(file)
        ax.legend(prop={'size':6}, bbox_to_anchor=(1.7, 1.1))
        ax.set_xlabel(r'wavelength ($\mu$m)')
        if xlog:
            ax.set_xscale('log')
        if ylog:
            ax.set_yscale('log')
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        plt.subplots_adjust(left=0.05, right=0.6, top=0.9, bottom=0.1)#, hspace=0.3, wspace=.3)
        for ifmt in fmt:
            fig.savefig('{0}.{1}'.format(outname,ifmt), transparent=True)
        plt.close()
        counter += 1
    print('# {0} file generated!'.format(outname))
    return


### MAIN ###
if __name__ == "__main__":
    pass
