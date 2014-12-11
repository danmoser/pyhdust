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
        (wl, flux) = linfit(self.wl, self.flux)
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
            (vels, flux) = linfit(vels[idx], flux[idx])
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
        if 'T' in dtobs:
            tobs, dtobs = dtobs.split('T')
            tobs = tobs.split(':')
            tobs = tobs[0]*3600+tobs[1]*60+tobs[2]
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
                            (vels, flux) = linfit(vels[idx], spdtb.flux[idx])
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

    Rydberg constant `R` manually adjusted to fit Halpha and Hbeta lines.
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
        return x, y
    else:
        yerr = yerr/np.average(new_y)
        return x, y, yerr

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
    (vels, flux) = linfit(vels, flux)
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

### MAIN ###
if __name__ == "__main__":
    pass
