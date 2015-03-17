#-*- coding:utf-8 -*-

"""
PyHdust *interftools* module: interferometry tools

`colors` keep the *amdlib* standard.


A biblioteca python XDRLIB eh MUITO lenta... Usa muitas listas!!!

>>> import xdrlib

A biblioteca PYDAP estah em desenvolvimento... Eh complicada de usar

>>> from pydap.model import *
>>> from pydap.xdr import DapUnpacker
>>> base_int = BaseType(name='base_int')
>>> base_float = BaseType(name='base_float', type=Float32)

Todas as leituras binarias baseiam-se no struct.

:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""
import os as _os
import struct as _struct
import numpy as _np
from glob import glob as _glob
import pyhdust.phc as _phc
import pyhdust.oifits as _oifits
from pyhdust.spectools import linfit as _linfit

try:
    import matplotlib.pyplot as _plt
    import matplotlib.ticker as _mtick
    import pyfits as _pyfits
except:
    print('# Warning! matplotlib and/or pyfits module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


colors = ["red","green","blue","black"]

### READ MAP ###
def log_transform(im):
    '''returns log(image) scaled to the interval [0,1]'''
    try:
        (min, max) = (_np.min(im[_np.where(im > 0)]), _np.max(im))
        if (max > min) and (max > 0):
            im = (_np.log( im.clip(min, max) )-_np.log(min)) / (_np.log(max)-\
            _np.log(min))
            idx = _np.where(im == 0)
            im[idx] = _np.NaN
            return im
    except:
        pass
    return im

def dat2png(file):
    """
    Save the image in the path of the .dat file

    | First: Run the IDL routine "export_merged_file.pro"
    | files = ['/data/hdust/runs/hdust/aeri/mod07/\
    | Ha_mod07_n01.0e12_1.1yr_a1.0_Tsh09000_t00.80_Rd030.0_Be_aeri_2014_SEI_00_00.dat']
    | for file in files:
    |     dat2png(file)
    """
    f0 = open(file)
    dim = f0.readline().split()
    dim = int(float(dim[-1]))
    f0.close()
    #
    img = _np.loadtxt(file, comments='%')
    img = img.reshape((dim,dim))
    #
    _plt.figure()
    #_plt.imshow(img, cmap=_plt.get_cmap('gist_heat'))
    _plt.imshow(log_transform(img), cmap=_plt.get_cmap('gist_heat'))
    _plt.savefig(file.replace('.dat','.png'), transparent=True)
    _plt.savefig(file.replace('.dat','.eps'), transparent=True)
    #
    return

def imshowl(img):
    """
    """
    _plt.clf()
    _plt.imshow(log_transform(img), cmap=_plt.get_cmap('gist_heat'))
    return

def readmap(file, quiet=False):
    """
    Read MAP or MAPS files.

    output = data, obslist, lbdc, Ra, xmax

    data(nimgs,nobs,nlbd,ny,nx,dfact)

    .map, dfact = 6

    .maps, dfact = 1
    """
    if file[-4:] == '.map':
        dfact = 6
    elif file[-5:] == '.maps':
        dfact = 1
    else:
        print('# ERROR: This is not a HDUST valid image!')
        return
    f = open(file).read()
    #
    ixdr=0
    nobs,lnum,nx,ny = _struct.unpack('>4l', f[ixdr:ixdr+4*4])
    ixdr+=4*4
    Ra,Rstar,Lratio,xmax = _struct.unpack('>4f', f[ixdr:ixdr+4*4])
    ixdr+=4*4
    #nf no IDL estah como DOUBLE, mas com certea eh float ou int...
    #caso contrario nm nao faz sentido.
    nf = _struct.unpack('>f', f[ixdr:ixdr+1*4])[0]
    ixdr+=1*4
    nm = _struct.unpack('>l', f[ixdr:ixdr+1*4])[0]
    ixdr+=1*4
    
    npxs = nm
    upck = '>{}f'.format(npxs)
    xmax = _np.array( _struct.unpack(upck, f[ixdr:ixdr+npxs*4]) )
    ixdr+=npxs*4
    
    if file[-4:] == '.map':
        #tmp must be a "small" array. Otherwise, a MemoryError will be raised
        npxs = nx*ny*lnum*nobs*nm
        tmp = _np.empty((npxs*dfact))
        upck = '>{}f'.format(npxs)
        for i in range(dfact):
            #skip first npxs... See below
            tmp[i*npxs:(i+1)*npxs] = _np.array( _struct.unpack(upck,\
            f[ixdr:ixdr+npxs*4]) )
            ixdr+=npxs*4
        data = _np.zeros((nm,nobs,lnum,ny,nx,dfact+1))
        for i in range(dfact):
            data[:,:,:,:,:,i+1] = tmp[i::dfact].reshape((nm,nobs,lnum,ny,nx))
        data[:,:,:,:,:,0] = data[:,:,:,:,:,1]+data[:,:,:,:,:,2]+\
        data[:,:,:,:,:,3]
        #
    elif file[-5:] == '.maps':
        npxs = dfact*nx*ny*lnum*nobs*nm
        upck = '>{}f'.format(npxs)
        data = _np.array( _struct.unpack(upck, f[ixdr:ixdr+npxs*4]) ).\
        reshape((nm,nobs,lnum,ny,nx))
        ixdr+=npxs*4

    npxs = 2*nobs
    upck = '>{}f'.format(npxs)
    obslist = _np.array( _struct.unpack(upck, f[ixdr:ixdr+npxs*4]) )
    ixdr+=npxs*4

    npxs = lnum+1
    upck = '>{}f'.format(npxs)
    lbdarr = _np.array( _struct.unpack(upck, f[ixdr:ixdr+npxs*4]) )
    ixdr+=npxs*4

    #this will check if the XDR is finished.
    if ixdr == len(f):
        if quiet == False:
            print('# XDR {} completely read!'.format(file))
    else:
        print('# Warning: XDR {} not completely read!'.format(file))
        print('# length difference is {}'.format( (len(f)-ixdr)/4 ) )
    
    #lbdarr tem lnum+1, pois reflete o INTERVALO de cada imagem.
    # Para termos o lambda central de cada imagem, fazemos o seguinte:
    lbdc = _np.zeros(lnum)
    for i in range(lnum):
        lbdc[i] = (lbdarr[i]+lbdarr[i+1])/2.

    return data, obslist, lbdc, Ra, xmax

def img2fits(img, lbd, Ra, xmax, dist, outname='model'):
    """ Export an image (e.g., data[0,0,0,:,:]) to the fits format.

    `lbd` is the wavelength of the image and must be in meters. """
    hdu = _pyfits.PrimaryHDU(img[::-1,:])
    hdulist = _pyfits.HDUList([hdu])
    pixsize = 2*xmax[0]/len(img)
    rad_per_pixel = _np.double(pixsize*_phc.Rsun.cgs/(dist*_phc.pc.cgs))#*60.*60.*1000.*180./_np.pi)
    hdulist[0].header['CDELT1'] = rad_per_pixel
    hdulist[0].header['CDELT2'] = rad_per_pixel
    hdulist[0].header['CDELT3'] = 1.
    hdulist[0].header['CRVAL1'] = 0.
    hdulist[0].header['CRVAL2'] = 0.
    hdulist[0].header['CRVAL3'] = lbd
    hdulist[0].header['NAXIS'] = 2
    hdulist[0].header['NAXIS1'] = len(img)
    hdulist[0].header['NAXIS2'] = len(img[0])
    hdulist[0].header['CPPIX1'] = len(img)/2
    hdulist[0].header['CPPIX2'] = len(img[0])/2
    outname = '{0}.fits'.format(outname.replace(".fits",""))
    hdu.writeto(outname, clobber=True)
    print('# Saved {0} !'.format(outname))
    return

def data2fitscube(data, obs, lbdc, Ra, xmax, dist, zoom=0, outname='model'):
    """ Export a set of images (e.g., data[zoom,obs,:,:,:]) to the fits cube
    format.

    `lbdc` is the wavelength array and must be in meters. """
    hdu = _pyfits.PrimaryHDU(data[zoom,obs,:,:,:])
    hdulist = _pyfits.HDUList([hdu])
    pixsize = 2*xmax[0]/len(data[zoom,obs,0,:,:])
    rad_per_pixel = _np.double(pixsize*_phc.Rsun.cgs/(dist*_phc.pc.cgs))#*60.*60.*1000.*180./_np.pi)
    hdulist[0].header['CDELT1'] = rad_per_pixel
    hdulist[0].header['CDELT2'] = rad_per_pixel
    hdulist[0].header['CDELT3'] = (lbdc[-1]-lbdc[0])/len(lbdc)
    hdulist[0].header['CRVAL1'] = 0.
    hdulist[0].header['CRVAL2'] = 0.
    hdulist[0].header['CRVAL3'] = lbdc[0]
    hdulist[0].header['NAXIS'] = 2
    hdulist[0].header['NAXIS1'] = len(data[zoom,obs,0,:,:])
    hdulist[0].header['NAXIS2'] = len(data[zoom,obs,0,:,:][0])
    hdulist[0].header['CPPIX1'] = len(data[zoom,obs,0,:,:])/2
    hdulist[0].header['CPPIX2'] = len(data[zoom,obs,0,:,:][0])/2
    outname = '{0}.fits'.format(outname.replace(".fits",""))
    hdu.writeto(outname, clobber=True)
    print('# Saved {0} !'.format(outname))
    return

def genSquare(size=64, halfside=16, center=(0,0)):
    """
    Generate a square inside a squa_re.

    If size is not even, the unit square will not be centered.

    center is the relative position
    """
    x = _np.zeros(size)
    y = x[:,_np.newaxis]
    #ABSOLUTE Center:
    #if center is None:
    #    x0 = y0 = size // 2
    #else:
    #    x0 = center[0]
    #    y0 = center[1]
    #Relative Center:
    x0 = size // 2 + center[0]
    y0 = size // 2 + center[1]
    img = x*y
    img[x0-halfside:x0+halfside,y0-halfside:y0+halfside] = 1
    return img

def genGaussian(size=64, sig=64/8, center=(0,0)):
    """
    Generate a square gaussian kernel (non-normalized).

    size is the length of a side of the square (pixels)

    center is the relative position
    """
    x = _np.arange(0, size, 1, float)
    y = x[:,_np.newaxis]
    #ABSOLUTE Center:
    #if center is None:
    #    x0 = y0 = size // 2
    #else:
    #    x0 = center[0]
    #    y0 = center[1]
    #Relative Center:
    x0 = size // 2 + center[0]
    y0 = size // 2 + center[1]
    #
    return _np.exp(-0.5 * ((x-x0)**2 + (y-y0)**2) / sig**2)

def setspacecoords(nx, ny, rad_per_pixel, xc=0., yc=0.):
    """
    return xx and yy, 2D physical coordinates in ANGULAR dimensions
    (unit = defined by 'rad_per_pixel', i.e., radians).
    The physical scale (or length) on both axis must be the same.

    xc and yc are the the center position in PIXELS
    """
    x = _np.arange(0.,nx)-(nx-1)/2.+xc
    xx = _np.repeat(x, ny).reshape(-1, ny).T*rad_per_pixel
    y = _np.arange(0.,ny)-(ny-1)/2.+yc
    yy = _np.repeat(y, nx).reshape(-1, nx)*rad_per_pixel
    return  xx, yy

def fastnumvis(img, lbd, Bproj, PA, rad_per_pixel, PAdisk=90.):
    """
    For a given image (in phys.units = rad_per_pixel) and a interf. setup,
        it returns the visibility and phase.

    PA and PAdisk in degrees.
    """
    PA = PA-PAdisk+90.
    idx = _np.where(img > 0)
    
    u = Bproj*_np.double(_np.sin(PA/_phc.ra2deg)/lbd)
    v = Bproj*_np.double(_np.cos(PA/_phc.ra2deg)/lbd)
    #print PA,phc.ra2deg,lbd,Bproj,v

    ny = len(img)
    nx = len(img[0])
    xx, yy = setspacecoords(nx, ny, rad_per_pixel, xc=0., yc=0.)

    arg = -2*_np.pi*(xx[idx]*u + yy[idx]*v)
    TF_z_re = _np.sum(img[idx]*_np.cos(arg))
    TF_z_im = _np.sum(img[idx]*_np.sin(arg))
    #print TF_z_re,TF_z_im

    TF_z = complex(TF_z_re, TF_z_im)
    TF_z0 = _np.sum(img[idx])

    complexVis= TF_z/TF_z0
    
    VisAmp = _np.abs(complexVis)
    VisPhase = _np.arctan2(complexVis.imag, complexVis.real)*_phc.ra2deg
    return complexVis, VisAmp, VisPhase

def fastnumvis3(img, lbd, Bprojs, PAs, rad_per_pixel, PAdisk=90.):
    """
    Call the routine fastnumvis for each of the 3 baselines available.
    """
    u1 = Bprojs[0]*_np.cos(PAs[0]*_np.pi/180.)
    u2 = Bprojs[1]*_np.cos(PAs[1]*_np.pi/180.)
    v1 = Bprojs[0]*_np.sin(PAs[0]*_np.pi/180.)
    v2 = Bprojs[1]*_np.sin(PAs[1]*_np.pi/180.)
    B3 = _np.sqrt( (u1+u2)**2+(v1+v2)**2 )
    PA3 = _np.arctan2( u1+u2, v1+v2 )*180/_np.pi

    cV1, VA1, VP1 = fastnumvis(img, lbd, Bprojs[0], PAs[0], rad_per_pixel, PAdisk=PAdisk)
    cV2, VA2, VP2 = fastnumvis(img, lbd, Bprojs[1], PAs[1], rad_per_pixel, PAdisk=PAdisk)
    cV3, VA3, VP3 = fastnumvis(img, lbd, B3, PA3, rad_per_pixel, PAdisk=PAdisk)

    complexVis = cV1*cV2*cV3.conjugate()

    VisAmp = _np.abs(complexVis)
    VisPhase = _np.arctan2(complexVis.imag, complexVis.real)*_phc.ra2deg
    return complexVis, VisAmp, VisPhase


def plot_pionier(oidata, ffile='last_run', fmt=['png'], legend=True, model=None,
    obs=None, dist=None):
    """  Standard observational log for PIONIER

    obs is a list
    dist is a number
    """
    fig = _plt.figure()#figsize=(5.6,8))
    alp = .75
    ms = 3 #markersize
    #xloc = _plt.MaxNLocator(6)
    #ax0 = display info
    #~ ax0 = fig.add_subplot(211)
    #~ ax0.axis('off')
    #~ hdrinfo = oidata.hdrinfo.returninfo()
    #~ for i in range(4):
        #~ ax0.text(0., .8-.2*i, hdrinfo[i])
    #ax1 = uvplane/lambda
    ax1 = fig.add_subplot(321)
    colorid = 0
    names = []
    for vis2 in oidata.vis2:
        ulbd = vis2.ucoord/vis2.wavelength.eff_wave
        vlbd = vis2.vcoord/vis2.wavelength.eff_wave
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        if label not in names:
            names+= [label]
        color = _phc.colors[names.index(label)]
        #~ colorid = _np.mod(colorid+1, len(phc.colors))
        #[u,-u] = W > E
        #[-u,u] = E < W
        ax1.plot([ulbd,-ulbd],[vlbd,-vlbd], '.', color=color)#label=label,
    ax1.get_xaxis().set_ticklabels([])
    #~ ax1.xaxis.tick_top()
    #~ ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    #~ names = list(_np.unique(names))
    ax1.set_ylabel(u'B$_{proj}$/$\lambda$')
    ax1.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    #~ if legend: ax1.legend(prop={'size':8},numpoints=1,bbox_to_anchor=(-0.25, 1.0))
    #ax2 = uvplane
    ax2 = fig.add_subplot(322)
    colorid = 0
    leg = []
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = _phc.colors[names.index(label)]
        #~ color = phc.colors[colorid]
        #~ colorid = _np.mod(colorid+1, len(phc.colors))
        #[u,-u] = W > E
        #[-u,u] = E < W
        if label not in leg:
            ax2.plot([u,-u],[v,-v], '.', label=label, color=color)
            leg.append(label)
        else:
            ax2.plot([u,-u],[v,-v], '.', color=color)#label=label, 
        #~ names.append(vis2.target.target)
    #~ names = list(_np.unique(names))
    ax2.xaxis.tick_top()
    ax2.set_ylabel(u'B$_{proj}$ (m)')
    ax2.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    if legend: ax2.legend(prop={'size':8},numpoints=1,bbox_to_anchor=(1.05, 1.0))
    #ax3 = VIS2 vs. B
    #ax4 = VIS2 vs. PA
    #~ names = []
    colorid = 0
    plotid = 323
    ax3 = fig.add_subplot(plotid)
    plotid = 324
    ax4 = fig.add_subplot(plotid)
    for vis2 in oidata.vis2:
        u = vis2.ucoord/vis2.wavelength.eff_wave
        v = vis2.vcoord/vis2.wavelength.eff_wave
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = _phc.colors[names.index(label)]
        #~ color = phc.colors[colorid]
        #~ colorid = _np.mod(colorid + 1, len(phc.colors))
        line = ax3.errorbar(_np.sqrt(u**2 + v**2), \
        vis2.vis2data, yerr=vis2.vis2err, color=color, fmt='o', markersize=ms)#, label=label)
        #_np.arctan(self.ucoord / self.vcoord) * 180.0 / _np.pi % 180.0
        PAobs = _np.arctan2(u,v)*180.0/_np.pi
        #~ idx = _np.where(PAobs < 0)
        #~ PAobs[idx] = PAobs[idx]+180
        line = ax4.errorbar(PAobs, \
        vis2.vis2data, yerr=vis2.vis2err, color=color, fmt='o', markersize=ms)#, label=label)
    Blim = ax3.get_xlim()
    if model != None:
        #~ res = 20
        #~ V2 = _np.empty(res)
        for mod in model:
            if obs == None:
                obs = [0]
            data, obslist, lbdc, Ra, xmax = readmap(mod)
            pixsize = 2*xmax[0]/len(data[0,0,0,:,:])
            rad_per_pixel = _np.double(pixsize*_phc.Rsun.cgs/(dist*_phc.pc.cgs))#*60.*60.*1000.*180./_np.pi)
            for vis2 in oidata.vis2:
                if (vis2.station[0] and vis2.station[1]):
                    label = vis2.station[0].sta_name + vis2.station[1].sta_name
                else:
                    label = 'unnamed'
                color = _phc.colors[names.index(label)]
                u = vis2.ucoord
                v = vis2.vcoord
                B = _np.sqrt(u**2 + v**2)
                PA = _np.arctan2(u,v)*180.0/_np.pi
                avlbd = _np.average(vis2.wavelength.eff_wave)
                V2 = []
                lbds = []
                for i in range(len(vis2.wavelength.eff_wave)):
                    for ob in obs:
                        lbd = vis2.wavelength.eff_wave[i]
                        j = list(lbdc*1e-6).index(_phc.find_nearest(lbdc*1e-6,lbd))
                        #~ print lbdc[j]*1e-6-lbd
                        lbcalc = lbdc[j]*1e-6
                        lbcalc = lbd
                        tmp, V, phvar = fastnumvis(data[ob,0,j,:,:], lbcalc, B, PA, rad_per_pixel, PAdisk=(216.9+90.))
                        V2 += [V**2]
                        lbds += [lbcalc]
                #~ ax3.plot(B/vis2.wavelength.eff_wave, V2, color='purple', alpha=.3)
                ax3.plot(B/lbds, V2, color='purple', alpha=.3)
                ax4.plot(_np.tile(PA, len(lbds)), V2, color='purple', alpha=.3, marker='s', markersize=ms)
                #~ lbd = phc.find_nearest(lbdc*1e-6,avlbd)
                #~ j = list(lbdc*1e-6).index(lbd)
                #~ Bs = _np.linspace(Blim[0]*lbd, Blim[1]*lbd, res)
                #~ for i in range(res):
                    #~ tmp, V2[i], phvar = fastnumvis(data[0,0,j,:,:], lbd, Bs[i], PA, rad_per_pixel, PAdisk=36.9+90)
                #~ ax3.plot(Bs/lbd, V2**2, color=color)
                #~ print Bs/lbd
                #~ print V2
                #~ a = raw_input('asdads')
    ax3.set_xlim(Blim)
    #~ ax1.xaxis.get_major_formatter().set_useOffset(False)
    ax3.xaxis.set_major_formatter(_mtick.FormatStrFormatter('%.2e'))
    #~ ax3.get_xaxis().set_visible(False)
    #~ ax4.get_xaxis().set_visible(False)
    ax3.get_xaxis().set_ticklabels([])
    PAlim = [-180,180]
    ax3.set_ylim([0,1.1])
    ax3.set_ylabel(u'$V$ $^2$')
    ax3.grid(b=True, linestyle=':', alpha=alp)
    ax4.set_ylim([0,1.1])
    ax4.get_xaxis().set_ticklabels([])
    ax4.set_ylabel(u'$V$ $^2$')
    ax4.grid(b=True, linestyle=':', alpha=alp)
    #ax5 = T3PHI vs. B
    #ax6 = T3PHI vs. PA
    #~ names = []
    colorid = 0
    plotid = 325
    ax5 = fig.add_subplot(plotid)
    plotid = 326
    ax6 = fig.add_subplot(plotid)
    ax5.plot(Blim, [0,0], ls='--')
    ax6.plot(PAlim, [0,0], ls='--')
    for t3 in oidata.t3:
        u1 = t3.u1coord
        v1 = t3.v1coord
        u2 = t3.u2coord
        v2 = t3.v2coord
        if _np.sqrt(u1**2 + v1**2) > _np.sqrt(u2**2 + v2**2):
            u = u1; v = v1
        else:
            u = u2; v = v2
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u,v)*180.0/_np.pi
        B = _np.tile(B, len(t3.wavelength.eff_wave))
        PA = _np.tile(PA, len(t3.wavelength.eff_wave))
        #~ idx = _np.where(PA < 0)
        #~ PA[idx] = PA[idx]+180
        #~ if (t3.station[0] and t3.station[1]):
            #~ label = t3.station[0].sta_name + t3.station[1].sta_name
        #~ else:
            #~ label = 'unnamed'
        #~ color = names.index(label)
        #~ color = phc.colors[colorid]
        #~ colorid = _np.mod(colorid + 1, len(phc.colors))
        color = 'Black'
        #~ y = _np.repeat(t3.t3phi, len(t3.wavelength.eff_wave))
        #~ yerr = _np.repeat(t3.t3phierr, len(t3.wavelength.eff_wave))
        y = t3.t3phi
        yerr = t3.t3phierr
        line = ax5.errorbar(B/t3.wavelength.eff_wave, y, yerr=yerr, color=color, fmt='o', markersize=ms)#, label=label)
        line = ax6.errorbar(PA, y, yerr=yerr, color=color, fmt='o', markersize=ms)#, label=label)
    if model != None:
        for mod in model:
            if obs == None:
                obs = [0]
            data, obslist, lbdc, Ra, xmax = readmap(mod)
            pixsize = 2*xmax[0]/len(data[0,0,0,:,:])
            rad_per_pixel = _np.double(pixsize*_phc.Rsun.cgs/(dist*_phc.pc.cgs))#*60.*60.*1000.*180./_np.pi)
            for t3 in oidata.t3:
                u1 = t3.u1coord
                v1 = t3.v1coord
                u2 = t3.u2coord
                v2 = t3.v2coord
                B = _np.append(_np.sqrt(u1**2 + v1**2), _np.sqrt(u2**2 + v2**2))
                PA = _np.append(_np.arctan2(u1,v1)*180.0/_np.pi, _np.arctan2(u2,v2)*180.0/_np.pi)
                Bmax = []
                PAmax = []
                t3m = []
                lbds = []
                for i in range(len(t3.wavelength.eff_wave)):
                    for ob in obs:
                        lbd = t3.wavelength.eff_wave[i]
                        j = list(lbdc*1e-6).index(_phc.find_nearest(lbdc*1e-6,lbd))
                        lcalc = lbdc[j]*1e-6
                        lcalc = lbd
                        tmp, V, phvar = fastnumvis3(data[ob,0,j,:,:], lcalc, B, PA, rad_per_pixel, PAdisk=(216.9+90.))
                        t3m += [phvar]
                        Bmax += [_np.max(B)]
                        if _np.max(B) == B[0]:
                            PAmax += [PA[0]]
                        else:
                            PAmax += [PA[1]]
                        lbds += [lcalc]
                Bmax = _np.array(Bmax); lbds = _np.array(lbds); t3m = _np.array(t3m); PAmax = _np.array(PAmax)
                ax5.plot(Bmax/lbds, t3m, color='purple', alpha=.9)
                ax6.plot(PAmax, t3m, color='purple', alpha=.6, marker='s', markersize=ms)
    #~ ax5.get_xaxis().set_ticklabels([])
    ax5.set_xlim(Blim)
    ax6.set_xlim(PAlim)
    ymax = _np.max(_np.abs(ax5.get_ylim()))
    ax5.set_ylim([-1.05*ymax,1.05*ymax])
    ax6.set_ylim([-1.05*ymax,1.05*ymax])
    ax5.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax5.set_ylabel(u'$\phi_{123}$ (deg.)')
    ax5.grid(b=True, linestyle=':', alpha=alp)
    ax6.set_xlabel(u'$PA$ (deg.)')
    ax6.set_ylabel(u'$\phi_{123}$ (deg.)')
    ax6.grid(b=True, linestyle=':', alpha=alp)
    #SAVING
    dir, name = _phc.trimpathname(ffile)
    name = _phc.rmext(name)
    #_plt.savefig('hdt/{}_{}.png'.format(hdrinfo[0], hdrinfo[2]), transparent=True)
    #_plt.locator_params(axis = 'x', nbins = 7)
    _plt.subplots_adjust(left=0.12, right=0.95, top=0.96, bottom=0.09, hspace=.009, wspace=.32)
    if not _os.path.exists('hdt'):
        os.system('mkdir hdt')
    for suf in fmt:
        _plt.savefig('hdt/{0}.{1}'.format(name,suf), transparent=True)
    _plt.close()
    return

def plot_pionier_res(oidata, model, ffile='last_run', fmt=['png'], legend=True, 
    obs=None, dist=42.75, quiet=True):
    """  Obs-Model comparison for PIONIER

    model, obs are lists
    dist is a number
    """
    fig = _plt.figure()#figsize=(5.6,8))
    alp = .75
    ms = 3 #markersize
    #xloc = _plt.MaxNLocator(6)
    #ax0 = display info
    #~ ax0 = fig.add_subplot(211)
    #~ ax0.axis('off')
    #~ hdrinfo = oidata.hdrinfo.returninfo()
    #~ for i in range(4):
        #~ ax0.text(0., .8-.2*i, hdrinfo[i])
    #ax1 = uvplane/lambda
    ax1 = fig.add_subplot(321)
    colorid = 0
    names = []
    for vis2 in oidata.vis2:
        ulbd = vis2.ucoord/vis2.wavelength.eff_wave
        vlbd = vis2.vcoord/vis2.wavelength.eff_wave
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        if label not in names:
            names+= [label]
        color = _phc.colors[names.index(label)]
        #~ colorid = _np.mod(colorid+1, len(phc.colors))
        #[u,-u] = W > E
        #[-u,u] = E < W
        ax1.plot([ulbd,-ulbd],[vlbd,-vlbd], '.', color=color)#label=label,
    ax1.get_xaxis().set_ticklabels([])
    #~ ax1.xaxis.tick_top()
    #~ ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    #~ names = list(_np.unique(names))
    ax1.set_ylabel(u'B$_{proj}$/$\lambda$')
    ax1.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    #~ if legend: ax1.legend(prop={'size':8},numpoints=1,bbox_to_anchor=(-0.25, 1.0))
    #ax2 = uvplane
    ax2 = fig.add_subplot(322)
    colorid = 0
    leg = []
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = _phc.colors[names.index(label)]
        #~ color = phc.colors[colorid]
        #~ colorid = _np.mod(colorid+1, len(phc.colors))
        #[u,-u] = W > E
        #[-u,u] = E < W
        if label not in leg:
            ax2.plot([u,-u],[v,-v], '.', label=label, color=color)
            leg.append(label)
        else:
            ax2.plot([u,-u],[v,-v], '.', color=color)#label=label, 
        #~ names.append(vis2.target.target)
    #~ names = list(_np.unique(names))
    ax2.xaxis.tick_top()
    ax2.set_ylabel(u'B$_{proj}$ (m)')
    ax2.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    if legend: ax2.legend(prop={'size':8},numpoints=1,bbox_to_anchor=(1.05, 1.0))
    #ax3 = VIS2 vs. B
    #ax4 = VIS2 vs. PA
    #~ names = []
    colorid = 0
    plotid = 323
    ax3 = fig.add_subplot(plotid)
    plotid = 324
    ax4 = fig.add_subplot(plotid)
    ax4.plot([0, 1e8], [0,0], ls='--')
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = _phc.colors[names.index(label)]
        #~ color = phc.colors[colorid]
        #~ colorid = _np.mod(colorid + 1, len(phc.colors))
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u, v)*180/_np.pi
        line = ax3.errorbar(B/vis2.wavelength.eff_wave, \
        vis2.vis2data, yerr=vis2.vis2err, color=color, fmt='o', markersize=ms)#, label=label)
        #_np.arctan(self.ucoord / self.vcoord) * 180.0 / _np.pi % 180.0
        if obs == None:
            obs = [0]
        for mod in model:
            data, obslist, lbdc, Ra, xmax = readmap(mod, quiet=quiet)
            pixsize = 2*xmax[0]/len(data[0,0,0,:,:])
            rad_per_pixel = _np.double(pixsize*_phc.Rsun.cgs/(dist*_phc.pc.cgs))#*60.*60.*1000.*180./_np.pi)
            V2 = []
            lbds = []
            for i in range(len(vis2.wavelength.eff_wave)):
                for ob in obs:
                    lbd = vis2.wavelength.eff_wave[i]
                    j = list(lbdc*1e-6).index(_phc.find_nearest(lbdc*1e-6,lbd))
                    #~ lbcalc = lbdc[j]*1e-6
                    lbcalc = lbd
                    tmp, V, phvar = fastnumvis(data[ob,0,j,:,:], lbcalc, B, PA, rad_per_pixel, PAdisk=(216.9+90.))
                    V2 += [V**2]
                    lbds += [lbcalc]
            lbds = _np.array(lbds)
            #~ V2 = _np.array(V2)+(.0875e-8*B/lbds-.026)
            V2 = _np.array(V2)+(.0875e-8*B/lbds-.026)
            ax3.plot(B/lbds, V2, color='purple', alpha=.3)
            ax4.plot(B/lbds, (vis2.vis2data-V2)/vis2.vis2err, color='black', markersize=ms, marker='o', ls='')#marker='s',
    Blim = ax3.get_xlim()
    ax4.set_xlim(Blim)
    #~ ax1.xaxis.get_major_formatter().set_useOffset(False)
    ax3.xaxis.set_major_formatter(_mtick.FormatStrFormatter('%.2e'))
    #~ ax3.get_xaxis().set_visible(False)
    #~ ax4.get_xaxis().set_visible(False)
    ax3.get_xaxis().set_ticklabels([])
    PAlim = [-180,180]
    ax3.set_ylim([0,1.1])
    ax3.set_ylabel(u'$V$ $^2$')
    ax3.grid(b=True, linestyle=':', alpha=alp)
    ax4.set_ylim([-9,9])
    ax4.get_xaxis().set_ticklabels([])
    ax4.set_ylabel(u'$V$ $^2$(data-mod)/err')
    ax4.grid(b=True, linestyle=':', alpha=alp)
    #ax5 = T3PHI vs. B
    #ax6 = T3PHI vs. PA
    #~ names = []
    colorid = 0
    plotid = 325
    ax5 = fig.add_subplot(plotid)
    plotid = 326
    ax6 = fig.add_subplot(plotid)
    ax5.plot(Blim, [0,0], ls='--')
    ax6.plot(Blim, [0,0], ls='--')
    for t3 in oidata.t3:
        u1 = t3.u1coord
        v1 = t3.v1coord
        u2 = t3.u2coord
        v2 = t3.v2coord
        if _np.sqrt(u1**2 + v1**2) > _np.sqrt(u2**2 + v2**2):
            u = u1; v = v1
        else:
            u = u2; v = v2
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u,v)*180.0/_np.pi
        B = _np.tile(B, len(t3.wavelength.eff_wave))
        PA = _np.tile(PA, len(t3.wavelength.eff_wave))
        #~ idx = _np.where(PA < 0)
        #~ PA[idx] = PA[idx]+180
        #~ if (t3.station[0] and t3.station[1]):
            #~ label = t3.station[0].sta_name + t3.station[1].sta_name
        #~ else:
            #~ label = 'unnamed'
        #~ color = names.index(label)
        #~ color = phc.colors[colorid]
        #~ colorid = _np.mod(colorid + 1, len(phc.colors))
        color = 'Black'
        #~ y = _np.repeat(t3.t3phi, len(t3.wavelength.eff_wave))
        #~ yerr = _np.repeat(t3.t3phierr, len(t3.wavelength.eff_wave))
        y = t3.t3phi
        yerr = t3.t3phierr
        line = ax5.errorbar(B/t3.wavelength.eff_wave, y, yerr=yerr, color=color, fmt='o', markersize=ms)#, label=label)
        if obs == None:
            obs = [0]
        for mod in model:
            data, obslist, lbdc, Ra, xmax = readmap(mod, quiet=True)
            pixsize = 2*xmax[0]/len(data[0,0,0,:,:])
            rad_per_pixel = _np.double(pixsize*_phc.Rsun.cgs/(dist*_phc.pc.cgs))#*60.*60.*1000.*180./_np.pi)
            B = _np.append(_np.sqrt(u1**2 + v1**2), _np.sqrt(u2**2 + v2**2))
            PA = _np.append(_np.arctan2(u1,v1)*180.0/_np.pi, _np.arctan2(u2,v2)*180.0/_np.pi)
            Bmax = []
            PAmax = []
            t3m = []
            lbds = []
            for i in range(len(t3.wavelength.eff_wave)):
                for ob in obs:
                    lbd = t3.wavelength.eff_wave[i]
                    j = list(lbdc*1e-6).index(_phc.find_nearest(lbdc*1e-6,lbd))
                    lcalc = lbdc[j]*1e-6
                    lcalc = lbd
                    tmp, V, phvar = fastnumvis3(data[ob,0,j,:,:], lcalc, B, PA, rad_per_pixel, PAdisk=(216.9+90.))
                    t3m += [phvar]
                    Bmax += [_np.max(B)]
                    if _np.max(B) == B[0]:
                        PAmax += [PA[0]]
                    else:
                        PAmax += [PA[1]]
                    lbds += [lcalc]
            Bmax = _np.array(Bmax); lbds = _np.array(lbds); t3m = _np.array(t3m); PAmax = _np.array(PAmax)
            ax5.plot(Bmax/lbds, t3m, color='purple', alpha=.9)
            ax6.plot(Bmax/lbds, (y-t3m)/yerr, color='black', markersize=ms, marker='o', ls='')#alpha=.6, 
    #~ ax5.get_xaxis().set_ticklabels([])
    #~ labels = ax5.get_yticks().tolist()
    #~ labels[-1] = ''
    #~ ax5.set_yticklabels(labels)
    ax5.set_xlim(Blim)
    ax6.set_xlim(Blim)
    ymax = _np.max(_np.abs(ax5.get_ylim()))
    ax5.set_ylim([-1.05*ymax,1.05*ymax])
    ax6.set_ylim([-9,9])
    ax5.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax5.set_ylabel(u'$\phi_{123}$ (deg.)')
    ax5.grid(b=True, linestyle=':', alpha=alp)
    ax6.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax6.set_ylabel(u'$\phi_{123}$(data-mod)/err')
    ax6.grid(b=True, linestyle=':', alpha=alp)
    #SAVING
    dir, name = _phc.trimpathname(ffile)
    name = _phc.rmext(name)
    #_plt.savefig('hdt/{}_{}.png'.format(hdrinfo[0], hdrinfo[2]), transparent=True)
    #_plt.locator_params(axis = 'x', nbins = 7)
    _plt.subplots_adjust(left=0.12, right=0.95, top=0.96, bottom=0.09, hspace=.009, wspace=.32)
    for suf in fmt:
        _plt.savefig('{0}_res.{1}'.format(name,suf), transparent=True)
    _plt.close()
    return
    

def plot_oifits(oidata, ffile='last_run', fmt=['png'], xrange=None, legend=True):
    """ Standard observational log for AMBER

    If the file starts with "PRODUCT_", it searchs for the specs in the "AVG"
    folder.

    (One could write this info into the fits file. Since I've only tested the
    reading features of the `oifits` routine, I prefered do it this way).
    """
    #If it is a PRODUCT ffile, tries to load the AVG spec.
    specfile = ''
    if ffile.find('_PRO/PRODUCT') > 0:
        dateobs = ffile[ffile.find('_20')+1:ffile.find('_20')+20]
        dateobs2 = _phc.strrep(dateobs, -3, ':')
        dateobs2 = _phc.strrep(dateobs2, -6, ':')
        specfile = ffile[:ffile.find('_PRO/PRODUCT_')].replace('_SPEC','')+\
        '/*{0}*_OIDATA_AVG.fits*'.format(dateobs)
        specfile = _glob(specfile)
        if len(specfile) == 0:
            specfile = ffile[:ffile.find('_PRO/PRODUCT_')].replace('_SPEC','')+\
        '/*{0}*_OIDATA_AVG.fits*'.format(dateobs2)
            specfile = _glob(specfile)
        if len(specfile) != 1:
            specfile = ''
            print('# ERROR! This is a PRO oifits and the AVG file was '+\
            'not found!')
            print(ffile[:ffile.find('_PRO/PRODUCT_')]+('/*{0}*_OIDATA_'+\
        'AVG.fits*').format(dateobs))
        else:
            specfile = specfile[0]
            specoidata = _oifits.open(specfile, quiet=True)
            #print('# {0} file read!!!'.format(specfile))
            oidata.amberspec = specoidata.amberspec
            spec = oidata.amberspec[0]
            vis = oidata.vis[0]
            if len(spec.wavelength.eff_wave) == len(vis.wavelength.eff_wave):
                for i in range(len(oidata.vis)):
                    oidata.amberspec[i].wavelength = oidata.vis[i].wavelength
            else:
                print('# ERROR! spec and vis sizes are different for {0}'.\
                format(_phc.trimpathname(ffile)[1]))
    #
    fig = _plt.figure(figsize=(5.6,8))
    alp = .75
    #xloc = _plt.MaxNLocator(6)
    #ax0 = display info
    ax0 = fig.add_subplot(521)
    ax0.axis('off')
    hdrinfo = oidata.hdrinfo.returninfo()
    for i in range(4):
        ax0.text(0., .8-.2*i, hdrinfo[i])
    #ax2 = uvplane.
    ax2 = fig.add_subplot(522)
    colorid = 0
    names = []
    for vis in oidata.vis:
        if xrange == None:
            xmin = None
            xmax = None
            xmin=_np.amin(_np.append(1e6*vis.wavelength.eff_wave[\
            _np.where(vis.flag == False)], xmax))
            xmax=_np.amax(_np.append(1e6*vis.wavelength.eff_wave[\
            _np.where(vis.flag == False)], xmin))
            xrange = (xmin,xmax)
        u = vis.ucoord
        v = vis.vcoord
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        colorid = _np.mod(colorid+1, len(colors))
        #[u,-u] = W > E
        #[-u,u] = E < W
        ax2.plot([-u,u],[v,-v], '.', label=label, color=color)
        ax2.xaxis.tick_top()
        names.append(vis.target.target)
    ax2.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    if legend: ax2.legend(prop={'size':8},numpoints=1,bbox_to_anchor=(-0.25, 1.0))
    #ax3 = dif.phases, RIGHT COLUMN
    plotid = 524
    names = []
    colorid = 0
    yrange = [0,0]
    for vis in oidata.vis:
        diff = _np.max(_np.abs(vis.visphi))
        if diff > yrange[1]:
            yrange = [-diff, diff]
    for vis in oidata.vis:
        ax3=fig.add_subplot(plotid)
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        colorid = _np.mod(colorid + 1, len(colors))
        line = ax3.errorbar(1e6*vis.wavelength.eff_wave, vis.visphi,\
        vis.visphierr, label=label, color=color)
        ax3.set_ylim(yrange)
        ax3.set_xlim(xrange)
        ax3.set_ylabel(u'$\phi$ (deg.)')
        names.append(vis.target.target)
        names = list(_np.unique(names))
        #title = names.pop()
        #for name in names:
        #    title += ', %s'%(name)
        #ax1.set_title(title)
        #ax1.set_ylabel('Differential phase')
        plotid+=2
        #ax3.get_xaxis().set_visible(False)
        ax3.get_xaxis().set_ticklabels([])
        _plt.grid(b=True, linestyle=':', alpha=alp)
    #ax5 = closure phases
    #plotid = (5,2,10)
    names = []
    colorid = 3
    for t3 in oidata.t3:
        ax5=fig.add_subplot(5,2,10)
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        line = ax5.errorbar(1e6*vis.wavelength.eff_wave, t3.t3phi,\
        t3.t3phierr, label=label, color=color)
        ax5.set_xlim(xrange)
        #diff = _np.max(_np.abs(t3.t3phi))
        #yrange2 = [-diff, diff]
        ax5.set_ylim(yrange)
        names.append(vis.target.target)
        names = list(_np.unique(names))
    ax5.set_ylabel(u'Closure $\phi$ (deg.)')
    ax5.set_xlabel('Wavelength ($\mu$m)')
    _plt.setp( ax5.xaxis.get_majorticklabels(), rotation=-35 )
    _plt.grid(b=True, linestyle=':', alpha=alp)
    #ax4 = visibilities, LEFT COLUMN
    plotid = 523
    names = []
    colorid = 0
    yrange = [1,0]
    for vis in oidata.vis2:
        yrange = [ _np.min([yrange[0],_np.min(vis.vis2data)]),\
        _np.max([yrange[1],_np.max(vis.vis2data)]) ]
    if yrange[1] > 1.1:
        yrange[1] = 1.1
    if yrange[0] < 0 or yrange[0] >= 1:
        yrange[0] = 0
    for vis in oidata.vis2:
        ax4=fig.add_subplot(plotid)
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        colorid = _np.mod(colorid + 1, len(colors))
        line = ax4.errorbar(1e6*vis.wavelength.eff_wave, vis.vis2data,\
        vis.vis2err, label=label, color=color)
        names.append(vis.target.target)
        names = list(_np.unique(names))
        ax4.set_ylim(yrange)
        ax4.set_xlim(xrange)
        ax4.set_ylabel(u'V$^{2}$')
        #title = names.pop()
        #for name in names:
        #    title += ', %s'%(name)
        #ax1.set_title(title)
        #ax1.set_ylabel('Differential phase')
        plotid+=2
        ax4.get_xaxis().set_ticklabels([])
        _plt.grid(b=True, linestyle=':', alpha=alp)
    #ax1 = Line profile
    if True:
        ax1 = fig.add_subplot(529)
        colorid = 0
        names = []
        for spec in oidata.amberspec:
            label = 'unnamed'
            color = colors[colorid]
            #x,y,yerr = linfit(1e6*spec.wavelength.eff_wave, spec.spectrum, yerr=spec.spectrumerr)
            #ax1.errorbar(x, y, yerr, label=label, color=color)
            x = 1e6*spec.wavelength.eff_wave
            y = _linfit(1e6*spec.wavelength.eff_wave, spec.spectrum)
            ax1.plot(x,y, label=label, color=color)
            colorid = _np.mod(colorid+1, len(colors))
            #[-u,u] = W > E
            #[u,-u] = E < W
        ax1.set_xlabel('Wavelength ($\mu$m)')
        ax1.set_xlim(xrange)
        ax1.set_ylabel('Norm. flux')
        _plt.setp( ax1.xaxis.get_majorticklabels(), rotation=-35 )
        _plt.grid(b=True, linestyle=':', alpha=alp)
    _phc.outfld()
    dir, name = _phc.trimpathname(ffile)
    name = _phc.rmext(name)
    #_plt.savefig('hdt/{}_{}.png'.format(hdrinfo[0], hdrinfo[2]), transparent=True)
    #_plt.locator_params(axis = 'x', nbins = 7)
    _plt.subplots_adjust(left=0.12, right=0.95, top=0.96, bottom=0.09, hspace=.009, wspace=.32)
    for suf in fmt:
        _plt.savefig('hdt/{0}.{1}'.format(name,suf), transparent=True)
    _plt.close()
    return

def genfinaloifits(oidata, ffile, xrange=None, legend=True):
    """ Standard observational log
    WITH the CORRECTED SPECTRUM
    """    
    print('# WAIT... Work in progress')
    return

def readesoquery(file):
    """ Read ESO query CSV ('utf-8-sig').

    There is a bug: the delimiter "," is used in the fields!!! This researches
    the line to replace the "," inside '"' symbols.
    """
    import codecs
    f0 = open(file)
    #DO NOT WORK
    #lines = f0.read().decode('utf-8-sig').encode('utf-8')
    lines = f0.readlines()
    f0.close()
    if lines[0].startswith(codecs.BOM_UTF8):
        lines[0] = lines[0].replace(codecs.BOM_UTF8, '', 1)
    outlines = []
    for i in range(len(lines)):
        if lines[i] != '\n' and lines[i] != '' and lines[i][0] != '#':
            k = lines[i].count('"')
            if k%2 == 1:
                print('# ERROR! Strange number os strings in line {} of {}'.format(\
                i, file))
                print lines[i]
                raise SystemExit(1)
            itmp = 0
            for l in range(k/2):
                i0 = lines[i][itmp:].find('"')+itmp
                i1 = i0+1
                i1 = lines[i][i1:].find('"')+i1
                itmp = i1+1
                lines[i] = lines[i][:i0]+lines[i][i0:i1+1].replace(',','_')+\
                    lines[i][i1+1:]
                #print l, lines[i], itmp, i0, i1
            outlines += [lines[i].replace('"','')]
    #print outlines[:2]
    f0 = open(file.replace('.csv','.txt'), 'w')
    f0.writelines(outlines)
    f0.close()
    return 

def checkESOdownload(path=None):
    """ check ESO download """
    if path == None:
        path = _os.getcwd()
    sh = _glob(path+'/*.sh')
    f0 = open(sh[0])
    lines = f0.readlines()
    f0.close()
    count = 0
    count = len(_glob(path+'/*.Z'))
    count+= len(_glob(path+'/*.txt'))
    count+= len(_glob(path+'/notused/*'))
    print(path,len(lines),count)
    return

def printinfo(file):
    """ Print OIFITS observational log, as appendix of Faes, D. M. (2015)

    DATE-OBS    MJD PA  B   PA  B   PA  B"""
    oidata = _oifits.open(file, quiet='True')
    info = list(oidata.hdrinfo.returninfo())
    info2 = []
    for vis in oidata.vis:
        info2+= ['{0:.1f}'.format(_np.sqrt(vis.ucoord**2 + vis.vcoord**2))]
        info2+= ['{0:.1f}'.format(_np.arctan(vis.ucoord / vis.vcoord) * 180.0 / _np.pi % 180.0)]
    return [info[2][:10], '{0:.7f}'.format(info[1])] + info2 


### MAIN ###
if __name__ == "__main__":
    pass
