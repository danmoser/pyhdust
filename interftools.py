#!/usr/bin/env python
#-*- coding:utf-8 -*-
#Modified by D. Moser in 2014-10-17

"""
INTERFEROMETRY tools

includes *readmap

About *readmap:
A biblioteca XDRLIB eh MUITO lenta... Usa muitas listas!!!

>>> import xdrlib

A biblioteca PYDAP estah em desenvolvimento... Eh complicada de usar

>>> from pydap.model import *
>>> from pydap.xdr import DapUnpacker
>>> base_int = BaseType(name='base_int')
>>> base_float = BaseType(name='base_float', type=Float32)

Todas as leituras binarias baseiam-se no struct
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import struct
import pyhdust.phc as phc
import pyhdust.oifits as oifits
import os
from glob import glob
from pyhdust.spectools import linfit
#from pyhdust import hdtpath

colors = ["red","green","blue","black"]

### READ MAP ###
def log_transform(im):
    '''returns log(image) scaled to the interval [0,1]'''
    try:
        (min, max) = (np.min(im[np.where(im > 0)]), np.max(im))
        if (max > min) and (max > 0):
            im = (np.log( im.clip(min, max) )-np.log(min)) / (np.log(max)-\
            np.log(min))
            idx = np.where(im == 0)
            im[idx] = np.NaN
            return im
    except:
        pass
    return im

def dat2png(file):
    """
    Save the image in the path of the .dat file

    First: Run the IDL routine "export_merged_file.pro"
    files = ['/data/hdust/runs/hdust/aeri/mod07/\
    Ha_mod07_n01.0e12_1.1yr_a1.0_Tsh09000_t00.80_Rd030.0_Be_aeri_2014_SEI_00_00.dat']
    for file in files:
        dat2png(file)
    """
    f0 = open(file)
    dim = f0.readline().split()
    dim = int(float(dim[-1]))
    f0.close()
    #
    img = np.loadtxt(file, comments='%')
    img = img.reshape((dim,dim))
    #
    plt.figure()
    #plt.imshow(img, cmap=plt.get_cmap('gist_heat'))
    plt.imshow(log_transform(img), cmap=plt.get_cmap('gist_heat'))
    plt.savefig(file.replace('.dat','.png'), transparent=True)
    plt.savefig(file.replace('.dat','.eps'), transparent=True)
    #
    return

def imshowl(img):
    """
    """
    plt.clf()
    plt.imshow(log_transform(img), cmap=plt.get_cmap('gist_heat'))
    return

def readmap(file):
    """
    MAP ou MAPS > Returns PNG !!!
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
    nobs,lnum,nx,ny = struct.unpack('>4l', f[ixdr:ixdr+4*4])
    ixdr+=4*4
    Ra,Rstar,Lratio,xmax = struct.unpack('>4f', f[ixdr:ixdr+4*4])
    ixdr+=4*4
    #nf no IDL estah como DOUBLE, mas com certea eh float ou int...
    #caso contrario nm nao faz sentido.
    nf = struct.unpack('>f', f[ixdr:ixdr+1*4])[0]
    ixdr+=1*4
    nm = struct.unpack('>l', f[ixdr:ixdr+1*4])[0]
    ixdr+=1*4
    
    npxs = nm
    upck = '>{}f'.format(npxs)
    xmax = np.array( struct.unpack(upck, f[ixdr:ixdr+npxs*4]) )
    ixdr+=npxs*4
    
    if file[-4:] == '.map':
        #tmp must be a "small" array. Otherwise, a MemoryError will be raised
        npxs = nx*ny*lnum*nobs*nm
        tmp = np.empty((npxs*dfact))
        upck = '>{}f'.format(npxs)
        for i in range(dfact):
            #skip first npxs... See below
            tmp[i*npxs:(i+1)*npxs] = np.array( struct.unpack(upck,\
            f[ixdr:ixdr+npxs*4]) )
            ixdr+=npxs*4
        data = np.zeros((nm,nobs,lnum,ny,nx,dfact+1))
        for i in range(dfact):
            data[:,:,:,:,:,i+1] = tmp[i::dfact].reshape((nm,nobs,lnum,ny,nx))
        data[:,:,:,:,:,0] = data[:,:,:,:,:,1]+data[:,:,:,:,:,2]+\
        data[:,:,:,:,:,3]
        #
    elif file[-5:] == '.maps':
        npxs = dfact*nx*ny*lnum*nobs*nm
        upck = '>{}f'.format(npxs)
        data = np.array( struct.unpack(upck, f[ixdr:ixdr+npxs*4]) ).\
        reshape((nm,nobs,lnum,ny,nx))
        ixdr+=npxs*4

    npxs = 2*nobs
    upck = '>{}f'.format(npxs)
    obslist = np.array( struct.unpack(upck, f[ixdr:ixdr+npxs*4]) )
    ixdr+=npxs*4

    npxs = lnum+1
    upck = '>{}f'.format(npxs)
    lbdarr = np.array( struct.unpack(upck, f[ixdr:ixdr+npxs*4]) )
    ixdr+=npxs*4

    #this will check if the XDR is finished.
    if ixdr == len(f):
        print('# XDR {} completely read!'.format(file))
    else:
        print('# Warning: XDR {} not completely read!'.format(file))
        print('# length difference is {}'.format( (len(f)-ixdr)/4 ) )
    
    #lbdarr tem lnum+1, pois reflete o INTERVALO de cada imagem.
    # Para termos o lambda central de cada imagem, fazemos o seguinte:
    lbdc = np.zeros(lnum)
    for i in range(lnum):
        lbdc[i] = (lbdarr[i]+lbdarr[i+1])/2.

    return data, obslist, lbdc, Ra, xmax

def genSquare(size=64, halfside=16, center=(0,0)):
    """
    Generate a square inside a square.

    If size is not even, the unit square will not be centered
    center is the relative position
    """
    x = np.zeros(size)
    y = x[:,np.newaxis]
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
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
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
    return np.exp(-0.5 * ((x-x0)**2 + (y-y0)**2) / sig**2)

def setspacecoords(nx, ny, rad_per_pixel, xc=0., yc=0.):
    """
    return xx and yy, 2D physical coordinates in ANGULAR dimensions
    (unit = defined by 'rad_per_pixel', i.e., radians).
    The physical scale (or length) on both axis must be the same.

    xc and yc are the the center position in PIXELS
    """
    x = np.arange(0.,nx)-(nx-1)/2.+xc
    xx = np.repeat(x, ny).reshape(-1, ny).T*rad_per_pixel
    y = np.arange(0.,ny)-(ny-1)/2.+yc
    yy = np.repeat(y, nx).reshape(-1, nx)*rad_per_pixel
    return  xx, yy

def fastnumvis(img, lbd, Bproj, PA, rad_per_pixel, PAdisk=90.):
    """
    For a given image (in phys.units = rad_per_pixel) and a interf. setup,
        it returns the visibility and phase.

    PA and PAdisk in degrees.
    """
    PA = PA+PAdisk
    idx = np.where(img > 0)
    
    u = Bproj*np.double(np.sin(PA/phc.ra2deg)/lbd)
    v = Bproj*np.double(np.cos(PA/phc.ra2deg)/lbd)
    #print PA,phc.ra2deg,lbd,Bproj,v

    ny = len(img)
    nx = len(img[0])
    xx, yy = setspacecoords(nx, ny, rad_per_pixel, xc=0., yc=0.)

    arg = -2*np.pi*(xx[idx]*u + yy[idx]*v)
    TF_z_re = np.sum(img[idx]*np.cos(arg))
    TF_z_im = np.sum(img[idx]*np.sin(arg))
    #print TF_z_re,TF_z_im

    TF_z = complex(TF_z_re, TF_z_im)
    TF_z0 = np.sum(img[idx])

    complexVis= TF_z/TF_z0
    
    VisAmp = np.abs(complexVis)
    VisPhase = np.arctan2(complexVis.imag, complexVis.real)*phc.ra2deg
    return complexVis, VisAmp, VisPhase

def fastnumvis3():
    """
    Call the routine fastnumvis for each of the 3 baselines available.
    """
    return

def fastnumvis4():
    """
    Call the routine fastnumvis for each of the 3 baselines available.
    """
    return

def plot_pionier(oidata, ffile='last_run', fmt=['png'], legend=True):
    """  Standard observational log for PIONIER """
    fig = plt.figure()#figsize=(5.6,8))
    alp = .75
    ms = 3 #markersize
    #xloc = plt.MaxNLocator(6)
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
        color = phc.colors[names.index(label)]
        #~ colorid = np.mod(colorid+1, len(phc.colors))
        #[u,-u] = W > E
        #[-u,u] = E < W
        ax1.plot([ulbd,-ulbd],[vlbd,-vlbd], '.', color=color)#label=label,
    ax1.get_xaxis().set_ticklabels([])
    #~ ax1.xaxis.tick_top()
    #~ ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    #~ names = list(np.unique(names))
    ax1.set_ylabel(u'B$_{proj}$/$\lambda$')
    ax1.axis('equal')
    plt.grid(b=True, linestyle=':', alpha=alp)
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
        color = phc.colors[names.index(label)]
        #~ color = phc.colors[colorid]
        #~ colorid = np.mod(colorid+1, len(phc.colors))
        #[u,-u] = W > E
        #[-u,u] = E < W
        if label not in leg:
            ax2.plot([u,-u],[v,-v], '.', label=label, color=color)
            leg.append(label)
        else:
            ax2.plot([u,-u],[v,-v], '.', color=color)#label=label, 
        #~ names.append(vis2.target.target)
    #~ names = list(np.unique(names))
    ax2.xaxis.tick_top()
    ax2.set_ylabel(u'B$_{proj}$ (m)')
    ax2.axis('equal')
    plt.grid(b=True, linestyle=':', alpha=alp)
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
        color = phc.colors[names.index(label)]
        #~ color = phc.colors[colorid]
        #~ colorid = np.mod(colorid + 1, len(phc.colors))
        line = ax3.errorbar(np.sqrt(u**2 + v**2), \
        vis2.vis2data, yerr=vis2.vis2err, color=color, fmt='o', markersize=ms)#, label=label)
        #np.arctan(self.ucoord / self.vcoord) * 180.0 / np.pi % 180.0
        PAobs = np.arctan2(u,v)*180.0/np.pi
        #~ idx = np.where(PAobs < 0)
        #~ PAobs[idx] = PAobs[idx]+180
        line = ax4.errorbar(PAobs, \
        vis2.vis2data, yerr=vis2.vis2err, color=color, fmt='o', markersize=ms)#, label=label)
    #~ ax1.xaxis.get_major_formatter().set_useOffset(False)
    ax3.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    #~ ax3.get_xaxis().set_visible(False)
    #~ ax4.get_xaxis().set_visible(False)
    ax3.get_xaxis().set_ticklabels([])
    Blim = ax3.get_xlim()
    PAlim = ax4.get_xlim()
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
        u1 = t3.u1coord/t3.wavelength.eff_wave
        v1 = t3.v1coord/t3.wavelength.eff_wave
        u2 = t3.u2coord/t3.wavelength.eff_wave
        v2 = t3.v2coord/t3.wavelength.eff_wave
        B = np.append(np.sqrt(u1**2 + v1**2), np.sqrt(u2**2 + v2**2))
        PA = np.append(np.arctan2(u1,v1)*180.0/np.pi, np.arctan2(u2,v2)*180.0/np.pi)
        #~ idx = np.where(PA < 0)
        #~ PA[idx] = PA[idx]+180
        #~ if (t3.station[0] and t3.station[1]):
            #~ label = t3.station[0].sta_name + t3.station[1].sta_name
        #~ else:
            #~ label = 'unnamed'
        #~ color = names.index(label)
        #~ color = phc.colors[colorid]
        #~ colorid = np.mod(colorid + 1, len(phc.colors))
        color = 'Black'
        y = np.repeat(t3.t3phi, len(B)/len(t3.wavelength.eff_wave))
        yerr = np.repeat(t3.t3phierr, len(B)/len(t3.wavelength.eff_wave))
        line = ax5.errorbar(B, y, yerr=yerr, color=color, fmt='o', markersize=ms)#, label=label)
        line = ax6.errorbar(PA, y, yerr=yerr, color=color, fmt='o', markersize=ms)#, label=label)
    #~ ax5.get_xaxis().set_ticklabels([])
    ax5.set_xlim(Blim)
    ax6.set_xlim(PAlim)
    ymax = np.max(np.abs(ax5.get_ylim()))
    ax5.set_ylim([-1.05*ymax,1.05*ymax])
    ax6.set_ylim([-1.05*ymax,1.05*ymax])
    ax5.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax5.set_ylabel(u'$\phi_{123}$ (deg.)')
    ax5.grid(b=True, linestyle=':', alpha=alp)
    ax6.set_xlabel(u'$PA$ (deg.)')
    ax6.set_ylabel(u'$\phi_{123}$ (deg.)')
    ax6.grid(b=True, linestyle=':', alpha=alp)
    #SAVING
    dir, name = phc.trimpathname(ffile)
    name = phc.rmext(name)
    #plt.savefig('hdt/{}_{}.png'.format(hdrinfo[0], hdrinfo[2]), transparent=True)
    #plt.locator_params(axis = 'x', nbins = 7)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.96, bottom=0.09, hspace=.009, wspace=.32)
    for suf in fmt:
        plt.savefig('{0}.{1}'.format(name,suf), transparent=True)
    plt.close()
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
        dateobs2 = phc.strrep(dateobs, -3, ':')
        dateobs2 = phc.strrep(dateobs2, -6, ':')
        specfile = ffile[:ffile.find('_PRO/PRODUCT_')].replace('_SPEC','')+\
        '/*{0}*_OIDATA_AVG.fits*'.format(dateobs)
        specfile = glob(specfile)
        if len(specfile) == 0:
            specfile = ffile[:ffile.find('_PRO/PRODUCT_')].replace('_SPEC','')+\
        '/*{0}*_OIDATA_AVG.fits*'.format(dateobs2)
            specfile = glob(specfile)
        if len(specfile) != 1:
            specfile = ''
            print('# ERROR! This is a PRO oifits and the AVG file was '+\
            'not found!')
            print(ffile[:ffile.find('_PRO/PRODUCT_')]+('/*{0}*_OIDATA_'+\
        'AVG.fits*').format(dateobs))
        else:
            specfile = specfile[0]
            specoidata = oifits.open(specfile, quiet=True)
            #print('# {0} file read!!!'.format(specfile))
            oidata.amberspec = specoidata.amberspec
            spec = oidata.amberspec[0]
            vis = oidata.vis[0]
            if len(spec.wavelength.eff_wave) == len(vis.wavelength.eff_wave):
                for i in range(len(oidata.vis)):
                    oidata.amberspec[i].wavelength = oidata.vis[i].wavelength
            else:
                print('# ERROR! spec and vis sizes are different for {0}'.\
                format(phc.trimpathname(ffile)[1]))
    #
    fig = plt.figure(figsize=(5.6,8))
    alp = .75
    #xloc = plt.MaxNLocator(6)
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
            xmin=np.amin(np.append(1e6*vis.wavelength.eff_wave[\
            np.where(vis.flag == False)], xmax))
            xmax=np.amax(np.append(1e6*vis.wavelength.eff_wave[\
            np.where(vis.flag == False)], xmin))
            xrange = (xmin,xmax)
        u = vis.ucoord
        v = vis.vcoord
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        colorid = np.mod(colorid+1, len(colors))
        #[u,-u] = W > E
        #[-u,u] = E < W
        ax2.plot([-u,u],[v,-v], '.', label=label, color=color)
        ax2.xaxis.tick_top()
        names.append(vis.target.target)
    ax2.axis('equal')
    plt.grid(b=True, linestyle=':', alpha=alp)
    if legend: ax2.legend(prop={'size':8},numpoints=1,bbox_to_anchor=(-0.25, 1.0))
    #ax3 = dif.phases, RIGHT COLUMN
    plotid = 524
    names = []
    colorid = 0
    yrange = [0,0]
    for vis in oidata.vis:
        diff = np.max(np.abs(vis.visphi))
        if diff > yrange[1]:
            yrange = [-diff, diff]
    for vis in oidata.vis:
        ax3=fig.add_subplot(plotid)
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        colorid = np.mod(colorid + 1, len(colors))
        line = ax3.errorbar(1e6*vis.wavelength.eff_wave, vis.visphi,\
        vis.visphierr, label=label, color=color)
        ax3.set_ylim(yrange)
        ax3.set_xlim(xrange)
        ax3.set_ylabel(u'$\phi$ (deg.)')
        names.append(vis.target.target)
        names = list(np.unique(names))
        #title = names.pop()
        #for name in names:
        #    title += ', %s'%(name)
        #ax1.set_title(title)
        #ax1.set_ylabel('Differential phase')
        plotid+=2
        #ax3.get_xaxis().set_visible(False)
        ax3.get_xaxis().set_ticklabels([])
        plt.grid(b=True, linestyle=':', alpha=alp)
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
        #diff = np.max(np.abs(t3.t3phi))
        #yrange2 = [-diff, diff]
        ax5.set_ylim(yrange)
        names.append(vis.target.target)
        names = list(np.unique(names))
    ax5.set_ylabel(u'Closure $\phi$ (deg.)')
    ax5.set_xlabel('Wavelength ($\mu$m)')
    plt.setp( ax5.xaxis.get_majorticklabels(), rotation=-35 )
    plt.grid(b=True, linestyle=':', alpha=alp)
    #ax4 = visibilities, LEFT COLUMN
    plotid = 523
    names = []
    colorid = 0
    yrange = [1,0]
    for vis in oidata.vis2:
        yrange = [ np.min([yrange[0],np.min(vis.vis2data)]),\
        np.max([yrange[1],np.max(vis.vis2data)]) ]
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
        colorid = np.mod(colorid + 1, len(colors))
        line = ax4.errorbar(1e6*vis.wavelength.eff_wave, vis.vis2data,\
        vis.vis2err, label=label, color=color)
        names.append(vis.target.target)
        names = list(np.unique(names))
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
        plt.grid(b=True, linestyle=':', alpha=alp)
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
            y = linfit(1e6*spec.wavelength.eff_wave, spec.spectrum)
            ax1.plot(x,y, label=label, color=color)
            colorid = np.mod(colorid+1, len(colors))
            #[-u,u] = W > E
            #[u,-u] = E < W
        ax1.set_xlabel('Wavelength ($\mu$m)')
        ax1.set_xlim(xrange)
        ax1.set_ylabel('Norm. flux')
        plt.setp( ax1.xaxis.get_majorticklabels(), rotation=-35 )
        plt.grid(b=True, linestyle=':', alpha=alp)
    phc.outfld()
    dir, name = phc.trimpathname(ffile)
    name = phc.rmext(name)
    #plt.savefig('hdt/{}_{}.png'.format(hdrinfo[0], hdrinfo[2]), transparent=True)
    #plt.locator_params(axis = 'x', nbins = 7)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.96, bottom=0.09, hspace=.009, wspace=.32)
    for suf in fmt:
        plt.savefig('hdt/{0}.{1}'.format(name,suf), transparent=True)
    plt.close()
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
        path = os.getcwd()
    sh = glob(path+'/*.sh')
    f0 = open(sh[0])
    lines = f0.readlines()
    f0.close()
    count = 0
    count = len(glob(path+'/*.Z'))
    count+= len(glob(path+'/*.txt'))
    count+= len(glob(path+'/notused/*'))
    print(path,len(lines),count)
    return

def printinfo(file):
    """ Print OIFITS observational log, as appendix of Faes, D. M. (2015)

    DATE-OBS    MJD PA  B   PA  B   PA  B"""
    oidata = oifits.open(file, quiet='True')
    info = list(oidata.hdrinfo.returninfo())
    info2 = []
    for vis in oidata.vis:
        info2+= ['{0:.1f}'.format(np.sqrt(vis.ucoord**2 + vis.vcoord**2))]
        info2+= ['{0:.1f}'.format(np.arctan(vis.ucoord / vis.vcoord) * 180.0 / np.pi % 180.0)]
    return [info[2][:10], '{0:.7f}'.format(info[1])] + info2 

### MAIN ###
if __name__ == "__main__":
    pass
