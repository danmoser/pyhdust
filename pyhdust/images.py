# -*- coding:utf-8 -*-

"""
PyHdust *images* module: auxiliary images module

# ramp yellow
sz = 100
r = _np.arange(sz*sz).reshape((sz,sz))*256./sz/sz
g = _np.arange(sz*sz).reshape((sz,sz))*256./sz/sz
b = _np.zeros((sz,sz))
yramp = _np.dstack((r,g,b)).astype(uint8)

img = PIL.Image.fromarray(yramp)
img.save('yrampRGB.jpg')

img = PIL.Image.fromarray(rgb2cmyk(yramp), mode='CMYK')
img.save('yrampCMYK.jpg')

# convert 'CMYK_color'
CM = PIL.Image.open('CMYK_color.jpg')
CM = CM.convert(mode='RGB')
CM.save('CMYK2rgbPIL.jpg')

img = PIL.Image.fromarray(rgb2cmyk(_np.asarray(CM)), mode='CMYK')
img.save('CMYK2rgb.jpg')

:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""
import os as _os
import numpy as _np
import time as _time
import pyhdust as _hdt
from pyhdust.hdrpil import hdr as _hdr

try:
    import PIL as _PIL
except:
    print('# Warning! PIL (Pillow) module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def rgb2cmyk(img):
    """ `img` must have (m,n,3) dimensions

    R' = R/255
    
    G' = G/255
    
    B' = B/255
    
    The black key (K) color is calculated from the red (R'), green (G') and blue (B') colors:
    
    K = min(R', G', B')
    
    The cyan color (C) is calculated from the red (R') and black (K) colors:
    
    C = (1-R'-K) / (1-K)
    
    The magenta color (M) is calculated from the green (G') and black (K) colors:
    
    M = (1-G'-K) / (1-K)
    
    The yellow color (Y) is calculated from the blue (B') and black (K) colors:
    
    Y = (1-B'-K) / (1-K)
    """
    r = img[...,0]/255.
    g = img[...,1]/255.
    b = img[...,2]/255.
    #~ 
    k = 1-_np.max((r, g, b), axis=0)
    c = (1-r-k)/(1-k)
    m = (1-g-k)/(1-k)
    y = (1-b-k)/(1-k)
    #~
    cmyk = _np.dstack((c,m,y,k))*255
    return cmyk.astype('uint8')


def doColorConv(lbd, flx, fdat=None):
    """ `fdat` is and (3,nlbd) matrix with the color space functions. If `None`,
        CIE 1931 RGB color space is used.

        INPUT = lambda (nm), flux (intensities) vectors

        OUTPUT = RGB values (0-255; uint8)
    """
    if fdat is None:
        fdat = _np.loadtxt(_hdt.hdtpath()+'/pyhdust/refs/rgb_eff.txt', unpack=True)
    idx = _np.where((lbd >= fdat[0, 0]) & (lbd <= fdat[0, -1]))
    lbd = lbd[idx]
    flx = flx[idx]
    out = _np.zeros(3)
    for i in range(0+1,3+1):
        finterp = _np.interp(lbd, fdat[0], fdat[i])
        out[i-1] = _np.trapz(finterp * flx, lbd)
    out = out*255./_np.max(out)
    return out.astype('uint8')


def doColorConvMaps(data, lbdc, fdat=None, di=0, oi=0, ii=0, outflt=False):
    """ `fdat` is and (3,nlbd) matrix with the color space functions. If `None`,
        CIE 1931 RGB color space is used.

        `di, oi, ii` are data[oi, ii, di]. `di` = 0 for *.maps files.

        INPUT = data (n,m,nlbdc), lambda (nlbdc; um unit)

        OUTPUT = RGB images (n,m,3; [0-255; uint8])
    """
    if fdat is None:
        fdat = _np.loadtxt(_hdt.hdtpath()+'/pyhdust/refs/rgb_eff.txt', unpack=True)
    #~ 
    lbd = lbdc*1e3
    idx = _np.where((lbd >= fdat[0, 0]) & (lbd <= fdat[0, -1]))
    lbd = lbd[idx]
    ndata = data[ii,oi,idx[0]]
    if len(_np.shape(ndata)) == 3:
        out = _np.empty(list(_np.shape(ndata[0,:,:]))+[3])
        ndata = _np.expand_dims(ndata,axis=-1)
        di = 0
    else:
        out = _np.empty(list(_np.shape(ndata[0,:,:,0]))+[3])
    ndata = ndata[...,di]
    ndata = _np.swapaxes(ndata, 0, 1)
    ndata = _np.swapaxes(ndata, 1, 2)
    rgblbd = _np.zeros((len(lbd),3))
    for i in range(3):
        rgblbd[:,i] = _np.interp(lbd, fdat[0], fdat[i+1])
        out[:,:,i] = _np.trapz( rgblbd[:,i] * ndata, lbd)
    #~ 
    RGB = out*255./_np.max(out)
    if not outflt:
        RGB = RGB.astype('uint8')
    return RGB


def rgb2png(RGB, savename=None, ext='png', size=None):
    """ Save a array(n,m,3; uint8) as regular image file. """
    if savename is None:
        savename = _time.strftime("%y%m%d-%H%M%S")
    #~ 
    img = _PIL.Image.fromarray(RGB.astype('uint8'))
    if size is not None:
        img = img.resize(size, _PIL.Image.ANTIALIAS)
    #~
    ext = ext.replace('.','')
    img.save(savename+'.'+ext)
    return


def doHDR(RGB, levels=6, folder='hdr', ext='png', strength=[1.],
        naturalness=[1.]):
    """ Do HDR of RGB images. """
    ext = ext.replace('.','')
    i = 0
    while _os.path.exists(folder) == True:
        folder+= str(i)
    _os.system('mkdir '+folder)
    for i in range(7):
        img = _np.array(RGB, dtype='float')*2**i
        img = _np.clip(img, 0, 255).astype('uint8')
        rgb2png(img, size=[1024,1024], savename='{}/hdr_{}'.format(folder,i),
            ext=ext)
    #~ 
    proj = _hdr(case=folder, img_type=ext, cur_dir='./')
    proj.get_hdr(strength=strength,naturalness=naturalness)
    return


def doBackground(foreimg, backimg, savename=None, pos='1', cut=0.5, ext='png',
    bga=.9):
    """ add ``backimg`` as background of ``foreimg``. ``*img`` are regular
    image files.

    ``pos`` equals the dialpad (without 0). """
    if savename is None:
        savename = _time.strftime("%y%m%d-%H%M%S")
    ext = ext.replace('.','')
    #~ 
    background = _PIL.Image.open(backimg)
    bwd = background.size[0]
    foreground = _PIL.Image.open(foreimg)
    fwd = foreground.size[0]
    #~ if fwd is not int(round(bwd/3.)):
        #~ fact = bwd/3./fwd
        #~ foreground = foreground.resize(np.round(foreground.size*fact).\
            #~ astype('int')), _PIL.Image.ANTIALIAS)
    #~ 
    aB = _PIL.Image.open(backimg)
    aF = _PIL.Image.open(foreimg)
    #~ 
    aB = _np.asarray(aB, dtype=float)
    aF = _np.asarray(aF, dtype=float)
    tF = _np.sum(aF, axis=2)
    tF = _np.clip(tF, 0, 3*256*cut)
    tF/= (3*256*cut)
    tB = -1*tF+bga
    tB = _np.clip(tB, 0, _np.max(tB))
    out = aF*tF[...,_np.newaxis] + aB[...,:3]*tB[...,_np.newaxis]
    img = _PIL.Image.fromarray(out.astype('uint8'))
    img.save(savename+'.'+ext)
    return
