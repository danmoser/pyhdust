# -*- coding:utf-8 -*-

"""PyHdust *images* module: auxiliary images module

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

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import os as _os
import numpy as _np
import time as _time
import pyhdust as _hdt
from pyhdust.hdrpil import hdr as _hdr
import pyhdust.phc as _phc
import warnings as _warn

try:
    import PIL as _PIL
except ImportError:
    _warn.warn('PIL (Pillow) module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def rgb2cmyk(img):
    """ `img` must have (m,n,3) dimensions

    R' = R/255

    G' = G/255

    B' = B/255

    The black key (K) color is calculated from the red (R'), green (G') and 
    blue (B') colors:

    K = min(R', G', B')

    The cyan color (C) is calculated from the red (R') and black (K) colors:

    C = (1-R'-K) / (1-K)

    The magenta color (M) is calculated from the green (G') and black (K) 
    colors:

    M = (1-G'-K) / (1-K)

    The yellow color (Y) is calculated from the blue (B') and black (K) colors:

    Y = (1-B'-K) / (1-K)
    """
    r = img[..., 0] / 255.
    g = img[..., 1] / 255.
    b = img[..., 2] / 255.
    # 
    k = 1 - _np.max((r, g, b), axis=0)
    c = (1 - r - k) / (1 - k)
    m = (1 - g - k) / (1 - k)
    y = (1 - b - k) / (1 - k)
    # 
    cmyk = _np.dstack((c, m, y, k)) * 255
    return cmyk.astype('uint8')


def doColorConv(lbd, flx, fdat=None):
    """ `fdat` is and (3,nlbd) matrix with the color space functions. If `None`,
        CIE 1931 RGB color space is used.

        INPUT = lambda (nm), flux (intensities) vectors

        OUTPUT = RGB values (0-255; uint8)
    """
    if fdat is None:
        fdat = _np.loadtxt(
            _hdt.hdtpath() + 'refs/rgb_eff.txt', unpack=True)
    idx = _np.where((lbd >= fdat[0, 0]) & (lbd <= fdat[0, -1]))
    lbd = lbd[idx]
    flx = flx[idx]
    out = _np.zeros(3)
    for i in range(0 + 1, 3 + 1):
        finterp = _np.interp(lbd, fdat[0], fdat[i])
        out[i - 1] = _np.trapz(finterp * flx, lbd)
    out = out * 255. / _np.max(out)
    return out.astype('uint8')


def doColorConvMaps(data, lbdc, fdat=None, di=0, oi=0, ii=0, outflt=False):
    """ `fdat` is and (3,nlbd) matrix with the color space functions. If `None`,
        CIE 1931 RGB color space is used.

        `di, oi, ii` are data[oi, ii, di]. `di` = 0 for *.maps files.

        outflt=True, the output is in float format.

        INPUT = data (n,m,nlbdc), lambda (nlbdc; um unit)

        OUTPUT = RGB images (n,m,3; [0-255; uint8])
    """
    if fdat is None:
        fdat = _np.loadtxt(
            _hdt.hdtpath() + 'refs/rgb_eff.txt', unpack=True)
    # 
    lbd = lbdc * 1e3
    idx = _np.where((lbd >= fdat[0, 0]) & (lbd <= fdat[0, -1]))
    lbd = lbd[idx]
    ndata = data[ii, oi, idx[0]]
    if len(_np.shape(ndata)) == 3:
        out = _np.empty(list(_np.shape(ndata[0, :, :])) + [3])
        ndata = _np.expand_dims(ndata, axis=-1)
        di = 0
    else:
        out = _np.empty(list(_np.shape(ndata[0, :, :, 0])) + [3])
    ndata = ndata[..., di]
    ndata = _np.swapaxes(ndata, 0, 1)
    ndata = _np.swapaxes(ndata, 1, 2)
    rgblbd = _np.zeros((len(lbd), 3))
    for i in range(3):
        rgblbd[:, i] = _np.interp(lbd, fdat[0], fdat[i + 1])
        out[:, :, i] = _np.trapz( rgblbd[:, i] * ndata, lbd)
    # 
    RGB = out * 255. / _np.max(out)
    if not outflt:
        RGB = RGB.astype('uint8')
    return RGB


def rgb2png(RGB, savename=None, ext='png', size=None):
    """ Save a array(n,m,3; uint8) as regular image file. """
    if savename is None:
        savename = _phc.dtflag()
    # 
    img = _PIL.Image.fromarray(RGB.astype('uint8'))
    if size is not None:
        img = img.resize(size, _PIL.Image.ANTIALIAS)
    # 
    ext = ext.replace('.', '')
    img.save(savename + '.' + ext)
    return


def doHDR(RGB, levels=6, folder='hdr', ext='png', strength=[1.],
        naturalness=[1.], size=None):
    """ Do HDR of RGB images. """
    ext = ext.replace('.', '')
    i = 0
    while _os.path.exists(folder) is True:
        folder += str(i)
    _os.system('mkdir ' + folder)
    for i in range(levels):
        img = _np.array(RGB, dtype='float') * 2**i
        img = _np.clip(img, 0, 255).astype('uint8')
        rgb2png(img, savename='{}/hdr_{}'.format(folder, i), size=size, 
            ext=ext)
    # 
    proj = _hdr(case=folder, img_type=ext, cur_dir='./')
    proj.get_hdr(strength=strength, naturalness=naturalness)
    return


def do_background(foreimg, backimg, savename=None, pos=[.1, .9], cut=0.5, 
    fmt=['png'], bga=.9, rotang=0., fsize=None):
    """ Add `backimg` as background (bg) of `foreimg` (fg). `*img` are regular
    image files (e.g., JPG or PNG).

    `pos` is the position of the central pixel of fg as fraction size of the 
    bg image. Example: 
    pos = [.1, .1] puts the foreground image in the lowe left position; 
    pos = [.5, .5] in the center. 

    `cut` = limit the brightest pixels. `cut == 1` says that only white points 
    will have no background (transparency. In other words, all fg will be 
    transparent as bg). 
    `cut == 0.5` says that every pixel with >= 50% of the white level will not  
    be transparent. `cut == 0` says that all black points will not be 
    transparent (i.e., no bg).

    `bga` = background strength

    METHOD: There the fg is brillhant (star), no bg. Where it is tenous, 
    (disk), it is applied a combination of bg + fg with the transparency 
    is proportiona to the fg brightness. 

    Is it important to note that even if a "black background" is applied, the 
    disk will be atenuated where it is below the cutting limit. 

    DETAILS: a combination of images, without saturation (standard level) is 
    equal to (img1 + img2 )/2.

    Problem 1: the star becomes transparent! Solution: we define a `cut`, where 
    pixels above a given level do not enter in the above sum, rather their 
    foreground values are kept.

    Problem 2: this procedure can create a big contrast between the star and 
    the disk, making the disk 50% fainter. Solucion: define the transparency 
    as function of the fg value (bigger the value, less transparent).

    Problem 3: if the star is seen edge-on, one can have bg where it is the 
    star! How to solve it? Analogous to the `cut`, if an mask can be created 
    if the photospheric image is provided (and then, no gb where px_val > 0).

    TODO: photospheric mask to cut.
    """
    if cut <= 0:
        cut = 1e-3
    if savename is None:
        savename = _time.strftime("%y%m%d-%H%M%S")
    # 
    # if fwd is not int(round(bwd/3.)):
        # fact = bwd/3./fwd
        # foreground = foreground.resize(np.round(foreground.size*fact).\
            # astype('int')), _PIL.Image.ANTIALIAS)
    # 
    aB = _PIL.Image.open(backimg)
    aF = _PIL.Image.open(foreimg)
    if fsize is not None:
        minbd = _np.min(aB.size)
        minfd = _np.min(aF.size)
        f = minbd*fsize/minfd
        fsize = ( _np.int(_np.round(f*aF.size[0])), _np.int(_np.round(
            f*aF.size[1])) )
        aF = aF.resize(fsize, _PIL.Image.ANTIALIAS)
    aB = _np.asarray(aB, dtype=float)
    aF = _np.asarray(aF, dtype=float)
    if rotang != 0:
        aF = rotate_cube(aF, rotang)
    aF = align_backfore(aB, aF, pos=pos)
    # tF = sum + clip in the limit (i.e., `cut` the brightest parts)
    # cut = mask to transparency (based on the values of the px, says if it 
    # will be transparent.
    tF = _np.sum(aF, axis=2)
    tF = _np.clip(tF, 0, 3 * 256 * cut)
    tF /= (3 * 256 * cut)
    # tB = inverse of tF, with the strength of `bga`
    tB = -1 * tF + bga
    tB = _np.clip(tB, 0, _np.max(tB))
    out = aF[..., :3] * tF[..., _np.newaxis] + \
        aB[..., :3] * tB[..., _np.newaxis]
    img = _PIL.Image.fromarray(out.astype('uint8'))
    for f in fmt:
        f = f.replace('.', '')
        img.save(savename + '.' + f)
    return


def align_backfore(bkimg, frimg, pos=[.1, .9]):
    """ doc """
    bkh, bkv = _np.shape(bkimg)[:2]
    f0h, f0v = _np.shape(frimg)[:2]
    frh, frv = _np.shape(bkimg)[:2]
    if len(_np.shape(frimg)) > 3:
        # TODO: Resige bkimg to the size of frimg
        print('# ERROR! foreground image wrong format (dimensions)!!')
        return None
    elif f0h == bkh and f0v == bkv and pos == [.5, .5]:
        # Nothing to do!
        print('# Nothing done!!')
        return frimg
    # Starting...
    if pos[0] < 0:
        pos[0] = 0
    if pos[1] > 1:
        pos[1] = 1
    # Horizontal
    ival = int(round(pos[0]*bkh))
    x0 = ival - f0h/2
    addx = 0
    if x0 < 0:
        x0 = 0
        addx += f0h/2 - ival
    elif bkh < x0 + f0h:
        addx += x0 + f0h - bkh   
    # vertical
    ival = int(round(pos[1]*bkv))
    y0 = ival - f0v/2
    addy = 0
    if y0 < 0:
        y0 = 0
        addy += f0v/2 - ival
    elif bkv < y0 + f0v:
        addy += y0 + f0v - bkv   

    if len(_np.shape(frimg)) == 3:
        d3 = _np.shape(frimg)[2]
        frtmp = _np.zeros((bkh+addx, bkv+addy, d3))
        # for i in range(d3):
        #     frtmp[x0:x0+f0h, y0:y0+f0v, i] = frimg[:, :, i]
    else:
        frtmp = _np.zeros((bkh+addx, bkv+addy))

    print(bkh, bkv, f0h, f0v, x0, y0, addx, addy, pos)
    frtmp[x0:x0+f0h, y0:y0+f0v] = frimg
    frnew = frtmp[:bkh, :bkv]
    # print frnew.shape, bkimg.shape
    return frnew


def rotate_coords(x, y, theta, ox, oy):
    """Rotate arrays of coordinates x and y by theta radians about the
    point (ox, oy).

    This routine was inspired on a http://codereview.stackexchange.com post.
    """
    s, c = _np.sin(theta), _np.cos(theta)
    x, y = _np.asarray(x) - ox, _np.asarray(y) - oy
    return x * c - y * s + ox, x * s + y * c + oy


def rotate_image(src, theta, ox=None, oy=None, fill=0):
    """Rotate the image src by theta radians about (ox, oy).
    Pixels in the result that don't correspond to pixels in src are
    replaced by the value fill.

    This routine was inspired on a http://codereview.stackexchange.com post.
    """
    # Images have origin at the top left, so negate the angle.
    theta = -theta

    # Dimensions of source image. Note that scipy.misc.imread loads
    # images in row-major order, so src.shape gives (height, width).
    sh, sw = src.shape
    if ox is None:
        ox = sw/2
    if oy is None:
        oy = sh/2

    # Rotated positions of the corners of the source image.
    cx, cy = rotate_coords([0, sw, sw, 0], [0, 0, sh, sh], theta, ox, oy)

    # Determine dimensions of destination image.
    dw, dh = (int(_np.ceil(c.max() - c.min())) for c in (cx, cy))

    # Coordinates of pixels in destination image.
    dx, dy = _np.meshgrid(_np.arange(dw), _np.arange(dh))

    # Corresponding coordinates in source image. Since we are
    # transforming dest-to-src here, the rotation is negated.
    sx, sy = rotate_coords(dx + cx.min(), dy + cy.min(), -theta, ox, oy)

    # Select nearest neighbour.
    sx, sy = sx.round().astype(int), sy.round().astype(int)

    # Mask for valid coordinates.
    mask = (0 <= sx) & (sx < sw) & (0 <= sy) & (sy < sh)

    # Create destination image.
    dest = _np.empty(shape=(dh, dw), dtype=src.dtype)

    # Copy valid coordinates from source image.
    dest[dy[mask], dx[mask]] = src[sy[mask], sx[mask]]

    # Fill invalid coordinates.
    dest[dy[~mask], dx[~mask]] = fill

    return dest


def rotate_cube(src, theta, ox=None, oy=None, fill=0):
    """ `src` must be a numpy array. """
    dims = src.shape
    if len(dims) == 2:
        print('# Warning! It is a image, not a cube!')
        return rotate_image(src, theta, ox, oy, fill)

    tmp = rotate_image(src[:, :, 0], theta, ox, oy, fill)
    dimsr = tmp.shape
    rot = _np.zeros((dimsr[0], dimsr[1], dims[2]))
    rot[:, :, 0] = tmp
    for i in range(1, dims[2]):
        rot[:, :, i] = rotate_image(src[:, :, i], theta, ox, oy, fill)
    return rot


def genSquare(size=64, halfside=16, center=(0, 0)):
    """
    Generate a square inside a square.

    If size is not even, the unit square will not be centered.

    center is the relative position
    """
    x = _np.zeros(size)
    y = x[:, _np.newaxis]
    # ABSOLUTE Center:
    # if center is None:
    #    x0 = y0 = size // 2
    # else:
    #    x0 = center[0]
    #    y0 = center[1]
    # Relative Center:
    x0 = size // 2 + center[0]
    y0 = size // 2 + center[1]
    img = x * y
    img[x0 - halfside:x0 + halfside, y0 - halfside:y0 + halfside] = 1
    return img


def genGaussian(size=64, sig=64 / 8, center=(0, 0)):
    """
    Generate a square gaussian kernel (non-normalized).

    `size` is the length of a side of the square (pixels)

    center is the relative position
    """
    x = _np.arange(0, size, 1, float)
    y = x[:, _np.newaxis]
    # ABSOLUTE Center:
    # if center is None:
    #    x0 = y0 = size // 2
    # else:
    #    x0 = center[0]
    #    y0 = center[1]
    # Relative Center:
    x0 = size // 2 + center[0]
    y0 = size // 2 + center[1]
    #
    return _np.exp(-0.5 * ((x - x0)**2 + (y - y0)**2) / sig**2)


# MAIN ###
if __name__ == "__main__":
    pass
