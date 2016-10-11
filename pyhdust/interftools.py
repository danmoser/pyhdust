# -*- coding:utf-8 -*-

"""PyHdust *interftools* module: interferometry tools

`colors` keep the *amdlib* standard.

A biblioteca python XDRLIB eh MUITO lenta... Usa muitas listas!!!

>>> import xdrlib

A biblioteca PYDAP estah em desenvolvimento... Eh complicada de usar

>>> from pydap.model import *
>>> from pydap.xdr import DapUnpacker
>>> base_int = BaseType(name='base_int')
>>> base_float = BaseType(name='base_float', type=Float32)

Todas as leituras binarias baseiam-se no struct.

- Rotines based on images (np.array matrix)
- Rotines based on oifiles
- Rotines based on HDUST (maps/data)
- Other?

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import os as _os
import struct as _struct
import numpy as _np
from glob import glob as _glob
import pyhdust.phc as _phc
import pyhdust.tabulate as _tab
import pyhdust.oifits as _oifits
import pyhdust.images as _img
from pyhdust.spectools import linfit as _linfit
import copy as _copy
from six import string_types as _strtypes
import warnings as _warn

try:
    import matplotlib as _mpl
    import matplotlib.pyplot as _plt
    import matplotlib.gridspec as _gridspec
    import matplotlib.ticker as _mtick
    from matplotlib.lines import Line2D as _Line2D
    import pyfits as _pyfits
except ImportError:
    _warn.warn('matplotlib, six and/or pyfits module not installed!!!')

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


colors = ["red", "green", "blue", "black"]


def azim_avg(img):
    """ Azimuthaly average a given image.

    .. note:: `round=0`. If it bigger than 0, then rounding problem appears.
    """
    img = img.copy()
    nrows, ncols = img.shape
    xx, yy = _np.meshgrid( *_np.ogrid[:ncols, :nrows] )
    xx = xx-(nrows-1)/2.
    yy = yy-(ncols-1)/2.
    rr = -_np.sqrt( xx**2 + yy**2 ).round(0)
    for r in _np.unique(rr):
        idx = _np.where(rr == r)
        img[idx] = _np.average(img[idx])
    return img


def log_transform(im):
    """ Returns log(image) scaled to the interval [0,1]
    """
    try:
        (imin, imax) = (_np.min(im[_np.where(im > 0)]), _np.max(im))
        if (imax > imin) and (imax > 0):
            im = (_np.log(im.clip(imin, imax))-_np.log(imin)) / (
                _np.log(imax)-_np.log(imin))
            idx = _np.where(im == 0)
            im[idx] = _np.NaN
            return im
    except:
        _warn.warn('# Warning! Wrong input at log_transform !!!')
    return im


def imshowl(img, cmap='gist_heat', origin='lower'):
    """ Plot the normalized image in log-scale.
    """
    _plt.clf()
    _plt.imshow(log_transform(img), cmap=_plt.get_cmap(cmap), origin=origin)
    return


def readmap(file, quiet=False):
    r"""
    Read *Hdust* MAP or MAPS files.

    `mapimg`: extract this component from the *.map* file.

        - 0 = total flux
        - 1 = transmitted flux
        - 2 = scattered flux
        - 3 = emitted flux
        - 4 = pol. *Q* flux
        - 5 = pol. *U* flux
        - 6 = pol. *V* flux

    OUTPUT = data, obslist, lbdc, Ra, xmax

        - data = image matrix
        - obslist = observers info (*i*, :math:`\phi`)
        - lbdc = central :math:`\lambda`
        - Ra = ?
        - xmax = image size in Rsun untis

    | data(nimgs,nobs,nlbd,ny,nx,dfact)
    | .map, dfact = 6
    | .maps, dfact = 1
    """
    if file[-4:] == '.map':
        dfact = 6
    elif file[-5:] == '.maps':
        dfact = 1
    else:
        print('# ERROR: This is not a HDUST valid image!')
        return
    f = open(file, 'rb').read()
    #
    ixdr = 0
    nobs, lnum, nx, ny = _struct.unpack('>4l', f[ixdr:ixdr + 4 * 4])
    ixdr += 4 * 4
    Ra, Rstar, Lratio, xmax = _struct.unpack('>4f', f[ixdr:ixdr + 4 * 4])
    ixdr += 4 * 4
    # nf no IDL estah como DOUBLE, mas com certea eh float ou int...
    # caso contrario nm nao faz sentido.
    nf = _struct.unpack('>f', f[ixdr:ixdr + 1 * 4])[0]
    ixdr += 1 * 4
    nm = _struct.unpack('>l', f[ixdr:ixdr + 1 * 4])[0]
    ixdr += 1 * 4

    npxs = nm
    upck = '>{}f'.format(npxs)
    xmax = _np.array( _struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) )
    ixdr += npxs * 4

    if file[-4:] == '.map':
        # tmp must be a "small" array. Otherwise, a MemoryError will be raised
        npxs = nx * ny * lnum * nobs * nm
        tmp = _np.empty((npxs * dfact))
        upck = '>{}f'.format(npxs)
        for i in range(dfact):
            # skip first npxs... See below
            tmp[i * npxs:(i + 1) * npxs] = _np.array( _struct.unpack(upck,
            f[ixdr:ixdr + npxs * 4]) )
            ixdr += npxs * 4
        data = _np.zeros((nm, nobs, lnum, ny, nx, dfact + 1))
        for i in range(dfact):
            data[:, :, :, :, :, i + 1] = tmp[i::dfact].reshape((
                nm, nobs, lnum, ny, nx))
        data[:, :, :, :, :, 0] = data[:, :, :, :, :, 1] +\
            data[:, :, :, :, :, 2] + data[:, :, :, :, :, 3]
        #
    elif file[-5:] == '.maps':
        npxs = dfact * nx * ny * lnum * nobs * nm
        upck = '>{}f'.format(npxs)
        data = _np.array( _struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) ).\
            reshape((nm, nobs, lnum, ny, nx))
        ixdr += npxs * 4

    npxs = 2 * nobs
    upck = '>{}f'.format(npxs)
    obslist = _np.array( _struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) )
    ixdr += npxs * 4

    npxs = lnum + 1
    upck = '>{}f'.format(npxs)
    lbdarr = _np.array( _struct.unpack(upck, f[ixdr:ixdr + npxs * 4]) )
    ixdr += npxs * 4

    # this will check if the XDR is finished.
    if ixdr == len(f):
        if not quiet:
            print('# XDR {} completely read!'.format(file))
    else:
        print('# Warning: XDR {} not completely read!'.format(file))
        print('# length difference is {}'.format( (len(f) - ixdr) / 4 ) )

    # lbdarr tem lnum+1, pois reflete o INTERVALO de cada imagem.
    # Para termos o lambda central de cada imagem, fazemos o seguinte:
    lbdc = _np.zeros(lnum)
    for i in range(lnum):
        lbdc[i] = (lbdarr[i] + lbdarr[i + 1]) / 2.

    return data, obslist, lbdc, Ra, xmax


def img2fits(img, lbd, xmax, dist, outname='model', rot=0., lum=0.,
    orient=0., coordsinf=None, deg=False, ulbd=''):
    r""" Export an image (e.g., data[0,0,0,:,:]) to the fits format.

    `lbd` is the wavelength value and the dimension is kept as it is. It must
    be in meters for JMMC softwares (ASPRO2/LITPRO).

    `ulbd` units of the lbd.

    `rot` = rotation angle to be applied to the images. 'x' and 'y' coordinate
    axes should be orientated with equatorial north corresponding to 'up' (and
    east == 'left'). Units according to `deg` bool.

    `orient` = orientation of the coordinate system. This is completely
    independent of the `rot` variable. Units according to `deg` bool.

    `lum` = luminosity given in Solar units.
    BUNIT sets the units of the image. Considering that HDUST images
    give the pixels counts as :math:`F_\lambda/F`, the same correction as done
    to BeAtlas is performed, and the final results are in 10^-17
    erg/s/cm^2/Ang. If `lum` = 0, no change is done.

    `coordsinf` = [RA,DEC], as ['21:51:12.055', '+28:51:38.72']

    `deg` = angles in degrees (instead of radians).

    Example: image at 21 cm, rotated 45 degrees, 2 AU long at 10 parsecs.

    .. code-block:: python

        img = np.arange(900).reshape((30,30))

        intt.img2fits(img, 21., [2*phc.au.cgs/phc.Rsun.cgs], 10, orient=45.,
        coordsinf=['21:51:12.055', '-28:51:38.72'], ulbd='cm', deg=True)

    .. image:: _static/modelfits.png
        :align: center

    .. code-block:: python

        intt.img2fits(img, 21., [2*phc.au.cgs/phc.Rsun.cgs], 10, rot=45.,
        coordsinf=['21:51:12.055', '-28:51:38.72'], ulbd='cm', deg=True,
        outname='model_rotated')

    .. image:: _static/modelrotfits.png
        :align: center
    """
    if deg:
        rot = rot * _np.pi / 180
        ucdelt = 'degrees'
        ushort = 'deg'
        orientr = orient * _np.pi / 180
    else:
        ucdelt = 'radians'
        ushort = 'rad'
        orientr = orient
        orient = orient * 180 / _np.pi
    #
    if rot != 0.:
        # print('# Total flux BEFORE rotation {0}'.format(_np.sum(img)))
        img = _img.rotate_image(img, rot, 0, 0, fill=0.)
        # print('# Total flux AFTER rotation {0}\n'.format(_np.sum(img)))
    if lum != 0:
        img = img * (lum * _phc.Lsun.cgs) / 4 / _np.pi / \
            (dist * _phc.pc.cgs)**2 * 1e-4 * 1e17
    #
    hdu = _pyfits.PrimaryHDU(img[::-1, :])
    hdulist = _pyfits.HDUList([hdu])
    pixsize = 2 * xmax[0] / len(img)
    ang_per_pixel = _np.double(pixsize * _phc.Rsun.cgs / (dist * _phc.pc.cgs))
    # *60.*60.*1000.*180./_np.pi)
    if deg:
        ang_per_pixel *= 180. / _np.pi
    #
    hdulist[0].header['CDELT1'] = (-ang_per_pixel, ucdelt)
    hdulist[0].header['CDELT2'] = (ang_per_pixel, ucdelt)
    hdulist[0].header['CDELT3'] = (1., ulbd)
    hdulist[0].header['CRVAL3'] = (lbd, ulbd)
    hdulist[0].header['NAXIS1'] = len(img)
    hdulist[0].header['NAXIS2'] = len(img[0])
    hdulist[0].header['CRPIX1'] = len(img) / 2.
    hdulist[0].header['CRPIX2'] = len(img[0]) / 2.
    hdulist[0].header['CPPIX1'] = len(img) / 2
    hdulist[0].header['CPPIX2'] = len(img[0]) / 2
    hdulist[0].header['CROTA2'] = (float('{0:.3f}'.format(orient)), 'degrees')
    if lum != 0:
        hdulist[0].header['BUNIT'] = '10^-17 erg/s/cm^2/Ang'
    if coordsinf is not None:
        hdulist[0].header['CTYPE1'] = 'RA---TAN'  # 'GLON-CAR'
        hdulist[0].header['CTYPE2'] = 'DEC--TAN'  # 'GLON-CAR'
        hdulist[0].header['RA'] = coordsinf[0]
        hdulist[0].header['DEC'] = coordsinf[1]
        # hdulist[0].header['RADECSYS'] = 'FK5'
        hdulist[0].header['EQUINOX'] = 2000.0
        hdulist[0].header['CRVAL1'] = (_phc.ra2degf(coordsinf[0]), 'degrees')
        hdulist[0].header['CRVAL2'] = (_phc.dec2degf(coordsinf[1]), 'degrees')
    else:
        hdulist[0].header['CRVAL1'] = 0.
        hdulist[0].header['CRVAL2'] = 0.
    hdulist[0].header['CROTA1'] = (0.000, 'degrees')
    hdulist[0].header['CD1_1'] = (-ang_per_pixel * _np.cos(orientr), ucdelt)
    hdulist[0].header['CD1_2'] = (-ang_per_pixel * _np.sin(orientr), ucdelt)
    hdulist[0].header['CD2_1'] = (ang_per_pixel * -_np.sin(orientr), ucdelt)
    hdulist[0].header['CD2_2'] = (ang_per_pixel * _np.cos(orientr), ucdelt)
    hdulist[0].header['CUNIT1'] = ushort
    hdulist[0].header['CUNIT2'] = ushort
    hdulist[0].header['LONPOLE'] = 180.000
    hdulist[0].header['LATPOLE'] = 0.000
    outname = '{0}.fits'.format(outname.replace(".fits", ""))
    hdu.writeto(outname, clobber=True)
    print('# Saved {0} !'.format(outname))
    return


def data2fitscube(data, obs, lbdc, xmax, dist, zoom=0, outname='model',
    orient=0., rot=0., lum=0., coordsinf=None, map=False, deg=False):
    r""" Export a set of images (e.g., data[zoom,obs,:,:,:]) to the fits cube
    format.

    `map` = if `data` is a *.map file, set it to True. Leave false to *.maps.

    `lbdc` is the wavelength array and the dimension is kept as it is. It must
    be in meters for JMMC softwares (ASPRO2/LITPRO).

    `rot` = rotation angle to be applied to the images. 'x' and 'y' coordinate
    axes should be orientated with equatorial north corresponding to 'up' (and
    east == 'left'). Units in Degrees.

    `orient` = orientation of the coordinate system. This is completely
    independent of the `rot` variable.

    `lum` = luminosity given in Solar units.
    BUNIT sets the units of the image. Considering that HDUST images
    give the pixels counts as :math:`F_\lambda/F`, the same correction as done
    to BeAtlas is performed, and the final results are in 10^-17
    erg/s/cm^2/Ang. If `lum` = 0, no change is done.

    `coordsinf` = [RA,DEC].
    Example: ['21:51:12.055', '+28:51:38.72']

    `deg` = angles in degrees (instead of radians).
    """
    if deg:
        rot = rot * _np.pi / 180
        ucdelt = 'degrees'
        ushort = 'deg'
        orientr = orient * _np.pi / 180
    else:
        ucdelt = 'radians'
        ushort = 'rad'
        orientr = orient
        orient = orientr * 180 / _np.pi
    ulbd = 'meters'
    #
    if not map:
        imgs = data[zoom, obs, :, ::-1, :]
    else:
        imgs = data[zoom, obs, :, ::-1, :, 0]
    #
    if rot == 0.:
        pass
    else:
        print('# ERROR! Rotation of CUBES not yet implemented!')
        raise SystemExit(1)
    if lum != 0:
        # iL = _phc.fltTxtOccur('L =', lines, seq=2)*_phc.Lsun.cgs
        imgs = imgs * (lum * _phc.Lsun.cgs) / 4 / _np.pi / (dist *
            _phc.pc.cgs)**2 * 1e-4 * 1e17
    #
    hdu = _pyfits.PrimaryHDU(imgs)
    hdulist = _pyfits.HDUList([hdu])
    pixsize = 2 * xmax[0] / len(data[zoom, obs, 0, :, :])
    ang_per_pixel = _np.double(pixsize * _phc.Rsun.cgs / (dist * _phc.pc.cgs))
    # *60.*60.*1000.*180./_np.pi)
    if deg:
        ang_per_pixel *= 180. / _np.pi
    #
    hdulist[0].header['CDELT1'] = (-ang_per_pixel, ucdelt)
    hdulist[0].header['CDELT2'] = (ang_per_pixel, ucdelt)
    hdulist[0].header['CDELT3'] = ((lbdc[-1] - lbdc[0]) / len(lbdc), ulbd)
    hdulist[0].header['CRVAL3'] = (lbdc[0], ulbd)
    hdulist[0].header['NAXIS1'] = len(imgs)
    hdulist[0].header['NAXIS2'] = len(imgs[0])
    hdulist[0].header['CRPIX1'] = len(imgs) / 2.
    hdulist[0].header['CRPIX2'] = len(imgs[0]) / 2.
    hdulist[0].header['CPPIX1'] = len(imgs) / 2
    hdulist[0].header['CPPIX2'] = len(imgs[0]) / 2
    hdulist[0].header['CROTA2'] = (float('{0:.3f}'.format(orient)), 'degrees')
    if lum != 0:
        hdulist[0].header['BUNIT'] = '10^-17 erg/s/cm^2/Ang'
    if coordsinf is not None:
        hdulist[0].header['CTYPE1'] = 'RA---TAN'  # 'GLON-CAR'
        hdulist[0].header['CTYPE2'] = 'DEC--TAN'  # 'GLON-CAR'
        hdulist[0].header['RA'] = coordsinf[0]
        hdulist[0].header['DEC'] = coordsinf[1]
        # hdulist[0].header['RADECSYS'] = 'FK5'
        hdulist[0].header['EQUINOX'] = 2000.0
        hdulist[0].header['CRVAL1'] = (_phc.ra2degf(coordsinf[0]), 'degrees')
        hdulist[0].header['CRVAL2'] = (_phc.dec2degf(coordsinf[1]), 'degrees')
    else:
        hdulist[0].header['CRVAL1'] = 0.
        hdulist[0].header['CRVAL2'] = 0.
    hdulist[0].header['CROTA1'] = (0.000, 'degrees')
    hdulist[0].header['CD1_1'] = (-ang_per_pixel * _np.cos(orientr), ucdelt)
    hdulist[0].header['CD1_2'] = (-ang_per_pixel * _np.sin(orientr), ucdelt)
    hdulist[0].header['CD2_1'] = (ang_per_pixel * -_np.sin(orientr), ucdelt)
    hdulist[0].header['CD2_2'] = (ang_per_pixel * _np.cos(orientr), ucdelt)
    hdulist[0].header['CUNIT1'] = ushort
    hdulist[0].header['CUNIT2'] = ushort
    hdulist[0].header['LONPOLE'] = 180.000
    hdulist[0].header['LATPOLE'] = 0.000
    outname = '{0}.fits'.format(outname.replace(".fits", ""))
    hdu.writeto(outname, clobber=True)
    print('# Saved {0} !'.format(outname))
    return


def setspacecoords(nx, ny, rad_per_pixel, xc=0., yc=0.):
    """
    return xx and yy, 2D physical coordinates in ANGULAR dimensions
    (unit = defined by 'rad_per_pixel', i.e., radians).
    The physical scale (or length) on both axis must be the same.

    xc and yc are the the center position in PIXELS
    """
    x = _np.arange(0., nx) - (nx - 1) / 2. + xc
    xx = _np.repeat(x, ny).reshape(-1, ny).T * rad_per_pixel
    y = _np.arange(0., ny) - (ny - 1) / 2. + yc
    yy = _np.repeat(y, nx).reshape(-1, nx) * rad_per_pixel
    return xx, yy


def fastnumvis(img, lbd, Bproj, PA, rad_per_pixel, PAdisk=90.):
    """
    For a given image (in phys.units = `rad_per_pixel`) and a interf. setup,
        it returns the visibility and phase.

    `PA` and `PAdisk` in degrees.

    `Bproj` and `lbd` must have the same units (m).

    output: complexVis, VisAmp, VisPhase
    """
    if lbd < 1e-6 or lbd > 4e-6:
        print('# Warning! *fastnumvis*(lbd) is {0:.1e} m!'.format(lbd))
    PA = PA - PAdisk + 90.
    idx = _np.where(img > 0)

    u = Bproj * _np.double(_np.sin(PA / _np.double(180. / _np.pi)) / lbd)
    v = Bproj * _np.double(_np.cos(PA / _np.double(180. / _np.pi)) / lbd)
    # print PA,phc.ra2deg,lbd,Bproj,v

    ny = len(img)
    nx = len(img[0])
    xx, yy = setspacecoords(nx, ny, rad_per_pixel, xc=0., yc=0.)

    arg = -2 * _np.pi * (xx[idx] * u + yy[idx] * v)
    TF_z_re = _np.sum(img[idx] * _np.cos(arg))
    TF_z_im = _np.sum(img[idx] * _np.sin(arg))
    # print TF_z_re,TF_z_im

    TF_z = complex(TF_z_re, TF_z_im)
    TF_z0 = _np.sum(img[idx])

    complexVis = TF_z / TF_z0

    VisAmp = _np.abs(complexVis)
    VisPhase = _np.arctan2(
        complexVis.imag, complexVis.real) * _np.double(180. / _np.pi)
    return complexVis, VisAmp, VisPhase


def fastnumvis3(img, lbd, Bprojs, PAs, rad_per_pixel, PAdisk=90.):
    """ For a given image (in phys.units = `rad_per_pixel`) and 3 telescopes
    interferometric setup, it returns the visibility and closure phase.

    Is 

    `PA` as array/list (two values) and `PAdisk` in degrees.

    `Bproj` as array/list (two values) and `lbd` must have the same units (m).

    output: complexVis, VisAmp, VisPhase (closure phase)
    """
    u1 = Bprojs[0] * _np.cos(PAs[0] * _np.pi / 180.)
    u2 = Bprojs[1] * _np.cos(PAs[1] * _np.pi / 180.)
    v1 = Bprojs[0] * _np.sin(PAs[0] * _np.pi / 180.)
    v2 = Bprojs[1] * _np.sin(PAs[1] * _np.pi / 180.)
    B3 = _np.sqrt( (u1 + u2)**2 + (v1 + v2)**2 )
    PA3 = _np.arctan2( u1 + u2, v1 + v2 ) * 180 / _np.pi

    cV1, VA1, VP1 = fastnumvis(
        img, lbd, Bprojs[0], PAs[0], rad_per_pixel, PAdisk=PAdisk)
    cV2, VA2, VP2 = fastnumvis(
        img, lbd, Bprojs[1], PAs[1], rad_per_pixel, PAdisk=PAdisk)
    cV3, VA3, VP3 = fastnumvis(img, lbd, B3, PA3, rad_per_pixel, PAdisk=PAdisk)

    complexVis = cV1 * cV2 * cV3.conjugate()

    VisAmp = _np.abs(complexVis)
    VisPhase = _np.arctan2(
        complexVis.imag, complexVis.real) * _np.double(180. / _np.pi)
    return complexVis, VisAmp, VisPhase


def reversePAs(PAs=[-180, 180]):
    """ Calculte the reverse of the PAs input range.

    OUTPUT: PArv (list, units deg)
    """ 
    PArv = [0, 0]
    if PAs[0] < 0:
        PArv[0] = PAs[0] + 180
    else:
        # PArv[0] == 0 must be -180!
        PArv[0] = PAs[0] - 180
    if PAs[1] <= 0:
        # PArv[1] == 0 must be +180!
        PArv[1] = PAs[1] + 180
    else:
        PArv[1] = PAs[1] - 180
    return PArv


class oiClass(object):
    """ PyHdust class for interferometric data. 
    """
    def __init__(self, oitype, oifile=None):
        self.oitype = oitype
        self.oifile = oifile

    def readvis2(self, oidata_vis2):
        tlist = []
        for v2i in oidata_vis2:
            n = len(v2i.wavelength.eff_wave)
            if (v2i.station[0] and v2i.station[1]):
                tlist.append([v2i.station[0].sta_name + 
                    v2i.station[1].sta_name]*n)
            else:
                tlist.append(['unnamed']*n)
            u = v2i.ucoord
            v = v2i.vcoord
            tlist.append([_np.sqrt(u**2 + v**2)]*n)
            tlist.append([_np.arctan2(u, v) * 180 / _np.pi]*n)
            tlist.extend([v2i.wavelength.eff_wave])
            tlist.extend([_np.array(v2i.vis2data)])
            tlist.extend([_np.array(v2i.vis2err)])
        # print tlist
        # tlist = _np.array(tlist)
        # print tlist[3::6]
        # raise SystemError
        self.v2t = _np.array(_phc.flatten(tlist[0::6]))
        self.v2B = _np.array(_phc.flatten(tlist[1::6])).astype(float)
        self.v2PA = _np.array(_phc.flatten(tlist[2::6])).astype(float)
        self.v2lb = _np.array(_phc.flatten(tlist[3::6])).astype(float)
        self.v2 = _np.array(_phc.flatten(tlist[4::6])).astype(float)
        self.v2e = _np.array(_phc.flatten(tlist[5::6])).astype(float)
        tlist = _np.unique(self.v2t)
        tlist2 = [_np.mean(self.v2B[_np.where(self.v2t == t)]) for t in tlist]
        idx = _np.argsort(tlist2)
        self.v2tavg = _np.array(tlist)[idx]
        self.v2Bavg = _np.array(tlist2)[idx]
        idx = _np.where(self.v2 >= 0)
        self.v2t = self.v2t[idx]
        self.v2B = self.v2B[idx]
        self.v2PA = self.v2PA[idx]
        self.v2lb = self.v2lb[idx]
        self.v2 = self.v2[idx]
        self.v2e = self.v2e[idx]
        return

    def readt3(self, oidata_t3):
        tlist = []
        for t3i in oidata_t3:
            n = len(t3i.wavelength.eff_wave)
            u1 = t3i.u1coord
            v1 = t3i.v1coord
            u2 = t3i.u2coord
            v2 = t3i.v2coord
            uvs = [[u1, v1], [u2, v2]]
            Bm = _np.sqrt((u2+u1)**2 + (v2+v1)**2)
            Pm = _np.arctan2(u2+u1, v2+v1)*180./_np.pi
            for uv in uvs:
                if _np.sqrt(uv[0]**2+uv[1]**2) > Bm:
                    Bm = _np.sqrt(uv[0]**2+uv[1]**2)
                    Pm = _np.arctan2(uv[0], uv[1])*180./_np.pi
            tlist.append([_np.sqrt(u1**2 + v1**2)]*n)
            tlist.append([_np.sqrt(u2**2 + v2**2)]*n)
            tlist.append([_np.arctan2(u1, v1) * 180 / _np.pi]*n)
            tlist.append([_np.arctan2(u2, v2) * 180 / _np.pi]*n)
            tlist.extend([t3i.wavelength.eff_wave])
            tlist.extend([_np.array(t3i.t3phi)])
            tlist.extend([_np.array(t3i.t3phierr)])
            tlist.append([Bm]*n)
            tlist.append([Pm]*n)
        # tlist = _np.array(tlist)
        self.t3B1 = _np.array(_phc.flatten(tlist[0::9]))
        self.t3B2 = _np.array(_phc.flatten(tlist[1::9]))
        self.t3PA1 = _np.array(_phc.flatten(tlist[2::9]))
        self.t3PA2 = _np.array(_phc.flatten(tlist[3::9]))
        self.t3lb = _np.array(_phc.flatten(tlist[4::9]))
        self.t3p = _np.array(_phc.flatten(tlist[5::9]))
        self.t3pe = _np.array(_phc.flatten(tlist[6::9]))
        self.t3Bm = _np.array(_phc.flatten(tlist[7::9]))
        self.t3Pm = _np.array(_phc.flatten(tlist[8::9]))
        return

    def readvis(self, oidata_vis):
        tlist = []
        self.idt = _np.empty(0)
        self.idB = _np.empty(0)
        self.idPA = _np.empty(0)
        for vi in oidata_vis:
            n = len(vi.wavelength.eff_wave)
            if (vi.station[0] and vi.station[1]):
                tlist.append([vi.station[0].sta_name + 
                    vi.station[1].sta_name]*n)
            else:
                tlist.append(['unnamed']*n)
            for st in vi.station:
                if st.sta_name not in self.idt:
                    self.idt = _np.append(self.idt, st.sta_name)
            u = vi.ucoord
            v = vi.vcoord
            B = _np.sqrt(u**2 + v**2)
            PA = _np.arctan2(u, v) * 180 / _np.pi
            if PA < 0:
                PA += 360
            if B not in self.idB:
                self.idB = _np.append(self.idB, B)
                self.idPA = _np.append(self.idPA, PA)
            tlist.append([B]*n)
            tlist.append([PA]*n)
            tlist.extend([vi.wavelength.eff_wave])
            tlist.extend([_np.array(vi.visamp)])
            tlist.extend([_np.array(vi.visamperr)])
            tlist.extend([_np.array(vi.visphi)])
            tlist.extend([_np.array(vi.visphierr)])
        self.v2t = _np.array(_phc.flatten(tlist[0::8]))
        self.vB = _np.array(_phc.flatten(tlist[1::8]))
        self.vPA = _np.array(_phc.flatten(tlist[2::8]))
        self.vlb = _np.array(_phc.flatten(tlist[3::8]))
        self.v = _np.array(_phc.flatten(tlist[4::8]))
        self.ver = _np.array(_phc.flatten(tlist[5::8]))
        self.vp = _np.array(_phc.flatten(tlist[6::8]))
        self.vpe = _np.array(_phc.flatten(tlist[7::8]))
        return

    def readhdrinfo(self, oidata_hdrinfo):
        self.dateobs = oidata_hdrinfo.dateobs
        self.target = oidata_hdrinfo.target
        self.datereduc = oidata_hdrinfo.datereduc
        self.mjd = oidata_hdrinfo.mjd
        return

    def clear(self):
        """ Clear observational info to be used by a model 
        """
        self.oifile = None
        if self.oitype.upper().startswith('PIO'):
            self.v2 = _np.empty(0)
            self.v2e = _np.empty(0)
            self.t3p = _np.empty(0)
            self.t3pe = _np.empty(0)
        elif self.oitype.upper().startswith('AMB'):
            self.v = _np.empty(0)
            self.ve = _np.empty(0)
            self.vp = _np.empty(0)
            self.vpe = _np.empty(0)
        return

    def add_model(self, imgs, rpx, PAdisk=90, modid=None, lbm=None):
        """ Fill the class with modeled values from observational setup already
        read.

        ``imgs`` is and matrix(lb, n, m), where lb is the lambda dimension (
        optional). ``rpx`` is radians_per_pixel.
        """
        self.clear()
        self.modid = modid
        if len(_np.shape(imgs)) == 3:
            iscube = True
        else:
            iscube = False
            img = imgs
        if self.oitype.upper().startswith('PIO'):
            if lbm is not None:
                self.v2lbm = _np.array(lbm)
                self.t3lbm = _np.array(lbm)
            # else:
            #     self.v2lbm = _np.array([1.66e-6])
            #     self.t3lbm = _np.array([1.66e-6])
            #
            tlist = []
            for i in range(len(self.v2B)):
                if iscube:
                    img = imgs[_phc.find_nearest(self.v2lbm, self.v2lb[i], 
                        idx=True)]
                _, v, _ = fastnumvis(img, self.v2lb[i], self.v2B[i], 
                    self.v2PA[i], rpx, PAdisk=PAdisk)
                tlist.extend([v**2])
            self.v2 = _np.array(tlist)
            self.v2e = _np.zeros(_np.shape(tlist))*_np.NaN
            #
            tlist = []
            for i in range(len(self.t3B1)):
                if iscube:
                    img = imgs[_phc.find_nearest(self.t3lbm, self.t3lb[i], 
                        idx=True)]
                _, _, p = fastnumvis3(img, self.t3lb[i], [self.t3B1[i], 
                    self.t3B2[i]], [self.t3PA1[i], self.t3PA2[i]], rpx,
                    PAdisk=PAdisk)
                tlist.extend([p])
            self.t3p = _np.array(tlist)
            self.t3pe = _np.zeros(_np.shape(tlist))*_np.NaN
        #
        elif self.oitype.upper().startswith('AMB'):
            if lbm is not None:
                self.vlbm = _np.array(lbm) 
            # else:
            #     self.vlbm = _np.array([2.21e-6])               
            self.v = _np.empty(0)
            self.vp = _np.empty(0)
        return


def readoifits(oifile, oitype='PIO', quiet=True):
    """ Read the oifits files/oidata according to its oitype.

    Possible oitypes=["PIO", "AMB"]
    """
    if isinstance(oifile, _strtypes):
        oidata = _oifits.open(oifile, quiet=quiet)
    else:
        oidata = _copy(oifile)
        oifile = None
    if oitype.upper().startswith('PIO'):
        oic = oiClass(oitype, oifile)
        oic.readvis2(oidata.vis2)
        oic.readt3(oidata.t3)
    elif oitype.upper().startswith('AMB'):
        oic = oiClass(oitype, oifile)
        oic.readvis(oidata.vis)
        oic.readhdrinfo(oidata.hdrinfo)
    else:
        raise NameError('Invalid `oitype` in intt.readoifits')
    return oic


def plot_pio_res(loic, oitype='PIO', outname=None, legend=True, fmt=['png'], 
    PAs=[-180, 180], lcor=None, lmk=None, legbbox=None, ard=None, rylim=None,
    v2leg=None):
    """ Plot a OIFITS files with VIS2 and CP information and RESIDUALS.

    :param loic: list of ``oiC``. Only the first ``loi`` will have the *uv* 
    plane coverage plotted.
    :param legbbox: legend(bbox_to_anchor) coordinates (list or tuple).

    If ``PAs != [-180, 180]``, then a limited range of PAs will be plotted.
    """
    if isinstance(loic, oiClass):
        loic = [loic]
    if lcor is None:
        lcor = [_phc.cycles(i) for i in range(len(loic))]
    if lmk is None:
        lmk = [_phc.cycles(i, 'mk') for i in range(len(loic))]
    fig = _plt.figure(figsize=1.5*_np.array([9/0.871875, 5]))
    lins, cols = (1, 8)
    gs0 = _gridspec.GridSpec(lins, cols)  # , width_ratios=[1]*cols, 
        # height_ratios=[1]*lins)
    gs0.update(wspace=4.)
    gs00 = _gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs0[0, 
        0:cols-2], wspace=0.1, hspace=0.01)
    gs01 = _gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0, 
        cols-2:], hspace=1.)
    axB1 = _plt.subplot(gs01[0, 0]) 
    axB2 = _plt.subplot(gs01[1, 0]) 
    axV = _plt.subplot(gs00[0, 0]) 
    axVr = _plt.subplot(gs00[1:, 0]) 
    axC = _plt.subplot(gs00[0, 1]) 
    axCr = _plt.subplot(gs00[1:, 1]) 
    print(get_ax_size(axB1, fig))
    for i in range(len(loic)):
        oic = loic[i]
        if i == 0:
            axB1 = uv_plot_ax(axB1, oic, PAs=PAs, legend=legend, lbd=False, 
                legbbox=legbbox, ard=ard)
            axB2 = uv_plot_ax(axB2, oic, PAs=PAs, legend=False, lbd=True, 
                ard=ard)
        axV = v2_plot_ax(axV, oic, PAs=PAs, cor=lcor[i], mk=lmk[i])
        axC = t3_plot_ax(axC, oic, PAs=PAs, cor=lcor[i], mk=lmk[i])
    axVr = v2res_plot_ax(axVr, loic, PAs=PAs, lcor=lcor, lmk=lmk, rylim=rylim)
    axCr = t3res_plot_ax(axCr, loic, PAs=PAs, lcor=lcor, lmk=lmk, rylim=rylim)
    if v2leg is not None:
        lines = []
        for i in range(len(loic)):
            lines.extend([_Line2D([], [], color=lcor[i], marker=lmk[i], 
                ls='', label=v2leg[i])])
        fs = dict(_mpl.rcParams.viewitems())['font.size']
        axV.legend(handles=lines, loc=1, fancybox=True, framealpha=0.5, 
            fontsize=fs-4, labelspacing=0.05)
    axC.yaxis.set_label_position("right")
    axC.yaxis.tick_right()
    axC.yaxis.set_ticks_position('both')
    axCr.yaxis.set_label_position("right")
    axCr.yaxis.tick_right()
    axCr.yaxis.set_ticks_position('both')
    axB1.yaxis.set_label_position("right")
    axB1.yaxis.tick_right()
    axB1.yaxis.set_ticks_position('both')
    axB2.yaxis.set_label_position("right")
    axB2.yaxis.tick_right()
    axB2.yaxis.set_ticks_position('both')
    axV.set_ylabel('V$^2$')
    axC.set_ylabel(r'$\phi_{123}$ (deg$.$)')
    axCr.set_ylabel(r'$\phi_{123}$(data-mod)/err')
    axVr.set_ylabel(r'$\phi_{123}$(data-mod)/err')
    axB1.set_xlabel(r'$B_{\rm proj}$ (m)')
    axB2.set_xlabel(r'$B_{\rm proj}/\lambda$')
    axCr.set_xlabel(r'$B_{\rm proj}/\lambda$')
    axVr.set_xlabel(r'$B_{\rm proj}/\lambda$')
    axV.set_xticklabels([])
    axC.set_xticklabels([])
    axVr.xaxis.set_major_locator(_mtick.MaxNLocator(nbins=5, prune='upper'))
    axCr.xaxis.set_major_locator(_mtick.MaxNLocator(nbins=5, prune='lower'))
    axVr.yaxis.set_major_locator(_mtick.MaxNLocator(nbins=8, prune='upper'))
    axCr.yaxis.set_major_locator(_mtick.MaxNLocator(nbins=8, prune='upper'))    
    axV.xaxis.set_major_locator(_mtick.MaxNLocator(nbins=5, prune='lower'))
    axC.xaxis.set_major_locator(_mtick.MaxNLocator(nbins=5, prune='lower'))
    axB1.xaxis.set_major_locator(_mtick.MaxNLocator(nbins=4))
    axB2.xaxis.set_major_locator(_mtick.MaxNLocator(nbins=4))
    axV.minorticks_on()
    axVr.minorticks_on()
    axC.minorticks_on()
    axCr.minorticks_on()
    if outname is None:
        outname = _phc.dtflag()
    _phc.savefig(fig, figname=outname, fmt=fmt)
    return


def get_ax_size(ax, fig):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height    


def plot_pionier(loic, oitype='PIO', outname=None, legend=True, fmt=['png'],
    PAs=[-180, 180], lcor=None, lmk=None):
    """ Plot a OIFITS files with VIS2 and CP information. 

    ``loic`` is a list of oic classes.

    Only the first ``loi`` will have the *uv* plane coverage plotted.

    If ``PAs != [-180, 180]``, then a limited range of PAs will be plotted.
    """
    if isinstance(loic, oiClass):
        loic = [loic]
    if lcor is None:
        lcor = [_phc.cycles(i) for i in range(len(loic))]
    if lmk is None:
        lmk = [_phc.cycles(i, 'mk') for i in range(len(loic))]
    lins, cols = (2, 5)
    fig = _plt.figure(figsize=3.*_np.array([cols, lins]))
    gs = _gridspec.GridSpec(lins, cols, width_ratios=[1]*lins, 
        height_ratios=[1]*cols)
    axB1 = _plt.subplot(gs[0, 4]) 
    axB2 = _plt.subplot(gs[1, 4]) 
    axV = _plt.subplot(gs[:, 0:2]) 
    axC = _plt.subplot(gs[:, 2:4])
    for i in range(len(loic)):
        oic = loic[i]
        if i == 0:
            axB1 = uv_plot_ax(axB1, oic, PAs=PAs, legend=legend, lbd=False)
            axB2 = uv_plot_ax(axB2, oic, PAs=PAs, legend=False, lbd=True)
        axV = v2_plot_ax(axV, oic, PAs=PAs, cor=lcor[i], mk=lmk[i])
        axC = t3_plot_ax(axC, oic, PAs=PAs, cor=lcor[i], mk=lmk[i])
    if outname is None:
        outname = _phc.dtflag()
    _phc.savefig(fig, figname=outname, fmt=fmt)
    return


def t3res_plot_ax(ax, loic, PAs=[-180, 180], lcor=None, lmk=None, rylim=None):
    """ Axis plot of t3 residuals

    :param loic: list of ``oiC`` models. The first one is considered as 
    observational data, and the others as models.
    """
    if lcor is None:
        lcor = [_phc.cycles(i) for i in range(len(loic))]
    if lmk is None:
        lmk = [_phc.cycles(i, 'mk') for i in range(len(loic))]
    PArv = reversePAs(PAs)
    #
    idx = _np.where( ((loic[0].t3Pm >= PAs[0]) & (loic[0].t3Pm <= PAs[1])) | (
        (loic[0].t3Pm >= PArv[0]) & (loic[0].t3Pm <= PArv[1])) )
    for j, oic in enumerate(loic[1:], 1):
        ax.scatter(loic[0].t3Bm/oic.t3lb[idx], (loic[0].t3p[idx] - 
            oic.t3p[idx]) / (loic[0].t3pe[idx]), marker=lmk[j], alpha=.7, 
            color=lcor[j])
    xlim = ax.get_xlim()
    ax.plot(xlim, [0, 0], color='gray', ls="-", zorder=0)
    ax.plot(xlim, [3, 3], color='k', ls=":", zorder=0)
    ax.plot(xlim, [-3, -3], color='k', ls=":", zorder=0)
    ax.set_xlim(xlim)
    if rylim is not None:
        ax.set_ylim(rylim)
    return ax


def v2res_plot_ax(ax, loic, PAs=[-180, 180], lcor=None, lmk=None, rylim=None):
    """ Axis plot of t3 residuals

    :param loic: list of ``oiC`` models. The first one is considered as 
    observational data, and the others as models.
    """
    if lcor is None:
        lcor = [_phc.cycles(i) for i in range(len(loic))]
    if lmk is None:
        lmk = [_phc.cycles(i, 'mk') for i in range(len(loic))]
    PArv = reversePAs(PAs)
    #
    idx = _np.where( ((loic[0].v2PA >= PAs[0]) & (loic[0].v2PA <= PAs[1])) | (
        (loic[0].v2PA >= PArv[0]) & (loic[0].v2PA <= PArv[1])) )
    for j, oic in enumerate(loic[1:], 1):
        ax.scatter(loic[0].v2B/oic.v2lb[idx], (loic[0].v2[idx] - oic.v2[idx]) /
            (loic[0].v2e[idx]), marker=lmk[j], alpha=.7, color=lcor[j])
    xlim = ax.get_xlim()
    ax.plot(xlim, [0, 0], color='gray', ls="-", zorder=0)
    ax.plot(xlim, [3, 3], color='k', ls=":", zorder=0)
    ax.plot(xlim, [-3, -3], color='k', ls=":", zorder=0)
    ax.set_xlim(xlim)
    if rylim is not None:
        ax.set_ylim(rylim)
    return ax


def t3_plot_ax(ax, oic, PAs=[-180, 180], cor='k', mk='o'):
    """ Axis plot of t3
    """
    PArv = reversePAs(PAs)
    #
    idx = _np.where( ((oic.t3Pm >= PAs[0]) & (oic.t3Pm <= PAs[1])) | (
        (oic.t3Pm >= PArv[0]) & (oic.t3Pm <= PArv[1])) )
    ax.errorbar(oic.t3Bm/oic.t3lb[idx], oic.t3p[idx], 
        yerr=oic.t3pe[idx], ls='', marker=mk, alpha=.7, color=cor,
        markeredgecolor=cor)
    xlim = ax.get_xlim()
    ax.plot(xlim, [0, 0], color='gray', ls="-", zorder=0)
    ax.set_xlim(xlim)
    return ax


def v2_plot_ax(ax, oic, PAs=[-180, 180], cor='k', mk='o'):
    """ Axis plot of v2 
    """
    PArv = reversePAs(PAs)
    #
    idx = _np.where( ((oic.v2PA >= PAs[0]) & (oic.v2PA <= PAs[1])) | (
        (oic.v2PA >= PArv[0]) & (oic.v2PA <= PArv[1])) )
    ax.errorbar(oic.v2B/oic.v2lb[idx], oic.v2[idx], 
        yerr=oic.v2e[idx], ls='', marker=mk, alpha=.7, color=cor,
        markeredgecolor=cor)
    return ax


def uv_plot_ax(ax, oic, PAs=[-180, 180], legend=True, lbd=False, legbbox=None,
    ard=None):
    """ Axis plot of the observational baselines (u, v) plane. 

    :param legbbox: legend(bbox_to_anchor) coordinates (list or tuple).
    :param ard: Is not ``None``, arrow PA (in deg).
    """
    PArv = reversePAs(PAs)
    #
    if oic.oitype.upper().startswith('PIO'):
        B = oic.v2B
        PA = oic.v2PA
        # t2Bavg = dict(zip(oic.v2tavg, oic.v2Bavg))
        # tlist = [t2Bavg[t] for t in oic.v2t]
        # colors = _phc.gradColor(tlist, min=oic.v2Bavg[0] - (
        #     oic.v2Bavg[-1]-oic.v2Bavg[0])*0.1, cmapn='viridis')
        t2Bavg = dict(zip(oic.v2tavg, range(len(oic.v2tavg))))
        tlist = [t2Bavg[t] for t in oic.v2t]
        colors = _phc.gradColor(tlist, cmapn='copper_r')
        lbds = oic.v2lb
    elif oic.oitype.upper().startswith('AMB'):
        B = oic.vB
        PA = oic.vPA
        lbds = oic.vlb
        colors = _phc.gradColor(oic.vB)
    if not lbd:
        lbds = _np.ones(len(B))
    #
    idx = _np.where( ((PA >= PAs[0]) & (PA <= PAs[1])) | (
        (PA >= PArv[0]) & (PA <= PArv[1])) )
    #
    u = B*_np.sin(PA*_np.pi/180.)/lbds
    v = B*_np.cos(PA*_np.pi/180.)/lbds
    ax.scatter(u[idx], v[idx], color=colors[idx], alpha=.7, zorder=10)
    ax.scatter(-u[idx], -v[idx], color=colors[idx], alpha=.7, zorder=10)
    if legend:
        # colors = _phc.gradColor(oic.v2Bavg, min=oic.v2Bavg[0] - (
        #     oic.v2Bavg[-1]-oic.v2Bavg[0])*0.1, cmapn='viridis')
        colors = _phc.gradColor(range(len(oic.v2Bavg)), cmapn='copper_r')
        lines = []
        for i in range(len(oic.v2Bavg)):
            lines.extend([_Line2D(range(1), range(1), color=colors[i], 
                marker='o', ls='')])
        fs = dict(_mpl.rcParams.viewitems())['font.size']-4
        if legbbox is None:
            ax.legend(lines, oic.v2tavg, numpoints=1, loc='best', 
                fancybox=True, framealpha=0.5, fontsize=fs, labelspacing=0.05)
        else:
            ax.legend(lines, oic.v2tavg, numpoints=1, loc='best', 
                fancybox=True, framealpha=0.5, fontsize=fs, labelspacing=0.05, 
                bbox_to_anchor=legbbox)
    if ard is not None:
        Bx = _np.max(B)*_np.sin(ard*_np.pi/180)/_np.max(lbds)
        By = _np.max(B)*_np.cos(ard*_np.pi/180)/_np.max(lbds)
        ax.arrow(Bx, -By, -2*Bx, 2*By, length_includes_head=True,
        head_width=0.1*_np.max(B)/_np.max(lbds), fc='gray', ec='gray',
        head_length=0.15*_np.max(B)/_np.max(lbds), zorder=1)
        ax.plot([-By, By], [-Bx, Bx], color='gray', ls='--')
    ax.minorticks_on()
    # ax.xaxis.tick_top()
    ax.axis('equal')
    # if lbd:
    #     ax.yaxis.set_major_formatter(_mtick.FormatStrFormatter('%.0e'))
    #     ax.xaxis.set_major_formatter(_mtick.FormatStrFormatter('%.0e'))
    # ax.locator_params(nbins=5)
    return ax


def plot_baseline(ax, oidata, colors, names, PAs=[-180, 180], PAsrev=False, 
    PArv=[0, 0], xlim=None, legend=True):
    leg = []
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        PA = _np.arctan2(u, v) * 180 / _np.pi
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[names.index(label)]
        # [u,-u] = W > E
        # [-u,u] = E < W
        if (PA >= PAs[0] and PA <= PAs[1]) or (PAsrev and 
            (PA >= PArv[0]) and PA <= PArv[1]):
            if label not in leg:
                ax.plot([u, -u], [v, -v], '.', label=label, color=color)
                leg.append(label)
            else:
                ax.plot([u, -u], [v, -v], '.', color=color)  # label=label,
        # names.append(vis2.target.target)
    f = 1.
    x = 84.06*f
    y = 111.96*f
    ax.arrow(-x, y, 2*x, -2*y, length_includes_head=True,
        head_width=0.1*111.96*2*f, 
        head_length=0.15*111.96*2*f, fc='gray', ec='gray', zorder=1)
    ax.plot([-y, y], [-x, x], color='gray', ls='--')
    # ax.plot([0], [0], 'o')
    # names = list(_np.unique(names))
    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel(u'B$_{proj}$ (m)')
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.axis('equal')
    # _plt.grid(b=True, linestyle=':', alpha=alp)
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        idx = [labels.index(i) for i in list(names)]
        handles = _np.array(handles)[idx]
        ax.legend(handles, names, prop={'size': 8}, numpoints=1, loc='best', 
            framealpha=0.5, fancybox=True, labelspacing=0.1)
    # ax.yaxis.set_major_formatter(_mtick.FormatStrFormatter('%.0e'))
    # ax.xaxis.set_major_formatter(_mtick.FormatStrFormatter('%.0e'))
    # ax.ticklabel_format(useOffset=False)  # style='sci' DO NOT WORK
    ax.locator_params(nbins=7)
    return ax


# Other?
def mapinterf(modf, im=0, obs=0, iflx=0, dist=10, PA=0., B=100., PAdisk=90.,
    quiet=False):
    """ Return Squared Visibilities (V2) and Diferential Phases (DP) for a given
    `hdust` map(s) file.

    If *.map file format, it takes `iflx` image layer.

    input: *.map(s) path (string), `dist` (float, parsecs)

    output: lbdc, V2, DP (float arrays) """
    data, obslist, lbdc, Ra, xmax = readmap(modf, quiet=quiet)
    pixsize = 2 * xmax[0] / _np.shape(data)[-1]
    rad_per_pixel = _np.double(pixsize * 6.96E10 / (dist * 3.08567758E18))
    npts = len(lbdc)
    V2 = _np.zeros(npts)
    DP = _np.zeros(npts)
    for i in range(npts):
        if len(_np.shape(data)) == 5:
            img = data[im, obs, i, ::-1, :]
        else:
            img = data[im, obs, i, :, :, iflx]
        tmp, V, DP[i] = fastnumvis(img, lbdc[i] * 1e-6, B, PA, rad_per_pixel,
            PAdisk=PAdisk)
        V2[i] = V**2
    # avg = (DP[0]+DP[-1])/2.
    ssize = int(.05 * len(DP))
    if ssize == 0:
        ssize = 1
    medx0, medx1 = _np.average(DP[:ssize]), _np.average(DP[-ssize:])
    avg = (medx0 + medx1) / 2.
    DP = DP - avg
    # DP = _spt.linfit(lbdc, DP)-1.
    # lbdc = _spt.air2vac(lbdc*1e4)*1e-4
    return lbdc, V2, DP


def datinterf(data, lbdc, xmax, im=0, obs=0, iflx=0, dist=10, PA=0., B=100., 
    PAdisk=90., quiet=False, normV2=False):
    """ Return Squared Visibilities (V2) and Diferential Phases (DP) for a given
    `hdust` data file.

    If *.map file format, it takes `iflx` image layer.

    input: data (np.ndarray), `dist` (float, parsecs)

    output: lbdc, V2, DP (float arrays) """
    pixsize = 2 * xmax[im] / _np.shape(data)[-1]
    rad_per_pixel = _np.double(pixsize * _phc.Rsun.cgs / (dist * _phc.pc.cgs))
    npts = len(lbdc)
    V2 = _np.zeros(npts)
    DP = _np.zeros(npts)
    for i in range(npts):
        if len(_np.shape(data)) == 5:
            img = data[im, obs, i, ::-1, :]
        else:
            img = data[im, obs, i, :, :, iflx]
        tmp, V, DP[i] = fastnumvis(img, lbdc[i] * 1e-6, B, PA, rad_per_pixel,
            PAdisk=PAdisk)
        V2[i] = V**2
    # avg = (DP[0]+DP[-1])/2.
    ssize = int(.05 * len(DP))
    if ssize == 0:
        ssize = 1
    medx0, medx1 = _np.average(DP[:ssize]), _np.average(DP[-ssize:])
    avg = (medx0 + medx1) / 2.
    DP = DP - avg
    # DP = _spt.linfit(lbdc, DP)-1.
    # lbdc = _spt.air2vac(lbdc*1e4)*1e-4
    if normV2:
        V2 = _linfit(lbdc, V2, ssize=0.1)
    return lbdc, V2, DP


def dataphot(data, lbdc, xmax, im=0, obs=0, iflx=0, dist=10, PA=0., B=100., 
    PAdisk=90., quiet=False, diff=True, phiunit=True): 
    """ Return photocenter positions for a given `hdust` data file.

    If *.map file format, it takes `iflx` image layer.

    input: data (np.ndarray), `dist` (float, parsecs)

    output: lbdc, V2, DP (float arrays) """
    pixsize = 2 * xmax[im] / _np.shape(data)[-1]
    rad_per_pixel = _np.double(pixsize * _phc.Rsun.cgs / (dist * _phc.pc.cgs))
    npts = len(lbdc)
    px = _np.zeros(npts)
    py = _np.zeros(npts)
    for i in range(npts):
        if len(_np.shape(data)) == 5:
            img = data[im, obs, i, ::-1, :]
        else:
            img = data[im, obs, i, :, :, iflx]
        img = _img.rotate_image(img, -(PA+90-PAdisk)*_np.pi/180)
        x, y = _np.shape(img)
        x = +_np.linspace(-(x-1)/2., (x-1)/2., x)
        y = +_np.linspace(-(y-1)/2., (y-1)/2., y)
        px[i] = _np.sum(_np.sum(img, axis=0)*x)/_np.sum(img)
        py[i] = _np.sum(_np.sum(img, axis=1)*y)/_np.sum(img)
    if diff:
        ssize = int(.05 * len(px))
        if ssize == 0:
            ssize = 1
        avg = _np.mean([_np.mean(px[:ssize]), _np.mean(px[-ssize:])])
        px -= avg
        ssize = int(.05 * len(py))
        if ssize == 0:
            ssize = 1
        avg = _np.mean([_np.mean(py[:ssize]), _np.mean(py[-ssize:])])
        py -= avg
    if phiunit:
        px *= -2*_np.pi*B/lbdc*1e6*rad_per_pixel*180/_np.pi
        py *= -2*_np.pi*B/lbdc*1e6*rad_per_pixel*180/_np.pi
    return lbdc, px, py


def plot_v2_img(ax, oidata, img, xlim=None, quiet=False, alp=.75):
    """ IMG from LitPRO

    TODO: scales
    """
    ms = 3
    ax.plot([0, 1e8], [0, 0], ls='--', zorder=0, color='k', alpha=alp)
    # Loop in vis2 creating the following vectors
    x, y, yerr, Blist, PAlist, lbd = ([], [], [], [], [], [])
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = _phc.cycles(list(oidata.vis2).index(vis2))
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u, v) * 180 / _np.pi
        # if (PA >= PAs[0] and PA <= PAs[1]) or (PAsrev and (PA >= PArv[0]) and 
        #     PA <= PArv[1]):
            # if bindata <= 0:
        lbd.extend(vis2.wavelength.eff_wave)
        x.extend(B/vis2.wavelength.eff_wave)
        y.extend(vis2.vis2data)
        yerr.extend(vis2.vis2err)
        Blist.extend([B]*len(vis2.wavelength.eff_wave))
        PAlist.extend([PA]*len(vis2.wavelength.eff_wave))
        ax.errorbar(B / vis2.wavelength.eff_wave, 
            vis2.vis2data, yerr=vis2.vis2err, fmt='o', markersize=ms, 
            color=color, alpha=alp)
    # # 
    fits = _pyfits.open(img)
    dx = fits[0].header['CDELT1']
    dy = fits[0].header['CDELT2']
    if dx != dy:
        print('# ERROR! Differential spatial scales for {0}'.format(img))
        return
    if fits[0].header['CUNIT1'].lower().find('deg') > -1:
        rad_per_pixel = dx*_phc.deg2rad()
    else:
        rad_per_pixel = dx
    data = fits[0].data
    data = data[:, :]
    V2 = []
    for i in range(len(x)):
        tmp, V, ph = fastnumvis(data, lbd[i], Blist[i], PAlist[i], 
            rad_per_pixel, PAdisk=90.)
        V2.append(V**2)
    ax.plot(x, V2, color='k', ls="", markersize=ms*2, marker="d")
    return ax, V2


def plot_t3_img(ax, oidata, img, xlim=None, quiet=False, alp=.75):
    """ IMG from LitPRO
    
    TODO: scales
    """
    ms = 3
    ax.plot([0, 1e8], [0, 0], ls='--', zorder=0, color='k', alpha=alp)
    # Loop in vis2 creating the following vectors
    x, y, yerr, Blist, PAlist, Bmax = ([], [], [], [], [], [])
    for t3 in oidata.t3:
        u1 = t3.u1coord
        v1 = t3.v1coord
        u2 = t3.u2coord
        v2 = t3.v2coord
        if _np.sqrt(u1**2 + v1**2) > _np.sqrt(u2**2 + v2**2):
            u = u1
            v = v1
        else:
            u = u2
            v = v2
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u, v) * 180.0 / _np.pi
        B = _np.tile(B, len(t3.wavelength.eff_wave))
        PA = _np.tile(PA, len(t3.wavelength.eff_wave))
        Bmax.extend(B)
        lbsz = len(t3.wavelength.eff_wave)
        x.extend(B/t3.wavelength.eff_wave)
        y.extend(t3.t3phi)
        yerr.extend(t3.t3phierr)
        B1, B2 = ( _np.sqrt(u1**2 + v1**2), _np.sqrt(u2**2 + v2**2) )
        Blist.extend([[B1, B2] for i in range(lbsz)])
        PA1, PA2 = ( _np.arctan2(u1, v1) * 180.0 / _np.pi, 
            _np.arctan2(u2, v2) * 180.0 / _np.pi )
        PAlist.extend([[PA1, PA2] for i in range(lbsz)])
        ax.errorbar(x, y, yerr=yerr, color='k', fmt='o', markersize=ms, 
            alpha=alp)
    # # 
    fits = _pyfits.open(img)
    dx = fits[0].header['CDELT1']
    dy = fits[0].header['CDELT2']
    if dx != dy:
        print('# ERROR! Differential spatial scales for {0}'.format(img))
        return
    if fits[0].header['CUNIT1'].lower().find('deg') > -1:
        rad_per_pixel = dx*_phc.deg2rad()
    else:
        rad_per_pixel = dx
    data = fits[0].data
    lblist = _np.array(Bmax) / _np.array(x) 
    t3 = []
    for i in range(len(x)):
        tmp, V, ph = fastnumvis3(data, lblist[i], Blist[i], PAlist[i], 
            rad_per_pixel, PAdisk=90.)
        t3.append(ph)
    ax.plot(x, t3, color='k', ls="", markersize=ms*2, marker="d")
    return ax, t3


def plot_v2_res(ax, ax2, oidata, colors, names, modfiles, obsdeg=None, 
    dist=None, xlim=None, PAs=[-180, 180], PAsrev=False, PArv=[0, 0], 
    bindata=0, quiet=False, alp=.75, printsum=False):
    r""" datas = models! 

    `bindata` refers to vis2.wavelength.eff_wave!! In other words: the binning
    only works on simultaneous observations with different :math:`\lambda`.
    """
    ms = 3.5
    ax2.plot([0, 1e8], [0, 0], ls='--', zorder=0, color='k', alpha=alp)
    # Loop in vis2 creating the following vectors
    x, y, yerr, Blist, PAlist = ([], [], [], [], [])
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[names.index(label)]
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u, v) * 180 / _np.pi
        if (PA >= PAs[0] and PA <= PAs[1]) or (PAsrev and (PA >= PArv[0]) and 
            PA <= PArv[1]):
            if bindata <= 0:
                x.extend(B / vis2.wavelength.eff_wave)
                y.extend(vis2.vis2data)
                yerr.extend(vis2.vis2err)
                Blist.extend([B] * len(vis2.wavelength.eff_wave))
                PAlist.extend([PA] * len(vis2.wavelength.eff_wave))
                ax.errorbar(B / vis2.wavelength.eff_wave, 
                    vis2.vis2data, yerr=vis2.vis2err, fmt='o', markersize=ms, 
                    color=color, alpha=alp, zorder=1)
            # Individual plot to keep each BASELINE color
            if bindata > 0:
                binned = _phc.bindata(B / vis2.wavelength.eff_wave, 
                    vis2.vis2data, bindata, yerr=vis2.vis2err)  # xlim=xlim)
                ax.errorbar(binned[0], binned[1], yerr=binned[2], 
                    fmt='o', markersize=ms, color=color, zorder=1)
    # Plotting models
    if obsdeg is None:
        obs = [0]
    else:
        obs = obsdeg
    mcolors = _phc.gradColor(_np.arange(len(modfiles) + 2), cmapn='Greens_r')
    mcolors = _np.array(mcolors)[1:-1]
    mcolors = ['dimgray', 'blue']
    mcolors = ['green', 'red']
    markers = ['*', 'd']
    # mcolors = ['red', 'green']
    for mod in modfiles:
        midx = modfiles.index(mod)
        data, obslist, lbdc, Ra, xmax = readmap(mod, quiet=quiet)
        data = data[..., ::-1, :, :]
        # It is just taking the length. Ignore the indexes
        pixsize = 2 * xmax[0] / len(data[0, 0, 0, :, :, 0])
        # *60.*60.*1000.*180./_np.pi)
        rad_per_pixel = _np.double(pixsize * _phc.Rsun.cgs / (dist * 
            _phc.pc.cgs))
        lblist = _np.array(Blist) / _np.array(x) 
        for deg in obs:
            V2 = []
            ob = list(obslist[::2]).index(_phc.find_nearest(obslist[::2], deg))
            if _np.abs(obslist[::2][ob] - deg) > 1:
                print('# The difference of angles is bigger than 1 deg!!!')
                print('# Input: {0:.1f}; Nearest model: {1:.1f}'.format(
                    deg, obslist[::2][ob]))
            for i in range(len(x)):
                j = list(lbdc * 1e-6).index(_phc.find_nearest(lbdc * 1e-6, 
                    lblist[i]))
                if _np.abs(lbdc[j]*1e-6 - lblist[i]) > 0.1*1e-6:
                    print('# The difference of lambd is bigger than 0.1 um!!!')
                    print('# Input: {0:.1f}; Nearest model: {1:.1f}'.format(
                        lblist[i], lbdc[j]*1e-6))
                # lbcalc = lbdc[j]*1e-6
                tmp, V, phvar = fastnumvis(data[0, ob, j, :, :, 0], lblist[i], 
                    Blist[i], PAlist[i], rad_per_pixel, PAdisk=(216.9 + 90.))
                V2 += [V**2]
            fig2 = _plt.figure()
            imshowl(data[0, ob, j, :, :, 0])
            _phc.savefig(fig2, figname='{0}_{1}'.format(_os.path.basename(mod), 
                ob))
            V2 = _np.array(V2) 
            if bindata <= 0:
                ax.plot(x, V2, color=mcolors[midx], alpha=alp, ls='',
                    markersize=ms,
                    marker=markers[midx])
                ax2.plot(x, (y - V2) / yerr, ls='', color=mcolors[midx],
                    markersize=ms, marker=markers[midx], 
                    zorder=1, alpha=alp)
                if printsum:
                    chi2 = _phc.chi2calc(V2, y, yerr)
                    ax2.text(0.4, 0.87, r'$\sum={0:.1f}$'.format(chi2), 
                        transform=ax2.transAxes)
            if printsum:
                v2sum = _np.sum((y - V2) / yerr)
                print(v2sum)
                ax2.text(0.4, 0.87, r'$\sum={0:.1f}$'.format(v2sum), 
                    transform=ax2.transAxes)
    # Blim = ax3.get_xlim()
    # ax4.set_xlim(Blim)
    # ax1.xaxis.get_major_formatter().set_useOffset(False)
    ax.xaxis.set_major_formatter(_mtick.FormatStrFormatter('%.2e'))
    # ax3.get_xaxis().set_visible(False)
    # ax4.get_xaxis().set_visible(False)
    ax.get_xaxis().set_ticklabels([])
    ax.set_ylim([0, 1.1])
    ax.set_ylabel(u'$V$ $^2$')
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.grid(b=True, linestyle=':', alpha=alp)
    ax2.set_ylim([-6, 6])
    ax2.get_xaxis().set_ticklabels([])
    ax2.set_ylabel(u'$V$ $^2$(data-mod)/err')
    if xlim is not None:
        ax2.set_xlim(xlim)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.grid(b=True, linestyle=':', alpha=alp)
    return ax, ax2


def plot_phi3_res(ax, ax2, oidata, colors, names, modfiles, obsdeg=None, 
    dist=None, xlim=None, PAs=[-180, 180], PAsrev=False, PArv=[0, 0], 
    bindata=0, quiet=False, alp=.75, printsum=False):
    r""" modfiles = models! 

    `bindata` refers to vis2.wavelength.eff_wave!! In other words: the binning
    only works on simultaneous observations with different :math:`\lambda`.
    """
    ms = 3.5
    ax.plot([0, 1e8], [0, 0], ls='--', zorder=0, color='k', alpha=alp)
    ax2.plot([0, 1e8], [0, 0], ls='--', zorder=0, color='k', alpha=alp)
    color = 'black'
    x, y, yerr, Blist, PAlist, Bmax = ([], [], [], [], [], [])
    for t3 in oidata.t3:
        u1 = t3.u1coord
        v1 = t3.v1coord
        u2 = t3.u2coord
        v2 = t3.v2coord
        # if _np.sqrt(u1**2 + v1**2) > _np.sqrt(u2**2 + v2**2):
        #     u = u1
        #     v = v1
        # else:
        #     u = u2
        #     v = v2
        # B = _np.sqrt(u**2 + v**2)
        # PA = _np.arctan2(u, v) * 180.0 / _np.pi
        Bcompar = _np.sqrt( [u1**2+v1**2, u2**2+v2**2, (u1+u2)**2+(v1+v2)**2])
        PAcompar = _np.array( [_np.arctan2(u1, v1), _np.arctan2(u2, v2), 
            _np.arctan2(u1+u2, v1+v2)] )*180/_np.pi
        idx = _np.where(Bcompar == _np.max(Bcompar))
        B = _np.tile(Bcompar[idx], len(t3.wavelength.eff_wave))
        PA = _np.tile(PAcompar[idx], len(t3.wavelength.eff_wave))
        if (PA[0] >= PAs[0] and PA[0] <= PAs[1]) or (PAsrev and PA[0] >= 
            PArv[0] and PA[0] <= PArv[1]):
                Bmax.extend(B)
                lbsz = len(t3.wavelength.eff_wave)
                x.extend(B/t3.wavelength.eff_wave)
                y.extend(t3.t3phi)
                yerr.extend(t3.t3phierr)
                B1, B2 = ( _np.sqrt(u1**2 + v1**2), _np.sqrt(u2**2 + v2**2) )
                Blist.extend([[B1, B2] for i in range(lbsz)])
                PA1, PA2 = ( _np.arctan2(u1, v1) * 180.0 / _np.pi, 
                    _np.arctan2(u2, v2) * 180.0 / _np.pi )
                PAlist.extend([[PA1, PA2] for i in range(lbsz)])
    if bindata <= 0:
        ax.errorbar(x, y, yerr=yerr, color=color, fmt='o', markersize=ms,
            zorder=1, alpha=alp)
    else:
        binned = _phc.bindata(B / t3.wavelength.eff_wave, t3.t3phi, 
            bindata, yerr=t3.t3phierr, xlim=xlim)
        ax.errorbar(binned[0], binned[1], yerr=binned[2], 
            color=color, fmt='o', markersize=ms, zorder=1, alpha=alp)
    mcolors = _phc.gradColor(_np.arange(len(modfiles)+2), cmapn='Greens_r')
    mcolors = _np.array(mcolors)[1:-1]
    mcolors = ['dimgray', 'blue']
    mcolors = ['green', 'red']
    markers = ['*', 'd']
    # mcolors = ['red', 'green']
    if obsdeg is None:
        obs = [0]
    else:
        obs = obsdeg
    for mod in modfiles:
        k = modfiles.index(mod)
        data, obslist, lbdc, Ra, xmax = readmap(mod, quiet=quiet)
        # It is just taking the length. Ignore the indexes
        pixsize = 2 * xmax[0] / len(data[0, 0, 0, :, :, 0])
        # *60.*60.*1000.*180./_np.pi)
        rad_per_pixel = _np.double( pixsize * _phc.Rsun.cgs / (dist * 
            _phc.pc.cgs))
        lblist = _np.array(Bmax) / _np.array(x) 
        for deg in obs:
            t3m = []
            ob = list(obslist[::2]).index(_phc.find_nearest(obslist[::2], deg))
            if _np.abs(obslist[::2][ob] - deg) > 1:
                print('# The difference of angles is bigger than 1 deg!!!')
                print('# Input: {0:.1f}; Nearest model: {1:.1f}'.format(
                    deg, obslist[::2][ob]))
            for i in range(len(x)):
                j = list(lbdc * 1e-6).index(_phc.find_nearest(lbdc * 1e-6, 
                    lblist[i]))
                if _np.abs(lbdc[j]*1e-6 - lblist[i]) > 0.1*1e-6:
                    print('# The difference of lambd is bigger than 0.1 um!!!')
                    print('# Input: {0:.1f}; Nearest model: {1:.1f}'.format(
                        lblist[i], lbdc[j]*1e-6))
                # lcalc = lbdc[j] * 1e-6
                tmp, V, phvar = fastnumvis3(data[0, ob, j, :, :, 0], lblist[i], 
                    Blist[i], PAlist[i], rad_per_pixel, PAdisk=(216.9 + 90.))
                t3m += [phvar]
            t3m = _np.array(t3m)
            if bindata <= 0:
                ax.plot(x, t3m, color=mcolors[k], markersize=ms, ls='',
                    marker=markers[k], alpha=alp)
                ax2.plot(x, (y - t3m) / yerr, color=mcolors[k],
                    markersize=ms, marker=markers[k], 
                    ls='', alpha=alp)  # alpha=.6,
                if printsum:
                    chi2 = _phc.chi2calc(t3m, y, yerr)
                    ax2.text(0.4, 0.87, r'$\sum={0:.1f}$'.format(chi2), 
                        transform=ax2.transAxes)
            else:
                binned = _phc.bindata(Bmax / lbds, t3m, bindata, xlim=xlim)
                ax.plot(binned[0], binned[1], color=mcolors[k], alpha=.9)
                binned = _phc.bindata(Bmax / lbds, (y - t3m) / yerr, bindata)
                ax2.plot(binned[0], binned[1],
                    color=mcolors[k], markersize=ms, marker='o', ls='')  
                    # alpha=.6,                
    # ax5.get_xaxis().set_ticklabels([])
    # labels = ax5.get_yticks().tolist()
    # labels[-1] = ''
    # ax5.set_yticklabels(labels)
    if xlim is not None:
        ax.set_xlim(xlim)
        ax.set_xlim(xlim)
    # ymax = _np.max(_np.abs(ax5.get_ylim()))
    # ax5.set_ylim([-1.05 * ymax, 1.05 * ymax])
    ax.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax.set_ylabel(u'$\phi_{123}$ (deg.)')
    ax.grid(b=True, linestyle=':', alpha=alp)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_ylim([-6, 6])
    ax2.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax2.set_ylabel(u'$\phi_{123}$(data-mod)/err')
    ax2.grid(b=True, linestyle=':', alpha=alp)
    return ax, ax2


class Disk(object):
    r""" To compute Vieira's models.

    `lbd` in cm (np.array. If None, default value from _phc.BBlbd), `barz2` is 
    the mean value of the square atomic number, `Td` isothermal disk 
    temperature, `fion` is the ionization fraction, mu (:math:`\mu`) is the 
    mean particle weight (mH units), `gb` is the...

    `Ms` is the stellar mass (Msun). """

    def __init__(self, lbd=[.001], barz2=1., Td=12120., fion=1., 
        mu=.5, gb=None, rho=2.8e-11):
        """ kappa class initialiser """
        self.lbd = _np.array(lbd)
        self.Td = Td
        self.barz2 = barz2
        self.fion = fion
        self.rho = rho
        self.mu = mu
        if gb is None:
            # 1e4 = cm to micron
            gb = _np.sum( _phc.gbf(self.Td, self.lbd * 1e4), axis=0 )
        self.gb = gb
        self.update()

    def update(self):
        """ Opacity expression from Brussaard & van de Hulst (1962) """
        self.kappa = 3.692e8 * (1 - _np.exp(-_phc.h.cgs * _phc.c.cgs / 
            self.lbd / _phc.kB.cgs / self.Td)) * self.barz2 * \
            self.Td**-.5 * (self.lbd / _phc.c.cgs)**3 * self.fion * \
            (self.rho / self.mu / _phc.mH.cgs)**2 * (self.gb)       
    # End of class


def I(Disk, Ms=7.7, Teff=20200., Rs=4.94, iang=0., bm2n=-5.5, fmin=5e-3, 
    Rmax=None, px=128, bartau=None):
    r""" I is the specific intensity of the star plus an isothermal disk. The 
    effects of limb-darkening, stellar rotation and circumstellar extinction 
    are neglected. Image constructed from LOWER origin.

    Disk: `lbd` is determined in Disk (lambda vector in cm). 

    `Teff` is the stellar effective temperature (K). All derived quantities are 
    based on this value, but the stellar emission, that is corrected by the 
    `fBBcor` function.

    `Rs`, the stellar (equatorial) radius (Rsun).

    `iang`, inclination angle (deg).

    `bm2n` is the :math:`-2n+\beta < 0` value. 

    `fmin`, minimum flux at the semi-major axis as fraction of the thick disk 
    flux.

    `Rmax`, maximum radius of the images.

    METHOD: If `Rmax` is None, then Rmax of the image is automatically 
    calculated where Athin is equal to `fmin * BBlbd(Tdisk)` at the longest 
    (vector last position) wavelength. Else, the `Rmax` value is taken.
    `Rmax` goes to all wavelengths. 

    `px`, image side in pixels (squared output).

    `bartau` = if None, 1.3 is considered. Otherwise, it is a vector with same 
    dimensions of `Disk.lbd`. See that `bartau` do not changes the 
    emission in Eq. 12(c, Vieira+2015).

    OUTPUT: (squared) images(len(lbd), px, px). The flux unit per pixel is 
    `BBlbd(Teff*fBBcor(Teff))` (cgs), in an area of pixelsize**2. 
    (`pixelsize` info is printed).
    """
    # """ The vertical optical depth. """ 
    # Sound speed for and ideal mono-particle gas (gamma=5/3.)
    csound = _np.sqrt(Disk.Td * _phc.kB.cgs * 5 / 3. / Disk.mu / _phc.mH.cgs)
    # H0 = csound*(G*Ms)**-.5 (A.5, Faes thesis)
    # H0 = csound*(_phc.G.cgs*Ms*_phc.Msun.cgs)**-.5
    H0 = csound * \
        _np.sqrt((_phc.Rsun.cgs * Rs)**3. / _phc.G.cgs / _phc.Msun.cgs / Ms)
    tau0 = _np.sqrt(_np.pi) * H0 * Disk.kappa

    fTd = Disk.Td / Teff
    lbd = Disk.lbd 
    iang = iang * _np.pi / 180.

    BBstar = _phc.BBlbd(Teff * _phc.fBBcor(Teff), lbd=lbd)
    F = _phc.BBlbd(Teff * fTd, lbd=lbd) / BBstar
    if lbd is None:
        bartau = _np.ones(len(lbd)) + .3
        lbd = _np.arange(1000, 10000, 100) * 1e-8  # Angs -> cm
    elif bartau is None:
        bartau = _np.ones(len(lbd)) + .3

    # barR in Stellar Radius
    barR = _np.exp(-_np.log(tau0 / bartau / _np.cos(iang)) / bm2n)
    # in barR units, ie, Rs
    if Rmax is None:
        Rmax = barR[-1] * _np.exp(_np.log(fmin / bartau[-1] / F[-1]) / bm2n)

    print('# Model info (cgs units):\n')
    print(_tab.tabulate(_np.array([lbd, Disk.kappa, tau0, barR]).T,
        headers=['lambda', 'kappa', 'tau0', 'barR (R*)']))
    print('\n# H0 = {0:.3e}'.format(H0))
    print('# Rmax = {0:.2f} R*'.format(Rmax))
    print('# pixelsize: {0:.3f} R*'.format(2*Rmax/px))
    # Create map structure
    xaxis = _np.linspace(-Rmax, Rmax, px)
    yaxis = xaxis / _np.cos(iang)
    X, Y = _np.meshgrid(xaxis, yaxis)
    X, Y2 = _np.meshgrid(xaxis, xaxis)
    R = _np.sqrt(X**2 + Y**2)
    Rrs = _np.sqrt(X**2 + Y2**2)

    if False:
    # if True:
        # Soh para o ultimo lambda...
        image = bartau[-1] * F[-1] * (R / barR[-1])**bm2n
        image[_np.where(R <= barR[-1])] = bartau[-1] * F[-1]
        image[_np.where(R <= 1)] = 1
        # Where Y>0 e R < Rstar, put the star
        image[_np.where((Rrs <= 1) & (Y2 > 0))] = 1
        # Where the disk is tenuous...
        idx = _np.where((Rrs <= 1) & (Y2 < 0) & (R >= barR[-1]))
        if len(idx) > 0:
            image[idx] = 1
    else:
        # Para todos os lambdas...
        image = _np.empty((len(lbd), px, px))
        for i in range(len(lbd)):
            # image[i] = bartau[i] * F[i] * (R / barR[i])**bm2n
            image[i] = F[i] * (R / barR[i])**bm2n
            image[i][_np.where(R <= barR[i])] = F[i]
            image[i][_np.where(R <= 1)] = 1
            # Where Y>0 e R < Rstar, put the star
            image[i][_np.where((Rrs <= 1) & (Y2 > 0))] = 1
            # Where the disk is tenuous...
            idx = _np.where((Rrs <= 1) & (Y2 < 0) & (R >= barR[i]))
            if len(idx) > 0:
                image[i][idx] = 1

    return image


# Old
def plot_pio_res_old(oidata, modellist, outname=None, fmt=['png'], legend=True,
    obsdeg=[60.6], distpc=42.75, quiet=False, xlim=None, bindata=0, 
    PAs=[-180, 180], PAsrev=True, shv2sum=False):
    """ Obs-Model comparison for PIONIER

    `legend`: ?

    `obsdeg`: ?

    `bindata`: ?

    `PArev`: Plot the reverse of the observed PAs

    `shv2sum`: ?
    """
    # Calculate the reverse PAs
    if PAsrev:
        PArv = [0, 0]
        if PAs[0] < 0:
            PArv[0] = PAs[0] + 180
        else:
            # PArv[0] == 0 must be -180!
            PArv[0] = PAs[0] - 180
        if PAs[1] <= 0:
            # PArv[1] == 0 must be +180!
            PArv[1] = PAs[1] + 180
        else:
            PArv[1] = PAs[1] - 180
    # Create array names and colors
    names = []
    Blist = []
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        if label not in names:
            names += [label]
            Blist.append([])
        i = names.index(label)
        Blist[i].extend([_np.sqrt(u**2 + v**2)])
    for i in range(len(Blist)):
        Blist[i] = _np.average(Blist[i])
    names = list(_np.array(names)[_np.argsort(Blist)])
    Blist.sort()
    colors = _phc.gradColor(Blist, min=Blist[0]-(Blist[-1]-Blist[0])*0.1, 
        cmapn='plasma_r')
    # color = _phc.colors[names.index(label)]
    # Create Plot
    # fig, axs = _plt.subplots(3, 2)
    # New plot creating
    fig = _plt.figure()
    axs = [[[], ]*2 for i in range(3)]
    ls, cs = (9, 2)
    axs[0][0] = _plt.subplot2grid((ls, cs), (0, 0), rowspan=2)
    axs[0][1] = _plt.subplot2grid((ls, cs), (0, 1), rowspan=2)
    axs[1][0] = _plt.subplot2grid((ls, cs), (2, 0), rowspan=4)
    axs[1][1] = _plt.subplot2grid((ls, cs), (2, 1), rowspan=4)
    axs[2][0] = _plt.subplot2grid((ls, cs), (6, 0), rowspan=3)
    axs[2][1] = _plt.subplot2grid((ls, cs), (6, 1), rowspan=3)
    # Plot of axis 1 = uv (B/lbd) map
    axs[0][0] = plot_uv_old(ax=axs[0][0], oidata=oidata, PAs=PAs, 
        PAsrev=PAsrev, 
        PArv=PArv, colors=colors, names=names)
    axs[0][1] = plot_baseline(ax=axs[0][1], oidata=oidata, PAs=PAs, 
        PAsrev=PAsrev, PArv=PArv, colors=colors, names=names, legend=legend)
    axs[1][0], axs[1][1] = plot_v2_res(axs[1][0], axs[1][1], oidata=oidata, 
        colors=colors, names=names, modfiles=modellist, obsdeg=obsdeg, 
        dist=distpc, xlim=xlim, PAs=PAs, PAsrev=PAsrev, PArv=PArv, bindata=0, 
        quiet=quiet, alp=.675, printsum=False)
    axs[2][0], axs[2][1] = plot_phi3_res(axs[2][0], axs[2][1], oidata=oidata, 
        colors=colors, names=names, modfiles=modellist, obsdeg=obsdeg, 
        dist=distpc, xlim=xlim, PAs=PAs, PAsrev=PAsrev, PArv=PArv, bindata=0, 
        quiet=quiet, alp=.675, printsum=False)
    # _plt.tight_layout()
    axs[1][1].set_xlim(axs[1][0].get_xlim())
    axs[2][0].set_xlim(axs[1][0].get_xlim())
    axs[2][1].set_xlim(axs[1][0].get_xlim())

    _plt.subplots_adjust(wspace=0.1)
    if outname is None:
        outname = _phc.dtflag()
    _phc.savefig(fig, figname=outname, fmt=fmt)
    return


def plot_uv_old(ax, oidata, colors, names, PAs=[-180, 180], PAsrev=False, 
    PArv=[0, 0], xlim=None):
    """ Plot uv map of a given oidata.vis """
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        ulbd = u / vis2.wavelength.eff_wave
        vlbd = v / vis2.wavelength.eff_wave
        PA = _np.arctan2(u, v) * 180 / _np.pi
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        # [u,-u] = W > E
        # [-u,u] = E < W
        # label=label,
        if (PA >= PAs[0] and PA <= PAs[1]) or (PAsrev and 
            (PA >= PArv[0]) and PA <= PArv[1]):
            ax.plot([ulbd, -ulbd], [vlbd, -vlbd], '.', 
                color=colors[names.index(label)])
    f = 1.4*5e7/112.
    x = 84.06*f
    y = 111.96*f
    ax.arrow(-x, y, 2*x, -2*y, length_includes_head=True,
        head_width=0.1*111.96*2*f, 
        head_length=0.15*111.96*2*f, fc='gray', ec='gray', zorder=1)
    ax.plot([-y, y], [-x, x], color='gray', ls='--')
    ax.set_ylabel(u'B$_{proj}$/$\lambda$')
    # ax2 = ax.twiny()
    if xlim is not None:
        ax.set_xlim(xlim)
    # ax2.set_xlim(ax.get_xlim())
    # ax.get_xaxis().set_ticklabels([])
    ax.xaxis.tick_top()
    ax.axis('equal')
    ax.yaxis.set_major_formatter(_mtick.FormatStrFormatter('%.0e'))
    ax.xaxis.set_major_formatter(_mtick.FormatStrFormatter('%.0e'))
    # _plt.grid(b=True, linestyle=':')
    ax.locator_params(nbins=5)
    return ax


def plot_pionier_old(oidata, ffile='last_run', fmt=['png'], legend=True,
    model=None, obs=None, dist=None):
    """  Standard observational log for PIONIER

    obs is a list
    dist is a number
    """
    fig = _plt.figure()  # figsize=(5.6,8))
    alp = .75
    ms = 3  # markersize
    # xloc = _plt.MaxNLocator(6)
    # ax0 = display info
    # ax0 = fig.add_subplot(211)
    # ax0.axis('off')
    # hdrinfo = oidata.hdrinfo.returninfo()
    # for i in range(4):
        # ax0.text(0., .8-.2*i, hdrinfo[i])
    # ax1 = uvplane/lambda
    ax1 = fig.add_subplot(321)
    colorid = 0
    names = []
    for vis2 in oidata.vis2:
        ulbd = vis2.ucoord / vis2.wavelength.eff_wave
        vlbd = vis2.vcoord / vis2.wavelength.eff_wave
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        if label not in names:
            names += [label]
        color = _phc.colors[names.index(label)]
        # colorid = _np.mod(colorid+1, len(phc.colors))
        # [u,-u] = W > E
        # [-u,u] = E < W
        # label=label,
        ax1.plot([ulbd, -ulbd], [vlbd, -vlbd], '.', color=color)
    ax1.get_xaxis().set_ticklabels([])
    # ax1.xaxis.tick_top()
    # ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    # names = list(_np.unique(names))
    ax1.set_ylabel(u'B$_{proj}$/$\lambda$')
    ax1.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    # if legend: ax1.legend(prop={'size':8},numpoints=1,bbox_to_
        # anchor=(-0.25, 1.0))
    # ax2 = uvplane
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
        # color = phc.colors[colorid]
        # colorid = _np.mod(colorid+1, len(phc.colors))
        # [u,-u] = W > E
        # [-u,u] = E < W
        if label not in leg:
            ax2.plot([u, -u], [v, -v], '.', label=label, color=color)
            leg.append(label)
        else:
            ax2.plot([u, -u], [v, -v], '.', color=color)  # label=label,
        # names.append(vis2.target.target)
    # names = list(_np.unique(names))
    ax2.xaxis.tick_top()
    ax2.set_ylabel(u'B$_{proj}$ (m)')
    ax2.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    if legend:
        ax2.legend(prop={'size': 8}, numpoints=1, bbox_to_anchor=(1.05, 1.0))
    # ax3 = VIS2 vs. B
    # ax4 = VIS2 vs. PA
    # names = []
    colorid = 0
    plotid = 323
    ax3 = fig.add_subplot(plotid)
    plotid = 324
    ax4 = fig.add_subplot(plotid)
    for vis2 in oidata.vis2:
        u = vis2.ucoord / vis2.wavelength.eff_wave
        v = vis2.vcoord / vis2.wavelength.eff_wave
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = _phc.colors[names.index(label)]
        # color = phc.colors[colorid]
        # colorid = _np.mod(colorid + 1, len(phc.colors))
        line = ax3.errorbar(_np.sqrt(u**2 + v**2),
        vis2.vis2data, yerr=vis2.vis2err, color=color, fmt='o', markersize=ms)
            # , label=label)
        # _np.arctan(self.ucoord / self.vcoord) * 180.0 / _np.pi % 180.0
        PAobs = _np.arctan2(u, v) * 180.0 / _np.pi
        # idx = _np.where(PAobs < 0)
        # PAobs[idx] = PAobs[idx]+180
        line = ax4.errorbar(PAobs,
        vis2.vis2data, yerr=vis2.vis2err, color=color, fmt='o', markersize=ms)
        # , label=label)
    Blim = ax3.get_xlim()
    if model is not None:
        # res = 20
        # V2 = _np.empty(res)
        for mod in model:
            if obs is None:
                obs = [0]
            data, obslist, lbdc, Ra, xmax = readmap(mod)
            pixsize = 2 * xmax[0] / len(data[0, 0, 0, :, :])
            # *60.*60.*1000.*180./_np.pi)
            rad_per_pixel = _np.double(
                pixsize * _phc.Rsun.cgs / (dist * _phc.pc.cgs))
            for vis2 in oidata.vis2:
                if (vis2.station[0] and vis2.station[1]):
                    label = vis2.station[0].sta_name + vis2.station[1].sta_name
                else:
                    label = 'unnamed'
                color = _phc.colors[names.index(label)]
                u = vis2.ucoord
                v = vis2.vcoord
                B = _np.sqrt(u**2 + v**2)
                PA = _np.arctan2(u, v) * 180.0 / _np.pi
                avlbd = _np.average(vis2.wavelength.eff_wave)
                V2 = []
                lbds = []
                for i in range(len(vis2.wavelength.eff_wave)):
                    for ob in obs:
                        lbd = vis2.wavelength.eff_wave[i]
                        j = list(lbdc * 1e-6).index(_phc.find_nearest(lbdc *
                                1e-6, lbd))
                        # print lbdc[j]*1e-6-lbd
                        lbcalc = lbdc[j] * 1e-6
                        lbcalc = lbd
                        tmp, V, phvar = fastnumvis(data[ob, 0, j, :, :],
                            lbcalc, B, PA, rad_per_pixel, PAdisk=(216.9 + 90.))
                        V2 += [V**2]
                        lbds += [lbcalc]
                # ax3.plot(B/vis2.wavelength.eff_wave, V2, color='purple',
                    # alpha=.3)
                ax3.plot(B / lbds, V2, color='purple', alpha=.3)
                ax4.plot(_np.tile(PA, len(lbds)), V2, color='purple', alpha=.3,
                    marker='s', markersize=ms)
                # lbd = phc.find_nearest(lbdc*1e-6,avlbd)
                # j = list(lbdc*1e-6).index(lbd)
                # Bs = _np.linspace(Blim[0]*lbd, Blim[1]*lbd, res)
                # for i in range(res):
                    # tmp, V2[i], phvar = fastnumvis(data[0,0,j,:,:], lbd,
                        # Bs[i], PA, rad_per_pixel, PAdisk=36.9+90)
                # ax3.plot(Bs/lbd, V2**2, color=color)
                # print Bs/lbd
                # print V2
                # a = raw_input('asdads')
    ax3.set_xlim(Blim)
    # ax1.xaxis.get_major_formatter().set_useOffset(False)
    ax3.xaxis.set_major_formatter(_mtick.FormatStrFormatter('%.2e'))
    # ax3.get_xaxis().set_visible(False)
    # ax4.get_xaxis().set_visible(False)
    ax3.get_xaxis().set_ticklabels([])
    PAlim = [-180, 180]
    ax3.set_ylim([0, 1.1])
    ax3.set_ylabel(u'$V$ $^2$')
    ax3.grid(b=True, linestyle=':', alpha=alp)
    ax4.set_ylim([0, 1.1])
    ax4.get_xaxis().set_ticklabels([])
    ax4.set_ylabel(u'$V$ $^2$')
    ax4.grid(b=True, linestyle=':', alpha=alp)
    # ax5 = T3PHI vs. B
    # ax6 = T3PHI vs. PA
    # names = []
    colorid = 0
    plotid = 325
    ax5 = fig.add_subplot(plotid)
    plotid = 326
    ax6 = fig.add_subplot(plotid)
    ax5.plot(Blim, [0, 0], ls='--')
    ax6.plot(PAlim, [0, 0], ls='--')
    for t3 in oidata.t3:
        u1 = t3.u1coord
        v1 = t3.v1coord
        u2 = t3.u2coord
        v2 = t3.v2coord
        if _np.sqrt(u1**2 + v1**2) > _np.sqrt(u2**2 + v2**2):
            u = u1
            v = v1
        else:
            u = u2
            v = v2
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u, v) * 180.0 / _np.pi
        B = _np.tile(B, len(t3.wavelength.eff_wave))
        PA = _np.tile(PA, len(t3.wavelength.eff_wave))
        # idx = _np.where(PA < 0)
        # PA[idx] = PA[idx]+180
        # if (t3.station[0] and t3.station[1]):
            # label = t3.station[0].sta_name + t3.station[1].sta_name
        # else:
            # label = 'unnamed'
        # color = names.index(label)
        # color = phc.colors[colorid]
        # colorid = _np.mod(colorid + 1, len(phc.colors))
        color = 'Black'
        # y = _np.repeat(t3.t3phi, len(t3.wavelength.eff_wave))
        # yerr = _np.repeat(t3.t3phierr, len(t3.wavelength.eff_wave))
        y = t3.t3phi
        yerr = t3.t3phierr
        line = ax5.errorbar(B / t3.wavelength.eff_wave, y, yerr=yerr,
            color=color, fmt='o', markersize=ms)  # , label=label)
        # , label=label)
        line = ax6.errorbar(
            PA, y, yerr=yerr, color=color, fmt='o', markersize=ms)
    if model is not None:
        for mod in model:
            if obs is None:
                obs = [0]
            data, obslist, lbdc, Ra, xmax = readmap(mod)
            pixsize = 2 * xmax[0] / len(data[0, 0, 0, :, :])
            # *60.*60.*1000.*180./_np.pi)
            rad_per_pixel = _np.double(
                pixsize * _phc.Rsun.cgs / (dist * _phc.pc.cgs))
            for t3 in oidata.t3:
                u1 = t3.u1coord
                v1 = t3.v1coord
                u2 = t3.u2coord
                v2 = t3.v2coord
                B = _np.append(
                    _np.sqrt(u1**2 + v1**2), _np.sqrt(u2**2 + v2**2))
                PA = _np.append(_np.arctan2(u1, v1) * 180.0 / _np.pi,
                    _np.arctan2(u2, v2) * 180.0 / _np.pi)
                Bmax = []
                PAmax = []
                t3m = []
                lbds = []
                for i in range(len(t3.wavelength.eff_wave)):
                    for ob in obs:
                        lbd = t3.wavelength.eff_wave[i]
                        j = list(lbdc * 1e-6).index(_phc.find_nearest(lbdc *
                            1e-6, lbd))
                        lcalc = lbdc[j] * 1e-6
                        lcalc = lbd
                        tmp, V, phvar = fastnumvis3(data[ob, 0, j, :, :],
                            lcalc, B, PA, rad_per_pixel, PAdisk=(216.9 + 90.))
                        t3m += [phvar]
                        Bmax += [_np.max(B)]
                        if _np.max(B) == B[0]:
                            PAmax += [PA[0]]
                        else:
                            PAmax += [PA[1]]
                        lbds += [lcalc]
                Bmax = _np.array(Bmax)
                lbds = _np.array(lbds)
                t3m = _np.array(t3m)
                PAmax = _np.array(PAmax)
                ax5.plot(Bmax / lbds, t3m, color='purple', alpha=.9)
                ax6.plot(PAmax, t3m, color='purple', alpha=.6, marker='s',
                    markersize=ms)
    # ax5.get_xaxis().set_ticklabels([])
    ax5.set_xlim(Blim)
    ax6.set_xlim(PAlim)
    ymax = _np.max(_np.abs(ax5.get_ylim()))
    ax5.set_ylim([-1.05 * ymax, 1.05 * ymax])
    ax6.set_ylim([-1.05 * ymax, 1.05 * ymax])
    ax5.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax5.set_ylabel(u'$\phi_{123}$ (deg.)')
    ax5.grid(b=True, linestyle=':', alpha=alp)
    ax6.set_xlabel(u'$PA$ (deg.)')
    ax6.set_ylabel(u'$\phi_{123}$ (deg.)')
    ax6.grid(b=True, linestyle=':', alpha=alp)
    # SAVING
    dir, name = _phc.trimpathname(ffile)
    name = _phc.rmext(name)
    # _plt.savefig('hdt/{}_{}.png'.format(hdrinfo[0], hdrinfo[2]),
        # transparent=True)
    # _plt.locator_params(axis = 'x', nbins = 7)
    _plt.subplots_adjust(left=0.12, right=0.95, top=0.96, bottom=0.09,
        hspace=.009, wspace=.32)
    if not _os.path.exists('hdt'):
        _os.system('mkdir hdt')
    for suf in fmt:
        _plt.savefig('hdt/{0}.{1}'.format(name, suf), transparent=True)
    _plt.close()
    return


def printinfo(file, extract=False):
    """ Print AMBER OIFITS observational info.

    If `extract` is False, output is:
    - DATE-OBS, MJD, Target, B1, PA1, B2, PA2, B3, PA3.

    If True, output is:
    - [DATE-OBS, MJD, Target, B1, PA1, B2, PA2, B3, PA3], WAVE, DPlist, V2list

    Where, DPlist = [DP1, eDP1, DP2, eDP2, DP3, eDP3] and V2list = [V2_1, 
    eV2_1, V2_2, eV2_2, V2_3, eV2_3].
    """
    oidata = _oifits.open(file, quiet='True')
    info = list(oidata.hdrinfo.returninfo())
    if extract:
        wav = []
        info2 = []
        info3 = []
        info4 = []
        for vis in oidata.vis:
            fact = 1.
            PA = _np.arctan2(vis.ucoord, vis.vcoord) * 180.0 / _np.pi % 360.0
            if PA > 180:
                PA = PA % 180
                fact = -1.
            info2 += ['{0:.1f}'.format(_np.sqrt(vis.ucoord**2 +
                                                vis.vcoord**2))]
            info2 += ['{0:.1f}'.format(PA)]
            wav += [1e6 * vis.wavelength.eff_wave]
            info3 += [fact * vis.visphi, vis.visphierr]
        for vis2 in oidata.vis2:
            info4 += [vis2.vis2data, vis2.vis2err]
        return [info[2][:10], '{0:.7f}'.format(info[1]), info[0]] + \
            info2, wav, info3, info4
    else:
        info2 = []
        for vis in oidata.vis:
            info2 += ['{0:.1f}'.format(_np.sqrt(vis.ucoord**2 + 
                vis.vcoord**2))]
            # info2+= ['{0:.1f}'.format(_np.arctan2(vis.ucoord , vis.vcoord) *
                # 180.0 / _np.pi % 180.0)]
            info2 += ['{0:.1f}'.format(_np.arctan2(vis.ucoord, vis.vcoord) *
                180.0 / _np.pi % 180.0)]
            # print vis.ucoord, vis.vcoord        
        return [info[2][:10], '{0:.7f}'.format(info[1]), info[0]] + info2


def plot_oifits(oidata, ffile='last_run', fmt=['png'], xlim=None,
    legend=True):
    """ Standard observational log for AMBER

    If the file starts with "PRODUCT_", it searchs for the specs in the "AVG"
    folder.

    (One could write this info into the fits file. Since I've only tested the
    reading features of the `oifits` routine, I prefered do it this way).
    """
    # If it is a PRODUCT ffile, tries to load the AVG spec.
    specfile = ''
    if ffile.find('_PRO/PRODUCT') > 0:
        dateobs = ffile[ffile.find('_20') + 1:ffile.find('_20') + 20]
        dateobs2 = _phc.strrep(dateobs, -3, ':')
        dateobs2 = _phc.strrep(dateobs2, -6, ':')
        specfile = ffile[:ffile.find('_PRO/PRODUCT_')].replace('_SPEC', '') +\
            '/*{0}*_OIDATA_AVG.fits*'.format(dateobs)
        specfile = _glob(specfile)
        if len(specfile) == 0:
            specfile = ffile[:ffile.find('_PRO/PRODUCT_')].replace('_SPEC',
                '') + '/*{0}*_OIDATA_AVG.fits*'.format(dateobs2)
            specfile = _glob(specfile)
        if len(specfile) != 1:
            specfile = ''
            print('# ERROR! This is a PRO oifits and the AVG file was ' +
            'not found!')
            print(ffile[:ffile.find('_PRO/PRODUCT_')] + ('/*{0}*_OIDATA_' +
        'AVG.fits*').format(dateobs))
        else:
            specfile = specfile[0]
            specoidata = _oifits.open(specfile, quiet=True)
            # print('# {0} file read!!!'.format(specfile))
            oidata.amberspec = specoidata.amberspec
            spec = oidata.amberspec[0]
            vis = oidata.vis[0]
            if len(spec.wavelength.eff_wave) == len(vis.wavelength.eff_wave):
                for i in range(len(oidata.vis)):
                    oidata.amberspec[i].wavelength = oidata.vis[i].wavelength
            else:
                print('# ERROR! spec and vis sizes are different for {0}'.
                format(_phc.trimpathname(ffile)[1]))
    #
    fig = _plt.figure(figsize=(5.6, 8))
    alp = .75
    # xloc = _plt.MaxNLocator(6)
    # ax0 = display info
    ax0 = fig.add_subplot(521)
    ax0.axis('off')
    hdrinfo = oidata.hdrinfo.returninfo()
    for i in range(4):
        ax0.text(0., .8 - .2 * i, hdrinfo[i])
    # ax2 = uvplane.
    ax2 = fig.add_subplot(522)
    colorid = 0
    names = []
    for vis in oidata.vis:
        if xlim is None:
            xmin = None
            xmax = None
            xmin = _np.amin(_np.append(1e6 * vis.wavelength.eff_wave[
                _np.where(vis.flag is False)], xmax))
            xmax = _np.amax(_np.append(1e6 * vis.wavelength.eff_wave[
                _np.where(vis.flag is False)], xmin))
            xlim = (xmin, xmax)
        u = vis.ucoord
        v = vis.vcoord
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        colorid = _np.mod(colorid + 1, len(colors))
        # [u,-u] = W > E
        # [-u,u] = E < W
        ax2.plot([-u, u], [v, -v], '.', label=label, color=color)
        ax2.xaxis.tick_top()
        names.append(vis.target.target)
    ax2.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    if legend:
        ax2.legend(prop={'size': 8}, numpoints=1, bbox_to_anchor=(-0.25, 1.0))
    # ax3 = dif.phases, RIGHT COLUMN
    plotid = (5,2,3)
    names = []
    colorid = 0
    yrange = [0, 0]
    for vis in oidata.vis:
        diff = _np.max(_np.abs(vis.visphi))
        if diff > yrange[1]:
            yrange = [-diff, diff]
    for vis in oidata.vis:
        ax3 = fig.add_subplot(*plotid)
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        colorid = _np.mod(colorid + 1, len(colors))
        line = ax3.errorbar(1e6 * vis.wavelength.eff_wave, vis.visphi,
            vis.visphierr, label=label, color=color)
        ax3.set_ylim(yrange)
        ax3.set_xlim(xlim)
        ax3.set_ylabel(u'$\phi$ (deg.)')
        names.append(vis.target.target)
        names = list(_np.unique(names))
        # title = names.pop()
        # for name in names:
        #    title += ', %s'%(name)
        # ax1.set_title(title)
        # ax1.set_ylabel('Differential phase')
        plotid += 2
        # ax3.get_xaxis().set_visible(False)
        ax3.get_xaxis().set_ticklabels([])
        _plt.grid(b=True, linestyle=':', alpha=alp)
    # ax5 = closure phases
    # plotid = (5,2,10)
    names = []
    colorid = 3
    for t3 in oidata.t3:
        ax5 = fig.add_subplot(5, 2, 10)
        if (t3.station[0] and t3.station[1] and  t3.station[2]):
            label = t3.station[0].sta_name + t3.station[1].sta_name + \
                t3.station[2].sta_name 
        else:
            label = 'unnamed'
        color = colors[colorid]
        line = ax5.errorbar(1e6 * t3.wavelength.eff_wave, t3.t3phi,
        t3.t3phierr, label=label, color=color)
        ax5.set_xlim(xlim)
        # diff = _np.max(_np.abs(t3.t3phi))
        # yrange2 = [-diff, diff]
        ax5.set_ylim(yrange)
        names.append(t3.target.target)
        names = list(_np.unique(names))
    ax5.set_ylabel(u'Closure $\phi$ (deg.)')
    ax5.set_xlabel('Wavelength ($\mu$m)')
    _plt.setp( ax5.xaxis.get_majorticklabels(), rotation=-35 )
    _plt.grid(b=True, linestyle=':', alpha=alp)
    # ax4 = visibilities, LEFT COLUMN
    plotid = (5,2,4)
    names = []
    colorid = 0
    yrange = [1, 0]
    for vis in oidata.vis2:
        yrange = [ _np.min([yrange[0], _np.min(vis.vis2data)]),
        _np.max([yrange[1], _np.max(vis.vis2data)]) ]
    if yrange[1] > 1.1:
        yrange[1] = 1.1
    if yrange[0] < 0 or yrange[0] >= 1:
        yrange[0] = 0
    for vis in oidata.vis2:
        ax4 = fig.add_subplot(*plotid)
        if (vis.station[0] and vis.station[1]):
            label = vis.station[0].sta_name + vis.station[1].sta_name
        else:
            label = 'unnamed'
        color = colors[colorid]
        colorid = _np.mod(colorid + 1, len(colors))
        line = ax4.errorbar(1e6 * vis.wavelength.eff_wave, vis.vis2data,
            vis.vis2err, label=label, color=color)
        names.append(vis.target.target)
        names = list(_np.unique(names))
        ax4.set_ylim(yrange)
        ax4.set_xlim(xlim)
        ax4.set_ylabel(u'V$^{2}$')
        # title = names.pop()
        # for name in names:
            # title += ', %s'%(name)
        # ax1.set_title(title)
        # ax1.set_ylabel('Differential phase')
        # plotid += 2
        ax4.get_xaxis().set_ticklabels([])
        _plt.grid(b=True, linestyle=':', alpha=alp)
    # ax1 = Line profile
    if True:
        ax1 = fig.add_subplot(5, 2, 9)
        colorid = 0
        names = []
        for spec in oidata.amberspec:
            label = 'unnamed'
            color = colors[colorid]
            # x,y,yerr = linfit(1e6*spec.wavelength.eff_wave, spec.spectrum,
                # yerr=spec.spectrumerr)
            # ax1.errorbar(x, y, yerr, label=label, color=color)
            x = 1e6 * spec.wavelength.eff_wave
            y = _linfit(1e6 * spec.wavelength.eff_wave, spec.spectrum)
            ax1.plot(x, y, label=label, color=color)
            colorid = _np.mod(colorid + 1, len(colors))
            # [-u,u] = W > E
            # [u,-u] = E < W
        ax1.set_xlabel('Wavelength ($\mu$m)')
        ax1.set_xlim(xlim)
        ax1.set_ylabel('Norm. flux')
        _plt.setp( ax1.xaxis.get_majorticklabels(), rotation=-35 )
        _plt.grid(b=True, linestyle=':', alpha=alp)
    _phc.outfld()
    dir, name = _phc.trimpathname(ffile)
    name = _phc.rmext(name)
    # _plt.savefig('hdt/{}_{}.png'.format(hdrinfo[0], hdrinfo[2]), 
        # transparent=True)
    # _plt.locator_params(axis = 'x', nbins = 7)
    _plt.subplots_adjust(
        left=0.12, right=0.95, top=0.96, bottom=0.09, hspace=.009, wspace=.32)
    for suf in fmt:
        _plt.savefig('hdt/{0}.{1}'.format(name, suf), transparent=True)
    _plt.close()
    return


def plot_pionier_res(oidata, model, outname=None, fmt=['png'], legend=True,
    obs=None, dist=42.75, quiet=True, xlim=None, bindata=0, 
    PArange=[-180, 180], PArr=True, shv2sum=False):
    """  Obs-Model comparison for PIONIER

    model, obs are lists
    dist is a number
    `PArange` is in deg., from -180 to 180.
    """
    # PArange = _np.array(PArange)*2*_np.pi/180
    v2sum = 0.
    if outname is None:
        outname = _phc.dtflag()
    if PArr:
        PArv = [0, 0]
        if PArange[0] < 0:
            PArv[0] = PArange[0] + 180
        else:
            # PArv[0] == 0 must be -180!
            PArv[0] = PArange[0] - 180
        if PArange[1] <= 0:
            # PArv[1] == 0 must be +180!
            PArv[1] = PArange[1] + 180
        else:
            PArv[1] = PArange[1] - 180
    fig = _plt.figure()  # figsize=(5.6,8))
    alp = .75
    ms = 3  # markersize
    # xloc = _plt.MaxNLocator(6)
    # ax0 = display info
    # ax0 = fig.add_subplot(211)
    # ax0.axis('off')
    # hdrinfo = oidata.hdrinfo.returninfo()
    # for i in range(4):
        # ax0.text(0., .8-.2*i, hdrinfo[i])
    # ax1 = uvplane/lambda
    ax1 = fig.add_subplot(321)
    colorid = 0
    names = []
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        ulbd = u / vis2.wavelength.eff_wave
        vlbd = v / vis2.wavelength.eff_wave
        PA = _np.arctan2(u, v) * 180 / _np.pi
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        if label not in names:
            names += [label]
        color = _phc.colors[names.index(label)]
        # colorid = _np.mod(colorid+1, len(phc.colors))
        # [u,-u] = W > E
        # [-u,u] = E < W
        # label=label,
        if (PA >= PArange[0] and PA <= PArange[1]) or (PArr and 
            (PA >= PArv[0]) and PA <= PArv[1]):
            ax1.plot([ulbd, -ulbd], [vlbd, -vlbd], '.', color=color)
    ax1.get_xaxis().set_ticklabels([])
    # ax1.xaxis.tick_top()
    # ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    # names = list(_np.unique(names))
    ax1.set_ylabel(u'B$_{proj}$/$\lambda$')
    ax1.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    # if legend: ax1.legend(prop={'size':8},numpoints=1,bbox_to_anchor=(-0.25,
        # 1.0))
    # ax2 = uvplane
    ax2 = fig.add_subplot(322)
    colorid = 0
    leg = []
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        PA = _np.arctan2(u, v) * 180 / _np.pi
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = _phc.colors[names.index(label)]
        # color = phc.colors[colorid]
        # colorid = _np.mod(colorid+1, len(phc.colors))
        # [u,-u] = W > E
        # [-u,u] = E < W
        if (PA >= PArange[0] and PA <= PArange[1]) or (PArr and 
            (PA >= PArv[0]) and PA <= PArv[1]):
            if label not in leg:
                ax2.plot([u, -u], [v, -v], '.', label=label, color=color)
                leg.append(label)
            else:
                ax2.plot([u, -u], [v, -v], '.', color=color)  # label=label,
        # names.append(vis2.target.target)
    # names = list(_np.unique(names))
    ax2.xaxis.tick_top()
    ax2.set_ylabel(u'B$_{proj}$ (m)')
    if xlim is not None:
        ax2.set_xlim(xlim)
    ax2.axis('equal')
    _plt.grid(b=True, linestyle=':', alpha=alp)
    if legend:
        ax2.legend(prop={'size': 8}, numpoints=1, bbox_to_anchor=(1.05, 1.0),
            framealpha=0.5, fancybox=True)
    # ax3 = VIS2 vs. B
    # ax4 = VIS2 vs. PA
    # names = []
    colorid = 0
    mcolors = ['black', 'red', 'green', 'blue']
    plotid = 323
    ax3 = fig.add_subplot(plotid)
    plotid = 324
    ax4 = fig.add_subplot(plotid)
    ax4.plot([0, 1e8], [0, 0], ls='--')
    print('ax4')
    for vis2 in oidata.vis2:
        u = vis2.ucoord
        v = vis2.vcoord
        if (vis2.station[0] and vis2.station[1]):
            label = vis2.station[0].sta_name + vis2.station[1].sta_name
        else:
            label = 'unnamed'
        color = _phc.colors[names.index(label)]
        # color = phc.colors[colorid]
        # colorid = _np.mod(colorid + 1, len(phc.colors))
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u, v) * 180 / _np.pi
        if (PA >= PArange[0] and PA <= PArange[1]) or (PArr and 
            (PA >= PArv[0]) and PA <= PArv[1]):
            if bindata <= 0:
                line = ax3.errorbar(B / vis2.wavelength.eff_wave, 
                    vis2.vis2data, yerr=vis2.vis2err, fmt='o', markersize=ms, 
                    color=color)
            else:
                binned = _phc.bindata(B / vis2.wavelength.eff_wave, 
                    vis2.vis2data, bindata, yerr=vis2.vis2err, xlim=xlim)
                line = ax3.errorbar(binned[0], binned[1], yerr=binned[2], 
                    fmt='o', markersize=ms, color=color)
            # , label=label)
        # _np.arctan(self.ucoord / self.vcoord) * 180.0 / _np.pi % 180.0
        if obs is None:
            obs = [0]
        for mod in model:
            k = model.index(mod)
            data, obslist, lbdc, Ra, xmax = readmap(mod, quiet=quiet)
            pixsize = 2 * xmax[0] / len(data[0, 0, 0, :, :, 0])
            # *60.*60.*1000.*180./_np.pi)
            rad_per_pixel = _np.double(
                pixsize * _phc.Rsun.cgs / (dist * _phc.pc.cgs))
            V2 = []
            lbds = []
            for i in range(len(vis2.wavelength.eff_wave)):
                for ob in obs:
                    lbd = vis2.wavelength.eff_wave[i]
                    j = list(
                        lbdc * 1e-6).index(_phc.find_nearest(lbdc * 1e-6, lbd))
                    # lbcalc = lbdc[j]*1e-6
                    lbcalc = lbd
                    tmp, V, phvar = fastnumvis(data[0, ob, j, :, :, 0], lbcalc, 
                        B, PA, rad_per_pixel, PAdisk=(216.9 + 90.))
                    V2 += [V**2]
                    lbds += [lbcalc]
            lbds = _np.array(lbds)
            # V2 = _np.array(V2)+(.0875e-8*B/lbds-.026)
            V2 = _np.array(V2) + (.0875e-8 * B / lbds - .026)
            if (PA >= PArange[0] and PA <= PArange[1]) or (PArr and 
                (PA >= PArv[0]) and PA <= PArv[1]):
                if bindata <= 0:
                    v2sum += (vis2.vis2data - V2) / vis2.vis2err
                    ax3.plot(B / lbds, V2, color=mcolors[k], alpha=.3)
                    ax4.plot(B / lbds, (vis2.vis2data - V2) / vis2.vis2err,
                        color=mcolors[k], markersize=ms, marker='o', ls='')
                else:
                    binned = _phc.bindata(B / lbds, V2, bindata, xlim=xlim)
                    ax3.plot(binned[0], binned[1], color=mcolors[k], alpha=.3)
                    binned = _phc.bindata(B / lbds, (vis2.vis2data - V2) / 
                        vis2.vis2err, bindata, xlim=xlim)
                    ax4.plot(binned[0], binned[1],
                        color=mcolors[k], markersize=ms, marker='o', ls='')
                # marker='s',
    Blim = ax3.get_xlim()
    ax4.set_xlim(Blim)
    # ax1.xaxis.get_major_formatter().set_useOffset(False)
    ax3.xaxis.set_major_formatter(_mtick.FormatStrFormatter('%.2e'))
    # ax3.get_xaxis().set_visible(False)
    # ax4.get_xaxis().set_visible(False)
    ax3.get_xaxis().set_ticklabels([])
    PAlim = [-180, 180]
    ax3.set_ylim([0, 1.1])
    ax3.set_ylabel(u'$V$ $^2$')
    if xlim is not None:
        ax3.set_xlim(xlim)
    ax3.grid(b=True, linestyle=':', alpha=alp)
    ax4.set_ylim([-6, 6])
    ax4.get_xaxis().set_ticklabels([])
    ax4.set_ylabel(u'$V$ $^2$(data-mod)/err')
    if xlim is not None:
        ax4.set_xlim(xlim)
    ax4.grid(b=True, linestyle=':', alpha=alp)
    if shv2sum:
        v2sum = _np.sum(v2sum)
        print(v2sum)
        ax4.text(0.4, 0.87, r'$\sum={0:.1f}$'.format(v2sum), transform=ax4.transAxes)
        # ax4.text(0.4, 0.8, 'ADS', transform=ax4.transAxes)
    # ax5 = T3PHI vs. B
    # ax6 = T3PHI vs. PA
    # names = []
    colorid = 0
    plotid = 325
    ax5 = fig.add_subplot(plotid)
    plotid = 326
    ax6 = fig.add_subplot(plotid)
    ax5.plot(Blim, [0, 0], ls='--')
    ax6.plot(Blim, [0, 0], ls='--')
    print('ax6')
    for t3 in oidata.t3:
        u1 = t3.u1coord
        v1 = t3.v1coord
        u2 = t3.u2coord
        v2 = t3.v2coord
        if _np.sqrt(u1**2 + v1**2) > _np.sqrt(u2**2 + v2**2):
            u = u1
            v = v1
        else:
            u = u2
            v = v2
        B = _np.sqrt(u**2 + v**2)
        PA = _np.arctan2(u, v) * 180.0 / _np.pi
        B = _np.tile(B, len(t3.wavelength.eff_wave))
        PA = _np.tile(PA, len(t3.wavelength.eff_wave))
        # idx = _np.where(PA < 0)
        # PA[idx] = PA[idx]+180
        # if (t3.station[0] and t3.station[1]):
            # label = t3.station[0].sta_name + t3.station[1].sta_name
        # else:
            # label = 'unnamed'
        # color = names.index(label)
        # color = phc.colors[colorid]
        # colorid = _np.mod(colorid + 1, len(phc.colors))
        color = 'Black'
        # y = _np.repeat(t3.t3phi, len(t3.wavelength.eff_wave))
        # yerr = _np.repeat(t3.t3phierr, len(t3.wavelength.eff_wave))
        y = t3.t3phi
        yerr = t3.t3phierr
        if (PA[0] >= PArange[0] and PA[0] <= PArange[1]) or (PArr and 
            PA[0] >= PArv[0] and PA[0] <= PArv[1]):
            if bindata <= 0:
                line = ax5.errorbar(B / t3.wavelength.eff_wave, y, yerr=yerr,
                    color=color, fmt='o', markersize=ms)
            else:
                binned = _phc.bindata(B / t3.wavelength.eff_wave, y, 
                    bindata, yerr=yerr, xlim=xlim)
                line = ax5.errorbar(binned[0], binned[1], yerr=binned[2], 
                    color=color, fmt='o', markersize=ms)
            # , label=label)
        if obs is None:
            obs = [0]
        for mod in model:
            k = model.index(mod)
            data, obslist, lbdc, Ra, xmax = readmap(mod, quiet=True)
            pixsize = 2 * xmax[0] / len(data[0, 0, 0, :, :, 0])
            # *60.*60.*1000.*180./_np.pi)
            rad_per_pixel = _np.double(
                pixsize * _phc.Rsun.cgs / (dist * _phc.pc.cgs))
            B = _np.append(_np.sqrt(u1**2 + v1**2), _np.sqrt(u2**2 + v2**2))
            PA = _np.append(_np.arctan2(u1, v1) * 180.0 / _np.pi, _np.arctan2(
                u2, v2) * 180.0 / _np.pi)
            Bmax = []
            PAmax = []
            t3m = []
            lbds = []
            for i in range(len(t3.wavelength.eff_wave)):
                for ob in obs:
                    lbd = t3.wavelength.eff_wave[i]
                    j = list(
                        lbdc * 1e-6).index(_phc.find_nearest(lbdc * 1e-6, lbd))
                    lcalc = lbdc[j] * 1e-6
                    lcalc = lbd
                    tmp, V, phvar = fastnumvis3(data[0, ob, j, :, :, 0], lcalc, 
                        B, PA, rad_per_pixel, PAdisk=(216.9 + 90.))
                    t3m += [phvar]
                    Bmax += [_np.max(B)]
                    if _np.max(B) == B[0]:
                        PAmax += [PA[0]]
                    else:
                        PAmax += [PA[1]]
                    lbds += [lcalc]
            Bmax = _np.array(Bmax)
            lbds = _np.array(lbds)
            t3m = _np.array(t3m)
            PAmax = _np.array(PAmax)
            # print PAmax
            if (PA[0] >= PArange[0] and PA[0] <= PArange[1]) or (PArr and 
                PA[0] >= PArv[0] and PA[0] <= PArv[1]):
                if bindata <= 0:
                    ax5.plot(Bmax / lbds, t3m, color=mcolors[k], alpha=.9)
                    ax6.plot(Bmax / lbds, (y - t3m) / yerr, color=mcolors[k],
                        markersize=ms, marker='o', ls='')  # alpha=.6,
                else:
                    binned = _phc.bindata(Bmax / lbds, t3m, bindata, xlim=xlim)
                    ax5.plot(binned[0], binned[1], color=mcolors[k], alpha=.9)
                    binned = _phc.bindata(Bmax / lbds, (y - t3m) / yerr, bindata)
                    ax6.plot(binned[0], binned[1],
                        color=mcolors[k], markersize=ms, marker='o', ls='')  
                        # alpha=.6,                
    # ax5.get_xaxis().set_ticklabels([])
    # labels = ax5.get_yticks().tolist()
    # labels[-1] = ''
    # ax5.set_yticklabels(labels)
    ax5.set_xlim(Blim)
    ax6.set_xlim(Blim)
    ymax = _np.max(_np.abs(ax5.get_ylim()))
    ax5.set_ylim([-1.05 * ymax, 1.05 * ymax])
    ax6.set_ylim([-6, 6])
    ax5.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax5.set_ylabel(u'$\phi_{123}$ (deg.)')
    if xlim is not None:
        ax5.set_xlim(xlim)
    ax5.grid(b=True, linestyle=':', alpha=alp)
    ax6.set_xlabel(u'B$_{proj}$/$\lambda$')
    ax6.set_ylabel(u'$\phi_{123}$(data-mod)/err')
    if xlim is not None:
        ax6.set_xlim(xlim)
    ax6.grid(b=True, linestyle=':', alpha=alp)
    # SAVING
    # dir, name = _phc.trimpathname(ffile)
    # name = _phc.rmext(name)
    # _plt.savefig('hdt/{}_{}.png'.format(hdrinfo[0], hdrinfo[2]),
        # transparent=True)
    # _plt.locator_params(axis = 'x', nbins = 7)
    _plt.subplots_adjust(
        left=0.12, right=0.95, top=0.96, bottom=0.09, hspace=.009, wspace=.32)
    print('asdhasduhasd')
    _phc.savefig(fig, figname=outname, fmt=fmt)
    return


if __name__ == "__main__":
    pass
