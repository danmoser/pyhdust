# -*- coding:utf-8 -*-

"""PyHdust *basoptics* module: Basic Optics functions.

Some refs:
- https://www.telescope-optics.net

Orion:
- Do = 114mm = 4.5 inch
- fo = 450mm  # f4
- De = 31.8mm = 1.25 inch
- fe = 17mm / 06mm
- mag = 26 / 75
- AFOV = 62. / 79.4
- pupil size = 4.25 mm / 1.5 mm
- True FoV = 2.34 deg / 1.06 deg
- max FoV = 4.1 deg

Celestron:
- Do = 130mm
- f0 = 650mm  # f5

:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""

import numpy as _np


def arctan(arg):
    """Output in DEGREES"""
    return _np.arctan(arg) * 180 / _np.pi


def sl_f(dobj, dimg):
    """Simple Lens

    return the (Effective) Focal Length (FL) based on the distance of the
    Object `dobj` and the Image `dimg`. The unit of the FL is the same as
    `dobj` and `dimg`.
    """
    return (dobj * dimg) / (dobj + dimg)


def sl_s(efl, d):
    """Simple Lens

    return the distance of the object/image based on the distance of the
    image/object `d` and the (Effective) Focal Length `efl`. The output unit
    is the same of `d` and the (Effective) Focal Length `efl`.
    """
    return (d * efl) / (d - efl)


def sl_mag(dobj, dimg):
    """Simple Lens

    return the magnification of the image (no-unit). Distance of the
    Object `dobj` and the Image `dimg` must be the same unit.
    """
    return -dimg / dobj


def st_mag(fo, fe):
    """Simple Telescope

    return the magnification of a simple telescope (no-unit). The Objective
    Focal Distance `fo` and the Eyepiece Focal Distance `fe` must be the same
    unit.
    """
    return fo / fe


def st_s(fo, fe):
    """Simple Telescope

    In order to find where the eye lens forms the exit pupil e', use the
    simple lens equation substituting `z = -(fo + fe)`. Output unit the same as
    input.
    """
    return fo / (fe * (fo + fe))


def st_maxfov(De, fo):
    """Simple telescope

    Returns the maximum telescope FoV in degrees. `De` (eyepiece diameter) and
    `fo` telescope focal lenght (same units).

    1 rad is equivalent do 57.3 degree.
    """
    return De / fo * 57.3


def st_fov(fo, fe, afov):
    """Simple Telescope

    Return the telescope FoV in the same units as the `afov`. `afov` is the
    Apparent Field of View of the eyepiece. Usually a spec of the component.

    For more info., see ISO 14132-1:2002 standard.
    """
    return afov * fe / fo


def st_exitpupil(fe, foDo):
    """Simple Telescope

    The average healthy human pupil normally opens about 2mm in daylight, and
    7mm in the dark. Exit pupil with 2 mm diameter is close to be optimal
    (in average) for small fine nebulae like tiny galaxies.

    At the exit pupil size larger than about 2mm in diameter, eye aberrations
    begin to dominate diffraction effect, increasing progressively with the
    pupil size. Thus, telescopic resolution is aberrations-limited for exits
    pupils larger than ~2mm, and diffraction-limited for smaller pupils.

    returns the Exit Pupil diameter size. `fe` focal length of the eyepiece,
    and `foDo` is the telescope f#. Output unit is the same as `fe`.
    """
    return fe / foDo


def st_airydiskep(exitpupil):
    """Simple Telescope

    returns the Airy Disk size in arcsec for a telescope with `exitpupil`
    diameter given in mm.
    """
    return 4.6 * 60 / exitpupil


def st_airydisk(a):
    """Simple telescope

    If a was the radius of the aperture in inches, and the wavelength of
    light was assumed to be 0.000022 inches (560 nm; the mean of visible
    wavelengths).

    returns the Airy Disk size in arcsec for a telescope with aperture `a`
    given in **mm**.
    """
    return 2.76 * 25.4 / a


def ep_afov(fe, De):
    """Eyepiece -- Simple lens

    return the apparent FoV (AFOV) in **degrees**. `fe` and `De` same units.
    """
    return arctan(De / fe)


def st_ppi(Do, M):
    """Simple telescope

    Returns the power per inch (PPI). `fo` in mm.

    The lower the PPI, the fainter the observable targets. In average seeing
    conditions, about 30 PPI as a practical maximum.
    """
    return M * 25.4 / Do


def ep_info():
    """Focal lengths of eyepieces are very short and will typically range
    anywhere between 3mm and about 40mm.

    The mostly used eyepiece barrel diameters for telescopes are:
    - 0.965 in. (24.5 mm)
    - 1.25 in. (31.75 mm)
    - 2 in. (50.8 mm)
    - 2.7 in. (68.58 mm)
    """
    return


def binocular_info():
    """Orion Explorer 10x50 binoculars have a magnification of 10x and an
    objective lens 50mm wide.

    Many binoculars are marked with a second set of numbers, close to one of
    the eyepieces and often below the magnification and objective lens size.
    For example, a Celestron UpClose G2 7x35 binoculars also have the
    specifications of 9.2/483 ft/161m. In this case, the first number is the
    field of view, in degrees. The second and third numbers indicate that you
    can observe entire object of 483 ft/161m long that is placed 1000m away.

    arctan(161/1000) = 9.2 [deg]

    Another example: The Orion Explorer 10x50 binocular has as second set of
    numbers 'Field 6.5^o 114m/1000m'.

    arctan(114/1000) = 6.5 [deg]
    """
    return


def st_platescale(fo):
    """Simple telescope

    d(theta)/d(size) = 1[angle]/f[length]

    output in arcsec/[length unit of `fo`]
    """
    return 206265.0 / fo


def tl_f(r1, r2, d, n):
    """Thick Lens

    return the (Effective) Focal Length (FL) based on the radius of curvature
    of the first `r1` and second `r2` surfaces, the full lens thickness `d` and
    the refraction index `n`. Same units as `r1`, `r2` and `d`.
    """
    return ((n - 1) * (1 / r1 - 1 / r2 + (n - 1) * d / (n * r1 * r2))) ** -1.0


def dl_f(f1, f2, d):
    """Double lens system

    return the (Effective) Focal Length (FL) of a set of two lenses.  Same
    units as `f1`, `f2` and `d`.

    Remember that there is a lot of difference between the EFL and the BFL. The
    EFL is built over the Principal Point, or the point whose angle builds the
    final image. The BFL is simply a length distance (no direct image property
    associated).
    """
    return f1 * f2 / (f1 + f2 - d)


def dl_mag(f1, f2, d, so):
    """Double lens system

    return the magnification of a set of two lenses. Same length units for
    `f1`, `f2`, `d` and `so` (position of the Object **from the first lens
    backwards**).

    For an object at the infinity, see "Simple Telescope".
    """
    return f1 * f2 / ((so - f1) * (d - f2 - (so * f1) / (so - f1)))


def dl_simg(f1, f2, d, so):
    """Double lens system

    return the position of the image of set of two lenses. Same length units
    for `f1`, `f2`, `d` and `so` (position of the Image **after the second lens
    onwards**).

    For an Object at the infinity, see "Simple Telescope".
    """
    return f2 * (d - (so * f1) / (so - f1)) / (d - f2 - (so * f1) / (so - f1))


def dl_bfl(f1, f2, d):
    """Double lens system

    return the Back Focal Length (BFL) of a set of two lenses.  Same
    units as `f1` and `d`.

    See dl_f().
    """
    # return dl_f(f1, f2, d)*(f1-d)/f1
    return (f1 - d) * f2 / (f1 + f2 - d)


def st_observablemag(Do, a=7, mag_n=0):
    """Simple Telescope

    Returns the Observable magnitude in telescope. `Do` is the telescope
    aperture **in mm** and `a=7`mm is the estimated eye pupil diameter.
    """
    return _np.log(Do**2 / a**2) + mag_n
