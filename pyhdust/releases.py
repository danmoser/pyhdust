# -*- coding:utf-8 -*-

"""PyHdust auxiliary module: PyHdust releases control.

History
============
v1.007 @ 2016-09-16
---------------------
- New ploting 2d routines

v1.005 @ 2016-08-17
---------------------
- New working generic interpolation sed2/XDR

v1.002 @ 2016-06-26
----------------------
- New folder structure

v0.999 @ 2016-06-14
----------------------
- Corrections on setup.py

v0.998 @ 2016-06-13
----------------------
- Corrections on setup.py

v0.997 @ 2016-06-06
----------------------
- Corrections on the documentation

v0.996 @ 2016-06-06
----------------------
- Setup.py and minor corrections

v0.995 @ 2016-04-30
----------------------
- Corrections on intt

v0.994 @ 2016-04-2x
----------------------
- New intt functions!

v0.993 @ 2016-04-30
----------------------
- New inp.makeCSGrid_bistabWind1Dust()

v.0970 @ 2016-01-xx
-----------------------
- New README.md
- phc.cycles (new)
- roadmap.TODO (new)
- hdt.plot_obs (new)
- bcd (new module)

v0.967 @ 2015-12-30
----------------------
- PEP8 Standardization
- intt.I: Vieira+2015 models

v0.966 @ 2015-12-16
----------------------
- Correction on setup.py

v0.965 @ 2015-12-16
----------------------
- hdt.mergesed2: SED was the first because, if present, the code will check if 
    other bands parameters are the same (i.e., observers, Rstar, Rwind).
- hdt.mergesed2: The criteria I elected for distinguishing between broad-band 
    and line (Sobolev 0/1) is the presence of "_SEI" extension in the filename, 
    assumed that the line rest wavelength is the BAND CENTER WAVELENGTH. There 
    is an option to the user manually put it.
- hdt.mergesed2: A new output format of the numbers was done.
- inp.makeSourceGrid: Function created as discussed
- inp.makeDiskGrid: Define convSig2Rho=True, then all the values of sig0vals 
will be considered as rho0 values.

v0.964 @ 2015-09-15
----------------------
- hdrpil module added
- readdust function implemented
- Other mirror changes

v0.956 @ 2015-05-08
----------------------
- Contributions from Bednarski to poltools+fieldstars
- Corrections to work on Windows (binary files)
- Triangle module added

v0.955 @ 2015-05-15
----------------------
- spt.fitzpatrick included
- intt.img2fits and intt.data2fitscube rotation correction
- Documentation of poltools.py (Portuguese)
- hdt.readSingleBe improved
- inp.makeDiskGrid with SingleBe option
- intt.img2fits debugged
- inp.makeInpJob BlueGene support added

v0.954 @ 2015-05-08
----------------------
- releases.py do automatically the docs' rsync
- hdt.plottemp improved
- General: "== False" replaced by "not"
- spt.*Plot* corrections
- spt.kuruczflux correction

v0.953 @ 2015-05-01
----------------------
- phc.bindata, now yerr is an optional array
- spt,dtb improvements save/load Class
- spt.shiftfits improved
- spt.plotSpecData improved
- spt.cardelli included

v0.952b @ 2015-04-28
---------------------
- Create of releases.py
-   Automatic update of setup.py and documentation files.
- Defined flux units at SED2 file manipulation
- Implemented Kurucz flux unity correction
- Removed *ra2deg* variable
- Updated phc.rot_stars (Beta(W))
- Created XDR BeAtlas

Previously
------------
- Added options `force` and `chknames` to genStdLog and genObjLog


:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def setRelease():
    """ Read the version values from __init__.py and write it to the setup.py
    and doc files. """
    import os
    import re
    from pyhdust import __version__, hdtpath
    #  
    f0 = open(os.path.join(os.path.split(hdtpath()[:-1])[0], 'setup.py'))
    lines = f0.readlines()
    f0.close()
    i = [lines.index(x) for x in lines if x.find('version=') > -1][0]
    verline = lines[i]
    if len(verline.split("'")) == 3:
        oldver = verline.split("'")[1]
    else:
        oldver = verline.split('"')[1]
    lines[i] = lines[i].replace(oldver, str(__version__))
    f0 = open(os.path.join(os.path.split(hdtpath()[:-1])[0], 'setup.py'), 'w')
    f0.writelines(lines)
    f0.close()
    print('# ../setup.py file updated!')
    #
    f0 = open(os.path.join(os.path.split(hdtpath()[:-1])[0], 'docs', 
        'index.rst'))
    lines = f0.readlines()
    f0.close()
    i = [lines.index(x) for x in lines if x.find('at **version') > -1][0]
    verline = lines[i]
    oldver = re.findall('[-+]?[0-9]*\.?[0-9]*\.?[0-9]+', verline)[0]
    lines[i] = lines[i].replace(oldver, str(__version__))
    f0 = open(os.path.join(os.path.split(hdtpath()[:-1])[0], 'docs', 
        'index.rst'), 'w')
    f0.writelines(lines)
    f0.close()    
    print('# docs/index.rst file updated!')
    os.chdir(os.path.join(os.path.split(hdtpath()[:-1])[0], 'docs'))
    os.system('make html')
    # os.system('rsync -rP _build/html/ astroweb:/www/moser/www/doc')
    # os.system('rsync -rP _build/html/ /data/Dropbox/Public/doc')
    # print('# From version {0} to {1}'.format(oldver, __version__))
    # os.system(
    # 'firefox https://dl.dropboxusercontent.com/u/6569986/doc/index.html &')
    # os.chdir('../..')
    # os.system('git push')
    #  os.system('disown')
    return


# MAIN ###
if __name__ == "__main__":
    pass
