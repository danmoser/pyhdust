# -*- coding:utf-8 -*-

"""
PyHdust auxiliary module: PyHdust releases control.

History
============
..
    v0.954 @ 2015-0
    ----------------------

v0.953 @ 2015-05-
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

To Be Done
============
General
---------
- pip install for pyhdust
- Transparent legend-box on the graphs
- Change scipy interpolate with numpy.interp (hdt.doFilterConv, others...)
- Define unique 'Planejamento_LNA' as MySQL (targets [coordinates], P and Ph0 [mags.], nights observed vs. bad weather, etc.)
- Standardization at CGS/SI units
- Finish help docs+translations to English
- \*Define `path` policy
- Replace if os.path.exists+os...mkdir with phc.outfld
- Remove *rmext*. Replace with os.* tool (Bebnarski)
- Remove *outfold* from spectools
- Check *makeInpJob* `scrid` variable
- Automatically Update of this list with TBD flag.

__init__
-----------
- update obs_calcs (spherical triangles)
- rotstars: add option for ellipsoidal star
- genlog: add check for *.sigma files
- merge: be ``smart'' and be independent of filters definitions (e.g., continuum or line based on the number of points). 

beatlas
----------
- Develop rotine correlations emcee (e.g., triangle.py)

input
-----------
- Implement Mdot11 option in makeDiskGrid

phc
---------
  
spec
-----------

:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def listFunctions():
    """ List all *PyHdust* functions. """ 
    import pyhdust as hdt
    print("# pyhdust")
    print("\n".join([x for x in dir(hdt) if x[0] != "_"]))
    
    import pyhdust.beatlas as bat
    print("# pyhdust.beatlas")
    print("\n".join([x for x in dir(bat) if x[0] != "_"]))
    
    import pyhdust.input as inp
    print("# pyhdust.input")
    print("\n".join([x for x in dir(inp) if x[0] != "_"]))
    
    import pyhdust.interftools as intt
    print("# pyhdust.interftools")
    print("\n".join([x for x in dir(intt) if x[0] != "_"]))
    
    import pyhdust.poltools as polt
    print("# pyhdust.poltools")
    print("\n".join([x for x in dir(polt) if x[0] != "_"]))
    
    import pyhdust.phc as phc
    print("# pyhdust.phc")
    print("\n".join([x for x in dir(phc) if x[0] != "_"]))
    
    import pyhdust.singscat as sst
    print("# pyhdust.singscat")
    print("\n".join([x for x in dir(sst) if x[0] != "_"]))
    
    import pyhdust.spectools as spt
    print("# pyhdust.spectools")
    print("\n".join([x for x in dir(spt) if x[0] != "_"]))

    return


def setRelease():
    """ Read the version values from __init__.py and write it to the setup.py
    and doc files. """
    import os
    import pyhdust.phc as phc
    from pyhdust import __version__, hdtpath
    #~ 
    f0 = open('{0}/setup.py'.format(hdtpath()))
    lines = f0.readlines()
    f0.close()
    i = [lines.index(x) for x in lines if x.find('version') > -1]
    i = i[0]
    oldver = phc.fltTxtOccur('version', [lines[i]], asstr=True)
    lines[i] = lines[i].replace(oldver, str(__version__))
    f0 = open('{0}/setup.py'.format(hdtpath()), 'w')
    f0.writelines(lines)
    f0.close()
    print('# ../setup.py file updated!')
    #~
    f0 = open('{0}/docs/index.rst'.format(hdtpath()))
    lines = f0.readlines()
    f0.close()
    i = [lines.index(x) for x in lines if x.find('at **version') > -1]
    i = i[0]
    oldver = phc.fltTxtOccur('version', [lines[i]], asstr=True)
    for i in range(len(lines)):
        if lines[i].find(oldver):
            lines[i] = lines[i].replace(oldver, str(__version__))
    f0 = open('{0}/docs/index.rst'.format(hdtpath()), 'w')
    f0.writelines(lines)
    f0.close()    
    print('# docs/index.rst file updated!')
    os.chdir('{0}/docs'.format(hdtpath()))
    os.system('make html')
    print('# From version {0} to {1}'.format(oldver, __version__))
    os.system('midori _build/html/index.html &')
    #~ os.system('disown')
    return


### MAIN ###
if __name__ == "__main__":
    pass
