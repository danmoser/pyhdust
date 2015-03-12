# -*- coding:utf-8 -*-

"""
PyHdust *beatlas* module: BeAtlas specific variables and functions.

:license: GNU GPL v3.0 (https://github.com/danmoser/pyhdust/blob/master/LICENSE)
"""
import os as _os
import numpy as _np
from glob import glob as _glob
import pyhdust.phc as _phc

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"
__date__ ='1 March 2015'
__version__ = '0.9'


vrots = [[259.759,354.834,417.792,464.549,483.847],\
     [252.050,346.163,406.388,449.818,468.126],\
     [245.127,336.834,399.983,448.076,467.806],\
     [239.522,329.496,388.734,432.532,450.806],\
     [234.301,321.139,379.297,423.241,441.122],\
     [228.538,313.797,370.343,412.488,429.914],\
     [219.126,299.656,354.547,395.821,413.008],\
     [211.544,288.840,341.081,380.426,396.978],\
     [203.438,279.328,328.666,365.697,380.660],\
     [197.823,268.964,316.901,353.568,368.506],\
     [192.620,262.688,308.208,341.963,356.410],\
     [187.003,255.125,299.737,332.511,346.043]]

obs = [1.1,1.2,1.3,1.4,1.45]

ms = [14.6, 12.5, 10.8, 9.6, 8.6, 7.7, 6.4, 5.5, 4.8, 4.2, 3.8,3.4]

Ms = _np.array([14.6, 12.5, 10.8, 9.6, 8.6, 7.7, 6.4, 5.5, 4.8, 4.2, 3.8, 3.4],\
dtype=str)

Tp11 = _np.array([28905.8,26945.8,25085.2,23629.3,22296.1,20919.7,\
18739.3,17063.8,15587.7,14300.3,13329.9,12307.1])

sig0 = _np.logspace(_np.log10(0.02),_np.log10(4.0),7)

Sig0 = ['{0:.2f}'.format(x) for x in sig0]

def rmMods(modn, Ms, clusters=['job']):
    """
    Remove the *.inp models of models `modn` according to the list structure
    below.

    | Masses list ans sig0 POSITION do be excluded
    | Ms = [
    | ['14.6', [0]],
    | ['12.5', [0,-1]],
    | ['10.8', [0,-1]],
    | ['09.6',  [0,-2,-1]],
    | ['08.6',  [0,-2,-1]],
    | ['07.7',  [0,-2,-1]],
    | ['06.4',  [0,-3,-2,-1]],
    | ['05.5',  [0,-3,-2,-1]],
    | ['04.8',  [-4,-3,-2,-1]],
    | ['04.2',  [-4,-3,-2,-1]],
    | ['03.8',  [-4,-3,-2,-1]],
    | ['03.4',  [-4,-3,-2,-1]],]

    INPUT: string, structured list

    OUTPUT: *files removed
    """
    #Create sig0 list
    sig0s = Sig0
    project = _phc.trimpathname(_os.getcwd())[1]
    for cl in clusters:
        file = open('{0}s/{0}s_{1}_mod{2}.sh'.format(cl, project, modn))
        lines = file.readlines()
        file.close()
        for item in Ms:
            M = item[0]
            exsig = item[1]
            for rm in exsig:
                _os.system('rm mod{0}/mod{0}*_sig{1}*_M{2}*.inp'.format(modn,
                sig0s[rm], M))
                print('# Deleted mod{0}/mod{0}*_sig{1}*_M{2}*.inp'.format(modn,
                sig0s[rm],M))
                _os.system('rm {3}s/mod{0}*_sig{1}*_M{2}*.{3}'.format(modn,
                sig0s[rm], M, cl))
                lines = [line for line in lines if (line.find('_sig{0}'.format(
                sig0s[rm]))==-1 or line.find('_M{0}'.format(M))==-1)]
		file = open('{0}s/{0}s_{1}_mod{2}.sh'.format(cl, project, modn), 'w')
		file.writelines(lines)
		file.close()
    #End prog
    return 


def readBAsed(xdrpath):
    """ Read the BeAtlas SED release.
    """
    #~ A lot of complicated code
    intervmods = _np.array([ [3.4, 14.6], [0., 1.], [0.01, 2.], [3., 4.5] ])
    lbdmods = _np.arange(10)
    modelos = _np.zeros((5, 10))
    return intervmods, lbdmods, modelos 

def interpolBA(params, modelos):
    """ Interpola os `modelos` para os parametros `params` """
    return 1

def breakJob(n, file):
	""" Break the jobs/jobs_Project_modn.sh into n files 
	../jobs_Project_modn_##.txt to be used with `dispara` """
	f0 = open(file)
	lines = f0.readlines()
	f0.close()
	lines.sort()
	lines = [line.replace('qsub ','') for line in lines]
	outname = _phc.trimpathname(file)[1].replace('.sh','')
	N = len(lines)
	for i in range(n):
		f0 = open('{0}_{1:02d}.txt'.format(outname, i), 'w')
		f0.writelines(lines[i*N/n:(i+1)*N/n])
		f0.close()
	print('# {0} files created!'.format(n))
	return

def correltable(pos):
    """ Create the correlation table of Domiciano de Souza+ 2014. """
    nwalkers = len(pos)
    ndim = len(pos[0])
    fig = _plt.figure()
    for i in range(ndim**2):
        ax = fig.add_subplot(ndim,ndim,i+1)
        if i+1 in [ 1+x*(ndim+1) for x in range(ndim) ]:
            ax.hist(pos[:,i/ndim], 20)
        else:
            ax.plot(pos[:,i/ndim], pos[:,i%ndim], 'o', markersize=2)
    _plt.savefig('correl.png', Transparent=True)
    _plt.close()
    print('# Figure "correl.png" saved!')
    return

### MAIN ###
if __name__ == "__main__":
    pass
