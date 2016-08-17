# -*- coding:utf-8 -*-

"""PyHdust *input* module: Hdust input tools.

Definitions on Classes:
- fname: is the original file for reference (kept as stated in the input). 
The following variables are defined based on the fname.
- proj: directory to the project (mandatory)
- modn: 'modn' structure
- suf: only the suffix.
- one use Class.set_fname(fname) to update everything.
- convention is to use "_" as separador. So, use "-" for a given **name**, as 
for the projects.

:co-author: Rodrigo Vieira
:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function
import os as _os
import numpy as _np
from glob import glob as _glob
from itertools import product as _product
import pyhdust.phc as _phc
import pyhdust.rotstars as _rot
import pyhdust as _hdt
from collections import OrderedDict as _OrderedDict
from collections import Mapping as _mapdict
from six import string_types as _string_types
import warnings as _warn

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


class Input(object):
    """ Input class; no fname (ie., specific suf/model defined) 

    - Criado no novo `make_inpjob`. 
    - Terminando verificação de existência prévia de *.inp e *.sed2
    - Incluído verificação do *_SEI.sed2
    """
    f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF_disco.txt'))
    mod = f0.readlines()
    f0.close()

    attl_to_chk = ['case1', 'case2', 'composition', 'controls', 'gridcells', 
        'observers', 'perturbations', 'source']

    def __init__(self, proj=None, modn='mod01'):
        self.proj = proj
        if proj is None:
            self.proj = _os.path.split(_os.getcwd())[1]
        self.set_modn(modn)
        return

    def set_modn(self, modn):
        self.modn = modn
        if isinstance(modn, _string_types):
            if modn.isdigit():
                modn = int(modn)
            elif not modn.startswith('mod'):
                raise ValueError('# ERROR! modn must startswith "mod"')
        if isinstance(modn, ( int, long )):
            if modn > 0 and modn < 1000:
                self.modn = 'mod{0:02d}'.format(int(modn))
            else:
                raise ValueError('# ERROR! modn must be <1000')

    def set_input(self, docasesl=None, case1='step1', case2='step1_refine', 
        case3l=None, imagesl=None, composition='pureH', controls='controls', 
        gridcells='grid', observers='observers', perturbations='', source='', 
        chkout=True, c1max=20, c2max=24, c3prefl=None, clustl=None, ncore=48, 
        walltime='72:00:00', email='$USER@localhost', touch=False, **kwargs):

        if isinstance(case3l, _string_types):
            self.case3l = [case3l]
        elif case3l is None:
            self.case3l = []
        else:
            self.case3l = case3l
        if isinstance(imagesl, _string_types):
            self.imagesl = [imagesl]
        elif imagesl is None:
            self.imagesl = []
        else:
            self.imagesl = imagesl

        # docases == list of numbers    
        if docasesl is None:
            self.docasesl = [1, 2, 3]
        else:
            self.docasesl = docasesl

        self.case1 = case1
        self.case2 = case2
        self.composition = composition
        self.controls = controls
        self.gridcells = gridcells
        self.observers = observers 
        self.perturbations = perturbations
        self.source = source

        self.chkout = chkout
        self.c1max = str(c1max)
        self.c2max = str(c2max)
        if isinstance(c3prefl, _string_types):
            self.c3prefl = [c3prefl]
        elif c3prefl is None:
            self.c3prefl = []
        else:
            self.c3prefl = c3prefl
        if isinstance(clustl, _string_types):
            self.clustl = [clustl]
        elif clustl is None:
            self.clustl = []
        else:
            self.clustl = clustl
        self.ncore = ncore
        self.walltime = walltime
        self.email = email
        self.touch = touch

        for it in kwargs.items():
            setattr(self, *it)
        return

    def write_inp(self):
        bdir = self.proj
        if bdir == _os.path.split(_os.getcwd())[1]:
            bdir = ""

        head = "PROJECT = {0}\nMODEL = {1}".format(_os.path.split(
            self.proj)[1], self.modn[3:])

        c1 = "\n\n"
        if (_os.path.exists(_os.path.join(bdir, self.modn, self.modn + 
            self.suf + self.c1max + '.temp')) and self.chkout) or (1 not in 
            self.docasesl):
                c1 += '! '
        c1 += ("SUFFIX='{suf}'  SIMULATION='{case1}' SOURCE='{source}' "
            "COMPOSITION='{composition}' GRIDCELLS='{gridcells}' "
            "CONTROLS='{controls}' ".format(suf=self.suf, case1=self.case1,
                source=self.source, composition=self.composition, 
                gridcells=self.gridcells, controls=self.controls))
        if len(self.perturbations) > 0:
            c1 += "PERTURBATIONS='{0}' ".format(self.perturbations)

        c2 = '\n\n'
        if (_os.path.exists(_os.path.join(bdir, self.modn, self.modn + 
            self.suf + self.c2max + '.temp')) and self.chkout) or (2 not in 
            self.docasesl):
                c2 += '! '
        c2 += ("SUFFIX='{suf}'  SIMULATION='{case2}' SOURCE='{source}' "
            "COMPOSITION='{composition}' GRIDCELLS='{gridcells}' "
            "CONTROLS='{controls}' ".format(suf=self.suf, case2=self.case2,
                source=self.source, composition=self.composition, 
                gridcells=self.gridcells, controls=self.controls))
        if len(self.perturbations) > 0:
            c2 += "PERTURBATIONS='{0}' ".format(self.perturbations)

        c3 = '\n\n'
        gen = (sim for sim in self.case3l if 3 in self.docasesl)
        for sim in gen:
            if self.chkout:
                pref = self.c3prefl[self.case3l.index(sim)]
                if not pref.endswith('_'):
                    pref = pref+'_'
                s2name = _os.path.join(bdir, self.modn, pref + self.modn + 
                    self.suf)
                if (_os.path.exists(s2name + '.sed2')) or (
                    _os.path.exists(s2name + '_SEI.sed2')):
                    c3 += '! '
            c3 += ("SUFFIX='{suf}'  SIMULATION='{sim}' SOURCE='{source}' "
                "COMPOSITION='{composition}' GRIDCELLS='{gridcells}' "
                "CONTROLS='{controls}' OBSERVERS='{obs}' ".format(suf=self.suf, 
                    sim=sim, source=self.source, 
                    composition=self.composition, gridcells=self.gridcells, 
                    controls=self.controls, obs=self.observers))
            if len(self.perturbations) > 0:
                c3 += "PERTURBATIONS='{0}' ".format(self.perturbations)
            img = self.imagesl[self.case3l.index(sim)]
            if len(img) > 0:
                c3 += "IMAGES='{0}' ".format(img)
            c3 += '\n'

        outname = _os.path.join(bdir, self.modn, self.modn + self.suf + '.inp')
        f0 = open(outname, 'w')
        f0.writelines(head+c1+c2+c3)
        f0.close()
        print('# {0} written!'.format(outname))
        return 

    def _hassiminp(self, inpname):
        inp = open(inpname).read().split('\n')
        sims = [1 for line in inp if (line.upper().find('SUFFIX')>= 0 and 
            not line.startswith('!'))]
        return bool(sims)

    def write_job(self):
        if self.touch is True:
            raise NotImplementedError('# touch option not available')

        bdir = self.proj
        if bdir == _os.path.split(_os.getcwd())[1]:
            bdir = ""
        outname = _os.path.join(bdir, self.modn, self.modn + self.suf + '.inp')
        intname = _os.path.join(self.proj, self.modn, 
            self.modn + self.suf + '.inp')

        if not self._hassiminp(outname):
            print('# Message: No simulations inside {0}'.format(outname))
            return

        for cl in self.clustl:
            if cl.lower().startswith('job'):
                cl = 'job'
                reffile = 'REF.job'
                reps = ( (3, 'hdust', 'hd_'+'_'.join((self.proj, self.modn))), 
                    (4, '128', self.ncore), (4, '36:00:00', self.walltime), 
                    (8, 'alexcarciofi@gmail.com', self.email), 
                    (11, 'hdust_bestar2.02.inp', intname),
                    (23, 'cd ${PBS_O_WORKDIR}', ('cd ${{PBS_O_WORKDIR}}\n'
                        'printf "{0}\\n\\n" >> output_${{PBS_JOBID}}\n'.format(
                            intname))),
                )
            elif cl.lower().startswith('oar'):
                cl = 'oar'
                reffile = 'REF.oar'
                reps = ( (1, 'hdust_dmf', 'hd_'+'_'.join((self.proj, 
                    self.modn))), 
                    (2, '12', int(round(self.ncore/6.))), 
                    (2, '24:00:00', self.walltime), 
                    (10, 'hdust_bestar2.02.inp', '{0}\n\nprintf "{0}\\n\\n"'
                        ' >> $OAR_JOB_ID'.format(intname)),
                )
            else:
                _warn.warn('# Warning! Option {0} for clusters ignored'.format(
                    cl))
                continue

            if not _os.path.exists(_os.path.join(bdir, cl)):
                _os.makedirs(_os.path.join(bdir, cl))
            clname = (_os.path.join(bdir, cl, self.modn + self.suf)+"."+cl)

            f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', reffile))
            modi = f0.readlines()
            f0.close()
            for r in reps:
                _phc.repl_fline_val(modi, *r)

            f0 = open(clname, 'w')
            print('# {0} written!'.format(clname))
            f0.writelines(modi)
            f0.close()
        return


class Disk(Input):
    """ Disk class """
    f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF_disco.txt'))
    mod = f0.readlines()
    f0.close()

    def __init__(self, proj=None, modn='mod01'):
        super(Disk, self).__init__(proj=proj, modn=modn)
        return

    def set_disk(self, renv=18.6, mh=1.5, ht=60., nr=2., dval=1e12, hseq=False, 
        alpha=.5, mu=.5, R0r=100, denstype=None, **kwargs):
        """ denstype = ['n0', 'sig0', 'mdot'] 

        ``sig0`` means g/cm2 instead of ``n0``"""
        self.renv = renv
        self.mh = mh
        self.ht = ht
        self.nr = nr
        self.dval = dval
        self.hseq = hseq
        self.alpha = alpha
        self.mu = mu
        self.R0r = R0r
        if isinstance(denstype, _string_types):
            if denstype.lower() in ['n0', 'sig0', 'mdot']:
                self.denstype = denstype
            denstype = None
        if denstype is None:
            _warn.warn('# Warning! Guessing the denstype...')
            self.denstype = 'n0'
            if dval < 1e-2:
                self.denstype = 'mdot'
            elif dval < 1e2:
                self.denstype = 'sig0'

        for it in kwargs.items():
            setattr(self, *it)
        return

    def set_suf(self, vfmt=''):
        """    vfmt = dict(zip(['renv', 'dval', 'ht', 'M', "n"], 
        ['{:04.1f}', '{:.1e}', '{:02.0f}', '{:04.1f}', '{:.1f}'])) """
        if isinstance(vfmt, _mapdict):
            self.vfmt = vfmt
        if hasattr(self, 'vfmt'):
            if isinstance(self.vfmt, _mapdict):    
                self.suf = '_' + '_'.join(it[0]+it[1].format(getattr(self, 
                    it[0])) for it in self.vfmt.items())
            else:
                delattr(self, vfmt)
                self.suf = vfmt
        else:
            self.suf = vfmt

    def __str__(self):
        return self.suf

    def makedisk(self):
        if not self.modn.startswith('mod'):
            raise ValueError('# ERROR! Invalid Disk.modn')
        modi = Disk.mod[:]

        if self.hseq is True and self.denstype is not 'mdot':
            print('# Warning! Your choice of HSEQ don\'t appears to be '
                'consistent')

        _ = _phc.repl_fline_val(modi, 13, '18.6', self.renv)
        if self.denstype != 'mdot':
            _ = _phc.repl_fline_val(modi, 20, '2.0', self.nr)
        else:
            _ = _phc.repl_fline_val(modi, 23, '1.', self.alpha)
            _ = _phc.repl_fline_val(modi, 24, '0.', self.R0r)
        if self.hseq:
            _ = _phc.repl_fline_val(modi, 25, '= 0', '= 1')
            _ = _phc.repl_fline_val(modi, 31, '0', '1')
            _ = _phc.repl_fline_val(modi, 36, '1.5', self.mh)
        else:
            _ = _phc.repl_fline_val(modi, 33, '1.5', self.mh)
        if float(self.ht) > 1000:
            _ = _phc.repl_fline_val(modi, 40, '18000.', self.ht)
        else:
            _ = _phc.repl_fline_val(modi, 38, '1', '2')
            if self.ht < 1:
                self.ht *= 100
            _ = _phc.repl_fline_val(modi, 43, '72.', self.ht)

        if self.denstype == 'sig0':
            raise ValueError('# ERROR! `sig0` is not a valid option yet')
        elif self.denstype == 'n0':
            _ = _phc.repl_fline_val(modi, 52, '2.35E13', self.dval)
        else:
            _ = _phc.repl_fline_val(modi, 49, '2', '3')
            _ = _phc.repl_fline_val(modi, 55, '1.E-9', self.dval)

        bdir = self.proj
        if bdir == _os.path.split(_os.getcwd())[1]:
            bdir = ""
        if not _os.path.exists(_os.path.join(bdir, self.modn)):
            _os.makedirs(_os.path.join(bdir, self.modn))
        self.fname = (_os.path.join(bdir, self.modn, self.modn + self.suf) + 
            ".txt")
        f0 = open(self.fname, 'w')
        print('# {0} written!'.format(self.fname))
        f0.writelines(modi)
        f0.close()

        return


class HdustMod(object):
    """ HdustMod doc 
    """
    def __init__(self, fname=None):
        self.fname = fname
        if fname is not None:
            self.set_fname(fname)
        return

    def set_fname(self, fname):
        self.fname = fname
        self.modn, self.suf = _os.path.split(fname)
        self.suf = _os.path.splitext(self.suf)[0]
        if (not self.modn.startswith('mod')) or (
            not self.modn.startswith('fullsed')):
            self.proj, self.modn = _os.path.split(self.modn)
        else:
            self.proj = ''
        if not self.modn.startswith('mod'):
            self.modn = ''
        if len(self.modn) == 0:
            self.modn = 'mod' + self.suf.split('mod')[-1].split('_')[0]
        self.suf = '_'+'_'.join(self.suf.split('mod')[-1].split('_')[1:])
        if hasattr(self, 'vdict'):
            if _os.path.splitext(fname)[1] == '.log':
                self._fill()
            else:
                self._fillfname()
        else:
            _warn.warn('Warning! No `vdict` found for HdustMod')

    def _fill(self):
        f0 = open(self.fname).read().split('\n')
        for it in self.vdict.items():
            setattr(self, it[0], _phc.fltTxtOccur(it[1], f0))
        return

    def _fillfname(self):
        for it in self.vdict.keys():
            setattr(self, it, _phc.fltTxtOccur(it, self.suf))
        return

    def get_nob(self):
        self.nob = int(_hdt.sed2info(self.fname)[1])
        return self.nob

    def get_lum(self):
        if _os.path.splitext(self.fname)[1] == '.log':
            logname = self.name
        else:
            bdir = self.proj
            if bdir == _os.path.split(_os.getcwd())[1]:
                bdir = ""
            logname = _os.path.join(bdir, self.modn, self.modn + self.suf + 
                '.log')
            if not _os.path.exists(logname):
                loglist = _glob(_os.path.join(bdir, self.modn, '*'+self.modn + 
                    self.suf + '*.log'))
                logname = [l for l in loglist if _os.path.getsize(l) > 500]
                if len(logname) == 0:
                    raise LookupError('Not found an equivalent to {0}'.format(
                        logname))
                logname = logname[0]
        #
        self.lum = _phc.fltTxtOccur('L =', open(logname).read().split('\n'))
        self.vrot = _phc.fltTxtOccur('Vrot ', open(logname).read().split('\n'))
        self.log = logname
        return self.lum

    def readfs2ob(self, pol=False):
        if not hasattr(self, 'nob'):
            self.get_nob()
        col = 3
        if pol:
            col = 7
        s2d = _hdt.readfullsed2(self.fname)
        self.arr_lbd = s2d[0, :, 2]
        self.arr_flx = s2d[:, :, col]
        self.obdeg = _np.round(_np.arccos(s2d[:, 0, 0])*180./_np.pi, 1)
        self.obs = s2d[:, 0, 0]
        return

    def __str__(self):
        if self.fname is None:
            return 'None'
        return '_'+'_'.join(it[0]+it[1].format(getattr(self, it[0])) for it in 
            self.vfmt.items())


class AeriMod(HdustMod):    
    """docstring for AeriMod"""
    vdict = _OrderedDict(zip(['dval', 'ht', "nr", 'renv', 'M'],
            ['n_0', 'Fraction', ' n ', 'R_env', ' M '])) 
    vfmt = _OrderedDict(zip(['dval', 'ht', "nr", 'renv', 'M'], 
        ['{:.1e}', '{:02.0f}', '{:.1f}', '{:04.1f}', '{:04.1f}']))

    def __init__(self, fname):
        super(AeriMod, self).__init__(fname)


class PhotMod(HdustMod):    
    """docstring for AeriMod"""
    vdict = _OrderedDict(zip(['W', 'M'], ['w =', ' M '])) 
    vfmt = _OrderedDict(zip(['W', 'M'], ['{:.3f}', '{:04.1f}']))

    def __init__(self, fname):
        super(PhotMod, self).__init__(fname)


def makeDiskGrid(modn='01', mhvals=[1.5], hvals=[60.], rdvals=[18.6], 
    mvals=[2.0], sig0vals=[1.], convSig2Rho=False, doFVDD=False, sBdays=None, 
    sBfiles=None, selsources='*', alpha=.5, mu=.5, R0r=300, Mdot11=False, 
    path=None):
    """
    | ###CONFIG. OPTIONS
    | #MODEL NUMBER
    | modn        = '02'
    | #The following filter will be applied to the SOURCE selection (string 
    | # fmt)
    | selsources = '*'
    |                              
    | #SUPERFICIAL DENSITY PROFILE EXPONENT
    | mvals       = [1.5,2.0,2.5,3.0]
    | #VERTICAL DENSITY PROFILE EXPONENT     
    | mhvals      = [1.5]
    | #FRACTION OF TEFF OF PRIMARY STAR
    | #This parameter sets if it you be FIXED to OB=1.1 case
    | hvals       = [72.]
    | #DISK RADIUS EQUATORIAL...    
    | rdvals      = [30.]
    | #SIGMA_0 VALUES
    | sig0vals    = _np.logspace(_np.log10(0.02),_np.log10(4.0),7)
    | 
    | #Do the Full VDD model for the corresponding sig0?
    | doFVDD = True 
    | alpha = 0.5
    | mu = 0.5
    | #WARNING: it only generates a single R0 value per batch. If you want to 
    | # change it, run it twice (or more)
    | R0r = 300
    | ###END CONFIG.
    |
    | Mdot11 = TODO
    |
    | WARNING: if convSig2Rho=True, the routine assumes ``sig0vals`` contains
    |  rho0vals. For that, (sBfiles==None or sBdays==None) and doFVDD==False.
    """
    # Consistency check
    if convSig2Rho and not ((sBfiles is None or sBdays is None) and doFVDD is 
        False):
        print('# ERROR! There is an inconsistency in the input.')
        print('# Check the *convSig2Rho*!')
    # 
    G = _phc.G.cgs
    Msun = _phc.Msun.cgs
    Rsun = _phc.Rsun.cgs
    kB = _phc.kB.cgs
    mH = _phc.mH.cgs
    yr = _phc.yr.cgs

    def doPL(prodI):
        '''
        Given a prodI (i.e., src,sig0,rd,h,m,mh), generates the Power-Law model 
        input
        '''
        src, sig0, rd, h, m, mh = prodI
        M, Req, Tp = _rot.readscr(src)
        Mstr = str(M)
        M *= Msun
        Req *= Rsun

        Th = h * Tp / 100.
        # a0 = (kB*h/100.*Tp/mu/mH)**.5
        a = (kB * Th / mu / mH)**.5
        n0 = (G * M / 2. / _np.pi)**.5 * sig0 / mu / mH / a / Req**1.5
        # Th = a**2*mu*mH/kB

        srcname = src.replace('source/', '').replace('.txt', '')

        wmod = mod[:]
        wmod[13] = wmod[13].replace('18.6', ('%.2f' % rd))
        wmod[20] = wmod[20].replace('2.0', ('%.2f' % m))
        wmod[33] = wmod[33].replace('1.5', ('%.2f' % mh))
        wmod[40] = wmod[40].replace('18000.', ('%.1f' % Th))      
        if convSig2Rho:
            wmod[52] = wmod[52].replace('2.35E13', ('%.2e' % sig0))    
            suffix = '_PLn{0:.1f}_rho{1:5.2e}_h{2:03.0f}_Rd{3:05.1f}_{4}'.\
                format((m + mh), sig0, h, rd, srcname)  
        else:
            wmod[52] = wmod[52].replace('2.35E13', ('%.2e' % n0))    
            suffix = '_PLn{0:.1f}_sig{1:.2f}_h{2:03.0f}_Rd{3:05.1f}_{4}'.\
                format((m + mh), sig0, h, rd, srcname)  

        f0 = open('mod' + modn + '/mod' + modn + suffix + '.txt', 'w')
        f0.writelines(wmod)
        f0.close()
        return

    def doMdot(prodI):
        '''
        Given a prodI (i.e., src,sig0,rd,h,m,mh), generates the full VDD model 
        input
        '''
        src, sig0, rd, h, m, mh = prodI
        M, Req, Tp = _rot.readscr(src)
        Mstr = str(M)
        M *= Msun
        Req *= Rsun

        Th = h * Tp / 100.
        a = (kB * Th / mu / mH)**.5
        # a0 = (kB*h/100*Tp/mu/mH)**.5
        # a = a0*Req0*Req**.25/Req/Req**.25

        R0 = R0r * Req
        Mdot = sig0 * Req**2 * 3 * _np.pi * alpha * \
            a**2 / (G * M * R0)**.5  # SI units
        Mdot = Mdot / Msun * yr
        # Th = a**2*mu*mH/kB

        srcname = src.replace('source/', '').replace('.txt', '')
        # suffix = '_NI_Mdot{:.1e}_Rd{:.1f}_R0{:.1f}_alp{:.1f}_h{:.1f}_{}'.\
        # format(Mdot,rd,R0/Req,alpha,h,srcname)
        suffix = '_NIa{0:.1f}_sig{1:.2f}_h{2:03.0f}_Rd{3:05.1f}_{4}'.format(
            alpha, sig0, h, rd, srcname)

        wmod = mod[:]
        wmod[13] = wmod[13].replace('18.6', ('%.2f' % rd))
        wmod[18] = wmod[18].replace('1', ('%d' % 2))
        wmod[23] = wmod[23].replace('1.', ('%.2f' % alpha)) 
        wmod[24] = wmod[24].replace('= 0.', ('= %.2f' % (R0 / Req)))
        wmod[25] = wmod[25].replace('= 0', ('= %d' % 1))
        wmod[31] = wmod[31].replace('0', ('%d' % 1))  
        wmod[40] = wmod[40].replace('18000.', '{0:.1f}'.format(Th)) 
        wmod[49] = wmod[49].replace('2', ('%d' % 3))            
        wmod[55] = wmod[55].replace('1.E-9', ('%.2e' % Mdot))

        f0 = open('mod' + modn + '/mod' + modn + suffix + '.txt', 'w')
        f0.writelines(wmod)
        f0.close()
        return

    def doSB(prodI, hseq=False):
        '''
        Given a prodI (i.e., sources,rdvals,hvals,mhvals,sBdays,sBfiles),
        generates the Single Be based model input
        '''

        src, rd, h, mh, day, sfile = prodI
        M, Req, Tp = _rot.readscr(src)
        Mstr = str(M)
        M *= Msun
        Req *= Rsun

        Th = h * Tp / 100.
        # a0 = (kB*h/100.*Tp/mu/mH)**.5
        a = (kB * Th / mu / mH)**.5
        # n0 = (G*M/2./_np.pi)**.5*sig0/mu/mH/a/Req**1.5
        # Th = a**2*mu*mH/kB

        srcname = src.replace('source/', '').replace('.txt', '')

        wmod = mod[:]
        wmod[13] = wmod[13].replace('18.6', ('%.2f' % rd))
        wmod[18] = wmod[18].replace('= 1', '= 4')
        wmod[28] = wmod[28].replace(
            'deltasco/Atsuo/1D/data/dSco_a035_01', (sfile))
        wmod[29] = wmod[29].replace('2.3', ('%.2f' % (day / 365.25)))
        if not hseq:
            wmod[33] = wmod[33].replace('1.5', ('%.2f' % mh))
            suffix = '_SB{0}_{1:.1f}d_h{2:03.0f}_Rd{3:05.1f}_{4}'.format(
                _phc.trimpathname(sfile)[1], day, h, rd, srcname)
        else:
            wmod[31] = wmod[31].replace('= 0', '= 1')
            wmod[36] = wmod[36].replace('1.5', ('%.2f' % mh))
            suffix = '_SB{0}_{1:.1f}d_hseq_Rd{2:05.1f}_{3}'.format(
                _phc.trimpathname(sfile)[1], day, rd, srcname)

        wmod[40] = wmod[40].replace('18000.', ('%.1f' % Th))      

        f0 = open('mod' + modn + '/mod' + modn + suffix + '.txt', 'w')
        f0.writelines(wmod)
        f0.close()
        return

    # TODO Setup Tpole = REF of a (scale height)    
    # Tps = dict(zip(Ms, Tp11))

    # PROGRAM BEGINS
    path0 = _os.getcwd()
    if path is not None:
        _os.chdir(path)
        if path[-1] != '/':
            path += '/'
    else:
        path = ''
    # Check modN folder
    if not _os.path.exists('mod{0}'.format(modn)):
        _os.system('mkdir mod{0}'.format(modn))

    # Select sources
    sources = _glob('source/' + selsources)

    # Load disk model
    f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF_disco.txt'))
    mod = f0.readlines()
    f0.close()

    if sBdays is None or sBfiles is None:
        for prodI in _product(sources, sig0vals, rdvals, hvals, mvals, mhvals):
            doPL(prodI)
            i = 0
            if doFVDD:
                i = 1
                doMdot(prodI)
        print('# {0:.0f} arquivos foram gerados !!'.format(len(sources) *
            len(sig0vals) * len(rdvals) * len(hvals) * (len(mvals) + i) * 
            len(mhvals)))
    else:
        for prodI in _product(sources, rdvals, hvals, mhvals, sBdays, sBfiles):
            doSB(prodI)
            i = 0
            if doFVDD:
                i = 1
                doSB(prodI, hseq=True)
        print('# {0:.0f} arquivos foram gerados !!'.format(len(sources) *
            len(rdvals) * len(hvals) * len(sBdays) * (len(mhvals) + i) * 
            len(sBfiles)))

    if path is not '':
        _os.chdir(path0)
    # END PROGRAM    
    return


def makeInpJob(modn='01', nodes=512, simulations=['SED'],
    docases=[1, 3], sim1='step1', sim2='step1_ref', composition='pureH',
    controls='controls', gridcells='grid', observers='observers',
    images=[''], perturbations=None, clusters=['job'], srcid='',
    walltime='24:00:00', wcheck=False, email='$USER@localhost', chkout=False,
    st1max=20, st1refmax=24, ctrM=False, touch=False, srcNf=None, path=None):
    """
    Create INP+JOB files to run Hdust.

    All SOURCE files must initiate by "Be_". Otherwise, the `makeInpJob` will
    not work. This is to satisfies the criteria of a specific disk model for
    each source star.

        ### Start edit here ###
        modn = '02'

        #clusters config
        # job = AlphaCrucis; oar = MesoCentre Licallo; ge = MesoCentre FRIPP
        clusters = ['job','oar','ge','bgq']
        clusters = ['oar']
        nodes    = 48
        #if wcheck == True, walltime will be AUTOMATICALLY estimated
        walltime = '3:00:00'
        wcheck   = True
        email    = 'user@gmail.com'

        #Check if the outputs already exist
        chkout = True
        #Above the values below, the step1 will be considered done!
        st1max = 26
        st1refmax = 30
        #Gera inp+job so' para o source com '1.45' no nome
        #Nao funciona caracteres especiais como * ou ?
        srcid       = '1.45'
        srcid       = ''
        #Se um dos 3 casos nao estiver presente, ele gera input comentado.
        docases = [1,2,3]
        #1 = step1  <> Gera inp+job so' para mod#/mod#.txt (SEM source, so disco)
        #habilita ADDSUFFIX; retira OBSERVERS e IMAGES
        sim1 = 'step1'
        #2 = step1_refine 
        sim2 = 'step1_refine'
        #3 = outros <> Gera inp+job so' para mod#/mod#SOURCE.txt (post-proc.)
        #retira ADDSUFFIX; adiciona OBSERVERS (e talvez IMAGES)
        simulations = ['sed','h','brg','halpha','uv']
        simulations = ['sed_sig','brg_M','halpha_M','uv','j','h','k','l','m','n',
         'q1','q2']
        simulations = ['SED','Ha']
        images      = ['','h','brg','halpha','uv']
        images        = simulations[:]
        composition = 'pureH'
        controls    = 'no_op'
        controls    = 'controls'
        ctrM        = False
        gridcells   = 'grid'
        observers   = 'obs'
        touch       = True
        ###stop edition here
    """
    def isFloat(x):
        try:
            a = float(x)
        except ValueError:
            return False
        else:
            return True

    def doCase1(inp, cases):
        case1 = inp[:]
        case1[0] = case1[0].replace('suffix', suf)
        case1[1] = case1[1].replace('pureH', composition)
        if ctrM:
            i = suf.find('_M')
            M = suf[i:i + 7]
            case1[2] = case1[2].replace('controls', controls + M)
        else:
            case1[2] = case1[2].replace('controls', controls)
        case1[3] = case1[3].replace('grid', gridcells)
        case1[4] = case1[4].replace('step1', sim1)
        case1[5] = case1[5].replace('source', src)
        if perturbations is not None:
            case1.append("PERTURBATIONS = '{0}'\n".format(perturbations))
        if 1 not in cases:
            for i in range(len(case1)):
                case1[i] = '! ' + case1[i]
        return case1

    def doCase2(inp, cases):
        case1 = inp[:]
        case1[0] = case1[0].replace('suffix', suf)
        case1[1] = case1[1].replace('pureH', composition)
        if ctrM:
            i = suf.find('_M')
            M = suf[i:i + 7]
            case1[2] = case1[2].replace('controls', controls + M)
        else:
            case1[2] = case1[2].replace('controls', controls)
        case1[3] = case1[3].replace('grid', gridcells)
        case1[4] = case1[4].replace('step1', sim2)
        case1[5] = case1[5].replace('source', src)
        if perturbations is not None:
            case1.append("PERTURBATIONS = '{0}'\n".format(perturbations))
        if 2 not in cases:
            for i in range(len(case1)):
                case1[i] = '! ' + case1[i]
        return case1

    def doCase3(inp, simchk):
        case3 = []
        for i in range(len(simulations)):
            case1 = inp[:]
            case1[0] = case1[0].replace('suffix', suf)
            case1[1] = case1[1].replace('pureH', composition)
            if ctrM:
                j = suf.find('_M')
                M = suf[j:j + 7]
                case1[2] = case1[2].replace('controls', controls + M)
            else:
                case1[2] = case1[2].replace('controls', controls)
            case1[3] = case1[3].replace('grid', gridcells)
            case1[5] = case1[5].replace('source', src)
            if simulations[i] == 'SED':
                sig = suf[suf.find('_sig') + 4:suf.find('_sig') + 8]
                if isFloat(sig) and srcNf[i]:
                    case1[4] = case1[4].replace(
                        'step1', 'SED_sig{0}'.format(sig))
                else:
                    case1[4] = case1[4].replace('step1', simulations[i])
            elif srcNf[i]:
                case1[4] = case1[4].replace('step1', '{0}_{1}'.format(
                    simulations[i], src))
            else:
                case1[4] = case1[4].replace('step1', simulations[i])
            case1.append("OBSERVERS   = '{0}'\n".format(observers))
            if images[i] != '':
                case1.append("IMAGES      = '{0}'\n".format(images[i]))
            if perturbations is not None:
                case1.append("PERTURBATIONS = '{0}'\n".format(perturbations))
            case1.append('\n')
            if not simchk[i]:
                for i in range(len(case1)):
                    case1[i] = '! ' + case1[i]
            case3 += case1
        return case3

    def doJobs(mod, sel, nodes, addtouch='\n'):
        # load Ref
        f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF.{0}'.format(sel)))
        wout = f0.readlines()
        f0.close()

        outname = mod[mod.find('/') + 1:].replace('txt', sel)

        f0 = open('{0}s/{0}s_{1}_mod{2}.sh'.format(sel, proj, modn), 'a')
        if sel == 'job':
            wout[4] = wout[4].replace('128', '{0}'.format(nodes))
            wout[4] = wout[4].replace('36:00:00', '{0}'.format(walltime))
            wout[8] = wout[8].replace(
                'alexcarciofi@gmail.com', '{0}'.format(email))
            wout[11] = wout[11].replace('hdust_bestar2.02.inp', '{0}/{1}'.
            format(proj, mod.replace('.txt', '.inp')))
            if touch:
                wout[24] = addtouch
                modchmod = _phc.trimpathname(mod)
                modchmod[1] = modchmod[1].replace('.txt', '*')
                wout[31] = ('chmod -f 664 {0}/{1}/*{2}\nchmod -f 664 log/*\n' 
                    'chmod -f 664 ../../tmp/*\n'.format(proj, *modchmod))
            f0.writelines('qsub {0}/{1}s/{2}\n'.format(proj, sel, outname))
        elif sel == 'oar':
            wout[2] = wout[2].replace('12', '{0}'.format(int(round(
                nodes/6.))))
            wout[2] = wout[2].replace('24:0:0', '{0}'.format(walltime))
            wout[10] = wout[10].replace('hdust_bestar2.02.inp', '{0}/{1}'.
            format(proj, mod.replace('.txt', '.inp')))
            f0.writelines(
                'chmod -f a+x {0}/{1}s/{2}\n'.format(proj, sel, outname))
            f0.writelines(
                'oarsub -S ./{0}/{1}s/{2}\n'.format(proj, sel, outname))      
        elif sel == 'ge':
            wout[3] = wout[3].replace('48', '{0}'.format(nodes))
            wout[4] = wout[4].replace('45:00:00', '{0}'.format(walltime))
            wout[7] = wout[7].replace('dmfaes@gmail.com', '{0}'.format(email))
            wout[11] = wout[11].replace('hdust_bestar2.02.inp', '{0}/{1}'.
            format(proj, mod.replace('.txt', '.inp')))
            f0.writelines(
                'qsub -P hdust {0}/{1}s/{2}\n'.format(proj, sel, outname))
        elif sel == 'bgq':
            wout[14] = wout[14].replace('512', '{0}'.format(nodes))
            nodes = int(nodes)
            if nodes % 512 != 0:
                nrsv = (nodes // 512 + 1) * 128
            else:
                nrsv = (nodes // 512) * 128
            wout[10] = wout[10].replace('128', '{0}'.format(nrsv))
            wout[5] = wout[5].replace('24:00:00', '{0}'.format(walltime))
            wout[14] = wout[14].replace('hdust_bestar2.02.inp', '{0}/{1}'.
            format(proj, mod.replace('.txt', '.inp')))
            f0.writelines('chmod -f +x {0}/{1}s/{2}\n'.format(proj, sel, 
                outname))
            f0.writelines(
                'llsubmit ./{0}/{1}s/{2}\n'.format(proj, sel, outname))
        f0.close()

        f0 = open('{0}s/{1}'.format(sel, outname), 'w')
        f0.writelines(wout)
        print('# Saved: {0}s/{1}'.format(sel, outname))
        f0.close()
        return

    # PROGRAM START
    if srcNf is None:
        srcNf = len(simulations) * [False]
    path0 = _os.getcwd()
    if path is not None:
        _os.chdir(path)
        if path[-1] != '/':
            path += '/'
    else:
        path = ''
    # obtain the actual directory
    proj = _os.getcwd() 
    proj = proj[proj.rfind('/') + 1:]

    # Folder's checks
    for sel in clusters:
        if not _os.path.exists('{0}s'.format(sel)):
            _os.system('mkdir {0}s'.format(sel))
        elif _os.path.exists('{0}s/{0}s_{1}_mod{2}.sh'.format(sel, proj, 
            modn)):
            _os.system('rm {0}s/{0}s_{1}_mod{2}.sh'.format(sel, proj, modn))

    # list of mods
    mods = _glob('mod{0}/mod{0}*.txt'.format(modn))

    # load REF_inp
    f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF_inp.txt'))
    inp = f0.readlines()
    f0.close()

    for mod in mods:
        # Write inps
        f0 = open(mod.replace('.txt', '.inp'), 'w')
        f0.writelines('PROJECT = {0}\nMODEL = {1}\n\n'.format(proj, modn))

        suf = mod[mod.find('_'):-4]
        if mod.find('Be_') > -1:
            src = mod[mod.find('Be_'):-4]
        else:
            src = _os.path.splitext(_os.path.basename(srcid))[0]
        # src = 'aeri_draft'
        # if src.find(srcid) == -1:
        #     print('# Script skipped by missing "Be" in source...')
        #     continue

        cases = docases[:]
        # Do the touch thing
        addtouch = '\n'
        addtouch += 'chmod -f 664 ../../tmp/*\nchmod -f 664 {0}/mod{1}/*\n'.\
            format(proj, modn)
        if touch and ( (1 in cases) or (2 in cases) ):
            addtouch += 'touch {0}/{1}\n'.format(proj,
                                                 mod.replace('.txt', '.log'))
        if touch and 3 in cases:
            for sim in simulations:
                addtouch += 'touch {0}/{1}\n'.format(proj, mod.replace('.txt', 
                    '.chk')).replace('mod{0}/mod{0}'.format(modn), 
                    'mod{0}/{1}_mod{0}'.format(modn, sim))
                addtouch += 'touch {0}/{1}\n'.format(proj, mod.replace('.txt', 
                    '.err')).replace('mod{0}/mod{0}'.format(modn), 
                    'mod{0}/{1}_mod{0}'.format(modn, sim))
                addtouch += 'touch {0}/{1}\n'.format(proj, mod.replace('.txt',
                    '.log')).replace('mod{0}/mod{0}'.format(modn), 
                    'mod{0}/{1}_mod{0}'.format(modn, sim))
                addtouch += 'touch {0}/{1}\n'.format(proj, mod.replace('.txt',
                    '_SEI.chk')).replace('mod{0}/mod{0}'.format(modn), 
                    'mod{0}/{1}_mod{0}'.format(modn, sim))
                addtouch += 'touch {0}/{1}\n'.format(proj, mod.replace('.txt',
                    '_SEI.err')).replace('mod{0}/mod{0}'.format(modn), 
                    'mod{0}/{1}_mod{0}'.format(modn, sim))
                addtouch += 'touch {0}/{1}\n'.format(proj, mod.replace('.txt', 
                    '_SEI.log')).replace('mod{0}/mod{0}'.format(modn), 
                    'mod{0}/{1}_mod{0}'.format(modn, sim))
                err90a = '{0}/{1}'.format(proj, mod.replace('.txt', '.err').
                    replace('mod{0}/mod{0}'.format(modn), 'mod{0}/{1}_mod{0}'.
                    format(modn, sim)))
                err90b = '{0}/{1}'.format(proj, mod.replace('.txt', 
                    '_SEI.err').replace('mod{0}/mod{0}'.format(modn), 
                    'mod{0}/{1}_mod{0}'.format(modn, sim)))
                addtouch += 'touch {0}\n'.format(err90a[:90])
                addtouch += 'touch {0}\n'.format(err90b[:90])
                addtouch += 'touch {0}\n'.format(
                    err90a[:90].replace(".err", ".chk").replace(".er", ".ch").
                    replace(".e", ".c"))
                addtouch += 'touch {0}\n'.format(
                    err90b[:90].replace(".err", ".chk").replace(".er", ".ch").
                    replace(".e", ".c"))
        modchmod = _phc.trimpathname(mod)
        modchmod[1] = modchmod[1].replace('.txt', '*')
        # addtouch += 'chmod 664 {0}/{1}/*{2}\n'.format(proj, *modchmod)

        # Set simulation check variable
        if 3 in cases:
            simchk = _np.ones(len(simulations), dtype=bool)
        else:    
            simchk = _np.zeros(len(simulations), dtype=bool)

        if _os.path.exists(mod.replace('.txt', '{0:02d}.temp'.format(st1max)))\
            and chkout and 1 in cases:
            cases.remove(1)
        case1 = doCase1(inp, cases)
        f0.writelines(case1 + ['\n'])

        if _os.path.exists(mod.replace('.txt', '{0:02d}.temp'.format(
            st1refmax))) and chkout and 2 in cases:
            cases.remove(2)    
        case2 = doCase2(inp, cases)
        f0.writelines(case2 + ['\n'])

        if chkout and 3 in cases:
            for i in range(len(simulations)):
                outs2a = 'mod{0}/{1}_mod{0}{2}.sed2'.format(
                    modn, simulations[i], suf)
                outs2b = 'mod{0}/{1}_mod{0}{2}_SEI.sed2'.format(
                    modn, simulations[i], suf)
                if _os.path.exists(outs2a) or _os.path.exists(outs2b):
                    simchk[i] = False
            if True not in simchk:
                cases.remove(3)
        case3 = doCase3(inp, simchk)
        f0.writelines(case3)
        f0.close()

        # Def automatic walltime:
        if wcheck:
            h = 0 
            if 1 in cases:
                h += 1
            if 2 in cases:
                h += 1
            idx = _np.where(simchk is True)
            if len(idx[0]) > 0:
                extra = 0 + 4 * len(idx[0])
                h = h + extra * 48 / nodes
            walltime = '{0}:0:0'.format(h)

        # Del old jobs
        for sel in clusters:
            outname = mod[mod.find('/') + 1:].replace('txt', sel)
            if _os.path.exists('{0}s/{1}'.format(sel, outname)):
                _os.system('rm {0}s/{1}'.format(sel, outname))

        # Write jobs (if necessary)
        if len(cases) > 0:
            for sel in clusters:
                doJobs(mod, sel, nodes, addtouch)

    if path is not '':
        _os.chdir(path0)
    # PROGRAM END
    return


def makeNoDiskGrid(modn, selsources, path=None):
    """
    #Create a model list with random disk parameters ("noCS" in filename)

    INPUT:  modn = '01'; selsources = '*' (filter that is applied to the SOURCE
    selection).

    OUTPUT: Files written
    """

    def doNoCS(src):
        '''
        Given a src, generates the noCS model input
        '''
        srcname = src.replace('source/', '').replace('.txt', '')
        suffix = '_noCS_{0}'.format(srcname)

        wmod = mod[:]
        # Remove a disk does not work:
        # wmod[9]=wmod[9].replace('1','0')
        wmod[13] = wmod[13].replace('18.6', '2.0')

        f0 = open('mod' + modn + '/mod' + modn + suffix + '.txt', 'w')
        f0.writelines(wmod)
        f0.close()
        return

    # PROGRAM BEGINS
    path0 = _os.getcwd()
    if path is not None:
        _os.chdir(path)
        if path[-1] != '/':
            path += '/'
    else:
        path = ''
    # Check modN folder
    if not _os.path.exists('mod{0}'.format(modn)):
        _os.system('mkdir mod{0}'.format(modn))

    # Select sources
    sources = _glob('source/' + selsources)

    # Load disk model
    f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF_disco.txt'))
    mod = f0.readlines()
    f0.close()

    for prodI in _product(sources):
        prodI = prodI[0]
        doNoCS(prodI)
    print('# {0:.0f} arquivos foram gerados !!'.format(len(sources)))
    if path is not "":
        _os.chdir(path0)
    # END PROGRAM
    return


def makeSimulLine(basesims, Rs, hwidth, Ms, Obs, lsuffix, path='source/'):
    """
    | vrots = [[167.023,229.187,271.072,301.299,313.702],
    |          [177.998,244.636,290.596,324.272,338.298],
    |          [192.612,267.017,318.288,355.320,370.638],
    |          [202.059,281.667,335.158,373.716,389.782],
    |          [209.244,292.409,358.626,410.439,430.844],
    |          [214.407,297.661,357.799,402.628,420.683]]

    | vrots = [[259.759,354.834,417.792,464.549,483.847],
    |          [252.050,346.163,406.388,449.818,468.126],
    |          [245.127,336.834,399.983,448.076,467.806],
    |          [239.522,329.496,388.734,432.532,450.806],
    |          [234.301,321.139,379.297,423.241,441.122],
    |          [228.538,313.797,370.343,412.488,429.914],
    |          [219.126,299.656,354.547,395.821,413.008],
    |          [211.544,288.840,341.081,380.426,396.978],
    |          [203.438,279.328,328.666,365.697,380.660],
    |          [197.823,268.964,316.901,353.568,368.506],
    |          [192.620,262.688,308.208,341.963,356.410],
    |          [187.003,255.125,299.737,332.511,346.043]]
    |         
    | basesims = ['simulation/Brg.txt','simulation/Ha.txt']
    | Rs = [12000, 20000]
    | 
    | Ms = [4.00,5.00,7.00,9.00,12.00,15.00]
    | Ms = [14.6, 12.5, 10.8, 9.6, 8.6, 7.7, 6.4, 5.5, 4.8, 4.2, 3.8, 3.4]
    | Obs = [1.1,1.2,1.3,1.4,1.45]
    | suffix = 'H0.30_Z0.014_bE_Ell'
    """
    c = _phc.c.cgs

    for prodI in _product(Ms, Obs, basesims, lsuffix):
        M, Ob, basesim, suffix = prodI

        f0 = open(basesim)
        mod = f0.readlines()
        f0.close()

        srcid = 'Be_M{0:05.2f}_ob{1:.2f}_{2}.txt'.format(M, Ob, suffix)

        # i = Ms.index(M)
        # j = Obs.index(Ob)
        k = basesims.index(basesim)
        R = Rs[k]
        nmod = mod[:]
        vrot = _rot.vrot_scr(path+srcid)
        vel = '{0:.1f}'.format(hwidth + vrot)
        nmod[103] = nmod[103].replace('1020.', vel)
        n = str(int(round(2 * (hwidth + vrot) * R / c * 1e5)))
        # if len(_np.shape(vrots)) == 2:
        #     vel = '{0:.1f}'.format(hwidth + vrots[i][j])
        #     nmod[103] = nmod[103].replace('1020.', vel)
        #     n = str(int(round(2 * (hwidth + vrots[i][j]) * R / c * 1e5)))
        # else:
        #     vel = '{0:.1f}'.format(hwidth + vrots[0])
        #     nmod[103] = nmod[103].replace('1020.', vel)
        #     n = str(int(round(2 * (hwidth + vrots[0]) * R / c * 1e5)))
        # print(srcid, n)
        nmod[100] = nmod[100].replace('100', n)

        f0 = open(
            basesim.replace('.txt', '_{0}'.format(srcid)), 'w')
        f0.writelines(nmod)
        f0.close()
    return


def makeSourceGrid(masses, rps, lums, Ws, betas, path=None):
    """ Arbitrary values (no stellar model info), for the SELF-CONSISTENT RIGID 
    ROTATOR setup.

    OUTPUT: The combination of the values of all lists in the *path/source/* 
    folder.

    Example: 2 of each parameters results in 32 models (2^5)
    """
    path0 = _os.getcwd()
    if path is not None:
        _os.chdir(path)
        if path[-1] != '/':
            path += '/'
    else:
        path = ''
    # 
    f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF_estrela.txt'))
    mod = f0.readlines()
    f0.close()
    if not _os.path.exists('{0}source'.format(path)):
        _os.system('mkdir {0}source'.format(path))
    # 
    for prod in _product(masses, rps, lums, Ws, betas):
        MI, Raio, Lum, W, Beta = prod
        # REGISTRA VALORES
        wmod = mod[:]        
        wmod[3] = wmod[3].replace('10.3065', ('%.2f' % MI))
        wmod[4] = wmod[4].replace('5.38462', ('%.2f' % Raio))
        wmod[5] = wmod[5].replace('0.775', ('%.4f' % W))            
        wmod[6] = wmod[6].replace('7500.', ('%.2f' % Lum))
        wmod[7] = wmod[7].replace('0.25', ('%.5f' % Beta))
        # 
        suffix = '_M{0:05.2f}_W{1:.2f}_b{2:.2f}_rp{3:05.2f}_L{4:05.0f}_Ell'. \
            format(MI, W, Beta, Raio, Lum)
        f0 = open('{0}source/Be'.format(path) + suffix + '.txt', 'w')
        f0.writelines(wmod)
        f0.close()
    #
    if path is not "":
        _os.chdir(path0)    
    return 


def makeStarGrid(oblats, Hfs, path=None):
    """
    Based on GENEVA models.

    | INPUT: oblats = [1.1,1.2,1.3,1.4,1.45] (example)
    | Hfs = [0.3] (example)

    Masses list a Z value are inside `geneve_par.pro` file.
    """
    path0 = _os.getcwd()
    if path is not None:
        _os.chdir(path)
        if path[-1] != '/':
            path += '/'
    else:
        path = ''
    if not _os.path.exists('stmodels'):
        _os.system('mkdir stmodels')
    try:
        runIDL = True 
        import pidly
    except ImportError:
        print('# This system do not have pIDLy installed...')
        runIDL = False

    if runIDL:
        key = _phc.user_input('# Do you want to run "geneve_par" (y/other):')
        if key != 'y':
            runIDL = False

    if runIDL:
        idl = pidly.IDL()
        propath = _os.path.join( _hdt.hdtpath(), 'refs' )
        idl('cd,"{0}"'.format(propath))
        idl('.r geneve_par')
        for ob in oblats:
            for H in Hfs:
                idl('geneve_par, {0}, {1}, /OBLAT,/makeeps'.format(ob, H))
                _os.system('mv {0}/geneve_lum.eps stmodels/geneve_lum_' +
                    '{1:.2f}_{2:.2f}.eps'.format(propath, ob, H))
                _os.system(
                    'mv {0}/geneve_rp.eps stmodels/geneve_rp_{1:.2f}_' +
                    '{2:.2f}.eps'.format(propath, ob, H))
                _os.system(
                    'mv {0}/geneve_par.txt stmodels/oblat{1}_h{2}.txt'.
                    format(propath, ob, H))
        idl.close()

    f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF_estrela.txt'))
    mod = f0.readlines()
    f0.close()

    if not _os.path.exists('source'):
        _os.system('mkdir source')

    for ob in oblats:                        
        for H in Hfs:  
            f0 = open('stmodels/oblat{0}_h{1}.txt'.format(ob, H))                     
            matriz = f0.readlines() 
            f0.close()
            Omega, W, Beta = map(float, matriz[1].split())

            m2 = []
            for i in range(4, len(matriz)):
                if len(matriz[i]) > 1:
                    m2 += [matriz[i].split()[1:]]
            matriz = _np.array(m2, dtype=float)
            M = matriz[:, 0]  # MASS (SOLAR MASSES)
            M = list(M)
            Rp = matriz[:, 1]  # POLAR RADIUS (SOLAR RADII)
            Rp = list(Rp)
            L = matriz[:, 2]  # LUMINOSITY (in solar lum.)
            L = list(L)
            Z = [0.014]  # METALLICITY(=Zsolar) (other options: 0.006, 0.002)

            print('Omega = ', Omega)
            print('W     = ', W)
            print('beta  = ', Beta)
            print('M     = ', M)
            print('Rp    = ', Rp)
            print('L     = ', L)
            print("%.0f arquivos gerados\n" % (len(M) * len(Hfs)))

            # DEFINE ALL INDEX
            for prodI in _product(M, Rp, L, Z):
                MI, RpI, LI, ZI = prodI
                a = M.index(MI)
                Raio = Rp[a]
                Lum = L[a]

                suffix = '_M{0:05.2f}_ob{1:.2f}_H{2:.2f}_Z{3}_bE_Ell'. \
                    format(MI, ob, H, ZI, Beta, RpI, LI)

                # REGISTRA VALORES
                wmod = mod[:]        
                wmod[3] = wmod[3].replace('10.3065', ('%.2f' % MI))
                wmod[4] = wmod[4].replace('5.38462', ('%.2f' % Raio))
                wmod[5] = wmod[5].replace('0.775', ('%.4f' % W))            
                wmod[6] = wmod[6].replace('7500.', ('%.2f' % Lum))
                wmod[7] = wmod[7].replace('0.25', ('%.5f' % Beta))

                f0 = open('source/Be' + suffix + '.txt', 'w')
                f0.writelines(wmod)
                f0.close()
    #
    if path is not "":
        _os.chdir(path0)
    return


def makeSimulDens(dbase, basesim):
    """
    Sets the SED simulations number of photos so that the signal/noise level
    is approximately constant at visible polarization.

    |dbase = _np.logspace(_np.log10(0.02),_np.log10(4.0),7)
    |basesim = 'simulation/sed.txt'
    """
    f0 = open(basesim)
    mod = f0.readlines()
    f0.close()
    # fact = 2.  Tempo execucao = d/1e13*fact
    # Nf0 = 500000000
    for d in dbase:
        srcid = 'sig{0:.2f}'.format(d)
        # alpha = .39794
        # beta =  13.87219
        alpha = 0.34588
        beta = 8.50927
        newd = int(10**(-alpha * _np.log10(d) + beta))
        print('{0}, N_f = {1:.2f}e+9'.format(srcid, newd / 1e9))
        nmod = mod[:]
        nmod[9] = nmod[9].replace('500000000', '{0}'.format(newd))
        f0 = open(basesim.replace('.txt', '_{0}.txt'.format(srcid)), 'w')
        f0.writelines(nmod)
        f0.close()
    # a = raw_input('asdads')
    return 


def makeCSGrid_bistabWind1Dust(modn='01', renv=[18.6], rcs=[5.], 
    dm_dOmega=[1.5e-5], v_0=[10.], v_inf=[1200.], beta_v=[2], A1=[49], 
    A2=[-0.7], A3=[0.], m_dth=[92.], grain_dust_ratio=[200.], grain_dens=[1.], 
    selsources='*', path=None):
    """ Great a CS Grid for HDUST for wind+dust shell (1) both in bi-stability

    Based on Carciofi+2010 and used in the B[e] grid between IAG/ON (Brazil).

    :param arg1: the first value
    :param arg2: the first value
    :type arg1: int, float,...
    :type arg2: int, float,...
    :returns: arg1/arg2 +arg3
    :rtype: int, float

    :Example:

    >>> import template
    >>> a = template.MainClass1()
    >>> a.function1(1,1,1)
    2

    .. note:: can be useful to emphasize important feature
    .. seealso:: :class:`MainClass2`
    .. warning:: arg2 must be non-zero.
    .. todo:: check that arg2 is non zero.
    """
    def dobistabWind1Dust(lpars):
        '''
        Given a proper list of parameters, do the thing...
        '''
        dm_d0, v0, vinf, beta, a1, a2, a3, mdth, gdratio, gdens, re, rd, src =\
            lpars

        srcname = _os.path.splitext(_os.path.basename(src))[0]
        # TODO
        suffix = ('_rd{0:02.0f}_a1+{1:02.0f}_m{2:03.0f}_md{3:.0e}_grh{4:04.1f}'
            '_gdr{5:03.0f}_v0+{6:02.0f}_vinf{7:04.0f}_a2+{8:04.1f}_b{9:03.1f}'
            '_a3+{10:03.1f}').format(rd, a1, mdth, dm_d0, gdens, gdratio, v0, 
            vinf, a2, beta, a3, srcname)

        wmod = mod[:]
        wmod = _phc.repl_fline_val(wmod, 13, '18.6', re)
        wmod = _phc.repl_fline_val(wmod, 14, '3.', rd)
        wmod = _phc.repl_fline_val(wmod, 19, '9.', a1)
        wmod = _phc.repl_fline_val(wmod, 20, '182.', mdth)
        wmod = _phc.repl_fline_val(wmod, 26, '39.', a1)
        wmod = _phc.repl_fline_val(wmod, 27, '182.', mdth)
        wmod = _phc.repl_fline_val(wmod, 35, '2.E-9', dm_d0)
        wmod = _phc.repl_fline_val(wmod, 30, '1.E-7', dm_d0)
        wmod = _phc.repl_fline_val(wmod, 31, '1.', gdens)
        wmod = _phc.repl_fline_val(wmod, 32, '200.', gdratio)
        wmod = _phc.repl_fline_val(wmod, 47, '10.', v0)
        wmod = _phc.repl_fline_val(wmod, 48, '400.', vinf)
        wmod = _phc.repl_fline_val(wmod, 49, '-0.7', a2)
        wmod = _phc.repl_fline_val(wmod, 50, '0.8', beta)
        wmod = _phc.repl_fline_val(wmod, 51, '2.75', a3)

        fmod = _os.path.join('mod'+modn, 'mod'+modn+suffix+'.txt')
        f0 = open(fmod, 'w')
        f0.writelines(wmod)
        f0.close()
        return

    # PROGRAM BEGINS
    path0 = _os.getcwd()
    if path is not None:
        _os.chdir(path)
    else:
        path = ''
    # Check modN folder
    if not _os.path.exists('mod{0}'.format(modn)):
        _os.system('mkdir mod{0}'.format(modn))

    # Select sources
    # print _os.path.join('source', selsources)
    sources = _glob(_os.path.join('source', selsources))

    # Load disk model
    f0 = open(_os.path.join(_hdt.hdtpath(), 'refs', 'REF_bistabWind1Dust.txt'))
    mod = f0.readlines()
    f0.close()

    for lpars in _product(dm_dOmega, v_0, v_inf, beta_v, A1, A2, A3, m_dth, 
        grain_dust_ratio, grain_dens, renv, rcs, sources):
        dobistabWind1Dust(lpars)

    if path is not "":
        _os.chdir(path0)
    # END PROGRAM
    return

# MAIN ###
if __name__ == "__main__":
    pass
