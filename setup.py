#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys


def check_py_ver():
    if (sys.version_info[0] > 2):
        return
    elif (sys.version_info[1] >= 7):
        return
    else: 
        _ = raw_input('Press Ctlr+Z to abort or ENTER to continue... ')
        return

print('# Warning! PyHdust in only supported for Python 2.7+\n'
    'Press Ctlr+Z to abort or ENTER to continue... ')
check_py_ver()

try:
    from setuptools import setup, find_packages
    errimport = False
except ImportError:
    from distutils.core import setup
    from glob import glob
    errimport = True

    def is_package(path):
        return ( os.path.isdir(path) and os.path.isfile(os.path.join(path, 
            '__init__.py')) )

    def find_packages(path=".", base="", exclude=[]):
        """ Find all packages in path """
        packages = {}
        lexc = []
        for item in exclude:
            lexc.extend(glob(item))
        for item in os.listdir(path):
            dir = os.path.join(path, item)
            if is_package( dir ) and item not in lexc:
                if base:
                    module_name = "%(base)s.%(item)s" % vars()
                else:
                    module_name = item
                packages[module_name] = dir
                packages.update(find_packages(dir, module_name))
        return packages


def rd(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r


def recfiles(chdir, path):
    out = []
    opath = os.getcwd()
    os.chdir(chdir)
    for root, dirs, fnames in os.walk(path):
        for f in fnames:
            out.append(os.path.join(root, f))
    os.chdir(opath)
    return out

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

setup(name='pyhdust',
    version='1.0.9',
    description=('Analysis tools for multi-technique astronomical data and '
        'hdust models'),
    url='http://j.mp/pyhdust',
    author='Daniel M. Faes',
    author_email='dmfaes@gmail.com',
    license='GNU GPLv3.0',      
    # packages=['pyhdust','pyhdust_refs'],
    scripts=[os.path.join('scripts', f) for f in os.listdir('scripts') if 
        os.path.splitext(f)[1] == '.py'],
    packages=find_packages(exclude=['build', 'docs', '*egg*', 'dist', 
        'scripts']),
    # include=recfiles('refs'),
    package_data={'pyhdust': recfiles('pyhdust', 'refs') + 
        ['LICENSE', 'README.rst']},
    include_package_data=True,
    # data_files = [('refs/*', 'stmodels/*')],
    # package_dir = {'../'},
    zip_safe=False,
    install_requires=['numpy >= 1.9.0', 'six'],
    long_description=rd('README.rst'),
    classifiers=[
        # "Development Status :: 4 - Beta",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later" + \
        " (GPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    keywords=['astronomy', 'astrophysics', 'science', 'hdust', 'Be stars',
        'spectroscopy', 'polarimetry', 'interferometry', 'radiative transfer',
        'optical-interferometry'],
)

if errimport:
    print('# You don\'t have "setuptools" installed!')
    print('# Because of this, you MAY need to copy the package data: \n')
    print('# Warning! The path can change according to your system')
    print('$ cp -r -f pyhdust ~/.local/lib/python2.7/site-packages/')
