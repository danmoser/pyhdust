#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Modified by D. Moser in 2015-03-17

import pyhdust

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def rd(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r


if __name__ == "__main__":
    setup(name='pyhdust',
    version=pyhdust.__version__,
    description='BeACoN’s Python tools for Hdust',
    url='http://astroweb.iag.usp.br/~moser/pyhdust/',
    author='Daniel M. Faes',
    author_email='dmfaes@gmail.com',
    license='GNU GPLv3.0',
    packages=['pyhdust'],
    zip_safe=False,
    install_requires=['numpy'],
    # ~ package_dir = {'../'},
    long_description=rd('README.md'),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        ],
    )
    
