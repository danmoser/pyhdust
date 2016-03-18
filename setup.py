#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Modified by D. Moser in 2015-04-15

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages


def rd(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r


if __name__ == "__main__":
    setup(name='pyhdust',
    version=0.97,
    description='BeACoN’s Python tools for Hdust',
    url='http://astroweb.iag.usp.br/~moser/doc/',
    author='Daniel M. Faes',
    author_email='dmfaes@gmail.com',
    license='GNU GPLv3.0',
    #~ packages=['pyhdust'],
    packages=find_packages(exclude=['build', 'docs', '*egg*', 'dist']),
    package_data={'pyhdust': ['../LICENSE', '../README.md', '../filters/*',
    '../refs/*', '../stmodels/*', '../stars/*']},
    #~ include_package_data=True,
    zip_safe=False,
    install_requires=['numpy'],
    #install_requires=['numpy >= 1.6.0'],
    #~ data_files = [('refs/*', 'stmodels/*')],
    #~ package_dir = {'../'},
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
    
