#!/usr/bin/env python
# -*- coding:utf-8 -*-

try:
    from setuptools import setup, find_packages
    errimport = False
except ImportError:
#~ if True:
    from distutils.core import setup
    import os
    from glob import glob
    errimport = True
    
    def is_package(path):
        return (
        os.path.isdir(path) and
        os.path.isfile(os.path.join(path, '__init__.py'))
        )

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


if __name__ == "__main__":
    setup(name='pyhdust',
    version='0.971',
    description='BeACoN’s Python tools for Hdust',
    url='https://dl.dropboxusercontent.com/u/6569986/doc/index.html',
    author='Daniel M. Faes',
    author_email='dmfaes@gmail.com',
    license='GNU GPLv3.0',
    #~ packages=['pyhdust','pyhdust_refs'],
    packages=find_packages(exclude=['build', 'docs', '*egg*', 'dist']),
    include_package_data=True,
    #~ include=['pyhdust_refs']),
    #~ package_data={'pyhdust':['*']},
    #~ , '../filters/*', '../refs/*', '../stmodels/*']},
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
        "License :: OSI Approved :: GNU General Public License v3 or later" + \
        " (GPLv3+)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        ],
    )

    if errimport:
        print('# You don\'t have "setuptools" installed!')
        print('# Because of this, you need to ADAPT and run this command: \n')
        print('# Warning! The cmd path MAY change according to your system')
        print('$ cp -r -f pyhdust ~/.local/lib/python2.7/site-packages/')
