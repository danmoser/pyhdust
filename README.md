pyhdust: Python tools for `hdust`
=============================================================

Full documentation at:
https://dl.dropboxusercontent.com/u/6569986/doc/index.html

About hdust` code: 
- http://adsabs.harvard.edu/abs/2006ApJ...639.1081C, 
- http://adsabs.harvard.edu/abs/2008ApJ...684.1374C


How to download the files:
================================
The package is available under **git** platform
(https://github.com/danmoser/pyhdust/). Dependencies (as filters
efficiency and stellar evolution models) are included in the subdirectories.

.. _https://github.com/danmoser/pyhdust/: https://github.com/danmoser/pyhdust/

git
-------

To download it:

.. code-block:: bash

    $ git clone https://github.com/danmoser/pyhdust.git

The command above saves the package in *pyhdust* subfolder. 

To **update** the package with *git*:

.. code-block:: bash

    $ git pull

zipball
---------

If *git* is not installed in your system, you can catch the current version
zipball:

.. code-block:: bash

    $ curl -L https://github.com/danmoser/pyhdust/zipball/master > pyhdust.zip
    $ unzip pyhdust.zip

pip
---

.. role:: strike

It is hoped that in a not too distant future pyHdust will be available in the
*pip* platform...


How to install
=================

In the downloaded pyHdust folder, type:

.. code-block:: bash

    $ python setup.py install

If your are not the root of the system, add the flag ``--user`` to the command
above.

**Alternatively**, to import the modules declare in your *~/.bashrc* file:

.. code-block:: bash
    
    export PYTHONPATH=$PYTHONPATH:~/pyhdust/

where *~/pyhdust/* is the package installation directory.

PyHdust requires the numpy module. Optionally, it makes use of: 

    - matplotlib
    - pyfits
    - emcee
    - pIDLy
    - scipy


How to use the tools
===============================================
To make use of all routines, the suggestion is to import them as follows:

.. code:: python

    import pyhdust as hdt
    import pyhdust.beatlas as bat
    import pyhdust.fieldstars as fls
    import pyhdust.input as inp
    import pyhdust.interftools as intt
    import pyhdust.poltools as polt
    import pyhdust.phc as phc
    import pyhdust.singscat as sst
    import pyhdust.spectools as spt


License
==========

The code is free, available under the terms of the `GNU GPL v3.0 license`_.

.. _GNU GPL v3.0 license: 
	https://github.com/danmoser/pyhdust/blob/master/LICENSE