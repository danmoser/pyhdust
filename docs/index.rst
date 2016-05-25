.. HDUST Python Tools documentation master file, created by
	sphinx-quickstart on Wed Oct 22 16:33:29 2014.
	You can adapt this file completely to your liking, but it should at least
	contain the root `toctree` directive.

Welcome to BeACoN's Python Tools (*PyHdust*) documentation!
=============================================================

Pyhdust is currently at **version 0.994**.

About `hdust` code: 

- http://adsabs.harvard.edu/abs/2006ApJ...639.1081C, 
- http://adsabs.harvard.edu/abs/2008ApJ...684.1374C

Package contents:

.. toctree::
    :maxdepth: 2

    pyhdust
    bcd
    beatlas
    fieldstars
    hdrpil
    images
    input
    interftools
    jdcal
    oifits
    phc
        phc_list
    poltools
    releases
    rotstars
    singscat
    spectools
    stats
    tabulate
    triangle

It should be independent of plataform (Linux, Mac, Windows) and compatible with any version of Python (superior to 2.7).


How to install/uninstall
================================
I strongly suggest users to use the ``pip`` plataform:

.. code-block:: bash

    pip install pyhdust

If your are not the root of the system, add the flag ``--user`` to the command above. Dependencies (as filters efficiency and stellar evolution models) are included in the subdirectories.

PyHdust requires the numpy module. Optionally, it makes use of: 

    - matplotlib
    - pyfits
    - emcee
    - pIDLy
    - scipy
    - pyqt_fit

To only **update** the package:

.. code-block:: bash

    pip install -U --no-deps pyhdust

``-U`` forces the upgrade and ``--no-deps`` do not reinstall the dependent packages. 

To uninstall it:

.. code-block:: bash

    pip uninstall pyhdust
 

Alternative downloads
---------------------------
git
^^^^
The package is available under **git** platform (https://github.com/danmoser/pyhdust/)

.. _https://github.com/danmoser/pyhdust/: https://github.com/danmoser/pyhdust/

.. code-block:: bash

    git clone https://github.com/danmoser/pyhdust.git

The command above saves the package in *pyhdust* subfolder. 

To **update** the package with *git*:

.. code-block:: bash

    git pull


zipball
^^^^^^^^
If neither *pip* or *git* are not installed in your system, you can catch the current version as a zipball:

.. code-block:: bash

    curl -L https://github.com/danmoser/pyhdust/zipball/master > pyhdust.zip
    unzip pyhdust.zip


Alternative installations
---------------------------
In the downloaded pyhdust folder, type:

.. code-block:: bash

    python setup.py install

If your are not the root of the system, add the flag ``--user`` to the command
above.

**Alternatively**, to import the modules declare in your *~/.bashrc* file:

.. code-block:: bash
    
    export PYTHONPATH=$PYTHONPATH:~/pyhdust/

where *~/pyhdust/* is the package installation directory.

.. note::

    To uninstall pyhdust from the alternative installations, using *git* you need to remove its folder manually. Using the zipball installation, you can use the following commands:

.. code-block:: bash

    python setup.py install --record files.txt
    cat files.txt | xargs rm -rf


How to use the modules
===============================================
To make use of all routines, the suggestion is to import them as follows:

.. code:: python

    import pyhdust as hdt
    import pyhdust.beatlas as bat
    import pyhdust.fieldstars as fls
    import pyhdust.images as img
    import pyhdust.input as inp
    import pyhdust.interftools as intt
    import pyhdust.jdcal as jdcal
    import pyhdust.poltools as polt
    import pyhdust.phc as phc
    impoty pyhdust.rotstars as rot
    import pyhdust.singscat as sst
    import pyhdust.spectools as spt
    import pyhdust.stats as stt


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* `Source code index`_

.. _`Source code index`: _modules/index.html


License
==========

The code is free, available under the terms of the `GNU GPL v3.0 license`_.

.. _GNU GPL v3.0 license: 
	https://github.com/danmoser/pyhdust/blob/master/LICENSE


.. raw:: html
    
    <a href="http://github.com/danmoser/pyhdust/" target="_blank"><img 
    style="position: fixed; top: 0; right: 0; border: 0;"  
    src="http://s3.amazonaws.com/github/ribbons/forkme_right_orange_ff7600.png" 
    alt="Fork me on GitHub"></a>
